/*
 *  @BEGIN LICENSE
 *
 *  Hilbert: a space for quantum chemistry plugins to Psi4
 *
 *  Copyright (c) 2020 by its authors (LICENSE).
 *
 *  The copyrights for code used from other parties are included in
 *  the corresponding files.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 *  @END LICENSE
 */

#include <tiledarray.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsio/psio.hpp>
#include <utility>

#include <psi4/libtrans/integraltransform.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libfunctional/superfunctional.h>
#include "psi4/libscf_solver/uhf.h"
#include <psi4/libmints/mintshelper.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/basisset.h>
#include <psi4/lib3index/dftensor.h>
#include <psi4/libqt/qt.h>

#include "cc_cavity/misc/ta_helper.h"
#include "misc/threeindexintegrals.h"

#include <mkl.h>
#include <omp.h>
#include "misc/blas.h"
#include "misc/hilbert_psifiles.h"
#include "polaritonic_scf/uhf.h"
#include <unistd.h>
#include <psi4/psifiles.h>
#include "cc_cavity/include/cc_cavity.h"

using namespace std;
using namespace TA;
using namespace psi;
using namespace TA_Helper;

namespace hilbert {

    CC_Cavity::CC_Cavity(std::shared_ptr<Wavefunction> reference_wavefunction, Options& options_):
            PolaritonicHF(std::move(reference_wavefunction), options_) {

        // set block_size for tiledarray
        TA_Helper::tile_size_ = (size_t) options_.get_int("TILE_SIZE");

        // adjustment of threaded mkl flags to work in parallel with mpi and omp
        // will default to mkl sequential if threads equal one
        int numThreads = Process::environment.get_n_threads();
        if ( numThreads > 1 ) {
            mkl_set_threading_layer(0);
            omp_set_num_threads(numThreads);
            mkl_set_num_threads(numThreads);
        } else {
            mkl_set_threading_layer(1);
            omp_set_num_threads(numThreads);
            mkl_set_num_threads(numThreads);
        }

        // if `has_photons` is true, but `n_photon_states` is not two, print warning
        if ( has_photon_ && n_photon_states_ != 2 )
            Printf("WARNING: CC Cavity should have two photon states if using bosonic operators\n");

        world_.gop.fence();
        common_init();
    }

    void CC_Cavity::common_init(){
        Printf("\n\n");
        Printf( "        *******************************************************\n");
        Printf( "        *                                                     *\n");
        Printf( "        *                      -----------                    *\n");
        Printf( "        *                       CC CAVITY                     *\n");
        Printf( "        *                      -----------                    *\n");
        Printf( "        *                                                     *\n");
        Printf( "        *******************************************************\n");
        Printf("\n\n");

        /// grab dimensions
        // ensure scf_type df
        if ( options_.get_str("SCF_TYPE") != "DF" && options_.get_str("SCF_TYPE") != "CD" ) {
            throw PsiException("CC Cavity only works with scf_type df or cd for now",__FILE__,__LINE__);
        }

        // ensure running_ in c1 symmetry
        if ( reference_wavefunction_->nirrep() > 1 ) {
            throw PsiException("CC Cavity only works with c1 symmetry for now.",__FILE__,__LINE__);
        }

        // grab the number of occupied and virtual orbitals
        o_ = nalpha_ + nbeta_;
        v_ = (nmo_-nalpha_) + (nmo_-nbeta_);

        // grab the number of occupied and virtual orbitals for each spin
        oa_ = nalpha_;
        ob_ = nbeta_;
        va_ = nmo_ - oa_;
        vb_ = nmo_ - ob_;

        // grab the number of spin orbitals
        size_t n  = 2L*(size_t)nmo_;
        ns_  = 2L*(size_t)nso_;

        if ( n != ns_ ) {
            throw PsiException("cc_cavity does not work with a basis that has linear dependencies", __FILE__, __LINE__);
        }

        // initialize diis solver
        int max_vecs = options_.get_int("DIIS_MAX_VECS");
        diis_ta = (std::shared_ptr<DIISTA>)(new DIISTA(max_vecs));

        /// cavity options

        // update cavity terms once more
        if ( n_photon_states_ > 1) {
            update_cavity_terms();
        }
        same_a_b_orbs_ = same_a_b_dens_ = false;

        cc_type_ = has_photon_ ? "QED-CC" : "CC";
        cc_type_ += "SD";
        if (include_t3_ || include_u3_) cc_type_ += "T";
        if (include_t4_ || include_u4_) cc_type_ += "Q";
        if (has_photon_) cc_type_ += "-1";

        // set coupling factors
        lambda_[0] = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
        lambda_[1] = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
        lambda_[2] = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

        /// initialize containers

        // initialize arbitrary index strings and identity tensors
        init_operators();

        // make containers for integrals
        init_integrals();

        // initialize orbital energies
        epsilon_ = (double*) calloc(2*nso_, sizeof(double));

        // transform integrals to MO basis
        transform_integrals(false);

    }

    void CC_Cavity::init_operators() {
        /// initialize arbitrary index strings
        idx_map_ = std::vector<std::string>(25);
        for (int i = 1; i < 25; i++) { // hard-coded for now. Unlikely to need more than 25 indices.

            // build string of indices for i'th rank tensor
            for (int j = 0; j < i; j++) idx_map_[i] += "x" + std::to_string(i) + ",";

            // remove last comma
            idx_map_[i].pop_back();
        }

        /// make identity tensors

        Id_blks_["a_o"] = makeTensor(world_, {oa_}, false);   Id_blks_["a_o"].fill(1.0);
        Id_blks_["a_v"] = makeTensor(world_, {va_}, false);   Id_blks_["a_v"].fill(1.0);
        Id_blks_["b_o"] = makeTensor(world_, {ob_}, false);   Id_blks_["b_o"].fill(1.0);
        Id_blks_["b_v"] = makeTensor(world_, {vb_}, false);   Id_blks_["b_v"].fill(1.0);

        Id_blks_["aa_oo"] = TiledArray::diagonal_array<TA::TArrayD>(world_, makeRange({oa_, oa_}), 1);
        Id_blks_["aa_vv"] = TiledArray::diagonal_array<TA::TArrayD>(world_, makeRange({va_, va_}), 1);
        Id_blks_["bb_oo"] = TiledArray::diagonal_array<TA::TArrayD>(world_, makeRange({ob_, ob_}), 1);
        Id_blks_["bb_vv"] = TiledArray::diagonal_array<TA::TArrayD>(world_, makeRange({vb_, vb_}), 1);

        Id_blks_["aaaa_oooo"]("p,q,r,s") = Id_blks_["aa_oo"]("p,r") * Id_blks_["aa_oo"]("q,s");
        Id_blks_["abab_oooo"]("p,q,r,s") = Id_blks_["aa_oo"]("p,r") * Id_blks_["bb_oo"]("q,s");
        Id_blks_["bbbb_oooo"]("p,q,r,s") = Id_blks_["bb_oo"]("p,r") * Id_blks_["bb_oo"]("q,s");

        /// initialize t1 and t2 amplitudes and residuals (every derived class must use these amplitudes and residuals)

        // t1
        amplitudes_["t1_aa"]  = makeTensor(world_, {va_, oa_}, true);
        amplitudes_["t1_bb"]  = makeTensor(world_, {vb_, ob_}, true);
        residuals_["t1_aa"]  = makeTensor(world_, {va_, oa_}, true);
        residuals_["t1_bb"]  = makeTensor(world_, {vb_, ob_}, true);

        // t2
        amplitudes_["t2_aaaa"]  = makeTensor(world_, {va_,va_, oa_,oa_}, true);
        amplitudes_["t2_abab"]  = makeTensor(world_, {va_,vb_, oa_,ob_}, true);
        amplitudes_["t2_bbbb"]  = makeTensor(world_, {vb_,vb_, ob_,ob_}, true);
        residuals_["t2_aaaa"]  = makeTensor(world_, {va_,va_, oa_,oa_}, true);
        residuals_["t2_abab"]  = makeTensor(world_, {va_,vb_, oa_,ob_}, true);
        residuals_["t2_bbbb"]  = makeTensor(world_, {vb_,vb_, ob_,ob_}, true);

        world_.gop.fence();
    }

    void CC_Cavity::init_integrals() {

        /// build MO transformation matrix
        size_t ns = ns_, oa = oa_, ob = ob_;

        // grab the MO coefficients
        double ** ca = Ca_->pointer();
        double ** cb = Cb_->pointer();

        C_blks_["a_o"] = makeTensor(world_, {ns_, oa_}, false);
        C_blks_["b_o"] = makeTensor(world_, {ns_, ob_}, false);
        C_blks_["a_v"] = makeTensor(world_, {ns_, va_}, false);
        C_blks_["b_v"] = makeTensor(world_, {ns_, vb_}, false);

        // copy the MO coefficients into the tensor
        size_t nso = nso_;
        C_blks_["a_o"].init_elements([ca, nso](auto& I) {
            if (I[0] < nso) return ca[I[0]][I[1]];
            else return 0.0;
        });
        C_blks_["b_o"].init_elements([cb, nso](auto& I) {
            if (I[0] >= nso) return cb[I[0]-nso][I[1]];
            else return 0.0;
        });
        C_blks_["a_v"].init_elements([ca, oa, nso](auto& I) {
            if (I[0] < nso) return ca[I[0]][I[1]+oa];
            else return 0.0;
        });
        C_blks_["b_v"].init_elements([cb, ob, nso](auto& I) {
            if (I[0] >= nso) return cb[I[0]-nso][I[1]+ob];
            else return 0.0;
        });

        /// core hamiltonian
        double ** h_p = reference_wavefunction_->H()->pointer();
        core_H_ = makeTensor(world_, {ns_,ns_}, false);
        core_H_.init_elements([h_p, nso](auto& I) {
            if (I[0] < nso && I[1] < nso) return h_p[I[0]][I[1]];
            else if (I[0] >= nso && I[1] >= nso) return h_p[I[0]-nso][I[1]-nso];
            else return 0.0;
        });

        if (has_photon_) {
            /// extra one-electron terms introduced by cavity
            // e-(n-<d>) contribution to dipole self energy
            double **scaled_e_n_dipole_squared_p = scaled_e_n_dipole_squared_->pointer();

            // one-electron part of electron-electron contribution to dipole self energy
            double **quadrupole_scaled_sum_p = quadrupole_scaled_sum_->pointer();

            oe_cavity_terms_ = makeTensor(world_, {ns_, ns_}, false);
            oe_cavity_terms_.init_elements([scaled_e_n_dipole_squared_p, quadrupole_scaled_sum_p, nso](auto &I) {
                if (I[0] < nso && I[1] < nso)
                    return scaled_e_n_dipole_squared_p[I[0]][I[1]] - quadrupole_scaled_sum_p[I[0]][I[1]];
                else if (I[0] >= nso && I[1] >= nso)
                    return scaled_e_n_dipole_squared_p[I[0] - nso][I[1] - nso] -
                           quadrupole_scaled_sum_p[I[0] - nso][I[1] - nso];
                else return 0.0;
            });
        }

        /// get number of auxiliary basis functions

        if ( options_.get_str("SCF_TYPE") == "DF" ) {
            // get auxiliary basis:
            std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_CC");

            // total number of auxiliary basis functions
            nQ_ = auxiliary->nbf();
        } else if ( options_.get_str("SCF_TYPE") == "CD" ) {

            std::shared_ptr<PSIO> psio(new PSIO());

            psio->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
            psio->read_entry(PSIF_DFSCF_BJ, "length", (char*)&nQ_, sizeof(long int));
            psio->close(PSIF_DFSCF_BJ,1);

        }

        /// initialize 3-index integral blocks

        /// fill 3-index integral blocks
        if ( options_.get_str("SCF_TYPE") == "DF" ) {
            // get primary/auxiliary basis:
            std::shared_ptr<BasisSet> primary = reference_wavefunction_->get_basisset("ORBITAL");
            std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_CC");

            nQ_ = auxiliary->nbf(); // total number of auxiliary basis functions

            // three-index integrals
            std::shared_ptr<DFTensor> DF (new DFTensor(primary,auxiliary,
                                                       Ca_,nalpha_,nso_-nalpha_,nalpha_,nso_-nalpha_,options_));

            std::shared_ptr<Matrix> Qso = DF->Qso();
            double ** Qso_p = Qso->pointer();
            Qso_ = makeTensor(world_, {nQ_, ns_, ns_}, false);
            world_.gop.fence();
            Qso_.init_elements([Qso_p, nso](auto &I) {
                size_t Q = I[0], mu = I[1], nu = I[2];
                if (mu < nso && nu < nso) {
                    return Qso_p[Q][mu*nso + nu];
                } else if (mu >= nso && nu >= nso) {
                    mu -= nso; nu -= nso;
                    return Qso_p[Q][mu*nso + nu];
                } else return 0.0;
            });
            world_.gop.fence();
            Qso.reset();
        } else if ( options_.get_str("SCF_TYPE") == "CD" ) {
            long int nQ = static_cast<size_t>(nQ_);
            ThreeIndexIntegrals(reference_wavefunction_,nQ,memory_);
            nQ_ = static_cast<size_t>(nQ);

            double * Qso_p = (double*)malloc(nso_*nso_*nQ_*sizeof(double));
            memset((void*)Qso_p,'\0',nso_*nso_*nQ_*sizeof(double));

            std::shared_ptr<PSIO> psio(new PSIO());
            psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
            psio->read_entry(PSIF_DCC_QSO,"(Q|mn) Integrals",(char*)Qso_p,sizeof(double)*nQ_ * nso_*nso_);
            psio->close(PSIF_DCC_QSO,1);

            Qso_ = makeTensor(world_, {nQ_, ns_, ns_}, false);
            world_.gop.fence();
            Qso_.init_elements([Qso_p, nso](auto &I) {
                size_t Q = I[0], mu = I[1], nu = I[2];
                if (mu < nso && nu < nso) {
                    return Qso_p[Q*nso*nso + mu*nso + nu];
                } else if (mu >= nso && nu >= nso) {
                    mu -= nso; nu -= nso;
                    return Qso_p[Q*nso*nso + mu*nso + nu];
                } else return 0.0;
            });
            world_.gop.fence();
            free(Qso_p);
        }
        world_.gop.fence();
    }

    void CC_Cavity::transform_integrals(bool use_t1) {
        has_t1_integrals_ = use_t1;

        if (!use_t1)
            return apply_transform(C_blks_, C_blks_);

        // copy C_blks_ to CL and CR
        TArrayMap CL, CR;
        for (auto& block : C_blks_) {
            CL[block.first] = C_blks_[block.first].clone();
            CR[block.first] = C_blks_[block.first].clone();
        }

        // grab references to t1 blocks
        TArrayD &t1_aa = amplitudes_["t1_aa"];
        TArrayD &t1_bb = amplitudes_["t1_bb"];

        // apply t1 to CL and CR
        CL["a_v"]("mu, a") -= C_blks_["a_o"]("mu, i") * t1_aa("a, i");
        CL["b_v"]("mu, a") -= C_blks_["b_o"]("mu, i") * t1_bb("a, i");
        CR["a_o"]("mu, i") += C_blks_["a_v"]("mu, a") * t1_aa("a, i");
        CR["b_o"]("mu, i") += C_blks_["b_v"]("mu, a") * t1_bb("a, i");

        apply_transform(CL, CR);

    }

    void CC_Cavity::print_dimensions() {

        /// compute the number of amplitudes
        singleDim_ = oa_ * va_
                     + ob_ * vb_;
        doubleDim_ = oa_ * oa_ * va_ * va_
                     + ob_ * ob_ * vb_ * vb_
                     + oa_ * ob_ * va_ * vb_;
        if (include_t3_ || include_u3_)
            tripleDim_ = oa_ * oa_ * oa_ * va_ * va_ * va_
                         + ob_ * ob_ * ob_ * vb_ * vb_ * vb_
                         + oa_ * oa_ * ob_ * va_ * va_ * vb_
                         + oa_ * ob_ * ob_ * va_ * vb_ * vb_;
        if (include_t4_ || include_u4_)
            quadDim_ = oa_ * oa_ * oa_ * oa_ * va_ * va_ * va_ * va_
                       + ob_ * ob_ * ob_ * ob_ * vb_ * vb_ * vb_ * vb_
                       + oa_ * oa_ * oa_ * ob_ * va_ * va_ * va_ * vb_
                       + oa_ * oa_ * ob_ * ob_ * va_ * va_ * vb_ * vb_
                       + oa_ * ob_ * ob_ * ob_ * va_ * vb_ * vb_ * vb_;

        /// evaluate the total number of amplitudes
        size_t ccamps_dim_ = singleDim_ + doubleDim_; // t1, t2
        if (include_t3_) ccamps_dim_ += tripleDim_; // t3
        if (include_t4_) ccamps_dim_ += quadDim_; // t4
        if (include_u0_) ccamps_dim_++; // u0
        if (include_u1_) ccamps_dim_ += singleDim_; // u1
        if (include_u2_) ccamps_dim_ += doubleDim_; // u2
        if (include_u3_) ccamps_dim_ += tripleDim_; // u3
        if (include_u4_) ccamps_dim_ += quadDim_; // u4

        /// print included amplitudes
        Printf("  Included amplitudes:\n");
        Printf(include_u1_ ? "    U0\n" : "");
        Printf("    %2s", "T1");
        Printf("    %2s\n", include_u1_ ? "U1" : "  ");
        Printf("    %2s", "T2");
        Printf("    %2s\n", include_u2_ ? "U2" : "  ");
        if (include_t3_ || include_u3_) {
            Printf("    %2s", include_t3_ ? "T3" : "  ");
            Printf("    %2s\n", include_u3_ ? "U3" : "  ");
        }
        if (include_t4_ || include_u4_) {
            Printf("    %2s", include_t4_ ? "T4" : "  ");
            Printf("    %2s\n", include_u4_ ? "U4" : "  ");
        }
        Printf("\n  Dimension of CC amplitudes: %d\n\n", ccamps_dim_);
        if (has_photon_){
            // calculate cavity volumes: lambda = (1/(4*pi)*V_cav)^(-1/2) -> V_cav = (4*pi/lambda^2)
            double lam_norm = 0;
            for (double lam : lambda_)
                lam_norm += lam * lam;

            double V_cav = 4.0 * M_PI / lam_norm;

            // convert from bohr^3 to nm^3 (https://physics.nist.gov/cgi-bin/cuu/Value?Abohrrada0)
            double bohr_to_nm = 0.052917721090380;
            V_cav *= bohr_to_nm * bohr_to_nm * bohr_to_nm;

            // print cavity parameters
            Printf("  V_cav (nm3) = %6.4lf\n", V_cav);
            Printf("  lambda      = x: %6.4lf y: %6.4lf z: %6.4lf\n", lambda_[0], lambda_[1], lambda_[2]);
            Printf("  hw          = x: %6.4lf y: %6.4lf z: %6.4lf\n",
                   cavity_frequency_[0], cavity_frequency_[1], cavity_frequency_[2]);
            Printf(
                    options_.get_bool("QED_USE_RELAXED_ORBITALS") ? "  Using relaxed orbitals.\n"
                                                                  : "  Using unrelaxed orbitals.\n"
            );

        }
    }

    double CC_Cavity::compute_energy() {

        // grab some input options_
        double e_convergence = options_.get_double("E_CONVERGENCE");
        double r_convergence = options_.get_double("R_CONVERGENCE");
        size_t maxiter          = options_.get_int("MAXITER");

        Printf("\n");
        Printf("    madness threads:                %5i\n", options_.get_int("MAD_NUM_THREADS"));
        Printf("    omp threads:                    %5i\n", Process::environment.get_n_threads());
        Printf("    No. basis functions:            %5i\n",nso_);
        Printf("    No. auxiliary basis functions:  %5i\n",nQ_);
        Printf("    No. alpha electrons:            %5i\n",nalpha_);
        Printf("    No. alpha virtuals:             %5i\n",va_);
        Printf("    No. beta electrons:             %5i\n",nbeta_);
        Printf("    No. beta virtuals:              %5i\n",vb_);
        Printf("    e_convergence:             %10.3le\n",e_convergence);
        Printf("    r_convergence:             %10.3le\n",r_convergence);
        Printf("    maxiter:                        %5i\n",maxiter);
        Printf("\n\n");

        // print out the dimensions of the amplitudes we will be using
        print_dimensions();

        // make containers for amplitudes from derived classes
        init_operators();

        // print header for iterations
        print_iter_header();

        t_ground.start();
        cc_energy_ = 0.0;
        cc_energy_ = cc_iterations();
        world_.gop.fence();
        t_ground.stop();

        Printf("\n");
        Printf("    %s iterations converged!\n", cc_type_.c_str());
        Printf("\n");
        Printf("    * %s total energy: %20.12lf\n", cc_type_.c_str(), cc_energy_);

        print_properties();
        Printf(
                   "  ==> Time spent in CC iterations <==  \n"
                   "                    Total: %s\n"
                   "      Residuals Equations: %s\n"
                   "        Update Amplitudes: %s\n"
                   "      Transform Integrals: %s\n",
                   t_ground.elapsed().c_str(),
                   t_resid.elapsed().c_str(),
                   t_ampUp.elapsed().c_str(),
                   t_transform.elapsed().c_str()
               );

        Printf("\n");

        Process::environment.globals["CC_CAVITY TOTAL ENERGY"] = cc_energy_;
        Process::environment.globals["CURRENT ENERGY"]     = cc_energy_;

        return cc_energy_;
    }

    double CC_Cavity::cc_iterations() {

        // grab some input options_
        double e_convergence = options_.get_double("E_CONVERGENCE");
        double r_convergence = options_.get_double("R_CONVERGENCE");
        size_t maxiter          = options_.get_int("MAXITER");

        double e_last  = 0.0;
        double dele    = 0.0;
        double tnorm   = 0.0;

        double energy = 0.0;
        size_t iter = 0;

        bool r_converged = false;

        do {
            e_last = energy; // save old energy

            t_resid.start();
            energy = build_residuals(); world_.gop.fence(); // build residuals and return total energy
            t_resid.stop();

            t_ampUp.start();
            extrapolate_amplitudes(); world_.gop.fence(); // update and extrapolate amplitudes
            tnorm = compute_residual_norms(); world_.gop.fence(); // compute residual norms
            t_ampUp.stop();

            t_transform.start();
            transform_integrals(true); world_.gop.fence(); // t1-transformation e(-T1) H e(T1)
            t_transform.stop();

            dele = (iter != 0) ? energy - e_last : 0.0; // calculate energy change (if not first iteration)

            print_iteration(iter, energy, dele, tnorm); // print iteration info


            r_converged = tnorm < r_convergence; // check residual convergence
            for (auto &resid: residuals_){ // check individual residual convergences
                if (!r_converged) break;
                if(!resid.second.is_initialized()) continue;
                if (sqrt(norm2(resid.second)) > r_convergence)
                    r_converged = false;
            }

            if (iter++ >= maxiter) break; // limit iterations
        } while (
                    fabs(dele) > e_convergence || // ensure energy convergence
                    tnorm > r_convergence // ensure residual convergence
                ); // end loop when all conditions are met

        if (iter >= maxiter) {
            throw std::runtime_error("CC iterations did not converge!");
        }

        return energy;
    }

    void CC_Cavity::print_iter_header() const {
        Printf("\n");
        Printf("    ==>  Begin %s iterations <==    \n", cc_type_.c_str());
        Printf("\n");
        Printf("%5s %16s %15s %15s  | %8s %8s",  "Iter","energy","dE","|dT|","|dT1|","|dT2|");
        Printf("\n");
    }

    void CC_Cavity::print_iteration(size_t iter, double energy, double dele, double tnorm) const {
        Printf("%5i %17.12lf %15.12lf %15.12lf | %-8.1e %-8.1e",iter,energy,dele,tnorm,resid_norms_.at("t1"),resid_norms_.at("t2"));
        Printf("\n");
    }

    void CC_Cavity::apply_transform(TArrayMap &CL, TArrayMap &CR) {
        /// transform integrals to new basis
        t_oei.start();
        build_oei(CL, CR); world_.gop.fence(); // build one-electron integrals in MO basis
        t_oei.stop();

        t_tei.start();
        build_tei(CL, CR); world_.gop.fence(); // build two-electron integrals to MO basis and make them antisymmetric
        t_tei.stop();

        /// update orbital energies
        t_ampUp.start();
        build_eps(); world_.gop.fence();
        t_ampUp.stop();

        world_.gop.fence();
    }

    void CC_Cavity::build_oei(TArrayMap &CL, TArrayMap &CR) {

        /// build fock matrix in AO basis

        TArrayD F = core_H_.clone();
        if (has_photon_) {
            F("mu,nu") += oe_cavity_terms_("mu,nu");
        }

        // build fock matrix blocks in MO basis
        F_blks_.clear();
        for (auto& l_blk : CL) {
            for (auto& r_blk : CR) {
                string blk = l_blk.first.substr(0, 1) + r_blk.first.substr(0, 1); // aa, ab, ba, bb
                blk += "_" + l_blk.first.substr(2, 1) + r_blk.first.substr(2, 1); // oo, ov, vo, vv
                F_blks_[blk]("p,q") = l_blk.second("mu,p") * F("mu,nu") * r_blk.second("nu,q");
            }
        }

        /// build dipole integrals

        if (has_photon_ || true) { // this does not work for ccsd yet, so always do it
            Dip_blks_.clear();
            size_t nso = nso_;
            for (int i = 0; i < 3; ++i) {
                double **dip_p = dipole_[i]->pointer(); // dx, dy, dz

                // get polarization
                string pol;
                if (i == 0) pol = "dx_";
                else if (i == 1) pol = "dy_";
                else if (i == 2) pol = "dz_";

                // build dipole matricies in AO basis
                TArrayD dip = makeTensor(world_, {ns_, ns_}, false);
                dip.init_elements([dip_p, nso](auto &I) {
                    if (I[0] < nso && I[1] < nso)
                        return dip_p[I[0]][I[1]];
                    else if (I[0] >= nso && I[1] >= nso)
                        return dip_p[I[0] - nso][I[1] - nso];
                    else return 0.0;
                });
                world_.gop.fence();

                // build blocks of dipole matricies in MO basis
                for (auto &l_blk: CL) {
                    for (auto &r_blk: CR) {
                        string blk = l_blk.first.substr(0, 1) + r_blk.first.substr(0, 1); // aa, ab, ba, bb
                        blk += "_" + l_blk.first.substr(2, 1) + r_blk.first.substr(2, 1); // oo, ov, vo, vv
                        Dip_blks_[pol + blk]("p,q") = l_blk.second("mu,p") * dip("mu,nu") * r_blk.second("nu,q");
                    }
                }
            }
        }
    }

    void CC_Cavity::build_tei(TArrayMap &CL, TArrayMap &CR) {

        /// transform Qso to MO basis

        TArrayMap Qmo_blks;
        for (auto& l_blk : CL) {
            for (auto& r_blk : CR) {
                string blk = l_blk.first.substr(0, 1) + r_blk.first.substr(0, 1); // aa, ab, ba, bb
                blk += "_" + l_blk.first.substr(2, 1) + r_blk.first.substr(2, 1); // oo, ov, vo, vv
                Qmo_blks[blk]("Q,p,q") = l_blk.second("mu,p") * Qso_("Q,mu,nu") * r_blk.second("nu,q");
            }
        }

        /// unpack Qmo into 4-index containers
        unpack_eris(Qmo_blks);
//        world_.gop.fence();

        /// build two-electron components of the Fock matrix

        // exchange:
        // sum k (q|rk) (q|ks)
        F_blks_["aa_oo"]("i, j") -= Qmo_blks["aa_oo"]("Q, i, k") * Qmo_blks["aa_oo"]("Q, k, j");
        F_blks_["aa_oo"]("i, j") -= Qmo_blks["ab_oo"]("Q, i, k") * Qmo_blks["ba_oo"]("Q, k, j");
        F_blks_["bb_oo"]("i, j") -= Qmo_blks["ba_oo"]("Q, i, k") * Qmo_blks["ab_oo"]("Q, k, j");
        F_blks_["bb_oo"]("i, j") -= Qmo_blks["bb_oo"]("Q, i, k") * Qmo_blks["bb_oo"]("Q, k, j");

        F_blks_["aa_ov"]("i, a") -= Qmo_blks["aa_oo"]("Q, i, k") * Qmo_blks["aa_ov"]("Q, k, a");
        F_blks_["aa_ov"]("i, a") -= Qmo_blks["ab_oo"]("Q, i, k") * Qmo_blks["ba_ov"]("Q, k, a");
        F_blks_["bb_ov"]("i, a") -= Qmo_blks["ba_oo"]("Q, i, k") * Qmo_blks["ab_ov"]("Q, k, a");
        F_blks_["bb_ov"]("i, a") -= Qmo_blks["bb_oo"]("Q, i, k") * Qmo_blks["bb_ov"]("Q, k, a");

        F_blks_["aa_vo"]("a, i") -= Qmo_blks["aa_vo"]("Q, a, k") * Qmo_blks["aa_oo"]("Q, k, i");
        F_blks_["aa_vo"]("a, i") -= Qmo_blks["ab_vo"]("Q, a, k") * Qmo_blks["ba_oo"]("Q, k, i");
        F_blks_["bb_vo"]("a, i") -= Qmo_blks["ba_vo"]("Q, a, k") * Qmo_blks["ab_oo"]("Q, k, i");
        F_blks_["bb_vo"]("a, i") -= Qmo_blks["bb_vo"]("Q, a, k") * Qmo_blks["bb_oo"]("Q, k, i");

        F_blks_["aa_vv"]("a, b") -= Qmo_blks["aa_vo"]("Q, a, k") * Qmo_blks["aa_ov"]("Q, k, b");
        F_blks_["aa_vv"]("a, b") -= Qmo_blks["ab_vo"]("Q, a, k") * Qmo_blks["ba_ov"]("Q, k, b");
        F_blks_["bb_vv"]("a, b") -= Qmo_blks["ba_vo"]("Q, a, k") * Qmo_blks["ab_ov"]("Q, k, b");
        F_blks_["bb_vv"]("a, b") -= Qmo_blks["bb_vo"]("Q, a, k") * Qmo_blks["bb_ov"]("Q, k, b");

        // coulomb
        // sum k (q|kk) (q|rs)
        TA::TArrayD Qkk;
        Qkk("Q")  = Qmo_blks["aa_oo"]("Q, k, l") * Id_blks_["aa_oo"]("k, l");
        Qkk("Q") += Qmo_blks["bb_oo"]("Q, k, l") * Id_blks_["bb_oo"]("k, l");

        F_blks_["aa_oo"]("i, j") += Qkk("Q") * Qmo_blks["aa_oo"]("Q, i, j");
        F_blks_["bb_oo"]("i, j") += Qkk("Q") * Qmo_blks["bb_oo"]("Q, i, j");

        F_blks_["aa_ov"]("i, a") += Qkk("Q") * Qmo_blks["aa_ov"]("Q, i, a");
        F_blks_["bb_ov"]("i, a") += Qkk("Q") * Qmo_blks["bb_ov"]("Q, i, a");

        F_blks_["aa_vo"]("a, i") += Qkk("Q") * Qmo_blks["aa_vo"]("Q, a, i");
        F_blks_["bb_vo"]("a, i") += Qkk("Q") * Qmo_blks["bb_vo"]("Q, a, i");

        F_blks_["aa_vv"]("a, b") += Qkk("Q") * Qmo_blks["aa_vv"]("Q, a, b");
        F_blks_["bb_vv"]("a, b") += Qkk("Q") * Qmo_blks["bb_vv"]("Q, a, b");

        // add dipole self energy contribution to Fock matrices
        if (has_photon_) {
            for (int i = 0; i < 3; ++i) {
                string pol; // polarization
                if (i == 0) pol = "dx_";
                else if (i == 1) pol = "dy_";
                else if (i == 2) pol = "dz_";

                double dp_oo = Dip_blks_[pol + "aa_oo"]("i, j").trace() + Dip_blks_[pol + "bb_oo"]("i, j").trace();

                // oo blocks
                F_blks_["aa_oo"]("i, j") += lambda_[i] * lambda_[i] * Dip_blks_[pol + "aa_oo"]("i, j") * dp_oo;
                F_blks_["aa_oo"]("i, j") -=
                        lambda_[i] * lambda_[i] * Dip_blks_[pol + "aa_oo"]("i, m") * Dip_blks_[pol + "aa_oo"]("m, j");

                F_blks_["bb_oo"]("i, j") += lambda_[i] * lambda_[i] * Dip_blks_[pol + "bb_oo"]("i, j") * dp_oo;
                F_blks_["bb_oo"]("i, j") -=
                        lambda_[i] * lambda_[i] * Dip_blks_[pol + "bb_oo"]("i, m") * Dip_blks_[pol + "bb_oo"]("m, j");

                // ov blocks
                F_blks_["aa_ov"]("i, a") += lambda_[i] * lambda_[i] * Dip_blks_[pol + "aa_ov"]("i, a") * dp_oo;
                F_blks_["aa_ov"]("i, a") -=
                        lambda_[i] * lambda_[i] * Dip_blks_[pol + "aa_oo"]("i, m") * Dip_blks_[pol + "aa_ov"]("m, a");

                F_blks_["bb_ov"]("i, a") += lambda_[i] * lambda_[i] * Dip_blks_[pol + "bb_ov"]("i, a") * dp_oo;
                F_blks_["bb_ov"]("i, a") -=
                        lambda_[i] * lambda_[i] * Dip_blks_[pol + "bb_oo"]("i, m") * Dip_blks_[pol + "bb_ov"]("m, a");

                // vo blocks
                F_blks_["aa_vo"]("a, i") += lambda_[i] * lambda_[i] * Dip_blks_[pol + "aa_vo"]("a, i") * dp_oo;
                F_blks_["aa_vo"]("a, i") -=
                        lambda_[i] * lambda_[i] * Dip_blks_[pol + "aa_vo"]("a, m") * Dip_blks_[pol + "aa_oo"]("m, i");

                F_blks_["bb_vo"]("a, i") += lambda_[i] * lambda_[i] * Dip_blks_[pol + "bb_vo"]("a, i") * dp_oo;
                F_blks_["bb_vo"]("a, i") -=
                        lambda_[i] * lambda_[i] * Dip_blks_[pol + "bb_vo"]("a, m") * Dip_blks_[pol + "bb_oo"]("m, i");

                // vv blocks
                F_blks_["aa_vv"]("a, b") += lambda_[i] * lambda_[i] * Dip_blks_[pol + "aa_vv"]("a, b") * dp_oo;
                F_blks_["aa_vv"]("a, b") -=
                        lambda_[i] * lambda_[i] * Dip_blks_[pol + "aa_vo"]("a, m") * Dip_blks_[pol + "aa_ov"]("m, b");

                F_blks_["bb_vv"]("a, b") += lambda_[i] * lambda_[i] * Dip_blks_[pol + "bb_vv"]("a, b") * dp_oo;
                F_blks_["bb_vv"]("a, b") -=
                        lambda_[i] * lambda_[i] * Dip_blks_[pol + "bb_vo"]("a, m") * Dip_blks_[pol + "bb_ov"]("m, b");
            }
        }
    }

    void CC_Cavity::unpack_eris(TArrayMap &Qmo_blks) {

        double lx2, ly2, lz2;
        if (has_photon_) {
            lx2 = lambda_[0] * lambda_[0];
            ly2 = lambda_[1] * lambda_[1];
            lz2 = lambda_[2] * lambda_[2];
        }

        static unordered_set<string> valid_names = {
            "aaaa_oooo", "aaaa_vvvv", "aaaa_oovv", "aaaa_vvoo", "aaaa_vovo", "aaaa_vooo", "aaaa_oovo", "aaaa_vovv", "aaaa_vvvo",
            "abab_oooo", "abab_vvvv", "abab_oovv", "abab_vvoo", "abab_vovo", "abab_vooo", "abab_oovo", "abab_vovv", "abab_vvvo",
            "bbbb_oooo", "bbbb_vvvv", "bbbb_oovv", "bbbb_vvoo", "bbbb_vovo", "bbbb_vooo", "bbbb_oovo", "bbbb_vovv", "bbbb_vvvo",

            // extra terms generated from parsed equations (helps optimization in TABuilder)
            "abba_oovo", "abba_vooo", "baab_vooo", "baba_vovo", "abba_vovo", "baab_vovo", "baab_vovv", "abba_vvvo"
        };

        // build string for eri spin
        static auto build_eri_spin = [](const string& colm_lspin, const string& colm_rspin) {
            string eri_spin(4, ' ');
            eri_spin[0] = colm_lspin[0]; eri_spin[1] = colm_rspin[0];
            eri_spin[2] = colm_lspin[1]; eri_spin[3] = colm_rspin[1];
            return eri_spin;
        };

        // build string for eri occupation
        static auto build_eri_ov = [](const string& colm_lov, const string& colm_rov) {
            string eri_ov(4, ' ');
            eri_ov[0] = colm_lov[0]; eri_ov[1] = colm_rov[0];
            eri_ov[2] = colm_lov[1]; eri_ov[3] = colm_rov[1];
            return eri_ov;
        };

        // check if eri spin is valid
        static auto is_valid_spin = [](const string& eri_spin) {
            size_t num_alpha = count(eri_spin.begin(), eri_spin.end(), 'a');
            size_t num_beta  = count(eri_spin.begin(), eri_spin.end(), 'b');
            return num_alpha % 2 == 0 && num_beta % 2 == 0;
        };

        // check if eri is in V_blks_
        auto &V_blks = V_blks_;
        auto is_in_V_blks = [&V_blks](const string& eri_name) {
            return V_blks.find(eri_name) != V_blks.end();
        };

        // check if eri name is valid
        static auto is_valid_name = [](const string& eri_name, const unordered_set<string>& ) {
            return valid_names.find(eri_name) != valid_names.end();
        };

        // build string for exchange name
        static auto build_exchange_name = [](const string& lspin, const string& rspin,
                                             const string& lov, const string& rov) {
            string exc_spin(2, ' ');
            string exc_ov(2, ' ');

            // Build spin and ov strings
            exc_spin[0] = lspin[0]; exc_spin[1] = rspin[1];
            exc_ov[0] =   lov[0];   exc_ov[1] =   rov[1];

            // Combine spin and orbital occupation strings
            string exchange_name = exc_spin + "_" + exc_ov;
            return exchange_name.substr(0, 5);
        };

        // remove eris from last build (if any)
        V_blks_.clear();

        // loop over all Qmo blocks and build electron repulsion integrals
        for (auto& Ql : Qmo_blks) {
            string colm_l = Ql.first;
            string colm_lspin = colm_l.substr(0, 2);
            string colm_lov = colm_l.substr(3, 2);

            for (auto& Qr : Qmo_blks) {
                string colm_r = Qr.first;
                string colm_rspin = colm_r.substr(0, 2);
                string colm_rov = colm_r.substr(3, 2);

                string eri_spin = build_eri_spin(colm_lspin, colm_rspin);
                string eri_ov = build_eri_ov(colm_lov, colm_rov);
                string eri_name = eri_spin + "_" + eri_ov;

                if (!is_valid_spin(eri_spin) || is_in_V_blks(eri_name) || !is_valid_name(eri_name, valid_names))
                    continue;

                string exc_l = build_exchange_name(colm_lspin, colm_rspin, colm_lov, colm_rov);
                string exc_r = build_exchange_name(colm_rspin, colm_lspin, colm_rov, colm_lov);

                TArrayD ijab, ibaj;
                ijab("i,j,a,b") = Ql.second("Q,i,j") * Qr.second("Q,a,b");
                ibaj("i,b,a,j") = Qmo_blks[exc_l]("Q,i,b") * Qmo_blks[exc_r]("Q,a,j");

                /// add dipole terms
                if (has_photon_) {
                    // add dipole terms for coulomb
                    ijab("i,j,a,b") += lx2 * Dip_blks_["dx_" + colm_l]("i,j") * Dip_blks_["dx_" + colm_r]("a,b");
                    ijab("i,j,a,b") += ly2 * Dip_blks_["dy_" + colm_l]("i,j") * Dip_blks_["dy_" + colm_r]("a,b");
                    ijab("i,j,a,b") += lz2 * Dip_blks_["dz_" + colm_l]("i,j") * Dip_blks_["dz_" + colm_r]("a,b");

                    // add dipole terms for exchange
                    ibaj("i,b,a,j") += lx2 * Dip_blks_["dx_"+exc_l]("i,b") * Dip_blks_["dx_"+exc_r]("a,j");
                    ibaj("i,b,a,j") += ly2 * Dip_blks_["dy_"+exc_l]("i,b") * Dip_blks_["dy_"+exc_r]("a,j");
                    ibaj("i,b,a,j") += lz2 * Dip_blks_["dz_"+exc_l]("i,b") * Dip_blks_["dz_"+exc_r]("a,j");
                }

                V_blks_[eri_name]("i,a,j,b") = ijab("i,j,a,b") - ibaj("i,b,a,j");
            }
        }
    }

    void CC_Cavity::build_eps() {

        double* eps = epsilon_;
        size_t oa = oa_, o = o_, va = va_;

        memset(eps, 0, sizeof(double)*ns_);

        // oo/aa block
        forall(F_blks_["aa_oo"], [eps,oa,o,va](auto &tile, auto &x){
            if (x[0] != x[1]) return;
            eps[x[0]] = tile[x];
        });

        // oo/bb block
        forall(F_blks_["bb_oo"], [eps,oa,o,va](auto &tile, auto &x){
            if (x[0] != x[1]) return;
            eps[x[0] + oa] = tile[x];
        });

        // vv/aa block
        forall(F_blks_["aa_vv"], [eps,oa,o,va](auto &tile, auto &x){
            if (x[0] != x[1]) return;
            eps[x[0] + o] = tile[x];
        });

        // vv/bb block
        forall(F_blks_["bb_vv"], [eps,oa,o,va](auto &tile, auto &x){
            if (x[0] != x[1]) return;
            eps[x[0] + o + va] = tile[x];
        });

        world_.gop.fence();
        world_.gop.reduce(eps, ns_, std::plus<>());

    }


    void CC_Cavity::update_amplitudes(){
        double *eps = epsilon_;
        size_t o = o_;
        size_t v = v_;
        size_t oa = oa_;
        size_t va = va_;

        /// dt = -residual / eps

        // t1
        forall(residuals_["t1_aa"],
                        [eps, o, oa, va](auto &tile, auto &x) {
                            tile[x] /= (eps[x[1]] - eps[x[0]+o]);
                        });

        forall(residuals_["t1_bb"],
                        [eps, o, oa, va](auto &tile, auto &x) {
                            tile[x] /= (eps[x[1]+oa] - eps[x[0]+o+va]);
                        });

        // t2
        forall(residuals_["t2_aaaa"],
                        [eps, o, oa, va](auto &tile, auto &x) {
                            double o_ep = eps[x[2]] + eps[x[3]],
                                    v_ep = eps[x[0]+o] + eps[x[1]+o];
                            tile[x] /= (o_ep - v_ep);
                        });
        forall(residuals_["t2_bbbb"],
                        [eps, o, oa, va](auto &tile, auto &x) {
                            double o_ep = eps[x[2]+oa] + eps[x[3]+oa],
                                    v_ep = eps[x[0]+o+va] + eps[x[1]+o+va];
                            tile[x] /= (o_ep - v_ep);
                        });
        forall(residuals_["t2_abab"],
                        [eps, o, oa, va](auto &tile, auto &x) {
                            double o_ep = eps[x[2]] + eps[x[3]+oa],
                                    v_ep = eps[x[0]+o] + eps[x[1]+o];
                            tile[x] /= (o_ep - v_ep);
                        });

        world_.gop.fence();

        /// update amplitudes according to t + dt = amplitude - residual / eps

        amplitudes_["t1_aa"](idx_map_[2]) += residuals_["t1_aa"](idx_map_[2]);
        amplitudes_["t1_bb"](idx_map_[2]) += residuals_["t1_bb"](idx_map_[2]);
        amplitudes_["t2_aaaa"](idx_map_[4]) += residuals_["t2_aaaa"](idx_map_[4]);
        amplitudes_["t2_abab"](idx_map_[4]) += residuals_["t2_abab"](idx_map_[4]);
        amplitudes_["t2_bbbb"](idx_map_[4]) += residuals_["t2_bbbb"](idx_map_[4]);

        world_.gop.fence();
    }

    void CC_Cavity::extrapolate_amplitudes() {

        // update amplitudes
        update_amplitudes();

        /// build vectors for DIIS

        // amplitudes
        std::vector<TA::TArrayD*> amp_vec, resid_vec;
        for (auto& amp : amplitudes_) {
            if(!amp.second.is_initialized()) continue;
            amp_vec.push_back(&(amp.second));
        }

        // residuals
        for (auto& resid : residuals_) {
            if(!resid.second.is_initialized()) continue;
            resid_vec.push_back(&(resid.second));
        }

        /// Perform DIIS extrapolation
        world_.gop.fence();
        diis_ta->WriteVector(amp_vec);
        diis_ta->WriteErrorVector(resid_vec);
        diis_ta->Extrapolate(amp_vec);
        world_.gop.fence();

        amp_vec.clear(); resid_vec.clear(); // clear vectors
    }

    double CC_Cavity::compute_residual_norms(bool return_tot) {

        resid_norms_["t1"] = sqrt(squared_norm(residuals_["t1_aa"])
                                  + squared_norm(residuals_["t1_bb"]));
        resid_norms_["t2"] = sqrt(squared_norm(residuals_["t2_aaaa"])
                                  + squared_norm(residuals_["t2_abab"])
                                  + squared_norm(residuals_["t2_bbbb"]));

        // total residual norm for all amplitudes (if requested)
        if (return_tot) {
            double norm = 0.0;
            for (auto & amp : residuals_) {
                if (!amp.second.is_initialized()) continue;
                norm += squared_norm(amp.second);
            }
            return sqrt(norm);
        }
        else return 0.0;
    }

    void CC_Cavity::print_properties() {
        // calculate norms of amplitudes
        map<string, double> amp_norms;
        double total_norm = 0.0;
        for (auto& amp : amplitudes_) {
            if(!amp.second.is_initialized()) continue;
            // name of amplitude is first two characters of key
            double amp_norm = squared_norm(amp.second);
            string amp_name = amp.first.substr(0,2);
            amp_norms[amp_name] += amp_norm;
            total_norm += amp_norm;
        }

        double nT1 = amp_norms["t1"];
        double nT2 = amp_norms["t2"];

        // print norms of cluster amplitudes
        outfile->Printf("\n\n   Norms of cluster amplitudes:");
        outfile->Printf("\n   ------------------------------");
        outfile->Printf("\n    T1: %15.12lf | %5.2f %%", sqrt(nT1), 100*nT1/total_norm);
        outfile->Printf("\n    T2: %15.12lf | %5.2f %%", sqrt(nT2), 100*nT2/total_norm);
        outfile->Printf("\n   ------------------------------");
        outfile->Printf("\n    Total: %15.12lf\n\n", sqrt(total_norm));

    }
}

