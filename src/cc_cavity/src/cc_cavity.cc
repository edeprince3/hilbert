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

#include "../misc/ta_helper.h"
#include "../misc/threeindexintegralsta.h"

#include <mkl.h>
#include <omp.h>
#include "../misc/qed_blas.h"
#include "../../misc/hilbert_psifiles.h"
#include "../../polaritonic_scf/uhf.h"
#include <unistd.h>
#include <psi4/psifiles.h>
#include "../include/cc_cavity.h"

using namespace std;
using namespace TA;
using namespace psi;
using namespace Helper;

namespace hilbert {

    CC_Cavity::CC_Cavity(std::shared_ptr<Wavefunction> reference_wavefunction, Options& options_):
            PolaritonicHF(std::move(reference_wavefunction), options_) {

        // set block_size for tiledarray
        HelperD::tile_size_ = (size_t) options_.get_int("TILE_SIZE");

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

        // grab dimensions
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

    }

    void CC_Cavity::init_operators() {
        /// initialize arbitrary index strings
        idx_map_ = std::vector<std::string>(25);
        idx_map_[0] = "";
        for (int i = 1; i < 25; i++) { // hard-coded for now. Unlikely to need more than 25 indices.

            // build string of indices for i'th rank tensor
            idx_map_[i] = "";
            for (int j = 0; j < i; j++) idx_map_[i] += "x" + std::to_string(i) + ",";

            // remove last comma
            idx_map_[i].pop_back();
        }

        /// make identity tensors

        Id_blks_["a_o"] = HelperD::makeTensor(world_, {oa_}, false);   Id_blks_["a_o"].fill(1.0);
        Id_blks_["a_v"] = HelperD::makeTensor(world_, {va_}, false);   Id_blks_["a_v"].fill(1.0);
        Id_blks_["b_o"] = HelperD::makeTensor(world_, {ob_}, false);   Id_blks_["b_o"].fill(1.0);
        Id_blks_["b_v"] = HelperD::makeTensor(world_, {vb_}, false);   Id_blks_["b_v"].fill(1.0);

        Id_blks_["aa_oo"] = TiledArray::diagonal_array<TA::TArrayD>(world_, HelperD::makeRange({oa_, oa_}), 1);
        Id_blks_["aa_vv"] = TiledArray::diagonal_array<TA::TArrayD>(world_, HelperD::makeRange({va_, va_}), 1);
        Id_blks_["bb_oo"] = TiledArray::diagonal_array<TA::TArrayD>(world_, HelperD::makeRange({ob_, ob_}), 1);
        Id_blks_["bb_vv"] = TiledArray::diagonal_array<TA::TArrayD>(world_, HelperD::makeRange({vb_, vb_}), 1);

        Id_blks_["aaaa_oooo"]("p,q,r,s") = Id_blks_["aa_oo"]("p,r") * Id_blks_["aa_oo"]("q,s");
        Id_blks_["abab_oooo"]("p,q,r,s") = Id_blks_["aa_oo"]("p,r") * Id_blks_["bb_oo"]("q,s");
        Id_blks_["bbbb_oooo"]("p,q,r,s") = Id_blks_["bb_oo"]("p,r") * Id_blks_["bb_oo"]("q,s");

        /// initialize t1 and t2 amplitudes and residuals (every derived class must use these amplitudes and residuals)

        // t1
        amplitudes_["t1_aa"]  = HelperD::makeTensor(world_, {va_, oa_}, true);
        amplitudes_["t1_bb"]  = HelperD::makeTensor(world_, {vb_, ob_}, true);
        residuals_["t1_aa"]  = HelperD::makeTensor(world_, {va_, oa_}, true);
        residuals_["t1_bb"]  = HelperD::makeTensor(world_, {vb_, ob_}, true);

        // t2
        amplitudes_["t2_aaaa"]  = HelperD::makeTensor(world_, {va_,va_, oa_,oa_}, true);
        amplitudes_["t2_abab"]  = HelperD::makeTensor(world_, {va_,vb_, oa_,ob_}, true);
        amplitudes_["t2_bbbb"]  = HelperD::makeTensor(world_, {vb_,vb_, ob_,ob_}, true);
        residuals_["t2_aaaa"]  = HelperD::makeTensor(world_, {va_,va_, oa_,oa_}, true);
        residuals_["t2_abab"]  = HelperD::makeTensor(world_, {va_,vb_, oa_,ob_}, true);
        residuals_["t2_bbbb"]  = HelperD::makeTensor(world_, {vb_,vb_, ob_,ob_}, true);
    }

    void CC_Cavity::init_integrals() {

        /// build MO transformation matrix
        size_t ns = ns_, oa = oa_, ob = ob_;

        // grab the MO coefficients
        double ** ca = Ca_->pointer();
        double ** cb = Cb_->pointer();

        C_blks_["a_o"] = HelperD::makeTensor(world_, {ns_, oa_}, false);
        C_blks_["b_o"] = HelperD::makeTensor(world_, {ns_, ob_}, false);
        C_blks_["a_v"] = HelperD::makeTensor(world_, {ns_, va_}, false);
        C_blks_["b_v"] = HelperD::makeTensor(world_, {ns_, vb_}, false);

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
        core_H_ = HelperD::makeTensor(world_, {ns_,ns_}, false);
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

            oe_cavity_terms_ = HelperD::makeTensor(world_, {ns_, ns_}, false);
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
        /// TODO: This is a hack; should be done in a memory efficient way (without breaking things)
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
            auto * tmp_so_p = (double*)calloc(nQ_ * ns_ * ns_, sizeof(double));
            for (size_t Q = 0; Q < nQ_; ++Q) {
                for (size_t mu = 0; mu < nso_; ++mu) {
                    for (size_t nu = 0; nu < nso_; ++nu) {
                        tmp_so_p[Q*ns_*ns_ + mu * ns_ + nu] = Qso_p[Q][mu*nso_ + nu];
                        tmp_so_p[Q*ns_*ns_ + (mu + nso_) * ns_ + (nu + nso_)] = Qso_p[Q][mu*nso_ + nu];
                    }
                }
            }
            world_.gop.fence();
            Qso.reset();
            Qso_ = HelperD::makeTensor(world_, {nQ_, ns_, ns_}, tmp_so_p);
            free(tmp_so_p);
        } else if ( options_.get_str("SCF_TYPE") == "CD" ) {
            double* Qso_p = ThreeIndexIntegrals(reference_wavefunction_,nQ_,memory_);
            auto * tmp_so_p = (double*)calloc(nQ_ * ns_ * ns_, sizeof(double));
            for (size_t Q = 0; Q < nQ_; ++Q) {
                for (size_t mu = 0; mu < nso_; ++mu) {
                    for (size_t nu = 0; nu < nso_; ++nu) {
                        tmp_so_p[Q*ns_*ns_ + mu * ns_ + nu] = Qso_p[Q * nso_ * nso_ + mu * nso_ + nu];
                        tmp_so_p[Q*ns_*ns_ + (mu + nso_) * ns_ + (nu + nso_)] = Qso_p[Q * nso_ * nso_ + mu * nso_ + nu];
                    }
                }
            }
            world_.gop.fence();
            free(Qso_p);
            Qso_ = HelperD::makeTensor(world_, {nQ_, ns_, ns_}, tmp_so_p);
            free(tmp_so_p);
        }
        world_.gop.fence();

        print_dimensions();

    }

    void CC_Cavity::transform_integrals(bool use_t1) {
        has_t1_integrals_ = use_t1;

        if (use_t1) {

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
        } else
            apply_transform(C_blks_, C_blks_);
    }

    double CC_Cavity::compute_energy() {

        // grab some input options_
        double e_convergence = options_.get_double("E_CONVERGENCE");
        double r_convergence = options_.get_double("R_CONVERGENCE");
        size_t maxiter          = options_.get_int("MAXITER");

        Printf("\n");
        Printf("    No. basis functions:            %5i\n",nso_);
        Printf("    No. auxiliary basis functions:  %5i\n",nQ_);
        Printf("    No. alpha electrons:            %5i\n",nalpha_);
        Printf("    No. alpha virtuals:             %5i\n",va_);
        Printf("    No. beta electrons:             %5i\n",nbeta_);
        Printf("    No. beta virtuals:              %5i\n",vb_);
        Printf("    e_convergence:             %10.3le\n",e_convergence);
        Printf("    r_convergence:             %10.3le\n",r_convergence);
        Printf("    maxiter:                        %5i\n",maxiter);
        Printf("\n");

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

        Process::environment.globals["UCCSD TOTAL ENERGY"] = cc_energy_;
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
            for (auto &amp: amplitudes_){ // check amplitudes convergence
                if (!r_converged) break;
                if(!amp.second.is_initialized()) continue;
                if (sqrt(norm2(amp.second)) > r_convergence) r_converged = true;
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

        if (has_photon_) {
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
                TArrayD dip = HelperD::makeTensor(world_, {ns_, ns_}, false);
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
        auto is_in_V_blks = [this](const string& eri_name) {
            return this->V_blks_.find(eri_name) != this->V_blks_.end();
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
        HelperD::forall(F_blks_["aa_oo"], [eps,oa,o,va](auto &tile, auto &x){
            if (x[0] != x[1]) return;
            eps[x[0]] = tile[x];
        });

        // oo/bb block
        HelperD::forall(F_blks_["bb_oo"], [eps,oa,o,va](auto &tile, auto &x){
            if (x[0] != x[1]) return;
            eps[x[0] + oa] = tile[x];
        });

        // vv/aa block
        HelperD::forall(F_blks_["aa_vv"], [eps,oa,o,va](auto &tile, auto &x){
            if (x[0] != x[1]) return;
            eps[x[0] + o] = tile[x];
        });

        // vv/bb block
        HelperD::forall(F_blks_["bb_vv"], [eps,oa,o,va](auto &tile, auto &x){
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
        HelperD::forall(residuals_["t1_aa"],
                        [eps, o, oa, va](auto &tile, auto &x) {
                            tile[x] /= (eps[x[1]] - eps[x[0]+o]);
                        });

        HelperD::forall(residuals_["t1_bb"],
                        [eps, o, oa, va](auto &tile, auto &x) {
                            tile[x] /= (eps[x[1]+oa] - eps[x[0]+o+va]);
                        });

        // t2
        HelperD::forall(residuals_["t2_aaaa"],
                        [eps, o, oa, va](auto &tile, auto &x) {
                            double o_ep = eps[x[2]] + eps[x[3]],
                                    v_ep = eps[x[0]+o] + eps[x[1]+o];
                            tile[x] /= (o_ep - v_ep);
                        });
        HelperD::forall(residuals_["t2_bbbb"],
                        [eps, o, oa, va](auto &tile, auto &x) {
                            double o_ep = eps[x[2]+oa] + eps[x[3]+oa],
                                    v_ep = eps[x[0]+o+va] + eps[x[1]+o+va];
                            tile[x] /= (o_ep - v_ep);
                        });
        HelperD::forall(residuals_["t2_abab"],
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

}