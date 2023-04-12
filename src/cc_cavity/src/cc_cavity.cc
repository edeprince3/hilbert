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
            tnorm = update_amplitudes(); world_.gop.fence(); // update amplitudes and return residual norm
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
}