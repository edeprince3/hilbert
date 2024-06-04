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
#include "cc_cavity/include/lambda_driver.h"
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

#include "cc_cavity/misc/ta_helper.hpp"
#include "misc/threeindexintegrals.h"

#include <mkl.h>
#include <omp.h>
#include "misc/blas.h"
#include "misc/hilbert_psifiles.h"
#include "polaritonic_scf/uhf.h"
#include <unistd.h>
#include <psi4/psifiles.h>
using namespace TA_Helper;

namespace hilbert {

    void LambdaDriver::init_operators() {

        // copy all amplitudes from T amplitudes and conjugate into L amplitudes
        for (const auto &t_pair: T_amplitudes_) {
            auto [t_name, t_amp] = t_pair; // extract name of amplitude and tensor

            // copy name of L amplitude from T
            std::string l_name = t_name;
            if (l_name[0] == 't') l_name[0] = 'l';
            if (l_name[0] == 'u') l_name[0] = 'm';

            // get index for T amplitude
            size_t rank = t_amp.trange().rank();
            std::string t_idxs = cc_wfn_->idx_map_[rank];
            size_t half_pos = t_idxs.size()/2;

            // conjugate index for L amplitude
            std::string l_idxs = t_idxs.substr(t_idxs.back() - half_pos, half_pos);
            l_idxs += ",";
            l_idxs += t_idxs.substr(t_idxs.front(), half_pos);

            // initialize L amplitude from T amplitude
            L_amplitudes_["l_name"](l_idxs) = t_amp(t_idxs);

        }

    }

    double LambdaDriver::compute_lambda() {
        // grab some input options_
        double e_convergence = options_.get_double("E_CONVERGENCE");
        double r_convergence = options_.get_double("R_CONVERGENCE");
        size_t maxiter = options_.get_int("MAXITER");

        // make containers for amplitudes from derived classes
        init_operators();

        // print header for iterations
        print_iter_header();

        t_ground.start();
        lam_energy_ = 0;
        lam_energy_ = lambda_iterations();
        world_.gop.fence();
        t_ground.stop();

        Printf("\n");
        Printf("    %s iterations converged!\n", cc_type_.c_str());
        Printf("\n");
        Printf("    * %s total energy: %20.12lf\n", cc_type_.c_str(), lam_energy_);

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

        Process::environment.globals["CC_CAVITY LAMBDA ENERGY"] = lam_energy_;
        return lam_energy_;
    }

    double LambdaDriver::lambda_iterations() {
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


        // ensure integrals are t1 transformed
        if (!cc_wfn_->has_t1_integrals_) cc_wfn_->transform_integrals(true);
        cc_wfn_->diis_ta->restart(); // restart DIIS

        // build shared intermediates
        build_intermediates();

        do {
            e_last = energy; // save old energy

            t_resid.start();
            energy = build_residuals(); world_.gop.fence(); // build residuals and return total energy
            t_resid.stop();

            t_ampUp.start();
            extrapolate_amplitudes(); world_.gop.fence(); // update and extrapolate amplitudes
            tnorm = compute_residual_norms(); world_.gop.fence(); // compute residual norms
            t_ampUp.stop();

            dele = (iter != 0) ? energy - e_last : 0.0; // calculate energy change (if not first iteration)

            print_iteration(iter, energy, dele, tnorm); // print iteration info

            r_converged = tnorm < r_convergence; // check residual convergence
            for (auto &amp: L_residuals_){ // check individual residual convergences
                if (!r_converged) break;
                if(!amp.second.is_initialized()) continue;
                if (sqrt(norm2(amp.second)) > r_convergence)
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

    void LambdaDriver::update_amplitudes() {

        /// dt = -residual / eps

        double *eps = cc_wfn_->epsilon_;
        double w0 = cc_wfn_->cavity_frequency_[2];
        size_t o =  o_;
        size_t v =  v_;
        size_t oa = oa_;
        size_t va = va_;

        // l1
        {
            forall(L_residuals_["l1_aa"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                tile[x] /= (eps[x[1]] - eps[x[0] + o]);
                            });

            forall(L_residuals_["l1_bb"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                tile[x] /= (eps[x[1] + oa] - eps[x[0] + o + va]);
                            });
        }

        // m1
        if (include_u1_) {
            forall(L_residuals_["m1_aa"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                tile[x] /= (eps[x[1]] - eps[x[0] + o] + w0);
                            });

            forall(L_residuals_["m1_bb"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                tile[x] /= (eps[x[1] + oa] - eps[x[0] + o + va] + w0);
                            });
        }

        // l2
        {
            // l2
            forall(L_residuals_["l2_aaaa"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                double o_ep = eps[x[2]] + eps[x[3]],
                                       v_ep = eps[x[0] + o] + eps[x[1] + o];
                                tile[x] /= (o_ep - v_ep);
                            });
            forall(L_residuals_["l2_bbbb"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                double o_ep = eps[x[2] + oa] + eps[x[3] + oa],
                                       v_ep = eps[x[0] + o + va] + eps[x[1] + o + va];
                                tile[x] /= (o_ep - v_ep);
                            });
            forall(L_residuals_["l2_abab"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                double o_ep = eps[x[2]] + eps[x[3] + oa],
                                       v_ep = eps[x[0] + o] + eps[x[1] + o] + w0;
                                tile[x] /= (o_ep - v_ep);
                            });
        }

        // m2
        if (include_u2_) {
            // m2
            forall(L_residuals_["m2_aaaa"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                double o_ep = eps[x[2]] + eps[x[3]],
                                       v_ep = eps[x[0] + o] + eps[x[1] + o] + w0;
                                tile[x] /= (o_ep - v_ep);
                            });
            forall(L_residuals_["u2_bbbb"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                double o_ep = eps[x[2] + oa] + eps[x[3] + oa],
                                       v_ep = eps[x[0] + o + va] + eps[x[1] + o + va] + w0;
                                tile[x] /= (o_ep - v_ep);
                            });
            forall(L_residuals_["u2_abab"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                double o_ep = eps[x[2]] + eps[x[3] + oa],
                                       v_ep = eps[x[0] + o] + eps[x[1] + o] + w0;
                                tile[x] /= (o_ep - v_ep);
                            });
        }


        world_.gop.fence();
        double rm0_ = L_scalar_resids_["m0"];
        double m0_ = L_scalar_resids_["m0"];

        // m0
        if ( include_u0_ ) {
            rm0_ /= -w0;
            L_amplitudes_["m0"] = makeTensor(world_, {1}, &m0_);
            L_residuals_["m0"] = makeTensor(world_, {1}, &rm0_);
        }

        /// update amplitudes according to u + du = amplitude - residual / (eps + w)
        if (include_u0_) L_amplitudes_["u0"](idx_map_[1]) -= L_residuals_["u0"](idx_map_[1]);
        if (include_u1_) {
            L_amplitudes_["u1_aa"](idx_map_[2]) -= L_residuals_["u1_aa"](idx_map_[2]);
            L_amplitudes_["u1_bb"](idx_map_[2]) -= L_residuals_["u1_bb"](idx_map_[2]);
        }
        if (include_u2_) {
            L_amplitudes_["u2_aaaa"](idx_map_[4]) -= L_residuals_["u2_aaaa"](idx_map_[4]);
            L_amplitudes_["u2_abab"](idx_map_[4]) -= L_residuals_["u2_abab"](idx_map_[4]);
            L_amplitudes_["u2_bbbb"](idx_map_[4]) -= L_residuals_["u2_bbbb"](idx_map_[4]);
        }

        world_.gop.fence();

    }


    void LambdaDriver::extrapolate_amplitudes() {
        // update amplitudes
        update_amplitudes();

        /// build vectors for DIIS

        // amplitudes
        std::vector<TA::TArrayD*> amp_vec, resid_vec;
        for (auto& amp : L_amplitudes_) {
            if(!amp.second.is_initialized()) continue;
            amp_vec.push_back(&(amp.second));
        }

        // residuals
        for (auto& resid : L_residuals_) {
            if(!resid.second.is_initialized()) continue;
            resid_vec.push_back(&(resid.second));
        }

        /// Perform DIIS extrapolation
        world_.gop.fence();
        cc_wfn_->diis_ta->WriteVector(amp_vec);
        cc_wfn_->diis_ta->WriteErrorVector(resid_vec);
        cc_wfn_->diis_ta->Extrapolate(amp_vec);
        world_.gop.fence();

        amp_vec.clear(); resid_vec.clear(); // clear vectors
    }

    double LambdaDriver::compute_residual_norms(bool return_tot) {

        /// residual norms for each amplitude

        // call parent function for t1 and t2

        // l residual norms
        if (include_u1_) {
            L_resid_norms_["l1"] = sqrt(squared_norm(L_residuals_["l1_aa"])
                                      + squared_norm(L_residuals_["l1_bb"]));
        }
        if (include_u2_) {
            L_resid_norms_["l2"] = sqrt(squared_norm(L_residuals_["l2_aaaa"])
                                      + squared_norm(L_residuals_["l2_abab"])
                                      + squared_norm(L_residuals_["l2_bbbb"]));
        }

        // m residual norms
        if (include_u0_) {
            double rm0 = L_scalar_resids_["m0"];
            L_resid_norms_["m0"] = sqrt(rm0 * rm0);
        }
        if (include_u1_) {
            L_resid_norms_["m1"] = sqrt(squared_norm(L_residuals_["m1_aa"])
                                      + squared_norm(L_residuals_["m1_bb"]));
        }
        if (include_u2_) {
            L_resid_norms_["m2"] = sqrt(squared_norm(L_residuals_["m2_aaaa"])
                                      + squared_norm(L_residuals_["m2_abab"])
                                      + squared_norm(L_residuals_["m2_bbbb"]));
        }

        // total residual norm for all amplitudes (if requested)
        if (return_tot) {
            double norm = 0.0;
            for (auto & resid : L_residuals_) {
                if (!resid.second.is_initialized()) continue;
                norm += squared_norm(resid.second);
            }
            return sqrt(norm);
        }
        else return 0.0;
    }
    
    void LambdaDriver::print_properties() {
        // calculate norms of amplitudes
        map<string, double> amp_norms;
        double total_norm = 0.0;
        for (auto& amp : L_amplitudes_) {
            if(!amp.second.is_initialized()) continue;
            // name of amplitude is first two characters of key
            double amp_norm = squared_norm(amp.second);
            string amp_name = amp.first.substr(0,2);
            amp_norms[amp_name] += amp_norm;
            total_norm += amp_norm;
        }

        double nL1 = amp_norms["l1"];
        double nL2 = amp_norms["l2"];
        double nM1 = amp_norms["m1"];
        double nM2 = amp_norms["m2"];

        double m0_ = L_scalar_amps_["m0"];

        // print norms of cluster amplitudes
        outfile->Printf("\n\n   Norms of lambda amplitudes:");
        outfile->Printf("\n   ------------------------------");
        outfile->Printf("\n    T1: %15.12lf | %5.2f %%", sqrt(nL1), 100*nL1/total_norm);
        outfile->Printf("\n    T2: %15.12lf | %5.2f %%", sqrt(nL2), 100*nL2/total_norm);
        if ( include_u0_ ) outfile->Printf("\n    U0: %15.12lf | %5.2f %%", m0_, 100.0*fabs(m0_*m0_)/total_norm);
        if ( include_u1_ ) outfile->Printf("\n    U1: %15.12lf | %5.2f %%", sqrt(nM1), 100.0*nM1/total_norm);
        if ( include_u2_ ) outfile->Printf("\n    U2: %15.12lf | %5.2f %%", sqrt(nM2), 100.0*nM2/total_norm);
        outfile->Printf("\n   ------------------------------");
        outfile->Printf("\n    Total: %15.12lf\n\n", sqrt(total_norm));

    }

    void LambdaDriver::print_iter_header() const {
        Printf("\n");
        Printf("    ==>  Begin %s lambda iterations <==    \n", cc_wfn_->cc_type_.c_str());
        Printf("\n");
        Printf("%5s %16s %15s %15s  | %8s %8s",  "Iter","energy","dE","|dL|","|dL1|","|dL2|");
        if (include_u0_) Printf(" %8s","|dM0|");
        if (include_u1_) Printf(" %8s","|dM1|");
        if (include_u2_) Printf(" %8s","|dM2|");
        Printf("\n");
    }

    void LambdaDriver::print_iteration(size_t iter, double energy, double dele, double tnorm) const {
        Printf("%5i %17.12lf %15.12lf %15.12lf | %-8.1e %-8.1e",iter,energy,dele,tnorm,L_resid_norms_.at("l1"),L_resid_norms_.at("l2"));
        if (include_u0_) Printf(" %-8.1e", L_resid_norms_.at("m0"));
        if (include_u1_) Printf(" %-8.1e", L_resid_norms_.at("m1"));
        if (include_u2_) Printf(" %-8.1e", L_resid_norms_.at("m2"));
        Printf("\n");
    }

} // cc_cavity
