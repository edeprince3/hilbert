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

#include "../include/eom_driver.h"

namespace hilbert {

    EOM_Driver::EOM_Driver(const shared_ptr<CC_Cavity> &cc_wfn, Options &options) :
            cc_wfn_(cc_wfn), options_(options), world_(TA::get_default_world()) {

        // type of eom-cc calculation (e.g. EOM-EE-CCSD)
        eom_type_ = "EOM-" + options_.get_str("EOM_TYPE") + "-" + cc_wfn_->cc_type_;

        // print header for eom-cc calculation
        print_banner();

    }

    void EOM_Driver::print_banner() const {
        size_t pad = 54;
        size_t pad_half = (pad - eom_type_.size()) / 2;
        size_t extra_space = (pad - eom_type_.size()) % 2; // Calculate if there's an extra space needed

        string pad_star_str = string(pad, '*');
        string pad_space_str = string(pad - 2, ' ');
        string pad_half_str = string(pad_half - 1, ' ');

        Printf("\n\n");
        Printf("        %s\n", pad_star_str.c_str());
        Printf("        *%s*\n", pad_space_str.c_str());
        Printf("        *%s*\n", pad_space_str.c_str());
        Printf("        *%s%s%s%*s*\n", pad_half_str.c_str(), eom_type_.c_str(), pad_half_str.c_str(), extra_space,
               ""); // Add the extra space if needed
        Printf("        *%s*\n", pad_space_str.c_str());
        Printf("        *%s*\n", pad_space_str.c_str());
        Printf("        %s\n", pad_star_str.c_str());
    }

    void EOM_Driver::compute_eom_energy() {

        // set dimension of the problem
        set_problem_size();

        if (build_hamiltonian_) M_ = N_;

        // number of desired roots:
        if (M_ == 0) {
            throw PsiException("Number of requested roots is zero. "
                               "Please set the number of roots to a positive integer.", __FILE__, __LINE__);
        }

        // if subspace size is larger than the dimension of the problem, set it to the dimension of the problem
        if (maxdim_ > N_) {
            maxdim_ = N_;
            Printf(
                    "Warning: Subspace size is larger than the dimension of the problem. Setting subspace size to %d\n",
                    N_
            );

            // if number of requested roots is larger than the dimension of the problem, set it to the dimension of the problem
            if (M_ > N_) {
                M_ = N_;
                Printf(
                        "Warning: Number of requested roots is larger than the dimension of the problem. "
                        "Setting number of roots to %d\n", N_
                );
            }
        }

        if (initdim_ > maxdim_) {
            initdim_ = maxdim_;
            Printf(
                    "Warning: Initial subspace size is larger than the maximum subspace size. "
                    "Setting initial subspace size to %d\n", maxdim_
            );
        }

        // print out EOM-CC parameters
        Printf("\n");
        Printf("  ==> %s Parameters <==\n", eom_type_.c_str());
        Printf("\n");
        if (build_hamiltonian_) {
            Printf("  Building Entire EOM Hamiltonian\n");
            Printf("  -- No. All States:        %d\n", N_);
        } else {
            Printf("  Non-Symmetric Davidson Algorithm for EOM Hamiltonian Subspace\n");
            Printf("  -- No. Requested Roots:   %d\n", M_);
            Printf("  -- No. Possible States:   %d\n", N_);
            Printf("  -- Initial Subspace Size: %d\n", initdim_);
            Printf("  -- Maximum Subspace Size: %d\n", maxdim_);
            Printf("  -- Max Iterations:        %d\n", eom_maxiter_);
            Printf("  -- Energy Convergence:    %.3e\n", eom_e_conv_);
            Printf("  -- Residual Convergence:  %.3e\n", eom_r_conv_);
            Printf("\n");
        }

        // compute the energy of the requested roots
        if (build_hamiltonian_) {
            // build the entire EOM Hamiltonian
            build_hamiltonian();
            return;
        }

        // build the guess subspace
        double *Hdiag = build_guess();
        eigvals_ = make_shared<Vector>(M_);
        eigvals_->zero();
        revec_ = make_shared<Matrix>(M_, N_);
        revec_->zero();
        levec_ = make_shared<Matrix>(M_, N_);
        levec_->zero();

        if (read_guess_) load_eigenvectors(load_id_); // read guess eigenvectors from file?


        Printf("\n\n  Building trial independent contractions across iterations...");
        common_timer_.start();
        build_common_ops();
        common_timer_.stop();
        Printf(" Done.\n      Finished in %s\n", common_timer_.elapsed().c_str());

        // solve the eigenvalue problem
        eigsolve_timer_.start();
        eigensolver_->solve(Hdiag,
                            N_,
                            M_,
                            eigvals_->pointer(),
                            revec_->pointer(),
                            levec_->pointer(),
                            [&](int N, int maxdim, int L, double **Q, double **sigmar, double **sigmal) {

                                //TODO: prevent each thread from each storing an instance of the davidson object
                                // this requires a redesign of the davidson object, which is currently not thread safe.
                                // Ideally, only one thread should enter the davidson object at a time and
                                // broadcast the results to the other threads.

                                eigsolve_timer_.stop();
                                build_sigma(L, Q, sigmar, sigmal);
                                eigsolve_timer_.start();

                            }, maxdim_, initdim_, eom_maxiter_, eom_e_conv_, eom_r_conv_, use_res_norm_, read_guess_,
                            eom_shift_);
        free(Hdiag); // free memory
        sigmaOps.clear(); // clear the sigma operators
        sigvec_blks_.clear(); // clear the sigma vectors
        evec_blks_.clear(); // clear the eigenvectors
        eigsolve_timer_.stop();

        // normalize the eigenvectors
        binormalize_states();

        // unpack the final eigenvectors
        unpack_eigenvectors();

        // save the eigenvectors to file if requested
        if (save_evecs_) save_eigenvectors();

        Printf("\n\n    ==>  Dominant Transitions:  <==    \n", eom_type_.c_str());

        // print out dominant transitions for each root
        print_dominant_transitions();

        // print the results
        Printf("\n");
        Printf("    ==>  %s energies:  <==    \n", eom_type_.c_str());
        print_eom_summary();

        // print out the timers
        print_timers();

    }

    void EOM_Driver::build_sigma(int L, double **Q, double **sigmar, double **sigmal) {

        pack_timer_.start();
        // initialize sigma vectors with zeros
        for (int i = 0; i < N_; i++) {
            for (int j = 0; j < maxdim_; j++) {
                memset(sigmar[i], 0, maxdim_ * sizeof(double));
            }
        }

        /// unpack the trial vectors from Q into the appropriate blocks of the left/right eigenvectors
        unpack_trial_vectors((size_t) L, Q);
        world_.gop.fence();
        pack_timer_.stop(false); // don't count this as a separate time

        build_timer_.start();
        /// build the sigma vectors
        build_Hc_cH((size_t) L);
        world_.gop.fence();
        build_timer_.stop();

        /// pack the sigma vectors into sigmar/sigmal
        pack_timer_.start();
        pack_sigma_vectors((size_t) L, sigmar, sigmal);
        world_.gop.fence();

        // clear the trial vectors and sigma vectors
        evec_blks_.clear();
        sigvec_blks_.clear();

        pack_timer_.stop();

        // normalize the eigenvectors
        binormalize_states();

        // print out a summary of the current iteration
        print_eom_summary();

        // print out the timers
        print_timers();
    }

    void EOM_Driver::binormalize_states() {
        double RValue1;
        double RValue2;
        double** rerp_ = revec_->pointer();
        double** relp_ = levec_->pointer();
        for (size_t i = 0; i < M_; i++) {
            RValue1 = C_DDOT(N_, rerp_[i], 1, rerp_[i], 1);
            for (size_t I = 0; I < N_; I++) {
                if (fabs(RValue1) < 1e-16) {
                    continue;
                }
                rerp_[i][I] /= sqrt(fabs(RValue1));
            }

            RValue2 = C_DDOT(N_, rerp_[i], 1, relp_[i], 1);
            if (fabs(RValue2) < 1e-16) {
                continue;
            }
            for (size_t I = 0; I < N_; I++) {
                relp_[i][I] /= fabs(RValue2);
            }
        }
    }

    void EOM_Driver::print_timers() const {

        // print out the time for each step
        auto total_time = (double) (unpack_timer_.get_runtime() + build_timer_.get_runtime() +
                                    pack_timer_.get_runtime());
        string total_time_str = Timer::format_time(total_time);

        Printf("\n");
        Printf("  Time Elapsed:           %s\n", total_time_str.c_str());
        Printf("  Iterations:             %d\n", build_timer_.num_calls());
        Printf("----------------------------------------\n");
        Printf("  Unpack/Pack Time:\n");
        Printf("    - Avg.  %s\n", pack_timer_.average_time().c_str());
        Printf("    - Total %s\n", pack_timer_.elapsed().c_str());
        Printf("  Eigen Solver Time:\n");
        Printf("    - Avg.  %s\n", eigsolve_timer_.average_time().c_str());
        Printf("    - Total %s\n", eigsolve_timer_.elapsed().c_str());
        Printf("  Sigma Build Time:\n");
        Printf("    - Avg.  %s\n", build_timer_.average_time().c_str());
        Printf("    - Total %s\n", build_timer_.elapsed().c_str());
        Printf("----------------------------------------\n");
        Printf("\n");
    }

    void EOM_Driver::save_eigenvectors() const {
        std::shared_ptr<PSIO> psio(new PSIO());
        Printf("Saving eigenvectors to file %s\n", ("evec_" + psio->getpid() + ".bin").c_str());
        world_.gop.fence();
        if (world_.rank() == 0) {
            revec_->save("revec_" + psio->getpid() + ".bin", false, false);
            levec_->save("levec_" + psio->getpid() + ".bin", false, false);
        }
        world_.gop.fence();
    }

    void EOM_Driver::load_eigenvectors(size_t pid) {
        Printf("Loading eigenvectors from file %s\n", ("evec_" + to_string(pid) + ".bin").c_str());
        // extract the PID from the filename
        world_.gop.fence();
        revec_->load("revec_" + to_string(pid) + ".bin");
        levec_->load("levec_" + to_string(pid) + ".bin");
        world_.gop.fence();
    }

} // cc_cavity