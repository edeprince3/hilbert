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

#include <psi4/libpsi4util/process.h>
#include "cc_cavity/include/eom_driver.h"

namespace hilbert {

    EOM_Driver::EOM_Driver(shared_ptr<CC_Cavity> &cc_wfn, Options &options) :
            cc_wfn_(cc_wfn), options_(options), world_(TA::get_default_world()) {

        // type of eom-cc calculation (e.g. EOM-EE-CCSD)
        eom_type_ = "EOM-" + options_.get_str("EOM_TYPE") + "-" + cc_wfn_->cc_type_;
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

        // print header for eom-cc calculation
        print_banner();

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

        // allocate memory for eigenvectors and eigenvalues
        eigvals_ = make_shared<Vector>(M_);
        revec_ = make_shared<Matrix>(M_, N_);
        levec_ = make_shared<Matrix>(M_, N_);
        eigvals_->zero();
        revec_->zero();
        levec_->zero();

        // compute the energy of the requested roots
        eigsolve_timer_.start();
        if (build_hamiltonian_) {
            // build the entire EOM Hamiltonian
            build_hamiltonian();
        } else {
            // build the EOM Hamiltonian subspace
            build_eom_subspace();
        }
        eigsolve_timer_.stop();

        // normalize the eigenvectors
        binormalize_states();

        // unpack the final eigenvectors into tiled array objects
        unpack_eigenvectors();

        // save the eigenvectors to file if requested
        if (save_evecs_) save_eigenvectors();

        // print out dominant transitions for each root
        transitions_summary();

        // print the results
        Printf("\n");
        Printf("    ==>  %s energies:  <==    \n", eom_type_.c_str());
        print_eom_summary();

        // print out the transition dipole moments

        // print out the timers
        print_timers();

    }

    void EOM_Driver::build_eom_subspace() {
        // build the guess subspace
        double *Hdiag = build_preconditioner();
        if (read_guess_)
            load_eigenvectors(load_id_); // read guess eigenvectors from file?


        Printf("\n\n  Building trial independent contractions across iterations...");
        common_timer_.start();
        build_common_ops();
        common_timer_.stop();
        Printf(" Done.\n      Finished in %s\n", common_timer_.elapsed().c_str());

        // set convergence criteria
        eigensolver_->eom_maxiter = eom_maxiter_;
        eigensolver_->eom_e_conv = eom_e_conv_;
        eigensolver_->eom_r_conv = eom_r_conv_;

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

                            }, maxdim_, initdim_, eom_r_conv_, use_res_norm_);

        // free memory

        free(Hdiag); // free the preconditioner
        sigmaOps.clear(); // clear the sigma operators
        reuse_tmps_.clear(); // clear the temporary operators
        sigvec_blks_.clear(); // clear the sigma vectors
        evec_blks_.clear(); // clear the eigenvectors
    }

    void EOM_Driver::build_sigma(int L, double **Q, double **sigmar, double **sigmal) {

        pack_timer_.start();

        // initialize sigma vectors with zeros
        memset(*sigmar, 0, maxdim_ * N_ * sizeof(double));
        memset(*sigmal, 0, maxdim_ * N_ * sizeof(double));

        auto Lsize = static_cast<size_t>(L);

        /// unpack the trial vectors from Q into the appropriate blocks of the left/right eigenvectors
        unpack_trial_vectors(Lsize, Q);
        world_.gop.fence();
        pack_timer_.stop(false); // don't count this as a separate instance of the timer

        build_timer_.start();
        /// build the sigma vectors
        build_Hc_cH(Lsize);
        world_.gop.fence();
        build_timer_.stop();

        /// pack the sigma vectors into sigmar/sigmal
        pack_timer_.start();
        pack_sigma_vectors(Lsize, sigmar, sigmal);
        world_.gop.fence();

        // clear the trial vectors and sigma vectors
        evec_blks_.clear();
        sigvec_blks_.clear();

        // normalize the eigenvectors
        binormalize_states();

        // print out a summary of the current iteration
        print_eom_summary();

        // print out the timers
        pack_timer_.stop();
        print_timers();

        // set environment variables using target state
        int target_state;

        vector<int> rdm_states = options_.get_int_vector("RDM_STATES");
        if (rdm_states.empty()) target_state = 0; // default to ground state
        else target_state = rdm_states.front(); // use the first state in the list

        Process::environment.globals["EOM TARGET ENERGY"] = eigvals_->get(target_state);
        Process::environment.globals["CURRENT ENERGY"]    = eigvals_->get(target_state);
    }

    void EOM_Driver::print_eom_summary() {

        // print out the header for the table
        print_eom_header();

        double last_en = 0.0;
        bool no_degeneracy = options_.get_bool("NO_DEGENERACY");

        ground_energy_ref_ = eigvals_->get(0);

        // loop over roots
        for (int i = 0; i < M_; i++) {
            // calculate norms of amplitudes
            double *norms = get_state_norms(i);

            // print energies and excitation energies
            double ee_energy = eigvals_->get(i);
            if (fabs(last_en - ee_energy) < 1e-10 && no_degeneracy)
                continue;
            last_en = ee_energy;

            Printf("%5d %20.12lf %17.12lf ", i, ee_energy,
                   ee_energy - ground_energy_ref_);

            // print out the norm of the amplitudes. Ignore if less than 1e-10
            for (int j = 0; j < nops_; j++) {
                if (fabs(norms[j]) > 1e-10) {
                    Printf("%13.10lf ", norms[j]);
                } else {
                    Printf("------------- ");
                }
            }
            Printf("\n");

            // free memory
            free(norms);
        }
    }

    void EOM_Driver::binormalize_states() {
        double RValue1, RValue2;
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
            if (fabs(RValue2) < 1e-16) continue;

            for (size_t I = 0; I < N_; I++) {
                relp_[i][I] /= fabs(RValue2);
            }
        }
    }

    void EOM_Driver::transitions_summary(){

        // print out the dominant transitions for each root if requested
        if (options_.get_bool("PRINT_TRANSITIONS"))
            Printf("\n\n    ==>  Dominant Transitions:  <==    \n", eom_type_.c_str());
        else return;

        for (size_t I = 0; I < M_; I++) {
            // get the dominant transitions for state I
            DominantTransitionsType dominant_transitions = find_dominant_transitions(I);
            size_t num_print = options_.get_int("NUM_PRINT_TRANSITIONS");
            // print the dominant transitions
            Printf("\n\n  --> Root %zu <--\n", I);

            // print top `n` transitions for each excitation block
            for(auto transitionBlock : dominant_transitions) {
                // grab first 5 values from priority queue
                std::deque<TransitionType> topTransitions;

                for (size_t i = 0; i < num_print; i++) {
                    if (!transitionBlock.second.empty()) {
                        topTransitions.push_back(transitionBlock.second.top());
                        transitionBlock.second.pop();
                    }
                }

                // print top 5 in following format:
                //     --> l1*r1_bbbb <--
                //          1: (  i  )->(  a  ), lr=(+-)0.1234, l=(+-)0.1234, r=(+-)0.1234
                //          2: (  i  )->(  a  ), lr=(+-)0.1234, l=(+-)0.1234, r=(+-)0.1234
                //          3: (  i  )->(  a  ), lr=(+-)0.1234, l=(+-)0.1234, r=(+-)0.1234
                //          4: (  i  )->(  a  ), lr=(+-)0.1234, l=(+-)0.1234, r=(+-)0.1234
                //          5: (  i  )->(  a  ), lr=(+-)0.1234, l=(+-)0.1234, r=(+-)0.1234
                //     --> l2*r2_aaaa <--
                //          1: (  i,  j  )->(  a,  b  ), lr=(+-)0.1234, l=(+-)0.1234, r=(+-)0.1234
                //          2: (  i,  j  )->(  a,  b  ), lr=(+-)0.1234, l=(+-)0.1234, r=(+-)0.1234
                //          3: (  i,  j  )->(  a,  b  ), lr=(+-)0.1234, l=(+-)0.1234, r=(+-)0.1234
                //          4: (  i,  j  )->(  a,  b  ), lr=(+-)0.1234, l=(+-)0.1234, r=(+-)0.1234
                //          5: (  i,  j  )->(  a,  b  ), lr=(+-)0.1234, l=(+-)0.1234, r=(+-)0.1234

                size_t top_size = topTransitions.size();
                if (top_size == 0) continue;

                Printf("\n    %s", transitionBlock.first.c_str());
                for (size_t i = 0; i < top_size; i++) {
                    TransitionType &tpair = topTransitions[i];

                    auto [lr, l, r, spin, labels] = tpair;
                    size_t size = labels.size(); // size of transitionBlock pairs
                    if (size == 0) { // ground state
                        Printf(": l*r = %15.12lf | l = %11.8lf, r = %11.8lf", i+1, lr, l, r);
                    } else if (size == 1) { // singles
                        Printf("\n        %3d: (+%2d%c), l*r = %15.12lf | l = %11.8lf, r = %11.8lf", i+1,
                               labels[0], spin[0], lr, l, r);
                    } else if (size == 2) { // singles
                        Printf("\n        %3d: (%2d%c ) -> (%2d%c ), l*r = %15.12lf | l = %11.8lf, r = %11.8lf", i+1,
                               labels[1], spin[1], labels[0], spin[0], lr, l, r);
                    } else if (size == 3) { // singles
                        Printf("\n        %3d: (%2d%c ) -> (%2d%c,+%2d%c ), l*r = %15.12lf | l = %11.8lf, r = %11.8lf", i+1,
                               labels[2], spin[2], labels[1], spin[1], labels[0], spin[0], lr, l, r);
                    } else if (size == 4) { // doubles
                        Printf("\n        %3d: (%2d%c,%2d%c ) -> (%2d%c,%2d%c ), l*r = %15.12lf | l = %11.8lf, r = %11.8lf", i+1,
                               labels[2], spin[2], labels[3], spin[3],
                               labels[0], spin[0], labels[1], spin[1], lr, l, r);
                    } else {
                        throw PsiException("Too many transitions for current EOM method. Update EOM_Driver::transitions_summary", __FILE__, __LINE__);
                    }
                }
            }

            // free memory
            dominant_transitions.clear();
        }
        Printf("\n\n");
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
        world_.gop.serial_invoke([=] () {
            revec_->save("revec_" + psio->getpid() + ".bin", false, false);
            levec_->save("levec_" + psio->getpid() + ".bin", false, false);
        });
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
