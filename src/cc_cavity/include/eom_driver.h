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

#ifndef CC_CAVITY_EOM_DRIVER_H
#define CC_CAVITY_EOM_DRIVER_H

#include "cc_cavity.h"

using namespace std;
using namespace TA;
using namespace psi;
using namespace Helper;

namespace hilbert {

    class EOM_Driver {
    public:

        EOM_Driver(const shared_ptr<CC_Cavity> &cc_wfn, Options & options);
        ~EOM_Driver() = default;

        // world object
        World & world_;

        // thread safe print function
        inline void Printf(const char *format, ...) const {
            if (world_.rank() == 0) {
                va_list argptr;
                va_start(argptr, format);
                char input[2048];
                vsprintf(input, format, argptr);
                va_end(argptr);
                outfile->Printf("%s", input);
            }
            world_.gop.fence();
        }

        /**
         * Calculate the index of an element in an upper triangular array from its row and column indices (i,j)
         * @param i row index
         * @param j column index
         * @param N leading dimension of the array
         * @return index of the element in the triangular array
         */
        inline static size_t sqr_2_tri_idx(size_t i, size_t j, size_t N) {
            return (2 * N - i - 3) * i / 2 + j - 1;
        }

        /** ------------------------ Common Attributes ------------------------ */

        // CC_Cavity object
        std::shared_ptr<CC_Cavity> cc_wfn_;
        Options & options_;


        // dimensions of the wavefunction
        size_t o_ = cc_wfn_->o_, oa_ = cc_wfn_->oa_, ob_ = cc_wfn_->ob_, // occupied
        v_ = cc_wfn_->v_, va_ = cc_wfn_->va_, vb_ = cc_wfn_->vb_, // virtual
        ns_ = cc_wfn_->ns_, // number of spin orbitals
        nQ_ = cc_wfn_->nQ_, // number of auxiliary basis functions
        nops_; // number of operators

        // dimensions of operators
        size_t singleDim_ = cc_wfn_->singleDim_,
                doubleDim_ = cc_wfn_->doubleDim_,
                tripleDim_ = cc_wfn_->tripleDim_,
                quadDim_ = cc_wfn_->quadDim_;

        // dimensions of the problem
        size_t dim_e_, dim_p_, N_; // number of electronic states, number of polaritonic states, total dimension
        size_t M_ = options_.get_int("NUMBER_ROOTS"); // number of desired roots
        size_t maxdim_ = M_ * options_.get_int("MAXDIM"); // maximum subspace size scaled by number of roots
        size_t initdim_ = M_ * options_.get_int("INDIM"); // initial subspace size scaled by number of roots

        // initialize included operator bools
        bool include_t3_ = cc_wfn_->include_t3_;
        bool include_t4_ = cc_wfn_->include_t4_;
        bool include_u0_ = cc_wfn_->include_u0_;
        bool include_u1_ = cc_wfn_->include_u1_;
        bool include_u2_ = cc_wfn_->include_u2_;
        bool include_u3_ = cc_wfn_->include_u3_;
        bool include_u4_ = cc_wfn_->include_u4_;

        /// eom-cc parameters

        bool build_hamiltonian_ = options_.get_bool("BUILD_HAMILTONIAN"); // build hamiltonian?
        double eom_e_conv_ = options_.get_double("EOM_E_CONV"); // energy convergence
        double eom_r_conv_ = options_.get_double("EOM_R_CONV"); // residual convergence
        size_t eom_maxiter_ = options_.get_int("EOM_MAXITER"); // maximum number of iterations
        bool read_guess_ = options_.get_int("LOAD_ID") != -1; // read guess from file?
        bool eom_ss_guess_ = options_.get_bool("EOM_SS_GUESS"); // use singles hamiltonian for guess?
        double eom_shift_ = options_.get_double("EOM_SHIFT"); // use shift?
//        bool use_res_norm_ = options_.get_bool("USE_RES_NORM"); // use residual norm?
        bool use_res_norm_ = true; // residual norm is always used to adjust the subspace size
        bool excited_transitions_ = options_.get_bool("EXCITED_PROPERTIES"); // compute excited state properties?

        bool save_evecs_ = options_.get_bool("SAVE_EVECS"); // save eigenvectors to file?
        size_t load_id_ = options_.get_int("LOAD_ID"); // id of file to load

        // eigenvalue solver (Non-symmetric Davidson, thanks to the wise Nam Vu)
        shared_ptr<Nonsym_DavidsonSolver_QED> eigensolver_ = make_shared<Nonsym_DavidsonSolver_QED>(Nonsym_DavidsonSolver_QED());

        /// eom-cc variables

        string eom_type_; // eom type
        vector<TArrayD> sigmaOps; // trial independent operators
        SharedVector eigvals_; // left/right eigenvalue blocks
        SharedMatrix revec_, levec_; // left/right eigenvalue blocks
        TArrayMap evec_blks_;   // left/right trial vector blocks
        TArrayMap sigvec_blks_;   // left/right trial eigenvalue*eigenvector blocks
        Timer common_timer_, build_timer_, unpack_timer_, pack_timer_, eigsolve_timer_; // timers

        /** ------------------------ Common Functions ------------------------ */

        /**
         * Print the banner for the EOM-CC calculation
         */
        virtual void print_banner() const;

        /**
         * Compute the EOM-CC energies
         */
        virtual void compute_eom_energy();

        /**
         * Prepare and process data for the EOM-CC sigma equations
         * @param L number of trial vectors
         * @param Q trial vectors
         * @param sigmar right sigma vectors
         * @param sigmal left sigma vectors
         */
        virtual void build_sigma(int L, double **Q, double **sigmar, double **sigmal);

        /**
         * Print the EOM-CC timers
         */
        virtual void print_timers() const;

        /**
         * Save the eigenvectors to file (not working)
         */
        virtual void save_eigenvectors() const;

        /**
         * Load the eigenvectors from file (not working)
         * @param pid process id
         */
        virtual void load_eigenvectors(size_t pid);

        /**
         * Binormalize the left and right eigenvectors
         */
        virtual void binormalize_states();

        /** ------------------------ Interface Functions ------------------------ */

        /**
         * Set the problem size for the EOM-CC calculation
         */
        virtual void set_problem_size() = 0;

        /**
         * Print the header for the EOM-CC summary of energies and properties
         */
        virtual void print_eom_header() const = 0;

        /**
         * Print the EOM-CC summary of energies and properties
         */
        virtual void print_eom_summary() const;

        /**
         * build the full hamiltonian (not implemented)
         */
        virtual void build_hamiltonian() = 0;

        /**
         * Build the preconditioner for the EOM-CC calculation (either uses orbital guess or singles guess)
         * @return guess vector
         */
        virtual double* build_preconditioner() = 0;

        /**
         * Build operators for the EOM-CC calculation that do not depend on the trial vectors
         */
        virtual void build_common_ops() = 0;

        /**
         * unpack the trial vectors from the array Q
         * @param L number of trial vectors
         * @param Q trial vectors
         */
        virtual void unpack_trial_vectors(size_t L, double **Q) = 0;

        /**
         * Solve the EOM-CC sigma equations for the trial vectors
         * @param L
         */
        virtual void build_Hc_cH(size_t L) = 0;

        /**
         * Pack the sigma vectors into the array sigmar
         * @param L number of trial vectors
         * @param sigmar right sigma vectors
         * @param sigmal left sigma vectors
         */
        virtual void pack_sigma_vectors(size_t L, double **sigmar, double **sigmal) = 0;

        /**
         * Compute the norms of each operator in the i'th eigenvector
         * @param i the i'th eigenvector
         * @return the norms of each operator in the i'th eigenvector
         */
        virtual double *get_state_norms(size_t i) const = 0;

        /**
         * unpack the final eigenvectors into TArrays
         */
        virtual void unpack_eigenvectors() = 0;


        // a map of the operator name to a priority queue of its dominant transitions with:
        //     the magnitude of the transition
        //     the spin of the transition
        //     the indicies of the transition
        typedef map<string, priority_queue<pair<double, pair<string, vector<size_t>>>>> DominantTransitionsType;

        /**
         * find the dominant transitions for a given state
         * @param I the state
         */
        virtual DominantTransitionsType find_dominant_transitions(size_t I) = 0;

        /**
         * Print the dominant transitions for all states
         */
        virtual void transitions_summary();
    };

} // cc_cavity

#endif //CC_CAVITY_EOM_DRIVER_H
