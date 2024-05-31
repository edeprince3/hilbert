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

#ifndef CC_CAVITY_LAMDA
#define CC_CAVITY_LAMDA

#include "tiledarray.h"
#include "cc_cavity/include/cc_cavity.h"
#include "polaritonic_scf/hf.h"
#include "misc/nonsym_davidson_solver.h"
#include "cc_cavity/misc/ta_helper.h"
#include <psi4/libpsio/psio.hpp>
#include "cc_cavity/misc/ta_diis.h"
#include "cc_cavity/misc/timer.h"

using namespace std;
using namespace TA;
using namespace psi;
using namespace Helper;

namespace hilbert {

    class LambdaDriver {

    public:
        LambdaDriver(shared_ptr<CC_Cavity> &cc_ref, Options & options) : cc_wfn_(cc_ref) {}
        ~LambdaDriver() = default;

        /**
         * @brief compute the lambda energy
         * @return the energy
         */
        double compute_lambda();
        double lam_energy_ = 0.0; // lambda energy

    protected:

        // thread safe print function
        inline void Printf(const char *format, ...) const {
            va_list argptr;
            va_start(argptr, format);
            char input[1024];
            vsprintf(input, format, argptr);
            va_end(argptr);
            world_.gop.serial_invoke(
                    [=]() {
                        outfile->Printf("%s", input);
                    }
            );
        }

        shared_ptr<CC_Cavity> cc_wfn_;

        World &world_ = TA::get_default_world();
        Options &options_ = cc_wfn_->options();

        TArrayMap &T_amplitudes_ = cc_wfn_->amplitudes_;
        TArrayMap L_amplitudes_;
        map<string, double> L_scalar_amps_; // scalar residuals

        TArrayMap L_residuals_;
        map<string, double> L_scalar_resids_; // scalar residuals
        map<string, double> L_resid_norms_;

        // dimensions of the wavefunction
        size_t o_ = cc_wfn_->o_, oa_ = cc_wfn_->oa_, ob_ = cc_wfn_->ob_, // occupied
        v_ = cc_wfn_->v_, va_ = cc_wfn_->va_, vb_ = cc_wfn_->vb_, // virtual
        ns_ = cc_wfn_->ns_, // number of spin orbitals
        nQ_ = cc_wfn_->nQ_; // number of auxiliary basis functions

        // dimensions of operators
        size_t singleDim_ = cc_wfn_->singleDim_,
                doubleDim_ = cc_wfn_->doubleDim_,
                tripleDim_ = cc_wfn_->tripleDim_,
                quadDim_ = cc_wfn_->quadDim_;
        vector<string> &idx_map_ = cc_wfn_->idx_map_;

        // cc_type
        string cc_type_ = cc_wfn_->cc_type_;

        // initialize included operator bools
        bool include_t3_ = cc_wfn_->include_t3_;
        bool include_t4_ = cc_wfn_->include_t4_;
        bool include_u0_ = cc_wfn_->include_u0_;
        bool include_u1_ = cc_wfn_->include_u1_;
        bool include_u2_ = cc_wfn_->include_u2_;
        bool include_u3_ = cc_wfn_->include_u3_;
        bool include_u4_ = cc_wfn_->include_u4_;

        /// timers
        Timer t_resid, t_ampUp, t_transform, t_oei, t_tei, t_ground;

        /**
         * @brief Initialize the CC operators
         */
        void init_operators();

        /**
         * @brief Print the header for the iteration
         */
        void print_iter_header() const;

        /**
         * @brief Print the iteration information
         * @param iter number of current iteration
         * @param energy current energy
         * @param dele energy change
         * @param tnorm current residual norm
         */
        void print_iteration(size_t iter, double energy, double dele, double tnorm) const;

        /**
         * @brief Build the residual vectors and compute the energy
         * @return the energy
         */
        double build_residuals();

        /**
         * @brief Build the intermediates for the lambda iterations
         */
        void build_intermediates();
        TArrayMap sharedOps;

        /**
         * @brief reset all tiles of a TA object
         * @param array the TA object
         */
        static void zero_tiles(TArrayD& array) {
            static TA::World& world = TA::get_default_world();
            array = TArrayD(world, array.trange()); array.fill(0.0);
        }

        /**
         * @brief Update the amplitudes
         */
        void update_amplitudes();

        /**
         * @brief Extrapolate the amplitudes using DIIS
         */
        void extrapolate_amplitudes();


        /**
         * @brief compute the residual norms for each operator
         * @return the total residual norm
         */
        double compute_residual_norms(bool return_tot = true);

        /**
         * @brief perform the lambda iterations
         * @return the current energy
         */
        double lambda_iterations();

        /**
         * @brief Print summary of amplitude properties
         */
        void print_properties();

    };

} // cc_cavity

#endif //CC_CAVITY_LAMDA
