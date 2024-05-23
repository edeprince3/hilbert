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

#ifndef CC_CAVITY_EOM_RDM_H
#define CC_CAVITY_EOM_RDM_H

#include "eom_driver.h"

namespace hilbert {
    class EOM_RDM {

    public:
        EOM_RDM(const shared_ptr<EOM_Driver> &eom_driver, Options &options);

        // initialize the eom driver
        shared_ptr<EOM_Driver> eom_driver_;
        Options &options_;

        // initialize the world
        World &world_ = eom_driver_->world_;

        // dimensions of the wavefunction
        size_t o_ = eom_driver_->o_, oa_ = eom_driver_->oa_, ob_ = eom_driver_->ob_, // occupied
        v_ = eom_driver_->v_, va_ = eom_driver_->va_, vb_ = eom_driver_->vb_, // virtual
        ns_ = eom_driver_->ns_, // number of spin orbitals
        nQ_ = eom_driver_->nQ_; // number of auxiliary basis functions

        size_t singleDim_ = eom_driver_->singleDim_,
                doubleDim_ = eom_driver_->doubleDim_;

        size_t dim_e_ = eom_driver_->dim_e_,
                dim_p_ = eom_driver_->dim_p_,
                M_ = eom_driver_->M_,
                N_ = eom_driver_->N_; // dimensions of the problem

        bool include_t3_ = eom_driver_->include_t3_,
                include_t4_ = eom_driver_->include_t4_,
                include_u0_ = eom_driver_->include_u0_,
                include_u1_ = eom_driver_->include_u1_,
                include_u2_ = eom_driver_->include_u2_,
                include_u3_ = eom_driver_->include_u3_,
                include_u4_ = eom_driver_->include_u4_;

        // initialize the rdm tensors
        TArrayMap RDM_blks_;

        // initialize array of properties
        TArrayMap properties_;

        /** ------------------------ Common Functions ------------------------ */

        /**
         * Print the oscillator strengths and x,y,z components of the transition dipole moment for each transition
         */
        virtual void print_oscillators();

        /**
             * Evaluate the oscillator strength and x,y,z components of the transition dipole moment for each transition
             * and store them in the properties_ map with the keys:
             *         "OSCILLATOR STRENGTHS", "X TRANSITION DIPOLES", "Y TRANSITION DIPOLES", "Z TRANSITION DIPOLES"
             */
        virtual void compute_oscillators();

        /** ------------------------ Virtual Functions ------------------------ */

        /**
         * Compute the 1-RDMs in MO basis
         */
        virtual void compute_eom_1rdm() = 0;

        /**
         * compute a 2-RDM in MO basis
         * @param rdm_states: the left and right states of the rdm to save
         */
        virtual void compute_eom_2rdm(vector<int> rdm_states) = 0;

        /**
         * save the density to wavefunction
         * @param rdm_states: the left and right states of the rdm to save
         */
        virtual void save_density(vector<int> rdm_states) = 0;

    };
}

#endif //CC_CAVITY_EOM_RDM_H
