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

#ifndef CC_CAVITY_QED_CC_H
#define CC_CAVITY_QED_CC_H

#include "../cc_cavity.h"
#include <tiledarray.h>
#include "../../../polaritonic_scf/hf.h"
#include "../../misc/nonsym_davidson_solver_qed.h"
#include "../../misc/ta_helper.h"
#include <psi4/libpsio/psio.hpp>
#include "../../misc/diis_qed.h"
#include "../../misc/timer.h"

using namespace std;
using namespace TA;
using namespace psi;
using namespace Helper;

namespace hilbert {

    class QED_CC : public CC_Cavity {

    public:
        QED_CC(const shared_ptr<Wavefunction> &reference_wavefunction, Options & options);
        ~QED_CC() override = default;

    protected:

        /**
         * @brief Initialize the CC operators
         */
        void init_operators() override;

        /**
         * @brief print the dimensions of the CC operators
         */
        void print_dimensions() override;

        /**
         * @brief Print the header for the iteration
         */
        void print_iter_header() const override;

        /**
         * @brief Print the iteration information
         * @param iter number of current iteration
         * @param energy current energy
         * @param dele energy change
         * @param tnorm current residual norm
         */
        void print_iteration(size_t iter, double energy, double dele, double tnorm) const override;

        /**
         * @brief Build the residual vectors and compute the energy
         * @return the energy
         */
        double build_residuals() override;

        /**
         * @brief Update the amplitudes
         */
        void update_amplitudes() override;

        /**
         * @brief compute the residual norms for each operator
         * @return the total residual norm
         */
        double compute_residual_norms() override;

        /**
         * @brief Print summary of amplitude properties
         */
        void print_properties() override;

    };

} // cc_cavity

#endif //CC_CAVITY_QED_CC_H
