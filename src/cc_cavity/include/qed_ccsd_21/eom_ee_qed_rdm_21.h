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

#ifndef CC_CAVITY_EOM_EE_QED_RDM_H
#define CC_CAVITY_EOM_EE_QED_RDM_H

#include "../ccsd/eom_ee_ccsd.h"
#include "../ccsd/eom_ee_rdm.h"

using namespace std;
using namespace psi;
using namespace TA;

namespace hilbert {
    class EOM_EE_QED_RDM_21 : public EOM_EE_RDM {

    public:
        EOM_EE_QED_RDM_21(const shared_ptr<EOM_Driver>& eom_driver, Options & options);
        ~EOM_EE_QED_RDM_21() override = default;

        /**
         * Compute the 1-RDMs in MO basis
         */
        void compute_eom_1rdm() override;

        /**
         * Compute the 2-RDMs in MO basis
         */
        void compute_eom_2rdm() override;
        void rdm2_21_1(), rdm2_21_2(), rdm2_21_3(), rdm2_21_4(), rdm2_21_5(), rdm2_21_6(), rdm2_21_7(), rdm2_21_8();

    };

} // cc_cavity

#endif //CC_CAVITY_EOM_EE_QED_RDM_H
