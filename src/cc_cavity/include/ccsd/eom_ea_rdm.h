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

#ifndef CC_CAVITY_EOM_EA_RDM_H
#define CC_CAVITY_EOM_EA_RDM_H

#include "../eom_rdm.h"
#include "eom_ea_ccsd.h"

using namespace std;
using namespace psi;
using namespace TA;

namespace hilbert {
    class EOM_EA_RDM : public EOM_RDM {

    public:
        EOM_EA_RDM(const shared_ptr<EOM_Driver>& eom_driver, Options & options);
        ~EOM_EA_RDM() override = default;

        /**
         * Compute the 1-RDMs in MO basis
         */
        void compute_eom_1rdm() override;

    };

} // cc_cavity

#endif //CC_CAVITY_EOM_EA_RDM_H
