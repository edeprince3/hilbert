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
        void init_integrals() override;

        void init_amplitudes() override;

        void print_iter_header() const override;

        void print_iteration(size_t iter, double energy, double dele, double tnorm) const override;

        void transform_integrals(bool use_t1) override;

        void build_oei(TArrayMap &CL, TArrayMap &CR) override;

        void build_tei(TArrayMap &CL, TArrayMap &CR) override;

        void build_eps() override;

        double build_residuals() override;

        double update_amplitudes() override;

        void print_properties() override;

        void unpack_eris(TArrayMap &Qmo_blks) override;
    };

} // cc_cavity

#endif //CC_CAVITY_QED_CC_H
