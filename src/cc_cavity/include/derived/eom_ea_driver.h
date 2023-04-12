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

#ifndef CC_CAVITY_EOM_EA_DRIVER_H
#define CC_CAVITY_EOM_EA_DRIVER_H

#include "../eom_driver.h"

using namespace std;
using namespace TA;
using namespace psi;
using namespace Helper;

namespace hilbert {

class EOM_EA_Driver : public EOM_Driver {
public:
    EOM_EA_Driver(const shared_ptr<CC_Cavity> &cc_wfn, Options & options);
    ~EOM_EA_Driver() = default;

    double old_e_dip_z_; // store the old dipole moment for backtracking
    double old_dipole_self_energy_; // store the old dipole_self_energy_ for backtracking

    void set_problem_size() override;

    void print_eom_summary() const override;

    void build_hamiltonian() override;
    double* build_guess() override;

    void build_common_ops() override;

    void unpack_trial_vectors(size_t L, double **Q) override;

    void pack_sigma_vectors(size_t L, double **sigmar, double **sigmal) override;

    void build_Hc_cH(size_t L) override;

    double *get_state_norms(size_t i) const override;

    void unpack_eigenvectors() override;

    void print_dominant_transitions() override;
};

} // cc_cavity

#endif //CC_CAVITY_EOM_EA_DRIVER_H
