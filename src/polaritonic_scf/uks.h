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

#ifndef POLARITONIC_UKS_H
#define POLARITONIC_UKS_H

#include "hf.h"

using namespace psi;

namespace hilbert{ 

class PolaritonicUKS: public PolaritonicHF {

  public:

    PolaritonicUKS(std::shared_ptr<Wavefunction> reference_wavefunction,Options & options);

    ~PolaritonicUKS();

    void common_init();

    double compute_energy();

  protected:

    std::shared_ptr<Matrix> Va_;
    std::shared_ptr<Matrix> Vb_;

    /// cavity Hamiltonian in basis of photon number states
    std::shared_ptr<Matrix> HCavity_x_;
    std::shared_ptr<Matrix> HCavity_y_;
    std::shared_ptr<Matrix> HCavity_z_;

    /// molecule->cavity interactoin Hamiltonian in basis of photon number states
    std::shared_ptr<Matrix> HCavityInteraction_x_;
    std::shared_ptr<Matrix> HCavityInteraction_y_;
    std::shared_ptr<Matrix> HCavityInteraction_z_;

    /// total cavity Hamiltonian in basis of photon number states
    std::shared_ptr<Matrix> HCavityTotal_x_;
    std::shared_ptr<Matrix> HCavityTotal_y_;
    std::shared_ptr<Matrix> HCavityTotal_z_;

    /// cavity dipole operator in basis of photon number states
    std::shared_ptr<Matrix> CavityDipole_x_;
    std::shared_ptr<Matrix> CavityDipole_y_;
    std::shared_ptr<Matrix> CavityDipole_z_;

    /// build and diagonalize cavity hamiltonian
    double build_cavity_hamiltonian();

    /// do print cavity properties?
    bool print_cavity_properties_ = false;

};

} // End namespaces

#endif 
