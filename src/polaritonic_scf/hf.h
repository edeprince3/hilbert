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

#ifndef POLARITONIC_HF_H
#define POLARITONIC_HF_H

#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsio/psio.hpp>

#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/molecule.h>
#include <psi4/libqt/qt.h>

using namespace psi;

namespace hilbert{ 

class PolaritonicHF: public Wavefunction {

  public:

    PolaritonicHF(std::shared_ptr<Wavefunction> reference_wavefunction,Options & options);

    ~PolaritonicHF();

    void common_init();

    virtual double compute_energy() {
        throw PsiException("compute_energy has not been implemented for this Polaritonic HF solver",__FILE__,__LINE__);
    }

  protected:

    /// nuclear repulsion energy
    double enuc_;

    /// the multiplicity
    int multiplicity_;

    // cavity related quantities: 

    std::vector< std::shared_ptr<Matrix> > dipole_;

    std::shared_ptr<Matrix> CavityDipolePotential_x_;
    std::shared_ptr<Matrix> CavityDipolePotential_y_;
    std::shared_ptr<Matrix> CavityDipolePotential_z_;
    std::shared_ptr<Matrix> CavityDipole_x_;
    std::shared_ptr<Matrix> CavityDipole_y_;
    std::shared_ptr<Matrix> CavityDipole_z_;
    std::shared_ptr<Matrix> HCavity_x_;
    std::shared_ptr<Matrix> HCavity_y_;
    std::shared_ptr<Matrix> HCavity_z_;
    std::shared_ptr<Matrix> HCavityTotal_x_;
    std::shared_ptr<Matrix> HCavityTotal_y_;
    std::shared_ptr<Matrix> HCavityTotal_z_;
    std::shared_ptr<Matrix> HCavityInteraction_x_;
    std::shared_ptr<Matrix> HCavityInteraction_y_;
    std::shared_ptr<Matrix> HCavityInteraction_z_;

    void build_cavity_hamiltonian();
    void dipole_potential_integrals();

    void InitializeCavity();

    double *cavity_e_, *cavity_tdm_, *cavity_coordinates_, *center_of_mass_;

    long int n_photon_states_;

    double nuc_dip_x_;
    double nuc_dip_y_; 
    double nuc_dip_z_;

    void initialize_cavity();

    bool print_cavity_properties_;

};

} // End namespaces

#endif 
