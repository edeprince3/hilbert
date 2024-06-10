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

    /// nuclear repulsion energy
    double enuc_;

    /// part of the dipole self energy: 1/2 lambda^2 (dn - <d>)^2 = 1/2 lambda^2 <de>^2
    double average_electric_dipole_self_energy_;

    /// parameters for the cavity
    double * cavity_frequency_;
    double * cavity_coupling_strength_;

    /// nuclear dipole moments
    double nuc_dip_x_;
    double nuc_dip_y_;
    double nuc_dip_z_;

    /// electronic dipole moments
    double e_dip_x_;
    double e_dip_y_;
    double e_dip_z_;

    /// total dipole moments
    double tot_dip_x_;
    double tot_dip_y_;
    double tot_dip_z_;

    /// the multiplicity
    int multiplicity_;

    std::vector< std::shared_ptr<Matrix> > dipole_;
    std::vector< std::shared_ptr<Matrix> > quadrupole_;
    std::shared_ptr<Matrix> dipole_scaled_sum_;
    std::shared_ptr<Matrix> quadrupole_scaled_sum_;
    std::shared_ptr<Matrix> scaled_e_n_dipole_squared_;

    long int n_photon_states_;

    /// do use coherent-state basis? default true
    bool use_coherent_state_basis_ = true;

  protected:
    /// evaluate orbital gradient
    std::shared_ptr<Matrix> OrbitalGradient(std::shared_ptr<Matrix> D,
                                            std::shared_ptr<Matrix> F,
                                            std::shared_ptr<Matrix> Shalf);

    /// set up the cavity
    void initialize_cavity();
    void update_cavity_terms();

    /// evaluate constant and one- and two-electron components of the dipole self energy
    void evaluate_dipole_self_energy();

    /// evaluate constant and one- and two-electron components of the dipole variance
    void evaluate_dipole_variance();

};

} // End namespaces

#endif 
