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

#ifndef ORBITAL_OPTIMIZER_H
#define ORBITAL_OPTIMIZER_H

#include <string>

#include <psi4/libplugin/plugin.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/matrix.h>

using namespace psi;

namespace hilbert{

class OrbitalOptimizer {

  public:

    OrbitalOptimizer(std::shared_ptr<Wavefunction> reference_wavefunction, Options & options);
    ~OrbitalOptimizer();

    /// take an orbital optimization step
    void optimize_orbitals(double * d2, double * d1, double * tei, double * oei, double * transformation_matrix);

    /// returns the orbital lagrangian
    void get_lagrangian(std::shared_ptr<Matrix> Lagrangian);

    /// returns the orbital hessian
    void get_hessian(std::shared_ptr<Matrix> Hessian);

    /// has the optimization converged?
    bool is_converged(){ return is_energy_converged_ && is_gradient_converged_; }

  protected:

    /// is the energy converged to e_convergence_?
    bool is_energy_converged_;

    /// is the gradient converged to g_convergence_?
    bool is_gradient_converged_;

    /// number of irreps
    int nirrep_;

    /// number of molecular orbitals
    int nmo_;

    /// number of frozen core orbitals
    int nfrzc_;

    /// number of restricted core orbitals
    int nrstc_;

    /// number of active orbitals
    int amo_;

    /// number of restricted virtual orbitals
    int nrstv_;

    /// number of frozen virtual orbitals
    int nfrzv_;

    /// active orbitals per irrep
    int * amopi_;

    /// restricted core orbitals per irrep
    int * rstcpi_;

    /// restricted virtual orbitals per irrep
    int * rstvpi_;

    /// molecular orbitals per irrep
    Dimension nmopi_;

    /// doubly occuiped  orbitals per irrep
    Dimension doccpi_;

    /// frozen core orbitals per irrep
    Dimension frzcpi_;

    /// frozen virtual orbitals per irrep
    Dimension frzvpi_;

    /// orbital symmetries, in energy order
    int * orbital_symmetries_;

    /// energy convergence
    double e_convergence_;

    /// orbital gradient convergence
    double g_convergence_;

    /// maximum number of iterations
    int maxiter_;

    /// do print iteration information?
    bool write_;

    /// do evaluate the exact diagonal hessian?
    bool exact_diagonal_hessian_;

    /// do rotate active-active orbital pairs?
    bool active_active_rotations_;

    /// algorithm
    std::string algorithm_;

};

}

#endif
