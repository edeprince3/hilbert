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
    void optimize_orbitals(double * d2, double * d1, double * tei, double * oei, double * transformation_matrix);
    void get_lagrangian(std::shared_ptr<Matrix> Lagrangian);
    void get_hessian(std::shared_ptr<Matrix> Hessian);
    bool is_converged(){ return is_energy_converged_ && is_gradient_converged_; }

  protected:

    bool is_energy_converged_;
    bool is_gradient_converged_;

    int nirrep_;

    int nmo_;
    int nfrzc_;
    int nrstc_;
    int amo_;
    int nrstv_;
    int nfrzv_;

    int * amopi_;
    int * rstcpi_;
    int * rstvpi_;

    Dimension nmopi_;
    Dimension doccpi_;
    Dimension frzcpi_;
    Dimension frzvpi_;

    int * orbital_symmetries_;

    /// options
    double e_convergence_;
    double g_convergence_;
    int maxiter_;
    bool write_;
    bool exact_diagonal_hessian_;
    std::string algorithm_;

};

}

#endif
