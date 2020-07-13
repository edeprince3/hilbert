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

#ifndef PYTHON_HELPERS_H
#define PYTHON_HELPERS_H

#include <v2rdm_doci/v2rdm_solver.h>
#include <doci/doci_solver.h>
#include <pp2rdm/pp2rdm_solver.h>

namespace hilbert{

class DOCIHelper{

  public:

    DOCIHelper(SharedWavefunction reference_wavefunction,Options & options);
    ~DOCIHelper();
    void common_init();
    double compute_energy();

  protected:

    /// the DOCISolver 
    std::shared_ptr<DOCISolver> doci;

};

class pp2RDMHelper{

  public:

    pp2RDMHelper(SharedWavefunction reference_wavefunction,Options & options);
    ~pp2RDMHelper();
    void common_init();
    double compute_energy();

  protected:

    /// the pp2RDMSolver 
    std::shared_ptr<pp2RDMSolver> pp2rdm;

};

class v2RDMHelper{

  public:

    v2RDMHelper(SharedWavefunction reference_wavefunction,Options & options);
    ~v2RDMHelper();
    void common_init();
    double compute_energy();

  protected:

    /// the v2RDMSolver 
    std::shared_ptr<v2RDMSolver> v2rdm;

};

}

#endif
