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

#include <v2rdm_doci/v2rdm_doci_solver.h>
#include <v2rdm_casscf/v2rdm_solver.h>
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

class v2RDM_DOCIHelper{

  public:

    v2RDM_DOCIHelper(SharedWavefunction reference_wavefunction,Options & options);
    ~v2RDM_DOCIHelper();
    void common_init();
    double compute_energy();

  protected:

    /// the v2RDMSolver 
    std::shared_ptr<v2RDM_DOCISolver> v2rdm_doci;

};

class v2RDMHelper{

  public:

    v2RDMHelper(SharedWavefunction reference_wavefunction,Options & options);
    ~v2RDMHelper();
    void common_init();

    double compute_energy();

    /// return spin-free one-particle density matrix. full space. sparse
    std::vector<opdm> get_opdm_sparse(std::string type);

    /// return spin-free two-particle density matrix. full space. sparse.
    std::vector<tpdm> get_tpdm_sparse(std::string type);

    /// return spin-free one-particle density matrix as shared matrix for python API
    std::shared_ptr<Matrix> get_opdm();

    /// return spin-free two-particle density matrix as shared matrix for python API
    std::shared_ptr<Matrix> get_tpdm();

    /// return subset of orbitals for python API
    /**!
     * Similar to wavefunction.Ca_subset(); however, this version knows about all of the CI
     * subspaces in the SO basis. We stick to the definitions used by DETCI for now.
     * @param  orbital_name fzc, drc, docc, act, ras1, ras2, ras3, ras4, pop, vir, fzv, drv, or all
     * @return C            Returns the appropriate orbitals in the SO basis.
     */
    std::shared_ptr<Matrix> get_orbitals(const std::string &orbital_name);

    /**!
     * Similar to wavefunction.Ca_subset(); however, this version knows about all of the CI
     * subspaces in the SO basis. We stick to the MCSCF definitions for now.
     * @param  orbital_name FZC, DRC, DOCC, ACT, RAS1, RAS2, RAS3, RAS4, POP, VIR, FZV, DRV, or ALL
     * @param  orbitals     SharedMatrix to set
     */
    void set_orbitals(const std::string &orbital_name, SharedMatrix orbitals);

  protected:

    /// the v2RDMSolver 
    std::shared_ptr<v2RDMSolver> v2rdm;

    /// Find out which orbitals belong where
    void orbital_locations(const std::string &orbital_name, int *start, int *end);

};

}

#endif
