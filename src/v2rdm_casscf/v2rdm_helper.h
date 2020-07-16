/*
 *@BEGIN LICENSE
 *
 * v2RDM-CASSCF, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * Copyright (c) 2014, The Florida State University. All rights reserved.
 *
 *@END LICENSE
 *
 */

#ifndef V2RDM_HELPER_H
#define V2RDM_HELPER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <psi4/libiwl/iwl.h>
#include <psi4/libplugin/plugin.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libqt/qt.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>

#include "v2rdm_solver.h"

namespace hilbert{

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

