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

#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libqt/qt.h>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libtrans/mospace.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/mintshelper.h>
#include <psi4/libmints/x2cint.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/matrix.h>

#include <psi4/psifiles.h>

#include "v2rdm_solver.h"

#include <misc/omp.h>

using namespace psi;

namespace hilbert{ 

SharedMatrix v2RDMSolver::GetOEI() {

    std::shared_ptr<MintsHelper> mints(new MintsHelper(reference_wavefunction_));
    T_ = (std::shared_ptr<Matrix> ) (new Matrix(mints->so_kinetic()));
    V_ = (std::shared_ptr<Matrix> ) (new Matrix(mints->so_potential()));

    //if (options_.get_str("RELATIVISTIC") == "X2C") {
    //    X2CInt x2cint;
    //    std::shared_ptr<BasisSet> basisset = reference_wavefunction_->get_basisset("ORBITAL"); 
    //    std::shared_ptr<BasisSet> rel_basisset = reference_wavefunction_->get_basisset("BASIS_RELATIVISTIC");
    //    x2cint.compute(basisset, rel_basisset, S_, T_, V_); 
    //}

    std::shared_ptr<Matrix> K1 (new Matrix(T_));
    K1->add(V_);
    K1->transform(Ca_);

    return K1;
}

}
