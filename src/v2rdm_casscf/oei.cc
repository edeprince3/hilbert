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

    std::shared_ptr<Matrix> K1 (new Matrix(T_));
    K1->add(V_);
    K1->transform(Ca_);

    return K1;
}

SharedMatrix v2RDMSolver::GetOEI_hubbard() {

    std::shared_ptr<Matrix> h (new Matrix(amo_,amo_));
    double ** h_p = h->pointer();
    double t = options_.get_double("HUBBARD_T");
    for (int i = 0; i < amo_; i++) {
        for (int j = 0; j < amo_; j++) {
            if ( abs(i-j) == 1 ) {
                h_p[i][j] = -t;
            }else if ( abs(i-j) == amo_ - 1 ) {
                h_p[i][j] = -t;
            }else {
                h_p[i][j] = 0.0;
            }
        }
    }

    return h;

}

}
