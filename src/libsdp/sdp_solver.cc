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

#include <math.h>

#include <psi4/libmints/vector.h>
#include <psi4/libmints/matrix.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libqt/qt.h>
#include <psi4/libpsi4util/PsiOutStream.h>

#include "sdp_solver.h"

using namespace psi;

namespace hilbert {

SDPSolver::SDPSolver(long int n_primal, long int n_dual, Options & options)
    :options_(options){

    n_primal_     = n_primal;
    n_dual_       = n_dual;
    mu_           = 0.1;
    primal_error_ = 0.0;
    dual_error_   = 0.0;
    oiter_        = 0;
    iiter_total_  = 0;
    oiter_time_   = 0.0;
    iiter_time_   = 0.0;

    y_   = (std::shared_ptr<Vector>)(new Vector(n_dual_));
    Au_  = (std::shared_ptr<Vector>)(new Vector(n_dual_));
    z_   = (std::shared_ptr<Vector>)(new Vector(n_primal_));
    ATu_ = (std::shared_ptr<Vector>)(new Vector(n_primal_));

    e_convergence_ = options_.get_double("E_CONVERGENCE");
    r_convergence_ = options_.get_double("R_CONVERGENCE");

    is_converged_ = false;
}

SDPSolver::~SDPSolver(){
}

}

