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

#include "lbfgs_helper.h"

#include <psi4/libqt/qt.h>
#include <psi4/libpsi4util/PsiOutStream.h>

using namespace psi;

namespace hilbert {

void lbfgs_error_check(int value) {

    //L-BFGS reaches convergence.
    //LBFGS_SUCCESS = 0,
    //LBFGS_CONVERGENCE = 0,
    //LBFGS_STOP,

    if (  value == 0 ) return;
    
    outfile->Printf("\n");
    outfile->Printf("    ==> WARNING <==\n");
    outfile->Printf("\n");
    outfile->Printf("    L-BFGS exited with an error:\n");
    outfile->Printf("\n");

    if ( value == (int)LBFGS_ALREADY_MINIMIZED) 
        outfile->Printf("        The initial variables already minimize the objective function.\n");
    if ( value == (int)LBFGSERR_UNKNOWNERROR) 
        outfile->Printf("        Unknown error.\n");
    if ( value == (int)LBFGSERR_LOGICERROR) 
        outfile->Printf("        Logic error.\n");
    if ( value == (int)LBFGSERR_OUTOFMEMORY) 
        outfile->Printf("        Insufficient memory.\n");
    if ( value == (int)LBFGSERR_CANCELED) 
        outfile->Printf("        The minimization process has been canceled.\n");
    if ( value == (int)LBFGSERR_INVALID_N) 
        outfile->Printf("        Invalid number of variables specified.\n");
    if ( value == (int)LBFGSERR_INVALID_N_SSE) 
        outfile->Printf("        Invalid number of variables (for SSE) specified.\n");
    if ( value == (int)LBFGSERR_INVALID_X_SSE) 
        outfile->Printf("        The array x must be aligned to 16 (for SSE).\n");
    if ( value == (int)LBFGSERR_INVALID_EPSILON) 
        outfile->Printf("        Invalid parameter lbfgs_parameter_t::epsilon specified.\n");
    if ( value == (int)LBFGSERR_INVALID_TESTPERIOD) 
        outfile->Printf("        Invalid parameter lbfgs_parameter_t::past specified.\n");
    if ( value == (int)LBFGSERR_INVALID_DELTA) 
        outfile->Printf("        Invalid parameter lbfgs_parameter_t::delta specified.\n");
    if ( value == (int)LBFGSERR_INVALID_LINESEARCH) 
        outfile->Printf("        Invalid parameter lbfgs_parameter_t::linesearch specified.\n");
    if ( value == (int)LBFGSERR_INVALID_MINSTEP) 
        outfile->Printf("        Invalid parameter lbfgs_parameter_t::max_step specified.\n");
    if ( value == (int)LBFGSERR_INVALID_MAXSTEP) 
        outfile->Printf("        Invalid parameter lbfgs_parameter_t::max_step specified.\n");
    if ( value == (int)LBFGSERR_INVALID_FTOL) 
        outfile->Printf("        Invalid parameter lbfgs_parameter_t::ftol specified.\n");
    if ( value == (int)LBFGSERR_INVALID_WOLFE) 
        outfile->Printf("        Invalid parameter lbfgs_parameter_t::wolfe specified.\n");
    if ( value == (int)LBFGSERR_INVALID_GTOL) 
        outfile->Printf("        Invalid parameter lbfgs_parameter_t::gtol specified.\n");
    if ( value == (int)LBFGSERR_INVALID_XTOL) 
        outfile->Printf("        Invalid parameter lbfgs_parameter_t::xtol specified.\n");
    if ( value == (int)LBFGSERR_INVALID_MAXLINESEARCH)
        outfile->Printf("        Invalid parameter lbfgs_parameter_t::max_linesearch specified.\n");
    if ( value == (int)LBFGSERR_INVALID_ORTHANTWISE)
        outfile->Printf("        Invalid parameter lbfgs_parameter_t::orthantwise_c specified.\n");
    if ( value == (int)LBFGSERR_INVALID_ORTHANTWISE_START)
        outfile->Printf("        Invalid parameter lbfgs_parameter_t::orthantwise_start specified.\n");
    if ( value == (int)LBFGSERR_INVALID_ORTHANTWISE_END)
        outfile->Printf("        Invalid parameter lbfgs_parameter_t::orthantwise_end specified.\n");
    if ( value == (int)LBFGSERR_OUTOFINTERVAL)
        outfile->Printf("        The line-search step went out of the interval of uncertainty.\n");
    if ( value == (int)LBFGSERR_INCORRECT_TMINMAX)
        outfile->Printf("        A logic error occurred; alternatively, the interval of uncertainty became too small.\n");
    if ( value == (int)LBFGSERR_ROUNDING_ERROR)
        outfile->Printf("        A rounding error occurred; alternatively, no line-search step satisfies the sufficient decrease and curvature conditions.\n");
    if ( value == (int)LBFGSERR_MINIMUMSTEP)
        outfile->Printf("        The line-search step became smaller than lbfgs_parameter_t::min_step.\n");
    if ( value == (int)LBFGSERR_MAXIMUMSTEP)
        outfile->Printf("        The line-search step became larger than lbfgs_parameter_t::max_step.\n");
    if ( value == (int)LBFGSERR_MAXIMUMLINESEARCH)
        outfile->Printf("        The line-search routine reaches the maximum number of evaluations.\n");
    if ( value == (int)LBFGSERR_MAXIMUMITERATION)
        outfile->Printf("        The algorithm routine reaches the maximum number of iterations.\n");
    if ( value == (int)LBFGSERR_WIDTHTOOSMALL)
        outfile->Printf("        Relative width of the interval of uncertainty is at most lbfgs_parameter_t::xtol.\n");
    if ( value == (int)LBFGSERR_INVALIDPARAMETERS)
        outfile->Printf("        A logic error (negative line-search step) occurred.\n");
    if ( value == (int)LBFGSERR_INCREASEGRADIENT)
        outfile->Printf("        The current search direction increases the objective function value.\n");

    outfile->Printf("\n");

    //exit(0);

}

}
