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

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <psi4/libplugin/plugin.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libqt/qt.h>

#include "cg_solver.h"

using namespace psi;

namespace hilbert{

CGSolver::CGSolver(long int n) {
    n_              = n;
    iter_           = 0;
    cg_max_iter_    = 10000;
    cg_convergence_ = 1e-9;
    p = SharedVector(new Vector(n));
    r = SharedVector(new Vector(n));
    //z = SharedVector(new Vector(n));
}
CGSolver::~CGSolver(){
}
void CGSolver::set_max_iter(int iter) {
    cg_max_iter_ = iter;
}
void CGSolver::set_convergence(double conv) {
    cg_convergence_ = conv;
}

void CGSolver::preconditioned_solve(long int n,
                    SharedVector Ap, 
                    SharedVector  x, 
                    SharedVector  b, 
                    SharedVector  precon, 
                    CallbackType function, void * data) {

    if ( n != n_ ) {
        throw PsiException("Warning: dimension does not match dimension from initialization",__FILE__,__LINE__);
    }

    double * p_p = p->pointer();
    double * r_p = r->pointer();
    double * z_p = z->pointer();

    double alpha = 0.0;
    double beta  = 0.0;

    // call some function to evaluate A.x.  Result in Ap
    function(n,Ap,x,data);

    double * b_p      = b->pointer();
    double * x_p      = x->pointer();
    double * Ap_p     = Ap->pointer();
    double * precon_p = precon->pointer();

    for (int i = 0; i < n; i++) {
        r_p[i] = b_p[i] - Ap_p[i];
        z_p[i] = precon_p[i] * r_p[i];
    }
    C_DCOPY(n,z_p,1,p_p,1);

    iter_ = 0;
    do {

        // call some function to evaluate A.p.  Result in Ap
        function(n,Ap,p,data);

        double rz  = C_DDOT(n_,r_p,1,z_p,1);
        double pap = C_DDOT(n_,p_p,1,Ap_p,1);
        double alpha = rz / pap;
        C_DAXPY(n_,alpha,p_p,1,x_p,1);
        C_DAXPY(n_,-alpha,Ap_p,1,r_p,1);

        // if r is sufficiently small, then exit loop
        double rrnew = C_DDOT(n_,r_p,1,r_p,1);
        double nrm = sqrt(rrnew);// / sqrt(n_);
        if ( nrm < cg_convergence_ ) break;

        for (int i = 0; i < n; i++) {
            z_p[i] = precon_p[i] * r_p[i];
        }
        double rznew  = C_DDOT(n_,r_p,1,z_p,1);
        double beta = rznew/rz;

        C_DSCAL(n_,beta,p_p,1);
        C_DAXPY(n_,1.0,z_p,1,p_p,1);

        iter_++;

    }while(iter_ < cg_max_iter_ );
}

void CGSolver::solve(long int n,
                    SharedVector Ap, 
                    SharedVector  x, 
                    SharedVector  b, 
                    CallbackType function, void * data) {

    if ( n != n_ ) {
        throw PsiException("Warning: dimension does not match dimension from initialization",__FILE__,__LINE__);
    }

    double * p_p = p->pointer();
    double * r_p = r->pointer();

    double alpha = 0.0;
    double beta  = 0.0;

    // call some function to evaluate A.x.  Result in Ap
    function(n,Ap,x,data);

    double * b_p  = b->pointer();
    double * x_p  = x->pointer();
    double * Ap_p = Ap->pointer();

    for (int i = 0; i < n; i++) {
        r_p[i] = b_p[i] - Ap_p[i];
    }

    C_DCOPY(n,r_p,1,p_p,1);

    iter_ = 0;
    do {

        // call some function to evaluate A.p.  Result in Ap
        function(n,Ap,p,data);

        double rr  = C_DDOT(n_,r_p,1,r_p,1);
        double pap = C_DDOT(n_,p_p,1,Ap_p,1);
        double alpha = rr / pap;
        C_DAXPY(n_,alpha,p_p,1,x_p,1);
        C_DAXPY(n_,-alpha,Ap_p,1,r_p,1);

        // if r is sufficiently small, then exit loop
        double rrnew = C_DDOT(n_,r_p,1,r_p,1);
        double nrm = sqrt(rrnew);// / sqrt(n_);
        double beta = rrnew/rr;
        if ( nrm < cg_convergence_ ) break;

        C_DSCAL(n_,beta,p_p,1);
        C_DAXPY(n_,1.0,r_p,1,p_p,1);

        iter_++;

    }while(iter_ < cg_max_iter_ );
}

int CGSolver::total_iterations() {
    return iter_;
}


}// end of namespace
