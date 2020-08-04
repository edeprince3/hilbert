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

#include <lbfgs.h>

#include "rrsdp_solver.h"

#include <misc/omp.h>
#include <misc/blas.h>

using namespace psi;
using namespace fnocc;

namespace hilbert {

typedef void (*RRSDPCallbackFunction)(std::shared_ptr<Vector>,std::shared_ptr<Vector>,void *);

// liblbfgs routines:
static lbfgsfloatval_t lbfgs_evaluate(void * instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step) {

    RRSDPSolver * sdp = reinterpret_cast<RRSDPSolver*>(instance);
    double f = sdp->evaluate_gradient(x,g);

    return f;
}

static int monitor_lbfgs_progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    )
{
    RRSDPSolver * sdp = reinterpret_cast<RRSDPSolver*>(instance);
    sdp->set_iiter(k);
    return 0;
}

RRSDPSolver::RRSDPSolver(long int n_primal, long int n_dual, Options & options)
    :options_(options){

    n_primal_     = n_primal;
    n_dual_       = n_dual;
    mu_           = 0.1;
    primal_error_ = 0.0;
    oiter_        = 0;
    iiter_        = 0;

    y_    = (std::shared_ptr<Vector>)(new Vector(n_dual_));
    Au_   = (std::shared_ptr<Vector>)(new Vector(n_dual_));
    ATu_  = (std::shared_ptr<Vector>)(new Vector(n_primal_));

    e_convergence_ = options_.get_double("E_CONVERGENCE");
    r_convergence_ = options_.get_double("R_CONVERGENCE");

    is_converged_ = false;
}

RRSDPSolver::~RRSDPSolver(){
}


void RRSDPSolver::solve(std::shared_ptr<Vector> x,   
                        std::shared_ptr<Vector> b, 
                        std::shared_ptr<Vector> c,
                        std::vector<int> primal_block_dim,
                        int maxiter,
                        RRSDPCallbackFunction evaluate_Au, 
                        RRSDPCallbackFunction evaluate_ATu, 
                        void * data){

    // class pointer to input data
    data_ = data;

    // class pointers to callback functions
    evaluate_Au_ = evaluate_Au;
    evaluate_ATu_ = evaluate_ATu;

    // copy block sizes
    primal_block_dim_.clear();
    for (int block = 0; block < primal_block_dim.size(); block++) {
        primal_block_dim_.push_back(primal_block_dim[block]);
    }

    // class pointers to input c,x,b
    c_ = c;
    x_ = x;
    b_ = b;

    // lbfgs container for primal solution vector, x
    lbfgsfloatval_t * lbfgs_vars  = lbfgs_malloc(n_primal_);
    C_DCOPY(n_primal_,x_->pointer(),1,lbfgs_vars,1);

    // initial energy   
    double energy =  C_DDOT(n_primal_,c_->pointer(),1,x_->pointer(),1);

    // this function can be called many times. don't forget to reset penalty parameter
    mu_ = 0.1;

    // the iterations
    outfile->Printf("\n");
    outfile->Printf("    initial primal energy: %20.12lf\n",energy);
    outfile->Printf("\n");
    outfile->Printf("           oiter");
    outfile->Printf("        iiter");
    outfile->Printf("            L");
    outfile->Printf("            E");
    outfile->Printf("           mu");
    outfile->Printf("     ||Ax-b||\n");

    int oiter_local = 0;

    double max_err = -999;
    do {

        // minimize lagrangian

        // initial objective function value default parameters
        lbfgsfloatval_t lagrangian = evaluate_gradient(x_->pointer(),ATu_->pointer());
        lbfgs_parameter_t param;
        lbfgs_parameter_init(&param);

        // adjust default parameters
        param.max_iterations = options_.get_int("MAXITER");
        //param.epsilon        = options_.get_double("LBFGS_CONVERGENCE");

        lbfgs(n_primal_,lbfgs_vars,&lagrangian,lbfgs_evaluate,monitor_lbfgs_progress,(void*)this,&param);

        // update lagrange multipliers and penalty parameter

        // build x = r.rT
        double * x_p = x_->pointer();
        double * r_p = lbfgs_vars;
        int off = 0;
        for (int block = 0; block < primal_block_dim_.size(); block++) {
            int n = primal_block_dim_[block];
            if ( n == 0 ) continue;
            F_DGEMM('n', 't', n, n, n, 1.0, r_p + off, n, r_p + off, n, 0.0, x_p + off, n);
            off += n*n;
        }

        // evaluate x^T.c
        double energy_primal = C_DDOT(n_primal_,x_->pointer(),1,c_->pointer(),1);

        // evaluate (Ax-b)
        evaluate_Au_(Au_,x_,data_);
        Au_->subtract(b_);
        primal_error_ = Au_->norm();

        double new_max_err = 0.0;
        int imax = -1;
        double * Au_p = Au_->pointer();
        double * y_p  = y_->pointer();
        for (int i = 0; i < n_dual_; i++) {
            if ( fabs(Au_p[i]) > new_max_err ) {
                new_max_err = fabs(Au_p[i]);
                imax = i;
             }
        }
        if ( new_max_err < 0.25 * max_err ){
            for (int i = 0; i < n_dual_; i++) {
                y_p[i] -= Au_p[i] / mu_;
            }
        }else{
            mu_ *= 0.1;
        }
        max_err = new_max_err;

        outfile->Printf("    %12i %12i %12.6lf %12.6lf %12.2le %12.3le\n",
                    oiter_,iiter_,lagrangian,energy_primal,mu_,primal_error_);

        iiter_total_ += iiter_;

        oiter_++;
        oiter_local++;

        if ( primal_error_ > r_convergence_  || fabs(energy - energy_primal) > e_convergence_ ) {
            is_converged_ = false;
        }else {
            is_converged_ = true;
        }

        // update energy
        energy = energy_primal;

        if ( oiter_local == maxiter ) break;

    }while( !is_converged_ );

    free(lbfgs_vars);

}

double RRSDPSolver::evaluate_gradient(const lbfgsfloatval_t * r, lbfgsfloatval_t * g) {

    // L = x^Tc - y^T (Ax - b) + 1/mu || Ax - b ||

    // pointers
    double * r_p   = (double*)r;
    double * x_p   = x_->pointer();
    double * c_p   = c_->pointer();

    // build x = r.rT
    int off = 0;
    for (int block = 0; block < primal_block_dim_.size(); block++) {
        int n = primal_block_dim_[block];
        if ( n == 0 ) continue;
        F_DGEMM('n', 't', n, n, n, 1.0, r_p + off, n, r_p + off, n, 0.0, x_p + off, n);
        off += n*n;
    }

    // evaluate primal energy
    double energy = C_DDOT(n_primal_,x_p,1,c_p,1);

    // evaluate (Ax-b)
    evaluate_Au_(Au_,x_,data_);
    Au_->subtract(b_);

    // evaluate sqrt(||Ax-b||)
    double nrm = Au_->norm();

    // evaluate lagrangian
    double lagrangian = energy - y_->vector_dot(Au_) + nrm*nrm/mu_;

    // dL/dR = 2( A^T [ 2/mu(Ax-b) - y] + c) . r

    // evaluate A^T (2/mu[Ax-b] - y)
    Au_->scale(2.0/mu_);
    Au_->subtract(y_);
    evaluate_ATu_(ATu_,Au_,data_);

    // add integrals for derivative of energy
    ATu_->add(c_);

    // evaluate gradient of lagrangian
    off = 0;
    for (int block = 0; block < primal_block_dim_.size(); block++) {
        int n = primal_block_dim_[block];
        if ( n == 0 ) continue;
        F_DGEMM('n', 'n', n, n, n, 2.0, ATu_->pointer() + off, n, r_p + off, n, 0.0, g + off, n);
        off += n*n;
    }

    return lagrangian;
}


}

