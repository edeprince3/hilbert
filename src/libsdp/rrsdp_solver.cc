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

#include "rrsdp_solver.h"

#include <misc/omp.h>
#include <misc/blas.h>
#include <misc/lbfgs_helper.h>

using namespace psi;
using namespace fnocc;

namespace hilbert {

// liblbfgs routines:
static lbfgsfloatval_t lbfgs_evaluate(void * instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step) {

    RRSDPSolver * sdp = reinterpret_cast<RRSDPSolver*>(instance);
    double f = sdp->evaluate_gradient_x(x,g);

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
    : SDPSolver(n_primal,n_dual,options) {

    iiter_ = 0;

    // lbfgs container auxiliary variables that define x
    lbfgs_vars_x_ = lbfgs_malloc(n_primal_);

    // seed auxiliary variables
    srand(0);
    for (int i = 0; i < n_primal_; i++) {
        lbfgs_vars_x_[i] = 2.0 * ( (double)rand()/RAND_MAX - 0.5 ) / 1000.0;
    }

    mu_ = 0.1;
    mu_reset_ = true;
    mu_scale_factor_ = 0.1;
}

RRSDPSolver::~RRSDPSolver(){

    free(lbfgs_vars_x_);

}

void RRSDPSolver::solve(std::shared_ptr<Vector> x,   
                        std::shared_ptr<Vector> b, 
                        std::shared_ptr<Vector> c,
                        std::vector<int> primal_block_dim,
                        int maxiter,
                        SDPCallbackFunction evaluate_Au, 
                        SDPCallbackFunction evaluate_ATu, 
                        void * data){

    // class pointer to input data
    data_ = data;

    // class pointers to callback functions
    evaluate_Au_  = evaluate_Au;
    evaluate_ATu_ = evaluate_ATu;

    // copy block sizes
    primal_block_dim_.clear();
    for (int block = 0; block < primal_block_dim.size(); block++) {
        primal_block_dim_.push_back(primal_block_dim[block]);
    }
    // if block ranks are not set, assign rank = dim
    if ( primal_block_rank_.size() != primal_block_dim_.size() ) {
        primal_block_rank_.clear();
        for (int block = 0; block < primal_block_dim.size(); block++) {
            primal_block_rank_.push_back(primal_block_dim[block]);
        }
    }

    // class pointers to input c,x,b
    c_ = c;
    x_ = x;
    b_ = b;

    build_x(lbfgs_vars_x_);

    // initial energy   
    double energy =  x_->vector_dot(c_);

    // this function can be called many times. don't forget to reset penalty parameter
    if ( mu_reset_ ) {
        mu_ = 0.1;
    }

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

    std::shared_ptr<Vector> tmp (new Vector(n_primal_));

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    param.max_iterations = options_.get_int("MAXITER");
    //param.epsilon        = options_.get_double("LBFGS_CONVERGENCE");

    do {

        // minimize lagrangian

        // initial objective function value default parameters
        lbfgsfloatval_t lagrangian = evaluate_gradient_x(lbfgs_vars_x_,tmp->pointer());

        if (oiter_ == 0) {
            param.epsilon = 0.01;
        }else {
            param.epsilon = 0.01 * primal_error_;
        }
        if ( param.epsilon < r_convergence_ ) {
            param.epsilon = r_convergence_;
        }
        int status = lbfgs(n_primal_,lbfgs_vars_x_,&lagrangian,lbfgs_evaluate,monitor_lbfgs_progress,(void*)this,&param);
        lbfgs_error_check(status);

        // update lagrange multipliers and penalty parameter

        // build x = r.rT
        build_x(lbfgs_vars_x_);

        // evaluate x^T.c
        double energy_primal = x_->vector_dot(c_);

        // evaluate (Ax-b)
        evaluate_Au_(Au_,x_,data_);
        Au_->subtract(b_);
        primal_error_ = Au_->norm();

        double new_max_err = 0.0;
        int imax = -1;
        double * Au_p = Au_->pointer();
        double * y_p  = y_->pointer();
        //for (int i = 0; i < n_dual_; i++) {
        //    if ( fabs(Au_p[i]) > new_max_err ) {
        //        new_max_err = fabs(Au_p[i]);
        //        imax = i;
        //     }
        //}
        //if ( new_max_err < 0.25 * max_err ){
            for (int i = 0; i < n_dual_; i++) {
                y_p[i] -= Au_p[i] / mu_;
            }
        //}else{
            mu_ *= mu_scale_factor_;
        //}
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

}

// build x = r.rT
void RRSDPSolver::build_x(double * r){

    double * x_p = x_->pointer();

    int off_nn = 0;
    int off_nm = 0;
    for (int block = 0; block < primal_block_dim_.size(); block++) {
        int n = primal_block_dim_[block];
        int m = primal_block_rank_[block];
        if ( n == 0 ) continue;
        F_DGEMM('n', 't', n, n, m, 1.0, r + off_nm, n, r + off_nm, n, 0.0, x_p + off_nn, n);
        off_nm += n*m;
        off_nn += n*n;
    }

}

double RRSDPSolver::evaluate_gradient_x(const lbfgsfloatval_t * r, lbfgsfloatval_t * g) {

    // L = x^Tc - y^T (Ax - b) + 1/mu || Ax - b ||

    // pointers
    double * r_p   = (double*)r;
    double * x_p   = x_->pointer();
    double * c_p   = c_->pointer();

    // build x = r.rT
    build_x(r_p);

    // evaluate primal energy
    double energy = x_->vector_dot(c_);

    // evaluate (Ax-b)
    evaluate_Au_(Au_,x_,data_);
    Au_->subtract(b_);

    // evaluate sqrt(||Ax-b||)
    double nrm = Au_->norm();

    // evaluate lagrangian
    double lagrangian = energy - y_->vector_dot(Au_) + nrm*nrm/(2.0 * mu_);

    // dL/dR = 2( A^T [ 2/mu(Ax-b) - y] + c) . r

    // evaluate A^T (1/mu[Ax-b] - y)
    Au_->scale(1.0/mu_);
    Au_->subtract(y_);
    evaluate_ATu_(ATu_,Au_,data_);

    // add integrals for derivative of energy
    ATu_->add(c_);

    double * ATu_p = ATu_->pointer();
    int off_nn = 0;
    for (int block = 0; block < primal_block_dim_.size(); block++) {
        int n = primal_block_dim_[block];
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                double dum_ij = ATu_p[i*n+j + off_nn];
                double dum_ji = ATu_p[j*n+i + off_nn];
                ATu_p[i*n+j + off_nn] = dum_ij + dum_ji;
                ATu_p[j*n+i + off_nn] = dum_ij + dum_ji;
            }
        }
        off_nn += n*n;
    }

    // evaluate gradient of lagrangian
    off_nn = 0;
    int off_nm = 0;
    for (int block = 0; block < primal_block_dim_.size(); block++) {
        int n = primal_block_dim_[block];
        int m = primal_block_rank_[block];
        if ( n == 0 ) continue;
        F_DGEMM('n', 'n', n, m, n, 1.0, ATu_->pointer() + off_nn, n, r_p + off_nm, n, 0.0, g + off_nm, n);
        off_nn += n*n;
        off_nm += n*m;
    }

    return lagrangian;
}


}

