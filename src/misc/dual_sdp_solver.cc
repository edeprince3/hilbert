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

#include "dual_sdp_solver.h"

#include <misc/cg_solver.h>
#include <misc/omp.h>
#include <misc/blas.h>

#include <lbfgs.h>

using namespace psi;
using namespace fnocc;

namespace hilbert {

// CG callback function
static void evaluate_cg_AATu(SharedVector Ax, SharedVector x, void * data) {

    // reinterpret void * as an instance of DualSDPSolver
    DualSDPSolver* sdp = reinterpret_cast<DualSDPSolver*>(data);
    sdp->evaluate_AATu(Ax, x);

}

void DualSDPSolver::evaluate_AATu(std::shared_ptr<Vector> AATu,std::shared_ptr<Vector> u) {
    AATu->zero();
    std::shared_ptr<Vector> ATu (new Vector(n_primal_));
    evaluate_ATu_(ATu,u,data_);
    evaluate_Au_(AATu,ATu,data_);
}

// liblbfgs routines:
static lbfgsfloatval_t lbfgs_evaluate_z(void * instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step) {

    DualSDPSolver * sdp = reinterpret_cast<DualSDPSolver*>(instance);
    double f = sdp->evaluate_gradient_z(x,g);

    return f;
}

static lbfgsfloatval_t lbfgs_evaluate_x(void * instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step) {

    DualSDPSolver * sdp = reinterpret_cast<DualSDPSolver*>(instance);
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
    DualSDPSolver * sdp = reinterpret_cast<DualSDPSolver*>(instance);
    sdp->set_iiter(k);
    return 0;
}

DualSDPSolver::DualSDPSolver(long int n_primal, long int n_dual, Options & options)
    : SDPSolver(n_primal,n_dual,options){

    mu_primal_    = 0.1;
    mu_xz_        = 0.1;
    iiter_        = 0;

    tmp_    = (std::shared_ptr<Vector>)(new Vector(n_primal_));
    cg_rhs_ = (std::shared_ptr<Vector>)(new Vector(n_dual_));

    // lbfgs container auxiliary variables that define z
    lbfgs_vars_z_ = lbfgs_malloc(n_primal_);

    // seed auxiliary variables
    srand(0);
    for (int i = 0; i < n_primal_; i++) {
        lbfgs_vars_z_[i] = 2.0 * ( (double)rand()/RAND_MAX - 0.5 ) / 1000.0;
    }

}

DualSDPSolver::~DualSDPSolver(){

    free(lbfgs_vars_z_);

}

void DualSDPSolver::solve_low_rank(std::shared_ptr<Vector> x,
                                   std::shared_ptr<Vector> b,
                                   std::shared_ptr<Vector> c,
                                   std::vector<int> primal_block_dim,
                                   std::vector<int> primal_block_rank,
                                   int maxiter,
                                   SDPCallbackFunction evaluate_Au,
                                   SDPCallbackFunction evaluate_ATu,
                                   void * data){
    // copy block ranks
    primal_block_rank_.clear();
    for (int block = 0; block < primal_block_rank.size(); block++) {
        primal_block_rank_.push_back(primal_block_rank[block]);
    }

    // call standard solve
    solve(x,b,c,primal_block_dim,maxiter,evaluate_Au,evaluate_ATu,data);
}

void DualSDPSolver::solve(std::shared_ptr<Vector> x,   
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
    evaluate_Au_ = evaluate_Au;
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

    // build initial z
    build_z(lbfgs_vars_z_);

    // class pointers to input c,x,b
    c_ = c;
    x_ = x;
    b_ = b;

    // this function can be called many times. don't forget to reset penalty parameter
    mu_primal_ = 0.1;
    mu_        = 0.1;
    mu_xz_     = 0.1;
    lag_xz_    = 1.0;

    // cg solver
    std::shared_ptr<CGSolver> cg (new CGSolver(n_dual_));
    cg->set_max_iter(options_.get_int("CG_MAXITER"));
    double cg_convergence = options_.get_double("CG_CONVERGENCE");
    cg->set_convergence(cg_convergence);

    // initial energy   
    double energy =  y_->vector_dot(b_);

    // the iterations
    outfile->Printf("\n");
    outfile->Printf("    initial dual energy: %20.12lf\n",energy);
    outfile->Printf("\n");
    outfile->Printf("           oiter");
    outfile->Printf("       y iter");
    outfile->Printf("       x iter");
    outfile->Printf("       z iter");
    outfile->Printf("        <x|z>");
    outfile->Printf("         E(p)");
    outfile->Printf("         E(d)");
    outfile->Printf("           mu");
    outfile->Printf("       eps(p)\n");
    outfile->Printf("        eps(d)\n");

    int oiter_local = 0;

    double max_err_primal = -999;
    double max_err_dual   = -999;
    double max_err_xz     = -999;

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);

    // adjust default lbfgs parameters
    param.max_iterations = options_.get_int("MAXITER");
    param.epsilon        = 1e-8; //options_.get_double("LBFGS_CONVERGENCE");

    do {

        // solve for y, given mu, x, and z

        double start = omp_get_wtime();

        // evaluate mu * (b - Ax) for CG
        evaluate_Au_(Au_, x_, data_);
        Au_->subtract(b_);
        Au_->scale(-mu_);

        // evaluate A(c-z) ( but don't overwrite c or z! )
        ATu_->copy(c_.get());
        ATu_->subtract(z_);
        evaluate_Au_(cg_rhs_, ATu_, data);

        // add tau*mu*(b-Ax) to A(c-z) and put result in cg_rhs_
        cg_rhs_->add(Au_);

        // solve CG problem (solve dL/dy = 0)
        cg->solve(Au_,y_,cg_rhs_,evaluate_cg_AATu,(void*)this);

        int cg_iter = cg->total_iterations();

        double energy_dual = y_->vector_dot(b_);

        double end = omp_get_wtime();

        iiter_time_  += end - start;
        iiter_total_ += cg_iter;

        // solve for z, given x, y, and mu

        start = omp_get_wtime();

        // initial objective function value default parameters
        lbfgsfloatval_t lag_z = evaluate_gradient_z(lbfgs_vars_z_,ATu_->pointer());
        lbfgs(n_primal_,lbfgs_vars_z_,&lag_z,lbfgs_evaluate_z,monitor_lbfgs_progress,(void*)this,&param);
        int z_iter = iiter_;

        // build z = r.rT
        build_z(lbfgs_vars_z_);

        // evaluate dual error

        // evaluate || A^T y - c + z||
        evaluate_ATu_(ATu_, y_, data_);
        ATu_->add(z_);
        ATu_->subtract(c_);
        dual_error_ = ATu_->norm();

        end = omp_get_wtime();

        oiter_time_ += end - start;

        // update penalty parameter (dual) and lagrange multipliers (x)

        double * err_p = ATu_->pointer();
        double * x_p = x_->pointer();
        if ( oiter_ % 20 == 0 && oiter_ > 0 ) {

/*
            // solve for x, given y, z, and mu
            lbfgsfloatval_t lag_x = evaluate_gradient_x(lbfgs_vars_x,ATu_->pointer());
            lbfgs(n_primal_,lbfgs_vars_x,&lag_x,lbfgs_evaluate_x,monitor_lbfgs_progress,(void*)this,&param);
            x_iter = iiter_;

            // build x = r.rT
            build_x(lbfgs_vars_x);

            mu_primal_ *= 0.1;
*/

            // update lagrange multipliers, x
            for (int i = 0; i < n_primal_; i++) {
                x_p[i] += err_p[i] / mu_;
            }
            //mu_ *= 0.1;

            double overlap = x_->vector_dot(z_);
            lag_xz_  += overlap / mu_xz_;
            mu_xz_ *= 0.9;

        }

        // evaluate primal energy
        double energy_primal = x_->vector_dot(c_);

        // evaluate primal energy: ||Ax - b||
        evaluate_Au_(Au_,x_,data_);
        Au_->subtract(b_);
        primal_error_ = Au_->norm();

        outfile->Printf("    %12i %12i %12i %12.6lf %12.6lf %12.6lf %12.2le %12.3le %12.3le\n",
                    oiter_,cg_iter, z_iter, x_->vector_dot(z_),energy_primal,energy_dual,mu_,primal_error_,dual_error_);

        oiter_++;
        oiter_local++;

        if ( dual_error_ > r_convergence_  || fabs(energy - energy_dual) > e_convergence_ ) {
            is_converged_ = false;
        }else {
            is_converged_ = true;
        }

        // update energy
        energy = energy_dual;

        if ( oiter_local == maxiter ) break;

    }while( !is_converged_ );

}

double DualSDPSolver::evaluate_gradient_z(const lbfgsfloatval_t * r, lbfgsfloatval_t * g) {

    // L = || c - ATy + z ||^2 + (x.z)^2, fixed y and z

    // pointers
    double * r_p   = (double*)r;
    double * x_p   = x_->pointer();
    double * z_p   = z_->pointer();
    double * ATu_p = ATu_->pointer();
    double * tmp_p = tmp_->pointer();

    // build z = r.rT
    build_z(r_p);

    // evaluate || A^T y - c + z||
    evaluate_ATu_(ATu_, y_, data_);
    ATu_->subtract(c_);
    ATu_->add(z_);
    double dual_error = ATu_->norm();

    // z.x = 0
    double xz = x_->vector_dot(z_);

    // evaluate lagrangian
    //double lagrangian = dual_error*dual_error/(2.0*mu_) + overlap*overlap/(2.0*mu_xz_);
    //double lagrangian = dual_error*dual_error + overlap*overlap;
    double lagrangian = 0.0;

    //lagrangian -= x_->vector_dot(ATu_);
    //lagrangian += dual_error*dual_error/(2.0 * mu_);
    lagrangian += dual_error*dual_error;

    //lagrangian -= xz * lag_xz_;
    //lagrangian += xz*xz/(2.0 * mu_xz_);
    lagrangian += xz*xz;

    int off_nn = 0;
    for (int block = 0; block < primal_block_dim_.size(); block++) {
        int n = primal_block_dim_[block];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                tmp_p[i*n+j + off_nn]  = ( ATu_p[i*n+j + off_nn] + ATu_p[j*n+i + off_nn] );
                tmp_p[i*n+j + off_nn] += (  x_p[i*n+j + off_nn] +  x_p[j*n+i + off_nn] ) * xz;

                // ||ATy - c + z||^2
                //tmp_p[i*n+j + off_nn]  = ( ATu_p[i*n+j + off_nn] + ATu_p[j*n+i + off_nn] ) / (2.0 * mu_);

                // -x (ATy - c + z)
                //tmp_p[i*n+j + off_nn] -= 0.5 * (  x_p[i*n+j + off_nn] +  x_p[j*n+i + off_nn] );

                // (xT.z)^2 
                //tmp_p[i*n+j + off_nn] += (  x_p[i*n+j + off_nn] +  x_p[j*n+i + off_nn] ) * xz / (2.0 * mu_xz_);

                // -lag xT.z
                //tmp_p[i*n+j + off_nn] -= 0.5 * (  x_p[i*n+j + off_nn] +  x_p[j*n+i + off_nn] ) * lag_xz_ ;

            }
        }
        off_nn += n*n;
    }

    off_nn = 0;
    int off_nm = 0;
    for (int block = 0; block < primal_block_dim_.size(); block++) {
        int n = primal_block_dim_[block];
        int m = primal_block_rank_[block];
        if ( n == 0 ) continue;
        F_DGEMM('n', 'n', n, m, n, 2.0, tmp_p + off_nn, n, r_p + off_nm, n, 0.0, g + off_nm, n);
        off_nn += n*n;
        off_nm += n*m;
    }

    return lagrangian;
}

double DualSDPSolver::evaluate_gradient_x(const lbfgsfloatval_t * r, lbfgsfloatval_t * g) {

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
    double lagrangian = energy - y_->vector_dot(Au_) + nrm*nrm/(2.0 * mu_primal_);

    // dL/dR = 2( A^T [ 2/mu(Ax-b) - y] + c) . r

    // evaluate A^T (1/mu[Ax-b] - y)
    Au_->scale(1.0/mu_primal_);
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

// build x = r.rT
void DualSDPSolver::build_x(double * r){

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

// build z = r.rT
void DualSDPSolver::build_z(double * r){

    double * z_p = z_->pointer();

    int off_nn = 0;
    int off_nm = 0;
    for (int block = 0; block < primal_block_dim_.size(); block++) {
        int n = primal_block_dim_[block];
        int m = primal_block_rank_[block];
        if ( n == 0 ) continue;
        F_DGEMM('n', 't', n, n, m, 1.0, r + off_nm, n, r + off_nm, n, 0.0, z_p + off_nn, n);
        off_nm += n*m;
        off_nn += n*n;
    }

}


}

