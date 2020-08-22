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

#include <misc/lbfgs_helper.h>
#include <misc/diis.h>

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
    tmp_->zero();
    evaluate_ATu_(tmp_,u,data_);
    evaluate_Au_(AATu,tmp_,data_);
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

    mu_primal_ = 0.1;
    mu_xz_     = 0.1;
    mu_        = 0.1;
    iiter_     = 0;

    xz_     = (std::shared_ptr<Vector>)(new Vector(n_primal_));
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

    // build initial z
    build_z(lbfgs_vars_z_);

    // class pointers to input c,x,b
    c_ = c;
    x_ = x;
    b_ = b;

    // this function can be called many times. don't forget to reset penalty parameter
    mu_primal_ = 0.1;
    mu_        = 0.1;
    //mu_xz_     = 0.1;
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
    outfile->Printf("       z iter");
    outfile->Printf("        <x|z>");
    outfile->Printf("         E(p)");
    outfile->Printf("         E(d)");
    outfile->Printf("           mu");
    outfile->Printf("       eps(p)");
    outfile->Printf("        eps(d)\n");

    int oiter_local = 0;

    std::shared_ptr<Vector> dy (new Vector(n_dual_));
    std::shared_ptr<Vector> dz (new Vector(n_primal_));
    std::shared_ptr<Vector> solution (new Vector(n_dual_ + 2 * n_primal_));
    std::shared_ptr<Vector> error    (new Vector(n_dual_ + 2 * n_primal_));
    std::shared_ptr<DIIS> diis (new DIIS(n_dual_ + 2 * n_primal_));

    double energy_primal = x_->vector_dot(c_);

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

        double cg_conv_i = options_.get_double("CG_CONVERGENCE");
        if (oiter_ == 0) {
            cg_conv_i = 0.01;
        }else {
            cg_conv_i = 0.01 * dual_error_;
        }
        if (cg_conv_i < cg_convergence) {
            cg_conv_i = cg_convergence;
        }
        cg->set_convergence(cg_conv_i);

        dy->copy(y_.get());

        // solve CG problem (solve dL/dy = 0)
        cg->solve(Au_,y_,cg_rhs_,evaluate_cg_AATu,(void*)this);

        dy->subtract(y_);
        dy->scale(-1.0);

        int cg_iter = cg->total_iterations();

        double energy_dual = y_->vector_dot(b_);

        double end = omp_get_wtime();

        iiter_time_  += end - start;
        iiter_total_ += cg_iter;

        // solve for z, given x, y, and mu

        start = omp_get_wtime();

        // initial objective function value default parameters

        // evaluate A^T y - c. don't overwrite ... ATu_ is used by lbfgs
        evaluate_ATu_(ATu_, y_, data_);
        ATu_->subtract(c_);

        // update z = r.rT by min ||A^Ty-c+z||^2 + |x^T z|^2

        lbfgs_parameter_t param;
        lbfgs_parameter_init(&param);

        // adjust default lbfgs parameters
        param.max_iterations = 100;//options_.get_int("MAXITER");

        if (oiter_ == 0) {
            param.epsilon = 0.01;
        }else {
            param.epsilon = 0.01 * dual_error_;
        }
        if ( param.epsilon < 0.01*r_convergence_ ) {
            param.epsilon = 0.01*r_convergence_;
        }

        C_DCOPY(n_primal_,(double*)lbfgs_vars_z_,1,dz->pointer(),1);

        lbfgsfloatval_t lag_z = evaluate_gradient_z(lbfgs_vars_z_,tmp_->pointer());
        int status = lbfgs(n_primal_,lbfgs_vars_z_,&lag_z,lbfgs_evaluate_z,monitor_lbfgs_progress,(void*)this,&param);
        lbfgs_error_check(status);
        int z_iter = iiter_;

        C_DAXPY(n_primal_,-1.0,(double*)lbfgs_vars_z_,1,dz->pointer(),1);
        dz->scale(-1.0);

        // build z = r.rT ... should have been built in lbfgs
        //build_z(lbfgs_vars_z_);

        // evaluate dual error

        // evaluate || A^T y - c + z||  ... A^Ty - c should already be precomputed
        //evaluate_ATu_(ATu_, y_, data_);
        //ATu_->subtract(c_);
        ATu_->add(z_);
        dual_error_ = ATu_->norm();

        end = omp_get_wtime();

        oiter_time_ += end - start;

        // update penalty parameter (dual) and lagrange multipliers (x)

        double * ATu_p = ATu_->pointer();
        double * x_p = x_->pointer();

        if ( oiter_ % 1 == 0 && oiter_ > 0 ) {

            // update lagrange multipliers, x
            ATu_->scale(1.0/mu_);
            x_->add(ATu_);

            //mu_ *= 0.98;

            // evaluate primal energy
            energy_primal = x_->vector_dot(c_);

            // evaluate primal energy: ||Ax - b||
            evaluate_Au_(Au_,x_,data_);
            Au_->subtract(b_);
            primal_error_ = Au_->norm();

        }

        //if ( oiter_ % 20 == 0 && oiter_ > 0 ) {
        //    mu_ *= primal_error_ / dual_error_;
        //}

        //build_xz(lbfgs_vars_z_);
        //double xz = xz_->norm();

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

double DualSDPSolver::evaluate_gradient_z_new(const lbfgsfloatval_t * r, lbfgsfloatval_t * g) {

    double lagrangian = 0.0;

    // L = || ATy - c + z ||^2 + (x.z)^2, fixed y and z

    // pointers
    double * r_p   = (double*)r;
    double * x_p   = x_->pointer();
    double * z_p   = z_->pointer();
    double * ATu_p = ATu_->pointer();
    double * tmp_p = tmp_->pointer();

    // build z = r.rT
    build_z(r_p);

    // gradient contribution from ||ATy - c + z||

    // evaluate || A^T y - c + z|| ... A^t - c has been precomputed in ATu
    ATu_->add(z_);
    double dual_error = ATu_->norm();

    lagrangian += dual_error*dual_error;

    #pragma omp parallel for schedule (dynamic)
    for (int block = 0; block < primal_block_dim_.size(); block++) {
        int n = primal_block_dim_[block];
        if ( n == 0 ) continue;
        int off_nn = 0;
        for (int myblock = 0; myblock < block; myblock++) {
            int myn = primal_block_dim_[myblock];
            off_nn += myn*myn;
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                tmp_p[i*n+j + off_nn]  = ( ATu_p[i*n+j + off_nn] + ATu_p[j*n+i + off_nn] );
            }
        }
    }

    // subtract z from ATy - c + z so we don't need to construct ATy - c on the next iteration
    ATu_->subtract(z_);

    gradient_contribution(2.0,0.0,tmp_p, r_p, g);

/*
    // now contributions to gradient from sum [ (x.z)_ij ]^2 

    // build x.z = x.r.rT
    build_xz(r_p);
    double xz = xz_->norm();

    lagrangian += xz*xz;

    // build x.r
    double * xr = tmp_->pointer();
    build_xr(xr,x_p,r_p);

    // use dr = x.z.x.r
    gradient_contribution(2.0,1.0,xz_->pointer(), xr, g);

    // build z.r
    double * zr = tmp_->pointer();
    build_xr(zr,z_p,r_p);

    // build x.z.r
    double * xzr = xz_->pointer();
    build_xr(xzr,x_p,zr);

    // use dr = x.z.r.x
    gradient_contribution(2.0,1.0,x_p, xzr, g);
*/

    // now contributions to gradient from sum [ (x.z)_ij - mu dij ]^2 

    // build x.z = x.r.rT
    build_xz(r_p);

    double * xz_p = xz_->pointer();
    for (int block = 0; block < primal_block_dim_.size(); block++) {
        int n = primal_block_dim_[block];
        if ( n == 0 ) continue;
        int off_nn = 0;
        for (int myblock = 0; myblock < block; myblock++) {
            int myn = primal_block_dim_[myblock];
            off_nn += myn*myn;
        }
        for (int i = 0; i < n; i++) {
            xz_p[i*n+i + off_nn] -= mu_xz_;
        }
    }
    double xz = xz_->norm();
    lagrangian += xz * xz;

    // build x.r
    double * xr = tmp_->pointer();
    build_xr(xr,x_p,r_p);

    // use dr = (x.z - mu I).x.r
    gradient_contribution(2.0,1.0,xz_->pointer(), xr, g);

    // build z.r
    double * zr = tmp_->pointer();
    build_xr(zr,z_p,r_p);

    // build x.z.r
    double * xzr = xz_->pointer();
    build_xr(xzr,x_p,zr);

    C_DAXPY(n_primal_,-mu_xz_,r_p,1,xzr,1);

    // use dr = (x.z.r - mu r).x
    gradient_contribution(2.0,1.0,x_p, xzr, g);

    return lagrangian;
}

void DualSDPSolver::gradient_contribution(double alpha, double beta, double * mat1, double * mat2, double * gradient) {

    #pragma omp parallel for schedule (dynamic)
    for (int block = 0; block < primal_block_dim_.size(); block++) {
        int n = primal_block_dim_[block];
        int m = primal_block_rank_[block];
        if ( n == 0 ) continue;
        int off_nn = 0;
        int off_nm = 0;
        for (int myblock = 0; myblock < block; myblock++) {
            int myn = primal_block_dim_[myblock];
            int mym = primal_block_rank_[myblock];
            off_nn += myn*myn;
            off_nm += myn*mym;
        }
        F_DGEMM('n', 'n', n, m, n, alpha, mat1 + off_nn, n, mat2 + off_nm, n, beta, gradient + off_nm, n);
    }
}

int DualSDPSolver::update_z(double * r) {

    // min || ATy - c + z ||^2 + (x.z)^2, fixed y and z

    // pointers
    double * r_p   = (double*)r;
    double * x_p   = x_->pointer();
    double * z_p   = z_->pointer();
    double * ATu_p = ATu_->pointer();
    double * tmp_p = tmp_->pointer();

    std::shared_ptr<Vector> denom (new Vector(n_primal_));
    std::shared_ptr<Vector> gradient (new Vector(n_primal_));

    int iter = 0;
    do {

        // build z = r.rT
        build_z(r_p);

        // evaluate || A^T y - c + z|| ... A^t - c has been precomputed in ATu
        ATu_->add(z_);
        double dual_error = ATu_->norm();

        // x^T.z = 0
        double xz = x_->vector_dot(z_);

        double lagrangian = 0.0;

        lagrangian += dual_error*dual_error;
        lagrangian += xz*xz;

        // d (x^T.z)^2 / dr = (x^T.z) (x + x^T).r 
        // +
        // d ATy... term
        double * d_p = denom->pointer();
        #pragma omp parallel for schedule (dynamic)
        for (int block = 0; block < primal_block_dim_.size(); block++) {
            int n = primal_block_dim_[block];
            if ( n == 0 ) continue;
            int off_nn = 0;
            for (int myblock = 0; myblock < block; myblock++) {
                int myn = primal_block_dim_[myblock];
                off_nn += myn*myn;
            }
            C_DCOPY(n*n,ATu_p + off_nn,1,tmp_p + off_nn,1);
            C_DAXPY(n*n,xz,x_p + off_nn,1,tmp_p + off_nn,1);
            for (int i = 0; i < n; i++) {
                C_DAXPY(n, 1.0, ATu_p + off_nn + i, n, tmp_p + off_nn + i*n, 1);
                C_DAXPY(n, xz,    x_p + off_nn + i, n, tmp_p + off_nn + i*n, 1);
                //for (int j = 0; j < n; j++) {
                //    //tmp_p[i*n+j + off_nn]  = ( ATu_p[i*n+j + off_nn] + ATu_p[j*n+i + off_nn] );
                //    //tmp_p[i*n+j + off_nn] += (  x_p[i*n+j + off_nn] +  x_p[j*n+i + off_nn] ) * xz;
                //    tmp_p[i*n+j + off_nn] += ATu_p[j*n+i + off_nn];
                //    tmp_p[i*n+j + off_nn] += x_p[j*n+i + off_nn] * xz;
                //}
                d_p[off_nn + i * n + i ] = tmp_p[off_nn + i * n + i ];
                tmp_p[off_nn + i * n + i ] = 0.0;
            }
        }

        double * g_p = gradient->pointer();
        gradient_contribution(1.0,0.0,tmp_p, r_p, g_p);

        C_DCOPY(n_primal_,g_p,1,r_p,1);

        #pragma omp parallel for schedule (dynamic)
        for (int block = 0; block < primal_block_dim_.size(); block++) {
            int n = primal_block_dim_[block];
            int m = primal_block_rank_[block];
            if ( n == 0 ) continue;
            int off_nn = 0;
            int off_nm = 0;
            for (int myblock = 0; myblock < block; myblock++) {
                int myn = primal_block_dim_[myblock];
                int mym = primal_block_rank_[myblock];
                off_nn += myn*myn;
                off_nm += myn*mym;
            }
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    r_p[i*n+j + off_nm] *= - 1.0 / d_p[j*n+j + off_nn];
                }
            }
        }

        // subtract z from ATy - c + z so we don't need to construct ATy - c on the next iteration
        ATu_->subtract(z_);

        iter++;

    }while(iter < 10);

    return iter;
}

double DualSDPSolver::evaluate_gradient_z(const lbfgsfloatval_t * r, lbfgsfloatval_t * g) {

    // min || ATy - c + z ||^2 + (x.z)^2, fixed y and z

    // pointers
    double * r_p   = (double*)r;
    double * x_p   = x_->pointer();
    double * z_p   = z_->pointer();
    double * ATu_p = ATu_->pointer();
    double * tmp_p = tmp_->pointer();

    // build z = r.rT
    build_z(r_p);

    // evaluate || A^T y - c + z|| ... A^t - c has been precomputed in ATu
    ATu_->add(z_);
    double dual_error = ATu_->norm();

    // x^T.z = 0
    double xz = x_->vector_dot(z_);

    double lagrangian = 0.0;

    lagrangian += dual_error*dual_error;
    lagrangian += xz*xz;

    // d (x^T.z)^2 / dr = (x^T.z) (x + x^T).r 
    // +
    // d ATy... term
    #pragma omp parallel for schedule (dynamic)
    for (int block = 0; block < primal_block_dim_.size(); block++) {
        int n = primal_block_dim_[block];
        if ( n == 0 ) continue;
        int off_nn = 0;
        for (int myblock = 0; myblock < block; myblock++) {
            int myn = primal_block_dim_[myblock];
            off_nn += myn*myn;
        }
        C_DCOPY(n*n,ATu_p + off_nn,1,tmp_p + off_nn,1);
        C_DAXPY(n*n,xz,x_p + off_nn,1,tmp_p + off_nn,1);
        for (int i = 0; i < n; i++) {
            C_DAXPY(n, 1.0, ATu_p + off_nn + i, n, tmp_p + off_nn + i*n, 1);
            C_DAXPY(n, xz,    x_p + off_nn + i, n, tmp_p + off_nn + i*n, 1);
            //for (int j = 0; j < n; j++) {
            //    //tmp_p[i*n+j + off_nn]  = ( ATu_p[i*n+j + off_nn] + ATu_p[j*n+i + off_nn] );
            //    //tmp_p[i*n+j + off_nn] += (  x_p[i*n+j + off_nn] +  x_p[j*n+i + off_nn] ) * xz;
            //    tmp_p[i*n+j + off_nn] += ATu_p[j*n+i + off_nn];
            //    tmp_p[i*n+j + off_nn] += x_p[j*n+i + off_nn] * xz;
            //}
        }
    }

    gradient_contribution(2.0,0.0,tmp_p, r_p, g);

    // subtract z from ATy - c + z so we don't need to construct ATy - c on the next iteration
    ATu_->subtract(z_);

    return lagrangian;
}

// build z = r.rT
void DualSDPSolver::build_z(double * r){

    double * z_p = z_->pointer();

    #pragma omp parallel for schedule (dynamic)
    for (int block = 0; block < primal_block_dim_.size(); block++) {
        int n = primal_block_dim_[block];
        int m = primal_block_rank_[block];
        if ( n == 0 ) continue;
        int off_nn = 0;
        int off_nm = 0;
        for (int myblock = 0; myblock < block; myblock++) {
            int myn = primal_block_dim_[myblock];
            int mym = primal_block_rank_[myblock];
            off_nn += myn*myn;
            off_nm += myn*mym;
        }
        F_DGEMM('n', 't', n, n, m, 1.0, r + off_nm, n, r + off_nm, n, 0.0, z_p + off_nn, n);
    }

}

// build xz = x.(r.rT)
void DualSDPSolver::build_xz(double * r) {

    double * xr = tmp_->pointer();

    build_xr(xr, x_->pointer(), r);

    double * xz_p = xz_->pointer();

    #pragma omp parallel for schedule (dynamic)
    for (int block = 0; block < primal_block_dim_.size(); block++) {
        int n = primal_block_dim_[block];
        int m = primal_block_rank_[block];
        if ( n == 0 ) continue;
        int off_nn = 0;
        int off_nm = 0;
        for (int myblock = 0; myblock < block; myblock++) {
            int myn = primal_block_dim_[myblock];
            int mym = primal_block_rank_[myblock];
            off_nn += myn*myn;
            off_nm += myn*mym;
        }

        // (x.z)^i_j = sum_a (x.R)^a_i R^a_j
        F_DGEMM('n', 't', n, n, m, 1.0, r + off_nm, n, xr + off_nm, n, 0.0, xz_p + off_nn, n);
    }

}

// build xr = x.r
void DualSDPSolver::build_xr(double * xr, double * x, double * r) {

    #pragma omp parallel for schedule (dynamic)
    for (int block = 0; block < primal_block_dim_.size(); block++) {
        int n = primal_block_dim_[block];
        int m = primal_block_rank_[block];
        if ( n == 0 ) continue;
        int off_nn = 0;
        int off_nm = 0;
        for (int myblock = 0; myblock < block; myblock++) {
            int myn = primal_block_dim_[myblock];
            int mym = primal_block_rank_[myblock];
            off_nn += myn*myn;
            off_nm += myn*mym;
        }

        // (x.R)^a_i = sum_k x^i_k R^a_k
        F_DGEMM('t', 'n', n, m, n, 1.0, x + off_nn, n, r + off_nm, n, 0.0, xr + off_nm, n);
    }

}


}

