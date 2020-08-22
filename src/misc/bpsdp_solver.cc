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

#include "bpsdp_solver.h"

#include <misc/cg_solver.h>
#include <misc/omp.h>
#include <misc/blas.h>

using namespace psi;
using namespace fnocc;

namespace hilbert {

// CG callback function
static void evaluate_cg_AATu(SharedVector Ax, SharedVector x, void * data) {

    // reinterpret void * as an instance of BPSDPSolver
    BPSDPSolver* sdp = reinterpret_cast<BPSDPSolver*>(data);
    sdp->evaluate_AATu(Ax, x);

}

void BPSDPSolver::evaluate_AATu(std::shared_ptr<Vector> AATu,std::shared_ptr<Vector> u) {
    AATu->zero();
    std::shared_ptr<Vector> ATu (new Vector(n_primal_));
    evaluate_ATu_(ATu,u,data_);
    evaluate_Au_(AATu,ATu,data_);
}

/*
// liblbfgs routines:
static lbfgsfloatval_t lbfgs_evaluate_z(void * instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step) {

    BPSDPSolver * sdp = reinterpret_cast<BPSDPSolver*>(instance);
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
    BPSDPSolver * sdp = reinterpret_cast<BPSDPSolver*>(instance);
    sdp->set_lbfgs_iter(k);
    return 0;
}
*/

BPSDPSolver::BPSDPSolver(long int n_primal, long int n_dual, Options & options)
    : SDPSolver(n_primal,n_dual,options) {

    cg_rhs_ = (std::shared_ptr<Vector>)(new Vector(n_dual_));
}

BPSDPSolver::~BPSDPSolver(){
}


void BPSDPSolver::solve(std::shared_ptr<Vector> x,   
                        std::shared_ptr<Vector> b, 
                        std::shared_ptr<Vector> c,
                        std::vector<int> primal_block_dim,
                        int maxiter,
                        SDPCallbackFunction evaluate_Au, 
                        SDPCallbackFunction evaluate_ATu, 
                        void * data){

/*
    // copy block sizes
    primal_block_dim_.clear();
    for (int block = 0; block < primal_block_dim.size(); block++) {
        primal_block_dim_.push_back(primal_block_dim[block]);
    }
    // copy block ranks (same as dimension for now)
    primal_block_rank_.clear();
    for (int block = 0; block < primal_block_dim.size(); block++) {
        primal_block_rank_.push_back(primal_block_dim[block]);
    }

    c_            = c;
    x_            = x;
*/

    data_         = data;
    evaluate_Au_  = evaluate_Au;
    evaluate_ATu_ = evaluate_ATu;

    // cg solver
    std::shared_ptr<CGSolver> cg (new CGSolver(n_dual_));
    cg->set_max_iter(options_.get_int("CG_MAXITER"));
    double cg_convergence = options_.get_int("CG_CONVERGENCE");
    cg->set_convergence(cg_convergence);

    // the iterations
    outfile->Printf("\n");
    outfile->Printf("    initial primal energy: %20.12lf\n",C_DDOT(n_primal_,c->pointer(),1,x->pointer(),1));
    outfile->Printf("\n");
    outfile->Printf("      oiter");
    outfile->Printf(" iiter");
    outfile->Printf("        E(p)");
    outfile->Printf("        E(d)");
    outfile->Printf("      E(gap)");
    outfile->Printf("      mu");
    outfile->Printf("     eps(p)");
    outfile->Printf("     eps(d)\n");

    double primal_dual_energy_gap = 0.0;

    int oiter_local = 0;

    do {

        double start = omp_get_wtime();

        // evaluate mu * (b - Ax) for CG
        evaluate_Au(Au_, x, data);
        Au_->subtract(b);
        Au_->scale(-mu_);

        // evaluate A(c-z) ( but don't overwrite c! )
        z_->scale(-1.0);
        z_->add(c);
        evaluate_Au(cg_rhs_, z_, data);

        // add tau*mu*(b-Ax) to A(c-z) and put result in cg_rhs_
        cg_rhs_->add(Au_);

        // set convergence for CG problem (step 1 in table 1 of PRL 106 083001)
        double cg_conv_i = options_.get_double("CG_CONVERGENCE");
        if (oiter_ == 0)
            cg_conv_i = 0.01;
        else
            cg_conv_i = (primal_error_ > dual_error_) ? 0.01 * dual_error_ : 0.01 * primal_error_;
        if (cg_conv_i < cg_convergence)
            cg_conv_i = cg_convergence;
        cg->set_convergence(cg_conv_i);

        // solve CG problem (step 1 in table 1 of PRL 106 083001)
        cg->solve(Au_,y_,cg_rhs_,evaluate_cg_AATu,(void*)this);
        int iiter = cg->total_iterations();

        double end = omp_get_wtime();

        iiter_time_  += end - start;
        iiter_total_ += iiter;

        start = omp_get_wtime();

        // update primal and dual solutions
        Update_xz(x, c, primal_block_dim, evaluate_ATu, data);

        end = omp_get_wtime();

        oiter_time_ += end - start;

        // update mu (step 3)

        // evaluate || A^T y - c + z||
        evaluate_ATu(ATu_, y_, data);
        ATu_->add(z_);
        ATu_->subtract(c);
        dual_error_ = ATu_->norm();

        // evaluate || Ax - b ||
        evaluate_Au(Au_, x, data);
        Au_->subtract(b);
        primal_error_ = Au_->norm();

        // compute current primal and dual energies
        double energy_primal = C_DDOT(n_primal_,c->pointer(),1,x->pointer(),1);
        double energy_dual   = C_DDOT(n_dual_,b->pointer(),1,y_->pointer(),1);

        primal_dual_energy_gap = fabs(energy_primal-energy_dual);

        outfile->Printf("      %5i %5i %11.6lf %11.6lf %11.6le %7.3lf %10.5le %10.5le\n",
                    oiter_,iiter,energy_primal,energy_dual,primal_dual_energy_gap,mu_,primal_error_,dual_error_);

        oiter_++;
        oiter_local++;

        // don't update mu every iteration
        if ( oiter_ % options_.get_int("MU_UPDATE_FREQUENCY") == 0 && oiter_ > 0 ){
            mu_ = mu_ * primal_error_ / dual_error_;
        }

        if ( primal_error_ > r_convergence_ || dual_error_ > r_convergence_  || primal_dual_energy_gap > e_convergence_ ) {
            is_converged_ = false;
        }else {
            is_converged_ = true;
        }

        if ( oiter_local == maxiter ) break;
    }while( !is_converged_ );

}

std::shared_ptr<Vector> BPSDPSolver::ATAx_minus_ATb(std::shared_ptr<Vector> x,
                                                    std::shared_ptr<Vector> b,
                                                    SDPCallbackFunction evaluate_Au,
                                                    SDPCallbackFunction evaluate_ATu,
                                                    void * data){

    std::shared_ptr<Vector> ret (new Vector(n_primal_));

    evaluate_Au(Au_, x, data);
    Au_->subtract(b);

    evaluate_ATu(ret,Au_,data);

    return ret;

}

std::shared_ptr<Vector> BPSDPSolver::ATy_plus_z_minus_c(std::shared_ptr<Vector> c,
                                                        SDPCallbackFunction evaluate_ATu,
                                                        void * data){

    std::shared_ptr<Vector> ret (new Vector(n_primal_));

    evaluate_ATu(ret,y_,data);
    ret->add(z_);
    ret->subtract(c);

    return ret;
}

// update x and z
void BPSDPSolver::Update_xz(std::shared_ptr<Vector> x, std::shared_ptr<Vector> c, std::vector<int> primal_block_dim, SDPCallbackFunction evaluate_ATu, void * data) {

    // evaluate M(mu*x + ATy - c)
    evaluate_ATu(ATu_, y_, data);
    ATu_->subtract(c);
    x->scale(mu_);
    ATu_->add(x);

    // loop over each block of x/z
    for (int i = 0; i < primal_block_dim.size(); i++) {
        if ( primal_block_dim[i] == 0 ) continue;
        int myoffset = 0;
        for (int j = 0; j < i; j++) {
            myoffset += primal_block_dim[j] * primal_block_dim[j];
        }

        SharedMatrix mat     (new Matrix(primal_block_dim[i],primal_block_dim[i]));
        SharedMatrix eigvec  (new Matrix(primal_block_dim[i],primal_block_dim[i]));
        SharedMatrix eigvec2 (new Matrix(primal_block_dim[i],primal_block_dim[i]));
        SharedVector eigval  (new Vector(primal_block_dim[i]));

        double ** mat_p = mat->pointer();
        double * A_p    = ATu_->pointer();

        for (int p = 0; p < primal_block_dim[i]; p++) {
            for (int q = p; q < primal_block_dim[i]; q++) {
                double dum = 0.5 * ( A_p[myoffset + p * primal_block_dim[i] + q] +
                                     A_p[myoffset + q * primal_block_dim[i] + p] );
                mat_p[p][q] = mat_p[q][p] = dum;
            }
        }

        mat->diagonalize(eigvec,eigval);

        // separate U+ and U-, transform back to nondiagonal basis

        double * eval_p   = eigval->pointer();
        double ** evec_p  = eigvec->pointer();
        double ** evec2_p = eigvec2->pointer();

        double * x_p      = x->pointer();
        double * z_p      = z_->pointer();

        // (+) part
        long int mydim = 0;
        for (long int j = 0; j < primal_block_dim[i]; j++) {
            if ( eval_p[j] > 0.0 ) {
                for (long int q = 0; q < primal_block_dim[i]; q++) {
                    mat_p[q][mydim]   = evec_p[q][j] * eval_p[j]/mu_;
                    evec2_p[q][mydim] = evec_p[q][j];
                }
                mydim++;
            }
        }
        F_DGEMM('t','n',primal_block_dim[i],primal_block_dim[i],mydim,1.0,&mat_p[0][0],primal_block_dim[i],&evec2_p[0][0],primal_block_dim[i],0.0,x_p+myoffset,primal_block_dim[i]);

        // (-) part
        mydim = 0;
        for (long int j = 0; j < primal_block_dim[i]; j++) {
            if ( eval_p[j] < 0.0 ) {
                for (long int q = 0; q < primal_block_dim[i]; q++) {
                    mat_p[q][mydim]   = -evec_p[q][j] * eval_p[j];
                    evec2_p[q][mydim] =  evec_p[q][j];
                }
                mydim++;
            }
        }
        F_DGEMM('t','n',primal_block_dim[i],primal_block_dim[i],mydim,1.0,&mat_p[0][0],primal_block_dim[i],&evec2_p[0][0],primal_block_dim[i],0.0,z_p+myoffset,primal_block_dim[i]);

    }

/*
printf("||z|| before %20.12lf\n",z_->norm());
    // lbfgs container for dual solution vector, z
    lbfgsfloatval_t * lbfgs_vars_z  = lbfgs_malloc(n_primal_);
    srand(0);
    for (int i = 0; i < n_primal_; i++) {
        lbfgs_vars_z[i] = 2.0 * ( (double)rand()/RAND_MAX - 0.5 ) / 1000.0;
    }
    build_z(lbfgs_vars_z);

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);

    // adjust default lbfgs parameters
    param.max_iterations = options_.get_int("MAXITER");
    //param.epsilon        = 1e-12;//options_.get_double("LBFGS_CONVERGENCE");

    // initial objective function value default parameters
    lbfgsfloatval_t lag_z = evaluate_gradient_z(lbfgs_vars_z,ATu_->pointer());
    lbfgs(n_primal_,lbfgs_vars_z,&lag_z,lbfgs_evaluate_z,monitor_lbfgs_progress,(void*)this,&param);

    // build z = r.rT
    build_z(lbfgs_vars_z);

    printf("hey %5i lbfgs iterations\n",lbfgs_iter_);fflush(stdout);

    free(lbfgs_vars_z);
*/

}

/*
double BPSDPSolver::evaluate_gradient_z(const lbfgsfloatval_t * r, lbfgsfloatval_t * g) {

    // L = || c - ATy + z ||^2 + || x.z ||^2, fixed y and z

    // pointers
    double * r_p   = (double*)r;
    double * x_p   = x_->pointer();
    double * z_p   = z_->pointer();
    double * ATu_p = ATu_->pointer();

    // build z = r.rT
    build_z(r_p);

    // evaluate || A^T y - c + z||
    evaluate_ATu_(ATu_, y_, data_);
    ATu_->subtract(c_);
    ATu_->add(z_);
    double dual_error = ATu_->norm();

    // z.x = 0
    double overlap = x_->vector_dot(z_);

    // evaluate lagrangian
    double lagrangian = dual_error*dual_error + overlap*overlap;

    int off_nn = 0;
    for (int block = 0; block < primal_block_dim_.size(); block++) {
        int n = primal_block_dim_[block];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                ATu_p[i*n+j + off_nn]  = ( ATu_p[i*n+j + off_nn] + ATu_p[j*n+i + off_nn] );
                ATu_p[i*n+j + off_nn] += (  x_p[i*n+j + off_nn] +  x_p[j*n+i + off_nn] ) * overlap;
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
        F_DGEMM('n', 'n', n, m, n, 2.0, ATu_p + off_nn, n, r_p + off_nm, n, 0.0, g + off_nm, n);
        off_nn += n*n;
        off_nm += n*m;
    }

    return lagrangian;
}

// build z = r.rT
void BPSDPSolver::build_z(double * r){

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
*/

}

