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

typedef void (*BPSDPCallbackFunction)(std::shared_ptr<Vector>,std::shared_ptr<Vector>,void *);

BPSDPSolver::BPSDPSolver(long int n_primal, long int n_dual, Options & options)
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

    cg_rhs_ = (std::shared_ptr<Vector>)(new Vector(n_dual_));

    e_convergence_ = options_.get_double("E_CONVERGENCE");
    r_convergence_ = options_.get_double("R_CONVERGENCE");

}

BPSDPSolver::~BPSDPSolver(){
}

void BPSDPSolver::solve(std::shared_ptr<Vector> x,   
                        std::shared_ptr<Vector> b, 
                        std::shared_ptr<Vector> c,
                        std::vector<int> primal_block_dim,
                        int maxiter,
                        BPSDPCallbackFunction evaluate_Au, 
                        BPSDPCallbackFunction evaluate_ATu, 
                        CGCallbackFunction evaluate_CG_LHS, 
                        void * data){

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
    outfile->Printf("      E gap)");
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
        cg->solve(Au_,y_,cg_rhs_,evaluate_CG_LHS,data);
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

        if ( oiter_local == maxiter ) break;

    }while( primal_error_ > r_convergence_ || dual_error_ > r_convergence_  || primal_dual_energy_gap > e_convergence_ );

}

// update x and z
void BPSDPSolver::Update_xz(std::shared_ptr<Vector> x, std::shared_ptr<Vector> c, std::vector<int> primal_block_dim, BPSDPCallbackFunction evaluate_ATu, void * data) {

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
}

}

