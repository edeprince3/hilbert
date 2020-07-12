/*
 *@BEGIN LICENSE
 *
 * pp2RDM, a plugin to:
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <psi4/libpsi4util/PsiOutStream.h>
#include <psi4/libtrans/integraltransform.h>

#include "../../external/install/include/lbfgs.h"

/*
//#include <lbfgs.h>
#include "/opt/local/include/lbfgs.h"
//#include "/usr/local/include/lbfgs.h"
*/

#include "pp2rdm_solver.h"

using namespace psi;
using namespace pp2rdm;

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

static lbfgsfloatval_t evaluate(void * instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step) {

    pp2RDMSolver * pp2RDM = reinterpret_cast<pp2RDMSolver*>(instance);

    // energy
    double energy = pp2RDM->evaluate_variational_energy();

    // gradient
    pp2RDM->evaluate_gradient(g);

    // it is possible that L-BFGS searches outside the domain of feasible amplitudes,
    // so, if the energy is not a number, just tell L-BFGS that it is unacceptably high

    if ( isnan(energy) ) energy = 9e9;

    // TODO: should probably do something to gradient in the event that isnan(energy).
    for (int i = 0; i < n; i++) {
        if ( isnan(g[i]) ) g[i] = 9e9;
    }

    return energy;
}

static int progress(
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
    pp2RDMSolver * sdp = reinterpret_cast<pp2RDMSolver*>(instance);
    sdp->set_number_of_lbfgs_iterations(k);
    return 0;
}

namespace psi { namespace pp2rdm {

void pp2RDMSolver::set_number_of_lbfgs_iterations(int iter) {

    ci_iter_ = iter;

}

double pp2RDMSolver::pp2rdm_lbfgs_iterations(int & ci_iter) {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    // build integrals and fock matrix
    setup_integrals();

    // check gradient for accuracy
    if ( options_.get_bool("CHECK_GRADIENT") ) {

        check_gradient();
        exit(0);

    }

    if ( options_.get_bool("CHECK_HESSIAN") ) {

        check_hessian();
        exit(0);

    }

    // allocate member for variables 
    lbfgsfloatval_t *  vars  = lbfgs_malloc(o*v);
    lbfgsfloatval_t * dvars  = lbfgs_malloc(o*v);

    memset((void*)vars, '\0',o*v*sizeof(double));

    // have t2_ point to vars.  reallocate t2_ at the end of this function.
    free(t2_);
    t2_ = vars;

    // restart from old solution?
    if ( restart_from_checkpoint_file_ ) {
        std::shared_ptr<PSIO> psio ( new PSIO() );
        psio->open(PSIF_V2RDM_CHECKPOINT,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_V2RDM_CHECKPOINT,"AMPLITUDES",(char*)t2_,o * v * sizeof(double));
        psio->close(PSIF_V2RDM_CHECKPOINT,1);
    }

    lbfgsfloatval_t f = evaluate_variational_energy();
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);

    // adjust default parameters
    //param.max_iterations = 100;
    //param.max_linesearch = 100;
    //param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
    //param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;
    //param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
    //param.gtol = 1e-3;  // default 0.9.  should be > ftol (1e-4) and < 1.0
    //param.epsilon = 1e-12; // default 1e-5
    //param.m = 8;

    // minimize energy
    //double * save = (double*)malloc(o*v*sizeof(double));
    //C_DCOPY(o*v,vars,1,save,1);
    int return_value = lbfgs(o*v,vars,&f,evaluate,progress,(void*)this,&param);
    lbfgs_error_check(return_value);
    //if ( return_value != 0 ) {
    //    C_DCOPY(o*v,save,1,vars,1);
    //}

    double en = f;

    // reallocate memory for amplitudes
    t2_ = (double*)malloc(o*v*sizeof(double));
    C_DCOPY(o*v,vars,1,t2_,1);

    free(vars);

    // save solution
    std::shared_ptr<PSIO> psio ( new PSIO() );
    psio->open(PSIF_V2RDM_CHECKPOINT,PSIO_OPEN_NEW);
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"AMPLITUDES",(char*)t2_,o * v * sizeof(double));
    psio->close(PSIF_V2RDM_CHECKPOINT,1);
    restart_from_checkpoint_file_ = true;

    return en;
}

// evaluate gradient of the energy
void pp2RDMSolver::evaluate_gradient(double * gradient) {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    memset((void*)gradient,'\0',o * v * sizeof(double));

    evaluate_residual(gradient);

    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            gradient[i * v + a ] += 2.0 * t2_[i * v + a] * ( fv_[a] - fo_[i] );
        }
    }

    C_DSCAL(o*v,2.0,gradient,1);

    if ( options_.get_str("P2RDM_TYPE") == "ACPF" || options_.get_str("P2RDM_TYPE") == "AQCC" ) {

        double fac = 0.0;
        if ( options_.get_str("P2RDM_TYPE") == "ACPF" ) {
            fac = 1.0 / o;
        }else if ( options_.get_str("P2RDM_TYPE") == "AQCC" ) {
            fac = 1.0 - (2.0 * o - 2.0) * (2.0 * o - 3.0) / (2.0 * o * (2.0 * o - 1.0));
        }

        double nrm = 1.0 + fac * C_DDOT(o * v, t2_, 1, t2_, 1);

        C_DSCAL(o*v,1.0/nrm,gradient,1);

        double en = evaluate_variational_energy();

        C_DAXPY(o*v,-2.0 * en / nrm * fac,t2_,1,gradient,1);

    }
}

double pp2RDMSolver::evaluate_variational_energy() {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    //memset((void*)r2_,'\0',o * v * sizeof(double));
    double * sigma = (double*)malloc(o*v*sizeof(double));
    memset((void*)sigma,'\0',o * v * sizeof(double));

    // construct c0
    Normalization();

    // evaluate <D|H-E0|0 + D>
    evaluate_sigma(sigma);
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            sigma[i * v + a ] += 2.0 * t2_[i * v + a] * ( fv_[a] - fo_[i] );
        }
    }

    // evaluate <D|H|0>
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            sigma[i * v + a ] += v_iaia_[i * v + a] * t0_[i * v + a];
        }
    }

    double en = C_DDOT(o*v,t2_,1,sigma,1);

    if ( options_.get_str("P2RDM_TYPE") == "ACPF" || options_.get_str("P2RDM_TYPE") == "AQCC" ) {

        double fac = 0.0;
        if ( options_.get_str("P2RDM_TYPE") == "ACPF" ) {
            fac = 1.0 / o;
        }else if ( options_.get_str("P2RDM_TYPE") == "AQCC" ) {
            fac = 1.0 - (2.0 * o - 2.0) * (2.0 * o - 3.0) / (2.0 * o * (2.0 * o - 1.0));
        }

        double nrm = 1.0 + fac * C_DDOT(o * v, t2_, 1, t2_, 1);
        en /= nrm;

    }
    free(sigma);

    return en;
}

void pp2RDMSolver::evaluate_a22(std::shared_ptr<Matrix> a22) {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    a22->zero();

    double **a22_p = a22->pointer();

    // diagonals (excluding c0 parts)
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            int ia = i * v + a;

            a22_p[ia][ia]  = 4.0 * (fv_[a] - fo_[i]);

            // -2 * (2 (ii|aa) - (ia|ai) c(ai)^2
            a22_p[ia][ia] -= 4.0 * ( 2.0 * v_iiaa_[i * v + a] - v_iaia_[i * v + a] );

            // c(ia) * c(ib) * (ab|ab)
            for (int b = 0; b < v; b++) {
                int ib = i * v + b;
                a22_p[ia][ib] += v_abab_[a * v + b];
                a22_p[ib][ia] += v_abab_[a * v + b];
            }

            // c(ia) * c(ja) * (ji|ji)
            for (int j = 0; j < o; j++) {
                int ja = j * v + a;
                a22_p[ia][ja] += v_ijij_[i * o + j];
                a22_p[ja][ia] += v_ijij_[i * o + j];
            }

            // all that remains are c0-dependent parts
            if ( options_.get_str("P2RDM_TYPE") == "K" ) {

                // diagonal contribution

                double dum = 0.0;
                //for (int k = 0; k < o; k++) {
                //    for (int c = 0; c < v; c++) {
                //        int kc = k * v + c;
                //        dum -= 2.0 * t2_[kc] / t0_[kc] * v_iaia_[kc] * ( (i==k) + (a==c) - (i==k)*(a==c) ); 
                //    }
                //}
                for (int c = 0; c < v; c++) {
                    int ic = i * v + c;
                    dum -= 2.0 * t2_[ic] / t0_[ic] * v_iaia_[ic];
                }
                for (int k = 0; k < o; k++) {
                    int ka = k * v + a;
                    dum -= 2.0 * t2_[ka] / t0_[ka] * v_iaia_[ka];
                }
                dum += 2.0 * t2_[ia] / t0_[ia] * v_iaia_[ia];

                a22_p[ia][ia] += dum;

                // off-diagonal contributions

                for (int j = 0; j < o; j++) {
                    for (int b = 0; b < v; b++) {
                        int jb = j * v + b;

                        double dum = 0.0;
                        //for (int k = 0; k < o; k++) {
                        //    for (int c = 0; c < v; c++) {
                        //        int kc = k * v + c;
                        //        dum -= 2.0 * t2_[kc] / t0_[kc] / t0_[kc] / t0_[kc] * v_iaia_[kc]
                        //             * ( (i==k) + (a==c) - (i==k)*(a==c) ) 
                        //             * ( (j==k) + (b==c) - (j==k)*(b==c) ) ;
                        //    }
                        //}
                        //a22_p[ia][jb] += t2_[ia] * t2_[jb] * dum;
                        //a22_p[ia][jb] -= 2.0 * t2_[ia] / t0_[jb] * v_iaia_[jb] * ( (i==j) + (a==b) - (i==j)*(a==b) );
                        //a22_p[ia][jb] -= 2.0 * t2_[jb] / t0_[ia] * v_iaia_[ia] * ( (i==j) + (a==b) - (i==j)*(a==b) );

                        // (i==k) * (b==c)
                        int ib = i * v + b;
                        dum -= 2.0 * t2_[ib] / t0_[ib] / t0_[ib] / t0_[ib] * v_iaia_[ib];

                        // (a==c) * (j==k) 
                        int ja = j * v + a;
                        dum -= 2.0 * t2_[ja] / t0_[ja] / t0_[ja] / t0_[ja] * v_iaia_[ja];

                        a22_p[ia][jb] += t2_[ia] * t2_[jb] * dum;


                    }
                }

                // unrolled

                // (i==j)
                for (int b = 0; b < v; b++) {

                    double dum = 0.0;

                    // (i==k) * (j==k)
                    for (int c = 0; c < v; c++) {
                        int ic = i * v + c;
                        dum -= 2.0 * t2_[ic] / t0_[ic] / t0_[ic] / t0_[ic] * v_iaia_[ic];
                    }

                    // (i==k) * (j==k) * (b==c)
                    int ib = i * v + b;
                    dum += 2.0 * t2_[ib] / t0_[ib] / t0_[ib] / t0_[ib] * v_iaia_[ib];

                    // (i==k) * (a==c) * (j==k)
                    dum += 2.0 * t2_[ia] / t0_[ia] / t0_[ia] / t0_[ia] * v_iaia_[ia];

                    a22_p[ia][ib] += t2_[ia] * t2_[ib] * dum;

                    // (i==j)
                    a22_p[ia][ib] -= 2.0 * t2_[ia] / t0_[ib] * v_iaia_[ib];
                    a22_p[ia][ib] -= 2.0 * t2_[ib] / t0_[ia] * v_iaia_[ia];

                }

                // (a==b)
                for (int j = 0; j < o; j++) {

                    double dum = 0.0;

                    // (a==c) * (b==c) 
                    for (int k = 0; k < o; k++) {
                        int ka = k * v + a;
                        dum -= 2.0 * t2_[ka] / t0_[ka] / t0_[ka] / t0_[ka] * v_iaia_[ka];
                    }

                    // (a==c) * (j==k) (b==c) 
                    int ja = j * v + a;
                    dum += 2.0 * t2_[ja] / t0_[ja] / t0_[ja] / t0_[ja] * v_iaia_[ja];

                    // (i==k) * (a==c) * (b==c)
                    dum += 2.0 * t2_[ia] / t0_[ia] / t0_[ia] / t0_[ia] * v_iaia_[ia];

                    a22_p[ia][ja] += t2_[ia] * t2_[ja] * dum;

                    // (a==b)
                    a22_p[ia][ja] -= 2.0 * t2_[ia] / t0_[ja] * v_iaia_[ja];
                    a22_p[ia][ja] -= 2.0 * t2_[ja] / t0_[ia] * v_iaia_[ia];

                }

                // more diagonal terms (i==j) (a==b)

                a22_p[ia][ia] += 2.0 * t2_[ia] / t0_[ia] * v_iaia_[ia];
                a22_p[ia][ia] += 2.0 * t2_[ia] / t0_[ia] * v_iaia_[ia];

                a22_p[ia][ia] -= 2.0 * t2_[ia] * t2_[ia] * t2_[ia] / t0_[ia] / t0_[ia] / t0_[ia] * v_iaia_[ia];

            }else if ( options_.get_str("P2RDM_TYPE") == "CEPA(1)" ) {
                // diagonal contribution
                double dum = 0.0;
                for (int k = 0; k < o; k++) {
                    for (int c = 0; c < v; c++) {
                        int kc = k * v + c;
                        dum -= 2.0 * t2_[kc] / t0_[kc] * v_iaia_[kc] * ( (i==k) );
                    }
                }
                a22_p[ia][ia] += dum;

                // off-diagonal contributions
                for (int j = 0; j < o; j++) {
                    for (int b = 0; b < v; b++) {
                        int jb = j * v + b;

                        double dum = 0.0;
                        for (int k = 0; k < o; k++) {
                            for (int c = 0; c < v; c++) {
                                int kc = k * v + c;
                                dum -= 2.0 * t2_[kc] / t0_[kc] / t0_[kc] / t0_[kc] * v_iaia_[kc]
                                     * ( (i==k) )
                                     * ( (j==k) );
                            }
                        }
                        a22_p[ia][jb] += t2_[ia] * t2_[jb] * dum;

                        a22_p[ia][jb] -= 2.0 * t2_[ia] / t0_[jb] * v_iaia_[jb] * ( (i==j) );
                        a22_p[ia][jb] -= 2.0 * t2_[jb] / t0_[ia] * v_iaia_[ia] * ( (i==j) );

                    }
                }

            }else if ( options_.get_str("P2RDM_TYPE") == "CID" ) {

                // diagonal contribution
                double dum = 0.0;
                for (int kc = 0; kc < o*v; kc++) {
                    dum -= 2.0 * t2_[kc] / t0_[kc] * v_iaia_[kc];
                }
                a22_p[ia][ia] += dum;

                // off-diagonal contributions
                for (int jb = 0; jb < o*v; jb++) {

                    double dum = 0.0;
                    for (int kc = 0; kc < o*v; kc++) {
                        dum -= 2.0 * t2_[kc] / t0_[kc] / t0_[kc] / t0_[kc] * v_iaia_[kc];
                    }
                    a22_p[ia][jb] += t2_[ia] * t2_[jb] * dum;

                    a22_p[ia][jb] -= 2.0 * t2_[ia] / t0_[jb] * v_iaia_[jb];
                    a22_p[ia][jb] -= 2.0 * t2_[jb] / t0_[ia] * v_iaia_[ia];

                }

            }else if ( options_.get_str("P2RDM_TYPE") == "CEPA(0)" ) {

                // nothing to add

            }else if ( options_.get_str("P2RDM_TYPE") == "ACPF" ) {

                throw PsiException("analytic hessians not yet implemented for p2rdm_type acpf",__FILE__,__LINE__);

            }else if ( options_.get_str("P2RDM_TYPE") == "AQCC" ) {

                throw PsiException("analytic hessians not yet implemented for p2rdm_type aqcc",__FILE__,__LINE__);

            }


        }
    }

}

void pp2RDMSolver::evaluate_numerical_hessian(std::shared_ptr<Matrix> hessian) {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    double * tmp_p1 = (double*)malloc(o*v*sizeof(double));
    double * tmp_p2 = (double*)malloc(o*v*sizeof(double));
    double * tmp_m1 = (double*)malloc(o*v*sizeof(double));
    double * tmp_m2 = (double*)malloc(o*v*sizeof(double));

    double **hp = hessian->pointer();

    double h = 1e-4;
    for (int i = 0; i < o*v; i++) {

        double save = t2_[i];

        t2_[i] = save + h;
        Normalization();
        evaluate_gradient(tmp_p1);

        t2_[i] = save + 2.0 * h;
        Normalization();
        evaluate_gradient(tmp_p2);

        t2_[i] = save - h;
        Normalization();
        evaluate_gradient(tmp_m1);

        t2_[i] = save - 2.0 * h;
        Normalization();
        evaluate_gradient(tmp_m2);

        
        for (int j = 0; j < o*v; j++) {
            hp[i][j] = (-tmp_p2[j] + 8.0 * tmp_p1[j] - 8.0 * tmp_m1[j] + tmp_m2[j]) / (12.0 * h);
        }

        t2_[i] = save;

    }

    free(tmp_p1);
    free(tmp_p2);
    free(tmp_m1);
    free(tmp_m2);

}

void pp2RDMSolver::check_gradient() {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    // random starting amplitudes on interval [-0.01:0.01]
    srand(time(0));
    for (int i = 0; i < o*v; i++) {
        t2_[i] = 2.0 * ( (double)rand()/RAND_MAX - 0.5 ) * 0.1;
    }
    Normalization();
    evaluate_gradient(r2_);
    double * temp = (double*)malloc(o*v*sizeof(double));
    C_DCOPY(o*v,r2_,1,temp,1);

    double h = 1e-4;
    for (int i = 0; i < o*v; i++) {

        double save = t2_[i];

        t2_[i] = save + h;
        double dump1 = evaluate_variational_energy();

        t2_[i] = save + 2.0 * h;
        double dump2 = evaluate_variational_energy();

        t2_[i] = save - h;
        double dumm1 = evaluate_variational_energy();

        t2_[i] = save - 2.0 * h;
        double dumm2 = evaluate_variational_energy();

        double myder = (-dump2 + 8.0 * dump1 - 8.0 * dumm1 + dumm2) / (12.0 * h);

        t2_[i] = save;

        printf("%5i %20.12le %20.12le %20.12le %20.12le %30.12lf %30.12lf %30.12lf\n",i,dump2,dump1,dumm1,dumm2,temp[i],myder,(temp[i]-myder)/temp[i] * 100.0);
    }
    free(temp);

}

void pp2RDMSolver::check_hessian() {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    // random starting amplitudes on interval [-0.01:0.01]
    srand(time(0));
    for (int i = 0; i < o*v; i++) {
        t2_[i] = 2.0 * ( (double)rand()/RAND_MAX - 0.5 ) * 0.05;
    }
    Normalization();

    std::shared_ptr<Matrix> hessian(new Matrix(o*v,o*v));
    evaluate_a22(hessian);
    double ** hp = hessian->pointer();

    std::shared_ptr<Matrix> approx_hessian(new Matrix(o*v,o*v));
    evaluate_numerical_hessian(approx_hessian);
    double ** ahp = approx_hessian->pointer();

    double h = 1e-4;
    double total = 0.0;
    for (int i = 0; i < o*v; i++) {
        for (int j = 0; j < o*v; j++) {
            printf("%5i %5i %20.12lf %20.12lf %20.12lf %20.12lf\n",i,j,hp[i][j],ahp[i][j],hp[i][j]-ahp[i][j],(hp[i][j]-ahp[i][j])/hp[i][j] * 100.0);

            total += (hp[i][j]-ahp[i][j]) * (hp[i][j]-ahp[i][j]);
        }

    }
    printf("|H - Happrox| = %20.12lf\n",sqrt(total));

}

double pp2RDMSolver::pp2rdm_newton_raphson_iterations(int & ci_iter) {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    // build integrals and fock matrix
    setup_integrals();

    // check gradient for accuracy
    if ( options_.get_bool("CHECK_GRADIENT") ) {

        check_gradient();
        exit(0);

    }

    if ( options_.get_bool("CHECK_HESSIAN") ) {

        check_hessian();
        exit(0);

    }

    // restart from old solution?
    if ( restart_from_checkpoint_file_ ) {
        std::shared_ptr<PSIO> psio ( new PSIO() );
        psio->open(PSIF_V2RDM_CHECKPOINT,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_V2RDM_CHECKPOINT,"AMPLITUDES",(char*)t2_,o * v * sizeof(double));
        psio->close(PSIF_V2RDM_CHECKPOINT,1);
    }

    std::shared_ptr<Matrix> hessian ( new Matrix(o*v,o*v) );
    double ** hp = hessian->pointer();

    ci_iter = 0;

    // start from pMP2 guess
    //evaluate_residual(r2_);
    //// update amplitudes
    //for (int i = 0; i < o; i++) {
    //    for (int a = 0; a < v; a++) {
    //        r2_[i * v + a] *= -0.5 / ( fv_[a] - fo_[i] );
    //    }
    //}
    //C_DCOPY(o * v, r2_, 1, t2_, 1);

    double energy = evaluate_variational_energy();
    evaluate_gradient(r2_);

    outfile->Printf("\n");
    outfile->Printf(
      "    ==> Begin p2RDM iterations <==\n");
    outfile->Printf("\n");
    outfile->Printf(
      "        iter          energy       d(Energy)          |g(T)|\n");

    double dE = 0.0;
    double rms_grad = 0.0;

    double * step = (double*)malloc(o*v*sizeof(double));
    do { 

        // d = -H^{-1}g
        evaluate_a22(hessian);

        hessian->invert();
        for (int i = 0; i < o*v; i++) {
            step[i] = 0.0;
            for (int j = 0; j < o*v; j++) {
                step[i] += hp[i][j] * r2_[j];
            }
        }
        C_DAXPY(o*v,-1.0,step,1,t2_,1);
        double new_energy = evaluate_variational_energy();
        double scale = 1.0;
        while ( isnan(new_energy) ) {
            C_DAXPY(o*v,scale,step,1,t2_,1);
            scale *= 0.5;
            outfile->Printf("\n");
            outfile->Printf("        ==> Warning! scaling NR step by %10.6lf <==\n",scale);
            outfile->Printf("\n");
            C_DAXPY(o*v,-scale,step,1,t2_,1);
            new_energy = evaluate_variational_energy();
            if ( scale < 1e-9 ) {
                throw PsiException("NR step size is unreasonably small.",__FILE__,__LINE__);
            }
        }

        evaluate_gradient(r2_);

        rms_grad = C_DNRM2(o*v,r2_,1);
        dE = fabs(new_energy - energy);
        energy = new_energy;

        outfile->Printf("       %5i %15.10lf %15.10lf %15.10lf\n",ci_iter,energy,dE,rms_grad);

        ci_iter++;
    }while(dE > e_convergence_ || rms_grad > r_convergence_);

    free(step);

    return energy;
}

}}// end of namespaces
