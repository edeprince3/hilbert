#include <psi4/libplugin/plugin.h>
#include <psi4/psi4-dec.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsi4util/PsiOutStream.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/molecule.h>
#include <psi4/libmints/basisset.h>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/psifiles.h>
#include <psi4/libqt/qt.h>

#include <misc/blas.h>
#include <misc/omp.h>

#include "orbital_optimizer.h"

using namespace psi;
using namespace fnocc;

namespace hilbert{

void OrbitalOptimizer::optimize_orbitals(double * d2, double * d1, double * tei, double * oei, double * transformation_matrix,\
                                         double& dE_orbopt, double& gnorm_orbopt, bool& converged_orbopt, int& iter_orbopt){

    Initiate_Optimizer(d2, d1, oei, transformation_matrix);
 
    // initialize iteration variables

    bool eval_G   = true; bool converged_e = false; bool converged_g = false;
    double r_current = 1.0e0; double r_new=-1.0e0;
    double final_E = 0.0e0; double initial_E = 0.0e0;
    double E_new; double E_init;
    double dE_quad; double dE; double dE_ratio;

    outfile->Printf("        iter       Energy        dE       |g|     r\n");

    iter_orbopt = 0;

    for ( int iter = 0; iter < maxiter_; iter++){

        if ( eval_G ) {

            E_init = Compute_EGH( d2, d1, tei, oei );

            if ( iter == 0 ) initial_E = E_init;

            dE_quad = Compute_Quadratic_dE();

            Precondition_Step();

            C_DSCAL(Nrot_,r_current,&kappa_current_[0],1);
            C_DCOPY(Nrot_,&kappa_current_[0],1,&kappa_save_[0],1);

        }

        ExponentiateKappa( kappa_current_ );

        TransformIntegrals( tei, oei, transformation_matrix);

        E_new    = ComputeE(d2, d1, tei, oei);
        dE       = E_new - E_init;

        outfile->Printf("        %4i %12.6f %9.2e %9.2e %5.3f\n",iter,E_new,dE,gradient_norm_,r_current);

        if ( gradient_norm_ < g_convergence_ ) converged_g = true;
        if ( fabs(dE) < e_convergence_ )       converged_e = true;

        iter_orbopt++;

        if ( converged_e || converged_g || iter == maxiter_ - 1 ) break;

        eval_G = Update_r(dE_quad, dE, r_current);

        if ( !eval_G ) Update_kappa(r_current);

    }

    if ( dE > 0.0e0 ) Restore_last_step(tei, oei, transformation_matrix, converged_e, converged_g, eval_G);

    final_E = ComputeE(d2, d1, tei,oei);

    converged_orbopt = false;
    if ( converged_e && converged_g ) converged_orbopt = true;
    dE_orbopt = final_E - initial_E;
    gnorm_orbopt =gradient_norm_;

    outfile->Printf("\n");

    Finalize_Optimizer(d2, d1, oei, transformation_matrix);

}

void OrbitalOptimizer::Update_kappa(double r_new){

    // at this point, it is assumed that
    // kappa_ref_     = -g/H               (preconditioned gradients)
    // kappa_current_ = r_old * kappa_ref_ (rejected step)
    // kappa_new_     = r_new * kappa_ref_ (next step)
    // and thus
    // U = e^(kappa_new_)*e^(-kappa_current_)

    C_DCOPY(Nrot_,&kappa_ref_[0],1,&kappa_new_[0],1);
    C_DSCAL(Nrot_,r_new,&kappa_new_[0],1);
    C_DAXPY(Nrot_,-1.0e0,&kappa_current_[0],1,&kappa_new_[0],1);

    C_DCOPY(Nrot_,&kappa_new_[0],1,&kappa_current_[0],1);

    C_DAXPY(Nrot_,1.0e0,&kappa_current_[0],1,&kappa_save_[0],1);

}

bool OrbitalOptimizer::Update_r(double dE_Quad, double dE, double& r){

    double r_inc_tol = 0.75e0;
    double r_dec_tol = 0.25e0;
    double r_inc_fac = 1.20e0;
    double r_dec_fac = 0.70e0;

    bool eval_G = true;

    if ( dE > 0.0e0 ) eval_G = false;

    double dE_ratio = dE/dE_Quad;

    if ( dE < 0.0e0 ){

        if ( dE_ratio > r_inc_tol ) r *= r_inc_fac;

        if ( dE_ratio < r_dec_tol ) r *= r_dec_fac;

    }

    else {

        r *= r_dec_fac;

    }

    return eval_G;

}

double OrbitalOptimizer::Compute_EGH(double * d2, double * d1, double * tei, double * oei){

    ComputeG( d2, d1, tei, oei );

    double E = ComputeE_2index( d1, oei );

//    ComputeG( d2, d1, tei, oei );

    ComputeH_diag( d2, d1, tei );

    return E;

}

double OrbitalOptimizer::Compute_Quadratic_dE(){

    double E = 2.0e0 * C_DDOT(Nrot_,orbital_gradient_,1,kappa_current_,1);

    for ( int pq =0; pq < Nrot_; pq++){

        E += kappa_current_[pq] * diagonal_orbital_hessian_[pq] * kappa_current_[pq];

    }

    E *= 0.5e0;

    return E;

}

void OrbitalOptimizer::Initiate_Optimizer(double * d2, double * d1, double * oei, double * transformation_matrix){

    // allocate mapping arrays
    Alloc_OindMap();

    // compute mapping arrays
    Determine_OindMap();

    // sort oei, d1, and d2 (needs maping arrays)
    int direction = 1;
    Sort( d2, d1, oei, transformation_matrix, direction );

    // allocate indexing arrays
    Alloc_Gradient();
    Alloc_Hessian();
    Determine_GradIndMap();

    // allocate transformation arrays
    Alloc_U();

    // allocate scratch arrays for optimizer
    Alloc_Optimizer();

}

void OrbitalOptimizer::Finalize_Optimizer(double * d2, double * d1, double * oei, double * transformation_matrix){

    // sort d2, d1, oei (needs mapping arrays)
    int direction = -1;
    Sort( d2, d1, oei, transformation_matrix, direction );

    // allocate mapping arrays
    Dealloc_OindMap();

    // allocate indexing arrays
    Dealloc_Gradient();
    Dealloc_Hessian();

    // allocate transformation arrays
    Dealloc_U();

    // allocate scratch arrays for optimizer
    Dealloc_Optimizer();

}

void OrbitalOptimizer::Restore_last_step( double * tei, double * oei, double * MOcoeff, bool converged_e, bool converged_g, bool eval_G){


    if ( converged_e || converged_g ) {

        C_DSCAL(Nrot_,-1.0e0,&kappa_current_[0],1);
        ExponentiateKappa( kappa_current_ );
        TransformIntegrals( tei, oei, MOcoeff);

    }

    if ( !eval_G ) {

        printf("heyho\n");

        C_DSCAL(Nrot_,-1.0e0,&kappa_save_[0],1);
        ExponentiateKappa( kappa_save_ );
        TransformIntegrals( tei, oei, MOcoeff );

    }

}

void OrbitalOptimizer::Precondition_Step(){

    for ( int pq = 0; pq < Nrot_; pq++){

        double Hval = diagonal_orbital_hessian_[pq];

        if ( Hval < 0.0e0  ) Hval = - Hval;
        if ( Hval < 1.0e-2 ) Hval = 1.0e-2;

        kappa_ref_[pq] = - orbital_gradient_[pq] / Hval;

    }

    C_DCOPY(Nrot_,&kappa_ref_[0],1,&kappa_current_[0],1);

}

void OrbitalOptimizer::Alloc_Optimizer(){

    kappa_current_ = (double *)malloc(Nrot_*sizeof(double));
    memset((void*)kappa_current_,'\0',Nrot_*sizeof(double));

    kappa_save_ = (double *)malloc(Nrot_*sizeof(double));
    memset((void*)kappa_save_,'\0',Nrot_*sizeof(double));

    kappa_ref_ = (double *)malloc(Nrot_*sizeof(double));
    memset((void*)kappa_ref_,'\0',Nrot_*sizeof(double));

    kappa_new_ = (double *)malloc(Nrot_*sizeof(double));
    memset((void*)kappa_new_,'\0',Nrot_*sizeof(double));

}

void OrbitalOptimizer::Dealloc_Optimizer(){

    free(kappa_current_);
    free(kappa_save_);
    free(kappa_ref_);
    free(kappa_new_);

}


}
