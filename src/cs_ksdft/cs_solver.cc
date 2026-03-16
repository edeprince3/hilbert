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

#include "cs_solver.h"

#include <misc/blas.h>
#include <misc/omp.h>

#include <lbfgs.h>
#include <random>

#include <psi4/libqt/qt.h>
#include <psi4/libpsio/psio.hpp>
#include "psi4/libpsi4util/process.h"
#include <psi4/liboptions/liboptions.h>
#include <psi4/psifiles.h>

// pybind
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace psi;
using namespace fnocc;

namespace hilbert{

static lbfgsfloatval_t lbfgs_evaluate(void * instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step) {

    cs_solver * cs = reinterpret_cast<cs_solver*>(instance);
    double f = cs->evaluate_augmented_lagrangian();
    C_DCOPY(n,cs->dvars_,1,g,1);

    return f;
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
    cs_solver * cs = reinterpret_cast<cs_solver*>(instance);
    cs->iiter_ = k;
    return 0;
}

cs_solver::cs_solver(Options & options):
    options_(options)
{
}

void cs_solver::solve(long int nalpha,
                      long int nbeta,
                      std::vector<double> T,
                      std::vector<double> super_phi,
                      std::vector<double> reference_rho,
                      std::vector<std::vector<double> > grid) {

    nalpha_ = nalpha;
    nbeta_ = nbeta;
    //eref_ = eref;
    T_.insert(T_.begin(), T.begin(), T.end());
    super_phi_.insert(super_phi_.begin(), super_phi.begin(), super_phi.end());
    temp_phi_.insert(temp_phi_.begin(), super_phi.begin(), super_phi.end());
    reference_rho_.insert(reference_rho_.begin(), reference_rho.begin(), reference_rho.end());
    grid_x_.insert(grid_x_.begin(), grid[0].begin(), grid[0].end());
    grid_y_.insert(grid_y_.begin(), grid[1].begin(), grid[1].end());
    grid_z_.insert(grid_z_.begin(), grid[2].begin(), grid[2].end());
    grid_w_.insert(grid_w_.begin(), grid[3].begin(), grid[3].end());

    nso_ = (long int)sqrt(T_.size());
    if ( nso_*nso_ != T_.size() ){
        outfile->Printf("\n");
        outfile->Printf("    error: could not determine nso\n");
        outfile->Printf("\n");
    }

    phi_points_ = reference_rho.size();

    // opdm
    D1a_     = (double*)malloc(nso_ * nso_ * sizeof(double));
    D1b_     = (double*)malloc(nso_ * nso_ * sizeof(double));

    // extra containers
    int dim = phi_points_ > nso_*nso_ ? phi_points_ : nso_*nso_;

    temp1_ = (double*)malloc(dim*sizeof(double));
    temp2_ = (double*)malloc(dim*sizeof(double));

    // options
    r_convergence_ = options_.get_double("R_CONVERGENCE");
    e_convergence_ = options_.get_double("E_CONVERGENCE");
    lbfgs_maxiter_ = options_.get_int("LBFGS_MAXITER");
    do_check_derivatives_ = options_.get_bool("CHECK_DERIVATIVES");

    constrained_search();
}


cs_solver::~cs_solver()
{
    grid_x_.clear();
    grid_y_.clear();
    grid_z_.clear();
    grid_w_.clear();

    free(D1a_);
    free(D1b_);

    free(temp1_);
    free(temp2_);

    free(vars_);
    free(dvars_);

    free(lambda_);
    free(c_);
    free(cval_);
    free(cerror_);
}

double cs_solver::constrained_search(){

    outfile->Printf("\n\n");
    outfile->Printf( "        ************************************************************************\n");
    outfile->Printf( "        *                                                                      *\n");
    outfile->Printf( "        *        Constrained Search Kohn-Sham Density Functional Theory        *\n");
    outfile->Printf( "        *                                                                      *\n");
    outfile->Printf( "        ************************************************************************\n");
    outfile->Printf("\n\n");

    // initialize penalty parameter 
    mu_ = 0.1;

    n_  = nso_ * (nalpha_ + nbeta_);

    std::string key = "DENSITY_CONSTRAINT";
    std::string value = options_.get_str(key);

    nconstraints_ = 0;
    nconstraints_ += nalpha_ * nalpha_; // orthonormalize alpha orbitals
    nconstraints_ += nbeta_  * nbeta_;  // orthonormalize beta orbitals
    // constrain rho
    if (value == "POINTWISE" ) {
        nconstraints_ += phi_points_;
    }else{
        nconstraints_ += 1;
    }

    outfile->Printf("\n");
    outfile->Printf("    ==> Optimization Details <==\n");
    outfile->Printf("\n");
    outfile->Printf("        r_convergence:                 %20.2le\n",r_convergence_);
    outfile->Printf("        e_convergence:                 %20.2le\n",e_convergence_);
    outfile->Printf("        lbfgs_maxiter:                 %20li\n",lbfgs_maxiter_);
    outfile->Printf("        total number of variables:     %20li\n",n_);
    outfile->Printf("        number of grid points:         %20li\n",phi_points_);
    outfile->Printf("        total number of constraints:   %20li\n",nconstraints_);
    outfile->Printf("        density constraint:            %20s\n",options_.get_str(key).c_str());
    //outfile->Printf("        reference kinetic energy:      %20.6lf\n",eref_);
    outfile->Printf("\n");

    lambda_ = (double*)malloc(nconstraints_*sizeof(double));
    c_      = (double*)malloc(nconstraints_*sizeof(double));
    cval_   = (double*)malloc(nconstraints_*sizeof(double));
    cerror_ = (double*)malloc(nconstraints_*sizeof(double));

    memset((void*)lambda_,'\0',nconstraints_*sizeof(double));
    memset((void*)c_,'\0',nconstraints_*sizeof(double));
    memset((void*)cval_,'\0',nconstraints_*sizeof(double));
    memset((void*)cerror_,'\0',nconstraints_*sizeof(double));

    vars_  = lbfgs_malloc(n_);
    dvars_ = lbfgs_malloc(n_);

    memset((void*)vars_,'\0',n_*sizeof(double));
    memset((void*)dvars_,'\0',n_*sizeof(double));

    // orbitals
    R1da_  = vars_;
    R1db_  = vars_ + nso_ * nalpha_;

    // gradient of objective function with respect to orbitals
    dR1da_ = dvars_;
    dR1db_ = dvars_ + nso_ * nalpha_;

    // set constraint values
    set_constraints();

    // random guess
    random_guess();

    // check derivatives:
    if ( do_check_derivatives_ ) {
        check_derivatives(vars_);
        return 0.0;
    }
    
    outfile->Printf("\n");
    outfile->Printf("    Guess energy: %20.12lf\n", kinetic_energy());
    outfile->Printf("\n");
    outfile->Printf("    =================================================================================\n");
    outfile->Printf("    oiter");
    outfile->Printf("  iiter");
    outfile->Printf("                              L");
    outfile->Printf("           E");
    outfile->Printf("   max error");
    outfile->Printf("    error norm\n");
    outfile->Printf("    ---------------------------------------------------------------------------------\n");

    // convergence
    double conv = 100.0;

    // outer iterations:
    int oiter = 0;

    // maximum violation in constraints
    double max = 1000.0;

    double current_energy = 0.0;
    double old_energy = 1.0;

    do {

        old_energy = current_energy;

        iiter_ = 0;

        lbfgsfloatval_t f = evaluate_augmented_lagrangian();
        lbfgs_parameter_t param;
        lbfgs_parameter_init(&param);

        // adjust default lbfgs parameters
        param.max_iterations = lbfgs_maxiter_;
        //param.epsilon = 1e-8;

        // minimize lagrangian function
        lbfgs(n_, vars_, &f, lbfgs_evaluate, progress, (void*)this, &param);

        f = evaluate_augmented_lagrangian();

        // update lagrange multipliers and penalty parameter
        double newmax = 0.0;
        int imax = -1;
        for (int i = 0; i < nconstraints_; i++) {
            if ( fabs(cerror_[i]) > newmax ) {
                newmax = fabs(cerror_[i]);
                imax = i;
             }
        }
        //if ( newmax < 0.25 * max && oiter > 0 ){
        //    f = evaluate_augmented_lagrangian();
        //    for (int i = 0; i < nconstraints_; i++) {
        //        lambda_[i] -= cerror_[i] / mu_;
        //    }
        //}else{
            f = evaluate_augmented_lagrangian();
            for (int i = 0; i < nconstraints_; i++) {
                //printf("%5i %20.12lf\n", i, cval_[i]);fflush(stdout);
                lambda_[i] -= cerror_[i] / mu_;
            }
            //mu_ *= 0.9;
            mu_ /= (10.0 + ( (double)rand()/RAND_MAX - 0.5 ) * 4.0);
        //}
        max = newmax;

        oiter++;
        current_energy = energy_;
        conv = C_DNRM2(nconstraints_,cerror_,1);

        build_d1();
        offset_ = 0;

        outfile->Printf("      %3i %6i %30.6lf %11.6lf %11.3le   %11.3le\n",oiter,iiter_,f,current_energy,max,conv);
        if (PyErr_CheckSignals() != 0) {
            throw pybind11::error_already_set();
        }
        fflush(stdout);

        // for correlated densities, this procedure might not converge tightly ... quit when penalty parameter gets small
        if ( mu_ < 1e-20 ) break;

    }while(fabs(current_energy - old_energy) > e_convergence_ || conv > r_convergence_);

    outfile->Printf("    =================================================================================\n");
    if ( mu_ >= 1e-20 ) {
        outfile->Printf("\n");
        outfile->Printf("    CS iterations converged!\n");
        outfile->Printf("\n");
    }else {
        outfile->Printf("\n");
        outfile->Printf("    CS iterations did not converge.\n");
        outfile->Printf("\n");
    }

    memset((void*)dvars_,'\0',n_*sizeof(double));
    evaluate_augmented_lagrangian();

    //outfile->Printf("        Reference kinetic energy: %20.12lf\n",eref_);
    outfile->Printf("        KS-DFT kinetic energy:    %20.12lf\n",energy_);
    outfile->Printf("\n");

    return energy_;
}

void cs_solver::set_constraints() {

    offset_ = 0;

    // orthonormalize alpha orbitals
    for (int i = 0; i < nalpha_; i++) {
        for (int j = 0; j < nalpha_; j++) {
            c_[offset_++] = (double)(i==j);
        }
    }
    // orthonormalize beta orbitals
    for (int i = 0; i < nbeta_; i++) {
        for (int j = 0; j < nbeta_; j++) {
            c_[offset_++] = (double)(i==j);
        }
    }

    // fixed rho
    std::string key = "DENSITY_CONSTRAINT";
    std::string value = options_.get_str(key);
    if (value == "POINTWISE" ) {
        for (int p = 0; p < phi_points_; p++) {
            c_[offset_++] += reference_rho_[p];
        }
    }else{
        c_[offset_++] = 0.0;
    }
}

void cs_solver::build_d1(){

    F_DGEMM('n', 't', nso_, nso_, nalpha_, 1.0, R1da_, nso_, R1da_, nso_, 0.0, D1a_, nso_);
    F_DGEMM('n', 't', nso_, nso_, nbeta_ , 1.0, R1db_, nso_, R1db_, nso_, 0.0, D1b_, nso_);
}

void cs_solver::evaluate_rho_pointwise(){

    // rho = rho,ref

    C_DCOPY(nso_*nso_,D1a_,1,temp1_,1);
    C_DAXPY(nso_*nso_,1.0,D1b_,1,temp1_,1);

    F_DGEMM('n','n',nso_,phi_points_,nso_,1.0,temp1_,nso_,&(super_phi_[0]),nso_,0.0,&(temp_phi_[0]),nso_);

    // evaluate rho at each grid point
    for (int p = 0; p < phi_points_; p++) {
        double dum = C_DDOT(nso_,&(temp_phi_[p*nso_]),1,&(super_phi_[p*nso_]),1);
        cval_[offset_+p] = dum;
    }

    // evaluate constraints

    C_DCOPY(phi_points_,cval_+offset_,1,cerror_+offset_,1);
    C_DAXPY(phi_points_,-1.0,c_+offset_,1,cerror_+offset_,1);

    // evaluate derivative of constraints
    for (int i = 0; i < nso_; i++) {
        for (int j = 0; j < nso_; j++) {
            double dum = 0.0;
            for (int p = 0; p < phi_points_; p++) {
                dum += (2.0/mu_ * cerror_[offset_+p] - lambda_[offset_+p]) * super_phi_[p*nso_ + i] * super_phi_[p*nso_ + j];
            }
            temp1_[i * nso_ + j] = 2.0 * dum;
        }
    }
    F_DGEMM('t', 'n', nso_, nalpha_, nso_, 1.0, temp1_, nso_, R1da_, nso_, 1.0, dR1da_, nso_);
    F_DGEMM('t', 'n', nso_, nbeta_, nso_,  1.0, temp1_, nso_, R1db_, nso_, 1.0, dR1db_, nso_);

    offset_ += phi_points_;
}

void cs_solver::evaluate_integral_rho_squared(){

    // int (p-p,ref)^2 dr = 0

    C_DCOPY(nso_*nso_,D1a_,1,temp1_,1);
    C_DAXPY(nso_*nso_,1.0,D1b_,1,temp1_,1);

    F_DGEMM('n','n',nso_,phi_points_,nso_,1.0,temp1_,nso_,&(super_phi_[0]),nso_,0.0,&(temp_phi_[0]),nso_);

    // evaluate int (rho - rho,ref)^2 dr
    cval_[offset_] = 0.0;
    for (int p = 0; p < phi_points_; p++) {
        // (rho - rho,ref ) ^2
        double dum = C_DDOT(nso_,&(temp_phi_[p*nso_]),1,&(super_phi_[p*nso_]),1) - reference_rho_[p];
        cval_[offset_] += dum*dum * grid_w_[p] ;
        // for intermediate, phi'(p,i) = grid_w_ [ rho(p) - reference_rho_(p) ] * phi(p,i)
        temp1_[p] = grid_w_[p] * dum;
    }

    // phi'(p,i) = grid_w_ [ rho(p) - reference_rho_(p) ] * phi(p,i)
    for (int p = 0; p < phi_points_; p++) {
        double dum = C_DDOT(nso_,&(temp_phi_[p*nso_]),1,&(super_phi_[p*nso_]),1) - reference_rho_[p];
        for (int i = 0; i < nso_; i++) {
            temp_phi_[p*nso_+i]  = temp1_[p] * super_phi_[p*nso_+i];
        }
    }

    // square root ...
    std::string key = "DENSITY_CONSTRAINT";
    if ( options_.get_str(key) == "L2_NORM" ) {
        cval_[offset_] = sqrt(cval_[offset_]);
    }

    // evaluate error in constraints
    cerror_[offset_] = cval_[offset_] - c_[offset_];

    // evaluate derivative of constraints ... note extra 1/2 / cval for square root
    F_DGEMM('n', 't', nso_, nso_, phi_points_,  4.0, &(temp_phi_[0]), nso_, &(super_phi_[0]), nso_, 0.0, temp1_, nso_);
    if ( options_.get_str(key) == "L2_NORM" ) {
        F_DGEMM('t', 'n', nso_, nalpha_, nso_, 0.5 / cval_[offset_] * (2.0/mu_ * cerror_[offset_] - lambda_[offset_]), temp1_, nso_, R1da_, nso_, 1.0, dR1da_, nso_);
        F_DGEMM('t', 'n', nso_, nbeta_, nso_,  0.5 / cval_[offset_] * (2.0/mu_ * cerror_[offset_] - lambda_[offset_]), temp1_, nso_, R1db_, nso_, 1.0, dR1db_, nso_);
    }else {
        F_DGEMM('t', 'n', nso_, nalpha_, nso_, (2.0/mu_ * cerror_[offset_] - lambda_[offset_]), temp1_, nso_, R1da_, nso_, 1.0, dR1da_, nso_);
        F_DGEMM('t', 'n', nso_, nbeta_, nso_, (2.0/mu_ * cerror_[offset_] - lambda_[offset_]), temp1_, nso_, R1db_, nso_, 1.0, dR1db_, nso_);
    }

    offset_++;
}

double cs_solver::evaluate_augmented_lagrangian(){

    memset((void*)dvars_,'\0',n_*sizeof(double));

    build_d1();
    energy_ = kinetic_energy();

    offset_ = 0;

    evaluate_orbital_orthonormality();

    std::string key = "DENSITY_CONSTRAINT";
    std::string value = options_.get_str(key);

    if (value == "POINTWISE" ) {
        evaluate_rho_pointwise();
    }else{
        evaluate_integral_rho_squared();
    }

    C_DCOPY(nconstraints_, cval_, 1, cerror_, 1);
    C_DAXPY(nconstraints_, -1.0, c_, 1, cerror_, 1);

    double f = energy_;
    for (int i = 0; i < nconstraints_; i++) {
         f += - lambda_[i] * cerror_[i] + 1.0 / mu_ * cerror_[i]*cerror_[i];
    }
    return f;
}

void cs_solver::evaluate_orbital_orthonormality() {

    // alpha
    for (int i = 0; i < nalpha_; i++) {
        for (int j = 0; j < nalpha_; j++) {

            double dum = C_DDOT(nso_, R1da_ + i * nso_, 1, R1da_ + j * nso_, 1);

            cval_[offset_ + i * nalpha_ + j]  = dum;
            cerror_[offset_ + i * nalpha_ + j] = cval_[offset_ + i * nalpha_ + j] - c_[offset_ + i * nalpha_ + j];
            temp1_[i * nalpha_ + j] = 2.0/mu_ *cerror_[offset_ + i * nalpha_ + j] - lambda_[offset_ + i * nalpha_ + j];
        }
    }
    F_DGEMM('n', 't', nso_, nalpha_, nalpha_, 1.0, R1da_, nso_, temp1_, nalpha_, 1.0, dR1da_, nso_);
    F_DGEMM('n', 'n', nso_, nalpha_, nalpha_, 1.0, R1da_, nso_, temp1_, nalpha_, 1.0, dR1da_, nso_);

    offset_ += nalpha_ * nalpha_;

    // beta
    for (int i = 0; i < nbeta_; i++) {
        for (int j = 0; j < nbeta_; j++) {

            double dum = C_DDOT(nso_, R1db_ + i * nso_ , 1, R1db_ + j * nso_, 1);

            cval_[offset_ + i * nbeta_ + j]  = dum;
            cerror_[offset_ + i * nbeta_ + j] = cval_[offset_ + i * nbeta_ + j] - c_[offset_ + i * nbeta_ + j];
            temp1_[i * nbeta_ + j] = 2.0/mu_ *cerror_[offset_ + i * nbeta_ + j] - lambda_[offset_ + i * nbeta_ + j];
        }
    }
    F_DGEMM('n', 't', nso_, nbeta_, nbeta_, 1.0, R1db_, nso_, temp1_, nbeta_, 1.0, dR1db_, nso_);
    F_DGEMM('n', 'n', nso_, nbeta_, nbeta_, 1.0, R1db_, nso_, temp1_, nbeta_, 1.0, dR1db_, nso_);

    offset_ += nbeta_ * nbeta_;
}

void cs_solver::random_guess(){

    // random guess for orbitals

    srand(time(0));
    memset((void*)R1da_,'\0',nso_*nalpha_*sizeof(double));
    memset((void*)R1db_,'\0',nso_*nbeta_*sizeof(double));
    //R1da_[0] = 1.0;
    //R1db_[0] = 1.0;

    for (int i = 0; i < nalpha_; i++) {
        for (int j = 0; j < nso_; j++) {
            double dum = ((double)rand()/RAND_MAX - 0.5 ) * 2.0 / 1e3;
            R1da_[i*nso_+j] = dum;
        }
    }

    for (int i = 0; i < nbeta_; i++) {
        for (int j = 0; j < nso_; j++) {
            double dum = ((double)rand()/RAND_MAX - 0.5 ) * 2.0 / 1e3;
            R1db_[i*nso_+j] = dum;
        }
    }

    build_d1();
}

double cs_solver::kinetic_energy(){

    // evaluate kinetic energy:

    double ena = C_DDOT(nso_*nso_, D1a_, 1, &(T_[0]), 1);
    double enb = C_DDOT(nso_*nso_, D1b_, 1, &(T_[0]), 1);
    
    double en = ena + enb;
    
    // derivative of kinetic energy:
    F_DGEMM('n', 'n', nso_, nalpha_, nso_, 2.0, &(T_[0]),nso_, R1da_, nso_, 1.0, dR1da_, nso_);
    F_DGEMM('n', 'n', nso_, nbeta_, nso_, 2.0, &(T_[0]), nso_, R1db_, nso_, 1.0, dR1db_, nso_);

    return en;
}

void cs_solver::check_derivatives(lbfgsfloatval_t * vars) {

    // displacement
    double h = 1e-4;

    evaluate_augmented_lagrangian();

    double * temp3 = (double*)malloc(n_*sizeof(double));
    C_DCOPY(n_,dvars_,1,temp3,1);

    for (int i = 0; i < n_; i++) {
        double save = vars[i];
        double refder  = temp3[i];
        vars[i] = save + h;
        double dump1 = evaluate_augmented_lagrangian();
        vars[i] = save + 2.0 * h;
        double dump2 = evaluate_augmented_lagrangian();
        vars[i] = save - h;
        double dumm1 = evaluate_augmented_lagrangian();
        vars[i] = save - 2.0 * h;
        double dumm2 = evaluate_augmented_lagrangian(); 
        double myder = (-dump2 + 8.0 * dump1 - 8.0 * dumm1 + dumm2) / (12.0 * h);
        vars[i] = save;
        outfile->Printf("%5i %20.12le %20.12le %20.12le %20.12le %30.12lf %30.12lf %30.12lf\n",i,dump2,dump1,dumm1,dumm2,refder,myder,(refder-myder)/refder * 100.0);
    }

    free(temp3);
}

}
