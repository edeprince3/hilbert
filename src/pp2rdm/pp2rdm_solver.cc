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

#include <psi4/libmints/writer.h>
#include <psi4/libmints/writer_file_prefix.h>

#include <psi4/libtrans/integraltransform.h>
#include <psi4/libtrans/mospace.h>

#include <psi4/libplugin/plugin.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libqt/qt.h>

#include <psi4/libpsi4util/process.h>
#include <psi4/libpsi4util/PsiOutStream.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libmints/wavefunction.h>
#include <psi4/psifiles.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libiwl/iwl.h>
#include <psi4/libmints/mintshelper.h>
#include <psi4/libmints/molecule.h>
#include <psi4/libmints/factory.h>
#include <psi4/libmints/basisset.h>
#include <psi4/libmints/local.h>

#include "diis.h"
#include "pp2rdm_solver.h"

#include <focas/focas_c_interface.h>
#include <misc/misc.h>
#include <misc/omp.h>
#include <misc/blas.h>

using namespace psi;
using namespace fnocc;

namespace psi{ namespace pp2rdm{

pp2RDMSolver::pp2RDMSolver(SharedWavefunction reference_wavefunction,Options & options):
    Wavefunction(options){
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

pp2RDMSolver::~pp2RDMSolver()

{
    free(orbital_lagrangian_);
}

void  pp2RDMSolver::common_init(){

    restart_from_checkpoint_file_ = false;

    if ( options_.get_str("P2RDM_TYPE") == "CCD" ) {
        if ( options_.get_str("P2RDM_ALGORITHM") != "PROJECTION" ) {
            throw PsiException("pccd only works with p2rdm_algorithm = projection",__FILE__,__LINE__);
        }
    }

    is_df_ = false;

    if ( options_.get_str("SCF_TYPE") == "DF" || options_.get_str("SCF_TYPE") == "CD" ) {
        is_df_ = true;
    }

    if ( !is_df_ ) {
        throw PsiException("plugin pp2rdm only works with scf_type = df or cd, for now",__FILE__,__LINE__);
    }

    shallow_copy(reference_wavefunction_);

    escf_     = reference_wavefunction_->energy();
    nalpha_   = reference_wavefunction_->nalpha();
    nbeta_    = reference_wavefunction_->nbeta();
    nalphapi_ = reference_wavefunction_->nalphapi();
    nbetapi_  = reference_wavefunction_->nbetapi();
    doccpi_   = reference_wavefunction_->doccpi();
    soccpi_   = reference_wavefunction_->soccpi();
    frzcpi_   = reference_wavefunction_->frzcpi();
    frzvpi_   = reference_wavefunction_->frzvpi();
    nmopi_    = reference_wavefunction_->nmopi();
    nirrep_   = reference_wavefunction_->nirrep();
    if ( nirrep_ > 1 ) {
        throw PsiException("plugin pp2rdm works only with c1 symmetry",__FILE__,__LINE__);
    }
    nso_      = reference_wavefunction_->nso();
    nmo_      = reference_wavefunction_->nmo();

    nsopi_    = reference_wavefunction_->nsopi();
    molecule_ = reference_wavefunction_->molecule();
    enuc_     = molecule_->nuclear_repulsion_energy({0.0,0.0,0.0});

    // were there linear dependencies in the primary basis set?
    if ( nmo_ != nso_ ) {
        throw PsiException("linear dependencies.  not sure what to do about that just yet",__FILE__,__LINE__);
    }

    AO2SO_ = SharedMatrix(reference_wavefunction_->aotoso());

    //Ca_  = (SharedMatrix)(new Matrix(reference_wavefunction_->Ca()));
    //Cb_  = (SharedMatrix)(new Matrix(reference_wavefunction_->Cb()));

    Ca_ = SharedMatrix(reference_wavefunction_->Ca());
    Cb_ = SharedMatrix(reference_wavefunction_->Cb());

    if ( options_.get_bool("LOCALIZE_ORBITALS") ) {
        // localize orbitals:
        std::shared_ptr<PMLocalizer> boys (new PMLocalizer(reference_wavefunction_->basisset(),reference_wavefunction_->Ca_subset("SO","OCC")));
        boys->localize();
        for (int mu = 0; mu < nso_; mu++) {
            for (int i = 0; i < nalpha_; i++) {
                Ca_->pointer()[mu][i] = boys->L()->pointer()[mu][i];
                Cb_->pointer()[mu][i] = boys->L()->pointer()[mu][i];
            }
        }
    }
    if ( options_.get_bool("LOCALIZE_VIRTUAL_ORBITALS") ) {
        // localize orbitals (virtual):
        std::shared_ptr<PMLocalizer> boys_vir (new PMLocalizer(reference_wavefunction_->basisset(),reference_wavefunction_->Ca_subset("SO","VIR")));
        boys_vir->localize();
        for (int mu = 0; mu < nso_; mu++) {
            for (int i = nalpha_; i < nso_; i++) {
                Ca_->pointer()[mu][i] = boys_vir->L()->pointer()[mu][i - nalpha_];
                Cb_->pointer()[mu][i] = boys_vir->L()->pointer()[mu][i - nalpha_];
            }
        }
    }

    if ( options_.get_bool("NOISY_ORBITALS") ) {

        // add noise
        // i' =  cos i + sin j
        // j' = -sin i + cos j
        srand(clock());

        int count = 0;
        do {

            int p = (int)( (double)rand()/RAND_MAX * nso_);
            int q = (int)( (double)rand()/RAND_MAX * nso_);
            if ( p == q ) continue;

            double angle = ( (double)rand()/RAND_MAX - 1.0 ) * 0.001;
            double cosine = cos(angle);
            double sine   = sin(angle);

            for (int mu = 0; mu < nso_; mu++) {

                double oldp = Ca_->pointer()[mu][p];
                double oldq = Ca_->pointer()[mu][q];

                double newp =  cosine * oldp +   sine * oldq;
                double newq =   -sine * oldp + cosine * oldq;

                Ca_->pointer()[mu][p] = newp;
                Ca_->pointer()[mu][q] = newq;

            }

            count++;

        }while(count < nso_*nso_*10);

        Cb_->copy(Ca_);

    }

    S_  = (SharedMatrix)(new Matrix(reference_wavefunction_->S()));

    Fa_  = SharedMatrix(reference_wavefunction_->Fa());
    Fb_  = SharedMatrix(reference_wavefunction_->Fb());

    Da_  = SharedMatrix(reference_wavefunction_->Da());
    Db_  = SharedMatrix(reference_wavefunction_->Db());
    
    // Lagrangian matrix
    Lagrangian_ = SharedMatrix(reference_wavefunction_->Lagrangian());

    epsilon_a_= SharedVector(new Vector(nirrep_, nmopi_));
    epsilon_a_->copy(reference_wavefunction_->epsilon_a().get());
    epsilon_b_= SharedVector(new Vector(nirrep_, nmopi_));
    epsilon_b_->copy(reference_wavefunction_->epsilon_b().get());

    // memory is from process::environment
    memory_ = Process::environment.get_memory();

    // set the wavefunction name
    name_ = options_.get_str("P2RDM_TYPE");
    if ( name_ == "K" ) {
        name_ = "p2RDM";
    }

    // pp2rdm convergence thresholds:
    r_convergence_  = options_.get_double("R_CONVERGENCE");
    e_convergence_  = options_.get_double("E_CONVERGENCE");

    maxiter_        = options_.get_int("MAXITER");

    // memory check happens here

    outfile->Printf("\n\n");
    outfile->Printf( "        ***************************************************\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *    pp2RDM                                       *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *    pair parametric 2RDM                         *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        ***************************************************\n");

    outfile->Printf("\n");
    outfile->Printf("    ==> Convergence parameters <==\n");
    outfile->Printf("\n");
    outfile->Printf("        r_convergence:                  %5.3le\n",r_convergence_);
    outfile->Printf("        e_convergence:                  %5.3le\n",e_convergence_);
    outfile->Printf("        maximum iterations:                 %5i\n",options_.get_int("MAXITER"));
    outfile->Printf("\n");

    outfile->Printf("    ==> Orbital optimization parameters <==\n");
    outfile->Printf("\n");
    outfile->Printf("        localize orbitals?                  %5s\n",options_.get_bool("LOCALIZE_ORBITALS") ? "yes" : "no");
    outfile->Printf("        optimize orbitals?                  %5s\n",options_.get_bool("OPTIMIZE_ORBITALS") ? "yes" : "no");
    outfile->Printf("        g_convergence:                  %5.3le\n",options_.get_double("ORBOPT_GRADIENT_CONVERGENCE"));
    outfile->Printf("        e_convergence:                  %5.3le\n",options_.get_double("ORBOPT_ENERGY_CONVERGENCE"));
    outfile->Printf("        maximum iterations:                 %5i\n",options_.get_int("ORBOPT_MAXITER"));
    //outfile->Printf("        maximum iterations:               dynamic\n");
    outfile->Printf("        exact diagonal Hessian:             %5s\n",options_.get_bool("ORBOPT_EXACT_DIAGONAL_HESSIAN") ? "true" : "false");
    outfile->Printf("        number of DIIS vectors:             %5i\n",options_.get_int("ORBOPT_NUM_DIIS_VECTORS"));
    outfile->Printf("        print iteration info:               %5s\n",options_.get_bool("ORBOPT_WRITE") ? "true" : "false");

    //outfile->Printf("\n");
    //outfile->Printf("  ==> Memory requirements <==\n");
    //outfile->Printf("\n");

    // if using 3-index integrals, transform them before allocating any memory integrals, transform
    if ( is_df_ ) {
        outfile->Printf("\n");
        outfile->Printf("    ==> Transform three-electron integrals <==\n");
        outfile->Printf("\n");

        double start = omp_get_wtime();

        ::ThreeIndexIntegrals(reference_wavefunction_,nQ_,memory_);

        Qmo_ = (double*)malloc(nmo_*(nmo_+1)/2*nQ_*sizeof(double));
        memset((void*)Qmo_,'\0',nmo_*(nmo_+1)/2*nQ_*sizeof(double));

        std::shared_ptr<PSIO> psio(new PSIO());
        psio->open(PSIF_DCC_QMO,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_QMO,"(Q|mn) Integrals",(char*)Qmo_,sizeof(double)*nQ_ * nmo_*(nmo_+1)/2);
        psio->close(PSIF_DCC_QMO,1);

        double end = omp_get_wtime();

        outfile->Printf("        Time for integral transformation:  %7.2lf s\n",end-start);
        outfile->Printf("\n");
    } else {
        // transform integrals
        outfile->Printf("    ==> Transform two-electron integrals <==\n");
        outfile->Printf("\n");

        double start = omp_get_wtime();
        std::vector<std::shared_ptr<MOSpace> > spaces;
        spaces.push_back(MOSpace::all);
        std::shared_ptr<IntegralTransform> ints(new IntegralTransform(reference_wavefunction_, spaces,
            IntegralTransform::TransformationType::Restricted, IntegralTransform::OutputType::IWLOnly,
            IntegralTransform::MOOrdering::PitzerOrder, IntegralTransform::FrozenOrbitals::None, false));
        ints->set_dpd_id(0);
        ints->set_keep_iwl_so_ints(true);
        ints->set_keep_dpd_so_ints(true);
        ints->initialize();
        ints->transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
        double end = omp_get_wtime();
        outfile->Printf("\n");
        outfile->Printf("        Time for integral transformation:  %7.2lf s\n",end-start);
        outfile->Printf("\n");

    }

    // one-electron hamiltonian:
    std::shared_ptr<MintsHelper> mints(new MintsHelper(reference_wavefunction_));
    std::shared_ptr<Matrix> Hcore = (std::shared_ptr<Matrix>)(new Matrix(mints->so_potential()));
    Hcore->add(mints->so_kinetic());
    Hcore->transform(Ca_);

    oei_dim_ = nmo_*(nmo_+1)/2;
    oei_ = (double*)malloc(oei_dim_*sizeof(double));
    for (long int i = 0; i < nmo_; i++) {
        for (long int j = i; j < nmo_; j++) {
                oei_[INDEX(i,j)] = Hcore->pointer()[i][j];
        }
    }

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    // integrals

    v_iiaa_ = (double*)malloc(o * v * sizeof(double)); // (ii|aa)
    v_iaia_ = (double*)malloc(o * v * sizeof(double)); // (ia|ai)
    v_abab_ = (double*)malloc(v * v * sizeof(double)); // (ab|ab)
    v_ijij_ = (double*)malloc(o * o * sizeof(double)); // (ij|ij)

    memset((void*)v_iiaa_,'\0',o * v * sizeof(double));
    memset((void*)v_iaia_,'\0',o * v * sizeof(double));
    memset((void*)v_abab_,'\0',v * v * sizeof(double));
    memset((void*)v_ijij_,'\0',o * o * sizeof(double));

    // fock matrix

    fo_ = (double*)malloc(o * sizeof(double));
    fv_ = (double*)malloc(v * sizeof(double));

    // amplitudes
    r2_ = (double*)malloc(o * v * sizeof(double)); // residual
    t2_ = (double*)malloc(o * v * sizeof(double)); // amplitudes
    z2_ = (double*)malloc(o * v * sizeof(double)); // z amplitudes (ccd only)
    t0_ = (double*)malloc(o * v * sizeof(double)); // normalization coefficients

    memset((void*)r2_,'\0',o * v * sizeof(double));
    memset((void*)z2_,'\0',o * v * sizeof(double));
    memset((void*)t2_,'\0',o * v * sizeof(double));
    memset((void*)t0_,'\0',o * v * sizeof(double));

    // stuff for orbital optimization:

    // tei
    tei_dim_ = (long int) nQ_ * (long int) nmo_ * ( (long int) ( nmo_ + 1L ) ) / 2L ;

    tei_ = Qmo_;

    // d2
    d2_dim_ = nmo_ * ( nmo_ + 1 ) / 2  * ( nmo_ * ( nmo_ + 1 ) / 2 + 1 ) / 2;
    d2_ = (double*)malloc(d2_dim_ * sizeof(double));
    memset((void*)d2_,'\0',d2_dim_ * sizeof(double));

    // mo-mo transformation matrix
    newMO_ = (SharedMatrix)(new Matrix(reference_wavefunction_->Ca()));
    newMO_->zero();
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < nmopi_[h]; i++) {
            newMO_->pointer(h)[i][i] = 1.0;
        }
    }

    orbopt_transformation_matrix_ = (double*)malloc(nmo_ * nmo_ * sizeof(double));
    memset((void*)orbopt_transformation_matrix_,'\0',nmo_ * nmo_ * sizeof(double));
    for (int i = 0; i < nmo_; i++) {
        orbopt_transformation_matrix_[i * nmo_ + i] = 1.0;
    }

    int nthread = omp_get_max_threads();

    orbopt_data_    = (double*)malloc(15*sizeof(double));
    orbopt_data_[0] = (double)nthread;
    orbopt_data_[1] = 1.0; // include active-active rotations
    orbopt_data_[2] = (double)nfrzc_; //(double)options_.get_int("ORBOPT_FROZEN_CORE");
    orbopt_data_[3] = (double)options_.get_double("ORBOPT_GRADIENT_CONVERGENCE");
    orbopt_data_[4] = (double)options_.get_double("ORBOPT_ENERGY_CONVERGENCE");
    orbopt_data_[5] = (double)(options_.get_bool("ORBOPT_WRITE") ? 1.0 : 0.0 );
    orbopt_data_[6] = (double)(options_.get_bool("ORBOPT_EXACT_DIAGONAL_HESSIAN") ? 1.0 : 0.0 );
    orbopt_data_[7] = (double)options_.get_int("ORBOPT_NUM_DIIS_VECTORS");
    orbopt_data_[8] = (double)options_.get_int("ORBOPT_MAXITER");
    orbopt_data_[9] = 1.0;   // use density fitting
    orbopt_data_[10] = 0.0;  // number of iterations (output)
    orbopt_data_[11] = 0.0;  // gradient norm (output)
    orbopt_data_[12] = 0.0;  // change in energy (output)
    orbopt_data_[13] = 0.0;  // converged?

    // orbital optimizatoin algorithm
    //orbopt_data_[14] = 1.0;
    //if      ( options_.get_str("ORBOPT_ALGORITHM") == "QUASI_NEWTON" )       orbopt_data_[14] = 0.0;
    //else if ( options_.get_str("ORBOPT_ALGORITHM") == "CONJUGATE_GRADIENT" ) orbopt_data_[14] = 1.0;
    //else if ( options_.get_str("ORBOPT_ALGORITHM") == "NEWTON_RAPHSON" )     orbopt_data_[14] = 2.0;
    // orbital optimizatoin algorithm

    // orbital optimization algorithm is handled automatically internally
    orbopt_data_[14] = 0.0;
    //orbopt_data_[14] = 3.0;
    //if      ( options_.get_str("ORBOPT_ALGORITHM") == "STEEPEST_DESCENT" ) orbopt_data_[14] = 0.0;
    //else if ( options_.get_str("ORBOPT_ALGORITHM") == "HESTENES_STIEFEL" ) orbopt_data_[14] = 1.0;
    //else if ( options_.get_str("ORBOPT_ALGORITHM") == "DAI_YUAN" )         orbopt_data_[14] = 2.0;
    //else if ( options_.get_str("ORBOPT_ALGORITHM") == "HAGER_ZHANG" )      orbopt_data_[14] = 3.0;
    //else if ( options_.get_str("ORBOPT_ALGORITHM") == "KOU_DAI" )          orbopt_data_[14] = 4.0;

    orbopt_converged_ = false;

    // don't change the length of this filename
    orbopt_outfile_ = (char*)malloc(120*sizeof(char));
    std::string filename = get_writer_file_prefix(reference_wavefunction_->molecule()->name()) + ".orbopt";
    strcpy(orbopt_outfile_,filename.c_str());
    if ( options_.get_bool("ORBOPT_WRITE") ) {
        FILE * fp = fopen(orbopt_outfile_,"w");
        fclose(fp);
    }

    same_a_b_orbs_ = true;
    same_a_b_dens_ = true;

//AED
    orbital_lagrangian_ = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
}

// compute the energy!
double pp2RDMSolver::compute_energy() {

    int iter = 0;

    d1_ = (double*)malloc( nmo_ * (nmo_+1)/2 * sizeof(double));
    d2ab_ = (double*)malloc(nmo_ * nmo_ * nmo_ * nmo_ * sizeof(double));
    d2aa_ = (double*)malloc(nmo_ * nmo_ * nmo_ * nmo_ * sizeof(double));

    //outfile->Printf("\n");
    outfile->Printf("    ==> Begin oo-p%s iterations <==\n",name_.c_str());
    outfile->Printf("\n");

    double dE = 0;
    double energy = -enuc_;
    do {

        ci_iter_ = 0;

        double ci_energy = 0.0;
        if ( options_.get_str("P2RDM_ALGORITHM") == "PROJECTION" ) {

            ci_energy = pp2rdm_iterations(ci_iter_);

            if ( options_.get_str("P2RDM_TYPE") == "CCD" ) {
                int z_iter = 0;
                z_iterations(z_iter);
            }

        }else if ( options_.get_str("P2RDM_ALGORITHM") == "LBFGS" ) {

            //throw PsiException("P2RDM_ALGORITHM LBFGS is currently disabled.",__FILE__,__LINE__);
            ci_energy = pp2rdm_lbfgs_iterations(ci_iter_);

        }else if ( options_.get_str("P2RDM_ALGORITHM") == "NEWTON_RAPHSON" ) {

            ci_energy = pp2rdm_newton_raphson_iterations(ci_iter_);

        }

        // extra logic because CCD RDMs are different
        double rdm_energy = 0.0;
        if ( options_.get_str("P2RDM_TYPE") != "CCD" ) {
            rdm_energy = BuildSeniorityZeroRDMs(true);
        }else {
            rdm_energy = BuildSeniorityZeroRDMs_pCCD(true);
        }

        double oo_energy;
        if ( options_.get_bool("OPTIMIZE_ORBITALS") ) {

            if ( iter < 50 ) orbopt_data_[14] = 3.0; // steepest descent
            else             orbopt_data_[14] = 3.0; // hager-zhang 

            oo_energy = RotateOrbitals();
            dE = fabs(energy - oo_energy);

            // increase number of orbital optimization steps
            orbopt_data_[8] += 1.0; // AED

        }else {
            oo_energy = rdm_energy;
            orbopt_converged_ = true;
            dE = 0.0;
        }

        energy = oo_energy;

        if ( options_.get_bool("OPTIMIZE_ORBITALS") ) {
            outfile->Printf("\n");
            outfile->Printf("    ==> oo-p%s summary <==\n",name_.c_str());
            outfile->Printf("\n");

            outfile->Printf("       iter");
            outfile->Printf("   iter(CI)");
            outfile->Printf("   iter(oo)");
            outfile->Printf("           E(CI)");
            outfile->Printf("           E(oo)");
            outfile->Printf("            |dE|\n");
            outfile->Printf("      %5i %9i %11i %15.10lf %15.10lf %15.10lf\n",iter,ci_iter_,(int)orbopt_data_[10],rdm_energy + enuc_,oo_energy + enuc_,dE);
        }

        if ( orbopt_converged_ && dE < e_convergence_ ) break;

        iter++;

    }while(iter < maxiter_); 

    if ( iter == maxiter_ ) {
        throw PsiException("oo-p2RDM did not converge.",__FILE__,__LINE__);
    }

    if ( options_.get_bool("OPTIMIZE_ORBITALS") ) {
        outfile->Printf("\n");
        outfile->Printf("      oo-p%s iterations converged!\n",name_.c_str());
        outfile->Printf("\n");

        outfile->Printf("    * oo-p%s total energy:              %20.12lf\n",name_.c_str(),energy+enuc_);
        outfile->Printf("\n");
    }else {
        outfile->Printf("\n");
        outfile->Printf("      p%s iterations converged!\n",name_.c_str());
        outfile->Printf("\n");

        outfile->Printf("    * p%s total energy:              %20.12lf\n",name_.c_str(),energy+enuc_);
        outfile->Printf("\n");
    }

    std::shared_ptr<Vector> NOs (new Vector("Natural Orbital Occupation Numbers (spin free)",nmo_));
    BuildSeniorityZeroRDMs(false);
    for (int i = 0; i < nmo_; i++) {
        NOs->pointer()[i] = d1_[INDEX(i,i)];
    }
    NOs->print();

    Process::environment.globals["CURRENT ENERGY"]    = energy + enuc_;
    Process::environment.globals["pp2RDM TOTAL ENERGY"] = energy + enuc_;

    // push final transformation matrix onto Ca_ and Cb_
    if ( options_.get_bool("OPTIMIZE_ORBITALS") ) {
        ::UpdateTransformationMatrix(reference_wavefunction_,newMO_,Ca_,Cb_,orbopt_transformation_matrix_);
    }

    // write tpdm to disk?
    if ( options_.get_bool("TPDM_WRITE_FULL") ) {
        WriteTPDM();
    }

    // print RDM elements to output file?

    if ( options_.get_bool("PRINT_RDMS") ) {
        // extra logic because CCD RDMs are different
        double rdm_energy = 0.0;
        if ( options_.get_str("P2RDM_TYPE") != "CCD" ) {
            rdm_energy = BuildSeniorityZeroRDMs(true);
        }else {
            rdm_energy = BuildSeniorityZeroRDMs_pCCD(true);
        }
        print_rdms();
    }
    if ( options_.get_bool("PRINT_INTEGRALS") ) {
        print_integrals();
    }

    return energy + enuc_;

}

void pp2RDMSolver::print_integrals() {

    // Gamma(pq,rs) = D2ab(pq,rs) + D2ba(pq,rs) + D2aa(pq,rs) + D2bb(pq,rs)
    outfile->Printf("\n");
    outfile->Printf("    ==> @TEI <==\n");
    outfile->Printf("\n");

    double e2 = 0.0; 

    for (int p = 0; p < nmo_; p++) {
        for (int q = 0; q < nmo_; q++) {
            for (int r = 0; r < nmo_; r++) {
                for (int s = 0; s < nmo_; s++) {
                    if ( p == q && r == s ) {
                        //double dum = 2.0 * d2ab_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s];
                        //dum       += 2.0 * d2aa_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s];
                        //outfile->Printf("%5i %5i %5i %5i %20.12lf\n",p,q,r,s,dum);

                        // 0.5 Gamma(pp,rr) (pr|pr)
                        int pr = INDEX(p,r);
                        double eint = C_DDOT(nQ_,Qmo_ + pr,nmo_*(nmo_+1)/2,Qmo_+pr,nmo_*(nmo_+1)/2);
                        //e2 += 0.5 * eint * dum;
                        outfile->Printf("%5i %5i %5i %5i %20.12lf\n",p,q,r,s,eint);

                    }else if ( p == r && q == s ) {
                        //double dum = 2.0 * d2ab_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s];
                        //dum       += 2.0 * d2aa_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s];
                        //outfile->Printf("%5i %5i %5i %5i %20.12lf\n",p,q,r,s,dum);

                        // 0.5 Gamma(pq,pq) (pp|qq)
                        int pp = INDEX(p,p);
                        int qq = INDEX(q,q);
                        double eint = C_DDOT(nQ_,Qmo_ + pp,nmo_*(nmo_+1)/2,Qmo_+qq,nmo_*(nmo_+1)/2);
                        //e2 += 0.5 * eint * dum;
                        outfile->Printf("%5i %5i %5i %5i %20.12lf\n",p,q,r,s,eint);

                    }else if ( p == s && q == r ) {
                        //double dum = 2.0 * d2ab_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s];
                        //dum       += 2.0 * d2aa_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s];
                        //outfile->Printf("%5i %5i %5i %5i %20.12lf\n",p,q,r,s,dum);

                        // 0.5 Gamma(pq,qp) (pq|qp)
                        int pq = INDEX(p,q);
                        double eint = C_DDOT(nQ_,Qmo_ + pq,nmo_*(nmo_+1)/2,Qmo_+pq,nmo_*(nmo_+1)/2);
                        //e2 += 0.5 * eint * dum;
                        outfile->Printf("%5i %5i %5i %5i %20.12lf\n",p,q,r,s,eint);

                    }
                }
            }
        }
    }
    outfile->Printf("\n");
    outfile->Printf("    @END\n");
    outfile->Printf("\n");


    outfile->Printf("\n");
    outfile->Printf("    ==> @OEI <==\n");
    outfile->Printf("\n");

    double e1 = 0.0; 

    for (int p = 0; p < nmo_; p++) {
        //double dum = d1_[INDEX(p,p)];
        //outfile->Printf("%5i %5i %20.12lf\n",p,p,dum);
        //e1 += dum * oei_[INDEX(p,p)];
        outfile->Printf("%5i %5i %20.12lf\n",p,p,oei_[INDEX(p,p)]);
    }
    outfile->Printf("\n");
    outfile->Printf("    @END\n");
    outfile->Printf("\n");

    outfile->Printf("\n");
    outfile->Printf("    ==> @ENUC %20.12lf <==\n",enuc_);
    outfile->Printf("\n");

    //printf("%20.12lf %20.12lf %20.12lf\n",e1,e2,e1+e2+enuc_);

}

void pp2RDMSolver::print_rdms() {

    // Gamma(pq,rs) = D2ab(pq,rs) + D2ba(pq,rs) + D2aa(pq,rs) + D2bb(pq,rs)
    outfile->Printf("\n");
    outfile->Printf("    ==> p%s @2RDM <==\n",name_.c_str());
    outfile->Printf("\n");

    double e2 = 0.0; 

    for (int p = 0; p < nmo_; p++) {
        for (int q = 0; q < nmo_; q++) {
            for (int r = 0; r < nmo_; r++) {
                for (int s = 0; s < nmo_; s++) {
                    if ( p == q && r == s ) {
                        double dum = 2.0 * d2ab_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s];
                        dum       += 2.0 * d2aa_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s];
                        outfile->Printf("%5i %5i %5i %5i %20.12lf\n",p,q,r,s,dum);

                        // 0.5 Gamma(pp,rr) (pr|pr)
                        int pr = INDEX(p,r);
                        double eint = C_DDOT(nQ_,Qmo_ + pr,nmo_*(nmo_+1)/2,Qmo_+pr,nmo_*(nmo_+1)/2);
                        e2 += 0.5 * eint * dum;

                    }else if ( p == r && q == s ) {
                        double dum = 2.0 * d2ab_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s];
                        dum       += 2.0 * d2aa_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s];
                        outfile->Printf("%5i %5i %5i %5i %20.12lf\n",p,q,r,s,dum);

                        // 0.5 Gamma(pq,pq) (pp|qq)
                        int pp = INDEX(p,p);
                        int qq = INDEX(q,q);
                        double eint = C_DDOT(nQ_,Qmo_ + pp,nmo_*(nmo_+1)/2,Qmo_+qq,nmo_*(nmo_+1)/2);
                        e2 += 0.5 * eint * dum;

                    }else if ( p == s && q == r ) {
                        double dum = 2.0 * d2ab_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s];
                        dum       += 2.0 * d2aa_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s];
                        outfile->Printf("%5i %5i %5i %5i %20.12lf\n",p,q,r,s,dum);

                        // 0.5 Gamma(pq,qp) (pq|qp)
                        int pq = INDEX(p,q);
                        double eint = C_DDOT(nQ_,Qmo_ + pq,nmo_*(nmo_+1)/2,Qmo_+pq,nmo_*(nmo_+1)/2);
                        e2 += 0.5 * eint * dum;

                    }
                }
            }
        }
    }
    outfile->Printf("\n");
    outfile->Printf("    @END\n");
    outfile->Printf("\n");


    outfile->Printf("\n");
    outfile->Printf("    ==> p%s @1RDM <==\n",name_.c_str());
    outfile->Printf("\n");

    double e1 = 0.0; 

    for (int p = 0; p < nmo_; p++) {
        double dum = d1_[INDEX(p,p)];
        outfile->Printf("%5i %5i %20.12lf\n",p,p,dum);
        e1 += dum * oei_[INDEX(p,p)];
    }
    outfile->Printf("\n");
    outfile->Printf("    @END\n");
    outfile->Printf("\n");

    //printf("%20.12lf %20.12lf %20.12lf\n",e1,e2,e1+e2+enuc_);

}

void pp2RDMSolver::PackSpatialDensity() {

    int n1 = nmo_;
    int n2 = nmo_ * nmo_ ;
    int n3 = nmo_ * nmo_ * nmo_;

    memset((void*)d2_,'\0',d2_dim_ * sizeof(double));

    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            for (int k = 0; k < nmo_; k++) {
                for (int l = 0; l < nmo_; l++) {

                    int id = INDEX(INDEX(i,k),INDEX(j,l));
                    double val = 0.0;

                    val += 0.5 * d2ab_[i * n3 + j * n2 + k * n1 + l];
                    val += 0.5 * d2ab_[k * n3 + j * n2 + i * n1 + l] * (1.0 - (double)(l==j));
                    val += 0.5 * d2ab_[i * n3 + l * n2 + k * n1 + j] * (1.0 - (double)(i==k));
                    val += 0.5 * d2ab_[k * n3 + l * n2 + i * n1 + j] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));

                    val += 0.5 * d2ab_[j * n3 + i * n2 + l * n1 + k];
                    val += 0.5 * d2ab_[j * n3 + k * n2 + l * n1 + i] * (1.0 - (double)(l==j));
                    val += 0.5 * d2ab_[l * n3 + i * n2 + j * n1 + k] * (1.0 - (double)(i==k));
                    val += 0.5 * d2ab_[l * n3 + k * n2 + j * n1 + i] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));

                    // aa / bb
                    if ( i != j && k != l ) {
                        int sij = 1;//( i < j ? 1 : -1 );
                        int skl = 1;//( k < l ? 1 : -1 );
                        val += 0.5 * sij * skl * d2aa_[i * n3 + j * n2 + k * n1 + l];
                        val += 0.5 * sij * skl * d2aa_[k * n3 + l * n2 + i * n1 + j] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));
                        val += 0.5 * sij * skl * d2aa_[i * n3 + j * n2 + k * n1 + l];
                        val += 0.5 * sij * skl * d2aa_[k * n3 + l * n2 + i * n1 + j] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));
                    }
                    if ( k != j && i != l ) {
                        int sij = 1;//( i < j ? 1 : -1 );
                        int skl = 1;//( k < l ? 1 : -1 );
                        val += 0.5 * sij * skl * d2aa_[k * n3 + j * n2 + i * n1 + l] * (1.0 - (double)(l==j));
                        val += 0.5 * sij * skl * d2aa_[i * n3 + l * n2 + k * n1 + j] * (1.0 - (double)(i==k));
                        val += 0.5 * sij * skl * d2aa_[k * n3 + j * n2 + i * n1 + l] * (1.0 - (double)(l==j));
                        val += 0.5 * sij * skl * d2aa_[i * n3 + l * n2 + k * n1 + j] * (1.0 - (double)(i==k));
                    }

                    // scale the off-diagonal elements
                    if ( INDEX(i,k) != INDEX(j,l) ) {
                        val *= 2.0;
                    }

                    d2_[id] = val;

                }
            }
        }
    }

}

double pp2RDMSolver::RotateOrbitals() {

    PackSpatialDensity();

    if ( orbopt_data_[8] > 0 ) {
        outfile->Printf("\n");
        outfile->Printf("    ==> Orbital Optimization <==\n");
        outfile->Printf("\n");
    }    

    int * symmetry_energy_order  = (int*)malloc(nmo_*sizeof(int));
    for (int i = 0; i < nmo_; i++) {
        symmetry_energy_order[i] = 1;
    }

    int nrstc = 0;
    int nrstv = 0;

    //printf("enter\n");fflush(stdout);
    OrbOpt(orbopt_transformation_matrix_,

          oei_, oei_dim_,

          tei_, tei_dim_,

          d1_, oei_dim_, d2_, d2_dim_,

          symmetry_energy_order, nrstc, nmo_, nrstv, nirrep_,

          orbopt_data_, orbopt_outfile_, orbital_lagrangian_);
    //printf("exit\n");fflush(stdout);

    free(symmetry_energy_order);

    if ( orbopt_data_[8] > 0 ) {
        //outfile->Printf("            Orbital Optimization %s in %3i iterations \n",(int)orbopt_data_[13] ? "converged" : "did not converge",(int)orbopt_data_[10]);
        //outfile->Printf("            Total energy change: %11.6le\n",orbopt_data_[12]);
        //outfile->Printf("            Final gradient norm: %11.6le\n",orbopt_data_[11]);
        //outfile->Printf("\n");

        if ( fabs(orbopt_data_[12]) < orbopt_data_[4] && fabs(orbopt_data_[11]) < orbopt_data_[3] ) {
            orbopt_converged_ = true;
        }
    }

    double e1 = 0.0;

    for (int i = 0; i < nmo_; i++) {
        e1 += d1_[INDEX(i,i)] * oei_[INDEX(i,i)];
    }

    double e2 = 0.0;

    // evaluating two-electron part of energy should be simple, since pair p2rdm looks like doci
    // D2(ij,ij) (ii|jj)
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            if ( i == j ) continue;
            int ii = INDEX(i,i);
            int jj = INDEX(j,j);
            double eint = C_DDOT(nQ_,Qmo_ + ii,nmo_*(nmo_+1)/2,Qmo_+jj,nmo_*(nmo_+1)/2);
            e2 += d2ab_[i * nmo_*nmo_*nmo_ + j * nmo_*nmo_ + i * nmo_ + j] * eint;
            e2 += d2aa_[i * nmo_*nmo_*nmo_ + j * nmo_*nmo_ + i * nmo_ + j] * eint;
        }
    }
    // D2(ij,ji) (ij|ji)
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            if ( i == j ) continue;
            int ij = INDEX(i,j);
            double eint = C_DDOT(nQ_,Qmo_ + ij,nmo_*(nmo_+1)/2,Qmo_+ij,nmo_*(nmo_+1)/2);
            e2 += d2aa_[i * nmo_*nmo_*nmo_ + j * nmo_*nmo_ + j * nmo_ + i] * eint;
        }
    }
    // D2(ii,jj) (ij|ij)
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            int ij = INDEX(i,j);
            double eint = C_DDOT(nQ_,Qmo_ + ij,nmo_*(nmo_+1)/2,Qmo_+ij,nmo_*(nmo_+1)/2);
            e2 += d2ab_[i * nmo_*nmo_*nmo_ + i * nmo_*nmo_ + j * nmo_ + j] * eint;
        }
    }

    //double * ints = (double*)malloc(nmo_*(nmo_+1)/2*nmo_*(nmo_+1)/2*sizeof(double));
    //F_DGEMM('n','t',nmo_*(nmo_+1)/2,nmo_*(nmo_+1)/2,nQ_,1.0,Qmo_,nmo_*(nmo_+1)/2,Qmo_,nmo_*(nmo_+1)/2,0.0,ints,nmo_*(nmo_+1)/2);
    //for (int i = 0; i < nmo_; i++) {
    //    for (int j = 0; j < nmo_; j++) {
    //        for (int k = 0; k < nmo_; k++) {
    //            for (int l = 0; l < nmo_; l++) {

    //                int ij = INDEX(i,j);
    //                int kl = INDEX(k,l);
    //                double eint = ints[ij*nmo_*(nmo_+1)/2+kl];
    //                e2 += d2ab_[i * nmo_*nmo_*nmo_ + k * nmo_*nmo_ + j * nmo_ + l] * eint;
    //                e2 += d2aa_[i * nmo_*nmo_*nmo_ + k * nmo_*nmo_ + j * nmo_ + l] * eint;
    //            }
    //        }
    //    }
    //}
    //free(ints);

    // d2 is being overwritten ...
    //PackSpatialDensity();

    if ( orbopt_data_[8] > 0 ) {
        outfile->Printf("        p2RDM one-electron energy = %20.12lf\n",e1);
        outfile->Printf("        p2RDM two-electron energy = %20.12lf\n",e2);
        outfile->Printf("        * p2RDM total energy      = %20.12lf\n",e1 + e2 + enuc_);
    }

    return e1+e2;

}

double pp2RDMSolver::z_iterations(int & iter) {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    // Initialize our diis extrapolation manager.  Note that,
    // in this implemenation, the DIIS managar allocates 
    // memory for two buffers of the size of the amplitudes
    std::shared_ptr<DIIS> diis (new DIIS(o*v));

    bool do_diis = true;//options_.get_bool("DIIS");

    // build integrals and fock matrix
    setup_integrals();

    double en = 0.0;

    outfile->Printf("\n");
    outfile->Printf(
      "    ==> Begin z iterations <==\n");
    outfile->Printf("\n");
    outfile->Printf(
      "        iter          energy       d(Energy)          |d(Z)|\n");

    do { 

        evaluate_residual_z(r2_);

        // update amplitudes
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                r2_[i * v + a] *= -0.5 / ( fv_[a] - fo_[i] );
            }
        }
        C_DCOPY(o * v, r2_, 1, t0_, 1);
        C_DAXPY(o * v, -1.0, z2_, 1, t0_, 1);  // error vector

        double nrm = C_DNRM2(o * v, t0_, 1);

        C_DCOPY(o * v, r2_, 1, z2_, 1); // new t2

        // The extrapolated parameters in DIIS for pp2rdm
        // are the amplitudes. This DIIS manager will write
        // the current amplitudes disk to avoid storing 
        // multiple copies in main memory.

        if ( do_diis ) {
            diis->WriteVector(z2_);
        }

        // The DIIS manager will write the current error vector to disk.
        if ( do_diis ) {
            diis->WriteErrorVector(t0_);
        }

        // The DIIS manager extrapolates the amplitudes, using
        // the amplitudes and error vectors generated in this
        // and previous iterations.
        if ( do_diis ) {
            diis->Extrapolate(z2_);
        }

        double dE = en;

        en = 0.0;
   
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                en += z2_[i * v + a] * v_iaia_[i * v + a];
            }  
        }  

        dE -= en;

        outfile->Printf("       %5i %15.10lf %15.10lf %15.10lf\n",iter,en,dE,nrm);

        iter++;

        if ( fabs(dE) < e_convergence_ && nrm < r_convergence_ ) break;

    }while(ci_iter_ < options_.get_int("MAXITER"));

    if (ci_iter_ == options_.get_int("MAXITER")){
        throw PsiException("  z iterations did not converge.",__FILE__,__LINE__);
    } 

    //printf("done var: %20.12lf proj: %20.12lf\n",evaluate_variational_energy(),evaluate_projection_energy());

    return en;

}

double pp2RDMSolver::pp2rdm_iterations(int & ci_iter_) {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    // Initialize our diis extrapolation manager.  Note that,
    // in this implemenation, the DIIS managar allocates 
    // memory for two buffers of the size of the amplitudes
    std::shared_ptr<DIIS> diis (new DIIS(o*v));

    bool do_diis = true;//options_.get_bool("DIIS");

    // build integrals and fock matrix
    setup_integrals();

    double en = 0.0;

    outfile->Printf("\n");
    outfile->Printf(
      "    ==> Begin p%s iterations <==\n",name_.c_str());
    outfile->Printf("\n");
    outfile->Printf(
      "        iter          energy       d(Energy)          |d(T)|\n");

    // restart from old solution?
    if ( restart_from_checkpoint_file_ ) {
        std::shared_ptr<PSIO> psio ( new PSIO() );
        psio->open(PSIF_V2RDM_CHECKPOINT,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_V2RDM_CHECKPOINT,"AMPLITUDES",(char*)t2_,o * v * sizeof(double));
        psio->close(PSIF_V2RDM_CHECKPOINT,1);
    }

    do { 

        evaluate_residual(r2_);

        // update amplitudes
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                r2_[i * v + a] *= -0.5 / ( fv_[a] - fo_[i] );
            }
        }
        C_DCOPY(o * v, r2_, 1, t0_, 1);
        C_DAXPY(o * v, -1.0, t2_, 1, t0_, 1);  // error vector

        double nrm = C_DNRM2(o * v, t0_, 1);

        C_DCOPY(o * v, r2_, 1, t2_, 1); // new t2

        // The extrapolated parameters in DIIS for pp2rdm
        // are the amplitudes. This DIIS manager will write
        // the current amplitudes disk to avoid storing 
        // multiple copies in main memory.

        if ( do_diis ) {
            diis->WriteVector(t2_);
        }

        // The DIIS manager will write the current error vector to disk.
        if ( do_diis ) {
            diis->WriteErrorVector(t0_);
        }

        // The DIIS manager extrapolates the amplitudes, using
        // the amplitudes and error vectors generated in this
        // and previous iterations.
        if ( do_diis ) {
            diis->Extrapolate(t2_);
        }

        double dE = en;

        en = evaluate_projection_energy();

        dE -= en;

        outfile->Printf("       %5i %15.10lf %15.10lf %15.10lf\n",ci_iter_,en,dE,nrm);

        ci_iter_++;

        if ( fabs(dE) < e_convergence_ && nrm < r_convergence_ ) break;

    }while(ci_iter_ < options_.get_int("MAXITER"));

    // restart from old solution?
    //if ( options_.get_str("P2RDM_TYPE") == "CCD" ) {
    //    restart_from_checkpoint_file_ = true;
    //}

    //if ( restart_from_checkpoint_file_ ) {
        std::shared_ptr<PSIO> psio ( new PSIO() );
        psio->open(PSIF_V2RDM_CHECKPOINT,PSIO_OPEN_NEW);
        psio->write_entry(PSIF_V2RDM_CHECKPOINT,"AMPLITUDES",(char*)t2_,o * v * sizeof(double));
        psio->close(PSIF_V2RDM_CHECKPOINT,1);
    //}

    if (ci_iter_ == options_.get_int("MAXITER")){
        throw PsiException(" pp2rdm  iterations did not converge.",__FILE__,__LINE__);
    } 

    //printf("done var: %20.12lf proj: %20.12lf\n",evaluate_variational_energy(),evaluate_projection_energy());

    return en;

}

// evaluate energy
double pp2RDMSolver::evaluate_projection_energy() {

    double en = 0.0;

    int o = nalpha_;
    int v = nmo_ - nalpha_;
    Normalization();
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            en += t2_[i * v + a] / t0_[i * v + a] * v_iaia_[i * v + a];
        }
    }

    return en;

}

// evaluate residual (z)
// see eq 19 of J. Chem. Phys. 141, 244104 (2014)
void pp2RDMSolver::evaluate_residual_z(double * residual) {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    memset((void*)residual,'\0',o * v * sizeof(double));

    // (ia|ia)
    C_DCOPY(o*v,v_iaia_,1,residual,1);


    double * VxT_v  = (double*)malloc(v * sizeof(double)); 
    double * VxT_o  = (double*)malloc(o * sizeof(double)); 
   
    // -2 ( (ja|ja) t(ja) + (ib|ib) t(ib) ) z(ia) 
    for (int a = 0; a < v; a++) {
        VxT_v[a] = -2.0 * C_DDOT(o, v_iaia_ + a, v, t2_ + a, v);
    }
    for (int i = 0; i < o; i++) {
        VxT_o[i] = -2.0 * C_DDOT(v, v_iaia_ + i*v, 1, t2_ + i*v, 1);
    }

    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {

            residual[i*v + a] += VxT_v[a] * z2_[i * v + a];
            residual[i*v + a] += VxT_o[i] * z2_[i * v + a];

        }
    }

    // -2 * (2 (ii|aa) - (ia|ai) - 2 (ia|ia) t(ia) ) z(ia)
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            residual[i * v + a] -= 2.0 * ( 2.0 * v_iiaa_[i * v + a] - v_iaia_[i * v + a] - v_iaia_[i * v + a] * t2_[i * v + a] ) * z2_[i * v + a];
        }
    }

    // r(ia) = (ab|ab) z(ib)
    F_DGEMM('n','n',v,o,v,1.0,v_abab_,v,z2_,v,1.0,residual,v);

    // r(ia) = z(ja) (ji|ji)
    F_DGEMM('n','n',v,o,o,1.0,z2_,v,v_ijij_,o,1.0,residual,v);

    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {

            double dum = 0.0;

            // -2 v(ia|ia) z(ja) t(ja) 
            for (int j = 0; j < o; j++) {
                dum += z2_[j * v + a] * t2_[j * v + a];
            }

            // -2 v(ia|ia) z(ib) t(ib) 
            for (int b = 0; b < v; b++) {
                dum += z2_[i * v + b] * t2_[i * v + b];
            }

            residual[i * v + a] -= 2.0 * v_iaia_[i * v + a] * dum;

            // t(jb) (ib|ib) z(ja) + t(jb) (ja|ja) z(ib)
            dum = 0.0;
            for (int j = 0; j < o; j++) {
                for (int b = 0; b < v; b++) {
                    dum += t2_[j * v + b] * v_iaia_[i * v + b] * z2_[j * v + a];
                    dum += t2_[j * v + b] * v_iaia_[j * v + a] * z2_[i * v + b];
                }
            }

            residual[i * v + a] += dum;

        }
    }

    free(VxT_v);
    free(VxT_o);

}
// evaluate residual
void pp2RDMSolver::evaluate_residual(double * residual) {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    memset((void*)residual,'\0',o * v * sizeof(double));

    // construct c0
    Normalization();

    evaluate_sigma(residual);

    if ( options_.get_str("P2RDM_TYPE") == "CCD" ) {

        double * VxT_v  = (double*)malloc(v * sizeof(double)); 
        double * VxT_o  = (double*)malloc(o * sizeof(double)); 
        double * VxT_oo = (double*)malloc(o * o * sizeof(double)); 
        
        for (int a = 0; a < v; a++) {
            VxT_v[a] = -2.0 * C_DDOT(o, v_iaia_ + a, v, t2_ + a, v);
        }
        for (int i = 0; i < o; i++) {
            VxT_o[i] = -2.0 * C_DDOT(v, v_iaia_ + i*v, 1, t2_ + i*v, 1);
        }

        // VxT_oo(i,j) = (jb|jb) t(i,b)
        F_DGEMM('t','n',o,o,v,1.0,v_iaia_,v,t2_,v,0.0,VxT_oo,o);

        // r2(i,a) += t(j,a) VxT_oo(i,j)
        F_DGEMM('n','n',v,o,o,1.0,t2_,v,VxT_oo,o,1.0,residual,v);

        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {

                double dum = 0.0;

                dum += VxT_v[a] * t2_[i * v + a];
                dum += VxT_o[i] * t2_[i * v + a];

                //for (int j = 0; j < o; j++) {
                //    dum -= 2.0 * v_iaia_[j * v + a] * t2_[j * v + a] * t2_[i * v + a];
                //}

                //for (int b = 0; b < v; b++) {
                //    dum -= 2.0 * v_iaia_[i * v + b] * t2_[i * v + b] * t2_[i * v + a];
                //}

                double t2 = t2_[i * v + a];
                dum += 2.0 * v_iaia_[i * v + a] * t2 * t2;

                //for (int j = 0; j < o; j++) {
                //    for (int b = 0; b < v; b++) {
                //        dum += v_iaia_[j * v + b] * t2_[j * v + a] * t2_[i * v + b];
                //    }
                //}

                residual[i * v + a] += dum;

            }
        }
        free(VxT_v);
        free(VxT_o);
        free(VxT_oo);
    }

    // last funny term for coupled-pair methods:
    if ( options_.get_str("P2RDM_TYPE") == "K" ) {
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                double dum = 0.0;
                for (int b = 0; b < v; b++) {
                    dum += t2_[i * v + b] / t0_[i * v + b] * v_iaia_[i * v + b];
                }
                for (int j = 0; j < o; j++) {
                    dum += t2_[j * v + a] / t0_[j * v + a] * v_iaia_[j * v + a];
                }
                dum -=  t2_[i * v + a] / t0_[i * v + a] * v_iaia_[i * v + a];

                residual[i * v + a] -= t2_[i * v + a] * dum;
            }
        }
    }else if ( options_.get_str("P2RDM_TYPE") == "CEPA(1)" ) {
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                double dum = 0.0;
                for (int b = 0; b < v; b++) {
                    dum += t2_[i * v + b] / t0_[i * v + b] * v_iaia_[i * v + b];
                }

                residual[i * v + a] -= t2_[i * v + a] * dum;
            }
        }
    }else if ( options_.get_str("P2RDM_TYPE") == "CID" ) {

        double dum = 0.0;
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                dum += t2_[i * v + a] / t0_[i * v + a] * v_iaia_[i * v + a];
            }
        }

        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                residual[i * v + a] -= t2_[i * v + a] * dum;
            }
        }

    }else if ( options_.get_str("P2RDM_TYPE") == "CEPA(0)" ) {

        // nothing to add

    }else if ( options_.get_str("P2RDM_TYPE") == "ACPF" ) {

        // nothing to add

    }else if ( options_.get_str("P2RDM_TYPE") == "AQCC" ) {

        // nothing to add

    }

}

void pp2RDMSolver::evaluate_sigma(double * sigma) {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    memset((void*)sigma,'\0',o*v*sizeof(double));

    // c0 * (ia|ia)
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            sigma[i * v + a] = v_iaia_[i * v + a] * t0_[i * v + a];
        }
    }

    // -2 * (2 (ii|aa) - (ia|ai) c(ai)
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            sigma[i * v + a] -= 2.0 * ( 2.0 * v_iiaa_[i * v + a] - v_iaia_[i * v + a] ) * t2_[i * v + a];
        }
    }

    // r(ia) = (ab|ab) c(ib)
    F_DGEMM('n','n',v,o,v,1.0,v_abab_,v,t2_,v,1.0,sigma,v);

    // r(ia) = c(ja) (ji|ji)
    F_DGEMM('n','n',v,o,o,1.0,t2_,v,v_ijij_,o,1.0,sigma,v);

    //for (int i = 0; i < o; i++) {
    //    for (int a = 0; a < v; a++) {
    //        double dum = 0.0;
    //        for (int b = 0; b < v; b++) {
    //            dum += v_abab_[a * v + b] * t2_[i * v + b];
    //        }
    //        for (int j = 0; j < o; j++) {
    //            dum += v_ijij_[i * o + j] * t2_[j * v + a];
    //        }
    //        sigma[i * v + a] += dum;
    //    }
    //}

}

// see eqs A8a - A11 of J. Chem. Phys. 141, 244104 (2014)
double pp2RDMSolver::BuildSeniorityZeroRDMs_pCCD(bool print){

    memset((void*)d1_,'\0', oei_dim_ *  sizeof(double));
    memset((void*)d2ab_,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)d2aa_,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    double energy = 0.0;

    // build opdm first:

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    double * xij = (double*)malloc(o*o*sizeof(double));
    double * xab = (double*)malloc(v*v*sizeof(double));

    memset((void*)xij,'\0',o*o*sizeof(double));
    memset((void*)xab,'\0',v*v*sizeof(double));

    // xji = [ia] * z[ja] (A9a)
    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            double dum = 0.0;
            for (int a = 0; a < v; a++) {
                dum += t2_[i * v + a] * z2_[j * v + a];
            }
            xij[j*o+i] = dum;
        }
    }

    // xba = t[ib] * z[ia] (A9b)
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            double dum = 0.0;
            for (int i = 0; i < o; i++) {
                dum += t2_[i * v + b] * z2_[i * v + a];
            }
            xab[b*v+a] = dum;
        }
    }

    // gamma_ij = 2 ( 1 - xij ) dij (A8a)
    for (int i = 0; i < o; i++) {
        d1_[INDEX(i,i)] = 2.0 * (1.0 - xij[i*o+i]);
    }
    // gamma_ab = 2 xab dab (A8b)
    for (int a = 0; a < v; a++) {
        d1_[INDEX(a+o,a+o)] = 2.0 * xab[a*v+a];
    }

    // TPDM:

    // xia = tib tja zjb (A11)
    double * xia = (double*)malloc(o*v*sizeof(double));
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            double dum = 0.0;
            for (int j = 0; j < o; j++) {
                for (int b = 0; b < v; b++) {
                    dum += t2_[i * v + b] * t2_[j * v + a] * z2_[j * v + b];
                }
            }
            xia[i * v + a] = dum;
        }
    }

    // note the expressions from J. Chem. Phys. 141, 244104 (2014) are spin-free

    // these four blocks only include ab contributions, so just divide by two

    // gamma(jj,ii) = 2 [ xji + dij ( 1 - 2xii) ]  (A10a)
    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            d2ab_[j*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+i] = xij[j * o + i] + (i==j) * (1.0 - 2.0 * xij[i * o + i]);
        }
    }
    // gamma(aa,ii) = 2 [tia + xia - 2tia (xaa + xii - tia zia)]
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            d2ab_[(a+o)*nmo_*nmo_*nmo_+(a+o)*nmo_*nmo_+i*nmo_+i] = t2_[i * v + a] + xia[i * v + a] 
                - 2.0 * t2_[i * v + a] * ( xab[a * v + a] + xij[i * o + i] - t2_[i * v + a] * z2_[i * v + a]);
        }
    }
    // gamma(ii,aa) = 2 zia
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            d2ab_[i*nmo_*nmo_*nmo_+i*nmo_*nmo_+(a+o)*nmo_+(a+o)] = z2_[i * v + a];
        }
    }
    // gamma(bb,aa) = 2 xba
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            d2ab_[(b+o)*nmo_*nmo_*nmo_+(b+o)*nmo_*nmo_+(a+o)*nmo_+(a+o)] = xab[b * v + a];
        }
    }

    // other blocks 

    // gamma(ab,ab) = 2 dab xaa ... actually already covered in gamma(aa,bb) block above
    //for (int a = 0; a < v; a++) {
    //    for (int b = 0; b < v; b++) {
    //        d2ab_[(a+o)*nmo_*nmo_*nmo_+(b+o)*nmo_*nmo_+(a+o)*nmo_+(b+o)] = (a==b) * xab[a * v + a];
    //    }
    //}

    // gamma(ij,ij) = 4.0 ( 1.0 - xii - xjj) + 2 * dij ( 3 xii - 1)
    // the case where i=j is covered in gamma(ii,jj), which leaves us
    // gamma(ij,ij) = 4.0 ( 1.0 - xii - xjj), for i!=j.  
    // i think for DOCI, d2ab(ij,ij) = d2aa(ij,ij), which means we can just drop the factor of four
    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            if ( i == j ) continue;
            d2ab_[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j] = 1.0 - xij[i * o + i] - xij[j * o + j];
        }
    }

    // gamma(ia,ia) = gamma(ai,ai) = 4 (xaa - tia zia)
    // i think for DOCI, d2ab(ia,ia) = d2aa(ia,ia), which means we can just drop the factor of four
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            double dum = xab[a * v + a] - t2_[i * v + a] * z2_[i * v + a];
            d2ab_[i*nmo_*nmo_*nmo_+(a+o)*nmo_*nmo_+i*nmo_+(a+o)] = dum;
            d2ab_[(a+o)*nmo_*nmo_*nmo_+i*nmo_*nmo_+(a+o)*nmo_+i] = dum;
        }
    }
    
    // symmetrize d2ab
    for (int p = 0; p < nmo_; p++) {
        for (int q = 0; q < nmo_; q++) {
            int pq = p * nmo_ + q;
            for (int r = 0; r < nmo_; r++) {
                for (int s = 0; s < nmo_; s++) {
                    int rs = r * nmo_ + s;
                    if ( pq > rs ) continue;
                    double dum = d2ab_[pq * nmo_ * nmo_ + rs] + d2ab_[rs * nmo_ * nmo_ + pq];
                    d2ab_[pq * nmo_ * nmo_ + rs] = 0.5 * dum;
                    d2ab_[rs * nmo_ * nmo_ + pq] = 0.5 * dum;
                }
            }
        }
    }

    // construct d2aa from d2ab:
    for (int p = 0; p < nmo_; p++) {
        for (int q = 0; q < nmo_; q++) {
            for (int r = 0; r < nmo_; r++) {
                for (int s = 0; s < nmo_; s++) {
                    d2aa_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s] = d2ab_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s]
                                                                 - d2ab_[q*nmo_*nmo_*nmo_+p*nmo_*nmo_+r*nmo_+s];

                }
            }
        }
    }

    // check trace
    double tra = 0.0;
    double trb = 0.0;
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            tra += d2aa_[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j];
            trb += d2ab_[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j];
        }
    }


    // check energy

    double e1 = 0.0;
    for (int i = 0; i < nmo_; i++) {
        e1 += d1_[INDEX(i,i)] * oei_[INDEX(i,i)];
    }

    // evaluating two-electron part of energy should be simple, since pair p2RDM looks like DOCI

    double e2 = 0.0;

    // D2(ij,ij) (ii|jj)
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            if ( i == j ) continue;
            int ii = INDEX(i,i);
            int jj = INDEX(j,j);
            double eint = C_DDOT(nQ_,Qmo_ + ii,nmo_*(nmo_+1)/2,Qmo_+jj,nmo_*(nmo_+1)/2);
            e2 += d2ab_[i * nmo_*nmo_*nmo_ + j * nmo_*nmo_ + i * nmo_ + j] * eint;
            e2 += d2aa_[i * nmo_*nmo_*nmo_ + j * nmo_*nmo_ + i * nmo_ + j] * eint;
        }
    }
    // D2(ij,ji) (ij|ji)
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            if ( i == j ) continue;
            int ij = INDEX(i,j);
            double eint = C_DDOT(nQ_,Qmo_ + ij,nmo_*(nmo_+1)/2,Qmo_+ij,nmo_*(nmo_+1)/2);
            e2 += d2aa_[i * nmo_*nmo_*nmo_ + j * nmo_*nmo_ + j * nmo_ + i] * eint;
        }
    }
    // D2(ii,jj) (ij|ij)
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            int ij = INDEX(i,j);
            double eint = C_DDOT(nQ_,Qmo_ + ij,nmo_*(nmo_+1)/2,Qmo_+ij,nmo_*(nmo_+1)/2);
            e2 += d2ab_[i * nmo_*nmo_*nmo_ + i * nmo_*nmo_ + j * nmo_ + j] * eint;
        }
    }

    double en = e1 + e2;
    if ( print ) {
        outfile->Printf("\n");
        outfile->Printf("        pCCD one-electron energy = %20.12lf\n",e1);
        outfile->Printf("        pCCD two-electron energy = %20.12lf\n",e2);
        outfile->Printf("        * pCCD total energy      = %20.12lf\n",e1 + e2 + enuc_);
    }

    free(xij);
    free(xab);
    free(xia);

    return e1 + e2;

}

double pp2RDMSolver::BuildSeniorityZeroRDMs(bool print){

    memset((void*)d1_,'\0', oei_dim_ *  sizeof(double));
    memset((void*)d2ab_,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)d2aa_,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    double energy = 0.0;

    // build opdm first:

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    double * Dij = (double*)malloc(o*o*sizeof(double));
    double * Dab = (double*)malloc(v*v*sizeof(double));

    memset((void*)Dij,'\0',o*o*sizeof(double));
    memset((void*)Dab,'\0',v*v*sizeof(double));

    // construct c0
    Normalization_RDMs();

    // Dij = - tb[abjk] * tb[abik]
    // Dij = - tb[aajk] * tb[aaik]
    // Dij = - tb[aajj] * tb[aaii] dij
    for (int i = 0; i < o; i++) {
        double dum = 0.0;
        for (int a = 0; a < v; a++) {
            dum -= t2_[i * v + a] * t2_[i * v + a];
        }
        Dij[i*o+i] = dum;
    }

    // Dab = tb[acij] * tb[bcij]
    // Dab = tb[acii] * tb[bcii]
    // Dab = tb[aaii] * tb[aaii] dab
    for (int a = 0; a < v; a++) {
        double dum = 0.0;
        for (int i = 0; i < o; i++) {
            dum += t2_[i * v + a] * t2_[i * v + a];
        }
        Dab[a*v+a] = dum;
    }

    // TPDM:

    // D2(ijkl): dik * djl + Dij[ik] * djl + Dij[jl] * dik + tb[abij] * tb[abkl]
    //
    //          notes: 
    //
    //                1. Dij is diagonal
    //                2. tb[abij] = tb[aaii] dab dij
    //
    // D2(ijij): 1 + Dij[ii] + Dij[jj] + tb[aaij] tb[aaij] dij
    // D2(iijj): tb[aaii] tb[aajj] 
    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            double dum = 
            d2ab_[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j] = 1.0 + Dij[i*o+i] + Dij[j*o+j];
        }
        double dum = 0.0;
        for (int a = 0; a < v; a++) {
            dum += t2_[i * v + a] * t2_[i * v + a];
        }
        d2ab_[i*nmo_*nmo_*nmo_+i*nmo_*nmo_+i*nmo_+i] += dum;
    }
    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            if ( i == j ) continue;
            double dum = 0.0;
            for (int a = 0; a < v; a++) {
                dum += t2_[i * v + a] * t2_[j * v + a];
            }
            d2ab_[i*nmo_*nmo_*nmo_+i*nmo_*nmo_+j*nmo_+j] += dum;
        }
    }

    // ijab / abij = t0(ijab) * t2(ijab)
    //             = t0(iiaa) * t2(iiaa) dij dab
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            double dum = t0_[i * v + a] * t2_[i * v + a];
            d2ab_[i*nmo_*nmo_*nmo_+i*nmo_*nmo_+(a+o)*nmo_+(a+o)] = dum;
            d2ab_[(a+o)*nmo_*nmo_*nmo_+(a+o)*nmo_*nmo_+i*nmo_+i] = dum;
        }
    }

    // iajb / aibj: dij Dab(ab) - tb(ackj) * tb(bcki)
    //              dij Dab(aa) dab - tb(aaii) * tb(aaii) dij dab
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            double dum = Dab[a*v+a] - t2_[i * v + a] * t2_[i * v + a];
            d2ab_[i*nmo_*nmo_*nmo_+(a+o)*nmo_*nmo_+i*nmo_+(a+o)] = dum;
            d2ab_[(a+o)*nmo_*nmo_*nmo_+i*nmo_*nmo_+(a+o)*nmo_+i] = dum;
        }
    }


    // abcd: tb[ijab] * tb[ijcd]
    //       tb[iiaa] * tb[iicc] dab dcd
    for (int a = 0; a < v; a++) {
        for (int c = 0; c < v; c++) {
            double dum = 0.0;
            for (int i = 0; i < o; i++) {
                dum += t2_[i * v + a] * t2_[i * v + c];
            }
            d2ab_[(a+o)*nmo_*nmo_*nmo_+(a+o)*nmo_*nmo_+(c+o)*nmo_+(c+o)] = dum;
        }
    }

    // construct d2aa from d2ab:
    for (int p = 0; p < nmo_; p++) {
        for (int q = 0; q < nmo_; q++) {
            for (int r = 0; r < nmo_; r++) {
                for (int s = 0; s < nmo_; s++) {
                    d2aa_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s] = d2ab_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s]
                                                                 - d2ab_[q*nmo_*nmo_*nmo_+p*nmo_*nmo_+r*nmo_+s];

                }
            }
        }
    }

    // check trace
    double tra = 0.0;
    double trb = 0.0;
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            tra += d2aa_[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j];
            trb += d2ab_[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j];
        }
    }

    // build d1 by contraction
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {

            double dum = 0.0;
            for (int p = 0; p < nmo_; p++) {
                dum += d2aa_[i*nmo_*nmo_*nmo_+p*nmo_*nmo_+j*nmo_+p];
                dum += d2ab_[i*nmo_*nmo_*nmo_+p*nmo_*nmo_+j*nmo_+p];
            }
            d1_[INDEX(i,j)] = 2.0 * dum / (2.0 * nalpha_-1.0);

        }
    }

    // check energy

    double e1 = 0.0;
    for (int i = 0; i < nmo_; i++) {
        e1 += d1_[INDEX(i,i)] * oei_[INDEX(i,i)];
    }

    // evaluating two-electron part of energy should be simple, since pair p2RDM looks like DOCI

    double e2 = 0.0;

    // D2(ij,ij) (ii|jj)
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            if ( i == j ) continue;
            int ii = INDEX(i,i);
            int jj = INDEX(j,j);
            double eint = C_DDOT(nQ_,Qmo_ + ii,nmo_*(nmo_+1)/2,Qmo_+jj,nmo_*(nmo_+1)/2);
            e2 += d2ab_[i * nmo_*nmo_*nmo_ + j * nmo_*nmo_ + i * nmo_ + j] * eint;
            e2 += d2aa_[i * nmo_*nmo_*nmo_ + j * nmo_*nmo_ + i * nmo_ + j] * eint;
        }
    }
    // D2(ij,ji) (ij|ji)
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            if ( i == j ) continue;
            int ij = INDEX(i,j);
            double eint = C_DDOT(nQ_,Qmo_ + ij,nmo_*(nmo_+1)/2,Qmo_+ij,nmo_*(nmo_+1)/2);
            e2 += d2aa_[i * nmo_*nmo_*nmo_ + j * nmo_*nmo_ + j * nmo_ + i] * eint;
        }
    }
    // D2(ii,jj) (ij|ij)
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            int ij = INDEX(i,j);
            double eint = C_DDOT(nQ_,Qmo_ + ij,nmo_*(nmo_+1)/2,Qmo_+ij,nmo_*(nmo_+1)/2);
            e2 += d2ab_[i * nmo_*nmo_*nmo_ + i * nmo_*nmo_ + j * nmo_ + j] * eint;
        }
    }

    double en = e1 + e2;
    if ( print ) {
        outfile->Printf("\n");
        outfile->Printf("        p2RDM one-electron energy = %20.12lf\n",e1);
        outfile->Printf("        p2RDM two-electron energy = %20.12lf\n",e2);
        outfile->Printf("        * p2RDM total energy      = %20.12lf\n",e1 + e2 + enuc_);
    }

    free(Dij);
    free(Dab);

    // need to re-scale amplitudes for ACPF and AQCC
    if ( options_.get_str("P2RDM_TYPE") == "ACPF" || options_.get_str("P2RDM_TYPE") == "AQCC" ) {

        double fac = 0.0;
        if ( options_.get_str("P2RDM_TYPE") == "ACPF" ) {
            fac = 1.0 / o;
        }else if ( options_.get_str("P2RDM_TYPE") == "AQCC" ) {
            fac = 1.0 - (2.0 * o - 2.0) * (2.0 * o - 3.0) / (2.0 * o * (2.0 * o - 1.0));
        }

        double dum = 1.0 + fac * C_DDOT(o * v, t2_, 1, t2_, 1);
        double c0 = sqrt(1.0 / dum);
        C_DSCAL(o * v, 1.0 / c0, t2_, 1);
    }

    return e1 + e2;

}

void pp2RDMSolver::Normalization() {

    int o = nalpha_;
    int v = nmo_ - nalpha_;
    
    if ( options_.get_str("P2RDM_TYPE") == "K" ) {

        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                double dum = 0.0;
                for (int b = 0; b < v; b++) {
                    double t2 = t2_[i * v + b];
                    dum += t2 * t2;
                }
                for (int j = 0; j < o; j++) {
                    double t2 = t2_[j * v + a];
                    dum += t2 * t2;
                }
                double t2 = t2_[i * v + a];
                dum -= t2 * t2;

                t0_[i * v + a] = sqrt(1.0 - dum);
            }
        }

    }else if ( options_.get_str("P2RDM_TYPE") == "CEPA(1)" ) {

        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                double dum = 0.0;
                for (int b = 0; b < v; b++) {
                    double t2 = t2_[i * v + b];
                    dum += t2 * t2;
                }

                t0_[i * v + a] = sqrt(1.0 - dum);
            }
        }

    }else if ( options_.get_str("P2RDM_TYPE") == "CID" ) {

        double dum = C_DDOT(o * v, t2_, 1, t2_, 1);
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                t0_[i * v + a] = sqrt(1.0 - dum);
            }
        }

    }else if ( options_.get_str("P2RDM_TYPE") == "CEPA(0)" ) {

        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                t0_[i * v + a] = 1.0;
            }
        }

    }else if ( options_.get_str("P2RDM_TYPE") == "ACPF" ) {

        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                t0_[i * v + a] = 1.0;
            }
        }

    }else if ( options_.get_str("P2RDM_TYPE") == "AQCC" ) {

        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                t0_[i * v + a] = 1.0;
            }
        }

    }else if ( options_.get_str("P2RDM_TYPE") == "CCD" ) {

        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                t0_[i * v + a] = 1.0;
            }
        }

    }
}

void pp2RDMSolver::Normalization_RDMs() {

    int o = nalpha_;
    int v = nmo_ - nalpha_;
    
    if ( options_.get_str("P2RDM_TYPE") == "K" ) {

        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                double dum = 0.0;
                for (int b = 0; b < v; b++) {
                    double t2 = t2_[i * v + b];
                    dum += t2 * t2;
                }
                for (int j = 0; j < o; j++) {
                    double t2 = t2_[j * v + a];
                    dum += t2 * t2;
                }
                double t2 = t2_[i * v + a];
                dum -= t2 * t2;

                t0_[i * v + a] = sqrt(1.0 - dum);
            }
        }

    }else if ( options_.get_str("P2RDM_TYPE") == "CEPA(1)" ) {

        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                double dum = 0.0;
                for (int b = 0; b < v; b++) {
                    double t2 = t2_[i * v + b];
                    dum += t2 * t2;
                }

                t0_[i * v + a] = sqrt(1.0 - dum);
            }
        }

    }else if ( options_.get_str("P2RDM_TYPE") == "CID" ) {

        double dum = C_DDOT(o * v, t2_, 1, t2_, 1);
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                t0_[i * v + a] = sqrt(1.0 - dum);
            }
        }

    }else if ( options_.get_str("P2RDM_TYPE") == "CEPA(0)" ) {

        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                t0_[i * v + a] = 1.0;
            }
        }

    }else if ( options_.get_str("P2RDM_TYPE") == "ACPF" || options_.get_str("P2RDM_TYPE") == "AQCC" ) {

        double fac = 0.0;
        if ( options_.get_str("P2RDM_TYPE") == "ACPF" ) {
            fac = 1.0 / o;
        }else if ( options_.get_str("P2RDM_TYPE") == "AQCC" ) {
            fac = 1.0 - (2.0 * o - 2.0) * (2.0 * o - 3.0) / (2.0 * o * (2.0 * o - 1.0));
        }

        double dum = 1.0 + fac * C_DDOT(o * v, t2_, 1, t2_, 1);
        double c0 = sqrt(1.0 / dum);
        C_DSCAL(o * v, c0, t2_, 1);
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                t0_[i * v + a] = c0;
            }
        }

    }
}

void pp2RDMSolver::setup_integrals() {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    // build required integral
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            v_iiaa_[i * v + a] = C_DDOT(nQ_,Qmo_ + INDEX(i,i),nmo_*(nmo_+1)/2,Qmo_+INDEX(a+o,a+o),nmo_*(nmo_+1)/2);
            v_iaia_[i * v + a] = C_DDOT(nQ_,Qmo_ + INDEX(i,a+o),nmo_*(nmo_+1)/2,Qmo_+INDEX(i,a+o),nmo_*(nmo_+1)/2);
        }
    }
    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            v_ijij_[i * o + j] = C_DDOT(nQ_,Qmo_ + INDEX(i,j),nmo_*(nmo_+1)/2,Qmo_+INDEX(i,j),nmo_*(nmo_+1)/2);
        }
    }
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            v_abab_[a * v + b] = C_DDOT(nQ_,Qmo_ + INDEX(a+o,b+o),nmo_*(nmo_+1)/2,Qmo_+INDEX(a+o,b+o),nmo_*(nmo_+1)/2);
        }
    }

    // fock matrix:

    for (int i = 0; i < o; i++) {
        double dum = oei_[INDEX(i,i)];
        for (int k = 0; k < o; k++) {
            dum += 2.0 * C_DDOT(nQ_,Qmo_ + INDEX(i,i),nmo_*(nmo_+1)/2,Qmo_+INDEX(k,k),nmo_*(nmo_+1)/2);
            dum -=       C_DDOT(nQ_,Qmo_ + INDEX(i,k),nmo_*(nmo_+1)/2,Qmo_+INDEX(k,i),nmo_*(nmo_+1)/2);
        }
        fo_[i] = dum;
    }
    for (int a = 0; a < v; a++) {
        double dum = oei_[INDEX(a+o,a+o)];
        for (int k = 0; k < o; k++) {
            dum += 2.0 * C_DDOT(nQ_,Qmo_ + INDEX(a+o,a+o),nmo_*(nmo_+1)/2,Qmo_+INDEX(k,k),nmo_*(nmo_+1)/2);
            dum -=       C_DDOT(nQ_,Qmo_ + INDEX(a+o,k),nmo_*(nmo_+1)/2,Qmo_+INDEX(k,a+o),nmo_*(nmo_+1)/2);
        }
        fv_[a] = dum;
    }

    escf_ = enuc_;
    for (int i = 0; i < o; i++) {
        escf_ += oei_[INDEX(i,i)] + fo_[i];
    }

    //outfile->Printf("hey new scf energy: %20.12lf\n",escf_);
}

}} //end namespaces
