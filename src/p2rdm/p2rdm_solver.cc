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

#include "p2rdm_solver.h"

#include <focas/focas_c_interface.h>
#include <misc/threeindexintegrals.h>
#include <misc/omp.h>
#include <misc/blas.h>
#include <misc/diis.h>

using namespace psi;
using namespace fnocc;

namespace hilbert{

p2RDMSolver::p2RDMSolver(SharedWavefunction reference_wavefunction,Options & options):
    Wavefunction(options){
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

p2RDMSolver::~p2RDMSolver()

{
    free(d1_);
    free(d2ab_);
    free(d2aa_);
    free(t2_);
    free(t0_);
    free(r2_);
    free(orbital_lagrangian_);
    free(Qoo_);
    free(Qov_);
    free(Qvo_);
    free(Qvv_);
}

void  p2RDMSolver::common_init(){

    if ( options_.get_str("P2RDM_TYPE") != "K" && options_.get_str("P2RDM_TYPE") != "CEPA(0)" && options_.get_str("P2RDM_TYPE") != "CID") {
        throw PsiException("invalid p2rdm type",__FILE__,__LINE__);
    }

    is_df_ = false;

    if ( options_.get_str("SCF_TYPE") == "DF" || options_.get_str("SCF_TYPE") == "CD" ) {
        is_df_ = true;
    }

    if ( !is_df_ ) {
        throw PsiException("plugin p2rdm only works with scf_type = df or cd, for now",__FILE__,__LINE__);
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
        throw PsiException("plugin p2rdm works only with c1 symmetry",__FILE__,__LINE__);
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

    // p2rdm convergence thresholds:
    r_convergence_  = options_.get_double("R_CONVERGENCE");
    e_convergence_  = options_.get_double("E_CONVERGENCE");

    maxiter_        = options_.get_int("MAXITER");

    // memory check happens here

    outfile->Printf("\n\n");
    outfile->Printf( "        ***************************************************\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *    p2RDM                                        *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *    parametric 2RDM                              *\n");
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
    outfile->Printf("        exact diagonal Hessian:             %5s\n",options_.get_bool("ORBOPT_EXACT_DIAGONAL_HESSIAN") ? "true" : "false");
    outfile->Printf("        number of DIIS vectors:             %5i\n",options_.get_int("ORBOPT_NUM_DIIS_VECTORS"));
    outfile->Printf("        print iteration info:               %5s\n",options_.get_bool("ORBOPT_WRITE") ? "true" : "false");

    // if using 3-index integrals, transform them before allocating any memory integrals, transform
    if ( is_df_ ) {
        outfile->Printf("\n");
        outfile->Printf("    ==> Transform three-electron integrals <==\n");
        outfile->Printf("\n");

        double start = omp_get_wtime();

        ThreeIndexIntegrals(reference_wavefunction_,nQ_,memory_);

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

    // amplitudes
    r2_ = (double*)malloc(o * o * v * v * sizeof(double)); // residual
    t2_ = (double*)malloc(o * o * v * v * sizeof(double)); // amplitudes
    t0_ = (double*)malloc(o * o * v * v * sizeof(double)); // normalization coefficients

    memset((void*)r2_,'\0',o * o * v * v * sizeof(double));
    memset((void*)t2_,'\0',o * o * v * v * sizeof(double));
    memset((void*)t0_,'\0',o * o * v * v * sizeof(double));

    // fock matrix
    foo_ = (double*)malloc(o*o*sizeof(double));
    fvv_ = (double*)malloc(v*v*sizeof(double));

    memset((void*)foo_,'\0',o * o * sizeof(double));
    memset((void*)fvv_,'\0',v * v * sizeof(double));

    // three-index integrals
    Qoo_ = (double*)malloc(o*o*nQ_*sizeof(double));
    Qov_ = (double*)malloc(o*v*nQ_*sizeof(double));
    Qvo_ = (double*)malloc(v*o*nQ_*sizeof(double));
    Qvv_ = (double*)malloc(v*v*nQ_*sizeof(double));

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

    orbital_lagrangian_ = (double*)malloc(nmo_*nmo_*sizeof(double));

    d1_ = (double*)malloc( nmo_ * (nmo_+1)/2 * sizeof(double));
    d2ab_ = (double*)malloc(nmo_ * nmo_ * nmo_ * nmo_ * sizeof(double));
    d2aa_ = (double*)malloc(nmo_ * nmo_ * nmo_ * nmo_ * sizeof(double));

}

// compute the energy!
double p2RDMSolver::compute_energy() {

    int iter = 0;

    //outfile->Printf("\n");
    outfile->Printf("    ==> Begin oo-%s iterations <==\n",name_.c_str());
    outfile->Printf("\n");

    double dE = 0;
    double energy = -enuc_;
    do {

        ci_iter_ = 0;

        double ci_energy = p2rdm_iterations(ci_iter_);

        double rdm_energy = build_rdms(true);

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
            outfile->Printf("    ==> oo-%s summary <==\n",name_.c_str());
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
        outfile->Printf("      oo-%s iterations converged!\n",name_.c_str());
        outfile->Printf("\n");

        outfile->Printf("    * oo-%s total energy:              %20.12lf\n",name_.c_str(),energy+enuc_);
        outfile->Printf("\n");
    }else {
        outfile->Printf("\n");
        outfile->Printf("      %s iterations converged!\n",name_.c_str());
        outfile->Printf("\n");

        outfile->Printf("    * %s total energy:              %20.12lf\n",name_.c_str(),energy+enuc_);
        outfile->Printf("\n");
    }

    // Natural Orbitals
    build_rdms(false);
    std::shared_ptr<Matrix> opdm (new Matrix(nmo_,nmo_));
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            opdm->pointer()[i][j] = d1_[INDEX(i,j)];
        }
    }
    std::shared_ptr<Vector> eigval (new Vector("Natural Orbital Occupation Numbers (spin free)",nmo_));
    std::shared_ptr<Matrix> eigvec (new Matrix(nmo_,nmo_));
    opdm->diagonalize(eigvec,eigval,descending);
    eigval->print();

    Process::environment.globals["CURRENT ENERGY"]    = energy + enuc_;
    Process::environment.globals["p2RDM TOTAL ENERGY"] = energy + enuc_;

    // push final transformation matrix onto Ca_ and Cb_
    if ( options_.get_bool("OPTIMIZE_ORBITALS") ) {
        UpdateTransformationMatrix(reference_wavefunction_,newMO_,Ca_,Cb_,orbopt_transformation_matrix_);
    }

    // write tpdm to disk?
    if ( options_.get_bool("TPDM_WRITE_FULL") ) {
        WriteTPDM();
    }

    return energy + enuc_;

}

double p2RDMSolver::update_amplitudes(bool do_diis, std::shared_ptr<DIIS> diis) {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    // df (ai|bj)

    double * integrals = (double*)malloc(o*o*v*v*sizeof(double));
    F_DGEMM('n', 't', o * v, o * v, nQ_, 1.0, Qvo_, o * v, Qvo_, o * v, 0.0, integrals, o * v);

    double * dt = (double*)malloc(o*o*v*v*sizeof(double));

    #pragma omp parallel for schedule(static)
    for (long int a = 0; a < v; a++) {

        double da = fvv_[a*v+a];

        for (long int b = 0; b < v; b++) {

            double dab = da + fvv_[b*v+b];

            for (long int i = 0; i < o; i++) {

                double dabi = dab - foo_[i*o+i];

                for (long int j = 0; j < o; j++) {

                    long int aibj = a * v * o * o + i * v * o + b * o + j;
                    long int abij = a * v * o * o + b * o * o + i * o + j;

                    double dabij = dabi - foo_[j*o+j];

                    dt[abij] = -(integrals[aibj] + r2_[abij]) / dabij;
                }
            }
        }
    }

    C_DAXPY(o*o*v*v, 1.0, dt, 1, t2_, 1);

    if ( do_diis ) {
        diis->WriteVector(t2_);
        diis->WriteErrorVector(dt);
        diis->Extrapolate(t2_);
    }

    free(integrals);
    free(dt);

    return C_DNRM2(o*o*v*v,dt,1);
}

double p2RDMSolver::evaluate_projection_energy() {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    Normalization();

    // df (ai|bj)

    double * integrals = (double*)malloc(o*o*v*v*sizeof(double));
    F_DGEMM('n', 't', o * v, o * v, nQ_, 1.0, Qov_, o * v, Qov_, o * v, 0.0, integrals, o * v);

    double energy = 0.0;

    for (long int a = 0; a < v; a++) {
        for (long int b = 0; b < v; b++) {
            for (long int i = 0; i < o; i++) {
                for (long int j = 0; j < o; j++) {
                    long int abij = a * v * o * o + b * o * o + i * o + j;
                    long int iajb = i * v * v * o + a * v * o + j * v + b;
                    long int jaib = j * v * v * o + a * v * o + i * v + b;
                    energy += (2.0 * integrals[iajb] - integrals[jaib]) * t2_[abij] / t0_[abij];
                }
            }
        }
    }

    return energy;
}

void p2RDMSolver::PackSpatialDensity() {

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

double p2RDMSolver::RotateOrbitals() {

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

    OrbOpt(orbopt_transformation_matrix_,

          oei_, oei_dim_,

          tei_, tei_dim_,

          d1_, oei_dim_, d2_, d2_dim_,

          symmetry_energy_order, nrstc, nmo_, nrstv, nirrep_,

          orbopt_data_, orbopt_outfile_, orbital_lagrangian_);

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

    // evaluate energy


    double e1 = 0.0;
    double e2 = 0.0;
    evaluate_rdm_energy(e1,e2);

    if ( orbopt_data_[8] > 0 ) {
        outfile->Printf("        p2RDM one-electron energy = %20.12lf\n",e1);
        outfile->Printf("        p2RDM two-electron energy = %20.12lf\n",e2);
        outfile->Printf("        * p2RDM total energy      = %20.12lf\n",e1 + e2 + enuc_);
    }

    return e1 + e2;

}

double p2RDMSolver::p2rdm_iterations(int & ci_iter_) {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    std::shared_ptr<DIIS> diis (new DIIS(o*o*v*v));

    bool do_diis = options_.get_bool("DIIS");

    // build integrals and fock matrix
    setup_integrals();

    double energy = 0.0;

    outfile->Printf("\n");
    outfile->Printf(
      "    ==> Begin %s iterations <==\n",name_.c_str());
    outfile->Printf("\n");
    outfile->Printf(
      "        iter          energy       d(Energy)          |d(T)|\n");

    do { 
        // evaluate residual
        evaluate_residual();

        // update amplitudes
        double nrm = update_amplitudes(do_diis, diis);

        double dE = energy;

        energy = evaluate_projection_energy();

        dE -= energy;

        outfile->Printf("       %5i %15.10lf %15.10lf %15.10lf\n",ci_iter_,energy,dE,nrm);

        ci_iter_++;

        if ( fabs(dE) < e_convergence_ && nrm < r_convergence_ ) break;

    }while(ci_iter_ < options_.get_int("MAXITER"));

    if (ci_iter_ == options_.get_int("MAXITER")){
        throw PsiException(" p2rdm  iterations did not converge.",__FILE__,__LINE__);
    } 

    return energy;

}


double p2RDMSolver::build_rdms(bool print){

    memset((void*)d1_,'\0', oei_dim_ *  sizeof(double));
    memset((void*)d2ab_,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)d2aa_,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    // build opdm first:
    int o = nalpha_;
    int v = nmo_ - nalpha_;

    double * Dij = (double*)malloc(o*o*sizeof(double));
    double * Dab = (double*)malloc(v*v*sizeof(double));

    memset((void*)Dij,'\0',o*o*sizeof(double));
    memset((void*)Dab,'\0',v*v*sizeof(double));

    double * ta = (double*)malloc(o*o*v*v*sizeof(double));
    double * tb = t2_;

    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            for (int a = 0; a < v; a++) {
                for (int b = 0; b < v; b++) {
                    ta[a*o*o*v + b*o*o + i*o + j] = tb[a*o*o*v + b*o*o + i*o + j] - tb[b*o*o*v + a*o*o + i*o + j];
                }
            }
        }
    }
    Normalization();

    // Dij
    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            double dum = 0.0;
            for (int k = 0; k < o; k++) {
                for (int a = 0; a < v; a++) {
                    for (int b = 0; b < v; b++) {
                        dum -= 0.5 * ta[a*o*o*v+b*o*o+j*o+k] * ta[a*o*o*v+b*o*o+i*o+k];
                        dum -= tb[a*o*o*v+b*o*o+j*o+k] * tb[a*o*o*v+b*o*o+i*o+k];
                    }
                }
            }
            Dij[i*o+j] = dum;
        }
    }
    // Dab
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            double dum = 0.0;
            for (int c = 0; c < v; c++) {
                for (int i = 0; i < o; i++) {
                    for (int j = 0; j < o; j++) {
                        dum += 0.5 * ta[a*o*o*v+c*o*o+i*o+j] * ta[b*o*o*v+c*o*o+i*o+j];
                        dum += tb[a*o*o*v+c*o*o+i*o+j] * tb[b*o*o*v+c*o*o+i*o+j];
                    }
                }
            }
            Dab[a*v+b] = dum;
        }
    }

    // TPDM:

    // ijkl:
    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            for (int k = 0; k < o; k++) {
                for (int l = 0; l < o; l++) {
                    double dumb = 0.0;

                    dumb += (i==k) * (j==l);
                    dumb += Dij[i*o+k] * (j==l) + Dij[j*o+l] * (i==k);

                    for (int a = 0; a < v; a++) {
                        for (int b = 0; b < v; b++) {
                            dumb += tb[a*o*o*v+b*o*o+i*o+j] * tb[a*o*o*v+b*o*o+k*o+l];
                        }
                    }
                    d2ab_[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l] = dumb;
                }
            }
        }
    }

    // ijab / abij
    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            for (int a = 0; a < v; a++) {
                for (int b = 0; b < v; b++) {
                    double dum = t0_[a*o*o*v+b*o*o+i*o+j] * tb[a*o*o*v+b*o*o+i*o+j];
                    d2ab_[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+(a+o)*nmo_+(b+o)] = dum;
                    d2ab_[(a+o)*nmo_*nmo_*nmo_+(b+o)*nmo_*nmo_+i*nmo_+j] = dum;
                }
            }
        }
    }

    // iajb / aibj:
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            for (int j = 0; j < o; j++) {
                for (int b = 0; b < v; b++) {

                    double dumb = (i==j) * Dab[a*v+b];
                    for (int k = 0; k < o; k++) {
                        for (int c = 0; c < v; c++) {
                            dumb -= tb[a*o*o*v+c*o*o+k*o+j] * tb[b*o*o*v+c*o*o+k*o+i];
                        }
                    }
                    d2ab_[i*nmo_*nmo_*nmo_+(a+o)*nmo_*nmo_+j*nmo_+(b+o)] = dumb;
                    d2ab_[(a+o)*nmo_*nmo_*nmo_+i*nmo_*nmo_+(b+o)*nmo_+j] = dumb;

                }
            }
        }
    }

    // iabj / aijb (this guy is a little weird for the ab block)
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            for (int j = 0; j < o; j++) {
                for (int b = 0; b < v; b++) {
                    double dumb = 0.0;
                    for (int k = 0; k < o; k++) {
                        for (int c = 0; c < v; c++) {
                            dumb += tb[a*o*o*v+c*o*o+j*o+k] * ta[b*o*o*v+c*o*o+i*o+k];
                            dumb += ta[a*o*o*v+c*o*o+j*o+k] * tb[b*o*o*v+c*o*o+i*o+k];
                        }
                    }
                    d2ab_[i*nmo_*nmo_*nmo_+(a+o)*nmo_*nmo_+(b+o)*nmo_+j] = dumb;
                    d2ab_[(a+o)*nmo_*nmo_*nmo_+i*nmo_*nmo_+j*nmo_+(b+o)] = dumb;
                }
            }
        }
    }

    F_DGEMM('t','n',v*v,v*v,o*o,1.0,tb,o*o,tb,o*o,0.0,d2aa_,v*v);
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            for (int c = 0; c < v; c++) {
                for (int d = 0; d < v; d++) {
                    d2ab_[(a+o)*nmo_*nmo_*nmo_+(b+o)*nmo_*nmo_+(c+o)*nmo_+(d+o)] = d2aa_[a*v*v*v+b*v*v+c*v+d];
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
    double e2 = 0.0;
    evaluate_rdm_energy(e1,e2);

    if ( print ) {
        outfile->Printf("\n");
        outfile->Printf("    p2RDM one-electron energy = %20.12lf\n",e1);
        outfile->Printf("    p2RDM two-electron energy = %20.12lf\n",e2);
        outfile->Printf("    * p2RDM total energy      = %20.12lf\n",e1 + e2 + enuc_);
        outfile->Printf("\n");
    }

    free(Dij);
    free(Dab);

    return e1 + e2;

}

void p2RDMSolver::evaluate_rdm_energy(double &e1, double &e2){

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    e1 = 0.0;
    e2 = 0.0;

    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            for (int k = 0; k < nmo_; k++) {
                for (int l = 0; l < nmo_; l++) {
                    double d2 = d2aa_[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l] + d2ab_[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];
                    int ik = INDEX(i,k);
                    int jl = INDEX(j,l);
                    double k2 = C_DDOT(nQ_,Qmo_ + ik,nmo_*(nmo_+1)/2,Qmo_+jl,nmo_*(nmo_+1)/2);
                    e2 += k2 * d2;
                    k2       = 1.0 / (2.0 * nalpha_ - 1.0) * ( (i==k)*oei_[INDEX(j,l)] + (j==l)*oei_[INDEX(i,k)]);
                    e1 += k2 * d2;
                }
            }
        }
    }

}

void p2RDMSolver::Normalization() {

    int o = nalpha_;
    int v = nmo_ - nalpha_;
    
    if ( options_.get_str("P2RDM_TYPE") == "K" ) {

        for (int i = 0; i < o; i++) {
            for (int j = 0; j < o; j++) {
                for (int a = 0; a < v; a++) {
                    for (int b = 0; b < v; b++) {
                        t0_[a*o*o*v + b*o*o + i*o + j] = 1.0;
                    }
                }
            }
        }

    }else if ( options_.get_str("P2RDM_TYPE") == "CID" ) {

        throw PsiException("implement me",__FILE__,__LINE__);

    }else if ( options_.get_str("P2RDM_TYPE") == "CEPA(0)" ) {

        for (int i = 0; i < o; i++) {
            for (int j = 0; j < o; j++) {
                for (int a = 0; a < v; a++) {
                    for (int b = 0; b < v; b++) {
                        t0_[a*o*o*v + b*o*o + i*o + j] = 1.0;
                    }
                }
            }
        }

    }

}

void p2RDMSolver::Normalization_RDMs() {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    throw PsiException("implement me",__FILE__,__LINE__);
    
    if ( options_.get_str("P2RDM_TYPE") == "K" ) {

        throw PsiException("implement me",__FILE__,__LINE__);

    }else if ( options_.get_str("P2RDM_TYPE") == "CID" ) {

        throw PsiException("implement me",__FILE__,__LINE__);

    }else if ( options_.get_str("P2RDM_TYPE") == "CEPA(0)" ) {

        for (int i = 0; i < o; i++) {
            for (int j = 0; j < o; j++) {
                for (int a = 0; a < v; a++) {
                    for (int b = 0; b < v; b++) {
                        t0_[a*o*o*v + b*o*o + i*o + j] = 1.0;
                    }
                }
            }
        }

    }

}

void p2RDMSolver::setup_integrals() {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    // unpack Qmo into Qoo, Qov, Qvo, Qvv

    int ntri = nmo_*(nmo_+1)/2;

    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            for (int Q = 0; Q < nQ_; Q++) {
                Qoo_[Q*o*o+i*o+j] = Qmo_[Q*ntri + INDEX(i,j)];
            }
        }
    }
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            for (int Q = 0; Q < nQ_; Q++) {
                Qov_[Q*o*v+i*v+a] = Qmo_[Q*ntri + INDEX(i,a+o)];
                Qvo_[Q*o*v+a*o+i] = Qmo_[Q*ntri + INDEX(i,a+o)];
            }
        }
    }
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            for (int Q = 0; Q < nQ_; Q++) {
                Qvv_[Q*v*v+a*v+b] = Qmo_[Q*ntri + INDEX(a+o,b+o)];
            }
        }
    }

    // fock matrix:

    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            double dum = oei_[INDEX(i,j)];
            for (int k = 0; k < o; k++) {
                dum += 2.0 * C_DDOT(nQ_,Qmo_ + INDEX(i,j),nmo_*(nmo_+1)/2,Qmo_+INDEX(k,k),nmo_*(nmo_+1)/2);
                dum -=       C_DDOT(nQ_,Qmo_ + INDEX(i,k),nmo_*(nmo_+1)/2,Qmo_+INDEX(k,j),nmo_*(nmo_+1)/2);
            }
            foo_[i * o + j] = dum;
        }
    }
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            double dum = oei_[INDEX(a+o,b+o)];
            for (int k = 0; k < o; k++) {
                dum += 2.0 * C_DDOT(nQ_,Qmo_ + INDEX(a+o,b+o),nmo_*(nmo_+1)/2,Qmo_+INDEX(k,k),nmo_*(nmo_+1)/2);
                dum -=       C_DDOT(nQ_,Qmo_ + INDEX(a+o,k),nmo_*(nmo_+1)/2,Qmo_+INDEX(k,b+o),nmo_*(nmo_+1)/2);
            }
            fvv_[a * v + b] = dum;
        }
    }

}

} //end namespaces
