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
#include "psi4/libmints/local.h"

#include "david.h"
#include "doci_solver.h"

#include <focas/focas_c_interface.h>
#include <misc/threeindexintegrals.h>
#include <misc/blas.h>
#include <misc/omp.h>
#include <focas/orbital_optimizer.h>

using namespace psi;

namespace hilbert{

static void evaluate_sigma(size_t N, size_t maxdim,double **sigma, double **b, void * data) {
    DOCISolver* doci = reinterpret_cast<DOCISolver*>(data);
    //doci->BuildSigma(N,maxdim,b,sigma);
    doci->BuildSigmaFast(N,maxdim,b,sigma);
}

DOCISolver::DOCISolver(SharedWavefunction reference_wavefunction,Options & options):
    Wavefunction(options){
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

DOCISolver::~DOCISolver()
{

    //free memory allocated by constructure and common_init()

    free(orbopt_outfile_);
    free(orbopt_data_);
    free(orbopt_transformation_matrix_);
    free(eint1_);
    free(eint2_);
    free(Qmo_);
    free(d2_);
    free(oei_);

}

void  DOCISolver::common_init(){

    orbopt_ = (std::shared_ptr<OrbitalOptimizer>)(new OrbitalOptimizer(reference_wavefunction_,options_));

    is_df_ = false;

    if ( options_.get_str("SCF_TYPE") == "DF" || options_.get_str("SCF_TYPE") == "CD" ) {
        is_df_ = true;
    }

    if ( !is_df_ ) {
        throw PsiException("plugin doci only works with scf_type = df or cd, for now",__FILE__,__LINE__);
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
        throw PsiException("plugin doci works only with c1 symmetry",__FILE__,__LINE__);
    }
    nso_      = reference_wavefunction_->nso();
    nmo_      = reference_wavefunction_->nmo();

    if ( nmo_ > 64 ) {
        throw PsiException("there is currently a hard limit of 64 orbitals for DOCI",__FILE__,__LINE__);
    }

    nsopi_    = reference_wavefunction_->nsopi();
    molecule_ = reference_wavefunction_->molecule();
    enuc_     = molecule_->nuclear_repulsion_energy({0.0,0.0,0.0});

    // were there linear dependencies in the primary basis set?
    if ( nmo_ != nso_ ) {
        throw PsiException("linear dependencies.  not sure what to do about that just yet",__FILE__,__LINE__);
    }

    AO2SO_ = SharedMatrix(reference_wavefunction_->aotoso());

    Ca_ = SharedMatrix(reference_wavefunction_->Ca());
    Cb_ = SharedMatrix(reference_wavefunction_->Cb());

    if ( options_.get_bool("LOCALIZE_ORBITALS") ) {
        // localize orbitals:
        std::shared_ptr<PMLocalizer> boys (new PMLocalizer(reference_wavefunction_->basisset(),reference_wavefunction_->Ca_subset("SO","OCC")));
        boys->localize();
        for (size_t mu = 0; mu < nso_; mu++) {
            for (size_t i = 0; i < nalpha_; i++) {
                Ca_->pointer()[mu][i] = boys->L()->pointer()[mu][i];
                Cb_->pointer()[mu][i] = boys->L()->pointer()[mu][i];
            }
        }
    }
    if ( options_.get_bool("LOCALIZE_VIRTUAL_ORBITALS") ) {
        // localize orbitals (virtual):
        std::shared_ptr<PMLocalizer> boys_vir (new PMLocalizer(reference_wavefunction_->basisset(),reference_wavefunction_->Ca_subset("SO","VIR")));
        boys_vir->localize();
        for (size_t mu = 0; mu < nso_; mu++) {
            for (size_t i = nalpha_; i < nso_; i++) {
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

        size_t count = 0;
        do {

            size_t p = (int)( (double)rand()/RAND_MAX * nso_);
            size_t q = (int)( (double)rand()/RAND_MAX * nso_);
            if ( p == q ) continue;

            double angle = ( (double)rand()/RAND_MAX - 1.0 ) * 0.001;
            double cosine = cos(angle);
            double sine   = sin(angle);

            for (size_t mu = 0; mu < nso_; mu++) {

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

    Fa_  = (SharedMatrix)(new Matrix(reference_wavefunction_->Fa()));
    Fb_  = (SharedMatrix)(new Matrix(reference_wavefunction_->Fb()));

    Da_  = (SharedMatrix)(new Matrix(reference_wavefunction_->Da()));
    Db_  = (SharedMatrix)(new Matrix(reference_wavefunction_->Db()));
    
    // Lagrangian matrix
    Lagrangian_ = SharedMatrix(reference_wavefunction_->Lagrangian());

    epsilon_a_= SharedVector(new Vector(nirrep_, nmopi_));
    epsilon_a_->copy(reference_wavefunction_->epsilon_a().get());
    epsilon_b_= SharedVector(new Vector(nirrep_, nmopi_));
    epsilon_b_->copy(reference_wavefunction_->epsilon_b().get());

    // memory is from process::environment
    memory_ = Process::environment.get_memory();

    // set the wavefunction name
    name_ = "DOCI";

    // doci convergence thresholds:
    r_convergence_  = options_.get_double("R_CONVERGENCE");
    e_convergence_  = options_.get_double("E_CONVERGENCE");

    maxiter_        = options_.get_int("MAXITER");

    // memory check happens here

    outfile->Printf("\n\n");
    outfile->Printf( "        ***************************************************\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *    DOCI                                         *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *    doubly occupied configuration interaction    *\n");
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
    for (size_t i = 0; i < nmo_; i++) {
        for (size_t j = i; j < nmo_; j++) {
                oei_[INDEX(i,j)] = Hcore->pointer()[i][j];
        }
    }

    // (ii|jj)
    eint1_ = (double*)malloc(oei_dim_*sizeof(double));
    memset((void*)eint1_,'\0',oei_dim_*sizeof(double));

    // (ij|ji)
    eint2_ = (double*)malloc(oei_dim_*sizeof(double));
    memset((void*)eint2_,'\0',oei_dim_*sizeof(double));

    // stuff for orbital optimization:

    // tei
    tei_dim_ = (size_t) nQ_ * (size_t) nmo_ * ( (size_t) ( nmo_ + 1L ) ) / 2L ;

    tei_ = Qmo_;

    // d2
    d2_dim_ = nmo_ * ( nmo_ + 1 ) / 2  * ( nmo_ * ( nmo_ + 1 ) / 2 + 1 ) / 2;
    d2_ = (double*)malloc(d2_dim_ * sizeof(double));
    memset((void*)d2_,'\0',d2_dim_ * sizeof(double));

    // mo-mo transformation matrix
    newMO_ = (SharedMatrix)(new Matrix(reference_wavefunction_->Ca()));
    newMO_->zero();
    for (size_t h = 0; h < nirrep_; h++) {
        for (size_t i = 0; i < nmopi_[h]; i++) {
            newMO_->pointer(h)[i][i] = 1.0;
        }
    }

    orbopt_transformation_matrix_ = (double*)malloc(nmo_ * nmo_ * sizeof(double));
    memset((void*)orbopt_transformation_matrix_,'\0',nmo_ * nmo_ * sizeof(double));
    for (size_t i = 0; i < nmo_; i++) {
        orbopt_transformation_matrix_[i * nmo_ + i] = 1.0;
    }

    size_t nthread = omp_get_max_threads();

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

}

size_t DOCISolver::GenerateConfigurations(size_t val) {

    size_t c = val & -val;
    size_t r = val + c;
    size_t next = (((r^val) >> 2) / c) | r;

    std::bitset<64> x(next);
    if ( x[nmo_] == 1 ) {
        return next;
    }else {
        configurations_[n_++] = x;
        return GenerateConfigurations(next);
    }

}

size_t popcount64d(uint64_t x)
{
    //return _mm_popcnt_u64(x);
    //return __builtin_popcountll(x);
    size_t count;
    for (count=0; x; ++count) {
        x &= x - 1;
        if ( count > 2 ) return count;
    }
    return count;
}

// Returns value of Binomial Coefficient C(n, k) 
size_t binomialCoeff(size_t n, size_t k) 
{ 
  // Base Cases 
  if (k==0 || k==n) 
    return 1; 
  
  // Recur 
  return  binomialCoeff(n-1, k-1) + binomialCoeff(n-1, k); 
} 


// compute the energy!
double DOCISolver::compute_energy() {

    n_ = binomialCoeff(nmo_,nalpha_);
    configurations_ = (std::bitset<64> * )malloc(n_ * sizeof(std::bitset<64>));

    outfile->Printf("\n");
    outfile->Printf("    Number of DOCI configurations: %10li\n",n_);
    outfile->Printf("\n");

    // diagonal of Hamiltonian
    Hdiag_ = (std::shared_ptr<Vector>)(new Vector(n_));

    n_ = 0;

    size_t nalpha_binary_rep = 0L;
    for (size_t i = 0; i < nalpha_; i++) {
        nalpha_binary_rep += (size_t)pow(2,i);
    }

    // generate unique configurations

    // reference
    std::bitset<64> x(nalpha_binary_rep);
    configurations_[n_++] = x;

    // the rest
    GenerateConfigurations( configurations_[0].to_ulong() );
    

    double * ci_wfn = (double*)malloc(n_*sizeof(double));
    memset((void*)ci_wfn,'\0',n_*sizeof(double));

    size_t iter = 0;

    // which configurations interact with which?
    configuration_map_ = (size_t**)malloc(n_*sizeof(size_t*));
    configuration_map_indices_ = (size_t***)malloc(n_*sizeof(size_t**));
    for (size_t i = 0; i < n_; i++) {
        configuration_map_[i] = (size_t*)malloc((nmo_-nalpha_)*nalpha_*sizeof(size_t));
        configuration_map_indices_[i] = (size_t**)malloc((nmo_-nalpha_)*nalpha_*sizeof(size_t*));
        for (size_t j = 0; j < (nmo_-nalpha_)*nalpha_; j++) {
            configuration_map_indices_[i][j] = (size_t*)malloc(2*sizeof(size_t));
        }
    }
   
    #pragma omp parallel for schedule (dynamic)
    for (size_t i = 0; i < n_; i++) {
        size_t count = 0;
        for (size_t j = 0; j < n_; j++) {
            if ( i == j ) continue;

            std::bitset<64> ixorj = (configurations_[i] ^ configurations_[j]);
            //if ( _mm_popcnt_u64( ixorj.to_ulong() ) != 2 ) continue;
            if ( __builtin_popcountll( ixorj.to_ulong() ) != 2 ) continue;
            //if ( popcount64d( ixorj.to_ulong() ) != 2 ) continue;

            configuration_map_[i][count] = j;

            size_t mybit_int = ixorj.to_ulong();
            size_t myj = __builtin_ffsl(mybit_int) - 1;
            size_t myi = 64 - __builtin_clzl((size_t)(mybit_int)) - 1;

            configuration_map_indices_[i][count][0] = myi;
            configuration_map_indices_[i][count][1] = myj;

            count++;
        }
    }

    d1_ = (double*)malloc( nmo_ * (nmo_+1)/2 * sizeof(double));
    d2ab_ = (double*)malloc(nmo_ * nmo_ * nmo_ * nmo_ * sizeof(double));
    d2aa_ = (double*)malloc(nmo_ * nmo_ * nmo_ * nmo_ * sizeof(double));

    //outfile->Printf("\n");
    outfile->Printf("    ==> Begin oo-DOCI iterations <==\n");
    outfile->Printf("\n");

    outfile->Printf("\n");
    outfile->Printf("       iter");
    outfile->Printf("   iter(CI)");
    outfile->Printf("   iter(oo)");
    outfile->Printf("                E(CI)");
    outfile->Printf("                E(oo)");
    outfile->Printf("                 |dE|\n");

    double dE = 0;
    double energy = -enuc_;
    size_t ci_iter = 0;
    do {

        double ci_energy = DiagonalizeHamiltonian(ci_wfn,ci_iter);
        //printf("wtf ci energy %20.12lf\n",ci_energy+enuc_);
        BuildRDMs(ci_wfn,false);

        double oo_energy;
        if ( options_.get_bool("OPTIMIZE_ORBITALS") ) {

            if ( iter < 50 ) orbopt_data_[14] = 0.0; // steepest descent
            else             orbopt_data_[14] = 3.0; // hager-zhang 

            oo_energy = RotateOrbitals();
            dE = fabs(energy - oo_energy);

            // increase number of orbital optimization steps 
            orbopt_data_[8] += 1.0; // AED

        }else {
            oo_energy = ci_energy;
            orbopt_converged_ = true;
            dE = 0.0;
        }

        energy = oo_energy;

        outfile->Printf("      %5i %9i %11i %20.12lf %20.12lf %20.12lf\n",iter,ci_iter,(int)orbopt_data_[10],ci_energy + enuc_,oo_energy + enuc_,dE);

        if ( orbopt_converged_ && dE < e_convergence_ ) break;

        iter++;

    }while(iter < maxiter_); 

    if ( iter == maxiter_ ) {
        throw PsiException("oo-DOCI did not converge.",__FILE__,__LINE__);
    }

    outfile->Printf("\n");
    outfile->Printf("      oo-DOCI iterations converged!\n");
    outfile->Printf("\n");

    outfile->Printf("    * DOCI total energy:              %20.12lf\n",energy+enuc_);
    outfile->Printf("\n");

    std::shared_ptr<Vector> NOs (new Vector("Natural Orbital Occupation Numbers (spin free)",nmo_));
    BuildRDMs(ci_wfn,true);
    for (size_t i = 0; i < nmo_; i++) {
        NOs->pointer()[i] = d1_[INDEX(i,i)];
    }
    NOs->print();

    Process::environment.globals["CURRENT ENERGY"]    = energy + enuc_;
    Process::environment.globals["DOCI TOTAL ENERGY"] = energy + enuc_;

    // push final transformation matrix onto Ca_ and Cb_
    if ( options_.get_bool("OPTIMIZE_ORBITALS") ) {
        UpdateTransformationMatrix(reference_wavefunction_,newMO_,Ca_,Cb_,orbopt_transformation_matrix_);
    }

    // write tpdm to disk?
    if ( options_.get_bool("TPDM_WRITE_FULL") ) {
        WriteTPDM();
    }

    if ( options_.get_bool("PRINT_RDMS") ) {
        print_rdms();
    }

    // free memory allocated in compute_energy()
    free(configurations_);
    free(ci_wfn);
    for (size_t i = 0; i < n_; i++) {
        free(configuration_map_[i]);
        for (size_t j = 0; j < (nmo_-nalpha_)*nalpha_; j++) {
            free(configuration_map_indices_[i][j]);
        }
        free(configuration_map_indices_[i]);
    }
    free(configuration_map_);
    free(configuration_map_indices_);
    free(d1_);
    free(d2ab_);
    free(d2aa_);

    return energy + enuc_;

}

void DOCISolver::print_rdms() {

    // Gamma(pq,rs) = D2ab(pq,rs) + D2ba(pq,rs) + D2aa(pq,rs) + D2bb(pq,rs)
    outfile->Printf("\n");
    outfile->Printf("    ==> DOCI @2RDM <==\n");
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
    outfile->Printf("    ==> DOCI @1RDM <==\n");
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

void DOCISolver::PackSpatialDensity() {

    size_t n1 = nmo_;
    size_t n2 = nmo_ * nmo_ ;
    size_t n3 = nmo_ * nmo_ * nmo_;

    memset((void*)d2_,'\0',d2_dim_ * sizeof(double));

    for (size_t i = 0; i < nmo_; i++) {
        for (size_t j = 0; j < nmo_; j++) {
            for (size_t k = 0; k < nmo_; k++) {
                for (size_t l = 0; l < nmo_; l++) {

                    size_t id = INDEX(INDEX(i,k),INDEX(j,l));
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

    //double energy = 0.0;

    //for (int i = 0; i < nmo_; i++) {
    //    energy += d1_[INDEX(i,i)] * oei_[INDEX(i,i)];
    //}

    //for (int i = 0; i < nmo_; i++) {
    //    for (int j = i; j < nmo_; j++) {
    //        for (int k = 0; k < nmo_; k++) {
    //            for (int l = k; l < nmo_; l++) {
    //                int ij = INDEX(i,j);
    //                int kl = INDEX(k,l);
    //                if ( ij > kl ) continue;
    //                double eint = C_DDOT(nQ_,Qmo_ + ij,nmo_*(nmo_+1)/2,Qmo_+kl,nmo_*(nmo_+1)/2);
    //                energy += d2_[INDEX(ij,kl)] * eint;
    //            }
    //        }
    //    }
    //}
    //printf("hey! %20.12lf\n",energy + enuc_);
    //exit(0);

}

double DOCISolver::RotateOrbitals() {

    PackSpatialDensity();

    if ( orbopt_data_[8] > 0 ) {
        //outfile->Printf("\n");
        //outfile->Printf("        ==> Orbital Optimization <==\n");
        //outfile->Printf("\n");
    }    

    double * X_ = (double*)malloc(nmo_*nmo_*sizeof(double));

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

          orbopt_data_, orbopt_outfile_, X_);
    //printf("exit\n");fflush(stdout);

    free(X_);
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

    for (size_t i = 0; i < nmo_; i++) {
        e1 += d1_[INDEX(i,i)] * oei_[INDEX(i,i)];
    }

    double e2 = 0.0;

    // D2(ij,ij) (ii|jj)
    for (size_t i = 0; i < nmo_; i++) {
        for (size_t j = 0; j < nmo_; j++) {
            if ( i == j ) continue;
            size_t ii = INDEX(i,i);
            size_t jj = INDEX(j,j);
            double eint = C_DDOT(nQ_,Qmo_ + ii,nmo_*(nmo_+1)/2,Qmo_+jj,nmo_*(nmo_+1)/2);
            e2 += d2ab_[i * nmo_*nmo_*nmo_ + j * nmo_*nmo_ + i * nmo_ + j] * eint;
            e2 += d2aa_[i * nmo_*nmo_*nmo_ + j * nmo_*nmo_ + i * nmo_ + j] * eint;
        }
    }
    // D2(ij,ji) (ij|ji)
    for (size_t i = 0; i < nmo_; i++) {
        for (size_t j = 0; j < nmo_; j++) {
            if ( i == j ) continue;
            size_t ij = INDEX(i,j);
            double eint = C_DDOT(nQ_,Qmo_ + ij,nmo_*(nmo_+1)/2,Qmo_+ij,nmo_*(nmo_+1)/2);
            e2 += d2aa_[i * nmo_*nmo_*nmo_ + j * nmo_*nmo_ + j * nmo_ + i] * eint;
        }
    }
    // D2(ii,jj) (ij|ij)
    for (size_t i = 0; i < nmo_; i++) {
        for (size_t j = 0; j < nmo_; j++) {
            size_t ij = INDEX(i,j);
            double eint = C_DDOT(nQ_,Qmo_ + ij,nmo_*(nmo_+1)/2,Qmo_+ij,nmo_*(nmo_+1)/2);
            e2 += d2ab_[i * nmo_*nmo_*nmo_ + i * nmo_*nmo_ + j * nmo_ + j] * eint;
        }
    }

    //printf("after oo            %20.12lf %20.12lf %20.12lf\n",e1+e2+enuc_,e1,e2);

    /*for (size_t i = 0; i < nmo_; i++) {
        for (size_t j = 0; j < nmo_; j++) {
            for (size_t k = 0; k < nmo_; k++) {
                for (size_t l = 0; l < nmo_; l++) {
                    size_t ij = INDEX(i,j);
                    size_t kl = INDEX(k,l);
                    double eint = C_DDOT(nQ_,Qmo_ + ij,nmo_*(nmo_+1)/2,Qmo_+kl,nmo_*(nmo_+1)/2);
                    e2 += d2ab_[i * nmo_*nmo_*nmo_ + k * nmo_*nmo_ + j * nmo_ + l] * eint;
                    e2 += d2aa_[i * nmo_*nmo_*nmo_ + k * nmo_*nmo_ + j * nmo_ + l] * eint;
                }
            }
        }
    }*/

    // d2 is being overwritten ...
    //PackSpatialDensity();

    return e1 + e2;

}

void DOCISolver::BuildSigmaFast(size_t N, size_t maxdim, double ** b, double ** sigma) {

    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < maxdim; j++) {
            sigma[i][j] = 0.0;
        }
    }

    double * Hdiag_p = Hdiag_->pointer();

    #pragma omp parallel for schedule (dynamic)
    for (size_t i = 0; i < n_; i++) {
        for (size_t jj = 0; jj < (nmo_-nalpha_)*nalpha_; jj++) {

            size_t j = configuration_map_[i][jj];

            size_t myi = configuration_map_indices_[i][jj][0];
            size_t myj = configuration_map_indices_[i][jj][1];

            double eint = eint2_[INDEX(myi,myj)];
            for (size_t k = 0; k < maxdim; k++) {
                sigma[i][k] += eint * b[k][j];
            }

        }

        for (size_t k = 0; k < maxdim; k++) {
            sigma[i][k] += Hdiag_p[i] * b[k][i];
        }

    }

}
void DOCISolver::BuildSigma(size_t N, size_t maxdim, double ** b, double ** sigma) {

    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < maxdim; j++) {
            sigma[i][j] = 0.0;
        }
    }

    double * Hdiag_p = Hdiag_->pointer();

    for (size_t i = 0; i < n_; i++) {
        for (size_t j = i + 1; j < n_; j++) {

            std::bitset<64> ixorj = (configurations_[i] ^ configurations_[j]);
            //if ( _mm_popcnt_u64( ixorj.to_ulong() ) != 2 ) continue;
            if ( __builtin_popcountll( ixorj.to_ulong() ) != 2 ) continue;
            //if ( popcount64d( ixorj.to_ulong() ) != 2 ) continue;

            size_t mybit_int = ixorj.to_ulong();
            size_t myj = __builtin_ffsl(mybit_int) - 1;
            size_t myi = 64 - __builtin_clzl((size_t)(mybit_int)) - 1;

            double eint = eint2_[INDEX(myi,myj)];
            for (size_t k = 0; k < maxdim; k++) {
                sigma[i][k] += eint * b[k][j];
                sigma[j][k] += eint * b[k][i];
            }

        }

        for (size_t k = 0; k < maxdim; k++) {
            sigma[i][k] += Hdiag_p[i] * b[k][i];
        }

    }

}

// return an element of the Hamiltonian.  The diagonal elements must be evaluated first.
double DOCISolver::HamiltonianElement(size_t i, size_t j) {

    double eint = 0.0;
    if ( i != j ) {
        std::bitset<64> ixorj = (configurations_[i] ^ configurations_[j]);
        if ( __builtin_popcountll( ixorj.to_ulong() ) != 2 ) return 0.0;
        
        size_t mybit_int = ixorj.to_ulong();
        size_t myj = __builtin_ffsl(mybit_int) - 1;
        size_t myi = 64 - __builtin_clzl(mybit_int) - 1;
        
        eint = eint2_[INDEX(myi,myj)];
    }else {
        eint = Hdiag_->pointer()[i];
    }
    return eint;
}

double DOCISolver::DiagonalizeHamiltonian(double * ci_wfn, size_t & ci_iter) {

    // build selected integrals that will contribute to H: (ii|jj), (ij|ji)
    for (size_t i = 0; i < nmo_; i++) {
        for (size_t j = i; j < nmo_; j++) {
            eint1_[INDEX(i,j)] = C_DDOT(nQ_,Qmo_ + INDEX(i,i),nmo_*(nmo_+1)/2,Qmo_+INDEX(j,j),nmo_*(nmo_+1)/2);
            eint2_[INDEX(i,j)] = C_DDOT(nQ_,Qmo_ + INDEX(i,j),nmo_*(nmo_+1)/2,Qmo_+INDEX(j,i),nmo_*(nmo_+1)/2);
        }
    }

    for (size_t i = 0; i < n_; i++) {

        double val = 0.0;
        for (size_t k = 0; k < nmo_; k++) {
            if ( !configurations_[i][k] )  continue;
            for (size_t l = 0; l < nmo_; l++) {
                if ( !configurations_[i][l] )  continue;

                val += 2.0 * eint1_[INDEX(k,l)] - eint2_[INDEX(k,l)];

            }
            val += 2.0 * oei_[INDEX(k,k)];
            
        }
        Hdiag_->pointer()[i] = val;

    }

    // diagonalize!
    std::shared_ptr<Matrix> eigvec (new Matrix(1,n_));
    std::shared_ptr<Vector> eigval (new Vector(n_));
    //david_in_core(H->pointer(),n_,1,eigval->pointer(),eigvec->pointer(),r_convergence_,0,ci_iter);

    //david_direct(Hdiag_->pointer(),n_,1,eigval->pointer(),eigvec->pointer(),r_convergence_,0,evaluate_sigma,ci_iter,(void*)this);
    david_direct_redo(Hdiag_->pointer(),n_,1,eigval->pointer(),eigvec->pointer(),r_convergence_,0,evaluate_sigma,ci_iter,(void*)this,options_.get_int("DAVIDSON_MAXDIM"));

    C_DCOPY(n_,eigvec->pointer()[0],1,ci_wfn,1);

    //printf("what is the norm: %20.12lf\n",C_DNRM2(n_,ci_wfn,1));fflush(stdout);

    //std::shared_ptr<Matrix> eigvec (new Matrix(n_,n_));
    //std::shared_ptr<Vector> eigval (new Vector(n_));
    //H->diagonalize(eigvec,eigval);
    //for (int i = 0; i < n_; i++) {
    //    ci_wfn[i] = eigvec->pointer()[i][0];
    //}

    return eigval->pointer()[0];

}

void DOCISolver::BuildRDMs(double * ci_wfn, bool print){

    memset((void*)d1_,'\0', oei_dim_ *  sizeof(double));
    memset((void*)d2ab_,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)d2aa_,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    double energy = 0.0;

    // build 2-RDM
    for (size_t i = 0; i < n_; i++) {
        for (size_t jj = 0; jj < (nmo_-nalpha_)*nalpha_; jj++) {

            size_t j = configuration_map_[i][jj];

            size_t myi = configuration_map_indices_[i][jj][0];
            size_t myj = configuration_map_indices_[i][jj][1];

            double cij = ci_wfn[i] * ci_wfn[j];

            d2ab_[myi * nmo_*nmo_*nmo_ + myi * nmo_*nmo_ + myj * nmo_ + myj] += cij;

            double eint = eint2_[INDEX(myi,myj)];
            energy += eint * cij;

        }

        double cii = ci_wfn[i] * ci_wfn[i];

        double val = 0.0;
        for (size_t k = 0; k < nmo_; k++) {
            if ( !configurations_[i][k] )  continue;
            for (size_t l = 0; l < nmo_; l++) {
                if ( !configurations_[i][l] )  continue;

                d2ab_[k * nmo_*nmo_*nmo_ + l * nmo_*nmo_ + k * nmo_ + l] += cii;
                d2aa_[k * nmo_*nmo_*nmo_ + l * nmo_*nmo_ + k * nmo_ + l] += cii;
                d2aa_[k * nmo_*nmo_*nmo_ + l * nmo_*nmo_ + l * nmo_ + k] -= cii;

            }
            d1_[INDEX(k,k)] += 2.0 * cii;

        }
        energy += Hdiag_->pointer()[i] * cii;

    }
    //printf("build rdm (CI way): %20.12lf\n",energy+enuc_);

    double e1 = 0.0;

    for (size_t i = 0; i < nmo_; i++) {
        e1 += d1_[INDEX(i,i)] * oei_[INDEX(i,i)];
    }

    double e2 = 0.0;

    // D2(ij,ij) (ii|jj)
    for (size_t i = 0; i < nmo_; i++) {
        for (size_t j = 0; j < nmo_; j++) {
            if ( i == j ) continue;
            size_t ii = INDEX(i,i);
            size_t jj = INDEX(j,j);
            double eint = C_DDOT(nQ_,Qmo_ + ii,nmo_*(nmo_+1)/2,Qmo_+jj,nmo_*(nmo_+1)/2);
            e2 += d2ab_[i * nmo_*nmo_*nmo_ + j * nmo_*nmo_ + i * nmo_ + j] * eint;
            e2 += d2aa_[i * nmo_*nmo_*nmo_ + j * nmo_*nmo_ + i * nmo_ + j] * eint;
        }
    }
    // D2(ij,ji) (ij|ji)
    for (size_t i = 0; i < nmo_; i++) {
        for (size_t j = 0; j < nmo_; j++) {
            if ( i == j ) continue;
            size_t ij = INDEX(i,j);
            double eint = C_DDOT(nQ_,Qmo_ + ij,nmo_*(nmo_+1)/2,Qmo_+ij,nmo_*(nmo_+1)/2);
            e2 += d2aa_[i * nmo_*nmo_*nmo_ + j * nmo_*nmo_ + j * nmo_ + i] * eint;
        }
    }
    // D2(ii,jj) (ij|ij)
    for (size_t i = 0; i < nmo_; i++) {
        for (size_t j = 0; j < nmo_; j++) {
            size_t ij = INDEX(i,j);
            double eint = C_DDOT(nQ_,Qmo_ + ij,nmo_*(nmo_+1)/2,Qmo_+ij,nmo_*(nmo_+1)/2);
            e2 += d2ab_[i * nmo_*nmo_*nmo_ + i * nmo_*nmo_ + j * nmo_ + j] * eint;
        }
    }

    if ( print ) {
        outfile->Printf("\n");
        outfile->Printf("        DOCI one-electron energy = %20.12lf\n",e1);
        outfile->Printf("        DOCI two-electron energy = %20.12lf\n",e2);
        outfile->Printf("        * DOCI total energy      = %20.12lf\n",e1 + e2 + enuc_);
        outfile->Printf("\n");
    }

}

} //end namespaces
