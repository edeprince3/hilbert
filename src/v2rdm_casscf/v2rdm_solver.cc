/*
 *@BEGIN LICENSE
 *
 * v2RDM-CASSCF, a plugin to:
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

#include <psi4/libmints/mintshelper.h>
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
#include <psi4/libmints/local.h>
#include <psi4/libmints/molecule.h>
#include <psi4/libmints/factory.h>
#include <psi4/libmints/basisset.h>

#include "v2rdm_solver.h"

#include <focas/focas_c_interface.h>
#include <misc/blas.h>
#include <misc/cg_solver.h>
#include <misc/bpsdp_solver.h>
#include <misc/threeindexintegrals.h>
#include <misc/omp.h>

using namespace psi;
using namespace fnocc;

extern "C" {
    void dgeev(char& JOBVL,char& JOBVR,long int& N,double* A,long int& LDA,double* WR,double* WI,
            double * VL,long int& LDVL,double* VR,long int& LDVR,double* WORK,long int& LWORK,long int& INFO);
};
inline void DGEEV(char& JOBVL,char& JOBVR,long int& N,double* A,long int& LDA,double* WR,double* WI,
            double * VL,long int& LDVL,double* VR,long int& LDVR,double* WORK,long int& LWORK,long int& INFO){
    dgeev(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO);
}


// diagonalize real, nonsymmetric matrix
void NonsymmetricEigenvalue(long int N, double * A, double * VL, double * VR, double * WR, double *WI){

    char JOBVL = 'V';
    char JOBVR = 'V';
    long int LDA  = N;
    long int LDVL = N;
    long int LDVR = N;
    long int LWORK = 4*N;
    double * WORK = (double*)malloc(LWORK*sizeof(double));
    long int INFO;

    DGEEV(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO);

    // kill complex eigenvalues
    for (int i = 0; i < N; i++) {
        if ( fabs(WI[i]) > 1e-6 ) {
            WR[i] = 0.0;
            WI[i] = 0.0;
        }
    }

    free(WORK);
}

namespace hilbert{

static void evaluate_Au(SharedVector Au, SharedVector u, void * data) {

    // reinterpret void * as an instance of v2RDMSolver
    v2RDMSolver* v2rdm = reinterpret_cast<v2RDMSolver*>(data);
    // call a function from class to evaluate Ax product:
    v2rdm->bpsdp_Au(Au,u);

}

static void evaluate_ATu(SharedVector ATu, SharedVector u, void * data) {

    // reinterpret void * as an instance of v2RDMSolver
    v2RDMSolver* v2rdm = reinterpret_cast<v2RDMSolver*>(data);
    // call a function from class to evaluate Ax product:
    v2rdm->bpsdp_ATu(ATu,u);

}


static void evaluate_cg_lhs(SharedVector Ax, SharedVector x, void * data) {

    // reinterpret void * as an instance of v2RDMSolver
    v2RDMSolver* v2rdm = reinterpret_cast<v2RDMSolver*>(data);

    Ax->zero();
    std::shared_ptr<Vector> ATy (new Vector(v2rdm->n_primal()));
    v2rdm->bpsdp_ATu(ATy,x);
    v2rdm->bpsdp_Au(Ax,ATy);

}

v2RDMSolver::v2RDMSolver(SharedWavefunction reference_wavefunction,Options & options):
    Wavefunction(options){
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

v2RDMSolver::~v2RDMSolver()
{
    free(tei_full_sym_);
    free(oei_full_sym_);
    free(d2_plus_core_sym_);
    free(d1_act_spatial_sym_);

    free(amopi_);
    free(rstcpi_);
    free(rstvpi_);
    free(d2aboff);
    free(d2aaoff);
    free(d2bboff);
    free(d1aoff);
    free(d1boff);
    free(q1aoff);
    free(q1boff);
    if ( constrain_q2_ ) {
        if ( !spin_adapt_q2_ ) {
            free(q2aboff);
            free(q2aaoff);
            free(q2bboff);
        }else {
            free(q2soff);
            free(q2toff);
            free(q2toff_p1);
            free(q2toff_m1);
        }
    }
    if ( constrain_g2_ ) {
        if ( ! spin_adapt_g2_ ) {
            free(g2aboff);
            free(g2baoff);
            free(g2aaoff);
        }else {
            free(g2soff);
            free(g2toff);
            free(g2toff_p1);
            free(g2toff_m1);
        }
    }
    if ( constrain_t1_ ) {
        free(t1aaboff);
        free(t1bbaoff);
        free(t1aaaoff);
        free(t1bbboff);
    }
    if ( constrain_t2_ ) {
        free(t2aaboff);
        free(t2bbaoff);
        free(t2aaaoff);
        free(t2bbboff);
    }
    if ( constrain_d3_ ) {
        free(d3aaaoff);
        free(d3bbboff);
        free(d3aaboff);
        free(d3bbaoff);
    }
    if ( constrain_d4_ ) {
        free(d4aaaaoff);
        free(d4aaaboff);
        free(d4aabboff);
        free(d4bbbaoff);
        free(d4bbbboff);
    }

    free(X_);

}

void  v2RDMSolver::common_init(){

    is_df_ = false;
    if ( options_.get_str("SCF_TYPE") == "DF" || options_.get_str("SCF_TYPE") == "CD" ) {
        is_df_ = true;
    }

    //if ( options_.get_bool("EXTENDED_KOOPMANS") && options_.get_bool("NAT_ORBS") ) {
    //    throw PsiException("EKT does not work with natural orbitals",__FILE__,__LINE__);
    //}
    //if ( options_.get_bool("EXTENDED_KOOPMANS") && options_.get_bool("FCIDUMP") ) {
    //    throw PsiException("EKT does not work with natural orbitals (triggered by FCIDUMP=true)",__FILE__,__LINE__);
    //}

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
    nso_      = reference_wavefunction_->nso();
    nmo_      = reference_wavefunction_->nmo();
    nsopi_    = reference_wavefunction_->nsopi();
    molecule_ = reference_wavefunction_->molecule();
    enuc_     = molecule_->nuclear_repulsion_energy({0.0,0.0,0.0});

    // need somewhere to store gradient, if required
    gradient_ =  reference_wavefunction_->matrix_factory()->create_shared_matrix("Total gradient", molecule_->natom(), 3);

    // restricted doubly occupied orbitals per irrep (optimized)
    rstcpi_   = (int*)malloc(nirrep_*sizeof(int));
    memset((void*)rstcpi_,'\0',nirrep_*sizeof(int));

    // restricted unoccupied occupied orbitals per irrep (optimized)
    rstvpi_   = (int*)malloc(nirrep_*sizeof(int));
    memset((void*)rstvpi_,'\0',nirrep_*sizeof(int));

    // active orbitals per irrep:
    amopi_    = (int*)malloc(nirrep_*sizeof(int));
    memset((void*)amopi_,'\0',nirrep_*sizeof(int));

    // multiplicity:
    multiplicity_ = Process::environment.molecule()->multiplicity();

    if (options_["FROZEN_DOCC"].has_changed()) {
        //gg need to take this out to allow frozen doubly occupied orbitals
        //gg-fc 
        //throw PsiException("FROZEN_DOCC is currently disabled.",__FILE__,__LINE__);

        if (options_["FROZEN_DOCC"].size() != nirrep_) {
            throw PsiException("The FROZEN_DOCC array has the wrong dimensions_",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            frzcpi_[h] = options_["FROZEN_DOCC"][h].to_double();
        }
    }
    if (options_["RESTRICTED_DOCC"].has_changed()) {
        if (options_["RESTRICTED_DOCC"].size() != nirrep_) {
            throw PsiException("The RESTRICTED_DOCC array has the wrong dimensions_",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            rstcpi_[h] = options_["RESTRICTED_DOCC"][h].to_double();
        }
    }
    if (options_["RESTRICTED_UOCC"].has_changed()) {
        if (options_["RESTRICTED_UOCC"].size() != nirrep_) {
            throw PsiException("The RESTRICTED_UOCC array has the wrong dimensions_",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            rstvpi_[h] = options_["RESTRICTED_UOCC"][h].to_double();
        }
    }
    if (options_["FROZEN_UOCC"].has_changed()) {

        //if ( !is_df_ ) {
        //    throw PsiException("FROZEN_UOCC is currently enabled only for SCF_TYPE CD and DF.",__FILE__,__LINE__);
        //}
        if (options_["FROZEN_UOCC"].size() != nirrep_) {
            throw PsiException("The FROZEN_UOCC array has the wrong dimensions_",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            frzvpi_[h] = options_["FROZEN_UOCC"][h].to_double();
        }
    }

    // user could specify active space with ACTIVE array
    if ( options_["ACTIVE"].has_changed() ) {
        //throw PsiException("The ACTIVE array is not yet enabled.",__FILE__,__LINE__);
        if (options_["ACTIVE"].size() != nirrep_) {
            throw PsiException("The ACTIVE array has the wrong dimensions_",__FILE__,__LINE__);
        }

        // warn user that active array takes precedence over restricted_uocc array
        if (options_["RESTRICTED_UOCC"].has_changed()) {
            outfile->Printf("\n");
            outfile->Printf("    <<< WARNING!! >>>\n");
            outfile->Printf("\n");
            outfile->Printf("    The ACTIVE array takes precedence over the RESTRICTED_UOCC array.\n");
            outfile->Printf("    Check below whether your active space was correctly specified.\n");
            outfile->Printf("\n");
        }

        // overwrite rstvpi_ array using the information in the frozen_docc,
        // restricted_docc, active, and frozen_virtual arrays.  start with nso total
        // orbitals and let the linear dependency check below adjust the spaces as needed
        for (int h = 0; h < nirrep_; h++) {
            amopi_[h]  = options_["ACTIVE"][h].to_double();
            rstvpi_[h] = nsopi_[h] - frzcpi_[h] - rstcpi_[h] - frzvpi_[h] - amopi_[h];
        }
    }

    // were there linear dependencies in the primary basis set?
    if ( nmo_ != nso_ ) {

        // which irreps lost orbitals?
        int * lost = (int*)malloc(nirrep_*sizeof(int));
        memset((void*)lost,'\0',nirrep_*sizeof(int));
        bool active_space_changed = false;
        for (int h = 0; h < factory_->nirrep(); h++){
            lost[h] = nsopi_[h] - nmopi_[h];
            if ( lost[h] > 0 ) {
                active_space_changed = true;
            }

            // eliminate frozen virtual orbitals first
            if ( frzvpi_[h] > 0 && lost[h] > 0 ) {
                frzvpi_[h] -= ( frzvpi_[h] < lost[h] ? frzvpi_[h] : lost[h] );
                lost[h]    -= ( frzvpi_[h] < lost[h] ? frzvpi_[h] : lost[h] );
            }
            // if necessary, eliminate restricted virtual orbitals next
            if ( rstvpi_[h] > 0 && lost[h] > 0 ) {
                rstvpi_[h] -= ( rstvpi_[h] < lost[h] ? rstvpi_[h] : lost[h] );
            }
        }
        if ( active_space_changed ) {
            outfile->Printf("\n");
            outfile->Printf("    <<< WARNING!! >>>\n");
            outfile->Printf("\n");
            outfile->Printf("    Your basis set may have linear dependencies.\n");
            outfile->Printf("    The number of restricted or frozen virtual orbitals per irrep may have changed.\n");
            outfile->Printf("\n");
            outfile->Printf("    No. orbitals removed per irrep: [");
            for (int h = 0; h < nirrep_; h++)
                outfile->Printf("%4i",nsopi_[h] - nmopi_[h]);
            outfile->Printf(" ]\n");
            //outfile->Printf("    No. frozen virtuals per irrep:  [");
            //for (int h = 0; h < nirrep_; h++)
            //    outfile->Printf("%4i",frzvpi_[h]);
            //outfile->Printf(" ]\n");
            //outfile->Printf("\n");
            outfile->Printf("    Check that your active space is still correct.\n");
            outfile->Printf("\n");
        }
    }

    AO2SO_ = SharedMatrix(reference_wavefunction_->aotoso());

    Ca_ = SharedMatrix(reference_wavefunction_->Ca());
    Cb_ = SharedMatrix(reference_wavefunction_->Cb());

    if ( options_.get_bool("LOCALIZE_ORBITALS") ) {
        // localize orbitals:
        std::shared_ptr<BoysLocalizer> boys (new BoysLocalizer(reference_wavefunction_->basisset(),reference_wavefunction_->Ca_subset("SO","OCC")));
        boys->localize();
        for (int mu = 0; mu < nso_; mu++) {
            for (int i = 0; i < nalpha_; i++) {
                Ca_->pointer()[mu][i] = boys->L()->pointer()[mu][i];
                Cb_->pointer()[mu][i] = boys->L()->pointer()[mu][i];
            }
        }
        // localize orbitals (virtual):
        std::shared_ptr<BoysLocalizer> boys_vir (new BoysLocalizer(reference_wavefunction_->basisset(),reference_wavefunction_->Ca_subset("SO","VIR")));
        boys_vir->localize();
        for (int mu = 0; mu < nso_; mu++) {
            for (int i = (int)nalpha_; i < nso_; i++) {
                Ca_->pointer()[mu][i] = boys_vir->L()->pointer()[mu][i - (int)nalpha_];
                Cb_->pointer()[mu][i] = boys_vir->L()->pointer()[mu][i - (int)nalpha_];
            }
        }
    }

    double fractional_charge = options_.get_double("FRACTIONAL_CHARGE");
    if ( fractional_charge > 0.0 ) { 
        nalpha_ -= options_.get_double("FRACTIONAL_CHARGE");
    }else {
        nbeta_ -= options_.get_double("FRACTIONAL_CHARGE");
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

    amo_      = 0;
    nfrzc_    = 0;
    nfrzv_    = 0;
    nrstc_    = 0;
    nrstv_    = 0;

    int ndocc = 0;
    int nvirt = 0;
    for (int h = 0; h < nirrep_; h++){
        nfrzc_   += frzcpi_[h];
        nrstc_   += rstcpi_[h];
        nrstv_   += rstvpi_[h];
        nfrzv_   += frzvpi_[h];
        amo_   += nmopi_[h]-frzcpi_[h]-rstcpi_[h]-rstvpi_[h]-frzvpi_[h];
        ndocc    += doccpi_[h];
        amopi_[h] = nmopi_[h]-frzcpi_[h]-rstcpi_[h]-rstvpi_[h]-frzvpi_[h];
    }

    int ndoccact = ndocc - nfrzc_ - nrstc_;
    nvirt    = amo_ - ndoccact;

    // sanity check for orbital occupancies:
    for (int h = 0; h < nirrep_; h++) {
        int tot = doccpi_[h] + soccpi_[h] + rstvpi_[h] + frzvpi_[h];
        if (doccpi_[h] + soccpi_[h] + rstvpi_[h] + frzvpi_[h] > nmopi_[h] ) {
            outfile->Printf("\n");
            outfile->Printf("    <<< WARNING >>> irrep %5i has too many orbitals:\n",h);
            outfile->Printf("\n");
            outfile->Printf("                    docc = %5i\n",doccpi_[h]);
            outfile->Printf("                    socc = %5i\n",soccpi_[h]);
            outfile->Printf("                    rstu = %5i\n",rstvpi_[h]);
            outfile->Printf("                    frzv = %5i\n",frzvpi_[h]);
            outfile->Printf("                    tot  = %5i\n",doccpi_[h] + soccpi_[h] + rstvpi_[h] + frzvpi_[h]);
            outfile->Printf("\n");
            outfile->Printf("                    total no. orbitals should be %5i\n",nmopi_[h]);
            outfile->Printf("\n");
            throw PsiException("at least one irrep has too many orbitals",__FILE__,__LINE__);
        }
        if (frzcpi_[h] + rstcpi_[h] > doccpi_[h] ) {
            outfile->Printf("\n");
            outfile->Printf("    <<< WARNING >>> irrep %5i has too many frozen and restricted core orbitals:\n",h);
            outfile->Printf("                    frzc = %5i\n",frzcpi_[h]);
            outfile->Printf("                    rstd = %5i\n",rstcpi_[h]);
            outfile->Printf("                    docc = %5i\n",doccpi_[h]);
            outfile->Printf("\n");
            throw PsiException("at least one irrep has too many frozen core orbitals",__FILE__,__LINE__);
        }
    }

    // memory is from process::environment
    memory_ = Process::environment.get_memory();
    // set the wavefunction name
    name_ = "V2RDM CASSCF";

    // pick conditions.  default is dqg
    constrain_q2_ = true;
    constrain_g2_ = true;
    constrain_t1_ = false;
    constrain_t2_ = false;
    constrain_d3_ = false;
    constrain_d4_ = false;
    if (options_.get_str("POSITIVITY")=="D") {
        constrain_q2_ = false;
        constrain_g2_ = false;
    }else if (options_.get_str("POSITIVITY")=="DQ") {
        constrain_q2_ = true;
        constrain_g2_ = false;
    }else if (options_.get_str("POSITIVITY")=="DG") {
        constrain_q2_ = false;
        constrain_g2_ = true;
    }else if (options_.get_str("POSITIVITY")=="DQGT1") {
        constrain_q2_ = true;
        constrain_g2_ = true;
        constrain_t1_ = true;
    }else if (options_.get_str("POSITIVITY")=="DQGT2") {
        constrain_q2_ = true;
        constrain_g2_ = true;
        constrain_t2_ = true;
    }else if (options_.get_str("POSITIVITY")=="DQGT1T2") {
        constrain_q2_ = true;
        constrain_g2_ = true;
        constrain_t1_ = true;
        constrain_t2_ = true;
    }else if (options_.get_str("POSITIVITY")=="DQGT") {
        constrain_q2_ = true;
        constrain_g2_ = true;
        constrain_t1_ = true;
        constrain_t2_ = true;
    }

    if ( options_.get_bool("CONSTRAIN_D3") ) {
        constrain_d3_ = true;
    }
    if ( options_.get_bool("CONSTRAIN_D4") ) {
        constrain_d4_ = true;
        if ( !constrain_d3_ ) {
            constrain_d3_ = true;
        }
    }

    spin_adapt_g2_  = options_.get_bool("SPIN_ADAPT_G2");
    spin_adapt_q2_  = options_.get_bool("SPIN_ADAPT_Q2");
    constrain_spin_ = options_.get_bool("CONSTRAIN_SPIN");

    if ( constrain_t1_ || constrain_t2_ ) {
        if (spin_adapt_g2_) {
            throw PsiException("If constraining T1/T2, G2 cannot currently be spin adapted.",__FILE__,__LINE__);
        }
        if (spin_adapt_q2_) {
            throw PsiException("If constraining T1/T2, Q2 cannot currently be spin adapted.",__FILE__,__LINE__);
        }
    }

    // build mapping arrays and determine the number of geminals per block
    BuildBasis();

    double ms = (multiplicity_ - 1.0)/2.0;
    if ( ms > 0 ) {
        if (spin_adapt_g2_) {
            throw PsiException("G2 not spin adapted for S = M != 0",__FILE__,__LINE__);
        }
        if (spin_adapt_q2_) {
            throw PsiException("Q2 not spin adapted for S = M != 0",__FILE__,__LINE__);
        }
    }

    // dimension of variable buffer (x)
    dimx_ = 0;
    for ( int h = 0; h < nirrep_; h++) {
        dimx_ += gems_ab[h]*gems_ab[h]; // D2ab
    }
    for ( int h = 0; h < nirrep_; h++) {
        dimx_ += gems_aa[h]*gems_aa[h]; // D2aa
    }
    for ( int h = 0; h < nirrep_; h++) {
        dimx_ += gems_aa[h]*gems_aa[h]; // D2bb
    }
    for ( int h = 0; h < nirrep_; h++) {
        dimx_ += amopi_[h]*amopi_[h]; // D1a
        dimx_ += amopi_[h]*amopi_[h]; // D1b
        dimx_ += amopi_[h]*amopi_[h]; // Q1b
        dimx_ += amopi_[h]*amopi_[h]; // Q1a
    }
    if ( constrain_spin_ && nalpha_ == nbeta_ ) {
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += gems_ab[h] * gems_ab[h]; // D200
        }
    }else if ( constrain_spin_ ) {
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += 4 * gems_ab[h] * gems_ab[h]; // D200
        }
    }
    if ( constrain_q2_ ) {
        if ( !spin_adapt_q2_ ) {
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_ab[h]*gems_ab[h]; // Q2ab
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_aa[h]*gems_aa[h]; // Q2aa
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_aa[h]*gems_aa[h]; // Q2bb
            }
        }else {
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_00[h]*gems_00[h]; // Q2s
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_aa[h]*gems_aa[h]; // Q2t
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_aa[h]*gems_aa[h]; // Q2t_p1
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_aa[h]*gems_aa[h]; // Q2t_m1
            }
        }
    }
    if ( constrain_g2_ ) {
        if ( !spin_adapt_g2_ ) {
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_ab[h]*gems_ab[h]; // G2ab
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_ab[h]*gems_ab[h]; // G2ba
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += 2*gems_ab[h]*2*gems_ab[h]; // G2aa/bb
            }
        }else {
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_ab[h]*gems_ab[h]; // G2s
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_ab[h]*gems_ab[h]; // G2t
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_ab[h]*gems_ab[h]; // G2t_p1
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_ab[h]*gems_ab[h]; // G2t_m1
            }
        }
    }
    if ( constrain_t1_ ) {
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aaa[h]*trip_aaa[h]; // T1aaa
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aaa[h]*trip_aaa[h]; // T1bbb
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aab[h]*trip_aab[h]; // T1aab
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aab[h]*trip_aab[h]; // T1bba
        }
    }
    if ( constrain_t2_ ) {
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += (trip_aba[h]+trip_aab[h])*(trip_aab[h]+trip_aba[h]); // T2aaa
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += (trip_aba[h]+trip_aab[h])*(trip_aab[h]+trip_aba[h]); // T2bbb
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aab[h]*trip_aab[h]; // T2aab
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aab[h]*trip_aab[h]; // T2bba
        }
    }
    if ( constrain_d3_ ) {
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aaa[h] * trip_aaa[h]; // D3aaa
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aaa[h] * trip_aaa[h]; // D3bbb
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aab[h]*trip_aab[h]; // D3aab
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aab[h]*trip_aab[h]; // D3bba
        }
    }
    if ( constrain_d4_ ) {
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += quartet_aaaa[h] * quartet_aaaa[h]; // D4aaaa
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += quartet_aaab[h] * quartet_aaab[h]; // D4aaab
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += quartet_aabb[h] * quartet_aabb[h]; // D4aabb
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += quartet_aaab[h] * quartet_aaab[h]; // D4bbba
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += quartet_aaaa[h] * quartet_aaaa[h]; // D4bbbb
        }
    }

    // offsets in x
    offset = 0;

    d2aboff = (int*)malloc(nirrep_*sizeof(int));
    d2aaoff = (int*)malloc(nirrep_*sizeof(int));
    d2bboff = (int*)malloc(nirrep_*sizeof(int));
    d200off = (int*)malloc(nirrep_*sizeof(int));
    for (int h = 0; h < nirrep_; h++) {
        d2aboff[h] = offset; offset += gems_ab[h]*gems_ab[h];
    }
    for (int h = 0; h < nirrep_; h++) {
        d2aaoff[h] = offset; offset += gems_aa[h]*gems_aa[h];
    }
    for (int h = 0; h < nirrep_; h++) {
        d2bboff[h] = offset; offset += gems_aa[h]*gems_aa[h];
    }
    if ( constrain_spin_ && nalpha_ == nbeta_ ) {
        for (int h = 0; h < nirrep_; h++) {
            d200off[h] = offset; offset += gems_ab[h]*gems_ab[h];
        }
    } else if ( constrain_spin_ ) {
        for (int h = 0; h < nirrep_; h++) {
            d200off[h] = offset; offset += 4*gems_ab[h]*gems_ab[h];
        }
    }

    d1aoff = (int*)malloc(nirrep_*sizeof(int));
    d1boff = (int*)malloc(nirrep_*sizeof(int));
    q1aoff = (int*)malloc(nirrep_*sizeof(int));
    q1boff = (int*)malloc(nirrep_*sizeof(int));
    for (int h = 0; h < nirrep_; h++) {
        d1aoff[h] = offset; offset += amopi_[h]*amopi_[h];
    }
    for (int h = 0; h < nirrep_; h++) {
        d1boff[h] = offset; offset += amopi_[h]*amopi_[h];
    }
    for (int h = 0; h < nirrep_; h++) {
        q1aoff[h] = offset; offset += amopi_[h]*amopi_[h];
    }
    for (int h = 0; h < nirrep_; h++) {
        q1boff[h] = offset; offset += amopi_[h]*amopi_[h];
    }

    if ( constrain_q2_ ) {
        if ( !spin_adapt_q2_ ) {
            q2aboff = (int*)malloc(nirrep_*sizeof(int));
            q2aaoff = (int*)malloc(nirrep_*sizeof(int));
            q2bboff = (int*)malloc(nirrep_*sizeof(int));
            for (int h = 0; h < nirrep_; h++) {
                q2aboff[h] = offset; offset += gems_ab[h]*gems_ab[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                q2aaoff[h] = offset; offset += gems_aa[h]*gems_aa[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                q2bboff[h] = offset; offset += gems_aa[h]*gems_aa[h];
            }
        }else {
            q2soff = (int*)malloc(nirrep_*sizeof(int));
            q2toff = (int*)malloc(nirrep_*sizeof(int));
            q2toff_p1 = (int*)malloc(nirrep_*sizeof(int));
            q2toff_m1 = (int*)malloc(nirrep_*sizeof(int));
            for (int h = 0; h < nirrep_; h++) {
                q2soff[h] = offset; offset += gems_00[h]*gems_00[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                q2toff[h] = offset; offset += gems_aa[h]*gems_aa[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                q2toff_p1[h] = offset; offset += gems_aa[h]*gems_aa[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                q2toff_m1[h] = offset; offset += gems_aa[h]*gems_aa[h];
            }
        }
    }

    if ( constrain_g2_ ) {
        if ( ! spin_adapt_g2_ ) {
            g2aboff = (int*)malloc(nirrep_*sizeof(int));
            g2baoff = (int*)malloc(nirrep_*sizeof(int));
            g2aaoff = (int*)malloc(nirrep_*sizeof(int));
            for (int h = 0; h < nirrep_; h++) {
                g2aboff[h] = offset; offset += gems_ab[h]*gems_ab[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                g2baoff[h] = offset; offset += gems_ab[h]*gems_ab[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                g2aaoff[h] = offset; offset += 2*gems_ab[h]*2*gems_ab[h];
            }
        }else {
            g2soff = (int*)malloc(nirrep_*sizeof(int));
            g2toff = (int*)malloc(nirrep_*sizeof(int));
            g2toff_p1 = (int*)malloc(nirrep_*sizeof(int));
            g2toff_m1 = (int*)malloc(nirrep_*sizeof(int));
            for (int h = 0; h < nirrep_; h++) {
                g2soff[h] = offset; offset += gems_ab[h]*gems_ab[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                g2toff[h] = offset; offset += gems_ab[h]*gems_ab[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                g2toff_p1[h] = offset; offset += gems_ab[h]*gems_ab[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                g2toff_m1[h] = offset; offset += gems_ab[h]*gems_ab[h];
            }
        }
    }

    if ( constrain_t1_ ) {
        t1aaboff = (int*)malloc(nirrep_*sizeof(int));
        t1bbaoff = (int*)malloc(nirrep_*sizeof(int));
        t1aaaoff = (int*)malloc(nirrep_*sizeof(int));
        t1bbboff = (int*)malloc(nirrep_*sizeof(int));
        for (int h = 0; h < nirrep_; h++) {
            t1aaaoff[h] = offset; offset += trip_aaa[h]*trip_aaa[h]; // T1aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            t1bbboff[h] = offset; offset += trip_aaa[h]*trip_aaa[h]; // T1bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            t1aaboff[h] = offset; offset += trip_aab[h]*trip_aab[h]; // T1aab
        }
        for (int h = 0; h < nirrep_; h++) {
            t1bbaoff[h] = offset; offset += trip_aab[h]*trip_aab[h]; // T1bba
        }
    }

    if ( constrain_t2_ ) {
        t2aaboff = (int*)malloc(nirrep_*sizeof(int));
        t2bbaoff = (int*)malloc(nirrep_*sizeof(int));
        t2aaaoff = (int*)malloc(nirrep_*sizeof(int));
        t2bbboff = (int*)malloc(nirrep_*sizeof(int));
        for (int h = 0; h < nirrep_; h++) {
            t2aaaoff[h] = offset; offset += (trip_aab[h]+trip_aba[h])*(trip_aab[h]+trip_aba[h]); // T2aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            t2bbboff[h] = offset; offset += (trip_aab[h]+trip_aba[h])*(trip_aab[h]+trip_aba[h]); // T2bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            t2aaboff[h] = offset; offset += trip_aab[h]*trip_aab[h]; // T2aab
        }
        for (int h = 0; h < nirrep_; h++) {
            t2bbaoff[h] = offset; offset += trip_aab[h]*trip_aab[h]; // T2bba
        }
    }
    if ( constrain_d3_ ) {
        d3aaaoff = (int*)malloc(nirrep_*sizeof(int));
        d3bbboff = (int*)malloc(nirrep_*sizeof(int));
        d3aaboff = (int*)malloc(nirrep_*sizeof(int));
        d3bbaoff = (int*)malloc(nirrep_*sizeof(int));
        for (int h = 0; h < nirrep_; h++) {
            d3aaaoff[h] = offset; offset += trip_aaa[h]*trip_aaa[h]; // D3aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            d3bbboff[h] = offset; offset += trip_aaa[h]*trip_aaa[h]; // D3bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            d3aaboff[h] = offset; offset += trip_aab[h]*trip_aab[h]; // D3aab
        }
        for (int h = 0; h < nirrep_; h++) {
            d3bbaoff[h] = offset; offset += trip_aab[h]*trip_aab[h]; // D3bba
        }
    }
    if ( constrain_d4_ ) {
        d4aaaaoff = (int*)malloc(nirrep_*sizeof(int));
        d4aaaboff = (int*)malloc(nirrep_*sizeof(int));
        d4aabboff = (int*)malloc(nirrep_*sizeof(int));
        d4bbbaoff = (int*)malloc(nirrep_*sizeof(int));
        d4bbbboff = (int*)malloc(nirrep_*sizeof(int));
        for (int h = 0; h < nirrep_; h++) {
            d4aaaaoff[h] = offset; offset += quartet_aaaa[h]*quartet_aaaa[h]; // D4aaaa
        }
        for (int h = 0; h < nirrep_; h++) {
            d4aaaboff[h] = offset; offset += quartet_aaab[h]*quartet_aaab[h]; // D4aaab
        }
        for (int h = 0; h < nirrep_; h++) {
            d4aabboff[h] = offset; offset += quartet_aabb[h]*quartet_aabb[h]; // D4aabb
        }
        for (int h = 0; h < nirrep_; h++) {
            d4bbbaoff[h] = offset; offset += quartet_aaab[h]*quartet_aaab[h]; // D4bbba
        }
        for (int h = 0; h < nirrep_; h++) {
            d4bbbboff[h] = offset; offset += quartet_aaaa[h]*quartet_aaaa[h]; // D4bbbb
        }
    }

    // constraints:
    nconstraints_ = 0;

    constrain_sz_ = options_.get_bool("CONSTRAIN_SZ");

    if ( constrain_sz_ ) {
        nconstraints_ += 1;                   // Tr(D2ab)
        nconstraints_ += 1;                   // Tr(D2aa)
        nconstraints_ += 1;                   // Tr(D2bb)
    }else{
        nconstraints_ += 1;                   // Tr(D2ab + D2aa + D2bb)
    }

    //for ( int h = 0; h < nirrep_; h++) {
    //    nconstraints_ += gems_ab[h]*gems_ab[h]; // D2ab hermiticity
    //    nconstraints_ += gems_aa[h]*gems_aa[h]; // D2aa hermiticity
    //    nconstraints_ += gems_aa[h]*gems_aa[h]; // D2bb hermiticity
    //}

    for ( int h = 0; h < nirrep_; h++) {
        nconstraints_ += amopi_[h]*amopi_[h]; // D1a <-> Q1a
    }
    for ( int h = 0; h < nirrep_; h++) {
        nconstraints_ += amopi_[h]*amopi_[h]; // D1b <-> Q1b
    }
    if ( constrain_sz_ ) {
        for ( int h = 0; h < nirrep_; h++) {
            nconstraints_ += amopi_[h]*amopi_[h]; // contract D2ab        -> D1 a
        }
        for ( int h = 0; h < nirrep_; h++) {
            nconstraints_ += amopi_[h]*amopi_[h]; // contract D2ab        -> D1 b
        }
        for ( int h = 0; h < nirrep_; h++) {
            nconstraints_ += amopi_[h]*amopi_[h]; // contract D2aa        -> D1 a
        }
        for ( int h = 0; h < nirrep_; h++) {
            nconstraints_ += amopi_[h]*amopi_[h]; // contract D2bb        -> D1 b
        }
    }else {
        for ( int h = 0; h < nirrep_; h++) {
            nconstraints_ += amopi_[h]*amopi_[h]; // contract D2aa + D2ab        -> D1 a
        }
        for ( int h = 0; h < nirrep_; h++) {
            nconstraints_ += amopi_[h]*amopi_[h]; // contract D2bb + D2ab        -> D1 b
        }
    }

    if ( constrain_spin_ ) {
        nconstraints_ += 1;               // spin
        // additional spin constraints for singlets:
        if ( nalpha_ == nbeta_ ) {
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += amopi_[h]*amopi_[h]; // D1a = D1b
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_aa[h]*gems_aa[h]; // D2aa = D2bb
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_aa[h]*gems_aa[h]; // D2aa[pq][rs] = 1/2(D2ab[pq][rs] - D2ab[pq][sr] - D2ab[qp][rs] + D2ab[qp][sr])
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_aa[h]*gems_aa[h]; // D2bb[pq][rs] = 1/2(D2ab[pq][rs] - D2ab[pq][sr] - D2ab[qp][rs] + D2ab[qp][sr])
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_ab[h]*gems_ab[h];  // D200[pq][rs] = 1/(2 sqrt(1+dpq)sqrt(1+drs))(D2ab[pq][rs] + D2ab[pq][sr] + D2ab[qp][rs] + D2ab[qp][sr])
            }
        }else { // nonsinglets
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += 4*gems_ab[h]*gems_ab[h]; // D200_0, D210_0, D201_0, D211_0
            }
        }
        // maximal spin constraints:
        if ( constrain_g2_ ) {
            nconstraints_ += gems_ab[0];
            nconstraints_ += gems_ab[0];
        }
    }

    if ( constrain_q2_ ) {
        if ( ! spin_adapt_q2_ ) {
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_ab[h]*gems_ab[h]; // Q2ab
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_aa[h]*gems_aa[h]; // Q2aa
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_aa[h]*gems_aa[h]; // Q2bb
            }
        }else {
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_00[h]*gems_00[h]; // Q2s
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_aa[h]*gems_aa[h]; // Q2t
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_aa[h]*gems_aa[h]; // Q2t_p1
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_aa[h]*gems_aa[h]; // Q2t_m1
            }
        }

    }
    if ( constrain_g2_ ) {
        if ( ! spin_adapt_g2_ ) {
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_ab[h]*gems_ab[h]; // G2ab
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_ab[h]*gems_ab[h]; // G2ba
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += 2*gems_ab[h]*2*gems_ab[h]; // G2aa
            }
        }else {
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_ab[h]*gems_ab[h]; // G2s
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_ab[h]*gems_ab[h]; // G2t
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_ab[h]*gems_ab[h]; // G2t_p1
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_ab[h]*gems_ab[h]; // G2t_m1
            }
        }
    }
    if ( constrain_t1_ ) {
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += trip_aaa[h]*trip_aaa[h]; // T1aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += trip_aaa[h]*trip_aaa[h]; // T1bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += trip_aab[h]*trip_aab[h]; // T1aab
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += trip_aab[h]*trip_aab[h]; // T1bba
        }
    }
    if ( constrain_t2_ ) {
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += (trip_aab[h]+trip_aba[h])*(trip_aab[h]+trip_aba[h]); // T2aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += (trip_aab[h]+trip_aba[h])*(trip_aab[h]+trip_aba[h]); // T2bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += trip_aab[h]*trip_aab[h]; // T2aab
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += trip_aab[h]*trip_aab[h]; // T2bba
        }
    }
    if ( constrain_d3_ ) {
        if ( nalpha_ - nrstc_ - nfrzc_ > 2 ) {
            for (int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_aa[h]*gems_aa[h]; // D3aaa -> D2aa
            }
        }
        if ( nbeta_ - nrstc_ - nfrzc_ > 2 ) {
            for (int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_aa[h]*gems_aa[h]; // D3bbb -> D2bb
            }
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += gems_aa[h]*gems_aa[h]; // D3aab -> D2aa
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += gems_aa[h]*gems_aa[h]; // D3bba -> D2bb
        }
        if ( nalpha_ - nrstc_ - nfrzc_ > 1 ) {
            for (int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_ab[h]*gems_ab[h]; // D3aab -> D2ab
            }
        }
        if ( nbeta_ - nrstc_ - nfrzc_ > 1 ) {
            for (int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_ab[h]*gems_ab[h]; // D3bba -> D2ab
            }
        }
        // additional spin constraints for singlets:
        if ( constrain_spin_ && nalpha_ == nbeta_ ) {
            for (int h = 0; h < nirrep_; h++) {
                nconstraints_ += trip_aab[h]*trip_aab[h]; // D3aab = D3bba
            }
            for (int h = 0; h < nirrep_; h++) {
                nconstraints_ += trip_aaa[h]*trip_aaa[h]; // D3aab -> D3aaa
                nconstraints_ += trip_aaa[h]*trip_aaa[h]; // D3bba -> D3bbb
            }
        }
    }
    if ( constrain_d4_ ) {
        if ( nalpha_ - nrstc_ - nfrzc_ > 3 ) {
            for (int h = 0; h < nirrep_; h++) {
                nconstraints_ += trip_aaa[h]*trip_aaa[h]; // D4aaaa -> D3aaa
            }
        }
        if ( nalpha_ - nrstc_ - nfrzc_ > 2 && nbeta_ - nrstc_ - nfrzc_ > 0 ) {
            for (int h = 0; h < nirrep_; h++) {
                nconstraints_ += trip_aaa[h]*trip_aaa[h]; // D4aaab -> D3aaa
            }
            for (int h = 0; h < nirrep_; h++) {
                nconstraints_ += trip_aab[h]*trip_aab[h]; // D4aaab -> D3aab
            }
        }
        if ( nalpha_ - nrstc_ - nfrzc_ > 1  && nbeta_ - nrstc_ - nfrzc_ > 1) {
            for (int h = 0; h < nirrep_; h++) {
                nconstraints_ += trip_aab[h]*trip_aab[h]; // D4aabb -> D3aab
            }
            for (int h = 0; h < nirrep_; h++) {
                nconstraints_ += trip_aab[h]*trip_aab[h]; // D4aabb -> D3abb
            }
        }
        if ( nalpha_ - nrstc_ - nfrzc_ > 1  && nbeta_ - nrstc_ - nfrzc_ > 1) {
            for (int h = 0; h < nirrep_; h++) {
                nconstraints_ += trip_aab[h]*trip_aab[h]; // D4aabb -> D3aab
            }
            for (int h = 0; h < nirrep_; h++) {
                nconstraints_ += trip_aab[h]*trip_aab[h]; // D4aabb -> D3abb
            }
        }
        if ( nalpha_ - nrstc_ - nfrzc_ > 0 && nbeta_ - nrstc_ - nfrzc_ > 2 ) {
            for (int h = 0; h < nirrep_; h++) {
                nconstraints_ += trip_aab[h]*trip_aab[h]; // D4abbb -> D3abb
            }
            for (int h = 0; h < nirrep_; h++) {
                nconstraints_ += trip_aaa[h]*trip_aaa[h]; // D4abbb -> D3bbb
            }
        }
        if ( nbeta_ - nrstc_ - nfrzc_ > 3 ) {
            for (int h = 0; h < nirrep_; h++) {
                nconstraints_ += trip_aaa[h]*trip_aaa[h]; // D4bbbb -> D3bbb
            }
        }
    }

    // list of dimensions_
    for (int h = 0; h < nirrep_; h++) {
        dimensions_.push_back(gems_ab[h]); // D2ab
    }
    for (int h = 0; h < nirrep_; h++) {
        dimensions_.push_back(gems_aa[h]); // D2aa
    }
    for (int h = 0; h < nirrep_; h++) {
        dimensions_.push_back(gems_aa[h]); // D2bb
    }
    if ( constrain_spin_ && nalpha_ == nbeta_ ) {
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(gems_ab[h]); // D200
        }
    }else if ( constrain_spin_ ) {
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(2*gems_ab[h]); // D200_0,D210_0,D201_0,D211_0
        }
    }
    for (int h = 0; h < nirrep_; h++) {
        dimensions_.push_back(amopi_[h]); // D1a
    }
    for (int h = 0; h < nirrep_; h++) {
        dimensions_.push_back(amopi_[h]); // D1b
    }
    for (int h = 0; h < nirrep_; h++) {
        dimensions_.push_back(amopi_[h]); // Q1a
    }
    for (int h = 0; h < nirrep_; h++) {
        dimensions_.push_back(amopi_[h]); // Q1b
    }
    if ( constrain_q2_ ) {
        if ( !spin_adapt_q2_ ) {
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_ab[h]); // Q2ab
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_aa[h]); // Q2aa
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_aa[h]); // Q2bb
            }
        }else {
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_00[h]); // Q2s
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_aa[h]); // Q2t
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_aa[h]); // Q2t_p1
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_aa[h]); // Q2t_m1
            }
        }
    }
    if ( constrain_g2_ ) {
        if ( !spin_adapt_g2_ ) {
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_ab[h]); // G2ab
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_ab[h]); // G2ba
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(2*gems_ab[h]); // G2aa
            }
        }else {
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_ab[h]); // G2s
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_ab[h]); // G2t
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_ab[h]); // G2t_p1
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_ab[h]); // G2t_m1
            }
        }
    }
    if ( constrain_t1_ ) {
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aaa[h]); // T1aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aaa[h]); // T1bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aab[h]); // T1aab
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aab[h]); // T1bba
        }
    }
    if ( constrain_t2_ ) {
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aab[h]+trip_aba[h]); // T2aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aab[h]+trip_aba[h]); // T2bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aab[h]); // T2aab
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aab[h]); // T2bba
        }
    }
    if ( constrain_d3_ ) {
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aaa[h]); // D3aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aaa[h]); // D3bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aab[h]); // D3aab
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aab[h]); // D3bba
        }
    }
    if ( constrain_d4_ ) {
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(quartet_aaaa[h]); // D4aaaa
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(quartet_aaab[h]); // D4aaab
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(quartet_aabb[h]); // D4aabb
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(quartet_aaab[h]); // D4bbba
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(quartet_aaaa[h]); // D4bbbb
        }
    }

    // memory check happens here

    outfile->Printf("\n\n");
    outfile->Printf( "        ****************************************************\n");
    outfile->Printf( "        *                                                  *\n");
    outfile->Printf( "        *    v2RDM-CASSCF                                  *\n");
    outfile->Printf( "        *                                                  *\n");
    outfile->Printf( "        *    A variational 2-RDM-driven approach to the    *\n");
    outfile->Printf( "        *    active space self-consistent field method     *\n");
    outfile->Printf( "        *                                                  *\n");
    outfile->Printf( "        ****************************************************\n");

    outfile->Printf("\n");
    outfile->Printf("\n");
    outfile->Printf("        The following papers should be cited when using v2RDM-CASSCF:\n");
    outfile->Printf("\n");
    outfile->Printf("        J. Fosso-Tande, D. R. Nascimento, and A. E. DePrince III,\n");
    outfile->Printf("        Mol. Phys. 114, 423-430 (2015).\n");
    outfile->Printf("\n");
    outfile->Printf("            URL: http://dx.doi.org/10.1080/00268976.2015.1078008\n");
    outfile->Printf("\n");
    outfile->Printf("        J. Fosso-Tande, T.-S. Nguyen, G. Gidofalvi, and\n");
    outfile->Printf("        A. E. DePrince III, J. Chem. Theory Comput. 12, 2260-2271 (2016).\n");
    outfile->Printf("\n");
    outfile->Printf("            URL: http://dx.doi.org/10.1021/acs.jctc.6b00190\n");
    outfile->Printf("\n");
    outfile->Printf("\n");

    outfile->Printf("\n");
    outfile->Printf("  ==> Convergence parameters <==\n");
    outfile->Printf("\n");
    outfile->Printf("        r_convergence:                      %5.3le\n",options_.get_double("R_CONVERGENCE"));
    outfile->Printf("        e_convergence:                      %5.3le\n",options_.get_double("E_CONVERGENCE"));
    outfile->Printf("        cg_convergence:                     %5.3le\n",options_.get_double("CG_CONVERGENCE"));
    outfile->Printf("        maxiter:                             %8i\n",options_.get_int("MAXITER"));
    outfile->Printf("        cg_maxiter:                          %8i\n",options_.get_int("CG_MAXITER"));
    outfile->Printf("\n");

    // print orbitals per irrep in each space
    outfile->Printf("  ==> Active space details <==\n");
    outfile->Printf("\n");
    outfile->Printf("        Number of frozen core orbitals:         %5i\n",nfrzc_);
    outfile->Printf("        Number of restricted occupied orbitals: %5i\n",nrstc_);
    outfile->Printf("        Number of active occupied orbitals:     %5i\n",ndoccact);
    outfile->Printf("        Number of active virtual orbitals:      %5i\n",nvirt);
    outfile->Printf("        Number of restricted virtual orbitals:  %5i\n",nrstv_);
    outfile->Printf("        Number of frozen virtual orbitals:      %5i\n",nfrzv_);
    outfile->Printf("\n");

    std::vector<std::string> labels = reference_wavefunction_->molecule()->irrep_labels();
    outfile->Printf("        Irrep:           ");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("%4s",labels[h].c_str());
        if ( h < nirrep_ - 1 ) {
            outfile->Printf(",");
        }
    }
    outfile->Printf(" \n");
    outfile->Printf(" \n");

    outfile->Printf("        frozen_docc     [");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("%4i",frzcpi_[h]);
        if ( h < nirrep_ - 1 ) {
            outfile->Printf(",");
        }
    }
    outfile->Printf(" ]\n");
    outfile->Printf("        restricted_docc [");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("%4i",rstcpi_[h]);
        if ( h < nirrep_ - 1 ) {
            outfile->Printf(",");
        }
    }
    outfile->Printf(" ]\n");
    outfile->Printf("        active          [");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("%4i",amopi_[h]);
        if ( h < nirrep_ - 1 ) {
            outfile->Printf(",");
        }
    }
    outfile->Printf(" ]\n");
    outfile->Printf("        restricted_uocc [");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("%4i",rstvpi_[h]);
        if ( h < nirrep_ - 1 ) {
            outfile->Printf(",");
        }
    }
    outfile->Printf(" ]\n");
    outfile->Printf("        frozen_uocc     [");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("%4i",frzvpi_[h]);
        if ( h < nirrep_ - 1 ) {
            outfile->Printf(",");
        }
    }
    outfile->Printf(" ]\n");
    outfile->Printf("\n");

    outfile->Printf("  ==> Orbital optimization parameters <==\n");
    outfile->Printf("\n");
    outfile->Printf("        1-step algorithm:                   %5s\n",options_.get_bool("ORBOPT_ONE_STEP") ? "true" : "false");
    outfile->Printf("        g_convergence:                      %5.3le\n",options_.get_double("ORBOPT_GRADIENT_CONVERGENCE"));
    outfile->Printf("        e_convergence:                      %5.3le\n",options_.get_double("ORBOPT_ENERGY_CONVERGENCE"));
    outfile->Printf("        maximum iterations:                 %5i\n",options_.get_int("ORBOPT_MAXITER"));
    outfile->Printf("        frequency:                          %5i\n",options_.get_int("ORBOPT_FREQUENCY"));
    outfile->Printf("        active-active rotations:            %5s\n",options_.get_bool("ORBOPT_ACTIVE_ACTIVE_ROTATIONS") ? "true" : "false");
    outfile->Printf("        exact diagonal Hessian:             %5s\n",options_.get_bool("ORBOPT_EXACT_DIAGONAL_HESSIAN") ? "true" : "false");
    outfile->Printf("        number of DIIS vectors:             %5i\n",options_.get_int("ORBOPT_NUM_DIIS_VECTORS"));
    outfile->Printf("        print iteration info:               %5s\n",options_.get_bool("ORBOPT_WRITE") ? "true" : "false");

    outfile->Printf("\n");
    outfile->Printf("  ==> Memory requirements <==\n");
    outfile->Printf("\n");
    int nd2   = 0;
    int ng2    = 0;
    int nt1    = 0;
    int nt2    = 0;
    int maxgem = 0;
    for (int h = 0; h < nirrep_; h++) {
        nd2 +=     gems_ab[h]*gems_ab[h];
        nd2 += 2 * gems_aa[h]*gems_aa[h];

        ng2 +=     gems_ab[h] * gems_ab[h]; // G2ab
        ng2 +=     gems_ab[h] * gems_ab[h]; // G2ba
        ng2 += 4 * gems_ab[h] * gems_ab[h]; // G2aa

        if ( gems_ab[h] > maxgem ) {
            maxgem = gems_ab[h];
        }
        if ( constrain_g2_ ) {
            if ( 2*gems_ab[h] > maxgem ) {
                maxgem = 2*gems_ab[h];
            }
        }

        if ( constrain_t1_ ) {
            nt1 += trip_aaa[h] * trip_aaa[h]; // T1aaa
            nt1 += trip_aaa[h] * trip_aaa[h]; // T1bbb
            nt1 += trip_aab[h] * trip_aab[h]; // T1aab
            nt1 += trip_aab[h] * trip_aab[h]; // T1bba
            if ( trip_aab[h] > maxgem ) {
                maxgem = trip_aab[h];
            }
        }

        if ( constrain_t2_ ) {
            nt2 += (trip_aab[h]+trip_aba[h]) * (trip_aab[h]+trip_aba[h]); // T2aaa
            nt2 += (trip_aab[h]+trip_aba[h]) * (trip_aab[h]+trip_aba[h]); // T2bbb
            nt2 += trip_aab[h] * trip_aab[h]; // T2aab
            nt2 += trip_aab[h] * trip_aab[h]; // T2bba

            if ( trip_aab[h]+trip_aaa[h] > maxgem ) {
                maxgem = trip_aab[h]+trip_aaa[h];
            }
        }

    }

    outfile->Printf("        D2:                       %7.2lf mb\n",nd2 * 8.0 / 1024.0 / 1024.0);
    if ( constrain_q2_ ) {
        outfile->Printf("        Q2:                       %7.2lf mb\n",nd2 * 8.0 / 1024.0 / 1024.0);
    }
    if ( constrain_g2_ ) {
        outfile->Printf("        G2:                       %7.2lf mb\n",ng2 * 8.0 / 1024.0 / 1024.0);
    }
    if ( constrain_d3_ ) {
        outfile->Printf("        D3:                       %7.2lf mb\n",nt1 * 8.0 / 1024.0 / 1024.0);
    }
    if ( constrain_d4_ ) {
        outfile->Printf("        D4:                       %7.2lf mb\n",nt1 * 8.0 / 1024.0 / 1024.0);
    }
    if ( constrain_t1_ ) {
        outfile->Printf("        T1:                       %7.2lf mb\n",nt1 * 8.0 / 1024.0 / 1024.0);
    }
    if ( constrain_t2_ ) {
        outfile->Printf("        T2:                       %7.2lf mb\n",nt2 * 8.0 / 1024.0 / 1024.0);
    }
    outfile->Printf("\n");

    // we have 4 arrays the size of x and 4 the size of y
    // in addition, we need to store 3 times whatever the largest
    // block of x is for the diagonalization step
    // integrals:
    //     K2a, K2b
    // casscf:
    //     4-index integrals (no permutational symmetry)
    //     3-index integrals

    double tot = 4.0*dimx_ + 4.0*nconstraints_ + 3.0*maxgem*maxgem;
    tot += nd2; // for K2a, K2b

    // for casscf, need d2 and 3- or 4-index integrals

    // storage requirements for full d2
    for (int h = 0; h < nirrep_; h++) {
        tot += gems_plus_core[h] * ( gems_plus_core[h] + 1 ) / 2;
    }
    if ( is_df_ ) {
        // storage requirements for df integrals
        nQ_ = Process::environment.globals["NAUX (SCF)"];
        if ( options_.get_str("SCF_TYPE") == "DF" ) {
            std::shared_ptr<BasisSet> primary = reference_wavefunction_->basisset();
            std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");

            nQ_ = auxiliary->nbf();
            Process::environment.globals["NAUX (SCF)"] = nQ_;
        }
        tot += (long int)nQ_*(long int)nmo_*((long int)nmo_+1)/2;
    }else {
        // storage requirements for four-index integrals
        for (int h = 0; h < nirrep_; h++) {
            tot += (long int)gems_full[h] * ( (long int)gems_full[h] + 1L ) / 2L;
        }
    }

    // memory available after allocating all we need for v2RDM-CASSCF
    available_memory_ = memory_ - tot * 8L;

    outfile->Printf("        Total number of variables:     %10i\n",dimx_);
    outfile->Printf("        Total number of constraints:   %10i\n",nconstraints_);
    outfile->Printf("        Total memory requirements:     %7.2lf mb\n",tot * 8.0 / 1024.0 / 1024.0);
    outfile->Printf("\n");

    if ( tot * 8.0 > (double)memory_ ) {
        outfile->Printf("\n");
        outfile->Printf("        Not enough memory!\n");
        outfile->Printf("\n");
        if ( !is_df_ ) {
            outfile->Printf("        Either increase the available memory by %7.2lf mb\n",(8.0 * tot - memory_)/1024.0/1024.0);
            outfile->Printf("        or try scf_type = df or scf_type = cd\n");

        }else {
            outfile->Printf("        Increase the available memory by %7.2lf mb.\n",(8.0 * tot - memory_)/1024.0/1024.0);
        }
        outfile->Printf("\n");
        throw PsiException("Not enough memory",__FILE__,__LINE__);
    }

    // mo-mo transformation matrix
    newMO_ = (SharedMatrix)(new Matrix(reference_wavefunction_->Ca()));
    newMO_->zero();
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < nmopi_[h]; i++) {
            newMO_->pointer(h)[i][i] = 1.0;
        }
    }

    orbopt_transformation_matrix_ = (double*)malloc((nmo_-nfrzc_-nfrzv_)*(nmo_-nfrzc_-nfrzv_)*sizeof(double));
    memset((void*)orbopt_transformation_matrix_,'\0',(nmo_-nfrzc_-nfrzv_)*(nmo_-nfrzc_-nfrzv_)*sizeof(double));
    for (int i = 0; i < nmo_-nfrzc_-nfrzv_; i++) {
        orbopt_transformation_matrix_[i*(nmo_-nfrzc_-nfrzv_)+i] = 1.0;
    }

    //  if restarting, need to grab Ca_ from disk before integral transformation
    // checkpoint file
    if ( options_.get_str("RESTART_FROM_CHECKPOINT_FILE") != "" ) {
        ReadOrbitalsFromCheckpointFile();
    }

    // if using 3-index integrals, transform them before allocating any memory integrals, transform
    if ( is_df_ ) {
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

        outfile->Printf("\n");
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

    // allocate vectors
    x      = SharedVector(new Vector("primal solution",dimx_));
    c      = SharedVector(new Vector("OEI and TEI",dimx_));
    b      = SharedVector(new Vector("constraints",nconstraints_));

    // input/output array for orbopt sweeps

    int nthread = omp_get_max_threads();

    orbopt_data_    = (double*)malloc(15*sizeof(double));
    orbopt_data_[0] = (double)nthread;
    orbopt_data_[1] = (double)(options_.get_bool("ORBOPT_ACTIVE_ACTIVE_ROTATIONS") ? 1.0 : 0.0 );
    orbopt_data_[2] = (double)nfrzc_; //(double)options_.get_int("ORBOPT_FROZEN_CORE");
    orbopt_data_[3] = (double)options_.get_double("ORBOPT_GRADIENT_CONVERGENCE");
    orbopt_data_[4] = (double)options_.get_double("ORBOPT_ENERGY_CONVERGENCE");
    orbopt_data_[5] = (double)(options_.get_bool("ORBOPT_WRITE") ? 1.0 : 0.0 );
    orbopt_data_[6] = (double)(options_.get_bool("ORBOPT_EXACT_DIAGONAL_HESSIAN") ? 1.0 : 0.0 );
    orbopt_data_[7] = (double)options_.get_int("ORBOPT_NUM_DIIS_VECTORS");
    orbopt_data_[8] = (double)options_.get_int("ORBOPT_MAXITER");
    orbopt_data_[9] = 0.0;
    if ( is_df_ ) {
      orbopt_data_[9] = 1.0;
    }
    orbopt_data_[10] = 0.0;  // number of iterations (output)
    orbopt_data_[11] = 0.0;  // gradient norm (output)
    orbopt_data_[12] = 0.0;  // change in energy (output)
    orbopt_data_[13] = 0.0;  // converged?

    // orbital optimizatoin algorithm
    orbopt_data_[14] = 0.0;
    if      ( options_.get_str("ORBOPT_ALGORITHM") == "QUASI_NEWTON" )       orbopt_data_[14] = 0.0;
    else if ( options_.get_str("ORBOPT_ALGORITHM") == "CONJUGATE_GRADIENT" ) orbopt_data_[14] = 1.0;
    else if ( options_.get_str("ORBOPT_ALGORITHM") == "NEWTON_RAPHSON" )     orbopt_data_[14] = 2.0;

    orbopt_converged_ = false;

    // don't change the length of this filename
    orbopt_outfile_ = (char*)malloc(120*sizeof(char));
    std::string filename = get_writer_file_prefix(reference_wavefunction_->molecule()->name()) + ".orbopt";
    strcpy(orbopt_outfile_,filename.c_str());
    if ( options_.get_bool("ORBOPT_WRITE") ) {
        FILE * fp = fopen(orbopt_outfile_,"w");
        fclose(fp);
    }

    // initialize timers and iteration counters
    orbopt_iter_total_ = 0;
    orbopt_time_       = 0.0;

    // allocate memory for orbital lagrangian (TODO: make these smaller)
    X_               = (double*)malloc(nmo_*nmo_*sizeof(double));

    // even if we use rhf/rohf reference, we need same_a_b_orbs_=false
    // to trigger the correct integral transformations in deriv.cc
    same_a_b_orbs_ = false;
    same_a_b_dens_ = false;

    // sdp solver
    sdp_ = (std::shared_ptr<BPSDPSolver>)(new BPSDPSolver(dimx_,nconstraints_,options_));
}

int v2RDMSolver::SymmetryPair(int i,int j) {
    return table[i*8+j];
}
int v2RDMSolver::TotalSym(int i,int j,int k, int l) {
    return SymmetryPair(SymmetryPair(symmetry[i],symmetry[j]),SymmetryPair(symmetry[k],symmetry[l]));
}

// compute the energy!
double v2RDMSolver::compute_energy() {

    double start_total_time = omp_get_wtime();

    // print guess orbitals in molden format
    if ( options_.get_bool("GUESS_ORBITALS_WRITE") ) {

        std::shared_ptr<MoldenWriter> molden(new MoldenWriter(reference_wavefunction_));
        std::shared_ptr<Vector> zero (new Vector("",nirrep_,nmopi_));
        zero->zero();
        std::string filename = get_writer_file_prefix(reference_wavefunction_->molecule()->name()) + ".guess.molden";

        Ca_ = SharedMatrix(reference_wavefunction_->Ca());
        Cb_ = SharedMatrix(reference_wavefunction_->Cb());

        SharedVector occupation_a= SharedVector(new Vector(nirrep_, nmopi_));
		SharedVector occupation_b= SharedVector(new Vector(nirrep_, nmopi_));

		for (int h = 0; h < nirrep_; h++) {
			for (int i = 0; i < nalphapi_[h]; i++) {
				occupation_a->set(h, i, 1.0);
			}
			for (int i = 0; i < nbetapi_[h]; i++) {
				occupation_b->set(h, i, 1.0);
			}
		}

        molden->write(filename, Ca_, Cb_, zero, zero, occupation_a, occupation_b, true);

    }

    // hartree-fock guess
    Guess();

    // checkpoint file
    if ( options_.get_str("RESTART_FROM_CHECKPOINT_FILE") != "" ) {
        ReadFromCheckpointFile();
    }

    // get integrals
    GetIntegrals();

    // generate constraint vector
    BuildConstraints();

    // iterate
    int orbopt_iter = 0;
    do {

        sdp_->solve(x, b, c, dimensions_, options_.get_int("ORBOPT_FREQUENCY"), evaluate_Au, evaluate_ATu, evaluate_cg_lhs, (void*)this);

        if ( options_.get_bool("OPTIMIZE_ORBITALS") ) {
    
            double start = omp_get_wtime();
            RotateOrbitals();
            double end = omp_get_wtime();
    
            orbopt_time_  += end - start;
            orbopt_iter_total_++;

            double energy_primal = C_DDOT(dimx_,c->pointer(),1,x->pointer(),1);
            outfile->Printf("            Total energy: %20.12lf\n",energy_primal+enuc_+efzc_);

        }else {
            orbopt_converged_ = true;
        }

        outfile->Printf("\n");

    }while( !orbopt_converged_ || !sdp_->is_converged() );

    outfile->Printf("\n");
    outfile->Printf("      v2RDM iterations converged!\n");
    outfile->Printf("\n");

    // evaluate spin squared
    double s2 = 0.0;
    double * x_p = x->pointer();
    for (int i = 0; i < amo_; i++){
        for (int j = 0; j < amo_; j++){
            int h = SymmetryPair(symmetry[i],symmetry[j]);
            int ij = ibas_ab_sym[h][i][j];
            int ji = ibas_ab_sym[h][j][i];
            s2 += x_p[d2aboff[h] + ij*gems_ab[h]+ji];
        }
    }
    double na = nalpha_ - nfrzc_ - nrstc_;
    double nb = nbeta_ - nfrzc_ - nrstc_;
    double ms = (multiplicity_ - 1.0)/2.0;

    // set energy for wavefunction
    double energy_primal = C_DDOT(dimx_,c->pointer(),1,x->pointer(),1);
    energy_ = energy_primal+enuc_+efzc_;

    // push final transformation matrix onto Ca_ and Cb_
    UpdateTransformationMatrix();

    energy_primal = C_DDOT(dimx_,c->pointer(),1,x->pointer(),1);

    // break down energy into one- and two-electron components

    double kinetic, potential, two_electron_energy;
    EnergyByComponent(kinetic,potential,two_electron_energy);

    outfile->Printf("      na:                        %20.6lf\n", na);
    outfile->Printf("      nb:                        %20.6lf\n", nb);
    outfile->Printf("      v2RDM total spin [S(S+1)]: %20.6lf\n", 0.5 * (na + nb) + ms*ms - s2);
    outfile->Printf("\n");
    outfile->Printf("      Nuclear Repulsion Energy:          %20.12lf\n",enuc_);
    outfile->Printf("      Two-Electron Energy:               %20.12lf\n",two_electron_energy);
    //outfile->Printf("      Two-Electron (active):             %20.12lf\n",active_two_electron_energy);
    outfile->Printf("      Kinetic Energy:                    %20.12lf\n",kinetic);
    outfile->Printf("      Electron-Nuclear Potential Energy: %20.12lf\n",potential);

    outfile->Printf("\n");
    outfile->Printf("    * v2RDM total energy:                %20.12lf\n",energy_primal+enuc_+efzc_);
    outfile->Printf("\n");

    Process::environment.globals["CURRENT ENERGY"]     = energy_primal+enuc_+efzc_;
    Process::environment.globals["v2RDM TOTAL ENERGY"] = energy_primal+enuc_+efzc_;

    if ( options_.get_bool("SEMICANONICALIZE_ORBITALS") ) {
        if ( options_.get_str("DERTYPE") == "FIRST" ) {
            outfile->Printf("\n");
            outfile->Printf("    <<< WARNING >>> skipping orbital semicanonicalization (blame gradients)\n");
            outfile->Printf("\n");
        }else {

            //printf("primal energy before semicanonicalization: %20.12lf\n",energy_primal+efzc_);
            orbopt_data_[8] = -2.0;
            RotateOrbitals();

            // push final transformation matrix onto Ca_ and Cb_
            UpdateTransformationMatrix();

        }
    } 

    // compute and natural orbitals and transform 1-RDM/2-RDM to the NO basis
    if ( options_.get_bool("NAT_ORBS") || options_.get_bool("FCIDUMP") || options_.get_bool("EXTENDED_KOOPMANS") ) {
        ComputeNaturalOrbitals();
    }
    if ( options_.get_bool("MOLDEN_WRITE") ) {
        WriteMoldenFile();
    }
    if ( options_.get_bool("EXTENDED_KOOPMANS") ) {
        ExtendedKoopmans();
    }
    if ( options_.get_bool("FCIDUMP") ) {
        FCIDUMP();
    }

    // compute and print natural orbital occupation numbers 
    PrintNaturalOrbitalOccupations();

    // push OPDM onto wavefunction object
    FinalizeOPDM();

    // write tpdm to disk?
    if ( options_.get_bool("TPDM_WRITE") ) {
        WriteActiveTPDM();
    }
    if ( options_.get_bool("TPDM_WRITE_FULL") ) {
        WriteTPDM();
        //ReadTPDM();
    }
    if ( options_.get_bool("TPDM_WRITE_SPIN_FREE") ) {
        WriteTPDMSpinFree();
        //ReadTPDM();
    }
    if ( options_.get_bool("OPDM_WRITE_FULL") ) {
        WriteOPDM();
    }
    // write 3-particle density matrix to disk?
    if ( options_.get_bool("3PDM_WRITE") && options_.get_bool("CONSTRAIN_D3")) {
        WriteActive3PDM();
        //Read3PDM();
    }

    // for derivatives:
    if ( options_.get_str("DERTYPE") == "FIRST" ) {

        if ( options_.get_bool("NAT_ORBS") || options_.get_bool("FCIDUMP") || options_.get_bool("EXTENDED_KOOPMANS") ) {
            throw PsiException("analytic gradients require nat_orbs false",__FILE__,__LINE__);
        }

        // write checkpoint file for next step in optimization
        WriteCheckpointFile();

        orbopt_data_[8] = -1.0;
        RotateOrbitals();

        // write 2-RDM in IWL format
        WriteTPDM_IWL();

        // push orbital lagrangian onto wave function
        OrbitalLagrangian();

    }else if ( options_.get_bool("WRITE_CHECKPOINT_FILE") ) {

        WriteCheckpointFile();

    }


    double end_total_time = omp_get_wtime();

    outfile->Printf("\n");
    outfile->Printf("  ==> Iteration count <==\n");
    outfile->Printf("\n");
    outfile->Printf("      Microiterations:            %12li\n",sdp_->iiter_total());
    outfile->Printf("      Macroiterations:            %12li\n",sdp_->oiter_total());
    outfile->Printf("      Orbital optimization steps: %12li\n",orbopt_iter_total_);
    outfile->Printf("\n");
    outfile->Printf("  ==> Wall time <==\n");
    outfile->Printf("\n");
    outfile->Printf("      Microiterations:            %12.2lf s\n",sdp_->iiter_time());
    outfile->Printf("      Macroiterations:            %12.2lf s\n",sdp_->oiter_time());
    outfile->Printf("      Orbital optimization:       %12.2lf s\n",orbopt_time_);
    outfile->Printf("      Total:                      %12.2lf s\n",end_total_time - start_total_time);
    outfile->Printf("\n");

    //CheckSpinStructure();

    return energy_primal + enuc_ + efzc_;
}

void v2RDMSolver::EnergyByComponent(double &kinetic, double &potential, double &two_electron_energy) {

    double * x_p = x->pointer();
    double * c_p = c->pointer();

    std::shared_ptr<MintsHelper> mints(new MintsHelper(reference_wavefunction_));
    std::shared_ptr<Matrix> myT (new Matrix(mints->so_kinetic()));
    std::shared_ptr<Matrix> myV (new Matrix(mints->so_potential()));

    myT->transform(Ca_);
    myV->transform(Ca_);
    kinetic = 0.0;
    potential = 0.0;

    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            int ii = i + rstcpi_[h] + frzcpi_[h];
            for (int j = 0; j < amopi_[h]; j++) {
                int jj = j + rstcpi_[h] + frzcpi_[h];

                if ( ii == jj ) continue;
                
                kinetic   += myT->pointer(h)[ii][jj] * x_p[d1aoff[h] + i * amopi_[h] + j];
                kinetic   += myT->pointer(h)[ii][jj] * x_p[d1boff[h] + i * amopi_[h] + j];
                
                potential += myV->pointer(h)[ii][jj] * x_p[d1aoff[h] + i * amopi_[h] + j];
                potential += myV->pointer(h)[ii][jj] * x_p[d1boff[h] + i * amopi_[h] + j];
            
            }

            kinetic   += myT->pointer(h)[ii][ii] * x_p[d1aoff[h] + i * amopi_[h] + i];
            kinetic   += myT->pointer(h)[ii][ii] * x_p[d1boff[h] + i * amopi_[h] + i];
            
            potential += myV->pointer(h)[ii][ii] * x_p[d1aoff[h] + i * amopi_[h] + i];
            potential += myV->pointer(h)[ii][ii] * x_p[d1boff[h] + i * amopi_[h] + i];
            
        }
    }

    two_electron_energy = efzc_;
    // remove one-electron parts of frozen core energy
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < rstcpi_[h] + frzcpi_[h]; i++) {
            two_electron_energy -= 2.0 * myT->pointer(h)[i][i];
            two_electron_energy -= 2.0 * myV->pointer(h)[i][i];
        }
    }
    double active_two_electron_energy = 0.0;
    for (int h = 0; h < nirrep_; h++) {
        active_two_electron_energy += C_DDOT(gems_ab[h] * gems_ab[h], x_p + d2aboff[h],1, c->pointer() + d2aboff[h], 1);
        active_two_electron_energy += C_DDOT(gems_aa[h] * gems_aa[h], x_p + d2aaoff[h],1, c->pointer() + d2aaoff[h], 1);
        active_two_electron_energy += C_DDOT(gems_aa[h] * gems_aa[h], x_p + d2bboff[h],1, c->pointer() + d2bboff[h], 1);
    }
    two_electron_energy += active_two_electron_energy;

    // core-active part of two-electron energy
    double core_active = 0.0;
    for (int h = 0; h < nirrep_; h++) {
        core_active += C_DDOT(amopi_[h] * amopi_[h], x_p + d1aoff[h],1, c->pointer() + d1aoff[h], 1);
        core_active += C_DDOT(amopi_[h] * amopi_[h], x_p + d1boff[h],1, c->pointer() + d1boff[h], 1);
    }
    core_active -= kinetic;
    core_active -= potential;

    two_electron_energy += core_active;

    // lastly, don't forget core contribution to kinetic and potential energy
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < rstcpi_[h] + frzcpi_[h]; i++) {
            kinetic   += 2.0 * myT->pointer(h)[i][i];
            potential += 2.0 * myV->pointer(h)[i][i];
        }
    }
}

void v2RDMSolver::CheckSpinStructure() {
    double * x_p = x->pointer();
    // D1a = D1b
    for ( int h = 0; h < nirrep_; h++) {
        C_DAXPY(amopi_[h]*amopi_[h],-1.0,x_p + d1boff[h],1,x_p + d1aoff[h],1);
        printf("d1 %20.12lf\n",C_DNRM2(amopi_[h]*amopi_[h],x_p + d1aoff[h],1));
    }
    // D2aa[pq][rs] = 1/2(D2ab[pq][rs] - D2ab[pq][sr])
    for ( int h = 0; h < nirrep_; h++) {
        double tot = 0.0;
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            int ijb = ibas_ab_sym[h][i][j];
            int jib = ibas_ab_sym[h][j][i];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                int klb = ibas_ab_sym[h][k][l];
                int lkb = ibas_ab_sym[h][l][k];
                double dum = x_p[d2aaoff[h] + ij*gems_aa[h]+kl];
                dum -= 0.5 * x_p[d2aboff[h] + ijb*gems_ab[h] + klb];
                dum += 0.5 * x_p[d2aboff[h] + jib*gems_ab[h] + klb];
                dum += 0.5 * x_p[d2aboff[h] + ijb*gems_ab[h] + lkb];
                dum -= 0.5 * x_p[d2aboff[h] + jib*gems_ab[h] + lkb];
                tot += dum*dum;
            }
        }
        printf("d2ab -> d2aa %20.12lf\n",sqrt(tot));
    }
    // D2bb[pq][rs] = 1/2(D2ab[pq][rs] - D2ab[pq][sr])
    for ( int h = 0; h < nirrep_; h++) {
        double tot = 0.0;
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            int ijb = ibas_ab_sym[h][i][j];
            int jib = ibas_ab_sym[h][j][i];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                int klb = ibas_ab_sym[h][k][l];
                int lkb = ibas_ab_sym[h][l][k];
                double dum = x_p[d2aaoff[h] + ij*gems_aa[h]+kl];
                dum -= 0.5 * x_p[d2aboff[h] + ijb*gems_ab[h] + klb];
                dum += 0.5 * x_p[d2aboff[h] + jib*gems_ab[h] + klb];
                dum += 0.5 * x_p[d2aboff[h] + ijb*gems_ab[h] + lkb];
                dum -= 0.5 * x_p[d2aboff[h] + jib*gems_ab[h] + lkb];
                tot += dum*dum;
            }
        }
        printf("d2ab -> d2bb %20.12lf\n",sqrt(tot));
    }
    // D2aa = D2bb
    for ( int h = 0; h < nirrep_; h++) {
        C_DAXPY(gems_aa[h]*gems_aa[h],-1.0,x_p + d2bboff[h],1,x_p + d2aaoff[h],1);
        printf("d2aa=d2bb %20.12lf\n",C_DNRM2(gems_aa[h]*gems_aa[h],x_p + d2aaoff[h],1));
    }


}

void v2RDMSolver::WriteMoldenFile() {

    // it is possible the 1-RDM is already in the NO basis:
    if ( options_.get_bool("NAT_ORBS") || options_.get_bool("FCIDUMP") || options_.get_bool("EXTENDED_KOOPMANS") ) {

        std::shared_ptr<Vector> eigval (new Vector("Natural Orbital Occupation Numbers (spin free)",nirrep_,nmopi_));
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < frzcpi_[h] + rstcpi_[h]; i++) {
                eigval->pointer(h)[i] = 1.0;
            }
            for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h] - rstvpi_[h] - frzvpi_[h]; i++) {
                eigval->pointer(h)[i]  = 0.5 * x->pointer()[d1aoff[h]+(i-rstcpi_[h]-frzcpi_[h])*amopi_[h]+(i-rstcpi_[h]-frzcpi_[h])];
                eigval->pointer(h)[i] += 0.5 * x->pointer()[d1aoff[h]+(i-rstcpi_[h]-frzcpi_[h])*amopi_[h]+(i-rstcpi_[h]-frzcpi_[h])];
            }
        }
        // Print a molden file
        if ( options_["RESTART_FROM_CHECKPOINT_FILE"].has_changed() ) {
            throw PsiException("printing orbitals is currently disabled when restarting v2rdm jobs.  i can't remember why, though... sorry!",__FILE__,__LINE__);
        }
        //std::shared_ptr<MoldenWriter> molden(new MoldenWriter((std::shared_ptr<Wavefunction>)this));
        std::shared_ptr<MoldenWriter> molden(new MoldenWriter(reference_wavefunction_));
        std::shared_ptr<Vector> zero (new Vector("",nirrep_,nmopi_));
        zero->zero();
        std::string filename = get_writer_file_prefix(reference_wavefunction_->molecule()->name()) + ".molden";
        molden->write(filename,Ca_,Ca_,zero, zero,eigval,eigval,true);

    // otherwise, we need to compute the natural orbitals:
    }else {
        std::shared_ptr<Matrix> D (new Matrix(nirrep_,nmopi_,nmopi_));
        std::shared_ptr<Matrix> eigvec (new Matrix(nirrep_,nmopi_,nmopi_));
        std::shared_ptr<Vector> eigval (new Vector("Natural Orbital Occupation Numbers (spin free)",nirrep_,nmopi_));
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < frzcpi_[h] + rstcpi_[h]; i++) {
                D->pointer(h)[i][i] = 1.0;
            }
            for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h] - rstvpi_[h] - frzvpi_[h]; i++) {
                for (int j = rstcpi_[h] + frzcpi_[h]; j < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; j++) {
                    D->pointer(h)[i][j]  = 0.5 * x->pointer()[d1aoff[h]+(i-rstcpi_[h]-frzcpi_[h])*amopi_[h]+(j-rstcpi_[h]-frzcpi_[h])];
                    D->pointer(h)[i][j] += 0.5 * x->pointer()[d1boff[h]+(i-rstcpi_[h]-frzcpi_[h])*amopi_[h]+(j-rstcpi_[h]-frzcpi_[h])];
                }
            }
        }
        std::shared_ptr<Matrix> saved ( new Matrix(D) );
        D->diagonalize(eigvec,eigval,descending);
        //eigval->print();

        std::shared_ptr<Matrix> Cno (new Matrix(Ca_));

        //Ca_->print();
        // build AO/NO transformation matrix 
        for (int h = 0; h < nirrep_; h++) {
            for (int mu = 0; mu < nsopi_[h]; mu++) {
                double *  temp = (double*)malloc(nmopi_[h]*sizeof(double));
                double ** cp   = Cno->pointer(h);
                double ** ep   = eigvec->pointer(h);
                for (int i = 0; i < nmopi_[h]; i++) {
                    double dum = 0.0;
                    for (int j = 0; j < nmopi_[h]; j++) {
                        dum += cp[mu][j] * ep[j][i];
                    }
                    temp[i] = dum;
                }
                for (int i = 0; i < nmopi_[h]; i++) {
                    cp[mu][i] = temp[i];
                }
                free(temp);
            }
        }

        // Print a molden file
        if ( options_["RESTART_FROM_CHECKPOINT_FILE"].has_changed() ) {
            throw PsiException("printing orbitals is currently disabled when restarting v2rdm jobs.  i can't remember why, though... sorry!",__FILE__,__LINE__);
        }
        //std::shared_ptr<MoldenWriter> molden(new MoldenWriter((std::shared_ptr<Wavefunction>)this));
        std::shared_ptr<MoldenWriter> molden(new MoldenWriter(reference_wavefunction_));
        std::shared_ptr<Vector> zero (new Vector("",nirrep_,nmopi_));
        zero->zero();
        std::string filename = get_writer_file_prefix(reference_wavefunction_->molecule()->name()) + ".molden";
        molden->write(filename,Cno,Cno,zero, zero,eigval,eigval,true);
    }

}

void v2RDMSolver::FinalizeOPDM() {
    // nee

    std::stringstream ss;
    ss << "v-2RDM";
    std::stringstream ss_a;
    std::stringstream ss_b;
    ss_a << ss.str() << " alpha";
    ss_b << ss.str() << " beta";
    SharedMatrix opdm_a(new Matrix(ss_a.str(), Ca_->colspi(), Ca_->colspi()));
    SharedMatrix opdm_b(new Matrix(ss_b.str(), Ca_->colspi(), Ca_->colspi()));

    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < rstcpi_[h]+frzcpi_[h]; i++) {
            opdm_a->pointer(h)[i][i] = 1.0;
        }
        for (int i = rstcpi_[h]+frzcpi_[h]; i < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; i++) {
            for (int j = rstcpi_[h]+frzcpi_[h]; j < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; j++) {
                opdm_a->pointer(h)[i][j] = x->pointer()[d1aoff[h]+(i-rstcpi_[h]-frzcpi_[h])*amopi_[h]+(j-rstcpi_[h]-frzcpi_[h])];
            }
        }
    }

    int symm = opdm_a->symmetry();

    double* temp = (double*)malloc(Ca_->max_ncol() * Ca_->max_nrow() * sizeof(double));

    Da_->zero();
    for (int h = 0; h < nirrep_; h++) {
        int nmol = Ca_->colspi()[h];
        int nmor = Ca_->colspi()[h^symm];
        int nsol = Ca_->rowspi()[h];
        int nsor = Ca_->rowspi()[h^symm];
        if (!nmol || !nmor || !nsol || !nsor) continue;
        double** Clp = Ca_->pointer(h);
        double** Crp = Ca_->pointer(h^symm);
        double** Dmop = opdm_a->pointer(h^symm);
        double** Dsop = Da_->pointer(h^symm);
        C_DGEMM('N','T',nmol,nsor,nmor,1.0,Dmop[0],nmor,Crp[0],nmor,0.0,temp,nsor);
        C_DGEMM('N','N',nsol,nsor,nmol,1.0,Clp[0],nmol,temp,nsor,0.0,Dsop[0],nsor);
    }

    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < rstcpi_[h]+frzcpi_[h]; i++) {
            opdm_b->pointer(h)[i][i] = 1.0;
        }
        for (int i = rstcpi_[h]+frzcpi_[h]; i < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; i++) {
            for (int j = rstcpi_[h]+frzcpi_[h]; j < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; j++) {
                opdm_b->pointer(h)[i][j] = x->pointer()[d1boff[h]+(i-rstcpi_[h]-frzcpi_[h])*amopi_[h]+(j-rstcpi_[h]-frzcpi_[h])];
            }
        }
    }
    Db_->zero();
    // hmm... the IntegralTransform type is Restricted, so only use Ca_ here.
    for (int h = 0; h < nirrep_; h++) {
        int nmol = Ca_->colspi()[h];
        int nmor = Ca_->colspi()[h^symm];
        int nsol = Ca_->rowspi()[h];
        int nsor = Ca_->rowspi()[h^symm];
        if (!nmol || !nmor || !nsol || !nsor) continue;
        double** Clp = Ca_->pointer(h);
        double** Crp = Ca_->pointer(h^symm);
        double** Dmop = opdm_b->pointer(h^symm);
        double** Dsop = Db_->pointer(h^symm);
        C_DGEMM('N','T',nmol,nsor,nmor,1.0,Dmop[0],nmor,Crp[0],nmor,0.0,temp,nsor);
        C_DGEMM('N','N',nsol,nsor,nmol,1.0,Clp[0],nmol,temp,nsor,0.0,Dsop[0],nsor);
    }

    free(temp);
}

void v2RDMSolver::Guess(){

    double* x_p = x->pointer();
    memset((void*)x_p,'\0',dimx_*sizeof(double));

    if ( options_.get_str("TPDM_GUESS") == "HF" ) {

        // Hartree-Fock guess for D2, D1, Q1, Q2, and G2

        // D2ab
        int poff1 = 0;
        for (int h1 = 0; h1 < nirrep_; h1++) {
            for (int i = 0; i < soccpi_[h1] + doccpi_[h1] - rstcpi_[h1] - frzcpi_[h1]; i++){
                int poff2 = 0;
                for (int h2 = 0; h2 < nirrep_; h2++) {
                    for (int j = 0; j < doccpi_[h2] - rstcpi_[h2] - frzcpi_[h2]; j++){
                        int ii = i + poff1;
                        int jj = j + poff2;
                        int h3 = SymmetryPair(symmetry[ii],symmetry[jj]);
                        int ij = ibas_ab_sym[h3][ii][jj];
                        x_p[d2aboff[h3] + ij*gems_ab[h3]+ij] = 1.0;
                    }
                    poff2   += nmopi_[h2] - rstcpi_[h2] - frzcpi_[h2] - rstvpi_[h2] - frzvpi_[h2];
                }
            }
            poff1   += nmopi_[h1] - rstcpi_[h1] - frzcpi_[h1] - rstvpi_[h1] - frzvpi_[h1];
        }

        // d2aa
        poff1 = 0;
        for (int h1 = 0; h1 < nirrep_; h1++) {
            for (int i = 0; i < soccpi_[h1] + doccpi_[h1] - rstcpi_[h1] - frzcpi_[h1]; i++){
                int poff2 = 0;
                for (int h2 = 0; h2 < nirrep_; h2++) {
                    for (int j = 0; j < soccpi_[h2] + doccpi_[h2] - rstcpi_[h2] - frzcpi_[h2]; j++){
                        int ii = i + poff1;
                        int jj = j + poff2;
                        if ( jj >= ii ) continue;
                        int h3 = SymmetryPair(symmetry[ii],symmetry[jj]);
                        int ij = ibas_aa_sym[h3][ii][jj];
                        x_p[d2aaoff[h3] + ij*gems_aa[h3]+ij] = 1.0;
                    }
                    poff2   += nmopi_[h2] - rstcpi_[h2] - frzcpi_[h2] - rstvpi_[h2] - frzvpi_[h2];
                }
            }
            poff1   += nmopi_[h1] - rstcpi_[h1] - frzcpi_[h1] - rstvpi_[h1] - frzvpi_[h1];
        }

        // d2bb
        poff1 = 0;
        for (int h1 = 0; h1 < nirrep_; h1++) {
            for (int i = 0; i < doccpi_[h1] - rstcpi_[h1] - frzcpi_[h1]; i++){
                int poff2 = 0;
                for (int h2 = 0; h2 < nirrep_; h2++) {
                    for (int j = 0; j < doccpi_[h2] - rstcpi_[h2] - frzcpi_[h2]; j++){
                        int ii = i + poff1;
                        int jj = j + poff2;
                        if ( jj >= ii ) continue;
                        int h3 = SymmetryPair(symmetry[ii],symmetry[jj]);
                        int ij = ibas_aa_sym[h3][ii][jj];
                        x_p[d2bboff[h3] + ij*gems_aa[h3]+ij] = 1.0;
                    }
                    poff2   += nmopi_[h2] - rstcpi_[h2] - frzcpi_[h2] - rstvpi_[h2] - frzvpi_[h2];
                }
            }
            poff1   += nmopi_[h1] - rstcpi_[h1] - frzcpi_[h1] - rstvpi_[h1] - frzvpi_[h1];
        }

        // D1
        for (int h = 0; h < nirrep_; h++) {
            for (int i = rstcpi_[h] + frzcpi_[h]; i < doccpi_[h]+soccpi_[h]; i++) {
                int ii = i - rstcpi_[h] - frzcpi_[h];
                x_p[d1aoff[h]+ii*amopi_[h]+ii] = 1.0;
            }
            for (int i = rstcpi_[h] + frzcpi_[h]; i < doccpi_[h]; i++) {
                int ii = i - rstcpi_[h] - frzcpi_[h];
                x_p[d1boff[h]+ii*amopi_[h]+ii] = 1.0;
            }
            // Q1
            for (int i = doccpi_[h]+soccpi_[h]; i < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; i++) {
                int ii = i - rstcpi_[h] - frzcpi_[h];
                x_p[q1aoff[h]+ii*amopi_[h]+ii] = 1.0;
            }
            for (int i = doccpi_[h]; i < nmopi_[h]-rstvpi_[h] - frzvpi_[h]; i++) {
                int ii = i - rstcpi_[h] - frzcpi_[h];
                x_p[q1boff[h]+ii*amopi_[h]+ii] = 1.0;
            }
        }
    }else { // random guess

        srand(0);
        for (int i = 0; i < dimx_; i++) {
            x_p[i] = ( (double)rand()/RAND_MAX - 1.0 ) * 2.0;
        }
        return;

    }

    if ( constrain_q2_ ) {
        if ( !spin_adapt_q2_) {
            Q2_constraints_guess(x);
        }else {
            Q2_constraints_guess_spin_adapted(x);
        }
    }

    if ( constrain_g2_ ) {
        if ( ! spin_adapt_g2_ ) {
            G2_constraints_guess(x);
        }else {
            G2_constraints_guess_spin_adapted(x);
        }
    }

    if ( constrain_t1_ ) {
        T1_constraints_guess(x);
    }
    if ( constrain_t2_ ) {
        T2_constraints_guess(x);
    }

}

void v2RDMSolver::BuildConstraints(){

    //constraint on the Trace of D2(s=0,ms=0)

    double na = nalpha_ - nfrzc_ - nrstc_;
    double nb = nbeta_ - nfrzc_ - nrstc_;
    double trdab = na * nb;

    //constraint on the Trace of D2(s=1,ms=0)
    double trdaa  = na*(na-1.0);
    double trdbb  = nb*(nb-1.0);
    double n = na + nb;
    double trd    = n*(n-1.0);

    b->zero();
    double* b_p = b->pointer();

    offset = 0;

    ///Trace of D2(s=0,ms=0) and D2(s=1,ms=0)
    if ( constrain_sz_ ) {
        b_p[offset++] = trdab;
        b_p[offset++] = trdaa;
        b_p[offset++] = trdbb;
    }else{
        b_p[offset++] = trd;
    }

    // d1 / q1 a
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                b_p[offset + i*amopi_[h]+j] = (double)(i==j);
            }
        }
        offset += amopi_[h]*amopi_[h];
    }

    // d1 / q1 b
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                b_p[offset + i*amopi_[h]+j] = (double)(i==j);
            }
        }
        offset += amopi_[h]*amopi_[h];
    }

    if ( constrain_sz_ ) {
        //contract D2ab -> D1a
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < amopi_[h]; i++){
                for(int j = 0; j < amopi_[h]; j++){
                    b_p[offset + i*amopi_[h]+j] = 0.0;
                }
            }
            offset += amopi_[h]*amopi_[h];
        }

        //contract D2ab -> D1b
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < amopi_[h]; i++){
                for(int j = 0; j < amopi_[h]; j++){
                    b_p[offset + i*amopi_[h]+j] = 0.0;
                }
            }
            offset += amopi_[h]*amopi_[h];
        }
        //contract D2aa -> D1a
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < amopi_[h]; i++){
                for(int j = 0; j < amopi_[h]; j++){
                    b_p[offset + i*amopi_[h]+j] = 0.0;
                }
            }
            offset += amopi_[h]*amopi_[h];
        }
        //contract D2bb -> D1b
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < amopi_[h]; i++){
                for(int j = 0; j < amopi_[h]; j++){
                    b_p[offset + i*amopi_[h]+j] = 0.0;
                }
            }
            offset += amopi_[h]*amopi_[h];
        }
    }else {
        //contract D2aa + D2ab -> D1a
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < amopi_[h]; i++){
                for(int j = 0; j < amopi_[h]; j++){
                    b_p[offset + i*amopi_[h]+j] = 0.0;
                }
            }
            offset += amopi_[h]*amopi_[h];
        }

        //contract D2bb + D2ab -> D1b
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < amopi_[h]; i++){
                for(int j = 0; j < amopi_[h]; j++){
                    b_p[offset + i*amopi_[h]+j] = 0.0;
                }
            }
            offset += amopi_[h]*amopi_[h];
        }
    }

    if ( constrain_spin_ ) {

        // funny ab trace with spin: N/2 + Ms^2 - S(S+1)
        double ms = (multiplicity_-1.0)/2.0;
        b_p[offset++] = (0.5 * (na + nb) + ms*ms - ms*(ms+1.0));

        // additional spin constraints for singlets:
        if ( nalpha_ == nbeta_ ) {
            for ( int h = 0; h < nirrep_; h++) {
                offset += amopi_[h]*amopi_[h]; // D1a = D1b
            }
            for ( int h = 0; h < nirrep_; h++) {
                offset += gems_aa[h]*gems_aa[h]; // D2aa = D2bb
            }
            for ( int h = 0; h < nirrep_; h++) {
                offset += gems_aa[h]*gems_aa[h]; // D2aa[pq][rs] = 1/2(D2ab[pq][rs] - D2ab[pq][sr] - D2ab[qp][rs] + D2ab[qp][sr])
            }
            for ( int h = 0; h < nirrep_; h++) {
                offset += gems_aa[h]*gems_aa[h]; // D2bb[pq][rs] = 1/2(D2ab[pq][rs] - D2ab[pq][sr] - D2ab[qp][rs] + D2ab[qp][sr]))
            }
            for ( int h = 0; h < nirrep_; h++) {
                offset += gems_ab[h]*gems_ab[h]; // D200[pq][rs] = 1/(sqrt(1+dpq)sqrt(1+drs))(D2ab[pq][rs] + D2ab[pq][sr] + D2ab[qp][rs] + D2ab[qp][sr])
            }
        }else { // nonsinglets
            for ( int h = 0; h < nirrep_; h++) {
                offset += 4*gems_ab[h]*gems_ab[h]; // D200_0, D210_0, D201_0, D211_0
            }
        }

        // maximal spin constraints:
        if ( constrain_g2_ ) {
            for (int i = 0; i < gems_ab[0]; i++){
                b_p[offset++] = 0.0;
            }
            for (int i = 0; i < gems_ab[0]; i++){
                b_p[offset++] = 0.0;
            }
        }

    }

    if ( constrain_q2_ ) {
        if ( !spin_adapt_q2_ ) {
            // map d2ab to q2ab
            for (int h = 0; h < nirrep_; h++) {
                for(int ij = 0; ij < gems_ab[h]; ij++){
                    int i = bas_ab_sym[h][ij][0];
                    int j = bas_ab_sym[h][ij][1];
                    for(int kl = 0; kl < gems_ab[h]; kl++){
                        int k = bas_ab_sym[h][kl][0];
                        int l = bas_ab_sym[h][kl][1];
                        b_p[offset + ij*gems_ab[h]+kl] = -(i==k)*(j==l);
                    }
                }
                offset += gems_ab[h]*gems_ab[h];
            }

            // map d2aa to q2aa
            for (int h = 0; h < nirrep_; h++) {
                for(int ij = 0; ij < gems_aa[h]; ij++){
                    int i = bas_aa_sym[h][ij][0];
                    int j = bas_aa_sym[h][ij][1];
                    for(int kl = 0; kl < gems_aa[h]; kl++){
                        int k = bas_aa_sym[h][kl][0];
                        int l = bas_aa_sym[h][kl][1];
                        b_p[offset + ij*gems_aa[h]+kl] = -(i==k)*(j==l) + (i==l)*(j==k);
                    }
                }
                offset += gems_aa[h]*gems_aa[h];
            }

            // map d2bb to q2bb
            for (int h = 0; h < nirrep_; h++) {
                for(int ij = 0; ij < gems_aa[h]; ij++){
                    int i = bas_aa_sym[h][ij][0];
                    int j = bas_aa_sym[h][ij][1];
                    for(int kl = 0; kl < gems_aa[h]; kl++){
                        int k = bas_aa_sym[h][kl][0];
                        int l = bas_aa_sym[h][kl][1];
                        b_p[offset + ij*gems_aa[h]+kl] = -(i==k)*(j==l) + (i==l)*(j==k);
                    }
                }
                offset += gems_aa[h]*gems_aa[h];
            }
        }else {
            // map d2ab to q2s
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_00[h]; i++){
                    for(int j = 0; j < gems_00[h]; j++){
                        b_p[offset + INDEX(i,j)] = 0.0;
                    }
                }
                offset += gems_00[h]*gems_00[h];
            }
            // map d2ab to q2t
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_aa[h]; i++){
                    for(int j = 0; j < gems_aa[h]; j++){
                        b_p[offset + INDEX(i,j)] = 0.0;
                    }
                }
                offset += gems_aa[h]*gems_aa[h];
            }
            // map d2aa to q2t_p1
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_aa[h]; i++){
                    for(int j = 0; j < gems_aa[h]; j++){
                        b_p[offset + INDEX(i,j)] = 0.0;
                    }
                }
                offset += gems_aa[h]*gems_aa[h];
            }
            // map d2bb to q2t_m1
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_aa[h]; i++){
                    for(int j = 0; j < gems_aa[h]; j++){
                        b_p[offset + INDEX(i,j)] = 0.0;
                    }
                }
                offset += gems_aa[h]*gems_aa[h];
            }
        }
    }

    if ( constrain_g2_ ) {
        if ( ! spin_adapt_g2_ ) {
            // map d2 and d1 to g2ab
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_ab[h]; i++){
                    for(int j = 0; j < gems_ab[h]; j++){
                        b_p[offset + i*gems_ab[h]+j] = 0.0;
                    }
                }
                offset += gems_ab[h]*gems_ab[h];
            }

            // map d2 and d1 to g2ba
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_ab[h]; i++){
                    for(int j = 0; j < gems_ab[h]; j++){
                        b_p[offset + i*gems_ab[h]+j] = 0.0;
                    }
                }
                offset += gems_ab[h]*gems_ab[h];
            }

            // map d2 and d1 to g2aa
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < 2*gems_ab[h]; i++){
                    for(int j = 0; j < 2*gems_ab[h]; j++){
                        b_p[offset + i*2*gems_ab[h]+j] = 0.0;
                    }
                }
                offset += 2*gems_ab[h]*2*gems_ab[h];
            }
        }else {
            // map d2 and d1 to g2s
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_ab[h]; i++){
                    for(int j = 0; j < gems_ab[h]; j++){
                        b_p[offset + i*gems_ab[h]+j] = 0.0;
                    }
                }
                offset += gems_ab[h]*gems_ab[h];
            }
            // map d2 and d1 to g2t
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_ab[h]; i++){
                    for(int j = 0; j < gems_ab[h]; j++){
                        b_p[offset + i*gems_ab[h]+j] = 0.0;
                    }
                }
                offset += gems_ab[h]*gems_ab[h];
            }
            // map d2 and d1 to g2t_p1
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_ab[h]; i++){
                    for(int j = 0; j < gems_ab[h]; j++){
                        b_p[offset + i*gems_ab[h]+j] = 0.0;
                    }
                }
                offset += gems_ab[h]*gems_ab[h];
            }
            // map d2 and d1 to g2t_m1
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_ab[h]; i++){
                    for(int j = 0; j < gems_ab[h]; j++){
                        b_p[offset + i*gems_ab[h]+j] = 0.0;
                    }
                }
                offset += gems_ab[h]*gems_ab[h];
            }
        }
    }

    if ( constrain_t1_ ) {
        // T1aaa
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < trip_aaa[h]; i++){
                for(int j = 0; j < trip_aaa[h]; j++){
                    b_p[offset + i*trip_aaa[h]+j] = 0.0;
                }
            }
            offset += trip_aaa[h]*trip_aaa[h];
        }
        // T1bbb
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < trip_aaa[h]; i++){
                for(int j = 0; j < trip_aaa[h]; j++){
                    b_p[offset + i*trip_aaa[h]+j] = 0.0;
                }
            }
            offset += trip_aaa[h]*trip_aaa[h];
        }
        // T1aab
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < trip_aab[h]; i++){
                for(int j = 0; j < trip_aab[h]; j++){
                    b_p[offset + i*trip_aab[h]+j] = 0.0;
                }
            }
            offset += trip_aab[h]*trip_aab[h];
        }
        // T1bba
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < trip_aab[h]; i++){
                for(int j = 0; j < trip_aab[h]; j++){
                    b_p[offset + i*trip_aab[h]+j] = 0.0;
                }
            }
            offset += trip_aab[h]*trip_aab[h];
        }
    }

    if ( constrain_t2_ ) {
        // T2aaa
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < trip_aba[h]+trip_aab[h]; i++){
                for(int j = 0; j < trip_aba[h] + trip_aab[h]; j++){
                    b_p[offset + i*(trip_aab[h]+trip_aba[h])+j] = 0.0;
                }
            }
            offset += (trip_aba[h]+trip_aab[h])*(trip_aba[h]+trip_aab[h]);
        }
        // T2bbb
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < trip_aba[h]+trip_aab[h]; i++){
                for(int j = 0; j < trip_aba[h] + trip_aab[h]; j++){
                    b_p[offset + i*(trip_aab[h]+trip_aba[h])+j] = 0.0;
                }
            }
            offset += (trip_aba[h]+trip_aab[h])*(trip_aba[h]+trip_aab[h]);
        }
        // T2aab
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < trip_aab[h]; i++){
                for(int j = 0; j < trip_aab[h]; j++){
                    b_p[offset + i*trip_aab[h]+j] = 0.0;
                }
            }
            offset += trip_aab[h]*trip_aab[h];
        }
        // T2bba
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < trip_aab[h]; i++){
                for(int j = 0; j < trip_aab[h]; j++){
                    b_p[offset + i*trip_aab[h]+j] = 0.0;
                }
            }
            offset += trip_aab[h]*trip_aab[h];
        }
    }
    if ( constrain_d3_ ) {
        if (  nalpha_ - nrstc_ - nfrzc_ > 2 ) {
            // D3aaa -> D2aa
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_aa[h]; i++){
                    for(int j = 0; j < gems_aa[h]; j++){
                        b_p[offset + i*gems_aa[h]+j] = 0.0;
                    }
                }
                offset += gems_aa[h]*gems_aa[h];
            }
        }
        if (  nbeta_ - nrstc_ - nfrzc_ > 2 ) {
            // D3bbb -> D2bb
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_aa[h]; i++){
                    for(int j = 0; j < gems_aa[h]; j++){
                        b_p[offset + i*gems_aa[h]+j] = 0.0;
                    }
                }
                offset += gems_aa[h]*gems_aa[h];
            }
        }
        // D3aab -> D2aa
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < gems_aa[h]; i++){
                for(int j = 0; j < gems_aa[h]; j++){
                    b_p[offset + i*gems_aa[h]+j] = 0.0;
                }
            }
            offset += gems_aa[h]*gems_aa[h];
        }
        // D3bba -> D2bb
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < gems_aa[h]; i++){
                for(int j = 0; j < gems_aa[h]; j++){
                    b_p[offset + i*gems_aa[h]+j] = 0.0;
                }
            }
            offset += gems_aa[h]*gems_aa[h];
        }
        if (  nalpha_ - nrstc_ - nfrzc_ > 1 ) {
            // D3aab -> D2ab
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_ab[h]; i++){
                    for(int j = 0; j < gems_ab[h]; j++){
                        b_p[offset + i*gems_ab[h]+j] = 0.0;
                    }
                }
                offset += gems_ab[h]*gems_ab[h];
            }
        }
        if (  nbeta_ - nrstc_ - nfrzc_ > 1 ) {
            // D3bba -> D2ab
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_ab[h]; i++){
                    for(int j = 0; j < gems_ab[h]; j++){
                        b_p[offset + i*gems_ab[h]+j] = 0.0;
                    }
                }
                offset += gems_ab[h]*gems_ab[h];
            }
        }
        // additional spin constraints for singlets:
        if ( constrain_spin_ && nalpha_ == nbeta_ ) {
            for (int h = 0; h < nirrep_; h++) {
                offset += trip_aab[h]*trip_aab[h]; // D3aab = D3bba
            }
            for (int h = 0; h < nirrep_; h++) {
                offset += trip_aaa[h]*trip_aaa[h]; // D3aab -> D3aaa
                offset += trip_aaa[h]*trip_aaa[h]; // D3bba -> D3bbb
            }
        }
    }
    if ( constrain_d4_ ) {
        if (  nalpha_ - nrstc_ - nfrzc_ > 3 ) {
            // D4aaaa -> D3aaa
            for (int h = 0; h < nirrep_; h++) {
                offset += trip_aaa[h]*trip_aaa[h];
            }
        }
        if (  nalpha_ - nrstc_ - nfrzc_ > 2 && nbeta_ - nrstc_ - nfrzc_ > 0) {
            // D4aaab -> D3aaa
            for (int h = 0; h < nirrep_; h++) {
                offset += trip_aaa[h]*trip_aaa[h];
            }
            // D4aaab -> D3aab
            for (int h = 0; h < nirrep_; h++) {
                offset += trip_aab[h]*trip_aab[h];
            }
        }
        if (  nalpha_ - nrstc_ - nfrzc_ > 1 && nbeta_ - nrstc_ - nfrzc_ > 1) {
            // D4aabb -> D3aab
            for (int h = 0; h < nirrep_; h++) {
                offset += trip_aab[h]*trip_aab[h];
            }
            // D4aabb -> D3abb
            for (int h = 0; h < nirrep_; h++) {
                offset += trip_aab[h]*trip_aab[h];
            }
        }
        if (  nalpha_ - nrstc_ - nfrzc_ > 0 && nbeta_ - nrstc_ - nfrzc_ > 2) {
            // D4abbb -> D3bbb
            for (int h = 0; h < nirrep_; h++) {
                offset += trip_aaa[h]*trip_aaa[h];
            }
            // D4abbb -> D3aab
            for (int h = 0; h < nirrep_; h++) {
                offset += trip_aab[h]*trip_aab[h];
            }
        }
        if (  nbeta_ - nrstc_ - nfrzc_ > 3) {
            // D4bbbb -> D3bbb
            for (int h = 0; h < nirrep_; h++) {
                offset += trip_aaa[h]*trip_aaa[h];
            }
        }
    }

}

///Build A dot u where u =[z,c]
void v2RDMSolver::bpsdp_Au(SharedVector A, SharedVector u){

    //A->zero();
    memset((void*)A->pointer(),'\0',nconstraints_*sizeof(double));

    offset = 0;
    D2_constraints_Au(A,u);

    if ( constrain_spin_ ) {
        Spin_constraints_Au(A,u);
    }

    if ( constrain_q2_ ) {
        if ( !spin_adapt_q2_ ) {
            Q2_constraints_Au(A,u);
        }else {
            Q2_constraints_Au_spin_adapted(A,u);
        }
    }

    if ( constrain_g2_ ) {
        if ( ! spin_adapt_g2_ ) {
            G2_constraints_Au(A,u);
        }else {
            G2_constraints_Au_spin_adapted(A,u);
        }
    }

    if ( constrain_t1_ ) {
        T1_constraints_Au(A,u);
    }

    if ( constrain_t2_ ) {
        //T2_constraints_Au(A,u);
        T2_constraints_Au_slow(A,u);
    }

    if ( constrain_d3_ ) {
        D3_constraints_Au(A,u);
    }

    if ( constrain_d4_ ) {
        D4_constraints_Au(A,u);
    }

} // end Au

///Build AT dot u where u =[z,c]
void v2RDMSolver::bpsdp_ATu(SharedVector A, SharedVector u){

    //A->zero();
    memset((void*)A->pointer(),'\0',dimx_*sizeof(double));

    offset = 0;
    D2_constraints_ATu(A,u);

    if ( constrain_spin_ ) {
        Spin_constraints_ATu(A,u);
    }

    if ( constrain_q2_ ) {
        if ( !spin_adapt_q2_ ) {
            Q2_constraints_ATu(A,u);
        }else {
            Q2_constraints_ATu_spin_adapted(A,u);
        }
    }

    if ( constrain_g2_ ) {
        if ( ! spin_adapt_g2_ ) {
            G2_constraints_ATu(A,u);
        }else {
            G2_constraints_ATu_spin_adapted(A,u);
        }
    }

    if ( constrain_t1_ ) {
        T1_constraints_ATu(A,u);
    }

    if ( constrain_t2_ ) {
        //T2_constraints_ATu(A,u);
        T2_constraints_ATu_slow(A,u);
    }
    if ( constrain_d3_ ) {
        D3_constraints_ATu(A,u);
    }

    if ( constrain_d4_ ) {
        D4_constraints_ATu(A,u);
    }

}//end ATu

// TODO: update remaining functions to use restricted vs frozen orbitals
void v2RDMSolver::UnpackDensityPlusCore() {

    memset((void*)d2_plus_core_sym_,'\0',d2_plus_core_dim_*sizeof(double));

    // D2 first
    double * x_p = x->pointer();
    // active active; active active
    int offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i            = bas_ab_sym[h][ij][0];
            int j            = bas_ab_sym[h][ij][1];
            int hi           = symmetry[i];
            int hj           = symmetry[j];
            int ifull        = full_basis[i];
            int jfull        = full_basis[j];
            int ij_ab        = ibas_ab_sym[h][i][j];
            int ji_ab        = ibas_ab_sym[h][j][i];
            int ij_aa        = ibas_aa_sym[h][i][j];

            for (int kl = 0; kl < gems_ab[h]; kl++) {
                int k          = bas_ab_sym[h][kl][0];
                int l          = bas_ab_sym[h][kl][1];
                int hk         = symmetry[k];
                int hl         = symmetry[l];
                int kfull      = full_basis[k];
                int lfull      = full_basis[l];
                int kl_ab      = ibas_ab_sym[h][k][l];
                int lk_ab      = ibas_ab_sym[h][l][k];
                int kl_aa      = ibas_aa_sym[h][k][l];

                //if ( i > k ) continue;
                //if ( j > l ) continue;

                int hik = SymmetryPair(hi,hk);

                int ik_full      = ibas_full_sym[hik][ifull][kfull];
                int jl_full      = ibas_full_sym[hik][jfull][lfull];

                //if ( ik_plus_core > jl_plus_core ) continue;

                int hkj = SymmetryPair(hk,hj);

                int kj_ab = ibas_ab_sym[hkj][k][j];
                int il_ab = ibas_ab_sym[hkj][i][l];

                int jk_ab = ibas_ab_sym[hkj][j][k];
                int li_ab = ibas_ab_sym[hkj][l][i];

                int kj_aa = ibas_aa_sym[hkj][k][j];
                int il_aa = ibas_aa_sym[hkj][i][l];

                offset = 0;
                for (int myh = 0; myh < hik; myh++) {
                    offset += gems_plus_core[myh] * ( gems_plus_core[myh] + 1 ) / 2;
                }
                int id = offset + INDEX(ik_full,jl_full);

                double val = 0.0;

                val += 0.5 * x_p[d2aboff[h]   + ij_ab*gems_ab[h]  + kl_ab];
                val += 0.5 * x_p[d2aboff[hkj] + kj_ab*gems_ab[hkj]+ il_ab] * (1.0 - (double)(l==j));
                val += 0.5 * x_p[d2aboff[hkj] + il_ab*gems_ab[hkj]+ kj_ab] * (1.0 - (double)(i==k));
                val += 0.5 * x_p[d2aboff[h]   + kl_ab*gems_ab[h]  + ij_ab] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));

                val += 0.5 * x_p[d2aboff[h]   + ji_ab*gems_ab[h]  + lk_ab];
                val += 0.5 * x_p[d2aboff[hkj] + jk_ab*gems_ab[hkj]+ li_ab] * (1.0 - (double)(l==j));
                val += 0.5 * x_p[d2aboff[hkj] + li_ab*gems_ab[hkj]+ jk_ab] * (1.0 - (double)(i==k));
                val += 0.5 * x_p[d2aboff[h]   + lk_ab*gems_ab[h]  + ji_ab] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));

                // aa / bb
                if ( i != j && k != l ) {
                    int sij = ( i < j ? 1 : -1 );
                    int skl = ( k < l ? 1 : -1 );
                    val += 0.5 * sij * skl * x_p[d2aaoff[h]   + ij_aa*gems_aa[h]  + kl_aa];
                    val += 0.5 * sij * skl * x_p[d2aaoff[h]   + kl_aa*gems_aa[h]  + ij_aa] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));
                    val += 0.5 * sij * skl * x_p[d2bboff[h]   + ij_aa*gems_aa[h]  + kl_aa];
                    val += 0.5 * sij * skl * x_p[d2bboff[h]   + kl_aa*gems_aa[h]  + ij_aa] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));
                }
                if ( k != j && i != l ) {
                    int skj = ( k < j ? 1 : -1 );
                    int sil = ( i < l ? 1 : -1 );
                    val += 0.5 * skj * sil * x_p[d2aaoff[hkj] + kj_aa*gems_aa[hkj]+ il_aa] * (1.0 - (double)(l==j));
                    val += 0.5 * skj * sil * x_p[d2aaoff[hkj] + il_aa*gems_aa[hkj]+ kj_aa] * (1.0 - (double)(i==k));
                    val += 0.5 * skj * sil * x_p[d2bboff[hkj] + kj_aa*gems_aa[hkj]+ il_aa] * (1.0 - (double)(l==j));
                    val += 0.5 * skj * sil * x_p[d2bboff[hkj] + il_aa*gems_aa[hkj]+ kj_aa] * (1.0 - (double)(i==k));
                }

                // scale the off-diagonal elements
                if ( ik_full != jl_full ) {
                    val *= 2.0;
                }
                d2_plus_core_sym_[id] = val;
            }
        }
    }

    // core core; core core
    double en = 0.0;
    for (int hi = 0; hi < nirrep_; hi++) {
        for (int i = 0; i < rstcpi_[hi] + frzcpi_[hi]; i++) {

            int ifull      = i + pitzer_offset_full[hi];

            for (int hj = 0; hj < nirrep_; hj++) {
                for (int j = 0; j < rstcpi_[hj] + frzcpi_[hj]; j++) {

                    int jfull      = j + pitzer_offset_full[hj];
                    int hij = SymmetryPair(hi,hj);

                    if ( ifull == jfull ) {

                        int iifull = ibas_full_sym[0][ifull][ifull];
                        int jjfull = ibas_full_sym[0][jfull][jfull];

                        d2_plus_core_sym_[INDEX(iifull,jjfull)] =  1.0;

                    }else {

                        int iifull = ibas_full_sym[0][ifull][ifull];
                        int jjfull = ibas_full_sym[0][jfull][jfull];

                        d2_plus_core_sym_[INDEX(iifull,jjfull)] =  4.0;

                        offset = 0;
                        for (int myh = 0; myh < hij; myh++) {
                            offset += gems_plus_core[myh] * ( gems_plus_core[myh] + 1 ) / 2;
                        }

                        int ijfull = ibas_full_sym[hij][ifull][jfull];

                        d2_plus_core_sym_[offset + INDEX(ijfull,ijfull)] = -2.0;
                    }
                }
            }
        }
    }

    // core active; core active
    for (int hi = 0; hi < nirrep_; hi++) {
        for (int i = 0; i < rstcpi_[hi] + frzcpi_[hi]; i++) {
            int ifull      = i + pitzer_offset_full[hi];
            int iifull     = ibas_full_sym[0][ifull][ifull];

            // D2(il; ij) ab, ba, aa, bb
            for (int hj = 0; hj < nirrep_; hj++) {
                for (int j = 0; j < amopi_[hj]; j++) {

                    //int jfull      = full_basis[j];
                    int jfull      = full_basis[j+pitzer_offset[hj]];

                    for (int l = j; l < amopi_[hj]; l++) {

                        //int lfull      = full_basis[l];
                        int lfull      = full_basis[l+pitzer_offset[hj]];

                        int jlfull = ibas_full_sym[0][jfull][lfull];

                        int id = INDEX(iifull,jlfull);

                        double val = 0.0;

                        // aa and bb pieces
                        if ( j == l ) {
                            val += 1.0 * x_p[d1aoff[hj]+j*amopi_[hj]+l];
                            val += 1.0 * x_p[d1boff[hj]+j*amopi_[hj]+l];
                        }else {
                            val += 2.0 * x_p[d1aoff[hj]+j*amopi_[hj]+l];
                            val += 2.0 * x_p[d1boff[hj]+j*amopi_[hj]+l];
                        }

                        // ab and ba pieces
                        if ( j == l ) {
                            val += 1.0 * x_p[d1aoff[hj]+j*amopi_[hj]+l];
                            val += 1.0 * x_p[d1boff[hj]+j*amopi_[hj]+l];
                        }else {
                            val += 2.0 * x_p[d1aoff[hj]+j*amopi_[hj]+l];
                            val += 2.0 * x_p[d1boff[hj]+j*amopi_[hj]+l];
                        }

                        d2_plus_core_sym_[id] = val;

                        // also (il|ji) with a minus sign
                        int hil = SymmetryPair(symmetry_full[ifull],symmetry_full[lfull]);
                        int hji = SymmetryPair(symmetry_full[ifull],symmetry_full[jfull]);

                        int ilfull = ibas_full_sym[hil][ifull][lfull];
                        int jifull = ibas_full_sym[hji][jfull][ifull];

                        //if ( il > ji ) continue;

                        val = 0.0;

                        // aa and bb pieces
                        if ( j == l ) {
                            val -= 1.0 * x_p[d1aoff[hj]+j*amopi_[hj]+l];
                            val -= 1.0 * x_p[d1boff[hj]+j*amopi_[hj]+l];
                        }else {
                            val -= 2.0 * x_p[d1aoff[hj]+j*amopi_[hj]+l];
                            val -= 2.0 * x_p[d1boff[hj]+j*amopi_[hj]+l];
                        }

                        offset = 0;
                        for (int myh = 0; myh < hil; myh++) {
                            offset += gems_plus_core[myh] * ( gems_plus_core[myh] + 1 ) / 2;
                        }

                        d2_plus_core_sym_[offset + INDEX(ilfull,jifull)] = val;

                    }
                }
            }
        }
    }

//gg -- changes relative to original code (see below)
//gg    1) only active orbitals referenced 
//gg    2) off-diagonal elements not scaled by 2

    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {

            for (int j = i; j < amopi_[h]; j++) {

                int id = offset + INDEX(i,j);

                d1_act_spatial_sym_[id]  = x_p[d1aoff[h] + i * amopi_[h] + j];
                d1_act_spatial_sym_[id] += x_p[d1boff[h] + i * amopi_[h] + j];

            }
        }
        offset += amopi_[h] * ( amopi_[h] + 1 ) / 2;
    }

//gg -- remove core-core d1 elements (only active referenced)
/*
    // now D1
    // active; active
    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {

            int iplus_core = i + rstcpi_[h] + frzcpi_[h];

            for (int j = i; j < amopi_[h]; j++) {

                int jplus_core = j + rstcpi_[h] + frzcpi_[h];

                int id = offset + INDEX(iplus_core,jplus_core);

                d1_act_spatial_sym_[id]  = x_p[d1aoff[h] + i * amopi_[h] + j];
                d1_act_spatial_sym_[id] += x_p[d1boff[h] + i * amopi_[h] + j];

                // scale off-diagonal elements
                if ( i != j ) {
                    d1_act_spatial_sym_[id] *= 2.0;
                }

            }
        }
        offset += (rstcpi_[h] + frzcpi_[h] + amopi_[h]) * ( rstcpi_[h] + frzcpi_[h] + amopi_[h] + 1 ) / 2;
    }

    // core; core;
    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < rstcpi_[h] + frzcpi_[h]; i++) {
            d1_act_spatial_sym_[offset + INDEX(i,i)] = 2.0;
        }
        offset += (rstcpi_[h] + frzcpi_[h] + amopi_[h]) * ( rstcpi_[h] + frzcpi_[h] + amopi_[h] + 1 ) / 2;
    }
*/

}

void v2RDMSolver::PackSpatialDensity() {

    memset((void*)d2_act_spatial_sym_,'\0',d2_act_spatial_dim_*sizeof(double));

    // D2 first
    double * x_p = x->pointer();
    // active active; active active
    int offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i            = bas_ab_sym[h][ij][0];
            int j            = bas_ab_sym[h][ij][1];
            int hi           = symmetry[i];
            int hj           = symmetry[j];
            int ij_ab        = ibas_ab_sym[h][i][j];
            int ji_ab        = ibas_ab_sym[h][j][i];
            int ij_aa        = ibas_aa_sym[h][i][j];

            for (int kl = 0; kl < gems_ab[h]; kl++) {
                int k          = bas_ab_sym[h][kl][0];
                int l          = bas_ab_sym[h][kl][1];
                int hk         = symmetry[k];
                int hl         = symmetry[l];
                int kl_ab      = ibas_ab_sym[h][k][l];
                int lk_ab      = ibas_ab_sym[h][l][k];
                int kl_aa      = ibas_aa_sym[h][k][l];

                int hik = SymmetryPair(hi,hk);

                int ik = ibas_00_sym[hik][i][k];
                int jl = ibas_00_sym[hik][j][l];

                int hkj = SymmetryPair(hk,hj);

                int kj_ab = ibas_ab_sym[hkj][k][j];
                int il_ab = ibas_ab_sym[hkj][i][l];

                int jk_ab = ibas_ab_sym[hkj][j][k];
                int li_ab = ibas_ab_sym[hkj][l][i];

                int kj_aa = ibas_aa_sym[hkj][k][j];
                int il_aa = ibas_aa_sym[hkj][i][l];

                offset = 0;
                for (int myh = 0; myh < hik; myh++) {
                    offset += gems_00[myh] * ( gems_00[myh] + 1 ) / 2;
//                    offset += gems_plus_core[myh] * ( gems_plus_core[myh] + 1 ) / 2;
                }
                int id = offset + INDEX(ik,jl);

                double val = 0.0;

                val += 0.5 * x_p[d2aboff[h]   + ij_ab*gems_ab[h]  + kl_ab];
                val += 0.5 * x_p[d2aboff[hkj] + kj_ab*gems_ab[hkj]+ il_ab] * (1.0 - (double)(l==j));
                val += 0.5 * x_p[d2aboff[hkj] + il_ab*gems_ab[hkj]+ kj_ab] * (1.0 - (double)(i==k));
                val += 0.5 * x_p[d2aboff[h]   + kl_ab*gems_ab[h]  + ij_ab] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));

                val += 0.5 * x_p[d2aboff[h]   + ji_ab*gems_ab[h]  + lk_ab];
                val += 0.5 * x_p[d2aboff[hkj] + jk_ab*gems_ab[hkj]+ li_ab] * (1.0 - (double)(l==j));
                val += 0.5 * x_p[d2aboff[hkj] + li_ab*gems_ab[hkj]+ jk_ab] * (1.0 - (double)(i==k));
                val += 0.5 * x_p[d2aboff[h]   + lk_ab*gems_ab[h]  + ji_ab] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));

                // aa / bb
                if ( i != j && k != l ) {
                    int sij = ( i < j ? 1 : -1 );
                    int skl = ( k < l ? 1 : -1 );
                    val += 0.5 * sij * skl * x_p[d2aaoff[h]   + ij_aa*gems_aa[h]  + kl_aa];
                    val += 0.5 * sij * skl * x_p[d2aaoff[h]   + kl_aa*gems_aa[h]  + ij_aa] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));
                    val += 0.5 * sij * skl * x_p[d2bboff[h]   + ij_aa*gems_aa[h]  + kl_aa];
                    val += 0.5 * sij * skl * x_p[d2bboff[h]   + kl_aa*gems_aa[h]  + ij_aa] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));
                }
                if ( k != j && i != l ) {
                    int skj = ( k < j ? 1 : -1 );
                    int sil = ( i < l ? 1 : -1 );
                    val += 0.5 * skj * sil * x_p[d2aaoff[hkj] + kj_aa*gems_aa[hkj]+ il_aa] * (1.0 - (double)(l==j));
                    val += 0.5 * skj * sil * x_p[d2aaoff[hkj] + il_aa*gems_aa[hkj]+ kj_aa] * (1.0 - (double)(i==k));
                    val += 0.5 * skj * sil * x_p[d2bboff[hkj] + kj_aa*gems_aa[hkj]+ il_aa] * (1.0 - (double)(l==j));
                    val += 0.5 * skj * sil * x_p[d2bboff[hkj] + il_aa*gems_aa[hkj]+ kj_aa] * (1.0 - (double)(i==k));
                }

                // scale the off-diagonal elements
                if ( ik != jl ) {
                    val *= 2.0;
                }
                d2_act_spatial_sym_[id] = val;
            }
        }
    }

//gg -- changes relative to original code (see below)
//gg    1) only active orbitals referenced 
//gg    2) off-diagonal elements not scaled by 2

    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {

            for (int j = i; j < amopi_[h]; j++) {

                int id = offset + INDEX(i,j);

                d1_act_spatial_sym_[id]  = x_p[d1aoff[h] + i * amopi_[h] + j];
                d1_act_spatial_sym_[id] += x_p[d1boff[h] + i * amopi_[h] + j];

            }
        }
        offset += amopi_[h] * ( amopi_[h] + 1 ) / 2;
    }

}

void v2RDMSolver::RotateOrbitals(){

    //UnpackDensityPlusCore();
    PackSpatialDensity();

    if ( orbopt_data_[8] > 0 ) {
        outfile->Printf("\n");
        outfile->Printf("        ==> Orbital Optimization <==\n");
        outfile->Printf("\n");
    }

    //int frzc = nfrzc_ + nrstc_;

    // notes for truly frozen core:
    //
    // 1.  symmetry_energy_order should start with first restricted orbital
    // 2.  orbopt_transformation_matrix_ should exclude frozen core

    // notes for truly frozen virtuals:
    // 1.  orbopt_transformation_matrix_ should exclude frozen virtuals
    // 2.  oei_full_dim_, tei_full_dim_ should exclude frozen virtuals
    // 3.  does symmetry_energy_order need to be the right length?

//gg -- added frzcpi_ to argument list
    //OrbOpt(orbopt_transformation_matrix_,
    //      oei_full_sym_,oei_full_dim_,tei_full_sym_,tei_full_dim_,
    //      d1_act_spatial_sym_,d1_act_spatial_dim_,d2_plus_core_sym_,d2_plus_core_dim_,
    //      symmetry_energy_order,frzcpi_,nrstc_,amo_,nrstv_,nirrep_,
    //      orbopt_data_,orbopt_outfile_);

    OrbOpt(orbopt_transformation_matrix_,
          oei_full_sym_,oei_full_dim_,tei_full_sym_,tei_full_dim_,
          d1_act_spatial_sym_,d1_act_spatial_dim_,d2_act_spatial_sym_,d2_act_spatial_dim_,
          symmetry_energy_order,nrstc_,amo_,nrstv_,nirrep_,
          orbopt_data_,orbopt_outfile_,X_);

    if ( orbopt_data_[8] > 0 ) {
        outfile->Printf("            Orbital Optimization %s in %3i iterations \n",(int)orbopt_data_[13] ? "converged" : "did not converge",(int)orbopt_data_[10]);
        outfile->Printf("            Total energy change: %11.6le\n",orbopt_data_[12]);
        outfile->Printf("            Final gradient norm: %11.6le\n",orbopt_data_[11]);
        outfile->Printf("\n");

        if ( fabs(orbopt_data_[12]) < orbopt_data_[4] && fabs(orbopt_data_[11]) < orbopt_data_[3] ) {
            orbopt_converged_ = true;
        }
    }

    RepackIntegrals();
}

} //end namespaces
