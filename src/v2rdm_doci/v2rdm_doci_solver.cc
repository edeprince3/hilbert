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
#include <psi4/libmints/molecule.h>
#include <psi4/libmints/factory.h>
#include <psi4/libmints/basisset.h>
#include <psi4/libmints/local.h>

#include "v2rdm_doci_solver.h"

#include <focas/focas_c_interface.h>
#include <misc/blas.h>
#include <misc/threeindexintegrals.h>
#include <misc/omp.h>
#include <misc/cg_solver.h>
#include <misc/bpsdp_solver.h>

using namespace psi;
using namespace fnocc;

namespace hilbert{

// AATy
static void evaluate_CG_LHS(SharedVector Ax, SharedVector x, void * data) {

    // reinterpret void * as an instance of v2RDM_DOCISolver
    v2RDM_DOCISolver* v2rdm_doci = reinterpret_cast<v2RDM_DOCISolver*>(data);

    Ax->zero();
    std::shared_ptr<Vector> ATy (new Vector(v2rdm_doci->n_primal()));
    v2rdm_doci->bpsdp_ATu(ATy,x);
    v2rdm_doci->bpsdp_Au(Ax,ATy);

}

static void evaluate_Au(std::shared_ptr<Vector> Au, std::shared_ptr<Vector> u, void * data) {

    // reinterpret void * as an instance of v2RDM_DOCISolver
    v2RDM_DOCISolver* v2rdm_doci = reinterpret_cast<v2RDM_DOCISolver*>(data);
    v2rdm_doci->bpsdp_Au(Au,u);

}

static void evaluate_ATu(std::shared_ptr<Vector> ATu, std::shared_ptr<Vector> u, void * data) {

    // reinterpret void * as an instance of v2RDM_DOCISolver
    v2RDM_DOCISolver* v2rdm_doci = reinterpret_cast<v2RDM_DOCISolver*>(data);
    v2rdm_doci->bpsdp_ATu(ATu,u);

}

v2RDM_DOCISolver::v2RDM_DOCISolver(SharedWavefunction reference_wavefunction,Options & options):
    Wavefunction(options){
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

v2RDM_DOCISolver::~v2RDM_DOCISolver()
{
    free(tei_full_sym_);
    free(oei_full_sym_);
    free(d2_plus_core_sym_);
    free(d1_act_spatial_sym_);

    free(amopi_);
    free(rstcpi_);
    free(rstvpi_);

    free(X_);

}

void  v2RDM_DOCISolver::common_init(){

    is_df_ = false;
    if ( options_.get_str("SCF_TYPE") == "DF" || options_.get_str("SCF_TYPE") == "CD" ) {
        is_df_ = true;
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
        throw PsiException("plugin v2rdm-doci works only with c1 symmetry",__FILE__,__LINE__);
    }
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
            throw PsiException("The FROZEN_DOCC array has the wrong dimensions",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            frzcpi_[h] = options_["FROZEN_DOCC"][h].to_double();
        }
    }
    if (options_["RESTRICTED_DOCC"].has_changed()) {
        if (options_["RESTRICTED_DOCC"].size() != nirrep_) {
            throw PsiException("The RESTRICTED_DOCC array has the wrong dimensions",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            rstcpi_[h] = options_["RESTRICTED_DOCC"][h].to_double();
        }
    }
    if (options_["RESTRICTED_UOCC"].has_changed()) {
        if (options_["RESTRICTED_UOCC"].size() != nirrep_) {
            throw PsiException("The RESTRICTED_UOCC array has the wrong dimensions",__FILE__,__LINE__);
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
            throw PsiException("The FROZEN_UOCC array has the wrong dimensions",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            frzvpi_[h] = options_["FROZEN_UOCC"][h].to_double();
        }
    }

    // user could specify active space with ACTIVE array
    if ( options_["ACTIVE"].has_changed() ) {
        //throw PsiException("The ACTIVE array is not yet enabled.",__FILE__,__LINE__);
        if (options_["ACTIVE"].size() != nirrep_) {
            throw PsiException("The ACTIVE array has the wrong dimensions",__FILE__,__LINE__);
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
    name_ = "V2RDM DOCI";

    // pick conditions.  default is dqg
    constrain_q2_ = true;
    constrain_g2_ = true;
    constrain_t1_ = false;
    constrain_t2_ = false;
    constrain_d3_ = false;

    if (options_.get_str("POSITIVITY")=="D") {

        constrain_q2_ = false;
        constrain_g2_ = false;

    }else if (options_.get_str("POSITIVITY")=="DQ") {

        constrain_q2_ = true;
        constrain_g2_ = false;

    }else if (options_.get_str("POSITIVITY")=="DG") {

        constrain_q2_ = false;
        constrain_g2_ = true;

    }else if (options_.get_str("POSITIVITY")=="DQG") {

        constrain_q2_ = true;
        constrain_g2_ = true;

    }else if (options_.get_str("POSITIVITY")=="DQGT2") {

        //throw PsiException("T2 constraints not implemented.",__FILE__,__LINE__);

        constrain_q2_ = true;
        constrain_g2_ = true;
        constrain_t2_ = true;

    }else if (options_.get_str("POSITIVITY")=="DQGT1") {

        //throw PsiException("T1 constraints not implemented.",__FILE__,__LINE__);

        constrain_q2_ = true;
        constrain_g2_ = true;
        constrain_t1_ = true;

    }else if (options_.get_str("POSITIVITY")=="DQGT1T2") {

        //throw PsiException("T2 constraints not implemented.",__FILE__,__LINE__);

        constrain_q2_ = true;
        constrain_g2_ = true;
        constrain_t1_ = true;
        constrain_t2_ = true;

    }else if (options_.get_str("POSITIVITY")=="DQGT2") {

        //throw PsiException("T2 constraints not implemented.",__FILE__,__LINE__);

        constrain_q2_ = true;
        constrain_g2_ = true;
        constrain_t1_ = true;
        constrain_t2_ = true;

    }else {

        throw PsiException("unknown positivity condition",__FILE__,__LINE__);

    }

    if ( options_.get_bool("CONSTRAIN_D3") ) {

        constrain_d3_ = true;

    }

    // build mapping arrays and determine the number of geminals per block
    BuildBasis();

    // dimension of variable buffer (x)

    dimx_ = 0;
    dimx_ += amo_ * amo_; // D2s2
    dimx_ += amo_ * amo_; // D2s0
    dimx_ += amo_;        // D1

    if ( constrain_q2_ ) {
        dimx_ += amo_ * (amo_-1)/2; // Q2s2
        dimx_ += amo_ * amo_; // Q2s0
    }

    if ( constrain_g2_ ) {
        for ( int i = 0; i < amo_; i++) {
            for ( int j = i+1; j < amo_; j++) {
                //if ( i == j ) continue;
                dimx_ += 2*2; // G2s2
            }
        }
        dimx_ += amo_*amo_; // G2s0
    }
    if ( constrain_d3_ ) {
        dimx_ += amo_ * amo_ * amo_;
        dimx_ += amo_ * (amo_ - 1) * (amo_ - 1);
        //for ( int i = 0; i < amo_; i++) {
        //    dimx_ += (amo_ - 1) * (amo_ - 1);
        //}
    }

    if ( constrain_t1_ ) {
        // poelmans thesis
        dimx_ += amo_ * (amo_ - 1) * (amo_ - 2) / 6; // T1(3)aaa
        dimx_ += amo_ * amo_ * (amo_ - 1) / 2;       // T1(3)aab
        dimx_ += amo_ * (amo_ - 1) * (amo_ - 1);
    }
    if ( constrain_t2_ ) {
        // poelmans thesis eq 2.208
        dimx_ += amo_ * (2 * amo_ - 1) * (2 * amo_ - 1);
        // jcp
        for ( int i = 0; i < amo_; i++) {
            for ( int k = i + 1; k < amo_; k++) {
                for ( int m = k + 1; m < amo_; m++) {
                    dimx_ += 3*3; // T2 #1
                }
            }
        }
    }

    // offsets in x

    offset = 0;

    d2s2off_ = offset; offset += amo_*amo_;
    d2s0off_ = offset; offset += amo_*amo_;

    d1off_   = offset; offset += amo_;

    if ( constrain_q2_ ) {

        q2s2off_ = offset; offset += amo_*(amo_-1)/2;
        q2s0off_ = offset; offset += amo_*amo_;

    }

    if ( constrain_g2_ ) {

            g2s2off_ = offset;
            for ( int i = 0; i < amo_; i++) {
                for ( int j = i+1; j < amo_; j++) {
                //if ( i == j ) continue;
                    offset += 2*2;
                }
            }

            g2s0off_ = offset; offset += amo_*amo_;
    }

    if ( constrain_d3_ ) {
        d3s3off_ = offset; offset += amo_ * amo_ * amo_;
        d3s1off_ = offset; offset += amo_ * (amo_ - 1) * (amo_ - 1);
    }

    if ( constrain_t1_ ) {
        t1s3aaaoff_ = offset; offset += amo_ * (amo_ - 1) * (amo_ - 2) / 6; 
        t1s3aaboff_ = offset; offset += amo_ * amo_ * (amo_ - 1) / 2;
        t1s1off_    = offset; offset += amo_ * (amo_ - 1) * (amo_ - 1);
    }
    if ( constrain_t2_ ) {
        t2s1off_ = offset; offset += amo_ * (2 * amo_ - 1) * (2 * amo_ - 1);
        t2off_ = offset;
        for ( int i = 0; i < amo_; i++) {
            for ( int k = i + 1; k < amo_; k++) {
                for ( int m = k + 1; m < amo_; m++) {
                    offset += 3*3; // T2 #1
                }
            }
        }
        //t2s3off_ = offset;
    }

    // constraints:
    nconstraints_ = 0;

    nconstraints_ += 1;             // Tr(D2s2)
    nconstraints_ += 1;             // Tr(D2s0)
    nconstraints_ += 1;             // Tr(D1a)

    nconstraints_ += amo_;          // D2s2 -> D1 a
    nconstraints_ += amo_;          // D2s2 -> D1 b
    nconstraints_ += amo_;          // D2s0 -> D1

    nconstraints_ += amo_ * amo_;          // D2s2 symmetric (with zero diagonal)

    if ( constrain_q2_ ) {
        nconstraints_ += amo_*(amo_-1)/2; // Q2s2
        nconstraints_ += amo_*amo_;     // Q2s0
    }

    if ( constrain_g2_ ) {
        for ( int i = 0; i < amo_; i++) {
            for ( int j = i+1; j < amo_; j++) {
                //if ( i == j ) continue;
                nconstraints_ += 2*2; // G2s2
            }
        }
        nconstraints_ += amo_*amo_; // G2s0
    }

    if ( constrain_d3_ ) {
        nconstraints_ += amo_ * (amo_ - 1); // D3s3 -> D2s2 contraction
        nconstraints_ += amo_ * (amo_ - 1); // D3s1 -> D2s2 contraction
        nconstraints_ += amo_ * (amo_ - 1); // D3s1 -> D2s2 contraction (again)
        nconstraints_ += amo_ * amo_;       // D3s1 -> D2s0 contraction
    }

    if ( constrain_t1_ ) {
        nconstraints_ += amo_ * (amo_ - 1) * (amo_ - 2) / 6;
        nconstraints_ += amo_ * amo_ * (amo_ - 1) / 2;
        nconstraints_ += amo_ * (amo_ - 1) * (amo_ - 1);
    }
    if ( constrain_t2_ ) {

        nconstraints_ += amo_ * (2 * amo_ - 1) * (2 * amo_ - 1);
        for ( int i = 0; i < amo_; i++) {
            for ( int k = i + 1; k < amo_; k++) {
                for ( int m = k + 1; m < amo_; m++) {
                    nconstraints_ += 3*3; // T2 #1
                }
            }
        }
    }

    // list of dimensions_
    for (int ij = 0; ij < amo_ * amo_; ij++) {
        dimensions_.push_back(1); // D2s2
    }
    dimensions_.push_back(amo_);  // D2s0

    for (int i = 0; i < amo_; i++) {
        dimensions_.push_back(1); // D1a
    }

    if ( constrain_q2_ ) {

        for (int ij = 0; ij < amo_ * (amo_-1)/2; ij++) {
            dimensions_.push_back(1); // Q2s2
        }
        dimensions_.push_back(amo_);  // Q2s0

    }

    if ( constrain_g2_ ) {
        for ( int i = 0; i < amo_; i++) {
            for ( int j = i+1; j < amo_; j++) {
                //if ( i == j ) continue;
                dimensions_.push_back(2); // G2s2
            }
        }
        dimensions_.push_back(amo_); // G2s0
    }

    if ( constrain_d3_ ) {
        for ( int i = 0; i < amo_ * amo_ * amo_; i++) {
            dimensions_.push_back(1); // D3s3
        }
        for ( int i = 0; i < amo_; i++) {
            dimensions_.push_back( amo_ - 1 ); // D3s1
        }
    }
    if ( constrain_t1_ ) {
        for ( int i = 0; i < amo_ * (amo_ - 1) * (amo_ - 2) / 6; i++) {
            dimensions_.push_back(1); // T1s3aaa
        }
        for ( int i = 0; i < amo_ * amo_ * (amo_ - 1) / 2; i++) {
            dimensions_.push_back(1); // T1s3aab
        }
        for ( int i = 0; i < amo_; i++) {
            dimensions_.push_back( amo_ - 1 ); // T1s1
        }
    }
    if ( constrain_t2_ ) {
        for ( int b = 0; b < amo_; b++) {
            dimensions_.push_back(2 * amo_ - 1);
        } 
        for ( int i = 0; i < amo_; i++) {
            for ( int k = i + 1; k < amo_; k++) {
                for ( int m = k + 1; m < amo_; m++) {
                    dimensions_.push_back(3); // T2 #1
                }
            }
        }
    }

    // v2rdm sdp convergence thresholds:
    r_convergence_  = options_.get_double("R_CONVERGENCE");
    e_convergence_  = options_.get_double("E_CONVERGENCE");
    maxiter_        = options_.get_int("MAXITER");

    // conjugate gradient solver thresholds:
    cg_convergence_ = options_.get_double("CG_CONVERGENCE");
    cg_maxiter_     = options_.get_double("CG_MAXITER");

    // memory check happens here

    outfile->Printf("\n\n");
    outfile->Printf( "        ***************************************************\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *    v2RDM-DOCI                                   *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *    A variational 2-RDM-driven approach to       *\n");
    outfile->Printf( "        *    doubly occupied configuration interaction    *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        ***************************************************\n");

    outfile->Printf("\n");
    outfile->Printf("  ==> Convergence parameters <==\n");
    outfile->Printf("\n");
    outfile->Printf("        r_convergence:                      %5.3le\n",r_convergence_);
    outfile->Printf("        e_convergence:                      %5.3le\n",e_convergence_);
    outfile->Printf("        cg_convergence:                     %5.3le\n",cg_convergence_);
    outfile->Printf("        maxiter:                             %8i\n",maxiter_);
    outfile->Printf("        cg_maxiter:                          %8i\n",cg_maxiter_);
    outfile->Printf("\n");

    // print orbitals per irrep in each space
    outfile->Printf("  ==> Active space details <==\n");
    outfile->Printf("\n");
    //outfile->Printf("        Freeze core orbitals?                   %5s\n",nfrzc_ > 0 ? "yes" : "no");
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
// gg
    outfile->Printf("        1-step algorithm:                   %5s\n","true");//options_.get_bool("ORBOPT_ONE_STEP") ? "true" : "false");
    outfile->Printf("        g_convergence:                  %5.3le\n",options_.get_double("ORBOPT_GRADIENT_CONVERGENCE"));
    outfile->Printf("        e_convergence:                  %5.3le\n",options_.get_double("ORBOPT_ENERGY_CONVERGENCE"));
    outfile->Printf("        maximum iterations:                 %5i\n",options_.get_int("ORBOPT_MAXITER"));
    outfile->Printf("        frequency:                          %5i\n",options_.get_int("ORBOPT_FREQUENCY"));
    outfile->Printf("        active-active rotations:            %5s\n","true");//options_.get_bool("ORBOPT_ACTIVE_ACTIVE_ROTATIONS") ? "true" : "false");
    outfile->Printf("        exact diagonal Hessian:             %5s\n",options_.get_bool("ORBOPT_EXACT_DIAGONAL_HESSIAN") ? "true" : "false");
    outfile->Printf("        number of DIIS vectors:             %5i\n",options_.get_int("ORBOPT_NUM_DIIS_VECTORS"));
    outfile->Printf("        print iteration info:               %5s\n",options_.get_bool("ORBOPT_WRITE") ? "true" : "false");
// gg

    outfile->Printf("\n");
    outfile->Printf("  ==> Memory requirements <==\n");
    outfile->Printf("\n");
    int nd2   = 0;
    int ng2    = 0;
    int maxgem = 0;
    // AED: some of this is still wrong (e.g. maxgem)
    for (int h = 0; h < nirrep_; h++) {
        nd2 +=     amopi_[h]*amopi_[h];
        nd2 +=     amopi_[h]*amopi_[h];

        ng2 +=     2*amopi_[h] * 2*amopi_[h]; // G2s2
        ng2 +=     amopi_[h]*amopi_[h]; // G2s0

        if ( gems_ab[h] > maxgem ) {
            maxgem = gems_ab[h];
        }
        if ( constrain_g2_ ) {
            if ( 2*gems_ab[h] > maxgem ) {
                maxgem = 2*gems_ab[h];
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

    // for doci, need d2 and 3- or 4-index integrals

    // storage requirements for full d2
    for (int h = 0; h < nirrep_; h++) {
        //tot += gems_plus_core[h] * ( gems_plus_core[h] + 1 ) / 2;
        tot += amopi_[h]*amopi_[h];
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

    // memory available after allocating all we need for v2RDM-DOCI
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
    //if ( options_.get_str("RESTART_FROM_CHECKPOINT_FILE") != "" ) {
    //    ReadOrbitalsFromCheckpointFile();
    //}

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
    orbopt_data_[14] = 3.0;
    if      ( options_.get_str("ORBOPT_ALGORITHM") == "STEEPEST_DESCENT" ) orbopt_data_[14] = 0.0;
    else if ( options_.get_str("ORBOPT_ALGORITHM") == "HESTENES_STIEFEL" ) orbopt_data_[14] = 1.0;
    else if ( options_.get_str("ORBOPT_ALGORITHM") == "DAI_YUAN" )         orbopt_data_[14] = 2.0;
    else if ( options_.get_str("ORBOPT_ALGORITHM") == "HAGER_ZHANG" )      orbopt_data_[14] = 3.0;
    else if ( options_.get_str("ORBOPT_ALGORITHM") == "KOU_DAI" )          orbopt_data_[14] = 4.0;

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

}

int v2RDM_DOCISolver::SymmetryPair(int i,int j) {
    return table[i*8+j];
}
int v2RDM_DOCISolver::TotalSym(int i,int j,int k, int l) {
    return SymmetryPair(SymmetryPair(symmetry[i],symmetry[j]),SymmetryPair(symmetry[k],symmetry[l]));
}

// compute the energy!
double v2RDM_DOCISolver::compute_energy() {

    double start_total_time = omp_get_wtime();

    // hartree-fock guess
    Guess();

    // checkpoint file
    //if ( options_.get_str("RESTART_FROM_CHECKPOINT_FILE") != "" ) {
    //    ReadFromCheckpointFile();
    //}

    // get integrals
    GetIntegrals();

    // generate constraint vector
    BuildConstraints();

    // sdp solver
    std::shared_ptr<BPSDPSolver> sdp (new BPSDPSolver(dimx_,nconstraints_,options_));

    // iterate
    int orbopt_iter = 0;
    do { 

        sdp->solve(x, b, c, dimensions_, options_.get_int("ORBOPT_FREQUENCY"), evaluate_Au, evaluate_ATu, evaluate_CG_LHS, (void*)this);

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

    }while( !orbopt_converged_);

    outfile->Printf("\n");
    outfile->Printf("      v2RDM-DOCI iterations converged!\n");
    outfile->Printf("\n");

    double energy_primal = C_DDOT(dimx_,c->pointer(),1,x->pointer(),1);
    outfile->Printf("    * v2RDM-DOCI total energy:        %20.12lf\n",energy_primal+enuc_+efzc_);
    outfile->Printf("\n");

    std::shared_ptr<Vector> NOs (new Vector("Natural Orbital Occupation Numbers (spin free)",nmo_));
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            int ii = i + pitzer_offset[h];
            NOs->pointer()[ii] = 2.0 * x->pointer()[d1off_ + ii];
        }
    }
    NOs->print();

    Process::environment.globals["CURRENT ENERGY"]     = energy_primal+enuc_+efzc_;
    Process::environment.globals["v2RDM-DOCI TOTAL ENERGY"] = energy_primal+enuc_+efzc_;

    // push final transformation matrix onto Ca_ and Cb_
    if ( options_.get_bool("OPTIMIZE_ORBITALS") ) {
        UpdateTransformationMatrix(reference_wavefunction_,newMO_,Ca_,Cb_,orbopt_transformation_matrix_);
    }

    if ( options_.get_bool("MOLDEN_WRITE") ) {
        WriteMoldenFile();
    }

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
            UpdateTransformationMatrix(reference_wavefunction_,newMO_,Ca_,Cb_,orbopt_transformation_matrix_);

            // transform D1, D2, D3 to semicanonical basis
            UpdatePrimal();
            //printf("primal energy after transformation:        %20.12lf\n",C_DDOT(dimx_,c->pointer(),1,x->pointer(),1)+efzc_);
    
        }
    } 

    if ( options_.get_bool("TPDM_WRITE_FULL") ) {
        WriteTPDM();
        //ReadTPDM();
    }

    // for derivatives:
    if ( options_.get_str("DERTYPE") == "FIRST" ) {

        throw PsiException("derivatives not yet enabled for v2rdm-doci",__FILE__,__LINE__);

        //// write checkpoint file for next step in optimization
        //WriteCheckpointFile();

        //orbopt_data_[8] = -1.0;
        //RotateOrbitals();

        //// write 2-RDM in IWL format
        //WriteTPDM_IWL();

        //// push orbital lagrangian onto wave function
        //OrbitalLagrangian();

    }

    double end_total_time = omp_get_wtime();

    outfile->Printf("\n");
    outfile->Printf("  ==> Iteration count <==\n");
    outfile->Printf("\n");
    outfile->Printf("      Microiterations:            %12li\n",sdp->iiter_total());
    outfile->Printf("      Macroiterations:            %12li\n",sdp->oiter_total());
    outfile->Printf("      Orbital optimization steps: %12li\n",orbopt_iter_total_);
    outfile->Printf("\n");
    outfile->Printf("  ==> Wall time <==\n");
    outfile->Printf("\n");
    outfile->Printf("      Microiterations:            %12.2lf s\n",sdp->iiter_time());
    outfile->Printf("      Macroiterations:            %12.2lf s\n",sdp->oiter_time());
    outfile->Printf("      Orbital optimization:       %12.2lf s\n",orbopt_time_);
    outfile->Printf("      Total:                      %12.2lf s\n",end_total_time - start_total_time);
    outfile->Printf("\n");

    return energy_primal + enuc_ + efzc_;
}

void v2RDM_DOCISolver::WriteMoldenFile() {

    throw PsiException("function WriteMoldenFile() is currently broken",__FILE__,__LINE__);

/*
    std::shared_ptr<Matrix> D (new Matrix(nirrep_,nmopi_,nmopi_));
    std::shared_ptr<Matrix> eigvec (new Matrix(nirrep_,nmopi_,nmopi_));
    std::shared_ptr<Vector> eigval (new Vector("Natural Orbital Occupation Numbers (spin free)",nirrep_,nmopi_));
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < frzcpi_[h] + rstcpi_[h]; i++) {
            D->pointer(h)[i][i] = 2.0;
        }
        for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h] - rstvpi_[h] - frzvpi_[h]; i++) {
            for (int j = rstcpi_[h] + frzcpi_[h]; j < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; j++) {
                D->pointer(h)[i][j]  = 2.0 * x->pointer()[d1off[h]+(i-rstcpi_[h]-frzcpi_[h])*amopi_[h]+(j-rstcpi_[h]-frzcpi_[h])];
            }
        }
    }
    std::shared_ptr<Matrix> saved ( new Matrix(D) );
    D->diagonalize(eigvec,eigval,descending);
    eigval->print();

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
        throw PsiException("printing orbitals is currently disabled when restarting v2rdm jobs.  sorry!",__FILE__,__LINE__);
    }
    //std::shared_ptr<MoldenWriter> molden(new MoldenWriter((std::shared_ptr<Wavefunction>)this));
    std::shared_ptr<MoldenWriter> molden(new MoldenWriter(reference_wavefunction_));
    std::shared_ptr<Vector> zero (new Vector("",nirrep_,nmopi_));
    zero->zero();
    std::string filename = get_writer_file_prefix(reference_wavefunction_->molecule()->name()) + ".molden";
    molden->write(filename,Cno,Cno,zero, zero,eigval,eigval,true);
*/

}

void v2RDM_DOCISolver::NaturalOrbitals() {

    throw PsiException("function NaturalOrbitals() is currently broken",__FILE__,__LINE__);

/*
    SharedVector eigval (new Vector("Natural Orbital Occupation Numbers",nirrep_,nmopi_));
    for (int h = 0; h < nirrep_; h++) {
        double * ep = eigval->pointer(h);
        for (int i = 0; i < frzcpi_[h] + rstcpi_[h]; i++) {
            ep[i] = 1.0;
        }
        for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h] - rstvpi_[h] - frzvpi_[h]; i++) {
            ep[i] = x->pointer()[d1off[h]+(i-rstcpi_[h]-frzcpi_[h])];
        }
    }
    eigval->print();
*/

}

void v2RDM_DOCISolver::MullikenPopulations() {

    throw PsiException("function NaturalOrbitals() is currently broken",__FILE__,__LINE__);

/*
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
                opdm_a->pointer(h)[i][j] = x->pointer()[d1off[h]+(i-rstcpi_[h]-frzcpi_[h])*amopi_[h]+(j-rstcpi_[h]-frzcpi_[h])];
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

    // DOCI: D1a = D1b
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < rstcpi_[h]+frzcpi_[h]; i++) {
            opdm_b->pointer(h)[i][i] = 1.0;
        }
        for (int i = rstcpi_[h]+frzcpi_[h]; i < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; i++) {
            for (int j = rstcpi_[h]+frzcpi_[h]; j < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; j++) {
                opdm_b->pointer(h)[i][j] = x->pointer()[d1off[h]+(i-rstcpi_[h]-frzcpi_[h])*amopi_[h]+(j-rstcpi_[h]-frzcpi_[h])];
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
*/

}

void v2RDM_DOCISolver::Guess(){

    double* x_p = x->pointer();

    memset((void*)x_p,'\0',dimx_*sizeof(double));

    srand(0);
    for (int i = 0; i < dimx_; i++) {
        x_p[i] = ( (double)rand()/RAND_MAX - 1.0 ) * 2.0;
    }

}

void v2RDM_DOCISolver::BuildConstraints(){

    //constraint on the Trace of D2(s=0,ms=0)

    b->zero();
    double* b_p = b->pointer();

    offset = 0;

    ///Trace of D2(seniority-2) and D2(seniority-0)

    b_p[offset++] = nalpha_ * ( nalpha_ - 1.0);
    b_p[offset++] = nalpha_;

    ///Trace of D1
    b_p[offset++] = nalpha_;

    //contract D2s2 -> D1a
    for(int i = 0; i < amo_; i++){
        b_p[offset + i] = 0.0;
    }
    offset += amo_;

    //contract D2s2 -> D1b
    for(int i = 0; i < amo_; i++){
        b_p[offset + i] = 0.0;
    }
    offset += amo_;

    //contract D2s0 -> D1a
    for(int i = 0; i < amo_; i++){
        b_p[offset + i] = 0.0;
    }
    offset += amo_;

    // D2_2 symmetric (with zero diagonal)
    for (int i = 0; i < amo_; i++){
        for (int j = 0; j < amo_; j++){
            b_p[offset++] = 0.0;
        }
    }

    if ( constrain_q2_ ) {

        // map d2s2 to q2s2 seniority-2
        int ij = 0;
        for(int i = 0; i < amo_; i++){
            for(int j = i + 1; j < amo_; j++){
                //if ( i == j ) continue;
                b_p[offset + ij] = -1.0;
                ij++;
            }
        }
        offset += amo_*(amo_-1)/2;

        // map d2s0 to q2s0 seniority-0
        for(int i = 0; i < amo_; i++){
            for(int j = 0; j < amo_; j++){
                b_p[offset + i*amo_+j] = -1.0 * (i==j);
            }
        }
        offset += amo_*amo_;
    }

    if ( constrain_g2_ ) {
 
        // map d2s2, d1 and d2s0 to g2s2
        for(int i = 0; i < amo_; i++){
            for(int j = i+1; j < amo_; j++){
                //if ( i == j ) continue;
                b_p[offset] = 0.0;
                b_p[offset + 1] = 0.0;
                b_p[offset + 2] = 0.0;
                b_p[offset + 3] = 0.0;

                offset += 2*2;
            }
        }

        // map d2s2 and d1 to g2s0
        for(int i = 0; i < amo_; i++){
            for(int j = 0; j < amo_; j++){
                b_p[offset + i*amo_+j] = 0.0;
            }
        }
        offset += amo_*amo_;

/*
        // map d2s2 and d1 to g2s0 (again)
        for(int i = 0; i < amo_; i++){
            for(int j = 0; j < amo_; j++){
                b_p[offset + i*amo_+j] = 0.0;
            }
        }
        offset += amo_*amo_;
*/
    }

    if ( constrain_d3_ ) {
        for ( int i = 0; i < amo_ * (amo_ - 1 ); i++) {
            b_p[offset++] = 0.0; //D3s3 -> D2s2
        }
        for ( int i = 0; i < amo_ * (amo_ - 1); i++) {
            b_p[offset++] = 0.0; //D3s1 -> D2s2
            b_p[offset++] = 0.0; //D3s1 -> D2s2 (again)
        }
        for ( int i = 0; i < amo_ * amo_; i++) {
            b_p[offset++] = 0.0; //D3s1 -> D2s0
        }
    }
    if ( constrain_t1_ ) { 
        for ( int i = 0; i < amo_ * (amo_ - 1) * (amo_ - 2) / 6; i++) {
            b_p[offset++] = 1.0; // T1s3aaa
        }
        for ( int i = 0; i < amo_ * amo_ * (amo_ - 1) / 2; i++) {
            b_p[offset++] = 1.0; // T1s3aab
        }
        for (int b = 0; b < amo_; b++) {
            for (int a = 0; a < amo_; a++) {

                if ( a == b ) continue;

                for (int c = 0; c < amo_; c++) {

                    if ( c == b ) continue;

                    int mya = a;
                    int myc = c;

                    if ( a > b ) mya--;
                    if ( c > b ) myc--;

                    b_p[offset++] = (double)(mya==myc); // T1s1
                }
            }
        }
    }
    if ( constrain_t2_ ) {
        for ( int b = 0; b < amo_ * (2 * amo_ - 1) * (2 * amo_ - 1); b++) {
            b_p[offset++] = 0.0;
        }
        //for ( int b = 0; b < amo_; b++) {
        //    for ( int ac = 0; ac < (2 * amo_ - 1) * (2 * amo_ - 1); ac++) {
        //        //b_p[offset + ac] = 0.0;
        //        b_p[offset++] = 0.0;
        //    }
        //    //b_p[offset + (2 * amo_ - 2) * (2 * amo_ - 1) + 2 * amo_ - 2] = 1.0;
        //    //offset += (2 * amo_ - 1) * (2 * amo_ - 1);
        //}
        for ( int i = 0; i < amo_; i++) {
            for ( int k = i + 1; k < amo_; k++) {
                for ( int m = k + 1; m < amo_; m++) {
                    b_p[offset + 0] = 0.0;
                    b_p[offset + 1] = 0.0;
                    b_p[offset + 2] = 0.0;
                    b_p[offset + 3] = 0.0;
                    b_p[offset + 4] = 0.0;
                    b_p[offset + 5] = 0.0;
                    b_p[offset + 6] = 0.0;
                    b_p[offset + 7] = 0.0;
                    b_p[offset + 8] = 0.0;

                    offset += 9;
                }
            }
        }
    }

}

///Build A dot u where u =[z,c]
void v2RDM_DOCISolver::bpsdp_Au(SharedVector A, SharedVector u){

    memset((void*)A->pointer(),'\0',nconstraints_*sizeof(double));

    offset = 0;
    D2_constraints_Au(A,u);

    if ( constrain_q2_ ) {
        Q2_constraints_Au(A,u);
    }

    if ( constrain_g2_ ) {
        G2_constraints_Au(A,u);
    }

    if ( constrain_d3_ ) {
        D3_constraints_Au(A,u);
    }

    if ( constrain_t1_ ) {
        T1_constraints_Au(A,u);
    }

    if ( constrain_t2_ ) {
        T2_constraints_Au(A,u);
    }

} // end Au

///Build AT dot u where u =[z,c]
void v2RDM_DOCISolver::bpsdp_ATu(SharedVector A, SharedVector u){

    memset((void*)A->pointer(),'\0',dimx_*sizeof(double));

    offset = 0;
    D2_constraints_ATu(A,u);

    if ( constrain_q2_ ) {
        Q2_constraints_ATu(A,u);
    }

    if ( constrain_g2_ ) {
        G2_constraints_ATu(A,u);
    }

    if ( constrain_d3_ ) {
        D3_constraints_ATu(A,u);
    }

    if ( constrain_t1_ ) {
        T1_constraints_ATu(A,u);
    }

    if ( constrain_t2_ ) {
        T2_constraints_ATu(A,u);
    }

}//end ATu

void v2RDM_DOCISolver::cg_Ax(long int N,SharedVector A,SharedVector ux){

//    A->zero();
//    bpsdp_ATu(ATy,ux);
//    bpsdp_Au(A,ATy);

}//end cg_Ax

void v2RDM_DOCISolver::PackSpatialDensity() {

    memset((void*)d2_act_spatial_sym_,'\0',d2_act_spatial_dim_*sizeof(double));

    double * D2aa = (double*)malloc(gems_aa[0]*gems_aa[0]*sizeof(double));
    double * D2ab = (double*)malloc(gems_ab[0]*gems_ab[0]*sizeof(double));
    memset((void*)D2aa,'\0',gems_aa[0]*gems_aa[0]*sizeof(double));
    memset((void*)D2ab,'\0',gems_ab[0]*gems_ab[0]*sizeof(double));

    // d2

    double * x_p = x->pointer();

    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            D2aa[ij*gems_aa[h] + ij] = x_p[d2s2off_ + i*amo_ + j]; //D2_2_[i*nmo_+j];
        }
    }

    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];
            if ( i == j ) continue;
            D2ab[ij*gems_ab[h] + ij]  = x_p[d2s2off_ + i*amo_ + j]; //D2_2_[i*nmo_+j];
        }
    }

    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            int ii = ibas_ab_sym[0][i][i];
            int jj = ibas_ab_sym[0][j][j];
            D2ab[ii*gems_ab[0] + jj] = x_p[d2s0off_ + i*amo_ + j]; //D2_0_[i*nmo_+j];
        }
    }

    // active active; active active
    int offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i            = bas_ab_sym[h][ij][0];
            int j            = bas_ab_sym[h][ij][1];
            int hi           = 0;
            int hj           = 0;
            int ij_ab        = ibas_ab_sym[h][i][j];
            int ji_ab        = ibas_ab_sym[h][j][i];
            int ij_aa        = ibas_aa_sym[h][i][j];

            for (int kl = 0; kl < gems_ab[h]; kl++) {
                int k          = bas_ab_sym[h][kl][0];
                int l          = bas_ab_sym[h][kl][1];
                int hk         = 0;
                int hl         = 0;
                int kl_ab      = ibas_ab_sym[h][k][l];
                int lk_ab      = ibas_ab_sym[h][l][k];
                int kl_aa      = ibas_aa_sym[h][k][l];

                int hik = 0;//SymmetryPair(hi,hk);

                int ik = ibas_00_sym[hik][i][k];
                int jl = ibas_00_sym[hik][j][l];

                int hkj = 0;//SymmetryPair(hk,hj);

                int kj_ab = ibas_ab_sym[hkj][k][j];
                int il_ab = ibas_ab_sym[hkj][i][l];

                int jk_ab = ibas_ab_sym[hkj][j][k];
                int li_ab = ibas_ab_sym[hkj][l][i];

                int kj_aa = ibas_aa_sym[hkj][k][j];
                int il_aa = ibas_aa_sym[hkj][i][l];

                int id = INDEX(ik,jl);

                double val = 0.0;

                val += 0.5 * D2ab[ij_ab*gems_ab[h]  + kl_ab];
                val += 0.5 * D2ab[kj_ab*gems_ab[hkj]+ il_ab] * (1.0 - (double)(l==j));
                val += 0.5 * D2ab[il_ab*gems_ab[hkj]+ kj_ab] * (1.0 - (double)(i==k));
                val += 0.5 * D2ab[kl_ab*gems_ab[h]  + ij_ab] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));

                val += 0.5 * D2ab[ji_ab*gems_ab[h]  + lk_ab];
                val += 0.5 * D2ab[jk_ab*gems_ab[hkj]+ li_ab] * (1.0 - (double)(l==j));
                val += 0.5 * D2ab[li_ab*gems_ab[hkj]+ jk_ab] * (1.0 - (double)(i==k));
                val += 0.5 * D2ab[lk_ab*gems_ab[h]  + ji_ab] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));

                // aa / bb
                if ( i != j && k != l ) {
                    int sij = ( i < j ? 1 : -1 );
                    int skl = ( k < l ? 1 : -1 );
                    val += 0.5 * sij * skl * D2aa[ij_aa*gems_aa[h]  + kl_aa];
                    val += 0.5 * sij * skl * D2aa[kl_aa*gems_aa[h]  + ij_aa] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));
                    val += 0.5 * sij * skl * D2aa[ij_aa*gems_aa[h]  + kl_aa];
                    val += 0.5 * sij * skl * D2aa[kl_aa*gems_aa[h]  + ij_aa] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));
                }
                if ( k != j && i != l ) {
                    int skj = ( k < j ? 1 : -1 );
                    int sil = ( i < l ? 1 : -1 );
                    val += 0.5 * skj * sil * D2aa[kj_aa*gems_aa[hkj]+ il_aa] * (1.0 - (double)(l==j));
                    val += 0.5 * skj * sil * D2aa[il_aa*gems_aa[hkj]+ kj_aa] * (1.0 - (double)(i==k));
                    val += 0.5 * skj * sil * D2aa[kj_aa*gems_aa[hkj]+ il_aa] * (1.0 - (double)(l==j));
                    val += 0.5 * skj * sil * D2aa[il_aa*gems_aa[hkj]+ kj_aa] * (1.0 - (double)(i==k));
                }

                // scale the off-diagonal elements
                if ( ik != jl ) {
                    val *= 2.0;
                }
                d2_act_spatial_sym_[id] = val;
            }
        }
    }

    free(D2aa);
    free(D2ab);

    // d1

    memset((void*)d1_act_spatial_sym_,'\0',d1_act_spatial_dim_*sizeof(double));

    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {

            int ii = i + pitzer_offset[h];
            d1_act_spatial_sym_[offset + INDEX(i,i)] = 2.0 * x_p[d1off_ + ii];

        }
        offset += amopi_[h] * ( amopi_[h] + 1 ) / 2;

    }

}

void v2RDM_DOCISolver::RotateOrbitals(){

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

    double * X_ = (double*)malloc(nmo_*nmo_*sizeof(double));

    int * symmetry_energy_order  = (int*)malloc(nmo_*sizeof(int));
    for (int i = 0; i < nmo_; i++) {
        symmetry_energy_order[i] = 1;
    }

    int nrstc = 0;
    int nrstv = 0;

    OrbOpt(orbopt_transformation_matrix_,
          oei_full_sym_,oei_full_dim_,tei_full_sym_,tei_full_dim_,
          d1_act_spatial_sym_,d1_act_spatial_dim_,d2_act_spatial_sym_,d2_act_spatial_dim_,
          symmetry_energy_order,nrstc,nmo_,nrstv,nirrep_,
          orbopt_data_,orbopt_outfile_,X_);

    free(X_);
    free(symmetry_energy_order);
    
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
