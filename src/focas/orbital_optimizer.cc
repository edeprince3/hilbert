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

#include <string>

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

#include "orbital_optimizer.h"

using namespace psi;
using namespace fnocc;

namespace hilbert{

OrbitalOptimizer::OrbitalOptimizer(std::shared_ptr<Wavefunction> reference_wavefunction, Options & options){

    is_energy_converged_ = false;
    is_gradient_converged_ = false;

    nirrep_ = reference_wavefunction->nirrep();
    nmo_    = reference_wavefunction->nmo();
    doccpi_ = reference_wavefunction->doccpi();
    frzcpi_ = reference_wavefunction->frzcpi();
    frzvpi_ = reference_wavefunction->frzvpi();
    nmopi_  = reference_wavefunction->nmopi();

    if ( reference_wavefunction->options().get_str("SCF_TYPE") == "DF" ) {
        std::shared_ptr<BasisSet> primary = reference_wavefunction->basisset();
        std::shared_ptr<BasisSet> auxiliary = reference_wavefunction->get_basisset("DF_BASIS_SCF");
        nQ_ = auxiliary->nbf();
    }else if ( reference_wavefunction->options().get_str("SCF_TYPE") == "CD" ) {
        std::shared_ptr<PSIO> psio(new PSIO());
        psio->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DFSCF_BJ, "length", (char*)&nQ_, sizeof(long int));
        psio->close(PSIF_DFSCF_BJ,1);
    }

    // restricted doubly occupied orbitals per irrep (optimized)
    rstcpi_   = (int*)malloc(nirrep_*sizeof(int));
    memset((void*)rstcpi_,'\0',nirrep_*sizeof(int));

    // restricted unoccupied occupied orbitals per irrep (optimized)
    rstvpi_   = (int*)malloc(nirrep_*sizeof(int));
    memset((void*)rstvpi_,'\0',nirrep_*sizeof(int));

    // active orbitals per irrep:
    amopi_    = (int*)malloc(nirrep_*sizeof(int));
    memset((void*)amopi_,'\0',nirrep_*sizeof(int));

    if (options["FROZEN_DOCC"].has_changed()) {
        if (options["FROZEN_DOCC"].size() != nirrep_) {
            throw PsiException("The FROZEN_DOCC array has the wrong dimensions_",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            frzcpi_[h] = options["FROZEN_DOCC"][h].to_double();
        }
    }
    if (options["RESTRICTED_DOCC"].has_changed()) {
        if (options["RESTRICTED_DOCC"].size() != nirrep_) {
            throw PsiException("The RESTRICTED_DOCC array has the wrong dimensions_",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            rstcpi_[h] = options["RESTRICTED_DOCC"][h].to_double();
        }
    }
    if (options["RESTRICTED_UOCC"].has_changed()) {
        if (options["RESTRICTED_UOCC"].size() != nirrep_) {
            throw PsiException("The RESTRICTED_UOCC array has the wrong dimensions_",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            rstvpi_[h] = options["RESTRICTED_UOCC"][h].to_double();
        }
    }
    if (options["FROZEN_UOCC"].has_changed()) {
        if (options["FROZEN_UOCC"].size() != nirrep_) {
            throw PsiException("The FROZEN_UOCC array has the wrong dimensions_",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            frzvpi_[h] = options["FROZEN_UOCC"][h].to_double();
        }
    }
    // user could specify active space with ACTIVE array
    if ( options["ACTIVE"].has_changed() ) {
        if (options["ACTIVE"].size() != nirrep_) {
            throw PsiException("The ACTIVE array has the wrong dimensions_",__FILE__,__LINE__);
        }

        // warn user that active array takes precedence over restricted_uocc array
        if (options["RESTRICTED_UOCC"].has_changed()) {
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
            amopi_[h]  = options["ACTIVE"][h].to_double();
            rstvpi_[h] = nmopi_[h] - frzcpi_[h] - rstcpi_[h] - frzvpi_[h] - amopi_[h];
        }
    }

    // now, total number of orbitals in each class:
    nfrzc_ = 0;
    nrstc_ = 0;
    amo_   = 0;
    nrstv_ = 0;
    nfrzv_ = 0;
    for (int h = 0; h < nirrep_; h++) {
        nfrzc_   += frzcpi_[h];
        nrstc_   += rstcpi_[h];
        nrstv_   += rstvpi_[h];
        nfrzv_   += frzvpi_[h];
        amo_     += nmopi_[h]-frzcpi_[h]-rstcpi_[h]-rstvpi_[h]-frzvpi_[h];
        amopi_[h] = nmopi_[h]-frzcpi_[h]-rstcpi_[h]-rstvpi_[h]-frzvpi_[h];
    }

    // lastly, determine orbital symmetries (in energy order)
    orbital_symmetries_ = (int*)malloc(nmo_*sizeof(int));
    memset((void*)orbital_symmetries_,'\0',nmo_*sizeof(int));

    bool * skip = (bool*)malloc(nmo_*sizeof(bool));
    memset((void*)skip,'\0',nmo_*sizeof(bool));

    // frozen core
    for (int i = 0; i < nfrzc_; i++){
        int me     = 0;
        double min = 1.0e99;
        int imin   = -999;
        int isym   = -999;
        for (int h = 0; h < nirrep_; h++) { 
            for (int j = 0; j < frzcpi_[h]; j++){
                if ( reference_wavefunction->epsilon_a()->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = reference_wavefunction->epsilon_a()->pointer(h)[j];
                        imin = me;
                        isym = h;
                    }
                }
                me++;
            }
            me += nmopi_[h] - frzcpi_[h] - frzvpi_[h];
        }
        if ( imin < 0 ) {
            throw PsiException("cannot determine energy-to-pitzer order (frzc)",__FILE__,__LINE__);
        }
        skip[imin] = true;
        orbital_symmetries_[i] = isym;
    }

    // restricted doubly occupied
    for (int i = nfrzc_; i < nrstc_ + nfrzc_; i++){

        int me     = 0;
        double min = 1.0e99;
        int imin   = -999;
        int isym   = -999;
        for (int h = 0; h < nirrep_; h++) {
            me += frzcpi_[h];
            for (int j = frzcpi_[h]; j < rstcpi_[h] + frzcpi_[h]; j++){
                if ( reference_wavefunction->epsilon_a()->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = reference_wavefunction->epsilon_a()->pointer(h)[j];
                        imin = me;
                        isym = h;
                    }
                }
                me++;
            }
            me += amopi_[h] + rstvpi_[h];// + frzvpi_[h];
        }
        if ( imin < 0 ) {
            throw PsiException("cannot determine energy-to-pitzer order (rstc)",__FILE__,__LINE__);
        }
        skip[imin] = true;
        orbital_symmetries_[i] = isym;
    }
    // active
    for (int i = nrstc_ + nfrzc_; i < amo_ + nrstc_ + nfrzc_; i++){

        int me     = 0;
        double min = 1.0e99;
        int imin   = -999;
        int isym   = -999;
        for (int h = 0; h < nirrep_; h++) {
            me += rstcpi_[h] + frzcpi_[h];
            for (int j = rstcpi_[h] + frzcpi_[h]; j < rstcpi_[h] + frzcpi_[h]+amopi_[h]; j++){
                if ( reference_wavefunction->epsilon_a()->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = reference_wavefunction->epsilon_a()->pointer(h)[j];
                        imin = me;
                        isym = h;
                    }
                }
                me++;
            }
            me += rstvpi_[h];// + frzvpi_[h];
        }
        if ( imin < 0 ) {
            throw PsiException("cannot determine energy-to-pitzer order (active)",__FILE__,__LINE__);
        }
        skip[imin] = true;
        orbital_symmetries_[i] = isym;
    }
    // restricted virtual
    for (int i = amo_ + nrstc_ + nfrzc_; i < nmo_ - nfrzv_ ; i++){

        int me     = 0;
        double min = 1.0e99;
        int imin   = -999;
        int isym   = -999;
        for (int h = 0; h < nirrep_; h++) {
            me += rstcpi_[h] + frzcpi_[h] + amopi_[h];
            for (int j = rstcpi_[h] + frzcpi_[h] + amopi_[h]; j < nmopi_[h] - frzvpi_[h] ; j++){
                if ( reference_wavefunction->epsilon_a()->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = reference_wavefunction->epsilon_a()->pointer(h)[j];
                        imin = me;
                        isym = h;
                    }
                }
                me++;
            }
            //me += frzvpi_[h];
        }
        if ( imin < 0 ) {
            throw PsiException("cannot determine energy-to-pitzer order (rstv)",__FILE__,__LINE__);
        }
        skip[imin] = true;
        orbital_symmetries_[i] = isym;
    }

    free(skip);

    // print orbitals per irrep in each space
    outfile->Printf("  ==> Active space details <==\n");
    outfile->Printf("\n");
    //outfile->Printf("        Freeze core orbitals?                   %5s\n",nfrzc_ > 0 ? "yes" : "no");
    outfile->Printf("        Number of frozen core orbitals:         %5i\n",nfrzc_);
    outfile->Printf("        Number of restricted occupied orbitals: %5i\n",nrstc_);
    outfile->Printf("        Number of active orbitals:              %5i\n",amo_);
    outfile->Printf("        Number of restricted virtual orbitals:  %5i\n",nrstv_);
    outfile->Printf("        Number of frozen virtual orbitals:      %5i\n",nfrzv_);
    outfile->Printf("\n");

    std::vector<std::string> labels = reference_wavefunction->molecule()->irrep_labels();
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

    // options
    // TODO: some options are missing:
    //    ORBOPT_ONE_STEP (should be handled outside of this class)
    //    ORBOPT_FREQUENCY (should be handled outside of this class)
    //    ORBOPT_NUM_DIIS_VECTORS (is diis implemented?)

    e_convergence_           = options.get_double("ORBOPT_ENERGY_CONVERGENCE");
    g_convergence_           = options.get_double("ORBOPT_GRADIENT_CONVERGENCE");
    active_active_rotations_ = options.get_bool("ORBOPT_ACTIVE_ACTIVE_ROTATIONS");
    maxiter_                 = options.get_int("ORBOPT_MAXITER");
    write_                   = options.get_bool("ORBOPT_WRITE");
    exact_diagonal_hessian_  = options.get_bool("ORBOPT_EXACT_DIAGONAL_HESSIAN");
    algorithm_               = options.get_str("ORBOPT_ALGORITHM");

    outfile->Printf("  ==> Orbital optimization parameters <==\n");
    outfile->Printf("\n");
    outfile->Printf("        g_convergence:               %12.3le\n",g_convergence_);
    outfile->Printf("        e_convergence:               %12.3le\n",e_convergence_);
    outfile->Printf("        maximum iterations:          %12i\n",maxiter_);
    outfile->Printf("        exact diagonal Hessian:      %12s\n",exact_diagonal_hessian_ ? "true" : "false");
    outfile->Printf("        print iteration info:        %12s\n",write_ ? "true" : "false");
    outfile->Printf("        algorithm:                   %12s\n",algorithm_.c_str());
    outfile->Printf("        active-active rotations:     %12s\n",active_active_rotations_ ? "true" : "false");
    outfile->Printf("\n");

}

OrbitalOptimizer::~OrbitalOptimizer(){

    free(rstcpi_);
    free(rstvpi_);
    free(amopi_);
    free(orbital_symmetries_);

}

void OrbitalOptimizer::get_hessian(std::shared_ptr<Matrix> Hessian){
    throw PsiException("get_hessian not yet implemented",__FILE__,__LINE__);
}

// GG functions

void OrbitalOptimizer::Determine_OindMap(){

    Qstride_ = nmo_ * ( nmo_ + 1 ) / 2;

    // Symmetry direct product table

    SymProd_[0*8+1] = SymProd_[1*8+0] = 1;
    SymProd_[0*8+2] = SymProd_[2*8+0] = 2;
    SymProd_[0*8+3] = SymProd_[3*8+0] = 3;
    SymProd_[0*8+4] = SymProd_[4*8+0] = 4;
    SymProd_[0*8+5] = SymProd_[5*8+0] = 5;
    SymProd_[0*8+6] = SymProd_[6*8+0] = 6;
    SymProd_[0*8+7] = SymProd_[7*8+0] = 7;
    SymProd_[1*8+2] = SymProd_[2*8+1] = 3;
    SymProd_[1*8+3] = SymProd_[3*8+1] = 2;
    SymProd_[1*8+4] = SymProd_[4*8+1] = 5;
    SymProd_[1*8+5] = SymProd_[5*8+1] = 4;
    SymProd_[1*8+6] = SymProd_[6*8+1] = 7;
    SymProd_[1*8+7] = SymProd_[7*8+1] = 6;
    SymProd_[2*8+3] = SymProd_[3*8+2] = 1;
    SymProd_[2*8+4] = SymProd_[4*8+2] = 6;
    SymProd_[2*8+5] = SymProd_[5*8+2] = 7;
    SymProd_[2*8+6] = SymProd_[6*8+2] = 4;
    SymProd_[2*8+7] = SymProd_[7*8+2] = 5;
    SymProd_[3*8+4] = SymProd_[4*8+3] = 7;
    SymProd_[3*8+5] = SymProd_[5*8+3] = 6;
    SymProd_[3*8+6] = SymProd_[6*8+3] = 5;
    SymProd_[3*8+7] = SymProd_[7*8+3] = 4;
    SymProd_[4*8+5] = SymProd_[5*8+4] = 1;
    SymProd_[4*8+6] = SymProd_[6*8+4] = 2;
    SymProd_[4*8+7] = SymProd_[7*8+4] = 3;
    SymProd_[5*8+6] = SymProd_[6*8+5] = 3;
    SymProd_[5*8+7] = SymProd_[7*8+5] = 2;
    SymProd_[6*8+7] = SymProd_[7*8+6] = 1;

    int * tmp_dim;
    tmp_dim = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)tmp_dim,'\0',nirrep_*sizeof(int));    

    // energy -> class & energy -> irrep & class -> irrep map

    int e2c_ind = 0;

    for (int h = 0; h < nirrep_; h++){

        for (int i=0;i<nrstc_;i++){

            if ( orbital_symmetries_[i] != h ) continue;

            OindMap_e2i_[i]       = tmp_dim[h];
            OindMap_e2c_[i]       = e2c_ind;
            OindMap_c2e_[e2c_ind] = i; 
            OindMap_c2i_[e2c_ind] = tmp_dim[h];

            tmp_dim[h]++;
            e2c_ind++;
               
        }

    }

    for (int h = 0; h < nirrep_; h++){

        for (int i=nrstc_; i < amo_ + nrstc_; i++){

            if ( orbital_symmetries_[i] != h ) continue;

            OindMap_e2i_[i] = tmp_dim[h];
            OindMap_e2c_[i] = e2c_ind;
            OindMap_c2e_[e2c_ind] = i;
            OindMap_c2i_[e2c_ind] = tmp_dim[h];

            tmp_dim[h]++;
            e2c_ind++;
               
        }

    }


    for (int h = 0;h < nirrep_; h++){

        for (int i=amo_+nrstc_; i < nrstv_ + amo_ + nrstc_; i++){

            if ( orbital_symmetries_[i] != h ) continue;

            OindMap_e2i_[i] = tmp_dim[h];
            OindMap_e2c_[i] = e2c_ind;
            OindMap_c2e_[e2c_ind] = i;
            OindMap_c2i_[e2c_ind] = tmp_dim[h];
              
            tmp_dim[h]++;
            e2c_ind++;

        }

    }

    // Geminal mappings

    memset((void*)tmp_dim,'\0',nirrep_*sizeof(int));

    for (int i = nrstc_; i < amo_ + nrstc_; i++){
 
        int h_i = orbital_symmetries_[i];

        for ( int j = nrstc_; j <= i; j++){

            int h_j = orbital_symmetries_[j];

            int h_ij = GemSym(h_i,h_j);

            Active_GindMap_e_[i * nmo_ + j] = tmp_dim[h_ij];
            Active_GindMap_e_[j * nmo_ + i] = tmp_dim[h_ij];

//            if ( h_ij == 0 ) printf("%i %i %i\n",i,j,tmp_dim[h_ij]);

            tmp_dim[h_ij]++;

        }

    }

    // Class index matrix

    int nmo_have = 0;

    for ( int h = 0; h < nirrep_; h++ ){
 
        first_index_[h][0] = nmo_have;
        last_index_[h][0]  = first_index_[h][0]+rstcpi_[h];
        nmo_have += rstcpi_[h];

    }
 
    for ( int h = 0; h < nirrep_; h++ ){

        first_index_[h][1] = nmo_have;
        last_index_[h][1]  = first_index_[h][1]+amopi_[h];
        nmo_have += amopi_[h];

    }

    for ( int h = 0; h < nirrep_; h++ ){

        first_index_[h][2] = nmo_have;
        last_index_[h][2]  = first_index_[h][2]+rstvpi_[h];
        nmo_have += rstvpi_[h];

    }

    int p_df = 0;

    for ( int h_p = 0; h_p < nirrep_; h_p++){

        for ( int p_class = 0; p_class < 3; p_class++){

            for ( int p_c = first_index_[h_p][p_class]; p_c < last_index_[h_p][p_class]; p_c++){

                OindMap_c2df_[p_c] = p_df;

                p_df++;

            }
     
        }
   
    }

    for ( int h = 1; h < nirrep_; h++){

        full_oe_offset_[h] = full_oe_offset_[h-1] + nmopi_[h-1] * nmopi_[h-1];
        lt_oe_offset_[h]   = lt_oe_offset_[h-1] + nmopi_[h-1] * ( nmopi_[h-1] + 1 ) / 2;    

    }

    for ( int h_ij = 0; h_ij < nirrep_; h_ij++){

        int ngem_tmp = 0;

        for ( int h_i=0; h_i<nirrep_; h_i++){

            int h_j = GemSym(h_ij,h_i); 
            if ( h_j > h_i ) continue;
            
            for ( int i_c = first_index_[h_i][1]; i_c < last_index_[h_i][1]; i_c++){

                int j_c_max = last_index_[h_j][1];
                if ( h_i == h_j ) j_c_max = i_c + 1;

                for ( int j_c = first_index_[h_j][1]; j_c < j_c_max; j_c++){

                    Active_GindMap_c_[i_c * nmo_ + j_c] = ngem_tmp;
                    Active_GindMap_c_[j_c * nmo_ + i_c] = ngem_tmp;                     

                    ngem_tmp++;

                }
 
            }

        }

    }

    memset((void*)tmp_dim,'\0',nirrep_*sizeof(int));

    for ( int p_c = 0; p_c < nmo_; p_c++){

        int h_p = orbital_symmetries_[OindMap_c2e_[p_c]];

        for ( int q_c = 0; q_c <= p_c ; q_c++){

            int h_q = orbital_symmetries_[OindMap_c2e_[q_c]];
            int h_pq = GemSym(h_p,h_q);

            Full_GindMap_c_[p_c * nmo_ + q_c] = tmp_dim[h_pq];
            Full_GindMap_c_[q_c * nmo_ + p_c] = tmp_dim[h_pq];

            tmp_dim[h_pq]++;

        }

    }

    // number ofgeminals

    for ( int h = 0; h< nirrep_; h++ ){

        for ( int h_i = 0; h_i < nirrep_; h_i++){

            int h_j = GemSym(h,h_i);

            if ( h_j < h_i ){

                doc_gempi_[h] += rstcpi_[h_i]*rstcpi_[h_j];
                act_gempi_[h] += amopi_[h_i]*amopi_[h_j];
                ext_gempi_[h] += rstvpi_[h_i]*rstvpi_[h_j];

            }

            else if ( h_j == h_i ){

                doc_gempi_[h] += rstcpi_[h_i]*(rstcpi_[h_i]+1)/2;
                act_gempi_[h] += amopi_[h_i]*(amopi_[h_i]+1)/2;
                ext_gempi_[h] += rstvpi_[h_i]*(rstvpi_[h_i]+1)/2;

            }

        }

        
    }

    for ( int h_tu =1; h_tu < nirrep_; h_tu++) d2_irrep_offset_[h_tu] = d2_irrep_offset_[h_tu-1] +\
                                                act_gempi_[h_tu-1] * ( act_gempi_[h_tu-1] + 1 ) / 2;

    nnz_Q_ = amopi_[0]*nmopi_[0];

    for ( int h = 1; h < nirrep_; h++){

        Q_offset_[h] = Q_offset_[h-1] + amopi_[h-1]*nmopi_[h-1];
        
        nnz_Q_ += amopi_[h]*nmopi_[h];
    }

    for ( int p_e = 0; p_e < nmo_; p_e++){

        int h_p = orbital_symmetries_[p_e];

        for ( int q_e = 0; q_e <= p_e; q_e++){

             int h_q = orbital_symmetries_[q_e];
             int h_pq = GemSym(h_p,h_q);

             full_gempi_[h_pq]++;

        }

    }

    free(tmp_dim); 

}

int OrbitalOptimizer::GemSym( int h_i, int h_j){
    return SymProd_[h_i*8 + h_j]; 
}

int OrbitalOptimizer::Full_ij_index( int row, int col, int nrow){
    return col * nrow + row;
}

int OrbitalOptimizer::Lt_ij_index( int row, int col){

    int index;

    if ( row > col) {

       index = row * ( row + 1 )/2 + col;        

    }
    else {   
 
       index = col * ( col + 1 )/2 + row;
 
    }

    return index;

}

int OrbitalOptimizer::Lt_ij_oe_index( int p_c, int q_c, int h){

    int p_i = OindMap_c2i_[p_c];
    int q_i = OindMap_c2i_[q_c];
   
    int index = Lt_ij_index(p_i, q_i) + lt_oe_offset_[h];

//    printf("input: %i %i %i  irrep: %i %i components: %i %i\n",h, p_c,q_c,p_i,q_i,Lt_ij_index(p_i, q_i),lt_oe_offset_[h]);

    return index;

}

int OrbitalOptimizer::Full_ij_oe_index( int row_c, int col_c, int h){

    int row_i = OindMap_c2i_[row_c];
    int col_i = OindMap_c2i_[col_c];

    int index = Full_ij_index(row_i, col_i,nmo_) + full_oe_offset_[h];

    return index;

}

int OrbitalOptimizer::df_ij_index( int p_c, int q_c ){

    int p_df = OindMap_c2df_[p_c];
    int q_df = OindMap_c2df_[q_c];

    int index = Lt_ij_index(p_df,q_df);

    return index;

}

int OrbitalOptimizer::Q_index( int t_c, int p_c, int h_p ){

    int t_red = t_c - first_index_[h_p][1]; 
    int p_i   = OindMap_c2i_[p_c] ;

    int index = Q_offset_[h_p] + p_i * amopi_[h_p] + t_red;

    return index;

}



int OrbitalOptimizer::Max_value( int * array, int dim){

    int maxval = 0;

    for ( int i = 0; i < dim; i++)

        if ( maxval < array[i] )
 
            maxval = array[i];

    return maxval;
}



void OrbitalOptimizer::Alloc_OindMap(){

    int nmo = nrstc_ + amo_ + nrstv_;

    SymProd_ = (int*)malloc(64*sizeof(int));
    memset((void*)SymProd_,'\0',64*sizeof(int));

    doc_gempi_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)doc_gempi_,'\0',nirrep_*sizeof(int));

    act_gempi_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)act_gempi_,'\0',nirrep_*sizeof(int));       

    ext_gempi_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)ext_gempi_,'\0',nirrep_*sizeof(int));

    full_gempi_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)full_gempi_,'\0',nirrep_*sizeof(int));

    lt_oe_offset_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)lt_oe_offset_,'\0',nirrep_*sizeof(int));

    full_oe_offset_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)full_oe_offset_,'\0',nirrep_*sizeof(int));    

    d2_irrep_offset_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)d2_irrep_offset_,'\0',nirrep_*sizeof(int));

    Q_offset_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)Q_offset_,'\0',nirrep_*sizeof(int));

    aa_ras1pi_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)aa_ras1pi_,'\0',nirrep_*sizeof(int));

    OindMap_e2c_ = (int *)malloc(nmo*sizeof(int)); 
    memset((void*)OindMap_e2c_,'\0',nmo*sizeof(int));

    OindMap_c2e_ = (int *)malloc(nmo*sizeof(int));
    memset((void*)OindMap_c2e_,'\0',nmo*sizeof(int));

    OindMap_e2i_ = (int *)malloc(nmo*sizeof(int));
    memset((void*)OindMap_e2i_,'\0',nmo*sizeof(int));

    OindMap_c2i_ = (int *)malloc(nmo*sizeof(int));
    memset((void*)OindMap_c2i_,'\0',nmo*sizeof(int));

    OindMap_c2df_ = (int *)malloc(nmo*sizeof(int));
    memset((void*)OindMap_c2df_,'\0',nmo*sizeof(int));

    Active_GindMap_e_ = (int *)malloc(nmo*nmo*sizeof(int));
    memset((void*)Active_GindMap_e_,'\0',nmo*nmo*sizeof(int));

    Active_GindMap_c_ = (int *)malloc(nmo*nmo*sizeof(int));
    memset((void*)Active_GindMap_c_,'\0',nmo*nmo*sizeof(int));

    Full_GindMap_c_ = (int *)malloc(nmo*nmo*sizeof(int));
    memset((void*)Full_GindMap_c_,'\0',nmo*nmo*sizeof(int));

    first_index_ = (int**)malloc(nirrep_*sizeof(int*));
    for (int i = 0; i < nirrep_; i++){
        first_index_[i] = (int*)malloc(3*sizeof(int));
    }

    last_index_ = (int**)malloc(nirrep_*sizeof(int*));
    for (int i = 0; i < nirrep_; i++){
        last_index_[i] = (int*)malloc(3*sizeof(int));
    }
};

void OrbitalOptimizer::Dealloc_OindMap(){

    free(SymProd_);
    free(doc_gempi_);
    free(act_gempi_);
    free(ext_gempi_);
    free(d2_irrep_offset_);
    free(full_gempi_);
    free(full_oe_offset_);
    free(Q_offset_);
    free(lt_oe_offset_);
    free(OindMap_e2c_);           
    free(OindMap_c2e_);
    free(OindMap_e2i_);
    free(OindMap_c2i_);
    free(OindMap_c2df_);
    free(Active_GindMap_c_);
    free(Active_GindMap_e_);
    free(Full_GindMap_c_);
    free(aa_ras1pi_); 

    for (int i = 0; i < nirrep_; i++){
        free(first_index_[i]);
    }
    free(first_index_);

    for (int i = 0; i < nirrep_; i++){
        free(last_index_[i]);
    }
    free(last_index_);


};

}// end of namespace
