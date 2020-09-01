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

#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libqt/qt.h>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libtrans/mospace.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/matrix.h>

#include"v2rdm_solver.h"

#include <misc/omp.h>

using namespace psi;

namespace hilbert{

void v2RDMSolver::build_mapping_arrays() {

    // product table:
    table = (int*)malloc(64*sizeof(int));
    memset((void*)table,'\0',64*sizeof(int));
    table[0*8+1] = table[1*8+0] = 1;
    table[0*8+2] = table[2*8+0] = 2;
    table[0*8+3] = table[3*8+0] = 3;
    table[0*8+4] = table[4*8+0] = 4;
    table[0*8+5] = table[5*8+0] = 5;
    table[0*8+6] = table[6*8+0] = 6;
    table[0*8+7] = table[7*8+0] = 7;
    table[1*8+2] = table[2*8+1] = 3;
    table[1*8+3] = table[3*8+1] = 2;
    table[1*8+4] = table[4*8+1] = 5;
    table[1*8+5] = table[5*8+1] = 4;
    table[1*8+6] = table[6*8+1] = 7;
    table[1*8+7] = table[7*8+1] = 6;
    table[2*8+3] = table[3*8+2] = 1;
    table[2*8+4] = table[4*8+2] = 6;
    table[2*8+5] = table[5*8+2] = 7;
    table[2*8+6] = table[6*8+2] = 4;
    table[2*8+7] = table[7*8+2] = 5;
    table[3*8+4] = table[4*8+3] = 7;
    table[3*8+5] = table[5*8+3] = 6;
    table[3*8+6] = table[6*8+3] = 5;
    table[3*8+7] = table[7*8+3] = 4;
    table[4*8+5] = table[5*8+4] = 1;
    table[4*8+6] = table[6*8+4] = 2;
    table[4*8+7] = table[7*8+4] = 3;
    table[5*8+6] = table[6*8+5] = 3;
    table[5*8+7] = table[7*8+5] = 2;
    table[6*8+7] = table[7*8+6] = 1;

    // orbitals are in pitzer order:
    symmetry               = (int*)malloc(nmo_*sizeof(int));
    symmetry_full          = (int*)malloc((nmo_-nfrzv_)*sizeof(int));
    //gg-fc modified dimension to include frozen core orbitals
    symmetry_energy_order  = (int*)malloc((nmo_-nfrzv_+nfrzc_)*sizeof(int));
    energy_to_pitzer_order = (int*)malloc((nmo_-nfrzv_)*sizeof(int));
    full_basis             = (int*)malloc((nmo_-nfrzv_)*sizeof(int));

    // including frozen virtuals
    symmetry_really_full               = (int*)malloc(nmo_*sizeof(int));
    energy_to_pitzer_order_really_full = (int*)malloc(nmo_*sizeof(int));

    memset((void*)symmetry,'\0',nmo_*sizeof(int));
    memset((void*)symmetry_full,'\0',(nmo_-nfrzv_)*sizeof(int));
    //gg-fc modified dimension to include frozen core orbitals
    memset((void*)symmetry_energy_order,'\0',(nmo_-nfrzv_+nfrzc_)*sizeof(int));
    memset((void*)energy_to_pitzer_order,'\0',(nmo_-nfrzv_)*sizeof(int));
    memset((void*)full_basis,'\0',(nmo_-nfrzv_)*sizeof(int));

    memset((void*)symmetry_really_full,'\0',nmo_*sizeof(int));
    memset((void*)energy_to_pitzer_order_really_full,'\0',nmo_*sizeof(int));

    int count = 0;
    int count_full = 0;

    // symmetry of ACTIVE orbitals
    for (int h = 0; h < nirrep_; h++) {
        count_full += rstcpi_[h] + frzcpi_[h];
        for (int norb = rstcpi_[h] + frzcpi_[h]; norb < nmopi_[h] - rstvpi_[h] - frzvpi_[h]; norb++){
            full_basis[count] = count_full++;
            symmetry[count++] = h;
        }
        count_full += rstvpi_[h]; // + frzvpi_[h];
    }

    // symmetry of all orbitals, except frozen virtuals
    count = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int norb = 0; norb < nmopi_[h] - frzvpi_[h]; norb++){
            symmetry_full[count++] = h;
        }
    }
    // ok, including frozen virtuals
    count = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int norb = 0; norb < nmopi_[h]; norb++){
            symmetry_really_full[count++] = h;
        }
    }

    int * skip = (int*)malloc(nmo_*sizeof(int));
    memset((void*)skip,'\0',nmo_*sizeof(int));

    // warning to future eugene:  it is possible that
    // this ordering will differ from that printed at the
    // end of the SCF routine if orbitals are truly
    // degenerate.  past eugene hasn't convinved himself
    // of whether or not this is actually a problem.

    // hey, past eugene! it turns out you were on to something!
    // when restarting jobs, if the SCF has degenerate orbitals,
    // sometimes their order can change.  How annoying!

    // TODO: the orbital ordering should be according to 
    // energy within each type of orbital

    // frozen core
    for (int i = 0; i < nfrzc_; i++){
        int me     = 0;
        double min = 1.0e99;
        int imin   = -999;
        int isym   = -999;
        for (int h = 0; h < nirrep_; h++) {
            for (int j = 0; j < frzcpi_[h]; j++){
                if ( epsilon_a_->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = epsilon_a_->pointer(h)[j];
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
        skip[imin] = 1;
        symmetry_energy_order[i] = isym + 1;
        energy_to_pitzer_order[i] = imin;
    }
    // active
    for (int i = nfrzc_; i < nrstc_ + nfrzc_; i++){

        int me     = 0;
        double min = 1.0e99;
        int imin   = -999;
        int isym   = -999;
        for (int h = 0; h < nirrep_; h++) {
            me += frzcpi_[h];
            for (int j = frzcpi_[h]; j < rstcpi_[h] + frzcpi_[h]; j++){
                if ( epsilon_a_->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = epsilon_a_->pointer(h)[j];
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
        skip[imin] = 1;
        symmetry_energy_order[i] = isym + 1;
        energy_to_pitzer_order[i] = imin;
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
                if ( epsilon_a_->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = epsilon_a_->pointer(h)[j];
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
        skip[imin] = 1;
        symmetry_energy_order[i] = isym + 1;
        energy_to_pitzer_order[i] = imin;
    }
    // restricted virtual
    //gg-fc modified to include frozen core orbitals
    for (int i = amo_ + nrstc_ + nfrzc_; i < nmo_ - nfrzv_ ; i++){

        int me     = 0;
        double min = 1.0e99;
        int imin   = -999;
        int isym   = -999;
        for (int h = 0; h < nirrep_; h++) {
            me += rstcpi_[h] + frzcpi_[h] + amopi_[h];
            for (int j = rstcpi_[h] + frzcpi_[h] + amopi_[h]; j < nmopi_[h] - frzvpi_[h] ; j++){
                if ( epsilon_a_->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = epsilon_a_->pointer(h)[j];
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
        skip[imin] = 1;
        symmetry_energy_order[i] = isym + 1;
        energy_to_pitzer_order[i] = imin;
    }
    // frozen virtual
    //for (int i = amo_ + nrstc_ + nfrzc_ + nrstv_; i < nmo_; i++){
    //    int me = 0;
    //    min = 1.0e99;
    //    for (int h = 0; h < nirrep_; h++) {
    //        me += rstcpi_[h] + frzcpi_[h] + amopi_[h] + rstvpi_[h];
    //        for (int j = rstcpi_[h] + frzcpi_[h] + amopi_[h] + rstvpi_[h]; j < nmopi_[h]; j++){
    //            if ( epsilon_a_->pointer(h)[j] < min ) {
    //                if ( !skip[me] ) {
    //                    min = epsilon_a_->pointer(h)[j];
    //                    imin = me;
    //                    isym = h;
    //                }
    //            }
    //            me++;
    //        }
    //    }
    //    skip[imin] = 1;
    //    symmetry_energy_order[i] = isym + 1;
    //    energy_to_pitzer_order[i] = imin;
    //}

    // again, including frozen virtuals:
    memset((void*)skip,'\0',nmo_*sizeof(int));

    // frozen core
    for (int i = 0; i < nfrzc_; i++){
        int me     = 0;
        double min = 1.0e99;
        int imin   = -999;
        int isym   = -999;
        for (int h = 0; h < nirrep_; h++) {
            for (int j = 0; j < frzcpi_[h]; j++){
                if ( epsilon_a_->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = epsilon_a_->pointer(h)[j];
                        imin = me;
                        isym = h;
                    }
                }
                me++;
            }
            me += nmopi_[h] - frzcpi_[h];// - frzvpi_[h];
        }
        if ( imin < 0 ) {
            throw PsiException("cannot determine energy-to-pitzer order (frzc)",__FILE__,__LINE__);
        }
        skip[imin] = 1;
        energy_to_pitzer_order_really_full[i] = imin;
    }
    // active
    for (int i = nfrzc_; i < nrstc_ + nfrzc_; i++){

        int me     = 0;
        double min = 1.0e99;
        int imin   = -999;
        int isym   = -999;
        for (int h = 0; h < nirrep_; h++) {
            me += frzcpi_[h];
            for (int j = frzcpi_[h]; j < rstcpi_[h] + frzcpi_[h]; j++){
                if ( epsilon_a_->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = epsilon_a_->pointer(h)[j];
                        imin = me;
                        isym = h;
                    }
                }
                me++;
            }
            me += amopi_[h] + rstvpi_[h] + frzvpi_[h];
        }
        if ( imin < 0 ) {
            throw PsiException("cannot determine energy-to-pitzer order (rstc)",__FILE__,__LINE__);
        }
        skip[imin] = 1;
        energy_to_pitzer_order_really_full[i] = imin;
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
                if ( epsilon_a_->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = epsilon_a_->pointer(h)[j];
                        imin = me;
                        isym = h;
                    }
                }
                me++;
            }
            me += rstvpi_[h] + frzvpi_[h];
        }
        if ( imin < 0 ) {
            throw PsiException("cannot determine energy-to-pitzer order (active)",__FILE__,__LINE__);
        }
        skip[imin] = 1;
        energy_to_pitzer_order_really_full[i] = imin;
    }
    // restricted virtual
    for (int i = amo_ + nrstc_ + nfrzc_; i < nmo_ - nfrzv_; i++){

        int me     = 0;
        double min = 1.0e99;
        int imin   = -999;
        int isym   = -999;
        for (int h = 0; h < nirrep_; h++) {
            me += rstcpi_[h] + frzcpi_[h] + amopi_[h];
            for (int j = rstcpi_[h] + frzcpi_[h] + amopi_[h]; j < nmopi_[h] - frzvpi_[h]; j++){
                if ( epsilon_a_->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = epsilon_a_->pointer(h)[j];
                        imin = me;
                        isym = h;
                    }
                }
                me++;
            }
            me += frzvpi_[h];
        }
        if ( imin < 0 ) {
            throw PsiException("cannot determine energy-to-pitzer order (rstv)",__FILE__,__LINE__);
        }
        skip[imin] = 1;
        energy_to_pitzer_order_really_full[i] = imin;
    }
    // frozen virtual
    for (int i = amo_ + nrstc_ + nfrzc_ + nrstv_; i < nmo_; i++){
        int me = 0;
        double min = 1.0e99;
        int imin = -999;
        int isym = -999;
        for (int h = 0; h < nirrep_; h++) {
            me += rstcpi_[h] + frzcpi_[h] + amopi_[h] + rstvpi_[h];
            for (int j = rstcpi_[h] + frzcpi_[h] + amopi_[h] + rstvpi_[h]; j < nmopi_[h]; j++){
                if ( epsilon_a_->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = epsilon_a_->pointer(h)[j];
                        imin = me;
                        isym = h;
                    }
                }
                me++;
            }
        }
        if ( imin < 0 ) {
            throw PsiException("cannot determine energy-to-pitzer order (frzv)",__FILE__,__LINE__);
        }
        skip[imin] = 1;
        energy_to_pitzer_order_really_full[i] = imin;
    }
    free(skip);

    pitzer_offset           = (int*)malloc(nirrep_*sizeof(int));
    pitzer_offset_full      = (int*)malloc(nirrep_*sizeof(int));
    count = 0;
    for (int h = 0; h < nirrep_; h++) {
        pitzer_offset[h] = count;
        count += amopi_[h]; 
    }
    count = 0;
    for (int h = 0; h < nirrep_; h++) {
        pitzer_offset_full[h] = count;
        count += nmopi_[h] - frzvpi_[h];
    }

    // geminals, by symmetry
    for (int h = 0; h < nirrep_; h++) {
        std::vector < std::pair<int,int> > mygems;
        for (int i = 0; i < amo_; i++) {
            for (int j = 0; j < amo_; j++) {
                int sym = SymmetryPair(symmetry[i],symmetry[j]);
                if (h==sym) {
                    mygems.push_back(std::make_pair(j,i));
                }

            }
        }
        gems.push_back(mygems);
    }

    bas_ab_sym           = (int***)malloc(nirrep_*sizeof(int**));
    bas_aa_sym           = (int***)malloc(nirrep_*sizeof(int**));
    bas_00_sym           = (int***)malloc(nirrep_*sizeof(int**));
    bas_full_sym         = (int***)malloc(nirrep_*sizeof(int**));
    bas_really_full_sym  = (int***)malloc(nirrep_*sizeof(int**));

    ibas_ab_sym          = (int***)malloc(nirrep_*sizeof(int**));
    ibas_aa_sym          = (int***)malloc(nirrep_*sizeof(int**));
    ibas_00_sym          = (int***)malloc(nirrep_*sizeof(int**));
    ibas_full_sym        = (int***)malloc(nirrep_*sizeof(int**));
    ibas_really_full_sym = (int***)malloc(nirrep_*sizeof(int**));

    gems_ab              = (int*)malloc(nirrep_*sizeof(int));
    gems_aa              = (int*)malloc(nirrep_*sizeof(int));
    gems_00              = (int*)malloc(nirrep_*sizeof(int));
    gems_full            = (int*)malloc(nirrep_*sizeof(int));
    gems_plus_core       = (int*)malloc(nirrep_*sizeof(int));

    for (int h = 0; h < nirrep_; h++) {

        ibas_ab_sym[h]          = (int**)malloc(amo_*sizeof(int*));
        ibas_aa_sym[h]          = (int**)malloc(amo_*sizeof(int*));
        ibas_00_sym[h]          = (int**)malloc(amo_*sizeof(int*));
        ibas_full_sym[h]        = (int**)malloc(nmo_*sizeof(int*));
        ibas_really_full_sym[h] = (int**)malloc(nmo_*sizeof(int*));

        bas_ab_sym[h]         = (int**)malloc(amo_*amo_*sizeof(int*));
        bas_aa_sym[h]         = (int**)malloc(amo_*amo_*sizeof(int*));
        bas_00_sym[h]         = (int**)malloc(amo_*amo_*sizeof(int*));
        bas_full_sym[h]       = (int**)malloc(nmo_*nmo_*sizeof(int*));
        bas_really_full_sym[h]       = (int**)malloc(nmo_*nmo_*sizeof(int*));

        // active space geminals
        for (int i = 0; i < amo_; i++) {
            ibas_ab_sym[h][i] = (int*)malloc(amo_*sizeof(int));
            ibas_aa_sym[h][i] = (int*)malloc(amo_*sizeof(int));
            ibas_00_sym[h][i] = (int*)malloc(amo_*sizeof(int));
            for (int j = 0; j < amo_; j++) {
                ibas_ab_sym[h][i][j] = -999;
                ibas_aa_sym[h][i][j] = -999;
                ibas_00_sym[h][i][j] = -999;
            }
        }
        for (int i = 0; i < amo_*amo_; i++) {
            bas_ab_sym[h][i] = (int*)malloc(2*sizeof(int));
            bas_aa_sym[h][i] = (int*)malloc(2*sizeof(int));
            bas_00_sym[h][i] = (int*)malloc(2*sizeof(int));
            for (int j = 0; j < 2; j++) {
                bas_ab_sym[h][i][j] = -999;
                bas_aa_sym[h][i][j] = -999;
                bas_00_sym[h][i][j] = -999;
            }
        }
        // full space geminals
        for (int i = 0; i < nmo_; i++) {
            ibas_full_sym[h][i]        = (int*)malloc(nmo_*sizeof(int));
            ibas_really_full_sym[h][i] = (int*)malloc(nmo_*sizeof(int));
            for (int j = 0; j < nmo_; j++) {
                ibas_full_sym[h][i][j]        = -999;
                ibas_really_full_sym[h][i][j] = -999;
            }
        }
        for (int i = 0; i < nmo_*nmo_; i++) {
            bas_full_sym[h][i]        = (int*)malloc(2*sizeof(int));
            bas_really_full_sym[h][i] = (int*)malloc(2*sizeof(int));
            for (int j = 0; j < 2; j++) {
                bas_full_sym[h][i][j]        = -999;
                bas_really_full_sym[h][i][j] = -999;
            }
        }

        // active space mappings:
        int count_ab = 0;
        int count_aa = 0;
        for (int n = 0; n < gems[h].size(); n++) {
            int i = gems[h][n].first;
            int j = gems[h][n].second;

            ibas_ab_sym[h][i][j] = n;
            bas_ab_sym[h][n][0]  = i;
            bas_ab_sym[h][n][1]  = j;
            count_ab++;

            if ( i >= j ) continue;

            ibas_aa_sym[h][i][j] = count_aa;
            ibas_aa_sym[h][j][i] = count_aa;
            bas_aa_sym[h][count_aa][0] = i;
            bas_aa_sym[h][count_aa][1] = j;
            count_aa++;
        }
        gems_ab[h] = count_ab;
        gems_aa[h] = count_aa;

    }

    // if restarting a job, need to read a few energy-order arrays from disk...
    /*if ( options_["RESTART_FROM_CHECKPOINT_FILE"].has_changed() ) {
        std::shared_ptr<PSIO> psio (new PSIO() );
        psio->open(PSIF_V2RDM_CHECKPOINT,PSIO_OPEN_OLD);

        // energy order to pitzer order mapping array
        psio->read_entry(PSIF_V2RDM_CHECKPOINT,"ENERGY_TO_PITZER_ORDER",
            (char*)energy_to_pitzer_order,(nmo_-nfrzv_)*sizeof(int));

        // energy order to pitzer order mapping array
        psio->read_entry(PSIF_V2RDM_CHECKPOINT,"ENERGY_TO_PITZER_ORDER_REALLY_FULL",
            (char*)energy_to_pitzer_order_really_full,nmo_*sizeof(int)); 

        psio->close(PSIF_V2RDM_CHECKPOINT,1);
    }*/

    // new way:
    memset((void*)gems_full,'\0',nirrep_*sizeof(int));
    memset((void*)gems_plus_core,'\0',nirrep_*sizeof(int));

    // everything except frozen virtuals
    for (int ieo = 0; ieo < nmo_ - nfrzv_; ieo++) {
        int ifull = energy_to_pitzer_order[ieo];
        int hi    = symmetry_full[ifull];
        //int i     = ifull - pitzer_offset_full[hi];
        for (int jeo = 0; jeo <= ieo; jeo++) {
            int jfull = energy_to_pitzer_order[jeo];
            int hj    = symmetry_full[jfull];
            //int j     = jfull - pitzer_offset_full[hj];

            int hij = SymmetryPair(hi,hj);
            ibas_full_sym[hij][ifull][jfull] = gems_full[hij];
            ibas_full_sym[hij][jfull][ifull] = gems_full[hij];
            bas_full_sym[hij][gems_full[hij]][0] = ifull;
            bas_full_sym[hij][gems_full[hij]][1] = jfull;
            gems_full[hij]++;
            if ( ieo < amo_ + nrstc_ + nfrzc_ && jeo < amo_ + nrstc_ + nfrzc_ ) {
                gems_plus_core[hij]++;
            }
        }
    }
    // including frozen virtuals
    int * gems_really_full = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)gems_really_full,'\0',nirrep_*sizeof(int));
    for (int ieo = 0; ieo < nmo_; ieo++) {
        int ifull = energy_to_pitzer_order_really_full[ieo];
        int hi    = symmetry_really_full[ifull];
        //int i     = ifull - pitzer_offset_full[hi];
        for (int jeo = 0; jeo <= ieo; jeo++) {
            int jfull = energy_to_pitzer_order_really_full[jeo];
            int hj    = symmetry_really_full[jfull];
            //int j     = jfull - pitzer_offset_full[hj];

            int hij = SymmetryPair(hi,hj);
            ibas_really_full_sym[hij][ifull][jfull] = gems_really_full[hij];
            ibas_really_full_sym[hij][jfull][ifull] = gems_really_full[hij];
            bas_really_full_sym[hij][gems_really_full[hij]][0] = ifull;
            bas_really_full_sym[hij][gems_really_full[hij]][1] = jfull;
            gems_really_full[hij]++;
        }
    }
    // active only
    int * pitzer_full_to_active = (int*)malloc(nmo_*sizeof(int));
    int off_full = 0;
    int off_act  = 0;
    for (int i = 0; i < nmo_; i++) {
        pitzer_full_to_active[i] = -999;
    }
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            pitzer_full_to_active[off_full + i + frzcpi_[h] + rstcpi_[h]] = off_act + i;
        }
        off_full += nmopi_[h] + frzvpi_[h];
        off_act  += amopi_[h];
    }

    memset((void*)gems_00,'\0',nirrep_*sizeof(int));
    for (int ieo = 0; ieo < nmo_ - nfrzv_; ieo++) {

        int ifull = energy_to_pitzer_order[ieo];
        int hi    = symmetry_full[ifull];
        int iact  = pitzer_full_to_active[ifull];
        if ( iact == -999 ) continue;

        for (int jeo = 0; jeo <= ieo; jeo++) {

            int jfull = energy_to_pitzer_order[jeo];
            int hj    = symmetry_full[jfull];
            int jact  = pitzer_full_to_active[jfull];
            if ( jact == -999 ) continue;

            int hij   = SymmetryPair(hi,hj);

            ibas_00_sym[hij][iact][jact] = gems_00[hij];
            ibas_00_sym[hij][jact][iact] = gems_00[hij];

            bas_00_sym[hij][gems_00[hij]][0] = iact;
            bas_00_sym[hij][gems_00[hij]][1] = jact;

            gems_00[hij]++;

        }
    }
    free(pitzer_full_to_active);

    if ( constrain_t1_ || constrain_t2_ || constrain_d3_ ) {
        // make all triplets
        for (int h = 0; h < nirrep_; h++) {
            std::vector < std::tuple<int,int,int> > mytrip;
            for (int i = 0; i < amo_; i++) {
                for (int j = 0; j < amo_; j++) {
                    int s1 = SymmetryPair(symmetry[i],symmetry[j]);
                    for (int k = 0; k < amo_; k++) {
                        int s2 = SymmetryPair(s1,symmetry[k]);
                        if (h==s2) {
                            mytrip.push_back(std::make_tuple(i,j,k));
                        }
                    }

                }
            }
            triplets.push_back(mytrip);
        }
        bas_aaa_sym  = (int***)malloc(nirrep_*sizeof(int**));
        bas_aab_sym  = (int***)malloc(nirrep_*sizeof(int**));
        bas_aba_sym  = (int***)malloc(nirrep_*sizeof(int**));
        ibas_aaa_sym = (int****)malloc(nirrep_*sizeof(int***));
        ibas_aab_sym = (int****)malloc(nirrep_*sizeof(int***));
        ibas_aba_sym = (int****)malloc(nirrep_*sizeof(int***));
        trip_aaa    = (int*)malloc(nirrep_*sizeof(int));
        trip_aab    = (int*)malloc(nirrep_*sizeof(int));
        trip_aba    = (int*)malloc(nirrep_*sizeof(int));
        for (int h = 0; h < nirrep_; h++) {
            ibas_aaa_sym[h] = (int***)malloc(amo_*sizeof(int**));
            ibas_aab_sym[h] = (int***)malloc(amo_*sizeof(int**));
            ibas_aba_sym[h] = (int***)malloc(amo_*sizeof(int**));
            bas_aaa_sym[h]  = (int**)malloc(amo_*amo_*amo_*sizeof(int*));
            bas_aab_sym[h]  = (int**)malloc(amo_*amo_*amo_*sizeof(int*));
            bas_aba_sym[h]  = (int**)malloc(amo_*amo_*amo_*sizeof(int*));
            for (int i = 0; i < amo_; i++) {
                ibas_aaa_sym[h][i] = (int**)malloc(amo_*sizeof(int*));
                ibas_aab_sym[h][i] = (int**)malloc(amo_*sizeof(int*));
                ibas_aba_sym[h][i] = (int**)malloc(amo_*sizeof(int*));
                for (int j = 0; j < amo_; j++) {
                    ibas_aaa_sym[h][i][j] = (int*)malloc(amo_*sizeof(int));
                    ibas_aab_sym[h][i][j] = (int*)malloc(amo_*sizeof(int));
                    ibas_aba_sym[h][i][j] = (int*)malloc(amo_*sizeof(int));
                    for (int k = 0; k < amo_; k++) {
                        ibas_aaa_sym[h][i][j][k] = -999;
                        ibas_aab_sym[h][i][j][k] = -999;
                        ibas_aba_sym[h][i][j][k] = -999;
                    }
                }
            }
            for (int i = 0; i < amo_*amo_*amo_; i++) {
                bas_aaa_sym[h][i] = (int*)malloc(3*sizeof(int));
                bas_aab_sym[h][i] = (int*)malloc(3*sizeof(int));
                bas_aba_sym[h][i] = (int*)malloc(3*sizeof(int));
                for (int j = 0; j < 3; j++) {
                    bas_aaa_sym[h][i][j] = -999;
                    bas_aab_sym[h][i][j] = -999;
                    bas_aba_sym[h][i][j] = -999;
                }
            }

            // mappings:
            int count_aaa = 0;
            int count_aab = 0;
            int count_aba = 0;
            for (int n = 0; n < triplets[h].size(); n++) {
                int i = std::get<0>(triplets[h][n]);
                int j = std::get<1>(triplets[h][n]);
                int k = std::get<2>(triplets[h][n]);

                ibas_aba_sym[h][i][j][k] = count_aba;
                bas_aba_sym[h][count_aba][0]  = i;
                bas_aba_sym[h][count_aba][1]  = j;
                bas_aba_sym[h][count_aba][2]  = k;
                count_aba++;

                if ( i >= j ) continue;

                ibas_aab_sym[h][i][j][k] = count_aab;
                ibas_aab_sym[h][j][i][k] = count_aab;
                bas_aab_sym[h][count_aab][0]  = i;
                bas_aab_sym[h][count_aab][1]  = j;
                bas_aab_sym[h][count_aab][2]  = k;
                count_aab++;

                if ( j >= k ) continue;

                ibas_aaa_sym[h][i][j][k] = count_aaa;
                ibas_aaa_sym[h][i][k][j] = count_aaa;
                ibas_aaa_sym[h][j][i][k] = count_aaa;
                ibas_aaa_sym[h][j][k][i] = count_aaa;
                ibas_aaa_sym[h][k][i][j] = count_aaa;
                ibas_aaa_sym[h][k][j][i] = count_aaa;
                bas_aaa_sym[h][count_aaa][0]  = i;
                bas_aaa_sym[h][count_aaa][1]  = j;
                bas_aaa_sym[h][count_aaa][2]  = k;
                count_aaa++;

            }
            trip_aaa[h] = count_aaa;
            trip_aab[h] = count_aab;
            trip_aba[h] = count_aba;
        }
    }
    if ( constrain_d4_ ) {
        // make all quartets
        for (int h = 0; h < nirrep_; h++) {
            std::vector < std::tuple<int,int,int,int> > myquartet;
            for (int i = 0; i < amo_; i++) {
                for (int j = 0; j < amo_; j++) {
                    int s1 = SymmetryPair(symmetry[i],symmetry[j]);
                    for (int k = 0; k < amo_; k++) {
                        int s2 = SymmetryPair(s1,symmetry[k]);
                        for (int l = 0; l < amo_; l++) {
                            int s3 = SymmetryPair(s2,symmetry[l]);
                            if (h==s3) {
                                myquartet.push_back(std::make_tuple(i,j,k,l));
                            }
                        }
                    }
                }
            }
            quartets.push_back(myquartet);
        }
        bas_aaaa_sym  = (int***)malloc(nirrep_*sizeof(int**));
        bas_aaab_sym  = (int***)malloc(nirrep_*sizeof(int**));
        bas_aabb_sym  = (int***)malloc(nirrep_*sizeof(int**));
        ibas_aaaa_sym = (int*****)malloc(nirrep_*sizeof(int****));
        ibas_aaab_sym = (int*****)malloc(nirrep_*sizeof(int****));
        ibas_aabb_sym = (int*****)malloc(nirrep_*sizeof(int****));
        quartet_aaaa  = (int*)malloc(nirrep_*sizeof(int));
        quartet_aaab  = (int*)malloc(nirrep_*sizeof(int));
        quartet_aabb  = (int*)malloc(nirrep_*sizeof(int));
        for (int h = 0; h < nirrep_; h++) {
            ibas_aaaa_sym[h] = (int****)malloc(amo_*sizeof(int***));
            ibas_aaab_sym[h] = (int****)malloc(amo_*sizeof(int***));
            ibas_aabb_sym[h] = (int****)malloc(amo_*sizeof(int***));
            bas_aaaa_sym[h]  = (int**)malloc(amo_*amo_*amo_*amo_*sizeof(int*));
            bas_aaab_sym[h]  = (int**)malloc(amo_*amo_*amo_*amo_*sizeof(int*));
            bas_aabb_sym[h]  = (int**)malloc(amo_*amo_*amo_*amo_*sizeof(int*));
            for (int i = 0; i < amo_; i++) {
                ibas_aaaa_sym[h][i] = (int***)malloc(amo_*sizeof(int**));
                ibas_aaab_sym[h][i] = (int***)malloc(amo_*sizeof(int**));
                ibas_aabb_sym[h][i] = (int***)malloc(amo_*sizeof(int**));
                for (int j = 0; j < amo_; j++) {
                    ibas_aaaa_sym[h][i][j] = (int**)malloc(amo_*sizeof(int*));
                    ibas_aaab_sym[h][i][j] = (int**)malloc(amo_*sizeof(int*));
                    ibas_aabb_sym[h][i][j] = (int**)malloc(amo_*sizeof(int*));
                    for (int k = 0; k < amo_; k++) {
                        ibas_aaaa_sym[h][i][j][k] = (int*)malloc(amo_*sizeof(int));
                        ibas_aaab_sym[h][i][j][k] = (int*)malloc(amo_*sizeof(int));
                        ibas_aabb_sym[h][i][j][k] = (int*)malloc(amo_*sizeof(int));
                        for (int l = 0; l < amo_; l++) {
                            ibas_aaaa_sym[h][i][j][k][l] = -999;
                            ibas_aaab_sym[h][i][j][k][l] = -999;
                            ibas_aabb_sym[h][i][j][k][l] = -999;
                        }
                    }
                }
            }
            for (int i = 0; i < amo_*amo_*amo_*amo_; i++) {
                bas_aaaa_sym[h][i] = (int*)malloc(4*sizeof(int));
                bas_aaab_sym[h][i] = (int*)malloc(4*sizeof(int));
                bas_aabb_sym[h][i] = (int*)malloc(4*sizeof(int));
                for (int j = 0; j < 4; j++) {
                    bas_aaaa_sym[h][i][j] = -999;
                    bas_aaab_sym[h][i][j] = -999;
                    bas_aabb_sym[h][i][j] = -999;
                    bas_aabb_sym[h][i][j] = -999;
                }
            }
            // mappings:
            int count_aaaa = 0;
            int count_aaab = 0;
            int count_aabb = 0;
            for (int n = 0; n < quartets[h].size(); n++) {
                int i = std::get<0>(quartets[h][n]);
                int j = std::get<1>(quartets[h][n]);
                int k = std::get<2>(quartets[h][n]);
                int l = std::get<3>(quartets[h][n]);

                if ( i < j && k < l ) {

                    ibas_aabb_sym[h][i][j][k][l] = count_aabb;
                    ibas_aabb_sym[h][j][i][k][l] = count_aabb;
                    ibas_aabb_sym[h][i][j][l][k] = count_aabb;
                    ibas_aabb_sym[h][j][i][l][k] = count_aabb;

                    bas_aabb_sym[h][count_aabb][0]  = i;
                    bas_aabb_sym[h][count_aabb][1]  = j;
                    bas_aabb_sym[h][count_aabb][2]  = k;
                    bas_aabb_sym[h][count_aabb][3]  = l;
                    count_aabb++;
                }
                if ( i < j && j < k ) {

                    ibas_aaab_sym[h][i][j][k][l] = count_aaab;
                    ibas_aaab_sym[h][i][k][j][l] = count_aaab;
                    ibas_aaab_sym[h][j][i][k][l] = count_aaab;
                    ibas_aaab_sym[h][j][k][i][l] = count_aaab;
                    ibas_aaab_sym[h][k][i][j][l] = count_aaab;
                    ibas_aaab_sym[h][k][j][i][l] = count_aaab;

                    bas_aaab_sym[h][count_aaab][0]  = i;
                    bas_aaab_sym[h][count_aaab][1]  = j;
                    bas_aaab_sym[h][count_aaab][2]  = k;
                    bas_aaab_sym[h][count_aaab][3]  = l;
                    count_aaab++;
                }
                if ( i < j && j < k  && k < l) {

                    ibas_aaaa_sym[h][i][j][k][l] = count_aaaa;
                    ibas_aaaa_sym[h][i][k][j][l] = count_aaaa;
                    ibas_aaaa_sym[h][j][i][k][l] = count_aaaa;
                    ibas_aaaa_sym[h][j][k][i][l] = count_aaaa;
                    ibas_aaaa_sym[h][k][i][j][l] = count_aaaa;
                    ibas_aaaa_sym[h][k][j][i][l] = count_aaaa;

                    ibas_aaaa_sym[h][i][j][l][k] = count_aaaa;
                    ibas_aaaa_sym[h][i][k][l][j] = count_aaaa;
                    ibas_aaaa_sym[h][j][i][l][k] = count_aaaa;
                    ibas_aaaa_sym[h][j][k][l][i] = count_aaaa;
                    ibas_aaaa_sym[h][k][i][l][j] = count_aaaa;
                    ibas_aaaa_sym[h][k][j][l][i] = count_aaaa;

                    ibas_aaaa_sym[h][i][l][j][k] = count_aaaa;
                    ibas_aaaa_sym[h][i][l][k][j] = count_aaaa;
                    ibas_aaaa_sym[h][j][l][i][k] = count_aaaa;
                    ibas_aaaa_sym[h][j][l][k][i] = count_aaaa;
                    ibas_aaaa_sym[h][k][l][i][j] = count_aaaa;
                    ibas_aaaa_sym[h][k][l][j][i] = count_aaaa;

                    ibas_aaaa_sym[h][l][i][j][k] = count_aaaa;
                    ibas_aaaa_sym[h][l][i][k][j] = count_aaaa;
                    ibas_aaaa_sym[h][l][j][i][k] = count_aaaa;
                    ibas_aaaa_sym[h][l][j][k][i] = count_aaaa;
                    ibas_aaaa_sym[h][l][k][i][j] = count_aaaa;
                    ibas_aaaa_sym[h][l][k][j][i] = count_aaaa;

                    bas_aaaa_sym[h][count_aaaa][0]  = i;
                    bas_aaaa_sym[h][count_aaaa][1]  = j;
                    bas_aaaa_sym[h][count_aaaa][2]  = k;
                    bas_aaaa_sym[h][count_aaaa][3]  = l;
                    count_aaaa++;
                }
            }
            quartet_aaaa[h] = count_aaaa;
            quartet_aaab[h] = count_aaab;
            quartet_aabb[h] = count_aabb;
        }
    }

    free(gems_really_full);
}


}
