/*
 *@BEGIN LICENSE
 *
 * v2RDM-DOCI, a plugin to:
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

#include<psi4/libtrans/integraltransform.h>
#include<psi4/libtrans/mospace.h>

#include<psi4/libmints/wavefunction.h>
//#include<psi4/libmints/mints.h>
#include<psi4/libmints/vector.h>
#include<psi4/libmints/matrix.h>
//#include<../bin/fnocc/blas.h>
#include<time.h>

#include"v2rdm_solver.h"

#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() ( (double)clock() / CLOCKS_PER_SEC )
    #define omp_get_max_threads() 1
#endif

using namespace psi;
//using namespace fnocc;

namespace psi{ namespace v2rdm_doci{


void v2RDMSolver::BuildBasis() {

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

            if ( i <= j ) continue;

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

    free(gems_really_full);
}


}}
