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
#include<psi4/libmints/vector.h>
#include<psi4/libmints/matrix.h>

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

void v2RDMSolver::GetIntegrals() {

    // one-electron integrals:  
    SharedMatrix K1 = GetOEI();

    // size of the tei buffer
    if ( is_df_ ) {

        // size of the 3-index integral buffer
        tei_full_dim_ = (long int) nQ_ * (long int) ( nmo_ - nfrzv_ ) * ( (long int) ( nmo_ - nfrzv_ ) + 1L ) / 2L ;

        // just point to 3-index integral buffer
        tei_full_sym_      = Qmo_;

    }else {

        // size of the 4-index integral buffer
        tei_full_dim_ = 0;
        for (int h = 0; h < nirrep_; h++) {
            tei_full_dim_ += (long int)gems_full[h] * ( (long int)gems_full[h] + 1L ) / 2L;
        }

        tei_full_sym_ = (double*)malloc(tei_full_dim_*sizeof(double));
        memset((void*)tei_full_sym_,'\0',tei_full_dim_*sizeof(double));

    }

    // size of d2, blocked by symmetry, including the core orbitals
    d2_plus_core_dim_ = 0;
    for (int h = 0; h < nirrep_; h++) {
        d2_plus_core_dim_ += (long int)gems_plus_core[h] * ( (long int)gems_plus_core[h] + 1L ) / 2L;
        //ionix apriori no tenemos nada en core, por tanto queda:
    }
    d2_plus_core_sym_  = (double*)malloc(d2_plus_core_dim_*sizeof(double));
    memset((void*)d2_plus_core_sym_,'\0',d2_plus_core_dim_*sizeof(double));

    d2_act_spatial_dim_ = 0;
    for (int h = 0; h < nirrep_; h++) {
        d2_act_spatial_dim_ += (long int)gems_00[h] * ( (long int)gems_00[h] + 1L ) / 2L;
    }
    d2_act_spatial_sym_  = (double*)malloc(d2_act_spatial_dim_*sizeof(double));
    memset((void*)d2_act_spatial_sym_,'\0',d2_act_spatial_dim_*sizeof(double));

    // allocate memory for oei tensor, blocked by symmetry, excluding frozen virtuals
    oei_full_dim_ = 0;
    for (int h = 0; h < nirrep_; h++) {
        oei_full_dim_ += nmopi_[h] * ( nmopi_[h] + 1 ) / 2;
    }


    // allocate memory for d1 tensor, blocked by symmetry, including the core orbitals
    //gg -- only active orbitals are stored now (old code below)
    d1_act_spatial_dim_ = 0;
    for ( int h = 0; h < nirrep_; h++) {
        d1_act_spatial_dim_ += amopi_[h] * ( amopi_[h] + 1 ) / 2;
    }

    oei_full_sym_ = (double*)malloc(oei_full_dim_*sizeof(double));
    memset((void*)oei_full_sym_,'\0',oei_full_dim_*sizeof(double));

    d1_act_spatial_sym_ = (double*)malloc(d1_act_spatial_dim_*sizeof(double));
    memset((void*)d1_act_spatial_sym_,'\0',d1_act_spatial_dim_*sizeof(double));

    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < nmopi_[h] - frzvpi_[h]; i++) {
            for (long int j = 0; j < nmopi_[h] - frzvpi_[h]; j++) {
                oei_full_sym_[offset + INDEX(i,j)] = K1->pointer(h)[i][j];
            }
        }
        offset += amopi_[h] * ( amopi_[h] + 1 ) / 2;
    }

    if ( !is_df_ ) {
        // read tei's from disk
        GetTEIFromDisk();
    }

    RepackIntegrals();

}

// repack rotated full-space integral into active-space integrals
void v2RDMSolver::RepackIntegrals(){ 

    FrozenCoreEnergy();

    double * c_p = c->pointer();
    // two-electron part
    long int na = nalpha_;

    for (int i = 0; i < amo_; i++) {

        int hi = symmetry[i];

        for (int j = 0; j < amo_; j++) {

            int hj  = symmetry[j];
            int hij = SymmetryPair(hi,hj);

            if ( i != j ){
                c_p[d2s2off_ + i*amo_+j]  = 2.0 * TEI(i,i,j,j,0);
                c_p[d2s2off_ + i*amo_+j] -=       TEI(i,j,i,j,hij);
                c_p[d2s2off_ + i*amo_+j] += 2.0/(nalpha_+nbeta_-1.0) * (oei_full_sym_[INDEX(i,i)] + oei_full_sym_[INDEX(j,j)]);
            }
            c_p[d2s0off_ + i*amo_+j]  = TEI(i,j,i,j,hij);
            c_p[d2s0off_ + i*amo_+j] += 2.0/(nalpha_+nbeta_-1.0) * (i==j) * oei_full_sym_[INDEX(i,i)];

        }
    }
}
void v2RDMSolver::FrozenCoreEnergy() {

    // if frozen core, adjust oei's and compute frozen core energy:
    efzc_ = 0.0;
    double * c_p = c->pointer();


    //offset = 0;
    //for (int h = 0; h < nirrep_; h++) {
    //    for (int i = 0; i < amopi_[h]; i++) {
    //        int ii = i + pitzer_offset[h];
    //        c_p[d1off_ + ii] = 2.0*oei_full_sym_[offset + INDEX(i,i)];
    //    }
    //    offset += nmopi_[h]*(nmopi_[h]+1)/2;
    //}

}

double v2RDMSolver::TEI(int i, int j, int k, int l, int h) {
    double dum = 0.0;

    if ( is_df_ ) {

        //dum = C_DDOT(nQ_,Qmo_ + nQ_*INDEX(i,j),1,Qmo_+nQ_*INDEX(k,l),1);
        dum = C_DDOT(nQ_,Qmo_ + INDEX(i,j),(nmo_-nfrzv_)*(nmo_-nfrzv_+1)/2,Qmo_+INDEX(k,l),(nmo_-nfrzv_)*(nmo_-nfrzv_+1)/2);

    }else {

        int myoff = 0;
        for (int myh = 0; myh < h; myh++) {
            myoff += (long int)gems_full[myh] * ( (long int)gems_full[myh] + 1L ) / 2L;
        }

        int ij    = ibas_full_sym[h][i][j];
        int kl    = ibas_full_sym[h][k][l];

        dum = tei_full_sym_[myoff + INDEX(ij,kl)];

    }

    return dum;
}



}}
