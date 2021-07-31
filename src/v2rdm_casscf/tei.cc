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

#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libqt/qt.h>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libtrans/mospace.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/matrix.h>

#include "v2rdm_solver.h"

#include <misc/omp.h>

using namespace psi;

namespace hilbert{

void v2RDMSolver::GetIntegrals() {

    
    // one-electron integrals:  
    SharedMatrix K1;
    if ( is_hubbard_ ) {
        K1 = GetOEI_hubbard();
    }else {
        K1 = GetOEI();
    }

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
        oei_full_dim_ += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    // allocate memory for d1 tensor, blocked by symmetry, including the core orbitals
    //gg -- only active orbitals are stored now (old code below)
    d1_act_spatial_dim_ = 0;
    for ( int h = 0; h < nirrep_; h++) {
        d1_act_spatial_dim_ += amopi_[h] * ( amopi_[h] + 1 ) / 2;
    }

/*
    d1_act_spatial_dim_ = 0;
    for ( int h = 0; h < nirrep_; h++) {
        d1_act_spatial_dim_ += (rstcpi_[h] + frzcpi_[h] + amopi_[h]) * ( rstcpi_[h] + frzcpi_[h] + amopi_[h] + 1 ) / 2;
    }
*/

    oei_full_sym_ = (double*)malloc(oei_full_dim_*sizeof(double));
    memset((void*)oei_full_sym_,'\0',oei_full_dim_*sizeof(double));

    d1_act_spatial_sym_ = (double*)malloc(d1_act_spatial_dim_*sizeof(double));
    memset((void*)d1_act_spatial_sym_,'\0',d1_act_spatial_dim_*sizeof(double));

    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < nmopi_[h] - frzvpi_[h]; i++) {
            for (long int j = i; j < nmopi_[h] - frzvpi_[h]; j++) {
                oei_full_sym_[offset + INDEX(i,j)] = K1->pointer(h)[i][j];
            }
        }
        offset += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    if ( !is_df_ && !is_hubbard_ ) {
        // read tei's from disk
        GetTEIFromDisk();
    }
    RepackIntegrals();

}

// repack rotated full-space integrals into active-space integrals
void v2RDMSolver::RepackIntegrals(){

    FrozenCoreEnergy();

    double * c_p = c->pointer();

    // two-electron part
    long int na = nalpha_ - nrstc_ - nfrzc_;
    long int nb = nbeta_ - nrstc_ - nfrzc_;
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (long int ij = 0; ij < gems_ab[h]; ij++) {
            long int i = bas_ab_sym[h][ij][0];
            long int j = bas_ab_sym[h][ij][1];

            long int ii = full_basis[i];
            long int jj = full_basis[j];

            for (long int kl = 0; kl < gems_ab[h]; kl++) {
                long int k = bas_ab_sym[h][kl][0];
                long int l = bas_ab_sym[h][kl][1];

                long int kk = full_basis[k];
                long int ll = full_basis[l];

                int hik = SymmetryPair(symmetry[i],symmetry[k]);

                c_p[d2aboff[h] + ij*gems_ab[h]+kl] = TEI(ii,kk,jj,ll,hik);

            }
        }
    }

    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (long int ij = 0; ij < gems_aa[h]; ij++) {
            long int i = bas_aa_sym[h][ij][0];
            long int j = bas_aa_sym[h][ij][1];

            long int ii = full_basis[i];
            long int jj = full_basis[j];

            for (long int kl = 0; kl < gems_aa[h]; kl++) {
                long int k = bas_aa_sym[h][kl][0];
                long int l = bas_aa_sym[h][kl][1];

                long int kk = full_basis[k];
                long int ll = full_basis[l];

                int hik = SymmetryPair(symmetry[i],symmetry[k]);
                int hil = SymmetryPair(symmetry[i],symmetry[l]);

                double dum1 = TEI(ii,kk,jj,ll,hik);
                double dum2 = TEI(ii,ll,jj,kk,hil);

                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]    = dum1 - dum2;
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]    = dum1 - dum2;
            }
        }
    }
}

void v2RDMSolver::FrozenCoreEnergy() {

    // if frozen core, adjust oei's and compute frozen core energy:
    efzc_ = 0.0;
    offset = 0;
    long int offset3 = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < rstcpi_[h] + frzcpi_[h]; i++) {

            int ifull = i + offset;

            efzc_ += 2.0 * oei_full_sym_[offset3 + INDEX(i,i)];

            long int offset2 = 0;
            for (int h2 = 0; h2 < nirrep_; h2++) {
                for (int j = 0; j < rstcpi_[h2] + frzcpi_[h2]; j++) {

                    int jfull = j + offset2;

                    int hij = SymmetryPair(h,h2);

                    double dum1 = TEI(ifull,ifull,jfull,jfull, 0);
                    double dum2 = TEI(ifull,jfull,ifull,jfull, hij);
                    efzc_ += 2.0 * dum1 - dum2;

                }
                offset2 += nmopi_[h2] - frzvpi_[h2];
            }
        }
        offset += nmopi_[h] - frzvpi_[h];
        offset3 += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    double * c_p = c->pointer();

    // adjust one-electron integrals for core repulsion contribution
    offset = 0;
    offset3 = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h] - rstvpi_[h] - frzvpi_[h]; i++) {

            int ifull = i + offset;

            for (int j = rstcpi_[h] + frzcpi_[h]; j < nmopi_[h] - rstvpi_[h] - frzvpi_[h]; j++) {

                int jfull = j + offset;

                double dum = 0.0;

                long int offset2 = 0;
                for (int h2 = 0; h2 < nirrep_; h2++) {
                    for (int k = 0; k < rstcpi_[h2] + frzcpi_[h2]; k++) {

                        int kfull = k + offset2;

                        int hik = SymmetryPair(h,h2);

                        double dum1 = TEI(ifull,jfull,kfull,kfull,0);
                        double dum2 = TEI(ifull,kfull,jfull,kfull,hik);
                        dum += 2.0 * dum1 - dum2;

                    }
                    offset2 += nmopi_[h2] - frzvpi_[h2];
                }
                c_p[d1aoff[h] + (i-rstcpi_[h]-frzcpi_[h])*amopi_[h] + (j-rstcpi_[h]-frzcpi_[h])] = oei_full_sym_[offset3+INDEX(i,j)];
                c_p[d1boff[h] + (i-rstcpi_[h]-frzcpi_[h])*amopi_[h] + (j-rstcpi_[h]-frzcpi_[h])] = oei_full_sym_[offset3+INDEX(i,j)];
                c_p[d1aoff[h] + (i-rstcpi_[h]-frzcpi_[h])*amopi_[h] + (j-rstcpi_[h]-frzcpi_[h])] += dum;
                c_p[d1boff[h] + (i-rstcpi_[h]-frzcpi_[h])*amopi_[h] + (j-rstcpi_[h]-frzcpi_[h])] += dum;
            }
        }
        offset += nmopi_[h] - frzvpi_[h];
        offset3 += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

}

double v2RDMSolver::TEI(int i, int j, int k, int l, int h) {
    double dum = 0.0;

    if ( is_df_ ) {

        //dum = C_DDOT(nQ_,Qmo_ + nQ_*INDEX(i,j),1,Qmo_+nQ_*INDEX(k,l),1);
        dum = C_DDOT(nQ_,Qmo_ + INDEX(i,j),(nmo_-nfrzv_)*(nmo_-nfrzv_+1)/2,Qmo_+INDEX(k,l),(nmo_-nfrzv_)*(nmo_-nfrzv_+1)/2);

    }else if ( is_hubbard_ ) {

        if ( i == j && k == l  && i == k ) {
            dum = options_.get_double("HUBBARD_U");
        }else {
            dum = 0.0;
        }

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

}
