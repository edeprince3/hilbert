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

#include "v2rdm_solver.h"

#include <misc/omp.h>

using namespace psi;

namespace hilbert{

// T1 portion of A.u 
void v2RDMSolver::T1_constraints_guess(double* u){

    // T1aab
    for (int h = 0; h < nirrep_; h++) {

        #pragma omp parallel for schedule (static)
        for (int ijk = 0; ijk < trip_aab[h]; ijk++) {

            int i = bas_aab_sym[h][ijk][0];
            int j = bas_aab_sym[h][ijk][1];
            int k = bas_aab_sym[h][ijk][2];

            for (int lmn = 0; lmn < trip_aab[h]; lmn++) {

                int l = bas_aab_sym[h][lmn][0];
                int m = bas_aab_sym[h][lmn][1];
                int n = bas_aab_sym[h][lmn][2];

                //double dum = -u[t1aaboff[h] + ijk*trip_aab[h]+lmn]; // - T1(ijk,lmn)
                double dum = 0.0;//-u[t1aaboff[h] + ijk*trip_aab[h]+lmn]; // - T1(ijk,lmn)

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    dum += u[q2aaoff[hij] + ij*gems_aa[hij] + lm];  // Q2(ij,lm) dkn
                }

                if ( j == l ) {
                    int hki = SymmetryPair(symmetry[k],symmetry[i]);
                    int nm = ibas_ab_sym[hki][m][n];
                    int ki = ibas_ab_sym[hki][i][k];
                    dum -= u[d2aboff[hki] + nm*gems_ab[hki] + ki];  // -D2(nm,ki) dlj
                }

                if ( l == i ) {
                    int hkj = SymmetryPair(symmetry[k],symmetry[j]);
                    int nm = ibas_ab_sym[hkj][m][n];
                    int kj = ibas_ab_sym[hkj][j][k];
                    dum += u[d2aboff[hkj] + nm*gems_ab[hkj] + kj];  // D2(nm,kj) dli
                }

                if ( j == m ) {
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int kl = ibas_ab_sym[hni][k][l];
                    dum -= u[g2baoff[hni] + ni*gems_ab[hni] + kl];  // -G2(ni,kl) djm
                    
                }

                if ( i == m ) {
                    int hkl = SymmetryPair(symmetry[k],symmetry[l]);
                    int nj = ibas_ab_sym[hkl][n][j];
                    int kl = ibas_ab_sym[hkl][k][l];
                    dum += u[g2baoff[hkl] + nj*gems_ab[hkl] + kl];  // G2(nj,kl) dim
                    
                }

                u[t1aaboff[h] + ijk*trip_aab[h]+lmn] = dum;

            }
        }
        offset += trip_aab[h]*trip_aab[h];

    }
    // T1bba
    for (int h = 0; h < nirrep_; h++) {

        #pragma omp parallel for schedule (static)
        for (int ijk = 0; ijk < trip_aab[h]; ijk++) {

            int i = bas_aab_sym[h][ijk][0];
            int j = bas_aab_sym[h][ijk][1];
            int k = bas_aab_sym[h][ijk][2];

            for (int lmn = 0; lmn < trip_aab[h]; lmn++) {

                int l = bas_aab_sym[h][lmn][0];
                int m = bas_aab_sym[h][lmn][1];
                int n = bas_aab_sym[h][lmn][2];

                double dum = 0.0;//-u[t1bbaoff[h] + ijk*trip_aab[h]+lmn]; // - T1(ijk,lmn)

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    dum += u[q2bboff[hij] + ij*gems_aa[hij] + lm];  // Q2(ij,lm) dkn
                }

                if ( j == l ) {
                    int hki = SymmetryPair(symmetry[k],symmetry[i]);
                    int nm = ibas_ab_sym[hki][n][m];
                    int ki = ibas_ab_sym[hki][k][i];
                    dum -= u[d2aboff[hki] + nm*gems_ab[hki] + ki];  // -D2(nm,ki) dlj
                }

                if ( l == i ) {
                    int hkj = SymmetryPair(symmetry[k],symmetry[j]);
                    int nm = ibas_ab_sym[hkj][n][m];
                    int kj = ibas_ab_sym[hkj][k][j];
                    dum += u[d2aboff[hkj] + nm*gems_ab[hkj] + kj];  // D2(nm,kj) dli
                }

                if ( j == m ) {
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int kl = ibas_ab_sym[hni][k][l];
                    dum -= u[g2aboff[hni] + ni*gems_ab[hni] + kl];  // -G2(ni,kl) djm
                    
                }

                if ( i == m ) {
                    int hkl = SymmetryPair(symmetry[k],symmetry[l]);
                    int nj = ibas_ab_sym[hkl][n][j];
                    int kl = ibas_ab_sym[hkl][k][l];
                    dum += u[g2aboff[hkl] + nj*gems_ab[hkl] + kl];  // G2(nj,kl) dim
                    
                }


                u[t1bbaoff[h] + ijk*trip_aab[h]+lmn] = dum; // - T1(ijk,lmn)

            }
        }
        offset += trip_aab[h]*trip_aab[h];

    }
    // T1aaa
    for (int h = 0; h < nirrep_; h++) {

        #pragma omp parallel for schedule (static)
        for (int ijk = 0; ijk < trip_aaa[h]; ijk++) {

            int i = bas_aaa_sym[h][ijk][0];
            int j = bas_aaa_sym[h][ijk][1];
            int k = bas_aaa_sym[h][ijk][2];

            for (int lmn = 0; lmn < trip_aaa[h]; lmn++) {

                int l = bas_aaa_sym[h][lmn][0];
                int m = bas_aaa_sym[h][lmn][1];
                int n = bas_aaa_sym[h][lmn][2];

                double dum = 0.0;//-u[t1aaaoff[h] + ijk*trip_aaa[h]+lmn]; // - T1(ijk,lmn)

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    dum += u[q2aaoff[hij] + ij*gems_aa[hij] + lm];  // Q2(ij,lm) dkn
                }

                if ( j == n ) {
                    int hik = SymmetryPair(symmetry[i],symmetry[k]);
                    int hlm = SymmetryPair(symmetry[l],symmetry[m]);
                    if ( hik == hlm ) {
                        int ik = ibas_aa_sym[hik][i][k];
                        int lm = ibas_aa_sym[hik][l][m];
                        dum -= u[q2aaoff[hik] + ik*gems_aa[hik] + lm];  // -Q2(ik,lm) djn
                    }
                }

                if ( i == n ) {
                    int hjk = SymmetryPair(symmetry[j],symmetry[k]);
                    int hlm = SymmetryPair(symmetry[l],symmetry[m]);
                    if ( hjk == hlm ) {
                        int jk = ibas_aa_sym[hjk][j][k];
                        int lm = ibas_aa_sym[hjk][l][m];
                        dum += u[q2aaoff[hjk] + jk*gems_aa[hjk] + lm];  // Q2(jk,lm) din
                    }
                }


                if ( l == k ) {
                    int hnm = SymmetryPair(symmetry[n],symmetry[m]);
                    int hji = SymmetryPair(symmetry[j],symmetry[i]);
                    if ( hji == hnm ) {
                        int nm = ibas_aa_sym[hji][n][m];
                        int ji = ibas_aa_sym[hji][j][i];
                        dum += u[d2aaoff[hji] + nm*gems_aa[hji] + ji];  // D2(nm,ji) dlk
                    }
                }

                if ( j == l ) {
                    int hki = SymmetryPair(symmetry[k],symmetry[i]);
                    int nm = ibas_aa_sym[hki][n][m];
                    int ki = ibas_aa_sym[hki][k][i];
                    dum -= u[d2aaoff[hki] + nm*gems_aa[hki] + ki];  // -D2(nm,ki) dlj
                }

                if ( l == i ) {
                    int hkj = SymmetryPair(symmetry[k],symmetry[j]);
                    int nm = ibas_aa_sym[hkj][n][m];
                    int kj = ibas_aa_sym[hkj][k][j];
                    dum += u[d2aaoff[hkj] + nm*gems_aa[hkj] + kj];  // D2(nm,kj) dli
                }

                if ( k == m ) {
                    if ( j == l ) {
                        int h2 = symmetry[n];
                        int nn = n - pitzer_offset[h2];
                        int ii = i - pitzer_offset[h2];
                        dum -= u[d1aoff[h2] + nn*amopi_[h2]+ii]; // - D1(n,i) djl dkm
                    }
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int jl = ibas_ab_sym[hni][j][l];
                    dum += u[g2aaoff[hni] + ni*2*gems_ab[hni] + jl];  // G2(ni,jl) dkm
                    
                }

                if ( j == m ) {
                    if ( k == l ) {
                        int h2 = symmetry[n];
                        int nn = n - pitzer_offset[h2];
                        int ii = i - pitzer_offset[h2];
                        dum += u[d1aoff[h2] + nn*amopi_[h2]+ii]; // D1(n,i) dkl djm
                    }
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int kl = ibas_ab_sym[hni][k][l];
                    dum -= u[g2aaoff[hni] + ni*2*gems_ab[hni] + kl];  // -G2(ni,kl) djm
                    
                }

                if ( i == m ) {
                    if ( k == l ) {
                        int h2 = symmetry[n];
                        int nn = n - pitzer_offset[h2];
                        int jj = j - pitzer_offset[h2];
                        dum -= u[d1aoff[h2] + nn*amopi_[h2]+jj]; // - D1(n,j) dkl dim
                    }
                    int hkl = SymmetryPair(symmetry[k],symmetry[l]);
                    int nj = ibas_ab_sym[hkl][n][j];
                    int kl = ibas_ab_sym[hkl][k][l];
                    dum += u[g2aaoff[hkl] + nj*2*gems_ab[hkl] + kl];  // G2(nj,kl) dim
                    
                }


                u[t1aaaoff[h] + ijk*trip_aaa[h]+lmn] = dum; // - T1(ijk,lmn)


            }
        }
        offset += trip_aaa[h]*trip_aaa[h];

    }
    // T1bbb
    for (int h = 0; h < nirrep_; h++) {

        #pragma omp parallel for schedule (static)
        for (int ijk = 0; ijk < trip_aaa[h]; ijk++) {

            int i = bas_aaa_sym[h][ijk][0];
            int j = bas_aaa_sym[h][ijk][1];
            int k = bas_aaa_sym[h][ijk][2];

            for (int lmn = 0; lmn < trip_aaa[h]; lmn++) {

                int l = bas_aaa_sym[h][lmn][0];
                int m = bas_aaa_sym[h][lmn][1];
                int n = bas_aaa_sym[h][lmn][2];

                double dum = 0.0;//-u[t1bbboff[h] + ijk*trip_aaa[h]+lmn]; // - T1(ijk,lmn)

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    dum += u[q2bboff[hij] + ij*gems_aa[hij] + lm];  // Q2(ij,lm) dkn
                }

                if ( j == n ) {
                    int hik = SymmetryPair(symmetry[i],symmetry[k]);
                    int hlm = SymmetryPair(symmetry[l],symmetry[m]);
                    if ( hik == hlm ) {
                        int ik = ibas_aa_sym[hik][i][k];
                        int lm = ibas_aa_sym[hik][l][m];
                        dum -= u[q2bboff[hik] + ik*gems_aa[hik] + lm];  // -Q2(ik,lm) djn
                    }
                }

                if ( i == n ) {
                    int hjk = SymmetryPair(symmetry[j],symmetry[k]);
                    int hlm = SymmetryPair(symmetry[l],symmetry[m]);
                    if ( hjk == hlm ) {
                        int jk = ibas_aa_sym[hjk][j][k];
                        int lm = ibas_aa_sym[hjk][l][m];
                        dum += u[q2bboff[hjk] + jk*gems_aa[hjk] + lm];  // Q2(jk,lm) din
                    }
                }


                if ( l == k ) {
                    int hnm = SymmetryPair(symmetry[n],symmetry[m]);
                    int hji = SymmetryPair(symmetry[j],symmetry[i]);
                    if ( hji == hnm ) {
                        int nm = ibas_aa_sym[hji][n][m];
                        int ji = ibas_aa_sym[hji][j][i];
                        dum += u[d2bboff[hji] + nm*gems_aa[hji] + ji];  // D2(nm,ji) dlk
                    }
                }

                if ( j == l ) {
                    int hki = SymmetryPair(symmetry[k],symmetry[i]);
                    int nm = ibas_aa_sym[hki][n][m];
                    int ki = ibas_aa_sym[hki][k][i];
                    dum -= u[d2bboff[hki] + nm*gems_aa[hki] + ki];  // -D2(nm,ki) dlj
                }

                if ( l == i ) {
                    int hkj = SymmetryPair(symmetry[k],symmetry[j]);
                    int nm = ibas_aa_sym[hkj][n][m];
                    int kj = ibas_aa_sym[hkj][k][j];
                    dum += u[d2bboff[hkj] + nm*gems_aa[hkj] + kj];  // D2(nm,kj) dli
                }

                if ( k == m ) {
                    if ( j == l ) {
                        int h2 = symmetry[n];
                        int nn = n - pitzer_offset[h2];
                        int ii = i - pitzer_offset[h2];
                        dum -= u[d1boff[h2] + nn*amopi_[h2]+ii]; // - D1(n,i) djl dkm
                    }
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int jl = ibas_ab_sym[hni][j][l];
                    dum += u[g2aaoff[hni] + (ni+gems_ab[hni])*2*gems_ab[hni] + (jl+gems_ab[hni])];  // G2(ni,jl) dkm
                    
                }

                if ( j == m ) {
                    if ( k == l ) {
                        int h2 = symmetry[n];
                        int nn = n - pitzer_offset[h2];
                        int ii = i - pitzer_offset[h2];
                        dum += u[d1boff[h2] + nn*amopi_[h2]+ii]; // D1(n,i) dkl djm
                    }
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int kl = ibas_ab_sym[hni][k][l];
                    dum -= u[g2aaoff[hni] + (ni+gems_ab[hni])*2*gems_ab[hni] + (kl+gems_ab[hni])];  // -G2(ni,kl) djm
                    
                }

                if ( i == m ) {
                    if ( k == l ) {
                        int h2 = symmetry[n];
                        int nn = n - pitzer_offset[h2];
                        int jj = j - pitzer_offset[h2];
                        dum -= u[d1boff[h2] + nn*amopi_[h2]+jj]; // - D1(n,j) dkl dim
                    }
                    int hkl = SymmetryPair(symmetry[k],symmetry[l]);
                    int nj = ibas_ab_sym[hkl][n][j];
                    int kl = ibas_ab_sym[hkl][k][l];
                    dum += u[g2aaoff[hkl] + (nj+gems_ab[hkl])*2*gems_ab[hkl] + (kl+gems_ab[hkl])];  // G2(nj,kl) dim
                    
                }


                u[t1bbboff[h] + ijk*trip_aaa[h]+lmn] = dum; // - T1(ijk,lmn)


            }
        }
        offset += trip_aaa[h]*trip_aaa[h];

    }

}

// T1 portion of A.u 
void v2RDMSolver::T1_constraints_Au(double* A,double* u){

    // T1aab
    for (int h = 0; h < nirrep_; h++) {

        #pragma omp parallel for schedule (static)
        for (int ijk = 0; ijk < trip_aab[h]; ijk++) {

            int i = bas_aab_sym[h][ijk][0];
            int j = bas_aab_sym[h][ijk][1];
            int k = bas_aab_sym[h][ijk][2];

            for (int lmn = 0; lmn < trip_aab[h]; lmn++) {

                int l = bas_aab_sym[h][lmn][0];
                int m = bas_aab_sym[h][lmn][1];
                int n = bas_aab_sym[h][lmn][2];

                double dum = -u[t1aaboff[h] + ijk*trip_aab[h]+lmn]; // - T1(ijk,lmn)

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    dum += u[q2aaoff[hij] + ij*gems_aa[hij] + lm];  // Q2(ij,lm) dkn
                }

                if ( j == l ) {
                    int hki = SymmetryPair(symmetry[k],symmetry[i]);
                    int nm = ibas_ab_sym[hki][m][n];
                    int ki = ibas_ab_sym[hki][i][k];
                    dum -= u[d2aboff[hki] + nm*gems_ab[hki] + ki];  // -D2(nm,ki) dlj
                }

                if ( l == i ) {
                    int hkj = SymmetryPair(symmetry[k],symmetry[j]);
                    int nm = ibas_ab_sym[hkj][m][n];
                    int kj = ibas_ab_sym[hkj][j][k];
                    dum += u[d2aboff[hkj] + nm*gems_ab[hkj] + kj];  // D2(nm,kj) dli
                }

                if ( j == m ) {
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int kl = ibas_ab_sym[hni][k][l];
                    dum -= u[g2baoff[hni] + ni*gems_ab[hni] + kl];  // -G2(ni,kl) djm
                    
                }

                if ( i == m ) {
                    int hkl = SymmetryPair(symmetry[k],symmetry[l]);
                    int nj = ibas_ab_sym[hkl][n][j];
                    int kl = ibas_ab_sym[hkl][k][l];
                    dum += u[g2baoff[hkl] + nj*gems_ab[hkl] + kl];  // G2(nj,kl) dim
                    
                }


                A[offset + ijk*trip_aab[h]+lmn] = dum;


            }
        }
        offset += trip_aab[h]*trip_aab[h];

    }
    // T1bba
    for (int h = 0; h < nirrep_; h++) {

        #pragma omp parallel for schedule (static)
        for (int ijk = 0; ijk < trip_aab[h]; ijk++) {

            int i = bas_aab_sym[h][ijk][0];
            int j = bas_aab_sym[h][ijk][1];
            int k = bas_aab_sym[h][ijk][2];

            for (int lmn = 0; lmn < trip_aab[h]; lmn++) {

                int l = bas_aab_sym[h][lmn][0];
                int m = bas_aab_sym[h][lmn][1];
                int n = bas_aab_sym[h][lmn][2];

                double dum = -u[t1bbaoff[h] + ijk*trip_aab[h]+lmn]; // - T1(ijk,lmn)

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    dum += u[q2bboff[hij] + ij*gems_aa[hij] + lm];  // Q2(ij,lm) dkn
                }

                if ( j == l ) {
                    int hki = SymmetryPair(symmetry[k],symmetry[i]);
                    int nm = ibas_ab_sym[hki][n][m];
                    int ki = ibas_ab_sym[hki][k][i];
                    dum -= u[d2aboff[hki] + nm*gems_ab[hki] + ki];  // -D2(nm,ki) dlj
                }

                if ( l == i ) {
                    int hkj = SymmetryPair(symmetry[k],symmetry[j]);
                    int nm = ibas_ab_sym[hkj][n][m];
                    int kj = ibas_ab_sym[hkj][k][j];
                    dum += u[d2aboff[hkj] + nm*gems_ab[hkj] + kj];  // D2(nm,kj) dli
                }

                if ( j == m ) {
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int kl = ibas_ab_sym[hni][k][l];
                    dum -= u[g2aboff[hni] + ni*gems_ab[hni] + kl];  // -G2(ni,kl) djm
                    
                }

                if ( i == m ) {
                    int hkl = SymmetryPair(symmetry[k],symmetry[l]);
                    int nj = ibas_ab_sym[hkl][n][j];
                    int kl = ibas_ab_sym[hkl][k][l];
                    dum += u[g2aboff[hkl] + nj*gems_ab[hkl] + kl];  // G2(nj,kl) dim
                    
                }


                A[offset + ijk*trip_aab[h]+lmn] = dum;


            }
        }
        offset += trip_aab[h]*trip_aab[h];

    }
    // T1aaa
    for (int h = 0; h < nirrep_; h++) {

        #pragma omp parallel for schedule (static)
        for (int ijk = 0; ijk < trip_aaa[h]; ijk++) {

            int i = bas_aaa_sym[h][ijk][0];
            int j = bas_aaa_sym[h][ijk][1];
            int k = bas_aaa_sym[h][ijk][2];

            for (int lmn = 0; lmn < trip_aaa[h]; lmn++) {

                int l = bas_aaa_sym[h][lmn][0];
                int m = bas_aaa_sym[h][lmn][1];
                int n = bas_aaa_sym[h][lmn][2];

                double dum = -u[t1aaaoff[h] + ijk*trip_aaa[h]+lmn]; // - T1(ijk,lmn)

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    dum += u[q2aaoff[hij] + ij*gems_aa[hij] + lm];  // Q2(ij,lm) dkn
                }

                if ( j == n ) {
                    int hik = SymmetryPair(symmetry[i],symmetry[k]);
                    int hlm = SymmetryPair(symmetry[l],symmetry[m]);
                    if ( hik == hlm ) {
                        int ik = ibas_aa_sym[hik][i][k];
                        int lm = ibas_aa_sym[hik][l][m];
                        dum -= u[q2aaoff[hik] + ik*gems_aa[hik] + lm];  // -Q2(ik,lm) djn
                    }
                }

                if ( i == n ) {
                    int hjk = SymmetryPair(symmetry[j],symmetry[k]);
                    int hlm = SymmetryPair(symmetry[l],symmetry[m]);
                    if ( hjk == hlm ) {
                        int jk = ibas_aa_sym[hjk][j][k];
                        int lm = ibas_aa_sym[hjk][l][m];
                        dum += u[q2aaoff[hjk] + jk*gems_aa[hjk] + lm];  // Q2(jk,lm) din
                    }
                }


                if ( l == k ) {
                    int hnm = SymmetryPair(symmetry[n],symmetry[m]);
                    int hji = SymmetryPair(symmetry[j],symmetry[i]);
                    if ( hji == hnm ) {
                        int nm = ibas_aa_sym[hji][n][m];
                        int ji = ibas_aa_sym[hji][j][i];
                        dum += u[d2aaoff[hji] + nm*gems_aa[hji] + ji];  // D2(nm,ji) dlk
                    }
                }

                if ( j == l ) {
                    int hki = SymmetryPair(symmetry[k],symmetry[i]);
                    int nm = ibas_aa_sym[hki][n][m];
                    int ki = ibas_aa_sym[hki][k][i];
                    dum -= u[d2aaoff[hki] + nm*gems_aa[hki] + ki];  // -D2(nm,ki) dlj
                }

                if ( l == i ) {
                    int hkj = SymmetryPair(symmetry[k],symmetry[j]);
                    int nm = ibas_aa_sym[hkj][n][m];
                    int kj = ibas_aa_sym[hkj][k][j];
                    dum += u[d2aaoff[hkj] + nm*gems_aa[hkj] + kj];  // D2(nm,kj) dli
                }

                if ( k == m ) {
                    if ( j == l ) {
                        int h2 = symmetry[n];
                        int nn = n - pitzer_offset[h2];
                        int ii = i - pitzer_offset[h2];
                        dum -= u[d1aoff[h2] + nn*amopi_[h2]+ii]; // - D1(n,i) djl dkm
                    }
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int jl = ibas_ab_sym[hni][j][l];
                    dum += u[g2aaoff[hni] + ni*2*gems_ab[hni] + jl];  // G2(ni,jl) dkm
                    
                }

                if ( j == m ) {
                    if ( k == l ) {
                        int h2 = symmetry[n];
                        int nn = n - pitzer_offset[h2];
                        int ii = i - pitzer_offset[h2];
                        dum += u[d1aoff[h2] + nn*amopi_[h2]+ii]; // D1(n,i) dkl djm
                    }
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int kl = ibas_ab_sym[hni][k][l];
                    dum -= u[g2aaoff[hni] + ni*2*gems_ab[hni] + kl];  // -G2(ni,kl) djm
                    
                }

                if ( i == m ) {
                    if ( k == l ) {
                        int h2 = symmetry[n];
                        int nn = n - pitzer_offset[h2];
                        int jj = j - pitzer_offset[h2];
                        dum -= u[d1aoff[h2] + nn*amopi_[h2]+jj]; // - D1(n,j) dkl dim
                    }
                    int hkl = SymmetryPair(symmetry[k],symmetry[l]);
                    int nj = ibas_ab_sym[hkl][n][j];
                    int kl = ibas_ab_sym[hkl][k][l];
                    dum += u[g2aaoff[hkl] + nj*2*gems_ab[hkl] + kl];  // G2(nj,kl) dim
                    
                }


                A[offset + ijk*trip_aaa[h]+lmn] = dum;


            }
        }
        offset += trip_aaa[h]*trip_aaa[h];

    }
    // T1bbb
    for (int h = 0; h < nirrep_; h++) {

        #pragma omp parallel for schedule (static)
        for (int ijk = 0; ijk < trip_aaa[h]; ijk++) {

            int i = bas_aaa_sym[h][ijk][0];
            int j = bas_aaa_sym[h][ijk][1];
            int k = bas_aaa_sym[h][ijk][2];

            for (int lmn = 0; lmn < trip_aaa[h]; lmn++) {

                int l = bas_aaa_sym[h][lmn][0];
                int m = bas_aaa_sym[h][lmn][1];
                int n = bas_aaa_sym[h][lmn][2];

                double dum = -u[t1bbboff[h] + ijk*trip_aaa[h]+lmn]; // - T1(ijk,lmn)

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    dum += u[q2bboff[hij] + ij*gems_aa[hij] + lm];  // Q2(ij,lm) dkn
                }

                if ( j == n ) {
                    int hik = SymmetryPair(symmetry[i],symmetry[k]);
                    int hlm = SymmetryPair(symmetry[l],symmetry[m]);
                    if ( hik == hlm ) {
                        int ik = ibas_aa_sym[hik][i][k];
                        int lm = ibas_aa_sym[hik][l][m];
                        dum -= u[q2bboff[hik] + ik*gems_aa[hik] + lm];  // -Q2(ik,lm) djn
                    }
                }

                if ( i == n ) {
                    int hjk = SymmetryPair(symmetry[j],symmetry[k]);
                    int hlm = SymmetryPair(symmetry[l],symmetry[m]);
                    if ( hjk == hlm ) {
                        int jk = ibas_aa_sym[hjk][j][k];
                        int lm = ibas_aa_sym[hjk][l][m];
                        dum += u[q2bboff[hjk] + jk*gems_aa[hjk] + lm];  // Q2(jk,lm) din
                    }
                }


                if ( l == k ) {
                    int hnm = SymmetryPair(symmetry[n],symmetry[m]);
                    int hji = SymmetryPair(symmetry[j],symmetry[i]);
                    if ( hji == hnm ) {
                        int nm = ibas_aa_sym[hji][n][m];
                        int ji = ibas_aa_sym[hji][j][i];
                        dum += u[d2bboff[hji] + nm*gems_aa[hji] + ji];  // D2(nm,ji) dlk
                    }
                }

                if ( j == l ) {
                    int hki = SymmetryPair(symmetry[k],symmetry[i]);
                    int nm = ibas_aa_sym[hki][n][m];
                    int ki = ibas_aa_sym[hki][k][i];
                    dum -= u[d2bboff[hki] + nm*gems_aa[hki] + ki];  // -D2(nm,ki) dlj
                }

                if ( l == i ) {
                    int hkj = SymmetryPair(symmetry[k],symmetry[j]);
                    int nm = ibas_aa_sym[hkj][n][m];
                    int kj = ibas_aa_sym[hkj][k][j];
                    dum += u[d2bboff[hkj] + nm*gems_aa[hkj] + kj];  // D2(nm,kj) dli
                }

                if ( k == m ) {
                    if ( j == l ) {
                        int h2 = symmetry[n];
                        int nn = n - pitzer_offset[h2];
                        int ii = i - pitzer_offset[h2];
                        dum -= u[d1boff[h2] + nn*amopi_[h2]+ii]; // - D1(n,i) djl dkm
                    }
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int jl = ibas_ab_sym[hni][j][l];
                    dum += u[g2aaoff[hni] + (ni+gems_ab[hni])*2*gems_ab[hni] + (jl+gems_ab[hni])];  // G2(ni,jl) dkm
                    
                }

                if ( j == m ) {
                    if ( k == l ) {
                        int h2 = symmetry[n];
                        int nn = n - pitzer_offset[h2];
                        int ii = i - pitzer_offset[h2];
                        dum += u[d1boff[h2] + nn*amopi_[h2]+ii]; // D1(n,i) dkl djm
                    }
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int kl = ibas_ab_sym[hni][k][l];
                    dum -= u[g2aaoff[hni] + (ni+gems_ab[hni])*2*gems_ab[hni] + (kl+gems_ab[hni])];  // -G2(ni,kl) djm
                    
                }

                if ( i == m ) {
                    if ( k == l ) {
                        int h2 = symmetry[n];
                        int nn = n - pitzer_offset[h2];
                        int jj = j - pitzer_offset[h2];
                        dum -= u[d1boff[h2] + nn*amopi_[h2]+jj]; // - D1(n,j) dkl dim
                    }
                    int hkl = SymmetryPair(symmetry[k],symmetry[l]);
                    int nj = ibas_ab_sym[hkl][n][j];
                    int kl = ibas_ab_sym[hkl][k][l];
                    dum += u[g2aaoff[hkl] + (nj+gems_ab[hkl])*2*gems_ab[hkl] + (kl+gems_ab[hkl])];  // G2(nj,kl) dim
                    
                }


                A[offset + ijk*trip_aaa[h]+lmn] = dum;


            }
        }
        offset += trip_aaa[h]*trip_aaa[h];

    }

}

// T1 portion of A^T.y 
void v2RDMSolver::T1_constraints_ATu(double* A,double* u){

    // T1aab
    for (int h = 0; h < nirrep_; h++) {

        for (int ijk = 0; ijk < trip_aab[h]; ijk++) {

            int i = bas_aab_sym[h][ijk][0];
            int j = bas_aab_sym[h][ijk][1];
            int k = bas_aab_sym[h][ijk][2];

            for (int lmn = 0; lmn < trip_aab[h]; lmn++) {

                int l = bas_aab_sym[h][lmn][0];
                int m = bas_aab_sym[h][lmn][1];
                int n = bas_aab_sym[h][lmn][2];

                double dum = u[offset + ijk*trip_aab[h]+lmn]; 

                A[t1aaboff[h] + ijk*trip_aab[h]+lmn] -= dum; // - T1(ijk,lmn)

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    A[q2aaoff[hij] + ij*gems_aa[hij] + lm] += dum;  // Q2(ij,lm) dkn
                }

                if ( j == l ) {
                    int hki = SymmetryPair(symmetry[k],symmetry[i]);
                    int nm = ibas_ab_sym[hki][m][n];
                    int ki = ibas_ab_sym[hki][i][k];
                    A[d2aboff[hki] + nm*gems_ab[hki] + ki] -= dum;  // -D2(nm,ki) dlj
                }

                if ( l == i ) {
                    int hkj = SymmetryPair(symmetry[k],symmetry[j]);
                    int nm = ibas_ab_sym[hkj][m][n];
                    int kj = ibas_ab_sym[hkj][j][k];
                    A[d2aboff[hkj] + nm*gems_ab[hkj] + kj] += dum;  // D2(nm,kj) dli
                }

                if ( j == m ) {
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int kl = ibas_ab_sym[hni][k][l];
                    A[g2baoff[hni] + ni*gems_ab[hni] + kl] -= dum;  // -G2(ni,kl) djm
                    
                }

                if ( i == m ) {
                    int hkl = SymmetryPair(symmetry[k],symmetry[l]);
                    int nj = ibas_ab_sym[hkl][n][j];
                    int kl = ibas_ab_sym[hkl][k][l];
                    A[g2baoff[hkl] + nj*gems_ab[hkl] + kl] += dum;  // G2(nj,kl) dim
                    
                }
            }
        }
        offset += trip_aab[h]*trip_aab[h];

    }
    // T1bba
    for (int h = 0; h < nirrep_; h++) {

        for (int ijk = 0; ijk < trip_aab[h]; ijk++) {

            int i = bas_aab_sym[h][ijk][0];
            int j = bas_aab_sym[h][ijk][1];
            int k = bas_aab_sym[h][ijk][2];

            for (int lmn = 0; lmn < trip_aab[h]; lmn++) {

                int l = bas_aab_sym[h][lmn][0];
                int m = bas_aab_sym[h][lmn][1];
                int n = bas_aab_sym[h][lmn][2];

                double dum = u[offset + ijk*trip_aab[h]+lmn]; 

                A[t1bbaoff[h] + ijk*trip_aab[h]+lmn] -= dum; // - T1(ijk,lmn)

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    A[q2bboff[hij] + ij*gems_aa[hij] + lm] += dum;  // Q2(ij,lm) dkn
                }

                if ( j == l ) {
                    int hki = SymmetryPair(symmetry[k],symmetry[i]);
                    int nm = ibas_ab_sym[hki][n][m];
                    int ki = ibas_ab_sym[hki][k][i];
                    A[d2aboff[hki] + nm*gems_ab[hki] + ki] -= dum;  // -D2(nm,ki) dlj
                }

                if ( l == i ) {
                    int hkj = SymmetryPair(symmetry[k],symmetry[j]);
                    int nm = ibas_ab_sym[hkj][n][m];
                    int kj = ibas_ab_sym[hkj][k][j];
                    A[d2aboff[hkj] + nm*gems_ab[hkj] + kj] += dum;  // D2(nm,kj) dli
                }

                if ( j == m ) {
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int kl = ibas_ab_sym[hni][k][l];
                    A[g2aboff[hni] + ni*gems_ab[hni] + kl] -= dum;  // -G2(ni,kl) djm
                    
                }

                if ( i == m ) {
                    int hkl = SymmetryPair(symmetry[k],symmetry[l]);
                    int nj = ibas_ab_sym[hkl][n][j];
                    int kl = ibas_ab_sym[hkl][k][l];
                    A[g2aboff[hkl] + nj*gems_ab[hkl] + kl] += dum;  // G2(nj,kl) dim
                    
                }
            }
        }
        offset += trip_aab[h]*trip_aab[h];

    }
    // T1aaa
    for (int h = 0; h < nirrep_; h++) {

        for (int ijk = 0; ijk < trip_aaa[h]; ijk++) {

            int i = bas_aaa_sym[h][ijk][0];
            int j = bas_aaa_sym[h][ijk][1];
            int k = bas_aaa_sym[h][ijk][2];

            for (int lmn = 0; lmn < trip_aaa[h]; lmn++) {

                int l = bas_aaa_sym[h][lmn][0];
                int m = bas_aaa_sym[h][lmn][1];
                int n = bas_aaa_sym[h][lmn][2];

                double dum = u[offset + ijk*trip_aaa[h] + lmn];

                A[t1aaaoff[h] + ijk*trip_aaa[h]+lmn] -= dum; // - T1(ijk,lmn)

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    A[q2aaoff[hij] + ij*gems_aa[hij] + lm] += dum;  // Q2(ij,lm) dkn
                }

                if ( j == n ) {
                    int hik = SymmetryPair(symmetry[i],symmetry[k]);
                    int hlm = SymmetryPair(symmetry[l],symmetry[m]);
                    if ( hik == hlm ) {
                        int ik = ibas_aa_sym[hik][i][k];
                        int lm = ibas_aa_sym[hik][l][m];
                        A[q2aaoff[hik] + ik*gems_aa[hik] + lm] -= dum;  // -Q2(ik,lm) djn
                    }
                }

                if ( i == n ) {
                    int hjk = SymmetryPair(symmetry[j],symmetry[k]);
                    int hlm = SymmetryPair(symmetry[l],symmetry[m]);
                    if ( hjk == hlm ) {
                        int jk = ibas_aa_sym[hjk][j][k];
                        int lm = ibas_aa_sym[hjk][l][m];
                        A[q2aaoff[hjk] + jk*gems_aa[hjk] + lm] += dum;  // Q2(jk,lm) din
                    }
                }


                if ( l == k ) {
                    int hji = SymmetryPair(symmetry[j],symmetry[i]);
                    int hnm = SymmetryPair(symmetry[n],symmetry[m]);
                    if ( hji == hnm ) {
                        int ji = ibas_aa_sym[hji][j][i];
                        int nm = ibas_aa_sym[hji][n][m];
                        A[d2aaoff[hji] + nm*gems_aa[hji] + ji] += dum;  // D2(nm,ji) dlk
                    }
                }

                if ( j == l ) {
                    int hki = SymmetryPair(symmetry[k],symmetry[i]);
                    int ki = ibas_aa_sym[hki][k][i];
                    int nm = ibas_aa_sym[hki][n][m];
                    A[d2aaoff[hki] + nm*gems_aa[hki] + ki] -= dum;  // -D2(nm,ki) dlj
                }

                if ( l == i ) {
                    int hkj = SymmetryPair(symmetry[k],symmetry[j]);
                    int kj = ibas_aa_sym[hkj][k][j];
                    int nm = ibas_aa_sym[hkj][n][m];
                    A[d2aaoff[hkj] + nm*gems_aa[hkj] + kj] += dum;  // D2(nm,kj) dli
                }

                if ( k == m ) {
                    if ( j == l ) {
                        int h2 = symmetry[n];
                        int nn = n - pitzer_offset[h2];
                        int ii = i - pitzer_offset[h2];
                        A[d1aoff[h2] + nn*amopi_[h2]+ii] -= dum; // - D1(n,i) djl dkm
                    }
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int jl = ibas_ab_sym[hni][j][l];
                    A[g2aaoff[hni] + ni*2*gems_ab[hni] + jl] += dum;  // G2(ni,jl) dkm
                    
                }

                if ( j == m ) {
                    if ( k == l ) {
                        int h2 = symmetry[n];
                        int nn = n - pitzer_offset[h2];
                        int ii = i - pitzer_offset[h2];
                        A[d1aoff[h2] + nn*amopi_[h2]+ii] += dum; // D1(n,i) dkl djm
                    }
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int kl = ibas_ab_sym[hni][k][l];
                    A[g2aaoff[hni] + ni*2*gems_ab[hni] + kl] -= dum;  // -G2(ni,kl) djm
                    
                }

                if ( i == m ) {
                    if ( k == l ) {
                        int h2 = symmetry[n];
                        int nn = n - pitzer_offset[h2];
                        int jj = j - pitzer_offset[h2];
                        A[d1aoff[h2] + nn*amopi_[h2]+jj] -= dum; // - D1(n,j) dkl dim
                    }
                    int hkl = SymmetryPair(symmetry[k],symmetry[l]);
                    int nj = ibas_ab_sym[hkl][n][j];
                    int kl = ibas_ab_sym[hkl][k][l];
                    A[g2aaoff[hkl] + nj*2*gems_ab[hkl] + kl] += dum;  // G2(nj,kl) dim
                    
                }

            }
        }
        offset += trip_aaa[h]*trip_aaa[h];

    }
    // T1bbb
    for (int h = 0; h < nirrep_; h++) {

        for (int ijk = 0; ijk < trip_aaa[h]; ijk++) {

            int i = bas_aaa_sym[h][ijk][0];
            int j = bas_aaa_sym[h][ijk][1];
            int k = bas_aaa_sym[h][ijk][2];

            for (int lmn = 0; lmn < trip_aaa[h]; lmn++) {

                int l = bas_aaa_sym[h][lmn][0];
                int m = bas_aaa_sym[h][lmn][1];
                int n = bas_aaa_sym[h][lmn][2];

                double dum = u[offset + ijk*trip_aaa[h] + lmn];

                A[t1bbboff[h] + ijk*trip_aaa[h]+lmn] -= dum; // - T1(ijk,lmn)

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    A[q2bboff[hij] + ij*gems_aa[hij] + lm] += dum;  // Q2(ij,lm) dkn
                }

                if ( j == n ) {
                    int hik = SymmetryPair(symmetry[i],symmetry[k]);
                    int hlm = SymmetryPair(symmetry[l],symmetry[m]);
                    if ( hik == hlm ) {
                        int ik = ibas_aa_sym[hik][i][k];
                        int lm = ibas_aa_sym[hik][l][m];
                        A[q2bboff[hik] + ik*gems_aa[hik] + lm] -= dum;  // -Q2(ik,lm) djn
                    }
                }

                if ( i == n ) {
                    int hjk = SymmetryPair(symmetry[j],symmetry[k]);
                    int hlm = SymmetryPair(symmetry[l],symmetry[m]);
                    if ( hjk == hlm ) {
                        int jk = ibas_aa_sym[hjk][j][k];
                        int lm = ibas_aa_sym[hjk][l][m];
                        A[q2bboff[hjk] + jk*gems_aa[hjk] + lm] += dum;  // Q2(jk,lm) din
                    }
                }


                if ( l == k ) {
                    int hji = SymmetryPair(symmetry[j],symmetry[i]);
                    int hnm = SymmetryPair(symmetry[n],symmetry[m]);
                    if ( hji == hnm ) {
                        int ji = ibas_aa_sym[hji][j][i];
                        int nm = ibas_aa_sym[hji][n][m];
                        A[d2bboff[hji] + nm*gems_aa[hji] + ji] += dum;  // D2(nm,ji) dlk
                    }
                }

                if ( j == l ) {
                    int hki = SymmetryPair(symmetry[k],symmetry[i]);
                    int ki = ibas_aa_sym[hki][k][i];
                    int nm = ibas_aa_sym[hki][n][m];
                    A[d2bboff[hki] + nm*gems_aa[hki] + ki] -= dum;  // -D2(nm,ki) dlj
                }

                if ( l == i ) {
                    int hkj = SymmetryPair(symmetry[k],symmetry[j]);
                    int kj = ibas_aa_sym[hkj][k][j];
                    int nm = ibas_aa_sym[hkj][n][m];
                    A[d2bboff[hkj] + nm*gems_aa[hkj] + kj] += dum;  // D2(nm,kj) dli
                }

                if ( k == m ) {
                    if ( j == l ) {
                        int h2 = symmetry[n];
                        int nn = n - pitzer_offset[h2];
                        int ii = i - pitzer_offset[h2];
                        A[d1boff[h2] + nn*amopi_[h2]+ii] -= dum; // - D1(n,i) djl dkm
                    }
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int jl = ibas_ab_sym[hni][j][l];
                    A[g2aaoff[hni] + (ni+gems_ab[hni])*2*gems_ab[hni] + (jl+gems_ab[hni])] += dum;  // G2(ni,jl) dkm
                    
                }

                if ( j == m ) {
                    if ( k == l ) {
                        int h2 = symmetry[n];
                        int nn = n - pitzer_offset[h2];
                        int ii = i - pitzer_offset[h2];
                        A[d1boff[h2] + nn*amopi_[h2]+ii] += dum; // D1(n,i) dkl djm
                    }
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int kl = ibas_ab_sym[hni][k][l];
                    A[g2aaoff[hni] + (ni+gems_ab[hni])*2*gems_ab[hni] + (kl+gems_ab[hni])] -= dum;  // -G2(ni,kl) djm
                    
                }

                if ( i == m ) {
                    if ( k == l ) {
                        int h2 = symmetry[n];
                        int nn = n - pitzer_offset[h2];
                        int jj = j - pitzer_offset[h2];
                        A[d1boff[h2] + nn*amopi_[h2]+jj] -= dum; // - D1(n,j) dkl dim
                    }
                    int hkl = SymmetryPair(symmetry[k],symmetry[l]);
                    int nj = ibas_ab_sym[hkl][n][j];
                    int kl = ibas_ab_sym[hkl][k][l];
                    A[g2aaoff[hkl] + (nj+gems_ab[hkl])*2*gems_ab[hkl] + (kl+gems_ab[hkl])] += dum;  // G2(nj,kl) dim
                    
                }

            }
        }
        offset += trip_aaa[h]*trip_aaa[h];

    }
}

} // end namespaces
