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

// F3 portion of A.u
void v2RDMSolver::F3_constraints_Au(double * A,double * u){

    // F3aab
    for (int h = 0; h < nirrep_; h++) {

        #pragma omp parallel for schedule (static)
        for (int ijk = 0; ijk < trip_aab[h]; ijk++) {

            int i = bas_aab_sym[h][ijk][0];
            int j = bas_aab_sym[h][ijk][1];
            int k = bas_aab_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aab[h]; lmn++) {

                int l = bas_aab_sym[h][lmn][0];
                int m = bas_aab_sym[h][lmn][1];
                int n = bas_aab_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int ijn = ibas_aab_sym[hijn][i][j][n];
                int lmk = ibas_aab_sym[hijn][l][m][k];

                double dum = -u[f3aaboff[h] + ijk*trip_aab[h]+lmn]; // - F3(ijk,lmn)
                dum       +=  u[d3aaboff[hijn] + ijn*trip_aab[hijn]+lmk]; // - D3(ijn,lmk)

                if ( j == m && i == l ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    dum += u[d1boff[h2] + nn*amopi_[h2]+kk]; // + D1(n,k) djm dil
                }
                if ( j == l && i == m ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    dum -= u[d1boff[h2] + nn*amopi_[h2]+kk]; // - D1(n,k) djl dim
                }

                if ( i == l ) {
                    int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                    int nj = ibas_ab_sym[hnj][j][n];
                    int km = ibas_ab_sym[hnj][m][k];
                    dum -= u[d2aboff[hnj] + nj*gems_ab[hnj]+km]; // - D2(nj,km) dil
                }
                if ( j == l ) {
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][i][n];
                    int km = ibas_ab_sym[hni][m][k];
                    dum += u[d2aboff[hni] + ni*gems_ab[hni]+km]; // D2(ni,km) djl
                }
                if ( i == m ) {
                    int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                    int nj = ibas_ab_sym[hnj][j][n];
                    int kl = ibas_ab_sym[hnj][l][k];
                    dum += u[d2aboff[hnj] + nj*gems_ab[hnj]+kl]; // D2(nj,kl) dim
                }
                if ( j == m ) {
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][i][n];
                    int kl = ibas_ab_sym[hni][l][k];
                    dum -= u[d2aboff[hni] + ni*gems_ab[hni]+kl]; // -D2(ni,kl) djm
                }

                A[offset + ijk*trip_aab[h]+lmn] = dum;

            }
        }
        offset += trip_aab[h]*trip_aab[h];
    }

    // F3bba
    for (int h = 0; h < nirrep_; h++) {

        #pragma omp parallel for schedule (static)
        for (int ijk = 0; ijk < trip_aab[h]; ijk++) {

            int i = bas_aab_sym[h][ijk][0];
            int j = bas_aab_sym[h][ijk][1];
            int k = bas_aab_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aab[h]; lmn++) {

                int l = bas_aab_sym[h][lmn][0];
                int m = bas_aab_sym[h][lmn][1];
                int n = bas_aab_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int ijn = ibas_aab_sym[hijn][i][j][n];
                int lmk = ibas_aab_sym[hijn][l][m][k];

                double dum = -u[f3bbaoff[h] + ijk*trip_aab[h]+lmn]; // - F3(ijk,lmn)
                dum       +=  u[d3bbaoff[hijn] + ijn*trip_aab[hijn]+lmk]; // - D3(ijn,lmk)

                if ( j == m && i == l ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    dum += u[d1aoff[h2] + nn*amopi_[h2]+kk]; // + D1(n,k) djm dil
                }
                if ( j == l && i == m ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    dum -= u[d1aoff[h2] + nn*amopi_[h2]+kk]; // - D1(n,k) djl dim
                }

                if ( i == l ) {
                    int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                    int nj = ibas_ab_sym[hnj][n][j];
                    int km = ibas_ab_sym[hnj][k][m];
                    dum -= u[d2aboff[hnj] + nj*gems_ab[hnj]+km]; // - D2(nj,km) dil
                }
                if ( j == l ) {
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int km = ibas_ab_sym[hni][k][m];
                    dum += u[d2aboff[hni] + ni*gems_ab[hni]+km]; // D2(ni,km) djl
                }
                if ( i == m ) {
                    int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                    int nj = ibas_ab_sym[hnj][n][j];
                    int kl = ibas_ab_sym[hnj][k][l];
                    dum += u[d2aboff[hnj] + nj*gems_ab[hnj]+kl]; // D2(nj,kl) dim
                }
                if ( j == m ) {
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int kl = ibas_ab_sym[hni][k][l];
                    dum -= u[d2aboff[hni] + ni*gems_ab[hni]+kl]; // -D2(ni,kl) djm
                }

                A[offset + ijk*trip_aab[h]+lmn] = dum;

            }
        }
        offset += trip_aab[h]*trip_aab[h];
    }

    // big block 1: F3aaa + F3abb
    for (int h = 0; h < nirrep_; h++) {

        // F3aaa/aaa
        #pragma omp parallel for schedule (static)
        for (int ijk = 0; ijk < trip_aab[h]; ijk++) {

            int i = bas_aab_sym[h][ijk][0];
            int j = bas_aab_sym[h][ijk][1];
            int k = bas_aab_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aab[h]; lmn++) {

                int l = bas_aab_sym[h][lmn][0];
                int m = bas_aab_sym[h][lmn][1];
                int n = bas_aab_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int id  = ijk*(trip_aab[h]+trip_aba[h])+lmn;

                double dum = -u[f3aaaoff[h] + id]; // - F3(ijk,lmn)

                if ( i != n && j != n && l != k && m != k ) {

                    int ijn = ibas_aaa_sym[hijn][i][j][n];
                    int lmk = ibas_aaa_sym[hijn][l][m][k];

                    // sign ... i < j, l < m guaranteed
                    //int sg = 1;
                    //if ( n < i && n < j ) {
                    //   // nothing
                    //}else if ( n > i && n > j ) {
                    //   // nothing
                    //}else if ( n > i && n < j ) {
                    //    sg = -sg;
                    //}else if ( n < i && n > j ) {
                    //    sg = -sg; // never will happen
                    //}
                    int sg = 1;
                    if ( i < n && n < j ) {
                        sg = -sg;
                    }
                    if ( l < k && k < m ) {
                        sg = -sg;
                    }

                    dum +=  sg * u[d3aaaoff[hijn] + ijn*trip_aaa[hijn]+lmk]; // - D3(ijn,lmk)

                }

                if ( j == m && i == l ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    dum += u[d1aoff[h2] + nn*amopi_[h2]+kk]; // + D1(n,k) djm dil
                }
                if ( j == l && i == m ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    dum -= u[d1aoff[h2] + nn*amopi_[h2]+kk]; // - D1(n,k) djl dim
                }

                if ( i == l ) {
                    if ( n != j && k != m ) {
                        int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                        int nj = ibas_aa_sym[hnj][n][j];
                        int km = ibas_aa_sym[hnj][k][m];

                        int s = 1;
                        if ( n > j ) s = -s;
                        if ( k > m ) s = -s;

                        dum -= s * u[d2aaoff[hnj] + nj*gems_aa[hnj]+km]; // - D2(nj,km) dil
                    }
                }
                if ( j == l ) {
                    if ( n != i && k != m ) {
                        int hni = SymmetryPair(symmetry[n],symmetry[i]);
                        int ni = ibas_aa_sym[hni][n][i];
                        int km = ibas_aa_sym[hni][k][m];

                        int s = 1;
                        if ( n > i ) s = -s;
                        if ( k > m ) s = -s;

                        dum += s * u[d2aaoff[hni] + ni*gems_aa[hni]+km]; // D2(ni,km) djl
                    }
                }
                if ( i == m ) {
                    if ( n != j && k != l ) {
                        int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                        int nj = ibas_aa_sym[hnj][n][j];
                        int kl = ibas_aa_sym[hnj][k][l];

                        int s = 1;
                        if ( n > j ) s = -s;
                        if ( k > l ) s = -s;

                        dum += s * u[d2aaoff[hnj] + nj*gems_aa[hnj]+kl]; // D2(nj,kl) dim
                    }
                }
                if ( j == m ) {
                    if ( n != i && k != l ) {
                        int hni = SymmetryPair(symmetry[n],symmetry[i]);
                        int ni = ibas_aa_sym[hni][n][i];
                        int kl = ibas_aa_sym[hni][k][l];

                        int s = 1;
                        if ( n > i ) s = -s;
                        if ( k > l ) s = -s;

                        dum -= s * u[d2aaoff[hni] + ni*gems_aa[hni]+kl]; // -D2(ni,kl) djm
                    }
                }

                A[offset + id] = dum;

            }
        }
        // F3aaa/abb
        #pragma omp parallel for schedule (static)
        for (int ijk = 0; ijk < trip_aab[h]; ijk++) {

            int i = bas_aab_sym[h][ijk][0];
            int j = bas_aab_sym[h][ijk][1];
            int k = bas_aab_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aba[h]; lmn++) {

                int l = bas_aba_sym[h][lmn][0];
                int m = bas_aba_sym[h][lmn][1];
                int n = bas_aba_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int id = ijk*(trip_aab[h]+trip_aba[h])+(lmn+trip_aab[h]);

                double dum = -u[f3aaaoff[h] + id]; // - F3(ijk,lmn)

                if ( i != j && l != k ) {

                    int ijn = ibas_aab_sym[hijn][i][j][n];
                    int lkm = ibas_aab_sym[hijn][l][k][m];

                    int sg = 1;
                    if ( l > k ) {
                        sg = -sg;
                    }

                    dum -=  sg * u[d3aaboff[hijn] + ijn*trip_aab[hijn]+lkm]; // + D3(ijn,lkm)
                }

                if ( i == l ) {
                    int hjn = SymmetryPair(symmetry[j],symmetry[n]);
                    int jn = ibas_ab_sym[hjn][j][n];
                    int km = ibas_ab_sym[hjn][k][m];
                    dum += u[d2aboff[hjn]+jn*gems_ab[hjn]+km]; // D2(jn,km) dil
                }
                if ( j == l ) {
                    int hin = SymmetryPair(symmetry[i],symmetry[n]);
                    int in = ibas_ab_sym[hin][i][n];
                    int km = ibas_ab_sym[hin][k][m];
                    dum -= u[d2aboff[hin]+in*gems_ab[hin]+km]; // -D2(in,km) djl
                }

                A[offset + id] = dum;

            }
        }

        // F3abb/aaa
        #pragma omp parallel for schedule (static)
        for (int ijk = 0; ijk < trip_aba[h]; ijk++) {

            int i = bas_aba_sym[h][ijk][0];
            int j = bas_aba_sym[h][ijk][1];
            int k = bas_aba_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aab[h]; lmn++) {

                int l = bas_aab_sym[h][lmn][0];
                int m = bas_aab_sym[h][lmn][1];
                int n = bas_aab_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int id = (ijk+trip_aab[h])*(trip_aab[h]+trip_aba[h])+lmn;

                double dum = -u[f3aaaoff[h] + id]; // - F3(ijk,lmn)

                if ( i != n && l != m ) {
                    
                    int inj = ibas_aab_sym[hijn][i][n][j];
                    int lmk = ibas_aab_sym[hijn][l][m][k];
                    
                    int sg = 1;
                    if ( i > n ) {
                        sg = -sg;
                    }

                    dum -=  sg * u[d3aaboff[hijn] + inj*trip_aab[hijn]+lmk]; // + D3(inj,lmk)
                }  

                if ( i == l ) {
                    int hjn = SymmetryPair(symmetry[j],symmetry[n]);
                    int jn = ibas_ab_sym[hjn][n][j];
                    int km = ibas_ab_sym[hjn][m][k];
                    dum += u[d2aboff[hjn]+jn*gems_ab[hjn]+km]; // D2(jn,km) dil
                }
                if ( i == m ) {
                    int hnj = SymmetryPair(symmetry[j],symmetry[n]);
                    int nj = ibas_ab_sym[hnj][n][j];
                    int lk = ibas_ab_sym[hnj][l][k];
                    dum -= u[d2aboff[hnj]+nj*gems_ab[hnj]+lk]; // -D2(in,km) djl
                }

                A[offset + id] = dum;

            }
        }

        // F3abb/abb
        #pragma omp parallel for schedule (static)
        for (int ijk = 0; ijk < trip_aba[h]; ijk++) {

            int i = bas_aba_sym[h][ijk][0];
            int j = bas_aba_sym[h][ijk][1];
            int k = bas_aba_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aba[h]; lmn++) {

                int l = bas_aba_sym[h][lmn][0];
                int m = bas_aba_sym[h][lmn][1];
                int n = bas_aba_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int id = (ijk+trip_aab[h])*(trip_aab[h]+trip_aba[h])+(lmn+trip_aab[h]);

                double dum = -u[f3aaaoff[h] + id]; // - F3(ijk,lmn)

                if ( j != n && m != k ) {
                    int jni = ibas_aab_sym[hijn][j][n][i];
                    int mkl = ibas_aab_sym[hijn][m][k][l];
                    int sg = 1;
                    if ( j > n ) sg = -sg;
                    if ( m > k ) sg = -sg;
                    dum += sg * u[d3bbaoff[hijn] + jni*trip_aab[hijn]+mkl]; // - D3(jni,mkl)
                }

                if ( j == m && i == l ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    dum += u[d1boff[h2] + nn*amopi_[h2]+kk]; // + D1(n,k) djm dil
                }

                if ( i == l ) {
                    if ( n != j && k != m ) {
                        int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                        int nj = ibas_aa_sym[hnj][n][j];
                        int km = ibas_aa_sym[hnj][k][m];

                        int s = 1;
                        if ( n > j ) s = -s;
                        if ( k > m ) s = -s;

                        dum -= s * u[d2bboff[hnj] + nj*gems_aa[hnj]+km]; // - D2(nj,km) dil
                    }
                }
                if ( j == m ) {
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][i][n];
                    int kl = ibas_ab_sym[hni][l][k];
                    dum -= u[d2aboff[hni] + ni*gems_ab[hni]+kl]; // -D2(ni,kl) djm
                }

                A[offset + id] = dum;

            }
        }
        offset += (trip_aba[h]+trip_aab[h])*(trip_aba[h]+trip_aab[h]);
    }

    // big block 2: F3bbb + F3baa
    for (int h = 0; h < nirrep_; h++) {

        // F3bbb/bbb
        #pragma omp parallel for schedule (static)
        for (int ijk = 0; ijk < trip_aab[h]; ijk++) {

            int i = bas_aab_sym[h][ijk][0];
            int j = bas_aab_sym[h][ijk][1];
            int k = bas_aab_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aab[h]; lmn++) {

                int l = bas_aab_sym[h][lmn][0];
                int m = bas_aab_sym[h][lmn][1];
                int n = bas_aab_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int id = ijk*(trip_aab[h]+trip_aba[h])+lmn;
                //int id = ijk*trip_aab[h]+lmn;

                double dum = -u[f3bbboff[h] + id]; // - F3(ijk,lmn)

                if ( i != n && j != n && l != k && m != k ) {

                    int ijn = ibas_aaa_sym[hijn][i][j][n];
                    int lmk = ibas_aaa_sym[hijn][l][m][k];
                    
                    // sign ... i < j, l < m guaranteed
                    //int sg = 1;
                    //if ( n < i && n < j ) {
                    //   // nothing
                    //}else if ( n > i && n > j ) {
                    //   // nothing
                    //}else if ( n > i && n < j ) {
                    //    sg = -sg;
                    //}else if ( n < i && n > j ) {
                    //    sg = -sg; // never will happen
                    //} 
                    int sg = 1;
                    if ( i < n && n < j ) {
                        sg = -sg;
                    }
                    if ( l < k && k < m ) {
                        sg = -sg;
                    }

                    dum +=  sg * u[d3bbboff[hijn] + ijn*trip_aaa[hijn]+lmk]; // - D3(ijn,lmk)

                }

                if ( j == m && i == l ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    dum += u[d1boff[h2] + nn*amopi_[h2]+kk]; // + D1(n,k) djm dil
                }
                if ( j == l && i == m ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    dum -= u[d1boff[h2] + nn*amopi_[h2]+kk]; // - D1(n,k) djl dim
                }

                if ( i == l ) {
                    if ( n != j && k != m ) {
                        int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                        int nj = ibas_aa_sym[hnj][n][j];
                        int km = ibas_aa_sym[hnj][k][m];

                        int s = 1;
                        if ( n > j ) s = -s;
                        if ( k > m ) s = -s;

                        dum -= s * u[d2bboff[hnj] + nj*gems_aa[hnj]+km]; // - D2(nj,km) dil
                    }
                }
                if ( j == l ) {
                    if ( n != i && k != m ) {
                        int hni = SymmetryPair(symmetry[n],symmetry[i]);
                        int ni = ibas_aa_sym[hni][n][i];
                        int km = ibas_aa_sym[hni][k][m];

                        int s = 1;
                        if ( n > i ) s = -s;
                        if ( k > m ) s = -s;

                        dum += s * u[d2bboff[hni] + ni*gems_aa[hni]+km]; // D2(ni,km) djl
                    }
                }
                if ( i == m ) {
                    if ( n != j && k != l ) {
                        int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                        int nj = ibas_aa_sym[hnj][n][j];
                        int kl = ibas_aa_sym[hnj][k][l];

                        int s = 1;
                        if ( n > j ) s = -s;
                        if ( k > l ) s = -s;

                        dum += s * u[d2bboff[hnj] + nj*gems_aa[hnj]+kl]; // D2(nj,kl) dim
                    }
                }
                if ( j == m ) {
                    if ( n != i && k != l ) {
                        int hni = SymmetryPair(symmetry[n],symmetry[i]);
                        int ni = ibas_aa_sym[hni][n][i];
                        int kl = ibas_aa_sym[hni][k][l];

                        int s = 1;
                        if ( n > i ) s = -s;
                        if ( k > l ) s = -s;

                        dum -= s * u[d2bboff[hni] + ni*gems_aa[hni]+kl]; // -D2(ni,kl) djm
                    }
                }

                A[offset + id] = dum;

            }
        }
        // F3bbb/baa
        #pragma omp parallel for schedule (static)
        for (int ijk = 0; ijk < trip_aab[h]; ijk++) {

            int i = bas_aab_sym[h][ijk][0];
            int j = bas_aab_sym[h][ijk][1];
            int k = bas_aab_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aba[h]; lmn++) {

                int l = bas_aba_sym[h][lmn][0];
                int m = bas_aba_sym[h][lmn][1];
                int n = bas_aba_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int id = ijk*(trip_aab[h]+trip_aba[h])+(lmn+trip_aab[h]);

                double dum = -u[f3bbboff[h] + id]; // - F3(ijk,lmn)

                if ( i != j && l != k ) {
                    int ijn = ibas_aab_sym[hijn][i][j][n];
                    int lkm = ibas_aab_sym[hijn][l][k][m];
                    int sg = 1;
                    if ( i > j ) sg = -sg;
                    if ( l > k ) sg = -sg;
                    dum -= sg * u[d3bbaoff[hijn] + ijn*trip_aab[hijn]+lkm]; // + D3(ijn,lkm)
                } 

                if ( i == l ) {
                    int hjn = SymmetryPair(symmetry[j],symmetry[n]);
                    int jn = ibas_ab_sym[hjn][n][j];
                    int km = ibas_ab_sym[hjn][m][k];
                    dum += u[d2aboff[hjn]+jn*gems_ab[hjn]+km]; // D2(jn,km) dil
                }
                if ( j == l ) {
                    int hin = SymmetryPair(symmetry[i],symmetry[n]);
                    int in = ibas_ab_sym[hin][n][i];
                    int km = ibas_ab_sym[hin][m][k];
                    dum -= u[d2aboff[hin]+in*gems_ab[hin]+km]; // -D2(in,km) djl
                }

                A[offset + id] = dum;

            }
        }

        // F3baa/bbb
        #pragma omp parallel for schedule (static)
        for (int ijk = 0; ijk < trip_aba[h]; ijk++) {

            int i = bas_aba_sym[h][ijk][0];
            int j = bas_aba_sym[h][ijk][1];
            int k = bas_aba_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aab[h]; lmn++) {

                int l = bas_aab_sym[h][lmn][0];
                int m = bas_aab_sym[h][lmn][1];
                int n = bas_aab_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int id = (ijk+trip_aab[h])*(trip_aab[h]+trip_aba[h])+lmn;

                double dum = -u[f3bbboff[h] + id]; // - F3(ijk,lmn)

                if ( i != n && l != m ) {
                    int inj = ibas_aab_sym[hijn][i][n][j];
                    int lmk = ibas_aab_sym[hijn][l][m][k];
                    int sg = 1;
                    if ( i > n ) sg = -sg;
                    if ( l > m ) sg = -sg;
                    dum -= sg * u[d3bbaoff[hijn] + inj*trip_aab[hijn]+lmk]; // + D3(inj,lmk)
                }

                if ( i == l ) {
                    int hjn = SymmetryPair(symmetry[j],symmetry[n]);
                    int jn = ibas_ab_sym[hjn][j][n];
                    int km = ibas_ab_sym[hjn][k][m];
                    dum += u[d2aboff[hjn]+jn*gems_ab[hjn]+km]; // D2(jn,km) dil
                }
                if ( i == m ) {
                    int hnj = SymmetryPair(symmetry[j],symmetry[n]);
                    int nj = ibas_ab_sym[hnj][j][n];
                    int lk = ibas_ab_sym[hnj][k][l];
                    dum -= u[d2aboff[hnj]+nj*gems_ab[hnj]+lk]; // -D2(in,km) djl
                }

                A[offset + id] = dum;

            }
        }
        // F3baa/baa
        #pragma omp parallel for schedule (static)
        for (int ijk = 0; ijk < trip_aba[h]; ijk++) {

            int i = bas_aba_sym[h][ijk][0];
            int j = bas_aba_sym[h][ijk][1];
            int k = bas_aba_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aba[h]; lmn++) {

                int l = bas_aba_sym[h][lmn][0];
                int m = bas_aba_sym[h][lmn][1];
                int n = bas_aba_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int id = (ijk+trip_aab[h])*(trip_aab[h]+trip_aba[h])+(lmn+trip_aab[h]);

                double dum = -u[f3bbboff[h] + id]; // - F3(ijk,lmn)

                if ( j != n && m != k ) {
                    int jni = ibas_aab_sym[hijn][j][n][i];
                    int mkl = ibas_aab_sym[hijn][m][k][l];
                    int sg = 1;
                    if ( j > n ) sg = -sg;
                    if ( m > k ) sg = -sg;
                    dum += sg * u[d3aaboff[hijn] + jni*trip_aab[hijn]+mkl]; // - D3(jni,mkl)
                }

                if ( j == m && i == l ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    dum += u[d1aoff[h2] + nn*amopi_[h2]+kk]; // + D1(n,k) djm dil
                }

                if ( i == l ) {
                    if ( n != j && k != m ) {
                        int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                        int nj = ibas_aa_sym[hnj][n][j];
                        int km = ibas_aa_sym[hnj][k][m];

                        int s = 1;
                        if ( n > j ) s = -s;
                        if ( k > m ) s = -s;

                        dum -= s * u[d2aaoff[hnj] + nj*gems_aa[hnj]+km]; // - D2(nj,km) dil
                    }
                }
                if ( j == m ) {
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int kl = ibas_ab_sym[hni][k][l];
                    dum -= u[d2aboff[hni] + ni*gems_ab[hni]+kl]; // -D2(ni,kl) djm
                }

                A[offset + id] = dum;

            }
        }
        offset += (trip_aba[h]+trip_aab[h])*(trip_aba[h]+trip_aab[h]);
    }

}

// F3 portion of AT.u
void v2RDMSolver::F3_constraints_ATu(double * A,double * u){

    // F3aab
    for (int h = 0; h < nirrep_; h++) {

        for (int ijk = 0; ijk < trip_aab[h]; ijk++) {

            int i = bas_aab_sym[h][ijk][0];
            int j = bas_aab_sym[h][ijk][1];
            int k = bas_aab_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aab[h]; lmn++) {

                int l = bas_aab_sym[h][lmn][0];
                int m = bas_aab_sym[h][lmn][1];
                int n = bas_aab_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int ijn = ibas_aab_sym[hijn][i][j][n];
                int lmk = ibas_aab_sym[hijn][l][m][k];

                double dum = u[offset + ijk*trip_aab[h]+lmn];

                A[f3aaboff[h] + ijk*trip_aab[h]+lmn] -= dum; // - F3(ijk,lmn)
                A[d3aaboff[hijn] + ijn*trip_aab[hijn]+lmk] += dum; // - D3(ijn,lmk)

                if ( j == m && i == l ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    A[d1boff[h2] + nn*amopi_[h2]+kk] += dum; // + D1(n,k) djm dil
                }   
                if ( j == l && i == m ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    A[d1boff[h2] + nn*amopi_[h2]+kk] -= dum; // - D1(n,k) djl dim
                }   
                
                if ( i == l ) {
                    int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                    int nj = ibas_ab_sym[hnj][j][n];
                    int km = ibas_ab_sym[hnj][m][k];
                    A[d2aboff[hnj] + nj*gems_ab[hnj]+km] -= dum; // - D2(nj,km) dil
                }   
                if ( j == l ) {
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][i][n];
                    int km = ibas_ab_sym[hni][m][k];
                    A[d2aboff[hni] + ni*gems_ab[hni]+km] += dum; // D2(ni,km) djl
                }   
                if ( i == m ) {
                    int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                    int nj = ibas_ab_sym[hnj][j][n];
                    int kl = ibas_ab_sym[hnj][l][k];
                    A[d2aboff[hnj] + nj*gems_ab[hnj]+kl] += dum; // D2(nj,kl) dim
                }   
                if ( j == m ) {
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][i][n];
                    int kl = ibas_ab_sym[hni][l][k];
                    A[d2aboff[hni] + ni*gems_ab[hni]+kl] -= dum; // -D2(ni,kl) djm
                }

            }
        }
        offset += trip_aab[h]*trip_aab[h];
    }

    // F3bba
    for (int h = 0; h < nirrep_; h++) {

        for (int ijk = 0; ijk < trip_aab[h]; ijk++) {

            int i = bas_aab_sym[h][ijk][0];
            int j = bas_aab_sym[h][ijk][1];
            int k = bas_aab_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aab[h]; lmn++) {

                int l = bas_aab_sym[h][lmn][0];
                int m = bas_aab_sym[h][lmn][1];
                int n = bas_aab_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int ijn = ibas_aab_sym[hijn][i][j][n];
                int lmk = ibas_aab_sym[hijn][l][m][k];

                double dum = u[offset + ijk*trip_aab[h]+lmn];

                A[f3bbaoff[h] + ijk*trip_aab[h]+lmn] -= dum; // - F3(ijk,lmn)
                A[d3bbaoff[hijn] + ijn*trip_aab[hijn]+lmk] += dum; // - D3(ijn,lmk)

                if ( j == m && i == l ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    A[d1aoff[h2] + nn*amopi_[h2]+kk] += dum; // + D1(n,k) djm dil
                }   
                if ( j == l && i == m ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    A[d1aoff[h2] + nn*amopi_[h2]+kk] -= dum; // - D1(n,k) djl dim
                }   
                
                if ( i == l ) {
                    int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                    int nj = ibas_ab_sym[hnj][n][j];
                    int km = ibas_ab_sym[hnj][k][m];
                    A[d2aboff[hnj] + nj*gems_ab[hnj]+km] -= dum; // - D2(nj,km) dil
                }   
                if ( j == l ) {
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int km = ibas_ab_sym[hni][k][m];
                    A[d2aboff[hni] + ni*gems_ab[hni]+km] += dum; // D2(ni,km) djl
                }   
                if ( i == m ) {
                    int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                    int nj = ibas_ab_sym[hnj][n][j];
                    int kl = ibas_ab_sym[hnj][k][l];
                    A[d2aboff[hnj] + nj*gems_ab[hnj]+kl] += dum; // D2(nj,kl) dim
                }   
                if ( j == m ) {
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int kl = ibas_ab_sym[hni][k][l];
                    A[d2aboff[hni] + ni*gems_ab[hni]+kl] -= dum; // -D2(ni,kl) djm
                }

            }
        }
        offset += trip_aab[h]*trip_aab[h];
    }

    // big block 1: F3aaa + F3abb
    for (int h = 0; h < nirrep_; h++) {

        // F3aaa/aaa
        for (int ijk = 0; ijk < trip_aab[h]; ijk++) {

            int i = bas_aab_sym[h][ijk][0];
            int j = bas_aab_sym[h][ijk][1];
            int k = bas_aab_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aab[h]; lmn++) {

                int l = bas_aab_sym[h][lmn][0];
                int m = bas_aab_sym[h][lmn][1];
                int n = bas_aab_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int id  = ijk*(trip_aab[h]+trip_aba[h])+lmn;

                double dum = u[offset + id];

                A[f3aaaoff[h] + id] -= dum; // - F3(ijk,lmn)

                if ( i != n && j != n && l != k && m != k ) {

                    int ijn = ibas_aaa_sym[hijn][i][j][n];
                    int lmk = ibas_aaa_sym[hijn][l][m][k];

                    // sign ... i < j, l < m guaranteed
                    //int sg = 1;
                    //if ( n < i && n < j ) {
                    //   // nothing
                    //}else if ( n > i && n > j ) {
                    //   // nothing
                    //}else if ( n > i && n < j ) {
                    //    sg = -sg;
                    //}else if ( n < i && n > j ) {
                    //    sg = -sg; // never will happen
                    //}
                    int sg = 1;
                    if ( i < n && n < j ) {
                        sg = -sg;
                    }
                    if ( l < k && k < m ) {
                        sg = -sg;
                    }

                    A[d3aaaoff[hijn] + ijn*trip_aaa[hijn]+lmk] += sg * dum; // - D3(ijn,lmk)

                }

                if ( j == m && i == l ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    A[d1aoff[h2] + nn*amopi_[h2]+kk] += dum; // + D1(n,k) djm dil
                }
                if ( j == l && i == m ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    A[d1aoff[h2] + nn*amopi_[h2]+kk] -= dum; // - D1(n,k) djl dim
                }

                if ( i == l ) {
                    if ( n != j && k != m ) {
                        int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                        int nj = ibas_aa_sym[hnj][n][j];
                        int km = ibas_aa_sym[hnj][k][m];
                
                        int s = 1;
                        if ( n > j ) s = -s;
                        if ( k > m ) s = -s;
                    
                        A[d2aaoff[hnj] + nj*gems_aa[hnj]+km] -= s * dum; // - D2(nj,km) dil
                    }
                }
                if ( j == l ) {
                    if ( n != i && k != m ) {
                        int hni = SymmetryPair(symmetry[n],symmetry[i]);
                        int ni = ibas_aa_sym[hni][n][i];
                        int km = ibas_aa_sym[hni][k][m];
                
                        int s = 1;
                        if ( n > i ) s = -s;
                        if ( k > m ) s = -s;
                        
                        A[d2aaoff[hni] + ni*gems_aa[hni]+km] += s * dum; // D2(ni,km) djl
                    }   
                }
                if ( i == m ) { 
                    if ( n != j && k != l ) {
                        int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                        int nj = ibas_aa_sym[hnj][n][j];
                        int kl = ibas_aa_sym[hnj][k][l];
                    
                        int s = 1;
                        if ( n > j ) s = -s;
                        if ( k > l ) s = -s;
                        
                        A[d2aaoff[hnj] + nj*gems_aa[hnj]+kl] += s * dum; // D2(nj,kl) dim
                    }   
                }
                if ( j == m ) { 
                    if ( n != i && k != l ) {
                        int hni = SymmetryPair(symmetry[n],symmetry[i]);
                        int ni = ibas_aa_sym[hni][n][i];
                        int kl = ibas_aa_sym[hni][k][l];
                    
                        int s = 1;
                        if ( n > i ) s = -s;
                        if ( k > l ) s = -s;
                        
                        A[d2aaoff[hni] + ni*gems_aa[hni]+kl] -= s * dum; // -D2(ni,kl) djm
                    }   
                }

            }
        }
        // F3aaa/abb
        for (int ijk = 0; ijk < trip_aab[h]; ijk++) {

            int i = bas_aab_sym[h][ijk][0];
            int j = bas_aab_sym[h][ijk][1];
            int k = bas_aab_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aba[h]; lmn++) {

                int l = bas_aba_sym[h][lmn][0];
                int m = bas_aba_sym[h][lmn][1];
                int n = bas_aba_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int id = ijk*(trip_aab[h]+trip_aba[h])+(lmn+trip_aab[h]);

                double dum = u[offset + id]; // - F3(ijk,lmn)

                A[f3aaaoff[h] + id] -= dum; // - F3(ijk,lmn)

                if ( i != j && l != k ) {

                    int ijn = ibas_aab_sym[hijn][i][j][n];
                    int lkm = ibas_aab_sym[hijn][l][k][m];

                    int sg = 1;
                    if ( k < l ) {
                        sg = -sg;
                    }

                    A[d3aaboff[hijn] + ijn*trip_aab[hijn]+lkm] -= sg * dum; // + D3(ijn,lkm)
                }

                if ( i == l ) {
                    int hjn = SymmetryPair(symmetry[j],symmetry[n]);
                    int jn = ibas_ab_sym[hjn][j][n];
                    int km = ibas_ab_sym[hjn][k][m];
                    A[d2aboff[hjn]+jn*gems_ab[hjn]+km] += dum; // D2(jn,km) dil
                }
                if ( j == l ) {
                    int hin = SymmetryPair(symmetry[i],symmetry[n]);
                    int in = ibas_ab_sym[hin][i][n];
                    int km = ibas_ab_sym[hin][k][m];
                    A[d2aboff[hin]+in*gems_ab[hin]+km] -= dum; // -D2(in,km) djl
                }

            }
        }

        // F3abb/aaa
        for (int ijk = 0; ijk < trip_aba[h]; ijk++) {

            int i = bas_aba_sym[h][ijk][0];
            int j = bas_aba_sym[h][ijk][1];
            int k = bas_aba_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aab[h]; lmn++) {

                int l = bas_aab_sym[h][lmn][0];
                int m = bas_aab_sym[h][lmn][1];
                int n = bas_aab_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int id = (ijk+trip_aab[h])*(trip_aab[h]+trip_aba[h])+lmn;

                double dum = u[offset + id];

                A[f3aaaoff[h] + id] -= dum; // - F3(ijk,lmn)

                if ( i != n && l != m ) {
                    
                    int inj = ibas_aab_sym[hijn][i][n][j];
                    int lmk = ibas_aab_sym[hijn][l][m][k];
                    
                    int sg = 1;
                    if ( n < i ) {
                        sg = -sg;
                    }

                    A[d3aaboff[hijn] + inj*trip_aab[hijn]+lmk] -= sg * dum; // + D3(inj,lmk)
                }  

                if ( i == l ) {
                    int hjn = SymmetryPair(symmetry[j],symmetry[n]);
                    int jn = ibas_ab_sym[hjn][n][j];
                    int km = ibas_ab_sym[hjn][m][k];
                    A[d2aboff[hjn]+jn*gems_ab[hjn]+km] += dum; // D2(jn,km) dil
                }
                if ( i == m ) {
                    int hnj = SymmetryPair(symmetry[j],symmetry[n]);
                    int nj = ibas_ab_sym[hnj][n][j];
                    int lk = ibas_ab_sym[hnj][l][k];
                    A[d2aboff[hnj]+nj*gems_ab[hnj]+lk] -= dum; // -D2(in,km) djl
                }

            }
        }

        // F3abb/abb
        for (int ijk = 0; ijk < trip_aba[h]; ijk++) {

            int i = bas_aba_sym[h][ijk][0];
            int j = bas_aba_sym[h][ijk][1];
            int k = bas_aba_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aba[h]; lmn++) {

                int l = bas_aba_sym[h][lmn][0];
                int m = bas_aba_sym[h][lmn][1];
                int n = bas_aba_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int id = (ijk+trip_aab[h])*(trip_aab[h]+trip_aba[h])+(lmn+trip_aab[h]);

                double dum = u[offset + id];

                A[f3aaaoff[h] + id] -= dum; // - F3(ijk,lmn)

                if ( j != n && m != k ) {
                    int jni = ibas_aab_sym[hijn][j][n][i];
                    int mkl = ibas_aab_sym[hijn][m][k][l];
                    int sg = 1;
                    if ( n < j ) sg = -sg;
                    if ( k < m ) sg = -sg;
                    A[d3bbaoff[hijn] + jni*trip_aab[hijn]+mkl] += sg * dum; // - D3(jni,mkl)
                }

                if ( j == m && i == l ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    A[d1boff[h2] + nn*amopi_[h2]+kk] += dum; // + D1(n,k) djm dil
                }

                if ( i == l ) {
                    if ( n != j && k != m ) {
                        int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                        int nj = ibas_aa_sym[hnj][n][j];
                        int km = ibas_aa_sym[hnj][k][m];

                        int s = 1;
                        if ( n > j ) s = -s;
                        if ( k > m ) s = -s;

                        A[d2bboff[hnj] + nj*gems_aa[hnj]+km] -= s * dum; // - D2(nj,km) dil
                    }
                }
                if ( j == m ) {
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][i][n];
                    int kl = ibas_ab_sym[hni][l][k];
                    A[d2aboff[hni] + ni*gems_ab[hni]+kl] -= dum; // -D2(ni,kl) djm
                }

            }
        }
        offset += (trip_aba[h]+trip_aab[h])*(trip_aba[h]+trip_aab[h]);
    }

    // big block 2: F3bbb + F3baa
    for (int h = 0; h < nirrep_; h++) {

        // F3bbb/bbb
        for (int ijk = 0; ijk < trip_aab[h]; ijk++) {

            int i = bas_aab_sym[h][ijk][0];
            int j = bas_aab_sym[h][ijk][1];
            int k = bas_aab_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aab[h]; lmn++) {

                int l = bas_aab_sym[h][lmn][0];
                int m = bas_aab_sym[h][lmn][1];
                int n = bas_aab_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int id = ijk*(trip_aab[h]+trip_aba[h])+lmn;
                //int id = ijk*trip_aab[h]+lmn;

                double dum = u[offset + id]; // - F3(ijk,lmn)

                A[f3bbboff[h] + id] -= dum; // - F3(ijk,lmn)

                if ( i != n && j != n && l != k && m != k ) {

                    int ijn = ibas_aaa_sym[hijn][i][j][n];
                    int lmk = ibas_aaa_sym[hijn][l][m][k];
                    
                    // sign ... i < j, l < m guaranteed
                    //int sg = 1;
                    //if ( n < i && n < j ) {
                    //   // nothing
                    //}else if ( n > i && n > j ) {
                    //   // nothing
                    //}else if ( n > i && n < j ) {
                    //    sg = -sg;
                    //}else if ( n < i && n > j ) {
                    //    sg = -sg; // never will happen
                    //} 
                    int sg = 1;
                    if ( i < n && n < j ) {
                        sg = -sg;
                    }
                    if ( l < k && k < m ) {
                        sg = -sg;
                    }

                    A[d3bbboff[hijn] + ijn*trip_aaa[hijn]+lmk] += sg * dum; // - D3(ijn,lmk)

                }

                if ( j == m && i == l ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    A[d1boff[h2] + nn*amopi_[h2]+kk] += dum; // + D1(n,k) djm dil
                }
                if ( j == l && i == m ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    A[d1boff[h2] + nn*amopi_[h2]+kk] -= dum; // - D1(n,k) djl dim
                }

                if ( i == l ) {
                    if ( n != j && k != m ) {
                        int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                        int nj = ibas_aa_sym[hnj][n][j];
                        int km = ibas_aa_sym[hnj][k][m];
                
                        int s = 1;
                        if ( n > j ) s = -s;
                        if ( k > m ) s = -s;
                    
                        A[d2bboff[hnj] + nj*gems_aa[hnj]+km] -= s * dum; // - D2(nj,km) dil
                    }
                }
                if ( j == l ) {
                    if ( n != i && k != m ) {
                        int hni = SymmetryPair(symmetry[n],symmetry[i]);
                        int ni = ibas_aa_sym[hni][n][i];
                        int km = ibas_aa_sym[hni][k][m];
                
                        int s = 1;
                        if ( n > i ) s = -s;
                        if ( k > m ) s = -s;
                        
                        A[d2bboff[hni] + ni*gems_aa[hni]+km] += s * dum; // D2(ni,km) djl
                    }   
                }
                if ( i == m ) { 
                    if ( n != j && k != l ) {
                        int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                        int nj = ibas_aa_sym[hnj][n][j];
                        int kl = ibas_aa_sym[hnj][k][l];
                    
                        int s = 1;
                        if ( n > j ) s = -s;
                        if ( k > l ) s = -s;
                        
                        A[d2bboff[hnj] + nj*gems_aa[hnj]+kl] += s * dum; // D2(nj,kl) dim
                    }   
                }
                if ( j == m ) { 
                    if ( n != i && k != l ) {
                        int hni = SymmetryPair(symmetry[n],symmetry[i]);
                        int ni = ibas_aa_sym[hni][n][i];
                        int kl = ibas_aa_sym[hni][k][l];
                    
                        int s = 1;
                        if ( n > i ) s = -s;
                        if ( k > l ) s = -s;
                        
                        A[d2bboff[hni] + ni*gems_aa[hni]+kl] -= s * dum; // -D2(ni,kl) djm
                    }   
                }
            }
        }
        // F3bbb/baa
        for (int ijk = 0; ijk < trip_aab[h]; ijk++) {

            int i = bas_aab_sym[h][ijk][0];
            int j = bas_aab_sym[h][ijk][1];
            int k = bas_aab_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aba[h]; lmn++) {

                int l = bas_aba_sym[h][lmn][0];
                int m = bas_aba_sym[h][lmn][1];
                int n = bas_aba_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int id = ijk*(trip_aab[h]+trip_aba[h])+(lmn+trip_aab[h]);

                double dum = u[offset + id]; // - F3(ijk,lmn)

                A[f3bbboff[h] + id] -= dum; // - F3(ijk,lmn)

                if ( i != j && l != k ) {
                    int ijn = ibas_aab_sym[hijn][i][j][n];
                    int lkm = ibas_aab_sym[hijn][l][k][m];
                    int sg = 1;
                    if ( i > j ) sg = -sg;
                    if ( l > k ) sg = -sg;
                    A[d3bbaoff[hijn] + ijn*trip_aab[hijn]+lkm] -= sg * dum; // + D3(ijn,lkm)
                } 

                if ( i == l ) {
                    int hjn = SymmetryPair(symmetry[j],symmetry[n]);
                    int jn = ibas_ab_sym[hjn][n][j];
                    int km = ibas_ab_sym[hjn][m][k];
                    A[d2aboff[hjn]+jn*gems_ab[hjn]+km] += dum; // D2(jn,km) dil
                }
                if ( j == l ) {
                    int hin = SymmetryPair(symmetry[i],symmetry[n]);
                    int in = ibas_ab_sym[hin][n][i];
                    int km = ibas_ab_sym[hin][m][k];
                    A[d2aboff[hin]+in*gems_ab[hin]+km] -= dum; // -D2(in,km) djl
                }

            }
        }

        // F3baa/bbb
        for (int ijk = 0; ijk < trip_aba[h]; ijk++) {

            int i = bas_aba_sym[h][ijk][0];
            int j = bas_aba_sym[h][ijk][1];
            int k = bas_aba_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aab[h]; lmn++) {

                int l = bas_aab_sym[h][lmn][0];
                int m = bas_aab_sym[h][lmn][1];
                int n = bas_aab_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int id = (ijk+trip_aab[h])*(trip_aab[h]+trip_aba[h])+lmn;

                double dum = u[offset + id]; 

                A[f3bbboff[h] + id] -= dum; // - F3(ijk,lmn)

                if ( i != n && l != m ) {
                    int inj = ibas_aab_sym[hijn][i][n][j];
                    int lmk = ibas_aab_sym[hijn][l][m][k];
                    int sg = 1;
                    if ( i > n ) sg = -sg;
                    if ( l > m ) sg = -sg;
                    A[d3bbaoff[hijn] + inj*trip_aab[hijn]+lmk] -= sg * dum; // + D3(inj,lmk)
                }

                if ( i == l ) {
                    int hjn = SymmetryPair(symmetry[j],symmetry[n]);
                    int jn = ibas_ab_sym[hjn][j][n];
                    int km = ibas_ab_sym[hjn][k][m];
                    A[d2aboff[hjn]+jn*gems_ab[hjn]+km] += dum; // D2(jn,km) dil
                }
                if ( i == m ) {
                    int hnj = SymmetryPair(symmetry[j],symmetry[n]);
                    int nj = ibas_ab_sym[hnj][j][n];
                    int lk = ibas_ab_sym[hnj][k][l];
                    A[d2aboff[hnj]+nj*gems_ab[hnj]+lk] -= dum; // -D2(in,km) djl
                }
            }
        }
        // F3baa/baa
        for (int ijk = 0; ijk < trip_aba[h]; ijk++) {

            int i = bas_aba_sym[h][ijk][0];
            int j = bas_aba_sym[h][ijk][1];
            int k = bas_aba_sym[h][ijk][2];

            int hi = symmetry[i];
            int hj = symmetry[j];
            int hk = symmetry[k];

            for (int lmn = 0; lmn < trip_aba[h]; lmn++) {

                int l = bas_aba_sym[h][lmn][0];
                int m = bas_aba_sym[h][lmn][1];
                int n = bas_aba_sym[h][lmn][2];

                int hl = symmetry[l];
                int hm = symmetry[m];
                int hn = symmetry[n];

                int hijn = (hi ^ hj) ^ hn;
                int hlmk = (hl ^ hm) ^ hk;

                int id = (ijk+trip_aab[h])*(trip_aab[h]+trip_aba[h])+(lmn+trip_aab[h]);

                double dum = u[offset + id]; // - F3(ijk,lmn)

                A[f3bbboff[h] + id] -= dum; // - F3(ijk,lmn)

                if ( j != n && m != k ) {
                    int jni = ibas_aab_sym[hijn][j][n][i];
                    int mkl = ibas_aab_sym[hijn][m][k][l];
                    int sg = 1;
                    if ( j > n ) sg = -sg;
                    if ( m > k ) sg = -sg;
                    A[d3aaboff[hijn] + jni*trip_aab[hijn]+mkl] += sg * dum; // - D3(jni,mkl)
                }

                if ( j == m && i == l ) {
                    int h2 = symmetry[k];
                    int kk = k - pitzer_offset[h2];
                    int nn = n - pitzer_offset[h2];
                    A[d1aoff[h2] + nn*amopi_[h2]+kk] += dum; // + D1(n,k) djm dil
                }

                if ( i == l ) {
                    if ( n != j && k != m ) {
                        int hnj = SymmetryPair(symmetry[n],symmetry[j]);
                        int nj = ibas_aa_sym[hnj][n][j];
                        int km = ibas_aa_sym[hnj][k][m];

                        int s = 1;
                        if ( n > j ) s = -s;
                        if ( k > m ) s = -s;

                        A[d2aaoff[hnj] + nj*gems_aa[hnj]+km] -= s * dum; // - D2(nj,km) dil
                    }
                }
                if ( j == m ) {
                    int hni = SymmetryPair(symmetry[n],symmetry[i]);
                    int ni = ibas_ab_sym[hni][n][i];
                    int kl = ibas_ab_sym[hni][k][l];
                    A[d2aboff[hni] + ni*gems_ab[hni]+kl] -= dum; // -D2(ni,kl) djm
                }

            }
        }
        offset += (trip_aba[h]+trip_aab[h])*(trip_aba[h]+trip_aab[h]);
    }

}

} // end namespaces
