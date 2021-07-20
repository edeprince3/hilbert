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

// E3 portion of A.u
void v2RDMSolver::E3_constraints_Au(double * A,double * u){

    // E3aab
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

                double dum = -u[e3aaboff[h] + ijk*trip_aab[h]+lmn]; // - E3(ijk,lmn)
                dum       -=  u[d3aaboff[hijn] + ijn*trip_aab[hijn]+lmk]; // - D3(ijn,lmk)

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    dum += u[d2aaoff[hij] + ij*gems_aa[hij]+lm]; // + D2(ij,lm) dkn
                }

                A[offset + ijk*trip_aab[h]+lmn] = dum;

            }
        }
        offset += trip_aab[h]*trip_aab[h];
    }

    // E3bba
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

                double dum = -u[e3bbaoff[h] + ijk*trip_aab[h]+lmn]; // - E3(ijk,lmn)
                dum       -=  u[d3bbaoff[hijn] + ijn*trip_aab[hijn]+lmk]; // - D3(ijn,lmk)

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    dum += u[d2bboff[hij] + ij*gems_aa[hij]+lm]; // + D2(ij,lm) dkn
                }

                A[offset + ijk*trip_aab[h]+lmn] = dum;

            }
        }
        offset += trip_aab[h]*trip_aab[h];
    }

    // big block 1: E3aaa + E3abb
    for (int h = 0; h < nirrep_; h++) {

        // E3aaa/aaa
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

                double dum = -u[e3aaaoff[h] + id]; // - E3(ijk,lmn)

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

                    dum -=  sg * u[d3aaaoff[hijn] + ijn*trip_aaa[hijn]+lmk]; // - D3(ijn,lmk)

                }

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    dum += u[d2aaoff[hij] + ij*gems_aa[hij]+lm]; // + D2(ij,lm) dkn
                }

                A[offset + id] = dum;

            }
        }
        // E3aaa/abb
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

                double dum = -u[e3aaaoff[h] + id]; // - E3(ijk,lmn)

                if ( i != j && l != k ) {

                    int ijn = ibas_aab_sym[hijn][i][j][n];
                    int lkm = ibas_aab_sym[hijn][l][k][m];

                    int sg = 1;
                    if ( l > k ) {
                        sg = -sg;
                    }

                    dum +=  sg * u[d3aaboff[hijn] + ijn*trip_aab[hijn]+lkm]; // + D3(ijn,lkm)
                }

                A[offset + id] = dum;

            }
        }

        // E3abb/aaa
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

                double dum = -u[e3aaaoff[h] + id]; // - E3(ijk,lmn)

                if ( i != n && l != m ) {
                    
                    int inj = ibas_aab_sym[hijn][i][n][j];
                    int lmk = ibas_aab_sym[hijn][l][m][k];
                    
                    int sg = 1;
                    if ( i > n ) {
                        sg = -sg;
                    }

                    dum +=  sg * u[d3aaboff[hijn] + inj*trip_aab[hijn]+lmk]; // + D3(inj,lmk)
                }  

                A[offset + id] = dum;

            }
        }

        // E3abb/abb
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

                double dum = -u[e3aaaoff[h] + id]; // - E3(ijk,lmn)

                if ( j != n && m != k ) {
                    int jni = ibas_aab_sym[hijn][j][n][i];
                    int mkl = ibas_aab_sym[hijn][m][k][l];
                    int sg = 1;
                    if ( j > n ) sg = -sg;
                    if ( m > k ) sg = -sg;
                    dum -= sg * u[d3bbaoff[hijn] + jni*trip_aab[hijn]+mkl]; // - D3(jni,mkl)
                }

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_ab_sym[hij][i][j];
                    int lm = ibas_ab_sym[hij][l][m];
                    dum += u[d2aboff[hij] + ij*gems_ab[hij]+lm]; // + D2(ij,lm) dkn
                }

                A[offset + id] = dum;

            }
        }
        offset += (trip_aba[h]+trip_aab[h])*(trip_aba[h]+trip_aab[h]);
    }

    // big block 2: E3bbb + E3baa
    for (int h = 0; h < nirrep_; h++) {

        // E3bbb/bbb
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

                double dum = -u[e3bbboff[h] + id]; // - E3(ijk,lmn)

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

                    dum -=  sg * u[d3bbboff[hijn] + ijn*trip_aaa[hijn]+lmk]; // - D3(ijn,lmk)

                }

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    dum += u[d2bboff[hij] + ij*gems_aa[hij]+lm]; // + D2(ij,lm) dkn
                }

                A[offset + id] = dum;

            }
        }
        // E3bbb/baa
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

                double dum = -u[e3bbboff[h] + id]; // - E3(ijk,lmn)

                if ( i != j && l != k ) {
                    int ijn = ibas_aab_sym[hijn][i][j][n];
                    int lkm = ibas_aab_sym[hijn][l][k][m];
                    int sg = 1;
                    if ( i > j ) sg = -sg;
                    if ( l > k ) sg = -sg;
                    dum += sg * u[d3bbaoff[hijn] + ijn*trip_aab[hijn]+lkm]; // + D3(ijn,lkm)
                } 

                A[offset + id] = dum;

            }
        }

        // E3baa/bbb
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

                double dum = -u[e3bbboff[h] + id]; // - E3(ijk,lmn)

                if ( i != n && l != m ) {
                    int inj = ibas_aab_sym[hijn][i][n][j];
                    int lmk = ibas_aab_sym[hijn][l][m][k];
                    int sg = 1;
                    if ( i > n ) sg = -sg;
                    if ( l > m ) sg = -sg;
                    dum += sg * u[d3bbaoff[hijn] + inj*trip_aab[hijn]+lmk]; // + D3(inj,lmk)
                }

                A[offset + id] = dum;

            }
        }
        // E3baa/baa
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

                double dum = -u[e3bbboff[h] + id]; // - E3(ijk,lmn)

                if ( j != n && m != k ) {
                    int jni = ibas_aab_sym[hijn][j][n][i];
                    int mkl = ibas_aab_sym[hijn][m][k][l];
                    int sg = 1;
                    if ( j > n ) sg = -sg;
                    if ( m > k ) sg = -sg;
                    dum -= sg * u[d3aaboff[hijn] + jni*trip_aab[hijn]+mkl]; // - D3(jni,mkl)
                }

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_ab_sym[hij][j][i];
                    int lm = ibas_ab_sym[hij][m][l];
                    dum += u[d2aboff[hij] + ij*gems_ab[hij]+lm]; // + D2(ij,lm) dkn
                }

                A[offset + id] = dum;

            }
        }
        offset += (trip_aba[h]+trip_aab[h])*(trip_aba[h]+trip_aab[h]);
    }

}

// E3 portion of AT.u
void v2RDMSolver::E3_constraints_ATu(double * A,double * u){

    // E3aab
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

                A[e3aaboff[h] + ijk*trip_aab[h]+lmn] -= dum; // - E3(ijk,lmn)
                A[d3aaboff[hijn] + ijn*trip_aab[hijn]+lmk] -= dum; // - D3(ijn,lmk)

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    A[d2aaoff[hij] + ij*gems_aa[hij]+lm] += dum; // + D2(ij,lm) dkn
                }

            }
        }
        offset += trip_aab[h]*trip_aab[h];
    }

    // E3bba
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

                A[e3bbaoff[h] + ijk*trip_aab[h]+lmn] -= dum; // - E3(ijk,lmn)
                A[d3bbaoff[hijn] + ijn*trip_aab[hijn]+lmk] -= dum; // - D3(ijn,lmk)

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    A[d2bboff[hij] + ij*gems_aa[hij]+lm] += dum; // + D2(ij,lm) dkn
                }

            }
        }
        offset += trip_aab[h]*trip_aab[h];
    }

    // big block 1: E3aaa + E3abb
    for (int h = 0; h < nirrep_; h++) {

        // E3aaa/aaa
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

                A[e3aaaoff[h] + id] -= dum; // - E3(ijk,lmn)

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

                    A[d3aaaoff[hijn] + ijn*trip_aaa[hijn]+lmk] -= sg * dum; // - D3(ijn,lmk)

                }

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    A[d2aaoff[hij] + ij*gems_aa[hij]+lm] += dum; // + D2(ij,lm) dkn
                }

            }
        }
        // E3aaa/abb
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

                double dum = u[offset + id]; // - E3(ijk,lmn)

                A[e3aaaoff[h] + id] -= dum; // - E3(ijk,lmn)

                if ( i != j && l != k ) {

                    int ijn = ibas_aab_sym[hijn][i][j][n];
                    int lkm = ibas_aab_sym[hijn][l][k][m];

                    int sg = 1;
                    if ( k < l ) {
                        sg = -sg;
                    }

                    A[d3aaboff[hijn] + ijn*trip_aab[hijn]+lkm] += sg * dum; // + D3(ijn,lkm)
                }

            }
        }

        // E3abb/aaa
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

                A[e3aaaoff[h] + id] -= dum; // - E3(ijk,lmn)

                if ( i != n && l != m ) {
                    
                    int inj = ibas_aab_sym[hijn][i][n][j];
                    int lmk = ibas_aab_sym[hijn][l][m][k];
                    
                    int sg = 1;
                    if ( n < i ) {
                        sg = -sg;
                    }

                    A[d3aaboff[hijn] + inj*trip_aab[hijn]+lmk] += sg * dum; // + D3(inj,lmk)
                }  

            }
        }

        // E3abb/abb
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

                A[e3aaaoff[h] + id] -= dum; // - E3(ijk,lmn)

                if ( j != n && m != k ) {
                    int jni = ibas_aab_sym[hijn][j][n][i];
                    int mkl = ibas_aab_sym[hijn][m][k][l];
                    int sg = 1;
                    if ( n < j ) sg = -sg;
                    if ( k < m ) sg = -sg;
                    A[d3bbaoff[hijn] + jni*trip_aab[hijn]+mkl] -= sg * dum; // - D3(jni,mkl)
                }

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_ab_sym[hij][i][j];
                    int lm = ibas_ab_sym[hij][l][m];
                    A[d2aboff[hij] + ij*gems_ab[hij]+lm] += dum; // + D2(ij,lm) dkn
                }

            }
        }
        offset += (trip_aba[h]+trip_aab[h])*(trip_aba[h]+trip_aab[h]);
    }

    // big block 2: E3bbb + E3baa
    for (int h = 0; h < nirrep_; h++) {

        // E3bbb/bbb
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

                double dum = u[offset + id]; // - E3(ijk,lmn)

                A[e3bbboff[h] + id] -= dum; // - E3(ijk,lmn)

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

                    A[d3bbboff[hijn] + ijn*trip_aaa[hijn]+lmk] -= sg * dum; // - D3(ijn,lmk)

                }

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_aa_sym[hij][i][j];
                    int lm = ibas_aa_sym[hij][l][m];
                    A[d2bboff[hij] + ij*gems_aa[hij]+lm] += dum; // + D2(ij,lm) dkn
                }

            }
        }
        // E3bbb/baa
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

                double dum = u[offset + id]; // - E3(ijk,lmn)

                A[e3bbboff[h] + id] -= dum; // - E3(ijk,lmn)

                if ( i != j && l != k ) {
                    int ijn = ibas_aab_sym[hijn][i][j][n];
                    int lkm = ibas_aab_sym[hijn][l][k][m];
                    int sg = 1;
                    if ( i > j ) sg = -sg;
                    if ( l > k ) sg = -sg;
                    A[d3bbaoff[hijn] + ijn*trip_aab[hijn]+lkm] += sg * dum; // + D3(ijn,lkm)
                } 

            }
        }

        // E3baa/bbb
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

                A[e3bbboff[h] + id] -= dum; // - E3(ijk,lmn)

                if ( i != n && l != m ) {
                    int inj = ibas_aab_sym[hijn][i][n][j];
                    int lmk = ibas_aab_sym[hijn][l][m][k];
                    int sg = 1;
                    if ( i > n ) sg = -sg;
                    if ( l > m ) sg = -sg;
                    A[d3bbaoff[hijn] + inj*trip_aab[hijn]+lmk] += sg * dum; // + D3(inj,lmk)
                }

            }
        }
        // E3baa/baa
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

                double dum = u[offset + id]; // - E3(ijk,lmn)

                A[e3bbboff[h] + id] -= dum; // - E3(ijk,lmn)

                if ( j != n && m != k ) {
                    int jni = ibas_aab_sym[hijn][j][n][i];
                    int mkl = ibas_aab_sym[hijn][m][k][l];
                    int sg = 1;
                    if ( j > n ) sg = -sg;
                    if ( m > k ) sg = -sg;
                    A[d3aaboff[hijn] + jni*trip_aab[hijn]+mkl] -= sg * dum; // - D3(jni,mkl)
                }

                if ( k == n ) {
                    int hij = SymmetryPair(symmetry[i],symmetry[j]);
                    int ij = ibas_ab_sym[hij][j][i];
                    int lm = ibas_ab_sym[hij][m][l];
                    A[d2aboff[hij] + ij*gems_ab[hij]+lm] += dum; // + D2(ij,lm) dkn
                }

            }
        }
        offset += (trip_aba[h]+trip_aab[h])*(trip_aba[h]+trip_aab[h]);
    }

}

} // end namespaces
