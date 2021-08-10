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

#include <psi4/liboptions/liboptions.h>
#include <psi4/libqt/qt.h>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libtrans/mospace.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/matrix.h>

#include "v2rdm_solver.h"

#include<misc/omp.h>

using namespace psi;

namespace hilbert{

// D4 portion of A.u 
void v2RDMSolver::D4_constraints_Au(double* A,double* u){

    int na = nalpha_ - nrstc_ - nfrzc_;
    int nb = nbeta_ - nrstc_ - nfrzc_;
    // D4aaaa -> D3aaa
    if ( na > 3 ) {
        for (int h = 0; h < nirrep_; h++) {
            for (int ijk = 0; ijk < trip_aaa[h]; ijk++) {
                int i = bas_aaa_sym[h][ijk][0];
                int j = bas_aaa_sym[h][ijk][1];
                int k = bas_aaa_sym[h][ijk][2];
                for (int lmn = 0; lmn < trip_aaa[h]; lmn++) {
                    int l = bas_aaa_sym[h][lmn][0];
                    int m = bas_aaa_sym[h][lmn][1];
                    int n = bas_aaa_sym[h][lmn][2];
                    double dum = (na - 3.0) * u[d3aaaoff[h] + ijk*trip_aaa[h] + lmn];
                    for ( int p = 0; p < amo_; p++) {
                        if ( i == p) continue;
                        if ( j == p) continue;
                        if ( k == p) continue;
                        if ( l == p) continue;
                        if ( m == p) continue;
                        if ( n == p) continue;
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijkp = ibas_aaaa_sym[h2][i][j][k][p];
                        int lmnp = ibas_aaaa_sym[h2][l][m][n][p];
                        int s = 1;
                        if ( p < i ) s = -s;
                        if ( p < j ) s = -s;
                        if ( p < k ) s = -s;
                        if ( p < l ) s = -s;
                        if ( p < m ) s = -s;
                        if ( p < n ) s = -s;
                        dum -= s * u[d4aaaaoff[h2] + ijkp*quartet_aaaa[h2]+lmnp];
                    }
                    A[offset + ijk*trip_aaa[h]+lmn] = dum;
                }
            }
            offset += trip_aaa[h] * trip_aaa[h];
        }
    }
    // D4aaab
    if ( na > 2 && nb > 0 ) {
        // D4aaab -> D3aaa
        for (int h = 0; h < nirrep_; h++) {
            for (int ijk = 0; ijk < trip_aaa[h]; ijk++) {
                int i = bas_aaa_sym[h][ijk][0];
                int j = bas_aaa_sym[h][ijk][1];
                int k = bas_aaa_sym[h][ijk][2];
                for (int lmn = 0; lmn < trip_aaa[h]; lmn++) {
                    int l = bas_aaa_sym[h][lmn][0];
                    int m = bas_aaa_sym[h][lmn][1];
                    int n = bas_aaa_sym[h][lmn][2];
                    double dum = nb * u[d3aaaoff[h] + ijk*trip_aaa[h] + lmn];
                    for ( int p = 0; p < amo_; p++) {
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijkp = ibas_aaab_sym[h2][i][j][k][p];
                        int lmnp = ibas_aaab_sym[h2][l][m][n][p];
                        dum -= u[d4aaaboff[h2] + ijkp*quartet_aaab[h2]+lmnp];
                    }
                    A[offset + ijk*trip_aaa[h]+lmn] = dum;
                }
            }
            offset += trip_aaa[h] * trip_aaa[h];
        }

        // D4aaab -> D3aab
        for (int h = 0; h < nirrep_; h++) {
            for (int ijk = 0; ijk < trip_aab[h]; ijk++) {
                int i = bas_aab_sym[h][ijk][0];
                int j = bas_aab_sym[h][ijk][1];
                int k = bas_aab_sym[h][ijk][2];
                for (int lmn = 0; lmn < trip_aab[h]; lmn++) {
                    int l = bas_aab_sym[h][lmn][0];
                    int m = bas_aab_sym[h][lmn][1];
                    int n = bas_aab_sym[h][lmn][2];
                    double dum = (na - 2.0) * u[d3aaboff[h] + ijk*trip_aab[h] + lmn];
                    for ( int p = 0; p < amo_; p++) {
                        if ( i == p) continue;
                        if ( j == p) continue;
                        if ( l == p) continue;
                        if ( m == p) continue;
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijkp = ibas_aaab_sym[h2][i][j][p][k];
                        int lmnp = ibas_aaab_sym[h2][l][m][p][n];
                        int s = 1;
                        if ( p < i ) s = -s;
                        if ( p < j ) s = -s;
                        if ( p < l ) s = -s;
                        if ( p < m ) s = -s;
                        dum -= s * u[d4aaaboff[h2] + ijkp*quartet_aaab[h2]+lmnp];
                    }
                    A[offset + ijk*trip_aab[h]+lmn] = dum;
                }
            }
            offset += trip_aab[h] * trip_aab[h];
        }
    }
    // D4aabb
    if ( na > 1 && nb > 1 ) {
        // D4aabb -> D3aab
        for (int h = 0; h < nirrep_; h++) {
            for (int ijk = 0; ijk < trip_aab[h]; ijk++) {
                int i = bas_aab_sym[h][ijk][0];
                int j = bas_aab_sym[h][ijk][1];
                int k = bas_aab_sym[h][ijk][2];
                for (int lmn = 0; lmn < trip_aab[h]; lmn++) {
                    int l = bas_aab_sym[h][lmn][0];
                    int m = bas_aab_sym[h][lmn][1];
                    int n = bas_aab_sym[h][lmn][2];
                    double dum = (nb-1.0) * u[d3aaboff[h] + ijk*trip_aab[h] + lmn];
                    for ( int p = 0; p < amo_; p++) {
                        int h2 = SymmetryPair(h,symmetry[p]);
                        if ( k == p ) continue;
                        if ( n == p ) continue;
                        int ijkp = ibas_aabb_sym[h2][i][j][k][p];
                        int lmnp = ibas_aabb_sym[h2][l][m][n][p];
                        int s = 1;
                        if ( p < k ) s = -s;
                        if ( p < n ) s = -s;
                        dum -= s * u[d4aabboff[h2] + ijkp*quartet_aabb[h2]+lmnp];
                    }
                    A[offset + ijk*trip_aab[h]+lmn] = dum;
                }
            }
            offset += trip_aab[h] * trip_aab[h];
        }
        // D4aabb -> D3bba
        for (int h = 0; h < nirrep_; h++) {
            for (int ijk = 0; ijk < trip_aab[h]; ijk++) {
                int i = bas_aab_sym[h][ijk][0];
                int j = bas_aab_sym[h][ijk][1];
                int k = bas_aab_sym[h][ijk][2];
                for (int lmn = 0; lmn < trip_aab[h]; lmn++) {
                    int l = bas_aab_sym[h][lmn][0];
                    int m = bas_aab_sym[h][lmn][1];
                    int n = bas_aab_sym[h][lmn][2];
                    double dum = (na-1.0) * u[d3bbaoff[h] + ijk*trip_aab[h] + lmn];
                    for ( int p = 0; p < amo_; p++) {
                        int h2 = SymmetryPair(h,symmetry[p]);
                        if ( k == p ) continue;
                        if ( n == p ) continue;
                        int ijkp = ibas_aabb_sym[h2][k][p][i][j];
                        int lmnp = ibas_aabb_sym[h2][n][p][l][m];
                        int s = 1;
                        if ( p < k ) s = -s;
                        if ( p < n ) s = -s;
                        dum -= s * u[d4aabboff[h2] + ijkp*quartet_aabb[h2]+lmnp];
                    }
                    A[offset + ijk*trip_aab[h]+lmn] = dum;
                }
            }
            offset += trip_aab[h] * trip_aab[h];
        }
    }


    // D4bbba
    if ( na > 0 && nb > 2 ) {
        // D4bbba -> D3bbb
        for (int h = 0; h < nirrep_; h++) {
            for (int ijk = 0; ijk < trip_aaa[h]; ijk++) {
                int i = bas_aaa_sym[h][ijk][0];
                int j = bas_aaa_sym[h][ijk][1];
                int k = bas_aaa_sym[h][ijk][2];
                for (int lmn = 0; lmn < trip_aaa[h]; lmn++) {
                    int l = bas_aaa_sym[h][lmn][0];
                    int m = bas_aaa_sym[h][lmn][1];
                    int n = bas_aaa_sym[h][lmn][2];
                    double dum = na * u[d3bbboff[h] + ijk*trip_aaa[h] + lmn];
                    for ( int p = 0; p < amo_; p++) {
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijkp = ibas_aaab_sym[h2][i][j][k][p];
                        int lmnp = ibas_aaab_sym[h2][l][m][n][p];
                        dum -= u[d4bbbaoff[h2] + ijkp*quartet_aaab[h2]+lmnp];
                    }
                    A[offset + ijk*trip_aaa[h]+lmn] = dum;
                }
            }
            offset += trip_aaa[h] * trip_aaa[h];
        }

        // D4bbba -> D3bba
        for (int h = 0; h < nirrep_; h++) {
            for (int ijk = 0; ijk < trip_aab[h]; ijk++) {
                int i = bas_aab_sym[h][ijk][0];
                int j = bas_aab_sym[h][ijk][1];
                int k = bas_aab_sym[h][ijk][2];
                for (int lmn = 0; lmn < trip_aab[h]; lmn++) {
                    int l = bas_aab_sym[h][lmn][0];
                    int m = bas_aab_sym[h][lmn][1];
                    int n = bas_aab_sym[h][lmn][2];
                    double dum = (nb - 2.0) * u[d3bbaoff[h] + ijk*trip_aab[h] + lmn];
                    for ( int p = 0; p < amo_; p++) {
                        if ( i == p) continue;
                        if ( j == p) continue;
                        if ( l == p) continue;
                        if ( m == p) continue;
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijkp = ibas_aaab_sym[h2][i][j][p][k];
                        int lmnp = ibas_aaab_sym[h2][l][m][p][n];
                        int s = 1;
                        if ( p < i ) s = -s;
                        if ( p < j ) s = -s;
                        if ( p < l ) s = -s;
                        if ( p < m ) s = -s;
                        dum -= s * u[d4bbbaoff[h2] + ijkp*quartet_aaab[h2]+lmnp];
                    }
                    A[offset + ijk*trip_aab[h]+lmn] = dum;
                }
            }
            offset += trip_aab[h] * trip_aab[h];
        }
    }

    // D4bbbb -> D3bbb
    if ( nb > 3 ) {
        for (int h = 0; h < nirrep_; h++) {
            for (int ijk = 0; ijk < trip_aaa[h]; ijk++) {
                int i = bas_aaa_sym[h][ijk][0];
                int j = bas_aaa_sym[h][ijk][1];
                int k = bas_aaa_sym[h][ijk][2];
                for (int lmn = 0; lmn < trip_aaa[h]; lmn++) {
                    int l = bas_aaa_sym[h][lmn][0];
                    int m = bas_aaa_sym[h][lmn][1];
                    int n = bas_aaa_sym[h][lmn][2];
                    double dum = (nb - 3.0) * u[d3bbboff[h] + ijk*trip_aaa[h] + lmn];
                    for ( int p = 0; p < amo_; p++) {
                        if ( i == p) continue;
                        if ( j == p) continue;
                        if ( k == p) continue;
                        if ( l == p) continue;
                        if ( m == p) continue;
                        if ( n == p) continue;
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijkp = ibas_aaaa_sym[h2][i][j][k][p];
                        int lmnp = ibas_aaaa_sym[h2][l][m][n][p];
                        int s = 1;
                        if ( p < i ) s = -s;
                        if ( p < j ) s = -s;
                        if ( p < k ) s = -s;
                        if ( p < l ) s = -s;
                        if ( p < m ) s = -s;
                        if ( p < n ) s = -s;
                        dum -= s * u[d4bbbboff[h2] + ijkp*quartet_aaaa[h2]+lmnp];
                    }
                    A[offset + ijk*trip_aaa[h]+lmn] = dum;
                }
            }
            offset += trip_aaa[h] * trip_aaa[h];
        }
    }

}

// D4 portion of A^T.y 
void v2RDMSolver::D4_constraints_ATu(double* A,double* u){

    int na = nalpha_ - nrstc_ - nfrzc_;
    int nb = nbeta_ - nrstc_ - nfrzc_;

    // D4aaaa -> D3aaa
    if ( na > 3 ) {
        for (int h = 0; h < nirrep_; h++) {
            for (int ijk = 0; ijk < trip_aaa[h]; ijk++) {
                int i = bas_aaa_sym[h][ijk][0];
                int j = bas_aaa_sym[h][ijk][1];
                int k = bas_aaa_sym[h][ijk][2];
                for (int lmn = 0; lmn < trip_aaa[h]; lmn++) {
                    int l = bas_aaa_sym[h][lmn][0];
                    int m = bas_aaa_sym[h][lmn][1];
                    int n = bas_aaa_sym[h][lmn][2];
                    double dum = u[offset + ijk*trip_aaa[h]+lmn];
                    A[d3aaaoff[h] + ijk*trip_aaa[h] + lmn] += (na - 3.0) * dum;
                    for ( int p = 0; p < amo_; p++) {
                        if ( i == p) continue;
                        if ( j == p) continue;
                        if ( k == p) continue;
                        if ( l == p) continue;
                        if ( m == p) continue;
                        if ( n == p) continue;
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijkp = ibas_aaaa_sym[h2][i][j][k][p];
                        int lmnp = ibas_aaaa_sym[h2][l][m][n][p];
                        int s = 1;
                        if ( p < i ) s = -s;
                        if ( p < j ) s = -s;
                        if ( p < k ) s = -s;
                        if ( p < l ) s = -s;
                        if ( p < m ) s = -s;
                        if ( p < n ) s = -s;
                        A[d4aaaaoff[h2] + ijkp*quartet_aaaa[h2]+lmnp] -= s * dum;
                    }
                }
            }
            offset += trip_aaa[h] * trip_aaa[h];
        }
    }

    // D4aaab
    if ( na > 2 && nb > 0 ) {
        // D4aaab -> D3aaa
        for (int h = 0; h < nirrep_; h++) {
            for (int ijk = 0; ijk < trip_aaa[h]; ijk++) {
                int i = bas_aaa_sym[h][ijk][0];
                int j = bas_aaa_sym[h][ijk][1];
                int k = bas_aaa_sym[h][ijk][2];
                for (int lmn = 0; lmn < trip_aaa[h]; lmn++) {
                    int l = bas_aaa_sym[h][lmn][0];
                    int m = bas_aaa_sym[h][lmn][1];
                    int n = bas_aaa_sym[h][lmn][2];
                    double dum = u[offset + ijk*trip_aaa[h]+lmn];

                    A[d3aaaoff[h] + ijk*trip_aaa[h] + lmn] += nb * dum;

                    for ( int p = 0; p < amo_; p++) {
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijkp = ibas_aaab_sym[h2][i][j][k][p];
                        int lmnp = ibas_aaab_sym[h2][l][m][n][p];
                        A[d4aaaboff[h2] + ijkp*quartet_aaab[h2]+lmnp] -= dum;
                    }
                }
            }
            offset += trip_aaa[h] * trip_aaa[h];
        }
        // D4aaab -> D3aab
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
                    A[d3aaboff[h] + ijk*trip_aab[h] + lmn] += (na - 2.0) * dum;
                    for ( int p = 0; p < amo_; p++) {
                        if ( i == p) continue;
                        if ( j == p) continue;
                        if ( l == p) continue;
                        if ( m == p) continue;
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijkp = ibas_aaab_sym[h2][i][j][p][k];
                        int lmnp = ibas_aaab_sym[h2][l][m][p][n];
                        int s = 1;
                        if ( p < i ) s = -s;
                        if ( p < j ) s = -s;
                        if ( p < l ) s = -s;
                        if ( p < m ) s = -s;
                        A[d4aaaboff[h2] + ijkp*quartet_aaab[h2]+lmnp] -= s * dum;
                    }
                }
            }
            offset += trip_aab[h] * trip_aab[h];
        }
    }

    // D4aabb
    if ( na > 1 && nb > 1 ) {
        // D4aabb -> D3aab
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
                    A[d3aaboff[h] + ijk*trip_aab[h] + lmn] += (nb - 1.0) * dum;
                    for ( int p = 0; p < amo_; p++) {
                        if ( k == p ) continue;
                        if ( n == p ) continue;
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijkp = ibas_aabb_sym[h2][i][j][k][p];
                        int lmnp = ibas_aabb_sym[h2][l][m][n][p];
                        int s = 1;
                        if ( p < k ) s = -s;
                        if ( p < n ) s = -s;
                        A[d4aabboff[h2] + ijkp*quartet_aabb[h2]+lmnp] -= s * dum;
                    }
                }
            }
            offset += trip_aab[h] * trip_aab[h];
        }
        // D4aabb -> D3bba
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
                    A[d3bbaoff[h] + ijk*trip_aab[h] + lmn] += (na - 1.0) * dum;
                    for ( int p = 0; p < amo_; p++) {
                        if ( k == p ) continue;
                        if ( n == p ) continue;
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijkp = ibas_aabb_sym[h2][k][p][i][j];
                        int lmnp = ibas_aabb_sym[h2][n][p][l][m];
                        int s = 1;
                        if ( p < k ) s = -s;
                        if ( p < n ) s = -s;
                        A[d4aabboff[h2] + ijkp*quartet_aabb[h2]+lmnp] -= s * dum;
                    }
                }
            }
            offset += trip_aab[h] * trip_aab[h];
        }
    }

    // D4bbba
    if ( na > 0 && nb > 2 ) {
        // D4bbba -> D3bbb
        for (int h = 0; h < nirrep_; h++) {
            for (int ijk = 0; ijk < trip_aaa[h]; ijk++) {
                int i = bas_aaa_sym[h][ijk][0];
                int j = bas_aaa_sym[h][ijk][1];
                int k = bas_aaa_sym[h][ijk][2];
                for (int lmn = 0; lmn < trip_aaa[h]; lmn++) {
                    int l = bas_aaa_sym[h][lmn][0];
                    int m = bas_aaa_sym[h][lmn][1];
                    int n = bas_aaa_sym[h][lmn][2];
                    double dum = u[offset + ijk*trip_aaa[h]+lmn];

                    A[d3bbboff[h] + ijk*trip_aaa[h] + lmn] += na * dum;

                    for ( int p = 0; p < amo_; p++) {
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijkp = ibas_aaab_sym[h2][i][j][k][p];
                        int lmnp = ibas_aaab_sym[h2][l][m][n][p];
                        A[d4bbbaoff[h2] + ijkp*quartet_aaab[h2]+lmnp] -= dum;
                    }
                }
            }
            offset += trip_aaa[h] * trip_aaa[h];
        }
        // D4bbba -> D3bba
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
                    A[d3bbaoff[h] + ijk*trip_aab[h] + lmn] += (nb - 2.0) * dum;
                    for ( int p = 0; p < amo_; p++) {
                        if ( i == p) continue;
                        if ( j == p) continue;
                        if ( l == p) continue;
                        if ( m == p) continue;
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijkp = ibas_aaab_sym[h2][i][j][p][k];
                        int lmnp = ibas_aaab_sym[h2][l][m][p][n];
                        int s = 1;
                        if ( p < i ) s = -s;
                        if ( p < j ) s = -s;
                        if ( p < l ) s = -s;
                        if ( p < m ) s = -s;
                        A[d4bbbaoff[h2] + ijkp*quartet_aaab[h2]+lmnp] -= s * dum;
                    }
                }
            }
            offset += trip_aab[h] * trip_aab[h];
        }
    }


    // D4bbbb -> D3bbb
    if ( nb > 3 ) {
        for (int h = 0; h < nirrep_; h++) {
            for (int ijk = 0; ijk < trip_aaa[h]; ijk++) {
                int i = bas_aaa_sym[h][ijk][0];
                int j = bas_aaa_sym[h][ijk][1];
                int k = bas_aaa_sym[h][ijk][2];
                for (int lmn = 0; lmn < trip_aaa[h]; lmn++) {
                    int l = bas_aaa_sym[h][lmn][0];
                    int m = bas_aaa_sym[h][lmn][1];
                    int n = bas_aaa_sym[h][lmn][2];
                    double dum = u[offset + ijk*trip_aaa[h]+lmn];
                    A[d3bbboff[h] + ijk*trip_aaa[h] + lmn] += (nb - 3.0) * dum;
                    for ( int p = 0; p < amo_; p++) {
                        if ( i == p) continue;
                        if ( j == p) continue;
                        if ( k == p) continue;
                        if ( l == p) continue;
                        if ( m == p) continue;
                        if ( n == p) continue;
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijkp = ibas_aaaa_sym[h2][i][j][k][p];
                        int lmnp = ibas_aaaa_sym[h2][l][m][n][p];
                        int s = 1;
                        if ( p < i ) s = -s;
                        if ( p < j ) s = -s;
                        if ( p < k ) s = -s;
                        if ( p < l ) s = -s;
                        if ( p < m ) s = -s;
                        if ( p < n ) s = -s;
                        A[d4bbbboff[h2] + ijkp*quartet_aaaa[h2]+lmnp] -= s * dum;
                    }
                }
            }
            offset += trip_aaa[h] * trip_aaa[h];
        }
    }
}

} // end namespaces
