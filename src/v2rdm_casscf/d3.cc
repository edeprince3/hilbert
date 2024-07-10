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

// D3 portion of A.u 
void v2RDMSolver::D3_constraints_Au(double * A,double * u){

    double na = nalpha_ - nrstc_ - nfrzc_;
    double nb = nbeta_ - nrstc_ - nfrzc_;

    if ( constrain_sz_ ) {
        // D3aaa -> D2aa
        //if ( na > 2 ) {
            for ( int h = 0; h < nirrep_; h++) {
                #pragma omp parallel for schedule (static)
                for ( int ij = 0; ij < gems_aa[h]; ij++) {
                    int i = bas_aa_sym[h][ij][0];
                    int j = bas_aa_sym[h][ij][1];
                    for ( int kl = 0; kl < gems_aa[h]; kl++) {
                        int k = bas_aa_sym[h][kl][0];
                        int l = bas_aa_sym[h][kl][1];
                        double dum = (na - 2.0) * u[d2aaoff[h] + ij*gems_aa[h] + kl];
                        for ( int p = 0; p < amo_; p++) {
                            if ( i == p || j == p ) continue;
                            if ( k == p || l == p ) continue;
                            int hp = symmetry[p];
                            int h2 = h ^ hp;
                            int ijp = ibas_aaa_sym[h2][i][j][p];
                            int klp = ibas_aaa_sym[h2][k][l][p];
                            int s = 1;
                            //if ( p < i ) s = -s;
                            //if ( p < j ) s = -s;
                            //if ( p < k ) s = -s;
                            //if ( p < l ) s = -s;
                            if ( i < p && p < j ) s = -s;
                            if ( k < p && p < l ) s = -s;
                            dum -= s * u[d3aaaoff[h2] + ijp*trip_aaa[h2]+klp];
                        }
                        A[offset + ij*gems_aa[h]+kl] = dum;
                    }
                }
                offset += gems_aa[h] * gems_aa[h];
            }
        //}
        //if ( nb > 2 ) {
            // D3bbb -> D2bb
            for ( int h = 0; h < nirrep_; h++) {
                #pragma omp parallel for schedule (static)
                for ( int ij = 0; ij < gems_aa[h]; ij++) {
                    int i = bas_aa_sym[h][ij][0];
                    int j = bas_aa_sym[h][ij][1];
                    for ( int kl = 0; kl < gems_aa[h]; kl++) {
                        int k = bas_aa_sym[h][kl][0];
                        int l = bas_aa_sym[h][kl][1];
                        double dum = (nb - 2.0) * u[d2bboff[h] + ij*gems_aa[h] + kl];
                        for ( int p = 0; p < amo_; p++) {
                            if ( i == p || j == p ) continue;
                            if ( k == p || l == p ) continue;
                            int hp = symmetry[p];
                            int h2 = h ^ hp;
                            int ijp = ibas_aaa_sym[h2][i][j][p];
                            int klp = ibas_aaa_sym[h2][k][l][p];
                            int s = 1;
                            //if ( p < i ) s = -s;
                            //if ( p < j ) s = -s;
                            //if ( p < k ) s = -s;
                            //if ( p < l ) s = -s;
                            if ( i < p && p < j ) s = -s;
                            if ( k < p && p < l ) s = -s;
                            dum -= s * u[d3bbboff[h2] + ijp*trip_aaa[h2]+klp];
                        }
                        A[offset + ij*gems_aa[h]+kl] = dum;
                    }
                }
                offset += gems_aa[h] * gems_aa[h];
            }
        //}
        // D3aab -> D2aa
        //if ( na > 1 ) {
            for ( int h = 0; h < nirrep_; h++) {
                #pragma omp parallel for schedule (static)
                for ( int ij = 0; ij < gems_aa[h]; ij++) {
                    int i = bas_aa_sym[h][ij][0];
                    int j = bas_aa_sym[h][ij][1];
                    for ( int kl = 0; kl < gems_aa[h]; kl++) {
                        int k = bas_aa_sym[h][kl][0];
                        int l = bas_aa_sym[h][kl][1];
                        double dum = nb * u[d2aaoff[h] + ij*gems_aa[h] + kl];
                        for ( int p = 0; p < amo_; p++) {
                            int hp = symmetry[p];
                            int h2 = h ^ hp;
                            int ijp = ibas_aab_sym[h2][i][j][p];
                            int klp = ibas_aab_sym[h2][k][l][p];
                            dum -= u[d3aaboff[h2] + ijp*trip_aab[h2]+klp];
                        }
                        A[offset + ij*gems_aa[h]+kl] = dum;
                    }
                }
                offset += gems_aa[h] * gems_aa[h];
            }
        //}
        //if ( nb > 1 ) {
            // D3bba -> D2bb
            for ( int h = 0; h < nirrep_; h++) {
                #pragma omp parallel for schedule (static)
                for ( int ij = 0; ij < gems_aa[h]; ij++) {
                    int i = bas_aa_sym[h][ij][0];
                    int j = bas_aa_sym[h][ij][1];
                    for ( int kl = 0; kl < gems_aa[h]; kl++) {
                        int k = bas_aa_sym[h][kl][0];
                        int l = bas_aa_sym[h][kl][1];
                        double dum = na * u[d2bboff[h] + ij*gems_aa[h] + kl];
                        for ( int p = 0; p < amo_; p++) {
                            int hp = symmetry[p];
                            int h2 = h ^ hp;
                            int ijp = ibas_aab_sym[h2][i][j][p];
                            int klp = ibas_aab_sym[h2][k][l][p];
                            dum -= u[d3bbaoff[h2] + ijp*trip_aab[h2]+klp];
                        }
                        A[offset + ij*gems_aa[h]+kl] = dum;
                    }
                }
                offset += gems_aa[h] * gems_aa[h];
            }
        //}
        //if ( na > 1 ) {
            // D3aab -> D2ab
            for ( int h = 0; h < nirrep_; h++) {
                #pragma omp parallel for schedule (static)
                for ( int ij = 0; ij < gems_ab[h]; ij++) {
                    int i = bas_ab_sym[h][ij][0];
                    int j = bas_ab_sym[h][ij][1];
                    for ( int kl = 0; kl < gems_ab[h]; kl++) {
                        int k = bas_ab_sym[h][kl][0];
                        int l = bas_ab_sym[h][kl][1];
                        double dum = (na - 1.0) * u[d2aboff[h] + ij*gems_ab[h] + kl];
                        for ( int p = 0; p < amo_; p++) {
                            if ( i == p) continue;
                            if ( k == p) continue;
                            int h2 = SymmetryPair(h,symmetry[p]);
                            int ijp = ibas_aab_sym[h2][i][p][j];
                            int klp = ibas_aab_sym[h2][k][p][l];
                            int s = 1;
                            if ( p < i ) s = -s;
                            if ( p < k ) s = -s;
                            dum -= s * u[d3aaboff[h2] + ijp*trip_aab[h2]+klp];
                        }
                        A[offset + ij*gems_ab[h]+kl] = dum;
                    }
                }
                offset += gems_ab[h] * gems_ab[h];
            }
        //}
        //if ( nb > 1 ) {
            // D3bba -> D2ab
            for ( int h = 0; h < nirrep_; h++) {
                #pragma omp parallel for schedule (static)
                for ( int ij = 0; ij < gems_ab[h]; ij++) {
                    int i = bas_ab_sym[h][ij][0];
                    int j = bas_ab_sym[h][ij][1];
                    for ( int kl = 0; kl < gems_ab[h]; kl++) {
                        int k = bas_ab_sym[h][kl][0];
                        int l = bas_ab_sym[h][kl][1];
                        double dum = (nb - 1.0) * u[d2aboff[h] + ij*gems_ab[h] + kl];
                        for ( int p = 0; p < amo_; p++) {
                            if ( j == p) continue;
                            if ( l == p) continue;
                            int h2 = SymmetryPair(h,symmetry[p]);
                            int ijp = ibas_aab_sym[h2][j][p][i];
                            int klp = ibas_aab_sym[h2][l][p][k];
                            int s = 1;
                            if ( p < j ) s = -s;
                            if ( p < l ) s = -s;
                            dum -= s * u[d3bbaoff[h2] + ijp*trip_aab[h2]+klp];
                        }
                        A[offset + ij*gems_ab[h]+kl] = dum;
                    }
                }
                offset += gems_ab[h] * gems_ab[h];
            }
        //}
    }else {
        // D3aaa + D3aab -> D2aa
        for ( int h = 0; h < nirrep_; h++) {
            #pragma omp parallel for schedule (static)
            for ( int ij = 0; ij < gems_aa[h]; ij++) {
                int i = bas_aa_sym[h][ij][0];
                int j = bas_aa_sym[h][ij][1];
                for ( int kl = 0; kl < gems_aa[h]; kl++) {
                    int k = bas_aa_sym[h][kl][0];
                    int l = bas_aa_sym[h][kl][1];
                    double dum = (na + nb - 2.0) * u[d2aaoff[h] + ij*gems_aa[h] + kl];
                    for ( int p = 0; p < amo_; p++) {
                        int hp = symmetry[p];
                        int h2 = h ^ hp;
                        int ijp = ibas_aab_sym[h2][i][j][p];
                        int klp = ibas_aab_sym[h2][k][l][p];
                        dum -= u[d3aaboff[h2] + ijp*trip_aab[h2]+klp];

                        if ( i == p || j == p ) continue;
                        if ( k == p || l == p ) continue;

                        ijp = ibas_aaa_sym[h2][i][j][p];
                        klp = ibas_aaa_sym[h2][k][l][p];
                        int s = 1;
                        //if ( p < i ) s = -s;
                        //if ( p < j ) s = -s;
                        //if ( p < k ) s = -s;
                        //if ( p < l ) s = -s;
                        if ( i < p && p < j ) s = -s;
                        if ( k < p && p < l ) s = -s;
                        dum -= s * u[d3aaaoff[h2] + ijp*trip_aaa[h2]+klp];
                    }
                    A[offset + ij*gems_aa[h]+kl] = dum;
                }
            }
            offset += gems_aa[h] * gems_aa[h];
        }
        // D3bbb + D3bba -> D2bb
        for ( int h = 0; h < nirrep_; h++) {
            #pragma omp parallel for schedule (static)
            for ( int ij = 0; ij < gems_aa[h]; ij++) {
                int i = bas_aa_sym[h][ij][0];
                int j = bas_aa_sym[h][ij][1];
                for ( int kl = 0; kl < gems_aa[h]; kl++) {
                    int k = bas_aa_sym[h][kl][0];
                    int l = bas_aa_sym[h][kl][1];
                    double dum = (na + nb - 2.0) * u[d2bboff[h] + ij*gems_aa[h] + kl];
                    for ( int p = 0; p < amo_; p++) {

                        int hp = symmetry[p];
                        int h2 = h ^ hp;
                        int ijp = ibas_aab_sym[h2][i][j][p];
                        int klp = ibas_aab_sym[h2][k][l][p];
                        dum -= u[d3bbaoff[h2] + ijp*trip_aab[h2]+klp];

                        if ( i == p || j == p ) continue;
                        if ( k == p || l == p ) continue;

                        ijp = ibas_aaa_sym[h2][i][j][p];
                        klp = ibas_aaa_sym[h2][k][l][p];
                        int s = 1;
                        //if ( p < i ) s = -s;
                        //if ( p < j ) s = -s;
                        //if ( p < k ) s = -s;
                        //if ( p < l ) s = -s;
                        if ( i < p && p < j ) s = -s;
                        if ( k < p && p < l ) s = -s;
                        dum -= s * u[d3bbboff[h2] + ijp*trip_aaa[h2]+klp];
                    }
                    A[offset + ij*gems_aa[h]+kl] = dum;
                }
            }
            offset += gems_aa[h] * gems_aa[h];
        }
        // D3abb + D3bba -> D2ab
        for ( int h = 0; h < nirrep_; h++) {
            #pragma omp parallel for schedule (static)
            for ( int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                for ( int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    double dum = (na + nb - 2.0) * u[d2aboff[h] + ij*gems_ab[h] + kl];
                    for ( int p = 0; p < amo_; p++) {
                        int hp = symmetry[p];
                        int h2 = h ^ hp;
                        if (i != p && k != p) {
                            int ijp = ibas_aab_sym[h2][i][p][j];
                            int klp = ibas_aab_sym[h2][k][p][l];
                            int s = 1;
                            if ( p < i ) s = -s;
                            if ( p < k ) s = -s;
                            dum -= s * u[d3aaboff[h2] + ijp*trip_aab[h2]+klp];
                        }
                        if ( j != p && l != p) {
                            int ijp = ibas_aab_sym[h2][j][p][i];
                            int klp = ibas_aab_sym[h2][l][p][k];
                            int s = 1;
                            if ( p < j ) s = -s;
                            if ( p < l ) s = -s;
                            dum -= s * u[d3bbaoff[h2] + ijp*trip_aab[h2]+klp];
                        }
                    }
                    A[offset + ij*gems_ab[h]+kl] = dum;
                }
            }
            offset += gems_ab[h] * gems_ab[h];
        }
    } 
                 
    // additional spin constraints for singlets:
    if ( constrain_spin_ && nalpha_ == nbeta_ ) {
        // D3aab = D3bba
        for ( int h = 0; h < nirrep_; h++) {
            C_DCOPY(trip_aab[h]*trip_aab[h],u + d3aaboff[h],1,A + offset,1);
            C_DAXPY(trip_aab[h]*trip_aab[h],-1.0,u + d3bbaoff[h],1,A + offset,1);
            offset += trip_aab[h]*trip_aab[h];
        }    
        // D3aaa <- D3aab
        for ( int h = 0; h < nirrep_; h++) {
            C_DCOPY(trip_aaa[h]*trip_aaa[h],u + d3aaaoff[h],1,A + offset,1);
            for (int pqr = 0; pqr < trip_aaa[h]; pqr++) {
                int p = bas_aaa_sym[h][pqr][0];
                int q = bas_aaa_sym[h][pqr][1];
                int r = bas_aaa_sym[h][pqr][2];
                int pqr_b = ibas_aab_sym[h][p][q][r];
                int prq_b = ibas_aab_sym[h][p][r][q];
                int qrp_b = ibas_aab_sym[h][q][r][p];
                for (int stv = 0; stv < trip_aaa[h]; stv++) {
                    int s = bas_aaa_sym[h][stv][0];
                    int t = bas_aaa_sym[h][stv][1];
                    int v = bas_aaa_sym[h][stv][2];
                    int stv_b = ibas_aab_sym[h][s][t][v];
                    int svt_b = ibas_aab_sym[h][s][v][t];
                    int tvs_b = ibas_aab_sym[h][t][v][s];
                    A[offset + pqr*trip_aaa[h] + stv] -= 1.0/3.0 * u[d3aaboff[h] + pqr_b * trip_aab[h] + stv_b];
                    A[offset + pqr*trip_aaa[h] + stv] += 1.0/3.0 * u[d3aaboff[h] + pqr_b * trip_aab[h] + svt_b];
                A[offset + pqr*trip_aaa[h] + stv] -= 1.0/3.0 * u[d3aaboff[h] + pqr_b * trip_aab[h] + tvs_b];

                A[offset + pqr*trip_aaa[h] + stv] += 1.0/3.0 * u[d3aaboff[h] + prq_b * trip_aab[h] + stv_b];
                A[offset + pqr*trip_aaa[h] + stv] -= 1.0/3.0 * u[d3aaboff[h] + prq_b * trip_aab[h] + svt_b];
                A[offset + pqr*trip_aaa[h] + stv] += 1.0/3.0 * u[d3aaboff[h] + prq_b * trip_aab[h] + tvs_b];

                A[offset + pqr*trip_aaa[h] + stv] -= 1.0/3.0 * u[d3aaboff[h] + qrp_b * trip_aab[h] + stv_b];
                A[offset + pqr*trip_aaa[h] + stv] += 1.0/3.0 * u[d3aaboff[h] + qrp_b * trip_aab[h] + svt_b];
                A[offset + pqr*trip_aaa[h] + stv] -= 1.0/3.0 * u[d3aaboff[h] + qrp_b * trip_aab[h] + tvs_b];
                }
            }
            offset += trip_aaa[h]*trip_aaa[h];
        }
        // D3bbb <- D3bba
        for ( int h = 0; h < nirrep_; h++) {
            C_DCOPY(trip_aaa[h]*trip_aaa[h],u + d3bbboff[h],1,A + offset,1);
            for (int pqr = 0; pqr < trip_aaa[h]; pqr++) {
                int p = bas_aaa_sym[h][pqr][0];
                int q = bas_aaa_sym[h][pqr][1];
                int r = bas_aaa_sym[h][pqr][2];
                int pqr_b = ibas_aab_sym[h][p][q][r];
                int prq_b = ibas_aab_sym[h][p][r][q];
                int qrp_b = ibas_aab_sym[h][q][r][p];
                for (int stv = 0; stv < trip_aaa[h]; stv++) {
                    int s = bas_aaa_sym[h][stv][0];
                    int t = bas_aaa_sym[h][stv][1];
                    int v = bas_aaa_sym[h][stv][2];
                    int stv_b = ibas_aab_sym[h][s][t][v];
                    int svt_b = ibas_aab_sym[h][s][v][t];
                    int tvs_b = ibas_aab_sym[h][t][v][s];
                    A[offset + pqr*trip_aaa[h] + stv] -= 1.0/3.0 * u[d3bbaoff[h] + pqr_b * trip_aab[h] + stv_b];
                    A[offset + pqr*trip_aaa[h] + stv] += 1.0/3.0 * u[d3bbaoff[h] + pqr_b * trip_aab[h] + svt_b];
                    A[offset + pqr*trip_aaa[h] + stv] -= 1.0/3.0 * u[d3bbaoff[h] + pqr_b * trip_aab[h] + tvs_b];

                    A[offset + pqr*trip_aaa[h] + stv] += 1.0/3.0 * u[d3bbaoff[h] + prq_b * trip_aab[h] + stv_b];
                    A[offset + pqr*trip_aaa[h] + stv] -= 1.0/3.0 * u[d3bbaoff[h] + prq_b * trip_aab[h] + svt_b];
                    A[offset + pqr*trip_aaa[h] + stv] += 1.0/3.0 * u[d3bbaoff[h] + prq_b * trip_aab[h] + tvs_b];

                    A[offset + pqr*trip_aaa[h] + stv] -= 1.0/3.0 * u[d3bbaoff[h] + qrp_b * trip_aab[h] + stv_b];
                    A[offset + pqr*trip_aaa[h] + stv] += 1.0/3.0 * u[d3bbaoff[h] + qrp_b * trip_aab[h] + svt_b];
                    A[offset + pqr*trip_aaa[h] + stv] -= 1.0/3.0 * u[d3bbaoff[h] + qrp_b * trip_aab[h] + tvs_b];
                }
            }
            offset += trip_aaa[h]*trip_aaa[h];
        }
    }
}

// D3 portion of A^T.y 
void v2RDMSolver::D3_constraints_ATu(double * A,double * u){

    double na = nalpha_ - nrstc_ - nfrzc_;
    double nb = nbeta_ - nrstc_ - nfrzc_;

    if ( constrain_sz_ ) {
        //if ( na > 2 ) {
            // D3aaa -> D2aa
            for ( int h = 0; h < nirrep_; h++) {
                for ( int ij = 0; ij < gems_aa[h]; ij++) {
                    int i = bas_aa_sym[h][ij][0];
                    int j = bas_aa_sym[h][ij][1];
                    for ( int kl = 0; kl < gems_aa[h]; kl++) {
                        int k = bas_aa_sym[h][kl][0];
                        int l = bas_aa_sym[h][kl][1];
                        double dum = u[offset + ij*gems_aa[h] + kl];
                        A[d2aaoff[h] + ij*gems_aa[h] + kl] += (na - 2.0) * dum;
                        for ( int p = 0; p < amo_; p++) {
                            if ( i == p || j == p ) continue;
                            if ( k == p || l == p ) continue;
                            int hp = symmetry[p];
                            int h2 = h ^ hp;
                            int ijp = ibas_aaa_sym[h2][i][j][p];
                            int klp = ibas_aaa_sym[h2][k][l][p];
                            int s = 1;
                            //if ( p < i ) s = -s;
                            //if ( p < j ) s = -s;
                            //if ( p < k ) s = -s;
                            //if ( p < l ) s = -s;
                            if ( i < p && p < j ) s = -s;
                            if ( k < p && p < l ) s = -s;
                            A[d3aaaoff[h2] + ijp*trip_aaa[h2]+klp] -= s * dum;
                        }
                    }
                }
                offset += gems_aa[h] * gems_aa[h];
            }
        //}
        //if ( nb > 2 ) {
            // D3bbb -> D2bb
            for ( int h = 0; h < nirrep_; h++) {
                for ( int ij = 0; ij < gems_aa[h]; ij++) {
                    int i = bas_aa_sym[h][ij][0];
                    int j = bas_aa_sym[h][ij][1];
                    for ( int kl = 0; kl < gems_aa[h]; kl++) {
                        int k = bas_aa_sym[h][kl][0];
                        int l = bas_aa_sym[h][kl][1];
                        double dum = u[offset + ij*gems_aa[h] + kl];
                        A[d2bboff[h] + ij*gems_aa[h] + kl] += (nb - 2.0) * dum;
                        for ( int p = 0; p < amo_; p++) {
                            if ( i == p || j == p ) continue;
                            if ( k == p || l == p ) continue;
                            int hp = symmetry[p];
                            int h2 = h ^ hp;
                            int ijp = ibas_aaa_sym[h2][i][j][p];
                            int klp = ibas_aaa_sym[h2][k][l][p];
                            int s = 1;
                            //if ( p < i ) s = -s;
                            //if ( p < j ) s = -s;
                            //if ( p < k ) s = -s;
                            //if ( p < l ) s = -s;
                            if ( i < p && p < j ) s = -s;
                            if ( k < p && p < l ) s = -s;
                            A[d3bbboff[h2] + ijp*trip_aaa[h2]+klp] -= s * dum;
                        }
                    }
                }
                offset += gems_aa[h] * gems_aa[h];
            }
        //}
        // D3aab -> D2aa
        //if ( na > 1 ) {
            for ( int h = 0; h < nirrep_; h++) {
                for ( int ij = 0; ij < gems_aa[h]; ij++) {
                    int i = bas_aa_sym[h][ij][0];
                    int j = bas_aa_sym[h][ij][1];
                    for ( int kl = 0; kl < gems_aa[h]; kl++) {
                        int k = bas_aa_sym[h][kl][0];
                        int l = bas_aa_sym[h][kl][1];
                        double dum = u[offset + ij*gems_aa[h] + kl];
                        A[d2aaoff[h] + ij*gems_aa[h] + kl] += nb * dum;
                        for ( int p = 0; p < amo_; p++) {
                            int hp = symmetry[p];
                            int h2 = h ^ hp;
                            int ijp = ibas_aab_sym[h2][i][j][p];
                            int klp = ibas_aab_sym[h2][k][l][p];
                            A[d3aaboff[h2] + ijp*trip_aab[h2]+klp] -= dum;
                        }
                    }
                }
                offset += gems_aa[h] * gems_aa[h];
            }
        //}
        //if ( nb > 1 ) {
            // D3bba -> D2bb
            for ( int h = 0; h < nirrep_; h++) {
                for ( int ij = 0; ij < gems_aa[h]; ij++) {
                    int i = bas_aa_sym[h][ij][0];
                    int j = bas_aa_sym[h][ij][1];
                    for ( int kl = 0; kl < gems_aa[h]; kl++) {
                        int k = bas_aa_sym[h][kl][0];
                        int l = bas_aa_sym[h][kl][1];
                        double dum = u[offset + ij*gems_aa[h] + kl];
                        A[d2bboff[h] + ij*gems_aa[h] + kl] += na * dum;
                        for ( int p = 0; p < amo_; p++) {
                            int hp = symmetry[p];
                            int h2 = h ^ hp;
                            int ijp = ibas_aab_sym[h2][i][j][p];
                            int klp = ibas_aab_sym[h2][k][l][p];
                            A[d3bbaoff[h2] + ijp*trip_aab[h2]+klp] -= dum;
                        }
                    }
                }
                offset += gems_aa[h] * gems_aa[h];
            }
        //}
        //if ( na > 1 ) {
            // D3aab -> D2ab
            for ( int h = 0; h < nirrep_; h++) {
                for ( int ij = 0; ij < gems_ab[h]; ij++) {
                    int i = bas_ab_sym[h][ij][0];
                    int j = bas_ab_sym[h][ij][1];
                    for ( int kl = 0; kl < gems_ab[h]; kl++) {
                        int k = bas_ab_sym[h][kl][0];
                        int l = bas_ab_sym[h][kl][1];
                        double dum = u[offset + ij*gems_ab[h] + kl];
                        A[d2aboff[h] + ij*gems_ab[h] + kl] += (na - 1.0) * dum;
                        for ( int p = 0; p < amo_; p++) {
                            if ( i == p) continue;
                            if ( k == p) continue;
                            int hp = symmetry[p];
                            int h2 = h ^ hp;
                            int ijp = ibas_aab_sym[h2][i][p][j];
                            int klp = ibas_aab_sym[h2][k][p][l];
                            int s = 1;
                            if ( p < i ) s = -s;
                            if ( p < k ) s = -s;
                            A[d3aaboff[h2] + ijp*trip_aab[h2]+klp] -= s * dum;
                        }
                    }
                }
                offset += gems_ab[h] * gems_ab[h];
            }
        //}
        //if ( nb > 1 ) {
            // D3bba -> D2ab
            for ( int h = 0; h < nirrep_; h++) {
                for ( int ij = 0; ij < gems_ab[h]; ij++) {
                    int i = bas_ab_sym[h][ij][0];
                    int j = bas_ab_sym[h][ij][1];
                    for ( int kl = 0; kl < gems_ab[h]; kl++) {
                        int k = bas_ab_sym[h][kl][0];
                        int l = bas_ab_sym[h][kl][1];
                        double dum = u[offset + ij*gems_ab[h] + kl];
                        A[d2aboff[h] + ij*gems_ab[h] + kl] += (nb - 1.0) * dum;
                        for ( int p = 0; p < amo_; p++) {
                            if ( j == p) continue;
                            if ( l == p) continue;
                            int hp = symmetry[p];
                            int h2 = h ^ hp;
                            int ijp = ibas_aab_sym[h2][j][p][i];
                            int klp = ibas_aab_sym[h2][l][p][k];
                            int s = 1;
                            if ( p < j ) s = -s;
                            if ( p < l ) s = -s;
                            A[d3bbaoff[h2] + ijp*trip_aab[h2]+klp] -= s * dum;
                        }
                    }
                }
                offset += gems_ab[h] * gems_ab[h];
            }
        //}
    }else {
        // D3aaa + D3aab -> D2aa
        for ( int h = 0; h < nirrep_; h++) {
            for ( int ij = 0; ij < gems_aa[h]; ij++) {
                int i = bas_aa_sym[h][ij][0];
                int j = bas_aa_sym[h][ij][1];
                for ( int kl = 0; kl < gems_aa[h]; kl++) {
                    int k = bas_aa_sym[h][kl][0];
                    int l = bas_aa_sym[h][kl][1];
                    double dum = u[offset + ij*gems_aa[h] + kl];
                    A[d2aaoff[h] + ij*gems_aa[h] + kl] += (na + nb - 2.0) * dum;
                    for ( int p = 0; p < amo_; p++) {

                        int hp = symmetry[p];
                        int h2 = h ^ hp;
                        int ijp = ibas_aab_sym[h2][i][j][p];
                        int klp = ibas_aab_sym[h2][k][l][p];
                        A[d3aaboff[h2] + ijp*trip_aab[h2]+klp] -= dum;

                        if ( i == p || j == p ) continue;
                        if ( k == p || l == p ) continue;

                        ijp = ibas_aaa_sym[h2][i][j][p];
                        klp = ibas_aaa_sym[h2][k][l][p];
                        int s = 1;
                        //if ( p < i ) s = -s;
                        //if ( p < j ) s = -s;
                        //if ( p < k ) s = -s;
                        //if ( p < l ) s = -s;
                        if ( i < p && p < j ) s = -s;
                        if ( k < p && p < l ) s = -s;
                        A[d3aaaoff[h2] + ijp*trip_aaa[h2]+klp] -= s * dum;
                    }
                }
            }
            offset += gems_aa[h] * gems_aa[h];
        }
        // D3bbb + D3bba -> D2bb
        for ( int h = 0; h < nirrep_; h++) {
            for ( int ij = 0; ij < gems_aa[h]; ij++) {
                int i = bas_aa_sym[h][ij][0];
                int j = bas_aa_sym[h][ij][1];
                for ( int kl = 0; kl < gems_aa[h]; kl++) {
                    int k = bas_aa_sym[h][kl][0];
                    int l = bas_aa_sym[h][kl][1];
                    double dum = u[offset + ij*gems_aa[h] + kl];
                    A[d2bboff[h] + ij*gems_aa[h] + kl] += (na + nb - 2.0) * dum;
                    for ( int p = 0; p < amo_; p++) {
                        int hp = symmetry[p];
                        int h2 = h ^ hp;
                        int ijp = ibas_aab_sym[h2][i][j][p];
                        int klp = ibas_aab_sym[h2][k][l][p];
                        A[d3bbaoff[h2] + ijp*trip_aab[h2]+klp] -= dum;

                        if ( i == p || j == p ) continue;
                        if ( k == p || l == p ) continue;

                        ijp = ibas_aaa_sym[h2][i][j][p];
                        klp = ibas_aaa_sym[h2][k][l][p];
                        int s = 1;
                        //if ( p < i ) s = -s;
                        //if ( p < j ) s = -s;
                        //if ( p < k ) s = -s;
                        //if ( p < l ) s = -s;
                        if ( i < p && p < j ) s = -s;
                        if ( k < p && p < l ) s = -s;
                        A[d3bbboff[h2] + ijp*trip_aaa[h2]+klp] -= s * dum;
                    }
                }
            }
            offset += gems_aa[h] * gems_aa[h];
        }
        // D3abb + D3bba -> D2ab
        for ( int h = 0; h < nirrep_; h++) {
            for ( int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                for ( int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    double dum = u[offset + ij*gems_ab[h] + kl];
                    A[d2aboff[h] + ij*gems_ab[h] + kl] += (na + nb - 2.0) * dum;
                    for ( int p = 0; p < amo_; p++) {
                        if (i != p && k != p) {
                            int hp = symmetry[p];
                            int h2 = h ^ hp;
                            int ijp = ibas_aab_sym[h2][i][p][j];
                            int klp = ibas_aab_sym[h2][k][p][l];
                            int s = 1;
                            if ( p < i ) s = -s;
                            if ( p < k ) s = -s;
                            A[d3aaboff[h2] + ijp*trip_aab[h2]+klp] -= s * dum;
                        }

                        if (j != p && l != p) {
                            int hp = symmetry[p];
                            int h2 = h ^ hp;
                            int ijp = ibas_aab_sym[h2][j][p][i];
                            int klp = ibas_aab_sym[h2][l][p][k];
                            int s = 1;
                            if ( p < j ) s = -s;
                            if ( p < l ) s = -s;
                            A[d3bbaoff[h2] + ijp*trip_aab[h2]+klp] -= s * dum;
                        }
                    }
                }
            }
            offset += gems_ab[h] * gems_ab[h];
        }
    }

    // additional spin constraints for singlets:
    if ( constrain_spin_ && nalpha_ == nbeta_ ) {
        // D3aab = D3bba
        for ( int h = 0; h < nirrep_; h++) {
            C_DAXPY(trip_aab[h]*trip_aab[h], 1.0,u + offset,1,A + d3aaboff[h],1);
            C_DAXPY(trip_aab[h]*trip_aab[h],-1.0,u + offset,1,A + d3bbaoff[h],1);
            offset += trip_aab[h]*trip_aab[h];
        }
        // D3aaa <- D3aab
        for ( int h = 0; h < nirrep_; h++) {
            C_DAXPY(trip_aaa[h]*trip_aaa[h],1.0,u + offset,1,A+d3aaaoff[h],1);
            for (int pqr = 0; pqr < trip_aaa[h]; pqr++) {
                int p = bas_aaa_sym[h][pqr][0];
                int q = bas_aaa_sym[h][pqr][1];
                int r = bas_aaa_sym[h][pqr][2];
                int pqr_b = ibas_aab_sym[h][p][q][r];
                int prq_b = ibas_aab_sym[h][p][r][q];
                int qrp_b = ibas_aab_sym[h][q][r][p];
                for (int stv = 0; stv < trip_aaa[h]; stv++) {
                    int s = bas_aaa_sym[h][stv][0];
                    int t = bas_aaa_sym[h][stv][1];
                    int v = bas_aaa_sym[h][stv][2];
                    int stv_b = ibas_aab_sym[h][s][t][v];
                    int svt_b = ibas_aab_sym[h][s][v][t];
                    int tvs_b = ibas_aab_sym[h][t][v][s];
                    A[d3aaboff[h] + pqr_b * trip_aab[h] + stv_b] -= 1.0/3.0 * u[offset + pqr*trip_aaa[h] + stv];
                    A[d3aaboff[h] + pqr_b * trip_aab[h] + svt_b] += 1.0/3.0 * u[offset + pqr*trip_aaa[h] + stv];
                    A[d3aaboff[h] + pqr_b * trip_aab[h] + tvs_b] -= 1.0/3.0 * u[offset + pqr*trip_aaa[h] + stv];
                                                                                                                   
                    A[d3aaboff[h] + prq_b * trip_aab[h] + stv_b] += 1.0/3.0 * u[offset + pqr*trip_aaa[h] + stv];
                    A[d3aaboff[h] + prq_b * trip_aab[h] + svt_b] -= 1.0/3.0 * u[offset + pqr*trip_aaa[h] + stv];
                    A[d3aaboff[h] + prq_b * trip_aab[h] + tvs_b] += 1.0/3.0 * u[offset + pqr*trip_aaa[h] + stv];
                                                                                                                   
                    A[d3aaboff[h] + qrp_b * trip_aab[h] + stv_b] -= 1.0/3.0 * u[offset + pqr*trip_aaa[h] + stv];
                    A[d3aaboff[h] + qrp_b * trip_aab[h] + svt_b] += 1.0/3.0 * u[offset + pqr*trip_aaa[h] + stv];
                    A[d3aaboff[h] + qrp_b * trip_aab[h] + tvs_b] -= 1.0/3.0 * u[offset + pqr*trip_aaa[h] + stv];
                }
            }
            offset += trip_aaa[h]*trip_aaa[h];
        }
        // D3bbb <- D3bba
        for ( int h = 0; h < nirrep_; h++) {
            C_DAXPY(trip_aaa[h]*trip_aaa[h],1.0,u + offset,1,A+d3bbboff[h],1);
            for (int pqr = 0; pqr < trip_aaa[h]; pqr++) {
                int p = bas_aaa_sym[h][pqr][0];
                int q = bas_aaa_sym[h][pqr][1];
                int r = bas_aaa_sym[h][pqr][2];
                int pqr_b = ibas_aab_sym[h][p][q][r];
                int prq_b = ibas_aab_sym[h][p][r][q];
                int qrp_b = ibas_aab_sym[h][q][r][p];
                for (int stv = 0; stv < trip_aaa[h]; stv++) {
                    int s = bas_aaa_sym[h][stv][0];
                    int t = bas_aaa_sym[h][stv][1];
                    int v = bas_aaa_sym[h][stv][2];
                    int stv_b = ibas_aab_sym[h][s][t][v];
                    int svt_b = ibas_aab_sym[h][s][v][t];
                    int tvs_b = ibas_aab_sym[h][t][v][s];
                    A[d3bbaoff[h] + pqr_b * trip_aab[h] + stv_b] -= 1.0/3.0 * u[offset + pqr*trip_aaa[h] + stv];
                    A[d3bbaoff[h] + pqr_b * trip_aab[h] + svt_b] += 1.0/3.0 * u[offset + pqr*trip_aaa[h] + stv];
                    A[d3bbaoff[h] + pqr_b * trip_aab[h] + tvs_b] -= 1.0/3.0 * u[offset + pqr*trip_aaa[h] + stv];
                                                                                                                   
                    A[d3bbaoff[h] + prq_b * trip_aab[h] + stv_b] += 1.0/3.0 * u[offset + pqr*trip_aaa[h] + stv];
                    A[d3bbaoff[h] + prq_b * trip_aab[h] + svt_b] -= 1.0/3.0 * u[offset + pqr*trip_aaa[h] + stv];
                    A[d3bbaoff[h] + prq_b * trip_aab[h] + tvs_b] += 1.0/3.0 * u[offset + pqr*trip_aaa[h] + stv];
                                                                                                                   
                    A[d3bbaoff[h] + qrp_b * trip_aab[h] + stv_b] -= 1.0/3.0 * u[offset + pqr*trip_aaa[h] + stv];
                    A[d3bbaoff[h] + qrp_b * trip_aab[h] + svt_b] += 1.0/3.0 * u[offset + pqr*trip_aaa[h] + stv];
                    A[d3bbaoff[h] + qrp_b * trip_aab[h] + tvs_b] -= 1.0/3.0 * u[offset + pqr*trip_aaa[h] + stv];
                }
            }
            offset += trip_aaa[h]*trip_aaa[h];
        }
    }
}

} // end namespaces
