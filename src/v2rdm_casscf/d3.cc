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

// D3 portion of A.u 
void v2RDMSolver::D3_constraints_Au(SharedVector A,SharedVector u){

    double * A_p = A->pointer();
    double * u_p = u->pointer();

    int na = nalpha_ - nrstc_ - nfrzc_;
    int nb = nbeta_ - nrstc_ - nfrzc_;

    // D3aaa -> D2aa
    if ( na > 2 ) {
        for ( int h = 0; h < nirrep_; h++) {
            #pragma omp parallel for schedule (static)
            for ( int ij = 0; ij < gems_aa[h]; ij++) {
                int i = bas_aa_sym[h][ij][0];
                int j = bas_aa_sym[h][ij][1];
                for ( int kl = 0; kl < gems_aa[h]; kl++) {
                    int k = bas_aa_sym[h][kl][0];
                    int l = bas_aa_sym[h][kl][1];
                    double dum = (na - 2.0) * u_p[d2aaoff[h] + ij*gems_aa[h] + kl];
                    for ( int p = 0; p < amo_; p++) {
                        if ( i == p || j == p ) continue;
                        if ( k == p || l == p ) continue;
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijp = ibas_aaa_sym[h2][i][j][p];
                        int klp = ibas_aaa_sym[h2][k][l][p];
                        int s = 1;
                        if ( p < i ) s = -s;
                        if ( p < j ) s = -s;
                        if ( p < k ) s = -s;
                        if ( p < l ) s = -s;
                        dum -= s * u_p[d3aaaoff[h2] + ijp*trip_aaa[h2]+klp];
                    }
                    A_p[offset + ij*gems_aa[h]+kl] = dum;
                }
            }
            offset += gems_aa[h] * gems_aa[h];
        }
    }
    if ( nb > 2 ) {
        // D3bbb -> D2bb
        for ( int h = 0; h < nirrep_; h++) {
            #pragma omp parallel for schedule (static)
            for ( int ij = 0; ij < gems_aa[h]; ij++) {
                int i = bas_aa_sym[h][ij][0];
                int j = bas_aa_sym[h][ij][1];
                for ( int kl = 0; kl < gems_aa[h]; kl++) {
                    int k = bas_aa_sym[h][kl][0];
                    int l = bas_aa_sym[h][kl][1];
                    double dum = (nb - 2.0) * u_p[d2bboff[h] + ij*gems_aa[h] + kl];
                    for ( int p = 0; p < amo_; p++) {
                        if ( i == p || j == p ) continue;
                        if ( k == p || l == p ) continue;
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijp = ibas_aaa_sym[h2][i][j][p];
                        int klp = ibas_aaa_sym[h2][k][l][p];
                        int s = 1;
                        if ( p < i ) s = -s;
                        if ( p < j ) s = -s;
                        if ( p < k ) s = -s;
                        if ( p < l ) s = -s;
                        dum -= s * u_p[d3bbboff[h2] + ijp*trip_aaa[h2]+klp];
                    }
                    A_p[offset + ij*gems_aa[h]+kl] = dum;
                }
            }
            offset += gems_aa[h] * gems_aa[h];
        }
    }
    // D3aab -> D2aa
    for ( int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for ( int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for ( int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double dum = nb * u_p[d2aaoff[h] + ij*gems_aa[h] + kl];
                for ( int p = 0; p < amo_; p++) {
                    int h2 = SymmetryPair(h,symmetry[p]);
                    int ijp = ibas_aab_sym[h2][i][j][p];
                    int klp = ibas_aab_sym[h2][k][l][p];
                    dum -= u_p[d3aaboff[h2] + ijp*trip_aab[h2]+klp];
                }
                A_p[offset + ij*gems_aa[h]+kl] = dum;
            }
        }
        offset += gems_aa[h] * gems_aa[h];
    }
    // D3bba -> D2bb
    for ( int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for ( int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for ( int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double dum = na * u_p[d2bboff[h] + ij*gems_aa[h] + kl];
                for ( int p = 0; p < amo_; p++) {
                    int h2 = SymmetryPair(h,symmetry[p]);
                    int ijp = ibas_aab_sym[h2][i][j][p];
                    int klp = ibas_aab_sym[h2][k][l][p];
                    dum -= u_p[d3bbaoff[h2] + ijp*trip_aab[h2]+klp];
                }
                A_p[offset + ij*gems_aa[h]+kl] = dum;
            }
        }
        offset += gems_aa[h] * gems_aa[h];
    }
    if ( na > 1 ) {
        // D3aab -> D2ab
        for ( int h = 0; h < nirrep_; h++) {
            #pragma omp parallel for schedule (static)
            for ( int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                for ( int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    double dum = (na - 1.0) * u_p[d2aboff[h] + ij*gems_ab[h] + kl];
                    for ( int p = 0; p < amo_; p++) {
                        if ( i == p) continue;
                        if ( k == p) continue;
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijp = ibas_aab_sym[h2][i][p][j];
                        int klp = ibas_aab_sym[h2][k][p][l];
                        int s = 1;
                        if ( p < i ) s = -s;
                        if ( p < k ) s = -s;
                        dum -= s * u_p[d3aaboff[h2] + ijp*trip_aab[h2]+klp];
                    }
                    A_p[offset + ij*gems_ab[h]+kl] = dum;
                }
            }
            offset += gems_ab[h] * gems_ab[h];
        }
    }
    if ( nb > 1 ) {
        // D3bba -> D2ab
        for ( int h = 0; h < nirrep_; h++) {
            #pragma omp parallel for schedule (static)
            for ( int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                for ( int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    double dum = (nb - 1.0) * u_p[d2aboff[h] + ij*gems_ab[h] + kl];
                    for ( int p = 0; p < amo_; p++) {
                        if ( j == p) continue;
                        if ( l == p) continue;
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijp = ibas_aab_sym[h2][j][p][i];
                        int klp = ibas_aab_sym[h2][l][p][k];
                        int s = 1;
                        if ( p < j ) s = -s;
                        if ( p < l ) s = -s;
                        dum -= s * u_p[d3bbaoff[h2] + ijp*trip_aab[h2]+klp];
                    }
                    A_p[offset + ij*gems_ab[h]+kl] = dum;
                }
            }
            offset += gems_ab[h] * gems_ab[h];
        }
    }

    // additional spin constraints for singlets:
    if ( constrain_spin_ && nalpha_ == nbeta_ ) {
        // D3aab = D3bba
        for ( int h = 0; h < nirrep_; h++) {
            C_DCOPY(trip_aab[h]*trip_aab[h],u_p + d3aaboff[h],1,A_p + offset,1);
            C_DAXPY(trip_aab[h]*trip_aab[h],-1.0,u_p + d3bbaoff[h],1,A_p + offset,1);
            offset += trip_aab[h]*trip_aab[h];
        }
        // D3aaa <- D3aab
        for ( int h = 0; h < nirrep_; h++) {
            C_DCOPY(trip_aaa[h]*trip_aaa[h],u_p + d3aaaoff[h],1,A_p + offset,1);
            for (int pqr = 0; pqr < trip_aaa[h]; pqr++) {
                int p = bas_aaa_sym[h][pqr][0];
                int q = bas_aaa_sym[h][pqr][1];
                int r = bas_aaa_sym[h][pqr][2];
                int pqr_b = ibas_aab_sym[h][p][q][r];
                int prq_b = ibas_aab_sym[h][p][r][q];
                int qrp_b = ibas_aab_sym[h][q][r][p];
                for (int stu = 0; stu < trip_aaa[h]; stu++) {
                    int s = bas_aaa_sym[h][stu][0];
                    int t = bas_aaa_sym[h][stu][1];
                    int u = bas_aaa_sym[h][stu][2];
                    int stu_b = ibas_aab_sym[h][s][t][u];
                    int sut_b = ibas_aab_sym[h][s][u][t];
                    int tus_b = ibas_aab_sym[h][t][u][s];
                    A_p[offset + pqr*trip_aaa[h] + stu] -= 1.0/3.0 * u_p[d3aaboff[h] + pqr_b * trip_aab[h] + stu_b];
                    A_p[offset + pqr*trip_aaa[h] + stu] += 1.0/3.0 * u_p[d3aaboff[h] + pqr_b * trip_aab[h] + sut_b];
                    A_p[offset + pqr*trip_aaa[h] + stu] -= 1.0/3.0 * u_p[d3aaboff[h] + pqr_b * trip_aab[h] + tus_b];

                    A_p[offset + pqr*trip_aaa[h] + stu] += 1.0/3.0 * u_p[d3aaboff[h] + prq_b * trip_aab[h] + stu_b];
                    A_p[offset + pqr*trip_aaa[h] + stu] -= 1.0/3.0 * u_p[d3aaboff[h] + prq_b * trip_aab[h] + sut_b];
                    A_p[offset + pqr*trip_aaa[h] + stu] += 1.0/3.0 * u_p[d3aaboff[h] + prq_b * trip_aab[h] + tus_b];

                    A_p[offset + pqr*trip_aaa[h] + stu] -= 1.0/3.0 * u_p[d3aaboff[h] + qrp_b * trip_aab[h] + stu_b];
                    A_p[offset + pqr*trip_aaa[h] + stu] += 1.0/3.0 * u_p[d3aaboff[h] + qrp_b * trip_aab[h] + sut_b];
                    A_p[offset + pqr*trip_aaa[h] + stu] -= 1.0/3.0 * u_p[d3aaboff[h] + qrp_b * trip_aab[h] + tus_b];
                }
            }
            offset += trip_aaa[h]*trip_aaa[h];
        }
        // D3bbb <- D3bba
        for ( int h = 0; h < nirrep_; h++) {
            C_DCOPY(trip_aaa[h]*trip_aaa[h],u_p + d3bbboff[h],1,A_p + offset,1);
            for (int pqr = 0; pqr < trip_aaa[h]; pqr++) {
                int p = bas_aaa_sym[h][pqr][0];
                int q = bas_aaa_sym[h][pqr][1];
                int r = bas_aaa_sym[h][pqr][2];
                int pqr_b = ibas_aab_sym[h][p][q][r];
                int prq_b = ibas_aab_sym[h][p][r][q];
                int qrp_b = ibas_aab_sym[h][q][r][p];
                for (int stu = 0; stu < trip_aaa[h]; stu++) {
                    int s = bas_aaa_sym[h][stu][0];
                    int t = bas_aaa_sym[h][stu][1];
                    int u = bas_aaa_sym[h][stu][2];
                    int stu_b = ibas_aab_sym[h][s][t][u];
                    int sut_b = ibas_aab_sym[h][s][u][t];
                    int tus_b = ibas_aab_sym[h][t][u][s];
                    A_p[offset + pqr*trip_aaa[h] + stu] -= 1.0/3.0 * u_p[d3bbaoff[h] + pqr_b * trip_aab[h] + stu_b];
                    A_p[offset + pqr*trip_aaa[h] + stu] += 1.0/3.0 * u_p[d3bbaoff[h] + pqr_b * trip_aab[h] + sut_b];
                    A_p[offset + pqr*trip_aaa[h] + stu] -= 1.0/3.0 * u_p[d3bbaoff[h] + pqr_b * trip_aab[h] + tus_b];

                    A_p[offset + pqr*trip_aaa[h] + stu] += 1.0/3.0 * u_p[d3bbaoff[h] + prq_b * trip_aab[h] + stu_b];
                    A_p[offset + pqr*trip_aaa[h] + stu] -= 1.0/3.0 * u_p[d3bbaoff[h] + prq_b * trip_aab[h] + sut_b];
                    A_p[offset + pqr*trip_aaa[h] + stu] += 1.0/3.0 * u_p[d3bbaoff[h] + prq_b * trip_aab[h] + tus_b];

                    A_p[offset + pqr*trip_aaa[h] + stu] -= 1.0/3.0 * u_p[d3bbaoff[h] + qrp_b * trip_aab[h] + stu_b];
                    A_p[offset + pqr*trip_aaa[h] + stu] += 1.0/3.0 * u_p[d3bbaoff[h] + qrp_b * trip_aab[h] + sut_b];
                    A_p[offset + pqr*trip_aaa[h] + stu] -= 1.0/3.0 * u_p[d3bbaoff[h] + qrp_b * trip_aab[h] + tus_b];
                }
            }
            offset += trip_aaa[h]*trip_aaa[h];
        }
    }

}

// D3 portion of A^T.y 
void v2RDMSolver::D3_constraints_ATu(SharedVector A,SharedVector u){

    double * A_p = A->pointer();
    double * u_p = u->pointer();

    int na = nalpha_ - nrstc_ - nfrzc_;
    int nb = nbeta_ - nrstc_ - nfrzc_;

    if ( na > 2 ) {
        // D3aaa -> D2aa
        for ( int h = 0; h < nirrep_; h++) {
            for ( int ij = 0; ij < gems_aa[h]; ij++) {
                int i = bas_aa_sym[h][ij][0];
                int j = bas_aa_sym[h][ij][1];
                for ( int kl = 0; kl < gems_aa[h]; kl++) {
                    int k = bas_aa_sym[h][kl][0];
                    int l = bas_aa_sym[h][kl][1];
                    double dum = u_p[offset + ij*gems_aa[h] + kl];
                    A_p[d2aaoff[h] + ij*gems_aa[h] + kl] += (na - 2.0) * dum;
                    for ( int p = 0; p < amo_; p++) {
                        if ( i == p || j == p ) continue;
                        if ( k == p || l == p ) continue;
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijp = ibas_aaa_sym[h2][i][j][p];
                        int klp = ibas_aaa_sym[h2][k][l][p];
                        int s = 1;
                        if ( p < i ) s = -s;
                        if ( p < j ) s = -s;
                        if ( p < k ) s = -s;
                        if ( p < l ) s = -s;
                        A_p[d3aaaoff[h2] + ijp*trip_aaa[h2]+klp] -= s * dum;
                    }
                }
            }
            offset += gems_aa[h] * gems_aa[h];
        }
    }
    if ( nb > 2 ) {
        // D3bbb -> D2bb
        for ( int h = 0; h < nirrep_; h++) {
            for ( int ij = 0; ij < gems_aa[h]; ij++) {
                int i = bas_aa_sym[h][ij][0];
                int j = bas_aa_sym[h][ij][1];
                for ( int kl = 0; kl < gems_aa[h]; kl++) {
                    int k = bas_aa_sym[h][kl][0];
                    int l = bas_aa_sym[h][kl][1];
                    double dum = u_p[offset + ij*gems_aa[h] + kl];
                    A_p[d2bboff[h] + ij*gems_aa[h] + kl] += (nb - 2.0) * dum;
                    for ( int p = 0; p < amo_; p++) {
                        if ( i == p || j == p ) continue;
                        if ( k == p || l == p ) continue;
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijp = ibas_aaa_sym[h2][i][j][p];
                        int klp = ibas_aaa_sym[h2][k][l][p];
                        int s = 1;
                        if ( p < i ) s = -s;
                        if ( p < j ) s = -s;
                        if ( p < k ) s = -s;
                        if ( p < l ) s = -s;
                        A_p[d3bbboff[h2] + ijp*trip_aaa[h2]+klp] -= s * dum;
                    }
                }
            }
            offset += gems_aa[h] * gems_aa[h];
        }
    }
    // D3aab -> D2aa
    for ( int h = 0; h < nirrep_; h++) {
        for ( int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for ( int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double dum = u_p[offset + ij*gems_aa[h] + kl];
                A_p[d2aaoff[h] + ij*gems_aa[h] + kl] += nb * dum;
                for ( int p = 0; p < amo_; p++) {
                    int h2 = SymmetryPair(h,symmetry[p]);
                    int ijp = ibas_aab_sym[h2][i][j][p];
                    int klp = ibas_aab_sym[h2][k][l][p];
                    A_p[d3aaboff[h2] + ijp*trip_aab[h2]+klp] -= dum;
                }
            }
        }
        offset += gems_aa[h] * gems_aa[h];
    }
    // D3bba -> D2bb
    for ( int h = 0; h < nirrep_; h++) {
        for ( int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for ( int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double dum = u_p[offset + ij*gems_aa[h] + kl];
                A_p[d2bboff[h] + ij*gems_aa[h] + kl] += na * dum;
                for ( int p = 0; p < amo_; p++) {
                    int h2 = SymmetryPair(h,symmetry[p]);
                    int ijp = ibas_aab_sym[h2][i][j][p];
                    int klp = ibas_aab_sym[h2][k][l][p];
                    A_p[d3bbaoff[h2] + ijp*trip_aab[h2]+klp] -= dum;
                }
            }
        }
        offset += gems_aa[h] * gems_aa[h];
    }
    if ( na > 1 ) {
        // D3aab -> D2ab
        for ( int h = 0; h < nirrep_; h++) {
            for ( int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                for ( int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    double dum = u_p[offset + ij*gems_ab[h] + kl];
                    A_p[d2aboff[h] + ij*gems_ab[h] + kl] += (na - 1.0) * dum;
                    for ( int p = 0; p < amo_; p++) {
                        if ( i == p) continue;
                        if ( k == p) continue;
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijp = ibas_aab_sym[h2][i][p][j];
                        int klp = ibas_aab_sym[h2][k][p][l];
                        int s = 1;
                        if ( p < i ) s = -s;
                        if ( p < k ) s = -s;
                        A_p[d3aaboff[h2] + ijp*trip_aab[h2]+klp] -= s * dum;
                    }
                }
            }
            offset += gems_ab[h] * gems_ab[h];
        }
    }
    if ( nb > 1 ) {
        // D3bba -> D2ab
        for ( int h = 0; h < nirrep_; h++) {
            for ( int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                for ( int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    double dum = u_p[offset + ij*gems_ab[h] + kl];
                    A_p[d2aboff[h] + ij*gems_ab[h] + kl] += (nb - 1.0) * dum;
                    for ( int p = 0; p < amo_; p++) {
                        if ( j == p) continue;
                        if ( l == p) continue;
                        int h2 = SymmetryPair(h,symmetry[p]);
                        int ijp = ibas_aab_sym[h2][j][p][i];
                        int klp = ibas_aab_sym[h2][l][p][k];
                        int s = 1;
                        if ( p < j ) s = -s;
                        if ( p < l ) s = -s;
                        A_p[d3bbaoff[h2] + ijp*trip_aab[h2]+klp] -= s * dum;
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
            C_DAXPY(trip_aab[h]*trip_aab[h], 1.0,u_p + offset,1,A_p + d3aaboff[h],1);
            C_DAXPY(trip_aab[h]*trip_aab[h],-1.0,u_p + offset,1,A_p + d3bbaoff[h],1);
            offset += trip_aab[h]*trip_aab[h];
        }
        // D3aaa <- D3aab
        for ( int h = 0; h < nirrep_; h++) {
            C_DAXPY(trip_aaa[h]*trip_aaa[h],1.0,u_p + offset,1,A_p+d3aaaoff[h],1);
            for (int pqr = 0; pqr < trip_aaa[h]; pqr++) {
                int p = bas_aaa_sym[h][pqr][0];
                int q = bas_aaa_sym[h][pqr][1];
                int r = bas_aaa_sym[h][pqr][2];
                int pqr_b = ibas_aab_sym[h][p][q][r];
                int prq_b = ibas_aab_sym[h][p][r][q];
                int qrp_b = ibas_aab_sym[h][q][r][p];
                for (int stu = 0; stu < trip_aaa[h]; stu++) {
                    int s = bas_aaa_sym[h][stu][0];
                    int t = bas_aaa_sym[h][stu][1];
                    int u = bas_aaa_sym[h][stu][2];
                    int stu_b = ibas_aab_sym[h][s][t][u];
                    int sut_b = ibas_aab_sym[h][s][u][t];
                    int tus_b = ibas_aab_sym[h][t][u][s];
                    A_p[d3aaboff[h] + pqr_b * trip_aab[h] + stu_b] -= 1.0/3.0 * u_p[offset + pqr*trip_aaa[h] + stu];
                    A_p[d3aaboff[h] + pqr_b * trip_aab[h] + sut_b] += 1.0/3.0 * u_p[offset + pqr*trip_aaa[h] + stu];
                    A_p[d3aaboff[h] + pqr_b * trip_aab[h] + tus_b] -= 1.0/3.0 * u_p[offset + pqr*trip_aaa[h] + stu];
                                                                                                                   
                    A_p[d3aaboff[h] + prq_b * trip_aab[h] + stu_b] += 1.0/3.0 * u_p[offset + pqr*trip_aaa[h] + stu];
                    A_p[d3aaboff[h] + prq_b * trip_aab[h] + sut_b] -= 1.0/3.0 * u_p[offset + pqr*trip_aaa[h] + stu];
                    A_p[d3aaboff[h] + prq_b * trip_aab[h] + tus_b] += 1.0/3.0 * u_p[offset + pqr*trip_aaa[h] + stu];
                                                                                                                   
                    A_p[d3aaboff[h] + qrp_b * trip_aab[h] + stu_b] -= 1.0/3.0 * u_p[offset + pqr*trip_aaa[h] + stu];
                    A_p[d3aaboff[h] + qrp_b * trip_aab[h] + sut_b] += 1.0/3.0 * u_p[offset + pqr*trip_aaa[h] + stu];
                    A_p[d3aaboff[h] + qrp_b * trip_aab[h] + tus_b] -= 1.0/3.0 * u_p[offset + pqr*trip_aaa[h] + stu];
                }
            }
            offset += trip_aaa[h]*trip_aaa[h];
        }
        // D3bbb <- D3bba
        for ( int h = 0; h < nirrep_; h++) {
            C_DAXPY(trip_aaa[h]*trip_aaa[h],1.0,u_p + offset,1,A_p+d3bbboff[h],1);
            for (int pqr = 0; pqr < trip_aaa[h]; pqr++) {
                int p = bas_aaa_sym[h][pqr][0];
                int q = bas_aaa_sym[h][pqr][1];
                int r = bas_aaa_sym[h][pqr][2];
                int pqr_b = ibas_aab_sym[h][p][q][r];
                int prq_b = ibas_aab_sym[h][p][r][q];
                int qrp_b = ibas_aab_sym[h][q][r][p];
                for (int stu = 0; stu < trip_aaa[h]; stu++) {
                    int s = bas_aaa_sym[h][stu][0];
                    int t = bas_aaa_sym[h][stu][1];
                    int u = bas_aaa_sym[h][stu][2];
                    int stu_b = ibas_aab_sym[h][s][t][u];
                    int sut_b = ibas_aab_sym[h][s][u][t];
                    int tus_b = ibas_aab_sym[h][t][u][s];
                    A_p[d3bbaoff[h] + pqr_b * trip_aab[h] + stu_b] -= 1.0/3.0 * u_p[offset + pqr*trip_aaa[h] + stu];
                    A_p[d3bbaoff[h] + pqr_b * trip_aab[h] + sut_b] += 1.0/3.0 * u_p[offset + pqr*trip_aaa[h] + stu];
                    A_p[d3bbaoff[h] + pqr_b * trip_aab[h] + tus_b] -= 1.0/3.0 * u_p[offset + pqr*trip_aaa[h] + stu];
                                                                                                                   
                    A_p[d3bbaoff[h] + prq_b * trip_aab[h] + stu_b] += 1.0/3.0 * u_p[offset + pqr*trip_aaa[h] + stu];
                    A_p[d3bbaoff[h] + prq_b * trip_aab[h] + sut_b] -= 1.0/3.0 * u_p[offset + pqr*trip_aaa[h] + stu];
                    A_p[d3bbaoff[h] + prq_b * trip_aab[h] + tus_b] += 1.0/3.0 * u_p[offset + pqr*trip_aaa[h] + stu];
                                                                                                                   
                    A_p[d3bbaoff[h] + qrp_b * trip_aab[h] + stu_b] -= 1.0/3.0 * u_p[offset + pqr*trip_aaa[h] + stu];
                    A_p[d3bbaoff[h] + qrp_b * trip_aab[h] + sut_b] += 1.0/3.0 * u_p[offset + pqr*trip_aaa[h] + stu];
                    A_p[d3bbaoff[h] + qrp_b * trip_aab[h] + tus_b] -= 1.0/3.0 * u_p[offset + pqr*trip_aaa[h] + stu];
                }
            }
            offset += trip_aaa[h]*trip_aaa[h];
        }
    }

}

} // end namespaces
