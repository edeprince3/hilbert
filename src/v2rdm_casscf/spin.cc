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


// Spin portion of A^T.y 
void v2RDMSolver::Spin_constraints_ATu(double* A,double* u){

    //if ( !constrain_g2_ ) {
        // spin 
        for (int i = 0; i < amo_; i++){
            for (int j = 0; j < amo_; j++){
                int h = SymmetryPair(symmetry[i],symmetry[j]);
                if ( gems_ab[h] == 0 ) continue;
                int ij = ibas_ab_sym[h][i][j];
                int ji = ibas_ab_sym[h][j][i];
                A[d2aboff[h] + ij*gems_ab[h]+ji] += u[offset];
            }
        }
        offset++;
    //}else {
    //    // spin (in terms of G2)
    //    double dum = u[offset];
    //    for (int h = 0; h < nirrep_; h++) {
    //        for (int i = 0; i < amopi_[h]; i++) {
    //            A[d1aoff[h] + i * amopi_[h] + i] += 0.5 * dum;
    //            A[d1boff[h] + i * amopi_[h] + i] -= 0.5 * dum;
    //        }
    //    }
    //    for (int i = 0; i < amo_; i++) {
    //        int ii = ibas_ab_sym[0][i][i];
    //        for (int j = 0; j < amo_; j++) {
    //            int jj = ibas_ab_sym[0][j][j];

    //            A[g2aaoff[0] + (ii             ) * 2*gems_ab[0] + jj             ] += 0.25 * dum;
    //            A[g2aaoff[0] + (ii             ) * 2*gems_ab[0] + jj + gems_ab[0]] -= 0.25 * dum;
    //            A[g2aaoff[0] + (ii + gems_ab[0]) * 2*gems_ab[0] + jj             ] -= 0.25 * dum;
    //            A[g2aaoff[0] + (ii + gems_ab[0]) * 2*gems_ab[0] + jj + gems_ab[0]] += 0.25 * dum;

    //            A[g2baoff[0] + ii * gems_ab[0] + jj] += dum;
    //        }
    //    }
    //    offset++;
    //}

    // additional spin constraints for singlets:
    if ( nalpha_ == nbeta_ ) {
        // D1a = D1b
        for ( int h = 0; h < nirrep_; h++) {
            C_DAXPY(amopi_[h]*amopi_[h], 1.0, u + offset, 1, A + d1aoff[h],1);
            C_DAXPY(amopi_[h]*amopi_[h],-1.0, u + offset, 1, A + d1boff[h],1);
            offset += amopi_[h]*amopi_[h];
        }
        // D2aa = D2bb
        for ( int h = 0; h < nirrep_; h++) {
            C_DAXPY(gems_aa[h]*gems_aa[h], 1.0, u + offset, 1, A + d2aaoff[h],1);
            C_DAXPY(gems_aa[h]*gems_aa[h],-1.0, u + offset, 1, A + d2bboff[h],1);
            offset += gems_aa[h]*gems_aa[h];
        }
        // D2aa[pq][rs] = 1/2(D2ab[pq][rs] - D2ab[pq][sr] - D2ab[qp][rs] + D2ab[qp][sr])
        for ( int h = 0; h < nirrep_; h++) {
            C_DAXPY(gems_aa[h]*gems_aa[h],1.0,u + offset,1,A + d2aaoff[h],1);
            for (int ij = 0; ij < gems_aa[h]; ij++) {
                int i = bas_aa_sym[h][ij][0]; 
                int j = bas_aa_sym[h][ij][1];
                int ijb = ibas_ab_sym[h][i][j];
                int jib = ibas_ab_sym[h][j][i];
                for (int kl = 0; kl < gems_aa[h]; kl++) {
                    int k = bas_aa_sym[h][kl][0]; 
                    int l = bas_aa_sym[h][kl][1];
                    int klb = ibas_ab_sym[h][k][l];
                    int lkb = ibas_ab_sym[h][l][k];
                    A[d2aboff[h] + ijb*gems_ab[h] + klb] -= 0.5 * u[offset + ij*gems_aa[h] + kl];
                    A[d2aboff[h] + jib*gems_ab[h] + klb] += 0.5 * u[offset + ij*gems_aa[h] + kl];
                    A[d2aboff[h] + ijb*gems_ab[h] + lkb] += 0.5 * u[offset + ij*gems_aa[h] + kl];
                    A[d2aboff[h] + jib*gems_ab[h] + lkb] -= 0.5 * u[offset + ij*gems_aa[h] + kl];
                }   
            }   
            offset += gems_aa[h]*gems_aa[h];
        }   
        // D2bb[pq][rs] = 1/2(D2ab[pq][rs] - D2ab[pq][sr] - D2ab[qp][rs] + D2ab[qp][sr])
        for ( int h = 0; h < nirrep_; h++) {
            C_DAXPY(gems_aa[h]*gems_aa[h],1.0,u + offset,1,A + d2bboff[h],1);
            for (int ij = 0; ij < gems_aa[h]; ij++) {
                int i = bas_aa_sym[h][ij][0];
                int j = bas_aa_sym[h][ij][1];
                int ijb = ibas_ab_sym[h][i][j];
                int jib = ibas_ab_sym[h][j][i];
                for (int kl = 0; kl < gems_aa[h]; kl++) {
                    int k = bas_aa_sym[h][kl][0];
                    int l = bas_aa_sym[h][kl][1];
                    int klb = ibas_ab_sym[h][k][l];
                    int lkb = ibas_ab_sym[h][l][k];
                    A[d2aboff[h] + ijb*gems_ab[h] + klb] -= 0.5 * u[offset + ij*gems_aa[h] + kl];
                    A[d2aboff[h] + jib*gems_ab[h] + klb] += 0.5 * u[offset + ij*gems_aa[h] + kl];
                    A[d2aboff[h] + ijb*gems_ab[h] + lkb] += 0.5 * u[offset + ij*gems_aa[h] + kl];
                    A[d2aboff[h] + jib*gems_ab[h] + lkb] -= 0.5 * u[offset + ij*gems_aa[h] + kl];
                }
            }
            offset += gems_aa[h]*gems_aa[h];
        }
        // D200 = 1/(2 sqrt(1+dpq)sqrt(1+drs)) ( D2ab[pq][rs] + D2ab[pq][sr] + D2ab[qp][rs] + D2ab[qp][sr] )
        for ( int h = 0; h < nirrep_; h++) {
            C_DAXPY(gems_ab[h]*gems_ab[h],1.0,u + offset,1,A + d200off[h],1);
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                double dij = ( i == j ) ? sqrt(2.0) : 1.0;
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    double dkl = ( k == l ) ? sqrt(2.0) : 1.0;
                    A[d2aboff[h] + ij*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u[offset + ij*gems_ab[h] + kl];
                    A[d2aboff[h] + ji*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u[offset + ij*gems_ab[h] + kl];
                    A[d2aboff[h] + ij*gems_ab[h] + lk] -= 0.5 / ( dij * dkl ) * u[offset + ij*gems_ab[h] + kl];
                    A[d2aboff[h] + ji*gems_ab[h] + lk] -= 0.5 / ( dij * dkl ) * u[offset + ij*gems_ab[h] + kl];
                }
            }
            offset += gems_ab[h]*gems_ab[h];
        }
        // D2ab[pq][rs] = D2ab[qp][sr]
        for ( int h = 0; h < nirrep_; h++) {
            C_DAXPY(gems_ab[h]*gems_ab[h],1.0,u + offset,1,A + d2aboff[h],1);
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    A[d2aboff[h] + ji*gems_ab[h] + lk] -= u[offset + ij*gems_ab[h] + kl];
                }
            }
            offset += gems_ab[h]*gems_ab[h];
        }
    }else { // nonsinglets ... big block

        for ( int h = 0; h < nirrep_; h++) {
            // D200
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                double dij = ( i == j ) ? sqrt(2.0) : 1.0;
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    double dkl = ( k == l ) ? sqrt(2.0) : 1.0;
                    A[d200off[h] + ij*2*gems_ab[h] + kl] += u[offset + ij*2*gems_ab[h] + kl];
                    A[d2aboff[h] + ij*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u[offset + ij*2*gems_ab[h] + kl];
                    A[d2aboff[h] + ji*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u[offset + ij*2*gems_ab[h] + kl];
                    A[d2aboff[h] + ij*gems_ab[h] + lk] -= 0.5 / ( dij * dkl ) * u[offset + ij*2*gems_ab[h] + kl];
                    A[d2aboff[h] + ji*gems_ab[h] + lk] -= 0.5 / ( dij * dkl ) * u[offset + ij*2*gems_ab[h] + kl];
                }
            }
            // D201
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                double dij = ( i == j ) ? sqrt(2.0) : 1.0;
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    A[d200off[h] + (ij)*2*gems_ab[h] + (kl+gems_ab[h])] += u[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])];
                    A[d2aboff[h] + ij*gems_ab[h] + kl] -= 0.5 / dij * u[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])];
                    A[d2aboff[h] + ij*gems_ab[h] + lk] += 0.5 / dij * u[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])];
                    A[d2aboff[h] + ji*gems_ab[h] + kl] -= 0.5 / dij * u[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])];
                    A[d2aboff[h] + ji*gems_ab[h] + lk] += 0.5 / dij * u[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])];
                }
            }
            // D210
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    double dkl = ( k == l ) ? sqrt(2.0) : 1.0;
                    A[d200off[h] + (ij+gems_ab[h])*2*gems_ab[h] + (kl)] += u[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)];
                    A[d2aboff[h] + ij*gems_ab[h] + kl] -= 0.5 / dkl * u[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)];
                    A[d2aboff[h] + ij*gems_ab[h] + lk] -= 0.5 / dkl * u[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)];
                    A[d2aboff[h] + ji*gems_ab[h] + kl] += 0.5 / dkl * u[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)];
                    A[d2aboff[h] + ji*gems_ab[h] + lk] += 0.5 / dkl * u[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)];
                }
            }
            // D211
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    A[d200off[h] + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])] += u[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])];
                    A[d2aboff[h] + ij*gems_ab[h] + kl] -= 0.5 * u[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])];
                    A[d2aboff[h] + ji*gems_ab[h] + kl] += 0.5 * u[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])];
                    A[d2aboff[h] + ij*gems_ab[h] + lk] += 0.5 * u[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])];
                    A[d2aboff[h] + ji*gems_ab[h] + lk] -= 0.5 * u[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])];
                }
            }
            offset += 4*gems_ab[h]*gems_ab[h];
        }
    }

    if ( constrain_g2_ ) {
        // maximal spin constraint (sum_i G2(kl,ii) = 0)
        for (int kl = 0; kl < gems_ab[0]; kl++) {
            double dum = u[offset + kl];
            for (int i = 0; i < amo_; i++) {
                int ii = ibas_ab_sym[0][i][i];
                A[g2baoff[0] + kl*gems_ab[0]+ii] += dum;
            }
        }
        offset += gems_ab[0];
        for (int kl = 0; kl < gems_ab[0]; kl++) {
            double dum = u[offset + kl];
            for (int i = 0; i < amo_; i++) {
                int ii = ibas_ab_sym[0][i][i];
                A[g2baoff[0] + ii*gems_ab[0]+kl] += dum;
            }
        }
        offset += gems_ab[0];
    }

}

// D2 portion of A.x (and D1/Q1)
void v2RDMSolver::Spin_constraints_Au(double* A,double* u){

    //if ( !constrain_g2_ ) {
        // spin
        double s2 = 0.0;
        for (int i = 0; i < amo_; i++){
            for (int j = 0; j < amo_; j++){
                int h = SymmetryPair(symmetry[i],symmetry[j]);
                if ( gems_ab[h] == 0 ) continue;
                int ij = ibas_ab_sym[h][i][j];
                int ji = ibas_ab_sym[h][j][i];
                s2 += u[d2aboff[h] + ij*gems_ab[h]+ji];
            }
        }
        A[offset] = s2;
        offset++;
    //}else {
    //    // spin (in terms of G2)
    //    double s2 = 0.0;
    //    for (int h = 0; h < nirrep_; h++) {
    //        for (int i = 0; i < amopi_[h]; i++) {
    //            s2 += 0.5 * u[d1aoff[h] + i * amopi_[h] + i];
    //            s2 -= 0.5 * u[d1boff[h] + i * amopi_[h] + i];
    //        }
    //    }
    //    for (int i = 0; i < amo_; i++) {
    //        int ii = ibas_ab_sym[0][i][i];
    //        for (int j = 0; j < amo_; j++) {
    //            int jj = ibas_ab_sym[0][j][j];

    //            s2 += 0.25 * u[g2aaoff[0] + (ii             ) * 2*gems_ab[0] + jj             ];
    //            s2 -= 0.25 * u[g2aaoff[0] + (ii             ) * 2*gems_ab[0] + jj + gems_ab[0]];
    //            s2 -= 0.25 * u[g2aaoff[0] + (ii + gems_ab[0]) * 2*gems_ab[0] + jj             ];
    //            s2 += 0.25 * u[g2aaoff[0] + (ii + gems_ab[0]) * 2*gems_ab[0] + jj + gems_ab[0]];

    //            s2 += u[g2baoff[0] + ii * gems_ab[0] + jj];
    //        }
    //    }
    //    A[offset] = s2;
    //    offset++;
    //}

    // additional spin constraints for singlets:
    if ( nalpha_ == nbeta_ ) {
        // D1a = D1b
        for ( int h = 0; h < nirrep_; h++) {
            C_DCOPY(amopi_[h]*amopi_[h],     u + d1aoff[h],1,A + offset,1);
            C_DAXPY(amopi_[h]*amopi_[h],-1.0,u + d1boff[h],1,A + offset,1);
            offset += amopi_[h]*amopi_[h]; 
        }
        // D2aa = D2bb
        for ( int h = 0; h < nirrep_; h++) {
            C_DCOPY(gems_aa[h]*gems_aa[h],     u + d2aaoff[h],1,A + offset,1);
            C_DAXPY(gems_aa[h]*gems_aa[h],-1.0,u + d2bboff[h],1,A + offset,1);
            offset += gems_aa[h]*gems_aa[h];
        }
        // D2aa[pq][rs] = 1/2(D2ab[pq][rs] - D2ab[pq][sr] - D2ab[qp][rs] + D2ab[qp][sr])
        for ( int h = 0; h < nirrep_; h++) {
            C_DCOPY(gems_aa[h]*gems_aa[h],u + d2aaoff[h],1,A + offset,1);
            for (int ij = 0; ij < gems_aa[h]; ij++) {
                int i = bas_aa_sym[h][ij][0];
                int j = bas_aa_sym[h][ij][1];
                int ijb = ibas_ab_sym[h][i][j];
                int jib = ibas_ab_sym[h][j][i];
                for (int kl = 0; kl < gems_aa[h]; kl++) {
                    int k = bas_aa_sym[h][kl][0];
                    int l = bas_aa_sym[h][kl][1];
                    int klb = ibas_ab_sym[h][k][l];
                    int lkb = ibas_ab_sym[h][l][k];
                    A[offset + ij*gems_aa[h] + kl] -= 0.5 * u[d2aboff[h] + ijb*gems_ab[h] + klb];
                    A[offset + ij*gems_aa[h] + kl] += 0.5 * u[d2aboff[h] + jib*gems_ab[h] + klb];
                    A[offset + ij*gems_aa[h] + kl] += 0.5 * u[d2aboff[h] + ijb*gems_ab[h] + lkb];
                    A[offset + ij*gems_aa[h] + kl] -= 0.5 * u[d2aboff[h] + jib*gems_ab[h] + lkb];
                }
            }
            offset += gems_aa[h]*gems_aa[h];
        }
        // D2bb[pq][rs] = 1/2(D2ab[pq][rs] - D2ab[pq][sr] - D2ab[qp][rs] + D2ab[qp][sr])
        for ( int h = 0; h < nirrep_; h++) {
            C_DCOPY(gems_aa[h]*gems_aa[h],u + d2bboff[h],1,A + offset,1);
            for (int ij = 0; ij < gems_aa[h]; ij++) {
                int i = bas_aa_sym[h][ij][0];
                int j = bas_aa_sym[h][ij][1];
                int ijb = ibas_ab_sym[h][i][j];
                int jib = ibas_ab_sym[h][j][i];
                for (int kl = 0; kl < gems_aa[h]; kl++) {
                    int k = bas_aa_sym[h][kl][0];
                    int l = bas_aa_sym[h][kl][1];
                    int klb = ibas_ab_sym[h][k][l];
                    int lkb = ibas_ab_sym[h][l][k];
                    A[offset + ij*gems_aa[h] + kl] -= 0.5 * u[d2aboff[h] + ijb*gems_ab[h] + klb];
                    A[offset + ij*gems_aa[h] + kl] += 0.5 * u[d2aboff[h] + jib*gems_ab[h] + klb];
                    A[offset + ij*gems_aa[h] + kl] += 0.5 * u[d2aboff[h] + ijb*gems_ab[h] + lkb];
                    A[offset + ij*gems_aa[h] + kl] -= 0.5 * u[d2aboff[h] + jib*gems_ab[h] + lkb];
                }
            }
            offset += gems_aa[h]*gems_aa[h];
        }
        // D200 = 1/(2 sqrt(1+dpq)sqrt(1+drs)) ( D2ab[pq][rs] + D2ab[pq][sr] + D2ab[qp][rs] + D2ab[qp][sr] )
        for ( int h = 0; h < nirrep_; h++) {
            C_DCOPY(gems_ab[h]*gems_ab[h],u + d200off[h],1,A + offset,1);
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                double dij = ( i == j ) ? sqrt(2.0) : 1.0;
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    double dkl = ( k == l ) ? sqrt(2.0) : 1.0;
                    A[offset + ij*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u[d2aboff[h] + ij*gems_ab[h] + kl];
                    A[offset + ij*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u[d2aboff[h] + ji*gems_ab[h] + kl];
                    A[offset + ij*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u[d2aboff[h] + ij*gems_ab[h] + lk];
                    A[offset + ij*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u[d2aboff[h] + ji*gems_ab[h] + lk];
                }
            }
            offset += gems_ab[h]*gems_ab[h];
        }
        // D2ab[pq][rs] = D2ab[qp][sr]
        for ( int h = 0; h < nirrep_; h++) {
            C_DAXPY(gems_ab[h]*gems_ab[h],1.0,u + d2aboff[h],1,A + offset,1);
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0]; 
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0]; 
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    A[offset + ij*gems_ab[h] + kl] -= u[d2aboff[h] + ji*gems_ab[h] + lk];
                }
            }
            offset += gems_ab[h]*gems_ab[h];
        }
    }else { // nonsinglets ... big block

        for ( int h = 0; h < nirrep_; h++) {
            // D200
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                double dij = ( i == j ) ? sqrt(2.0) : 1.0;
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    double dkl = ( k == l ) ? sqrt(2.0) : 1.0;
                    A[offset + ij*2*gems_ab[h] + kl] += u[d200off[h] + ij*2*gems_ab[h] + kl];
                    A[offset + ij*2*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u[d2aboff[h] + ij*gems_ab[h] + kl];
                    A[offset + ij*2*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u[d2aboff[h] + ji*gems_ab[h] + kl];
                    A[offset + ij*2*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u[d2aboff[h] + ij*gems_ab[h] + lk];
                    A[offset + ij*2*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u[d2aboff[h] + ji*gems_ab[h] + lk];
                }
            }
            // D201
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                double dij = ( i == j ) ? sqrt(2.0) : 1.0;
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    A[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])] += u[d200off[h] + (ij)*2*gems_ab[h] + (kl+gems_ab[h])];
                    A[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])] -= 0.5 / dij * u[d2aboff[h] + ij*gems_ab[h] + kl];
                    A[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])] += 0.5 / dij * u[d2aboff[h] + ij*gems_ab[h] + lk];
                    A[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])] -= 0.5 / dij * u[d2aboff[h] + ji*gems_ab[h] + kl];
                    A[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])] += 0.5 / dij * u[d2aboff[h] + ji*gems_ab[h] + lk];
                }
            }
            // D210
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    double dkl = ( k == l ) ? sqrt(2.0) : 1.0;
                    A[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)] += u[d200off[h] + (ij+gems_ab[h])*2*gems_ab[h] + (kl)];
                    A[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)] -= 0.5 / dkl * u[d2aboff[h] + ij*gems_ab[h] + kl];
                    A[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)] -= 0.5 / dkl * u[d2aboff[h] + ij*gems_ab[h] + lk];
                    A[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)] += 0.5 / dkl * u[d2aboff[h] + ji*gems_ab[h] + kl];
                    A[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)] += 0.5 / dkl * u[d2aboff[h] + ji*gems_ab[h] + lk];
                }
            }
            // D211
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    A[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])] += u[d200off[h] + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])];
                    A[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])] -= 0.5 * u[d2aboff[h] + ij*gems_ab[h] + kl];
                    A[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])] += 0.5 * u[d2aboff[h] + ji*gems_ab[h] + kl];
                    A[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])] += 0.5 * u[d2aboff[h] + ij*gems_ab[h] + lk];
                    A[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])] -= 0.5 * u[d2aboff[h] + ji*gems_ab[h] + lk];
                }
            }
            offset += 4*gems_ab[h]*gems_ab[h];
        }
    }

    if ( constrain_g2_ ) {
        // maximal spin constraint (sum_i G2(kl,ii) = 0)
        for (int kl = 0; kl < gems_ab[0]; kl++) {
            double dum = 0.0; 
            for (int i = 0; i < amo_; i++) {
                int ii = ibas_ab_sym[0][i][i];
                dum += u[g2baoff[0] + kl*gems_ab[0]+ii];
            }
            A[offset + kl] = dum;
        }
        offset += gems_ab[0];
        for (int kl = 0; kl < gems_ab[0]; kl++) {
            double dum = 0.0;
            for (int i = 0; i < amo_; i++) {
                int ii = ibas_ab_sym[0][i][i];
                dum += u[g2baoff[0] + ii*gems_ab[0]+kl];
            }
            A[offset + kl] = dum;
        }
        offset += gems_ab[0];
    }

}

}
