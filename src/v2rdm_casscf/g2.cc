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

void v2RDMSolver::G2_constraints_guess_spin_adapted(SharedVector u){

    double * u_p = u->pointer();

    // G200
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = 0.0;

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum       +=  u_p[d1aoff[h3] + ii*amopi_[h3]+kk] * 0.5; //   D1(i,k) djl
                    dum       +=  u_p[d1boff[h3] + ii*amopi_[h3]+kk] * 0.5; //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );


                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];
                    dum       -=  u_p[d2aaoff[h2] + ild*gems_aa[h2]+kjd] * sil * skj * 0.5; // -D2aa(il,kj)
                    dum       -=  u_p[d2bboff[h2] + ild*gems_aa[h2]+kjd] * sil * skj * 0.5; // -D2bb(il,kj)

                }

                int ild = ibas_ab_sym[h2][i][l];
                int jkd = ibas_ab_sym[h2][j][k];

                dum       +=  u_p[d2aboff[h2] + ild*gems_ab[h2]+jkd] * 0.5; // D2ab(il,jk)

                int lid = ibas_ab_sym[h2][l][i];
                int kjd = ibas_ab_sym[h2][k][j];

                dum       +=  u_p[d2aboff[h2] + lid*gems_ab[h2]+kjd] * 0.5; // D2ab(li,kj)

                u_p[g2soff[h] + ijg*gems_ab[h]+klg] = dum;
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }
    // G210
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = 0.0;

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum       +=  u_p[d1aoff[h3] + ii*amopi_[h3]+kk] * 0.5; //   D1(i,k) djl
                    dum       +=  u_p[d1boff[h3] + ii*amopi_[h3]+kk] * 0.5; //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];
                    dum       -=  u_p[d2aaoff[h2] + ild*gems_aa[h2]+kjd] * sil * skj * 0.5; // -D2aa(il,kj)
                    dum       -=  u_p[d2bboff[h2] + ild*gems_aa[h2]+kjd] * sil * skj * 0.5; // -D2bb(il,kj)
                }

                int ild = ibas_ab_sym[h2][i][l];
                int jkd = ibas_ab_sym[h2][j][k];

                dum       -=  u_p[d2aboff[h2] + ild*gems_ab[h2]+jkd] * 0.5; // D2ab(il,jk)

                int lid = ibas_ab_sym[h2][l][i];
                int kjd = ibas_ab_sym[h2][k][j];

                dum       -=  u_p[d2aboff[h2] + lid*gems_ab[h2]+kjd] * 0.5; // D2ab(li,kj)

                u_p[g2toff[h] + ijg*gems_ab[h]+klg] = dum;
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }
       
    // G211 constraints:
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = 0.0;

                if ( j==l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum   +=  u_p[d1aoff[h3] + ii*amopi_[h3]+kk];      //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                int ild = ibas_ab_sym[h2][i][l];
                int kjd = ibas_ab_sym[h2][k][j];

                dum       -=  u_p[d2aboff[h2] + ild*gems_ab[h2]+kjd];   // - D2ab(il,kj)

                u_p[g2toff_p1[h] + ijg*gems_ab[h]+klg] = dum;
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }

    // G21-1 constraints:
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = 0.0;

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum   +=  u_p[d1boff[h3] + ii*amopi_[h3]+kk];      //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                int ild = ibas_ab_sym[h2][l][i];
                int kjd = ibas_ab_sym[h2][j][k];

                dum       -=  u_p[d2aboff[h2] + ild*gems_ab[h2]+kjd];   // - D2ab(il,kj)

                u_p[g2toff_m1[h] + ijg*gems_ab[h]+klg] = dum;
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }
}
void v2RDMSolver::G2_constraints_Au_spin_adapted(SharedVector A,SharedVector u){

    double * A_p = A->pointer();
    double * u_p = u->pointer();

    // G200
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = -u_p[g2soff[h] + ijg*gems_ab[h]+klg];          // - G2s(ij,kl)

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum       +=  u_p[d1aoff[h3] + ii*amopi_[h3]+kk] * 0.5; //   D1(i,k) djl
                    dum       +=  u_p[d1boff[h3] + ii*amopi_[h3]+kk] * 0.5; //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                //int ils = ibas_00_sym[h2][i][l];
                //int jks = ibas_00_sym[h2][j][k];

                //dum       +=  u_p[d2soff[h2] + INDEX(ils,jks)] * 0.5; //   D2s(li,kj)

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );


                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];
                    dum       -=  u_p[d2aaoff[h2] + ild*gems_aa[h2]+kjd] * sil * skj * 0.5; // -D2aa(il,kj)
                    dum       -=  u_p[d2bboff[h2] + ild*gems_aa[h2]+kjd] * sil * skj * 0.5; // -D2bb(il,kj)

                    //int ilt = ibas_aa_sym[h2][i][l];
                    //int kjt = ibas_aa_sym[h2][k][j];
                    //dum       -=  u_p[d2toff[h2]    + INDEX(ilt,kjt)] * sil * skj * 0.5; //   D210(il,kj)
                    //dum       -=  u_p[d2toff_p1[h2] + INDEX(ilt,kjt)] * sil * skj * 0.5; //   D211(il,kj)
                    //dum       -=  u_p[d2toff_m1[h2] + INDEX(ilt,kjt)] * sil * skj * 0.5; //   D21-1(il,kj)

                }

                int ild = ibas_ab_sym[h2][i][l];
                int jkd = ibas_ab_sym[h2][j][k];

                dum       +=  u_p[d2aboff[h2] + ild*gems_ab[h2]+jkd] * 0.5; // D2ab(il,jk)

                int lid = ibas_ab_sym[h2][l][i];
                int kjd = ibas_ab_sym[h2][k][j];

                dum       +=  u_p[d2aboff[h2] + lid*gems_ab[h2]+kjd] * 0.5; // D2ab(li,kj)

                A_p[offset + ijg*gems_ab[h]+klg] = dum;
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }
    // G210
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = -u_p[g2toff[h] + ijg*gems_ab[h]+klg];          // - G2t(ij,kl)

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum       +=  u_p[d1aoff[h3] + ii*amopi_[h3]+kk] * 0.5; //   D1(i,k) djl
                    dum       +=  u_p[d1boff[h3] + ii*amopi_[h3]+kk] * 0.5; //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                //int ils = ibas_00_sym[h2][i][l];
                //int jks = ibas_00_sym[h2][j][k];

                //dum       -=  u_p[d2soff[h2] + INDEX(ils,jks)] * 0.5; //   D2s(li,kj)

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];
                    dum       -=  u_p[d2aaoff[h2] + ild*gems_aa[h2]+kjd] * sil * skj * 0.5; // -D2aa(il,kj)
                    dum       -=  u_p[d2bboff[h2] + ild*gems_aa[h2]+kjd] * sil * skj * 0.5; // -D2bb(il,kj)

                    //int ilt = ibas_aa_sym[h2][i][l];
                    //int kjt = ibas_aa_sym[h2][k][j];
                    //dum       +=  u_p[d2toff[h2]    + INDEX(ilt,kjt)] * sil * skj * 0.5; //   D210(il,kj)
                    //dum       -=  u_p[d2toff_p1[h2] + INDEX(ilt,kjt)] * sil * skj * 0.5; //   D211(il,kj)
                    //dum       -=  u_p[d2toff_m1[h2] + INDEX(ilt,kjt)] * sil * skj * 0.5; //   D21-1(il,kj)

                }

                int ild = ibas_ab_sym[h2][i][l];
                int jkd = ibas_ab_sym[h2][j][k];

                dum       -=  u_p[d2aboff[h2] + ild*gems_ab[h2]+jkd] * 0.5; // D2ab(il,jk)

                int lid = ibas_ab_sym[h2][l][i];
                int kjd = ibas_ab_sym[h2][k][j];

                dum       -=  u_p[d2aboff[h2] + lid*gems_ab[h2]+kjd] * 0.5; // D2ab(li,kj)

                A_p[offset + ijg*gems_ab[h]+klg] = dum;
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }
       
    // G211 constraints:
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = -u_p[g2toff_p1[h] + ijg*gems_ab[h]+klg];    // - G2ab(ij,kl)

                if ( j==l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum   +=  u_p[d1aoff[h3] + ii*amopi_[h3]+kk];      //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                int ild = ibas_ab_sym[h2][i][l];
                int kjd = ibas_ab_sym[h2][k][j];

                dum       -=  u_p[d2aboff[h2] + ild*gems_ab[h2]+kjd];   // - D2ab(il,kj)

                //int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                //int ils = ibas_00_sym[h2][i][l];
                //int kjs = ibas_00_sym[h2][k][j];
                //dum       -=  u_p[d2soff[h2] + INDEX(ils,kjs)] * 0.5; //   D2s(li,kj)

                //if ( i != l && k != j ) {

                //    int sil = ( i < l ? 1 : -1 );
                //    int skj = ( k < j ? 1 : -1 );

                //    int ilt = ibas_aa_sym[h2][i][l];
                //    int kjt = ibas_aa_sym[h2][k][j];

                //    dum       -=  u_p[d2toff_p1[h2] + INDEX(ilt,kjt)] * sil * skj * 0.5; //   D211(il,kj)

                //}

                A_p[offset + ijg*gems_ab[h]+klg] = dum;
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }
    // G21-1 constraints:
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = -u_p[g2toff_m1[h] + ijg*gems_ab[h]+klg];    // - G2ab(ij,kl)

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum   +=  u_p[d1boff[h3] + ii*amopi_[h3]+kk];      //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                int ild = ibas_ab_sym[h2][l][i];
                int kjd = ibas_ab_sym[h2][j][k];

                dum       -=  u_p[d2aboff[h2] + ild*gems_ab[h2]+kjd];   // - D2ab(il,kj)

                //int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                //int ils = ibas_00_sym[h2][i][l];
                //int kjs = ibas_00_sym[h2][k][j];
                //dum       -=  u_p[d2soff[h2] + INDEX(ils,kjs)] * 0.5; //   D2s(li,kj)

                //if ( i != l && k != j ) {

                //    int sil = ( i < l ? 1 : -1 );
                //    int skj = ( k < j ? 1 : -1 );

                //    int ilt = ibas_aa_sym[h2][i][l];
                //    int kjt = ibas_aa_sym[h2][k][j];

                //    dum       -=  u_p[d2toff_m1[h2] + INDEX(ilt,kjt)] * sil * skj * 0.5; //   D211(il,kj)

                //}

                A_p[offset + ijg*gems_ab[h]+klg] = dum;
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }
    // contraction condition:
    /*int na = nalpha_ - nrstc_ - nfrzc_;
    int nb = nbeta_ - nrstc_ - nfrzc_;
    int ms = (multiplicity_ - 1)/2;
    for (int h = 0; h < nirrep_; h++) {
        for (int klg = 0; klg < gems_ab[h]; klg++) {

            int k = bas_ab_sym[h][klg][0];
            int l = bas_ab_sym[h][klg][1];

            double dum = 0.0;

            for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

                int i = bas_ab_sym[h][ijg][0];
                int j = bas_ab_sym[h][ijg][1];
                if ( i != j ) continue;


                dum += (na+nb)*u_p[g2toff[h] + ijg*gems_ab[h] + klg];
                dum += -2.0 * ms * u_p[g2soff[h] + ijg*gems_ab[h] + klg];

            }

            A_p[offset + k*amo_+l] = dum;

        }
    }
    offset += amo_*amo_;*/
    
    // maximal spin constraint:
    /*for (int h = 0; h < nirrep_; h++) {
        for (int klg = 0; klg < gems_ab[h]; klg++) {

            int k = bas_ab_sym[h][klg][0];
            int l = bas_ab_sym[h][klg][1];

            double maxs = 0.0;

            for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

                int i = bas_ab_sym[h][ijg][0];
                int j = bas_ab_sym[h][ijg][1];
                if ( i != j ) continue;

                maxs += u_p[g2toff_p1[h] + ijg*gems_ab[h] + klg];

            }

printf("%20.12lf\n",maxs);
            A_p[offset + k*amo_+l] = maxs;

        }
    }
    offset += amo_*amo_;
    for (int h = 0; h < nirrep_; h++) {
        for (int klg = 0; klg < gems_ab[h]; klg++) {

            int k = bas_ab_sym[h][klg][0];
            int l = bas_ab_sym[h][klg][1];

            double maxs = 0.0;

            for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

                int i = bas_ab_sym[h][ijg][0];
                int j = bas_ab_sym[h][ijg][1];
                if ( i != j ) continue;

                maxs += u_p[g2toff_p1[h] + klg*gems_ab[h] + ijg];

            }

            A_p[offset + k*amo_+l] = maxs;

        }
    }
    offset += amo_*amo_;*/

}
// G2 portion of A^T.y (spin adapted)
void v2RDMSolver::G2_constraints_ATu_spin_adapted(SharedVector A,SharedVector u){

    double * A_p = A->pointer();
    double * u_p = u->pointer();

    // G200
    for (int h = 0; h < nirrep_; h++) {
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = u_p[offset + ijg*gems_ab[h]+klg];

                A_p[g2soff[h] + ijg*gems_ab[h]+klg] -= dum;

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    A_p[d1aoff[h3] + ii*amopi_[h3]+kk]         += 0.5 * dum;
                    A_p[d1boff[h3] + ii*amopi_[h3]+kk]         += 0.5 * dum;
                }

                //int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                //int ils = ibas_00_sym[h2][i][l];
                //int jks = ibas_00_sym[h2][j][k];
                //A_p[d2soff[h2] + INDEX(ils,jks)] += dum * 0.5;

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];

                    A_p[d2aaoff[h2] + ild*gems_aa[h2]+kjd] -= 0.5 * dum * sil * skj;
                    A_p[d2bboff[h2] + ild*gems_aa[h2]+kjd] -= 0.5 * dum * sil * skj;

                    //int ilt = ibas_aa_sym[h2][i][l];
                    //int kjt = ibas_aa_sym[h2][k][j];
                    //A_p[d2toff[h2]    + INDEX(ilt,kjt)] -= dum * sil * skj * 0.5; // 10
                    //A_p[d2toff_p1[h2] + INDEX(ilt,kjt)] -= dum * sil * skj * 0.5; // 11
                    //A_p[d2toff_m1[h2] + INDEX(ilt,kjt)] -= dum * sil * skj * 0.5; // 1-1
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                int ild = ibas_ab_sym[h2][i][l];
                int jkd = ibas_ab_sym[h2][j][k];

                A_p[d2aboff[h2] + ild*gems_ab[h2]+jkd] += 0.5 * dum;

                int lid = ibas_ab_sym[h2][l][i];
                int kjd = ibas_ab_sym[h2][k][j];

                A_p[d2aboff[h2] + lid*gems_ab[h2]+kjd] += 0.5 * dum;
            }
        }

        offset += gems_ab[h]*gems_ab[h];
    }
    // G210
    for (int h = 0; h < nirrep_; h++) {
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = u_p[offset + ijg*gems_ab[h]+klg];

                A_p[g2toff[h] + ijg*gems_ab[h]+klg] -= dum;

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    A_p[d1aoff[h3] + ii*amopi_[h3]+kk]         += 0.5 * dum;
                    A_p[d1boff[h3] + ii*amopi_[h3]+kk]         += 0.5 * dum;
                }

                //int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                //int ils = ibas_00_sym[h2][i][l];
                //int jks = ibas_00_sym[h2][j][k];

                //A_p[d2soff[h2] + INDEX(ils,jks)] -= dum * 0.5;

                //if ( i != l && k != j ) {

                //    int sil = ( i < l ? 1 : -1 );
                //    int skj = ( k < j ? 1 : -1 );

                //    int ilt = ibas_aa_sym[h2][i][l];
                //    int kjt = ibas_aa_sym[h2][k][j];

                //    A_p[d2toff[h2]    + INDEX(ilt,kjt)] += dum * sil * skj * 0.5; // 10
                //    A_p[d2toff_p1[h2] + INDEX(ilt,kjt)] -= dum * sil * skj * 0.5; // 11
                //    A_p[d2toff_m1[h2] + INDEX(ilt,kjt)] -= dum * sil * skj * 0.5; // 1-1
                //}

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];

                    A_p[d2aaoff[h2] + ild*gems_aa[h2]+kjd] -= 0.5 * dum * sil * skj;
                    A_p[d2bboff[h2] + ild*gems_aa[h2]+kjd] -= 0.5 * dum * sil * skj;
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                int ild = ibas_ab_sym[h2][i][l];
                int jkd = ibas_ab_sym[h2][j][k];

                A_p[d2aboff[h2] + ild*gems_ab[h2]+jkd] -= 0.5 * dum;

                int lid = ibas_ab_sym[h2][l][i];
                int kjd = ibas_ab_sym[h2][k][j];

                A_p[d2aboff[h2] + lid*gems_ab[h2]+kjd] -= 0.5 * dum;
            }
        }

        offset += gems_ab[h]*gems_ab[h];
    }
    // G211 constraints:
    for (int h = 0; h < nirrep_; h++) {
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = u_p[offset + ijg*gems_ab[h]+klg];

                A_p[g2toff_p1[h] + ijg*gems_ab[h]+klg] -= dum;    // - G2ab(ij,kl)

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    A_p[d1aoff[h3] + ii*amopi_[h3]+kk] += dum;      //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                int ild = ibas_ab_sym[h2][i][l];
                int kjd = ibas_ab_sym[h2][k][j];

                A_p[d2aboff[h2] + ild*gems_ab[h2]+kjd] -= dum;   // - D2ab(il,kj)

                //int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                //int ils = ibas_00_sym[h2][i][l];
                //int kjs = ibas_00_sym[h2][k][j];
                //A_p[d2soff[h2] + INDEX(ils,kjs)]             -= dum * 0.5;

                //if ( i != l && k != j ) {

                //    int sil = ( i < l ? 1 : -1 );
                //    int skj = ( k < j ? 1 : -1 );

                //    int ilt = ibas_aa_sym[h2][i][l];
                //    int kjt = ibas_aa_sym[h2][k][j];

                //    A_p[d2toff_p1[h2] + INDEX(ilt,kjt)] -= dum * sil * skj * 0.5;
                //}
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }
    // G21-1 constraints:
    for (int h = 0; h < nirrep_; h++) {
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = u_p[offset + ijg*gems_ab[h]+klg];

                A_p[g2toff_m1[h] + ijg*gems_ab[h]+klg] -= dum;    // - G2ab(ij,kl)

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    A_p[d1boff[h3] + ii*amopi_[h3]+kk] += dum;      //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                int ild = ibas_ab_sym[h2][l][i];
                int kjd = ibas_ab_sym[h2][j][k];

                A_p[d2aboff[h2] + ild*gems_ab[h2]+kjd] -= dum;   // - D2ab(il,kj)

                //int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                //int ils = ibas_00_sym[h2][i][l];
                //int kjs = ibas_00_sym[h2][k][j];
                //A_p[d2soff[h2] + INDEX(ils,kjs)]             -= dum * 0.5;

                //if ( i != l && k != j ) {

                //    int sil = ( i < l ? 1 : -1 );
                //    int skj = ( k < j ? 1 : -1 );

                //    int ilt = ibas_aa_sym[h2][i][l];
                //    int kjt = ibas_aa_sym[h2][k][j];

                //    A_p[d2toff_m1[h2] + INDEX(ilt,kjt)] -= dum * sil * skj * 0.5;
                //}
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }
    // contraction condition:
    /*int na = nalpha_ - nrstc_ - nfrzc_;
    int nb = nbeta_ - nrstc_ - nfrzc_;
    int ms = (multiplicity_ - 1)/2;
    for (int h = 0; h < nirrep_; h++) {
        for (int klg = 0; klg < gems_ab[h]; klg++) {
    
            int k = bas_ab_sym[h][klg][0];
            int l = bas_ab_sym[h][klg][1];

            for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

                int i = bas_ab_sym[h][ijg][0];
                int j = bas_ab_sym[h][ijg][1];
                if ( i != j ) continue;


                A_p[g2toff[h] + ijg*gems_ab[h] + klg] += (na+nb)*u_p[offset + k*amo_+l];
                A_p[g2soff[h] + ijg*gems_ab[h] + klg] += -2.0 * ms * u_p[offset + k*amo_+l];

            }

        }
    }
    offset += amo_*amo_;*/
    // maximal spin constraint:
    /*for (int h = 0; h < nirrep_; h++) {
        for (int klg = 0; klg < gems_ab[h]; klg++) {

            int k = bas_ab_sym[h][klg][0];
            int l = bas_ab_sym[h][klg][1];

            double dum = u_p[offset + k * amo_ + l];

            for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

                int i = bas_ab_sym[h][ijg][0];
                int j = bas_ab_sym[h][ijg][1];
                if ( i != j ) continue;

                A_p[g2toff_p1[h] + ijg*gems_ab[h] + klg] += dum;
            }

        }
    }
    offset += amo_*amo_;
    for (int h = 0; h < nirrep_; h++) {
        for (int klg = 0; klg < gems_ab[h]; klg++) {

            int k = bas_ab_sym[h][klg][0];
            int l = bas_ab_sym[h][klg][1];

            double dum = u_p[offset + k * amo_ + l];

            for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

                int i = bas_ab_sym[h][ijg][0];
                int j = bas_ab_sym[h][ijg][1];
                if ( i != j ) continue;

                A_p[g2toff_p1[h] + klg*gems_ab[h] + ijg] += dum;
            }

        }
    }
    offset += amo_*amo_;*/
}

void v2RDMSolver::G2_constraints_guess(SharedVector u){

    double* u_p = u->pointer();

    // G2ab constraints:
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = 0.0;

                if ( j==l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum   +=  u_p[d1aoff[h3] + ii*amopi_[h3]+kk];      //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                int ild = ibas_ab_sym[h2][i][l];
                int kjd = ibas_ab_sym[h2][k][j];

                dum       -=  u_p[d2aboff[h2] + ild*gems_ab[h2]+kjd];   // - D2ab(il,kj)

                u_p[g2aboff[h] + ijg*gems_ab[h]+klg] = dum;
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }
    // G2ba constraints:
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = 0.0;

                if ( j==l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum       +=  u_p[d1boff[h3] + ii*amopi_[h3]+kk];      //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                int lid = ibas_ab_sym[h2][l][i];
                int jkd = ibas_ab_sym[h2][j][k];

                dum    -=  u_p[d2aboff[h2] + lid*gems_ab[h2]+jkd];       //   -D2ab(li,jk)

                u_p[g2baoff[h] + ijg*gems_ab[h]+klg] = dum;
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }
    // G2aaaa / G2aabb / G2bbaa / G2bbbb
    for (int h = 0; h < nirrep_; h++) {
        // G2aaaa
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                double dum = 0.0;

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum       +=  u_p[d1aoff[h3] + ii*amopi_[h3]+kk];    //   D1(i,k) djl
                }

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];

                    dum       -=  u_p[d2aaoff[h2] + ild*gems_aa[h2]+kjd] * sil * skj; // -D2aa(il,kj)

                }

                u_p[g2aaoff[h] + ijg*2*gems_ab[h]+klg] = dum;
            }
        }
        // G2bbbb
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                double dum = 0.0;

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum       +=  u_p[d1boff[h3] + ii*amopi_[h3]+kk];    //   D1(i,k) djl
                }

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];

                    dum       -=  u_p[d2bboff[h2] + ild*gems_aa[h2]+kjd] * sil * skj; // -D2bb(il,kj)

                }

                u_p[g2aaoff[h] + (gems_ab[h] + ijg)*2*gems_ab[h] + (gems_ab[h] + klg)] = dum;

            }
        }
        // G2aabb
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                double dum = 0.0;

                int ild = ibas_ab_sym[h2][i][l];
                int jkd = ibas_ab_sym[h2][j][k];

                dum       +=  u_p[d2aboff[h2] + ild*gems_ab[h2]+jkd]; // D2ab(il,jk)

                u_p[g2aaoff[h] + (ijg)*2*gems_ab[h] + (gems_ab[h] + klg)] = dum;

            }
        }
        // G2bbaa
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                double dum = 0.0;

                int lid = ibas_ab_sym[h2][l][i];
                int kjd = ibas_ab_sym[h2][k][j];

                dum       +=  u_p[d2aboff[h2] + lid*gems_ab[h2]+kjd]; // D2ab(li,kj)

                u_p[g2aaoff[h] + (gems_ab[h] + ijg)*2*gems_ab[h] + (klg)] = dum;
            }
        }
        offset += 2*gems_ab[h]*2*gems_ab[h];
    }
}

// G2 portion of A.x (with symmetry)
void v2RDMSolver::G2_constraints_Au(SharedVector A,SharedVector u){
    double* A_p = A->pointer();
    double* u_p = u->pointer();

    // G2ab constraints:
// heyheyhey
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];


                double dum = -u_p[g2aboff[h] + ijg*gems_ab[h]+klg];    // - G2ab(ij,kl)

                if ( j==l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum   +=  u_p[d1aoff[h3] + ii*amopi_[h3]+kk];      //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                int ild = ibas_ab_sym[h2][i][l];
                int kjd = ibas_ab_sym[h2][k][j];

                dum       -=  u_p[d2aboff[h2] + ild*gems_ab[h2]+kjd];   // - D2ab(il,kj)

                A_p[offset + ijg*gems_ab[h]+klg] = dum;
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }
    // G2ba constraints:
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = -u_p[g2baoff[h] + ijg*gems_ab[h]+klg];        // - G2ba(ij,kl)

                if ( j==l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum       +=  u_p[d1boff[h3] + ii*amopi_[h3]+kk];      //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                int lid = ibas_ab_sym[h2][l][i];
                int jkd = ibas_ab_sym[h2][j][k];

                dum    -=  u_p[d2aboff[h2] + lid*gems_ab[h2]+jkd];       //   -D2ab(li,jk)

                A_p[offset + ijg*gems_ab[h]+klg] = dum;
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }
    // G2aaaa / G2aabb / G2bbaa / G2bbbb
    for (int h = 0; h < nirrep_; h++) {
        // G2aaaa
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                double dum = -u_p[g2aaoff[h] + ijg*2*gems_ab[h]+klg];       // - G2aaaa(ij,kl)
                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum       +=  u_p[d1aoff[h3] + ii*amopi_[h3]+kk];    //   D1(i,k) djl
                }

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];

                    dum       -=  u_p[d2aaoff[h2] + ild*gems_aa[h2]+kjd] * sil * skj; // -D2aa(il,kj)

                }

                A_p[offset + ijg*2*gems_ab[h]+klg] = dum;
            }
        }
        // G2bbbb
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                double dum = -u_p[g2aaoff[h] + (gems_ab[h] + ijg)*2*gems_ab[h] + (gems_ab[h] + klg)]; // - G2bbbb(ij,kl)
                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum       +=  u_p[d1boff[h3] + ii*amopi_[h3]+kk];    //   D1(i,k) djl
                }

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];

                    dum       -=  u_p[d2bboff[h2] + ild*gems_aa[h2]+kjd] * sil * skj; // -D2bb(il,kj)

                }

                A_p[offset + (gems_ab[h] + ijg)*2*gems_ab[h] + (gems_ab[h] + klg)] = dum;
            }
        }
        // G2aabb
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                double dum = -u_p[g2aaoff[h] + (ijg)*2*gems_ab[h] + (gems_ab[h] + klg)];       // - G2aabb(ij,kl)

                int ild = ibas_ab_sym[h2][i][l];
                int jkd = ibas_ab_sym[h2][j][k];

                dum       +=  u_p[d2aboff[h2] + ild*gems_ab[h2]+jkd]; // D2ab(il,jk)

                A_p[offset + (ijg)*2*gems_ab[h] + (gems_ab[h] + klg)] = dum;
            }
        }
        // G2bbaa
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                double dum = -u_p[g2aaoff[h] + (gems_ab[h] + ijg)*2*gems_ab[h] + (klg)];       // - G2bbaa(ij,kl)

                int lid = ibas_ab_sym[h2][l][i];
                int kjd = ibas_ab_sym[h2][k][j];

                dum       +=  u_p[d2aboff[h2] + lid*gems_ab[h2]+kjd]; // D2ab(li,kj)

                A_p[offset + (gems_ab[h] + ijg)*2*gems_ab[h] + (klg)] = dum;
            }
        }
        offset += 2*gems_ab[h]*2*gems_ab[h];
    }

}

// G2 portion of A^T.y (with symmetry)
void v2RDMSolver::G2_constraints_ATu(SharedVector A,SharedVector u){

    double* A_p = A->pointer();
    double* u_p = u->pointer();

    // G2ab constraints:
// heyheyhey
    for (int h = 0; h < nirrep_; h++) {
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = u_p[offset + ijg*gems_ab[h]+klg];

                A_p[g2aboff[h] + ijg*gems_ab[h]+klg] -= dum;    // - G2ab(ij,kl)

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    A_p[d1aoff[h3] + ii*amopi_[h3]+kk] += dum;      //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                int ild = ibas_ab_sym[h2][i][l];
                int kjd = ibas_ab_sym[h2][k][j];

                A_p[d2aboff[h2] + ild*gems_ab[h2]+kjd] -= dum;   // - D2ab(il,kj)
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }
    // G2ba constraints:
    for (int h = 0; h < nirrep_; h++) {
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = u_p[offset + ijg*gems_ab[h]+klg];

                A_p[g2baoff[h] + ijg*gems_ab[h]+klg]  -= dum;

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    A_p[d1boff[h3] + ii*amopi_[h3]+kk]      += dum;
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                int lid = ibas_ab_sym[h2][l][i];
                int jkd = ibas_ab_sym[h2][j][k];

                A_p[d2aboff[h2] + lid*gems_ab[h2]+jkd] -= dum;
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }
    // G2aaaa / G2aabb / G2bbaa / G2bbbb constraints:
    for (int h = 0; h < nirrep_; h++) {
        // G2aaaa
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = u_p[offset + ijg*2*gems_ab[h]+klg];

                A_p[g2aaoff[h] + ijg*2*gems_ab[h]+klg] -= dum;

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    A_p[d1aoff[h3] + ii*amopi_[h3]+kk]         += dum;
                }

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];

                    A_p[d2aaoff[h2] + ild*gems_aa[h2]+kjd] -= dum * sil * skj;
                }
           }
        }
        // G2bbbb
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = u_p[offset + (gems_ab[h] + ijg)*2*gems_ab[h]+(gems_ab[h] + klg)];

                A_p[g2aaoff[h] + (gems_ab[h] + ijg)*2*gems_ab[h]+(gems_ab[h] + klg)] -= dum;

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    A_p[d1boff[h3] + ii*amopi_[h3]+kk]         += dum;
                }

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];

                    A_p[d2bboff[h2] + ild*gems_aa[h2]+kjd] -= dum * sil * skj;
                }
            }
        }
        // G2aabb
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = u_p[offset + ijg*2*gems_ab[h]+(klg + gems_ab[h])];

                A_p[g2aaoff[h] + ijg*2*gems_ab[h]+(klg + gems_ab[h])] -= dum;

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                int ild = ibas_ab_sym[h2][i][l];
                int jkd = ibas_ab_sym[h2][j][k];

                A_p[d2aboff[h2] + ild*gems_ab[h2]+jkd] += dum;
            }
        }
        // G2bbaa
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = u_p[offset + (ijg + gems_ab[h])*2*gems_ab[h]+klg];

                A_p[g2aaoff[h] + (ijg + gems_ab[h])*2*gems_ab[h]+klg] -= dum;

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                int lid = ibas_ab_sym[h2][l][i];
                int kjd = ibas_ab_sym[h2][k][j];

                A_p[d2aboff[h2] + lid*gems_ab[h2]+kjd] += dum;
            }
        }

        offset += 2*gems_ab[h]*2*gems_ab[h];
    }

}

} // end namespaces
