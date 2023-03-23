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

void v2RDMSolver::G2_constraints_guess(double* u){

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
                    dum   +=  u[d1aoff[h3] + ii*amopi_[h3]+kk];      //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                int ild = ibas_ab_sym[h2][i][l];
                int kjd = ibas_ab_sym[h2][k][j];

                dum       -=  u[d2aboff[h2] + ild*gems_ab[h2]+kjd];   // - D2ab(il,kj)

                u[g2aboff[h] + ijg*gems_ab[h]+klg] = dum;
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
                    dum       +=  u[d1boff[h3] + ii*amopi_[h3]+kk];      //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                int lid = ibas_ab_sym[h2][l][i];
                int jkd = ibas_ab_sym[h2][j][k];

                dum    -=  u[d2aboff[h2] + lid*gems_ab[h2]+jkd];       //   -D2ab(li,jk)

                u[g2baoff[h] + ijg*gems_ab[h]+klg] = dum;
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
                    dum       +=  u[d1aoff[h3] + ii*amopi_[h3]+kk];    //   D1(i,k) djl
                }

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];

                    dum       -=  u[d2aaoff[h2] + ild*gems_aa[h2]+kjd] * sil * skj; // -D2aa(il,kj)

                }

                u[g2aaoff[h] + ijg*2*gems_ab[h]+klg] = dum;
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
                    dum       +=  u[d1boff[h3] + ii*amopi_[h3]+kk];    //   D1(i,k) djl
                }

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];

                    dum       -=  u[d2bboff[h2] + ild*gems_aa[h2]+kjd] * sil * skj; // -D2bb(il,kj)

                }

                u[g2aaoff[h] + (gems_ab[h] + ijg)*2*gems_ab[h] + (gems_ab[h] + klg)] = dum;

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

                dum       +=  u[d2aboff[h2] + ild*gems_ab[h2]+jkd]; // D2ab(il,jk)

                u[g2aaoff[h] + (ijg)*2*gems_ab[h] + (gems_ab[h] + klg)] = dum;

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

                dum       +=  u[d2aboff[h2] + lid*gems_ab[h2]+kjd]; // D2ab(li,kj)

                u[g2aaoff[h] + (gems_ab[h] + ijg)*2*gems_ab[h] + (klg)] = dum;
            }
        }
        offset += 2*gems_ab[h]*2*gems_ab[h];
    }
}

// G2 portion of A.x (with symmetry)
void v2RDMSolver::G2_constraints_Au(double* A,double* u){

    // G2ab constraints:
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];


                double dum = -u[g2aboff[h] + ijg*gems_ab[h]+klg];    // - G2ab(ij,kl)

                if ( j==l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum   +=  u[d1aoff[h3] + ii*amopi_[h3]+kk];      //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                int ild = ibas_ab_sym[h2][i][l];
                int kjd = ibas_ab_sym[h2][k][j];

                dum       -=  u[d2aboff[h2] + ild*gems_ab[h2]+kjd];   // - D2ab(il,kj)

                A[offset + ijg*gems_ab[h]+klg] = dum;
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

                double dum = -u[g2baoff[h] + ijg*gems_ab[h]+klg];        // - G2ba(ij,kl)

                if ( j==l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum       +=  u[d1boff[h3] + ii*amopi_[h3]+kk];      //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                int lid = ibas_ab_sym[h2][l][i];
                int jkd = ibas_ab_sym[h2][j][k];

                dum    -=  u[d2aboff[h2] + lid*gems_ab[h2]+jkd];       //   -D2ab(li,jk)

                A[offset + ijg*gems_ab[h]+klg] = dum;
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

                double dum = -u[g2aaoff[h] + ijg*2*gems_ab[h]+klg];       // - G2aaaa(ij,kl)
                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum       +=  u[d1aoff[h3] + ii*amopi_[h3]+kk];    //   D1(i,k) djl
                }

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];

                    dum       -=  u[d2aaoff[h2] + ild*gems_aa[h2]+kjd] * sil * skj; // -D2aa(il,kj)

                }

                A[offset + ijg*2*gems_ab[h]+klg] = dum;
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

                double dum = -u[g2aaoff[h] + (gems_ab[h] + ijg)*2*gems_ab[h] + (gems_ab[h] + klg)]; // - G2bbbb(ij,kl)
                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    dum       +=  u[d1boff[h3] + ii*amopi_[h3]+kk];    //   D1(i,k) djl
                }

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];

                    dum       -=  u[d2bboff[h2] + ild*gems_aa[h2]+kjd] * sil * skj; // -D2bb(il,kj)

                }

                A[offset + (gems_ab[h] + ijg)*2*gems_ab[h] + (gems_ab[h] + klg)] = dum;
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

                double dum = -u[g2aaoff[h] + (ijg)*2*gems_ab[h] + (gems_ab[h] + klg)];       // - G2aabb(ij,kl)

                int ild = ibas_ab_sym[h2][i][l];
                int jkd = ibas_ab_sym[h2][j][k];

                dum       +=  u[d2aboff[h2] + ild*gems_ab[h2]+jkd]; // D2ab(il,jk)

                A[offset + (ijg)*2*gems_ab[h] + (gems_ab[h] + klg)] = dum;
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

                double dum = -u[g2aaoff[h] + (gems_ab[h] + ijg)*2*gems_ab[h] + (klg)];       // - G2bbaa(ij,kl)

                int lid = ibas_ab_sym[h2][l][i];
                int kjd = ibas_ab_sym[h2][k][j];

                dum       +=  u[d2aboff[h2] + lid*gems_ab[h2]+kjd]; // D2ab(li,kj)

                A[offset + (gems_ab[h] + ijg)*2*gems_ab[h] + (klg)] = dum;
            }
        }
        offset += 2*gems_ab[h]*2*gems_ab[h];
    }

}

// G2 portion of A^T.y (with symmetry)
void v2RDMSolver::G2_constraints_ATu(double* A,double* u){

    // G2ab constraints:
    for (int h = 0; h < nirrep_; h++) {
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = u[offset + ijg*gems_ab[h]+klg];

                A[g2aboff[h] + ijg*gems_ab[h]+klg] -= dum;    // - G2ab(ij,kl)

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    A[d1aoff[h3] + ii*amopi_[h3]+kk] += dum;      //   D1(i,k) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                int ild = ibas_ab_sym[h2][i][l];
                int kjd = ibas_ab_sym[h2][k][j];

                A[d2aboff[h2] + ild*gems_ab[h2]+kjd] -= dum;   // - D2ab(il,kj)
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

                double dum = u[offset + ijg*gems_ab[h]+klg];

                A[g2baoff[h] + ijg*gems_ab[h]+klg]  -= dum;

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    A[d1boff[h3] + ii*amopi_[h3]+kk]      += dum;
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                int lid = ibas_ab_sym[h2][l][i];
                int jkd = ibas_ab_sym[h2][j][k];

                A[d2aboff[h2] + lid*gems_ab[h2]+jkd] -= dum;
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

                double dum = u[offset + ijg*2*gems_ab[h]+klg];

                A[g2aaoff[h] + ijg*2*gems_ab[h]+klg] -= dum;

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    A[d1aoff[h3] + ii*amopi_[h3]+kk]         += dum;
                }

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];

                    A[d2aaoff[h2] + ild*gems_aa[h2]+kjd] -= dum * sil * skj;
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

                double dum = u[offset + (gems_ab[h] + ijg)*2*gems_ab[h]+(gems_ab[h] + klg)];

                A[g2aaoff[h] + (gems_ab[h] + ijg)*2*gems_ab[h]+(gems_ab[h] + klg)] -= dum;

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    A[d1boff[h3] + ii*amopi_[h3]+kk]         += dum;
                }

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];

                    A[d2bboff[h2] + ild*gems_aa[h2]+kjd] -= dum * sil * skj;
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

                double dum = u[offset + ijg*2*gems_ab[h]+(klg + gems_ab[h])];

                A[g2aaoff[h] + ijg*2*gems_ab[h]+(klg + gems_ab[h])] -= dum;

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                int ild = ibas_ab_sym[h2][i][l];
                int jkd = ibas_ab_sym[h2][j][k];

                A[d2aboff[h2] + ild*gems_ab[h2]+jkd] += dum;
            }
        }
        // G2bbaa
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = u[offset + (ijg + gems_ab[h])*2*gems_ab[h]+klg];

                A[g2aaoff[h] + (ijg + gems_ab[h])*2*gems_ab[h]+klg] -= dum;

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                int lid = ibas_ab_sym[h2][l][i];
                int kjd = ibas_ab_sym[h2][k][j];

                A[d2aboff[h2] + lid*gems_ab[h2]+kjd] += dum;
            }
        }

        offset += 2*gems_ab[h]*2*gems_ab[h];
    }

}

// G2 portion of A.u for the sos problem, which looks an awful lot like AT.u for the primal problem 
void v2RDMSolver::G2_constraints_Au_sos(double* A,double* u){

    // diagonal contribution to H (E)

    double na = nalpha_ - nrstc_ - nfrzc_;
    double nb = nbeta_ - nrstc_ - nfrzc_;
    double n = na + nb;

    // ab
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            A[d2aboff[h] + ij * gems_ab[h] + ij] -= u[sos_e2ab_off];
        }
    }
    // aa, bb
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            A[d2aaoff[h] + ij * gems_aa[h] + ij] -= u[sos_e2aa_off];
            A[d2bboff[h] + ij * gems_aa[h] + ij] -= u[sos_e2bb_off];
        }
    }

    for (int h = 0; h < nirrep_; h++) {
        // D2ab contribution to H:
        C_DAXPY(gems_ab[h] * gems_ab[h], 1.0, u + d2aboff[h], 1, A + d2aboff[h], 1);

        // D2aa contribution to H:
        C_DAXPY(gems_aa[h] * gems_aa[h], 1.0, u + d2aaoff[h], 1, A + d2aaoff[h], 1);

        // D2bb contribution to H:
        C_DAXPY(gems_aa[h] * gems_aa[h], 1.0, u + d2bboff[h], 1, A + d2bboff[h], 1);
    }

    // G2ab contributions to H:
    for (int h = 0; h < nirrep_; h++) {
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = u[g2aboff[h] + ijg*gems_ab[h]+klg];        // G2ba(ij,kl)

                if ( j==l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    A[d1aoff[h3] + ii*amopi_[h3]+kk] += dum;   //   H(i,k) += G2ab(ij,kl) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                int ild = ibas_ab_sym[h2][i][l];
                int kjd = ibas_ab_sym[h2][k][j];

                A[d2aboff[h2] + ild*gems_ab[h2]+kjd] -= dum;   // D2ab(il,kj) -= G2ab(ij,kl)
            }
        }
    }

    // G2ba contributions to H:
    for (int h = 0; h < nirrep_; h++) {
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                double dum = u[g2baoff[h] + ijg*gems_ab[h]+klg];        // G2ba(ij,kl)

                if ( j==l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    A[d1boff[h3] + ii*amopi_[h3]+kk] += dum;      //   H1(i,k) += G2ba(ij,kl) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                int lid = ibas_ab_sym[h2][l][i];
                int jkd = ibas_ab_sym[h2][j][k];

                A[d2aboff[h2] + lid*gems_ab[h2]+jkd] -= dum;       //   H2(li,jk) -= G2ba(ij,kl)
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }

    // G2aaaa / G2aabb / G2bbaa / G2bbbb contributions to H
    for (int h = 0; h < nirrep_; h++) {
        // G2aaaa
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                double dum = u[g2aaoff[h] + ijg*2*gems_ab[h]+klg];       // G2aaaa(ij,kl)

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    A[d1aoff[h3] + ii*amopi_[h3]+kk] += dum;    //   H1(i,k) += G2aaaa(ij,kl) djl
                }

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];

                    A[d2aaoff[h2] + ild*gems_aa[h2]+kjd] -= dum * sil * skj; // H2aa(il,kj) -= G2aaaa(ij,kl)

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

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                double dum = u[g2aaoff[h] + (gems_ab[h] + ijg)*2*gems_ab[h] + (gems_ab[h] + klg)]; // G2bbbb(ij,kl)

                if ( j == l ) {
                    int h3 = symmetry[i];
                    int ii = i - pitzer_offset[h3];
                    int kk = k - pitzer_offset[h3];
                    A[d1boff[h3] + ii*amopi_[h3]+kk] += dum;    //   H1(i,k) += G2bbbb(ij,kl) djl
                }

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];

                    A[d2bboff[h2] + ild*gems_aa[h2]+kjd] -= dum * sil * skj; // H2bb(il,kj) -= G2bbbb(ij,kl)

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

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                double dum = u[g2aaoff[h] + (ijg)*2*gems_ab[h] + (gems_ab[h] + klg)];       // G2aabb(ij,kl)

                int ild = ibas_ab_sym[h2][i][l];
                int jkd = ibas_ab_sym[h2][j][k];

                A[d2aboff[h2] + ild*gems_ab[h2]+jkd] += dum; // H2ab(il,jk) += G2aabb(ij,kl)
            }
        }

        // G2bbaa
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                double dum = u[g2aaoff[h] + (gems_ab[h] + ijg)*2*gems_ab[h] + (klg)]; // G2bbaa(ij,kl)

                int lid = ibas_ab_sym[h2][l][i];
                int kjd = ibas_ab_sym[h2][k][j];

                A[d2aboff[h2] + lid*gems_ab[h2]+kjd] += dum; // H2ab(li,kj) += G2bbaa(ij,kl)
            }
        }
    }

    // adjust A such that d1 elements are mapped onto d2
    for (int h = 0; h < nirrep_; h++) {
        for (long int ij = 0; ij < gems_ab[h]; ij++) {
            long int i = bas_ab_sym[h][ij][0];
            long int j = bas_ab_sym[h][ij][1];

            long int ii = full_basis[i];
            long int jj = full_basis[j];

            int hi = symmetry[i];
            int hj = symmetry[j];

            for (long int kl = 0; kl < gems_ab[h]; kl++) {
                long int k = bas_ab_sym[h][kl][0];
                long int l = bas_ab_sym[h][kl][1];

                long int kk = full_basis[k];
                long int ll = full_basis[l];

                int hik = SymmetryPair(symmetry[i],symmetry[k]);

                int hk = symmetry[k];
                int hl = symmetry[l];

                int iii = i - pitzer_offset[hi];
                int kkk = k - pitzer_offset[hk];
                int jjj = j - pitzer_offset[hj];
                int lll = l - pitzer_offset[hl];


                if ( j == l ) {
                    A[d2aboff[h] + ij*gems_ab[h]+kl] += 1.0 / (na + nb - 1) * A[d1aoff[hi] + iii * amopi_[hi] + kkk];
                }
                if ( i == k ) {
                    A[d2aboff[h] + ij*gems_ab[h]+kl] += 1.0 / (na + nb - 1) * A[d1boff[hj] + jjj * amopi_[hj] + lll];
                }

            }
        }
    }

    for (int h = 0; h < nirrep_; h++) {
        for (long int ij = 0; ij < gems_aa[h]; ij++) {
            long int i = bas_aa_sym[h][ij][0];
            long int j = bas_aa_sym[h][ij][1];

            long int ii = full_basis[i];
            long int jj = full_basis[j];

            int hi = symmetry[i];
            int hj = symmetry[j];

            for (long int kl = 0; kl < gems_aa[h]; kl++) {
                long int k = bas_aa_sym[h][kl][0];
                long int l = bas_aa_sym[h][kl][1];

                long int kk = full_basis[k];
                long int ll = full_basis[l];

                int hk = symmetry[k];
                int hl = symmetry[l];

                int hik = SymmetryPair(symmetry[i],symmetry[k]);

                long int ijb = ibas_ab_sym[h][i][j];
                long int jib = ibas_ab_sym[h][j][i];

                long int klb = ibas_ab_sym[h][k][l];
                long int lkb = ibas_ab_sym[h][l][k];

                int iii = i - pitzer_offset[hi];
                int kkk = k - pitzer_offset[hk];
                int jjj = j - pitzer_offset[hj];
                int lll = l - pitzer_offset[hl];

                if ( j == l ) {
                    A[d2aaoff[h] + ij*gems_aa[h]+kl] += 1.0 / (na + nb - 1) * A[d1aoff[hi] + iii * amopi_[hi] + kkk];
                    A[d2bboff[h] + ij*gems_aa[h]+kl] += 1.0 / (na + nb - 1) * A[d1boff[hi] + iii * amopi_[hi] + kkk];
                }
                if ( i == k ) {
                    A[d2aaoff[h] + ij*gems_aa[h]+kl] += 1.0 / (na + nb - 1) * A[d1aoff[hj] + jjj * amopi_[hj] + lll];
                    A[d2bboff[h] + ij*gems_aa[h]+kl] += 1.0 / (na + nb - 1) * A[d1boff[hj] + jjj * amopi_[hj] + lll];
                }
                if ( j == k ) {
                    A[d2aaoff[h] + ij*gems_aa[h]+kl] -= 1.0 / (na + nb - 1) * A[d1aoff[hi] + iii * amopi_[hi] + lll];
                    A[d2bboff[h] + ij*gems_aa[h]+kl] -= 1.0 / (na + nb - 1) * A[d1boff[hi] + iii * amopi_[hi] + lll];
                }
                if ( i == l ) {
                    A[d2aaoff[h] + ij*gems_aa[h]+kl] -= 1.0 / (na + nb - 1) * A[d1aoff[hj] + jjj * amopi_[hj] + kkk];
                    A[d2bboff[h] + ij*gems_aa[h]+kl] -= 1.0 / (na + nb - 1) * A[d1boff[hj] + jjj * amopi_[hj] + kkk];
                }

            }
        }
    }

    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < amopi_[h] * amopi_[h]; ij++) {
            A[d1aoff[h] + ij] = 0.0;
            A[d1boff[h] + ij] = 0.0;
        }
    }

}

// G2 portion of AT.u for the sos problem, which looks an awful lot like A.u for the primal problem
void v2RDMSolver::G2_constraints_ATu_sos(double* A,double* u){

    // diagonal contribution to H (E)

    double na = nalpha_ - nrstc_ - nfrzc_;
    double nb = nbeta_ - nrstc_ - nfrzc_;
    double n = na + nb;

    // ab
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            A[sos_e2ab_off] -= u[d2aboff[h] + ij * gems_ab[h] + ij];
        }
    }

    // aa, bb
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            A[sos_e2aa_off] -= u[d2aaoff[h] + ij * gems_aa[h] + ij];
            A[sos_e2bb_off] -= u[d2bboff[h] + ij * gems_aa[h] + ij];
        }
    }

    for (int h = 0; h < nirrep_; h++) {
        // D2ab contribution to H:
        C_DAXPY(gems_ab[h] * gems_ab[h], 1.0, u + d2aboff[h], 1, A + d2aboff[h], 1);

        // D2aa contribution to H:
        C_DAXPY(gems_aa[h] * gems_aa[h], 1.0, u + d2aaoff[h], 1, A + d2aaoff[h], 1);

        // D2bb contribution to H:
        C_DAXPY(gems_aa[h] * gems_aa[h], 1.0, u + d2bboff[h], 1, A + d2bboff[h], 1);

    }

    // adjust u such that d2 elements are mapped onto d1

    int poff = 0;

    // contraction: D2aa + D2ab -> D1 a
    poff = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++){
            for (int j = 0; j < amopi_[h]; j++){
                double sum = 0.0; 
                int ii  = i + poff;
                int jj  = j + poff;
                for(int k = 0; k < amo_; k++){
                    int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                    int ik = ibas_ab_sym[h2][ii][k];
                    int jk = ibas_ab_sym[h2][jj][k];
                    sum += u[d2aboff[h2] + ik*gems_ab[h2]+jk];
                }
                for(int k = 0; k < amo_; k++){
                    if( ii==k || jj==k ) continue;
                    int h2   = SymmetryPair(symmetry[ii],symmetry[k]);
                    int ik  = ibas_aa_sym[h2][ii][k];
                    int jk  = ibas_aa_sym[h2][jj][k];
                    int sik = ( ii < k ) ? 1 : -1;
                    int sjk = ( jj < k ) ? 1 : -1;
                    sum += sik*sjk*u[d2aaoff[h2] + ik*gems_aa[h2]+jk];
                }
                u[d1aoff[h] + i*amopi_[h]+j] = sum / (n - 1.0);
            }
        }
        offset += amopi_[h]*amopi_[h];
        poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
    }

    // contraction: D2bb + D2ab -> D1 b
    poff = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++){
            for (int j = 0; j < amopi_[h]; j++){
                double sum = 0.0;
                int ii  = i + poff;
                int jj  = j + poff;
                for(int k = 0; k < amo_; k++){
                    int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                    int ik = ibas_ab_sym[h2][k][ii];
                    int jk = ibas_ab_sym[h2][k][jj];
                    sum += u[d2aboff[h2] + ik*gems_ab[h2]+jk];
                }
                for(int k = 0; k < amo_; k++){
                    if( ii==k || jj==k ) continue;
                    int h2   = SymmetryPair(symmetry[ii],symmetry[k]);
                    int ik  = ibas_aa_sym[h2][ii][k];
                    int jk  = ibas_aa_sym[h2][jj][k];
                    int sik = ( ii < k ) ? 1 : -1;
                    int sjk = ( jj < k ) ? 1 : -1;
                    sum += sik*sjk*u[d2bboff[h2] + ik*gems_aa[h2]+jk];
                }
                u[d1boff[h] + i*amopi_[h]+j] = sum / (n - 1.0);
            }
        }
        offset += amopi_[h]*amopi_[h];
        poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
    }

    // G2ab contributions to H:
    for (int h = 0; h < nirrep_; h++) {
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
                    dum += u[d1aoff[h3] + ii*amopi_[h3]+kk];   //   H(i,k) += G2ab(ij,kl) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                int ild = ibas_ab_sym[h2][i][l];
                int kjd = ibas_ab_sym[h2][k][j];

                dum -= u[d2aboff[h2] + ild*gems_ab[h2]+kjd];   // D2ab(il,kj) -= G2ab(ij,kl)

                A[g2aboff[h] + ijg*gems_ab[h]+klg] += dum;        // G2ba(ij,kl)
            }
        }
    }

    // G2ba contributions to H:
    for (int h = 0; h < nirrep_; h++) {
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
                    dum += u[d1boff[h3] + ii*amopi_[h3]+kk];      //   H1(i,k) += G2ba(ij,kl) djl
                }

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                int lid = ibas_ab_sym[h2][l][i];
                int jkd = ibas_ab_sym[h2][j][k];

                dum -= u[d2aboff[h2] + lid*gems_ab[h2]+jkd];       //   H2(li,jk) -= G2ba(ij,kl)

                A[g2baoff[h] + ijg*gems_ab[h]+klg] += dum;         // G2ba(ij,kl)
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }

    // G2aaaa / G2aabb / G2bbaa / G2bbbb contributions to H
    for (int h = 0; h < nirrep_; h++) {

        // G2aaaa
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
                    dum += u[d1aoff[h3] + ii*amopi_[h3]+kk];    //   H1(i,k) += G2aaaa(ij,kl) djl
                }

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];

                    dum -= u[d2aaoff[h2] + ild*gems_aa[h2]+kjd] * sil * skj; // H2aa(il,kj) -= G2aaaa(ij,kl)

                }

                A[g2aaoff[h] + ijg*2*gems_ab[h]+klg] += dum;       // G2aaaa(ij,kl)
            }
        }

        // G2bbbb
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
                    dum += u[d1boff[h3] + ii*amopi_[h3]+kk];    //   H1(i,k) += G2bbbb(ij,kl) djl
                }

                if ( i != l && k != j ) {

                    int sil = ( i < l ? 1 : -1 );
                    int skj = ( k < j ? 1 : -1 );

                    int ild = ibas_aa_sym[h2][i][l];
                    int kjd = ibas_aa_sym[h2][k][j];

                    dum -= u[d2bboff[h2] + ild*gems_aa[h2]+kjd] * sil * skj; // H2bb(il,kj) -= G2bbbb(ij,kl)

                }

                A[g2aaoff[h] + (gems_ab[h] + ijg)*2*gems_ab[h] + (gems_ab[h] + klg)] += dum; // G2bbbb(ij,kl)
            }
        }

        // G2aabb
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                int ild = ibas_ab_sym[h2][i][l];
                int jkd = ibas_ab_sym[h2][j][k];

                A[g2aaoff[h] + (ijg)*2*gems_ab[h] + (gems_ab[h] + klg)] += u[d2aboff[h2] + ild*gems_ab[h2]+jkd]; // H2ab(il,jk) += G2aabb(ij,kl)
            }
        }

        // G2bbaa
        for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

            int i = bas_ab_sym[h][ijg][0];
            int j = bas_ab_sym[h][ijg][1];

            for (int klg = 0; klg < gems_ab[h]; klg++) {

                int k = bas_ab_sym[h][klg][0];
                int l = bas_ab_sym[h][klg][1];

                int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                int lid = ibas_ab_sym[h2][l][i];
                int kjd = ibas_ab_sym[h2][k][j];

                A[g2aaoff[h] + (gems_ab[h] + ijg)*2*gems_ab[h] + (klg)] += u[d2aboff[h2] + lid*gems_ab[h2]+kjd]; // H2ab(li,kj) += G2bbaa(ij,kl)
            }
        }
    }

    // re-zero d1 part of u
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < amopi_[h] * amopi_[h]; ij++) {
            u[d1aoff[h] + ij] = 0.0;
            u[d1boff[h] + ij] = 0.0;
        }
    }

}

} // end namespaces
