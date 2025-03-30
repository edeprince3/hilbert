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

#include "v2rdm_doci_solver.h"

#include <misc/omp.h>

using namespace psi;

namespace hilbert{

// T2 portion of A.x (with symmetry)
void v2RDM_DOCISolver::T2_constraints_Au(double* A,double* u){

    // Rubio-Garcia et. al JCTC (2018): Eq. 33
    //
    //
    //  D(ac)            | -dac Pi(ab)                            | D(ab)
    //  -dac Pi(ba)      | Pi(ac) + dac ( -D(ba) - D(ab) + p(b) ) | Pi(ab)
    //  D(bc)            | Pi(bc)                                 | p(b)
    //
    // 

    // 
    //  D(jk)            | -djk Pi(ik)                            | D(ik)
    //  -djk Pi(ki)      | Pi(jk) + djk ( -D(ik) - D(ki) + p(i) ) | Pi(ik)
    //  D(ji)            | Pi(ji)                                 | p(i)
    // 
    int aux = 0;
    for (int i = 0; i < amo_; i++) {
        // first block: D(jk)
        for (int j = 0; j < amo_; j++) {
            if ( i == j ) continue;
            int myj = j;
            if ( j > i ) myj--;
            for (int k = 0; k < amo_; k++) {
                if ( i == k ) continue;
                int myk = k;
                if ( k > i ) myk--;
                int id = myk * ( 2 * amo_ - 1 ) + myj;

                double dum = -u[t2s1off_ + id + aux];
                if ( j == k ) {
                    dum += u[d1off_ + j];
                }else {
                    dum += 0.5 * u[d2s2off_ + j * amo_ + k];
                    dum += 0.5 * u[d2s2off_ + k * amo_ + j];
                }

                A[offset + id] = dum;
            }
        }
        // second block: -djk Pi(ik) 
        for (int j = 0; j < amo_; j++) {
            if ( i == j ) continue;
            int myj = j;
            if ( j > i ) myj--;
            for (int k = 0; k < amo_; k++) {
                if ( i == k ) continue;
                int myk = k;
                if ( k > i ) myk--;
                int id = myk * ( 2 * amo_ - 1 ) + myj + amo_ - 1;

                double dum = -u[t2s1off_ + id + aux];
                if ( j == k ) {
                    dum -= 0.5 * u[d2s0off_ + i * amo_ + k];
                    dum -= 0.5 * u[d2s0off_ + k * amo_ + i];
                }

                A[offset + id] = dum;
            }
        }
        // third block: D(ik)
        for (int k = 0; k < amo_; k++) {
            if ( i == k ) continue;
            int myk = k;
            if ( k > i ) myk--;
            int id = myk * ( 2 * amo_ - 1 ) + 2 * amo_ - 2;

            double dum = -u[t2s1off_ + id + aux];
            dum += 0.5 * u[d2s2off_ + i * amo_ + k];
            dum += 0.5 * u[d2s2off_ + k * amo_ + i];

            A[offset + id] = dum;
        }
        // fourth block: -djk Pi(ki) 
        for (int j = 0; j < amo_; j++) {
            if ( i == j ) continue;
            int myj = j;
            if ( j > i ) myj--;
            for (int k = 0; k < amo_; k++) {
                if ( i == k ) continue;
                int myk = k;
                if ( k > i ) myk--;
                int id = (myk + amo_ - 1) * ( 2 * amo_ - 1 ) + myj;

                double dum = -u[t2s1off_ + id + aux];
                if ( j == k ) {
                    dum -= 0.5 * u[d2s0off_ + i * amo_ + k];
                    dum -= 0.5 * u[d2s0off_ + k * amo_ + i];
                }

                A[offset + id] = dum;
            }
        }
        // fifth block: Pi(jk) + djk ( -D(ik) - D(ki) + p(i) )
        for (int j = 0; j < amo_; j++) {
            if ( i == j ) continue;
            int myj = j;
            if ( j > i ) myj--;
            for (int k = 0; k < amo_; k++) {
                if ( i == k ) continue;
                int myk = k;
                if ( k > i ) myk--;
                int id = (myk + amo_ - 1) * ( 2 * amo_ - 1 ) + myj + amo_ - 1;

                double dum = -u[t2s1off_ + id + aux];
                if ( j == k ) {
                    dum += u[d1off_ + i];
                    dum -= u[d2s2off_ + i * amo_ + k];
                    dum -= u[d2s2off_ + k * amo_ + i];
                }
                dum += 0.5 * u[d2s0off_ + j * amo_ + k];
                dum += 0.5 * u[d2s0off_ + k * amo_ + j];

                A[offset + id] = dum;
            }
        }
        // sixth block: Pi(ik)
        for (int k = 0; k < amo_; k++) {
            if ( i == k ) continue;
            int myk = k;
            if ( k > i ) myk--;
            int id = (myk + amo_ - 1) * ( 2 * amo_ - 1 ) + 2 * amo_ - 2;

            double dum = -u[t2s1off_ + id + aux];
            dum += 0.5 * u[d2s0off_ + i * amo_ + k];
            dum += 0.5 * u[d2s0off_ + k * amo_ + i];

            A[offset + id] = dum;
        }
        // seventh block: D(ji)
        for (int j = 0; j < amo_; j++) {
            if ( i == j ) continue;
            int myj = j;
            if ( j > i ) myj--;
            int id = (2 * amo_ - 2) * ( 2 * amo_ - 1 ) + myj;

            double dum = -u[t2s1off_ + id + aux];
            dum += 0.5 * u[d2s2off_ + j * amo_ + i];
            dum += 0.5 * u[d2s2off_ + i * amo_ + j];

            A[offset + id] = dum;
        }
        // eighth block: Pi(ji)
        for (int j = 0; j < amo_; j++) {
            if ( i == j ) continue;
            int myj = j;
            if ( j > i ) myj--;
            int id = (2 * amo_ - 2) * ( 2 * amo_ - 1 ) + myj + amo_ - 1;

            double dum = -u[t2s1off_ + id + aux];
            dum += 0.5 * u[d2s0off_ + j * amo_ + i];
            dum += 0.5 * u[d2s0off_ + i * amo_ + j];

            A[offset + id] = dum;
        }
        // ninth block: p(i)
        int id = (2 * amo_ - 2) * ( 2 * amo_ - 1 ) + 2 * amo_ - 2;

        double dum = -u[t2s1off_ + id + aux];
        dum += u[d1off_ + i];

        A[offset + id] = dum;


        offset += ( 2 * amo_ - 1 ) * ( 2 * amo_ - 1 );
        aux    += ( 2 * amo_ - 1 ) * ( 2 * amo_ - 1 );
    }

/*
    int aux = 0;
    for ( int b = 0; b < amo_; b++) {

        // first block: D(ac) = dac p(c) + D(ac)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            for ( int c = 0; c < amo_; c++) {
                if ( c == b ) continue;
                int myc = c;
                if ( c > b ) myc--;
                int id = mya * ( 2 * amo_ - 1 ) + myc;

                double dum = -u[t2s1off_ + id + aux];
                if ( a == c ) {
                    dum += u[d1off_ + c];
                } else {
                    dum += 0.5 * u[d2s2off_ + a * amo_ + c];
                    dum += 0.5 * u[d2s2off_ + c * amo_ + a];
                }

                A[offset + id] = dum;
            }
        }

        // second block: -dac Pi(ab)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            for ( int c = 0; c < amo_; c++) {
                if ( c == b ) continue;
                int myc = c;
                if ( c > b ) myc--;
                int id = mya * ( 2 * amo_ - 1 ) + myc + amo_ - 1;

                double dum = -u[t2s1off_ + id + aux];
                if ( a == c ) {
                    dum -= 0.5 * u[d2s0off_ + a * amo_ + b];
                    dum -= 0.5 * u[d2s0off_ + b * amo_ + a];
                }

                A[offset + id] = dum;
            }
        }

        // third block: D(ab) 
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            int id = mya * ( 2 * amo_ - 1 ) + 2 * amo_ - 2;

            double dum = -u[t2s1off_ + id + aux];
            dum += 0.5 * u[d2s2off_ + a*amo_ + b];
            dum += 0.5 * u[d2s2off_ + b*amo_ + a];

            A[offset + id] = dum;
        }

        // fourth block: -dac Pi(ab)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            for ( int c = 0; c < amo_; c++) {
                if ( c == b ) continue;
                int myc = c;
                if ( c > b ) myc--;
                int id = (mya + amo_ - 1) * ( 2 * amo_ - 1 ) + myc;

                double dum = -u[t2s1off_ + id + aux];
                if ( a == c ) {
                    dum -= 0.5 * u[d2s0off_ + b * amo_ + a];
                    dum -= 0.5 * u[d2s0off_ + a * amo_ + b];
                }

                A[offset + id] = dum;
            }
        }

        // fifth block:  Pi(ac) + dac ( -D(ba) - D(ab) + p(b) )
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            for ( int c = 0; c < amo_; c++) {
                if ( c == b ) continue;
                int myc = c;
                if ( c > b ) myc--;
                int id = (mya + amo_ - 1) * ( 2 * amo_ - 1 ) + myc + amo_ - 1;

                double dum = -u[t2s1off_ + id + aux];
                if ( a == c ) {
                    dum += u[d1off_ + b];
                    dum -= u[d2s2off_ + b * amo_ + a];
                    dum -= u[d2s2off_ + a * amo_ + b];
                }
                dum += 0.5 * u[d2s0off_ + a * amo_ + c];
                dum += 0.5 * u[d2s0off_ + c * amo_ + a];

                A[offset + id] = dum;
            }
        }

        // sixth block: Pi(ab)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            int id = (mya + amo_ - 1) * ( 2 * amo_ - 1 ) + 2 * amo_ - 2;

            double dum = -u[t2s1off_ + id + aux];
            dum       += 0.5 * u[d2s0off_ + a*amo_ + b];
            dum       += 0.5 * u[d2s0off_ + b*amo_ + a];
            A[offset + id] = dum;
        }

        // seventh block: D(bc)
        for ( int c = 0; c < amo_; c++) {
            if ( c == b ) continue;
            int myc = c;
            if ( c > b ) myc--;
            int id = (2 * amo_ - 2) * ( 2 * amo_ - 1 ) + myc;

            double dum = -u[t2s1off_ + id + aux];
            dum       +=  0.5 * u[d2s2off_ + b*amo_ + c];
            dum       +=  0.5 * u[d2s2off_ + c*amo_ + b];

            A[offset + id] = dum;
        }

        // eigth block: Pi(bc)
        for ( int c = 0; c < amo_; c++) {
            if ( c == b ) continue;
            int myc = c;
            if ( c > b ) myc--;
            int id = (2 * amo_ - 2) * ( 2 * amo_ - 1 ) + myc + amo_ - 1;

            double dum = -u[t2s1off_ + id + aux];
            dum       += 0.5 * u[d2s0off_ + b * amo_ + c];
            dum       += 0.5 * u[d2s0off_ + c * amo_ + b];

            A[offset + id] = dum;
        }

        // ninth block: p(b)
        int id = (2 * amo_ - 2) * ( 2 * amo_ - 1 ) + 2 * amo_ - 2;

        double dum = -u[t2s1off_ + id + aux];
        dum       +=  u[d1off_ + b];

        A[offset + id] = dum;

        offset += ( 2 * amo_ - 1 ) * ( 2 * amo_ - 1 );
        aux    += ( 2 * amo_ - 1 ) * ( 2 * amo_ - 1 );

    }
*/
/*
    // Poelmans' thesis, Eq (2.208) - note, AED symmetrized all D2/Pi terms
    //
    //
    //  dac p(c) + D(ac) | dac Pi(ab)                      | p(a) 
    //  dac Pi(bc)       | dac ( p(b) - 2 D(bc) ) - Pi(ac) | 0
    //  p(c)             | 0                               | 1
    //
    //

    int aux = 0;
    for ( int b = 0; b < amo_; b++) {

        // first block: dac p(c) + D(ac)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            for ( int c = 0; c < amo_; c++) {
                if ( c == b ) continue;
                int myc = c;
                if ( c > b ) myc--;
                int id = mya * ( 2 * amo_ - 1 ) + myc;

                double dum = -u[t2s1off_ + id + aux];
                if ( a == c ) {
                    dum += u[d1off_ + c];
                } else {
                    dum += 0.5 * u[d2s2off_ + a * amo_ + c];
                    dum += 0.5 * u[d2s2off_ + c * amo_ + a];
                }

                A[offset + id] = dum;
            }
        }

        // second block: dac Pi(ab)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            for ( int c = 0; c < amo_; c++) {
                if ( c == b ) continue;
                int myc = c;
                if ( c > b ) myc--;
                int id = mya * ( 2 * amo_ - 1 ) + myc + amo_ - 1;

                double dum = -u[t2s1off_ + id + aux];
                if ( a == c ) {
                    dum += 0.5 * u[d2s0off_ + a * amo_ + b];
                    dum += 0.5 * u[d2s0off_ + b * amo_ + a];
                }

                A[offset + id] = dum;
            }
        }

        // third block: p(a)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            int id = mya * ( 2 * amo_ - 1 ) + 2 * amo_ - 2;

            double dum = -u[t2s1off_ + id + aux];
            dum += u[d1off_ + a];

            A[offset + id] = dum;
        }

        // fourth block: dac Pi(bc)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            for ( int c = 0; c < amo_; c++) {
                if ( c == b ) continue;
                int myc = c;
                if ( c > b ) myc--;
                int id = (mya + amo_ - 1) * ( 2 * amo_ - 1 ) + myc;

                double dum = -u[t2s1off_ + id + aux];
                if ( a == c ) {
                    dum += 0.5 * u[d2s0off_ + b * amo_ + c];
                    dum += 0.5 * u[d2s0off_ + c * amo_ + b];
                }

                A[offset + id] = dum;
            }
        }

        // fifth block: dac ( p(b) - 2 D(bc) ) - Pi(ac)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            for ( int c = 0; c < amo_; c++) {
                if ( c == b ) continue;
                int myc = c;
                if ( c > b ) myc--;
                int id = (mya + amo_ - 1) * ( 2 * amo_ - 1 ) + myc + amo_ - 1;

                double dum = -u[t2s1off_ + id + aux];
                if ( a == c ) {
                    dum += u[d1off_ + b];
                    dum -= u[d2s2off_ + b * amo_ + c];
                    dum -= u[d2s2off_ + c * amo_ + b];
                }
                dum -= 0.5 * u[d2s0off_ + a * amo_ + c];
                dum -= 0.5 * u[d2s0off_ + c * amo_ + a];

                A[offset + id] = dum;
            }
        }

        // sixth block: 0
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            int id = (mya + amo_ - 1) * ( 2 * amo_ - 1 ) + 2 * amo_ - 2;

            double dum = -u[t2s1off_ + id + aux];
            A[offset + id] = dum;
        }

        // seventh block: p(c)
        for ( int c = 0; c < amo_; c++) {
            if ( c == b ) continue;
            int myc = c;
            if ( c > b ) myc--;
            int id = (2 * amo_ - 2) * ( 2 * amo_ - 1 ) + myc;

            double dum = -u[t2s1off_ + id + aux];
            dum       +=  u[d1off_ + c];

            A[offset + id] = dum;
        }

        // eigth block: 0
        for ( int c = 0; c < amo_; c++) {
            if ( c == b ) continue;
            int myc = c;
            if ( c > b ) myc--;
            int id = (2 * amo_ - 2) * ( 2 * amo_ - 1 ) + myc + amo_ - 1;

            double dum = -u[t2s1off_ + id + aux];

            A[offset + id] = dum;
        }

        // ninth block: 1
        int id = (2 * amo_ - 2) * ( 2 * amo_ - 1 ) + 2 * amo_ - 2;

        double dum = -u[t2s1off_ + id + aux];
        A[offset + id] = dum;

        offset += ( 2 * amo_ - 1 ) * ( 2 * amo_ - 1 );
        aux    += ( 2 * amo_ - 1 ) * ( 2 * amo_ - 1 );

    }
*/

    // jcp: first T2 constraint in appendix
    // note: AED symmetrized terms

    aux = 0;
    for (int i = 0; i < amo_; i++) {
        for (int k = i + 1; k < amo_; k++) {
            for (int m = k + 1; m < amo_; m++) {

                // 0,0

                double dum = -u[t2off_ + aux + 0];      // - T2(ikm;0,0)

                dum += 0.5 * u[d2s2off_ + i*amo_+ k];   // + 0.5 * D2s2(ik)
                dum += 0.5 * u[d2s2off_ + k*amo_+ i];   // + 0.5 * D2s2(ki)

                dum -= 0.5 * u[d2s2off_ + i*amo_+ m];   // - 0.5 * D2s2(im)
                dum -= 0.5 * u[d2s2off_ + m*amo_+ i];   // - 0.5 * D2s2(mi)

                dum -= 0.5 * u[d2s2off_ + k*amo_+ m];   // - 0.5 * D2s2(km)
                dum -= 0.5 * u[d2s2off_ + m*amo_+ k];   // - 0.5 * D2s2(mk)

                dum += u[d1off_ + m];                   // + D1(m)

                A[offset + 0] = dum;

                // 0,1

                dum = -u[t2off_ + aux + 1];      // - T2(ikm;0,1)

                dum += 0.5 * u[d2s0off_ + m*amo_+ k];   // + 0.5 * D2s0(mk)
                dum += 0.5 * u[d2s0off_ + k*amo_+ m];   // + 0.5 * D2s0(km)

                A[offset + 1] = dum;

                // 0,2

                dum = -u[t2off_ + aux + 2];      // - T2(ikm;0,2)

                dum += 0.5 * u[d2s0off_ + m*amo_+ i];   // + 0.5 * D2s0(mi)
                dum += 0.5 * u[d2s0off_ + i*amo_+ m];   // + 0.5 * D2s0(im)

                A[offset + 2] = dum;

                // 1,0

                dum = -u[t2off_ + aux + 3];      // - T2(ikm;1,0)

                dum += 0.5 * u[d2s0off_ + k*amo_+ m];   // + 0.5 * D2s0(km)
                dum += 0.5 * u[d2s0off_ + m*amo_+ k];   // + 0.5 * D2s0(mk)

                A[offset + 3] = dum;

                // 1,1

                dum = -u[t2off_ + aux + 4];      // - T2(ikm;1,1)

                dum -= 0.5 * u[d2s2off_ + i*amo_+ k];   // - 0.5 * D2s2(ik)
                dum -= 0.5 * u[d2s2off_ + k*amo_+ i];   // - 0.5 * D2s2(ki)

                dum += 0.5 * u[d2s2off_ + i*amo_+ m];   // + 0.5 * D2s2(im)
                dum += 0.5 * u[d2s2off_ + m*amo_+ i];   // + 0.5 * D2s2(mi)

                dum -= 0.5 * u[d2s2off_ + m*amo_+ k];   // - 0.5 * D2s2(mk)
                dum -= 0.5 * u[d2s2off_ + k*amo_+ m];   // - 0.5 * D2s2(km)

                dum += u[d1off_ + k];                   // + D1(m)

                A[offset + 4] = dum;

                // 1,2

                dum = -u[t2off_ + aux + 5];      // - T2(ikm;1,2)

                dum += 0.5 * u[d2s0off_ + k*amo_+ i];   // + 0.5 * D2s0(ki)
                dum += 0.5 * u[d2s0off_ + i*amo_+ k];   // + 0.5 * D2s0(ik)

                A[offset + 5] = dum;

                // 2,0

                dum = -u[t2off_ + aux + 6];      // - T2(ikm;2,0)

                dum += 0.5 * u[d2s0off_ + i*amo_+ m];   // + 0.5 * D2s0(im)
                dum += 0.5 * u[d2s0off_ + m*amo_+ i];   // + 0.5 * D2s0(mi)

                A[offset + 6] = dum;

                // 2,1

                dum = -u[t2off_ + aux + 7];      // - T2(ikm;2,1)

                dum += 0.5 * u[d2s0off_ + i*amo_+ k];   // + 0.5 * D2s0(ik)
                dum += 0.5 * u[d2s0off_ + k*amo_+ i];   // + 0.5 * D2s0(ki)

                A[offset + 7] = dum;

                // 2,2

                dum = -u[t2off_ + aux + 8];      // - T2(ikm;2,2)

                dum -= 0.5 * u[d2s2off_ + i*amo_+ k];   // - 0.5 * D2s2(ik)
                dum -= 0.5 * u[d2s2off_ + k*amo_+ i];   // - 0.5 * D2s2(ki)

                dum -= 0.5 * u[d2s2off_ + i*amo_+ m];   // - 0.5 * D2s2(im)
                dum -= 0.5 * u[d2s2off_ + m*amo_+ i];   // - 0.5 * D2s2(mi)

                dum += 0.5 * u[d2s2off_ + k*amo_+ m];   // + 0.5 * D2s2(km)
                dum += 0.5 * u[d2s2off_ + m*amo_+ k];   // + 0.5 * D2s2(mk)

                dum += u[d1off_ + i];                   // + D1(m)

                A[offset + 8] = dum;

                aux    += 9;
                offset += 9;

            }
        }
    }

}

// T2 portion of A^T.y (with symmetry)
void v2RDM_DOCISolver::T2_constraints_ATu(double* A,double* u){

    // Rubio-Garcia et. al JCTC (2018): Eq. 33
    //
    //
    //  D(ac)            | -dac Pi(ab)                            | D(ab)
    //  -dac Pi(ba)      | Pi(ac) + dac ( -D(ba) - D(ab) + p(b) ) | Pi(ab)
    //  D(bc)            | Pi(bc)                                 | p(b)
    //
    //  (for all b; no restriction on a/c ... but JCP has restriction a!=b and c!=b)

    // 
    //  D(jk)            | -djk Pi(ik)                            | D(ik)
    //  -djk Pi(ki)      | Pi(jk) + djk ( -D(ik) - D(ki) + p(i) ) | Pi(ik)
    //  D(ji)            | Pi(ji)                                 | p(i)
    // 
    int aux = 0;
    for (int i = 0; i < amo_; i++) {
        // first block: D(jk)
        for (int j = 0; j < amo_; j++) {
            if ( i == j ) continue;
            int myj = j;
            if ( j > i ) myj--;
            for (int k = 0; k < amo_; k++) {
                if ( i == k ) continue;
                int myk = k;
                if ( k > i ) myk--;
                int id = myk * ( 2 * amo_ - 1 ) + myj;

                double dum = u[offset + id];

                A[t2s1off_ + id + aux] -= dum;
                if ( j == k ) {
                    A[d1off_ + j] += dum;
                }else {
                    A[d2s2off_ + j * amo_ + k] += 0.5 * dum;
                    A[d2s2off_ + k * amo_ + j] += 0.5 * dum;
                }

            }
        }
        // second block: -djk Pi(ik) 
        for (int j = 0; j < amo_; j++) {
            if ( i == j ) continue;
            int myj = j;
            if ( j > i ) myj--;
            for (int k = 0; k < amo_; k++) {
                if ( i == k ) continue;
                int myk = k;
                if ( k > i ) myk--;
                int id = myk * ( 2 * amo_ - 1 ) + myj + amo_ - 1;

                double dum = u[offset + id];

                A[t2s1off_ + id + aux] -= dum;
                if ( j == k ) {
                    A[d2s0off_ + i * amo_ + k] -= 0.5 * dum;
                    A[d2s0off_ + k * amo_ + i] -= 0.5 * dum;
                }

            }
        }
        // third block: D(ik)
        for (int k = 0; k < amo_; k++) {
            if ( i == k ) continue;
            int myk = k;
            if ( k > i ) myk--;
            int id = myk * ( 2 * amo_ - 1 ) + 2 * amo_ - 2;

            double dum = u[offset + id];

            A[t2s1off_ + id + aux] -= dum;
            A[d2s2off_ + i * amo_ + k] += 0.5 * dum;
            A[d2s2off_ + k * amo_ + i] += 0.5 * dum;

        }
        // fourth block: -djk Pi(ki) 
        for (int j = 0; j < amo_; j++) {
            if ( i == j ) continue;
            int myj = j;
            if ( j > i ) myj--;
            for (int k = 0; k < amo_; k++) {
                if ( i == k ) continue;
                int myk = k;
                if ( k > i ) myk--;
                int id = (myk + amo_ - 1) * ( 2 * amo_ - 1 ) + myj;

                double dum = u[offset + id];

                A[t2s1off_ + id + aux] -= dum;
                if ( j == k ) {
                    A[d2s0off_ + i * amo_ + k] -= 0.5 * dum;
                    A[d2s0off_ + k * amo_ + i] -= 0.5 * dum;
                }

            }
        }
        // fifth block: Pi(jk) + djk ( -D(ik) - D(ki) + p(i) )
        for (int j = 0; j < amo_; j++) {
            if ( i == j ) continue;
            int myj = j;
            if ( j > i ) myj--;
            for (int k = 0; k < amo_; k++) {
                if ( i == k ) continue;
                int myk = k;
                if ( k > i ) myk--;
                int id = (myk + amo_ - 1) * ( 2 * amo_ - 1 ) + myj + amo_ - 1;

                double dum = u[offset + id];

                A[t2s1off_ + id + aux] -= dum;
                if ( j == k ) {
                    A[d1off_ + i] += dum;
                    A[d2s2off_ + i * amo_ + k] -= dum;
                    A[d2s2off_ + k * amo_ + i] -= dum;
                }
                A[d2s0off_ + j * amo_ + k] += 0.5 * dum;
                A[d2s0off_ + k * amo_ + j] += 0.5 * dum;

            }
        }
        // sixth block: Pi(ik)
        for (int k = 0; k < amo_; k++) {
            if ( i == k ) continue;
            int myk = k;
            if ( k > i ) myk--;
            int id = (myk + amo_ - 1) * ( 2 * amo_ - 1 ) + 2 * amo_ - 2;

            double dum = u[offset + id];

            A[t2s1off_ + id + aux] -= dum;
            A[d2s0off_ + i * amo_ + k] += 0.5 * dum;
            A[d2s0off_ + k * amo_ + i] += 0.5 * dum;

        }
        // seventh block: D(ji)
        for (int j = 0; j < amo_; j++) {
            if ( i == j ) continue;
            int myj = j;
            if ( j > i ) myj--;
            int id = (2 * amo_ - 2) * ( 2 * amo_ - 1 ) + myj;

            double dum = u[offset + id];

            A[t2s1off_ + id + aux] -= dum;
            A[d2s2off_ + j * amo_ + i] += 0.5 * dum;
            A[d2s2off_ + i * amo_ + j] += 0.5 * dum;

        }
        // eighth block: Pi(ji)
        for (int j = 0; j < amo_; j++) {
            if ( i == j ) continue;
            int myj = j;
            if ( j > i ) myj--;
            int id = (2 * amo_ - 2) * ( 2 * amo_ - 1 ) + myj + amo_ - 1;

            double dum = u[offset + id];

            A[t2s1off_ + id + aux] -= dum;
            A[d2s0off_ + j * amo_ + i] += 0.5 * dum;
            A[d2s0off_ + i * amo_ + j] += 0.5 * dum;

        }
        // ninth block: p(i)
        int id = (2 * amo_ - 2) * ( 2 * amo_ - 1 ) + 2 * amo_ - 2;

        double dum = u[offset + id];

        A[t2s1off_ + id + aux] -= dum;
        A[d1off_ + i] += dum;

        offset += ( 2 * amo_ - 1 ) * ( 2 * amo_ - 1 );
        aux    += ( 2 * amo_ - 1 ) * ( 2 * amo_ - 1 );
    }

/*
    int aux = 0;
    for ( int b = 0; b < amo_; b++) {

        // first block: dac p(c) + D(ac)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            for ( int c = 0; c < amo_; c++) {
                if ( c == b ) continue;
                int myc = c;
                if ( c > b ) myc--;
                int id = mya * ( 2 * amo_ - 1 ) + myc;

                double dum = u[offset + id];

                A[t2s1off_ + id + aux] -= dum;
                if ( a == c ) {
                    A[d1off_ + c] += dum;
                } else {
                    A[d2s2off_ + a * amo_ + c] += 0.5 * dum;
                    A[d2s2off_ + c * amo_ + a] += 0.5 * dum;
                }
            }
        }

        // second block: -dac Pi(ab)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            for ( int c = 0; c < amo_; c++) {
                if ( c == b ) continue;
                int myc = c;
                if ( c > b ) myc--;
                int id = mya * ( 2 * amo_ - 1 ) + myc + amo_ - 1;

                double dum = u[offset + id];

                A[t2s1off_ + id + aux] -= dum;;
                if ( a == c ) {
                    A[d2s0off_ + a * amo_ + b] -= 0.5 * dum;
                    A[d2s0off_ + b * amo_ + a] -= 0.5 * dum;
                }
            }
        }

        // third block: D(ab)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            int id = mya * ( 2 * amo_ - 1 ) + 2 * amo_ - 2;

            double dum = u[offset + id];

            A[t2s1off_ + id + aux] -= dum;
            A[d2s2off_ + a * amo_ + b] += 0.5 * dum;
            A[d2s2off_ + b * amo_ + a] += 0.5 * dum;
        }

        // fourth block: -dac Pi(ba)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            for ( int c = 0; c < amo_; c++) {
                if ( c == b ) continue;
                int myc = c;
                if ( c > b ) myc--;
                int id = (mya + amo_ - 1) * ( 2 * amo_ - 1 ) + myc;

                double dum = u[offset + id];

                A[t2s1off_ + id + aux] -= dum;
                if ( a == c ) {
                    A[d2s0off_ + b * amo_ + a] -= 0.5 * dum;
                    A[d2s0off_ + a * amo_ + b] -= 0.5 * dum;
                }
            }
        }

        // fifth block: Pi(ac) + dac ( -D(ba) - D(ab) + p(b) )
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            for ( int c = 0; c < amo_; c++) {
                if ( c == b ) continue;
                int myc = c;
                if ( c > b ) myc--;
                int id = (mya + amo_ - 1) * ( 2 * amo_ - 1 ) + myc + amo_ - 1;

                double dum = u[offset + id];

                A[t2s1off_ + id + aux] -= dum;
                if ( a == c ) {
                    A[d1off_ + b]              += dum;
                    A[d2s2off_ + b * amo_ + a] -= dum;
                    A[d2s2off_ + a * amo_ + b] -= dum;
                }
                A[d2s0off_ + a * amo_ + c] += 0.5 * dum;
                A[d2s0off_ + c * amo_ + a] += 0.5 * dum;
            }
        }

        // sixth block: Pi(ab)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            int id = (mya + amo_ - 1) * ( 2 * amo_ - 1 ) + 2 * amo_ - 2;

            double dum = u[offset + id];

            A[t2s1off_ + id + aux] -= dum;
            A[d2s0off_ + a * amo_ + b] += 0.5 * dum;
            A[d2s0off_ + b * amo_ + a] += 0.5 * dum;
        }

        // seventh block: D(bc)
        for ( int c = 0; c < amo_; c++) {
            if ( c == b ) continue;
            int myc = c;
            if ( c > b ) myc--;
            int id = (2 * amo_ - 2) * ( 2 * amo_ - 1 ) + myc;

            double dum = u[offset + id];

            A[t2s1off_ + id + aux] -= dum;
            A[d2s2off_ + b * amo_ + c] += 0.5 * dum;
            A[d2s2off_ + c * amo_ + b] += 0.5 * dum;
        }

        // eigth block: Pi(bc)
        for ( int c = 0; c < amo_; c++) {
            if ( c == b ) continue;
            int myc = c;
            if ( c > b ) myc--;
            int id = (2 * amo_ - 2) * ( 2 * amo_ - 1 ) + myc + amo_ - 1;

            double dum = u[offset + id];

            A[t2s1off_ + id + aux] -= dum;
            A[d2s0off_ + b * amo_ + c] -= 0.5 * dum;
            A[d2s0off_ + c * amo_ + b] -= 0.5 * dum;
        }

        // ninth block: p(b)
        int id = (2 * amo_ - 2) * ( 2 * amo_ - 1 ) + 2 * amo_ - 2;

        double dum = u[offset + id];

        A[t2s1off_ + id + aux] -= dum;
        A[d1off_ + b]          += dum;

        offset += ( 2 * amo_ - 1 ) * ( 2 * amo_ - 1 );
        aux    += ( 2 * amo_ - 1 ) * ( 2 * amo_ - 1 );

    }
*/

    // Poelmans' thesis, Eq (2.208) - note, AED symmetrized all D2/Pi terms
    //
    //
    //  dac p(c) + D(ac) | dac Pi(ab)                      | p(a) 
    //  dac Pi(bc)       | dac ( p(b) - 2 D(bc) ) - Pi(ac) | 0
    //  p(c)             | 0                               | 1
    //
    //

/*
    int aux = 0;
    for ( int b = 0; b < amo_; b++) {

        // first block: dac p(c) + D(ac)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            for ( int c = 0; c < amo_; c++) {
                if ( c == b ) continue;
                int myc = c;
                if ( c > b ) myc--;
                int id = mya * ( 2 * amo_ - 1 ) + myc;

                double dum = u[offset + id];

                A[t2s1off_ + id + aux] -= dum;
                if ( a == c ) {
                    A[d1off_ + c] += dum;
                } else {
                    A[d2s2off_ + a * amo_ + c] += 0.5 * dum;
                    A[d2s2off_ + c * amo_ + a] += 0.5 * dum;
                }
            }
        }

        // second block: dac Pi(ab)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            for ( int c = 0; c < amo_; c++) {
                if ( c == b ) continue;
                int myc = c;
                if ( c > b ) myc--;
                int id = mya * ( 2 * amo_ - 1 ) + myc + amo_ - 1;

                double dum = u[offset + id];

                A[t2s1off_ + id + aux] -= dum;;
                if ( a == c ) {
                    A[d2s0off_ + a * amo_ + b] += 0.5 * dum;
                    A[d2s0off_ + b * amo_ + a] += 0.5 * dum;
                }
            }
        }

        // third block: p(a)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            int id = mya * ( 2 * amo_ - 1 ) + 2 * amo_ - 2;

            double dum = u[offset + id];

            A[t2s1off_ + id + aux] -= dum;
            A[d1off_ + a] += dum;
        }

        // fourth block: dac Pi(bc)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            for ( int c = 0; c < amo_; c++) {
                if ( c == b ) continue;
                int myc = c;
                if ( c > b ) myc--;
                int id = (mya + amo_ - 1) * ( 2 * amo_ - 1 ) + myc;

                double dum = u[offset + id];

                A[t2s1off_ + id + aux] -= dum;
                if ( a == c ) {
                    A[d2s0off_ + b * amo_ + c] += 0.5 * dum;
                    A[d2s0off_ + c * amo_ + b] += 0.5 * dum;
                }
            }
        }

        // fifth block: dac ( p(b) - 2 D(bc) ) - Pi(ac)
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            for ( int c = 0; c < amo_; c++) {
                if ( c == b ) continue;
                int myc = c;
                if ( c > b ) myc--;
                int id = (mya + amo_ - 1) * ( 2 * amo_ - 1 ) + myc + amo_ - 1;

                double dum = u[offset + id];

                A[t2s1off_ + id + aux] -= dum;
                if ( a == c ) {
                    A[d1off_ + b] += dum;
                    A[d2s2off_ + b * amo_ + c] -= dum;
                    A[d2s2off_ + c * amo_ + b] -= dum;
                }
                A[d2s0off_ + a * amo_ + c] -= 0.5 * dum;
                A[d2s0off_ + c * amo_ + a] -= 0.5 * dum;
            }
        }

        // sixth block: 0
        for ( int a = 0; a < amo_; a++) {
            if ( a == b ) continue;
            int mya = a;
            if ( a > b ) mya--;
            int id = (mya + amo_ - 1) * ( 2 * amo_ - 1 ) + 2 * amo_ - 2;

            double dum = u[offset + id];

            A[t2s1off_ + id + aux] -= dum;
        }

        // seventh block: p(c)
        for ( int c = 0; c < amo_; c++) {
            if ( c == b ) continue;
            int myc = c;
            if ( c > b ) myc--;
            int id = (2 * amo_ - 2) * ( 2 * amo_ - 1 ) + myc;

            double dum = u[offset + id];

            A[t2s1off_ + id + aux] -= dum;
            A[d1off_ + c]    += dum;
        }

        // eigth block: 0
        for ( int c = 0; c < amo_; c++) {
            if ( c == b ) continue;
            int myc = c;
            if ( c > b ) myc--;
            int id = (2 * amo_ - 2) * ( 2 * amo_ - 1 ) + myc + amo_ - 1;

            double dum = u[offset + id];

            A[t2s1off_ + id + aux] -= dum;
        }

        // ninth block: 1
        int id = (2 * amo_ - 2) * ( 2 * amo_ - 1 ) + 2 * amo_ - 2;

        double dum = u[offset + id];

        A[t2s1off_ + id + aux] -= dum;

        offset += ( 2 * amo_ - 1 ) * ( 2 * amo_ - 1 );
        aux    += ( 2 * amo_ - 1 ) * ( 2 * amo_ - 1 );

    }
*/

    // jcp: first T2 constraint in appendix
    // note: AED symmetrized terms

    aux = 0;
    for (int i = 0; i < amo_; i++) {
        for (int k = i + 1; k < amo_; k++) {
            for (int m = k + 1; m < amo_; m++) {

                // 0,0

                double dum = u[offset + 0];

                A[t2off_ + aux + 0] -= dum;      // - T2(ikm;0,0)

                A[d2s2off_ + i*amo_+ k] += 0.5 * dum;   // + 0.5 * D2s2(ik)
                A[d2s2off_ + k*amo_+ i] += 0.5 * dum;   // + 0.5 * D2s2(ki)

                A[d2s2off_ + i*amo_+ m] -= 0.5 * dum;   // - 0.5 * D2s2(im)
                A[d2s2off_ + m*amo_+ i] -= 0.5 * dum;   // - 0.5 * D2s2(mi)

                A[d2s2off_ + k*amo_+ m] -= 0.5 * dum;   // - 0.5 * D2s2(km)
                A[d2s2off_ + m*amo_+ k] -= 0.5 * dum;   // - 0.5 * D2s2(mk)

                A[d1off_ + m] += dum;                   // + D1(m)

                // 0,1

                dum = u[offset + 1];

                A[t2off_ + aux + 1] -= dum;      // - T2(ikm;0,1)

                A[d2s0off_ + m*amo_+ k] += 0.5 * dum;   // + 0.5 * D2s0(mk)
                A[d2s0off_ + k*amo_+ m] += 0.5 * dum;   // + 0.5 * D2s0(km)

                // 0,2

                dum = u[offset + 2];

                A[t2off_ + aux + 2] -= dum;      // - T2(ikm;0,2)

                A[d2s0off_ + m*amo_+ i] += 0.5 * dum;   // + 0.5 * D2s0(mi)
                A[d2s0off_ + i*amo_+ m] += 0.5 * dum;   // + 0.5 * D2s0(im)

                // 1,0

                dum = u[offset + 3];      // - T2(ikm;1,0)

                A[t2off_ + aux + 3] -= dum;      // - T2(ikm;1,0)

                A[d2s0off_ + k*amo_+ m] += 0.5 * dum;   // + 0.5 * D2s0(km)
                A[d2s0off_ + m*amo_+ k] += 0.5 * dum;   // + 0.5 * D2s0(mk)

                // 1,1

                dum = u[offset + 4];

                A[t2off_ + aux + 4] -= dum;      // - T2(ikm;1,1)

                A[d2s2off_ + i*amo_+ k] -= 0.5 * dum;   // - 0.5 * D2s2(ik)
                A[d2s2off_ + k*amo_+ i] -= 0.5 * dum;   // - 0.5 * D2s2(ki)

                A[d2s2off_ + i*amo_+ m] += 0.5 * dum;   // + 0.5 * D2s2(im)
                A[d2s2off_ + m*amo_+ i] += 0.5 * dum;   // + 0.5 * D2s2(mi)

                A[d2s2off_ + k*amo_+ m] -= 0.5 * dum;   // - 0.5 * D2s2(km)
                A[d2s2off_ + m*amo_+ k] -= 0.5 * dum;   // - 0.5 * D2s2(mk)

                A[d1off_ + k] += dum;                   // + D1(m)

                // 1,2

                dum = u[offset + 5];

                A[t2off_ + aux + 5] -= dum;      // - T2(ikm;1,2)

                A[d2s0off_ + k*amo_+ i] += 0.5 * dum;   // + 0.5 * D2s0(ki)
                A[d2s0off_ + i*amo_+ k] += 0.5 * dum;   // + 0.5 * D2s0(ik)

                // 2,0

                dum = u[offset + 6];

                A[t2off_ + aux + 6] -= dum;      // - T2(ikm;2,0)

                A[d2s0off_ + i*amo_+ m] += 0.5 * dum;   // + 0.5 * D2s0(im)
                A[d2s0off_ + m*amo_+ i] += 0.5 * dum;   // + 0.5 * D2s0(mi)

                // 2,1

                dum = u[offset + 7];      // - T2(ikm;2,1)

                A[t2off_ + aux + 7] -= dum;      // - T2(ikm;2,0)

                A[d2s0off_ + i*amo_+ k] += 0.5 * dum;   // + 0.5 * D2s0(ik)
                A[d2s0off_ + k*amo_+ i] += 0.5 * dum;   // + 0.5 * D2s0(ki)

                // 2,2

                dum = u[offset + 8];      // - T2(ikm;2,2)

                A[t2off_ + aux + 8] -= dum;      // - T2(ikm;2,2)

                A[d2s2off_ + i*amo_+ k] -= 0.5 * dum;   // - 0.5 * D2s2(ik)
                A[d2s2off_ + k*amo_+ i] -= 0.5 * dum;   // - 0.5 * D2s2(ki)

                A[d2s2off_ + i*amo_+ m] -= 0.5 * dum;   // - 0.5 * D2s2(im)
                A[d2s2off_ + m*amo_+ i] -= 0.5 * dum;   // - 0.5 * D2s2(mi)

                A[d2s2off_ + k*amo_+ m] += 0.5 * dum;   // + 0.5 * D2s2(km)
                A[d2s2off_ + m*amo_+ k] += 0.5 * dum;   // + 0.5 * D2s2(mk)

                A[d1off_ + i] += dum;                   // + D1(m)

                aux    += 9;
                offset += 9;

            }
        }
    }

}

} // end namespaces
