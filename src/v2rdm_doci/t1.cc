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

namespace psi{ namespace v2rdm_doci{

// T1 portion of A.x (with symmetry)
void v2RDMSolver::T1_constraints_Au(SharedVector A,SharedVector u){

    double* A_p = A->pointer();
    double* u_p = u->pointer();

    // seniority 3, aaa
    int abc = 0;
    for (int a = 0; a < amo_; a++) {
        for (int b = a + 1; b < amo_; b++) {
            for (int c = b + 1; c < amo_; c++) {

                double dum = u_p[t1s3aaaoff_ + abc];        // T1(abc)

                dum += u_p[d1off_ + a];                  // + D1(a)
                dum += u_p[d1off_ + b];                  // + D1(b)
                dum += u_p[d1off_ + c];                  // + D1(c)
                dum -= 0.5 * u_p[d2s2off_ + a*amo_ + b]; // - 0.5 * D2s2(ab)
                dum -= 0.5 * u_p[d2s2off_ + b*amo_ + a]; // - 0.5 * D2s2(ba)
                dum -= 0.5 * u_p[d2s2off_ + b*amo_ + c]; // - 0.5 * D2s2(bc)
                dum -= 0.5 * u_p[d2s2off_ + c*amo_ + b]; // - 0.5 * D2s2(cb)
                dum -= 0.5 * u_p[d2s2off_ + a*amo_ + c]; // - 0.5 * D2s2(ac)
                dum -= 0.5 * u_p[d2s2off_ + c*amo_ + a]; // - 0.5 * D2s2(ca)

                A_p[offset++] = dum;

                abc++;

            }
        }
    }

    // seniority 3, aab
    abc = 0;
    for (int a = 0; a < amo_; a++) {
        for (int b = a + 1; b < amo_; b++) {
            for (int c = 0; c < amo_; c++) {

                double dum = u_p[t1s3aaboff_ + abc];        // T1(abc)

                dum += u_p[d1off_ + a];                  // + D1(a)
                dum += u_p[d1off_ + b];                  // + D1(b)
                dum += u_p[d1off_ + c];                  // + D1(c)
                dum -= 0.5 * u_p[d2s2off_ + a*amo_ + b]; // - 0.5 * D2s2(ab)
                dum -= 0.5 * u_p[d2s2off_ + b*amo_ + a]; // - 0.5 * D2s2(ba)
                if ( b != c ) {
                    dum -= 0.5 * u_p[d2s2off_ + b*amo_ + c]; // - 0.5 * D2s2(bc)
                    dum -= 0.5 * u_p[d2s2off_ + c*amo_ + b]; // - 0.5 * D2s2(cb)
                }else {
                    dum -= 0.5 * u_p[d2s0off_ + b*amo_ + c]; // - 0.5 * D2s0(bc)
                    dum -= 0.5 * u_p[d2s0off_ + c*amo_ + b]; // - 0.5 * D2s0(cb)
                }
                if ( a != c ) {
                    dum -= 0.5 * u_p[d2s2off_ + a*amo_ + c]; // - 0.5 * D2s2(ac)
                    dum -= 0.5 * u_p[d2s2off_ + c*amo_ + a]; // - 0.5 * D2s2(ca)
                }else {
                    dum -= 0.5 * u_p[d2s0off_ + a*amo_ + c]; // - 0.5 * D2s0(ac)
                    dum -= 0.5 * u_p[d2s0off_ + c*amo_ + a]; // - 0.5 * D2s0(ca)
                }

                A_p[offset++] = dum;

                abc++;

            }
        }
    }

    // seniority 1
    for (int b = 0; b < amo_; b++) {
        for (int a = 0; a < amo_; a++) {

            if ( a == b ) continue;

            for (int c = 0; c < amo_; c++) {

                if ( c == b ) continue;

                int mya = a;
                int myc = c;

                if ( a > b ) mya--;
                if ( c > b ) myc--;

                double dum = u_p[t1s1off_ + b*(amo_-1)*(amo_-1) + mya*(amo_-1) + myc]; // T1^b_ac

                if ( a == c ) {

                    dum += 2.0 * u_p[d1off_ + a];        // + 2 D1(a)
                    dum +=       u_p[d1off_ + b];        // + D1(b)
                    dum -= u_p[d2s2off_ + a*amo_ + b];   // - D2s2(ab)
                    dum -= u_p[d2s2off_ + b*amo_ + a];   // - D2s2(ba)

                }
                dum -= 0.5 * u_p[d2s0off_ + a*amo_ + c]; // - 0.5 * D2s0(ac)
                dum -= 0.5 * u_p[d2s0off_ + c*amo_ + a]; // - 0.5 * D2s0(ca)

                A_p[offset++] = dum;

            }
        }
    }

}

// T1 portion of A^T.y (with symmetry)
void v2RDMSolver::T1_constraints_ATu(SharedVector A,SharedVector u){

    double* A_p = A->pointer();
    double* u_p = u->pointer();

    // seniority 3, aaa
    int abc = 0;
    for (int a = 0; a < amo_; a++) {
        for (int b = a + 1; b < amo_; b++) {
            for (int c = b + 1; c < amo_; c++) {

                double dum = u_p[offset++];

                A_p[t1s3aaaoff_ + abc] += dum;              // T1(abc)

                A_p[d1off_ + a] += dum;                  // + D1(a)
                A_p[d1off_ + b] += dum;                  // + D1(b)
                A_p[d1off_ + c] += dum;                  // + D1(c)
                A_p[d2s2off_ + a*amo_ + b] -= 0.5 * dum; // - 0.5 * D2s2(ab)
                A_p[d2s2off_ + b*amo_ + a] -= 0.5 * dum; // - 0.5 * D2s2(ba)
                A_p[d2s2off_ + b*amo_ + c] -= 0.5 * dum; // - 0.5 * D2s2(bc)
                A_p[d2s2off_ + c*amo_ + b] -= 0.5 * dum; // - 0.5 * D2s2(cb)
                A_p[d2s2off_ + a*amo_ + c] -= 0.5 * dum; // - 0.5 * D2s2(ac)
                A_p[d2s2off_ + c*amo_ + a] -= 0.5 * dum; // - 0.5 * D2s2(ca)

                abc++;

            }
        }
    }

    // seniority 3, aab
    abc = 0;
    for (int a = 0; a < amo_; a++) {
        for (int b = a + 1; b < amo_; b++) {
            for (int c = 0; c < amo_; c++) {

                double dum = u_p[offset++];

                A_p[t1s3aaboff_ + abc] += dum;        // T1(abc)

                A_p[d1off_ + a] += dum;                  // + D1(a)
                A_p[d1off_ + b] += dum;                  // + D1(b)
                A_p[d1off_ + c] += dum;                  // + D1(c)
                A_p[d2s2off_ + a*amo_ + b] -= 0.5 * dum; // - 0.5 * D2s2(ab)
                A_p[d2s2off_ + b*amo_ + a] -= 0.5 * dum; // - 0.5 * D2s2(ba)
                if ( b != c ) {
                    A_p[d2s2off_ + b*amo_ + c] -= 0.5 * dum; // - 0.5 * D2s2(bc)
                    A_p[d2s2off_ + c*amo_ + b] -= 0.5 * dum; // - 0.5 * D2s2(cb)
                }else {
                    A_p[d2s0off_ + b*amo_ + c] -= 0.5 * dum; // - 0.5 * D2s0(bc)
                    A_p[d2s0off_ + c*amo_ + b] -= 0.5 * dum; // - 0.5 * D2s0(cb)
                }
                if ( a != c ) {
                    A_p[d2s2off_ + a*amo_ + c] -= 0.5 * dum; // - 0.5 * D2s2(ac)
                    A_p[d2s2off_ + c*amo_ + a] -= 0.5 * dum; // - 0.5 * D2s2(ca)
                }else {
                    A_p[d2s0off_ + a*amo_ + c] -= 0.5 * dum; // - 0.5 * D2s0(ac)
                    A_p[d2s0off_ + c*amo_ + a] -= 0.5 * dum; // - 0.5 * D2s0(ca)
                }

                abc++;

            }
        }
    }

    for (int b = 0; b < amo_; b++) {
        for (int a = 0; a < amo_; a++) {

            if ( a == b ) continue;

            for (int c = 0; c < amo_; c++) {

                if ( c == b ) continue;

                int mya = a;
                int myc = c;

                if ( a > b ) mya--;
                if ( c > b ) myc--;

                double dum = u_p[offset++];

                A_p[t1s1off_ + b*(amo_-1)*(amo_-1) + mya*(amo_-1) + myc] += dum; // T1^b_ac

                if ( a == c ) {
                    
                    A_p[d1off_ + a] += 2.0 * dum;        // + 2 D1(a)
                    A_p[d1off_ + b] +=       dum;        // + D1(b)
                    A_p[d2s2off_ + a*amo_ + b] -= dum;   // - D2s2(ab)
                    A_p[d2s2off_ + b*amo_ + a] -= dum;   // - D2s2(ba)
                
                }
                A_p[d2s0off_ + a*amo_ + c] -= 0.5 * dum; // - 0.5 * D2s0(ac)
                A_p[d2s0off_ + c*amo_ + a] -= 0.5 * dum; // - 0.5 * D2s0(ca)

            }
        }
    }
}

}}
