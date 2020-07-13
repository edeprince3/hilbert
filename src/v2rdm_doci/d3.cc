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

#include "v2rdm_solver.h"

using namespace psi;

namespace psi{ namespace v2rdm_doci{


// D3 portion of A^T.y 
void v2RDMSolver::D3_constraints_ATu(SharedVector A,SharedVector u){
    double* A_p = A->pointer();
    double* u_p = u->pointer();

    // seniority 3: sum_c Dabc = (N/2 - 2) Dab
    // note: AED symmetrized
    for (int a = 0; a < amo_; a++) {
        for (int b = 0; b < amo_; b++) {
            if ( a == b ) continue;
            double dum = u_p[offset];
            for (int c = 0; c < amo_; c++) {
                if ( a == c || b == c) continue;
                A_p[d3s3off_ + c * amo_ * amo_ + a * amo_ + b] += 1.0 / 6.0 * dum;
                A_p[d3s3off_ + c * amo_ * amo_ + b * amo_ + a] += 1.0 / 6.0 * dum;
                A_p[d3s3off_ + a * amo_ * amo_ + c * amo_ + b] += 1.0 / 6.0 * dum;
                A_p[d3s3off_ + b * amo_ * amo_ + c * amo_ + a] += 1.0 / 6.0 * dum;
                A_p[d3s3off_ + a * amo_ * amo_ + b * amo_ + c] += 1.0 / 6.0 * dum;
                A_p[d3s3off_ + b * amo_ * amo_ + a * amo_ + c] += 1.0 / 6.0 * dum;
            }
            A_p[d2s2off_ + a * amo_ + b] -= ( 0.5 * (nalpha_ + nbeta_) - 2.0 ) * dum;
            offset++;
        }
    }

    // seniority 1: Pi^b_aa = Pi^a_bb = Dab = Dba
    for (int a = 0; a < amo_; a++) {
        for (int b = 0; b < amo_; b++) {
            if ( a == b ) continue;

            int myb = b;
            if ( b > a ) myb--;
            //int mya = a;
            //if ( a > b ) mya--;

            double dum = u_p[offset];

            A_p[d3s1off_ + a * (amo_-1) * (amo_-1) + myb * (amo_-1) + myb] += dum;
            //A_p[d3s1off_ + b * (amo_-1) * (amo_-1) + mya * (amo_-1) + mya] += dum;
            A_p[d2s2off_ + a * amo_ + b] -= 0.5 * dum;
            A_p[d2s2off_ + b * amo_ + a] -= 0.5 * dum;

            offset++;

        }
    }

    // seniority 1: Pi^b_aa = Pi^a_bb = Dab = Dba (again)
    for (int a = 0; a < amo_; a++) {
        for (int b = 0; b < amo_; b++) {
            if ( a == b ) continue;

            int mya = a;
            if ( a > b ) mya--;

            double dum = u_p[offset];

            A_p[d3s1off_ + b * (amo_-1) * (amo_-1) + mya * (amo_-1) + mya] += dum;
            A_p[d2s2off_ + a * amo_ + b] -= 0.5 * dum;
            A_p[d2s2off_ + b * amo_ + a] -= 0.5 * dum;

            offset++;

        }
    }

    // seniority 1: sum_b Pi^b_ac = (N/2 - 1) Pi_ac
    for (int a = 0; a < amo_; a++) {
        for (int c = 0; c < amo_; c++) {

            double dum = u_p[offset];

            for (int b = 0; b < amo_; b++) {

                if ( b == a || b == c ) continue;

                int mya = a;
                int myc = c;

                if ( a > b ) mya--;
                if ( c > b ) myc--;

                A_p[d3s1off_ + b * (amo_-1) * (amo_-1) + mya * (amo_-1) + myc] += 0.5 * dum;
                A_p[d3s1off_ + b * (amo_-1) * (amo_-1) + myc * (amo_-1) + mya] += 0.5 * dum;
            }

            A_p[d2s0off_ + a * amo_ + c] -= 0.5 * ( 0.5 * (nalpha_ + nbeta_) - 1.0 ) * dum;
            A_p[d2s0off_ + c * amo_ + a] -= 0.5 * ( 0.5 * (nalpha_ + nbeta_) - 1.0 ) * dum;

            offset++;
        }
    }

}

// D3 portion of A.x (and D1/Q1) 
void v2RDMSolver::D3_constraints_Au(SharedVector A,SharedVector u){

    double* A_p = A->pointer();
    double* u_p = u->pointer();

    // seniority 3: sum_c Dabc = (N/2 - 2) Dab
    // note: AED symmetrized
    for (int a = 0; a < amo_; a++) {
        for (int b = 0; b < amo_; b++) {
            if ( a == b ) continue;
            double dum = 0.0;
            for (int c = 0; c < amo_; c++) {
                if ( a == c || b == c) continue;
                dum += u_p[d3s3off_ + c * amo_ * amo_ + a * amo_ + b];
                dum += u_p[d3s3off_ + c * amo_ * amo_ + b * amo_ + a];
                dum += u_p[d3s3off_ + a * amo_ * amo_ + c * amo_ + b];
                dum += u_p[d3s3off_ + b * amo_ * amo_ + c * amo_ + a];
                dum += u_p[d3s3off_ + a * amo_ * amo_ + b * amo_ + c];
                dum += u_p[d3s3off_ + b * amo_ * amo_ + a * amo_ + c];
            }
            A_p[offset++] = 1.0 / 6.0 * dum - ( 0.5 * (nalpha_ + nbeta_) - 2.0 ) * u_p[d2s2off_ + a * amo_ + b];
        }
    }

    // seniority 1: Pi^b_aa = Pi^a_bb = Dab = Dba
    for (int a = 0; a < amo_; a++) {
        for (int b = 0; b < amo_; b++) {
            if ( a == b ) continue;

            int myb = b;
            if ( b > a ) myb--;
            //int mya = a;
            //if ( a > b ) mya--;

            double dum = 0.0;

            dum += u_p[d3s1off_ + a * (amo_-1) * (amo_-1) + myb * (amo_-1) + myb];
            //dum += u_p[d3s1off_ + b * (amo_-1) * (amo_-1) + mya * (amo_-1) + mya];
            dum -= 0.5 * u_p[d2s2off_ + a * amo_ + b];
            dum -= 0.5 * u_p[d2s2off_ + b * amo_ + a];

            A_p[offset++] = dum;
        }
    }

    // seniority 1: Pi^b_aa = Pi^a_bb = Dab = Dba (again)
    for (int a = 0; a < amo_; a++) {
        for (int b = 0; b < amo_; b++) {
            if ( a == b ) continue;

            int mya = a;
            if ( a > b ) mya--;

            double dum = 0.0;

            dum += u_p[d3s1off_ + b * (amo_-1) * (amo_-1) + mya * (amo_-1) + mya];
            dum -= 0.5 * u_p[d2s2off_ + a * amo_ + b];
            dum -= 0.5 * u_p[d2s2off_ + b * amo_ + a];

            A_p[offset++] = dum;
        }
    }

    // seniority 1: sum_b Pi^b_ac = (N/2 - 1) Pi_ac
    for (int a = 0; a < amo_; a++) {
        for (int c = 0; c < amo_; c++) {

            double dum = 0.0;

            for (int b = 0; b < amo_; b++) {

                if ( b == a || b == c ) continue;

                int mya = a;
                int myc = c;

                if ( a > b ) mya--;
                if ( c > b ) myc--;

                dum += 0.5 * u_p[d3s1off_ + b * (amo_-1) * (amo_-1) + mya * (amo_-1) + myc];
                dum += 0.5 * u_p[d3s1off_ + b * (amo_-1) * (amo_-1) + myc * (amo_-1) + mya];
            }

            dum -= 0.5 * ( 0.5 * (nalpha_ + nbeta_) - 1.0 ) * u_p[d2s0off_ + a * amo_ + c];
            dum -= 0.5 * ( 0.5 * (nalpha_ + nbeta_) - 1.0 ) * u_p[d2s0off_ + c * amo_ + a];

            A_p[offset++] = dum;
        }
    }

}

}}
