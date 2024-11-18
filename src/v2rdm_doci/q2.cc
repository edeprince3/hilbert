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

// Q2 portion of A.x (with symmetry)
void v2RDM_DOCISolver::Q2_constraints_Au(double* A,double* u){

    // map D2s2 to Q2s2 seniority-2
    //C_DCOPY(amo_*amo_,u + d2s2off_,1,A + offset,1);      // + D2(ij,ij)
    //C_DAXPY(amo_*amo_,-1.0,u + q2s2off_,1,A + offset,1); // - Q2(ij,ij)

    // note: AED symmetrized
    int ij = 0;
    for (int i = 0; i < amo_; i++) {
        for (int j = i + 1; j < amo_; j++) {
            //if ( i == j ) continue;
            double dum = 0.0;

            dum += 0.5 * u[d2s2off_ + i*amo_ + j]; //  D2(ij,ij)
            dum += 0.5 * u[d2s2off_ + j*amo_ + i]; //  D2(ij,ij)
            dum -= u[q2s2off_ + ij];         // -Q2(ij,ij)
            dum -= u[d1off_ + i];            // -D1(i)
            dum -= u[d1off_ + j];            // -D1(j)

            A[offset + ij] = dum;

            ij++;
        }
    }
    offset += amo_*(amo_-1)/2;

    // map D2s0 to Q2s0 seniority-0
    C_DCOPY(amo_*amo_,u + d2s0off_,1,A + offset,1);      // + D2(ii,jj)
    C_DAXPY(amo_*amo_,-1.0,u + q2s0off_,1,A + offset,1); // - Q2(ii,jj)

    #pragma omp parallel for schedule (static)
    for (int i = 0; i < amo_; i++) {
        double dum = 0.0;
        dum        -=  u[d1off_ + i];  // -D1(i) 
        dum        -=  u[d1off_ + i];  // -D1(j) 
        A[offset + i*amo_ + i] += dum;
    }
    offset += amo_*amo_;

}

// Q2 portion of A^T.y (with symmetry)
void v2RDM_DOCISolver::Q2_constraints_ATu(double* A,double* u){

    // map D2s2 to Q2s2
    //C_DAXPY(amo_*amo_, 1.0,u + offset,1,A + d2s2off_,1); // + D2(ij,ij)
    //C_DAXPY(amo_*amo_,-1.0,u + offset,1,A + q2s2off_,1); // - Q2(ii,jj)

    // note: AED symmetrized
    int ij = 0;
    for (int i = 0; i < amo_; i++) {
        for (int j = i + 1; j < amo_; j++) {
            //if ( i == j ) continue;
            A[d2s2off_ + i*amo_ + j] += 0.5 * u[offset + ij]; // + D2(ij,ij)
            A[d2s2off_ + j*amo_ + i] += 0.5 * u[offset + ij]; // + D2(ij,ij)
            A[q2s2off_ + ij]         -= u[offset + ij]; // - Q2(ij,ij)
            A[d1off_ + i]            -= u[offset + ij]; // -D1(i)
            A[d1off_ + j]            -= u[offset + ij]; // -D1(j)

            ij++;
        }
    }
    offset += amo_*(amo_-1)/2;

    // map D2s0 to Q2s0
    C_DAXPY(amo_*amo_, 1.0,u + offset,1,A + d2s0off_,1); // + D2(ii,jj)
    C_DAXPY(amo_*amo_,-1.0,u + offset,1,A + q2s0off_,1); // - Q2(ii,jj)

    for (int i = 0; i < amo_; i++) {
        double val = u[offset + i*amo_ + i];
        A[d1off_  + i] -= val; // -D1(i)
        A[d1off_  + i] -= val; // -D1(i)
    }
    offset += amo_*amo_;

}

} // end namespaces
