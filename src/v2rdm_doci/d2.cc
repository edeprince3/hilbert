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

#include"v2rdm_doci_solver.h"

#include <misc/omp.h>

using namespace psi;

namespace hilbert{

// D2 portion of A^T.y ( and D1 / Q1 ) 
void v2RDM_DOCISolver::D2_constraints_ATu(double* A,double* u){

    // Traces
    // Tr(D2s2) seniority-2
    for (int i = 0; i < amo_; i++){
        for (int j = 0; j < amo_; j++){
            if ( i==j ) continue;
            A[d2s2off_ + i*amo_+j] += u[offset];
        }
    }
    offset++;

    // Tr(D2s0) seniority-0
    for (int i = 0; i < amo_; i++){
        A[d2s0off_ + i*amo_+i] += u[offset];
    }
    offset++;

    // Tr(D1a)
    for (int i = 0; i < amo_; i++){
        A[d1off_ + i] += u[offset];
    }
    offset++;

    int na = nalpha_;

    // contraction: D2s2 -> D1 a seniority-2
    for (int i = 0; i < amo_; i++){
        A[d1off_ + i] += (na-1.0)*u[offset + i];
        for (int j = 0; j < amo_; j++){
            if (i==j) continue; 
            A[d2s2off_ + i*amo_ + j] -= u[offset + i];
        }
    }
    offset += amo_;

    // contraction: D2s2 -> D1 b seniority-2
    for (int i = 0; i < amo_; i++){
        A[d1off_ + i] += (na-1.0)*u[offset + i];
        for (int j = 0; j < amo_; j++){
            if (i==j) continue; 
            A[d2s2off_ + j*amo_ + i] -= u[offset + i];
        }
    }
    offset += amo_;

    //contract D2s0 -> D1 a seniority-0
    for (int i = 0; i < amo_; i++){
        A[d1off_ + i]          += u[offset + i];
        A[d2s0off_ + i*amo_+i] -= u[offset + i];
    }
    offset += amo_;

    // D2_2 symmetric (with zero diagonal)
    for (int i = 0; i < amo_; i++){
        for (int j = 0; j < amo_; j++){
            if ( i == j ) {
                A[d2s2off_ + i*amo_ + i] += u[offset+i*amo_+i];
            }else {
                A[d2s2off_ + i*amo_ + j] += u[offset+i*amo_+j];
                A[d2s2off_ + j*amo_ + i] -= u[offset+i*amo_+j];
            }
        }
    }
    offset += amo_*amo_;

}

// D2 portion of A.x (and D1/Q1) 
void v2RDM_DOCISolver::D2_constraints_Au(double* A,double* u){

    // Traces
    // Tr(D2s2) seniority-2
    double sums2 =0.0;
    for (int i = 0; i < amo_; i++){
        for (int j = 0; j < amo_; j++){
            if ( i==j ) continue;
            sums2 += u[d2s2off_ + i*amo_ + j];
        }
    }
    A[offset] = sums2;
    offset++;

    // Tr(D2s0) seniority-0
    double sums0 =0.0;
    for (int i = 0; i < amo_; i++){
        sums0 += u[d2s0off_ + i*amo_ + i];
    }
    A[offset] = sums0;
    offset++;

    // Tr(D1a)
    double sums1 =0.0;
    for (int i = 0; i < amo_; i++){
        sums1 += u[d1off_ + i];
    }
    A[offset] = sums1;
    offset++;

    int na = nalpha_;

    // contraction: D2s2 -> D1 a seniority-2
    for (int i = 0; i < amo_; i++){
        double sum = (na-1.0)*u[d1off_ + i];
        for (int j = 0; j < amo_; j++){
            if (i==j) continue; 
            sum -= u[d2s2off_ + i*amo_ + j];
        }
        A[offset + i] = sum;
    }
    offset += amo_;

    // contraction: D2s2 -> D1 b seniority-2
    for (int i = 0; i < amo_; i++){
        double sum = (na-1.0)*u[d1off_ + i];
        for (int j = 0; j < amo_; j++){
            if (i==j) continue; 
            sum -= u[d2s2off_ + j*amo_ + i];
        }
        A[offset + i] = sum;
    }
    offset += amo_;

    //contract D2s0 -> D1 a seniority-0
    for (int i = 0; i < amo_; i++){
        A[offset+i] = u[d1off_ + i] - u[d2s0off_ + i*amo_+i];
    }
    offset += amo_;

    // D2_2 symmetric (with zero diagonal)
    for (int i = 0; i < amo_; i++){
        for (int j = 0; j < amo_; j++){
            if ( i == j ) {
                A[offset+i*amo_+i] = u[d2s2off_ + i*amo_ + i];
            }else {
                A[offset+i*amo_+j] = u[d2s2off_ + i*amo_ + j] - u[d2s2off_ + j*amo_+i];
            }
        }
    }

    offset += amo_*amo_;
}

}
