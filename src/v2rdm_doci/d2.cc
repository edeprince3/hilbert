/*
 *@BEGIN LICENSE
 *
 * v2RDM-DOCI, a plugin to:
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

#include"v2rdm_solver.h"

#include <misc/omp.h>

using namespace psi;

namespace psi{ namespace v2rdm_doci{

// D2 portion of A^T.y ( and D1 / Q1 ) 
void v2RDMSolver::D2_constraints_ATu(SharedVector A,SharedVector u){
    double* A_p = A->pointer();
    double* u_p = u->pointer();

    // Traces
    // Tr(D2s2) seniority-2
    for (int i = 0; i < amo_; i++){
        for (int j = 0; j < amo_; j++){
            if ( i==j ) continue;
            A_p[d2s2off_ + i*amo_+j] += u_p[offset];
        }
    }
    offset++;

    // Tr(D2s0) seniority-0
    for (int i = 0; i < amo_; i++){
        A_p[d2s0off_ + i*amo_+i] += u_p[offset];
    }
    offset++;

    // Tr(D1a)
    for (int i = 0; i < amo_; i++){
        A_p[d1off_ + i] += u_p[offset];
    }
    offset++;

    int na = nalpha_;

    // contraction: D2s2 -> D1 a seniority-2
    for (int i = 0; i < amo_; i++){
        A_p[d1off_ + i] += (na-1.0)*u_p[offset + i];
        for (int j = 0; j < amo_; j++){
            if (i==j) continue; 
            A_p[d2s2off_ + i*amo_ + j] -= u_p[offset + i];
        }
    }
    offset += amo_;

    // contraction: D2s2 -> D1 b seniority-2
    for (int i = 0; i < amo_; i++){
        A_p[d1off_ + i] += (na-1.0)*u_p[offset + i];
        for (int j = 0; j < amo_; j++){
            if (i==j) continue; 
            A_p[d2s2off_ + j*amo_ + i] -= u_p[offset + i];
        }
    }
    offset += amo_;

    //contract D2s0 -> D1 a seniority-0
    for (int i = 0; i < amo_; i++){
        A_p[d1off_ + i]          += u_p[offset + i];
        A_p[d2s0off_ + i*amo_+i] -= u_p[offset + i];
    }
    offset += amo_;

    // D2_2 symmetric (with zero diagonal)
    for (int i = 0; i < amo_; i++){
        for (int j = 0; j < amo_; j++){
            if ( i == j ) {
                A_p[d2s2off_ + i*amo_ + i] += u_p[offset+i*amo_+i];
            }else {
                A_p[d2s2off_ + i*amo_ + j] += u_p[offset+i*amo_+j];
                A_p[d2s2off_ + j*amo_ + i] -= u_p[offset+i*amo_+j];
            }
        }
    }
    offset += amo_*amo_;

}

// D2 portion of A.x (and D1/Q1) 
void v2RDMSolver::D2_constraints_Au(SharedVector A,SharedVector u){

    double* A_p = A->pointer();
    double* u_p = u->pointer();

    // Traces
    // Tr(D2s2) seniority-2
    double sums2 =0.0;
    for (int i = 0; i < amo_; i++){
        for (int j = 0; j < amo_; j++){
            if ( i==j ) continue;
            sums2 += u_p[d2s2off_ + i*amo_ + j];
        }
    }
    A_p[offset] = sums2;
    offset++;

    // Tr(D2s0) seniority-0
    double sums0 =0.0;
    for (int i = 0; i < amo_; i++){
        sums0 += u_p[d2s0off_ + i*amo_ + i];
    }
    A_p[offset] = sums0;
    offset++;

    // Tr(D1a)
    double sums1 =0.0;
    for (int i = 0; i < amo_; i++){
        sums1 += u_p[d1off_ + i];
    }
    A_p[offset] = sums1;
    offset++;

    int na = nalpha_;

    // contraction: D2s2 -> D1 a seniority-2
    for (int i = 0; i < amo_; i++){
        double sum = (na-1.0)*u_p[d1off_ + i];
        for (int j = 0; j < amo_; j++){
            if (i==j) continue; 
            sum -= u_p[d2s2off_ + i*amo_ + j];
        }
        A_p[offset + i] = sum;
    }
    offset += amo_;

    // contraction: D2s2 -> D1 b seniority-2
    for (int i = 0; i < amo_; i++){
        double sum = (na-1.0)*u_p[d1off_ + i];
        for (int j = 0; j < amo_; j++){
            if (i==j) continue; 
            sum -= u_p[d2s2off_ + j*amo_ + i];
        }
        A_p[offset + i] = sum;
    }
    offset += amo_;

    //contract D2s0 -> D1 a seniority-0
    for (int i = 0; i < amo_; i++){
        A_p[offset+i] = u_p[d1off_ + i] - u_p[d2s0off_ + i*amo_+i];
    }
    offset += amo_;

    // D2_2 symmetric (with zero diagonal)
    for (int i = 0; i < amo_; i++){
        for (int j = 0; j < amo_; j++){
            if ( i == j ) {
                A_p[offset+i*amo_+i] = u_p[d2s2off_ + i*amo_ + i];
            }else {
                A_p[offset+i*amo_+j] = u_p[d2s2off_ + i*amo_ + j] - u_p[d2s2off_ + j*amo_+i];
            }
        }
    }

    offset += amo_*amo_;
}

}}
