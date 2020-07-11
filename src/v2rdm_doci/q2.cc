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

#include<psi4/libtrans/integraltransform.h>
#include<psi4/libtrans/mospace.h>

#include<psi4/libmints/wavefunction.h>
//#include<psi4/libmints/mints.h>
#include<psi4/libmints/vector.h>
#include<psi4/libmints/matrix.h>
//#include<../bin/fnocc/blas.h>
#include<time.h>

#include"v2rdm_solver.h"

#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() ( (double)clock() / CLOCKS_PER_SEC )
    #define omp_get_max_threads() 1
#endif

//using namespace boost;
using namespace psi;
//using namespace fnocc;

namespace psi{ namespace v2rdm_doci{

// Q2 portion of A.x (with symmetry)
void v2RDMSolver::Q2_constraints_Au(SharedVector A,SharedVector u){

    double * A_p = A->pointer();
    double * u_p = u->pointer();

    // map D2s2 to Q2s2 seniority-2
    //C_DCOPY(amo_*amo_,u_p + d2s2off_,1,A_p + offset,1);      // + D2(ij,ij)
    //C_DAXPY(amo_*amo_,-1.0,u_p + q2s2off_,1,A_p + offset,1); // - Q2(ij,ij)

    // note: AED symmetrized
    int ij = 0;
    for (int i = 0; i < amo_; i++) {
        for (int j = i + 1; j < amo_; j++) {
            //if ( i == j ) continue;
            double dum = 0.0;

            dum += 0.5 * u_p[d2s2off_ + i*amo_ + j]; //  D2(ij,ij)
            dum += 0.5 * u_p[d2s2off_ + j*amo_ + i]; //  D2(ij,ij)
            dum -= u_p[q2s2off_ + ij];         // -Q2(ij,ij)
            dum -= u_p[d1off_ + i];            // -D1(i)
            dum -= u_p[d1off_ + j];            // -D1(j)

            A_p[offset + ij] = dum;

            ij++;
        }
    }
    offset += amo_*(amo_-1)/2;

    // map D2s0 to Q2s0 seniority-0
    C_DCOPY(amo_*amo_,u_p + d2s0off_,1,A_p + offset,1);      // + D2(ii,jj)
    C_DAXPY(amo_*amo_,-1.0,u_p + q2s0off_,1,A_p + offset,1); // - Q2(ii,jj)

    #pragma omp parallel for schedule (static)
    for (int i = 0; i < amo_; i++) {
        double dum = 0.0;
        dum        -=  u_p[d1off_ + i];  // -D1(i) 
        dum        -=  u_p[d1off_ + i];  // -D1(j) 
        A_p[offset + i*amo_ + i] += dum;
    }
    offset += amo_*amo_;

}

// Q2 portion of A^T.y (with symmetry)
void v2RDMSolver::Q2_constraints_ATu(SharedVector A,SharedVector u){

    double * A_p = A->pointer();
    double * u_p = u->pointer();

    // map D2s2 to Q2s2
    //C_DAXPY(amo_*amo_, 1.0,u_p + offset,1,A_p + d2s2off_,1); // + D2(ij,ij)
    //C_DAXPY(amo_*amo_,-1.0,u_p + offset,1,A_p + q2s2off_,1); // - Q2(ii,jj)

    // note: AED symmetrized
    int ij = 0;
    for (int i = 0; i < amo_; i++) {
        for (int j = i + 1; j < amo_; j++) {
            //if ( i == j ) continue;
            A_p[d2s2off_ + i*amo_ + j] += 0.5 * u_p[offset + ij]; // + D2(ij,ij)
            A_p[d2s2off_ + j*amo_ + i] += 0.5 * u_p[offset + ij]; // + D2(ij,ij)
            A_p[q2s2off_ + ij]         -= u_p[offset + ij]; // - Q2(ij,ij)
            A_p[d1off_ + i]            -= u_p[offset + ij]; // -D1(i)
            A_p[d1off_ + j]            -= u_p[offset + ij]; // -D1(j)

            ij++;
        }
    }
    offset += amo_*(amo_-1)/2;

    // map D2s0 to Q2s0
    C_DAXPY(amo_*amo_, 1.0,u_p + offset,1,A_p + d2s0off_,1); // + D2(ii,jj)
    C_DAXPY(amo_*amo_,-1.0,u_p + offset,1,A_p + q2s0off_,1); // - Q2(ii,jj)

    for (int i = 0; i < amo_; i++) {
        double val = u_p[offset + i*amo_ + i];
        A_p[d1off_  + i] -= val; // -D1(i)
        A_p[d1off_  + i] -= val; // -D1(i)
    }
    offset += amo_*amo_;

}

}} // end namespaces
