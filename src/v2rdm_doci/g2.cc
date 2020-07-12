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
#include<psi4/libmints/vector.h>
#include<psi4/libmints/matrix.h>


#include"v2rdm_solver.h"

#include <misc/omp.h>

using namespace psi;

namespace psi{ namespace v2rdm_doci{

// G2 portion of A.x (with symmetry)
void v2RDMSolver::G2_constraints_Au(SharedVector A,SharedVector u){

    double* A_p = A->pointer();
    double* u_p = u->pointer();

    int aux = 0;
    // Eq (35) ... note AED symmetried
    for (int i = 0; i < amo_; i++) {
        for (int j = i + 1; j < amo_; j++) {
            //if ( i == j ) continue;

            double dum;

            dum    = -u_p[g2s2off_ + aux + 0];     // - G2s2(ij,ij)
            dum   +=  u_p[d1off_ + i];             // + D1(i)
            dum   -=  0.5 * u_p[d2s2off_ + i*amo_+ j];   // - D2s2(ij,ij)
            dum   -=  0.5 * u_p[d2s2off_ + j*amo_+ i];   // - D2s2(ij,ij)

            A_p[offset + 0] = dum;

            dum    = -u_p[g2s2off_ + aux + 1];     // - G2s2(ij,ij)
            dum   +=  0.5 * u_p[d2s0off_ + i*amo_+ j];   // + D2s0(ij,ij)
            dum   +=  0.5 * u_p[d2s0off_ + j*amo_+ i];   // + D2s0(ij,ij)

            A_p[offset + 1] = dum;

            dum    = -u_p[g2s2off_ + aux + 2];     // - G2s2(ij,ij)
            dum   +=  0.5 * u_p[d2s0off_ + i*amo_+ j];   // + D2s0(ij,ij)
            dum   +=  0.5 * u_p[d2s0off_ + j*amo_+ i];   // + D2s0(ij,ij)

            A_p[offset + 2] = dum;

            dum    = -u_p[g2s2off_ + aux + 3];     // - G2s2(ij,ij)
            dum   +=  u_p[d1off_ + j];             // + D1(i)
            dum   -=  0.5 * u_p[d2s2off_ + i*amo_+ j];   // - D2s2(ij,ij)
            dum   -=  0.5 * u_p[d2s2off_ + j*amo_+ i];   // - D2s2(ij,ij)

            A_p[offset + 3] = dum;

            offset += 2*2;
            aux += 2*2;
        }
    }

    // G2s0 (note AED symmetrized):
    for (int i = 0; i < amo_; i++) {
        for (int j = 0; j < amo_; j++) {
            double dum = -u_p[g2s0off_ + i*amo_+j];
            if ( i != j ) {
                dum += 0.5 * u_p[d2s2off_ + i*amo_+j];
                dum += 0.5 * u_p[d2s2off_ + j*amo_+i];
            }else {
                dum += u_p[d1off_ + i];
            }

            A_p[offset + i*amo_+j] = dum;

        }
    }
    offset += amo_*amo_;

}

// G2 portion of A^T.y (with symmetry)
void v2RDMSolver::G2_constraints_ATu(SharedVector A,SharedVector u){

    double* A_p = A->pointer();
    double* u_p = u->pointer();

    // Eq (35) ... note AED symmetrized
    int aux = 0;
    for (int i = 0; i < amo_; i++) {
        for (int j = i + 1; j < amo_; j++) {
            //if ( i == j ) continue;

            double dum = u_p[offset + 0];

            A_p[g2s2off_ + aux + 0]   -= dum;
            A_p[d1off_ + i]           += dum;
            A_p[d2s2off_ + i*amo_+j]  -= 0.5 * dum;
            A_p[d2s2off_ + j*amo_+i]  -= 0.5 * dum;

            dum = u_p[offset + 1];

            A_p[g2s2off_ + aux + 1]   -= dum;
            A_p[d2s0off_ + i*amo_+ j] += 0.5 * dum;
            A_p[d2s0off_ + j*amo_+ i] += 0.5 * dum;

            dum = u_p[offset + 2];

            A_p[g2s2off_ + aux + 2]   -= dum;
            A_p[d2s0off_ + i*amo_+ j] += 0.5 * dum;
            A_p[d2s0off_ + j*amo_+ i] += 0.5 * dum;

            dum = u_p[offset + 3];
            A_p[g2s2off_ + aux + 3]   -= dum;
            A_p[d1off_ + j]           += dum;
            A_p[d2s2off_ + j*amo_+i]  -= 0.5 * dum;
            A_p[d2s2off_ + i*amo_+j]  -= 0.5 * dum;

            aux    += 2*2;
            offset += 2*2;
        }
    }

    // G2s0 ... note AED symmetrized:
    for (int i = 0; i < amo_; i++) {
        for (int j = 0; j < amo_; j++) {
            double dum = u_p[offset + i*amo_+j];
            A_p[g2s0off_ + i*amo_+j] -= dum;
            if ( i != j ) {
                A_p[d2s2off_ + i*amo_+j] += 0.5 * dum;
                A_p[d2s2off_ + j*amo_+i] += 0.5 * dum;
            }else {
                A_p[d1off_ + i] += dum;
            }

        }
    }
    offset += amo_*amo_;

}

}} // end namespaces
