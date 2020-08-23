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

#include <psi4/libpsi4util/PsiOutStream.h>

#include"v2rdm_solver.h"

using namespace psi;

namespace hilbert{


// portion of A^T.y corresponding to generalized paulit constraints
void v2RDMSolver::Generalized_Pauli_3_6_constraints_ATu(SharedVector A,SharedVector u, int state){

    double* A_p = A->pointer();
    double* u_p = u->pointer();
    double ** orb_p = NatOrbs_[state]->pointer();

    // l1 + l6 = 1
    double dum1 = u_p[ offset++ ];

    // l2 + l5 = 1
    double dum2 = u_p[ offset++ ];

    // l3 + l4 = 1
    double dum3 = u_p[ offset++ ];

    // l4 - l5 - l6 <= 0
    double dum4 = u_p[ offset++ ];


    Generalized_Pauli_ATu_term(dum1,orb_p,A_p,d1aoff,d1boff,1);
    Generalized_Pauli_ATu_term(dum2,orb_p,A_p,d1aoff,d1boff,2);
    Generalized_Pauli_ATu_term(dum3,orb_p,A_p,d1aoff,d1boff,3);

    Generalized_Pauli_ATu_term(dum3+dum4,orb_p,A_p,d1aoff,d1boff,4);
    Generalized_Pauli_ATu_term(dum2-dum4,orb_p,A_p,d1aoff,d1boff,5);
    Generalized_Pauli_ATu_term(dum1-dum4,orb_p,A_p,d1aoff,d1boff,6);

    A_p[ gpcoff[state][0] ]  = 0.0;
    A_p[ gpcoff[state][1] ]  = 0.0;
    A_p[ gpcoff[state][2] ]  = 0.0;
    A_p[ gpcoff[state][3] ] += dum4;

}

// portion of A.x corresponding to generalized paulit constraints
void v2RDMSolver::Generalized_Pauli_3_6_constraints_Au(SharedVector A,SharedVector u, int state){

    int saveoff = offset;

    double* A_p = A->pointer();
    double* u_p = u->pointer();
    double ** orb_p = NatOrbs_[state]->pointer();

    double l1 = Generalized_Pauli_Au_term(orb_p,u_p,d1aoff,d1boff,1);
    double l2 = Generalized_Pauli_Au_term(orb_p,u_p,d1aoff,d1boff,2);
    double l3 = Generalized_Pauli_Au_term(orb_p,u_p,d1aoff,d1boff,3);
    double l4 = Generalized_Pauli_Au_term(orb_p,u_p,d1aoff,d1boff,4);
    double l5 = Generalized_Pauli_Au_term(orb_p,u_p,d1aoff,d1boff,5);
    double l6 = Generalized_Pauli_Au_term(orb_p,u_p,d1aoff,d1boff,6);

    // l1 + l6 = 1
    A_p[offset++] = l1 + l6;

    // l2 + l5 = 1
    A_p[offset++] = l2 + l5;

    // l3 + l4 = 1
    A_p[offset++] = l3 + l4;

    // l4 - l5 - l6 <= 0
    A_p[offset++] = u_p[ gpcoff[state][3] ] + l4 - l5 - l6;

    if ( print_gpc_error_ ) {
        outfile->Printf("\n");        outfile->Printf("    ==> Generalized Pauli Constraint Errors <===\n");
        outfile->Printf("\n");
        for (int i = saveoff; i < saveoff+n_gpc_/n_gpc_states_; i++) {
            outfile->Printf("    %5i %20.12lf %20.12lf %5s\n",i,A_p[i],b->pointer()[i],A_p[i] <= b->pointer()[i] ? "" : "XXX");
        }
    }
}

}
