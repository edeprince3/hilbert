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

#include <psi4/libpsi4util/PsiOutStream.h>

#include"v2rdm_solver.h"

using namespace psi;

namespace hilbert{


// portion of A^T.y corresponding to generalized paulit constraints
void v2RDMSolver::Generalized_Pauli_3_6_constraints_ATu(double* A,double* u, int state){

    double ** orb_p = NatOrbs_[state]->pointer();

    // l1 + l6 = 1
    double dum1 = u[ offset++ ];

    // l2 + l5 = 1
    double dum2 = u[ offset++ ];

    // l3 + l4 = 1
    double dum3 = u[ offset++ ];

    // l4 - l5 - l6 <= 0
    double dum4 = u[ offset++ ];


    Generalized_Pauli_ATu_term(gpc_rdm_sign_a_[state], gpc_rdm_sign_b_[state],dum1 / gpc_rdm_nrm_[state],orb_p,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],1);
    Generalized_Pauli_ATu_term(gpc_rdm_sign_a_[state], gpc_rdm_sign_b_[state],dum2 / gpc_rdm_nrm_[state],orb_p,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],2);
    Generalized_Pauli_ATu_term(gpc_rdm_sign_a_[state], gpc_rdm_sign_b_[state],dum3 / gpc_rdm_nrm_[state],orb_p,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],3);

    Generalized_Pauli_ATu_term(gpc_rdm_sign_a_[state], gpc_rdm_sign_b_[state],(dum3+dum4) / gpc_rdm_nrm_[state],orb_p,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],4);
    Generalized_Pauli_ATu_term(gpc_rdm_sign_a_[state], gpc_rdm_sign_b_[state],(dum2-dum4) / gpc_rdm_nrm_[state],orb_p,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],5);
    Generalized_Pauli_ATu_term(gpc_rdm_sign_a_[state], gpc_rdm_sign_b_[state],(dum1-dum4) / gpc_rdm_nrm_[state],orb_p,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],6);

    A[ gpcoff[state][0] ]  = 0.0;
    A[ gpcoff[state][1] ]  = 0.0;
    A[ gpcoff[state][2] ]  = 0.0;
    A[ gpcoff[state][3] ] += dum4;

}

// portion of A.x corresponding to generalized paulit constraints
void v2RDMSolver::Generalized_Pauli_3_6_constraints_Au(double* A,double* u, int state){

    int saveoff = offset;

    double ** orb_p = NatOrbs_[state]->pointer();

    double l1 = Generalized_Pauli_Au_term(orb_p,u,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],1,gpc_rdm_nrm_[state],gpc_rdm_sign_a_[state], gpc_rdm_sign_b_[state]);
    double l2 = Generalized_Pauli_Au_term(orb_p,u,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],2,gpc_rdm_nrm_[state],gpc_rdm_sign_a_[state], gpc_rdm_sign_b_[state]);
    double l3 = Generalized_Pauli_Au_term(orb_p,u,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],3,gpc_rdm_nrm_[state],gpc_rdm_sign_a_[state], gpc_rdm_sign_b_[state]);
    double l4 = Generalized_Pauli_Au_term(orb_p,u,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],4,gpc_rdm_nrm_[state],gpc_rdm_sign_a_[state], gpc_rdm_sign_b_[state]);
    double l5 = Generalized_Pauli_Au_term(orb_p,u,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],5,gpc_rdm_nrm_[state],gpc_rdm_sign_a_[state], gpc_rdm_sign_b_[state]);
    double l6 = Generalized_Pauli_Au_term(orb_p,u,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],6,gpc_rdm_nrm_[state],gpc_rdm_sign_a_[state], gpc_rdm_sign_b_[state]);

    // l1 + l6 = 1
    A[offset++] = l1 + l6;

    // l2 + l5 = 1
    A[offset++] = l2 + l5;

    // l3 + l4 = 1
    A[offset++] = l3 + l4;

    // l4 - l5 - l6 <= 0
    A[offset++] = u[ gpcoff[state][3] ] + l4 - l5 - l6;

    if ( print_gpc_error_ ) {
        outfile->Printf("\n");        outfile->Printf("    ==> Generalized Pauli Constraint Errors <===\n");
        outfile->Printf("\n");
        for (int i = saveoff; i < saveoff+n_gpc_[state]; i++) {
            outfile->Printf("    %5i %20.12lf %20.12lf %5s\n",i,A[i],b->pointer()[i],A[i] <= b->pointer()[i] ? "" : "XXX");
        }
    }
}

}
