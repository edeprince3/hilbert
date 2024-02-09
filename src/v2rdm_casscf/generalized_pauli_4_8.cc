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
void v2RDMSolver::Generalized_Pauli_4_8_constraints_ATu(double* A,double* u, int state){

    double ** orb_p = NatOrbs_[state]->pointer();

    //int * x1aoff;
    //int * x1boff;

    //x1aoff = d1aoff;
    //x1boff = d1boff;

    int off = gpcoff[state][0];

    //GP_N_8_ATu(u[offset++],off,A,x1aoff,x1boff,orb_p,-1,1,0,0,0,0,0,0); //=0,   ########################
    //GP_N_8_ATu(u[offset++],off,A,x1aoff,x1boff,orb_p,0,-1,1,0,0,0,0,0); //=0,   ##                    ##
    //GP_N_8_ATu(u[offset++],off,A,x1aoff,x1boff,orb_p,0,0,-1,1,0,0,0,0); //=0,   ##                    ##
    //GP_N_8_ATu(u[offset++],off,A,x1aoff,x1boff,orb_p,0,0,0,-1,1,0,0,0); //=0,   ##    ordering for    ##
    //GP_N_8_ATu(u[offset++],off,A,x1aoff,x1boff,orb_p,0,0,0,0,-1,1,0,0); //=0,   ##    lambda[i]'s     ##
    //GP_N_8_ATu(u[offset++],off,A,x1aoff,x1boff,orb_p,0,0,0,0,0,-1,1,0); //=0,   ##                    ##
    //GP_N_8_ATu(u[offset++],off,A,x1aoff,x1boff,orb_p,0,0,0,0,0,0,-1,1); //=0,   ##                    ##
    //GP_N_8_ATu(u[offset++],off,A,x1aoff,x1boff,orb_p,0,0,0,0,0,0,0,-1); //=0,   ########################

    //GP_N_8_ATu(u[offset++],off,A,x1aoff,x1boff,orb_p,1,0,0,0,0,0,0,0); //=1,    ##  Pauli inequality:lambda[1); //=1 ##

    GP_N_8_ATu(gpc_rdm_sign_a_[state],gpc_rdm_sign_b_[state],gpc_rdm_nrm_[state],u[offset++],off,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,0,0,-1,0,-1,-1,0); //=0,
    GP_N_8_ATu(gpc_rdm_sign_a_[state],gpc_rdm_sign_b_[state],gpc_rdm_nrm_[state],u[offset++],off,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,0,0,-1,-1,0,0,-1); //=0,
    GP_N_8_ATu(gpc_rdm_sign_a_[state],gpc_rdm_sign_b_[state],gpc_rdm_nrm_[state],u[offset++],off,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,0,-1,0,0,-1,0,-1); //=0,
    GP_N_8_ATu(gpc_rdm_sign_a_[state],gpc_rdm_sign_b_[state],gpc_rdm_nrm_[state],u[offset++],off,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,-1,0,0,0,0,-1,-1); //=0,
    GP_N_8_ATu(gpc_rdm_sign_a_[state],gpc_rdm_sign_b_[state],gpc_rdm_nrm_[state],u[offset++],off,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,0,1,0,-1,0,-1,0,-1); //=0,
    GP_N_8_ATu(gpc_rdm_sign_a_[state],gpc_rdm_sign_b_[state],gpc_rdm_nrm_[state],u[offset++],off,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,0,0,1,-1,0,0,-1,-1); //=0,
    GP_N_8_ATu(gpc_rdm_sign_a_[state],gpc_rdm_sign_b_[state],gpc_rdm_nrm_[state],u[offset++],off,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,0,0,0,0,1,-1,-1,-1); //=0,

    GP_N_8_ATu(gpc_rdm_sign_a_[state],gpc_rdm_sign_b_[state],gpc_rdm_nrm_[state],u[offset++],off,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,1,1,-1,0,0,0,0); //=2,
    GP_N_8_ATu(gpc_rdm_sign_a_[state],gpc_rdm_sign_b_[state],gpc_rdm_nrm_[state],u[offset++],off,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,1,0,0,0,0,1,-1); //=2,
    GP_N_8_ATu(gpc_rdm_sign_a_[state],gpc_rdm_sign_b_[state],gpc_rdm_nrm_[state],u[offset++],off,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,1,0,0,1,-1,0,0); //=2,
    GP_N_8_ATu(gpc_rdm_sign_a_[state],gpc_rdm_sign_b_[state],gpc_rdm_nrm_[state],u[offset++],off,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,0,1,0,1,0,-1,0); //=2,
    GP_N_8_ATu(gpc_rdm_sign_a_[state],gpc_rdm_sign_b_[state],gpc_rdm_nrm_[state],u[offset++],off,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,0,1,0,0,1,0,-1); //=2,
    GP_N_8_ATu(gpc_rdm_sign_a_[state],gpc_rdm_sign_b_[state],gpc_rdm_nrm_[state],u[offset++],off,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,0,0,1,1,0,0,-1); //=2,
    GP_N_8_ATu(gpc_rdm_sign_a_[state],gpc_rdm_sign_b_[state],gpc_rdm_nrm_[state],u[offset++],off,A,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,0,1,1,0,1,0,0,-1); //=2

}

// portion of A.x corresponding to generalized paulit constraints
void v2RDMSolver::Generalized_Pauli_4_8_constraints_Au(double* A,double* u,int state){

    double ** orb_p = NatOrbs_[state]->pointer();

    //int * x1aoff;
    //int * x1boff;

    //x1aoff = d1aoff;
    //x1boff = d1boff;

    double * eigvals = (double*)malloc(8*sizeof(double));
    for (int i = 0; i < 8; i++) {
        eigvals[i] = Generalized_Pauli_Au_term(orb_p,u,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],i+1,gpc_rdm_nrm_[state],gpc_rdm_sign_a_[state],gpc_rdm_sign_b_[state]);
    }

    int off = gpcoff[state][0];
    int saveoff = offset;

    //A[offset++] = GP_N_8_Au(off,u,x1aoff,x1boff,orb_p,eigvals,-1,1,0,0,0,0,0,0); //=0,   ########################
    //A[offset++] = GP_N_8_Au(off,u,x1aoff,x1boff,orb_p,eigvals,0,-1,1,0,0,0,0,0); //=0,   ##                    ##
    //A[offset++] = GP_N_8_Au(off,u,x1aoff,x1boff,orb_p,eigvals,0,0,-1,1,0,0,0,0); //=0,   ##                    ##
    //A[offset++] = GP_N_8_Au(off,u,x1aoff,x1boff,orb_p,eigvals,0,0,0,-1,1,0,0,0); //=0,   ##    ordering for    ##
    //A[offset++] = GP_N_8_Au(off,u,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,-1,1,0,0); //=0,   ##    lambda[i]'s     ##
    //A[offset++] = GP_N_8_Au(off,u,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,0,-1,1,0); //=0,   ##                    ##
    //A[offset++] = GP_N_8_Au(off,u,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,0,0,-1,1); //=0,   ##                    ##
    //A[offset++] = GP_N_8_Au(off,u,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,0,0,0,-1); //=0,   ########################
 
    //A[offset++] = GP_N_8_Au(off,u,x1aoff,x1boff,orb_p,eigvals,1,0,0,0,0,0,0,0); //=1,    ##  Pauli inequality:lambda[1); //=1 ##
 
    A[offset++] = GP_N_8_Au(off,u,orb_p,eigvals,1,0,0,-1,0,-1,-1,0); //=0,
    A[offset++] = GP_N_8_Au(off,u,orb_p,eigvals,1,0,0,-1,-1,0,0,-1); //=0,
    A[offset++] = GP_N_8_Au(off,u,orb_p,eigvals,1,0,-1,0,0,-1,0,-1); //=0,
    A[offset++] = GP_N_8_Au(off,u,orb_p,eigvals,1,-1,0,0,0,0,-1,-1); //=0,
    A[offset++] = GP_N_8_Au(off,u,orb_p,eigvals,0,1,0,-1,0,-1,0,-1); //=0,
    A[offset++] = GP_N_8_Au(off,u,orb_p,eigvals,0,0,1,-1,0,0,-1,-1); //=0,
    A[offset++] = GP_N_8_Au(off,u,orb_p,eigvals,0,0,0,0,1,-1,-1,-1); //=0,
 
    A[offset++] = GP_N_8_Au(off,u,orb_p,eigvals,1,1,1,-1,0,0,0,0); //=2,
    A[offset++] = GP_N_8_Au(off,u,orb_p,eigvals,1,1,0,0,0,0,1,-1); //=2,
    A[offset++] = GP_N_8_Au(off,u,orb_p,eigvals,1,1,0,0,1,-1,0,0); //=2,
    A[offset++] = GP_N_8_Au(off,u,orb_p,eigvals,1,0,1,0,1,0,-1,0); //=2,
    A[offset++] = GP_N_8_Au(off,u,orb_p,eigvals,1,0,1,0,0,1,0,-1); //=2,
    A[offset++] = GP_N_8_Au(off,u,orb_p,eigvals,1,0,0,1,1,0,0,-1); //=2,
    A[offset++] = GP_N_8_Au(off,u,orb_p,eigvals,0,1,1,0,1,0,0,-1); //=2

    if ( print_gpc_error_ ) {
        outfile->Printf("\n");        outfile->Printf("    ==> Generalized Pauli Constraint Errors <===\n");
        outfile->Printf("\n");
        for (int i = saveoff; i < saveoff+n_gpc_[state]; i++) {
            outfile->Printf("    %5i %20.12lf %20.12lf %5s\n",i,A[i],b->pointer()[i],A[i] <= b->pointer()[i] ? "" : "XXX");
        }
    }

    free(eigvals);
}

}
