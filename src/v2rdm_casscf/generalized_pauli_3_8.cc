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
 * MERCHANTABILITY or FITNESS FORPARTICULAR PURPOSE.  See the
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
void v2RDMSolver::Generalized_Pauli_3_8_constraints_ATu(SharedVector A,SharedVector u, int state){

    double* A_p = A->pointer();
    double* u_p = u->pointer();
    double ** orb_p = NatOrbs_[state]->pointer();

    //int * x1aoff;
    //int * x1boff;

    //if ( gpc_[state] == GeneralizedPauli_5_8 ) {
    //    x1aoff = q1aoff;
    //    x1boff = q1boff;
    //}else {
    //    x1aoff = d1aoff;
    //    x1boff = d1boff;
    //}

    int off = gpcoff[state][0];

    //GP_N_8_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-1, 1, 0, 0, 0, 0, 0, 0); //= 0,    ##################
    //GP_N_8_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0, -1, 1, 0, 0, 0, 0, 0); //= 0,    ##              ##
    //GP_N_8_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0, 0, -1, 1, 0, 0, 0, 0); //= 0,    ##  ordering    ##
    //GP_N_8_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0, 0, 0, -1, 1, 0, 0, 0); //= 0,    ##    for       ##
    //GP_N_8_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0, 0, 0, 0, -1, 1, 0, 0); //= 0,    ## lambda[i]'s  ##
    //GP_N_8_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0, 0, 0, 0, 0, -1, 1, 0); //= 0,    ##              ##
    //GP_N_8_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0, 0, 0, 0, 0, 0, -1, 1); //= 0,    ##              ##
    //GP_N_8_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1, 1, 1, 1, 1, 1, 1, -7); //= 3,    ##################

    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,3, -1, -1, -1, -1, -1, -1, 3); //= 1,  ## Extended Pauli inequality: lambda[1]+lambda[8]<=1 ##

    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-1, 1, 1, 1, 1, -1, -1, -1); //= 1, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1, 1, -1, -1, 1, 1, -1, -1); //= 1, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1, 1, -1, 1, -1, -1, 1, -1); //= 1, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1, -1, 1, 1, -1, 1, -1, -1); //= 1, 

    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,2, 1, -2, -1, 0, -1, 0, 1); //= 1, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,2, -1, 0, -1, 0, 1, -2, 1); //= 1, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,0, 0, 1, 2, -2, -1, -1, 1); //= 1, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1, 2, -2, 0, -1, -1, 0, 1); //= 1, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,2, -1, 0, 1, -2, -1, 0, 1); //= 1, 

    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,5, 5, -7, -3, -3, 1, 1, 1); //= 3, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,5, -3, -3, 1, 1, 5, -7, 1); //= 3, 

    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,5, 1, -3, 1, -3, 1, -3, 1); //= 3, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1, 1, 1, 5, -3, -3, -3, 1); //= 3, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1, 5, -3, 1, 1, -3, -3, 1); //= 3, 


    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9, 1, -7, -7, -7, 1, 1, 9); //= 3, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9, -7, -7, 1, 1, 1, -7, 9); //= 3, 

    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7, -1, -1, -1, -1, 7, -9, -1); //= 5, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7, -1, -1, 7, -9, -1, -1, -1); //= 5, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7, 7, -9, -1, -1, -1, -1, -1); //= 5,
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-1, -1, 7, 7, -1, -1, -9, -1); //= 5,  
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-1, 7, -1, 7, -1, -9, -1, -1); //= 5, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-1, 7, -1, -1, 7, -1, -9, -1); //= 5, 

    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-3, 5, 5, 13, -11, -3, -11, 5); //= 7, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,5, 13, -11, 5, -11, -3, -3, 5); //= 7, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,5, -3, 5, 13, -11, -11, -3, 5); //= 7, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,5, 13, -11, -3, 5, -11, -3, 5); //= 7, 

    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19, 11, -21, -13, -5, -5, 3, 11); //= 9, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19, -13, -5, -5, 3, 11, -21, 11); //= 9, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,11, 19, -21, -5, -13, -5, 3, 11); //= 9, 
    GP_N_8_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-5, 3, 11, 19, -21, -13, -5, 11); //= 9]:
}

// portion of A.x corresponding to generalized paulit constraints
void v2RDMSolver::Generalized_Pauli_3_8_constraints_Au(SharedVector A,SharedVector u, int state){

    int saveoff = offset;

    double* A_p = A->pointer();
    double* u_p = u->pointer();
    double ** orb_p = NatOrbs_[state]->pointer();

    //int * x1aoff;
    //int * x1boff;

    //if ( gpc_[state] == GeneralizedPauli_5_8 ) {
    //    x1aoff = q1aoff;
    //    x1boff = q1boff;
    //}else {
    //    x1aoff = d1aoff;
    //    x1boff = d1boff;
    //}

    double * eigvals = (double*)malloc(8*sizeof(double));
    for (int i = 0; i < 8; i++) {
        eigvals[i] = Generalized_Pauli_Au_term(orb_p,u_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],i+1,gpc_rdm_nrm_[state]);
//printf("%5i %5i %20.12lf\n",state,i,eigvals[i]);
    }

    int off = gpcoff[state][0];

    //A_p[offset++] = GP_N_8_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-1, 1, 0, 0, 0, 0, 0, 0); //= 0,    ##################
    //A_p[offset++] = GP_N_8_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0, -1, 1, 0, 0, 0, 0, 0); //= 0,    ##              ##
    //A_p[offset++] = GP_N_8_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0, 0, -1, 1, 0, 0, 0, 0); //= 0,    ##  ordering    ##
    //A_p[offset++] = GP_N_8_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0, 0, 0, -1, 1, 0, 0, 0); //= 0,    ##    for       ##
    //A_p[offset++] = GP_N_8_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0, 0, 0, 0, -1, 1, 0, 0); //= 0,    ## lambda[i]'s  ##
    //A_p[offset++] = GP_N_8_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0, 0, 0, 0, 0, -1, 1, 0); //= 0,    ##              ##
    //A_p[offset++] = GP_N_8_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0, 0, 0, 0, 0, 0, -1, 1); //= 0,    ##              ##
    //A_p[offset++] = GP_N_8_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1, 1, 1, 1, 1, 1, 1, -7); //= 3,    ##################
               
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,3, -1, -1, -1, -1, -1, -1, 3); //= 1,  ## Extended Pauli inequality: lambda[1]+lambda[8]<=1 ##
               
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,-1, 1, 1, 1, 1, -1, -1, -1); //= 1, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,1, 1, -1, -1, 1, 1, -1, -1); //= 1, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,1, 1, -1, 1, -1, -1, 1, -1); //= 1, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,1, -1, 1, 1, -1, 1, -1, -1); //= 1, 
               
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,2, 1, -2, -1, 0, -1, 0, 1); //= 1, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,2, -1, 0, -1, 0, 1, -2, 1); //= 1, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,0, 0, 1, 2, -2, -1, -1, 1); //= 1, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,1, 2, -2, 0, -1, -1, 0, 1); //= 1, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,2, -1, 0, 1, -2, -1, 0, 1); //= 1, 
               
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,5, 5, -7, -3, -3, 1, 1, 1); //= 3, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,5, -3, -3, 1, 1, 5, -7, 1); //= 3, 
               
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,5, 1, -3, 1, -3, 1, -3, 1); //= 3, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,1, 1, 1, 5, -3, -3, -3, 1); //= 3, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,1, 5, -3, 1, 1, -3, -3, 1); //= 3, 
               
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,9, 1, -7, -7, -7, 1, 1, 9); //= 3, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,9, -7, -7, 1, 1, 1, -7, 9); //= 3, 
               
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,7, -1, -1, -1, -1, 7, -9, -1); //= 5, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,7, -1, -1, 7, -9, -1, -1, -1); //= 5, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,7, 7, -9, -1, -1, -1, -1, -1); //= 5,
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,-1, -1, 7, 7, -1, -1, -9, -1); //= 5,  
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,-1, 7, -1, 7, -1, -9, -1, -1); //= 5, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,-1, 7, -1, -1, 7, -1, -9, -1); //= 5, 
               
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,-3, 5, 5, 13, -11, -3, -11, 5); //= 7, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,5, 13, -11, 5, -11, -3, -3, 5); //= 7, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,5, -3, 5, 13, -11, -11, -3, 5); //= 7, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,5, 13, -11, -3, 5, -11, -3, 5); //= 7, 
               
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,19, 11, -21, -13, -5, -5, 3, 11); //= 9, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,19, -13, -5, -5, 3, 11, -21, 11); //= 9, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,11, 19, -21, -5, -13, -5, 3, 11); //= 9, 
    A_p[offset++] = GP_N_8_Au(off,u_p,orb_p,eigvals,-5, 3, 11, 19, -21, -13, -5, 11); //= 9]:

    if ( print_gpc_error_ ) {
        outfile->Printf("\n");        outfile->Printf("    ==> Generalized Pauli Constraint Errors <===\n");
        outfile->Printf("\n");
        for (int i = saveoff; i < saveoff+n_gpc_/n_gpc_states_; i++) {
            outfile->Printf("    %5i %20.12lf %20.12lf %5s\n",i,A_p[i],b->pointer()[i],A_p[i] <= b->pointer()[i] ? "" : "XXX");
        }
    }

    free(eigvals);
}

void v2RDMSolver::GP_N_8_ATu(double rdm_nrm, double dum,int & off, double * A, int *** map_a, int *** map_b, double ** orbs, 
        int d1, int d2, int d3, int d4, int d5, int d6, int d7,int d8) {
    A[off++] = dum;
    Generalized_Pauli_ATu_term(dum * d1 / rdm_nrm,orbs,A,map_a,map_b,1);
    Generalized_Pauli_ATu_term(dum * d2 / rdm_nrm,orbs,A,map_a,map_b,2);
    Generalized_Pauli_ATu_term(dum * d3 / rdm_nrm,orbs,A,map_a,map_b,3);
    Generalized_Pauli_ATu_term(dum * d4 / rdm_nrm,orbs,A,map_a,map_b,4);
    Generalized_Pauli_ATu_term(dum * d5 / rdm_nrm,orbs,A,map_a,map_b,5);
    Generalized_Pauli_ATu_term(dum * d6 / rdm_nrm,orbs,A,map_a,map_b,6);
    Generalized_Pauli_ATu_term(dum * d7 / rdm_nrm,orbs,A,map_a,map_b,7);
    Generalized_Pauli_ATu_term(dum * d8 / rdm_nrm,orbs,A,map_a,map_b,8);

}

double v2RDMSolver::GP_N_8_Au(int & off, double * u, double ** orbs, double * eigvals,
        int d1, int d2, int d3, int d4, int d5, int d6, int d7,int d8) {

    double dum = u[off++]
               + d1  * eigvals[0]
               + d2  * eigvals[1]
               + d3  * eigvals[2]
               + d4  * eigvals[3]
               + d5  * eigvals[4]
               + d6  * eigvals[5]
               + d7  * eigvals[6]
               + d8  * eigvals[7];
    return dum;
}

}
