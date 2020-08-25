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
void v2RDMSolver::Generalized_Pauli_3_10_constraints_ATu(SharedVector A,SharedVector u, int state){

    double* A_p = A->pointer();
    double* u_p = u->pointer();
    double ** orb_p = NatOrbs_[state]->pointer();

    //int * x1aoff;
    //int * x1boff;

    //if ( gpc_[state] == GeneralizedPauli_7_10 ) {
    //    x1aoff = q1aoff;
    //    x1boff = q1boff;
    //}else {
    //    x1aoff = d1aoff;
    //    x1boff = d1boff;
    //}

    int off = gpcoff[state][0];

    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-1,1,0,0,0,0,0,0,0,0);//=0      #################
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,0,0,0,0,0,-1,1,0);//=0     ##             ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,0,0,-1,1,0,0,0,0);//=0     ##             ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,0,-1,1,0,0,0,0,0);//=0     ##  ordering   ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,0,0,0,-1,1,0,0,0);//=0     ##     for     ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,-1,1,0,0,0,0,0,0);//=0     ## lambda[i]'s ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,-1,1,0,0,0,0,0,0,0);//=0     ##             ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,0,0,0,0,-1,1,0,0);//=0     ##             ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,0,0,0,0,0,0,-1,1);//=0     ##             ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,1,1,1,1,1,1,1,1,-9);//=3     #################

    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,1,-1,1,-1,-1,0,0,0,0);//=1

    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,4,-1,-1,-1,-1,-1,-1,-1,-1,4);//=2  ## Extended Pauli inequality:lambda[1]+lambda[10);//=1 ##

    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,6,-4,-4,-4,1,1,1,1,-4,6);//=3
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,6,1,-4,-4,-4,-4,1,1,1,6);//=3

    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,5,3,-5,-3,-1,-1,-1,-1,1,3);//=3
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,5,-3,-1,-1,-1,-1,1,3,-5,3);//=3

    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,3,-2,3,3,-2,3,-2,-2,-2,-2);//=4
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,3,3,-2,3,-2,-2,3,-2,-2,-2);//=4
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,3,3,-2,-2,3,3,-2,-2,-2,-2);//=4
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-2,3,3,3,3,-2,-2,-2,-2,-2);//=4

    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,-3,2,7,-8,2,-3,-3,-3,2);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,-3,2,2,-3,7,-8,-3,-3,2);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,2,-3,7,7,-3,2,-8,-3,-3,2);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,2,-3,-3,-3,2,2,7,-8,-3);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,-3,-3,2,2,7,-8,2,-3,-3);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,-3,-3,2,2,2,-3,7,-8,-3);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,2,-3,-3,2,7,-8,-3,-3,2);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,2,7,-3,-3,7,2,-8,-3,-3,2);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-3,7,2,2,7,-3,-8,-3,-3,2);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-3,2,7,7,2,-3,-8,-3,-3,2);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-3,7,2,7,2,-8,-3,-3,-3,2);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,2,-3,7,-8,-3,2,-3,-3,2);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,2,7,-3,7,-3,-8,2,-3,-3,2);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,2,-3,7,-8,-3,-3,2,2,-3);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,2,7,-3,7,-3,-8,-3,2,2,-3);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,2,7,-3,7,-8,-3,2,-3,2,-3);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,7,-8,2,-3,-3,2,-3,-3,2);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,7,-8,2,-3,-3,-3,2,2,-3);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-3,7,2,2,7,-8,-3,-3,2,-3);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-3,2,7,7,2,-8,-3,-3,2,-3);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,7,-8,-3,2,2,-3,-3,-3,2);//=6
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,7,-8,-3,-3,2,2,2,-3,-3);//=6

    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,4,9,-6,-1,4,-6,-1,4,-6,-1);//=7
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,4,-1,4,9,-6,-6,-1,4,-6,-1);//=7
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-1,4,4,9,-6,-1,-6,4,-6,-1);//=7
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,4,9,-6,4,-6,-1,-1,4,-6,-1);//=7

    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-1,9,-1,9,-1,-11,-1,-1,-1,-1);//=7
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-1,9,-1,-1,9,-1,-11,-1,-1,-1);//=7
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,-1,-1,9,-11,-1,-1,-1,-1,-1);//=7
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,-1,-1,-1,-1,9,-11,-1,-1,-1);//=7
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-1,-1,9,9,-1,-1,-11,-1,-1,-1);//=7
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,9,-11,-1,-1,-1,-1,-1,-1,-1);//=7
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,-1,-1,-1,-1,-1,-1,9,-11,-1);//=7

    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,13,-7,-7,-7,3,3,3,13,-17,3);//=9
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,13,13,-17,-7,-7,-7,3,3,3,3);//=9

    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,13,3,-7,3,-7,3,-7,3,-7,3);//=9
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,3,13,-7,3,3,-7,-7,3,-7,3);//=9
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,3,3,3,13,-7,-7,-7,3,-7,3);//=9

    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,17,-13,-3,7,-13,-3,-3,-3,7);//=11
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,-3,7,17,-13,-13,-3,-3,-3,7);//=11
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-3,7,7,17,-13,-3,-13,-3,-3,7);//=11

    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,14,-11,-1,-6,-1,4,9,-11,-6);//=12
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-1,4,9,14,-11,-6,-1,9,-11,-6);//=12
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,14,-6,4,9,-11,9,-11,-6,-1,-1);//=12
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,-6,9,14,-11,4,-11,-6,-1,-1);//=12
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,14,9,-11,-6,-1,-1,4,9,-11,-6);//=12
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,14,-6,-1,-1,4,9,-11,9,-11,-6);//=12
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,14,9,-11,-6,4,9,-11,-6,-1,-1);//=12
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,14,-11,-6,9,4,-11,-6,-1,-1);//=12
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-6,14,4,9,9,-11,-11,-6,-1,-1);//=12
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-6,9,9,14,4,-11,-11,-6,-1,-1);//=12
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,14,-11,9,-11,-6,-1,4,-1,-6);//=12
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,14,-11,9,-11,-6,4,-1,-6,-1);//=12
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,14,9,-11,9,-11,-6,4,-6,-1,-1);//=12
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,9,-6,14,-11,-11,4,-6,-1,-1);//=12
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,14,-11,9,-6,-11,4,-6,-1,-1);//=12
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,14,9,-11,9,-11,-6,-1,-1,4,-6);//=12
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,9,-6,14,-11,-11,-1,-1,4,-6);//=12
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,14,-11,9,-6,-11,-1,-1,4,-6);//=12

    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,14,-21,-11,-6,-1,-6,-1,4,9);//=12	
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,-11,-6,-1,-6,-1,4,14,-21,9);//=12	

    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,-9,1,-9,1,11,-19,-9,1,11);//=13
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,-9,1,11,-19,-9,1,-9,1,11);//=13
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,11,-19,-9,1,-9,1,-9,1,11);//=13
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,1,11,21,-19,-9,-9,-9,1,11);//=13
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,-9,1,-9,1,-9,1,11,-19,11);//=13

    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,13,23,-17,3,-7,-7,3,13,-17,-7);//=19
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,3,3,13,23,-17,-7,-7,13,-17,-7);//=19
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,23,-7,3,13,-17,13,-17,-7,3,-7);//=19
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,23,-7,3,13,-17,-7,3,13,-17,-7);//=19
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,13,-7,13,23,-17,3,-17,-7,3,-7);//=19
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,23,13,-17,-7,3,-7,3,13,-17,-7);//=19
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,23,-7,3,-7,3,13,-17,13,-17,-7);//=19
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,23,13,-17,-7,3,13,-17,-7,3,-7);//=19
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,13,23,-17,-7,13,3,-17,-7,3,-7);//=19
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-7,23,3,13,13,-17,-17,-7,3,-7);//=19
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-7,13,13,23,3,-17,-17,-7,3,-7);//=19
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,13,23,-17,13,-17,-7,3,3,-7,-7);//=19
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,23,13,-17,13,-17,-7,3,-7,3,-7);//=19
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,13,13,-7,23,-17,-17,3,-7,3,-7);//=19
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,13,23,-17,13,-7,-17,3,-7,3,-7);//=19

    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,37,-23,-13,-3,-13,-3,7,17,-33,27);//=21	
    GP_N_10_ATu(gpc_rdm_nrm_[state],u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,37,17,-33,-23,-13,-3,-13,-3,7,27);//=21	
}

// portion of A.x corresponding to generalized paulit constraints
void v2RDMSolver::Generalized_Pauli_3_10_constraints_Au(SharedVector A,SharedVector u, int state){

    int saveoff = offset;

    double* A_p = A->pointer();
    double* u_p = u->pointer();
    double ** orb_p = NatOrbs_[state]->pointer();

    //int * x1aoff;
    //int * x1boff;
    
    //if ( gpc_[state] == GeneralizedPauli_7_10 ) {
    //    x1aoff = q1aoff;
    //    x1boff = q1boff;
    //}else {
    //    x1aoff = d1aoff;
    //    x1boff = d1boff;
    //}

    double * eigvals = (double*)malloc(10*sizeof(double));
    for (int i = 0; i < 10; i++) {
        //eigvals[i] = Generalized_Pauli_Au_term(orb_p,u_p,x1aoff,x1boff,i+1);
        eigvals[i] = Generalized_Pauli_Au_term(orb_p,u_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],i+1, gpc_rdm_nrm_[state]);
    }

    int off = gpcoff[state][0];

    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-1,1,0,0,0,0,0,0,0,0);//=0      #################
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,0,0,0,-1,1,0);//=0     ##             ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,-1,1,0,0,0,0);//=0     ##             ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,0,-1,1,0,0,0,0,0);//=0     ##  ordering   ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,0,-1,1,0,0,0);//=0     ##     for     ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,-1,1,0,0,0,0,0,0);//=0     ## lambda[i]'s ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,-1,1,0,0,0,0,0,0,0);//=0     ##             ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,0,0,-1,1,0,0);//=0     ##             ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,0,0,0,0,-1,1);//=0     ##             ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,1,1,1,1,1,1,1,1,-9);//=3     #################

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,1,1,-1,1,-1,-1,0,0,0,0);//=1

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,4,-1,-1,-1,-1,-1,-1,-1,-1,4);//=2  ## Extended Pauli inequality:lambda[1]+lambda[10);//=1 ##

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,6,-4,-4,-4,1,1,1,1,-4,6);//=3
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,6,1,-4,-4,-4,-4,1,1,1,6);//=3

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,5,3,-5,-3,-1,-1,-1,-1,1,3);//=3
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,5,-3,-1,-1,-1,-1,1,3,-5,3);//=3

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,3,-2,3,3,-2,3,-2,-2,-2,-2);//=4
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,3,3,-2,3,-2,-2,3,-2,-2,-2);//=4
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,3,3,-2,-2,3,3,-2,-2,-2,-2);//=4
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,-2,3,3,3,3,-2,-2,-2,-2,-2);//=4

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,-3,2,7,-8,2,-3,-3,-3,2);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,-3,2,2,-3,7,-8,-3,-3,2);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,2,-3,7,7,-3,2,-8,-3,-3,2);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,2,-3,-3,-3,2,2,7,-8,-3);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,-3,-3,2,2,7,-8,2,-3,-3);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,-3,-3,2,2,2,-3,7,-8,-3);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,2,-3,-3,2,7,-8,-3,-3,2);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,2,7,-3,-3,7,2,-8,-3,-3,2);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,-3,7,2,2,7,-3,-8,-3,-3,2);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,-3,2,7,7,2,-3,-8,-3,-3,2);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,-3,7,2,7,2,-8,-3,-3,-3,2);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,2,-3,7,-8,-3,2,-3,-3,2);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,2,7,-3,7,-3,-8,2,-3,-3,2);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,2,-3,7,-8,-3,-3,2,2,-3);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,2,7,-3,7,-3,-8,-3,2,2,-3);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,2,7,-3,7,-8,-3,2,-3,2,-3);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,7,-8,2,-3,-3,2,-3,-3,2);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,7,-8,2,-3,-3,-3,2,2,-3);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,-3,7,2,2,7,-8,-3,-3,2,-3);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,-3,2,7,7,2,-8,-3,-3,2,-3);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,7,-8,-3,2,2,-3,-3,-3,2);//=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,7,-8,-3,-3,2,2,2,-3,-3);//=6

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,4,9,-6,-1,4,-6,-1,4,-6,-1);//=7
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,4,-1,4,9,-6,-6,-1,4,-6,-1);//=7
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,-1,4,4,9,-6,-1,-6,4,-6,-1);//=7
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,4,9,-6,4,-6,-1,-1,4,-6,-1);//=7

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,-1,9,-1,9,-1,-11,-1,-1,-1,-1);//=7
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,-1,9,-1,-1,9,-1,-11,-1,-1,-1);//=7
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,-1,-1,9,-11,-1,-1,-1,-1,-1);//=7
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,-1,-1,-1,-1,9,-11,-1,-1,-1);//=7
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,-1,-1,9,9,-1,-1,-11,-1,-1,-1);//=7
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,9,-11,-1,-1,-1,-1,-1,-1,-1);//=7
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,-1,-1,-1,-1,-1,-1,9,-11,-1);//=7

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,13,-7,-7,-7,3,3,3,13,-17,3);//=9
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,13,13,-17,-7,-7,-7,3,3,3,3);//=9

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,13,3,-7,3,-7,3,-7,3,-7,3);//=9
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,3,13,-7,3,3,-7,-7,3,-7,3);//=9
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,3,3,3,13,-7,-7,-7,3,-7,3);//=9

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,17,-13,-3,7,-13,-3,-3,-3,7);//=11
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,-3,7,17,-13,-13,-3,-3,-3,7);//=11
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,-3,7,7,17,-13,-3,-13,-3,-3,7);//=11

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,14,-11,-1,-6,-1,4,9,-11,-6);//=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,-1,4,9,14,-11,-6,-1,9,-11,-6);//=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,14,-6,4,9,-11,9,-11,-6,-1,-1);//=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,-6,9,14,-11,4,-11,-6,-1,-1);//=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,14,9,-11,-6,-1,-1,4,9,-11,-6);//=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,14,-6,-1,-1,4,9,-11,9,-11,-6);//=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,14,9,-11,-6,4,9,-11,-6,-1,-1);//=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,14,-11,-6,9,4,-11,-6,-1,-1);//=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,-6,14,4,9,9,-11,-11,-6,-1,-1);//=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,-6,9,9,14,4,-11,-11,-6,-1,-1);//=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,14,-11,9,-11,-6,-1,4,-1,-6);//=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,14,-11,9,-11,-6,4,-1,-6,-1);//=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,14,9,-11,9,-11,-6,4,-6,-1,-1);//=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,9,-6,14,-11,-11,4,-6,-1,-1);//=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,14,-11,9,-6,-11,4,-6,-1,-1);//=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,14,9,-11,9,-11,-6,-1,-1,4,-6);//=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,9,-6,14,-11,-11,-1,-1,4,-6);//=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,14,-11,9,-6,-11,-1,-1,4,-6);//=12

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,14,-21,-11,-6,-1,-6,-1,4,9);//=12	
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,-11,-6,-1,-6,-1,4,14,-21,9);//=12	

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,-9,1,-9,1,11,-19,-9,1,11);//=13
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,-9,1,11,-19,-9,1,-9,1,11);//=13
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,11,-19,-9,1,-9,1,-9,1,11);//=13
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,1,1,11,21,-19,-9,-9,-9,1,11);//=13
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,-9,1,-9,1,-9,1,11,-19,11);//=13

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,13,23,-17,3,-7,-7,3,13,-17,-7);//=19
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,3,3,13,23,-17,-7,-7,13,-17,-7);//=19
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,23,-7,3,13,-17,13,-17,-7,3,-7);//=19
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,23,-7,3,13,-17,-7,3,13,-17,-7);//=19
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,13,-7,13,23,-17,3,-17,-7,3,-7);//=19
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,23,13,-17,-7,3,-7,3,13,-17,-7);//=19
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,23,-7,3,-7,3,13,-17,13,-17,-7);//=19
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,23,13,-17,-7,3,13,-17,-7,3,-7);//=19
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,13,23,-17,-7,13,3,-17,-7,3,-7);//=19
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,-7,23,3,13,13,-17,-17,-7,3,-7);//=19
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,-7,13,13,23,3,-17,-17,-7,3,-7);//=19
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,13,23,-17,13,-17,-7,3,3,-7,-7);//=19
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,23,13,-17,13,-17,-7,3,-7,3,-7);//=19
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,13,13,-7,23,-17,-17,3,-7,3,-7);//=19
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,13,23,-17,13,-7,-17,3,-7,3,-7);//=19

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,37,-23,-13,-3,-13,-3,7,17,-33,27);//=21	
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,37,17,-33,-23,-13,-3,-13,-3,7,27);//=21	

    if ( print_gpc_error_ ) {
        outfile->Printf("\n");
        outfile->Printf("    ==> Generalized Pauli Constraint Errors <===\n");
        outfile->Printf("\n");
        for (int i = saveoff; i < saveoff+n_gpc_/n_gpc_states_; i++) {
            outfile->Printf("    %5i %20.12lf %20.12lf %5s\n",i,A_p[i],b->pointer()[i],A_p[i] <= b->pointer()[i] ? "" : "XXX");
        }
    }

}

}
