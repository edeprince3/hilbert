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

// N:=[
//[0,0,-1,1,0,0,0,0,0,0]<=0      #################
//,[0,0,0,0,-1,1,0,0,0,0]<=0     ##             ##
//,[0,0,0,0,0,-1,1,0,0,0]<=0     ##             ##
//,[0,0,0,0,0,0,0,0,-1,1]<=0     ##  ordering   ##
//,[0,0,0,0,0,0,-1,1,0,0]<=0     ##     for     ##
//,[0,-1,1,0,0,0,0,0,0,0]<=0     ## lambda[i]'s ## 
//,[-1,1,0,0,0,0,0,0,0,0]<=0     ##             ##
//,[0,0,0,0,0,0,0,-1,1,0]<=0     ##             ##
//,[0,0,0,-1,1,0,0,0,0,0]<=0     ##             ##
//,[1,1,1,1,1,1,1,1,1,-9]<=4     #################
//
//,[9,-1,-1,-1,-1,-1,-1,-1,-1,-1]<=6  ##  Pauli inequality:lambda[1]<=1 ##
//
//,[1,1,1,-1,1,-1,-1,-1,1,-1]<=2
//
//,[1,1,2,-1,2,-2,-1,-2,0,0]<=3
//,[1,2,1,-1,2,-2,-2,-1,0,0]<=3
//,[2,1,2,-2,1,-2,-1,-1,0,0]<=3
//,[1,2,2,-2,1,-1,-2,-1,0,0]<=3
//
//,[1,2,2,-1,3,-3,-2,-2,0,0]<=4
//,[2,2,3,-3,1,-2,-2,-1,0,0]<=4
//
//,[5,1,3,-3,-1,1,-3,-1,1,-3]<=6
//,[3,1,1,1,5,-3,-3,-3,-1,-1]<=6
//,[5,1,-1,1,3,-3,-3,-1,1,-3]<=6
//,[5,1,3,-3,1,-3,-1,-1,1,-3]<=6
//,[3,5,1,-3,1,-3,-1,-1,1,-3]<=6
//,[3,5,-1,-1,1,-3,1,-3,1,-3]<=6
//,[3,5,1,-3,-1,-1,1,-3,1,-3]<=6
//,[5,-1,1,1,3,-3,-1,-3,1,-3]<=6
//,[3,1,5,-3,1,1,-3,-3,-1,-1]<=6
//,[3,5,1,-3,1,-3,1,-3,-1,-1]<=6
//
//,[4,-1,4,-1,4,-1,-6,-1,-1,-1]<=6
//,[4,-1,-1,4,4,-1,-1,-6,-1,-1]<=6
//,[4,4,-1,-1,-1,-1,-1,-1,4,-6]<=6
//,[4,4,-1,-1,-1,-1,4,-6,-1,-1]<=6
//,[-1,4,4,-1,4,-1,-1,-6,-1,-1]<=6
//,[4,-1,4,-1,-1,4,-1,-6,-1,-1]<=6
//,[4,4,-1,-1,4,-6,-1,-1,-1,-1]<=6
//,[4,4,4,-6,-1,-1,-1,-1,-1,-1]<=6
//
//,[7,2,2,-3,2,-3,-3,2,-3,-3]<=8
//,[7,-3,2,2,2,2,-3,-3,-3,-3]<=8
//,[7,2,2,-3,-3,2,2,-3,-3,-3]<=8
//,[7,2,-3,2,2,-3,2,-3,-3,-3]<=8
//,[2,2,2,2,7,-3,-3,-3,-3,-3]<=8
//,[2,2,7,-3,2,2,-3,-3,-3,-3]<=8
//,[2,7,2,-3,2,-3,2,-3,-3,-3]<=8
//
//,[6,1,1,1,6,-4,-4,-4,1,-4]<=9
//,[6,6,1,-4,1,-4,1,-4,1,-4]<=9
//,[6,1,6,-4,1,1,-4,-4,1,-4]<=9
//
//,[9,4,4,-6,4,-6,-6,-1,-1,-1]<=11
//,[4,4,4,-1,9,-6,-6,-6,-1,-1]<=11
//,[4,9,4,-6,4,-6,-1,-6,-1,-1]<=11
//,[4,4,9,-6,4,-1,-6,-6,-1,-1]<=11
//
//,[3,5,7,-3,9,-9,-5,-7,-1,1]<=12
//,[3,5,7,-3,9,-9,-7,-5,1,-1]<=12
//,[5,7,9,-9,3,-5,-7,-3,-1,1]<=12
//,[5,7,9,-9,3,-7,-5,-3,1,-1]<=12
//,[7,5,9,-9,3,-7,-5,-3,-1,1]<=12
//,[3,7,5,-3,9,-9,-7,-5,-1,1]<=12
//
//,[8,-2,3,3,8,-7,-7,-2,3,-7]<=12
//
//,[11,1,21,-9,1,11,-19,-9,-9,1]<=24
//,[1,11,21,-9,11,1,-19,-9,-9,1]<=24
//,[21,11,-9,-9,1,1,11,-19,1,-9]<=24
//,[21,11,-9,-9,1,1,1,-9,11,-19]<=24
//,[21,11,1,-9,-9,-9,1,1,11,-19]<=24
//,[21,1,11,-9,11,-9,-19,-9,1,1]<=24
//,[21,1,11,-9,11,-19,-9,1,-9,1]<=24
//,[21,-9,1,11,11,1,-19,-9,-9,1]<=24
//,[21,-9,11,1,1,11,-19,-9,-9,1]<=24
//,[21,11,1,-9,11,-19,-9,-9,1,1]<=24
//,[21,11,11,-19,-9,-9,1,1,1,-9]<=24
//,[11,1,1,11,21,-9,-19,-9,-9,1]<=24
//,[11,1,1,11,21,-19,-9,-9,1,-9]<=24
//,[21,11,11,-19,1,-9,-9,-9,1,1]<=24
//,[1,11,11,1,21,-9,-19,-9,-9,1]<=24
//,[11,1,11,1,21,-19,-9,-9,-9,1]<=24
//,[1,11,11,1,21,-19,-9,-9,1,-9]<=24
//
//,[6,6,11,1,21,-14,-14,-9,-4,-4]<=24
//,[11,1,6,6,21,-14,-14,-9,-4,-4]<=24
//,[21,-4,1,6,11,-14,-9,-4,6,-14]<=24
//,[21,6,11,-14,-4,-9,-4,1,6,-14]<=24
//,[21,11,-9,-4,-4,1,6,-14,6,-14]<=24
//,[21,11,6,-14,-9,-4,-4,1,6,-14]<=24
//,[21,6,11,-14,6,-14,-9,-4,1,-4]<=24
//,[21,6,6,-9,11,-14,-14,-4,-4,1]<=24
//,[21,6,11,-14,6,-9,-14,-4,-4,1]<=24
//,[21,11,6,-14,6,-14,-9,-4,-4,1]<=24
//
//,[19,29,9,-21,9,-21,-1,-11,-11,-1]<=36
//,[19,29,-1,-11,9,-21,9,-21,-11,-1]<=36
//,[9,29,9,-11,19,-21,-1,-21,-11,-1]<=36
//,[19,-1,29,-11,9,9,-21,-21,-11,-1]<=36
//,[19,9,29,-21,-1,9,-11,-21,-11,-1]<=36
//,[9,9,29,-11,19,-1,-21,-21,-11,-1]<=36
//,[9,19,29,-21,9,-1,-11,-21,-11,-1]<=36
//,[19,29,9,-21,-1,-11,9,-21,-11,-1]<=36
//,[9,29,19,-21,9,-11,-1,-21,-11,-1]<=36
//,[29,9,-11,9,19,-21,-1,-21,-11,-1]<=36
//,[29,9,19,-21,-11,9,-1,-21,-11,-1]<=36
//,[29,19,-11,-1,9,-21,9,-21,-11,-1]<=36
//,[29,19,-11,-1,9,-21,-11,-1,9,-21]<=36
//,[29,19,9,-21,-11,-1,9,-21,-11,-1]<=36
//,[9,19,29,-21,9,-11,-1,-21,-1,-11]<=36
//,[19,29,9,-21,-11,-1,-1,-11,9,-21]<=36
//,[19,29,9,-21,-11,-1,9,-21,-1,-11]<=36
//,[19,29,-11,-1,9,-21,9,-21,-1,-11]<=36
//,[19,29,-11,-1,-1,-11,9,-21,9,-21]<=36
//,[19,29,-11,-1,9,-21,-1,-11,9,-21]<=36
//,[19,29,9,-21,-1,-11,-11,-1,9,-21]<=36
//,[19,29,-1,-11,9,-21,-11,-1,9,-21]<=36
//,[19,29,-1,-11,-11,-1,9,-21,9,-21]<=36
//,[29,19,-11,-1,-11,-1,9,-21,9,-21]<=36
//,[29,19,9,-21,-11,-1,-11,-1,9,-21]<=36
//,[29,-11,9,9,19,-1,-21,-21,-11,-1]<=36
//,[29,-11,19,-1,9,9,-21,-21,-11,-1]<=36
//,[19,29,9,-21,9,-21,-11,-1,-1,-11]<=36
//,[19,9,29,-21,9,-1,-21,-11,-11,-1]<=36
//,[29,9,19,-21,9,-11,-21,-1,-11,-1]<=36
//,[29,9,9,-11,19,-21,-21,-1,-11,-1]<=36
//,[29,19,9,-21,9,-21,-11,-1,-11,-1]<=36
//,[19,9,-1,9,29,-21,-11,-21,-11,-1]<=36
//,[9,19,9,-1,29,-21,-11,-21,-11,-1]<=36
//,[19,-1,9,9,29,-11,-21,-21,-11,-1]<=36
//,[9,9,19,-1,29,-11,-21,-21,-11,-1]<=36
//,[19,9,9,-1,29,-21,-21,-11,-11,-1]<=36
//,[19,9,-1,9,29,-21,-21,-11,-1,-11]<=36
//,[29,9,19,-21,9,-21,-11,-1,-1,-11]<=36
//,[19,-1,9,9,29,-21,-11,-21,-1,-11]<=36
//,[19,9,29,-21,-1,9,-21,-11,-1,-11]<=36
//,[29,-1,-1,9,19,-21,-11,-11,9,-21]<=36
//,[29,9,19,-21,-1,-11,-11,-1,9,-21]<=36
//,[9,19,9,-1,29,-21,-21,-11,-1,-11]<=36
//,[9,9,19,-1,29,-21,-11,-21,-1,-11]<=36
//,[9,19,29,-21,9,-1,-21,-11,-1,-11]<=36
//
//,[13,23,33,-27,13,-7,-27,-17,-7,3]<=42
//,[13,13,23,-7,33,-27,-17,-27,-7,3]<=42
//,[13,13,23,-7,33,-27,-27,-17,3,-7]<=42
//,[33,13,23,-27,13,-27,-17,-7,-7,3]<=42
//,[13,23,13,-7,33,-27,-27,-17,-7,3]<=42]:


// portion of A^T.y corresponding to generalized paulit constraints
void v2RDMSolver::Generalized_Pauli_4_10_constraints_ATu(SharedVector A,SharedVector u, int state){

    double* A_p = A->pointer();
    double* u_p = u->pointer();
    double ** orb_p = NatOrbs_[state]->pointer();

    //int * x1aoff;
    //int * x1boff;

    //if ( gpc_[state] == GeneralizedPauli_6_10 ) {
    //    x1aoff = q1aoff;
    //    x1boff = q1boff;
    //}else {
    //    x1aoff = d1aoff;
    //    x1boff = d1boff;
    //}

    int off = gpcoff[state][0];

    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,-1,1,0,0,0,0,0,0);//<=0      #################
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,0,0,-1,1,0,0,0,0);//<=0     ##             ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,0,0,0,-1,1,0,0,0);//<=0     ##             ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,0,0,0,0,0,0,-1,1);//<=0     ##  ordering   ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,0,0,0,0,-1,1,0,0);//<=0     ##     for     ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,-1,1,0,0,0,0,0,0,0);//<=0     ## lambda[i]'s ## 
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-1,1,0,0,0,0,0,0,0,0);//<=0     ##             ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,0,0,0,0,0,-1,1,0);//<=0     ##             ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,0,-1,1,0,0,0,0,0);//<=0     ##             ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,1,1,1,1,1,1,1,1,-9);//<=4     #################
         
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,9,-1,-1,-1,-1,-1,-1,-1,-1,-1);//6  ##  Pauli inequality:lambda[1]<=1 ##
        
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,1,1,-1,1,-1,-1,-1,1,-1);//2
       
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,1,2,-1,2,-2,-1,-2,0,0);//<=3
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,2,1,-1,2,-2,-2,-1,0,0);//<=3
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,2,1,2,-2,1,-2,-1,-1,0,0);//<=3
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,2,2,-2,1,-1,-2,-1,0,0);//<=3
      
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,2,2,-1,3,-3,-2,-2,0,0);//<=4
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,2,2,3,-3,1,-2,-2,-1,0,0);//<=4
     
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,5,1,3,-3,-1,1,-3,-1,1,-3);//<=6
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,3,1,1,1,5,-3,-3,-3,-1,-1);//<=6
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,5,1,-1,1,3,-3,-3,-1,1,-3);//<=6
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,5,1,3,-3,1,-3,-1,-1,1,-3);//<=6
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,3,5,1,-3,1,-3,-1,-1,1,-3);//<=6
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,3,5,-1,-1,1,-3,1,-3,1,-3);//<=6
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,3,5,1,-3,-1,-1,1,-3,1,-3);//<=6
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,5,-1,1,1,3,-3,-1,-3,1,-3);//<=6
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,3,1,5,-3,1,1,-3,-3,-1,-1);//<=6
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,3,5,1,-3,1,-3,1,-3,-1,-1);//<=6
    
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,4,-1,4,-1,4,-1,-6,-1,-1,-1);//<=6
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,4,-1,-1,4,4,-1,-1,-6,-1,-1);//<=6
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,4,4,-1,-1,-1,-1,-1,-1,4,-6);//<=6
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,4,4,-1,-1,-1,-1,4,-6,-1,-1);//<=6
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,-1,4,4,-1,4,-1,-1,-6,-1,-1);//<=6
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,4,-1,4,-1,-1,4,-1,-6,-1,-1);//<=6
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,4,4,-1,-1,4,-6,-1,-1,-1,-1);//<=6
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,4,4,4,-6,-1,-1,-1,-1,-1,-1);//<=6
   
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,2,2,-3,2,-3,-3,2,-3,-3);//<=8
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,-3,2,2,2,2,-3,-3,-3,-3);//<=8
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,2,2,-3,-3,2,2,-3,-3,-3);//<=8
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,2,-3,2,2,-3,2,-3,-3,-3);//<=8
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,2,2,2,2,7,-3,-3,-3,-3,-3);//<=8
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,2,2,7,-3,2,2,-3,-3,-3,-3);//<=8
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,2,7,2,-3,2,-3,2,-3,-3,-3);//<=8
  
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,6,1,1,1,6,-4,-4,-4,1,-4);//<=9
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,6,6,1,-4,1,-4,1,-4,1,-4);//<=9
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,6,1,6,-4,1,1,-4,-4,1,-4);//<=9
 
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,4,4,-6,4,-6,-6,-1,-1,-1);//<=11
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,4,4,4,-1,9,-6,-6,-6,-1,-1);//<=11
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,4,9,4,-6,4,-6,-1,-6,-1,-1);//<=11
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,4,4,9,-6,4,-1,-6,-6,-1,-1);//<=11

    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,3,5,7,-3,9,-9,-5,-7,-1,1);//<=12
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,3,5,7,-3,9,-9,-7,-5,1,-1);//<=12
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,5,7,9,-9,3,-5,-7,-3,-1,1);//<=12
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,5,7,9,-9,3,-7,-5,-3,1,-1);//<=12
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,7,5,9,-9,3,-7,-5,-3,-1,1);//<=12
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,3,7,5,-3,9,-9,-7,-5,-1,1);//<=12

    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,8,-2,3,3,8,-7,-7,-2,3,-7);//12

    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,11,1,21,-9,1,11,-19,-9,-9,1);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,11,21,-9,11,1,-19,-9,-9,1);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,11,-9,-9,1,1,11,-19,1,-9);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,11,-9,-9,1,1,1,-9,11,-19);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,11,1,-9,-9,-9,1,1,11,-19);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,1,11,-9,11,-9,-19,-9,1,1);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,1,11,-9,11,-19,-9,1,-9,1);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,-9,1,11,11,1,-19,-9,-9,1);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,-9,11,1,1,11,-19,-9,-9,1);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,11,1,-9,11,-19,-9,-9,1,1);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,11,11,-19,-9,-9,1,1,1,-9);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,11,1,1,11,21,-9,-19,-9,-9,1);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,11,1,1,11,21,-19,-9,-9,1,-9);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,11,11,-19,1,-9,-9,-9,1,1);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,11,11,1,21,-9,-19,-9,-9,1);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,11,1,11,1,21,-19,-9,-9,-9,1);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,1,11,11,1,21,-19,-9,-9,1,-9);//<=24

    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,6,6,11,1,21,-14,-14,-9,-4,-4);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,11,1,6,6,21,-14,-14,-9,-4,-4);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,-4,1,6,11,-14,-9,-4,6,-14);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,6,11,-14,-4,-9,-4,1,6,-14);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,11,-9,-4,-4,1,6,-14,6,-14);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,11,6,-14,-9,-4,-4,1,6,-14);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,6,11,-14,6,-14,-9,-4,1,-4);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,6,6,-9,11,-14,-14,-4,-4,1);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,6,11,-14,6,-9,-14,-4,-4,1);//<=24
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,21,11,6,-14,6,-14,-9,-4,-4,1);//<=24

    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,29,9,-21,9,-21,-1,-11,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,29,-1,-11,9,-21,9,-21,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,29,9,-11,19,-21,-1,-21,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,-1,29,-11,9,9,-21,-21,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,9,29,-21,-1,9,-11,-21,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,9,29,-11,19,-1,-21,-21,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,19,29,-21,9,-1,-11,-21,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,29,9,-21,-1,-11,9,-21,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,29,19,-21,9,-11,-1,-21,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,29,9,-11,9,19,-21,-1,-21,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,29,9,19,-21,-11,9,-1,-21,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,29,19,-11,-1,9,-21,9,-21,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,29,19,-11,-1,9,-21,-11,-1,9,-21);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,29,19,9,-21,-11,-1,9,-21,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,19,29,-21,9,-11,-1,-21,-1,-11);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,29,9,-21,-11,-1,-1,-11,9,-21);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,29,9,-21,-11,-1,9,-21,-1,-11);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,29,-11,-1,9,-21,9,-21,-1,-11);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,29,-11,-1,-1,-11,9,-21,9,-21);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,29,-11,-1,9,-21,-1,-11,9,-21);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,29,9,-21,-1,-11,-11,-1,9,-21);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,29,-1,-11,9,-21,-11,-1,9,-21);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,29,-1,-11,-11,-1,9,-21,9,-21);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,29,19,-11,-1,-11,-1,9,-21,9,-21);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,29,19,9,-21,-11,-1,-11,-1,9,-21);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,29,-11,9,9,19,-1,-21,-21,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,29,-11,19,-1,9,9,-21,-21,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,29,9,-21,9,-21,-11,-1,-1,-11);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,9,29,-21,9,-1,-21,-11,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,29,9,19,-21,9,-11,-21,-1,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,29,9,9,-11,19,-21,-21,-1,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,29,19,9,-21,9,-21,-11,-1,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,9,-1,9,29,-21,-11,-21,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,19,9,-1,29,-21,-11,-21,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,-1,9,9,29,-11,-21,-21,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,9,19,-1,29,-11,-21,-21,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,9,9,-1,29,-21,-21,-11,-11,-1);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,9,-1,9,29,-21,-21,-11,-1,-11);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,29,9,19,-21,9,-21,-11,-1,-1,-11);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,-1,9,9,29,-21,-11,-21,-1,-11);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,19,9,29,-21,-1,9,-21,-11,-1,-11);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,29,-1,-1,9,19,-21,-11,-11,9,-21);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,29,9,19,-21,-1,-11,-11,-1,9,-21);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,19,9,-1,29,-21,-21,-11,-1,-11);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,9,19,-1,29,-21,-11,-21,-1,-11);//<=36
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,9,19,29,-21,9,-1,-21,-11,-1,-11);//<=36

    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,13,23,33,-27,13,-7,-27,-17,-7,3);//<=42
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,13,13,23,-7,33,-27,-17,-27,-7,3);//<=42
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,13,13,23,-7,33,-27,-27,-17,3,-7);//<=42
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,33,13,23,-27,13,-27,-17,-7,-7,3);//<=42
    GP_N_10_ATu(u_p[offset++],off,A_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],orb_p,13,23,13,-7,33,-27,-27,-17,-7,3);//<=42]:
}

// portion of A.x corresponding to generalized paulit constraints
void v2RDMSolver::Generalized_Pauli_4_10_constraints_Au(SharedVector A,SharedVector u, int state){

    int saveoff = offset;

    double* A_p = A->pointer();
    double* u_p = u->pointer();
    double ** orb_p = NatOrbs_[state]->pointer();

    //int * x1aoff;
    //int * x1boff;

    //if ( gpc_[state] == GeneralizedPauli_6_10 ) {
    //    x1aoff = q1aoff;
    //    x1boff = q1boff;
    //}else {
    //    x1aoff = d1aoff;
    //    x1boff = d1boff;
    //}

    double * eigvals = (double*)malloc(10*sizeof(double));
    for (int i = 0; i < 10; i++) {
        eigvals[i] = Generalized_Pauli_Au_term(orb_p,u_p,gpc_rdm_map_a_[state],gpc_rdm_map_b_[state],i+1);
    }

    int off = gpcoff[state][0];
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,-1,1,0,0,0,0,0,0);//<=0      #################
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,-1,1,0,0,0,0);//<=0     ##             ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,0,-1,1,0,0,0);//<=0     ##             ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,0,0,0,0,-1,1);//<=0     ##  ordering   ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,0,0,-1,1,0,0);//<=0     ##     for     ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,-1,1,0,0,0,0,0,0,0);//<=0     ## lambda[i]'s ## 
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-1,1,0,0,0,0,0,0,0,0);//<=0     ##             ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,0,0,0,-1,1,0);//<=0     ##             ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,0,-1,1,0,0,0,0,0);//<=0     ##             ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,1,1,1,1,1,1,1,1,-9);//<=4     #################

    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,9,-1,-1,-1,-1,-1,-1,-1,-1,-1);//6  ##  Pauli inequality:lambda[1]<=1 ##

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,1,1,1,-1,1,-1,-1,-1,1,-1);//2

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,1,1,2,-1,2,-2,-1,-2,0,0);//<=3
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,1,2,1,-1,2,-2,-2,-1,0,0);//<=3
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,2,1,2,-2,1,-2,-1,-1,0,0);//<=3
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,1,2,2,-2,1,-1,-2,-1,0,0);//<=3

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,1,2,2,-1,3,-3,-2,-2,0,0);//<=4
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,2,2,3,-3,1,-2,-2,-1,0,0);//<=4

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,5,1,3,-3,-1,1,-3,-1,1,-3);//<=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,3,1,1,1,5,-3,-3,-3,-1,-1);//<=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,5,1,-1,1,3,-3,-3,-1,1,-3);//<=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,5,1,3,-3,1,-3,-1,-1,1,-3);//<=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,3,5,1,-3,1,-3,-1,-1,1,-3);//<=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,3,5,-1,-1,1,-3,1,-3,1,-3);//<=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,3,5,1,-3,-1,-1,1,-3,1,-3);//<=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,5,-1,1,1,3,-3,-1,-3,1,-3);//<=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,3,1,5,-3,1,1,-3,-3,-1,-1);//<=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,3,5,1,-3,1,-3,1,-3,-1,-1);//<=6

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,4,-1,4,-1,4,-1,-6,-1,-1,-1);//<=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,4,-1,-1,4,4,-1,-1,-6,-1,-1);//<=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,4,4,-1,-1,-1,-1,-1,-1,4,-6);//<=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,4,4,-1,-1,-1,-1,4,-6,-1,-1);//<=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,-1,4,4,-1,4,-1,-1,-6,-1,-1);//<=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,4,-1,4,-1,-1,4,-1,-6,-1,-1);//<=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,4,4,-1,-1,4,-6,-1,-1,-1,-1);//<=6
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,4,4,4,-6,-1,-1,-1,-1,-1,-1);//<=6

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,2,2,-3,2,-3,-3,2,-3,-3);//<=8
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,-3,2,2,2,2,-3,-3,-3,-3);//<=8
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,2,2,-3,-3,2,2,-3,-3,-3);//<=8
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,2,-3,2,2,-3,2,-3,-3,-3);//<=8
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,2,2,2,2,7,-3,-3,-3,-3,-3);//<=8
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,2,2,7,-3,2,2,-3,-3,-3,-3);//<=8
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,2,7,2,-3,2,-3,2,-3,-3,-3);//<=8

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,6,1,1,1,6,-4,-4,-4,1,-4);//<=9
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,6,6,1,-4,1,-4,1,-4,1,-4);//<=9
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,6,1,6,-4,1,1,-4,-4,1,-4);//<=9

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,4,4,-6,4,-6,-6,-1,-1,-1);//<=11
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,4,4,4,-1,9,-6,-6,-6,-1,-1);//<=11
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,4,9,4,-6,4,-6,-1,-6,-1,-1);//<=11
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,4,4,9,-6,4,-1,-6,-6,-1,-1);//<=11

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,3,5,7,-3,9,-9,-5,-7,-1,1);//<=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,3,5,7,-3,9,-9,-7,-5,1,-1);//<=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,5,7,9,-9,3,-5,-7,-3,-1,1);//<=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,5,7,9,-9,3,-7,-5,-3,1,-1);//<=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,7,5,9,-9,3,-7,-5,-3,-1,1);//<=12
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,3,7,5,-3,9,-9,-7,-5,-1,1);//<=12

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,8,-2,3,3,8,-7,-7,-2,3,-7);//12

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,11,1,21,-9,1,11,-19,-9,-9,1);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,1,11,21,-9,11,1,-19,-9,-9,1);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,11,-9,-9,1,1,11,-19,1,-9);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,11,-9,-9,1,1,1,-9,11,-19);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,11,1,-9,-9,-9,1,1,11,-19);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,1,11,-9,11,-9,-19,-9,1,1);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,1,11,-9,11,-19,-9,1,-9,1);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,-9,1,11,11,1,-19,-9,-9,1);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,-9,11,1,1,11,-19,-9,-9,1);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,11,1,-9,11,-19,-9,-9,1,1);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,11,11,-19,-9,-9,1,1,1,-9);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,11,1,1,11,21,-9,-19,-9,-9,1);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,11,1,1,11,21,-19,-9,-9,1,-9);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,11,11,-19,1,-9,-9,-9,1,1);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,1,11,11,1,21,-9,-19,-9,-9,1);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,11,1,11,1,21,-19,-9,-9,-9,1);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,1,11,11,1,21,-19,-9,-9,1,-9);//<=24

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,6,6,11,1,21,-14,-14,-9,-4,-4);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,11,1,6,6,21,-14,-14,-9,-4,-4);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,-4,1,6,11,-14,-9,-4,6,-14);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,6,11,-14,-4,-9,-4,1,6,-14);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,11,-9,-4,-4,1,6,-14,6,-14);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,11,6,-14,-9,-4,-4,1,6,-14);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,6,11,-14,6,-14,-9,-4,1,-4);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,6,6,-9,11,-14,-14,-4,-4,1);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,6,11,-14,6,-9,-14,-4,-4,1);//<=24
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,21,11,6,-14,6,-14,-9,-4,-4,1);//<=24

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,29,9,-21,9,-21,-1,-11,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,29,-1,-11,9,-21,9,-21,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,29,9,-11,19,-21,-1,-21,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,-1,29,-11,9,9,-21,-21,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,9,29,-21,-1,9,-11,-21,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,9,29,-11,19,-1,-21,-21,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,19,29,-21,9,-1,-11,-21,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,29,9,-21,-1,-11,9,-21,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,29,19,-21,9,-11,-1,-21,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,29,9,-11,9,19,-21,-1,-21,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,29,9,19,-21,-11,9,-1,-21,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,29,19,-11,-1,9,-21,9,-21,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,29,19,-11,-1,9,-21,-11,-1,9,-21);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,29,19,9,-21,-11,-1,9,-21,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,19,29,-21,9,-11,-1,-21,-1,-11);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,29,9,-21,-11,-1,-1,-11,9,-21);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,29,9,-21,-11,-1,9,-21,-1,-11);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,29,-11,-1,9,-21,9,-21,-1,-11);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,29,-11,-1,-1,-11,9,-21,9,-21);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,29,-11,-1,9,-21,-1,-11,9,-21);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,29,9,-21,-1,-11,-11,-1,9,-21);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,29,-1,-11,9,-21,-11,-1,9,-21);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,29,-1,-11,-11,-1,9,-21,9,-21);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,29,19,-11,-1,-11,-1,9,-21,9,-21);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,29,19,9,-21,-11,-1,-11,-1,9,-21);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,29,-11,9,9,19,-1,-21,-21,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,29,-11,19,-1,9,9,-21,-21,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,29,9,-21,9,-21,-11,-1,-1,-11);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,9,29,-21,9,-1,-21,-11,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,29,9,19,-21,9,-11,-21,-1,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,29,9,9,-11,19,-21,-21,-1,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,29,19,9,-21,9,-21,-11,-1,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,9,-1,9,29,-21,-11,-21,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,19,9,-1,29,-21,-11,-21,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,-1,9,9,29,-11,-21,-21,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,9,19,-1,29,-11,-21,-21,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,9,9,-1,29,-21,-21,-11,-11,-1);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,9,-1,9,29,-21,-21,-11,-1,-11);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,29,9,19,-21,9,-21,-11,-1,-1,-11);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,-1,9,9,29,-21,-11,-21,-1,-11);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,19,9,29,-21,-1,9,-21,-11,-1,-11);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,29,-1,-1,9,19,-21,-11,-11,9,-21);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,29,9,19,-21,-1,-11,-11,-1,9,-21);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,19,9,-1,29,-21,-21,-11,-1,-11);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,9,19,-1,29,-21,-11,-21,-1,-11);//<=36
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,9,19,29,-21,9,-1,-21,-11,-1,-11);//<=36

    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,13,23,33,-27,13,-7,-27,-17,-7,3);//<=42
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,13,13,23,-7,33,-27,-17,-27,-7,3);//<=42
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,13,13,23,-7,33,-27,-27,-17,3,-7);//<=42
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,33,13,23,-27,13,-27,-17,-7,-7,3);//<=42
    A_p[offset++] = GP_N_10_Au(off,u_p,orb_p,eigvals,13,23,13,-7,33,-27,-27,-17,-7,3);//<=42]:

    if ( print_gpc_error_ ) {
        outfile->Printf("\n");        outfile->Printf("    ==> Generalized Pauli Constraint Errors <===\n");
        outfile->Printf("\n");
        for (int i = saveoff; i < saveoff+n_gpc_/n_gpc_states_; i++) {
            outfile->Printf("    %5i %20.12lf %20.12lf %5s\n",i,A_p[i],b->pointer()[i],A_p[i] <= b->pointer()[i] ? "" : "XXX");
        }
    }

    free(eigvals);
}

void v2RDMSolver::GP_N_10_ATu(double dum,int & off, double * A, int *** map_a, int *** map_b, double ** orbs, 
        int d1, int d2, int d3, int d4, int d5, int d6, int d7,int d8, int d9, int d10) {
    A[off++] = dum;
    Generalized_Pauli_ATu_term(dum *  d1,orbs,A,map_a,map_b,1);
    Generalized_Pauli_ATu_term(dum *  d2,orbs,A,map_a,map_b,2);
    Generalized_Pauli_ATu_term(dum *  d3,orbs,A,map_a,map_b,3);
    Generalized_Pauli_ATu_term(dum *  d4,orbs,A,map_a,map_b,4);
    Generalized_Pauli_ATu_term(dum *  d5,orbs,A,map_a,map_b,5);
    Generalized_Pauli_ATu_term(dum *  d6,orbs,A,map_a,map_b,6);
    Generalized_Pauli_ATu_term(dum *  d7,orbs,A,map_a,map_b,7);
    Generalized_Pauli_ATu_term(dum *  d8,orbs,A,map_a,map_b,8);
    Generalized_Pauli_ATu_term(dum *  d9,orbs,A,map_a,map_b,9);
    Generalized_Pauli_ATu_term(dum * d10,orbs,A,map_a,map_b,10);

}

double v2RDMSolver::GP_N_10_Au(int & off, double * u, double ** orbs, double * eigvals,
        int d1, int d2, int d3, int d4, int d5, int d6, int d7,int d8, int d9, int d10) {

    double dum = u[off++]
               + d1  * eigvals[0]
               + d2  * eigvals[1]
               + d3  * eigvals[2]
               + d4  * eigvals[3]
               + d5  * eigvals[4]
               + d6  * eigvals[5]
               + d7  * eigvals[6]
               + d8  * eigvals[7]
               + d9  * eigvals[8]
               + d10 * eigvals[9];
    return dum;
}

}
