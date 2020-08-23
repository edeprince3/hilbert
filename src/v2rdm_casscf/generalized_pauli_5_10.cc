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
void v2RDMSolver::Generalized_Pauli_5_10_constraints_ATu(SharedVector A,SharedVector u, int state){

    double* A_p = A->pointer();
    double* u_p = u->pointer();
    double ** orb_p = NatOrbs_[state]->pointer();

    int * x1aoff = d1aoff;
    int * x1boff = d1boff;

    int off = gpcoff[state][0];

    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-1,1,0,0,0,0,0,0,0,0);//<=0,    #################
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,-1,1,0,0,0,0,0,0,0);//<=0,    ##             ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,-1,1,0,0,0,0,0,0);//<=0,    ##             ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,0,-1,1,0,0,0,0,0);//<=0,    ##  ordering   ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,0,0,-1,1,0,0,0,0);//<=0,    ##     for     ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,0,0,0,-1,1,0,0,0);//<=0,    ## lambda[i]'s ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,0,0,0,0,-1,1,0,0);//<=0,    ##             ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,0,0,0,0,0,-1,1,0);//<=0,    ##             ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,0,0,0,0,0,0,0,-1,1);//<=0,    ##             ##
    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,1,1,1,1,1,1,1,1,-9);//<=5,    #################

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,1,1,1,-1,0,0,0,-1,-1);//<=2,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,1,0,0,0,1,-1,-1,-1,0);//<=2,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,1,1,0,0,1,-1,0,-1,-1);//<=2,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,1,1,0,-1,0,-1,0,-1,0);//<=2,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,1,0,1,-1,0,0,-1,-1,0);//<=2,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,0,1,0,1,0,1,0,-1,-1,-1);//<=2,

    //GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,9,-1,-1,-1,-1,-1,-1,-1,-1,-1);//<=5,   ##  Pauli inequality:lambda[1]<=1 ##

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-3,2,2,2,2,2,2,-3,-3,-3);//<=5,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,3,3,-2,-2,-2,-2,-2,-2,3);//<=5,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,3,5,-3,-1,-3,-1,1,-5,1);//<=7,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,3,5,1,-5,-3,-1,-3,-1,1);//<=7,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-1,1,3,5,-1,1,3,-5,-3,-3);//<=7,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,3,5,-3,-1,1,-5,-3,-1,1);//<=7,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-1,5,-1,1,3,1,3,-5,-3,-3);//<=7,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-1,5,-1,1,1,3,3,-3,-5,-3);//<=7,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-1,1,3,1,3,5,-1,-5,-3,-3);//<=7,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,5,3,1,-5,-3,-3,-1,-1,1);//<=7,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-1,1,1,3,3,5,-1,-3,-5,-3);//<=7,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,5,3,-3,-3,-1,-1,1,-5,1);//<=7,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,3,1,5,-1,5,-1,-5,-5,-3);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,5,-1,5,-1,3,1,-5,-5,-3);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,5,-1,5,-3,1,1,-5,-5,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,5,1,5,-5,1,-1,-5,-3,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,5,-1,3,1,5,-1,-5,-5,-3);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,5,5,1,-5,-1,-3,1,-5,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,5,5,1,-5,1,-5,-1,-3,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,5,1,5,-5,-1,1,-3,-5,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,3,5,-1,1,5,-5,-1,-3,-5);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,5,3,-1,1,5,-5,-1,-5,-3);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,3,5,1,-1,5,-5,-1,-5,-3);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,5,-1,1,1,5,-3,-5,-5,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,5,1,1,-1,5,-5,-5,-3,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,5,5,-1,-3,1,-5,1,-5,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,5,5,-1,-1,3,-5,1,-5,-3);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,5,3,1,5,-5,-1,1,-5,-3,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,5,5,3,-5,-1,-1,1,-5,-3);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,3,5,5,-5,1,-1,-1,-5,-3);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,5,1,-1,1,5,-5,-3,-5,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,5,3,5,-5,-1,1,-1,-5,-3);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,5,3,5,-3,-1,1,-5,1,-5,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,5,-1,5,-1,1,3,-5,-3,-5);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,5,3,-1,1,1,5,-5,-3,-5,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,5,3,5,-5,-1,-1,1,-3,-5);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,3,5,5,-5,-1,1,-1,-3,-5);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,5,3,1,-1,1,5,-5,-5,-3,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,5,1,3,5,-5,1,-3,-1,-5,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,5,1,3,-1,5,-5,-3,-1,-5);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,5,1,3,5,-5,1,-1,-5,-3,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,3,5,1,-1,5,-5,-3,-1,-5);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,5,3,5,-5,1,-3,-1,-1,-5);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,5,1,1,3,-1,5,-5,-3,-5,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,5,1,3,1,-1,5,-5,-5,-3,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,3,5,5,-5,1,-1,-3,-1,-5);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,5,3,5,1,-5,-3,-1,1,-5,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,5,-1,1,3,5,-1,-5,-3,-5);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,1,3,5,-1,5,-1,-5,-3,-5);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,5,3,5,1,-5,1,-5,-3,-1,-1);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,1,3,5,-1,5,-3,-1,-5,-5);//<=9,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,5,5,1,3,-5,1,-5,-3,-1,-1);//<=9,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,5,7,9,1,-9,-7,-5,-3,-1,3);//<=13,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,5,7,9,-7,-5,-3,-1,1,-9,3);//<=13,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-3,9,-1,1,3,5,7,-9,-7,-5);//<=13,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-3,1,3,5,7,9,-1,-9,-7,-5);//<=13,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,4,9,-1,9,-6,-1,4,-11,-6,-1);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-1,4,9,-1,4,9,-6,-1,-11,-6);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,9,9,4,-6,-6,-1,-1,4,-11,-1);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-1,4,9,9,-6,-1,-1,4,-6,-11);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,9,-1,4,4,-1,9,-11,-6,-6,-1);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-1,4,9,9,-6,-1,4,-1,-11,-6);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-1,9,4,9,-6,-1,-1,4,-11,-6);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,9,4,-1,-1,4,9,-11,-6,-6,-1);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,9,9,4,4,-11,-6,-6,-1,-1,-1);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-1,4,4,9,-1,9,-6,-1,-11,-6);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,4,9,-1,4,-1,9,-11,-6,-6,-1);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,4,9,-1,-1,4,9,-6,-11,-6,-1);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-1,9,-1,4,4,9,-1,-11,-6,-6);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-1,4,4,9,-1,9,-1,-11,-6,-6);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-1,9,-1,9,-1,4,4,-11,-6,-6);//<=15,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,11,-4,1,1,6,6,-4,-9,-9);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,6,11,-4,1,6,-9,1,-9,-4);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,6,11,1,6,-9,-4,1,-9,-4,1);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,11,6,-4,1,1,6,-9,-9,-4,1);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,6,6,11,-9,1,-4,-4,1,-9);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,6,11,1,-4,1,6,-9,-9,-4,1);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,6,11,-4,1,1,6,-9,-4,-9,1);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,6,6,11,-9,-4,1,1,-4,-9);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,6,11,1,6,-9,1,-9,-4,-4,1);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,6,6,11,-9,1,-4,1,-9,-4);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,6,11,6,-9,-4,1,1,-9,-4);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,6,6,11,-4,-4,1,-9,1,-9,1);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,1,1,1,6,6,11,-4,-4,-9,-9);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,6,6,11,1,-9,1,-9,-4,-4,1);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,6,6,11,1,-9,-4,-4,1,-9,1);//<=15,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-3,7,7,7,-3,7,-3,-3,-3,-13);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,7,-3,7,7,-3,7,-3,-3,-13,-3);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,7,7,-3,-3,7,7,-3,-3,-13,-3);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,7,7,7,-3,-3,7,-13,-3,-3,-3);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,7,7,-3,7,-3,-3,7,-3,-13,-3);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,7,7,-3,7,-3,7,-3,-13,-3,-3);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,7,7,7,7,-13,-3,-3,-3,-3,-3);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,7,7,7,-3,-3,-3,-3,7,-13,-3);//<=15,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,13,3,3,-7,3,-7,-7,3,-7);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,13,3,3,3,-7,3,-7,-7,-7,3);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,13,3,3,-7,-7,3,3,-7,-7);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,13,3,-7,3,3,-7,3,-7,-7);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,3,3,3,3,13,-7,-7,-7,-7);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,3,13,3,-7,3,-7,3,-7,-7);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,13,-7,3,3,3,3,-7,-7,-7);//<=15,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,3,3,13,-7,3,3,-7,-7,-7);//<=15,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-2,8,3,8,3,13,-2,-12,-12,-7);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-2,13,-2,8,3,8,3,-12,-12,-7);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-2,8,3,13,-2,8,3,-12,-12,-7);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-2,13,-2,3,8,3,8,-7,-12,-12);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-2,3,8,8,3,13,-7,-2,-12,-12);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-2,8,3,3,8,13,-2,-12,-7,-12);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-2,8,3,13,-2,3,8,-12,-7,-12);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-2,3,8,3,8,13,-2,-7,-12,-12);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-2,3,8,13,-2,3,8,-7,-12,-12);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-2,3,8,13,-2,8,3,-12,-7,-12);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-2,13,-2,8,3,3,8,-12,-7,-12);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-2,13,-2,3,8,8,3,-12,-7,-12);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-2,8,8,13,-7,-2,3,3,-12,-12);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-2,3,8,8,3,13,-2,-12,-7,-12);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,8,8,-2,3,3,13,-12,-12,-7,-2);//<=20,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,12,7,12,-8,-3,-3,-8,2,-13,2);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,12,7,12,-3,-8,-8,-3,2,-13,2);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,12,12,7,-8,-3,-8,-3,2,-13,2);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,12,12,2,7,-13,-3,-8,-8,-3,2);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,12,7,12,2,-13,-8,-3,-3,-8,2);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,12,7,12,-8,-3,2,-13,-3,-8,2);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,12,12,7,2,-13,-8,-3,-8,-3,2);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,12,12,7,-8,-3,2,-13,-8,-3,2);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,12,7,12,-3,-8,2,-13,-8,-3,2);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,12,12,-3,-3,2,7,-13,-8,-8,2);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,12,7,12,2,-13,-3,-8,-8,-3,2);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,7,12,12,-3,-8,2,-13,-3,-8,2);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,7,12,12,2,-13,-3,-8,-3,-8,2);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,7,12,12,-3,-8,-3,-8,2,-13,2);//<=20,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,2,7,12,12,-13,-3,-3,2,-8,-8);//<=20,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,3,8,13,18,-17,-7,-2,3,-12,-7);//<=25,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,13,18,-7,-2,3,8,-17,-12,-7,3);//<=25,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,13,18,8,-12,-7,-7,-2,3,-17,3);//<=25,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,13,18,8,3,-17,-12,-7,-7,-2,3);//<=25,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,13,18,3,8,-17,-7,-12,-7,-2,3);//<=25,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,7,12,-3,2,7,17,-18,-13,-8,-3);//<=25,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-3,7,12,17,-8,-3,2,7,-18,-13);//<=25,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-3,2,7,12,7,17,-8,-3,-18,-13);//<=25,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-3,17,-3,2,7,7,12,-8,-18,-13);//<=25,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-3,2,7,7,12,17,-3,-8,-18,-13);//<=25,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,17,17,7,-3,-13,-13,-13,-3,-3,7);//<=25,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,17,17,7,-13,-13,-3,-3,-3,-13,7);//<=25,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-7,13,3,3,3,13,13,-7,-17,-17);//<=25,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-7,3,3,13,13,13,3,-7,-17,-17);//<=25,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-3,27,-8,2,7,12,17,-23,-18,-13);//<=35,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-3,2,7,12,17,27,-8,-23,-18,-13);//<=35,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,13,18,23,-17,-12,-7,-2,8,-27,3);//<=35,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,13,18,23,8,-27,-17,-12,-7,-2,3);//<=35,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-7,3,13,13,23,33,-7,-27,-27,-17);//<=45,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,-7,33,-7,3,13,13,23,-27,-27,-17);//<=45,

    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,17,27,27,7,-33,-23,-13,-13,-3,7);//<=45,
    GP_N_10_ATu(u_p[offset++],off,A_p,x1aoff,x1boff,orb_p,17,27,27,-23,-13,-13,-3,7,-33,7);//<=45,

}

// portion of A.x corresponding to generalized paulit constraints
void v2RDMSolver::Generalized_Pauli_5_10_constraints_Au(SharedVector A,SharedVector u, int state){

    int saveoff = offset;

    double* A_p = A->pointer();
    double* u_p = u->pointer();
    double ** orb_p = NatOrbs_[state]->pointer();

    int * x1aoff = d1aoff;
    int * x1boff = d1boff;

    double * eigvals = (double*)malloc(10*sizeof(double));
    for (int i = 0; i < 10; i++) {
        eigvals[i] = Generalized_Pauli_Au_term(orb_p,u_p,x1aoff,x1boff,i+1);
    }

    int off = gpcoff[state][0];
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-1,1,0,0,0,0,0,0,0,0);//<=0,    #################
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,-1,1,0,0,0,0,0,0,0);//<=0,    ##             ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,-1,1,0,0,0,0,0,0);//<=0,    ##             ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,0,-1,1,0,0,0,0,0);//<=0,    ##  ordering   ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,-1,1,0,0,0,0);//<=0,    ##     for     ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,0,-1,1,0,0,0);//<=0,    ## lambda[i]'s ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,0,0,-1,1,0,0);//<=0,    ##             ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,0,0,0,-1,1,0);//<=0,    ##             ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,0,0,0,0,0,0,0,-1,1);//<=0,    ##             ##
    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,1,1,1,1,1,1,1,1,-9);//<=5,    #################

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,1,1,1,-1,0,0,0,-1,-1);//<=2,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,1,0,0,0,1,-1,-1,-1,0);//<=2,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,1,1,0,0,1,-1,0,-1,-1);//<=2,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,1,1,0,-1,0,-1,0,-1,0);//<=2,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,1,0,1,-1,0,0,-1,-1,0);//<=2,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,0,1,0,1,0,1,0,-1,-1,-1);//<=2,

    //A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,9,-1,-1,-1,-1,-1,-1,-1,-1,-1);//<=5,   ##  Pauli inequality:lambda[1]<=1 ##

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-3,2,2,2,2,2,2,-3,-3,-3);//<=5,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,3,3,-2,-2,-2,-2,-2,-2,3);//<=5,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,3,5,-3,-1,-3,-1,1,-5,1);//<=7,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,3,5,1,-5,-3,-1,-3,-1,1);//<=7,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-1,1,3,5,-1,1,3,-5,-3,-3);//<=7,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,3,5,-3,-1,1,-5,-3,-1,1);//<=7,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-1,5,-1,1,3,1,3,-5,-3,-3);//<=7,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-1,5,-1,1,1,3,3,-3,-5,-3);//<=7,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-1,1,3,1,3,5,-1,-5,-3,-3);//<=7,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,5,3,1,-5,-3,-3,-1,-1,1);//<=7,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-1,1,1,3,3,5,-1,-3,-5,-3);//<=7,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,5,3,-3,-3,-1,-1,1,-5,1);//<=7,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,3,1,5,-1,5,-1,-5,-5,-3);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,5,-1,5,-1,3,1,-5,-5,-3);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,5,-1,5,-3,1,1,-5,-5,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,5,1,5,-5,1,-1,-5,-3,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,5,-1,3,1,5,-1,-5,-5,-3);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,5,5,1,-5,-1,-3,1,-5,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,5,5,1,-5,1,-5,-1,-3,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,5,1,5,-5,-1,1,-3,-5,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,3,5,-1,1,5,-5,-1,-3,-5);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,5,3,-1,1,5,-5,-1,-5,-3);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,3,5,1,-1,5,-5,-1,-5,-3);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,5,-1,1,1,5,-3,-5,-5,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,5,1,1,-1,5,-5,-5,-3,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,5,5,-1,-3,1,-5,1,-5,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,5,5,-1,-1,3,-5,1,-5,-3);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,5,3,1,5,-5,-1,1,-5,-3,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,5,5,3,-5,-1,-1,1,-5,-3);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,3,5,5,-5,1,-1,-1,-5,-3);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,5,1,-1,1,5,-5,-3,-5,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,5,3,5,-5,-1,1,-1,-5,-3);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,5,3,5,-3,-1,1,-5,1,-5,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,5,-1,5,-1,1,3,-5,-3,-5);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,5,3,-1,1,1,5,-5,-3,-5,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,5,3,5,-5,-1,-1,1,-3,-5);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,3,5,5,-5,-1,1,-1,-3,-5);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,5,3,1,-1,1,5,-5,-5,-3,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,5,1,3,5,-5,1,-3,-1,-5,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,5,1,3,-1,5,-5,-3,-1,-5);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,5,1,3,5,-5,1,-1,-5,-3,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,3,5,1,-1,5,-5,-3,-1,-5);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,5,3,5,-5,1,-3,-1,-1,-5);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,5,1,1,3,-1,5,-5,-3,-5,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,5,1,3,1,-1,5,-5,-5,-3,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,3,5,5,-5,1,-1,-3,-1,-5);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,5,3,5,1,-5,-3,-1,1,-5,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,5,-1,1,3,5,-1,-5,-3,-5);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,1,3,5,-1,5,-1,-5,-3,-5);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,5,3,5,1,-5,1,-5,-3,-1,-1);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,1,3,5,-1,5,-3,-1,-5,-5);//<=9,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,5,5,1,3,-5,1,-5,-3,-1,-1);//<=9,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,5,7,9,1,-9,-7,-5,-3,-1,3);//<=13,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,5,7,9,-7,-5,-3,-1,1,-9,3);//<=13,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-3,9,-1,1,3,5,7,-9,-7,-5);//<=13,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-3,1,3,5,7,9,-1,-9,-7,-5);//<=13,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,4,9,-1,9,-6,-1,4,-11,-6,-1);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-1,4,9,-1,4,9,-6,-1,-11,-6);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,9,9,4,-6,-6,-1,-1,4,-11,-1);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-1,4,9,9,-6,-1,-1,4,-6,-11);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,9,-1,4,4,-1,9,-11,-6,-6,-1);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-1,4,9,9,-6,-1,4,-1,-11,-6);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-1,9,4,9,-6,-1,-1,4,-11,-6);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,9,4,-1,-1,4,9,-11,-6,-6,-1);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,9,9,4,4,-11,-6,-6,-1,-1,-1);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-1,4,4,9,-1,9,-6,-1,-11,-6);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,4,9,-1,4,-1,9,-11,-6,-6,-1);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,4,9,-1,-1,4,9,-6,-11,-6,-1);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-1,9,-1,4,4,9,-1,-11,-6,-6);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-1,4,4,9,-1,9,-1,-11,-6,-6);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-1,9,-1,9,-1,4,4,-11,-6,-6);//<=15,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,11,-4,1,1,6,6,-4,-9,-9);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,6,11,-4,1,6,-9,1,-9,-4);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,6,11,1,6,-9,-4,1,-9,-4,1);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,11,6,-4,1,1,6,-9,-9,-4,1);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,6,6,11,-9,1,-4,-4,1,-9);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,6,11,1,-4,1,6,-9,-9,-4,1);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,6,11,-4,1,1,6,-9,-4,-9,1);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,6,6,11,-9,-4,1,1,-4,-9);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,6,11,1,6,-9,1,-9,-4,-4,1);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,6,6,11,-9,1,-4,1,-9,-4);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,6,11,6,-9,-4,1,1,-9,-4);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,6,6,11,-4,-4,1,-9,1,-9,1);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,1,1,1,6,6,11,-4,-4,-9,-9);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,6,6,11,1,-9,1,-9,-4,-4,1);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,6,6,11,1,-9,-4,-4,1,-9,1);//<=15,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-3,7,7,7,-3,7,-3,-3,-3,-13);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,7,-3,7,7,-3,7,-3,-3,-13,-3);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,7,7,-3,-3,7,7,-3,-3,-13,-3);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,7,7,7,-3,-3,7,-13,-3,-3,-3);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,7,7,-3,7,-3,-3,7,-3,-13,-3);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,7,7,-3,7,-3,7,-3,-13,-3,-3);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,7,7,7,7,-13,-3,-3,-3,-3,-3);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,7,7,7,-3,-3,-3,-3,7,-13,-3);//<=15,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,13,3,3,-7,3,-7,-7,3,-7);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,13,3,3,3,-7,3,-7,-7,-7,3);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,13,3,3,-7,-7,3,3,-7,-7);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,13,3,-7,3,3,-7,3,-7,-7);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,3,3,3,3,13,-7,-7,-7,-7);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,3,13,3,-7,3,-7,3,-7,-7);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,13,-7,3,3,3,3,-7,-7,-7);//<=15,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,3,3,13,-7,3,3,-7,-7,-7);//<=15,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-2,8,3,8,3,13,-2,-12,-12,-7);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-2,13,-2,8,3,8,3,-12,-12,-7);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-2,8,3,13,-2,8,3,-12,-12,-7);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-2,13,-2,3,8,3,8,-7,-12,-12);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-2,3,8,8,3,13,-7,-2,-12,-12);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-2,8,3,3,8,13,-2,-12,-7,-12);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-2,8,3,13,-2,3,8,-12,-7,-12);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-2,3,8,3,8,13,-2,-7,-12,-12);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-2,3,8,13,-2,3,8,-7,-12,-12);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-2,3,8,13,-2,8,3,-12,-7,-12);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-2,13,-2,8,3,3,8,-12,-7,-12);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-2,13,-2,3,8,8,3,-12,-7,-12);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-2,8,8,13,-7,-2,3,3,-12,-12);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-2,3,8,8,3,13,-2,-12,-7,-12);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,8,8,-2,3,3,13,-12,-12,-7,-2);//<=20,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,12,7,12,-8,-3,-3,-8,2,-13,2);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,12,7,12,-3,-8,-8,-3,2,-13,2);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,12,12,7,-8,-3,-8,-3,2,-13,2);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,12,12,2,7,-13,-3,-8,-8,-3,2);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,12,7,12,2,-13,-8,-3,-3,-8,2);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,12,7,12,-8,-3,2,-13,-3,-8,2);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,12,12,7,2,-13,-8,-3,-8,-3,2);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,12,12,7,-8,-3,2,-13,-8,-3,2);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,12,7,12,-3,-8,2,-13,-8,-3,2);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,12,12,-3,-3,2,7,-13,-8,-8,2);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,12,7,12,2,-13,-3,-8,-8,-3,2);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,7,12,12,-3,-8,2,-13,-3,-8,2);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,7,12,12,2,-13,-3,-8,-3,-8,2);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,7,12,12,-3,-8,-3,-8,2,-13,2);//<=20,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,2,7,12,12,-13,-3,-3,2,-8,-8);//<=20,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,3,8,13,18,-17,-7,-2,3,-12,-7);//<=25,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,13,18,-7,-2,3,8,-17,-12,-7,3);//<=25,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,13,18,8,-12,-7,-7,-2,3,-17,3);//<=25,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,13,18,8,3,-17,-12,-7,-7,-2,3);//<=25,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,13,18,3,8,-17,-7,-12,-7,-2,3);//<=25,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,7,12,-3,2,7,17,-18,-13,-8,-3);//<=25,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-3,7,12,17,-8,-3,2,7,-18,-13);//<=25,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-3,2,7,12,7,17,-8,-3,-18,-13);//<=25,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-3,17,-3,2,7,7,12,-8,-18,-13);//<=25,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-3,2,7,7,12,17,-3,-8,-18,-13);//<=25,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,17,17,7,-3,-13,-13,-13,-3,-3,7);//<=25,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,17,17,7,-13,-13,-3,-3,-3,-13,7);//<=25,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-7,13,3,3,3,13,13,-7,-17,-17);//<=25,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-7,3,3,13,13,13,3,-7,-17,-17);//<=25,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-3,27,-8,2,7,12,17,-23,-18,-13);//<=35,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-3,2,7,12,17,27,-8,-23,-18,-13);//<=35,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,13,18,23,-17,-12,-7,-2,8,-27,3);//<=35,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,13,18,23,8,-27,-17,-12,-7,-2,3);//<=35,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-7,3,13,13,23,33,-7,-27,-27,-17);//<=45,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,-7,33,-7,3,13,13,23,-27,-27,-17);//<=45,

    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,17,27,27,7,-33,-23,-13,-13,-3,7);//<=45,
    A_p[offset++] = GP_N_10_Au(off,u_p,x1aoff,x1boff,orb_p,eigvals,17,27,27,-23,-13,-13,-3,7,-33,7);//<=45,

    if ( print_gpc_error_ ) {
        outfile->Printf("\n");        outfile->Printf("    ==> Generalized Pauli Constraint Errors <===\n");
        outfile->Printf("\n");
        for (int i = saveoff; i < saveoff+n_gpc_/n_gpc_states_; i++) {
            outfile->Printf("    %5i %20.12lf %20.12lf %5s\n",i,A_p[i],b->pointer()[i],A_p[i] <= b->pointer()[i] ? "" : "XXX");
        }
    }

    free(eigvals);
}

}
