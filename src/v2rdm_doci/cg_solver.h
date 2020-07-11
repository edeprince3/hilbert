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

#ifndef CG_SOLVER_H
#define CG_SOLVER_H

#include<psi4/libmints/vector.h>


namespace psi{ 

typedef void (*CallbackType)(long int,SharedVector,SharedVector,void *);  

class CGSolver {
public:

    CGSolver(long int n);
    ~CGSolver();
    void preconditioned_solve(long int n,
               SharedVector Ap,
               SharedVector  x,
               SharedVector  b,
               SharedVector  precon,
               CallbackType function, void * data);
    void solve(long int n,
               SharedVector Ap,
               SharedVector  x,
               SharedVector  b,
               CallbackType function, void * data);

    int total_iterations();
    void set_max_iter(int iter);
    void set_convergence(double conv);

private:

    int    n_;
    int    iter_;
    int    cg_max_iter_;
    double cg_convergence_;
    SharedVector p;
    SharedVector r;
    SharedVector z;

};

} // end of namespace

#endif
