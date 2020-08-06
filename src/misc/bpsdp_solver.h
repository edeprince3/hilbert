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

#ifndef BPSDP_SOLVER_H
#define BPSDP_SOLVER_H

#include <psi4/libmints/vector.h>
#include <psi4/liboptions/liboptions.h>

#include <misc/sdp_solver.h>

using namespace psi;

namespace hilbert {

class BPSDPSolver: public SDPSolver{

  public:

    /// BPSDPSolver constructor
    BPSDPSolver(long int n_primal, long int n_dual, Options & options);

    /// BPSDPSolver destructor
    ~BPSDPSolver();

    /// solve the sdp problem
    void solve(std::shared_ptr<Vector> x, 
               std::shared_ptr<Vector> b, 
               std::shared_ptr<Vector> c,
               std::vector<int> primal_block_dim,
               int maxiter,
               SDPCallbackFunction evaluate_Au,
               SDPCallbackFunction evaluate_ATu,
               void * data);

    void evaluate_AATu(std::shared_ptr<Vector> AATu,std::shared_ptr<Vector> u);

/*
    /// evaluate gradient of (x.z)^2 + ||ATy-c+z||^2
    double evaluate_gradient_z(const lbfgsfloatval_t * r, lbfgsfloatval_t * g);
    void set_lbfgs_iter(int iter) { lbfgs_iter_ = iter; }
*/

    // augmented primal error vector
    std::shared_ptr<Vector> ATAx_minus_ATb(std::shared_ptr<Vector> x,
                                           std::shared_ptr<Vector> b,
                                           SDPCallbackFunction evaluate_Au,
                                           SDPCallbackFunction evaluate_ATu,
                                           void * data);

    // dual error vector
    std::shared_ptr<Vector> ATy_plus_z_minus_c(std::shared_ptr<Vector> c,
                                               SDPCallbackFunction evaluate_ATu,
                                               void * data);

  protected:

    /// the right-hand side of AATy = A(c - z) + mu(b - Ax) 
    std::shared_ptr<Vector> cg_rhs_;

    /// update the primal (x) and dual (z) solutions
    void Update_xz(std::shared_ptr<Vector> x, 
                   std::shared_ptr<Vector> c, 
                   std::vector<int> primal_block_dim, 
                   SDPCallbackFunction evaluate_ATu,
                   void * data);

};

}

#endif
