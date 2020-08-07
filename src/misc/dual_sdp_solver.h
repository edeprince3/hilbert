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

#ifndef DUAL_SDP_SOLVER_H
#define DUAL_SDP_SOLVER_H

#include <psi4/libmints/vector.h>
#include <psi4/liboptions/liboptions.h>

#include <lbfgs.h>

#include <misc/sdp_solver.h>

using namespace psi;

namespace hilbert {

class DualSDPSolver: public SDPSolver {

  public:

    /// DualSDPSolver constructor
    DualSDPSolver(long int n_primal, long int n_dual, Options & options);

    /// DualSDPSolver destructor
    ~DualSDPSolver();

    /// solve the sdp problem
    void solve(std::shared_ptr<Vector> x, 
               std::shared_ptr<Vector> b, 
               std::shared_ptr<Vector> c,
               std::vector<int> primal_block_dim,
               int maxiter,
               SDPCallbackFunction evaluate_Au,
               SDPCallbackFunction evaluate_ATu,
               void * data);

    /// solve the sdp problem (low rank)
    void solve_low_rank(std::shared_ptr<Vector> x, 
                        std::shared_ptr<Vector> b, 
                        std::shared_ptr<Vector> c,
                        std::vector<int> primal_block_dim,
                        std::vector<int> primal_block_rank,
                        int maxiter,
                        SDPCallbackFunction evaluate_Au,
                        SDPCallbackFunction evaluate_ATu,
                        void * data);

    double evaluate_gradient_x(const lbfgsfloatval_t * r, lbfgsfloatval_t * g);
    double evaluate_gradient_z(const lbfgsfloatval_t * r, lbfgsfloatval_t * g);

    void evaluate_AATu(std::shared_ptr<Vector> AATu,std::shared_ptr<Vector> u);

    void set_iiter(int iiter) { iiter_ = iiter; }

  protected:

    /// pointer to input data
    void * data_;

    /// container for auxiliary variables
    lbfgsfloatval_t * lbfgs_vars_z_;

    /// the number of inner (cg) iterations
    int iiter_;

    /// the penalty parameter (primal)
    double mu_primal_;

    /// the penalty parameter (dual)
    double mu_;

    /// the penalty parameter ( x.z )
    double mu_xz_;

    /// the lagrange multiplier for x.z
    double lag_xz_;

    /// pointer to the input c vector
    std::shared_ptr<Vector> c_;

    /// pointer to the input x vector
    std::shared_ptr<Vector> x_;

    /// pointer to the input b vector
    std::shared_ptr<Vector> b_;

    /// the right-hand side of AATy = A(c - z) + mu(b - Ax) 
    std::shared_ptr<Vector> cg_rhs_;

    /// temporary container 
    std::shared_ptr<Vector> tmp_;

    /// copy of Au callback function
    SDPCallbackFunction evaluate_Au_;

    /// copy of ATu callback function
    SDPCallbackFunction evaluate_ATu_;

    /// copy of list of block sizes
    std::vector<int> primal_block_dim_;

    /// copy of list of block ranks
    std::vector<int> primal_block_rank_;

    /// build x from auxiliary parameters
    void build_x(double * r);

    /// build z from auxiliary parameters
    void build_z(double * r);

};

}

#endif
