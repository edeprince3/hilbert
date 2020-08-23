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

#ifndef RRSDP_SOLVER_H
#define RRSDP_SOLVER_H

#include <psi4/libmints/vector.h>
#include <psi4/liboptions/liboptions.h>

#include <lbfgs.h>

#include <misc/sdp_solver.h>

using namespace psi;

namespace hilbert {

class RRSDPSolver: public SDPSolver {

  public:

    /// RRSDPSolver constructor
    RRSDPSolver(long int n_primal, long int n_dual, Options & options);

    /// RRSDPSolver destructor
    ~RRSDPSolver();

    /// solve the sdp problem
    void solve(std::shared_ptr<Vector> x, 
               std::shared_ptr<Vector> b, 
               std::shared_ptr<Vector> c,
               std::vector<int> primal_block_dim,
               int maxiter,
               SDPCallbackFunction evaluate_Au,
               SDPCallbackFunction evaluate_ATu,
               void * data);

    double evaluate_gradient_x(const lbfgsfloatval_t * r, lbfgsfloatval_t * g);

    void set_iiter(int iiter) { iiter_ = iiter; }

    void set_mu_scale_factor(double mu_scale_factor) { mu_scale_factor_ = mu_scale_factor; }

    void set_mu_reset(bool mu_reset) { mu_reset_ = mu_reset; }

  protected:

    /// do reset penalty parameter each time solve() is called?
    bool mu_reset_;

    /// scaling factor for updating penalty parameter
    double mu_scale_factor_;

    /// pointer to input data
    void * data_;

    /// container for auxiliary variables
    lbfgsfloatval_t * lbfgs_vars_x_;

    /// copy of Au callback function
    SDPCallbackFunction evaluate_Au_;

    /// copy of ATu callback function
    SDPCallbackFunction evaluate_ATu_;

    /// copy of list of block sizes
    std::vector<int> primal_block_dim_;

    /// copy of list of block ranks
    std::vector<int> primal_block_rank_;

    /// the number of inner (lbfgs) iterations
    int iiter_;

    /// pointer to the input c vector
    std::shared_ptr<Vector> c_;

    /// pointer to the input x vector
    std::shared_ptr<Vector> x_;

    /// pointer to the input b vector
    std::shared_ptr<Vector> b_;

    /// build x from auxiliary parameters
    void build_x(double * r);

};

}

#endif
