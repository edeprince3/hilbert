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

#include <misc/cg_solver.h>

using namespace psi;

namespace hilbert {

typedef void (*BPSDPCallbackFunction)(std::shared_ptr<Vector>,std::shared_ptr<Vector>,void *);

class BPSDPSolver{

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
               BPSDPCallbackFunction evaluate_Au,
               BPSDPCallbackFunction evaluate_ATu,
               CGCallbackFunction evaluate_CG_LHS,
               void * data);

    int iiter_total() { return iiter_total_; }
    int oiter_total() { return oiter_; }
    double iiter_time() { return iiter_time_; }
    double oiter_time() { return oiter_time_; }

    void set_mu(double mu) { mu_ = mu; }
    void set_y(std::shared_ptr<Vector> y) { y_->copy(y.get()); }
    void set_z(std::shared_ptr<Vector> z) { z_->copy(z.get()); }

    double get_mu() { return mu_; }
    std::shared_ptr<Vector> get_y() { return y_; }
    std::shared_ptr<Vector> get_z() { return z_; }

    bool is_converged(){ return is_converged_; }

  protected:

    /// Options object
    Options& options_;

    /// the error in the primal constraints
    double primal_error_;

    /// the error in the dual constraints
    double dual_error_;

    /// is the solver converged?
    bool is_converged_;

    /// the number of outer iterations
    int oiter_;

    /// the total number of inner iterations
    int iiter_total_;

    /// the total time taken by the inner iterations
    double iiter_time_;

    /// the total time taken by the outer iterations
    double oiter_time_;

    /// the penalty parameter
    double mu_;

    /// the dual solution vector
    std::shared_ptr<Vector> y_;

    /// the second dual solution vector
    std::shared_ptr<Vector> z_;

    /// the right-hand side of AATy = A(c - z) + mu(b - Ax) 
    std::shared_ptr<Vector> cg_rhs_;

    /// temporary container the size of Au
    std::shared_ptr<Vector> Au_;

    /// temporary container the size of ATu
    std::shared_ptr<Vector> ATu_;

    /// the requested convergence of the primal dual gap
    double e_convergence_;

    /// the requested convergence of the primal and dual errors
    double r_convergence_;

    /// the dimension of the primal vector
    long int n_primal_;

    /// the dimension of the dual vector
    long int n_dual_;

    /// update the primal (x) and dual (z) solutions
    void Update_xz(std::shared_ptr<Vector> x, 
                   std::shared_ptr<Vector> c, 
                   std::vector<int> primal_block_dim, 
                   BPSDPCallbackFunction evaluate_ATu,
                   void * data);

};

}

#endif
