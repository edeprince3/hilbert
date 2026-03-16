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

#ifndef CS_SOLVER_H
#define CS_SOLVER_H

#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsio/psio.hpp>

#include <lbfgs.h>

#include <vector>

using namespace psi;

namespace hilbert{

class cs_solver {

  public:

    cs_solver(Options & options);
    ~cs_solver();

    void solve(long int nalpha,
               long int nbeta,
               //double eref,
               std::vector<double> T,
               std::vector<double> super_phi,
               std::vector<double> reference_rho,
               std::vector<std::vector<double> > grid);

    // minimize non-interacting kinetic energy s.t. density recovers the reference density
    double constrained_search();

    // evaluate augmented lagrangian function and gradient
    double evaluate_augmented_lagrangian();

    // pointer to variables (public for returning Ca/Cb)
    double * vars_;

    // pointer to gradient of variables (public for lbfgs)
    double * dvars_;

    // number of lbfgs iterations (public for lbfgs)
    int iiter_;

  protected:

    /// options dictionary
    Options & options_;

    /// kinetic energy
    std::vector<double> T_;

    /// number of alpha electrons
    long int nalpha_;

    /// number of beta electrons
    long int nbeta_;

    /// number of orbitals
    long int nso_;

    /// convergence in the constraints
    double r_convergence_;

    /// convergence in the spatial density
    double rho_convergence_;

    /// convergence in the non-interacting kinetic energy
    double e_convergence_;

    /// maximum number of lbfgs iterations
    long int lbfgs_maxiter_;

    /// do check correctness of derivatives?
    bool do_check_derivatives_;

    /// number of variables
    long int n_;

    /// number of constraints
    long int nconstraints_;

    /// non-interacting kinetic energy
    double energy_;

    /// trace of d1a
    double trd1a_;

    /// trace of d1b
    double trd1b_;

    /// number of grid points_
    long int phi_points_;

    /// phi matrix
    std::vector<double> super_phi_;

    /// additional buffer of the size of the phi matrix
    std::vector<double> temp_phi_;

    /// reference spatial density
    std::vector<double> reference_rho_;

    /// grid x values
    std::vector<double> grid_x_;

    /// grid y values
    std::vector<double> grid_y_;

    /// grid z values
    std::vector<double> grid_z_;

    /// grid weights
    std::vector<double> grid_w_;

    double * D1a_, *R1da_, * dR1da_;
    double * D1b_, *R1db_, * dR1db_;

    double * temp1_, * temp2_;

    /// lagrange multipliers
    double * lambda_;

    /// penalty parameter for augmented lagrangian
    double mu_;

    /// target constraint values
    double * c_;

    /// actual constraint values
    double * cval_;

    /// errors in constraints
    double * cerror_;

    /// offset in constraint list
    int offset_;

    /// reference kinetic energy
    //double eref_;

    /// set constraint values
    void set_constraints();

    /// choose starting solutions
    void random_guess();

    /// build d1
    void build_d1();

    /// build density, evaluate constraint int [rho(r)-rho0(r)]^2 dr = 0
    void evaluate_integral_rho_squared();

    /// build density, evaluate constraint rho(r)-rho0(r) for all r
    void evaluate_rho_pointwise();

    /// evaluate the energy and the gradient of the energy
    double kinetic_energy();

    /// evaluate orbital orthonormality conditions
    void evaluate_orbital_orthonormality();

    /// numerical check of derivatives of augmented lagrangian function
    void check_derivatives(lbfgsfloatval_t * vars);
};

}


#endif
