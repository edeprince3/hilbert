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

#ifndef P2RDM_SOLVER_H
#define P2RDM_SOLVER_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <psi4/libiwl/iwl.h>
#include <psi4/libplugin/plugin.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libqt/qt.h>

#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>

#include <misc/hilbert_psifiles.h>
#include <misc/diis.h>

using namespace psi;

namespace hilbert{

class p2RDMSolver: public Wavefunction{
  public:
    p2RDMSolver(SharedWavefunction reference_wavefunction,Options & options);
    ~p2RDMSolver();
    void common_init();
    double compute_energy();
    double compute_excited_energy();

    /// evaluate gradient of the energy (public so lbfgs can call it)
    void evaluate_gradient(double * gradient);

    /// evaluate ( 1/2 x ) the orbital gradient of the energy (public just because)
    void evaluate_orbital_gradient(std::shared_ptr<Matrix> gradient);

    /// evaluate amplitude-amplitude block of the hessian of the energy (public just because)
    void evaluate_a22(std::shared_ptr<Matrix> a22);

    /// evaluate amplitude-amplitude block of the hessian of the energy hessian of the energy numerically (public just because)
    void evaluate_numerical_hessian(std::shared_ptr<Matrix> hessian);

    /// evaluate orbital-amplitude block of the hessian of the energy (public just because)
    void evaluate_a21_b21(std::shared_ptr<Matrix> a21, std::shared_ptr<Matrix> b21);

    /// evaluate orbital-orbital block of the hessian of the energy (public just because)
    void evaluate_a11_b11(std::shared_ptr<Matrix> a11,std::shared_ptr<Matrix> b11);

    /// evaluate metric matrix
    void evaluate_metric_matrix(std::shared_ptr<Matrix> metric);

    /// evaluate one permutation of hessian element pq,rs
    double evaluate_hessian_pqrs(int p, int q, int r, int s);

    /// evaluate the variational energy (public so lbfgs can call it)
    double evaluate_variational_energy();

    /// set number of ci iterations (public so lbfgs can call it)
    void set_number_of_lbfgs_iterations(int iter);

  protected:

    virtual bool same_a_b_orbs() const { return same_a_b_orbs_; }
    virtual bool same_a_b_dens() const { return same_a_b_dens_; }

    /// orbital lagrangian
    double * orbital_lagrangian_;

    /// evaluate the projection energy
    double evaluate_projection_energy(double * t2);

    /// evaluate residual vector
    void evaluate_residual(double * residual, double * t2);

    /// update amplitudes
    double update_amplitudes(double * residual, double * t2, bool do_diis, std::shared_ptr<DIIS> diis);

    /// 1-RDM
    double * d1_;

    /// 2-RDM (ab)
    double * d2ab_;

    /// 2-RDM (aa)
    double * d2aa_;

    /// number of ci iterations
    int ci_iter_;

    /// solve p2rdm projection equations
    double p2rdm_iterations(int & iter);

    /// construct c0
    void Normalization();

    /// construct c0 (for constructing RDMs)
    void Normalization_RDMs();

    /// construct 1- and 2-RDM
    double BuildRDMs(bool print);

    /// optimize orbitals
    double RotateOrbitals();

    /// are we using density fitting?
    bool is_df_;

    /// convergence in energy
    double e_convergence_;

    /// convergence in amplitudes
    double r_convergence_;

    /// three-index integrals in the MO basis
    double * Qmo_;
    double * Qoo_;
    double * Qov_;
    double * Qvo_;
    double * Qvv_;

    /// number of auxiliary basis functions
    long int nQ_;

    /// SCF energy
    double escf_;

    /// nuclear repulsion energy
    double enuc_;

    /// full space of integrals for MO gradient / Hessian, blocked by symmetry
    double * tei_;
    long int tei_dim_;

    /// fock matrix (occupied-occupied)
    double * foo_;

    /// fock matrix (virtual-virtual)
    double * fvv_;

    /// amplitudes
    double * t2_;

    /// residual
    double * r2_;

    /// normalization coefficients
    double * t0_;

    double * oei_;
    int oei_dim_;
    int d2_dim_;
    double * d2_;

    /// pack active-space spin-blocked density into spatial density
    void PackSpatialDensity();

    /// orbital optimization stuff:

    double * orbopt_transformation_matrix_;
    double * orbopt_data_;
    char * orbopt_outfile_;
    bool orbopt_converged_;

    /// total number of orbital optimization
    long int orbopt_iter_total_;

    /// wall time for orbital optimization
    double orbopt_time_;

    /// maximum number of macroiterations
    int maxiter_;

    /// write full 2RDM to disk
    void WriteTPDM();

    /// mo-mo transformation matrix
    SharedMatrix newMO_;

    /// build integrals and fock matrix
    void setup_integrals();

};

}
#endif

