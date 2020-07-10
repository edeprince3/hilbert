/*
 *@BEGIN LICENSE
 *
 * pp2RDM, a plugin to:
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

#ifndef PP2RDM_SOLVER_H
#define PP2RDM_SOLVER_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include<bitset>
#include <psi4/libiwl/iwl.h>
#include <psi4/libplugin/plugin.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libqt/qt.h>

#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>

// TODO: move to psifiles.h
#define PSIF_DCC_QMO          268
#define PSIF_V2RDM_CHECKPOINT 269
#define PSIF_V2RDM_D2AA       270
#define PSIF_V2RDM_D2AB       271
#define PSIF_V2RDM_D2BB       272
#define PSIF_V2RDM_D3AAA      273
#define PSIF_V2RDM_D3AAB      274
#define PSIF_V2RDM_D3BBA      275
#define PSIF_V2RDM_D3BBB      276
#define PSIF_V2RDM_D1A        277
#define PSIF_V2RDM_D1B        278

namespace psi{ namespace pp2rdm{

class pp2RDMSolver: public Wavefunction{
  public:
    pp2RDMSolver(SharedWavefunction reference_wavefunction,Options & options);
    ~pp2RDMSolver();
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
    double evaluate_projection_energy();

    /// evaluate sigma vector
    void evaluate_sigma(double * sigma);

    /// evaluate residual vector
    void evaluate_residual(double * residual);

    /// evaluate residual vector for z equations
    void evaluate_residual_z(double * residual);

    /// 1-RDM
    double * d1_;

    /// 2-RDM (ab)
    double * d2ab_;

    /// 2-RDM (aa)
    double * d2aa_;

    /// number of ci iterations
    int ci_iter_;

    /// solve pp2rdm projection equations
    double pp2rdm_iterations(int & iter);

    /// solve z projection equations
    double z_iterations(int & iter);

    /// minimize pp2rdm energy functional using lbfgs
    double pp2rdm_lbfgs_iterations(int & iter);

    /// minimize pp2rdm energy functional using newton-raphson
    double pp2rdm_newton_raphson_iterations(int & iter);

    /// construct c0
    void Normalization();

    /// construct c0 (for constructing RDMs)
    void Normalization_RDMs();

    /// construct 1- and 2-RDM (assuming seniority zero structure)
    double BuildSeniorityZeroRDMs(bool print);

    /// construct 1- and 2-RDM for pCCD
    double BuildSeniorityZeroRDMs_pCCD(bool print);

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

    /// number of auxiliary basis functions
    long int nQ_;

    /// SCF energy
    double escf_;

    /// nuclear repulsion energy
    double enuc_;

    /// full space of integrals for MO gradient / Hessian, blocked by symmetry
    double * tei_;
    long int tei_dim_;

    /// fock (o)
    double * fo_;

    /// fock (v)
    double * fv_;

    /// (ii|aa)
    double * v_iiaa_;

    /// (ia|ai)
    double * v_iaia_;

    /// (ij|ij)
    double * v_ijij_;

    /// (ab|ab)
    double * v_abab_;

    /// amplitudes
    double * t2_;

    /// z-amplitudes
    double * z2_;

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

    /// do restart from a checkpoint file?
    bool restart_from_checkpoint_file_;

    /// build integrals and fock matrix
    void setup_integrals();

    /// check analytic gradient
    void check_gradient();

    /// check analytic hessian
    void check_hessian();

    /// print 1- and 2-RDM elements
    void print_rdms();

    /// print 1- and 2-electron integrals
    void print_integrals();
    ///Build hessian for HF (test)
    double * Fij;
    double * Fab;
    double *  tei_4index_;
    void BuildFock();
    std::shared_ptr<Matrix> mux_a;
    std::shared_ptr<Matrix> mux_b;
    std::shared_ptr<Matrix> muy_a;
    std::shared_ptr<Matrix> muz_a;
    /// construct hf rdms
    double BuildHFRDMs(bool print);
    ///sorting eigenvectors
    void Sort(double *x, double *d, int n, int isort);
    ///spin-free 2rdm
    double * d2t_;

};

}}
#endif

