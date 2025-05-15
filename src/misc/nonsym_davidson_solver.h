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
#include <functional>
#ifndef NONSYM_DAVID_H
#define NONSYM_DAVID_H

#include <misc/blas.h>

typedef std::function<void(int,int,int,double**,double**,double**)>BuildSigma;
typedef std::function<void(size_t,size_t,size_t,double**,double**,double**)>BuildSigma2;

namespace psi{

// sigma-vector build function
//typedef void (*SigmaBuildFunction)(size_t N, size_t maxdim,double **sigma, double **b, void * data);

/**
 *  diagonalize general matrix and keep eigenvectors
 *
 */
void NonSymmetricEigenvalueEigenvector(long int dim, double **M, double * eigval, double * wi, double ** el, double ** er);
void GeneralizedEigenvalueProblem(long int N, double **A, double **B, double **vr, double *eig, bool * energy_is_real);

void qr_orthogonalization(long int L, long int M, long int N, double **A);


class Nonsym_DavidsonSolver{

  public:
      size_t eom_maxiter = 100; // maximum number of iterations
      double eom_e_conv = 1e-8; // convergence threshold for energy
      double eom_r_conv = 1e-8; // convergence threshold for residual

      Nonsym_DavidsonSolver();
      ~Nonsym_DavidsonSolver();
      void solve(double *Adiag, int N, int M, double *reval, double **rer, double **rel, BuildSigma build_sigma, int maxdim, size_t init_dim, double residual_norm, bool use_residual_norm);
      void real_generalized_eigenvalue_problem(double *Hdiag, double *Sdiag, size_t N, size_t M, double *reval, double **rer, BuildSigma2 build_sigma2, size_t maxdim,size_t init_dim, double residual_norm, bool use_residual_norm);
      void eigenvectors_swap(double **L, double **R, int l, int n, int m);
      void schmidt_biorthogonal(double **A, double **B, size_t m, size_t n);
  protected:
      int writeVectors(int N, int M, int maxdim, double *reval, double **rer, double **rel,
                     double *const *Clp, double *const *Crp, double *const *Qp,
                     const double *lambdap, const double *lambdaip, int L, bool doPrint) const;

};

}

#endif
