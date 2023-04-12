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

typedef std::function<void(int,int,int,double**,double**,double**)>BuildSigma;
namespace psi{

// sigma-vector build function
//typedef void (*SigmaBuildFunction)(size_t N, size_t maxdim,double **sigma, double **b, void * data);

class Nonsym_DavidsonSolver_QED{

  public:
      Nonsym_DavidsonSolver_QED();
      ~Nonsym_DavidsonSolver_QED();
      SharedMatrix solve(double *Adiag, size_t N, size_t M, double *reval, double **rer, double **rel,
                         const BuildSigma& build_sigma, size_t maxdim, size_t init_dim, size_t maxiter, double residual_error,
                         double residual_norm, bool use_residual_norm, bool read_guess = false, double shift=0.0);
      void eigenvectors_swap(double **L, double **R, int l, int n, int m);
      void schmidt_biorthogonal(double **A, double **B, size_t m, size_t n);
      static void NonSymmetricEigenvalueEigenvector(long int dim, double **M, double * eigval, double * wi, double ** el, double ** er, bool applyShift = false, double shift = 0);
      void qr_orthogonalization(long int L, long int M, long int N, double **A);
  protected:

    int writeVectors(int N, int M, int maxdim, double *reval, double **rer, double **rel,
                     double *const *Clp, double *const *Crp, double *const *Qp,
                     const double *lambdap, const double *lambdaip, int L, bool doPrint) const;
};

}

#endif
