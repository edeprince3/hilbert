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

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <psi4/libciomr/libciomr.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>
#include <psi4/libplugin/plugin.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libqt/qt.h>
#include <psi4/libpsi4util/process.h>
#include <psi4/libpsi4util/PsiOutStream.h>
#include <psi4/libpsio/psio.hpp>
#include <complex>

#include "nonsym_davidson_solver_qed.h"
#include "qed_blas.h"
#include "omp.h"

using namespace std;
using namespace psi;

#define BIGNUM 1E100
#define MAXIT 1000

namespace psi {
    Nonsym_DavidsonSolver_QED::Nonsym_DavidsonSolver_QED() {
    }

    Nonsym_DavidsonSolver_QED::~Nonsym_DavidsonSolver_QED() {
    }

    SharedMatrix Nonsym_DavidsonSolver_QED::solve(double *Adiag, size_t N, size_t M, double *reval, double **rer, double **rel,
                                              const BuildSigma &build_sigma, size_t maxdim, size_t init_dim, size_t maxiter,
                                              double residual_error,
                                              double residual_norm, bool use_residual_norm, bool read_guess, double shift) {
        int min_pos;
        double minimum;

        bool applyShift = fabs(shift) > 1e-6;

        std::shared_ptr<Matrix> B(new Matrix(maxdim, maxdim));//Q^T(A-B)(A+B)Q
        std::shared_ptr<Matrix> Q(new Matrix(maxdim, N));
        std::shared_ptr<Matrix> Qn(new Matrix(maxdim, N));
        std::shared_ptr<Matrix> sigmar(new Matrix(N, maxdim));//(A-B)(A+B)q
        std::shared_ptr<Matrix> sigmal(new Matrix(N, maxdim));//(A-B)(A+B)q
        std::shared_ptr<Matrix> Cl(new Matrix(maxdim, maxdim));
        std::shared_ptr<Matrix> Cr(new Matrix(maxdim, maxdim));
        std::shared_ptr<Matrix> Rl(new Matrix(M + 1, N));
        std::shared_ptr<Matrix> Rr(new Matrix(M + 1, N));
        std::shared_ptr<Vector> lambda(new Vector(maxdim));
        std::shared_ptr<Vector> lambdai(new Vector(maxdim));
        std::shared_ptr<Vector> lambda_old(new Vector(M + 1));
        std::shared_ptr<Vector> lambdai_old(new Vector(M + 1));

        double **Clp = Cl->pointer();
        double **Crp = Cr->pointer();
        double **Qp = Q->pointer();
        double **Bp = B->pointer();
        double **sigmarp = sigmar->pointer();
        double **sigmalp = sigmal->pointer();
        double **Qnp = Qn->pointer();
        double **Rlp = Rl->pointer();
        double **Rrp = Rr->pointer();
        double *lambdap = lambda->pointer();
        double *lambdaip = lambdai->pointer();
        double *lambda_oldp = lambda_old->pointer();
        double *lambdai_oldp = lambdai_old->pointer();
        //Adiag guess
        auto *Adiag2 = (double *) malloc(N * sizeof(double));
        C_DCOPY(N, Adiag, 1, Adiag2, 1);
        if (!applyShift) {
            //Construct guess
            for (int i = 0; i < init_dim; i++) {
                minimum = Adiag2[0];
                min_pos = 0;
                for (int j = 1; j < N; j++) {
                    if (Adiag2[j] < minimum) {
                        minimum = Adiag2[j];
                        min_pos = j;
                    }
                }
                Qp[i][min_pos] = 1.0;
                Adiag2[min_pos] = BIGNUM;
//            lambdap[i] = minimum;
            }
        } else {
            for (int i = 0; i < N; i++) {
                Adiag2[i] = fabs(Adiag2[i] - shift);
            }
            //use unit vectors as guess
            for (int i = 0; i < init_dim; i++) {
                minimum = Adiag2[0];
                min_pos = 0;
                for (int j = 1; j < N; j++) {
                    if (Adiag2[j] < minimum) {
                        minimum = Adiag2[j];
                        min_pos = j;
                    }
                }
                Qp[i][min_pos] = 1.0;
                Adiag2[min_pos] = BIGNUM;
//            lambdap[i] = minimum;
            }
        }
        free(Adiag2);
        int L = (int)init_dim;

        if (read_guess) {
            //use current eigenvectors as guess
            for (size_t i = 0; i < M; i++) {
                for (size_t k = 0; k < N; k++) {
                    Qnp[2 * i][k] = rer[i][k];
                    Qnp[2 * i + 1][k] = rel[i][k];
                }
            }
        }

        int iter = 0;
        int maxit = maxiter;
        bool convergence = false;
        outfile->Printf("\n\n");
        outfile->Printf("  ==> Non-Symmetric Davidson Solver <==\n\n");

        while (iter < maxit) {


            //this code only checks convergence of the real part of the eigenvalues
            double start = omp_get_wtime();
            for (int i = 0; i < M + 1; i++) {
                lambda_oldp[i] = lambdap[i];
                lambdai_oldp[i] = lambdaip[i];
            }

            //make all elements of eigenvectors and eigenvalues zero
            for (int i = 0; i < maxdim; i++) {
                lambdap[i] = 0.0;
                lambdaip[i] = 0.0;
                for (int j = 0; j < maxdim; j++) {
                    Clp[i][j] = 0.0;
                    Crp[i][j] = 0.0;

                }
            }


            //make all elements of QAQ^T zero
            for (int i = 0; i < maxdim; i++) {
                for (int j = 0; j < maxdim; j++) {
                    Bp[i][j] = 0.0;
                }
            }
            double end = omp_get_wtime();

            build_sigma(N, maxdim, L, Qp, sigmarp, sigmalp);
            start = omp_get_wtime();

            //Complex Value test
            /*std::shared_ptr<Matrix> C(new Matrix(N, N));
            double **Cp = C->pointer();
            for (int i = 0; i < N; i++) {
                Cp[0][i] = 0.5;
                Cp[N - 1][i] = 0.5;
                Cp[i][0] = 0.4;
                Cp[i][N - 1] = 0.4;
                Cp[i][i] = 0.2 * i + 1.0;
            }
            for (int i = 0; i < N - 1; i++) {
                Cp[i][i + 1] = -1.0;
                Cp[i + 1][i] = 1.0;
            }

            for (int i = 0; i < L; i++) {
                for (int j = 0; j < N; j++) {
                    double dum1 = 0.0;
                    double dum2 = 0.0;
                    for (int k = 0; k < N; k++) {
                        dum1 += Cp[j][k] * Qp[i][k];
                        dum2 += Cp[k][j] * Qp[i][k];
                    }
                    sigmarp[j][i] = dum1;
                    sigmalp[j][i] = dum2;
                }
            }*/

            //projection of A onto subspace
            for (int i = 0; i < L; i++) {
                for (int j = 0; j < L; j++) {
                    double dum1 = 0.0;
                    for (int k = 0; k < N; k++) {
                        dum1 += sigmarp[k][j] * Qp[i][k];
                    }
                    Bp[i][j] = dum1;
                }
            }
            bool ascending = true;
            NonSymmetricEigenvalueEigenvector(L, Bp, lambdap, lambdaip, Clp, Crp, ascending, applyShift, shift);

            for (size_t k = 0; k < M + 1; k++) {
                for (size_t n = 0; n < N; n++) {
                    Rrp[k][n] = 0.0;
                    Rlp[k][n] = 0.0;
                }
            }
            //residual
            outfile->Printf("%5s %20s %20s\n", "State", "Real Part", "Imaginary Part");
            size_t complex_root = 0;
            for (size_t k = 0; k < M; k++) {
                if (fabs(lambdaip[k]) < 1e-16) {
                    for (size_t n = 0; n < N; n++) {
                        double dum1 = 0.0;
                        double dum2 = 0.0;
                        for (size_t m = 0; m < L; m++) {
                            dum1 += (sigmarp[n][m] - lambdap[k] * Qp[m][n]) * Crp[k][m];
                            dum2 += (sigmalp[n][m] - lambdap[k] * Qp[m][n]) * Clp[k][m];

                        }
                        Rrp[k][n] = dum1;
                        Rlp[k][n] = dum2;
                    }
                    outfile->Printf("%5d %20.12lf --------------------\n", k, lambdap[k], lambdaip[k]);
                } else {
                    complex_root++;
                    if (complex_root % 2 == 0) continue;

                    double re = lambdap[k];
                    double im = lambdaip[k];
                    outfile->Printf("%5d %20.12lf %20.12lf\n", k, lambdap[k], lambdaip[k]);
                    for (size_t n = 0; n < N; n++) {
                        double res_rr = 0.0;
                        double res_ri = 0.0;
                        double res_lr = 0.0;
                        double res_li = 0.0;
                        for (size_t m = 0; m < L; m++) {
                            // f_R = AQR-wQR
                            res_rr += sigmarp[n][m] * Crp[k][m] - (re * Crp[k][m] - im * Crp[k + 1][m]) * Qp[m][n];
                            res_ri += sigmarp[n][m] * Crp[k + 1][m] - (im * Crp[k][m] + re * Crp[k + 1][m]) * Qp[m][n];
                            // f_L = A^T QL - (w*) QL
                            res_lr += sigmalp[n][m] * Clp[k][m] - (re * Clp[k][m] + im * Clp[k + 1][m]) * Qp[m][n];
                            res_li += sigmalp[n][m] * Clp[k + 1][m] - (-im * Clp[k][m] + re * Clp[k + 1][m]) * Qp[m][n];
                        }
                        Rrp[k][n] = res_rr;
                        Rrp[k + 1][n] = res_ri;
                        Rlp[k][n] = res_lr;
                        Rlp[k + 1][n] = res_li;
                    }
                }
            }

            size_t Mn;
            if (complex_root % 2 != 0) {
                Mn = M + 1;
            } else Mn = M;
            vector<int> index;
            int check_residual = 0;
            int conv = 0;
            auto residual_norms = (double *) calloc((M + 1) * N, sizeof(double));
            outfile->Printf("residual norm, dimension %d\n", L);
            for (int k = 0; k < Mn; k++) {
                double normvalr = C_DDOT(N, Rrp[k], 1, Rrp[k], 1);
                double normvall = C_DDOT(N, Rlp[k], 1, Rlp[k], 1);
                double error = lambdap[k] - lambda_oldp[k];
                normvalr = sqrt(normvalr);
                normvall = sqrt(normvall);
                if ((normvalr < residual_norm) && (normvall < residual_norm) &&
                    ((fabs(error) < residual_error) || applyShift)) {
                    conv++;
                } else {
                    //find rows of residual that have norm larger than some threshold
                    index.push_back(k);
                }
                outfile->Printf("residual norm%d %20.12lf %20.12lf\n", k, normvalr, normvall);

                for (int n = 0; n < N; n++) {
                    residual_norms[k*N+n] = normvalr*normvall;
                }
            }

            outfile->Printf("check_convergence %d\n", conv);
            outfile->Printf("energy error, dimension %d\n", L);
            for (int k = 0; k < Mn; k++) {
                double error = lambdap[k] - lambda_oldp[k];
                outfile->Printf("eig_error%d %20.12lf \n", k, error);
            }
            if (conv == Mn) {
                convergence = true;
                outfile->Printf("energy and eigenvectors converged\n");
                break;
            }

            //preconditioner
            complex_root = 0;
            for (size_t k = 0; k < M; k++) {
                if (fabs(lambdaip[k]) < 1e-16) {
                    for (size_t n = 0; n < N; n++) {
                        double dum = lambdap[k] - Adiag[n];
                        if (fabs(dum) > 1e-6) {
                            Rrp[k][n] /= dum;
                            Rlp[k][n] /= dum;
                        } else {
                            Rrp[k][n] = 0.0;
                            Rlp[k][n] = 0.0;
                        }
                    }
                    //outfile->Printf("%20.12lf %20.12lf\n", lambdap[k], lambdaip[k]);
                } else {
                    complex_root++;
                    if (complex_root % 2 == 0) continue;
                    double re = lambdap[k];
                    double im = lambdaip[k];
                    //outfile->Printf("%20.12lf %20.12lf\n", re,im);
                    for (size_t n = 0; n < N; n++) {
                        double md = (re - Adiag[n]) * (re - Adiag[n]) + im * im;
                        double res_rr = ((re - Adiag[n]) * Rrp[k][n] + im * Rrp[k + 1][n]) / md;
                        double res_lr = ((re - Adiag[n]) * Rlp[k][n] - im * Rlp[k + 1][n]) / md;
                        double res_ri = ((re - Adiag[n]) * Rrp[k + 1][n] - im * Rrp[k][n]) / md;
                        double res_li = ((re - Adiag[n]) * Rlp[k + 1][n] + im * Rlp[k][n]) / md;
                        Rrp[k][n] = res_rr;
                        Rlp[k][n] = res_lr;
                        Rrp[k + 1][n] = res_ri;
                        Rlp[k + 1][n] = res_li;
                    }
                }
            }
            if (use_residual_norm) {
                outfile->Printf(
                        "\n Use residual norm to control the number of correction vectors added to new subspace\n");
                outfile->Printf("dimension %d, number of complex roots %d\n", L, complex_root);
                if (complex_root % 2 != 0) {
                    Mn = M + 1;
                    outfile->Printf("need %d vectors (for imaginary part)\n", Mn);
                } else Mn = M;
                if (maxdim - L < 2 * index.size()) {
                    outfile->Printf("Subspace too large: maxdim = %d, L = %d\n", maxdim, L);
                    outfile->Printf("Collapsing eigenvectors.\n");
                    for (size_t i = 0; i < maxdim; i++) {
                        for (size_t k = 0; k < N; k++) {
                            Qnp[i][k] = 0.0;
                        }
                    }
                    for (size_t i = 0; i < Mn; i++) {
                        for (size_t k = 0; k < N; k++) {
                            for (size_t j = 0; j < L; j++) {
                                Qnp[2 * i][k] += Clp[i][j] * Qp[j][k];
                                Qnp[2 * i + 1][k] += Crp[i][j] * Qp[j][k];
                            }
                        }
                    }
                    for (size_t i = 0; i < index.size(); i++) {
                        for (size_t k = 0; k < N; k++) {
                            Qnp[2 * i + 2 * Mn][k] = Rlp[index[i]][k];
                            Qnp[2 * i + 1 + 2 * Mn][k] = Rrp[index[i]][k];
                        }
                    }
                    //zero current subspace bases?
                    for (size_t i = 0; i < maxdim; i++) {
                        for (size_t k = 0; k < N; k++) {
                            Qp[i][k] = 0.0;
                        }
                    }

                    //Lapack Householder
                    //copy new vectors into place
                    for (size_t i = 0; i < 2 * (Mn + index.size()); i++) {
                        for (size_t k = 0; k < N; k++) {
                            Qp[i][k] = Qnp[i][k];
                        }
                    }

                    L = 2 * (Mn + index.size());
                    qr_orthogonalization(maxdim, N, L, Qp);
                } else {
                    for (size_t i = 0; i < maxdim; i++) {
                        for (size_t k = 0; k < N; k++) {
                            Qnp[i][k] = 0.0;
                        }
                    }
                    for (size_t k = 0; k < index.size(); k++) {
                        for (size_t n = 0; n < N; n++) {
                            Qnp[2 * k][n] = Rlp[index[k]][n];
                            Qnp[2 * k + 1][n] = Rrp[index[k]][n];
                        }
                    }

                    //Lapack Householder
                    for (size_t k = 0; k < index.size(); k++) {
                        for (size_t n = 0; n < N; n++) {
                            Qp[L + k][n] = Rlp[index[k]][n];
                            Qp[L + index.size() + k][n] = Rrp[index[k]][n];
                        }
                    }
                    L = L + 2 * index.size();
                    qr_orthogonalization(maxdim, N, L, Qp);
                }
            } else {
                if (complex_root % 2 != 0) {
                    Mn = M + 1;
                    outfile->Printf("need %d vectors (for imaginary part)\n", Mn);
                } else Mn = M;


                if (maxdim - L < 2 * Mn) {
                    outfile->Printf("Subspace too large: maxdim = %d, L = %d\n", maxdim, L);
                    outfile->Printf("Collapsing eigenvectors.\n");
                    //zero current subspace bases?
                    for (size_t i = 0; i < maxdim; i++) {
                        for (size_t k = 0; k < N; k++) {
                            Qnp[i][k] = 0.0;
                        }
                    }

                    for (size_t i = 0; i < Mn; i++) {
                        for (size_t k = 0; k < N; k++) {
                            for (size_t j = 0; j < L; j++) {
                                Qnp[2 * i][k] += Clp[i][j] * Qp[j][k];
                                Qnp[2 * i + 1][k] += Crp[i][j] * Qp[j][k];
                            }
                            Qnp[2 * i + 2 * Mn][k] = Rlp[i][k];
                            Qnp[2 * i + 1 + 2 * Mn][k] = Rrp[i][k];
                        }
                    }
                    for (size_t i = 0; i < maxdim; i++) {
                        for (size_t k = 0; k < N; k++) {
                            Qp[i][k] = 0.0;
                        }
                    }

                    //Lapack Householder
                    // copy new vectors into place
                    for (size_t i = 0; i < 4 * Mn; i++) {
                        for (size_t k = 0; k < N; k++) {
                            Qp[i][k] = Qnp[i][k];
                        }
                    }
                    // Q->print();

                    L = 4 * Mn;
                    qr_orthogonalization(maxdim, N, L, Qp);
                } else {
                    for (size_t i = 0; i < maxdim; i++) {
                        for (size_t k = 0; k < N; k++) {
                            Qnp[i][k] = 0.0;
                        }
                    }
                    for (size_t k = 0; k < Mn; k++) {
                        for (size_t n = 0; n < N; n++) {
                            Qnp[2 * k][n] = Rlp[k][n];
                            Qnp[2 * k + 1][n] = Rrp[k][n];
                        }
                    }

                    //Lapack Householder
                    for (size_t k = 0; k < Mn; k++) {
                        for (size_t n = 0; n < N; n++) {
                            Qp[L + k][n] = Rrp[k][n];
                            Qp[L + Mn + k][n] = Rlp[k][n];
                        }
                    }
                    L = L + 2 * Mn;
                    qr_orthogonalization(maxdim, N, L, Qp);
                    //Q->print();
                }
            }

            outfile->Printf("\nIter %5d / %-5d\n", ++iter, maxiter);
            writeVectors(N, M, maxdim, reval, rer, rel, Clp, Crp, Qp, lambdap, lambdaip, L, false);
        }
        outfile->Printf("\n");
        if (!convergence) {
            outfile->Printf("Davidson procedure fails after %d iterations, try increase maxdim", maxit);
        }
        outfile->Printf("\n");

        //get the eigenvectors from desired roots; even if the Davidson procedure fails

        outfile->Printf("\n\n############### Final Eigenvalues ###############\n");
        outfile->Printf("%5s %20s %20s\n", "State", "Real Part", "Imaginary Part");
        writeVectors((int)N, (int)M, maxdim, reval, rer, rel, Clp, Crp, Qp, lambdap, lambdaip, L, true);

        return Q;
    }

    /**
     * @brief Write the eigenvectors to container
     *
     * @param N Dimension of the matrix
     * @param M Number of roots
     * @param maxdim Maximum dimension of the subspace
     * @param reval Right eigenvalues
     * @param rer Right eigenvectors
     * @param rel Left eigenvectors
     * @param Clp Coefficients of left eigenvectors
     * @param Crp Coefficients of right eigenvectors
     * @param Qp Subspace vectors
     * @param lambdap Real part of eigenvalues
     * @param lambdaip Imaginary part of eigenvalues
     * @param L Current dimension of the subspace
     * @param print Whether to print the eigenvectors
     */
    int Nonsym_DavidsonSolver_QED::writeVectors(int N, int M, int maxdim, double *reval, double **rer, double **rel,
                                            double *const *Clp, double *const *Crp, double *const *Qp,
                                            const double *lambdap, const double *lambdaip, int L, bool doPrint) const {
        int M_i = M;
        int num_i = 0;
        for (int k = 0; k < M && M <= maxdim; ++k) {
            bool isComplex = false;
            if (fabs(lambdaip[k]) >= 1e-16) {
                ++M;
                ++k;
                ++num_i;
                isComplex = true;
            }
            int l = k - num_i;
            if (isComplex && doPrint) outfile->Printf("%5d %20.12lf %20.12lf\n", l, lambdap[k - 1], lambdaip[k - 1]);
            if (doPrint) outfile->Printf("%5d %20.12lf --------------------\n", l, lambdap[k], lambdaip[k]);
            if (!isComplex) {
                reval[l] = lambdap[k];
                for (int n = 0; n < N; n++) {
                    rer[l][n] = 0.0;
                    rel[l][n] = 0.0;
                    for (int m = 0; m < L; m++) {
                        rer[l][n] += Qp[m][n] * Crp[k][m];
                        rel[l][n] += Qp[m][n] * Clp[k][m];
                    }
                }
            } else {
//                std::complex evalConj1(lambdap[k - 1], lambdaip[k - 1]);
//                std::complex evalConj2(lambdap[k], lambdaip[k]);

//                // project out the real part of the eigenvalues by rotating the magnitude of the complex vector
//                // to the real axis with the nearest sign based on the phase factor.
//                std::complex md = sqrt(evalConj1 * evalConj2);
                reval[l] = lambdap[k]; //copysign(md.real(), evalConj1.real());
                for (int n = 0; n < N; n++) {
                    rer[l][n] = 0.0;
                    rel[l][n] = 0.0;
                    std::complex<double> evec_r;
                    std::complex<double> evec_l;
                    for (int m = 0; m < L; m++) {
                        std::complex dumr(Qp[m][n] * Crp[k - 1][m] * lambdap[k - 1],
                                          Qp[m][n] * Crp[k - 1][m] * lambdaip[k - 1]);
                        std::complex duml(Qp[m][n] * Clp[k][m] * lambdap[k], Qp[m][n] * Clp[k][m] * lambdaip[k]);

                        evec_r += dumr;
                        evec_l += duml;
                    }

                    // Project out the real part of the eigenvectors by rotating the complex vector's magnitude
                    // to the real axis, using the nearest sign based on the phase factor.
                    // This preserves the eigenvectors' sign without requiring phase correction.
                    // note: this will cause disagreements with traditional methods that only use the real part.
                    // Both methods are valid for the purposes of this code.
                    std::complex<double> dumr_mag = sqrt(conj(evec_r) * evec_r);
                    std::complex<double> duml_mag = sqrt(conj(evec_l) * evec_l);
                    rer[l][n] =  copysign(dumr_mag.real(), evec_r.real());
                    rel[l][n] = -copysign(duml_mag.real(), (evec_l).real());
                }
            }
        }
        if (M > maxdim){
            outfile->Printf("WARNING: %d complex eigenvalues encountered. Subspace is insufficient to store them and provide the requested roots.\n", num_i);
        }
        M = M_i;
        return M;
    }

    void Nonsym_DavidsonSolver_QED::eigenvectors_swap(double **L, double **R, int l, int n, int m) {
        //this function is tested with simple systems, further tests may be required
        //l=row l of matrix of eigenvectors where degeneracy happens
        //n=number of degeneracy
        //m=length of each eigenvectors
        //selection sort
        double Smax, S, t;
        int i, j, jmax;
        for (j = 0; j < n - 1; j++) {
            Smax = C_DDOT(m, L[l + j], 1, R[l + j], 1);
            jmax = j;
            // find maximum eigenvalue
            for (i = j + 1; i < n; i++) {
                S = C_DDOT(m, L[l + j], 1, R[l + i], 1);
                outfile->Printf("%d %d %20.12lf %20.12lf \n", j, i, Smax, S);
                if (fabs(Smax) - fabs(S) < 0e0) {
                    jmax = i;
                    Smax = S;
                }
            }
            if (jmax != j) {
                //  outfile->Printf("detect wrong order\n");
                for (i = 0; i < m; i++) {
                    t = R[l + j][i];
                    R[l + j][i] = R[l + jmax][i];
                    R[l + jmax][i] = t;
                }
            }
        }
    }

    void Nonsym_DavidsonSolver_QED::schmidt_biorthogonal(double **A, double **B, size_t m, size_t n) {
        //assuming that A and B have m linear independent row vectors and A[i]*B[i] is non-zero, this will biorthogonalize the rows of A and B
        double RValue1;
        double RValue2;
        for (size_t i = 0; i < m; i++) {
            RValue1 = C_DDOT(n, A[i], 1, B[i], 1);
            for (size_t I = 0; I < n; I++) {
                if (RValue1 < 0.0) {
                    B[i][I] /= -sqrt(fabs(RValue1));
                    A[i][I] /= sqrt(fabs(RValue1));
                } else {
                    B[i][I] /= sqrt(RValue1);
                    A[i][I] /= sqrt(RValue1);
                }

                // B[i][I] /= RValue1;
            }

            for (size_t j = i + 1; j < m; j++) {
                RValue1 = C_DDOT(n, A[i], 1, B[j], 1);
                RValue2 = C_DDOT(n, B[i], 1, A[j], 1);
                for (size_t I = 0; I < n; I++) {
                    B[j][I] -= RValue1 * B[i][I];
                    A[j][I] -= RValue2 * A[i][I];
                }
            }
        }

    }

/**
  * diagonalize general matrix and keep eigenvectors
  */

    void Nonsym_DavidsonSolver_QED::NonSymmetricEigenvalueEigenvector(long int dim, double **M, double *eigval, double *wi,
                                                                  double **el, double **er, bool ascending, bool applyShift,
                                                                  double shift) {
        long int info;
        char vl = 'V';
        char vr = 'V';
        long int n = dim;
        long int lwork = 4 * n;
        auto *work = (doublereal *) malloc(lwork * sizeof(doublereal));
        //double * wi = (double*)malloc(n*sizeof(double));
        //copy the non-zero part of the matrix we need to solve for eigenvalue to another mini matrix
        std::shared_ptr<Matrix> As(new Matrix(dim, dim));
        std::shared_ptr<Matrix> ers(new Matrix(dim, dim));
        std::shared_ptr<Matrix> els(new Matrix(dim, dim));
        std::shared_ptr<Vector> eigr(new Vector(dim));
        std::shared_ptr<Vector> eigi(new Vector(dim));
        double **Asp = As->pointer();
        double **elsp = els->pointer();
        double **ersp = ers->pointer();
        double *eigrp = eigr->pointer();
        double *eigip = eigi->pointer();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Asp[j][i] = M[i][j];
            }
        }


        psi::fnocc::DGEEV(vl, vr, n, Asp[0], n, eigrp, eigip, elsp[0], n, ersp[0], n, work, lwork, info);

        // perform argument sort by pairing eigenvalues with index
        auto *eig_pair = (std::pair<std::complex<double>, size_t> *) malloc(
                n * sizeof(std::pair<std::complex<double>, size_t>));

        for (size_t i = 0; i < n; i++) {
            if (applyShift) {
                eig_pair[i].first = std::complex<double>(fabs(eigrp[i] - shift), eigip[i]);
            } else {
                eig_pair[i].first = std::complex<double>(eigrp[i], eigip[i]);
            }
            eig_pair[i].second = i;
        }

        // sort by eigenvalue
        std::sort(eig_pair, eig_pair + n,
            [ascending](std::pair<std::complex<double>, size_t> &l, std::pair<std::complex<double>, size_t> &r) {
                if (ascending) {
                    bool moveOrder = l.first.real() < r.first.real();
                    if (!moveOrder) {
                        if (fabs(l.first.real() - r.first.real()) <= 1e-16) {
                            return l.first.imag() > r.first.imag();
                        }
                    }
                    return moveOrder;
                } else {
                    bool moveOrder = l.first.real() > r.first.real();
                    if (!moveOrder) {
                        if (fabs(l.first.real() - r.first.real()) <= 1e-16) {
                            return l.first.imag() < r.first.imag();
                        }
                    }
                    return moveOrder;
                }
            }
        );

        // apply sorted indices to reorder eigensystem
        for (size_t i = 0; i < n; ++i) {

            // get ordered index
            size_t min_i = eig_pair[i].second;

            // reorder eigenvalue
            eigval[i] = eigrp[min_i];
            wi[i] = eigip[min_i];

            // reorder eigenvector
            for (size_t j = 0; j < n; j++) {
                er[i][j] = ersp[min_i][j];
                el[i][j] = elsp[min_i][j];
            }

            if (fabs(wi[i]) >= 1e-16) {
                size_t min_i1 = eig_pair[i+1].second;
                eigval[i+1] = eigrp[min_i1];
                wi[i + 1] = eigip[min_i1];

                // reorder eigenvector
                for (size_t j = 0; j < n; j++) {
                    er[i + 1][j] = ersp[min_i1][j];
                    el[i + 1][j] = elsp[min_i1][j];
                }

                // skip the next index
                ++i;
            }
        }
        free(work);
    }


    void Nonsym_DavidsonSolver_QED::qr_orthogonalization(long int L, long int M, long int N, double **A) {

        //this will orthonormalize N first rows of a matrix A with dimension L x M (N <= L <= M).
        //copy nonzero part of A to a smaller maxtrix before orthogonalization
        std::shared_ptr<Matrix> At(new Matrix(N, M));
        std::shared_ptr<Vector> tau(new Vector(N));
        double **Atp = At->pointer();
        //provide variable L here just in case I need to change type of A later
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                Atp[i][j] = A[i][j];
            }
        }

        double *taup = tau->pointer();
        long int lwork = 8 * N;
        double *work = (doublereal *) malloc(lwork * sizeof(doublereal));
        long int info1;
        long int info2;
        psi::fnocc::dgeqrf(M, N, Atp[0], M, taup, work, lwork, info1);
        //outfile->Printf("info1 %d\n",info1);
        psi::fnocc::dorgqr(M, N, N, Atp[0], M, taup, work, lwork, info2);
        //outfile->Printf("info2 %d\n",info2);

        //copy orthonormal rows back to A
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                A[i][j] = Atp[i][j];
            }
        }

        free(work);
    }
}
