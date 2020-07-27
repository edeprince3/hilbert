/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
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

#include "davidson_solver.h"

using namespace psi;

namespace hilbert{

#define BIGNUM 1E100
#define MAXIT 1000

DavidsonSolver::DavidsonSolver(){
}
DavidsonSolver::~DavidsonSolver(){
}

int DavidsonSolver::solve(double *Adiag, int N, int M, double *eps, double **v, double cutoff, int print, 
        SigmaBuildFunction build_sigma, HamiltonianElementFunction hamiltonian_element, size_t & iter, void * data, int maxdim) {

    int i, j, k, L, I;
    double minimum;
    int min_pos, numf, converged, skip_check;
    int init_dim;
    double norm, denom, diff;

    //int maxdim = 20 * M;

    // current set of guess vectors, stored by row
    std::shared_ptr<Matrix> b ( new Matrix(maxdim,N) );
    // guess vectors formed from old vectors, stored by row
    std::shared_ptr<Matrix> bnew ( new Matrix(maxdim,N) );
    // sigma vectors, stored by column
    std::shared_ptr<Matrix> sigma ( new Matrix(N,maxdim) );
    // Davidson mini-Hamiltonian
    std::shared_ptr<Matrix> G ( new Matrix(maxdim,maxdim) );
    // residual eigenvectors, stored by row
    std::shared_ptr<Matrix> f ( new Matrix(maxdim,N) );
   // eigenvectors of G
    std::shared_ptr<Matrix> alpha ( new Matrix(maxdim,maxdim) );
    // eigenvalues of G
    std::shared_ptr<Vector> lambda ( new Vector(maxdim) );
    // approximate roots from previous iteration
    std::shared_ptr<Vector> lambda_old ( new Vector(maxdim) );

    double ** b_p         = b->pointer();
    double ** bnew_p      = bnew->pointer();
    double ** sigma_p     = sigma->pointer();
    double ** G_p         = G->pointer();
    double ** f_p         = f->pointer();
    double ** alpha_p     = alpha->pointer();
    double * lambda_p     = lambda->pointer();
    double * lambda_old_p = lambda_old->pointer();

    double * Adiag2 = (double*)malloc(N*sizeof(double));
    C_DCOPY(N,Adiag,1,Adiag2,1);

    int * small2big = (int*)malloc(19*M*sizeof(int));
    memset((void*)small2big,'\0',19*M*sizeof(int));

    int smart_guess = 1;
    if ( hamiltonian_element == NULL ) {
        smart_guess = 0;
    }

    if(N > maxdim-M) init_dim = maxdim-M;
    else             init_dim =        M;

    if (smart_guess) { /* Use eigenvectors of a sub-matrix as initial guesses */

        //Adiag = init_array(N);
        //for (i = 0; i < N; i++) {
        //    Adiag[i] = A[i][i];
        //}

        for (i = 0; i < init_dim; i++) {
            minimum = Adiag2[0];
            min_pos = 0;
            for (j = 1; j < N; j++)
                if (Adiag2[j] < minimum) {
                    minimum = Adiag2[j];
                    min_pos = j;
                    small2big[i] = j;
                }

            Adiag2[min_pos] = BIGNUM;
            lambda_old_p[i] = minimum;
        }

        for (i = 0; i < init_dim; i++) {
            for (j = 0; j < init_dim; j++) G_p[i][j] = hamiltonian_element(small2big[i],small2big[j],data);
        }

        sq_rsp(init_dim, init_dim, G_p, lambda_p, 1, alpha_p, 1e-12);

        for (i = 0; i < init_dim; i++) {
            for (j = 0; j < init_dim; j++) {
                b_p[i][small2big[j]] = alpha_p[j][i];
            }
        }

    } else { /* Use unit vectors as initial guesses */
        //Adiag = init_array(N);
        //for (i = 0; i < N; i++) {
        //    Adiag[i] = A[i][i];
        //}
        for (i = 0; i < init_dim; i++) {
            minimum = Adiag2[0];
            min_pos = 0;
            for (j = 1; j < N; j++){
                if (Adiag2[j] < minimum) {
                    minimum = Adiag2[j];
                    min_pos = j;
                }
            }

            b_p[i][min_pos] = 1.0;
            Adiag2[min_pos] = BIGNUM;
            lambda_old_p[i] = minimum;
        }

        //free(Adiag);
    }

    free(Adiag2);
    free(small2big);

    L = init_dim;
    iter = 0;
    converged = 0;
    //conv = init_int_array(M); /* boolean array for convergence of each
    //                             root */
    int * conv = (int*)malloc(M*sizeof(int));
    memset((void*)conv,'\0',M*sizeof(int));

    while (converged < M && iter < MAXIT) {
        skip_check = 0;

        // use call-back function to get sigma vector:
        build_sigma(N,L,sigma_p,b_p,data);
        C_DGEMM('n', 'n', L, L, N, 1.0, &(b_p[0][0]), N, &(sigma_p[0][0]), maxdim, 0.0, &(G_p[0][0]), maxdim);

        /* diagonalize mini-matrix */
        sq_rsp(L, L, G_p, lambda_p, 1, alpha_p, 1e-12);

        /* form preconditioned residue vectors */
        for (k = 0; k < M; k++) {
            for (I = 0; I < N; I++) {
                f_p[k][I] = 0.0;
                for (i = 0; i < L; i++) {
                    f_p[k][I] += alpha_p[i][k] * (sigma_p[I][i] - lambda_p[k] * b_p[i][I]);
                }
                denom = lambda_p[k] - Adiag[I];//A[I][I];
                if (std::fabs(denom) > 1e-6) {
                    f_p[k][I] /= denom;
                }else {
                    f_p[k][I] = 0.0;
                }
            }
        }

        /* normalize each residual */
        for (k = 0; k < M; k++) {
            norm = 0.0;
            for (I = 0; I < N; I++) {
                norm += f_p[k][I] * f_p[k][I];
            }
            norm = sqrt(norm);
            for (I = 0; I < N; I++) {
                if (norm > 1e-6) {
                    f_p[k][I] /= norm;
                }else {
                    f_p[k][I] = 0.0;
                }
            }
        }

        /* schmidt orthogonalize the f[k] against the set of b[i] and add
           new vectors */
        for (k = 0, numf = 0; k < M; k++) {
            if (schmidt_add(b_p, L, N, f_p[k])) {
                L++;
                numf++;
            }
        }

        /* If L is close to maxdim, collapse to one guess per root */
        if (maxdim - L < M) {
            if (print) {
                printf("Subspace too large: maxdim = %d, L = %d\n", maxdim, L);
                printf("Collapsing eigenvectors.\n");
            }
            for (i = 0; i < M; i++) {
                memset((void *)bnew_p[i], 0, N * sizeof(double));
                for (j = 0; j < L; j++) {
                    for (k = 0; k < N; k++) {
                        bnew_p[i][k] += alpha_p[j][i] * b_p[j][k];
                    }
                }
            }

            /* copy new vectors into place */
            for (i = 0; i < M; i++)
                for (k = 0; k < N; k++) b_p[i][k] = bnew_p[i][k];

            skip_check = 1;

            L = M;
        }

        /* check convergence on all roots */
        if (!skip_check) {
            converged = 0;
            //zero_int_array(conv, M);
            memset((void*)conv,'\0',M*sizeof(int));
            if (print) {
                printf("Root      Eigenvalue       Delta  Converged?\n");
                printf("---- -------------------- ------- ----------\n");
            }
            for (k = 0; k < M; k++) {
                diff = std::fabs(lambda_p[k] - lambda_old_p[k]);
                if (diff < cutoff) {
                    conv[k] = 1;
                    converged++;
                }
                lambda_old_p[k] = lambda_p[k];
                if (print) {
                    printf("%3d  %20.14f %4.3e    %1s\n", k, lambda_p[k], diff, conv[k] == 1 ? "Y" : "N");
                }
            }
        }

        iter++;
    }

    /* generate final eigenvalues and eigenvectors */
    if (converged == M) {
        for (i = 0; i < M; i++) {
            eps[i] = lambda_p[i];
            for (j = 0; j < L; j++) {
                double dum = 0.0;
                for (I = 0; I < N; I++) {
                    v[i][I] += alpha_p[j][i] * b_p[j][I];
                }
            }
            // ensure ci vector is normalized?
            norm = 0.0;
            for(int I = 0; I < N; I++) {
                norm += v[i][I] * v[i][I];
            }
            norm = sqrt(norm);
            for(size_t I = 0; I < N; I++) {
               v[i][I] /= norm;
            }
            //printf("hey is this thing normalized? %5i %20.12lf\n",i,norm);fflush(stdout);
            if ( fabs(1.0 - norm) > 1e-6 ) {
                throw PsiException("davidson CI vector is not normalized. try increasing maxdim?",__FILE__,__LINE__);
            }
        }
        if (print) printf("Davidson algorithm converged in %zu iterations.\n", iter);
    }

    free(conv);
    //free_block(b);
    //free_block(bnew);
    //free_block(sigma);
    //free_block(G);
    //free_block(f);
    //free_block(alpha);
    //free(lambda);
    //free(lambda_old);

    return converged;
}

}
