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

#include "nonsym_davidson_solver.h"
#include "blas.h"

using namespace std;
using namespace psi;

#define BIGNUM 1E100
#define SMALLNUM -1E100
#define MAXIT 1000

namespace psi{
void GeneralizedEigenvalueProblem(long int N, double **A, double **B, double **vr, double *eig, bool * energy_is_real){

    char JOBVR = 'V';
    char JOBVL = 'V';
    std::shared_ptr<Matrix> As (new Matrix(N,N));
    std::shared_ptr<Matrix> Bs (new Matrix(N,N));
    double** Asp = As->pointer();
    double** Bsp = Bs->pointer();
    for(int i = 0; i < N; i++) {
       for(int j = 0; j < N; j++) {
          Asp[i][j] = A[i][j];
          Bsp[i][j] = B[i][j];
       }
    }

    long int LDA  = N;
    long int LDB  = N;
    long int LDVR = N;
    long int LDVL = N;
    long int LWORK = 8*N;
    std::shared_ptr<Matrix> VRp (new Matrix(N,N));
    std::shared_ptr<Matrix> VLp (new Matrix(N,N));
    std::shared_ptr<Vector> ALPHARp (new Vector(N));
    std::shared_ptr<Vector> ALPHAIp (new Vector(N));
    std::shared_ptr<Vector> BETAp   (new Vector(N));
    //std::shared_ptr<Vector> WORK   (new Vector(N));
    double* ALPHAR = ALPHARp->pointer();
    double* ALPHAI = ALPHAIp->pointer();
    double* BETA = BETAp->pointer();
    double** VR = VRp->pointer();
    double** VL = VLp->pointer();

    double*WORK=(double*)malloc(LWORK*sizeof(double));
    memset((void*)WORK,'\0',LWORK*sizeof(double));
    integer INFO=0;
    //INFO = C_DGGEV(JOBVL,JOBVR,N,Asp[0],LDA,Bsp[0],LDB,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK);
    INFO = C_DGGEV(JOBVL,JOBVR,N,Asp[0],LDA,Bsp[0],LDB,ALPHAR,ALPHAI,BETA,VL[0],LDVL,VR[0],LDVR,WORK,LWORK);
    //DGGEV(JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO);
    outfile->Printf("\n");
    //outfile->Printf("    ==> positive excitation energies <==\n");
    outfile->Printf("    ==> excitation energies <==\n");
    outfile->Printf("\n");

    outfile->Printf("    ");
    outfile->Printf("             Re(eig)");
    outfile->Printf("             Im(eig)");
    outfile->Printf("\n");
    outfile->Printf("\n");
    for (long int i = 0; i < N; i++) {
        outfile->Printf("    %20.12lf %20.12lf\n", ALPHAR[i]/BETA[i],/*27.2114*/  ALPHAI[i]/BETA[i]);
       if ( fabs(ALPHAI[i]/BETA[i]) > 1e-6 ) {
           //eig[i] = 0.0;
           eig[i] = ALPHAR[i] / BETA[i];
           outfile->Printf("<<<warning>>> excitation energy is complex: %20.12lf + %20.12lf I\n",eig[i], ALPHAI[i]/BETA[i]);
           energy_is_real[i] = false;
       }else {
           eig[i] = ALPHAR[i] / BETA[i];
           energy_is_real[i] = true;
       }
        //if ( eig[i] > 0.0 && fabs(ALPHAI[i]/BETA[i]) < 1e-6) {
        //    outfile->Printf("%20.12lf %20.12lf %20.12lf\n",eig[i], eig[i] * 27.21138, ALPHAI[i]/BETA[i]);
        //}else {
        //    outfile->Printf("%20.12lf %20.12lf %20.12lf (negative)\n",eig[i], eig[i] * 27.21138, ALPHAI[i]/BETA[i]);
        //}
    }
    outfile->Printf("\n");
    outfile->Printf("info: %5i\n",INFO);
  //stable selection sort
    std::shared_ptr<Vector> evr_temp (new Vector(N));
    double* evrp = evr_temp->pointer();
    size_t i, j, jmin;

    for (j = 0; j < N-1; j++) {
        jmin = j;
        // find index of minimum eigenvalue from j to n
        for (i = j+1; i < N; i++) {
            if (eig[jmin] < eig[i] ) {
                jmin = i;
            }
        }
        //save current eigenvalue and eigenvectors associated with index jmin 
        double eigr = eig[jmin];
        for (i = 0; i < N; i++) {
            evrp[i] = VR[jmin][i];
        }
        //shift values
        while (jmin > j) {
            eig[jmin] = eig[jmin-1];
            for (i = 0; i < N; i++) {
                VR[jmin][i] = VR[jmin-1][i];
            }
            jmin--;
        }
        eig[j] = eigr;
        for (i = 0; i < N; i++) {
            VR[j][i] = evrp[i];
        }
     }


    //C_DCOPY(N*N,VR,1,Bsp[0],1);
    //C_DCOPY(N*N,VL,1,Asp[0],1);
    for(int i = 0; i<N;i++) {
       for(int j = 0; j<N;j++) {
       vr[i][j] = VR[i][j];
       }
    }
    //free(VR);
    //free(VL);
    free(WORK);
    //free(ALPHAR);
    //free(ALPHAI);
    //free(BETA);
}
void NonSymmetricEigenvalueEigenvector(long int dim, double **M, double * eigval, double * wi, double ** el, double ** er) {
    long int info;
    char vl = 'V';
    char vr = 'V';
    long int n = dim; 
    long int lwork = 4*n; 
    double * work  = (doublereal*)malloc(lwork*sizeof(doublereal));
    //double * wi = (double*)malloc(n*sizeof(double));
    //copy the non-zero part of the matrix we need to solve for eigenvalue to another mini matrix
    std::shared_ptr<Matrix> As (new Matrix(dim,dim));
    std::shared_ptr<Matrix> ers (new Matrix(dim,dim));
    std::shared_ptr<Matrix> els (new Matrix(dim,dim));
    std::shared_ptr<Vector> eigr (new Vector(dim));
    std::shared_ptr<Vector> eigi (new Vector(dim));
    double** Asp = As->pointer();
    double** elsp = els->pointer();
    double** ersp = ers->pointer();
    double* eigrp = eigr->pointer();
    double* eigip = eigi->pointer();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
             Asp[j][i]=M[i][j];
        }
    }   
            
    psi::fnocc::DGEEV(vl,vr,n,Asp[0],n,eigrp,eigip,elsp[0],n,ersp[0],n,work,lwork,info);

    // perform argument sort by pairing eigenvalues with index
    typedef std::pair<std::complex<double>, size_t> complex_w_idx;
    auto *eig_pair = (complex_w_idx *) malloc(n * sizeof(complex_w_idx));

    for (size_t i = 0; i < n; i++) {
        eig_pair[i].first = {eigrp[i], eigip[i]};
        eig_pair[i].second = i;
    }

    // sort by eigenvalue
    std::stable_sort(eig_pair, eig_pair + n,
        [](const complex_w_idx &l, const complex_w_idx &r) {
            // sort by real part (lowest first)
            std::complex<double> l_val = l.first;
            std::complex<double> r_val = r.first;
            bool correctOrder = l_val.real() < r_val.real();

            if (!correctOrder) {
                // sort by imaginary part (positive first) if same real part
                bool same_real_part = fabs(l_val.real() - r_val.real()) <= 1e-16;
                if (same_real_part) {
                    return l_val.imag() > r_val.imag();
                }
            }

            // return true if l should come before r
            return correctOrder;
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
    }

     free(work);

}

void qr_orthogonalization(long int L, long int M, long int N, double **A) {
   //this will orthonormalize N first rows of a matrix A with dimension L x M (N <= L <= M).
   //copy nonzero part of A to a smaller maxtrix before orthogonalization
   std::shared_ptr<Matrix> At (new Matrix(N,M));
   std::shared_ptr<Vector> tau (new Vector(N));
   double** Atp = At->pointer();
   //provide variable L here just in case I need to change type of A later
   for (int i = 0; i < N; i++) {
       for (int j = 0; j < M; j++) {
           Atp[i][j] = A[i][j];	 
       }
   }
   //outfile->Printf("At before orthogonalization");
   // At->print();
   //double * tau= (double*)malloc(N*sizeof(double));
   double * taup= tau->pointer();
   long int lwork = 8*N;
   double * work1  = (doublereal*)malloc(lwork*sizeof(doublereal));
   double * work2  = (doublereal*)malloc(lwork*sizeof(doublereal));
   long int info1;
   long int info2;
   psi::fnocc::DGEQRF( M, N, Atp[0], M, taup, work1, lwork, info1);
   //tau->print();
   //At->print();
   outfile->Printf("info1 %d\n",info1);
   psi::fnocc::DORGQR(M,N,N,Atp[0],M,taup,work2,lwork,info2);
   outfile->Printf("info2 %d\n",info2);
   //copy orthonormal rows back to A
   for (int i = 0; i < N; i++) {
       for (int j = 0; j < M; j++) {
           A[i][j] = Atp[i][j];	 
       }
   }
   //outfile->Printf("At after orthogonalization");
   //At->print();
   free(work1);
   free(work2);
}


Nonsym_DavidsonSolver::Nonsym_DavidsonSolver(){
}
Nonsym_DavidsonSolver::~Nonsym_DavidsonSolver(){
}

void Nonsym_DavidsonSolver::solve(double *Adiag, int N, int M, double *reval, double **rer, double **rel, BuildSigma build_sigma, int maxdim, size_t init_dim, double residual_norm, bool use_residual_norm) {
    int min_pos;
    double minimum;
    //int N = o*v;
    //int M =10;//number of roots
    //int maxdim = 18 * M;//max size of the subspace
    // init_dim= N;
 
    std::shared_ptr<Matrix> B ( new Matrix(maxdim,maxdim) );//Q^T(A-B)(A+B)Q
    std::shared_ptr<Matrix> Q ( new Matrix(maxdim,N) );
    std::shared_ptr<Matrix> Qn ( new Matrix(maxdim,N) );
    std::shared_ptr<Matrix> sigmar ( new Matrix(N,maxdim) );//(A-B)(A+B)q
    std::shared_ptr<Matrix> sigmal ( new Matrix(N,maxdim) );//(A-B)(A+B)q
    std::shared_ptr<Matrix> Cl ( new Matrix(maxdim,maxdim) );
    std::shared_ptr<Matrix> Cr ( new Matrix(maxdim,maxdim) );
    std::shared_ptr<Matrix> Rl ( new Matrix(M+1,N) );
    std::shared_ptr<Matrix> Rr ( new Matrix(M+1,N) );
    std::shared_ptr<Vector> lambda ( new Vector(maxdim) );
    std::shared_ptr<Vector> lambdai ( new Vector(maxdim) );
    std::shared_ptr<Vector> lambda_old ( new Vector(M+1) );
    std::shared_ptr<Vector> lambdai_old ( new Vector(M+1) );
    
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
    double * Adiag2 = (double*)malloc(N*sizeof(double));
    C_DCOPY(N,Adiag,1,Adiag2,1);
    //use unit vectors as guess
    for (int i = 0; i < init_dim; i++) {
        minimum = Adiag2[0];
        min_pos = 0;
        for (int j = 1; j < N; j++){
            if (Adiag2[j] < minimum) {
                minimum = Adiag2[j];
                min_pos = j;
            }
        }
        Qp[i][min_pos] = 1.0;
        Adiag2[min_pos] = BIGNUM;
    }
    free(Adiag2);
    int L = init_dim;

    int iter = 0;
    int maxit = eom_maxiter;
    bool convergence = false;
    outfile->Printf("\n");
    while (iter < maxit) {
      //this code only checks convergence of the real part of the eigenvalues
      for (int i = 0; i < M+1; i++) {
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
      build_sigma(N,maxdim,L,Qp,sigmarp,sigmalp);
      
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
      NonSymmetricEigenvalueEigenvector(L, Bp, lambdap, lambdaip, Clp,Crp);

     for (size_t k = 0; k < M+1; k++) {
          for (size_t n = 0; n < N; n++) {
              Rrp[k][n]= 0.0;
              Rlp[k][n]= 0.0;
          }
      }
      //residual
      size_t complex_root = 0;
      for (size_t k = 0; k < M; k++) {
          if (fabs(lambdaip[k]) < 1e-16) {
             for (size_t n = 0; n < N; n++) {
                 double dum1 = 0.0;
                 double dum2 = 0.0;
                 for (size_t m = 0; m < L; m++) {
                    dum1 += (sigmarp[n][m] - lambdap[k]  * Qp[m][n]) * Crp[k][m];
                    dum2 += (sigmalp[n][m] - lambdap[k]  * Qp[m][n]) * Clp[k][m];

                 }
                   Rrp[k][n]= dum1;
                   Rlp[k][n]= dum2;
             }
          }
          else {
             complex_root++;
             if (complex_root%2==0) continue;
	     double re = lambdap[k];
	     double im = lambdaip[k];
             for (size_t n = 0; n < N; n++) {
                 double res_rr = 0.0;
                 double res_ri = 0.0;
                 double res_lr = 0.0;
                 double res_li = 0.0;
                 for (size_t m = 0; m < L; m++) {
		     // f_R = AQR-wQR
                     res_rr += sigmarp[n][m] * Crp[k  ][m] - (re * Crp[k][m] - im * Crp[k+1][m]) * Qp[m][n];
                     res_ri += sigmarp[n][m] * Crp[k+1][m] - (im * Crp[k][m] + re * Crp[k+1][m]) * Qp[m][n];
		     // f_L = A^T QL - (w*) QL
                     res_lr += sigmalp[n][m] * Clp[k  ][m] - (re * Clp[k][m] + im * Clp[k+1][m]) * Qp[m][n];
                     res_li += sigmalp[n][m] * Clp[k+1][m] - (-im * Clp[k][m] + re * Clp[k+1][m]) * Qp[m][n];
                 }
                   Rrp[k  ][n] = res_rr;
                   Rrp[k+1][n] = res_ri;
                   Rlp[k  ][n] = res_lr;
                   Rlp[k+1][n] = res_li;
             }
          }
      }
       
      size_t Mn;
      if (complex_root%2!=0) {
         Mn =M+1;
      }
      else Mn=M;
      vector<int> index;
      int check_residual = 0;
      int conv = 0;
      outfile->Printf("residual norm, dimension %d\n", L);
      for (int k = 0; k < Mn; k++) {
          double normvalr = C_DDOT(N, Rrp[k], 1, Rrp[k], 1);
          double normvall = C_DDOT(N, Rlp[k], 1, Rlp[k], 1);
          normvalr = sqrt(normvalr);
          normvall = sqrt(normvall);

          double e_error = lambdap[k]-lambda_oldp[k];
          bool e_conv = fabs(e_error) < eom_r_conv;

          // residual norms are less than the threshold
          bool r_resid_conv = normvalr < residual_norm;
          bool l_resid_conv = normvall < residual_norm;

          // or residual norms are less than convergence criterion
          r_resid_conv |= normvall < eom_r_conv;
          l_resid_conv |= normvalr < eom_r_conv;

          if (r_resid_conv && l_resid_conv && e_conv) {
                 conv++; // count the number of converged roots
          }
          else {
             //find rows of residual that have norm larger than some threshold
             index.push_back(k);
          }
          outfile->Printf("residual norm%d %20.12lf %20.12lf\n", k,normvalr,normvall);
      }

      outfile->Printf("check_convergence %d\n", conv);
      outfile->Printf("energy error, dimension %d\n", L);
      for (int k = 0; k < Mn; k++) {
          double error = lambdap[k]-lambda_oldp[k];
           outfile->Printf("eig_error%d %20.12lf %20.12lf \n", k,error,lambdap[k]);
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
                 double dum = lambdap[k]-Adiag[n];
                 if (fabs(dum) > 1e-6) {
                    Rrp[k][n] /= dum;
                    Rlp[k][n] /= dum;
                 }else {
                    Rrp[k][n]= 0.0;
                    Rlp[k][n]= 0.0;
                 }
             }
             //outfile->Printf("%20.12lf %20.12lf\n", lambdap[k], lambdaip[k]);
          }
          else {
             complex_root++;
             if (complex_root%2==0) continue;
             double re = lambdap[k];
             double im = lambdaip[k];
             //outfile->Printf("%20.12lf %20.12lf\n", re,im);
             for (size_t n = 0; n < N; n++) {
                 double md = (re-Adiag[n]) * (re-Adiag[n]) + im * im;
                 double res_rr = ((re-Adiag[n]) * Rrp[k][n] + im            * Rrp[k+1][n])/md;
                 double res_ri = (-im           * Rrp[k][n] + (re-Adiag[n]) * Rrp[k+1][n])/md;
                 double res_lr = ((re-Adiag[n]) * Rlp[k][n] - im            * Rlp[k+1][n])/md;
                 double res_li = ( im           * Rlp[k][n] + (re-Adiag[n]) * Rlp[k+1][n])/md;
                 Rrp[k  ][n] = res_rr;
                 Rrp[k+1][n] = res_ri;
                 Rlp[k  ][n] = res_lr;
                 Rlp[k+1][n] = res_li;
             }
          }
      }
      if (use_residual_norm==true) {
          outfile->Printf("\n Use residual norm to control the numer of correction vectors added to new subspace\n");
          outfile->Printf("dimension %d, number of complex roots %d\n",L,complex_root);
          if (complex_root%2!=0) {
             Mn = M+1;
             outfile->Printf("need %d vectors (for imaginary part)\n", Mn);
          }
          else Mn = M;
          if (maxdim - L < 2 *index.size()) {
              outfile->Printf("Subspace too large: maxdim = %d, L = %d\n", maxdim, L);
              printf("Collapsing eigenvectors.\n");
              for (size_t i = 0; i < maxdim; i++) {
                  for (size_t k = 0; k < N; k++) {
                      Qnp[i][k] = 0.0;
                  }
              }
              for (size_t i = 0; i < Mn; i++) {
                  for (size_t k = 0; k < N; k++) {
                      for (size_t j = 0; j < L; j++) {
                          Qnp[2*i][k] += Clp[i][j] * Qp[j][k];
                          Qnp[2*i+1][k] += Crp[i][j] * Qp[j][k];
                      }
                  }
              }
              for (size_t i = 0; i < index.size(); i++) {
                  for (size_t k = 0; k < N; k++) {
                      Qnp[2*i+2*Mn][k] = Rlp[index[i]][k];
                      Qnp[2*i+1+2*Mn][k] = Rrp[index[i]][k];
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
              for (size_t i = 0; i < 2 * (Mn +index.size()) ; i++) {
                  for (size_t k = 0; k < N; k++) {
                      Qp[i][k] = Qnp[i][k];
                  }
              }

              L = 2 * (Mn +index.size());
              qr_orthogonalization(maxdim,N,L,Qp);
          }else {
              for (size_t i = 0; i < maxdim; i++) {
                  for (size_t k = 0; k < N; k++) {
                      Qnp[i][k] = 0.0;
                  }
              }
              for (size_t k = 0; k < index.size(); k++) {
                  for (size_t n = 0; n < N; n++) {
                      Qnp[2*k][n]= Rlp[index[k]][n];
                      Qnp[2*k+1][n]= Rrp[index[k]][n];
                  }
              }
	     //Lapack Householder
              for (size_t k = 0; k < index.size(); k++) {
                  for (size_t n = 0; n < N; n++) {
                      Qp[L+k][n]= Rlp[index[k]][n];
                      Qp[L+index.size()+k][n]= Rrp[index[k]][n];
                  }
              }
              L=L + 2*index.size();
              qr_orthogonalization(maxdim,N,L,Qp);
          }
      }else {
          if (complex_root%2!=0) {
             Mn =M+1;
             outfile->Printf("need %d vectors (for imaginary part)\n", Mn);
          }
          else Mn=M;

          if (maxdim - L < 2 * Mn) {
             outfile->Printf("Subspace too large: maxdim = %d, L = %d\n", maxdim, L);
             printf("Collapsing eigenvectors.\n");
             //zero current subspace bases?
             for (size_t i = 0; i < maxdim; i++) {
                 for (size_t k = 0; k < N; k++) {
                     Qnp[i][k] = 0.0;
                 }
             }

             for (size_t i = 0; i < Mn; i++) {
                 for (size_t k = 0; k < N; k++) {
                     for (size_t j = 0; j < L; j++) {
                         Qnp[2*i][k] += Clp[i][j] * Qp[j][k];
                         Qnp[2*i+1][k] += Crp[i][j] * Qp[j][k];
                     }
                         Qnp[2*i+2*Mn][k] = Rlp[i][k];
                         Qnp[2*i+1+2*Mn][k] = Rrp[i][k];
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
             qr_orthogonalization(maxdim,N,L,Qp);
          }else {
             for (size_t i = 0; i < maxdim; i++) {
                 for (size_t k = 0; k < N; k++) {
                     Qnp[i][k] = 0.0;
                 }
             }
             for (size_t k = 0; k < Mn; k++) {
                 for (size_t n = 0; n < N; n++) {
                     Qnp[2*k][n]= Rlp[k][n];
                     Qnp[2*k+1][n]= Rrp[k][n];
                 }
             }
	         //Lapack Householder
             for (size_t k = 0; k < Mn; k++) {
                 for (size_t n = 0; n < N; n++) {
                     Qp[L+k][n]= Rrp[k][n];
                     Qp[L+Mn+k][n]= Rlp[k][n];
                 }
             }
             L=L + 2*Mn;
             qr_orthogonalization(maxdim,N,L,Qp);
             //Q->print();
          }
      }

      //get the eigenvectors from desired roots
      outfile->Printf("\nIter %5d / %-5d\n", ++iter, maxit);
      writeVectors(N, M, maxdim, reval, rer, rel, Clp, Crp, Qp, lambdap, lambdaip, L, false);

    }
    outfile->Printf("\n");
    if (convergence == false) {
        outfile->Printf("Davidson procedure fails after %d iterations, try increase maxdim", maxit);
        exit(0);
    }
    outfile->Printf("\n");

}

int Nonsym_DavidsonSolver::writeVectors(int N, int M, int maxdim, double *reval, double **rer, double **rel,
                                        double *const *Clp, double *const *Crp, double *const *Qp,
                                        const double *lambdap, const double *lambdaip, int L, bool doPrint) const {
    int M_i = M;
    int num_i = 0;
    for (int k = 0; k < M && M <= maxdim; ++k) {
        bool isComplex = fabs(lambdaip[k]) >= 1e-16;
        if (isComplex) {
            // move to the next eigenvalue for complex conjugate pairs
            ++M;
            ++k;
            ++num_i;
            isComplex = true;
        }
        int l = k - num_i;
        if (isComplex && doPrint) outfile->Printf("%5d %20.12lf %20.12lf\n", l, lambdap[k - 1], lambdaip[k - 1]);
        if (doPrint) outfile->Printf("%5d %20.12lf --------------------\n", l, lambdap[k], lambdaip[k]);
        if (!isComplex) {
            reval[l] = lambdap[k]; // eigenvalues
            for (int n = 0; n < N; n++) {
                rer[l][n] = 0.0;
                rel[l][n] = 0.0;

                // eigenvectors
                for (int m = 0; m < L; m++) {
                    rer[l][n] += Qp[m][n] * Crp[k][m];
                    rel[l][n] += Qp[m][n] * Clp[k][m];
                }
            }
        } else {
            reval[l] = lambdap[k]; // get real part of eigenvalues
            for (int n = 0; n < N; n++) {
                rer[l][n] = 0.0;
                rel[l][n] = 0.0;
                std::complex<double> evec_r;
                std::complex<double> evec_l;

                // get complex eigenvectors
                for (int m = 0; m < L; m++) {
                    std::complex dumr(Qp[m][n] * Crp[k - 1][m] * lambdap[k - 1],
                                      Qp[m][n] * Crp[k - 1][m] * lambdaip[k - 1]);
                    std::complex duml(Qp[m][n] * Clp[k][m] * lambdap[k], Qp[m][n] * Clp[k][m] * lambdaip[k]);

                    evec_r += dumr;
                    evec_l += duml;
                }

                // get real components of eigenvectors
                rer[l][n] = evec_r.real();
                rel[l][n] = evec_l.real();
            }
        }
    }
    if (M > maxdim){
        outfile->Printf("WARNING: %d complex eigenvalues encountered. Subspace is insufficient to store them and provide the requested roots.\n", num_i);
    }
    M = M_i;
    return M;
}

void Nonsym_DavidsonSolver::real_generalized_eigenvalue_problem(double *Hdiag, double *Sdiag, size_t N, size_t M, double *reval, double **rer, BuildSigma2 build_sigma2, size_t maxdim,size_t init_dim, double residual_norm, bool use_residual_norm) {
    double minimum;
    int min_pos;
    std::shared_ptr<Matrix> mh ( new Matrix(maxdim,maxdim) );//Q^T*H* Q
    std::shared_ptr<Matrix> ms ( new Matrix(maxdim,maxdim) );//Q^T*S*Q
    std::shared_ptr<Matrix> Q ( new Matrix(maxdim,N) );
    std::shared_ptr<Matrix> Qn ( new Matrix(maxdim,N) );
    std::shared_ptr<Matrix> vh ( new Matrix(N,maxdim) );//Hq
    std::shared_ptr<Matrix> vs ( new Matrix(N,maxdim) );//Sq
    std::shared_ptr<Matrix> Cr ( new Matrix(maxdim,maxdim) );//eigenvector in subspace
    std::shared_ptr<Matrix> Rr ( new Matrix(M,N) );
    //std::shared_ptr<Matrix> Rr ( new Matrix(2*M,N) );
    std::shared_ptr<Matrix> QC ( new Matrix(M,N) );//Q*C
    std::shared_ptr<Matrix> pQC ( new Matrix(M,N) );//apply precondtioner on Q*C
    std::shared_ptr<Vector> epsilon ( new Vector(M) );
    std::shared_ptr<Vector> lambda ( new Vector(maxdim) );
    std::shared_ptr<Vector> lambda_old ( new Vector(M) );

    double **mhp = mh->pointer();
    double **msp = ms->pointer();
    double **vhp = vh->pointer();
    double **vsp = vs->pointer();
    double **Qp = Q->pointer();
    double **Qnp = Qn->pointer();
    double **Crp = Cr->pointer();
    double **Rrp = Rr->pointer();
    double **QCp = QC->pointer();
    double **pQCp = pQC->pointer();
    double *epsilonp = epsilon->pointer();
    double *lambdap = lambda->pointer();
    double *lambda_oldp = lambda_old->pointer();

    //maybe maxdim <=N/2
    double * Adiag2 = (double*)malloc(N*sizeof(double));
    for (size_t m = 0; m < N; m++) {
           // Adiag[m] = pabp[m][m];
            Adiag2[m] = Sdiag[m]/Hdiag[m];
    }
    //use unit vectors as guess
    for (size_t i = 0; i < init_dim; i++) {
           minimum = Adiag2[0];
           min_pos = 0;
           for (size_t j = 1; j < N; j++){
               if (Adiag2[j] > minimum) {
                   minimum = Adiag2[j];
                   min_pos = j;
               }
           }

           Qp[i][min_pos] = 1.0;
           Adiag2[min_pos] = SMALLNUM;
           //Adiag2[min_pos] = BIGNUM;
          // lambda_oldp[i] = minimum;
    }
    size_t L = init_dim;
    //for (int i = 0; i < L; i++) {
    //    for (int m = 0; m < N/2; m++) {
    //        double dum =Qp[i][m];
    //       // Qp[i][m] = 0;
    //        Qp[i+L][m+N/2]=dum;
    //    }
    //}
    //L=2*init_dim;
    //srand(time(0));
    //for (size_t i = 0; i < init_dim; i++) {
    //    for (size_t j = 1; j < N; j++){
    //   Qp[i][j] = ((double)rand()/RAND_MAX);
    //    }
    //}
    //qr_orthogonalization(maxdim,N,L,Qp);

    //free(Adiag);
    free(Adiag2);
    int iter = 0;
    int maxit = 1000;
    bool convergence = false;
    outfile->Printf("\n");
    while (iter < maxit) {
// outfile->Printf("test orthogonality of basis Q\n");
//    for (int m = 0; m < L; m++) {
//        for (int n = 0; n < L; n++) {
//        double a = C_DDOT(N, Qp[m], 1, Qp[n], 1);
//          outfile->Printf("%d %d %20.12lf\n",m, n, a);
//        }
//    }
    for (size_t i = 0; i < M; i++) {
        lambda_oldp[i] = lambdap[i];
    }

    //make all elements of eigenvectors and eigenvalues zero
    for (size_t i = 0; i < maxdim; i++) {
         lambdap[i] = 0.0;
         for (size_t j = 0; j < maxdim; j++) {
             Crp[i][j] = 0.0;
         }
    }
    //make all elements of sigma zero
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < maxdim; j++) {
            vhp[i][j] = 0.0;
            vsp[i][j] = 0.0;

        }
    }

    //make all elements of H and S projection zero
    for (size_t i = 0; i < maxdim; i++) {
        for (size_t j = 0; j < maxdim; j++) {
            mhp[i][j] = 0.0;
            msp[i][j] = 0.0;

        }
    }
    //range of H and S. This should be done in some sigma build function
    build_sigma2(N,maxdim,L,Qp,vhp,vsp);
    //projection of H and S onto subspace
    for (size_t i = 0; i < L; i++) {
        for (size_t j = 0; j < L; j++) {
            double dum1 = 0.0;
            double dum2 = 0.0;
            for (int k = 0; k < N; k++) {
                dum1 += vhp[k][j] * Qp[i][k];
                dum2 += vsp[k][j] * Qp[i][k];
            }
            mhp[i][j] = dum1;
            msp[i][j] = dum2;
        }
    }
    //ms->print();
    bool * energy_is_real = (bool*)malloc(L*sizeof(bool));
    GeneralizedEigenvalueProblem(L,msp,mhp,Crp,lambdap,energy_is_real);
    //lambda->print();
    outfile->Printf("inverse eig, dimension %d\n",L);
    for (size_t k = 0; k < M; k++) {
        outfile->Printf("inverse eig%d %20.12lf\n",k,27.21138 * 1.0/lambdap[k]);
    }
    for (size_t k = 0; k < M; k++) {
        epsilonp[k]= 0.0;
        for (size_t n = 0; n < N; n++) {
            Rrp[k][n]= 0.0;
            QCp[k][n]= 0.0;
            pQCp[k][n]= 0.0;
        }
    }
    size_t n_positive = 0;
    for (size_t i = 0; i < L; i++) {
        if (lambdap[i]<1e-16) continue;
        n_positive++;
        //outfile->Printf("%20.12lf\n",lambdap[i]);
        for (size_t n = 0; n < N; n++) {
            for (size_t m = 0; m < L; m++) {
                 Rrp[n_positive-1][n] += vsp[n][m] * Crp[i][m] - lambdap[i] * vhp[n][m] *Crp[i][m];
            }
        }
        if (n_positive ==M) break;
    }

    //Rr->print();
    size_t Mn = M;
    vector<int> index;
    int check_residual = 0;
    int conv = 0;
    outfile->Printf("residual norm, dimension %d\n", L);
    for (int k = 0; k < Mn; k++) {
        double normvalr = C_DDOT(N, Rrp[k], 1, Rrp[k], 1);
        double error = lambdap[k]-lambda_oldp[k];
        normvalr = sqrt(normvalr);
        if ((normvalr < residual_norm) &&((fabs(error) < 1e-8)) ) {
           conv++;
        }
        else { 
        //find rows of residual that have norm larger than some threshold
           index.push_back(k);
        }
        outfile->Printf("residual norm%d %20.12lf \n", k,normvalr);
    }
    //outfile->Printf("size of index %d\n", index.size());
    if (index.size() >0) {
       for (auto it = index.begin(); it != index.end(); ++it) {
            // outfile->Printf("%d\n", *it);
       }
    } 
    outfile->Printf("check_convergence %d\n", conv);
    outfile->Printf("eig_error, dimension %d\n", L);
    for (int k = 0; k < Mn; k++) {
        double error = lambdap[k]-lambda_oldp[k];
         outfile->Printf("eig_error%d %30.16lf\n ", k,error);
         //if(fabs(error)<1e-8) outfile->Printf("check");
         //outfile->Printf("\n");
    }
    
    n_positive = 0;
    for (size_t i = 0; i < L; i++) {
        if (lambdap[i]<1e-16) continue;
        n_positive++;
        //outfile->Printf("%20.12lf\n",lambdap[i]);
        for (size_t n = 0; n < N; n++) {
            for (size_t m = 0; m < L; m++) {
                 QCp[n_positive-1][n] += Qp[m][n] * Crp[i][m];
            }
        }
        if (n_positive ==M) break;
    }
    if (conv == Mn) {
       convergence = true;
       outfile->Printf("energy and eigenvectors converged\n");
       break;
    }
    n_positive = 0;
    for (size_t i = 0; i < L; i++) {
        if (lambdap[i]<1e-16) continue;
        n_positive++;
        for (size_t n = 0; n < N; n++) {
            double dum = Sdiag[n] - lambdap[i] * Hdiag[n];
        //outfile->Printf("dum%20.12lf\n",dum);
            if (fabs(dum) >1e-10) {
               Rrp[n_positive-1][n] /= dum;
               pQCp[n_positive-1][n] = QCp[n_positive-1][n]/dum;
            }else {
        //outfile->Printf("%d %d dum%20.12lf\n",i,n,dum);
               Rrp[n_positive-1][n] = 0.0;
               pQCp[n_positive-1][n] = 0.0;
            }

        }
        if (n_positive ==M) break;
    }

    for (int k = 0; k < M; k++) {
        double a1 = C_DDOT(N, QCp[k], 1, Rrp[k], 1);
        double a2 = C_DDOT(N, QCp[k], 1, pQCp[k], 1);
        //outfile->Printf("a1%20.12lf,a2%20.12lf\n",a1,a2);
        if (fabs(a2) >1e-10) {
           epsilonp[k] = a1/a2;
        }else {
           epsilonp[k] = 0.0;
        }
        //outfile->Printf("epsilon%20.12lf\n",epsilonp[k]);
    }
    for (int k = 0; k < M; k++) {
        for (int n = 0; n < N; n++) {
            double dum = epsilonp[k] * pQCp[k][n] - Rrp[k][n];
            Rrp[k][n] = dum;
        }
    }
    Mn =M;
    if (use_residual_norm==true) {
          outfile->Printf("\n Use residual norm to control the numer of correction vectors added to new subspace\n");
          if (maxdim - L < index.size()) {
              outfile->Printf("Subspace too large: maxdim = %d, L = %d\n", maxdim, L);
              printf("Collapsing eigenvectors.\n");
              for (size_t i = 0; i < maxdim; i++) {
                  for (size_t k = 0; k < N; k++) {
                      Qnp[i][k] = 0.0;
                  }
              }
	      //add current eigenvectors
              for (size_t i = 0; i < Mn; i++) {
                  for (size_t k = 0; k < N; k++) {
		      Qnp[i][k] = QCp[i][k];
                  }
              }
	      //add current residuals 
              for (size_t i = 0; i < index.size(); i++) {
                  for (size_t k = 0; k < N; k++) {
                      Qnp[i+Mn][k] = Rrp[index[i]][k];
                  }
              }
              //zero current subspace bases?
              for (size_t i = 0; i < maxdim; i++) {
                  for (size_t k = 0; k < N; k++) {
                      Qp[i][k] = 0.0;
                  }
              }
	      /*
	      //Gram-Schmidt 
              double dotval, normval;
              size_t k,i, I;
              L = 0;
              for (k = 0; k < Mn +index.size(); k++) {
                  if (L>0) {
                     // Q->print();
                      for (i = 0; i < L; i++) {
                          dotval = C_DDOT(N, Qp[i], 1, Qnp[k], 1);
                          for (I = 0; I < N; I++) Qnp[k][I] -= dotval * Qp[i][I];
                      }
		      //reorthogonalization
                      for (i = 0; i < L; i++) {
                          dotval = C_DDOT(N, Qp[i], 1, Qnp[k], 1);
                          for (I = 0; I < N; I++) Qnp[k][I] -= dotval * Qp[i][I];
                      }
                  }
                  normval = C_DDOT(N, Qnp[k], 1, Qnp[k], 1);
                  normval = sqrt(normval);
                  //outfile->Printf("trial vector norm%30.18lf\n",normval);
                  if (normval > 1e-20) {
                      for (I = 0; I < N; I++) {
                          Qp[L][I] = Qnp[k][I] / normval;
                      }
                      L++;
                      //outfile->Printf("check orthogonality1\n");
                      for (int i = 0; i < L; i++) {
                          for (int j = 0; j < L; j++) {
                              double a = C_DDOT(N, Qp[i], 1, Qp[j], 1);
                              //outfile->Printf("%d %d %30.16lf",i,j,a);
                              if ((i!=j) && (fabs(a)>1e-12)) outfile->Printf(" detect linear dependency\n");
                              //outfile->Printf("\n");
                          }
                      }
                  }
              }
	      */
	      //Lapack Householder
              //copy new vectors into place 
              for (size_t i = 0; i < Mn +index.size(); i++) {
                  for (size_t k = 0; k < N; k++) {
                      Qp[i][k] = Qnp[i][k];
                  }
              }

              L = Mn +index.size();
              qr_orthogonalization(maxdim,N,L,Qp);
          }else {
              for (size_t i = 0; i < maxdim; i++) {
                  for (size_t k = 0; k < N; k++) {
                      Qnp[i][k] = 0.0;
                  }
              }
              for (size_t k = 0; k < index.size(); k++) {
                  for (size_t n = 0; n < N; n++) {
                      Qnp[k][n]= Rrp[index[k]][n];
                  }
              }
	      /*
	      //Gram-Schmidt 
              double dotval, normval;
              size_t k,i, I;
              for (k = 0; k < index.size(); k++) {
                  for (i = 0; i < L; i++) {
                      dotval = C_DDOT(N, Qp[i], 1, Qnp[k], 1);
                      for (I = 0; I < N; I++) Qnp[k][I] -= dotval * Qp[i][I];
                  }
                  //reorthogonalization
                  for (i = 0; i < L; i++) {
                      dotval = C_DDOT(N, Qp[i], 1, Qnp[k], 1);
                      for (I = 0; I < N; I++) Qnp[k][I] -= dotval * Qp[i][I];
                  }
                  normval = C_DDOT(N, Qnp[k], 1, Qnp[k], 1);
                  normval = sqrt(normval);
                  //outfile->Printf("trial vector norm2%30.18lf\n",normval);
                  if (normval > 1e-20) {
                      for (I = 0; I < N; I++) {
                          Qp[L][I] = Qnp[k][I] / normval;
                      }
                      L++;
                      //outfile->Printf("check orthogonality2\n");
                      for (int i = 0; i < L; i++) {
                          for (int j = 0; j < L; j++) {
                              double a = C_DDOT(N, Qp[i], 1, Qp[j], 1);
                              //outfile->Printf("%d %d %30.16lf",i,j,a);
                              if (i!=j && fabs(a)>1e-12) outfile->Printf(" detect linear dependency\n");
                              //outfile->Printf("\n");
                          }
                      }
                  }
              }
	      */
	     //Lapack Householder
              for (size_t k = 0; k < index.size(); k++) {
                  for (size_t n = 0; n < N; n++) {
                      Qp[L+k][n]= Rrp[index[k]][n];
                  }
              }
              L = L + index.size();
              qr_orthogonalization(maxdim,N,L,Qp);
          }
      }else {

          if (maxdim - L < Mn) {
             outfile->Printf("Subspace too large: maxdim = %d, L = %d\n", maxdim, L);
             printf("Collapsing eigenvectors.\n");
             //zero current subspace bases?
             for (size_t i = 0; i < maxdim; i++) {
                 for (size_t k = 0; k < N; k++) {
                     Qnp[i][k] = 0.0;
                 }
             }

             for (size_t i = 0; i < Mn; i++) {
                 for (size_t k = 0; k < N; k++) {
		     Qnp[i][k] = QCp[i][k];
                     Qnp[i+Mn][k] = Rrp[i][k];
                 }
             }
             for (size_t i = 0; i < maxdim; i++) {
                 for (size_t k = 0; k < N; k++) {
                     Qp[i][k] = 0.0;
                 }
             }
             /*
	     //Gram-Schmidt 
             double dotval, normval;
             size_t k,i, I;
             L = 0;
             for (k = 0; k < 2 * Mn; k++) {
                 if (L>0) {
                    // Q->print();
                     for (i = 0; i < L; i++) {
                         dotval = C_DDOT(N, Qp[i], 1, Qnp[k], 1);
                         for (I = 0; I < N; I++) Qnp[k][I] -= dotval * Qp[i][I];
                     }
                     //reorthogonalization
                     for (i = 0; i < L; i++) {
                         dotval = C_DDOT(N, Qp[i], 1, Qnp[k], 1);
                         for (I = 0; I < N; I++) Qnp[k][I] -= dotval * Qp[i][I];
                     }
                 }
                 normval = C_DDOT(N, Qnp[k], 1, Qnp[k], 1);
                 normval = sqrt(normval);
                 //outfile->Printf("trial vector norm%30.18lf\n",normval);
                 if (normval > 1e-20) {
                     for (I = 0; I < N; I++) {
                         Qp[L][I] = Qnp[k][I] / normval;
                     }
                     L++;
                     //outfile->Printf("check orthogonality1\n");
                     for (int i = 0; i < L; i++) {
                         for (int j = 0; j < L; j++) {
                             double a = C_DDOT(N, Qp[i], 1, Qp[j], 1);
                             //outfile->Printf("%d %d %30.16lf",i,j,a);
                             if ((i!=j) && (fabs(a)>1e-12)) outfile->Printf(" detect linear dependency\n");
                             //outfile->Printf("\n");
                         }
                     }
                 }
             }
	     */
	     //Lapack Householder
             // copy new vectors into place 
             for (size_t i = 0; i < 2 * Mn; i++) {
                 for (size_t k = 0; k < N; k++) {
                     Qp[i][k] = Qnp[i][k];
             }
             }
              // Q->print();

             L = 2 * Mn;
             qr_orthogonalization(maxdim,N,L,Qp);
          }else {
             for (size_t i = 0; i < maxdim; i++) {
                 for (size_t k = 0; k < N; k++) {
                     Qnp[i][k] = 0.0;
                 }
             }
             for (size_t k = 0; k < Mn; k++) {
                 for (size_t n = 0; n < N; n++) {
                     Qnp[k][n]= Rrp[k][n];
                 }
             }
	     /*
	     //Gram-Schmidt 
             double dotval, normval;
             size_t k,i, I;
             for (k = 0; k < M; k++) {
                 //outfile->Printf("L %d\n",L);
                 for (i = 0; i < L; i++) {
                     dotval = C_DDOT(N, Qp[i], 1, Qnp[k], 1);
                     for (I = 0; I < N; I++) Qnp[k][I] -= dotval * Qp[i][I];
                 }
                 //reorthogonalization
                 for (i = 0; i < L; i++) {
                     dotval = C_DDOT(N, Qp[i], 1, Qnp[k], 1);
                     for (I = 0; I < N; I++) Qnp[k][I] -= dotval * Qp[i][I];
                 }
                 normval = C_DDOT(N, Qnp[k], 1, Qnp[k], 1);
                 normval = sqrt(normval);
                 //outfile->Printf("trial vector norm%30.18lf\n",normval);
                 if (normval > 1e-20) {
                     for (I = 0; I < N; I++) {
                         Qp[L][I] = Qnp[k][I] / normval;
                     }
                     L++;
                     //outfile->Printf("check orthogonality1\n");
                     for (int i = 0; i < L; i++) {
                         for (int j = 0; j < L; j++) {
                             double a = C_DDOT(N, Qp[i], 1, Qp[j], 1);
                             //outfile->Printf("%d %d %30.16lf",i,j,a);
                             if (i!=j && fabs(a)>1e-12) outfile->Printf(" detect linear dependency\n");
                             //outfile->Printf("\n");
                         }
                     }
                 }
             }
             */
	     //Lapack Householder
             for (size_t k = 0; k < Mn; k++) {
                 for (size_t n = 0; n < N; n++) {
                     Qp[L+k][n]= Rrp[k][n];
                 }
             }
             L=L + Mn;
             qr_orthogonalization(maxdim,N,L,Qp);
             //Q->print();
          }
      }
      iter++;
      //outfile->Printf("%d\n",L);
    free(energy_is_real);
    }
    if (convergence == false) {
            outfile->Printf("Davidson procedure fails after %d iterations, try increase maxdim", maxit);
            exit(0);
    }
    //get the eigenvectors 
    for (size_t k = 0; k < M; k++) {
         reval[k]=1.0/lambdap[k];
         for (size_t n = 0; n < N; n++) {
             rer[k][n] = QCp[k][n];
         }
    }

}
void Nonsym_DavidsonSolver::eigenvectors_swap(double **L, double **R, int l, int n, int m) {
    //this function is tested with simple systems, further tests may be required
    //l=row l of matrix of eigenvectors where degeneracy happens
    //n=number of degeneracy
    //m=length of each eigenvectors
    //selection sort
    double Smax, S,t;
    int i, j, jmax;
    for (j=0; j<n-1; j++) {
        Smax = C_DDOT(m, L[l+j], 1, R[l+j], 1);
        jmax = j; 
        // find maximum eigenvalue
        for (i=j+1; i<n; i++) {
            S = C_DDOT(m, L[l+j], 1, R[l+i], 1);
           outfile->Printf("%d %d %20.12lf %20.12lf \n",j,i,Smax,S);
            if (fabs(Smax) - fabs(S) < 0e0) {
                jmax = i; Smax = S;
            }
        }
        if (jmax != j) {
         //  outfile->Printf("detect wrong order\n");
           for (i=0; i<m; i++) {
               t = R[l+j][i]; R[l+j][i] = R[l+jmax][i]; R[l+jmax][i] = t;
           }
        }
    }
}
void Nonsym_DavidsonSolver::schmidt_biorthogonal(double **A, double **B, size_t m, size_t n) {
    //assuming that A and B have m linear independent row vectors and A[i]*B[i] is non-zero, this will biorthogonalize the rows of A and B
    double RValue1;
    double RValue2;
    for (size_t i = 0; i < m; i++) {
        RValue1 = C_DDOT(n, A[i], 1, B[i], 1);
        if (fabs(RValue1) < 1e-6) {
                outfile->Printf ("vector A[%d] and  B[%d] are already orthogonal.", i,i);//possibly need to choose new B[i] from other rows so that A[i] *B[i] \neq 0
                exit(0);
        }
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

        for (size_t j = i+1; j < m; j++) {
            RValue1 = C_DDOT(n, A[i], 1, B[j], 1);
            RValue2 = C_DDOT(n, B[i], 1, A[j], 1);
            for (size_t I = 0; I < n; I++) {
                B[j][I] -= RValue1 * B[i][I];
                A[j][I] -= RValue2 * A[i][I];
            }
        }
    }

}
}
