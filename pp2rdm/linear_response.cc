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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pp2rdm_solver.h"

#include <psi4/libpsi4util/PsiOutStream.h>

#include <psi4/libtrans/integraltransform.h>
#include<psi4/libmints/mintshelper.h>
#include "blas.h"
#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() ( (double)clock() / CLOCKS_PER_SEC )
    #define omp_get_max_threads() 1
#endif

using namespace psi;
using namespace pp2rdm;
using namespace fnocc;

namespace psi { namespace pp2rdm {


double pp2RDMSolver::evaluate_hessian_pqrs(int p, int q, int r, int s) {

    double val = 0.0;

    int n1 = nmo_;
    int n2 = nmo_ * nmo_;
    int n3 = nmo_ * nmo_  * nmo_;
    int ntri = nmo_*(nmo_+1)/2;

    //one-electron terms
    val += -oei_[INDEX( s , p )] * d1_[INDEX( q , r )] - oei_[INDEX( q , r )] * d1_[INDEX( s , p )];

    if ( q == r ) {
        for (int u = 0; u < nmo_; u++) {
            val += 0.5 * (oei_[INDEX( u , p )] * d1_[INDEX( s , u )] + oei_[INDEX( s , u )] * d1_[INDEX( u , p )]);
        }
    }
    if ( p == s ) {
        for (int u = 0; u < nmo_; u++) {
            val += 0.5 * (oei_[INDEX( u , r )] * d1_[INDEX( q , u )] + oei_[INDEX( q , u )] * d1_[INDEX( u , r )]);
        }
    }

    //two-electron terms
    if ( q == r ) {
        for (int t = 0; t < nmo_; t++) {
            for (int u = 0; u < nmo_; u++) {
                for (int v = 0; v < nmo_; v++) {
                    val += 0.5 * tei_4index_[INDEX( u , p ) * ntri + INDEX( v , t )] * d2t_[ s * n3 + t * n2 + u * n1 + v ];
                    val += 0.5 * tei_4index_[INDEX( s , u ) * ntri + INDEX( t , v )] * d2t_[ u * n3 + v * n2 + p * n1 + t ];
                }
            }
        }
    }
    if ( p == s ) {
        for (int t = 0; t < nmo_; t++) {
            for (int u = 0; u < nmo_; u++) {
                for (int v = 0; v < nmo_; v++) {
                    val += 0.5 * tei_4index_[INDEX( q , u ) * ntri + INDEX( t , v )] * d2t_[ u * n3 + v * n2 + r * n1 + t ];
                    val += 0.5 * tei_4index_[INDEX( u , r ) * ntri + INDEX( v , t )] * d2t_[ q * n3 + t * n2 + u * n1 + v ];
                }
            }
        }
    }
    for (int u = 0; u < nmo_; u++) {
        for (int v = 0; v < nmo_; v++) {
            val += tei_4index_[INDEX( u , p ) * ntri + INDEX( v , r )] * d2t_[ q * n3 + s * n2 + u * n1 + v ];
            val += tei_4index_[INDEX( q , u ) * ntri + INDEX( s , v )] * d2t_[ u * n3 + v * n2 + p * n1 + r ];
        }
    }
    for (int u = 0; u < nmo_; u++) {
        for (int t = 0; t < nmo_; t++) {
            val -= tei_4index_[INDEX( s , p ) * ntri + INDEX( t , u )] * d2t_[ q * n3 + u * n2 + r * n1 + t ];
            val -= tei_4index_[INDEX( t , p ) * ntri + INDEX( s , u )] * d2t_[ q * n3 + u * n2 + t * n1 + r ];
            val -= tei_4index_[INDEX( q , r ) * ntri + INDEX( u , t )] * d2t_[ s * n3 + t * n2 + p * n1 + u ];
            val -= tei_4index_[INDEX( q , t ) * ntri + INDEX( u , r )] * d2t_[ t * n3 + s * n2 + p * n1 + u ];
        }
    }

    return val;
}

void pp2RDMSolver::evaluate_a11_b11(std::shared_ptr<Matrix> a11, std::shared_ptr<Matrix> b11) {

    //BuildSeniorityZeroRDMs(false);
    BuildHFRDMs(false);
    int o = nalpha_;
    int v = nmo_ - nalpha_;

    //spin-free 2RDM
    memset((void*)d2t_,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    for (int p = 0; p < nmo_; p++) {
        for (int q = 0; q < nmo_; q++) {
            for (int r = 0; r < nmo_; r++) {
                for (int s = 0; s < nmo_; s++) {
                     d2t_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s] = 2.0 * d2ab_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s]
                                                                 + 2.0 * d2aa_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s];
                }
            }
        }
    }

    double ** ap = a11->pointer();
    double ** bp = b11->pointer();
    for (int i = 0; i < o; i++) { 
        for (int a = 0; a < v; a++) { 
            for (int j = 0; j < o; j++) { 
                for (int b = 0; b < v; b++) { 
                    ap[i*v+a][j*v+b] =  evaluate_hessian_pqrs(i,a+o,b+o,j);
                    bp[i*v+a][j*v+b] = -evaluate_hessian_pqrs(i,a+o,j,b+o);
                }
            }
        }
    }

/*
    int n3=nmo_*nmo_*nmo_;
    int n2=nmo_*nmo_;
    int n1=nmo_;
    //eq. A6 from J. Chem. Phys. 141, 244104 (2014)
    outfile->Printf("Electronic Hessian submatrix A from RDMs\n");
    for (int a = 0; a < v; a++) {
        for (int i = 0; i < o; i++) {
            for (int b = 0; b < v; b++) {
                for (int j = 0; j < o; j++) {

                    // int ai = (2*i*nmo_-i*i+2*(a+o)-3*i-2)/2;    
                    // int bj = (2*j*nmo_-j*j+2*(b+o)-3*j-2)/2;    
                    // hessian->pointer()[a*o+i][b*o+j]= orbital_hessian[INDEX(ai,bj)];
                    //hessian->pointer()[a*o+i][b*o+j]= orbital_hessian[INDEX(INDEX(a,i),INDEX(b,j))];
                    // hessian->pointer()[a*o+i][b*o+j]= 1.0;
                    double dum1 =0.0;
                    double dum = -oei_[INDEX(j,i)] * d1_[INDEX(b+o,a+o)] - oei_[INDEX(a+o,b+o)] * d1_[INDEX(i,j)];
                    //one-electron
                    for (int u = 0; u < nmo_; u++) {
                        dum += 0.5 * (oei_[INDEX(u,i  )] * d1_[INDEX(u,  j)] + oei_[INDEX(j  ,u)] * d1_[INDEX(i  ,u)]) * (a==b);
                        dum += 0.5 * (oei_[INDEX(u,b+o)] * d1_[INDEX(u,a+o)] + oei_[INDEX(a+o,u)] * d1_[INDEX(b+o,u)]) * (i==j);
                        //two-electron
                        for (int w = 0; w < nmo_; w++) {
                            dum += tei_4index_[INDEX(u  ,i  )*nmo_*(nmo_+1)/2+INDEX(w,b+o)] * d2t_[    u*n3+    w*n2+(a+o)*n1+j];
                            dum += tei_4index_[INDEX(a+o,u  )*nmo_*(nmo_+1)/2+INDEX(j,w  )] * d2t_[    i*n3+(b+o)*n2+    u*n1+w];
                            dum -= tei_4index_[INDEX(j  ,i  )*nmo_*(nmo_+1)/2+INDEX(w,u  )] * d2t_[(b+o)*n3+    w*n2+(a+o)*n1+u];
                            dum -= tei_4index_[INDEX(w  ,i  )*nmo_*(nmo_+1)/2+INDEX(j,u  )] * d2t_[    w*n3+(b+o)*n2+(a+o)*n1+u];
                            dum -= tei_4index_[INDEX(a+o,b+o)*nmo_*(nmo_+1)/2+INDEX(u,w  )] * d2t_[    i*n3+    u*n2+    j*n1+w];
                            dum -= tei_4index_[INDEX(a+o,w  )*nmo_*(nmo_+1)/2+INDEX(u,b+o)] * d2t_[    i*n3+    u*n2+    w*n1+j];
                            for (int t = 0; t < nmo_; t++) {
                                dum += 0.5 * tei_4index_[INDEX(u  ,i  )*nmo_*(nmo_+1)/2+INDEX(w,t)] * d2t_[u*n3+w*n2+j*n1+t] * (a==b);
                                dum += 0.5 * tei_4index_[INDEX(j  ,u  )*nmo_*(nmo_+1)/2+INDEX(t,w)] * d2t_[i*n3+t*n2+u*n1+w] * (a==b);
                                dum += 0.5 * tei_4index_[INDEX(a+o,u  )*nmo_*(nmo_+1)/2+INDEX(t,w)] * d2t_[(b+o)*n3+t*n2+u*n1+w] * (i==j);
                                dum += 0.5 * tei_4index_[INDEX(u  ,b+o)*nmo_*(nmo_+1)/2+INDEX(w,t)] * d2t_[u*n3+w*n2+(a+o)*n1+t] * (i==j);
                            }
                        }
                    }
                        a11->pointer()[i*v+a][j*v+b]= -dum/2.0;//normalization
                        outfile->Printf("%20.12lf", a11->pointer()[i*v+a][j*v+b]);
                }
            }
            //outfile->Printf("\n");
        }
    }

    outfile->Printf("Electronic Hessian submatrix B from RDMs\n");
    for (int a = 0; a < v; a++) {
        for (int i = 0; i < o; i++) {
            for (int b = 0; b < v; b++) {
                for (int j = 0; j < o; j++) {
                    //one-electron
                    double dum = -oei_[INDEX(b+o,i)]*d1_[INDEX(j,a+o)] - oei_[INDEX(a+o,j)]*d1_[INDEX(i,b+o)];
                    for (int u = 0; u < nmo_; u++) {
                        dum += 0.5 * (oei_[INDEX(u,i)] * d1_[INDEX(u,b+o)] + oei_[INDEX(b+o,u)] * d1_[INDEX(i,u)]) * (a+o==j);
                        dum += 0.5 * (oei_[INDEX(u,j)] * d1_[INDEX(u,a+o)] + oei_[INDEX(a+o,u)] * d1_[INDEX(j,u)])*(i==b+o);
                        //two-electron
                        for (int w = 0; w < nmo_; w++) {
                            dum += tei_4index_[INDEX(u,i)*nmo_*(nmo_+1)/2+INDEX(w,j)] * d2t_[u*n3+w*n2+(a+o)*n1+b+o];
                            dum += tei_4index_[INDEX(a+o,u)*nmo_*(nmo_+1)/2+INDEX(b+o,w)] * d2t_[i*n3+j*n2+u*n1+w];
                            dum -= tei_4index_[INDEX(b+o,i)*nmo_*(nmo_+1)/2+INDEX(w,u)] * d2t_[j*n3+w*n2+(a+o)*n1+u];
                            dum -= tei_4index_[INDEX(w,i)*nmo_*(nmo_+1)/2+INDEX(b+o,u)] * d2t_[w*n3+j*n2+(a+o)*n1+u];
                            dum -= tei_4index_[INDEX(a+o,j)*nmo_*(nmo_+1)/2+INDEX(u,w)] * d2t_[i*n3+u*n2+(b+o)*n1+w];
                            dum -= tei_4index_[INDEX(a+o,w)*nmo_*(nmo_+1)/2+INDEX(u,j)] * d2t_[i*n3+u*n2+w*n1+(b+o)];
                            for (int t = 0; t < nmo_; t++) { 
                                dum += 0.5 * tei_4index_[INDEX(u,i)*nmo_*(nmo_+1)/2+INDEX(w,t)] * d2t_[u*n3+w*n2+(b+o)*n1+t] * (a+o==j);
                                dum += 0.5 * tei_4index_[INDEX(b+o,u)*nmo_*(nmo_+1)/2+INDEX(t,w)] * d2t_[i*n3+t*n2+u*n1+w] * (a+o==j);
                                dum += 0.5 * tei_4index_[INDEX(a+o,u)*nmo_*(nmo_+1)/2+INDEX(t,w)] * d2t_[j*n3+t*n2+u*n1+w] * (i==b+o);
                                dum += 0.5 * tei_4index_[INDEX(u,j)*nmo_*(nmo_+1)/2+INDEX(w,t)] * d2t_[u*n3+w*n2+(a+o)*n1+t] * (i==b+o);
                            }
                        }
                    }
                        b11->pointer()[i*v+a][j*v+b]= -dum/2.0;
                    outfile->Printf("%20.12lf", b11->pointer()[i*v+a][j*v+b]);
                }
            }
            outfile->Printf("\n");
        }
    }
            outfile->Printf("(A-B) * 4\n");
    for (int a = 0; a < v; a++) {
        for (int i = 0; i < o; i++) {
            for (int b = 0; b < v; b++) {
                for (int j = 0; j < o; j++) {
                    double dum=a11->pointer()[i*v+a][j*v+b]- b11->pointer()[i*v+a][j*v+b];
                    outfile->Printf("%20.12lf", dum*4);
                }
            }
            outfile->Printf("\n");
            
        }
    }
*/

    // compare greg's hessian to A11+B11

    PackSpatialDensity();

    int * symmetry_energy_order  = (int*)malloc(nmo_*sizeof(int));
    for (int i = 0; i < nmo_; i++) {
        symmetry_energy_order[i] = 1;
    }

    int nrstc = 0;
    int nrstv = 0;

    int save_val = orbopt_data_[8];
    orbopt_data_[8] = -3;

    int ntri = nmo_*(nmo_-1)/2;
    double * orbital_hessian = (double*)malloc(ntri*(ntri+1)/2*sizeof(double));
    memset((void*)orbital_hessian,'\0',ntri*(ntri+1)/2*sizeof(double));

    //this require the stop statement in line 713 of focas_hessian.F90 to be uncommented, otherwise the program will stop
    OrbOpt(orbopt_transformation_matrix_,
          oei_, oei_dim_,
          tei_, tei_dim_,
          d1_, oei_dim_, d2_, d2_dim_,
          symmetry_energy_order, nrstc, nmo_, nrstv, nirrep_,
          orbopt_data_, orbopt_outfile_, orbital_hessian);

    free(symmetry_energy_order);

    int ** idx_map = (int**)malloc(nmo_*sizeof(int*));
    for (int p = 0; p < nmo_; p++) {
        idx_map[p] = (int*)malloc(nmo_*sizeof(int));
    }
    int idx = 0;
    for (int p = 0; p < nmo_; p++) {
        for (int q = p+1; q < nmo_; q++) {
            idx_map[p][q] = idx;
            idx_map[q][p] = idx;
            idx++;
        }
    }

    for (int p = 0; p < nmo_; p++) {
        for (int q = p+1; q < nmo_; q++) {

            int pq = idx_map[p][q];

            for (int r = 0; r < nmo_; r++) {
                //int smin = r+1;
                //if ( r == p ) smin = q;
                for (int s = r+1; s < nmo_; s++) {

                    int rs = idx_map[r][s];

                    if ( pq != rs ) continue;
                    if ( pq < rs ) continue;
                    //if ( pq > rs ) continue;

                    double val1 = -evaluate_hessian_pqrs(p,q,s,r);
                    double val2 =  evaluate_hessian_pqrs(p,q,r,s);
                    double val3 =  evaluate_hessian_pqrs(q,p,s,r);
                    double val4 = -evaluate_hessian_pqrs(q,p,r,s);

                    double val5 = -evaluate_hessian_pqrs(s,r,p,q);
                    double val6 =  evaluate_hessian_pqrs(r,s,p,q);
                    double val7 =  evaluate_hessian_pqrs(s,r,q,p);
                    double val8 = -evaluate_hessian_pqrs(r,s,q,p);

                    double hess  = val1 + val2 + val3 + val4;
                    double hess2 = val5 + val6 + val7 + val8;

                    //double ref  = orbital_hessian[count++];
                    //double ref  = orbital_hessian[INDEX(pq,rs)];
                    double ref  = orbital_hessian[INDEX(pq,rs)];

                    printf("%5i %5i %5i %5i %20.12lf %20.12lf %20.12lf %20.12lf\n",p,q,r,s,ref,hess,hess2,ref-hess);
                    fflush(stdout);
                }
            }
        }
    }

    free(orbital_hessian);
    orbopt_data_[8] = save_val;

    free(d2t_);
}

void pp2RDMSolver::evaluate_a21_b21(std::shared_ptr<Matrix> a21, std::shared_ptr<Matrix> b21) {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    double * tmp_p1 = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * tmp_p2 = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * tmp_m1 = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * tmp_m2 = (double*)malloc(nmo_*nmo_*sizeof(double));

    double **a21_p = a21->pointer();
    double **b21_p = b21->pointer();

    std::shared_ptr<Matrix> orbital_gradient (new Matrix(nmo_,nmo_));
    double ** orbital_gradient_p = orbital_gradient->pointer();

    double h = 1e-4;
    for (int ia = 0; ia < o*v; ia++) {

        double save = t2_[ia];

        t2_[ia] = save + h;
        double dum = BuildSeniorityZeroRDMs(false);
        evaluate_orbital_gradient(orbital_gradient);
        C_DCOPY(nmo_*nmo_,orbital_gradient_p[0],1,tmp_p1,1);

        t2_[ia] = save + 2.0 * h;
        dum = BuildSeniorityZeroRDMs(false);
        evaluate_orbital_gradient(orbital_gradient);
        C_DCOPY(nmo_*nmo_,orbital_gradient_p[0],1,tmp_p2,1);

        t2_[ia] = save - h;
        dum = BuildSeniorityZeroRDMs(false);
        evaluate_orbital_gradient(orbital_gradient);
        C_DCOPY(nmo_*nmo_,orbital_gradient_p[0],1,tmp_m1,1);

        t2_[ia] = save - 2.0 * h;
        dum = BuildSeniorityZeroRDMs(false);
        evaluate_orbital_gradient(orbital_gradient);
        C_DCOPY(nmo_*nmo_,orbital_gradient_p[0],1,tmp_m2,1);

        for (int j = 0; j < o; j++) {
            for (int b = 0; b < v; b++) {
                a21_p[ia][j*v+b] = (       - tmp_p2[j*nmo_+(b+o)] 
                                     + 8.0 * tmp_p1[j*nmo_+(b+o)] 
                                     - 8.0 * tmp_m1[j*nmo_+(b+o)] 
                                           + tmp_m2[j*nmo_+(b+o)]) / (12.0 * h);
            }
        }
        for (int j = 0; j < o; j++) {
            for (int b = 0; b < v; b++) {
                b21_p[ia][j*v+b] = -(       - tmp_p2[(b+o)*nmo_+j] 
                                      + 8.0 * tmp_p1[(b+o)*nmo_+j] 
                                      - 8.0 * tmp_m1[(b+o)*nmo_+j] 
                                            + tmp_m2[(b+o)*nmo_+j]) / (12.0 * h);
            }
        }

        t2_[ia] = save;

    }

    free(tmp_p1);
    free(tmp_p2);
    free(tmp_m1);
    free(tmp_m2);

}

// evaluates orbital gradient (actually 1/2 of orbital gradient)
// g(pq) = <0| [E_{pq}^- , H ] |0> = 2 <0| [E_{pq} , H ] |0>
// i am only returning <0| [ E_{pq} , H ] |0>, though
void pp2RDMSolver::evaluate_orbital_gradient(std::shared_ptr<Matrix> gradient) {

    gradient->zero();
    double ** gradient_p = gradient->pointer();

    BuildSeniorityZeroRDMs(false);

    for (int p = 0; p < nmo_; p++) {
        for (int q = 0; q < nmo_; q++) {
            double dum = 0.0;
            for (int r = 0; r < nmo_; r++) {
                dum += oei_[INDEX(r,p)] * d1_[INDEX(q,r)];
                dum -= oei_[INDEX(q,r)] * d1_[INDEX(r,p)];
            }
            gradient_p[p][q] = dum;
        }
    }
    /*for (int p = 0; p < nmo_; p++) {
        for (int r = 0; r < nmo_; r++) {
            for (int s = 0; s < nmo_; s++) {
                for (int t = 0; t < nmo_; t++) {
                    int rp = INDEX(r,p);
                    int st = INDEX(s,t);
                    double v_rspt = C_DDOT(nQ_,Qmo_ + rp,nmo_*(nmo_+1)/2,Qmo_+st,nmo_*(nmo_+1)/2);
                    for (int q = 0; q < nmo_; q++) {
                        double d2ab = d2ab_[q*nmo_*nmo_*nmo_ + t*nmo_*nmo_ + r*nmo_ + s];
                        double d2aa = d2aa_[q*nmo_*nmo_*nmo_ + t*nmo_*nmo_ + r*nmo_ + s];
                        double d2   = 2.0 * (d2ab + d2aa);
                        gradient_p[p][q] += v_rspt * d2;
                        gradient_p[q][p] -= v_rspt * d2;
                    }
                }
            }
        }
    }*/
    for (int p = 0; p < nmo_; p++) {
        for (int q = 0; q < nmo_; q++) {
            for (int t = 0; t < nmo_; t++) {
                int tp = INDEX(t,p);
                int tq = INDEX(t,q);
                int qt = INDEX(q,t);
                int qp = INDEX(q,p);
                int tt = INDEX(t,t);
                double v_ttpq = C_DDOT(nQ_,Qmo_ + tp,nmo_*(nmo_+1)/2,Qmo_+tq,nmo_*(nmo_+1)/2);
                double v_qtpt = C_DDOT(nQ_,Qmo_ + qp,nmo_*(nmo_+1)/2,Qmo_+tt,nmo_*(nmo_+1)/2);

                double d2ab = d2ab_[q*nmo_*nmo_*nmo_ + t*nmo_*nmo_ + q*nmo_ + t];
                double d2aa = d2aa_[q*nmo_*nmo_*nmo_ + t*nmo_*nmo_ + q*nmo_ + t];
                double d2   = 2.0 * (d2ab + d2aa);
                gradient_p[p][q] += v_qtpt * d2;
                gradient_p[q][p] -= v_qtpt * d2;

                if ( q != t ) {
                    d2ab = d2ab_[q*nmo_*nmo_*nmo_ + q*nmo_*nmo_ + t*nmo_ + t];
                    d2aa = d2aa_[q*nmo_*nmo_*nmo_ + q*nmo_*nmo_ + t*nmo_ + t];
                    d2   = 2.0 * (d2ab + d2aa);
                    gradient_p[p][q] += v_ttpq * d2;
                    gradient_p[q][p] -= v_ttpq * d2;

                    d2ab = d2ab_[q*nmo_*nmo_*nmo_ + t*nmo_*nmo_ + t*nmo_ + q];
                    d2aa = d2aa_[q*nmo_*nmo_*nmo_ + t*nmo_*nmo_ + t*nmo_ + q];
                    d2   = 2.0 * (d2ab + d2aa);
                    gradient_p[p][q] += v_ttpq * d2;
                    gradient_p[q][p] -= v_ttpq * d2;
                }

            }
        }
    }

    //for (int q = 0; q < nmo_; q++) {
    //    for (int r = 0; r < nmo_; r++) {
    //        for (int s = 0; s < nmo_; s++) {
    //            for (int t = 0; t < nmo_; t++) {
    //                int qr = INDEX(q,r);
    //                int ts = INDEX(t,s);
    //                double v_qtrs = C_DDOT(nQ_,Qmo_ + qr,nmo_*(nmo_+1)/2,Qmo_+ts,nmo_*(nmo_+1)/2);
    //                for (int p = 0; p < nmo_; p++) {
    //                    double d2ab = d2ab_[r*nmo_*nmo_*nmo_ + s*nmo_*nmo_ + p*nmo_ + t];
    //                    double d2aa = d2aa_[r*nmo_*nmo_*nmo_ + s*nmo_*nmo_ + p*nmo_ + t];
    //                    //gradient_p[p][q] -= 2.0 * v_qtrs * (d2aa + d2ab);
    //                }
    //            }
    //        }
    //    }
    //}

    //double save_value = orbopt_data_[8];
    //orbopt_data_[8] = -1.0;
    //RotateOrbitals();
    //orbopt_data_[8] = save_value;

    //for (int p = 0; p < nmo_; p++) {
    //    for (int q = 0; q < nmo_; q++) {
    //        if ( p == q ) continue;
    //        printf("%5i %5i %20.12lf %20.12lf %20.12lf\n",p,q,orbital_lagrangian_[p*nmo_+q],
    //            gradient_p[q][p],orbital_lagrangian_[p*nmo_+q]-gradient_p[q][p]);
    //    }
    //}

    //double nrm = 0.0;
    //for (int p = 0; p < nmo_; p++) {
    //    for (int q = 0; q < p; q++) {
    //        nrm += 4.0 * gradient_p[p][q] * gradient_p[p][q];
    //    }
    //}
    //printf("%20.12lf %20.12lf\n",orbopt_data_[11],sqrt(nrm));

}
void pp2RDMSolver::evaluate_metric_matrix(std::shared_ptr<Matrix> metric) {
    // need to test again, I built this long ago
    int o = nalpha_;
    int v = nmo_ - nalpha_;
    // orbital-orbital block (M_(1-1)): S11_(ia,ia) = D_ii-D_aa=2-2*sum t_ib*t_ib-2*sum t_aj*t_aj
    double * S = (double*)malloc(o*v*o*v*sizeof(double));
    double * Dij = (double*)malloc(o*o*sizeof(double));
    double * Dab = (double*)malloc(v*v*sizeof(double));
    memset((void*)S,'\0',o*v*o*v*sizeof(double));
    memset((void*)Dij,'\0',o*o*sizeof(double));
    memset((void*)Dab,'\0',v*v*sizeof(double));
    for (int i = 0; i < o; i++) {
        double dum = 0.0;
        for (int a = 0; a < v; a++) {
            dum += t2_[i*v+a]*t2_[i*v+a];
            }
        Dij[i*o+i]=1.0-1.0*dum;
    }

    for (int a = 0; a < v; a++) {
        double dum = 0.0;
        for (int i = 0; i < o; i++) {
            dum += t2_[i*v+a]*t2_[i*v+a];
        }
        Dab[a*o+a]=1.0*dum;
    }

    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            int ia = i*v+a;
            S[ia*o*v+ia]= Dij[i*o+i]-Dab[a*v+a];
        }
    }
//amplitude-amplitude block is an identiy matrix M_(2-2)=<HF|[a2+,a2]|HF> 
// M
    double **M = metric->pointer();

    // S(ia,jb) = dab gamma_ij - dij gamma_ab
    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            for (int a = 0; a < v; a++) {
                for (int b = 0; b < v; b++) {
                    M[i*v + a][j*v + b]             =   (a==b) * d1_[INDEX(i,j)] - (i==j) * d1_[INDEX(a+o,b+o)];
                    M[o*v + i*v + a][o*v + j*v + b] = - (a==b) * d1_[INDEX(i,j)] + (i==j) * d1_[INDEX(a+o,b+o)];
                }
            }
        }
    }

    //for (int p = 0; p < 2*o*v; p++) {
    //    for (int q = 0; q < 2*o*v; q++) {
    //        if (p < o*v && q < o*v) M[p][q] = S[p*o*v+q];
    //        else if (p >= o*v && q >= o*v) M[p][q] = (p==q);
    //        else M[p][q] = 0.0;
    //    }
    //}
    free(S);
    free(Dij);
    free(Dab);
}

double pp2RDMSolver::BuildHFRDMs(bool print){

    memset((void*)d1_,'\0', oei_dim_ *  sizeof(double));
    memset((void*)d2ab_,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)d2aa_,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
   // memset((void*)d2t_,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    double energy = 0.0;

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            d2ab_[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j] = 1.0;
        }
    }

 // construct d2aa from d2ab:
    for (int p = 0; p < nmo_; p++) {
        for (int q = 0; q < nmo_; q++) {
            for (int r = 0; r < nmo_; r++) {
                for (int s = 0; s < nmo_; s++) {
                    d2aa_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s] = d2ab_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s]
                                                                 - d2ab_[q*nmo_*nmo_*nmo_+p*nmo_*nmo_+r*nmo_+s];

                }
            }
        }
    }

    // check trace
    double tra = 0.0;
    double trb = 0.0;
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            tra += d2aa_[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j];
            trb += d2ab_[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j];
        }
    }

    // build d1 by contraction
    for (int i = 0; i < o; i++) {
        d1_[INDEX(i,i)] = 2.0;
    }
    //for (int i = 0; i < nmo_; i++) {
    //    for (int j = 0; j < nmo_; j++) {

    //        double dum = 0.0;
    //        for (int p = 0; p < nmo_; p++) {
    //            dum += d2aa_[i*nmo_*nmo_*nmo_+p*nmo_*nmo_+j*nmo_+p];
    //            dum += d2ab_[i*nmo_*nmo_*nmo_+p*nmo_*nmo_+j*nmo_+p];
    //        }
    //        d1_[INDEX(i,j)] = 2.0 * dum / (2.0 * nalpha_-1.0);
    //    //    outfile->Printf("%20.12lf, %d %d ",d1_[INDEX(i,j)], i, j);

    //    }
    //  //  outfile->Printf("\n");
    //}

    // check energy

    double e1 = 0.0;
    for (int i = 0; i < nmo_; i++) {
        e1 += d1_[INDEX(i,i)] * oei_[INDEX(i,i)];
    }


    double e2 = 0.0;

    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            e2 += d2ab_[i * nmo_*nmo_*nmo_ + j * nmo_*nmo_ + i * nmo_ + j] * tei_4index_[INDEX(i,i)*nmo_*(nmo_+1)/2+INDEX(j,j)] ;
            e2 += d2aa_[i * nmo_*nmo_*nmo_ + j * nmo_*nmo_ + i * nmo_ + j] * tei_4index_[INDEX(i,i)*nmo_*(nmo_+1)/2+INDEX(j,j)] ;
            e2 += d2aa_[i * nmo_*nmo_*nmo_ + j * nmo_*nmo_ + j * nmo_ + i] * tei_4index_[INDEX(i,j)*nmo_*(nmo_+1)/2+INDEX(j,i)] ;
        }
    }
   //spin-free 2RDM
  /* for (int p = 0; p < nmo_; p++) {
        for (int q = 0; q < nmo_; q++) {
            for (int r = 0; r < nmo_; r++) {
                for (int s = 0; s < nmo_; s++) {
                    d2t_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s] = 2.0 * d2ab_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s]
                                                                + 2.0 * d2aa_[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+r*nmo_+s];

                }
            }
        }
    }*/
    double en = e1 + e2;
    if ( print ) {
        outfile->Printf("\n");
        outfile->Printf("        HF one-electron energy = %20.12lf\n",e1);
        outfile->Printf("        HF two-electron energy = %20.12lf\n",e2);
        outfile->Printf("        * HF total energy      = %20.12lf\n",e1 + e2 + enuc_);fflush(stdout);
    }
      return e1 + e2;

}

double pp2RDMSolver::compute_excited_energy() {
    d1_ = (double*)malloc( nmo_ * (nmo_+1)/2 * sizeof(double));
    d2ab_ = (double*)malloc(nmo_ * nmo_ * nmo_ * nmo_ * sizeof(double));
    d2aa_ = (double*)malloc(nmo_ * nmo_ * nmo_ * nmo_ * sizeof(double));
    d2t_ = (double*)malloc(nmo_ * nmo_ * nmo_ * nmo_ * sizeof(double));
    //LR for HF    
    // grab the dipole integrals from MintsHelper:
    std::shared_ptr<MintsHelper> mints (new MintsHelper(reference_wavefunction_));
    std::vector<std::shared_ptr<Matrix> > dipole_ = mints->so_dipole();
    mux_a = (std::shared_ptr<Matrix>)(new Matrix(dipole_[0]));
    mux_b = (std::shared_ptr<Matrix>)(new Matrix(dipole_[0]));
    muy_a = (std::shared_ptr<Matrix>)(new Matrix(dipole_[1]));
    muz_a = (std::shared_ptr<Matrix>)(new Matrix(dipole_[2]));
    mux_a->transform(Ca_);
    muy_a->transform(Ca_);
    muz_a->transform(Ca_);
    int o = nalpha_;
    int v = nmo_ - nalpha_;
    Fij   = (double*)malloc(o*o*sizeof(double));
    Fab   = (double*)malloc(v*v*sizeof(double));
    tei_4index_ = (double*)malloc(nmo_*(nmo_+1)/2*nmo_*(nmo_+1)/2*sizeof(double));
    memset((void*)Fij,'\0',o*o*sizeof(double));
    memset((void*)Fab,'\0',v*v*sizeof(double));
    memset((void*)tei_4index_,'\0',nmo_*(nmo_+1)/2*nmo_*(nmo_+1)/2*sizeof(double));
    BuildFock();

    //build orbital-orbital Hessian
    std::shared_ptr<Matrix> A11 ( new Matrix(o*v,o*v) );
    std::shared_ptr<Matrix> B11 ( new Matrix(o*v,o*v) );
    evaluate_a11_b11(A11,B11);

   /* double * A = (double*)malloc(o*v*o*v*sizeof(double));
    double * B = (double*)malloc(o*v*o*v*sizeof(double));
    memset((void*)A,'\0',o*v*o*v*sizeof(double));
    memset((void*)B,'\0',o*v*o*v*sizeof(double));
   */
    /*
    outfile->Printf("\n");
    outfile->Printf("Electronic Hessian submatrix A\n");
    for (int a = 0; a < v; a++) {
        for (int i = 0; i < o; i++) {
            for (int b = 0; b < v; b++) {
                for (int j = 0; j < o; j++) {
                     int ia = i*v+a;
                     int jb = j*v+b;
                     A[ia * o * v + jb] = Fab[a * v + b] * (i == j) - Fij[i * o + j] * (a == b) - tei_4index_[INDEX(a+o,b+o)*nmo_*(nmo_+1)/2+INDEX(j,i)] + 2.0 * tei_4index_[INDEX(a+o,i)*nmo_*(nmo_+1)/2+INDEX(j,b+o)] ;
                    // outfile->Printf("\t%d%d%d%d: %20.12lf",a+o,i,b+o,j,A[ai * o * v + bj]);
                     outfile->Printf("%20.12lf",A[ia * o * v + jb]);
                }
            }
            outfile->Printf("\n");
        }
    }
            outfile->Printf("\n");
    outfile->Printf("Electronic Hessian submatrix B\n");
    for (int a = 0; a < v; a++) {
        for (int i = 0; i < o; i++) {
            for (int b = 0; b < v; b++) {
                for (int j = 0; j < o; j++) {
                     int ia = i*v+a;
                     int jb = j*v+b;
                     B[ia * o * v + jb] =  - 2.0 * tei_4index_[INDEX(a+o,i)*nmo_*(nmo_+1)/2+INDEX(b+o,j)] +  tei_4index_[INDEX(a+o,j)*nmo_*(nmo_+1)/2+INDEX(b+o,i)] ;
                    // outfile->Printf("\t%d%d%d%d: %20.12lf",a+o,i,b+o,j,B[ai * o * v + bj]);
                     outfile->Printf("%20.12lf",B[ia * o * v + jb]);
                }
            }
            outfile->Printf("\n");
        }
    }

    */
    long int totdim = 2*o*v;
    double * EH = (double*)malloc(totdim*totdim*sizeof(double));
    double * cc   = (double*)malloc(totdim*totdim*sizeof(double));
    double * eig  = (double*)malloc(totdim*sizeof(double));

    memset((void*)EH,'\0',totdim*totdim*sizeof(double));
    memset((void*)cc,'\0',totdim*totdim*sizeof(double));
    memset((void*)eig,'\0',totdim*sizeof(double));

    // evaluate metric matrix
    std::shared_ptr<Matrix> S (new Matrix(2*o*v,2*o*v));
    evaluate_metric_matrix(S);

    for (int p = 0; p < o*v; p++) {
        for (int q = 0; q < o*v; q++) {
            EH[p*totdim+q]            =  A11->pointer()[p][q];
            EH[(p+o*v)*totdim+q]      =  B11->pointer()[p][q];
            EH[(p+o*v)*totdim+(q+o*v)]=  A11->pointer()[p][q];
            EH[p*totdim+(q+o*v)]      =  B11->pointer()[p][q];

        }
    }

    int info = 0;
    std::shared_ptr<Matrix>  Xr (new Matrix(totdim,totdim/2));
    std::shared_ptr<Matrix>  Yr (new Matrix(totdim,totdim/2));
    std::shared_ptr<Matrix>  Xl (new Matrix(totdim,totdim/2));
    std::shared_ptr<Matrix>  Yl (new Matrix(totdim,totdim/2));
    bool prune = false;
    bool * energy_is_real = (bool*)malloc(totdim*sizeof(bool));
    if (prune) {
          info = SymmetricGeneralizedEigenvalueProblem(totdim,EH,S->pointer()[0],cc,eig);
    }else {
          GeneralizedEigenvalueProblem(totdim,EH,S->pointer()[0],cc,eig,energy_is_real);
    }
    Sort(S->pointer()[0], eig, totdim, 1);
    for (int p = 0; p < totdim; p++) {
        for (int a = 0; a < v; a++) {
            for (int i = 0; i < o; i++) { 
                Xr->pointer()[p][i*v+a] = S->pointer()[p][i*v+a      ]; 
                Yr->pointer()[p][i*v+a] = S->pointer()[p][i*v+a + o*v];
                Xl->pointer()[p][i*v+a] = EH[p*totdim + i*v+a]; 
                Yl->pointer()[p][i*v+a] = EH[p*totdim + i*v+a + o*v];
            }
       }
   }

   Xr->print();
   Yr->print();
  //  Xl->print();
  // Yl->print();
   outfile->Printf("state           eig                 X:            Y:         Z:             oscillator strengh\n");
   for (long int state = 0; state < totdim; state++) {
       if ( eig[state] > 0.0 ) {
          if ( prune ) eig[state] = 1.0 / eig[state];
          else if ( !energy_is_real[state] ) {
                  continue;
                  }
                // avoid considering nonsensical states
          if ( fabs(eig[state] ) > 10000.0) continue;
          double dum = 0.0;
          // normalization
          for (int a = 0; a < v; a++) {
              for (int i = 0; i < o; i++) {
                  dum += Xr->pointer()[state][i*v+a] * Xr->pointer()[state][i*v+a];
                  dum -= Yr->pointer()[state][i*v+a] * Yr->pointer()[state][i*v+a];
               }
          }
         //outfile->Printf("dum %20.12lf\n ", dum);
          dum=sqrt(dum);
          for (int a = 0; a < v; a++) {
              for (int i = 0; i < o; i++) {
                  Xr->pointer()[state][i*v+a] /= dum;
                  Yr->pointer()[state][i*v+a] /= dum;

              }
          }
          double dumx = 0.0;
          double dumy = 0.0;
          double dumz = 0.0;

         // outfile->Printf("state %5i\n ", state);
          for (int a = 0; a < v; a++) {
              for (int i = 0; i < o; i++) {
         // outfile->Printf("a %5i i%5i\n ", a,i);
                  dumx +=sqrt(2.0) * (Xr->pointer()[state][i*v+a]-Yr->pointer()[state][i*v+a]) * mux_a->pointer()[a+o][i];
                  dumy +=sqrt(2.0) * (Xr->pointer()[state][i*v+a]-Yr->pointer()[state][i*v+a]) * muy_a->pointer()[a+o][i];
                  dumz +=sqrt(2.0) * (Xr->pointer()[state][i*v+a]-Yr->pointer()[state][i*v+a]) * muz_a->pointer()[a+o][i];

              }
          }
           double val = 2./3. * eig[state] * (dumx*dumx+dumy*dumy+dumz*dumz);
         outfile->Printf("%5i %20.12lf    %10.6lf   %10.6lf   %10.6lf    %20.12lf\n",state, 27.21138*  eig[state], dumx, dumy, dumz,val);

   }
}
    free(energy_is_real);
    free(Fij);
    free(Fab);
    free(tei_4index_);
    //free(S);
    free(EH);
   // free(A);
   // free(B);
    return 1;
}

void pp2RDMSolver::BuildFock() {

    int o = nalpha_;
    int v = nmo_ - nalpha_;
    // Fij
    for (long int i = 0; i < o; i++) {
        for (long int j = 0; j < o; j++) {
            double dum = oei_[INDEX(i,j)];

            for (int k = 0; k < o; k++) {



                double coulomb = C_DDOT(nQ_,Qmo_ + INDEX(i,j),nmo_*(nmo_+1)/2,Qmo_+INDEX(k,k),nmo_*(nmo_+1)/2);

                double exchange = C_DDOT(nQ_,Qmo_ + INDEX(i,k),nmo_*(nmo_+1)/2,Qmo_+INDEX(k,j),nmo_*(nmo_+1)/2);

                dum += 2.0 * coulomb - exchange;

            }
            Fij[i * o + j] = dum;
        }
    }

    // Fab
    for (long int a = 0; a < v; a++) {
        for (long int b = 0; b < v; b++) {

            double dum = oei_[INDEX(a+o,b+o)];

            for (int k = 0; k < o; k++) {

                double coulomb = C_DDOT(nQ_,Qmo_ + INDEX(a+o,b+o),nmo_*(nmo_+1)/2,Qmo_+INDEX(k,k),nmo_*(nmo_+1)/2);

                double exchange = C_DDOT(nQ_,Qmo_ + INDEX(a+o,k),nmo_*(nmo_+1)/2,Qmo_+INDEX(k,b+o),nmo_*(nmo_+1)/2);

                dum += 2.0 * coulomb - exchange;
            }
            Fab[a * v + b] = dum;
        }
    }


   // construct 4-index integrals
  //     F_DGEMM('t','n',nmo_*(nmo_+1)/2,nmo_*(nmo_+1)/2,nQ_,1.0,Qmo_,nQ_,Qmo_,nQ_,0.0,tei_4index_,nmo_*(nmo_+1)/2);
   F_DGEMM('n','t',nmo_*(nmo_+1)/2,nmo_*(nmo_+1)/2,nQ_,1.0,Qmo_,nmo_*(nmo_+1)/2,Qmo_,nmo_*(nmo_+1)/2,0.0,tei_4index_,nmo_*(nmo_+1)/2);
}

void pp2RDMSolver::Sort(double *x, double *d, int n, int isort) {
//selection sort
    double dmax, t;
    int i, j, jmax;
    if (isort == 0) return;
    for (j=0; j<n-1; j++) {
        jmax = j; dmax = d[j];
        // find maximum eigenvalue
        for (i=j+1; i<n; i++) {
            if (isort * (dmax - d[i]) > 0e0) {
                jmax = i; dmax = d[i];
            }
        }
        if (jmax != j) {
           d[jmax] = d[j]; d[j] = dmax; // swap current component with maximum
           for (i=0; i<n; i++) {
               t = x[j*n+i]; x[j*n+i] = x[jmax*n+i]; x[jmax*n+i] = t;
           }
        }
     }
}

}}// end of namespaces
