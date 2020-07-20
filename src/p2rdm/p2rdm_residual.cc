/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include <psi4/psi4-dec.h>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libqt/qt.h>
#include <psi4/libmints/basisset.h>
#include <psi4/lib3index/3index.h>
#include <psi4/libciomr/libciomr.h>

#include <misc/blas.h>
#include <misc/omp.h>

#include "p2rdm_solver.h"

using namespace psi;
using namespace fnocc;

namespace hilbert {

void p2RDMSolver::evaluate_residual() {

    int o = nalpha_;
    int v = nmo_ - nalpha_;

    // need several buffers the size of t2 plus some other stuff

    long int                dim = 2L*v*v*v;
    if (2*nQ_*o*v     > dim) dim = 2*nQ_*o*v;
    if (o*o*v*v      > dim) dim = o*o*v*v;
    if (nQ_*v*v       > dim) dim = nQ_*v*v;
    if (nQ_*nso_*nso_ > dim) dim = nQ_*nso_*nso_;

    long int                  tempvdim = o*o*v*v;
    if (nQ_*o*v > tempvdim)    tempvdim = nQ_*o*v;
    if (nso_*nso_ > tempvdim) tempvdim = nso_*nso_;

    double * integrals = (double *)malloc(dim * sizeof(double));
    double * tempt     = (double *)malloc((o*(o+1)*v*(v+1)) * sizeof(double));
    double * tempv     = (double *)malloc(tempvdim * sizeof(double));
    double * Abij      = (double *)malloc(o*(o+1)/2*v * sizeof(double));
    double * Sbij      = (double *)malloc(o*(o+1)/2*v * sizeof(double));

    // construct c0
    //Normalization();

    // C2 = -1/2 t(bc,kj) (ki|ac)
    //      +    t(bc,ki) (kj|ac) 
    F_DGEMM('n','t',v*v,o*o,nQ_,1.0,Qvv_,v*v,Qoo_,o*o,0.0,tempv,v*v);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int i = 0; i < o; i++) {
            for (int k = 0; k < o; k++) {
                for (int c = 0; c < v; c++) {
                    integrals[a*o*o*v+i*o*v+k*v+c] = tempv[k*o*v*v+i*v*v+a*v+c];
                }
            }
        }
    }
    #pragma omp parallel for schedule (static)
    for (int b = 0; b < v; b++) {
        for (int j = 0; j < o; j++) {
            for (int k = 0; k < o; k++) {
                for (int c = 0; c < v; c++) {
                    tempt[b*o*o*v+j*o*v+k*v+c] = t2_[b*o*o*v+c*o*o+k*o+j];
                }
            }
        }
    }
    F_DGEMM('t','n',o*v,o*v,o*v,-1.0,integrals,o*v,tempt,o*v,0.0,tempv,o*v);

    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    r2_[a*o*o*v+b*o*o+i*o+j] = 0.5 * tempv[b*o*o*v+j*o*v+a*o+i] + tempv[b*o*o*v+i*o*v+a*o+j];
                }
            }
        }
    }

    // D2: 1/2 U(b,c,j,k) L(a,i,k,c)
    F_DGEMM('n','t',o*v,o*v,nQ_,2.0,Qov_,o*v,Qvo_,o*v,0.0,integrals,o*v);
    F_DGEMM('n','t',o*o,v*v,nQ_,-1.0,Qoo_,o*o,Qvv_,v*v,0.0,tempv,o*o);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int i = 0; i < o; i++) {
            for (int k = 0; k < o; k++) {
                for (int c = 0; c < v; c++) {
                    integrals[a*o*o*v+i*o*v+k*v+c] += tempv[a*o*o*v+c*o*o+k*o+i];
                }
            }
        }
    }
    #pragma omp parallel for schedule (static)
    for (int k = 0; k < o; k++) {
        for (int c = 0; c < v; c++) {
            for (int b = 0; b < v; b++) {
                for (int j = 0; j < o; j++) {
                    tempt[k*o*v*v+c*o*v+b*o+j] = 2.0 * t2_[b*o*o*v+c*o*o+j*o+k] - t2_[b*o*o*v+c*o*o+k*o+j];
                }
            }
        }
    }
    F_DGEMM('n','n',o*v,o*v,o*v,0.5,tempt,o*v,integrals,o*v,0.0,tempv,o*v);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    r2_[a*o*o*v+b*o*o+i*o+j] += tempv[a*o*o*v+i*o*v+b*o+j];
                }
            }
        }
    }

    // E2 a: t(ac,ij) F(bc)
    #pragma omp parallel for schedule (static)
    for (int c = 0; c < v; c++) {
        for (int a = 0; a < v; a++) {
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    tempt[c*o*o*v+a*o*o+i*o+j] = t2_[a*o*o*v+c*o*o+i*o+j];
                }
            }
        }
    }
    F_DGEMM('n','n',o*o*v,v,v,1.0,tempt,o*o*v,fvv_,v,0.0,tempv,o*o*v);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    r2_[a*o*o*v+b*o*o+i*o+j] += tempv[b*o*o*v+a*o*o+i*o+j];
                }
            }
        }
    }

    // E2 b: -t(a,b,i,k) F(kj)
    F_DGEMM('n','n',o,o*v*v,o,-1.0,foo_,o,t2_,o,1.0,r2_,o);

    // R2 = R2 + P(ia,jb) R2
    C_DCOPY(o*o*v*v,r2_,1,integrals,1);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    r2_[a*o*o*v+b*o*o+i*o+j] += integrals[b*o*o*v+a*o*o+j*o+i];
                }
            }
        }
    }

    // B2 = t(ab,kl) [ (ki|lj) + t(cd,ij) (kc|ld) ]
    F_DGEMM('n','t',o*o,o*o,nQ_,1.0,Qoo_,o*o,Qoo_,o*o,0.0,integrals,o*o);
    #pragma omp parallel for schedule (static)
    for (int k = 0; k < o; k++) {
        for (int i = 0; i < o; i++) {
            for (int l = 0; l < o; l++) {
                for (int j = 0; j < o; j++) {
                    tempt[k*o*o*o+l*o*o+i*o+j] = integrals[k*o*o*o+i*o*o+l*o+j];
                }
            }
        }
    }
    F_DGEMM('n','n',o*o,v*v,o*o,1.0,tempt,o*o,t2_,o*o,1.0,r2_,o*o);

    // (ac|bd) t2(cd,ij)
    long int oov = o * o * v;
    long int oo = o * o;
    long int otri = o * (o + 1) / 2;
    long int vtri = v * (v + 1) / 2;

    #pragma omp parallel for schedule(static)
    for (long int i = 0; i < o; i++) {
        for (long int j = i; j < o; j++) {
            long int ij = INDEX(i, j);
            for (long int a = 0; a < v; a++) {
                for (long int b = a; b < v; b++) {
                    tempt[INDEX(a, b) * otri + ij] =
                        (t2_[a * oov + b * oo + i * o + j] + t2_[b * oov + a * oo + i * o + j]);
                    tempt[INDEX(a, b) * otri + ij + vtri * otri] =
                        (t2_[a * oov + b * oo + i * o + j] - t2_[b * oov + a * oo + i * o + j]);
                }
                tempt[INDEX(a, a) * otri + ij] = t2_[a * oov + a * oo + i * o + j];
            }
        }
    }

    double *Vcdb = integrals;
    double *Vm = integrals + v * v * v;
    double *Vp = Vm;

    // qvv transpose
    #pragma omp parallel for schedule(static)
    for (int q = 0; q < nQ_; q++) {
        C_DCOPY(v * v, Qvv_ + q * v * v, 1, integrals + q, nQ_);
    }
    C_DCOPY(nQ_ * v * v, integrals, 1, Qvv_, 1);

    for (long int a = 0; a < v; a++) {
        int nb = v - a;
        F_DGEMM('t', 'n', v, v * nb, nQ_, 1.0, Qvv_ + a * v * nQ_, nQ_, Qvv_ + a * v * nQ_, nQ_, 0.0, Vcdb, v);

        #pragma omp parallel for schedule(static)
        for (long int b = a; b < v; b++) {
            long int cd = 0;
            long int ind1 = (b - a) * vtri;
            long int ind2 = (b - a) * v * v;
            long int v1, v2;
            for (long int c = 0; c < v; c++) {
                for (long int d = 0; d <= c; d++) {
                    Vp[ind1 + cd] = Vcdb[ind2 + d * v + c] + Vcdb[ind2 + c * v + d];
                    cd++;
                }
            }
        }
        F_DGEMM('n', 'n', otri, nb, vtri, 0.5, tempt, otri, Vp, vtri, 0.0, Abij, otri);
        #pragma omp parallel for schedule(static)
        for (long int b = a; b < v; b++) {
            long int cd = 0;
            long int ind1 = (b - a) * vtri;
            long int ind2 = (b - a) * v * v;
            long int v1, v2;
            for (long int c = 0; c < v; c++) {
                for (long int d = 0; d <= c; d++) {
                    Vm[ind1 + cd] = Vcdb[ind2 + d * v + c] - Vcdb[ind2 + c * v + d];
                    cd++;
                }
            }
        }
        F_DGEMM('n', 'n', otri, nb, vtri, 0.5, tempt + otri * vtri, otri, Vm, vtri, 0.0, Sbij, otri);

        // contribute to residual
        #pragma omp parallel for schedule(static)
        for (long int b = a; b < v; b++) {
            for (long int i = 0; i < o; i++) {
                for (long int j = 0; j < o; j++) {
                    int sg = (i > j) ? 1 : -1;
                    r2_[a * oo * v + b * oo + i * o + j] +=
                        Abij[(b - a) * otri + INDEX(i, j)] + sg * Sbij[(b - a) * otri + INDEX(i, j)];
                    if (a != b) {
                        r2_[b * oov + a * oo + i * o + j] +=
                            Abij[(b - a) * otri + INDEX(i, j)] - sg * Sbij[(b - a) * otri + INDEX(i, j)];
                    }
                }
            }
        }
    }

    // qvv un-transpose
    #pragma omp parallel for schedule (static)
    for (int q = 0; q < nQ_; q++) {
        C_DCOPY(v*v,Qvv_+q,nQ_,integrals+q*v*v,1);
    }
    C_DCOPY(nQ_*v*v,integrals,1,Qvv_,1);

    free(integrals);
    free(tempt);
    free(tempv);
    free(Abij);
    free(Sbij);

/*
    // last funny term for coupled-pair methods:
    if ( options_.get_str("P2RDM_TYPE") == "K" ) {

        throw PsiException("implement me",__FILE__,__LINE__);

    }else if ( options_.get_str("P2RDM_TYPE") == "CID" ) {

        throw PsiException("implement me",__FILE__,__LINE__);

    }else if ( options_.get_str("P2RDM_TYPE") == "CEPA(0)" ) {

        // nothing to add

    }
*/

}

}// end of namespaces
