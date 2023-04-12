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

#ifndef BLAS_H
#define BLAS_H
#ifdef __GNUC__
// Avoid tons of warnings with 'more than one instance of overloaded function'
#pragma GCC system_header
#endif

/**
 * fortran-ordered blas routines
 */

#include "qed_blas_mangle.h"

typedef long int myinteger;
typedef double doublereal;

namespace psi{ namespace fnocc{

/**
 * fortran-ordered dgemv
 */
void F_DGEMV(char trans,myinteger m,myinteger n,doublereal alpha,doublereal*A,myinteger lda,
            doublereal*X,myinteger incx,doublereal beta,doublereal*Y,myinteger incy);
/**
 * fortran-ordered dgemm
 */
void F_DGEMM(char transa,char transb, myinteger m, myinteger n, myinteger k,
            doublereal alpha,doublereal*A,myinteger lda,doublereal*B,myinteger ldb,
            doublereal beta,doublereal*C,myinteger ldc);

/**
 * name mangling for fortran-ordered dgemv
 */
extern "C" {
    void dgemv(char&trans,myinteger&m,myinteger&n,doublereal&alpha,doublereal*A,myinteger&lda,
            doublereal*X,myinteger&incx,doublereal&beta,doublereal*Y,myinteger&incy);
};
inline void DGEMV(char&trans,myinteger&m,myinteger&n,doublereal&alpha,doublereal*A,myinteger&lda,
            doublereal*X,myinteger&incx,doublereal&beta,doublereal*Y,myinteger&incy){
    dgemv(trans,m,n,alpha,A,lda,X,incx,beta,Y,incy);
}
/**
 * name mangling for fortran-ordered dgemm
 */
extern "C" {
    void dgemm(char&transa,char&transb,myinteger&m,myinteger&n,myinteger&k,
         doublereal&alpha,doublereal*A,myinteger&lda,doublereal*B,myinteger&ldb,
         doublereal&beta,doublereal*C,myinteger&ldc);
};
inline void DGEMM(char&transa,char&transb,myinteger&m,myinteger&n,myinteger&k,
         doublereal&alpha,doublereal*A,myinteger&lda,doublereal*B,myinteger&ldb,
         doublereal&beta,doublereal*C,myinteger&ldc)
{
    dgemm(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
};
/**
 * name mangling dcopy
 */
extern "C" {
    void dcopy(myinteger&n,doublereal*dx,myinteger&incx,doublereal*dy,
         myinteger&incy);
};
inline void DCOPY(myinteger&n,doublereal*dx,myinteger&incx,doublereal*dy,
            myinteger&incy){
    dcopy(n,dx,incx,dy,incy);
}
/**
 * name mangling dnrm2
 */
extern"C"{
    double dnrm2(myinteger&N,doublereal*X,myinteger&INCX);
};
inline double DNRM2(myinteger&N,doublereal*X,myinteger&INCX){
    return dnrm2(N,X,INCX);
};
/**
 * name mangling dgesv
 */
extern"C" {
    void dgesv(myinteger &N,myinteger &NRHS,doublereal*A,myinteger &LDA,myinteger*IPIV,doublereal*B,myinteger &LDB,myinteger &INFO);
};
inline void DGESV(myinteger &N,myinteger &NRHS,doublereal*A,myinteger &LDA,myinteger*IPIV,doublereal*B,myinteger &LDB,myinteger &INFO){
    dgesv(N,NRHS,A,LDA,IPIV,B,LDB,INFO);
};
/**
 * name mangling ddot
 */
extern "C" {
    double ddot(myinteger&n,doublereal*dx,myinteger&incx,doublereal*dy,myinteger&incy);
};
inline double DDOT(myinteger&n,doublereal*dx,myinteger&incx,doublereal*dy,myinteger&incy){
    return ddot(n,dx,incx,dy,incy);
}


/**
 * diagonalize a real symmetric matrix
 */
void Diagonalize(myinteger N,doublereal*A,doublereal*W);
/**
 * name mangling dsyev
 */
extern "C" {
    void dsyev(char&JOBZ,char&UPLO,myinteger&N,doublereal*A,myinteger&LDA,doublereal*W,doublereal*WORK,myinteger&LWORK,myinteger&INFO);
};
inline void DSYEV(char&JOBZ,char&UPLO,myinteger&N,doublereal*A,myinteger&LDA,doublereal*W,doublereal*WORK,myinteger&LWORK,myinteger&INFO){
    dsyev(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO);
}
/**
 * name mangling dspev
 */
extern "C" {
    void dspev(char&JOBZ,char&UPLO,myinteger&N,doublereal*AP,doublereal*W,doublereal*Z,myinteger&LDZ,doublereal*WORK,myinteger&INFO);
};
inline void DSPEV(char&JOBZ,char&UPLO,myinteger&N,doublereal*AP,doublereal*W,doublereal*Z,myinteger&LDZ,doublereal*WORK,myinteger&INFO){
    dspev(JOBZ,UPLO,N,AP,W,Z,LDZ,WORK,INFO);
}

/**
 *  General SVD
 */
void SVD(myinteger M,myinteger N,doublereal*A,doublereal*U,doublereal*VT,doublereal*S);
/**
 * name mangling dgesvd
 */
extern "C" {
    void dgesvd(char&JOBU,char&JOBVT,myinteger&M,myinteger&N,doublereal*A,myinteger&LDA,doublereal*S,doublereal*U,myinteger&LDU,doublereal*VT,myinteger&LDVT,doublereal*WORK,myinteger&LWORK,myinteger&INFO);
};
inline void DGESVD(char&JOBU,char&JOBVT,myinteger&M,myinteger&N,doublereal*A,myinteger&LDA,doublereal*S,doublereal*U,myinteger&LDU,doublereal*VT,myinteger&LDVT,doublereal*WORK,myinteger&LWORK,myinteger&INFO){
    dgesvd(JOBU,JOBVT,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,INFO);
}


/**
 * name manging dgeev
 */
extern "C" {
    void dgeev(char &jobvl,char &jobvr,myinteger &n,doublereal *a,myinteger &lda,
                    doublereal *wr,doublereal *wi, doublereal *vl,myinteger &ldvl,doublereal *vr,
                    myinteger &ldvr,doublereal * work,myinteger &lwork,myinteger &info);
};
inline void DGEEV(char &jobvl,char &jobvr,myinteger &n,doublereal*a,myinteger &lda,
                  doublereal*wr,doublereal*wi, doublereal*vl,myinteger&ldvl,doublereal*vr,
                  myinteger &ldvr,doublereal * work,myinteger &lwork,myinteger&info)
{
    dgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info);
}
/**
 *  * name mangling dggev
 *   */
extern "C" {
    void dgeqrf(myinteger &m,myinteger &n,doublereal *a,myinteger &lda,doublereal *tau,
                doublereal * work,myinteger &lwork,myinteger &info);
};
extern "C" {
    void dorgqr(myinteger &m,myinteger &n,myinteger &k,doublereal *a,myinteger &lda,doublereal *tau,
                doublereal * work,myinteger &lwork,myinteger &info);
};



/**
 * diagonalize a real symmetric packed matrix
 */
void Diagonalize2(myinteger N,doublereal*AP,doublereal*W,doublereal*Z);

/**
  * diagonalize general matrix and keep eigenvectors
  */
inline void Diagonalize3(double*M,long int dim, double*eigval,double*eigvec,double*leigvec, bool doSort = true, double* wi = nullptr){
    myinteger info;
    char vl = 'V';
    char vr = 'V';
    myinteger n = dim;
    myinteger lwork = 4*n;
    auto * work  = (doublereal*)malloc(lwork*sizeof(doublereal));

    bool keepwi = wi != nullptr;
    if(!keepwi) wi = (double *) malloc(n * sizeof(double));

    //double * eigvecl = (double*)malloc(n*n*sizeof(double));
    //memset((void*)eigvecl,'\0',n*n*sizeof(double));
    DGEEV(vl,vr,n,M,n,eigval,wi,leigvec,n,eigvec,n,work,lwork,info);

    auto * tempv = (double *) malloc(n * sizeof(double));
    double * tempvi;
    if(keepwi) tempvi = (double *) malloc(n * sizeof(double));
    auto * templ = (double*)malloc(n*n*sizeof(double));
    auto * tempr = (double*)malloc(n*n*sizeof(double));

    // sort eigenvalues and eigenvectors
    int count = 0;
    int*skip=(int*)malloc(n*sizeof(int));
    for (int i = 0; i < n; i++) skip[i] = 0;

    if( doSort ){
        for (int i = 0; i < n; i++) {
            double min = 1e99;
            int mini   = -1;
            for (int j = 0; j < n; j++) {
                if (skip[j]) continue;
                if (eigval[j] < min) {
                    min  = eigval[j];
                    mini = j;
                }
            }
            skip[mini] = 1;
            tempv[i] = min;
            if(keepwi) tempvi[i] = wi[mini];
            for (int j = 0; j < n; j++) {
                tempr[i*n+j] = eigvec[mini*n+j];
                templ[i*n+j] = leigvec[mini*n+j];
            }
        }
    }

    myinteger dx = 1;
    myinteger nn = n*n;
    if( doSort ){
        DCOPY(n,tempv,dx,eigval,dx);
        if(keepwi) DCOPY(n, tempvi, dx, wi, dx);
        DCOPY(nn,tempr,dx,eigvec,dx);
        DCOPY(nn,templ,dx,leigvec,dx);
    }
    //free(eigvecl);
    free(templ);
    free(tempr);
    free(tempv);
    free(keepwi ? tempvi : wi);
    free(skip);
    free(work);
}

}}

#endif
