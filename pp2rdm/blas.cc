#ifndef BLAS_EIG_H
#define BLAS_EIG_H

#include"fortran.h"
#include "blas.h"

#include <psi4/libqt/qt.h>

#include<math.h>

#include <psi4/libpsi4util/process.h>
#include <psi4/libpsi4util/PsiOutStream.h>

#include<psi4/libmints/vector.h>
#include<psi4/libmints/matrix.h>

typedef long int integer;
typedef double doublereal;


namespace psi{ namespace pp2rdm{

/**
  * diagonalize general matrix and keep eigenvectors
  */
void NonSymmetricEigenvalueEigenvector(long int dim, double * M, double * eigval, double * el, double * er) {
  long int info;
  char vl = 'V';
  char vr = 'V';
  long int n = dim;
  long int lwork = 4*n;
  double * work  = (doublereal*)malloc(lwork*sizeof(doublereal));
  double * wi = (double*)malloc(n*sizeof(double));

 // DGEEV(&vl,&vr,&n,M,&n,eigval,wi,el,&n,er,&n,work,&lwork,&info);
  DGEEV(vl,vr,n,M,n,eigval,wi,el,n,er,n,work,lwork,info);

  // sort eigenvalues and eigenvectors
  //int count = 0;
  //int*skip=(int*)malloc(n*sizeof(int));
  //for (int i = 0; i < n; i++) skip[i] = 0;
  //for (int i = 0; i < n; i++) {
  //    double min = 1e99;
  //    int mini   = -1;
  //    for (int j = 0; j < n; j++) {
  //        if (skip[j]) continue;
  //        if (eigval[j] < min) {
  //            min  = eigval[j];
  //            mini = j;
  //        }
  //    }
  //    skip[mini] = 1;
  //    wi[i] = min;
  //}
  //C_DCOPY(n,wi,1,eigval,1);

  //free(skip);

  free(wi);
  free(work);
}

/**
  * diagonalize general matrix, don't compute eigenvectors
  */
void GeneralizedEigenvalueProblem(long int N, double * A, double * B, double * c, double * eig, bool * energy_is_real){
    char JOBVR = 'V';
    char JOBVL = 'V';

    long int LDA  = N;
    long int LDB  = N;
    long int LDVR = N;
    long int LDVL = N;
    long int LWORK = 8*N;

    double*WORK=(double*)malloc(LWORK*sizeof(double));
    double*ALPHAR=(double*)malloc(N*sizeof(double));
    double*ALPHAI=(double*)malloc(N*sizeof(double));
    double*BETA=(double*)malloc(N*sizeof(double));
    double*VR=(double*)malloc(N*N*sizeof(double));
    double*VL=(double*)malloc(N*N*sizeof(double));
    std::shared_ptr<Matrix> Al (new Matrix(N,N));
//    outfile->Printf("\n");
//    outfile->Printf("Before solving eq, the matrix A is\n");
//    for(int i = 0; i<N;i++){
//       for(int j = 0; j<N;j++){
//       Al->pointer()[i][j] = A[i*N+j];
//       }
//    }
//    Al->print();
//     outfile->Printf("\n");
//    std::shared_ptr<Matrix> Bl (new Matrix(N,N));
//    outfile->Printf("\n");
//    outfile->Printf("Before solving eq, the matrix B is\n");
//    for(int i = 0; i<N;i++){
//       for(int j = 0; j<N;j++){
//       Bl->pointer()[i][j] = B[i*N+j];
//       }
//    }
//    Bl->print();
//     outfile->Printf("\n");
    integer INFO=0;
    INFO = C_DGGEV(JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK);
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
       // eig[i] = ALPHAR[i] / BETA[i];
       // alphai->pointer()[i]=ALPHAI[i];
      //  beta->pointer()[i]=BETA[i];
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
    C_DCOPY(N*N,VR,1,B,1);
    C_DCOPY(N*N,VL,1,A,1);
    outfile->Printf("\n");
    outfile->Printf("After solving eigenvalue problem, the matrix A is\n");
    std::shared_ptr<Matrix> An (new Matrix(N,N));
    for(int i = 0; i<N;i++){
       for(int j = 0; j<N;j++){
       An->pointer()[i][j] = A[i*N+j];
       }
    }
    An->print();

    outfile->Printf("\n");
    outfile->Printf("After solving eigenvalue problem, the matrix B is\n");
    std::shared_ptr<Matrix> Bn (new Matrix(N,N));
    for(int i = 0; i<N;i++){
       for(int j = 0; j<N;j++){
       Bn->pointer()[i][j] = B[i*N+j];
       }
    }
    Bn->print();

    free(VR);
    free(VL);
    free(WORK);
    free(ALPHAR);
    free(ALPHAI);
    free(BETA);
}

/**
 * name mangling dggev
 */
extern "C" {
    void dsygv(long int &ITYPE,char&JOBZ,char&UPLO, long int &N, double*A,long int&LDA,double*B,long int &LDB,
        double *W, double*WORK,long int&LWORK,long int&INFO);
};
inline void DSYGV(long int &ITYPE,char&JOBZ,char&UPLO, long int &N, double*A,long int&LDA,double*B,long int &LDB,
        double *W, double*WORK,long int&LWORK,long int&INFO){
    dsygv(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,INFO);
}

int SymmetricGeneralizedEigenvalueProblem(long int N, double * A, double * B, double * c, double * eig){
    long int ITYPE = 1; // Ax = l Bx (actually solve Bx = 1/l Ax)
    char JOBZ = 'V';
    char UPLO = 'U';
    long int LDA  = N;
    long int LDB  = N;

    long int LWORK = 3*N-1;
    double*WORK=(double*)malloc(LWORK*sizeof(double));

    integer INFO=0;
    C_DSYGV(ITYPE,JOBZ,UPLO,N,B,LDA,A,LDB,eig,WORK,LWORK);
    return (int)INFO;

    free(WORK);
}

}} // end of namespace

#endif
