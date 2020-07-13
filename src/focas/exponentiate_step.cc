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

#include"v2rdm_solver.h"

#include <misc/blas.h>
#include <misc/omp.h>

#define TINY 1.0e-16

using namespace psi;
using namespace fnocc;

namespace psi{ namespace v2rdm_casscf{

void v2RDMSolver::exponentiate_step(double * X) { 

    int maxnmopi = 0;
    for (int h = 0; h < nirrep_; h++) {
        int mynmopi = rstcpi_[h] + amopi_[h] + rstvpi_[h];
        if ( mynmopi > maxnmopi ) maxnmopi = mynmopi;
    }

    double * mat  = (double *)malloc(maxnmopi * maxnmopi * sizeof(double));
    double * mat2 = (double *)malloc(maxnmopi * maxnmopi * sizeof(double));
    double * cosK = (double *)malloc(maxnmopi * maxnmopi * sizeof(double));
    double * sinK = (double *)malloc(maxnmopi * maxnmopi * sizeof(double));
    double * evec = (double *)malloc(maxnmopi * maxnmopi * sizeof(double));
    double * eval = (double *)malloc(maxnmopi * sizeof(double));

    double * expK = (double *)malloc(nmo_ * nmo_ * sizeof(double));
    memset((void*)expK,'\0',nmo_ * nmo_ * sizeof(double));

    int lwork     = (1 + 6*maxnmopi + 2*maxnmopi * maxnmopi);
    double * work = (double*)malloc(lwork * sizeof(double));

    int liwork    = 3 + 5 * maxnmopi;
    int * iwork   = (int*)malloc(liwork * sizeof(int));

    // loop over target irrep
    for (int h_want = 0; h_want < nirrep_; h_want++){

        if ( nmopi_[h_want] == 0 ) continue;

        //memset((void*)mat,'\0',nmopi_[h_want] * nmopi_[h_want] * sizeof(double));
        memset((void*)mat,'\0',maxnmopi * maxnmopi * sizeof(double));

        int ind = 0;
        for (int h = 0; h < nirrep_; h++){

            for (int i = 0; i < rstcpi_[h]; i++){

               int j_off = rstcpi_[h_want];

               // d-a pairs
               for (int j = 0; j < amopi_[h]; j++){

                   if ( h_want == h ) {

                       mat[(i      )*nmopi_[h] + (j+j_off)] =  X[ind];
                       mat[(j+j_off)*nmopi_[h] + (i      )] = -X[ind];

                   }

                   ind++;

               }

               j_off += amopi_[h_want];

               // d-e pairs
               for (int j =0; j < rstvpi_[h]; j++){

                   if ( h_want == h ) {

                       mat[(i      )*nmopi_[h] + (j+j_off)] =  X[ind];
                       mat[(j+j_off)*nmopi_[h] + (i      )] = -X[ind];
 
                   }

                   ind++;

               }

            }

        }

        for (int h = 0; h < nirrep_; h++){

            int i_off = rstcpi_[h_want];

            for (int i = 0; i < amopi_[h]; i++){

                int j_off = rstcpi_[h_want];

                if ( options_.get_bool("ORBOPT_ACTIVE_ACTIVE_ROTATIONS") ) {

                    // a-a pairs

                    for (int j = i + 1; j < amopi_[h]; j++){

                        if ( h_want == h ){

                            mat[(i+i_off)*nmopi_[h] + (j+j_off)] =  X[ind];
                            mat[(j+j_off)*nmopi_[h] + (i+i_off)] = -X[ind];

                        }

                        ind++;

                    }

                }

                j_off += amopi_[h_want];

                // a-e pairs

                for (int j = 0; j<rstvpi_[h]; j++){

                    if ( h_want == h ){

                        mat[(i+i_off)*nmopi_[h] + (j+j_off)] =  X[ind];
                        mat[(j+j_off)*nmopi_[h] + (i+i_off)] = -X[ind];

                    }

                    ind++;

                }

            }

        }
    
        // square of step matrix
        F_DGEMM('n','n',nmopi_[h_want],nmopi_[h_want],nmopi_[h_want],1.0,mat,nmopi_[h_want],mat,nmopi_[h_want],0.0,evec,nmopi_[h_want]);

        // diagonalize square of step matrix
        C_DSYEVD('V', 'U', nmopi_[h_want], evec, nmopi_[h_want], eval, work, lwork, iwork, liwork);

        // d = sqrt(-lambda)
        for (int i = 0; i < nmopi_[h_want]; i++) {
            if ( eval[i] > 0.0 ) {
                eval[i] = 0.0; // K*K should be negative definite
            }else {
                eval[i] = sqrt(-eval[i]);
            }
        }

        // U = exp(K) = X cos(d) X^T + K X d^-1 sin(d) X^T

        // cos(d) X^T
        for (int p = 0; p < nmopi_[h_want]; p++) {
            for (int j = 0; j < nmopi_[h_want]; j++) {
                mat2[p * nmopi_[h_want] + j] = evec[p * nmopi_[h_want] + j] * cos(eval[p]);
            }
        }

        // X cos(d) X^T
        F_DGEMM('n', 't', nmopi_[h_want], nmopi_[h_want], nmopi_[h_want], 1.0, mat2, nmopi_[h_want], evec, nmopi_[h_want], 0.0, cosK, nmopi_[h_want]);

        // d^-1 sin(d) X^T
        for (int p = 0; p < nmopi_[h_want]; p++) {
            for (int j = 0; j < nmopi_[h_want]; j++) {
                if ( eval[p] != 0.0 ) {
                    mat2[p * nmopi_[h_want] + j] = evec[p * nmopi_[h_want] + j] * sin(eval[p])/eval[p];
                }else {
                    mat2[p * nmopi_[h_want] + j] = evec[p * nmopi_[h_want] + j];
                }
            }
        }

        // X d^-1 sin(d) X^T
        F_DGEMM('n', 't', nmopi_[h_want], nmopi_[h_want], nmopi_[h_want], 1.0, mat2, nmopi_[h_want], evec, nmopi_[h_want], 0.0, sinK, nmopi_[h_want]);

        // K X d^-1 sin(d) X^T
        F_DGEMM('t', 'n', nmopi_[h_want], nmopi_[h_want], nmopi_[h_want], 1.0, sinK, nmopi_[h_want], mat, nmopi_[h_want], 1.0, cosK, nmopi_[h_want]);

        //printf("C++ U for irrep %5i nmo: %5i\n",h_want+1,nmopi_[h_want]);
        //for (int i = 0; i < nmopi_[h_want]; i++) {
        //    for (int j = 0; j <nmopi_[h_want]; j++) {
        //        printf("%5i %5i %20.12le\n",i+1,j+1,cosK[i*nmopi_[h_want]+j]);
        //        fflush(stdout);
        //    }
        //}

        // now, put this block of exponential of K should be in expK (actually, we need to transpose it ...)
        int expK_off = 0;
        for (int h = 0; h < h_want; h++) {
            expK_off += nmopi_[h]*nmopi_[h];
        }
        for (int i = 0; i < nmopi_[h_want]; i++) {
            for (int j = 0; j < nmopi_[h_want]; j++) {
                expK[ expK_off + j * nmopi_[h_want] + i ] = cosK[i * nmopi_[h_want] + j];
            }
        }
        //C_DCOPY(nmopi_[h_want]*nmopi_[h_want],cosK,1,expK + expK_off,1);

    }

    // copy exp(X) into X
    C_DCOPY(nmo_*nmo_,expK,1,X,1);

    free(work);
    free(iwork);

    free(mat);
    free(mat2);

    free(cosK);
    free(sinK);
    free(expK);

    free(evec);
    free(eval);

}


}} // end of namespaces
