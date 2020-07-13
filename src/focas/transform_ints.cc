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

#include <psi4/libtrans/integraltransform.h>

#include "transform_ints.h"

#include <misc/blas.h>

using namespace psi;
using namespace fnocc;

void transform_ints_driver(double * int1, double * int2, double * U, int *nmopi, int *frzvpi, int nirrep, int nQ){

  int max_num_threads = 1;

  int * U_offset;
  int * nmo_offset;

  // figure out total number of MOs, maximum number of MOs in an irrep as well as the nnumber of LT geminals

  int nmo_tot   = 0;
  int max_nmopi = 0;

  for (int sym_R = 0; sym_R < nirrep; sym_R++){

      nmo_tot += nmopi[sym_R];
      if ( max_nmopi > nmopi[sym_R] ) continue;
      max_nmopi = nmopi[sym_R];

  }

  int ngem_tot_lt = nmo_tot * ( nmo_tot + 1 ) / 2;

  // allocate memory on CPU

  U_offset      = (int*) malloc(nirrep*sizeof(int));
  nmo_offset    = (int*) malloc(nirrep*sizeof(int));

  // figure out offset arrays for U and MO indeces

  U_offset[0]   = 0;
  nmo_offset[0] = 0;

  for ( int sym_L = 1; sym_L < nirrep; sym_L++ ){

      U_offset[sym_L]   = U_offset[sym_L - 1] + nmopi[sym_L - 1] * nmopi[sym_L - 1];
      nmo_offset[sym_L] = nmo_offset[sym_L - 1] + nmopi[sym_L - 1];

  }

  // 3-index transformation
  transform_3index_tei( int2, U, nmopi, U_offset, nmo_offset, nirrep, nQ, 
                        max_nmopi, nmo_tot, ngem_tot_lt, max_num_threads);

  // 2-index transformation
  transform_oei( int1, U, nmopi, frzvpi, U_offset, nmo_offset, nirrep,
                 max_nmopi, nmo_tot, max_num_threads);

  // deallocate memory on CPU

  free(U_offset);
  free(nmo_offset);

}

void transform_3index_tei(double *int2, double *U, int *nmopi, int* U_offset, 
                                   int* nmo_offset, int nirrep, int nQ, int max_nmopi, 
                                   int nmo_tot, int ngem_tot_lt, int max_num_threads){

  long int int_ind_off, int_ind, int_tmp_off;

  double* int2_tmp1;
  double* int2_tmp2;

  int2_tmp1 = (double*)malloc(max_num_threads*max_nmopi*max_nmopi*sizeof(double));
  int2_tmp2 = (double*)malloc(max_num_threads*max_nmopi*max_nmopi*sizeof(double));

  // transform integrals

  for ( int Q = 0; Q < nQ; Q++){

      int_tmp_off = 0; //(long int) Q * (long int) max_nmopi * (long int) max_nmopi;

      int_ind_off = (long int) ngem_tot_lt * (long int) Q;  
  
      // unpack integrals for this Q  (L > R)

      for (int sym_L = 0; sym_L < nirrep; sym_L++){

          // *************
          // SYM_R < SYM_L
          // *************

          for (int sym_R = 0; sym_R < sym_L; sym_R++){

              // UNPACK 

              for (int L = 0; L < nmopi[sym_L]; L++){

                  int_ind = int_ind_off + (long int) INDEX(L+nmo_offset[sym_L],nmo_offset[sym_R]);

                  for (int R = 0; R < nmopi[sym_R]; R++){

                      int2_tmp1[L*nmopi[sym_R] + R] = int2[int_ind];
                      int_ind++;
  
                  }

              }

              // TRANSFORM

              // note: everything is backwards, greg. :(

              // I''(ij) = U(pi).I(pq).U(qj)
              //         = U(pi).I'(pj)

              // I'(pj) = I(pq).U(qj) = U(qj).I(pq)
              F_DGEMM('n','n',nmopi[sym_R],nmopi[sym_L],nmopi[sym_R],1.0,U + U_offset[sym_R],nmopi[sym_R],int2_tmp1,nmopi[sym_R],0.0,int2_tmp2,nmopi[sym_R]);
              // or is it
              // I'(p*sym[R] + j) = I(q*sym[L] + p).U(j*sym[R] + q) = U(j*sym[R] + q).I(q*sym[L] + p)
              //F_DGEMM('t','t',nmopi[sym_R],nmopi[sym_L],nmopi[sym_R],1.0,U + U_offset[sym_R],nmopi[sym_R],int2_tmp1,nmopi[sym_L],0.0,int2_tmp2,nmopi[sym_R]);

              // I''(ij) = U(pi).I(pq).U(qj)
              //         = I'(pj).U(pi)
              F_DGEMM('n','t',nmopi[sym_R],nmopi[sym_L],nmopi[sym_L],1.0,int2_tmp2,nmopi[sym_R],U + U_offset[sym_L],nmopi[sym_L],0.0,int2_tmp1,nmopi[sym_R]);
              // or is it
              //         = I'(pj).U(ip)
              //F_DGEMM('n','n',nmopi[sym_R],nmopi[sym_L],nmopi[sym_L],1.0,int2_tmp2,nmopi[sym_R],U + U_offset[sym_L],nmopi[sym_L],0.0,int2_tmp1,nmopi[sym_R]);

              // REPACK

              for (int L = 0; L < nmopi[sym_L]; L++){

                  int_ind = int_ind_off + (long int) INDEX(L+nmo_offset[sym_L],nmo_offset[sym_R]);

                  for (int R = 0; R < nmopi[sym_R]; R++){
  
                      int2[int_ind] = int2_tmp1[L*nmopi[sym_R] + R];
                      int_ind++;

                  }

              }
 
          }

          // **************
          // SYM_R == SYM_L
          // **************

          // UNPACK

          for (int L = 0; L < nmopi[sym_L]; L++){

              int_ind = int_ind_off + (long int) INDEX(L+nmo_offset[sym_L],nmo_offset[sym_L]);

              for (int R = 0; R < L; R++){

                  int2_tmp1[L*nmopi[sym_L] + R] = int2[int_ind];
                  int2_tmp1[R*nmopi[sym_L] + L] = int2[int_ind];

                  int_ind++;
  
              }

              int2_tmp1[L*nmopi[sym_L] + L] = int2[int_ind];

          }

          // TRANSFORM

          // I'(pj) = I(pq).U(qj) = U(qj).I(pq)
          F_DGEMM('n','n',nmopi[sym_L],nmopi[sym_L],nmopi[sym_L],1.0,U + U_offset[sym_L],nmopi[sym_L],int2_tmp1,nmopi[sym_L],0.0,int2_tmp2,nmopi[sym_L]);

          // I''(ij) = U(pi).I(pq).U(qj)
          //         = I'(pj).U(pi)
          F_DGEMM('n','t',nmopi[sym_L],nmopi[sym_L],nmopi[sym_L],1.0,int2_tmp2,nmopi[sym_L],U + U_offset[sym_L],nmopi[sym_L],0.0,int2_tmp1,nmopi[sym_L]);

          // UNPACK

          for (int L = 0; L < nmopi[sym_L]; L++){

              int_ind = int_ind_off + (long int) INDEX(L+nmo_offset[sym_L],nmo_offset[sym_L]);

              for (int R = 0; R < L; R++){

                  int2[int_ind] = int2_tmp1[L*nmopi[sym_L] + R];
                  int_ind++;

              }

              int2[int_ind] = int2_tmp1[L * nmopi[sym_L] + L];

          }

      }

  }
 
  free(int2_tmp1);
  free(int2_tmp2);

}

void transform_oei(double *int1, double *U, int *nmopi, int *frzvpi, int* U_offset, 
                      int* nmo_offset, int nirrep, int max_nmopi, 
                      int nmo_tot, int max_num_threads){

  double * int1_tmp1 = (double*)malloc(max_num_threads*max_nmopi*max_nmopi*sizeof(double));
  double * int1_tmp2 = (double*)malloc(max_num_threads*max_nmopi*max_nmopi*sizeof(double));

  // transform integrals

  for (int h = 0; h < nirrep; h++){

      // UNPACK
      int dim = nmopi[h] - frzvpi[h];

      for (int L = 0; L < dim; L++){

          int myoff = 0;
          for (int h2 = 0; h2 < h; h2++) {
              myoff +=  (nmopi[h2] - frzvpi[h2])* (nmopi[h2] - frzvpi[h2] + 1) / 2;
          }

          for (int R = 0; R < L; R++){
              int1_tmp1[L*dim + R] = int1[myoff + INDEX(L,R)];
              int1_tmp1[R*dim + L] = int1[myoff + INDEX(L,R)];
          }
          int1_tmp1[L*dim + L] = int1[myoff + INDEX(L,L)];

      }

      // TRANSFORM

      // I'(pj) = I(pq).U(qj) = U(qj).I(pq)
      F_DGEMM('n','n',dim,dim,dim,1.0,U + U_offset[h],dim,int1_tmp1,dim,0.0,int1_tmp2,dim);

      // I''(ij) = U(pi).I(pq).U(qj)
      //         = I'(pj).U(pi)
      F_DGEMM('n','t',dim,dim,dim,1.0,int1_tmp2,dim,U + U_offset[h],dim,0.0,int1_tmp1,dim);

      // REPACK

      for (int L = 0; L < dim; L++){

          int myoff = 0;
          for (int h2 = 0; h2 < h; h2++) {
              myoff +=  (nmopi[h2] - frzvpi[h2])* (nmopi[h2] - frzvpi[h2] + 1) / 2;
          }   

          for (int R = 0; R < L; R++){
              int1[myoff + INDEX(L,R)] = int1_tmp1[L*dim + R];
          }   
          int1[myoff + INDEX(L,L)] = int1_tmp1[L*dim + L];

      }

  }
 
  free(int1_tmp1);
  free(int1_tmp2);

}
