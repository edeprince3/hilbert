#include "gpu_transform_3index_teint.h"

void transform_3index_teints_driver_(double *int2_, double *U_, int *nmopi_, int* nirrep_, int* nQ_){

  int* U_offset;
  int* nmo_offset;
 
  double *int2_tmp1;
  double *int2_tmp2;

  int nirrep = *nirrep_;
  int nQ     = *nQ_;

  int nmo_tot   = 0;
  int max_nmopi = 0;

  for (int sym_R = 0; sym_R < nirrep; sym_R++){

    nmo_tot += nmopi_[sym_R];

    if ( max_nmopi > nmopi_[sym_R] ) continue;

    max_nmopi = nmopi_[sym_R];

  }

  double max_gpu_mem_mb = 1.0e3;

  int ngem_tot_lt = nmo_tot * ( nmo_tot + 1 ) / 2;

  int max_tile_size=max_numQ_per_pass(nQ, nmo_tot, max_nmopi, nirrep, max_gpu_mem_mb, nmopi_);
 
  int num_tiles= nQ / max_tile_size;

  // allocate and figure out offset matrix for indexing elements of U

  U_offset      = (int*) malloc(nirrep*sizeof(int));

  nmo_offset    = (int*) malloc(nirrep*sizeof(int));

  U_offset[0]   = 0;

  nmo_offset[0] = 0;

  for ( int sym_L = 1; sym_L < nirrep; sym_L++ ){

    U_offset[sym_L]   = U_offset[sym_L - 1] + nmopi_[sym_L - 1] * nmopi_[sym_L - 1];

    nmo_offset[sym_L] = nmo_offset[sym_L - 1] + nmopi_[sym_L - 1];

  }

  int2_tmp1     = (double*)malloc(max_tile_size * max_nmopi * max_nmopi * sizeof(double));

  int2_tmp2     = (double*)malloc(max_tile_size * max_nmopi * max_nmopi * sizeof(double));


  for ( int tile = 0 ; tile < num_tiles; tile++){

    // first Q index is tile*max_tile_size, and max_tile_size Qs are copied

    transform_3index_teints_block( &int2_[tile*max_tile_size*ngem_tot_lt], int2_tmp1, int2_tmp2, U_, nmopi_, U_offset, nmo_offset, nirrep, max_tile_size, max_nmopi, nmo_tot, ngem_tot_lt);
  }

  if ( nQ > num_tiles * max_tile_size){

    int tile_size = nQ  - num_tiles*max_tile_size;

    // first Q index is num_tiles*max_tile_size and tile_size Qs are copied

    transform_3index_teints_block( &int2_[num_tiles*max_tile_size*ngem_tot_lt], int2_tmp1, int2_tmp2, U_, nmopi_, U_offset, nmo_offset, nirrep, tile_size, max_nmopi, nmo_tot, ngem_tot_lt);

  }

  free(int2_tmp1);

  free(int2_tmp2);

  free(U_offset);

  free(nmo_offset);

}

void transform_3index_teints_block(double *int2_, double *int2_tmp1, double *int2_tmp2, double *U_, int *nmopi_, int* U_offset, int* nmo_offset, int nirrep, int nQ, int max_nmopi, int nmo_tot, int ngem_tot_lt){

  long int int_ind_off, int_ind, int_tmp_off;

  // transform integrals

  for ( int Q = 0; Q < nQ; Q++){

    int_tmp_off = (long int) Q * (long int) max_nmopi * (long int) max_nmopi;

    int_ind_off = (long int) ngem_tot_lt * (long int) Q;  
  
    // unpack integrals for this Q  (L > R)

    for (int sym_L = 0; sym_L < nirrep; sym_L++){

      // *************
      // SYM_R < SYM_L
      // *************

      for (int sym_R = 0; sym_R < sym_L; sym_R++){

        // UNPACK 

        for (int L = 0; L < nmopi_[sym_L]; L++){

          int_ind = int_ind_off + (long int) gamma_mn(L+nmo_offset[sym_L],nmo_offset[sym_R]);

          for (int R = 0; R < nmopi_[sym_R]; R++){

            int LR = (long int) ind_ij( L, R, nmopi_[sym_L] ) + int_tmp_off;

            int2_tmp1[LR] = int2_[int_ind];

            int_ind++;

          }

        }

        // TRANSFORM R

        for ( int L = 0; L < nmopi_[sym_L]; L++){

          for ( int R = 0; R < nmopi_[sym_R]; R++){

            double tmp = 0.0;

            for ( int K = 0; K < nmopi_[sym_R]; K++){

              long int LK = (long int) ind_ij( L, K, nmopi_[sym_L] ) + int_tmp_off;

              int KR = ind_ij( K, R, nmopi_[sym_R] ) + U_offset[sym_R];

              tmp += int2_tmp1[LK] * U_[KR];

            }

            long int LR = (long int) ind_ij( L, R, nmopi_[sym_L] ) + int_tmp_off;

            int2_tmp2[LR] = tmp;

          }

        }

        // TRANSFORM L

        for ( int R = 0; R < nmopi_[sym_R]; R++){

          for ( int L = 0; L < nmopi_[sym_L]; L++){

            double tmp = 0.0;

            for ( int K = 0; K < nmopi_[sym_L]; K++){

              int KL = ind_ij(K, L, nmopi_[sym_L]) + U_offset[sym_L];

              long int KR = (long int) ind_ij(K, R, nmopi_[sym_L]) + int_tmp_off;

              tmp += U_[KL] * int2_tmp2[KR];

            }

            long int LR = (long int) ind_ij(L, R, nmopi_[sym_L]) + int_tmp_off;

            int2_tmp1[LR] = tmp;

          }

        }

        // REPACK

        for (int L = 0; L < nmopi_[sym_L]; L++){

          int_ind = int_ind_off + (long int) gamma_mn(L+nmo_offset[sym_L],nmo_offset[sym_R]);

          for (int R = 0; R < nmopi_[sym_R]; R++){

            long int LR = (long int) ind_ij( L, R, nmopi_[sym_L] ) + int_tmp_off;

            int2_[int_ind] = int2_tmp1[LR];

            int_ind++;

          }

        }

      }

      // **************
      // SYM_R == SYM_L
      // **************

      // UNPACK

      for (int L = 0; L < nmopi_[sym_L]; L++){

        int_ind = int_ind_off + (long int) gamma_mn(L+nmo_offset[sym_L],nmo_offset[sym_L]);

        for (int R = 0; R < L; R++){

          long int LR = (long int) ind_ij( L, R, nmopi_[sym_L] ) + int_tmp_off;

          int2_tmp1[LR] = int2_[int_ind];

          LR = (long int) ind_ij( R, L, nmopi_[sym_L] ) + int_tmp_off;

          int2_tmp1[LR] = int2_[int_ind];
 
          int_ind++;

        }

        long int LL = (long int) ind_ij( L, L, nmopi_[sym_L] ) + int_tmp_off;

        int2_tmp1[LL] = int2_[int_ind];

      }

      // TRANSFORM R

      for ( int L = 0; L < nmopi_[sym_L]; L++){

        for ( int R = 0; R < nmopi_[sym_L]; R++){

          double tmp = 0.0;

          for ( int K = 0; K < nmopi_[sym_L]; K++){

            long int LK = (long int) ind_ij( L, K, nmopi_[sym_L] ) + int_tmp_off;

            int KR = ind_ij( K, R, nmopi_[sym_L] ) + U_offset[sym_L];

            tmp += int2_tmp1[LK] * U_[KR];

          }

          long int LR = (long int) ind_ij( L, R, nmopi_[sym_L] ) + int_tmp_off;

          int2_tmp2[LR] = tmp;

        }

      }

      // TRANSFORM L

      for ( int R = 0; R < nmopi_[sym_L]; R++){

        for ( int L = 0; L < nmopi_[sym_L]; L++){

          double tmp = 0.0;

          for ( int K = 0; K < nmopi_[sym_L]; K++){

            int KL = ind_ij(K, L, nmopi_[sym_L]) + U_offset[sym_L];

            long int KR = (long int) ind_ij(K, R, nmopi_[sym_L]) + int_tmp_off;

            tmp += U_[KL] * int2_tmp2[KR];

          }

          long int LR = (long int) ind_ij(L, R, nmopi_[sym_L]) + int_tmp_off;

          int2_tmp1[LR] = tmp;

        }

      }

      // UNPACK

      for (int L = 0; L < nmopi_[sym_L]; L++){

        int_ind = int_ind_off + (long int) gamma_mn(L+nmo_offset[sym_L],nmo_offset[sym_L]);

        for (int R = 0; R < L; R++){

          long int LR = (long int) ind_ij( L, R, nmopi_[sym_L] ) + int_tmp_off;

          int2_[int_ind] = int2_tmp1[LR];

          int_ind++;

        }

        long int LL = (long int) ind_ij( L, L, nmopi_[sym_L] ) + int_tmp_off;

        int2_[int_ind] = int2_tmp1[LL];

      }

    }

  }
 
}

int gamma_mn(int gamma_m, int gamma_n){

  if ( gamma_m >= gamma_n ){

    return gamma_m * ( gamma_m + 1 ) / 2 + gamma_n;

  }

  return gamma_n * ( gamma_n + 1 ) / 2 + gamma_m;

}

int ind_ij(int i, int j, int dim_i){

  return j * dim_i + i;

}

int max_numQ_per_pass(int nQ, int nmo_tot, int nirrep, int max_nmopi, double max_gpu_mem_mb, int *nmopi){

  //number of doubles for a single instance of: 

  //   intermediate matrices transformation (int2_tmp1 and int2_tmp2)
  //   packed integrals

  double mem_total  = 2.0 * (double) max_nmopi * (double) max_nmopi;
  mem_total        += (double) nmo_tot * ( (double) nmo_tot + 1.0 ) / 2.0;

  // this is the total memory (in mb) for auxiliary arrays for a single GPU block
  mem_total *= (double) sizeof(double) / ( 1024.0 * 1024.0 );
 
  int num_Q = max_gpu_mem_mb / mem_total;

  if ( num_Q < 1 ) num_Q = 1;

  if ( num_Q > nQ ) num_Q = nQ;
  
  return num_Q;

}
