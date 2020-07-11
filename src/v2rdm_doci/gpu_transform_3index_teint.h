#include<stdio.h>
#include<stdlib.h>

extern "C" {
  void transform_3index_teints_driver_(double *int2_, double *U_, int *nmopi_, int* nirrep_, int* nQ_);
}

void transform_3index_teints_block(double *int2_, double *int2_tmp1, double *int2_tmp2, double *U_, int *nmopi_, int *U_offset, int *nmo_offset, int nirrep_, int nQ_, int max_nmopi, int nmo_tot, int ngem_tot_lt);

int gamma_mn(int gamma_m, int gamma_n);

int ind_ij(int i, int j, int dim_i);

int max_numQ_per_pass(int nQ, int nmo_tot, int max_nmopi, int nirrep, double max_gpu_mem_mb, int *nmopi);
