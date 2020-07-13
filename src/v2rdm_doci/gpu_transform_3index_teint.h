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

#include <stdio.h>
#include <stdlib.h>

extern "C" {
  void transform_3index_teints_driver_(double *int2_, double *U_, int *nmopi_, int* nirrep_, int* nQ_);
}

void transform_3index_teints_block(double *int2_, double *int2_tmp1, double *int2_tmp2, double *U_, int *nmopi_, int *U_offset, int *nmo_offset, int nirrep_, int nQ_, int max_nmopi, int nmo_tot, int ngem_tot_lt);

int gamma_mn(int gamma_m, int gamma_n);

int ind_ij(int i, int j, int dim_i);

int max_numQ_per_pass(int nQ, int nmo_tot, int max_nmopi, int nirrep, double max_gpu_mem_mb, int *nmopi);
