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

#ifndef transform_ints_h
#define transform_ints_h

#include<stdlib.h>
#include<stdio.h>

namespace hilbert{

void transform_ints_driver(double * int1, double *int2, double *U, int *nmopi, int *frzvpi, int nirrep, int nQ);

void transform_3index_tei(double *int2, double *U, int *nmopi, int* U_offset,
                          int* nmo_offset, int nirrep, int nQ, int max_nmopi,
                          int nmo_tot, int ngem_tot_lt, int max_num_threads);

void transform_oei(double *int1, double *U, int *nmopi, int * frzvpi, int* U_offset,
                   int* nmo_offset, int nirrep, int max_nmopi,
                   int nmo_tot, int max_num_threads);

}

#endif
