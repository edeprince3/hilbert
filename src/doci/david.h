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

#ifndef DAVID_H
#define DAVID_H

namespace hilbert{

// callback function
typedef void (*CallbackType)(size_t N, size_t maxdim,double **sigma, double **b, void * data);

int david_direct_redo(double *Adiag, int N, int M, double *eps, double **v, double cutoff, int print, CallbackType function, size_t & iter, void * data, int maxdim);

size_t david_direct(double *Adiag, size_t N, size_t M, double *eps, double **v, double cutoff,size_t print, CallbackType function, size_t & iter, void * data);
size_t david_in_core(double **A, size_t N, size_t M, double *eps, double **v,double cutoff, size_t print, size_t & iter);

}

#endif
