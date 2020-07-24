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

#ifndef LEGENDRE_H
#define LEGENDRE_H

class Legendre {
    public:
       void legendre_compute_glr ( int n, double x[], double w[] );
       void legendre_compute_glr0 ( int n, double *p, double *pp );
       void legendre_compute_glr1 ( int n, double *roots, double *ders );
       void legendre_compute_glr2 ( double p, int n, double *roots, double *ders );
       void rescale ( double a, double b, int n, double x[], double w[] );
       double rk2_leg ( double t, double tn, double x, int n );
       double ts_mult ( double *u, double h, int n );
};

#endif
