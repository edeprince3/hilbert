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
#include <math.h>
#include "legendre.h"

/******************************************************************************/

void Legendre::legendre_compute_glr ( int n, double x[], double w[] )

/******************************************************************************/
/*
 *   Purpose:
 *
 *   LEGENDRE_COMPUTE_GLR: Legendre quadrature by the Glaser-Liu-Rokhlin method.
 *
 *   Licensing:
 *
 *   This code is distributed under the GNU LGPL license. 
 *
 *   Modified:
 *
 *   19 October 2009
 *
 *   Author:
 *
 *   Original C version by Nick Hale.
 *   This C version by John Burkardt.
 *
 *   Reference:
 *
 *   Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
 *   A fast algorithm for the calculation of the roots of special functions, 
 *   SIAM Journal on Scientific Computing,
 *   Volume 29, Number 4, pages 1420-1438, 2007.
 *
 *   Parameters:
 *
 *   Input, int N, the order.
 *
 *   Output, double X[N], the abscissas.
 *
 *   Output, double W[N], the weights.
 *                                                             */
{
  double p, pp, w_sum;
/*
 *   Get the value and derivative of the N-th Legendre polynomial at 0.
 *   */
  legendre_compute_glr0 ( n, &p, &pp );
/*
 *   Either zero is a root, or we have to call a function to find the first root.
 *   */
  if ( n % 2 == 1 ){
    x[(n-1)/2] = p;
    w[(n-1)/2] = pp;
  }else{
    legendre_compute_glr2 ( p, n, &x[n/2], &w[n/2] );
  }
/*
 *   Get the complete set of roots and derivatives.
 *   */
  legendre_compute_glr1 ( n, x, w );
/*
 *   Compute the weights.
 *   */
  for (int i = 0; i < n; i++ ){
    w[i] = 2.0 / ( 1.0 - x[i] ) / ( 1.0 + x[i] ) / w[i] / w[i];
  }
  w_sum = 0.0;
  for (int i = 0; i < n; i++ ){
    w_sum = w_sum + w[i];
  }
  for (int i = 0; i < n; i++ ){
    w[i] = 2.0 * w[i] / w_sum;
  }
  return;
}

/******************************************************************************/

void Legendre::legendre_compute_glr0 ( int n, double *p, double *pp )

/******************************************************************************/
/*
 *   Purpose:
 *
 *   LEGENDRE_COMPUTE_GLR0 gets a starting value for the fast algorithm.
 *
 *   Licensing:
 *
 *   This code is distributed under the GNU LGPL license. 
 *
 *   Modified:
 *
 *   19 October 2009
 *
 *   Author:
 *
 *   Original C version by Nick Hale.
 *   This C version by John Burkardt.
 *
 *   Reference:
 *
 *   Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
 *   A fast algorithm for the calculation of the roots of special functions, 
 *   SIAM Journal on Scientific Computing,
 *   Volume 29, Number 4, pages 1420-1438, 2007.
 *
 *   Parameters:
 *   Input, int N, the order of the Legendre polynomial.
 *   Output, double *P, *PP, the value of the N-th Legendre polynomial
 *   and its derivative at 0.
 *                                                             */
{
  double dk;
  double pm1 = 1.0;
  double pm2 = 0.0, ppm1 = 0.0, ppm2 = 0.0;

  for (int k = 0; k < n; k++ ){
    dk = ( double ) k;
    *p = - dk * pm2 / ( dk + 1.0 );
    *pp = ( ( 2.0 * dk + 1.0 ) * pm1 - dk * ppm2 ) / ( dk + 1.0 );
    pm2 = pm1;
    pm1 = *p;
    ppm2 = ppm1;
    ppm1 = *pp;
  }
  return;
}

/******************************************************************************/

void Legendre::legendre_compute_glr1 ( int n, double *x, double *ders )

/******************************************************************************/
/*
 *   Purpose:
 *
 *   LEGENDRE_COMPUTE_GLR1 gets the complete set of Legendre points and weights.
 *
 *   Discussion:
 *
 *   This routine requires that a starting estimate be provided for one
 *   root and its derivative.  This information will be stored in entry
 *   (N+1)/2 if N is odd, or N/2 if N is even, of ROOTS and DERS.
 *
 *   Licensing:
 *
 *   This code is distributed under the GNU LGPL license. 
 *
 *   Modified:
 *
 *   19 October 2009
 *
 *   Author:
 *
 *   Original C version by Nick Hale.
 *   This C version by John Burkardt.
 *
 *   Reference:
 *
 *   Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
 *   A fast algorithm for the calculation of the roots of special functions, 
 *   SIAM Journal on Scientific Computing,
 *   Volume 29, Number 4, pages 1420-1438, 2007.
 *
 *   Parameters:
 *
 *   Input, int N, the order of the Legendre polynomial.
 *
 *   Input/output, double X[N].  On input, a starting value
 *   has been set in one entry.  On output, the roots of the Legendre 
 *   polynomial.
 *
 *   Input/output, double DERS[N].  On input, a starting value
 *   has been set in one entry.  On output, the derivatives of the Legendre 
 *   polynomial at the zeros.
 *
 *   Local Parameters:
 *
 *   Local, int M, the number of terms in the Taylor expansion.
 *                                                                                                 */
{
  double dk, dn, h;
  int n2, s;
  int m = 30;
  double *u, *up, xp;
  double pi = M_PI;

  if ( n % 2 == 1 ){
    n2 = ( n - 1 ) / 2;
    s = 1;
  }else{
    n2 = n / 2;
    s = 0;
  }

  u = ( double * ) malloc ( ( m + 2 ) * sizeof ( double ) );
  up = ( double * ) malloc ( ( m + 1 ) * sizeof ( double ) );
  dn = ( double ) n;

  for (int j = n2; j < n - 1; j++ ){
    xp = x[j];

    h = rk2_leg (pi/2.0, -pi/2.0, xp, n ) - xp;

    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = ders[j];

    up[0] = 0.0;
    up[1] = u[2];
    for (int k = 0; k <= m - 2; k++ ){
      dk = ( double ) k;

      u[k+3] =
      (
        2.0 * xp * ( dk + 1.0 ) * u[k+2]
        + ( dk * ( dk + 1.0 ) - dn * ( dn + 1.0 ) ) * u[k+1] / ( dk + 1.0 )
      ) / ( 1.0 - xp ) / ( 1.0 + xp ) / ( dk + 2.0 );

      up[k+2] = ( dk + 2.0 ) * u[k+3];
    }
    for (int l = 0; l < 5; l++ ){
      h = h - ts_mult ( u, h, m ) / ts_mult ( up, h, m-1 );
    }

    x[j+1] = xp + h;
    ders[j+1] = ts_mult ( up, h, m-1 );
  }

  free ( u );
  free ( up );

  for (int k = 0; k < n2 + s; k++ ){
    x[k] = - x[n-k-1];
    ders[k] = ders[n-k-1];
  }
  return;
}

/******************************************************************************/

void Legendre::legendre_compute_glr2 ( double pn0, int n, double *x1,  double *d1 )

/******************************************************************************/
/*
 *   Purpose:
 *
 *   LEGENDRE_COMPUTE_GLR2 finds the first real root.
 *
 *   Discussion:
 *
 *   This routine is only called if N is even.
 *
 *   Thanks to Morten Welinder, for pointing out a typographical error
 *   in indexing, 17 May 2013.
 *
 *   Licensing:
 *
 *   This code is distributed under the GNU LGPL license. 
 *
 *   Modified:
 *
 *   17 May 2013
 *
 *   Author:
 *
 *   Original C version by Nick Hale.
 *   This C version by John Burkardt.
 *
 *   Reference:
 *
 *   Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
 *   A fast algorithm for the calculation of the roots of special functions, 
 *   SIAM Journal on Scientific Computing,
 *   Volume 29, Number 4, pages 1420-1438, 2007.
 *
 *   Parameters:
 *
 *   Input, double PN0, the value of the N-th Legendre polynomial at 0.
 *
 *   Input, int N, the order of the Legendre polynomial.
 *
 *   Output, double *X1, the first real root.
 *
 *   Output, double *D1, the derivative at X1.
 *
 *   Local Parameters:
 *
 *   Local, int M, the number of terms in the Taylor expansion.
 *                                                                                     */
{
  double dk, dn, t, *u, *up;
  int m = 30;
  double pi = M_PI;
  t = 0.0;
  *x1 = rk2_leg ( t, -pi/2.0, 0.0, n );

  u = ( double * ) malloc ( ( m + 2 ) * sizeof ( double ) );
  up = ( double * ) malloc ( ( m + 1 ) * sizeof ( double ) );

  dn = ( double ) n;
/*
 *   U[0] and UP[0] are never used.
 *     U[M+1] is set, but not used, and UP[M] is set and not used.
 *       What gives?
 *       */
  u[0] = 0.0;
  u[1] = pn0;

  up[0] = 0.0;

  for (int k = 0; k <= m - 2; k = k + 2 ){
    dk = ( double ) k;

    u[k+2] = 0.0;
    u[k+3] = ( dk * ( dk + 1.0 ) - dn * ( dn + 1.0 ) ) * u[k+1]
      / ( dk + 1.0 ) / ( dk + 2.0 );

    up[k+1] = 0.0;
    up[k+2] = ( dk + 2.0 ) * u[k+3];
  }

  for (int l = 0; l < 5; l++ ){
    *x1 = *x1 - ts_mult ( u, *x1, m ) / ts_mult ( up, *x1, m-1 );
  }
  *d1 = ts_mult ( up, *x1, m-1 );

  free ( u );
  free ( up) ;

  return;
}

/******************************************************************************/

void Legendre::rescale ( double a, double b, int n, double x[], double w[] ){

/******************************************************************************/
/*
 *   Purpose:
 *
 *   RESCALE rescales a Legendre quadrature rule from [-1,+1] to [A,B].
 *
 *   Licensing:
 *
 *   This code is distributed under the GNU LGPL license. 
 *
 *   Modified:
 *
 *   22 October 2009
 *
 *   Author:
 *
 *   Original MATLAB version by Nick Hale.
 *   C version by John Burkardt.
 *
 *   Reference:
 *
 *   Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
 *   A fast algorithm for the calculation of the roots of special functions, 
 *   SIAM Journal on Scientific Computing,
 *   Volume 29, Number 4, pages 1420-1438, 2007.
 *
 *   Parameters:
 *
 *   Input, double A, B, the endpoints of the new interval.
 *
 *   Input, int N, the order.
 *
 *   Input/output, double X[N], on input, the abscissas for [-1,+1].
 *   On output, the abscissas for [A,B].
 *
 *   Input/output, double W[N], on input, the weights for [-1,+1].
 *   On output, the weights for [A,B].
 *                                                                         */
  for (int i = 0; i < n; i++ ){
    x[i] = ( ( a + b ) + ( b - a ) * x[i] ) / 2.0;
  }
  for (int i = 0; i < n; i++ ){
    w[i] = ( b - a ) * w[i] / 2.0;
  }
  return;
}
/******************************************************************************/

double Legendre::rk2_leg ( double t1, double t2, double x, int n )

/******************************************************************************/
/*
 *   Purpose:
 *
 *   RK2_LEG advances the value of X(T) using a Runge-Kutta method.
 *
 *   Licensing:
 *
 *   This code is distributed under the GNU LGPL license. 
 *
 *   Modified:
 *
 *   22 October 2009
 *
 *   Author:
 *
 *   Original C version by Nick Hale.
 *   This C version by John Burkardt.
 *
 *   Parameters:
 *
 *   Input, double T1, T2, the range of the integration interval.
 *
 *   Input, double X, the value of X at T1.
 *
 *   Input, int N, the number of steps to take.
 *
 *   Output, double RK2_LEG, the value of X at T2.
 *                                               */
{
  double f, h, k1, k2, snn1, t;
  int m = 10;

  h = ( t2 - t1 ) / ( double ) m;
  snn1 = sqrt ( ( double ) ( n * ( n + 1 ) ) );

  t = t1;

  for (int j = 0; j < m; j++ ){
    f = ( 1.0 - x ) * ( 1.0 + x );
    k1 = - h * f / ( snn1 * sqrt ( f ) - 0.5 * x * sin ( 2.0 * t ) );
    x = x + k1;
    t = t + h;
    f = ( 1.0 - x ) * ( 1.0 + x );
    k2 = - h * f / ( snn1 * sqrt ( f ) - 0.5 * x * sin ( 2.0 * t ) );
    x = x + 0.5 * ( k2 - k1 );
  }
  return x;
}

/******************************************************************************/

double Legendre::ts_mult ( double *u, double h, int n )

/******************************************************************************/
/*
 *   Purpose:
 *
 *   TS_MULT evaluates a polynomial.
 *
 *   Discussion:
 *
 *   TS_MULT = U[1] + U[2] * H + ... + U[N] * H^(N-1).
 *
 *   Licensing:
 *
 *   This code is distributed under the GNU LGPL license. 
 *
 *   Modified:
 *
 *   17 May 2013
 *
 *   Author:
 *
 *   Original C version by Nick Hale.
 *   This C version by John Burkardt.
 *
 *   Parameters:
 *
 *   Input, double U[N+1], the polynomial coefficients.
 *   U[0] is ignored.
 *
 *   Input, double H, the polynomial argument.
 *
 *   Input, int N, the number of terms to compute.
 *
 *   Output, double TS_MULT, the value of the polynomial.
 *                                                         */
{
  double hk, ts;

  ts = 0.0;
  hk = 1.0;
  for (int k = 1; k<= n; k++ ){
    ts = ts + u[k] * hk;
    hk = hk * h;
  }
  return ts;
}

