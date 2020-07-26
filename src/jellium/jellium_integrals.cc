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

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

#include <psi4/psi4-dec.h>
#include <psi4/libpsi4util/PsiOutStream.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/matrix.h>

#include "jellium_integrals.h"
#include "legendre.h"

#include <misc/omp.h>

using namespace psi;

namespace hilbert{

JelliumIntegrals::JelliumIntegrals(Options & options):
        options_(options)
{

    outfile->Printf("\n");
    outfile->Printf("    ==> Jellium Integral Construction <==\n");
    outfile->Printf("\n");

    orbitalMax = options.get_int("N_BASIS_FUNCTIONS");
    electrons = options.get_int("N_ELECTRONS");
    symmetry = true;//options.get_bool("SYMMETRY");
    compute_integrals();

    //x1 = (int *)malloc(3*sizeof(int));
    //x2 = (int *)malloc(3*sizeof(int));
    //y1 = (int *)malloc(3*sizeof(int));
    //y2 = (int *)malloc(3*sizeof(int));
    //z1 = (int *)malloc(3*sizeof(int));
    //z2 = (int *)malloc(3*sizeof(int));
}

// free memory here
JelliumIntegrals::~JelliumIntegrals()
{
  //free(x1);
  //free(x2);
  //free(y1);
  //free(y2);
  //free(z1);
  //free(z2);
}

void JelliumIntegrals::compute_integrals() {

        //printf ( "\n" );
        //printf ( "LEGENDRE_RULE_FAST:\n" );
        //printf ( "  Normal end of execution.\n" );

        //printf ( "\n" );
        double a = 0.0;
        double b = 1.0;
        int n = options_.get_int("N_GRID_POINTS");
        double *x;
        int *mu, *nu, *sig, *lam;

        x   = (double *)malloc(n*sizeof(double));
        w   = (double *)malloc(n*sizeof(double));
        grid_points = x; 

        sig  = (int*)malloc(3*sizeof(int));
        lam  = (int*)malloc(3*sizeof(int));

        nmax=30;
        //TODO: make this vector irrep from the get go
        std::shared_ptr<Vector> ORBE = std::shared_ptr<Vector>( new Vector(3*nmax*nmax*nmax));//VEC_INT(3*nmax*nmax*nmax);
        MO  = MAT_INT(3*nmax*nmax*nmax,3);
        OrderPsis3D(nmax, ORBE->pointer(), MO);

        //set up symmetry
        if(symmetry){
            Orderirrep(nmax, ORBE->pointer(), MO, electrons);
        }else{
            nirrep_ = 1;
            nsopi_ = (int*)malloc(sizeof(int));
            doccpi_ = (int*)malloc(sizeof(int));
            //TODO if not restricted this wont work
            nsopi_[0] = orbitalMax;
            doccpi_[0] = electrons/2;
        }

        Legendre tmp;
        //  Constructe grid and weights, store them to the vectors x and w, respectively.
        //  This is one of John Burkhardt's library functions
        tmp.legendre_compute_glr(n, x, w);

        // Scale the grid to start at value a and end at value b. 
        // We want our integration range to go from 0 to 1, so a = 0, b = 1
        // This is also one of John Burkhardt's library functions
        tmp.rescale( a, b, n, x, w);

        // build g tensor g[npq] * w[n]
        outfile->Printf("\n");
        outfile->Printf("    build g tensor................"); fflush(stdout);
        g_tensor = std::shared_ptr<Vector>( new Vector(n * 2 * (nmax + 1) * 2 * (nmax + 1)));
        for (int pt = 0; pt < n; pt++) {
                double xval = x[pt];
                for (int p = 0; p <= nmax*2; p++) {
                        for (int q = 0; q <= nmax*2; q++) {
                                g_tensor->pointer()[pt*2*nmax*2*nmax+p*2*nmax+q] = g_pq(p, q, xval) * w[pt];
                        }
                }
        }
        outfile->Printf("done.\n");

        // build sqrt(x*x+y*y+z*z)
        outfile->Printf("    build sqrt tensor............."); fflush(stdout);
        sqrt_tensor = std::shared_ptr<Vector>(new Vector(n*n*n));
        double * s_p = sqrt_tensor->pointer();
        for (int i = 0; i < n; i++) {
                double xval = x[i];
                for (int j = 0; j < n; j++) {
                        double yval = x[j];
                        for (int k = 0; k < n; k++) {
                                double zval = x[k];
                                double val = sqrt(xval*xval+yval*yval+zval*zval);
                                s_p[i*n*n + j*n + k] = 1.0/val;
                        }
                }
        }
        outfile->Printf("done.\n");


        unsigned long start_pq = clock();
        // now, compute (P|Q)
        outfile->Printf("    build (P|Q)..................."); fflush(stdout);
        PQmap = (int ***)malloc((2*nmax+1)*sizeof(int**));
        for (int i = 0; i < 2*nmax+1; i++) {
                PQmap[i] = (int **)malloc((2*nmax+1)*sizeof(int*));
                for (int j = 0; j < 2*nmax+1; j++) {
                        PQmap[i][j] = (int *)malloc((2*nmax+1)*sizeof(int));
                        for (int k = 0; k < 2*nmax+1; k++) {
                                PQmap[i][j][k] = 999;
                        }
                }
        }
        int Pdim = 0;
        for (int px = 0; px < 2*nmax+1; px++) {
                for (int py = 0; py < 2*nmax+1; py++) {
                        for (int pz = 0; pz < 2*nmax+1; pz++) {
                                PQmap[px][py][pz] = Pdim;
                                Pdim++;
                        }
                }
        }
        //printf("1, -2, -3, 4, 2, 0\n"); 
        //pq_int(orbitalMax, x, w, 1, -2, -3, 4, 2, 0);

        PQ = std::shared_ptr<Matrix>(new Matrix(Pdim,Pdim));
        double ** PQ_p = PQ->pointer();
        Ke = std::shared_ptr<Matrix>(new Matrix(nirrep_,nsopi_,nsopi_));
        NucAttrac = std::shared_ptr<Matrix>(new Matrix(nirrep_,nsopi_,nsopi_));

        #pragma omp parallel
        {
        #pragma omp for schedule(dynamic) nowait 
                for (int px = 0; px < 2*nmax+1; px++) {
                        for (int qx = px; qx < 2*nmax+1; qx++) {

                                int pq_x = px*(2*nmax+1) + qx;

                                for (int py = 0; py < 2*nmax+1; py++) {
                                        for (int qy = py; qy < 2*nmax+1; qy++) {

                                                int pq_y = py*(2*nmax+1) + qy;
                                                if ( pq_x > pq_y ) continue;

                                                for (int pz = 0; pz < 2*nmax+1; pz++) {
                                                        for (int qz = pz; qz < 2*nmax+1; qz++) {

                                                                int pq_z = pz*(2*nmax+1) + qz;
                                                                if ( pq_y > pq_z ) continue;

                                                                //if ( P > Q ) continue;
                                                                if((px+qx)%2==0 && (py+qy)%2==0 && (pz+qz)%2==0){
                                                                        double dum = pq_int_new(n, px, py, pz, qx, qy, qz);
                                                                        int P,Q;

                                                                        // start 
                                                                        P = PQmap[px][py][pz];
                                                                        Q = PQmap[qx][qy][qz];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qx][py][pz];
                                                                        Q = PQmap[px][qy][qz];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[px][qy][pz];
                                                                        Q = PQmap[qx][py][qz];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qx][qy][pz];
                                                                        Q = PQmap[px][py][qz];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[px][py][qz];
                                                                        Q = PQmap[qx][qy][pz];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qx][py][qz];
                                                                        Q = PQmap[px][qy][pz];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[px][qy][qz];
                                                                        Q = PQmap[qx][py][pz];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qx][qy][qz];
                                                                        Q = PQmap[px][py][pz];
                                                                        PQ_p[P][Q] = dum;

                                                                        // pxqx - pyqy

                                                                        P = PQmap[py][px][pz];
                                                                        Q = PQmap[qy][qx][qz];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[py][qx][pz];
                                                                        Q = PQmap[qy][px][qz];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qy][px][pz];
                                                                        Q = PQmap[py][qx][qz];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qy][qx][pz];
                                                                        Q = PQmap[py][px][qz];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[py][px][qz];
                                                                        Q = PQmap[qy][qx][pz];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[py][qx][qz];
                                                                        Q = PQmap[qy][px][pz];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qy][px][qz];
                                                                        Q = PQmap[py][qx][pz];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qy][qx][qz];
                                                                        Q = PQmap[py][px][pz];
                                                                        PQ_p[P][Q] = dum;

                                                                        // now begins pxqx < pyqy < pzqz

                                                                        P = PQmap[pz][px][py];
                                                                        Q = PQmap[qz][qx][qy];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[pz][qx][py];
                                                                        Q = PQmap[qz][px][qy];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[pz][px][qy];
                                                                        Q = PQmap[qz][qx][py];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[pz][qx][qy];
                                                                        Q = PQmap[qz][px][py];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qz][px][py];
                                                                        Q = PQmap[pz][qx][qy];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qz][qx][py];
                                                                        Q = PQmap[pz][px][qy];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qz][px][qy];
                                                                        Q = PQmap[pz][qx][py];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qz][qx][qy];
                                                                        Q = PQmap[pz][px][py];
                                                                        PQ_p[P][Q] = dum;

                                                                        // pxqx - pyqy

                                                                        P = PQmap[pz][py][px];
                                                                        Q = PQmap[qz][qy][qx];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[pz][py][qx];
                                                                        Q = PQmap[qz][qy][px];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[pz][qy][px];
                                                                        Q = PQmap[qz][py][qx];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[pz][qy][qx];
                                                                        Q = PQmap[qz][py][px];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qz][py][px];
                                                                        Q = PQmap[pz][qy][qx];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qz][py][qx];
                                                                        Q = PQmap[pz][qy][px];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qz][qy][px];
                                                                        Q = PQmap[pz][py][qx];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qz][qy][qx];
                                                                        Q = PQmap[pz][py][px];
                                                                        PQ_p[P][Q] = dum;

                                                                        // now begins last set of 16

                                                                        P = PQmap[px][pz][py];
                                                                        Q = PQmap[qx][qz][qy];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qx][pz][py];
                                                                        Q = PQmap[px][qz][qy];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[px][pz][qy];
                                                                        Q = PQmap[qx][qz][py];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qx][pz][qy];
                                                                        Q = PQmap[px][qz][py];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[px][qz][py];
                                                                        Q = PQmap[qx][pz][qy];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qx][qz][py];
                                                                        Q = PQmap[px][pz][qy];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[px][qz][qy];
                                                                        Q = PQmap[qx][pz][py];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qx][qz][qy];
                                                                        Q = PQmap[px][pz][py];
                                                                        PQ_p[P][Q] = dum;

                                                                        // pxqx - pyqy

                                                                        P = PQmap[py][pz][px];
                                                                        Q = PQmap[qy][qz][qx];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[py][pz][qx];
                                                                        Q = PQmap[qy][qz][px];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qy][pz][px];
                                                                        Q = PQmap[py][qz][qx];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qy][pz][qx];
                                                                        Q = PQmap[py][qz][px];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[py][qz][px];
                                                                        Q = PQmap[qy][pz][qx];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[py][qz][qx];
                                                                        Q = PQmap[qy][pz][px];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qy][qz][px];
                                                                        Q = PQmap[py][pz][qx];
                                                                        PQ_p[P][Q] = dum;

                                                                        P = PQmap[qy][qz][qx];
                                                                        Q = PQmap[py][pz][px];
                                                                        PQ_p[P][Q] = dum;

                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
        unsigned long end_pq = clock();
        outfile->Printf("done.\n");fflush(stdout);
        outfile->Printf("\n");
        outfile->Printf("    time for (P|Q) construction:                %6.1f s\n",(double)(end_pq-start_pq)/CLOCKS_PER_SEC); fflush(stdout);
        //printf("    time for (P|Q) construction:                %6.1f s\n",(double)(end_pq-start_pq)/CLOCKS_PER_SEC);
        //printf("value of 000|000 %20.20f\n",PQ_p[PQmap[0][0][0]][PQmap[0][0][0]]);
        outfile->Printf("\n");
        
        //testing a small pq setup
        int offset = 0;
        for(int i = 0; i <=2*nmax; i++){
           for(int j = i; j <=2*nmax; j++){
               if((i+j)%2==0){
                  offset++;
               }
           }
        }
        int count = 0;
        int count2 = 0;
        int P = 0;
        int Q = 0;
        //printf("%d \n",offset*offset*offset);

        //Creating small PQ
        //TODO disabled for now due to symmetry assumtions made
        PQ_small = (double*)malloc(offset*offset*offset*sizeof(double));
        for(int px = 0; px <= 2*nmax; px++){
           for(int qx = px; qx <= 2*nmax; qx++){
              for(int py = 0; py <= 2*nmax; py++){
                 for(int qy = py; qy <= 2*nmax; qy++){
                    for(int pz = 0; pz <= 2*nmax; pz++){
                       for(int qz = pz; qz <= 2*nmax; qz++){
                          //count2++;
                          //if(px == px_t && py == py_t && pz == pz_t && qx == qx_t && qy == qy_t && qz == qz_t){

                          //    printf("%d \n",count);}
                          if((px+qx)%2==0 && (py+qy)%2==0 && (pz+qz)%2==0){
                              P = PQmap[px][py][pz];
                              Q = PQmap[qx][qy][qz];
                              PQ_small[count] = PQ_p[P][Q];
                              count++;
                          }
                       }
                    }
                 }
              }
           }
        }
        offset_pq = 0;
        for(int i = 0; i <=2*nmax; i++){
           for(int j = i; j <=2*nmax; j++){
               if((i+j)%2==0){
                  offset_pq++;
               }
           }
        }
        //outfile->Printf("canonical integrals");

        // Four nested loops to compute lower triange of electron repulsion integrals - roughly have of the non-unique integrals
        // will not be computed, but this is still not exploiting symmetry fully
        outfile->Printf("    build potential integrals.....");fflush(stdout);
        unsigned long start = clock();
        #pragma omp parallel
        {
        #pragma omp for schedule(dynamic) nowait
                for(int h = 0; h < nirrep_;h++){
                        int offset = 0;
                        double** Ke_p = Ke->pointer(h);
                        for(int i = 0; i < h; i ++){
                                offset += nsopi_[i];
                        }
                        double** Nu_p = NucAttrac->pointer(h);
                        for (int i=0; i< nsopi_[h]; i++) {
                                int* mu;
                                int* nu;
                                mu = MO[i+offset];

                                for (int j=i; j< nsopi_[h]; j++) {
                                        nu = MO[j+offset];

                                        // Kinetic Energy Integrals - already computed and stored in ORBE vector    
                                        if (i==j) { 
                                                Ke_p[i][j] = 0.5*ORBE->pointer()[i+offset];
                                        }
                                        // Nuclear-attraction Integrals
                                        double dum = Vab_Int_new(n, x, w, mu, nu);
                                        Nu_p[i][j] = dum;
                                        Nu_p[j][i] = dum;

                                }
                        }
                        //offset += nsopi_[h];
                }
        }
        unsigned long end = clock();
        outfile->Printf("done.\n");fflush(stdout);
        outfile->Printf("\n");
        outfile->Printf("    time for potential integral construction:   %6.1f s\n",(double)(end-start)/CLOCKS_PER_SEC); fflush(stdout);
        outfile->Printf("\n");
        // Compute self energy
        selfval = E0_Int(n, x, w);
        //PQ = NULL;
        sqrt_tensor = NULL;
        for(int j = 0; j < 2*nmax; j++){
        int pqx = 0;
        int tmp = 1;
        for(int i = 0; i < j; i++){
            pqx += nmax+tmp;
            if(i%2==0){
               tmp--;
            }
        }
        if(j%2==0){
           //printf("test: %d\n",(nmax*j+j-(j*j)/4));
           //printf("j: %d pqx: %d\n",j,pqx);
        }
        if(j%2==1){
           //printf("test: %d\n",(nmax*j)-(j*j)/4+3/4*j+1/2);
           //printf("j: %d pqx: %d\n",j,pqx);
        }
        }

}

//Commented out code broken by removing symmetry
double JelliumIntegrals::ERI(int a, int b, int c, int d){

    if((MO[a][0]+MO[b][0]+MO[c][0]+MO[d][0])%2==1){
      return 0.0;
    }
    
    //return ERI_unrolled_test(MO[a], MO[b], MO[c], MO[d], PQ->pointer(), PQmap);
    return ERI_unrolled_new(MO[a], MO[b], MO[c], MO[d], PQ->pointer(), PQmap);
}

double JelliumIntegrals::g_pq(int p, int q, double x) {
    int d = abs(p-q);
    double pi = M_PI;
    ////if(q < 0 || p < 0){
    ////   return 0;
    ////}
    //if (p == q && p == 0) {
    //  return 1.0 - x;
    //}
    //else if ( p == q && p > 0 ) {
    //  return (1.0 - x)*cos(p*pi*x)/2.0 - sin(p*pi*x)/(2*p*pi);
    //}
    //else if ( (d % 2)==0 && d!=0) {
    //  return (q*sin(q*pi*x) - p*sin(p*pi*x))/((p*p-q*q)*pi);
    //}
    //else 
    //  return 0.0;

    if ( p == q ) {
        if ( p == 0 ) {
            return 1.0 - x;
        }else {
            return (1.0 - x)*cos(p*pi*x)/2.0 - sin(p*pi*x)/(2*p*pi);
        }
    }else if ( (d % 2)==0 && d != 0 ) {
        return (q*sin(q*pi*x) - p*sin(p*pi*x))/((p*p-q*q)*pi);
    }   
    return 0.0;
}


// From Eq. 3.6 in the Peter Gill paper, 
// the Vab integrals are -1/pi^3 \int (phi_a phi_b )/(|r1 - r2|) dr1 dr2
// This essentially leads to (p|q) integrals with 
// px = a_x - b_x
// py = a_y - b_y
// pz = a_z - b_z
// qx = a_x + b_x
// qy = a_y + b_y
// qz = a_z + b_z
// Note the expansion of the trigonetric identity:
// Cos[px x1] Cos[py y1] Cos[pz z1] - Cos[qx x1] Cos[py y1] Cos[pz z1] - 
// Cos[px x1] Cos[qy y1] Cos[pz z1] + Cos[qx x1] Cos[qy y1] Cos[pz z1] -
// Cos[px x1] Cos[py y1] Cos[qz z1] + 
// Cos[qx x1] Cos[py y1] Cos[qz z1] + Cos[px x1] Cos[qy y1] Cos[qz z1] -
// Cos[qx x1] Cos[qy y1] Cos[qz z1]
// In order to be consistent with the defintiion of the (p|q) integrals, 
// the term Cos[px x1] Cos[py y1] Cos[pz z1] -> Cos[px x1] Cos[py y1] Cos[pz z1] Cos[0 x2] Cos[0 y2] Cos[0 z2]
// In terms of how the pq_int function is called for the above integral, it should be
// pq_int(dim, xa, w, px, py, pz, 0, 0, 0)

double JelliumIntegrals::Vab_Int_new(int dim, double *xa, double *w, int *a, int *b) {

    int px, py, pz, qx, qy, qz;
    double Vab;
    px = abs(a[0] - b[0]);
    py = abs(a[1] - b[1]);
    pz = abs(a[2] - b[2]);
    

    qx = a[0] + b[0];
    qy = a[1] + b[1];
    qz = a[2] + b[2];
    if(px%2==1 || py %2==1 || pz%2==1 || qx%2==1 || qy%2==1 || qz%2==1){
       return 0.0;
    }
    Vab = 0.0;
    int P, Q;
    double ** PQ_p = PQ->pointer(); 

    //    Cos[px x1] Cos[py y1] Cos[pz z1]
    //Vab  += pq_int_new(dim, px, py, pz,  0,  0,  0);
    P = PQmap[px][py][pz];
    Q = PQmap[0][0][0];
    Vab += PQ_p[P][Q];
         
    // -  Cos[qx x1] Cos[py y1] Cos[pz z1]
    //Vab  -= pq_int_new(dim,  0, py, pz, qx,  0,  0);
    P = PQmap[0][py][pz];
    Q = PQmap[qx][0][0];
    Vab -= PQ_p[P][Q];
         
    // -  Cos[px x1] Cos[qy y1] Cos[pz z1]
    //Vab  -= pq_int_new(dim, px,  0, pz,  0, qy,  0);
    P = PQmap[px][0][pz];
    Q = PQmap[0][qy][0];
    Vab -= PQ_p[P][Q];
         
    // +  Cos[qx x1] Cos[qy y1] Cos[pz z1]
    //Vab  += pq_int_new(dim,  0,  0, pz, qx, qy,  0);
    P = PQmap[0][0][pz];
    Q = PQmap[qx][qy][0];
    Vab += PQ_p[P][Q];
         
    // - Cos[px x1] Cos[py y1] Cos[qz z1]  
    //Vab  -= pq_int_new(dim, px, py,  0,  0,  0, qz);
    P = PQmap[px][py][0];
    Q = PQmap[0][0][qz];
    Vab -= PQ_p[P][Q];
         
    // + Cos[qx x1] Cos[py y1] Cos[qz z1] 
    //Vab  += pq_int_new(dim,  0, py,  0, qx,  0, qz);
    P = PQmap[0][py][0];
    Q = PQmap[qx][0][qz];
    Vab += PQ_p[P][Q];
         
    // +  Cos[px x1] Cos[qy y1] Cos[qz z1]
    //Vab  += pq_int_new(dim, px,  0,  0,  0, qy, qz);
    P = PQmap[px][0][0];
    Q = PQmap[0][qy][qz];
    Vab += PQ_p[P][Q];
         
    // - Cos[qx x1] Cos[qy y1] Cos[qz z1]
    //Vab  -= pq_int_new(dim,  0,  0,  0, qx, qy, qz);
    P = PQmap[0][0][0];
    Q = PQmap[qx][qy][qz];
    Vab -= PQ_p[P][Q];

    return -Vab;

}

// the E integral is 1/pi^6 \int 1/(|r1 - r2|) dr1 dr2
// which is equivalent to 
// 1/pi^6 \int cos(0 x1) cos(0 y1) cos(0 z1) cos(0 x2) cos(0 y2) cos(0 z2)/|r1-r2| dr1 dr2
//
double JelliumIntegrals::E0_Int(int dim, double *xa, double *w) {
  return  pq_int(dim, xa, w, 0, 0, 0, 0, 0, 0);
}

// new on 7/24/20
double JelliumIntegrals::ERI_unrolled_new(int * a, int * b, int * c, int * d, double ** PQ, int *** PQmap) {

  //x1[0] = ax-bx, x1[1] = ax+bx
  int* x1 = (int *)malloc(3*sizeof(int));
  int* x2 = (int *)malloc(3*sizeof(int));
  int* y1 = (int *)malloc(3*sizeof(int));
  int* y2 = (int *)malloc(3*sizeof(int));
  int* z1 = (int *)malloc(3*sizeof(int));
  int* z2 = (int *)malloc(3*sizeof(int));
 
  x1[0] = abs(a[0] - b[0]);
  y1[0] = abs(a[1] - b[1]);
  z1[0] = abs(a[2] - b[2]);

  x1[1] = a[0] + b[0];
  y1[1] = a[1] + b[1];
  z1[1] = a[2] + b[2];

  //x1[0] abs(= cx-dx, x1)[1] = cx+dx
  x2[0] = abs(c[0] - d[0]);
  y2[0] = abs(c[1] - d[1]);
  z2[0] = abs(c[2] - d[2]);

  x2[1] = c[0] + d[0];
  y2[1] = c[1] + d[1];
  z2[1] = c[2] + d[2];
  // Generate all combinations of phi_a phi_b phi_c phi_d in expanded cosine form

  double eri_val = 0.0;

  int Q;
  int P;

  // 0,0,0
  Q = PQmap[ x2[0] ][ y2[0] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  eri_val -= PQ[P][Q];


  // 1,0,0

  Q = PQmap[ x2[1] ][ y2[0] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  eri_val += PQ[P][Q];

  // 0,1,0

  Q = PQmap[ x2[0] ][ y2[1] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  eri_val += PQ[P][Q];


  // 1,1,0
  Q = PQmap[ x2[1] ][ y2[1] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  // 0,0,1

  Q = PQmap[ x2[0] ][ y2[0] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  eri_val += PQ[P][Q];

  // 1,0,1

  Q = PQmap[ x2[1] ][ y2[0] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  // 0,1,1

  Q = PQmap[ x2[0] ][ y2[1] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  // 1,1,1

  Q = PQmap[ x2[1] ][ y2[1] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  eri_val += PQ[P][Q];

  free(x1);
  free(x2);
  free(y1);
  free(y2);
  free(z1);
  free(z2);

  return eri_val;

} 

double JelliumIntegrals::ERI_unrolled_test(int * a, int * b, int * c, int * d, double ** PQ, int *** PQmap) {

  //x1[0] = ax-bx, x1[1] = ax+bx
  int* x1 = (int *)malloc(3*sizeof(int));
  int* x2 = (int *)malloc(3*sizeof(int));
  int* y1 = (int *)malloc(3*sizeof(int));
  int* y2 = (int *)malloc(3*sizeof(int));
  int* z1 = (int *)malloc(3*sizeof(int));
  int* z2 = (int *)malloc(3*sizeof(int));
 
  x1[0] = abs(a[0] - b[0]);
  y1[0] = abs(a[1] - b[1]);
  z1[0] = abs(a[2] - b[2]);

  x1[1] = a[0] + b[0];
  y1[1] = a[1] + b[1];
  z1[1] = a[2] + b[2];

  //x1[0] abs(= cx-dx, x1)[1] = cx+dx
  x2[0] = abs(c[0] - d[0]);
  y2[0] = abs(c[1] - d[1]);
  z2[0] = abs(c[2] - d[2]);

  x2[1] = c[0] + d[0];
  y2[1] = c[1] + d[1];
  z2[1] = c[2] + d[2];
  // Generate all combinations of phi_a phi_b phi_c phi_d in expanded cosine form

  double eri_val = 0.0;

  int Q = PQmap[ x2[0] ][ y2[0] ][ z2[0] ];
  int P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];

  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  Q = PQmap[ x2[1] ][ y2[0] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  eri_val += PQ[P][Q];

  Q = PQmap[ x2[0] ][ y2[0] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  eri_val += PQ[P][Q];

  Q = PQmap[ x2[1] ][ y2[0] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  Q = PQmap[ x2[0] ][ y2[1] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  eri_val += PQ[P][Q];

  Q = PQmap[ x2[1] ][ y2[1] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  Q = PQmap[ x2[0] ][ y2[1] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  Q = PQmap[ x2[1] ][ y2[1] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  eri_val += PQ[P][Q];

  Q = PQmap[ x2[0] ][ y2[0] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  eri_val += PQ[P][Q];

  Q = PQmap[ x2[1] ][ y2[0] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  Q = PQmap[ x2[0] ][ y2[0] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  Q = PQmap[ x2[1] ][ y2[0] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  eri_val += PQ[P][Q];

  Q = PQmap[ x2[0] ][ y2[1] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  Q = PQmap[ x2[1] ][ y2[1] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  eri_val += PQ[P][Q];

  Q = PQmap[ x2[0] ][ y2[1] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  eri_val += PQ[P][Q];

  Q = PQmap[ x2[1] ][ y2[1] ][ z2[0] ];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  Q = PQmap[ x2[0] ][ y2[0] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  eri_val += PQ[P][Q];

  Q = PQmap[ x2[1] ][ y2[0] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  Q = PQmap[ x2[0] ][ y2[0] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  Q = PQmap[ x2[1] ][ y2[0] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  eri_val += PQ[P][Q];

  Q = PQmap[ x2[0] ][ y2[1] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  Q = PQmap[ x2[1] ][ y2[1] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  eri_val += PQ[P][Q];

  Q = PQmap[ x2[0] ][ y2[1] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  eri_val += PQ[P][Q];

  Q = PQmap[ x2[1] ][ y2[1] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  eri_val -= PQ[P][Q];

  Q = PQmap[ x2[0] ][ y2[0] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  Q = PQmap[ x2[1] ][ y2[0] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  eri_val += PQ[P][Q];

  Q = PQmap[ x2[0] ][ y2[0] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  eri_val += PQ[P][Q];

  Q = PQmap[ x2[1] ][ y2[0] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  Q = PQmap[ x2[0] ][ y2[1] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  eri_val += PQ[P][Q];

  Q = PQmap[ x2[1] ][ y2[1] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  Q = PQmap[ x2[0] ][ y2[1] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  eri_val += PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  Q = PQmap[ x2[1] ][ y2[1] ][ z2[1] ];

  P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  eri_val -= PQ[P][Q];

  P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  eri_val += PQ[P][Q];

  free(x1);
  free(x2);
  free(y1);
  free(y2);
  free(z1);
  free(z2);

  return eri_val;

} 

//not used for now as this uses small_pq which does not work yet for the symmetry disabled version
double JelliumIntegrals::ERI_unrolled(int * a, int * b, int * c, int * d) {
  //x1[0] = ax-bx, x1[1] = ax+bx
  int* x1 = (int *)malloc(2*sizeof(int));
  int* x2 = (int *)malloc(2*sizeof(int));
  int* y1 = (int *)malloc(2*sizeof(int));
  int* y2 = (int *)malloc(2*sizeof(int));
  
int* z1 = (int *)malloc(2*sizeof(int));
  int* z2 = (int *)malloc(2*sizeof(int));
 
  
  x1[0] = abs(a[0] - b[0]);
  y1[0] = abs(a[1] - b[1]);
  z1[0] = abs(a[2] - b[2]);

  x1[1] = a[0] + b[0];
  y1[1] = a[1] + b[1];
  z1[1] = a[2] + b[2];

  //x1[0] abs(= cx-dx, x1)[1] = cx+dx
  x2[0] = abs(c[0] - d[0]);
  y2[0] = abs(c[1] - d[1]);
  z2[0] = abs(c[2] - d[2]);

  x2[1] = c[0] + d[0];
  y2[1] = c[1] + d[1];
  z2[1] = c[2] + d[2];
  // Generate all combinations of phi_a phi_b phi_c phi_d in expanded cosine form

  //double //eri_val = 0.0;
  double eri_val_test = 0.0;

  //int //Q = PQmap[ x2[0] ][ y2[0] ][ z2[0] ];

  //int //P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[0],y2[0],z2[0],x1[0],y1[0],z1[0]);

//P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[0],y2[0],z2[0],x1[1],y1[0],z1[0]);

//Q = PQmap[ x2[1] ][ y2[0] ][ z2[0] ];

//P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[1],y2[0],z2[0],x1[0],y1[0],z1[0]);

//P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[1],y2[0],z2[0],x1[1],y1[0],z1[0]);

//Q = PQmap[ x2[0] ][ y2[0] ][ z2[0] ];

//P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[0],y2[0],z2[0],x1[0],y1[1],z1[0]);

//P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[0],y2[0],z2[0],x1[1],y1[1],z1[0]);

//Q = PQmap[ x2[1] ][ y2[0] ][ z2[0] ];

//P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[1],y2[0],z2[0],x1[0],y1[1],z1[0]);

//P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[1],y2[0],z2[0],x1[1],y1[1],z1[0]);

//Q = PQmap[ x2[0] ][ y2[1] ][ z2[0] ];

//P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[0],y2[1],z2[0],x1[0],y1[0],z1[0]);

//P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[0],y2[1],z2[0],x1[1],y1[0],z1[0]);

//Q = PQmap[ x2[1] ][ y2[1] ][ z2[0] ];

//P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[1],y2[1],z2[0],x1[0],y1[0],z1[0]);

//P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[1],y2[1],z2[0],x1[1],y1[0],z1[0]);

//Q = PQmap[ x2[0] ][ y2[1] ][ z2[0] ];

//P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[0],y2[1],z2[0],x1[0],y1[1],z1[0]);

//P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[0],y2[1],z2[0],x1[1],y1[1],z1[0]);

//Q = PQmap[ x2[1] ][ y2[1] ][ z2[0] ];

//P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[1],y2[1],z2[0],x1[0],y1[1],z1[0]);

//P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[1],y2[1],z2[0],x1[1],y1[1],z1[0]);

//Q = PQmap[ x2[0] ][ y2[0] ][ z2[0] ];

//P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[0],y2[0],z2[0],x1[0],y1[0],z1[1]);

//P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[0],y2[0],z2[0],x1[1],y1[0],z1[1]);

//Q = PQmap[ x2[1] ][ y2[0] ][ z2[0] ];

//P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[1],y2[0],z2[0],x1[0],y1[0],z1[1]);

//P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[1],y2[0],z2[0],x1[1],y1[0],z1[1]);

//Q = PQmap[ x2[0] ][ y2[0] ][ z2[0] ];

//P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[0],y2[0],z2[0],x1[0],y1[1],z1[1]);

//P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[0],y2[0],z2[0],x1[1],y1[1],z1[1]);

//Q = PQmap[ x2[1] ][ y2[0] ][ z2[0] ];

//P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[1],y2[0],z2[0],x1[0],y1[1],z1[1]);

//P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[1],y2[0],z2[0],x1[1],y1[1],z1[1]);

//Q = PQmap[ x2[0] ][ y2[1] ][ z2[0] ];

//P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[0],y2[1],z2[0],x1[0],y1[0],z1[1]);

//P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[0],y2[1],z2[0],x1[1],y1[0],z1[1]);

//Q = PQmap[ x2[1] ][ y2[1] ][ z2[0] ];

//P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[1],y2[1],z2[0],x1[0],y1[0],z1[1]);

//P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[1],y2[1],z2[0],x1[1],y1[0],z1[1]);

//Q = PQmap[ x2[0] ][ y2[1] ][ z2[0] ];

//P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[0],y2[1],z2[0],x1[0],y1[1],z1[1]);

//P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[0],y2[1],z2[0],x1[1],y1[1],z1[1]);

//Q = PQmap[ x2[1] ][ y2[1] ][ z2[0] ];

//P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[1],y2[1],z2[0],x1[0],y1[1],z1[1]);

//P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[1],y2[1],z2[0],x1[1],y1[1],z1[1]);

//Q = PQmap[ x2[0] ][ y2[0] ][ z2[1] ];

//P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[0],y2[0],z2[1],x1[0],y1[0],z1[0]);

//P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[0],y2[0],z2[1],x1[1],y1[0],z1[0]);

//Q = PQmap[ x2[1] ][ y2[0] ][ z2[1] ];

//P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[1],y2[0],z2[1],x1[0],y1[0],z1[0]);

//P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[1],y2[0],z2[1],x1[1],y1[0],z1[0]);

//Q = PQmap[ x2[0] ][ y2[0] ][ z2[1] ];

//P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[0],y2[0],z2[1],x1[0],y1[1],z1[0]);

//P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[0],y2[0],z2[1],x1[1],y1[1],z1[0]);

//Q = PQmap[ x2[1] ][ y2[0] ][ z2[1] ];

//P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[1],y2[0],z2[1],x1[0],y1[1],z1[0]);

//P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[1],y2[0],z2[1],x1[1],y1[1],z1[0]);

//Q = PQmap[ x2[0] ][ y2[1] ][ z2[1] ];

//P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[0],y2[1],z2[1],x1[0],y1[0],z1[0]);

//P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[0],y2[1],z2[1],x1[1],y1[0],z1[0]);

//Q = PQmap[ x2[1] ][ y2[1] ][ z2[1] ];

//P = PQmap[ x1[0] ][ y1[0] ][ z1[0] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[1],y2[1],z2[1],x1[0],y1[0],z1[0]);

//P = PQmap[ x1[1] ][ y1[0] ][ z1[0] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[1],y2[1],z2[1],x1[1],y1[0],z1[0]);

//Q = PQmap[ x2[0] ][ y2[1] ][ z2[1] ];

//P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[0],y2[1],z2[1],x1[0],y1[1],z1[0]);

//P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[0],y2[1],z2[1],x1[1],y1[1],z1[0]);

//Q = PQmap[ x2[1] ][ y2[1] ][ z2[1] ];

//P = PQmap[ x1[0] ][ y1[1] ][ z1[0] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[1],y2[1],z2[1],x1[0],y1[1],z1[0]);

//P = PQmap[ x1[1] ][ y1[1] ][ z1[0] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[1],y2[1],z2[1],x1[1],y1[1],z1[0]);

//Q = PQmap[ x2[0] ][ y2[0] ][ z2[1] ];

//P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[0],y2[0],z2[1],x1[0],y1[0],z1[1]);

//P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[0],y2[0],z2[1],x1[1],y1[0],z1[1]);

//Q = PQmap[ x2[1] ][ y2[0] ][ z2[1] ];

//P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[1],y2[0],z2[1],x1[0],y1[0],z1[1]);

//P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[1],y2[0],z2[1],x1[1],y1[0],z1[1]);

//Q = PQmap[ x2[0] ][ y2[0] ][ z2[1] ];

//P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[0],y2[0],z2[1],x1[0],y1[1],z1[1]);

//P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[0],y2[0],z2[1],x1[1],y1[1],z1[1]);

//Q = PQmap[ x2[1] ][ y2[0] ][ z2[1] ];

//P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[1],y2[0],z2[1],x1[0],y1[1],z1[1]);

//P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[1],y2[0],z2[1],x1[1],y1[1],z1[1]);

//Q = PQmap[ x2[0] ][ y2[1] ][ z2[1] ];

//P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[0],y2[1],z2[1],x1[0],y1[0],z1[1]);

//P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[0],y2[1],z2[1],x1[1],y1[0],z1[1]);

//Q = PQmap[ x2[1] ][ y2[1] ][ z2[1] ];

//P = PQmap[ x1[0] ][ y1[0] ][ z1[1] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[1],y2[1],z2[1],x1[0],y1[0],z1[1]);

//P = PQmap[ x1[1] ][ y1[0] ][ z1[1] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[1],y2[1],z2[1],x1[1],y1[0],z1[1]);

//Q = PQmap[ x2[0] ][ y2[1] ][ z2[1] ];

//P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[0],y2[1],z2[1],x1[0],y1[1],z1[1]);

//P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[0],y2[1],z2[1],x1[1],y1[1],z1[1]);

//Q = PQmap[ x2[1] ][ y2[1] ][ z2[1] ];

//P = PQmap[ x1[0] ][ y1[1] ][ z1[1] ];
  //eri_val -= PQ[P][Q];
  eri_val_test -= smallpq(x2[1],y2[1],z2[1],x1[0],y1[1],z1[1]);

//  P = PQmap[ x1[1] ][ y1[1] ][ z1[1] ];
  //eri_val += PQ[P][Q];
  eri_val_test += smallpq(x2[1],y2[1],z2[1],x1[1],y1[1],z1[1]);

  free(x1);
  free(x2);
  free(y1);
  free(y2);
  free(z1);
  free(z2);

  //if(eri_val != eri_val_test){
  //   printf("not the same\n");
  //   printf("true %f\n",eri_val);
  //   printf("mine %f\n",eri_val_test);
  ////   //printf("px: %d py: %d pz: %d qx: %d qy: %d qz: %d \n",x2[0],y2[0],z2[0],x1[0],y1[0],z1[0]);
  //}
  return eri_val_test;

}

// This function implements Eq. 4.7 and 4.8 in Peter Gills paper on 2-electrons in a cube
// Gauss-Legendre quadrature is used for the 3d integral on the range 0->1 for x, y, and z
// int dim is the number of points on this grid, double *xa is a vector containing the actual points on this grid, and
// double *w is a vector containing the weights associated with this grid (analogous to differential length elements
// in rectangle rule integration).
// double px, py, pz, qx, qy, qz has the same interpretation as it does in the Gill paper.
double JelliumIntegrals::pq_int(int dim, double *xa, double *w, int px, int py, int pz, int qx, int qy, int qz) {
  double sum = 0.;
  double num, denom, x, y, z, dx, dy, dz, gx, gy, gz;
  double pi = M_PI;
  if (px<0 || qx<0 || py<0 || qy<0 || pz<0 || qz<0) {
         return 0.;
  } else {
    for (int i=0; i<dim; i++) {
        x = xa[i];
        dx = w[i];
        gx = g_pq(px, qx, x);
        for (int j=0; j<dim; j++) {
            y = xa[j];
            dy = w[j];
            gy = g_pq(py, qy, y);
            for (int k=0; k<dim; k++) {
                z = xa[k];
                dz = w[k];
                gz = g_pq(pz, qz, z);
                num = gx*gy*gz;
                denom = sqrt(x*x+y*y+z*z);
                sum += (num/denom)*dx*dy*dz;
                 
        //printf("  sum %f  i %d  j %d  k %d denom %f x %f dx %f y %f dy %f z %f dz %f gx %f gy %f gz %f\n",sum, i, j, k,denom,x,dx,y,dy,z,dz,gx,gy,gz);
            }
        }
    }
    return (8./pi)*sum;
  }
}
double JelliumIntegrals::pq_int_new(int dim, int px, int py, int pz, int qx, int qy, int qz) {
    //if (px<0 || qx<0 || py<0 || qy<0 || pz<0 || qz<0){
    //    return 0.;
    //}
    double * s_p = sqrt_tensor->pointer();
    double * g_p = g_tensor->pointer();
    double sum = 0.0;
    double tmp_gxgygz = 0.0;
    //double tmp_gxgy = 0.0;
    for (int i = 0; i < dim; i++){
        double gx = g_p[i * 2 * nmax * 2 * nmax + px * 2 * nmax + qx];
        for (int j = 0; j < dim; j++){
            double gxgy = gx*g_p[j * 2 * nmax * 2 * nmax + py * 2 * nmax + qy];
            for (int k = 0; k < dim; k++){
                double gxgygz = g_p[k * 2 * nmax * 2 * nmax + pz * 2 * nmax + qz];
                tmp_gxgygz += gxgygz * s_p[i*dim*dim + j*dim + k];
                //printf("  sum %f  x %f  y %f  z %f\n",sum, x, y, z);
            }
            sum = tmp_gxgygz * gxgy + sum;
            tmp_gxgygz = 0.0;
        }
    }
    return 8.0 * sum / M_PI;
}

/* 
/ Function to Determine Energy Calculation
/ Take function and loop through n to keep all atomic orbital energy levels.
*/
void JelliumIntegrals::OrderPsis3D(int &norbs, double *E, int **MO) {

    int c, d, i, j, k, l, idx;
    double swap;
    int **N;
    int cond, Ecur;
    N = MAT_INT(2*(norbs+1)*(norbs+1)*(norbs+1),3);
  
    // Going to start with i and j=0 because that's the beginning
    // of the array... nx=i+1, ny=j+1
  
    for ( i=0; i<norbs; i++) {
      for ( j=0; j<norbs; j++) {
        for ( k=0; k<norbs; k++) {
  
          idx = i*norbs*norbs+j*norbs+k;
          // l is related to nx^2 + ny^2 + nz^2, aka, (i+1)^2 + (j+1)^2 (k+1)^2
          l = (i+1)*(i+1) + (j+1)*(j+1) + (k+1)*(k+1);
          E[idx] = l;
          // index is and energy is ...
          //outfile->Printf("  Index is %i and Energy[%i,%i,%i] is %i\n",idx,i+1,j+1,k+1,l);
          // element N[k][0] is nx = i+1
          N[l][0] = i+1;
          // element N[k][1] is ny = j+1
          N[l][1] = j+1;
          // element N[k][2] is nz = k+1
          N[l][2] = k+1;
  
        }
      }
    }
  
    for ( c = 0 ; c < ( norbs*norbs*norbs-1 ); c++){
      for (d = 0 ; d < norbs*norbs*norbs - c - 1; d++){
        if (E[d] > E[d+1]) /* For decreasing order use < */
        {
          swap       = E[d];
          E[d]   = E[d+1];
          E[d+1] = swap;
        }
      }
    }
    // print all energy values
    //for (int i=0; i<(norbs*norbs*norbs); i++) {
    //  outfile->Printf(" E[%i] is %f \n",i,E[i]);
    //}
    c=0;
    do {
      Ecur = E[c];
      i=0;
      do {
        i++;
        j=0;
        do {
          j++;
          k=0;
          do {
            k++;
            cond=Ecur-(i*i+j*j+k*k);
            if (cond==0) {
              MO[c][0] = i;
              MO[c][1] = j;
              MO[c][2] = k;
              c++;
            }
          }while( Ecur==E[c] && k<norbs);
        }while( Ecur==E[c] && j<norbs);
      }while (Ecur==E[c] && i<norbs);
    }while(c<norbs*norbs*norbs);
    for ( c = 0 ; c < ( norbs*norbs*norbs ); c++){
       E[c] = E[c];///(length*length);
    }
  
    // reset nmax to be actual maximum necessary to consider, given orbitalMax
  
    //outfile->Printf(" exit successful \n");
  
      for (i=0; i<(norbs*norbs*norbs); i++) {
       // outfile->Printf("  Psi( %i , %i, %i ) %i\n",MO[i][0],MO[i][1],MO[i][2],E[i]);
      }

    int new_nmax = 0;
    for (int i = 0; i < orbitalMax; i++) {
        if ( MO[i][0] > new_nmax ) new_nmax = MO[i][0];
        if ( MO[i][1] > new_nmax ) new_nmax = MO[i][1];
        if ( MO[i][2] > new_nmax ) new_nmax = MO[i][2];
    }
    //printf("%5i\n",nmax);
    //printf("%5i\n",new_nmax);
    //exit(0);
    norbs = new_nmax;
}

//create irreps
void JelliumIntegrals::Orderirrep(int &norbs, double *E, int **MO, int electrons) {
    int eee = 0;
    int eeo = 0;
    int eoe = 0;
    int eoo = 0;
    int oee = 0;
    int oeo = 0;
    int ooe = 0;
    int ooo = 0;
    for (int i = 0; i < orbitalMax; i++) {
        if ( MO[i][0]%2==0 ){
            if ( MO[i][1]%2==0 ){
                if( MO[i][2]%2==0 ){
                   eee++;
                } else {
                   eeo++;
                }
            }else{
                if( MO[i][2]%2==0 ){
                   eoe++;
                } else {
                   eoo++;
                }
            }
        } else {
            if ( MO[i][1]%2==0 ){
                if( MO[i][2]%2==0 ){
                   oee++;
                } else {
                   oeo++;
                }
            }else{
                if( MO[i][2]%2==0 ){
                   ooe++;
                } else {
                   ooo++;
                }
            } 
        }
    } 
    //printf("eee: %d \neeo: %d \neoe: %d \neoo: %d\n oee: %d\n oeo: %d\n ooe: %d\n ooo: %d\n",eee,eeo,eoe,eoo,oee,oeo,ooe,ooo);
    nirrep_ = 8;
    nsopi_ = (int*)malloc(nirrep_*sizeof(int));
    nsopi_[0] = eee;
    nsopi_[1] = eeo;
    nsopi_[2] = eoe;
    nsopi_[3] = eoo;
    nsopi_[4] = oee;
    nsopi_[5] = oeo;
    nsopi_[6] = ooe;
    nsopi_[7] = ooo;
    int tmp = 0;
    double tmp_energy = 0;
    double max_energy = E[electrons/2];
    int* tmp_swap = (int*)malloc(3*sizeof(int));
    for(int i = 0; i < orbitalMax; i++){
        if(MO[i][0]%2==0 && MO[i][1]%2==0 && MO[i][2]%2==0){
           tmp_swap[0]=MO[tmp][0];
           tmp_swap[1]=MO[tmp][1];
           tmp_swap[2]=MO[tmp][2];
           tmp_energy = E[tmp];
           MO[tmp][0]=MO[i][0];
           MO[tmp][1]=MO[i][1];
           MO[tmp][2]=MO[i][2];
           E[tmp] = E[i];
           tmp++;
           MO[i][0]=tmp_swap[0];
           MO[i][1]=tmp_swap[1];
           MO[i][2]=tmp_swap[2];
           E[i]=tmp_energy;
        }
   
    }
    for(int i = 0; i < orbitalMax; i++){
        if(MO[i][0]%2==0 && MO[i][1]%2==0 && MO[i][2]%2==1){
           tmp_swap[0]=MO[tmp][0];
           tmp_swap[1]=MO[tmp][1];
           tmp_swap[2]=MO[tmp][2];
           tmp_energy = E[tmp];
           MO[tmp][0]=MO[i][0];
           MO[tmp][1]=MO[i][1];
           MO[tmp][2]=MO[i][2];
           E[tmp] = E[i];
           tmp++;
           MO[i][0]=tmp_swap[0];
           MO[i][1]=tmp_swap[1];
           MO[i][2]=tmp_swap[2];
           E[i]=tmp_energy;
        }
   
    }
    for(int i = 0; i < orbitalMax; i++){
        if(MO[i][0]%2==0 && MO[i][1]%2==1 && MO[i][2]%2==0){
           tmp_swap[0]=MO[tmp][0];
           tmp_swap[1]=MO[tmp][1];
           tmp_swap[2]=MO[tmp][2];
           tmp_energy = E[tmp];
           MO[tmp][0]=MO[i][0];
           MO[tmp][1]=MO[i][1];
           MO[tmp][2]=MO[i][2];
           E[tmp] = E[i];
           tmp++;
           MO[i][0]=tmp_swap[0];
           MO[i][1]=tmp_swap[1];
           MO[i][2]=tmp_swap[2];
           E[i]=tmp_energy;
        }
   
    }
    for(int i = 0; i < orbitalMax; i++){
        if(MO[i][0]%2==0 && MO[i][1]%2==1 && MO[i][2]%2==1){
           tmp_swap[0]=MO[tmp][0];
           tmp_swap[1]=MO[tmp][1];
           tmp_swap[2]=MO[tmp][2];
           tmp_energy = E[tmp];
           MO[tmp][0]=MO[i][0];
           MO[tmp][1]=MO[i][1];
           MO[tmp][2]=MO[i][2];
           E[tmp] = E[i];
           tmp++;
           MO[i][0]=tmp_swap[0];
           MO[i][1]=tmp_swap[1];
           MO[i][2]=tmp_swap[2];
           E[i]=tmp_energy;
        }
   
    }
    for(int i = 0; i < orbitalMax; i++){
        if(MO[i][0]%2==1 && MO[i][1]%2==0 && MO[i][2]%2==0){
           tmp_swap[0]=MO[tmp][0];
           tmp_swap[1]=MO[tmp][1];
           tmp_swap[2]=MO[tmp][2];
           tmp_energy = E[tmp];
           MO[tmp][0]=MO[i][0];
           MO[tmp][1]=MO[i][1];
           MO[tmp][2]=MO[i][2];
           E[tmp] = E[i];
           tmp++;
           MO[i][0]=tmp_swap[0];
           MO[i][1]=tmp_swap[1];
           MO[i][2]=tmp_swap[2];
           E[i]=tmp_energy;
        }
   
    }
    for(int i = 0; i < orbitalMax; i++){
        if(MO[i][0]%2==1 && MO[i][1]%2==0 && MO[i][2]%2==1){
           tmp_swap[0]=MO[tmp][0];
           tmp_swap[1]=MO[tmp][1];
           tmp_swap[2]=MO[tmp][2];
           tmp_energy = E[tmp];
           MO[tmp][0]=MO[i][0];
           MO[tmp][1]=MO[i][1];
           MO[tmp][2]=MO[i][2];
           E[tmp] = E[i];
           tmp++;
           MO[i][0]=tmp_swap[0];
           MO[i][1]=tmp_swap[1];
           MO[i][2]=tmp_swap[2];
           E[i]=tmp_energy;
        }
   
    }
    for(int i = 0; i < orbitalMax; i++){
        if(MO[i][0]%2==1 && MO[i][1]%2==1 && MO[i][2]%2==0){
           tmp_swap[0]=MO[tmp][0];
           tmp_swap[1]=MO[tmp][1];
           tmp_swap[2]=MO[tmp][2];
           tmp_energy = E[tmp];
           MO[tmp][0]=MO[i][0];
           MO[tmp][1]=MO[i][1];
           MO[tmp][2]=MO[i][2];
           E[tmp] = E[i];
           tmp++;
           MO[i][0]=tmp_swap[0];
           MO[i][1]=tmp_swap[1];
           MO[i][2]=tmp_swap[2];
           E[i]=tmp_energy;
        }
   
    }
    doccpi_ = (int*)malloc(nirrep_*sizeof(int));
    int offsetJ = 0;
    int ecounter = 0;
    for(int i = 0; i < nirrep_; i++){
        doccpi_[i]=0;
    }
    for(int i = 0; i < nirrep_; i++){
        for(int j = 0; j < nsopi_[i]; j++){
            if(E[offsetJ + j] < max_energy && ecounter < electrons/2){
                doccpi_[i]++;
                ecounter++;
            }
        }
        offsetJ += nsopi_[i];
    }
    offsetJ = 0;
    for(int i = 0; i < nirrep_; i++){
        for(int j = 0; j < nsopi_[i]; j++){
            if(E[offsetJ + j] == max_energy && ecounter < electrons/2){
                doccpi_[i]++;
                ecounter++;
            }
        }
        offsetJ += nsopi_[i];
    }
    //printf("max energy%f\n",max_energy);
    for(int i = 0; i < nirrep_; i++){
        //printf("electrons[%d] with <= max energy %d\n",i,doccpi_[i]);
    }
    offsetJ = 0;
    for(int i = 0; i < nirrep_; i++){
       for(int j = 0; j < nsopi_[i]; j++){
          for(int k = j+1; k < nsopi_[i]; k++){
             if(E[j+offsetJ]>E[k+offsetJ]){
                tmp_swap[0] = MO[j+offsetJ][0];     
                tmp_swap[1] = MO[j+offsetJ][1];     
                tmp_swap[2] = MO[j+offsetJ][2];
                tmp_energy = E[j+offsetJ];
                MO[j+offsetJ][0] = MO[k+offsetJ][0];    
                MO[j+offsetJ][1] = MO[k+offsetJ][1];    
                MO[j+offsetJ][2] = MO[k+offsetJ][2];
                E[j+offsetJ] = E[k+offsetJ];
                MO[k+offsetJ][0] = tmp_swap[0];    
                MO[k+offsetJ][1] = tmp_swap[1];    
                MO[k+offsetJ][2] = tmp_swap[2]; 
                E[k+offsetJ] = tmp_energy; 
             }
          }
       }     
       offsetJ += nsopi_[i];
    }
    for(int i = 0; i < orbitalMax; i++){
       //printf("MO[%d][0]:%d\tMO[%d][1]:%d\tMO[%d][2]:%d energy: %f\n",i,MO[i][0],i,MO[i][1],i,MO[i][2],E[i]);
    }
}
/*
void JelliumIntegrals::Orderirrep(int &norbs, double *E, int **MO, int electrons) {
    int ee = 0;
    int eo = 0;
    int oe = 0;
    int oo = 0;
    for (int i = 0; i < orbitalMax; i++) {
            if ( MO[i][1]%2==0 ){
                if( MO[i][2]%2==0 ){
                   ee++;
                } else {
                   eo++;
                }
            }else{
                if( MO[i][2]%2==0 ){
                   oe++;
                } else {
                   oo++;
                }
            }
    } 
    //printf("eee: %d \neeo: %d \neoe: %d \neoo: %d\n oee: %d\n oeo: %d\n ooe: %d\n ooo: %d\n",eee,eeo,eoe,eoo,oee,oeo,ooe,ooo);
    nirrep_ = 4;
    nsopi_ = (int*)malloc(nirrep_*sizeof(int));
    nsopi_[0] = ee;
    nsopi_[1] = eo;
    nsopi_[2] = oe;
    nsopi_[3] = oo;
    int tmp = 0;
    double tmp_energy = 0;
    double max_energy = E[electrons/2];
    int* tmp_swap = (int*)malloc(3*sizeof(int));
    for(int i = 0; i < orbitalMax; i++){
        if(MO[i][1]%2==0 && MO[i][2]%2==0){
           tmp_swap[0]=MO[tmp][0];
           tmp_swap[1]=MO[tmp][1];
           tmp_swap[2]=MO[tmp][2];
           tmp_energy = E[tmp];
           MO[tmp][0]=MO[i][0];
           MO[tmp][1]=MO[i][1];
           MO[tmp][2]=MO[i][2];
           E[tmp] = E[i];
           tmp++;
           MO[i][0]=tmp_swap[0];
           MO[i][1]=tmp_swap[1];
           MO[i][2]=tmp_swap[2];
           E[i]=tmp_energy;
        }
   
    }
    for(int i = 0; i < orbitalMax; i++){
        if(MO[i][1]%2==0 && MO[i][2]%2==1){
           tmp_swap[0]=MO[tmp][0];
           tmp_swap[1]=MO[tmp][1];
           tmp_swap[2]=MO[tmp][2];
           tmp_energy = E[tmp];
           MO[tmp][0]=MO[i][0];
           MO[tmp][1]=MO[i][1];
           MO[tmp][2]=MO[i][2];
           E[tmp] = E[i];
           tmp++;
           MO[i][0]=tmp_swap[0];
           MO[i][1]=tmp_swap[1];
           MO[i][2]=tmp_swap[2];
           E[i]=tmp_energy;
        }
   
    }
    for(int i = 0; i < orbitalMax; i++){
        if(MO[i][1]%2==1 && MO[i][2]%2==0){
           tmp_swap[0]=MO[tmp][0];
           tmp_swap[1]=MO[tmp][1];
           tmp_swap[2]=MO[tmp][2];
           tmp_energy = E[tmp];
           MO[tmp][0]=MO[i][0];
           MO[tmp][1]=MO[i][1];
           MO[tmp][2]=MO[i][2];
           E[tmp] = E[i];
           tmp++;
           MO[i][0]=tmp_swap[0];
           MO[i][1]=tmp_swap[1];
           MO[i][2]=tmp_swap[2];
           E[i]=tmp_energy;
        }
   
    }
    for(int i = 0; i < orbitalMax; i++){
        if(MO[i][1]%2==1 && MO[i][2]%2==1){
           tmp_swap[0]=MO[tmp][0];
           tmp_swap[1]=MO[tmp][1];
           tmp_swap[2]=MO[tmp][2];
           tmp_energy = E[tmp];
           MO[tmp][0]=MO[i][0];
           MO[tmp][1]=MO[i][1];
           MO[tmp][2]=MO[i][2];
           E[tmp] = E[i];
           tmp++;
           MO[i][0]=tmp_swap[0];
           MO[i][1]=tmp_swap[1];
           MO[i][2]=tmp_swap[2];
           E[i]=tmp_energy;
        }
   
    }
    doccpi_ = (int*)malloc(nirrep_*sizeof(int));
    int offsetJ = 0;
    int ecounter = 0;
    for(int i = 0; i < nirrep_; i++){
        doccpi_[i]=0;
    }
    for(int i = 0; i < nirrep_; i++){
        for(int j = 0; j < nsopi_[i]; j++){
            if(E[offsetJ + j] < max_energy && ecounter < electrons/2){
                doccpi_[i]++;
                ecounter++;
            }
        }
        offsetJ += nsopi_[i];
    }
    offsetJ = 0;
    for(int i = 0; i < nirrep_; i++){
        for(int j = 0; j < nsopi_[i]; j++){
            if(E[offsetJ + j] == max_energy && ecounter < electrons/2){
                doccpi_[i]++;
                ecounter++;
            }
        }
        offsetJ += nsopi_[i];
    }
    //printf("max energy%f\n",max_energy);
    //for(int i = 0; i < nirrep_; i++){
    //    printf("electrons[%d] with <= max energy %d\n",i,doccpi_[i]);
    //}
    offsetJ = 0;
    for(int i = 0; i < nirrep_; i++){
       for(int j = 0; j < nsopi_[i]; j++){
          for(int k = j+1; k < nsopi_[i]; k++){
             if(E[j+offsetJ]>E[k+offsetJ]){
                tmp_swap[0] = MO[j+offsetJ][0];     
                tmp_swap[1] = MO[j+offsetJ][1];     
                tmp_swap[2] = MO[j+offsetJ][2];
                tmp_energy = E[j+offsetJ];
                MO[j+offsetJ][0] = MO[k+offsetJ][0];    
                MO[j+offsetJ][1] = MO[k+offsetJ][1];    
                MO[j+offsetJ][2] = MO[k+offsetJ][2];
                E[j+offsetJ] = E[k+offsetJ];
                MO[k+offsetJ][0] = tmp_swap[0];    
                MO[k+offsetJ][1] = tmp_swap[1];    
                MO[k+offsetJ][2] = tmp_swap[2]; 
                E[k+offsetJ] = tmp_energy; 
             }
          }
       }     
       offsetJ += nsopi_[i];
    }
    //for(int i = 0; i < orbitalMax; i++){
    //   printf("MO[%d][0]:%d\tMO[%d][1]:%d\tMO[%d][2]:%d energy: %f\n",i,MO[i][0],i,MO[i][1],i,MO[i][2],E[i]);
    //}
}
*/

int ** JelliumIntegrals::MAT_INT(int dim1, int dim2){
  int i,j;
  int **M;
  M = (int **)malloc(dim1*sizeof(int *));
  if (M==NULL) {
     outfile->Printf("\n\nMAT_INT: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim1; i++){
      M[i] = (int *)malloc(dim2*sizeof(int));
      if (M[i]==NULL) {
         outfile->Printf("\n\nMAT_INT: Memory allocation error\n\n");
         exit(0);
      }
      for (j=0; j<dim2; j++){
          M[i][j] = 0;
      }
  }
  return M;
}

int JelliumIntegrals::get_nmax(){
    return nmax;
}

//to get the 1 electron integral result from the small pq structure
double JelliumIntegrals::smallpq(int px, int py, int pz, int qx, int qy, int qz){
        int px_t = px;
        int py_t = py;
        int pz_t = pz;
        int qx_t = qx;
        int qy_t = qy;
        int qz_t = qz;
        if(qx_t < px_t){
           px_t = qx_t;
           qx_t = px;
        }
        if(qy_t < py_t){
           py_t = qy_t;
           qy_t = py;
        }
        if(qz_t < pz_t){
           pz_t = qz_t;
           qz_t = pz;
        }
        //int test = 0;
        int pqx = 0;
        //for(int i = 0; i < px_t; i++){
        //   for(int j = i; j <= 2*nmax; j++){
        //       if((i+j)%2==0){
        //       pqx++;
        //       }
        //   }
        //}
        //pqx = 0;
        //int tmp = 1;
        if(px_t%2==0){
           //pqx = (nmax*px_t+(1-px_t/2)*px_t/2)+px_t/2;
           pqx = (nmax*px_t+px_t-(px_t*px_t)/4);
        }
        else{
           pqx = (nmax*px_t)+(1-px_t/2)*(px_t-1)/2+1;
        }
        //for(int i = 0; i < px_t; i++){
        //    pqx += nmax+tmp;
        //    if(i%2==0){
        //       tmp--; 
        //    }
        //}
        //printf("pqx: %d 2*nmax: %d px_t: %d\n",pqx,2*nmax,px_t);
        pqx += (qx_t-px_t)/2;
        //int pqy = 0;
        //for(int i = 0; i < py_t; i++){
        //   for(int j = i; j <= 2*nmax; j++){
        //       if((i+j)%2==0){
        //       pqy++;
        //       }
        //   }
        //}
        //tmp = 1;
        int pqy = 0;
        if(py_t%2==0){
           //pqy = (nmax*py_t+(1-py_t/2)*py_t/2)+py_t/2;
           pqy = (nmax*py_t+py_t-(py_t*py_t)/4);
        }
        else{
           pqy = (nmax*py_t)+(1-py_t/2)*(py_t-1)/2+1;
        }
        //for(int i = 0; i < py_t; i++){
        //    pqy += nmax+tmp;
        //    if(i%2==0){
        //       tmp--; 
        //    }
        //}
        pqy += (qy_t-py_t)/2;
        int pqz = 0;
        //for(int i = 0; i < pz_t; i++){
        //   for(int j = i; j <= 2*nmax; j++){
        //       if((i+j)%2==0){
        //       pqz++;
        //       }
        //   }
        //}
        //tmp = 1;
        if(pz_t%2==0){
           //pqz = (nmax*pz_t+(1-pz_t/2)*pz_t/2)+pz_t/2;
           pqz = (nmax*pz_t+pz_t-(pz_t*pz_t)/4);
        }
        else{
           pqz = (nmax*pz_t)+(1-pz_t/2)*(pz_t-1)/2+1;
        }
        //for(int i = 0; i < pz_t; i++){
        //    pqz += nmax+tmp;
        //    if(i%2==0){
        //       tmp--; 
        //    }
        //}
        pqz += (qz_t-pz_t)/2;
        //test = pqx*offset_pq*offset_pq + pqy*offset_pq + pqz;
        return PQ_small[pqx*offset_pq*offset_pq + pqy*offset_pq + pqz];
}

} // end of namespaces
