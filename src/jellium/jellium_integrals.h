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

#ifndef JELLIUM_INTEGRALS_H
#define JELLIUM_INTEGRALS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unordered_map>

#include <psi4/liboptions/liboptions.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>

using namespace psi;

namespace hilbert{ 

class JelliumIntegrals{

  public:
    JelliumIntegrals(Options & options);
    ~JelliumIntegrals();

    int get_nmax();
    int * nsopi(){ return nsopi_;}
    int * doccpi(){ return doccpi_;}
    int nirrep(){ return nirrep_;}

    double ERI(int a, int b, int c, int d);

    //Integral matricies
    std::shared_ptr<Matrix> NucAttrac;
    std::shared_ptr<Matrix> Ke;
    std::shared_ptr<Matrix> PQ;
    double selfval = 0.0;
    double* grid_points;
    std::shared_ptr<Vector> sqrt_tensor;
    std::shared_ptr<Vector> g_tensor;
    int ** MO;
    double* w;
    double pq_int_new(int dim, int px, int py, int pz, int qx, int qy, int qz);

  private:
    int nirrep_;
    int* nsopi_;
    int* doccpi_;
    int offset_pq;
    double ERI_unrolled_test(int * a, int * b, int * c, int * d, double ** PQ, int *** PQmap);

    //int * x1, * x2, * y1, * y2, * z1, * z2;
    double ERI_unrolled_new(int * a, int * b, int * c, int * d, double ** PQ, int *** PQmap);
    bool symmetry;
    int get_pq(int px, int py, int pz, int qx, int qy, int qz);
    double * PQ_small;
    double smallpq(int px, int py, int pz, int qx, int qy, int qz);
    double * px1;
    void compute_integrals();
    void Orderirrep(int &norbs, double *E, int **MO, int electrons);
    double n_order_;
    int electrons;
    int *** PQmap;
    int orbitalMax;
    int nmax = 0;
    double pi;
    /// Options object
    Options & options_;
    //  Electron integral functions
    double g_pq(int p, int q, double r);
    double pq_int(int dim, double *x, double *w, int px, int py, int pz, int qx, int qy, int qz);
    double E0_Int(int dim, double *xa, double *w);
    double Vab_Int(int dim, double *xa, double *w, int *a, int *b);
    double Vab_Int_new(int dim, double *xa, double *w, int *a, int *b);

    double ERI_unrolled(int * a, int * b, int * c, int * d);
    
    void OrderPsis3D(int &norbs, double *E, int **MO);
    int **MAT_INT(int dim1, int dim2);
};

}

#endif
