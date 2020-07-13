/*
 *@BEGIN LICENSE
 *
 * DOCI, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * Copyright (c) 2014, The Florida State University. All rights reserved.
 * 
 *@END LICENSE
 *
 */

#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libqt/qt.h>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libtrans/mospace.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/matrix.h>

using namespace psi;

namespace psi{

// update Ca/Cb matrices and repack energy-order transformation matrix as pitzer order
void UpdateTransformationMatrix(std::shared_ptr<Wavefunction> ref, std::shared_ptr<Matrix> newMO, 
        std::shared_ptr<Matrix> Ca, std::shared_ptr<Matrix> Cb, double * orbopt_transformation_matrix) {

    int nmo    = ref->nmo();
    int nirrep = ref->nirrep();

    SharedMatrix temp ( new Matrix(newMO) );
    // repack energy-order transformation matrix in pitzer order
    for (int ieo = 0; ieo < nmo; ieo++) {

        int ifull = ieo;
        int hi    = 0;
        int i     = ifull;

        double ** tp = temp->pointer(hi);
        double ** np = newMO->pointer(hi);
            
        for (int jeo = 0; jeo < nmo; jeo++) {

            int jfull = jeo;
            int hj    = 0;
            int j     = jfull;

            if ( hi != hj ) continue;

            double dum = 0.0;

            for (int keo = 0; keo < nmo; keo++) {

                int kfull = keo;
                int hk    = 0;
                int k     = kfull;

                if ( hi != hk ) continue;

                dum += np[k][j] * orbopt_transformation_matrix[ieo * nmo + keo];
            }
            
            //tp[i][j] = dum;
            tp[i][j] = orbopt_transformation_matrix[ieo * nmo + jeo];
            
        }
    }

    // reset energy-order transformation matrix:
    memset((void*)orbopt_transformation_matrix,'\0',nmo*nmo*sizeof(double));
    for (int i = 0; i < nmo; i++) {
        orbopt_transformation_matrix[i*nmo+i] = 1.0;
    }

    // update so/mo coefficient matrix (only need Ca):
    for (int h = 0; h < nirrep; h++) {

        double ** tp = temp->pointer(h);
        double **cap = Ca->pointer(h);
        double **cbp = Cb->pointer(h);

        for (int mu = 0; mu < ref->nsopi()[h]; mu++) {

            double * temp3 = (double*)malloc(ref->nmopi()[h] * sizeof(double));

            for (int i = 0; i < ref->nmopi()[h]; i++) {
                double dum = 0.0;
                for (int j = 0; j < ref->nmopi()[h]; j++) {
                    dum += cap[mu][j] * tp[i][j];
                }
                temp3[i] = dum;
            }
            for (int i = 0; i < ref->nmopi()[h]; i++) {
                cap[mu][i] = temp3[i];
                cbp[mu][i] = temp3[i];
            }
            free(temp3);
        }
    }
    //newMO->print();

    SharedMatrix temp2 (new Matrix(temp));
    temp2->gemm(false,false,1.0,newMO,temp,0.0);
    newMO->copy(temp2);

}

}
