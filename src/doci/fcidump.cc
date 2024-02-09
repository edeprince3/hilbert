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

#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libmints/mintshelper.h>

#include <psi4/libtrans/integraltransform.h>
#include <psi4/libtrans/mospace.h>

#include <psi4/libmints/writer.h>
#include <psi4/libmints/writer_file_prefix.h>
#include <psi4/libmints/molecule.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/matrix.h>

#include "doci_solver.h"

#include <misc/omp.h>

using namespace psi;

namespace hilbert{


// 
// dump 1-/2-electron integrals and 1-/2-RDM to disk
// 
void DOCISolver::FCIDUMP() {

    std::string filename = get_writer_file_prefix(reference_wavefunction_->molecule()->name());
    std::string rdm_filename = filename + ".rdm";
    std::string int_filename = filename + ".int";

    FILE * int_fp = fopen(int_filename.c_str(),"w");
    FILE * rdm_fp = fopen(rdm_filename.c_str(),"w");

    fprintf(int_fp,"%5i\n",nmo_);
    fprintf(rdm_fp,"%5i\n",nmo_);

    int zero = 0;

    // two-electron integrals
    for (int p = 0; p < nmo_; p++) {
        for (int q = p; q < nmo_; q++) {
            long int pq = INDEX(p,q);
            for (int r = 0; r < nmo_; r++) {
                for (int s = r; s < nmo_; s++) {
                    long int rs = INDEX(r,s);
                    if ( pq > rs ) continue;
                    double dum = C_DDOT(nQ_,Qmo_ + pq,nmo_*(nmo_+1)/2,Qmo_+rs,nmo_*(nmo_+1)/2);
                    if ( fabs(dum) > 1e-12 ) {
                        fprintf(int_fp,"%20.12le %5i %5i %5i %5i\n",dum,p+1,q+1,r+1,s+1);
                    }
                }
            }
        }
    }

    // one-electron integrals

    std::shared_ptr<MintsHelper> mints(new MintsHelper(reference_wavefunction_));
    std::shared_ptr<Matrix> T (new Matrix(mints->so_kinetic()));
    std::shared_ptr<Matrix> V (new Matrix(mints->so_potential()));

    T->transform(Ca_);
    V->transform(Ca_);

    double ** Tp = T->pointer();
    double ** Vp = V->pointer();
    for (int p = 0; p < nmo_; p++) {
        for (int q = 0; q < nmo_; q++) {
            double dum = Tp[p][q] + Vp[p][q];
            if ( fabs(dum) > 1e-12 ) {
                fprintf(int_fp,"%20.12le %5i %5i %5i %5i\n",dum,p+1,q+1,0,0);
            }
        }
    }
    fprintf(int_fp,"%20.12le %5i %5i %5i %5i\n",enuc_,0,0,0,0);
    fclose(int_fp);

    // two-electron rdm

    size_t n1 = nmo_;
    size_t n2 = nmo_ * nmo_ ;
    size_t n3 = nmo_ * nmo_ * nmo_;

    for (size_t i = 0; i < nmo_; i++) {
        for (size_t j = 0; j < nmo_; j++) {
            for (size_t k = 0; k < nmo_; k++) {
                for (size_t l = 0; l < nmo_; l++) {

                    double dum = d2ab_[i * n3 + j * n2 + k * n1 + l] 
                               + d2ab_[j * n3 + i * n2 + l * n1 + k];

                    if ( i != j && k != l ) {

                        dum += 2.0 * d2aa_[i * n3 + j * n2 + k * n1 + l];

                    }

                    if ( fabs(dum) > 1e-12 ) {
                        fprintf(rdm_fp,"%20.12le %5zu %5zu %5zu %5zu\n",dum,i+1,j+1,k+1,l+1);
                    }
                }
            }
        }
    }

    // one-electron rdm
    for (int i = 0; i < nmo_; i++) {

        double dum = d1_[INDEX(i,i)];

        fprintf(rdm_fp,"%20.12le %5i %5i %5i %5i\n",dum,i+1,i+1,0,0);

    }

    fclose(rdm_fp);


}


}
