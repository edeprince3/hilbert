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

#include "ta_diis.h"

#include <psi4/libqt/qt.h>
#include <psi4/psi4-dec.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/psifiles.h>

#include <string.h>

using namespace psi;

namespace hilbert {

DIISTA::DIISTA(int n) : hist(), histerr() {

    maxdiis_           = n;
    diis_iter_         = 0;
    replace_diis_iter_ = 1;
    diisvec_           = (double*)malloc(sizeof(double)*(maxdiis_+1));
    errmtx             = (double*)malloc(sizeof(double)*maxdiis_*maxdiis_);
    memset((void*)errmtx,'\0',maxdiis_*maxdiis_*sizeof(double));


}
DIISTA::~DIISTA() {
    free(diisvec_);
    free(errmtx);
}

// reset diis solver (without freeing memory).
void DIISTA::restart() {
    diis_iter_         = 0;
    replace_diis_iter_ = 1;
    hist.clear();
    histerr.clear();
    memset((void*)errmtx,'\0',maxdiis_*maxdiis_*sizeof(double));
}

std::string DIISTA::get_index(const TA::TArrayD& array){
    size_t rank = array.trange().rank();

    std::string idx = "";
    for (int i = 0; i < rank-1; i++) {
        idx += "i" + std::to_string(i) + ",";
    }
    idx += "i" + std::to_string(rank-1);
    return idx;
}

// Store current solution vector to DIIS history
void DIISTA::WriteVector(const std::vector<TA::TArrayD*> &vector){
    std::vector<TA::TArrayD> copyVec;
    for (const auto &item: vector){
        TA::TArrayD temp;
        std::string idx = get_index(*item);
        temp(idx) = (*item)(idx);
        copyVec.push_back(temp);
    }
    if ( diis_iter_ <= maxdiis_ ){
        hist.push_back(copyVec);
    }
    else{
        hist[replace_diis_iter_] = copyVec;
    }
}

// Store current error vector to DIIS history
void DIISTA::WriteErrorVector(const std::vector<TA::TArrayD*>& vector){
    std::vector<TA::TArrayD> copyVec;
    for (const auto &item: vector){
        TA::TArrayD temp;
        std::string idx = get_index(*item);
        temp(idx) = (*item)(idx);
        copyVec.push_back(temp);
    }
    if ( diis_iter_ <= maxdiis_ ){
        histerr.push_back(copyVec);
    }
    else{
        histerr[replace_diis_iter_] = copyVec;
    }
}

// Perform DIIS extrapolation.
void DIISTA::Extrapolate(std::vector<TA::TArrayD*> &vector){

    if ( diis_iter_ > 1 ) {
        // Compute coefficients for the extrapolation
        DIISCoefficients( diis_iter_ < maxdiis_ ? diis_iter_ : maxdiis_ );
        int max = diis_iter_;
        if (max > maxdiis_) max = maxdiis_;

        //zero out amplitudes in vector to write in extrapolated amplitudes
        for (auto & vec : vector) {
            std::string idx = get_index(*vec);
            (*vec)(idx) = 0.0 * (*vec)(idx);
        }

        for (int j = 1; j <= max; j++){
            // Accumulate extrapolated vector.
            for (int k = 0; k < vector.size(); k++) {
                // adds next term in extrapolation to the amplitude.
                std::string idx = get_index(*vector[k]);
                (*(vector[k]))(idx) += diisvec_[j-1] * hist[j][k](idx);
            }
        }
    }

    if (diis_iter_ <= maxdiis_){
        diis_iter_++;
    }
    else {
        // If we already have maxdiis_ vectors, choose the one with
        // the largest error as the one to replace.
        int jmax   = 0;
        double max = -1.0e99;
        for (int j = 0; j < maxdiis_; j++){
            double nrm = diis_norm(j+1);
            if ( nrm > max ) {
                max  = nrm;
                jmax = j+1;
            }
        }
        replace_diis_iter_ = jmax;
    }
}

// Evaluate extrapolation coefficients for DIIS.
void DIISTA::DIISCoefficients(int nvec){

    // Allocate memory for small matrices/vectors.
    int * ipiv    = (int*)malloc((nvec+1)*sizeof(int));
    auto * A    = (double*)malloc(sizeof(double)*(nvec+1)*(nvec+1));
    auto * B    = (double*)malloc(sizeof(double)*(nvec+1));
    memset((void*)A,'\0',(nvec+1)*(nvec+1)*sizeof(double));
    memset((void*)B,'\0',(nvec+1)*sizeof(double));
    B[nvec] = -1.;

    // Reshape the error matrix, in case its dimension is less than maxdiis_.
    for (int i = 0; i < nvec; i++){
        for (int j = 0; j < nvec; j++){
            A[i*(nvec+1)+j] = errmtx[i*maxdiis_+j];
        }
    }

    if (nvec <= 3) {
        // At early iterations, just build the whole matrix.
        for (int i = 0; i < nvec; i++) {
            for (int j = i+1; j < nvec; j++){
                double sum  = diis_dot(i+1, j+1);
                A[i*(nvec+1)+j] = sum;
                A[j*(nvec+1)+i] = sum;
            }
            double sum  = diis_dot(i+1, i+1);
            A[i*(nvec+1)+i] = sum;
        }
    }else {
        // At later iterations, don't build the whote matrix.
        // Just replace one row/column.

        // Which row/column will be replaced?
        int i = nvec < maxdiis_ ? nvec - 1 : replace_diis_iter_-1;

        for (int j = 0; j < nvec; j++){
            double sum  = diis_dot(i+1, j+1);
            A[i*(nvec+1)+j] = sum;
            A[j*(nvec+1)+i] = sum;
        }
    }

    int j = nvec;
    for (int i = 0; i < (nvec+1); i++){
        A[j*(nvec+1)+i] = -1.0;
        A[i*(nvec+1)+j] = -1.0;
    }
    A[(nvec+1)*(nvec+1)-1] = 0.0;

    // Write error matrix for next iteration.
    //memset((void*)errmtx,'\0',maxdiis_*maxdiis_*sizeof(double));
    for (int i = 0; i < nvec; i++){
        for (int k = 0; k < nvec; k++){
            errmtx[i*maxdiis_ + k] = A[i * (nvec + 1) + k];
        }
    }

    // Solve the set of linear equations for the extrapolation coefficients
    int nrhs = 1;
    int lda  = nvec+1;
    int ldb  = nvec+1;
    C_DGESV(nvec+1,nrhs,A,lda,ipiv,B,ldb);
    C_DCOPY(nvec,B,1,diisvec_,1);

    free(A);
    free(B);
    free(ipiv);
}

double DIISTA::diis_dot(int i, int j){
    std::vector<TA::TArrayD> evec_i(histerr[i]);
    std::vector<TA::TArrayD> evec_j(histerr[j]);

    double sum = 0.0;
    //Loop over tiled array amplitudes
    for (int k = 0; k < evec_i.size(); k++) {
        //Compute dot product between the amplitude in corresponding error vectors
        std::string idx = get_index(evec_i[k]);
        sum += evec_i[k](idx).dot(evec_j[k](idx));
    }

    return sum;
}

double DIISTA::diis_norm(int i){
    std::vector<TA::TArrayD> evec_i(histerr[i]);

    double norm = 0.0;
    for (auto & evec : evec_i) {
        std::string idx = get_index(evec);
        norm += evec(idx).dot(evec(idx));
    }

    return sqrt(norm);
}

} // end of namespaces
