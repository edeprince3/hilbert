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

#include <psi4/libpsi4util/process.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/basisset.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libmints/sieve.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>
#include <psi4/psifiles.h>
#include <psi4/libtrans/integraltransform.h>
#include "threeindexintegralsta.h"

#include "qed_blas.h"
#include "../../misc/hilbert_psifiles.h"

using namespace psi;
using namespace fnocc;

namespace hilbert {

double* ThreeIndexIntegrals(std::shared_ptr<Wavefunction> ref, size_t &nQ, long int memory) {

    int nmo    = ref->nmo();
    int nso    = ref->nso();
    int nirrep = ref->nirrep();


    std::shared_ptr<BasisSet> basisset = ref->basisset();

    // get ntri from sieve
    std::shared_ptr<ERISieve> sieve (new ERISieve(basisset, ref->options().get_double("INTS_TOLERANCE")));
    const std::vector<std::pair<int, int> >& function_pairs = sieve->function_pairs();
    long int ntri = function_pairs.size();

    // read integrals that were written to disk in the scf
    std::shared_ptr<PSIO> psio(new PSIO());

    if ( ref->options().get_str("SCF_TYPE") == "CD" ) {
        long int nQ_tmp;
        psio->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DFSCF_BJ, "length", (char*)&nQ_tmp, sizeof(long int));
        psio->close(PSIF_DFSCF_BJ,1);
        nQ = (size_t) nQ_tmp;
    } else {
        throw PsiException("This file should not be called if scf_type != CD",__FILE__,__LINE__);
    }

    // 100 mb extra to account for all mapping arrays already 
    // allocated. this should be WAY more than necessary.
    long int extra = 100 * 1024 * 1024;  
    long int ndoubles = (memory-extra) / 8;

    // orbitals will end up in energy order.  
    // we will want them in pitzer.  for sorting: 
    long int * reorder  = (long int*)malloc(nmo*sizeof(long int));
    long int * sym      = (long int*)malloc(nmo*sizeof(long int));
    bool * skip    = (bool*)malloc(nmo*sizeof(bool));

    for (long int i = 0; i < nmo; i++) {
        skip[i] = false;
    }
    for (long int i = 0; i < nmo; i++) {
        double min   = 1.e99;
        long int count    = 0;
        long int minj     = -999;
        long int minh     = -999;
        long int mincount = -999;
        for (int h = 0; h < nirrep; h++) {
            for (long int j = 0; j < ref->nmopi()[h]; j++) {
                if ( skip[count+j] ) continue;
                if ( ref->epsilon_a()->pointer(h)[j] < min ) {
                    min      = ref->epsilon_a()->pointer(h)[j];
                    mincount = count;
                    minj     = j;
                    minh     = h;
                }
            }
            count += ref->nmopi()[h];
        }
        skip[mincount + minj]     = true;
        reorder[i]                = minj;
        sym[i]                    = minh;
    }

    // how many rows of (Q|mn) can we read in at once?
    if ( ndoubles < nso*nso ) {
        throw PsiException("holy moses, we can't fit nso^2 doubles in memory.  increase memory!",__FILE__,__LINE__);
    }

    long int nrows = 1;
    long int rowsize = nQ;
    while ( rowsize*nso*nso*2 > ndoubles ) {
        nrows++;
        rowsize = nQ / nrows;
        if (nrows * rowsize < nQ) rowsize++;
        if (rowsize == 1) break;
    }
    long int lastrowsize = nQ - (nrows - 1L) * rowsize;
    long int * rowdims = new long int [nrows];
    for (int i = 0; i < nrows-1; i++) rowdims[i] = rowsize;
    rowdims[nrows-1] = lastrowsize;

    double * tmp1 = (double*)malloc(rowdims[0]*nso*nso*sizeof(double));
    double * tmp2 = (double*)malloc(rowdims[0]*nso*nso*sizeof(double));

    // read integrals from SCF and unpack them
    psio_address addr  = PSIO_ZERO;
    psio->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
    psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);

    memset((void*)tmp1,'\0',nso*nso*rowdims[0]*sizeof(double));
    for (long int row = 0; row < nrows; row++) {

        // read
        psio->read(PSIF_DFSCF_BJ, "(Q|mn) Integrals", (char*) tmp2, sizeof(double) * ntri * rowdims[row],addr,&addr);

        // unpack
        #pragma omp parallel for schedule (static)
        for (long int Q = 0; Q < rowdims[row]; Q++) {
            for (long int mn = 0; mn < ntri; mn++) {

                long int m = function_pairs[mn].first;
                long int n = function_pairs[mn].second;

                tmp1[Q*nso*nso+m*nso+n] = tmp2[Q*ntri+mn];
                tmp1[Q*nso*nso+n*nso+m] = tmp2[Q*ntri+mn];
            }
        }
    }
    psio->close(PSIF_DFSCF_BJ,1);

    delete[] rowdims;
    free(reorder);
    free(skip);
    free(sym);
    free(tmp2);
    return tmp1;

}


}
