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

#include "doci_solver.h"

#include <psi4/libpsi4util/process.h>
#include <psi4/libmints/basisset.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libmints/sieve.h>
#include <psi4/libmints/matrix.h>
#include <psi4/psifiles.h>
#include <psi4/libtrans/integraltransform.h>

#include "blas.h"

using namespace psi;
using namespace fnocc;

namespace psi { namespace doci {

void DOCISolver::ThreeIndexIntegrals() {

    basisset_ = reference_wavefunction_->basisset();

    // get ntri from sieve
    std::shared_ptr<ERISieve> sieve (new ERISieve(basisset_, options_.get_double("INTS_TOLERANCE")));
    const std::vector<std::pair<int, int> >& function_pairs = sieve->function_pairs();
    long int ntri = function_pairs.size();

    // read integrals that were written to disk in the scf
    //nQ_ = Process::environment.globals["NAUX (SCF)"];

    std::shared_ptr<PSIO> psio(new PSIO());
    //psio->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
    //psio->read_entry(PSIF_DFSCF_BJ, "length", (char*)&nQ_, sizeof(long int));
    //psio->close(PSIF_DFSCF_BJ,1);

    if ( options_.get_str("SCF_TYPE") == "DF" ) {
        std::shared_ptr<BasisSet> primary = reference_wavefunction_->basisset(); 
        std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");

        nQ_ = auxiliary->nbf();
        Process::environment.globals["NAUX (SCF)"] = nQ_;
    }else if ( options_.get_str("SCF_TYPE") == "CD" ) {
        psio->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DFSCF_BJ, "length", (char*)&nQ_, sizeof(long int));
        psio->close(PSIF_DFSCF_BJ,1);
    }

    // 100 mb extra to account for all mapping arrays already 
    // allocated. this should be WAY more than necessary.
    long int extra = 100 * 1024 * 1024;  
    long int ndoubles = (memory_-extra) / 8;

    // orbitals will end up in energy order.  
    // we will want them in pitzer.  for sorting: 
    long int * reorder  = (long int*)malloc(nmo_*sizeof(long int));
    long int * sym      = (long int*)malloc(nmo_*sizeof(long int));
    bool * skip    = (bool*)malloc(nmo_*sizeof(bool));

    for (long int i = 0; i < nmo_; i++) {
        skip[i] = false;
    }
    for (long int i = 0; i < nmo_; i++) {
        double min   = 1.e99;
        long int count    = 0;
        long int minj     = -999;
        long int minh     = -999;
        long int mincount = -999;
        for (int h = 0; h < nirrep_; h++) {
            for (long int j = 0; j < nmopi_[h]; j++) {
                if ( skip[count+j] ) continue;
                if ( epsilon_a_->pointer(h)[j] < min ) {
                    min      = epsilon_a_->pointer(h)[j];
                    mincount = count;
                    minj     = j;
                    minh     = h;
                }
            }
            count += nmopi_[h];
        }
        skip[mincount + minj]     = true;
        reorder[i]                = minj;
        sym[i]                    = minh;
    }

    // how many rows of (Q|mn) can we read in at once?
    if ( ndoubles < nso_*nso_ ) {
        throw PsiException("holy moses, we can't fit nso^2 doubles in memory.  increase memory!",__FILE__,__LINE__);
    }

    long int nrows = 1;
    long int rowsize = nQ_;
    while ( rowsize*nso_*nso_*2 > ndoubles ) {
        nrows++;
        rowsize = nQ_ / nrows;
        if (nrows * rowsize < nQ_) rowsize++;
        if (rowsize == 1) break;
    }
    long int lastrowsize = nQ_ - (nrows - 1L) * rowsize;
    long int * rowdims = new long int [nrows];
    for (int i = 0; i < nrows-1; i++) rowdims[i] = rowsize;
    rowdims[nrows-1] = lastrowsize;

    double * tmp1 = (double*)malloc(rowdims[0]*nso_*nso_*sizeof(double));
    double * tmp2 = (double*)malloc(rowdims[0]*nso_*nso_*sizeof(double));

    long int nn1mo = nmo_*(nmo_+1)/2;

    psio->open(PSIF_DCC_QSO,PSIO_OPEN_NEW);
    psio->open(PSIF_DCC_QMO,PSIO_OPEN_NEW);
    psio_address addr  = PSIO_ZERO;
    psio_address addr2 = PSIO_ZERO;
    for (long int row = 0; row < nrows; row++) {
        psio->write(PSIF_DCC_QSO, "(Q|mn) Integrals", (char*) tmp1, sizeof(double) * rowdims[row] * nso_ * nso_,addr,&addr);
        psio->write(PSIF_DCC_QMO, "(Q|mn) Integrals", (char*) tmp1, sizeof(double) * rowdims[row] * nn1mo,addr2,&addr2);
    }
    psio->close(PSIF_DCC_QSO,1);
    psio->close(PSIF_DCC_QMO,1);

    // read integrals from SCF and unpack them
    addr  = PSIO_ZERO;
    addr2 = PSIO_ZERO;
    psio->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
    psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);

    memset((void*)tmp1,'\0',nso_*nso_*rowdims[0]*sizeof(double));
    for (long int row = 0; row < nrows; row++) {

        // read
        psio->read(PSIF_DFSCF_BJ, "(Q|mn) Integrals", (char*) tmp2, sizeof(double) * ntri * rowdims[row],addr,&addr);

        // unpack
        #pragma omp parallel for schedule (static)
        for (long int Q = 0; Q < rowdims[row]; Q++) {
            for (long int mn = 0; mn < ntri; mn++) {

                long int m = function_pairs[mn].first;
                long int n = function_pairs[mn].second;

                tmp1[Q*nso_*nso_+m*nso_+n] = tmp2[Q*ntri+mn];
                tmp1[Q*nso_*nso_+n*nso_+m] = tmp2[Q*ntri+mn];
            }
        }

        // write
        psio->write(PSIF_DCC_QSO, "(Q|mn) Integrals", (char*) tmp1, sizeof(double) * nso_*nso_ * rowdims[row],addr2,&addr2);
    }
    psio->close(PSIF_DFSCF_BJ,1);

    // AO->MO transformation matrix:
    SharedMatrix myCa (new Matrix(reference_wavefunction_->Ca_subset("AO","ALL")));

    // transform first index:
    addr  = PSIO_ZERO;
    addr2 = PSIO_ZERO;
    for (long int row = 0; row < nrows; row++) {
        // read
        psio->read(PSIF_DCC_QSO, "(Q|mn) Integrals", (char*) tmp1, sizeof(double) * nso_*nso_ * rowdims[row],addr,&addr);

        // transform first index:
        F_DGEMM('n','n',nmo_,nso_*rowdims[row],nso_,1.0,&(myCa->pointer()[0][0]),nmo_,tmp1,nso_,0.0,tmp2,nmo_);

        // sort
        #pragma omp parallel for schedule (static)
        for (long int Q = 0; Q < rowdims[row]; Q++) {
            for (long int i = 0; i < nmo_; i++) {
                for (long int m = 0; m < nso_; m++) {
                    tmp1[Q*nso_*nmo_+i*nso_+m] = tmp2[Q*nso_*nmo_+m*nmo_+i];
                }
            }
        }

        // write
        psio->write(PSIF_DCC_QSO, "(Q|mn) Integrals", (char*) tmp1, sizeof(double) * nso_*nmo_ * rowdims[row],addr2,&addr2);
    }
    // transform second index:
    addr  = PSIO_ZERO;
    addr2 = PSIO_ZERO;
    psio->open(PSIF_DCC_QMO,PSIO_OPEN_OLD);
    for (long int row = 0; row < nrows; row++) {
        // read
        psio->read(PSIF_DCC_QSO, "(Q|mn) Integrals", (char*) tmp1, sizeof(double) * nso_*nmo_ * rowdims[row],addr,&addr);

        // transform second index:
        F_DGEMM('n','n',nmo_,nmo_*rowdims[row],nso_,1.0,&(myCa->pointer()[0][0]),nmo_,tmp1,nso_,0.0,tmp2,nmo_);

        // sort orbitals into pitzer order
        #pragma omp parallel for schedule (static)
        for (long int Q = 0; Q < rowdims[row]; Q++) {
            for (long int m = 0; m < nmo_; m++) {
                int hm = sym[m];
                long int offm = 0;
                for (int h = 0; h < hm; h++) {
                    offm += nmopi_[h] - frzvpi_[h];
                }
                if ( reorder[m] >= nmopi_[hm] - frzvpi_[hm] ) continue;
                long int mm = reorder[m] + offm;
                for (long int n = 0; n < nmo_; n++) {
                    int hn = sym[n];
                    long int offn = 0;
                    for (int h = 0; h < hn; h++) {
                        offn += nmopi_[h] - frzvpi_[h];
                    }
                    if ( reorder[n] >= nmopi_[hn] - frzvpi_[hn] ) continue;
                    long int nn = reorder[n] + offn;
                    tmp1[Q*nn1mo+INDEX(mm,nn)] = tmp2[Q*nmo_*nmo_+m*nmo_+n];
                }
            }
        }

        // write
        psio->write(PSIF_DCC_QMO, "(Q|mn) Integrals", (char*) tmp1, sizeof(double) * nn1mo * rowdims[row],addr2,&addr2);
    }
    psio->close(PSIF_DCC_QMO,1);
    psio->close(PSIF_DCC_QSO,1);

    delete[] rowdims;

    //F_DGEMM('t','t',nso_*nQ_,nso_,nso_,1.0,tmp1,nso_,&(myCa->pointer()[0][0]),nso_,0.0,tmp2,nso_*nQ_);
    //F_DGEMM('t','t',nso_*nQ_,nso_,nso_,1.0,tmp2,nso_,&(myCa->pointer()[0][0]),nso_,0.0,tmp1,nso_*nQ_);

    free(reorder);
    free(skip);
    free(sym);
    free(tmp2);
    free(tmp1);

    Qmo_ = (double*)malloc(nn1mo*nQ_*sizeof(double));
    memset((void*)Qmo_,'\0',nn1mo*nQ_*sizeof(double));
    psio->open(PSIF_DCC_QMO,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_QMO,"(Q|mn) Integrals",(char*)Qmo_,sizeof(double)*nQ_ * nn1mo);
    psio->close(PSIF_DCC_QMO,1);

    // the following code would transpose the three-index integrals (Q|mn) -> (mn|Q)

/*
    // with 3-index integrals in memory, how much else can we hold?
    ndoubles -= nn1fv*nQ_;

    nrows = 1;
    rowsize = nQ_;
    while ( rowsize*nn1fv > ndoubles ) {
        nrows++;
        rowsize = nQ_ / nrows;
        if (nrows * rowsize < nQ_) rowsize++;
        if (rowsize == 1) break;
    }

    lastrowsize = nQ_ - (nrows - 1L) * rowsize;
    long int * rowdims2 = new long int [nrows];
    for (long int i = 0; i < nrows-1; i++) rowdims2[i] = rowsize;
    rowdims2[nrows-1] = lastrowsize;

    double * tmp3 = (double*)malloc(rowdims2[0] * nn1fv*sizeof(double));
    memset((void*)tmp3,'\0',rowdims2[0] * nn1fv*sizeof(double));
    addr = PSIO_ZERO;
    psio->open(PSIF_DCC_QMO,PSIO_OPEN_OLD);

    long int totalQ = 0;
    for (long int row = 0; row < nrows; row++) {
        psio->read(PSIF_DCC_QMO,"(Q|mn) Integrals",(char*)tmp3,sizeof(double)*rowdims2[row] * nn1fv,addr,&addr);
        for (long int Q = 0; Q < rowdims2[row]; Q++) {
            for (long int mn = 0; mn < nn1fv; mn++) {
                Qmo_[mn*nQ_+totalQ] = tmp3[Q*nn1fv+mn];
            }
            totalQ++;
        }
    }
    psio->close(PSIF_DCC_QMO,1);

    delete rowdims2;
    free(tmp3);
*/
}


}}
