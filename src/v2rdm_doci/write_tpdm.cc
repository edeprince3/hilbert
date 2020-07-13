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
#include <psi4/psifiles.h>
#include <psi4/libiwl/iwl.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libmints/mintshelper.h>

#include "v2rdm_solver.h"

using namespace psi;

namespace hilbert{

struct tpdm {
    int i;
    int j;
    int k;
    int l;
    double val;
};

void v2RDMSolver::WriteTPDM(){

    memset((void*)d2_act_spatial_sym_,'\0',d2_act_spatial_dim_*sizeof(double));

    double * D2aa = (double*)malloc(gems_aa[0]*gems_aa[0]*sizeof(double));
    double * D2ab = (double*)malloc(gems_ab[0]*gems_ab[0]*sizeof(double));
    memset((void*)D2aa,'\0',gems_aa[0]*gems_aa[0]*sizeof(double));
    memset((void*)D2ab,'\0',gems_ab[0]*gems_ab[0]*sizeof(double));

    // d2

    double * x_p = x->pointer();

    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            D2aa[ij*gems_aa[h] + ij] = x_p[d2s2off_ + i*amo_ + j]; //D2_2_[i*nmo_+j];
        }
    }

    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];
            if ( i == j ) continue;
            D2ab[ij*gems_ab[h] + ij]  = x_p[d2s2off_ + i*amo_ + j]; //D2_2_[i*nmo_+j];
        }
    }

    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            int ii = ibas_ab_sym[0][i][i];
            int jj = ibas_ab_sym[0][j][j];
            D2ab[ii*gems_ab[0] + jj] = x_p[d2s0off_ + i*amo_ + j]; //D2_0_[i*nmo_+j];
        }
    }

    std::shared_ptr<PSIO> psio (new PSIO());

    psio->open(PSIF_V2RDM_D2AA,PSIO_OPEN_NEW);
    psio->open(PSIF_V2RDM_D2BB,PSIO_OPEN_NEW);
    psio->open(PSIF_V2RDM_D2AB,PSIO_OPEN_NEW);

    psio_address addr_aa = PSIO_ZERO;
    psio_address addr_bb = PSIO_ZERO;
    psio_address addr_ab = PSIO_ZERO;

    long int countaa = 0;
    long int countbb = 0;
    long int countab = 0;

    // active-active part

    for (int h = 0; h < nirrep_; h++) {

        for (int ij = 0; ij < gems_ab[h]; ij++) {

            int i     = bas_ab_sym[h][ij][0];
            int j     = bas_ab_sym[h][ij][1];
            int ifull = i;
            int jfull = j;

            for (int kl = 0; kl < gems_ab[h]; kl++) {

                int k     = bas_ab_sym[h][kl][0];
                int l     = bas_ab_sym[h][kl][1];
                int kfull = k;
                int lfull = l;

                double valab = D2ab[ij*gems_ab[h] + kl];

                tpdm d2;
                d2.i   = ifull;
                d2.j   = jfull;
                d2.k   = kfull;
                d2.l   = lfull;
                d2.val = valab;
                psio->write(PSIF_V2RDM_D2AB,"D2ab",(char*)&d2,sizeof(tpdm),addr_ab,&addr_ab);
                countab++;

                if ( i != j && k != l ) {

                    int ija = ibas_aa_sym[h][i][j];
                    int kla = ibas_aa_sym[h][k][l];

                    int sij = i < j ? 1 : -1;
                    int skl = k < l ? 1 : -1;

                    double valaa = sij * skl * D2aa[ija*gems_aa[h] + kla];
                    double valbb = sij * skl * D2aa[ija*gems_aa[h] + kla];

                    d2.val = valaa;
                    psio->write(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);
                    countaa++;

                    d2.val = valbb;
                    psio->write(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);
                    countbb++;

                }
            }
        }
    }

    // write the number of entries in each file
    psio->write_entry(PSIF_V2RDM_D2AA,"length",(char*)&countaa,sizeof(long int));
    psio->write_entry(PSIF_V2RDM_D2BB,"length",(char*)&countbb,sizeof(long int));
    psio->write_entry(PSIF_V2RDM_D2AB,"length",(char*)&countab,sizeof(long int));

    // it might be nice for post CASSCF codes to know what orbitals are active:

    psio->write_entry(PSIF_V2RDM_D2AA,"NUMBER ACTIVE ORBITALS",(char*)&nmo_,sizeof(int));
    psio->write_entry(PSIF_V2RDM_D2BB,"NUMBER ACTIVE ORBITALS",(char*)&nmo_,sizeof(int));
    psio->write_entry(PSIF_V2RDM_D2AB,"NUMBER ACTIVE ORBITALS",(char*)&nmo_,sizeof(int));

    addr_aa = PSIO_ZERO;
    addr_bb = PSIO_ZERO;
    addr_ab = PSIO_ZERO;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < nmopi_[h]; i++) {
            int ifull = i;
            psio->write(PSIF_V2RDM_D2AA,"ACTIVE ORBITALS",(char*)&ifull,sizeof(int),addr_aa,&addr_aa);
            psio->write(PSIF_V2RDM_D2BB,"ACTIVE ORBITALS",(char*)&ifull,sizeof(int),addr_bb,&addr_bb);
            psio->write(PSIF_V2RDM_D2AB,"ACTIVE ORBITALS",(char*)&ifull,sizeof(int),addr_ab,&addr_ab);
        }
    }
    
    // it might be nice for post CASSCF codes to know what orbitals are inactive:

    int inact = 0;
    psio->write_entry(PSIF_V2RDM_D2AA,"NUMBER INACTIVE ORBITALS",(char*)&inact,sizeof(int));
    psio->write_entry(PSIF_V2RDM_D2BB,"NUMBER INACTIVE ORBITALS",(char*)&inact,sizeof(int));
    psio->write_entry(PSIF_V2RDM_D2AB,"NUMBER INACTIVE ORBITALS",(char*)&inact,sizeof(int));

    // close files
    psio->close(PSIF_V2RDM_D2AA,1);
    psio->close(PSIF_V2RDM_D2BB,1);
    psio->close(PSIF_V2RDM_D2AB,1);

    free(D2aa);
    free(D2ab);

}

} //end namespaces


