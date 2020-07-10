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
#include <psi4/psifiles.h>
#include <psi4/libiwl/iwl.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libmints/mintshelper.h>

#include "doci_solver.h"

using namespace psi;

namespace psi{namespace doci{

struct tpdm {
    int i;
    int j;
    int k;
    int l;
    double val;
};

void DOCISolver::WriteTPDM(){

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

    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            for (int k = 0; k < nmo_; k++) {
                for (int l = 0; l < nmo_; l++) {

                    double valab = d2ab_[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];

                    tpdm d2;
                    d2.i   = i;
                    d2.j   = j;
                    d2.k   = k;
                    d2.l   = l;
                    d2.val = valab;
                    psio->write(PSIF_V2RDM_D2AB,"D2ab",(char*)&d2,sizeof(tpdm),addr_ab,&addr_ab);
                    countab++;

                    if ( i != j && k != l ) {

                        double valaa = d2aa_[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];
                        double valbb = d2aa_[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];

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

}

}} //end namespaces


