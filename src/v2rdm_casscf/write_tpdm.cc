/*
 *@BEGIN LICENSE
 *
 * v2RDM-CASSCF, a plugin to:
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

#include "v2rdm_solver.h"

using namespace psi;

namespace hilbert{

void v2RDMSolver::WriteTPDM(){

    double * x_p = x->pointer();

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
            int ifull = full_basis[i];
            int jfull = full_basis[j];

            for (int kl = 0; kl < gems_ab[h]; kl++) {

                int k     = bas_ab_sym[h][kl][0];
                int l     = bas_ab_sym[h][kl][1];
                int kfull = full_basis[k];
                int lfull = full_basis[l];

                double valab = x_p[d2aboff[h] + ij*gems_ab[h] + kl];

                tpdm d2;
                d2.i   = ifull;
                d2.j   = jfull;
                d2.k   = kfull;
                d2.l   = lfull;
                d2.value = valab;
                psio->write(PSIF_V2RDM_D2AB,"D2ab",(char*)&d2,sizeof(tpdm),addr_ab,&addr_ab);
                countab++;

                if ( i != j && k != l ) {

                    int ija = ibas_aa_sym[h][i][j];
                    int kla = ibas_aa_sym[h][k][l];

                    int sij = i < j ? 1 : -1;
                    int skl = k < l ? 1 : -1;

                    double valaa = sij * skl * x_p[d2aaoff[h] + ija*gems_aa[h] + kla];
                    double valbb = sij * skl * x_p[d2bboff[h] + ija*gems_aa[h] + kla];

                    d2.value = valaa;
                    psio->write(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);
                    countaa++;

                    d2.value = valbb;
                    psio->write(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);
                    countbb++;

                }
            }
        }
    }

    // core-core
    for (int hi = 0; hi < nirrep_; hi++) {

        for (int i = 0; i < rstcpi_[hi] + frzcpi_[hi]; i++) {

            int ifull      = i + pitzer_offset_full[hi];

            for (int hj = 0; hj < nirrep_; hj++) {

                for (int j = 0; j < rstcpi_[hj] + frzcpi_[hj]; j++) {

                    int jfull      = j + pitzer_offset_full[hj];

                    tpdm d2;

                    d2.i   = ifull;
                    d2.j   = jfull;
                    d2.k   = ifull;
                    d2.l   = jfull;

                    d2.value = 1.0;
                    psio->write(PSIF_V2RDM_D2AB,"D2ab",(char*)&d2,sizeof(tpdm),addr_ab,&addr_ab);
                    countab++;

                    if ( ifull != jfull ) {

                        psio->write(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);
                        countaa++;

                        psio->write(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);
                        countbb++;

                        // ij;ji

                        d2.value = -1.0;
                        d2.k   = jfull;
                        d2.l   = ifull;

                        psio->write(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);
                        countaa++;

                        psio->write(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);
                        countbb++;

                    }
                }
            }
        }
    }

    // core active; core active
    for (int hi = 0; hi < nirrep_; hi++) {

        for (int i = 0; i < rstcpi_[hi] + frzcpi_[hi]; i++) {

            int ifull      = i + pitzer_offset_full[hi];

            // D2(ij; il) 
            for (int hj = 0; hj < nirrep_; hj++) {

                for (int j = 0; j < amopi_[hj]; j++) {

                    int jfull      = full_basis[j+pitzer_offset[hj]];

                    for (int l = 0; l < amopi_[hj]; l++) {

                        int lfull      = full_basis[l+pitzer_offset[hj]];

                        // aa and bb pieces
                        double valaa = x_p[d1aoff[hj]+j*amopi_[hj]+l];
                        double valbb = x_p[d1boff[hj]+j*amopi_[hj]+l];

                        tpdm d2;

                        // ij;il
                        d2.i   = ifull;
                        d2.j   = jfull;
                        d2.k   = ifull;
                        d2.l   = lfull;

                        d2.value = valaa;
                        psio->write(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);
                        countaa++;

                        d2.value = valbb;
                        psio->write(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);
                        countbb++;

                        // ij;li
                        d2.k   = lfull;
                        d2.l   = ifull;

                        d2.value = -valaa;
                        psio->write(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);
                        countaa++;

                        d2.value = -valbb;
                        psio->write(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);
                        countbb++;

                        // ji;li
                        d2.i   = jfull;
                        d2.j   = ifull;

                        d2.value = valaa;
                        psio->write(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);
                        countaa++;

                        d2.value = valbb;
                        psio->write(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);
                        countbb++;

                        // ji;il
                        d2.k   = ifull;
                        d2.l   = lfull;

                        d2.value = -valaa;
                        psio->write(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);
                        countaa++;

                        d2.value = -valbb;
                        psio->write(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);
                        countbb++;


                        // ab (ij;il) and ba (ji;li) pieces
                        double valab = x_p[d1boff[hj]+j*amopi_[hj]+l];
                        double valba = x_p[d1aoff[hj]+j*amopi_[hj]+l];

                        // ij;il
                        d2.i   = ifull;
                        d2.j   = jfull;

                        d2.value = valab;
                        psio->write(PSIF_V2RDM_D2AB,"D2ab",(char*)&d2,sizeof(tpdm),addr_ab,&addr_ab);
                        countab++;

                        // ji;li
                        d2.i   = jfull;
                        d2.j   = ifull;
                        d2.k   = lfull;
                        d2.l   = ifull;

                        d2.value = valba;
                        psio->write(PSIF_V2RDM_D2AB,"D2ab",(char*)&d2,sizeof(tpdm),addr_ab,&addr_ab);
                        countab++;

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

    psio->write_entry(PSIF_V2RDM_D2AA,"NUMBER ACTIVE ORBITALS",(char*)&amo_,sizeof(int));
    psio->write_entry(PSIF_V2RDM_D2BB,"NUMBER ACTIVE ORBITALS",(char*)&amo_,sizeof(int));
    psio->write_entry(PSIF_V2RDM_D2AB,"NUMBER ACTIVE ORBITALS",(char*)&amo_,sizeof(int));

    addr_aa = PSIO_ZERO;
    addr_bb = PSIO_ZERO;
    addr_ab = PSIO_ZERO;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            int ifull = full_basis[i + pitzer_offset[h]];
            psio->write(PSIF_V2RDM_D2AA,"ACTIVE ORBITALS",(char*)&ifull,sizeof(int),addr_aa,&addr_aa);
            psio->write(PSIF_V2RDM_D2BB,"ACTIVE ORBITALS",(char*)&ifull,sizeof(int),addr_bb,&addr_bb);
            psio->write(PSIF_V2RDM_D2AB,"ACTIVE ORBITALS",(char*)&ifull,sizeof(int),addr_ab,&addr_ab);
        }
    }
    
    // it might be nice for post CASSCF codes to know what orbitals are inactive:

    int inact = nfrzc_ + nrstc_;
    psio->write_entry(PSIF_V2RDM_D2AA,"NUMBER INACTIVE ORBITALS",(char*)&inact,sizeof(int));
    psio->write_entry(PSIF_V2RDM_D2BB,"NUMBER INACTIVE ORBITALS",(char*)&inact,sizeof(int));
    psio->write_entry(PSIF_V2RDM_D2AB,"NUMBER INACTIVE ORBITALS",(char*)&inact,sizeof(int));

    addr_aa = PSIO_ZERO;
    addr_bb = PSIO_ZERO;
    addr_ab = PSIO_ZERO;
    for (int hi = 0; hi < nirrep_; hi++) {
        for (int i = 0; i < rstcpi_[hi] + frzcpi_[hi]; i++) {
            int ifull = i + pitzer_offset_full[hi];
            psio->write(PSIF_V2RDM_D2AA,"INACTIVE ORBITALS",(char*)&ifull,sizeof(int),addr_aa,&addr_aa);
            psio->write(PSIF_V2RDM_D2BB,"INACTIVE ORBITALS",(char*)&ifull,sizeof(int),addr_bb,&addr_bb);
            psio->write(PSIF_V2RDM_D2AB,"INACTIVE ORBITALS",(char*)&ifull,sizeof(int),addr_ab,&addr_ab);
        }
    }

    // close files
    psio->close(PSIF_V2RDM_D2AA,1);
    psio->close(PSIF_V2RDM_D2BB,1);
    psio->close(PSIF_V2RDM_D2AB,1);

}
void v2RDMSolver::WriteTPDMSpinFree(){

    double * x_p = x->pointer();

    std::shared_ptr<PSIO> psio (new PSIO());

    psio->open(PSIF_V2RDM_D2_SPIN_FREE,PSIO_OPEN_NEW);

    psio_address addr = PSIO_ZERO;

    long int count = 0;

    // active-active part

    for (int h = 0; h < nirrep_; h++) {

        for (int ij_ab = 0; ij_ab < gems_ab[h]; ij_ab++) {

            int i     = bas_ab_sym[h][ij_ab][0];
            int j     = bas_ab_sym[h][ij_ab][1];

            int ij_aa = ibas_aa_sym[h][i][j];
            int ji_ab = ibas_ab_sym[h][j][i];

            int ifull = full_basis[i];
            int jfull = full_basis[j];

            for (int kl_ab = 0; kl_ab < gems_ab[h]; kl_ab++) {

                int k     = bas_ab_sym[h][kl_ab][0];
                int l     = bas_ab_sym[h][kl_ab][1];

                int kl_aa = ibas_aa_sym[h][k][l];
                int lk_ab = ibas_ab_sym[h][l][k];

                int kfull = full_basis[k];
                int lfull = full_basis[l];

                double val = 0.0;
                
                val += x_p[d2aboff[h] + ij_ab*gems_ab[h] + kl_ab];
                val += x_p[d2aboff[h] + ji_ab*gems_ab[h] + lk_ab];

                if ( i != j && k != l ) {
                    int sg = 1;
                    if ( i > j ) sg = -1;
                    if ( k > l ) sg = -1;
                    val += sg * x_p[d2aaoff[h] + ij_aa*gems_aa[h] + kl_aa];
                    val += sg * x_p[d2bboff[h] + ij_aa*gems_aa[h] + kl_aa];
                }

                tpdm d2;
                d2.i   = ifull;
                d2.j   = jfull;
                d2.k   = kfull;
                d2.l   = lfull;
                d2.value = val;
                psio->write(PSIF_V2RDM_D2_SPIN_FREE,"D2",(char*)&d2,sizeof(tpdm),addr,&addr);
                count++;

            }
        }
    }

    // core-core
    for (int hi = 0; hi < nirrep_; hi++) {

        for (int i = 0; i < rstcpi_[hi] + frzcpi_[hi]; i++) {

            int ifull      = i + pitzer_offset_full[hi];

            for (int hj = 0; hj < nirrep_; hj++) {

                for (int j = 0; j < rstcpi_[hj] + frzcpi_[hj]; j++) {

                    int jfull      = j + pitzer_offset_full[hj];

                    tpdm d2;

                    d2.i   = ifull;
                    d2.j   = jfull;
                    d2.k   = ifull;
                    d2.l   = jfull;

                    d2.value = 2.0;
                    if ( ifull != jfull ) {
                        d2.value += 2.0;
                    }

                    psio->write(PSIF_V2RDM_D2_SPIN_FREE,"D2",(char*)&d2,sizeof(tpdm),addr,&addr);
                    count++;

                    if ( ifull != jfull ) {

                        // ij;ji

                        d2.value = -2.0;
                        d2.k   = jfull;
                        d2.l   = ifull;

                        psio->write(PSIF_V2RDM_D2_SPIN_FREE,"D2",(char*)&d2,sizeof(tpdm),addr,&addr);
                        count++;

                    }
                }
            }
        }
    }

    // core active; core active
    for (int hi = 0; hi < nirrep_; hi++) {

        for (int i = 0; i < rstcpi_[hi] + frzcpi_[hi]; i++) {
            int ifull      = i + pitzer_offset_full[hi];

            // D2(ij; il) 
            for (int hj = 0; hj < nirrep_; hj++) {

                for (int j = 0; j < amopi_[hj]; j++) {

                    int jfull      = full_basis[j+pitzer_offset[hj]];

                    for (int l = 0; l < amopi_[hj]; l++) {

                        int lfull      = full_basis[l+pitzer_offset[hj]];

                        // aa and bb pieces
                        double valaa = x_p[d1aoff[hj]+j*amopi_[hj]+l];
                        double valbb = x_p[d1boff[hj]+j*amopi_[hj]+l];

                        // ab (ij;il) and ba (ji;li) pieces
                        double valab = x_p[d1boff[hj]+j*amopi_[hj]+l];
                        double valba = x_p[d1aoff[hj]+j*amopi_[hj]+l];

                        tpdm d2;

                        // ij;il
                        d2.i   = ifull;
                        d2.j   = jfull;
                        d2.k   = ifull;
                        d2.l   = lfull;

                        d2.value = valaa + valbb + valab + valba;
                        psio->write(PSIF_V2RDM_D2_SPIN_FREE,"D2",(char*)&d2,sizeof(tpdm),addr,&addr);
                        count++;

                        // ij;li
                        d2.k   = lfull;
                        d2.l   = ifull;

                        d2.value = - valaa - valbb;
                        psio->write(PSIF_V2RDM_D2_SPIN_FREE,"D2",(char*)&d2,sizeof(tpdm),addr,&addr);
                        count++;

                        // ji;li
                        d2.i   = jfull;
                        d2.j   = ifull;

                        d2.value = valaa + valbb + valab + valba;
                        psio->write(PSIF_V2RDM_D2_SPIN_FREE,"D2",(char*)&d2,sizeof(tpdm),addr,&addr);
                        count++;

                        // ji;il
                        d2.k   = ifull;
                        d2.l   = lfull;

                        d2.value = - valaa - valbb;
                        psio->write(PSIF_V2RDM_D2_SPIN_FREE,"D2",(char*)&d2,sizeof(tpdm),addr,&addr);
                        count++;

                    }
                }
            }
        }
    }

    // write the number of entries in each file
    psio->write_entry(PSIF_V2RDM_D2_SPIN_FREE,"length",(char*)&count,sizeof(long int));

    // it might be nice for post CASSCF codes to know what orbitals are active:

    psio->write_entry(PSIF_V2RDM_D2_SPIN_FREE,"NUMBER ACTIVE ORBITALS",(char*)&amo_,sizeof(int));

    addr = PSIO_ZERO;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            int ifull = full_basis[i + pitzer_offset[h]];
            psio->write(PSIF_V2RDM_D2_SPIN_FREE,"ACTIVE ORBITALS",(char*)&ifull,sizeof(int),addr,&addr);
        }
    }
    
    // it might be nice for post CASSCF codes to know what orbitals are inactive:

    int inact = nfrzc_ + nrstc_;
    psio->write_entry(PSIF_V2RDM_D2_SPIN_FREE,"NUMBER INACTIVE ORBITALS",(char*)&inact,sizeof(int));

    addr = PSIO_ZERO;
    for (int hi = 0; hi < nirrep_; hi++) {
        for (int i = 0; i < rstcpi_[hi] + frzcpi_[hi]; i++) {
            int ifull = i + pitzer_offset_full[hi];
            psio->write(PSIF_V2RDM_D2_SPIN_FREE,"INACTIVE ORBITALS",(char*)&ifull,sizeof(int),addr,&addr);
        }
    }

    // close file
    psio->close(PSIF_V2RDM_D2_SPIN_FREE,1);

}

void v2RDMSolver::WriteActiveTPDM(){

    double * x_p = x->pointer();

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
            int ifull = full_basis[i];
            int jfull = full_basis[j];

            for (int kl = 0; kl < gems_ab[h]; kl++) {

                int k     = bas_ab_sym[h][kl][0];
                int l     = bas_ab_sym[h][kl][1];
                int kfull = full_basis[k];
                int lfull = full_basis[l];

                double valab = x_p[d2aboff[h] + ij*gems_ab[h] + kl];

                tpdm d2;
                d2.i   = ifull;
                d2.j   = jfull;
                d2.k   = kfull;
                d2.l   = lfull;
                d2.value = valab;
                psio->write(PSIF_V2RDM_D2AB,"D2ab",(char*)&d2,sizeof(tpdm),addr_ab,&addr_ab);
                countab++;

                if ( i != j && k != l ) {

                    int ija = ibas_aa_sym[h][i][j];
                    int kla = ibas_aa_sym[h][k][l];

                    int sij = i < j ? 1 : -1;
                    int skl = k < l ? 1 : -1;

                    double valaa = sij * skl * x_p[d2aaoff[h] + ija*gems_aa[h] + kla];
                    double valbb = sij * skl * x_p[d2bboff[h] + ija*gems_aa[h] + kla];

                    d2.value = valaa;
                    psio->write(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);
                    countaa++;

                    d2.value = valbb;
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

    // close files
    psio->close(PSIF_V2RDM_D2AA,1);
    psio->close(PSIF_V2RDM_D2BB,1);
    psio->close(PSIF_V2RDM_D2AB,1);

}

void v2RDMSolver::ReadTPDM(){

    std::shared_ptr<PSIO> psio (new PSIO());

    if ( !psio->exists(PSIF_V2RDM_D2AB) ) return;
    if ( !psio->exists(PSIF_V2RDM_D2AA) ) return;
    if ( !psio->exists(PSIF_V2RDM_D2BB) ) return;

    //Ca_->print();

    double * D2aa = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D2bb = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D2ab = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));

    memset((void*)D2aa,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2bb,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2ab,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    psio_address addr_aa = PSIO_ZERO;
    psio_address addr_bb = PSIO_ZERO;
    psio_address addr_ab = PSIO_ZERO;

    // ab
    psio->open(PSIF_V2RDM_D2AB,PSIO_OPEN_OLD);

    long int nab;
    psio->read_entry(PSIF_V2RDM_D2AB,"length",(char*)&nab,sizeof(long int));

    for (int n = 0; n < nab; n++) {
        tpdm d2;
        psio->read(PSIF_V2RDM_D2AB,"D2ab",(char*)&d2,sizeof(tpdm),addr_ab,&addr_ab);
        int i = d2.i;
        int j = d2.j;
        int k = d2.k;
        int l = d2.l;
        long int id = i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l;
        D2ab[id] = d2.value;
    }
    psio->close(PSIF_V2RDM_D2AB,1);

    // aa
    psio->open(PSIF_V2RDM_D2AA,PSIO_OPEN_OLD);

    long int naa;
    psio->read_entry(PSIF_V2RDM_D2AA,"length",(char*)&naa,sizeof(long int));

    for (int n = 0; n < naa; n++) {
        tpdm d2;
        psio->read(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);
        int i = d2.i;
        int j = d2.j;
        int k = d2.k;
        int l = d2.l;
        long int id = i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l;
        D2aa[id] = d2.value;
    }
    psio->close(PSIF_V2RDM_D2AA,1);

    // bb
    psio->open(PSIF_V2RDM_D2BB,PSIO_OPEN_OLD);

    long int nbb;
    psio->read_entry(PSIF_V2RDM_D2BB,"length",(char*)&nbb,sizeof(long int));

    for (int n = 0; n < nbb; n++) {
        tpdm d2;
        psio->read(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);
        int i = d2.i;
        int j = d2.j;
        int k = d2.k;
        int l = d2.l;
        long int id = i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l;
        D2bb[id] = d2.value;
    }
    psio->close(PSIF_V2RDM_D2BB,1);

    // check traces:
    double traa = 0.0;
    double trbb = 0.0;
    double trab = 0.0;
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            traa += D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j];
            trbb += D2bb[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j];
            trab += D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j];
        }
    }
    printf("  tr(d2aa) = %20.12lf\n",traa); fflush(stdout);
    printf("  tr(d2bb) = %20.12lf\n",trbb); fflush(stdout);
    printf("  tr(d2ab) = %20.12lf\n",trab); fflush(stdout);

    double * Da = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * Db = (double*)malloc(nmo_*nmo_*sizeof(double));

    memset((void*)Da,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)Db,'\0',nmo_*nmo_*sizeof(double));

    double tra = 0.0;
    double trb = 0.0;

    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {

            double duma = 0.0;
            double dumb = 0.0;
            for (int k = 0; k < nmo_; k++) {
                duma += D2ab[i*nmo_*nmo_*nmo_+k*nmo_*nmo_+j*nmo_+k];
                duma += D2aa[i*nmo_*nmo_*nmo_+k*nmo_*nmo_+j*nmo_+k];

                dumb += D2ab[k*nmo_*nmo_*nmo_+i*nmo_*nmo_+k*nmo_+j];
                dumb += D2bb[i*nmo_*nmo_*nmo_+k*nmo_*nmo_+j*nmo_+k];
            }
            Da[i*nmo_+j] = 1.0/(nalpha_+nbeta_-1.0) * duma;
            Db[i*nmo_+j] = 1.0/(nalpha_+nbeta_-1.0) * dumb;

            if ( i == j ) {
                tra += Da[i*nmo_+j];
                trb += Db[i*nmo_+j];
            }

        }
    }

    printf("  tr(da) = %20.12lf\n",tra); fflush(stdout);
    printf("  tr(db) = %20.12lf\n",trb); fflush(stdout);

    // check energy:

    double en2 = 0.0;
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            for (int k = 0; k < nmo_; k++) {
                for (int l = 0; l < nmo_; l++) {

                    double eri = C_DDOT(nQ_,Qmo_ + INDEX(i,k),(nmo_-nfrzv_)*(nmo_-nfrzv_+1)/2,Qmo_+INDEX(j,l),(nmo_-nfrzv_)*(nmo_-nfrzv_+1)/2);
                    
                    en2 +=       eri * D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];
                    en2 += 0.5 * eri * D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];
                    en2 += 0.5 * eri * D2bb[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];

                }
            }
        }
    }

    std::shared_ptr<MintsHelper> mints(new MintsHelper(reference_wavefunction_));
    SharedMatrix K1 (new Matrix(mints->so_potential()));
    K1->add(mints->so_kinetic());
    K1->transform(Ca_);

    double en1 = 0.0;

    long int offset = 0;
    long int offset2 = 0;
    for (int h = 0; h < nirrep_; h++) {

        for (int i = 0; i < nmopi_[h]; i++) {

            int ifull = i + offset;

            for (int j = 0; j < nmopi_[h]; j++) {

                int jfull = j + offset;


                en1 += oei_full_sym_[offset2 + INDEX(i,j)] * Da[ifull*nmo_+jfull];
                en1 += oei_full_sym_[offset2 + INDEX(i,j)] * Db[ifull*nmo_+jfull];
                //en1 += K1->pointer(h)[i][j] * Da[ifull*nmo_+jfull];
                //en1 += K1->pointer(h)[i][j] * Db[ifull*nmo_+jfull];
            }
        }

        offset  += nmopi_[h] - frzvpi_[h];
        offset2 += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;

    }

    printf("%20.12lf %20.12lf %20.12lf %20.12lf\n",en1,en2,enuc_,en1+en2+enuc_);

    free(D2aa);
    free(D2bb);
    free(D2ab);
    free(Da);
    free(Db);
}


} //end namespaces


