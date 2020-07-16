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

#ifndef _python_api_h_
#define _python_api_h_

#include "v2rdm_solver.h"

#include <psi4/libpsi4util/process.h>
#include <psi4/libmints/wavefunction.h>

using namespace psi;


namespace hilbert {

std::vector<opdm> v2RDMSolver::get_opdm_sparse(std::string type) {

    if ( type != "A" && type != "B" && type != "SUM" ) {
        throw PsiException("invalid opdm type.",__FILE__,__LINE__);
    }

    std::vector<opdm> my_opdm;

    double * x_p = x->pointer();

    // active-active part
    for (int h = 0; h < nirrep_; h++) {

        for (int i = 0; i < amopi_[h]; i++) {

            int ifull = full_basis[i + pitzer_offset[h]];

            for (int j = 0; j < amopi_[h]; j++) {

                int jfull = full_basis[j + pitzer_offset[h]];

                double vala = x_p[d1aoff[h] + i*amopi_[h] + j];
                double valb = x_p[d1boff[h] + i*amopi_[h] + j];

                opdm d1;

                d1.i   = ifull;
                d1.j   = jfull;

                if ( type == "A" ) {
                    d1.value = vala;
                }else if ( type == "B") {
                    d1.value = valb;
                }else if ( type == "SUM") {
                    d1.value = vala + valb;
                }

                my_opdm.push_back(d1);

            }
        }
    }

    // core-core
    for (int hi = 0; hi < nirrep_; hi++) {

        for (int i = 0; i < rstcpi_[hi] + frzcpi_[hi]; i++) {

            int ifull      = i + pitzer_offset_full[hi];

            opdm d1;

            d1.i   = ifull;
            d1.j   = ifull;

            if ( type == "A" ) {
                d1.value = 1.0;
            }else if ( type == "B") {
                d1.value = 1.0;
            }else if ( type == "SUM") {
                d1.value = 2.0;
            }

            my_opdm.push_back(d1);

        }
    }

    return my_opdm;

}

std::vector<tpdm> v2RDMSolver::get_tpdm_sparse(std::string type) {

    if ( type != "AA" && type != "BB"  && type != "AB" &&
         type != "BA" && type != "SUM" ) {
        throw PsiException("invalid tpdm type.",__FILE__,__LINE__);
    }

    std::vector<tpdm> my_tpdm;

    double * x_p = x->pointer();

    // active-active part
    for (int h = 0; h < nirrep_; h++) {

        for (int ij = 0; ij < gems_ab[h]; ij++) {

            int i     = bas_ab_sym[h][ij][0];
            int j     = bas_ab_sym[h][ij][1];
            int ji    = ibas_ab_sym[h][j][i];
            int ifull = full_basis[i];
            int jfull = full_basis[j];

            for (int kl = 0; kl < gems_ab[h]; kl++) {

                int k     = bas_ab_sym[h][kl][0];
                int l     = bas_ab_sym[h][kl][1];
                int lk    = ibas_ab_sym[h][l][k];
                int kfull = full_basis[k];
                int lfull = full_basis[l];

                double valab = x_p[d2aboff[h] + ij*gems_ab[h] + kl];
                double valba = x_p[d2aboff[h] + ji*gems_ab[h] + lk];

                tpdm d2;
                d2.i   = ifull;
                d2.j   = jfull;
                d2.k   = kfull;
                d2.l   = lfull;

                double valaa = 0.0;
                double valbb = 0.0;
                if ( i != j && k != l ) {

                    int ija = ibas_aa_sym[h][i][j];
                    int kla = ibas_aa_sym[h][k][l];

                    int sij = i < j ? 1 : -1;
                    int skl = k < l ? 1 : -1;

                    valaa = sij * skl * x_p[d2aaoff[h] + ija*gems_aa[h] + kla];
                    valbb = sij * skl * x_p[d2bboff[h] + ija*gems_aa[h] + kla];

                }

                if ( type == "AA" ) {
                    d2.value = valaa;
                }else if ( type == "BB" ) {
                    d2.value = valbb;
                }else if ( type == "AB" ) {
                    d2.value = valab;
                }else if ( type == "BA" ) {
                    d2.value = valba;
                }else if ( type == "SUM" ) {
                    d2.value = valaa + valbb + valab + valba;
                }

                my_tpdm.push_back(d2);

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

                    // ij;ij

                    double valab = 1.0;
                    double valba = 1.0;
                    double valaa = 0.0;
                    double valbb = 0.0;

                    if ( ifull != jfull ) {

                        valaa = 1.0;
                        valbb = 1.0;

                    }

                    if ( type == "AA" ) {
                        d2.value = valaa;
                    }else if ( type == "BB" ) {
                        d2.value = valbb;
                    }else if ( type == "AB" ) {
                        d2.value = valab;
                    }else if ( type == "BA" ) {
                        d2.value = valba;
                    }else if ( type == "SUM" ) {
                        d2.value = valaa + valbb + valab + valba;
                    }

                    my_tpdm.push_back(d2);

                    // ij;ji
                    if ( ifull != jfull ) {

                        d2.k   = jfull;
                        d2.l   = ifull;
                        valaa  = -1.0;
                        valbb  = -1.0;
                        if ( type == "AA" ) {
                            d2.value = valaa;
                        }else if ( type == "BB" ) {
                            d2.value = valbb;
                        }else if ( type == "AB" ) { // skip
                            d2.value = 0.0;
                        }else if ( type == "BA" ) { // skip
                            d2.value = 0.0;
                        }else if ( type == "SUM" ) {
                            d2.value = valaa + valbb;
                        }

                        my_tpdm.push_back(d2);
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

                        tpdm d2;

                        // ij;il (aa, bb, ab, ba)
                        d2.i   = ifull;
                        d2.j   = jfull;
                        d2.k   = ifull;
                        d2.l   = lfull;

                        // aa, bb, ab, ba pieces
                        double vala = x_p[d1aoff[hj]+j*amopi_[hj]+l];
                        double valb = x_p[d1boff[hj]+j*amopi_[hj]+l];

                        if ( type == "AA" ) {
                            d2.value = vala;
                        }else if ( type == "BB" ) {
                            d2.value = valb;
                        }else if ( type == "AB" ) {
                            d2.value = valb;
                        }else if ( type == "BA" ) {
                            d2.value = vala;
                        }else if ( type == "SUM" ) {
                            d2.value = 2.0 * vala + 2.0 * valb;
                        }

                        my_tpdm.push_back(d2);

                        // ji;li (aa, bb, ab, ba)
                        d2.i   = jfull;
                        d2.j   = ifull;
                        d2.k   = lfull;
                        d2.l   = ifull;

                        if ( type == "AA" ) {
                            d2.value = vala;
                        }else if ( type == "BB" ) {
                            d2.value = valb;
                        }else if ( type == "AB" ) {
                            d2.value = vala;
                        }else if ( type == "BA" ) {
                            d2.value = valb;
                        }else if ( type == "SUM" ) {
                            d2.value = 2.0 * vala + 2.0 * valb;
                        }

                        my_tpdm.push_back(d2);

                        // ij;li
                        d2.i   = ifull;
                        d2.j   = jfull;
                        d2.k   = lfull;
                        d2.l   = ifull;

                        if ( type == "AA" ) {
                            d2.value = -vala;
                        }else if ( type == "BB" ) {
                            d2.value = -valb;
                        }else if ( type == "AB" ) { // skip
                            d2.value = 0.0;
                        }else if ( type == "BA" ) { // skip
                            d2.value = 0.0;
                        }else if ( type == "SUM" ) {
                            d2.value = -vala - valb;
                        }

                        my_tpdm.push_back(d2);

                        // ji;il
                        d2.i   = jfull;
                        d2.j   = ifull;
                        d2.k   = ifull;
                        d2.l   = lfull;

                        if ( type == "AA" ) {
                            d2.value = -vala;
                        }else if ( type == "BB" ) {
                            d2.value = -valb;
                        }else if ( type == "AB" ) { // skip
                            d2.value = 0.0;
                        }else if ( type == "BA" ) { // skip
                            d2.value = 0.0;
                        }else if ( type == "SUM" ) {
                            d2.value = -vala - valb;
                        }

                        my_tpdm.push_back(d2);

                    }
                }
            }
        }
    }

    return my_tpdm;

}


std::shared_ptr<Matrix> v2RDMSolver::get_opdm() {

    Dimension amopi = Dimension(nirrep_);
    for (int h = 0; h < nirrep_; h++) {
        amopi[h] = amopi_[h];
    }
    std::shared_ptr<Matrix> opdm(new Matrix("OPDM",amopi,amopi));

    double * x_p = x->pointer();
    for (int h = 0; h < nirrep_; h++) {
        double ** opdm_p = opdm->pointer(h);
        for (int i = 0; i < amopi_[h]; i++) {
            for (int j = 0; j < amopi_[h]; j++) {
                opdm_p[i][j]  = x_p[d1aoff[h]+i*amopi_[h]+j];
                opdm_p[i][j] += x_p[d1aoff[h]+i*amopi_[h]+j];
            }
        }
    }

    return opdm;
}

std::shared_ptr<Matrix> v2RDMSolver::get_tpdm() {

    if ( nirrep_ > 1 ) {
        throw PsiException("not sure how CIWavefunction stores TPDM with symmetry",__FILE__,__LINE__);
    }
    //Dimension amopi = Dimension(nirrep_);
    //for (int h = 0; h < nirrep_; h++) {
    //    amopi[h] = amopi_[h];
    //}
    std::shared_ptr<Matrix> tpdm(new Matrix("TPDM",amo_*amo_,amo_*amo_));

    double ** tpdm_p = tpdm->pointer();
    double * x_p = x->pointer();
    for (int i = 0; i < amo_; i++) {
        for (int j = 0; j < amo_; j++) {
            int ij_ab = ibas_ab_sym[0][i][j];
            int ji_ab = ibas_ab_sym[0][j][i];
            int ij_aa = ibas_aa_sym[0][i][j];
            for (int k = 0; k < amo_; k++) {
                for (int l = 0; l < amo_; l++) {
                    int kl_ab = ibas_ab_sym[0][k][l];
                    int lk_ab = ibas_ab_sym[0][l][k];
                    int il_ab = ibas_ab_sym[0][i][l];
                    int kj_ab = ibas_ab_sym[0][k][j];
                    int li_ab = ibas_ab_sym[0][l][i];
                    int jk_ab = ibas_ab_sym[0][j][k];

                    int kl_aa = ibas_aa_sym[0][k][l];
                    int il_aa = ibas_aa_sym[0][i][l];
                    int kj_aa = ibas_aa_sym[0][k][j];

                    double dum_ikjl = x_p[d2aboff[0] + ij_ab * gems_ab[0] + kl_ab];
                    dum_ikjl       += x_p[d2aboff[0] + ji_ab * gems_ab[0] + lk_ab];

                    double dum_iklj = x_p[d2aboff[0] + il_ab * gems_ab[0] + kj_ab];
                    dum_iklj       += x_p[d2aboff[0] + li_ab * gems_ab[0] + jk_ab];

                    if ( i != j && k != l ) {
                        int sg = 1;
                        if ( i > j ) sg = -sg;
                        if ( k > l ) sg = -sg;
                        dum_ikjl += sg * x_p[d2aaoff[0] + ij_aa * gems_aa[0] + kl_aa];
                        dum_ikjl += sg * x_p[d2bboff[0] + ij_aa * gems_aa[0] + kl_aa];
                    }
                    if ( i != l && k != j ) {
                        int sg = 1;
                        if ( i > l ) sg = -sg;
                        if ( k > j ) sg = -sg;
                        dum_ikjl += sg * x_p[d2aaoff[0] + il_aa * gems_aa[0] + kj_aa];
                        dum_ikjl += sg * x_p[d2bboff[0] + il_aa * gems_aa[0] + kj_aa];
                    }

                    //tpdm_p[i*amo_+j][k*amo_+l] =  dum;
                    tpdm_p[i*amo_+k][j*amo_+l] =  0.5 * (dum_ikjl + dum_iklj);
                    //tpdm_p[i*amo_+j][k*amo_+l] =  0.5 * (dum_ijkl + dum_ijlk);
                }
            }
        }
    }
    return tpdm;
}

void v2RDMSolver::orbital_locations(const std::string& orbitals, int* start, int* end) {
    if (orbitals == "FROZEN_DOCC") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = 0;
            end[h] = frzcpi_[h];
        }
    } else if (orbitals == "RESTRICTED_DOCC") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = frzcpi_[h];
            end[h] = frzcpi_[h] + rstcpi_[h];
        }
    } else if (orbitals == "DOCC") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = 0;
            end[h] = frzcpi_[h] + rstcpi_[h];
        }
    } else if (orbitals == "ACTIVE") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = frzcpi_[h] + rstcpi_[h];
            end[h] = frzcpi_[h] + rstcpi_[h] + amopi_[h];
        }
    } else if (orbitals == "RESTRICTED_UOCC") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = frzcpi_[h] + rstcpi_[h] + amopi_[h];
            end[h] = frzcpi_[h] + rstcpi_[h] + amopi_[h] + rstvpi_[h];
        }
    } else if (orbitals == "FROZEN_UOCC") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = frzcpi_[h] + rstcpi_[h] + amopi_[h] + rstvpi_[h];
            end[h] = frzcpi_[h] + rstcpi_[h] + amopi_[h] + rstvpi_[h] + frzvpi_[h];
        }
    } else if (orbitals == "VIRTUAL") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = frzcpi_[h] + rstcpi_[h] + amopi_[h];
            end[h] = frzcpi_[h] + rstcpi_[h] + amopi_[h] + rstvpi_[h] + frzvpi_[h];
        }

    } else if (orbitals == "ALL") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = 0;
            end[h] = nmopi_[h];
        }
    } else if (orbitals == "ROT") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = frzcpi_[h];
            end[h] = nmopi_[h] - frzvpi_[h];
        }
    } else {
        throw PSIEXCEPTION(
            "v2RDMSolver: Orbital subset is not defined, should be FROZEN_DOCC, "
            "RESTRICTED_DOCC, DOCC, ACTIVE, RESTRICTED_UOCC, FROZEN_UOCC, "
            "VIRTUAL, ROT, or ALL");
    }
}

std::shared_ptr<Matrix> v2RDMSolver::get_orbitals(const std::string& orbital_name) {
    /// Figure out orbital positions
    auto* start = new int[nirrep_];
    auto* end = new int[nirrep_];

    orbital_locations(orbital_name, start, end);

    auto* spread = new int[nirrep_];
    for (int h = 0; h < nirrep_; h++) {
        spread[h] = end[h] - start[h];
    }

    /// Fill desired orbitals
    auto retC = std::make_shared<Matrix>("C " + orbital_name, nirrep_, nsopi_, spread);
    for (int h = 0; h < nirrep_; h++) {
        for (int i = start[h], pos = 0; i < end[h]; i++, pos++) {
            C_DCOPY(nsopi_[h], &Ca_->pointer(h)[0][i], nmopi_[h], &retC->pointer(h)[0][pos], spread[h]);
        }
    }

    /// Cleanup
    delete[] start;
    delete[] end;
    delete[] spread;

    return retC;
}

void v2RDMSolver::set_orbitals(const std::string& orbital_name, SharedMatrix orbitals) {
    /// Figure out orbital positions
    auto* start = new int[nirrep_];
    auto* end = new int[nirrep_];

    orbital_locations(orbital_name, start, end);

    auto* spread = new int[nirrep_];
    for (int h = 0; h < nirrep_; h++) {
        spread[h] = end[h] - start[h];
    }

    /// Fill desired orbitals
    for (int h = 0; h < nirrep_; h++) {
        for (int i = start[h], pos = 0; i < end[h]; i++, pos++) {
            C_DCOPY(nsopi_[h], &orbitals->pointer(h)[0][pos], spread[h], &Ca_->pointer(h)[0][i], nmopi_[h]);
        }
    }

    /// Cleanup
    delete[] start;
    delete[] end;
    delete[] spread;
}


} // End namespaces

#endif
