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
#include <psi4/libiwl/iwl.hpp>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libdpd/dpd.h>
#include <psi4/libqt/qt.h>
#include <psi4/libdpd/dpd.h>

#define EXTERN

#include "v2rdm_solver.h"

using namespace psi;

namespace hilbert{

// function to write TPDM to disk in MO basis for subsequent use in deriv()
void v2RDMSolver::WriteTPDM_IWL(){

    double * x_p = x->pointer();

    std::shared_ptr<PSIO> psio (new PSIO());

    IWL d2aa(psio.get(), PSIF_MO_AA_TPDM, 1.0e-14, 0, 0);
    IWL d2ab(psio.get(), PSIF_MO_AB_TPDM, 1.0e-14, 0, 0);
    IWL d2bb(psio.get(), PSIF_MO_BB_TPDM, 1.0e-14, 0, 0);

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

                d2ab.write_value(ifull, kfull, jfull, lfull, valab, 0, "NULL", 0);

                // NOTE (AED): I can't figure out why I need a factor of 1/4 here.  1/2 makes sense,
                // but I need the 1/4 to get the two-electron part of the gradient correct for 
                // an exact case (FH/STO-3G, a two-hole system)
                if ( i != j && k != l ) {

                    int ija = ibas_aa_sym[h][i][j];
                    int kla = ibas_aa_sym[h][k][l];

                    int sij = i < j ? 1 : -1;
                    int skl = k < l ? 1 : -1;

                    double valaa = 0.25 * sij * skl * x_p[d2aaoff[h] + ija*gems_aa[h] + kla];
                    double valbb = 0.25 * sij * skl * x_p[d2bboff[h] + ija*gems_aa[h] + kla];

                    d2aa.write_value(ifull, kfull, jfull, lfull, valaa, 0, "NULL", 0);
                    d2bb.write_value(ifull, kfull, jfull, lfull, valbb, 0, "NULL", 0);
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

                    d2ab.write_value(ifull, ifull, jfull, jfull, 1.0, 0, "NULL", 0);

                    // NOTE (AED): I can't figure out why I need a factor of 1/4 here.  1/2 makes sense,
                    // but I need the 1/4 to get the two-electron part of the gradient correct for 
                    // Hartree-Fock.
                    if ( ifull != jfull ) {

                        d2aa.write_value(ifull, ifull, jfull, jfull, 0.25, 0, "NULL", 0);
                        d2bb.write_value(ifull, ifull, jfull, jfull, 0.25, 0, "NULL", 0);

                        // ij;ji

                        d2aa.write_value(ifull, jfull, jfull, ifull, -0.25, 0, "NULL", 0);
                        d2bb.write_value(ifull, jfull, jfull, ifull, -0.25, 0, "NULL", 0);

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

                        // NOTE (AED): I haven't verified this yet, but I'm assuming I will need the 
                        // same factor of 1/4 for the aa/bb blocks here as I needed above.

                        double valaa = 0.25 * x_p[d1aoff[hj]+j*amopi_[hj]+l];
                        double valbb = 0.25 * x_p[d1boff[hj]+j*amopi_[hj]+l];

                        //// ij;il

                        d2aa.write_value(ifull, ifull, jfull, lfull, valaa, 0, "NULL", 0);
                        d2bb.write_value(ifull, ifull, jfull, lfull, valbb, 0, "NULL", 0);

                        // ij;li
                        d2aa.write_value(ifull, lfull, jfull, ifull, -valaa, 0, "NULL", 0);
                        d2bb.write_value(ifull, lfull, jfull, ifull, -valbb, 0, "NULL", 0);

                        // ji;li
                        d2aa.write_value(jfull, lfull, ifull, ifull, valaa, 0, "NULL", 0);
                        d2bb.write_value(jfull, lfull, ifull, ifull, valbb, 0, "NULL", 0);

                        // ji;il
                        d2aa.write_value(jfull, ifull, ifull, lfull, -valaa, 0, "NULL", 0);
                        d2bb.write_value(jfull, ifull, ifull, lfull, -valbb, 0, "NULL", 0);


                        // ab (ij;il) and ba (ji;li) pieces
                        double valab = x_p[d1boff[hj]+j*amopi_[hj]+l];
                        double valba = x_p[d1aoff[hj]+j*amopi_[hj]+l];

                        // ij;il
                        d2ab.write_value(ifull, ifull, jfull, lfull, valab, 0, "NULL", 0);

                        // ji;li
                        d2ab.write_value(jfull, lfull, ifull, ifull, valba, 0, "NULL", 0);

                    }
                }
            }
        }
    }

    d2aa.flush(1);
    d2bb.flush(1);
    d2ab.flush(1);

    d2aa.set_keep_flag(1);
    d2bb.set_keep_flag(1);
    d2ab.set_keep_flag(1);

    d2aa.close();
    d2bb.close();
    d2ab.close();

}

} //end namespaces
