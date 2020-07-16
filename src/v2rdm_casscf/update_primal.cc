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
#include <psi4/liboptions/liboptions.h>
#include <psi4/libqt/qt.h>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libtrans/mospace.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/matrix.h>

#include "v2rdm_solver.h"

#include <misc/omp.h>

using namespace psi;

namespace hilbert{

// transform primal solute to new basis after semicanonical orbital transformation
// note, we're only transforming the D1/D2/D3
void v2RDMSolver::UpdatePrimal() {

    // D1a, D1b
    double * z_p  = z->pointer();
    double * x_p  = x->pointer();
    double * tmp_p = ATy->pointer();
    for (int h = 0; h < nirrep_; h++) {
        double ** t_p = newMO_->pointer(h);
        for (int i = 0; i < amopi_[h]; i++) {
            for (int j = 0; j < amopi_[h]; j++) {
                double duma = 0.0;
                double dumb = 0.0;
                for (int k = 0; k < amopi_[h]; k++) {
                    duma += x_p[d1aoff[h] + i*amopi_[h] + k] * t_p[j+frzcpi_[h]+rstcpi_[h]][k+frzcpi_[h]+rstcpi_[h]];
                    dumb += x_p[d1boff[h] + i*amopi_[h] + k] * t_p[j+frzcpi_[h]+rstcpi_[h]][k+frzcpi_[h]+rstcpi_[h]];
                }
                tmp_p[d1aoff[h] + i*amopi_[h] + j] = duma;
                tmp_p[d1boff[h] + i*amopi_[h] + j] = dumb;
            }
        }
        for (int i = 0; i < amopi_[h]; i++) {
            for (int j = 0; j < amopi_[h]; j++) {
                double duma = 0.0;
                double dumb = 0.0;
                for (int k = 0; k < amopi_[h]; k++) {
                    duma += tmp_p[d1aoff[h] + k*amopi_[h] + i] * t_p[j+frzcpi_[h]+rstcpi_[h]][k+frzcpi_[h]+rstcpi_[h]];
                    dumb += tmp_p[d1boff[h] + k*amopi_[h] + i] * t_p[j+frzcpi_[h]+rstcpi_[h]][k+frzcpi_[h]+rstcpi_[h]];
                }
                x_p[d1aoff[h] + j*amopi_[h] + i] = duma;
                x_p[d1boff[h] + j*amopi_[h] + i] = dumb;
            }
        }
    }

    // transform D2ab
    TransformFourIndex(x_p+d2aboff[0],tmp_p+d2aboff[0],newMO_);

    // unpack D2aa block and copy into z
    memset((void*)z_p,'\0',dimx_*sizeof(double));
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            int ijb = ibas_ab_sym[h][i][j];
            int jib = ibas_ab_sym[h][j][i];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                int klb = ibas_ab_sym[h][k][l];
                int lkb = ibas_ab_sym[h][l][k];
                double dum = x_p[d2aaoff[h] + ij*gems_aa[h] + kl];
                z_p[d2aboff[h] + ijb*gems_ab[h] + klb] =  dum;
                z_p[d2aboff[h] + ijb*gems_ab[h] + lkb] = -dum;
                z_p[d2aboff[h] + jib*gems_ab[h] + klb] = -dum;
                z_p[d2aboff[h] + jib*gems_ab[h] + lkb] =  dum;
            }
        }
    }

    // transform D2aa
    TransformFourIndex(z_p+d2aboff[0],tmp_p+d2aboff[0],newMO_);

    // repack D2aa
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            int ijb = ibas_ab_sym[h][i][j];
            int jib = ibas_ab_sym[h][j][i];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                int klb = ibas_ab_sym[h][k][l];
                int lkb = ibas_ab_sym[h][l][k];
                x_p[d2aaoff[h] + ij*gems_aa[h] + kl] = z_p[d2aboff[h] + ijb*gems_ab[h] + klb];
            }
        }
    }

    // unpack D2bb block and copy into z
    memset((void*)z_p,'\0',dimx_*sizeof(double));
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            int ijb = ibas_ab_sym[h][i][j];
            int jib = ibas_ab_sym[h][j][i];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                int klb = ibas_ab_sym[h][k][l];
                int lkb = ibas_ab_sym[h][l][k];
                double dum = x_p[d2bboff[h] + ij*gems_aa[h] + kl];
                z_p[d2aboff[h] + ijb*gems_ab[h] + klb] =  dum;
                z_p[d2aboff[h] + ijb*gems_ab[h] + lkb] = -dum;
                z_p[d2aboff[h] + jib*gems_ab[h] + klb] = -dum;
                z_p[d2aboff[h] + jib*gems_ab[h] + lkb] =  dum;
            }
        }
    }

    // transform D2bb
    TransformFourIndex(z_p+d2aboff[0],tmp_p+d2aboff[0],newMO_);

    // repack D2bb
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            int ijb = ibas_ab_sym[h][i][j];
            int jib = ibas_ab_sym[h][j][i];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                int klb = ibas_ab_sym[h][k][l];
                int lkb = ibas_ab_sym[h][l][k];
                x_p[d2bboff[h] + ij*gems_aa[h] + kl] = z_p[d2aboff[h] + ijb*gems_ab[h] + klb];
            }
        }
    }

}

void v2RDMSolver::TransformFourIndex(double * inout, double * tmp, SharedMatrix trans) {

    // transform fourth index 
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            for (int kl = 0; kl < gems_ab[h]; kl++) {
                int k = bas_ab_sym[h][kl][0];
                int l = bas_ab_sym[h][kl][1];
                double dum = 0.0;
                int hl = symmetry[l];
                double ** t_p = trans->pointer(hl);
                for (int pp = 0; pp < amopi_[hl]; pp++) {

                    int p = pp + pitzer_offset[hl];
                    int kp = ibas_ab_sym[h][k][p];

                    dum += inout[d2aboff[h] + ij*gems_ab[h] + kp] * t_p[l-pitzer_offset[hl] + frzcpi_[hl]+rstcpi_[hl]][pp + frzcpi_[hl]+rstcpi_[hl]];

                }
                tmp[d2aboff[h] + ij*gems_ab[h] + kl] = dum;
            }
        }
    }
    // transform third index
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            for (int kl = 0; kl < gems_ab[h]; kl++) {
                int k = bas_ab_sym[h][kl][0];
                int l = bas_ab_sym[h][kl][1];
                double dum = 0.0;
                int hk = symmetry[k];
                double ** t_p = trans->pointer(hk);
                for (int pp = 0; pp < amopi_[hk]; pp++) {

                    int p = pp + pitzer_offset[hk];
                    int pl = ibas_ab_sym[h][p][l];

                    dum += tmp[d2aboff[h] + ij*gems_ab[h] + pl] * t_p[k-pitzer_offset[hk] + frzcpi_[hk]+rstcpi_[hk]][pp + frzcpi_[hk]+rstcpi_[hk]];

                }
                inout[d2aboff[h] + ij*gems_ab[h] + kl] = dum;
            }
        }
    }
    // transform second index
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            for (int kl = 0; kl < gems_ab[h]; kl++) {
                int k = bas_ab_sym[h][kl][0];
                int l = bas_ab_sym[h][kl][1];
                double dum = 0.0;
                int hl = symmetry[l];
                double ** t_p = trans->pointer(hl);
                for (int pp = 0; pp < amopi_[hl]; pp++) {

                    int p = pp + pitzer_offset[hl];
                    int kp = ibas_ab_sym[h][k][p];

                    dum += inout[d2aboff[h] + kp*gems_ab[h] + ij] * t_p[l-pitzer_offset[hl] + frzcpi_[hl]+rstcpi_[hl]][pp + frzcpi_[hl]+rstcpi_[hl]];

                }
                tmp[d2aboff[h] + kl*gems_ab[h] + ij] = dum;
            }
        }
    }
    // transform first index
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            for (int kl = 0; kl < gems_ab[h]; kl++) {
                int k = bas_ab_sym[h][kl][0];
                int l = bas_ab_sym[h][kl][1];
                double dum = 0.0;
                int hk = symmetry[k];
                double ** t_p = trans->pointer(hk);
                for (int pp = 0; pp < amopi_[hk]; pp++) {

                    int p = pp + pitzer_offset[hk];
                    int pl = ibas_ab_sym[h][p][l];

                    dum += tmp[d2aboff[h] + pl*gems_ab[h] + ij] * t_p[k-pitzer_offset[hk] + frzcpi_[hk]+rstcpi_[hk]][pp + frzcpi_[hk]+rstcpi_[hk]];

                }
                inout[d2aboff[h] + kl*gems_ab[h] + ij] = dum;
            }
        }
    }
}



}
