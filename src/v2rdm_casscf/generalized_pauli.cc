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

#include "v2rdm_solver.h"
#include <misc/omp.h>

using namespace psi;

namespace hilbert {


// portion of A^T.y corresponding to generalized pauli constraints
void v2RDMSolver::Generalized_Pauli_constraints_ATu(double* A,double* u, int state){

    if ( gpc_[state] == GeneralizedPauli_3_8 ) {

        Generalized_Pauli_3_8_constraints_ATu(A,u,state);

    }else if ( gpc_[state] == GeneralizedPauli_4_8 ) {

        Generalized_Pauli_4_8_constraints_ATu(A,u,state);

    }else if ( gpc_[state] == GeneralizedPauli_5_8 ) {

        Generalized_Pauli_3_8_constraints_ATu(A,u,state);

    }else if ( gpc_[state] == GeneralizedPauli_3_6 ) {

        Generalized_Pauli_3_6_constraints_ATu(A,u,state);

    }else if ( gpc_[state] == GeneralizedPauli_4_10 ) {

        Generalized_Pauli_4_10_constraints_ATu(A,u,state);

    }else if ( gpc_[state] == GeneralizedPauli_6_10 ) {

        Generalized_Pauli_4_10_constraints_ATu(A,u,state);

    }else if ( gpc_[state] == GeneralizedPauli_5_10 ) {

        Generalized_Pauli_5_10_constraints_ATu(A,u,state);

    }else if ( gpc_[state] == GeneralizedPauli_3_10 ) {

        Generalized_Pauli_3_10_constraints_ATu(A,u,state);

    }else if ( gpc_[state] == GeneralizedPauli_7_10 ) {

        Generalized_Pauli_3_10_constraints_ATu(A,u,state);

    }

}

// portion of A.x corresponding to generalized pauli constraints
void v2RDMSolver::Generalized_Pauli_constraints_Au(double* A,double* u, int state){

    if ( gpc_[state] == GeneralizedPauli_3_8 ) {

        Generalized_Pauli_3_8_constraints_Au(A,u,state);

    }else if ( gpc_[state] == GeneralizedPauli_4_8 ) {

        Generalized_Pauli_4_8_constraints_Au(A,u,state);

    }else if ( gpc_[state] == GeneralizedPauli_5_8 ) {

        Generalized_Pauli_3_8_constraints_Au(A,u,state);

    }else if ( gpc_[state] == GeneralizedPauli_3_6 ) {

        Generalized_Pauli_3_6_constraints_Au(A,u,state);

    }else if ( gpc_[state] == GeneralizedPauli_4_10 ) {

        Generalized_Pauli_4_10_constraints_Au(A,u,state);

    }else if ( gpc_[state] == GeneralizedPauli_6_10 ) {

        Generalized_Pauli_4_10_constraints_Au(A,u,state);

    }else if ( gpc_[state] == GeneralizedPauli_5_10 ) {

        Generalized_Pauli_5_10_constraints_Au(A,u,state);

    }else if ( gpc_[state] == GeneralizedPauli_3_10 ) {

        Generalized_Pauli_3_10_constraints_Au(A,u,state);

    }else if ( gpc_[state] == GeneralizedPauli_7_10 ) {

        Generalized_Pauli_3_10_constraints_Au(A,u,state);

    }

}

void v2RDMSolver::Generalized_Pauli_ATu_term(int *** sign_a, int *** sign_b, double val, double ** orbs,double * A,int *** map_a, int *** map_b,int index) {

    for (int h = 0; h < nirrep_; h++) {
        for (int ii = 0; ii < amopi_[h]; ii++) {
            int i = ii + pitzer_offset[h];
            for (int jj = 0; jj < amopi_[h]; jj++) {
                int j = jj + pitzer_offset[h];

                int id_a = map_a[h][ii][jj];
                int id_b = map_b[h][ii][jj];

                int sg_a = sign_a[h][ii][jj];
                int sg_b = sign_b[h][ii][jj];

                if ( id_a >= 0 ) {
                    A[ id_a ] += val * orbs[i     ][index - 1] * orbs[j     ][index - 1] * sg_a;
                }
                if ( id_b >= 0 ) {
                    A[ id_b ] += val * orbs[i+amo_][index - 1] * orbs[j+amo_][index - 1] * sg_b;
                }
            }
        }
    }

}

double v2RDMSolver::Generalized_Pauli_Au_term(double ** orbs,double * u,int *** map_a, int *** map_b,int index, double rdm_nrm, int *** sign_a, int *** sign_b) {

    double dum = 0.0;
    for (int h = 0; h < nirrep_; h++) {
        for (int ii = 0; ii < amopi_[h]; ii++) {
            int i = ii + pitzer_offset[h];
            for (int jj = 0; jj < amopi_[h]; jj++) {
                int j = jj + pitzer_offset[h];

                int id_a = map_a[h][ii][jj];
                int id_b = map_b[h][ii][jj];

                int sg_a = sign_a[h][ii][jj];
                int sg_b = sign_b[h][ii][jj];

                if ( id_a >= 0 ) {
                    dum += orbs[i     ][index - 1] * orbs[j     ][index - 1] * u[ id_a ] * sg_a;
                }
                if ( id_b >= 0 ) {
                    dum += orbs[i+amo_][index - 1] * orbs[j+amo_][index - 1] * u[ id_b ] * sg_b;
                }

            }
        }
    }

    return dum / rdm_nrm;
}

void v2RDMSolver::set_gpc_rdm_nrm() {

    gpc_rdm_nrm_.clear();

    if ( constrain_gpc_1rdm_ ) {
        gpc_rdm_nrm_.push_back(1.0);
    }

    if ( constrain_gpc_2rdm_ ) {
        throw PsiException("GPCs may not yet be applied to the 2RDM",__FILE__,__LINE__);
    }

}

void v2RDMSolver::SortedNaturalOrbitals(int state) {

    std::shared_ptr<Matrix> temp (new Matrix(2*amo_,2*amo_));

    //int * x1aoff;
    //int * x1boff;

    //if ( gpc_[state] == GeneralizedPauli_5_8 || gpc_[state] == GeneralizedPauli_6_10 || gpc_[state] == GeneralizedPauli_7_10 ) {
    //    x1aoff = q1aoff;
    //    x1boff = q1boff;
    //}else {
    //    x1aoff = d1aoff;
    //    x1boff = d1boff;
    //}

    std::shared_ptr<Matrix> Da (new Matrix(nirrep_,amopi_,amopi_));
    std::shared_ptr<Matrix> eigveca (new Matrix(nirrep_,amopi_,amopi_));
    std::shared_ptr<Vector> eigvala (new Vector("Natural Orbital Occupation Numbers (alpha)",nirrep_,amopi_));
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h] ; i++) {
            for (int j = 0; j < amopi_[h]; j++) {

                int id_a = gpc_rdm_map_a_[state][h][i][j];

                if ( id_a < 0 ) {
                    Da->pointer(h)[i][j] = 0.0;
                }else {
                    Da->pointer(h)[i][j] = x->pointer()[id_a] / gpc_rdm_nrm_[state] * gpc_rdm_sign_a_[state][h][i][j];
                }

            }
        }
    }
    Da->diagonalize(eigveca,eigvala,descending);

    std::shared_ptr<Matrix> Db (new Matrix(nirrep_,amopi_,amopi_));
    std::shared_ptr<Matrix> eigvecb (new Matrix(nirrep_,amopi_,amopi_));
    std::shared_ptr<Vector> eigvalb (new Vector("Natural Orbital Occupation Numbers (beta)",nirrep_,amopi_));
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h] ; i++) {
            for (int j = 0; j < amopi_[h]; j++) {

                int id_b = gpc_rdm_map_b_[state][h][i][j];

                if ( id_b < 0 ) {
                    Db->pointer(h)[i][j] = 0.0;
                }else {
                    Db->pointer(h)[i][j] = x->pointer()[id_b] / gpc_rdm_nrm_[state] * gpc_rdm_sign_b_[state][h][i][j];
                }

            }
        }
    }
    Db->diagonalize(eigvecb,eigvalb,descending);

    // TEST
    if ( print_gpc_error_ ) {
        eigvala->print();
        eigvalb->print();
    }

    // sort!
    int * skipa = (int*)malloc(amo_*sizeof(int));
    int * skipb = (int*)malloc(amo_*sizeof(int));
    memset((void*)skipa,'\0',amo_*sizeof(int));
    memset((void*)skipb,'\0',amo_*sizeof(int));
    double ** orb_p = temp->pointer();

    for (int i = 0; i < 2*amo_; i++) {
        bool is_alpha = true;
        int  imax = 0;
        double max = -999;
        for (int h = 0; h < nirrep_; h++) {
            for (int j = 0; j < amopi_[h]; j++) {
                int jj = j + pitzer_offset[h];
                if ( skipa[jj] ) continue;
                if ( eigvala->pointer(h)[j] > max ) {
                    imax = jj;
                    max = eigvala->pointer(h)[j];
                    is_alpha = true;
                }
            }
        }
        for (int h = 0; h < nirrep_; h++) {
            for (int j = 0; j < amopi_[h]; j++) {
                int jj = j + pitzer_offset[h];
                if ( skipb[jj] ) continue;
                if ( eigvalb->pointer(h)[j] > max ) {
                    imax = jj;
                    max = eigvalb->pointer(h)[j];
                    is_alpha = false;
                }
            }
        }
        //printf("%s %i\n",is_alpha ? "alpha" : "beta", imax);
        if ( is_alpha ) {
            int h = symmetry[imax];
            double ** e_p = eigveca->pointer(h);
            skipa[imax] = 1;
            for (int j = 0; j < amopi_[h]; j++) {
                int jj = j + pitzer_offset[h];
                orb_p[jj][i] = e_p[j][imax-pitzer_offset[h]];
            }
        }else {
            int h = symmetry[imax];
            double ** e_p = eigvecb->pointer(h);
            skipb[imax] = 1;
            for (int j = 0; j < amopi_[h]; j++) {
                int jj = j + pitzer_offset[h];
                orb_p[jj+amo_][i] = e_p[j][imax-pitzer_offset[h]];
            }
        }
    }

    free(skipa);
    free(skipb);

    // check phase:
    double **old = NatOrbs_[state]->pointer();
    for (int i = 0; i < 2*amo_; i++) {
        double max = -999;
        int jmax   = 0;
        for (int j = 0; j < 2*amo_; j++) {
            if ( fabs(orb_p[j][i]) > max ) {
                max = fabs(orb_p[j][i]);
                jmax = j;
            }
        }
        int sg = 1;
        if ( orb_p[jmax][i] * old[jmax][i] < 0.0 ) {
            sg = -1;
        }
        for (int j = 0; j < 2*amo_; j++) {
             old[j][i] = sg * orb_p[j][i];
        }
    }

    //NatOrbs_->print();
}

}
