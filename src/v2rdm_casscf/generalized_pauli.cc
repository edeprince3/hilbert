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
void v2RDMSolver::Generalized_Pauli_constraints_ATu(SharedVector A,SharedVector u){

    if ( gpconstraint_ == GeneralizedPauli_3_8 ) {

        Generalized_Pauli_3_8_constraints_ATu(A,u);

    }else if ( gpconstraint_ == GeneralizedPauli_4_8 ) {

        Generalized_Pauli_4_8_constraints_ATu(A,u);

    }else if ( gpconstraint_ == GeneralizedPauli_5_8 ) {

        Generalized_Pauli_3_8_constraints_ATu(A,u);

    }else if ( gpconstraint_ == GeneralizedPauli_3_6 ) {

        Generalized_Pauli_3_6_constraints_ATu(A,u);

    }else if ( gpconstraint_ == GeneralizedPauli_4_10 ) {

        Generalized_Pauli_4_10_constraints_ATu(A,u);

    }else if ( gpconstraint_ == GeneralizedPauli_6_10 ) {

        Generalized_Pauli_4_10_constraints_ATu(A,u);

    }else if ( gpconstraint_ == GeneralizedPauli_5_10 ) {

        Generalized_Pauli_5_10_constraints_ATu(A,u);

    }else if ( gpconstraint_ == GeneralizedPauli_3_10 ) {

        Generalized_Pauli_3_10_constraints_ATu(A,u);

    }else if ( gpconstraint_ == GeneralizedPauli_7_10 ) {

        Generalized_Pauli_3_10_constraints_ATu(A,u);

    }

}

// portion of A.x corresponding to generalized pauli constraints
void v2RDMSolver::Generalized_Pauli_constraints_Au(SharedVector A,SharedVector u){

    if ( gpconstraint_ == GeneralizedPauli_3_8 ) {

        Generalized_Pauli_3_8_constraints_Au(A,u);

    }else if ( gpconstraint_ == GeneralizedPauli_4_8 ) {

        Generalized_Pauli_4_8_constraints_Au(A,u);

    }else if ( gpconstraint_ == GeneralizedPauli_5_8 ) {

        Generalized_Pauli_3_8_constraints_Au(A,u);

    }else if ( gpconstraint_ == GeneralizedPauli_3_6 ) {

        Generalized_Pauli_3_6_constraints_Au(A,u);

    }else if ( gpconstraint_ == GeneralizedPauli_4_10 ) {

        Generalized_Pauli_4_10_constraints_Au(A,u);

    }else if ( gpconstraint_ == GeneralizedPauli_6_10 ) {

        Generalized_Pauli_4_10_constraints_Au(A,u);

    }else if ( gpconstraint_ == GeneralizedPauli_5_10 ) {

        Generalized_Pauli_5_10_constraints_Au(A,u);

    }else if ( gpconstraint_ == GeneralizedPauli_3_10 ) {

        Generalized_Pauli_3_10_constraints_Au(A,u);

    }else if ( gpconstraint_ == GeneralizedPauli_7_10 ) {

        Generalized_Pauli_3_10_constraints_Au(A,u);

    }

}

void v2RDMSolver::Generalized_Pauli_ATu_term(double val, double ** orbs,double * A,int * offa, int * offb,int index) {
    for (int i = 0; i < amo_; i++) {
        int hi = symmetry[i];
        int ii = i - pitzer_offset[hi];
        for (int j = 0; j < amo_; j++) {
            int hj = symmetry[j];
            if ( hi != hj ) continue;
            int jj = j - pitzer_offset[hj];

            A[offa[hi] + ii*amopi_[hi] + jj] += val * orbs[i     ][index - 1] * orbs[j     ][index - 1];
            A[offb[hi] + ii*amopi_[hi] + jj] += val * orbs[i+amo_][index - 1] * orbs[j+amo_][index - 1];

        }
    }
}

double v2RDMSolver::Generalized_Pauli_Au_term(double ** orbs,double * u,int * offa, int * offb,int index) {
    double dum = 0.0;
    for (int i = 0; i < amo_; i++) {
        int hi = symmetry[i];
        int ii = i - pitzer_offset[hi];
        for (int j = 0; j < amo_; j++) {
            int hj = symmetry[j];
            if ( hi != hj ) continue;
            int jj = j - pitzer_offset[hj];
            dum += orbs[i     ][index - 1] * orbs[j     ][index - 1] * u[offa[hi] + ii*amopi_[hi] + jj];
            dum += orbs[i+amo_][index - 1] * orbs[j+amo_][index - 1] * u[offb[hi] + ii*amopi_[hi] + jj];
        }
    }
    return dum;
}

void v2RDMSolver::SortedNaturalOrbitals() {

    //NatOrbs_->zero();
    std::shared_ptr<Matrix> temp (new Matrix(2*amo_,2*amo_));

    int * x1aoff;
    int * x1boff;

    if ( gpconstraint_ == GeneralizedPauli_5_8 || gpconstraint_ == GeneralizedPauli_6_10 || gpconstraint_ == GeneralizedPauli_7_10 ) {
        x1aoff = q1aoff;
        x1boff = q1boff;
    }else {
        x1aoff = d1aoff;
        x1boff = d1boff;
    }

    std::shared_ptr<Matrix> Da (new Matrix(nirrep_,amopi_,amopi_));
    std::shared_ptr<Matrix> eigveca (new Matrix(nirrep_,amopi_,amopi_));
    std::shared_ptr<Vector> eigvala (new Vector("Natural Orbital Occupation Numbers (alpha)",nirrep_,amopi_));
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h] ; i++) {
            for (int j = 0; j < amopi_[h]; j++) {
                Da->pointer(h)[i][j] = x->pointer()[x1aoff[h]+i*amopi_[h]+j];
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
                Db->pointer(h)[i][j] = x->pointer()[x1boff[h]+i*amopi_[h]+j];
            }
        }
    }
    Db->diagonalize(eigvecb,eigvalb,descending);

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
    double **old = NatOrbs_->pointer();
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
