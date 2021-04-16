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


// D2 portion of A^T.y ( and D1 / Q1 )
void v2RDMSolver::D2_constraints_ATu(double* A,double* u){

    // Traces
    if ( constrain_sz_ ) {
        // Tr(D2ab)
        for (int i = 0; i < amo_; i++){
            for (int j = 0; j < amo_; j++){
                int h = SymmetryPair(symmetry[i],symmetry[j]);
                int ij = ibas_ab_sym[h][i][j];
                if ( gems_ab[h] == 0 ) continue;
                A[d2aboff[h] + ij*gems_ab[h]+ij] += u[offset];
            }
        }
        offset++;

        // Tr(D2aa)
        for (int i = 0; i < amo_; i++){
            for (int j = 0; j < amo_; j++){
                if ( i==j ) continue;
                int h = SymmetryPair(symmetry[i],symmetry[j]);
                if ( gems_aa[h] == 0 ) continue;
                int ij = ibas_aa_sym[h][i][j];
                A[d2aaoff[h]+ij*gems_aa[h]+ij] += u[offset];
            }
        }
        offset++;
        // Tr(D2bb)
        for (int i = 0; i < amo_; i++){
            for (int j = 0; j < amo_; j++){
                if ( i==j ) continue;
                int h = SymmetryPair(symmetry[i],symmetry[j]);
                if ( gems_aa[h] == 0 ) continue;
                int ij = ibas_aa_sym[h][i][j];
                A[d2bboff[h]+ij*gems_aa[h]+ij] += u[offset];
            }
        }
        offset++;
    }else {
        // Tr(D2ab + D2aa + D2bb)
        for (int i = 0; i < amo_; i++){
            for (int j = 0; j < amo_; j++){
                int h = SymmetryPair(symmetry[i],symmetry[j]);

                if ( gems_ab[h] == 0 ) continue;
                int ij = ibas_ab_sym[h][i][j];
                A[d2aboff[h] + ij*gems_ab[h]+ij] += 2.0 * u[offset];

                if ( i == j ) continue;   
                if ( gems_aa[h] == 0 ) continue;
                ij = ibas_aa_sym[h][i][j];
                A[d2aaoff[h]+ij*gems_aa[h]+ij] += u[offset];
                A[d2bboff[h]+ij*gems_aa[h]+ij] += u[offset];
            }
        }
        offset++;
    }

    // hermiticity (necessary for dual SDP?)
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                double dum = u[offset];
                A[d2aaoff[h] + ij*gems_aa[h] + kl] += dum;
                A[d2aaoff[h] + kl*gems_aa[h] + ij] -= dum;
                offset++;
            }
        }
    }
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                double dum = u[offset];
                A[d2bboff[h] + ij*gems_aa[h] + kl] += dum;
                A[d2bboff[h] + kl*gems_aa[h] + ij] -= dum;
                offset++;
            }
        }
    }
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            for (int kl = 0; kl < gems_ab[h]; kl++) {
                double dum = u[offset];
                A[d2aboff[h] + ij*gems_ab[h] + kl] += dum;
                A[d2aboff[h] + kl*gems_ab[h] + ij] -= dum;
                offset++;
            }
        }
    }

    // d1 / q1 a
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                double dum = u[offset + i*amopi_[h]+j];
                A[d1aoff[h] + j*amopi_[h]+i] += dum;
                A[q1aoff[h] + i*amopi_[h]+j] += dum;
            }
        }
        offset += amopi_[h]*amopi_[h];
    }

    // d1 / q1 b
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                double dum = u[offset + i*amopi_[h]+j];
                A[d1boff[h] + j*amopi_[h]+i] += dum;
                A[q1boff[h] + i*amopi_[h]+j] += dum;
            }
        }
        offset += amopi_[h]*amopi_[h];
    }

    double na = nalpha_ - nrstc_ - nfrzc_;
    double nb = nbeta_ - nrstc_ - nfrzc_;

    int poff = 0;

    if ( !constrain_sz_ ) {

        // contraction: D2aa + D2ab -> D1 a
        poff = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < amopi_[h]; i++){
                for (int j = 0; j < amopi_[h]; j++){
                    A[d1aoff[h] + i*amopi_[h]+j] += (na + nb - 1.0) * u[offset + i*amopi_[h]+j];
                    int ii = i + poff;
                    int jj = j + poff;
                    for (int k = 0; k < amo_; k++){
                        int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                        int ik = ibas_ab_sym[h2][ii][k];
                        int jk = ibas_ab_sym[h2][jj][k];
                        A[d2aboff[h2] + ik*gems_ab[h2]+jk] -= u[offset + i*amopi_[h]+j];
                    }
                    for(int k =0; k < amo_; k++){
                        if( ii==k || jj==k )continue;
                        int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                        int ik = ibas_aa_sym[h2][ii][k];
                        int jk = ibas_aa_sym[h2][jj][k];
                        int sik = ( ii < k ? 1 : -1);
                        int sjk = ( jj < k ? 1 : -1);
                        A[d2aaoff[h2] + ik*gems_aa[h2]+jk] -= sik*sjk*u[offset + i*amopi_[h]+j];
                    }
                }
            }
            offset += amopi_[h]*amopi_[h];
            poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
        }

        // contraction: D2bb + D2ab -> D1 b
        poff = 0;
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < amopi_[h]; i++){
                for(int j = 0; j < amopi_[h]; j++){
                    A[d1boff[h] + i*amopi_[h]+j] += (na + nb - 1.0) * u[offset + i*amopi_[h]+j];
                    int ii = i + poff;
                    int jj = j + poff;
                    for(int k = 0; k < amo_; k++){
                        int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                        int ik = ibas_ab_sym[h2][k][ii];
                        int jk = ibas_ab_sym[h2][k][jj];
                        A[d2aboff[h2] + ik*gems_ab[h2]+jk] -= u[offset + i*amopi_[h]+j];
                    }
                    for(int k =0; k < amo_; k++){
                        if( ii==k || jj==k )continue;
                        int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                        int ik = ibas_aa_sym[h2][ii][k];
                        int jk = ibas_aa_sym[h2][jj][k];
                        int sik = ( ii < k ? 1 : -1);
                        int sjk = ( jj < k ? 1 : -1);
                        A[d2bboff[h2] + ik*gems_aa[h2]+jk] -= sik*sjk*u[offset + i*amopi_[h]+j];
                    }
                }
            }
            offset += amopi_[h]*amopi_[h];
            poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
        }
    }else {
        // contraction: D2ab -> D1 a
        poff = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < amopi_[h]; i++){
                for (int j = 0; j < amopi_[h]; j++){
                    A[d1aoff[h] + i*amopi_[h]+j] += nb * u[offset + i*amopi_[h]+j];
                    int ii = i + poff;
                    int jj = j + poff;
                    for (int k = 0; k < amo_; k++){
                        int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                        int ik = ibas_ab_sym[h2][ii][k];
                        int jk = ibas_ab_sym[h2][jj][k];
                        A[d2aboff[h2] + ik*gems_ab[h2]+jk] -= u[offset + i*amopi_[h]+j];
                    }
                }
            }
            offset += amopi_[h]*amopi_[h];
            poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
        }

        // contraction: D2ab -> D1 b
        poff = 0;
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < amopi_[h]; i++){
                for(int j = 0; j < amopi_[h]; j++){
                    A[d1boff[h] + i*amopi_[h]+j] += na * u[offset + i*amopi_[h]+j];
                    int ii = i + poff;
                    int jj = j + poff;
                    for(int k = 0; k < amo_; k++){
                        int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                        int ik = ibas_ab_sym[h2][k][ii];
                        int jk = ibas_ab_sym[h2][k][jj];
                        A[d2aboff[h2] + ik*gems_ab[h2]+jk] -= u[offset + i*amopi_[h]+j];
                    }
                }
            }
            offset += amopi_[h]*amopi_[h];
            poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
        }

        //contract D2aa -> D1 a
        poff = 0;
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < amopi_[h]; i++){
                for(int j = 0; j < amopi_[h]; j++){
                    A[d1aoff[h] + i*amopi_[h]+j] += (na - 1.0) * u[offset + i*amopi_[h]+j];
                    int ii = i + poff;
                    int jj = j + poff;
                    for(int k =0; k < amo_; k++){
                        if( ii==k || jj==k )continue;
                        int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                        int ik = ibas_aa_sym[h2][ii][k];
                        int jk = ibas_aa_sym[h2][jj][k];
                        int sik = ( ii < k ? 1 : -1);
                        int sjk = ( jj < k ? 1 : -1);
                        A[d2aaoff[h2] + ik*gems_aa[h2]+jk] -= sik*sjk*u[offset + i*amopi_[h]+j];
                    }
                }
            }
            offset += amopi_[h]*amopi_[h];
            poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
        }

        //contract D2bb -> D1 b
        poff = 0;
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < amopi_[h]; i++){
                for(int j = 0; j < amopi_[h]; j++){
                    A[d1boff[h] + i*amopi_[h]+j] += (nb - 1.0) * u[offset + i*amopi_[h]+j];
                    int ii = i + poff;
                    int jj = j + poff;
                    for(int k =0; k < amo_; k++){
                        if( ii==k || jj==k )continue;
                        int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                        int ik = ibas_aa_sym[h2][ii][k];
                        int jk = ibas_aa_sym[h2][jj][k];
                        int sik = ( ii < k ? 1 : -1);
                        int sjk = ( jj < k ? 1 : -1);
                        A[d2bboff[h2] + ik*gems_aa[h2]+jk] -= sik*sjk*u[offset + i*amopi_[h]+j];
                    }
                }
            }
            offset += amopi_[h]*amopi_[h];
            poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
        }
    }

}

// D2 portion of A.x (and D1/Q1)
void v2RDMSolver::D2_constraints_Au(double* A,double* u){

    if ( constrain_sz_ ) {
        // Traces
        // Tr(D2ab)
        double sumab =0.0;
        for (int i = 0; i < amo_; i++){
            for (int j = 0; j < amo_; j++){
                int h = SymmetryPair(symmetry[i],symmetry[j]);
                if ( gems_ab[h] == 0 ) continue;
                int ij = ibas_ab_sym[h][i][j];
                sumab += u[d2aboff[h] + ij*gems_ab[h]+ij];
            }
        }
        A[offset] = sumab;
        offset++;

        // Tr(D2aa)
        double sumaa =0.0;
        for (int i = 0; i < amo_; i++){
            for (int j = 0; j < amo_; j++){
                if ( i==j ) continue;
                int h = SymmetryPair(symmetry[i],symmetry[j]);
                if ( gems_aa[h] == 0 ) continue;
                int ij = ibas_aa_sym[h][i][j];
                sumaa += u[d2aaoff[h] + ij*gems_aa[h]+ij];
            }

        }
        A[offset] = sumaa;
        offset++;

        // Tr(D2bb)
        double sumbb =0.0;
        for (int i = 0; i < amo_; i++){
            for (int j = 0; j < amo_; j++){
                if ( i==j ) continue;
                int h = SymmetryPair(symmetry[i],symmetry[j]);
                if ( gems_aa[h] == 0 ) continue;
                int ij = ibas_aa_sym[h][i][j];
                sumbb += u[d2bboff[h] + ij*gems_aa[h]+ij];
            }

        }
        A[offset] = sumbb;
        offset++;
    }else {
        // Tr(D2ab + D2aa + D2bb)
        double sum =0.0;
        for (int i = 0; i < amo_; i++){
            for (int j = 0; j < amo_; j++){
                int h = SymmetryPair(symmetry[i],symmetry[j]);

                if ( gems_ab[h] == 0 ) continue;
                int ij = ibas_ab_sym[h][i][j];
                sum += 2.0 * u[d2aboff[h] + ij*gems_ab[h]+ij];

                if ( i==j ) continue;

                if ( gems_aa[h] == 0 ) continue;
                ij = ibas_aa_sym[h][i][j];
                sum += u[d2aaoff[h] + ij*gems_aa[h]+ij];
                sum += u[d2bboff[h] + ij*gems_aa[h]+ij];
            }

        }
        A[offset] = sum;
        offset++;
    }

    // hermiticity (necessary for dual SDP?)
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                double dum = u[d2aaoff[h] + ij*gems_aa[h] + kl];
                dum       -= u[d2aaoff[h] + kl*gems_aa[h] + ij];
                A[offset] = dum;
                offset++;
            }
        }
    }
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                double dum = u[d2bboff[h] + ij*gems_aa[h] + kl];
                dum       -= u[d2bboff[h] + kl*gems_aa[h] + ij];
                A[offset] = dum;
                offset++;
            }
        }
    }
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            for (int kl = 0; kl < gems_ab[h]; kl++) {
                double dum = u[d2aboff[h] + ij*gems_ab[h] + kl];
                dum       -= u[d2aboff[h] + kl*gems_ab[h] + ij];
                A[offset] = dum;
                offset++;
            }
        }
    }

    // d1 / q1 a
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                A[offset+i*amopi_[h]+j] = u[d1aoff[h]+j*amopi_[h]+i] + u[q1aoff[h]+i*amopi_[h]+j];
            }
        }
        offset += amopi_[h]*amopi_[h];
    }

    // d1 / q1 b
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                A[offset+i*amopi_[h]+j] = u[d1boff[h]+j*amopi_[h]+i] + u[q1boff[h]+i*amopi_[h]+j];
            }
        }
        offset += amopi_[h]*amopi_[h];
    }

    double na = nalpha_ - nrstc_ - nfrzc_;
    double nb = nbeta_ - nrstc_ - nfrzc_;
    double n = na + nb;

    int poff = 0;

    if ( !constrain_sz_ ) {
        // contraction: D2aa + D2ab -> D1 a
        poff = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < amopi_[h]; i++){
                for (int j = 0; j < amopi_[h]; j++){
                    double sum = (n - 1.0) * u[d1aoff[h] + i*amopi_[h]+j];
                    int ii  = i + poff;
                    int jj  = j + poff;
                    for(int k = 0; k < amo_; k++){
                        int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                        int ik = ibas_ab_sym[h2][ii][k];
                        int jk = ibas_ab_sym[h2][jj][k];
                        sum -= u[d2aboff[h2] + ik*gems_ab[h2]+jk];
                    }
                    for(int k = 0; k < amo_; k++){
                        if( ii==k || jj==k ) continue;
                        int h2   = SymmetryPair(symmetry[ii],symmetry[k]);
                        int ik  = ibas_aa_sym[h2][ii][k];
                        int jk  = ibas_aa_sym[h2][jj][k];
                        int sik = ( ii < k ) ? 1 : -1;
                        int sjk = ( jj < k ) ? 1 : -1;
                        sum -= sik*sjk*u[d2aaoff[h2] + ik*gems_aa[h2]+jk];
                    }
                    A[offset + i*amopi_[h]+j] = sum;
                }
            }
            offset += amopi_[h]*amopi_[h];
            poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
        }

        // contraction: D2bb + D2ab -> D1 b
        poff = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < amopi_[h]; i++){
                for (int j = 0; j < amopi_[h]; j++){
                    double sum = (n - 1.0) * u[d1boff[h] + i*amopi_[h]+j];
                    int ii  = i + poff;
                    int jj  = j + poff;
                    for(int k = 0; k < amo_; k++){
                        int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                        int ik = ibas_ab_sym[h2][k][ii];
                        int jk = ibas_ab_sym[h2][k][jj];
                        sum -= u[d2aboff[h2] + ik*gems_ab[h2]+jk];
                    }
                    for(int k = 0; k < amo_; k++){
                        if( ii==k || jj==k ) continue;
                        int h2   = SymmetryPair(symmetry[ii],symmetry[k]);
                        int ik  = ibas_aa_sym[h2][ii][k];
                        int jk  = ibas_aa_sym[h2][jj][k];
                        int sik = ( ii < k ) ? 1 : -1;
                        int sjk = ( jj < k ) ? 1 : -1;
                        sum -= sik*sjk*u[d2bboff[h2] + ik*gems_aa[h2]+jk];
                    }
                    A[offset + i*amopi_[h]+j] = sum;
                }
            }
            offset += amopi_[h]*amopi_[h];
            poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
        }
    }else {
        // contraction: D2ab -> D1 a
        poff = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < amopi_[h]; i++){
                for (int j = 0; j < amopi_[h]; j++){
                    double sum = nb * u[d1aoff[h] + i*amopi_[h]+j];
                    int ii  = i + poff;
                    int jj  = j + poff;
                    for(int k = 0; k < amo_; k++){
                        int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                        int ik = ibas_ab_sym[h2][ii][k];
                        int jk = ibas_ab_sym[h2][jj][k];
                        sum -= u[d2aboff[h2] + ik*gems_ab[h2]+jk];
                    }
                    A[offset + i*amopi_[h]+j] = sum;
                }
            }
            offset += amopi_[h]*amopi_[h];
            poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
        }

        // contraction: D2ab -> D1 b
        poff = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < amopi_[h]; i++){
                for (int j = 0; j < amopi_[h]; j++){
                    double sum = na * u[d1boff[h] + i*amopi_[h]+j];
                    int ii  = i + poff;
                    int jj  = j + poff;
                    for(int k = 0; k < amo_; k++){
                        int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                        int ik = ibas_ab_sym[h2][k][ii];
                        int jk = ibas_ab_sym[h2][k][jj];
                        sum -= u[d2aboff[h2] + ik*gems_ab[h2]+jk];
                    }
                    A[offset + i*amopi_[h]+j] = sum;
                }
            }
            offset += amopi_[h]*amopi_[h];
            poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
        }

        //contract D2aa -> D1 a
        poff = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < amopi_[h]; i++){
                for (int j = 0; j < amopi_[h]; j++){
                    double sum = (na - 1.0) * u[d1aoff[h] + i*amopi_[h]+j];
                    int ii  = i + poff;
                    int jj  = j + poff;
                    for(int k = 0; k < amo_; k++){
                        if( ii==k || jj==k ) continue;
                        int h2   = SymmetryPair(symmetry[ii],symmetry[k]);
                        int ik  = ibas_aa_sym[h2][ii][k];
                        int jk  = ibas_aa_sym[h2][jj][k];
                        int sik = ( ii < k ) ? 1 : -1;
                        int sjk = ( jj < k ) ? 1 : -1;
                        sum -= sik*sjk*u[d2aaoff[h2] + ik*gems_aa[h2]+jk];
                    }
                    A[offset+i*amopi_[h]+j] = sum;
                }
            }
            offset += amopi_[h]*amopi_[h];
            poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
        }

        //contract D2bb -> D1 b
        poff = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < amopi_[h]; i++){
                for (int j = 0; j < amopi_[h]; j++){
                    double sum = (nb - 1.0) * u[d1boff[h] + i*amopi_[h]+j];
                    int ii  = i + poff;
                    int jj  = j + poff;
                    for(int k = 0; k < amo_; k++){
                        if( ii==k || jj==k ) continue;
                        int h2   = SymmetryPair(symmetry[ii],symmetry[k]);
                        int ik  = ibas_aa_sym[h2][ii][k];
                        int jk  = ibas_aa_sym[h2][jj][k];
                        int sik = ( ii < k ) ? 1 : -1;
                        int sjk = ( jj < k ) ? 1 : -1;
                        sum -= sik*sjk*u[d2bboff[h2] + ik*gems_aa[h2]+jk];
                    }
                    A[offset+i*amopi_[h]+j] = sum;
                }
            }
            offset += amopi_[h]*amopi_[h];
            poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
        }
    }

}

}
