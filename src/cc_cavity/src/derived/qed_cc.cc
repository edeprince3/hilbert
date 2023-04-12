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

#include "../../include/derived/qed_cc.h"
#include <tiledarray.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsio/psio.hpp>
#include <utility>

#include <psi4/libtrans/integraltransform.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libfunctional/superfunctional.h>
#include "psi4/libscf_solver/uhf.h"
#include <psi4/libmints/mintshelper.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/basisset.h>
#include <psi4/lib3index/dftensor.h>
#include <psi4/libqt/qt.h>

#include "../../misc/ta_helper.h"
#include "../../misc/threeindexintegralsta.h"

#include <mkl.h>
#include <omp.h>
#include "../../misc/qed_blas.h"
#include "../../../misc/hilbert_psifiles.h"
#include "../../../polaritonic_scf/uhf.h"
#include <unistd.h>
#include <psi4/psifiles.h>

namespace hilbert {

    QED_CC::QED_CC(const shared_ptr<Wavefunction> &reference_wavefunction, Options &options) :
                             CC_Cavity(reference_wavefunction, options) {
        // update cavity terms once more
        if ( n_photon_states_ > 1) {
            update_cavity_terms();
        }
        same_a_b_orbs_ = same_a_b_dens_ = false;
        Printf("\n  Spin-tracing enabled.\n");
        Printf(
                options_.get_bool("QED_USE_RELAXED_ORBITALS") ? "  Using relaxed orbitals.\n"
                                                              : "  Using unrelaxed orbitals.\n"
              );
        Printf("\n");

        cc_type_ = has_photon_ ? "QED-CC" : "CC";
        cc_type_ += "SD";
        if (include_t3_ || include_u3_) cc_type_ += "T";
        if (include_t4_ || include_u4_) cc_type_ += "Q";
        if (has_photon_) cc_type_ += "-1";

        // set coupling factors
        lambda_[0] = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
        lambda_[1] = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
        lambda_[2] = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

        // make containers for amplitudes
        init_operators();

        // make containers for integrals
        init_integrals();

        // initialize orbital energies
        epsilon_ = (double*) calloc(2*nso_, sizeof(double));

        // transform integrals to MO basis
        transform_integrals(false);
    }

    void QED_CC::init_operators() {

        // call base class for t1 and t2, plus some other stuff
        CC_Cavity::init_operators(); 

        /// initialize amplitude and residual blocks

        // t3
        if (include_t3_) {
            amplitudes_["t3_aaaaaa"] = HelperD::makeTensor(world_, {va_,va_,va_, oa_,oa_,oa_}, true);
            amplitudes_["t3_aabaab"] = HelperD::makeTensor(world_, {va_,va_,vb_, oa_,oa_,ob_}, true);
            amplitudes_["t3_abbabb"] = HelperD::makeTensor(world_, {va_,vb_,vb_, oa_,ob_,ob_}, true);
            amplitudes_["t3_bbbbbb"] = HelperD::makeTensor(world_, {vb_,vb_,vb_, ob_,ob_,ob_}, true);
            residuals_["t3_aaaaaa"] = HelperD::makeTensor(world_, {va_,va_,va_, oa_,oa_,oa_}, true);
            residuals_["t3_aabaab"] = HelperD::makeTensor(world_, {va_,va_,vb_, oa_,oa_,ob_}, true);
            residuals_["t3_abbabb"] = HelperD::makeTensor(world_, {va_,vb_,vb_, oa_,ob_,ob_}, true);
            residuals_["t3_bbbbbb"] = HelperD::makeTensor(world_, {vb_,vb_,vb_, ob_,ob_,ob_}, true);
        }

        // t4
        if (include_t4_) {
            amplitudes_["t4_aaaaaaaa"] = HelperD::makeTensor(world_, {va_,va_,va_,va_, oa_,oa_,oa_,oa_}, true);
            amplitudes_["t4_aaabaaab"] = HelperD::makeTensor(world_, {va_,va_,va_,vb_, oa_,oa_,oa_,ob_}, true);
            amplitudes_["t4_aabbaabb"] = HelperD::makeTensor(world_, {va_,va_,vb_,vb_, oa_,oa_,ob_,ob_}, true);
            amplitudes_["t4_abbbabbb"] = HelperD::makeTensor(world_, {va_,vb_,vb_,vb_, oa_,ob_,ob_,ob_}, true);
            amplitudes_["t4_bbbbbbbb"] = HelperD::makeTensor(world_, {vb_,vb_,vb_,vb_, ob_,ob_,ob_,ob_}, true);
            residuals_["t4_aaaaaaaa"] = HelperD::makeTensor(world_, {va_,va_,va_,va_, oa_,oa_,oa_,oa_}, true);
            residuals_["t4_aaabaaab"] = HelperD::makeTensor(world_, {va_,va_,va_,vb_, oa_,oa_,oa_,ob_}, true);
            residuals_["t4_aabbaabb"] = HelperD::makeTensor(world_, {va_,va_,vb_,vb_, oa_,oa_,ob_,ob_}, true);
            residuals_["t4_abbbabbb"] = HelperD::makeTensor(world_, {va_,vb_,vb_,vb_, oa_,ob_,ob_,ob_}, true);
            residuals_["t4_bbbbbbbb"] = HelperD::makeTensor(world_, {vb_,vb_,vb_,vb_, ob_,ob_,ob_,ob_}, true);
        }

        // u0
        scalar_amps_["u0"] = 0.0;
        scalar_resids_["u0"] = 0.0;

        // u1
        if (include_u1_) {
            amplitudes_["u1_aa"]  = HelperD::makeTensor(world_, {va_, oa_}, true);
            amplitudes_["u1_bb"]  = HelperD::makeTensor(world_, {vb_, ob_}, true);
            residuals_["u1_aa"]  = HelperD::makeTensor(world_, {va_, oa_}, true);
            residuals_["u1_bb"]  = HelperD::makeTensor(world_, {vb_, ob_}, true);
        }

        // u2
        if (include_u2_) {
            amplitudes_["u2_aaaa"]  = HelperD::makeTensor(world_, {va_,va_, oa_,oa_}, true);
            amplitudes_["u2_abab"]  = HelperD::makeTensor(world_, {va_,vb_, oa_,ob_}, true);
            amplitudes_["u2_bbbb"]  = HelperD::makeTensor(world_, {vb_,vb_, ob_,ob_}, true);
            residuals_["u2_aaaa"]  = HelperD::makeTensor(world_, {va_,va_, oa_,oa_}, true);
            residuals_["u2_abab"]  = HelperD::makeTensor(world_, {va_,vb_, oa_,ob_}, true);
            residuals_["u2_bbbb"]  = HelperD::makeTensor(world_, {vb_,vb_, ob_,ob_}, true);
        }

        // u3
        if (include_u3_) {
            amplitudes_["u3_aaaaaa"] = HelperD::makeTensor(world_, {va_,va_,va_, oa_,oa_,oa_}, true);
            amplitudes_["u3_aabaab"] = HelperD::makeTensor(world_, {va_,va_,vb_, oa_,oa_,ob_}, true);
            amplitudes_["u3_abbabb"] = HelperD::makeTensor(world_, {va_,vb_,vb_, oa_,ob_,ob_}, true);
            amplitudes_["u3_bbbbbb"] = HelperD::makeTensor(world_, {vb_,vb_,vb_, ob_,ob_,ob_}, true);
            residuals_["u3_aaaaaa"] = HelperD::makeTensor(world_, {va_,va_,va_, oa_,oa_,oa_}, true);
            residuals_["u3_aabaab"] = HelperD::makeTensor(world_, {va_,va_,vb_, oa_,oa_,ob_}, true);
            residuals_["u3_abbabb"] = HelperD::makeTensor(world_, {va_,vb_,vb_, oa_,ob_,ob_}, true);
            residuals_["u3_bbbbbb"] = HelperD::makeTensor(world_, {vb_,vb_,vb_, ob_,ob_,ob_}, true);
        }

        // u4
        if (include_u4_) {
            amplitudes_["u4_aaaaaaaa"] = HelperD::makeTensor(world_, {va_,va_,va_,va_, oa_,oa_,oa_,oa_}, true);
            amplitudes_["u4_aaabaaab"] = HelperD::makeTensor(world_, {va_,va_,va_,vb_, oa_,oa_,oa_,ob_}, true);
            amplitudes_["u4_aabbaabb"] = HelperD::makeTensor(world_, {va_,va_,vb_,vb_, oa_,oa_,ob_,ob_}, true);
            amplitudes_["u4_abbbabbb"] = HelperD::makeTensor(world_, {va_,vb_,vb_,vb_, oa_,ob_,ob_,ob_}, true);
            amplitudes_["u4_bbbbbbbb"] = HelperD::makeTensor(world_, {vb_,vb_,vb_,vb_, ob_,ob_,ob_,ob_}, true);
            residuals_["u4_aaaaaaaa"] = HelperD::makeTensor(world_, {va_,va_,va_,va_, oa_,oa_,oa_,oa_}, true);
            residuals_["u4_aaabaaab"] = HelperD::makeTensor(world_, {va_,va_,va_,vb_, oa_,oa_,oa_,ob_}, true);
            residuals_["u4_aabbaabb"] = HelperD::makeTensor(world_, {va_,va_,vb_,vb_, oa_,oa_,ob_,ob_}, true);
            residuals_["u4_abbbabbb"] = HelperD::makeTensor(world_, {va_,vb_,vb_,vb_, oa_,ob_,ob_,ob_}, true);
            residuals_["u4_bbbbbbbb"] = HelperD::makeTensor(world_, {vb_,vb_,vb_,vb_, ob_,ob_,ob_,ob_}, true);
        }
    }

    void QED_CC::print_dimensions() {
        
        /// compute the number of amplitudes
        singleDim_ = oa_ * va_
                     + ob_ * vb_;
        doubleDim_ = oa_ * oa_ * va_ * va_
                     + ob_ * ob_ * vb_ * vb_
                     + oa_ * ob_ * va_ * vb_;
        if (include_t3_ || include_u3_)
            tripleDim_ = oa_ * oa_ * oa_ * va_ * va_ * va_
                         + ob_ * ob_ * ob_ * vb_ * vb_ * vb_
                         + oa_ * oa_ * ob_ * va_ * va_ * vb_
                         + oa_ * ob_ * ob_ * va_ * vb_ * vb_;
        if (include_t4_ || include_u4_)
            quadDim_ = oa_ * oa_ * oa_ * oa_ * va_ * va_ * va_ * va_
                       + ob_ * ob_ * ob_ * ob_ * vb_ * vb_ * vb_ * vb_
                       + oa_ * oa_ * oa_ * ob_ * va_ * va_ * va_ * vb_
                       + oa_ * oa_ * ob_ * ob_ * va_ * va_ * vb_ * vb_
                       + oa_ * ob_ * ob_ * ob_ * va_ * vb_ * vb_ * vb_;

        /// evaluate the total number of amplitudes
        size_t ccamps_dim_ = singleDim_ + doubleDim_; // t1, t2
        if (include_t3_) ccamps_dim_ += tripleDim_; // t3
        if (include_t4_) ccamps_dim_ += quadDim_; // t4
        if (include_u0_) ccamps_dim_++; // u0
        if (include_u1_) ccamps_dim_ += singleDim_; // u1
        if (include_u2_) ccamps_dim_ += doubleDim_; // u2
        if (include_u3_) ccamps_dim_ += tripleDim_; // u3
        if (include_u4_) ccamps_dim_ += quadDim_; // u4

        /// print included amplitudes
        Printf("  Included amplitudes:\n");
        Printf(include_u1_ ? "    U0\n" : "");
        Printf("    %2s", "T1");
        Printf("    %2s\n", include_u1_ ? "U1" : "  ");
        Printf("    %2s", "T2");
        Printf("    %2s\n", include_u2_ ? "U2" : "  ");
        if (include_t3_ || include_u3_) {
            Printf("    %2s", include_t3_ ? "T3" : "  ");
            Printf("    %2s\n", include_u3_ ? "U3" : "  ");
        }
        if (include_t4_ || include_u4_) {
            Printf("    %2s", include_t4_ ? "T4" : "  ");
            Printf("    %2s\n", include_u4_ ? "U4" : "  ");
        }
        Printf("\n  Dimension of CC amplitudes: %d\n\n", ccamps_dim_);
        if (has_photon_){
            // calculate cavity volumes: lambda = (1/(4*pi)*V_cav)^(-1/2) -> V_cav = (4*pi/lambda^2)
            double lam_norm = 0;
            for (double lam : lambda_)
                lam_norm += lam * lam;

            double V_cav = 4.0 * M_PI / lam_norm;

            // convert from bohr^3 to nm^3 (https://physics.nist.gov/cgi-bin/cuu/Value?Abohrrada0)
            double bohr_to_nm = 0.052917721090380;
            V_cav *= bohr_to_nm * bohr_to_nm * bohr_to_nm;

            // print cavity parameters
            Printf("  V_cav (nm3) = %6.4lf\n", V_cav);
            Printf("  lambda      = x: %6.4lf y: %6.4lf z: %6.4lf\n", lambda_[0], lambda_[1], lambda_[2]);
            Printf("  hw          = x: %6.4lf y: %6.4lf z: %6.4lf\n",
                   cavity_frequency_[0], cavity_frequency_[1], cavity_frequency_[2]);

        }
    }

    void QED_CC::update_amplitudes() {

        CC_Cavity::update_amplitudes(); // call parent function for t1 and t2

        /// dt = -residual / eps

        double *eps = epsilon_;
        size_t o = o_;
        size_t v = v_;
        size_t oa = oa_;
        size_t va = va_;

        // t3
        if (include_t3_) {
            HelperD::forall(residuals_["t3_aaaaaa"],
                            [eps, o, oa, va](auto &tile, auto &x) {
                                double o_ep = eps[x[3]] + eps[x[4]] + eps[x[5]],
                                        v_ep = eps[x[0]+o] + eps[x[1]+o] + eps[x[2]+o];
                                tile[x] /= (o_ep - v_ep);
                            });
            HelperD::forall(residuals_["t3_aabaab"],
                            [eps, o, oa, va](auto &tile, auto &x) {
                                double o_ep = eps[x[3]] + eps[x[4]] + eps[x[5]+oa],
                                        v_ep = eps[x[0]+o] + eps[x[1]+o] + eps[x[2]+o+va];
                                tile[x] /= (o_ep - v_ep);
                            });
            HelperD::forall(residuals_["t3_abbabb"],
                            [eps, o, oa, va](auto &tile, auto &x) {
                                double o_ep = eps[x[3]] + eps[x[4]+oa] + eps[x[5]+oa],
                                        v_ep = eps[x[0]+o] + eps[x[1]+o+va] + eps[x[2]+o+va];
                                tile[x] /= (o_ep - v_ep);
                            });
            HelperD::forall(residuals_["t3_bbbbbb"],
                            [eps, o, oa, va](auto &tile, auto &x) {
                                double o_ep = eps[x[3]+oa] + eps[x[4]+oa] + eps[x[5]+oa],
                                        v_ep = eps[x[0]+o+va] + eps[x[1]+o+va] + eps[x[2]+o+va];
                                tile[x] /= (o_ep - v_ep);
                            });
        }

        // t4
        if (include_t4_) {
            HelperD::forall(residuals_["t4_aaaaaaaa"],
                            [eps, o, oa, va](auto &tile, auto &x) {
                                double o_ep = eps[x[4]] + eps[x[5]] + eps[x[6]] + eps[x[7]],
                                        v_ep = eps[x[0]+o] + eps[x[1]+o] + eps[x[2]+o] + eps[x[3]+o];
                                tile[x] /= (o_ep - v_ep);
                            });
            HelperD::forall(residuals_["t4_aaabaaab"],
                            [eps, o, oa, va](auto &tile, auto &x) {
                                double o_ep = eps[x[4]] + eps[x[5]] + eps[x[6]] + eps[x[7]+oa],
                                        v_ep = eps[x[0]+o] + eps[x[1]+o] + eps[x[2]+o] + eps[x[3]+o+va];
                                tile[x] /= (o_ep - v_ep);
                            });
            HelperD::forall(residuals_["t4_aabbaabb"],
                            [eps, o, oa, va](auto &tile, auto &x) {
                                double o_ep = eps[x[4]] + eps[x[5]] + eps[x[6]+oa] + eps[x[7]+oa],
                                        v_ep = eps[x[0]+o] + eps[x[1]+o] + eps[x[2]+o+va] + eps[x[3]+o+va];
                                tile[x] /= (o_ep - v_ep);
                            });
            HelperD::forall(residuals_["t4_abbbabbb"],
                            [eps, o, oa, va](auto &tile, auto &x) {
                                double o_ep = eps[x[4]] + eps[x[5]+oa] + eps[x[6]+oa] + eps[x[7]+oa],
                                        v_ep = eps[x[0]+o] + eps[x[1]+o+va] + eps[x[2]+o+va] + eps[x[3]+o+va];
                                tile[x] /= (o_ep - v_ep);
                            });
            HelperD::forall(residuals_["t4_bbbbbbbb"],
                            [eps, o, oa, va](auto &tile, auto &x) {
                                double o_ep = eps[x[4]+oa] + eps[x[5]+oa] + eps[x[6]+oa] + eps[x[7]+oa],
                                        v_ep = eps[x[0]+o+va] + eps[x[1]+o+va] + eps[x[2]+o+va] + eps[x[3]+o+va];
                                tile[x] /= (o_ep - v_ep);
                            });
        }

        /// du = -residual / (eps + w)
        double w0 = cavity_frequency_[2];

        // u1
        if (include_u1_) {
            HelperD::forall(residuals_["u1_aa"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                tile[x] /= (eps[x[1]] - eps[x[0] + o]);
                            });

            HelperD::forall(residuals_["u1_bb"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                tile[x] /= (eps[x[1] + oa] - eps[x[0] + o + va]);
                            });
        }

        // u2
        if (include_u2_) {
            // u2
            HelperD::forall(residuals_["u2_aaaa"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                double o_ep = eps[x[2]] + eps[x[3]],
                                        v_ep = eps[x[0] + o] + eps[x[1] + o] + w0;
                                tile[x] /= (o_ep - v_ep);
                            });
            HelperD::forall(residuals_["u2_bbbb"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                double o_ep = eps[x[2] + oa] + eps[x[3] + oa],
                                        v_ep = eps[x[0] + o + va] + eps[x[1] + o + va] + w0;
                                tile[x] /= (o_ep - v_ep);
                            });
            HelperD::forall(residuals_["u2_abab"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                double o_ep = eps[x[2]] + eps[x[3] + oa],
                                        v_ep = eps[x[0] + o] + eps[x[1] + o] + w0;
                                tile[x] /= (o_ep - v_ep);
                            });
        }

        // u3
        if (include_u3_) {
            HelperD::forall(residuals_["u3_aaaaaa"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                double o_ep = eps[x[3]] + eps[x[4]] + eps[x[5]],
                                        v_ep = eps[x[0]+o] + eps[x[1]+o] + eps[x[2]+o] + w0;
                                tile[x] /= (o_ep - v_ep);
                            });
            HelperD::forall(residuals_["u3_aabaab"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                double o_ep = eps[x[3]] + eps[x[4]] + eps[x[5]+oa],
                                        v_ep = eps[x[0]+o] + eps[x[1]+o] + eps[x[2]+o+va] + w0;
                                tile[x] /= (o_ep - v_ep);
                            });
            HelperD::forall(residuals_["u3_abbabb"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                double o_ep = eps[x[3]] + eps[x[4]+oa] + eps[x[5]+oa],
                                        v_ep = eps[x[0]+o] + eps[x[1]+o+va] + eps[x[2]+o+va] + w0;
                                tile[x] /= (o_ep - v_ep);
                            });
            HelperD::forall(residuals_["u3_bbbbbb"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                double o_ep = eps[x[3]+oa] + eps[x[4]+oa] + eps[x[5]+oa],
                                        v_ep = eps[x[0]+o+va] + eps[x[1]+o+va] + eps[x[2]+o+va] + w0;
                                tile[x] /= (o_ep - v_ep);
                            });
        }

        // u4
        if (include_u4_) {
            HelperD::forall(residuals_["u4_aaaaaaaa"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                double o_ep = eps[x[4]] + eps[x[5]] + eps[x[6]] + eps[x[7]],
                                        v_ep = eps[x[0]+o] + eps[x[1]+o] + eps[x[2]+o] + eps[x[3]+o] + w0;
                                tile[x] /= (o_ep - v_ep);
                            });
            HelperD::forall(residuals_["u4_aaabaaab"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                double o_ep = eps[x[4]] + eps[x[5]] + eps[x[6]] + eps[x[7]+oa],
                                        v_ep = eps[x[0]+o] + eps[x[1]+o] + eps[x[2]+o] + eps[x[3]+o+va] + w0;
                                tile[x] /= (o_ep - v_ep);
                            });
            HelperD::forall(residuals_["u4_aabbaabb"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                double o_ep = eps[x[4]] + eps[x[5]] + eps[x[6]+oa] + eps[x[7]+oa],
                                        v_ep = eps[x[0]+o] + eps[x[1]+o] + eps[x[2]+o+va] + eps[x[3]+o+va] + w0;
                                tile[x] /= (o_ep - v_ep);
                            });
            HelperD::forall(residuals_["u4_abbbabbb"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                double o_ep = eps[x[4]] + eps[x[5]+oa] + eps[x[6]+oa] + eps[x[7]+oa],
                                        v_ep = eps[x[0]+o] + eps[x[1]+o+va] + eps[x[2]+o+va] + eps[x[3]+o+va] + w0;
                                tile[x] /= (o_ep - v_ep);
                            });
            HelperD::forall(residuals_["u4_bbbbbbbb"],
                            [eps, o, oa, va, w0](auto &tile, auto &x) {
                                double o_ep = eps[x[4]+oa] + eps[x[5]+oa] + eps[x[6]+oa] + eps[x[7]+oa],
                                        v_ep = eps[x[0]+o+va] + eps[x[1]+o+va] + eps[x[2]+o+va] + eps[x[3]+o+va] + w0;
                                tile[x] /= (o_ep - v_ep);
                            });
        }

        world_.gop.fence();
        double ru0_ = scalar_resids_["u0"];
        double u0_ = scalar_amps_["u0"];

        // u0
        if ( include_u0_ ) {
            ru0_ /= -w0;
            amplitudes_["u0"] = HelperD::makeTensor(world_, {1}, &u0_);
             residuals_["u0"] = HelperD::makeTensor(world_, {1}, &ru0_);
        }

        /// update amplitudes according to t + dt = amplitude - residual / eps
        if (include_t3_) {
            amplitudes_["t3_aaaaaa"](idx_map_[6]) += residuals_["t3_aaaaaa"](idx_map_[6]);
            amplitudes_["t3_aabaab"](idx_map_[6]) += residuals_["t3_aabaab"](idx_map_[6]);
            amplitudes_["t3_abbabb"](idx_map_[6]) += residuals_["t3_abbabb"](idx_map_[6]);
            amplitudes_["t3_bbbbbb"](idx_map_[6]) += residuals_["t3_bbbbbb"](idx_map_[6]);
        }
        if (include_t4_) {
            amplitudes_["t4_aaaaaaaa"](idx_map_[8]) += residuals_["t4_aaaaaaaa"](idx_map_[8]);
            amplitudes_["t4_aaabaaab"](idx_map_[8]) += residuals_["t4_aaabaaab"](idx_map_[8]);
            amplitudes_["t4_aabbaabb"](idx_map_[8]) += residuals_["t4_aabbaabb"](idx_map_[8]);
            amplitudes_["t4_abbbabbb"](idx_map_[8]) += residuals_["t4_abbbabbb"](idx_map_[8]);
            amplitudes_["t4_bbbbbbbb"](idx_map_[8]) += residuals_["t4_bbbbbbbb"](idx_map_[8]);
        }

        /// update amplitudes according to u + du = amplitude - residual / (eps + w)
        if (include_u0_) amplitudes_["u0"](idx_map_[1]) += residuals_["u0"](idx_map_[1]);
        if (include_u1_) {
            amplitudes_["u1_aa"](idx_map_[2]) += residuals_["u1_aa"](idx_map_[2]);
            amplitudes_["u1_bb"](idx_map_[2]) += residuals_["u1_bb"](idx_map_[2]);
        }
        if (include_u2_) {
            amplitudes_["u2_aaaa"](idx_map_[4]) += residuals_["u2_aaaa"](idx_map_[4]);
            amplitudes_["u2_abab"](idx_map_[4]) += residuals_["u2_abab"](idx_map_[4]);
            amplitudes_["u2_bbbb"](idx_map_[4]) += residuals_["u2_bbbb"](idx_map_[4]);
        }
        if (include_u3_) {
            amplitudes_["u3_aaaaaa"](idx_map_[6]) += residuals_["u3_aaaaaa"](idx_map_[6]);
            amplitudes_["u3_aabaab"](idx_map_[6]) += residuals_["u3_aabaab"](idx_map_[6]);
            amplitudes_["u3_abbabb"](idx_map_[6]) += residuals_["u3_abbabb"](idx_map_[6]);
            amplitudes_["u3_bbbbbb"](idx_map_[6]) += residuals_["u3_bbbbbb"](idx_map_[6]);
        }
        if (include_u4_) {
            amplitudes_["u4_aaaaaaaa"](idx_map_[8]) += residuals_["u4_aaaaaaaa"](idx_map_[8]);
            amplitudes_["u4_aaabaaab"](idx_map_[8]) += residuals_["u4_aaabaaab"](idx_map_[8]);
            amplitudes_["u4_aabbaabb"](idx_map_[8]) += residuals_["u4_aabbaabb"](idx_map_[8]);
            amplitudes_["u4_abbbabbb"](idx_map_[8]) += residuals_["u4_abbbabbb"](idx_map_[8]);
            amplitudes_["u4_bbbbbbbb"](idx_map_[8]) += residuals_["u4_bbbbbbbb"](idx_map_[8]);
        }

        world_.gop.fence();

    }

    double QED_CC::compute_residual_norms() {

        if(include_u0_) {
            double& u0 = scalar_amps_["u0"];
            world_.gop.fence();
            HelperD::forall(amplitudes_["u0"], [&u0](auto &tile, auto &x){
               u0 = tile[x];
            });
            world_.gop.fence();
        }


        /// residual norms for each amplitude

        resid_norms_["t1"] = sqrt(squared_norm(residuals_["t1_aa"])
                                  + squared_norm(residuals_["t1_bb"]));
        resid_norms_["t2"] = sqrt(squared_norm(residuals_["t2_aaaa"])
                                  + squared_norm(residuals_["t2_abab"])
                                  + squared_norm(residuals_["t2_bbbb"]));
        if (include_t3_) {
            resid_norms_["t3"] = sqrt(squared_norm(residuals_["t3_aaaaaa"])
                                      + squared_norm(residuals_["t3_aabaab"])
                                      + squared_norm(residuals_["t3_abbabb"])
                                      + squared_norm(residuals_["t3_bbbbbb"]));
        }
        if (include_t4_) {
            resid_norms_["t4"] = sqrt(squared_norm(residuals_["t4_aaaaaaaa"])
                                      + squared_norm(residuals_["t4_aaabaaab"])
                                      + squared_norm(residuals_["t4_aabbaabb"])
                                      + squared_norm(residuals_["t4_abbbabbb"])
                                      + squared_norm(residuals_["t4_bbbbbbbb"]));
        }

        // u residual norms
        if (include_u0_) {
            double ru0 = scalar_resids_["u0"];
            resid_norms_["u0"] = sqrt(ru0 * ru0);
        }
        if (include_u1_) {
            resid_norms_["u1"] = sqrt(squared_norm(residuals_["u1_aa"])
                                      + squared_norm(residuals_["u1_bb"]));
        }
        if (include_u2_) {
            resid_norms_["u2"] = sqrt(squared_norm(residuals_["u2_aaaa"])
                                      + squared_norm(residuals_["u2_abab"])
                                      + squared_norm(residuals_["u2_bbbb"]));
        }
        if (include_u3_) {
            resid_norms_["u3"] = sqrt(squared_norm(residuals_["u3_aaaaaa"])
                                      + squared_norm(residuals_["u3_aabaab"])
                                      + squared_norm(residuals_["u3_abbabb"])
                                      + squared_norm(residuals_["u3_bbbbbb"]));
        }
        if (include_u4_) {
            resid_norms_["u4"] = sqrt(squared_norm(residuals_["u4_aaaaaaaa"])
                                      + squared_norm(residuals_["u4_aaabaaab"])
                                      + squared_norm(residuals_["u4_aabbaabb"])
                                      + squared_norm(residuals_["u4_abbbabbb"])
                                      + squared_norm(residuals_["u4_bbbbbbbb"]));
        }


        // total residual norm
        double norm = 0.0;
        for (auto & amp : residuals_) {
            if (!amp.second.is_initialized()) continue;
            norm += squared_norm(amp.second);
        }

        return sqrt(norm);
    }
    
    void QED_CC::print_properties() {
        // calculate norms of amplitudes
        map<string, double> amp_norms;
        double total_norm = 0.0;
        for (auto& amp : amplitudes_) {
            if(!amp.second.is_initialized()) continue;
            // name of amplitude is first two characters of key
            double amp_norm = squared_norm(amp.second);
            string amp_name = amp.first.substr(0,2);
            amp_norms[amp_name] += amp_norm;
            total_norm += amp_norm;
        }

        double nT1 = amp_norms["t1"];
        double nT2 = amp_norms["t2"];
        double nT3 = amp_norms["t3"];
        double nT4 = amp_norms["t4"];
        double nU1 = amp_norms["u1"];
        double nU2 = amp_norms["u2"];
        double nU3 = amp_norms["u3"];
        double nU4 = amp_norms["u4"];

        double u0_ = scalar_amps_["u0"];

        // print norms of cluster amplitudes
        outfile->Printf("\n\n   Norms of cluster amplitudes:");
        outfile->Printf("\n   ----------------------------");
        outfile->Printf("\n    T1: %15.12lf | %5.2f %%"
                        "\n    T2: %15.12lf | %5.2f %%",
                        sqrt(nT1), 100*nT1/total_norm,
                        sqrt(nT2), 100*nT2/total_norm);
        if ( include_t3_ ) outfile->Printf("\n    T3: %15.12lf | %5.2f %%", sqrt(nT3), 100.0*nT3/total_norm);
        if ( include_u0_ ) outfile->Printf("\n    U0: %15.12lf | %5.2f %%", u0_, 100.0*fabs(u0_*u0_)/total_norm);
        if ( include_u1_ ) outfile->Printf("\n    U1: %15.12lf | %5.2f %%", sqrt(nU1), 100.0*nU1/total_norm);
        if ( include_u2_ ) outfile->Printf("\n    U2: %15.12lf | %5.2f %%", sqrt(nU2), 100.0*nU2/total_norm);
        if ( include_u3_ ) outfile->Printf("\n    U3: %15.12lf | %5.2f %%", sqrt(nU3), 100.0*nU3/total_norm);
        if ( include_u4_ ) outfile->Printf("\n    U4: %15.12lf | %5.2f %%", sqrt(nU4), 100.0*nU4/total_norm);
        outfile->Printf("\n   ----------------------------");
        outfile->Printf("\n    Total: %15.12lf\n\n", sqrt(total_norm));

    }

    void QED_CC::print_iter_header() const {
        Printf("\n");
        Printf("    ==>  Begin %s iterations <==    \n", cc_type_.c_str());
        Printf("\n");
        Printf("%5s %16s %15s %15s | %8s %8s",  "Iter","energy","dE","|dT|","|dT1|","|dT2|");
        if (include_t3_) Printf(" %8s","|dT3|");
        if (include_u0_) Printf(" %8s","|dU0|");
        if (include_u1_) Printf(" %8s","|dU1|");
        if (include_u2_) Printf(" %8s","|dU2|");
        if (include_u3_) Printf(" %8s","|dU3|");
        if (include_u4_) Printf(" %8s","|dU4|");
        Printf("\n");
    }

    void QED_CC::print_iteration(size_t iter, double energy, double dele, double tnorm) const {
        Printf("%5i %17.12lf %15.12lf %15.12lf | %-8.1e %-8.1e",iter,energy,dele,tnorm,resid_norms_.at("t1"),resid_norms_.at("t2"));
        if (include_t3_) Printf(" %-8.1e", resid_norms_.at("t3"));
        if (include_u0_) Printf(" %-8.1e", resid_norms_.at("u0"));
        if (include_u1_) Printf(" %-8.1e", resid_norms_.at("u1"));
        if (include_u2_) Printf(" %-8.1e", resid_norms_.at("u2"));
        if (include_u3_) Printf(" %-8.1e", resid_norms_.at("u3"));
        if (include_u4_) Printf(" %-8.1e", resid_norms_.at("u4"));
        Printf("\n");
    }
    
    void QED_CC::unpack_eris(TArrayMap &Qmo_blks) {



        // TODO: We should use Qmo_blks directly like the code commented out below
        //       but this is not working for some reason. The rest of the code is a messy workaround, but it works.
        /*for (auto& Ql : Qmo_blks) {
            for (auto& Qr : Qmo_blks) {

                string l_ov = Ql.first.substr(3, 2);
                string l_spin = Ql.first.substr(0, 2);
                string r_ov = Qr.first.substr(3, 2);
                string r_spin = Qr.first.substr(0, 2);

                // if exchange is not possible, skip
                if (l_ov[1] != r_ov[1]) continue;
                if (l_spin[1] != r_spin[1]) continue;

                // spin must not be odd. check for odd 'a' or 'b' spin
                size_t num_a = 0, num_b = 0; // too lazy to write a function
                if (l_spin[0] == 'a') num_a++; if (l_spin[0] == 'b') num_b++;
                if (l_spin[1] == 'a') num_a++; if (l_spin[1] == 'b') num_b++;
                if (r_spin[0] == 'a') num_a++; if (r_spin[0] == 'b') num_b++;
                if (r_spin[1] == 'a') num_a++; if (r_spin[1] == 'b') num_b++;
                if (num_a % 2 != 0 || num_b % 2 != 0) continue;


                // build 4-index container

                TArrayD pqrs;
                TArrayD psrq;
                pqrs("p,q,r,s") = Ql.second("Q,p,q") * Qr.second("Q,r,s");
                psrq("p,s,r,q") = Ql.second("Q,p,s") * Qr.second("Q,r,q");

                // add dipole integrals
                if (l_spin[0] == l_spin[1] && r_spin[0] == r_spin[1] &&
                    l_ov[0] == l_ov[1] && r_ov[0] == r_ov[1]) {

                    string l_dip = Ql.first;
                    string r_dip = Qr.first;
                    pqrs("p,q,r,s") += lambda_[0]*lambda_[0] * Dip_blks_["dx_"+l_dip]("p,q") * Dip_blks_["dx_"+r_dip]("r,s");
                    pqrs("p,q,r,s") += lambda_[1]*lambda_[1] * Dip_blks_["dy_"+l_dip]("p,q") * Dip_blks_["dy_"+r_dip]("r,s");
                    pqrs("p,q,r,s") += lambda_[2]*lambda_[2] * Dip_blks_["dz_"+l_dip]("p,q") * Dip_blks_["dz_"+r_dip]("r,s");
                }
                if (l_spin[0] == r_spin[1] && l_spin[1] == r_spin[0] &&
                    l_ov[0] == r_ov[1] && l_ov[1] == r_ov[0]) {

                    char l_dip[5] = {l_spin[0], r_spin[1], '_', l_ov[0], r_ov[1]};
                    char r_dip[5] = {l_spin[1], r_spin[0], '_', l_ov[1], r_ov[0]};
                    string l_str(l_dip);
                    string r_str(r_dip);

                    // l_str is supposed to be something like aa_ov, but for some reason gdb thinks it is aa_ov\026aa_vo
                    // how in the world is that possible? (break qed_cc_spin.cc:40
                    // anyway, this is a hack to fix it
                    if (l_str.size() > 5) l_str = l_str.substr(0, 5);
                    if (r_str.size() > 5) r_str = r_str.substr(0, 5);

                    psrq("p,s,r,q") += lambda_[0]*lambda_[0] * Dip_blks_["dx_"+l_str]("p,s") * Dip_blks_["dx_"+r_str]("r,q");
                    psrq("p,s,r,q") += lambda_[1]*lambda_[1] * Dip_blks_["dy_"+l_str]("p,s") * Dip_blks_["dy_"+r_str]("r,q");
                    psrq("p,s,r,q") += lambda_[2]*lambda_[2] * Dip_blks_["dz_"+l_str]("p,s") * Dip_blks_["dz_"+r_str]("r,q");
                }

                // build electron repulsion integrals
                char v_blk[9] = {l_spin[0], r_spin[1], l_spin[0], r_spin[1], '_', l_ov[0], r_ov[0], l_ov[1], r_ov[1]};
                string v_blk_str(v_blk);
                if (v_blk_str.size() > 9) v_blk_str = v_blk_str.substr(0, 9);
                V_blks_[v_blk_str]("p,s,q,r") = pqrs("p,q,r,s") - psrq("p,s,r,q");
//                Printf(Ql.first.c_str());
//                Printf("\n");
//                Printf(Qr.first.c_str());
//                Printf("\n");
//                Printf(v_blk_str.c_str());
//                Printf("\n\n");
            }
        }
//        exit(0);*/



        /*******************************************************************************************
         * Compute non-antisymmetric 4-electron integrals
         *******************************************************************************************/
        /// (note: not all of these will be used. Some containers can be reused, and some contractions are not allowed)

        double lambda_x = lambda_[0];
        double lambda_y = lambda_[1];
        double lambda_z = lambda_[2];

        double lx2 = lambda_x * lambda_x;
        double ly2 = lambda_y * lambda_y;
        double lz2 = lambda_z * lambda_z;

        // aaaa
        TArrayD ikjl_aaaa, adbc_aaaa, aibj_aaaa, ibaj_aaaa, jika_aaaa, acib_aaaa,
                iljk_aaaa, iajb_aaaa, ajbi_aaaa, ijak_aaaa, jaki_aaaa, aibc_aaaa,
                acbd_aaaa, ibja_aaaa, ijab_aaaa, ikaj_aaaa, abic_aaaa, acbi_aaaa;
        // aaab
        TArrayD ikjl_aaab, adbc_aaab, aibj_aaab, ibaj_aaab, jika_aaab, acib_aaab,
                iljk_aaab, iajb_aaab, ajbi_aaab, ijak_aaab, jaki_aaab, aibc_aaab,
                acbd_aaab, ibja_aaab, ijab_aaab, ikaj_aaab, abic_aaab, acbi_aaab;
        // aaba
        TArrayD ikjl_aaba, adbc_aaba, aibj_aaba, ibaj_aaba, jika_aaba, acib_aaba,
                iljk_aaba, iajb_aaba, ajbi_aaba, ijak_aaba, jaki_aaba, aibc_aaba,
                acbd_aaba, ibja_aaba, ijab_aaba, ikaj_aaba, abic_aaba, acbi_aaba;
        // aabb
        TArrayD ikjl_aabb, adbc_aabb, aibj_aabb, ibaj_aabb, jika_aabb, acib_aabb,
                iljk_aabb, iajb_aabb, ajbi_aabb, ijak_aabb, jaki_aabb, aibc_aabb,
                acbd_aabb, ibja_aabb, ijab_aabb, ikaj_aabb, abic_aabb, acbi_aabb;
        // abaa
        TArrayD ikjl_abaa, adbc_abaa, aibj_abaa, ibaj_abaa, jika_abaa, acib_abaa,
                iljk_abaa, iajb_abaa, ajbi_abaa, ijak_abaa, jaki_abaa, aibc_abaa,
                acbd_abaa, ibja_abaa, ijab_abaa, ikaj_abaa, abic_abaa, acbi_abaa;
        // abab
        TArrayD ikjl_abab, adbc_abab, aibj_abab, ibaj_abab, jika_abab, acib_abab,
                iljk_abab, iajb_abab, ajbi_abab, ijak_abab, jaki_abab, aibc_abab,
                acbd_abab, ibja_abab, ijab_abab, ikaj_abab, abic_abab, acbi_abab;
        // abba
        TArrayD ikjl_abba, adbc_abba, aibj_abba, ibaj_abba, jika_abba, acib_abba,
                iljk_abba, iajb_abba, ajbi_abba, ijak_abba, jaki_abba, aibc_abba,
                acbd_abba, ibja_abba, ijab_abba, ikaj_abba, abic_abba, acbi_abba;
        // abbb
        TArrayD ikjl_abbb, adbc_abbb, aibj_abbb, ibaj_abbb, jika_abbb, acib_abbb,
                iljk_abbb, iajb_abbb, ajbi_abbb, ijak_abbb, jaki_abbb, aibc_abbb,
                acbd_abbb, ibja_abbb, ijab_abbb, ikaj_abbb, abic_abbb, acbi_abbb;
        // baaa
        TArrayD ikjl_baaa, adbc_baaa, aibj_baaa, ibaj_baaa, jika_baaa, acib_baaa,
                iljk_baaa, iajb_baaa, ajbi_baaa, ijak_baaa, jaki_baaa, aibc_baaa,
                acbd_baaa, ibja_baaa, ijab_baaa, ikaj_baaa, abic_baaa, acbi_baaa;
        // baab
        TArrayD ikjl_baab, adbc_baab, aibj_baab, ibaj_baab, jika_baab, acib_baab,
                iljk_baab, iajb_baab, ajbi_baab, ijak_baab, jaki_baab, aibc_baab,
                acbd_baab, ibja_baab, ijab_baab, ikaj_baab, abic_baab, acbi_baab;
        // baba
        TArrayD ikjl_baba, adbc_baba, aibj_baba, ibaj_baba, jika_baba, acib_baba,
                iljk_baba, iajb_baba, ajbi_baba, ijak_baba, jaki_baba, aibc_baba,
                acbd_baba, ibja_baba, ijab_baba, ikaj_baba, abic_baba, acbi_baba;
        // babb
        TArrayD ikjl_babb, adbc_babb, aibj_babb, ibaj_babb, jika_babb, acib_babb,
                iljk_babb, iajb_babb, ajbi_babb, ijak_babb, jaki_babb, aibc_babb,
                acbd_babb, ibja_babb, ijab_babb, ikaj_babb, abic_babb, acbi_babb;
        // bbaa
        TArrayD ikjl_bbaa, adbc_bbaa, aibj_bbaa, ibaj_bbaa, jika_bbaa, acib_bbaa,
                iljk_bbaa, iajb_bbaa, ajbi_bbaa, ijak_bbaa, jaki_bbaa, aibc_bbaa,
                acbd_bbaa, ibja_bbaa, ijab_bbaa, ikaj_bbaa, abic_bbaa, acbi_bbaa;
        // bbab
        TArrayD ikjl_bbab, adbc_bbab, aibj_bbab, ibaj_bbab, jika_bbab, acib_bbab,
                iljk_bbab, iajb_bbab, ajbi_bbab, ijak_bbab, jaki_bbab, aibc_bbab,
                acbd_bbab, ibja_bbab, ijab_bbab, ikaj_bbab, abic_bbab, acbi_bbab;
        // bbba
        TArrayD ikjl_bbba, adbc_bbba, aibj_bbba, ibaj_bbba, jika_bbba, acib_bbba,
                iljk_bbba, iajb_bbba, ajbi_bbba, ijak_bbba, jaki_bbba, aibc_bbba,
                acbd_bbba, ibja_bbba, ijab_bbba, ikaj_bbba, abic_bbba, acbi_bbba;
        // bbbb
        TArrayD ikjl_bbbb, adbc_bbbb, aibj_bbbb, ibaj_bbbb, jika_bbbb, acib_bbbb,
                iljk_bbbb, iajb_bbbb, ajbi_bbbb, ijak_bbbb, jaki_bbbb, aibc_bbbb,
                acbd_bbbb, ibja_bbbb, ijab_bbbb, ikaj_bbbb, abic_bbbb, acbi_bbbb;

        /*******************************************************************************************
         * <ij|kl> = (ik|jl) - (il|jk)
         *******************************************************************************************/

        bool pure_photon_integrals_ = false; // TODO: remove this. Not needed anymore.

        /// aaaa - aaaa
        if (!pure_photon_integrals_) ikjl_aaaa("i, k, j, l") = Qmo_blks["aa_oo"]("Q, i, k") * Qmo_blks["aa_oo"]("Q, j, l");
        else ikjl_aaaa = HelperD::makeTensor(world_, {oa_, oa_, oa_, oa_}, true);

        // add dipole terms
        ikjl_aaaa("i, k, j, l") += lx2 * Dip_blks_["dx_aa_oo"]("i, k") * Dip_blks_["dx_aa_oo"]("j, l");
        ikjl_aaaa("i, k, j, l") += ly2 * Dip_blks_["dy_aa_oo"]("i, k") * Dip_blks_["dy_aa_oo"]("j, l");
        ikjl_aaaa("i, k, j, l") += lz2 * Dip_blks_["dz_aa_oo"]("i, k") * Dip_blks_["dz_aa_oo"]("j, l");

        // make antisymmetrized
        V_blks_["aaaa_oooo"]("i, j, k, l") = ikjl_aaaa("i, k, j, l") - ikjl_aaaa("i, l, j, k");

        /// aabb - abba
        if (!pure_photon_integrals_) {
            ikjl_aabb("i, k, j, l") = Qmo_blks["aa_oo"]("Q, i, k") * Qmo_blks["bb_oo"]("Q, j, l");
            iljk_abba("i, l, j, k") = Qmo_blks["ab_oo"]("Q, i, l") * Qmo_blks["ba_oo"]("Q, j, k");
        } else {
            ikjl_aabb = HelperD::makeTensor(world_, {oa_, oa_, ob_, ob_}, true);
            iljk_abba = HelperD::makeTensor(world_, {oa_, ob_, ob_, oa_}, true);
        }

        // add dipole terms
        ikjl_aabb("i, k, j, l") += lx2 * Dip_blks_["dx_aa_oo"]("i, k") * Dip_blks_["dx_bb_oo"]("j, l");
        ikjl_aabb("i, k, j, l") += ly2 * Dip_blks_["dy_aa_oo"]("i, k") * Dip_blks_["dy_bb_oo"]("j, l");
        ikjl_aabb("i, k, j, l") += lz2 * Dip_blks_["dz_aa_oo"]("i, k") * Dip_blks_["dz_bb_oo"]("j, l");

        // make antisymmetrized
        V_blks_["abab_oooo"]("i, j, k, l") = ikjl_aabb("i, k, j, l") - iljk_abba("i, l, j, k");

        /// bbbb - bbbb
        if (!pure_photon_integrals_) ikjl_bbbb("i, k, j, l") = Qmo_blks["bb_oo"]("Q, i, k") * Qmo_blks["bb_oo"]("Q, j, l");
        else ikjl_bbbb = HelperD::makeTensor(world_, {ob_, ob_, ob_, ob_}, true);

        // add dipole terms
        ikjl_bbbb("i, k, j, l") += lx2 * Dip_blks_["dx_bb_oo"]("i, k") * Dip_blks_["dx_bb_oo"]("j, l");
        ikjl_bbbb("i, k, j, l") += ly2 * Dip_blks_["dy_bb_oo"]("i, k") * Dip_blks_["dy_bb_oo"]("j, l");
        ikjl_bbbb("i, k, j, l") += lz2 * Dip_blks_["dz_bb_oo"]("i, k") * Dip_blks_["dz_bb_oo"]("j, l");

        // make antisymmetrized
        V_blks_["bbbb_oooo"]("i, j, k, l") = ikjl_bbbb("i, k, j, l") - ikjl_bbbb("i, l, j, k");

        /*******************************************************************************************
         * <ab||cd> = (ac|bd) - (ad|bc)
         *******************************************************************************************/

        /// aaaa - aaaa
        if (!pure_photon_integrals_) acbd_aaaa("a, c, b, d") = Qmo_blks["aa_vv"]("Q, a, c") * Qmo_blks["aa_vv"]("Q, b, d");
        else acbd_aaaa = HelperD::makeTensor(world_, {va_, va_, va_, va_}, true);

        // add dipole terms
        acbd_aaaa("a, c, b, d") += lx2 * Dip_blks_["dx_aa_vv"]("a, c") * Dip_blks_["dx_aa_vv"]("b, d");
        acbd_aaaa("a, c, b, d") += ly2 * Dip_blks_["dy_aa_vv"]("a, c") * Dip_blks_["dy_aa_vv"]("b, d");
        acbd_aaaa("a, c, b, d") += lz2 * Dip_blks_["dz_aa_vv"]("a, c") * Dip_blks_["dz_aa_vv"]("b, d");

        // make antisymmetrized
        V_blks_["aaaa_vvvv"]("a, b, c, d") = acbd_aaaa("a, c, b, d") - acbd_aaaa("a, d, b, c");

        /// aabb - abba
        if (!pure_photon_integrals_) {
            acbd_aabb("a, c, b, d") = Qmo_blks["aa_vv"]("Q, a, c") * Qmo_blks["bb_vv"]("Q, b, d");
            adbc_abba("a, d, b, c") = Qmo_blks["ab_vv"]("Q, a, d") * Qmo_blks["ba_vv"]("Q, b, c");
        } else {
            acbd_aabb = HelperD::makeTensor(world_, {va_, va_, vb_, vb_}, true);
            adbc_abba = HelperD::makeTensor(world_, {va_, vb_, vb_, va_}, true);
        }

        // add dipole terms
        acbd_aabb("a, c, b, d") += lx2 * Dip_blks_["dx_aa_vv"]("a, c") * Dip_blks_["dx_bb_vv"]("b, d");
        acbd_aabb("a, c, b, d") += ly2 * Dip_blks_["dy_aa_vv"]("a, c") * Dip_blks_["dy_bb_vv"]("b, d");
        acbd_aabb("a, c, b, d") += lz2 * Dip_blks_["dz_aa_vv"]("a, c") * Dip_blks_["dz_bb_vv"]("b, d");

        // make antisymmetrized
        V_blks_["abab_vvvv"]("a, b, c, d") = acbd_aabb("a, c, b, d") - adbc_abba("a, d, b, c");

        /// bbbb - bbbb
        if (!pure_photon_integrals_) acbd_bbbb("a, c, b, d") = Qmo_blks["bb_vv"]("Q, a, c") * Qmo_blks["bb_vv"]("Q, b, d");
        else acbd_bbbb = HelperD::makeTensor(world_, {vb_, vb_, vb_, vb_}, true);

        // add dipole terms
        acbd_bbbb("a, c, b, d") += lx2 * Dip_blks_["dx_bb_vv"]("a, c") * Dip_blks_["dx_bb_vv"]("b, d");
        acbd_bbbb("a, c, b, d") += ly2 * Dip_blks_["dy_bb_vv"]("a, c") * Dip_blks_["dy_bb_vv"]("b, d");
        acbd_bbbb("a, c, b, d") += lz2 * Dip_blks_["dz_bb_vv"]("a, c") * Dip_blks_["dz_bb_vv"]("b, d");

        // make antisymmetrized
        V_blks_["bbbb_vvvv"]("a, b, c, d") = acbd_bbbb("a, c, b, d") - acbd_bbbb("a, d, b, c");

        /***************************************************************************
         * <ij||ab> = (ia|jb) - (ib|ja)
         **************************************************************************/

        /// aaaa - aaaa
        if (!pure_photon_integrals_) iajb_aaaa("i, a, j, b") = Qmo_blks["aa_ov"]("Q, i, a") * Qmo_blks["aa_ov"]("Q, j, b");
        else iajb_aaaa = HelperD::makeTensor(world_, {oa_, va_, oa_, va_}, true);

        // add dipole terms
        iajb_aaaa("i, a, j, b") += lx2 * Dip_blks_["dx_aa_ov"]("i, a") * Dip_blks_["dx_aa_ov"]("j, b");
        iajb_aaaa("i, a, j, b") += ly2 * Dip_blks_["dy_aa_ov"]("i, a") * Dip_blks_["dy_aa_ov"]("j, b");
        iajb_aaaa("i, a, j, b") += lz2 * Dip_blks_["dz_aa_ov"]("i, a") * Dip_blks_["dz_aa_ov"]("j, b");

        // make antisymmetrized
        V_blks_["aaaa_oovv"]("i, j, a, b") = iajb_aaaa("i, a, j, b") - iajb_aaaa("i, b, j, a");

        /// aabb - abba
        if (!pure_photon_integrals_) {
            iajb_aabb("i, a, j, b") = Qmo_blks["aa_ov"]("Q, i, a") * Qmo_blks["bb_ov"]("Q, j, b");
            ibja_abba("i, b, j, a") = Qmo_blks["ab_ov"]("Q, i, b") * Qmo_blks["ba_ov"]("Q, j, a");
        } else {
            iajb_aabb = HelperD::makeTensor(world_, {oa_, va_, ob_, vb_}, true);
            ibja_abba = HelperD::makeTensor(world_, {oa_, vb_, ob_, va_}, true);
        }

        // add dipole terms
        iajb_aabb("i, a, j, b") += lx2 * Dip_blks_["dx_aa_ov"]("i, a") * Dip_blks_["dx_bb_ov"]("j, b");
        iajb_aabb("i, a, j, b") += ly2 * Dip_blks_["dy_aa_ov"]("i, a") * Dip_blks_["dy_bb_ov"]("j, b");
        iajb_aabb("i, a, j, b") += lz2 * Dip_blks_["dz_aa_ov"]("i, a") * Dip_blks_["dz_bb_ov"]("j, b");

        // make antisymmetrized
        V_blks_["abab_oovv"]("i, j, a, b") = iajb_aabb("i, a, j, b") - ibja_abba("i, b, j, a");

        /// bbbb - bbbb
        if (!pure_photon_integrals_) iajb_bbbb("i, a, j, b") = Qmo_blks["bb_ov"]("Q, i, a") * Qmo_blks["bb_ov"]("Q, j, b");
        else iajb_bbbb = HelperD::makeTensor(world_, {ob_, vb_, ob_, vb_}, true);

        // add dipole terms
        iajb_bbbb("i, a, j, b") += lx2 * Dip_blks_["dx_bb_ov"]("i, a") * Dip_blks_["dx_bb_ov"]("j, b");
        iajb_bbbb("i, a, j, b") += ly2 * Dip_blks_["dy_bb_ov"]("i, a") * Dip_blks_["dy_bb_ov"]("j, b");
        iajb_bbbb("i, a, j, b") += lz2 * Dip_blks_["dz_bb_ov"]("i, a") * Dip_blks_["dz_bb_ov"]("j, b");

        // make antisymmetrized
        V_blks_["bbbb_oovv"]("i, j, a, b") = iajb_bbbb("i, a, j, b") - iajb_bbbb("i, b, j, a");

        /***************************************************************************
         * <ab||ij> = (ai|bj) - (aj|bi)
         **************************************************************************/

        /// aaaa - aaaa
        if (!pure_photon_integrals_) aibj_aaaa("a, i, b, j") = Qmo_blks["aa_vo"]("Q, a, i") * Qmo_blks["aa_vo"]("Q, b, j");
        else aibj_aaaa = HelperD::makeTensor(world_, {va_, oa_, va_, oa_}, true);

        // add dipole terms
        aibj_aaaa("a, i, b, j") += lx2 * Dip_blks_["dx_aa_vo"]("a, i") * Dip_blks_["dx_aa_vo"]("b, j");
        aibj_aaaa("a, i, b, j") += ly2 * Dip_blks_["dy_aa_vo"]("a, i") * Dip_blks_["dy_aa_vo"]("b, j");
        aibj_aaaa("a, i, b, j") += lz2 * Dip_blks_["dz_aa_vo"]("a, i") * Dip_blks_["dz_aa_vo"]("b, j");

        // make antisymmetrized
        V_blks_["aaaa_vvoo"]("a, b, i, j") = aibj_aaaa("a, i, b, j") - aibj_aaaa("a, j, b, i");

        /// aabb - abba
        if (!pure_photon_integrals_) {
            aibj_aabb("a, i, b, j") = Qmo_blks["aa_vo"]("Q, a, i") * Qmo_blks["bb_vo"]("Q, b, j");
            ajbi_abba("a, j, b, i") = Qmo_blks["ab_vo"]("Q, a, j") * Qmo_blks["ba_vo"]("Q, b, i");
        } else {
            aibj_aabb = HelperD::makeTensor(world_, {va_, oa_, vb_, ob_}, true);
            ajbi_abba = HelperD::makeTensor(world_, {va_, ob_, vb_, oa_}, true);
        }

        // add dipole terms
        aibj_aabb("a, i, b, j") += lx2 * Dip_blks_["dx_aa_vo"]("a, i") * Dip_blks_["dx_bb_vo"]("b, j");
        aibj_aabb("a, i, b, j") += ly2 * Dip_blks_["dy_aa_vo"]("a, i") * Dip_blks_["dy_bb_vo"]("b, j");
        aibj_aabb("a, i, b, j") += lz2 * Dip_blks_["dz_aa_vo"]("a, i") * Dip_blks_["dz_bb_vo"]("b, j");

        // make antisymmetrized
        V_blks_["abab_vvoo"]("a, b, i, j") = aibj_aabb("a, i, b, j") - ajbi_abba("a, j, b, i");

        /// bbbb - bbbb

        if (!pure_photon_integrals_) aibj_bbbb("a, i, b, j") = Qmo_blks["bb_vo"]("Q, a, i") * Qmo_blks["bb_vo"]("Q, b, j");
        else aibj_bbbb = HelperD::makeTensor(world_, {vb_, ob_, vb_, ob_}, true);

        // add dipole terms
        aibj_bbbb("a, i, b, j") += lx2 * Dip_blks_["dx_bb_vo"]("a, i") * Dip_blks_["dx_bb_vo"]("b, j");
        aibj_bbbb("a, i, b, j") += ly2 * Dip_blks_["dy_bb_vo"]("a, i") * Dip_blks_["dy_bb_vo"]("b, j");
        aibj_bbbb("a, i, b, j") += lz2 * Dip_blks_["dz_bb_vo"]("a, i") * Dip_blks_["dz_bb_vo"]("b, j");

        // make antisymmetrized
        V_blks_["bbbb_vvoo"]("a, b, i, j") = aibj_bbbb("a, i, b, j") - aibj_bbbb("a, j, b, i");

        /***************************************************************************
         * <ia||jb> = (ij|ab) - (ib|aj)
         **************************************************************************/

        /// aaaa - aaaa
        if (!pure_photon_integrals_) {
            ijab_aaaa("i, j, a, b") = Qmo_blks["aa_oo"]("Q, i, j") * Qmo_blks["aa_vv"]("Q, a, b");
            ibaj_aaaa("i, b, a, j") = Qmo_blks["aa_ov"]("Q, i, b") * Qmo_blks["aa_vo"]("Q, a, j");
        } else {
            ijab_aaaa = HelperD::makeTensor(world_, {oa_, oa_, va_, va_}, true);
            ibaj_aaaa = HelperD::makeTensor(world_, {oa_, va_, va_, oa_}, true);
        }

        // add dipole terms
        ijab_aaaa("i, j, a, b") += lx2 * Dip_blks_["dx_aa_oo"]("i, j") * Dip_blks_["dx_aa_vv"]("a, b");
        ijab_aaaa("i, j, a, b") += ly2 * Dip_blks_["dy_aa_oo"]("i, j") * Dip_blks_["dy_aa_vv"]("a, b");
        ijab_aaaa("i, j, a, b") += lz2 * Dip_blks_["dz_aa_oo"]("i, j") * Dip_blks_["dz_aa_vv"]("a, b");

        ibaj_aaaa("i, b, a, j") += lx2 * Dip_blks_["dx_aa_ov"]("i, b") * Dip_blks_["dx_aa_vo"]("a, j");
        ibaj_aaaa("i, b, a, j") += ly2 * Dip_blks_["dy_aa_ov"]("i, b") * Dip_blks_["dy_aa_vo"]("a, j");
        ibaj_aaaa("i, b, a, j") += lz2 * Dip_blks_["dz_aa_ov"]("i, b") * Dip_blks_["dz_aa_vo"]("a, j");

        // make antisymmetrized
        V_blks_["aaaa_ovov"]("i, a, j, b") = ijab_aaaa("i, j, a, b") - ibaj_aaaa("i, b, a, j");

        /// aabb - abba
        if (!pure_photon_integrals_) {
            ijab_aabb("i, j, a, b") = Qmo_blks["aa_oo"]("Q, i, j") * Qmo_blks["bb_vv"]("Q, a, b");
            ibaj_abba("i, b, a, j") = Qmo_blks["ab_ov"]("Q, i, b") * Qmo_blks["ba_vo"]("Q, a, j");
        } else {
            ijab_aabb = HelperD::makeTensor(world_, {oa_, oa_, vb_, vb_}, true);
            ibaj_abba = HelperD::makeTensor(world_, {oa_, vb_, vb_, oa_}, true);
        }

        // add dipole terms
        ijab_aabb("i, j, a, b") += lx2 * Dip_blks_["dx_aa_oo"]("i, j") * Dip_blks_["dx_bb_vv"]("a, b");
        ijab_aabb("i, j, a, b") += ly2 * Dip_blks_["dy_aa_oo"]("i, j") * Dip_blks_["dy_bb_vv"]("a, b");
        ijab_aabb("i, j, a, b") += lz2 * Dip_blks_["dz_aa_oo"]("i, j") * Dip_blks_["dz_bb_vv"]("a, b");

        // make antisymmetrized
        V_blks_["abab_ovov"]("i, a, j, b") = ijab_aabb("i, j, a, b") - ibaj_abba("i, b, a, j");

        /// abba - aabb
        if (!pure_photon_integrals_) {
            ijab_abba("i, j, a, b") = Qmo_blks["ab_oo"]("Q, i, j") * Qmo_blks["ba_vv"]("Q, a, b");
            ibaj_aabb("i, b, a, j") = Qmo_blks["aa_ov"]("Q, i, b") * Qmo_blks["bb_vo"]("Q, a, j");
        } else {
            ijab_abba = HelperD::makeTensor(world_, {oa_, ob_, vb_, va_}, true);
            ibaj_aabb = HelperD::makeTensor(world_, {oa_, va_, vb_, ob_}, true);
        }

        // add dipole terms
        ibaj_aabb("i, b, a, j") += lx2 * Dip_blks_["dx_aa_ov"]("i, b") * Dip_blks_["dx_bb_vo"]("a, j");
        ibaj_aabb("i, b, a, j") += ly2 * Dip_blks_["dy_aa_ov"]("i, b") * Dip_blks_["dy_bb_vo"]("a, j");
        ibaj_aabb("i, b, a, j") += lz2 * Dip_blks_["dz_aa_ov"]("i, b") * Dip_blks_["dz_bb_vo"]("a, j");

        // make antisymmetrized
        V_blks_["abba_ovov"]("i, a, j, b") = ijab_abba("i, j, a, b") - ibaj_aabb("i, b, a, j");

        /// baab - bbaa
        if (!pure_photon_integrals_) {
            ijab_baab("i, j, a, b") = Qmo_blks["ba_oo"]("Q, i, j") * Qmo_blks["ab_vv"]("Q, a, b");
            ibaj_bbaa("i, b, a, j") = Qmo_blks["bb_ov"]("Q, i, b") * Qmo_blks["aa_vo"]("Q, a, j");
        } else {
            ijab_baab = HelperD::makeTensor(world_, {ob_, oa_, va_, vb_}, true);
            ibaj_bbaa = HelperD::makeTensor(world_, {ob_, vb_, va_, oa_}, true);
        }

        // add dipole terms
        ibaj_bbaa("i, b, a, j") += lx2 * Dip_blks_["dx_bb_ov"]("i, b") * Dip_blks_["dx_aa_vo"]("a, j");
        ibaj_bbaa("i, b, a, j") += ly2 * Dip_blks_["dy_bb_ov"]("i, b") * Dip_blks_["dy_aa_vo"]("a, j");
        ibaj_bbaa("i, b, a, j") += lz2 * Dip_blks_["dz_bb_ov"]("i, b") * Dip_blks_["dz_aa_vo"]("a, j");

        // make antisymmetrized
        V_blks_["baab_ovov"]("i, a, j, b") = ijab_baab("i, j, a, b") - ibaj_bbaa("i, b, a, j");

        /// bbaa - baab
        if (!pure_photon_integrals_) {
            ijab_bbaa("i, j, a, b") = Qmo_blks["bb_oo"]("Q, i, j") * Qmo_blks["aa_vv"]("Q, a, b");
            ibaj_baab("i, b, a, j") = Qmo_blks["ba_ov"]("Q, i, b") * Qmo_blks["ab_vo"]("Q, a, j");
        } else {
            ijab_bbaa = HelperD::makeTensor(world_, {ob_, ob_, va_, va_}, true);
            ibaj_baab = HelperD::makeTensor(world_, {ob_, va_, va_, ob_}, true);
        }

        // add dipole terms
        ijab_bbaa("i, j, a, b") += lx2 * Dip_blks_["dx_bb_oo"]("i, j") * Dip_blks_["dx_aa_vv"]("a, b");
        ijab_bbaa("i, j, a, b") += ly2 * Dip_blks_["dy_bb_oo"]("i, j") * Dip_blks_["dy_aa_vv"]("a, b");
        ijab_bbaa("i, j, a, b") += lz2 * Dip_blks_["dz_bb_oo"]("i, j") * Dip_blks_["dz_aa_vv"]("a, b");

        // make antisymmetrized
        V_blks_["baba_ovov"]("i, a, j, b") = ijab_bbaa("i, j, a, b") - ibaj_baab("i, b, a, j");

        /// bbbb - bbbb
        if (!pure_photon_integrals_) {
            ijab_bbbb("i, j, a, b") = Qmo_blks["bb_oo"]("Q, i, j") * Qmo_blks["bb_vv"]("Q, a, b");
            ibaj_bbbb("i, b, a, j") = Qmo_blks["bb_ov"]("Q, i, b") * Qmo_blks["bb_vo"]("Q, a, j");
        } else {
            ijab_bbbb = HelperD::makeTensor(world_, {ob_, ob_, vb_, vb_}, true);
            ibaj_bbbb = HelperD::makeTensor(world_, {ob_, vb_, vb_, ob_}, true);
        }

        // add dipole terms
        ijab_bbbb("i, j, a, b") += lx2 * Dip_blks_["dx_bb_oo"]("i, j") * Dip_blks_["dx_bb_vv"]("a, b");
        ijab_bbbb("i, j, a, b") += ly2 * Dip_blks_["dy_bb_oo"]("i, j") * Dip_blks_["dy_bb_vv"]("a, b");
        ijab_bbbb("i, j, a, b") += lz2 * Dip_blks_["dz_bb_oo"]("i, j") * Dip_blks_["dz_bb_vv"]("a, b");

        ibaj_bbbb("i, b, a, j") += lx2 * Dip_blks_["dx_bb_ov"]("i, b") * Dip_blks_["dx_bb_vo"]("a, j");
        ibaj_bbbb("i, b, a, j") += ly2 * Dip_blks_["dy_bb_ov"]("i, b") * Dip_blks_["dy_bb_vo"]("a, j");
        ibaj_bbbb("i, b, a, j") += lz2 * Dip_blks_["dz_bb_ov"]("i, b") * Dip_blks_["dz_bb_vo"]("a, j");

        // make antisymmetrized
        V_blks_["bbbb_ovov"]("i, a, j, b") = ijab_bbbb("i, j, a, b") - ibaj_bbbb("i, b, a, j");

        /*******************************************************************************************
         * <ia||jk> = (ij|ak) - (ik|aj)
         * i,j,k,l are occupied orbitals (o); a,b,c,d are virtual orbitals (v)
         *******************************************************************************************/

        /// aaaa - aaaa
        if (!pure_photon_integrals_) ijak_aaaa("i, j, a, k") = Qmo_blks["aa_oo"]("Q, i, j") * Qmo_blks["aa_vo"]("Q, a, k");
        else ijak_aaaa = HelperD::makeTensor(world_, {oa_, oa_, va_, oa_}, true);

        // add dipole terms
        ijak_aaaa("i, j, a, k") += lx2 * Dip_blks_["dx_aa_oo"]("i, j") * Dip_blks_["dx_aa_vo"]("a, k");
        ijak_aaaa("i, j, a, k") += ly2 * Dip_blks_["dy_aa_oo"]("i, j") * Dip_blks_["dy_aa_vo"]("a, k");
        ijak_aaaa("i, j, a, k") += lz2 * Dip_blks_["dz_aa_oo"]("i, j") * Dip_blks_["dz_aa_vo"]("a, k");

        // make antisymmetrized
        V_blks_["aaaa_ovoo"]("i, a, j, k") = ijak_aaaa("i, j, a, k") - ijak_aaaa("i, k, a, j");

        /// aabb - abba
        if (!pure_photon_integrals_) {
            ijak_aabb("i, j, a, k") = Qmo_blks["aa_oo"]("Q, i, j") * Qmo_blks["bb_vo"]("Q, a, k");
            ikaj_abba("i, k, a, j") = Qmo_blks["ab_oo"]("Q, i, k") * Qmo_blks["ba_vo"]("Q, a, j");
        } else {
            ijak_aabb = HelperD::makeTensor(world_, {oa_, oa_, vb_, ob_}, true);
            ikaj_abba = HelperD::makeTensor(world_, {oa_, ob_, vb_, oa_}, true);
        }

        // add dipole terms
        ijak_aabb("i, j, a, k") += lx2 * Dip_blks_["dx_aa_oo"]("i, j") * Dip_blks_["dx_bb_vo"]("a, k");
        ijak_aabb("i, j, a, k") += ly2 * Dip_blks_["dy_aa_oo"]("i, j") * Dip_blks_["dy_bb_vo"]("a, k");
        ijak_aabb("i, j, a, k") += lz2 * Dip_blks_["dz_aa_oo"]("i, j") * Dip_blks_["dz_bb_vo"]("a, k");

        // make antisymmetrized
        V_blks_["abab_ovoo"]("i, a, j, k") = ijak_aabb("i, j, a, k") - ikaj_abba("i, k, a, j");

        /// bbaa - baab
        if (!pure_photon_integrals_) {
            ijak_bbaa("i, j, a, k") = Qmo_blks["bb_oo"]("Q, i, j") * Qmo_blks["aa_vo"]("Q, a, k");
            ikaj_baab("i, k, a, j") = Qmo_blks["ba_oo"]("Q, i, k") * Qmo_blks["ab_vo"]("Q, a, j");
        } else {
            ijak_bbaa = HelperD::makeTensor(world_, {ob_, ob_, va_, oa_}, true);
            ikaj_baab = HelperD::makeTensor(world_, {ob_, oa_, va_, ob_}, true);
        }

        // add dipole terms
        ijak_bbaa("i, j, a, k") += lx2 * Dip_blks_["dx_bb_oo"]("i, j") * Dip_blks_["dx_aa_vo"]("a, k");
        ijak_bbaa("i, j, a, k") += ly2 * Dip_blks_["dy_bb_oo"]("i, j") * Dip_blks_["dy_aa_vo"]("a, k");
        ijak_bbaa("i, j, a, k") += lz2 * Dip_blks_["dz_bb_oo"]("i, j") * Dip_blks_["dz_aa_vo"]("a, k");

        // make antisymmetrized
        V_blks_["baba_ovoo"]("i, a, j, k") = ijak_bbaa("i, j, a, k") - ikaj_baab("i, k, a, j");

        /// bbbb - bbbb
        if (!pure_photon_integrals_) ijak_bbbb("i, j, a, k") = Qmo_blks["bb_oo"]("Q, i, j") * Qmo_blks["bb_vo"]("Q, a, k");
        else ijak_bbbb = HelperD::makeTensor(world_, {ob_, ob_, vb_, ob_}, true);

        // add dipole terms
        ijak_bbbb("i, j, a, k") += lx2 * Dip_blks_["dx_bb_oo"]("i, j") * Dip_blks_["dx_bb_vo"]("a, k");
        ijak_bbbb("i, j, a, k") += ly2 * Dip_blks_["dy_bb_oo"]("i, j") * Dip_blks_["dy_bb_vo"]("a, k");
        ijak_bbbb("i, j, a, k") += lz2 * Dip_blks_["dz_bb_oo"]("i, j") * Dip_blks_["dz_bb_vo"]("a, k");

        // make antisymmetrized
        V_blks_["bbbb_ovoo"]("i, a, j, k") = ijak_bbbb("i, j, a, k") - ijak_bbbb("i, k, a, j");

        /*******************************************************************************************
         * <jk||ia> = (ji|ka) - (ja|ki)
         * i,j,k,l are occupied orbitals (o); a,b,c,d are virtual orbitals (v)
         *******************************************************************************************/

        /// aaaa - aaaa
        if (!pure_photon_integrals_) {
            jika_aaaa("j, i, k, a") = Qmo_blks["aa_oo"]("Q, j, i") * Qmo_blks["aa_ov"]("Q, k, a");
            jaki_aaaa("j, a, k, i") = Qmo_blks["aa_ov"]("Q, j, a") * Qmo_blks["aa_oo"]("Q, k, i");
        } else {
            jika_aaaa = HelperD::makeTensor(world_, {oa_, oa_, oa_, va_}, true);
            jaki_aaaa = HelperD::makeTensor(world_, {oa_, va_, oa_, oa_}, true);
        }

        // add dipole terms
        jika_aaaa("j, i, k, a") += lx2 * Dip_blks_["dx_aa_oo"]("j, i") * Dip_blks_["dx_aa_ov"]("k, a");
        jika_aaaa("j, i, k, a") += ly2 * Dip_blks_["dy_aa_oo"]("j, i") * Dip_blks_["dy_aa_ov"]("k, a");
        jika_aaaa("j, i, k, a") += lz2 * Dip_blks_["dz_aa_oo"]("j, i") * Dip_blks_["dz_aa_ov"]("k, a");

        jaki_aaaa("j, a, k, i") += lx2 * Dip_blks_["dx_aa_ov"]("j, a") * Dip_blks_["dx_aa_oo"]("k, i");
        jaki_aaaa("j, a, k, i") += ly2 * Dip_blks_["dy_aa_ov"]("j, a") * Dip_blks_["dy_aa_oo"]("k, i");
        jaki_aaaa("j, a, k, i") += lz2 * Dip_blks_["dz_aa_ov"]("j, a") * Dip_blks_["dz_aa_oo"]("k, i");

        // make antisymmetrized
        V_blks_["aaaa_ooov"]("j, k, i, a") = jika_aaaa("j, i, k, a") - jaki_aaaa("j, a, k, i");

        /// aabb - abba
        if (!pure_photon_integrals_) {
            jika_aabb("j, i, k, a") = Qmo_blks["aa_oo"]("Q, j, i") * Qmo_blks["bb_ov"]("Q, k, a");
            jaki_abba("j, a, k, i") = Qmo_blks["ab_ov"]("Q, j, a") * Qmo_blks["ba_oo"]("Q, k, i");
        } else {
            jika_aabb = HelperD::makeTensor(world_, {oa_, oa_, ob_, vb_}, true);
            jaki_abba = HelperD::makeTensor(world_, {oa_, vb_, ob_, oa_}, true);
        }

        // add dipole terms
        jika_aabb("j, i, k, a") += lx2 * Dip_blks_["dx_aa_oo"]("j, i") * Dip_blks_["dx_bb_ov"]("k, a");
        jika_aabb("j, i, k, a") += ly2 * Dip_blks_["dy_aa_oo"]("j, i") * Dip_blks_["dy_bb_ov"]("k, a");
        jika_aabb("j, i, k, a") += lz2 * Dip_blks_["dz_aa_oo"]("j, i") * Dip_blks_["dz_bb_ov"]("k, a");

        // make antisymmetrized
        V_blks_["abab_ooov"]("j, k, i, a") = jika_aabb("j, i, k, a") - jaki_abba("j, a, k, i");

        /// bbaa - baab
        if (!pure_photon_integrals_) {
            jika_bbaa("j, i, k, a") = Qmo_blks["bb_oo"]("Q, j, i") * Qmo_blks["aa_ov"]("Q, k, a");
            jaki_baab("j, a, k, i") = Qmo_blks["ba_ov"]("Q, j, a") * Qmo_blks["ab_oo"]("Q, k, i");
        } else {
            jika_bbaa = HelperD::makeTensor(world_, {ob_, ob_, oa_, va_}, true);
            jaki_baab = HelperD::makeTensor(world_, {ob_, va_, oa_, ob_}, true);
        }

        // add dipole terms
        jika_bbaa("j, i, k, a") += lx2 * Dip_blks_["dx_bb_oo"]("j, i") * Dip_blks_["dx_aa_ov"]("k, a");
        jika_bbaa("j, i, k, a") += ly2 * Dip_blks_["dy_bb_oo"]("j, i") * Dip_blks_["dy_aa_ov"]("k, a");
        jika_bbaa("j, i, k, a") += lz2 * Dip_blks_["dz_bb_oo"]("j, i") * Dip_blks_["dz_aa_ov"]("k, a");

        // make antisymmetrized
        V_blks_["baba_ooov"]("j, k, i, a") = jika_bbaa("j, i, k, a") - jaki_baab("j, a, k, i");

        /// bbbb - bbbb
        if (!pure_photon_integrals_) {
            jika_bbbb("j, i, k, a") = Qmo_blks["bb_oo"]("Q, j, i") * Qmo_blks["bb_ov"]("Q, k, a");
            jaki_bbbb("j, a, k, i") = Qmo_blks["bb_ov"]("Q, j, a") * Qmo_blks["bb_oo"]("Q, k, i");
        } else {
            jika_bbbb = HelperD::makeTensor(world_, {ob_, ob_, ob_, vb_}, true);
            jaki_bbbb = HelperD::makeTensor(world_, {ob_, vb_, ob_, ob_}, true);
        }

        // add dipole terms
        jika_bbbb("j, i, k, a") += lx2 * Dip_blks_["dx_bb_oo"]("j, i") * Dip_blks_["dx_bb_ov"]("k, a");
        jika_bbbb("j, i, k, a") += ly2 * Dip_blks_["dy_bb_oo"]("j, i") * Dip_blks_["dy_bb_ov"]("k, a");
        jika_bbbb("j, i, k, a") += lz2 * Dip_blks_["dz_bb_oo"]("j, i") * Dip_blks_["dz_bb_ov"]("k, a");

        jaki_bbbb("j, a, k, i") += lx2 * Dip_blks_["dx_bb_ov"]("j, a") * Dip_blks_["dx_bb_oo"]("k, i");
        jaki_bbbb("j, a, k, i") += ly2 * Dip_blks_["dy_bb_ov"]("j, a") * Dip_blks_["dy_bb_oo"]("k, i");
        jaki_bbbb("j, a, k, i") += lz2 * Dip_blks_["dz_bb_ov"]("j, a") * Dip_blks_["dz_bb_oo"]("k, i");

        // make antisymmetrized
        V_blks_["bbbb_ooov"]("j, k, i, a") = jika_bbbb("j, i, k, a") - jaki_bbbb("j, a, k, i");

        /*******************************************************************************************
         * <ai||bc> = (ab|ic) - (ac|ib)
         * i,j,k,l are occupied orbitals (o); a,b,c,d are virtual orbitals (v)
         *******************************************************************************************/

        /// aaaa - aaaa
        if (!pure_photon_integrals_) abic_aaaa("a, b, i, c") = Qmo_blks["aa_vv"]("Q, a, b") * Qmo_blks["aa_ov"]("Q, i, c");
        else abic_aaaa = HelperD::makeTensor(world_, {va_, va_, oa_, va_}, true);

        // add dipole terms
        abic_aaaa("a, b, i, c") += lx2 * Dip_blks_["dx_aa_vv"]("a, b") * Dip_blks_["dx_aa_ov"]("i, c");
        abic_aaaa("a, b, i, c") += ly2 * Dip_blks_["dy_aa_vv"]("a, b") * Dip_blks_["dy_aa_ov"]("i, c");
        abic_aaaa("a, b, i, c") += lz2 * Dip_blks_["dz_aa_vv"]("a, b") * Dip_blks_["dz_aa_ov"]("i, c");

        // make antisymmetrized
        V_blks_["aaaa_vovv"]("a, i, b, c") = abic_aaaa("a, b, i, c") - abic_aaaa("a, c, i, b");

        /// aabb - abba
        if (!pure_photon_integrals_) {
            abic_aabb("a, b, i, c") = Qmo_blks["aa_vv"]("Q, a, b") * Qmo_blks["bb_ov"]("Q, i, c");
            acib_abba("a, c, i, b") = Qmo_blks["ab_vv"]("Q, a, c") * Qmo_blks["ba_ov"]("Q, i, b");
        } else {
            abic_aabb = HelperD::makeTensor(world_, {va_, va_, ob_, vb_}, true);
            acib_abba = HelperD::makeTensor(world_, {va_, vb_, ob_, va_}, true);
        }

        // add dipole terms
        abic_aabb("a, b, i, c") += lx2 * Dip_blks_["dx_aa_vv"]("a, b") * Dip_blks_["dx_bb_ov"]("i, c");
        abic_aabb("a, b, i, c") += ly2 * Dip_blks_["dy_aa_vv"]("a, b") * Dip_blks_["dy_bb_ov"]("i, c");
        abic_aabb("a, b, i, c") += lz2 * Dip_blks_["dz_aa_vv"]("a, b") * Dip_blks_["dz_bb_ov"]("i, c");

        // make antisymmetrized
        V_blks_["abab_vovv"]("a, i, b, c") = abic_aabb("a, b, i, c") - acib_abba("a, c, i, b");

        /// bbaa - baab
        if (!pure_photon_integrals_) {
            abic_bbaa("a, b, i, c") = Qmo_blks["bb_vv"]("Q, a, b") * Qmo_blks["aa_ov"]("Q, i, c");
            acib_baab("a, c, i, b") = Qmo_blks["ba_vv"]("Q, a, c") * Qmo_blks["ab_ov"]("Q, i, b");
        } else {
            abic_bbaa = HelperD::makeTensor(world_, {vb_, vb_, oa_, va_}, true);
            acib_baab = HelperD::makeTensor(world_, {vb_, va_, oa_, vb_}, true);
        }

        // add dipole terms
        abic_bbaa("a, b, i, c") += lx2 * Dip_blks_["dx_bb_vv"]("a, b") * Dip_blks_["dx_aa_ov"]("i, c");
        abic_bbaa("a, b, i, c") += ly2 * Dip_blks_["dy_bb_vv"]("a, b") * Dip_blks_["dy_aa_ov"]("i, c");
        abic_bbaa("a, b, i, c") += lz2 * Dip_blks_["dz_bb_vv"]("a, b") * Dip_blks_["dz_aa_ov"]("i, c");

        // make antisymmetrized
        V_blks_["baba_vovv"]("a, i, b, c") = abic_bbaa("a, b, i, c") - acib_baab("a, c, i, b");

        /// bbbb - bbbb
        if (!pure_photon_integrals_) abic_bbbb("a, b, i, c") = Qmo_blks["bb_vv"]("Q, a, b") * Qmo_blks["bb_ov"]("Q, i, c");
        else abic_bbbb = HelperD::makeTensor(world_, {vb_, vb_, ob_, vb_}, true);

        // add dipole terms
        abic_bbbb("a, b, i, c") += lx2 * Dip_blks_["dx_bb_vv"]("a, b") * Dip_blks_["dx_bb_ov"]("i, c");
        abic_bbbb("a, b, i, c") += ly2 * Dip_blks_["dy_bb_vv"]("a, b") * Dip_blks_["dy_bb_ov"]("i, c");
        abic_bbbb("a, b, i, c") += lz2 * Dip_blks_["dz_bb_vv"]("a, b") * Dip_blks_["dz_bb_ov"]("i, c");

        // make antisymmetrized
        V_blks_["bbbb_vovv"]("a, i, b, c") = abic_bbbb("a, b, i, c") - abic_bbbb("a, c, i, b");


        /*******************************************************************************************
         * <ab||ic> = (ai|bc) - (ac|bi)
         *******************************************************************************************/

        /// aaaa - aaaa
        if (!pure_photon_integrals_) {
            aibc_aaaa("a, i, b, c") = Qmo_blks["aa_vo"]("Q, a, i") * Qmo_blks["aa_vv"]("Q, b, c");
            acbi_aaaa("a, c, b, i") = Qmo_blks["aa_vv"]("Q, a, c") * Qmo_blks["aa_vo"]("Q, b, i");
        } else {
            aibc_aaaa = HelperD::makeTensor(world_, {va_, oa_, va_, va_}, true);
            acbi_aaaa = HelperD::makeTensor(world_, {va_, va_, va_, oa_}, true);
        }

        // add dipole terms
        aibc_aaaa("a, i, b, c") += lx2 * Dip_blks_["dx_aa_vo"]("a, i") * Dip_blks_["dx_aa_vv"]("b, c");
        aibc_aaaa("a, i, b, c") += ly2 * Dip_blks_["dy_aa_vo"]("a, i") * Dip_blks_["dy_aa_vv"]("b, c");
        aibc_aaaa("a, i, b, c") += lz2 * Dip_blks_["dz_aa_vo"]("a, i") * Dip_blks_["dz_aa_vv"]("b, c");

        acbi_aaaa("a, c, b, i") += lx2 * Dip_blks_["dx_aa_vv"]("a, c") * Dip_blks_["dx_aa_vo"]("b, i");
        acbi_aaaa("a, c, b, i") += ly2 * Dip_blks_["dy_aa_vv"]("a, c") * Dip_blks_["dy_aa_vo"]("b, i");
        acbi_aaaa("a, c, b, i") += lz2 * Dip_blks_["dz_aa_vv"]("a, c") * Dip_blks_["dz_aa_vo"]("b, i");

        // make antisymmetrized
        V_blks_["aaaa_vvov"]("a, b, i, c") = aibc_aaaa("a, i, b, c") - acbi_aaaa("a, c, b, i");

        /// aabb - abba
        if (!pure_photon_integrals_) {
            aibc_aabb("a, i, b, c") = Qmo_blks["aa_vo"]("Q, a, i") * Qmo_blks["bb_vv"]("Q, b, c");
            acbi_abba("a, c, b, i") = Qmo_blks["ab_vv"]("Q, a, c") * Qmo_blks["ba_vo"]("Q, b, i");
        } else {
            aibc_aabb = HelperD::makeTensor(world_, {va_, oa_, vb_, vb_}, true);
            acbi_abba = HelperD::makeTensor(world_, {va_, vb_, vb_, oa_}, true);
        }

        // add dipole terms
        aibc_aabb("a, i, b, c") += lx2 * Dip_blks_["dx_aa_vo"]("a, i") * Dip_blks_["dx_bb_vv"]("b, c");
        aibc_aabb("a, i, b, c") += ly2 * Dip_blks_["dy_aa_vo"]("a, i") * Dip_blks_["dy_bb_vv"]("b, c");
        aibc_aabb("a, i, b, c") += lz2 * Dip_blks_["dz_aa_vo"]("a, i") * Dip_blks_["dz_bb_vv"]("b, c");

        // make antisymmetrized
        V_blks_["abab_vvov"]("a, b, i, c") = aibc_aabb("a, i, b, c") - acbi_abba("a, c, b, i");

        /// bbaa - baab
        if (!pure_photon_integrals_) {
            aibc_bbaa("a, i, b, c") = Qmo_blks["bb_vo"]("Q, a, i") * Qmo_blks["aa_vv"]("Q, b, c");
            acbi_baab("a, c, b, i") = Qmo_blks["ba_vv"]("Q, a, c") * Qmo_blks["ab_vo"]("Q, b, i");
        } else {
            aibc_bbaa = HelperD::makeTensor(world_, {vb_, ob_, va_, va_}, true);
            acbi_baab = HelperD::makeTensor(world_, {vb_, va_, va_, ob_}, true);
        }

        // add dipole terms
        aibc_bbaa("a, i, b, c") += lx2 * Dip_blks_["dx_bb_vo"]("a, i") * Dip_blks_["dx_aa_vv"]("b, c");
        aibc_bbaa("a, i, b, c") += ly2 * Dip_blks_["dy_bb_vo"]("a, i") * Dip_blks_["dy_aa_vv"]("b, c");
        aibc_bbaa("a, i, b, c") += lz2 * Dip_blks_["dz_bb_vo"]("a, i") * Dip_blks_["dz_aa_vv"]("b, c");

        // make antisymmetrized
        V_blks_["baba_vvov"]("a, b, i, c") = aibc_bbaa("a, i, b, c") - acbi_baab("a, c, b, i");

        /// bbbb - bbbb
        if (!pure_photon_integrals_) {
            aibc_bbbb("a, i, b, c") = Qmo_blks["bb_vo"]("Q, a, i") * Qmo_blks["bb_vv"]("Q, b, c");
            acbi_bbbb("a, c, b, i") = Qmo_blks["bb_vv"]("Q, a, c") * Qmo_blks["bb_vo"]("Q, b, i");
        } else {
            aibc_bbbb = HelperD::makeTensor(world_, {vb_, ob_, vb_, vb_}, true);
            acbi_bbbb = HelperD::makeTensor(world_, {vb_, vb_, vb_, ob_}, true);
        }

        // add dipole terms
        aibc_bbbb("a, i, b, c") += lx2 * Dip_blks_["dx_bb_vo"]("a, i") * Dip_blks_["dx_bb_vv"]("b, c");
        aibc_bbbb("a, i, b, c") += ly2 * Dip_blks_["dy_bb_vo"]("a, i") * Dip_blks_["dy_bb_vv"]("b, c");
        aibc_bbbb("a, i, b, c") += lz2 * Dip_blks_["dz_bb_vo"]("a, i") * Dip_blks_["dz_bb_vv"]("b, c");

        acbi_bbbb("a, c, b, i") += lx2 * Dip_blks_["dx_bb_vv"]("a, c") * Dip_blks_["dx_bb_vo"]("b, i");
        acbi_bbbb("a, c, b, i") += ly2 * Dip_blks_["dy_bb_vv"]("a, c") * Dip_blks_["dy_bb_vo"]("b, i");
        acbi_bbbb("a, c, b, i") += lz2 * Dip_blks_["dz_bb_vv"]("a, c") * Dip_blks_["dz_bb_vo"]("b, i");

        // make antisymmetrized
        V_blks_["bbbb_vvov"]("a, b, i, c") = aibc_bbbb("a, i, b, c") - acbi_bbbb("a, c, b, i");


        // make new ordering
        V_blks_["aaaa_oovo"]("p,q,r,s") =  1.0 * V_blks_["aaaa_ooov"]("q,p,s,r");
        V_blks_["abba_oovo"]("p,q,r,s") = -1.0 * V_blks_["abab_ooov"]("p,q,s,r");
        V_blks_["abab_oovo"]("p,q,r,s") =  1.0 * V_blks_["baba_ooov"]("q,p,s,r");
        V_blks_["bbbb_oovo"]("p,q,r,s") =  1.0 * V_blks_["bbbb_ooov"]("q,p,s,r");

        V_blks_["aaaa_vooo"]("p,q,r,s") =  1.0 * V_blks_["aaaa_ovoo"]("q,p,s,r");
        V_blks_["abba_vooo"]("p,q,r,s") = -1.0 * V_blks_["baba_ovoo"]("q,p,r,s");
        V_blks_["abab_vooo"]("p,q,r,s") =  1.0 * V_blks_["baba_ovoo"]("q,p,s,r");
        V_blks_["baab_vooo"]("p,q,r,s") = -1.0 * V_blks_["abab_ovoo"]("q,p,r,s");
        V_blks_["bbbb_vooo"]("p,q,r,s") =  1.0 * V_blks_["bbbb_ovoo"]("q,p,s,r");

        V_blks_["aaaa_vovo"]("p,q,r,s") =  1.0 * V_blks_["aaaa_ovov"]("q,p,s,r");
        V_blks_["abab_vovo"]("p,q,r,s") =  1.0 * V_blks_["baba_ovov"]("q,p,s,r");
        V_blks_["baba_vovo"]("p,q,r,s") =  1.0 * V_blks_["abab_ovov"]("q,p,s,r");
        V_blks_["abba_vovo"]("p,q,r,s") =  1.0 * V_blks_["baab_ovov"]("q,p,s,r");
        V_blks_["baab_vovo"]("p,q,r,s") =  1.0 * V_blks_["abba_ovov"]("q,p,s,r");
        V_blks_["bbbb_vovo"]("p,q,r,s") =  1.0 * V_blks_["bbbb_ovov"]("q,p,s,r");

        V_blks_["baab_vovv"]("p,q,r,s") = -1.0 * V_blks_["baba_vovv"]("p,q,s,r");

        V_blks_["aaaa_vvvo"]("p,q,r,s") =  1.0 * V_blks_["aaaa_vvov"]("q,p,s,r");
        V_blks_["abba_vvvo"]("p,q,r,s") = -1.0 * V_blks_["abab_vvov"]("p,q,s,r");
        V_blks_["abab_vvvo"]("p,q,r,s") =  1.0 * V_blks_["baba_vvov"]("q,p,s,r");
        V_blks_["bbbb_vvvo"]("p,q,r,s") =  1.0 * V_blks_["bbbb_vvov"]("q,p,s,r");

        /*******************************************************************************************
         * All electron integrals are now computed.
         *******************************************************************************************/
    }


} // cc_cavity
