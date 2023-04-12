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

} // cc_cavity
