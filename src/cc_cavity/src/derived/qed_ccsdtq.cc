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

#include "../../include/derived/qed_ccsdtq.h"
#include <tiledarray.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libfunctional/superfunctional.h>
#include <psi4/libmints/mintshelper.h>
#include <psi4/libmints/basisset.h>
#include "../../misc/ta_helper.h"
#include <omp.h>

namespace hilbert {

    QED_CCSDTQ::QED_CCSDTQ(const shared_ptr<Wavefunction> &reference_wavefunction, Options &options) :
                             QED_CCSDT(reference_wavefunction, options) {
    }

    void QED_CCSDTQ::init_operators() {

        // call base class for t1, t2, u0, u1, u2, and other operators
        QED_CCSDT::init_operators();

        /// initialize amplitude and residual blocks

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

    void QED_CCSDTQ::update_amplitudes() {

        QED_CCSDT::update_amplitudes(); // call parent function for t1 and t2

        /// dt = -residual / eps

        double *eps = epsilon_;
        size_t o = o_;
        size_t v = v_;
        size_t oa = oa_;
        size_t va = va_;

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

        /// update amplitudes according to t + dt = amplitude - residual / eps
        if (include_t4_) {
            amplitudes_["t4_aaaaaaaa"](idx_map_[8]) += residuals_["t4_aaaaaaaa"](idx_map_[8]);
            amplitudes_["t4_aaabaaab"](idx_map_[8]) += residuals_["t4_aaabaaab"](idx_map_[8]);
            amplitudes_["t4_aabbaabb"](idx_map_[8]) += residuals_["t4_aabbaabb"](idx_map_[8]);
            amplitudes_["t4_abbbabbb"](idx_map_[8]) += residuals_["t4_abbbabbb"](idx_map_[8]);
            amplitudes_["t4_bbbbbbbb"](idx_map_[8]) += residuals_["t4_bbbbbbbb"](idx_map_[8]);
        }

        /// update amplitudes according to u + du = amplitude - residual / (eps + w)
        if (include_u4_) {
            amplitudes_["u4_aaaaaaaa"](idx_map_[8]) += residuals_["u4_aaaaaaaa"](idx_map_[8]);
            amplitudes_["u4_aaabaaab"](idx_map_[8]) += residuals_["u4_aaabaaab"](idx_map_[8]);
            amplitudes_["u4_aabbaabb"](idx_map_[8]) += residuals_["u4_aabbaabb"](idx_map_[8]);
            amplitudes_["u4_abbbabbb"](idx_map_[8]) += residuals_["u4_abbbabbb"](idx_map_[8]);
            amplitudes_["u4_bbbbbbbb"](idx_map_[8]) += residuals_["u4_bbbbbbbb"](idx_map_[8]);
        }

        world_.gop.fence();

    }

    double QED_CCSDTQ::compute_residual_norms(bool return_tot) {

        QED_CCSDT::compute_residual_norms(false);

        /// residual norms for each amplitude
        if (include_t4_) {
            resid_norms_["t4"] = sqrt(squared_norm(residuals_["t4_aaaaaaaa"])
                                      + squared_norm(residuals_["t4_aaabaaab"])
                                      + squared_norm(residuals_["t4_aabbaabb"])
                                      + squared_norm(residuals_["t4_abbbabbb"])
                                      + squared_norm(residuals_["t4_bbbbbbbb"]));
        }
        if (include_u4_) {
            resid_norms_["u4"] = sqrt(squared_norm(residuals_["u4_aaaaaaaa"])
                                      + squared_norm(residuals_["u4_aaabaaab"])
                                      + squared_norm(residuals_["u4_aabbaabb"])
                                      + squared_norm(residuals_["u4_abbbabbb"])
                                      + squared_norm(residuals_["u4_bbbbbbbb"]));
        }

        // total residual norm for all amplitudes (if requested)
        if (return_tot) {
            double norm = 0.0;
            for (auto & amp : residuals_) {
                if (!amp.second.is_initialized()) continue;
                norm += squared_norm(amp.second);
            }
            return sqrt(norm);
        }
        else return 0.0;
    }

    double QED_CCSDTQ::build_residuals(){
        throw PsiException("QED_CCSDTQ::build_residuals() not implemented", __FILE__, __LINE__);
    }

    void QED_CCSDTQ::print_properties() {
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
        if ( include_t4_ ) outfile->Printf("\n    T4: %15.12lf | %5.2f %%", sqrt(nT4), 100.0*nT4/total_norm);
        if ( include_u0_ ) outfile->Printf("\n    U0: %15.12lf | %5.2f %%", u0_, 100.0*fabs(u0_*u0_)/total_norm);
        if ( include_u1_ ) outfile->Printf("\n    U1: %15.12lf | %5.2f %%", sqrt(nU1), 100.0*nU1/total_norm);
        if ( include_u2_ ) outfile->Printf("\n    U2: %15.12lf | %5.2f %%", sqrt(nU2), 100.0*nU2/total_norm);
        if ( include_u3_ ) outfile->Printf("\n    U3: %15.12lf | %5.2f %%", sqrt(nU3), 100.0*nU3/total_norm);
        if ( include_u4_ ) outfile->Printf("\n    U4: %15.12lf | %5.2f %%", sqrt(nU4), 100.0*nU4/total_norm);
        outfile->Printf("\n   ----------------------------");
        outfile->Printf("\n    Total: %15.12lf\n\n", sqrt(total_norm));

    }

    void QED_CCSDTQ::print_iter_header() const {
        Printf("\n");
        Printf("    ==>  Begin %s iterations <==    \n", cc_type_.c_str());
        Printf("\n");
        Printf("%5s %16s %15s %15s | %8s %8s",  "Iter","energy","dE","|dT|","|dT1|","|dT2|");
        if (include_t3_) Printf(" %8s","|dT3|");
        if (include_t4_) Printf(" %8s","|dT4|");
        if (include_u0_) Printf(" %8s","|dU0|");
        if (include_u1_) Printf(" %8s","|dU1|");
        if (include_u2_) Printf(" %8s","|dU2|");
        if (include_u3_) Printf(" %8s","|dU3|");
        if (include_u4_) Printf(" %8s","|dU4|");
        Printf("\n");
    }

    void QED_CCSDTQ::print_iteration(size_t iter, double energy, double dele, double tnorm) const {
        Printf("%5i %17.12lf %15.12lf %15.12lf | %-8.1e %-8.1e",iter,energy,dele,tnorm,resid_norms_.at("t1"),resid_norms_.at("t2"));
        if (include_t3_) Printf(" %-8.1e", resid_norms_.at("t3"));
        if (include_t4_) Printf(" %-8.1e", resid_norms_.at("t4"));
        if (include_u0_) Printf(" %-8.1e", resid_norms_.at("u0"));
        if (include_u1_) Printf(" %-8.1e", resid_norms_.at("u1"));
        if (include_u2_) Printf(" %-8.1e", resid_norms_.at("u2"));
        if (include_u3_) Printf(" %-8.1e", resid_norms_.at("u3"));
        if (include_u4_) Printf(" %-8.1e", resid_norms_.at("u4"));
        Printf("\n");
    }

} // cc_cavity
