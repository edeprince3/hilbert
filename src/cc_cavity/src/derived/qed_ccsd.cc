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

#include "../../include/derived/qed_ccsd.h"
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

    QED_CCSD::QED_CCSD(const shared_ptr<Wavefunction> &reference_wavefunction, Options &options) :
                             CC_Cavity(reference_wavefunction, options) {
    }


    void QED_CCSD::init_operators() {

        // call base class for t1 and t2, plus some other stuff
        CC_Cavity::init_operators();

        /// initialize amplitude and residual blocks

        // u0
        if (include_u0_){
            amplitudes_["u0"] = HelperD::makeTensor(world_, {1}, true);
            residuals_["u0"] = HelperD::makeTensor(world_, {1}, true);

            // should remove this and just use amplitudes_["u0"]. I need both for now.
            scalar_amps_["u0"] = 0.0;
            scalar_resids_["u0"] = 0.0;
        }


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

    }

    void QED_CCSD::update_amplitudes() {

        CC_Cavity::update_amplitudes(); // call parent function for t1 and t2

        /// dt = -residual / eps

        double *eps = epsilon_;
        size_t o = o_;
        size_t v = v_;
        size_t oa = oa_;
        size_t va = va_;

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


        world_.gop.fence();
        double ru0_ = scalar_resids_["u0"];
        double u0_ = scalar_amps_["u0"];

        // u0
        if ( include_u0_ ) {
            ru0_ /= -w0;
            amplitudes_["u0"] = HelperD::makeTensor(world_, {1}, &u0_);
             residuals_["u0"] = HelperD::makeTensor(world_, {1}, &ru0_);
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

        world_.gop.fence();

    }

    double QED_CCSD::compute_residual_norms(bool return_tot) {

        /// residual norms for each amplitude

        // call parent function for t1 and t2
        CC_Cavity::compute_residual_norms();

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
    
    void QED_CCSD::print_properties() {
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
        double nU1 = amp_norms["u1"];
        double nU2 = amp_norms["u2"];

        double u0_ = scalar_amps_["u0"];

        // print norms of cluster amplitudes
        outfile->Printf("\n\n   Norms of cluster amplitudes:");
        outfile->Printf("\n   ------------------------------");
        outfile->Printf("\n    T1: %15.12lf | %5.2f %%"
                        "\n    T2: %15.12lf | %5.2f %%",
                        sqrt(nT1), 100*nT1/total_norm,
                        sqrt(nT2), 100*nT2/total_norm);
        if ( include_u0_ ) outfile->Printf("\n    U0: %15.12lf | %5.2f %%", u0_, 100.0*fabs(u0_*u0_)/total_norm);
        if ( include_u1_ ) outfile->Printf("\n    U1: %15.12lf | %5.2f %%", sqrt(nU1), 100.0*nU1/total_norm);
        if ( include_u2_ ) outfile->Printf("\n    U2: %15.12lf | %5.2f %%", sqrt(nU2), 100.0*nU2/total_norm);
        outfile->Printf("\n   ------------------------------");
        outfile->Printf("\n    Total: %15.12lf\n\n", sqrt(total_norm));

    }

    void QED_CCSD::print_iter_header() const {
        Printf("\n");
        Printf("    ==>  Begin %s iterations <==    \n", cc_type_.c_str());
        Printf("\n");
        Printf("%5s %16s %15s %15s  | %8s %8s",  "Iter","energy","dE","|dT|","|dT1|","|dT2|");
        if (include_u0_) Printf(" %8s","|dU0|");
        if (include_u1_) Printf(" %8s","|dU1|");
        if (include_u2_) Printf(" %8s","|dU2|");
        Printf("\n");
    }

    void QED_CCSD::print_iteration(size_t iter, double energy, double dele, double tnorm) const {
        Printf("%5i %17.12lf %15.12lf %15.12lf | %-8.1e %-8.1e",iter,energy,dele,tnorm,resid_norms_.at("t1"),resid_norms_.at("t2"));
        if (include_u0_) Printf(" %-8.1e", resid_norms_.at("u0"));
        if (include_u1_) Printf(" %-8.1e", resid_norms_.at("u1"));
        if (include_u2_) Printf(" %-8.1e", resid_norms_.at("u2"));
        Printf("\n");
    }

} // cc_cavity
