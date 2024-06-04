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

#include "cc_cavity/include/derived/ccsd.h"
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

#include "cc_cavity/misc/ta_helper.hpp"

#include <mkl.h>
#include <omp.h>
#include "misc/blas.h"
#include "misc/hilbert_psifiles.h"
#include "polaritonic_scf/uhf.h"
#include <unistd.h>
#include <psi4/psifiles.h>

using namespace TA_Helper;

namespace hilbert {

    CCSD::CCSD(const shared_ptr<Wavefunction> &reference_wavefunction, Options &options) :
                             CC_Cavity(reference_wavefunction, options) {
    }


    void CCSD::init_operators() {

        // call base class for building delta operators
        CC_Cavity::init_operators();

        /// initialize t1 / t2 amplitudes and residuals

        // t1
        amplitudes_["t1_aa"] = makeTensor(world_, {va_, oa_}, true);
        amplitudes_["t1_bb"] = makeTensor(world_, {vb_, ob_}, true);
        residuals_["t1_aa"]  = makeTensor(world_, {va_, oa_}, true);
        residuals_["t1_bb"]  = makeTensor(world_, {vb_, ob_}, true);

        // t2
        amplitudes_["t2_aaaa"] = makeTensor(world_, {va_,va_, oa_,oa_}, true);
        amplitudes_["t2_abab"] = makeTensor(world_, {va_,vb_, oa_,ob_}, true);
        amplitudes_["t2_bbbb"] = makeTensor(world_, {vb_,vb_, ob_,ob_}, true);
        residuals_["t2_aaaa"]  = makeTensor(world_, {va_,va_, oa_,oa_}, true);
        residuals_["t2_abab"]  = makeTensor(world_, {va_,vb_, oa_,ob_}, true);
        residuals_["t2_bbbb"]  = makeTensor(world_, {vb_,vb_, ob_,ob_}, true);

    }

    void CCSD::update_amplitudes() {

        double *eps = epsilon_;
        size_t o = o_;
        size_t v = v_;
        size_t oa = oa_;
        size_t va = va_;

        /// dt = -residual / eps

        // t1
        forall(residuals_["t1_aa"],
               [eps, o, oa, va](auto &tile, auto &x) {
                   tile[x] /= (eps[x[1]] - eps[x[0]+o]);
               });

        forall(residuals_["t1_bb"],
               [eps, o, oa, va](auto &tile, auto &x) {
                   tile[x] /= (eps[x[1]+oa] - eps[x[0]+o+va]);
               });

        // t2
        forall(residuals_["t2_aaaa"],
               [eps, o, oa, va](auto &tile, auto &x) {
                   double o_ep = eps[x[2]] + eps[x[3]],
                           v_ep = eps[x[0]+o] + eps[x[1]+o];
                   tile[x] /= (o_ep - v_ep);
               });
        forall(residuals_["t2_bbbb"],
               [eps, o, oa, va](auto &tile, auto &x) {
                   double o_ep = eps[x[2]+oa] + eps[x[3]+oa],
                           v_ep = eps[x[0]+o+va] + eps[x[1]+o+va];
                   tile[x] /= (o_ep - v_ep);
               });
        forall(residuals_["t2_abab"],
               [eps, o, oa, va](auto &tile, auto &x) {
                   double o_ep = eps[x[2]] + eps[x[3]+oa],
                           v_ep = eps[x[0]+o] + eps[x[1]+o];
                   tile[x] /= (o_ep - v_ep);
               });

        world_.gop.fence();

        /// update amplitudes according to t + dt = amplitude - residual / eps

        amplitudes_["t1_aa"]("a,i")   += residuals_["t1_aa"]("a,i");
        amplitudes_["t1_bb"]("a,i")   += residuals_["t1_bb"]("a,i");
        amplitudes_["t2_aaaa"]("a,b,i,j") += residuals_["t2_aaaa"]("a,b,i,j");
        amplitudes_["t2_abab"]("a,b,i,j") += residuals_["t2_abab"]("a,b,i,j");
        amplitudes_["t2_bbbb"]("a,b,i,j") += residuals_["t2_bbbb"]("a,b,i,j");

        world_.gop.fence();

    }

    double CCSD::compute_residual_norms(bool return_tot) {

        /// residual norms for each amplitude

        resid_norms_["t1"] = sqrt(squared_norm(residuals_["t1_aa"])
                                  + squared_norm(residuals_["t1_bb"]));
        resid_norms_["t2"] = sqrt(squared_norm(residuals_["t2_aaaa"])
                                  + squared_norm(residuals_["t2_abab"])
                                  + squared_norm(residuals_["t2_bbbb"]));

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
    
    void CCSD::print_properties() {
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

        // print norms of cluster amplitudes
        outfile->Printf("\n\n   Norms of cluster amplitudes:");
        outfile->Printf("\n   ------------------------------");
        outfile->Printf("\n    T1: %15.12lf | %5.2f %%", sqrt(nT1), 100*nT1/total_norm);
        outfile->Printf("\n    T2: %15.12lf | %5.2f %%", sqrt(nT2), 100*nT2/total_norm);
        outfile->Printf("\n   ------------------------------");
        outfile->Printf("\n    Total: %15.12lf\n\n", sqrt(total_norm));
    }

    void CCSD::print_iter_header() const {
        Printf("\n");
        Printf("    ==>  Begin %s iterations <==    \n", cc_type_.c_str());
        Printf("\n");
        Printf("%5s %16s %15s %15s  | %8s %8s",  "Iter","energy","dE","|dT|","|dT1|","|dT2|\n");
    }

    void CCSD::print_iteration(size_t iter, double energy, double dele, double tnorm) const {
        Printf("%5i %17.12lf %15.12lf %15.12lf | %-8.1e %-8.1e",iter,energy,dele,tnorm,
               resid_norms_.at("t1"),
               resid_norms_.at("t2"));
        Printf("\n");
    }

} // cc_cavity
