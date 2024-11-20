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

#include "cc_cavity/include/qed_ccsd_21/qed_ccsd_21.h"
#include <psi4/libfunctional/superfunctional.h>
#include <psi4/libmints/basisset.h>
#include <psi4/libmints/mintshelper.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/psi4-dec.h>
#include <unistd.h>
#include "cc_cavity/misc/ta_helper.h"
#include "tiledarray.h"

using namespace TA_Helper;

namespace hilbert {

    QED_CCSD_21::QED_CCSD_21(const shared_ptr<Wavefunction> &reference_wavefunction, Options &options, map<string,bool> &includes) :
                             CCSD(reference_wavefunction, options, includes) {
    }


    void QED_CCSD_21::init_operators() {

        // call base class for t1 and t2, plus some other stuff
        CCSD::init_operators();

        /// initialize amplitude and residual blocks

        // t0_1
        amplitudes_["t0_1"] = makeTensor(world_, {1});
        residuals_["t0_1"]  = makeTensor(world_, {1});
        idxs_["t0_1"] = "0";

        // t1_1
        amplitudes_["t1_1_aa"]  = makeTensor(world_, {va_, oa_});
        amplitudes_["t1_1_bb"]  = makeTensor(world_, {vb_, ob_});
        residuals_["t1_1_aa"]   = makeTensor(world_, {va_, oa_});
        residuals_["t1_1_bb"]   = makeTensor(world_, {vb_, ob_});
        idxs_["t1_1_aa"] = "a,i"; idxs_["t1_1_bb"] = "a,i";

        // t2_1
        amplitudes_["t2_1_aaaa"]  = makeTensor(world_, {va_,va_, oa_,oa_});
        amplitudes_["t2_1_abab"]  = makeTensor(world_, {va_,vb_, oa_,ob_});
        amplitudes_["t2_1_bbbb"]  = makeTensor(world_, {vb_,vb_, ob_,ob_});
        residuals_["t2_1_aaaa"]   = makeTensor(world_, {va_,va_, oa_,oa_});
        residuals_["t2_1_abab"]   = makeTensor(world_, {va_,vb_, oa_,ob_});
        residuals_["t2_1_bbbb"]   = makeTensor(world_, {vb_,vb_, ob_,ob_});
        idxs_["t2_1_aaaa"] = "a,b,i,j"; idxs_["t2_1_abab"] = "a,b,i,j"; idxs_["t2_1_bbbb"] = "a,b,i,j";

    }

    void QED_CCSD_21::update_residuals() {

        CCSD::update_residuals(); // call parent function for t1 and t2

        /// dt = -residual / eps

        double *eps = epsilon_;
        size_t o = o_;
        size_t oa = oa_;
        size_t va = va_;

        /// du = -residual / (eps + w)
        double w0 = cavity_frequency_[2];

        // t0_1
        foreach_inplace(residuals_["t0_1"], [w0](auto &tile) {
            for (auto &x : tile.range()) {
                if (fabs(w0) > 1.0e-12)
                    tile[x] /= w0;
                else tile[x] = 0.0;
            }
        });

        // t1_1
        foreach_inplace(residuals_["t1_1_aa"], [eps, o, w0](auto &tile) {
                            for (auto &x : tile.range()) {
                                size_t a = x[0] + o, i = x[1];
                                tile[x] /= eps[a] + w0 - eps[i];
                            }
                        });

        foreach_inplace(residuals_["t1_1_bb"], [eps, o, oa, va, w0](auto &tile) {
                            for (auto &x : tile.range()) {
                                size_t a = x[0] + o + va, i = x[1] + oa;
                                tile[x] /= eps[a] + w0 - eps[i];
                            }
                        });
        // t2_1
        foreach_inplace(residuals_["t2_1_aaaa"], [eps, o, w0](auto &tile) {
                            for (auto &x : tile.range()) {
                                size_t a = x[0] + o, b = x[1] + o, i = x[2], j = x[3];
                                tile[x] /= eps[a] + eps[b] + w0 - eps[i] - eps[j];
                            }
                        });
        foreach_inplace(residuals_["t2_1_bbbb"], [eps, o, oa, va, w0](auto &tile) {
                            for (auto &x : tile.range()) {
                                size_t a = x[0] + o+va, b = x[1] + o+va, i = x[2] + oa, j = x[3] + oa;
                                tile[x] /= eps[a] + eps[b] + w0 - eps[i] - eps[j];
                            }
                        });
        foreach_inplace(residuals_["t2_1_abab"], [eps, o, oa, va, w0](auto &tile) {
                            for (auto &x : tile.range()) {
                                size_t a = x[0] + o, b = x[1] + o+va, i = x[2], j = x[3] + oa;
                                tile[x] /= eps[a] + eps[b] + w0 - eps[i] - eps[j];
                            }
                        });

        world_.gop.fence();

    }

    double QED_CCSD_21::compute_residual_norms(bool return_tot) {

        /// residual norms for each amplitude

        // call parent function for t1 and t2
        CCSD::compute_residual_norms(false);

        // t_1 residual norms
        if (includes_.at("t0_1")) {
            resid_norms_["t0_1"] = sqrt(squared_norm(residuals_["t0_1"]));
        }
        if (includes_.at("t1_1")) {
            resid_norms_["t1_1"] = sqrt(squared_norm(residuals_["t1_1_aa"])
                                      + squared_norm(residuals_["t1_1_bb"]));
        }
        if (includes_.at("t2_1")) {
            resid_norms_["t2_1"] = sqrt(squared_norm(residuals_["t2_1_aaaa"])
                                      + squared_norm(residuals_["t2_1_abab"])
                                      + squared_norm(residuals_["t2_1_bbbb"]));
        }

        // total residual norm for all amplitudes (if requested)
        if (return_tot) {
            double norm = 0.0;
            for (auto & [name, amp] : residuals_) {
                if (!amp.is_initialized()) continue;
                norm += squared_norm(amp);
            }
            return sqrt(norm);
        }

        return 0.0;
    }
    
    void QED_CCSD_21::print_properties() {
        // calculate norms of amplitudes
        map<string, double> amp_norms;
        double total_norm = 0.0;
        for (auto& [name, amp] : amplitudes_) {
            if(!amp.is_initialized()) continue;
            // name of amplitude are all characters before the last '_a' or '_b'
            auto pos = name.find("_a");
            if (pos == string::npos)
                 pos = name.find("_b");
            if (pos == string::npos) {
                pos = name.size();
            }

            string amp_name = name.substr(0, pos);
            double amp_norm = squared_norm(amp);
            amp_norms[amp_name] += amp_norm;
            total_norm += amp_norm;
        }

        double nT1   = amp_norms["t1"];
        double nT2   = amp_norms["t2"];
        double nT01 = amp_norms["t0_1"];
        double nT11  = amp_norms["t1_1"];
        double nT21  = amp_norms["t2_1"];

        // print norms of cluster amplitudes
        Printf("\n\n   Norms of cluster amplitudes:");
        Printf("\n   ------------------------------");
        Printf("\n     T1: %15.12lf | %5.2f %%", sqrt(nT1), 100*nT1/total_norm);
        Printf("\n     T2: %15.12lf | %5.2f %%", sqrt(nT2), 100*nT2/total_norm);
        Printf("\n   T0,1: %15.12lf | %5.2f %%", sqrt(nT01), 100.0*nT01/total_norm);
        Printf("\n   T1,1: %15.12lf | %5.2f %%", sqrt(nT11), 100.0*nT11/total_norm);
        Printf("\n   T2,1: %15.12lf | %5.2f %%", sqrt(nT21), 100.0*nT21/total_norm);
        Printf("\n   ------------------------------");
        Printf("\n    Total: %15.12lf\n\n", sqrt(total_norm));
    }

    void QED_CCSD_21::print_iter_header() const {
        Printf("\n");
        Printf("    ==>  Begin %s iterations <==    \n", cc_type_.c_str());
        Printf("\n");
        Printf("%5s %16s %15s %15s  | %7s%7s ",  "Iter","energy","dE","|dT|","|dT1|","|dT2|");
        Printf("%7s","|dT0,1|");
        Printf("%7s","|dT1,1|");
        Printf("%7s","|dT2,1|");
        Printf("\n");
    }

    void QED_CCSD_21::print_iteration(size_t iter, double energy, double dele, double tnorm) const {
        Printf("%5i %17.12lf %15.12lf %15.12lf | %7.0e%7.0e",iter,energy,dele,tnorm,resid_norms_.at("t1"),resid_norms_.at("t2"));
        Printf("%7.0e", resid_norms_.at("t0_1"));
        Printf("%7.0e", resid_norms_.at("t1_1"));
        Printf("%7.0e", resid_norms_.at("t2_1"));
        Printf("\n");
    }

    TArrayMap QED_CCSD_21::effective_dipole() {
        // build effective dipole integrals

        // Get cavity information
        double w0 = cavity_frequency_[2];
        double coupling_factor_z = w0 * cavity_coupling_strength_[2];

        // build dipole integrals
        TArrayMap dp;
        dp["aa_oo"]("i, j") = coupling_factor_z * Dip_blks_["dz_aa_oo"]("i, j");
        dp["aa_ov"]("i, a") = coupling_factor_z * Dip_blks_["dz_aa_ov"]("i, a");
        dp["aa_vo"]("a, i") = coupling_factor_z * Dip_blks_["dz_aa_vo"]("a, i");
        dp["aa_vv"]("a, b") = coupling_factor_z * Dip_blks_["dz_aa_vv"]("a, b");

        dp["bb_oo"]("i, j") = coupling_factor_z * Dip_blks_["dz_bb_oo"]("i, j");
        dp["bb_ov"]("i, a") = coupling_factor_z * Dip_blks_["dz_bb_ov"]("i, a");
        dp["bb_vo"]("a, i") = coupling_factor_z * Dip_blks_["dz_bb_vo"]("a, i");
        dp["bb_vv"]("a, b") = coupling_factor_z * Dip_blks_["dz_bb_vv"]("a, b");

        return dp;
    }

    double QED_CCSD_21::build_residuals() {
        // t1 transformation of integrals
        if ( !has_t1_integrals_ ) transform_integrals(true);

        // zero out residual tensors
        for (auto &[name, resid] : residuals_) {
            zero_tiles(resid);
        }

        // initialize coherent state amplitudes to zero
        for (auto &[name, resid] : residuals_) {
            tmps_["c" + name] = resid.clone();
            zero_tiles(tmps_["c" + name]);
        }

        /// initialize doubles

        // extract ground state + hw amplitude
        double t0_1;
        foreach_inplace(amplitudes_["t0_1"], [&t0_1](auto &tile){
            for(auto &x : tile.range())
                t0_1 = tile[x];
        });

        // initialize scalars
        scalars_["energy"]  = 0.0;
        scalars_["cenergy"] = 0.0;
        scalars_["rt0_1"]   = 0.0;
        scalars_["crt0_1"]  = 0.0;
        scalars_["t0_1"]    = t0_1;
        scalars_["w0"]      = cavity_frequency_[2];

        // compute residuals
        resid_21_1();
        resid_21_2();
        resid_21_3();
        resid_21_4();
        world_.gop.fence();

        // extract scalars
        double &energy  = scalars_["energy"];
        double &cenergy = scalars_["cenergy"];
        double &rt0_1   = scalars_["rt0_1"];
        double &crt0_1  = scalars_["crt0_1"];

        // process residuals
        double coherent_scalar = cavity_frequency_[2] * cavity_coupling_strength_[2];
        if ( options_.get_bool("QED_USE_RELAXED_ORBITALS"))
             coherent_scalar *=    e_dip_z_;
        else coherent_scalar *= -nuc_dip_z_;

        // add coherent state basis terms
        energy += coherent_scalar * cenergy;
        rt0_1  += coherent_scalar * crt0_1;
        for (auto &[name, resid] : residuals_) {
            if (name == "t0_1") continue; // handled separately

            // add coherent state terms to residuals
            residuals_[name](idxs_[name]) += coherent_scalar * tmps_["c" + name](idxs_[name]);
        }
        residuals_["t0_1"] = makeTensor(world_, {1}, &rt0_1);

        // clear temporary arrays
        tmps_.clear();

        world_.gop.fence();
        return energy + enuc_ + average_electric_dipole_self_energy_;

    }


} // cc_cavity
