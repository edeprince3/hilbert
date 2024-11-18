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

#include "cc_cavity/include/qed_ccsd_21/eom_ee_qed_ccsd_21.h"
#include "cc_cavity/include/qed_ccsd_21/qed_ccsd_21.h"

void hilbert::EOM_EE_QED_CCSD_21::sigma_ee_21_7() {

    // Get cavity information
    double w0 = cc_wfn_->cavity_frequency_[2];
    double coupling_factor_z = w0 * cc_wfn_->cavity_coupling_strength_[2];

    // get effective dipole integrals
    TArrayMap dp = reinterpret_pointer_cast<QED_CCSD_21>(cc_wfn_)->effective_dipole();

    // get integrals
    TArrayMap &eri = cc_wfn_->V_blks_;
    TArrayMap &f = cc_wfn_->F_blks_;
    TArrayMap &Id = cc_wfn_->Id_blks_;


    /// extract amplitudes

    // extract 0-body amplitudes
    double t0_1;
    foreach_inplace(cc_wfn_->amplitudes_["t0_1"], [&t0_1](auto &tile){
        for(auto &x : tile.range())
            t0_1 = tile[x];
    });

    // extract 1-body amplitudes
    std::map<std::string, TA::TArrayD> t1_1 {
            {"aa", cc_wfn_->amplitudes_["t1_1_aa"]},
            {"bb", cc_wfn_->amplitudes_["t1_1_bb"]}
    };

    // extract 2-body amplitudes
    std::map<std::string, TA::TArrayD> t2 {
            {"aaaa", cc_wfn_->amplitudes_["t2_aaaa"]},
            {"abab", cc_wfn_->amplitudes_["t2_abab"]},
            {"bbbb", cc_wfn_->amplitudes_["t2_bbbb"]}
    };
    std::map<std::string, TA::TArrayD> t2_1 {
            {"aaaa", cc_wfn_->amplitudes_["t2_1_aaaa"]},
            {"abab", cc_wfn_->amplitudes_["t2_1_abab"]},
            {"bbbb", cc_wfn_->amplitudes_["t2_1_bbbb"]}
    };

    world_.gop.fence();

    /// unpack left/right operators

    TArrayMap r1, r2, r1_1, r2_1,
            l1, l2, l1_1, l2_1;

    TArrayD &r0 = evec_blks_["r0"];
    r1["aa"] = evec_blks_["r1_aa"];
    r1["bb"] = evec_blks_["r1_bb"];
    r2["aaaa"] = evec_blks_["r2_aaaa"];
    r2["abab"] = evec_blks_["r2_abab"];
    r2["bbbb"] = evec_blks_["r2_bbbb"];

    TArrayD &r0_1 = evec_blks_["r0_1"];
    r1_1["aa"] = evec_blks_["r1_1_aa"];
    r1_1["bb"] = evec_blks_["r1_1_bb"];
    r2_1["aaaa"] = evec_blks_["r2_1_aaaa"];
    r2_1["abab"] = evec_blks_["r2_1_abab"];
    r2_1["bbbb"] = evec_blks_["r2_1_bbbb"];

    TArrayD &l0 = evec_blks_["l0"];
    l1["aa"] = evec_blks_["l1_aa"];
    l1["bb"] = evec_blks_["l1_bb"];
    l2["aaaa"] = evec_blks_["l2_aaaa"];
    l2["abab"] = evec_blks_["l2_abab"];
    l2["bbbb"] = evec_blks_["l2_bbbb"];

    TArrayD &l0_1 = evec_blks_["l0_1"];
    l1_1["aa"] = evec_blks_["l1_1_aa"];
    l1_1["bb"] = evec_blks_["l1_1_bb"];
    l2_1["aaaa"] = evec_blks_["l2_1_aaaa"];
    l2_1["abab"] = evec_blks_["l2_1_abab"];
    l2_1["bbbb"] = evec_blks_["l2_1_bbbb"];

    /// extract the sigma vectors

    TArrayD &sigmar0 = sigvec_blks_["sigmar0"];
    TArrayD &sigmal0 = sigvec_blks_["sigmal0"];

    TArrayD &sigmar1_aa = sigvec_blks_["sigmar1_aa"];
    TArrayD &sigmar1_bb = sigvec_blks_["sigmar1_bb"];
    TArrayD &sigmal1_aa = sigvec_blks_["sigmal1_aa"];
    TArrayD &sigmal1_bb = sigvec_blks_["sigmal1_bb"];

    TArrayD &sigmar2_aaaa = sigvec_blks_["sigmar2_aaaa"];
    TArrayD &sigmar2_abab = sigvec_blks_["sigmar2_abab"];
    TArrayD &sigmar2_bbbb = sigvec_blks_["sigmar2_bbbb"];
    TArrayD &sigmal2_aaaa = sigvec_blks_["sigmal2_aaaa"];
    TArrayD &sigmal2_abab = sigvec_blks_["sigmal2_abab"];
    TArrayD &sigmal2_bbbb = sigvec_blks_["sigmal2_bbbb"];

    TArrayD &sigmar0_1      = sigvec_blks_["sigmar0_1"];
    TArrayD &sigmal0_1      = sigvec_blks_["sigmal0_1"];
    TArrayD &sigmar1_1_aa   = sigvec_blks_["sigmar1_1_aa"];
    TArrayD &sigmar1_1_bb   = sigvec_blks_["sigmar1_1_bb"];
    TArrayD &sigmal1_1_aa   = sigvec_blks_["sigmal1_1_aa"];
    TArrayD &sigmal1_1_bb   = sigvec_blks_["sigmal1_1_bb"];
    TArrayD &sigmar2_1_aaaa = sigvec_blks_["sigmar2_1_aaaa"];
    TArrayD &sigmar2_1_abab = sigvec_blks_["sigmar2_1_abab"];
    TArrayD &sigmar2_1_bbbb = sigvec_blks_["sigmar2_1_bbbb"];
    TArrayD &sigmal2_1_aaaa = sigvec_blks_["sigmal2_1_aaaa"];
    TArrayD &sigmal2_1_abab = sigvec_blks_["sigmal2_1_abab"];
    TArrayD &sigmal2_1_bbbb = sigvec_blks_["sigmal2_1_bbbb"];


    /// extract coherent state basis operators

    TArrayD &csigmar0        = tmps_["csigmar0"];
    TArrayD &csigmal0        = tmps_["csigmal0"];
    TArrayD &csigmar0_1      = tmps_["csigmar0_1"];
    TArrayD &csigmal0_1      = tmps_["csigmal0_1"];

    TArrayD &csigmar1_aa     = tmps_["csigmar1_aa"];
    TArrayD &csigmar1_bb     = tmps_["csigmar1_bb"];
    TArrayD &csigmar1_1_aa   = tmps_["csigmar1_1_aa"];
    TArrayD &csigmar1_1_bb   = tmps_["csigmar1_1_bb"];
    TArrayD &csigmal1_aa     = tmps_["csigmal1_aa"];
    TArrayD &csigmal1_bb     = tmps_["csigmal1_bb"];
    TArrayD &csigmal1_1_aa   = tmps_["csigmal1_1_aa"];
    TArrayD &csigmal1_1_bb   = tmps_["csigmal1_1_bb"];

    TArrayD &csigmar2_aaaa   = tmps_["csigmar2_aaaa"];
    TArrayD &csigmar2_abab   = tmps_["csigmar2_abab"];
    TArrayD &csigmar2_bbbb   = tmps_["csigmar2_bbbb"];
    TArrayD &csigmar2_1_aaaa = tmps_["csigmar2_1_aaaa"];
    TArrayD &csigmar2_1_abab = tmps_["csigmar2_1_abab"];
    TArrayD &csigmar2_1_bbbb = tmps_["csigmar2_1_bbbb"];
    TArrayD &csigmal2_aaaa   = tmps_["csigmal2_aaaa"];
    TArrayD &csigmal2_abab   = tmps_["csigmal2_abab"];
    TArrayD &csigmal2_bbbb   = tmps_["csigmal2_bbbb"];
    TArrayD &csigmal2_1_aaaa = tmps_["csigmal2_1_aaaa"];
    TArrayD &csigmal2_1_abab = tmps_["csigmal2_1_abab"];
    TArrayD &csigmal2_1_bbbb = tmps_["csigmal2_1_bbbb"];


    // sigmar2_abab  = -1.00 d-_aa(a,c) r1_aa(c,i) t1_1_bb(b,j)
    //              += -1.00 <l,k||c,j>_abab r1_aa(c,l) t2_abab(a,b,i,k)
    //              += -1.00 <l,k||c,j>_bbbb r1_bb(c,l) t2_abab(a,b,i,k)
    //              += -1.00 d-_bb(k,c) r1_bb(c,k) t2_1_abab(a,b,i,j)
    //              += -1.00 d-_aa(k,c) r1_aa(c,k) t2_1_abab(a,b,i,j)
    //              += +0.50 <l,k||c,d>_aaaa r2_aaaa(c,d,i,l) t2_abab(a,b,k,j)
    //              += -1.00 f_aa(k,c) r1_aa(c,i) t2_abab(a,b,k,j)
    //              += -0.50 <k,l||c,d>_abab r2_abab(c,d,i,l) t2_abab(a,b,k,j)
    //              += -0.50 <k,l||d,c>_abab r2_abab(d,c,i,l) t2_abab(a,b,k,j)
    //              += -1.00 d-_aa(k,c) r1_aa(a,i) t2_abab(c,b,k,j) t0_1
    //              += +1.00 d-_bb(k,c) r1_aa(a,i) t2_bbbb(c,b,j,k) t0_1
    //              += -1.00 <k,a||c,d>_aaaa r1_aa(d,k) t2_abab(c,b,i,j)
    //              += +1.00 <a,k||c,d>_abab r1_bb(d,k) t2_abab(c,b,i,j)
    //              += +1.00 <l,k||d,c>_abab r2_abab(d,b,l,j) t2_abab(a,c,i,k)
    //              += -1.00 d-_bb(b,c) r1_bb(c,j) t1_1_aa(a,i)
    //              += +0.50 <l,k||c,d>_aaaa r2_aaaa(d,a,l,k) t2_abab(c,b,i,j)
    //              += -0.50 <l,k||c,d>_abab r2_abab(a,d,l,k) t2_abab(c,b,i,j)
    //              += -0.50 <k,l||c,d>_abab r2_abab(a,d,k,l) t2_abab(c,b,i,j)
    //              += +1.00 d-_aa(k,c) r1_aa(c,i) t2_abab(a,b,k,j) t0_1
    //              += +1.00 d-_aa(k,c) r1_aa(a,k) t2_abab(c,b,i,j) t0_1
    //              += +1.00 d-_bb(k,c) r1_bb(c,j) t2_abab(a,b,i,k) t0_1
    //              += +1.00 d-_bb(k,c) r1_bb(b,k) t2_abab(a,c,i,j) t0_1
    //              += +1.00 d-_aa(k,i) r1_aa(a,k) t1_1_bb(b,j)
    //              += +0.25 <l,k||c,d>_abab r2_abab(c,d,i,j) t2_abab(a,b,l,k)
    //              += +0.25 <k,l||c,d>_abab r2_abab(c,d,i,j) t2_abab(a,b,k,l)
    //              += +0.25 <l,k||d,c>_abab r2_abab(d,c,i,j) t2_abab(a,b,l,k)
    //              += +0.25 <k,l||d,c>_abab r2_abab(d,c,i,j) t2_abab(a,b,k,l)
    //              += +1.00 d-_aa(k,c) r1_1_aa(c,i) t2_abab(a,b,k,j)
    //              += +1.00 d-_aa(k,c) r2_aaaa(c,a,i,k) t1_1_bb(b,j)
    //              += -1.00 d-_bb(k,c) r2_abab(a,c,i,k) t1_1_bb(b,j)
    //              += +1.00 <k,l||d,c>_abab r2_abab(d,b,i,l) t2_abab(a,c,k,j)
    //              += +0.50 <l,k||c,d>_bbbb r2_bbbb(d,b,l,k) t2_abab(a,c,i,j)
    //              += -0.50 <l,k||d,c>_abab r2_abab(d,b,l,k) t2_abab(a,c,i,j)
    //              += -0.50 <k,l||d,c>_abab r2_abab(d,b,k,l) t2_abab(a,c,i,j)
    //              += +1.00 d-_bb(k,j) r1_bb(b,k) t1_1_aa(a,i)
    //              += +1.00 d-_aa(k,c) r1_aa(c,i) t2_1_abab(a,b,k,j)
    //              += -1.00 d-_aa(k,c) r2_abab(c,b,k,j) t1_1_aa(a,i)
    //              += +1.00 d-_bb(k,c) r2_bbbb(c,b,j,k) t1_1_aa(a,i)
    //              += +1.00 d-_bb(k,c) r2_abab(a,c,i,j) t1_1_bb(b,k)
    //              += -1.00 <k,l||i,c>_abab r1_bb(c,l) t2_abab(a,b,k,j)
    //              += -1.00 <l,k||c,i>_aaaa r1_aa(c,l) t2_abab(a,b,k,j)
    //              += -1.00 f_bb(k,c) r1_bb(c,j) t2_abab(a,b,i,k)
    //              += +0.50 <l,k||c,d>_bbbb r2_bbbb(c,d,j,l) t2_abab(a,b,i,k)
    //              += -0.50 <l,k||c,d>_abab r2_abab(c,d,l,j) t2_abab(a,b,i,k)
    //              += -0.50 <l,k||d,c>_abab r2_abab(d,c,l,j) t2_abab(a,b,i,k)
    //              += +1.00 d-_aa(k,c) r2_abab(c,b,i,j) t1_1_aa(a,k)
    //              += +1.00 <k,b||d,c>_abab r1_aa(d,k) t2_abab(a,c,i,j)
    //              += -1.00 <k,b||c,d>_bbbb r1_bb(d,k) t2_abab(a,c,i,j)
    //              += +1.00 d-_bb(k,c) r1_bb(c,j) t2_1_abab(a,b,i,k)
    //              += +1.00 d-_bb(k,c) r1_1_bb(c,j) t2_abab(a,b,i,k)
    //              += -0.50 <k,b||c,d>_abab r1_aa(a,k) t2_abab(c,d,i,j)
    //              += -0.50 <k,b||d,c>_abab r1_aa(a,k) t2_abab(d,c,i,j)
    //              += -1.00 f_aa(k,c) r1_aa(a,k) t2_abab(c,b,i,j)
    //              += +1.00 <a,b||c,j>_abab r1_aa(c,i)
    //              += -0.50 <l,k||d,c>_abab r2_abab(d,b,i,j) t2_abab(a,c,l,k)
    //              += -0.50 <k,l||d,c>_abab r2_abab(d,b,i,j) t2_abab(a,c,k,l)
    //              += +1.00 f_bb(b,j) r1_aa(a,i)
    //              += -1.00 d-_bb(b,c) r2_abab(a,c,i,j) t0_1
    //              += -1.00 f_bb(k,c) r1_bb(b,k) t2_abab(a,c,i,j)
    //              += -0.50 <a,k||c,d>_abab r1_bb(b,k) t2_abab(c,d,i,j)
    //              += -0.50 <a,k||d,c>_abab r1_bb(b,k) t2_abab(d,c,i,j)
    //              += +1.00 d-_bb(k,c) r1_1_aa(a,i) t2_bbbb(c,b,j,k)
    //              += -1.00 d-_bb(b,c) r0_1 t2_abab(a,c,i,j)
    //              += +1.00 <l,k||i,c>_abab r1_aa(a,l) t2_bbbb(c,b,j,k)
    //              += +1.00 d-_aa(k,i) r1_bb(b,j) t1_1_aa(a,k)
    //              += +1.00 d-_aa(k,c) r1_bb(b,j) t2_1_aaaa(c,a,i,k)
    //              += -1.00 d-_aa(a,c) r1_bb(b,j) t1_1_aa(c,i)
    //              += -1.00 d-_bb(k,c) r1_bb(b,j) t2_1_abab(a,c,i,k)
    //              += +1.00 f_aa(a,i) r1_bb(b,j)
    //              += +1.00 <k,b||c,d>_bbbb r1_bb(d,j) t2_abab(a,c,i,k)
    //              += +0.50 <l,k||i,c>_abab r1_bb(c,j) t2_abab(a,b,l,k)
    //              += +0.50 <k,l||i,c>_abab r1_bb(c,j) t2_abab(a,b,k,l)
    //              += +0.50 <l,k||c,j>_abab r1_aa(c,i) t2_abab(a,b,l,k)
    //              += +0.50 <k,l||c,j>_abab r1_aa(c,i) t2_abab(a,b,k,l)
    //              += +1.00 <l,k||c,d>_bbbb r2_abab(a,d,i,l) t2_bbbb(c,b,j,k)
    //              += +1.00 <k,l||c,d>_abab r2_abab(a,d,i,l) t2_abab(c,b,k,j)
    //              += -1.00 d-_aa(k,c) r2_abab(a,b,i,j) t1_1_aa(c,k)
    //              += -1.00 d-_bb(b,j) r1_1_aa(a,i)
    //              += +1.00 <l,k||c,j>_abab r1_aa(a,l) t2_abab(c,b,i,k)
    //              += +1.00 d-_aa(k,i) r2_abab(a,b,k,j) t0_1
    //              += +1.00 <k,l||i,c>_abab r1_bb(b,l) t2_abab(a,c,k,j)
    //              += +1.00 <a,b||i,c>_abab r1_bb(c,j)
    //              += -0.50 <l,k||c,d>_aaaa r2_abab(a,b,l,j) t2_aaaa(c,d,i,k)
    //              += -1.00 d-_bb(k,k) r1_bb(b,j) t1_1_aa(a,i)
    //              += -1.00 d-_aa(k,k) r1_bb(b,j) t1_1_aa(a,i)
    //              += +1.00 f_aa(a,c) r2_abab(c,b,i,j)
    //              += -1.00 d-_bb(k,k) r2_1_abab(a,b,i,j)
    //              += -1.00 d-_bb(k,c) r1_1_bb(b,j) t2_abab(a,c,i,k)
    //              += +1.00 <l,k||d,c>_abab r2_aaaa(d,a,i,l) t2_bbbb(c,b,j,k)
    //              += +1.00 d-_bb(k,j) r2_abab(a,b,i,k) t0_1
    //              += +0.25 <l,k||c,d>_abab r2_abab(a,b,l,k) t2_abab(c,d,i,j)
    //              += +0.25 <l,k||d,c>_abab r2_abab(a,b,l,k) t2_abab(d,c,i,j)
    //              += +0.25 <k,l||c,d>_abab r2_abab(a,b,k,l) t2_abab(c,d,i,j)
    //              += +0.25 <k,l||d,c>_abab r2_abab(a,b,k,l) t2_abab(d,c,i,j)
    //              += -1.00 d-_aa(a,c) r0_1 t2_abab(c,b,i,j)
    //              += +1.00 d-_bb(k,j) r0_1 t2_abab(a,b,i,k)
    //              += -0.50 <l,k||c,d>_aaaa r2_abab(d,b,i,j) t2_aaaa(c,a,l,k)
    //              += +1.00 <k,a||c,i>_aaaa r2_abab(c,b,k,j)
    //              += +1.00 <l,k||c,j>_bbbb r1_bb(b,l) t2_abab(a,c,i,k)
    //              += +1.00 <l,k||c,d>_bbbb r2_bbbb(d,b,j,l) t2_abab(a,c,i,k)
    //              += +1.00 <l,k||c,d>_aaaa r2_abab(d,b,l,j) t2_aaaa(c,a,i,k)
    //              += +1.00 <k,l||c,j>_abab r1_bb(b,l) t2_aaaa(c,a,i,k)
    //              += -1.00 <a,k||i,c>_abab r2_bbbb(c,b,j,k)
    //              += +1.00 <l,k||c,i>_aaaa r1_aa(a,l) t2_abab(c,b,k,j)
    //              += +1.00 f_bb(b,c) r2_abab(a,c,i,j)
    //              += -1.00 <a,k||c,j>_abab r2_abab(c,b,i,k)
    //              += +1.00 d-_aa(k,c) r1_1_bb(b,j) t2_aaaa(c,a,i,k)
    //              += -1.00 d-_bb(k,c) r2_abab(a,b,i,j) t1_1_bb(c,k)
    //              += -1.00 d-_aa(k,k) r2_1_abab(a,b,i,j)
    //              += -0.50 <l,k||c,d>_abab r2_abab(a,d,i,j) t2_abab(c,b,l,k)
    //              += -0.50 <k,l||c,d>_abab r2_abab(a,d,i,j) t2_abab(c,b,k,l)
    //              += -1.00 d-_bb(k,c) r1_bb(b,j) t2_abab(a,c,i,k) t0_1
    //              += -0.50 <k,l||c,d>_abab r2_abab(a,b,i,l) t2_abab(c,d,k,j)
    //              += -0.50 <k,l||d,c>_abab r2_abab(a,b,i,l) t2_abab(d,c,k,j)
    //              += -1.00 d-_bb(b,c) r1_aa(a,i) t1_1_bb(c,j)
    //              += +1.00 d-_bb(k,j) r1_aa(a,i) t1_1_bb(b,k)
    //              += -1.00 d-_aa(k,c) r1_aa(a,i) t2_1_abab(c,b,k,j)
    //              += +1.00 d-_bb(k,c) r1_aa(a,i) t2_1_bbbb(c,b,j,k)
    //              += +0.50 <l,k||i,j>_abab r2_abab(a,b,l,k)
    //              += +0.50 <k,l||i,j>_abab r2_abab(a,b,k,l)
    //              += -1.00 d-_aa(k,k) r2_abab(a,b,i,j) t0_1
    //              += +1.00 f_aa(j,j) l1_1_bb(i,a)
    //              += -0.50 <k,j||k,j>_abab l1_1_bb(i,a)
    //              += -0.50 <j,k||j,k>_abab l1_1_bb(i,a)
    //              += +0.25 <k,j||b,c>_aaaa t2_aaaa(b,c,k,j) l1_1_bb(i,a)
    //              += -0.50 <k,j||k,j>_bbbb l1_1_bb(i,a)
    //              += -0.50 <k,j||k,j>_aaaa l1_1_bb(i,a)
    //              += +0.25 <k,j||b,c>_bbbb t2_bbbb(b,c,k,j) l1_1_bb(i,a)
    //              += +1.00 f_bb(j,j) l1_1_bb(i,a)
    //              += -1.00 d-_bb(j,j) t0_1 l1_1_bb(i,a)
    //              += +0.25 <k,j||b,c>_abab t2_abab(b,c,k,j) l1_1_bb(i,a)
    //              += +0.25 <j,k||b,c>_abab t2_abab(b,c,j,k) l1_1_bb(i,a)
    //              += +0.25 <k,j||c,b>_abab t2_abab(c,b,k,j) l1_1_bb(i,a)
    //              += +0.25 <j,k||c,b>_abab t2_abab(c,b,j,k) l1_1_bb(i,a)
    //              += -1.00 d-_aa(a,i) r1_1_bb(b,j)
    //              += -0.50 <l,k||c,d>_abab r2_abab(a,b,l,j) t2_abab(c,d,i,k)
    //              += -0.50 <l,k||d,c>_abab r2_abab(a,b,l,j) t2_abab(d,c,i,k)
    //              += -1.00 <k,b||d,c>_abab r1_aa(d,i) t2_abab(a,c,k,j)
    //              += +1.00 <l,k||c,d>_aaaa r2_aaaa(d,a,i,l) t2_abab(c,b,k,j)
    //              += +0.50 <a,b||c,d>_abab r2_abab(c,d,i,j)
    //              += +0.50 <a,b||d,c>_abab r2_abab(d,c,i,j)
    //              += +1.00 d-_bb(k,j) r2_1_abab(a,b,i,k)
    //              += -1.00 <k,b||i,c>_abab r2_abab(a,c,k,j)
    //              += -1.00 <k,b||c,j>_abab r2_aaaa(c,a,i,k)
    //              += -0.50 <k,a||c,d>_aaaa r1_bb(b,j) t2_aaaa(c,d,i,k)
    //              += -0.50 <l,k||i,c>_abab r1_bb(b,j) t2_abab(a,c,l,k)
    //              += -0.50 <k,l||i,c>_abab r1_bb(b,j) t2_abab(a,c,k,l)
    //              += -0.50 <l,k||c,i>_aaaa r1_bb(b,j) t2_aaaa(c,a,l,k)
    //              += -1.00 f_aa(k,c) r1_bb(b,j) t2_aaaa(c,a,i,k)
    //              += +0.50 <a,k||c,d>_abab r1_bb(b,j) t2_abab(c,d,i,k)
    //              += +0.50 <a,k||d,c>_abab r1_bb(b,j) t2_abab(d,c,i,k)
    //              += -1.00 d-_aa(a,i) r1_bb(b,j) t0_1
    //              += +1.00 f_bb(k,c) r1_bb(b,j) t2_abab(a,c,i,k)
    //              += -1.00 d-_aa(a,c) r2_abab(c,b,i,j) t0_1
    //              += -1.00 <a,k||c,d>_abab r1_bb(d,j) t2_abab(c,b,i,k)
    //              += -1.00 <a,k||i,j>_abab r1_bb(b,k)
    //              += -1.00 <a,k||d,c>_abab r1_aa(d,i) t2_bbbb(c,b,j,k)
    //              += -1.00 f_aa(k,i) r2_abab(a,b,k,j)
    //              += -0.50 <l,k||c,d>_bbbb r2_abab(a,d,i,j) t2_bbbb(c,b,l,k)
    //              += +1.00 <k,b||c,j>_bbbb r2_abab(a,c,i,k)
    //              += -1.00 f_bb(k,j) r2_abab(a,b,i,k)
    //              += -1.00 d-_aa(k,k) r1_aa(a,i) t1_1_bb(b,j)
    //              += -1.00 d-_bb(k,k) r1_aa(a,i) t1_1_bb(b,j)
    //              += -0.50 <l,k||c,d>_bbbb r2_abab(a,b,i,l) t2_bbbb(c,d,j,k)
    //              += +1.00 <l,k||c,d>_abab r2_abab(a,d,l,j) t2_abab(c,b,i,k)
    //              += -1.00 d-_aa(k,c) r1_1_aa(a,i) t2_abab(c,b,k,j)
    //              += -1.00 d-_bb(b,c) r2_1_abab(a,c,i,j)
    //              += -1.00 d-_aa(a,c) r2_1_abab(c,b,i,j)
    //              += +1.00 <k,a||c,d>_aaaa r1_aa(d,i) t2_abab(c,b,k,j)
    //              += -1.00 <k,b||i,j>_abab r1_aa(a,k)
    //              += +1.00 d-_aa(k,i) r2_1_abab(a,b,k,j)
    //              += +1.00 <k,l||c,d>_abab r2_bbbb(d,b,j,l) t2_aaaa(c,a,i,k)
    //              += -1.00 <k,b||c,d>_abab r1_bb(d,j) t2_aaaa(c,a,i,k)
    //              += -1.00 f_bb(k,c) r1_aa(a,i) t2_bbbb(c,b,j,k)
    //              += -1.00 d-_bb(b,j) r1_aa(a,i) t0_1
    //              += -0.50 <l,k||c,j>_abab r1_aa(a,i) t2_abab(c,b,l,k)
    //              += -0.50 <k,l||c,j>_abab r1_aa(a,i) t2_abab(c,b,k,l)
    //              += -0.50 <l,k||c,j>_bbbb r1_aa(a,i) t2_bbbb(c,b,l,k)
    //              += +1.00 f_aa(k,c) r1_aa(a,i) t2_abab(c,b,k,j)
    //              += +0.50 <k,b||c,d>_abab r1_aa(a,i) t2_abab(c,d,k,j)
    //              += +0.50 <k,b||d,c>_abab r1_aa(a,i) t2_abab(d,c,k,j)
    //              += -0.50 <k,b||c,d>_bbbb r1_aa(a,i) t2_bbbb(c,d,j,k)
    //              += +1.00 d-_aa(k,c) r1_bb(b,j) t2_aaaa(c,a,i,k) t0_1
    sigmar2_abab("R,a,b,i,j")  = -1.00 * tmps_["280_bbaa_Lvovo"]("R,b,j,a,i");
    tmps_["280_bbaa_Lvovo"].~TArrayD();

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["285_abab_Lvvoo"]("R,a,b,i,j")  = r0_1("R") * t2_1["abab"]("a,b,i,j");

    // csigmar2_1_abab += +1.00 r0_1 t2_1_abab(a,b,i,j)
    csigmar2_1_abab("R,a,b,i,j") += tmps_["285_abab_Lvvoo"]("R,a,b,i,j");

    // flops: o2v2L1  = o3v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["289_aabb_Lovvo"]("R,i,a,b,j")  = r1["aa"]("R,a,l") * reused_["190_baba_vooo"]("b,i,j,l");
    tmps_["289_aabb_Lovvo"]("R,i,a,b,j") += r0_1("R") * reused_["72_aabb_ovvo"]("i,a,b,j");
    tmps_["289_aabb_Lovvo"]("R,i,a,b,j") += r2["abab"]("R,a,b,l,j") * reused_["120_aa_oo"]("l,i");
    tmps_["289_aabb_Lovvo"]("R,i,a,b,j") += r1_1["aa"]("R,a,l") * reused_["83_baba_vooo"]("b,i,j,l");
    tmps_["289_aabb_Lovvo"]("R,i,a,b,j") += r2["abab"]("R,a,b,i,k") * reused_["38_bb_oo"]("k,j");
    tmps_["289_aabb_Lovvo"]("R,i,a,b,j") += r1_1["bb"]("R,b,k") * reused_["77_aabb_vooo"]("a,i,j,k");
    tmps_["289_aabb_Lovvo"]("R,i,a,b,j") += r1["bb"]("R,b,k") * reused_["208_aabb_vooo"]("a,i,j,k");

    // sigmar2_abab += +1.00 d-_aa(k,i) r0_1 t2_abab(a,b,k,j)
    //              += +1.00 d-_aa(k,c) r1_aa(a,k) t2_1_abab(c,b,i,j)
    //              += +1.00 d-_aa(k,c) r2_abab(a,b,k,j) t1_1_aa(c,i)
    //              += +1.00 d-_aa(k,c) r1_1_aa(a,k) t2_abab(c,b,i,j)
    //              += +1.00 d-_bb(k,c) r2_abab(a,b,i,k) t1_1_bb(c,j)
    //              += +1.00 d-_bb(k,c) r1_1_bb(b,k) t2_abab(a,c,i,j)
    //              += +1.00 d-_bb(k,c) r1_bb(b,k) t2_1_abab(a,c,i,j)
    sigmar2_abab("R,a,b,i,j") += tmps_["289_aabb_Lovvo"]("R,i,a,b,j");

    // flops: o1v1L1  = o1v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["292_bb_Lvo"]("R,a,i")  = dp["bb_vv"]("a,b") * r1_1["bb"]("R,b,i");

    // sigmar1_bb += -1.00 d-_bb(a,b) r1_1_bb(b,i)
    sigmar1_bb("R,a,i") -= tmps_["292_bb_Lvo"]("R,a,i");

    // flops: o1v1L1  = o2v1L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["293_bb_Lov"]("R,i,a")  = dp["bb_oo"]("j,i") * r1_1["bb"]("R,a,j");

    // sigmar1_bb += +1.00 d-_bb(j,i) r1_1_bb(a,j)
    sigmar1_bb("R,a,i") += tmps_["293_bb_Lov"]("R,i,a");

    // flops: o1v1L1  = o2v2L1 o2v2L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1
    tmps_["294_bb_Lvo"]("R,a,i")  = -1.00 * dp["bb_ov"]("k,c") * r2_1["bbbb"]("R,c,a,i,k");
    tmps_["294_bb_Lvo"]("R,a,i") += dp["aa_ov"]("j,b") * r2_1["abab"]("R,b,a,j,i");

    // sigmar1_bb += -1.00 d-_aa(j,b) r2_1_abab(b,a,j,i)
    //            += +1.00 d-_bb(j,b) r2_1_bbbb(b,a,i,j)
    sigmar1_bb("R,a,i") -= tmps_["294_bb_Lvo"]("R,a,i");

    // flops: o1v1L1  = o2v1L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["295_bb_Lvo"]("R,a,i")  = t1_1["bb"]("a,j") * tmps_["73_bb_Loo"]("R,i,j");
    tmps_["73_bb_Loo"].~TArrayD();

    // sigmar1_bb += +1.00 d-_bb(j,b) r1_bb(b,i) t1_1_bb(a,j)
    sigmar1_bb("R,a,i") += tmps_["295_bb_Lvo"]("R,a,i");

    // flops: o3v1L1  = o3v2L1 o3v2L1 o3v2L1 o3v1L1 o4v2L1 o3v1L1 o3v2L1 o3v1L1 o3v2L1 o3v1L1 o4v2L1 o3v1L1 o4v2L1 o3v1L1 o3v3L1 o3v1L1 o3v2L1 o3v1L1 o3v2L1 o3v1L1 o3v1L1
    //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1
    tmps_["291_aabb_Lvooo"]("R,a,i,j,k")  = t2["abab"]("a,b,i,j") * tmps_["68_bb_Lov"]("R,k,b");
    tmps_["291_aabb_Lvooo"]("R,a,i,j,k") -= eri["abba_vovo"]("a,k,b,i") * r1["bb"]("R,b,j");
    tmps_["291_aabb_Lvooo"]("R,a,i,j,k") += eri["abab_vovo"]("a,k,e,j") * r1["aa"]("R,e,i");
    tmps_["291_aabb_Lvooo"]("R,a,i,j,k") -= eri["abab_oovo"]("l,k,e,j") * r2["aaaa"]("R,e,a,i,l");
    tmps_["291_aabb_Lvooo"]("R,a,i,j,k") -= 2.00 * r2_1["abab"]("R,a,b,i,j") * dp["bb_ov"]("k,b");
    tmps_["291_aabb_Lvooo"]("R,a,i,j,k") += r2["abab"]("R,a,b,i,j") * f["bb_ov"]("k,b");
    tmps_["291_aabb_Lvooo"]("R,a,i,j,k") -= eri["bbbb_oovo"]("k,m,b,j") * r2["abab"]("R,a,b,i,m");
    tmps_["291_aabb_Lvooo"]("R,a,i,j,k") += eri["abba_oovo"]("l,k,b,i") * r2["abab"]("R,a,b,l,j");
    tmps_["291_aabb_Lvooo"]("R,a,i,j,k") += eri["abab_vovv"]("a,k,e,d") * r2["abab"]("R,e,d,i,j");
    tmps_["291_aabb_Lvooo"]("R,a,i,j,k") -= r2["abab"]("R,a,b,i,j") * reused_["12_bb_ov"]("k,b");
    tmps_["291_aabb_Lvooo"]("R,a,i,j,k") += t2["abab"]("a,b,i,j") * tmps_["151_bb_Lov"]("R,k,b");
    tmps_["151_bb_Lov"].~TArrayD();
    tmps_["68_bb_Lov"].~TArrayD();

    // flops: o3v1L1  = o3v2L1
    //  mems: o3v1L1  = o3v1L1
    tmps_["290_baba_Lvooo"]("R,a,i,j,k")  = t2["abab"]("b,a,i,j") * tmps_["144_aa_Lov"]("R,k,b");
    tmps_["144_aa_Lov"].~TArrayD();

    // flops: o3v1L1  = o3v2L1 o3v2L1 o3v2L1 o3v1L1 o3v1L1
    //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1
    tmps_["288_aabb_Lvooo"]("R,a,i,j,k")  = -1.00 * reused_["56_aabb_voov"]("a,i,j,b") * r1["bb"]("R,b,k");
    tmps_["288_aabb_Lvooo"]("R,a,i,j,k") += reused_["6_aabb_voov"]("a,i,j,b") * r1["bb"]("R,b,k");
    tmps_["288_aabb_Lvooo"]("R,a,i,j,k") += reused_["191_abba_voov"]("a,k,j,c") * r1["aa"]("R,c,i");

    // flops: o3v1L1  = o3v2L1 o3v2L1 o3v2L1 o3v1L1 o3v1L1
    //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1
    tmps_["287_baab_Lvooo"]("R,a,i,j,k")  = reused_["1_bbaa_voov"]("a,k,j,b") * r1["aa"]("R,b,i");
    tmps_["287_baab_Lvooo"]("R,a,i,j,k") += reused_["78_baab_voov"]("a,i,j,c") * r1["bb"]("R,c,k");
    tmps_["287_baab_Lvooo"]("R,a,i,j,k") -= reused_["27_bbaa_voov"]("a,k,j,b") * r1["aa"]("R,b,i");

    // flops: o3v1L1  = o3v2L1 o3v3L1 o3v2L1 o3v1L1 o4v2L1 o3v1L1 o3v2L1 o3v1L1 o3v2L1 o3v1L1 o4v2L1 o3v1L1 o4v1L1 o3v1L1 o3v1L1 o4v2L1 o3v1L1 o4v1L1 o3v1L1 o3v2L1 o3v1L1
    //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1
    tmps_["286_baab_Lvooo"]("R,a,i,j,k")  = -1.00 * eri["baba_vovo"]("a,i,c,j") * r1["bb"]("R,c,k");
    tmps_["286_baab_Lvooo"]("R,a,i,j,k") += eri["baab_vovv"]("a,i,b,d") * r2["abab"]("R,b,d,j,k");
    tmps_["286_baab_Lvooo"]("R,a,i,j,k") += r2["abab"]("R,b,a,j,k") * reused_["13_aa_ov"]("i,b");
    tmps_["286_baab_Lvooo"]("R,a,i,j,k") += eri["abab_oovo"]("i,l,b,k") * r2["abab"]("R,b,a,j,l");
    tmps_["286_baab_Lvooo"]("R,a,i,j,k") -= r2["abab"]("R,b,a,j,k") * f["aa_ov"]("i,b");
    tmps_["286_baab_Lvooo"]("R,a,i,j,k") += eri["baab_vovo"]("a,i,b,k") * r1["aa"]("R,b,j");
    tmps_["286_baab_Lvooo"]("R,a,i,j,k") -= eri["abba_oovo"]("i,l,c,j") * r2["bbbb"]("R,c,a,k,l");
    tmps_["286_baab_Lvooo"]("R,a,i,j,k") += r1["bb"]("R,a,l") * reused_["80_abab_oooo"]("i,l,j,k");
    tmps_["286_baab_Lvooo"]("R,a,i,j,k") += eri["aaaa_oovo"]("i,m,b,j") * r2["abab"]("R,b,a,m,k");
    tmps_["286_baab_Lvooo"]("R,a,i,j,k") += eri["abab_oooo"]("i,l,j,k") * r1["bb"]("R,a,l");
    tmps_["286_baab_Lvooo"]("R,a,i,j,k") += 2.00 * dp["aa_ov"]("i,b") * r2_1["abab"]("R,b,a,j,k");

    // flops: o2v2L1  = o3v3L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["284_bbbb_Lvoov"]("R,a,i,j,b")  = r2_1["abab"]("R,c,a,k,i") * eri["abab_oovv"]("k,j,c,b");

    // flops: o2v2L1  = o3v3L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["283_baab_Lvoov"]("R,a,i,j,b")  = r2_1["abab"]("R,c,a,i,k") * eri["abab_oovv"]("j,k,c,b");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["281_bb_Lvo"]("R,a,i")  = -1.00 * dp["bb_ov"]("j,b") * r2_1["bbbb"]("R,b,a,i,j");
    tmps_["281_bb_Lvo"].~TArrayD();

    // flops: o4v0L1  = o4v2L1
    //  mems: o4v0L1  = o4v0L1
    tmps_["282_abab_Loooo"]("R,i,j,k,l")  = r2_1["abab"]("R,a,b,i,j") * eri["abab_oovv"]("k,l,a,b");

    // flops: o2v2L1  = o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o3v2L1 o2v2L1 o2v4L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o0v0L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o4v2L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o0v0L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o4v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o0v0L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o0v0L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o0v0L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o0v0L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j")  = -1.00 * t1_1["aa"]("b,j") * tmps_["107_bb_Lvo"]("R,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= tmps_["147_aa_Lvv"]("R,b,c") * t2["abab"]("c,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= t2["abab"]("b,e,k,i") * tmps_["283_baab_Lvoov"]("R,a,j,k,e");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += t2_1["abab"]("b,e,j,i") * tmps_["65_bb_Lvv"]("R,a,e");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= 2.00 * t2_1["abab"]("b,a,k,i") * tmps_["34_aa_Loo"]("R,j,k");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= t2["abab"]("b,a,k,i") * tmps_["129_aa_Loo"]("R,k,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += tmps_["162_aa_Lvo"]("R,b,j") * t1_1["bb"]("a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += t2["abab"]("b,e,j,i") * tmps_["154_bb_Lvv"]("R,a,e");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += scalars_["1"] * tmps_["285_abab_Lvvoo"]("R,b,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= t1_1["aa"]("b,k") * tmps_["286_baab_Lvooo"]("R,a,k,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += dp["aa_vo"]("b,j") * tmps_["108_bb_Lvo"]("R,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["58_aa_vo"]("b,j") * tmps_["108_bb_Lvo"]("R,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= t2["abab"]("b,e,j,m") * tmps_["284_bbbb_Lvoov"]("R,a,i,m,e");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= t2_1["abab"]("b,a,l,m") * tmps_["279_abab_Loooo"]("R,j,i,l,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += t1_1["aa"]("b,j") * tmps_["294_bb_Lvo"]("R,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += t2_1["abab"]("b,a,j,m") * tmps_["54_bb_Loo"]("R,m,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += tmps_["146_aa_Lvv"]("R,b,c") * t2["abab"]("c,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= t0_1 * tmps_["289_aabb_Lovvo"]("R,j,b,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r1["aa"]("R,b,j") * reused_["200_bb_vo"]("a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= t1_1["aa"]("b,l") * tmps_["287_baab_Lvooo"]("R,a,j,l,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= t2["abab"]("b,e,j,i") * tmps_["153_bb_Lvv"]("R,a,e");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= t2_1["abab"]("b,e,k,i") * tmps_["278_baab_Lvoov"]("R,a,j,k,e");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= t2_1["abab"]("b,a,k,i") * tmps_["37_aa_Loo"]("R,k,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += 0.50 * t2_1["abab"]("b,a,k,i") * tmps_["36_aa_Loo"]("R,j,k");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= 2.00 * t2_1["abab"]("b,a,j,m") * tmps_["49_bb_Loo"]("R,i,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += 0.50 * tmps_["135_aa_Lvv"]("R,b,c") * t2_1["abab"]("c,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= tmps_["134_aa_Lvv"]("R,b,c") * t2_1["abab"]("c,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= t1_1["bb"]("a,i") * tmps_["165_aa_Lov"]("R,j,b");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r0_1("R") * reused_["240_abab_vvoo"]("b,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["55_aaaa_voov"]("b,j,l,d") * r2_1["abab"]("R,d,a,l,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r2["abab"]("R,b,a,l,i") * reused_["125_aa_oo"]("l,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= eri["abab_vvvv"]("b,a,c,f") * r2_1["abab"]("R,c,f,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= 0.50 * reused_["87_aa_vo"]("b,j") * r1_1["bb"]("R,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["58_aa_vo"]("b,j") * r1["bb"]("R,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += dp["bb_vo"]("a,i") * r1["aa"]("R,b,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += 2.00 * r1_1["aa"]("R,b,j") * reused_["88_bb_vo"]("a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= reused_["221_abab_vooo"]("b,m,j,i") * r1["bb"]("R,a,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r2["aaaa"]("R,d,b,j,k") * reused_["105_baab_vovo"]("a,k,d,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r2_1["abab"]("R,b,f,j,i") * reused_["3_bb_vv"]("f,a");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r1_1["aa"]("R,b,k") * reused_["175_baab_vooo"]("a,k,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r1_1["aa"]("R,b,j") * reused_["84_bb_vo"]("a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r2_1["abab"]("R,b,a,j,n") * reused_["30_bb_oo"]("n,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["242_abab_vovo"]("b,m,d,i") * r2["abab"]("R,d,a,j,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= f["aa_vv"]("b,c") * r2_1["abab"]("R,c,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += 2.00 * scalars_["4"] * r2_1["abab"]("R,b,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r1_1["aa"]("R,b,l") * reused_["233_aabb_oovo"]("l,j,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r1["aa"]("R,b,l") * reused_["198_baab_vooo"]("a,j,l,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["122_aa_vv"]("b,d") * r2["abab"]("R,d,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["227_aabb_vooo"]("b,j,i,n") * r1["bb"]("R,a,n");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += scalars_["2"] * r0_1("R") * t2_1["abab"]("b,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += 0.50 * r0_1("R") * reused_["237_baab_vvoo"]("a,b,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["152_aa_vo"]("b,j") * r1_1["bb"]("R,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += scalars_["5"] * r2_1["abab"]("R,b,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r2["abab"]("R,b,a,l,i") * reused_["124_aa_oo"]("l,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["230_aabb_vvvo"]("b,d,a,i") * r1_1["aa"]("R,d,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= 2.00 * r1_1["bb"]("R,a,i") * reused_["89_aa_ov"]("j,b");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= eri["baab_vooo"]("a,k,j,i") * r1_1["aa"]("R,b,k");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["57_aa_vv"]("d,b") * r2_1["abab"]("R,d,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r2["abab"]("R,b,a,k,i") * reused_["121_aa_oo"]("j,k");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= reused_["114_aabb_voov"]("b,j,n,f") * r2["bbbb"]("R,f,a,i,n");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= 2.00 * r1_1["aa"]("R,b,k") * reused_["190_baba_vooo"]("a,j,i,k");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["167_aabb_vvvo"]("b,d,a,i") * r1_1["aa"]("R,d,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r2["abab"]("R,b,f,j,i") * reused_["108_bb_vv"]("a,f");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= w0 * r2_1["abab"]("R,b,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= eri["abab_vvoo"]("b,a,j,i") * r0_1("R");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r2["aaaa"]("R,d,b,j,l") * reused_["101_bbaa_voov"]("a,i,l,d");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r2_1["abab"]("R,b,f,l,i") * reused_["78_baab_voov"]("a,j,l,f");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["56_aabb_voov"]("b,j,n,f") * r2_1["bbbb"]("R,f,a,i,n");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r2_1["abab"]("R,b,f,j,n") * reused_["28_bbbb_voov"]("a,i,n,f");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r1["aa"]("R,b,l") * reused_["238_aabb_oovo"]("l,j,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += eri["baba_vovo"]("a,k,e,j") * r2_1["abab"]("R,b,e,k,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= reused_["6_aabb_voov"]("b,j,n,f") * r2_1["bbbb"]("R,f,a,i,n");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r0_1("R") * reused_["95_aabb_voov"]("b,j,i,a");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r0_1("R") * reused_["93_aabb_vovo"]("b,j,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r2["abab"]("R,b,f,k,i") * reused_["243_baba_vovo"]("a,k,f,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= reused_["54_aa_vo"]("b,j") * r1["bb"]("R,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r1_1["aa"]("R,b,j") * reused_["91_bb_vo"]("a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r2_1["abab"]("R,b,e,j,i") * reused_["34_bb_vv"]("a,e");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["224_abba_vooo"]("b,m,i,j") * r1["bb"]("R,a,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= f["bb_vo"]("a,i") * r1_1["aa"]("R,b,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= f["bb_vv"]("a,e") * r2_1["abab"]("R,b,e,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= reused_["117_aaaa_vovo"]("b,k,d,j") * r2["abab"]("R,d,a,k,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += f["aa_oo"]("k,j") * r2_1["abab"]("R,b,a,k,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["209_abba_vvvo"]("b,f,a,j") * r1_1["bb"]("R,f,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["223_abba_vvvo"]("b,f,a,j") * r1["bb"]("R,f,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r1["aa"]("R,b,j") * reused_["32_bb_vo"]("a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += eri["abba_vvvo"]("b,a,e,j") * r1_1["bb"]("R,e,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r2["abab"]("R,b,f,l,i") * reused_["241_baab_voov"]("a,j,l,f");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r2["abab"]("R,b,a,j,n") * reused_["48_bb_oo"]("n,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= eri["abab_oooo"]("l,m,j,i") * r2_1["abab"]("R,b,a,l,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= 0.50 * r2["abab"]("R,b,f,j,i") * reused_["107_bb_vv"]("a,f");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r2["abab"]("R,b,f,j,n") * reused_["104_bbbb_voov"]("a,i,n,f");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r1_1["aa"]("R,b,l") * reused_["173_bbaa_vooo"]("a,i,l,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= dp["aa_oo"]("k,j") * r2["abab"]("R,b,a,k,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r1["aa"]("R,b,j") * reused_["33_bb_vo"]("a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= 2.00 * r2_1["abab"]("R,b,a,j,m") * reused_["38_bb_oo"]("m,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= 0.50 * r2_1["abab"]("R,b,a,l,i") * reused_["15_aa_oo"]("j,l");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["90_aa_vo"]("b,j") * r1_1["bb"]("R,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r2["abab"]("R,b,a,l,m") * reused_["245_abab_oooo"]("j,i,l,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += t0_1 * r0_1("R") * reused_["73_abab_vvoo"]("b,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["59_aa_vv"]("b,c") * r2_1["abab"]("R,c,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["115_aabb_voov"]("b,j,n,f") * r2["bbbb"]("R,f,a,i,n");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["225_aabb_vooo"]("b,j,i,m") * r1_1["bb"]("R,a,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= eri["baab_vovo"]("a,k,c,i") * r2_1["aaaa"]("R,c,b,j,k");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r2_1["aaaa"]("R,d,b,j,l") * reused_["27_bbaa_voov"]("a,i,l,d");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r2["abab"]("R,b,f,j,m") * reused_["106_bbbb_vovo"]("a,m,f,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r2_1["abab"]("R,b,a,l,m") * reused_["80_abab_oooo"]("l,m,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["239_aabb_vooo"]("b,j,n,i") * r1["bb"]("R,a,n");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r1_1["aa"]("R,c,j") * reused_["169_abab_vovv"]("c,i,b,a");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["113_aaaa_voov"]("b,j,l,d") * r2["abab"]("R,d,a,l,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += 2.00 * scalars_["3"] * r2_1["abab"]("R,b,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += dp["bb_vv"]("a,e") * r2["abab"]("R,b,e,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["118_abba_vovo"]("b,m,f,j") * r2["bbbb"]("R,f,a,i,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += dp["aa_vo"]("b,j") * r1["bb"]("R,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r1["aa"]("R,b,l") * reused_["197_baba_vooo"]("a,j,i,l");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= eri["abab_vvvo"]("b,a,c,i") * r1_1["aa"]("R,c,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r2["abab"]("R,b,f,j,n") * reused_["102_bbbb_voov"]("a,i,n,f");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r2_1["aaaa"]("R,d,b,j,l") * reused_["1_bbaa_voov"]("a,i,l,d");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += eri["bbbb_vovo"]("a,m,e,i") * r2_1["abab"]("R,b,e,j,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r1_1["aa"]("R,b,l") * reused_["171_abba_oovo"]("l,i,a,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r2_1["abab"]("R,b,a,l,i") * reused_["16_aa_oo"]("l,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r2_1["abab"]("R,b,a,k,i") * reused_["17_aa_oo"]("k,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r2_1["abab"]("R,b,a,j,m") * reused_["35_bb_oo"]("m,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r0_1("R") * reused_["98_abba_vvoo"]("b,a,i,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= 0.50 * r2["abab"]("R,b,a,j,n") * reused_["47_bb_oo"]("i,n");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= eri["abba_vovo"]("b,m,e,j") * r2_1["bbbb"]("R,e,a,i,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r1["aa"]("R,c,j") * reused_["193_abab_vovv"]("c,i,b,a");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r1["bb"]("R,a,n") * reused_["226_baab_oovo"]("n,j,b,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r1["aa"]("R,b,k") * reused_["194_baab_vooo"]("a,k,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r2["abab"]("R,b,a,l,m") * reused_["244_abab_oooo"]("l,m,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= 0.50 * r1["aa"]("R,b,j") * reused_["199_bb_ov"]("i,a");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["235_aabb_vvvo"]("b,d,a,i") * r1["aa"]("R,d,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += eri["abab_vooo"]("b,m,j,i") * r1_1["bb"]("R,a,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["215_baab_vovv"]("e,j,b,a") * r1_1["bb"]("R,e,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["236_bbaa_vvvo"]("a,f,b,j") * r1["bb"]("R,f,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r1["aa"]("R,b,k") * reused_["195_baba_vooo"]("a,j,i,k");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= 0.50 * r2_1["abab"]("R,b,a,j,n") * reused_["31_bb_oo"]("i,n");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r2["aaaa"]("R,d,b,j,l") * reused_["103_bbaa_voov"]("a,i,l,d");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r2["abab"]("R,b,a,j,m") * reused_["40_bb_oo"]("m,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r1_1["bb"]("R,a,n") * reused_["232_bbaa_oovo"]("n,i,b,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= scalars_["6"] * r2["abab"]("R,b,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= f["aa_vo"]("b,j") * r1_1["bb"]("R,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= reused_["211_bbaa_vvvo"]("a,f,b,j") * r1["bb"]("R,f,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= 0.50 * r2_1["abab"]("R,b,f,j,i") * reused_["29_bb_vv"]("a,f");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= reused_["168_abba_vovv"]("b,i,a,d") * r1_1["aa"]("R,d,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += scalars_["1"] * r2["abab"]("R,b,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r1_1["bb"]("R,a,n") * reused_["217_bbaa_oovo"]("n,i,b,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= 0.50 * reused_["23_aa_vv"]("b,d") * r2_1["abab"]("R,d,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= reused_["123_aa_vv"]("b,d") * r2["abab"]("R,d,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["216_baab_vovv"]("e,j,b,a") * r1["bb"]("R,e,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["153_aa_vo"]("b,j") * r1["bb"]("R,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += eri["aaaa_vovo"]("b,k,c,j") * r2_1["abab"]("R,c,a,k,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += scalars_["2"] * r2["abab"]("R,b,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= dp["bb_oo"]("m,i") * r2["abab"]("R,b,a,j,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += eri["abab_vovo"]("b,m,c,i") * r2_1["abab"]("R,c,a,j,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["219_abba_vooo"]("b,i,n,j") * r1_1["bb"]("R,a,n");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["196_aabb_vvvo"]("b,d,a,i") * r1["aa"]("R,d,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["231_bbaa_vvvo"]("a,f,b,j") * r1_1["bb"]("R,f,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r0_1("R") * reused_["99_abab_vvoo"]("b,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += dp["aa_vv"]("b,c") * r2["abab"]("R,c,a,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += f["bb_oo"]("m,i") * r2_1["abab"]("R,b,a,j,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += t0_1 * r0_1("R") * reused_["74_baab_vvoo"]("a,b,j,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= 2.00 * r2_1["abab"]("R,b,a,k,i") * reused_["120_aa_oo"]("k,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= reused_["210_bbaa_vvvo"]("a,f,b,j") * r1_1["bb"]("R,f,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= 2.00 * reused_["208_aabb_vooo"]("b,j,i,m") * r1_1["bb"]("R,a,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r2_1["abab"]("R,b,f,j,n") * reused_["26_bbbb_voov"]("a,i,n,f");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r0_1("R") * reused_["81_abba_vovo"]("b,i,a,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= reused_["192_abba_vovv"]("b,i,a,d") * r1["aa"]("R,d,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += t1_1["aa"]("b,j") * tmps_["292_bb_Lvo"]("R,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += 0.50 * t2_1["abab"]("b,e,j,i") * tmps_["64_bb_Lvv"]("R,a,e");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += t2["abab"]("b,a,j,m") * tmps_["57_bb_Loo"]("R,i,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= t2_1["abab"]("b,a,j,m") * tmps_["53_bb_Loo"]("R,i,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= tmps_["288_aabb_Lvooo"]("R,b,j,n,i") * t1_1["bb"]("a,n");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += t2_1["abab"]("b,a,j,m") * tmps_["56_bb_Loo"]("R,i,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += t2["abab"]("b,a,j,m") * tmps_["55_bb_Loo"]("R,m,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= reused_["157_aa_vo"]("b,j") * r1["bb"]("R,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += tmps_["122_aa_Lvo"]("R,b,j") * reused_["32_bb_vo"]("a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r1_1["aa"]("R,b,j") * reused_["180_bb_vo"]("a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += t2_1["abab"]("b,a,j,i") * tmps_["43_L"]("R");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += tmps_["122_aa_Lvo"]("R,b,j") * dp["bb_vo"]("a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= t2_1["abab"]("b,e,j,m") * tmps_["277_bbbb_Lvoov"]("R,a,i,m,e");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += t2["abab"]("b,a,k,i") * tmps_["128_aa_Loo"]("R,j,k");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= tmps_["122_aa_Lvo"]("R,b,j") * reused_["33_bb_vo"]("a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= t2["abab"]("b,a,l,m") * tmps_["282_abab_Loooo"]("R,j,i,l,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= t1_1["aa"]("b,j") * tmps_["293_bb_Lov"]("R,i,a");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= t2_1["abab"]("b,a,k,i") * tmps_["121_aa_Loo"]("R,j,k");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= reused_["54_aa_vo"]("b,j") * tmps_["108_bb_Lvo"]("R,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += t1_1["aa"]("b,k") * tmps_["290_baba_Lvooo"]("R,a,j,i,k");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += t1_1["bb"]("a,m") * tmps_["291_aabb_Lvooo"]("R,b,j,i,m");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= t1_1["aa"]("b,j") * tmps_["295_bb_Lvo"]("R,a,i");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= t1_1["bb"]("a,i") * tmps_["166_aa_Lvo"]("R,b,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") -= r1_1["bb"]("R,a,i") * reused_["156_aa_vo"]("b,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += r1["bb"]("R,a,i") * reused_["155_aa_vo"]("b,j");
    tmps_["296_bbaa_Lvovo"]("R,a,i,b,j") += reused_["234_bb_vo"]("a,i") * tmps_["120_aa_Lvo"]("R,b,j");
    tmps_["291_aabb_Lvooo"].~TArrayD();
    tmps_["290_baba_Lvooo"].~TArrayD();
    tmps_["289_aabb_Lovvo"].~TArrayD();
    tmps_["288_aabb_Lvooo"].~TArrayD();
    tmps_["287_baab_Lvooo"].~TArrayD();
    tmps_["286_baab_Lvooo"].~TArrayD();
    tmps_["285_abab_Lvvoo"].~TArrayD();
    tmps_["284_bbbb_Lvoov"].~TArrayD();
    tmps_["283_baab_Lvoov"].~TArrayD();
    tmps_["282_abab_Loooo"].~TArrayD();
    tmps_["279_abab_Loooo"].~TArrayD();
    tmps_["278_baab_Lvoov"].~TArrayD();
    tmps_["277_bbbb_Lvoov"].~TArrayD();
    tmps_["166_aa_Lvo"].~TArrayD();
    tmps_["165_aa_Lov"].~TArrayD();
    tmps_["162_aa_Lvo"].~TArrayD();
    tmps_["154_bb_Lvv"].~TArrayD();
    tmps_["153_bb_Lvv"].~TArrayD();
    tmps_["147_aa_Lvv"].~TArrayD();
    tmps_["146_aa_Lvv"].~TArrayD();
    tmps_["135_aa_Lvv"].~TArrayD();
    tmps_["134_aa_Lvv"].~TArrayD();
    tmps_["129_aa_Loo"].~TArrayD();
    tmps_["128_aa_Loo"].~TArrayD();
    tmps_["122_aa_Lvo"].~TArrayD();
    tmps_["121_aa_Loo"].~TArrayD();
    tmps_["120_aa_Lvo"].~TArrayD();
    tmps_["65_bb_Lvv"].~TArrayD();
    tmps_["64_bb_Lvv"].~TArrayD();
    tmps_["57_bb_Loo"].~TArrayD();
    tmps_["56_bb_Loo"].~TArrayD();
    tmps_["55_bb_Loo"].~TArrayD();
    tmps_["54_bb_Loo"].~TArrayD();
    tmps_["53_bb_Loo"].~TArrayD();
    tmps_["49_bb_Loo"].~TArrayD();
    tmps_["43_L"].~TArrayD();
    tmps_["37_aa_Loo"].~TArrayD();
    tmps_["36_aa_Loo"].~TArrayD();
    tmps_["34_aa_Loo"].~TArrayD();

    // sigmar2_1_abab += +1.00 <l,k||c,d>_aaaa r1_bb(b,j) t2_aaaa(c,a,i,k) t1_1_aa(d,l)
    //                += +1.00 d-_aa(k,c) r1_1_bb(b,j) t2_aaaa(c,a,i,k) t0_1
    //                += +1.00 <l,k||c,d>_bbbb r1_bb(d,l) t2_abab(a,c,i,j) t1_1_bb(b,k)
    //                += -1.00 <a,k||i,c>_abab r1_bb(c,j) t1_1_bb(b,k)
    //                += -1.00 <a,k||c,j>_abab r1_aa(c,i) t1_1_bb(b,k)
    //                += +1.00 <l,k||c,j>_abab r2_aaaa(c,a,i,l) t1_1_bb(b,k)
    //                += +2.00 d-_bb(k,c) r2_1_abab(a,c,i,j) t1_1_bb(b,k)
    //                += -1.00 f_bb(k,c) r2_abab(a,c,i,j) t1_1_bb(b,k)
    //                += -1.00 <l,k||c,j>_bbbb r2_abab(a,c,i,l) t1_1_bb(b,k)
    //                += +1.00 <l,k||i,c>_abab r2_abab(a,c,l,j) t1_1_bb(b,k)
    //                += -0.50 <a,k||c,d>_abab r2_abab(c,d,i,j) t1_1_bb(b,k)
    //                += -0.50 <a,k||d,c>_abab r2_abab(d,c,i,j) t1_1_bb(b,k)
    //                += +1.00 d-_bb(k,c) r2_abab(a,c,i,j) t0_1 t1_1_bb(b,k)
    //                += -1.00 <l,k||d,c>_abab r1_aa(d,l) t2_abab(a,c,i,j) t1_1_bb(b,k)
    //                += -1.00 <k,l||c,d>_abab r1_bb(d,l) t2_abab(c,b,i,j) t1_1_aa(a,k)
    //                += +1.00 <l,k||c,d>_aaaa r1_aa(d,l) t2_abab(c,b,i,j) t1_1_aa(a,k)
    //                += +1.00 d-_bb(k,c) r1_bb(b,k) t1_1_aa(a,i) t1_1_bb(c,j)
    //                += +1.00 <a,k||c,d>_abab r1_1_bb(d,k) t2_abab(c,b,i,j)
    //                += -1.00 <k,a||c,d>_aaaa r1_1_aa(d,k) t2_abab(c,b,i,j)
    //                += +1.00 <k,l||d,c>_abab r2_1_abab(d,b,i,l) t2_abab(a,c,k,j)
    //                += +1.00 <k,b||d,c>_abab r1_aa(d,k) t2_1_abab(a,c,i,j)
    //                += -1.00 <k,b||c,d>_bbbb r1_bb(d,k) t2_1_abab(a,c,i,j)
    //                += +2.00 d-_aa(k,c) r1_1_aa(c,i) t2_1_abab(a,b,k,j)
    //                += -1.00 <l,k||c,i>_aaaa r1_1_aa(c,l) t2_abab(a,b,k,j)
    //                += -1.00 <k,l||i,c>_abab r1_1_bb(c,l) t2_abab(a,b,k,j)
    //                += -1.00 d-_bb(k,c) r2_1_abab(a,c,i,k) t1_1_bb(b,j)
    //                += +1.00 d-_aa(k,c) r2_1_aaaa(c,a,i,k) t1_1_bb(b,j)
    //                += -0.50 <l,k||d,c>_abab r2_1_abab(d,b,l,k) t2_abab(a,c,i,j)
    //                += -0.50 <k,l||d,c>_abab r2_1_abab(d,b,k,l) t2_abab(a,c,i,j)
    //                += +0.50 <l,k||c,d>_bbbb r2_1_bbbb(d,b,l,k) t2_abab(a,c,i,j)
    //                += -1.00 d-_aa(k,k) r0_1 t2_1_abab(a,b,i,j)
    //                += -0.50 <k,b||c,d>_abab r2_abab(c,d,i,j) t1_1_aa(a,k)
    //                += -0.50 <k,b||d,c>_abab r2_abab(d,c,i,j) t1_1_aa(a,k)
    //                += +1.00 d-_aa(k,c) r2_abab(c,b,i,j) t0_1 t1_1_aa(a,k)
    //                += +1.00 <k,l||c,j>_abab r2_abab(c,b,i,l) t1_1_aa(a,k)
    //                += -1.00 f_aa(k,c) r2_abab(c,b,i,j) t1_1_aa(a,k)
    //                += -1.00 <k,b||c,j>_abab r1_aa(c,i) t1_1_aa(a,k)
    //                += +1.00 <k,l||i,c>_abab r2_bbbb(c,b,j,l) t1_1_aa(a,k)
    //                += +0.50 <k,l||c,d>_abab r1_bb(b,l) t2_abab(c,d,i,j) t1_1_aa(a,k)
    //                += +0.50 <k,l||d,c>_abab r1_bb(b,l) t2_abab(d,c,i,j) t1_1_aa(a,k)
    //                += -1.00 <k,b||i,c>_abab r1_bb(c,j) t1_1_aa(a,k)
    //                += -1.00 <l,k||c,i>_aaaa r2_abab(c,b,l,j) t1_1_aa(a,k)
    //                += +1.00 <k,l||i,j>_abab r1_bb(b,l) t1_1_aa(a,k)
    //                += +2.00 d-_aa(k,c) r2_1_abab(c,b,i,j) t1_1_aa(a,k)
    //                += -1.00 d-_aa(a,i) r0_1 t1_1_bb(b,j)
    //                += -1.00 d-_bb(k,c) r0_1 t2_abab(a,c,i,k) t1_1_bb(b,j)
    //                += +1.00 <l,k||d,c>_abab r2_1_abab(d,b,l,j) t2_abab(a,c,i,k)
    //                += +0.25 <l,k||c,d>_abab r2_abab(c,d,i,j) t2_1_abab(a,b,l,k)
    //                += +0.25 <k,l||c,d>_abab r2_abab(c,d,i,j) t2_1_abab(a,b,k,l)
    //                += +0.25 <l,k||d,c>_abab r2_abab(d,c,i,j) t2_1_abab(a,b,l,k)
    //                += +0.25 <k,l||d,c>_abab r2_abab(d,c,i,j) t2_1_abab(a,b,k,l)
    //                += -1.00 d-_aa(k,c) r2_1_abab(c,b,k,j) t1_1_aa(a,i)
    //                += +1.00 d-_bb(k,c) r2_1_bbbb(c,b,j,k) t1_1_aa(a,i)
    //                += -1.00 <l,k||c,j>_abab r1_aa(c,l) t2_1_abab(a,b,i,k)
    //                += -1.00 <l,k||c,j>_bbbb r1_bb(c,l) t2_1_abab(a,b,i,k)
    //                += -0.50 <l,k||c,d>_abab r2_1_abab(a,d,l,k) t2_abab(c,b,i,j)
    //                += -0.50 <k,l||c,d>_abab r2_1_abab(a,d,k,l) t2_abab(c,b,i,j)
    //                += +0.50 <l,k||c,d>_aaaa r2_1_aaaa(d,a,l,k) t2_abab(c,b,i,j)
    //                += +1.00 d-_aa(k,i) r0_1 t2_abab(a,b,k,j) t0_1
    //                += +1.00 d-_aa(k,c) r1_aa(a,k) t0_1 t2_1_abab(c,b,i,j)
    //                += +1.00 d-_aa(k,c) r2_abab(a,b,k,j) t0_1 t1_1_aa(c,i)
    //                += +1.00 d-_aa(k,c) r1_1_aa(a,k) t2_abab(c,b,i,j) t0_1
    //                += +1.00 d-_bb(k,c) r2_abab(a,b,i,k) t0_1 t1_1_bb(c,j)
    //                += +1.00 d-_bb(k,c) r1_1_bb(b,k) t2_abab(a,c,i,j) t0_1
    //                += +1.00 d-_bb(k,c) r1_bb(b,k) t0_1 t2_1_abab(a,c,i,j)
    //                += +1.00 <l,k||c,d>_bbbb r1_aa(a,i) t2_bbbb(c,b,j,k) t1_1_bb(d,l)
    //                += -1.00 <l,k||c,d>_aaaa r1_aa(a,i) t2_abab(c,b,k,j) t1_1_aa(d,l)
    //                += +1.00 <l,k||c,d>_abab r1_bb(d,j) t2_abab(c,b,i,k) t1_1_aa(a,l)
    //                += +1.00 <l,k||c,d>_aaaa r1_aa(d,i) t2_abab(c,b,k,j) t1_1_aa(a,l)
    //                += +1.00 <l,k||d,c>_abab r1_aa(d,i) t2_bbbb(c,b,j,k) t1_1_aa(a,l)
    //                += -1.00 <k,b||c,d>_bbbb r1_1_bb(d,k) t2_abab(a,c,i,j)
    //                += +1.00 <k,b||d,c>_abab r1_1_aa(d,k) t2_abab(a,c,i,j)
    //                += +1.00 <k,l||d,c>_abab r2_abab(d,b,i,l) t2_1_abab(a,c,k,j)
    //                += -1.00 <k,l||i,c>_abab r1_bb(c,l) t2_1_abab(a,b,k,j)
    //                += -1.00 <l,k||c,i>_aaaa r1_aa(c,l) t2_1_abab(a,b,k,j)
    //                += +0.50 <l,k||c,d>_aaaa r2_aaaa(c,d,i,l) t2_1_abab(a,b,k,j)
    //                += -1.00 f_aa(k,c) r1_aa(c,i) t2_1_abab(a,b,k,j)
    //                += -0.50 <k,l||c,d>_abab r2_abab(c,d,i,l) t2_1_abab(a,b,k,j)
    //                += -0.50 <k,l||d,c>_abab r2_abab(d,c,i,l) t2_1_abab(a,b,k,j)
    //                += +2.00 d-_bb(k,c) r1_1_bb(c,j) t2_1_abab(a,b,i,k)
    //                += +0.50 <l,k||c,d>_aaaa r2_aaaa(d,a,l,k) t2_1_abab(c,b,i,j)
    //                += -0.50 <l,k||c,d>_abab r2_abab(a,d,l,k) t2_1_abab(c,b,i,j)
    //                += -0.50 <k,l||c,d>_abab r2_abab(a,d,k,l) t2_1_abab(c,b,i,j)
    //                += -1.00 <k,a||c,d>_aaaa r1_aa(d,k) t2_1_abab(c,b,i,j)
    //                += +1.00 <a,k||c,d>_abab r1_bb(d,k) t2_1_abab(c,b,i,j)
    //                += +1.00 d-_aa(k,i) r1_1_aa(a,k) t1_1_bb(b,j)
    //                += -0.50 <l,k||d,c>_abab r0_1 t2_abab(a,c,l,k) t2_abab(d,b,i,j)
    //                += -0.50 <k,l||d,c>_abab r0_1 t2_abab(a,c,k,l) t2_abab(d,b,i,j)
    //                += +1.00 <l,k||c,d>_aaaa r2_1_abab(d,b,l,j) t2_aaaa(c,a,i,k)
    //                += -1.00 <l,k||i,c>_abab r2_abab(a,b,l,j) t1_1_bb(c,k)
    //                += -0.50 <l,k||c,d>_abab r2_abab(a,b,l,j) t2_1_abab(c,d,i,k)
    //                += -0.50 <l,k||d,c>_abab r2_abab(a,b,l,j) t2_1_abab(d,c,i,k)
    //                += +0.50 <a,b||c,d>_abab r2_1_abab(c,d,i,j)
    //                += +0.50 <a,b||d,c>_abab r2_1_abab(d,c,i,j)
    //                += -0.50 <k,a||c,d>_aaaa r1_1_bb(b,j) t2_aaaa(c,d,i,k)
    //                += -0.50 <l,k||i,c>_abab r1_1_bb(b,j) t2_abab(a,c,l,k)
    //                += -0.50 <k,l||i,c>_abab r1_1_bb(b,j) t2_abab(a,c,k,l)
    //                += -0.50 <l,k||c,i>_aaaa r1_1_bb(b,j) t2_aaaa(c,a,l,k)
    //                += -1.00 f_aa(k,c) r1_1_bb(b,j) t2_aaaa(c,a,i,k)
    //                += +0.50 <a,k||c,d>_abab r1_1_bb(b,j) t2_abab(c,d,i,k)
    //                += +0.50 <a,k||d,c>_abab r1_1_bb(b,j) t2_abab(d,c,i,k)
    //                += -1.00 d-_aa(a,i) r1_1_bb(b,j) t0_1
    //                += +1.00 f_bb(k,c) r1_1_bb(b,j) t2_abab(a,c,i,k)
    //                += -1.00 d+_bb(k,c) r1_bb(b,j) t2_abab(a,c,i,k)
    //                += -1.00 d+_bb(b,j) r1_aa(a,i)
    //                += -2.00 d-_bb(b,c) r1_1_aa(a,i) t1_1_bb(c,j)
    //                += +2.00 d-_bb(k,j) r1_1_aa(a,i) t1_1_bb(b,k)
    //                += -2.00 d-_aa(k,c) r1_1_aa(a,i) t2_1_abab(c,b,k,j)
    //                += +2.00 d-_bb(k,c) r1_1_aa(a,i) t2_1_bbbb(c,b,j,k)
    //                += -1.00 <a,k||i,c>_abab r1_bb(b,k) t1_1_bb(c,j)
    //                += -1.00 <k,b||d,c>_abab r2_aaaa(d,a,i,k) t1_1_bb(c,j)
    //                += -0.50 <l,k||c,d>_abab r2_1_abab(a,d,i,j) t2_abab(c,b,l,k)
    //                += -0.50 <k,l||c,d>_abab r2_1_abab(a,d,i,j) t2_abab(c,b,k,l)
    //                += -0.50 <k,b||c,d>_abab r1_1_aa(a,k) t2_abab(c,d,i,j)
    //                += -0.50 <k,b||d,c>_abab r1_1_aa(a,k) t2_abab(d,c,i,j)
    //                += -1.00 f_aa(k,c) r1_1_aa(a,k) t2_abab(c,b,i,j)
    //                += -1.00 f_bb(k,c) r1_1_aa(a,i) t2_bbbb(c,b,j,k)
    //                += -1.00 d-_bb(b,j) r1_1_aa(a,i) t0_1
    //                += -0.50 <l,k||c,j>_abab r1_1_aa(a,i) t2_abab(c,b,l,k)
    //                += -0.50 <k,l||c,j>_abab r1_1_aa(a,i) t2_abab(c,b,k,l)
    //                += -0.50 <l,k||c,j>_bbbb r1_1_aa(a,i) t2_bbbb(c,b,l,k)
    //                += +1.00 f_aa(k,c) r1_1_aa(a,i) t2_abab(c,b,k,j)
    //                += +0.50 <k,b||c,d>_abab r1_1_aa(a,i) t2_abab(c,d,k,j)
    //                += +0.50 <k,b||d,c>_abab r1_1_aa(a,i) t2_abab(d,c,k,j)
    //                += -0.50 <k,b||c,d>_bbbb r1_1_aa(a,i) t2_bbbb(c,d,j,k)
    //                += -0.50 <k,l||c,d>_abab r2_1_abab(a,b,i,l) t2_abab(c,d,k,j)
    //                += -0.50 <k,l||d,c>_abab r2_1_abab(a,b,i,l) t2_abab(d,c,k,j)
    //                += -1.00 <a,k||d,c>_abab r2_abab(d,b,i,k) t1_1_bb(c,j)
    //                += +1.00 f_aa(a,c) r2_1_abab(c,b,i,j)
    //                += -2.00 d-_aa(k,c) r2_1_abab(a,b,i,j) t1_1_aa(c,k)
    //                += +1.00 <l,k||c,i>_aaaa r1_1_aa(a,l) t2_abab(c,b,k,j)
    //                += +1.00 <l,k||c,d>_abab r1_aa(a,l) t2_abab(c,b,i,k) t1_1_bb(d,j)
    //                += +1.00 <l,k||i,c>_abab r1_aa(a,l) t2_1_bbbb(c,b,j,k)
    //                += +1.00 <k,a||c,d>_aaaa r2_abab(d,b,i,j) t1_1_aa(c,k)
    //                += -0.50 <l,k||c,d>_aaaa r2_abab(d,b,i,j) t2_1_aaaa(c,a,l,k)
    //                += -1.00 <k,l||d,c>_abab r1_bb(b,l) t2_abab(a,c,i,j) t1_1_aa(d,k)
    //                += +1.00 <k,l||c,j>_abab r1_bb(b,l) t2_1_aaaa(c,a,i,k)
    //                += +1.00 <k,l||d,c>_abab r1_bb(b,l) t2_abab(a,c,k,j) t1_1_aa(d,i)
    //                += -1.00 d-_bb(k,k) r0_1 t2_1_abab(a,b,i,j)
    //                += +0.50 <l,k||c,d>_bbbb r0_1 t2_abab(a,c,i,j) t2_bbbb(d,b,l,k)
    //                += -1.00 d-_bb(k,c) r1_1_bb(b,j) t2_abab(a,c,i,k) t0_1
    //                += -1.00 d-_aa(k,k) r2_1_abab(a,b,i,j) t0_1
    //                += +1.00 f_aa(k,k) l2_1_aaaa(i,j,a,b)
    //                += -0.50 <l,k||l,k>_abab l2_1_aaaa(i,j,a,b)
    //                += -0.50 <k,l||k,l>_abab l2_1_aaaa(i,j,a,b)
    //                += +0.25 <l,k||c,d>_aaaa t2_aaaa(c,d,l,k) l2_1_aaaa(i,j,a,b)
    //                += -0.50 <l,k||l,k>_bbbb l2_1_aaaa(i,j,a,b)
    //                += -0.50 <l,k||l,k>_aaaa l2_1_aaaa(i,j,a,b)
    //                += +0.25 <l,k||c,d>_bbbb t2_bbbb(c,d,l,k) l2_1_aaaa(i,j,a,b)
    //                += +1.00 f_bb(k,k) l2_1_aaaa(i,j,a,b)
    //                += -1.00 d-_bb(k,k) t0_1 l2_1_aaaa(i,j,a,b)
    //                += +0.25 <l,k||c,d>_abab t2_abab(c,d,l,k) l2_1_aaaa(i,j,a,b)
    //                += +0.25 <k,l||c,d>_abab t2_abab(c,d,k,l) l2_1_aaaa(i,j,a,b)
    //                += +0.25 <l,k||d,c>_abab t2_abab(d,c,l,k) l2_1_aaaa(i,j,a,b)
    //                += +0.25 <k,l||d,c>_abab t2_abab(d,c,k,l) l2_1_aaaa(i,j,a,b)
    //                += +1.00 <l,k||c,i>_aaaa r2_abab(a,b,l,j) t1_1_aa(c,k)
    //                += -0.50 <l,k||c,d>_aaaa r2_abab(a,b,l,j) t2_1_aaaa(c,d,i,k)
    //                += +1.00 <k,a||c,d>_aaaa r1_1_aa(d,i) t2_abab(c,b,k,j)
    //                += +2.00 d-_aa(k,i) r1_1_bb(b,j) t1_1_aa(a,k)
    //                += +2.00 d-_aa(k,c) r1_1_bb(b,j) t2_1_aaaa(c,a,i,k)
    //                += -2.00 d-_aa(a,c) r1_1_bb(b,j) t1_1_aa(c,i)
    //                += -2.00 d-_bb(k,c) r1_1_bb(b,j) t2_1_abab(a,c,i,k)
    //                += -1.00 <k,b||i,j>_abab r1_1_aa(a,k)
    //                += -0.50 <l,k||d,c>_abab r2_1_abab(d,b,i,j) t2_abab(a,c,l,k)
    //                += -0.50 <k,l||d,c>_abab r2_1_abab(d,b,i,j) t2_abab(a,c,k,l)
    //                += -1.00 f_aa(k,c) r2_abab(a,b,k,j) t1_1_aa(c,i)
    //                += +1.00 <k,l||c,d>_abab r2_bbbb(d,b,j,l) t2_1_aaaa(c,a,i,k)
    //                += +2.00 d-_aa(k,c) r1_1_aa(a,k) t2_1_abab(c,b,i,j)
    //                += -1.00 <a,k||d,c>_abab r1_1_aa(d,i) t2_bbbb(c,b,j,k)
    //                += +1.00 <k,b||c,d>_abab r2_abab(a,d,i,j) t1_1_aa(c,k)
    //                += -0.50 <l,k||c,d>_abab r2_abab(a,d,i,j) t2_1_abab(c,b,l,k)
    //                += -0.50 <k,l||c,d>_abab r2_abab(a,d,i,j) t2_1_abab(c,b,k,l)
    //                += +1.00 r2_1_abab(a,b,i,j) w0
    //                += +1.00 <a,b||i,j>_abab r0_1
    //                += +1.00 <l,k||c,d>_aaaa r2_aaaa(d,a,i,l) t2_1_abab(c,b,k,j)
    //                += +1.00 <l,k||c,d>_abab r2_1_abab(a,d,l,j) t2_abab(c,b,i,k)
    //                += +1.00 <l,k||c,d>_bbbb r2_1_bbbb(d,b,j,l) t2_abab(a,c,i,k)
    //                += +1.00 <l,k||c,d>_bbbb r2_1_abab(a,d,i,l) t2_bbbb(c,b,j,k)
    //                += +1.00 <l,k||c,i>_aaaa r1_aa(a,l) t2_1_abab(c,b,k,j)
    //                += -1.00 <l,k||c,d>_aaaa r1_aa(a,l) t2_abab(c,b,i,j) t1_1_aa(d,k)
    //                += +1.00 <l,k||c,d>_aaaa r1_aa(a,l) t2_abab(c,b,k,j) t1_1_aa(d,i)
    //                += -1.00 <k,b||i,c>_abab r2_1_abab(a,c,k,j)
    //                += +1.00 <k,l||c,d>_abab r2_1_bbbb(d,b,j,l) t2_aaaa(c,a,i,k)
    //                += -0.50 <l,k||d,c>_abab r0_1 t2_abab(a,c,i,j) t2_abab(d,b,l,k)
    //                += -0.50 <k,l||d,c>_abab r0_1 t2_abab(a,c,i,j) t2_abab(d,b,k,l)
    //                += +1.00 <k,b||c,j>_bbbb r0_1 t2_abab(a,c,i,k)
    //                += -1.00 <k,b||c,j>_abab r0_1 t2_aaaa(c,a,i,k)
    //                += +1.00 f_bb(b,c) r0_1 t2_abab(a,c,i,j)
    //                += -2.00 d-_bb(b,c) r0_1 t2_1_abab(a,c,i,j)
    //                += +1.00 <k,l||c,d>_abab r0_1 t2_aaaa(c,a,i,k) t2_bbbb(d,b,j,l)
    //                += +2.00 d-_aa(k,c) r0_1 t2_abab(c,b,i,j) t1_1_aa(a,k)
    //                += +1.00 <l,k||c,d>_bbbb r0_1 t2_abab(a,c,i,k) t2_bbbb(d,b,j,l)
    //                += +1.00 <l,k||c,d>_aaaa r0_1 t2_aaaa(c,a,i,k) t2_abab(d,b,l,j)
    //                += -1.00 <k,b||c,d>_abab r2_abab(a,d,k,j) t1_1_aa(c,i)
    //                += +1.00 d+_aa(k,c) r1_bb(b,j) t2_aaaa(c,a,i,k)
    //                += -1.00 d-_aa(k,k) r1_1_aa(a,i) t1_1_bb(b,j)
    //                += -1.00 d-_bb(k,k) r1_1_aa(a,i) t1_1_bb(b,j)
    //                += -1.00 d-_bb(b,c) r2_1_abab(a,c,i,j) t0_1
    //                += -1.00 <a,k||c,j>_abab r1_bb(b,k) t1_1_aa(c,i)
    //                += -1.00 f_bb(k,c) r1_bb(b,k) t2_1_abab(a,c,i,j)
    //                += -0.50 <a,k||c,d>_abab r1_bb(b,k) t2_1_abab(c,d,i,j)
    //                += -0.50 <a,k||d,c>_abab r1_bb(b,k) t2_1_abab(d,c,i,j)
    //                += +1.00 f_bb(b,j) r1_1_aa(a,i)
    //                += +1.00 f_bb(b,c) r2_1_abab(a,c,i,j)
    //                += -1.00 <k,a||c,d>_aaaa r2_abab(d,b,k,j) t1_1_aa(c,i)
    //                += -1.00 f_aa(k,i) r2_1_abab(a,b,k,j)
    //                += -1.00 <a,k||c,d>_abab r1_1_bb(d,j) t2_abab(c,b,i,k)
    //                += -1.00 <a,k||c,d>_abab r1_bb(d,j) t2_1_abab(c,b,i,k)
    //                += +1.00 <a,b||c,d>_abab r1_bb(d,j) t1_1_aa(c,i)
    //                += +0.50 <l,k||c,d>_abab r1_bb(d,j) t2_abab(a,b,l,k) t1_1_aa(c,i)
    //                += +0.50 <k,l||c,d>_abab r1_bb(d,j) t2_abab(a,b,k,l) t1_1_aa(c,i)
    //                += -1.00 d+_aa(k,c) r1_aa(a,i) t2_abab(c,b,k,j)
    //                += +1.00 <a,b||i,c>_abab r1_1_bb(c,j)
    //                += +1.00 <l,k||c,d>_abab r2_abab(a,d,l,j) t2_1_abab(c,b,i,k)
    //                += -1.00 <k,l||c,j>_abab r2_abab(a,b,i,l) t1_1_aa(c,k)
    //                += -0.50 <k,l||c,d>_abab r2_abab(a,b,i,l) t2_1_abab(c,d,k,j)
    //                += -0.50 <k,l||d,c>_abab r2_abab(a,b,i,l) t2_1_abab(d,c,k,j)
    //                += +0.50 <l,k||i,j>_abab r2_1_abab(a,b,l,k)
    //                += +0.50 <k,l||i,j>_abab r2_1_abab(a,b,k,l)
    //                += -0.50 <l,k||c,d>_bbbb r2_abab(a,d,i,j) t2_1_bbbb(c,b,l,k)
    //                += +1.00 <k,b||c,d>_bbbb r2_abab(a,d,i,j) t1_1_bb(c,k)
    //                += +1.00 <l,k||c,d>_bbbb r2_abab(a,d,i,l) t2_1_bbbb(c,b,j,k)
    //                += +1.00 <l,k||i,c>_abab r1_1_aa(a,l) t2_bbbb(c,b,j,k)
    //                += +1.00 d+_aa(k,i) r2_abab(a,b,k,j)
    //                += +1.00 d+_bb(k,c) r1_aa(a,i) t2_bbbb(c,b,j,k)
    //                += +2.00 d-_bb(k,c) r2_1_abab(a,b,i,k) t1_1_bb(c,j)
    //                += -0.50 <l,k||c,d>_aaaa r2_1_abab(a,b,l,j) t2_aaaa(c,d,i,k)
    //                += -1.00 d-_bb(k,k) r1_1_bb(b,j) t1_1_aa(a,i)
    //                += -1.00 d-_aa(k,k) r1_1_bb(b,j) t1_1_aa(a,i)
    //                += +0.25 <l,k||c,d>_abab r2_abab(a,b,l,k) t2_1_abab(c,d,i,j)
    //                += +0.25 <l,k||d,c>_abab r2_abab(a,b,l,k) t2_1_abab(d,c,i,j)
    //                += +0.25 <k,l||c,d>_abab r2_abab(a,b,k,l) t2_1_abab(c,d,i,j)
    //                += +0.25 <k,l||d,c>_abab r2_abab(a,b,k,l) t2_1_abab(d,c,i,j)
    //                += +0.50 <l,k||c,j>_abab r2_abab(a,b,l,k) t1_1_aa(c,i)
    //                += +0.50 <k,l||c,j>_abab r2_abab(a,b,k,l) t1_1_aa(c,i)
    //                += -1.00 d-_aa(a,c) r0_1 t2_abab(c,b,i,j) t0_1
    //                += +1.00 d-_bb(k,j) r0_1 t2_abab(a,b,i,k) t0_1
    //                += -1.00 d-_aa(a,c) r2_1_abab(c,b,i,j) t0_1
    //                += +1.00 <l,k||c,d>_bbbb r2_bbbb(d,b,j,l) t2_1_abab(a,c,i,k)
    //                += -1.00 f_bb(k,c) r1_1_bb(b,k) t2_abab(a,c,i,j)
    //                += -0.50 <a,k||c,d>_abab r1_1_bb(b,k) t2_abab(c,d,i,j)
    //                += -0.50 <a,k||d,c>_abab r1_1_bb(b,k) t2_abab(d,c,i,j)
    //                += -1.00 <k,b||c,j>_abab r2_1_aaaa(c,a,i,k)
    //                += +1.00 <l,k||c,d>_aaaa r2_1_aaaa(d,a,i,l) t2_abab(c,b,k,j)
    //                += -1.00 <k,b||c,d>_bbbb r2_abab(a,d,i,k) t1_1_bb(c,j)
    //                += +0.25 <l,k||c,d>_abab r2_1_abab(a,b,l,k) t2_abab(c,d,i,j)
    //                += +0.25 <l,k||d,c>_abab r2_1_abab(a,b,l,k) t2_abab(d,c,i,j)
    //                += +0.25 <k,l||c,d>_abab r2_1_abab(a,b,k,l) t2_abab(c,d,i,j)
    //                += +0.25 <k,l||d,c>_abab r2_1_abab(a,b,k,l) t2_abab(d,c,i,j)
    //                += +1.00 <l,k||c,d>_bbbb r1_bb(b,l) t2_abab(a,c,i,k) t1_1_bb(d,j)
    //                += -1.00 <l,k||c,d>_bbbb r1_bb(b,l) t2_abab(a,c,i,j) t1_1_bb(d,k)
    //                += +1.00 <l,k||c,j>_bbbb r1_bb(b,l) t2_1_abab(a,c,i,k)
    //                += +0.50 <l,k||c,j>_abab r1_1_aa(c,i) t2_abab(a,b,l,k)
    //                += +0.50 <k,l||c,j>_abab r1_1_aa(c,i) t2_abab(a,b,k,l)
    //                += +1.00 <l,k||c,d>_aaaa r2_abab(d,b,l,j) t2_1_aaaa(c,a,i,k)
    //                += -2.00 d-_bb(k,c) r2_1_abab(a,b,i,j) t1_1_bb(c,k)
    //                += -1.00 d+_bb(b,c) r2_abab(a,c,i,j)
    //                += -1.00 <a,k||c,d>_abab r2_bbbb(d,b,j,k) t1_1_aa(c,i)
    //                += -1.00 d+_aa(a,i) r1_bb(b,j)
    //                += -1.00 <l,k||c,d>_abab r1_aa(a,l) t2_abab(c,b,i,j) t1_1_bb(d,k)
    //                += +1.00 <l,k||c,j>_abab r1_aa(a,l) t2_1_abab(c,b,i,k)
    //                += +1.00 <l,k||i,j>_abab r1_aa(a,l) t1_1_bb(b,k)
    //                += +0.50 <l,k||c,d>_abab r1_aa(a,l) t2_abab(c,d,i,j) t1_1_bb(b,k)
    //                += +0.50 <l,k||d,c>_abab r1_aa(a,l) t2_abab(d,c,i,j) t1_1_bb(b,k)
    //                += +1.00 <l,k||d,c>_abab r1_aa(a,l) t2_bbbb(c,b,j,k) t1_1_aa(d,i)
    //                += +1.00 <a,b||c,j>_abab r1_1_aa(c,i)
    //                += +1.00 <k,l||c,d>_abab r2_abab(a,d,i,l) t2_1_abab(c,b,k,j)
    //                += +1.00 <l,k||d,c>_abab r2_1_aaaa(d,a,i,l) t2_bbbb(c,b,j,k)
    //                += +1.00 <k,b||c,j>_bbbb r2_1_abab(a,c,i,k)
    //                += +1.00 <l,k||c,j>_abab r1_1_aa(a,l) t2_abab(c,b,i,k)
    //                += -0.50 <l,k||c,d>_abab r2_1_abab(a,b,l,j) t2_abab(c,d,i,k)
    //                += -0.50 <l,k||d,c>_abab r2_1_abab(a,b,l,j) t2_abab(d,c,i,k)
    //                += +1.00 d-_aa(k,i) r2_1_abab(a,b,k,j) t0_1
    //                += +1.00 d-_bb(k,j) r2_1_abab(a,b,i,k) t0_1
    //                += -0.50 <l,k||c,d>_abab r0_1 t2_abab(a,b,l,j) t2_abab(c,d,i,k)
    //                += -0.50 <l,k||d,c>_abab r0_1 t2_abab(a,b,l,j) t2_abab(d,c,i,k)
    //                += +1.00 <k,a||c,i>_aaaa r0_1 t2_abab(c,b,k,j)
    //                += +2.00 d-_aa(k,i) r0_1 t2_1_abab(a,b,k,j)
    //                += +2.00 d-_bb(k,c) r0_1 t2_abab(a,b,i,k) t1_1_bb(c,j)
    //                += -1.00 <a,k||i,c>_abab r0_1 t2_bbbb(c,b,j,k)
    //                += -1.00 f_aa(k,i) r0_1 t2_abab(a,b,k,j)
    //                += +1.00 <l,k||d,c>_abab r0_1 t2_abab(a,c,i,k) t2_abab(d,b,l,j)
    //                += -0.50 <l,k||c,d>_aaaa r0_1 t2_abab(a,b,l,j) t2_aaaa(c,d,i,k)
    //                += -0.50 <l,k||c,d>_bbbb r2_abab(a,b,i,l) t2_1_bbbb(c,d,j,k)
    //                += +1.00 <l,k||c,j>_bbbb r2_abab(a,b,i,l) t1_1_bb(c,k)
    //                += -1.00 <a,k||i,c>_abab r2_1_bbbb(c,b,j,k)
    //                += +0.50 <l,k||c,j>_abab r1_aa(c,i) t2_1_abab(a,b,l,k)
    //                += +0.50 <k,l||c,j>_abab r1_aa(c,i) t2_1_abab(a,b,k,l)
    //                += +1.00 <k,l||i,c>_abab r1_bb(b,l) t2_1_abab(a,c,k,j)
    //                += +1.00 <k,l||c,d>_abab r1_bb(b,l) t2_aaaa(c,a,i,k) t1_1_bb(d,j)
    //                += -1.00 <k,b||i,c>_abab r1_aa(a,k) t1_1_bb(c,j)
    //                += +0.50 <l,k||i,c>_abab r2_abab(a,b,l,k) t1_1_bb(c,j)
    //                += +0.50 <k,l||i,c>_abab r2_abab(a,b,k,l) t1_1_bb(c,j)
    //                += -0.50 <l,k||c,d>_bbbb r1_aa(a,i) t2_bbbb(c,d,j,k) t1_1_bb(b,l)
    //                += -1.00 <l,k||d,c>_abab r1_aa(a,i) t2_bbbb(c,b,j,k) t1_1_aa(d,l)
    //                += -0.50 <l,k||c,d>_abab r1_aa(a,i) t2_abab(c,b,l,k) t1_1_bb(d,j)
    //                += -0.50 <k,l||c,d>_abab r1_aa(a,i) t2_abab(c,b,k,l) t1_1_bb(d,j)
    //                += +1.00 f_aa(k,c) r1_aa(a,i) t2_1_abab(c,b,k,j)
    //                += -1.00 f_bb(k,j) r1_aa(a,i) t1_1_bb(b,k)
    //                += -0.50 <k,b||c,d>_bbbb r1_aa(a,i) t2_1_bbbb(c,d,j,k)
    //                += +1.00 <k,b||c,j>_abab r1_aa(a,i) t1_1_aa(c,k)
    //                += -1.00 f_bb(k,c) r1_aa(a,i) t2_1_bbbb(c,b,j,k)
    //                += +1.00 f_bb(b,c) r1_aa(a,i) t1_1_bb(c,j)
    //                += +0.50 <k,b||c,d>_abab r1_aa(a,i) t2_1_abab(c,d,k,j)
    //                += +0.50 <k,b||d,c>_abab r1_aa(a,i) t2_1_abab(d,c,k,j)
    //                += +1.00 <k,b||c,j>_bbbb r1_aa(a,i) t1_1_bb(c,k)
    //                += -0.50 <l,k||c,j>_bbbb r1_aa(a,i) t2_1_bbbb(c,b,l,k)
    //                += -1.00 d-_bb(k,c) r1_aa(a,i) t1_1_bb(b,j) t1_1_bb(c,k)
    //                += -1.00 d-_aa(k,c) r1_aa(a,i) t1_1_bb(b,j) t1_1_aa(c,k)
    //                += +1.00 r1_aa(a,i) t1_1_bb(b,j) w0
    //                += -0.50 <l,k||c,j>_abab r1_aa(a,i) t2_1_abab(c,b,l,k)
    //                += -0.50 <k,l||c,j>_abab r1_aa(a,i) t2_1_abab(c,b,k,l)
    //                += -0.50 <k,l||c,d>_abab r1_aa(a,i) t2_abab(c,d,k,j) t1_1_bb(b,l)
    //                += -0.50 <k,l||d,c>_abab r1_aa(a,i) t2_abab(d,c,k,j) t1_1_bb(b,l)
    //                += +2.00 d-_bb(k,c) r1_aa(a,i) t1_1_bb(b,k) t1_1_bb(c,j)
    //                += -1.00 d-_bb(b,c) r1_aa(a,i) t0_1 t1_1_bb(c,j)
    //                += +1.00 d-_bb(k,j) r1_aa(a,i) t0_1 t1_1_bb(b,k)
    //                += -1.00 d-_aa(k,c) r1_aa(a,i) t0_1 t2_1_abab(c,b,k,j)
    //                += +1.00 d-_bb(k,c) r1_aa(a,i) t0_1 t2_1_bbbb(c,b,j,k)
    //                += -0.50 <l,k||c,d>_bbbb r1_aa(a,i) t2_bbbb(c,b,l,k) t1_1_bb(d,j)
    //                += +1.00 <k,l||c,d>_abab r1_aa(a,i) t2_abab(c,b,k,j) t1_1_bb(d,l)
    //                += +1.00 <k,a||c,d>_aaaa r1_aa(d,i) t2_1_abab(c,b,k,j)
    //                += -1.00 <a,k||i,j>_abab r1_1_bb(b,k)
    //                += +0.50 <l,k||i,c>_abab r1_1_bb(c,j) t2_abab(a,b,l,k)
    //                += +0.50 <k,l||i,c>_abab r1_1_bb(c,j) t2_abab(a,b,k,l)
    //                += +1.00 <k,b||c,d>_bbbb r1_bb(d,j) t2_1_abab(a,c,i,k)
    //                += -1.00 f_aa(k,c) r1_aa(a,k) t2_1_abab(c,b,i,j)
    //                += -1.00 <k,b||c,j>_abab r1_aa(a,k) t1_1_aa(c,i)
    //                += -0.50 <k,b||c,d>_abab r1_aa(a,k) t2_1_abab(c,d,i,j)
    //                += -0.50 <k,b||d,c>_abab r1_aa(a,k) t2_1_abab(d,c,i,j)
    //                += -0.50 <l,k||c,d>_bbbb r2_1_abab(a,b,i,l) t2_bbbb(c,d,j,k)
    //                += +1.00 <l,k||d,c>_abab r2_aaaa(d,a,i,l) t2_1_bbbb(c,b,j,k)
    //                += -1.00 f_bb(k,c) r2_abab(a,b,i,k) t1_1_bb(c,j)
    //                += +1.00 <l,k||c,j>_bbbb r1_1_bb(b,l) t2_abab(a,c,i,k)
    //                += +1.00 r2_abab(a,b,i,j) t0_1 w0
    //                += +0.25 <l,k||c,d>_aaaa r2_abab(a,b,i,j) t2_1_aaaa(c,d,l,k)
    //                += +0.25 <l,k||c,d>_abab r2_abab(a,b,i,j) t2_1_abab(c,d,l,k)
    //                += +0.25 <k,l||c,d>_abab r2_abab(a,b,i,j) t2_1_abab(c,d,k,l)
    //                += +0.25 <l,k||d,c>_abab r2_abab(a,b,i,j) t2_1_abab(d,c,l,k)
    //                += +0.25 <k,l||d,c>_abab r2_abab(a,b,i,j) t2_1_abab(d,c,k,l)
    //                += -1.00 d-_aa(k,c) r2_abab(a,b,i,j) t0_1 t1_1_aa(c,k)
    //                += +1.00 f_bb(k,c) r2_abab(a,b,i,j) t1_1_bb(c,k)
    //                += -1.00 d-_bb(k,c) r2_abab(a,b,i,j) t0_1 t1_1_bb(c,k)
    //                += +0.25 <l,k||c,d>_bbbb r2_abab(a,b,i,j) t2_1_bbbb(c,d,l,k)
    //                += +1.00 f_aa(k,c) r2_abab(a,b,i,j) t1_1_aa(c,k)
    //                += +1.00 f_aa(a,i) r1_1_bb(b,j)
    //                += -1.00 <k,b||c,d>_abab r1_bb(d,j) t2_1_aaaa(c,a,i,k)
    //                += -0.50 <l,k||c,d>_bbbb r2_1_abab(a,d,i,j) t2_bbbb(c,b,l,k)
    //                += -1.00 <k,b||d,c>_abab r1_1_aa(d,i) t2_abab(a,c,k,j)
    //                += -1.00 d+_aa(k,k) r2_abab(a,b,i,j)
    //                += +1.00 <k,l||c,j>_abab r1_1_bb(b,l) t2_aaaa(c,a,i,k)
    //                += -0.50 <l,k||c,d>_aaaa r2_1_abab(d,b,i,j) t2_aaaa(c,a,l,k)
    //                += +1.00 <a,k||d,c>_abab r2_abab(d,b,i,j) t1_1_bb(c,k)
    //                += -0.50 <l,k||d,c>_abab r2_abab(d,b,i,j) t2_1_abab(a,c,l,k)
    //                += -0.50 <k,l||d,c>_abab r2_abab(d,b,i,j) t2_1_abab(a,c,k,l)
    //                += +0.50 <l,k||i,c>_abab r1_bb(c,j) t2_1_abab(a,b,l,k)
    //                += +0.50 <k,l||i,c>_abab r1_bb(c,j) t2_1_abab(a,b,k,l)
    //                += -0.50 <l,k||d,c>_abab r1_bb(b,j) t2_abab(a,c,l,k) t1_1_aa(d,i)
    //                += -0.50 <k,l||d,c>_abab r1_bb(b,j) t2_abab(a,c,k,l) t1_1_aa(d,i)
    //                += +1.00 d-_aa(k,i) r1_bb(b,j) t0_1 t1_1_aa(a,k)
    //                += +1.00 d-_aa(k,c) r1_bb(b,j) t0_1 t2_1_aaaa(c,a,i,k)
    //                += -1.00 d-_aa(a,c) r1_bb(b,j) t0_1 t1_1_aa(c,i)
    //                += -1.00 d-_bb(k,c) r1_bb(b,j) t0_1 t2_1_abab(a,c,i,k)
    //                += -1.00 <k,l||c,d>_abab r1_bb(b,j) t2_aaaa(c,a,i,k) t1_1_bb(d,l)
    //                += +2.00 d-_aa(k,c) r1_bb(b,j) t1_1_aa(a,k) t1_1_aa(c,i)
    //                += -0.50 <l,k||c,d>_abab r1_bb(b,j) t2_abab(c,d,i,k) t1_1_aa(a,l)
    //                += -0.50 <l,k||d,c>_abab r1_bb(b,j) t2_abab(d,c,i,k) t1_1_aa(a,l)
    //                += -0.50 <l,k||c,d>_aaaa r1_bb(b,j) t2_aaaa(c,a,l,k) t1_1_aa(d,i)
    //                += +1.00 r1_bb(b,j) t1_1_aa(a,i) w0
    //                += -1.00 f_aa(k,c) r1_bb(b,j) t2_1_aaaa(c,a,i,k)
    //                += -1.00 d-_bb(k,c) r1_bb(b,j) t1_1_aa(a,i) t1_1_bb(c,k)
    //                += -0.50 <l,k||i,c>_abab r1_bb(b,j) t2_1_abab(a,c,l,k)
    //                += -0.50 <k,l||i,c>_abab r1_bb(b,j) t2_1_abab(a,c,k,l)
    //                += -1.00 d-_aa(k,c) r1_bb(b,j) t1_1_aa(a,i) t1_1_aa(c,k)
    //                += +1.00 <k,a||c,i>_aaaa r1_bb(b,j) t1_1_aa(c,k)
    //                += +1.00 f_bb(k,c) r1_bb(b,j) t2_1_abab(a,c,i,k)
    //                += -1.00 f_aa(k,i) r1_bb(b,j) t1_1_aa(a,k)
    //                += -0.50 <k,a||c,d>_aaaa r1_bb(b,j) t2_1_aaaa(c,d,i,k)
    //                += +1.00 <a,k||i,c>_abab r1_bb(b,j) t1_1_bb(c,k)
    //                += -0.50 <l,k||c,i>_aaaa r1_bb(b,j) t2_1_aaaa(c,a,l,k)
    //                += +0.50 <a,k||c,d>_abab r1_bb(b,j) t2_1_abab(c,d,i,k)
    //                += +0.50 <a,k||d,c>_abab r1_bb(b,j) t2_1_abab(d,c,i,k)
    //                += +1.00 f_aa(a,c) r1_bb(b,j) t1_1_aa(c,i)
    //                += -0.50 <l,k||c,d>_aaaa r1_bb(b,j) t2_aaaa(c,d,i,k) t1_1_aa(a,l)
    //                += +1.00 <k,a||c,i>_aaaa r2_1_abab(c,b,k,j)
    //                += -1.00 d+_bb(k,k) r2_abab(a,b,i,j)
    //                += +1.00 d+_bb(k,j) r2_abab(a,b,i,k)
    //                += -1.00 <a,k||c,j>_abab r2_1_abab(c,b,i,k)
    //                += +1.00 <k,l||i,c>_abab r1_1_bb(b,l) t2_abab(a,c,k,j)
    //                += -1.00 <a,k||d,c>_abab r1_aa(d,i) t2_1_bbbb(c,b,j,k)
    //                += +0.50 <l,k||d,c>_abab r1_aa(d,i) t2_abab(a,b,l,k) t1_1_bb(c,j)
    //                += +0.50 <k,l||d,c>_abab r1_aa(d,i) t2_abab(a,b,k,l) t1_1_bb(c,j)
    //                += +1.00 <a,b||d,c>_abab r1_aa(d,i) t1_1_bb(c,j)
    //                += +1.00 <k,b||c,d>_bbbb r1_1_bb(d,j) t2_abab(a,c,i,k)
    //                += -0.50 <k,l||c,d>_abab r0_1 t2_abab(a,b,i,l) t2_abab(c,d,k,j)
    //                += -0.50 <k,l||d,c>_abab r0_1 t2_abab(a,b,i,l) t2_abab(d,c,k,j)
    //                += +0.50 <l,k||i,j>_abab r0_1 t2_abab(a,b,l,k)
    //                += +0.50 <k,l||i,j>_abab r0_1 t2_abab(a,b,k,l)
    //                += +0.50 <a,b||c,d>_abab r0_1 t2_abab(c,d,i,j)
    //                += +0.50 <a,b||d,c>_abab r0_1 t2_abab(d,c,i,j)
    //                += +1.00 f_aa(a,c) r0_1 t2_abab(c,b,i,j)
    //                += +2.00 d-_bb(k,j) r0_1 t2_1_abab(a,b,i,k)
    //                += -1.00 f_bb(k,j) r0_1 t2_abab(a,b,i,k)
    //                += +2.00 d-_aa(k,c) r0_1 t2_abab(a,b,k,j) t1_1_aa(c,i)
    //                += -0.50 <l,k||c,d>_bbbb r0_1 t2_abab(a,b,i,l) t2_bbbb(c,d,j,k)
    //                += -2.00 d-_aa(a,c) r0_1 t2_1_abab(c,b,i,j)
    //                += -1.00 <a,k||c,j>_abab r0_1 t2_abab(c,b,i,k)
    //                += +1.00 <k,l||d,c>_abab r0_1 t2_abab(a,c,k,j) t2_abab(d,b,i,l)
    //                += +2.00 d-_bb(k,c) r0_1 t2_abab(a,c,i,j) t1_1_bb(b,k)
    //                += +0.25 <l,k||c,d>_abab r0_1 t2_abab(a,b,l,k) t2_abab(c,d,i,j)
    //                += +0.25 <l,k||d,c>_abab r0_1 t2_abab(a,b,l,k) t2_abab(d,c,i,j)
    //                += +0.25 <k,l||c,d>_abab r0_1 t2_abab(a,b,k,l) t2_abab(c,d,i,j)
    //                += +0.25 <k,l||d,c>_abab r0_1 t2_abab(a,b,k,l) t2_abab(d,c,i,j)
    //                += -0.50 <l,k||c,d>_aaaa r0_1 t2_aaaa(c,a,l,k) t2_abab(d,b,i,j)
    //                += -1.00 d+_aa(a,c) r2_abab(c,b,i,j)
    //                += -1.00 f_bb(k,j) r2_1_abab(a,b,i,k)
    //                += -1.00 d-_bb(b,c) r0_1 t2_abab(a,c,i,j) t0_1
    //                += +2.00 d-_aa(k,c) r2_1_abab(a,b,k,j) t1_1_aa(c,i)
    //                += -1.00 <k,b||c,d>_abab r1_1_bb(d,j) t2_aaaa(c,a,i,k)
    //                += +2.00 d-_bb(k,c) r1_1_bb(b,k) t2_1_abab(a,c,i,j)
    //                += +1.00 <k,l||c,d>_abab r2_1_abab(a,d,i,l) t2_abab(c,b,k,j)
    //                += -1.00 <k,b||i,c>_abab r0_1 t2_abab(a,c,k,j)
    //                += -1.00 <k,b||d,c>_abab r1_aa(d,i) t2_1_abab(a,c,k,j)
    //                += -1.00 d-_bb(b,c) r1_1_bb(c,j) t1_1_aa(a,i)
    //                += +0.50 <l,k||c,d>_bbbb r2_bbbb(d,b,l,k) t2_1_abab(a,c,i,j)
    //                += -0.50 <l,k||d,c>_abab r2_abab(d,b,l,k) t2_1_abab(a,c,i,j)
    //                += -0.50 <k,l||d,c>_abab r2_abab(d,b,k,l) t2_1_abab(a,c,i,j)
    //                += -1.00 f_bb(k,c) r1_1_bb(c,j) t2_abab(a,b,i,k)
    //                += -1.00 <l,k||c,d>_abab r1_bb(d,j) t2_abab(a,b,i,k) t1_1_aa(c,l)
    //                += -1.00 <l,k||c,d>_bbbb r1_bb(d,j) t2_abab(a,b,i,k) t1_1_bb(c,l)
    //                += +1.00 d-_bb(k,c) r1_1_bb(c,j) t2_abab(a,b,i,k) t0_1
    //                += -0.50 <l,k||c,d>_abab r2_1_abab(c,d,l,j) t2_abab(a,b,i,k)
    //                += -0.50 <l,k||d,c>_abab r2_1_abab(d,c,l,j) t2_abab(a,b,i,k)
    //                += +1.00 <l,k||c,d>_bbbb r1_bb(d,l) t2_abab(a,b,i,k) t1_1_bb(c,j)
    //                += +0.50 <l,k||c,d>_bbbb r2_1_bbbb(c,d,j,l) t2_abab(a,b,i,k)
    //                += -1.00 <l,k||d,c>_abab r1_aa(d,l) t2_abab(a,b,i,k) t1_1_bb(c,j)
    //                += +1.00 d-_bb(k,c) r1_bb(c,j) t0_1 t2_1_abab(a,b,i,k)
    //                += +1.00 <k,l||c,d>_abab r1_bb(d,j) t2_aaaa(c,a,i,k) t1_1_bb(b,l)
    //                += +1.00 <k,l||d,c>_abab r1_aa(d,i) t2_abab(a,c,k,j) t1_1_bb(b,l)
    //                += +1.00 <l,k||c,d>_bbbb r1_bb(d,j) t2_abab(a,c,i,k) t1_1_bb(b,l)
    //                += -1.00 f_bb(k,c) r1_bb(c,j) t2_1_abab(a,b,i,k)
    //                += +0.50 <l,k||c,d>_bbbb r2_bbbb(c,d,j,l) t2_1_abab(a,b,i,k)
    //                += -0.50 <l,k||c,d>_abab r2_abab(c,d,l,j) t2_1_abab(a,b,i,k)
    //                += -0.50 <l,k||d,c>_abab r2_abab(d,c,l,j) t2_1_abab(a,b,i,k)
    //                += -1.00 <l,k||c,j>_abab r1_1_aa(c,l) t2_abab(a,b,i,k)
    //                += -1.00 <l,k||c,j>_bbbb r1_1_bb(c,l) t2_abab(a,b,i,k)
    //                += +1.00 <l,k||d,c>_abab r1_bb(b,j) t2_abab(a,c,i,k) t1_1_aa(d,l)
    //                += -1.00 <l,k||c,d>_bbbb r1_bb(b,j) t2_abab(a,c,i,k) t1_1_bb(d,l)
    //                += -1.00 d-_aa(k,c) r0_1 t2_abab(c,b,k,j) t1_1_aa(a,i)
    //                += -1.00 d-_aa(k,c) r1_1_aa(a,i) t2_abab(c,b,k,j) t0_1
    //                += +1.00 d-_bb(k,c) r1_1_aa(a,i) t2_bbbb(c,b,j,k) t0_1
    //                += -1.00 d-_bb(k,c) r1_1_bb(c,k) t2_1_abab(a,b,i,j)
    //                += -1.00 d-_aa(k,c) r1_1_aa(c,k) t2_1_abab(a,b,i,j)
    //                += -1.00 d-_bb(b,j) r0_1 t1_1_aa(a,i)
    //                += +1.00 <l,k||d,c>_abab r2_abab(d,b,l,j) t2_1_abab(a,c,i,k)
    //                += -0.50 <k,l||c,d>_abab r2_1_abab(c,d,i,l) t2_abab(a,b,k,j)
    //                += -0.50 <k,l||d,c>_abab r2_1_abab(d,c,i,l) t2_abab(a,b,k,j)
    //                += -1.00 <l,k||c,d>_aaaa r1_aa(d,i) t2_abab(a,b,k,j) t1_1_aa(c,l)
    //                += +0.50 <l,k||c,d>_aaaa r2_1_aaaa(c,d,i,l) t2_abab(a,b,k,j)
    //                += +1.00 d-_aa(k,c) r1_1_aa(c,i) t2_abab(a,b,k,j) t0_1
    //                += -1.00 <k,l||c,d>_abab r1_bb(d,l) t2_abab(a,b,k,j) t1_1_aa(c,i)
    //                += -1.00 f_aa(k,c) r1_1_aa(c,i) t2_abab(a,b,k,j)
    //                += -1.00 <k,l||d,c>_abab r1_aa(d,i) t2_abab(a,b,k,j) t1_1_bb(c,l)
    //                += +1.00 <l,k||c,d>_aaaa r1_aa(d,l) t2_abab(a,b,k,j) t1_1_aa(c,i)
    //                += +1.00 d-_bb(k,c) r0_1 t2_bbbb(c,b,j,k) t1_1_aa(a,i)
    //                += +0.25 <l,k||c,d>_abab r2_1_abab(c,d,i,j) t2_abab(a,b,l,k)
    //                += +0.25 <k,l||c,d>_abab r2_1_abab(c,d,i,j) t2_abab(a,b,k,l)
    //                += +0.25 <l,k||d,c>_abab r2_1_abab(d,c,i,j) t2_abab(a,b,l,k)
    //                += +0.25 <k,l||d,c>_abab r2_1_abab(d,c,i,j) t2_abab(a,b,k,l)
    //                += +1.00 d-_bb(k,j) r1_1_bb(b,k) t1_1_aa(a,i)
    //                += +1.00 d-_aa(k,c) r1_aa(c,i) t0_1 t2_1_abab(a,b,k,j)
    //                += +1.00 d-_aa(k,c) r0_1 t2_aaaa(c,a,i,k) t1_1_bb(b,j)
    //                += +1.00 d-_bb(k,c) r1_bb(c,j) t1_1_aa(a,i) t1_1_bb(b,k)
    //                += +1.00 d-_aa(k,c) r1_aa(c,i) t1_1_aa(a,k) t1_1_bb(b,j)
    //                += -1.00 d-_aa(a,c) r1_1_aa(c,i) t1_1_bb(b,j)
    //                += +1.00 d-_aa(k,c) r1_aa(a,k) t1_1_bb(b,j) t1_1_aa(c,i)
    sigmar2_1_abab("R,a,b,i,j") -= tmps_["296_bbaa_Lvovo"]("R,b,j,a,i");
    tmps_["296_bbaa_Lvovo"].~TArrayD();
}
