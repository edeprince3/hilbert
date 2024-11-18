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

#include "cc_cavity/include/qed_ccsd_21/eom_ea_qed_ccsd_21.h"
#include "cc_cavity/include/qed_ccsd_21/qed_ccsd_21.h"

void hilbert::EOM_EA_QED_CCSD::sigma_ea_21_4() {

    // Get cavity information
    double w0 = cc_wfn_->cavity_frequency_[2];
    double coupling_factor_z = w0 * cc_wfn_->cavity_coupling_strength_[2];

    // get effective dipole integrals
    TArrayMap dp = reinterpret_pointer_cast<QED_CCSD>(cc_wfn_)->effective_dipole();

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

    /// reference right operators
    TArrayMap r1 = {
            {"a", evec_blks_["r1_a"]},
            {"b", evec_blks_["r1_b"]}
    };
    TArrayMap r2 = {
            {"aaa", evec_blks_["r2_aaa"]},
            {"abb", evec_blks_["r2_abb"]},
    };
    TArrayMap r1_1 = {
            {"a", evec_blks_["r1_1_a"]},
            {"b", evec_blks_["r1_1_b"]}
    };
    TArrayMap r2_1 = {
            {"aaa", evec_blks_["r2_1_aaa"]},
            {"abb", evec_blks_["r2_1_abb"]},
    };

    /// reference left operators

    TArrayMap l1 = {
            {"a", evec_blks_["l1_a"]},
            {"b", evec_blks_["l1_b"]}
    };
    TArrayMap l2 = {
            {"aaa", evec_blks_["l2_aaa"]},
            {"bab", evec_blks_["l2_bab"]},
    };
    TArrayMap l1_1 = {
            {"a", evec_blks_["l1_1_a"]},
            {"b", evec_blks_["l1_1_b"]}
    };
    TArrayMap l2_1 = {
            {"aaa", evec_blks_["l2_1_aaa"]},
            {"bab", evec_blks_["l2_1_bab"]},
    };

    /// reference sigma vectors

    // electronic

    auto &sigmar1_a = sigvec_blks_["sigmar1_a"];
    auto &sigmal1_a = sigvec_blks_["sigmal1_a"];

    auto &sigmar2_aaa = sigvec_blks_["sigmar2_aaa"];
    auto &sigmar2_abb = sigvec_blks_["sigmar2_abb"];
    auto &sigmal2_aaa = sigvec_blks_["sigmal2_aaa"];
    auto &sigmal2_abb = sigvec_blks_["sigmal2_abb"];

    // photonic

    auto &sigmar1_1_a = sigvec_blks_["sigmar1_1_a"];
    auto &sigmal1_1_a = sigvec_blks_["sigmal1_1_a"];

    auto &sigmar2_1_aaa = sigvec_blks_["sigmar2_1_aaa"];
    auto &sigmar2_1_abb = sigvec_blks_["sigmar2_1_abb"];
    auto &sigmal2_1_aaa = sigvec_blks_["sigmal2_1_aaa"];
    auto &sigmal2_1_abb = sigvec_blks_["sigmal2_1_abb"];


    /// extract coherent state basis operators

    // electronic

    auto &csigmar1_a = tmps_["csigmar1_a"];
    auto &csigmal1_a = tmps_["csigmal1_a"];

    auto &csigmar2_aaa = tmps_["csigmar2_aaa"];
    auto &csigmar2_abb = tmps_["csigmar2_abb"];
    auto &csigmal2_aaa = tmps_["csigmal2_aaa"];
    auto &csigmal2_abb = tmps_["csigmal2_abb"];

    // photonic

    auto &csigmar1_1_a = tmps_["csigmar1_1_a"];
    auto &csigmal1_1_a = tmps_["csigmal1_1_a"];

    auto &csigmar2_1_aaa = tmps_["csigmar2_1_aaa"];
    auto &csigmar2_1_abb = tmps_["csigmar2_1_abb"];
    auto &csigmal2_1_aaa = tmps_["csigmal2_1_aaa"];
    auto &csigmal2_1_abb = tmps_["csigmal2_1_abb"];


    // sigmal2_1_abb += +1.00 d-_bb(i,b) t1_1_aa(c,j) l2_1_aaa(j,c,a)
    //               += -0.50 <k,i||a,b>_abab t2_abab(d,c,k,j) l2_1_bab(j,d,c)
    //               += -0.50 <k,i||a,b>_abab t2_abab(c,d,k,j) l2_1_bab(j,c,d)
    //               += -1.00 f_bb(i,j) l2_1_bab(j,a,b)
    //               += -2.00 d-_aa(j,c) t1_1_aa(c,j) l2_1_bab(i,a,b)
    //               += -1.00 d-_aa(j,j) l2_bab(i,a,b)
    //               += -1.00 <c,i||j,b>_abab l2_1_aaa(j,c,a)
    //               += +1.00 <i,k||d,b>_bbbb t2_abab(c,d,j,k) l2_1_aaa(j,c,a)
    //               += -1.00 d-_bb(i,b) t0_1 l1_1_a(a)
    //               += +1.00 <k,i||d,b>_abab t2_aaaa(d,c,j,k) l2_1_aaa(j,c,a)
    //               += +1.00 <i,k||d,b>_bbbb t2_bbbb(d,c,j,k) l2_1_bab(j,a,c)
    //               += -0.50 <k,j||d,a>_aaaa t2_aaaa(d,c,k,j) l2_1_bab(i,c,b)
    //               += -1.00 d-_bb(c,b) l2_bab(i,a,c)
    //               += -1.00 d-_bb(c,b) t0_1 l2_1_bab(i,a,c)
    //               += -0.50 <k,j||a,d>_abab t2_abab(c,d,k,j) l2_1_bab(i,c,b)
    //               += -0.50 <j,k||a,d>_abab t2_abab(c,d,j,k) l2_1_bab(i,c,b)
    //               += -0.50 <k,j||d,b>_abab t2_abab(d,c,k,j) l2_1_bab(i,a,c)
    //               += -0.50 <j,k||d,b>_abab t2_abab(d,c,j,k) l2_1_bab(i,a,c)
    //               += -2.00 d-_bb(j,c) t1_1_bb(c,j) l2_1_bab(i,a,b)
    //               += +1.00 f_bb(i,b) l1_1_a(a)
    //               += +1.00 <k,i||d,b>_abab t2_abab(d,c,k,j) l2_1_bab(j,a,c)
    //               += -1.00 d-_aa(c,a) t0_1 l2_1_bab(i,c,b)
    //               += -0.50 <i,k||c,d>_bbbb t2_bbbb(c,d,j,k) l2_1_bab(j,a,b)
    //               += +1.00 f_bb(c,b) l2_1_bab(i,a,c)
    //               += -1.00 d-_bb(j,j) l2_bab(i,a,b)
    //               += +1.00 <i,c||b,j>_bbbb l2_1_bab(j,a,c)
    //               += +1.00 l2_1_bab(i,a,b) w0
    //               += -1.00 d-_bb(i,b) l1_a(a)
    //               += +0.50 <d,c||a,b>_abab l2_1_bab(i,d,c)
    //               += +0.50 <c,d||a,b>_abab l2_1_bab(i,c,d)
    //               += +1.00 d-_bb(i,j) l2_bab(j,a,b)
    //               += -0.50 <k,i||c,d>_abab t2_abab(c,d,k,j) l2_1_bab(j,a,b)
    //               += -0.50 <k,i||d,c>_abab t2_abab(d,c,k,j) l2_1_bab(j,a,b)
    //               += -0.50 <k,j||k,j>_abab l2_1_bab(i,a,b)
    //               += -0.50 <j,k||j,k>_abab l2_1_bab(i,a,b)
    //               += +0.25 <k,j||c,d>_bbbb t2_bbbb(c,d,k,j) l2_1_bab(i,a,b)
    //               += +1.00 f_bb(j,j) l2_1_bab(i,a,b)
    //               += +0.25 <k,j||c,d>_aaaa t2_aaaa(c,d,k,j) l2_1_bab(i,a,b)
    //               += -1.00 d-_aa(j,j) t0_1 l2_1_bab(i,a,b)
    //               += +1.00 f_aa(j,j) l2_1_bab(i,a,b)
    //               += -0.50 <k,j||k,j>_bbbb l2_1_bab(i,a,b)
    //               += -0.50 <k,j||k,j>_aaaa l2_1_bab(i,a,b)
    //               += +0.25 <k,j||c,d>_abab t2_abab(c,d,k,j) l2_1_bab(i,a,b)
    //               += +0.25 <j,k||c,d>_abab t2_abab(c,d,j,k) l2_1_bab(i,a,b)
    //               += +0.25 <k,j||d,c>_abab t2_abab(d,c,k,j) l2_1_bab(i,a,b)
    //               += +0.25 <j,k||d,c>_abab t2_abab(d,c,j,k) l2_1_bab(i,a,b)
    //               += -1.00 d-_bb(j,j) t0_1 l2_1_bab(i,a,b)
    //               += +1.00 d-_bb(i,j) t0_1 l2_1_bab(j,a,b)
    //               += +1.00 f_aa(c,a) l2_1_bab(i,c,b)
    //               += +1.00 <c,i||a,b>_abab l1_1_a(c)
    //               += -1.00 d-_aa(c,a) l2_bab(i,c,b)
    //               += -0.50 <k,j||d,b>_bbbb t2_bbbb(d,c,k,j) l2_1_bab(i,a,c)
    //               += +0.25 <k,j||a,b>_abab t2_abab(d,c,k,j) l2_1_bab(i,d,c)
    //               += +0.25 <j,k||a,b>_abab t2_abab(d,c,j,k) l2_1_bab(i,d,c)
    //               += +0.25 <k,j||a,b>_abab t2_abab(c,d,k,j) l2_1_bab(i,c,d)
    //               += +0.25 <j,k||a,b>_abab t2_abab(c,d,j,k) l2_1_bab(i,c,d)
    //               += +2.00 d-_aa(j,a) t1_1_aa(c,j) l2_1_bab(i,c,b)
    //               += +0.50 <k,i||a,b>_abab t2_aaaa(d,c,j,k) l2_1_aaa(j,d,c)
    //               += +2.00 d-_bb(j,b) t1_1_bb(c,j) l2_1_bab(i,a,c)
    //               += +1.00 <k,i||a,d>_abab t2_abab(c,d,k,j) l2_1_bab(j,c,b)
    //               += -1.00 d-_bb(i,b) t1_1_bb(c,j) l2_1_bab(j,a,c)
    sigmal2_1_abb("L,a,b,i") += tmps_["98_bba_Lovv"]("L,i,b,a");
    tmps_["98_bba_Lovv"].~TArrayD();

    // flops: o0v1L1  = o1v2L1
    //  mems: o0v1L1  = o0v1L1
    tmps_["99_a_Lv"]("R,a")  = -1.00 * dp["bb_ov"]("i,b") * r2_1["abb"]("R,a,b,i");
    tmps_["99_a_Lv"].~TArrayD();

    // flops: o0v1L1  = o0v2L1
    //  mems: o0v1L1  = o0v1L1
    tmps_["100_a_Lv"]("R,a")  = -1.00 * dp["aa_vv"]("a,b") * r1_1["a"]("R,b");
    tmps_["100_a_Lv"].~TArrayD();

    // flops: o1v2L1  = o2v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["101_abb_Lvvo"]("R,a,b,i")  = r2["abb"]("R,a,b,j") * reused_["59_bb_oo"]("j,i");

    // sigmar2_abb += +1.00 d-_bb(j,c) r2_abb(a,b,j) t1_1_bb(c,i)
    sigmar2_abb("R,a,b,i") += tmps_["101_abb_Lvvo"]("R,a,b,i");

    // flops: o0v1L1  = o1v2L1 o1v2L1 o0v1L1
    //  mems: o0v1L1  = o0v1L1 o0v1L1 o0v1L1
    tmps_["110_a_Lv"]("R,a")  = -1.00 * dp["bb_ov"]("j,c") * r2_1["abb"]("R,a,c,j");
    tmps_["110_a_Lv"]("R,a") += dp["aa_ov"]("i,b") * r2_1["aaa"]("R,b,a,i");

    // sigmar1_a += +1.00 d-_aa(i,b) r2_1_aaa(b,a,i)
    //           += -1.00 d-_bb(i,b) r2_1_abb(a,b,i)
    sigmar1_a("R,a") += tmps_["110_a_Lv"]("R,a");

    // flops: o0v1L1  = o0v2L1 o1v1L1 o0v1L1
    //  mems: o0v1L1  = o0v1L1 o0v1L1 o0v1L1
    tmps_["111_a_Lv"]("R,a")  = -1.00 * dp["aa_vv"]("a,b") * r1_1["a"]("R,b");
    tmps_["111_a_Lv"]("R,a") += t1_1["aa"]("a,i") * tmps_["76_a_Lo"]("R,i");
    tmps_["76_a_Lo"].~TArrayD();

    // sigmar1_a += +1.00 d-_aa(i,b) r1_a(b) t1_1_aa(a,i)
    //           += -1.00 d-_aa(a,b) r1_1_a(b)
    sigmar1_a("R,a") += tmps_["111_a_Lv"]("R,a");

    // flops: o2v1L1  = o2v2L1 o2v3L1 o2v2L1 o2v1L1 o2v1L1 o2v2L1 o2v1L1 o3v2L1 o2v1L1 o2v2L1 o2v1L1
    //  mems: o2v1L1  = o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1
    tmps_["107_bba_Lvoo"]("R,a,i,j")  = 0.50 * eri["baab_vovo"]("a,j,b,i") * r1["a"]("R,b");
    tmps_["107_bba_Lvoo"]("R,a,i,j") += 0.50 * eri["baab_vovv"]("a,j,b,c") * r2["abb"]("R,b,c,i");
    tmps_["107_bba_Lvoo"]("R,a,i,j") += dp["aa_ov"]("j,b") * r2_1["abb"]("R,b,a,i");
    tmps_["107_bba_Lvoo"]("R,a,i,j") -= 0.50 * r2["abb"]("R,b,a,i") * f["aa_ov"]("j,b");
    tmps_["107_bba_Lvoo"]("R,a,i,j") += 0.50 * r2["abb"]("R,b,a,k") * eri["abab_oovo"]("j,k,b,i");
    tmps_["107_bba_Lvoo"]("R,a,i,j") += 0.50 * r2["abb"]("R,b,a,i") * reused_["4_aa_ov"]("j,b");

    // flops: o1v2L1  = o2v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["109_abb_Lvvo"]("R,a,b,i")  = -4.00 * t1_1["aa"]("a,j") * tmps_["107_bba_Lvoo"]("R,b,i,j");
    tmps_["107_bba_Lvoo"].~TArrayD();

    // flops: o2v1L1  = o2v2L1 o2v2L1 o2v1L1 o2v3L1 o2v1L1 o2v2L1 o2v1L1 o3v2L1 o2v1L1 o3v2L1 o2v1L1 o2v2L1 o2v1L1
    //  mems: o2v1L1  = o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1
    tmps_["108_abb_Lvoo"]("R,a,i,j")  = -0.50 * r2["abb"]("R,a,b,i") * f["bb_ov"]("j,b");
    tmps_["108_abb_Lvoo"]("R,a,i,j") += r2_1["abb"]("R,a,b,i") * dp["bb_ov"]("j,b");
    tmps_["108_abb_Lvoo"]("R,a,i,j") -= 0.50 * eri["abab_vovv"]("a,j,c,d") * r2["abb"]("R,c,d,i");
    tmps_["108_abb_Lvoo"]("R,a,i,j") -= 0.50 * eri["abab_vovo"]("a,j,c,i") * r1["a"]("R,c");
    tmps_["108_abb_Lvoo"]("R,a,i,j") += 0.50 * r2["abb"]("R,a,b,l") * eri["bbbb_oovo"]("j,l,b,i");
    tmps_["108_abb_Lvoo"]("R,a,i,j") += 0.50 * r2["aaa"]("R,c,a,k") * eri["abab_oovo"]("k,j,c,i");
    tmps_["108_abb_Lvoo"]("R,a,i,j") += 0.50 * r2["abb"]("R,a,b,i") * reused_["5_bb_ov"]("j,b");

    // flops: o2v1L1  = o2v2 o2v2L1
    //  mems: o2v1L1  = o2v2 o2v1L1
    tmps_["106_bba_Lvoo"]("R,a,i,j")  = (reused_["48_bbaa_voov"]("a,i,j,b") + -1.00 * reused_["50_bbaa_voov"]("a,i,j,b")) * r1["a"]("R,b");

    // flops: o1v2L1  = o2v3L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["105_bab_Lvov"]("R,a,i,b")  = r2_1["abb"]("R,c,a,j") * eri["abab_oovv"]("i,j,c,b");

    // flops: o2v1L1  = o2v2L1
    //  mems: o2v1L1  = o2v1L1
    tmps_["104_abb_Lvoo"]("R,a,i,j")  = r1["a"]("R,b") * reused_["45_abba_voov"]("a,i,j,b");

    // flops: o3v0L1  = o3v2L1
    //  mems: o3v0L1  = o3v0L1
    tmps_["103_abb_Looo"]("R,i,j,k")  = eri["abab_oovv"]("i,j,a,b") * r2_1["abb"]("R,a,b,k");

    // flops: o1v0L1  = o1v1L1
    //  mems: o1v0L1  = o1v0L1
    tmps_["102_a_Lo"]("R,i")  = r1_1["a"]("R,a") * reused_["4_aa_ov"]("i,a");

    // flops: o1v2L1  = o2v2L1 o2v2L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o3v2L1 o1v2L1 o1v3L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v4L1 o1v2L1 o2v3L1 o1v2L1 o2v2L1 o1v2L1 o1v3L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o3v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["112_bab_Lvvo"]("R,a,b,i")  = -0.50 * t2_1["abab"]("b,a,k,i") * tmps_["26_a_Lo"]("R,k");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= t1_1["aa"]("b,m") * tmps_["106_bba_Lvoo"]("R,a,i,m");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += t2_1["abab"]("b,d,k,i") * tmps_["80_bab_Lvov"]("R,a,k,d");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += tmps_["111_a_Lv"]("R,b") * t1_1["bb"]("a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += t2["abab"]("b,d,k,i") * tmps_["105_bab_Lvov"]("R,a,k,d");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += t2_1["abab"]("b,a,k,i") * tmps_["25_a_Lo"]("R,k");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += tmps_["110_a_Lv"]("R,b") * t1_1["bb"]("a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += 2.00 * t2_1["abab"]("b,a,k,i") * tmps_["24_a_Lo"]("R,k");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += t2["abab"]("b,a,m,l") * tmps_["103_abb_Looo"]("R,m,l,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2["abb"]("R,b,d,i") * dp["bb_vv"]("a,d");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2["abb"]("R,b,f,i") * reused_["94_bb_vv"]("a,f");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += r1_1["a"]("R,e") * reused_["30_abab_vovv"]("e,i,b,a");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += r1_1["a"]("R,c") * reused_["29_baab_vvvo"]("a,c,b,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2_1["aaa"]("R,c,b,m") * reused_["48_bbaa_voov"]("a,i,m,c");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += 0.50 * r2_1["abb"]("R,c,a,i") * reused_["1_aa_vv"]("b,c");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += r1["a"]("R,e") * reused_["46_abab_vovv"]("e,i,b,a");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2["abb"]("R,b,a,l") * reused_["87_bb_oo"]("i,l");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += r2_1["abb"]("R,b,a,l") * reused_["67_bb_oo"]("l,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += 0.50 * r2_1["abb"]("R,b,f,i") * reused_["56_bb_vv"]("a,f");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2_1["abb"]("R,c,a,i") * reused_["2_aa_vv"]("b,c");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2_1["abb"]("R,b,a,l") * f["bb_oo"]("l,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= scalars_["5"] * r2_1["abb"]("R,b,a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= dp["aa_vv"]("b,e") * r2["abb"]("R,e,a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += r1_1["a"]("R,b") * f["bb_vo"]("a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += w0 * r2_1["abb"]("R,b,a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2_1["abb"]("R,b,d,l") * eri["bbbb_vovo"]("a,l,d,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2["abb"]("R,b,f,j") * reused_["92_bbbb_voov"]("a,i,j,f");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += r2["abb"]("R,b,f,j") * reused_["80_bbbb_voov"]("a,i,j,f");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += r1_1["a"]("R,b") * reused_["37_bb_vo"]("a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2["abb"]("R,c,a,i") * reused_["14_aa_vv"]("b,c");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += eri["abab_vvvv"]("b,a,e,f") * r2_1["abb"]("R,e,f,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += r2_1["aaa"]("R,c,b,m") * reused_["50_bbaa_voov"]("a,i,m,c");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += r2["abb"]("R,b,a,l") * dp["bb_oo"]("l,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2_1["abb"]("R,b,f,i") * reused_["55_bb_vv"]("a,f");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2_1["abb"]("R,b,f,j") * reused_["53_bbbb_voov"]("a,i,j,f");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2["abb"]("R,c,a,l") * reused_["83_abab_vovo"]("b,l,c,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += f["aa_vv"]("b,e") * r2_1["abb"]("R,e,a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += 0.50 * r2_1["abb"]("R,b,a,j") * reused_["58_bb_oo"]("i,j");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2["aaa"]("R,c,b,m") * reused_["91_bbaa_voov"]("a,i,m,c");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += r2_1["aaa"]("R,e,b,k") * eri["baab_vovo"]("a,k,e,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= 2.00 * scalars_["4"] * r2_1["abb"]("R,b,a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += r2["aaa"]("R,c,b,m") * reused_["72_bbaa_voov"]("a,i,m,c");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r1_1["a"]("R,c") * reused_["28_aabb_vvvo"]("b,c,a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= eri["abab_vovo"]("b,l,e,i") * r2_1["abb"]("R,e,a,l");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2["abb"]("R,b,f,i") * reused_["88_bb_vv"]("a,f");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= scalars_["1"] * r2["abb"]("R,b,a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += r2["aaa"]("R,c,b,k") * reused_["73_baab_vovo"]("a,k,c,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= 2.00 * scalars_["3"] * r2_1["abb"]("R,b,a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += r2_1["abb"]("R,b,f,j") * reused_["51_bbbb_voov"]("a,i,j,f");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += r2["abb"]("R,b,f,l") * reused_["93_bbbb_vovo"]("a,l,f,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += 0.50 * r1["a"]("R,b") * reused_["63_bb_ov"]("i,a");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2["abb"]("R,c,a,i") * reused_["13_aa_vv"]("b,c");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += 2.00 * r1_1["a"]("R,b") * reused_["35_bb_ov"]("i,a");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += r1["a"]("R,c") * reused_["60_abab_vovv"]("c,i,b,a");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += 2.00 * r2_1["abb"]("R,b,a,l") * reused_["59_bb_oo"]("l,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += scalars_["6"] * r2["abb"]("R,b,a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2_1["abb"]("R,b,d,i") * reused_["66_bb_vv"]("a,d");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r1["a"]("R,b") * dp["bb_vo"]("a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += r1["a"]("R,c") * reused_["43_abba_vovv"]("b,i,a,c");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2_1["abb"]("R,b,a,j") * reused_["57_bb_oo"]("j,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2["abb"]("R,b,a,j") * reused_["89_bb_oo"]("j,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += r2_1["abb"]("R,b,d,i") * f["bb_vv"]("a,d");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += r1["a"]("R,b") * reused_["33_bb_vo"]("a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += eri["abab_vvvo"]("b,a,e,i") * r1_1["a"]("R,e");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2_1["abb"]("R,e,a,i") * reused_["3_aa_vv"]("b,e");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r1["a"]("R,c") * reused_["90_aabb_vvvo"]("b,c,a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r1_1["a"]("R,c") * reused_["65_aabb_vvvo"]("b,c,a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= scalars_["2"] * r2["abb"]("R,b,a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= r2["abb"]("R,b,a,j") * reused_["95_bb_oo"]("j,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += t2_1["abab"]("b,a,m,l") * tmps_["81_abb_Looo"]("R,m,l,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += t2["abab"]("b,a,k,i") * tmps_["47_a_Lo"]("R,k");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += 2.00 * tmps_["108_abb_Lvoo"]("R,b,i,l") * t1_1["bb"]("a,l");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += t0_1 * tmps_["101_abb_Lvvo"]("R,b,a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += t2["abab"]("b,a,k,i") * tmps_["102_a_Lo"]("R,k");
    tmps_["112_bab_Lvvo"]("R,a,b,i") -= 0.50 * tmps_["109_abb_Lvvo"]("R,b,a,i");
    tmps_["112_bab_Lvvo"]("R,a,b,i") += t1_1["bb"]("a,j") * tmps_["104_abb_Lvoo"]("R,b,i,j");
    tmps_["109_abb_Lvvo"].~TArrayD();
    tmps_["108_abb_Lvoo"].~TArrayD();
    tmps_["106_bba_Lvoo"].~TArrayD();
    tmps_["105_bab_Lvov"].~TArrayD();
    tmps_["104_abb_Lvoo"].~TArrayD();
    tmps_["103_abb_Looo"].~TArrayD();
    tmps_["102_a_Lo"].~TArrayD();
    tmps_["101_abb_Lvvo"].~TArrayD();
    tmps_["81_abb_Looo"].~TArrayD();
    tmps_["80_bab_Lvov"].~TArrayD();
    tmps_["47_a_Lo"].~TArrayD();
    tmps_["26_a_Lo"].~TArrayD();
    tmps_["25_a_Lo"].~TArrayD();
    tmps_["24_a_Lo"].~TArrayD();

    // sigmar2_1_abb += +1.00 <j,k||d,c>_abab r1_a(d) t2_abab(a,c,j,i) t1_1_bb(b,k)
    //               += +0.50 <k,j||c,d>_aaaa r2_aaa(c,d,k) t2_1_abab(a,b,j,i)
    //               += -0.50 <j,k||c,d>_abab r2_abb(c,d,k) t2_1_abab(a,b,j,i)
    //               += -0.50 <j,k||d,c>_abab r2_abb(d,c,k) t2_1_abab(a,b,j,i)
    //               += -1.00 f_aa(j,c) r1_a(c) t2_1_abab(a,b,j,i)
    //               += +1.00 <k,j||c,d>_aaaa r1_a(d) t2_abab(c,b,j,i) t1_1_aa(a,k)
    //               += +1.00 <k,j||d,c>_abab r1_a(d) t2_bbbb(c,b,i,j) t1_1_aa(a,k)
    //               += +1.00 <j,k||d,c>_abab r2_abb(d,b,k) t2_1_abab(a,c,j,i)
    //               += +1.00 d-_aa(j,c) r1_a(c) t1_1_aa(a,j) t1_1_bb(b,i)
    //               += -1.00 d-_aa(a,c) r1_1_a(c) t1_1_bb(b,i)
    //               += +1.00 <j,k||d,c>_abab r2_1_abb(d,b,k) t2_abab(a,c,j,i)
    //               += +1.00 d-_aa(j,c) r1_a(c) t0_1 t2_1_abab(a,b,j,i)
    //               += +1.00 d-_aa(j,c) r2_1_aaa(c,a,j) t1_1_bb(b,i)
    //               += -1.00 d-_bb(j,c) r2_1_abb(a,c,j) t1_1_bb(b,i)
    //               += +2.00 d-_aa(j,c) r1_1_a(c) t2_1_abab(a,b,j,i)
    //               += +0.25 <k,j||c,d>_abab r2_1_abb(c,d,i) t2_abab(a,b,k,j)
    //               += +0.25 <j,k||c,d>_abab r2_1_abb(c,d,i) t2_abab(a,b,j,k)
    //               += +0.25 <k,j||d,c>_abab r2_1_abb(d,c,i) t2_abab(a,b,k,j)
    //               += +0.25 <j,k||d,c>_abab r2_1_abb(d,c,i) t2_abab(a,b,j,k)
    //               += -1.00 d+_bb(b,c) r2_abb(a,c,i)
    //               += +1.00 <j,b||c,d>_bbbb r2_abb(a,d,i) t1_1_bb(c,j)
    //               += -0.50 <k,j||c,d>_bbbb r2_abb(a,d,i) t2_1_bbbb(c,b,k,j)
    //               += +0.50 <k,j||c,i>_abab r1_1_a(c) t2_abab(a,b,k,j)
    //               += +0.50 <j,k||c,i>_abab r1_1_a(c) t2_abab(a,b,j,k)
    //               += -1.00 <j,b||d,c>_abab r1_1_a(d) t2_abab(a,c,j,i)
    //               += +1.00 <k,j||c,d>_aaaa r2_1_aaa(d,a,k) t2_abab(c,b,j,i)
    //               += -0.50 <k,j||c,d>_aaaa r2_1_abb(d,b,i) t2_aaaa(c,a,k,j)
    //               += +0.50 <k,j||c,i>_abab r1_a(c) t2_1_abab(a,b,k,j)
    //               += +0.50 <j,k||c,i>_abab r1_a(c) t2_1_abab(a,b,j,k)
    //               += -1.00 f_bb(j,c) r2_abb(a,b,j) t1_1_bb(c,i)
    //               += +1.00 d-_bb(j,i) r2_1_abb(a,b,j) t0_1
    //               += -0.50 <k,j||c,d>_bbbb r2_1_abb(a,d,i) t2_bbbb(c,b,k,j)
    //               += -0.50 <k,j||d,c>_abab r2_1_abb(d,b,i) t2_abab(a,c,k,j)
    //               += -0.50 <j,k||d,c>_abab r2_1_abb(d,b,i) t2_abab(a,c,j,k)
    //               += -1.00 f_bb(j,i) r2_1_abb(a,b,j)
    //               += -0.50 <k,j||k,j>_abab r2_1_abb(a,b,i)
    //               += -0.50 <j,k||j,k>_abab r2_1_abb(a,b,i)
    //               += +0.25 <k,j||c,d>_bbbb r2_1_abb(a,b,i) t2_bbbb(c,d,k,j)
    //               += +1.00 f_bb(j,j) r2_1_abb(a,b,i)
    //               += +0.25 <k,j||c,d>_aaaa r2_1_abb(a,b,i) t2_aaaa(c,d,k,j)
    //               += -1.00 d-_aa(j,j) r2_1_abb(a,b,i) t0_1
    //               += +1.00 f_aa(j,j) r2_1_abb(a,b,i)
    //               += -0.50 <k,j||k,j>_bbbb r2_1_abb(a,b,i)
    //               += -0.50 <k,j||k,j>_aaaa r2_1_abb(a,b,i)
    //               += +0.25 <k,j||c,d>_abab r2_1_abb(a,b,i) t2_abab(c,d,k,j)
    //               += +0.25 <j,k||c,d>_abab r2_1_abb(a,b,i) t2_abab(c,d,j,k)
    //               += +0.25 <k,j||d,c>_abab r2_1_abb(a,b,i) t2_abab(d,c,k,j)
    //               += +0.25 <j,k||d,c>_abab r2_1_abb(a,b,i) t2_abab(d,c,j,k)
    //               += -1.00 d-_bb(j,j) r2_1_abb(a,b,i) t0_1
    //               += -1.00 d+_aa(a,c) r2_abb(c,b,i)
    //               += +1.00 f_bb(b,i) r1_1_a(a)
    //               += +1.00 r2_1_abb(a,b,i) w0
    //               += +1.00 <j,b||c,i>_bbbb r2_1_abb(a,c,j)
    //               += +1.00 <k,j||c,d>_bbbb r2_abb(a,d,k) t2_1_bbbb(c,b,i,j)
    //               += +1.00 <j,k||c,d>_abab r2_abb(a,d,k) t2_1_abab(c,b,j,i)
    //               += +1.00 d-_bb(j,c) r1_1_a(a) t2_bbbb(c,b,i,j) t0_1
    //               += -1.00 d-_aa(j,c) r1_1_a(a) t2_abab(c,b,j,i) t0_1
    //               += -1.00 d-_bb(b,i) r1_1_a(a) t0_1
    //               += -1.00 d-_aa(j,j) r1_1_a(a) t1_1_bb(b,i)
    //               += -1.00 d-_bb(j,j) r1_1_a(a) t1_1_bb(b,i)
    //               += -0.50 <k,j||c,i>_abab r1_1_a(a) t2_abab(c,b,k,j)
    //               += -0.50 <j,k||c,i>_abab r1_1_a(a) t2_abab(c,b,j,k)
    //               += -0.50 <j,b||c,d>_bbbb r1_1_a(a) t2_bbbb(c,d,i,j)
    //               += +1.00 f_aa(j,c) r1_1_a(a) t2_abab(c,b,j,i)
    //               += +0.50 <j,b||c,d>_abab r1_1_a(a) t2_abab(c,d,j,i)
    //               += +0.50 <j,b||d,c>_abab r1_1_a(a) t2_abab(d,c,j,i)
    //               += -0.50 <k,j||c,i>_bbbb r1_1_a(a) t2_bbbb(c,b,k,j)
    //               += -1.00 f_bb(j,c) r1_1_a(a) t2_bbbb(c,b,i,j)
    //               += +1.00 <j,a||c,d>_aaaa r2_abb(d,b,i) t1_1_aa(c,j)
    //               += -0.50 <k,j||c,d>_aaaa r2_abb(d,b,i) t2_1_aaaa(c,a,k,j)
    //               += +0.50 <a,b||c,d>_abab r2_1_abb(c,d,i)
    //               += +0.50 <a,b||d,c>_abab r2_1_abb(d,c,i)
    //               += +1.00 <k,j||d,c>_abab r2_1_aaa(d,a,k) t2_bbbb(c,b,i,j)
    //               += +1.00 d+_bb(j,i) r2_abb(a,b,j)
    //               += -0.50 <k,j||c,d>_abab r2_1_abb(a,d,i) t2_abab(c,b,k,j)
    //               += -0.50 <j,k||c,d>_abab r2_1_abb(a,d,i) t2_abab(c,b,j,k)
    //               += +1.00 <k,j||c,d>_bbbb r2_1_abb(a,d,k) t2_bbbb(c,b,i,j)
    //               += -1.00 <a,j||d,c>_abab r2_abb(d,b,j) t1_1_bb(c,i)
    //               += +1.00 f_aa(a,c) r2_1_abb(c,b,i)
    //               += -0.50 <k,j||c,d>_bbbb r2_1_abb(a,b,k) t2_bbbb(c,d,i,j)
    //               += +1.00 <k,j||c,d>_aaaa r2_aaa(d,a,k) t2_1_abab(c,b,j,i)
    //               += -1.00 <j,b||c,i>_abab r2_1_aaa(c,a,j)
    //               += -2.00 d-_aa(j,c) r2_1_abb(a,b,i) t1_1_aa(c,j)
    //               += +1.00 <k,j||d,c>_abab r2_aaa(d,a,k) t2_1_bbbb(c,b,i,j)
    //               += -1.00 <a,j||d,c>_abab r1_1_a(d) t2_bbbb(c,b,i,j)
    //               += -1.00 <a,j||c,i>_abab r2_1_abb(c,b,j)
    //               += -0.50 <k,j||c,d>_abab r2_abb(a,d,i) t2_1_abab(c,b,k,j)
    //               += -0.50 <j,k||c,d>_abab r2_abb(a,d,i) t2_1_abab(c,b,j,k)
    //               += +1.00 <j,b||c,d>_abab r2_abb(a,d,i) t1_1_aa(c,j)
    //               += -1.00 d+_bb(j,j) r2_abb(a,b,i)
    //               += -1.00 <j,b||d,c>_abab r2_aaa(d,a,j) t1_1_bb(c,i)
    //               += -2.00 d-_bb(j,c) r2_1_abb(a,b,i) t1_1_bb(c,j)
    //               += +1.00 <j,k||c,d>_abab r2_1_abb(a,d,k) t2_abab(c,b,j,i)
    //               += -1.00 <j,b||c,d>_bbbb r2_abb(a,d,j) t1_1_bb(c,i)
    //               += -0.50 <k,j||c,d>_bbbb r1_a(a) t2_bbbb(c,d,i,j) t1_1_bb(b,k)
    //               += +1.00 r1_a(a) t1_1_bb(b,i) w0
    //               += -0.50 <k,j||c,i>_abab r1_a(a) t2_1_abab(c,b,k,j)
    //               += -0.50 <j,k||c,i>_abab r1_a(a) t2_1_abab(c,b,j,k)
    //               += -1.00 d-_bb(j,c) r1_a(a) t1_1_bb(b,i) t1_1_bb(c,j)
    //               += +1.00 f_aa(j,c) r1_a(a) t2_1_abab(c,b,j,i)
    //               += -1.00 f_bb(j,i) r1_a(a) t1_1_bb(b,j)
    //               += +1.00 <j,b||c,i>_bbbb r1_a(a) t1_1_bb(c,j)
    //               += -1.00 d-_aa(j,c) r1_a(a) t1_1_bb(b,i) t1_1_aa(c,j)
    //               += -1.00 f_bb(j,c) r1_a(a) t2_1_bbbb(c,b,i,j)
    //               += -0.50 <j,b||c,d>_bbbb r1_a(a) t2_1_bbbb(c,d,i,j)
    //               += +1.00 f_bb(b,c) r1_a(a) t1_1_bb(c,i)
    //               += -0.50 <k,j||c,i>_bbbb r1_a(a) t2_1_bbbb(c,b,k,j)
    //               += +0.50 <j,b||c,d>_abab r1_a(a) t2_1_abab(c,d,j,i)
    //               += +0.50 <j,b||d,c>_abab r1_a(a) t2_1_abab(d,c,j,i)
    //               += +1.00 <j,b||c,i>_abab r1_a(a) t1_1_aa(c,j)
    //               += -0.50 <k,j||c,d>_bbbb r1_a(a) t2_bbbb(c,b,k,j) t1_1_bb(d,i)
    //               += -1.00 <k,j||c,d>_aaaa r1_a(a) t2_abab(c,b,j,i) t1_1_aa(d,k)
    //               += +1.00 <k,j||c,d>_bbbb r1_a(a) t2_bbbb(c,b,i,j) t1_1_bb(d,k)
    //               += +2.00 d-_bb(j,c) r1_a(a) t1_1_bb(b,j) t1_1_bb(c,i)
    //               += -1.00 <k,j||d,c>_abab r1_a(a) t2_bbbb(c,b,i,j) t1_1_aa(d,k)
    //               += -0.50 <j,k||c,d>_abab r1_a(a) t2_abab(c,d,j,i) t1_1_bb(b,k)
    //               += -0.50 <j,k||d,c>_abab r1_a(a) t2_abab(d,c,j,i) t1_1_bb(b,k)
    //               += +1.00 <j,k||c,d>_abab r1_a(a) t2_abab(c,b,j,i) t1_1_bb(d,k)
    //               += -0.50 <k,j||c,d>_abab r1_a(a) t2_abab(c,b,k,j) t1_1_bb(d,i)
    //               += -0.50 <j,k||c,d>_abab r1_a(a) t2_abab(c,b,j,k) t1_1_bb(d,i)
    //               += +1.00 d-_bb(j,i) r1_a(a) t0_1 t1_1_bb(b,j)
    //               += -1.00 d-_aa(j,c) r1_a(a) t0_1 t2_1_abab(c,b,j,i)
    //               += +1.00 d-_bb(j,c) r1_a(a) t0_1 t2_1_bbbb(c,b,i,j)
    //               += -1.00 d-_bb(b,c) r1_a(a) t0_1 t1_1_bb(c,i)
    //               += -0.50 <k,j||d,c>_abab r2_abb(d,b,i) t2_1_abab(a,c,k,j)
    //               += -0.50 <j,k||d,c>_abab r2_abb(d,b,i) t2_1_abab(a,c,j,k)
    //               += +1.00 <a,j||d,c>_abab r2_abb(d,b,i) t1_1_bb(c,j)
    //               += +2.00 d-_bb(j,i) r1_1_a(a) t1_1_bb(b,j)
    //               += -2.00 d-_aa(j,c) r1_1_a(a) t2_1_abab(c,b,j,i)
    //               += +2.00 d-_bb(j,c) r1_1_a(a) t2_1_bbbb(c,b,i,j)
    //               += -2.00 d-_bb(b,c) r1_1_a(a) t1_1_bb(c,i)
    //               += +0.50 <k,j||d,c>_abab r1_a(d) t2_abab(a,b,k,j) t1_1_bb(c,i)
    //               += +0.50 <j,k||d,c>_abab r1_a(d) t2_abab(a,b,j,k) t1_1_bb(c,i)
    //               += -1.00 <a,j||d,c>_abab r1_a(d) t2_1_bbbb(c,b,i,j)
    //               += +1.00 <a,b||d,c>_abab r1_a(d) t1_1_bb(c,i)
    //               += +2.00 d-_bb(j,c) r2_1_abb(a,b,j) t1_1_bb(c,i)
    //               += +1.00 r2_abb(a,b,i) t0_1 w0
    //               += +0.25 <k,j||c,d>_abab r2_abb(a,b,i) t2_1_abab(c,d,k,j)
    //               += +0.25 <j,k||c,d>_abab r2_abb(a,b,i) t2_1_abab(c,d,j,k)
    //               += +0.25 <k,j||d,c>_abab r2_abb(a,b,i) t2_1_abab(d,c,k,j)
    //               += +0.25 <j,k||d,c>_abab r2_abb(a,b,i) t2_1_abab(d,c,j,k)
    //               += -1.00 d-_bb(j,c) r2_abb(a,b,i) t0_1 t1_1_bb(c,j)
    //               += +0.25 <k,j||c,d>_bbbb r2_abb(a,b,i) t2_1_bbbb(c,d,k,j)
    //               += -1.00 d-_aa(j,c) r2_abb(a,b,i) t0_1 t1_1_aa(c,j)
    //               += +1.00 f_bb(j,c) r2_abb(a,b,i) t1_1_bb(c,j)
    //               += +0.25 <k,j||c,d>_aaaa r2_abb(a,b,i) t2_1_aaaa(c,d,k,j)
    //               += +1.00 f_aa(j,c) r2_abb(a,b,i) t1_1_aa(c,j)
    //               += -1.00 d-_bb(b,c) r2_1_abb(a,c,i) t0_1
    //               += -1.00 d+_bb(b,i) r1_a(a)
    //               += -1.00 <j,b||d,c>_abab r1_a(d) t2_1_abab(a,c,j,i)
    //               += -0.50 <j,k||c,d>_abab r2_1_abb(a,b,k) t2_abab(c,d,j,i)
    //               += -0.50 <j,k||d,c>_abab r2_1_abb(a,b,k) t2_abab(d,c,j,i)
    //               += -1.00 <j,k||c,i>_abab r2_abb(a,b,k) t1_1_aa(c,j)
    //               += -0.50 <j,k||c,d>_abab r2_abb(a,b,k) t2_1_abab(c,d,j,i)
    //               += -0.50 <j,k||d,c>_abab r2_abb(a,b,k) t2_1_abab(d,c,j,i)
    //               += +1.00 f_bb(b,c) r2_1_abb(a,c,i)
    //               += +1.00 d+_bb(j,c) r1_a(a) t2_bbbb(c,b,i,j)
    //               += -1.00 d+_aa(j,c) r1_a(a) t2_abab(c,b,j,i)
    //               += +1.00 <a,b||c,i>_abab r1_1_a(c)
    //               += -1.00 d-_aa(a,c) r2_1_abb(c,b,i) t0_1
    //               += +1.00 <j,a||c,d>_aaaa r1_a(d) t2_1_abab(c,b,j,i)
    //               += +1.00 <j,a||c,d>_aaaa r1_1_a(d) t2_abab(c,b,j,i)
    //               += -1.00 d+_aa(j,j) r2_abb(a,b,i)
    //               += +1.00 <k,j||c,i>_bbbb r2_abb(a,b,k) t1_1_bb(c,j)
    //               += -0.50 <k,j||c,d>_bbbb r2_abb(a,b,k) t2_1_bbbb(c,d,i,j)
    //               += +0.25 <k,j||c,d>_abab r2_abb(c,d,i) t2_1_abab(a,b,k,j)
    //               += +0.25 <j,k||c,d>_abab r2_abb(c,d,i) t2_1_abab(a,b,j,k)
    //               += +0.25 <k,j||d,c>_abab r2_abb(d,c,i) t2_1_abab(a,b,k,j)
    //               += +0.25 <j,k||d,c>_abab r2_abb(d,c,i) t2_1_abab(a,b,j,k)
    //               += -1.00 <k,j||c,d>_aaaa r1_a(d) t2_abab(a,b,j,i) t1_1_aa(c,k)
    //               += -1.00 f_aa(j,c) r1_1_a(c) t2_abab(a,b,j,i)
    //               += +0.50 <k,j||c,d>_aaaa r2_1_aaa(c,d,k) t2_abab(a,b,j,i)
    //               += -1.00 <j,k||d,c>_abab r1_a(d) t2_abab(a,b,j,i) t1_1_bb(c,k)
    //               += -0.50 <j,k||c,d>_abab r2_1_abb(c,d,k) t2_abab(a,b,j,i)
    //               += -0.50 <j,k||d,c>_abab r2_1_abb(d,c,k) t2_abab(a,b,j,i)
    //               += +2.00 d-_bb(j,c) r2_1_abb(a,c,i) t1_1_bb(b,j)
    //               += -1.00 f_bb(j,c) r2_abb(a,c,i) t1_1_bb(b,j)
    //               += -0.50 <a,j||c,d>_abab r2_abb(c,d,i) t1_1_bb(b,j)
    //               += -0.50 <a,j||d,c>_abab r2_abb(d,c,i) t1_1_bb(b,j)
    //               += -1.00 <a,j||c,i>_abab r1_a(c) t1_1_bb(b,j)
    //               += -1.00 <k,j||c,i>_bbbb r2_abb(a,c,k) t1_1_bb(b,j)
    //               += +1.00 <k,j||c,i>_abab r2_aaa(c,a,k) t1_1_bb(b,j)
    //               += +1.00 d-_bb(j,c) r2_abb(a,c,i) t0_1 t1_1_bb(b,j)
    //               += +1.00 d-_bb(j,c) r2_abb(a,b,j) t0_1 t1_1_bb(c,i)
    //               += +1.00 d-_aa(j,c) r1_1_a(c) t2_abab(a,b,j,i) t0_1
    //               += +2.00 d-_aa(j,c) r2_1_abb(c,b,i) t1_1_aa(a,j)
    //               += -0.50 <j,b||c,d>_abab r2_abb(c,d,i) t1_1_aa(a,j)
    //               += -0.50 <j,b||d,c>_abab r2_abb(d,c,i) t1_1_aa(a,j)
    //               += -1.00 <j,b||c,i>_abab r1_a(c) t1_1_aa(a,j)
    //               += -1.00 f_aa(j,c) r2_abb(c,b,i) t1_1_aa(a,j)
    //               += +1.00 <j,k||c,i>_abab r2_abb(c,b,k) t1_1_aa(a,j)
    //               += +1.00 d-_aa(j,c) r2_abb(c,b,i) t0_1 t1_1_aa(a,j)
    sigmar2_1_abb("R,a,b,i") += tmps_["112_bab_Lvvo"]("R,b,a,i");
    tmps_["112_bab_Lvvo"].~TArrayD();

    // flops: o2v1L1  = o2v2L1 o2v3L1 o2v1L1
    //  mems: o2v1L1  = o2v1L1 o2v1L1 o2v1L1
    tmps_["116_aaa_Lvoo"]("R,a,i,j")  = 2.00 * eri["aaaa_vovo"]("a,i,b,j") * r1["a"]("R,b");
    tmps_["116_aaa_Lvoo"]("R,a,i,j") += eri["aaaa_vovv"]("a,i,b,c") * r2["aaa"]("R,b,c,j");

    // flops: o2v1L1  = o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v1L1 o2v1L1 o2v1L1 o2v2L1 o2v1L1
    //  mems: o2v1L1  = o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1 o2v1L1
    tmps_["115_aaa_Lvoo"]("R,a,i,j")  = -1.00 * r2["aaa"]("R,b,a,k") * eri["aaaa_oovo"]("i,k,b,j");
    tmps_["115_aaa_Lvoo"]("R,a,i,j") += r2["aaa"]("R,b,a,j") * f["aa_ov"]("i,b");
    tmps_["115_aaa_Lvoo"]("R,a,i,j") += r2["abb"]("R,a,c,l") * eri["abba_oovo"]("i,l,c,j");
    tmps_["115_aaa_Lvoo"]("R,a,i,j") -= 2.00 * r2_1["aaa"]("R,b,a,j") * dp["aa_ov"]("i,b");
    tmps_["115_aaa_Lvoo"]("R,a,i,j") -= r2["aaa"]("R,b,a,j") * reused_["4_aa_ov"]("i,b");

    // flops: o2v1L1  = o2v2 o2v2L1
    //  mems: o2v1L1  = o2v2 o2v1L1
    tmps_["114_aaa_Lvoo"]("R,a,i,j")  = (reused_["44_aaaa_voov"]("a,i,j,b") + -1.00 * reused_["49_aaaa_voov"]("a,i,j,b")) * r1["a"]("R,b");

    // flops: o1v2L1  = o2v3L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["113_abb_Lvov"]("R,a,i,b")  = r2_1["aaa"]("R,c,a,j") * eri["abab_oovv"]("j,i,c,b");

    // flops: o1v2L1  = o2v3L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o2v3L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["117_aaa_Lvvo"]("R,a,b,i")  = -1.00 * tmps_["72_abb_Lvov"]("R,a,j,c") * t2_1["abab"]("b,c,i,j");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") += t1_1["aa"]("a,k") * tmps_["114_aaa_Lvoo"]("R,b,i,k");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") += tmps_["115_aaa_Lvoo"]("R,a,l,i") * t1_1["aa"]("b,l");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") -= eri["abba_vovo"]("b,j,c,i") * r2_1["abb"]("R,a,c,j");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") -= r1["a"]("R,d") * reused_["42_aaaa_vvvo"]("b,d,a,i");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") -= r1_1["a"]("R,a") * reused_["36_aa_vo"]("b,i");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") -= r2["aaa"]("R,d,a,l") * reused_["99_aaaa_vovo"]("b,l,d,i");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") += f["aa_vo"]("b,i") * r1_1["a"]("R,a");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") -= 2.00 * r1_1["a"]("R,a") * reused_["34_aa_vo"]("b,i");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") -= r2["abb"]("R,a,f,m") * reused_["79_aabb_voov"]("b,i,m,f");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") += r2["aaa"]("R,d,a,i") * reused_["14_aa_vv"]("b,d");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") -= r1_1["a"]("R,d") * reused_["64_aaaa_vvvo"]("b,d,a,i");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") += r2_1["aaa"]("R,e,a,i") * reused_["3_aa_vv"]("b,e");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") += r2_1["aaa"]("R,d,a,i") * reused_["2_aa_vv"]("b,d");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") += r2["abb"]("R,a,f,m") * reused_["98_aabb_voov"]("b,i,m,f");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") += r2["abb"]("R,a,f,j") * reused_["84_abba_vovo"]("b,j,f,i");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") -= 0.50 * r2_1["aaa"]("R,d,a,i") * reused_["1_aa_vv"]("b,d");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") += dp["aa_vv"]("b,e") * r2["aaa"]("R,e,a,i");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") -= r1_1["a"]("R,d") * reused_["26_aaaa_vvvo"]("b,d,a,i");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") -= f["aa_vv"]("b,e") * r2_1["aaa"]("R,e,a,i");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") += r2_1["abb"]("R,a,f,m") * reused_["52_aabb_voov"]("b,i,m,f");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") -= r2_1["abb"]("R,a,f,m") * reused_["47_aabb_voov"]("b,i,m,f");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") += eri["aaaa_vovo"]("b,l,e,i") * r2_1["aaa"]("R,e,a,l");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") -= r1["a"]("R,d") * reused_["96_aaaa_vvvo"]("b,d,a,i");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") -= r1["a"]("R,a") * reused_["62_aa_vo"]("b,i");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") -= r1["a"]("R,a") * reused_["32_aa_vo"]("b,i");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") += r2["aaa"]("R,d,a,k") * reused_["97_aaaa_voov"]("b,i,k,d");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") += r2_1["aaa"]("R,d,a,k") * reused_["44_aaaa_voov"]("b,i,k,d");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") -= dp["aa_vo"]("b,i") * r1["a"]("R,a");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") += r2["aaa"]("R,d,a,i") * reused_["13_aa_vv"]("b,d");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") -= t1_1["aa"]("a,i") * tmps_["111_a_Lv"]("R,b");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") += 0.50 * t1_1["aa"]("a,l") * tmps_["116_aaa_Lvoo"]("R,b,l,i");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") -= t2["abab"]("b,c,i,j") * tmps_["113_abb_Lvov"]("R,a,j,c");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") += t1_1["aa"]("b,i") * tmps_["110_a_Lv"]("R,a");
    tmps_["117_aaa_Lvvo"]("R,a,b,i") += r1["a"]("R,a") * reused_["100_aa_vo"]("b,i");
    tmps_["116_aaa_Lvoo"].~TArrayD();
    tmps_["115_aaa_Lvoo"].~TArrayD();
    tmps_["114_aaa_Lvoo"].~TArrayD();
    tmps_["113_abb_Lvov"].~TArrayD();
    tmps_["111_a_Lv"].~TArrayD();
    tmps_["110_a_Lv"].~TArrayD();
    tmps_["72_abb_Lvov"].~TArrayD();

    // sigmar2_1_aaa += -1.00 P(a,b) <k,j||d,c>_abab r1_a(b) t2_abab(a,c,i,j) t1_1_aa(d,k)
    //               += +1.00 P(a,b) <k,j||c,d>_aaaa r1_a(d) t2_aaaa(c,a,i,j) t1_1_aa(b,k)
    //               += +1.00 P(a,b) <k,j||d,c>_abab r1_a(d) t2_abab(a,c,i,j) t1_1_aa(b,k)
    //               += +1.00 P(a,b) <k,j||d,c>_abab r2_aaa(d,b,k) t2_1_abab(a,c,i,j)
    //               += +1.00 P(a,b) <j,k||i,c>_abab r2_abb(b,c,k) t1_1_aa(a,j)
    //               += +2.00 P(a,b) d-_aa(j,c) r2_1_aaa(c,b,i) t1_1_aa(a,j)
    //               += -1.00 P(a,b) f_aa(j,c) r2_aaa(c,b,i) t1_1_aa(a,j)
    //               += -1.00 P(a,b) <k,j||c,i>_aaaa r2_aaa(c,b,k) t1_1_aa(a,j)
    //               += +1.00 P(a,b) d-_aa(j,c) r2_aaa(c,b,i) t0_1 t1_1_aa(a,j)
    //               += -1.00 P(a,b) <a,j||i,c>_abab r2_1_abb(b,c,j)
    //               += +1.00 P(a,b) <a,j||d,c>_abab r1_a(d) t2_1_abab(b,c,i,j)
    //               += +1.00 P(a,b) d-_aa(a,i) r1_1_a(b) t0_1
    //               += +0.50 P(a,b) <k,j||c,i>_aaaa r1_1_a(b) t2_aaaa(c,a,k,j)
    //               += +0.50 P(a,b) <k,j||i,c>_abab r1_1_a(b) t2_abab(a,c,k,j)
    //               += +0.50 P(a,b) <j,k||i,c>_abab r1_1_a(b) t2_abab(a,c,j,k)
    //               += +1.00 P(a,b) f_aa(j,c) r1_1_a(b) t2_aaaa(c,a,i,j)
    //               += -1.00 P(a,b) f_bb(j,c) r1_1_a(b) t2_abab(a,c,i,j)
    //               += -0.50 P(a,b) <a,j||c,d>_abab r1_1_a(b) t2_abab(c,d,i,j)
    //               += -0.50 P(a,b) <a,j||d,c>_abab r1_1_a(b) t2_abab(d,c,i,j)
    //               += +1.00 P(a,b) d-_bb(j,j) r1_1_a(b) t1_1_aa(a,i)
    //               += +1.00 P(a,b) d-_aa(j,j) r1_1_a(b) t1_1_aa(a,i)
    //               += +0.50 P(a,b) <j,a||c,d>_aaaa r1_1_a(b) t2_aaaa(c,d,i,j)
    //               += +1.00 P(a,b) d-_bb(j,c) r1_1_a(b) t2_abab(a,c,i,j) t0_1
    //               += -1.00 P(a,b) d-_aa(j,c) r1_1_a(b) t2_aaaa(c,a,i,j) t0_1
    //               += -1.00 P(a,b) <j,a||c,d>_aaaa r2_aaa(d,b,j) t1_1_aa(c,i)
    //               += -1.00 P(a,b) f_aa(a,i) r1_1_a(b)
    //               += +2.00 P(a,b) d-_aa(a,c) r1_1_a(b) t1_1_aa(c,i)
    //               += +2.00 P(a,b) d-_bb(j,c) r1_1_a(b) t2_1_abab(a,c,i,j)
    //               += -2.00 P(a,b) d-_aa(j,c) r1_1_a(b) t2_1_aaaa(c,a,i,j)
    //               += -2.00 P(a,b) d-_aa(j,i) r1_1_a(b) t1_1_aa(a,j)
    //               += +1.00 P(a,b) <j,k||c,d>_abab r2_abb(b,d,k) t2_1_aaaa(c,a,i,j)
    //               += +1.00 P(a,b) <j,a||c,d>_aaaa r2_aaa(d,b,i) t1_1_aa(c,j)
    //               += -0.50 P(a,b) <k,j||c,d>_aaaa r2_aaa(d,b,i) t2_1_aaaa(c,a,k,j)
    //               += -1.00 P(a,b) <j,a||c,d>_aaaa r1_1_a(d) t2_aaaa(c,b,i,j)
    //               += -1.00 P(a,b) d-_aa(a,c) r2_1_aaa(c,b,i) t0_1
    //               += -0.50 P(a,b) <k,j||d,c>_abab r2_1_aaa(d,b,i) t2_abab(a,c,k,j)
    //               += -0.50 P(a,b) <j,k||d,c>_abab r2_1_aaa(d,b,i) t2_abab(a,c,j,k)
    //               += +1.00 P(a,b) <k,j||c,d>_bbbb r2_abb(b,d,k) t2_1_abab(a,c,i,j)
    //               += -1.00 P(a,b) <a,j||c,d>_abab r2_abb(b,d,j) t1_1_aa(c,i)
    //               += -0.50 P(a,b) <k,j||c,d>_aaaa r2_1_aaa(d,b,i) t2_aaaa(c,a,k,j)
    //               += -1.00 P(a,b) d+_aa(a,c) r2_aaa(c,b,i)
    //               += +1.00 P(a,b) <a,j||d,c>_abab r1_1_a(d) t2_abab(b,c,i,j)
    //               += +1.00 P(a,b) f_aa(a,c) r2_1_aaa(c,b,i)
    //               += +1.00 P(a,b) <k,j||c,d>_bbbb r2_1_abb(b,d,k) t2_abab(a,c,i,j)
    //               += +1.00 P(a,b) <j,k||c,d>_abab r2_1_abb(b,d,k) t2_aaaa(c,a,i,j)
    //               += +1.00 P(a,b) <j,a||c,i>_aaaa r2_1_aaa(c,b,j)
    //               += -1.00 P(a,b) <j,a||c,d>_aaaa r1_a(d) t2_1_aaaa(c,b,i,j)
    //               += +0.50 P(a,b) <k,j||c,d>_abab r1_a(b) t2_abab(c,d,i,j) t1_1_aa(a,k)
    //               += +0.50 P(a,b) <k,j||d,c>_abab r1_a(b) t2_abab(d,c,i,j) t1_1_aa(a,k)
    //               += +1.00 P(a,b) <k,j||c,d>_bbbb r1_a(b) t2_abab(a,c,i,j) t1_1_bb(d,k)
    //               += +0.50 P(a,b) <k,j||c,d>_aaaa r1_a(b) t2_aaaa(c,d,i,j) t1_1_aa(a,k)
    //               += +0.50 P(a,b) <k,j||c,d>_aaaa r1_a(b) t2_aaaa(c,a,k,j) t1_1_aa(d,i)
    //               += +1.00 P(a,b) <j,k||c,d>_abab r1_a(b) t2_aaaa(c,a,i,j) t1_1_bb(d,k)
    //               += +0.50 P(a,b) <k,j||d,c>_abab r1_a(b) t2_abab(a,c,k,j) t1_1_aa(d,i)
    //               += +0.50 P(a,b) <j,k||d,c>_abab r1_a(b) t2_abab(a,c,j,k) t1_1_aa(d,i)
    //               += -1.00 P(a,b) <k,j||c,d>_aaaa r1_a(b) t2_aaaa(c,a,i,j) t1_1_aa(d,k)
    //               += -2.00 P(a,b) d-_aa(j,c) r1_a(b) t1_1_aa(a,j) t1_1_aa(c,i)
    //               += +1.00 P(a,b) d-_aa(a,c) r1_a(b) t0_1 t1_1_aa(c,i)
    //               += +1.00 P(a,b) d-_bb(j,c) r1_a(b) t0_1 t2_1_abab(a,c,i,j)
    //               += -1.00 P(a,b) d-_aa(j,c) r1_a(b) t0_1 t2_1_aaaa(c,a,i,j)
    //               += -1.00 P(a,b) d-_aa(j,i) r1_a(b) t0_1 t1_1_aa(a,j)
    //               += +1.00 P(a,b) f_aa(j,i) r1_a(b) t1_1_aa(a,j)
    //               += -1.00 P(a,b) f_aa(a,c) r1_a(b) t1_1_aa(c,i)
    //               += +1.00 P(a,b) f_aa(j,c) r1_a(b) t2_1_aaaa(c,a,i,j)
    //               += -1.00 P(a,b) <a,j||i,c>_abab r1_a(b) t1_1_bb(c,j)
    //               += -0.50 P(a,b) <a,j||c,d>_abab r1_a(b) t2_1_abab(c,d,i,j)
    //               += -0.50 P(a,b) <a,j||d,c>_abab r1_a(b) t2_1_abab(d,c,i,j)
    //               += +0.50 P(a,b) <k,j||c,i>_aaaa r1_a(b) t2_1_aaaa(c,a,k,j)
    //               += -1.00 P(a,b) <j,a||c,i>_aaaa r1_a(b) t1_1_aa(c,j)
    //               += -1.00 P(a,b) f_bb(j,c) r1_a(b) t2_1_abab(a,c,i,j)
    //               += +0.50 P(a,b) <k,j||i,c>_abab r1_a(b) t2_1_abab(a,c,k,j)
    //               += +0.50 P(a,b) <j,k||i,c>_abab r1_a(b) t2_1_abab(a,c,j,k)
    //               += +0.50 P(a,b) <j,a||c,d>_aaaa r1_a(b) t2_1_aaaa(c,d,i,j)
    //               += -1.00 P(a,b) r1_a(b) t1_1_aa(a,i) w0
    //               += +1.00 P(a,b) d-_aa(j,c) r1_a(b) t1_1_aa(a,i) t1_1_aa(c,j)
    //               += +1.00 P(a,b) d-_bb(j,c) r1_a(b) t1_1_aa(a,i) t1_1_bb(c,j)
    //               += +1.00 P(a,b) d+_bb(j,c) r1_a(b) t2_abab(a,c,i,j)
    //               += -1.00 P(a,b) d+_aa(j,c) r1_a(b) t2_aaaa(c,a,i,j)
    //               += +1.00 P(a,b) <k,j||c,d>_aaaa r2_aaa(d,b,k) t2_1_aaaa(c,a,i,j)
    //               += +1.00 P(a,b) <k,j||c,d>_aaaa r2_1_aaa(d,b,k) t2_aaaa(c,a,i,j)
    //               += +1.00 P(a,b) d+_aa(a,i) r1_a(b)
    //               += -0.50 P(a,b) <k,j||d,c>_abab r2_aaa(d,b,i) t2_1_abab(a,c,k,j)
    //               += -0.50 P(a,b) <j,k||d,c>_abab r2_aaa(d,b,i) t2_1_abab(a,c,j,k)
    //               += +1.00 P(a,b) <a,j||d,c>_abab r2_aaa(d,b,i) t1_1_bb(c,j)
    //               += +1.00 P(a,b) d-_aa(j,c) r1_a(c) t1_1_aa(a,j) t1_1_aa(b,i)
    //               += -1.00 P(a,b) d-_aa(a,c) r1_1_a(c) t1_1_aa(b,i)
    //               += +0.50 P(a,b) <j,a||c,d>_aaaa r2_aaa(c,d,i) t1_1_aa(b,j)
    //               += +1.00 P(a,b) <j,a||c,i>_aaaa r1_a(c) t1_1_aa(b,j)
    //               += +1.00 P(a,b) <k,j||d,c>_abab r2_1_aaa(d,b,k) t2_abab(a,c,i,j)
    //               += -1.00 P(a,b) d-_aa(j,c) r2_1_aaa(c,b,j) t1_1_aa(a,i)
    //               += +1.00 P(a,b) d-_bb(j,c) r2_1_abb(b,c,j) t1_1_aa(a,i)
    sigmar2_1_aaa("R,a,b,i") -= tmps_["117_aaa_Lvvo"]("R,b,a,i");
    sigmar2_1_aaa("R,a,b,i") += tmps_["117_aaa_Lvvo"]("R,a,b,i");
    tmps_["117_aaa_Lvvo"].~TArrayD();

    // flops: o1v2L1  = o0v1L1 o1v2L1 o2v3L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1
    //  mems: o1v2L1  = o0v1L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["118_aaa_Lvov"]("L,a,i,b")  = (tmps_["97_a_Lv"]("L,a") + -1.00 * tmps_["53_a_Lv"]("L,a")) * dp["aa_ov"]("i,b");
    tmps_["118_aaa_Lvov"]("L,a,i,b") -= eri["abab_oovv"]("i,j,b,c") * tmps_["66_abb_Lvvo"]("L,a,c,j");
    tmps_["118_aaa_Lvov"]("L,a,i,b") -= 2.00 * dp["aa_ov"]("k,b") * tmps_["61_aaa_Lovo"]("L,i,a,k");
    tmps_["118_aaa_Lvov"]("L,a,i,b") -= l1_1["a"]("L,a") * reused_["4_aa_ov"]("i,b");
    tmps_["118_aaa_Lvov"]("L,a,i,b") -= l1["a"]("L,a") * dp["aa_ov"]("i,b");
    tmps_["118_aaa_Lvov"]("L,a,i,b") += l1_1["a"]("L,a") * f["aa_ov"]("i,b");
    tmps_["118_aaa_Lvov"]("L,a,i,b") += l2_1["aaa"]("L,k,d,a") * reused_["70_aaaa_voov"]("d,k,i,b");
    tmps_["118_aaa_Lvov"]("L,a,i,b") += 0.50 * l2_1["aaa"]("L,i,d,a") * reused_["31_aa_vv"]("b,d");
    tmps_["118_aaa_Lvov"]("L,a,i,b") += dp["aa_vv"]("d,b") * l2["aaa"]("L,i,d,a");
    tmps_["118_aaa_Lvov"]("L,a,i,b") -= f["aa_vv"]("d,b") * l2_1["aaa"]("L,i,d,a");
    tmps_["118_aaa_Lvov"]("L,a,i,b") += l2_1["aaa"]("L,i,d,a") * reused_["3_aa_vv"]("d,b");
    tmps_["118_aaa_Lvov"]("L,a,i,b") -= l2_1["bab"]("L,l,a,e") * reused_["50_bbaa_voov"]("e,l,i,b");
    tmps_["118_aaa_Lvov"]("L,a,i,b") += l2_1["aaa"]("L,i,d,a") * reused_["2_aa_vv"]("d,b");
    tmps_["118_aaa_Lvov"]("L,a,i,b") += l2_1["bab"]("L,l,a,e") * reused_["71_bbaa_voov"]("e,l,i,b");
    tmps_["118_aaa_Lvov"]("L,a,i,b") -= l2_1["bab"]("L,l,a,e") * eri["baab_vovo"]("e,i,b,l");
    tmps_["118_aaa_Lvov"]("L,a,i,b") += l2_1["aaa"]("L,k,d,a") * eri["aaaa_vovo"]("d,i,b,k");
    tmps_["97_a_Lv"].~TArrayD();
    tmps_["66_abb_Lvvo"].~TArrayD();
    tmps_["61_aaa_Lovo"].~TArrayD();
    tmps_["53_a_Lv"].~TArrayD();

    // sigmal2_1_aaa += -1.00 P(a,b) d-_aa(i,a) t1_1_aa(c,j) l2_1_aaa(j,c,b)
    //               += +1.00 P(a,b) d-_aa(i,a) t1_1_bb(c,j) l2_1_bab(j,b,c)
    //               += +1.00 P(a,b) <i,k||a,d>_abab t2_abab(c,d,j,k) l2_1_aaa(j,c,b)
    //               += +2.00 P(a,b) d-_aa(j,a) t1_1_aa(c,j) l2_1_aaa(i,c,b)
    //               += +1.00 P(a,b) d-_aa(i,a) t0_1 l1_1_a(b)
    //               += +1.00 P(a,b) d-_aa(i,a) l1_a(b)
    //               += -1.00 P(a,b) f_aa(i,a) l1_1_a(b)
    //               += +1.00 P(a,b) <i,k||d,a>_aaaa t2_aaaa(d,c,j,k) l2_1_aaa(j,c,b)
    //               += -0.50 P(a,b) <k,j||d,a>_aaaa t2_aaaa(d,c,k,j) l2_1_aaa(i,c,b)
    //               += -1.00 P(a,b) d-_aa(c,a) l2_aaa(i,c,b)
    //               += +1.00 P(a,b) f_aa(c,a) l2_1_aaa(i,c,b)
    //               += -1.00 P(a,b) d-_aa(c,a) t0_1 l2_1_aaa(i,c,b)
    //               += +1.00 P(a,b) <i,k||a,d>_abab t2_bbbb(d,c,j,k) l2_1_bab(j,b,c)
    //               += -0.50 P(a,b) <k,j||a,d>_abab t2_abab(c,d,k,j) l2_1_aaa(i,c,b)
    //               += -0.50 P(a,b) <j,k||a,d>_abab t2_abab(c,d,j,k) l2_1_aaa(i,c,b)
    //               += +1.00 P(a,b) <i,k||d,a>_aaaa t2_abab(d,c,k,j) l2_1_bab(j,b,c)
    //               += -1.00 P(a,b) <i,c||a,j>_abab l2_1_bab(j,b,c)
    //               += +1.00 P(a,b) <i,c||a,j>_aaaa l2_1_aaa(j,c,b)
    sigmal2_1_aaa("L,a,b,i") -= tmps_["118_aaa_Lvov"]("L,b,i,a");
    sigmal2_1_aaa("L,a,b,i") += tmps_["118_aaa_Lvov"]("L,a,i,b");
    tmps_["118_aaa_Lvov"].~TArrayD();
}