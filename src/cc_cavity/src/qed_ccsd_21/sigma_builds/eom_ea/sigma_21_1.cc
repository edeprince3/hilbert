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

void hilbert::EOM_EA_QED_CCSD::sigma_ea_21_1() {

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


    // flops: o1v2L1  = o1v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["1_aaa_Lvov"]("R,a,i,b")  = t1_1["aa"]("a,i") * r1["a"]("R,b");

    // csigmar2_aaa  = -1.00 P(a,b) r1_a(b) t1_1_aa(a,i)
    csigmar2_aaa("R,a,b,i")  = -1.00 * tmps_["1_aaa_Lvov"]("R,a,i,b");
    csigmar2_aaa("R,a,b,i") += tmps_["1_aaa_Lvov"]("R,b,i,a");
    tmps_["1_aaa_Lvov"].~TArrayD();

    // flops: o1v2L1  = o1v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["2_aaa_Lvov"]("R,a,i,b")  = t1_1["aa"]("a,i") * r1_1["a"]("R,b");

    // csigmar2_1_aaa  = -1.00 P(a,b) r1_1_a(b) t1_1_aa(a,i)
    csigmar2_1_aaa("R,a,b,i")  = -1.00 * tmps_["2_aaa_Lvov"]("R,a,i,b");
    csigmar2_1_aaa("R,a,b,i") += tmps_["2_aaa_Lvov"]("R,b,i,a");
    tmps_["2_aaa_Lvov"].~TArrayD();

    // flops: o0v1L1  = o1v3L1 o0v2L1 o1v2L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v2L1 o0v1L1 o0v1L1 o0v1L1 o1v2L1 o0v1L1 o0v2L1 o0v1L1 o1v2L1 o0v1L1 o0v2L1 o0v1L1 o0v1L1 o0v1L1 o1v2L1 o0v1L1 o1v3L1 o0v1L1
    //  mems: o0v1L1  = o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1
    tmps_["3_a_Lv"]("R,a")  = -0.50 * eri["aaaa_vovv"]("a,i,b,d") * r2["aaa"]("R,b,d,i");
    tmps_["3_a_Lv"]("R,a") -= f["aa_vv"]("a,b") * r1["a"]("R,b");
    tmps_["3_a_Lv"]("R,a") += r2["abb"]("R,a,e,j") * reused_["5_bb_ov"]("j,e");
    tmps_["3_a_Lv"]("R,a") += scalars_["4"] * r1["a"]("R,a");
    tmps_["3_a_Lv"]("R,a") += scalars_["5"] * r1["a"]("R,a");
    tmps_["3_a_Lv"]("R,a") += scalars_["2"] * r1_1["a"]("R,a");
    tmps_["3_a_Lv"]("R,a") += r1["a"]("R,b") * reused_["3_aa_vv"]("a,b");
    tmps_["3_a_Lv"]("R,a") += scalars_["3"] * r1["a"]("R,a");
    tmps_["3_a_Lv"]("R,a") -= r2["aaa"]("R,b,a,i") * reused_["4_aa_ov"]("i,b");
    tmps_["3_a_Lv"]("R,a") += r1["a"]("R,d") * reused_["2_aa_vv"]("a,d");
    tmps_["3_a_Lv"]("R,a") -= f["bb_ov"]("j,e") * r2["abb"]("R,a,e,j");
    tmps_["3_a_Lv"]("R,a") -= 0.50 * r1["a"]("R,d") * reused_["1_aa_vv"]("a,d");
    tmps_["3_a_Lv"]("R,a") += scalars_["1"] * r1_1["a"]("R,a");
    tmps_["3_a_Lv"]("R,a") += f["aa_ov"]("i,b") * r2["aaa"]("R,b,a,i");
    tmps_["3_a_Lv"]("R,a") -= eri["abab_vovv"]("a,j,b,c") * r2["abb"]("R,b,c,j");

    // sigmar1_a  = -1.00 d-_aa(i,b) r1_a(a) t1_1_aa(b,i)
    //           += -1.00 d-_bb(i,b) r2_abb(a,b,i) t0_1
    //           += -0.50 <j,i||j,i>_abab r1_a(a)
    //           += -0.50 <i,j||i,j>_abab r1_a(a)
    //           += +0.25 <j,i||b,c>_bbbb r1_a(a) t2_bbbb(b,c,j,i)
    //           += +1.00 f_bb(i,i) r1_a(a)
    //           += +0.25 <j,i||b,c>_aaaa r1_a(a) t2_aaaa(b,c,j,i)
    //           += -1.00 d-_aa(i,i) r1_a(a) t0_1
    //           += +1.00 f_aa(i,i) r1_a(a)
    //           += -0.50 <j,i||j,i>_bbbb r1_a(a)
    //           += -0.50 <j,i||j,i>_aaaa r1_a(a)
    //           += +0.25 <j,i||b,c>_abab r1_a(a) t2_abab(b,c,j,i)
    //           += +0.25 <i,j||b,c>_abab r1_a(a) t2_abab(b,c,i,j)
    //           += +0.25 <j,i||c,b>_abab r1_a(a) t2_abab(c,b,j,i)
    //           += +0.25 <i,j||c,b>_abab r1_a(a) t2_abab(c,b,i,j)
    //           += -1.00 d-_bb(i,i) r1_a(a) t0_1
    //           += +1.00 f_aa(a,b) r1_a(b)
    //           += -1.00 d-_aa(i,i) r1_1_a(a)
    //           += -0.50 <i,a||b,c>_aaaa r2_aaa(b,c,i)
    //           += -1.00 d-_aa(a,b) r1_a(b) t0_1
    //           += -1.00 d-_bb(i,b) r1_a(a) t1_1_bb(b,i)
    //           += +1.00 d-_aa(i,b) r2_aaa(b,a,i) t0_1
    //           += -0.50 <j,i||c,b>_abab r1_a(c) t2_abab(a,b,j,i)
    //           += -0.50 <i,j||c,b>_abab r1_a(c) t2_abab(a,b,i,j)
    //           += +1.00 f_bb(i,b) r2_abb(a,b,i)
    //           += -0.50 <j,i||b,c>_aaaa r1_a(c) t2_aaaa(b,a,j,i)
    //           += -1.00 d-_bb(i,i) r1_1_a(a)
    //           += -1.00 f_aa(i,b) r2_aaa(b,a,i)
    //           += +0.50 <a,i||b,c>_abab r2_abb(b,c,i)
    //           += +0.50 <a,i||c,b>_abab r2_abb(c,b,i)
    sigmar1_a("R,a")  = -1.00 * tmps_["3_a_Lv"]("R,a");
    tmps_["3_a_Lv"].~TArrayD();

    // flops: o0v1L1  = o1v3L1
    //  mems: o0v1L1  = o0v1L1
    tmps_["4_a_Lv"]("R,a")  = -1.00 * eri["abab_vovv"]("a,i,b,c") * r2["abb"]("R,b,c,i");
    tmps_["4_a_Lv"].~TArrayD();

    // flops: o0v1L1  = o0v1L1 o0v1L1
    //  mems: o0v1L1  = o0v1L1 o0v1L1
    tmps_["5_a_Lv"]("R,a")  = r1_1["a"]("R,a");
    tmps_["5_a_Lv"]("R,a") += t0_1 * r1["a"]("R,a");

    // csigmar1_a  = +1.00 r1_a(a) t0_1
    //            += +1.00 r1_1_a(a)
    csigmar1_a("R,a")  = tmps_["5_a_Lv"]("R,a");
    tmps_["5_a_Lv"].~TArrayD();

    // flops: o0v1L1  = o0v1L1 o0v1L1
    //  mems: o0v1L1  = o0v1L1 o0v1L1
    tmps_["6_a_Lv"]("R,a")  = r1["a"]("R,a");
    tmps_["6_a_Lv"]("R,a") += t0_1 * r1_1["a"]("R,a");

    // csigmar1_1_a  = +1.00 r1_1_a(a) t0_1
    //              += +1.00 r1_a(a)
    csigmar1_1_a("R,a")  = tmps_["6_a_Lv"]("R,a");
    tmps_["6_a_Lv"].~TArrayD();

    // flops: o0v1L1  = o0v1L1 o0v1L1
    //  mems: o0v1L1  = o0v1L1 o0v1L1
    tmps_["7_a_Lv"]("L,a")  = l1_1["a"]("L,a");
    tmps_["7_a_Lv"]("L,a") += t0_1 * l1["a"]("L,a");

    // csigmal1_a  = +1.00 t0_1 l1_a(a)
    //            += +1.00 l1_1_a(a)
    csigmal1_a("L,a")  = tmps_["7_a_Lv"]("L,a");
    tmps_["7_a_Lv"].~TArrayD();

    // flops: o0v1L1  = o0v1L1 o0v1L1
    //  mems: o0v1L1  = o0v1L1 o0v1L1
    tmps_["8_a_Lv"]("L,a")  = l1["a"]("L,a");
    tmps_["8_a_Lv"]("L,a") += t0_1 * l1_1["a"]("L,a");

    // csigmal1_1_a  = +1.00 l1_a(a)
    //              += +1.00 t0_1 l1_1_a(a)
    csigmal1_1_a("L,a")  = tmps_["8_a_Lv"]("L,a");
    tmps_["8_a_Lv"].~TArrayD();

    // flops: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["9_abb_Lvvo"]("R,a,b,i")  = r1["a"]("R,a") * t1_1["bb"]("b,i");
    tmps_["9_abb_Lvvo"]("R,a,b,i") += r2_1["abb"]("R,a,b,i");
    tmps_["9_abb_Lvvo"]("R,a,b,i") += t0_1 * r2["abb"]("R,a,b,i");

    // csigmar2_abb  = +1.00 r1_a(a) t1_1_bb(b,i)
    //              += +1.00 r2_1_abb(a,b,i)
    //              += +1.00 r2_abb(a,b,i) t0_1
    csigmar2_abb("R,a,b,i")  = tmps_["9_abb_Lvvo"]("R,a,b,i");
    tmps_["9_abb_Lvvo"].~TArrayD();

    // flops: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["10_abb_Lvvo"]("R,a,b,i")  = r1_1["a"]("R,a") * t1_1["bb"]("b,i");
    tmps_["10_abb_Lvvo"]("R,a,b,i") += r2["abb"]("R,a,b,i");
    tmps_["10_abb_Lvvo"]("R,a,b,i") += t0_1 * r2_1["abb"]("R,a,b,i");

    // csigmar2_1_abb  = +1.00 r1_1_a(a) t1_1_bb(b,i)
    //                += +1.00 r2_abb(a,b,i)
    //                += +1.00 r2_1_abb(a,b,i) t0_1
    csigmar2_1_abb("R,a,b,i")  = tmps_["10_abb_Lvvo"]("R,a,b,i");
    tmps_["10_abb_Lvvo"].~TArrayD();

    // flops: o1v2L1  = o1v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1
    tmps_["11_bab_Lovv"]("L,i,a,b")  = l2_1["bab"]("L,i,a,b");
    tmps_["11_bab_Lovv"]("L,i,a,b") += t0_1 * l2["bab"]("L,i,a,b");

    // csigmal2_abb  = +1.00 t0_1 l2_bab(i,a,b)
    //              += +1.00 l2_1_bab(i,a,b)
    csigmal2_abb("L,a,b,i")  = tmps_["11_bab_Lovv"]("L,i,a,b");
    tmps_["11_bab_Lovv"].~TArrayD();

    // flops: o1v2L1  = o1v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1
    tmps_["12_bab_Lovv"]("L,i,a,b")  = l2["bab"]("L,i,a,b");
    tmps_["12_bab_Lovv"]("L,i,a,b") += t0_1 * l2_1["bab"]("L,i,a,b");

    // csigmal2_1_abb  = +1.00 l2_bab(i,a,b)
    //                += +1.00 t0_1 l2_1_bab(i,a,b)
    csigmal2_1_abb("L,a,b,i")  = tmps_["12_bab_Lovv"]("L,i,a,b");
    tmps_["12_bab_Lovv"].~TArrayD();

    // flops: o1v2L1  = o1v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1
    tmps_["13_aaa_Lvvo"]("R,a,b,i")  = r2_1["aaa"]("R,a,b,i");
    tmps_["13_aaa_Lvvo"]("R,a,b,i") += t0_1 * r2["aaa"]("R,a,b,i");

    // csigmar2_aaa += +1.00 r2_aaa(a,b,i) t0_1
    //              += +1.00 r2_1_aaa(a,b,i)
    csigmar2_aaa("R,a,b,i") += tmps_["13_aaa_Lvvo"]("R,a,b,i");
    tmps_["13_aaa_Lvvo"].~TArrayD();

    // flops: o1v2L1  = o1v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1
    tmps_["14_aaa_Lvvo"]("R,a,b,i")  = r2["aaa"]("R,a,b,i");
    tmps_["14_aaa_Lvvo"]("R,a,b,i") += t0_1 * r2_1["aaa"]("R,a,b,i");

    // csigmar2_1_aaa += +1.00 r2_1_aaa(a,b,i) t0_1
    //                += +1.00 r2_aaa(a,b,i)
    csigmar2_1_aaa("R,a,b,i") += tmps_["14_aaa_Lvvo"]("R,a,b,i");
    tmps_["14_aaa_Lvvo"].~TArrayD();

    // flops: o1v2L1  = o1v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1
    tmps_["15_aaa_Lovv"]("L,i,a,b")  = l2_1["aaa"]("L,i,a,b");
    tmps_["15_aaa_Lovv"]("L,i,a,b") += t0_1 * l2["aaa"]("L,i,a,b");

    // csigmal2_aaa  = +1.00 l2_1_aaa(i,a,b)
    //              += +1.00 t0_1 l2_aaa(i,a,b)
    csigmal2_aaa("L,a,b,i")  = tmps_["15_aaa_Lovv"]("L,i,a,b");
    tmps_["15_aaa_Lovv"].~TArrayD();

    // flops: o1v2L1  = o1v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1
    tmps_["16_aaa_Lovv"]("L,i,a,b")  = l2["aaa"]("L,i,a,b");
    tmps_["16_aaa_Lovv"]("L,i,a,b") += t0_1 * l2_1["aaa"]("L,i,a,b");

    // csigmal2_1_aaa  = +1.00 l2_aaa(i,a,b)
    //                += +1.00 t0_1 l2_1_aaa(i,a,b)
    csigmal2_1_aaa("L,a,b,i")  = tmps_["16_aaa_Lovv"]("L,i,a,b");
    tmps_["16_aaa_Lovv"].~TArrayD();

    // flops: o1v0L1  = o2v2L1
    //  mems: o1v0L1  = o1v0L1
    tmps_["20_a_Lo"]("L,i")  = l2_1["bab"]("L,j,a,b") * t2["abab"]("a,b,i,j");

    // flops: o1v0L1  = o2v2L1
    //  mems: o1v0L1  = o1v0L1
    tmps_["19_a_Lo"]("L,i")  = l2_1["aaa"]("L,j,a,b") * t2["aaaa"]("a,b,j,i");

    // flops: o1v2L1  = o2v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["18_aaa_Lvvo"]("L,a,b,i")  = 0.50 * l2_1["aaa"]("L,j,a,b") * reused_["6_aa_oo"]("i,j");

    // flops: o3v0L1  = o3v2L1
    //  mems: o3v0L1  = o3v0L1
    tmps_["17_aaa_Looo"]("L,i,j,k")  = l2_1["aaa"]("L,i,a,b") * t2["aaaa"]("a,b,j,k");

    // flops: o1v2L1  = o2v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v4L1 o1v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o3v2L1 o1v2L1 o2v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["21_aaa_Lovv"]("L,i,a,b")  = -0.50 * l2["aaa"]("L,k,a,b") * dp["aa_oo"]("i,k");
    tmps_["21_aaa_Lovv"]("L,i,a,b") -= 0.50 * l2_1["aaa"]("L,k,a,b") * reused_["8_aa_oo"]("i,k");
    tmps_["21_aaa_Lovv"]("L,i,a,b") += scalars_["3"] * l2_1["aaa"]("L,i,a,b");
    tmps_["21_aaa_Lovv"]("L,i,a,b") += 0.50 * scalars_["5"] * l2_1["aaa"]("L,i,a,b");
    tmps_["21_aaa_Lovv"]("L,i,a,b") += 0.25 * eri["aaaa_vvvv"]("d,c,a,b") * l2_1["aaa"]("L,i,c,d");
    tmps_["21_aaa_Lovv"]("L,i,a,b") += 0.50 * scalars_["1"] * l2["aaa"]("L,i,a,b");
    tmps_["21_aaa_Lovv"]("L,i,a,b") += 0.25 * l2_1["aaa"]("L,k,a,b") * reused_["7_aa_oo"]("i,k");
    tmps_["21_aaa_Lovv"]("L,i,a,b") += 0.50 * scalars_["2"] * l2["aaa"]("L,i,a,b");
    tmps_["21_aaa_Lovv"]("L,i,a,b") -= 0.50 * eri["aaaa_vovv"]("d,i,a,b") * l1_1["a"]("L,d");
    tmps_["21_aaa_Lovv"]("L,i,a,b") += 0.50 * f["aa_oo"]("i,k") * l2_1["aaa"]("L,k,a,b");
    tmps_["21_aaa_Lovv"]("L,i,a,b") += scalars_["4"] * l2_1["aaa"]("L,i,a,b");
    tmps_["21_aaa_Lovv"]("L,i,a,b") -= 0.50 * w0 * l2_1["aaa"]("L,i,a,b");
    tmps_["21_aaa_Lovv"]("L,i,a,b") += tmps_["18_aaa_Lvvo"]("L,a,b,i");
    tmps_["21_aaa_Lovv"]("L,i,a,b") -= 0.50 * eri["aaaa_oovv"]("i,j,a,b") * tmps_["20_a_Lo"]("L,j");
    tmps_["21_aaa_Lovv"]("L,i,a,b") += 0.1250 * eri["aaaa_oovv"]("k,j,a,b") * tmps_["17_aaa_Looo"]("L,i,j,k");
    tmps_["21_aaa_Lovv"]("L,i,a,b") += 0.25 * eri["aaaa_oovv"]("i,j,a,b") * tmps_["19_a_Lo"]("L,j");
    tmps_["18_aaa_Lvvo"].~TArrayD();
    tmps_["17_aaa_Looo"].~TArrayD();

    // sigmal2_1_aaa  = -2.00 d-_bb(j,c) t1_1_bb(c,j) l2_1_aaa(i,a,b)
    //               += +1.00 d-_aa(i,j) t0_1 l2_1_aaa(j,a,b)
    //               += -0.50 <k,j||k,j>_abab l2_1_aaa(i,a,b)
    //               += -0.50 <j,k||j,k>_abab l2_1_aaa(i,a,b)
    //               += +0.25 <k,j||c,d>_bbbb t2_bbbb(c,d,k,j) l2_1_aaa(i,a,b)
    //               += +1.00 f_bb(j,j) l2_1_aaa(i,a,b)
    //               += +0.25 <k,j||c,d>_aaaa t2_aaaa(c,d,k,j) l2_1_aaa(i,a,b)
    //               += -1.00 d-_aa(j,j) t0_1 l2_1_aaa(i,a,b)
    //               += +1.00 f_aa(j,j) l2_1_aaa(i,a,b)
    //               += -0.50 <k,j||k,j>_bbbb l2_1_aaa(i,a,b)
    //               += -0.50 <k,j||k,j>_aaaa l2_1_aaa(i,a,b)
    //               += +0.25 <k,j||c,d>_abab t2_abab(c,d,k,j) l2_1_aaa(i,a,b)
    //               += +0.25 <j,k||c,d>_abab t2_abab(c,d,j,k) l2_1_aaa(i,a,b)
    //               += +0.25 <k,j||d,c>_abab t2_abab(d,c,k,j) l2_1_aaa(i,a,b)
    //               += +0.25 <j,k||d,c>_abab t2_abab(d,c,j,k) l2_1_aaa(i,a,b)
    //               += -1.00 d-_bb(j,j) t0_1 l2_1_aaa(i,a,b)
    //               += +0.50 <d,c||a,b>_aaaa l2_1_aaa(i,d,c)
    //               += -1.00 d-_bb(j,j) l2_aaa(i,a,b)
    //               += -0.50 <i,k||c,d>_aaaa t2_aaaa(c,d,j,k) l2_1_aaa(j,a,b)
    //               += +1.00 d-_aa(i,j) l2_aaa(j,a,b)
    //               += -1.00 d-_aa(j,j) l2_aaa(i,a,b)
    //               += -1.00 <i,c||a,b>_aaaa l1_1_a(c)
    //               += -1.00 f_aa(i,j) l2_1_aaa(j,a,b)
    //               += -2.00 d-_aa(j,c) t1_1_aa(c,j) l2_1_aaa(i,a,b)
    //               += +1.00 l2_1_aaa(i,a,b) w0
    //               += -0.50 <i,k||c,d>_abab t2_abab(c,d,j,k) l2_1_aaa(j,a,b)
    //               += -0.50 <i,k||d,c>_abab t2_abab(d,c,j,k) l2_1_aaa(j,a,b)
    //               += +0.50 <i,k||a,b>_aaaa t2_abab(d,c,k,j) l2_1_bab(j,d,c)
    //               += +0.50 <i,k||a,b>_aaaa t2_abab(c,d,k,j) l2_1_bab(j,c,d)
    //               += +0.25 <k,j||a,b>_aaaa t2_aaaa(d,c,k,j) l2_1_aaa(i,d,c)
    //               += -0.50 <i,k||a,b>_aaaa t2_aaaa(d,c,j,k) l2_1_aaa(j,d,c)
    sigmal2_1_aaa("L,a,b,i")  = -2.00 * tmps_["21_aaa_Lovv"]("L,i,a,b");
    tmps_["21_aaa_Lovv"].~TArrayD();

    // flops: o1v0L1  = o2v2L1 o2v2L1 o1v0L1 o1v1L1 o1v0L1
    //  mems: o1v0L1  = o1v0L1 o1v0L1 o1v0L1 o1v0L1 o1v0L1
    tmps_["26_a_Lo"]("R,i")  = 2.00 * eri["abab_oovv"]("i,k,a,c") * r2["abb"]("R,a,c,k");
    tmps_["26_a_Lo"]("R,i") += eri["aaaa_oovv"]("i,j,a,b") * r2["aaa"]("R,a,b,j");
    tmps_["26_a_Lo"]("R,i") += 2.00 * f["aa_ov"]("i,a") * r1["a"]("R,a");

    // flops: o1v0L1  = o1v1L1
    //  mems: o1v0L1  = o1v0L1
    tmps_["25_a_Lo"]("R,i")  = r1["a"]("R,a") * reused_["4_aa_ov"]("i,a");

    // flops: o1v0L1  = o1v1L1
    //  mems: o1v0L1  = o1v0L1
    tmps_["24_a_Lo"]("R,i")  = dp["aa_ov"]("i,a") * r1_1["a"]("R,a");

    // flops: o1v2L1  = o2v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["22_aaa_Lovv"]("L,i,a,b")  = 0.25 * eri["aaaa_oovv"]("i,j,a,b") * tmps_["19_a_Lo"]("L,j");
    tmps_["22_aaa_Lovv"].~TArrayD();

    // flops: o0v1L1  = o0v1L1
    //  mems: o0v1L1  = o0v1L1
    tmps_["23_a_Lv"]("R,a")  = -1.00 * scalars_["1"] * r1["a"]("R,a");

    // flops: o0v1L1  = o1v1L1 o1v2L1 o0v2L1 o0v1L1 o1v2L1 o0v1L1 o0v2L1 o0v1L1 o0v1L1 o0v1L1 o0v2L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v2L1 o0v1L1 o1v2L1 o0v1L1 o0v1L1 o0v1L1 o1v2L1 o0v1L1 o1v2L1 o0v1L1 o0v2L1 o0v1L1 o0v1L1 o0v1L1 o1v2L1 o0v1L1 o1v3L1 o0v1L1 o1v2L1 o0v1L1 o1v2L1 o0v1L1 o1v3L1 o0v1L1 o0v2L1 o0v1L1 o0v1L1 o1v1L1 o0v1L1 o1v1L1 o0v1L1 o0v1L1
    //  mems: o0v1L1  = o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1
    tmps_["27_a_Lv"]("R,a")  = -0.25 * t1_1["aa"]("a,i") * tmps_["26_a_Lo"]("R,i");
    tmps_["27_a_Lv"]("R,a") += 0.50 * r2["abb"]("R,a,d,k") * reused_["12_bb_ov"]("k,d");
    tmps_["27_a_Lv"]("R,a") -= 0.50 * r1["a"]("R,c") * reused_["14_aa_vv"]("a,c");
    tmps_["27_a_Lv"]("R,a") -= 0.50 * r2["aaa"]("R,c,a,j") * reused_["9_aa_ov"]("j,c");
    tmps_["27_a_Lv"]("R,a") -= 0.50 * r1_1["a"]("R,c") * reused_["2_aa_vv"]("a,c");
    tmps_["27_a_Lv"]("R,a") += 0.50 * scalars_["6"] * r1["a"]("R,a");
    tmps_["27_a_Lv"]("R,a") -= 0.50 * r1["a"]("R,c") * reused_["13_aa_vv"]("a,c");
    tmps_["27_a_Lv"]("R,a") += 0.50 * w0 * r1_1["a"]("R,a");
    tmps_["27_a_Lv"]("R,a") -= scalars_["3"] * r1_1["a"]("R,a");
    tmps_["27_a_Lv"]("R,a") -= 0.50 * scalars_["5"] * r1_1["a"]("R,a");
    tmps_["27_a_Lv"]("R,a") += 0.50 * f["aa_vv"]("a,b") * r1_1["a"]("R,b");
    tmps_["27_a_Lv"]("R,a") += 0.50 * f["bb_ov"]("l,e") * r2_1["abb"]("R,a,e,l");
    tmps_["27_a_Lv"]("R,a") -= 0.50 * scalars_["2"] * r1["a"]("R,a");
    tmps_["27_a_Lv"]("R,a") -= 0.50 * r2_1["abb"]("R,a,e,l") * reused_["5_bb_ov"]("l,e");
    tmps_["27_a_Lv"]("R,a") += 0.50 * r2["abb"]("R,a,d,k") * reused_["10_bb_ov"]("k,d");
    tmps_["27_a_Lv"]("R,a") += 0.25 * r1_1["a"]("R,c") * reused_["1_aa_vv"]("a,c");
    tmps_["27_a_Lv"]("R,a") -= scalars_["4"] * r1_1["a"]("R,a");
    tmps_["27_a_Lv"]("R,a") += 0.50 * r2_1["aaa"]("R,b,a,i") * reused_["4_aa_ov"]("i,b");
    tmps_["27_a_Lv"]("R,a") += 0.50 * eri["abab_vovv"]("a,l,b,d") * r2_1["abb"]("R,b,d,l");
    tmps_["27_a_Lv"]("R,a") -= 0.50 * r2["aaa"]("R,c,a,j") * reused_["11_aa_ov"]("j,c");
    tmps_["27_a_Lv"]("R,a") -= 0.50 * f["aa_ov"]("i,b") * r2_1["aaa"]("R,b,a,i");
    tmps_["27_a_Lv"]("R,a") += 0.25 * eri["aaaa_vovv"]("a,i,b,c") * r2_1["aaa"]("R,b,c,i");
    tmps_["27_a_Lv"]("R,a") -= 0.50 * r1_1["a"]("R,b") * reused_["3_aa_vv"]("a,b");
    tmps_["27_a_Lv"]("R,a") += 0.50 * tmps_["23_a_Lv"]("R,a");
    tmps_["27_a_Lv"]("R,a") += t1_1["aa"]("a,i") * tmps_["24_a_Lo"]("R,i");
    tmps_["27_a_Lv"]("R,a") += 0.50 * t1_1["aa"]("a,i") * tmps_["25_a_Lo"]("R,i");
    tmps_["23_a_Lv"].~TArrayD();

    // sigmar1_1_a  = +2.00 d-_aa(i,b) r1_1_a(b) t1_1_aa(a,i)
    //             += +1.00 <i,j||b,c>_abab r2_abb(a,c,j) t1_1_aa(b,i)
    //             += +1.00 <i,a||b,c>_aaaa r1_a(c) t1_1_aa(b,i)
    //             += -0.50 <j,i||b,c>_aaaa r1_a(c) t2_1_aaaa(b,a,j,i)
    //             += +1.00 <j,i||b,c>_aaaa r2_aaa(c,a,j) t1_1_aa(b,i)
    //             += -0.50 <j,i||c,b>_abab r1_1_a(c) t2_abab(a,b,j,i)
    //             += -0.50 <i,j||c,b>_abab r1_1_a(c) t2_abab(a,b,i,j)
    //             += +1.00 r1_a(a) t0_1 w0
    //             += +0.25 <j,i||b,c>_abab r1_a(a) t2_1_abab(b,c,j,i)
    //             += +0.25 <i,j||b,c>_abab r1_a(a) t2_1_abab(b,c,i,j)
    //             += +0.25 <j,i||c,b>_abab r1_a(a) t2_1_abab(c,b,j,i)
    //             += +0.25 <i,j||c,b>_abab r1_a(a) t2_1_abab(c,b,i,j)
    //             += -1.00 d-_bb(i,b) r1_a(a) t0_1 t1_1_bb(b,i)
    //             += +0.25 <j,i||b,c>_bbbb r1_a(a) t2_1_bbbb(b,c,j,i)
    //             += -1.00 d-_aa(i,b) r1_a(a) t0_1 t1_1_aa(b,i)
    //             += +1.00 f_bb(i,b) r1_a(a) t1_1_bb(b,i)
    //             += +0.25 <j,i||b,c>_aaaa r1_a(a) t2_1_aaaa(b,c,j,i)
    //             += +1.00 f_aa(i,b) r1_a(a) t1_1_aa(b,i)
    //             += -0.50 <j,i||c,b>_abab r1_a(c) t2_1_abab(a,b,j,i)
    //             += -0.50 <i,j||c,b>_abab r1_a(c) t2_1_abab(a,b,i,j)
    //             += +1.00 <a,i||c,b>_abab r1_a(c) t1_1_bb(b,i)
    //             += +1.00 r1_1_a(a) w0
    //             += -2.00 d-_bb(i,b) r1_1_a(a) t1_1_bb(b,i)
    //             += -0.50 <j,i||j,i>_abab r1_1_a(a)
    //             += -0.50 <i,j||i,j>_abab r1_1_a(a)
    //             += +0.25 <j,i||b,c>_bbbb r1_1_a(a) t2_bbbb(b,c,j,i)
    //             += +1.00 f_bb(i,i) r1_1_a(a)
    //             += +0.25 <j,i||b,c>_aaaa r1_1_a(a) t2_aaaa(b,c,j,i)
    //             += -1.00 d-_aa(i,i) r1_1_a(a) t0_1
    //             += +1.00 f_aa(i,i) r1_1_a(a)
    //             += -0.50 <j,i||j,i>_bbbb r1_1_a(a)
    //             += -0.50 <j,i||j,i>_aaaa r1_1_a(a)
    //             += +0.25 <j,i||b,c>_abab r1_1_a(a) t2_abab(b,c,j,i)
    //             += +0.25 <i,j||b,c>_abab r1_1_a(a) t2_abab(b,c,i,j)
    //             += +0.25 <j,i||c,b>_abab r1_1_a(a) t2_abab(c,b,j,i)
    //             += +0.25 <i,j||c,b>_abab r1_1_a(a) t2_abab(c,b,i,j)
    //             += -1.00 d-_bb(i,i) r1_1_a(a) t0_1
    //             += +1.00 f_aa(a,b) r1_1_a(b)
    //             += +1.00 f_bb(i,b) r2_1_abb(a,b,i)
    //             += -1.00 d+_aa(i,i) r1_a(a)
    //             += -1.00 d-_bb(i,b) r2_1_abb(a,b,i) t0_1
    //             += -1.00 <j,i||b,c>_bbbb r2_abb(a,c,j) t1_1_bb(b,i)
    //             += -0.50 <j,i||b,c>_aaaa r1_1_a(c) t2_aaaa(b,a,j,i)
    //             += -2.00 d-_aa(i,b) r1_1_a(a) t1_1_aa(b,i)
    //             += +1.00 d-_aa(i,b) r2_1_aaa(b,a,i) t0_1
    //             += +0.50 <a,i||b,c>_abab r2_1_abb(b,c,i)
    //             += +0.50 <a,i||c,b>_abab r2_1_abb(c,b,i)
    //             += -1.00 <j,i||c,b>_abab r2_aaa(c,a,j) t1_1_bb(b,i)
    //             += -1.00 f_aa(i,b) r2_1_aaa(b,a,i)
    //             += -0.50 <i,a||b,c>_aaaa r2_1_aaa(b,c,i)
    //             += -1.00 d-_aa(a,b) r1_1_a(b) t0_1
    //             += -1.00 d+_bb(i,i) r1_a(a)
    //             += +1.00 d-_aa(i,b) r1_a(b) t0_1 t1_1_aa(a,i)
    //             += +0.50 <j,i||b,c>_aaaa r2_aaa(b,c,j) t1_1_aa(a,i)
    //             += -0.50 <i,j||b,c>_abab r2_abb(b,c,j) t1_1_aa(a,i)
    //             += -0.50 <i,j||c,b>_abab r2_abb(c,b,j) t1_1_aa(a,i)
    //             += -1.00 f_aa(i,b) r1_a(b) t1_1_aa(a,i)
    sigmar1_1_a("R,a")  = 2.00 * tmps_["27_a_Lv"]("R,a");
    tmps_["27_a_Lv"].~TArrayD();

    // flops: o0v1L1  = o1v1L1
    //  mems: o0v1L1  = o0v1L1
    tmps_["28_a_Lv"]("R,a")  = -0.25 * t1_1["aa"]("a,i") * tmps_["26_a_Lo"]("R,i");
    tmps_["28_a_Lv"].~TArrayD();

    // flops: o1v2L1  = o1v1L1 o2v2L1
    //  mems: o1v2L1  = o1v0L1 o1v2L1
    tmps_["29_aaa_Lvvo"]("R,a,b,i")  = dp["aa_ov"]("j,c") * r1["a"]("R,c") * t2["aaaa"]("a,b,i,j");

    // sigmar2_1_aaa  = -1.00 d+_aa(j,c) r1_a(c) t2_aaaa(a,b,i,j)
    sigmar2_1_aaa("R,a,b,i")  = -1.00 * tmps_["29_aaa_Lvvo"]("R,a,b,i");

    // flops: o3v0L1  = o3v2L1
    //  mems: o3v0L1  = o3v0L1
    tmps_["30_aaa_Looo"]("R,i,j,k")  = r2["aaa"]("R,a,b,i") * eri["aaaa_oovv"]("j,k,a,b");

    // flops: o1v2L1  = o1v4L1 o2v2L1 o2v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o3v2L1 o1v2L1 o2v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["31_aaa_Lvvo"]("R,a,b,i")  = -0.50 * eri["aaaa_vvvv"]("a,b,c,d") * r2["aaa"]("R,c,d,i");
    tmps_["31_aaa_Lvvo"]("R,a,b,i") -= 0.50 * r2["aaa"]("R,a,b,k") * reused_["16_aa_oo"]("i,k");
    tmps_["31_aaa_Lvvo"]("R,a,b,i") += f["aa_oo"]("j,i") * r2["aaa"]("R,a,b,j");
    tmps_["31_aaa_Lvvo"]("R,a,b,i") += r2["aaa"]("R,a,b,k") * reused_["6_aa_oo"]("k,i");
    tmps_["31_aaa_Lvvo"]("R,a,b,i") += scalars_["1"] * r2_1["aaa"]("R,a,b,i");
    tmps_["31_aaa_Lvvo"]("R,a,b,i") += scalars_["5"] * r2["aaa"]("R,a,b,i");
    tmps_["31_aaa_Lvvo"]("R,a,b,i") -= r2["aaa"]("R,a,b,j") * reused_["8_aa_oo"]("j,i");
    tmps_["31_aaa_Lvvo"]("R,a,b,i") += scalars_["3"] * r2["aaa"]("R,a,b,i");
    tmps_["31_aaa_Lvvo"]("R,a,b,i") += 0.50 * r1["a"]("R,c") * reused_["15_aaaa_vovv"]("c,i,a,b");
    tmps_["31_aaa_Lvvo"]("R,a,b,i") -= eri["aaaa_vvvo"]("a,b,c,i") * r1["a"]("R,c");
    tmps_["31_aaa_Lvvo"]("R,a,b,i") += scalars_["4"] * r2["aaa"]("R,a,b,i");
    tmps_["31_aaa_Lvvo"]("R,a,b,i") -= dp["aa_oo"]("j,i") * r2_1["aaa"]("R,a,b,j");
    tmps_["31_aaa_Lvvo"]("R,a,b,i") += scalars_["2"] * r2_1["aaa"]("R,a,b,i");
    tmps_["31_aaa_Lvvo"]("R,a,b,i") += t0_1 * tmps_["29_aaa_Lvvo"]("R,a,b,i");
    tmps_["31_aaa_Lvvo"]("R,a,b,i") += 0.25 * t2["aaaa"]("a,b,k,j") * tmps_["30_aaa_Looo"]("R,i,j,k");
    tmps_["31_aaa_Lvvo"]("R,a,b,i") -= 0.50 * t2["aaaa"]("a,b,i,j") * tmps_["26_a_Lo"]("R,j");
    tmps_["29_aaa_Lvvo"].~TArrayD();

    // sigmar2_aaa  = -1.00 d-_aa(j,c) r1_a(c) t2_aaaa(a,b,i,j) t0_1
    //             += -1.00 d-_bb(j,j) r2_1_aaa(a,b,i)
    //             += -0.50 <k,j||c,d>_abab r2_aaa(a,b,k) t2_abab(c,d,i,j)
    //             += -0.50 <k,j||d,c>_abab r2_aaa(a,b,k) t2_abab(d,c,i,j)
    //             += -1.00 f_aa(j,i) r2_aaa(a,b,j)
    //             += -0.50 <k,j||c,d>_aaaa r2_aaa(a,b,k) t2_aaaa(c,d,i,j)
    //             += -0.50 <k,j||k,j>_abab r2_aaa(a,b,i)
    //             += -0.50 <j,k||j,k>_abab r2_aaa(a,b,i)
    //             += +0.25 <k,j||c,d>_bbbb r2_aaa(a,b,i) t2_bbbb(c,d,k,j)
    //             += +1.00 f_bb(j,j) r2_aaa(a,b,i)
    //             += +0.25 <k,j||c,d>_aaaa r2_aaa(a,b,i) t2_aaaa(c,d,k,j)
    //             += -1.00 d-_aa(j,j) r2_aaa(a,b,i) t0_1
    //             += +1.00 f_aa(j,j) r2_aaa(a,b,i)
    //             += -0.50 <k,j||k,j>_bbbb r2_aaa(a,b,i)
    //             += -0.50 <k,j||k,j>_aaaa r2_aaa(a,b,i)
    //             += +0.25 <k,j||c,d>_abab r2_aaa(a,b,i) t2_abab(c,d,k,j)
    //             += +0.25 <j,k||c,d>_abab r2_aaa(a,b,i) t2_abab(c,d,j,k)
    //             += +0.25 <k,j||d,c>_abab r2_aaa(a,b,i) t2_abab(d,c,k,j)
    //             += +0.25 <j,k||d,c>_abab r2_aaa(a,b,i) t2_abab(d,c,j,k)
    //             += -1.00 d-_bb(j,j) r2_aaa(a,b,i) t0_1
    //             += +1.00 d-_aa(j,i) r2_aaa(a,b,j) t0_1
    //             += -1.00 d-_bb(j,c) r2_aaa(a,b,i) t1_1_bb(c,j)
    //             += +0.50 <k,j||c,i>_aaaa r1_a(c) t2_aaaa(a,b,k,j)
    //             += +1.00 <a,b||c,i>_aaaa r1_a(c)
    //             += +0.50 <a,b||c,d>_aaaa r2_aaa(c,d,i)
    //             += -1.00 d-_aa(j,c) r2_aaa(a,b,i) t1_1_aa(c,j)
    //             += +1.00 d-_aa(j,i) r2_1_aaa(a,b,j)
    //             += -1.00 d-_aa(j,j) r2_1_aaa(a,b,i)
    //             += +0.25 <k,j||c,d>_aaaa r2_aaa(c,d,i) t2_aaaa(a,b,k,j)
    //             += -0.50 <k,j||c,d>_aaaa r2_aaa(c,d,k) t2_aaaa(a,b,i,j)
    //             += +0.50 <j,k||c,d>_abab r2_abb(c,d,k) t2_aaaa(a,b,i,j)
    //             += +0.50 <j,k||d,c>_abab r2_abb(d,c,k) t2_aaaa(a,b,i,j)
    //             += +1.00 f_aa(j,c) r1_a(c) t2_aaaa(a,b,i,j)
    sigmar2_aaa("R,a,b,i")  = -1.00 * tmps_["31_aaa_Lvvo"]("R,a,b,i");
    tmps_["31_aaa_Lvvo"].~TArrayD();

    // flops: o1v2L1  = o2v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["32_aaa_Lvvo"]("R,a,b,i")  = -0.50 * t2["aaaa"]("a,b,i,j") * tmps_["26_a_Lo"]("R,j");
    tmps_["32_aaa_Lvvo"].~TArrayD();

    // flops: o1v0L1  = o2v2L1
    //  mems: o1v0L1  = o1v0L1
    tmps_["33_a_Lo"]("L,i")  = l2_1["bab"]("L,j,a,b") * t2_1["abab"]("a,b,i,j");

    // flops: o1v2L1  = o2v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["34_aaa_Lovv"]("L,i,a,b")  = -4.00 * eri["aaaa_oovv"]("i,j,a,b") * tmps_["33_a_Lo"]("L,j");
    tmps_["34_aaa_Lovv"].~TArrayD();

    // flops: o1v2L1  = o2v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["35_aaa_Lvvo"]("L,a,b,i")  = l2_1["aaa"]("L,j,a,b") * reused_["18_aa_oo"]("i,j");

    // sigmal2_1_aaa += +2.00 d-_aa(i,c) t1_1_aa(c,j) l2_1_aaa(j,a,b)
    sigmal2_1_aaa("L,a,b,i") += 2.00 * tmps_["35_aaa_Lvvo"]("L,a,b,i");

    // flops: o1v0L1  = o1v1L1
    //  mems: o1v0L1  = o1v0L1
    tmps_["42_a_Lo"]("L,i")  = l1_1["a"]("L,a") * t1_1["aa"]("a,i");

    // flops: o1v0L1  = o2v2L1
    //  mems: o1v0L1  = o1v0L1
    tmps_["41_a_Lo"]("L,i")  = l2["bab"]("L,j,a,b") * t2["abab"]("a,b,i,j");

    // flops: o2v1L1  = o2v2L1
    //  mems: o2v1L1  = o2v1L1
    tmps_["40_aaa_Lovo"]("L,i,a,j")  = l2_1["aaa"]("L,i,a,b") * t1_1["aa"]("b,j");

    // flops: o1v0L1  = o2v2L1
    //  mems: o1v0L1  = o1v0L1
    tmps_["39_a_Lo"]("L,i")  = l2["aaa"]("L,j,a,b") * t2["aaaa"]("a,b,j,i");

    // flops: o1v0L1  = o2v2L1
    //  mems: o1v0L1  = o1v0L1
    tmps_["38_a_Lo"]("L,i")  = l2_1["aaa"]("L,j,a,b") * t2_1["aaaa"]("a,b,j,i");

    // flops: o1v2L1  = o1v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["37_aaa_Lovv"]("L,i,a,b")  = -1.00 * scalars_["4"] * l2["aaa"]("L,i,a,b");

    // flops: o3v0L1  = o3v2L1 o3v2L1 o3v0L1
    //  mems: o3v0L1  = o3v0L1 o3v0L1 o3v0L1
    tmps_["36_aaa_Looo"]("L,i,j,k")  = l2["aaa"]("L,i,a,b") * t2["aaaa"]("a,b,j,k");
    tmps_["36_aaa_Looo"]("L,i,j,k") += l2_1["aaa"]("L,i,a,b") * t2_1["aaaa"]("a,b,j,k");

    // flops: o1v2L1  = o1v2L1 o2v2L1 o2v3L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o1v3L1 o1v2L1 o1v4L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o3v2L1 o1v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["43_aaa_Lvvo"]("L,a,b,i")  = -4.00 * t0_1 * tmps_["35_aaa_Lvvo"]("L,a,b,i");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") -= 4.00 * eri["aaaa_oovv"]("i,j,a,b") * tmps_["41_a_Lo"]("L,j");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") += 4.00 * eri["aaaa_vovv"]("c,l,a,b") * tmps_["40_aaa_Lovo"]("L,i,c,l");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") -= 4.00 * scalars_["6"] * l2_1["aaa"]("L,i,a,b");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") -= 4.00 * l2_1["aaa"]("L,l,a,b") * dp["aa_oo"]("i,l");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") += 4.00 * scalars_["3"] * l2["aaa"]("L,i,a,b");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") += 2.00 * l2["aaa"]("L,l,a,b") * reused_["7_aa_oo"]("i,l");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") += 4.00 * l2["aaa"]("L,l,a,b") * f["aa_oo"]("i,l");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") += 4.00 * l2["aaa"]("L,l,a,b") * reused_["6_aa_oo"]("i,l");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") -= 4.00 * eri["aaaa_vovv"]("e,i,a,b") * l1["a"]("L,e");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") += 2.00 * eri["aaaa_vvvv"]("e,c,a,b") * l2["aaa"]("L,i,c,e");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") += 4.00 * scalars_["1"] * l2_1["aaa"]("L,i,a,b");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") += 4.00 * scalars_["5"] * l2["aaa"]("L,i,a,b");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") += 4.00 * scalars_["2"] * l2_1["aaa"]("L,i,a,b");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") += 4.00 * l2_1["aaa"]("L,l,a,b") * reused_["19_aa_oo"]("i,l");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") -= 4.00 * l2["aaa"]("L,l,a,b") * reused_["8_aa_oo"]("i,l");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") -= 4.00 * l2_1["aaa"]("L,l,a,b") * reused_["20_aa_oo"]("i,l");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") -= 4.00 * l2_1["aaa"]("L,l,a,b") * reused_["17_aa_oo"]("i,l");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") -= 4.00 * l2["aaa"]("L,l,a,b") * reused_["18_aa_oo"]("i,l");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") -= 4.00 * tmps_["37_aaa_Lovv"]("L,i,a,b");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") += eri["aaaa_oovv"]("l,j,a,b") * tmps_["36_aaa_Looo"]("L,i,j,l");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") += 2.00 * eri["aaaa_oovv"]("i,j,a,b") * tmps_["39_a_Lo"]("L,j");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") -= 4.00 * eri["aaaa_oovv"]("i,l,a,b") * tmps_["42_a_Lo"]("L,l");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") += 2.00 * eri["aaaa_oovv"]("i,j,a,b") * tmps_["38_a_Lo"]("L,j");
    tmps_["43_aaa_Lvvo"]("L,a,b,i") -= 4.00 * eri["aaaa_oovv"]("i,j,a,b") * tmps_["33_a_Lo"]("L,j");
    tmps_["37_aaa_Lovv"].~TArrayD();
    tmps_["36_aaa_Looo"].~TArrayD();
    tmps_["35_aaa_Lvvo"].~TArrayD();

    // sigmal2_aaa  = +0.25 <k,j||a,b>_aaaa t2_aaaa(d,c,k,j) l2_aaa(i,d,c)
    //             += +0.25 <k,j||a,b>_aaaa t2_1_aaaa(d,c,k,j) l2_1_aaa(i,d,c)
    //             += +1.00 t0_1 l2_1_aaa(i,a,b) w0
    //             += +0.25 <k,j||c,d>_abab t2_1_abab(c,d,k,j) l2_1_aaa(i,a,b)
    //             += +0.25 <j,k||c,d>_abab t2_1_abab(c,d,j,k) l2_1_aaa(i,a,b)
    //             += +0.25 <k,j||d,c>_abab t2_1_abab(d,c,k,j) l2_1_aaa(i,a,b)
    //             += +0.25 <j,k||d,c>_abab t2_1_abab(d,c,j,k) l2_1_aaa(i,a,b)
    //             += -1.00 d-_bb(j,c) t0_1 t1_1_bb(c,j) l2_1_aaa(i,a,b)
    //             += +0.25 <k,j||c,d>_bbbb t2_1_bbbb(c,d,k,j) l2_1_aaa(i,a,b)
    //             += -1.00 d-_aa(j,c) t0_1 t1_1_aa(c,j) l2_1_aaa(i,a,b)
    //             += +1.00 f_bb(j,c) t1_1_bb(c,j) l2_1_aaa(i,a,b)
    //             += +0.25 <k,j||c,d>_aaaa t2_1_aaaa(c,d,k,j) l2_1_aaa(i,a,b)
    //             += +1.00 f_aa(j,c) t1_1_aa(c,j) l2_1_aaa(i,a,b)
    //             += +1.00 d+_aa(i,j) l2_1_aaa(j,a,b)
    //             += -1.00 d-_bb(j,c) t1_1_bb(c,j) l2_aaa(i,a,b)
    //             += -0.50 <i,k||c,d>_aaaa t2_aaaa(c,d,j,k) l2_aaa(j,a,b)
    //             += -1.00 f_aa(i,j) l2_aaa(j,a,b)
    //             += -0.50 <i,k||c,d>_abab t2_abab(c,d,j,k) l2_aaa(j,a,b)
    //             += -0.50 <i,k||d,c>_abab t2_abab(d,c,j,k) l2_aaa(j,a,b)
    //             += -1.00 <i,c||a,b>_aaaa l1_a(c)
    //             += +0.50 <d,c||a,b>_aaaa l2_aaa(i,d,c)
    //             += -1.00 d+_bb(j,j) l2_1_aaa(i,a,b)
    //             += -0.50 <k,j||k,j>_abab l2_aaa(i,a,b)
    //             += -0.50 <j,k||j,k>_abab l2_aaa(i,a,b)
    //             += +0.25 <k,j||c,d>_bbbb t2_bbbb(c,d,k,j) l2_aaa(i,a,b)
    //             += +1.00 f_bb(j,j) l2_aaa(i,a,b)
    //             += +0.25 <k,j||c,d>_aaaa t2_aaaa(c,d,k,j) l2_aaa(i,a,b)
    //             += -1.00 d-_aa(j,j) t0_1 l2_aaa(i,a,b)
    //             += +1.00 f_aa(j,j) l2_aaa(i,a,b)
    //             += -0.50 <k,j||k,j>_bbbb l2_aaa(i,a,b)
    //             += -0.50 <k,j||k,j>_aaaa l2_aaa(i,a,b)
    //             += +0.25 <k,j||c,d>_abab t2_abab(c,d,k,j) l2_aaa(i,a,b)
    //             += +0.25 <j,k||c,d>_abab t2_abab(c,d,j,k) l2_aaa(i,a,b)
    //             += +0.25 <k,j||d,c>_abab t2_abab(d,c,k,j) l2_aaa(i,a,b)
    //             += +0.25 <j,k||d,c>_abab t2_abab(d,c,j,k) l2_aaa(i,a,b)
    //             += -1.00 d-_bb(j,j) t0_1 l2_aaa(i,a,b)
    //             += -1.00 d+_aa(j,j) l2_1_aaa(i,a,b)
    //             += -1.00 f_aa(i,c) t1_1_aa(c,j) l2_1_aaa(j,a,b)
    //             += +1.00 d-_aa(i,j) t0_1 l2_aaa(j,a,b)
    //             += -1.00 <i,k||j,c>_abab t1_1_bb(c,k) l2_1_aaa(j,a,b)
    //             += -0.50 <i,k||c,d>_abab t2_1_abab(c,d,j,k) l2_1_aaa(j,a,b)
    //             += -0.50 <i,k||d,c>_abab t2_1_abab(d,c,j,k) l2_1_aaa(j,a,b)
    //             += +1.00 <i,k||c,j>_aaaa t1_1_aa(c,k) l2_1_aaa(j,a,b)
    //             += -0.50 <i,k||c,d>_aaaa t2_1_aaaa(c,d,j,k) l2_1_aaa(j,a,b)
    //             += +1.00 d-_aa(i,c) t1_1_aa(c,j) l2_aaa(j,a,b)
    //             += -1.00 d-_aa(j,c) t1_1_aa(c,j) l2_aaa(i,a,b)
    //             += +1.00 <j,d||a,b>_aaaa t1_1_aa(c,j) l2_1_aaa(i,d,c)
    //             += +0.50 <i,k||a,b>_aaaa t2_abab(d,c,k,j) l2_bab(j,d,c)
    //             += +0.50 <i,k||a,b>_aaaa t2_abab(c,d,k,j) l2_bab(j,c,d)
    //             += -0.50 <i,k||a,b>_aaaa t2_aaaa(d,c,j,k) l2_aaa(j,d,c)
    //             += +1.00 d-_aa(i,c) t0_1 t1_1_aa(c,j) l2_1_aaa(j,a,b)
    //             += +1.00 <i,j||a,b>_aaaa t1_1_aa(c,j) l1_1_a(c)
    //             += -0.50 <i,k||a,b>_aaaa t2_1_aaaa(d,c,j,k) l2_1_aaa(j,d,c)
    //             += +0.50 <i,k||a,b>_aaaa t2_1_abab(d,c,k,j) l2_1_bab(j,d,c)
    //             += +0.50 <i,k||a,b>_aaaa t2_1_abab(c,d,k,j) l2_1_bab(j,c,d)
    sigmal2_aaa("L,a,b,i")  = -0.25 * tmps_["43_aaa_Lvvo"]("L,a,b,i");
    tmps_["43_aaa_Lvvo"].~TArrayD();

    // flops: o1v2L1  = o2v2L1 o1v1L1 o2v2L1 o1v2L1 o1v1L1 o2v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v0L1 o1v2L1 o1v2L1 o1v0L1 o1v2L1 o1v2L1
    tmps_["45_aaa_Lvvo"]("R,a,b,i")  = -1.00 * r2["aaa"]("R,a,b,j") * reused_["18_aa_oo"]("j,i");
    tmps_["45_aaa_Lvvo"]("R,a,b,i") += dp["aa_ov"]("j,c") * r1["a"]("R,c") * t2_1["aaaa"]("a,b,i,j");
    tmps_["45_aaa_Lvvo"]("R,a,b,i") += dp["aa_ov"]("j,c") * r1_1["a"]("R,c") * t2["aaaa"]("a,b,i,j");

    // sigmar2_aaa += -1.00 d-_aa(j,c) r1_a(c) t2_1_aaaa(a,b,i,j)
    //             += +1.00 d-_aa(j,c) r2_aaa(a,b,j) t1_1_aa(c,i)
    //             += -1.00 d-_aa(j,c) r1_1_a(c) t2_aaaa(a,b,i,j)
    sigmar2_aaa("R,a,b,i") -= tmps_["45_aaa_Lvvo"]("R,a,b,i");

    // flops: o1v0L1  = o2v2L1 o1v1L1 o1v1L1 o1v0L1 o1v0L1 o1v1L1 o1v0L1 o2v2L1 o1v0L1
    //  mems: o1v0L1  = o1v0L1 o1v0L1 o1v0L1 o1v0L1 o1v0L1 o1v0L1 o1v0L1 o1v0L1 o1v0L1
    tmps_["47_a_Lo"]("R,i")  = -0.50 * eri["aaaa_oovv"]("i,k,a,c") * r2_1["aaa"]("R,a,c,k");
    tmps_["47_a_Lo"]("R,i") -= f["aa_ov"]("i,a") * r1_1["a"]("R,a");
    tmps_["47_a_Lo"]("R,i") += r1["a"]("R,c") * reused_["24_aa_ov"]("i,c");
    tmps_["47_a_Lo"]("R,i") -= r1["a"]("R,c") * reused_["11_aa_ov"]("i,c");
    tmps_["47_a_Lo"]("R,i") -= eri["abab_oovv"]("i,j,a,b") * r2_1["abb"]("R,a,b,j");

    // flops: o1v2L1  = o2v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["46_aaa_Lvvo"]("R,a,b,i")  = -0.50 * r2["aaa"]("R,a,b,j") * reused_["22_aa_oo"]("i,j");

    // flops: o3v0L1  = o3v2L1
    //  mems: o3v0L1  = o3v0L1
    tmps_["44_aaa_Looo"]("R,i,j,k")  = r2_1["aaa"]("R,a,b,i") * eri["aaaa_oovv"]("j,k,a,b");

    // flops: o1v2L1  = o2v2L1 o3v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o1v4L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o3v2L1 o1v2L1 o2v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["48_aaa_Lvvo"]("R,a,b,i")  = -2.00 * t2_1["aaaa"]("a,b,i,k") * tmps_["26_a_Lo"]("R,k");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") += t2["aaaa"]("a,b,j,k") * tmps_["44_aaa_Looo"]("R,i,k,j");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") += 8.00 * t2_1["aaaa"]("a,b,i,k") * tmps_["24_a_Lo"]("R,k");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") += 4.00 * r2_1["aaa"]("R,a,b,j") * reused_["6_aa_oo"]("j,i");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") += 8.00 * scalars_["3"] * r2_1["aaa"]("R,a,b,i");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") += 8.00 * scalars_["4"] * r2_1["aaa"]("R,a,b,i");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") += 4.00 * r1["a"]("R,d") * reused_["21_aaaa_vvvo"]("a,b,d,i");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") += 4.00 * scalars_["5"] * r2_1["aaa"]("R,a,b,i");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") += 2.00 * r1["a"]("R,c") * reused_["23_aaaa_vovv"]("c,i,a,b");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") -= 4.00 * r2["aaa"]("R,a,b,j") * reused_["20_aa_oo"]("j,i");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") -= 4.00 * scalars_["6"] * r2["aaa"]("R,a,b,i");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") -= 4.00 * w0 * r2_1["aaa"]("R,a,b,i");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") += 2.00 * r1_1["a"]("R,c") * reused_["15_aaaa_vovv"]("c,i,a,b");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") -= 4.00 * eri["aaaa_vvvo"]("a,b,c,i") * r1_1["a"]("R,c");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") += 4.00 * scalars_["1"] * r2["aaa"]("R,a,b,i");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") -= 8.00 * r2_1["aaa"]("R,a,b,k") * reused_["18_aa_oo"]("k,i");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") -= 2.00 * eri["aaaa_vvvv"]("a,b,c,d") * r2_1["aaa"]("R,c,d,i");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") += 4.00 * r2_1["aaa"]("R,a,b,k") * f["aa_oo"]("k,i");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") += 4.00 * r2["aaa"]("R,a,b,k") * reused_["19_aa_oo"]("k,i");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") += 4.00 * scalars_["2"] * r2["aaa"]("R,a,b,i");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") -= 4.00 * r2["aaa"]("R,a,b,k") * dp["aa_oo"]("k,i");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") -= 2.00 * r2_1["aaa"]("R,a,b,j") * reused_["16_aa_oo"]("i,j");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") -= 4.00 * r2_1["aaa"]("R,a,b,k") * reused_["8_aa_oo"]("k,i");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") += 4.00 * tmps_["46_aaa_Lvvo"]("R,a,b,i");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") += 4.00 * t0_1 * tmps_["45_aaa_Lvvo"]("R,a,b,i");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") += t2_1["aaaa"]("a,b,j,k") * tmps_["30_aaa_Looo"]("R,i,k,j");
    tmps_["48_aaa_Lvvo"]("R,a,b,i") += 4.00 * t2["aaaa"]("a,b,i,k") * tmps_["47_a_Lo"]("R,k");
    tmps_["46_aaa_Lvvo"].~TArrayD();
    tmps_["45_aaa_Lvvo"].~TArrayD();
    tmps_["44_aaa_Looo"].~TArrayD();
    tmps_["30_aaa_Looo"].~TArrayD();
}