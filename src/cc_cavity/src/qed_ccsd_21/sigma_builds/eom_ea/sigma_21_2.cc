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

void hilbert::EOM_EA_QED_CCSD_21::sigma_ea_21_2() {

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


    // sigmar2_1_aaa += +0.25 <k,j||c,d>_aaaa r2_1_aaa(c,d,i) t2_aaaa(a,b,k,j)
    //               += -0.50 <k,j||c,d>_aaaa r2_aaa(c,d,k) t2_1_aaaa(a,b,i,j)
    //               += +0.50 <j,k||c,d>_abab r2_abb(c,d,k) t2_1_aaaa(a,b,i,j)
    //               += +0.50 <j,k||d,c>_abab r2_abb(d,c,k) t2_1_aaaa(a,b,i,j)
    //               += +1.00 f_aa(j,c) r1_a(c) t2_1_aaaa(a,b,i,j)
    //               += -2.00 d-_aa(j,c) r1_1_a(c) t2_1_aaaa(a,b,i,j)
    //               += -0.50 <k,j||c,d>_abab r2_1_aaa(a,b,k) t2_abab(c,d,i,j)
    //               += -0.50 <k,j||d,c>_abab r2_1_aaa(a,b,k) t2_abab(d,c,i,j)
    //               += -2.00 d-_bb(j,c) r2_1_aaa(a,b,i) t1_1_bb(c,j)
    //               += -2.00 d-_aa(j,c) r2_1_aaa(a,b,i) t1_1_aa(c,j)
    //               += -1.00 <a,b||c,d>_aaaa r1_a(d) t1_1_aa(c,i)
    //               += -0.50 <k,j||c,d>_aaaa r1_a(d) t2_aaaa(a,b,k,j) t1_1_aa(c,i)
    //               += -0.50 <k,j||k,j>_abab r2_1_aaa(a,b,i)
    //               += -0.50 <j,k||j,k>_abab r2_1_aaa(a,b,i)
    //               += +0.25 <k,j||c,d>_bbbb r2_1_aaa(a,b,i) t2_bbbb(c,d,k,j)
    //               += +1.00 f_bb(j,j) r2_1_aaa(a,b,i)
    //               += +0.25 <k,j||c,d>_aaaa r2_1_aaa(a,b,i) t2_aaaa(c,d,k,j)
    //               += -1.00 d-_aa(j,j) r2_1_aaa(a,b,i) t0_1
    //               += +1.00 f_aa(j,j) r2_1_aaa(a,b,i)
    //               += -0.50 <k,j||k,j>_bbbb r2_1_aaa(a,b,i)
    //               += -0.50 <k,j||k,j>_aaaa r2_1_aaa(a,b,i)
    //               += +0.25 <k,j||c,d>_abab r2_1_aaa(a,b,i) t2_abab(c,d,k,j)
    //               += +0.25 <j,k||c,d>_abab r2_1_aaa(a,b,i) t2_abab(c,d,j,k)
    //               += +0.25 <k,j||d,c>_abab r2_1_aaa(a,b,i) t2_abab(d,c,k,j)
    //               += +0.25 <j,k||d,c>_abab r2_1_aaa(a,b,i) t2_abab(d,c,j,k)
    //               += -1.00 d-_bb(j,j) r2_1_aaa(a,b,i) t0_1
    //               += +0.50 <k,j||c,i>_aaaa r1_a(c) t2_1_aaaa(a,b,k,j)
    //               += -1.00 <k,j||i,c>_abab r2_aaa(a,b,k) t1_1_bb(c,j)
    //               += -0.50 <k,j||c,d>_abab r2_aaa(a,b,k) t2_1_abab(c,d,i,j)
    //               += -0.50 <k,j||d,c>_abab r2_aaa(a,b,k) t2_1_abab(d,c,i,j)
    //               += +1.00 r2_aaa(a,b,i) t0_1 w0
    //               += +0.25 <k,j||c,d>_abab r2_aaa(a,b,i) t2_1_abab(c,d,k,j)
    //               += +0.25 <j,k||c,d>_abab r2_aaa(a,b,i) t2_1_abab(c,d,j,k)
    //               += +0.25 <k,j||d,c>_abab r2_aaa(a,b,i) t2_1_abab(d,c,k,j)
    //               += +0.25 <j,k||d,c>_abab r2_aaa(a,b,i) t2_1_abab(d,c,j,k)
    //               += -1.00 d-_bb(j,c) r2_aaa(a,b,i) t0_1 t1_1_bb(c,j)
    //               += +0.25 <k,j||c,d>_bbbb r2_aaa(a,b,i) t2_1_bbbb(c,d,k,j)
    //               += -1.00 d-_aa(j,c) r2_aaa(a,b,i) t0_1 t1_1_aa(c,j)
    //               += +1.00 f_bb(j,c) r2_aaa(a,b,i) t1_1_bb(c,j)
    //               += +0.25 <k,j||c,d>_aaaa r2_aaa(a,b,i) t2_1_aaaa(c,d,k,j)
    //               += +1.00 f_aa(j,c) r2_aaa(a,b,i) t1_1_aa(c,j)
    //               += +1.00 r2_1_aaa(a,b,i) w0
    //               += +0.50 <k,j||c,i>_aaaa r1_1_a(c) t2_aaaa(a,b,k,j)
    //               += +1.00 <a,b||c,i>_aaaa r1_1_a(c)
    //               += -1.00 d+_bb(j,j) r2_aaa(a,b,i)
    //               += +2.00 d-_aa(j,c) r2_1_aaa(a,b,j) t1_1_aa(c,i)
    //               += +0.50 <a,b||c,d>_aaaa r2_1_aaa(c,d,i)
    //               += -1.00 f_aa(j,i) r2_1_aaa(a,b,j)
    //               += -1.00 f_aa(j,c) r2_aaa(a,b,j) t1_1_aa(c,i)
    //               += -1.00 d+_aa(j,j) r2_aaa(a,b,i)
    //               += +1.00 d+_aa(j,i) r2_aaa(a,b,j)
    //               += -0.50 <k,j||c,d>_aaaa r2_1_aaa(a,b,k) t2_aaaa(c,d,i,j)
    //               += +1.00 d-_aa(j,i) r2_1_aaa(a,b,j) t0_1
    //               += -0.50 <k,j||c,d>_aaaa r2_aaa(a,b,k) t2_1_aaaa(c,d,i,j)
    //               += +1.00 <k,j||c,i>_aaaa r2_aaa(a,b,k) t1_1_aa(c,j)
    //               += -1.00 d-_aa(j,c) r1_a(c) t0_1 t2_1_aaaa(a,b,i,j)
    //               += +1.00 d-_aa(j,c) r2_aaa(a,b,j) t0_1 t1_1_aa(c,i)
    //               += -1.00 d-_aa(j,c) r1_1_a(c) t2_aaaa(a,b,i,j) t0_1
    //               += +0.25 <k,j||c,d>_aaaa r2_aaa(c,d,i) t2_1_aaaa(a,b,k,j)
    //               += +1.00 <k,j||c,d>_aaaa r1_a(d) t2_aaaa(a,b,i,j) t1_1_aa(c,k)
    //               += +1.00 f_aa(j,c) r1_1_a(c) t2_aaaa(a,b,i,j)
    //               += -0.50 <k,j||c,d>_aaaa r2_1_aaa(c,d,k) t2_aaaa(a,b,i,j)
    //               += +1.00 <j,k||d,c>_abab r1_a(d) t2_aaaa(a,b,i,j) t1_1_bb(c,k)
    //               += +0.50 <j,k||c,d>_abab r2_1_abb(c,d,k) t2_aaaa(a,b,i,j)
    //               += +0.50 <j,k||d,c>_abab r2_1_abb(d,c,k) t2_aaaa(a,b,i,j)
    sigmar2_1_aaa("R,a,b,i") -= 0.25 * tmps_["48_aaa_Lvvo"]("R,a,b,i");
    tmps_["48_aaa_Lvvo"].~TArrayD();

    // flops: o0v1L1  = o1v2L1
    //  mems: o0v1L1  = o0v1L1
    tmps_["53_a_Lv"]("L,a")  = l2_1["bab"]("L,i,a,b") * t1_1["bb"]("b,i");

    // csigmal1_1_a += +1.00 t1_1_bb(b,i) l2_1_bab(i,a,b)
    csigmal1_1_a("L,a") += tmps_["53_a_Lv"]("L,a");

    // flops: o0v1L1  = o1v0L1 o1v1L1
    //  mems: o0v1L1  = o1v0L1 o0v1L1
    tmps_["54_a_Lv"]("L,a")  = (tmps_["19_a_Lo"]("L,i") + -2.00 * tmps_["20_a_Lo"]("L,i")) * dp["aa_ov"]("i,a");

    // sigmal1_a  = -0.50 d+_aa(j,a) t2_aaaa(c,b,i,j) l2_1_aaa(i,c,b)
    //           += +0.50 d+_aa(j,a) t2_abab(c,b,j,i) l2_1_bab(i,c,b)
    //           += +0.50 d+_aa(j,a) t2_abab(b,c,j,i) l2_1_bab(i,b,c)
    sigmal1_a("L,a")  = -0.50 * tmps_["54_a_Lv"]("L,a");

    // flops: o0v1L1  = o0v1L1
    //  mems: o0v1L1  = o0v1L1
    tmps_["51_a_Lv"]("L,a")  = -1.00 * scalars_["1"] * l1["a"]("L,a");

    // flops: o1v2L1  = o2v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["49_aaa_Lvvo"]("R,a,b,i")  = 4.00 * t2["aaaa"]("a,b,i,j") * tmps_["47_a_Lo"]("R,j");
    tmps_["49_aaa_Lvvo"].~TArrayD();

    // flops: o0v1L1  = o1v2L1
    //  mems: o0v1L1  = o0v1L1
    tmps_["50_a_Lv"]("L,a")  = l2_1["aaa"]("L,i,a,b") * t1_1["aa"]("b,i");

    // flops: o0v1L1  = o1v3L1 o1v2L1 o0v1L1 o1v2L1 o0v1L1 o1v3L1 o0v1L1 o0v1L1 o1v2L1 o0v1L1 o1v3L1 o0v1L1 o1v2L1 o0v1L1 o0v1L1 o0v1L1 o1v2L1 o0v1L1 o0v1L1 o1v3L1 o0v1L1 o1v3L1 o0v1L1 o1v3L1 o0v1L1 o0v2L1 o0v1L1 o0v2L1 o0v1L1 o1v3L1 o0v1L1 o0v1L1 o0v1L1 o1v2L1 o0v1L1 o1v2L1 o0v1L1 o1v3L1 o0v1L1 o0v2L1 o0v1L1 o1v2L1 o0v1L1 o1v3L1 o0v1L1 o0v2L1 o0v1L1 o0v2L1 o0v1L1 o1v2L1 o0v1L1 o0v2L1 o0v1L1 o0v1L1 o0v1L1 o1v2L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1
    //  mems: o0v1L1  = o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1
    tmps_["52_a_Lv"]("L,a")  = 0.25 * l2_1["aaa"]("L,k,e,b") * reused_["15_aaaa_vovv"]("a,k,e,b");
    tmps_["52_a_Lv"]("L,a") += f["aa_vo"]("b,k") * l2_1["aaa"]("L,k,b,a");
    tmps_["52_a_Lv"]("L,a") += scalars_["5"] * l1_1["a"]("L,a");
    tmps_["52_a_Lv"]("L,a") -= l2_1["bab"]("L,i,a,c") * reused_["37_bb_vo"]("c,i");
    tmps_["52_a_Lv"]("L,a") -= l2_1["bab"]("L,i,e,c") * reused_["27_aabb_vvvo"]("e,a,c,i");
    tmps_["52_a_Lv"]("L,a") -= dp["aa_vo"]("b,k") * l2["aaa"]("L,k,b,a");
    tmps_["52_a_Lv"]("L,a") -= l2_1["bab"]("L,i,e,c") * reused_["30_abab_vovv"]("a,i,e,c");
    tmps_["52_a_Lv"]("L,a") -= 2.00 * l2_1["bab"]("L,i,a,c") * reused_["35_bb_ov"]("i,c");
    tmps_["52_a_Lv"]("L,a") += 2.00 * scalars_["4"] * l1_1["a"]("L,a");
    tmps_["52_a_Lv"]("L,a") -= l2["aaa"]("L,k,b,a") * reused_["32_aa_vo"]("b,k");
    tmps_["52_a_Lv"]("L,a") -= l2_1["aaa"]("L,k,e,b") * reused_["26_aaaa_vvvo"]("e,a,b,k");
    tmps_["52_a_Lv"]("L,a") -= eri["abab_vvvo"]("e,c,a,i") * l2_1["bab"]("L,i,e,c");
    tmps_["52_a_Lv"]("L,a") -= l2_1["bab"]("L,i,b,d") * reused_["29_baab_vvvo"]("d,a,b,i");
    tmps_["52_a_Lv"]("L,a") += dp["aa_vv"]("b,a") * l1["a"]("L,b");
    tmps_["52_a_Lv"]("L,a") += l1_1["a"]("L,b") * reused_["3_aa_vv"]("b,a");
    tmps_["52_a_Lv"]("L,a") += 0.50 * eri["aaaa_vvvo"]("b,e,a,k") * l2_1["aaa"]("L,k,e,b");
    tmps_["52_a_Lv"]("L,a") -= w0 * l1_1["a"]("L,a");
    tmps_["52_a_Lv"]("L,a") -= l2_1["aaa"]("L,k,b,a") * reused_["36_aa_vo"]("b,k");
    tmps_["52_a_Lv"]("L,a") -= l2["bab"]("L,i,a,c") * reused_["33_bb_vo"]("c,i");
    tmps_["52_a_Lv"]("L,a") += l2_1["aaa"]("L,k,e,b") * reused_["25_aaaa_vovv"]("b,k,e,a");
    tmps_["52_a_Lv"]("L,a") += 0.50 * l1_1["a"]("L,b") * reused_["31_aa_vv"]("a,b");
    tmps_["52_a_Lv"]("L,a") -= f["bb_vo"]("c,i") * l2_1["bab"]("L,i,a,c");
    tmps_["52_a_Lv"]("L,a") += l2_1["bab"]("L,i,e,c") * reused_["28_aabb_vvvo"]("e,a,c,i");
    tmps_["52_a_Lv"]("L,a") += dp["aa_vv"]("e,a") * tmps_["50_a_Lv"]("L,e");
    tmps_["52_a_Lv"]("L,a") -= f["aa_vv"]("b,a") * l1_1["a"]("L,b");
    tmps_["52_a_Lv"]("L,a") += dp["bb_vo"]("c,i") * l2["bab"]("L,i,a,c");
    tmps_["52_a_Lv"]("L,a") += l1_1["a"]("L,b") * reused_["2_aa_vv"]("b,a");
    tmps_["52_a_Lv"]("L,a") += scalars_["2"] * l1["a"]("L,a");
    tmps_["52_a_Lv"]("L,a") -= 2.00 * l2_1["aaa"]("L,k,b,a") * reused_["34_aa_vo"]("b,k");
    tmps_["52_a_Lv"]("L,a") += 2.00 * scalars_["3"] * l1_1["a"]("L,a");
    tmps_["52_a_Lv"]("L,a") -= tmps_["51_a_Lv"]("L,a");
    tmps_["51_a_Lv"].~TArrayD();
    tmps_["50_a_Lv"].~TArrayD();

    // flops: o0v1L1  = o1v0L1 o1v1L1 o0v2L1 o0v1L1 o1v1L1 o0v1L1 o1v1L1 o0v1L1 o1v1L1 o0v1L1 o1v1L1 o0v1L1 o1v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1
    //  mems: o0v1L1  = o1v0L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1
    tmps_["55_a_Lv"]("L,a")  = (tmps_["38_a_Lo"]("L,i") + 0.50 * tmps_["39_a_Lo"]("L,i")) * dp["aa_ov"]("i,a");
    tmps_["55_a_Lv"]("L,a") += dp["aa_vv"]("b,a") * tmps_["53_a_Lv"]("L,b");
    tmps_["55_a_Lv"]("L,a") -= 2.00 * dp["aa_ov"]("k,a") * tmps_["42_a_Lo"]("L,k");
    tmps_["55_a_Lv"]("L,a") -= dp["aa_ov"]("i,a") * tmps_["41_a_Lo"]("L,i");
    tmps_["55_a_Lv"]("L,a") -= 2.00 * dp["aa_ov"]("i,a") * tmps_["33_a_Lo"]("L,i");
    tmps_["55_a_Lv"]("L,a") += f["aa_ov"]("i,a") * tmps_["20_a_Lo"]("L,i");
    tmps_["55_a_Lv"]("L,a") -= 0.50 * f["aa_ov"]("i,a") * tmps_["19_a_Lo"]("L,i");
    tmps_["55_a_Lv"]("L,a") += tmps_["52_a_Lv"]("L,a");
    tmps_["55_a_Lv"]("L,a") += 0.50 * t0_1 * tmps_["54_a_Lv"]("L,a");
    tmps_["52_a_Lv"].~TArrayD();

    // sigmal1_1_a  = -1.00 d-_aa(j,a) t2_1_aaaa(c,b,i,j) l2_1_aaa(i,c,b)
    //             += -0.50 d-_aa(j,a) t2_aaaa(c,b,i,j) l2_aaa(i,c,b)
    //             += -1.00 d-_aa(c,a) t1_1_bb(b,i) l2_1_bab(i,c,b)
    //             += +2.00 d-_aa(i,a) t1_1_aa(b,i) l1_1_a(b)
    //             += +0.50 d-_aa(j,a) t2_abab(c,b,j,i) l2_bab(i,c,b)
    //             += +0.50 d-_aa(j,a) t2_abab(b,c,j,i) l2_bab(i,b,c)
    //             += +1.00 d-_aa(j,a) t2_1_abab(c,b,j,i) l2_1_bab(i,c,b)
    //             += +1.00 d-_aa(j,a) t2_1_abab(b,c,j,i) l2_1_bab(i,b,c)
    //             += -0.50 f_aa(j,a) t2_abab(c,b,j,i) l2_1_bab(i,c,b)
    //             += -0.50 f_aa(j,a) t2_abab(b,c,j,i) l2_1_bab(i,b,c)
    //             += +0.50 f_aa(j,a) t2_aaaa(c,b,i,j) l2_1_aaa(i,c,b)
    //             += +1.00 d-_bb(j,c) t2_bbbb(c,b,i,j) t0_1 l2_1_bab(i,a,b)
    //             += -1.00 d-_aa(j,c) t2_abab(c,b,j,i) t0_1 l2_1_bab(i,a,b)
    //             += -1.00 d-_bb(b,i) t0_1 l2_1_bab(i,a,b)
    //             += -1.00 d-_aa(j,j) t1_1_bb(b,i) l2_1_bab(i,a,b)
    //             += -1.00 d-_bb(j,j) t1_1_bb(b,i) l2_1_bab(i,a,b)
    //             += -0.50 <k,j||c,i>_abab t2_abab(c,b,k,j) l2_1_bab(i,a,b)
    //             += -0.50 <j,k||c,i>_abab t2_abab(c,b,j,k) l2_1_bab(i,a,b)
    //             += -0.50 <j,b||c,d>_bbbb t2_bbbb(c,d,i,j) l2_1_bab(i,a,b)
    //             += +1.00 f_aa(j,c) t2_abab(c,b,j,i) l2_1_bab(i,a,b)
    //             += +0.50 <j,b||c,d>_abab t2_abab(c,d,j,i) l2_1_bab(i,a,b)
    //             += +0.50 <j,b||d,c>_abab t2_abab(d,c,j,i) l2_1_bab(i,a,b)
    //             += -0.50 <k,j||c,i>_bbbb t2_bbbb(c,b,k,j) l2_1_bab(i,a,b)
    //             += -1.00 f_bb(j,c) t2_bbbb(c,b,i,j) l2_1_bab(i,a,b)
    //             += -0.50 <j,i||j,i>_abab l1_1_a(a)
    //             += -0.50 <i,j||i,j>_abab l1_1_a(a)
    //             += +0.25 <j,i||b,c>_bbbb t2_bbbb(b,c,j,i) l1_1_a(a)
    //             += +1.00 f_bb(i,i) l1_1_a(a)
    //             += +0.25 <j,i||b,c>_aaaa t2_aaaa(b,c,j,i) l1_1_a(a)
    //             += -1.00 d-_aa(i,i) t0_1 l1_1_a(a)
    //             += +1.00 f_aa(i,i) l1_1_a(a)
    //             += -0.50 <j,i||j,i>_bbbb l1_1_a(a)
    //             += -0.50 <j,i||j,i>_aaaa l1_1_a(a)
    //             += +0.25 <j,i||b,c>_abab t2_abab(b,c,j,i) l1_1_a(a)
    //             += +0.25 <i,j||b,c>_abab t2_abab(b,c,i,j) l1_1_a(a)
    //             += +0.25 <j,i||c,b>_abab t2_abab(c,b,j,i) l1_1_a(a)
    //             += +0.25 <i,j||c,b>_abab t2_abab(c,b,i,j) l1_1_a(a)
    //             += -1.00 d-_bb(i,i) t0_1 l1_1_a(a)
    //             += +1.00 <j,c||d,a>_aaaa t2_abab(d,b,j,i) l2_1_bab(i,c,b)
    //             += -1.00 f_aa(b,i) l2_1_aaa(i,b,a)
    //             += +1.00 d-_aa(b,i) l2_aaa(i,b,a)
    //             += +0.25 <k,j||a,i>_abab t2_abab(c,b,k,j) l2_1_bab(i,c,b)
    //             += +0.25 <j,k||a,i>_abab t2_abab(c,b,j,k) l2_1_bab(i,c,b)
    //             += +0.25 <k,j||a,i>_abab t2_abab(b,c,k,j) l2_1_bab(i,b,c)
    //             += +0.25 <j,k||a,i>_abab t2_abab(b,c,j,k) l2_1_bab(i,b,c)
    //             += +2.00 d-_bb(j,i) t1_1_bb(b,j) l2_1_bab(i,a,b)
    //             += -2.00 d-_aa(j,c) t2_1_abab(c,b,j,i) l2_1_bab(i,a,b)
    //             += +2.00 d-_bb(j,c) t2_1_bbbb(c,b,i,j) l2_1_bab(i,a,b)
    //             += -2.00 d-_bb(b,c) t1_1_bb(c,i) l2_1_bab(i,a,b)
    //             += -2.00 d-_aa(i,b) t1_1_aa(b,i) l1_1_a(a)
    //             += +1.00 d-_bb(j,c) t2_abab(b,c,i,j) l2_aaa(i,b,a)
    //             += -1.00 d-_aa(j,c) t2_aaaa(c,b,i,j) l2_aaa(i,b,a)
    //             += +0.25 <k,j||a,i>_aaaa t2_aaaa(c,b,k,j) l2_1_aaa(i,c,b)
    //             += +1.00 <c,j||a,d>_abab t2_abab(b,d,i,j) l2_1_aaa(i,c,b)
    //             += +0.50 <c,b||a,i>_abab l2_1_bab(i,c,b)
    //             += +0.50 <b,c||a,i>_abab l2_1_bab(i,b,c)
    //             += -1.00 <j,c||a,d>_abab t2_abab(b,d,j,i) l2_1_bab(i,b,c)
    //             += -1.00 d-_aa(b,a) l1_a(b)
    //             += -1.00 d-_aa(b,a) t0_1 l1_1_a(b)
    //             += +0.50 <c,b||a,i>_aaaa l2_1_aaa(i,c,b)
    //             += +1.00 l1_1_a(a) w0
    //             += +1.00 d-_aa(b,i) t0_1 l2_1_aaa(i,b,a)
    //             += +0.50 <k,j||c,i>_aaaa t2_aaaa(c,b,k,j) l2_1_aaa(i,b,a)
    //             += +0.50 <k,j||i,c>_abab t2_abab(b,c,k,j) l2_1_aaa(i,b,a)
    //             += +0.50 <j,k||i,c>_abab t2_abab(b,c,j,k) l2_1_aaa(i,b,a)
    //             += +1.00 f_aa(j,c) t2_aaaa(c,b,i,j) l2_1_aaa(i,b,a)
    //             += -1.00 f_bb(j,c) t2_abab(b,c,i,j) l2_1_aaa(i,b,a)
    //             += -0.50 <b,j||c,d>_abab t2_abab(c,d,i,j) l2_1_aaa(i,b,a)
    //             += -0.50 <b,j||d,c>_abab t2_abab(d,c,i,j) l2_1_aaa(i,b,a)
    //             += +1.00 d-_bb(j,j) t1_1_aa(b,i) l2_1_aaa(i,b,a)
    //             += +1.00 d-_aa(j,j) t1_1_aa(b,i) l2_1_aaa(i,b,a)
    //             += +0.50 <j,b||c,d>_aaaa t2_aaaa(c,d,i,j) l2_1_aaa(i,b,a)
    //             += +1.00 d-_bb(j,c) t2_abab(b,c,i,j) t0_1 l2_1_aaa(i,b,a)
    //             += -1.00 d-_aa(j,c) t2_aaaa(c,b,i,j) t0_1 l2_1_aaa(i,b,a)
    //             += +1.00 d-_bb(j,c) t2_bbbb(c,b,i,j) l2_bab(i,a,b)
    //             += -1.00 d-_aa(j,c) t2_abab(c,b,j,i) l2_bab(i,a,b)
    //             += -1.00 <j,c||d,a>_aaaa t2_aaaa(d,b,i,j) l2_1_aaa(i,c,b)
    //             += -0.50 <j,i||c,a>_aaaa t2_aaaa(c,b,j,i) l1_1_a(b)
    //             += +1.00 f_bb(b,i) l2_1_bab(i,a,b)
    //             += -1.00 <c,j||a,d>_abab t2_bbbb(d,b,i,j) l2_1_bab(i,c,b)
    //             += -1.00 d-_aa(c,a) t1_1_aa(b,i) l2_1_aaa(i,c,b)
    //             += +1.00 f_aa(b,a) l1_1_a(b)
    //             += -1.00 d-_bb(b,i) l2_bab(i,a,b)
    //             += -0.50 <j,i||a,c>_abab t2_abab(b,c,j,i) l1_1_a(b)
    //             += -0.50 <i,j||a,c>_abab t2_abab(b,c,i,j) l1_1_a(b)
    //             += -1.00 d-_aa(i,i) l1_a(a)
    //             += +2.00 d-_aa(b,c) t1_1_aa(c,i) l2_1_aaa(i,b,a)
    //             += +2.00 d-_bb(j,c) t2_1_abab(b,c,i,j) l2_1_aaa(i,b,a)
    //             += -2.00 d-_aa(j,c) t2_1_aaaa(c,b,i,j) l2_1_aaa(i,b,a)
    //             += -2.00 d-_aa(j,i) t1_1_aa(b,j) l2_1_aaa(i,b,a)
    //             += -2.00 d-_bb(i,b) t1_1_bb(b,i) l1_1_a(a)
    //             += -1.00 d-_bb(i,i) l1_a(a)
    //             += -0.50 d-_aa(j,a) t2_aaaa(c,b,i,j) t0_1 l2_1_aaa(i,c,b)
    //             += +0.50 d-_aa(j,a) t2_abab(c,b,j,i) t0_1 l2_1_bab(i,c,b)
    //             += +0.50 d-_aa(j,a) t2_abab(b,c,j,i) t0_1 l2_1_bab(i,b,c)
    sigmal1_1_a("L,a")  = -1.00 * tmps_["55_a_Lv"]("L,a");
    tmps_["55_a_Lv"].~TArrayD();

    // flops: o0v1L1  = o0v1L1
    //  mems: o0v1L1  = o0v1L1
    tmps_["56_a_Lv"]("L,a")  = 0.50 * t0_1 * tmps_["54_a_Lv"]("L,a");
    tmps_["56_a_Lv"].~TArrayD();
    tmps_["54_a_Lv"].~TArrayD();

    // flops: o1v2L1  = o2v3L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["57_bba_Lvov"]("L,a,i,b")  = l2_1["bab"]("L,j,c,a") * eri["abab_vovo"]("c,i,b,j");

    // sigmal2_1_abb  = -1.00 <c,i||a,j>_abab l2_1_bab(j,c,b)
    sigmal2_1_abb("L,a,b,i")  = -1.00 * tmps_["57_bba_Lvov"]("L,b,i,a");

    // flops: o0v1L1  = o1v2L1
    //  mems: o0v1L1  = o0v1L1
    tmps_["67_a_Lv"]("L,a")  = l2["bab"]("L,i,a,b") * t1_1["bb"]("b,i");

    // csigmal1_a += +1.00 t1_1_bb(b,i) l2_bab(i,a,b)
    csigmal1_a("L,a") += tmps_["67_a_Lv"]("L,a");

    // flops: o1v2L1  = o2v3L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["66_abb_Lvvo"]("L,a,b,i")  = l2_1["aaa"]("L,j,c,a") * t2["abab"]("c,b,j,i");

    // flops: o2v1L1  = o2v2L1
    //  mems: o2v1L1  = o2v1L1
    tmps_["65_bab_Lovo"]("L,i,a,j")  = l2_1["bab"]("L,i,a,b") * t1_1["bb"]("b,j");

    // flops: o0v1L1  = o1v1L1
    //  mems: o0v1L1  = o0v1L1
    tmps_["64_a_Lv"]("L,a")  = 0.50 * f["aa_ov"]("i,a") * tmps_["39_a_Lo"]("L,i");

    // flops: o2v1L1  = o2v2L1
    //  mems: o2v1L1  = o2v1L1
    tmps_["62_abb_Loov"]("L,i,j,a")  = t1_1["aa"]("b,i") * l2_1["bab"]("L,j,b,a");

    // flops: o2v1L1  = o2v2L1
    //  mems: o2v1L1  = o2v1L1
    tmps_["61_aaa_Lovo"]("L,i,a,j")  = l2_1["aaa"]("L,i,b,a") * t1_1["aa"]("b,j");

    // flops: o1v0L1  = o2v1L1 o2v1L1 o1v0L1
    //  mems: o1v0L1  = o1v0L1 o1v0L1 o1v0L1
    tmps_["63_a_Lo"]("L,i")  = t1_1["aa"]("b,k") * tmps_["61_aaa_Lovo"]("L,k,b,i");
    tmps_["63_a_Lo"]("L,i") += t1_1["bb"]("a,j") * tmps_["62_abb_Loov"]("L,i,j,a");

    // flops: o1v0L1  = o2v2L1 o2v2L1 o1v0L1
    //  mems: o1v0L1  = o1v0L1 o1v0L1 o1v0L1
    tmps_["60_a_Lo"]("L,i")  = -0.50 * l2["aaa"]("L,k,a,c") * t2_1["aaaa"]("a,c,k,i");
    tmps_["60_a_Lo"]("L,i") += l2["bab"]("L,j,a,b") * t2_1["abab"]("a,b,i,j");

    // flops: o0v1L1  = o1v2L1
    //  mems: o0v1L1  = o0v1L1
    tmps_["59_a_Lv"]("L,a")  = l2["aaa"]("L,i,a,b") * t1_1["aa"]("b,i");

    // flops: o1v0L1  = o1v1L1
    //  mems: o1v0L1  = o1v0L1
    tmps_["58_a_Lo"]("L,i")  = l1["a"]("L,a") * t1_1["aa"]("a,i");

    // flops: o0v1L1  = o1v3L1 o1v3L1 o1v2L1 o1v2L1 o1v3L1 o1v3L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o1v2L1 o0v1L1 o0v1L1 o1v3L1 o0v1L1 o0v2L1 o0v1L1 o0v2L1 o0v1L1 o0v1L1 o0v1L1 o1v2L1 o0v1L1 o1v3L1 o0v1L1 o1v2L1 o0v1L1 o1v3L1 o0v1L1 o1v3L1 o0v1L1 o0v1L1 o0v1L1 o1v2L1 o0v1L1 o1v3L1 o0v1L1 o0v2L1 o0v1L1 o1v2L1 o0v1L1 o1v3L1 o0v1L1 o1v3L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o1v2L1 o0v1L1 o1v2L1 o0v1L1 o1v3L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o1v2L1 o0v1L1 o0v2L1 o0v1L1 o1v3L1 o0v1L1 o0v2L1 o0v1L1 o0v2L1 o0v1L1 o1v3L1 o0v1L1 o1v3L1 o0v1L1 o1v2L1 o0v1L1 o0v2L1 o0v1L1 o1v2L1 o0v1L1 o1v3L1 o0v1L1 o0v2L1 o0v1L1 o0v2L1 o0v1L1 o2v2L1 o0v1L1 o1v1L1 o0v1L1 o2v2L1 o0v1L1 o1v1L1 o0v1L1 o1v1L1 o0v1L1 o2v2L1 o0v1L1 o1v1L1 o0v1L1 o1v1L1 o0v1L1 o1v2L1 o0v1L1 o1v1L1 o0v1L1 o1v1L1 o0v1L1 o2v2L1 o0v1L1 o2v2L1 o0v1L1 o1v1L1 o0v1L1 o1v1L1 o0v1L1 o1v2L1 o0v1L1 o1v1L1 o0v1L1 o1v1L1 o0v1L1 o1v1L1 o0v1L1 o2v2L1 o0v1L1 o2v2L1 o0v1L1 o1v1L1 o0v1L1 o1v1L1 o0v1L1 o1v1L1 o0v1L1 o0v1L1 o1v1L1 o0v1L1
    //  mems: o0v1L1  = o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1
    tmps_["68_a_Lv"]("L,a")  = -0.25 * l2["aaa"]("L,k,d,c") * reused_["15_aaaa_vovv"]("a,k,d,c");
    tmps_["68_a_Lv"]("L,a") -= 0.50 * l2_1["aaa"]("L,k,d,c") * reused_["41_aaaa_vvvo"]("c,d,a,k");
    tmps_["68_a_Lv"]("L,a") += f["bb_vo"]("b,j") * l2["bab"]("L,j,a,b");
    tmps_["68_a_Lv"]("L,a") += l2["aaa"]("L,k,c,a") * reused_["36_aa_vo"]("c,k");
    tmps_["68_a_Lv"]("L,a") += l2["bab"]("L,j,d,b") * reused_["30_abab_vovv"]("a,j,d,b");
    tmps_["68_a_Lv"]("L,a") += l2_1["bab"]("L,j,d,b") * reused_["38_aabb_vvvo"]("d,a,b,j");
    tmps_["68_a_Lv"]("L,a") += l2["bab"]("L,j,a,b") * reused_["35_bb_ov"]("j,b");
    tmps_["68_a_Lv"]("L,a") -= 0.25 * l2_1["aaa"]("L,k,d,c") * reused_["23_aaaa_vovv"]("a,k,d,c");
    tmps_["68_a_Lv"]("L,a") -= 0.50 * l1["a"]("L,c") * reused_["31_aa_vv"]("a,c");
    tmps_["68_a_Lv"]("L,a") -= l1["a"]("L,c") * reused_["2_aa_vv"]("c,a");
    tmps_["68_a_Lv"]("L,a") -= scalars_["1"] * l1_1["a"]("L,a");
    tmps_["68_a_Lv"]("L,a") += l2_1["aaa"]("L,k,c,a") * reused_["32_aa_vo"]("c,k");
    tmps_["68_a_Lv"]("L,a") += l2["aaa"]("L,k,d,c") * reused_["26_aaaa_vvvo"]("d,a,c,k");
    tmps_["68_a_Lv"]("L,a") += l2["bab"]("L,j,a,b") * reused_["37_bb_vo"]("b,j");
    tmps_["68_a_Lv"]("L,a") += l2["bab"]("L,j,d,b") * reused_["27_aabb_vvvo"]("d,a,b,j");
    tmps_["68_a_Lv"]("L,a") += l2["bab"]("L,j,c,g") * reused_["29_baab_vvvo"]("g,a,c,j");
    tmps_["68_a_Lv"]("L,a") -= scalars_["3"] * l1["a"]("L,a");
    tmps_["68_a_Lv"]("L,a") -= dp["bb_vo"]("b,j") * l2_1["bab"]("L,j,a,b");
    tmps_["68_a_Lv"]("L,a") += eri["abab_vvvo"]("d,b,a,j") * l2["bab"]("L,j,d,b");
    tmps_["68_a_Lv"]("L,a") -= dp["aa_vv"]("c,a") * l1_1["a"]("L,c");
    tmps_["68_a_Lv"]("L,a") += l2["aaa"]("L,k,c,a") * reused_["34_aa_vo"]("c,k");
    tmps_["68_a_Lv"]("L,a") -= l2["bab"]("L,j,d,b") * reused_["28_aabb_vvvo"]("d,a,b,j");
    tmps_["68_a_Lv"]("L,a") += l2_1["aaa"]("L,k,d,c") * reused_["42_aaaa_vvvo"]("d,a,c,k");
    tmps_["68_a_Lv"]("L,a") -= scalars_["5"] * l1["a"]("L,a");
    tmps_["68_a_Lv"]("L,a") += scalars_["6"] * l1_1["a"]("L,a");
    tmps_["68_a_Lv"]("L,a") -= f["aa_vo"]("c,k") * l2["aaa"]("L,k,c,a");
    tmps_["68_a_Lv"]("L,a") += 0.50 * l2_1["bab"]("L,j,a,b") * reused_["63_bb_ov"]("j,b");
    tmps_["68_a_Lv"]("L,a") += l2_1["bab"]("L,j,d,b") * reused_["46_abab_vovv"]("a,j,d,b");
    tmps_["68_a_Lv"]("L,a") -= scalars_["2"] * l1_1["a"]("L,a");
    tmps_["68_a_Lv"]("L,a") -= scalars_["4"] * l1["a"]("L,a");
    tmps_["68_a_Lv"]("L,a") += dp["aa_vo"]("c,k") * l2_1["aaa"]("L,k,c,a");
    tmps_["68_a_Lv"]("L,a") -= l1_1["a"]("L,c") * reused_["13_aa_vv"]("c,a");
    tmps_["68_a_Lv"]("L,a") -= 0.50 * eri["aaaa_vvvo"]("c,d,a,k") * l2["aaa"]("L,k,d,c");
    tmps_["68_a_Lv"]("L,a") += f["aa_vv"]("c,a") * l1["a"]("L,c");
    tmps_["68_a_Lv"]("L,a") -= l1["a"]("L,c") * reused_["3_aa_vv"]("c,a");
    tmps_["68_a_Lv"]("L,a") += l2_1["bab"]("L,j,c,g") * reused_["43_abba_vovv"]("c,j,g,a");
    tmps_["68_a_Lv"]("L,a") += l2_1["bab"]("L,j,d,b") * reused_["60_abab_vovv"]("a,j,d,b");
    tmps_["68_a_Lv"]("L,a") += l2_1["bab"]("L,j,a,b") * reused_["33_bb_vo"]("b,j");
    tmps_["68_a_Lv"]("L,a") -= dp["aa_vv"]("d,a") * tmps_["59_a_Lv"]("L,d");
    tmps_["68_a_Lv"]("L,a") += l2_1["aaa"]("L,k,c,a") * reused_["62_aa_vo"]("c,k");
    tmps_["68_a_Lv"]("L,a") -= l2["aaa"]("L,k,d,c") * reused_["25_aaaa_vovv"]("c,k,d,a");
    tmps_["68_a_Lv"]("L,a") += l1_1["a"]("L,c") * reused_["61_aa_vv"]("c,a");
    tmps_["68_a_Lv"]("L,a") -= dp["aa_vv"]("d,a") * tmps_["67_a_Lv"]("L,d");
    tmps_["68_a_Lv"]("L,a") += reused_["50_bbaa_voov"]("g,j,m,a") * tmps_["62_abb_Loov"]("L,m,j,g");
    tmps_["68_a_Lv"]("L,a") += dp["aa_ov"]("i,a") * tmps_["60_a_Lo"]("L,i");
    tmps_["68_a_Lv"]("L,a") += reused_["39_aaaa_voov"]("d,k,m,a") * tmps_["40_aaa_Lovo"]("L,k,d,m");
    tmps_["68_a_Lv"]("L,a") -= 0.50 * reused_["4_aa_ov"]("i,a") * tmps_["39_a_Lo"]("L,i");
    tmps_["68_a_Lv"]("L,a") -= reused_["11_aa_ov"]("i,a") * tmps_["20_a_Lo"]("L,i");
    tmps_["68_a_Lv"]("L,a") += eri["baab_vovo"]("g,i,a,j") * tmps_["62_abb_Loov"]("L,i,j,g");
    tmps_["68_a_Lv"]("L,a") += reused_["4_aa_ov"]("i,a") * tmps_["41_a_Lo"]("L,i");
    tmps_["68_a_Lv"]("L,a") += reused_["4_aa_ov"]("k,a") * tmps_["42_a_Lo"]("L,k");
    tmps_["68_a_Lv"]("L,a") -= t1_1["bb"]("b,n") * tmps_["57_bba_Lvov"]("L,b,n,a");
    tmps_["68_a_Lv"]("L,a") += 0.50 * reused_["11_aa_ov"]("i,a") * tmps_["19_a_Lo"]("L,i");
    tmps_["68_a_Lv"]("L,a") -= 0.50 * reused_["4_aa_ov"]("i,a") * tmps_["38_a_Lo"]("L,i");
    tmps_["68_a_Lv"]("L,a") += reused_["40_bbaa_voov"]("g,j,m,a") * tmps_["62_abb_Loov"]("L,m,j,g");
    tmps_["68_a_Lv"]("L,a") += reused_["49_aaaa_voov"]("d,k,m,a") * tmps_["40_aaa_Lovo"]("L,k,d,m");
    tmps_["68_a_Lv"]("L,a") += dp["aa_ov"]("k,a") * tmps_["58_a_Lo"]("L,k");
    tmps_["68_a_Lv"]("L,a") -= f["aa_ov"]("i,a") * tmps_["33_a_Lo"]("L,i");
    tmps_["68_a_Lv"]("L,a") -= reused_["12_bb_ov"]("n,g") * tmps_["66_abb_Lvvo"]("L,a,g,n");
    tmps_["68_a_Lv"]("L,a") += 0.50 * reused_["54_aa_ov"]("i,a") * tmps_["19_a_Lo"]("L,i");
    tmps_["68_a_Lv"]("L,a") += reused_["4_aa_ov"]("i,a") * tmps_["33_a_Lo"]("L,i");
    tmps_["68_a_Lv"]("L,a") -= f["aa_ov"]("i,a") * tmps_["41_a_Lo"]("L,i");
    tmps_["68_a_Lv"]("L,a") -= eri["aaaa_vovo"]("d,i,a,k") * tmps_["40_aaa_Lovo"]("L,k,d,i");
    tmps_["68_a_Lv"]("L,a") += reused_["45_abba_voov"]("d,j,l,a") * tmps_["65_bab_Lovo"]("L,j,d,l");
    tmps_["68_a_Lv"]("L,a") += 0.50 * f["aa_ov"]("i,a") * tmps_["38_a_Lo"]("L,i");
    tmps_["68_a_Lv"]("L,a") -= reused_["54_aa_ov"]("i,a") * tmps_["20_a_Lo"]("L,i");
    tmps_["68_a_Lv"]("L,a") -= f["aa_ov"]("k,a") * tmps_["42_a_Lo"]("L,k");
    tmps_["68_a_Lv"]("L,a") += tmps_["64_a_Lv"]("L,a");
    tmps_["68_a_Lv"]("L,a") += dp["aa_ov"]("i,a") * tmps_["63_a_Lo"]("L,i");
    tmps_["64_a_Lv"].~TArrayD();
    tmps_["63_a_Lo"].~TArrayD();
    tmps_["60_a_Lo"].~TArrayD();
    tmps_["59_a_Lv"].~TArrayD();
    tmps_["58_a_Lo"].~TArrayD();
    tmps_["57_bba_Lvov"].~TArrayD();
    tmps_["40_aaa_Lovo"].~TArrayD();

    // sigmal1_a += +1.00 d-_aa(j,a) t1_1_bb(b,i) t1_1_aa(c,j) l2_1_bab(i,c,b)
    //           += +1.00 d-_aa(j,a) t1_1_aa(b,i) t1_1_aa(c,j) l2_1_aaa(i,c,b)
    //           += +1.00 <j,c||d,a>_aaaa t2_1_abab(d,b,j,i) l2_1_bab(i,c,b)
    //           += +0.25 <k,j||a,i>_abab t2_abab(c,b,k,j) l2_bab(i,c,b)
    //           += +0.25 <j,k||a,i>_abab t2_abab(c,b,j,k) l2_bab(i,c,b)
    //           += +0.25 <k,j||a,i>_abab t2_abab(b,c,k,j) l2_bab(i,b,c)
    //           += +0.25 <j,k||a,i>_abab t2_abab(b,c,j,k) l2_bab(i,b,c)
    //           += +1.00 d-_aa(b,i) t0_1 l2_aaa(i,b,a)
    //           += +0.50 <k,j||c,i>_aaaa t2_aaaa(c,b,k,j) l2_aaa(i,b,a)
    //           += +0.50 <k,j||i,c>_abab t2_abab(b,c,k,j) l2_aaa(i,b,a)
    //           += +0.50 <j,k||i,c>_abab t2_abab(b,c,j,k) l2_aaa(i,b,a)
    //           += +1.00 f_aa(j,c) t2_aaaa(c,b,i,j) l2_aaa(i,b,a)
    //           += -1.00 f_bb(j,c) t2_abab(b,c,i,j) l2_aaa(i,b,a)
    //           += -0.50 <b,j||c,d>_abab t2_abab(c,d,i,j) l2_aaa(i,b,a)
    //           += -0.50 <b,j||d,c>_abab t2_abab(d,c,i,j) l2_aaa(i,b,a)
    //           += +1.00 d-_bb(j,j) t1_1_aa(b,i) l2_aaa(i,b,a)
    //           += +1.00 d-_aa(j,j) t1_1_aa(b,i) l2_aaa(i,b,a)
    //           += +0.50 <j,b||c,d>_aaaa t2_aaaa(c,d,i,j) l2_aaa(i,b,a)
    //           += +1.00 d-_bb(j,c) t2_abab(b,c,i,j) t0_1 l2_aaa(i,b,a)
    //           += -1.00 d-_aa(j,c) t2_aaaa(c,b,i,j) t0_1 l2_aaa(i,b,a)
    //           += +1.00 f_bb(b,i) l2_bab(i,a,b)
    //           += -0.50 <c,b||d,a>_aaaa t1_1_aa(d,i) l2_1_aaa(i,c,b)
    //           += -1.00 <j,c||d,a>_aaaa t2_1_aaaa(d,b,i,j) l2_1_aaa(i,c,b)
    //           += -0.25 <k,j||d,a>_aaaa t2_aaaa(c,b,k,j) t1_1_aa(d,i) l2_1_aaa(i,c,b)
    //           += +1.00 d-_bb(j,i) t1_1_bb(b,j) l2_bab(i,a,b)
    //           += -1.00 d-_aa(j,c) t2_1_abab(c,b,j,i) l2_bab(i,a,b)
    //           += +1.00 d-_bb(j,c) t2_1_bbbb(c,b,i,j) l2_bab(i,a,b)
    //           += -1.00 d-_bb(b,c) t1_1_bb(c,i) l2_bab(i,a,b)
    //           += +0.25 <k,j||a,i>_aaaa t2_aaaa(c,b,k,j) l2_aaa(i,c,b)
    //           += +0.25 <k,j||a,i>_aaaa t2_1_aaaa(c,b,k,j) l2_1_aaa(i,c,b)
    //           += -0.50 <j,i||c,a>_aaaa t2_aaaa(c,b,j,i) l1_a(b)
    //           += -0.50 <j,i||a,c>_abab t2_abab(b,c,j,i) l1_a(b)
    //           += -0.50 <i,j||a,c>_abab t2_abab(b,c,i,j) l1_a(b)
    //           += -1.00 d+_bb(i,i) l1_1_a(a)
    //           += +1.00 d+_bb(j,c) t2_abab(b,c,i,j) l2_1_aaa(i,b,a)
    //           += -1.00 d+_aa(j,c) t2_aaaa(c,b,i,j) l2_1_aaa(i,b,a)
    //           += +1.00 <c,j||a,d>_abab t2_abab(b,d,i,j) l2_aaa(i,c,b)
    //           += +1.00 d-_bb(j,c) t2_bbbb(c,b,i,j) t0_1 l2_bab(i,a,b)
    //           += -1.00 d-_aa(j,c) t2_abab(c,b,j,i) t0_1 l2_bab(i,a,b)
    //           += -1.00 d-_bb(b,i) t0_1 l2_bab(i,a,b)
    //           += -1.00 d-_aa(j,j) t1_1_bb(b,i) l2_bab(i,a,b)
    //           += -1.00 d-_bb(j,j) t1_1_bb(b,i) l2_bab(i,a,b)
    //           += -0.50 <k,j||c,i>_abab t2_abab(c,b,k,j) l2_bab(i,a,b)
    //           += -0.50 <j,k||c,i>_abab t2_abab(c,b,j,k) l2_bab(i,a,b)
    //           += -0.50 <j,b||c,d>_bbbb t2_bbbb(c,d,i,j) l2_bab(i,a,b)
    //           += +1.00 f_aa(j,c) t2_abab(c,b,j,i) l2_bab(i,a,b)
    //           += +0.50 <j,b||c,d>_abab t2_abab(c,d,j,i) l2_bab(i,a,b)
    //           += +0.50 <j,b||d,c>_abab t2_abab(d,c,j,i) l2_bab(i,a,b)
    //           += -0.50 <k,j||c,i>_bbbb t2_bbbb(c,b,k,j) l2_bab(i,a,b)
    //           += -1.00 f_bb(j,c) t2_bbbb(c,b,i,j) l2_bab(i,a,b)
    //           += +1.00 <j,c||d,a>_aaaa t2_abab(d,b,j,i) l2_bab(i,c,b)
    //           += -1.00 <j,c||a,d>_abab t2_abab(b,d,j,i) l2_bab(i,b,c)
    //           += -1.00 d-_bb(i,b) t1_1_bb(b,i) l1_a(a)
    //           += -1.00 d+_bb(b,i) l2_1_bab(i,a,b)
    //           += +0.50 <c,b||a,i>_abab l2_bab(i,c,b)
    //           += +0.50 <b,c||a,i>_abab l2_bab(i,b,c)
    //           += -1.00 d+_aa(b,a) l1_1_a(b)
    //           += +1.00 d-_aa(b,c) t1_1_aa(c,i) l2_aaa(i,b,a)
    //           += +1.00 d-_bb(j,c) t2_1_abab(b,c,i,j) l2_aaa(i,b,a)
    //           += -1.00 d-_aa(j,c) t2_1_aaaa(c,b,i,j) l2_aaa(i,b,a)
    //           += -1.00 d-_aa(j,i) t1_1_aa(b,j) l2_aaa(i,b,a)
    //           += -1.00 <c,j||a,d>_abab t2_bbbb(d,b,i,j) l2_bab(i,c,b)
    //           += +1.00 <c,j||a,d>_abab t2_1_abab(b,d,i,j) l2_1_aaa(i,c,b)
    //           += -0.50 <j,i||j,i>_abab l1_a(a)
    //           += -0.50 <i,j||i,j>_abab l1_a(a)
    //           += +0.25 <j,i||b,c>_bbbb t2_bbbb(b,c,j,i) l1_a(a)
    //           += +1.00 f_bb(i,i) l1_a(a)
    //           += +0.25 <j,i||b,c>_aaaa t2_aaaa(b,c,j,i) l1_a(a)
    //           += -1.00 d-_aa(i,i) t0_1 l1_a(a)
    //           += +1.00 f_aa(i,i) l1_a(a)
    //           += -0.50 <j,i||j,i>_bbbb l1_a(a)
    //           += -0.50 <j,i||j,i>_aaaa l1_a(a)
    //           += +0.25 <j,i||b,c>_abab t2_abab(b,c,j,i) l1_a(a)
    //           += +0.25 <i,j||b,c>_abab t2_abab(b,c,i,j) l1_a(a)
    //           += +0.25 <j,i||c,b>_abab t2_abab(c,b,j,i) l1_a(a)
    //           += +0.25 <i,j||c,b>_abab t2_abab(c,b,i,j) l1_a(a)
    //           += -1.00 d-_bb(i,i) t0_1 l1_a(a)
    //           += +1.00 t0_1 l1_1_a(a) w0
    //           += +0.25 <j,i||b,c>_abab t2_1_abab(b,c,j,i) l1_1_a(a)
    //           += +0.25 <i,j||b,c>_abab t2_1_abab(b,c,i,j) l1_1_a(a)
    //           += +0.25 <j,i||c,b>_abab t2_1_abab(c,b,j,i) l1_1_a(a)
    //           += +0.25 <i,j||c,b>_abab t2_1_abab(c,b,i,j) l1_1_a(a)
    //           += -1.00 d-_bb(i,b) t0_1 t1_1_bb(b,i) l1_1_a(a)
    //           += +0.25 <j,i||b,c>_bbbb t2_1_bbbb(b,c,j,i) l1_1_a(a)
    //           += -1.00 d-_aa(i,b) t0_1 t1_1_aa(b,i) l1_1_a(a)
    //           += +1.00 f_bb(i,b) t1_1_bb(b,i) l1_1_a(a)
    //           += +0.25 <j,i||b,c>_aaaa t2_1_aaaa(b,c,j,i) l1_1_a(a)
    //           += +1.00 f_aa(i,b) t1_1_aa(b,i) l1_1_a(a)
    //           += -1.00 f_aa(b,i) l2_aaa(i,b,a)
    //           += -0.50 <k,j||c,d>_bbbb t2_bbbb(c,d,i,j) t1_1_bb(b,k) l2_1_bab(i,a,b)
    //           += +1.00 t1_1_bb(b,i) l2_1_bab(i,a,b) w0
    //           += -0.50 <k,j||c,i>_abab t2_1_abab(c,b,k,j) l2_1_bab(i,a,b)
    //           += -0.50 <j,k||c,i>_abab t2_1_abab(c,b,j,k) l2_1_bab(i,a,b)
    //           += -1.00 d-_bb(j,c) t1_1_bb(b,i) t1_1_bb(c,j) l2_1_bab(i,a,b)
    //           += +1.00 f_aa(j,c) t2_1_abab(c,b,j,i) l2_1_bab(i,a,b)
    //           += -1.00 f_bb(j,i) t1_1_bb(b,j) l2_1_bab(i,a,b)
    //           += +1.00 <j,b||c,i>_bbbb t1_1_bb(c,j) l2_1_bab(i,a,b)
    //           += -1.00 d-_aa(j,c) t1_1_bb(b,i) t1_1_aa(c,j) l2_1_bab(i,a,b)
    //           += -1.00 f_bb(j,c) t2_1_bbbb(c,b,i,j) l2_1_bab(i,a,b)
    //           += -0.50 <j,b||c,d>_bbbb t2_1_bbbb(c,d,i,j) l2_1_bab(i,a,b)
    //           += +1.00 f_bb(b,c) t1_1_bb(c,i) l2_1_bab(i,a,b)
    //           += -0.50 <k,j||c,i>_bbbb t2_1_bbbb(c,b,k,j) l2_1_bab(i,a,b)
    //           += +0.50 <j,b||c,d>_abab t2_1_abab(c,d,j,i) l2_1_bab(i,a,b)
    //           += +0.50 <j,b||d,c>_abab t2_1_abab(d,c,j,i) l2_1_bab(i,a,b)
    //           += +1.00 <j,b||c,i>_abab t1_1_aa(c,j) l2_1_bab(i,a,b)
    //           += -0.50 <k,j||c,d>_bbbb t2_bbbb(c,b,k,j) t1_1_bb(d,i) l2_1_bab(i,a,b)
    //           += -1.00 <k,j||c,d>_aaaa t2_abab(c,b,j,i) t1_1_aa(d,k) l2_1_bab(i,a,b)
    //           += +1.00 <k,j||c,d>_bbbb t2_bbbb(c,b,i,j) t1_1_bb(d,k) l2_1_bab(i,a,b)
    //           += +2.00 d-_bb(j,c) t1_1_bb(b,j) t1_1_bb(c,i) l2_1_bab(i,a,b)
    //           += -1.00 <k,j||d,c>_abab t2_bbbb(c,b,i,j) t1_1_aa(d,k) l2_1_bab(i,a,b)
    //           += -0.50 <j,k||c,d>_abab t2_abab(c,d,j,i) t1_1_bb(b,k) l2_1_bab(i,a,b)
    //           += -0.50 <j,k||d,c>_abab t2_abab(d,c,j,i) t1_1_bb(b,k) l2_1_bab(i,a,b)
    //           += +1.00 <j,k||c,d>_abab t2_abab(c,b,j,i) t1_1_bb(d,k) l2_1_bab(i,a,b)
    //           += -0.50 <k,j||c,d>_abab t2_abab(c,b,k,j) t1_1_bb(d,i) l2_1_bab(i,a,b)
    //           += -0.50 <j,k||c,d>_abab t2_abab(c,b,j,k) t1_1_bb(d,i) l2_1_bab(i,a,b)
    //           += +1.00 d-_bb(j,i) t0_1 t1_1_bb(b,j) l2_1_bab(i,a,b)
    //           += -1.00 d-_aa(j,c) t0_1 t2_1_abab(c,b,j,i) l2_1_bab(i,a,b)
    //           += +1.00 d-_bb(j,c) t0_1 t2_1_bbbb(c,b,i,j) l2_1_bab(i,a,b)
    //           += -1.00 d-_bb(b,c) t0_1 t1_1_bb(c,i) l2_1_bab(i,a,b)
    //           += +0.25 <k,j||a,i>_abab t2_1_abab(c,b,k,j) l2_1_bab(i,c,b)
    //           += +0.25 <j,k||a,i>_abab t2_1_abab(c,b,j,k) l2_1_bab(i,c,b)
    //           += +0.25 <k,j||a,i>_abab t2_1_abab(b,c,k,j) l2_1_bab(i,b,c)
    //           += +0.25 <j,k||a,i>_abab t2_1_abab(b,c,j,k) l2_1_bab(i,b,c)
    //           += -1.00 d+_aa(i,i) l1_1_a(a)
    //           += -1.00 d-_aa(i,b) t1_1_aa(b,i) l1_a(a)
    //           += +1.00 d+_aa(b,i) l2_1_aaa(i,b,a)
    //           += -0.50 <j,i||a,c>_abab t2_1_abab(b,c,j,i) l1_1_a(b)
    //           += -0.50 <i,j||a,c>_abab t2_1_abab(b,c,i,j) l1_1_a(b)
    //           += +1.00 <b,i||a,c>_abab t1_1_bb(c,i) l1_1_a(b)
    //           += +0.50 <c,b||a,i>_aaaa l2_aaa(i,c,b)
    //           += +1.00 f_aa(b,a) l1_a(b)
    //           += -1.00 d-_aa(b,a) t0_1 l1_a(b)
    //           += -1.00 <j,c||a,d>_abab t2_1_abab(b,d,j,i) l2_1_bab(i,b,c)
    //           += +0.25 <k,j||a,d>_abab t2_abab(c,b,k,j) t1_1_bb(d,i) l2_1_bab(i,c,b)
    //           += +0.25 <j,k||a,d>_abab t2_abab(c,b,j,k) t1_1_bb(d,i) l2_1_bab(i,c,b)
    //           += +0.25 <k,j||a,d>_abab t2_abab(b,c,k,j) t1_1_bb(d,i) l2_1_bab(i,b,c)
    //           += +0.25 <j,k||a,d>_abab t2_abab(b,c,j,k) t1_1_bb(d,i) l2_1_bab(i,b,c)
    //           += -1.00 <c,j||a,d>_abab t2_1_bbbb(d,b,i,j) l2_1_bab(i,c,b)
    //           += +0.50 <c,b||a,d>_abab t1_1_bb(d,i) l2_1_bab(i,c,b)
    //           += +0.50 <b,c||a,d>_abab t1_1_bb(d,i) l2_1_bab(i,b,c)
    //           += +1.00 d+_bb(j,c) t2_bbbb(c,b,i,j) l2_1_bab(i,a,b)
    //           += -1.00 d+_aa(j,c) t2_abab(c,b,j,i) l2_1_bab(i,a,b)
    //           += -1.00 d-_aa(c,a) t1_1_aa(b,i) l2_aaa(i,c,b)
    //           += +0.50 <k,j||c,d>_abab t2_abab(c,d,i,j) t1_1_aa(b,k) l2_1_aaa(i,b,a)
    //           += +0.50 <k,j||d,c>_abab t2_abab(d,c,i,j) t1_1_aa(b,k) l2_1_aaa(i,b,a)
    //           += +1.00 <k,j||c,d>_bbbb t2_abab(b,c,i,j) t1_1_bb(d,k) l2_1_aaa(i,b,a)
    //           += +0.50 <k,j||c,d>_aaaa t2_aaaa(c,d,i,j) t1_1_aa(b,k) l2_1_aaa(i,b,a)
    //           += +0.50 <k,j||c,d>_aaaa t2_aaaa(c,b,k,j) t1_1_aa(d,i) l2_1_aaa(i,b,a)
    //           += +1.00 <j,k||c,d>_abab t2_aaaa(c,b,i,j) t1_1_bb(d,k) l2_1_aaa(i,b,a)
    //           += +0.50 <k,j||d,c>_abab t2_abab(b,c,k,j) t1_1_aa(d,i) l2_1_aaa(i,b,a)
    //           += +0.50 <j,k||d,c>_abab t2_abab(b,c,j,k) t1_1_aa(d,i) l2_1_aaa(i,b,a)
    //           += -1.00 <k,j||c,d>_aaaa t2_aaaa(c,b,i,j) t1_1_aa(d,k) l2_1_aaa(i,b,a)
    //           += -2.00 d-_aa(j,c) t1_1_aa(b,j) t1_1_aa(c,i) l2_1_aaa(i,b,a)
    //           += +1.00 d-_aa(b,c) t0_1 t1_1_aa(c,i) l2_1_aaa(i,b,a)
    //           += +1.00 d-_bb(j,c) t0_1 t2_1_abab(b,c,i,j) l2_1_aaa(i,b,a)
    //           += -1.00 d-_aa(j,c) t0_1 t2_1_aaaa(c,b,i,j) l2_1_aaa(i,b,a)
    //           += -1.00 d-_aa(j,i) t0_1 t1_1_aa(b,j) l2_1_aaa(i,b,a)
    //           += +1.00 f_aa(j,i) t1_1_aa(b,j) l2_1_aaa(i,b,a)
    //           += -1.00 f_aa(b,c) t1_1_aa(c,i) l2_1_aaa(i,b,a)
    //           += +1.00 f_aa(j,c) t2_1_aaaa(c,b,i,j) l2_1_aaa(i,b,a)
    //           += -1.00 <b,j||i,c>_abab t1_1_bb(c,j) l2_1_aaa(i,b,a)
    //           += -0.50 <b,j||c,d>_abab t2_1_abab(c,d,i,j) l2_1_aaa(i,b,a)
    //           += -0.50 <b,j||d,c>_abab t2_1_abab(d,c,i,j) l2_1_aaa(i,b,a)
    //           += +0.50 <k,j||c,i>_aaaa t2_1_aaaa(c,b,k,j) l2_1_aaa(i,b,a)
    //           += -1.00 <j,b||c,i>_aaaa t1_1_aa(c,j) l2_1_aaa(i,b,a)
    //           += -1.00 f_bb(j,c) t2_1_abab(b,c,i,j) l2_1_aaa(i,b,a)
    //           += +0.50 <k,j||i,c>_abab t2_1_abab(b,c,k,j) l2_1_aaa(i,b,a)
    //           += +0.50 <j,k||i,c>_abab t2_1_abab(b,c,j,k) l2_1_aaa(i,b,a)
    //           += +0.50 <j,b||c,d>_aaaa t2_1_aaaa(c,d,i,j) l2_1_aaa(i,b,a)
    //           += -1.00 t1_1_aa(b,i) l2_1_aaa(i,b,a) w0
    //           += +1.00 d-_aa(j,c) t1_1_aa(b,i) t1_1_aa(c,j) l2_1_aaa(i,b,a)
    //           += +1.00 d-_bb(j,c) t1_1_aa(b,i) t1_1_bb(c,j) l2_1_aaa(i,b,a)
    //           += -1.00 <j,c||d,a>_aaaa t2_aaaa(d,b,i,j) l2_aaa(i,c,b)
    //           += +1.00 <i,b||c,a>_aaaa t1_1_aa(c,i) l1_1_a(b)
    //           += -0.50 <j,i||c,a>_aaaa t2_1_aaaa(c,b,j,i) l1_1_a(b)
    //           += -1.00 d-_aa(c,a) t1_1_bb(b,i) l2_bab(i,c,b)
    //           += +1.00 <k,j||a,d>_abab t2_bbbb(d,c,i,j) t1_1_aa(b,k) l2_1_bab(i,b,c)
    //           += +0.50 d-_aa(j,a) t2_1_abab(c,b,j,i) l2_bab(i,c,b)
    //           += +0.50 d-_aa(j,a) t2_1_abab(b,c,j,i) l2_bab(i,b,c)
    //           += -0.50 d-_aa(j,a) t2_1_aaaa(c,b,i,j) l2_aaa(i,c,b)
    //           += +1.00 <k,j||d,a>_aaaa t2_aaaa(d,c,i,j) t1_1_aa(b,k) l2_1_aaa(i,c,b)
    //           += -0.50 d-_aa(j,a) t2_aaaa(c,b,i,j) t0_1 l2_aaa(i,c,b)
    //           += -0.50 <j,k||a,d>_abab t2_abab(c,b,j,i) t1_1_bb(d,k) l2_1_bab(i,c,b)
    //           += -0.50 <j,k||a,d>_abab t2_abab(b,c,j,i) t1_1_bb(d,k) l2_1_bab(i,b,c)
    //           += -1.00 <j,c||a,i>_abab t1_1_aa(b,j) l2_1_bab(i,b,c)
    //           += +0.50 d-_aa(j,a) t2_abab(c,b,j,i) t0_1 l2_bab(i,c,b)
    //           += +0.50 d-_aa(j,a) t2_abab(b,c,j,i) t0_1 l2_bab(i,b,c)
    //           += +1.00 d-_aa(i,a) t0_1 t1_1_aa(b,i) l1_1_a(b)
    //           += -1.00 <c,j||a,i>_abab t1_1_bb(b,j) l2_1_bab(i,c,b)
    //           += +0.50 <j,k||a,d>_abab t2_aaaa(c,b,i,j) t1_1_bb(d,k) l2_1_aaa(i,c,b)
    //           += -0.50 d-_aa(j,a) t0_1 t2_1_aaaa(c,b,i,j) l2_1_aaa(i,c,b)
    //           += +1.00 <k,j||d,a>_aaaa t2_abab(d,c,j,i) t1_1_aa(b,k) l2_1_bab(i,b,c)
    //           += +1.00 <k,j||a,d>_abab t2_abab(c,d,i,j) t1_1_aa(b,k) l2_1_aaa(i,c,b)
    //           += +1.00 d-_aa(i,a) t1_1_aa(b,i) l1_a(b)
    //           += -0.50 f_aa(j,a) t2_1_abab(c,b,j,i) l2_1_bab(i,c,b)
    //           += -0.50 f_aa(j,a) t2_1_abab(b,c,j,i) l2_1_bab(i,b,c)
    //           += -1.00 <k,j||d,c>_abab t2_abab(b,c,i,j) t1_1_aa(d,k) l2_1_aaa(i,b,a)
    //           += +0.50 <k,j||d,a>_aaaa t2_aaaa(c,b,i,j) t1_1_aa(d,k) l2_1_aaa(i,c,b)
    //           += +0.50 d-_aa(j,a) t0_1 t2_1_abab(c,b,j,i) l2_1_bab(i,c,b)
    //           += +0.50 d-_aa(j,a) t0_1 t2_1_abab(b,c,j,i) l2_1_bab(i,b,c)
    //           += -0.50 f_aa(j,a) t2_abab(c,b,j,i) l2_bab(i,c,b)
    //           += -0.50 f_aa(j,a) t2_abab(b,c,j,i) l2_bab(i,b,c)
    //           += +1.00 <j,c||a,i>_aaaa t1_1_aa(b,j) l2_1_aaa(i,c,b)
    //           += +1.00 <j,k||a,d>_abab t2_abab(c,d,j,i) t1_1_bb(b,k) l2_1_bab(i,c,b)
    //           += +0.50 f_aa(j,a) t2_1_aaaa(c,b,i,j) l2_1_aaa(i,c,b)
    //           += -0.50 <k,j||d,a>_aaaa t2_abab(c,b,j,i) t1_1_aa(d,k) l2_1_bab(i,c,b)
    //           += -0.50 <k,j||d,a>_aaaa t2_abab(b,c,j,i) t1_1_aa(d,k) l2_1_bab(i,b,c)
    //           += -1.00 f_aa(i,a) t1_1_aa(b,i) l1_1_a(b)
    //           += +0.50 f_aa(j,a) t2_aaaa(c,b,i,j) l2_aaa(i,c,b)
    sigmal1_a("L,a") += tmps_["68_a_Lv"]("L,a");
    tmps_["68_a_Lv"].~TArrayD();

    // flops: o0v1L1  = o0v2L1
    //  mems: o0v1L1  = o0v1L1
    tmps_["73_a_Lv"]("R,a")  = dp["aa_vv"]("a,b") * r1["a"]("R,b");

    // sigmar1_1_a += -1.00 d+_aa(a,b) r1_a(b)
    sigmar1_1_a("R,a") -= tmps_["73_a_Lv"]("R,a");

    // flops: o0v1L1  = o1v2L1 o1v2L1 o0v1L1
    //  mems: o0v1L1  = o0v1L1 o0v1L1 o0v1L1
    tmps_["74_a_Lv"]("R,a")  = -1.00 * dp["bb_ov"]("j,c") * r2["abb"]("R,a,c,j");
    tmps_["74_a_Lv"]("R,a") += dp["aa_ov"]("i,b") * r2["aaa"]("R,b,a,i");

    // sigmar1_1_a += +1.00 d+_aa(i,b) r2_aaa(b,a,i)
    //             += -1.00 d+_bb(i,b) r2_abb(a,b,i)
    sigmar1_1_a("R,a") += tmps_["74_a_Lv"]("R,a");

    // flops: o1v2L1  = o2v3L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["72_abb_Lvov"]("R,a,i,b")  = r2["aaa"]("R,c,a,j") * eri["abab_oovv"]("j,i,c,b");

    // flops: o1v2L1  = o1v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["71_aaa_Lvvo"]("R,a,b,i")  = -1.00 * r1_1["a"]("R,a") * reused_["32_aa_vo"]("b,i");

    // flops: o0v1L1  = o1v2L1
    //  mems: o0v1L1  = o0v1L1
    tmps_["69_a_Lv"]("R,a")  = -1.00 * dp["bb_ov"]("i,b") * r2["abb"]("R,a,b,i");
    tmps_["69_a_Lv"].~TArrayD();

    // flops: o2v1L1  = o2v2L1
    //  mems: o2v1L1  = o2v1L1
    tmps_["70_aaa_Lvoo"]("R,a,i,j")  = r2["aaa"]("R,b,a,i") * dp["aa_ov"]("j,b");
}