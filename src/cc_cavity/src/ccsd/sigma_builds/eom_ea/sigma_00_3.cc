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

#include "cc_cavity/include/ccsd/eom_ea_ccsd.h"

void hilbert::EOM_EA_CCSD::sigma_ea_00_3() {

    // get reference to electronic integrals
    TArrayMap &V_blks_ = cc_wfn_->V_blks_;
    TArrayMap &F_blks_ = cc_wfn_->F_blks_;
    TArrayMap &Id_blks_ = cc_wfn_->Id_blks_;
    TArrayMap &eri = cc_wfn_->V_blks_;
    TArrayMap &f = cc_wfn_->F_blks_;
    TArrayMap &Id = cc_wfn_->Id_blks_;

    TArrayMap t2 = {
            {"aaaa", cc_wfn_->amplitudes_["t2_aaaa"]},
            {"abab", cc_wfn_->amplitudes_["t2_abab"]},
            {"bbbb", cc_wfn_->amplitudes_["t2_bbbb"]}
    };

    /// reference right operators
    TArrayMap r1 = {
            {"a", evec_blks_["r1_a"]},
            {"b", evec_blks_["r1_b"]}
    };

    TArrayMap r2 = {
            {"aaa", evec_blks_["r2_aaa"]},
            {"abb", evec_blks_["r2_abb"]},
            {"aba", evec_blks_["r2_aba"]},
            {"bbb", evec_blks_["r2_bbb"]}
    };

    /// reference left operators

    TArrayMap l1 = {
            {"a", evec_blks_["l1_a"]},
            {"b", evec_blks_["l1_b"]}
    };
    TArrayMap l2 = {
            {"aaa", evec_blks_["l2_aaa"]},
            {"bab", evec_blks_["l2_bab"]},
            {"aab", evec_blks_["l2_aab"]},
            {"bbb", evec_blks_["l2_bbb"]}
    };

    /// reference sigma vectors
    auto &sigmar1_a = sigvec_blks_["sigmar1_a"];
    auto &sigmar1_b = sigvec_blks_["sigmar1_b"];
    auto &sigmal1_a = sigvec_blks_["sigmal1_a"];
    auto &sigmal1_b = sigvec_blks_["sigmal1_b"];

    auto &sigmar2_aaa = sigvec_blks_["sigmar2_aaa"];
    auto &sigmar2_abb = sigvec_blks_["sigmar2_abb"];
    auto &sigmal2_aaa = sigvec_blks_["sigmal2_aaa"];
    auto &sigmal2_abb = sigvec_blks_["sigmal2_abb"];


    // sigmal1_b  = +0.50 f_bb(j,a) t2_abab(c,b,i,j) l2_aab(i,c,b)
    //           += +0.50 f_bb(j,a) t2_abab(b,c,i,j) l2_aab(i,b,c)
    //           += +0.50 f_bb(j,a) t2_bbbb(c,b,i,j) l2_bbb(i,c,b)
    //           += -1.00 f_aa(b,i) l2_aab(i,b,a)
    //           += -0.50 <j,i||c,a>_bbbb t2_bbbb(c,b,j,i) l1_b(b)
    //           += +0.50 <k,j||c,i>_aaaa t2_aaaa(c,b,k,j) l2_aab(i,b,a)
    //           += -0.50 <b,j||c,d>_abab t2_abab(c,d,i,j) l2_aab(i,b,a)
    //           += -0.50 <b,j||d,c>_abab t2_abab(d,c,i,j) l2_aab(i,b,a)
    //           += +0.50 <j,b||c,d>_aaaa t2_aaaa(c,d,i,j) l2_aab(i,b,a)
    //           += +1.00 f_aa(j,c) t2_aaaa(c,b,i,j) l2_aab(i,b,a)
    //           += +0.50 <k,j||i,c>_abab t2_abab(b,c,k,j) l2_aab(i,b,a)
    //           += +0.50 <j,k||i,c>_abab t2_abab(b,c,j,k) l2_aab(i,b,a)
    //           += -1.00 f_bb(j,c) t2_abab(b,c,i,j) l2_aab(i,b,a)
    //           += +1.00 f_bb(b,a) l1_b(b)
    //           += -0.50 <j,i||j,i>_bbbb l1_b(a)
    //           += -0.50 <j,i||j,i>_aaaa l1_b(a)
    //           += -0.50 <j,i||j,i>_abab l1_b(a)
    //           += -0.50 <i,j||i,j>_abab l1_b(a)
    //           += +0.25 <j,i||b,c>_bbbb t2_bbbb(b,c,j,i) l1_b(a)
    //           += +0.25 <j,i||b,c>_abab t2_abab(b,c,j,i) l1_b(a)
    //           += +0.25 <i,j||b,c>_abab t2_abab(b,c,i,j) l1_b(a)
    //           += +0.25 <j,i||c,b>_abab t2_abab(c,b,j,i) l1_b(a)
    //           += +0.25 <i,j||c,b>_abab t2_abab(c,b,i,j) l1_b(a)
    //           += +1.00 f_bb(i,i) l1_b(a)
    //           += +0.25 <j,i||b,c>_aaaa t2_aaaa(b,c,j,i) l1_b(a)
    //           += +1.00 f_aa(i,i) l1_b(a)
    //           += -1.00 f_bb(b,i) l2_bbb(i,b,a)
    //           += -0.50 <c,b||i,a>_abab l2_aab(i,c,b)
    //           += -0.50 <b,c||i,a>_abab l2_aab(i,b,c)
    //           += +1.00 <j,c||d,a>_abab t2_aaaa(d,b,i,j) l2_aab(i,b,c)
    //           += -1.00 <j,c||d,a>_bbbb t2_abab(b,d,i,j) l2_aab(i,b,c)
    //           += +1.00 <j,c||d,a>_abab t2_abab(d,b,j,i) l2_bbb(i,c,b)
    //           += +0.25 <k,j||a,i>_bbbb t2_bbbb(c,b,k,j) l2_bbb(i,c,b)
    //           += -1.00 <j,c||d,a>_bbbb t2_bbbb(d,b,i,j) l2_bbb(i,c,b)
    //           += -0.50 <j,i||c,a>_abab t2_abab(c,b,j,i) l1_b(b)
    //           += -0.50 <i,j||c,a>_abab t2_abab(c,b,i,j) l1_b(b)
    //           += +0.50 <c,b||a,i>_bbbb l2_bbb(i,c,b)
    //           += +1.00 <c,j||d,a>_abab t2_abab(d,b,i,j) l2_aab(i,c,b)
    //           += -0.25 <k,j||i,a>_abab t2_abab(c,b,k,j) l2_aab(i,c,b)
    //           += -0.25 <j,k||i,a>_abab t2_abab(c,b,j,k) l2_aab(i,c,b)
    //           += -0.25 <k,j||i,a>_abab t2_abab(b,c,k,j) l2_aab(i,b,c)
    //           += -0.25 <j,k||i,a>_abab t2_abab(b,c,j,k) l2_aab(i,b,c)
    //           += -1.00 f_aa(j,c) t2_abab(c,b,j,i) l2_bbb(i,b,a)
    //           += +0.50 <j,b||c,d>_bbbb t2_bbbb(c,d,i,j) l2_bbb(i,b,a)
    //           += -0.50 <j,b||c,d>_abab t2_abab(c,d,j,i) l2_bbb(i,b,a)
    //           += -0.50 <j,b||d,c>_abab t2_abab(d,c,j,i) l2_bbb(i,b,a)
    //           += +1.00 f_bb(j,c) t2_bbbb(c,b,i,j) l2_bbb(i,b,a)
    //           += +0.50 <k,j||c,i>_bbbb t2_bbbb(c,b,k,j) l2_bbb(i,b,a)
    //           += +0.50 <k,j||c,i>_abab t2_abab(c,b,k,j) l2_bbb(i,b,a)
    //           += +0.50 <j,k||c,i>_abab t2_abab(c,b,j,k) l2_bbb(i,b,a)
    sigmal1_b("L,a")  = tmps_["16_b_Lv"]("L,a");
    tmps_["16_b_Lv"].~TArrayD();

    // flops: o1v2L1  = o1v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["19_aaa_Lvov"]("R,a,i,b")  = -1.00 * f["aa_vo"]("a,i") * r1["a"]("R,b");

    // flops: o0v1L1  = o1v3L1 o0v2L1 o0v2L1 o1v2L1 o0v1L1 o1v2L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o1v2L1 o0v1L1 o0v1L1 o1v3L1 o0v1L1 o1v3L1 o0v1L1 o0v2L1 o0v1L1 o1v3L1 o0v1L1 o1v3L1 o0v1L1 o1v2L1 o0v1L1
    //  mems: o0v1L1  = o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1
    tmps_["17_b_Lv"]("L,a")  = eri["abba_vvvo"]("c,b,a,j") * l2["aab"]("L,j,c,b");
    tmps_["17_b_Lv"]("L,a") += f["bb_vv"]("b,a") * l1["b"]("L,b");
    tmps_["17_b_Lv"]("L,a") -= 0.50 * l1["b"]("L,b") * reused_["22_bb_vv"]("a,b");
    tmps_["17_b_Lv"]("L,a") -= f["aa_vo"]("f,j") * l2["aab"]("L,j,f,a");
    tmps_["17_b_Lv"]("L,a") -= 0.50 * l2["aab"]("L,j,f,a") * reused_["13_aa_ov"]("j,f");
    tmps_["17_b_Lv"]("L,a") -= 0.50 * scalars_["1"] * l1["b"]("L,a");
    tmps_["17_b_Lv"]("L,a") -= f["bb_vo"]("b,i") * l2["bbb"]("L,i,b,a");
    tmps_["17_b_Lv"]("L,a") -= l2["aab"]("L,j,f,d") * reused_["24_aabb_vovv"]("f,j,d,a");
    tmps_["17_b_Lv"]("L,a") -= l2["bbb"]("L,i,d,b") * reused_["25_bbbb_vovv"]("b,i,d,a");
    tmps_["17_b_Lv"]("L,a") -= l1["b"]("L,b") * reused_["1_bb_vv"]("a,b");
    tmps_["17_b_Lv"]("L,a") -= 0.50 * eri["bbbb_vvvo"]("b,d,a,i") * l2["bbb"]("L,i,d,b");
    tmps_["17_b_Lv"]("L,a") += l2["aab"]("L,j,c,b") * reused_["23_abba_vvvo"]("c,a,b,j");
    tmps_["17_b_Lv"]("L,a") -= l2["bbb"]("L,i,b,a") * reused_["14_bb_vo"]("b,i");
    tmps_["17_b_Lv"].~TArrayD();

    // flops: o1v2L1  = o2v3L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["18_abb_Lvov"]("R,a,i,b")  = eri["abab_oovv"]("j,i,c,b") * r2["aaa"]("R,c,a,j");

    // flops: o1v2L1  = o2v3L1 o1v3L1 o1v2L1 o2v3L1 o1v3L1 o1v3L1 o2v3L1 o1v3L1 o1v3L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["20_aaa_Lvov"]("R,a,i,b")  = eri["abba_vovo"]("a,j,c,i") * r2["abb"]("R,b,c,j");
    tmps_["20_aaa_Lvov"]("R,a,i,b") += r1["a"]("R,d") * reused_["26_aaaa_vvvo"]("a,d,b,i");
    tmps_["20_aaa_Lvov"]("R,a,i,b") -= 0.50 * r1["a"]("R,b") * reused_["13_aa_ov"]("i,a");
    tmps_["20_aaa_Lvov"]("R,a,i,b") -= eri["aaaa_vovo"]("a,l,e,i") * r2["aaa"]("R,e,b,l");
    tmps_["20_aaa_Lvov"]("R,a,i,b") -= reused_["4_aa_vv"]("a,d") * r2["aaa"]("R,d,b,i");
    tmps_["20_aaa_Lvov"]("R,a,i,b") += f["aa_vv"]("a,e") * r2["aaa"]("R,e,b,i");
    tmps_["20_aaa_Lvov"]("R,a,i,b") -= reused_["27_aaaa_voov"]("a,i,k,d") * r2["aaa"]("R,d,b,k");
    tmps_["20_aaa_Lvov"]("R,a,i,b") += 0.50 * reused_["3_aa_vv"]("a,d") * r2["aaa"]("R,d,b,i");
    tmps_["20_aaa_Lvov"]("R,a,i,b") += r1["a"]("R,d") * reused_["6_aaaa_vvvo"]("a,d,b,i");
    tmps_["20_aaa_Lvov"]("R,a,i,b") += reused_["29_aabb_voov"]("a,i,m,f") * r2["abb"]("R,b,f,m");
    tmps_["20_aaa_Lvov"]("R,a,i,b") -= reused_["28_aabb_voov"]("a,i,m,f") * r2["abb"]("R,b,f,m");
    tmps_["20_aaa_Lvov"]("R,a,i,b") += tmps_["19_aaa_Lvov"]("R,a,i,b");
    tmps_["20_aaa_Lvov"]("R,a,i,b") += t2["abab"]("a,c,i,j") * tmps_["18_abb_Lvov"]("R,b,j,c");
    tmps_["19_aaa_Lvov"].~TArrayD();
    tmps_["18_abb_Lvov"].~TArrayD();

    // sigmar2_aaa  = +1.00 P(a,b) <k,j||d,c>_abab r2_aaa(d,b,k) t2_abab(a,c,i,j)
    //             += +1.00 P(a,b) <a,j||d,c>_abab r1_a(d) t2_abab(b,c,i,j)
    //             += -0.50 P(a,b) <k,j||c,d>_aaaa r2_aaa(d,b,i) t2_aaaa(c,a,k,j)
    //             += +1.00 P(a,b) <k,j||c,d>_aaaa r2_aaa(d,b,k) t2_aaaa(c,a,i,j)
    //             += +1.00 P(a,b) <j,k||c,d>_abab r2_abb(b,d,k) t2_aaaa(c,a,i,j)
    //             += +1.00 P(a,b) f_aa(a,c) r2_aaa(c,b,i)
    //             += -0.50 P(a,b) <k,j||d,c>_abab r2_aaa(d,b,i) t2_abab(a,c,k,j)
    //             += -0.50 P(a,b) <j,k||d,c>_abab r2_aaa(d,b,i) t2_abab(a,c,j,k)
    //             += +1.00 P(a,b) <j,a||c,i>_aaaa r2_aaa(c,b,j)
    //             += +1.00 P(a,b) <k,j||c,d>_bbbb r2_abb(b,d,k) t2_abab(a,c,i,j)
    //             += +0.50 P(a,b) <k,j||c,i>_aaaa r1_a(b) t2_aaaa(c,a,k,j)
    //             += -0.50 P(a,b) <a,j||c,d>_abab r1_a(b) t2_abab(c,d,i,j)
    //             += -0.50 P(a,b) <a,j||d,c>_abab r1_a(b) t2_abab(d,c,i,j)
    //             += +0.50 P(a,b) <j,a||c,d>_aaaa r1_a(b) t2_aaaa(c,d,i,j)
    //             += +1.00 P(a,b) f_aa(j,c) r1_a(b) t2_aaaa(c,a,i,j)
    //             += +0.50 P(a,b) <k,j||i,c>_abab r1_a(b) t2_abab(a,c,k,j)
    //             += +0.50 P(a,b) <j,k||i,c>_abab r1_a(b) t2_abab(a,c,j,k)
    //             += -1.00 P(a,b) f_bb(j,c) r1_a(b) t2_abab(a,c,i,j)
    //             += -1.00 P(a,b) <j,a||c,d>_aaaa r1_a(d) t2_aaaa(c,b,i,j)
    //             += -1.00 P(a,b) <a,j||i,c>_abab r2_abb(b,c,j)
    //             += -1.00 P(a,b) f_aa(a,i) r1_a(b)
    sigmar2_aaa("R,a,b,i")  = tmps_["20_aaa_Lvov"]("R,a,i,b");
    sigmar2_aaa("R,a,b,i") -= tmps_["20_aaa_Lvov"]("R,b,i,a");
    tmps_["20_aaa_Lvov"].~TArrayD();

    // flops: o3v0L1  = o3v2L1
    //  mems: o3v0L1  = o3v0L1
    tmps_["21_aaa_Looo"]("L,i,j,k")  = l2["aaa"]("L,i,a,b") * t2["aaaa"]("a,b,j,k");

    // flops: o1v2L1  = o3v2L1 o2v2L1 o2v2L1 o1v4L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["22_aaa_Lovv"]("L,i,a,b")  = 0.50 * eri["aaaa_oovv"]("k,j,a,b") * tmps_["21_aaa_Looo"]("L,i,j,k");
    tmps_["22_aaa_Lovv"]("L,i,a,b") += 2.00 * l2["aaa"]("L,k,a,b") * f["aa_oo"]("i,k");
    tmps_["22_aaa_Lovv"]("L,i,a,b") += 2.00 * l2["aaa"]("L,k,a,b") * reused_["31_aa_oo"]("k,i");
    tmps_["22_aaa_Lovv"]("L,i,a,b") += eri["aaaa_vvvv"]("d,c,a,b") * l2["aaa"]("L,i,c,d");
    tmps_["22_aaa_Lovv"]("L,i,a,b") -= 2.00 * eri["aaaa_vovv"]("d,i,a,b") * l1["a"]("L,d");
    tmps_["22_aaa_Lovv"]("L,i,a,b") += scalars_["1"] * l2["aaa"]("L,i,a,b");
    tmps_["22_aaa_Lovv"]("L,i,a,b") += l2["aaa"]("L,k,a,b") * reused_["30_aa_oo"]("i,k");
    tmps_["22_aaa_Lovv"]("L,i,a,b") += eri["aaaa_oovv"]("i,j,a,b") * tmps_["5_a_Lo"]("L,j");

    // sigmal2_aaa  = -0.50 <i,k||a,b>_aaaa t2_aaaa(d,c,j,k) l2_aaa(j,d,c)
    //             += +0.50 <i,k||a,b>_aaaa t2_abab(d,c,k,j) l2_bab(j,d,c)
    //             += +0.50 <i,k||a,b>_aaaa t2_abab(c,d,k,j) l2_bab(j,c,d)
    //             += +0.50 <d,c||a,b>_aaaa l2_aaa(i,d,c)
    //             += -0.50 <i,k||c,d>_abab t2_abab(c,d,j,k) l2_aaa(j,a,b)
    //             += -0.50 <i,k||d,c>_abab t2_abab(d,c,j,k) l2_aaa(j,a,b)
    //             += -1.00 <i,c||a,b>_aaaa l1_a(c)
    //             += -0.50 <k,j||k,j>_bbbb l2_aaa(i,a,b)
    //             += -0.50 <k,j||k,j>_aaaa l2_aaa(i,a,b)
    //             += -0.50 <k,j||k,j>_abab l2_aaa(i,a,b)
    //             += -0.50 <j,k||j,k>_abab l2_aaa(i,a,b)
    //             += +0.25 <k,j||c,d>_bbbb t2_bbbb(c,d,k,j) l2_aaa(i,a,b)
    //             += +0.25 <k,j||c,d>_abab t2_abab(c,d,k,j) l2_aaa(i,a,b)
    //             += +0.25 <j,k||c,d>_abab t2_abab(c,d,j,k) l2_aaa(i,a,b)
    //             += +0.25 <k,j||d,c>_abab t2_abab(d,c,k,j) l2_aaa(i,a,b)
    //             += +0.25 <j,k||d,c>_abab t2_abab(d,c,j,k) l2_aaa(i,a,b)
    //             += +1.00 f_bb(j,j) l2_aaa(i,a,b)
    //             += +0.25 <k,j||c,d>_aaaa t2_aaaa(c,d,k,j) l2_aaa(i,a,b)
    //             += +1.00 f_aa(j,j) l2_aaa(i,a,b)
    //             += -1.00 f_aa(i,j) l2_aaa(j,a,b)
    //             += -0.50 <i,k||c,d>_aaaa t2_aaaa(c,d,j,k) l2_aaa(j,a,b)
    //             += +0.25 <k,j||a,b>_aaaa t2_aaaa(d,c,k,j) l2_aaa(i,d,c)
    sigmal2_aaa("L,a,b,i")  = -0.50 * tmps_["22_aaa_Lovv"]("L,i,a,b");
    tmps_["22_aaa_Lovv"].~TArrayD();

    // flops: o1v2L1  = o1v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["25_aaa_Lvvo"]("R,a,b,i")  = -0.50 * scalars_["1"] * r2["aaa"]("R,a,b,i");
}