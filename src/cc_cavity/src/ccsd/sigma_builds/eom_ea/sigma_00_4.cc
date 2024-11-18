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

void hilbert::EOM_EA_CCSD::sigma_ea_00_4() {

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


    // flops: o1v2L1  = o3v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["23_aaa_Lvvo"]("L,a,b,i")  = 0.50 * eri["aaaa_oovv"]("j,k,a,b") * tmps_["21_aaa_Looo"]("L,i,k,j");
    tmps_["23_aaa_Lvvo"].~TArrayD();
    tmps_["21_aaa_Looo"].~TArrayD();

    // flops: o3v0L1  = o3v2L1
    //  mems: o3v0L1  = o3v0L1
    tmps_["24_aaa_Looo"]("R,i,j,k")  = eri["aaaa_oovv"]("j,k,a,b") * r2["aaa"]("R,a,b,i");

    // flops: o1v2L1  = o3v2L1 o2v2L1 o1v2L1 o2v2L1 o1v4L1 o2v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["26_aaa_Lvvo"]("R,a,b,i")  = -0.25 * t2["aaaa"]("a,b,k,j") * tmps_["24_aaa_Looo"]("R,i,j,k");
    tmps_["26_aaa_Lvvo"]("R,a,b,i") += t2["aaaa"]("a,b,i,j") * tmps_["12_a_Lo"]("R,j");
    tmps_["26_aaa_Lvvo"]("R,a,b,i") -= r2["aaa"]("R,a,b,j") * f["aa_oo"]("j,i");
    tmps_["26_aaa_Lvvo"]("R,a,b,i") += 0.50 * eri["aaaa_vvvv"]("a,b,c,d") * r2["aaa"]("R,c,d,i");
    tmps_["26_aaa_Lvvo"]("R,a,b,i") += 0.50 * r2["aaa"]("R,a,b,k") * reused_["32_aa_oo"]("i,k");
    tmps_["26_aaa_Lvvo"]("R,a,b,i") += eri["aaaa_vvvo"]("a,b,c,i") * r1["a"]("R,c");
    tmps_["26_aaa_Lvvo"]("R,a,b,i") -= r2["aaa"]("R,a,b,k") * reused_["31_aa_oo"]("i,k");
    tmps_["26_aaa_Lvvo"]("R,a,b,i") -= 0.50 * r1["a"]("R,c") * reused_["10_aaaa_vovv"]("c,i,a,b");
    tmps_["26_aaa_Lvvo"]("R,a,b,i") += tmps_["25_aaa_Lvvo"]("R,a,b,i");
    tmps_["25_aaa_Lvvo"].~TArrayD();
    tmps_["24_aaa_Looo"].~TArrayD();
    tmps_["12_a_Lo"].~TArrayD();

    // sigmar2_aaa += +1.00 f_aa(j,c) r1_a(c) t2_aaaa(a,b,i,j)
    //             += +0.50 <j,k||c,d>_abab r2_abb(c,d,k) t2_aaaa(a,b,i,j)
    //             += +0.50 <j,k||d,c>_abab r2_abb(d,c,k) t2_aaaa(a,b,i,j)
    //             += -0.50 <k,j||c,d>_aaaa r2_aaa(c,d,k) t2_aaaa(a,b,i,j)
    //             += +0.25 <k,j||c,d>_aaaa r2_aaa(c,d,i) t2_aaaa(a,b,k,j)
    //             += +1.00 <a,b||c,i>_aaaa r1_a(c)
    //             += -0.50 <k,j||c,d>_aaaa r2_aaa(a,b,k) t2_aaaa(c,d,i,j)
    //             += +0.50 <a,b||c,d>_aaaa r2_aaa(c,d,i)
    //             += -1.00 f_aa(j,i) r2_aaa(a,b,j)
    //             += -0.50 <k,j||c,d>_abab r2_aaa(a,b,k) t2_abab(c,d,i,j)
    //             += -0.50 <k,j||d,c>_abab r2_aaa(a,b,k) t2_abab(d,c,i,j)
    //             += +0.50 <k,j||c,i>_aaaa r1_a(c) t2_aaaa(a,b,k,j)
    //             += -0.50 <k,j||k,j>_bbbb r2_aaa(a,b,i)
    //             += -0.50 <k,j||k,j>_aaaa r2_aaa(a,b,i)
    //             += -0.50 <k,j||k,j>_abab r2_aaa(a,b,i)
    //             += -0.50 <j,k||j,k>_abab r2_aaa(a,b,i)
    //             += +0.25 <k,j||c,d>_bbbb r2_aaa(a,b,i) t2_bbbb(c,d,k,j)
    //             += +0.25 <k,j||c,d>_abab r2_aaa(a,b,i) t2_abab(c,d,k,j)
    //             += +0.25 <j,k||c,d>_abab r2_aaa(a,b,i) t2_abab(c,d,j,k)
    //             += +0.25 <k,j||d,c>_abab r2_aaa(a,b,i) t2_abab(d,c,k,j)
    //             += +0.25 <j,k||d,c>_abab r2_aaa(a,b,i) t2_abab(d,c,j,k)
    //             += +1.00 f_bb(j,j) r2_aaa(a,b,i)
    //             += +0.25 <k,j||c,d>_aaaa r2_aaa(a,b,i) t2_aaaa(c,d,k,j)
    //             += +1.00 f_aa(j,j) r2_aaa(a,b,i)
    sigmar2_aaa("R,a,b,i") += tmps_["26_aaa_Lvvo"]("R,a,b,i");
    tmps_["26_aaa_Lvvo"].~TArrayD();

    // flops: o3v0L1  = o3v2L1
    //  mems: o3v0L1  = o3v0L1
    tmps_["27_abb_Looo"]("L,i,j,k")  = l2["bab"]("L,k,a,b") * t2["abab"]("a,b,i,j");

    // flops: o1v2L1  = o3v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["30_abb_Lvvo"]("L,a,b,i")  = 2.00 * eri["abab_oovv"]("j,k,a,b") * tmps_["27_abb_Looo"]("L,j,k,i");
    tmps_["27_abb_Looo"].~TArrayD();

    // flops: o1v2L1  = o1v3L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["29_abb_Lvov"]("L,a,i,b")  = -2.00 * f["aa_vv"]("c,a") * l2["bab"]("L,i,c,b");

    // flops: o1v2L1  = o2v3L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["28_bba_Lvvo"]("L,a,b,i")  = l2["bab"]("L,j,c,a") * t2["abab"]("c,b,i,j");

    // flops: o1v2L1  = o2v2L1 o2v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v4L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o1v3L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["31_bab_Lovv"]("L,i,a,b")  = 0.50 * eri["abab_oovv"]("j,i,a,b") * tmps_["5_a_Lo"]("L,j");
    tmps_["31_bab_Lovv"]("L,i,a,b") -= 0.50 * l2["bab"]("L,k,a,b") * reused_["35_bb_oo"]("i,k");
    tmps_["31_bab_Lovv"]("L,i,a,b") += l2["bab"]("L,k,a,f") * reused_["19_bbbb_voov"]("f,k,i,b");
    tmps_["31_bab_Lovv"]("L,i,a,b") -= 0.50 * l2["bab"]("L,i,a,f") * reused_["22_bb_vv"]("b,f");
    tmps_["31_bab_Lovv"]("L,i,a,b") -= l2["bab"]("L,i,a,f") * reused_["1_bb_vv"]("b,f");
    tmps_["31_bab_Lovv"]("L,i,a,b") -= l2["bab"]("L,i,d,b") * reused_["4_aa_vv"]("d,a");
    tmps_["31_bab_Lovv"]("L,i,a,b") += f["bb_vv"]("f,b") * l2["bab"]("L,i,a,f");
    tmps_["31_bab_Lovv"]("L,i,a,b") += l2["aaa"]("L,l,d,a") * reused_["29_aabb_voov"]("d,l,i,b");
    tmps_["31_bab_Lovv"]("L,i,a,b") -= l2["bab"]("L,k,a,f") * eri["bbbb_vovo"]("f,i,b,k");
    tmps_["31_bab_Lovv"]("L,i,a,b") -= l2["bab"]("L,k,d,b") * eri["abab_vovo"]("d,i,a,k");
    tmps_["31_bab_Lovv"]("L,i,a,b") -= l2["bab"]("L,k,a,b") * f["bb_oo"]("i,k");
    tmps_["31_bab_Lovv"]("L,i,a,b") += l1["a"]("L,a") * f["bb_ov"]("i,b");
    tmps_["31_bab_Lovv"]("L,i,a,b") -= 0.50 * scalars_["1"] * l2["bab"]("L,i,a,b");
    tmps_["31_bab_Lovv"]("L,i,a,b") += eri["abab_vvvv"]("e,f,a,b") * l2["bab"]("L,i,e,f");
    tmps_["31_bab_Lovv"]("L,i,a,b") -= l2["bab"]("L,k,a,f") * reused_["34_bbbb_voov"]("f,k,i,b");
    tmps_["31_bab_Lovv"]("L,i,a,b") += l2["aaa"]("L,l,d,a") * eri["abba_vovo"]("d,i,b,l");
    tmps_["31_bab_Lovv"]("L,i,a,b") += eri["abab_vovv"]("d,i,a,b") * l1["a"]("L,d");
    tmps_["31_bab_Lovv"]("L,i,a,b") -= l2["bab"]("L,k,a,b") * reused_["21_bb_oo"]("k,i");
    tmps_["31_bab_Lovv"]("L,i,a,b") -= 0.50 * l2["bab"]("L,i,d,b") * reused_["12_aa_vv"]("a,d");
    tmps_["31_bab_Lovv"]("L,i,a,b") -= l2["aaa"]("L,l,d,a") * reused_["33_aabb_voov"]("d,l,i,b");
    tmps_["31_bab_Lovv"]("L,i,a,b") -= 0.50 * tmps_["29_abb_Lvov"]("L,a,i,b");
    tmps_["31_bab_Lovv"]("L,i,a,b") += 0.50 * tmps_["30_abb_Lvvo"]("L,a,b,i");
    tmps_["31_bab_Lovv"]("L,i,a,b") += eri["abab_oovv"]("j,i,a,c") * tmps_["28_bba_Lvvo"]("L,b,c,j");
    tmps_["30_abb_Lvvo"].~TArrayD();
    tmps_["29_abb_Lvov"].~TArrayD();
    tmps_["28_bba_Lvvo"].~TArrayD();
    tmps_["5_a_Lo"].~TArrayD();

    // sigmal2_abb  = +1.00 <k,i||a,d>_abab t2_abab(c,d,k,j) l2_bab(j,c,b)
    //             += +0.50 <k,i||a,b>_abab t2_aaaa(d,c,j,k) l2_aaa(j,d,c)
    //             += -0.50 <k,i||a,b>_abab t2_abab(d,c,k,j) l2_bab(j,d,c)
    //             += -0.50 <k,i||a,b>_abab t2_abab(c,d,k,j) l2_bab(j,c,d)
    //             += -0.50 <i,k||c,d>_bbbb t2_bbbb(c,d,j,k) l2_bab(j,a,b)
    //             += +1.00 <k,i||d,b>_abab t2_abab(d,c,k,j) l2_bab(j,a,c)
    //             += -0.50 <k,j||d,b>_bbbb t2_bbbb(d,c,k,j) l2_bab(i,a,c)
    //             += -0.50 <k,j||d,b>_abab t2_abab(d,c,k,j) l2_bab(i,a,c)
    //             += -0.50 <j,k||d,b>_abab t2_abab(d,c,j,k) l2_bab(i,a,c)
    //             += -0.50 <k,j||a,d>_abab t2_abab(c,d,k,j) l2_bab(i,c,b)
    //             += -0.50 <j,k||a,d>_abab t2_abab(c,d,j,k) l2_bab(i,c,b)
    //             += +1.00 f_bb(c,b) l2_bab(i,a,c)
    //             += +1.00 <k,i||d,b>_abab t2_aaaa(d,c,j,k) l2_aaa(j,c,a)
    //             += +1.00 <i,c||b,j>_bbbb l2_bab(j,a,c)
    //             += -1.00 <c,i||a,j>_abab l2_bab(j,c,b)
    //             += -1.00 f_bb(i,j) l2_bab(j,a,b)
    //             += +1.00 f_bb(i,b) l1_a(a)
    //             += -0.50 <k,j||k,j>_bbbb l2_bab(i,a,b)
    //             += -0.50 <k,j||k,j>_aaaa l2_bab(i,a,b)
    //             += -0.50 <k,j||k,j>_abab l2_bab(i,a,b)
    //             += -0.50 <j,k||j,k>_abab l2_bab(i,a,b)
    //             += +0.25 <k,j||c,d>_bbbb t2_bbbb(c,d,k,j) l2_bab(i,a,b)
    //             += +0.25 <k,j||c,d>_abab t2_abab(c,d,k,j) l2_bab(i,a,b)
    //             += +0.25 <j,k||c,d>_abab t2_abab(c,d,j,k) l2_bab(i,a,b)
    //             += +0.25 <k,j||d,c>_abab t2_abab(d,c,k,j) l2_bab(i,a,b)
    //             += +0.25 <j,k||d,c>_abab t2_abab(d,c,j,k) l2_bab(i,a,b)
    //             += +1.00 f_bb(j,j) l2_bab(i,a,b)
    //             += +0.25 <k,j||c,d>_aaaa t2_aaaa(c,d,k,j) l2_bab(i,a,b)
    //             += +1.00 f_aa(j,j) l2_bab(i,a,b)
    //             += +0.50 <d,c||a,b>_abab l2_bab(i,d,c)
    //             += +0.50 <c,d||a,b>_abab l2_bab(i,c,d)
    //             += +1.00 <i,k||d,b>_bbbb t2_bbbb(d,c,j,k) l2_bab(j,a,c)
    //             += -1.00 <c,i||j,b>_abab l2_aaa(j,c,a)
    //             += +1.00 <c,i||a,b>_abab l1_a(c)
    //             += -0.50 <k,i||c,d>_abab t2_abab(c,d,k,j) l2_bab(j,a,b)
    //             += -0.50 <k,i||d,c>_abab t2_abab(d,c,k,j) l2_bab(j,a,b)
    //             += -0.50 <k,j||d,a>_aaaa t2_aaaa(d,c,k,j) l2_bab(i,c,b)
    //             += +1.00 <i,k||d,b>_bbbb t2_abab(c,d,j,k) l2_aaa(j,c,a)
    //             += +1.00 f_aa(c,a) l2_bab(i,c,b)
    //             += +0.25 <k,j||a,b>_abab t2_abab(d,c,k,j) l2_bab(i,d,c)
    //             += +0.25 <j,k||a,b>_abab t2_abab(d,c,j,k) l2_bab(i,d,c)
    //             += +0.25 <k,j||a,b>_abab t2_abab(c,d,k,j) l2_bab(i,c,d)
    //             += +0.25 <j,k||a,b>_abab t2_abab(c,d,j,k) l2_bab(i,c,d)
    sigmal2_abb("L,a,b,i")  = tmps_["31_bab_Lovv"]("L,i,a,b");
    tmps_["31_bab_Lovv"].~TArrayD();

    // flops: o1v2L1  = o2v3L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["32_abb_Lvvo"]("L,a,b,i")  = l2["aaa"]("L,j,c,a") * t2["abab"]("c,b,j,i");

    // flops: o1v2L1  = o0v2 o1v3L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1
    //  mems: o1v2L1  = o0v2 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["33_aaa_Lovv"]("L,i,a,b")  = -0.50 * (reused_["12_aa_vv"]("a,d") + -2.00 * f["aa_vv"]("d,a")) * l2["aaa"]("L,i,d,b");
    tmps_["33_aaa_Lovv"]("L,i,a,b") -= l2["aaa"]("L,k,d,b") * eri["aaaa_vovo"]("d,i,a,k");
    tmps_["33_aaa_Lovv"]("L,i,a,b") -= l2["aaa"]("L,k,d,b") * reused_["36_aaaa_voov"]("d,k,i,a");
    tmps_["33_aaa_Lovv"]("L,i,a,b") -= l2["bab"]("L,l,b,e") * reused_["37_bbaa_voov"]("e,l,i,a");
    tmps_["33_aaa_Lovv"]("L,i,a,b") -= l1["a"]("L,b") * f["aa_ov"]("i,a");
    tmps_["33_aaa_Lovv"]("L,i,a,b") += l2["bab"]("L,l,b,e") * reused_["20_bbaa_voov"]("e,l,i,a");
    tmps_["33_aaa_Lovv"]("L,i,a,b") -= l2["aaa"]("L,i,d,b") * reused_["4_aa_vv"]("d,a");
    tmps_["33_aaa_Lovv"]("L,i,a,b") += l2["bab"]("L,l,b,e") * eri["baab_vovo"]("e,i,a,l");
    tmps_["33_aaa_Lovv"]("L,i,a,b") += eri["abab_oovv"]("i,j,a,c") * tmps_["32_abb_Lvvo"]("L,b,c,j");
    tmps_["32_abb_Lvvo"].~TArrayD();

    // sigmal2_aaa += +1.00 P(a,b) <i,k||a,d>_abab t2_abab(c,d,j,k) l2_aaa(j,c,b)
    //             += -0.50 P(a,b) <k,j||d,a>_aaaa t2_aaaa(d,c,k,j) l2_aaa(i,c,b)
    //             += +1.00 P(a,b) f_aa(c,a) l2_aaa(i,c,b)
    //             += +1.00 P(a,b) <i,c||a,j>_aaaa l2_aaa(j,c,b)
    //             += +1.00 P(a,b) <i,k||d,a>_aaaa t2_aaaa(d,c,j,k) l2_aaa(j,c,b)
    //             += +1.00 P(a,b) <i,k||d,a>_aaaa t2_abab(d,c,k,j) l2_bab(j,b,c)
    //             += -1.00 P(a,b) f_aa(i,a) l1_a(b)
    //             += +1.00 P(a,b) <i,k||a,d>_abab t2_bbbb(d,c,j,k) l2_bab(j,b,c)
    //             += -0.50 P(a,b) <k,j||a,d>_abab t2_abab(c,d,k,j) l2_aaa(i,c,b)
    //             += -0.50 P(a,b) <j,k||a,d>_abab t2_abab(c,d,j,k) l2_aaa(i,c,b)
    //             += -1.00 P(a,b) <i,c||a,j>_abab l2_bab(j,b,c)
    sigmal2_aaa("L,a,b,i") += tmps_["33_aaa_Lovv"]("L,i,a,b");
    sigmal2_aaa("L,a,b,i") -= tmps_["33_aaa_Lovv"]("L,i,b,a");
    tmps_["33_aaa_Lovv"].~TArrayD();
}