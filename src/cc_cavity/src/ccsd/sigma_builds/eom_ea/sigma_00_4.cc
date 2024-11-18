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
    };

    TArrayMap r2 = {
            {"aaa", evec_blks_["r2_aaa"]},
            {"abb", evec_blks_["r2_abb"]},
    };

    /// reference left operators

    TArrayMap l1 = {
            {"a", evec_blks_["l1_a"]},
    };
    TArrayMap l2 = {
            {"aaa", evec_blks_["l2_aaa"]},
            {"bab", evec_blks_["l2_bab"]},
    };

    /// reference sigma vectors
    auto &sigmar1_a = sigvec_blks_["sigmar1_a"];
    auto &sigmal1_a = sigvec_blks_["sigmal1_a"];

    auto &sigmar2_aaa = sigvec_blks_["sigmar2_aaa"];
    auto &sigmar2_abb = sigvec_blks_["sigmar2_abb"];
    auto &sigmal2_aaa = sigvec_blks_["sigmal2_aaa"];
    auto &sigmal2_abb = sigvec_blks_["sigmal2_abb"];


    // sigmar2_aaa += +1.00 f_aa(j,c) r1_a(c) t2_aaaa(a,b,i,j)
    //             += +0.50 <j,k||c,d>_abab r2_abb(c,d,k) t2_aaaa(a,b,i,j)
    //             += +0.50 <j,k||d,c>_abab r2_abb(d,c,k) t2_aaaa(a,b,i,j)
    //             += -0.50 <k,j||c,d>_aaaa r2_aaa(c,d,k) t2_aaaa(a,b,i,j)
    //             += +0.25 <k,j||c,d>_aaaa r2_aaa(c,d,i) t2_aaaa(a,b,k,j)
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
    //             += -1.00 f_aa(j,i) r2_aaa(a,b,j)
    //             += +0.50 <a,b||c,d>_aaaa r2_aaa(c,d,i)
    //             += -0.50 <k,j||c,d>_aaaa r2_aaa(a,b,k) t2_aaaa(c,d,i,j)
    //             += +1.00 <a,b||c,i>_aaaa r1_a(c)
    //             += -0.50 <k,j||c,d>_abab r2_aaa(a,b,k) t2_abab(c,d,i,j)
    //             += -0.50 <k,j||d,c>_abab r2_aaa(a,b,k) t2_abab(d,c,i,j)
    //             += +0.50 <k,j||c,i>_aaaa r1_a(c) t2_aaaa(a,b,k,j)
    sigmar2_aaa("R,a,b,i") += tmps_["17_aaa_Lvvo"]("R,a,b,i");
    tmps_["17_aaa_Lvvo"].~TArrayD();

    // flops: o1v2L1  = o2v3L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["18_abb_Lvvo"]("L,a,b,i")  = l2["aaa"]("L,j,c,a") * t2["abab"]("c,b,j,i");

    // flops: o1v2L1  = o2v3L1 o2v3L1 o1v3L1 o1v2L1 o2v3L1 o2v3L1 o2v3L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["19_aaa_Lovv"]("L,i,a,b")  = -2.00 * l2["bab"]("L,k,a,d") * eri["baab_vovo"]("d,i,b,k");
    tmps_["19_aaa_Lovv"]("L,i,a,b") -= 2.00 * l2["bab"]("L,k,a,d") * reused_["19_bbaa_voov"]("d,k,i,b");
    tmps_["19_aaa_Lovv"]("L,i,a,b") -= 2.00 * f["aa_vv"]("e,b") * l2["aaa"]("L,i,e,a");
    tmps_["19_aaa_Lovv"]("L,i,a,b") += 2.00 * f["aa_ov"]("i,b") * l1["a"]("L,a");
    tmps_["19_aaa_Lovv"]("L,i,a,b") += 2.00 * l2["bab"]("L,k,a,d") * reused_["30_bbaa_voov"]("d,k,i,b");
    tmps_["19_aaa_Lovv"]("L,i,a,b") += 2.00 * eri["aaaa_vovo"]("e,i,b,l") * l2["aaa"]("L,l,e,a");
    tmps_["19_aaa_Lovv"]("L,i,a,b") += 2.00 * l2["aaa"]("L,l,e,a") * reused_["29_aaaa_voov"]("e,l,i,b");
    tmps_["19_aaa_Lovv"]("L,i,a,b") += reused_["10_aa_vv"]("b,e") * l2["aaa"]("L,i,e,a");
    tmps_["19_aaa_Lovv"]("L,i,a,b") += 2.00 * l2["aaa"]("L,i,e,a") * reused_["2_aa_vv"]("e,b");
    tmps_["19_aaa_Lovv"]("L,i,a,b") -= 2.00 * eri["abab_oovv"]("i,j,b,c") * tmps_["18_abb_Lvvo"]("L,a,c,j");
    tmps_["18_abb_Lvvo"].~TArrayD();

    // sigmal2_aaa += -0.50 P(a,b) <k,j||d,a>_aaaa t2_aaaa(d,c,k,j) l2_aaa(i,c,b)
    //             += +1.00 P(a,b) <i,k||d,a>_aaaa t2_aaaa(d,c,j,k) l2_aaa(j,c,b)
    //             += +1.00 P(a,b) <i,c||a,j>_aaaa l2_aaa(j,c,b)
    //             += +1.00 P(a,b) <i,k||d,a>_aaaa t2_abab(d,c,k,j) l2_bab(j,b,c)
    //             += -1.00 P(a,b) f_aa(i,a) l1_a(b)
    //             += +1.00 P(a,b) f_aa(c,a) l2_aaa(i,c,b)
    //             += +1.00 P(a,b) <i,k||a,d>_abab t2_bbbb(d,c,j,k) l2_bab(j,b,c)
    //             += -0.50 P(a,b) <k,j||a,d>_abab t2_abab(c,d,k,j) l2_aaa(i,c,b)
    //             += -0.50 P(a,b) <j,k||a,d>_abab t2_abab(c,d,j,k) l2_aaa(i,c,b)
    //             += -1.00 P(a,b) <i,c||a,j>_abab l2_bab(j,b,c)
    //             += +1.00 P(a,b) <i,k||a,d>_abab t2_abab(c,d,j,k) l2_aaa(j,c,b)
    sigmal2_aaa("L,a,b,i") -= 0.50 * tmps_["19_aaa_Lovv"]("L,i,b,a");
    sigmal2_aaa("L,a,b,i") += 0.50 * tmps_["19_aaa_Lovv"]("L,i,a,b");
    tmps_["19_aaa_Lovv"].~TArrayD();

    // flops: o1v2L1  = o2v3L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["21_bba_Lvvo"]("L,a,b,i")  = l2["bab"]("L,j,c,a") * t2["abab"]("c,b,i,j");

    // flops: o3v0L1  = o3v2L1
    //  mems: o3v0L1  = o3v0L1
    tmps_["20_abb_Looo"]("L,i,j,k")  = l2["bab"]("L,k,a,b") * t2["abab"]("a,b,i,j");

    // flops: o1v2L1  = o2v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v4L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o3v2L1 o1v2L1 o2v3L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["22_bab_Lovv"]("L,i,a,b")  = -0.50 * l2["bab"]("L,k,a,b") * reused_["34_bb_oo"]("i,k");
    tmps_["22_bab_Lovv"]("L,i,a,b") += l2["bab"]("L,k,a,f") * reused_["18_bbbb_voov"]("f,k,i,b");
    tmps_["22_bab_Lovv"]("L,i,a,b") -= 0.50 * reused_["33_bb_vv"]("b,f") * l2["bab"]("L,i,a,f");
    tmps_["22_bab_Lovv"]("L,i,a,b") += f["bb_vv"]("f,b") * l2["bab"]("L,i,a,f");
    tmps_["22_bab_Lovv"]("L,i,a,b") -= reused_["20_bb_vv"]("b,f") * l2["bab"]("L,i,a,f");
    tmps_["22_bab_Lovv"]("L,i,a,b") += l2["aaa"]("L,l,d,a") * reused_["25_aabb_voov"]("d,l,i,b");
    tmps_["22_bab_Lovv"]("L,i,a,b") -= reused_["2_aa_vv"]("d,a") * l2["bab"]("L,i,d,b");
    tmps_["22_bab_Lovv"]("L,i,a,b") += f["aa_vv"]("d,a") * l2["bab"]("L,i,d,b");
    tmps_["22_bab_Lovv"]("L,i,a,b") -= eri["abab_vovo"]("d,i,a,k") * l2["bab"]("L,k,d,b");
    tmps_["22_bab_Lovv"]("L,i,a,b") -= eri["bbbb_vovo"]("f,i,b,k") * l2["bab"]("L,k,a,f");
    tmps_["22_bab_Lovv"]("L,i,a,b") -= f["bb_oo"]("i,k") * l2["bab"]("L,k,a,b");
    tmps_["22_bab_Lovv"]("L,i,a,b") += f["bb_ov"]("i,b") * l1["a"]("L,a");
    tmps_["22_bab_Lovv"]("L,i,a,b") -= 0.50 * scalars_["1"] * l2["bab"]("L,i,a,b");
    tmps_["22_bab_Lovv"]("L,i,a,b") += eri["abab_vvvv"]("e,f,a,b") * l2["bab"]("L,i,e,f");
    tmps_["22_bab_Lovv"]("L,i,a,b") -= l2["bab"]("L,k,a,f") * reused_["32_bbbb_voov"]("f,k,i,b");
    tmps_["22_bab_Lovv"]("L,i,a,b") += eri["abba_vovo"]("d,i,b,l") * l2["aaa"]("L,l,d,a");
    tmps_["22_bab_Lovv"]("L,i,a,b") += eri["abab_vovv"]("d,i,a,b") * l1["a"]("L,d");
    tmps_["22_bab_Lovv"]("L,i,a,b") -= l2["bab"]("L,k,a,b") * reused_["21_bb_oo"]("k,i");
    tmps_["22_bab_Lovv"]("L,i,a,b") -= l2["aaa"]("L,l,d,a") * reused_["31_aabb_voov"]("d,l,i,b");
    tmps_["22_bab_Lovv"]("L,i,a,b") -= 0.50 * reused_["10_aa_vv"]("a,d") * l2["bab"]("L,i,d,b");
    tmps_["22_bab_Lovv"]("L,i,a,b") += 0.50 * eri["abab_oovv"]("j,i,a,b") * tmps_["3_a_Lo"]("L,j");
    tmps_["22_bab_Lovv"]("L,i,a,b") += eri["abab_oovv"]("j,k,a,b") * tmps_["20_abb_Looo"]("L,j,k,i");
    tmps_["22_bab_Lovv"]("L,i,a,b") += eri["abab_oovv"]("j,i,a,c") * tmps_["21_bba_Lvvo"]("L,b,c,j");
    tmps_["21_bba_Lvvo"].~TArrayD();
    tmps_["20_abb_Looo"].~TArrayD();
    tmps_["3_a_Lo"].~TArrayD();

    // sigmal2_abb  = +1.00 <k,i||a,d>_abab t2_abab(c,d,k,j) l2_bab(j,c,b)
    //             += +1.00 <k,i||d,b>_abab t2_abab(d,c,k,j) l2_bab(j,a,c)
    //             += -0.50 <i,k||c,d>_bbbb t2_bbbb(c,d,j,k) l2_bab(j,a,b)
    //             += -0.50 <k,j||d,b>_bbbb t2_bbbb(d,c,k,j) l2_bab(i,a,c)
    //             += +1.00 f_bb(c,b) l2_bab(i,a,c)
    //             += -0.50 <k,j||d,b>_abab t2_abab(d,c,k,j) l2_bab(i,a,c)
    //             += -0.50 <j,k||d,b>_abab t2_abab(d,c,j,k) l2_bab(i,a,c)
    //             += +1.00 <k,i||d,b>_abab t2_aaaa(d,c,j,k) l2_aaa(j,c,a)
    //             += -0.50 <k,j||a,d>_abab t2_abab(c,d,k,j) l2_bab(i,c,b)
    //             += -0.50 <j,k||a,d>_abab t2_abab(c,d,j,k) l2_bab(i,c,b)
    //             += +1.00 f_aa(c,a) l2_bab(i,c,b)
    //             += -1.00 <c,i||a,j>_abab l2_bab(j,c,b)
    //             += +1.00 <i,c||b,j>_bbbb l2_bab(j,a,c)
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
    //             += +1.00 <i,k||d,b>_bbbb t2_abab(c,d,j,k) l2_aaa(j,c,a)
    //             += -0.50 <k,j||d,a>_aaaa t2_aaaa(d,c,k,j) l2_bab(i,c,b)
    //             += +0.50 <k,i||a,b>_abab t2_aaaa(d,c,j,k) l2_aaa(j,d,c)
    //             += -0.50 <k,i||a,b>_abab t2_abab(d,c,k,j) l2_bab(j,d,c)
    //             += -0.50 <k,i||a,b>_abab t2_abab(c,d,k,j) l2_bab(j,c,d)
    //             += +0.25 <k,j||a,b>_abab t2_abab(d,c,k,j) l2_bab(i,d,c)
    //             += +0.25 <j,k||a,b>_abab t2_abab(d,c,j,k) l2_bab(i,d,c)
    //             += +0.25 <k,j||a,b>_abab t2_abab(c,d,k,j) l2_bab(i,c,d)
    //             += +0.25 <j,k||a,b>_abab t2_abab(c,d,j,k) l2_bab(i,c,d)
    sigmal2_abb("L,a,b,i")  = tmps_["22_bab_Lovv"]("L,i,a,b");
    tmps_["22_bab_Lovv"].~TArrayD();
}