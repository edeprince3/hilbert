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


    // flops: o1v2L1  = o2v3L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["11_abb_Lvov"]("R,a,i,b")  = eri["abab_oovv"]("j,i,c,b") * r2["aaa"]("R,c,a,j");

    // flops: o1v2L1  = o1v3L1 o2v3L1 o2v3L1 o1v3L1 o1v3L1 o1v2L1 o2v3L1 o1v3L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["13_aaa_Lvov"]("R,a,i,b")  = -1.00 * reused_["2_aa_vv"]("a,d") * r2["aaa"]("R,d,b,i");
    tmps_["13_aaa_Lvov"]("R,a,i,b") -= eri["aaaa_vovo"]("a,l,e,i") * r2["aaa"]("R,e,b,l");
    tmps_["13_aaa_Lvov"]("R,a,i,b") -= reused_["24_aabb_voov"]("a,i,m,f") * r2["abb"]("R,b,f,m");
    tmps_["13_aaa_Lvov"]("R,a,i,b") += f["aa_vv"]("a,e") * r2["aaa"]("R,e,b,i");
    tmps_["13_aaa_Lvov"]("R,a,i,b") += r1["a"]("R,d") * reused_["22_aaaa_vvvo"]("a,d,b,i");
    tmps_["13_aaa_Lvov"]("R,a,i,b") -= 0.50 * r1["a"]("R,b") * reused_["11_aa_ov"]("i,a");
    tmps_["13_aaa_Lvov"]("R,a,i,b") -= reused_["23_aaaa_voov"]("a,i,k,d") * r2["aaa"]("R,d,b,k");
    tmps_["13_aaa_Lvov"]("R,a,i,b") += 0.50 * reused_["1_aa_vv"]("a,d") * r2["aaa"]("R,d,b,i");
    tmps_["13_aaa_Lvov"]("R,a,i,b") += r1["a"]("R,d") * reused_["4_aaaa_vvvo"]("a,d,b,i");
    tmps_["13_aaa_Lvov"]("R,a,i,b") += reused_["25_aabb_voov"]("a,i,m,f") * r2["abb"]("R,b,f,m");
    tmps_["13_aaa_Lvov"]("R,a,i,b") += eri["abba_vovo"]("a,j,c,i") * r2["abb"]("R,b,c,j");
    tmps_["13_aaa_Lvov"]("R,a,i,b") += tmps_["12_aaa_Lvov"]("R,a,i,b");
    tmps_["13_aaa_Lvov"]("R,a,i,b") += t2["abab"]("a,c,i,j") * tmps_["11_abb_Lvov"]("R,b,j,c");
    tmps_["12_aaa_Lvov"].~TArrayD();
    tmps_["11_abb_Lvov"].~TArrayD();

    // sigmar2_aaa  = +1.00 P(a,b) <k,j||d,c>_abab r2_aaa(d,b,k) t2_abab(a,c,i,j)
    //             += +1.00 P(a,b) <a,j||d,c>_abab r1_a(d) t2_abab(b,c,i,j)
    //             += -0.50 P(a,b) <k,j||c,d>_aaaa r2_aaa(d,b,i) t2_aaaa(c,a,k,j)
    //             += +1.00 P(a,b) <k,j||c,d>_aaaa r2_aaa(d,b,k) t2_aaaa(c,a,i,j)
    //             += +0.50 P(a,b) <k,j||c,i>_aaaa r1_a(b) t2_aaaa(c,a,k,j)
    //             += -0.50 P(a,b) <a,j||c,d>_abab r1_a(b) t2_abab(c,d,i,j)
    //             += -0.50 P(a,b) <a,j||d,c>_abab r1_a(b) t2_abab(d,c,i,j)
    //             += +0.50 P(a,b) <j,a||c,d>_aaaa r1_a(b) t2_aaaa(c,d,i,j)
    //             += +1.00 P(a,b) f_aa(j,c) r1_a(b) t2_aaaa(c,a,i,j)
    //             += +0.50 P(a,b) <k,j||i,c>_abab r1_a(b) t2_abab(a,c,k,j)
    //             += +0.50 P(a,b) <j,k||i,c>_abab r1_a(b) t2_abab(a,c,j,k)
    //             += -1.00 P(a,b) f_bb(j,c) r1_a(b) t2_abab(a,c,i,j)
    //             += -1.00 P(a,b) <j,a||c,d>_aaaa r1_a(d) t2_aaaa(c,b,i,j)
    //             += +1.00 P(a,b) f_aa(a,c) r2_aaa(c,b,i)
    //             += +1.00 P(a,b) <k,j||c,d>_bbbb r2_abb(b,d,k) t2_abab(a,c,i,j)
    //             += +1.00 P(a,b) <j,k||c,d>_abab r2_abb(b,d,k) t2_aaaa(c,a,i,j)
    //             += -1.00 P(a,b) <a,j||i,c>_abab r2_abb(b,c,j)
    //             += +1.00 P(a,b) <j,a||c,i>_aaaa r2_aaa(c,b,j)
    //             += -0.50 P(a,b) <k,j||d,c>_abab r2_aaa(d,b,i) t2_abab(a,c,k,j)
    //             += -0.50 P(a,b) <j,k||d,c>_abab r2_aaa(d,b,i) t2_abab(a,c,j,k)
    //             += -1.00 P(a,b) f_aa(a,i) r1_a(b)
    sigmar2_aaa("R,a,b,i")  = tmps_["13_aaa_Lvov"]("R,a,i,b");
    sigmar2_aaa("R,a,b,i") -= tmps_["13_aaa_Lvov"]("R,b,i,a");
    tmps_["13_aaa_Lvov"].~TArrayD();

    // flops: o3v0L1  = o3v2L1
    //  mems: o3v0L1  = o3v0L1
    tmps_["14_aaa_Looo"]("L,i,j,k")  = l2["aaa"]("L,i,a,b") * t2["aaaa"]("a,b,j,k");

    // flops: o1v2L1  = o3v2L1 o2v2L1 o2v2L1 o1v4L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["15_aaa_Lvvo"]("L,a,b,i")  = 0.50 * eri["aaaa_oovv"]("k,j,a,b") * tmps_["14_aaa_Looo"]("L,i,j,k");
    tmps_["15_aaa_Lvvo"]("L,a,b,i") += 2.00 * l2["aaa"]("L,k,a,b") * f["aa_oo"]("i,k");
    tmps_["15_aaa_Lvvo"]("L,a,b,i") += 2.00 * l2["aaa"]("L,k,a,b") * reused_["27_aa_oo"]("k,i");
    tmps_["15_aaa_Lvvo"]("L,a,b,i") += eri["aaaa_vvvv"]("d,c,a,b") * l2["aaa"]("L,i,c,d");
    tmps_["15_aaa_Lvvo"]("L,a,b,i") += scalars_["1"] * l2["aaa"]("L,i,a,b");
    tmps_["15_aaa_Lvvo"]("L,a,b,i") -= 2.00 * eri["aaaa_vovv"]("d,i,a,b") * l1["a"]("L,d");
    tmps_["15_aaa_Lvvo"]("L,a,b,i") += l2["aaa"]("L,k,a,b") * reused_["26_aa_oo"]("i,k");
    tmps_["15_aaa_Lvvo"]("L,a,b,i") += eri["aaaa_oovv"]("i,j,a,b") * tmps_["3_a_Lo"]("L,j");
    tmps_["14_aaa_Looo"].~TArrayD();

    // sigmal2_aaa  = +0.50 <d,c||a,b>_aaaa l2_aaa(i,d,c)
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
    //             += -0.50 <i,k||c,d>_abab t2_abab(c,d,j,k) l2_aaa(j,a,b)
    //             += -0.50 <i,k||d,c>_abab t2_abab(d,c,j,k) l2_aaa(j,a,b)
    //             += -1.00 <i,c||a,b>_aaaa l1_a(c)
    //             += -1.00 f_aa(i,j) l2_aaa(j,a,b)
    //             += -0.50 <i,k||c,d>_aaaa t2_aaaa(c,d,j,k) l2_aaa(j,a,b)
    //             += +0.25 <k,j||a,b>_aaaa t2_aaaa(d,c,k,j) l2_aaa(i,d,c)
    //             += -0.50 <i,k||a,b>_aaaa t2_aaaa(d,c,j,k) l2_aaa(j,d,c)
    //             += +0.50 <i,k||a,b>_aaaa t2_abab(d,c,k,j) l2_bab(j,d,c)
    //             += +0.50 <i,k||a,b>_aaaa t2_abab(c,d,k,j) l2_bab(j,c,d)
    sigmal2_aaa("L,a,b,i")  = -0.50 * tmps_["15_aaa_Lvvo"]("L,a,b,i");
    tmps_["15_aaa_Lvvo"].~TArrayD();

    // flops: o3v0L1  = o3v2L1
    //  mems: o3v0L1  = o3v0L1
    tmps_["16_aaa_Looo"]("R,i,j,k")  = eri["aaaa_oovv"]("j,k,a,b") * r2["aaa"]("R,a,b,i");

    // flops: o1v2L1  = o3v2L1 o2v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o1v4L1 o1v2L1 o2v2L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["17_aaa_Lvvo"]("R,a,b,i")  = -0.25 * t2["aaaa"]("a,b,k,j") * tmps_["16_aaa_Looo"]("R,i,j,k");
    tmps_["17_aaa_Lvvo"]("R,a,b,i") += t2["aaaa"]("a,b,i,j") * tmps_["9_a_Lo"]("R,j");
    tmps_["17_aaa_Lvvo"]("R,a,b,i") -= 0.50 * scalars_["1"] * r2["aaa"]("R,a,b,i");
    tmps_["17_aaa_Lvvo"]("R,a,b,i") -= r2["aaa"]("R,a,b,j") * f["aa_oo"]("j,i");
    tmps_["17_aaa_Lvvo"]("R,a,b,i") += 0.50 * eri["aaaa_vvvv"]("a,b,c,d") * r2["aaa"]("R,c,d,i");
    tmps_["17_aaa_Lvvo"]("R,a,b,i") += 0.50 * r2["aaa"]("R,a,b,k") * reused_["28_aa_oo"]("i,k");
    tmps_["17_aaa_Lvvo"]("R,a,b,i") += eri["aaaa_vvvo"]("a,b,c,i") * r1["a"]("R,c");
    tmps_["17_aaa_Lvvo"]("R,a,b,i") -= r2["aaa"]("R,a,b,k") * reused_["27_aa_oo"]("i,k");
    tmps_["17_aaa_Lvvo"]("R,a,b,i") -= 0.50 * r1["a"]("R,c") * reused_["8_aaaa_vovv"]("c,i,a,b");
    tmps_["16_aaa_Looo"].~TArrayD();
    tmps_["9_a_Lo"].~TArrayD();
}