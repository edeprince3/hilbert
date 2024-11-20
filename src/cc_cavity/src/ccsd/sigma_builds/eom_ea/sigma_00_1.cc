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

void hilbert::EOM_EA_CCSD::sigma_ea_00_1() {

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

    // flops: o0v1L1  = o1v3L1 o0v2L1 o0v1L1 o0v1L1 o1v3L1 o0v1L1 o1v2L1 o0v1L1 o0v2L1 o0v1L1 o0v1L1 o1v2L1 o0v1L1 o0v2L1 o0v1L1
    //  mems: o0v1L1  = o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1
    tmps_["1_a_Lv"]("R,a")  = -1.00 * eri["aaaa_vovv"]("a,i,b,c") * r2["aaa"]("R,b,c,i");
    tmps_["1_a_Lv"]("R,a") -= r1["a"]("R,c") * reused_["1_aa_vv"]("a,c");
    tmps_["1_a_Lv"]("R,a") += scalars_["1"] * r1["a"]("R,a");
    tmps_["1_a_Lv"]("R,a") -= 2.00 * eri["abab_vovv"]("a,j,b,e") * r2["abb"]("R,b,e,j");
    tmps_["1_a_Lv"]("R,a") -= 2.00 * f["bb_ov"]("j,d") * r2["abb"]("R,a,d,j");
    tmps_["1_a_Lv"]("R,a") += 2.00 * r1["a"]("R,c") * reused_["2_aa_vv"]("a,c");
    tmps_["1_a_Lv"]("R,a") += 2.00 * f["aa_ov"]("i,b") * r2["aaa"]("R,b,a,i");
    tmps_["1_a_Lv"]("R,a") -= 2.00 * f["aa_vv"]("a,b") * r1["a"]("R,b");

    // sigmar1_a  = -0.50 <j,i||j,i>_bbbb r1_a(a)
    //           += -0.50 <j,i||j,i>_aaaa r1_a(a)
    //           += -0.50 <j,i||j,i>_abab r1_a(a)
    //           += -0.50 <i,j||i,j>_abab r1_a(a)
    //           += +0.25 <j,i||b,c>_bbbb r1_a(a) t2_bbbb(b,c,j,i)
    //           += +0.25 <j,i||b,c>_abab r1_a(a) t2_abab(b,c,j,i)
    //           += +0.25 <i,j||b,c>_abab r1_a(a) t2_abab(b,c,i,j)
    //           += +0.25 <j,i||c,b>_abab r1_a(a) t2_abab(c,b,j,i)
    //           += +0.25 <i,j||c,b>_abab r1_a(a) t2_abab(c,b,i,j)
    //           += +1.00 f_bb(i,i) r1_a(a)
    //           += +0.25 <j,i||b,c>_aaaa r1_a(a) t2_aaaa(b,c,j,i)
    //           += +1.00 f_aa(i,i) r1_a(a)
    //           += -0.50 <j,i||b,c>_aaaa r1_a(c) t2_aaaa(b,a,j,i)
    //           += +0.50 <a,i||b,c>_abab r2_abb(b,c,i)
    //           += +0.50 <a,i||c,b>_abab r2_abb(c,b,i)
    //           += +1.00 f_bb(i,b) r2_abb(a,b,i)
    //           += -0.50 <j,i||c,b>_abab r1_a(c) t2_abab(a,b,j,i)
    //           += -0.50 <i,j||c,b>_abab r1_a(c) t2_abab(a,b,i,j)
    //           += -0.50 <i,a||b,c>_aaaa r2_aaa(b,c,i)
    //           += -1.00 f_aa(i,b) r2_aaa(b,a,i)
    //           += +1.00 f_aa(a,b) r1_a(b)
    sigmar1_a("R,a")  = -0.50 * tmps_["1_a_Lv"]("R,a");
    tmps_["1_a_Lv"].~TArrayD();

    // flops: o0v1L1  = o1v3L1
    //  mems: o0v1L1  = o0v1L1
    tmps_["4_a_Lv"]("L,a")  = -4.00 * l2["bab"]("L,i,b,c") * reused_["7_baab_vvvo"]("c,a,b,i");

    // flops: o0v1L1  = o0v2L1
    //  mems: o0v1L1  = o0v1L1
    tmps_["2_a_Lv"]("R,a")  = -2.00 * f["aa_vv"]("a,b") * r1["a"]("R,b");
    tmps_["2_a_Lv"].~TArrayD();

    // flops: o1v0L1  = o2v2L1 o2v2L1 o1v0L1
    //  mems: o1v0L1  = o1v0L1 o1v0L1 o1v0L1
    tmps_["3_a_Lo"]("L,i")  = -2.00 * l2["bab"]("L,k,a,c") * t2["abab"]("a,c,i,k");
    tmps_["3_a_Lo"]("L,i") += l2["aaa"]("L,j,a,b") * t2["aaaa"]("a,b,j,i");

    // flops: o0v1L1  = o1v3 o1v3 o1v3L1 o1v3L1 o0v1L1 o0v2L1 o0v1L1 o0v2L1 o0v1L1 o1v3L1 o0v1L1 o1v3L1 o0v1L1 o0v1L1 o0v1L1 o1v3L1 o0v1L1 o1v2L1 o0v1L1 o1v2L1 o0v1L1 o1v2L1 o0v1L1 o1v2L1 o0v1L1 o1v3L1 o0v1L1 o0v2L1 o0v1L1 o0v1L1 o1v1L1 o0v1L1
    //  mems: o0v1L1  = o1v3 o1v3 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1 o0v1L1
    tmps_["5_a_Lv"]("L,a")  = -0.50 * (4.00 * reused_["3_aaaa_vovv"]("c,j,b,a") + reused_["8_aaaa_vovv"]("a,j,b,c") + 2.00 * eri["aaaa_vvvo"]("c,b,a,j")) * l2["aaa"]("L,j,b,c");
    tmps_["5_a_Lv"]("L,a") += 2.00 * eri["abab_vvvo"]("b,d,a,k") * l2["bab"]("L,k,b,d");
    tmps_["5_a_Lv"]("L,a") -= 2.00 * l1["a"]("L,c") * reused_["2_aa_vv"]("c,a");
    tmps_["5_a_Lv"]("L,a") -= l1["a"]("L,c") * reused_["10_aa_vv"]("a,c");
    tmps_["5_a_Lv"]("L,a") += 2.00 * l2["bab"]("L,k,b,d") * reused_["5_aabb_vvvo"]("b,a,d,k");
    tmps_["5_a_Lv"]("L,a") += 2.00 * l2["aaa"]("L,j,b,c") * reused_["4_aaaa_vvvo"]("b,a,c,j");
    tmps_["5_a_Lv"]("L,a") -= scalars_["1"] * l1["a"]("L,a");
    tmps_["5_a_Lv"]("L,a") += 2.00 * l2["bab"]("L,k,b,d") * reused_["9_abab_vovv"]("a,k,b,d");
    tmps_["5_a_Lv"]("L,a") += 2.00 * l2["bab"]("L,k,a,d") * reused_["12_bb_vo"]("d,k");
    tmps_["5_a_Lv"]("L,a") -= 2.00 * f["aa_vo"]("c,j") * l2["aaa"]("L,j,c,a");
    tmps_["5_a_Lv"]("L,a") += 2.00 * f["bb_vo"]("d,k") * l2["bab"]("L,k,a,d");
    tmps_["5_a_Lv"]("L,a") -= l2["aaa"]("L,j,c,a") * reused_["11_aa_ov"]("j,c");
    tmps_["5_a_Lv"]("L,a") -= 2.00 * l2["bab"]("L,k,b,d") * reused_["6_aabb_vvvo"]("b,a,d,k");
    tmps_["5_a_Lv"]("L,a") += 2.00 * f["aa_vv"]("c,a") * l1["a"]("L,c");
    tmps_["5_a_Lv"]("L,a") -= 0.50 * tmps_["4_a_Lv"]("L,a");
    tmps_["5_a_Lv"]("L,a") += f["aa_ov"]("i,a") * tmps_["3_a_Lo"]("L,i");

    // sigmal1_a  = +0.50 f_aa(j,a) t2_aaaa(c,b,i,j) l2_aaa(i,c,b)
    //           += -0.50 f_aa(j,a) t2_abab(c,b,j,i) l2_bab(i,c,b)
    //           += -0.50 f_aa(j,a) t2_abab(b,c,j,i) l2_bab(i,b,c)
    //           += +0.25 <k,j||a,i>_aaaa t2_aaaa(c,b,k,j) l2_aaa(i,c,b)
    //           += -1.00 <j,c||d,a>_aaaa t2_aaaa(d,b,i,j) l2_aaa(i,c,b)
    //           += +0.50 <c,b||a,i>_aaaa l2_aaa(i,c,b)
    //           += +0.50 <c,b||a,i>_abab l2_bab(i,c,b)
    //           += +0.50 <b,c||a,i>_abab l2_bab(i,b,c)
    //           += -0.50 <j,i||a,c>_abab t2_abab(b,c,j,i) l1_a(b)
    //           += -0.50 <i,j||a,c>_abab t2_abab(b,c,i,j) l1_a(b)
    //           += -0.50 <j,i||c,a>_aaaa t2_aaaa(c,b,j,i) l1_a(b)
    //           += +1.00 <j,c||d,a>_aaaa t2_abab(d,b,j,i) l2_bab(i,c,b)
    //           += +1.00 <c,j||a,d>_abab t2_abab(b,d,i,j) l2_aaa(i,c,b)
    //           += -0.50 <j,i||j,i>_bbbb l1_a(a)
    //           += -0.50 <j,i||j,i>_aaaa l1_a(a)
    //           += -0.50 <j,i||j,i>_abab l1_a(a)
    //           += -0.50 <i,j||i,j>_abab l1_a(a)
    //           += +0.25 <j,i||b,c>_bbbb t2_bbbb(b,c,j,i) l1_a(a)
    //           += +0.25 <j,i||b,c>_abab t2_abab(b,c,j,i) l1_a(a)
    //           += +0.25 <i,j||b,c>_abab t2_abab(b,c,i,j) l1_a(a)
    //           += +0.25 <j,i||c,b>_abab t2_abab(c,b,j,i) l1_a(a)
    //           += +0.25 <i,j||c,b>_abab t2_abab(c,b,i,j) l1_a(a)
    //           += +1.00 f_bb(i,i) l1_a(a)
    //           += +0.25 <j,i||b,c>_aaaa t2_aaaa(b,c,j,i) l1_a(a)
    //           += +1.00 f_aa(i,i) l1_a(a)
    //           += +0.25 <k,j||a,i>_abab t2_abab(c,b,k,j) l2_bab(i,c,b)
    //           += +0.25 <j,k||a,i>_abab t2_abab(c,b,j,k) l2_bab(i,c,b)
    //           += +0.25 <k,j||a,i>_abab t2_abab(b,c,k,j) l2_bab(i,b,c)
    //           += +0.25 <j,k||a,i>_abab t2_abab(b,c,j,k) l2_bab(i,b,c)
    //           += +1.00 f_aa(j,c) t2_abab(c,b,j,i) l2_bab(i,a,b)
    //           += -0.50 <j,b||c,d>_bbbb t2_bbbb(c,d,i,j) l2_bab(i,a,b)
    //           += -0.50 <k,j||c,i>_abab t2_abab(c,b,k,j) l2_bab(i,a,b)
    //           += -0.50 <j,k||c,i>_abab t2_abab(c,b,j,k) l2_bab(i,a,b)
    //           += -1.00 f_bb(j,c) t2_bbbb(c,b,i,j) l2_bab(i,a,b)
    //           += -0.50 <k,j||c,i>_bbbb t2_bbbb(c,b,k,j) l2_bab(i,a,b)
    //           += +0.50 <j,b||c,d>_abab t2_abab(c,d,j,i) l2_bab(i,a,b)
    //           += +0.50 <j,b||d,c>_abab t2_abab(d,c,j,i) l2_bab(i,a,b)
    //           += -1.00 f_aa(b,i) l2_aaa(i,b,a)
    //           += +1.00 f_bb(b,i) l2_bab(i,a,b)
    //           += +0.50 <k,j||c,i>_aaaa t2_aaaa(c,b,k,j) l2_aaa(i,b,a)
    //           += -0.50 <b,j||c,d>_abab t2_abab(c,d,i,j) l2_aaa(i,b,a)
    //           += -0.50 <b,j||d,c>_abab t2_abab(d,c,i,j) l2_aaa(i,b,a)
    //           += +0.50 <j,b||c,d>_aaaa t2_aaaa(c,d,i,j) l2_aaa(i,b,a)
    //           += +1.00 f_aa(j,c) t2_aaaa(c,b,i,j) l2_aaa(i,b,a)
    //           += +0.50 <k,j||i,c>_abab t2_abab(b,c,k,j) l2_aaa(i,b,a)
    //           += +0.50 <j,k||i,c>_abab t2_abab(b,c,j,k) l2_aaa(i,b,a)
    //           += -1.00 f_bb(j,c) t2_abab(b,c,i,j) l2_aaa(i,b,a)
    //           += -1.00 <c,j||a,d>_abab t2_bbbb(d,b,i,j) l2_bab(i,c,b)
    //           += +1.00 f_aa(b,a) l1_a(b)
    //           += -1.00 <j,c||a,d>_abab t2_abab(b,d,j,i) l2_bab(i,b,c)
    sigmal1_a("L,a")  = 0.50 * tmps_["5_a_Lv"]("L,a");
    tmps_["5_a_Lv"].~TArrayD();

    // flops: o1v0L1  = o2v2L1 o2v2L1 o1v1L1 o1v0L1 o1v0L1
    //  mems: o1v0L1  = o1v0L1 o1v0L1 o1v0L1 o1v0L1 o1v0L1
    tmps_["9_a_Lo"]("R,i")  = 0.50 * eri["aaaa_oovv"]("i,k,a,c") * r2["aaa"]("R,a,c,k");
    tmps_["9_a_Lo"]("R,i") += eri["abab_oovv"]("i,j,a,b") * r2["abb"]("R,a,b,j");
    tmps_["9_a_Lo"]("R,i") += f["aa_ov"]("i,a") * r1["a"]("R,a");
}