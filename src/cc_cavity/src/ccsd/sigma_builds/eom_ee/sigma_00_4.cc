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

#include "cc_cavity/include/ccsd/eom_ee_ccsd.h"

void hilbert::EOM_EE_CCSD::sigma_ee_00_4() {

    // get integrals
    TArrayMap &eri = cc_wfn_->V_blks_;
    TArrayMap &f = cc_wfn_->F_blks_;
    TArrayMap &Id = cc_wfn_->Id_blks_;


    // extract amplitudes
    TArrayMap t2 = {
            {"aaaa", cc_wfn_->amplitudes_["t2_aaaa"]},
            {"abab", cc_wfn_->amplitudes_["t2_abab"]},
            {"bbbb", cc_wfn_->amplitudes_["t2_bbbb"]}
    };

    /// extract right operators
    TArrayD &r0 = evec_blks_["r0"];
    TArrayMap r1 = {
            {"aa", evec_blks_["r1_aa"]},
            {"bb", evec_blks_["r1_bb"]}
    };

    TArrayMap r2 = {
            {"aaaa", evec_blks_["r2_aaaa"]},
            {"abab", evec_blks_["r2_abab"]},
            {"bbbb", evec_blks_["r2_bbbb"]}
    };

    /// extract left operators
    TArrayD &l0 = evec_blks_["l0"];
    TArrayMap l1 = {
            {"aa", evec_blks_["l1_aa"]},
            {"bb", evec_blks_["l1_bb"]}
    };
    TArrayMap l2 = {
            {"aaaa", evec_blks_["l2_aaaa"]},
            {"abab", evec_blks_["l2_abab"]},
            {"bbbb", evec_blks_["l2_bbbb"]}
    };

    /// extract sigma vectors
    auto &sigmar0 = sigvec_blks_["sigmar0"];
    auto &sigmal0 = sigvec_blks_["sigmal0"];

    auto &sigmar1_aa = sigvec_blks_["sigmar1_aa"];
    auto &sigmar1_bb = sigvec_blks_["sigmar1_bb"];
    auto &sigmal1_aa = sigvec_blks_["sigmal1_aa"];
    auto &sigmal1_bb = sigvec_blks_["sigmal1_bb"];

    auto &sigmar2_aaaa = sigvec_blks_["sigmar2_aaaa"];
    auto &sigmar2_abab = sigvec_blks_["sigmar2_abab"];
    auto &sigmar2_bbbb = sigvec_blks_["sigmar2_bbbb"];
    auto &sigmal2_aaaa = sigvec_blks_["sigmal2_aaaa"];
    auto &sigmal2_abab = sigvec_blks_["sigmal2_abab"];
    auto &sigmal2_bbbb = sigvec_blks_["sigmal2_bbbb"];


    // flops: o1v1L1  = o3v1L1 o2v1L1 o1v1L1 o1v3L1 o1v1L1 o1v3L1 o1v1L1 o3v1L1 o1v1L1 o3v1L1 o1v1L1 o1v3L1 o1v1L1 o1v3L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o2v1L1 o1v1L1 o1v2L1 o1v1L1 o3v2L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o1v2L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o2v1L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o2v1L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o1v1L1 o3v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    tmps_["13_bb_Lov"]("L,i,a")  = eri["abba_oovo"]("l,i,a,k") * tmps_["8_aa_Loo"]("L,k,l");
    tmps_["13_bb_Lov"]("L,i,a") -= f["bb_ov"]("k,a") * tmps_["7_bb_Loo"]("L,i,k");
    tmps_["13_bb_Lov"]("L,i,a") -= eri["bbbb_vovv"]("c,i,a,d") * tmps_["3_bb_Lvv"]("L,c,d");
    tmps_["13_bb_Lov"]("L,i,a") += eri["abab_vovv"]("c,i,d,a") * tmps_["2_aa_Lvv"]("L,c,d");
    tmps_["13_bb_Lov"]("L,i,a") -= 0.50 * eri["bbbb_oovo"]("i,l,a,k") * tmps_["9_bb_Loo"]("L,k,l");
    tmps_["13_bb_Lov"]("L,i,a") -= eri["bbbb_oovo"]("i,l,a,k") * tmps_["7_bb_Loo"]("L,k,l");
    tmps_["13_bb_Lov"]("L,i,a") += 0.50 * eri["abab_vovv"]("c,i,d,a") * tmps_["1_aa_Lvv"]("L,c,d");
    tmps_["13_bb_Lov"]("L,i,a") -= 0.50 * eri["bbbb_vovv"]("c,i,a,d") * tmps_["4_bb_Lvv"]("L,c,d");
    tmps_["13_bb_Lov"]("L,i,a") += 0.50 * f["bb_ov"]("k,a") * tmps_["12_bb_Loo"]("L,i,k");
    tmps_["13_bb_Lov"]("L,i,a") -= l1["aa"]("L,j,b") * reused_["7_aabb_voov"]("b,j,i,a");
    tmps_["13_bb_Lov"]("L,i,a") -= reused_["45_baab_vovv"]("a,j,c,b") * l2["abab"]("L,j,i,c,b");
    tmps_["13_bb_Lov"]("L,i,a") += l2["abab"]("L,j,i,b,a") * reused_["82_aa_vo"]("b,j");
    tmps_["13_bb_Lov"]("L,i,a") -= l2["abab"]("L,j,k,b,a") * reused_["93_abab_vooo"]("b,i,j,k");
    tmps_["13_bb_Lov"]("L,i,a") -= 0.50 * l1["bb"]("L,j,a") * reused_["87_bb_oo"]("i,j");
    tmps_["13_bb_Lov"]("L,i,a") -= 0.50 * reused_["55_bb_vv"]("a,b") * l1["bb"]("L,i,b");
    tmps_["13_bb_Lov"]("L,i,a") += l2["abab"]("L,j,k,b,a") * reused_["33_aabb_vooo"]("b,j,i,k");
    tmps_["13_bb_Lov"]("L,i,a") -= 0.25 * scalars_["9"] * l1["bb"]("L,i,a");
    tmps_["13_bb_Lov"]("L,i,a") -= f["bb_vo"]("b,j") * l2["bbbb"]("L,i,j,b,a");
    tmps_["13_bb_Lov"]("L,i,a") -= 0.50 * eri["bbbb_vvvo"]("b,c,a,j") * l2["bbbb"]("L,i,j,c,b");
    tmps_["13_bb_Lov"]("L,i,a") += l2["abab"]("L,j,k,b,a") * reused_["38_bbaa_oovo"]("i,k,b,j");
    tmps_["13_bb_Lov"]("L,i,a") -= reused_["1_bbbb_vvvo"]("c,a,b,j") * l2["bbbb"]("L,i,j,c,b");
    tmps_["13_bb_Lov"]("L,i,a") -= l2["bbbb"]("L,j,k,b,a") * reused_["66_bbbb_vooo"]("b,j,i,k");
    tmps_["13_bb_Lov"]("L,i,a") -= reused_["59_bb_vv"]("a,b") * l1["bb"]("L,i,b");
    tmps_["13_bb_Lov"]("L,i,a") -= reused_["14_bbbb_vovv"]("b,j,c,a") * l2["bbbb"]("L,i,j,c,b");
    tmps_["13_bb_Lov"]("L,i,a") -= l2["bbbb"]("L,j,k,b,a") * reused_["35_bbbb_oovo"]("i,k,b,j");
    tmps_["13_bb_Lov"]("L,i,a") -= l1["bb"]("L,j,a") * reused_["50_bb_oo"]("j,i");
    tmps_["13_bb_Lov"]("L,i,a") += reused_["18_bbaa_vvvo"]("c,a,b,j") * l2["abab"]("L,j,i,b,c");
    tmps_["13_bb_Lov"]("L,i,a") -= eri["abba_vovo"]("b,i,a,j") * l1["aa"]("L,j,b");
    tmps_["13_bb_Lov"]("L,i,a") += reused_["16_aabb_vovv"]("b,j,c,a") * l2["abab"]("L,j,i,b,c");
    tmps_["13_bb_Lov"]("L,i,a") -= l2["abab"]("L,k,j,b,a") * reused_["31_baab_oovo"]("i,k,b,j");
    tmps_["13_bb_Lov"]("L,i,a") += l1["bb"]("L,j,b") * reused_["4_bbbb_voov"]("b,j,i,a");
    tmps_["13_bb_Lov"]("L,i,a") += l1["aa"]("L,j,b") * reused_["46_aabb_voov"]("b,j,i,a");
    tmps_["13_bb_Lov"]("L,i,a") -= 0.25 * reused_["39_bbbb_vovv"]("a,j,c,b") * l2["bbbb"]("L,i,j,c,b");
    tmps_["13_bb_Lov"]("L,i,a") -= eri["bbbb_vovo"]("b,i,a,j") * l1["bb"]("L,j,b");
    tmps_["13_bb_Lov"]("L,i,a") += l2["bbbb"]("L,i,j,b,a") * reused_["78_bb_vo"]("b,j");
    tmps_["13_bb_Lov"]("L,i,a") += f["bb_vv"]("b,a") * l1["bb"]("L,i,b");
    tmps_["13_bb_Lov"]("L,i,a") -= l1["bb"]("L,j,b") * reused_["41_bbbb_voov"]("b,j,i,a");
    tmps_["13_bb_Lov"]("L,i,a") += f["aa_vo"]("b,j") * l2["abab"]("L,j,i,b,a");
    tmps_["13_bb_Lov"]("L,i,a") -= eri["abba_vvvo"]("c,b,a,j") * l2["abab"]("L,j,i,c,b");
    tmps_["13_bb_Lov"]("L,i,a") -= reused_["17_abba_vvvo"]("c,a,b,j") * l2["abab"]("L,j,i,c,b");
    tmps_["13_bb_Lov"]("L,i,a") -= 0.50 * eri["bbbb_vooo"]("b,i,j,k") * l2["bbbb"]("L,j,k,b,a");
    tmps_["13_bb_Lov"]("L,i,a") -= f["bb_oo"]("i,j") * l1["bb"]("L,j,a");
    tmps_["13_bb_Lov"]("L,i,a") -= 0.25 * l2["bbbb"]("L,j,k,b,a") * reused_["72_bbbb_vooo"]("b,i,j,k");
    tmps_["13_bb_Lov"]("L,i,a") -= eri["abab_vooo"]("b,i,j,k") * l2["abab"]("L,j,k,b,a");
    tmps_["13_bb_Lov"]("L,i,a") += 0.50 * eri["abba_oovo"]("l,i,a,k") * tmps_["6_aa_Loo"]("L,k,l");
    tmps_["12_bb_Loo"].~TArrayD();
    tmps_["9_bb_Loo"].~TArrayD();
    tmps_["8_aa_Loo"].~TArrayD();
    tmps_["7_bb_Loo"].~TArrayD();
    tmps_["6_aa_Loo"].~TArrayD();
    tmps_["4_bb_Lvv"].~TArrayD();
    tmps_["3_bb_Lvv"].~TArrayD();
    tmps_["2_aa_Lvv"].~TArrayD();
    tmps_["1_aa_Lvv"].~TArrayD();

    // sigmal1_bb  = -0.50 <l,i||k,a>_abab t2_abab(c,b,l,j) l2_abab(k,j,c,b)
    //            += -0.50 <l,i||k,a>_abab t2_abab(b,c,l,j) l2_abab(k,j,b,c)
    //            += -0.50 f_bb(k,a) t2_abab(c,b,j,k) l2_abab(j,i,c,b)
    //            += -0.50 f_bb(k,a) t2_abab(b,c,j,k) l2_abab(j,i,b,c)
    //            += -0.50 <i,c||d,a>_bbbb t2_abab(b,d,j,k) l2_abab(j,k,b,c)
    //            += -0.50 <i,c||d,a>_bbbb t2_abab(b,d,k,j) l2_abab(k,j,b,c)
    //            += +0.50 <c,i||d,a>_abab t2_abab(d,b,j,k) l2_abab(j,k,c,b)
    //            += +0.50 <c,i||d,a>_abab t2_abab(d,b,k,j) l2_abab(k,j,c,b)
    //            += -0.50 <i,l||a,k>_bbbb t2_bbbb(c,b,j,l) l2_bbbb(j,k,c,b)
    //            += -0.50 <i,l||a,k>_bbbb t2_abab(c,b,j,l) l2_abab(j,k,c,b)
    //            += -0.50 <i,l||a,k>_bbbb t2_abab(b,c,j,l) l2_abab(j,k,b,c)
    //            += +0.50 <c,i||d,a>_abab t2_aaaa(d,b,j,k) l2_aaaa(j,k,c,b)
    //            += -0.50 <i,c||d,a>_bbbb t2_bbbb(d,b,j,k) l2_bbbb(j,k,c,b)
    //            += +0.50 f_bb(k,a) t2_bbbb(c,b,j,k) l2_bbbb(i,j,c,b)
    //            += -1.00 <k,i||c,a>_abab t2_aaaa(c,b,j,k) l1_aa(j,b)
    //            += +0.25 <l,k||j,a>_abab t2_abab(c,b,l,k) l2_abab(j,i,c,b)
    //            += +0.25 <k,l||j,a>_abab t2_abab(c,b,k,l) l2_abab(j,i,c,b)
    //            += +0.25 <l,k||j,a>_abab t2_abab(b,c,l,k) l2_abab(j,i,b,c)
    //            += +0.25 <k,l||j,a>_abab t2_abab(b,c,k,l) l2_abab(j,i,b,c)
    //            += +1.00 f_bb(k,c) t2_abab(b,c,j,k) l2_abab(j,i,b,a)
    //            += -0.50 <l,k||j,c>_abab t2_abab(b,c,l,k) l2_abab(j,i,b,a)
    //            += -0.50 <k,l||j,c>_abab t2_abab(b,c,k,l) l2_abab(j,i,b,a)
    //            += +0.50 <b,k||c,d>_abab t2_abab(c,d,j,k) l2_abab(j,i,b,a)
    //            += +0.50 <b,k||d,c>_abab t2_abab(d,c,j,k) l2_abab(j,i,b,a)
    //            += -1.00 f_aa(k,c) t2_aaaa(c,b,j,k) l2_abab(j,i,b,a)
    //            += -0.50 <l,k||c,j>_aaaa t2_aaaa(c,b,l,k) l2_abab(j,i,b,a)
    //            += -0.50 <k,b||c,d>_aaaa t2_aaaa(c,d,j,k) l2_abab(j,i,b,a)
    //            += -0.25 <b,i||c,d>_abab t2_abab(c,d,j,k) l2_abab(j,k,b,a)
    //            += -0.25 <b,i||c,d>_abab t2_abab(c,d,k,j) l2_abab(k,j,b,a)
    //            += -0.25 <b,i||d,c>_abab t2_abab(d,c,j,k) l2_abab(j,k,b,a)
    //            += -0.25 <b,i||d,c>_abab t2_abab(d,c,k,j) l2_abab(k,j,b,a)
    //            += -0.50 f_bb(i,c) t2_abab(b,c,j,k) l2_abab(j,k,b,a)
    //            += -0.50 f_bb(i,c) t2_abab(b,c,k,j) l2_abab(k,j,b,a)
    //            += -0.50 <i,k||b,c>_bbbb t2_bbbb(b,c,j,k) l1_bb(j,a)
    //            += -0.50 <k,j||c,a>_bbbb t2_bbbb(c,b,k,j) l1_bb(i,b)
    //            += +1.00 <i,l||c,k>_bbbb t2_abab(b,c,j,l) l2_abab(j,k,b,a)
    //            += +0.25 <k,j||b,c>_bbbb t2_bbbb(b,c,k,j) l1_bb(i,a)
    //            += +1.00 f_aa(j,j) l1_bb(i,a)
    //            += +1.00 f_bb(j,j) l1_bb(i,a)
    //            += -0.50 <k,j||k,j>_aaaa l1_bb(i,a)
    //            += -0.50 <k,j||k,j>_bbbb l1_bb(i,a)
    //            += +0.25 <k,j||b,c>_aaaa t2_aaaa(b,c,k,j) l1_bb(i,a)
    //            += -0.50 <k,j||k,j>_abab l1_bb(i,a)
    //            += -0.50 <j,k||j,k>_abab l1_bb(i,a)
    //            += +0.25 <k,j||b,c>_abab t2_abab(b,c,k,j) l1_bb(i,a)
    //            += +0.25 <j,k||b,c>_abab t2_abab(b,c,j,k) l1_bb(i,a)
    //            += +0.25 <k,j||c,b>_abab t2_abab(c,b,k,j) l1_bb(i,a)
    //            += +0.25 <j,k||c,b>_abab t2_abab(c,b,j,k) l1_bb(i,a)
    //            += -1.00 f_bb(b,j) l2_bbbb(i,j,b,a)
    //            += +0.50 <c,b||a,j>_bbbb l2_bbbb(i,j,c,b)
    //            += +1.00 <l,i||c,k>_abab t2_aaaa(c,b,j,l) l2_abab(j,k,b,a)
    //            += +1.00 <k,c||d,a>_abab t2_abab(d,b,k,j) l2_bbbb(i,j,c,b)
    //            += -1.00 <i,l||c,k>_bbbb t2_bbbb(c,b,j,l) l2_bbbb(j,k,b,a)
    //            += -0.50 <k,j||c,a>_abab t2_abab(c,b,k,j) l1_bb(i,b)
    //            += -0.50 <j,k||c,a>_abab t2_abab(c,b,j,k) l1_bb(i,b)
    //            += -1.00 <k,c||d,a>_bbbb t2_bbbb(d,b,j,k) l2_bbbb(i,j,c,b)
    //            += -1.00 <l,i||c,k>_abab t2_abab(c,b,l,j) l2_bbbb(j,k,b,a)
    //            += -0.50 <k,i||b,c>_abab t2_abab(b,c,k,j) l1_bb(j,a)
    //            += -0.50 <k,i||c,b>_abab t2_abab(c,b,k,j) l1_bb(j,a)
    //            += -1.00 <k,c||d,a>_abab t2_aaaa(d,b,j,k) l2_abab(j,i,b,c)
    //            += +1.00 <b,i||j,a>_abab l1_aa(j,b)
    //            += +1.00 <k,c||d,a>_bbbb t2_abab(b,d,j,k) l2_abab(j,i,b,c)
    //            += +1.00 <l,i||k,c>_abab t2_abab(b,c,l,j) l2_abab(k,j,b,a)
    //            += +1.00 <k,i||c,a>_abab t2_abab(c,b,k,j) l1_bb(j,b)
    //            += -1.00 <i,k||c,a>_bbbb t2_abab(b,c,j,k) l1_aa(j,b)
    //            += +0.25 <l,k||a,j>_bbbb t2_bbbb(c,b,l,k) l2_bbbb(i,j,c,b)
    //            += +1.00 <i,b||a,j>_bbbb l1_bb(j,b)
    //            += +1.00 f_bb(k,c) t2_bbbb(c,b,j,k) l2_bbbb(i,j,b,a)
    //            += +0.50 <l,k||c,j>_bbbb t2_bbbb(c,b,l,k) l2_bbbb(i,j,b,a)
    //            += +0.50 <l,k||c,j>_abab t2_abab(c,b,l,k) l2_bbbb(i,j,b,a)
    //            += +0.50 <k,l||c,j>_abab t2_abab(c,b,k,l) l2_bbbb(i,j,b,a)
    //            += -1.00 f_aa(k,c) t2_abab(c,b,k,j) l2_bbbb(i,j,b,a)
    //            += -0.50 <k,b||c,d>_abab t2_abab(c,d,k,j) l2_bbbb(i,j,b,a)
    //            += -0.50 <k,b||d,c>_abab t2_abab(d,c,k,j) l2_bbbb(i,j,b,a)
    //            += +0.50 <k,b||c,d>_bbbb t2_bbbb(c,d,j,k) l2_bbbb(i,j,b,a)
    //            += +1.00 f_bb(b,a) l1_bb(i,b)
    //            += +1.00 <i,k||c,a>_bbbb t2_bbbb(c,b,j,k) l1_bb(j,b)
    //            += +1.00 f_aa(b,j) l2_abab(j,i,b,a)
    //            += +0.50 <c,b||j,a>_abab l2_abab(j,i,c,b)
    //            += +0.50 <b,c||j,a>_abab l2_abab(j,i,b,c)
    //            += -1.00 <c,k||d,a>_abab t2_abab(d,b,j,k) l2_abab(j,i,c,b)
    //            += +0.50 <i,b||j,k>_bbbb l2_bbbb(j,k,b,a)
    //            += -1.00 f_bb(i,j) l1_bb(j,a)
    //            += +0.25 <i,b||c,d>_bbbb t2_bbbb(c,d,j,k) l2_bbbb(j,k,b,a)
    //            += +0.50 f_bb(i,c) t2_bbbb(c,b,j,k) l2_bbbb(j,k,b,a)
    //            += -0.50 <b,i||j,k>_abab l2_abab(j,k,b,a)
    //            += -0.50 <b,i||k,j>_abab l2_abab(k,j,b,a)
    //            += -0.50 <l,i||k,a>_abab t2_aaaa(c,b,j,l) l2_aaaa(j,k,c,b)
    sigmal1_bb("L,a,i")  = tmps_["13_bb_Lov"]("L,i,a");
    tmps_["13_bb_Lov"].~TArrayD();

    // flops: o2v0L1  = o3v2L1 o3v2L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
    tmps_["14_aa_Loo"]("L,i,l")  = l2["abab"]("L,i,k,d,c") * t2["abab"]("d,c,l,k");
    tmps_["14_aa_Loo"]("L,i,l") -= 0.50 * l2["aaaa"]("L,i,k,d,c") * t2["aaaa"]("d,c,k,l");

    // sigmal2_aaaa += +0.50 P(i,j) <j,l||a,b>_aaaa t2_abab(d,c,l,k) l2_abab(i,k,d,c)
    //              += +0.50 P(i,j) <j,l||a,b>_aaaa t2_abab(c,d,l,k) l2_abab(i,k,c,d)
    //              += -0.50 P(i,j) <j,l||a,b>_aaaa t2_aaaa(d,c,k,l) l2_aaaa(i,k,d,c)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j")  = eri["aaaa_oovv"]("j,l,a,b") * tmps_["14_aa_Loo"]("L,i,l");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,a,b,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmal2_abab += -0.50 <l,j||a,b>_abab t2_abab(d,c,l,k) l2_abab(i,k,d,c)
    //              += -0.50 <l,j||a,b>_abab t2_abab(c,d,l,k) l2_abab(i,k,c,d)
    //              += +0.50 <l,j||a,b>_abab t2_aaaa(d,c,k,l) l2_aaaa(i,k,d,c)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= eri["abab_oovv"]("l,j,a,b") * tmps_["14_aa_Loo"]("L,i,l");
    tmps_["14_aa_Loo"].~TArrayD();

    // flops: o0v2L1  = o2v3L1 o2v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["15_aa_Lvv"]("L,a,d")  = l2["abab"]("L,k,l,a,c") * t2["abab"]("d,c,k,l");
    tmps_["15_aa_Lvv"]("L,a,d") -= 0.50 * l2["aaaa"]("L,k,l,c,a") * t2["aaaa"]("d,c,k,l");

    // sigmal2_aaaa += +0.50 P(a,b) <i,j||d,a>_aaaa t2_abab(d,c,k,l) l2_abab(k,l,b,c)
    //              += +0.50 P(a,b) <i,j||d,a>_aaaa t2_abab(d,c,l,k) l2_abab(l,k,b,c)
    //              += -0.50 P(a,b) <i,j||d,a>_aaaa t2_aaaa(d,c,k,l) l2_aaaa(k,l,c,b)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j")  = eri["aaaa_oovv"]("i,j,a,d") * tmps_["15_aa_Lvv"]("L,b,d");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,b,a,i,j");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmal2_abab += -0.50 <i,j||d,b>_abab t2_abab(d,c,k,l) l2_abab(k,l,a,c)
    //              += -0.50 <i,j||d,b>_abab t2_abab(d,c,l,k) l2_abab(l,k,a,c)
    //              += +0.50 <i,j||d,b>_abab t2_aaaa(d,c,k,l) l2_aaaa(k,l,c,a)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= eri["abab_oovv"]("i,j,d,b") * tmps_["15_aa_Lvv"]("L,a,d");
    tmps_["15_aa_Lvv"].~TArrayD();

    // flops: o0v2L1  = o2v3L1 o2v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["16_bb_Lvv"]("L,b,d")  = l2["bbbb"]("L,k,l,c,b") * t2["bbbb"]("d,c,k,l");
    tmps_["16_bb_Lvv"]("L,b,d") -= 2.00 * l2["abab"]("L,k,l,c,b") * t2["abab"]("c,d,k,l");

    // sigmal2_abab += +0.50 <i,j||a,d>_abab t2_bbbb(d,c,k,l) l2_bbbb(k,l,c,b)
    //              += -0.50 <i,j||a,d>_abab t2_abab(c,d,k,l) l2_abab(k,l,c,b)
    //              += -0.50 <i,j||a,d>_abab t2_abab(c,d,l,k) l2_abab(l,k,c,b)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") += 0.50 * eri["abab_oovv"]("i,j,a,d") * tmps_["16_bb_Lvv"]("L,b,d");

    // sigmal2_bbbb += -0.50 P(a,b) <i,j||d,a>_bbbb t2_bbbb(d,c,k,l) l2_bbbb(k,l,c,b)
    //              += +0.50 P(a,b) <i,j||d,a>_bbbb t2_abab(c,d,k,l) l2_abab(k,l,c,b)
    //              += +0.50 P(a,b) <i,j||d,a>_bbbb t2_abab(c,d,l,k) l2_abab(l,k,c,b)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j")  = 0.50 * eri["bbbb_oovv"]("i,j,a,d") * tmps_["16_bb_Lvv"]("L,b,d");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,b,a,i,j");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();
    tmps_["16_bb_Lvv"].~TArrayD();

    // flops: o2v0L1  = o3v2L1 o3v2L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
    tmps_["17_bb_Loo"]("L,j,l")  = l2["bbbb"]("L,j,k,d,c") * t2["bbbb"]("d,c,k,l");
    tmps_["17_bb_Loo"]("L,j,l") -= 2.00 * l2["abab"]("L,k,j,d,c") * t2["abab"]("d,c,k,l");

    // sigmal2_abab += +0.50 <i,l||a,b>_abab t2_bbbb(d,c,k,l) l2_bbbb(j,k,d,c)
    //              += -0.50 <i,l||a,b>_abab t2_abab(d,c,k,l) l2_abab(k,j,d,c)
    //              += -0.50 <i,l||a,b>_abab t2_abab(c,d,k,l) l2_abab(k,j,c,d)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") += 0.50 * eri["abab_oovv"]("i,l,a,b") * tmps_["17_bb_Loo"]("L,j,l");

    // sigmal2_bbbb += -0.50 P(i,j) <j,l||a,b>_bbbb t2_bbbb(d,c,k,l) l2_bbbb(i,k,d,c)
    //              += +0.50 P(i,j) <j,l||a,b>_bbbb t2_abab(d,c,k,l) l2_abab(k,i,d,c)
    //              += +0.50 P(i,j) <j,l||a,b>_bbbb t2_abab(c,d,k,l) l2_abab(k,i,c,d)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j")  = 0.50 * eri["bbbb_oovv"]("j,l,a,b") * tmps_["17_bb_Loo"]("L,i,l");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,a,b,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();
    tmps_["17_bb_Loo"].~TArrayD();

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["18_bb_Lov"]("R,j,b")  = eri["abab_oovv"]("k,j,c,b") * r1["aa"]("R,c,k");

    // flops: o1v1L1  = o2v2L1 o2v1L1 o1v1L1 o1v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o1v1L1 o1v1L1 o2v3L1 o1v1L1 o1v2L1 o1v1L1 o2v2L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    tmps_["19_aa_Lvo"]("R,a,i")  = f["aa_ov"]("j,b") * r2["aaaa"]("R,b,a,i,j");
    tmps_["19_aa_Lvo"]("R,a,i") += r1["aa"]("R,a,k") * reused_["92_aa_oo"]("i,k");
    tmps_["19_aa_Lvo"]("R,a,i") += r1["aa"]("R,c,i") * reused_["61_aa_vv"]("a,c");
    tmps_["19_aa_Lvo"]("R,a,i") += eri["abba_vovo"]("a,j,b,i") * r1["bb"]("R,b,j");
    tmps_["19_aa_Lvo"]("R,a,i") += r1["aa"]("R,c,k") * reused_["12_aaaa_voov"]("a,i,k,c");
    tmps_["19_aa_Lvo"]("R,a,i") -= f["aa_vv"]("a,b") * r1["aa"]("R,b,i");
    tmps_["19_aa_Lvo"]("R,a,i") -= 0.50 * eri["aaaa_oovo"]("j,k,b,i") * r2["aaaa"]("R,b,a,k,j");
    tmps_["19_aa_Lvo"]("R,a,i") += r1["bb"]("R,c,k") * reused_["7_aabb_voov"]("a,i,k,c");
    tmps_["19_aa_Lvo"]("R,a,i") -= r1["bb"]("R,c,k") * reused_["47_aabb_voov"]("a,i,k,c");
    tmps_["19_aa_Lvo"]("R,a,i") += f["aa_oo"]("j,i") * r1["aa"]("R,a,j");
    tmps_["19_aa_Lvo"]("R,a,i") += eri["aaaa_vovo"]("a,j,b,i") * r1["aa"]("R,b,j");
    tmps_["19_aa_Lvo"]("R,a,i") -= 0.50 * eri["aaaa_vovv"]("a,j,b,c") * r2["aaaa"]("R,b,c,i,j");
    tmps_["19_aa_Lvo"]("R,a,i") -= 0.50 * r1["aa"]("R,a,k") * reused_["53_aa_oo"]("i,k");
    tmps_["19_aa_Lvo"]("R,a,i") -= f["bb_ov"]("j,b") * r2["abab"]("R,a,b,i,j");
    tmps_["19_aa_Lvo"]("R,a,i") -= eri["abba_oovo"]("k,j,b,i") * r2["abab"]("R,a,b,k,j");
    tmps_["19_aa_Lvo"]("R,a,i") += 0.25 * scalars_["9"] * r1["aa"]("R,a,i");
    tmps_["19_aa_Lvo"]("R,a,i") -= eri["abab_vovv"]("a,j,b,c") * r2["abab"]("R,b,c,i,j");
    tmps_["19_aa_Lvo"]("R,a,i") -= 0.50 * r1["aa"]("R,c,i") * reused_["63_aa_vv"]("a,c");
    tmps_["19_aa_Lvo"]("R,a,i") -= t2["abab"]("a,b,i,j") * tmps_["18_bb_Lov"]("R,j,b");
    tmps_["18_bb_Lov"].~TArrayD();

    // sigmar1_aa  = -1.00 f_aa(j,b) r2_aaaa(b,a,i,j)
    //            += -0.50 <k,j||b,c>_abab r1_aa(a,k) t2_abab(b,c,i,j)
    //            += -0.50 <k,j||c,b>_abab r1_aa(a,k) t2_abab(c,b,i,j)
    //            += -0.50 <k,j||c,b>_abab r1_aa(c,i) t2_abab(a,b,k,j)
    //            += -0.50 <j,k||c,b>_abab r1_aa(c,i) t2_abab(a,b,j,k)
    //            += +1.00 <a,j||i,b>_abab r1_bb(b,j)
    //            += +1.00 <k,j||b,c>_aaaa r1_aa(c,k) t2_aaaa(b,a,i,j)
    //            += +1.00 f_aa(a,b) r1_aa(b,i)
    //            += -0.50 <k,j||b,i>_aaaa r2_aaaa(b,a,k,j)
    //            += -1.00 <j,k||b,c>_abab r1_bb(c,k) t2_aaaa(b,a,i,j)
    //            += -1.00 <k,j||b,c>_bbbb r1_bb(c,k) t2_abab(a,b,i,j)
    //            += -1.00 f_aa(j,i) r1_aa(a,j)
    //            += +1.00 <j,a||b,i>_aaaa r1_aa(b,j)
    //            += -0.50 <j,a||b,c>_aaaa r2_aaaa(b,c,i,j)
    //            += -0.50 <k,j||b,c>_aaaa r1_aa(a,k) t2_aaaa(b,c,i,j)
    //            += +1.00 f_bb(j,b) r2_abab(a,b,i,j)
    //            += -0.50 <k,j||i,b>_abab r2_abab(a,b,k,j)
    //            += -0.50 <j,k||i,b>_abab r2_abab(a,b,j,k)
    //            += +0.25 <k,j||b,c>_bbbb r1_aa(a,i) t2_bbbb(b,c,k,j)
    //            += +1.00 f_aa(j,j) r1_aa(a,i)
    //            += +1.00 f_bb(j,j) r1_aa(a,i)
    //            += -0.50 <k,j||k,j>_aaaa r1_aa(a,i)
    //            += -0.50 <k,j||k,j>_bbbb r1_aa(a,i)
    //            += +0.25 <k,j||b,c>_aaaa r1_aa(a,i) t2_aaaa(b,c,k,j)
    //            += -0.50 <k,j||k,j>_abab r1_aa(a,i)
    //            += -0.50 <j,k||j,k>_abab r1_aa(a,i)
    //            += +0.25 <k,j||b,c>_abab r1_aa(a,i) t2_abab(b,c,k,j)
    //            += +0.25 <j,k||b,c>_abab r1_aa(a,i) t2_abab(b,c,j,k)
    //            += +0.25 <k,j||c,b>_abab r1_aa(a,i) t2_abab(c,b,k,j)
    //            += +0.25 <j,k||c,b>_abab r1_aa(a,i) t2_abab(c,b,j,k)
    //            += +0.50 <a,j||b,c>_abab r2_abab(b,c,i,j)
    //            += +0.50 <a,j||c,b>_abab r2_abab(c,b,i,j)
    //            += -0.50 <k,j||b,c>_aaaa r1_aa(c,i) t2_aaaa(b,a,k,j)
    //            += +1.00 <k,j||c,b>_abab r1_aa(c,k) t2_abab(a,b,i,j)
    sigmar1_aa("R,a,i")  = -1.00 * tmps_["19_aa_Lvo"]("R,a,i");
    tmps_["19_aa_Lvo"].~TArrayD();

    // flops: o1v1L1  = o3v2L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o3v2L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    tmps_["20_bb_Lov"]("R,i,a")  = r2["bbbb"]("R,b,a,k,j") * eri["bbbb_oovo"]("j,k,b,i");
    tmps_["20_bb_Lov"]("R,i,a") -= 2.00 * f["bb_ov"]("j,b") * r2["bbbb"]("R,b,a,i,j");
    tmps_["20_bb_Lov"]("R,i,a") += eri["bbbb_vovv"]("a,j,b,c") * r2["bbbb"]("R,b,c,i,j");
    tmps_["20_bb_Lov"]("R,i,a") -= 2.00 * eri["bbbb_vovo"]("a,j,b,i") * r1["bb"]("R,b,j");
    tmps_["20_bb_Lov"]("R,i,a") -= 2.00 * r1["bb"]("R,a,j") * f["bb_oo"]("j,i");
    tmps_["20_bb_Lov"]("R,i,a") += 2.00 * r1["bb"]("R,c,k") * reused_["4_bbbb_voov"]("a,i,k,c");
    tmps_["20_bb_Lov"]("R,i,a") -= 2.00 * r1["bb"]("R,c,k") * reused_["42_bbbb_voov"]("a,i,k,c");
    tmps_["20_bb_Lov"]("R,i,a") -= 2.00 * r1["bb"]("R,c,i") * reused_["59_bb_vv"]("c,a");
    tmps_["20_bb_Lov"]("R,i,a") -= 0.50 * scalars_["9"] * r1["bb"]("R,a,i");
    tmps_["20_bb_Lov"]("R,i,a") += 2.00 * f["aa_ov"]("j,b") * r2["abab"]("R,b,a,j,i");
    tmps_["20_bb_Lov"]("R,i,a") += r1["bb"]("R,c,i") * reused_["56_bb_vv"]("a,c");
    tmps_["20_bb_Lov"]("R,i,a") += 2.00 * r1["aa"]("R,c,k") * reused_["6_bbaa_voov"]("a,i,k,c");
    tmps_["20_bb_Lov"]("R,i,a") -= 2.00 * eri["baab_vovv"]("a,j,b,c") * r2["abab"]("R,b,c,j,i");
    tmps_["20_bb_Lov"]("R,i,a") -= 2.00 * eri["baab_vovo"]("a,j,b,i") * r1["aa"]("R,b,j");
    tmps_["20_bb_Lov"]("R,i,a") -= 2.00 * r1["bb"]("R,a,k") * reused_["50_bb_oo"]("i,k");
    tmps_["20_bb_Lov"]("R,i,a") += r1["bb"]("R,a,k") * reused_["88_bb_oo"]("i,k");
    tmps_["20_bb_Lov"]("R,i,a") -= 2.00 * r1["aa"]("R,c,k") * reused_["15_bbaa_voov"]("a,i,k,c");
    tmps_["20_bb_Lov"]("R,i,a") += 2.00 * f["bb_vv"]("a,b") * r1["bb"]("R,b,i");
    tmps_["20_bb_Lov"]("R,i,a") -= 2.00 * eri["abab_oovo"]("k,j,b,i") * r2["abab"]("R,b,a,k,j");

    // sigmar1_bb  = -0.50 <k,j||b,i>_bbbb r2_bbbb(b,a,k,j)
    //            += -1.00 f_bb(j,b) r2_bbbb(b,a,i,j)
    //            += -0.50 <j,a||b,c>_bbbb r2_bbbb(b,c,i,j)
    //            += +1.00 <j,a||b,i>_bbbb r1_bb(b,j)
    //            += -1.00 f_bb(j,i) r1_bb(a,j)
    //            += +1.00 <j,k||b,c>_abab r1_bb(c,k) t2_abab(b,a,j,i)
    //            += +1.00 <k,j||b,c>_bbbb r1_bb(c,k) t2_bbbb(b,a,i,j)
    //            += -0.50 <k,j||b,c>_abab r1_bb(c,i) t2_abab(b,a,k,j)
    //            += -0.50 <j,k||b,c>_abab r1_bb(c,i) t2_abab(b,a,j,k)
    //            += +0.25 <k,j||b,c>_bbbb r1_bb(a,i) t2_bbbb(b,c,k,j)
    //            += +1.00 f_aa(j,j) r1_bb(a,i)
    //            += +1.00 f_bb(j,j) r1_bb(a,i)
    //            += -0.50 <k,j||k,j>_aaaa r1_bb(a,i)
    //            += -0.50 <k,j||k,j>_bbbb r1_bb(a,i)
    //            += +0.25 <k,j||b,c>_aaaa r1_bb(a,i) t2_aaaa(b,c,k,j)
    //            += -0.50 <k,j||k,j>_abab r1_bb(a,i)
    //            += -0.50 <j,k||j,k>_abab r1_bb(a,i)
    //            += +0.25 <k,j||b,c>_abab r1_bb(a,i) t2_abab(b,c,k,j)
    //            += +0.25 <j,k||b,c>_abab r1_bb(a,i) t2_abab(b,c,j,k)
    //            += +0.25 <k,j||c,b>_abab r1_bb(a,i) t2_abab(c,b,k,j)
    //            += +0.25 <j,k||c,b>_abab r1_bb(a,i) t2_abab(c,b,j,k)
    //            += +1.00 f_aa(j,b) r2_abab(b,a,j,i)
    //            += -0.50 <k,j||b,c>_bbbb r1_bb(c,i) t2_bbbb(b,a,k,j)
    //            += -1.00 <k,j||b,c>_aaaa r1_aa(c,k) t2_abab(b,a,j,i)
    //            += +0.50 <j,a||b,c>_abab r2_abab(b,c,j,i)
    //            += +0.50 <j,a||c,b>_abab r2_abab(c,b,j,i)
    //            += +1.00 <j,a||b,i>_abab r1_aa(b,j)
    //            += -0.50 <j,k||b,c>_abab r1_bb(a,k) t2_abab(b,c,j,i)
    //            += -0.50 <j,k||c,b>_abab r1_bb(a,k) t2_abab(c,b,j,i)
    //            += -0.50 <k,j||b,c>_bbbb r1_bb(a,k) t2_bbbb(b,c,i,j)
    //            += -1.00 <k,j||c,b>_abab r1_aa(c,k) t2_bbbb(b,a,i,j)
    //            += +1.00 f_bb(a,b) r1_bb(b,i)
    //            += -0.50 <k,j||b,i>_abab r2_abab(b,a,k,j)
    //            += -0.50 <j,k||b,i>_abab r2_abab(b,a,j,k)
    sigmar1_bb("R,a,i")  = 0.50 * tmps_["20_bb_Lov"]("R,i,a");
    tmps_["20_bb_Lov"].~TArrayD();

    // flops: o2v0L1  = o3v1L1 o3v1L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
    tmps_["21_aa_Loo"]("R,k,i")  = eri["abba_oovo"]("k,l,c,i") * r1["bb"]("R,c,l");
    tmps_["21_aa_Loo"]("R,k,i") += eri["aaaa_oovo"]("k,l,c,i") * r1["aa"]("R,c,l");

    // sigmar2_aaaa += -1.00 P(i,j) <k,l||j,c>_abab r1_bb(c,l) t2_aaaa(a,b,i,k)
    //              += -1.00 P(i,j) <l,k||c,j>_aaaa r1_aa(c,l) t2_aaaa(a,b,i,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = t2["aaaa"]("a,b,i,k") * tmps_["21_aa_Loo"]("R,k,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,a,b,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_abab += -1.00 <k,l||i,c>_abab r1_bb(c,l) t2_abab(a,b,k,j)
    //              += -1.00 <l,k||c,i>_aaaa r1_aa(c,l) t2_abab(a,b,k,j)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += t2["abab"]("a,b,k,j") * tmps_["21_aa_Loo"]("R,k,i");
    tmps_["21_aa_Loo"].~TArrayD();

    // flops: o0v2L1  = o2v3L1 o2v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["22_aa_Lvv"]("R,a,c")  = eri["abab_oovv"]("l,k,c,d") * r2["abab"]("R,a,d,l,k");
    tmps_["22_aa_Lvv"]("R,a,c") += 0.50 * eri["aaaa_oovv"]("k,l,c,d") * r2["aaaa"]("R,d,a,l,k");

    // sigmar2_aaaa += +0.50 P(a,b) <l,k||c,d>_abab r2_abab(b,d,l,k) t2_aaaa(c,a,i,j)
    //              += +0.50 P(a,b) <k,l||c,d>_abab r2_abab(b,d,k,l) t2_aaaa(c,a,i,j)
    //              += -0.50 P(a,b) <l,k||c,d>_aaaa r2_aaaa(d,b,l,k) t2_aaaa(c,a,i,j)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = t2["aaaa"]("c,a,i,j") * tmps_["22_aa_Lvv"]("R,b,c");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,b,a,i,j");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_abab += -0.50 <l,k||c,d>_abab r2_abab(a,d,l,k) t2_abab(c,b,i,j)
    //              += -0.50 <k,l||c,d>_abab r2_abab(a,d,k,l) t2_abab(c,b,i,j)
    //              += +0.50 <l,k||c,d>_aaaa r2_aaaa(d,a,l,k) t2_abab(c,b,i,j)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= t2["abab"]("c,b,i,j") * tmps_["22_aa_Lvv"]("R,a,c");
    tmps_["22_aa_Lvv"].~TArrayD();

    // flops: o2v0L1  = o2v1L1 o3v2L1 o2v0L1 o3v2L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1
    tmps_["23_aa_Loo"]("R,i,k")  = f["aa_ov"]("k,c") * r1["aa"]("R,c,i");
    tmps_["23_aa_Loo"]("R,i,k") += 0.50 * eri["aaaa_oovv"]("k,l,c,d") * r2["aaaa"]("R,c,d,i,l");
    tmps_["23_aa_Loo"]("R,i,k") += eri["abab_oovv"]("k,l,c,d") * r2["abab"]("R,c,d,i,l");

    // sigmar2_aaaa += +1.00 P(i,j) f_aa(k,c) r1_aa(c,i) t2_aaaa(a,b,j,k)
    //              += -0.50 P(i,j) <l,k||c,d>_aaaa r2_aaaa(c,d,i,l) t2_aaaa(a,b,j,k)
    //              += +0.50 P(i,j) <k,l||c,d>_abab r2_abab(c,d,i,l) t2_aaaa(a,b,j,k)
    //              += +0.50 P(i,j) <k,l||d,c>_abab r2_abab(d,c,i,l) t2_aaaa(a,b,j,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = t2["aaaa"]("a,b,j,k") * tmps_["23_aa_Loo"]("R,i,k");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,a,b,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_abab += -1.00 f_aa(k,c) r1_aa(c,i) t2_abab(a,b,k,j)
    //              += +0.50 <l,k||c,d>_aaaa r2_aaaa(c,d,i,l) t2_abab(a,b,k,j)
    //              += -0.50 <k,l||c,d>_abab r2_abab(c,d,i,l) t2_abab(a,b,k,j)
    //              += -0.50 <k,l||d,c>_abab r2_abab(d,c,i,l) t2_abab(a,b,k,j)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= t2["abab"]("a,b,k,j") * tmps_["23_aa_Loo"]("R,i,k");
    tmps_["23_aa_Loo"].~TArrayD();

    // flops: o0v2L1  = o1v3L1 o1v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["24_aa_Lvv"]("R,a,c")  = eri["abab_vovv"]("a,k,c,d") * r1["bb"]("R,d,k");
    tmps_["24_aa_Lvv"]("R,a,c") += eri["aaaa_vovv"]("a,k,c,d") * r1["aa"]("R,d,k");

    // sigmar2_aaaa += +1.00 P(a,b) <a,k||c,d>_abab r1_bb(d,k) t2_aaaa(c,b,i,j)
    //              += -1.00 P(a,b) <k,a||c,d>_aaaa r1_aa(d,k) t2_aaaa(c,b,i,j)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = t2["aaaa"]("c,b,i,j") * tmps_["24_aa_Lvv"]("R,a,c");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,b,a,i,j");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_abab += +1.00 <a,k||c,d>_abab r1_bb(d,k) t2_abab(c,b,i,j)
    //              += -1.00 <k,a||c,d>_aaaa r1_aa(d,k) t2_abab(c,b,i,j)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += t2["abab"]("c,b,i,j") * tmps_["24_aa_Lvv"]("R,a,c");
    tmps_["24_aa_Lvv"].~TArrayD();

    // flops: o2v0L1  = o3v1L1 o3v1L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
    tmps_["25_bb_Loo"]("R,k,j")  = eri["abab_oovo"]("l,k,c,j") * r1["aa"]("R,c,l");
    tmps_["25_bb_Loo"]("R,k,j") -= eri["bbbb_oovo"]("k,l,c,j") * r1["bb"]("R,c,l");

    // sigmar2_abab += -1.00 <l,k||c,j>_abab r1_aa(c,l) t2_abab(a,b,i,k)
    //              += -1.00 <l,k||c,j>_bbbb r1_bb(c,l) t2_abab(a,b,i,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= t2["abab"]("a,b,i,k") * tmps_["25_bb_Loo"]("R,k,j");

    // sigmar2_bbbb += -1.00 P(i,j) <l,k||c,j>_abab r1_aa(c,l) t2_bbbb(a,b,i,k)
    //              += -1.00 P(i,j) <l,k||c,j>_bbbb r1_bb(c,l) t2_bbbb(a,b,i,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = t2["bbbb"]("a,b,i,k") * tmps_["25_bb_Loo"]("R,k,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();
    tmps_["25_bb_Loo"].~TArrayD();

    // flops: o0v2L1  = o1v3L1 o1v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["26_bb_Lvv"]("R,b,c")  = eri["bbbb_vovv"]("b,k,c,d") * r1["bb"]("R,d,k");
    tmps_["26_bb_Lvv"]("R,b,c") -= eri["baab_vovv"]("b,k,d,c") * r1["aa"]("R,d,k");

    // sigmar2_abab += -1.00 <k,b||c,d>_bbbb r1_bb(d,k) t2_abab(a,c,i,j)
    //              += +1.00 <k,b||d,c>_abab r1_aa(d,k) t2_abab(a,c,i,j)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += t2["abab"]("a,c,i,j") * tmps_["26_bb_Lvv"]("R,b,c");

    // sigmar2_bbbb += -1.00 P(a,b) <k,a||c,d>_bbbb r1_bb(d,k) t2_bbbb(c,b,i,j)
    //              += +1.00 P(a,b) <k,a||d,c>_abab r1_aa(d,k) t2_bbbb(c,b,i,j)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = t2["bbbb"]("c,b,i,j") * tmps_["26_bb_Lvv"]("R,a,c");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,b,a,i,j");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();
    tmps_["26_bb_Lvv"].~TArrayD();

    // flops: o2v0L1  = o3v2L1 o2v1L1 o2v0L1 o3v2L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1
    tmps_["27_bb_Loo"]("R,j,k")  = eri["bbbb_oovv"]("k,l,c,d") * r2["bbbb"]("R,c,d,j,l");
    tmps_["27_bb_Loo"]("R,j,k") += 2.00 * f["bb_ov"]("k,c") * r1["bb"]("R,c,j");
    tmps_["27_bb_Loo"]("R,j,k") += 2.00 * eri["abab_oovv"]("l,k,c,d") * r2["abab"]("R,c,d,l,j");

    // sigmar2_abab += +0.50 <l,k||c,d>_bbbb r2_bbbb(c,d,j,l) t2_abab(a,b,i,k)
    //              += -1.00 f_bb(k,c) r1_bb(c,j) t2_abab(a,b,i,k)
    //              += -0.50 <l,k||c,d>_abab r2_abab(c,d,l,j) t2_abab(a,b,i,k)
    //              += -0.50 <l,k||d,c>_abab r2_abab(d,c,l,j) t2_abab(a,b,i,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= 0.50 * t2["abab"]("a,b,i,k") * tmps_["27_bb_Loo"]("R,j,k");

    // sigmar2_bbbb += -0.50 P(i,j) <l,k||c,d>_bbbb r2_bbbb(c,d,i,l) t2_bbbb(a,b,j,k)
    //              += +1.00 P(i,j) f_bb(k,c) r1_bb(c,i) t2_bbbb(a,b,j,k)
    //              += +0.50 P(i,j) <l,k||c,d>_abab r2_abab(c,d,l,i) t2_bbbb(a,b,j,k)
    //              += +0.50 P(i,j) <l,k||d,c>_abab r2_abab(d,c,l,i) t2_bbbb(a,b,j,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = 0.50 * t2["bbbb"]("a,b,j,k") * tmps_["27_bb_Loo"]("R,i,k");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();
    tmps_["27_bb_Loo"].~TArrayD();

    // flops: o0v2L1  = o2v3L1 o2v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["28_bb_Lvv"]("R,b,c")  = eri["abab_oovv"]("l,k,d,c") * r2["abab"]("R,d,b,l,k");
    tmps_["28_bb_Lvv"]("R,b,c") += 0.50 * eri["bbbb_oovv"]("k,l,c,d") * r2["bbbb"]("R,d,b,l,k");

    // sigmar2_abab += -0.50 <l,k||d,c>_abab r2_abab(d,b,l,k) t2_abab(a,c,i,j)
    //              += -0.50 <k,l||d,c>_abab r2_abab(d,b,k,l) t2_abab(a,c,i,j)
    //              += +0.50 <l,k||c,d>_bbbb r2_bbbb(d,b,l,k) t2_abab(a,c,i,j)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= t2["abab"]("a,c,i,j") * tmps_["28_bb_Lvv"]("R,b,c");

    // sigmar2_bbbb += +0.50 P(a,b) <l,k||d,c>_abab r2_abab(d,b,l,k) t2_bbbb(c,a,i,j)
    //              += +0.50 P(a,b) <k,l||d,c>_abab r2_abab(d,b,k,l) t2_bbbb(c,a,i,j)
    //              += -0.50 P(a,b) <l,k||c,d>_bbbb r2_bbbb(d,b,l,k) t2_bbbb(c,a,i,j)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = t2["bbbb"]("c,a,i,j") * tmps_["28_bb_Lvv"]("R,b,c");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,b,a,i,j");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();
    tmps_["28_bb_Lvv"].~TArrayD();

}