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

void hilbert::EOM_EE_CCSD::sigma_ee_00_3() {

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


    // sigmar2_abab += +1.00 <k,l||c,d>_abab r2_abab(a,d,i,l) t2_abab(c,b,k,j)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += r2["abab"]("R,a,d,i,l") * reused_["4_bbbb_voov"]("b,j,l,d");

    // sigmar2_abab += +1.00 <k,l||c,d>_abab r2_bbbb(d,b,j,l) t2_aaaa(c,a,i,k)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += r2["bbbb"]("R,d,b,j,l") * reused_["7_aabb_voov"]("a,i,l,d");

    // sigmar2_abab += +1.00 <l,k||c,d>_bbbb r2_abab(a,d,i,l) t2_bbbb(c,b,j,k)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= r2["abab"]("R,a,d,i,l") * reused_["42_bbbb_voov"]("b,j,l,d");

    // sigmar2_abab += +1.00 <l,k||c,d>_bbbb r2_bbbb(d,b,j,l) t2_abab(a,c,i,k)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= r2["bbbb"]("R,d,b,j,l") * reused_["47_aabb_voov"]("a,i,l,d");

    // sigmar2_abab += +1.00 <k,a||c,i>_aaaa r2_abab(c,b,k,j)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= eri["aaaa_vovo"]("a,k,c,i") * r2["abab"]("R,c,b,k,j");

    // sigmar2_abab += -1.00 <k,b||c,j>_abab r2_aaaa(c,a,i,k)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += r2["aaaa"]("R,c,a,i,k") * eri["baab_vovo"]("b,k,c,j");

    // sigmar2_abab += +1.00 <l,k||c,d>_aaaa r2_aaaa(d,a,i,l) t2_abab(c,b,k,j)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= r2["aaaa"]("R,d,a,i,l") * reused_["6_bbaa_voov"]("b,j,l,d");

    // sigmar2_abab += +1.00 <l,k||c,d>_aaaa r2_abab(d,b,l,j) t2_aaaa(c,a,i,k)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= r2["abab"]("R,d,b,l,j") * reused_["12_aaaa_voov"]("a,i,l,d");

    // sigmar2_abab += +1.00 <l,k||d,c>_abab r2_aaaa(d,a,i,l) t2_bbbb(c,b,j,k)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += r2["aaaa"]("R,d,a,i,l") * reused_["15_bbaa_voov"]("b,j,l,d");

    // sigmar2_abab += -1.00 <a,k||c,j>_abab r2_abab(c,b,i,k)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= eri["abab_vovo"]("a,k,c,j") * r2["abab"]("R,c,b,i,k");

    // sigmar2_abab += +0.50 <a,b||c,d>_abab r2_abab(c,d,i,j)
    //              += +0.50 <a,b||d,c>_abab r2_abab(d,c,i,j)
    // flops: o2v2L1 += o2v4L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += eri["abab_vvvv"]("a,b,c,d") * r2["abab"]("R,c,d,i,j");

    // sigmar2_abab += +0.25 <l,k||c,d>_abab r2_abab(c,d,i,j) t2_abab(a,b,l,k)
    //              += +0.25 <k,l||c,d>_abab r2_abab(c,d,i,j) t2_abab(a,b,k,l)
    //              += +0.25 <l,k||d,c>_abab r2_abab(d,c,i,j) t2_abab(a,b,l,k)
    //              += +0.25 <k,l||d,c>_abab r2_abab(d,c,i,j) t2_abab(a,b,k,l)
    // flops: o2v2L1 += o4v2L1 o4v2L1
    //  mems: o2v2L1 += o4v0L1 o2v2L1
    sigmar2_abab("R,a,b,i,j") += eri["abab_oovv"]("l,k,c,d") * r2["abab"]("R,c,d,i,j") * t2["abab"]("a,b,l,k");

    // sigmar2_abab += +1.00 <k,l||d,c>_abab r2_abab(d,b,i,l) t2_abab(a,c,k,j)
    // flops: o2v2L1 += o3v3L1 o3v3L1
    //  mems: o2v2L1 += o2v2L1 o2v2L1
    sigmar2_abab("R,a,b,i,j") += eri["abab_oovv"]("k,l,d,c") * r2["abab"]("R,d,b,i,l") * t2["abab"]("a,c,k,j");

    // sigmar2_abab += +1.00 <l,k||d,c>_abab r2_abab(d,b,l,j) t2_abab(a,c,i,k)
    // flops: o2v2L1 += o3v3L1 o3v3L1
    //  mems: o2v2L1 += o2v2L1 o2v2L1
    sigmar2_abab("R,a,b,i,j") += eri["abab_oovv"]("l,k,d,c") * r2["abab"]("R,d,b,l,j") * t2["abab"]("a,c,i,k");

    // sigmar2_bbbb += +0.25 <l,k||c,d>_bbbb r2_bbbb(a,b,i,j) t2_bbbb(c,d,l,k)
    //              += +1.00 f_aa(k,k) r2_bbbb(a,b,i,j)
    //              += +1.00 f_bb(k,k) r2_bbbb(a,b,i,j)
    //              += -0.50 <l,k||l,k>_aaaa r2_bbbb(a,b,i,j)
    //              += -0.50 <l,k||l,k>_bbbb r2_bbbb(a,b,i,j)
    //              += +0.25 <l,k||c,d>_aaaa r2_bbbb(a,b,i,j) t2_aaaa(c,d,l,k)
    //              += -0.50 <l,k||l,k>_abab r2_bbbb(a,b,i,j)
    //              += -0.50 <k,l||k,l>_abab r2_bbbb(a,b,i,j)
    //              += +0.25 <l,k||c,d>_abab r2_bbbb(a,b,i,j) t2_abab(c,d,l,k)
    //              += +0.25 <k,l||c,d>_abab r2_bbbb(a,b,i,j) t2_abab(c,d,k,l)
    //              += +0.25 <l,k||d,c>_abab r2_bbbb(a,b,i,j) t2_abab(d,c,l,k)
    //              += +0.25 <k,l||d,c>_abab r2_bbbb(a,b,i,j) t2_abab(d,c,k,l)
    // flops: o2v2L1 += o2v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_bbbb("R,a,b,i,j") -= 0.25 * scalars_["9"] * r2["bbbb"]("R,a,b,i,j");

    // sigmar2_bbbb += +1.00 P(i,j) P(a,b) f_bb(k,c) r1_bb(b,i) t2_bbbb(c,a,j,k)
    //              += +0.50 P(i,j) P(a,b) <l,k||c,j>_bbbb r1_bb(b,i) t2_bbbb(c,a,l,k)
    //              += +0.50 P(i,j) P(a,b) <l,k||c,j>_abab r1_bb(b,i) t2_abab(c,a,l,k)
    //              += +0.50 P(i,j) P(a,b) <k,l||c,j>_abab r1_bb(b,i) t2_abab(c,a,k,l)
    //              += -1.00 P(i,j) P(a,b) f_aa(k,c) r1_bb(b,i) t2_abab(c,a,k,j)
    //              += -0.50 P(i,j) P(a,b) <k,a||c,d>_abab r1_bb(b,i) t2_abab(c,d,k,j)
    //              += -0.50 P(i,j) P(a,b) <k,a||d,c>_abab r1_bb(b,i) t2_abab(d,c,k,j)
    //              += +0.50 P(i,j) P(a,b) <k,a||c,d>_bbbb r1_bb(b,i) t2_bbbb(c,d,j,k)
    // flops: o2v2L1 += o2v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = r1["bb"]("R,b,i") * reused_["78_bb_vo"]("a,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,j,i");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,b,a,i,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += -1.00 P(i,j) f_bb(k,j) r2_bbbb(a,b,i,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = r2["bbbb"]("R,a,b,i,k") * f["bb_oo"]("k,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += +1.00 P(a,b) <k,a||i,j>_bbbb r1_bb(b,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = eri["bbbb_vooo"]("a,k,i,j") * r1["bb"]("R,b,k");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,b,a,i,j");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += -0.50 P(i,j) <k,l||c,d>_abab r2_bbbb(a,b,i,l) t2_abab(c,d,k,j)
    //              += -0.50 P(i,j) <k,l||d,c>_abab r2_bbbb(a,b,i,l) t2_abab(d,c,k,j)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = r2["bbbb"]("R,a,b,i,l") * reused_["50_bb_oo"]("j,l");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += +0.50 P(a,b) <k,a||c,d>_bbbb r1_bb(b,k) t2_bbbb(c,d,i,j)
    //              += +1.00 P(a,b) f_bb(k,c) r1_bb(b,k) t2_bbbb(c,a,i,j)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = 0.50 * r1["bb"]("R,b,k") * reused_["72_bbbb_vooo"]("a,k,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,b,a,i,j");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += -1.00 P(i,j) P(a,b) <k,l||c,j>_abab r1_bb(b,l) t2_abab(c,a,k,i)
    //              += -1.00 P(i,j) P(a,b) <l,k||c,j>_bbbb r1_bb(b,l) t2_bbbb(c,a,i,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = r1["bb"]("R,b,l") * reused_["85_bbbb_oovo"]("l,j,a,i");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,j,i");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,b,a,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += -0.50 P(i,j) <l,k||c,d>_bbbb r2_bbbb(a,b,i,l) t2_bbbb(c,d,j,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = 0.50 * r2["bbbb"]("R,a,b,i,l") * reused_["88_bb_oo"]("j,l");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += +1.00 P(i,j) <a,b||c,j>_bbbb r1_bb(c,i)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = eri["bbbb_vvvo"]("a,b,c,j") * r1["bb"]("R,c,i");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += +1.00 P(a,b) f_bb(a,c) r2_bbbb(c,b,i,j)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = f["bb_vv"]("a,c") * r2["bbbb"]("R,c,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,b,a,i,j");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += -1.00 P(i,j) P(a,b) <k,a||c,d>_bbbb r1_bb(d,i) t2_bbbb(c,b,j,k)
    //              += +1.00 P(i,j) P(a,b) <k,a||c,d>_abab r1_bb(d,i) t2_abab(c,b,k,j)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = r1["bb"]("R,d,i") * reused_["3_bbbb_vvvo"]("a,d,b,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,j,i");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,b,a,i,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += +0.50 P(i,j) <l,k||c,j>_bbbb r1_bb(c,i) t2_bbbb(a,b,l,k)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = 0.50 * r1["bb"]("R,c,i") * reused_["39_bbbb_vovv"]("c,j,a,b");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += -0.50 P(a,b) <l,k||c,d>_bbbb r2_bbbb(d,b,i,j) t2_bbbb(c,a,l,k)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = 0.50 * r2["bbbb"]("R,d,b,i,j") * reused_["56_bb_vv"]("a,d");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,b,a,i,j");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += -0.50 P(a,b) <l,k||c,d>_abab r2_bbbb(d,b,i,j) t2_abab(c,a,l,k)
    //              += -0.50 P(a,b) <k,l||c,d>_abab r2_bbbb(d,b,i,j) t2_abab(c,a,k,l)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = r2["bbbb"]("R,d,b,i,j") * reused_["59_bb_vv"]("d,a");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,b,a,i,j");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += +0.50 <l,k||i,j>_bbbb r2_bbbb(a,b,l,k)
    // flops: o2v2L1 += o4v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_bbbb("R,a,b,i,j") -= 0.50 * r2["bbbb"]("R,a,b,l,k") * eri["bbbb_oooo"]("k,l,i,j");

    // sigmar2_bbbb += +0.25 <l,k||c,d>_bbbb r2_bbbb(a,b,l,k) t2_bbbb(c,d,i,j)
    // flops: o2v2L1 += o4v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_bbbb("R,a,b,i,j") -= 0.25 * r2["bbbb"]("R,a,b,l,k") * reused_["68_bbbb_oooo"]("i,j,k,l");

    // sigmar2_bbbb += +1.00 P(i,j) P(a,b) <k,a||c,j>_bbbb r2_bbbb(c,b,i,k)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = eri["bbbb_vovo"]("a,k,c,j") * r2["bbbb"]("R,c,b,i,k");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,j,i");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,b,a,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += +1.00 P(i,j) P(a,b) <k,l||c,d>_abab r2_bbbb(d,b,i,l) t2_abab(c,a,k,j)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = r2["bbbb"]("R,d,b,i,l") * reused_["4_bbbb_voov"]("a,j,l,d");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,j,i");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,b,a,i,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += +1.00 P(i,j) P(a,b) <l,k||c,d>_bbbb r2_bbbb(d,b,i,l) t2_bbbb(c,a,j,k)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = r2["bbbb"]("R,d,b,i,l") * reused_["42_bbbb_voov"]("a,j,l,d");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,j,i");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,b,a,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += -1.00 P(i,j) P(a,b) <k,a||c,j>_abab r2_abab(c,b,k,i)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = eri["baab_vovo"]("a,k,c,j") * r2["abab"]("R,c,b,k,i");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,j,i");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,b,a,i,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += +1.00 P(i,j) P(a,b) <l,k||c,d>_aaaa r2_abab(d,b,l,i) t2_abab(c,a,k,j)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = r2["abab"]("R,d,b,l,i") * reused_["6_bbaa_voov"]("a,j,l,d");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,j,i");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,b,a,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += +1.00 P(i,j) P(a,b) <l,k||d,c>_abab r2_abab(d,b,l,i) t2_bbbb(c,a,j,k)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = r2["abab"]("R,d,b,l,i") * reused_["15_bbaa_voov"]("a,j,l,d");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,a,b,j,i");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,b,a,i,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_bbbb += +0.50 <a,b||c,d>_bbbb r2_bbbb(c,d,i,j)
    // flops: o2v2L1 += o2v4L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_bbbb("R,a,b,i,j") += 0.50 * eri["bbbb_vvvv"]("a,b,c,d") * r2["bbbb"]("R,c,d,i,j");

    // sigmar2_bbbb += +0.25 <l,k||c,d>_bbbb r2_bbbb(c,d,i,j) t2_bbbb(a,b,l,k)
    // flops: o2v2L1 += o4v2L1 o4v2L1
    //  mems: o2v2L1 += o4v0L1 o2v2L1
    sigmar2_bbbb("R,a,b,i,j") -= 0.25 * eri["bbbb_oovv"]("k,l,c,d") * r2["bbbb"]("R,c,d,i,j") * t2["bbbb"]("a,b,l,k");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["10_bb_Lvo"]("L,c,k")  = l1["aa"]("L,j,b") * t2["abab"]("b,c,j,k");

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["9_bb_Loo"]("L,k,l")  = l2["bbbb"]("L,j,k,c,b") * t2["bbbb"]("c,b,j,l");

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["8_aa_Loo"]("L,i,k")  = l2["abab"]("L,i,j,c,b") * t2["abab"]("c,b,k,j");

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["7_bb_Loo"]("L,k,l")  = l2["abab"]("L,j,k,c,b") * t2["abab"]("c,b,j,l");

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["6_aa_Loo"]("L,k,l")  = l2["aaaa"]("L,j,k,c,b") * t2["aaaa"]("c,b,j,l");

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["5_aa_Loo"]("L,i,k")  = l2["aaaa"]("L,i,j,c,b") * t2["aaaa"]("c,b,j,k");

    // flops: o0v2L1  = o2v3L1
    //  mems: o0v2L1  = o0v2L1
    tmps_["4_bb_Lvv"]("L,c,d")  = l2["bbbb"]("L,j,k,c,b") * t2["bbbb"]("d,b,j,k");

    // flops: o0v2L1  = o2v3L1
    //  mems: o0v2L1  = o0v2L1
    tmps_["3_bb_Lvv"]("L,c,d")  = l2["abab"]("L,j,k,b,c") * t2["abab"]("b,d,j,k");

    // flops: o0v2L1  = o2v3L1
    //  mems: o0v2L1  = o0v2L1
    tmps_["2_aa_Lvv"]("L,c,d")  = l2["abab"]("L,j,k,c,b") * t2["abab"]("d,b,j,k");

    // flops: o0v2L1  = o2v3L1
    //  mems: o0v2L1  = o0v2L1
    tmps_["1_aa_Lvv"]("L,c,d")  = l2["aaaa"]("L,j,k,c,b") * t2["aaaa"]("d,b,j,k");

    // flops: o1v1L1  = o3v1L1 o1v3L1 o1v1L1 o2v2L1 o1v1L1 o1v3L1 o1v1L1 o3v1L1 o1v1L1 o2v1L1 o1v1L1 o3v2L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o3v2L1 o1v1L1 o2v3L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o2v3L1 o1v1L1 o2v3L1 o1v1L1 o2v3L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o1v2L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o2v1L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o1v1L1 o1v1L1 o1v2L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o3v1L1 o1v1L1 o3v1L1 o1v1L1 o1v3L1 o1v1L1 o1v3L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    tmps_["11_aa_Lov"]("L,i,a")  = eri["aaaa_oovo"]("i,l,a,k") * tmps_["8_aa_Loo"]("L,k,l");
    tmps_["11_aa_Lov"]("L,i,a") += 0.50 * eri["aaaa_vovv"]("c,i,a,d") * tmps_["1_aa_Lvv"]("L,c,d");
    tmps_["11_aa_Lov"]("L,i,a") -= eri["abab_oovv"]("i,k,a,c") * tmps_["10_bb_Lvo"]("L,c,k");
    tmps_["11_aa_Lov"]("L,i,a") += eri["aaaa_vovv"]("c,i,a,d") * tmps_["2_aa_Lvv"]("L,c,d");
    tmps_["11_aa_Lov"]("L,i,a") += eri["abab_oovo"]("i,l,a,k") * tmps_["7_bb_Loo"]("L,k,l");
    tmps_["11_aa_Lov"]("L,i,a") -= 0.50 * f["aa_ov"]("k,a") * tmps_["5_aa_Loo"]("L,i,k");
    tmps_["11_aa_Lov"]("L,i,a") -= 0.50 * l2["aaaa"]("L,j,k,b,a") * reused_["89_aaaa_vooo"]("b,j,k,i");
    tmps_["11_aa_Lov"]("L,i,a") += l1["aa"]("L,j,a") * reused_["92_aa_oo"]("j,i");
    tmps_["11_aa_Lov"]("L,i,a") += 0.50 * l1["aa"]("L,j,a") * reused_["52_aa_oo"]("i,j");
    tmps_["11_aa_Lov"]("L,i,a") += l1["bb"]("L,j,b") * reused_["15_bbaa_voov"]("b,j,i,a");
    tmps_["11_aa_Lov"]("L,i,a") += l2["abab"]("L,k,j,a,b") * reused_["65_aabb_oovo"]("i,k,b,j");
    tmps_["11_aa_Lov"]("L,i,a") += l2["aaaa"]("L,j,k,b,a") * reused_["36_aaaa_vooo"]("b,j,i,k");
    tmps_["11_aa_Lov"]("L,i,a") += f["aa_vo"]("b,j") * l2["aaaa"]("L,i,j,b,a");
    tmps_["11_aa_Lov"]("L,i,a") -= f["aa_vv"]("b,a") * l1["aa"]("L,i,b");
    tmps_["11_aa_Lov"]("L,i,a") -= l2["aaaa"]("L,j,k,b,a") * reused_["29_aaaa_oovo"]("i,k,b,j");
    tmps_["11_aa_Lov"]("L,i,a") += reused_["22_aabb_vvvo"]("c,a,b,j") * l2["abab"]("L,i,j,c,b");
    tmps_["11_aa_Lov"]("L,i,a") -= reused_["8_abab_vovv"]("a,j,c,b") * l2["abab"]("L,i,j,c,b");
    tmps_["11_aa_Lov"]("L,i,a") += l2["abab"]("L,i,j,a,b") * reused_["78_bb_vo"]("b,j");
    tmps_["11_aa_Lov"]("L,i,a") += 0.50 * eri["aaaa_vvvo"]("b,c,a,j") * l2["aaaa"]("L,i,j,c,b");
    tmps_["11_aa_Lov"]("L,i,a") -= reused_["69_aaaa_vvvo"]("c,a,b,j") * l2["aaaa"]("L,i,j,c,b");
    tmps_["11_aa_Lov"]("L,i,a") += 0.25 * reused_["10_aaaa_vovv"]("a,j,c,b") * l2["aaaa"]("L,i,j,c,b");
    tmps_["11_aa_Lov"]("L,i,a") -= eri["abab_vvvo"]("c,b,a,j") * l2["abab"]("L,i,j,c,b");
    tmps_["11_aa_Lov"]("L,i,a") += reused_["70_aaaa_vovv"]("b,j,c,a") * l2["aaaa"]("L,i,j,c,b");
    tmps_["11_aa_Lov"]("L,i,a") += eri["baab_vovo"]("b,i,a,j") * l1["bb"]("L,j,b");
    tmps_["11_aa_Lov"]("L,i,a") += l1["aa"]("L,j,b") * reused_["13_aaaa_voov"]("b,j,i,a");
    tmps_["11_aa_Lov"]("L,i,a") -= l2["abab"]("L,k,j,a,b") * reused_["27_bbaa_vooo"]("b,j,i,k");
    tmps_["11_aa_Lov"]("L,i,a") += eri["aaaa_vovo"]("b,i,a,j") * l1["aa"]("L,j,b");
    tmps_["11_aa_Lov"]("L,i,a") -= l2["abab"]("L,j,k,a,b") * reused_["91_baab_vooo"]("b,i,j,k");
    tmps_["11_aa_Lov"]("L,i,a") += reused_["61_aa_vv"]("b,a") * l1["aa"]("L,i,b");
    tmps_["11_aa_Lov"]("L,i,a") -= l2["abab"]("L,j,k,a,b") * reused_["32_abba_oovo"]("i,k,b,j");
    tmps_["11_aa_Lov"]("L,i,a") -= eri["baab_vooo"]("b,i,j,k") * l2["abab"]("L,j,k,a,b");
    tmps_["11_aa_Lov"]("L,i,a") += l2["aaaa"]("L,i,j,b,a") * reused_["82_aa_vo"]("b,j");
    tmps_["11_aa_Lov"]("L,i,a") -= l1["bb"]("L,j,b") * reused_["5_bbaa_voov"]("b,j,i,a");
    tmps_["11_aa_Lov"]("L,i,a") -= reused_["23_baab_vvvo"]("c,a,b,j") * l2["abab"]("L,i,j,b,c");
    tmps_["11_aa_Lov"]("L,i,a") += f["aa_oo"]("i,j") * l1["aa"]("L,j,a");
    tmps_["11_aa_Lov"]("L,i,a") -= reused_["25_aabb_vvvo"]("c,a,b,j") * l2["abab"]("L,i,j,c,b");
    tmps_["11_aa_Lov"]("L,i,a") -= f["bb_vo"]("b,j") * l2["abab"]("L,i,j,a,b");
    tmps_["11_aa_Lov"]("L,i,a") += 0.50 * eri["aaaa_vooo"]("b,i,j,k") * l2["aaaa"]("L,j,k,b,a");
    tmps_["11_aa_Lov"]("L,i,a") += 0.25 * scalars_["9"] * l1["aa"]("L,i,a");
    tmps_["11_aa_Lov"]("L,i,a") += 0.50 * reused_["62_aa_vv"]("a,b") * l1["aa"]("L,i,b");
    tmps_["11_aa_Lov"]("L,i,a") += f["aa_ov"]("k,a") * tmps_["8_aa_Loo"]("L,i,k");
    tmps_["11_aa_Lov"]("L,i,a") += 0.50 * eri["abab_oovo"]("i,l,a,k") * tmps_["9_bb_Loo"]("L,k,l");
    tmps_["11_aa_Lov"]("L,i,a") += 0.50 * eri["aaaa_oovo"]("i,l,a,k") * tmps_["6_aa_Loo"]("L,k,l");
    tmps_["11_aa_Lov"]("L,i,a") += 0.50 * eri["baab_vovv"]("c,i,a,d") * tmps_["4_bb_Lvv"]("L,c,d");
    tmps_["11_aa_Lov"]("L,i,a") += eri["baab_vovv"]("c,i,a,d") * tmps_["3_bb_Lvv"]("L,c,d");
    tmps_["10_bb_Lvo"].~TArrayD();
    tmps_["5_aa_Loo"].~TArrayD();

    // sigmal1_aa  = -0.50 <i,l||a,k>_aaaa t2_abab(c,b,l,j) l2_abab(k,j,c,b)
    //            += -0.50 <i,l||a,k>_aaaa t2_abab(b,c,l,j) l2_abab(k,j,b,c)
    //            += -0.50 <i,c||d,a>_aaaa t2_aaaa(d,b,j,k) l2_aaaa(j,k,c,b)
    //            += +1.00 <i,k||a,c>_abab t2_abab(b,c,j,k) l1_aa(j,b)
    //            += -0.50 <i,c||d,a>_aaaa t2_abab(d,b,j,k) l2_abab(j,k,c,b)
    //            += -0.50 <i,c||d,a>_aaaa t2_abab(d,b,k,j) l2_abab(k,j,c,b)
    //            += -0.50 <i,l||a,k>_abab t2_abab(c,b,j,l) l2_abab(j,k,c,b)
    //            += -0.50 <i,l||a,k>_abab t2_abab(b,c,j,l) l2_abab(j,k,b,c)
    //            += +0.50 f_aa(k,a) t2_aaaa(c,b,j,k) l2_aaaa(i,j,c,b)
    //            += +0.50 f_aa(i,c) t2_aaaa(c,b,j,k) l2_aaaa(j,k,b,a)
    //            += +0.25 <i,b||c,d>_aaaa t2_aaaa(c,d,j,k) l2_aaaa(j,k,b,a)
    //            += -0.50 <i,k||b,c>_abab t2_abab(b,c,j,k) l1_aa(j,a)
    //            += -0.50 <i,k||c,b>_abab t2_abab(c,b,j,k) l1_aa(j,a)
    //            += -0.50 <i,k||b,c>_aaaa t2_aaaa(b,c,j,k) l1_aa(j,a)
    //            += -1.00 <i,k||a,c>_abab t2_bbbb(c,b,j,k) l1_bb(j,b)
    //            += +1.00 <i,l||k,c>_abab t2_bbbb(c,b,j,l) l2_abab(k,j,a,b)
    //            += -1.00 <i,l||c,k>_aaaa t2_aaaa(c,b,j,l) l2_aaaa(j,k,b,a)
    //            += -1.00 f_aa(b,j) l2_aaaa(i,j,b,a)
    //            += +1.00 f_aa(b,a) l1_aa(i,b)
    //            += -1.00 <i,l||k,c>_abab t2_abab(b,c,j,l) l2_aaaa(j,k,b,a)
    //            += -1.00 <c,k||a,d>_abab t2_bbbb(d,b,j,k) l2_abab(i,j,c,b)
    //            += +0.25 <l,k||a,j>_abab t2_abab(c,b,l,k) l2_abab(i,j,c,b)
    //            += +0.25 <k,l||a,j>_abab t2_abab(c,b,k,l) l2_abab(i,j,c,b)
    //            += +0.25 <l,k||a,j>_abab t2_abab(b,c,l,k) l2_abab(i,j,b,c)
    //            += +0.25 <k,l||a,j>_abab t2_abab(b,c,k,l) l2_abab(i,j,b,c)
    //            += -1.00 f_bb(k,c) t2_bbbb(c,b,j,k) l2_abab(i,j,a,b)
    //            += -0.50 <l,k||c,j>_bbbb t2_bbbb(c,b,l,k) l2_abab(i,j,a,b)
    //            += -0.50 <l,k||c,j>_abab t2_abab(c,b,l,k) l2_abab(i,j,a,b)
    //            += -0.50 <k,l||c,j>_abab t2_abab(c,b,k,l) l2_abab(i,j,a,b)
    //            += +1.00 f_aa(k,c) t2_abab(c,b,k,j) l2_abab(i,j,a,b)
    //            += +0.50 <k,b||c,d>_abab t2_abab(c,d,k,j) l2_abab(i,j,a,b)
    //            += +0.50 <k,b||d,c>_abab t2_abab(d,c,k,j) l2_abab(i,j,a,b)
    //            += -0.50 <k,b||c,d>_bbbb t2_bbbb(c,d,j,k) l2_abab(i,j,a,b)
    //            += +0.50 <c,b||a,j>_aaaa l2_aaaa(i,j,c,b)
    //            += +1.00 <c,k||a,d>_abab t2_abab(b,d,j,k) l2_aaaa(i,j,c,b)
    //            += +0.25 <l,k||a,j>_aaaa t2_aaaa(c,b,l,k) l2_aaaa(i,j,c,b)
    //            += +0.50 <c,b||a,j>_abab l2_abab(i,j,c,b)
    //            += +0.50 <b,c||a,j>_abab l2_abab(i,j,b,c)
    //            += -1.00 <k,c||d,a>_aaaa t2_aaaa(d,b,j,k) l2_aaaa(i,j,c,b)
    //            += +1.00 <i,b||a,j>_abab l1_bb(j,b)
    //            += +1.00 <i,k||c,a>_aaaa t2_aaaa(c,b,j,k) l1_aa(j,b)
    //            += +1.00 <i,l||c,k>_aaaa t2_abab(c,b,l,j) l2_abab(k,j,a,b)
    //            += +1.00 <i,b||a,j>_aaaa l1_aa(j,b)
    //            += -0.25 <i,b||c,d>_abab t2_abab(c,d,j,k) l2_abab(j,k,a,b)
    //            += -0.25 <i,b||c,d>_abab t2_abab(c,d,k,j) l2_abab(k,j,a,b)
    //            += -0.25 <i,b||d,c>_abab t2_abab(d,c,j,k) l2_abab(j,k,a,b)
    //            += -0.25 <i,b||d,c>_abab t2_abab(d,c,k,j) l2_abab(k,j,a,b)
    //            += -0.50 f_aa(i,c) t2_abab(c,b,j,k) l2_abab(j,k,a,b)
    //            += -0.50 f_aa(i,c) t2_abab(c,b,k,j) l2_abab(k,j,a,b)
    //            += -0.50 <k,j||a,c>_abab t2_abab(b,c,k,j) l1_aa(i,b)
    //            += -0.50 <j,k||a,c>_abab t2_abab(b,c,j,k) l1_aa(i,b)
    //            += +1.00 <i,l||c,k>_abab t2_abab(c,b,j,l) l2_abab(j,k,a,b)
    //            += -0.50 <i,b||j,k>_abab l2_abab(j,k,a,b)
    //            += -0.50 <i,b||k,j>_abab l2_abab(k,j,a,b)
    //            += -1.00 f_bb(k,c) t2_abab(b,c,j,k) l2_aaaa(i,j,b,a)
    //            += +0.50 <l,k||j,c>_abab t2_abab(b,c,l,k) l2_aaaa(i,j,b,a)
    //            += +0.50 <k,l||j,c>_abab t2_abab(b,c,k,l) l2_aaaa(i,j,b,a)
    //            += -0.50 <b,k||c,d>_abab t2_abab(c,d,j,k) l2_aaaa(i,j,b,a)
    //            += -0.50 <b,k||d,c>_abab t2_abab(d,c,j,k) l2_aaaa(i,j,b,a)
    //            += +1.00 f_aa(k,c) t2_aaaa(c,b,j,k) l2_aaaa(i,j,b,a)
    //            += +0.50 <l,k||c,j>_aaaa t2_aaaa(c,b,l,k) l2_aaaa(i,j,b,a)
    //            += +0.50 <k,b||c,d>_aaaa t2_aaaa(c,d,j,k) l2_aaaa(i,j,b,a)
    //            += -1.00 <i,k||c,a>_aaaa t2_abab(c,b,k,j) l1_bb(j,b)
    //            += -1.00 <k,c||a,d>_abab t2_abab(b,d,k,j) l2_abab(i,j,b,c)
    //            += -1.00 f_aa(i,j) l1_aa(j,a)
    //            += +1.00 <k,c||d,a>_aaaa t2_abab(d,b,k,j) l2_abab(i,j,c,b)
    //            += +1.00 f_bb(b,j) l2_abab(i,j,a,b)
    //            += +0.50 <i,b||j,k>_aaaa l2_aaaa(j,k,b,a)
    //            += +0.25 <k,j||b,c>_bbbb t2_bbbb(b,c,k,j) l1_aa(i,a)
    //            += +1.00 f_aa(j,j) l1_aa(i,a)
    //            += +1.00 f_bb(j,j) l1_aa(i,a)
    //            += -0.50 <k,j||k,j>_aaaa l1_aa(i,a)
    //            += -0.50 <k,j||k,j>_bbbb l1_aa(i,a)
    //            += +0.25 <k,j||b,c>_aaaa t2_aaaa(b,c,k,j) l1_aa(i,a)
    //            += -0.50 <k,j||k,j>_abab l1_aa(i,a)
    //            += -0.50 <j,k||j,k>_abab l1_aa(i,a)
    //            += +0.25 <k,j||b,c>_abab t2_abab(b,c,k,j) l1_aa(i,a)
    //            += +0.25 <j,k||b,c>_abab t2_abab(b,c,j,k) l1_aa(i,a)
    //            += +0.25 <k,j||c,b>_abab t2_abab(c,b,k,j) l1_aa(i,a)
    //            += +0.25 <j,k||c,b>_abab t2_abab(c,b,j,k) l1_aa(i,a)
    //            += -0.50 <k,j||c,a>_aaaa t2_aaaa(c,b,k,j) l1_aa(i,b)
    //            += -0.50 f_aa(k,a) t2_abab(c,b,k,j) l2_abab(i,j,c,b)
    //            += -0.50 f_aa(k,a) t2_abab(b,c,k,j) l2_abab(i,j,b,c)
    //            += -0.50 <i,l||a,k>_abab t2_bbbb(c,b,j,l) l2_bbbb(j,k,c,b)
    //            += -0.50 <i,l||a,k>_aaaa t2_aaaa(c,b,j,l) l2_aaaa(j,k,c,b)
    //            += +0.50 <i,c||a,d>_abab t2_bbbb(d,b,j,k) l2_bbbb(j,k,c,b)
    //            += +0.50 <i,c||a,d>_abab t2_abab(b,d,j,k) l2_abab(j,k,b,c)
    //            += +0.50 <i,c||a,d>_abab t2_abab(b,d,k,j) l2_abab(k,j,b,c)
    sigmal1_aa("L,a,i")  = -1.00 * tmps_["11_aa_Lov"]("L,i,a");
    tmps_["11_aa_Lov"].~TArrayD();

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["12_bb_Loo"]("L,i,k")  = l2["bbbb"]("L,i,j,c,b") * t2["bbbb"]("c,b,j,k");
}