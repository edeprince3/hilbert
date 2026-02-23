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

void hilbert::EOM_EE_CCSD::sigma_ee_00_2() {

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


    // sigmal2_bbbb += +0.25 <i,j||c,d>_bbbb t2_bbbb(c,d,k,l) l2_bbbb(k,l,a,b)
    // flops: o2v2L1 += o4v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_bbbb("L,a,b,i,j") += 0.25 * l2["bbbb"]("L,k,l,a,b") * reused_["68_bbbb_oooo"]("k,l,i,j");

    // sigmal2_bbbb += +1.00 P(i,j) P(a,b) <j,c||a,k>_bbbb l2_bbbb(i,k,c,b)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j")  = l2["bbbb"]("L,i,k,c,b") * eri["bbbb_vovo"]("c,j,a,k");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,a,b,j,i");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,b,a,i,j");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,b,a,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmal2_bbbb += +1.00 P(i,j) P(a,b) <l,j||d,a>_abab t2_abab(d,c,l,k) l2_bbbb(i,k,c,b)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j")  = l2["bbbb"]("L,i,k,c,b") * reused_["4_bbbb_voov"]("c,k,j,a");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,a,b,j,i");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,b,a,i,j");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,b,a,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmal2_bbbb += +1.00 P(i,j) P(a,b) <j,l||d,a>_bbbb t2_bbbb(d,c,k,l) l2_bbbb(i,k,c,b)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j")  = l2["bbbb"]("L,i,k,c,b") * reused_["41_bbbb_voov"]("c,k,j,a");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,a,b,j,i");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,b,a,i,j");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,b,a,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmal2_bbbb += -1.00 P(i,j) P(a,b) <c,j||k,a>_abab l2_abab(k,i,c,b)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j")  = l2["abab"]("L,k,i,c,b") * eri["abba_vovo"]("c,j,a,k");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,a,b,j,i");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,b,a,i,j");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,b,a,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmal2_bbbb += +1.00 P(i,j) P(a,b) <l,j||d,a>_abab t2_aaaa(d,c,k,l) l2_abab(k,i,c,b)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j")  = l2["abab"]("L,k,i,c,b") * reused_["7_aabb_voov"]("c,k,j,a");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,a,b,j,i");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,b,a,i,j");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,b,a,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmal2_bbbb += +1.00 P(i,j) P(a,b) <j,l||d,a>_bbbb t2_abab(c,d,k,l) l2_abab(k,i,c,b)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j")  = l2["abab"]("L,k,i,c,b") * reused_["46_aabb_voov"]("c,k,j,a");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,a,b,j,i");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,b,a,i,j");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,b,a,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmal2_bbbb += +0.50 <d,c||a,b>_bbbb l2_bbbb(i,j,d,c)
    // flops: o2v2L1 += o2v4L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_bbbb("L,a,b,i,j") -= 0.50 * eri["bbbb_vvvv"]("c,d,a,b") * l2["bbbb"]("L,i,j,d,c");

    // sigmal2_bbbb += +0.25 <l,k||a,b>_bbbb t2_bbbb(d,c,l,k) l2_bbbb(i,j,d,c)
    // flops: o2v2L1 += o4v2L1 o4v2L1
    //  mems: o2v2L1 += o4v0L1 o2v2L1
    sigmal2_bbbb("L,a,b,i,j") -= 0.25 * l2["bbbb"]("L,i,j,d,c") * t2["bbbb"]("d,c,l,k") * eri["bbbb_oovv"]("k,l,a,b");

    // sigmar2_aaaa += +0.25 <l,k||c,d>_bbbb r2_aaaa(a,b,i,j) t2_bbbb(c,d,l,k)
    //              += +1.00 f_aa(k,k) r2_aaaa(a,b,i,j)
    //              += +1.00 f_bb(k,k) r2_aaaa(a,b,i,j)
    //              += -0.50 <l,k||l,k>_aaaa r2_aaaa(a,b,i,j)
    //              += -0.50 <l,k||l,k>_bbbb r2_aaaa(a,b,i,j)
    //              += +0.25 <l,k||c,d>_aaaa r2_aaaa(a,b,i,j) t2_aaaa(c,d,l,k)
    //              += -0.50 <l,k||l,k>_abab r2_aaaa(a,b,i,j)
    //              += -0.50 <k,l||k,l>_abab r2_aaaa(a,b,i,j)
    //              += +0.25 <l,k||c,d>_abab r2_aaaa(a,b,i,j) t2_abab(c,d,l,k)
    //              += +0.25 <k,l||c,d>_abab r2_aaaa(a,b,i,j) t2_abab(c,d,k,l)
    //              += +0.25 <l,k||d,c>_abab r2_aaaa(a,b,i,j) t2_abab(d,c,l,k)
    //              += +0.25 <k,l||d,c>_abab r2_aaaa(a,b,i,j) t2_abab(d,c,k,l)
    // flops: o2v2L1 += o2v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_aaaa("R,a,b,i,j") -= 0.25 * scalars_["9"] * r2["aaaa"]("R,a,b,i,j");

    // sigmar2_aaaa += -1.00 P(i,j) P(a,b) f_bb(k,c) r1_aa(b,i) t2_abab(a,c,j,k)
    //              += +0.50 P(i,j) P(a,b) <l,k||j,c>_abab r1_aa(b,i) t2_abab(a,c,l,k)
    //              += +0.50 P(i,j) P(a,b) <k,l||j,c>_abab r1_aa(b,i) t2_abab(a,c,k,l)
    //              += -0.50 P(i,j) P(a,b) <a,k||c,d>_abab r1_aa(b,i) t2_abab(c,d,j,k)
    //              += -0.50 P(i,j) P(a,b) <a,k||d,c>_abab r1_aa(b,i) t2_abab(d,c,j,k)
    //              += +1.00 P(i,j) P(a,b) f_aa(k,c) r1_aa(b,i) t2_aaaa(c,a,j,k)
    //              += +0.50 P(i,j) P(a,b) <l,k||c,j>_aaaa r1_aa(b,i) t2_aaaa(c,a,l,k)
    //              += +0.50 P(i,j) P(a,b) <k,a||c,d>_aaaa r1_aa(b,i) t2_aaaa(c,d,j,k)
    // flops: o2v2L1 += o2v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = r1["aa"]("R,b,i") * reused_["82_aa_vo"]("a,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,j,i");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,b,a,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_aaaa += +1.00 P(a,b) <k,a||i,j>_aaaa r1_aa(b,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = eri["aaaa_vooo"]("a,k,i,j") * r1["aa"]("R,b,k");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,b,a,i,j");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_aaaa += -1.00 P(i,j) f_aa(k,j) r2_aaaa(a,b,i,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = r2["aaaa"]("R,a,b,i,k") * f["aa_oo"]("k,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_aaaa += -0.50 P(i,j) <l,k||c,d>_aaaa r2_aaaa(a,b,i,l) t2_aaaa(c,d,j,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = 0.50 * r2["aaaa"]("R,a,b,i,l") * reused_["53_aa_oo"]("j,l");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,a,b,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_aaaa += -1.00 P(i,j) P(a,b) <l,k||c,j>_aaaa r1_aa(b,l) t2_aaaa(c,a,i,k)
    //              += -1.00 P(i,j) P(a,b) <l,k||j,c>_abab r1_aa(b,l) t2_abab(a,c,i,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = r1["aa"]("R,b,l") * reused_["86_aaaa_oovo"]("l,j,a,i");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,a,b,j,i");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,b,a,i,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_aaaa += +1.00 P(a,b) f_aa(k,c) r1_aa(b,k) t2_aaaa(c,a,i,j)
    //              += +0.50 P(a,b) <k,a||c,d>_aaaa r1_aa(b,k) t2_aaaa(c,d,i,j)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = r1["aa"]("R,b,k") * reused_["89_aaaa_vooo"]("a,i,j,k");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,b,a,i,j");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_aaaa += -0.50 P(i,j) <l,k||c,d>_abab r2_aaaa(a,b,i,l) t2_abab(c,d,j,k)
    //              += -0.50 P(i,j) <l,k||d,c>_abab r2_aaaa(a,b,i,l) t2_abab(d,c,j,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = r2["aaaa"]("R,a,b,i,l") * reused_["92_aa_oo"]("j,l");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_aaaa += +1.00 P(i,j) <a,b||c,j>_aaaa r1_aa(c,i)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = eri["aaaa_vvvo"]("a,b,c,j") * r1["aa"]("R,c,i");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,a,b,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_aaaa += +1.00 P(a,b) f_aa(a,c) r2_aaaa(c,b,i,j)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = f["aa_vv"]("a,c") * r2["aaaa"]("R,c,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,b,a,i,j");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_aaaa += +0.50 P(i,j) <l,k||c,j>_aaaa r1_aa(c,i) t2_aaaa(a,b,l,k)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = 0.50 * r1["aa"]("R,c,i") * reused_["10_aaaa_vovv"]("c,j,a,b");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_aaaa += -0.50 P(a,b) <l,k||d,c>_abab r2_aaaa(d,b,i,j) t2_abab(a,c,l,k)
    //              += -0.50 P(a,b) <k,l||d,c>_abab r2_aaaa(d,b,i,j) t2_abab(a,c,k,l)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = r2["aaaa"]("R,d,b,i,j") * reused_["61_aa_vv"]("a,d");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,b,a,i,j");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_aaaa += -0.50 P(a,b) <l,k||c,d>_aaaa r2_aaaa(d,b,i,j) t2_aaaa(c,a,l,k)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = 0.50 * r2["aaaa"]("R,d,b,i,j") * reused_["63_aa_vv"]("a,d");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,b,a,i,j");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_aaaa += +1.00 P(i,j) P(a,b) <a,k||d,c>_abab r1_aa(d,i) t2_abab(b,c,j,k)
    //              += -1.00 P(i,j) P(a,b) <k,a||c,d>_aaaa r1_aa(d,i) t2_aaaa(c,b,j,k)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = r1["aa"]("R,d,i") * reused_["74_aaaa_vvvo"]("a,d,b,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,a,b,j,i");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,b,a,i,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_aaaa += +0.50 <l,k||i,j>_aaaa r2_aaaa(a,b,l,k)
    // flops: o2v2L1 += o4v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_aaaa("R,a,b,i,j") -= 0.50 * r2["aaaa"]("R,a,b,l,k") * eri["aaaa_oooo"]("k,l,i,j");

    // sigmar2_aaaa += +0.25 <l,k||c,d>_aaaa r2_aaaa(a,b,l,k) t2_aaaa(c,d,i,j)
    // flops: o2v2L1 += o4v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_aaaa("R,a,b,i,j") -= 0.25 * r2["aaaa"]("R,a,b,l,k") * reused_["26_aaaa_oooo"]("i,j,k,l");

    // sigmar2_aaaa += -1.00 P(i,j) P(a,b) <a,k||j,c>_abab r2_abab(b,c,i,k)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = eri["abba_vovo"]("a,k,c,j") * r2["abab"]("R,b,c,i,k");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,a,b,j,i");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,b,a,i,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_aaaa += +1.00 P(i,j) P(a,b) <k,l||c,d>_abab r2_abab(b,d,i,l) t2_aaaa(c,a,j,k)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = r2["abab"]("R,b,d,i,l") * reused_["7_aabb_voov"]("a,j,l,d");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,a,b,j,i");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,b,a,i,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_aaaa += +1.00 P(i,j) P(a,b) <l,k||c,d>_bbbb r2_abab(b,d,i,l) t2_abab(a,c,j,k)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = r2["abab"]("R,b,d,i,l") * reused_["47_aabb_voov"]("a,j,l,d");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,j,i");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,b,a,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_aaaa += +1.00 P(i,j) P(a,b) <k,a||c,j>_aaaa r2_aaaa(c,b,i,k)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = eri["aaaa_vovo"]("a,k,c,j") * r2["aaaa"]("R,c,b,i,k");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,j,i");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,b,a,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_aaaa += +1.00 P(i,j) P(a,b) <l,k||c,d>_aaaa r2_aaaa(d,b,i,l) t2_aaaa(c,a,j,k)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = r2["aaaa"]("R,d,b,i,l") * reused_["12_aaaa_voov"]("a,j,l,d");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,j,i");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,b,a,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_aaaa += +0.50 <a,b||c,d>_aaaa r2_aaaa(c,d,i,j)
    // flops: o2v2L1 += o2v4L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_aaaa("R,a,b,i,j") += 0.50 * eri["aaaa_vvvv"]("a,b,c,d") * r2["aaaa"]("R,c,d,i,j");

    // sigmar2_aaaa += +0.25 <l,k||c,d>_aaaa r2_aaaa(c,d,i,j) t2_aaaa(a,b,l,k)
    // flops: o2v2L1 += o4v2L1 o4v2L1
    //  mems: o2v2L1 += o4v0L1 o2v2L1
    sigmar2_aaaa("R,a,b,i,j") -= 0.25 * eri["aaaa_oovv"]("k,l,c,d") * r2["aaaa"]("R,c,d,i,j") * t2["aaaa"]("a,b,l,k");

    // sigmar2_aaaa += +1.00 P(i,j) P(a,b) <l,k||d,c>_abab r2_aaaa(d,b,i,l) t2_abab(a,c,j,k)
    // flops: o2v2L1 += o3v3L1 o3v3L1
    //  mems: o2v2L1 += o2v2L1 o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = eri["abab_oovv"]("l,k,d,c") * r2["aaaa"]("R,d,b,i,l") * t2["abab"]("a,c,j,k");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,a,b,j,i");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,b,a,i,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_abab += +1.00 f_aa(a,i) r1_bb(b,j)
    // flops: o2v2L1 += o2v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += f["aa_vo"]("a,i") * r1["bb"]("R,b,j");

    // sigmar2_abab += +0.25 <l,k||c,d>_bbbb r2_abab(a,b,i,j) t2_bbbb(c,d,l,k)
    //              += +1.00 f_aa(k,k) r2_abab(a,b,i,j)
    //              += +1.00 f_bb(k,k) r2_abab(a,b,i,j)
    //              += -0.50 <l,k||l,k>_aaaa r2_abab(a,b,i,j)
    //              += -0.50 <l,k||l,k>_bbbb r2_abab(a,b,i,j)
    //              += +0.25 <l,k||c,d>_aaaa r2_abab(a,b,i,j) t2_aaaa(c,d,l,k)
    //              += -0.50 <l,k||l,k>_abab r2_abab(a,b,i,j)
    //              += -0.50 <k,l||k,l>_abab r2_abab(a,b,i,j)
    //              += +0.25 <l,k||c,d>_abab r2_abab(a,b,i,j) t2_abab(c,d,l,k)
    //              += +0.25 <k,l||c,d>_abab r2_abab(a,b,i,j) t2_abab(c,d,k,l)
    //              += +0.25 <l,k||d,c>_abab r2_abab(a,b,i,j) t2_abab(d,c,l,k)
    //              += +0.25 <k,l||d,c>_abab r2_abab(a,b,i,j) t2_abab(d,c,k,l)
    // flops: o2v2L1 += o2v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= 0.25 * scalars_["9"] * r2["abab"]("R,a,b,i,j");

    // sigmar2_abab += -1.00 f_bb(k,c) r1_aa(a,i) t2_bbbb(c,b,j,k)
    //              += -0.50 <l,k||c,j>_bbbb r1_aa(a,i) t2_bbbb(c,b,l,k)
    //              += -0.50 <l,k||c,j>_abab r1_aa(a,i) t2_abab(c,b,l,k)
    //              += -0.50 <k,l||c,j>_abab r1_aa(a,i) t2_abab(c,b,k,l)
    //              += +1.00 f_aa(k,c) r1_aa(a,i) t2_abab(c,b,k,j)
    //              += +0.50 <k,b||c,d>_abab r1_aa(a,i) t2_abab(c,d,k,j)
    //              += +0.50 <k,b||d,c>_abab r1_aa(a,i) t2_abab(d,c,k,j)
    //              += -0.50 <k,b||c,d>_bbbb r1_aa(a,i) t2_bbbb(c,d,j,k)
    // flops: o2v2L1 += o2v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= r1["aa"]("R,a,i") * reused_["78_bb_vo"]("b,j");

    // sigmar2_abab += +1.00 f_bb(k,c) r1_bb(b,j) t2_abab(a,c,i,k)
    //              += -0.50 <l,k||i,c>_abab r1_bb(b,j) t2_abab(a,c,l,k)
    //              += -0.50 <k,l||i,c>_abab r1_bb(b,j) t2_abab(a,c,k,l)
    //              += +0.50 <a,k||c,d>_abab r1_bb(b,j) t2_abab(c,d,i,k)
    //              += +0.50 <a,k||d,c>_abab r1_bb(b,j) t2_abab(d,c,i,k)
    //              += -1.00 f_aa(k,c) r1_bb(b,j) t2_aaaa(c,a,i,k)
    //              += -0.50 <l,k||c,i>_aaaa r1_bb(b,j) t2_aaaa(c,a,l,k)
    //              += -0.50 <k,a||c,d>_aaaa r1_bb(b,j) t2_aaaa(c,d,i,k)
    // flops: o2v2L1 += o2v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += r1["bb"]("R,b,j") * reused_["82_aa_vo"]("a,i");

    // sigmar2_abab += -1.00 f_aa(k,i) r2_abab(a,b,k,j)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= r2["abab"]("R,a,b,k,j") * f["aa_oo"]("k,i");

    // sigmar2_abab += -1.00 <k,b||i,j>_abab r1_aa(a,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += r1["aa"]("R,a,k") * eri["baab_vooo"]("b,k,i,j");

    // sigmar2_abab += -0.50 <l,k||c,d>_aaaa r2_abab(a,b,l,j) t2_aaaa(c,d,i,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += 0.50 * r2["abab"]("R,a,b,l,j") * reused_["53_aa_oo"]("i,l");

    // sigmar2_abab += +1.00 <l,k||i,c>_abab r1_aa(a,l) t2_bbbb(c,b,j,k)
    //              += +1.00 <l,k||c,i>_aaaa r1_aa(a,l) t2_abab(c,b,k,j)
    //              += +1.00 <l,k||c,j>_abab r1_aa(a,l) t2_abab(c,b,i,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= r1["aa"]("R,a,l") * reused_["83_aabb_oovo"]("l,i,b,j");

    // sigmar2_abab += -0.50 <k,b||c,d>_abab r1_aa(a,k) t2_abab(c,d,i,j)
    //              += -0.50 <k,b||d,c>_abab r1_aa(a,k) t2_abab(d,c,i,j)
    //              += -1.00 f_aa(k,c) r1_aa(a,k) t2_abab(c,b,i,j)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += r1["aa"]("R,a,k") * reused_["91_baab_vooo"]("b,k,i,j");

    // sigmar2_abab += -0.50 <l,k||c,d>_abab r2_abab(a,b,l,j) t2_abab(c,d,i,k)
    //              += -0.50 <l,k||d,c>_abab r2_abab(a,b,l,j) t2_abab(d,c,i,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= r2["abab"]("R,a,b,l,j") * reused_["92_aa_oo"]("i,l");

    // sigmar2_abab += -1.00 f_bb(k,j) r2_abab(a,b,i,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= r2["abab"]("R,a,b,i,k") * f["bb_oo"]("k,j");

    // sigmar2_abab += -1.00 <a,k||i,j>_abab r1_bb(b,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= eri["abab_vooo"]("a,k,i,j") * r1["bb"]("R,b,k");

    // sigmar2_abab += -0.50 <k,l||c,d>_abab r2_abab(a,b,i,l) t2_abab(c,d,k,j)
    //              += -0.50 <k,l||d,c>_abab r2_abab(a,b,i,l) t2_abab(d,c,k,j)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= r2["abab"]("R,a,b,i,l") * reused_["50_bb_oo"]("j,l");

    // sigmar2_abab += +1.00 <l,k||c,j>_bbbb r1_bb(b,l) t2_abab(a,c,i,k)
    //              += +1.00 <k,l||c,j>_abab r1_bb(b,l) t2_aaaa(c,a,i,k)
    //              += +1.00 <k,l||i,c>_abab r1_bb(b,l) t2_abab(a,c,k,j)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= r1["bb"]("R,b,l") * reused_["84_bbaa_oovo"]("l,j,a,i");

    // sigmar2_abab += -0.50 <l,k||c,d>_bbbb r2_abab(a,b,i,l) t2_bbbb(c,d,j,k)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += 0.50 * r2["abab"]("R,a,b,i,l") * reused_["88_bb_oo"]("j,l");

    // sigmar2_abab += -0.50 <a,k||c,d>_abab r1_bb(b,k) t2_abab(c,d,i,j)
    //              += -0.50 <a,k||d,c>_abab r1_bb(b,k) t2_abab(d,c,i,j)
    //              += -1.00 f_bb(k,c) r1_bb(b,k) t2_abab(a,c,i,j)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= r1["bb"]("R,b,k") * reused_["93_abab_vooo"]("a,k,i,j");

    // sigmar2_abab += +1.00 <a,b||i,c>_abab r1_bb(c,j)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= eri["abba_vvvo"]("a,b,c,i") * r1["bb"]("R,c,j");

    // sigmar2_abab += +1.00 f_bb(b,c) r2_abab(a,c,i,j)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += r2["abab"]("R,a,c,i,j") * f["bb_vv"]("b,c");

    // sigmar2_abab += -1.00 <k,b||c,d>_abab r1_bb(d,j) t2_aaaa(c,a,i,k)
    //              += +1.00 <k,b||c,d>_bbbb r1_bb(d,j) t2_abab(a,c,i,k)
    //              += -1.00 <a,k||c,d>_abab r1_bb(d,j) t2_abab(c,b,i,k)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += r1["bb"]("R,d,j") * reused_["20_bbaa_vvvo"]("b,d,a,i");

    // sigmar2_abab += +0.50 <l,k||i,c>_abab r1_bb(c,j) t2_abab(a,b,l,k)
    //              += +0.50 <k,l||i,c>_abab r1_bb(c,j) t2_abab(a,b,k,l)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= r1["bb"]("R,c,j") * reused_["45_baab_vovv"]("c,i,a,b");

    // sigmar2_abab += -0.50 <l,k||c,d>_bbbb r2_abab(a,d,i,j) t2_bbbb(c,b,l,k)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += 0.50 * r2["abab"]("R,a,d,i,j") * reused_["56_bb_vv"]("b,d");

    // sigmar2_abab += -0.50 <l,k||c,d>_abab r2_abab(a,d,i,j) t2_abab(c,b,l,k)
    //              += -0.50 <k,l||c,d>_abab r2_abab(a,d,i,j) t2_abab(c,b,k,l)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= r2["abab"]("R,a,d,i,j") * reused_["59_bb_vv"]("d,b");

    // sigmar2_abab += +1.00 <a,b||c,j>_abab r1_aa(c,i)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += eri["abab_vvvo"]("a,b,c,j") * r1["aa"]("R,c,i");

    // sigmar2_abab += +1.00 f_aa(a,c) r2_abab(c,b,i,j)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += f["aa_vv"]("a,c") * r2["abab"]("R,c,b,i,j");

    // sigmar2_abab += +0.50 <l,k||c,j>_abab r1_aa(c,i) t2_abab(a,b,l,k)
    //              += +0.50 <k,l||c,j>_abab r1_aa(c,i) t2_abab(a,b,k,l)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += r1["aa"]("R,c,i") * reused_["8_abab_vovv"]("c,j,a,b");

    // sigmar2_abab += -1.00 <k,b||d,c>_abab r1_aa(d,i) t2_abab(a,c,k,j)
    //              += +1.00 <k,a||c,d>_aaaa r1_aa(d,i) t2_abab(c,b,k,j)
    //              += -1.00 <a,k||d,c>_abab r1_aa(d,i) t2_bbbb(c,b,j,k)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += r1["aa"]("R,d,i") * reused_["24_baab_vvvo"]("b,d,a,j");

    // sigmar2_abab += -0.50 <l,k||d,c>_abab r2_abab(d,b,i,j) t2_abab(a,c,l,k)
    //              += -0.50 <k,l||d,c>_abab r2_abab(d,b,i,j) t2_abab(a,c,k,l)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= r2["abab"]("R,d,b,i,j") * reused_["61_aa_vv"]("a,d");

    // sigmar2_abab += -0.50 <l,k||c,d>_aaaa r2_abab(d,b,i,j) t2_aaaa(c,a,l,k)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += 0.50 * r2["abab"]("R,d,b,i,j") * reused_["63_aa_vv"]("a,d");

    // sigmar2_abab += +0.50 <l,k||i,j>_abab r2_abab(a,b,l,k)
    //              += +0.50 <k,l||i,j>_abab r2_abab(a,b,k,l)
    // flops: o2v2L1 += o4v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += r2["abab"]("R,a,b,l,k") * eri["abab_oooo"]("l,k,i,j");

    // sigmar2_abab += +0.25 <l,k||c,d>_abab r2_abab(a,b,l,k) t2_abab(c,d,i,j)
    //              += +0.25 <l,k||d,c>_abab r2_abab(a,b,l,k) t2_abab(d,c,i,j)
    //              += +0.25 <k,l||c,d>_abab r2_abab(a,b,k,l) t2_abab(c,d,i,j)
    //              += +0.25 <k,l||d,c>_abab r2_abab(a,b,k,l) t2_abab(d,c,i,j)
    // flops: o2v2L1 += o4v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += r2["abab"]("R,a,b,l,k") * reused_["30_abab_oooo"]("l,k,i,j");

    // sigmar2_abab += -1.00 <k,b||i,c>_abab r2_abab(a,c,k,j)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= r2["abab"]("R,a,c,k,j") * eri["baba_vovo"]("b,k,c,i");

    // sigmar2_abab += +1.00 <l,k||c,d>_abab r2_abab(a,d,l,j) t2_abab(c,b,i,k)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += r2["abab"]("R,a,d,l,j") * reused_["44_baab_voov"]("b,i,l,d");

    // sigmar2_abab += +1.00 <k,b||c,j>_bbbb r2_abab(a,c,i,k)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") -= r2["abab"]("R,a,c,i,k") * eri["bbbb_vovo"]("b,k,c,j");

    // sigmar2_abab += -1.00 <a,k||i,c>_abab r2_bbbb(c,b,j,k)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmar2_abab("R,a,b,i,j") += eri["abba_vovo"]("a,k,c,i") * r2["bbbb"]("R,c,b,j,k");
}