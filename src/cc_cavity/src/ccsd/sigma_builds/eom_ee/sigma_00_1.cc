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

void hilbert::EOM_EE_CCSD::sigma_ee_00_1() {

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

    // sigmar2_bbbb  = -1.00 P(i,j) P(a,b) f_bb(a,j) r1_bb(b,i)
    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j")  = f["bb_vo"]("a,j") * r1["bb"]("R,b,i");
    sigmar2_bbbb("R,a,b,i,j")  = -1.00 * tmps_["perm_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,a,b,j,i");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("R,b,a,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmal2_abab  = +1.00 f_aa(i,a) l1_bb(j,b)
    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    sigmal2_abab("L,a,b,i,j")  = f["aa_ov"]("i,a") * l1["bb"]("L,j,b");

    // sigmal2_bbbb  = -1.00 P(i,j) P(a,b) f_bb(j,a) l1_bb(i,b)
    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j")  = l1["bb"]("L,i,b") * f["bb_ov"]("j,a");
    sigmal2_bbbb("L,a,b,i,j")  = -1.00 * tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,a,b,j,i");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,b,a,i,j");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,b,a,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmar2_aaaa  = -1.00 P(i,j) P(a,b) f_aa(a,j) r1_aa(b,i)
    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j")  = f["aa_vo"]("a,j") * r1["aa"]("R,b,i");
    sigmar2_aaaa("R,a,b,i,j")  = -1.00 * tmps_["perm_aaaa_Lvvoo"]("R,a,b,i,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,a,b,j,i");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("R,b,a,i,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("R,b,a,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmar2_abab  = +1.00 f_bb(b,j) r1_aa(a,i)
    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    sigmar2_abab("R,a,b,i,j")  = r1["aa"]("R,a,i") * f["bb_vo"]("b,j");

    // sigmal2_aaaa  = -1.00 P(i,j) P(a,b) f_aa(j,a) l1_aa(i,b)
    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j")  = l1["aa"]("L,i,b") * f["aa_ov"]("j,a");
    sigmal2_aaaa("L,a,b,i,j")  = -1.00 * tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,a,b,j,i");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,b,a,i,j");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,b,a,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmal2_aaaa += +0.25 <l,k||c,d>_bbbb t2_bbbb(c,d,l,k) l2_aaaa(i,j,a,b)
    //              += +1.00 f_aa(k,k) l2_aaaa(i,j,a,b)
    //              += +1.00 f_bb(k,k) l2_aaaa(i,j,a,b)
    //              += -0.50 <l,k||l,k>_aaaa l2_aaaa(i,j,a,b)
    //              += -0.50 <l,k||l,k>_bbbb l2_aaaa(i,j,a,b)
    //              += +0.25 <l,k||c,d>_aaaa t2_aaaa(c,d,l,k) l2_aaaa(i,j,a,b)
    //              += -0.50 <l,k||l,k>_abab l2_aaaa(i,j,a,b)
    //              += -0.50 <k,l||k,l>_abab l2_aaaa(i,j,a,b)
    //              += +0.25 <l,k||c,d>_abab t2_abab(c,d,l,k) l2_aaaa(i,j,a,b)
    //              += +0.25 <k,l||c,d>_abab t2_abab(c,d,k,l) l2_aaaa(i,j,a,b)
    //              += +0.25 <l,k||d,c>_abab t2_abab(d,c,l,k) l2_aaaa(i,j,a,b)
    //              += +0.25 <k,l||d,c>_abab t2_abab(d,c,k,l) l2_aaaa(i,j,a,b)
    // flops: o2v2L1 += o2v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_aaaa("L,a,b,i,j") -= 0.25 * scalars_["9"] * l2["aaaa"]("L,i,j,a,b");

    // sigmal2_aaaa += -1.00 P(a,b) <i,j||a,k>_aaaa l1_aa(k,b)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j")  = l1["aa"]("L,k,b") * eri["aaaa_oovo"]("i,j,a,k");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,b,a,i,j");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmal2_aaaa += -1.00 P(i,j) f_aa(j,k) l2_aaaa(i,k,a,b)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j")  = l2["aaaa"]("L,i,k,a,b") * f["aa_oo"]("j,k");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,a,b,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmal2_aaaa += -0.50 P(i,j) <j,l||c,d>_aaaa t2_aaaa(c,d,k,l) l2_aaaa(i,k,a,b)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j")  = 0.50 * l2["aaaa"]("L,i,k,a,b") * reused_["52_aa_oo"]("j,k");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,a,b,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmal2_aaaa += -0.50 P(i,j) <j,l||c,d>_abab t2_abab(c,d,k,l) l2_aaaa(i,k,a,b)
    //              += -0.50 P(i,j) <j,l||d,c>_abab t2_abab(d,c,k,l) l2_aaaa(i,k,a,b)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j")  = l2["aaaa"]("L,i,k,a,b") * reused_["92_aa_oo"]("k,j");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,a,b,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmal2_aaaa += -1.00 P(i,j) <j,c||a,b>_aaaa l1_aa(i,c)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j")  = l1["aa"]("L,i,c") * eri["aaaa_vovv"]("c,j,a,b");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,a,b,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmal2_aaaa += +1.00 P(a,b) f_aa(c,a) l2_aaaa(i,j,c,b)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j")  = f["aa_vv"]("c,a") * l2["aaaa"]("L,i,j,c,b");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,b,a,i,j");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmal2_aaaa += -0.50 P(a,b) <l,k||a,d>_abab t2_abab(c,d,l,k) l2_aaaa(i,j,c,b)
    //              += -0.50 P(a,b) <k,l||a,d>_abab t2_abab(c,d,k,l) l2_aaaa(i,j,c,b)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j")  = l2["aaaa"]("L,i,j,c,b") * reused_["61_aa_vv"]("c,a");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,b,a,i,j");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmal2_aaaa += -0.50 P(a,b) <l,k||d,a>_aaaa t2_aaaa(d,c,l,k) l2_aaaa(i,j,c,b)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j")  = 0.50 * l2["aaaa"]("L,i,j,c,b") * reused_["62_aa_vv"]("a,c");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,b,a,i,j");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmal2_aaaa += +0.50 <i,j||k,l>_aaaa l2_aaaa(k,l,a,b)
    // flops: o2v2L1 += o4v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_aaaa("L,a,b,i,j") += 0.50 * l2["aaaa"]("L,k,l,a,b") * eri["aaaa_oooo"]("i,j,k,l");

    // sigmal2_aaaa += +0.25 <i,j||c,d>_aaaa t2_aaaa(c,d,k,l) l2_aaaa(k,l,a,b)
    // flops: o2v2L1 += o4v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_aaaa("L,a,b,i,j") += 0.25 * l2["aaaa"]("L,k,l,a,b") * reused_["26_aaaa_oooo"]("k,l,i,j");

    // sigmal2_aaaa += -1.00 P(i,j) P(a,b) <j,c||a,k>_abab l2_abab(i,k,b,c)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j")  = l2["abab"]("L,i,k,b,c") * eri["baab_vovo"]("c,j,a,k");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,a,b,j,i");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,b,a,i,j");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,b,a,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmal2_aaaa += +1.00 P(i,j) P(a,b) <j,l||d,a>_aaaa t2_abab(d,c,l,k) l2_abab(i,k,b,c)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j")  = l2["abab"]("L,i,k,b,c") * reused_["5_bbaa_voov"]("c,k,j,a");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,a,b,j,i");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,b,a,i,j");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,b,a,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmal2_aaaa += +1.00 P(i,j) P(a,b) <j,l||a,d>_abab t2_bbbb(d,c,k,l) l2_abab(i,k,b,c)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j")  = l2["abab"]("L,i,k,b,c") * reused_["15_bbaa_voov"]("c,k,j,a");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,a,b,j,i");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,b,a,i,j");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,b,a,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmal2_aaaa += +1.00 P(i,j) P(a,b) <j,c||a,k>_aaaa l2_aaaa(i,k,c,b)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j")  = l2["aaaa"]("L,i,k,c,b") * eri["aaaa_vovo"]("c,j,a,k");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,a,b,j,i");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,b,a,i,j");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,b,a,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmal2_aaaa += +1.00 P(i,j) P(a,b) <j,l||d,a>_aaaa t2_aaaa(d,c,k,l) l2_aaaa(i,k,c,b)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j")  = l2["aaaa"]("L,i,k,c,b") * reused_["13_aaaa_voov"]("c,k,j,a");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,a,b,j,i");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,b,a,i,j");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,b,a,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmal2_aaaa += +0.50 <d,c||a,b>_aaaa l2_aaaa(i,j,d,c)
    // flops: o2v2L1 += o2v4L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_aaaa("L,a,b,i,j") -= 0.50 * eri["aaaa_vvvv"]("c,d,a,b") * l2["aaaa"]("L,i,j,d,c");

    // sigmal2_aaaa += +0.25 <l,k||a,b>_aaaa t2_aaaa(d,c,l,k) l2_aaaa(i,j,d,c)
    // flops: o2v2L1 += o4v2L1 o4v2L1
    //  mems: o2v2L1 += o4v0L1 o2v2L1
    sigmal2_aaaa("L,a,b,i,j") -= 0.25 * l2["aaaa"]("L,i,j,d,c") * t2["aaaa"]("d,c,l,k") * eri["aaaa_oovv"]("k,l,a,b");

    // sigmal2_aaaa += +1.00 P(i,j) P(a,b) <j,l||a,d>_abab t2_abab(c,d,k,l) l2_aaaa(i,k,c,b)
    // flops: o2v2L1 += o3v3L1 o3v3L1
    //  mems: o2v2L1 += o2v2L1 o2v2L1
    tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j")  = l2["aaaa"]("L,i,k,c,b") * t2["abab"]("c,d,k,l") * eri["abab_oovv"]("j,l,a,d");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,a,b,i,j");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,a,b,j,i");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["perm_aaaa_Lvvoo"]("L,b,a,i,j");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["perm_aaaa_Lvvoo"]("L,b,a,j,i");
    tmps_["perm_aaaa_Lvvoo"].~TArrayD();

    // sigmal2_abab += +1.00 f_bb(j,b) l1_aa(i,a)
    // flops: o2v2L1 += o2v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") += l1["aa"]("L,i,a") * f["bb_ov"]("j,b");

    // sigmal2_abab += +0.25 <l,k||c,d>_bbbb t2_bbbb(c,d,l,k) l2_abab(i,j,a,b)
    //              += +1.00 f_aa(k,k) l2_abab(i,j,a,b)
    //              += +1.00 f_bb(k,k) l2_abab(i,j,a,b)
    //              += -0.50 <l,k||l,k>_aaaa l2_abab(i,j,a,b)
    //              += -0.50 <l,k||l,k>_bbbb l2_abab(i,j,a,b)
    //              += +0.25 <l,k||c,d>_aaaa t2_aaaa(c,d,l,k) l2_abab(i,j,a,b)
    //              += -0.50 <l,k||l,k>_abab l2_abab(i,j,a,b)
    //              += -0.50 <k,l||k,l>_abab l2_abab(i,j,a,b)
    //              += +0.25 <l,k||c,d>_abab t2_abab(c,d,l,k) l2_abab(i,j,a,b)
    //              += +0.25 <k,l||c,d>_abab t2_abab(c,d,k,l) l2_abab(i,j,a,b)
    //              += +0.25 <l,k||d,c>_abab t2_abab(d,c,l,k) l2_abab(i,j,a,b)
    //              += +0.25 <k,l||d,c>_abab t2_abab(d,c,k,l) l2_abab(i,j,a,b)
    // flops: o2v2L1 += o2v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= 0.25 * scalars_["9"] * l2["abab"]("L,i,j,a,b");

    // sigmal2_abab += -1.00 <i,j||k,b>_abab l1_aa(k,a)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") += l1["aa"]("L,k,a") * eri["abba_oovo"]("i,j,b,k");

    // sigmal2_abab += -1.00 f_aa(i,k) l2_abab(k,j,a,b)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= f["aa_oo"]("i,k") * l2["abab"]("L,k,j,a,b");

    // sigmal2_abab += -0.50 <i,l||c,d>_aaaa t2_aaaa(c,d,k,l) l2_abab(k,j,a,b)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= 0.50 * l2["abab"]("L,k,j,a,b") * reused_["52_aa_oo"]("i,k");

    // sigmal2_abab += -0.50 <i,l||c,d>_abab t2_abab(c,d,k,l) l2_abab(k,j,a,b)
    //              += -0.50 <i,l||d,c>_abab t2_abab(d,c,k,l) l2_abab(k,j,a,b)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= l2["abab"]("L,k,j,a,b") * reused_["92_aa_oo"]("k,i");

    // sigmal2_abab += -1.00 f_bb(j,k) l2_abab(i,k,a,b)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= l2["abab"]("L,i,k,a,b") * f["bb_oo"]("j,k");

    // sigmal2_abab += -1.00 <i,j||a,k>_abab l1_bb(k,b)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= l1["bb"]("L,k,b") * eri["abab_oovo"]("i,j,a,k");

    // sigmal2_abab += -0.50 <l,j||c,d>_abab t2_abab(c,d,l,k) l2_abab(i,k,a,b)
    //              += -0.50 <l,j||d,c>_abab t2_abab(d,c,l,k) l2_abab(i,k,a,b)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= l2["abab"]("L,i,k,a,b") * reused_["50_bb_oo"]("k,j");

    // sigmal2_abab += -0.50 <j,l||c,d>_bbbb t2_bbbb(c,d,k,l) l2_abab(i,k,a,b)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= 0.50 * l2["abab"]("L,i,k,a,b") * reused_["87_bb_oo"]("j,k");

    // sigmal2_abab += +1.00 f_bb(c,b) l2_abab(i,j,a,c)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") += f["bb_vv"]("c,b") * l2["abab"]("L,i,j,a,c");

    // sigmal2_abab += +1.00 <i,c||a,b>_abab l1_bb(j,c)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= eri["baab_vovv"]("c,i,a,b") * l1["bb"]("L,j,c");

    // sigmal2_abab += -0.50 <l,k||d,b>_bbbb t2_bbbb(d,c,l,k) l2_abab(i,j,a,c)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= 0.50 * l2["abab"]("L,i,j,a,c") * reused_["55_bb_vv"]("b,c");

    // sigmal2_abab += -0.50 <l,k||d,b>_abab t2_abab(d,c,l,k) l2_abab(i,j,a,c)
    //              += -0.50 <k,l||d,b>_abab t2_abab(d,c,k,l) l2_abab(i,j,a,c)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= l2["abab"]("L,i,j,a,c") * reused_["59_bb_vv"]("b,c");

    // sigmal2_abab += +1.00 f_aa(c,a) l2_abab(i,j,c,b)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") += f["aa_vv"]("c,a") * l2["abab"]("L,i,j,c,b");

    // sigmal2_abab += +1.00 <c,j||a,b>_abab l1_aa(i,c)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") += l1["aa"]("L,i,c") * eri["abab_vovv"]("c,j,a,b");

    // sigmal2_abab += -0.50 <l,k||a,d>_abab t2_abab(c,d,l,k) l2_abab(i,j,c,b)
    //              += -0.50 <k,l||a,d>_abab t2_abab(c,d,k,l) l2_abab(i,j,c,b)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= l2["abab"]("L,i,j,c,b") * reused_["61_aa_vv"]("c,a");

    // sigmal2_abab += -0.50 <l,k||d,a>_aaaa t2_aaaa(d,c,l,k) l2_abab(i,j,c,b)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= 0.50 * l2["abab"]("L,i,j,c,b") * reused_["62_aa_vv"]("a,c");

    // sigmal2_abab += +0.50 <i,j||k,l>_abab l2_abab(k,l,a,b)
    //              += +0.50 <i,j||l,k>_abab l2_abab(l,k,a,b)
    // flops: o2v2L1 += o4v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") += l2["abab"]("L,k,l,a,b") * eri["abab_oooo"]("i,j,k,l");

    // sigmal2_abab += +0.25 <i,j||c,d>_abab t2_abab(c,d,k,l) l2_abab(k,l,a,b)
    //              += +0.25 <i,j||c,d>_abab t2_abab(c,d,l,k) l2_abab(l,k,a,b)
    //              += +0.25 <i,j||d,c>_abab t2_abab(d,c,k,l) l2_abab(k,l,a,b)
    //              += +0.25 <i,j||d,c>_abab t2_abab(d,c,l,k) l2_abab(l,k,a,b)
    // flops: o2v2L1 += o4v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") += l2["abab"]("L,k,l,a,b") * reused_["30_abab_oooo"]("i,j,k,l");

    // sigmal2_abab += -1.00 <i,c||k,b>_abab l2_abab(k,j,a,c)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= eri["baba_vovo"]("c,i,b,k") * l2["abab"]("L,k,j,a,c");

    // sigmal2_abab += +1.00 <i,l||d,b>_abab t2_abab(d,c,k,l) l2_abab(k,j,a,c)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") += l2["abab"]("L,k,j,a,c") * reused_["44_baab_voov"]("c,k,i,b");

    // sigmal2_abab += +1.00 <j,c||b,k>_bbbb l2_abab(i,k,a,c)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= l2["abab"]("L,i,k,a,c") * eri["bbbb_vovo"]("c,j,b,k");

    // sigmal2_abab += -1.00 <i,c||a,k>_abab l2_bbbb(j,k,c,b)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") += eri["baab_vovo"]("c,i,a,k") * l2["bbbb"]("L,j,k,c,b");

    // sigmal2_abab += +1.00 <l,j||d,b>_abab t2_abab(d,c,l,k) l2_abab(i,k,a,c)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") += l2["abab"]("L,i,k,a,c") * reused_["4_bbbb_voov"]("c,k,j,b");

    // sigmal2_abab += +1.00 <i,l||d,a>_aaaa t2_abab(d,c,l,k) l2_bbbb(j,k,c,b)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= l2["bbbb"]("L,j,k,c,b") * reused_["5_bbaa_voov"]("c,k,i,a");

    // sigmal2_abab += +1.00 <i,l||a,d>_abab t2_bbbb(d,c,k,l) l2_bbbb(j,k,c,b)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") += l2["bbbb"]("L,j,k,c,b") * reused_["15_bbaa_voov"]("c,k,i,a");

    // sigmal2_abab += +1.00 <j,l||d,b>_bbbb t2_bbbb(d,c,k,l) l2_abab(i,k,a,c)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= l2["abab"]("L,i,k,a,c") * reused_["41_bbbb_voov"]("c,k,j,b");

    // sigmal2_abab += +1.00 <i,c||a,k>_aaaa l2_abab(k,j,c,b)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= eri["aaaa_vovo"]("c,i,a,k") * l2["abab"]("L,k,j,c,b");

    // sigmal2_abab += -1.00 <c,j||k,b>_abab l2_aaaa(i,k,c,a)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") += l2["aaaa"]("L,i,k,c,a") * eri["abba_vovo"]("c,j,b,k");

    // sigmal2_abab += +1.00 <l,j||d,b>_abab t2_aaaa(d,c,k,l) l2_aaaa(i,k,c,a)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") += l2["aaaa"]("L,i,k,c,a") * reused_["7_aabb_voov"]("c,k,j,b");

    // sigmal2_abab += +1.00 <i,l||d,a>_aaaa t2_aaaa(d,c,k,l) l2_abab(k,j,c,b)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= l2["abab"]("L,k,j,c,b") * reused_["13_aaaa_voov"]("c,k,i,a");

    // sigmal2_abab += +1.00 <j,l||d,b>_bbbb t2_abab(c,d,k,l) l2_aaaa(i,k,c,a)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= l2["aaaa"]("L,i,k,c,a") * reused_["46_aabb_voov"]("c,k,j,b");

    // sigmal2_abab += -1.00 <c,j||a,k>_abab l2_abab(i,k,c,b)
    // flops: o2v2L1 += o3v3L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") -= l2["abab"]("L,i,k,c,b") * eri["abab_vovo"]("c,j,a,k");

    // sigmal2_abab += +0.50 <d,c||a,b>_abab l2_abab(i,j,d,c)
    //              += +0.50 <c,d||a,b>_abab l2_abab(i,j,c,d)
    // flops: o2v2L1 += o2v4L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_abab("L,a,b,i,j") += eri["abab_vvvv"]("d,c,a,b") * l2["abab"]("L,i,j,d,c");

    // sigmal2_abab += +0.25 <l,k||a,b>_abab t2_abab(d,c,l,k) l2_abab(i,j,d,c)
    //              += +0.25 <k,l||a,b>_abab t2_abab(d,c,k,l) l2_abab(i,j,d,c)
    //              += +0.25 <l,k||a,b>_abab t2_abab(c,d,l,k) l2_abab(i,j,c,d)
    //              += +0.25 <k,l||a,b>_abab t2_abab(c,d,k,l) l2_abab(i,j,c,d)
    // flops: o2v2L1 += o4v2L1 o4v2L1
    //  mems: o2v2L1 += o4v0L1 o2v2L1
    sigmal2_abab("L,a,b,i,j") += l2["abab"]("L,i,j,d,c") * t2["abab"]("d,c,l,k") * eri["abab_oovv"]("l,k,a,b");

    // sigmal2_abab += +1.00 <l,j||a,d>_abab t2_abab(c,d,l,k) l2_abab(i,k,c,b)
    // flops: o2v2L1 += o3v3L1 o3v3L1
    //  mems: o2v2L1 += o2v2L1 o2v2L1
    sigmal2_abab("L,a,b,i,j") += l2["abab"]("L,i,k,c,b") * t2["abab"]("c,d,l,k") * eri["abab_oovv"]("l,j,a,d");

    // sigmal2_abab += +1.00 <i,l||a,d>_abab t2_abab(c,d,k,l) l2_abab(k,j,c,b)
    // flops: o2v2L1 += o3v3L1 o3v3L1
    //  mems: o2v2L1 += o2v2L1 o2v2L1
    sigmal2_abab("L,a,b,i,j") += l2["abab"]("L,k,j,c,b") * t2["abab"]("c,d,k,l") * eri["abab_oovv"]("i,l,a,d");

    // sigmal2_bbbb += +0.25 <l,k||c,d>_bbbb t2_bbbb(c,d,l,k) l2_bbbb(i,j,a,b)
    //              += +1.00 f_aa(k,k) l2_bbbb(i,j,a,b)
    //              += +1.00 f_bb(k,k) l2_bbbb(i,j,a,b)
    //              += -0.50 <l,k||l,k>_aaaa l2_bbbb(i,j,a,b)
    //              += -0.50 <l,k||l,k>_bbbb l2_bbbb(i,j,a,b)
    //              += +0.25 <l,k||c,d>_aaaa t2_aaaa(c,d,l,k) l2_bbbb(i,j,a,b)
    //              += -0.50 <l,k||l,k>_abab l2_bbbb(i,j,a,b)
    //              += -0.50 <k,l||k,l>_abab l2_bbbb(i,j,a,b)
    //              += +0.25 <l,k||c,d>_abab t2_abab(c,d,l,k) l2_bbbb(i,j,a,b)
    //              += +0.25 <k,l||c,d>_abab t2_abab(c,d,k,l) l2_bbbb(i,j,a,b)
    //              += +0.25 <l,k||d,c>_abab t2_abab(d,c,l,k) l2_bbbb(i,j,a,b)
    //              += +0.25 <k,l||d,c>_abab t2_abab(d,c,k,l) l2_bbbb(i,j,a,b)
    // flops: o2v2L1 += o2v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_bbbb("L,a,b,i,j") -= 0.25 * scalars_["9"] * l2["bbbb"]("L,i,j,a,b");

    // sigmal2_bbbb += -1.00 P(a,b) <i,j||a,k>_bbbb l1_bb(k,b)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j")  = l1["bb"]("L,k,b") * eri["bbbb_oovo"]("i,j,a,k");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,b,a,i,j");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmal2_bbbb += -1.00 P(i,j) f_bb(j,k) l2_bbbb(i,k,a,b)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j")  = l2["bbbb"]("L,i,k,a,b") * f["bb_oo"]("j,k");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,a,b,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmal2_bbbb += -0.50 P(i,j) <l,j||c,d>_abab t2_abab(c,d,l,k) l2_bbbb(i,k,a,b)
    //              += -0.50 P(i,j) <l,j||d,c>_abab t2_abab(d,c,l,k) l2_bbbb(i,k,a,b)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j")  = l2["bbbb"]("L,i,k,a,b") * reused_["50_bb_oo"]("k,j");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,a,b,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmal2_bbbb += -0.50 P(i,j) <j,l||c,d>_bbbb t2_bbbb(c,d,k,l) l2_bbbb(i,k,a,b)
    // flops: o2v2L1 += o3v2L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j")  = 0.50 * l2["bbbb"]("L,i,k,a,b") * reused_["87_bb_oo"]("j,k");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,a,b,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmal2_bbbb += -1.00 P(i,j) <j,c||a,b>_bbbb l1_bb(i,c)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j")  = l1["bb"]("L,i,c") * eri["bbbb_vovv"]("c,j,a,b");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,a,b,j,i");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmal2_bbbb += +1.00 P(a,b) f_bb(c,a) l2_bbbb(i,j,c,b)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j")  = f["bb_vv"]("c,a") * l2["bbbb"]("L,i,j,c,b");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,b,a,i,j");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmal2_bbbb += -0.50 P(a,b) <l,k||d,a>_bbbb t2_bbbb(d,c,l,k) l2_bbbb(i,j,c,b)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j")  = 0.50 * l2["bbbb"]("L,i,j,c,b") * reused_["55_bb_vv"]("a,c");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,b,a,i,j");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmal2_bbbb += -0.50 P(a,b) <l,k||d,a>_abab t2_abab(d,c,l,k) l2_bbbb(i,j,c,b)
    //              += -0.50 P(a,b) <k,l||d,a>_abab t2_abab(d,c,k,l) l2_bbbb(i,j,c,b)
    // flops: o2v2L1 += o2v3L1
    //  mems: o2v2L1 += o2v2L1
    tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j")  = l2["bbbb"]("L,i,j,c,b") * reused_["59_bb_vv"]("a,c");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["perm_bbbb_Lvvoo"]("L,a,b,i,j");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["perm_bbbb_Lvvoo"]("L,b,a,i,j");
    tmps_["perm_bbbb_Lvvoo"].~TArrayD();

    // sigmal2_bbbb += +0.50 <i,j||k,l>_bbbb l2_bbbb(k,l,a,b)
    // flops: o2v2L1 += o4v2L1
    //  mems: o2v2L1 += o2v2L1
    sigmal2_bbbb("L,a,b,i,j") += 0.50 * l2["bbbb"]("L,k,l,a,b") * eri["bbbb_oooo"]("i,j,k,l");
}