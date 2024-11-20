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

#include "cc_cavity/include/qed_ccsd_21/eom_ee_qed_ccsd_21.h"
#include "cc_cavity/include/qed_ccsd_21/qed_ccsd_21.h"

void hilbert::EOM_EE_QED_CCSD_21::sigma_ee_21_2() {

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

    /// unpack left/right operators

    TArrayMap r1, r2, r1_1, r2_1,
            l1, l2, l1_1, l2_1;

    TArrayD &r0 = evec_blks_["r0"];
    r1["aa"] = evec_blks_["r1_aa"];
    r1["bb"] = evec_blks_["r1_bb"];
    r2["aaaa"] = evec_blks_["r2_aaaa"];
    r2["abab"] = evec_blks_["r2_abab"];
    r2["bbbb"] = evec_blks_["r2_bbbb"];

    TArrayD &r0_1 = evec_blks_["r0_1"];
    r1_1["aa"] = evec_blks_["r1_1_aa"];
    r1_1["bb"] = evec_blks_["r1_1_bb"];
    r2_1["aaaa"] = evec_blks_["r2_1_aaaa"];
    r2_1["abab"] = evec_blks_["r2_1_abab"];
    r2_1["bbbb"] = evec_blks_["r2_1_bbbb"];

    TArrayD &l0 = evec_blks_["l0"];
    l1["aa"] = evec_blks_["l1_aa"];
    l1["bb"] = evec_blks_["l1_bb"];
    l2["aaaa"] = evec_blks_["l2_aaaa"];
    l2["abab"] = evec_blks_["l2_abab"];
    l2["bbbb"] = evec_blks_["l2_bbbb"];

    TArrayD &l0_1 = evec_blks_["l0_1"];
    l1_1["aa"] = evec_blks_["l1_1_aa"];
    l1_1["bb"] = evec_blks_["l1_1_bb"];
    l2_1["aaaa"] = evec_blks_["l2_1_aaaa"];
    l2_1["abab"] = evec_blks_["l2_1_abab"];
    l2_1["bbbb"] = evec_blks_["l2_1_bbbb"];

    /// extract the sigma vectors

    TArrayD &sigmar0 = sigvec_blks_["sigmar0"];
    TArrayD &sigmal0 = sigvec_blks_["sigmal0"];

    TArrayD &sigmar1_aa = sigvec_blks_["sigmar1_aa"];
    TArrayD &sigmar1_bb = sigvec_blks_["sigmar1_bb"];
    TArrayD &sigmal1_aa = sigvec_blks_["sigmal1_aa"];
    TArrayD &sigmal1_bb = sigvec_blks_["sigmal1_bb"];

    TArrayD &sigmar2_aaaa = sigvec_blks_["sigmar2_aaaa"];
    TArrayD &sigmar2_abab = sigvec_blks_["sigmar2_abab"];
    TArrayD &sigmar2_bbbb = sigvec_blks_["sigmar2_bbbb"];
    TArrayD &sigmal2_aaaa = sigvec_blks_["sigmal2_aaaa"];
    TArrayD &sigmal2_abab = sigvec_blks_["sigmal2_abab"];
    TArrayD &sigmal2_bbbb = sigvec_blks_["sigmal2_bbbb"];

    TArrayD &sigmar0_1      = sigvec_blks_["sigmar0_1"];
    TArrayD &sigmal0_1      = sigvec_blks_["sigmal0_1"];
    TArrayD &sigmar1_1_aa   = sigvec_blks_["sigmar1_1_aa"];
    TArrayD &sigmar1_1_bb   = sigvec_blks_["sigmar1_1_bb"];
    TArrayD &sigmal1_1_aa   = sigvec_blks_["sigmal1_1_aa"];
    TArrayD &sigmal1_1_bb   = sigvec_blks_["sigmal1_1_bb"];
    TArrayD &sigmar2_1_aaaa = sigvec_blks_["sigmar2_1_aaaa"];
    TArrayD &sigmar2_1_abab = sigvec_blks_["sigmar2_1_abab"];
    TArrayD &sigmar2_1_bbbb = sigvec_blks_["sigmar2_1_bbbb"];
    TArrayD &sigmal2_1_aaaa = sigvec_blks_["sigmal2_1_aaaa"];
    TArrayD &sigmal2_1_abab = sigvec_blks_["sigmal2_1_abab"];
    TArrayD &sigmal2_1_bbbb = sigvec_blks_["sigmal2_1_bbbb"];


    /// extract coherent state basis operators

    TArrayD &csigmar0        = tmps_["csigmar0"];
    TArrayD &csigmal0        = tmps_["csigmal0"];
    TArrayD &csigmar0_1      = tmps_["csigmar0_1"];
    TArrayD &csigmal0_1      = tmps_["csigmal0_1"];

    TArrayD &csigmar1_aa     = tmps_["csigmar1_aa"];
    TArrayD &csigmar1_bb     = tmps_["csigmar1_bb"];
    TArrayD &csigmar1_1_aa   = tmps_["csigmar1_1_aa"];
    TArrayD &csigmar1_1_bb   = tmps_["csigmar1_1_bb"];
    TArrayD &csigmal1_aa     = tmps_["csigmal1_aa"];
    TArrayD &csigmal1_bb     = tmps_["csigmal1_bb"];
    TArrayD &csigmal1_1_aa   = tmps_["csigmal1_1_aa"];
    TArrayD &csigmal1_1_bb   = tmps_["csigmal1_1_bb"];

    TArrayD &csigmar2_aaaa   = tmps_["csigmar2_aaaa"];
    TArrayD &csigmar2_abab   = tmps_["csigmar2_abab"];
    TArrayD &csigmar2_bbbb   = tmps_["csigmar2_bbbb"];
    TArrayD &csigmar2_1_aaaa = tmps_["csigmar2_1_aaaa"];
    TArrayD &csigmar2_1_abab = tmps_["csigmar2_1_abab"];
    TArrayD &csigmar2_1_bbbb = tmps_["csigmar2_1_bbbb"];
    TArrayD &csigmal2_aaaa   = tmps_["csigmal2_aaaa"];
    TArrayD &csigmal2_abab   = tmps_["csigmal2_abab"];
    TArrayD &csigmal2_bbbb   = tmps_["csigmal2_bbbb"];
    TArrayD &csigmal2_1_aaaa = tmps_["csigmal2_1_aaaa"];
    TArrayD &csigmal2_1_abab = tmps_["csigmal2_1_abab"];
    TArrayD &csigmal2_1_bbbb = tmps_["csigmal2_1_bbbb"];


    // flops: o2v2L1  = o3v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["67_bbbb_Lvvoo"]("R,a,b,i,j")  = -1.00 * t1_1["bb"]("a,k") * tmps_["63_bbbb_Lvooo"]("R,b,i,j,k");
    tmps_["67_bbbb_Lvvoo"].~TArrayD();
    tmps_["63_bbbb_Lvooo"].~TArrayD();

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["68_bb_Lov"]("R,i,a")  = eri["abab_oovv"]("j,i,b,a") * r1["aa"]("R,b,j");

    // flops: o1v1L1  = o1v1L1 o2v1L1 o1v2L1 o3v2L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v2L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o2v3L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v2L1 o1v1L1 o1v2L1 o1v1L1 o2v1L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o2v3L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    tmps_["69_aa_Lov"]("R,i,a")  = -1.00 * dp["aa_vo"]("a,i") * r0_1("R");
    tmps_["69_aa_Lov"]("R,i,a") -= r1["aa"]("R,a,m") * f["aa_oo"]("m,i");
    tmps_["69_aa_Lov"]("R,i,a") -= r1["aa"]("R,c,i") * reused_["57_aa_vv"]("c,a");
    tmps_["69_aa_Lov"]("R,i,a") += r2["abab"]("R,a,b,k,j") * eri["abba_oovo"]("k,j,b,i");
    tmps_["69_aa_Lov"]("R,i,a") -= scalars_["2"] * r1_1["aa"]("R,a,i");
    tmps_["69_aa_Lov"]("R,i,a") += f["bb_ov"]("j,b") * r2["abab"]("R,a,b,i,j");
    tmps_["69_aa_Lov"]("R,i,a") += 0.50 * r1["aa"]("R,c,i") * reused_["23_aa_vv"]("a,c");
    tmps_["69_aa_Lov"]("R,i,a") -= scalars_["5"] * r1["aa"]("R,a,i");
    tmps_["69_aa_Lov"]("R,i,a") += r0_1("R") * reused_["54_aa_vo"]("a,i");
    tmps_["69_aa_Lov"]("R,i,a") -= eri["abba_vovo"]("a,j,b,i") * r1["bb"]("R,b,j");
    tmps_["69_aa_Lov"]("R,i,a") += r2["aaaa"]("R,e,a,i,m") * reused_["13_aa_ov"]("m,e");
    tmps_["69_aa_Lov"]("R,i,a") -= r1["bb"]("R,d,l") * reused_["6_aabb_voov"]("a,i,l,d");
    tmps_["69_aa_Lov"]("R,i,a") -= eri["aaaa_vovo"]("a,m,e,i") * r1["aa"]("R,e,m");
    tmps_["69_aa_Lov"]("R,i,a") -= r0_1("R") * reused_["58_aa_vo"]("a,i");
    tmps_["69_aa_Lov"]("R,i,a") -= r1["aa"]("R,c,k") * reused_["55_aaaa_voov"]("a,i,k,c");
    tmps_["69_aa_Lov"]("R,i,a") -= f["aa_ov"]("m,e") * r2["aaaa"]("R,e,a,i,m");
    tmps_["69_aa_Lov"]("R,i,a") += 0.50 * r1["aa"]("R,a,k") * reused_["15_aa_oo"]("i,k");
    tmps_["69_aa_Lov"]("R,i,a") += eri["abab_vovv"]("a,j,e,d") * r2["abab"]("R,e,d,i,j");
    tmps_["69_aa_Lov"]("R,i,a") -= scalars_["1"] * r1_1["aa"]("R,a,i");
    tmps_["69_aa_Lov"]("R,i,a") -= scalars_["3"] * r1["aa"]("R,a,i");
    tmps_["69_aa_Lov"]("R,i,a") += f["aa_vv"]("a,e") * r1["aa"]("R,e,i");
    tmps_["69_aa_Lov"]("R,i,a") -= r1["aa"]("R,e,i") * reused_["59_aa_vv"]("a,e");
    tmps_["69_aa_Lov"]("R,i,a") -= r1["aa"]("R,a,k") * reused_["16_aa_oo"]("k,i");
    tmps_["69_aa_Lov"]("R,i,a") -= scalars_["4"] * r1["aa"]("R,a,i");
    tmps_["69_aa_Lov"]("R,i,a") -= r2["abab"]("R,a,b,i,j") * reused_["12_bb_ov"]("j,b");
    tmps_["69_aa_Lov"]("R,i,a") += 0.50 * eri["aaaa_oovo"]("m,k,e,i") * r2["aaaa"]("R,e,a,k,m");
    tmps_["69_aa_Lov"]("R,i,a") += 0.50 * eri["aaaa_vovv"]("a,m,e,c") * r2["aaaa"]("R,e,c,i,m");
    tmps_["69_aa_Lov"]("R,i,a") += r1["aa"]("R,a,m") * reused_["17_aa_oo"]("m,i");
    tmps_["69_aa_Lov"]("R,i,a") += r1["bb"]("R,d,l") * reused_["56_aabb_voov"]("a,i,l,d");
    tmps_["69_aa_Lov"]("R,i,a") -= t1_1["aa"]("a,i") * tmps_["45_L"]("R");
    tmps_["69_aa_Lov"]("R,i,a") += t2["abab"]("a,b,i,j") * tmps_["68_bb_Lov"]("R,j,b");

    // sigmar1_aa  = -0.50 <k,j||i,b>_abab r2_abab(a,b,k,j)
    //            += -0.50 <j,k||i,b>_abab r2_abab(a,b,j,k)
    //            += -0.50 <k,j||c,b>_abab r1_aa(c,i) t2_abab(a,b,k,j)
    //            += -0.50 <j,k||c,b>_abab r1_aa(c,i) t2_abab(a,b,j,k)
    //            += -1.00 f_aa(j,i) r1_aa(a,j)
    //            += -1.00 d-_bb(j,j) r1_1_aa(a,i)
    //            += +1.00 f_bb(j,b) r2_abab(a,b,i,j)
    //            += -1.00 d-_aa(a,i) r0_1
    //            += -0.50 <k,j||b,c>_aaaa r1_aa(c,i) t2_aaaa(b,a,k,j)
    //            += -1.00 d-_aa(j,j) r1_aa(a,i) t0_1
    //              += +1.00 f_aa(k,k) l2_aaaa(i,j,a,b)
    //              += -0.50 <l,k||l,k>_abab l2_aaaa(i,j,a,b)
    //              += -0.50 <k,l||k,l>_abab l2_aaaa(i,j,a,b)
    //              += +0.25 <l,k||c,d>_aaaa t2_aaaa(c,d,l,k) l2_aaaa(i,j,a,b)
    //              += -0.50 <l,k||l,k>_bbbb l2_aaaa(i,j,a,b)
    //              += -0.50 <l,k||l,k>_aaaa l2_aaaa(i,j,a,b)
    //              += +0.25 <l,k||c,d>_bbbb t2_bbbb(c,d,l,k) l2_aaaa(i,j,a,b)
    //              += +1.00 f_bb(k,k) l2_aaaa(i,j,a,b)
    //              += -1.00 d-_bb(k,k) t0_1 l2_aaaa(i,j,a,b)
    //              += +0.25 <l,k||c,d>_abab t2_abab(c,d,l,k) l2_aaaa(i,j,a,b)
    //              += +0.25 <k,l||c,d>_abab t2_abab(c,d,k,l) l2_aaaa(i,j,a,b)
    //              += +0.25 <l,k||d,c>_abab t2_abab(d,c,l,k) l2_aaaa(i,j,a,b)
    //              += +0.25 <k,l||d,c>_abab t2_abab(d,c,k,l) l2_aaaa(i,j,a,b)
    //            += +1.00 d-_aa(j,b) r0_1 t2_aaaa(b,a,i,j)
    //            += +1.00 <a,j||i,b>_abab r1_bb(b,j)
    //            += +1.00 d-_aa(j,b) r2_aaaa(b,a,i,j) t0_1
    //            += -1.00 <j,k||b,c>_abab r1_bb(c,k) t2_aaaa(b,a,i,j)
    //            += +1.00 <j,a||b,i>_aaaa r1_aa(b,j)
    //            += -1.00 d-_bb(j,b) r0_1 t2_abab(a,b,i,j)
    //            += +1.00 <k,j||b,c>_aaaa r1_aa(c,k) t2_aaaa(b,a,i,j)
    //            += -1.00 f_aa(j,b) r2_aaaa(b,a,i,j)
    //            += -0.50 <k,j||b,c>_aaaa r1_aa(a,k) t2_aaaa(b,c,i,j)
    //            += +0.50 <a,j||b,c>_abab r2_abab(b,c,i,j)
    //            += +0.50 <a,j||c,b>_abab r2_abab(c,b,i,j)
    //            += -1.00 d-_aa(j,j) r1_1_aa(a,i)
    //            += -1.00 d-_bb(j,b) r1_aa(a,i) t1_1_bb(b,j)
    //            += +1.00 f_aa(a,b) r1_aa(b,i)
    //            += -1.00 d-_aa(a,b) r1_aa(b,i) t0_1
    //            += -0.50 <k,j||b,c>_abab r1_aa(a,k) t2_abab(b,c,i,j)
    //            += -0.50 <k,j||c,b>_abab r1_aa(a,k) t2_abab(c,b,i,j)
    //            += -1.00 d-_aa(j,b) r1_aa(a,i) t1_1_aa(b,j)
    //            += -1.00 d-_bb(j,b) r2_abab(a,b,i,j) t0_1
    //            += -0.50 <k,j||b,i>_aaaa r2_aaaa(b,a,k,j)
    //            += -0.50 <j,a||b,c>_aaaa r2_aaaa(b,c,i,j)
    //            += +1.00 d-_aa(j,i) r1_aa(a,j) t0_1
    //            += -1.00 <k,j||b,c>_bbbb r1_bb(c,k) t2_abab(a,b,i,j)
    //            += -1.00 d-_bb(j,b) r1_bb(b,j) t1_1_aa(a,i)
    //            += -1.00 d-_aa(j,b) r1_aa(b,j) t1_1_aa(a,i)
    //            += +1.00 <k,j||c,b>_abab r1_aa(c,k) t2_abab(a,b,i,j)
    sigmar1_aa("R,a,i")  = tmps_["69_aa_Lov"]("R,i,a");
    tmps_["69_aa_Lov"].~TArrayD();

    // flops: o2v0L1  = o2v1L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["70_bb_Loo"]("R,i,j")  = dp["bb_ov"]("i,a") * r1["bb"]("R,a,j");
    tmps_["70_bb_Loo"].~TArrayD();

    // flops: o2v2L1  = o2v1L1 o3v2L1
    //  mems: o2v2L1  = o2v0L1 o2v2L1
    tmps_["71_bbbb_Lovvo"]("R,i,a,b,j")  = dp["bb_ov"]("k,c") * r1["bb"]("R,c,i") * t2["bbbb"]("a,b,j,k");

    // sigmar2_1_bbbb += -1.00 P(i,j) d+_bb(k,c) r1_bb(c,i) t2_bbbb(a,b,j,k)
    sigmar2_1_bbbb("R,a,b,i,j") -= tmps_["71_bbbb_Lovvo"]("R,i,a,b,j");
    sigmar2_1_bbbb("R,a,b,i,j") += tmps_["71_bbbb_Lovvo"]("R,j,a,b,i");

    // flops: o2v0L1  = o2v1L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["73_bb_Loo"]("R,i,j")  = r1["bb"]("R,a,i") * dp["bb_ov"]("j,a");

    // flops: o2v2L1  = o3v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["72_bbbb_Lvvoo"]("R,a,b,i,j")  = 0.50 * r2["bbbb"]("R,a,b,i,k") * reused_["31_bb_oo"]("j,k");

    // flops: o2v2L1  = o3v2L1 o2v3L1 o3v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["74_bbbb_Lvvoo"]("R,a,b,i,j")  = -1.00 * r2["bbbb"]("R,a,b,i,k") * f["bb_oo"]("k,j");
    tmps_["74_bbbb_Lvvoo"]("R,a,b,i,j") -= 0.50 * r1["bb"]("R,e,i") * reused_["41_bbbb_vovv"]("e,j,a,b");
    tmps_["74_bbbb_Lvvoo"]("R,a,b,i,j") -= r2["bbbb"]("R,a,b,i,m") * reused_["30_bb_oo"]("m,j");
    tmps_["74_bbbb_Lvvoo"]("R,a,b,i,j") += r2["bbbb"]("R,a,b,i,k") * reused_["35_bb_oo"]("k,j");
    tmps_["74_bbbb_Lvvoo"]("R,a,b,i,j") += r2_1["bbbb"]("R,a,b,i,k") * dp["bb_oo"]("k,j");
    tmps_["74_bbbb_Lvvoo"]("R,a,b,i,j") += eri["bbbb_vvvo"]("a,b,e,j") * r1["bb"]("R,e,i");
    tmps_["74_bbbb_Lvvoo"]("R,a,b,i,j") += tmps_["72_bbbb_Lvvoo"]("R,a,b,i,j");
    tmps_["74_bbbb_Lvvoo"]("R,a,b,i,j") -= t2_1["bbbb"]("a,b,j,k") * tmps_["73_bb_Loo"]("R,i,k");
    tmps_["74_bbbb_Lvvoo"]("R,a,b,i,j") -= t2["bbbb"]("a,b,j,k") * tmps_["49_bb_Loo"]("R,i,k");
    tmps_["74_bbbb_Lvvoo"]("R,a,b,i,j") += t2["bbbb"]("a,b,j,k") * tmps_["56_bb_Loo"]("R,i,k");
    tmps_["74_bbbb_Lvvoo"]("R,a,b,i,j") -= t0_1 * tmps_["71_bbbb_Lovvo"]("R,i,a,b,j");
    tmps_["74_bbbb_Lvvoo"]("R,a,b,i,j") -= t2["bbbb"]("a,b,i,k") * tmps_["54_bb_Loo"]("R,k,j");
    tmps_["72_bbbb_Lvvoo"].~TArrayD();
    tmps_["71_bbbb_Lovvo"].~TArrayD();

    // sigmar2_bbbb += +1.00 P(i,j) d-_bb(k,j) r2_bbbb(a,b,i,k) t0_1
    //              += -0.50 P(i,j) <k,l||c,d>_abab r2_bbbb(a,b,i,l) t2_abab(c,d,k,j)
    //              += -0.50 P(i,j) <k,l||d,c>_abab r2_bbbb(a,b,i,l) t2_abab(d,c,k,j)
    //              += +1.00 P(i,j) d-_bb(k,j) r2_1_bbbb(a,b,i,k)
    //              += +1.00 P(i,j) <a,b||c,j>_bbbb r1_bb(c,i)
    //              += +0.50 P(i,j) <l,k||c,j>_bbbb r1_bb(c,i) t2_bbbb(a,b,l,k)
    //              += -1.00 P(i,j) f_bb(k,j) r2_bbbb(a,b,i,k)
    //              += -0.50 P(i,j) <l,k||c,d>_bbbb r2_bbbb(a,b,i,l) t2_bbbb(c,d,j,k)
    //              += -1.00 P(i,j) d-_bb(k,c) r1_bb(c,i) t2_1_bbbb(a,b,j,k)
    //              += -1.00 P(i,j) d-_bb(k,c) r1_1_bb(c,i) t2_bbbb(a,b,j,k)
    //              += +1.00 P(i,j) f_bb(k,c) r1_bb(c,i) t2_bbbb(a,b,j,k)
    //              += -0.50 P(i,j) <l,k||c,d>_bbbb r2_bbbb(c,d,i,l) t2_bbbb(a,b,j,k)
    //              += +0.50 P(i,j) <l,k||c,d>_abab r2_abab(c,d,l,i) t2_bbbb(a,b,j,k)
    //              += +0.50 P(i,j) <l,k||d,c>_abab r2_abab(d,c,l,i) t2_bbbb(a,b,j,k)
    //              += -1.00 P(i,j) d-_bb(k,c) r1_bb(c,i) t2_bbbb(a,b,j,k) t0_1
    //              += -1.00 P(i,j) <l,k||c,j>_abab r1_aa(c,l) t2_bbbb(a,b,i,k)
    //              += -1.00 P(i,j) <l,k||c,j>_bbbb r1_bb(c,l) t2_bbbb(a,b,i,k)
    sigmar2_bbbb("R,a,b,i,j") += tmps_["74_bbbb_Lvvoo"]("R,a,b,i,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["74_bbbb_Lvvoo"]("R,a,b,j,i");
    tmps_["74_bbbb_Lvvoo"].~TArrayD();

    // flops: o3v1L1  = o3v2L1
    //  mems: o3v1L1  = o3v1L1
    tmps_["78_aaaa_Loovo"]("L,i,j,a,k")  = l2_1["aaaa"]("L,i,j,a,b") * t1_1["aa"]("b,k");

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["77_aaaa_Loovv"]("L,i,j,a,b")  = -4.00 * scalars_["1"] * l2_1["aaaa"]("L,i,j,a,b");

    // flops: o2v2L1  = o3v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["75_bbbb_Lvvoo"]("R,a,b,i,j")  = -1.00 * t2["bbbb"]("a,b,i,k") * tmps_["54_bb_Loo"]("R,k,j");
    tmps_["75_bbbb_Lvvoo"].~TArrayD();

    // flops: o4v0L1  = o4v2L1 o4v2L1 o4v0L1
    //  mems: o4v0L1  = o4v0L1 o4v0L1 o4v0L1
    tmps_["76_aaaa_Loooo"]("L,i,j,k,l")  = l2["aaaa"]("L,i,j,a,b") * t2["aaaa"]("a,b,k,l");
    tmps_["76_aaaa_Loooo"]("L,i,j,k,l") += l2_1["aaaa"]("L,i,j,a,b") * t2_1["aaaa"]("a,b,k,l");

    // flops: o2v2L1  = o4v2L1 o2v4L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["79_aaaa_Lvvoo"]("L,a,b,i,j")  = -1.00 * eri["aaaa_oovv"]("k,l,a,b") * tmps_["76_aaaa_Loooo"]("L,i,j,l,k");
    tmps_["79_aaaa_Lvvoo"]("L,a,b,i,j") -= 2.00 * eri["aaaa_vvvv"]("d,c,a,b") * l2["aaaa"]("L,i,j,c,d");
    tmps_["79_aaaa_Lvvoo"]("L,a,b,i,j") -= 4.00 * scalars_["4"] * l2["aaaa"]("L,i,j,a,b");
    tmps_["79_aaaa_Lvvoo"]("L,a,b,i,j") += 4.00 * scalars_["6"] * l2_1["aaaa"]("L,i,j,a,b");
    tmps_["79_aaaa_Lvvoo"]("L,a,b,i,j") += l2_1["aaaa"]("L,k,l,a,b") * reused_["22_aaaa_oooo"]("k,l,i,j");
    tmps_["79_aaaa_Lvvoo"]("L,a,b,i,j") += 4.00 * l2_1["aaaa"]("L,k,l,a,b") * reused_["60_aaaa_oooo"]("i,j,l,k");
    tmps_["79_aaaa_Lvvoo"]("L,a,b,i,j") -= 4.00 * scalars_["3"] * l2["aaaa"]("L,i,j,a,b");
    tmps_["79_aaaa_Lvvoo"]("L,a,b,i,j") += 2.00 * eri["aaaa_oooo"]("i,j,k,l") * l2["aaaa"]("L,k,l,a,b");
    tmps_["79_aaaa_Lvvoo"]("L,a,b,i,j") -= 4.00 * scalars_["5"] * l2["aaaa"]("L,i,j,a,b");
    tmps_["79_aaaa_Lvvoo"]("L,a,b,i,j") -= 4.00 * scalars_["2"] * l2_1["aaaa"]("L,i,j,a,b");
    tmps_["79_aaaa_Lvvoo"]("L,a,b,i,j") += l2["aaaa"]("L,k,l,a,b") * reused_["21_aaaa_oooo"]("i,j,k,l");
    tmps_["79_aaaa_Lvvoo"]("L,a,b,i,j") += tmps_["77_aaaa_Loovv"]("L,i,j,a,b");
    tmps_["79_aaaa_Lvvoo"]("L,a,b,i,j") -= 4.00 * eri["aaaa_vovv"]("c,k,a,b") * tmps_["78_aaaa_Loovo"]("L,i,j,c,k");
    tmps_["77_aaaa_Loovv"].~TArrayD();
    tmps_["76_aaaa_Loooo"].~TArrayD();

    // sigmal2_aaaa  = +0.25 <i,j||c,d>_aaaa t2_1_aaaa(c,d,k,l) l2_1_aaaa(k,l,a,b)
    //              += +1.00 t0_1 l2_1_aaaa(i,j,a,b) w0
    //              += +0.25 <l,k||c,d>_aaaa t2_1_aaaa(c,d,l,k) l2_1_aaaa(i,j,a,b)
    //              += +0.25 <l,k||c,d>_abab t2_1_abab(c,d,l,k) l2_1_aaaa(i,j,a,b)
    //              += +0.25 <k,l||c,d>_abab t2_1_abab(c,d,k,l) l2_1_aaaa(i,j,a,b)
    //              += +0.25 <l,k||d,c>_abab t2_1_abab(d,c,l,k) l2_1_aaaa(i,j,a,b)
    //              += +0.25 <k,l||d,c>_abab t2_1_abab(d,c,k,l) l2_1_aaaa(i,j,a,b)
    //              += -1.00 d-_aa(k,c) t0_1 t1_1_aa(c,k) l2_1_aaaa(i,j,a,b)
    //              += +1.00 f_bb(k,c) t1_1_bb(c,k) l2_1_aaaa(i,j,a,b)
    //              += -1.00 d-_bb(k,c) t0_1 t1_1_bb(c,k) l2_1_aaaa(i,j,a,b)
    //              += +0.25 <l,k||c,d>_bbbb t2_1_bbbb(c,d,l,k) l2_1_aaaa(i,j,a,b)
    //              += +1.00 f_aa(k,c) t1_1_aa(c,k) l2_1_aaaa(i,j,a,b)
    //              += -1.00 d-_aa(k,c) t1_1_aa(c,k) l2_aaaa(i,j,a,b)
    //              += +0.50 <d,c||a,b>_aaaa l2_aaaa(i,j,d,c)
    //              += +1.00 <i,j||c,l>_aaaa t1_1_aa(c,k) l2_1_aaaa(k,l,a,b)
    //              += -1.00 d-_bb(k,c) t1_1_bb(c,k) l2_aaaa(i,j,a,b)
    //              += +0.50 <i,j||k,l>_aaaa l2_aaaa(k,l,a,b)
    //              += -1.00 d-_aa(k,k) t0_1 l2_aaaa(i,j,a,b)
    //            += +1.00 f_aa(j,j) r1_aa(a,i)
    //            += -0.50 <k,j||k,j>_abab r1_aa(a,i)
    //            += -0.50 <j,k||j,k>_abab r1_aa(a,i)
    //            += +0.25 <k,j||b,c>_aaaa r1_aa(a,i) t2_aaaa(b,c,k,j)
    //            += -0.50 <k,j||k,j>_bbbb r1_aa(a,i)
    //            += -0.50 <k,j||k,j>_aaaa r1_aa(a,i)
    //            += +0.25 <k,j||b,c>_bbbb r1_aa(a,i) t2_bbbb(b,c,k,j)
    //            += +1.00 f_bb(j,j) r1_aa(a,i)
    //            += -1.00 d-_bb(j,j) r1_aa(a,i) t0_1
    //            += +0.25 <k,j||b,c>_abab r1_aa(a,i) t2_abab(b,c,k,j)
    //            += +0.25 <j,k||b,c>_abab r1_aa(a,i) t2_abab(b,c,j,k)
    //            += +0.25 <k,j||c,b>_abab r1_aa(a,i) t2_abab(c,b,k,j)
    //            += +0.25 <j,k||c,b>_abab r1_aa(a,i) t2_abab(c,b,j,k)
    //              += -1.00 d+_bb(k,k) l2_1_aaaa(i,j,a,b)
    //              += +0.25 <i,j||c,d>_aaaa t2_aaaa(c,d,k,l) l2_aaaa(k,l,a,b)
    //              += -1.00 d+_aa(k,k) l2_1_aaaa(i,j,a,b)
    //              += +0.25 <l,k||a,b>_aaaa t2_1_aaaa(d,c,l,k) l2_1_aaaa(i,j,d,c)
    //              += +0.25 <l,k||a,b>_aaaa t2_aaaa(d,c,l,k) l2_aaaa(i,j,d,c)
    //              += +1.00 <k,d||a,b>_aaaa t1_1_aa(c,k) l2_1_aaaa(i,j,d,c)
    sigmal2_aaaa("L,a,b,i,j")  = 0.25 * tmps_["79_aaaa_Lvvoo"]("L,a,b,i,j");
    tmps_["79_aaaa_Lvvoo"].~TArrayD();

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["81_aaaa_Loovv"]("L,i,j,a,b")  = -2.00 * scalars_["3"] * l2_1["aaaa"]("L,i,j,a,b");

    // flops: o4v0L1  = o4v2L1
    //  mems: o4v0L1  = o4v0L1
    tmps_["80_aaaa_Loooo"]("L,i,j,k,l")  = l2_1["aaaa"]("L,i,j,a,b") * t2["aaaa"]("a,b,k,l");

    // flops: o2v2L1  = o4v2L1 o2v4L1 0 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 0 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["82_aaaa_Loovv"]("L,i,j,a,b")  = -0.25 * eri["aaaa_oovv"]("k,l,a,b") * tmps_["80_aaaa_Loooo"]("L,i,j,l,k");
    tmps_["82_aaaa_Loovv"]("L,i,j,a,b") -= 0.50 * eri["aaaa_vvvv"]("d,c,a,b") * l2_1["aaaa"]("L,i,j,c,d");
    tmps_["82_aaaa_Loovv"]("L,i,j,a,b") += (w0 + -1.00 * scalars_["5"]) * l2_1["aaaa"]("L,i,j,a,b");
    tmps_["82_aaaa_Loovv"]("L,i,j,a,b") += eri["aaaa_oovv"]("i,j,a,b") * l0_1("L");
    tmps_["82_aaaa_Loovv"]("L,i,j,a,b") -= scalars_["2"] * l2["aaaa"]("L,i,j,a,b");
    tmps_["82_aaaa_Loovv"]("L,i,j,a,b") -= scalars_["1"] * l2["aaaa"]("L,i,j,a,b");
    tmps_["82_aaaa_Loovv"]("L,i,j,a,b") -= 2.00 * scalars_["4"] * l2_1["aaaa"]("L,i,j,a,b");
    tmps_["82_aaaa_Loovv"]("L,i,j,a,b") += 0.25 * l2_1["aaaa"]("L,k,l,a,b") * reused_["21_aaaa_oooo"]("i,j,k,l");
    tmps_["82_aaaa_Loovv"]("L,i,j,a,b") += 0.50 * l2_1["aaaa"]("L,k,l,a,b") * eri["aaaa_oooo"]("i,j,k,l");
    tmps_["82_aaaa_Loovv"]("L,i,j,a,b") += tmps_["81_aaaa_Loovv"]("L,i,j,a,b");
    tmps_["81_aaaa_Loovv"].~TArrayD();
    tmps_["80_aaaa_Loooo"].~TArrayD();

    // sigmal2_1_aaaa  = +1.00 l2_1_aaaa(i,j,a,b) w0
    //                += -1.00 d-_aa(k,k) t0_1 l2_1_aaaa(i,j,a,b)
    //                += +1.00 f_aa(k,k) r2_1_abab(a,b,i,j)
    //                += -0.50 <l,k||l,k>_abab r2_1_abab(a,b,i,j)
    //                += -0.50 <k,l||k,l>_abab r2_1_abab(a,b,i,j)
    //                += +0.25 <l,k||c,d>_aaaa r2_1_abab(a,b,i,j) t2_aaaa(c,d,l,k)
    //                += -0.50 <l,k||l,k>_bbbb r2_1_abab(a,b,i,j)
    //                += -0.50 <l,k||l,k>_aaaa r2_1_abab(a,b,i,j)
    //                += +0.25 <l,k||c,d>_bbbb r2_1_abab(a,b,i,j) t2_bbbb(c,d,l,k)
    //                += +1.00 f_bb(k,k) r2_1_abab(a,b,i,j)
    //                += -1.00 d-_bb(k,k) r2_1_abab(a,b,i,j) t0_1
    //                += +0.25 <l,k||c,d>_abab r2_1_abab(a,b,i,j) t2_abab(c,d,l,k)
    //                += +0.25 <k,l||c,d>_abab r2_1_abab(a,b,i,j) t2_abab(c,d,k,l)
    //                += +0.25 <l,k||d,c>_abab r2_1_abab(a,b,i,j) t2_abab(d,c,l,k)
    //                += +0.25 <k,l||d,c>_abab r2_1_abab(a,b,i,j) t2_abab(d,c,k,l)
    //                += +1.00 <i,j||a,b>_aaaa l0_1
    //                += +0.50 <d,c||a,b>_aaaa l2_1_aaaa(i,j,d,c)
    //                += -1.00 d-_bb(k,k) l2_aaaa(i,j,a,b)
    //                += -1.00 d-_aa(k,k) l2_aaaa(i,j,a,b)
    //                += -2.00 d-_aa(k,c) t1_1_aa(c,k) l2_1_aaaa(i,j,a,b)
    //                += +0.25 <i,j||c,d>_aaaa t2_aaaa(c,d,k,l) l2_1_aaaa(k,l,a,b)
    //                += +0.50 <i,j||k,l>_aaaa l2_1_aaaa(k,l,a,b)
    //                += -2.00 d-_bb(k,c) t1_1_bb(c,k) l2_1_aaaa(i,j,a,b)
    //                += +0.25 <l,k||a,b>_aaaa t2_aaaa(d,c,l,k) l2_1_aaaa(i,j,d,c)
    sigmal2_1_aaaa("L,a,b,i,j")  = tmps_["82_aaaa_Loovv"]("L,i,j,a,b");
    tmps_["82_aaaa_Loovv"].~TArrayD();

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["83_aaaa_Lvvoo"]("R,a,b,i,j")  = 2.00 * t2_1["aaaa"]("a,b,i,j") * tmps_["45_L"]("R");
    tmps_["83_aaaa_Lvvoo"].~TArrayD();

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["84_aaaa_Lvvoo"]("R,a,b,i,j")  = 2.00 * scalars_["4"] * r2["aaaa"]("R,a,b,i,j");

    // flops: o2v2L1  = o2v4L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["85_aaaa_Loovv"]("R,i,j,a,b")  = -1.00 * eri["aaaa_vvvv"]("a,b,c,d") * r2["aaaa"]("R,c,d,i,j");
    tmps_["85_aaaa_Loovv"]("R,i,j,a,b") += 2.00 * scalars_["3"] * r2["aaaa"]("R,a,b,i,j");
    tmps_["85_aaaa_Loovv"]("R,i,j,a,b") += eri["aaaa_oooo"]("l,k,i,j") * r2["aaaa"]("R,a,b,k,l");
    tmps_["85_aaaa_Loovv"]("R,i,j,a,b") += 2.00 * scalars_["1"] * r2_1["aaaa"]("R,a,b,i,j");
    tmps_["85_aaaa_Loovv"]("R,i,j,a,b") += 0.50 * r2["aaaa"]("R,a,b,k,l") * reused_["21_aaaa_oooo"]("l,k,i,j");
    tmps_["85_aaaa_Loovv"]("R,i,j,a,b") += 2.00 * scalars_["2"] * r2_1["aaaa"]("R,a,b,i,j");
    tmps_["85_aaaa_Loovv"]("R,i,j,a,b") += 2.00 * scalars_["5"] * r2["aaaa"]("R,a,b,i,j");
    tmps_["85_aaaa_Loovv"]("R,i,j,a,b") += tmps_["84_aaaa_Lvvoo"]("R,a,b,i,j");
    tmps_["85_aaaa_Loovv"]("R,i,j,a,b") += 0.50 * t2["aaaa"]("a,b,k,l") * tmps_["42_aaaa_Loooo"]("R,i,j,l,k");
    tmps_["85_aaaa_Loovv"]("R,i,j,a,b") += 2.00 * t2_1["aaaa"]("a,b,i,j") * tmps_["45_L"]("R");
    tmps_["84_aaaa_Lvvoo"].~TArrayD();
    tmps_["42_aaaa_Loooo"].~TArrayD();

    // sigmar2_aaaa += +0.50 <l,k||i,j>_aaaa r2_aaaa(a,b,l,k)
    //              += -1.00 d-_bb(k,c) r2_aaaa(a,b,i,j) t1_1_bb(c,k)
    //              += +0.50 <a,b||c,d>_aaaa r2_aaaa(c,d,i,j)
    //              += -1.00 d-_aa(k,k) r2_1_aaaa(a,b,i,j)
    //              += +0.25 <l,k||c,d>_aaaa r2_aaaa(a,b,l,k) t2_aaaa(c,d,i,j)
    //              += -1.00 d-_bb(k,k) r2_1_aaaa(a,b,i,j)
    //              += -1.00 d-_aa(k,k) r2_aaaa(a,b,i,j) t0_1
    //            += +1.00 f_aa(j,j) l1_aa(i,a)
    //            += -0.50 <k,j||k,j>_abab l1_aa(i,a)
    //            += -0.50 <j,k||j,k>_abab l1_aa(i,a)
    //            += +0.25 <k,j||b,c>_aaaa t2_aaaa(b,c,k,j) l1_aa(i,a)
    //            += -0.50 <k,j||k,j>_bbbb l1_aa(i,a)
    //            += -0.50 <k,j||k,j>_aaaa l1_aa(i,a)
    //            += +0.25 <k,j||b,c>_bbbb t2_bbbb(b,c,k,j) l1_aa(i,a)
    //            += +1.00 f_bb(j,j) l1_aa(i,a)
    //            += -1.00 d-_bb(j,j) t0_1 l1_aa(i,a)
    //            += +0.25 <k,j||b,c>_abab t2_abab(b,c,k,j) l1_aa(i,a)
    //            += +0.25 <j,k||b,c>_abab t2_abab(b,c,j,k) l1_aa(i,a)
    //            += +0.25 <k,j||c,b>_abab t2_abab(c,b,k,j) l1_aa(i,a)
    //            += +0.25 <j,k||c,b>_abab t2_abab(c,b,j,k) l1_aa(i,a)
    //              += -1.00 d-_aa(k,c) r2_aaaa(a,b,i,j) t1_1_aa(c,k)
    //              += +0.25 <l,k||c,d>_aaaa r2_aaaa(c,d,i,j) t2_aaaa(a,b,l,k)
    //              += -1.00 d-_bb(k,c) r1_bb(c,k) t2_1_aaaa(a,b,i,j)
    //              += -1.00 d-_aa(k,c) r1_aa(c,k) t2_1_aaaa(a,b,i,j)
    sigmar2_aaaa("R,a,b,i,j") -= 0.50 * tmps_["85_aaaa_Loovv"]("R,i,j,a,b");
    tmps_["85_aaaa_Loovv"].~TArrayD();

    // flops: o0v0L1  = o2v2L1
    //  mems: o0v0L1  = o0v0L1
    tmps_["86_L"]("L")  = 4.00 * l2_1["abab"]("L,i,j,a,b") * t2_1["abab"]("a,b,i,j");
    tmps_["86_L"].~TArrayD();

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["87_bb_Lov"]("L,i,a")  = dp["bb_vo"]("b,j") * l2_1["bbbb"]("L,i,j,b,a");

    // sigmal1_bb  = +1.00 d+_bb(b,j) l2_1_bbbb(i,j,b,a)
    sigmal1_bb("L,a,i")  = tmps_["87_bb_Lov"]("L,i,a");

    // flops: o1v1L1  = o1v1 o2v2L1
    //  mems: o1v1L1  = o1v1 o1v1L1
    tmps_["97_aa_Lov"]("L,i,a")  = (dp["aa_vo"]("b,j") + reused_["58_aa_vo"]("b,j")) * l2_1["aaaa"]("L,i,j,b,a");

    // sigmal1_aa  = +1.00 d+_aa(b,j) l2_1_aaaa(i,j,b,a)
    //            += +1.00 d+_bb(k,c) t2_abab(b,c,j,k) l2_1_aaaa(i,j,b,a)
    sigmal1_aa("L,a,i")  = tmps_["97_aa_Lov"]("L,i,a");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["100_aa_Lov"]("L,i,a")  = l2_1["aaaa"]("L,i,j,b,a") * reused_["54_aa_vo"]("b,j");

    // sigmal1_aa += -1.00 d+_aa(k,c) t2_aaaa(c,b,j,k) l2_1_aaaa(i,j,b,a)
    sigmal1_aa("L,a,i") -= tmps_["100_aa_Lov"]("L,i,a");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["101_bb_Lov"]("L,i,a")  = l2_1["abab"]("L,j,i,b,a") * t1_1["aa"]("b,j");

    // csigmal1_1_bb += +1.00 t1_1_aa(b,j) l2_1_abab(j,i,b,a)
    csigmal1_1_bb("L,a,i") += tmps_["101_bb_Lov"]("L,i,a");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["102_aa_Lov"]("L,i,a")  = l2_1["abab"]("L,i,j,a,b") * t1_1["bb"]("b,j");

    // csigmal1_1_aa += +1.00 t1_1_bb(b,j) l2_1_abab(i,j,a,b)
    csigmal1_1_aa("L,a,i") += tmps_["102_aa_Lov"]("L,i,a");

    // flops: o0v0L1  = o2v2L1 o2v2L1 o2v2L1 o0v0L1 o0v0L1
    //  mems: o0v0L1  = o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1
    tmps_["103_L"]("L")  = 4.00 * l2_1["abab"]("L,k,j,c,b") * t2_1["abab"]("c,b,k,j");
    tmps_["103_L"]("L") += l2_1["aaaa"]("L,k,l,c,d") * t2_1["aaaa"]("c,d,k,l");
    tmps_["103_L"]("L") += l2_1["bbbb"]("L,i,j,a,b") * t2_1["bbbb"]("a,b,i,j");

    // csigmal0_1 += +0.25 t2_1_bbbb(b,a,i,j) l2_1_bbbb(i,j,b,a)
    //            += +0.25 t2_1_aaaa(b,a,i,j) l2_1_aaaa(i,j,b,a)
    //            += +0.25 t2_1_abab(b,a,i,j) l2_1_abab(i,j,b,a)
    //            += +0.25 t2_1_abab(b,a,j,i) l2_1_abab(j,i,b,a)
    //            += +0.25 t2_1_abab(a,b,i,j) l2_1_abab(i,j,a,b)
    //            += +0.25 t2_1_abab(a,b,j,i) l2_1_abab(j,i,a,b)
    csigmal0_1("L") += 0.25 * tmps_["103_L"]("L");

    // flops: o1v1L1  = o1v1 o2v2L1
    //  mems: o1v1L1  = o1v1 o1v1L1
    tmps_["104_bb_Lov"]("L,i,a")  = (reused_["33_bb_vo"]("b,j") + -1.00 * reused_["32_bb_vo"]("b,j")) * l2_1["bbbb"]("L,i,j,b,a");

    // sigmal1_bb += -1.00 d+_bb(k,c) t2_bbbb(c,b,j,k) l2_1_bbbb(i,j,b,a)
    //            += +1.00 d+_aa(k,c) t2_abab(c,b,k,j) l2_1_bbbb(i,j,b,a)
    sigmal1_bb("L,a,i") -= tmps_["104_bb_Lov"]("L,i,a");

    // flops: o0v2L1  = o2v3L1
    //  mems: o0v2L1  = o0v2L1
    tmps_["99_aa_Lvv"]("L,a,b")  = l2_1["abab"]("L,i,j,a,c") * t2["abab"]("b,c,i,j");

    // flops: o0v0L1  = o1v1L1 o2v2L1 o0v0L1
    //  mems: o0v0L1  = o0v0L1 o0v0L1 o0v0L1
    tmps_["98_L"]("L")  = -2.00 * l1_1["aa"]("L,i,b") * reused_["54_aa_vo"]("b,i");
    tmps_["98_L"]("L") += l2_1["aaaa"]("L,i,j,a,b") * reused_["76_aaaa_vvoo"]("a,b,i,j");

    // flops: o0v0L1  = o2v2L1
    //  mems: o0v0L1  = o0v0L1
    tmps_["96_L"]("L")  = l2_1["aaaa"]("L,i,j,a,b") * reused_["75_aaaa_ovvo"]("j,a,b,i");

    // flops: o0v0L1  = o2v2L1
    //  mems: o0v0L1  = o0v0L1
    tmps_["95_L"]("L")  = l2_1["abab"]("L,i,j,a,b") * reused_["74_baab_vvoo"]("b,a,i,j");

    // flops: o0v0L1  = o2v2L1
    //  mems: o0v0L1  = o0v0L1
    tmps_["94_L"]("L")  = l2_1["abab"]("L,i,j,a,b") * reused_["73_abab_vvoo"]("a,b,i,j");

    // flops: o0v0L1  = o2v2L1
    //  mems: o0v0L1  = o0v0L1
    tmps_["93_L"]("L")  = l2_1["abab"]("L,i,j,a,b") * reused_["72_aabb_ovvo"]("i,a,b,j");

    // flops: o0v0L1  = o2v2L1
    //  mems: o0v0L1  = o0v0L1
    tmps_["92_L"]("L")  = l2_1["bbbb"]("L,i,j,a,b") * reused_["39_bbbb_ovvo"]("j,a,b,i");

    // flops: o0v0L1  = o2v2L1
    //  mems: o0v0L1  = o0v0L1
    tmps_["91_L"]("L")  = l2_1["bbbb"]("L,i,j,a,b") * reused_["52_bbbb_vvoo"]("a,b,i,j");

    // flops: o0v0L1  = o1v1L1
    //  mems: o0v0L1  = o0v0L1
    tmps_["90_L"]("L")  = l1_1["aa"]("L,i,a") * reused_["58_aa_vo"]("a,i");

    // flops: o0v0L1  = o1v1L1
    //  mems: o0v0L1  = o0v0L1
    tmps_["89_L"]("L")  = l1_1["bb"]("L,i,a") * reused_["32_bb_vo"]("a,i");

    // flops: o0v0L1  = o1v1L1
    //  mems: o0v0L1  = o0v0L1
    tmps_["88_L"]("L")  = l1_1["bb"]("L,i,a") * reused_["33_bb_vo"]("a,i");

    // flops: o0v0L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o1v1L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o1v1L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o1v1L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o1v1L1 o0v0L1 o2v2L1 o0v0L1 o1v1L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o0v0L1 o0v0L1 o2v2L1 o0v0L1 o1v1L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o1v1L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o1v1L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o1v1L1 o0v0L1 o2v2L1 o0v0L1 o0v0L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o0v0L1 o0v0L1 o2v2L1 o0v0L1 o1v1L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o1v1L1 o0v0L1 o0v0L1 o0v0L1 o1v1L1 o0v0L1 o1v1L1 o0v0L1 o2v2L1 o0v0L1 o1v1L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o1v1L1 o1v1L1 o0v0L1 o0v0L1 o0v2L1 o1v1L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o1v1L1 o0v0L1 o1v1L1 o0v0L1 o0v0L1 o1v1L1 o0v0L1 o1v1L1 o0v0L1 o1v1L1 o0v0L1 o1v1L1 o0v0L1 o0v0L1 o1v1L1 o0v0L1
    //  mems: o0v0L1  = o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1
    tmps_["105_L"]("L")  = -0.25 * l2_1["aaaa"]("L,j,k,c,b") * reused_["24_aaaa_vvoo"]("c,b,j,k");
    tmps_["105_L"]("L") -= 0.25 * l2_1["bbbb"]("L,l,i,a,d") * reused_["49_bbbb_ovvo"]("i,a,d,l");
    tmps_["105_L"]("L") -= 0.50 * l2_1["aaaa"]("L,j,k,c,b") * reused_["7_aaaa_vovo"]("b,j,c,k");
    tmps_["105_L"]("L") -= 0.50 * l2_1["aaaa"]("L,j,k,c,b") * reused_["7_aaaa_vovo"]("c,k,b,j");
    tmps_["105_L"]("L") += l2_1["aaaa"]("L,j,k,c,b") * reused_["97_aaaa_vvoo"]("b,c,j,k");
    tmps_["105_L"]("L") += t0_1 * tmps_["90_L"]("L");
    tmps_["105_L"]("L") += t0_1 * tmps_["89_L"]("L");
    tmps_["105_L"]("L") += 0.06250 * l2_1["aaaa"]("L,j,k,c,b") * reused_["25_aaaa_oovv"]("j,k,c,b");
    tmps_["105_L"]("L") += 0.50 * l2["bbbb"]("L,l,i,a,d") * reused_["52_bbbb_vvoo"]("a,d,l,i");
    tmps_["105_L"]("L") += 0.25 * eri["bbbb_vvoo"]("d,a,l,i") * l2_1["bbbb"]("L,l,i,a,d");
    tmps_["105_L"]("L") += 0.50 * l2_1["bbbb"]("L,l,i,a,d") * reused_["68_bbbb_ovvo"]("i,a,d,l");
    tmps_["105_L"]("L") += 0.25 * eri["aaaa_vvoo"]("b,c,j,k") * l2_1["aaaa"]("L,j,k,c,b");
    tmps_["105_L"]("L") -= 0.50 * l2["aaaa"]("L,j,k,c,b") * reused_["75_aaaa_ovvo"]("k,c,b,j");
    tmps_["105_L"]("L") -= l2_1["abab"]("L,n,i,c,d") * reused_["65_aabb_ovvo"]("n,c,d,i");
    tmps_["105_L"]("L") -= f["aa_vo"]("b,j") * l1_1["aa"]("L,j,b");
    tmps_["105_L"]("L") += l2_1["abab"]("L,j,i,c,d") * reused_["63_baba_ovvo"]("i,c,d,j");
    tmps_["105_L"]("L") -= l2_1["aaaa"]("L,j,k,c,b") * reused_["85_aaaa_vovo"]("c,k,b,j");
    tmps_["105_L"]("L") -= f["bb_vo"]("d,l") * l1_1["bb"]("L,l,d");
    tmps_["105_L"]("L") += scalars_["5"] * l0_1("L");
    tmps_["105_L"]("L") += t0_1 * tmps_["94_L"]("L");
    tmps_["105_L"]("L") += 2.00 * scalars_["4"] * l0_1("L");
    tmps_["105_L"]("L") += 0.1250 * l2_1["aaaa"]("L,j,k,c,b") * reused_["20_aaaa_vvoo"]("b,c,j,k");
    tmps_["105_L"]("L") -= 0.50 * l2_1["bbbb"]("L,i,m,a,d") * reused_["67_bbbb_ovvo"]("m,a,d,i");
    tmps_["105_L"]("L") += 0.1250 * l2_1["bbbb"]("L,l,i,a,d") * reused_["79_bbbb_vvoo"]("d,a,l,i");
    tmps_["105_L"]("L") += 0.50 * l2_1["aaaa"]("L,j,k,c,b") * reused_["5_aaaa_vvoo"]("c,b,j,k");
    tmps_["105_L"]("L") += l2_1["abab"]("L,j,i,b,a") * reused_["95_aabb_voov"]("b,j,i,a");
    tmps_["105_L"]("L") += l1_1["aa"]("L,j,b") * reused_["90_aa_vo"]("b,j");
    tmps_["105_L"]("L") -= 0.25 * l2_1["aaaa"]("L,j,k,c,b") * reused_["96_aaaa_ovvo"]("k,c,b,j");
    tmps_["105_L"]("L") -= eri["abab_vvoo"]("c,d,j,i") * l2_1["abab"]("L,j,i,c,d");
    tmps_["105_L"]("L") += dp["bb_vo"]("d,l") * l1["bb"]("L,l,d");
    tmps_["105_L"]("L") += l2_1["abab"]("L,k,l,c,d") * reused_["66_aabb_ovvo"]("k,c,d,l");
    tmps_["105_L"]("L") += l1["aa"]("L,j,b") * reused_["58_aa_vo"]("b,j");
    tmps_["105_L"]("L") += l2_1["abab"]("L,j,i,c,d") * reused_["99_abab_vvoo"]("c,d,j,i");
    tmps_["105_L"]("L") += l2_1["abab"]("L,k,l,c,d") * reused_["98_abba_vvoo"]("c,d,l,k");
    tmps_["105_L"]("L") += 0.50 * l2_1["abab"]("L,j,i,b,a") * reused_["71_aabb_vovo"]("b,j,a,i");
    tmps_["105_L"]("L") -= w0 * l0_1("L");
    tmps_["105_L"]("L") += 0.50 * l2_1["abab"]("L,k,l,c,d") * reused_["93_aabb_vovo"]("c,k,d,l");
    tmps_["105_L"]("L") -= 2.00 * l1_1["aa"]("L,j,b") * reused_["89_aa_ov"]("j,b");
    tmps_["105_L"]("L") -= 0.50 * l2_1["bbbb"]("L,l,i,a,d") * reused_["2_bbbb_vovo"]("d,l,a,i");
    tmps_["105_L"]("L") -= 0.50 * l2["bbbb"]("L,l,i,a,d") * reused_["39_bbbb_ovvo"]("i,a,d,l");
    tmps_["105_L"]("L") += l1_1["bb"]("L,l,d") * reused_["84_bb_vo"]("d,l");
    tmps_["105_L"]("L") -= l2["abab"]("L,k,l,c,d") * reused_["72_aabb_ovvo"]("k,c,d,l");
    tmps_["105_L"]("L") += 0.50 * l2["aaaa"]("L,j,k,c,b") * reused_["76_aaaa_vvoo"]("c,b,j,k");
    tmps_["105_L"]("L") -= l2_1["abab"]("L,k,m,c,d") * reused_["64_baba_ovvo"]("m,c,d,k");
    tmps_["105_L"]("L") += l2_1["bbbb"]("L,l,i,a,d") * reused_["100_bbbb_voov"]("a,l,i,d");
    tmps_["105_L"]("L") -= 0.50 * l2_1["bbbb"]("L,l,i,a,d") * reused_["2_bbbb_vovo"]("a,i,d,l");
    tmps_["105_L"]("L") += l1["bb"]("L,l,d") * reused_["32_bb_vo"]("d,l");
    tmps_["105_L"]("L") += 0.50 * t0_1 * tmps_["91_L"]("L");
    tmps_["105_L"]("L") -= 0.50 * t0_1 * tmps_["92_L"]("L");
    tmps_["105_L"]("L") -= l1["aa"]("L,j,b") * reused_["54_aa_vo"]("b,j");
    tmps_["105_L"]("L") += 0.50 * l2_1["aaaa"]("L,j,k,c,b") * reused_["62_aaaa_ovvo"]("k,c,b,j");
    tmps_["105_L"]("L") -= t0_1 * tmps_["88_L"]("L");
    tmps_["105_L"]("L") -= 0.25 * l2_1["bbbb"]("L,l,i,a,d") * reused_["92_bbbb_voov"]("d,l,i,a");
    tmps_["105_L"]("L") += l2["abab"]("L,j,i,c,d") * reused_["73_abab_vvoo"]("c,d,j,i");
    tmps_["105_L"]("L") += t0_1 * tmps_["95_L"]("L");
    tmps_["105_L"]("L") += l2_1["abab"]("L,k,l,b,a") * reused_["81_abba_vovo"]("b,l,a,k");
    tmps_["105_L"]("L") += dp["aa_vo"]("b,j") * l1["aa"]("L,j,b");
    tmps_["105_L"]("L") -= 0.50 * t0_1 * tmps_["96_L"]("L");
    tmps_["105_L"]("L") -= t0_1 * tmps_["93_L"]("L");
    tmps_["105_L"]("L") -= 0.50 * l2_1["aaaa"]("L,k,n,c,b") * reused_["61_aaaa_ovvo"]("n,c,b,k");
    tmps_["105_L"]("L") += 0.06250 * l2_1["bbbb"]("L,l,i,a,d") * reused_["94_bbbb_oovv"]("l,i,a,d");
    tmps_["105_L"]("L") -= l1["bb"]("L,l,d") * reused_["33_bb_vo"]("d,l");
    tmps_["105_L"]("L") += 2.00 * scalars_["3"] * l0_1("L");
    tmps_["105_L"]("L") += 2.00 * l1_1["bb"]("L,l,d") * reused_["88_bb_vo"]("d,l");
    tmps_["105_L"]("L") += l1_1["bb"]("L,l,d") * reused_["91_bb_vo"]("d,l");
    tmps_["105_L"]("L") += 0.50 * l2_1["bbbb"]("L,l,i,a,d") * reused_["4_bbbb_voov"]("d,l,i,a");
    tmps_["105_L"]("L") -= 0.50 * l1_1["aa"]("L,j,b") * reused_["87_aa_vo"]("b,j");
    tmps_["105_L"]("L") += l2_1["bbbb"]("L,l,i,a,d") * reused_["86_bbbb_vovo"]("a,i,d,l");
    tmps_["105_L"]("L") += l2["abab"]("L,j,i,b,a") * reused_["74_baab_vvoo"]("a,b,j,i");
    tmps_["105_L"]("L") -= reused_["54_aa_vo"]("c,k") * tmps_["102_aa_Lov"]("L,k,c");
    tmps_["105_L"]("L") -= t1_1["bb"]("d,l") * tmps_["87_bb_Lov"]("L,l,d");
    tmps_["105_L"]("L") += 0.25 * scalars_["1"] * tmps_["103_L"]("L");
    tmps_["105_L"]("L") += 0.50 * t0_1 * tmps_["98_L"]("L");
    tmps_["105_L"]("L") += reused_["57_aa_vv"]("h,c") * tmps_["99_aa_Lvv"]("L,c,h");
    tmps_["105_L"]("L") += reused_["58_aa_vo"]("c,k") * tmps_["102_aa_Lov"]("L,k,c");
    tmps_["105_L"]("L") += 0.25 * scalars_["2"] * tmps_["103_L"]("L");
    tmps_["105_L"]("L") += t1_1["aa"]("b,j") * tmps_["100_aa_Lov"]("L,j,b");
    tmps_["105_L"]("L") += reused_["32_bb_vo"]("a,i") * tmps_["101_bb_Lov"]("L,i,a");
    tmps_["105_L"]("L") += t1_1["bb"]("d,l") * tmps_["104_bb_Lov"]("L,l,d");
    tmps_["105_L"]("L") -= reused_["33_bb_vo"]("a,i") * tmps_["101_bb_Lov"]("L,i,a");
    tmps_["105_L"]("L") -= t1_1["aa"]("b,j") * tmps_["97_aa_Lov"]("L,j,b");
    tmps_["105_L"]("L") += dp["aa_vo"]("c,k") * tmps_["102_aa_Lov"]("L,k,c");
    tmps_["105_L"]("L") += dp["bb_vo"]("a,i") * tmps_["101_bb_Lov"]("L,i,a");
    tmps_["98_L"].~TArrayD();
    tmps_["97_aa_Lov"].~TArrayD();
    tmps_["96_L"].~TArrayD();
    tmps_["95_L"].~TArrayD();
    tmps_["94_L"].~TArrayD();
    tmps_["93_L"].~TArrayD();
    tmps_["92_L"].~TArrayD();
    tmps_["91_L"].~TArrayD();
    tmps_["90_L"].~TArrayD();
    tmps_["89_L"].~TArrayD();
    tmps_["88_L"].~TArrayD();
    tmps_["87_bb_Lov"].~TArrayD();

    // sigmal0_1  = -0.25 <l,k||d,c>_abab t2_abab(d,a,i,j) t2_abab(b,c,l,k) l2_1_abab(i,j,b,a)
    //           += -0.25 <k,l||d,c>_abab t2_abab(d,a,i,j) t2_abab(b,c,k,l) l2_1_abab(i,j,b,a)
    //           += -0.25 <l,k||d,c>_abab t2_abab(d,a,j,i) t2_abab(b,c,l,k) l2_1_abab(j,i,b,a)
    //           += -0.25 <k,l||d,c>_abab t2_abab(d,a,j,i) t2_abab(b,c,k,l) l2_1_abab(j,i,b,a)
    //           += -1.00 d-_bb(k,c) t2_abab(b,c,j,k) t1_1_bb(a,i) l2_1_abab(j,i,b,a)
    //           += -0.50 d-_aa(b,c) t2_aaaa(c,a,i,j) t0_1 l2_1_aaaa(i,j,b,a)
    //           += +1.00 d-_aa(j,b) t2_aaaa(b,a,i,j) t0_1 l1_1_aa(i,a)
    //           += -0.25 d-_aa(k,k) t2_1_bbbb(b,a,i,j) l2_1_bbbb(i,j,b,a)
    //           += -0.25 d-_aa(k,k) t2_1_aaaa(b,a,i,j) l2_1_aaaa(i,j,b,a)
    //           += -0.25 d-_aa(k,k) t2_1_abab(b,a,i,j) l2_1_abab(i,j,b,a)
    //           += -0.25 d-_aa(k,k) t2_1_abab(b,a,j,i) l2_1_abab(j,i,b,a)
    //           += -0.25 d-_aa(k,k) t2_1_abab(a,b,i,j) l2_1_abab(i,j,a,b)
    //           += -0.25 d-_aa(k,k) t2_1_abab(a,b,j,i) l2_1_abab(j,i,a,b)
    //           += +1.00 d-_bb(b,j) t1_1_bb(a,i) l2_1_bbbb(i,j,b,a)
    //           += -0.25 d-_bb(k,k) t2_1_bbbb(b,a,i,j) l2_1_bbbb(i,j,b,a)
    //           += -0.25 d-_bb(k,k) t2_1_aaaa(b,a,i,j) l2_1_aaaa(i,j,b,a)
    //           += -0.25 d-_bb(k,k) t2_1_abab(b,a,i,j) l2_1_abab(i,j,b,a)
    //           += -0.25 d-_bb(k,k) t2_1_abab(b,a,j,i) l2_1_abab(j,i,b,a)
    //           += -0.25 d-_bb(k,k) t2_1_abab(a,b,i,j) l2_1_abab(i,j,a,b)
    //           += -0.25 d-_bb(k,k) t2_1_abab(a,b,j,i) l2_1_abab(j,i,a,b)
    //           += -1.00 d-_aa(k,c) t2_aaaa(c,b,j,k) t1_1_aa(a,i) l2_1_aaaa(i,j,b,a)
    //           += -1.00 d-_aa(k,c) t2_abab(c,b,k,j) t1_1_aa(a,i) l2_1_abab(i,j,a,b)
    //           += +1.00 d-_aa(k,c) t2_aaaa(c,b,j,k) t1_1_bb(a,i) l2_1_abab(j,i,b,a)
    //           += -1.00 d-_bb(k,c) t2_bbbb(c,b,j,k) t1_1_bb(a,i) l2_1_bbbb(i,j,b,a)
    //           += +1.00 d-_aa(k,c) t2_abab(c,b,k,j) t1_1_bb(a,i) l2_1_bbbb(i,j,b,a)
    //           += +1.00 d-_bb(k,c) t2_bbbb(c,b,j,k) t1_1_aa(a,i) l2_1_abab(i,j,a,b)
    //           += +1.00 d-_aa(b,j) t1_1_aa(a,i) l2_1_aaaa(i,j,b,a)
    //           += +1.00 d-_bb(k,c) t2_abab(b,c,j,k) t1_1_aa(a,i) l2_1_aaaa(i,j,b,a)
    //           += -1.00 d-_aa(b,j) t1_1_bb(a,i) l2_1_abab(j,i,b,a)
    //           += -1.00 d-_bb(j,b) t2_abab(a,b,i,j) t0_1 l1_1_aa(i,a)
    //           += -1.00 d-_aa(k,c) t2_aaaa(c,b,i,j) t1_1_aa(a,k) l2_1_aaaa(i,j,b,a)
    //           += -1.00 d-_aa(b,c) t2_1_aaaa(c,a,i,j) l2_1_aaaa(i,j,b,a)
    //           += +0.50 f_aa(b,c) t2_aaaa(c,a,i,j) l2_1_aaaa(i,j,b,a)
    //           += +0.50 <l,k||d,c>_abab t2_aaaa(d,a,i,l) t2_abab(b,c,j,k) l2_1_aaaa(i,j,b,a)
    //           += -1.00 d-_aa(j,b) t2_abab(b,a,j,i) t0_1 l1_1_bb(i,a)
    //           += +0.06250 <l,k||c,d>_aaaa t2_aaaa(b,a,l,k) t2_aaaa(c,d,i,j) l2_1_aaaa(i,j,b,a)
    //           += +0.1250 <l,k||i,j>_aaaa t2_aaaa(b,a,l,k) l2_1_aaaa(i,j,b,a)
    //           += -0.50 d-_bb(b,c) t2_bbbb(c,a,i,j) l2_bbbb(i,j,b,a)
    //           += +0.50 <k,l||c,d>_abab t2_abab(a,d,i,l) t2_aaaa(c,b,j,k) l2_1_aaaa(i,j,b,a)
    //           += +0.25 <b,a||i,j>_bbbb l2_1_bbbb(i,j,b,a)
    //           += -0.50 <k,l||j,l>_bbbb t2_bbbb(b,a,i,k) l2_1_bbbb(i,j,b,a)
    //           += +0.25 <b,a||i,j>_aaaa l2_1_aaaa(i,j,b,a)
    //           += -0.25 <l,k||c,d>_bbbb t2_bbbb(b,a,i,l) t2_bbbb(c,d,j,k) l2_1_bbbb(i,j,b,a)
    //           += -0.50 f_bb(k,j) t2_bbbb(b,a,i,k) l2_1_bbbb(i,j,b,a)
    //           += -1.00 d-_bb(k,c) t2_bbbb(b,a,j,k) t1_1_bb(c,i) l2_1_bbbb(i,j,b,a)
    //           += +1.00 d-_bb(k,j) t2_1_bbbb(b,a,i,k) l2_1_bbbb(i,j,b,a)
    //           += +0.50 <l,k||c,d>_aaaa t2_abab(d,a,l,i) t2_abab(c,b,k,j) l2_1_bbbb(i,j,b,a)
    //           += +0.50 <l,k||c,d>_bbbb t2_bbbb(d,a,i,l) t2_bbbb(c,b,j,k) l2_1_bbbb(i,j,b,a)
    //           += -0.25 <k,l||c,d>_abab t2_bbbb(b,a,i,l) t2_abab(c,d,k,j) l2_1_bbbb(i,j,b,a)
    //           += -0.25 <k,l||d,c>_abab t2_bbbb(b,a,i,l) t2_abab(d,c,k,j) l2_1_bbbb(i,j,b,a)
    //           += +0.50 d-_aa(k,j) t2_aaaa(b,a,i,k) l2_aaaa(i,j,b,a)
    //           += +0.50 <l,i||k,i>_aaaa t2_abab(b,a,l,j) l2_1_abab(k,j,b,a)
    //           += +0.50 <l,i||k,i>_aaaa t2_abab(a,b,l,j) l2_1_abab(k,j,a,b)
    //           += +1.00 f_aa(a,i) l1_1_aa(i,a)
    //           += -0.50 <k,l||j,l>_bbbb t2_abab(b,a,i,k) l2_1_abab(i,j,b,a)
    //           += -0.50 <k,l||j,l>_bbbb t2_abab(a,b,i,k) l2_1_abab(i,j,a,b)
    //           += -1.00 <b,k||j,c>_abab t2_abab(a,c,i,k) l2_1_aaaa(i,j,b,a)
    //           += +1.00 <k,b||c,j>_aaaa t2_aaaa(c,a,i,k) l2_1_aaaa(i,j,b,a)
    //           += +1.00 f_bb(a,i) l1_1_bb(i,a)
    //           += -0.25 <l,k||c,d>_aaaa t2_aaaa(d,a,i,j) t2_aaaa(c,b,l,k) l2_1_aaaa(i,j,b,a)
    //           += -1.00 d-_aa(i,i) t0_1 l0_1
    //           += +1.00 f_aa(i,i) r0_1
    //           += -0.50 <j,i||j,i>_abab r0_1
    //           += -0.50 <i,j||i,j>_abab r0_1
    //           += +0.25 <j,i||a,b>_aaaa r0_1 t2_aaaa(a,b,j,i)
    //           += -0.50 <j,i||j,i>_bbbb r0_1
    //           += -0.50 <j,i||j,i>_aaaa r0_1
    //           += +0.25 <j,i||a,b>_bbbb r0_1 t2_bbbb(a,b,j,i)
    //           += +1.00 f_bb(i,i) r0_1
    //           += -1.00 d-_bb(i,i) r0_1 t0_1
    //           += +0.25 <j,i||a,b>_abab r0_1 t2_abab(a,b,j,i)
    //           += +0.25 <i,j||a,b>_abab r0_1 t2_abab(a,b,i,j)
    //           += +0.25 <j,i||b,a>_abab r0_1 t2_abab(b,a,j,i)
    //           += +0.25 <i,j||b,a>_abab r0_1 t2_abab(b,a,i,j)
    //           += -0.50 d-_aa(b,c) t2_abab(c,a,i,j) t0_1 l2_1_abab(i,j,b,a)
    //           += -0.50 d-_aa(b,c) t2_abab(c,a,j,i) t0_1 l2_1_abab(j,i,b,a)
    //           += +0.50 d-_bb(k,j) t2_abab(b,a,i,k) t0_1 l2_1_abab(i,j,b,a)
    //           += +0.50 d-_bb(k,j) t2_abab(a,b,i,k) t0_1 l2_1_abab(i,j,a,b)
    //           += -2.00 d-_aa(i,a) t1_1_aa(a,i) l0_1
    //           += +0.1250 <b,a||c,d>_aaaa t2_aaaa(c,d,i,j) l2_1_aaaa(i,j,b,a)
    //           += +0.50 <l,i||k,i>_bbbb t2_bbbb(b,a,j,l) l2_1_bbbb(j,k,b,a)
    //           += +0.1250 <b,a||c,d>_bbbb t2_bbbb(c,d,i,j) l2_1_bbbb(i,j,b,a)
    //           += -0.25 <l,k||d,c>_abab t2_aaaa(d,a,i,j) t2_abab(b,c,l,k) l2_1_aaaa(i,j,b,a)
    //           += -0.25 <k,l||d,c>_abab t2_aaaa(d,a,i,j) t2_abab(b,c,k,l) l2_1_aaaa(i,j,b,a)
    //           += -0.25 <l,k||c,d>_abab t2_abab(a,d,i,j) t2_abab(c,b,l,k) l2_1_abab(i,j,a,b)
    //           += -0.25 <k,l||c,d>_abab t2_abab(a,d,i,j) t2_abab(c,b,k,l) l2_1_abab(i,j,a,b)
    //           += -0.25 <l,k||c,d>_abab t2_abab(a,d,j,i) t2_abab(c,b,l,k) l2_1_abab(j,i,a,b)
    //           += -0.25 <k,l||c,d>_abab t2_abab(a,d,j,i) t2_abab(c,b,k,l) l2_1_abab(j,i,a,b)
    //           += +1.00 <k,b||c,j>_bbbb t2_abab(a,c,i,k) l2_1_abab(i,j,a,b)
    //           += -1.00 <k,b||c,j>_abab t2_aaaa(c,a,i,k) l2_1_abab(i,j,a,b)
    //           += +0.50 f_bb(b,c) t2_abab(a,c,i,j) l2_1_abab(i,j,a,b)
    //           += +0.50 f_bb(b,c) t2_abab(a,c,j,i) l2_1_abab(j,i,a,b)
    //           += -1.00 d-_bb(b,c) t2_1_abab(a,c,i,j) l2_1_abab(i,j,a,b)
    //           += -1.00 d-_bb(b,c) t2_1_abab(a,c,j,i) l2_1_abab(j,i,a,b)
    //           += +0.50 <l,k||d,c>_abab t2_aaaa(d,a,i,l) t2_bbbb(c,b,j,k) l2_1_abab(i,j,a,b)
    //           += +0.50 <k,l||c,d>_abab t2_bbbb(d,a,i,l) t2_aaaa(c,b,j,k) l2_1_abab(j,i,b,a)
    //           += +1.00 d-_aa(k,c) t2_abab(c,b,i,j) t1_1_aa(a,k) l2_1_abab(i,j,a,b)
    //           += +1.00 d-_aa(k,c) t2_abab(c,b,j,i) t1_1_aa(a,k) l2_1_abab(j,i,a,b)
    //           += -1.00 d-_bb(j,j) t1_1_aa(a,i) l1_1_aa(i,a)
    //           += -1.00 d-_aa(j,j) t1_1_aa(a,i) l1_1_aa(i,a)
    //           += -0.25 <l,k||c,d>_aaaa t2_aaaa(b,a,i,l) t2_aaaa(c,d,j,k) l2_1_aaaa(i,j,b,a)
    //           += +0.50 <l,k||c,d>_bbbb t2_abab(a,d,i,l) t2_abab(b,c,j,k) l2_1_aaaa(i,j,b,a)
    //           += -0.25 <l,k||c,d>_abab t2_aaaa(b,a,i,l) t2_abab(c,d,j,k) l2_1_aaaa(i,j,b,a)
    //           += -0.25 <l,k||d,c>_abab t2_aaaa(b,a,i,l) t2_abab(d,c,j,k) l2_1_aaaa(i,j,b,a)
    //           += -1.00 d-_aa(k,c) t2_aaaa(b,a,j,k) t1_1_aa(c,i) l2_1_aaaa(i,j,b,a)
    //           += +1.00 d-_aa(k,j) t2_1_aaaa(b,a,i,k) l2_1_aaaa(i,j,b,a)
    //           += -0.50 f_aa(k,j) t2_aaaa(b,a,i,k) l2_1_aaaa(i,j,b,a)
    //           += +0.50 <l,k||c,d>_aaaa t2_aaaa(d,a,i,l) t2_aaaa(c,b,j,k) l2_1_aaaa(i,j,b,a)
    //           += +0.25 <b,a||i,j>_abab l2_1_abab(i,j,b,a)
    //           += +0.25 <a,b||i,j>_abab l2_1_abab(i,j,a,b)
    //           += +0.25 <b,a||j,i>_abab l2_1_abab(j,i,b,a)
    //           += +0.25 <a,b||j,i>_abab l2_1_abab(j,i,a,b)
    //           += -1.00 d-_bb(a,i) l1_bb(i,a)
    //           += -0.50 <k,l||j,l>_aaaa t2_abab(b,a,k,i) l2_1_abab(j,i,b,a)
    //           += -0.50 <k,l||j,l>_aaaa t2_abab(a,b,k,i) l2_1_abab(j,i,a,b)
    //           += -1.00 d-_bb(j,b) t2_abab(a,b,i,j) l1_aa(i,a)
    //           += -0.25 <k,l||c,d>_abab t2_abab(b,a,i,l) t2_abab(c,d,k,j) l2_1_abab(i,j,b,a)
    //           += -0.25 <k,l||d,c>_abab t2_abab(b,a,i,l) t2_abab(d,c,k,j) l2_1_abab(i,j,b,a)
    //           += -0.25 <k,l||c,d>_abab t2_abab(a,b,i,l) t2_abab(c,d,k,j) l2_1_abab(i,j,a,b)
    //           += -0.25 <k,l||d,c>_abab t2_abab(a,b,i,l) t2_abab(d,c,k,j) l2_1_abab(i,j,a,b)
    //           += +0.1250 <l,k||i,j>_abab t2_abab(b,a,l,k) l2_1_abab(i,j,b,a)
    //           += +0.1250 <l,k||j,i>_abab t2_abab(b,a,l,k) l2_1_abab(j,i,b,a)
    //           += +0.1250 <k,l||i,j>_abab t2_abab(b,a,k,l) l2_1_abab(i,j,b,a)
    //           += +0.1250 <k,l||j,i>_abab t2_abab(b,a,k,l) l2_1_abab(j,i,b,a)
    //           += +0.1250 <l,k||i,j>_abab t2_abab(a,b,l,k) l2_1_abab(i,j,a,b)
    //           += +0.1250 <l,k||j,i>_abab t2_abab(a,b,l,k) l2_1_abab(j,i,a,b)
    //           += +0.1250 <k,l||i,j>_abab t2_abab(a,b,k,l) l2_1_abab(i,j,a,b)
    //           += +0.1250 <k,l||j,i>_abab t2_abab(a,b,k,l) l2_1_abab(j,i,a,b)
    //           += +0.1250 <b,a||c,d>_abab t2_abab(c,d,i,j) l2_1_abab(i,j,b,a)
    //           += +0.1250 <a,b||c,d>_abab t2_abab(c,d,i,j) l2_1_abab(i,j,a,b)
    //           += +0.1250 <b,a||c,d>_abab t2_abab(c,d,j,i) l2_1_abab(j,i,b,a)
    //           += +0.1250 <a,b||c,d>_abab t2_abab(c,d,j,i) l2_1_abab(j,i,a,b)
    //           += +0.1250 <b,a||d,c>_abab t2_abab(d,c,i,j) l2_1_abab(i,j,b,a)
    //           += +0.1250 <a,b||d,c>_abab t2_abab(d,c,i,j) l2_1_abab(i,j,a,b)
    //           += +0.1250 <b,a||d,c>_abab t2_abab(d,c,j,i) l2_1_abab(j,i,b,a)
    //           += +0.1250 <a,b||d,c>_abab t2_abab(d,c,j,i) l2_1_abab(j,i,a,b)
    //           += +0.50 f_aa(b,c) t2_abab(c,a,i,j) l2_1_abab(i,j,b,a)
    //           += +0.50 f_aa(b,c) t2_abab(c,a,j,i) l2_1_abab(j,i,b,a)
    //           += +1.00 d-_bb(k,j) t2_1_abab(b,a,i,k) l2_1_abab(i,j,b,a)
    //           += +1.00 d-_bb(k,j) t2_1_abab(a,b,i,k) l2_1_abab(i,j,a,b)
    //           += -0.50 f_bb(k,j) t2_abab(b,a,i,k) l2_1_abab(i,j,b,a)
    //           += -0.50 f_bb(k,j) t2_abab(a,b,i,k) l2_1_abab(i,j,a,b)
    //           += +1.00 d-_aa(k,c) t2_abab(b,a,k,j) t1_1_aa(c,i) l2_1_abab(i,j,b,a)
    //           += +1.00 d-_aa(k,c) t2_abab(a,b,k,j) t1_1_aa(c,i) l2_1_abab(i,j,a,b)
    //           += -0.25 <l,k||c,d>_bbbb t2_abab(b,a,i,l) t2_bbbb(c,d,j,k) l2_1_abab(i,j,b,a)
    //           += -0.25 <l,k||c,d>_bbbb t2_abab(a,b,i,l) t2_bbbb(c,d,j,k) l2_1_abab(i,j,a,b)
    //           += -1.00 d-_aa(b,c) t2_1_abab(c,a,i,j) l2_1_abab(i,j,b,a)
    //           += -1.00 d-_aa(b,c) t2_1_abab(c,a,j,i) l2_1_abab(j,i,b,a)
    //           += -1.00 <b,k||c,j>_abab t2_abab(c,a,i,k) l2_1_abab(i,j,b,a)
    //           += +0.50 <k,l||d,c>_abab t2_abab(d,a,i,l) t2_abab(b,c,k,j) l2_1_abab(i,j,b,a)
    //           += +0.50 <l,k||c,d>_abab t2_abab(a,d,l,i) t2_abab(c,b,j,k) l2_1_abab(j,i,a,b)
    //           += +1.00 d-_bb(k,c) t2_abab(b,c,i,j) t1_1_bb(a,k) l2_1_abab(i,j,b,a)
    //           += +1.00 d-_bb(k,c) t2_abab(b,c,j,i) t1_1_bb(a,k) l2_1_abab(j,i,b,a)
    //           += +0.06250 <l,k||c,d>_abab t2_abab(b,a,l,k) t2_abab(c,d,i,j) l2_1_abab(i,j,b,a)
    //           += +0.06250 <l,k||c,d>_abab t2_abab(b,a,l,k) t2_abab(c,d,j,i) l2_1_abab(j,i,b,a)
    //           += +0.06250 <l,k||d,c>_abab t2_abab(b,a,l,k) t2_abab(d,c,i,j) l2_1_abab(i,j,b,a)
    //           += +0.06250 <l,k||d,c>_abab t2_abab(b,a,l,k) t2_abab(d,c,j,i) l2_1_abab(j,i,b,a)
    //           += +0.06250 <k,l||c,d>_abab t2_abab(b,a,k,l) t2_abab(c,d,i,j) l2_1_abab(i,j,b,a)
    //           += +0.06250 <k,l||c,d>_abab t2_abab(b,a,k,l) t2_abab(c,d,j,i) l2_1_abab(j,i,b,a)
    //           += +0.06250 <k,l||d,c>_abab t2_abab(b,a,k,l) t2_abab(d,c,i,j) l2_1_abab(i,j,b,a)
    //           += +0.06250 <k,l||d,c>_abab t2_abab(b,a,k,l) t2_abab(d,c,j,i) l2_1_abab(j,i,b,a)
    //           += +0.06250 <l,k||c,d>_abab t2_abab(a,b,l,k) t2_abab(c,d,i,j) l2_1_abab(i,j,a,b)
    //           += +0.06250 <l,k||c,d>_abab t2_abab(a,b,l,k) t2_abab(c,d,j,i) l2_1_abab(j,i,a,b)
    //           += +0.06250 <l,k||d,c>_abab t2_abab(a,b,l,k) t2_abab(d,c,i,j) l2_1_abab(i,j,a,b)
    //           += +0.06250 <l,k||d,c>_abab t2_abab(a,b,l,k) t2_abab(d,c,j,i) l2_1_abab(j,i,a,b)
    //           += +0.06250 <k,l||c,d>_abab t2_abab(a,b,k,l) t2_abab(c,d,i,j) l2_1_abab(i,j,a,b)
    //           += +0.06250 <k,l||c,d>_abab t2_abab(a,b,k,l) t2_abab(c,d,j,i) l2_1_abab(j,i,a,b)
    //           += +0.06250 <k,l||d,c>_abab t2_abab(a,b,k,l) t2_abab(d,c,i,j) l2_1_abab(i,j,a,b)
    //           += +0.06250 <k,l||d,c>_abab t2_abab(a,b,k,l) t2_abab(d,c,j,i) l2_1_abab(j,i,a,b)
    //           += -0.25 <l,k||c,d>_aaaa t2_abab(d,a,i,j) t2_aaaa(c,b,l,k) l2_1_abab(i,j,b,a)
    //           += -0.25 <l,k||c,d>_aaaa t2_abab(d,a,j,i) t2_aaaa(c,b,l,k) l2_1_abab(j,i,b,a)
    //           += -0.25 <l,k||c,d>_abab t2_abab(b,a,l,i) t2_abab(c,d,j,k) l2_1_abab(j,i,b,a)
    //           += -0.25 <l,k||d,c>_abab t2_abab(b,a,l,i) t2_abab(d,c,j,k) l2_1_abab(j,i,b,a)
    //           += -0.25 <l,k||c,d>_abab t2_abab(a,b,l,i) t2_abab(c,d,j,k) l2_1_abab(j,i,a,b)
    //           += -0.25 <l,k||d,c>_abab t2_abab(a,b,l,i) t2_abab(d,c,j,k) l2_1_abab(j,i,a,b)
    //           += +1.00 <k,b||c,j>_aaaa t2_abab(c,a,k,i) l2_1_abab(j,i,b,a)
    //           += +1.00 d-_aa(k,j) t2_1_abab(b,a,k,i) l2_1_abab(j,i,b,a)
    //           += +1.00 d-_aa(k,j) t2_1_abab(a,b,k,i) l2_1_abab(j,i,a,b)
    //           += +1.00 d-_bb(k,c) t2_abab(b,a,j,k) t1_1_bb(c,i) l2_1_abab(j,i,b,a)
    //           += +1.00 d-_bb(k,c) t2_abab(a,b,j,k) t1_1_bb(c,i) l2_1_abab(j,i,a,b)
    //           += -1.00 <b,k||j,c>_abab t2_bbbb(c,a,i,k) l2_1_abab(j,i,b,a)
    //           += -0.50 f_aa(k,j) t2_abab(b,a,k,i) l2_1_abab(j,i,b,a)
    //           += -0.50 f_aa(k,j) t2_abab(a,b,k,i) l2_1_abab(j,i,a,b)
    //           += +0.50 <l,k||d,c>_abab t2_abab(d,a,l,i) t2_abab(b,c,j,k) l2_1_abab(j,i,b,a)
    //           += +0.50 <k,l||c,d>_abab t2_abab(a,d,i,l) t2_abab(c,b,k,j) l2_1_abab(i,j,a,b)
    //           += -0.25 <l,k||c,d>_aaaa t2_abab(b,a,l,i) t2_aaaa(c,d,j,k) l2_1_abab(j,i,b,a)
    //           += -0.25 <l,k||c,d>_aaaa t2_abab(a,b,l,i) t2_aaaa(c,d,j,k) l2_1_abab(j,i,a,b)
    //           += +0.50 <l,k||c,d>_aaaa t2_aaaa(d,a,i,l) t2_abab(c,b,k,j) l2_1_abab(i,j,a,b)
    //           += -0.25 <l,k||c,d>_bbbb t2_abab(a,d,i,j) t2_bbbb(c,b,l,k) l2_1_abab(i,j,a,b)
    //           += -0.25 <l,k||c,d>_bbbb t2_abab(a,d,j,i) t2_bbbb(c,b,l,k) l2_1_abab(j,i,a,b)
    //           += +0.50 <l,k||c,d>_bbbb t2_abab(a,d,i,l) t2_bbbb(c,b,j,k) l2_1_abab(i,j,a,b)
    //           += +1.00 l0_1 w0
    //           += +0.50 <l,k||c,d>_bbbb t2_bbbb(d,a,i,l) t2_abab(b,c,j,k) l2_1_abab(j,i,b,a)
    //           += +0.50 <l,k||c,d>_aaaa t2_abab(d,a,l,i) t2_aaaa(c,b,j,k) l2_1_abab(j,i,b,a)
    //           += +2.00 d-_aa(j,i) t1_1_aa(a,j) l1_1_aa(i,a)
    //           += +2.00 d-_aa(j,b) t2_1_aaaa(b,a,i,j) l1_1_aa(i,a)
    //           += -2.00 d-_aa(a,b) t1_1_aa(b,i) l1_1_aa(i,a)
    //           += -2.00 d-_bb(j,b) t2_1_abab(a,b,i,j) l1_1_aa(i,a)
    //           += +0.50 <l,k||d,c>_abab t2_abab(d,a,l,i) t2_bbbb(c,b,j,k) l2_1_bbbb(i,j,b,a)
    //           += +0.50 d-_bb(k,j) t2_bbbb(b,a,i,k) l2_bbbb(i,j,b,a)
    //           += -1.00 f_bb(j,b) t2_bbbb(b,a,i,j) l1_1_bb(i,a)
    //           += -1.00 d-_bb(a,i) t0_1 l1_1_bb(i,a)
    //           += -0.50 <k,j||b,i>_abab t2_abab(b,a,k,j) l1_1_bb(i,a)
    //           += -0.50 <j,k||b,i>_abab t2_abab(b,a,j,k) l1_1_bb(i,a)
    //           += -0.50 <k,j||b,i>_bbbb t2_bbbb(b,a,k,j) l1_1_bb(i,a)
    //           += +1.00 f_aa(j,b) t2_abab(b,a,j,i) l1_1_bb(i,a)
    //           += +0.50 <j,a||b,c>_abab t2_abab(b,c,j,i) l1_1_bb(i,a)
    //           += +0.50 <j,a||c,b>_abab t2_abab(c,b,j,i) l1_1_bb(i,a)
    //           += -0.50 <j,a||b,c>_bbbb t2_bbbb(b,c,i,j) l1_1_bb(i,a)
    //           += +0.50 d-_aa(k,j) t2_abab(b,a,k,i) l2_abab(j,i,b,a)
    //           += +0.50 d-_aa(k,j) t2_abab(a,b,k,i) l2_abab(j,i,a,b)
    //           += -0.50 d-_aa(b,c) t2_aaaa(c,a,i,j) l2_aaaa(i,j,b,a)
    //           += +0.50 <l,i||k,i>_bbbb t2_abab(b,a,j,l) l2_1_abab(j,k,b,a)
    //           += +0.50 <l,i||k,i>_bbbb t2_abab(a,b,j,l) l2_1_abab(j,k,a,b)
    //           += -1.00 d-_bb(k,c) t2_bbbb(c,b,i,j) t1_1_bb(a,k) l2_1_bbbb(i,j,b,a)
    //           += -1.00 d-_bb(b,c) t2_1_bbbb(c,a,i,j) l2_1_bbbb(i,j,b,a)
    //           += +0.50 f_bb(b,c) t2_bbbb(c,a,i,j) l2_1_bbbb(i,j,b,a)
    //           += +0.50 <k,l||c,d>_abab t2_bbbb(d,a,i,l) t2_abab(c,b,k,j) l2_1_bbbb(i,j,b,a)
    //           += -1.00 d-_aa(j,b) t2_abab(b,a,j,i) l1_bb(i,a)
    //           += -0.50 d-_bb(b,c) t2_bbbb(c,a,i,j) t0_1 l2_1_bbbb(i,j,b,a)
    //           += +0.50 d-_bb(k,j) t2_bbbb(b,a,i,k) t0_1 l2_1_bbbb(i,j,b,a)
    //           += +1.00 d-_aa(j,b) t2_aaaa(b,a,i,j) l1_aa(i,a)
    //           += -0.50 <k,l||j,l>_aaaa t2_aaaa(b,a,i,k) l2_1_aaaa(i,j,b,a)
    //           += +1.00 d-_bb(j,b) t2_bbbb(b,a,i,j) t0_1 l1_1_bb(i,a)
    //           += -0.25 <l,k||c,d>_bbbb t2_bbbb(d,a,i,j) t2_bbbb(c,b,l,k) l2_1_bbbb(i,j,b,a)
    //           += -0.50 d-_aa(b,c) t2_abab(c,a,i,j) l2_abab(i,j,b,a)
    //           += -0.50 d-_aa(b,c) t2_abab(c,a,j,i) l2_abab(j,i,b,a)
    //           += +0.50 d-_bb(k,j) t2_abab(b,a,i,k) l2_abab(i,j,b,a)
    //           += +0.50 d-_bb(k,j) t2_abab(a,b,i,k) l2_abab(i,j,a,b)
    //           += -0.50 d-_bb(b,c) t2_abab(a,c,i,j) t0_1 l2_1_abab(i,j,a,b)
    //           += -0.50 d-_bb(b,c) t2_abab(a,c,j,i) t0_1 l2_1_abab(j,i,a,b)
    //           += -1.00 <k,b||j,c>_abab t2_abab(a,c,k,i) l2_1_abab(j,i,a,b)
    //           += -1.00 d-_aa(a,i) l1_aa(i,a)
    //           += +0.50 d-_aa(k,j) t2_aaaa(b,a,i,k) t0_1 l2_1_aaaa(i,j,b,a)
    //           += +0.50 d-_aa(k,j) t2_abab(b,a,k,i) t0_1 l2_1_abab(j,i,b,a)
    //           += +0.50 d-_aa(k,j) t2_abab(a,b,k,i) t0_1 l2_1_abab(j,i,a,b)
    //           += +0.50 <l,i||k,i>_aaaa t2_aaaa(b,a,j,l) l2_1_aaaa(j,k,b,a)
    //           += +0.06250 <l,k||c,d>_bbbb t2_bbbb(b,a,l,k) t2_bbbb(c,d,i,j) l2_1_bbbb(i,j,b,a)
    //           += +0.1250 <l,k||i,j>_bbbb t2_bbbb(b,a,l,k) l2_1_bbbb(i,j,b,a)
    //           += +1.00 d-_bb(j,b) t2_bbbb(b,a,i,j) l1_bb(i,a)
    //           += -2.00 d-_bb(i,a) t1_1_bb(a,i) l0_1
    //           += -2.00 d-_bb(a,b) t1_1_bb(b,i) l1_1_bb(i,a)
    //           += +2.00 d-_bb(j,i) t1_1_bb(a,j) l1_1_bb(i,a)
    //           += -2.00 d-_aa(j,b) t2_1_abab(b,a,j,i) l1_1_bb(i,a)
    //           += +2.00 d-_bb(j,b) t2_1_bbbb(b,a,i,j) l1_1_bb(i,a)
    //           += -1.00 d-_aa(j,j) t1_1_bb(a,i) l1_1_bb(i,a)
    //           += -1.00 d-_bb(j,j) t1_1_bb(a,i) l1_1_bb(i,a)
    //           += -0.25 <l,k||c,d>_abab t2_bbbb(d,a,i,j) t2_abab(c,b,l,k) l2_1_bbbb(i,j,b,a)
    //           += -0.25 <k,l||c,d>_abab t2_bbbb(d,a,i,j) t2_abab(c,b,k,l) l2_1_bbbb(i,j,b,a)
    //           += -0.50 <j,a||b,c>_aaaa t2_aaaa(b,c,i,j) l1_1_aa(i,a)
    //           += -0.50 <k,j||i,b>_abab t2_abab(a,b,k,j) l1_1_aa(i,a)
    //           += -0.50 <j,k||i,b>_abab t2_abab(a,b,j,k) l1_1_aa(i,a)
    //           += -0.50 <k,j||b,i>_aaaa t2_aaaa(b,a,k,j) l1_1_aa(i,a)
    //           += -1.00 f_aa(j,b) t2_aaaa(b,a,i,j) l1_1_aa(i,a)
    //           += +0.50 <a,j||b,c>_abab t2_abab(b,c,i,j) l1_1_aa(i,a)
    //           += +0.50 <a,j||c,b>_abab t2_abab(c,b,i,j) l1_1_aa(i,a)
    //           += -1.00 d-_aa(a,i) t0_1 l1_1_aa(i,a)
    //           += +1.00 f_bb(j,b) t2_abab(a,b,i,j) l1_1_aa(i,a)
    //           += +1.00 <k,b||c,j>_bbbb t2_bbbb(c,a,i,k) l2_1_bbbb(i,j,b,a)
    //           += -1.00 <k,b||c,j>_abab t2_abab(c,a,k,i) l2_1_bbbb(i,j,b,a)
    //           += -0.50 d-_bb(b,c) t2_abab(a,c,i,j) l2_abab(i,j,a,b)
    //           += -0.50 d-_bb(b,c) t2_abab(a,c,j,i) l2_abab(j,i,a,b)
    //           += -1.00 d-_bb(b,j) t1_1_aa(a,i) l2_1_abab(i,j,a,b)
    sigmal0_1("L")  = -1.00 * tmps_["105_L"]("L");
    tmps_["105_L"].~TArrayD();

    // flops: o1v1L1  = o2v1L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["107_bb_Lvo"]("R,a,i")  = r1["bb"]("R,a,j") * reused_["38_bb_oo"]("j,i");

    // sigmar1_bb += +1.00 d-_bb(j,b) r1_bb(a,j) t1_1_bb(b,i)
    sigmar1_bb("R,a,i") += tmps_["107_bb_Lvo"]("R,a,i");

    // flops: o1v1L1  = o1v1L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["108_bb_Lvo"]("R,a,i")  = r0_1("R") * t1_1["bb"]("a,i");

    // csigmar1_1_bb += +1.00 r0_1 t1_1_bb(a,i)
    csigmar1_1_bb("R,a,i") += tmps_["108_bb_Lvo"]("R,a,i");

    // flops: o1v1L1  = o1v1L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["106_bb_Lvo"]("R,a,i")  = 2.00 * scalars_["6"] * r1["bb"]("R,a,i");

    // flops: o1v1L1  = o2v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o1v1L1 o1v1L1 o1v2L1 o1v1L1 o2v2L1 o1v1L1 o0v0L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o0v0L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o1v1L1 o1v1L1 o3v2L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o1v2L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o1v2L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o2v3L1 o1v1L1 o2v1L1 o1v1L1 o1v2L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o2v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o0v0L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o0v0L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    tmps_["109_bb_Lvo"]("R,a,i")  = -0.50 * r1_1["bb"]("R,a,n") * reused_["31_bb_oo"]("i,n");
    tmps_["109_bb_Lvo"]("R,a,i") -= 0.50 * eri["bbbb_oovo"]("l,n,c,i") * r2_1["bbbb"]("R,c,a,n,l");
    tmps_["109_bb_Lvo"]("R,a,i") += r1_1["bb"]("R,b,n") * reused_["28_bbbb_voov"]("a,i,n,b");
    tmps_["109_bb_Lvo"]("R,a,i") += f["bb_oo"]("l,i") * r1_1["bb"]("R,a,l");
    tmps_["109_bb_Lvo"]("R,a,i") += 2.00 * scalars_["4"] * r1_1["bb"]("R,a,i");
    tmps_["109_bb_Lvo"]("R,a,i") += reused_["108_bb_vv"]("a,b") * r1["bb"]("R,b,i");
    tmps_["109_bb_Lvo"]("R,a,i") -= r1_1["aa"]("R,e,o") * reused_["27_bbaa_voov"]("a,i,o,e");
    tmps_["109_bb_Lvo"]("R,a,i") -= t0_1 * r0_1("R") * reused_["33_bb_vo"]("a,i");
    tmps_["109_bb_Lvo"]("R,a,i") -= w0 * r1_1["bb"]("R,a,i");
    tmps_["109_bb_Lvo"]("R,a,i") += reused_["3_bb_vv"]("b,a") * r1_1["bb"]("R,b,i");
    tmps_["109_bb_Lvo"]("R,a,i") += r1_1["aa"]("R,e,o") * reused_["1_bbaa_voov"]("a,i,o,e");
    tmps_["109_bb_Lvo"]("R,a,i") -= r1["bb"]("R,b,n") * reused_["102_bbbb_voov"]("a,i,n,b");
    tmps_["109_bb_Lvo"]("R,a,i") += t0_1 * r0_1("R") * reused_["32_bb_vo"]("a,i");
    tmps_["109_bb_Lvo"]("R,a,i") += 2.00 * r0_1("R") * reused_["88_bb_vo"]("a,i");
    tmps_["109_bb_Lvo"]("R,a,i") -= 2.00 * r1_1["bb"]("R,a,l") * reused_["38_bb_oo"]("l,i");
    tmps_["109_bb_Lvo"]("R,a,i") += eri["bbbb_vovo"]("a,l,c,i") * r1_1["bb"]("R,c,l");
    tmps_["109_bb_Lvo"]("R,a,i") -= 0.50 * r1["bb"]("R,a,n") * reused_["47_bb_oo"]("i,n");
    tmps_["109_bb_Lvo"]("R,a,i") += scalars_["2"] * r1["bb"]("R,a,i");
    tmps_["109_bb_Lvo"]("R,a,i") += 0.50 * r2["bbbb"]("R,b,a,n,l") * reused_["44_bbbb_ooov"]("i,l,n,b");
    tmps_["109_bb_Lvo"]("R,a,i") += scalars_["5"] * r1_1["bb"]("R,a,i");
    tmps_["109_bb_Lvo"]("R,a,i") -= r1_1["bb"]("R,a,l") * reused_["35_bb_oo"]("l,i");
    tmps_["109_bb_Lvo"]("R,a,i") += r1["bb"]("R,a,l") * reused_["40_bb_oo"]("l,i");
    tmps_["109_bb_Lvo"]("R,a,i") += r0_1("R") * reused_["84_bb_vo"]("a,i");
    tmps_["109_bb_Lvo"]("R,a,i") += eri["baab_vovo"]("a,j,d,i") * r1_1["aa"]("R,d,j");
    tmps_["109_bb_Lvo"]("R,a,i") -= f["bb_vo"]("a,i") * r0_1("R");
    tmps_["109_bb_Lvo"]("R,a,i") -= 0.50 * reused_["29_bb_vv"]("a,b") * r1_1["bb"]("R,b,i");
    tmps_["109_bb_Lvo"]("R,a,i") += r1["bb"]("R,a,n") * reused_["48_bb_oo"]("n,i");
    tmps_["109_bb_Lvo"]("R,a,i") += r1["aa"]("R,e,j") * reused_["105_baab_vovo"]("a,j,e,i");
    tmps_["109_bb_Lvo"]("R,a,i") -= r1_1["bb"]("R,b,n") * reused_["26_bbbb_voov"]("a,i,n,b");
    tmps_["109_bb_Lvo"]("R,a,i") += r2["bbbb"]("R,b,a,i,n") * reused_["11_bb_ov"]("n,b");
    tmps_["109_bb_Lvo"]("R,a,i") += 2.00 * scalars_["3"] * r1_1["bb"]("R,a,i");
    tmps_["109_bb_Lvo"]("R,a,i") += r1["aa"]("R,e,o") * reused_["103_bbaa_voov"]("a,i,o,e");
    tmps_["109_bb_Lvo"]("R,a,i") += scalars_["1"] * r1["bb"]("R,a,i");
    tmps_["109_bb_Lvo"]("R,a,i") -= 0.50 * reused_["107_bb_vv"]("a,b") * r1["bb"]("R,b,i");
    tmps_["109_bb_Lvo"]("R,a,i") -= r1["aa"]("R,e,o") * reused_["101_bbaa_voov"]("a,i,o,e");
    tmps_["109_bb_Lvo"]("R,a,i") += eri["baab_vovv"]("a,j,d,b") * r2_1["abab"]("R,d,b,j,i");
    tmps_["109_bb_Lvo"]("R,a,i") += eri["abab_oovo"]("o,l,d,i") * r2_1["abab"]("R,d,a,o,l");
    tmps_["109_bb_Lvo"]("R,a,i") -= r2["abab"]("R,e,a,o,i") * reused_["10_aa_ov"]("o,e");
    tmps_["109_bb_Lvo"]("R,a,i") -= r2_1["bbbb"]("R,c,a,i,l") * reused_["12_bb_ov"]("l,c");
    tmps_["109_bb_Lvo"]("R,a,i") -= f["aa_ov"]("j,d") * r2_1["abab"]("R,d,a,j,i");
    tmps_["109_bb_Lvo"]("R,a,i") -= r2["abab"]("R,e,a,o,i") * reused_["9_aa_ov"]("o,e");
    tmps_["109_bb_Lvo"]("R,a,i") += r2["bbbb"]("R,b,a,i,n") * reused_["8_bb_ov"]("n,b");
    tmps_["109_bb_Lvo"]("R,a,i") += f["bb_ov"]("l,c") * r2_1["bbbb"]("R,c,a,i,l");
    tmps_["109_bb_Lvo"]("R,a,i") += r1["bb"]("R,b,n") * reused_["104_bbbb_voov"]("a,i,n,b");
    tmps_["109_bb_Lvo"]("R,a,i") += reused_["34_bb_vv"]("a,c") * r1_1["bb"]("R,c,i");
    tmps_["109_bb_Lvo"]("R,a,i") += r2_1["abab"]("R,d,a,j,i") * reused_["13_aa_ov"]("j,d");
    tmps_["109_bb_Lvo"]("R,a,i") += r2["abab"]("R,e,a,o,l") * reused_["43_abab_oovo"]("o,l,e,i");
    tmps_["109_bb_Lvo"]("R,a,i") -= 0.50 * eri["bbbb_vovv"]("a,l,c,b") * r2_1["bbbb"]("R,c,b,i,l");
    tmps_["109_bb_Lvo"]("R,a,i") += r1_1["bb"]("R,a,n") * reused_["30_bb_oo"]("n,i");
    tmps_["109_bb_Lvo"]("R,a,i") -= f["bb_vv"]("a,c") * r1_1["bb"]("R,c,i");
    tmps_["109_bb_Lvo"]("R,a,i") -= r1["bb"]("R,b,l") * reused_["106_bbbb_vovo"]("a,l,b,i");
    tmps_["109_bb_Lvo"]("R,a,i") -= 0.50 * tmps_["106_bb_Lvo"]("R,a,i");
    tmps_["109_bb_Lvo"]("R,a,i") -= t0_1 * tmps_["107_bb_Lvo"]("R,a,i");
    tmps_["109_bb_Lvo"]("R,a,i") -= t1_1["bb"]("a,l") * tmps_["53_bb_Loo"]("R,i,l");
    tmps_["109_bb_Lvo"]("R,a,i") -= 2.00 * t1_1["bb"]("a,l") * tmps_["49_bb_Loo"]("R,i,l");
    tmps_["109_bb_Lvo"]("R,a,i") += t1_1["bb"]("a,i") * tmps_["43_L"]("R");
    tmps_["109_bb_Lvo"]("R,a,i") += t1_1["bb"]("a,l") * tmps_["54_bb_Loo"]("R,l,i");
    tmps_["109_bb_Lvo"]("R,a,i") += t1_1["bb"]("a,l") * tmps_["56_bb_Loo"]("R,i,l");
    tmps_["109_bb_Lvo"]("R,a,i") += scalars_["2"] * tmps_["108_bb_Lvo"]("R,a,i");
    tmps_["109_bb_Lvo"]("R,a,i") += scalars_["1"] * tmps_["108_bb_Lvo"]("R,a,i");
    tmps_["106_bb_Lvo"].~TArrayD();
}
