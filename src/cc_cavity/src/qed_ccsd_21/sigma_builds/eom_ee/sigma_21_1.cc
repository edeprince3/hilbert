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

void hilbert::EOM_EE_QED_CCSD::sigma_ee_21_1() {

    // Get cavity information
    double w0 = cc_wfn_->cavity_frequency_[2];
    double coupling_factor_z = w0 * cc_wfn_->cavity_coupling_strength_[2];

    // get effective dipole integrals
    TArrayMap dp = reinterpret_pointer_cast<QED_CCSD>(cc_wfn_)->effective_dipole();

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


    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["1_bbbb_Lvovo"]("R,a,i,b,j")  = r0_1("R") * reused_["2_bbbb_vovo"]("a,i,b,j");

    // sigmar2_1_bbbb  = +1.00 P(i,j) <k,l||c,d>_abab r0_1 t2_abab(c,a,k,j) t2_bbbb(d,b,i,l)
    sigmar2_1_bbbb("R,a,b,i,j")  = tmps_["1_bbbb_Lvovo"]("R,a,j,b,i");
    sigmar2_1_bbbb("R,a,b,i,j") -= tmps_["1_bbbb_Lvovo"]("R,a,i,b,j");

    // sigmar2_1_bbbb += +1.00 P(i,j) <l,k||d,c>_abab r0_1 t2_bbbb(c,a,j,k) t2_abab(d,b,l,i)
    sigmar2_1_bbbb("R,a,b,i,j") += tmps_["1_bbbb_Lvovo"]("R,b,i,a,j");
    sigmar2_1_bbbb("R,a,b,i,j") -= tmps_["1_bbbb_Lvovo"]("R,b,j,a,i");
    tmps_["1_bbbb_Lvovo"].~TArrayD();

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["2_bbbb_Lvoov"]("R,a,i,j,b")  = r0_1("R") * reused_["4_bbbb_voov"]("a,i,j,b");

    // sigmar2_1_bbbb += -0.50 <l,k||c,d>_abab r0_1 t2_abab(c,a,l,k) t2_bbbb(d,b,i,j)
    //                += -0.50 <k,l||c,d>_abab r0_1 t2_abab(c,a,k,l) t2_bbbb(d,b,i,j)
    sigmar2_1_bbbb("R,a,b,i,j") -= tmps_["2_bbbb_Lvoov"]("R,b,i,j,a");

    // sigmar2_1_bbbb += +0.50 <l,k||d,c>_abab r0_1 t2_bbbb(c,a,i,j) t2_abab(d,b,l,k)
    //                += +0.50 <k,l||d,c>_abab r0_1 t2_bbbb(c,a,i,j) t2_abab(d,b,k,l)
    sigmar2_1_bbbb("R,a,b,i,j") += tmps_["2_bbbb_Lvoov"]("R,a,i,j,b");
    tmps_["2_bbbb_Lvoov"].~TArrayD();

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["3_aaaa_Lvvoo"]("R,a,b,i,j")  = r0_1("R") * reused_["5_aaaa_vvoo"]("a,b,i,j");

    // sigmar2_1_aaaa  = +0.50 <l,k||c,d>_abab r0_1 t2_aaaa(c,a,i,j) t2_abab(b,d,l,k)
    //                += +0.50 <k,l||c,d>_abab r0_1 t2_aaaa(c,a,i,j) t2_abab(b,d,k,l)
    sigmar2_1_aaaa("R,a,b,i,j")  = tmps_["3_aaaa_Lvvoo"]("R,b,a,i,j");

    // sigmar2_1_aaaa += -0.50 <l,k||d,c>_abab r0_1 t2_abab(a,c,l,k) t2_aaaa(d,b,i,j)
    //                += -0.50 <k,l||d,c>_abab r0_1 t2_abab(a,c,k,l) t2_aaaa(d,b,i,j)
    sigmar2_1_aaaa("R,a,b,i,j") -= tmps_["3_aaaa_Lvvoo"]("R,a,b,i,j");
    tmps_["3_aaaa_Lvvoo"].~TArrayD();

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["4_aaaa_Lvovo"]("R,a,i,b,j")  = r0_1("R") * reused_["7_aaaa_vovo"]("a,i,b,j");

    // sigmar2_1_aaaa += +1.00 P(i,j) <l,k||d,c>_abab r0_1 t2_abab(a,c,j,k) t2_aaaa(d,b,i,l)
    sigmar2_1_aaaa("R,a,b,i,j") += tmps_["4_aaaa_Lvovo"]("R,a,j,b,i");
    sigmar2_1_aaaa("R,a,b,i,j") -= tmps_["4_aaaa_Lvovo"]("R,a,i,b,j");

    // sigmar2_1_aaaa += +1.00 P(i,j) <k,l||c,d>_abab r0_1 t2_aaaa(c,a,j,k) t2_abab(b,d,i,l)
    sigmar2_1_aaaa("R,a,b,i,j") += tmps_["4_aaaa_Lvovo"]("R,b,i,a,j");
    sigmar2_1_aaaa("R,a,b,i,j") -= tmps_["4_aaaa_Lvovo"]("R,b,j,a,i");
    tmps_["4_aaaa_Lvovo"].~TArrayD();

    // flops: o0v0L1  = o0v0L1
    //  mems: o0v0L1  = o0v0L1
    tmps_["5_L"]("R")  = t0_1 * r0_1("R");

    // csigmar0_1  = +1.00 r0_1 t0_1
    csigmar0_1("R")  = tmps_["5_L"]("R");
    tmps_["5_L"].~TArrayD();

    // flops: o0v0L1  = o0v0L1
    //  mems: o0v0L1  = o0v0L1
    tmps_["6_L"]("L")  = t0_1 * l0_1("L");

    // csigmal0_1  = +1.00 t0_1 l0_1
    csigmal0_1("L")  = tmps_["6_L"]("L");
    tmps_["6_L"].~TArrayD();

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["7_bbbb_Lvovo"]("R,a,i,b,j")  = r1["bb"]("R,a,i") * t1_1["bb"]("b,j");

    // csigmar2_bbbb  = -1.00 P(i,j) P(a,b) r1_bb(b,i) t1_1_bb(a,j)
    csigmar2_bbbb("R,a,b,i,j")  = -1.00 * tmps_["7_bbbb_Lvovo"]("R,b,i,a,j");
    csigmar2_bbbb("R,a,b,i,j") += tmps_["7_bbbb_Lvovo"]("R,b,j,a,i");
    csigmar2_bbbb("R,a,b,i,j") += tmps_["7_bbbb_Lvovo"]("R,a,i,b,j");
    csigmar2_bbbb("R,a,b,i,j") -= tmps_["7_bbbb_Lvovo"]("R,a,j,b,i");
    tmps_["7_bbbb_Lvovo"].~TArrayD();

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["8_bbbb_Lvovo"]("R,a,i,b,j")  = r1_1["bb"]("R,a,i") * t1_1["bb"]("b,j");

    // csigmar2_1_bbbb  = -1.00 P(i,j) P(a,b) r1_1_bb(b,i) t1_1_bb(a,j)
    csigmar2_1_bbbb("R,a,b,i,j")  = -1.00 * tmps_["8_bbbb_Lvovo"]("R,b,i,a,j");
    csigmar2_1_bbbb("R,a,b,i,j") += tmps_["8_bbbb_Lvovo"]("R,b,j,a,i");
    csigmar2_1_bbbb("R,a,b,i,j") += tmps_["8_bbbb_Lvovo"]("R,a,i,b,j");
    csigmar2_1_bbbb("R,a,b,i,j") -= tmps_["8_bbbb_Lvovo"]("R,a,j,b,i");
    tmps_["8_bbbb_Lvovo"].~TArrayD();

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["9_aaaa_Lvovo"]("R,a,i,b,j")  = r1["aa"]("R,a,i") * t1_1["aa"]("b,j");

    // csigmar2_aaaa  = -1.00 P(i,j) P(a,b) r1_aa(b,i) t1_1_aa(a,j)
    csigmar2_aaaa("R,a,b,i,j")  = -1.00 * tmps_["9_aaaa_Lvovo"]("R,b,i,a,j");
    csigmar2_aaaa("R,a,b,i,j") += tmps_["9_aaaa_Lvovo"]("R,b,j,a,i");
    csigmar2_aaaa("R,a,b,i,j") += tmps_["9_aaaa_Lvovo"]("R,a,i,b,j");
    csigmar2_aaaa("R,a,b,i,j") -= tmps_["9_aaaa_Lvovo"]("R,a,j,b,i");
    tmps_["9_aaaa_Lvovo"].~TArrayD();

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["10_aaaa_Lvovo"]("R,a,i,b,j")  = r1_1["aa"]("R,a,i") * t1_1["aa"]("b,j");

    // csigmar2_1_aaaa  = -1.00 P(i,j) P(a,b) r1_1_aa(b,i) t1_1_aa(a,j)
    csigmar2_1_aaaa("R,a,b,i,j")  = -1.00 * tmps_["10_aaaa_Lvovo"]("R,b,i,a,j");
    csigmar2_1_aaaa("R,a,b,i,j") += tmps_["10_aaaa_Lvovo"]("R,b,j,a,i");
    csigmar2_1_aaaa("R,a,b,i,j") += tmps_["10_aaaa_Lvovo"]("R,a,i,b,j");
    csigmar2_1_aaaa("R,a,b,i,j") -= tmps_["10_aaaa_Lvovo"]("R,a,j,b,i");
    tmps_["10_aaaa_Lvovo"].~TArrayD();

    // flops: o0v0L1  = o2v2L1 o1v1L1 o0v0L1 o0v0L1 o1v1L1 o1v1L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o1v1L1 o0v0L1 o0v0L1 o0v0L1 o1v1L1 o0v0L1 o0v0L1 o1v1L1 o0v0L1 o2v2L1 o0v0L1 o2v2L1 o0v0L1 o1v1L1 o0v0L1 o1v1L1 o0v0L1 o0v0L1 o0v0L1
    //  mems: o0v0L1  = o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1
    tmps_["11_L"]("R")  = -0.25 * eri["aaaa_oovv"]("j,i,b,a") * r2_1["aaaa"]("R,b,a,i,j");
    tmps_["11_L"]("R") -= r1_1["aa"]("R,b,j") * reused_["13_aa_ov"]("j,b");
    tmps_["11_L"]("R") -= 2.00 * scalars_["3"] * r0_1("R");
    tmps_["11_L"]("R") -= 2.00 * scalars_["4"] * r0_1("R");
    tmps_["11_L"]("R") += f["bb_ov"]("k,c") * r1_1["bb"]("R,c,k");
    tmps_["11_L"]("R") += r1["bb"]("R,d,l") * reused_["8_bb_ov"]("l,d");
    tmps_["11_L"]("R") += r1["bb"]("R,d,l") * reused_["11_bb_ov"]("l,d");
    tmps_["11_L"]("R") -= scalars_["5"] * r0_1("R");
    tmps_["11_L"]("R") -= r1_1["bb"]("R,c,k") * reused_["12_bb_ov"]("k,c");
    tmps_["11_L"]("R") += r1["aa"]("R,a,i") * reused_["10_aa_ov"]("i,a");
    tmps_["11_L"]("R") -= 0.25 * eri["bbbb_oovv"]("k,l,c,d") * r2_1["bbbb"]("R,c,d,l,k");
    tmps_["11_L"]("R") += eri["abab_oovv"]("i,k,b,d") * r2_1["abab"]("R,b,d,i,k");
    tmps_["11_L"]("R") += f["aa_ov"]("j,b") * r1_1["aa"]("R,b,j");
    tmps_["11_L"]("R") += r1["aa"]("R,a,i") * reused_["9_aa_ov"]("i,a");
    tmps_["11_L"]("R") += w0 * r0_1("R");

    // sigmar0_1  = +1.00 f_bb(i,a) r1_1_bb(a,i)
    //           += +1.00 <i,j||a,b>_abab r1_bb(b,j) t1_1_aa(a,i)
    //           += -2.00 d-_aa(i,a) r0_1 t1_1_aa(a,i)
    //           += -2.00 d-_bb(i,a) r0_1 t1_1_bb(a,i)
    //           += -1.00 d-_aa(i,a) r1_1_aa(a,i) t0_1
    //           += -1.00 <j,i||a,b>_bbbb r1_bb(b,j) t1_1_bb(a,i)
    //           += -1.00 d-_aa(i,i) r0_1 t0_1
    //           += +1.00 f_aa(i,i) l0_1
    //           += -0.50 <j,i||j,i>_abab l0_1
    //           += -0.50 <i,j||i,j>_abab l0_1
    //           += +0.25 <j,i||a,b>_aaaa t2_aaaa(a,b,j,i) l0_1
    //           += -0.50 <j,i||j,i>_bbbb l0_1
    //           += -0.50 <j,i||j,i>_aaaa l0_1
    //           += +0.25 <j,i||a,b>_bbbb t2_bbbb(a,b,j,i) l0_1
    //           += +1.00 f_bb(i,i) l0_1
    //           += -1.00 d-_bb(i,i) t0_1 l0_1
    //           += +0.25 <j,i||a,b>_abab t2_abab(a,b,j,i) l0_1
    //           += +0.25 <i,j||a,b>_abab t2_abab(a,b,i,j) l0_1
    //           += +0.25 <j,i||b,a>_abab t2_abab(b,a,j,i) l0_1
    //           += +0.25 <i,j||b,a>_abab t2_abab(b,a,i,j) l0_1
    //           += -1.00 d-_bb(i,a) r1_1_bb(a,i) t0_1
    //           += +0.25 <j,i||a,b>_aaaa r2_1_aaaa(a,b,j,i)
    //           += -1.00 <j,i||a,b>_aaaa r1_aa(b,j) t1_1_aa(a,i)
    //           += +0.25 <j,i||a,b>_bbbb r2_1_bbbb(a,b,j,i)
    //           += +0.25 <j,i||a,b>_abab r2_1_abab(a,b,j,i)
    //           += +0.25 <i,j||a,b>_abab r2_1_abab(a,b,i,j)
    //           += +0.25 <j,i||b,a>_abab r2_1_abab(b,a,j,i)
    //           += +0.25 <i,j||b,a>_abab r2_1_abab(b,a,i,j)
    //           += +1.00 f_aa(i,a) r1_1_aa(a,i)
    //           += +1.00 <j,i||b,a>_abab r1_aa(b,j) t1_1_bb(a,i)
    //           += +1.00 r0_1 w0
    sigmar0_1("R")  = tmps_["11_L"]("R");
    tmps_["11_L"].~TArrayD();

    // flops: o1v1L1  = o1v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1
    tmps_["12_bb_Lvo"]("R,a,i")  = r1["bb"]("R,a,i");
    tmps_["12_bb_Lvo"]("R,a,i") += t0_1 * r1_1["bb"]("R,a,i");

    // csigmar1_1_bb  = +1.00 r1_1_bb(a,i) t0_1
    //               += +1.00 r1_bb(a,i)
    csigmar1_1_bb("R,a,i")  = tmps_["12_bb_Lvo"]("R,a,i");
    tmps_["12_bb_Lvo"].~TArrayD();

    // flops: o1v1L1  = o1v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1
    tmps_["13_bb_Lvo"]("R,a,i")  = r1_1["bb"]("R,a,i");
    tmps_["13_bb_Lvo"]("R,a,i") += t0_1 * r1["bb"]("R,a,i");

    // csigmar1_bb  = +1.00 r1_bb(a,i) t0_1
    //             += +1.00 r1_1_bb(a,i)
    csigmar1_bb("R,a,i")  = tmps_["13_bb_Lvo"]("R,a,i");
    tmps_["13_bb_Lvo"].~TArrayD();

    // flops: o1v1L1  = o1v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1
    tmps_["14_bb_Lov"]("L,i,a")  = l1_1["bb"]("L,i,a");
    tmps_["14_bb_Lov"]("L,i,a") += t0_1 * l1["bb"]("L,i,a");

    // csigmal1_bb  = +1.00 l1_1_bb(i,a)
    //             += +1.00 t0_1 l1_bb(i,a)
    csigmal1_bb("L,a,i")  = tmps_["14_bb_Lov"]("L,i,a");
    tmps_["14_bb_Lov"].~TArrayD();

    // flops: o1v1L1  = o1v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1
    tmps_["15_bb_Lov"]("L,i,a")  = l1["bb"]("L,i,a");
    tmps_["15_bb_Lov"]("L,i,a") += t0_1 * l1_1["bb"]("L,i,a");

    // csigmal1_1_bb  = +1.00 l1_bb(i,a)
    //               += +1.00 t0_1 l1_1_bb(i,a)
    csigmal1_1_bb("L,a,i")  = tmps_["15_bb_Lov"]("L,i,a");
    tmps_["15_bb_Lov"].~TArrayD();

    // flops: o1v1L1  = o1v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1
    tmps_["16_aa_Lvo"]("R,a,i")  = r1["aa"]("R,a,i");
    tmps_["16_aa_Lvo"]("R,a,i") += t0_1 * r1_1["aa"]("R,a,i");

    // csigmar1_1_aa  = +1.00 r1_1_aa(a,i) t0_1
    //               += +1.00 r1_aa(a,i)
    csigmar1_1_aa("R,a,i")  = tmps_["16_aa_Lvo"]("R,a,i");
    tmps_["16_aa_Lvo"].~TArrayD();

    // flops: o1v1L1  = o1v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1
    tmps_["17_aa_Lvo"]("R,a,i")  = r1_1["aa"]("R,a,i");
    tmps_["17_aa_Lvo"]("R,a,i") += t0_1 * r1["aa"]("R,a,i");

    // csigmar1_aa  = +1.00 r1_aa(a,i) t0_1
    //             += +1.00 r1_1_aa(a,i)
    csigmar1_aa("R,a,i")  = tmps_["17_aa_Lvo"]("R,a,i");
    tmps_["17_aa_Lvo"].~TArrayD();

    // flops: o1v1L1  = o1v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1
    tmps_["18_aa_Lov"]("L,i,a")  = l1_1["aa"]("L,i,a");
    tmps_["18_aa_Lov"]("L,i,a") += t0_1 * l1["aa"]("L,i,a");

    // csigmal1_aa  = +1.00 l1_1_aa(i,a)
    //             += +1.00 t0_1 l1_aa(i,a)
    csigmal1_aa("L,a,i")  = tmps_["18_aa_Lov"]("L,i,a");
    tmps_["18_aa_Lov"].~TArrayD();

    // flops: o1v1L1  = o1v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1
    tmps_["19_aa_Lov"]("L,i,a")  = l1["aa"]("L,i,a");
    tmps_["19_aa_Lov"]("L,i,a") += t0_1 * l1_1["aa"]("L,i,a");

    // csigmal1_1_aa  = +1.00 t0_1 l1_1_aa(i,a)
    //               += +1.00 l1_aa(i,a)
    csigmal1_1_aa("L,a,i")  = tmps_["19_aa_Lov"]("L,i,a");
    tmps_["19_aa_Lov"].~TArrayD();

    // flops: o2v2L1  = o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1
    tmps_["20_bbbb_Lvvoo"]("R,a,b,i,j")  = r2_1["bbbb"]("R,a,b,i,j");
    tmps_["20_bbbb_Lvvoo"]("R,a,b,i,j") += t0_1 * r2["bbbb"]("R,a,b,i,j");

    // csigmar2_bbbb += +1.00 r2_1_bbbb(a,b,i,j)
    //               += +1.00 r2_bbbb(a,b,i,j) t0_1
    csigmar2_bbbb("R,a,b,i,j") += tmps_["20_bbbb_Lvvoo"]("R,a,b,i,j");
    tmps_["20_bbbb_Lvvoo"].~TArrayD();

    // flops: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["21_bbbb_Lvvoo"]("R,a,b,i,j")  = r0_1("R") * t2_1["bbbb"]("a,b,i,j");
    tmps_["21_bbbb_Lvvoo"]("R,a,b,i,j") += r2["bbbb"]("R,a,b,i,j");
    tmps_["21_bbbb_Lvvoo"]("R,a,b,i,j") += t0_1 * r2_1["bbbb"]("R,a,b,i,j");

    // csigmar2_1_bbbb += +1.00 r0_1 t2_1_bbbb(a,b,i,j)
    //                 += +1.00 r2_bbbb(a,b,i,j)
    //                 += +1.00 r2_1_bbbb(a,b,i,j) t0_1
    csigmar2_1_bbbb("R,a,b,i,j") += tmps_["21_bbbb_Lvvoo"]("R,a,b,i,j");
    tmps_["21_bbbb_Lvvoo"].~TArrayD();

    // flops: o2v2L1  = o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1
    tmps_["22_bbbb_Loovv"]("L,i,j,a,b")  = l2_1["bbbb"]("L,i,j,a,b");
    tmps_["22_bbbb_Loovv"]("L,i,j,a,b") += t0_1 * l2["bbbb"]("L,i,j,a,b");

    // csigmal2_bbbb  = +1.00 t0_1 l2_bbbb(i,j,a,b)
    //               += +1.00 l2_1_bbbb(i,j,a,b)
    csigmal2_bbbb("L,a,b,i,j")  = tmps_["22_bbbb_Loovv"]("L,i,j,a,b");
    tmps_["22_bbbb_Loovv"].~TArrayD();

    // flops: o2v2L1  = o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1
    tmps_["23_bbbb_Loovv"]("L,i,j,a,b")  = l2["bbbb"]("L,i,j,a,b");
    tmps_["23_bbbb_Loovv"]("L,i,j,a,b") += t0_1 * l2_1["bbbb"]("L,i,j,a,b");

    // csigmal2_1_bbbb  = +1.00 l2_bbbb(i,j,a,b)
    //                 += +1.00 t0_1 l2_1_bbbb(i,j,a,b)
    csigmal2_1_bbbb("L,a,b,i,j")  = tmps_["23_bbbb_Loovv"]("L,i,j,a,b");
    tmps_["23_bbbb_Loovv"].~TArrayD();

    // flops: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["24_abab_Lvvoo"]("R,a,b,i,j")  = r1["bb"]("R,b,j") * t1_1["aa"]("a,i");
    tmps_["24_abab_Lvvoo"]("R,a,b,i,j") += r2_1["abab"]("R,a,b,i,j");
    tmps_["24_abab_Lvvoo"]("R,a,b,i,j") += t0_1 * r2["abab"]("R,a,b,i,j");
    tmps_["24_abab_Lvvoo"]("R,a,b,i,j") += r1["aa"]("R,a,i") * t1_1["bb"]("b,j");

    // csigmar2_abab  = +1.00 r2_abab(a,b,i,j) t0_1
    //               += +1.00 r2_1_abab(a,b,i,j)
    //               += +1.00 r1_aa(a,i) t1_1_bb(b,j)
    //               += +1.00 r1_bb(b,j) t1_1_aa(a,i)
    csigmar2_abab("R,a,b,i,j")  = tmps_["24_abab_Lvvoo"]("R,a,b,i,j");
    tmps_["24_abab_Lvvoo"].~TArrayD();

    // flops: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["25_abab_Lvvoo"]("R,a,b,i,j")  = r1_1["bb"]("R,b,j") * t1_1["aa"]("a,i");
    tmps_["25_abab_Lvvoo"]("R,a,b,i,j") += r2["abab"]("R,a,b,i,j");
    tmps_["25_abab_Lvvoo"]("R,a,b,i,j") += t0_1 * r2_1["abab"]("R,a,b,i,j");
    tmps_["25_abab_Lvvoo"]("R,a,b,i,j") += r1_1["aa"]("R,a,i") * t1_1["bb"]("b,j");

    // csigmar2_1_abab  = +1.00 r2_1_abab(a,b,i,j) t0_1
    //                 += +1.00 r2_abab(a,b,i,j)
    //                 += +1.00 r1_1_aa(a,i) t1_1_bb(b,j)
    //                 += +1.00 r1_1_bb(b,j) t1_1_aa(a,i)
    csigmar2_1_abab("R,a,b,i,j")  = tmps_["25_abab_Lvvoo"]("R,a,b,i,j");
    tmps_["25_abab_Lvvoo"].~TArrayD();

    // flops: o2v2L1  = o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1
    tmps_["26_abab_Loovv"]("L,i,j,a,b")  = l2_1["abab"]("L,i,j,a,b");
    tmps_["26_abab_Loovv"]("L,i,j,a,b") += t0_1 * l2["abab"]("L,i,j,a,b");

    // csigmal2_abab  = +1.00 l2_1_abab(i,j,a,b)
    //               += +1.00 t0_1 l2_abab(i,j,a,b)
    csigmal2_abab("L,a,b,i,j")  = tmps_["26_abab_Loovv"]("L,i,j,a,b");
    tmps_["26_abab_Loovv"].~TArrayD();

    // flops: o2v2L1  = o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1
    tmps_["27_abab_Loovv"]("L,i,j,a,b")  = l2["abab"]("L,i,j,a,b");
    tmps_["27_abab_Loovv"]("L,i,j,a,b") += t0_1 * l2_1["abab"]("L,i,j,a,b");

    // csigmal2_1_abab  = +1.00 l2_abab(i,j,a,b)
    //                 += +1.00 t0_1 l2_1_abab(i,j,a,b)
    csigmal2_1_abab("L,a,b,i,j")  = tmps_["27_abab_Loovv"]("L,i,j,a,b");
    tmps_["27_abab_Loovv"].~TArrayD();

    // flops: o2v2L1  = o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1
    tmps_["28_aaaa_Lvvoo"]("R,a,b,i,j")  = r2_1["aaaa"]("R,a,b,i,j");
    tmps_["28_aaaa_Lvvoo"]("R,a,b,i,j") += t0_1 * r2["aaaa"]("R,a,b,i,j");

    // csigmar2_aaaa += +1.00 r2_aaaa(a,b,i,j) t0_1
    //               += +1.00 r2_1_aaaa(a,b,i,j)
    csigmar2_aaaa("R,a,b,i,j") += tmps_["28_aaaa_Lvvoo"]("R,a,b,i,j");
    tmps_["28_aaaa_Lvvoo"].~TArrayD();

    // flops: o2v2L1  = o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1
    tmps_["29_aaaa_Lvvoo"]("R,a,b,i,j")  = r2["aaaa"]("R,a,b,i,j");
    tmps_["29_aaaa_Lvvoo"]("R,a,b,i,j") += t0_1 * r2_1["aaaa"]("R,a,b,i,j");

    // csigmar2_1_aaaa += +1.00 r2_1_aaaa(a,b,i,j) t0_1
    //                 += +1.00 r2_aaaa(a,b,i,j)
    csigmar2_1_aaaa("R,a,b,i,j") += tmps_["29_aaaa_Lvvoo"]("R,a,b,i,j");
    tmps_["29_aaaa_Lvvoo"].~TArrayD();

    // flops: o2v2L1  = o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1
    tmps_["30_aaaa_Loovv"]("L,i,j,a,b")  = l2_1["aaaa"]("L,i,j,a,b");
    tmps_["30_aaaa_Loovv"]("L,i,j,a,b") += t0_1 * l2["aaaa"]("L,i,j,a,b");

    // csigmal2_aaaa  = +1.00 t0_1 l2_aaaa(i,j,a,b)
    //               += +1.00 l2_1_aaaa(i,j,a,b)
    csigmal2_aaaa("L,a,b,i,j")  = tmps_["30_aaaa_Loovv"]("L,i,j,a,b");
    tmps_["30_aaaa_Loovv"].~TArrayD();

    // flops: o2v2L1  = o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1
    tmps_["31_aaaa_Loovv"]("L,i,j,a,b")  = l2["aaaa"]("L,i,j,a,b");
    tmps_["31_aaaa_Loovv"]("L,i,j,a,b") += t0_1 * l2_1["aaaa"]("L,i,j,a,b");

    // csigmal2_1_aaaa  = +1.00 t0_1 l2_1_aaaa(i,j,a,b)
    //                 += +1.00 l2_aaaa(i,j,a,b)
    csigmal2_1_aaaa("L,a,b,i,j")  = tmps_["31_aaaa_Loovv"]("L,i,j,a,b");
    tmps_["31_aaaa_Loovv"].~TArrayD();

    // flops: o2v0L1  = o2v1L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["32_aa_Loo"]("R,i,j")  = dp["aa_ov"]("i,a") * r1["aa"]("R,a,j");
    tmps_["32_aa_Loo"].~TArrayD();

    // flops: o2v2L1  = o2v1L1 o3v2L1
    //  mems: o2v2L1  = o2v0L1 o2v2L1
    tmps_["33_aaaa_Lovvo"]("R,i,a,b,j")  = dp["aa_ov"]("k,c") * r1["aa"]("R,c,i") * t2["aaaa"]("a,b,j,k");

    // sigmar2_1_aaaa += -1.00 P(i,j) d+_aa(k,c) r1_aa(c,i) t2_aaaa(a,b,j,k)
    sigmar2_1_aaaa("R,a,b,i,j") -= tmps_["33_aaaa_Lovvo"]("R,i,a,b,j");
    sigmar2_1_aaaa("R,a,b,i,j") += tmps_["33_aaaa_Lovvo"]("R,j,a,b,i");

    // flops: o2v0L1  = o3v1L1 o3v1L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
    tmps_["37_aa_Loo"]("R,i,j")  = eri["aaaa_oovo"]("i,l,b,j") * r1["aa"]("R,b,l");
    tmps_["37_aa_Loo"]("R,i,j") += eri["abba_oovo"]("i,k,a,j") * r1["bb"]("R,a,k");

    // flops: o2v0L1  = o3v2L1 o2v1L1 o3v2L1 o2v0L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1
    tmps_["36_aa_Loo"]("R,i,j")  = r2["abab"]("R,a,c,i,l") * reused_["129_abab_oovv"]("j,l,a,c");
    tmps_["36_aa_Loo"]("R,i,j") += 2.00 * f["aa_ov"]("j,a") * r1["aa"]("R,a,i");
    tmps_["36_aa_Loo"]("R,i,j") += eri["aaaa_oovv"]("j,k,a,b") * r2["aaaa"]("R,a,b,i,k");

    // flops: o2v0L1  = o2v1L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["35_aa_Loo"]("R,i,j")  = r1["aa"]("R,a,i") * dp["aa_ov"]("j,a");

    // flops: o2v0L1  = o2v1L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["34_aa_Loo"]("R,i,j")  = r1_1["aa"]("R,a,i") * dp["aa_ov"]("j,a");

    // flops: o2v2L1  = o3v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["38_aaaa_Lovvo"]("R,i,a,b,j")  = -0.50 * t2["aaaa"]("a,b,j,k") * tmps_["36_aa_Loo"]("R,i,k");
    tmps_["38_aaaa_Lovvo"]("R,i,a,b,j") += t0_1 * tmps_["33_aaaa_Lovvo"]("R,i,a,b,j");
    tmps_["38_aaaa_Lovvo"]("R,i,a,b,j") -= t2["aaaa"]("a,b,i,k") * tmps_["37_aa_Loo"]("R,k,j");
    tmps_["38_aaaa_Lovvo"]("R,i,a,b,j") += t2["aaaa"]("a,b,j,k") * tmps_["34_aa_Loo"]("R,i,k");
    tmps_["38_aaaa_Lovvo"]("R,i,a,b,j") += t2_1["aaaa"]("a,b,j,k") * tmps_["35_aa_Loo"]("R,i,k");
    tmps_["38_aaaa_Lovvo"]("R,i,a,b,j") -= dp["aa_oo"]("k,j") * r2_1["aaaa"]("R,a,b,i,k");
    tmps_["38_aaaa_Lovvo"]("R,i,a,b,j") += 0.50 * r1["aa"]("R,c,i") * reused_["14_aaaa_vovv"]("c,j,a,b");
    tmps_["38_aaaa_Lovvo"]("R,i,a,b,j") -= 0.50 * r2["aaaa"]("R,a,b,i,l") * reused_["15_aa_oo"]("j,l");
    tmps_["38_aaaa_Lovvo"]("R,i,a,b,j") -= r2["aaaa"]("R,a,b,i,k") * reused_["17_aa_oo"]("k,j");
    tmps_["38_aaaa_Lovvo"]("R,i,a,b,j") += f["aa_oo"]("k,j") * r2["aaaa"]("R,a,b,i,k");
    tmps_["38_aaaa_Lovvo"]("R,i,a,b,j") += r2["aaaa"]("R,a,b,i,l") * reused_["16_aa_oo"]("l,j");
    tmps_["38_aaaa_Lovvo"]("R,i,a,b,j") -= eri["aaaa_vvvo"]("a,b,c,j") * r1["aa"]("R,c,i");
    tmps_["33_aaaa_Lovvo"].~TArrayD();

    // sigmar2_aaaa  = -1.00 P(i,j) d-_aa(k,c) r1_aa(c,i) t2_aaaa(a,b,j,k) t0_1
    //              += -0.50 P(i,j) <l,k||c,d>_aaaa r2_aaaa(c,d,i,l) t2_aaaa(a,b,j,k)
    //              += +1.00 P(i,j) f_aa(k,c) r1_aa(c,i) t2_aaaa(a,b,j,k)
    //              += +0.50 P(i,j) <k,l||c,d>_abab r2_abab(c,d,i,l) t2_aaaa(a,b,j,k)
    //              += +0.50 P(i,j) <k,l||d,c>_abab r2_abab(d,c,i,l) t2_aaaa(a,b,j,k)
    //              += -1.00 P(i,j) <k,l||j,c>_abab r1_bb(c,l) t2_aaaa(a,b,i,k)
    //              += -1.00 P(i,j) <l,k||c,j>_aaaa r1_aa(c,l) t2_aaaa(a,b,i,k)
    //              += -1.00 P(i,j) d-_aa(k,c) r1_1_aa(c,i) t2_aaaa(a,b,j,k)
    //              += -1.00 P(i,j) d-_aa(k,c) r1_aa(c,i) t2_1_aaaa(a,b,j,k)
    //              += +1.00 P(i,j) d-_aa(k,j) r2_1_aaaa(a,b,i,k)
    //              += +0.50 P(i,j) <l,k||c,j>_aaaa r1_aa(c,i) t2_aaaa(a,b,l,k)
    //              += -0.50 P(i,j) <l,k||c,d>_aaaa r2_aaaa(a,b,i,l) t2_aaaa(c,d,j,k)
    //              += +1.00 P(i,j) d-_aa(k,j) r2_aaaa(a,b,i,k) t0_1
    //              += -1.00 P(i,j) f_aa(k,j) r2_aaaa(a,b,i,k)
    //              += -0.50 P(i,j) <l,k||c,d>_abab r2_aaaa(a,b,i,l) t2_abab(c,d,j,k)
    //              += -0.50 P(i,j) <l,k||d,c>_abab r2_aaaa(a,b,i,l) t2_abab(d,c,j,k)
    //              += +1.00 P(i,j) <a,b||c,j>_aaaa r1_aa(c,i)
    sigmar2_aaaa("R,a,b,i,j")  = -1.00 * tmps_["38_aaaa_Lovvo"]("R,i,a,b,j");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["38_aaaa_Lovvo"]("R,j,a,b,i");
    tmps_["38_aaaa_Lovvo"].~TArrayD();

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["40_aaaa_Lvvoo"]("R,a,b,i,j")  = r0_1("R") * t2_1["aaaa"]("a,b,i,j");

    // csigmar2_1_aaaa += +1.00 r0_1 t2_1_aaaa(a,b,i,j)
    csigmar2_1_aaaa("R,a,b,i,j") += tmps_["40_aaaa_Lvvoo"]("R,a,b,i,j");

    // flops: o0v0L1  = o1v1L1 o1v1L1 o0v0L1
    //  mems: o0v0L1  = o0v0L1 o0v0L1 o0v0L1
    tmps_["43_L"]("R")  = dp["aa_ov"]("j,b") * r1_1["aa"]("R,b,j");
    tmps_["43_L"]("R") += dp["bb_ov"]("i,a") * r1_1["bb"]("R,a,i");

    // flops: o4v0L1  = o4v2L1
    //  mems: o4v0L1  = o4v0L1
    tmps_["42_aaaa_Loooo"]("R,i,j,k,l")  = r2["aaaa"]("R,a,b,i,j") * eri["aaaa_oovv"]("k,l,a,b");

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["41_aaaa_Lvvoo"]("R,a,b,i,j")  = -4.00 * scalars_["6"] * r2["aaaa"]("R,a,b,i,j");

    // flops: o4v0L1  = o4v2L1
    //  mems: o4v0L1  = o4v0L1
    tmps_["39_aaaa_Loooo"]("R,i,j,k,l")  = r2_1["aaaa"]("R,a,b,i,j") * eri["aaaa_oovv"]("k,l,a,b");

    // flops: o2v2L1  = 0 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v4L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = 0 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j")  = (scalars_["2"] + scalars_["1"]) * tmps_["40_aaaa_Lvvoo"]("R,a,b,i,j");
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j") += t2_1["aaaa"]("a,b,i,j") * tmps_["43_L"]("R");
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j") += 0.25 * t2["aaaa"]("a,b,k,l") * tmps_["39_aaaa_Loooo"]("R,i,j,l,k");
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j") += 0.25 * t2_1["aaaa"]("a,b,k,l") * tmps_["42_aaaa_Loooo"]("R,i,j,l,k");
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j") += 0.25 * (reused_["25_aaaa_oovv"]("i,j,a,b") + -4.00 * eri["aaaa_vvoo"]("a,b,i,j")) * r0_1("R");
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j") += 0.25 * r2["aaaa"]("R,a,b,k,l") * reused_["22_aaaa_oooo"]("i,j,l,k");
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j") += scalars_["5"] * r2_1["aaaa"]("R,a,b,i,j");
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j") += 2.00 * scalars_["3"] * r2_1["aaaa"]("R,a,b,i,j");
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j") += scalars_["1"] * r2["aaaa"]("R,a,b,i,j");
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j") -= 0.50 * r0_1("R") * reused_["24_aaaa_vvoo"]("a,b,i,j");
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j") -= 0.50 * r0_1("R") * reused_["20_aaaa_vvoo"]("a,b,i,j");
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j") += 0.50 * eri["aaaa_oooo"]("l,k,i,j") * r2_1["aaaa"]("R,a,b,k,l");
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j") += 2.00 * scalars_["4"] * r2_1["aaaa"]("R,a,b,i,j");
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j") += scalars_["2"] * r2["aaaa"]("R,a,b,i,j");
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j") -= w0 * r2_1["aaaa"]("R,a,b,i,j");
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j") -= 0.50 * eri["aaaa_vvvv"]("a,b,c,d") * r2_1["aaaa"]("R,c,d,i,j");
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j") += 0.25 * r2_1["aaaa"]("R,a,b,k,l") * reused_["21_aaaa_oooo"]("l,k,i,j");
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j") -= 0.50 * r0_1("R") * reused_["19_aaaa_voov"]("a,i,j,b");
    tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j") += 0.25 * tmps_["41_aaaa_Lvvoo"]("R,a,b,i,j");
    tmps_["41_aaaa_Lvvoo"].~TArrayD();
    tmps_["40_aaaa_Lvvoo"].~TArrayD();
    tmps_["39_aaaa_Loooo"].~TArrayD();

    // sigmar2_1_aaaa += -1.00 d-_bb(k,k) r0_1 t2_1_aaaa(a,b,i,j)
    //                += -1.00 d-_aa(k,k) r0_1 t2_1_aaaa(a,b,i,j)
    //                += -1.00 d-_bb(k,c) r1_1_bb(c,k) t2_1_aaaa(a,b,i,j)
    //                += -1.00 d-_aa(k,c) r1_1_aa(c,k) t2_1_aaaa(a,b,i,j)
    //                += +0.25 <l,k||c,d>_aaaa r2_1_aaaa(c,d,i,j) t2_aaaa(a,b,l,k)
    //                += +0.25 <l,k||c,d>_aaaa r2_aaaa(c,d,i,j) t2_1_aaaa(a,b,l,k)
    //                += +0.25 <l,k||c,d>_aaaa r0_1 t2_aaaa(a,b,l,k) t2_aaaa(c,d,i,j)
    //                += +0.50 <l,k||i,j>_aaaa r0_1 t2_aaaa(a,b,l,k)
    //                += +1.00 <a,b||i,j>_aaaa r0_1
    //                += +0.25 <l,k||c,d>_aaaa r2_aaaa(a,b,l,k) t2_1_aaaa(c,d,i,j)
    //                += -1.00 d-_aa(k,k) r2_1_aaaa(a,b,i,j) t0_1
    //                += +1.00 f_aa(k,k) l2_1_abab(i,j,a,b)
    //                += -0.50 <l,k||l,k>_abab l2_1_abab(i,j,a,b)
    //                += -0.50 <k,l||k,l>_abab l2_1_abab(i,j,a,b)
    //                += +0.25 <l,k||c,d>_aaaa t2_aaaa(c,d,l,k) l2_1_abab(i,j,a,b)
    //                += -0.50 <l,k||l,k>_bbbb l2_1_abab(i,j,a,b)
    //                += -0.50 <l,k||l,k>_aaaa l2_1_abab(i,j,a,b)
    //                += +0.25 <l,k||c,d>_bbbb t2_bbbb(c,d,l,k) l2_1_abab(i,j,a,b)
    //                += +1.00 f_bb(k,k) l2_1_abab(i,j,a,b)
    //                += -1.00 d-_bb(k,k) t0_1 l2_1_abab(i,j,a,b)
    //                += +0.25 <l,k||c,d>_abab t2_abab(c,d,l,k) l2_1_abab(i,j,a,b)
    //                += +0.25 <k,l||c,d>_abab t2_abab(c,d,k,l) l2_1_abab(i,j,a,b)
    //                += +0.25 <l,k||d,c>_abab t2_abab(d,c,l,k) l2_1_abab(i,j,a,b)
    //                += +0.25 <k,l||d,c>_abab t2_abab(d,c,k,l) l2_1_abab(i,j,a,b)
    //                += -2.00 d-_bb(k,c) r2_1_aaaa(a,b,i,j) t1_1_bb(c,k)
    //                += -1.00 d+_aa(k,k) r2_aaaa(a,b,i,j)
    //                += -0.50 <l,k||c,d>_aaaa r0_1 t2_aaaa(c,a,l,k) t2_aaaa(d,b,i,j)
    //                += +0.50 <a,b||c,d>_aaaa r0_1 t2_aaaa(c,d,i,j)
    //                += +0.50 <l,k||i,j>_aaaa r2_1_aaaa(a,b,l,k)
    //                += -2.00 d-_aa(k,c) r2_1_aaaa(a,b,i,j) t1_1_aa(c,k)
    //                += -1.00 d+_bb(k,k) r2_aaaa(a,b,i,j)
    //                += +1.00 r2_1_aaaa(a,b,i,j) w0
    //                += +0.50 <a,b||c,d>_aaaa r2_1_aaaa(c,d,i,j)
    //                += +0.25 <l,k||c,d>_aaaa r2_1_aaaa(a,b,l,k) t2_aaaa(c,d,i,j)
    //                += -0.50 <l,k||c,d>_aaaa r0_1 t2_aaaa(c,a,i,j) t2_aaaa(d,b,l,k)
    //                += +1.00 r2_aaaa(a,b,i,j) t0_1 w0
    //                += +0.25 <l,k||c,d>_aaaa r2_aaaa(a,b,i,j) t2_1_aaaa(c,d,l,k)
    //                += +0.25 <l,k||c,d>_abab r2_aaaa(a,b,i,j) t2_1_abab(c,d,l,k)
    //                += +0.25 <k,l||c,d>_abab r2_aaaa(a,b,i,j) t2_1_abab(c,d,k,l)
    //                += +0.25 <l,k||d,c>_abab r2_aaaa(a,b,i,j) t2_1_abab(d,c,l,k)
    //                += +0.25 <k,l||d,c>_abab r2_aaaa(a,b,i,j) t2_1_abab(d,c,k,l)
    //                += -1.00 d-_aa(k,c) r2_aaaa(a,b,i,j) t0_1 t1_1_aa(c,k)
    //                += +1.00 f_bb(k,c) r2_aaaa(a,b,i,j) t1_1_bb(c,k)
    //                += -1.00 d-_bb(k,c) r2_aaaa(a,b,i,j) t0_1 t1_1_bb(c,k)
    //                += +0.25 <l,k||c,d>_bbbb r2_aaaa(a,b,i,j) t2_1_bbbb(c,d,l,k)
    //                += +1.00 f_aa(k,c) r2_aaaa(a,b,i,j) t1_1_aa(c,k)
    sigmar2_1_aaaa("R,a,b,i,j") -= tmps_["44_aaaa_Lvvoo"]("R,a,b,i,j");
    tmps_["44_aaaa_Lvvoo"].~TArrayD();

    // flops: o0v0L1  = o1v1L1 o1v1L1 o0v0L1
    //  mems: o0v0L1  = o0v0L1 o0v0L1 o0v0L1
    tmps_["45_L"]("R")  = dp["aa_ov"]("j,b") * r1["aa"]("R,b,j");
    tmps_["45_L"]("R") += dp["bb_ov"]("i,a") * r1["bb"]("R,a,i");

    // sigmar0_1 += -1.00 d+_bb(i,a) r1_bb(a,i)
    //           += -1.00 d+_aa(i,a) r1_aa(a,i)
    sigmar0_1("R") -= tmps_["45_L"]("R");

    // flops: o1v1L1  = o1v1L1 o3v2L1 o2v3L1 o2v1L1 o1v2L1 o2v2L1 o1v2L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o2v3L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    tmps_["46_bb_Lvo"]("R,a,i")  = t1_1["bb"]("a,i") * tmps_["45_L"]("R");
    tmps_["46_bb_Lvo"]("R,a,i") -= 0.50 * eri["bbbb_oovo"]("j,k,b,i") * r2["bbbb"]("R,b,a,k,j");
    tmps_["46_bb_Lvo"]("R,a,i") -= 0.50 * eri["bbbb_vovv"]("a,j,b,e") * r2["bbbb"]("R,b,e,i,j");
    tmps_["46_bb_Lvo"]("R,a,i") -= 0.50 * r1["bb"]("R,a,k") * reused_["31_bb_oo"]("i,k");
    tmps_["46_bb_Lvo"]("R,a,i") -= 0.50 * reused_["29_bb_vv"]("a,e") * r1["bb"]("R,e,i");
    tmps_["46_bb_Lvo"]("R,a,i") += r1["bb"]("R,e,k") * reused_["28_bbbb_voov"]("a,i,k,e");
    tmps_["46_bb_Lvo"]("R,a,i") += r1["bb"]("R,b,i") * reused_["34_bb_vv"]("a,b");
    tmps_["46_bb_Lvo"]("R,a,i") += scalars_["4"] * r1["bb"]("R,a,i");
    tmps_["46_bb_Lvo"]("R,a,i") += scalars_["5"] * r1["bb"]("R,a,i");
    tmps_["46_bb_Lvo"]("R,a,i") += dp["bb_vo"]("a,i") * r0_1("R");
    tmps_["46_bb_Lvo"]("R,a,i") += eri["baab_vovv"]("a,l,c,e") * r2["abab"]("R,c,e,l,i");
    tmps_["46_bb_Lvo"]("R,a,i") -= r1["bb"]("R,a,j") * reused_["35_bb_oo"]("j,i");
    tmps_["46_bb_Lvo"]("R,a,i") -= r1["bb"]("R,e,k") * reused_["26_bbbb_voov"]("a,i,k,e");
    tmps_["46_bb_Lvo"]("R,a,i") -= r2["bbbb"]("R,b,a,i,j") * reused_["12_bb_ov"]("j,b");
    tmps_["46_bb_Lvo"]("R,a,i") += r1["bb"]("R,a,k") * reused_["30_bb_oo"]("k,i");
    tmps_["46_bb_Lvo"]("R,a,i") += scalars_["1"] * r1_1["bb"]("R,a,i");
    tmps_["46_bb_Lvo"]("R,a,i") += r0_1("R") * reused_["32_bb_vo"]("a,i");
    tmps_["46_bb_Lvo"]("R,a,i") += eri["baab_vovo"]("a,l,c,i") * r1["aa"]("R,c,l");
    tmps_["46_bb_Lvo"]("R,a,i") += f["bb_oo"]("j,i") * r1["bb"]("R,a,j");
    tmps_["46_bb_Lvo"]("R,a,i") += r2["abab"]("R,c,a,l,i") * reused_["13_aa_ov"]("l,c");
    tmps_["46_bb_Lvo"]("R,a,i") -= f["bb_vv"]("a,b") * r1["bb"]("R,b,i");
    tmps_["46_bb_Lvo"]("R,a,i") -= r0_1("R") * reused_["33_bb_vo"]("a,i");
    tmps_["46_bb_Lvo"]("R,a,i") += f["bb_ov"]("j,b") * r2["bbbb"]("R,b,a,i,j");
    tmps_["46_bb_Lvo"]("R,a,i") += r1["bb"]("R,e,i") * reused_["3_bb_vv"]("e,a");
    tmps_["46_bb_Lvo"]("R,a,i") += r1["aa"]("R,d,m") * reused_["1_bbaa_voov"]("a,i,m,d");
    tmps_["46_bb_Lvo"]("R,a,i") -= r1["aa"]("R,d,m") * reused_["27_bbaa_voov"]("a,i,m,d");
    tmps_["46_bb_Lvo"]("R,a,i") += scalars_["2"] * r1_1["bb"]("R,a,i");
    tmps_["46_bb_Lvo"]("R,a,i") -= f["aa_ov"]("l,c") * r2["abab"]("R,c,a,l,i");
    tmps_["46_bb_Lvo"]("R,a,i") += eri["bbbb_vovo"]("a,j,b,i") * r1["bb"]("R,b,j");
    tmps_["46_bb_Lvo"]("R,a,i") += eri["abab_oovo"]("m,j,c,i") * r2["abab"]("R,c,a,m,j");
    tmps_["46_bb_Lvo"]("R,a,i") += scalars_["3"] * r1["bb"]("R,a,i");

    // sigmar1_bb  = -1.00 d-_bb(a,b) r1_bb(b,i) t0_1
    //            += +1.00 <k,j||b,c>_bbbb r1_bb(c,k) t2_bbbb(b,a,i,j)
    //            += -1.00 d-_aa(j,b) r1_bb(a,i) t1_1_aa(b,j)
    //            += -1.00 d-_aa(j,j) r1_bb(a,i) t0_1
    //                += +1.00 f_aa(k,k) l2_1_bbbb(i,j,a,b)
    //                += -0.50 <l,k||l,k>_abab l2_1_bbbb(i,j,a,b)
    //                += -0.50 <k,l||k,l>_abab l2_1_bbbb(i,j,a,b)
    //                += +0.25 <l,k||c,d>_aaaa t2_aaaa(c,d,l,k) l2_1_bbbb(i,j,a,b)
    //                += -0.50 <l,k||l,k>_bbbb l2_1_bbbb(i,j,a,b)
    //                += -0.50 <l,k||l,k>_aaaa l2_1_bbbb(i,j,a,b)
    //                += +0.25 <l,k||c,d>_bbbb t2_bbbb(c,d,l,k) l2_1_bbbb(i,j,a,b)
    //                += +1.00 f_bb(k,k) l2_1_bbbb(i,j,a,b)
    //                += -1.00 d-_bb(k,k) t0_1 l2_1_bbbb(i,j,a,b)
    //                += +0.25 <l,k||c,d>_abab t2_abab(c,d,l,k) l2_1_bbbb(i,j,a,b)
    //                += +0.25 <k,l||c,d>_abab t2_abab(c,d,k,l) l2_1_bbbb(i,j,a,b)
    //                += +0.25 <l,k||d,c>_abab t2_abab(d,c,l,k) l2_1_bbbb(i,j,a,b)
    //                += +0.25 <k,l||d,c>_abab t2_abab(d,c,k,l) l2_1_bbbb(i,j,a,b)
    //            += -0.50 <k,j||b,c>_bbbb r1_bb(c,i) t2_bbbb(b,a,k,j)
    //            += -1.00 d-_bb(a,i) r0_1
    //            += +0.50 <j,a||b,c>_abab r2_abab(b,c,j,i)
    //            += +0.50 <j,a||c,b>_abab r2_abab(c,b,j,i)
    //            += +1.00 d-_bb(j,i) r1_bb(a,j) t0_1
    //            += +1.00 <j,k||b,c>_abab r1_bb(c,k) t2_abab(b,a,j,i)
    //            += +1.00 d-_bb(j,b) r2_bbbb(b,a,i,j) t0_1
    //            += -0.50 <k,j||b,c>_bbbb r1_bb(a,k) t2_bbbb(b,c,i,j)
    //            += -0.50 <j,a||b,c>_bbbb r2_bbbb(b,c,i,j)
    //            += -0.50 <j,k||b,c>_abab r1_bb(a,k) t2_abab(b,c,j,i)
    //            += -0.50 <j,k||c,b>_abab r1_bb(a,k) t2_abab(c,b,j,i)
    //            += -1.00 d-_aa(j,j) r1_1_bb(a,i)
    //            += -1.00 d-_aa(j,b) r0_1 t2_abab(b,a,j,i)
    //            += +1.00 <j,a||b,i>_abab r1_aa(b,j)
    //            += -1.00 f_bb(j,i) r1_bb(a,j)
    //            += -1.00 d-_aa(j,b) r2_abab(b,a,j,i) t0_1
    //            += +1.00 f_bb(a,b) r1_bb(b,i)
    //            += +1.00 d-_bb(j,b) r0_1 t2_bbbb(b,a,i,j)
    //            += -1.00 f_bb(j,b) r2_bbbb(b,a,i,j)
    //            += -0.50 <k,j||b,c>_abab r1_bb(c,i) t2_abab(b,a,k,j)
    //            += -0.50 <j,k||b,c>_abab r1_bb(c,i) t2_abab(b,a,j,k)
    //            += -1.00 <k,j||c,b>_abab r1_aa(c,k) t2_bbbb(b,a,i,j)
    //            += -1.00 <k,j||b,c>_aaaa r1_aa(c,k) t2_abab(b,a,j,i)
    //            += -1.00 d-_bb(j,j) r1_1_bb(a,i)
    //            += +1.00 f_aa(j,b) r2_abab(b,a,j,i)
    //            += +1.00 <j,a||b,i>_bbbb r1_bb(b,j)
    //            += -0.50 <k,j||b,i>_abab r2_abab(b,a,k,j)
    //            += -0.50 <j,k||b,i>_abab r2_abab(b,a,j,k)
    //            += -1.00 d-_bb(j,b) r1_bb(a,i) t1_1_bb(b,j)
    //            += -0.50 <k,j||b,i>_bbbb r2_bbbb(b,a,k,j)
    //            += -1.00 d-_bb(j,b) r1_bb(b,j) t1_1_bb(a,i)
    //            += -1.00 d-_aa(j,b) r1_aa(b,j) t1_1_bb(a,i)
    sigmar1_bb("R,a,i")  = -1.00 * tmps_["46_bb_Lvo"]("R,a,i");
    tmps_["46_bb_Lvo"].~TArrayD();

    // flops: o4v0L1  = o4v2L1
    //  mems: o4v0L1  = o4v0L1
    tmps_["47_bbbb_Loooo"]("R,i,j,k,l")  = r2["bbbb"]("R,a,b,i,j") * eri["bbbb_oovv"]("k,l,a,b");

    // flops: o2v2L1  = o2v4L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["48_bbbb_Lvvoo"]("R,a,b,i,j")  = -0.50 * eri["bbbb_vvvv"]("a,b,c,d") * r2["bbbb"]("R,c,d,i,j");
    tmps_["48_bbbb_Lvvoo"]("R,a,b,i,j") += 0.25 * r2["bbbb"]("R,a,b,k,l") * reused_["36_bbbb_oooo"]("l,k,i,j");
    tmps_["48_bbbb_Lvvoo"]("R,a,b,i,j") += scalars_["4"] * r2["bbbb"]("R,a,b,i,j");
    tmps_["48_bbbb_Lvvoo"]("R,a,b,i,j") += scalars_["1"] * r2_1["bbbb"]("R,a,b,i,j");
    tmps_["48_bbbb_Lvvoo"]("R,a,b,i,j") += scalars_["3"] * r2["bbbb"]("R,a,b,i,j");
    tmps_["48_bbbb_Lvvoo"]("R,a,b,i,j") += scalars_["5"] * r2["bbbb"]("R,a,b,i,j");
    tmps_["48_bbbb_Lvvoo"]("R,a,b,i,j") += 0.50 * eri["bbbb_oooo"]("l,k,i,j") * r2["bbbb"]("R,a,b,k,l");
    tmps_["48_bbbb_Lvvoo"]("R,a,b,i,j") += scalars_["2"] * r2_1["bbbb"]("R,a,b,i,j");
    tmps_["48_bbbb_Lvvoo"]("R,a,b,i,j") += t2_1["bbbb"]("a,b,i,j") * tmps_["45_L"]("R");
    tmps_["48_bbbb_Lvvoo"]("R,a,b,i,j") += 0.25 * t2["bbbb"]("a,b,k,l") * tmps_["47_bbbb_Loooo"]("R,i,j,l,k");

    // sigmar2_bbbb  = -1.00 d-_aa(k,c) r2_bbbb(a,b,i,j) t1_1_aa(c,k)
    //              += +0.25 <l,k||c,d>_bbbb r2_bbbb(a,b,l,k) t2_bbbb(c,d,i,j)
    //              += -1.00 d-_aa(k,k) r2_1_bbbb(a,b,i,j)
    //              += -1.00 d-_bb(k,c) r2_bbbb(a,b,i,j) t1_1_bb(c,k)
    //              += -1.00 d-_aa(k,k) r2_bbbb(a,b,i,j) t0_1
    //              += +1.00 f_aa(j,j) l1_1_aa(i,a)
    //              += -0.50 <k,j||k,j>_abab l1_1_aa(i,a)
    //              += -0.50 <j,k||j,k>_abab l1_1_aa(i,a)
    //              += +0.25 <k,j||b,c>_aaaa t2_aaaa(b,c,k,j) l1_1_aa(i,a)
    //              += -0.50 <k,j||k,j>_bbbb l1_1_aa(i,a)
    //              += -0.50 <k,j||k,j>_aaaa l1_1_aa(i,a)
    //              += +0.25 <k,j||b,c>_bbbb t2_bbbb(b,c,k,j) l1_1_aa(i,a)
    //              += +1.00 f_bb(j,j) l1_1_aa(i,a)
    //              += -1.00 d-_bb(j,j) t0_1 l1_1_aa(i,a)
    //              += +0.25 <k,j||b,c>_abab t2_abab(b,c,k,j) l1_1_aa(i,a)
    //              += +0.25 <j,k||b,c>_abab t2_abab(b,c,j,k) l1_1_aa(i,a)
    //              += +0.25 <k,j||c,b>_abab t2_abab(c,b,k,j) l1_1_aa(i,a)
    //              += +0.25 <j,k||c,b>_abab t2_abab(c,b,j,k) l1_1_aa(i,a)
    //              += +0.50 <a,b||c,d>_bbbb r2_bbbb(c,d,i,j)
    //              += +0.50 <l,k||i,j>_bbbb r2_bbbb(a,b,l,k)
    //              += -1.00 d-_bb(k,k) r2_1_bbbb(a,b,i,j)
    //              += -1.00 d-_bb(k,c) r1_bb(c,k) t2_1_bbbb(a,b,i,j)
    //              += -1.00 d-_aa(k,c) r1_aa(c,k) t2_1_bbbb(a,b,i,j)
    //              += +0.25 <l,k||c,d>_bbbb r2_bbbb(c,d,i,j) t2_bbbb(a,b,l,k)
    sigmar2_bbbb("R,a,b,i,j")  = -1.00 * tmps_["48_bbbb_Lvvoo"]("R,a,b,i,j");
    tmps_["48_bbbb_Lvvoo"].~TArrayD();

    // flops: o2v0L1  = o2v1L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["49_bb_Loo"]("R,i,j")  = r1_1["bb"]("R,a,i") * dp["bb_ov"]("j,a");

    // flops: o2v2L1  = o3v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["50_bbbb_Lvvoo"]("R,a,b,i,j")  = -2.00 * t2_1["bbbb"]("a,b,i,k") * tmps_["49_bb_Loo"]("R,j,k");
    tmps_["50_bbbb_Lvvoo"].~TArrayD();

    // flops: o2v2L1  = o2v2L1 o3v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1
    tmps_["51_bbbb_Lovvo"]("R,i,a,b,j")  = r0_1("R") * reused_["39_bbbb_ovvo"]("i,a,b,j");
    tmps_["51_bbbb_Lovvo"]("R,i,a,b,j") += r2["bbbb"]("R,a,b,j,k") * reused_["38_bb_oo"]("k,i");

    // sigmar2_bbbb += +1.00 P(i,j) d-_bb(k,j) r0_1 t2_bbbb(a,b,i,k)
    //              += +1.00 P(i,j) d-_bb(k,c) r2_bbbb(a,b,i,k) t1_1_bb(c,j)
    sigmar2_bbbb("R,a,b,i,j") += tmps_["51_bbbb_Lovvo"]("R,j,a,b,i");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["51_bbbb_Lovvo"]("R,i,a,b,j");

    // flops: o2v0L1  = o2v1L1 o2v1L1 o2v1L1 o2v0L1 o2v0L1 o2v1L1 o2v0L1 o3v2L1 o2v0L1 o3v1L1 o2v0L1 o3v2L1 o2v0L1 o3v1L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1
    tmps_["57_bb_Loo"]("R,i,j")  = -1.00 * r1["bb"]("R,c,i") * reused_["46_bb_ov"]("j,c");
    tmps_["57_bb_Loo"]("R,i,j") += r1["bb"]("R,c,i") * reused_["8_bb_ov"]("j,c");
    tmps_["57_bb_Loo"]("R,i,j") += f["bb_ov"]("j,b") * r1_1["bb"]("R,b,i");
    tmps_["57_bb_Loo"]("R,i,j") -= r1_1["bb"]("R,b,i") * reused_["12_bb_ov"]("j,b");
    tmps_["57_bb_Loo"]("R,i,j") += eri["abab_oovv"]("k,j,d,c") * r2_1["abab"]("R,d,c,k,i");
    tmps_["57_bb_Loo"]("R,i,j") += r1["bb"]("R,c,l") * reused_["44_bbbb_ooov"]("i,j,l,c");
    tmps_["57_bb_Loo"]("R,i,j") += 0.50 * eri["bbbb_oovv"]("j,l,b,c") * r2_1["bbbb"]("R,b,c,i,l");
    tmps_["57_bb_Loo"]("R,i,j") += r1["aa"]("R,a,k") * reused_["43_abab_oovo"]("k,j,a,i");

    // flops: o2v0L1  = o3v2L1 o2v1L1 o2v0L1 o3v2L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1
    tmps_["56_bb_Loo"]("R,i,j")  = 0.50 * r2["bbbb"]("R,c,b,i,l") * eri["bbbb_oovv"]("j,l,c,b");
    tmps_["56_bb_Loo"]("R,i,j") += r1["bb"]("R,c,i") * f["bb_ov"]("j,c");
    tmps_["56_bb_Loo"]("R,i,j") += r2["abab"]("R,a,b,k,i") * eri["abab_oovv"]("k,j,a,b");

    // flops: o2v0L1  = o3v1L1 o3v1L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
    tmps_["55_bb_Loo"]("R,i,j")  = -1.00 * eri["bbbb_oovo"]("i,l,b,j") * r1_1["bb"]("R,b,l");
    tmps_["55_bb_Loo"]("R,i,j") += eri["abab_oovo"]("k,i,a,j") * r1_1["aa"]("R,a,k");

    // flops: o2v0L1  = o3v1L1 o3v1L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
    tmps_["54_bb_Loo"]("R,i,j")  = -1.00 * eri["bbbb_oovo"]("i,l,b,j") * r1["bb"]("R,b,l");
    tmps_["54_bb_Loo"]("R,i,j") += eri["abab_oovo"]("k,i,a,j") * r1["aa"]("R,a,k");

    // flops: o2v0L1  = o2v1L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["53_bb_Loo"]("R,i,j")  = r1["bb"]("R,a,i") * reused_["12_bb_ov"]("j,a");

    // flops: o2v2L1  = o3v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["52_bbbb_Lvvoo"]("R,a,b,i,j")  = -1.00 * r2["bbbb"]("R,a,b,i,k") * reused_["40_bb_oo"]("k,j");

    // flops: o2v2L1  = o4v2L1 o3v2L1 o3v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j")  = -0.50 * r2["bbbb"]("R,a,b,m,k") * reused_["45_bbbb_oooo"]("k,m,i,j");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") -= r2_1["bbbb"]("R,a,b,j,m") * reused_["30_bb_oo"]("m,i");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") += 2.00 * r2_1["bbbb"]("R,a,b,j,k") * reused_["38_bb_oo"]("k,i");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") += eri["bbbb_vvvo"]("a,b,c,i") * r1_1["bb"]("R,c,j");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") -= r2["bbbb"]("R,a,b,j,m") * reused_["48_bb_oo"]("m,i");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") += 0.50 * r2_1["bbbb"]("R,a,b,j,m") * reused_["31_bb_oo"]("i,m");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") += 0.50 * r2["bbbb"]("R,a,b,j,m") * reused_["47_bb_oo"]("i,m");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") -= 0.50 * reused_["41_bbbb_vovv"]("c,i,a,b") * r1_1["bb"]("R,c,j");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") += 0.50 * r0_1("R") * reused_["49_bbbb_ovvo"]("i,a,b,j");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") -= r2_1["bbbb"]("R,a,b,j,k") * f["bb_oo"]("k,i");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") += r2_1["bbbb"]("R,a,b,j,k") * reused_["35_bb_oo"]("k,i");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") += r2["bbbb"]("R,a,b,j,k") * dp["bb_oo"]("k,i");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") -= 0.50 * reused_["42_bbbb_vovv"]("c,i,a,b") * r1["bb"]("R,c,j");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") -= r1["bb"]("R,e,j") * reused_["37_bbbb_vvvo"]("a,b,e,i");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") += tmps_["52_bbbb_Lvvoo"]("R,a,b,j,i");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") += t2["bbbb"]("a,b,i,k") * tmps_["57_bb_Loo"]("R,j,k");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") -= t2_1["bbbb"]("a,b,i,k") * tmps_["53_bb_Loo"]("R,j,k");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") += t0_1 * tmps_["51_bbbb_Lovvo"]("R,i,a,b,j");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") -= t2_1["bbbb"]("a,b,j,k") * tmps_["54_bb_Loo"]("R,k,i");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") += t2_1["bbbb"]("a,b,i,k") * tmps_["56_bb_Loo"]("R,j,k");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") -= t2["bbbb"]("a,b,j,k") * tmps_["55_bb_Loo"]("R,k,i");
    tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j") -= 2.00 * t2_1["bbbb"]("a,b,i,k") * tmps_["49_bb_Loo"]("R,j,k");
    tmps_["52_bbbb_Lvvoo"].~TArrayD();
    tmps_["51_bbbb_Lovvo"].~TArrayD();

    // sigmar2_1_bbbb += +1.00 P(i,j) <a,b||c,j>_bbbb r1_1_bb(c,i)
    //                += +2.00 P(i,j) d-_bb(k,c) r2_1_bbbb(a,b,i,k) t1_1_bb(c,j)
    //                += -0.50 P(i,j) <k,l||c,d>_abab r2_1_bbbb(a,b,i,l) t2_abab(c,d,k,j)
    //                += -0.50 P(i,j) <k,l||d,c>_abab r2_1_bbbb(a,b,i,l) t2_abab(d,c,k,j)
    //                += +0.50 P(i,j) <l,k||c,j>_bbbb r2_bbbb(a,b,l,k) t1_1_bb(c,i)
    //                += -1.00 P(i,j) <k,l||c,j>_abab r2_bbbb(a,b,i,l) t1_1_aa(c,k)
    //                += -0.50 P(i,j) <k,l||c,d>_abab r2_bbbb(a,b,i,l) t2_1_abab(c,d,k,j)
    //                += -0.50 P(i,j) <k,l||d,c>_abab r2_bbbb(a,b,i,l) t2_1_abab(d,c,k,j)
    //                += -0.50 P(i,j) <l,k||c,d>_bbbb r2_1_bbbb(a,b,i,l) t2_bbbb(c,d,j,k)
    //                += -0.50 P(i,j) <l,k||c,d>_bbbb r2_bbbb(a,b,i,l) t2_1_bbbb(c,d,j,k)
    //                += +1.00 P(i,j) <l,k||c,j>_bbbb r2_bbbb(a,b,i,l) t1_1_bb(c,k)
    //                += +0.50 P(i,j) <l,k||c,j>_bbbb r1_1_bb(c,i) t2_bbbb(a,b,l,k)
    //                += -0.50 P(i,j) <l,k||c,d>_bbbb r0_1 t2_bbbb(a,b,i,l) t2_bbbb(c,d,j,k)
    //                += -1.00 P(i,j) f_bb(k,j) r0_1 t2_bbbb(a,b,i,k)
    //                += -2.00 P(i,j) d-_bb(k,c) r0_1 t2_bbbb(a,b,j,k) t1_1_bb(c,i)
    //                += +2.00 P(i,j) d-_bb(k,j) r0_1 t2_1_bbbb(a,b,i,k)
    //                += +1.00 P(i,j) <l,k||c,d>_aaaa r0_1 t2_abab(c,a,k,j) t2_abab(d,b,l,i)
    //                += +1.00 P(i,j) <l,k||c,d>_bbbb r0_1 t2_bbbb(c,a,j,k) t2_bbbb(d,b,i,l)
    //                += -0.50 P(i,j) <k,l||c,d>_abab r0_1 t2_bbbb(a,b,i,l) t2_abab(c,d,k,j)
    //                += -0.50 P(i,j) <k,l||d,c>_abab r0_1 t2_bbbb(a,b,i,l) t2_abab(d,c,k,j)
    //                += -1.00 P(i,j) f_bb(k,j) r2_1_bbbb(a,b,i,k)
    //                += +1.00 P(i,j) d-_bb(k,j) r2_1_bbbb(a,b,i,k) t0_1
    //                += +1.00 P(i,j) d+_bb(k,j) r2_bbbb(a,b,i,k)
    //                += +0.50 P(i,j) <l,k||c,j>_bbbb r1_bb(c,i) t2_1_bbbb(a,b,l,k)
    //                += -1.00 P(i,j) <a,b||c,d>_bbbb r1_bb(d,i) t1_1_bb(c,j)
    //                += -0.50 P(i,j) <l,k||c,d>_bbbb r1_bb(d,i) t2_bbbb(a,b,l,k) t1_1_bb(c,j)
    //                += -1.00 P(i,j) f_bb(k,c) r2_bbbb(a,b,i,k) t1_1_bb(c,j)
    //                += +1.00 P(i,j) f_bb(k,c) r1_1_bb(c,i) t2_bbbb(a,b,j,k)
    //                += +1.00 P(i,j) <l,k||c,d>_abab r1_bb(d,i) t2_bbbb(a,b,j,k) t1_1_aa(c,l)
    //                += +1.00 P(i,j) <l,k||c,d>_bbbb r1_bb(d,i) t2_bbbb(a,b,j,k) t1_1_bb(c,l)
    //                += -1.00 P(i,j) d-_bb(k,c) r1_1_bb(c,i) t2_bbbb(a,b,j,k) t0_1
    //                += +0.50 P(i,j) <l,k||c,d>_abab r2_1_abab(c,d,l,i) t2_bbbb(a,b,j,k)
    //                += +0.50 P(i,j) <l,k||d,c>_abab r2_1_abab(d,c,l,i) t2_bbbb(a,b,j,k)
    //                += -1.00 P(i,j) <l,k||c,d>_bbbb r1_bb(d,l) t2_bbbb(a,b,j,k) t1_1_bb(c,i)
    //                += -0.50 P(i,j) <l,k||c,d>_bbbb r2_1_bbbb(c,d,i,l) t2_bbbb(a,b,j,k)
    //                += +1.00 P(i,j) <l,k||d,c>_abab r1_aa(d,l) t2_bbbb(a,b,j,k) t1_1_bb(c,i)
    //                += -1.00 P(i,j) d-_bb(k,c) r1_bb(c,i) t0_1 t2_1_bbbb(a,b,j,k)
    //                += +1.00 P(i,j) d-_bb(k,j) r0_1 t2_bbbb(a,b,i,k) t0_1
    //                += +1.00 P(i,j) d-_bb(k,c) r2_bbbb(a,b,i,k) t0_1 t1_1_bb(c,j)
    //                += -1.00 P(i,j) <l,k||c,j>_abab r1_aa(c,l) t2_1_bbbb(a,b,i,k)
    //                += -1.00 P(i,j) <l,k||c,j>_bbbb r1_bb(c,l) t2_1_bbbb(a,b,i,k)
    //                += +1.00 P(i,j) f_bb(k,c) r1_bb(c,i) t2_1_bbbb(a,b,j,k)
    //                += -0.50 P(i,j) <l,k||c,d>_bbbb r2_bbbb(c,d,i,l) t2_1_bbbb(a,b,j,k)
    //                += +0.50 P(i,j) <l,k||c,d>_abab r2_abab(c,d,l,i) t2_1_bbbb(a,b,j,k)
    //                += +0.50 P(i,j) <l,k||d,c>_abab r2_abab(d,c,l,i) t2_1_bbbb(a,b,j,k)
    //                += -1.00 P(i,j) <l,k||c,j>_abab r1_1_aa(c,l) t2_bbbb(a,b,i,k)
    //                += -1.00 P(i,j) <l,k||c,j>_bbbb r1_1_bb(c,l) t2_bbbb(a,b,i,k)
    //                += -2.00 P(i,j) d-_bb(k,c) r1_1_bb(c,i) t2_1_bbbb(a,b,j,k)
    sigmar2_1_bbbb("R,a,b,i,j") += tmps_["58_bbbb_Lvvoo"]("R,a,b,j,i");
    sigmar2_1_bbbb("R,a,b,i,j") -= tmps_["58_bbbb_Lvvoo"]("R,a,b,i,j");
    tmps_["58_bbbb_Lvvoo"].~TArrayD();

    // flops: o3v1L1  = o3v2L1
    //  mems: o3v1L1  = o3v1L1
    tmps_["60_bbbb_Loovo"]("L,i,j,a,k")  = l2_1["bbbb"]("L,i,j,a,b") * t1_1["bb"]("b,k");

    // flops: o4v0L1  = o4v2L1 o4v2L1 o4v0L1
    //  mems: o4v0L1  = o4v0L1 o4v0L1 o4v0L1
    tmps_["59_bbbb_Loooo"]("L,i,j,k,l")  = l2["bbbb"]("L,i,j,a,b") * t2["bbbb"]("a,b,k,l");
    tmps_["59_bbbb_Loooo"]("L,i,j,k,l") += l2_1["bbbb"]("L,i,j,a,b") * t2_1["bbbb"]("a,b,k,l");

    // flops: o2v2L1  = o4v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v4L1 o2v2L1 o4v2L1 o2v2L1 o3v3L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["61_bbbb_Loovv"]("L,i,j,a,b")  = -0.25 * l2["bbbb"]("L,k,l,a,b") * reused_["36_bbbb_oooo"]("i,j,k,l");
    tmps_["61_bbbb_Loovv"]("L,i,j,a,b") += scalars_["3"] * l2["bbbb"]("L,i,j,a,b");
    tmps_["61_bbbb_Loovv"]("L,i,j,a,b") -= l2_1["bbbb"]("L,k,l,a,b") * reused_["45_bbbb_oooo"]("i,j,l,k");
    tmps_["61_bbbb_Loovv"]("L,i,j,a,b") += scalars_["1"] * l2_1["bbbb"]("L,i,j,a,b");
    tmps_["61_bbbb_Loovv"]("L,i,j,a,b") -= scalars_["6"] * l2_1["bbbb"]("L,i,j,a,b");
    tmps_["61_bbbb_Loovv"]("L,i,j,a,b") -= 0.50 * eri["bbbb_oooo"]("i,j,k,l") * l2["bbbb"]("L,k,l,a,b");
    tmps_["61_bbbb_Loovv"]("L,i,j,a,b") += scalars_["4"] * l2["bbbb"]("L,i,j,a,b");
    tmps_["61_bbbb_Loovv"]("L,i,j,a,b") += scalars_["5"] * l2["bbbb"]("L,i,j,a,b");
    tmps_["61_bbbb_Loovv"]("L,i,j,a,b") -= 0.25 * l2_1["bbbb"]("L,k,l,a,b") * reused_["50_bbbb_oooo"]("k,l,i,j");
    tmps_["61_bbbb_Loovv"]("L,i,j,a,b") += scalars_["2"] * l2_1["bbbb"]("L,i,j,a,b");
    tmps_["61_bbbb_Loovv"]("L,i,j,a,b") += 0.50 * eri["bbbb_vvvv"]("d,c,a,b") * l2["bbbb"]("L,i,j,c,d");
    tmps_["61_bbbb_Loovv"]("L,i,j,a,b") += 0.25 * eri["bbbb_oovv"]("k,l,a,b") * tmps_["59_bbbb_Loooo"]("L,i,j,l,k");
    tmps_["61_bbbb_Loovv"]("L,i,j,a,b") += eri["bbbb_vovv"]("c,k,a,b") * tmps_["60_bbbb_Loovo"]("L,i,j,c,k");
    tmps_["59_bbbb_Loooo"].~TArrayD();

    // sigmal2_bbbb  = -1.00 d-_bb(k,c) t1_1_bb(c,k) l2_bbbb(i,j,a,b)
    //              += +0.25 <i,j||c,d>_bbbb t2_bbbb(c,d,k,l) l2_bbbb(k,l,a,b)
    //              += +1.00 <i,j||c,l>_bbbb t1_1_bb(c,k) l2_1_bbbb(k,l,a,b)
    //              += -1.00 d+_aa(k,k) l2_1_bbbb(i,j,a,b)
    //              += +1.00 t0_1 l2_1_bbbb(i,j,a,b) w0
    //              += +0.25 <l,k||c,d>_aaaa t2_1_aaaa(c,d,l,k) l2_1_bbbb(i,j,a,b)
    //              += +0.25 <l,k||c,d>_abab t2_1_abab(c,d,l,k) l2_1_bbbb(i,j,a,b)
    //              += +0.25 <k,l||c,d>_abab t2_1_abab(c,d,k,l) l2_1_bbbb(i,j,a,b)
    //              += +0.25 <l,k||d,c>_abab t2_1_abab(d,c,l,k) l2_1_bbbb(i,j,a,b)
    //              += +0.25 <k,l||d,c>_abab t2_1_abab(d,c,k,l) l2_1_bbbb(i,j,a,b)
    //              += -1.00 d-_aa(k,c) t0_1 t1_1_aa(c,k) l2_1_bbbb(i,j,a,b)
    //              += +1.00 f_bb(k,c) t1_1_bb(c,k) l2_1_bbbb(i,j,a,b)
    //              += -1.00 d-_bb(k,c) t0_1 t1_1_bb(c,k) l2_1_bbbb(i,j,a,b)
    //              += +0.25 <l,k||c,d>_bbbb t2_1_bbbb(c,d,l,k) l2_1_bbbb(i,j,a,b)
    //              += +1.00 f_aa(k,c) t1_1_aa(c,k) l2_1_bbbb(i,j,a,b)
    //              += +0.50 <i,j||k,l>_bbbb l2_bbbb(k,l,a,b)
    //              += -1.00 d-_aa(k,c) t1_1_aa(c,k) l2_bbbb(i,j,a,b)
    //              += -1.00 d-_aa(k,k) t0_1 l2_bbbb(i,j,a,b)
    //              += +1.00 f_aa(j,j) r1_1_aa(a,i)
    //              += -0.50 <k,j||k,j>_abab r1_1_aa(a,i)
    //              += -0.50 <j,k||j,k>_abab r1_1_aa(a,i)
    //              += +0.25 <k,j||b,c>_aaaa r1_1_aa(a,i) t2_aaaa(b,c,k,j)
    //              += -0.50 <k,j||k,j>_bbbb r1_1_aa(a,i)
    //              += -0.50 <k,j||k,j>_aaaa r1_1_aa(a,i)
    //              += +0.25 <k,j||b,c>_bbbb r1_1_aa(a,i) t2_bbbb(b,c,k,j)
    //              += +1.00 f_bb(j,j) r1_1_aa(a,i)
    //              += -1.00 d-_bb(j,j) r1_1_aa(a,i) t0_1
    //              += +0.25 <k,j||b,c>_abab r1_1_aa(a,i) t2_abab(b,c,k,j)
    //              += +0.25 <j,k||b,c>_abab r1_1_aa(a,i) t2_abab(b,c,j,k)
    //              += +0.25 <k,j||c,b>_abab r1_1_aa(a,i) t2_abab(c,b,k,j)
    //              += +0.25 <j,k||c,b>_abab r1_1_aa(a,i) t2_abab(c,b,j,k)
    //              += +0.25 <i,j||c,d>_bbbb t2_1_bbbb(c,d,k,l) l2_1_bbbb(k,l,a,b)
    //              += -1.00 d+_bb(k,k) l2_1_bbbb(i,j,a,b)
    //              += +0.50 <d,c||a,b>_bbbb l2_bbbb(i,j,d,c)
    //              += +0.25 <l,k||a,b>_bbbb t2_1_bbbb(d,c,l,k) l2_1_bbbb(i,j,d,c)
    //              += +0.25 <l,k||a,b>_bbbb t2_bbbb(d,c,l,k) l2_bbbb(i,j,d,c)
    //              += +1.00 <k,d||a,b>_bbbb t1_1_bb(c,k) l2_1_bbbb(i,j,d,c)
    sigmal2_bbbb("L,a,b,i,j")  = -1.00 * tmps_["61_bbbb_Loovv"]("L,i,j,a,b");
    tmps_["61_bbbb_Loovv"].~TArrayD();

    // flops: o2v2L1  = o3v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["62_bbbb_Lvoov"]("R,a,i,j,b")  = reused_["51_bbbb_vooo"]("a,i,j,k") * r1["bb"]("R,b,k");

    // sigmar2_1_bbbb += -1.00 P(a,b) d+_bb(k,c) r1_bb(b,k) t2_bbbb(c,a,i,j)
    sigmar2_1_bbbb("R,a,b,i,j") -= tmps_["62_bbbb_Lvoov"]("R,a,i,j,b");
    sigmar2_1_bbbb("R,a,b,i,j") += tmps_["62_bbbb_Lvoov"]("R,b,i,j,a");

    // flops: o0v2L1  = o1v3L1 o1v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["65_bb_Lvv"]("R,a,b")  = -1.00 * eri["bbbb_vovv"]("a,j,b,d") * r1["bb"]("R,d,j");
    tmps_["65_bb_Lvv"]("R,a,b") += eri["baab_vovv"]("a,i,c,b") * r1["aa"]("R,c,i");

    // flops: o0v2L1  = o2v3L1 o2v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["64_bb_Lvv"]("R,a,b")  = eri["bbbb_oovv"]("j,i,b,c") * r2["bbbb"]("R,c,a,i,j");
    tmps_["64_bb_Lvv"]("R,a,b") += r2["abab"]("R,d,a,k,j") * reused_["129_abab_oovv"]("k,j,d,b");

    // flops: o3v1L1  = o3v2L1
    //  mems: o3v1L1  = o3v1L1
    tmps_["63_bbbb_Lvooo"]("R,a,i,j,k")  = r2["bbbb"]("R,b,a,i,j") * dp["bb_ov"]("k,b");

    // flops: o2v2L1  = o2v3L1 o2v3L1 o2v3L1 o2v2L1 o3v2L1 o3v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["66_bbbb_Lvoov"]("R,a,i,j,b")  = -0.50 * reused_["29_bb_vv"]("a,d") * r2["bbbb"]("R,d,b,i,j");
    tmps_["66_bbbb_Lvoov"]("R,a,i,j,b") -= f["bb_vv"]("a,c") * r2["bbbb"]("R,c,b,i,j");
    tmps_["66_bbbb_Lvoov"]("R,a,i,j,b") += dp["bb_vv"]("a,c") * r2_1["bbbb"]("R,c,b,i,j");
    tmps_["66_bbbb_Lvoov"]("R,a,i,j,b") += r0_1("R") * reused_["52_bbbb_vvoo"]("a,b,i,j");
    tmps_["66_bbbb_Lvoov"]("R,a,i,j,b") -= r1["bb"]("R,b,k") * reused_["53_bbbb_vooo"]("a,i,j,k");
    tmps_["66_bbbb_Lvoov"]("R,a,i,j,b") += eri["bbbb_vooo"]("a,k,i,j") * r1["bb"]("R,b,k");
    tmps_["66_bbbb_Lvoov"]("R,a,i,j,b") += r2["bbbb"]("R,c,b,i,j") * reused_["34_bb_vv"]("a,c");
    tmps_["66_bbbb_Lvoov"]("R,a,i,j,b") += r2["bbbb"]("R,d,b,i,j") * reused_["3_bb_vv"]("d,a");
    tmps_["66_bbbb_Lvoov"]("R,a,i,j,b") += tmps_["65_bb_Lvv"]("R,a,c") * t2["bbbb"]("c,b,i,j");
    tmps_["66_bbbb_Lvoov"]("R,a,i,j,b") += t0_1 * tmps_["62_bbbb_Lvoov"]("R,a,i,j,b");
    tmps_["66_bbbb_Lvoov"]("R,a,i,j,b") -= 0.50 * t2["bbbb"]("c,a,i,j") * tmps_["64_bb_Lvv"]("R,b,c");
    tmps_["66_bbbb_Lvoov"]("R,a,i,j,b") -= t1_1["bb"]("a,k") * tmps_["63_bbbb_Lvooo"]("R,b,i,j,k");
    tmps_["62_bbbb_Lvoov"].~TArrayD();

    // sigmar2_bbbb += +1.00 P(a,b) <k,a||i,j>_bbbb r1_bb(b,k)
    //              += -1.00 P(a,b) d-_bb(a,c) r2_bbbb(c,b,i,j) t0_1
    //              += +1.00 P(a,b) f_bb(k,c) r1_bb(b,k) t2_bbbb(c,a,i,j)
    //              += +0.50 P(a,b) <k,a||c,d>_bbbb r1_bb(b,k) t2_bbbb(c,d,i,j)
    //              += -1.00 P(a,b) d-_bb(a,c) r0_1 t2_bbbb(c,b,i,j)
    //              += -0.50 P(a,b) <l,k||c,d>_abab r2_bbbb(d,b,i,j) t2_abab(c,a,l,k)
    //              += -0.50 P(a,b) <k,l||c,d>_abab r2_bbbb(d,b,i,j) t2_abab(c,a,k,l)
    //              += -1.00 P(a,b) d-_bb(a,c) r2_1_bbbb(c,b,i,j)
    //              += +1.00 P(a,b) f_bb(a,c) r2_bbbb(c,b,i,j)
    //              += -0.50 P(a,b) <l,k||c,d>_bbbb r2_bbbb(d,b,i,j) t2_bbbb(c,a,l,k)
    //              += +1.00 P(a,b) <k,a||d,c>_abab r1_aa(d,k) t2_bbbb(c,b,i,j)
    //              += -1.00 P(a,b) <k,a||c,d>_bbbb r1_bb(d,k) t2_bbbb(c,b,i,j)
    //              += -1.00 P(a,b) d-_bb(k,c) r1_bb(b,k) t2_bbbb(c,a,i,j) t0_1
    //              += -0.50 P(a,b) <l,k||c,d>_bbbb r2_bbbb(d,b,l,k) t2_bbbb(c,a,i,j)
    //              += +0.50 P(a,b) <l,k||d,c>_abab r2_abab(d,b,l,k) t2_bbbb(c,a,i,j)
    //              += +0.50 P(a,b) <k,l||d,c>_abab r2_abab(d,b,k,l) t2_bbbb(c,a,i,j)
    //              += +1.00 P(a,b) d-_bb(k,c) r2_bbbb(c,b,i,j) t1_1_bb(a,k)
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["66_bbbb_Lvoov"]("R,a,i,j,b");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["66_bbbb_Lvoov"]("R,b,i,j,a");
    tmps_["66_bbbb_Lvoov"].~TArrayD();
}
