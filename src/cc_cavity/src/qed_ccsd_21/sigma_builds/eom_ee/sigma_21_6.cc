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

void hilbert::EOM_EE_QED_CCSD::sigma_ee_21_6() {

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


    // flops: o1v1L1  = o1v1 o2v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o3v2L1 o1v1L1 o3v1L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o3v1L1 o1v1L1 o2v1L1 o1v1L1 o3v2L1 o1v1L1 o2v1L1 o1v1L1 o1v3L1 o1v1L1 o3v1L1 o1v1L1 o2v2L1 o1v1L1 o1v3L1 o1v1L1 o3v1L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o2v1L1 o1v1L1 o1v3L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o3v1L1 o1v1L1 o1v1L1 o1v1L1 o3v1L1 o1v1L1 o3v1L1 o1v1L1 o2v1L1 o1v1L1 o1v2L1 o1v1L1 o2v2L1 o1v1L1 o3v1L1 o1v1L1 o1v1L1 o1v2L1 o1v1L1 o2v3L1 o1v1L1 o2v3L1 o1v1L1 o1v1L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o1v2L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o2v3L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o1v2L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o2v1L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o2v1L1 o1v1L1 o3v2L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o1v1L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o2v1L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o1v1L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o1v3L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o3v1L1 o1v1L1 o3v1L1 o1v1L1 o3v1L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o1v1L1 o2v2L1 o2v1L1 o1v1L1 o1v1L1 o1v1L1 o4v1L1 o1v1L1 o4v1L1 o1v1L1 o2v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    tmps_["262_bb_Lov"]("L,i,a")  = -0.50 * (dp["bb_ov"]("j,a") + -1.00 * reused_["8_bb_ov"]("j,a")) * tmps_["179_bb_Loo"]("L,i,j");
    tmps_["262_bb_Lov"]("L,i,a") += t1_1["bb"]("f,n") * tmps_["221_bbbb_Lovvo"]("L,n,f,a,i");
    tmps_["262_bb_Lov"]("L,i,a") += dp["bb_ov"]("n,a") * tmps_["250_bb_Loo"]("L,i,n");
    tmps_["262_bb_Lov"]("L,i,a") += reused_["201_aabb_voov"]("c,l,m,a") * tmps_["169_abab_Loovo"]("L,l,i,c,m");
    tmps_["262_bb_Lov"]("L,i,a") -= reused_["111_bbbb_oovo"]("i,m,a,n") * tmps_["178_bb_Loo"]("L,n,m");
    tmps_["262_bb_Lov"]("L,i,a") -= 0.50 * t0_1 * tmps_["253_bb_Lvo"]("L,a,i");
    tmps_["262_bb_Lov"]("L,i,a") -= reused_["148_bb_ov"]("j,a") * tmps_["178_bb_Loo"]("L,i,j");
    tmps_["262_bb_Lov"]("L,i,a") += reused_["148_bb_ov"]("j,g") * tmps_["261_bbbb_Lvoov"]("L,g,j,i,a");
    tmps_["262_bb_Lov"]("L,i,a") += 0.50 * reused_["111_bbbb_oovo"]("i,m,a,n") * tmps_["179_bb_Loo"]("L,n,m");
    tmps_["262_bb_Lov"]("L,i,a") += 0.50 * reused_["148_bb_ov"]("j,a") * tmps_["179_bb_Loo"]("L,i,j");
    tmps_["262_bb_Lov"]("L,i,a") += reused_["202_bbbb_voov"]("g,n,m,a") * tmps_["60_bbbb_Loovo"]("L,i,n,g,m");
    tmps_["262_bb_Lov"]("L,i,a") -= f["bb_ov"]("n,a") * tmps_["222_bb_Loo"]("L,i,n");
    tmps_["262_bb_Lov"]("L,i,a") += 0.50 * eri["abab_vovv"]("c,i,e,a") * tmps_["240_aa_Lvv"]("L,c,e");
    tmps_["262_bb_Lov"]("L,i,a") -= eri["bbbb_oovo"]("i,m,a,j") * tmps_["217_bb_Loo"]("L,j,m");
    tmps_["262_bb_Lov"]("L,i,a") -= l2_1["abab"]("L,l,i,b,a") * reused_["155_aa_vo"]("b,l");
    tmps_["262_bb_Lov"]("L,i,a") -= eri["bbbb_vovv"]("g,i,a,d") * tmps_["114_bb_Lvv"]("L,g,d");
    tmps_["262_bb_Lov"]("L,i,a") -= 0.50 * eri["bbbb_oovo"]("i,m,a,j") * tmps_["243_bb_Loo"]("L,j,m");
    tmps_["262_bb_Lov"]("L,i,a") += reused_["26_bbbb_voov"]("g,n,m,a") * tmps_["60_bbbb_Loovo"]("L,i,n,g,m");
    tmps_["262_bb_Lov"]("L,i,a") += eri["abba_vovo"]("c,j,a,l") * tmps_["169_abab_Loovo"]("L,l,i,c,j");
    tmps_["262_bb_Lov"]("L,i,a") += reused_["6_aabb_voov"]("c,l,m,a") * tmps_["169_abab_Loovo"]("L,l,i,c,m");
    tmps_["262_bb_Lov"]("L,i,a") -= reused_["8_bb_ov"]("j,a") * tmps_["178_bb_Loo"]("L,i,j");
    tmps_["262_bb_Lov"]("L,i,a") -= 0.50 * eri["bbbb_vovv"]("g,i,a,d") * tmps_["241_bb_Lvv"]("L,g,d");
    tmps_["262_bb_Lov"]("L,i,a") -= f["bb_ov"]("j,a") * tmps_["218_bb_Loo"]("L,i,j");
    tmps_["262_bb_Lov"]("L,i,a") += dp["bb_ov"]("j,a") * tmps_["178_bb_Loo"]("L,i,j");
    tmps_["262_bb_Lov"]("L,i,a") += reused_["12_bb_ov"]("j,a") * tmps_["218_bb_Loo"]("L,i,j");
    tmps_["262_bb_Lov"]("L,i,a") += reused_["12_bb_ov"]("n,a") * tmps_["222_bb_Loo"]("L,i,n");
    tmps_["262_bb_Lov"]("L,i,a") -= eri["bbbb_oovo"]("i,j,a,n") * tmps_["222_bb_Loo"]("L,n,j");
    tmps_["262_bb_Lov"]("L,i,a") -= t0_1 * tmps_["254_bb_Lov"]("L,i,a");
    tmps_["262_bb_Lov"]("L,i,a") += eri["abba_oovo"]("o,i,a,l") * tmps_["209_aa_Loo"]("L,l,o");
    tmps_["262_bb_Lov"]("L,i,a") -= reused_["116_abba_oovo"]("k,i,a,l") * tmps_["205_aa_Loo"]("L,l,k");
    tmps_["262_bb_Lov"]("L,i,a") -= f["bb_ov"]("j,a") * tmps_["217_bb_Loo"]("L,i,j");
    tmps_["262_bb_Lov"]("L,i,a") -= dp["bb_vv"]("g,a") * tmps_["182_bb_Lov"]("L,i,g");
    tmps_["262_bb_Lov"]("L,i,a") += reused_["8_bb_ov"]("j,g") * tmps_["261_bbbb_Lvoov"]("L,g,j,i,a");
    tmps_["262_bb_Lov"]("L,i,a") += eri["abba_oovo"]("k,i,a,o") * tmps_["206_aa_Loo"]("L,o,k");
    tmps_["262_bb_Lov"]("L,i,a") -= dp["bb_ov"]("i,a") * l0_1("L");
    tmps_["262_bb_Lov"]("L,i,a") -= reused_["3_bb_vv"]("a,f") * l1["bb"]("L,i,f");
    tmps_["262_bb_Lov"]("L,i,a") -= 0.25 * reused_["41_bbbb_vovv"]("a,n,g,f") * l2["bbbb"]("L,i,n,g,f");
    tmps_["262_bb_Lov"]("L,i,a") += reused_["211_bbaa_vvvo"]("g,a,b,l") * l2_1["abab"]("L,l,i,b,g");
    tmps_["262_bb_Lov"]("L,i,a") -= scalars_["5"] * l1["bb"]("L,i,a");
    tmps_["262_bb_Lov"]("L,i,a") -= l2["abab"]("L,l,j,b,a") * reused_["225_aabb_vooo"]("b,l,j,i");
    tmps_["262_bb_Lov"]("L,i,a") += l1_1["bb"]("L,n,f") * reused_["102_bbbb_voov"]("f,n,i,a");
    tmps_["262_bb_Lov"]("L,i,a") -= eri["bbbb_vovo"]("f,i,a,n") * l1["bb"]("L,n,f");
    tmps_["262_bb_Lov"]("L,i,a") -= scalars_["3"] * l1["bb"]("L,i,a");
    tmps_["262_bb_Lov"]("L,i,a") -= dp["bb_vv"]("f,a") * l1_1["bb"]("L,i,f");
    tmps_["262_bb_Lov"]("L,i,a") -= reused_["215_baab_vovv"]("a,l,c,f") * l2["abab"]("L,l,i,c,f");
    tmps_["262_bb_Lov"]("L,i,a") -= 0.50 * l2["bbbb"]("L,n,j,f,a") * reused_["138_bbbb_vooo"]("f,n,j,i");
    tmps_["262_bb_Lov"]("L,i,a") -= l2["bbbb"]("L,n,j,f,a") * reused_["178_bbbb_oovo"]("i,j,f,n");
    tmps_["262_bb_Lov"]("L,i,a") -= reused_["216_baab_vovv"]("a,l,c,f") * l2_1["abab"]("L,l,i,c,f");
    tmps_["262_bb_Lov"]("L,i,a") += l1["bb"]("L,n,a") * reused_["35_bb_oo"]("i,n");
    tmps_["262_bb_Lov"]("L,i,a") += tmps_["251_bb_Lov"]("L,j,a") * dp["bb_oo"]("i,j");
    tmps_["262_bb_Lov"]("L,i,a") -= l1_1["bb"]("L,n,f") * reused_["164_bbbb_voov"]("f,n,i,a");
    tmps_["262_bb_Lov"]("L,i,a") -= l2_1["abab"]("L,l,i,b,a") * reused_["58_aa_vo"]("b,l");
    tmps_["262_bb_Lov"]("L,i,a") -= l1_1["aa"]("L,l,b") * reused_["114_aabb_voov"]("b,l,i,a");
    tmps_["262_bb_Lov"]("L,i,a") += 0.50 * l2_1["bbbb"]("L,n,j,f,a") * reused_["139_bbbb_vooo"]("f,n,j,i");
    tmps_["262_bb_Lov"]("L,i,a") += 0.50 * l2_1["bbbb"]("L,n,j,f,a") * reused_["140_bbbb_vooo"]("f,n,j,i");
    tmps_["262_bb_Lov"]("L,i,a") += 0.50 * l2["abab"]("L,l,i,b,a") * reused_["87_aa_vo"]("b,l");
    tmps_["262_bb_Lov"]("L,i,a") -= l2_1["abab"]("L,l,i,b,a") * reused_["153_aa_vo"]("b,l");
    tmps_["262_bb_Lov"]("L,i,a") -= 0.50 * l2_1["bbbb"]("L,i,n,f,a") * reused_["199_bb_ov"]("n,f");
    tmps_["262_bb_Lov"]("L,i,a") += reused_["112_bb_vv"]("f,a") * l1_1["bb"]("L,i,f");
    tmps_["262_bb_Lov"]("L,i,a") += l2_1["bbbb"]("L,n,j,f,a") * reused_["228_bbbb_vooo"]("f,j,i,n");
    tmps_["262_bb_Lov"]("L,i,a") += l1_1["aa"]("L,l,b") * reused_["162_aabb_voov"]("b,l,i,a");
    tmps_["262_bb_Lov"]("L,i,a") -= l2_1["bbbb"]("L,n,j,f,a") * reused_["222_bbbb_vooo"]("f,i,j,n");
    tmps_["262_bb_Lov"]("L,i,a") -= reused_["108_bb_vv"]("f,a") * l1_1["bb"]("L,i,f");
    tmps_["262_bb_Lov"]("L,i,a") -= l1_1["bb"]("L,n,a") * reused_["48_bb_oo"]("i,n");
    tmps_["262_bb_Lov"]("L,i,a") -= l2["abab"]("L,l,i,b,a") * reused_["152_aa_vo"]("b,l");
    tmps_["262_bb_Lov"]("L,i,a") -= f["bb_oo"]("i,n") * l1["bb"]("L,n,a");
    tmps_["262_bb_Lov"]("L,i,a") += l2["abab"]("L,l,j,b,a") * reused_["218_aabb_vooo"]("b,l,i,j");
    tmps_["262_bb_Lov"]("L,i,a") -= l2_1["abab"]("L,o,n,b,a") * reused_["226_baab_oovo"]("i,o,b,n");
    tmps_["262_bb_Lov"]("L,i,a") -= 0.25 * reused_["42_bbbb_vovv"]("a,n,g,f") * l2_1["bbbb"]("L,i,n,g,f");
    tmps_["262_bb_Lov"]("L,i,a") -= l1["aa"]("L,l,b") * reused_["6_aabb_voov"]("b,l,i,a");
    tmps_["262_bb_Lov"]("L,i,a") += reused_["213_aabb_vovv"]("b,l,g,a") * l2["abab"]("L,l,i,b,g");
    tmps_["262_bb_Lov"]("L,i,a") -= reused_["214_bbbb_vovv"]("f,n,g,a") * l2["bbbb"]("L,i,n,g,f");
    tmps_["262_bb_Lov"]("L,i,a") += l1_1["aa"]("L,l,b") * reused_["118_abba_vovo"]("b,i,a,l");
    tmps_["262_bb_Lov"]("L,i,a") -= 0.50 * reused_["109_bb_vv"]("a,f") * l1["bb"]("L,i,f");
    tmps_["262_bb_Lov"]("L,i,a") += l1["bb"]("L,n,a") * reused_["38_bb_oo"]("i,n");
    tmps_["262_bb_Lov"]("L,i,a") -= l2_1["abab"]("L,l,j,b,a") * reused_["227_aabb_vooo"]("b,l,j,i");
    tmps_["262_bb_Lov"]("L,i,a") += l2["abab"]("L,l,j,b,a") * reused_["208_aabb_vooo"]("b,l,j,i");
    tmps_["262_bb_Lov"]("L,i,a") -= 0.50 * l1["bb"]("L,n,a") * reused_["161_bb_oo"]("i,n");
    tmps_["262_bb_Lov"]("L,i,a") += l2_1["abab"]("L,o,n,b,a") * reused_["221_abab_vooo"]("b,i,o,n");
    tmps_["262_bb_Lov"]("L,i,a") -= scalars_["2"] * l1_1["bb"]("L,i,a");
    tmps_["262_bb_Lov"]("L,i,a") -= scalars_["4"] * l1["bb"]("L,i,a");
    tmps_["262_bb_Lov"]("L,i,a") += f["aa_vo"]("b,l") * l2["abab"]("L,l,i,b,a");
    tmps_["262_bb_Lov"]("L,i,a") += l0_1("L") * reused_["148_bb_ov"]("i,a");
    tmps_["262_bb_Lov"]("L,i,a") -= reused_["223_abba_vvvo"]("c,a,f,l") * l2_1["abab"]("L,l,i,c,f");
    tmps_["262_bb_Lov"]("L,i,a") += 0.50 * l2_1["bbbb"]("L,n,j,f,a") * reused_["207_bbbb_vooo"]("f,n,j,i");
    tmps_["262_bb_Lov"]("L,i,a") -= l2_1["abab"]("L,l,j,b,a") * reused_["224_abba_vooo"]("b,i,j,l");
    tmps_["262_bb_Lov"]("L,i,a") -= reused_["176_bbbb_vvvo"]("g,a,f,n") * l2["bbbb"]("L,i,n,g,f");
    tmps_["262_bb_Lov"]("L,i,a") += l2["bbbb"]("L,i,n,f,a") * reused_["91_bb_vo"]("f,n");
    tmps_["262_bb_Lov"]("L,i,a") -= 0.50 * eri["bbbb_vvvo"]("f,g,a,n") * l2["bbbb"]("L,i,n,g,f");
    tmps_["262_bb_Lov"]("L,i,a") += l2["bbbb"]("L,i,n,f,a") * reused_["84_bb_vo"]("f,n");
    tmps_["262_bb_Lov"]("L,i,a") += f["bb_vv"]("f,a") * l1["bb"]("L,i,f");
    tmps_["262_bb_Lov"]("L,i,a") -= scalars_["1"] * l1_1["bb"]("L,i,a");
    tmps_["262_bb_Lov"]("L,i,a") -= l2_1["abab"]("L,l,j,b,a") * reused_["206_aabb_vooo"]("b,l,j,i");
    tmps_["262_bb_Lov"]("L,i,a") -= l2["abab"]("L,l,i,b,a") * reused_["90_aa_vo"]("b,l");
    tmps_["262_bb_Lov"]("L,i,a") -= reused_["212_bbbb_vvvo"]("g,a,f,n") * l2_1["bbbb"]("L,i,n,g,f");
    tmps_["262_bb_Lov"]("L,i,a") += reused_["203_aabb_vovv"]("b,l,g,a") * l2_1["abab"]("L,l,i,b,g");
    tmps_["262_bb_Lov"]("L,i,a") += l2["abab"]("L,l,i,b,a") * reused_["89_aa_ov"]("l,b");
    tmps_["262_bb_Lov"]("L,i,a") -= l1_1["bb"]("L,n,a") * reused_["40_bb_oo"]("i,n");
    tmps_["262_bb_Lov"]("L,i,a") -= eri["abba_vvvo"]("c,f,a,l") * l2["abab"]("L,l,i,c,f");
    tmps_["262_bb_Lov"]("L,i,a") -= 0.50 * eri["bbbb_vooo"]("f,i,n,j") * l2["bbbb"]("L,n,j,f,a");
    tmps_["262_bb_Lov"]("L,i,a") += l1["bb"]("L,n,f") * reused_["26_bbbb_voov"]("f,n,i,a");
    tmps_["262_bb_Lov"]("L,i,a") -= eri["abab_vooo"]("b,i,l,j") * l2["abab"]("L,l,j,b,a");
    tmps_["262_bb_Lov"]("L,i,a") += 0.50 * l2["bbbb"]("L,n,j,f,a") * reused_["53_bbbb_vooo"]("f,n,j,i");
    tmps_["262_bb_Lov"]("L,i,a") -= l1["bb"]("L,n,f") * reused_["163_bbbb_voov"]("f,n,i,a");
    tmps_["262_bb_Lov"]("L,i,a") += scalars_["6"] * l1_1["bb"]("L,i,a");
    tmps_["262_bb_Lov"]("L,i,a") -= l1["bb"]("L,n,a") * reused_["30_bb_oo"]("i,n");
    tmps_["262_bb_Lov"]("L,i,a") += dp["bb_oo"]("i,n") * l1_1["bb"]("L,n,a");
    tmps_["262_bb_Lov"]("L,i,a") -= reused_["209_abba_vvvo"]("c,a,f,l") * l2["abab"]("L,l,i,c,f");
    tmps_["262_bb_Lov"]("L,i,a") -= f["bb_vo"]("f,n") * l2["bbbb"]("L,i,n,f,a");
    tmps_["262_bb_Lov"]("L,i,a") += reused_["210_bbaa_vvvo"]("g,a,b,l") * l2["abab"]("L,l,i,b,g");
    tmps_["262_bb_Lov"]("L,i,a") -= l2["bbbb"]("L,n,j,f,a") * reused_["220_bbbb_vooo"]("f,n,i,j");
    tmps_["262_bb_Lov"]("L,i,a") -= dp["aa_vo"]("b,l") * l2_1["abab"]("L,l,i,b,a");
    tmps_["262_bb_Lov"]("L,i,a") -= tmps_["252_bb_Lov"]("L,i,g") * dp["bb_vv"]("g,a");
    tmps_["262_bb_Lov"]("L,i,a") += l1_1["bb"]("L,n,a") * reused_["181_bb_oo"]("i,n");
    tmps_["262_bb_Lov"]("L,i,a") += l2["abab"]("L,l,j,b,a") * reused_["217_bbaa_oovo"]("i,j,b,l");
    tmps_["262_bb_Lov"]("L,i,a") -= eri["abba_vovo"]("b,i,a,l") * l1["aa"]("L,l,b");
    tmps_["262_bb_Lov"]("L,i,a") -= reused_["34_bb_vv"]("f,a") * l1["bb"]("L,i,f");
    tmps_["262_bb_Lov"]("L,i,a") += l2_1["abab"]("L,o,n,b,a") * reused_["204_bbaa_oovo"]("i,n,b,o");
    tmps_["262_bb_Lov"]("L,i,a") += l1["aa"]("L,l,b") * reused_["70_aabb_voov"]("b,l,i,a");
    tmps_["262_bb_Lov"]("L,i,a") -= 0.25 * reused_["205_bbbb_vovv"]("a,n,g,f") * l2_1["bbbb"]("L,i,n,g,f");
    tmps_["262_bb_Lov"]("L,i,a") += l2["bbbb"]("L,i,n,f,a") * reused_["88_bb_vo"]("f,n");
    tmps_["262_bb_Lov"]("L,i,a") -= l2["abab"]("L,o,n,b,a") * reused_["219_abba_vooo"]("b,n,i,o");
    tmps_["262_bb_Lov"]("L,i,a") -= tmps_["256_bb_Lov"]("L,i,a");
    tmps_["262_bb_Lov"]("L,i,a") += reused_["78_baab_voov"]("g,l,k,a") * tmps_["181_abba_Loovo"]("L,l,i,g,k");
    tmps_["262_bb_Lov"]("L,i,a") -= eri["bbbb_vovo"]("g,j,a,n") * tmps_["60_bbbb_Loovo"]("L,i,n,g,j");
    tmps_["262_bb_Lov"]("L,i,a") += eri["abab_vovv"]("c,i,e,a") * tmps_["140_aa_Lvv"]("L,c,e");
    tmps_["262_bb_Lov"]("L,i,a") -= dp["bb_ov"]("i,a") * tmps_["244_L"]("L");
    tmps_["262_bb_Lov"]("L,i,a") += l2_1["bbbb"]("L,i,n,f,a") * reused_["200_bb_vo"]("f,n");
    tmps_["262_bb_Lov"]("L,i,a") -= eri["baba_vovo"]("g,o,a,l") * tmps_["181_abba_Loovo"]("L,l,i,g,o");
    tmps_["262_bb_Lov"]("L,i,a") += t1_1["aa"]("b,l") * tmps_["260_aabb_Lovvo"]("L,l,b,a,i");
    tmps_["262_bb_Lov"]("L,i,a") += 0.50 * eri["abba_oovo"]("k,i,a,o") * tmps_["242_aa_Loo"]("L,o,k");
    tmps_["262_bb_Lov"]("L,i,a") -= eri["bbbb_oovo"]("i,m,a,j") * tmps_["218_bb_Loo"]("L,j,m");
    tmps_["262_bb_Lov"]("L,i,a") += eri["abba_oovo"]("k,i,a,o") * tmps_["207_aa_Loo"]("L,o,k");
    tmps_["262_bb_Lov"]("L,i,a") += 0.50 * f["bb_ov"]("j,a") * tmps_["219_bb_Loo"]("L,i,j");
    tmps_["262_bb_Lov"]("L,i,a") -= 0.50 * reused_["12_bb_ov"]("j,a") * tmps_["220_bb_Loo"]("L,i,j");
    tmps_["262_bb_Lov"]("L,i,a") -= 0.50 * reused_["12_bb_ov"]("j,a") * tmps_["219_bb_Loo"]("L,i,j");
    tmps_["262_bb_Lov"]("L,i,a") += reused_["12_bb_ov"]("j,a") * tmps_["217_bb_Loo"]("L,i,j");
    tmps_["262_bb_Lov"]("L,i,a") += dp["bb_ov"]("j,a") * tmps_["255_bb_Loo"]("L,i,j");
    tmps_["262_bb_Lov"]("L,i,a") += 0.50 * f["bb_ov"]("j,a") * tmps_["220_bb_Loo"]("L,i,j");
    tmps_["262_bb_Lov"]("L,i,a") -= 0.50 * tmps_["258_bb_Lov"]("L,i,a");
    tmps_["262_bb_Lov"]("L,i,a") -= eri["bbbb_oovv"]("i,m,a,d") * tmps_["246_bb_Lvo"]("L,d,m");
    tmps_["262_bb_Lov"]("L,i,a") += dp["bb_ov"]("j,a") * tmps_["257_bb_Loo"]("L,i,j");
    tmps_["262_bb_Lov"]("L,i,a") += tmps_["259_bb_Lov"]("L,i,a");
    tmps_["262_bb_Lov"]("L,i,a") += eri["abab_oooo"]("k,i,l,j") * tmps_["181_abba_Loovo"]("L,l,j,a,k");
    tmps_["262_bb_Lov"]("L,i,a") += reused_["80_abab_oooo"]("k,i,l,j") * tmps_["181_abba_Loovo"]("L,l,j,a,k");
    tmps_["262_bb_Lov"]("L,i,a") += dp["bb_oo"]("i,j") * tmps_["182_bb_Lov"]("L,j,a");
    tmps_["259_bb_Lov"].~TArrayD();
    tmps_["258_bb_Lov"].~TArrayD();
    tmps_["257_bb_Loo"].~TArrayD();
    tmps_["256_bb_Lov"].~TArrayD();
    tmps_["255_bb_Loo"].~TArrayD();
    tmps_["254_bb_Lov"].~TArrayD();
    tmps_["253_bb_Lvo"].~TArrayD();
    tmps_["252_bb_Lov"].~TArrayD();
    tmps_["251_bb_Lov"].~TArrayD();
    tmps_["250_bb_Loo"].~TArrayD();
    tmps_["246_bb_Lvo"].~TArrayD();
    tmps_["244_L"].~TArrayD();
    tmps_["243_bb_Loo"].~TArrayD();
    tmps_["242_aa_Loo"].~TArrayD();
    tmps_["241_bb_Lvv"].~TArrayD();
    tmps_["240_aa_Lvv"].~TArrayD();
    tmps_["221_bbbb_Lovvo"].~TArrayD();
    tmps_["60_bbbb_Loovo"].~TArrayD();

    // sigmal1_bb += +0.50 <l,i||j,k>_abab t1_1_aa(b,l) l2_1_abab(j,k,b,a)
    //            += +0.50 <l,i||k,j>_abab t1_1_aa(b,l) l2_1_abab(k,j,b,a)
    //            += +1.00 d-_bb(k,a) t1_1_bb(b,j) t1_1_bb(c,k) l2_1_bbbb(i,j,c,b)
    //            += +1.00 d-_bb(k,a) t1_1_aa(b,j) t1_1_bb(c,k) l2_1_abab(j,i,b,c)
    //            += +0.50 <i,l||d,a>_bbbb t2_abab(c,d,j,k) t1_1_bb(b,l) l2_1_abab(j,k,c,b)
    //            += +0.50 <i,l||d,a>_bbbb t2_abab(c,d,k,j) t1_1_bb(b,l) l2_1_abab(k,j,c,b)
    //            += -0.50 <i,l||d,a>_bbbb t2_bbbb(d,c,j,k) t1_1_bb(b,l) l2_1_bbbb(j,k,c,b)
    //            += -0.50 d+_bb(k,a) t2_bbbb(c,b,j,k) l2_1_bbbb(i,j,c,b)
    //            += +0.50 <l,k||d,a>_abab t2_bbbb(c,b,j,k) t1_1_aa(d,l) l2_1_bbbb(i,j,c,b)
    //            += +1.00 d-_bb(i,c) t1_1_bb(b,j) t1_1_bb(c,k) l2_1_bbbb(j,k,b,a)
    //            += +1.00 d-_bb(j,a) t1_1_bb(b,j) l1_bb(i,b)
    //            += +1.00 <l,k||d,a>_bbbb t2_abab(c,d,j,k) t1_1_bb(b,l) l2_1_abab(j,i,c,b)
    //            += +0.50 <i,l||d,a>_bbbb t2_abab(c,b,k,l) t1_1_bb(d,j) l2_1_abab(k,j,c,b)
    //            += +0.50 <i,l||d,a>_bbbb t2_abab(b,c,k,l) t1_1_bb(d,j) l2_1_abab(k,j,b,c)
    //            += -0.50 d-_bb(i,c) t0_1 t2_1_bbbb(c,b,j,k) l2_1_bbbb(j,k,b,a)
    //            += +0.50 d-_bb(i,c) t0_1 t2_1_abab(b,c,j,k) l2_1_abab(j,k,b,a)
    //            += +0.50 d-_bb(i,c) t0_1 t2_1_abab(b,c,k,j) l2_1_abab(k,j,b,a)
    //            += +1.00 d-_bb(i,b) t0_1 t1_1_bb(b,j) l1_1_bb(j,a)
    //            += -0.50 <l,k||d,a>_bbbb t2_abab(c,b,j,k) t1_1_bb(d,l) l2_1_abab(j,i,c,b)
    //            += -0.50 <l,k||d,a>_bbbb t2_abab(b,c,j,k) t1_1_bb(d,l) l2_1_abab(j,i,b,c)
    //            += -1.00 <l,k||c,d>_bbbb t2_abab(b,c,j,k) t1_1_bb(d,l) l2_1_abab(j,i,b,a)
    //            += -0.50 <i,l||d,a>_bbbb t2_bbbb(c,b,k,l) t1_1_bb(d,j) l2_1_bbbb(j,k,c,b)
    //            += +0.50 <l,k||d,a>_bbbb t2_bbbb(c,b,j,k) t1_1_bb(d,l) l2_1_bbbb(i,j,c,b)
    //            += +1.00 <l,k||d,a>_bbbb t2_bbbb(d,c,j,k) t1_1_bb(b,l) l2_1_bbbb(i,j,c,b)
    //            += -1.00 f_bb(j,a) t1_1_bb(b,j) l1_1_bb(i,b)
    //            += +0.50 <c,i||d,a>_abab t2_aaaa(d,b,j,k) l2_aaaa(j,k,c,b)
    //            += +0.50 <c,i||d,a>_abab t2_1_aaaa(d,b,j,k) l2_1_aaaa(j,k,c,b)
    //            += -0.50 <i,l||a,k>_bbbb t2_1_abab(c,b,j,l) l2_1_abab(j,k,c,b)
    //            += -0.50 <i,l||a,k>_bbbb t2_1_abab(b,c,j,l) l2_1_abab(j,k,b,c)
    //            += +1.00 <l,k||c,d>_aaaa t2_aaaa(c,b,j,k) t1_1_aa(d,l) l2_1_abab(j,i,b,a)
    //            += -0.50 <i,c||d,a>_bbbb t2_1_abab(b,d,j,k) l2_1_abab(j,k,b,c)
    //            += -0.50 <i,c||d,a>_bbbb t2_1_abab(b,d,k,j) l2_1_abab(k,j,b,c)
    //            += -0.50 <i,c||d,a>_bbbb t2_abab(b,d,j,k) l2_abab(j,k,b,c)
    //            += -0.50 <i,c||d,a>_bbbb t2_abab(b,d,k,j) l2_abab(k,j,b,c)
    //            += -0.50 <i,l||a,k>_bbbb t2_1_bbbb(c,b,j,l) l2_1_bbbb(j,k,c,b)
    //            += -0.50 <i,l||a,k>_bbbb t2_bbbb(c,b,j,l) l2_bbbb(j,k,c,b)
    //            += +1.00 <k,l||d,a>_abab t2_abab(d,c,k,j) t1_1_bb(b,l) l2_1_bbbb(i,j,c,b)
    //            += -1.00 <c,k||j,a>_abab t1_1_bb(b,k) l2_1_abab(j,i,c,b)
    //            += +1.00 <k,l||d,a>_abab t2_aaaa(d,c,j,k) t1_1_bb(b,l) l2_1_abab(j,i,c,b)
    //            += -0.50 <l,k||d,a>_abab t2_abab(c,b,j,k) t1_1_aa(d,l) l2_1_abab(j,i,c,b)
    //            += -0.50 <l,k||d,a>_abab t2_abab(b,c,j,k) t1_1_aa(d,l) l2_1_abab(j,i,b,c)
    //            += -0.50 <i,c||d,a>_bbbb t2_1_bbbb(d,b,j,k) l2_1_bbbb(j,k,c,b)
    //            += -0.50 <i,c||d,a>_bbbb t2_bbbb(d,b,j,k) l2_bbbb(j,k,c,b)
    //            += -0.50 f_bb(k,a) t2_abab(c,b,j,k) l2_abab(j,i,c,b)
    //            += -0.50 f_bb(k,a) t2_abab(b,c,j,k) l2_abab(j,i,b,c)
    //            += +0.50 d+_bb(k,a) t2_abab(c,b,j,k) l2_1_abab(j,i,c,b)
    //            += +0.50 d+_bb(k,a) t2_abab(b,c,j,k) l2_1_abab(j,i,b,c)
    //            += +0.50 d-_bb(k,a) t2_abab(c,b,j,k) t0_1 l2_abab(j,i,c,b)
    //            += +0.50 d-_bb(k,a) t2_abab(b,c,j,k) t0_1 l2_abab(j,i,b,c)
    //            += +1.00 d-_bb(j,a) t0_1 t1_1_bb(b,j) l1_1_bb(i,b)
    //            += -1.00 <i,k||a,j>_bbbb t1_1_bb(b,k) l1_1_bb(j,b)
    //            += -1.00 d-_bb(k,c) t2_bbbb(c,b,j,k) t0_1 l2_bbbb(i,j,b,a)
    //            += +0.50 d-_bb(i,c) t2_abab(b,c,j,k) t0_1 l2_abab(j,k,b,a)
    //            += +0.50 d-_bb(i,c) t2_abab(b,c,k,j) t0_1 l2_abab(k,j,b,a)
    //            += +1.00 d-_aa(k,c) t2_aaaa(c,b,j,k) t0_1 l2_abab(j,i,b,a)
    //            += +1.00 d-_aa(k,c) t2_abab(c,b,k,j) t0_1 l2_bbbb(i,j,b,a)
    //            += -0.50 d-_bb(i,c) t2_bbbb(c,b,j,k) t0_1 l2_bbbb(j,k,b,a)
    //            += -1.00 <k,i||j,a>_abab t1_1_aa(b,k) l1_1_aa(j,b)
    //            += -0.50 <l,i||d,a>_abab t2_abab(c,b,l,k) t1_1_aa(d,j) l2_1_abab(j,k,c,b)
    //            += -0.50 <l,i||d,a>_abab t2_abab(b,c,l,k) t1_1_aa(d,j) l2_1_abab(j,k,b,c)
    //            += -0.50 f_bb(k,a) t2_1_abab(c,b,j,k) l2_1_abab(j,i,c,b)
    //            += -0.50 f_bb(k,a) t2_1_abab(b,c,j,k) l2_1_abab(j,i,b,c)
    //            += -1.00 d-_bb(c,a) t1_1_aa(b,j) l2_abab(j,i,b,c)
    //            += +1.00 <l,k||d,c>_abab t2_abab(b,c,j,k) t1_1_aa(d,l) l2_1_abab(j,i,b,a)
    //            += -0.50 <l,i||k,a>_abab t2_abab(c,b,l,j) l2_abab(k,j,c,b)
    //            += -0.50 <l,i||k,a>_abab t2_abab(b,c,l,j) l2_abab(k,j,b,c)
    //            += -1.00 d+_bb(i,a) l0_1
    //            += -0.50 <k,j||c,a>_abab t2_abab(c,b,k,j) l1_bb(i,b)
    //            += -0.50 <j,k||c,a>_abab t2_abab(c,b,j,k) l1_bb(i,b)
    //            += +0.25 <l,k||a,j>_bbbb t2_bbbb(c,b,l,k) l2_bbbb(i,j,c,b)
    //            += -1.00 <k,c||d,a>_abab t2_1_aaaa(d,b,j,k) l2_1_abab(j,i,b,c)
    //            += -1.00 d-_aa(j,j) t0_1 l1_bb(i,a)
    //                += +1.00 f_aa(k,k) r2_1_bbbb(a,b,i,j)
    //                += -0.50 <l,k||l,k>_abab r2_1_bbbb(a,b,i,j)
    //                += -0.50 <k,l||k,l>_abab r2_1_bbbb(a,b,i,j)
    //                += +0.25 <l,k||c,d>_aaaa r2_1_bbbb(a,b,i,j) t2_aaaa(c,d,l,k)
    //                += -0.50 <l,k||l,k>_bbbb r2_1_bbbb(a,b,i,j)
    //                += -0.50 <l,k||l,k>_aaaa r2_1_bbbb(a,b,i,j)
    //                += +0.25 <l,k||c,d>_bbbb r2_1_bbbb(a,b,i,j) t2_bbbb(c,d,l,k)
    //                += +1.00 f_bb(k,k) r2_1_bbbb(a,b,i,j)
    //                += -1.00 d-_bb(k,k) r2_1_bbbb(a,b,i,j) t0_1
    //                += +0.25 <l,k||c,d>_abab r2_1_bbbb(a,b,i,j) t2_abab(c,d,l,k)
    //                += +0.25 <k,l||c,d>_abab r2_1_bbbb(a,b,i,j) t2_abab(c,d,k,l)
    //                += +0.25 <l,k||d,c>_abab r2_1_bbbb(a,b,i,j) t2_abab(d,c,l,k)
    //                += +0.25 <k,l||d,c>_abab r2_1_bbbb(a,b,i,j) t2_abab(d,c,k,l)
    //            += -0.50 f_bb(i,c) t2_abab(b,c,j,k) l2_abab(j,k,b,a)
    //            += -0.50 f_bb(i,c) t2_abab(b,c,k,j) l2_abab(k,j,b,a)
    //            += -0.25 <b,i||c,d>_abab t2_abab(c,d,j,k) l2_abab(j,k,b,a)
    //            += -0.25 <b,i||c,d>_abab t2_abab(c,d,k,j) l2_abab(k,j,b,a)
    //            += -0.25 <b,i||d,c>_abab t2_abab(d,c,j,k) l2_abab(j,k,b,a)
    //            += -0.25 <b,i||d,c>_abab t2_abab(d,c,k,j) l2_abab(k,j,b,a)
    //            += +1.00 <k,i||c,a>_abab t2_1_abab(c,b,k,j) l1_1_bb(j,b)
    //            += +1.00 <i,b||a,j>_bbbb l1_bb(j,b)
    //            += -1.00 d-_bb(j,b) t1_1_bb(b,j) l1_bb(i,a)
    //            += -1.00 d+_bb(b,a) l1_1_bb(i,b)
    //            += +0.25 <l,k||j,a>_abab t2_abab(c,b,l,k) l2_abab(j,i,c,b)
    //            += +0.25 <k,l||j,a>_abab t2_abab(c,b,k,l) l2_abab(j,i,c,b)
    //            += +0.25 <l,k||j,a>_abab t2_abab(b,c,l,k) l2_abab(j,i,b,c)
    //            += +0.25 <k,l||j,a>_abab t2_abab(b,c,k,l) l2_abab(j,i,b,c)
    //            += -0.50 d-_bb(i,c) t2_1_bbbb(c,b,j,k) l2_bbbb(j,k,b,a)
    //            += -1.00 <l,i||c,k>_abab t2_abab(c,b,l,j) l2_bbbb(j,k,b,a)
    //            += +0.25 <l,k||j,a>_abab t2_1_abab(c,b,l,k) l2_1_abab(j,i,c,b)
    //            += +0.25 <k,l||j,a>_abab t2_1_abab(c,b,k,l) l2_1_abab(j,i,c,b)
    //            += +0.25 <l,k||j,a>_abab t2_1_abab(b,c,l,k) l2_1_abab(j,i,b,c)
    //            += +0.25 <k,l||j,a>_abab t2_1_abab(b,c,k,l) l2_1_abab(j,i,b,c)
    //            += +1.00 d-_bb(i,j) t0_1 l1_bb(j,a)
    //            += +1.00 d-_bb(i,k) t1_1_bb(b,j) l2_bbbb(j,k,b,a)
    //            += +1.00 <i,k||c,a>_bbbb t2_1_bbbb(c,b,j,k) l1_1_bb(j,b)
    //            += -1.00 <i,b||c,a>_bbbb t1_1_bb(c,j) l1_1_bb(j,b)
    //            += -1.00 d+_bb(k,c) t2_abab(b,c,j,k) l2_1_abab(j,i,b,a)
    //            += -1.00 <k,i||c,a>_abab t2_1_aaaa(c,b,j,k) l1_1_aa(j,b)
    //            += +0.50 f_bb(i,c) t2_1_bbbb(c,b,j,k) l2_1_bbbb(j,k,b,a)
    //            += +0.25 <i,b||c,d>_bbbb t2_1_bbbb(c,d,j,k) l2_1_bbbb(j,k,b,a)
    //            += +0.50 <l,i||d,c>_abab t2_bbbb(c,b,j,k) t1_1_aa(d,l) l2_1_bbbb(j,k,b,a)
    //            += -0.50 <k,b||c,d>_aaaa t2_aaaa(c,d,j,k) l2_abab(j,i,b,a)
    //            += -0.50 <l,k||j,c>_abab t2_abab(b,c,l,k) l2_abab(j,i,b,a)
    //            += -0.50 <k,l||j,c>_abab t2_abab(b,c,k,l) l2_abab(j,i,b,a)
    //            += -0.50 <l,k||c,j>_aaaa t2_aaaa(c,b,l,k) l2_abab(j,i,b,a)
    //            += -1.00 f_aa(k,c) t2_aaaa(c,b,j,k) l2_abab(j,i,b,a)
    //            += +0.50 <b,k||c,d>_abab t2_abab(c,d,j,k) l2_abab(j,i,b,a)
    //            += +0.50 <b,k||d,c>_abab t2_abab(d,c,j,k) l2_abab(j,i,b,a)
    //            += -1.00 d-_aa(b,j) t0_1 l2_abab(j,i,b,a)
    //            += +1.00 f_bb(k,c) t2_abab(b,c,j,k) l2_abab(j,i,b,a)
    //            += -0.50 <l,k||d,c>_abab t2_abab(b,c,l,k) t1_1_aa(d,j) l2_1_abab(j,i,b,a)
    //            += -0.50 <k,l||d,c>_abab t2_abab(b,c,k,l) t1_1_aa(d,j) l2_1_abab(j,i,b,a)
    //            += +1.00 d-_aa(k,j) t0_1 t1_1_aa(b,k) l2_1_abab(j,i,b,a)
    //            += +1.00 d-_aa(k,c) t0_1 t2_1_aaaa(c,b,j,k) l2_1_abab(j,i,b,a)
    //            += -1.00 d-_aa(b,c) t0_1 t1_1_aa(c,j) l2_1_abab(j,i,b,a)
    //            += -1.00 d-_bb(k,c) t0_1 t2_1_abab(b,c,j,k) l2_1_abab(j,i,b,a)
    //            += -1.00 <k,l||c,d>_abab t2_aaaa(c,b,j,k) t1_1_bb(d,l) l2_1_abab(j,i,b,a)
    //            += +2.00 d-_aa(k,c) t1_1_aa(b,k) t1_1_aa(c,j) l2_1_abab(j,i,b,a)
    //            += -0.50 <l,k||c,d>_abab t2_abab(c,d,j,k) t1_1_aa(b,l) l2_1_abab(j,i,b,a)
    //            += -0.50 <l,k||d,c>_abab t2_abab(d,c,j,k) t1_1_aa(b,l) l2_1_abab(j,i,b,a)
    //            += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,b,l,k) t1_1_aa(d,j) l2_1_abab(j,i,b,a)
    //            += +1.00 t1_1_aa(b,j) l2_1_abab(j,i,b,a) w0
    //            += -1.00 f_aa(k,c) t2_1_aaaa(c,b,j,k) l2_1_abab(j,i,b,a)
    //            += -1.00 d-_bb(k,c) t1_1_aa(b,j) t1_1_bb(c,k) l2_1_abab(j,i,b,a)
    //            += -0.50 <l,k||j,c>_abab t2_1_abab(b,c,l,k) l2_1_abab(j,i,b,a)
    //            += -0.50 <k,l||j,c>_abab t2_1_abab(b,c,k,l) l2_1_abab(j,i,b,a)
    //            += -1.00 d-_aa(k,c) t1_1_aa(b,j) t1_1_aa(c,k) l2_1_abab(j,i,b,a)
    //            += +1.00 <k,b||c,j>_aaaa t1_1_aa(c,k) l2_1_abab(j,i,b,a)
    //            += +1.00 f_bb(k,c) t2_1_abab(b,c,j,k) l2_1_abab(j,i,b,a)
    //            += -1.00 f_aa(k,j) t1_1_aa(b,k) l2_1_abab(j,i,b,a)
    //            += -0.50 <k,b||c,d>_aaaa t2_1_aaaa(c,d,j,k) l2_1_abab(j,i,b,a)
    //            += +1.00 <b,k||j,c>_abab t1_1_bb(c,k) l2_1_abab(j,i,b,a)
    //            += -0.50 <l,k||c,j>_aaaa t2_1_aaaa(c,b,l,k) l2_1_abab(j,i,b,a)
    //            += +0.50 <b,k||c,d>_abab t2_1_abab(c,d,j,k) l2_1_abab(j,i,b,a)
    //            += +0.50 <b,k||d,c>_abab t2_1_abab(d,c,j,k) l2_1_abab(j,i,b,a)
    //            += +1.00 f_aa(b,c) t1_1_aa(c,j) l2_1_abab(j,i,b,a)
    //            += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,d,j,k) t1_1_aa(b,l) l2_1_abab(j,i,b,a)
    //            += +0.50 <l,k||c,d>_bbbb t2_bbbb(c,d,j,k) t1_1_bb(b,l) l2_1_bbbb(i,j,b,a)
    //            += +1.00 <l,k||d,c>_abab t2_bbbb(c,b,j,k) t1_1_aa(d,l) l2_1_bbbb(i,j,b,a)
    //            += +0.50 <l,k||c,d>_abab t2_abab(c,b,l,k) t1_1_bb(d,j) l2_1_bbbb(i,j,b,a)
    //            += +0.50 <k,l||c,d>_abab t2_abab(c,b,k,l) t1_1_bb(d,j) l2_1_bbbb(i,j,b,a)
    //            += -1.00 f_aa(k,c) t2_1_abab(c,b,k,j) l2_1_bbbb(i,j,b,a)
    //            += +1.00 f_bb(k,j) t1_1_bb(b,k) l2_1_bbbb(i,j,b,a)
    //            += +0.50 <k,b||c,d>_bbbb t2_1_bbbb(c,d,j,k) l2_1_bbbb(i,j,b,a)
    //            += -1.00 <k,b||c,j>_abab t1_1_aa(c,k) l2_1_bbbb(i,j,b,a)
    //            += +1.00 f_bb(k,c) t2_1_bbbb(c,b,j,k) l2_1_bbbb(i,j,b,a)
    //            += -1.00 f_bb(b,c) t1_1_bb(c,j) l2_1_bbbb(i,j,b,a)
    //            += -0.50 <k,b||c,d>_abab t2_1_abab(c,d,k,j) l2_1_bbbb(i,j,b,a)
    //            += -0.50 <k,b||d,c>_abab t2_1_abab(d,c,k,j) l2_1_bbbb(i,j,b,a)
    //            += -1.00 <k,b||c,j>_bbbb t1_1_bb(c,k) l2_1_bbbb(i,j,b,a)
    //            += +0.50 <l,k||c,j>_bbbb t2_1_bbbb(c,b,l,k) l2_1_bbbb(i,j,b,a)
    //            += +1.00 d-_bb(k,c) t1_1_bb(b,j) t1_1_bb(c,k) l2_1_bbbb(i,j,b,a)
    //            += +1.00 d-_aa(k,c) t1_1_bb(b,j) t1_1_aa(c,k) l2_1_bbbb(i,j,b,a)
    //            += -1.00 t1_1_bb(b,j) l2_1_bbbb(i,j,b,a) w0
    //            += +0.50 <l,k||c,j>_abab t2_1_abab(c,b,l,k) l2_1_bbbb(i,j,b,a)
    //            += +0.50 <k,l||c,j>_abab t2_1_abab(c,b,k,l) l2_1_bbbb(i,j,b,a)
    //            += +0.50 <k,l||c,d>_abab t2_abab(c,d,k,j) t1_1_bb(b,l) l2_1_bbbb(i,j,b,a)
    //            += +0.50 <k,l||d,c>_abab t2_abab(d,c,k,j) t1_1_bb(b,l) l2_1_bbbb(i,j,b,a)
    //            += -2.00 d-_bb(k,c) t1_1_bb(b,k) t1_1_bb(c,j) l2_1_bbbb(i,j,b,a)
    //            += +1.00 d-_bb(b,c) t0_1 t1_1_bb(c,j) l2_1_bbbb(i,j,b,a)
    //            += -1.00 d-_bb(k,j) t0_1 t1_1_bb(b,k) l2_1_bbbb(i,j,b,a)
    //            += +1.00 d-_aa(k,c) t0_1 t2_1_abab(c,b,k,j) l2_1_bbbb(i,j,b,a)
    //            += -1.00 d-_bb(k,c) t0_1 t2_1_bbbb(c,b,j,k) l2_1_bbbb(i,j,b,a)
    //            += +0.50 <l,k||c,d>_bbbb t2_bbbb(c,b,l,k) t1_1_bb(d,j) l2_1_bbbb(i,j,b,a)
    //            += -1.00 <k,l||c,d>_abab t2_abab(c,b,k,j) t1_1_bb(d,l) l2_1_bbbb(i,j,b,a)
    //            += +1.00 <j,b||c,a>_bbbb t1_1_bb(c,j) l1_1_bb(i,b)
    //            += -0.50 <k,j||c,a>_bbbb t2_1_bbbb(c,b,k,j) l1_1_bb(i,b)
    //            += +1.00 <l,i||c,d>_abab t2_abab(c,b,l,k) t1_1_bb(d,j) l2_1_bbbb(j,k,b,a)
    //            += -1.00 <l,i||c,k>_abab t2_1_abab(c,b,l,j) l2_1_bbbb(j,k,b,a)
    //            += -1.00 <i,k||c,a>_bbbb t2_1_abab(b,c,j,k) l1_1_aa(j,b)
    //            += +1.00 <i,b||c,k>_bbbb t1_1_bb(c,j) l2_1_bbbb(j,k,b,a)
    //            += +1.00 <j,b||c,a>_abab t1_1_aa(c,j) l1_1_bb(i,b)
    //            += -0.50 <k,j||c,a>_abab t2_1_abab(c,b,k,j) l1_1_bb(i,b)
    //            += -0.50 <j,k||c,a>_abab t2_1_abab(c,b,j,k) l1_1_bb(i,b)
    //            += -1.00 <k,i||b,j>_abab t1_1_aa(b,k) l1_1_bb(j,a)
    //            += -0.50 <k,i||b,c>_abab t2_1_abab(b,c,k,j) l1_1_bb(j,a)
    //            += -0.50 <k,i||c,b>_abab t2_1_abab(c,b,k,j) l1_1_bb(j,a)
    //            += -1.00 d-_bb(k,c) t2_abab(b,c,j,k) t0_1 l2_abab(j,i,b,a)
    //            += -1.00 f_bb(i,j) l1_bb(j,a)
    //            += +1.00 <i,l||c,k>_bbbb t2_abab(b,c,j,l) l2_abab(j,k,b,a)
    //            += +1.00 <l,i||k,c>_abab t2_1_abab(b,c,l,j) l2_1_abab(k,j,b,a)
    //            += +1.00 <l,i||c,d>_abab t2_aaaa(c,b,k,l) t1_1_bb(d,j) l2_1_abab(k,j,b,a)
    //            += +0.25 <l,k||a,j>_bbbb t2_1_bbbb(c,b,l,k) l2_1_bbbb(i,j,c,b)
    //            += -1.00 <k,i||c,a>_abab t2_aaaa(c,b,j,k) l1_aa(j,b)
    //            += +1.00 <k,c||d,a>_bbbb t2_abab(b,d,j,k) l2_abab(j,i,b,c)
    //            += -1.00 <k,c||d,a>_bbbb t2_bbbb(d,b,j,k) l2_bbbb(i,j,c,b)
    //            += +1.00 <b,i||c,a>_abab t1_1_aa(c,j) l1_1_aa(j,b)
    //            += -0.50 <k,j||c,a>_bbbb t2_bbbb(c,b,k,j) l1_bb(i,b)
    //            += +1.00 d-_bb(i,b) t1_1_bb(b,j) l1_bb(j,a)
    //            += -0.50 <l,i||d,c>_abab t2_abab(b,c,j,k) t1_1_aa(d,l) l2_1_abab(j,k,b,a)
    //            += -0.50 <l,i||d,c>_abab t2_abab(b,c,k,j) t1_1_aa(d,l) l2_1_abab(k,j,b,a)
    //            += +1.00 <l,i||c,k>_abab t2_1_aaaa(c,b,j,l) l2_1_abab(j,k,b,a)
    //            += +1.00 <l,i||d,c>_abab t2_abab(b,c,l,k) t1_1_aa(d,j) l2_1_abab(j,k,b,a)
    //            += +0.50 d-_bb(i,c) t2_1_abab(b,c,j,k) l2_abab(j,k,b,a)
    //            += +0.50 d-_bb(i,c) t2_1_abab(b,c,k,j) l2_abab(k,j,b,a)
    //            += -0.50 <i,k||b,c>_bbbb t2_bbbb(b,c,j,k) l1_bb(j,a)
    //            += -1.00 <b,i||k,c>_abab t1_1_bb(c,j) l2_1_abab(k,j,b,a)
    //            += -1.00 d+_bb(j,j) l1_1_bb(i,a)
    //            += -1.00 d-_aa(j,b) t1_1_aa(b,j) l1_bb(i,a)
    //            += +1.00 f_aa(b,j) l2_abab(j,i,b,a)
    //            += -1.00 <i,j||b,a>_bbbb t1_1_bb(b,j) l0_1
    //            += -1.00 <c,k||d,a>_abab t2_1_abab(d,b,j,k) l2_1_abab(j,i,c,b)
    //            += +0.50 <c,b||d,a>_abab t1_1_aa(d,j) l2_1_abab(j,i,c,b)
    //            += +0.50 <b,c||d,a>_abab t1_1_aa(d,j) l2_1_abab(j,i,b,c)
    //            += +0.25 <l,k||d,a>_abab t2_abab(c,b,l,k) t1_1_aa(d,j) l2_1_abab(j,i,c,b)
    //            += +0.25 <k,l||d,a>_abab t2_abab(c,b,k,l) t1_1_aa(d,j) l2_1_abab(j,i,c,b)
    //            += +0.25 <l,k||d,a>_abab t2_abab(b,c,l,k) t1_1_aa(d,j) l2_1_abab(j,i,b,c)
    //            += +0.25 <k,l||d,a>_abab t2_abab(b,c,k,l) t1_1_aa(d,j) l2_1_abab(j,i,b,c)
    //            += +0.50 <i,l||c,d>_bbbb t2_bbbb(c,b,j,k) t1_1_bb(d,l) l2_1_bbbb(j,k,b,a)
    //            += +1.00 <i,l||c,d>_bbbb t2_bbbb(c,b,k,l) t1_1_bb(d,j) l2_1_bbbb(j,k,b,a)
    //            += -1.00 <i,l||c,k>_bbbb t2_1_bbbb(c,b,j,l) l2_1_bbbb(j,k,b,a)
    //            += -0.50 <i,l||j,k>_bbbb t1_1_bb(b,l) l2_1_bbbb(j,k,b,a)
    //            += -0.25 <i,l||c,d>_bbbb t2_bbbb(c,d,j,k) t1_1_bb(b,l) l2_1_bbbb(j,k,b,a)
    //            += -1.00 <b,i||c,k>_abab t1_1_aa(c,j) l2_1_abab(j,k,b,a)
    //            += -0.50 f_bb(i,c) t2_1_abab(b,c,j,k) l2_1_abab(j,k,b,a)
    //            += -0.50 f_bb(i,c) t2_1_abab(b,c,k,j) l2_1_abab(k,j,b,a)
    //            += -0.25 <b,i||c,d>_abab t2_1_abab(c,d,j,k) l2_1_abab(j,k,b,a)
    //            += -0.25 <b,i||c,d>_abab t2_1_abab(c,d,k,j) l2_1_abab(k,j,b,a)
    //            += -0.25 <b,i||d,c>_abab t2_1_abab(d,c,j,k) l2_1_abab(j,k,b,a)
    //            += -0.25 <b,i||d,c>_abab t2_1_abab(d,c,k,j) l2_1_abab(k,j,b,a)
    //            += +1.00 <k,c||d,a>_abab t2_abab(d,b,k,j) l2_bbbb(i,j,c,b)
    //            += +1.00 d-_aa(k,k) t1_1_bb(b,j) l2_bbbb(i,j,b,a)
    //            += +1.00 d-_bb(k,k) t1_1_bb(b,j) l2_bbbb(i,j,b,a)
    //            += +0.50 <c,b||a,j>_bbbb l2_bbbb(i,j,c,b)
    //            += +1.00 f_bb(k,c) t2_bbbb(c,b,j,k) l2_bbbb(i,j,b,a)
    //            += +1.00 d-_bb(b,j) t0_1 l2_bbbb(i,j,b,a)
    //            += +0.50 <l,k||c,j>_abab t2_abab(c,b,l,k) l2_bbbb(i,j,b,a)
    //            += +0.50 <k,l||c,j>_abab t2_abab(c,b,k,l) l2_bbbb(i,j,b,a)
    //            += +0.50 <l,k||c,j>_bbbb t2_bbbb(c,b,l,k) l2_bbbb(i,j,b,a)
    //            += -1.00 f_aa(k,c) t2_abab(c,b,k,j) l2_bbbb(i,j,b,a)
    //            += -0.50 <k,b||c,d>_abab t2_abab(c,d,k,j) l2_bbbb(i,j,b,a)
    //            += -0.50 <k,b||d,c>_abab t2_abab(d,c,k,j) l2_bbbb(i,j,b,a)
    //            += +0.50 <k,b||c,d>_bbbb t2_bbbb(c,d,j,k) l2_bbbb(i,j,b,a)
    //            += +1.00 f_bb(b,a) l1_bb(i,b)
    //            += -1.00 d+_aa(j,j) l1_1_bb(i,a)
    //            += -0.50 <i,l||c,d>_bbbb t2_abab(b,c,j,k) t1_1_bb(d,l) l2_1_abab(j,k,b,a)
    //            += -0.50 <i,l||c,d>_bbbb t2_abab(b,c,k,j) t1_1_bb(d,l) l2_1_abab(k,j,b,a)
    //            += +1.00 <i,l||c,k>_bbbb t2_1_abab(b,c,j,l) l2_1_abab(j,k,b,a)
    //            += -1.00 d-_bb(k,k) t1_1_aa(b,j) l2_abab(j,i,b,a)
    //            += -1.00 d-_aa(k,k) t1_1_aa(b,j) l2_abab(j,i,b,a)
    //            += +1.00 <k,c||d,a>_abab t2_1_abab(d,b,k,j) l2_1_bbbb(i,j,c,b)
    //            += +1.00 <k,c||d,a>_bbbb t2_1_abab(b,d,j,k) l2_1_abab(j,i,b,c)
    //            += +1.00 d-_aa(k,j) t1_1_aa(b,k) l2_abab(j,i,b,a)
    //            += +1.00 d-_aa(k,c) t2_1_aaaa(c,b,j,k) l2_abab(j,i,b,a)
    //            += -1.00 d-_aa(b,c) t1_1_aa(c,j) l2_abab(j,i,b,a)
    //            += -1.00 d-_bb(k,c) t2_1_abab(b,c,j,k) l2_abab(j,i,b,a)
    //            += -1.00 f_bb(i,b) t1_1_bb(b,j) l1_1_bb(j,a)
    //            += +0.50 <c,b||j,a>_abab l2_abab(j,i,c,b)
    //            += +0.50 <b,c||j,a>_abab l2_abab(j,i,b,c)
    //            += +0.50 <i,b||j,k>_bbbb l2_bbbb(j,k,b,a)
    //            += +1.00 <k,i||c,a>_abab t2_abab(c,b,k,j) l1_bb(j,b)
    //            += -0.50 <b,i||j,k>_abab l2_abab(j,k,b,a)
    //            += -0.50 <b,i||k,j>_abab l2_abab(k,j,b,a)
    //            += +0.50 f_bb(i,c) t2_bbbb(c,b,j,k) l2_bbbb(j,k,b,a)
    //            += +0.25 <i,b||c,d>_bbbb t2_bbbb(c,d,j,k) l2_bbbb(j,k,b,a)
    //            += +1.00 <i,k||c,a>_bbbb t2_bbbb(c,b,j,k) l1_bb(j,b)
    //            += +1.00 t0_1 l1_1_bb(i,a) w0
    //            += +0.25 <k,j||b,c>_aaaa t2_1_aaaa(b,c,k,j) l1_1_bb(i,a)
    //            += +0.25 <k,j||b,c>_abab t2_1_abab(b,c,k,j) l1_1_bb(i,a)
    //            += +0.25 <j,k||b,c>_abab t2_1_abab(b,c,j,k) l1_1_bb(i,a)
    //            += +0.25 <k,j||c,b>_abab t2_1_abab(c,b,k,j) l1_1_bb(i,a)
    //            += +0.25 <j,k||c,b>_abab t2_1_abab(c,b,j,k) l1_1_bb(i,a)
    //            += -1.00 d-_aa(j,b) t0_1 t1_1_aa(b,j) l1_1_bb(i,a)
    //            += +1.00 f_bb(j,b) t1_1_bb(b,j) l1_1_bb(i,a)
    //            += -1.00 d-_bb(j,b) t0_1 t1_1_bb(b,j) l1_1_bb(i,a)
    //            += +0.25 <k,j||b,c>_bbbb t2_1_bbbb(b,c,k,j) l1_1_bb(i,a)
    //            += +1.00 f_aa(j,b) t1_1_aa(b,j) l1_1_bb(i,a)
    //            += -0.50 <k,i||b,c>_abab t2_abab(b,c,k,j) l1_bb(j,a)
    //            += -0.50 <k,i||c,b>_abab t2_abab(c,b,k,j) l1_bb(j,a)
    //            += +1.00 d+_bb(i,j) l1_1_bb(j,a)
    //            += -1.00 <c,k||d,a>_abab t2_abab(d,b,j,k) l2_abab(j,i,c,b)
    //            += -1.00 f_bb(b,j) l2_bbbb(i,j,b,a)
    //            += -1.00 <k,c||d,a>_abab t2_aaaa(d,b,j,k) l2_abab(j,i,b,c)
    //            += -1.00 <i,l||c,k>_bbbb t2_bbbb(c,b,j,l) l2_bbbb(j,k,b,a)
    //            += -1.00 d+_aa(b,j) l2_1_abab(j,i,b,a)
    //            += -1.00 d-_bb(c,a) t1_1_bb(b,j) l2_bbbb(i,j,c,b)
    //            += +1.00 <i,k||b,j>_bbbb t1_1_bb(b,k) l1_1_bb(j,a)
    //            += -0.50 <i,k||b,c>_bbbb t2_1_bbbb(b,c,j,k) l1_1_bb(j,a)
    //            += +1.00 <l,i||c,k>_abab t2_aaaa(c,b,j,l) l2_abab(j,k,b,a)
    //            += +1.00 <b,i||j,a>_abab l1_aa(j,b)
    //            += -1.00 d-_bb(b,a) t0_1 l1_bb(i,b)
    //            += +1.00 <i,l||c,d>_bbbb t2_abab(b,c,k,l) t1_1_bb(d,j) l2_1_abab(k,j,b,a)
    //            += -1.00 <i,k||c,a>_bbbb t2_abab(b,c,j,k) l1_aa(j,b)
    //            += -0.25 <l,k||d,a>_bbbb t2_bbbb(c,b,l,k) t1_1_bb(d,j) l2_1_bbbb(i,j,c,b)
    //            += -1.00 <k,c||d,a>_bbbb t2_1_bbbb(d,b,j,k) l2_1_bbbb(i,j,c,b)
    //            += -0.50 <c,b||d,a>_bbbb t1_1_bb(d,j) l2_1_bbbb(i,j,c,b)
    //            += +1.00 d-_bb(b,c) t1_1_bb(c,j) l2_bbbb(i,j,b,a)
    //            += -1.00 d-_bb(k,j) t1_1_bb(b,k) l2_bbbb(i,j,b,a)
    //            += +1.00 d-_aa(k,c) t2_1_abab(c,b,k,j) l2_bbbb(i,j,b,a)
    //            += -1.00 d-_bb(k,c) t2_1_bbbb(c,b,j,k) l2_bbbb(i,j,b,a)
    //            += +1.00 <l,i||k,c>_abab t2_abab(b,c,l,j) l2_abab(k,j,b,a)
    //            += +1.00 <j,i||b,a>_abab t1_1_aa(b,j) l0_1
    //            += +1.00 <l,k||d,a>_abab t2_abab(d,c,j,k) t1_1_aa(b,l) l2_1_abab(j,i,b,c)
    //            += +1.00 <k,c||a,j>_bbbb t1_1_bb(b,k) l2_1_bbbb(i,j,c,b)
    //            += +0.50 <c,i||d,a>_abab t2_abab(d,b,j,k) l2_abab(j,k,c,b)
    //            += +0.50 <c,i||d,a>_abab t2_abab(d,b,k,j) l2_abab(k,j,c,b)
    //            += +0.50 <c,i||d,a>_abab t2_1_abab(d,b,j,k) l2_1_abab(j,k,c,b)
    //            += +0.50 <c,i||d,a>_abab t2_1_abab(d,b,k,j) l2_1_abab(k,j,c,b)
    //            += -1.00 d-_bb(i,a) t1_1_bb(b,j) l1_bb(j,b)
    //            += -1.00 d-_bb(i,a) t1_1_aa(b,j) l1_aa(j,b)
    //            += -0.25 d-_bb(i,a) t2_1_bbbb(c,b,j,k) l2_bbbb(j,k,c,b)
    //            += -0.25 d-_bb(i,a) t2_1_aaaa(c,b,j,k) l2_aaaa(j,k,c,b)
    //            += -0.25 d-_bb(i,a) t2_1_abab(c,b,j,k) l2_abab(j,k,c,b)
    //            += -0.25 d-_bb(i,a) t2_1_abab(c,b,k,j) l2_abab(k,j,c,b)
    //            += -0.25 d-_bb(i,a) t2_1_abab(b,c,j,k) l2_abab(j,k,b,c)
    //            += -0.25 d-_bb(i,a) t2_1_abab(b,c,k,j) l2_abab(k,j,b,c)
    //            += -1.00 <l,k||c,d>_bbbb t2_bbbb(c,b,j,k) t1_1_bb(d,l) l2_1_bbbb(i,j,b,a)
    //            += +1.00 <l,k||c,d>_aaaa t2_abab(c,b,k,j) t1_1_aa(d,l) l2_1_bbbb(i,j,b,a)
    //            += -1.00 <k,c||j,a>_abab t1_1_aa(b,k) l2_1_abab(j,i,b,c)
    //            += +1.00 d-_bb(i,c) t1_1_aa(b,j) t1_1_bb(c,k) l2_1_abab(j,k,b,a)
    //            += -0.50 <l,i||k,a>_abab t2_aaaa(c,b,j,l) l2_aaaa(j,k,c,b)
    //            += -0.50 <l,i||k,a>_abab t2_1_aaaa(c,b,j,l) l2_1_aaaa(j,k,c,b)
    //            += -0.50 <i,l||a,k>_bbbb t2_abab(c,b,j,l) l2_abab(j,k,c,b)
    //            += -0.50 <i,l||a,k>_bbbb t2_abab(b,c,j,l) l2_abab(j,k,b,c)
    //            += -0.50 <l,i||k,a>_abab t2_1_abab(c,b,l,j) l2_1_abab(k,j,c,b)
    //            += -0.50 <l,i||k,a>_abab t2_1_abab(b,c,l,j) l2_1_abab(k,j,b,c)
    //            += +0.50 f_bb(k,a) t2_1_bbbb(c,b,j,k) l2_1_bbbb(i,j,c,b)
    //            += -0.50 d-_bb(k,a) t2_bbbb(c,b,j,k) t0_1 l2_bbbb(i,j,c,b)
    //            += -0.50 d-_bb(k,a) t0_1 t2_1_bbbb(c,b,j,k) l2_1_bbbb(i,j,c,b)
    //            += +0.50 d-_bb(k,a) t0_1 t2_1_abab(c,b,j,k) l2_1_abab(j,i,c,b)
    //            += +0.50 d-_bb(k,a) t0_1 t2_1_abab(b,c,j,k) l2_1_abab(j,i,b,c)
    //            += +0.50 d-_bb(k,a) t2_1_abab(c,b,j,k) l2_abab(j,i,c,b)
    //            += +0.50 d-_bb(k,a) t2_1_abab(b,c,j,k) l2_abab(j,i,b,c)
    //            += -0.50 d-_bb(k,a) t2_1_bbbb(c,b,j,k) l2_bbbb(i,j,c,b)
    //            += +0.50 f_bb(k,a) t2_bbbb(c,b,j,k) l2_bbbb(i,j,c,b)
    //            += +0.50 <l,i||d,a>_abab t2_aaaa(c,b,k,l) t1_1_aa(d,j) l2_1_aaaa(j,k,c,b)
    //            += -0.50 <l,i||d,a>_abab t2_abab(d,c,j,k) t1_1_aa(b,l) l2_1_abab(j,k,b,c)
    //            += -0.50 <l,i||d,a>_abab t2_abab(d,c,k,j) t1_1_aa(b,l) l2_1_abab(k,j,b,c)
    //            += +0.50 <l,i||d,a>_abab t2_aaaa(d,c,j,k) t1_1_aa(b,l) l2_1_aaaa(j,k,c,b)
    //            += +0.25 <l,i||c,d>_abab t2_abab(c,d,j,k) t1_1_aa(b,l) l2_1_abab(j,k,b,a)
    //            += +0.25 <l,i||c,d>_abab t2_abab(c,d,k,j) t1_1_aa(b,l) l2_1_abab(k,j,b,a)
    //            += +0.25 <l,i||d,c>_abab t2_abab(d,c,j,k) t1_1_aa(b,l) l2_1_abab(j,k,b,a)
    //            += +0.25 <l,i||d,c>_abab t2_abab(d,c,k,j) t1_1_aa(b,l) l2_1_abab(k,j,b,a)
    //            += +1.00 d-_bb(i,k) t1_1_aa(b,j) l2_abab(j,k,b,a)
    sigmal1_bb("L,a,i") += tmps_["262_bb_Lov"]("L,i,a");
    tmps_["262_bb_Lov"].~TArrayD();

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["263_bb_Lov"]("L,i,a")  = l2_1["bbbb"]("L,i,j,b,a") * t1_1["bb"]("b,j");

    // csigmal1_1_bb += -1.00 t1_1_bb(b,j) l2_1_bbbb(i,j,b,a)
    csigmal1_1_bb("L,a,i") -= tmps_["263_bb_Lov"]("L,i,a");

    // flops: o2v2L1  = o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["264_bbbb_Lovov"]("L,i,a,j,b")  = -1.00 * dp["bb_ov"]("i,a") * tmps_["263_bb_Lov"]("L,j,b");
    tmps_["264_bbbb_Lovov"]("L,i,a,j,b") -= eri["bbbb_vovo"]("d,i,a,l") * l2_1["bbbb"]("L,j,l,d,b");
    tmps_["264_bbbb_Lovov"]("L,i,a,j,b") -= f["bb_ov"]("i,a") * l1_1["bb"]("L,j,b");
    tmps_["264_bbbb_Lovov"]("L,i,a,j,b") += l2_1["bbbb"]("L,j,l,d,b") * reused_["26_bbbb_voov"]("d,l,i,a");
    tmps_["264_bbbb_Lovov"]("L,i,a,j,b") += eri["abba_vovo"]("c,i,a,k") * l2_1["abab"]("L,k,j,c,b");
    tmps_["264_bbbb_Lovov"]("L,i,a,j,b") += dp["bb_ov"]("i,a") * l1["bb"]("L,j,b");
    tmps_["264_bbbb_Lovov"]("L,i,a,j,b") += l1_1["bb"]("L,j,b") * reused_["12_bb_ov"]("i,a");
    tmps_["264_bbbb_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,k,j,c,b") * reused_["6_aabb_voov"]("c,k,i,a");
    tmps_["264_bbbb_Lovov"]("L,i,a,j,b") -= l2_1["bbbb"]("L,j,l,d,b") * reused_["163_bbbb_voov"]("d,l,i,a");
    tmps_["264_bbbb_Lovov"]("L,i,a,j,b") -= l2_1["abab"]("L,k,j,c,b") * reused_["70_aabb_voov"]("c,k,i,a");
    tmps_["264_bbbb_Lovov"]("L,i,a,j,b") += dp["bb_ov"]("i,a") * tmps_["101_bb_Lov"]("L,j,b");

    // sigmal2_1_bbbb += +1.00 P(i,j) P(a,b) d-_bb(j,a) t1_1_aa(c,k) l2_1_abab(k,i,c,b)
    //                += +1.00 P(i,j) P(a,b) <l,j||d,a>_abab t2_abab(d,c,l,k) l2_1_bbbb(i,k,c,b)
    //                += -1.00 P(i,j) P(a,b) f_bb(j,a) l1_1_bb(i,b)
    //                += -1.00 P(i,j) P(a,b) <c,j||k,a>_abab l2_1_abab(k,i,c,b)
    //                += +1.00 P(i,j) P(a,b) d-_bb(j,a) l1_bb(i,b)
    //                += +1.00 P(i,j) P(a,b) d-_bb(j,a) t0_1 l1_1_bb(i,b)
    //                += +1.00 P(i,j) P(a,b) <l,j||d,a>_abab t2_aaaa(d,c,k,l) l2_1_abab(k,i,c,b)
    //                += +1.00 P(i,j) P(a,b) <j,c||a,k>_bbbb l2_1_bbbb(i,k,c,b)
    //                += +1.00 P(i,j) P(a,b) <j,l||d,a>_bbbb t2_bbbb(d,c,k,l) l2_1_bbbb(i,k,c,b)
    //                += +1.00 P(i,j) P(a,b) <j,l||d,a>_bbbb t2_abab(c,d,k,l) l2_1_abab(k,i,c,b)
    //                += -1.00 P(i,j) P(a,b) d-_bb(j,a) t1_1_bb(c,k) l2_1_bbbb(i,k,c,b)
    sigmal2_1_bbbb("L,a,b,i,j") += tmps_["264_bbbb_Lovov"]("L,j,a,i,b");
    sigmal2_1_bbbb("L,a,b,i,j") -= tmps_["264_bbbb_Lovov"]("L,i,a,j,b");
    sigmal2_1_bbbb("L,a,b,i,j") -= tmps_["264_bbbb_Lovov"]("L,j,b,i,a");
    sigmal2_1_bbbb("L,a,b,i,j") += tmps_["264_bbbb_Lovov"]("L,i,b,j,a");
    tmps_["264_bbbb_Lovov"].~TArrayD();

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["265_bbbb_Lovov"]("L,i,a,j,b")  = -1.00 * dp["bb_ov"]("i,a") * tmps_["263_bb_Lov"]("L,j,b");
    tmps_["265_bbbb_Lovov"].~TArrayD();

    // flops: o2v2L1  = o0v2L1 o2v3L1 o3v2L1 o0v2 o2v3L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o0v2L1 o2v2L1 o2v2L1 o0v2 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["266_bbbb_Lvoov"]("L,a,i,j,b")  = (tmps_["238_bb_Lvv"]("L,a,c") + -2.00 * tmps_["199_bb_Lvv"]("L,a,c")) * eri["bbbb_oovv"]("i,j,b,c");
    tmps_["266_bbbb_Lvoov"]("L,a,i,j,b") += 4.00 * dp["bb_ov"]("k,b") * tmps_["113_bbbb_Loovo"]("L,i,j,a,k");
    tmps_["266_bbbb_Lvoov"]("L,a,i,j,b") -= 2.00 * (reused_["34_bb_vv"]("d,b") + 0.50 * reused_["109_bb_vv"]("b,d")) * l2_1["bbbb"]("L,i,j,d,a");
    tmps_["266_bbbb_Lvoov"]("L,a,i,j,b") -= 2.00 * dp["bb_vv"]("d,b") * l2["bbbb"]("L,i,j,d,a");
    tmps_["266_bbbb_Lvoov"]("L,a,i,j,b") -= 2.00 * eri["bbbb_oovo"]("i,j,b,k") * l1_1["bb"]("L,k,a");
    tmps_["266_bbbb_Lvoov"]("L,a,i,j,b") -= 2.00 * reused_["3_bb_vv"]("b,d") * l2_1["bbbb"]("L,i,j,d,a");
    tmps_["266_bbbb_Lvoov"]("L,a,i,j,b") += 2.00 * f["bb_vv"]("d,b") * l2_1["bbbb"]("L,i,j,d,a");

    // sigmal2_1_bbbb += -0.50 P(a,b) <i,j||d,a>_bbbb t2_bbbb(d,c,k,l) l2_1_bbbb(k,l,c,b)
    //                += +0.50 P(a,b) <i,j||d,a>_bbbb t2_abab(c,d,k,l) l2_1_abab(k,l,c,b)
    //                += +0.50 P(a,b) <i,j||d,a>_bbbb t2_abab(c,d,l,k) l2_1_abab(l,k,c,b)
    //                += +2.00 P(a,b) d-_bb(k,a) t1_1_bb(c,k) l2_1_bbbb(i,j,c,b)
    //                += -1.00 P(a,b) d-_bb(c,a) t0_1 l2_1_bbbb(i,j,c,b)
    //                += -0.50 P(a,b) <l,k||d,a>_bbbb t2_bbbb(d,c,l,k) l2_1_bbbb(i,j,c,b)
    //                += -1.00 P(a,b) d-_bb(c,a) l2_bbbb(i,j,c,b)
    //                += -1.00 P(a,b) <i,j||a,k>_bbbb l1_1_bb(k,b)
    //                += -0.50 P(a,b) <l,k||d,a>_abab t2_abab(d,c,l,k) l2_1_bbbb(i,j,c,b)
    //                += -0.50 P(a,b) <k,l||d,a>_abab t2_abab(d,c,k,l) l2_1_bbbb(i,j,c,b)
    //                += +1.00 P(a,b) f_bb(c,a) l2_1_bbbb(i,j,c,b)
    sigmal2_1_bbbb("L,a,b,i,j") += 0.50 * tmps_["266_bbbb_Lvoov"]("L,b,i,j,a");
    sigmal2_1_bbbb("L,a,b,i,j") -= 0.50 * tmps_["266_bbbb_Lvoov"]("L,a,i,j,b");
    tmps_["266_bbbb_Lvoov"].~TArrayD();

    // flops: o1v1L1  = o3v2L1 o2v2L1 o3v2L1 o1v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    tmps_["270_bb_Lvo"]("L,a,i")  = -2.00 * l2_1["abab"]("L,l,k,c,a") * reused_["77_aabb_vooo"]("c,l,k,i");
    tmps_["270_bb_Lvo"]("L,a,i") -= 2.00 * l2_1["abab"]("L,l,i,c,a") * reused_["54_aa_vo"]("c,l");
    tmps_["270_bb_Lvo"]("L,a,i") += l2_1["bbbb"]("L,j,k,b,a") * reused_["51_bbbb_vooo"]("b,j,k,i");

    // sigmal1_bb += -0.50 d+_bb(i,c) t2_bbbb(c,b,j,k) l2_1_bbbb(j,k,b,a)
    //            += +1.00 d+_aa(k,c) t2_aaaa(c,b,j,k) l2_1_abab(j,i,b,a)
    //            += +0.50 d+_bb(i,c) t2_abab(b,c,j,k) l2_1_abab(j,k,b,a)
    //            += +0.50 d+_bb(i,c) t2_abab(b,c,k,j) l2_1_abab(k,j,b,a)
    sigmal1_bb("L,a,i") -= 0.50 * tmps_["270_bb_Lvo"]("L,a,i");

    // flops: o1v1L1  = o2v1L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["271_bb_Lvo"]("L,a,i")  = -0.50 * f["bb_ov"]("j,a") * tmps_["179_bb_Loo"]("L,i,j");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["269_bb_Lov"]("L,i,a")  = l2_1["bbbb"]("L,j,i,b,a") * t1_1["bb"]("b,j");

    // flops: o1v1L1  = o3v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["267_bb_Lvo"]("L,a,i")  = -2.00 * l2_1["abab"]("L,j,k,b,a") * reused_["77_aabb_vooo"]("b,j,k,i");
    tmps_["267_bb_Lvo"].~TArrayD();

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["268_bb_Lov"]("L,i,a")  = l2_1["bbbb"]("L,i,j,a,b") * t1_1["bb"]("b,j");

    // flops: o1v1L1  = o2v1L1 o2v1L1 o1v3L1 o2v1L1 o3v1L1 o2v1L1 o1v3L1 o1v2L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v3L1 o1v1L1 o3v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o3v1L1 o1v1L1 o1v2L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o1v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o1v1L1 o1v1L1 o1v2L1 o1v1L1 o2v3L1 o1v1L1 o2v1L1 o1v1L1 o3v2L1 o1v1L1 o2v3L1 o1v1L1 o1v1L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o2v1L1 o1v1L1 o3v2L1 o1v1L1 o1v1L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o1v1L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o2v3L1 o1v1L1 o1v1L1 o1v3L1 o1v1L1 o3v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    tmps_["272_bb_Lvo"]("L,a,i")  = -1.00 * dp["bb_ov"]("j,a") * tmps_["218_bb_Loo"]("L,i,j");
    tmps_["272_bb_Lvo"]("L,a,i") -= reused_["12_bb_ov"]("j,a") * tmps_["178_bb_Loo"]("L,i,j");
    tmps_["272_bb_Lvo"]("L,a,i") -= eri["abab_vovv"]("f,i,g,a") * tmps_["99_aa_Lvv"]("L,f,g");
    tmps_["272_bb_Lvo"]("L,a,i") -= 2.00 * dp["bb_ov"]("n,a") * tmps_["222_bb_Loo"]("L,i,n");
    tmps_["272_bb_Lvo"]("L,a,i") += 0.50 * eri["bbbb_oovo"]("i,o,a,j") * tmps_["208_bb_Loo"]("L,j,o");
    tmps_["272_bb_Lvo"]("L,a,i") += dp["bb_ov"]("j,a") * tmps_["219_bb_Loo"]("L,i,j");
    tmps_["272_bb_Lvo"]("L,a,i") -= 0.50 * eri["abab_vovv"]("f,i,g,a") * tmps_["198_aa_Lvv"]("L,f,g");
    tmps_["272_bb_Lvo"]("L,a,i") += dp["bb_vv"]("c,a") * tmps_["101_bb_Lov"]("L,i,c");
    tmps_["272_bb_Lvo"]("L,a,i") += t0_1 * tmps_["104_bb_Lov"]("L,i,a");
    tmps_["272_bb_Lvo"]("L,a,i") += dp["bb_ov"]("i,a") * tmps_["210_L"]("L");
    tmps_["272_bb_Lvo"]("L,a,i") += eri["bbbb_vovv"]("c,i,a,e") * tmps_["199_bb_Lvv"]("L,c,e");
    tmps_["272_bb_Lvo"]("L,a,i") += eri["bbbb_oovo"]("i,o,a,j") * tmps_["178_bb_Loo"]("L,j,o");
    tmps_["272_bb_Lvo"]("L,a,i") += f["bb_ov"]("j,a") * tmps_["178_bb_Loo"]("L,i,j");
    tmps_["272_bb_Lvo"]("L,a,i") += 0.25 * dp["bb_ov"]("i,a") * tmps_["103_L"]("L");
    tmps_["272_bb_Lvo"]("L,a,i") -= 2.00 * dp["bb_ov"]("j,a") * tmps_["217_bb_Loo"]("L,i,j");
    tmps_["272_bb_Lvo"]("L,a,i") += 0.50 * t0_1 * tmps_["270_bb_Lvo"]("L,a,i");
    tmps_["272_bb_Lvo"]("L,a,i") -= eri["abba_oovo"]("l,i,a,m") * tmps_["205_aa_Loo"]("L,m,l");
    tmps_["272_bb_Lvo"]("L,a,i") += dp["bb_vv"]("c,a") * tmps_["268_bb_Lov"]("L,i,c");
    tmps_["272_bb_Lvo"]("L,a,i") += 0.50 * dp["bb_ov"]("j,a") * tmps_["220_bb_Loo"]("L,i,j");
    tmps_["272_bb_Lvo"]("L,a,i") += 0.50 * reused_["12_bb_ov"]("j,a") * tmps_["179_bb_Loo"]("L,i,j");
    tmps_["272_bb_Lvo"]("L,a,i") -= 2.00 * l2_1["bbbb"]("L,i,n,d,a") * reused_["88_bb_vo"]("d,n");
    tmps_["272_bb_Lvo"]("L,a,i") -= l2_1["abab"]("L,k,j,b,a") * reused_["217_bbaa_oovo"]("i,j,b,k");
    tmps_["272_bb_Lvo"]("L,a,i") -= l1_1["aa"]("L,k,b") * reused_["70_aabb_voov"]("b,k,i,a");
    tmps_["272_bb_Lvo"]("L,a,i") += l2_1["abab"]("L,k,i,b,a") * reused_["152_aa_vo"]("b,k");
    tmps_["272_bb_Lvo"]("L,a,i") += f["bb_vo"]("d,n") * l2_1["bbbb"]("L,i,n,d,a");
    tmps_["272_bb_Lvo"]("L,a,i") += scalars_["2"] * l1["bb"]("L,i,a");
    tmps_["272_bb_Lvo"]("L,a,i") -= f["bb_vv"]("d,a") * l1_1["bb"]("L,i,d");
    tmps_["272_bb_Lvo"]("L,a,i") += eri["bbbb_vovo"]("d,i,a,n") * l1_1["bb"]("L,n,d");
    tmps_["272_bb_Lvo"]("L,a,i") -= f["aa_vo"]("b,k") * l2_1["abab"]("L,k,i,b,a");
    tmps_["272_bb_Lvo"]("L,a,i") -= 0.50 * l2_1["bbbb"]("L,n,j,d,a") * reused_["53_bbbb_vooo"]("d,n,j,i");
    tmps_["272_bb_Lvo"]("L,a,i") -= w0 * l1_1["bb"]("L,i,a");
    tmps_["272_bb_Lvo"]("L,a,i") += reused_["34_bb_vv"]("d,a") * l1_1["bb"]("L,i,d");
    tmps_["272_bb_Lvo"]("L,a,i") -= reused_["210_bbaa_vvvo"]("c,a,b,k") * l2_1["abab"]("L,k,i,b,c");
    tmps_["272_bb_Lvo"]("L,a,i") += f["bb_oo"]("i,n") * l1_1["bb"]("L,n,a");
    tmps_["272_bb_Lvo"]("L,a,i") += l2_1["bbbb"]("L,n,j,d,a") * reused_["220_bbbb_vooo"]("d,n,i,j");
    tmps_["272_bb_Lvo"]("L,a,i") -= reused_["213_aabb_vovv"]("b,k,c,a") * l2_1["abab"]("L,k,i,b,c");
    tmps_["272_bb_Lvo"]("L,a,i") += l0_1("L") * reused_["12_bb_ov"]("i,a");
    tmps_["272_bb_Lvo"]("L,a,i") += 0.25 * reused_["41_bbbb_vovv"]("a,n,c,d") * l2_1["bbbb"]("L,i,n,c,d");
    tmps_["272_bb_Lvo"]("L,a,i") -= l2_1["bbbb"]("L,i,n,d,a") * reused_["91_bb_vo"]("d,n");
    tmps_["272_bb_Lvo"]("L,a,i") += reused_["209_abba_vvvo"]("f,a,d,k") * l2_1["abab"]("L,k,i,f,d");
    tmps_["272_bb_Lvo"]("L,a,i") += l1_1["bb"]("L,n,a") * reused_["30_bb_oo"]("i,n");
    tmps_["272_bb_Lvo"]("L,a,i") += eri["abab_vooo"]("b,i,k,j") * l2_1["abab"]("L,k,j,b,a");
    tmps_["272_bb_Lvo"]("L,a,i") += scalars_["1"] * l1["bb"]("L,i,a");
    tmps_["272_bb_Lvo"]("L,a,i") += reused_["176_bbbb_vvvo"]("c,a,d,n") * l2_1["bbbb"]("L,i,n,c,d");
    tmps_["272_bb_Lvo"]("L,a,i") += l2_1["abab"]("L,k,j,b,a") * reused_["225_aabb_vooo"]("b,k,j,i");
    tmps_["272_bb_Lvo"]("L,a,i") += l2_1["abab"]("L,m,n,b,a") * reused_["219_abba_vooo"]("b,n,i,m");
    tmps_["272_bb_Lvo"]("L,a,i") += l2["abab"]("L,k,i,b,a") * reused_["58_aa_vo"]("b,k");
    tmps_["272_bb_Lvo"]("L,a,i") += 2.00 * scalars_["4"] * l1_1["bb"]("L,i,a");
    tmps_["272_bb_Lvo"]("L,a,i") += 0.50 * l1_1["bb"]("L,n,a") * reused_["161_bb_oo"]("i,n");
    tmps_["272_bb_Lvo"]("L,a,i") += l1_1["aa"]("L,k,b") * reused_["6_aabb_voov"]("b,k,i,a");
    tmps_["272_bb_Lvo"]("L,a,i") -= l2_1["bbbb"]("L,i,n,d,a") * reused_["84_bb_vo"]("d,n");
    tmps_["272_bb_Lvo"]("L,a,i") += l2_1["abab"]("L,k,i,b,a") * reused_["90_aa_vo"]("b,k");
    tmps_["272_bb_Lvo"]("L,a,i") += dp["bb_vv"]("d,a") * l1["bb"]("L,i,d");
    tmps_["272_bb_Lvo"]("L,a,i") += reused_["215_baab_vovv"]("a,k,f,d") * l2_1["abab"]("L,k,i,f,d");
    tmps_["272_bb_Lvo"]("L,a,i") -= l2_1["abab"]("L,k,j,b,a") * reused_["218_aabb_vooo"]("b,k,i,j");
    tmps_["272_bb_Lvo"]("L,a,i") += scalars_["5"] * l1_1["bb"]("L,i,a");
    tmps_["272_bb_Lvo"]("L,a,i") += reused_["214_bbbb_vovv"]("d,n,c,a") * l2_1["bbbb"]("L,i,n,c,d");
    tmps_["272_bb_Lvo"]("L,a,i") += l2_1["bbbb"]("L,n,j,d,a") * reused_["178_bbbb_oovo"]("i,j,d,n");
    tmps_["272_bb_Lvo"]("L,a,i") -= l1_1["bb"]("L,n,a") * reused_["35_bb_oo"]("i,n");
    tmps_["272_bb_Lvo"]("L,a,i") -= 0.50 * l2_1["abab"]("L,k,i,b,a") * reused_["87_aa_vo"]("b,k");
    tmps_["272_bb_Lvo"]("L,a,i") -= f["bb_ov"]("i,a") * l0_1("L");
    tmps_["272_bb_Lvo"]("L,a,i") += dp["aa_vo"]("b,k") * l2["abab"]("L,k,i,b,a");
    tmps_["272_bb_Lvo"]("L,a,i") += l1_1["bb"]("L,n,d") * reused_["163_bbbb_voov"]("d,n,i,a");
    tmps_["272_bb_Lvo"]("L,a,i") -= tmps_["269_bb_Lov"]("L,j,a") * dp["bb_oo"]("i,j");
    tmps_["272_bb_Lvo"]("L,a,i") += eri["abba_vvvo"]("f,d,a,k") * l2_1["abab"]("L,k,i,f,d");
    tmps_["272_bb_Lvo"]("L,a,i") += eri["abba_vovo"]("b,i,a,k") * l1_1["aa"]("L,k,b");
    tmps_["272_bb_Lvo"]("L,a,i") -= dp["bb_vo"]("d,n") * l2["bbbb"]("L,i,n,d,a");
    tmps_["272_bb_Lvo"]("L,a,i") += 0.50 * reused_["109_bb_vv"]("a,d") * l1_1["bb"]("L,i,d");
    tmps_["272_bb_Lvo"]("L,a,i") += 0.50 * eri["bbbb_vooo"]("d,i,n,j") * l2_1["bbbb"]("L,n,j,d,a");
    tmps_["272_bb_Lvo"]("L,a,i") -= 2.00 * l2_1["abab"]("L,k,i,b,a") * reused_["89_aa_ov"]("k,b");
    tmps_["272_bb_Lvo"]("L,a,i") += reused_["3_bb_vv"]("a,d") * l1_1["bb"]("L,i,d");
    tmps_["272_bb_Lvo"]("L,a,i") += 2.00 * scalars_["3"] * l1_1["bb"]("L,i,a");
    tmps_["272_bb_Lvo"]("L,a,i") -= l1_1["bb"]("L,n,d") * reused_["26_bbbb_voov"]("d,n,i,a");
    tmps_["272_bb_Lvo"]("L,a,i") -= dp["bb_oo"]("i,n") * l1["bb"]("L,n,a");
    tmps_["272_bb_Lvo"]("L,a,i") += 0.50 * eri["bbbb_vvvo"]("d,c,a,n") * l2_1["bbbb"]("L,i,n,c,d");
    tmps_["272_bb_Lvo"]("L,a,i") += 0.50 * eri["bbbb_vovv"]("c,i,a,e") * tmps_["200_bb_Lvv"]("L,c,e");
    tmps_["272_bb_Lvo"]("L,a,i") -= 0.50 * eri["abba_oovo"]("l,i,a,m") * tmps_["204_aa_Loo"]("L,m,l");
    tmps_["272_bb_Lvo"]("L,a,i") += tmps_["271_bb_Lvo"]("L,a,i");
    tmps_["272_bb_Lvo"]("L,a,i") -= dp["bb_oo"]("i,j") * tmps_["101_bb_Lov"]("L,j,a");
    tmps_["271_bb_Lvo"].~TArrayD();
    tmps_["270_bb_Lvo"].~TArrayD();
    tmps_["269_bb_Lov"].~TArrayD();
    tmps_["268_bb_Lov"].~TArrayD();
    tmps_["210_L"].~TArrayD();
    tmps_["208_bb_Loo"].~TArrayD();
    tmps_["204_aa_Loo"].~TArrayD();
    tmps_["200_bb_Lvv"].~TArrayD();
    tmps_["198_aa_Lvv"].~TArrayD();
    tmps_["104_bb_Lov"].~TArrayD();
    tmps_["103_L"].~TArrayD();

    // sigmal1_1_bb += -1.00 d-_bb(c,a) t1_1_aa(b,j) l2_1_abab(j,i,b,c)
    //              += +0.50 <c,i||d,a>_abab t2_aaaa(d,b,j,k) l2_1_aaaa(j,k,c,b)
    //              += -1.00 d-_bb(k,c) t2_bbbb(c,b,j,k) t0_1 l2_1_bbbb(i,j,b,a)
    //              += +1.00 d-_aa(k,c) t2_abab(c,b,k,j) t0_1 l2_1_bbbb(i,j,b,a)
    //              += -1.00 d-_bb(k,a) t2_1_bbbb(c,b,j,k) l2_1_bbbb(i,j,c,b)
    //              += -0.50 <i,l||a,k>_bbbb t2_bbbb(c,b,j,l) l2_1_bbbb(j,k,c,b)
    //              += +2.00 d-_bb(j,a) t1_1_bb(b,j) l1_1_bb(i,b)
    //              += +0.50 <c,i||d,a>_abab t2_abab(d,b,j,k) l2_1_abab(j,k,c,b)
    //              += +0.50 <c,i||d,a>_abab t2_abab(d,b,k,j) l2_1_abab(k,j,c,b)
    //              += -1.00 d-_bb(i,a) t1_1_bb(b,j) l1_1_bb(j,b)
    //              += -1.00 d-_bb(i,a) t1_1_aa(b,j) l1_1_aa(j,b)
    //              += -0.50 <i,c||d,a>_bbbb t2_abab(b,d,j,k) l2_1_abab(j,k,b,c)
    //              += -0.50 <i,c||d,a>_bbbb t2_abab(b,d,k,j) l2_1_abab(k,j,b,c)
    //              += -0.50 <i,l||a,k>_bbbb t2_abab(c,b,j,l) l2_1_abab(j,k,c,b)
    //              += -0.50 <i,l||a,k>_bbbb t2_abab(b,c,j,l) l2_1_abab(j,k,b,c)
    //              += +0.50 d-_bb(k,a) t2_abab(c,b,j,k) t0_1 l2_1_abab(j,i,c,b)
    //              += +0.50 d-_bb(k,a) t2_abab(b,c,j,k) t0_1 l2_1_abab(j,i,b,c)
    //              += -0.50 f_bb(k,a) t2_abab(c,b,j,k) l2_1_abab(j,i,c,b)
    //              += -0.50 f_bb(k,a) t2_abab(b,c,j,k) l2_1_abab(j,i,b,c)
    //              += -0.25 d-_bb(i,a) t2_1_bbbb(c,b,j,k) l2_1_bbbb(j,k,c,b)
    //              += -0.25 d-_bb(i,a) t2_1_aaaa(c,b,j,k) l2_1_aaaa(j,k,c,b)
    //              += -0.25 d-_bb(i,a) t2_1_abab(c,b,j,k) l2_1_abab(j,k,c,b)
    //              += -0.25 d-_bb(i,a) t2_1_abab(c,b,k,j) l2_1_abab(k,j,c,b)
    //              += -0.25 d-_bb(i,a) t2_1_abab(b,c,j,k) l2_1_abab(j,k,b,c)
    //              += -0.25 d-_bb(i,a) t2_1_abab(b,c,k,j) l2_1_abab(k,j,b,c)
    //              += +1.00 d-_bb(k,a) t2_1_abab(c,b,j,k) l2_1_abab(j,i,c,b)
    //              += +1.00 d-_bb(k,a) t2_1_abab(b,c,j,k) l2_1_abab(j,i,b,c)
    //              += -0.50 d-_bb(i,c) t2_bbbb(c,b,j,k) t0_1 l2_1_bbbb(j,k,b,a)
    //              += +1.00 d-_aa(k,c) t2_aaaa(c,b,j,k) t0_1 l2_1_abab(j,i,b,a)
    //              += +0.50 d-_bb(i,c) t2_abab(b,c,j,k) t0_1 l2_1_abab(j,k,b,a)
    //              += +0.50 d-_bb(i,c) t2_abab(b,c,k,j) t0_1 l2_1_abab(k,j,b,a)
    //              += +0.50 d-_bb(k,a) t2_abab(c,b,j,k) l2_abab(j,i,c,b)
    //              += +0.50 d-_bb(k,a) t2_abab(b,c,j,k) l2_abab(j,i,b,c)
    //              += -0.50 <l,i||k,a>_abab t2_abab(c,b,l,j) l2_1_abab(k,j,c,b)
    //              += -0.50 <l,i||k,a>_abab t2_abab(b,c,l,j) l2_1_abab(k,j,b,c)
    //              += -1.00 d-_bb(c,a) t1_1_bb(b,j) l2_1_bbbb(i,j,c,b)
    //              += -0.50 d-_bb(k,a) t2_bbbb(c,b,j,k) l2_bbbb(i,j,c,b)
    //              += -0.50 d-_bb(k,a) t2_bbbb(c,b,j,k) t0_1 l2_1_bbbb(i,j,c,b)
    //              += +2.00 d-_bb(b,c) t1_1_bb(c,j) l2_1_bbbb(i,j,b,a)
    //              += -2.00 d-_bb(k,j) t1_1_bb(b,k) l2_1_bbbb(i,j,b,a)
    //              += +2.00 d-_aa(k,c) t2_1_abab(c,b,k,j) l2_1_bbbb(i,j,b,a)
    //              += -2.00 d-_bb(k,c) t2_1_bbbb(c,b,j,k) l2_1_bbbb(i,j,b,a)
    //              += +1.00 <l,i||c,k>_abab t2_aaaa(c,b,j,l) l2_1_abab(j,k,b,a)
    //              += -1.00 <i,k||c,a>_bbbb t2_abab(b,c,j,k) l1_1_aa(j,b)
    //              += -1.00 d-_bb(k,c) t2_abab(b,c,j,k) t0_1 l2_1_abab(j,i,b,a)
    //              += -1.00 f_bb(b,j) l2_1_bbbb(i,j,b,a)
    //              += -1.00 d-_bb(j,j) l1_bb(i,a)
    //              += +1.00 f_bb(b,a) l1_1_bb(i,b)
    //              += +1.00 <i,b||a,j>_bbbb l1_1_bb(j,b)
    //              += +1.00 f_aa(b,j) l2_1_abab(j,i,b,a)
    //              += +0.50 f_bb(i,c) t2_bbbb(c,b,j,k) l2_1_bbbb(j,k,b,a)
    //              += +0.25 <i,b||c,d>_bbbb t2_bbbb(c,d,j,k) l2_1_bbbb(j,k,b,a)
    //              += +1.00 l1_1_bb(i,a) w0
    //              += -1.00 d-_bb(b,a) t0_1 l1_1_bb(i,b)
    //              += -1.00 <k,c||d,a>_abab t2_aaaa(d,b,j,k) l2_1_abab(j,i,b,c)
    //              += -1.00 f_bb(i,j) l1_1_bb(j,a)
    //              += -1.00 <i,l||c,k>_bbbb t2_bbbb(c,b,j,l) l2_1_bbbb(j,k,b,a)
    //              += +1.00 <k,c||d,a>_bbbb t2_abab(b,d,j,k) l2_1_abab(j,i,b,c)
    //              += -1.00 d-_bb(i,a) t0_1 l0_1
    //              += +0.25 <l,k||a,j>_bbbb t2_bbbb(c,b,l,k) l2_1_bbbb(i,j,c,b)
    //              += +1.00 d-_aa(k,k) t1_1_bb(b,j) l2_1_bbbb(i,j,b,a)
    //              += +1.00 d-_bb(k,k) t1_1_bb(b,j) l2_1_bbbb(i,j,b,a)
    //              += -1.00 <c,k||d,a>_abab t2_abab(d,b,j,k) l2_1_abab(j,i,c,b)
    //              += -0.50 <k,i||b,c>_abab t2_abab(b,c,k,j) l1_1_bb(j,a)
    //              += -0.50 <k,i||c,b>_abab t2_abab(c,b,k,j) l1_1_bb(j,a)
    //              += -0.50 <b,i||j,k>_abab l2_1_abab(j,k,b,a)
    //              += -0.50 <b,i||k,j>_abab l2_1_abab(k,j,b,a)
    //              += -1.00 d-_aa(j,j) l1_bb(i,a)
    //              += +1.00 <k,c||d,a>_abab t2_abab(d,b,k,j) l2_1_bbbb(i,j,c,b)
    //              += -0.50 f_bb(i,c) t2_abab(b,c,j,k) l2_1_abab(j,k,b,a)
    //              += -0.50 f_bb(i,c) t2_abab(b,c,k,j) l2_1_abab(k,j,b,a)
    //              += -0.25 <b,i||c,d>_abab t2_abab(c,d,j,k) l2_1_abab(j,k,b,a)
    //              += -0.25 <b,i||c,d>_abab t2_abab(c,d,k,j) l2_1_abab(k,j,b,a)
    //              += -0.25 <b,i||d,c>_abab t2_abab(d,c,j,k) l2_1_abab(j,k,b,a)
    //              += -0.25 <b,i||d,c>_abab t2_abab(d,c,k,j) l2_1_abab(k,j,b,a)
    //              += +1.00 <l,i||k,c>_abab t2_abab(b,c,l,j) l2_1_abab(k,j,b,a)
    //              += -1.00 d-_bb(k,c) t2_abab(b,c,j,k) l2_abab(j,i,b,a)
    //              += -2.00 d-_aa(j,b) t1_1_aa(b,j) l1_1_bb(i,a)
    //              += -0.50 <i,k||b,c>_bbbb t2_bbbb(b,c,j,k) l1_1_bb(j,a)
    //              += -1.00 <k,i||c,a>_abab t2_aaaa(c,b,j,k) l1_1_aa(j,b)
    //              += +1.00 f_bb(k,c) t2_bbbb(c,b,j,k) l2_1_bbbb(i,j,b,a)
    //              += +1.00 d-_bb(b,j) t0_1 l2_1_bbbb(i,j,b,a)
    //              += +0.50 <l,k||c,j>_abab t2_abab(c,b,l,k) l2_1_bbbb(i,j,b,a)
    //              += +0.50 <k,l||c,j>_abab t2_abab(c,b,k,l) l2_1_bbbb(i,j,b,a)
    //              += +0.50 <l,k||c,j>_bbbb t2_bbbb(c,b,l,k) l2_1_bbbb(i,j,b,a)
    //              += -1.00 f_aa(k,c) t2_abab(c,b,k,j) l2_1_bbbb(i,j,b,a)
    //              += -0.50 <k,b||c,d>_abab t2_abab(c,d,k,j) l2_1_bbbb(i,j,b,a)
    //              += -0.50 <k,b||d,c>_abab t2_abab(d,c,k,j) l2_1_bbbb(i,j,b,a)
    //              += +0.50 <k,b||c,d>_bbbb t2_bbbb(c,d,j,k) l2_1_bbbb(i,j,b,a)
    //              += -1.00 d-_bb(k,k) t1_1_aa(b,j) l2_1_abab(j,i,b,a)
    //              += -1.00 d-_aa(k,k) t1_1_aa(b,j) l2_1_abab(j,i,b,a)
    //              += -1.00 d-_bb(b,a) l1_bb(i,b)
    //              += +0.25 <l,k||j,a>_abab t2_abab(c,b,l,k) l2_1_abab(j,i,c,b)
    //              += +0.25 <k,l||j,a>_abab t2_abab(c,b,k,l) l2_1_abab(j,i,c,b)
    //              += +0.25 <l,k||j,a>_abab t2_abab(b,c,l,k) l2_1_abab(j,i,b,c)
    //              += +0.25 <k,l||j,a>_abab t2_abab(b,c,k,l) l2_1_abab(j,i,b,c)
    //              += +1.00 <i,l||c,k>_bbbb t2_abab(b,c,j,l) l2_1_abab(j,k,b,a)
    //              += -1.00 d-_aa(j,j) t0_1 l1_1_bb(i,a)
    //              += +1.00 f_aa(k,k) r2_abab(a,b,i,j)
    //              += -0.50 <l,k||l,k>_abab r2_abab(a,b,i,j)
    //              += -0.50 <k,l||k,l>_abab r2_abab(a,b,i,j)
    //              += +0.25 <l,k||c,d>_aaaa r2_abab(a,b,i,j) t2_aaaa(c,d,l,k)
    //              += -0.50 <l,k||l,k>_bbbb r2_abab(a,b,i,j)
    //              += -0.50 <l,k||l,k>_aaaa r2_abab(a,b,i,j)
    //              += +0.25 <l,k||c,d>_bbbb r2_abab(a,b,i,j) t2_bbbb(c,d,l,k)
    //              += +1.00 f_bb(k,k) r2_abab(a,b,i,j)
    //              += -1.00 d-_bb(k,k) r2_abab(a,b,i,j) t0_1
    //              += +0.25 <l,k||c,d>_abab r2_abab(a,b,i,j) t2_abab(c,d,l,k)
    //              += +0.25 <k,l||c,d>_abab r2_abab(a,b,i,j) t2_abab(c,d,k,l)
    //              += +0.25 <l,k||d,c>_abab r2_abab(a,b,i,j) t2_abab(d,c,l,k)
    //              += +0.25 <k,l||d,c>_abab r2_abab(a,b,i,j) t2_abab(d,c,k,l)
    //              += -1.00 <k,c||d,a>_bbbb t2_bbbb(d,b,j,k) l2_1_bbbb(i,j,c,b)
    //              += -1.00 <l,i||c,k>_abab t2_abab(c,b,l,j) l2_1_bbbb(j,k,b,a)
    //              += +1.00 d-_bb(i,j) t0_1 l1_1_bb(j,a)
    //              += -0.50 <k,b||c,d>_aaaa t2_aaaa(c,d,j,k) l2_1_abab(j,i,b,a)
    //              += -0.50 <l,k||j,c>_abab t2_abab(b,c,l,k) l2_1_abab(j,i,b,a)
    //              += -0.50 <k,l||j,c>_abab t2_abab(b,c,k,l) l2_1_abab(j,i,b,a)
    //              += -0.50 <l,k||c,j>_aaaa t2_aaaa(c,b,l,k) l2_1_abab(j,i,b,a)
    //              += -1.00 f_aa(k,c) t2_aaaa(c,b,j,k) l2_1_abab(j,i,b,a)
    //              += +0.50 <b,k||c,d>_abab t2_abab(c,d,j,k) l2_1_abab(j,i,b,a)
    //              += +0.50 <b,k||d,c>_abab t2_abab(d,c,j,k) l2_1_abab(j,i,b,a)
    //              += -1.00 d-_aa(b,j) t0_1 l2_1_abab(j,i,b,a)
    //              += +1.00 f_bb(k,c) t2_abab(b,c,j,k) l2_1_abab(j,i,b,a)
    //              += +1.00 f_bb(i,a) l0_1
    //              += -1.00 d-_aa(b,j) l2_abab(j,i,b,a)
    //              += +1.00 <i,k||c,a>_bbbb t2_bbbb(c,b,j,k) l1_1_bb(j,b)
    //              += +1.00 d-_bb(i,k) t1_1_bb(b,j) l2_1_bbbb(j,k,b,a)
    //              += +0.50 <c,b||j,a>_abab l2_1_abab(j,i,c,b)
    //              += +0.50 <b,c||j,a>_abab l2_1_abab(j,i,b,c)
    //              += +1.00 <b,i||j,a>_abab l1_1_aa(j,b)
    //              += +1.00 d-_bb(b,j) l2_bbbb(i,j,b,a)
    //              += -0.50 <k,j||c,a>_bbbb t2_bbbb(c,b,k,j) l1_1_bb(i,b)
    //              += +0.50 <i,b||j,k>_bbbb l2_1_bbbb(j,k,b,a)
    //              += +2.00 d-_aa(k,j) t1_1_aa(b,k) l2_1_abab(j,i,b,a)
    //              += +2.00 d-_aa(k,c) t2_1_aaaa(c,b,j,k) l2_1_abab(j,i,b,a)
    //              += -2.00 d-_aa(b,c) t1_1_aa(c,j) l2_1_abab(j,i,b,a)
    //              += -2.00 d-_bb(k,c) t2_1_abab(b,c,j,k) l2_1_abab(j,i,b,a)
    //              += -0.50 <k,j||c,a>_abab t2_abab(c,b,k,j) l1_1_bb(i,b)
    //              += -0.50 <j,k||c,a>_abab t2_abab(c,b,j,k) l1_1_bb(i,b)
    //              += -2.00 d-_bb(j,b) t1_1_bb(b,j) l1_1_bb(i,a)
    //              += +1.00 <k,i||c,a>_abab t2_abab(c,b,k,j) l1_1_bb(j,b)
    //              += +1.00 d-_bb(i,j) l1_bb(j,a)
    //              += +0.50 <c,b||a,j>_bbbb l2_1_bbbb(i,j,c,b)
    //              += -0.50 <i,c||d,a>_bbbb t2_bbbb(d,b,j,k) l2_1_bbbb(j,k,c,b)
    //              += -0.50 <l,i||k,a>_abab t2_aaaa(c,b,j,l) l2_1_aaaa(j,k,c,b)
    //              += +0.50 f_bb(k,a) t2_bbbb(c,b,j,k) l2_1_bbbb(i,j,c,b)
    //              += +1.00 d-_bb(i,k) t1_1_aa(b,j) l2_1_abab(j,k,b,a)
    sigmal1_1_bb("L,a,i") -= tmps_["272_bb_Lvo"]("L,a,i");
    tmps_["272_bb_Lvo"].~TArrayD();

    // flops: o2v2L1  = o2v1L1 o3v2L1 o3v2L1 o2v2L1 o2v1L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1
    //  mems: o2v2L1  = o2v0L1 o2v2L1 o2v2L1 o2v2L1 o2v0L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["275_aabb_Lovvo"]("R,i,a,b,j")  = dp["aa_ov"]("l,d") * r1["aa"]("R,d,i") * t2["abab"]("a,b,l,j");
    tmps_["275_aabb_Lovvo"]("R,i,a,b,j") += r1["aa"]("R,a,l") * reused_["83_baba_vooo"]("b,i,j,l");
    tmps_["275_aabb_Lovvo"]("R,i,a,b,j") += dp["bb_ov"]("k,c") * r1["bb"]("R,c,j") * t2["abab"]("a,b,i,k");
    tmps_["275_aabb_Lovvo"]("R,i,a,b,j") += reused_["77_aabb_vooo"]("a,i,j,k") * r1["bb"]("R,b,k");

    // sigmar2_1_abab  = +1.00 d+_aa(k,c) r1_aa(c,i) t2_abab(a,b,k,j)
    //                += +1.00 d+_aa(k,c) r1_aa(a,k) t2_abab(c,b,i,j)
    //                += +1.00 d+_bb(k,c) r1_bb(c,j) t2_abab(a,b,i,k)
    //                += +1.00 d+_bb(k,c) r1_bb(b,k) t2_abab(a,c,i,j)
    sigmar2_1_abab("R,a,b,i,j")  = tmps_["275_aabb_Lovvo"]("R,i,a,b,j");

    // flops: o4v0L1  = o4v2L1
    //  mems: o4v0L1  = o4v0L1
    tmps_["279_abab_Loooo"]("R,i,j,k,l")  = r2["abab"]("R,a,b,i,j") * eri["abab_oovv"]("k,l,a,b");

    // flops: o2v2L1  = o3v3L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["278_baab_Lvoov"]("R,a,i,j,b")  = r2["abab"]("R,c,a,i,k") * eri["abab_oovv"]("j,k,c,b");

    // flops: o2v2L1  = o3v3L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["277_bbbb_Lvoov"]("R,a,i,j,b")  = r2["abab"]("R,c,a,k,i") * eri["abab_oovv"]("k,j,c,b");

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["276_aabb_Lvovo"]("R,a,i,b,j")  = -1.00 * r1["aa"]("R,a,i") * reused_["84_bb_vo"]("b,j");

    // flops: o3v1L1  = o3v2L1
    //  mems: o3v1L1  = o3v1L1
    tmps_["274_aabb_Lvooo"]("R,a,i,j,k")  = r2["abab"]("R,a,b,i,j") * dp["bb_ov"]("k,b");

    // flops: o3v1L1  = o3v2L1
    //  mems: o3v1L1  = o3v1L1
    tmps_["273_baba_Lvooo"]("R,a,i,j,k")  = r2["abab"]("R,b,a,i,j") * dp["aa_ov"]("k,b");

    // flops: o2v2L1  = o3v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o2v4L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o3v2L1 o4v2L1 o2v2L1 o2v2L1 o3v3L1 o2v3L1 o3v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j")  = -1.00 * r1["aa"]("R,b,k") * reused_["175_baab_vooo"]("a,k,j,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= eri["abab_vvvo"]("b,a,c,i") * r1["aa"]("R,c,j");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += reused_["57_aa_vv"]("f,b") * r2["abab"]("R,f,a,j,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= f["bb_vo"]("a,i") * r1["aa"]("R,b,j");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += r2["abab"]("R,b,d,j,i") * reused_["34_bb_vv"]("a,d");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += reused_["225_aabb_vooo"]("b,j,i,l") * r1["bb"]("R,a,l");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= r1_1["aa"]("R,b,j") * reused_["33_bb_vo"]("a,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += r0_1("R") * reused_["74_baab_vvoo"]("a,b,j,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += r1["aa"]("R,b,n") * reused_["173_bbaa_vooo"]("a,i,n,j");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= r1["bb"]("R,a,i") * reused_["89_aa_ov"]("j,b");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= f["aa_vo"]("b,j") * r1["bb"]("R,a,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += reused_["231_bbaa_vvvo"]("a,e,b,j") * r1["bb"]("R,e,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += reused_["215_baab_vovv"]("d,j,b,a") * r1["bb"]("R,d,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= r1["aa"]("R,c,j") * reused_["169_abab_vovv"]("c,i,b,a");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += r2["abab"]("R,b,e,j,m") * reused_["28_bbbb_voov"]("a,i,m,e");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= r2["abab"]("R,b,e,j,m") * reused_["26_bbbb_voov"]("a,i,m,e");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += scalars_["4"] * r2["abab"]("R,b,a,j,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += dp["bb_vo"]("a,i") * r1_1["aa"]("R,b,j");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= r1["aa"]("R,b,n") * reused_["171_abba_oovo"]("n,i,a,j");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= r2["abab"]("R,b,a,k,i") * reused_["17_aa_oo"]("k,j");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += reused_["219_abba_vooo"]("b,i,m,j") * r1["bb"]("R,a,m");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += eri["abba_vvvo"]("b,a,d,j") * r1["bb"]("R,d,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= 0.50 * r2["abab"]("R,b,a,n,i") * reused_["15_aa_oo"]("j,n");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += reused_["90_aa_vo"]("b,j") * r1["bb"]("R,a,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= f["aa_vv"]("b,c") * r2["abab"]("R,c,a,j,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += scalars_["2"] * r2_1["abab"]("R,b,a,j,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += reused_["58_aa_vo"]("b,j") * r1_1["bb"]("R,a,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= r2["aaaa"]("R,f,b,j,n") * reused_["1_bbaa_voov"]("a,i,n,f");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= r2["abab"]("R,b,a,j,l") * reused_["35_bb_oo"]("l,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= r2["abab"]("R,b,a,n,l") * reused_["80_abab_oooo"]("n,l,j,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += r0_1("R") * reused_["73_abab_vvoo"]("b,a,j,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= 0.50 * reused_["23_aa_vv"]("b,f") * r2["abab"]("R,f,a,j,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += eri["aaaa_vovo"]("b,k,c,j") * r2["abab"]("R,c,a,k,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += r1["bb"]("R,a,m") * reused_["232_bbaa_oovo"]("m,i,b,j");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += reused_["56_aabb_voov"]("b,j,m,e") * r2["bbbb"]("R,e,a,i,m");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += reused_["55_aaaa_voov"]("b,j,n,f") * r2["abab"]("R,f,a,n,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= r1["bb"]("R,a,m") * reused_["217_bbaa_oovo"]("m,i,b,j");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= eri["abba_vovo"]("b,l,d,j") * r2["bbbb"]("R,d,a,i,l");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += r1["aa"]("R,b,n") * reused_["233_aabb_oovo"]("n,j,a,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= f["bb_vv"]("a,d") * r2["abab"]("R,b,d,j,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += eri["abab_vovo"]("b,l,c,i") * r2["abab"]("R,c,a,j,l");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= reused_["54_aa_vo"]("b,j") * r1_1["bb"]("R,a,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += scalars_["3"] * r2["abab"]("R,b,a,j,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += scalars_["1"] * r2_1["abab"]("R,b,a,j,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += r2["abab"]("R,b,e,j,i") * reused_["3_bb_vv"]("e,a");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += reused_["152_aa_vo"]("b,j") * r1["bb"]("R,a,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += r2["abab"]("R,b,a,j,m") * reused_["30_bb_oo"]("m,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += r1["aa"]("R,b,j") * reused_["88_bb_vo"]("a,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= eri["abab_oooo"]("n,l,j,i") * r2["abab"]("R,b,a,n,l");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += scalars_["5"] * r2["abab"]("R,b,a,j,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += dp["aa_vo"]("b,j") * r1_1["bb"]("R,a,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += r2["abab"]("R,b,a,n,i") * reused_["16_aa_oo"]("n,j");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= reused_["168_abba_vovv"]("b,i,a,f") * r1["aa"]("R,f,j");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += r2["aaaa"]("R,f,b,j,n") * reused_["27_bbaa_voov"]("a,i,n,f");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= eri["abab_vvvv"]("b,a,c,e") * r2["abab"]("R,c,e,j,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= dp["bb_oo"]("l,i") * r2_1["abab"]("R,b,a,j,l");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += eri["baba_vovo"]("a,k,d,j") * r2["abab"]("R,b,d,k,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= eri["baab_vovo"]("a,k,c,i") * r2["aaaa"]("R,c,b,j,k");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= 0.50 * reused_["87_aa_vo"]("b,j") * r1["bb"]("R,a,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += reused_["59_aa_vv"]("b,c") * r2["abab"]("R,c,a,j,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += reused_["209_abba_vvvo"]("b,e,a,j") * r1["bb"]("R,e,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += eri["abab_vooo"]("b,l,j,i") * r1["bb"]("R,a,l");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += reused_["167_aabb_vvvo"]("b,f,a,i") * r1["aa"]("R,f,j");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += f["aa_oo"]("k,j") * r2["abab"]("R,b,a,k,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= 0.50 * r2["abab"]("R,b,e,j,i") * reused_["29_bb_vv"]("a,e");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += eri["bbbb_vovo"]("a,l,d,i") * r2["abab"]("R,b,d,j,l");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += f["bb_oo"]("l,i") * r2["abab"]("R,b,a,j,l");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += r1["aa"]("R,b,j") * reused_["91_bb_vo"]("a,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= 0.50 * r2["abab"]("R,b,a,j,m") * reused_["31_bb_oo"]("i,m");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= r2["abab"]("R,b,e,n,i") * reused_["78_baab_voov"]("a,j,n,e");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += r1_1["aa"]("R,b,j") * reused_["32_bb_vo"]("a,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += dp["bb_vv"]("a,d") * r2_1["abab"]("R,b,d,j,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += dp["aa_vv"]("b,c") * r2_1["abab"]("R,c,a,j,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += reused_["230_aabb_vvvo"]("b,f,a,i") * r1["aa"]("R,f,j");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= eri["baab_vooo"]("a,k,j,i") * r1["aa"]("R,b,k");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= dp["aa_oo"]("k,j") * r2_1["abab"]("R,b,a,k,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= reused_["6_aabb_voov"]("b,j,m,e") * r2["bbbb"]("R,e,a,i,m");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= reused_["210_bbaa_vvvo"]("a,e,b,j") * r1["bb"]("R,e,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= tmps_["276_aabb_Lvovo"]("R,b,j,a,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= t1_1["aa"]("b,k") * tmps_["273_baba_Lvooo"]("R,a,j,i,k");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= t2["abab"]("b,a,n,l") * tmps_["279_abab_Loooo"]("R,j,i,n,l");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= t1_1["bb"]("a,i") * tmps_["192_aa_Lov"]("R,j,b");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= t0_1 * tmps_["275_aabb_Lovvo"]("R,j,b,a,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= t2["abab"]("b,d,j,l") * tmps_["277_bbbb_Lvoov"]("R,a,i,l,d");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= t2["abab"]("c,a,j,i") * tmps_["134_aa_Lvv"]("R,b,c");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += 0.50 * t2["abab"]("b,a,k,i") * tmps_["36_aa_Loo"]("R,j,k");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += t2["abab"]("b,a,j,l") * tmps_["54_bb_Loo"]("R,l,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += t2_1["abab"]("b,a,j,i") * tmps_["45_L"]("R");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += r1["aa"]("R,b,j") * reused_["180_bb_vo"]("a,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += t1_1["aa"]("b,j") * tmps_["213_bb_Lvo"]("R,a,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += 0.50 * t2["abab"]("c,a,j,i") * tmps_["135_aa_Lvv"]("R,b,c");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= t2["abab"]("b,a,k,i") * tmps_["34_aa_Loo"]("R,j,k");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= tmps_["189_aa_Lvo"]("R,b,j") * t1_1["bb"]("a,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= t2["abab"]("b,d,k,i") * tmps_["278_baab_Lvoov"]("R,a,j,k,d");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += 0.50 * t2["abab"]("b,d,j,i") * tmps_["64_bb_Lvv"]("R,a,d");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= t1_1["aa"]("b,j") * tmps_["214_bb_Lov"]("R,i,a");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= t2_1["abab"]("b,a,k,i") * tmps_["35_aa_Loo"]("R,j,k");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += t1_1["aa"]("b,j") * tmps_["215_bb_Lvo"]("R,a,i");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= tmps_["274_aabb_Lvooo"]("R,b,j,i,l") * t1_1["bb"]("a,l");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= t2["abab"]("b,a,k,i") * tmps_["37_aa_Loo"]("R,k,j");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += t2["abab"]("b,a,j,l") * tmps_["56_bb_Loo"]("R,i,l");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += t2["abab"]("b,d,j,i") * tmps_["65_bb_Lvv"]("R,a,d");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= t2_1["abab"]("b,a,j,l") * tmps_["73_bb_Loo"]("R,i,l");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") -= t2["abab"]("b,a,j,l") * tmps_["49_bb_Loo"]("R,i,l");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += t1_1["bb"]("a,i") * tmps_["191_aa_Lvo"]("R,b,j");
    tmps_["280_bbaa_Lvovo"]("R,a,i,b,j") += r1["bb"]("R,a,i") * reused_["229_aa_vo"]("b,j");
    tmps_["276_aabb_Lvovo"].~TArrayD();
    tmps_["275_aabb_Lovvo"].~TArrayD();
    tmps_["274_aabb_Lvooo"].~TArrayD();
    tmps_["273_baba_Lvooo"].~TArrayD();
    tmps_["215_bb_Lvo"].~TArrayD();
    tmps_["214_bb_Lov"].~TArrayD();
    tmps_["213_bb_Lvo"].~TArrayD();
    tmps_["192_aa_Lov"].~TArrayD();
    tmps_["191_aa_Lvo"].~TArrayD();
    tmps_["189_aa_Lvo"].~TArrayD();
    tmps_["45_L"].~TArrayD();
    tmps_["35_aa_Loo"].~TArrayD();
}
