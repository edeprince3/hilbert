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

void hilbert::EOM_EE_QED_CCSD_21::sigma_ee_21_8() {

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


    // flops: o2v2L1  = o3v3L1 o3v3L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1
    tmps_["301_bbbb_Lvoov"]("L,a,i,j,b")  = t2["abab"]("c,a,k,i") * l2["abab"]("L,k,j,c,b");
    tmps_["301_bbbb_Lvoov"]("L,a,i,j,b") += t2_1["abab"]("c,a,k,i") * l2_1["abab"]("L,k,j,c,b");

    // flops: o2v2L1  = o3v3L1 o3v3L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1
    tmps_["300_baab_Lvoov"]("L,a,i,j,b")  = t2["abab"]("c,a,i,k") * l2["abab"]("L,j,k,c,b");
    tmps_["300_baab_Lvoov"]("L,a,i,j,b") += t2_1["abab"]("c,a,i,k") * l2_1["abab"]("L,j,k,c,b");

    // flops: o4v0L1  = o4v2L1 o4v2L1 o4v0L1
    //  mems: o4v0L1  = o4v0L1 o4v0L1 o4v0L1
    tmps_["299_abab_Loooo"]("L,i,j,k,l")  = l2["abab"]("L,i,j,a,b") * t2["abab"]("a,b,k,l");
    tmps_["299_abab_Loooo"]("L,i,j,k,l") += l2_1["abab"]("L,i,j,a,b") * t2_1["abab"]("a,b,k,l");

    // flops: o3v1L1  = o3v2L1
    //  mems: o3v1L1  = o3v1L1
    tmps_["298_abab_Loovo"]("L,i,j,a,k")  = l2["abab"]("L,i,j,a,b") * t1_1["bb"]("b,k");

    // flops: o3v1L1  = o3v2L1
    //  mems: o3v1L1  = o3v1L1
    tmps_["297_abba_Loovo"]("L,i,j,a,k")  = l2["abab"]("L,i,j,b,a") * t1_1["aa"]("b,k");

    // flops: o2v2L1  = o2v2L1 o3v2L1 o4v2L1 o2v0L1 o3v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o4v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o4v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v4L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o4v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o4v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o4v2L1 o2v2L1 o4v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v0L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b")  = -1.00 * tmps_["170_aa_Lov"]("L,j,b") * dp["bb_ov"]("i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += tmps_["209_aa_Loo"]("L,j,m") * eri["abab_oovv"]("m,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += tmps_["138_aaaa_Loovo"]("L,j,m,b,k") * eri["abba_oovo"]("k,i,a,m");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= 0.50 * (tmps_["219_bb_Loo"]("L,i,n") + -2.00 * tmps_["218_bb_Loo"]("L,i,n")) * eri["abab_oovv"]("j,n,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= eri["abab_oovv"]("j,n,b,f") * tmps_["301_bbbb_Lvoov"]("L,f,n,i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= tmps_["300_baab_Lvoov"]("L,f,k,j,a") * eri["abab_oovv"]("k,i,b,f");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += eri["abab_oovv"]("j,l,b,a") * tmps_["222_bb_Loo"]("L,i,l");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += tmps_["114_bb_Lvv"]("L,a,f") * eri["abab_oovv"]("j,i,b,f");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= t0_1 * tmps_["260_aabb_Lovvo"]("L,j,b,a,i");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += eri["abab_oovv"]("j,n,b,a") * tmps_["217_bb_Loo"]("L,i,n");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= eri["abab_oovv"]("k,l,b,a") * tmps_["299_abab_Loooo"]("L,j,i,k,l");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= 0.50 * eri["abab_oovv"]("j,n,b,a") * tmps_["220_bb_Loo"]("L,i,n");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += dp["aa_ov"]("j,b") * tmps_["182_bb_Lov"]("L,i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= reused_["12_bb_ov"]("l,a") * tmps_["169_abab_Loovo"]("L,j,i,b,l");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += eri["abab_vovv"]("c,l,b,a") * tmps_["169_abab_Loovo"]("L,j,i,c,l");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += f["bb_ov"]("l,a") * tmps_["169_abab_Loovo"]("L,j,i,b,l");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= dp["aa_ov"]("m,b") * tmps_["297_abba_Loovo"]("L,j,i,a,m");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= t0_1 * tmps_["239_abab_Loovv"]("L,j,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= dp["bb_ov"]("l,a") * tmps_["298_abab_Loovo"]("L,j,i,b,l");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= dp["aa_ov"]("j,b") * tmps_["183_bb_Lov"]("L,i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= l2["aaaa"]("L,j,m,e,b") * reused_["6_aabb_voov"]("e,m,i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += l2["abab"]("L,j,l,b,a") * reused_["30_bb_oo"]("i,l");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += l1["aa"]("L,j,b") * reused_["12_bb_ov"]("i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += eri["baba_vovo"]("d,j,a,m") * l2["abab"]("L,m,i,b,d");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += scalars_["4"] * l2["abab"]("L,j,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += 0.50 * reused_["109_bb_vv"]("a,d") * l2["abab"]("L,j,i,b,d");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= l2["abab"]("L,m,n,b,a") * reused_["80_abab_oooo"]("j,i,m,n");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += dp["aa_ov"]("j,b") * l1_1["bb"]("L,i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= eri["abab_vovv"]("e,i,b,a") * l1["aa"]("L,j,e");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += l1_1["aa"]("L,m,b") * reused_["116_abba_oovo"]("j,i,a,m");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += eri["abab_oovo"]("j,i,b,l") * l1["bb"]("L,l,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += reused_["3_bb_vv"]("a,d") * l2["abab"]("L,j,i,b,d");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += reused_["69_aaaa_voov"]("e,m,j,b") * l2["abab"]("L,m,i,e,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= l2["abab"]("L,j,l,b,d") * reused_["26_bbbb_voov"]("d,l,i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= scalars_["6"] * l2_1["abab"]("L,j,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += dp["aa_vv"]("e,b") * l2_1["abab"]("L,j,i,e,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += scalars_["1"] * l2_1["abab"]("L,j,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += l1_1["bb"]("L,l,a") * reused_["43_abab_oovo"]("j,i,b,l");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += reused_["57_aa_vv"]("b,e") * l2["abab"]("L,j,i,e,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += reused_["159_bbaa_voov"]("d,l,j,b") * l2_1["bbbb"]("L,i,l,d,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= f["bb_ov"]("i,a") * l1["aa"]("L,j,b");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += scalars_["5"] * l2["abab"]("L,j,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += reused_["59_aa_vv"]("e,b") * l2["abab"]("L,j,i,e,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= reused_["182_aa_oo"]("j,m") * l2_1["abab"]("L,m,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= l1_1["aa"]("L,j,b") * reused_["148_bb_ov"]("i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= reused_["17_aa_oo"]("j,m") * l2["abab"]("L,m,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= reused_["147_aa_ov"]("j,b") * l1_1["bb"]("L,i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= reused_["9_aa_ov"]("j,b") * l1_1["bb"]("L,i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= eri["abba_vovo"]("e,i,a,m") * l2["aaaa"]("L,j,m,e,b");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += reused_["121_aa_oo"]("m,j") * l2_1["abab"]("L,m,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,j,l,b,a") * reused_["48_bb_oo"]("i,l");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= reused_["132_aa_vv"]("e,b") * l2_1["abab"]("L,j,i,e,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= f["aa_ov"]("j,b") * l1["bb"]("L,i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += 0.50 * reused_["174_aa_oo"]("j,m") * l2["abab"]("L,m,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += 0.50 * reused_["18_aa_vv"]("b,e") * l2["abab"]("L,j,i,e,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += reused_["158_bbaa_voov"]("d,l,j,b") * l2["bbbb"]("L,i,l,d,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= l2["abab"]("L,j,l,b,a") * reused_["35_bb_oo"]("i,l");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,k,l,b,a") * reused_["244_abab_oooo"]("j,i,k,l");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= reused_["123_aa_vv"]("e,b") * l2_1["abab"]("L,j,i,e,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += dp["bb_ov"]("i,a") * l1_1["aa"]("L,j,b");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= l2_1["abab"]("L,m,n,b,a") * reused_["245_abab_oooo"]("m,n,j,i");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= eri["abba_oovo"]("j,i,a,m") * l1["aa"]("L,m,b");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= reused_["78_baab_voov"]("d,m,j,a") * l2["abab"]("L,m,i,b,d");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += l2_1["aaaa"]("L,j,m,e,b") * reused_["162_aabb_voov"]("e,m,i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= l1_1["aa"]("L,j,b") * reused_["8_bb_ov"]("i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += scalars_["2"] * l2_1["abab"]("L,j,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += 0.50 * l2["abab"]("L,j,l,b,a") * reused_["161_bb_oo"]("i,l");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += dp["bb_vv"]("d,a") * l2_1["abab"]("L,j,i,b,d");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= reused_["1_bbaa_voov"]("d,l,j,b") * l2["bbbb"]("L,i,l,d,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= reused_["105_baab_vovo"]("d,j,b,l") * l2_1["bbbb"]("L,i,l,d,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= eri["abab_oooo"]("j,i,m,n") * l2["abab"]("L,m,n,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += reused_["13_aa_ov"]("j,b") * l1["bb"]("L,i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += f["bb_oo"]("i,l") * l2["abab"]("L,j,l,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += reused_["34_bb_vv"]("d,a") * l2["abab"]("L,j,i,b,d");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += eri["baab_vovv"]("d,j,b,a") * l1["bb"]("L,i,d");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= reused_["103_bbaa_voov"]("d,l,j,b") * l2_1["bbbb"]("L,i,l,d,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,j,l,b,a") * reused_["40_bb_oo"]("i,l");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += eri["aaaa_vovo"]("e,j,b,m") * l2["abab"]("L,m,i,e,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,j,l,e,a") * reused_["242_abab_vovo"]("e,i,b,l");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += reused_["108_bb_vv"]("d,a") * l2_1["abab"]("L,j,i,b,d");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= dp["bb_oo"]("i,l") * l2_1["abab"]("L,j,l,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= eri["baab_vovo"]("d,j,b,l") * l2["bbbb"]("L,i,l,d,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= dp["aa_oo"]("j,m") * l2_1["abab"]("L,m,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,j,l,b,d") * reused_["164_bbbb_voov"]("d,l,i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= f["aa_vv"]("e,b") * l2["abab"]("L,j,i,e,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= reused_["112_bb_vv"]("d,a") * l2_1["abab"]("L,j,i,b,d");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += f["aa_oo"]("j,m") * l2["abab"]("L,m,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= reused_["243_baba_vovo"]("d,j,a,m") * l2_1["abab"]("L,m,i,b,d");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += l2_1["aaaa"]("L,j,m,e,b") * reused_["118_abba_vovo"]("e,i,a,m");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= eri["abab_vvvv"]("c,d,b,a") * l2["abab"]("L,j,i,c,d");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= f["bb_vv"]("d,a") * l2["abab"]("L,j,i,b,d");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= l2_1["abab"]("L,j,l,b,a") * reused_["181_bb_oo"]("i,l");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= l2_1["abab"]("L,j,l,b,d") * reused_["102_bbbb_voov"]("d,l,i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += reused_["16_aa_oo"]("j,m") * l2["abab"]("L,m,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= l2_1["aaaa"]("L,j,m,e,b") * reused_["114_aabb_voov"]("e,m,i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= reused_["241_baab_voov"]("d,m,j,a") * l2_1["abab"]("L,m,i,b,d");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= l2["abab"]("L,j,l,b,a") * reused_["38_bb_oo"]("i,l");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += l2["abab"]("L,j,l,b,d") * reused_["163_bbbb_voov"]("d,l,i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += scalars_["3"] * l2["abab"]("L,j,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += eri["abab_vovo"]("e,i,b,l") * l2["abab"]("L,j,l,e,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += l2["aaaa"]("L,j,m,e,b") * reused_["70_aabb_voov"]("e,m,i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += eri["bbbb_vovo"]("d,i,a,l") * l2["abab"]("L,j,l,b,d");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= reused_["125_aa_oo"]("j,m") * l2_1["abab"]("L,m,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += reused_["160_aaaa_vovo"]("e,j,b,m") * l2_1["abab"]("L,m,i,e,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= reused_["120_aa_oo"]("j,m") * l2["abab"]("L,m,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= reused_["13_aa_ov"]("m,b") * tmps_["181_abba_Loovo"]("L,j,i,a,m");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += eri["aaaa_oovo"]("j,k,b,m") * tmps_["181_abba_Loovo"]("L,m,i,a,k");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= 0.50 * tmps_["115_bb_Lvv"]("L,a,f") * eri["abab_oovv"]("j,i,b,f");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += tmps_["140_aa_Lvv"]("L,b,c") * eri["abab_oovv"]("j,i,c,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= eri["abab_oovo"]("j,n,b,l") * tmps_["113_bbbb_Loovo"]("L,i,l,a,n");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += f["aa_ov"]("m,b") * tmps_["181_abba_Loovo"]("L,j,i,a,m");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= 0.50 * tmps_["139_aa_Lvv"]("L,b,c") * eri["abab_oovv"]("j,i,c,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= eri["baab_vovv"]("f,m,b,a") * tmps_["181_abba_Loovo"]("L,j,i,f,m");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += eri["abba_oovo"]("j,n,a,m") * tmps_["169_abab_Loovo"]("L,m,i,b,n");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += tmps_["169_abab_Loovo"]("L,j,l,b,n") * eri["bbbb_oovo"]("i,n,a,l");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += tmps_["207_aa_Loo"]("L,j,k") * eri["abab_oovv"]("k,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += tmps_["202_aa_Loo"]("L,j,k") * reused_["246_abab_oovv"]("k,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") -= tmps_["181_abba_Loovo"]("L,j,l,a,k") * eri["abab_oovo"]("k,i,b,l");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += tmps_["171_aa_Lov"]("L,j,b") * dp["bb_ov"]("i,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += tmps_["206_aa_Loo"]("L,j,k") * eri["abab_oovv"]("k,i,b,a");
    tmps_["302_bbaa_Lovov"]("L,i,a,j,b") += tmps_["203_aa_Loo"]("L,j,k") * reused_["246_abab_oovv"]("k,i,b,a");
    tmps_["301_bbbb_Lvoov"].~TArrayD();
    tmps_["300_baab_Lvoov"].~TArrayD();
    tmps_["299_abab_Loooo"].~TArrayD();
    tmps_["298_abab_Loovo"].~TArrayD();
    tmps_["297_abba_Loovo"].~TArrayD();
    tmps_["260_aabb_Lovvo"].~TArrayD();
    tmps_["239_abab_Loovv"].~TArrayD();
    tmps_["222_bb_Loo"].~TArrayD();
    tmps_["220_bb_Loo"].~TArrayD();
    tmps_["219_bb_Loo"].~TArrayD();
    tmps_["218_bb_Loo"].~TArrayD();
    tmps_["217_bb_Loo"].~TArrayD();
    tmps_["209_aa_Loo"].~TArrayD();
    tmps_["207_aa_Loo"].~TArrayD();
    tmps_["206_aa_Loo"].~TArrayD();
    tmps_["203_aa_Loo"].~TArrayD();
    tmps_["202_aa_Loo"].~TArrayD();
    tmps_["183_bb_Lov"].~TArrayD();
    tmps_["182_bb_Lov"].~TArrayD();
    tmps_["171_aa_Lov"].~TArrayD();
    tmps_["170_aa_Lov"].~TArrayD();
    tmps_["140_aa_Lvv"].~TArrayD();
    tmps_["139_aa_Lvv"].~TArrayD();
    tmps_["138_aaaa_Loovo"].~TArrayD();
    tmps_["115_bb_Lvv"].~TArrayD();
    tmps_["114_bb_Lvv"].~TArrayD();
    tmps_["113_bbbb_Loovo"].~TArrayD();

    // sigmal2_abab  = -1.00 <j,l||b,k>_bbbb t1_1_bb(c,l) l2_1_abab(i,k,a,c)
    //              += +0.50 <i,l||a,b>_abab t2_1_bbbb(d,c,k,l) l2_1_bbbb(j,k,d,c)
    //              += -0.50 <i,l||a,b>_abab t2_abab(d,c,k,l) l2_abab(k,j,d,c)
    //              += -0.50 <i,l||a,b>_abab t2_abab(c,d,k,l) l2_abab(k,j,c,d)
    //              += +1.00 <i,l||a,d>_abab t2_1_abab(c,d,k,l) l2_1_abab(k,j,c,b)
    //              += +1.00 <i,l||a,d>_abab t2_abab(c,d,k,l) l2_abab(k,j,c,b)
    //              += +1.00 <l,j||a,d>_abab t2_1_abab(c,d,l,k) l2_1_abab(i,k,c,b)
    //              += +1.00 <l,j||a,d>_abab t2_abab(c,d,l,k) l2_abab(i,k,c,b)
    //              += -1.00 <i,k||a,b>_abab t1_1_bb(c,k) l1_1_bb(j,c)
    //              += -0.50 <i,j||a,d>_abab t2_1_abab(c,d,k,l) l2_1_abab(k,l,c,b)
    //              += -0.50 <i,j||a,d>_abab t2_1_abab(c,d,l,k) l2_1_abab(l,k,c,b)
    //              += -0.50 <i,j||a,d>_abab t2_abab(c,d,k,l) l2_abab(k,l,c,b)
    //              += -0.50 <i,j||a,d>_abab t2_abab(c,d,l,k) l2_abab(l,k,c,b)
    //              += +1.00 d-_bb(j,c) t0_1 t1_1_bb(c,k) l2_1_abab(i,k,a,b)
    //              += -0.50 <i,l||a,b>_abab t2_1_abab(d,c,k,l) l2_1_abab(k,j,d,c)
    //              += -0.50 <i,l||a,b>_abab t2_1_abab(c,d,k,l) l2_1_abab(k,j,c,d)
    //              += +0.25 <l,k||a,b>_abab t2_abab(d,c,l,k) l2_abab(i,j,d,c)
    //              += +0.25 <k,l||a,b>_abab t2_abab(d,c,k,l) l2_abab(i,j,d,c)
    //              += +0.25 <l,k||a,b>_abab t2_abab(c,d,l,k) l2_abab(i,j,c,d)
    //              += +0.25 <k,l||a,b>_abab t2_abab(c,d,k,l) l2_abab(i,j,c,d)
    //              += +0.25 <l,k||a,b>_abab t2_1_abab(d,c,l,k) l2_1_abab(i,j,d,c)
    //              += +0.25 <k,l||a,b>_abab t2_1_abab(d,c,k,l) l2_1_abab(i,j,d,c)
    //              += +0.25 <l,k||a,b>_abab t2_1_abab(c,d,l,k) l2_1_abab(i,j,c,d)
    //              += +0.25 <k,l||a,b>_abab t2_1_abab(c,d,k,l) l2_1_abab(i,j,c,d)
    //              += +0.50 <i,l||a,b>_abab t2_bbbb(d,c,k,l) l2_bbbb(j,k,d,c)
    //              += -1.00 d-_aa(i,a) t1_1_aa(c,k) l2_abab(k,j,c,b)
    //              += +1.00 d-_bb(k,b) t0_1 t1_1_bb(c,k) l2_1_abab(i,j,a,c)
    //              += -1.00 <d,k||a,b>_abab t1_1_bb(c,k) l2_1_abab(i,j,d,c)
    //              += -1.00 f_bb(k,b) t1_1_bb(c,k) l2_1_abab(i,j,a,c)
    //              += +1.00 d-_aa(k,a) t1_1_aa(c,k) l2_abab(i,j,c,b)
    //              += +1.00 d-_aa(i,c) t0_1 t1_1_aa(c,k) l2_1_abab(k,j,a,b)
    //              += +1.00 d-_bb(k,b) t1_1_bb(c,k) l2_abab(i,j,a,c)
    //              += +1.00 d-_aa(i,a) t1_1_bb(c,k) l2_bbbb(j,k,c,b)
    //              += +1.00 <l,j||d,b>_abab t2_aaaa(d,c,k,l) l2_aaaa(i,k,c,a)
    //              += -0.50 <l,j||c,d>_abab t2_abab(c,d,l,k) l2_abab(i,k,a,b)
    //              += -0.50 <l,j||d,c>_abab t2_abab(d,c,l,k) l2_abab(i,k,a,b)
    //              += -1.00 d-_bb(j,b) t0_1 l1_aa(i,a)
    //              += -1.00 <i,c||k,b>_abab l2_abab(k,j,a,c)
    //              += -1.00 d-_aa(k,c) t1_1_aa(c,k) l2_abab(i,j,a,b)
    //              += -0.50 <l,k||d,b>_bbbb t2_bbbb(d,c,l,k) l2_abab(i,j,a,c)
    //              += +0.25 <i,j||c,d>_abab t2_abab(c,d,k,l) l2_abab(k,l,a,b)
    //              += +0.25 <i,j||c,d>_abab t2_abab(c,d,l,k) l2_abab(l,k,a,b)
    //              += +0.25 <i,j||d,c>_abab t2_abab(d,c,k,l) l2_abab(k,l,a,b)
    //              += +0.25 <i,j||d,c>_abab t2_abab(d,c,l,k) l2_abab(l,k,a,b)
    //              += -1.00 d+_aa(i,a) l1_1_bb(j,b)
    //              += +1.00 <c,j||a,b>_abab l1_aa(i,c)
    //              += -1.00 <i,j||c,b>_abab t1_1_aa(c,k) l1_1_aa(k,a)
    //              += -1.00 <i,j||a,k>_abab l1_bb(k,b)
    //              += -0.50 <l,k||d,b>_abab t2_abab(d,c,l,k) l2_abab(i,j,a,c)
    //              += -0.50 <k,l||d,b>_abab t2_abab(d,c,k,l) l2_abab(i,j,a,c)
    //              += +1.00 <i,l||d,a>_aaaa t2_aaaa(d,c,k,l) l2_abab(k,j,c,b)
    //              += +1.00 <l,j||d,b>_abab t2_abab(d,c,l,k) l2_abab(i,k,a,c)
    //              += +1.00 t0_1 l2_1_abab(i,j,a,b) w0
    //              += +0.25 <l,k||c,d>_aaaa t2_1_aaaa(c,d,l,k) l2_1_abab(i,j,a,b)
    //              += +0.25 <l,k||c,d>_abab t2_1_abab(c,d,l,k) l2_1_abab(i,j,a,b)
    //              += +0.25 <k,l||c,d>_abab t2_1_abab(c,d,k,l) l2_1_abab(i,j,a,b)
    //              += +0.25 <l,k||d,c>_abab t2_1_abab(d,c,l,k) l2_1_abab(i,j,a,b)
    //              += +0.25 <k,l||d,c>_abab t2_1_abab(d,c,k,l) l2_1_abab(i,j,a,b)
    //              += -1.00 d-_aa(k,c) t0_1 t1_1_aa(c,k) l2_1_abab(i,j,a,b)
    //              += +1.00 f_bb(k,c) t1_1_bb(c,k) l2_1_abab(i,j,a,b)
    //              += -1.00 d-_bb(k,c) t0_1 t1_1_bb(c,k) l2_1_abab(i,j,a,b)
    //              += +0.25 <l,k||c,d>_bbbb t2_1_bbbb(c,d,l,k) l2_1_abab(i,j,a,b)
    //              += +1.00 f_aa(k,c) t1_1_aa(c,k) l2_1_abab(i,j,a,b)
    //              += -1.00 d+_aa(c,a) l2_1_abab(i,j,c,b)
    //              += -1.00 d+_aa(k,k) l2_1_abab(i,j,a,b)
    //              += -1.00 <i,j||a,c>_abab t1_1_bb(c,k) l1_1_bb(k,b)
    //              += -0.50 <l,k||a,d>_abab t2_abab(c,d,l,k) l2_abab(i,j,c,b)
    //              += -0.50 <k,l||a,d>_abab t2_abab(c,d,k,l) l2_abab(i,j,c,b)
    //              += +1.00 <i,l||d,a>_aaaa t2_1_abab(d,c,l,k) l2_1_bbbb(j,k,c,b)
    //              += +1.00 f_bb(j,b) l1_aa(i,a)
    //              += -1.00 d-_aa(k,k) t0_1 l2_abab(i,j,a,b)
    //              += +1.00 f_aa(j,j) r1_1_bb(a,i)
    //              += -0.50 <k,j||k,j>_abab r1_1_bb(a,i)
    //              += -0.50 <j,k||j,k>_abab r1_1_bb(a,i)
    //              += +0.25 <k,j||b,c>_aaaa r1_1_bb(a,i) t2_aaaa(b,c,k,j)
    //              += -0.50 <k,j||k,j>_bbbb r1_1_bb(a,i)
    //              += -0.50 <k,j||k,j>_aaaa r1_1_bb(a,i)
    //              += +0.25 <k,j||b,c>_bbbb r1_1_bb(a,i) t2_bbbb(b,c,k,j)
    //              += +1.00 f_bb(j,j) r1_1_bb(a,i)
    //              += -1.00 d-_bb(j,j) r1_1_bb(a,i) t0_1
    //              += +0.25 <k,j||b,c>_abab r1_1_bb(a,i) t2_abab(b,c,k,j)
    //              += +0.25 <j,k||b,c>_abab r1_1_bb(a,i) t2_abab(b,c,j,k)
    //              += +0.25 <k,j||c,b>_abab r1_1_bb(a,i) t2_abab(c,b,k,j)
    //              += +0.25 <j,k||c,b>_abab r1_1_bb(a,i) t2_abab(c,b,j,k)
    //              += -1.00 d-_aa(c,a) t0_1 l2_abab(i,j,c,b)
    //              += +1.00 <i,l||c,k>_aaaa t1_1_aa(c,l) l2_1_abab(k,j,a,b)
    //              += -0.50 <i,l||c,d>_aaaa t2_1_aaaa(c,d,k,l) l2_1_abab(k,j,a,b)
    //              += -1.00 <j,k||c,b>_bbbb t1_1_bb(c,k) l1_1_aa(i,a)
    //              += +1.00 d-_aa(i,k) t0_1 l2_abab(k,j,a,b)
    //              += -1.00 <i,k||c,a>_aaaa t1_1_aa(c,k) l1_1_bb(j,b)
    //              += +1.00 <i,k||a,c>_abab t1_1_bb(c,k) l1_1_bb(j,b)
    //              += -1.00 <c,j||k,b>_abab l2_aaaa(i,k,c,a)
    //              += -1.00 f_aa(i,c) t1_1_aa(c,k) l2_1_abab(k,j,a,b)
    //              += -1.00 <l,j||c,k>_abab t1_1_aa(c,l) l2_1_abab(i,k,a,b)
    //              += -0.50 <l,j||c,d>_abab t2_1_abab(c,d,l,k) l2_1_abab(i,k,a,b)
    //              += -0.50 <l,j||d,c>_abab t2_1_abab(d,c,l,k) l2_1_abab(i,k,a,b)
    //              += +1.00 <k,c||d,a>_aaaa t1_1_aa(d,k) l2_1_abab(i,j,c,b)
    //              += -0.50 <l,k||d,a>_aaaa t2_1_aaaa(d,c,l,k) l2_1_abab(i,j,c,b)
    //              += +1.00 f_aa(i,a) l1_bb(j,b)
    //              += -0.50 <i,l||c,d>_aaaa t2_aaaa(c,d,k,l) l2_abab(k,j,a,b)
    //              += -0.50 <l,k||d,a>_aaaa t2_aaaa(d,c,l,k) l2_abab(i,j,c,b)
    //              += +1.00 <i,l||d,a>_aaaa t2_abab(d,c,l,k) l2_bbbb(j,k,c,b)
    //              += +1.00 d-_bb(j,k) t0_1 l2_abab(i,k,a,b)
    //              += +1.00 <i,j||l,c>_abab t1_1_bb(c,k) l2_1_abab(l,k,a,b)
    //              += +1.00 <c,k||a,d>_abab t1_1_bb(d,k) l2_1_abab(i,j,c,b)
    //              += -0.50 <l,k||a,d>_abab t2_1_abab(c,d,l,k) l2_1_abab(i,j,c,b)
    //              += -0.50 <k,l||a,d>_abab t2_1_abab(c,d,k,l) l2_1_abab(i,j,c,b)
    //              += -1.00 d+_bb(j,b) l1_1_aa(i,a)
    //              += +0.25 <i,j||c,d>_abab t2_1_abab(c,d,k,l) l2_1_abab(k,l,a,b)
    //              += +0.25 <i,j||c,d>_abab t2_1_abab(c,d,l,k) l2_1_abab(l,k,a,b)
    //              += +0.25 <i,j||d,c>_abab t2_1_abab(d,c,k,l) l2_1_abab(k,l,a,b)
    //              += +0.25 <i,j||d,c>_abab t2_1_abab(d,c,l,k) l2_1_abab(l,k,a,b)
    //              += +1.00 <i,j||c,l>_abab t1_1_aa(c,k) l2_1_abab(k,l,a,b)
    //              += -1.00 <i,j||k,b>_abab l1_aa(k,a)
    //              += +1.00 <i,l||d,b>_abab t2_abab(d,c,k,l) l2_abab(k,j,a,c)
    //              += +1.00 <j,l||d,b>_bbbb t2_1_abab(c,d,k,l) l2_1_aaaa(i,k,c,a)
    //              += +1.00 <k,j||c,b>_abab t1_1_aa(c,k) l1_1_aa(i,a)
    //              += -1.00 d+_bb(k,k) l2_1_abab(i,j,a,b)
    //              += -0.50 <j,l||c,d>_bbbb t2_bbbb(c,d,k,l) l2_abab(i,k,a,b)
    //              += -1.00 d+_bb(c,b) l2_1_abab(i,j,a,c)
    //              += +1.00 <i,l||a,d>_abab t2_bbbb(d,c,k,l) l2_bbbb(j,k,c,b)
    //              += -1.00 <i,c||a,d>_abab t1_1_bb(d,k) l2_1_bbbb(j,k,c,b)
    //              += +0.50 <i,j||k,l>_abab l2_abab(k,l,a,b)
    //              += +0.50 <i,j||l,k>_abab l2_abab(l,k,a,b)
    //              += -1.00 d-_aa(i,a) t0_1 l1_bb(j,b)
    //              += -1.00 f_bb(j,k) l2_abab(i,k,a,b)
    //              += -1.00 d-_bb(c,b) t0_1 l2_abab(i,j,a,c)
    //              += +1.00 <i,c||a,b>_abab l1_bb(j,c)
    //              += +1.00 <i,l||a,d>_abab t2_1_bbbb(d,c,k,l) l2_1_bbbb(j,k,c,b)
    //              += -1.00 f_bb(j,c) t1_1_bb(c,k) l2_1_abab(i,k,a,b)
    //              += +1.00 <i,c||a,k>_aaaa l2_abab(k,j,c,b)
    //              += -1.00 <c,j||a,d>_abab t1_1_bb(d,k) l2_1_abab(i,k,c,b)
    //              += +1.00 <k,c||d,b>_abab t1_1_aa(d,k) l2_1_abab(i,j,a,c)
    //              += -0.50 <l,k||d,b>_abab t2_1_abab(d,c,l,k) l2_1_abab(i,j,a,c)
    //              += -0.50 <k,l||d,b>_abab t2_1_abab(d,c,k,l) l2_1_abab(i,j,a,c)
    //              += +1.00 d+_bb(j,k) l2_1_abab(i,k,a,b)
    //              += -1.00 <i,c||a,k>_abab l2_bbbb(j,k,c,b)
    //              += +1.00 d+_aa(i,k) l2_1_abab(k,j,a,b)
    //              += +1.00 <j,l||d,b>_bbbb t2_1_bbbb(d,c,k,l) l2_1_abab(i,k,a,c)
    //              += -1.00 <j,c||d,b>_bbbb t1_1_bb(d,k) l2_1_abab(i,k,a,c)
    //              += +1.00 f_aa(c,a) l2_abab(i,j,c,b)
    //              += +1.00 <k,c||d,b>_bbbb t1_1_bb(d,k) l2_1_abab(i,j,a,c)
    //              += -0.50 <l,k||d,b>_bbbb t2_1_bbbb(d,c,l,k) l2_1_abab(i,j,a,c)
    //              += -1.00 f_aa(i,k) l2_abab(k,j,a,b)
    //              += -1.00 <i,c||d,b>_abab t1_1_aa(d,k) l2_1_abab(k,j,a,c)
    //              += -1.00 <c,j||d,b>_abab t1_1_aa(d,k) l2_1_aaaa(i,k,c,a)
    //              += +0.50 <d,c||a,b>_abab l2_abab(i,j,d,c)
    //              += +0.50 <c,d||a,b>_abab l2_abab(i,j,c,d)
    //              += +1.00 f_bb(c,b) l2_abab(i,j,a,c)
    //              += +1.00 <j,l||c,k>_bbbb t1_1_bb(c,l) l2_1_abab(i,k,a,b)
    //              += -0.50 <j,l||c,d>_bbbb t2_1_bbbb(c,d,k,l) l2_1_abab(i,k,a,b)
    //              += +1.00 <l,j||d,b>_abab t2_1_abab(d,c,l,k) l2_1_abab(i,k,a,c)
    //              += -0.50 <i,l||c,d>_abab t2_abab(c,d,k,l) l2_abab(k,j,a,b)
    //              += -0.50 <i,l||d,c>_abab t2_abab(d,c,k,l) l2_abab(k,j,a,b)
    //              += +1.00 <l,j||d,b>_abab t2_1_aaaa(d,c,k,l) l2_1_aaaa(i,k,c,a)
    //              += +1.00 <i,l||d,b>_abab t2_1_abab(d,c,k,l) l2_1_abab(k,j,a,c)
    //              += +1.00 d-_bb(j,c) t1_1_bb(c,k) l2_abab(i,k,a,b)
    //              += +1.00 <j,l||d,b>_bbbb t2_bbbb(d,c,k,l) l2_abab(i,k,a,c)
    //              += -1.00 d-_bb(k,c) t1_1_bb(c,k) l2_abab(i,j,a,b)
    //              += -1.00 <c,j||a,k>_abab l2_abab(i,k,c,b)
    //              += +1.00 <j,l||d,b>_bbbb t2_abab(c,d,k,l) l2_aaaa(i,k,c,a)
    //              += +1.00 <j,c||b,k>_bbbb l2_abab(i,k,a,c)
    //              += -1.00 <i,l||k,c>_abab t1_1_bb(c,l) l2_1_abab(k,j,a,b)
    //              += -0.50 <i,l||c,d>_abab t2_1_abab(c,d,k,l) l2_1_abab(k,j,a,b)
    //              += -0.50 <i,l||d,c>_abab t2_1_abab(d,c,k,l) l2_1_abab(k,j,a,b)
    //              += -1.00 <i,c||d,a>_aaaa t1_1_aa(d,k) l2_1_abab(k,j,c,b)
    //              += +1.00 <i,l||d,a>_aaaa t2_1_aaaa(d,c,k,l) l2_1_abab(k,j,c,b)
    //              += +1.00 d-_aa(i,c) t1_1_aa(c,k) l2_abab(k,j,a,b)
    //              += +1.00 d-_aa(k,a) t0_1 t1_1_aa(c,k) l2_1_abab(i,j,c,b)
    //              += -1.00 <i,l||a,k>_aaaa t1_1_aa(c,l) l2_1_abab(k,j,c,b)
    //              += +0.50 <i,j||a,d>_abab t2_1_bbbb(d,c,k,l) l2_1_bbbb(k,l,c,b)
    //              += +0.50 <i,j||a,d>_abab t2_bbbb(d,c,k,l) l2_bbbb(k,l,c,b)
    //              += -0.50 <i,j||d,b>_abab t2_abab(d,c,k,l) l2_abab(k,l,a,c)
    //              += -0.50 <i,j||d,b>_abab t2_abab(d,c,l,k) l2_abab(l,k,a,c)
    //              += -0.50 <i,j||d,b>_abab t2_1_abab(d,c,k,l) l2_1_abab(k,l,a,c)
    //              += -0.50 <i,j||d,b>_abab t2_1_abab(d,c,l,k) l2_1_abab(l,k,a,c)
    //              += +1.00 <i,l||a,k>_abab t1_1_bb(c,l) l2_1_bbbb(j,k,c,b)
    //              += -1.00 f_aa(k,a) t1_1_aa(c,k) l2_1_abab(i,j,c,b)
    //              += +0.50 <i,j||d,b>_abab t2_aaaa(d,c,k,l) l2_aaaa(k,l,c,a)
    //              += +0.50 <i,j||d,b>_abab t2_1_aaaa(d,c,k,l) l2_1_aaaa(k,l,c,a)
    //              += -1.00 <k,d||a,b>_abab t1_1_aa(c,k) l2_1_abab(i,j,c,d)
    //              += +1.00 <i,l||k,b>_abab t1_1_bb(c,l) l2_1_abab(k,j,a,c)
    //              += -0.50 <l,j||a,b>_abab t2_1_abab(d,c,l,k) l2_1_abab(i,k,d,c)
    //              += -0.50 <l,j||a,b>_abab t2_1_abab(c,d,l,k) l2_1_abab(i,k,c,d)
    //              += +0.50 <l,j||a,b>_abab t2_1_aaaa(d,c,k,l) l2_1_aaaa(i,k,d,c)
    //              += +1.00 <l,j||k,b>_abab t1_1_aa(c,l) l2_1_aaaa(i,k,c,a)
    //              += -1.00 <k,j||a,b>_abab t1_1_aa(c,k) l1_1_aa(i,c)
    //              += +1.00 d-_bb(j,b) t1_1_aa(c,k) l2_aaaa(i,k,c,a)
    //              += +1.00 <l,j||a,k>_abab t1_1_aa(c,l) l2_1_abab(i,k,c,b)
    //              += -1.00 d-_bb(j,b) t1_1_bb(c,k) l2_abab(i,k,a,c)
    //              += -0.50 <l,j||a,b>_abab t2_abab(d,c,l,k) l2_abab(i,k,d,c)
    //              += -0.50 <l,j||a,b>_abab t2_abab(c,d,l,k) l2_abab(i,k,c,d)
    //              += +0.50 <l,j||a,b>_abab t2_aaaa(d,c,k,l) l2_aaaa(i,k,d,c)
    sigmal2_abab("L,a,b,i,j")  = -1.00 * tmps_["302_bbaa_Lovov"]("L,j,b,i,a");
    tmps_["302_bbaa_Lovov"].~TArrayD();

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["306_bbbb_Lvovo"]("R,a,i,b,j")  = -1.00 * t1_1["bb"]("a,i") * tmps_["107_bb_Lvo"]("R,b,j");
    tmps_["107_bb_Lvo"].~TArrayD();

    // flops: o3v1L1  = o4v2L1 o4v2L1 o3v1L1
    //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1
    tmps_["305_bbbb_Loovo"]("R,i,j,a,k")  = eri["abab_oovo"]("m,i,c,j") * r2["abab"]("R,c,a,m,k");
    tmps_["305_bbbb_Loovo"]("R,i,j,a,k") += eri["bbbb_oovo"]("i,l,b,j") * r2["bbbb"]("R,b,a,k,l");

    // flops: o3v1L1  = o2v2 o3v2L1
    //  mems: o3v1L1  = o2v2 o3v1L1
    tmps_["304_bbbb_Lvooo"]("R,a,i,j,k")  = (reused_["26_bbbb_voov"]("a,i,j,b") + -1.00 * reused_["28_bbbb_voov"]("a,i,j,b")) * r1["bb"]("R,b,k");

    // flops: o3v1L1  = o3v2L1
    //  mems: o3v1L1  = o3v1L1
    tmps_["303_bbbb_Lvooo"]("R,a,i,j,k")  = eri["bbbb_vovo"]("a,i,b,j") * r1["bb"]("R,b,k");

    // flops: o2v2L1  = o1v1L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1
    //  mems: o2v2L1  = o1v1L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j")  = (tmps_["295_bb_Lvo"]("R,a,i") + -1.00 * tmps_["292_bb_Lvo"]("R,a,i")) * t1_1["bb"]("b,j");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") -= reused_["33_bb_vo"]("a,j") * tmps_["108_bb_Lvo"]("R,b,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += t1_1["bb"]("a,k") * tmps_["305_bbbb_Loovo"]("R,k,j,b,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += tmps_["304_bbbb_Lvooo"]("R,a,j,l,i") * t1_1["bb"]("b,l");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += reused_["32_bb_vo"]("a,j") * tmps_["108_bb_Lvo"]("R,b,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += t1_1["bb"]("a,j") * tmps_["294_bb_Lvo"]("R,b,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += dp["bb_vo"]("a,j") * tmps_["108_bb_Lvo"]("R,b,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += reused_["180_bb_vo"]("a,j") * r1_1["bb"]("R,b,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += t1_1["bb"]("a,i") * tmps_["293_bb_Lov"]("R,j,b");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += reused_["1_bbaa_voov"]("a,j,n,f") * r2_1["abab"]("R,f,b,n,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += reused_["105_baab_vovo"]("a,m,f,j") * r2["abab"]("R,f,b,m,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") -= reused_["27_bbaa_voov"]("a,j,n,f") * r2_1["abab"]("R,f,b,n,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") -= reused_["222_bbbb_vooo"]("a,k,j,i") * r1["bb"]("R,b,k");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += reused_["106_bbbb_vovo"]("a,k,d,j") * r2["bbbb"]("R,d,b,i,k");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") -= reused_["212_bbbb_vvvo"]("a,d,b,j") * r1["bb"]("R,d,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += reused_["32_bb_vo"]("a,j") * r1["bb"]("R,b,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") -= reused_["104_bbbb_voov"]("a,j,l,d") * r2["bbbb"]("R,d,b,i,l");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += reused_["84_bb_vo"]("a,j") * r1_1["bb"]("R,b,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += reused_["102_bbbb_voov"]("a,j,l,d") * r2["bbbb"]("R,d,b,i,l");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += 2.00 * reused_["88_bb_vo"]("a,j") * r1_1["bb"]("R,b,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += reused_["247_bbbb_vvvo"]("a,d,b,j") * r1["bb"]("R,d,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") -= r0_1("R") * reused_["86_bbbb_vovo"]("a,j,b,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += eri["baab_vovo"]("a,m,e,j") * r2_1["abab"]("R,e,b,m,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += reused_["103_bbaa_voov"]("a,j,n,f") * r2["abab"]("R,f,b,n,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") -= f["bb_vo"]("a,j") * r1_1["bb"]("R,b,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") -= reused_["248_bbbb_vooo"]("a,j,l,i") * r1["bb"]("R,b,l");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") -= reused_["33_bb_vo"]("a,j") * r1["bb"]("R,b,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += r1_1["bb"]("R,b,l") * reused_["179_bbbb_oovo"]("l,j,a,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") -= eri["bbbb_vovo"]("a,k,c,j") * r2_1["bbbb"]("R,c,b,i,k");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") -= 0.50 * r1["bb"]("R,b,i") * reused_["199_bb_ov"]("j,a");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") -= reused_["176_bbbb_vvvo"]("a,d,b,j") * r1_1["bb"]("R,d,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += reused_["228_bbbb_vooo"]("a,j,l,i") * r1["bb"]("R,b,l");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") -= reused_["28_bbbb_voov"]("a,j,l,d") * r2_1["bbbb"]("R,d,b,i,l");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") -= reused_["101_bbaa_voov"]("a,j,n,f") * r2["abab"]("R,f,b,n,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") -= r1_1["bb"]("R,b,l") * reused_["178_bbbb_oovo"]("l,j,a,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += reused_["26_bbbb_voov"]("a,j,l,d") * r2_1["bbbb"]("R,d,b,i,l");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += reused_["177_bbbb_vvvo"]("a,d,b,j") * r1_1["bb"]("R,d,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += dp["bb_vo"]("a,j") * r1["bb"]("R,b,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += reused_["91_bb_vo"]("a,j") * r1_1["bb"]("R,b,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") += reused_["200_bb_vo"]("a,j") * r1["bb"]("R,b,i");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") -= tmps_["306_bbbb_Lvovo"]("R,a,i,b,j");
    tmps_["307_bbbb_Lvovo"]("R,a,i,b,j") -= t1_1["bb"]("b,k") * tmps_["303_bbbb_Lvooo"]("R,a,k,j,i");
    tmps_["306_bbbb_Lvovo"].~TArrayD();
    tmps_["305_bbbb_Loovo"].~TArrayD();
    tmps_["304_bbbb_Lvooo"].~TArrayD();
    tmps_["303_bbbb_Lvooo"].~TArrayD();
    tmps_["295_bb_Lvo"].~TArrayD();
    tmps_["294_bb_Lvo"].~TArrayD();
    tmps_["293_bb_Lov"].~TArrayD();
    tmps_["292_bb_Lvo"].~TArrayD();
    tmps_["108_bb_Lvo"].~TArrayD();

    // sigmar2_1_bbbb += +1.00 P(i,j) P(a,b) d-_bb(k,c) r1_bb(c,i) t1_1_bb(a,k) t1_1_bb(b,j)
    //                += -1.00 P(i,j) P(a,b) d-_bb(a,c) r1_1_bb(c,i) t1_1_bb(b,j)
    //                += -1.00 P(i,j) P(a,b) d-_bb(k,c) r0_1 t2_bbbb(c,a,j,k) t1_1_bb(b,i)
    //                += -1.00 P(i,j) P(a,b) <l,k||c,j>_bbbb r2_bbbb(c,b,i,l) t1_1_bb(a,k)
    //                += +1.00 P(i,j) P(a,b) <l,k||c,j>_abab r2_abab(c,b,l,i) t1_1_bb(a,k)
    //                += +1.00 P(i,j) P(a,b) <k,l||c,d>_abab r1_bb(d,i) t2_abab(c,a,k,j) t1_1_bb(b,l)
    //                += +1.00 P(i,j) P(a,b) <l,k||c,d>_bbbb r1_bb(d,i) t2_bbbb(c,a,j,k) t1_1_bb(b,l)
    //                += +1.00 P(i,j) P(a,b) d-_aa(k,c) r0_1 t2_abab(c,a,k,j) t1_1_bb(b,i)
    //                += +1.00 P(i,j) P(a,b) d-_aa(k,c) r2_1_abab(c,b,k,i) t1_1_bb(a,j)
    //                += -1.00 P(i,j) P(a,b) d-_bb(k,c) r2_1_bbbb(c,b,i,k) t1_1_bb(a,j)
    //                += +1.00 P(i,j) P(a,b) d-_bb(a,j) r0_1 t1_1_bb(b,i)
    //                += +1.00 P(i,j) P(a,b) d-_aa(k,c) r1_1_bb(b,i) t2_abab(c,a,k,j) t0_1
    //                += -1.00 P(i,j) P(a,b) d-_bb(k,c) r1_1_bb(b,i) t2_bbbb(c,a,j,k) t0_1
    //                += +1.00 P(i,j) P(a,b) d-_bb(k,j) r1_1_bb(b,k) t1_1_bb(a,i)
    //                += +1.00 P(i,j) P(a,b) <l,k||d,c>_abab r2_1_abab(d,b,l,i) t2_bbbb(c,a,j,k)
    //                += -1.00 P(i,j) P(a,b) <k,a||d,c>_abab r2_abab(d,b,k,i) t1_1_bb(c,j)
    //                += +1.00 P(i,j) P(a,b) <l,k||c,d>_aaaa r2_1_abab(d,b,l,i) t2_abab(c,a,k,j)
    //                += +1.00 P(i,j) P(a,b) <k,a||c,j>_bbbb r1_bb(b,k) t1_1_bb(c,i)
    //                += -1.00 P(i,j) P(a,b) <k,a||c,d>_bbbb r2_bbbb(d,b,i,k) t1_1_bb(c,j)
    //                += +1.00 P(i,j) P(a,b) <k,a||c,d>_abab r1_bb(d,i) t2_1_abab(c,b,k,j)
    //                += +1.00 P(i,j) P(a,b) d+_aa(k,c) r1_bb(b,i) t2_abab(c,a,k,j)
    //                += +1.00 P(i,j) P(a,b) <l,k||c,d>_bbbb r2_bbbb(d,b,i,l) t2_1_bbbb(c,a,j,k)
    //                += +1.00 P(i,j) P(a,b) f_bb(k,c) r1_1_bb(b,i) t2_bbbb(c,a,j,k)
    //                += +1.00 P(i,j) P(a,b) d-_bb(a,j) r1_1_bb(b,i) t0_1
    //                += +0.50 P(i,j) P(a,b) <l,k||c,j>_abab r1_1_bb(b,i) t2_abab(c,a,l,k)
    //                += +0.50 P(i,j) P(a,b) <k,l||c,j>_abab r1_1_bb(b,i) t2_abab(c,a,k,l)
    //                += +0.50 P(i,j) P(a,b) <l,k||c,j>_bbbb r1_1_bb(b,i) t2_bbbb(c,a,l,k)
    //                += -1.00 P(i,j) P(a,b) f_aa(k,c) r1_1_bb(b,i) t2_abab(c,a,k,j)
    //                += -0.50 P(i,j) P(a,b) <k,a||c,d>_abab r1_1_bb(b,i) t2_abab(c,d,k,j)
    //                += -0.50 P(i,j) P(a,b) <k,a||d,c>_abab r1_1_bb(b,i) t2_abab(d,c,k,j)
    //                += +0.50 P(i,j) P(a,b) <k,a||c,d>_bbbb r1_1_bb(b,i) t2_bbbb(c,d,j,k)
    //                += +1.00 P(i,j) P(a,b) <k,l||c,d>_abab r2_bbbb(d,b,i,l) t2_1_abab(c,a,k,j)
    //                += +2.00 P(i,j) P(a,b) d-_bb(a,c) r1_1_bb(b,i) t1_1_bb(c,j)
    //                += -2.00 P(i,j) P(a,b) d-_bb(k,j) r1_1_bb(b,i) t1_1_bb(a,k)
    //                += +2.00 P(i,j) P(a,b) d-_aa(k,c) r1_1_bb(b,i) t2_1_abab(c,a,k,j)
    //                += -2.00 P(i,j) P(a,b) d-_bb(k,c) r1_1_bb(b,i) t2_1_bbbb(c,a,j,k)
    //                += -1.00 P(i,j) P(a,b) <k,a||c,d>_bbbb r1_bb(d,i) t2_1_bbbb(c,b,j,k)
    //                += +1.00 P(i,j) P(a,b) <k,a||c,j>_bbbb r0_1 t2_bbbb(c,b,i,k)
    //                += -1.00 P(i,j) P(a,b) <k,a||c,j>_abab r0_1 t2_abab(c,b,k,i)
    //                += -1.00 P(i,j) P(a,b) <k,a||c,j>_abab r2_1_abab(c,b,k,i)
    //                += +1.00 P(i,j) P(a,b) <l,k||d,c>_abab r2_abab(d,b,l,i) t2_1_bbbb(c,a,j,k)
    //                += -1.00 P(i,j) P(a,b) f_bb(a,j) r1_1_bb(b,i)
    //                += +1.00 P(i,j) P(a,b) <l,k||c,d>_bbbb r1_bb(b,l) t2_bbbb(c,a,j,k) t1_1_bb(d,i)
    //                += -1.00 P(i,j) P(a,b) <l,k||c,j>_bbbb r1_bb(b,l) t2_1_bbbb(c,a,i,k)
    //                += -1.00 P(i,j) P(a,b) d+_bb(k,c) r1_bb(b,i) t2_bbbb(c,a,j,k)
    //                += -1.00 P(i,j) P(a,b) <l,k||c,j>_bbbb r1_1_bb(b,l) t2_bbbb(c,a,i,k)
    //                += +1.00 P(i,j) P(a,b) <k,a||c,j>_bbbb r2_1_bbbb(c,b,i,k)
    //                += +0.50 P(i,j) P(a,b) <l,k||c,d>_bbbb r1_bb(b,i) t2_bbbb(c,d,j,k) t1_1_bb(a,l)
    //                += +1.00 P(i,j) P(a,b) <l,k||d,c>_abab r1_bb(b,i) t2_bbbb(c,a,j,k) t1_1_aa(d,l)
    //                += +0.50 P(i,j) P(a,b) <l,k||c,d>_abab r1_bb(b,i) t2_abab(c,a,l,k) t1_1_bb(d,j)
    //                += +0.50 P(i,j) P(a,b) <k,l||c,d>_abab r1_bb(b,i) t2_abab(c,a,k,l) t1_1_bb(d,j)
    //                += -1.00 P(i,j) P(a,b) f_aa(k,c) r1_bb(b,i) t2_1_abab(c,a,k,j)
    //                += +1.00 P(i,j) P(a,b) f_bb(k,j) r1_bb(b,i) t1_1_bb(a,k)
    //                += +0.50 P(i,j) P(a,b) <k,a||c,d>_bbbb r1_bb(b,i) t2_1_bbbb(c,d,j,k)
    //                += -1.00 P(i,j) P(a,b) <k,a||c,j>_abab r1_bb(b,i) t1_1_aa(c,k)
    //                += +1.00 P(i,j) P(a,b) f_bb(k,c) r1_bb(b,i) t2_1_bbbb(c,a,j,k)
    //                += -1.00 P(i,j) P(a,b) f_bb(a,c) r1_bb(b,i) t1_1_bb(c,j)
    //                += -0.50 P(i,j) P(a,b) <k,a||c,d>_abab r1_bb(b,i) t2_1_abab(c,d,k,j)
    //                += -0.50 P(i,j) P(a,b) <k,a||d,c>_abab r1_bb(b,i) t2_1_abab(d,c,k,j)
    //                += -1.00 P(i,j) P(a,b) <k,a||c,j>_bbbb r1_bb(b,i) t1_1_bb(c,k)
    //                += +0.50 P(i,j) P(a,b) <l,k||c,j>_bbbb r1_bb(b,i) t2_1_bbbb(c,a,l,k)
    //                += +1.00 P(i,j) P(a,b) d-_bb(k,c) r1_bb(b,i) t1_1_bb(a,j) t1_1_bb(c,k)
    //                += +1.00 P(i,j) P(a,b) d-_aa(k,c) r1_bb(b,i) t1_1_bb(a,j) t1_1_aa(c,k)
    //                += -1.00 P(i,j) P(a,b) r1_bb(b,i) t1_1_bb(a,j) w0
    //                += +0.50 P(i,j) P(a,b) <l,k||c,j>_abab r1_bb(b,i) t2_1_abab(c,a,l,k)
    //                += +0.50 P(i,j) P(a,b) <k,l||c,j>_abab r1_bb(b,i) t2_1_abab(c,a,k,l)
    //                += +0.50 P(i,j) P(a,b) <k,l||c,d>_abab r1_bb(b,i) t2_abab(c,d,k,j) t1_1_bb(a,l)
    //                += +0.50 P(i,j) P(a,b) <k,l||d,c>_abab r1_bb(b,i) t2_abab(d,c,k,j) t1_1_bb(a,l)
    //                += -2.00 P(i,j) P(a,b) d-_bb(k,c) r1_bb(b,i) t1_1_bb(a,k) t1_1_bb(c,j)
    //                += +1.00 P(i,j) P(a,b) d-_bb(a,c) r1_bb(b,i) t0_1 t1_1_bb(c,j)
    //                += -1.00 P(i,j) P(a,b) d-_bb(k,j) r1_bb(b,i) t0_1 t1_1_bb(a,k)
    //                += +1.00 P(i,j) P(a,b) d-_aa(k,c) r1_bb(b,i) t0_1 t2_1_abab(c,a,k,j)
    //                += -1.00 P(i,j) P(a,b) d-_bb(k,c) r1_bb(b,i) t0_1 t2_1_bbbb(c,a,j,k)
    //                += +0.50 P(i,j) P(a,b) <l,k||c,d>_bbbb r1_bb(b,i) t2_bbbb(c,a,l,k) t1_1_bb(d,j)
    //                += -1.00 P(i,j) P(a,b) <k,l||c,d>_abab r1_bb(b,i) t2_abab(c,a,k,j) t1_1_bb(d,l)
    //                += +1.00 P(i,j) P(a,b) <k,a||c,d>_abab r1_1_bb(d,i) t2_abab(c,b,k,j)
    //                += +1.00 P(i,j) P(a,b) <k,l||c,d>_abab r1_bb(b,l) t2_abab(c,a,k,j) t1_1_bb(d,i)
    //                += -1.00 P(i,j) P(a,b) <k,l||c,j>_abab r1_bb(b,l) t2_1_abab(c,a,k,i)
    //                += +1.00 P(i,j) P(a,b) <l,k||c,d>_bbbb r2_1_bbbb(d,b,i,l) t2_bbbb(c,a,j,k)
    //                += +1.00 P(i,j) P(a,b) <l,k||c,d>_aaaa r2_abab(d,b,l,i) t2_1_abab(c,a,k,j)
    //                += -1.00 P(i,j) P(a,b) <k,l||c,j>_abab r1_1_bb(b,l) t2_abab(c,a,k,i)
    //                += +1.00 P(i,j) P(a,b) <k,l||c,d>_abab r2_1_bbbb(d,b,i,l) t2_abab(c,a,k,j)
    //                += -1.00 P(i,j) P(a,b) <k,a||c,d>_bbbb r1_1_bb(d,i) t2_bbbb(c,b,j,k)
    //                += +1.00 P(i,j) P(a,b) d+_bb(a,j) r1_bb(b,i)
    //                += +1.00 P(i,j) P(a,b) d-_aa(k,k) r1_1_bb(b,i) t1_1_bb(a,j)
    //                += +1.00 P(i,j) P(a,b) d-_bb(k,k) r1_1_bb(b,i) t1_1_bb(a,j)
    //                += -1.00 P(i,j) P(a,b) <l,k||c,d>_bbbb r1_bb(b,i) t2_bbbb(c,a,j,k) t1_1_bb(d,l)
    //                += +1.00 P(i,j) P(a,b) <l,k||c,d>_aaaa r1_bb(b,i) t2_abab(c,a,k,j) t1_1_aa(d,l)
    //                += +1.00 P(i,j) P(a,b) d-_bb(k,c) r1_bb(b,k) t1_1_bb(a,i) t1_1_bb(c,j)
    //                += +1.00 P(i,j) P(a,b) <k,a||c,j>_bbbb r1_bb(c,i) t1_1_bb(b,k)
    sigmar2_1_bbbb("R,a,b,i,j") += tmps_["307_bbbb_Lvovo"]("R,a,i,b,j");
    sigmar2_1_bbbb("R,a,b,i,j") -= tmps_["307_bbbb_Lvovo"]("R,a,j,b,i");
    sigmar2_1_bbbb("R,a,b,i,j") -= tmps_["307_bbbb_Lvovo"]("R,b,i,a,j");
    sigmar2_1_bbbb("R,a,b,i,j") += tmps_["307_bbbb_Lvovo"]("R,b,j,a,i");
    tmps_["307_bbbb_Lvovo"].~TArrayD();

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["310_aabb_Lovov"]("L,i,a,j,b")  = -1.00 * dp["aa_ov"]("i,a") * tmps_["263_bb_Lov"]("L,j,b");
    tmps_["263_bb_Lov"].~TArrayD();

    // flops: o2v2L1  = o3v3L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["309_baab_Lvoov"]("L,a,i,j,b")  = t2["abab"]("c,a,i,k") * l2_1["abab"]("L,j,k,c,b");

    // flops: o4v0L1  = o4v2L1
    //  mems: o4v0L1  = o4v0L1
    tmps_["308_abab_Loooo"]("L,i,j,k,l")  = l2_1["abab"]("L,i,j,a,b") * t2["abab"]("a,b,k,l");

    // flops: o2v2L1  = o3v2L1 o2v3L1 o3v2L1 o3v3L1 o4v2L1 o2v4L1 o2v3L1 o3v2L1 o3v2L1 o3v2L1 o3v3L1 o3v2L1 o2v3L1 o3v2L1 o2v2L1 o3v3L1 o2v3L1 o2v3L1 o3v3L1 o2v2L1 o2v3L1 o4v2L1 o3v2L1 o3v3L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o4v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b")  = -0.50 * tmps_["201_aa_Loo"]("L,j,k") * eri["abab_oovv"]("k,i,b,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= 0.50 * tmps_["186_aa_Lvv"]("L,b,c") * eri["abab_oovv"]("j,i,c,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= 0.50 * eri["abab_oovv"]("j,m,b,a") * tmps_["179_bb_Loo"]("L,i,m");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= eri["abab_oovv"]("j,m,b,e") * tmps_["261_bbbb_Lvoov"]("L,e,m,i,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= eri["abab_oovv"]("k,n,b,a") * tmps_["308_abab_Loooo"]("L,j,i,k,n");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= eri["abab_vvvv"]("c,f,b,a") * l2_1["abab"]("L,j,i,c,f");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= f["aa_vv"]("d,b") * l2_1["abab"]("L,j,i,d,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= l2["abab"]("L,j,n,b,a") * dp["bb_oo"]("i,n");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += l1_1["bb"]("L,n,a") * eri["abab_oovo"]("j,i,b,n");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += 0.50 * l2_1["abab"]("L,j,n,b,a") * reused_["161_bb_oo"]("i,n");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= l2_1["aaaa"]("L,j,l,d,b") * eri["abba_vovo"]("d,i,a,l");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= dp["aa_oo"]("j,l") * l2["abab"]("L,l,i,b,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,j,i,d,a") * reused_["59_aa_vv"]("d,b");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,j,n,b,a") * f["bb_oo"]("i,n");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= l1_1["aa"]("L,j,b") * f["bb_ov"]("i,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= l2_1["bbbb"]("L,i,n,f,a") * reused_["1_bbaa_voov"]("f,n,j,b");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += 0.50 * l2_1["abab"]("L,j,i,b,f") * reused_["109_bb_vv"]("a,f");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += eri["baab_vovv"]("f,j,b,a") * l1_1["bb"]("L,i,f");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,j,n,b,f") * reused_["163_bbbb_voov"]("f,n,i,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= eri["abab_oovv"]("j,i,b,a") * l0_1("L");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= l1_1["aa"]("L,j,d") * eri["abab_vovv"]("d,i,b,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= l2_1["abab"]("L,l,m,b,a") * eri["abab_oooo"]("j,i,l,m");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= l1_1["aa"]("L,l,b") * eri["abba_oovo"]("j,i,a,l");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= eri["baab_vovo"]("f,j,b,n") * l2_1["bbbb"]("L,i,n,f,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= l2_1["abab"]("L,j,n,b,a") * reused_["35_bb_oo"]("i,n");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += l1_1["aa"]("L,j,b") * reused_["12_bb_ov"]("i,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,j,i,d,a") * reused_["57_aa_vv"]("b,d");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,j,i,b,f") * reused_["34_bb_vv"]("f,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += eri["baba_vovo"]("f,j,a,l") * l2_1["abab"]("L,l,i,b,f");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= l2_1["abab"]("L,j,n,b,f") * reused_["26_bbbb_voov"]("f,n,i,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= l2_1["aaaa"]("L,j,l,d,b") * reused_["6_aabb_voov"]("d,l,i,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += f["aa_oo"]("j,l") * l2_1["abab"]("L,l,i,b,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += scalars_["1"] * l2["abab"]("L,j,i,b,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= w0 * l2_1["abab"]("L,j,i,b,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,j,i,b,f") * reused_["3_bb_vv"]("a,f");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += 0.50 * l2_1["abab"]("L,l,i,b,a") * reused_["174_aa_oo"]("j,l");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,j,n,b,f") * eri["bbbb_vovo"]("f,i,a,n");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,l,i,d,a") * reused_["69_aaaa_voov"]("d,l,j,b");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += scalars_["5"] * l2_1["abab"]("L,j,i,b,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += l1["aa"]("L,j,b") * dp["bb_ov"]("i,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,j,n,b,a") * reused_["30_bb_oo"]("i,n");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += 2.00 * scalars_["3"] * l2_1["abab"]("L,j,i,b,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= l2_1["abab"]("L,l,i,b,f") * reused_["78_baab_voov"]("f,l,j,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,l,i,b,a") * reused_["16_aa_oo"]("j,l");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += dp["bb_vv"]("f,a") * l2["abab"]("L,j,i,b,f");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += eri["aaaa_vovo"]("d,j,b,l") * l2_1["abab"]("L,l,i,d,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= l2_1["abab"]("L,l,i,b,a") * reused_["17_aa_oo"]("j,l");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += scalars_["2"] * l2["abab"]("L,j,i,b,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += 0.50 * l2_1["abab"]("L,j,i,d,a") * reused_["18_aa_vv"]("b,d");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += 2.00 * scalars_["4"] * l2_1["abab"]("L,j,i,b,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += dp["aa_ov"]("j,b") * l1["bb"]("L,i,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= f["aa_ov"]("j,b") * l1_1["bb"]("L,i,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += reused_["13_aa_ov"]("j,b") * l1_1["bb"]("L,i,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += l2_1["aaaa"]("L,j,l,d,b") * reused_["70_aabb_voov"]("d,l,i,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= f["bb_vv"]("f,a") * l2_1["abab"]("L,j,i,b,f");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += l2_1["bbbb"]("L,i,n,f,a") * reused_["158_bbaa_voov"]("f,n,j,b");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= l2_1["abab"]("L,l,m,b,a") * reused_["80_abab_oooo"]("j,i,l,m");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += dp["aa_vv"]("d,b") * l2["abab"]("L,j,i,d,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= 2.00 * dp["aa_ov"]("l,b") * tmps_["181_abba_Loovo"]("L,j,i,a,l");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= 2.00 * dp["bb_ov"]("n,a") * tmps_["169_abab_Loovo"]("L,j,i,b,n");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += eri["abab_oovv"]("j,m,b,a") * tmps_["178_bb_Loo"]("L,i,m");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += dp["aa_ov"]("j,b") * tmps_["101_bb_Lov"]("L,i,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= eri["abab_oovv"]("k,i,b,e") * tmps_["309_baab_Lvoov"]("L,e,k,j,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += tmps_["310_aabb_Lovov"]("L,j,b,i,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += tmps_["99_aa_Lvv"]("L,b,c") * eri["abab_oovv"]("j,i,c,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += tmps_["102_aa_Lov"]("L,j,b") * dp["bb_ov"]("i,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= 0.50 * tmps_["238_bb_Lvv"]("L,a,e") * eri["abab_oovv"]("j,i,b,e");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += tmps_["205_aa_Loo"]("L,j,k") * eri["abab_oovv"]("k,i,b,a");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") += tmps_["199_bb_Lvv"]("L,a,e") * eri["abab_oovv"]("j,i,b,e");
    tmps_["311_bbaa_Lovov"]("L,i,a,j,b") -= tmps_["173_aa_Lov"]("L,j,b") * dp["bb_ov"]("i,a");
    tmps_["310_aabb_Lovov"].~TArrayD();
    tmps_["309_baab_Lvoov"].~TArrayD();
    tmps_["308_abab_Loooo"].~TArrayD();
    tmps_["261_bbbb_Lvoov"].~TArrayD();
    tmps_["238_bb_Lvv"].~TArrayD();
    tmps_["199_bb_Lvv"].~TArrayD();
    tmps_["186_aa_Lvv"].~TArrayD();
    tmps_["181_abba_Loovo"].~TArrayD();
    tmps_["179_bb_Loo"].~TArrayD();
    tmps_["178_bb_Loo"].~TArrayD();
    tmps_["173_aa_Lov"].~TArrayD();
    tmps_["169_abab_Loovo"].~TArrayD();
    tmps_["102_aa_Lov"].~TArrayD();
    tmps_["101_bb_Lov"].~TArrayD();
    tmps_["99_aa_Lvv"].~TArrayD();

    // sigmal2_1_abab += -1.00 d-_bb(j,b) t1_1_bb(c,k) l2_1_abab(i,k,a,c)
    //                += -0.50 <i,j||d,b>_abab t2_abab(d,c,k,l) l2_1_abab(k,l,a,c)
    //                += -0.50 <i,j||d,b>_abab t2_abab(d,c,l,k) l2_1_abab(l,k,a,c)
    //                += +0.50 <i,j||a,d>_abab t2_bbbb(d,c,k,l) l2_1_bbbb(k,l,c,b)
    //                += -0.50 <l,j||a,b>_abab t2_abab(d,c,l,k) l2_1_abab(i,k,d,c)
    //                += -0.50 <l,j||a,b>_abab t2_abab(c,d,l,k) l2_1_abab(i,k,c,d)
    //                += -0.50 <i,j||a,d>_abab t2_abab(c,d,k,l) l2_1_abab(k,l,c,b)
    //                += -0.50 <i,j||a,d>_abab t2_abab(c,d,l,k) l2_1_abab(l,k,c,b)
    //                += -1.00 d-_bb(j,b) t0_1 l1_1_aa(i,a)
    //                += -0.50 <l,k||a,d>_abab t2_abab(c,d,l,k) l2_1_abab(i,j,c,b)
    //                += -0.50 <k,l||a,d>_abab t2_abab(c,d,k,l) l2_1_abab(i,j,c,b)
    //                += +1.00 d-_bb(j,k) t0_1 l2_1_abab(i,k,a,b)
    //                += -1.00 d-_bb(c,b) t0_1 l2_1_abab(i,j,a,c)
    //                += -1.00 <i,c||a,k>_abab l2_1_bbbb(j,k,c,b)
    //                += -1.00 <i,c||k,b>_abab l2_1_abab(k,j,a,c)
    //                += +1.00 <l,j||d,b>_abab t2_abab(d,c,l,k) l2_1_abab(i,k,a,c)
    //                += +1.00 <l,j||d,b>_abab t2_aaaa(d,c,k,l) l2_1_aaaa(i,k,c,a)
    //                += -1.00 f_aa(i,k) l2_1_abab(k,j,a,b)
    //                += -1.00 <i,j||k,b>_abab l1_1_aa(k,a)
    //                += +0.50 <i,j||k,l>_abab l2_1_abab(k,l,a,b)
    //                += +0.50 <i,j||l,k>_abab l2_1_abab(l,k,a,b)
    //                += -1.00 d-_aa(k,k) l2_abab(i,j,a,b)
    //                += +1.00 l2_1_abab(i,j,a,b) w0
    //                += +1.00 <c,j||a,b>_abab l1_1_aa(i,c)
    //                += -0.50 <l,k||d,b>_abab t2_abab(d,c,l,k) l2_1_abab(i,j,a,c)
    //                += -0.50 <k,l||d,b>_abab t2_abab(d,c,k,l) l2_1_abab(i,j,a,c)
    //                += +1.00 <i,j||a,b>_abab l0_1
    //                += -0.50 <i,l||c,d>_aaaa t2_aaaa(c,d,k,l) l2_1_abab(k,j,a,b)
    //                += +1.00 <j,l||d,b>_bbbb t2_bbbb(d,c,k,l) l2_1_abab(i,k,a,c)
    //                += +1.00 <i,c||a,b>_abab l1_1_bb(j,c)
    //                += +1.00 <j,c||b,k>_bbbb l2_1_abab(i,k,a,c)
    //                += -0.50 <l,k||d,b>_bbbb t2_bbbb(d,c,l,k) l2_1_abab(i,j,a,c)
    //                += +1.00 <i,l||a,d>_abab t2_bbbb(d,c,k,l) l2_1_bbbb(j,k,c,b)
    //                += +1.00 f_bb(j,b) l1_1_aa(i,a)
    //                += +1.00 <i,l||d,a>_aaaa t2_aaaa(d,c,k,l) l2_1_abab(k,j,c,b)
    //                += -1.00 f_bb(j,k) l2_1_abab(i,k,a,b)
    //                += -1.00 d-_aa(c,a) t0_1 l2_1_abab(i,j,c,b)
    //                += +1.00 d-_aa(i,k) l2_abab(k,j,a,b)
    //                += -1.00 <c,j||k,b>_abab l2_1_aaaa(i,k,c,a)
    //                += -0.50 <j,l||c,d>_bbbb t2_bbbb(c,d,k,l) l2_1_abab(i,k,a,b)
    //                += -1.00 d-_aa(k,k) t0_1 l2_1_abab(i,j,a,b)
    //                += +1.00 f_aa(k,k) r2_1_aaaa(a,b,i,j)
    //                += -0.50 <l,k||l,k>_abab r2_1_aaaa(a,b,i,j)
    //                += -0.50 <k,l||k,l>_abab r2_1_aaaa(a,b,i,j)
    //                += +0.25 <l,k||c,d>_aaaa r2_1_aaaa(a,b,i,j) t2_aaaa(c,d,l,k)
    //                += -0.50 <l,k||l,k>_bbbb r2_1_aaaa(a,b,i,j)
    //                += -0.50 <l,k||l,k>_aaaa r2_1_aaaa(a,b,i,j)
    //                += +0.25 <l,k||c,d>_bbbb r2_1_aaaa(a,b,i,j) t2_bbbb(c,d,l,k)
    //                += +1.00 f_bb(k,k) r2_1_aaaa(a,b,i,j)
    //                += -1.00 d-_bb(k,k) r2_1_aaaa(a,b,i,j) t0_1
    //                += +0.25 <l,k||c,d>_abab r2_1_aaaa(a,b,i,j) t2_abab(c,d,l,k)
    //                += +0.25 <k,l||c,d>_abab r2_1_aaaa(a,b,i,j) t2_abab(c,d,k,l)
    //                += +0.25 <l,k||d,c>_abab r2_1_aaaa(a,b,i,j) t2_abab(d,c,l,k)
    //                += +0.25 <k,l||d,c>_abab r2_1_aaaa(a,b,i,j) t2_abab(d,c,k,l)
    //                += -1.00 <i,j||a,k>_abab l1_1_bb(k,b)
    //                += +1.00 d-_bb(j,k) l2_abab(i,k,a,b)
    //                += -1.00 d-_bb(j,b) l1_aa(i,a)
    //                += -0.50 <l,j||c,d>_abab t2_abab(c,d,l,k) l2_1_abab(i,k,a,b)
    //                += -0.50 <l,j||d,c>_abab t2_abab(d,c,l,k) l2_1_abab(i,k,a,b)
    //                += -2.00 d-_bb(k,c) t1_1_bb(c,k) l2_1_abab(i,j,a,b)
    //                += +1.00 <i,l||d,b>_abab t2_abab(d,c,k,l) l2_1_abab(k,j,a,c)
    //                += -0.50 <i,l||c,d>_abab t2_abab(c,d,k,l) l2_1_abab(k,j,a,b)
    //                += -0.50 <i,l||d,c>_abab t2_abab(d,c,k,l) l2_1_abab(k,j,a,b)
    //                += +1.00 f_aa(c,a) l2_1_abab(i,j,c,b)
    //                += -1.00 d-_bb(c,b) l2_abab(i,j,a,c)
    //                += +1.00 <i,c||a,k>_aaaa l2_1_abab(k,j,c,b)
    //                += +1.00 d-_aa(i,k) t0_1 l2_1_abab(k,j,a,b)
    //                += -1.00 d-_bb(k,k) l2_abab(i,j,a,b)
    //                += -0.50 <l,k||d,a>_aaaa t2_aaaa(d,c,l,k) l2_1_abab(i,j,c,b)
    //                += -2.00 d-_aa(k,c) t1_1_aa(c,k) l2_1_abab(i,j,a,b)
    //                += -1.00 d-_aa(i,a) l1_bb(j,b)
    //                += +1.00 f_aa(i,a) l1_1_bb(j,b)
    //                += -1.00 d-_aa(i,a) t0_1 l1_1_bb(j,b)
    //                += +1.00 <j,l||d,b>_bbbb t2_abab(c,d,k,l) l2_1_aaaa(i,k,c,a)
    //                += +1.00 f_bb(c,b) l2_1_abab(i,j,a,c)
    //                += +0.50 <d,c||a,b>_abab l2_1_abab(i,j,d,c)
    //                += +0.50 <c,d||a,b>_abab l2_1_abab(i,j,c,d)
    //                += +1.00 <i,l||d,a>_aaaa t2_abab(d,c,l,k) l2_1_bbbb(j,k,c,b)
    //                += +0.25 <i,j||c,d>_abab t2_abab(c,d,k,l) l2_1_abab(k,l,a,b)
    //                += +0.25 <i,j||c,d>_abab t2_abab(c,d,l,k) l2_1_abab(l,k,a,b)
    //                += +0.25 <i,j||d,c>_abab t2_abab(d,c,k,l) l2_1_abab(k,l,a,b)
    //                += +0.25 <i,j||d,c>_abab t2_abab(d,c,l,k) l2_1_abab(l,k,a,b)
    //                += -1.00 d-_aa(c,a) l2_abab(i,j,c,b)
    //                += +2.00 d-_aa(k,a) t1_1_aa(c,k) l2_1_abab(i,j,c,b)
    //                += +2.00 d-_bb(k,b) t1_1_bb(c,k) l2_1_abab(i,j,a,c)
    //                += -0.50 <i,l||a,b>_abab t2_abab(d,c,k,l) l2_1_abab(k,j,d,c)
    //                += -0.50 <i,l||a,b>_abab t2_abab(c,d,k,l) l2_1_abab(k,j,c,d)
    //                += +0.25 <l,k||a,b>_abab t2_abab(d,c,l,k) l2_1_abab(i,j,d,c)
    //                += +0.25 <k,l||a,b>_abab t2_abab(d,c,k,l) l2_1_abab(i,j,d,c)
    //                += +0.25 <l,k||a,b>_abab t2_abab(c,d,l,k) l2_1_abab(i,j,c,d)
    //                += +0.25 <k,l||a,b>_abab t2_abab(c,d,k,l) l2_1_abab(i,j,c,d)
    //                += -1.00 d-_aa(i,a) t1_1_aa(c,k) l2_1_abab(k,j,c,b)
    //                += +1.00 <i,l||a,d>_abab t2_abab(c,d,k,l) l2_1_abab(k,j,c,b)
    //                += +1.00 <l,j||a,d>_abab t2_abab(c,d,l,k) l2_1_abab(i,k,c,b)
    //                += +0.50 <i,l||a,b>_abab t2_bbbb(d,c,k,l) l2_1_bbbb(j,k,d,c)
    //                += +1.00 d-_aa(i,a) t1_1_bb(c,k) l2_1_bbbb(j,k,c,b)
    //                += +1.00 d-_bb(j,b) t1_1_aa(c,k) l2_1_aaaa(i,k,c,a)
    //                += +0.50 <i,j||d,b>_abab t2_aaaa(d,c,k,l) l2_1_aaaa(k,l,c,a)
    //                += +0.50 <l,j||a,b>_abab t2_aaaa(d,c,k,l) l2_1_aaaa(i,k,d,c)
    sigmal2_1_abab("L,a,b,i,j") -= tmps_["311_bbaa_Lovov"]("L,j,b,i,a");
    tmps_["311_bbaa_Lovov"].~TArrayD();

    // flops: o2v2L1  = o2v0L1 o3v2L1
    //  mems: o2v2L1  = o2v0L1 o2v2L1
    tmps_["312_aaaa_Loovv"]("L,i,j,a,b")  = (tmps_["205_aa_Loo"]("L,i,k") + -0.50 * tmps_["201_aa_Loo"]("L,i,k")) * eri["aaaa_oovv"]("j,k,a,b");
    tmps_["312_aaaa_Loovv"].~TArrayD();

    // flops: o2v2L1  = o2v0L1 o3v2L1 o2v0 o3v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v0L1 o2v2L1 o2v0 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["313_aaaa_Loovv"]("L,i,j,a,b")  = (tmps_["205_aa_Loo"]("L,i,k") + -0.50 * tmps_["201_aa_Loo"]("L,i,k")) * eri["aaaa_oovv"]("j,k,a,b");
    tmps_["313_aaaa_Loovv"]("L,i,j,a,b") -= (reused_["16_aa_oo"]("j,m") + 0.50 * reused_["174_aa_oo"]("j,m")) * l2_1["aaaa"]("L,i,m,a,b");
    tmps_["313_aaaa_Loovv"]("L,i,j,a,b") += eri["aaaa_vovv"]("e,j,a,b") * l1_1["aa"]("L,i,e");
    tmps_["313_aaaa_Loovv"]("L,i,j,a,b") += dp["aa_oo"]("j,m") * l2["aaaa"]("L,i,m,a,b");
    tmps_["313_aaaa_Loovv"]("L,i,j,a,b") += l2_1["aaaa"]("L,i,m,a,b") * reused_["17_aa_oo"]("j,m");
    tmps_["313_aaaa_Loovv"]("L,i,j,a,b") -= f["aa_oo"]("j,m") * l2_1["aaaa"]("L,i,m,a,b");
    tmps_["205_aa_Loo"].~TArrayD();
    tmps_["201_aa_Loo"].~TArrayD();

    // sigmal2_1_aaaa += +0.50 P(i,j) <j,l||a,b>_aaaa t2_abab(d,c,l,k) l2_1_abab(i,k,d,c)
    //                += +0.50 P(i,j) <j,l||a,b>_aaaa t2_abab(c,d,l,k) l2_1_abab(i,k,c,d)
    //                += -0.50 P(i,j) <j,l||a,b>_aaaa t2_aaaa(d,c,k,l) l2_1_aaaa(i,k,d,c)
    //                += -0.50 P(i,j) <j,l||c,d>_abab t2_abab(c,d,k,l) l2_1_aaaa(i,k,a,b)
    //                += -0.50 P(i,j) <j,l||d,c>_abab t2_abab(d,c,k,l) l2_1_aaaa(i,k,a,b)
    //                += -0.50 P(i,j) <j,l||c,d>_aaaa t2_aaaa(c,d,k,l) l2_1_aaaa(i,k,a,b)
    //                += -1.00 P(i,j) <j,c||a,b>_aaaa l1_1_aa(i,c)
    //                += +1.00 P(i,j) d-_aa(j,k) l2_aaaa(i,k,a,b)
    //                += +1.00 P(i,j) d-_aa(j,k) t0_1 l2_1_aaaa(i,k,a,b)
    //                += -1.00 P(i,j) f_aa(j,k) l2_1_aaaa(i,k,a,b)
    sigmal2_1_aaaa("L,a,b,i,j") += tmps_["313_aaaa_Loovv"]("L,i,j,a,b");
    sigmal2_1_aaaa("L,a,b,i,j") -= tmps_["313_aaaa_Loovv"]("L,j,i,a,b");
    tmps_["313_aaaa_Loovv"].~TArrayD();
}
