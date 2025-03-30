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

void hilbert::EOM_EE_QED_CCSD_21::sigma_ee_21_4() {

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


    // flops: o1v1L1  = o2v2L1 o2v2L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1
    tmps_["162_aa_Lvo"]("R,a,i")  = -1.00 * dp["aa_ov"]("k,c") * r2_1["aaaa"]("R,c,a,i,k");
    tmps_["162_aa_Lvo"]("R,a,i") += dp["bb_ov"]("j,b") * r2_1["abab"]("R,a,b,i,j");

    // sigmar1_aa += -1.00 d-_bb(j,b) r2_1_abab(a,b,i,j)
    //            += +1.00 d-_aa(j,b) r2_1_aaaa(b,a,i,j)
    sigmar1_aa("R,a,i") -= tmps_["162_aa_Lvo"]("R,a,i");

    // flops: o1v1L1  = o2v1L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["165_aa_Lov"]("R,i,a")  = dp["aa_oo"]("j,i") * r1_1["aa"]("R,a,j");

    // sigmar1_aa += +1.00 d-_aa(j,i) r1_1_aa(a,j)
    sigmar1_aa("R,a,i") += tmps_["165_aa_Lov"]("R,i,a");

    // flops: o1v1L1  = o1v2L1 o2v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1
    tmps_["166_aa_Lvo"]("R,a,i")  = -1.00 * dp["aa_vv"]("a,b") * r1_1["aa"]("R,b,i");
    tmps_["166_aa_Lvo"]("R,a,i") += t1_1["aa"]("a,j") * tmps_["35_aa_Loo"]("R,i,j");

    // sigmar1_aa += +1.00 d-_aa(j,b) r1_aa(b,i) t1_1_aa(a,j)
    //            += -1.00 d-_aa(a,b) r1_1_aa(b,i)
    sigmar1_aa("R,a,i") += tmps_["166_aa_Lvo"]("R,a,i");

    // flops: o2v2L1  = o3v3L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["164_aabb_Lvoov"]("R,a,i,j,b")  = r2["aaaa"]("R,c,a,i,k") * eri["abab_oovv"]("k,j,c,b");

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["163_aaaa_Lvovo"]("R,a,i,b,j")  = -1.00 * t1_1["aa"]("a,i") * tmps_["162_aa_Lvo"]("R,b,j");

    // flops: o3v1L1  = o4v2L1 o4v2L1 o3v1L1
    //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1
    tmps_["161_aaaa_Loovo"]("R,i,j,a,k")  = r2["aaaa"]("R,b,a,k,l") * eri["aaaa_oovo"]("i,l,b,j");
    tmps_["161_aaaa_Loovo"]("R,i,j,a,k") -= eri["abba_oovo"]("i,m,c,j") * r2["abab"]("R,a,c,k,m");

    // flops: o3v1L1  = o2v2 o3v2L1
    //  mems: o3v1L1  = o2v2 o3v1L1
    tmps_["160_aaaa_Lvooo"]("R,a,i,j,k")  = (reused_["143_aaaa_voov"]("a,i,j,b") + -1.00 * reused_["55_aaaa_voov"]("a,i,j,b")) * r1["aa"]("R,b,k");

    // flops: o2v2L1  = o3v3L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["159_aabb_Lvoov"]("R,a,i,j,b")  = r2_1["aaaa"]("R,c,a,i,k") * eri["abab_oovv"]("k,j,c,b");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["156_aa_Lvo"]("R,a,i")  = -1.00 * dp["aa_ov"]("j,b") * r2_1["aaaa"]("R,b,a,i,j");
    tmps_["156_aa_Lvo"].~TArrayD();

    // flops: o1v1L1  = o1v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["157_aa_Lvo"]("R,a,i")  = -1.00 * dp["aa_vv"]("a,b") * r1_1["aa"]("R,b,i");
    tmps_["157_aa_Lvo"].~TArrayD();

    // flops: o3v1L1  = o3v2L1
    //  mems: o3v1L1  = o3v1L1
    tmps_["158_aaaa_Lvooo"]("R,a,i,j,k")  = eri["aaaa_vovo"]("a,i,b,j") * r1["aa"]("R,b,k");

    // flops: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j")  = -1.00 * r1["aa"]("R,a,i") * reused_["155_aa_vo"]("b,j");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= tmps_["166_aa_Lvo"]("R,b,i") * t1_1["aa"]("a,j");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") += reused_["157_aa_vo"]("b,j") * r1["aa"]("R,a,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= reused_["58_aa_vo"]("b,j") * tmps_["122_aa_Lvo"]("R,a,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= t1_1["aa"]("b,i") * tmps_["120_aa_Lvo"]("R,a,j");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= dp["aa_vo"]("b,j") * tmps_["122_aa_Lvo"]("R,a,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= t1_1["aa"]("b,k") * tmps_["161_aaaa_Loovo"]("R,k,j,a,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= t1_1["aa"]("b,i") * tmps_["165_aa_Lov"]("R,j,a");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= tmps_["160_aaaa_Lvooo"]("R,b,j,l,i") * t1_1["aa"]("a,l");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= r1_1["aa"]("R,a,l") * reused_["150_aaaa_oovo"]("l,j,b,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") += 2.00 * r1_1["aa"]("R,a,i") * reused_["89_aa_ov"]("j,b");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") += reused_["151_aaaa_vooo"]("b,k,j,i") * r1["aa"]("R,a,k");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") += reused_["118_abba_vovo"]("b,m,f,j") * r2["abab"]("R,a,f,i,m");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= reused_["117_aaaa_vovo"]("b,k,d,j") * r2["aaaa"]("R,d,a,i,k");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= reused_["114_aabb_voov"]("b,j,n,f") * r2["abab"]("R,a,f,i,n");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") += reused_["113_aaaa_voov"]("b,j,l,d") * r2["aaaa"]("R,d,a,i,l");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= reused_["144_aaaa_vvvo"]("b,d,a,j") * r1_1["aa"]("R,d,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= reused_["146_aaaa_vvvo"]("b,d,a,j") * r1["aa"]("R,d,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= dp["aa_vo"]("b,j") * r1["aa"]("R,a,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") += reused_["55_aaaa_voov"]("b,j,l,d") * r2_1["aaaa"]("R,d,a,i,l");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= reused_["141_aaaa_vvvo"]("b,d,a,j") * r1["aa"]("R,d,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= r1_1["aa"]("R,a,l") * reused_["149_aaaa_oovo"]("l,j,b,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") += 0.50 * reused_["87_aa_vo"]("b,j") * r1_1["aa"]("R,a,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= r1["aa"]("R,a,l") * reused_["142_aaaa_oovo"]("l,j,b,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") += reused_["56_aabb_voov"]("b,j,n,f") * r2_1["abab"]("R,a,f,i,n");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") += reused_["54_aa_vo"]("b,j") * r1["aa"]("R,a,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= reused_["145_aaaa_vvvo"]("b,d,a,j") * r1_1["aa"]("R,d,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= reused_["90_aa_vo"]("b,j") * r1_1["aa"]("R,a,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= reused_["6_aabb_voov"]("b,j,n,f") * r2_1["abab"]("R,a,f,i,n");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= reused_["153_aa_vo"]("b,j") * r1["aa"]("R,a,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= eri["abba_vovo"]("b,m,e,j") * r2_1["abab"]("R,a,e,i,m");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= reused_["58_aa_vo"]("b,j") * r1["aa"]("R,a,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") += f["aa_vo"]("b,j") * r1_1["aa"]("R,a,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= reused_["154_aaaa_vooo"]("b,j,l,i") * r1["aa"]("R,a,l");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") += eri["aaaa_vovo"]("b,k,c,j") * r2_1["aaaa"]("R,c,a,i,k");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") += reused_["115_aabb_voov"]("b,j,n,f") * r2["abab"]("R,a,f,i,n");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= reused_["152_aa_vo"]("b,j") * r1_1["aa"]("R,a,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= r0_1("R") * reused_["85_aaaa_vovo"]("b,j,a,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") += reused_["54_aa_vo"]("b,j") * tmps_["122_aa_Lvo"]("R,a,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= t2_1["abab"]("b,e,j,m") * tmps_["164_aabb_Lvoov"]("R,a,i,m,e");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") -= t2["abab"]("b,e,j,m") * tmps_["159_aabb_Lvoov"]("R,a,i,m,e");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") += tmps_["163_aaaa_Lvovo"]("R,b,j,a,i");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") += r1_1["aa"]("R,a,i") * reused_["156_aa_vo"]("b,j");
    tmps_["167_aaaa_Lvovo"]("R,a,i,b,j") += t1_1["aa"]("a,k") * tmps_["158_aaaa_Lvooo"]("R,b,k,j,i");
    tmps_["163_aaaa_Lvovo"].~TArrayD();
    tmps_["161_aaaa_Loovo"].~TArrayD();
    tmps_["160_aaaa_Lvooo"].~TArrayD();
    tmps_["159_aabb_Lvoov"].~TArrayD();
    tmps_["158_aaaa_Lvooo"].~TArrayD();

    // sigmar2_1_aaaa += -1.00 P(i,j) P(a,b) d-_aa(k,c) r1_1_aa(b,i) t2_aaaa(c,a,j,k) t0_1
    //                += +1.00 P(i,j) P(a,b) d-_aa(k,c) r1_aa(c,i) t1_1_aa(a,k) t1_1_aa(b,j)
    //                += -1.00 P(i,j) P(a,b) d-_aa(a,c) r1_1_aa(c,i) t1_1_aa(b,j)
    //                += -1.00 P(i,j) P(a,b) <l,k||d,c>_abab r1_aa(b,i) t2_abab(a,c,j,k) t1_1_aa(d,l)
    //                += +1.00 P(i,j) P(a,b) <l,k||c,d>_bbbb r1_aa(b,i) t2_abab(a,c,j,k) t1_1_bb(d,l)
    //                += +1.00 P(i,j) P(a,b) d-_bb(k,c) r0_1 t2_abab(a,c,j,k) t1_1_aa(b,i)
    //                += +1.00 P(i,j) P(a,b) d-_aa(k,c) r1_aa(b,k) t1_1_aa(a,i) t1_1_aa(c,j)
    //                += +1.00 P(i,j) P(a,b) d-_aa(a,j) r0_1 t1_1_aa(b,i)
    //                += -1.00 P(i,j) P(a,b) <l,k||c,j>_aaaa r2_aaaa(c,b,i,l) t1_1_aa(a,k)
    //                += +1.00 P(i,j) P(a,b) <k,l||j,c>_abab r2_abab(b,c,i,l) t1_1_aa(a,k)
    //                += +1.00 P(i,j) P(a,b) d-_aa(k,j) r1_1_aa(b,k) t1_1_aa(a,i)
    //                += +1.00 P(i,j) P(a,b) <l,k||d,c>_abab r1_aa(d,i) t2_abab(a,c,j,k) t1_1_aa(b,l)
    //                += +1.00 P(i,j) P(a,b) <l,k||c,d>_aaaa r1_aa(d,i) t2_aaaa(c,a,j,k) t1_1_aa(b,l)
    //                += -1.00 P(i,j) P(a,b) <l,k||j,c>_abab r1_1_aa(b,l) t2_abab(a,c,i,k)
    //                += -2.00 P(i,j) P(a,b) d-_aa(k,j) r1_1_aa(b,i) t1_1_aa(a,k)
    //                += -2.00 P(i,j) P(a,b) d-_aa(k,c) r1_1_aa(b,i) t2_1_aaaa(c,a,j,k)
    //                += +2.00 P(i,j) P(a,b) d-_aa(a,c) r1_1_aa(b,i) t1_1_aa(c,j)
    //                += +2.00 P(i,j) P(a,b) d-_bb(k,c) r1_1_aa(b,i) t2_1_abab(a,c,j,k)
    //                += +1.00 P(i,j) P(a,b) <k,a||c,j>_aaaa r1_aa(b,k) t1_1_aa(c,i)
    //                += -1.00 P(i,j) P(a,b) <a,k||c,d>_abab r2_abab(b,d,i,k) t1_1_aa(c,j)
    //                += -1.00 P(i,j) P(a,b) <k,a||c,d>_aaaa r2_aaaa(d,b,i,k) t1_1_aa(c,j)
    //                += +1.00 P(i,j) P(a,b) <k,l||c,d>_abab r2_abab(b,d,i,l) t2_1_aaaa(c,a,j,k)
    //                += +1.00 P(i,j) P(a,b) <l,k||c,d>_aaaa r2_aaaa(d,b,i,l) t2_1_aaaa(c,a,j,k)
    //                += -1.00 P(i,j) P(a,b) <k,a||c,d>_aaaa r1_1_aa(d,i) t2_aaaa(c,b,j,k)
    //                += +1.00 P(i,j) P(a,b) <a,k||d,c>_abab r1_aa(d,i) t2_1_abab(b,c,j,k)
    //                += +1.00 P(i,j) P(a,b) d+_aa(a,j) r1_aa(b,i)
    //                += +1.00 P(i,j) P(a,b) <l,k||c,d>_aaaa r2_1_aaaa(d,b,i,l) t2_aaaa(c,a,j,k)
    //                += -1.00 P(i,j) P(a,b) <k,a||c,d>_aaaa r1_aa(d,i) t2_1_aaaa(c,b,j,k)
    //                += -1.00 P(i,j) P(a,b) <l,k||c,j>_aaaa r1_1_aa(b,l) t2_aaaa(c,a,i,k)
    //                += +0.50 P(i,j) P(a,b) <k,a||c,d>_aaaa r1_1_aa(b,i) t2_aaaa(c,d,j,k)
    //                += +0.50 P(i,j) P(a,b) <l,k||j,c>_abab r1_1_aa(b,i) t2_abab(a,c,l,k)
    //                += +0.50 P(i,j) P(a,b) <k,l||j,c>_abab r1_1_aa(b,i) t2_abab(a,c,k,l)
    //                += +0.50 P(i,j) P(a,b) <l,k||c,j>_aaaa r1_1_aa(b,i) t2_aaaa(c,a,l,k)
    //                += +1.00 P(i,j) P(a,b) f_aa(k,c) r1_1_aa(b,i) t2_aaaa(c,a,j,k)
    //                += -0.50 P(i,j) P(a,b) <a,k||c,d>_abab r1_1_aa(b,i) t2_abab(c,d,j,k)
    //                += -0.50 P(i,j) P(a,b) <a,k||d,c>_abab r1_1_aa(b,i) t2_abab(d,c,j,k)
    //                += +1.00 P(i,j) P(a,b) d-_aa(a,j) r1_1_aa(b,i) t0_1
    //                += -1.00 P(i,j) P(a,b) f_bb(k,c) r1_1_aa(b,i) t2_abab(a,c,j,k)
    //                += -1.00 P(i,j) P(a,b) <l,k||c,j>_aaaa r1_aa(b,l) t2_1_aaaa(c,a,i,k)
    //                += +1.00 P(i,j) P(a,b) <l,k||c,d>_aaaa r1_aa(b,l) t2_aaaa(c,a,j,k) t1_1_aa(d,i)
    //                += +1.00 P(i,j) P(a,b) <l,k||c,d>_bbbb r2_1_abab(b,d,i,l) t2_abab(a,c,j,k)
    //                += -1.00 P(i,j) P(a,b) d+_aa(k,c) r1_aa(b,i) t2_aaaa(c,a,j,k)
    //                += +1.00 P(i,j) P(a,b) <a,k||d,c>_abab r1_1_aa(d,i) t2_abab(b,c,j,k)
    //                += +1.00 P(i,j) P(a,b) d-_bb(k,k) r1_1_aa(b,i) t1_1_aa(a,j)
    //                += +1.00 P(i,j) P(a,b) d-_aa(k,k) r1_1_aa(b,i) t1_1_aa(a,j)
    //                += +1.00 P(i,j) P(a,b) <k,l||c,d>_abab r2_1_abab(b,d,i,l) t2_aaaa(c,a,j,k)
    //                += +0.50 P(i,j) P(a,b) <l,k||d,c>_abab r1_aa(b,i) t2_abab(a,c,l,k) t1_1_aa(d,j)
    //                += +0.50 P(i,j) P(a,b) <k,l||d,c>_abab r1_aa(b,i) t2_abab(a,c,k,l) t1_1_aa(d,j)
    //                += -1.00 P(i,j) P(a,b) d-_aa(k,j) r1_aa(b,i) t0_1 t1_1_aa(a,k)
    //                += -1.00 P(i,j) P(a,b) d-_aa(k,c) r1_aa(b,i) t0_1 t2_1_aaaa(c,a,j,k)
    //                += +1.00 P(i,j) P(a,b) d-_aa(a,c) r1_aa(b,i) t0_1 t1_1_aa(c,j)
    //                += +1.00 P(i,j) P(a,b) d-_bb(k,c) r1_aa(b,i) t0_1 t2_1_abab(a,c,j,k)
    //                += +1.00 P(i,j) P(a,b) <k,l||c,d>_abab r1_aa(b,i) t2_aaaa(c,a,j,k) t1_1_bb(d,l)
    //                += -2.00 P(i,j) P(a,b) d-_aa(k,c) r1_aa(b,i) t1_1_aa(a,k) t1_1_aa(c,j)
    //                += +0.50 P(i,j) P(a,b) <l,k||c,d>_abab r1_aa(b,i) t2_abab(c,d,j,k) t1_1_aa(a,l)
    //                += +0.50 P(i,j) P(a,b) <l,k||d,c>_abab r1_aa(b,i) t2_abab(d,c,j,k) t1_1_aa(a,l)
    //                += +0.50 P(i,j) P(a,b) <l,k||c,d>_aaaa r1_aa(b,i) t2_aaaa(c,a,l,k) t1_1_aa(d,j)
    //                += -1.00 P(i,j) P(a,b) r1_aa(b,i) t1_1_aa(a,j) w0
    //                += +1.00 P(i,j) P(a,b) f_aa(k,c) r1_aa(b,i) t2_1_aaaa(c,a,j,k)
    //                += +1.00 P(i,j) P(a,b) d-_bb(k,c) r1_aa(b,i) t1_1_aa(a,j) t1_1_bb(c,k)
    //                += +0.50 P(i,j) P(a,b) <l,k||j,c>_abab r1_aa(b,i) t2_1_abab(a,c,l,k)
    //                += +0.50 P(i,j) P(a,b) <k,l||j,c>_abab r1_aa(b,i) t2_1_abab(a,c,k,l)
    //                += +1.00 P(i,j) P(a,b) d-_aa(k,c) r1_aa(b,i) t1_1_aa(a,j) t1_1_aa(c,k)
    //                += -1.00 P(i,j) P(a,b) <k,a||c,j>_aaaa r1_aa(b,i) t1_1_aa(c,k)
    //                += -1.00 P(i,j) P(a,b) f_bb(k,c) r1_aa(b,i) t2_1_abab(a,c,j,k)
    //                += +1.00 P(i,j) P(a,b) f_aa(k,j) r1_aa(b,i) t1_1_aa(a,k)
    //                += +0.50 P(i,j) P(a,b) <k,a||c,d>_aaaa r1_aa(b,i) t2_1_aaaa(c,d,j,k)
    //                += -1.00 P(i,j) P(a,b) <a,k||j,c>_abab r1_aa(b,i) t1_1_bb(c,k)
    //                += +0.50 P(i,j) P(a,b) <l,k||c,j>_aaaa r1_aa(b,i) t2_1_aaaa(c,a,l,k)
    //                += -0.50 P(i,j) P(a,b) <a,k||c,d>_abab r1_aa(b,i) t2_1_abab(c,d,j,k)
    //                += -0.50 P(i,j) P(a,b) <a,k||d,c>_abab r1_aa(b,i) t2_1_abab(d,c,j,k)
    //                += -1.00 P(i,j) P(a,b) f_aa(a,c) r1_aa(b,i) t1_1_aa(c,j)
    //                += +0.50 P(i,j) P(a,b) <l,k||c,d>_aaaa r1_aa(b,i) t2_aaaa(c,d,j,k) t1_1_aa(a,l)
    //                += -1.00 P(i,j) P(a,b) <a,k||j,c>_abab r2_1_abab(b,c,i,k)
    //                += +1.00 P(i,j) P(a,b) d+_bb(k,c) r1_aa(b,i) t2_abab(a,c,j,k)
    //                += -1.00 P(i,j) P(a,b) f_aa(a,j) r1_1_aa(b,i)
    //                += +1.00 P(i,j) P(a,b) <l,k||d,c>_abab r1_aa(b,l) t2_abab(a,c,j,k) t1_1_aa(d,i)
    //                += -1.00 P(i,j) P(a,b) <l,k||j,c>_abab r1_aa(b,l) t2_1_abab(a,c,i,k)
    //                += +1.00 P(i,j) P(a,b) <k,a||c,j>_aaaa r2_1_aaaa(c,b,i,k)
    //                += +1.00 P(i,j) P(a,b) <l,k||c,d>_bbbb r2_abab(b,d,i,l) t2_1_abab(a,c,j,k)
    //                += +1.00 P(i,j) P(a,b) d-_bb(k,c) r1_1_aa(b,i) t2_abab(a,c,j,k) t0_1
    //                += -1.00 P(i,j) P(a,b) <a,k||j,c>_abab r0_1 t2_abab(b,c,i,k)
    //                += +1.00 P(i,j) P(a,b) <k,a||c,j>_aaaa r0_1 t2_aaaa(c,b,i,k)
    //                += -1.00 P(i,j) P(a,b) d-_aa(k,c) r0_1 t2_aaaa(c,a,j,k) t1_1_aa(b,i)
    //                += +1.00 P(i,j) P(a,b) <l,k||d,c>_abab r2_aaaa(d,b,i,l) t2_1_abab(a,c,j,k)
    //                += +1.00 P(i,j) P(a,b) <l,k||d,c>_abab r2_1_aaaa(d,b,i,l) t2_abab(a,c,j,k)
    //                += +1.00 P(i,j) P(a,b) d-_bb(k,c) r2_1_abab(b,c,i,k) t1_1_aa(a,j)
    //                += -1.00 P(i,j) P(a,b) d-_aa(k,c) r2_1_aaaa(c,b,i,k) t1_1_aa(a,j)
    //                += -1.00 P(i,j) P(a,b) <l,k||c,d>_aaaa r1_aa(b,i) t2_aaaa(c,a,j,k) t1_1_aa(d,l)
    //                += +1.00 P(i,j) P(a,b) <k,a||c,j>_aaaa r1_aa(c,i) t1_1_aa(b,k)
    sigmar2_1_aaaa("R,a,b,i,j") -= tmps_["167_aaaa_Lvovo"]("R,b,i,a,j");
    sigmar2_1_aaaa("R,a,b,i,j") += tmps_["167_aaaa_Lvovo"]("R,b,j,a,i");
    sigmar2_1_aaaa("R,a,b,i,j") += tmps_["167_aaaa_Lvovo"]("R,a,i,b,j");
    sigmar2_1_aaaa("R,a,b,i,j") -= tmps_["167_aaaa_Lvovo"]("R,a,j,b,i");
    tmps_["167_aaaa_Lvovo"].~TArrayD();

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["170_aa_Lov"]("L,i,a")  = l2["aaaa"]("L,i,j,b,a") * t1_1["aa"]("b,j");

    // csigmal1_aa += -1.00 t1_1_aa(b,j) l2_aaaa(i,j,b,a)
    csigmal1_aa("L,a,i") -= tmps_["170_aa_Lov"]("L,i,a");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["171_aa_Lov"]("L,i,a")  = l2["abab"]("L,i,j,a,b") * t1_1["bb"]("b,j");

    // csigmal1_aa += +1.00 t1_1_bb(b,j) l2_abab(i,j,a,b)
    csigmal1_aa("L,a,i") += tmps_["171_aa_Lov"]("L,i,a");

    // flops: o3v1L1  = o3v2L1
    //  mems: o3v1L1  = o3v1L1
    tmps_["169_abab_Loovo"]("L,i,j,a,k")  = l2_1["abab"]("L,i,j,a,b") * t1_1["bb"]("b,k");

    // flops: o2v2L1  = o3v3L1 o3v3L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1
    tmps_["168_bbaa_Lvoov"]("L,a,i,j,b")  = t2["abab"]("c,a,k,i") * l2["aaaa"]("L,j,k,c,b");
    tmps_["168_bbaa_Lvoov"]("L,a,i,j,b") += t2_1["abab"]("c,a,k,i") * l2_1["aaaa"]("L,j,k,c,b");

    // flops: o2v2L1  = o2v2L1 o3v3L1 o2v2L1 o2v2L1 o3v3L1 o3v3L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b")  = -1.00 * dp["aa_ov"]("i,a") * tmps_["170_aa_Lov"]("L,j,b");
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b") -= eri["aaaa_vovo"]("d,i,a,l") * l2["aaaa"]("L,j,l,d,b");
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b") -= f["aa_ov"]("i,a") * l1["aa"]("L,j,b");
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b") -= l1_1["aa"]("L,j,b") * reused_["9_aa_ov"]("i,a");
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b") -= l2["abab"]("L,j,m,b,c") * reused_["158_bbaa_voov"]("c,m,i,a");
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b") -= reused_["69_aaaa_voov"]("d,l,i,a") * l2["aaaa"]("L,j,l,d,b");
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b") += reused_["105_baab_vovo"]("c,i,a,m") * l2_1["abab"]("L,j,m,b,c");
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b") += eri["baab_vovo"]("c,i,a,m") * l2["abab"]("L,j,m,b,c");
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b") += dp["aa_ov"]("i,a") * l1_1["aa"]("L,j,b");
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b") += l1["aa"]("L,j,b") * reused_["13_aa_ov"]("i,a");
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b") += l2["abab"]("L,j,m,b,c") * reused_["1_bbaa_voov"]("c,m,i,a");
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b") -= l2_1["aaaa"]("L,j,l,d,b") * reused_["160_aaaa_vovo"]("d,i,a,l");
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b") -= l2_1["abab"]("L,j,m,b,c") * reused_["159_bbaa_voov"]("c,m,i,a");
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,j,m,b,c") * reused_["103_bbaa_voov"]("c,m,i,a");
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b") -= l1_1["aa"]("L,j,b") * reused_["147_aa_ov"]("i,a");
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b") += eri["abab_oovv"]("i,n,a,e") * tmps_["168_bbaa_Lvoov"]("L,e,n,j,b");
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b") += eri["abab_oovo"]("i,n,a,m") * tmps_["169_abab_Loovo"]("L,j,m,b,n");
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b") += dp["aa_ov"]("i,a") * tmps_["171_aa_Lov"]("L,j,b");
    tmps_["172_aaaa_Lovov"]("L,i,a,j,b") -= eri["aaaa_oovo"]("i,k,a,l") * tmps_["138_aaaa_Loovo"]("L,j,l,b,k");
    tmps_["168_bbaa_Lvoov"].~TArrayD();

    // sigmal2_aaaa += +1.00 P(i,j) P(a,b) <j,l||a,k>_abab t1_1_bb(c,l) l2_1_abab(i,k,b,c)
    //              += -1.00 P(i,j) P(a,b) <j,c||a,d>_abab t1_1_bb(d,k) l2_1_abab(i,k,b,c)
    //              += +1.00 P(i,j) P(a,b) <j,l||d,a>_aaaa t2_aaaa(d,c,k,l) l2_aaaa(i,k,c,b)
    //              += -1.00 P(i,j) P(a,b) <j,c||a,k>_abab l2_abab(i,k,b,c)
    //              += +1.00 P(i,j) P(a,b) d+_aa(j,a) l1_1_aa(i,b)
    //              += +1.00 P(i,j) P(a,b) <j,l||d,a>_aaaa t2_abab(d,c,l,k) l2_abab(i,k,b,c)
    //              += +1.00 P(i,j) P(a,b) d-_aa(j,a) t0_1 l1_aa(i,b)
    //              += -1.00 P(i,j) P(a,b) <j,k||a,c>_abab t1_1_bb(c,k) l1_1_aa(i,b)
    //              += +1.00 P(i,j) P(a,b) <j,l||a,d>_abab t2_bbbb(d,c,k,l) l2_abab(i,k,b,c)
    //              += -1.00 P(i,j) P(a,b) <j,c||d,a>_aaaa t1_1_aa(d,k) l2_1_aaaa(i,k,c,b)
    //              += +1.00 P(i,j) P(a,b) <j,l||d,a>_aaaa t2_1_aaaa(d,c,k,l) l2_1_aaaa(i,k,c,b)
    //              += +1.00 P(i,j) P(a,b) <j,l||d,a>_aaaa t2_1_abab(d,c,l,k) l2_1_abab(i,k,b,c)
    //              += +1.00 P(i,j) P(a,b) <j,l||a,d>_abab t2_1_bbbb(d,c,k,l) l2_1_abab(i,k,b,c)
    //              += -1.00 P(i,j) P(a,b) f_aa(j,a) l1_aa(i,b)
    //              += +1.00 P(i,j) P(a,b) <j,k||c,a>_aaaa t1_1_aa(c,k) l1_1_aa(i,b)
    //              += +1.00 P(i,j) P(a,b) <j,c||a,k>_aaaa l2_aaaa(i,k,c,b)
    //              += +1.00 P(i,j) P(a,b) <j,l||a,d>_abab t2_1_abab(c,d,k,l) l2_1_aaaa(i,k,c,b)
    //              += +1.00 P(i,j) P(a,b) <j,l||a,d>_abab t2_abab(c,d,k,l) l2_aaaa(i,k,c,b)
    //              += -1.00 P(i,j) P(a,b) d-_aa(j,a) t1_1_aa(c,k) l2_aaaa(i,k,c,b)
    //              += +1.00 P(i,j) P(a,b) d-_aa(j,a) t1_1_bb(c,k) l2_abab(i,k,b,c)
    //              += -1.00 P(i,j) P(a,b) <j,l||a,k>_aaaa t1_1_aa(c,l) l2_1_aaaa(i,k,c,b)
    sigmal2_aaaa("L,a,b,i,j") += tmps_["172_aaaa_Lovov"]("L,j,a,i,b");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["172_aaaa_Lovov"]("L,i,a,j,b");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["172_aaaa_Lovov"]("L,j,b,i,a");
    sigmal2_aaaa("L,a,b,i,j") += tmps_["172_aaaa_Lovov"]("L,i,b,j,a");
    tmps_["172_aaaa_Lovov"].~TArrayD();

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["173_aa_Lov"]("L,i,a")  = l2_1["aaaa"]("L,i,j,b,a") * t1_1["aa"]("b,j");

    // csigmal1_1_aa += -1.00 t1_1_aa(b,j) l2_1_aaaa(i,j,b,a)
    csigmal1_1_aa("L,a,i") -= tmps_["173_aa_Lov"]("L,i,a");

    // flops: o2v2L1  = o3v3L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["174_bbaa_Lvoov"]("L,a,i,j,b")  = t2["abab"]("c,a,k,i") * l2_1["aaaa"]("L,j,k,c,b");

    // flops: o2v2L1  = o1v1L1 o2v2L1 o2v2L1 o2v2 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o1v1L1 o2v2L1 o2v2L1 o2v2 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["175_aaaa_Lovov"]("L,i,a,j,b")  = (tmps_["102_aa_Lov"]("L,i,a") + -1.00 * tmps_["173_aa_Lov"]("L,i,a")) * dp["aa_ov"]("j,b");
    tmps_["175_aaaa_Lovov"]("L,i,a,j,b") -= l1_1["aa"]("L,i,a") * f["aa_ov"]("j,b");
    tmps_["175_aaaa_Lovov"]("L,i,a,j,b") += (eri["baab_vovo"]("d,j,b,l") + -1.00 * reused_["158_bbaa_voov"]("d,l,j,b")) * l2_1["abab"]("L,i,l,a,d");
    tmps_["175_aaaa_Lovov"]("L,i,a,j,b") += l1["aa"]("L,i,a") * dp["aa_ov"]("j,b");
    tmps_["175_aaaa_Lovov"]("L,i,a,j,b") += l1_1["aa"]("L,i,a") * reused_["13_aa_ov"]("j,b");
    tmps_["175_aaaa_Lovov"]("L,i,a,j,b") -= l2_1["aaaa"]("L,i,m,e,a") * reused_["69_aaaa_voov"]("e,m,j,b");
    tmps_["175_aaaa_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,i,l,a,d") * reused_["1_bbaa_voov"]("d,l,j,b");
    tmps_["175_aaaa_Lovov"]("L,i,a,j,b") -= l2_1["aaaa"]("L,i,m,e,a") * eri["aaaa_vovo"]("e,j,b,m");
    tmps_["175_aaaa_Lovov"]("L,i,a,j,b") += eri["abab_oovv"]("j,k,b,c") * tmps_["174_bbaa_Lvoov"]("L,c,k,i,a");

    // sigmal2_1_aaaa += +1.00 P(i,j) P(a,b) d-_aa(j,a) t1_1_bb(c,k) l2_1_abab(i,k,b,c)
    //                += -1.00 P(i,j) P(a,b) d-_aa(j,a) t1_1_aa(c,k) l2_1_aaaa(i,k,c,b)
    //                += +1.00 P(i,j) P(a,b) <j,l||a,d>_abab t2_abab(c,d,k,l) l2_1_aaaa(i,k,c,b)
    //                += -1.00 P(i,j) P(a,b) <j,c||a,k>_abab l2_1_abab(i,k,b,c)
    //                += +1.00 P(i,j) P(a,b) <j,l||d,a>_aaaa t2_abab(d,c,l,k) l2_1_abab(i,k,b,c)
    //                += -1.00 P(i,j) P(a,b) f_aa(j,a) l1_1_aa(i,b)
    //                += +1.00 P(i,j) P(a,b) d-_aa(j,a) l1_aa(i,b)
    //                += +1.00 P(i,j) P(a,b) d-_aa(j,a) t0_1 l1_1_aa(i,b)
    //                += +1.00 P(i,j) P(a,b) <j,l||d,a>_aaaa t2_aaaa(d,c,k,l) l2_1_aaaa(i,k,c,b)
    //                += +1.00 P(i,j) P(a,b) <j,l||a,d>_abab t2_bbbb(d,c,k,l) l2_1_abab(i,k,b,c)
    //                += +1.00 P(i,j) P(a,b) <j,c||a,k>_aaaa l2_1_aaaa(i,k,c,b)
    sigmal2_1_aaaa("L,a,b,i,j") += tmps_["175_aaaa_Lovov"]("L,i,b,j,a");
    sigmal2_1_aaaa("L,a,b,i,j") -= tmps_["175_aaaa_Lovov"]("L,j,b,i,a");
    sigmal2_1_aaaa("L,a,b,i,j") -= tmps_["175_aaaa_Lovov"]("L,i,a,j,b");
    sigmal2_1_aaaa("L,a,b,i,j") += tmps_["175_aaaa_Lovov"]("L,j,a,i,b");
    tmps_["175_aaaa_Lovov"].~TArrayD();

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["179_bb_Loo"]("L,i,j")  = l2_1["bbbb"]("L,i,k,a,b") * t2["bbbb"]("a,b,k,j");

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["178_bb_Loo"]("L,i,j")  = l2_1["abab"]("L,k,i,a,b") * t2["abab"]("a,b,k,j");

    // flops: o2v2L1  = o1v1L1 o2v2L1
    //  mems: o2v2L1  = o1v1L1 o2v2L1
    tmps_["176_aaaa_Lovov"]("L,i,a,j,b")  = (tmps_["102_aa_Lov"]("L,i,a") + -1.00 * tmps_["173_aa_Lov"]("L,i,a")) * dp["aa_ov"]("j,b");
    tmps_["176_aaaa_Lovov"].~TArrayD();

    // flops: o2v2L1  = o3v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["177_bbbb_Lovvo"]("L,i,a,b,j")  = -1.00 * l2_1["bbbb"]("L,i,k,a,b") * reused_["35_bb_oo"]("j,k");

    // flops: o2v2L1  = o3v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["180_bbbb_Lovvo"]("L,i,a,b,j")  = -1.00 * l2["bbbb"]("L,i,m,a,b") * dp["bb_oo"]("j,m");
    tmps_["180_bbbb_Lovvo"]("L,i,a,b,j") += l2_1["bbbb"]("L,i,m,a,b") * reused_["30_bb_oo"]("j,m");
    tmps_["180_bbbb_Lovvo"]("L,i,a,b,j") += 0.50 * l2_1["bbbb"]("L,i,m,a,b") * reused_["161_bb_oo"]("j,m");
    tmps_["180_bbbb_Lovvo"]("L,i,a,b,j") -= eri["bbbb_vovv"]("d,j,a,b") * l1_1["bb"]("L,i,d");
    tmps_["180_bbbb_Lovvo"]("L,i,a,b,j") += f["bb_oo"]("j,m") * l2_1["bbbb"]("L,i,m,a,b");
    tmps_["180_bbbb_Lovvo"]("L,i,a,b,j") += tmps_["177_bbbb_Lovvo"]("L,i,a,b,j");
    tmps_["180_bbbb_Lovvo"]("L,i,a,b,j") += 0.50 * eri["bbbb_oovv"]("j,k,a,b") * tmps_["179_bb_Loo"]("L,i,k");
    tmps_["180_bbbb_Lovvo"]("L,i,a,b,j") -= eri["bbbb_oovv"]("j,k,a,b") * tmps_["178_bb_Loo"]("L,i,k");
    tmps_["177_bbbb_Lovvo"].~TArrayD();

    // sigmal2_1_bbbb += -0.50 P(i,j) <l,j||c,d>_abab t2_abab(c,d,l,k) l2_1_bbbb(i,k,a,b)
    //                += -0.50 P(i,j) <l,j||d,c>_abab t2_abab(d,c,l,k) l2_1_bbbb(i,k,a,b)
    //                += +1.00 P(i,j) d-_bb(j,k) l2_bbbb(i,k,a,b)
    //                += -0.50 P(i,j) <j,l||c,d>_bbbb t2_bbbb(c,d,k,l) l2_1_bbbb(i,k,a,b)
    //                += -1.00 P(i,j) <j,c||a,b>_bbbb l1_1_bb(i,c)
    //                += -1.00 P(i,j) f_bb(j,k) l2_1_bbbb(i,k,a,b)
    //                += +1.00 P(i,j) d-_bb(j,k) t0_1 l2_1_bbbb(i,k,a,b)
    //                += -0.50 P(i,j) <j,l||a,b>_bbbb t2_bbbb(d,c,k,l) l2_1_bbbb(i,k,d,c)
    //                += +0.50 P(i,j) <j,l||a,b>_bbbb t2_abab(d,c,k,l) l2_1_abab(k,i,d,c)
    //                += +0.50 P(i,j) <j,l||a,b>_bbbb t2_abab(c,d,k,l) l2_1_abab(k,i,c,d)
    sigmal2_1_bbbb("L,a,b,i,j") -= tmps_["180_bbbb_Lovvo"]("L,i,a,b,j");
    sigmal2_1_bbbb("L,a,b,i,j") += tmps_["180_bbbb_Lovvo"]("L,j,a,b,i");
    tmps_["180_bbbb_Lovvo"].~TArrayD();

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["182_bb_Lov"]("L,i,a")  = l2["abab"]("L,j,i,b,a") * t1_1["aa"]("b,j");

    // csigmal1_bb += +1.00 t1_1_aa(b,j) l2_abab(j,i,b,a)
    csigmal1_bb("L,a,i") += tmps_["182_bb_Lov"]("L,i,a");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["183_bb_Lov"]("L,i,a")  = l2["bbbb"]("L,i,j,b,a") * t1_1["bb"]("b,j");

    // csigmal1_bb += -1.00 t1_1_bb(b,j) l2_bbbb(i,j,b,a)
    csigmal1_bb("L,a,i") -= tmps_["183_bb_Lov"]("L,i,a");

    // flops: o3v1L1  = o3v2L1
    //  mems: o3v1L1  = o3v1L1
    tmps_["181_abba_Loovo"]("L,i,j,a,k")  = l2_1["abab"]("L,i,j,b,a") * t1_1["aa"]("b,k");

    // flops: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b")  = -1.00 * l1["bb"]("L,i,a") * f["bb_ov"]("j,b");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") -= l1_1["bb"]("L,i,a") * reused_["8_bb_ov"]("j,b");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") -= l1_1["bb"]("L,i,a") * reused_["148_bb_ov"]("j,b");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") -= l2["bbbb"]("L,i,n,d,a") * reused_["163_bbbb_voov"]("d,n,j,b");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") += l1["bb"]("L,i,a") * reused_["12_bb_ov"]("j,b");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") += l2["bbbb"]("L,i,n,d,a") * reused_["26_bbbb_voov"]("d,n,j,b");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") -= l2_1["abab"]("L,l,i,c,a") * reused_["162_aabb_voov"]("c,l,j,b");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") += l2_1["bbbb"]("L,i,n,d,a") * reused_["102_bbbb_voov"]("d,n,j,b");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") += l1_1["bb"]("L,i,a") * dp["bb_ov"]("j,b");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") += l2["abab"]("L,l,i,c,a") * reused_["6_aabb_voov"]("c,l,j,b");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") += eri["abba_vovo"]("c,j,b,l") * l2["abab"]("L,l,i,c,a");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") += l2_1["abab"]("L,l,i,c,a") * reused_["114_aabb_voov"]("c,l,j,b");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") -= l2_1["bbbb"]("L,i,n,d,a") * reused_["164_bbbb_voov"]("d,n,j,b");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") -= l2_1["abab"]("L,l,i,c,a") * reused_["118_abba_vovo"]("c,j,b,l");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") -= l2["abab"]("L,l,i,c,a") * reused_["70_aabb_voov"]("c,l,j,b");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") -= l2["bbbb"]("L,i,n,d,a") * eri["bbbb_vovo"]("d,j,b,n");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") -= tmps_["113_bbbb_Loovo"]("L,i,n,a,m") * eri["bbbb_oovo"]("j,m,b,n");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") -= tmps_["183_bb_Lov"]("L,i,a") * dp["bb_ov"]("j,b");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") += tmps_["182_bb_Lov"]("L,i,a") * dp["bb_ov"]("j,b");
    tmps_["184_bbbb_Lovov"]("L,i,a,j,b") -= tmps_["181_abba_Loovo"]("L,l,i,a,k") * eri["abba_oovo"]("k,j,b,l");

    // sigmal2_bbbb += +1.00 P(i,j) P(a,b) d-_bb(j,a) t0_1 l1_bb(i,b)
    //              += +1.00 P(i,j) P(a,b) <j,l||d,a>_bbbb t2_bbbb(d,c,k,l) l2_bbbb(i,k,c,b)
    //              += +1.00 P(i,j) P(a,b) <j,k||c,a>_bbbb t1_1_bb(c,k) l1_1_bb(i,b)
    //              += +1.00 P(i,j) P(a,b) <l,j||d,a>_abab t2_abab(d,c,l,k) l2_bbbb(i,k,c,b)
    //              += -1.00 P(i,j) P(a,b) <k,j||c,a>_abab t1_1_aa(c,k) l1_1_bb(i,b)
    //              += +1.00 P(i,j) P(a,b) <j,l||d,a>_bbbb t2_1_abab(c,d,k,l) l2_1_abab(k,i,c,b)
    //              += +1.00 P(i,j) P(a,b) <l,j||d,a>_abab t2_1_abab(d,c,l,k) l2_1_bbbb(i,k,c,b)
    //              += +1.00 P(i,j) P(a,b) d+_bb(j,a) l1_1_bb(i,b)
    //              += -1.00 P(i,j) P(a,b) f_bb(j,a) l1_bb(i,b)
    //              += +1.00 P(i,j) P(a,b) <l,j||d,a>_abab t2_aaaa(d,c,k,l) l2_abab(k,i,c,b)
    //              += -1.00 P(i,j) P(a,b) <c,j||k,a>_abab l2_abab(k,i,c,b)
    //              += +1.00 P(i,j) P(a,b) <l,j||d,a>_abab t2_1_aaaa(d,c,k,l) l2_1_abab(k,i,c,b)
    //              += +1.00 P(i,j) P(a,b) <j,l||d,a>_bbbb t2_1_bbbb(d,c,k,l) l2_1_bbbb(i,k,c,b)
    //              += -1.00 P(i,j) P(a,b) <j,c||d,a>_bbbb t1_1_bb(d,k) l2_1_bbbb(i,k,c,b)
    //              += -1.00 P(i,j) P(a,b) <c,j||d,a>_abab t1_1_aa(d,k) l2_1_abab(k,i,c,b)
    //              += +1.00 P(i,j) P(a,b) <j,l||d,a>_bbbb t2_abab(c,d,k,l) l2_abab(k,i,c,b)
    //              += +1.00 P(i,j) P(a,b) <j,c||a,k>_bbbb l2_bbbb(i,k,c,b)
    //              += -1.00 P(i,j) P(a,b) <j,l||a,k>_bbbb t1_1_bb(c,l) l2_1_bbbb(i,k,c,b)
    //              += -1.00 P(i,j) P(a,b) d-_bb(j,a) t1_1_bb(c,k) l2_bbbb(i,k,c,b)
    //              += +1.00 P(i,j) P(a,b) d-_bb(j,a) t1_1_aa(c,k) l2_abab(k,i,c,b)
    //              += +1.00 P(i,j) P(a,b) <l,j||k,a>_abab t1_1_aa(c,l) l2_1_abab(k,i,c,b)
    sigmal2_bbbb("L,a,b,i,j") += tmps_["184_bbbb_Lovov"]("L,i,b,j,a");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["184_bbbb_Lovov"]("L,j,b,i,a");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["184_bbbb_Lovov"]("L,i,a,j,b");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["184_bbbb_Lovov"]("L,j,a,i,b");
    tmps_["184_bbbb_Lovov"].~TArrayD();

    // flops: o0v2L1  = o2v3L1
    //  mems: o0v2L1  = o0v2L1
    tmps_["186_aa_Lvv"]("L,a,b")  = l2_1["aaaa"]("L,i,j,c,a") * t2["aaaa"]("b,c,i,j");

    // flops: o2v2L1  = o3v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["185_aaaa_Lvoov"]("L,a,i,j,b")  = -4.00 * dp["aa_ov"]("k,a") * tmps_["138_aaaa_Loovo"]("L,i,j,b,k");

    // flops: o2v2L1  = o0v2 o0v2 o0v2 o2v3L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1
    //  mems: o2v2L1  = o0v2 o0v2 o0v2 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["187_aaaa_Lvoov"]("L,a,i,j,b")  = (2.00 * reused_["57_aa_vv"]("a,e") + reused_["18_aa_vv"]("a,e") + -2.00 * f["aa_vv"]("e,a") + 2.00 * reused_["59_aa_vv"]("e,a")) * l2_1["aaaa"]("L,i,j,e,b");
    tmps_["187_aaaa_Lvoov"]("L,a,i,j,b") += 2.00 * eri["aaaa_oovo"]("i,j,a,k") * l1_1["aa"]("L,k,b");
    tmps_["187_aaaa_Lvoov"]("L,a,i,j,b") += 2.00 * dp["aa_vv"]("e,a") * l2["aaaa"]("L,i,j,e,b");
    tmps_["187_aaaa_Lvoov"]("L,a,i,j,b") += tmps_["185_aaaa_Lvoov"]("L,a,i,j,b");
    tmps_["187_aaaa_Lvoov"]("L,a,i,j,b") -= eri["aaaa_oovv"]("i,j,a,c") * tmps_["186_aa_Lvv"]("L,b,c");
    tmps_["187_aaaa_Lvoov"]("L,a,i,j,b") += 2.00 * eri["aaaa_oovv"]("i,j,a,c") * tmps_["99_aa_Lvv"]("L,b,c");
    tmps_["185_aaaa_Lvoov"].~TArrayD();

    // sigmal2_1_aaaa += -0.50 P(a,b) <l,k||d,a>_aaaa t2_aaaa(d,c,l,k) l2_1_aaaa(i,j,c,b)
    //                += +1.00 P(a,b) f_aa(c,a) l2_1_aaaa(i,j,c,b)
    //                += -0.50 P(a,b) <l,k||a,d>_abab t2_abab(c,d,l,k) l2_1_aaaa(i,j,c,b)
    //                += -0.50 P(a,b) <k,l||a,d>_abab t2_abab(c,d,k,l) l2_1_aaaa(i,j,c,b)
    //                += -1.00 P(a,b) d-_aa(c,a) t0_1 l2_1_aaaa(i,j,c,b)
    //                += -1.00 P(a,b) <i,j||a,k>_aaaa l1_1_aa(k,b)
    //                += -1.00 P(a,b) d-_aa(c,a) l2_aaaa(i,j,c,b)
    //                += +2.00 P(a,b) d-_aa(k,a) t1_1_aa(c,k) l2_1_aaaa(i,j,c,b)
    //                += -0.50 P(a,b) <i,j||d,a>_aaaa t2_aaaa(d,c,k,l) l2_1_aaaa(k,l,c,b)
    //                += +0.50 P(a,b) <i,j||d,a>_aaaa t2_abab(d,c,k,l) l2_1_abab(k,l,b,c)
    //                += +0.50 P(a,b) <i,j||d,a>_aaaa t2_abab(d,c,l,k) l2_1_abab(l,k,b,c)
    sigmal2_1_aaaa("L,a,b,i,j") -= 0.50 * tmps_["187_aaaa_Lvoov"]("L,a,i,j,b");
    sigmal2_1_aaaa("L,a,b,i,j") += 0.50 * tmps_["187_aaaa_Lvoov"]("L,b,i,j,a");
    tmps_["187_aaaa_Lvoov"].~TArrayD();

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["188_aa_Lvo"]("R,a,i")  = -1.00 * dp["bb_ov"]("j,b") * r2["abab"]("R,a,b,i,j");
    tmps_["188_aa_Lvo"].~TArrayD();

    // flops: o1v1L1  = o2v2L1 o2v2L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1
    tmps_["189_aa_Lvo"]("R,a,i")  = -1.00 * dp["bb_ov"]("k,c") * r2["abab"]("R,a,c,i,k");
    tmps_["189_aa_Lvo"]("R,a,i") += dp["aa_ov"]("j,b") * r2["aaaa"]("R,b,a,i,j");

    // sigmar1_1_aa += +1.00 d+_aa(j,b) r2_aaaa(b,a,i,j)
    //              += -1.00 d+_bb(j,b) r2_abab(a,b,i,j)
    sigmar1_1_aa("R,a,i") += tmps_["189_aa_Lvo"]("R,a,i");

    // flops: o1v1L1  = o1v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["191_aa_Lvo"]("R,a,i")  = dp["aa_vv"]("a,b") * r1["aa"]("R,b,i");

    // sigmar1_1_aa += -1.00 d+_aa(a,b) r1_aa(b,i)
    sigmar1_1_aa("R,a,i") -= tmps_["191_aa_Lvo"]("R,a,i");

    // flops: o1v1L1  = o2v1L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["192_aa_Lov"]("R,i,a")  = dp["aa_oo"]("j,i") * r1["aa"]("R,a,j");

    // sigmar1_1_aa += +1.00 d+_aa(j,i) r1_aa(a,j)
    sigmar1_1_aa("R,a,i") += tmps_["192_aa_Lov"]("R,i,a");

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["190_aaaa_Lvovo"]("R,a,i,b,j")  = -1.00 * t1_1["aa"]("a,i") * tmps_["189_aa_Lvo"]("R,b,j");

    // flops: o2v2L1  = o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j")  = eri["abba_vovo"]("a,l,d,i") * r2["abab"]("R,b,d,j,l");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") -= f["aa_vo"]("a,i") * r1["aa"]("R,b,j");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") -= r1["aa"]("R,b,j") * reused_["89_aa_ov"]("i,a");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") -= reused_["55_aaaa_voov"]("a,i,m,e") * r2["aaaa"]("R,e,b,j,m");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") += reused_["144_aaaa_vvvo"]("a,e,b,i") * r1["aa"]("R,e,j");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") += reused_["145_aaaa_vvvo"]("a,e,b,i") * r1["aa"]("R,e,j");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") += r1["aa"]("R,b,m") * reused_["149_aaaa_oovo"]("m,i,a,j");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") += reused_["152_aa_vo"]("a,i") * r1["aa"]("R,b,j");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") -= eri["aaaa_vovo"]("a,k,c,i") * r2["aaaa"]("R,c,b,j,k");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") += dp["aa_vo"]("a,i") * r1_1["aa"]("R,b,j");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") -= reused_["54_aa_vo"]("a,i") * r1_1["aa"]("R,b,j");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") += reused_["6_aabb_voov"]("a,i,n,f") * r2["abab"]("R,b,f,j,n");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") += r1["aa"]("R,b,m") * reused_["150_aaaa_oovo"]("m,i,a,j");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") -= 0.50 * reused_["87_aa_vo"]("a,i") * r1["aa"]("R,b,j");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") += reused_["58_aa_vo"]("a,i") * r1_1["aa"]("R,b,j");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") += reused_["90_aa_vo"]("a,i") * r1["aa"]("R,b,j");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") -= reused_["56_aabb_voov"]("a,i,n,f") * r2["abab"]("R,b,f,j,n");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") += t2["abab"]("a,d,i,l") * tmps_["164_aabb_Lvoov"]("R,b,j,l,d");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") += t1_1["aa"]("a,j") * tmps_["192_aa_Lov"]("R,i,b");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") += tmps_["190_aaaa_Lvovo"]("R,a,i,b,j");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") -= t1_1["aa"]("b,i") * tmps_["191_aa_Lvo"]("R,a,j");
    tmps_["193_aaaa_Lvovo"]("R,a,i,b,j") += r1["aa"]("R,b,j") * reused_["229_aa_vo"]("a,i");
    tmps_["190_aaaa_Lvovo"].~TArrayD();
    tmps_["164_aabb_Lvoov"].~TArrayD();

    // sigmar2_aaaa += -1.00 P(i,j) P(a,b) <a,k||j,c>_abab r2_abab(b,c,i,k)
    //              += -1.00 P(i,j) P(a,b) f_aa(a,j) r1_aa(b,i)
    //              += -1.00 P(i,j) P(a,b) d-_aa(k,j) r1_aa(b,i) t1_1_aa(a,k)
    //              += -1.00 P(i,j) P(a,b) d-_aa(k,c) r1_aa(b,i) t2_1_aaaa(c,a,j,k)
    //              += +1.00 P(i,j) P(a,b) d-_aa(a,c) r1_aa(b,i) t1_1_aa(c,j)
    //              += +1.00 P(i,j) P(a,b) d-_bb(k,c) r1_aa(b,i) t2_1_abab(a,c,j,k)
    //              += +1.00 P(i,j) P(a,b) <l,k||c,d>_aaaa r2_aaaa(d,b,i,l) t2_aaaa(c,a,j,k)
    //              += -1.00 P(i,j) P(a,b) <k,a||c,d>_aaaa r1_aa(d,i) t2_aaaa(c,b,j,k)
    //              += +1.00 P(i,j) P(a,b) <a,k||d,c>_abab r1_aa(d,i) t2_abab(b,c,j,k)
    //              += -1.00 P(i,j) P(a,b) <l,k||c,j>_aaaa r1_aa(b,l) t2_aaaa(c,a,i,k)
    //              += +1.00 P(i,j) P(a,b) d-_bb(k,c) r1_aa(b,i) t2_abab(a,c,j,k) t0_1
    //              += +1.00 P(i,j) P(a,b) <k,a||c,j>_aaaa r2_aaaa(c,b,i,k)
    //              += +1.00 P(i,j) P(a,b) d-_aa(a,j) r1_1_aa(b,i)
    //              += -1.00 P(i,j) P(a,b) d-_aa(k,c) r1_1_aa(b,i) t2_aaaa(c,a,j,k)
    //              += +1.00 P(i,j) P(a,b) <k,l||c,d>_abab r2_abab(b,d,i,l) t2_aaaa(c,a,j,k)
    //              += -1.00 P(i,j) P(a,b) <l,k||j,c>_abab r1_aa(b,l) t2_abab(a,c,i,k)
    //              += +0.50 P(i,j) P(a,b) <k,a||c,d>_aaaa r1_aa(b,i) t2_aaaa(c,d,j,k)
    //              += +0.50 P(i,j) P(a,b) <l,k||j,c>_abab r1_aa(b,i) t2_abab(a,c,l,k)
    //              += +0.50 P(i,j) P(a,b) <k,l||j,c>_abab r1_aa(b,i) t2_abab(a,c,k,l)
    //              += +0.50 P(i,j) P(a,b) <l,k||c,j>_aaaa r1_aa(b,i) t2_aaaa(c,a,l,k)
    //              += +1.00 P(i,j) P(a,b) f_aa(k,c) r1_aa(b,i) t2_aaaa(c,a,j,k)
    //              += -0.50 P(i,j) P(a,b) <a,k||c,d>_abab r1_aa(b,i) t2_abab(c,d,j,k)
    //              += -0.50 P(i,j) P(a,b) <a,k||d,c>_abab r1_aa(b,i) t2_abab(d,c,j,k)
    //              += +1.00 P(i,j) P(a,b) d-_aa(a,j) r1_aa(b,i) t0_1
    //              += -1.00 P(i,j) P(a,b) f_bb(k,c) r1_aa(b,i) t2_abab(a,c,j,k)
    //              += +1.00 P(i,j) P(a,b) d-_bb(k,c) r1_1_aa(b,i) t2_abab(a,c,j,k)
    //              += +1.00 P(i,j) P(a,b) d-_bb(k,k) r1_aa(b,i) t1_1_aa(a,j)
    //              += +1.00 P(i,j) P(a,b) d-_aa(k,k) r1_aa(b,i) t1_1_aa(a,j)
    //              += +1.00 P(i,j) P(a,b) <l,k||c,d>_bbbb r2_abab(b,d,i,l) t2_abab(a,c,j,k)
    //              += +1.00 P(i,j) P(a,b) <l,k||d,c>_abab r2_aaaa(d,b,i,l) t2_abab(a,c,j,k)
    //              += +1.00 P(i,j) P(a,b) d-_aa(k,j) r1_aa(b,k) t1_1_aa(a,i)
    //              += -1.00 P(i,j) P(a,b) d-_aa(k,c) r2_aaaa(c,b,i,k) t1_1_aa(a,j)
    //              += +1.00 P(i,j) P(a,b) d-_bb(k,c) r2_abab(b,c,i,k) t1_1_aa(a,j)
    //              += -1.00 P(i,j) P(a,b) d-_aa(a,c) r1_aa(c,i) t1_1_aa(b,j)
    //              += -1.00 P(i,j) P(a,b) d-_aa(k,c) r1_aa(b,i) t2_aaaa(c,a,j,k) t0_1
    sigmar2_aaaa("R,a,b,i,j") += tmps_["193_aaaa_Lvovo"]("R,a,j,b,i");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["193_aaaa_Lvovo"]("R,a,i,b,j");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["193_aaaa_Lvovo"]("R,b,j,a,i");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["193_aaaa_Lvovo"]("R,b,i,a,j");
    tmps_["193_aaaa_Lvovo"].~TArrayD();

    // flops: o1v1L1  = o3v2L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    tmps_["197_aa_Lov"]("L,i,a")  = -0.50 * l2_1["aaaa"]("L,j,m,e,a") * reused_["82_aaaa_vooo"]("e,j,m,i");
    tmps_["197_aa_Lov"]("L,i,a") += l2_1["abab"]("L,i,l,a,b") * reused_["33_bb_vo"]("b,l");
    tmps_["197_aa_Lov"]("L,i,a") -= l2_1["abab"]("L,i,l,a,b") * reused_["32_bb_vo"]("b,l");
    tmps_["197_aa_Lov"]("L,i,a") += l2_1["abab"]("L,j,k,a,b") * reused_["83_baba_vooo"]("b,j,k,i");

    // sigmal1_aa += +1.00 d+_bb(k,c) t2_bbbb(c,b,j,k) l2_1_abab(i,j,a,b)
    //            += -0.50 d+_aa(i,c) t2_aaaa(c,b,j,k) l2_1_aaaa(j,k,b,a)
    //            += -1.00 d+_aa(k,c) t2_abab(c,b,k,j) l2_1_abab(i,j,a,b)
    //            += +0.50 d+_aa(i,c) t2_abab(c,b,j,k) l2_1_abab(j,k,a,b)
    //            += +0.50 d+_aa(i,c) t2_abab(c,b,k,j) l2_1_abab(k,j,a,b)
    sigmal1_aa("L,a,i") += tmps_["197_aa_Lov"]("L,i,a");

    // flops: o0v0L1  = o1v1L1 o1v1L1 o0v0L1
    //  mems: o0v0L1  = o0v0L1 o0v0L1 o0v0L1
    tmps_["210_L"]("L")  = l1_1["aa"]("L,j,b") * t1_1["aa"]("b,j");
    tmps_["210_L"]("L") += l1_1["bb"]("L,i,a") * t1_1["bb"]("a,i");

    // csigmal0_1 += +1.00 t1_1_bb(a,i) l1_1_bb(i,a)
    //            += +1.00 t1_1_aa(a,i) l1_1_aa(i,a)
    csigmal0_1("L") += tmps_["210_L"]("L");

    // flops: o2v0L1  = o2v1L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["209_aa_Loo"]("L,i,j")  = l1_1["aa"]("L,i,a") * t1_1["aa"]("a,j");

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["208_bb_Loo"]("L,i,j")  = l2_1["bbbb"]("L,k,i,a,b") * t2["bbbb"]("a,b,k,j");

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["207_aa_Loo"]("L,i,j")  = l2_1["abab"]("L,i,k,a,b") * t2_1["abab"]("a,b,j,k");

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["206_aa_Loo"]("L,i,j")  = l2["abab"]("L,i,k,a,b") * t2["abab"]("a,b,j,k");

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["205_aa_Loo"]("L,i,j")  = l2_1["abab"]("L,i,k,a,b") * t2["abab"]("a,b,j,k");

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["204_aa_Loo"]("L,i,j")  = l2_1["aaaa"]("L,k,i,a,b") * t2["aaaa"]("a,b,k,j");

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["203_aa_Loo"]("L,i,j")  = l2["aaaa"]("L,i,k,a,b") * t2["aaaa"]("a,b,k,j");

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["202_aa_Loo"]("L,i,j")  = l2_1["aaaa"]("L,i,k,a,b") * t2_1["aaaa"]("a,b,k,j");

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["201_aa_Loo"]("L,i,j")  = l2_1["aaaa"]("L,i,k,a,b") * t2["aaaa"]("a,b,k,j");

    // flops: o0v2L1  = o2v3L1
    //  mems: o0v2L1  = o0v2L1
    tmps_["200_bb_Lvv"]("L,a,b")  = l2_1["bbbb"]("L,i,j,a,c") * t2["bbbb"]("b,c,i,j");

    // flops: o0v2L1  = o2v3L1
    //  mems: o0v2L1  = o0v2L1
    tmps_["199_bb_Lvv"]("L,a,b")  = l2_1["abab"]("L,i,j,c,a") * t2["abab"]("c,b,i,j");

    // flops: o0v2L1  = o2v3L1
    //  mems: o0v2L1  = o0v2L1
    tmps_["198_aa_Lvv"]("L,a,b")  = l2_1["aaaa"]("L,i,j,a,c") * t2["aaaa"]("b,c,i,j");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["196_aa_Lov"]("L,i,a")  = l2_1["aaaa"]("L,j,i,b,a") * t1_1["aa"]("b,j");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["195_aa_Lov"]("L,i,a")  = l2_1["aaaa"]("L,i,j,a,b") * t1_1["aa"]("b,j");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["194_bb_Lvo"]("L,a,i")  = l1_1["aa"]("L,j,b") * t2["abab"]("b,a,j,i");

    // flops: o1v1L1  = o2v1L1 o1v1L1 o2v1L1 o2v1L1 o1v1L1 o1v1L1 o3v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o3v1L1 o1v1L1 o1v3L1 o1v1L1 o2v1L1 o1v1L1 o1v3L1 o1v1L1 o3v1L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o1v3L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o3v1L1 o1v1L1 o3v2L1 o2v3L1 o2v1L1 o2v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o1v2L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o2v3L1 o1v1L1 o2v1L1 o1v1L1 o3v2L1 o1v1L1 o1v2L1 o1v1L1 o1v1L1 o1v1L1 o1v2L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o1v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o2v3L1 o1v1L1 o1v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o2v3L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o1v2L1 o1v1L1 o1v1L1 o1v1L1 o1v3L1 o1v1L1 o2v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    tmps_["211_aa_Lvo"]("L,a,i")  = -0.50 * f["aa_ov"]("j,a") * tmps_["201_aa_Loo"]("L,i,j");
    tmps_["211_aa_Lvo"]("L,a,i") += 0.25 * dp["aa_ov"]("i,a") * tmps_["103_L"]("L");
    tmps_["211_aa_Lvo"]("L,a,i") += 0.50 * reused_["13_aa_ov"]("j,a") * tmps_["201_aa_Loo"]("L,i,j");
    tmps_["211_aa_Lvo"]("L,a,i") += f["aa_ov"]("j,a") * tmps_["205_aa_Loo"]("L,i,j");
    tmps_["211_aa_Lvo"]("L,a,i") += 0.50 * eri["abab_oovo"]("i,o,a,l") * tmps_["208_bb_Loo"]("L,l,o");
    tmps_["211_aa_Lvo"]("L,a,i") -= reused_["13_aa_ov"]("j,a") * tmps_["205_aa_Loo"]("L,i,j");
    tmps_["211_aa_Lvo"]("L,a,i") -= eri["abab_oovv"]("i,l,a,b") * tmps_["194_bb_Lvo"]("L,b,l");
    tmps_["211_aa_Lvo"]("L,a,i") -= 2.00 * dp["aa_ov"]("k,a") * tmps_["209_aa_Loo"]("L,i,k");
    tmps_["211_aa_Lvo"]("L,a,i") += eri["abab_oovo"]("i,o,a,l") * tmps_["178_bb_Loo"]("L,l,o");
    tmps_["211_aa_Lvo"]("L,a,i") += 0.50 * eri["baab_vovv"]("b,i,a,c") * tmps_["200_bb_Lvv"]("L,b,c");
    tmps_["211_aa_Lvo"]("L,a,i") -= dp["aa_ov"]("j,a") * tmps_["206_aa_Loo"]("L,i,j");
    tmps_["211_aa_Lvo"]("L,a,i") += eri["aaaa_vovv"]("f,i,a,g") * tmps_["99_aa_Lvv"]("L,f,g");
    tmps_["211_aa_Lvo"]("L,a,i") += eri["aaaa_oovo"]("i,n,a,j") * tmps_["205_aa_Loo"]("L,j,n");
    tmps_["211_aa_Lvo"]("L,a,i") += t0_1 * tmps_["100_aa_Lov"]("L,i,a");
    tmps_["211_aa_Lvo"]("L,a,i") += dp["aa_ov"]("j,a") * tmps_["202_aa_Loo"]("L,i,j");
    tmps_["211_aa_Lvo"]("L,a,i") += 0.50 * eri["aaaa_vovv"]("f,i,a,g") * tmps_["198_aa_Lvv"]("L,f,g");
    tmps_["211_aa_Lvo"]("L,a,i") -= t0_1 * tmps_["197_aa_Lov"]("L,i,a");
    tmps_["211_aa_Lvo"]("L,a,i") -= 2.00 * dp["aa_ov"]("j,a") * tmps_["207_aa_Loo"]("L,i,j");
    tmps_["211_aa_Lvo"]("L,a,i") += 0.50 * eri["aaaa_oovo"]("i,n,a,j") * tmps_["204_aa_Loo"]("L,j,n");
    tmps_["211_aa_Lvo"]("L,a,i") -= 0.50 * l2_1["aaaa"]("L,k,j,d,a") * reused_["130_aaaa_vooo"]("d,k,j,i");
    tmps_["211_aa_Lvo"]("L,a,i") -= eri["abab_vvvo"]("f,e,a,m") * l2_1["abab"]("L,i,m,f,e");
    tmps_["211_aa_Lvo"]("L,a,i") -= l1["aa"]("L,k,a") * dp["aa_oo"]("i,k");
    tmps_["211_aa_Lvo"]("L,a,i") += 0.50 * l1_1["aa"]("L,k,a") * reused_["174_aa_oo"]("i,k");
    tmps_["211_aa_Lvo"]("L,a,i") += 2.00 * scalars_["3"] * l1_1["aa"]("L,i,a");
    tmps_["211_aa_Lvo"]("L,a,i") += scalars_["2"] * l1["aa"]("L,i,a");
    tmps_["211_aa_Lvo"]("L,a,i") -= l2_1["abab"]("L,i,m,d,b") * reused_["168_abba_vovv"]("d,m,b,a");
    tmps_["211_aa_Lvo"]("L,a,i") -= l2_1["aaaa"]("L,k,j,d,a") * reused_["150_aaaa_oovo"]("i,j,d,k");
    tmps_["211_aa_Lvo"]("L,a,i") += dp["aa_vv"]("f,a") * tmps_["195_aa_Lov"]("L,i,f");
    tmps_["211_aa_Lvo"]("L,a,i") -= dp["aa_vo"]("d,k") * l2["aaaa"]("L,i,k,d,a");
    tmps_["211_aa_Lvo"]("L,a,i") += 0.50 * l2_1["aaaa"]("L,i,k,d,a") * reused_["87_aa_vo"]("d,k");
    tmps_["211_aa_Lvo"]("L,a,i") += eri["baab_vovo"]("e,i,a,m") * l1_1["bb"]("L,m,e");
    tmps_["211_aa_Lvo"]("L,a,i") += l2_1["abab"]("L,j,m,a,e") * reused_["173_bbaa_vooo"]("e,m,i,j");
    tmps_["211_aa_Lvo"]("L,a,i") += reused_["167_aabb_vvvo"]("f,a,e,m") * l2_1["abab"]("L,i,m,f,e");
    tmps_["211_aa_Lvo"]("L,a,i") += l1_1["aa"]("L,k,a") * reused_["16_aa_oo"]("i,k");
    tmps_["211_aa_Lvo"]("L,a,i") += l2_1["aaaa"]("L,k,j,d,a") * reused_["170_aaaa_vooo"]("d,k,i,j");
    tmps_["211_aa_Lvo"]("L,a,i") += 0.50 * reused_["18_aa_vv"]("a,d") * l1_1["aa"]("L,i,d");
    tmps_["211_aa_Lvo"]("L,a,i") -= w0 * l1_1["aa"]("L,i,a");
    tmps_["211_aa_Lvo"]("L,a,i") += reused_["59_aa_vv"]("d,a") * l1_1["aa"]("L,i,d");
    tmps_["211_aa_Lvo"]("L,a,i") -= tmps_["196_aa_Lov"]("L,j,a") * dp["aa_oo"]("i,j");
    tmps_["211_aa_Lvo"]("L,a,i") += f["aa_oo"]("i,k") * l1_1["aa"]("L,k,a");
    tmps_["211_aa_Lvo"]("L,a,i") -= l1_1["aa"]("L,k,a") * reused_["17_aa_oo"]("i,k");
    tmps_["211_aa_Lvo"]("L,a,i") += l0_1("L") * reused_["13_aa_ov"]("i,a");
    tmps_["211_aa_Lvo"]("L,a,i") -= f["aa_ov"]("i,a") * l0_1("L");
    tmps_["211_aa_Lvo"]("L,a,i") += 2.00 * l2_1["abab"]("L,i,m,a,e") * reused_["88_bb_vo"]("e,m");
    tmps_["211_aa_Lvo"]("L,a,i") += scalars_["1"] * l1["aa"]("L,i,a");
    tmps_["211_aa_Lvo"]("L,a,i") -= l2_1["aaaa"]("L,i,k,d,a") * reused_["90_aa_vo"]("d,k");
    tmps_["211_aa_Lvo"]("L,a,i") -= l2_1["abab"]("L,i,m,f,e") * reused_["169_abab_vovv"]("a,m,f,e");
    tmps_["211_aa_Lvo"]("L,a,i") -= eri["baab_vooo"]("e,i,k,l") * l2_1["abab"]("L,k,l,a,e");
    tmps_["211_aa_Lvo"]("L,a,i") += 0.50 * eri["aaaa_vooo"]("d,i,k,j") * l2_1["aaaa"]("L,k,j,d,a");
    tmps_["211_aa_Lvo"]("L,a,i") += l1_1["aa"]("L,i,d") * reused_["57_aa_vv"]("a,d");
    tmps_["211_aa_Lvo"]("L,a,i") += l1_1["aa"]("L,k,d") * reused_["69_aaaa_voov"]("d,k,i,a");
    tmps_["211_aa_Lvo"]("L,a,i") -= l2_1["aaaa"]("L,i,k,d,a") * reused_["152_aa_vo"]("d,k");
    tmps_["211_aa_Lvo"]("L,a,i") += dp["aa_vv"]("d,a") * l1["aa"]("L,i,d");
    tmps_["211_aa_Lvo"]("L,a,i") -= l1_1["bb"]("L,m,e") * reused_["158_bbaa_voov"]("e,m,i,a");
    tmps_["211_aa_Lvo"]("L,a,i") -= l2_1["aaaa"]("L,i,k,f,d") * reused_["145_aaaa_vvvo"]("f,a,d,k");
    tmps_["211_aa_Lvo"]("L,a,i") -= l2_1["abab"]("L,k,l,a,e") * reused_["171_abba_oovo"]("i,l,e,k");
    tmps_["211_aa_Lvo"]("L,a,i") += 0.50 * eri["aaaa_vvvo"]("d,f,a,k") * l2_1["aaaa"]("L,i,k,f,d");
    tmps_["211_aa_Lvo"]("L,a,i") -= f["aa_vv"]("d,a") * l1_1["aa"]("L,i,d");
    tmps_["211_aa_Lvo"]("L,a,i") -= f["bb_vo"]("e,m") * l2_1["abab"]("L,i,m,a,e");
    tmps_["211_aa_Lvo"]("L,a,i") += f["aa_vo"]("d,k") * l2_1["aaaa"]("L,i,k,d,a");
    tmps_["211_aa_Lvo"]("L,a,i") += l2_1["abab"]("L,i,m,a,e") * reused_["84_bb_vo"]("e,m");
    tmps_["211_aa_Lvo"]("L,a,i") += l2_1["abab"]("L,i,m,a,e") * reused_["91_bb_vo"]("e,m");
    tmps_["211_aa_Lvo"]("L,a,i") -= l2_1["abab"]("L,k,l,a,e") * reused_["175_baab_vooo"]("e,i,k,l");
    tmps_["211_aa_Lvo"]("L,a,i") -= l2["aaaa"]("L,i,k,d,a") * reused_["58_aa_vo"]("d,k");
    tmps_["211_aa_Lvo"]("L,a,i") += scalars_["5"] * l1_1["aa"]("L,i,a");
    tmps_["211_aa_Lvo"]("L,a,i") += 0.25 * l2_1["aaaa"]("L,i,k,f,d") * reused_["14_aaaa_vovv"]("a,k,f,d");
    tmps_["211_aa_Lvo"]("L,a,i") += 2.00 * l2_1["aaaa"]("L,i,k,d,a") * reused_["89_aa_ov"]("k,d");
    tmps_["211_aa_Lvo"]("L,a,i") += l1_1["bb"]("L,m,e") * reused_["1_bbaa_voov"]("e,m,i,a");
    tmps_["211_aa_Lvo"]("L,a,i") -= l2_1["abab"]("L,j,m,a,e") * reused_["172_bbaa_vooo"]("e,m,i,j");
    tmps_["211_aa_Lvo"]("L,a,i") += l2_1["aaaa"]("L,i,k,f,d") * reused_["165_aaaa_vovv"]("d,k,f,a");
    tmps_["211_aa_Lvo"]("L,a,i") += eri["aaaa_vovo"]("d,i,a,k") * l1_1["aa"]("L,k,d");
    tmps_["211_aa_Lvo"]("L,a,i") += dp["bb_vo"]("e,m") * l2["abab"]("L,i,m,a,e");
    tmps_["211_aa_Lvo"]("L,a,i") += 2.00 * scalars_["4"] * l1_1["aa"]("L,i,a");
    tmps_["211_aa_Lvo"]("L,a,i") -= l2_1["abab"]("L,i,m,f,e") * reused_["166_aabb_vvvo"]("f,a,e,m");
    tmps_["211_aa_Lvo"]("L,a,i") += 0.50 * dp["aa_ov"]("j,a") * tmps_["203_aa_Loo"]("L,i,j");
    tmps_["211_aa_Lvo"]("L,a,i") += dp["aa_vv"]("f,a") * tmps_["102_aa_Lov"]("L,i,f");
    tmps_["211_aa_Lvo"]("L,a,i") += dp["aa_ov"]("i,a") * tmps_["210_L"]("L");
    tmps_["211_aa_Lvo"]("L,a,i") += eri["baab_vovv"]("b,i,a,c") * tmps_["199_bb_Lvv"]("L,b,c");
    tmps_["211_aa_Lvo"]("L,a,i") -= dp["aa_oo"]("i,j") * tmps_["102_aa_Lov"]("L,j,a");
    tmps_["197_aa_Lov"].~TArrayD();
    tmps_["196_aa_Lov"].~TArrayD();
    tmps_["195_aa_Lov"].~TArrayD();
    tmps_["194_bb_Lvo"].~TArrayD();
    tmps_["100_aa_Lov"].~TArrayD();

    // sigmal1_1_aa  = -0.50 f_aa(k,a) t2_abab(c,b,k,j) l2_1_abab(i,j,c,b)
    //              += -0.50 f_aa(k,a) t2_abab(b,c,k,j) l2_1_abab(i,j,b,c)
    //              += -0.50 d-_aa(k,a) t2_aaaa(c,b,j,k) t0_1 l2_1_aaaa(i,j,c,b)
    //              += -0.25 d-_aa(i,a) t2_1_bbbb(c,b,j,k) l2_1_bbbb(j,k,c,b)
    //              += -0.25 d-_aa(i,a) t2_1_aaaa(c,b,j,k) l2_1_aaaa(j,k,c,b)
    //              += -0.25 d-_aa(i,a) t2_1_abab(c,b,j,k) l2_1_abab(j,k,c,b)
    //              += -0.25 d-_aa(i,a) t2_1_abab(c,b,k,j) l2_1_abab(k,j,c,b)
    //              += -0.25 d-_aa(i,a) t2_1_abab(b,c,j,k) l2_1_abab(j,k,b,c)
    //              += -0.25 d-_aa(i,a) t2_1_abab(b,c,k,j) l2_1_abab(k,j,b,c)
    //              += -0.50 <i,l||a,k>_abab t2_bbbb(c,b,j,l) l2_1_bbbb(j,k,c,b)
    //              += +0.50 f_aa(k,a) t2_aaaa(c,b,j,k) l2_1_aaaa(i,j,c,b)
    //              += +0.50 d-_aa(k,a) t2_abab(c,b,k,j) t0_1 l2_1_abab(i,j,c,b)
    //              += +0.50 d-_aa(k,a) t2_abab(b,c,k,j) t0_1 l2_1_abab(i,j,b,c)
    //              += +1.00 <i,k||a,c>_abab t2_abab(b,c,j,k) l1_1_aa(j,b)
    //              += +2.00 d-_aa(j,a) t1_1_aa(b,j) l1_1_aa(i,b)
    //              += -0.50 <i,l||a,k>_abab t2_abab(c,b,j,l) l2_1_abab(j,k,c,b)
    //              += -0.50 <i,l||a,k>_abab t2_abab(b,c,j,l) l2_1_abab(j,k,b,c)
    //              += +0.50 <i,c||a,d>_abab t2_bbbb(d,b,j,k) l2_1_bbbb(j,k,c,b)
    //              += +0.50 d-_aa(k,a) t2_abab(c,b,k,j) l2_abab(i,j,c,b)
    //              += +0.50 d-_aa(k,a) t2_abab(b,c,k,j) l2_abab(i,j,b,c)
    //              += -0.50 <i,c||d,a>_aaaa t2_abab(d,b,j,k) l2_1_abab(j,k,c,b)
    //              += -0.50 <i,c||d,a>_aaaa t2_abab(d,b,k,j) l2_1_abab(k,j,c,b)
    //              += -0.50 <i,l||a,k>_aaaa t2_abab(c,b,l,j) l2_1_abab(k,j,c,b)
    //              += -0.50 <i,l||a,k>_aaaa t2_abab(b,c,l,j) l2_1_abab(k,j,b,c)
    //              += -1.00 d-_aa(k,c) t2_aaaa(c,b,j,k) t0_1 l2_1_aaaa(i,j,b,a)
    //              += -1.00 d-_aa(k,a) t2_1_aaaa(c,b,j,k) l2_1_aaaa(i,j,c,b)
    //              += -0.50 <i,c||d,a>_aaaa t2_aaaa(d,b,j,k) l2_1_aaaa(j,k,c,b)
    //              += +1.00 d-_bb(k,c) t2_bbbb(c,b,j,k) t0_1 l2_1_abab(i,j,a,b)
    //              += -0.50 d-_aa(i,c) t2_aaaa(c,b,j,k) t0_1 l2_1_aaaa(j,k,b,a)
    //              += -1.00 d-_aa(k,c) t2_abab(c,b,k,j) t0_1 l2_1_abab(i,j,a,b)
    //              += +0.50 d-_aa(i,c) t2_abab(c,b,j,k) t0_1 l2_1_abab(j,k,a,b)
    //              += +0.50 d-_aa(i,c) t2_abab(c,b,k,j) t0_1 l2_1_abab(k,j,a,b)
    //              += +1.00 d-_aa(k,a) t2_1_abab(c,b,k,j) l2_1_abab(i,j,c,b)
    //              += +1.00 d-_aa(k,a) t2_1_abab(b,c,k,j) l2_1_abab(i,j,b,c)
    //              += -0.50 <i,l||a,k>_aaaa t2_aaaa(c,b,j,l) l2_1_aaaa(j,k,c,b)
    //              += -1.00 d-_bb(j,j) l1_aa(i,a)
    //              += -2.00 d-_bb(j,b) t1_1_bb(b,j) l1_1_aa(i,a)
    //              += -0.50 <i,k||b,c>_aaaa t2_aaaa(b,c,j,k) l1_1_aa(j,a)
    //              += +1.00 d-_aa(i,j) l1_aa(j,a)
    //              += -1.00 <k,c||a,d>_abab t2_abab(b,d,k,j) l2_1_abab(i,j,b,c)
    //              += -1.00 <i,l||k,c>_abab t2_abab(b,c,j,l) l2_1_aaaa(j,k,b,a)
    //              += -1.00 d-_aa(c,a) t1_1_aa(b,j) l2_1_aaaa(i,j,c,b)
    //              += +0.50 <c,b||a,j>_abab l2_1_abab(i,j,c,b)
    //              += +0.50 <b,c||a,j>_abab l2_1_abab(i,j,b,c)
    //              += +1.00 d-_aa(b,j) l2_aaaa(i,j,b,a)
    //              += +0.50 <k,b||c,d>_aaaa t2_aaaa(c,d,j,k) l2_1_aaaa(i,j,b,a)
    //              += +0.50 <l,k||j,c>_abab t2_abab(b,c,l,k) l2_1_aaaa(i,j,b,a)
    //              += +0.50 <k,l||j,c>_abab t2_abab(b,c,k,l) l2_1_aaaa(i,j,b,a)
    //              += +0.50 <l,k||c,j>_aaaa t2_aaaa(c,b,l,k) l2_1_aaaa(i,j,b,a)
    //              += +1.00 f_aa(k,c) t2_aaaa(c,b,j,k) l2_1_aaaa(i,j,b,a)
    //              += -0.50 <b,k||c,d>_abab t2_abab(c,d,j,k) l2_1_aaaa(i,j,b,a)
    //              += -0.50 <b,k||d,c>_abab t2_abab(d,c,j,k) l2_1_aaaa(i,j,b,a)
    //              += +1.00 d-_aa(b,j) t0_1 l2_1_aaaa(i,j,b,a)
    //              += -1.00 f_bb(k,c) t2_abab(b,c,j,k) l2_1_aaaa(i,j,b,a)
    //              += +1.00 <i,b||a,j>_abab l1_1_bb(j,b)
    //              += +1.00 <i,l||k,c>_abab t2_bbbb(c,b,j,l) l2_1_abab(k,j,a,b)
    //              += -1.00 <c,k||a,d>_abab t2_bbbb(d,b,j,k) l2_1_abab(i,j,c,b)
    //              += -0.50 <i,k||b,c>_abab t2_abab(b,c,j,k) l1_1_aa(j,a)
    //              += -0.50 <i,k||c,b>_abab t2_abab(c,b,j,k) l1_1_aa(j,a)
    //              += -1.00 <i,l||c,k>_aaaa t2_aaaa(c,b,j,l) l2_1_aaaa(j,k,b,a)
    //              += -0.50 <k,j||c,a>_aaaa t2_aaaa(c,b,k,j) l1_1_aa(i,b)
    //              += +1.00 l1_1_aa(i,a) w0
    //              += -1.00 d-_aa(b,a) t0_1 l1_1_aa(i,b)
    //              += +1.00 d-_aa(i,k) t1_1_aa(b,j) l2_1_aaaa(j,k,b,a)
    //              += -1.00 f_aa(i,j) l1_1_aa(j,a)
    //              += +1.00 d-_aa(i,j) t0_1 l1_1_aa(j,a)
    //              += -1.00 d-_aa(i,a) t0_1 l0_1
    //              += +1.00 f_aa(i,a) l0_1
    //              += -2.00 d-_bb(b,c) t1_1_bb(c,j) l2_1_abab(i,j,a,b)
    //              += +2.00 d-_bb(k,j) t1_1_bb(b,k) l2_1_abab(i,j,a,b)
    //              += -2.00 d-_aa(k,c) t2_1_abab(c,b,k,j) l2_1_abab(i,j,a,b)
    //              += +2.00 d-_bb(k,c) t2_1_bbbb(c,b,j,k) l2_1_abab(i,j,a,b)
    //              += -1.00 d-_aa(j,j) l1_aa(i,a)
    //              += +1.00 d-_bb(k,k) t1_1_aa(b,j) l2_1_aaaa(i,j,b,a)
    //              += +1.00 d-_aa(k,k) t1_1_aa(b,j) l2_1_aaaa(i,j,b,a)
    //              += +0.25 <l,k||a,j>_abab t2_abab(c,b,l,k) l2_1_abab(i,j,c,b)
    //              += +0.25 <k,l||a,j>_abab t2_abab(c,b,k,l) l2_1_abab(i,j,c,b)
    //              += +0.25 <l,k||a,j>_abab t2_abab(b,c,l,k) l2_1_abab(i,j,b,c)
    //              += +0.25 <k,l||a,j>_abab t2_abab(b,c,k,l) l2_1_abab(i,j,b,c)
    //              += -0.50 <i,b||j,k>_abab l2_1_abab(j,k,a,b)
    //              += -0.50 <i,b||k,j>_abab l2_1_abab(k,j,a,b)
    //              += +0.50 <i,b||j,k>_aaaa l2_1_aaaa(j,k,b,a)
    //              += -0.50 <k,j||a,c>_abab t2_abab(b,c,k,j) l1_1_aa(i,b)
    //              += -0.50 <j,k||a,c>_abab t2_abab(b,c,j,k) l1_1_aa(i,b)
    //              += +1.00 <i,k||c,a>_aaaa t2_aaaa(c,b,j,k) l1_1_aa(j,b)
    //              += +1.00 d-_bb(k,c) t2_abab(b,c,j,k) t0_1 l2_1_aaaa(i,j,b,a)
    //              += -1.00 d-_aa(b,a) l1_aa(i,b)
    //              += -1.00 <i,k||c,a>_aaaa t2_abab(c,b,k,j) l1_1_bb(j,b)
    //              += +1.00 <c,k||a,d>_abab t2_abab(b,d,j,k) l2_1_aaaa(i,j,c,b)
    //              += +1.00 <i,l||c,k>_abab t2_abab(c,b,j,l) l2_1_abab(j,k,a,b)
    //              += +0.50 <c,b||a,j>_aaaa l2_1_aaaa(i,j,c,b)
    //              += +1.00 f_aa(b,a) l1_1_aa(i,b)
    //              += +1.00 f_bb(b,j) l2_1_abab(i,j,a,b)
    //              += -1.00 f_aa(b,j) l2_1_aaaa(i,j,b,a)
    //              += -1.00 f_bb(k,c) t2_bbbb(c,b,j,k) l2_1_abab(i,j,a,b)
    //              += -1.00 d-_bb(b,j) t0_1 l2_1_abab(i,j,a,b)
    //              += -0.50 <l,k||c,j>_abab t2_abab(c,b,l,k) l2_1_abab(i,j,a,b)
    //              += -0.50 <k,l||c,j>_abab t2_abab(c,b,k,l) l2_1_abab(i,j,a,b)
    //              += -0.50 <l,k||c,j>_bbbb t2_bbbb(c,b,l,k) l2_1_abab(i,j,a,b)
    //              += +1.00 f_aa(k,c) t2_abab(c,b,k,j) l2_1_abab(i,j,a,b)
    //              += +0.50 <k,b||c,d>_abab t2_abab(c,d,k,j) l2_1_abab(i,j,a,b)
    //              += +0.50 <k,b||d,c>_abab t2_abab(d,c,k,j) l2_1_abab(i,j,a,b)
    //              += -0.50 <k,b||c,d>_bbbb t2_bbbb(c,d,j,k) l2_1_abab(i,j,a,b)
    //              += -1.00 d-_aa(k,k) t1_1_bb(b,j) l2_1_abab(i,j,a,b)
    //              += -1.00 d-_bb(k,k) t1_1_bb(b,j) l2_1_abab(i,j,a,b)
    //              += -0.25 <i,b||c,d>_abab t2_abab(c,d,j,k) l2_1_abab(j,k,a,b)
    //              += -0.25 <i,b||c,d>_abab t2_abab(c,d,k,j) l2_1_abab(k,j,a,b)
    //              += -0.25 <i,b||d,c>_abab t2_abab(d,c,j,k) l2_1_abab(j,k,a,b)
    //              += -0.25 <i,b||d,c>_abab t2_abab(d,c,k,j) l2_1_abab(k,j,a,b)
    //              += -0.50 f_aa(i,c) t2_abab(c,b,j,k) l2_1_abab(j,k,a,b)
    //              += -0.50 f_aa(i,c) t2_abab(c,b,k,j) l2_1_abab(k,j,a,b)
    //              += +1.00 d-_bb(k,c) t2_abab(b,c,j,k) l2_aaaa(i,j,b,a)
    //              += -1.00 d-_aa(j,j) t0_1 l1_1_aa(i,a)
    //              += +1.00 f_aa(k,k) r2_bbbb(a,b,i,j)
    //              += -0.50 <l,k||l,k>_abab r2_bbbb(a,b,i,j)
    //              += -0.50 <k,l||k,l>_abab r2_bbbb(a,b,i,j)
    //              += +0.25 <l,k||c,d>_aaaa r2_bbbb(a,b,i,j) t2_aaaa(c,d,l,k)
    //              += -0.50 <l,k||l,k>_bbbb r2_bbbb(a,b,i,j)
    //              += -0.50 <l,k||l,k>_aaaa r2_bbbb(a,b,i,j)
    //              += +0.25 <l,k||c,d>_bbbb r2_bbbb(a,b,i,j) t2_bbbb(c,d,l,k)
    //              += +1.00 f_bb(k,k) r2_bbbb(a,b,i,j)
    //              += -1.00 d-_bb(k,k) r2_bbbb(a,b,i,j) t0_1
    //              += +0.25 <l,k||c,d>_abab r2_bbbb(a,b,i,j) t2_abab(c,d,l,k)
    //              += +0.25 <k,l||c,d>_abab r2_bbbb(a,b,i,j) t2_abab(c,d,k,l)
    //              += +0.25 <l,k||d,c>_abab r2_bbbb(a,b,i,j) t2_abab(d,c,l,k)
    //              += +0.25 <k,l||d,c>_abab r2_bbbb(a,b,i,j) t2_abab(d,c,k,l)
    //              += +0.25 <l,k||a,j>_aaaa t2_aaaa(c,b,l,k) l2_1_aaaa(i,j,c,b)
    //              += -2.00 d-_aa(k,j) t1_1_aa(b,k) l2_1_aaaa(i,j,b,a)
    //              += -2.00 d-_aa(k,c) t2_1_aaaa(c,b,j,k) l2_1_aaaa(i,j,b,a)
    //              += +2.00 d-_aa(b,c) t1_1_aa(c,j) l2_1_aaaa(i,j,b,a)
    //              += +2.00 d-_bb(k,c) t2_1_abab(b,c,j,k) l2_1_aaaa(i,j,b,a)
    //              += -1.00 <i,k||a,c>_abab t2_bbbb(c,b,j,k) l1_1_bb(j,b)
    //              += +1.00 <i,l||c,k>_aaaa t2_abab(c,b,l,j) l2_1_abab(k,j,a,b)
    //              += +0.50 f_aa(i,c) t2_aaaa(c,b,j,k) l2_1_aaaa(j,k,b,a)
    //              += +0.25 <i,b||c,d>_aaaa t2_aaaa(c,d,j,k) l2_1_aaaa(j,k,b,a)
    //              += -1.00 <k,c||d,a>_aaaa t2_aaaa(d,b,j,k) l2_1_aaaa(i,j,c,b)
    //              += +1.00 <i,b||a,j>_aaaa l1_1_aa(j,b)
    //              += -1.00 d-_bb(b,j) l2_abab(i,j,a,b)
    //              += -2.00 d-_aa(j,b) t1_1_aa(b,j) l1_1_aa(i,a)
    //              += +1.00 <k,c||d,a>_aaaa t2_abab(d,b,k,j) l2_1_abab(i,j,c,b)
    //              += -0.50 d-_aa(k,a) t2_aaaa(c,b,j,k) l2_aaaa(i,j,c,b)
    //              += -1.00 d-_aa(c,a) t1_1_bb(b,j) l2_1_abab(i,j,c,b)
    //              += -1.00 d-_aa(i,a) t1_1_bb(b,j) l1_1_bb(j,b)
    //              += -1.00 d-_aa(i,a) t1_1_aa(b,j) l1_1_aa(j,b)
    //              += +0.50 <i,c||a,d>_abab t2_abab(b,d,j,k) l2_1_abab(j,k,b,c)
    //              += +0.50 <i,c||a,d>_abab t2_abab(b,d,k,j) l2_1_abab(k,j,b,c)
    //              += +1.00 d-_aa(i,k) t1_1_bb(b,j) l2_1_abab(k,j,a,b)
    sigmal1_1_aa("L,a,i")  = -1.00 * tmps_["211_aa_Lvo"]("L,a,i");
    tmps_["211_aa_Lvo"].~TArrayD();

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["212_bb_Lvo"]("R,a,i")  = -1.00 * dp["bb_ov"]("j,b") * r2["bbbb"]("R,b,a,i,j");
    tmps_["212_bb_Lvo"].~TArrayD();

    // flops: o1v1L1  = o1v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["213_bb_Lvo"]("R,a,i")  = dp["bb_vv"]("a,b") * r1["bb"]("R,b,i");

    // sigmar1_1_bb += -1.00 d+_bb(a,b) r1_bb(b,i)
    sigmar1_1_bb("R,a,i") -= tmps_["213_bb_Lvo"]("R,a,i");

    // flops: o1v1L1  = o2v1L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["214_bb_Lov"]("R,i,a")  = dp["bb_oo"]("j,i") * r1["bb"]("R,a,j");

    // sigmar1_1_bb += +1.00 d+_bb(j,i) r1_bb(a,j)
    sigmar1_1_bb("R,a,i") += tmps_["214_bb_Lov"]("R,i,a");

    // flops: o1v1L1  = o2v2L1 o2v2L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1
    tmps_["215_bb_Lvo"]("R,a,i")  = -1.00 * dp["bb_ov"]("k,c") * r2["bbbb"]("R,c,a,i,k");
    tmps_["215_bb_Lvo"]("R,a,i") += dp["aa_ov"]("j,b") * r2["abab"]("R,b,a,j,i");

    // sigmar1_1_bb += -1.00 d+_aa(j,b) r2_abab(b,a,j,i)
    //              += +1.00 d+_bb(j,b) r2_bbbb(b,a,i,j)
    sigmar1_1_bb("R,a,i") -= tmps_["215_bb_Lvo"]("R,a,i");
}
