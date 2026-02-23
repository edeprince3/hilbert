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

void hilbert::EOM_EE_QED_CCSD_21::sigma_ee_21_5() {

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


    // flops: o2v2L1  = o3v3L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j")  = -1.00 * eri["bbbb_vovo"]("a,k,c,i") * r2["bbbb"]("R,c,b,j,k");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") -= r2["bbbb"]("R,d,b,j,m") * reused_["28_bbbb_voov"]("a,i,m,d");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") += dp["bb_vo"]("a,i") * r1_1["bb"]("R,b,j");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") += r1["bb"]("R,b,j") * reused_["91_bb_vo"]("a,i");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") -= r1["bb"]("R,b,m") * reused_["178_bbbb_oovo"]("m,i,a,j");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") += r1["bb"]("R,d,j") * reused_["177_bbbb_vvvo"]("a,d,b,i");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") += r1["bb"]("R,b,j") * reused_["84_bb_vo"]("a,i");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") += eri["baab_vovo"]("a,l,e,i") * r2["abab"]("R,e,b,l,j");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") -= r1_1["bb"]("R,b,j") * reused_["33_bb_vo"]("a,i");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") += r1["bb"]("R,b,m") * reused_["179_bbbb_oovo"]("m,i,a,j");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") -= r2["abab"]("R,f,b,n,j") * reused_["27_bbaa_voov"]("a,i,n,f");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") += r2["abab"]("R,f,b,n,j") * reused_["1_bbaa_voov"]("a,i,n,f");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") += r2["bbbb"]("R,d,b,j,m") * reused_["26_bbbb_voov"]("a,i,m,d");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") += r1_1["bb"]("R,b,j") * reused_["32_bb_vo"]("a,i");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") += r1["bb"]("R,b,j") * reused_["88_bb_vo"]("a,i");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") -= f["bb_vo"]("a,i") * r1["bb"]("R,b,j");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") -= r1["bb"]("R,d,j") * reused_["176_bbbb_vvvo"]("a,d,b,i");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") += reused_["180_bb_vo"]("a,i") * r1["bb"]("R,b,j");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") += t1_1["bb"]("a,i") * tmps_["215_bb_Lvo"]("R,b,j");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") += t1_1["bb"]("a,j") * tmps_["214_bb_Lov"]("R,i,b");
    tmps_["216_bbbb_Lvovo"]("R,a,i,b,j") += tmps_["213_bb_Lvo"]("R,a,j") * reused_["234_bb_vo"]("b,i");

    // sigmar2_bbbb += +1.00 P(i,j) P(a,b) d-_aa(k,c) r2_abab(c,b,k,i) t1_1_bb(a,j)
    //              += -1.00 P(i,j) P(a,b) d-_bb(k,c) r2_bbbb(c,b,i,k) t1_1_bb(a,j)
    //              += +1.00 P(i,j) P(a,b) d-_aa(k,c) r1_bb(b,i) t2_abab(c,a,k,j) t0_1
    //              += -1.00 P(i,j) P(a,b) d-_bb(k,c) r1_bb(b,i) t2_bbbb(c,a,j,k) t0_1
    //              += +1.00 P(i,j) P(a,b) d-_bb(k,j) r1_bb(b,k) t1_1_bb(a,i)
    //              += +1.00 P(i,j) P(a,b) d-_bb(a,j) r1_1_bb(b,i)
    //              += +1.00 P(i,j) P(a,b) <l,k||c,d>_bbbb r2_bbbb(d,b,i,l) t2_bbbb(c,a,j,k)
    //              += +1.00 P(i,j) P(a,b) <k,a||c,j>_bbbb r2_bbbb(c,b,i,k)
    //              += +1.00 P(i,j) P(a,b) d-_aa(k,k) r1_bb(b,i) t1_1_bb(a,j)
    //              += +1.00 P(i,j) P(a,b) d-_bb(k,k) r1_bb(b,i) t1_1_bb(a,j)
    //              += -1.00 P(i,j) P(a,b) <k,l||c,j>_abab r1_bb(b,l) t2_abab(c,a,k,i)
    //              += -1.00 P(i,j) P(a,b) <k,a||c,d>_bbbb r1_bb(d,i) t2_bbbb(c,b,j,k)
    //              += +1.00 P(i,j) P(a,b) f_bb(k,c) r1_bb(b,i) t2_bbbb(c,a,j,k)
    //              += +1.00 P(i,j) P(a,b) d-_bb(a,j) r1_bb(b,i) t0_1
    //              += +0.50 P(i,j) P(a,b) <l,k||c,j>_abab r1_bb(b,i) t2_abab(c,a,l,k)
    //              += +0.50 P(i,j) P(a,b) <k,l||c,j>_abab r1_bb(b,i) t2_abab(c,a,k,l)
    //              += +0.50 P(i,j) P(a,b) <l,k||c,j>_bbbb r1_bb(b,i) t2_bbbb(c,a,l,k)
    //              += -1.00 P(i,j) P(a,b) f_aa(k,c) r1_bb(b,i) t2_abab(c,a,k,j)
    //              += -0.50 P(i,j) P(a,b) <k,a||c,d>_abab r1_bb(b,i) t2_abab(c,d,k,j)
    //              += -0.50 P(i,j) P(a,b) <k,a||d,c>_abab r1_bb(b,i) t2_abab(d,c,k,j)
    //              += +0.50 P(i,j) P(a,b) <k,a||c,d>_bbbb r1_bb(b,i) t2_bbbb(c,d,j,k)
    //              += -1.00 P(i,j) P(a,b) <k,a||c,j>_abab r2_abab(c,b,k,i)
    //              += -1.00 P(i,j) P(a,b) d-_bb(k,c) r1_1_bb(b,i) t2_bbbb(c,a,j,k)
    //              += -1.00 P(i,j) P(a,b) <l,k||c,j>_bbbb r1_bb(b,l) t2_bbbb(c,a,i,k)
    //              += +1.00 P(i,j) P(a,b) <l,k||c,d>_aaaa r2_abab(d,b,l,i) t2_abab(c,a,k,j)
    //              += +1.00 P(i,j) P(a,b) <l,k||d,c>_abab r2_abab(d,b,l,i) t2_bbbb(c,a,j,k)
    //              += +1.00 P(i,j) P(a,b) <k,l||c,d>_abab r2_bbbb(d,b,i,l) t2_abab(c,a,k,j)
    //              += +1.00 P(i,j) P(a,b) d-_aa(k,c) r1_1_bb(b,i) t2_abab(c,a,k,j)
    //              += +1.00 P(i,j) P(a,b) d-_bb(a,c) r1_bb(b,i) t1_1_bb(c,j)
    //              += -1.00 P(i,j) P(a,b) d-_bb(k,j) r1_bb(b,i) t1_1_bb(a,k)
    //              += +1.00 P(i,j) P(a,b) d-_aa(k,c) r1_bb(b,i) t2_1_abab(c,a,k,j)
    //              += -1.00 P(i,j) P(a,b) d-_bb(k,c) r1_bb(b,i) t2_1_bbbb(c,a,j,k)
    //              += -1.00 P(i,j) P(a,b) f_bb(a,j) r1_bb(b,i)
    //              += +1.00 P(i,j) P(a,b) <k,a||c,d>_abab r1_bb(d,i) t2_abab(c,b,k,j)
    //              += -1.00 P(i,j) P(a,b) d-_bb(a,c) r1_bb(c,i) t1_1_bb(b,j)
    sigmar2_bbbb("R,a,b,i,j") += tmps_["216_bbbb_Lvovo"]("R,a,j,b,i");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["216_bbbb_Lvovo"]("R,a,i,b,j");
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["216_bbbb_Lvovo"]("R,b,j,a,i");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["216_bbbb_Lvovo"]("R,b,i,a,j");
    tmps_["216_bbbb_Lvovo"].~TArrayD();

    // flops: o2v2L1  = o3v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["221_bbbb_Lovvo"]("L,i,a,b,j")  = l2_1["bbbb"]("L,i,k,a,b") * reused_["38_bb_oo"]("j,k");

    // sigmal2_1_bbbb += +2.00 P(i,j) d-_bb(j,c) t1_1_bb(c,k) l2_1_bbbb(i,k,a,b)
    sigmal2_1_bbbb("L,a,b,i,j") += 2.00 * tmps_["221_bbbb_Lovvo"]("L,i,a,b,j");
    sigmal2_1_bbbb("L,a,b,i,j") -= 2.00 * tmps_["221_bbbb_Lovvo"]("L,j,a,b,i");

    // flops: o2v0L1  = o2v1L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["222_bb_Loo"]("L,i,j")  = l1_1["bb"]("L,i,a") * t1_1["bb"]("a,j");

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["220_bb_Loo"]("L,i,j")  = l2["bbbb"]("L,i,k,a,b") * t2["bbbb"]("a,b,k,j");

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["219_bb_Loo"]("L,i,j")  = l2_1["bbbb"]("L,i,k,a,b") * t2_1["bbbb"]("a,b,k,j");

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["218_bb_Loo"]("L,i,j")  = l2["abab"]("L,k,i,a,b") * t2["abab"]("a,b,k,j");

    // flops: o2v0L1  = o3v2L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["217_bb_Loo"]("L,i,j")  = l2_1["abab"]("L,k,i,a,b") * t2_1["abab"]("a,b,k,j");

    // flops: o2v2L1  = o3v2L1 o3v2L1 o3v2L1 o2v3L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["223_bbbb_Lovvo"]("L,i,a,b,j")  = -0.50 * l2["bbbb"]("L,i,m,a,b") * reused_["161_bb_oo"]("j,m");
    tmps_["223_bbbb_Lovvo"]("L,i,a,b,j") -= l2["bbbb"]("L,i,m,a,b") * f["bb_oo"]("j,m");
    tmps_["223_bbbb_Lovvo"]("L,i,a,b,j") -= l2["bbbb"]("L,i,m,a,b") * reused_["30_bb_oo"]("j,m");
    tmps_["223_bbbb_Lovvo"]("L,i,a,b,j") += l1["bb"]("L,i,d") * eri["bbbb_vovv"]("d,j,a,b");
    tmps_["223_bbbb_Lovvo"]("L,i,a,b,j") += l2_1["bbbb"]("L,i,m,a,b") * dp["bb_oo"]("j,m");
    tmps_["223_bbbb_Lovvo"]("L,i,a,b,j") += l2_1["bbbb"]("L,i,m,a,b") * reused_["181_bb_oo"]("j,m");
    tmps_["223_bbbb_Lovvo"]("L,i,a,b,j") += l2["bbbb"]("L,i,m,a,b") * reused_["38_bb_oo"]("j,m");
    tmps_["223_bbbb_Lovvo"]("L,i,a,b,j") += l2["bbbb"]("L,i,m,a,b") * reused_["35_bb_oo"]("j,m");
    tmps_["223_bbbb_Lovvo"]("L,i,a,b,j") -= l2_1["bbbb"]("L,i,m,a,b") * reused_["48_bb_oo"]("j,m");
    tmps_["223_bbbb_Lovvo"]("L,i,a,b,j") -= l2_1["bbbb"]("L,i,m,a,b") * reused_["40_bb_oo"]("j,m");
    tmps_["223_bbbb_Lovvo"]("L,i,a,b,j") += t0_1 * tmps_["221_bbbb_Lovvo"]("L,i,a,b,j");
    tmps_["223_bbbb_Lovvo"]("L,i,a,b,j") += eri["bbbb_oovv"]("j,m,a,b") * tmps_["222_bb_Loo"]("L,i,m");
    tmps_["223_bbbb_Lovvo"]("L,i,a,b,j") += eri["bbbb_oovv"]("j,k,a,b") * tmps_["218_bb_Loo"]("L,i,k");
    tmps_["223_bbbb_Lovvo"]("L,i,a,b,j") -= 0.50 * tmps_["220_bb_Loo"]("L,i,k") * eri["bbbb_oovv"]("j,k,a,b");
    tmps_["223_bbbb_Lovvo"]("L,i,a,b,j") -= 0.50 * tmps_["219_bb_Loo"]("L,i,k") * eri["bbbb_oovv"]("j,k,a,b");
    tmps_["223_bbbb_Lovvo"]("L,i,a,b,j") += tmps_["217_bb_Loo"]("L,i,k") * eri["bbbb_oovv"]("j,k,a,b");

    // sigmal2_bbbb += +1.00 P(i,j) d-_bb(j,c) t0_1 t1_1_bb(c,k) l2_1_bbbb(i,k,a,b)
    //              += +1.00 P(i,j) d+_bb(j,k) l2_1_bbbb(i,k,a,b)
    //              += -1.00 P(i,j) <j,c||a,b>_bbbb l1_bb(i,c)
    //              += -0.50 P(i,j) <l,j||c,d>_abab t2_abab(c,d,l,k) l2_bbbb(i,k,a,b)
    //              += -0.50 P(i,j) <l,j||d,c>_abab t2_abab(d,c,l,k) l2_bbbb(i,k,a,b)
    //              += -1.00 P(i,j) f_bb(j,k) l2_bbbb(i,k,a,b)
    //              += +1.00 P(i,j) <j,l||c,k>_bbbb t1_1_bb(c,l) l2_1_bbbb(i,k,a,b)
    //              += -0.50 P(i,j) <j,l||c,d>_bbbb t2_1_bbbb(c,d,k,l) l2_1_bbbb(i,k,a,b)
    //              += +1.00 P(i,j) d-_bb(j,c) t1_1_bb(c,k) l2_bbbb(i,k,a,b)
    //              += +1.00 P(i,j) d-_bb(j,k) t0_1 l2_bbbb(i,k,a,b)
    //              += -1.00 P(i,j) <l,j||c,k>_abab t1_1_aa(c,l) l2_1_bbbb(i,k,a,b)
    //              += -0.50 P(i,j) <l,j||c,d>_abab t2_1_abab(c,d,l,k) l2_1_bbbb(i,k,a,b)
    //              += -0.50 P(i,j) <l,j||d,c>_abab t2_1_abab(d,c,l,k) l2_1_bbbb(i,k,a,b)
    //              += -0.50 P(i,j) <j,l||c,d>_bbbb t2_bbbb(c,d,k,l) l2_bbbb(i,k,a,b)
    //              += -1.00 P(i,j) f_bb(j,c) t1_1_bb(c,k) l2_1_bbbb(i,k,a,b)
    //              += +1.00 P(i,j) <j,k||a,b>_bbbb t1_1_bb(c,k) l1_1_bb(i,c)
    //              += +0.50 P(i,j) <j,l||a,b>_bbbb t2_abab(d,c,k,l) l2_abab(k,i,d,c)
    //              += +0.50 P(i,j) <j,l||a,b>_bbbb t2_abab(c,d,k,l) l2_abab(k,i,c,d)
    //              += -0.50 P(i,j) <j,l||a,b>_bbbb t2_bbbb(d,c,k,l) l2_bbbb(i,k,d,c)
    //              += -0.50 P(i,j) <j,l||a,b>_bbbb t2_1_bbbb(d,c,k,l) l2_1_bbbb(i,k,d,c)
    //              += +0.50 P(i,j) <j,l||a,b>_bbbb t2_1_abab(d,c,k,l) l2_1_abab(k,i,d,c)
    //              += +0.50 P(i,j) <j,l||a,b>_bbbb t2_1_abab(c,d,k,l) l2_1_abab(k,i,c,d)
    sigmal2_bbbb("L,a,b,i,j") += tmps_["223_bbbb_Lovvo"]("L,i,a,b,j");
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["223_bbbb_Lovvo"]("L,j,a,b,i");
    tmps_["223_bbbb_Lovvo"].~TArrayD();

    // flops: o2v2L1  = o3v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["224_aaaa_Lovvo"]("L,i,a,b,j")  = l2_1["aaaa"]("L,i,k,a,b") * reused_["120_aa_oo"]("j,k");

    // sigmal2_1_aaaa += +2.00 P(i,j) d-_aa(j,c) t1_1_aa(c,k) l2_1_aaaa(i,k,a,b)
    sigmal2_1_aaaa("L,a,b,i,j") += 2.00 * tmps_["224_aaaa_Lovvo"]("L,i,a,b,j");
    sigmal2_1_aaaa("L,a,b,i,j") -= 2.00 * tmps_["224_aaaa_Lovvo"]("L,j,a,b,i");

    // flops: o2v2L1  = o3v2L1 o3v2L1 o3v2L1 o2v0 o3v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v0 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["225_aaaa_Loovv"]("L,i,j,a,b")  = -0.50 * eri["aaaa_oovv"]("i,k,a,b") * tmps_["202_aa_Loo"]("L,j,k");
    tmps_["225_aaaa_Loovv"]("L,i,j,a,b") -= 0.50 * eri["aaaa_oovv"]("i,k,a,b") * tmps_["203_aa_Loo"]("L,j,k");
    tmps_["225_aaaa_Loovv"]("L,i,j,a,b") += eri["aaaa_oovv"]("i,k,a,b") * tmps_["206_aa_Loo"]("L,j,k");
    tmps_["225_aaaa_Loovv"]("L,i,j,a,b") += (reused_["120_aa_oo"]("i,l") + -0.50 * reused_["174_aa_oo"]("i,l")) * l2["aaaa"]("L,j,l,a,b");
    tmps_["225_aaaa_Loovv"]("L,i,j,a,b") += eri["aaaa_vovv"]("e,i,a,b") * l1["aa"]("L,j,e");
    tmps_["225_aaaa_Loovv"]("L,i,j,a,b") -= l2_1["aaaa"]("L,j,l,a,b") * reused_["121_aa_oo"]("l,i");
    tmps_["225_aaaa_Loovv"]("L,i,j,a,b") -= l2["aaaa"]("L,j,l,a,b") * reused_["16_aa_oo"]("i,l");
    tmps_["225_aaaa_Loovv"]("L,i,j,a,b") -= f["aa_oo"]("i,l") * l2["aaaa"]("L,j,l,a,b");
    tmps_["225_aaaa_Loovv"]("L,i,j,a,b") += dp["aa_oo"]("i,l") * l2_1["aaaa"]("L,j,l,a,b");
    tmps_["225_aaaa_Loovv"]("L,i,j,a,b") += l2["aaaa"]("L,j,l,a,b") * reused_["17_aa_oo"]("i,l");
    tmps_["225_aaaa_Loovv"]("L,i,j,a,b") += l2_1["aaaa"]("L,j,l,a,b") * reused_["125_aa_oo"]("i,l");
    tmps_["225_aaaa_Loovv"]("L,i,j,a,b") += l2_1["aaaa"]("L,j,l,a,b") * reused_["182_aa_oo"]("i,l");
    tmps_["225_aaaa_Loovv"]("L,i,j,a,b") += t0_1 * tmps_["224_aaaa_Lovvo"]("L,j,a,b,i");
    tmps_["225_aaaa_Loovv"]("L,i,j,a,b") += eri["aaaa_oovv"]("i,l,a,b") * tmps_["209_aa_Loo"]("L,j,l");
    tmps_["225_aaaa_Loovv"]("L,i,j,a,b") += eri["aaaa_oovv"]("i,k,a,b") * tmps_["207_aa_Loo"]("L,j,k");

    // sigmal2_aaaa += +1.00 P(i,j) d-_aa(j,c) t1_1_aa(c,k) l2_aaaa(i,k,a,b)
    //              += -0.50 P(i,j) <j,l||c,d>_aaaa t2_aaaa(c,d,k,l) l2_aaaa(i,k,a,b)
    //              += -1.00 P(i,j) <j,c||a,b>_aaaa l1_aa(i,c)
    //              += -1.00 P(i,j) f_aa(j,c) t1_1_aa(c,k) l2_1_aaaa(i,k,a,b)
    //              += -0.50 P(i,j) <j,l||c,d>_abab t2_abab(c,d,k,l) l2_aaaa(i,k,a,b)
    //              += -0.50 P(i,j) <j,l||d,c>_abab t2_abab(d,c,k,l) l2_aaaa(i,k,a,b)
    //              += -1.00 P(i,j) f_aa(j,k) l2_aaaa(i,k,a,b)
    //              += +1.00 P(i,j) d+_aa(j,k) l2_1_aaaa(i,k,a,b)
    //              += +1.00 P(i,j) d-_aa(j,k) t0_1 l2_aaaa(i,k,a,b)
    //              += -1.00 P(i,j) <j,l||k,c>_abab t1_1_bb(c,l) l2_1_aaaa(i,k,a,b)
    //              += -0.50 P(i,j) <j,l||c,d>_abab t2_1_abab(c,d,k,l) l2_1_aaaa(i,k,a,b)
    //              += -0.50 P(i,j) <j,l||d,c>_abab t2_1_abab(d,c,k,l) l2_1_aaaa(i,k,a,b)
    //              += +1.00 P(i,j) <j,l||c,k>_aaaa t1_1_aa(c,l) l2_1_aaaa(i,k,a,b)
    //              += -0.50 P(i,j) <j,l||c,d>_aaaa t2_1_aaaa(c,d,k,l) l2_1_aaaa(i,k,a,b)
    //              += +1.00 P(i,j) d-_aa(j,c) t0_1 t1_1_aa(c,k) l2_1_aaaa(i,k,a,b)
    //              += +0.50 P(i,j) <j,l||a,b>_aaaa t2_abab(d,c,l,k) l2_abab(i,k,d,c)
    //              += +0.50 P(i,j) <j,l||a,b>_aaaa t2_abab(c,d,l,k) l2_abab(i,k,c,d)
    //              += -0.50 P(i,j) <j,l||a,b>_aaaa t2_aaaa(d,c,k,l) l2_aaaa(i,k,d,c)
    //              += -0.50 P(i,j) <j,l||a,b>_aaaa t2_1_aaaa(d,c,k,l) l2_1_aaaa(i,k,d,c)
    //              += +1.00 P(i,j) <j,k||a,b>_aaaa t1_1_aa(c,k) l1_1_aa(i,c)
    //              += +0.50 P(i,j) <j,l||a,b>_aaaa t2_1_abab(d,c,l,k) l2_1_abab(i,k,d,c)
    //              += +0.50 P(i,j) <j,l||a,b>_aaaa t2_1_abab(c,d,l,k) l2_1_abab(i,k,c,d)
    sigmal2_aaaa("L,a,b,i,j") += tmps_["225_aaaa_Loovv"]("L,j,i,a,b");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["225_aaaa_Loovv"]("L,i,j,a,b");
    tmps_["225_aaaa_Loovv"].~TArrayD();

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["226_aa_Lov"]("L,i,a")  = 2.00 * l2["aaaa"]("L,i,j,b,a") * reused_["54_aa_vo"]("b,j");
    tmps_["226_aa_Lov"].~TArrayD();

    // flops: o1v1L1  = o3v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["227_aa_Lvo"]("L,a,i")  = -0.50 * l2_1["aaaa"]("L,j,k,b,a") * reused_["134_aaaa_vooo"]("b,j,k,i");
    tmps_["227_aa_Lvo"].~TArrayD();

    // flops: o2v2L1  = o3v3L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["228_baab_Lovov"]("L,i,a,j,b")  = eri["abab_vovo"]("c,i,a,k") * l2_1["abab"]("L,j,k,c,b");

    // sigmal2_1_abab  = -1.00 <c,j||a,k>_abab l2_1_abab(i,k,c,b)
    sigmal2_1_abab("L,a,b,i,j")  = -1.00 * tmps_["228_baab_Lovov"]("L,j,a,i,b");

    // flops: o1v1L1  = o3v2L1 o2v2L1 o3v2L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    tmps_["232_aa_Lvo"]("L,a,i")  = -2.00 * l2["abab"]("L,j,m,a,c") * reused_["83_baba_vooo"]("c,j,m,i");
    tmps_["232_aa_Lvo"]("L,a,i") -= 2.00 * l2["abab"]("L,i,k,a,c") * reused_["33_bb_vo"]("c,k");
    tmps_["232_aa_Lvo"]("L,a,i") += l2["aaaa"]("L,j,l,b,a") * reused_["82_aaaa_vooo"]("b,j,l,i");
    tmps_["232_aa_Lvo"]("L,a,i") += 2.00 * l2["abab"]("L,i,k,a,c") * reused_["32_bb_vo"]("c,k");
    tmps_["232_aa_Lvo"]("L,a,i") += 2.00 * l2["aaaa"]("L,i,j,b,a") * reused_["54_aa_vo"]("b,j");

    // sigmal1_1_aa += -0.50 d-_aa(i,c) t2_aaaa(c,b,j,k) l2_aaaa(j,k,b,a)
    //              += +1.00 d-_bb(k,c) t2_bbbb(c,b,j,k) l2_abab(i,j,a,b)
    //              += +0.50 d-_aa(i,c) t2_abab(c,b,j,k) l2_abab(j,k,a,b)
    //              += +0.50 d-_aa(i,c) t2_abab(c,b,k,j) l2_abab(k,j,a,b)
    //              += -1.00 d-_aa(k,c) t2_abab(c,b,k,j) l2_abab(i,j,a,b)
    //              += -1.00 d-_aa(k,c) t2_aaaa(c,b,j,k) l2_aaaa(i,j,b,a)
    sigmal1_1_aa("L,a,i") -= 0.50 * tmps_["232_aa_Lvo"]("L,a,i");

    // flops: o1v1L1  = o3v2L1 o2v1L1 o3v2L1 o1v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    tmps_["235_aa_Lvo"]("L,a,i")  = -0.50 * l2_1["aaaa"]("L,j,l,c,a") * reused_["134_aaaa_vooo"]("c,j,l,i");
    tmps_["235_aa_Lvo"]("L,a,i") += l1_1["aa"]("L,j,a") * reused_["120_aa_oo"]("i,j");
    tmps_["235_aa_Lvo"]("L,a,i") += l2_1["abab"]("L,j,k,a,b") * reused_["190_baba_vooo"]("b,j,k,i");

    // sigmal1_1_aa += +1.00 d-_aa(i,c) t2_1_abab(c,b,j,k) l2_1_abab(j,k,a,b)
    //              += +1.00 d-_aa(i,c) t2_1_abab(c,b,k,j) l2_1_abab(k,j,a,b)
    //              += +2.00 d-_aa(i,b) t1_1_aa(b,j) l1_1_aa(j,a)
    //              += -1.00 d-_aa(i,c) t2_1_aaaa(c,b,j,k) l2_1_aaaa(j,k,b,a)
    sigmal1_1_aa("L,a,i") += 2.00 * tmps_["235_aa_Lvo"]("L,a,i");

    // flops: o2v2L1  = o3v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["239_abab_Loovv"]("L,i,j,a,b")  = reused_["120_aa_oo"]("i,k") * l2_1["abab"]("L,k,j,a,b");

    // sigmal2_1_abab += +2.00 d-_aa(i,c) t1_1_aa(c,k) l2_1_abab(k,j,a,b)
    sigmal2_1_abab("L,a,b,i,j") += 2.00 * tmps_["239_abab_Loovv"]("L,i,j,a,b");

    // flops: o0v2L1  = o2v3L1
    //  mems: o0v2L1  = o0v2L1
    tmps_["238_bb_Lvv"]("L,a,b")  = l2_1["bbbb"]("L,i,j,c,a") * t2["bbbb"]("b,c,i,j");

    // flops: o1v1L1  = o0v2L1 o1v2L1
    //  mems: o1v1L1  = o0v2L1 o1v1L1
    tmps_["246_bb_Lvo"]("L,a,i")  = (tmps_["199_bb_Lvv"]("L,b,a") + -0.50 * tmps_["238_bb_Lvv"]("L,b,a")) * t1_1["bb"]("b,i");

    // flops: o1v1L1  = o0v2L1 o1v2L1
    //  mems: o1v1L1  = o0v2L1 o1v1L1
    tmps_["245_aa_Lvo"]("L,a,i")  = (tmps_["99_aa_Lvv"]("L,b,a") + -0.50 * tmps_["186_aa_Lvv"]("L,b,a")) * t1_1["aa"]("b,i");

    // flops: o0v0L1  = o2v2L1 o2v2L1 o1v1L1 o1v1L1 o0v0L1 o0v0L1 o0v0L1 o2v2L1 o0v0L1
    //  mems: o0v0L1  = o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1 o0v0L1
    tmps_["244_L"]("L")  = 0.25 * l2["aaaa"]("L,i,k,a,c") * t2_1["aaaa"]("a,c,i,k");
    tmps_["244_L"]("L") += 0.25 * l2["bbbb"]("L,l,j,d,b") * t2_1["bbbb"]("d,b,l,j");
    tmps_["244_L"]("L") += l1["aa"]("L,i,c") * t1_1["aa"]("c,i");
    tmps_["244_L"]("L") += l1["bb"]("L,l,b") * t1_1["bb"]("b,l");
    tmps_["244_L"]("L") += l2["abab"]("L,i,j,a,b") * t2_1["abab"]("a,b,i,j");

    // flops: o2v0L1  = o3v2L1 o3v2L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
    tmps_["243_bb_Loo"]("L,i,j")  = l2["bbbb"]("L,k,i,a,b") * t2["bbbb"]("a,b,k,j");
    tmps_["243_bb_Loo"]("L,i,j") += l2_1["bbbb"]("L,k,i,a,b") * t2_1["bbbb"]("a,b,k,j");

    // flops: o2v0L1  = o3v2L1 o3v2L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
    tmps_["242_aa_Loo"]("L,i,j")  = l2["aaaa"]("L,k,i,a,b") * t2["aaaa"]("a,b,k,j");
    tmps_["242_aa_Loo"]("L,i,j") += l2_1["aaaa"]("L,k,i,a,b") * t2_1["aaaa"]("a,b,k,j");

    // flops: o0v2L1  = o2v3L1 o2v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["241_bb_Lvv"]("L,a,b")  = l2["bbbb"]("L,i,j,a,c") * t2["bbbb"]("b,c,i,j");
    tmps_["241_bb_Lvv"]("L,a,b") += l2_1["bbbb"]("L,i,j,a,c") * t2_1["bbbb"]("b,c,i,j");

    // flops: o0v2L1  = o2v3L1 o2v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["240_aa_Lvv"]("L,a,b")  = l2["aaaa"]("L,i,j,a,c") * t2["aaaa"]("b,c,i,j");
    tmps_["240_aa_Lvv"]("L,a,b") += l2_1["aaaa"]("L,i,j,a,c") * t2_1["aaaa"]("b,c,i,j");

    // flops: o1v1L1  = o2v1L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["237_aa_Lvo"]("L,a,i")  = -1.00 * dp["aa_ov"]("j,a") * tmps_["205_aa_Loo"]("L,i,j");

    // flops: o2v0L1  = o3v1L1 o3v1L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
    tmps_["236_aa_Loo"]("L,i,j")  = t1_1["aa"]("b,l") * tmps_["138_aaaa_Loovo"]("L,i,l,b,j");
    tmps_["236_aa_Loo"]("L,i,j") += t1_1["bb"]("a,k") * tmps_["181_abba_Loovo"]("L,i,k,a,j");

    // flops: o1v1L1  = o2v2L1 o2v2L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1
    tmps_["234_bb_Lvo"]("L,a,i")  = l1["aa"]("L,j,b") * t2["abab"]("b,a,j,i");
    tmps_["234_bb_Lvo"]("L,a,i") += l1_1["aa"]("L,j,b") * t2_1["abab"]("b,a,j,i");

    // flops: o2v0L1  = o3v2L1 o3v2L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
    tmps_["233_aa_Loo"]("L,i,j")  = -2.00 * l2["abab"]("L,i,l,a,c") * t2_1["abab"]("a,c,j,l");
    tmps_["233_aa_Loo"]("L,i,j") += l2["aaaa"]("L,i,k,a,b") * t2_1["aaaa"]("a,b,k,j");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["231_aa_Lov"]("L,i,a")  = l2["aaaa"]("L,j,i,b,a") * t1_1["aa"]("b,j");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["230_aa_Lov"]("L,i,a")  = l2["aaaa"]("L,i,j,a,b") * t1_1["aa"]("b,j");

    // flops: o2v0L1  = o2v1L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["229_aa_Loo"]("L,i,j")  = l1["aa"]("L,i,a") * t1_1["aa"]("a,j");

    // flops: o1v1L1  = o2v1L1 o3v2L1 o3v2L1 o2v3L1 o2v2L1 o2v1L1 o2v1L1 o2v2 o2v2L1 o1v1L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o2v3L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o1v1L1 o1v1L1 o3v2L1 o1v1L1 o2v3L1 o1v1L1 o1v1L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o2v3L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o3v2L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o3v2L1 o1v1L1 o2v1L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o1v1L1 o1v1L1 o3v2L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o1v1L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v2L1 o1v1L1 o2v3L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o2v3L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o1v2L1 o1v1L1 o1v2L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o2v1L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v3L1 o1v1L1 o2v1L1 o1v1L1 o2v3L1 o1v1L1 o2v3L1 o1v1L1 o2v1L1 o1v1L1 o1v2L1 o1v1L1 o2v2L1 o1v1L1 o3v1L1 o1v1L1 o3v2L1 o1v1L1 o2v1L1 o1v1L1 o1v2L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o3v1L1 o1v1L1 o3v1L1 o1v1L1 o2v2L1 o1v1L1 o3v1L1 o1v1L1 o3v2L1 o1v1L1 o1v3L1 o1v1L1 o3v1L1 o1v1L1 o1v3L1 o1v1L1 o2v1L1 o1v1L1 o3v1L1 o1v1L1 o2v1L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o3v1L1 o1v1L1 o2v2L1 o1v1L1 o3v1L1 o1v1L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o1v3L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o3v1L1 o1v1L1 o2v1L1 o1v1L1 o3v1L1 o1v1L1 o2v1L1 o1v1L1 o1v2L1 o1v1L1 o3v2L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o3v1L1 o1v1L1 o3v2L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o3v2L1 o1v1L1 o2v1L1 o1v1L1 o3v1L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o1v3L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o2v2 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    tmps_["247_aa_Lov"]("L,i,a")  = -0.50 * f["aa_ov"]("k,a") * tmps_["202_aa_Loo"]("L,i,k");
    tmps_["247_aa_Lov"]("L,i,a") -= 0.50 * l2["aaaa"]("L,n,k,e,a") * reused_["130_aaaa_vooo"]("e,n,k,i");
    tmps_["247_aa_Lov"]("L,i,a") -= 0.50 * l2_1["aaaa"]("L,n,k,e,a") * reused_["189_aaaa_vooo"]("e,n,k,i");
    tmps_["247_aa_Lov"]("L,i,a") -= eri["abab_vvvo"]("f,d,a,m") * l2["abab"]("L,i,m,f,d");
    tmps_["247_aa_Lov"]("L,i,a") -= 0.50 * l2_1["abab"]("L,i,m,a,d") * reused_["199_bb_ov"]("m,d");
    tmps_["247_aa_Lov"]("L,i,a") -= l1_1["aa"]("L,n,a") * reused_["125_aa_oo"]("i,n");
    tmps_["247_aa_Lov"]("L,i,a") += l1["aa"]("L,n,a") * reused_["16_aa_oo"]("i,n");
    tmps_["247_aa_Lov"]("L,i,a") += (reused_["69_aaaa_voov"]("e,n,i,a") + eri["aaaa_vovo"]("e,i,a,n")) * l1["aa"]("L,n,e");
    tmps_["247_aa_Lov"]("L,i,a") += 0.25 * l2_1["aaaa"]("L,i,n,f,e") * reused_["127_aaaa_vovv"]("a,n,f,e");
    tmps_["247_aa_Lov"]("L,i,a") -= l2["aaaa"]("L,i,n,e,a") * reused_["90_aa_vo"]("e,n");
    tmps_["247_aa_Lov"]("L,i,a") += l2["aaaa"]("L,i,n,f,e") * reused_["165_aaaa_vovv"]("e,n,f,a");
    tmps_["247_aa_Lov"]("L,i,a") += l2["abab"]("L,i,m,f,d") * reused_["167_aabb_vvvo"]("f,a,d,m");
    tmps_["247_aa_Lov"]("L,i,a") += scalars_["4"] * l1["aa"]("L,i,a");
    tmps_["247_aa_Lov"]("L,i,a") += dp["bb_vo"]("d,m") * l2_1["abab"]("L,i,m,a,d");
    tmps_["247_aa_Lov"]("L,i,a") += dp["aa_vv"]("e,a") * l1_1["aa"]("L,i,e");
    tmps_["247_aa_Lov"]("L,i,a") += scalars_["3"] * l1["aa"]("L,i,a");
    tmps_["247_aa_Lov"]("L,i,a") -= l2["abab"]("L,n,o,a,d") * eri["baab_vooo"]("d,i,n,o");
    tmps_["247_aa_Lov"]("L,i,a") -= l2["abab"]("L,i,m,e,g") * reused_["168_abba_vovv"]("e,m,g,a");
    tmps_["247_aa_Lov"]("L,i,a") += dp["aa_ov"]("i,a") * l0_1("L");
    tmps_["247_aa_Lov"]("L,i,a") += 0.50 * l2["aaaa"]("L,n,k,e,a") * eri["aaaa_vooo"]("e,i,n,k");
    tmps_["247_aa_Lov"]("L,i,a") += l1["bb"]("L,m,d") * reused_["1_bbaa_voov"]("d,m,i,a");
    tmps_["247_aa_Lov"]("L,i,a") += 0.25 * l2["aaaa"]("L,i,n,f,e") * reused_["14_aaaa_vovv"]("a,n,f,e");
    tmps_["247_aa_Lov"]("L,i,a") += 0.25 * l2_1["aaaa"]("L,n,k,e,a") * reused_["135_aaaa_vooo"]("e,i,n,k");
    tmps_["247_aa_Lov"]("L,i,a") += l2["abab"]("L,i,m,a,d") * reused_["88_bb_vo"]("d,m");
    tmps_["247_aa_Lov"]("L,i,a") -= l2_1["abab"]("L,i,m,e,g") * reused_["192_abba_vovv"]("e,m,g,a");
    tmps_["247_aa_Lov"]("L,i,a") -= l2_1["abab"]("L,i,m,f,d") * reused_["193_abab_vovv"]("a,m,f,d");
    tmps_["247_aa_Lov"]("L,i,a") -= f["bb_vo"]("d,m") * l2["abab"]("L,i,m,a,d");
    tmps_["247_aa_Lov"]("L,i,a") += 0.50 * eri["aaaa_vvvo"]("e,f,a,n") * l2["aaaa"]("L,i,n,f,e");
    tmps_["247_aa_Lov"]("L,i,a") += f["aa_vo"]("e,n") * l2["aaaa"]("L,i,n,e,a");
    tmps_["247_aa_Lov"]("L,i,a") += scalars_["1"] * l1_1["aa"]("L,i,a");
    tmps_["247_aa_Lov"]("L,i,a") -= l0_1("L") * reused_["9_aa_ov"]("i,a");
    tmps_["247_aa_Lov"]("L,i,a") += 0.50 * l2["aaaa"]("L,n,k,e,a") * reused_["134_aaaa_vooo"]("e,n,k,i");
    tmps_["247_aa_Lov"]("L,i,a") -= l0_1("L") * reused_["147_aa_ov"]("i,a");
    tmps_["247_aa_Lov"]("L,i,a") += l1_1["aa"]("L,n,a") * reused_["121_aa_oo"]("n,i");
    tmps_["247_aa_Lov"]("L,i,a") -= l2["abab"]("L,n,o,a,d") * reused_["175_baab_vooo"]("d,i,n,o");
    tmps_["247_aa_Lov"]("L,i,a") -= l1["aa"]("L,n,a") * reused_["120_aa_oo"]("i,n");
    tmps_["247_aa_Lov"]("L,i,a") += l2_1["abab"]("L,n,o,a,d") * reused_["197_baba_vooo"]("d,n,o,i");
    tmps_["247_aa_Lov"]("L,i,a") -= l2["abab"]("L,k,m,a,d") * reused_["172_bbaa_vooo"]("d,m,i,k");
    tmps_["247_aa_Lov"]("L,i,a") += l2["abab"]("L,k,m,a,d") * reused_["173_bbaa_vooo"]("d,m,i,k");
    tmps_["247_aa_Lov"]("L,i,a") += l2["aaaa"]("L,n,k,e,a") * reused_["170_aaaa_vooo"]("e,n,i,k");
    tmps_["247_aa_Lov"]("L,i,a") -= scalars_["6"] * l1_1["aa"]("L,i,a");
    tmps_["247_aa_Lov"]("L,i,a") -= l2_1["aaaa"]("L,n,k,e,a") * reused_["154_aaaa_vooo"]("e,k,i,n");
    tmps_["247_aa_Lov"]("L,i,a") += scalars_["2"] * l1_1["aa"]("L,i,a");
    tmps_["247_aa_Lov"]("L,i,a") -= dp["aa_oo"]("i,n") * l1_1["aa"]("L,n,a");
    tmps_["247_aa_Lov"]("L,i,a") += l1_1["bb"]("L,m,d") * reused_["105_baab_vovo"]("d,i,a,m");
    tmps_["247_aa_Lov"]("L,i,a") += 0.50 * l2["aaaa"]("L,i,n,e,a") * reused_["87_aa_vo"]("e,n");
    tmps_["247_aa_Lov"]("L,i,a") -= l2_1["abab"]("L,i,m,f,d") * reused_["185_aabb_vvvo"]("f,a,d,m");
    tmps_["247_aa_Lov"]("L,i,a") -= l2["aaaa"]("L,n,k,e,a") * reused_["150_aaaa_oovo"]("i,k,e,n");
    tmps_["247_aa_Lov"]("L,i,a") -= l2_1["abab"]("L,k,m,a,d") * reused_["198_baab_vooo"]("d,k,i,m");
    tmps_["247_aa_Lov"]("L,i,a") += scalars_["5"] * l1["aa"]("L,i,a");
    tmps_["247_aa_Lov"]("L,i,a") -= l2_1["abab"]("L,n,o,a,d") * reused_["188_aabb_oovo"]("i,n,d,o");
    tmps_["247_aa_Lov"]("L,i,a") -= l1["bb"]("L,m,d") * reused_["158_bbaa_voov"]("d,m,i,a");
    tmps_["247_aa_Lov"]("L,i,a") -= f["aa_vv"]("e,a") * l1["aa"]("L,i,e");
    tmps_["247_aa_Lov"]("L,i,a") -= l2["abab"]("L,i,m,f,d") * reused_["169_abab_vovv"]("a,m,f,d");
    tmps_["247_aa_Lov"]("L,i,a") -= l1_1["aa"]("L,n,a") * reused_["182_aa_oo"]("i,n");
    tmps_["247_aa_Lov"]("L,i,a") -= l2_1["aaaa"]("L,i,n,e,a") * reused_["153_aa_vo"]("e,n");
    tmps_["247_aa_Lov"]("L,i,a") -= l1_1["bb"]("L,m,d") * reused_["159_bbaa_voov"]("d,m,i,a");
    tmps_["247_aa_Lov"]("L,i,a") -= l1_1["aa"]("L,i,e") * reused_["123_aa_vv"]("e,a");
    tmps_["247_aa_Lov"]("L,i,a") += l2["aaaa"]("L,i,n,e,a") * reused_["89_aa_ov"]("n,e");
    tmps_["247_aa_Lov"]("L,i,a") += l2["abab"]("L,i,m,a,d") * reused_["91_bb_vo"]("d,m");
    tmps_["247_aa_Lov"]("L,i,a") += l1_1["bb"]("L,m,d") * reused_["103_bbaa_voov"]("d,m,i,a");
    tmps_["247_aa_Lov"]("L,i,a") -= l2_1["aaaa"]("L,i,n,f,e") * reused_["146_aaaa_vvvo"]("f,a,e,n");
    tmps_["247_aa_Lov"]("L,i,a") += l2_1["aaaa"]("L,i,n,f,e") * reused_["187_aaaa_vovv"]("e,n,f,a");
    tmps_["247_aa_Lov"]("L,i,a") += l2_1["abab"]("L,k,m,a,d") * reused_["194_baab_vooo"]("d,i,k,m");
    tmps_["247_aa_Lov"]("L,i,a") -= l2["abab"]("L,n,o,a,d") * reused_["171_abba_oovo"]("i,o,d,n");
    tmps_["247_aa_Lov"]("L,i,a") -= l2["aaaa"]("L,i,n,e,a") * reused_["152_aa_vo"]("e,n");
    tmps_["247_aa_Lov"]("L,i,a") += l1["aa"]("L,i,e") * reused_["59_aa_vv"]("e,a");
    tmps_["247_aa_Lov"]("L,i,a") += l1["aa"]("L,i,e") * reused_["57_aa_vv"]("a,e");
    tmps_["247_aa_Lov"]("L,i,a") += 0.50 * l1["aa"]("L,i,e") * reused_["18_aa_vv"]("a,e");
    tmps_["247_aa_Lov"]("L,i,a") -= l2_1["abab"]("L,k,m,a,d") * reused_["186_bbaa_vooo"]("d,m,i,k");
    tmps_["247_aa_Lov"]("L,i,a") += l2_1["abab"]("L,n,o,a,d") * reused_["195_baba_vooo"]("d,n,o,i");
    tmps_["247_aa_Lov"]("L,i,a") -= l2["abab"]("L,n,o,a,d") * reused_["190_baba_vooo"]("d,n,o,i");
    tmps_["247_aa_Lov"]("L,i,a") += l2_1["aaaa"]("L,n,k,e,a") * reused_["151_aaaa_vooo"]("e,i,k,n");
    tmps_["247_aa_Lov"]("L,i,a") += f["aa_oo"]("i,n") * l1["aa"]("L,n,a");
    tmps_["247_aa_Lov"]("L,i,a") -= 0.50 * l2_1["aaaa"]("L,n,k,e,a") * reused_["136_aaaa_vooo"]("e,n,k,i");
    tmps_["247_aa_Lov"]("L,i,a") += eri["baab_vovo"]("d,i,a,m") * l1["bb"]("L,m,d");
    tmps_["247_aa_Lov"]("L,i,a") -= l1["aa"]("L,n,a") * reused_["17_aa_oo"]("i,n");
    tmps_["247_aa_Lov"]("L,i,a") += l2["abab"]("L,i,m,a,d") * reused_["84_bb_vo"]("d,m");
    tmps_["247_aa_Lov"]("L,i,a") += l1_1["aa"]("L,n,e") * reused_["160_aaaa_vovo"]("e,i,a,n");
    tmps_["247_aa_Lov"]("L,i,a") -= l2["abab"]("L,i,m,f,d") * reused_["166_aabb_vvvo"]("f,a,d,m");
    tmps_["247_aa_Lov"]("L,i,a") += 0.50 * l1["aa"]("L,n,a") * reused_["174_aa_oo"]("i,n");
    tmps_["247_aa_Lov"]("L,i,a") -= l2["aaaa"]("L,i,n,f,e") * reused_["145_aaaa_vvvo"]("f,a,e,n");
    tmps_["247_aa_Lov"]("L,i,a") += l2_1["abab"]("L,i,m,f,d") * reused_["196_aabb_vvvo"]("f,a,d,m");
    tmps_["247_aa_Lov"]("L,i,a") -= dp["aa_oo"]("i,k") * tmps_["231_aa_Lov"]("L,k,a");
    tmps_["247_aa_Lov"]("L,i,a") -= l1_1["aa"]("L,i,e") * reused_["132_aa_vv"]("e,a");
    tmps_["247_aa_Lov"]("L,i,a") += reused_["148_bb_ov"]("o,g") * tmps_["174_bbaa_Lvoov"]("L,g,o,i,a");
    tmps_["247_aa_Lov"]("L,i,a") += 0.50 * eri["aaaa_oovo"]("i,l,a,k") * tmps_["242_aa_Loo"]("L,k,l");
    tmps_["247_aa_Lov"]("L,i,a") -= reused_["191_abba_voov"]("f,m,j,a") * tmps_["169_abab_Loovo"]("L,i,m,f,j");
    tmps_["247_aa_Lov"]("L,i,a") += 0.50 * reused_["13_aa_ov"]("k,a") * tmps_["202_aa_Loo"]("L,i,k");
    tmps_["247_aa_Lov"]("L,i,a") += dp["aa_vv"]("f,a") * tmps_["230_aa_Lov"]("L,i,f");
    tmps_["247_aa_Lov"]("L,i,a") -= t1_1["bb"]("d,m") * tmps_["239_abab_Loovv"]("L,i,m,a,d");
    tmps_["247_aa_Lov"]("L,i,a") += 0.50 * dp["aa_ov"]("k,a") * tmps_["201_aa_Loo"]("L,i,k");
    tmps_["247_aa_Lov"]("L,i,a") += eri["abab_oovo"]("i,j,a,o") * tmps_["218_bb_Loo"]("L,o,j");
    tmps_["247_aa_Lov"]("L,i,a") += eri["abab_oovo"]("i,o,a,m") * tmps_["222_bb_Loo"]("L,m,o");
    tmps_["247_aa_Lov"]("L,i,a") -= l2_1["aaaa"]("L,i,n,e,a") * reused_["155_aa_vo"]("e,n");
    tmps_["247_aa_Lov"]("L,i,a") += eri["aaaa_oovo"]("i,l,a,k") * tmps_["207_aa_Loo"]("L,k,l");
    tmps_["247_aa_Lov"]("L,i,a") -= reused_["143_aaaa_voov"]("f,n,l,a") * tmps_["78_aaaa_Loovo"]("L,i,n,f,l");
    tmps_["247_aa_Lov"]("L,i,a") += eri["aaaa_vovv"]("f,i,a,c") * tmps_["140_aa_Lvv"]("L,f,c");
    tmps_["247_aa_Lov"]("L,i,a") -= 0.50 * reused_["131_aaaa_oovo"]("i,l,a,n") * tmps_["201_aa_Loo"]("L,n,l");
    tmps_["247_aa_Lov"]("L,i,a") += eri["baab_vovv"]("g,i,a,b") * tmps_["114_bb_Lvv"]("L,g,b");
    tmps_["247_aa_Lov"]("L,i,a") -= reused_["13_aa_ov"]("k,a") * tmps_["206_aa_Loo"]("L,i,k");
    tmps_["247_aa_Lov"]("L,i,a") -= 0.50 * reused_["43_abab_oovo"]("i,j,a,m") * tmps_["179_bb_Loo"]("L,m,j");
    tmps_["247_aa_Lov"]("L,i,a") += 0.50 * dp["aa_ov"]("k,a") * tmps_["233_aa_Loo"]("L,i,k");
    tmps_["247_aa_Lov"]("L,i,a") += 0.50 * t0_1 * tmps_["232_aa_Lvo"]("L,a,i");
    tmps_["247_aa_Lov"]("L,i,a") += reused_["147_aa_ov"]("k,a") * tmps_["205_aa_Loo"]("L,i,k");
    tmps_["247_aa_Lov"]("L,i,a") += 0.50 * eri["abab_oovo"]("i,j,a,o") * tmps_["243_bb_Loo"]("L,o,j");
    tmps_["247_aa_Lov"]("L,i,a") += reused_["8_bb_ov"]("o,g") * tmps_["174_bbaa_Lvoov"]("L,g,o,i,a");
    tmps_["247_aa_Lov"]("L,i,a") += eri["aaaa_oovo"]("i,l,a,k") * tmps_["206_aa_Loo"]("L,k,l");
    tmps_["247_aa_Lov"]("L,i,a") += eri["aaaa_vovo"]("f,k,a,n") * tmps_["78_aaaa_Loovo"]("L,i,n,f,k");
    tmps_["247_aa_Lov"]("L,i,a") += l2_1["abab"]("L,i,m,a,d") * reused_["200_bb_vo"]("d,m");
    tmps_["247_aa_Lov"]("L,i,a") += 0.50 * eri["baab_vovv"]("g,i,a,b") * tmps_["241_bb_Lvv"]("L,g,b");
    tmps_["247_aa_Lov"]("L,i,a") -= t0_1 * tmps_["235_aa_Lvo"]("L,a,i");
    tmps_["247_aa_Lov"]("L,i,a") -= dp["aa_ov"]("n,a") * tmps_["229_aa_Loo"]("L,i,n");
    tmps_["247_aa_Lov"]("L,i,a") -= reused_["13_aa_ov"]("n,a") * tmps_["209_aa_Loo"]("L,i,n");
    tmps_["247_aa_Lov"]("L,i,a") += eri["aaaa_oovo"]("i,k,a,n") * tmps_["209_aa_Loo"]("L,n,k");
    tmps_["247_aa_Lov"]("L,i,a") -= 0.50 * reused_["9_aa_ov"]("k,a") * tmps_["201_aa_Loo"]("L,i,k");
    tmps_["247_aa_Lov"]("L,i,a") += eri["abab_oovo"]("i,j,a,o") * tmps_["217_bb_Loo"]("L,o,j");
    tmps_["247_aa_Lov"]("L,i,a") += f["aa_ov"]("n,a") * tmps_["209_aa_Loo"]("L,i,n");
    tmps_["247_aa_Lov"]("L,i,a") += dp["aa_vv"]("f,a") * tmps_["171_aa_Lov"]("L,i,f");
    tmps_["247_aa_Lov"]("L,i,a") -= reused_["183_aaaa_voov"]("f,n,l,a") * tmps_["78_aaaa_Loovo"]("L,i,n,f,l");
    tmps_["247_aa_Lov"]("L,i,a") += 0.50 * reused_["13_aa_ov"]("k,a") * tmps_["203_aa_Loo"]("L,i,k");
    tmps_["247_aa_Lov"]("L,i,a") -= reused_["13_aa_ov"]("k,a") * tmps_["207_aa_Loo"]("L,i,k");
    tmps_["247_aa_Lov"]("L,i,a") += reused_["9_aa_ov"]("k,a") * tmps_["205_aa_Loo"]("L,i,k");
    tmps_["247_aa_Lov"]("L,i,a") += t1_1["bb"]("d,o") * tmps_["228_baab_Lovov"]("L,o,a,i,d");
    tmps_["247_aa_Lov"]("L,i,a") += reused_["43_abab_oovo"]("i,j,a,m") * tmps_["178_bb_Loo"]("L,m,j");
    tmps_["247_aa_Lov"]("L,i,a") -= reused_["1_bbaa_voov"]("g,m,l,a") * tmps_["181_abba_Loovo"]("L,i,m,g,l");
    tmps_["247_aa_Lov"]("L,i,a") += dp["aa_ov"]("i,a") * tmps_["244_L"]("L");
    tmps_["247_aa_Lov"]("L,i,a") += f["aa_ov"]("k,a") * tmps_["206_aa_Loo"]("L,i,k");
    tmps_["247_aa_Lov"]("L,i,a") -= 0.50 * reused_["147_aa_ov"]("k,a") * tmps_["201_aa_Loo"]("L,i,k");
    tmps_["247_aa_Lov"]("L,i,a") -= reused_["184_bbaa_voov"]("g,m,l,a") * tmps_["181_abba_Loovo"]("L,i,m,g,l");
    tmps_["247_aa_Lov"]("L,i,a") += f["aa_ov"]("k,a") * tmps_["207_aa_Loo"]("L,i,k");
    tmps_["247_aa_Lov"]("L,i,a") += reused_["131_aaaa_oovo"]("i,l,a,n") * tmps_["205_aa_Loo"]("L,n,l");
    tmps_["247_aa_Lov"]("L,i,a") -= eri["abab_oovv"]("i,o,a,g") * tmps_["234_bb_Lvo"]("L,g,o");
    tmps_["247_aa_Lov"]("L,i,a") -= eri["baab_vovo"]("g,k,a,m") * tmps_["181_abba_Loovo"]("L,i,m,g,k");
    tmps_["247_aa_Lov"]("L,i,a") += 0.50 * eri["aaaa_vovv"]("f,i,a,c") * tmps_["240_aa_Lvv"]("L,f,c");
    tmps_["247_aa_Lov"]("L,i,a") -= t1_1["aa"]("e,n") * tmps_["224_aaaa_Lovvo"]("L,n,e,a,i");
    tmps_["247_aa_Lov"]("L,i,a") -= 0.50 * f["aa_ov"]("k,a") * tmps_["203_aa_Loo"]("L,i,k");
    tmps_["247_aa_Lov"]("L,i,a") += tmps_["237_aa_Lvo"]("L,a,i");
    tmps_["247_aa_Lov"]("L,i,a") -= dp["aa_ov"]("k,a") * tmps_["236_aa_Loo"]("L,i,k");
    tmps_["247_aa_Lov"]("L,i,a") += eri["aaaa_oovv"]("i,l,a,c") * tmps_["245_aa_Lvo"]("L,c,l");
    tmps_["247_aa_Lov"]("L,i,a") += eri["abab_oovv"]("i,j,a,b") * tmps_["246_bb_Lvo"]("L,b,j");
    tmps_["247_aa_Lov"]("L,i,a") -= dp["aa_oo"]("i,k") * tmps_["171_aa_Lov"]("L,k,a");
    tmps_["237_aa_Lvo"].~TArrayD();
    tmps_["236_aa_Loo"].~TArrayD();
    tmps_["235_aa_Lvo"].~TArrayD();
    tmps_["234_bb_Lvo"].~TArrayD();
    tmps_["233_aa_Loo"].~TArrayD();
    tmps_["232_aa_Lvo"].~TArrayD();
    tmps_["231_aa_Lov"].~TArrayD();
    tmps_["230_aa_Lov"].~TArrayD();
    tmps_["229_aa_Loo"].~TArrayD();
    tmps_["228_baab_Lovov"].~TArrayD();
    tmps_["224_aaaa_Lovvo"].~TArrayD();
    tmps_["174_bbaa_Lvoov"].~TArrayD();
    tmps_["78_aaaa_Loovo"].~TArrayD();

    // sigmal1_aa += +1.00 <l,k||c,d>_bbbb t2_abab(b,c,j,k) t1_1_bb(d,l) l2_1_aaaa(i,j,b,a)
    //            += +1.00 <i,k||c,a>_aaaa t2_aaaa(c,b,j,k) l1_aa(j,b)
    //            += +1.00 <i,b||a,j>_aaaa l1_aa(j,b)
    //            += -0.50 <i,k||b,c>_abab t2_abab(b,c,j,k) l1_aa(j,a)
    //            += -0.50 <i,k||c,b>_abab t2_abab(c,b,j,k) l1_aa(j,a)
    //            += -1.00 <i,k||j,b>_abab t1_1_bb(b,k) l1_1_aa(j,a)
    //            += -0.50 <i,k||b,c>_abab t2_1_abab(b,c,j,k) l1_1_aa(j,a)
    //            += -0.50 <i,k||c,b>_abab t2_1_abab(c,b,j,k) l1_1_aa(j,a)
    //            += +0.25 <l,k||a,j>_aaaa t2_1_aaaa(c,b,l,k) l2_1_aaaa(i,j,c,b)
    //            += +1.00 d-_bb(k,k) t1_1_aa(b,j) l2_aaaa(i,j,b,a)
    //            += +1.00 d-_aa(k,k) t1_1_aa(b,j) l2_aaaa(i,j,b,a)
    //            += -1.00 <k,c||d,a>_aaaa t2_aaaa(d,b,j,k) l2_aaaa(i,j,c,b)
    //            += -1.00 <c,k||a,d>_abab t2_bbbb(d,b,j,k) l2_abab(i,j,c,b)
    //            += -1.00 d-_aa(j,b) t1_1_aa(b,j) l1_aa(i,a)
    //            += -1.00 d+_bb(b,j) l2_1_abab(i,j,a,b)
    //            += -1.00 d+_aa(b,a) l1_1_aa(i,b)
    //            += -1.00 d-_bb(j,b) t1_1_bb(b,j) l1_aa(i,a)
    //            += -0.50 <i,b||j,k>_abab l2_abab(j,k,a,b)
    //            += -0.50 <i,b||k,j>_abab l2_abab(k,j,a,b)
    //            += -1.00 <k,c||a,d>_abab t2_abab(b,d,k,j) l2_abab(i,j,b,c)
    //            += -1.00 d+_aa(i,a) l0_1
    //            += +0.50 <i,b||j,k>_aaaa l2_aaaa(j,k,b,a)
    //            += -1.00 <i,k||a,c>_abab t2_bbbb(c,b,j,k) l1_bb(j,b)
    //            += +0.25 <l,k||a,j>_aaaa t2_aaaa(c,b,l,k) l2_aaaa(i,j,c,b)
    //            += +0.25 <i,b||c,d>_aaaa t2_1_aaaa(c,d,j,k) l2_1_aaaa(j,k,b,a)
    //            += +0.50 f_aa(i,c) t2_1_aaaa(c,b,j,k) l2_1_aaaa(j,k,b,a)
    //            += -0.50 <l,k||c,d>_bbbb t2_bbbb(c,d,j,k) t1_1_bb(b,l) l2_1_abab(i,j,a,b)
    //            += -1.00 <l,k||d,c>_abab t2_bbbb(c,b,j,k) t1_1_aa(d,l) l2_1_abab(i,j,a,b)
    //            += -0.50 <l,k||c,d>_abab t2_abab(c,b,l,k) t1_1_bb(d,j) l2_1_abab(i,j,a,b)
    //            += -0.50 <k,l||c,d>_abab t2_abab(c,b,k,l) t1_1_bb(d,j) l2_1_abab(i,j,a,b)
    //            += +1.00 f_aa(k,c) t2_1_abab(c,b,k,j) l2_1_abab(i,j,a,b)
    //            += -1.00 f_bb(k,j) t1_1_bb(b,k) l2_1_abab(i,j,a,b)
    //            += -0.50 <k,b||c,d>_bbbb t2_1_bbbb(c,d,j,k) l2_1_abab(i,j,a,b)
    //            += +1.00 <k,b||c,j>_abab t1_1_aa(c,k) l2_1_abab(i,j,a,b)
    //            += -1.00 f_bb(k,c) t2_1_bbbb(c,b,j,k) l2_1_abab(i,j,a,b)
    //            += +1.00 f_bb(b,c) t1_1_bb(c,j) l2_1_abab(i,j,a,b)
    //            += +0.50 <k,b||c,d>_abab t2_1_abab(c,d,k,j) l2_1_abab(i,j,a,b)
    //            += +0.50 <k,b||d,c>_abab t2_1_abab(d,c,k,j) l2_1_abab(i,j,a,b)
    //            += +1.00 <k,b||c,j>_bbbb t1_1_bb(c,k) l2_1_abab(i,j,a,b)
    //            += -0.50 <l,k||c,j>_bbbb t2_1_bbbb(c,b,l,k) l2_1_abab(i,j,a,b)
    //            += -1.00 d-_bb(k,c) t1_1_bb(b,j) t1_1_bb(c,k) l2_1_abab(i,j,a,b)
    //            += -1.00 d-_aa(k,c) t1_1_bb(b,j) t1_1_aa(c,k) l2_1_abab(i,j,a,b)
    //            += +1.00 t1_1_bb(b,j) l2_1_abab(i,j,a,b) w0
    //            += -0.50 <l,k||c,j>_abab t2_1_abab(c,b,l,k) l2_1_abab(i,j,a,b)
    //            += -0.50 <k,l||c,j>_abab t2_1_abab(c,b,k,l) l2_1_abab(i,j,a,b)
    //            += -0.50 <k,l||c,d>_abab t2_abab(c,d,k,j) t1_1_bb(b,l) l2_1_abab(i,j,a,b)
    //            += -0.50 <k,l||d,c>_abab t2_abab(d,c,k,j) t1_1_bb(b,l) l2_1_abab(i,j,a,b)
    //            += +2.00 d-_bb(k,c) t1_1_bb(b,k) t1_1_bb(c,j) l2_1_abab(i,j,a,b)
    //            += -1.00 d-_bb(b,c) t0_1 t1_1_bb(c,j) l2_1_abab(i,j,a,b)
    //            += +1.00 d-_bb(k,j) t0_1 t1_1_bb(b,k) l2_1_abab(i,j,a,b)
    //            += -1.00 d-_aa(k,c) t0_1 t2_1_abab(c,b,k,j) l2_1_abab(i,j,a,b)
    //            += +1.00 d-_bb(k,c) t0_1 t2_1_bbbb(c,b,j,k) l2_1_abab(i,j,a,b)
    //            += -0.50 <l,k||c,d>_bbbb t2_bbbb(c,b,l,k) t1_1_bb(d,j) l2_1_abab(i,j,a,b)
    //            += +1.00 <k,l||c,d>_abab t2_abab(c,b,k,j) t1_1_bb(d,l) l2_1_abab(i,j,a,b)
    //            += +0.50 <c,b||a,j>_abab l2_abab(i,j,c,b)
    //            += +0.50 <b,c||a,j>_abab l2_abab(i,j,b,c)
    //            += -1.00 d-_bb(b,c) t1_1_bb(c,j) l2_abab(i,j,a,b)
    //            += +1.00 d-_bb(k,j) t1_1_bb(b,k) l2_abab(i,j,a,b)
    //            += -1.00 d-_aa(k,c) t2_1_abab(c,b,k,j) l2_abab(i,j,a,b)
    //            += +1.00 d-_bb(k,c) t2_1_bbbb(c,b,j,k) l2_abab(i,j,a,b)
    //            += -1.00 <k,c||a,d>_abab t2_1_abab(b,d,k,j) l2_1_abab(i,j,b,c)
    //            += +0.25 <l,k||a,j>_abab t2_1_abab(c,b,l,k) l2_1_abab(i,j,c,b)
    //            += +0.25 <k,l||a,j>_abab t2_1_abab(c,b,k,l) l2_1_abab(i,j,c,b)
    //            += +0.25 <l,k||a,j>_abab t2_1_abab(b,c,l,k) l2_1_abab(i,j,b,c)
    //            += +0.25 <k,l||a,j>_abab t2_1_abab(b,c,k,l) l2_1_abab(i,j,b,c)
    //            += +0.50 <i,l||c,d>_aaaa t2_aaaa(c,b,j,k) t1_1_aa(d,l) l2_1_aaaa(j,k,b,a)
    //            += -1.00 <i,l||c,k>_aaaa t2_1_aaaa(c,b,j,l) l2_1_aaaa(j,k,b,a)
    //            += -0.50 <i,l||j,k>_aaaa t1_1_aa(b,l) l2_1_aaaa(j,k,b,a)
    //            += +1.00 <i,l||c,d>_aaaa t2_aaaa(c,b,k,l) t1_1_aa(d,j) l2_1_aaaa(j,k,b,a)
    //            += -0.25 <i,l||c,d>_aaaa t2_aaaa(c,d,j,k) t1_1_aa(b,l) l2_1_aaaa(j,k,b,a)
    //            += +1.00 f_bb(b,j) l2_abab(i,j,a,b)
    //            += +0.50 <c,b||a,j>_aaaa l2_aaaa(i,j,c,b)
    //            += -1.00 f_aa(b,j) l2_aaaa(i,j,b,a)
    //            += -1.00 d+_aa(j,j) l1_1_aa(i,a)
    //            += +1.00 <i,j||a,b>_abab t1_1_bb(b,j) l0_1
    //            += -0.50 d-_aa(i,c) t2_1_aaaa(c,b,j,k) l2_aaaa(j,k,b,a)
    //            += -1.00 <i,j||b,a>_aaaa t1_1_aa(b,j) l0_1
    //            += -1.00 f_aa(i,b) t1_1_aa(b,j) l1_1_aa(j,a)
    //            += -0.25 <i,b||c,d>_abab t2_abab(c,d,j,k) l2_abab(j,k,a,b)
    //            += -0.25 <i,b||c,d>_abab t2_abab(c,d,k,j) l2_abab(k,j,a,b)
    //            += -0.25 <i,b||d,c>_abab t2_abab(d,c,j,k) l2_abab(j,k,a,b)
    //            += -0.25 <i,b||d,c>_abab t2_abab(d,c,k,j) l2_abab(k,j,a,b)
    //            += -0.50 f_aa(i,c) t2_abab(c,b,j,k) l2_abab(j,k,a,b)
    //            += -0.50 f_aa(i,c) t2_abab(c,b,k,j) l2_abab(k,j,a,b)
    //            += +1.00 d-_aa(i,b) t1_1_aa(b,j) l1_aa(j,a)
    //            += -0.50 <i,l||c,d>_abab t2_abab(c,b,j,k) t1_1_bb(d,l) l2_1_abab(j,k,a,b)
    //            += -0.50 <i,l||c,d>_abab t2_abab(c,b,k,j) t1_1_bb(d,l) l2_1_abab(k,j,a,b)
    //            += +1.00 <i,l||c,k>_abab t2_1_abab(c,b,j,l) l2_1_abab(j,k,a,b)
    //            += +0.50 <i,l||j,k>_abab t1_1_bb(b,l) l2_1_abab(j,k,a,b)
    //            += +0.50 <i,l||k,j>_abab t1_1_bb(b,l) l2_1_abab(k,j,a,b)
    //            += +0.25 <i,l||c,d>_abab t2_abab(c,d,j,k) t1_1_bb(b,l) l2_1_abab(j,k,a,b)
    //            += +0.25 <i,l||c,d>_abab t2_abab(c,d,k,j) t1_1_bb(b,l) l2_1_abab(k,j,a,b)
    //            += +0.25 <i,l||d,c>_abab t2_abab(d,c,j,k) t1_1_bb(b,l) l2_1_abab(j,k,a,b)
    //            += +0.25 <i,l||d,c>_abab t2_abab(d,c,k,j) t1_1_bb(b,l) l2_1_abab(k,j,a,b)
    //            += +1.00 <i,l||d,c>_abab t2_bbbb(c,b,k,l) t1_1_aa(d,j) l2_1_abab(j,k,a,b)
    //            += +1.00 <i,l||c,k>_aaaa t2_abab(c,b,l,j) l2_abab(k,j,a,b)
    //            += +1.00 <i,l||k,c>_abab t2_bbbb(c,b,j,l) l2_abab(k,j,a,b)
    //            += -1.00 <i,l||c,k>_aaaa t2_aaaa(c,b,j,l) l2_aaaa(j,k,b,a)
    //            += +1.00 t0_1 l1_1_aa(i,a) w0
    //            += +0.25 <k,j||b,c>_aaaa t2_1_aaaa(b,c,k,j) l1_1_aa(i,a)
    //            += +0.25 <k,j||b,c>_abab t2_1_abab(b,c,k,j) l1_1_aa(i,a)
    //            += +0.25 <j,k||b,c>_abab t2_1_abab(b,c,j,k) l1_1_aa(i,a)
    //            += +0.25 <k,j||c,b>_abab t2_1_abab(c,b,k,j) l1_1_aa(i,a)
    //            += +0.25 <j,k||c,b>_abab t2_1_abab(c,b,j,k) l1_1_aa(i,a)
    //            += -1.00 d-_aa(j,b) t0_1 t1_1_aa(b,j) l1_1_aa(i,a)
    //            += +1.00 f_bb(j,b) t1_1_bb(b,j) l1_1_aa(i,a)
    //            += -1.00 d-_bb(j,b) t0_1 t1_1_bb(b,j) l1_1_aa(i,a)
    //            += +0.25 <k,j||b,c>_bbbb t2_1_bbbb(b,c,k,j) l1_1_aa(i,a)
    //            += +1.00 f_aa(j,b) t1_1_aa(b,j) l1_1_aa(i,a)
    //            += +1.00 <i,l||d,c>_abab t2_abab(b,c,k,l) t1_1_aa(d,j) l2_1_aaaa(j,k,b,a)
    //            += -1.00 <i,l||k,c>_abab t2_1_abab(b,c,j,l) l2_1_aaaa(j,k,b,a)
    //            += -1.00 d+_bb(j,j) l1_1_aa(i,a)
    //            += +1.00 d+_aa(i,j) l1_1_aa(j,a)
    //            += +1.00 <i,b||a,c>_abab t1_1_bb(c,j) l1_1_bb(j,b)
    //            += +0.50 <k,b||c,d>_aaaa t2_aaaa(c,d,j,k) l2_aaaa(i,j,b,a)
    //            += +0.50 <l,k||j,c>_abab t2_abab(b,c,l,k) l2_aaaa(i,j,b,a)
    //            += +0.50 <k,l||j,c>_abab t2_abab(b,c,k,l) l2_aaaa(i,j,b,a)
    //            += +0.50 <l,k||c,j>_aaaa t2_aaaa(c,b,l,k) l2_aaaa(i,j,b,a)
    //            += +1.00 f_aa(k,c) t2_aaaa(c,b,j,k) l2_aaaa(i,j,b,a)
    //            += -0.50 <b,k||c,d>_abab t2_abab(c,d,j,k) l2_aaaa(i,j,b,a)
    //            += -0.50 <b,k||d,c>_abab t2_abab(d,c,j,k) l2_aaaa(i,j,b,a)
    //            += +1.00 d-_aa(b,j) t0_1 l2_aaaa(i,j,b,a)
    //            += -1.00 f_bb(k,c) t2_abab(b,c,j,k) l2_aaaa(i,j,b,a)
    //            += +1.00 <k,c||d,a>_aaaa t2_1_abab(d,b,k,j) l2_1_abab(i,j,c,b)
    //            += -1.00 <i,l||k,c>_abab t2_abab(b,c,j,l) l2_aaaa(j,k,b,a)
    //            += +1.00 <i,l||c,d>_abab t2_abab(c,b,k,l) t1_1_bb(d,j) l2_1_abab(k,j,a,b)
    //            += +1.00 <i,l||k,c>_abab t2_1_bbbb(c,b,j,l) l2_1_abab(k,j,a,b)
    //            += -1.00 d-_aa(j,j) t0_1 l1_aa(i,a)
    //              += +1.00 f_aa(k,k) r2_aaaa(a,b,i,j)
    //              += -0.50 <l,k||l,k>_abab r2_aaaa(a,b,i,j)
    //              += -0.50 <k,l||k,l>_abab r2_aaaa(a,b,i,j)
    //              += +0.25 <l,k||c,d>_aaaa r2_aaaa(a,b,i,j) t2_aaaa(c,d,l,k)
    //              += -0.50 <l,k||l,k>_bbbb r2_aaaa(a,b,i,j)
    //              += -0.50 <l,k||l,k>_aaaa r2_aaaa(a,b,i,j)
    //              += +0.25 <l,k||c,d>_bbbb r2_aaaa(a,b,i,j) t2_bbbb(c,d,l,k)
    //              += +1.00 f_bb(k,k) r2_aaaa(a,b,i,j)
    //              += -1.00 d-_bb(k,k) r2_aaaa(a,b,i,j) t0_1
    //              += +0.25 <l,k||c,d>_abab r2_aaaa(a,b,i,j) t2_abab(c,d,l,k)
    //              += +0.25 <k,l||c,d>_abab r2_aaaa(a,b,i,j) t2_abab(c,d,k,l)
    //              += +0.25 <l,k||d,c>_abab r2_aaaa(a,b,i,j) t2_abab(d,c,l,k)
    //              += +0.25 <k,l||d,c>_abab r2_aaaa(a,b,i,j) t2_abab(d,c,k,l)
    //            += +1.00 <i,l||c,d>_aaaa t2_abab(c,b,l,k) t1_1_aa(d,j) l2_1_abab(j,k,a,b)
    //            += -0.50 <i,l||c,d>_aaaa t2_abab(c,b,j,k) t1_1_aa(d,l) l2_1_abab(j,k,a,b)
    //            += -0.50 <i,l||c,d>_aaaa t2_abab(c,b,k,j) t1_1_aa(d,l) l2_1_abab(k,j,a,b)
    //            += -1.00 <i,k||c,a>_aaaa t2_abab(c,b,k,j) l1_bb(j,b)
    //            += +0.50 f_aa(i,c) t2_aaaa(c,b,j,k) l2_aaaa(j,k,b,a)
    //            += +0.25 <i,b||c,d>_aaaa t2_aaaa(c,d,j,k) l2_aaaa(j,k,b,a)
    //            += +1.00 f_aa(b,a) l1_aa(i,b)
    //            += +0.25 <l,k||a,j>_abab t2_abab(c,b,l,k) l2_abab(i,j,c,b)
    //            += +0.25 <k,l||a,j>_abab t2_abab(c,b,k,l) l2_abab(i,j,c,b)
    //            += +0.25 <l,k||a,j>_abab t2_abab(b,c,l,k) l2_abab(i,j,b,c)
    //            += +0.25 <k,l||a,j>_abab t2_abab(b,c,k,l) l2_abab(i,j,b,c)
    //            += +1.00 <i,k||b,j>_aaaa t1_1_aa(b,k) l1_1_aa(j,a)
    //            += -0.50 <i,k||b,c>_aaaa t2_1_aaaa(b,c,j,k) l1_1_aa(j,a)
    //            += +0.50 <l,k||d,c>_abab t2_abab(b,c,l,k) t1_1_aa(d,j) l2_1_aaaa(i,j,b,a)
    //            += +0.50 <k,l||d,c>_abab t2_abab(b,c,k,l) t1_1_aa(d,j) l2_1_aaaa(i,j,b,a)
    //            += -1.00 d-_aa(k,j) t0_1 t1_1_aa(b,k) l2_1_aaaa(i,j,b,a)
    //            += -1.00 d-_aa(k,c) t0_1 t2_1_aaaa(c,b,j,k) l2_1_aaaa(i,j,b,a)
    //            += +1.00 d-_aa(b,c) t0_1 t1_1_aa(c,j) l2_1_aaaa(i,j,b,a)
    //            += +1.00 d-_bb(k,c) t0_1 t2_1_abab(b,c,j,k) l2_1_aaaa(i,j,b,a)
    //            += +1.00 <k,l||c,d>_abab t2_aaaa(c,b,j,k) t1_1_bb(d,l) l2_1_aaaa(i,j,b,a)
    //            += -2.00 d-_aa(k,c) t1_1_aa(b,k) t1_1_aa(c,j) l2_1_aaaa(i,j,b,a)
    //            += +0.50 <l,k||c,d>_abab t2_abab(c,d,j,k) t1_1_aa(b,l) l2_1_aaaa(i,j,b,a)
    //            += +0.50 <l,k||d,c>_abab t2_abab(d,c,j,k) t1_1_aa(b,l) l2_1_aaaa(i,j,b,a)
    //            += +0.50 <l,k||c,d>_aaaa t2_aaaa(c,b,l,k) t1_1_aa(d,j) l2_1_aaaa(i,j,b,a)
    //            += -1.00 t1_1_aa(b,j) l2_1_aaaa(i,j,b,a) w0
    //            += +1.00 f_aa(k,c) t2_1_aaaa(c,b,j,k) l2_1_aaaa(i,j,b,a)
    //            += +1.00 d-_bb(k,c) t1_1_aa(b,j) t1_1_bb(c,k) l2_1_aaaa(i,j,b,a)
    //            += +0.50 <l,k||j,c>_abab t2_1_abab(b,c,l,k) l2_1_aaaa(i,j,b,a)
    //            += +0.50 <k,l||j,c>_abab t2_1_abab(b,c,k,l) l2_1_aaaa(i,j,b,a)
    //            += +1.00 d-_aa(k,c) t1_1_aa(b,j) t1_1_aa(c,k) l2_1_aaaa(i,j,b,a)
    //            += -1.00 <k,b||c,j>_aaaa t1_1_aa(c,k) l2_1_aaaa(i,j,b,a)
    //            += -1.00 f_bb(k,c) t2_1_abab(b,c,j,k) l2_1_aaaa(i,j,b,a)
    //            += +1.00 f_aa(k,j) t1_1_aa(b,k) l2_1_aaaa(i,j,b,a)
    //            += +0.50 <k,b||c,d>_aaaa t2_1_aaaa(c,d,j,k) l2_1_aaaa(i,j,b,a)
    //            += -1.00 <b,k||j,c>_abab t1_1_bb(c,k) l2_1_aaaa(i,j,b,a)
    //            += +0.50 <l,k||c,j>_aaaa t2_1_aaaa(c,b,l,k) l2_1_aaaa(i,j,b,a)
    //            += -0.50 <b,k||c,d>_abab t2_1_abab(c,d,j,k) l2_1_aaaa(i,j,b,a)
    //            += -0.50 <b,k||d,c>_abab t2_1_abab(d,c,j,k) l2_1_aaaa(i,j,b,a)
    //            += -1.00 f_aa(b,c) t1_1_aa(c,j) l2_1_aaaa(i,j,b,a)
    //            += +0.50 <l,k||c,d>_aaaa t2_aaaa(c,d,j,k) t1_1_aa(b,l) l2_1_aaaa(i,j,b,a)
    //            += -1.00 <i,k||c,a>_aaaa t2_1_abab(c,b,k,j) l1_1_bb(j,b)
    //            += +1.00 <b,j||a,c>_abab t1_1_bb(c,j) l1_1_aa(i,b)
    //            += -0.50 <k,j||a,c>_abab t2_1_abab(b,c,k,j) l1_1_aa(i,b)
    //            += -0.50 <j,k||a,c>_abab t2_1_abab(b,c,j,k) l1_1_aa(i,b)
    //            += -1.00 d-_aa(k,j) t1_1_aa(b,k) l2_aaaa(i,j,b,a)
    //            += -1.00 d-_aa(k,c) t2_1_aaaa(c,b,j,k) l2_aaaa(i,j,b,a)
    //            += +1.00 d-_aa(b,c) t1_1_aa(c,j) l2_aaaa(i,j,b,a)
    //            += +1.00 d-_bb(k,c) t2_1_abab(b,c,j,k) l2_aaaa(i,j,b,a)
    //            += -1.00 d-_aa(k,k) t1_1_bb(b,j) l2_abab(i,j,a,b)
    //            += -1.00 d-_bb(k,k) t1_1_bb(b,j) l2_abab(i,j,a,b)
    //            += -1.00 <i,k||a,c>_abab t2_1_bbbb(c,b,j,k) l1_1_bb(j,b)
    //            += +1.00 <c,k||a,d>_abab t2_1_abab(b,d,j,k) l2_1_aaaa(i,j,c,b)
    //            += -1.00 <k,c||d,a>_aaaa t2_1_aaaa(d,b,j,k) l2_1_aaaa(i,j,c,b)
    //            += -0.25 <l,k||d,a>_aaaa t2_aaaa(c,b,l,k) t1_1_aa(d,j) l2_1_aaaa(i,j,c,b)
    //            += -0.50 <c,b||d,a>_aaaa t1_1_aa(d,j) l2_1_aaaa(i,j,c,b)
    //            += -1.00 <i,b||k,c>_abab t1_1_bb(c,j) l2_1_abab(k,j,a,b)
    //            += +1.00 <i,l||c,k>_abab t2_abab(c,b,j,l) l2_abab(j,k,a,b)
    //            += +1.00 d-_bb(k,c) t2_abab(b,c,j,k) t0_1 l2_aaaa(i,j,b,a)
    //            += -1.00 d-_aa(b,a) t0_1 l1_aa(i,b)
    //            += -0.50 <k,j||a,c>_abab t2_abab(b,c,k,j) l1_aa(i,b)
    //            += -0.50 <j,k||a,c>_abab t2_abab(b,c,j,k) l1_aa(i,b)
    //            += -0.50 <k,j||c,a>_aaaa t2_aaaa(c,b,k,j) l1_aa(i,b)
    //            += +1.00 <i,l||c,k>_aaaa t2_1_abab(c,b,l,j) l2_1_abab(k,j,a,b)
    //            += -0.50 f_aa(i,c) t2_1_abab(c,b,j,k) l2_1_abab(j,k,a,b)
    //            += -0.50 f_aa(i,c) t2_1_abab(c,b,k,j) l2_1_abab(k,j,a,b)
    //            += -1.00 <i,b||c,k>_abab t1_1_aa(c,j) l2_1_abab(j,k,a,b)
    //            += -0.25 <i,b||c,d>_abab t2_1_abab(c,d,j,k) l2_1_abab(j,k,a,b)
    //            += -0.25 <i,b||c,d>_abab t2_1_abab(c,d,k,j) l2_1_abab(k,j,a,b)
    //            += -0.25 <i,b||d,c>_abab t2_1_abab(d,c,j,k) l2_1_abab(j,k,a,b)
    //            += -0.25 <i,b||d,c>_abab t2_1_abab(d,c,k,j) l2_1_abab(k,j,a,b)
    //            += +0.50 d-_aa(i,c) t2_1_abab(c,b,j,k) l2_abab(j,k,a,b)
    //            += +0.50 d-_aa(i,c) t2_1_abab(c,b,k,j) l2_abab(k,j,a,b)
    //            += +1.00 <i,b||c,k>_aaaa t1_1_aa(c,j) l2_1_aaaa(j,k,b,a)
    //            += -1.00 f_aa(i,j) l1_aa(j,a)
    //            += +0.50 <i,l||c,d>_abab t2_aaaa(c,b,j,k) t1_1_bb(d,l) l2_1_aaaa(j,k,b,a)
    //            += +1.00 <i,b||a,j>_abab l1_bb(j,b)
    //            += +1.00 d-_aa(i,j) t0_1 l1_aa(j,a)
    //            += -1.00 f_bb(k,c) t2_bbbb(c,b,j,k) l2_abab(i,j,a,b)
    //            += -1.00 d-_bb(b,j) t0_1 l2_abab(i,j,a,b)
    //            += -0.50 <l,k||c,j>_abab t2_abab(c,b,l,k) l2_abab(i,j,a,b)
    //            += -0.50 <k,l||c,j>_abab t2_abab(c,b,k,l) l2_abab(i,j,a,b)
    //            += -0.50 <l,k||c,j>_bbbb t2_bbbb(c,b,l,k) l2_abab(i,j,a,b)
    //            += +1.00 f_aa(k,c) t2_abab(c,b,k,j) l2_abab(i,j,a,b)
    //            += +0.50 <k,b||c,d>_abab t2_abab(c,d,k,j) l2_abab(i,j,a,b)
    //            += +0.50 <k,b||d,c>_abab t2_abab(d,c,k,j) l2_abab(i,j,a,b)
    //            += -0.50 <k,b||c,d>_bbbb t2_bbbb(c,d,j,k) l2_abab(i,j,a,b)
    //            += -1.00 <i,b||c,a>_aaaa t1_1_aa(c,j) l1_1_aa(j,b)
    //            += +1.00 <i,k||c,a>_aaaa t2_1_aaaa(c,b,j,k) l1_1_aa(j,b)
    //            += +1.00 <k,c||d,a>_aaaa t2_abab(d,b,k,j) l2_abab(i,j,c,b)
    //            += -0.50 <i,k||b,c>_aaaa t2_aaaa(b,c,j,k) l1_aa(j,a)
    //            += +1.00 <c,k||a,d>_abab t2_abab(b,d,j,k) l2_aaaa(i,j,c,b)
    //            += -1.00 <c,k||a,d>_abab t2_1_bbbb(d,b,j,k) l2_1_abab(i,j,c,b)
    //            += +0.25 <l,k||a,d>_abab t2_abab(c,b,l,k) t1_1_bb(d,j) l2_1_abab(i,j,c,b)
    //            += +0.25 <k,l||a,d>_abab t2_abab(c,b,k,l) t1_1_bb(d,j) l2_1_abab(i,j,c,b)
    //            += +0.25 <l,k||a,d>_abab t2_abab(b,c,l,k) t1_1_bb(d,j) l2_1_abab(i,j,b,c)
    //            += +0.25 <k,l||a,d>_abab t2_abab(b,c,k,l) t1_1_bb(d,j) l2_1_abab(i,j,b,c)
    //            += +0.50 <c,b||a,d>_abab t1_1_bb(d,j) l2_1_abab(i,j,c,b)
    //            += +0.50 <b,c||a,d>_abab t1_1_bb(d,j) l2_1_abab(i,j,b,c)
    //            += +1.00 d-_aa(i,k) t1_1_aa(b,j) l2_aaaa(j,k,b,a)
    //            += +1.00 <j,b||c,a>_aaaa t1_1_aa(c,j) l1_1_aa(i,b)
    //            += -0.50 <k,j||c,a>_aaaa t2_1_aaaa(c,b,k,j) l1_1_aa(i,b)
    //            += -0.50 <i,l||a,k>_aaaa t2_aaaa(c,b,j,l) l2_aaaa(j,k,c,b)
    //            += -0.50 <i,l||a,k>_aaaa t2_1_aaaa(c,b,j,l) l2_1_aaaa(j,k,c,b)
    //            += +1.00 <k,l||a,d>_abab t2_abab(c,d,k,j) t1_1_bb(b,l) l2_1_abab(i,j,c,b)
    //            += -0.50 d-_aa(k,a) t0_1 t2_1_aaaa(c,b,j,k) l2_1_aaaa(i,j,c,b)
    //            += -1.00 d-_aa(c,a) t1_1_aa(b,j) l2_aaaa(i,j,c,b)
    //            += +1.00 d-_aa(i,c) t1_1_bb(b,j) t1_1_aa(c,k) l2_1_abab(k,j,a,b)
    //            += -0.50 d+_aa(k,a) t2_aaaa(c,b,j,k) l2_1_aaaa(i,j,c,b)
    //            += -0.50 <i,l||a,k>_abab t2_abab(c,b,j,l) l2_abab(j,k,c,b)
    //            += -0.50 <i,l||a,k>_abab t2_abab(b,c,j,l) l2_abab(j,k,b,c)
    //            += -1.00 <i,k||a,j>_abab t1_1_bb(b,k) l1_1_bb(j,b)
    //            += -1.00 <l,k||c,d>_aaaa t2_aaaa(c,b,j,k) t1_1_aa(d,l) l2_1_aaaa(i,j,b,a)
    //            += -0.50 <i,l||a,k>_aaaa t2_1_abab(c,b,l,j) l2_1_abab(k,j,c,b)
    //            += -0.50 <i,l||a,k>_aaaa t2_1_abab(b,c,l,j) l2_1_abab(k,j,b,c)
    //            += +1.00 <l,k||a,d>_abab t2_abab(c,d,j,k) t1_1_aa(b,l) l2_1_aaaa(i,j,c,b)
    //            += -0.50 <i,c||d,a>_aaaa t2_abab(d,b,j,k) l2_abab(j,k,c,b)
    //            += -0.50 <i,c||d,a>_aaaa t2_abab(d,b,k,j) l2_abab(k,j,c,b)
    //            += -0.50 <i,c||d,a>_aaaa t2_1_abab(d,b,j,k) l2_1_abab(j,k,c,b)
    //            += -0.50 <i,c||d,a>_aaaa t2_1_abab(d,b,k,j) l2_1_abab(k,j,c,b)
    //            += -0.50 <i,l||d,a>_aaaa t2_aaaa(c,b,k,l) t1_1_aa(d,j) l2_1_aaaa(j,k,c,b)
    //            += +0.50 <i,c||a,d>_abab t2_1_abab(b,d,j,k) l2_1_abab(j,k,b,c)
    //            += +0.50 <i,c||a,d>_abab t2_1_abab(b,d,k,j) l2_1_abab(k,j,b,c)
    //            += +0.50 <i,c||a,d>_abab t2_abab(b,d,j,k) l2_abab(j,k,b,c)
    //            += +0.50 <i,c||a,d>_abab t2_abab(b,d,k,j) l2_abab(k,j,b,c)
    //            += +0.50 d-_aa(k,a) t2_abab(c,b,k,j) t0_1 l2_abab(i,j,c,b)
    //            += +0.50 d-_aa(k,a) t2_abab(b,c,k,j) t0_1 l2_abab(i,j,b,c)
    //            += +0.50 <i,l||a,d>_abab t2_bbbb(c,b,k,l) t1_1_bb(d,j) l2_1_bbbb(j,k,c,b)
    //            += -0.50 d-_aa(k,a) t2_1_aaaa(c,b,j,k) l2_aaaa(i,j,c,b)
    //            += +0.50 d-_aa(k,a) t2_1_abab(c,b,k,j) l2_abab(i,j,c,b)
    //            += +0.50 d-_aa(k,a) t2_1_abab(b,c,k,j) l2_abab(i,j,b,c)
    //            += -0.50 d-_aa(i,c) t2_aaaa(c,b,j,k) t0_1 l2_aaaa(j,k,b,a)
    //            += +1.00 d-_bb(k,c) t2_bbbb(c,b,j,k) t0_1 l2_abab(i,j,a,b)
    //            += +0.50 d-_aa(i,c) t2_abab(c,b,j,k) t0_1 l2_abab(j,k,a,b)
    //            += +0.50 d-_aa(i,c) t2_abab(c,b,k,j) t0_1 l2_abab(k,j,a,b)
    //            += -1.00 d-_aa(k,c) t2_abab(c,b,k,j) t0_1 l2_abab(i,j,a,b)
    //            += -1.00 d-_aa(k,c) t2_aaaa(c,b,j,k) t0_1 l2_aaaa(i,j,b,a)
    //            += -0.50 <l,k||d,a>_aaaa t2_abab(c,b,k,j) t1_1_aa(d,l) l2_1_abab(i,j,c,b)
    //            += -0.50 <l,k||d,a>_aaaa t2_abab(b,c,k,j) t1_1_aa(d,l) l2_1_abab(i,j,b,c)
    //            += -0.50 <i,l||a,k>_abab t2_1_bbbb(c,b,j,l) l2_1_bbbb(j,k,c,b)
    //            += -0.50 <i,l||a,k>_abab t2_bbbb(c,b,j,l) l2_bbbb(j,k,c,b)
    //            += -1.00 <l,k||d,c>_abab t2_abab(b,c,j,k) t1_1_aa(d,l) l2_1_aaaa(i,j,b,a)
    //            += -0.50 <i,l||a,k>_aaaa t2_abab(c,b,l,j) l2_abab(k,j,c,b)
    //            += -0.50 <i,l||a,k>_aaaa t2_abab(b,c,l,j) l2_abab(k,j,b,c)
    //            += +0.50 f_aa(k,a) t2_1_aaaa(c,b,j,k) l2_1_aaaa(i,j,c,b)
    //            += +1.00 <k,c||a,j>_aaaa t1_1_aa(b,k) l2_1_aaaa(i,j,c,b)
    //            += +1.00 <l,k||c,d>_bbbb t2_bbbb(c,b,j,k) t1_1_bb(d,l) l2_1_abab(i,j,a,b)
    //            += -1.00 <l,k||c,d>_aaaa t2_abab(c,b,k,j) t1_1_aa(d,l) l2_1_abab(i,j,a,b)
    //            += +0.50 <i,c||a,d>_abab t2_1_bbbb(d,b,j,k) l2_1_bbbb(j,k,c,b)
    //            += +0.50 <i,c||a,d>_abab t2_bbbb(d,b,j,k) l2_bbbb(j,k,c,b)
    //            += +0.50 d-_aa(i,c) t0_1 t2_1_abab(c,b,j,k) l2_1_abab(j,k,a,b)
    //            += +0.50 d-_aa(i,c) t0_1 t2_1_abab(c,b,k,j) l2_1_abab(k,j,a,b)
    //            += +1.00 d-_aa(i,b) t0_1 t1_1_aa(b,j) l1_1_aa(j,a)
    //            += -0.50 d-_aa(i,c) t0_1 t2_1_aaaa(c,b,j,k) l2_1_aaaa(j,k,b,a)
    //            += +1.00 d-_aa(j,a) t1_1_aa(b,j) l1_aa(i,b)
    //            += +1.00 d-_aa(j,a) t0_1 t1_1_aa(b,j) l1_1_aa(i,b)
    //            += -1.00 <i,k||a,j>_aaaa t1_1_aa(b,k) l1_1_aa(j,b)
    //            += +0.50 <k,l||a,d>_abab t2_aaaa(c,b,j,k) t1_1_bb(d,l) l2_1_aaaa(i,j,c,b)
    //            += -0.50 <i,l||a,k>_abab t2_1_abab(c,b,j,l) l2_1_abab(j,k,c,b)
    //            += -0.50 <i,l||a,k>_abab t2_1_abab(b,c,j,l) l2_1_abab(j,k,b,c)
    //            += -1.00 f_aa(j,a) t1_1_aa(b,j) l1_1_aa(i,b)
    //            += -1.00 d-_aa(c,a) t1_1_bb(b,j) l2_abab(i,j,c,b)
    //            += +1.00 <l,k||d,a>_aaaa t2_aaaa(d,c,j,k) t1_1_aa(b,l) l2_1_aaaa(i,j,c,b)
    //            += -0.50 d-_aa(k,a) t2_aaaa(c,b,j,k) t0_1 l2_aaaa(i,j,c,b)
    //            += +0.50 d-_aa(k,a) t0_1 t2_1_abab(c,b,k,j) l2_1_abab(i,j,c,b)
    //            += +0.50 d-_aa(k,a) t0_1 t2_1_abab(b,c,k,j) l2_1_abab(i,j,b,c)
    //            += -0.50 <k,l||a,d>_abab t2_abab(c,b,k,j) t1_1_bb(d,l) l2_1_abab(i,j,c,b)
    //            += -0.50 <k,l||a,d>_abab t2_abab(b,c,k,j) t1_1_bb(d,l) l2_1_abab(i,j,b,c)
    //            += -1.00 <c,k||a,j>_abab t1_1_bb(b,k) l2_1_abab(i,j,c,b)
    //            += -0.50 <i,l||a,d>_abab t2_abab(c,b,k,l) t1_1_bb(d,j) l2_1_abab(k,j,c,b)
    //            += -0.50 <i,l||a,d>_abab t2_abab(b,c,k,l) t1_1_bb(d,j) l2_1_abab(k,j,b,c)
    //            += +1.00 <l,k||a,d>_abab t2_bbbb(d,c,j,k) t1_1_aa(b,l) l2_1_abab(i,j,b,c)
    //            += -1.00 d-_aa(i,a) t1_1_bb(b,j) l1_bb(j,b)
    //            += -1.00 d-_aa(i,a) t1_1_aa(b,j) l1_aa(j,b)
    //            += -0.25 d-_aa(i,a) t2_1_bbbb(c,b,j,k) l2_bbbb(j,k,c,b)
    //            += -0.25 d-_aa(i,a) t2_1_aaaa(c,b,j,k) l2_aaaa(j,k,c,b)
    //            += -0.25 d-_aa(i,a) t2_1_abab(c,b,j,k) l2_abab(j,k,c,b)
    //            += -0.25 d-_aa(i,a) t2_1_abab(c,b,k,j) l2_abab(k,j,c,b)
    //            += -0.25 d-_aa(i,a) t2_1_abab(b,c,j,k) l2_abab(j,k,b,c)
    //            += -0.25 d-_aa(i,a) t2_1_abab(b,c,k,j) l2_abab(k,j,b,c)
    //            += -0.50 f_aa(k,a) t2_abab(c,b,k,j) l2_abab(i,j,c,b)
    //            += -0.50 f_aa(k,a) t2_abab(b,c,k,j) l2_abab(i,j,b,c)
    //            += +0.50 <l,k||d,a>_aaaa t2_aaaa(c,b,j,k) t1_1_aa(d,l) l2_1_aaaa(i,j,c,b)
    //            += +1.00 <l,k||d,a>_aaaa t2_abab(d,c,k,j) t1_1_aa(b,l) l2_1_abab(i,j,b,c)
    //            += -0.50 f_aa(k,a) t2_1_abab(c,b,k,j) l2_1_abab(i,j,c,b)
    //            += -0.50 f_aa(k,a) t2_1_abab(b,c,k,j) l2_1_abab(i,j,b,c)
    //            += +0.50 <i,l||d,a>_aaaa t2_abab(c,b,l,k) t1_1_aa(d,j) l2_1_abab(j,k,c,b)
    //            += +0.50 <i,l||d,a>_aaaa t2_abab(b,c,l,k) t1_1_aa(d,j) l2_1_abab(j,k,b,c)
    //            += +1.00 <i,k||a,c>_abab t2_1_abab(b,c,j,k) l1_1_aa(j,b)
    //            += +1.00 <i,k||a,c>_abab t2_abab(b,c,j,k) l1_aa(j,b)
    //            += -1.00 <k,c||a,j>_abab t1_1_aa(b,k) l2_1_abab(i,j,b,c)
    //            += -0.50 <i,c||d,a>_aaaa t2_aaaa(d,b,j,k) l2_aaaa(j,k,c,b)
    //            += -0.50 <i,c||d,a>_aaaa t2_1_aaaa(d,b,j,k) l2_1_aaaa(j,k,c,b)
    //            += +1.00 d-_aa(i,c) t1_1_aa(b,j) t1_1_aa(c,k) l2_1_aaaa(j,k,b,a)
    //            += +0.50 f_aa(k,a) t2_aaaa(c,b,j,k) l2_aaaa(i,j,c,b)
    //            += +0.50 d+_aa(k,a) t2_abab(c,b,k,j) l2_1_abab(i,j,c,b)
    //            += +0.50 d+_aa(k,a) t2_abab(b,c,k,j) l2_1_abab(i,j,b,c)
    //            += +1.00 d-_aa(k,a) t1_1_bb(b,j) t1_1_aa(c,k) l2_1_abab(i,j,c,b)
    //            += +1.00 d-_aa(k,a) t1_1_aa(b,j) t1_1_aa(c,k) l2_1_aaaa(i,j,c,b)
    //            += +0.50 <i,l||d,a>_aaaa t2_abab(d,c,j,k) t1_1_aa(b,l) l2_1_abab(j,k,b,c)
    //            += +0.50 <i,l||d,a>_aaaa t2_abab(d,c,k,j) t1_1_aa(b,l) l2_1_abab(k,j,b,c)
    //            += -0.50 <i,l||d,a>_aaaa t2_aaaa(d,c,j,k) t1_1_aa(b,l) l2_1_aaaa(j,k,c,b)
    //            += -0.50 <i,l||a,d>_abab t2_abab(c,d,j,k) t1_1_bb(b,l) l2_1_abab(j,k,c,b)
    //            += -0.50 <i,l||a,d>_abab t2_abab(c,d,k,j) t1_1_bb(b,l) l2_1_abab(k,j,c,b)
    //            += +0.50 <i,l||a,d>_abab t2_bbbb(d,c,j,k) t1_1_bb(b,l) l2_1_bbbb(j,k,c,b)
    //            += +1.00 d-_aa(i,k) t1_1_bb(b,j) l2_abab(k,j,a,b)
    sigmal1_aa("L,a,i") -= tmps_["247_aa_Lov"]("L,i,a");
    tmps_["247_aa_Lov"].~TArrayD();

    // flops: o1v1L1  = o2v1L1 o3v2L1 o3v2L1 o1v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    tmps_["253_bb_Lvo"]("L,a,i")  = -2.00 * l1_1["bb"]("L,j,a") * reused_["38_bb_oo"]("i,j");
    tmps_["253_bb_Lvo"]("L,a,i") -= 2.00 * l2_1["abab"]("L,l,k,d,a") * reused_["208_aabb_vooo"]("d,l,k,i");
    tmps_["253_bb_Lvo"]("L,a,i") += l2_1["bbbb"]("L,j,k,b,a") * reused_["138_bbbb_vooo"]("b,j,k,i");

    // sigmal1_1_bb  = -1.00 d-_bb(i,c) t2_1_bbbb(c,b,j,k) l2_1_bbbb(j,k,b,a)
    //              += +1.00 d-_bb(i,c) t2_1_abab(b,c,j,k) l2_1_abab(j,k,b,a)
    //              += +1.00 d-_bb(i,c) t2_1_abab(b,c,k,j) l2_1_abab(k,j,b,a)
    //              += +2.00 d-_bb(i,b) t1_1_bb(b,j) l1_1_bb(j,a)
    sigmal1_1_bb("L,a,i")  = -1.00 * tmps_["253_bb_Lvo"]("L,a,i");

    // flops: o1v1L1  = o3v2L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    tmps_["254_bb_Lov"]("L,i,a")  = -1.00 * l2["abab"]("L,l,k,d,a") * reused_["77_aabb_vooo"]("d,l,k,i");
    tmps_["254_bb_Lov"]("L,i,a") += l2["bbbb"]("L,i,j,b,a") * reused_["33_bb_vo"]("b,j");
    tmps_["254_bb_Lov"]("L,i,a") -= l2["abab"]("L,l,i,d,a") * reused_["54_aa_vo"]("d,l");
    tmps_["254_bb_Lov"]("L,i,a") -= l2["bbbb"]("L,i,j,b,a") * reused_["32_bb_vo"]("b,j");
    tmps_["254_bb_Lov"]("L,i,a") += 0.50 * l2["bbbb"]("L,j,k,b,a") * reused_["51_bbbb_vooo"]("b,j,k,i");

    // sigmal1_1_bb += -1.00 d-_bb(k,c) t2_bbbb(c,b,j,k) l2_bbbb(i,j,b,a)
    //              += +0.50 d-_bb(i,c) t2_abab(b,c,j,k) l2_abab(j,k,b,a)
    //              += +0.50 d-_bb(i,c) t2_abab(b,c,k,j) l2_abab(k,j,b,a)
    //              += +1.00 d-_aa(k,c) t2_aaaa(c,b,j,k) l2_abab(j,i,b,a)
    //              += +1.00 d-_aa(k,c) t2_abab(c,b,k,j) l2_bbbb(i,j,b,a)
    //              += -0.50 d-_bb(i,c) t2_bbbb(c,b,j,k) l2_bbbb(j,k,b,a)
    sigmal1_1_bb("L,a,i") -= tmps_["254_bb_Lov"]("L,i,a");

    // flops: o2v2L1  = o3v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["260_aabb_Lovvo"]("L,i,a,b,j")  = l2_1["abab"]("L,i,k,a,b") * reused_["38_bb_oo"]("j,k");

    // sigmal2_1_abab += +2.00 d-_bb(j,c) t1_1_bb(c,k) l2_1_abab(i,k,a,b)
    sigmal2_1_abab("L,a,b,i,j") += 2.00 * tmps_["260_aabb_Lovvo"]("L,i,a,b,j");

    // flops: o2v2L1  = o3v3L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["261_bbbb_Lvoov"]("L,a,i,j,b")  = t2["abab"]("c,a,k,i") * l2_1["abab"]("L,k,j,c,b");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["259_bb_Lov"]("L,i,a")  = -1.00 * eri["abab_oovv"]("j,i,b,a") * tmps_["245_aa_Lvo"]("L,b,j");
    tmps_["245_aa_Lvo"].~TArrayD();

    // flops: o1v1L1  = o3v1L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["258_bb_Lov"]("L,i,a")  = -1.00 * reused_["116_abba_oovo"]("j,i,a,k") * tmps_["201_aa_Loo"]("L,k,j");

    // flops: o2v0L1  = o3v1L1 o3v1L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
    tmps_["257_bb_Loo"]("L,i,j")  = t1_1["aa"]("b,l") * tmps_["169_abab_Loovo"]("L,l,i,b,j");
    tmps_["257_bb_Loo"]("L,i,j") += t1_1["bb"]("a,k") * tmps_["113_bbbb_Loovo"]("L,i,k,a,j");

    // flops: o1v1L1  = o1v1L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["256_bb_Lov"]("L,i,a")  = -1.00 * l0_1("L") * reused_["8_bb_ov"]("i,a");

    // flops: o2v0L1  = o3v2L1 o3v2L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
    tmps_["255_bb_Loo"]("L,i,j")  = -0.50 * l2["bbbb"]("L,i,l,c,b") * t2_1["bbbb"]("c,b,l,j");
    tmps_["255_bb_Loo"]("L,i,j") += l2["abab"]("L,k,i,a,b") * t2_1["abab"]("a,b,k,j");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["252_bb_Lov"]("L,i,a")  = l2["bbbb"]("L,i,j,a,b") * t1_1["bb"]("b,j");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["251_bb_Lov"]("L,i,a")  = l2["bbbb"]("L,j,i,b,a") * t1_1["bb"]("b,j");

    // flops: o1v1L1  = o2v1L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["248_bb_Lvo"]("L,a,i")  = -2.00 * l1_1["bb"]("L,j,a") * reused_["38_bb_oo"]("i,j");
    tmps_["248_bb_Lvo"].~TArrayD();

    // flops: o1v1L1  = o3v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["249_bb_Lvo"]("L,a,i")  = 0.50 * l2["bbbb"]("L,j,k,b,a") * reused_["51_bbbb_vooo"]("b,j,k,i");
    tmps_["249_bb_Lvo"].~TArrayD();

    // flops: o2v0L1  = o2v1L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["250_bb_Loo"]("L,i,j")  = l1["bb"]("L,i,a") * t1_1["bb"]("a,j");
}
