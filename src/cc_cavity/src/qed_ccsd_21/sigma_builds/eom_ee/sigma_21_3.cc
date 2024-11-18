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

void hilbert::EOM_EE_QED_CCSD::sigma_ee_21_3() {

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


    // sigmar1_1_bb  = -1.00 d-_bb(j,b) r1_1_bb(b,j) t1_1_bb(a,i)
    //              += -1.00 d-_aa(j,b) r1_1_aa(b,j) t1_1_bb(a,i)
    //              += +2.00 d-_bb(j,b) r1_1_bb(b,i) t1_1_bb(a,j)
    //              += -1.00 <k,j||b,i>_abab r1_aa(b,k) t1_1_bb(a,j)
    //              += -1.00 <k,j||b,i>_bbbb r1_bb(b,k) t1_1_bb(a,j)
    //              += -1.00 f_bb(j,b) r1_bb(b,i) t1_1_bb(a,j)
    //              += +0.50 <k,j||b,c>_bbbb r2_bbbb(b,c,i,k) t1_1_bb(a,j)
    //              += -0.50 <k,j||b,c>_abab r2_abab(b,c,k,i) t1_1_bb(a,j)
    //              += -0.50 <k,j||c,b>_abab r2_abab(c,b,k,i) t1_1_bb(a,j)
    //              += +1.00 d-_bb(j,b) r1_bb(b,i) t0_1 t1_1_bb(a,j)
    //              += +1.00 d-_bb(j,b) r1_bb(a,j) t0_1 t1_1_bb(b,i)
    //              += -1.00 d-_bb(j,j) r0_1 t1_1_bb(a,i)
    //              += -1.00 d-_aa(j,j) r0_1 t1_1_bb(a,i)
    //              += -0.50 <k,j||b,c>_bbbb r1_1_bb(a,k) t2_bbbb(b,c,i,j)
    //              += -0.50 <k,j||b,i>_bbbb r2_1_bbbb(b,a,k,j)
    //              += +1.00 <k,j||b,c>_bbbb r1_1_bb(c,k) t2_bbbb(b,a,i,j)
    //              += -1.00 f_bb(j,i) r1_1_bb(a,j)
    //              += -2.00 d-_aa(j,b) r1_1_bb(a,i) t1_1_aa(b,j)
    //              += +1.00 <j,a||b,c>_abab r1_bb(c,i) t1_1_aa(b,j)
    //              += -0.50 <k,j||b,c>_abab r1_bb(c,i) t2_1_abab(b,a,k,j)
    //              += -0.50 <j,k||b,c>_abab r1_bb(c,i) t2_1_abab(b,a,j,k)
    //              += -1.00 <k,j||b,c>_aaaa r1_1_aa(c,k) t2_abab(b,a,j,i)
    //              += +1.00 d-_bb(j,b) r0_1 t2_bbbb(b,a,i,j) t0_1
    //              += +1.00 r1_1_bb(a,i) w0
    //              += -0.50 <k,j||b,c>_abab r1_1_bb(c,i) t2_abab(b,a,k,j)
    //              += -0.50 <j,k||b,c>_abab r1_1_bb(c,i) t2_abab(b,a,j,k)
    //              += -1.00 <k,j||c,b>_abab r1_1_aa(c,k) t2_bbbb(b,a,i,j)
    //              += +1.00 <j,k||b,c>_abab r1_bb(c,k) t2_1_abab(b,a,j,i)
    //              += -1.00 d-_aa(j,b) r0_1 t2_abab(b,a,j,i) t0_1
    //              += -2.00 d-_bb(a,b) r0_1 t1_1_bb(b,i)
    //              += +2.00 d-_bb(j,i) r0_1 t1_1_bb(a,j)
    //              += -2.00 d-_aa(j,b) r0_1 t2_1_abab(b,a,j,i)
    //              += +2.00 d-_bb(j,b) r0_1 t2_1_bbbb(b,a,i,j)
    //              += +2.00 d-_bb(j,b) r1_1_bb(a,j) t1_1_bb(b,i)
    //              += +1.00 <j,a||b,i>_bbbb r1_1_bb(b,j)
    //              += -0.50 <k,j||b,c>_bbbb r1_bb(a,k) t2_1_bbbb(b,c,i,j)
    //              += +1.00 <k,j||b,i>_bbbb r1_bb(a,k) t1_1_bb(b,j)
    //              += -1.00 d+_bb(j,j) r1_bb(a,i)
    //              += +0.50 <k,j||b,c>_bbbb r2_bbbb(c,a,k,j) t1_1_bb(b,i)
    //              += -1.00 d-_aa(j,j) r1_1_bb(a,i) t0_1
    //              += +1.00 f_aa(k,k) l2_abab(i,j,a,b)
    //              += -0.50 <l,k||l,k>_abab l2_abab(i,j,a,b)
    //              += -0.50 <k,l||k,l>_abab l2_abab(i,j,a,b)
    //              += +0.25 <l,k||c,d>_aaaa t2_aaaa(c,d,l,k) l2_abab(i,j,a,b)
    //              += -0.50 <l,k||l,k>_bbbb l2_abab(i,j,a,b)
    //              += -0.50 <l,k||l,k>_aaaa l2_abab(i,j,a,b)
    //              += +0.25 <l,k||c,d>_bbbb t2_bbbb(c,d,l,k) l2_abab(i,j,a,b)
    //              += +1.00 f_bb(k,k) l2_abab(i,j,a,b)
    //              += -1.00 d-_bb(k,k) t0_1 l2_abab(i,j,a,b)
    //              += +0.25 <l,k||c,d>_abab t2_abab(c,d,l,k) l2_abab(i,j,a,b)
    //              += +0.25 <k,l||c,d>_abab t2_abab(c,d,k,l) l2_abab(i,j,a,b)
    //              += +0.25 <l,k||d,c>_abab t2_abab(d,c,l,k) l2_abab(i,j,a,b)
    //              += +0.25 <k,l||d,c>_abab t2_abab(d,c,k,l) l2_abab(i,j,a,b)
    //              += +1.00 d-_bb(j,i) r1_1_bb(a,j) t0_1
    //              += -1.00 f_bb(j,b) r1_bb(a,j) t1_1_bb(b,i)
    //              += -1.00 f_bb(j,b) r0_1 t2_bbbb(b,a,i,j)
    //              += -1.00 d-_bb(a,i) r0_1 t0_1
    //              += -0.50 <k,j||b,i>_abab r0_1 t2_abab(b,a,k,j)
    //              += -0.50 <j,k||b,i>_abab r0_1 t2_abab(b,a,j,k)
    //              += -0.50 <k,j||b,i>_bbbb r0_1 t2_bbbb(b,a,k,j)
    //              += +1.00 f_aa(j,b) r0_1 t2_abab(b,a,j,i)
    //              += +0.50 <j,a||b,c>_abab r0_1 t2_abab(b,c,j,i)
    //              += +0.50 <j,a||c,b>_abab r0_1 t2_abab(c,b,j,i)
    //              += -0.50 <j,a||b,c>_bbbb r0_1 t2_bbbb(b,c,i,j)
    //              += +1.00 <j,a||b,i>_abab r1_1_aa(b,j)
    //              += +1.00 f_bb(a,i) r0_1
    //              += -0.50 <k,j||b,c>_bbbb r1_1_bb(c,i) t2_bbbb(b,a,k,j)
    //              += -1.00 <j,k||b,i>_abab r1_bb(a,k) t1_1_aa(b,j)
    //              += -0.50 <j,k||b,c>_abab r1_bb(a,k) t2_1_abab(b,c,j,i)
    //              += -0.50 <j,k||c,b>_abab r1_bb(a,k) t2_1_abab(c,b,j,i)
    //              += +1.00 <j,a||c,b>_abab r1_aa(c,j) t1_1_bb(b,i)
    //              += +1.00 <j,k||b,c>_abab r1_1_bb(c,k) t2_abab(b,a,j,i)
    //              += +1.00 <k,j||b,c>_bbbb r2_bbbb(c,a,i,k) t1_1_bb(b,j)
    //              += -2.00 d-_bb(j,b) r1_1_bb(a,i) t1_1_bb(b,j)
    //              += -1.00 <k,j||c,b>_abab r1_aa(c,k) t2_1_bbbb(b,a,i,j)
    //              += -1.00 d+_aa(j,j) r1_bb(a,i)
    //              += -0.50 <k,j||b,c>_bbbb r1_bb(c,i) t2_1_bbbb(b,a,k,j)
    //              += +1.00 <j,a||b,c>_bbbb r1_bb(c,i) t1_1_bb(b,j)
    //              += -1.00 <k,j||b,c>_aaaa r1_aa(c,k) t2_1_abab(b,a,j,i)
    //              += +0.50 <j,a||b,c>_abab r2_1_abab(b,c,j,i)
    //              += +0.50 <j,a||c,b>_abab r2_1_abab(c,b,j,i)
    //              += -0.50 <k,j||b,i>_abab r2_1_abab(b,a,k,j)
    //              += -0.50 <j,k||b,i>_abab r2_1_abab(b,a,j,k)
    //              += -1.00 <k,j||b,c>_aaaa r2_abab(c,a,k,i) t1_1_aa(b,j)
    //              += +1.00 d-_bb(j,b) r2_1_bbbb(b,a,i,j) t0_1
    //              += +1.00 f_aa(j,b) r2_1_abab(b,a,j,i)
    //              += +1.00 <k,j||c,b>_abab r2_abab(c,a,k,i) t1_1_bb(b,j)
    //              += -1.00 <j,k||b,c>_abab r2_bbbb(c,a,i,k) t1_1_aa(b,j)
    //              += -1.00 f_bb(j,b) r2_1_bbbb(b,a,i,j)
    //              += +1.00 <k,j||b,c>_bbbb r1_bb(c,k) t2_1_bbbb(b,a,i,j)
    //              += -1.00 d-_bb(a,b) r1_1_bb(b,i) t0_1
    //              += -1.00 d-_aa(j,b) r2_1_abab(b,a,j,i) t0_1
    //              += -0.50 <k,j||c,b>_abab r2_abab(c,a,k,j) t1_1_bb(b,i)
    //              += -0.50 <j,k||c,b>_abab r2_abab(c,a,j,k) t1_1_bb(b,i)
    //              += -0.50 <j,a||b,c>_bbbb r2_1_bbbb(b,c,i,j)
    //              += -0.50 <j,k||b,c>_abab r1_1_bb(a,k) t2_abab(b,c,j,i)
    //              += -0.50 <j,k||c,b>_abab r1_1_bb(a,k) t2_abab(c,b,j,i)
    //              += +1.00 f_bb(a,b) r1_1_bb(b,i)
    //              += -1.00 <j,a||b,c>_bbbb r1_bb(c,j) t1_1_bb(b,i)
    //              += +1.00 r1_bb(a,i) t0_1 w0
    //              += +0.25 <k,j||b,c>_aaaa r1_bb(a,i) t2_1_aaaa(b,c,k,j)
    //              += +0.25 <k,j||b,c>_abab r1_bb(a,i) t2_1_abab(b,c,k,j)
    //              += +0.25 <j,k||b,c>_abab r1_bb(a,i) t2_1_abab(b,c,j,k)
    //              += +0.25 <k,j||c,b>_abab r1_bb(a,i) t2_1_abab(c,b,k,j)
    //              += +0.25 <j,k||c,b>_abab r1_bb(a,i) t2_1_abab(c,b,j,k)
    //              += -1.00 d-_aa(j,b) r1_bb(a,i) t0_1 t1_1_aa(b,j)
    //              += +1.00 f_bb(j,b) r1_bb(a,i) t1_1_bb(b,j)
    //              += -1.00 d-_bb(j,b) r1_bb(a,i) t0_1 t1_1_bb(b,j)
    //              += +0.25 <k,j||b,c>_bbbb r1_bb(a,i) t2_1_bbbb(b,c,k,j)
    //              += +1.00 f_aa(j,b) r1_bb(a,i) t1_1_aa(b,j)
    sigmar1_1_bb("R,a,i")  = -1.00 * tmps_["109_bb_Lvo"]("R,a,i");
    tmps_["109_bb_Lvo"].~TArrayD();

    // flops: o4v0L1  = o4v2L1
    //  mems: o4v0L1  = o4v0L1
    tmps_["110_bbbb_Loooo"]("R,i,j,k,l")  = r2_1["bbbb"]("R,a,b,i,j") * eri["bbbb_oovv"]("k,l,a,b");

    // flops: o2v2L1  = o4v2L1 o2v2 o2v2L1 o0v0L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o0v0L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o2v4L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2 o2v2L1 o0v0L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o0v0L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j")  = 0.25 * t2["bbbb"]("a,b,k,l") * tmps_["110_bbbb_Loooo"]("R,i,j,l,k");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") += 0.25 * (reused_["94_bbbb_oovv"]("i,j,a,b") + -2.00 * reused_["79_bbbb_vvoo"]("a,b,i,j")) * r0_1("R");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") += scalars_["2"] * r0_1("R") * t2_1["bbbb"]("a,b,i,j");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") -= 0.50 * r0_1("R") * reused_["92_bbbb_voov"]("b,i,j,a");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") += 0.50 * eri["bbbb_oooo"]("l,k,i,j") * r2_1["bbbb"]("R,a,b,k,l");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") += scalars_["1"] * r2["bbbb"]("R,a,b,i,j");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") += scalars_["2"] * r2["bbbb"]("R,a,b,i,j");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") -= scalars_["6"] * r2["bbbb"]("R,a,b,i,j");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") -= eri["bbbb_vvoo"]("a,b,i,j") * r0_1("R");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") += 0.25 * r2_1["bbbb"]("R,a,b,k,l") * reused_["36_bbbb_oooo"]("l,k,i,j");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") += scalars_["5"] * r2_1["bbbb"]("R,a,b,i,j");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") += 2.00 * scalars_["3"] * r2_1["bbbb"]("R,a,b,i,j");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") -= 0.50 * r0_1("R") * reused_["110_bbbb_voov"]("a,i,j,b");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") -= w0 * r2_1["bbbb"]("R,a,b,i,j");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") += 2.00 * scalars_["4"] * r2_1["bbbb"]("R,a,b,i,j");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") += scalars_["1"] * r0_1("R") * t2_1["bbbb"]("a,b,i,j");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") += 0.25 * r2["bbbb"]("R,a,b,k,l") * reused_["50_bbbb_oooo"]("i,j,l,k");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") -= 0.50 * eri["bbbb_vvvv"]("a,b,c,d") * r2_1["bbbb"]("R,c,d,i,j");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") += t2_1["bbbb"]("a,b,i,j") * tmps_["43_L"]("R");
    tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j") += 0.25 * t2_1["bbbb"]("a,b,k,l") * tmps_["47_bbbb_Loooo"]("R,i,j,l,k");
    tmps_["110_bbbb_Loooo"].~TArrayD();
    tmps_["47_bbbb_Loooo"].~TArrayD();

    // sigmar2_1_bbbb += -1.00 d-_bb(k,c) r1_1_bb(c,k) t2_1_bbbb(a,b,i,j)
    //                += -1.00 d-_aa(k,c) r1_1_aa(c,k) t2_1_bbbb(a,b,i,j)
    //                += +0.25 <l,k||c,d>_bbbb r0_1 t2_bbbb(a,b,l,k) t2_bbbb(c,d,i,j)
    //                += +0.50 <l,k||i,j>_bbbb r0_1 t2_bbbb(a,b,l,k)
    //                += +0.50 <a,b||c,d>_bbbb r0_1 t2_bbbb(c,d,i,j)
    //                += -1.00 d-_bb(k,k) r0_1 t2_1_bbbb(a,b,i,j)
    //                += -0.50 <l,k||c,d>_bbbb r0_1 t2_bbbb(c,a,l,k) t2_bbbb(d,b,i,j)
    //                += +0.50 <l,k||i,j>_bbbb r2_1_bbbb(a,b,l,k)
    //                += -1.00 d+_aa(k,k) r2_bbbb(a,b,i,j)
    //                += -1.00 d+_bb(k,k) r2_bbbb(a,b,i,j)
    //                += +1.00 r2_bbbb(a,b,i,j) t0_1 w0
    //                += +0.25 <l,k||c,d>_aaaa r2_bbbb(a,b,i,j) t2_1_aaaa(c,d,l,k)
    //                += +0.25 <l,k||c,d>_abab r2_bbbb(a,b,i,j) t2_1_abab(c,d,l,k)
    //                += +0.25 <k,l||c,d>_abab r2_bbbb(a,b,i,j) t2_1_abab(c,d,k,l)
    //                += +0.25 <l,k||d,c>_abab r2_bbbb(a,b,i,j) t2_1_abab(d,c,l,k)
    //                += +0.25 <k,l||d,c>_abab r2_bbbb(a,b,i,j) t2_1_abab(d,c,k,l)
    //                += -1.00 d-_aa(k,c) r2_bbbb(a,b,i,j) t0_1 t1_1_aa(c,k)
    //                += +1.00 f_bb(k,c) r2_bbbb(a,b,i,j) t1_1_bb(c,k)
    //                += -1.00 d-_bb(k,c) r2_bbbb(a,b,i,j) t0_1 t1_1_bb(c,k)
    //                += +0.25 <l,k||c,d>_bbbb r2_bbbb(a,b,i,j) t2_1_bbbb(c,d,l,k)
    //                += +1.00 f_aa(k,c) r2_bbbb(a,b,i,j) t1_1_aa(c,k)
    //                += +1.00 <a,b||i,j>_bbbb r0_1
    //                += +0.25 <l,k||c,d>_bbbb r2_1_bbbb(a,b,l,k) t2_bbbb(c,d,i,j)
    //                += -1.00 d-_aa(k,k) r2_1_bbbb(a,b,i,j) t0_1
    //            += +1.00 f_aa(j,j) l1_bb(i,a)
    //            += -0.50 <k,j||k,j>_abab l1_bb(i,a)
    //            += -0.50 <j,k||j,k>_abab l1_bb(i,a)
    //            += +0.25 <k,j||b,c>_aaaa t2_aaaa(b,c,k,j) l1_bb(i,a)
    //            += -0.50 <k,j||k,j>_bbbb l1_bb(i,a)
    //            += -0.50 <k,j||k,j>_aaaa l1_bb(i,a)
    //            += +0.25 <k,j||b,c>_bbbb t2_bbbb(b,c,k,j) l1_bb(i,a)
    //            += +1.00 f_bb(j,j) l1_bb(i,a)
    //            += -1.00 d-_bb(j,j) t0_1 l1_bb(i,a)
    //            += +0.25 <k,j||b,c>_abab t2_abab(b,c,k,j) l1_bb(i,a)
    //            += +0.25 <j,k||b,c>_abab t2_abab(b,c,j,k) l1_bb(i,a)
    //            += +0.25 <k,j||c,b>_abab t2_abab(c,b,k,j) l1_bb(i,a)
    //            += +0.25 <j,k||c,b>_abab t2_abab(c,b,j,k) l1_bb(i,a)
    //                += -2.00 d-_bb(k,c) r2_1_bbbb(a,b,i,j) t1_1_bb(c,k)
    //                += -0.50 <l,k||c,d>_bbbb r0_1 t2_bbbb(c,a,i,j) t2_bbbb(d,b,l,k)
    //                += +1.00 r2_1_bbbb(a,b,i,j) w0
    //                += -2.00 d-_aa(k,c) r2_1_bbbb(a,b,i,j) t1_1_aa(c,k)
    //                += -1.00 d-_aa(k,k) r0_1 t2_1_bbbb(a,b,i,j)
    //                += +0.25 <l,k||c,d>_bbbb r2_bbbb(a,b,l,k) t2_1_bbbb(c,d,i,j)
    //                += +0.50 <a,b||c,d>_bbbb r2_1_bbbb(c,d,i,j)
    //                += +0.25 <l,k||c,d>_bbbb r2_1_bbbb(c,d,i,j) t2_bbbb(a,b,l,k)
    //                += +0.25 <l,k||c,d>_bbbb r2_bbbb(c,d,i,j) t2_1_bbbb(a,b,l,k)
    sigmar2_1_bbbb("R,a,b,i,j") -= tmps_["111_bbbb_Lvvoo"]("R,a,b,i,j");
    tmps_["111_bbbb_Lvvoo"].~TArrayD();

    // flops: o0v2L1  = o2v3L1 o2v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["115_bb_Lvv"]("L,a,b")  = l2["bbbb"]("L,i,j,c,a") * t2["bbbb"]("b,c,i,j");
    tmps_["115_bb_Lvv"]("L,a,b") += l2_1["bbbb"]("L,i,j,c,a") * t2_1["bbbb"]("b,c,i,j");

    // flops: o0v2L1  = o2v3L1 o2v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["114_bb_Lvv"]("L,a,b")  = l2["abab"]("L,i,j,c,a") * t2["abab"]("c,b,i,j");
    tmps_["114_bb_Lvv"]("L,a,b") += l2_1["abab"]("L,i,j,c,a") * t2_1["abab"]("c,b,i,j");

    // flops: o3v1L1  = o3v2L1
    //  mems: o3v1L1  = o3v1L1
    tmps_["113_bbbb_Loovo"]("L,i,j,a,k")  = l2_1["bbbb"]("L,i,j,b,a") * t1_1["bb"]("b,k");

    // flops: o3v1L1  = o3v2L1
    //  mems: o3v1L1  = o3v1L1
    tmps_["112_bbbb_Loovo"]("L,i,j,a,k")  = l2["bbbb"]("L,i,j,b,a") * t1_1["bb"]("b,k");

    // flops: o2v2L1  = o2v3L1 o3v2L1 o3v2L1 o2v3L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o3v2L1 o3v2L1 o2v3L1 o2v3L1 o2v3L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["116_bbbb_Lvoov"]("L,a,i,j,b")  = -0.50 * tmps_["115_bb_Lvv"]("L,a,c") * eri["bbbb_oovv"]("i,j,b,c");
    tmps_["116_bbbb_Lvoov"]("L,a,i,j,b") -= dp["bb_ov"]("k,b") * tmps_["112_bbbb_Loovo"]("L,i,j,a,k");
    tmps_["116_bbbb_Lvoov"]("L,a,i,j,b") += f["bb_ov"]("k,b") * tmps_["113_bbbb_Loovo"]("L,i,j,a,k");
    tmps_["116_bbbb_Lvoov"]("L,a,i,j,b") += eri["bbbb_oovv"]("i,j,b,c") * tmps_["114_bb_Lvv"]("L,a,c");
    tmps_["116_bbbb_Lvoov"]("L,a,i,j,b") -= reused_["12_bb_ov"]("k,b") * tmps_["113_bbbb_Loovo"]("L,i,j,a,k");
    tmps_["116_bbbb_Lvoov"]("L,a,i,j,b") += dp["bb_vv"]("d,b") * l2_1["bbbb"]("L,i,j,d,a");
    tmps_["116_bbbb_Lvoov"]("L,a,i,j,b") += l1["bb"]("L,k,a") * eri["bbbb_oovo"]("i,j,b,k");
    tmps_["116_bbbb_Lvoov"]("L,a,i,j,b") += l1_1["bb"]("L,k,a") * reused_["111_bbbb_oovo"]("i,j,b,k");
    tmps_["116_bbbb_Lvoov"]("L,a,i,j,b") -= f["bb_vv"]("d,b") * l2["bbbb"]("L,i,j,d,a");
    tmps_["116_bbbb_Lvoov"]("L,a,i,j,b") += l2["bbbb"]("L,i,j,d,a") * reused_["34_bb_vv"]("d,b");
    tmps_["116_bbbb_Lvoov"]("L,a,i,j,b") += l2_1["bbbb"]("L,i,j,d,a") * reused_["108_bb_vv"]("d,b");
    tmps_["116_bbbb_Lvoov"]("L,a,i,j,b") += l2["bbbb"]("L,i,j,d,a") * reused_["3_bb_vv"]("b,d");
    tmps_["116_bbbb_Lvoov"]("L,a,i,j,b") += 0.50 * l2["bbbb"]("L,i,j,d,a") * reused_["109_bb_vv"]("b,d");
    tmps_["116_bbbb_Lvoov"]("L,a,i,j,b") -= l2_1["bbbb"]("L,i,j,d,a") * reused_["112_bb_vv"]("d,b");
    tmps_["112_bbbb_Loovo"].~TArrayD();

    // sigmal2_bbbb += +0.50 P(a,b) <i,j||d,a>_bbbb t2_1_abab(c,d,k,l) l2_1_abab(k,l,c,b)
    //              += +0.50 P(a,b) <i,j||d,a>_bbbb t2_1_abab(c,d,l,k) l2_1_abab(l,k,c,b)
    //              += +0.50 P(a,b) <i,j||d,a>_bbbb t2_abab(c,d,k,l) l2_abab(k,l,c,b)
    //              += +0.50 P(a,b) <i,j||d,a>_bbbb t2_abab(c,d,l,k) l2_abab(l,k,c,b)
    //              += -1.00 P(a,b) f_bb(k,a) t1_1_bb(c,k) l2_1_bbbb(i,j,c,b)
    //              += +1.00 P(a,b) d-_bb(k,a) t1_1_bb(c,k) l2_bbbb(i,j,c,b)
    //              += +1.00 P(a,b) d-_bb(k,a) t0_1 t1_1_bb(c,k) l2_1_bbbb(i,j,c,b)
    //              += -1.00 P(a,b) d-_bb(c,a) t0_1 l2_bbbb(i,j,c,b)
    //              += +1.00 P(a,b) <k,c||d,a>_abab t1_1_aa(d,k) l2_1_bbbb(i,j,c,b)
    //              += -0.50 P(a,b) <l,k||d,a>_abab t2_1_abab(d,c,l,k) l2_1_bbbb(i,j,c,b)
    //              += -0.50 P(a,b) <k,l||d,a>_abab t2_1_abab(d,c,k,l) l2_1_bbbb(i,j,c,b)
    //              += +1.00 P(a,b) f_bb(c,a) l2_bbbb(i,j,c,b)
    //              += -0.50 P(a,b) <l,k||d,a>_abab t2_abab(d,c,l,k) l2_bbbb(i,j,c,b)
    //              += -0.50 P(a,b) <k,l||d,a>_abab t2_abab(d,c,k,l) l2_bbbb(i,j,c,b)
    //              += -0.50 P(a,b) <l,k||d,a>_bbbb t2_bbbb(d,c,l,k) l2_bbbb(i,j,c,b)
    //              += +1.00 P(a,b) <k,c||d,a>_bbbb t1_1_bb(d,k) l2_1_bbbb(i,j,c,b)
    //              += -0.50 P(a,b) <l,k||d,a>_bbbb t2_1_bbbb(d,c,l,k) l2_1_bbbb(i,j,c,b)
    //              += +1.00 P(a,b) <i,j||c,a>_bbbb t1_1_bb(c,k) l1_1_bb(k,b)
    //              += -1.00 P(a,b) <i,j||a,k>_bbbb l1_bb(k,b)
    //              += -1.00 P(a,b) d+_bb(c,a) l2_1_bbbb(i,j,c,b)
    //              += -0.50 P(a,b) <i,j||d,a>_bbbb t2_1_bbbb(d,c,k,l) l2_1_bbbb(k,l,c,b)
    //              += -0.50 P(a,b) <i,j||d,a>_bbbb t2_bbbb(d,c,k,l) l2_bbbb(k,l,c,b)
    sigmal2_bbbb("L,a,b,i,j") -= tmps_["116_bbbb_Lvoov"]("L,b,i,j,a");
    sigmal2_bbbb("L,a,b,i,j") += tmps_["116_bbbb_Lvoov"]("L,a,i,j,b");
    tmps_["116_bbbb_Lvoov"].~TArrayD();

    // flops: o4v0L1  = o4v2L1
    //  mems: o4v0L1  = o4v0L1
    tmps_["117_bbbb_Loooo"]("L,i,j,k,l")  = l2_1["bbbb"]("L,i,j,a,b") * t2["bbbb"]("a,b,k,l");

    // flops: o2v2L1  = o4v2L1 0 o2v2L1 o4v2L1 o2v2L1 o2v2L1 o2v2L1 o4v2L1 o2v2L1 o2v4L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 0 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["118_bbbb_Lvvoo"]("L,a,b,i,j")  = eri["bbbb_oovv"]("k,l,a,b") * tmps_["117_bbbb_Loooo"]("L,i,j,l,k");
    tmps_["118_bbbb_Lvvoo"]("L,a,b,i,j") += 4.00 * (scalars_["5"] + -1.00 * w0) * l2_1["bbbb"]("L,i,j,a,b");
    tmps_["118_bbbb_Lvvoo"]("L,a,b,i,j") -= l2_1["bbbb"]("L,k,l,a,b") * reused_["36_bbbb_oooo"]("i,j,k,l");
    tmps_["118_bbbb_Lvvoo"]("L,a,b,i,j") += 4.00 * scalars_["2"] * l2["bbbb"]("L,i,j,a,b");
    tmps_["118_bbbb_Lvvoo"]("L,a,b,i,j") -= 2.00 * eri["bbbb_oooo"]("i,j,k,l") * l2_1["bbbb"]("L,k,l,a,b");
    tmps_["118_bbbb_Lvvoo"]("L,a,b,i,j") += 2.00 * eri["bbbb_vvvv"]("d,c,a,b") * l2_1["bbbb"]("L,i,j,c,d");
    tmps_["118_bbbb_Lvvoo"]("L,a,b,i,j") -= 4.00 * eri["bbbb_oovv"]("i,j,a,b") * l0_1("L");
    tmps_["118_bbbb_Lvvoo"]("L,a,b,i,j") += 8.00 * scalars_["4"] * l2_1["bbbb"]("L,i,j,a,b");
    tmps_["118_bbbb_Lvvoo"]("L,a,b,i,j") += 8.00 * scalars_["3"] * l2_1["bbbb"]("L,i,j,a,b");
    tmps_["118_bbbb_Lvvoo"]("L,a,b,i,j") += 4.00 * scalars_["1"] * l2["bbbb"]("L,i,j,a,b");
    tmps_["117_bbbb_Loooo"].~TArrayD();

    // sigmal2_1_bbbb  = +0.25 <l,k||a,b>_bbbb t2_bbbb(d,c,l,k) l2_1_bbbb(i,j,d,c)
    //                += -1.00 d-_aa(k,k) t0_1 l2_1_bbbb(i,j,a,b)
    //            += +1.00 f_aa(j,j) r1_bb(a,i)
    //            += -0.50 <k,j||k,j>_abab r1_bb(a,i)
    //            += -0.50 <j,k||j,k>_abab r1_bb(a,i)
    //            += +0.25 <k,j||b,c>_aaaa r1_bb(a,i) t2_aaaa(b,c,k,j)
    //            += -0.50 <k,j||k,j>_bbbb r1_bb(a,i)
    //            += -0.50 <k,j||k,j>_aaaa r1_bb(a,i)
    //            += +0.25 <k,j||b,c>_bbbb r1_bb(a,i) t2_bbbb(b,c,k,j)
    //            += +1.00 f_bb(j,j) r1_bb(a,i)
    //            += -1.00 d-_bb(j,j) r1_bb(a,i) t0_1
    //            += +0.25 <k,j||b,c>_abab r1_bb(a,i) t2_abab(b,c,k,j)
    //            += +0.25 <j,k||b,c>_abab r1_bb(a,i) t2_abab(b,c,j,k)
    //            += +0.25 <k,j||c,b>_abab r1_bb(a,i) t2_abab(c,b,k,j)
    //            += +0.25 <j,k||c,b>_abab r1_bb(a,i) t2_abab(c,b,j,k)
    //                += +1.00 l2_1_bbbb(i,j,a,b) w0
    //                += +0.25 <i,j||c,d>_bbbb t2_bbbb(c,d,k,l) l2_1_bbbb(k,l,a,b)
    //                += -1.00 d-_bb(k,k) l2_bbbb(i,j,a,b)
    //                += +0.50 <i,j||k,l>_bbbb l2_1_bbbb(k,l,a,b)
    //                += +0.50 <d,c||a,b>_bbbb l2_1_bbbb(i,j,d,c)
    //                += +1.00 <i,j||a,b>_bbbb l0_1
    //                += -2.00 d-_aa(k,c) t1_1_aa(c,k) l2_1_bbbb(i,j,a,b)
    //                += -2.00 d-_bb(k,c) t1_1_bb(c,k) l2_1_bbbb(i,j,a,b)
    //                += -1.00 d-_aa(k,k) l2_bbbb(i,j,a,b)
    sigmal2_1_bbbb("L,a,b,i,j")  = -0.25 * tmps_["118_bbbb_Lvvoo"]("L,a,b,i,j");
    tmps_["118_bbbb_Lvvoo"].~TArrayD();

    // flops: o1v1L1  = o2v1L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["120_aa_Lvo"]("R,a,i")  = r1["aa"]("R,a,j") * reused_["120_aa_oo"]("j,i");

    // sigmar1_aa += +1.00 d-_aa(j,b) r1_aa(a,j) t1_1_aa(b,i)
    sigmar1_aa("R,a,i") += tmps_["120_aa_Lvo"]("R,a,i");

    // flops: o1v1L1  = o1v1L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["122_aa_Lvo"]("R,a,i")  = r0_1("R") * t1_1["aa"]("a,i");

    // csigmar1_1_aa += +1.00 r0_1 t1_1_aa(a,i)
    csigmar1_1_aa("R,a,i") += tmps_["122_aa_Lvo"]("R,a,i");

    // flops: o2v0L1  = o2v1L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["121_aa_Loo"]("R,i,j")  = r1["aa"]("R,a,i") * reused_["13_aa_ov"]("j,a");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["119_bb_Lov"]("R,i,a")  = eri["abab_oovv"]("j,i,b,a") * r1_1["aa"]("R,b,j");

    // flops: o1v1L1  = o2v3L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o2v3L1 o2v2L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v2L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o1v2L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o1v2L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o3v2L1 o1v1L1 o3v2L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o0v0L1 o1v1L1 o1v1L1 o0v0L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o0v0L1 o1v1L1 o1v1L1 o0v0L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
    tmps_["123_aa_Lvo"]("R,a,i")  = -0.50 * eri["aaaa_vovv"]("a,j,d,c") * r2_1["aaaa"]("R,d,c,i,j");
    tmps_["123_aa_Lvo"]("R,a,i") -= f["aa_vo"]("a,i") * r0_1("R");
    tmps_["123_aa_Lvo"]("R,a,i") -= f["bb_ov"]("k,b") * r2_1["abab"]("R,a,b,i,k");
    tmps_["123_aa_Lvo"]("R,a,i") += 2.00 * scalars_["3"] * r1_1["aa"]("R,a,i");
    tmps_["123_aa_Lvo"]("R,a,i") -= 0.50 * r1_1["aa"]("R,c,i") * reused_["23_aa_vv"]("a,c");
    tmps_["123_aa_Lvo"]("R,a,i") -= eri["abab_vovv"]("a,k,d,e") * r2_1["abab"]("R,d,e,i,k");
    tmps_["123_aa_Lvo"]("R,a,i") -= r2_1["aaaa"]("R,d,a,i,j") * reused_["13_aa_ov"]("j,d");
    tmps_["123_aa_Lvo"]("R,a,i") += r2["aaaa"]("R,c,a,i,m") * reused_["9_aa_ov"]("m,c");
    tmps_["123_aa_Lvo"]("R,a,i") += r1["aa"]("R,a,j") * reused_["121_aa_oo"]("i,j");
    tmps_["123_aa_Lvo"]("R,a,i") -= scalars_["6"] * r1["aa"]("R,a,i");
    tmps_["123_aa_Lvo"]("R,a,i") += 2.00 * scalars_["4"] * r1_1["aa"]("R,a,i");
    tmps_["123_aa_Lvo"]("R,a,i") += r1["aa"]("R,c,i") * reused_["122_aa_vv"]("a,c");
    tmps_["123_aa_Lvo"]("R,a,i") += scalars_["5"] * r1_1["aa"]("R,a,i");
    tmps_["123_aa_Lvo"]("R,a,i") += r1_1["aa"]("R,c,m") * reused_["55_aaaa_voov"]("a,i,m,c");
    tmps_["123_aa_Lvo"]("R,a,i") -= r1["aa"]("R,c,i") * reused_["123_aa_vv"]("a,c");
    tmps_["123_aa_Lvo"]("R,a,i") += eri["abba_vovo"]("a,k,b,i") * r1_1["bb"]("R,b,k");
    tmps_["123_aa_Lvo"]("R,a,i") -= f["aa_vv"]("a,d") * r1_1["aa"]("R,d,i");
    tmps_["123_aa_Lvo"]("R,a,i") += scalars_["2"] * r1["aa"]("R,a,i");
    tmps_["123_aa_Lvo"]("R,a,i") += r1_1["bb"]("R,e,n") * reused_["6_aabb_voov"]("a,i,n,e");
    tmps_["123_aa_Lvo"]("R,a,i") += r1["aa"]("R,c,m") * reused_["113_aaaa_voov"]("a,i,m,c");
    tmps_["123_aa_Lvo"]("R,a,i") += r2["abab"]("R,a,e,m,k") * reused_["116_abba_oovo"]("m,k,e,i");
    tmps_["123_aa_Lvo"]("R,a,i") += r1_1["aa"]("R,c,i") * reused_["57_aa_vv"]("c,a");
    tmps_["123_aa_Lvo"]("R,a,i") += r1_1["aa"]("R,a,j") * f["aa_oo"]("j,i");
    tmps_["123_aa_Lvo"]("R,a,i") += r2_1["abab"]("R,a,b,i,k") * reused_["12_bb_ov"]("k,b");
    tmps_["123_aa_Lvo"]("R,a,i") -= 2.00 * r0_1("R") * reused_["89_aa_ov"]("i,a");
    tmps_["123_aa_Lvo"]("R,a,i") -= r2["abab"]("R,a,e,i,n") * reused_["8_bb_ov"]("n,e");
    tmps_["123_aa_Lvo"]("R,a,i") -= r1["aa"]("R,c,j") * reused_["117_aaaa_vovo"]("a,j,c,i");
    tmps_["123_aa_Lvo"]("R,a,i") -= r1_1["aa"]("R,a,j") * reused_["17_aa_oo"]("j,i");
    tmps_["123_aa_Lvo"]("R,a,i") += r1_1["aa"]("R,d,i") * reused_["59_aa_vv"]("a,d");
    tmps_["123_aa_Lvo"]("R,a,i") -= 0.50 * r1_1["aa"]("R,a,m") * reused_["15_aa_oo"]("i,m");
    tmps_["123_aa_Lvo"]("R,a,i") += r1["aa"]("R,a,m") * reused_["124_aa_oo"]("m,i");
    tmps_["123_aa_Lvo"]("R,a,i") += scalars_["1"] * r1["aa"]("R,a,i");
    tmps_["123_aa_Lvo"]("R,a,i") -= w0 * r1_1["aa"]("R,a,i");
    tmps_["123_aa_Lvo"]("R,a,i") -= 2.00 * r1_1["aa"]("R,a,j") * reused_["120_aa_oo"]("j,i");
    tmps_["123_aa_Lvo"]("R,a,i") -= 0.50 * r2_1["aaaa"]("R,d,a,m,j") * eri["aaaa_oovo"]("j,m,d,i");
    tmps_["123_aa_Lvo"]("R,a,i") -= r1["bb"]("R,e,k") * reused_["118_abba_vovo"]("a,k,e,i");
    tmps_["123_aa_Lvo"]("R,a,i") += r1_1["aa"]("R,a,m") * reused_["16_aa_oo"]("m,i");
    tmps_["123_aa_Lvo"]("R,a,i") += eri["aaaa_vovo"]("a,j,d,i") * r1_1["aa"]("R,d,j");
    tmps_["123_aa_Lvo"]("R,a,i") += r2["aaaa"]("R,c,a,i,m") * reused_["10_aa_ov"]("m,c");
    tmps_["123_aa_Lvo"]("R,a,i") -= eri["abba_oovo"]("m,k,b,i") * r2_1["abab"]("R,a,b,m,k");
    tmps_["123_aa_Lvo"]("R,a,i") += 0.50 * r2["aaaa"]("R,c,a,m,j") * reused_["119_aaaa_ooov"]("i,j,m,c");
    tmps_["123_aa_Lvo"]("R,a,i") -= r1_1["bb"]("R,e,n") * reused_["56_aabb_voov"]("a,i,n,e");
    tmps_["123_aa_Lvo"]("R,a,i") -= r1["bb"]("R,e,n") * reused_["115_aabb_voov"]("a,i,n,e");
    tmps_["123_aa_Lvo"]("R,a,i") -= 0.50 * r0_1("R") * reused_["87_aa_vo"]("a,i");
    tmps_["123_aa_Lvo"]("R,a,i") -= r1["aa"]("R,a,m") * reused_["125_aa_oo"]("m,i");
    tmps_["123_aa_Lvo"]("R,a,i") += r1["bb"]("R,e,n") * reused_["114_aabb_voov"]("a,i,n,e");
    tmps_["123_aa_Lvo"]("R,a,i") -= t0_1 * r0_1("R") * reused_["54_aa_vo"]("a,i");
    tmps_["123_aa_Lvo"]("R,a,i") += t0_1 * r0_1("R") * reused_["58_aa_vo"]("a,i");
    tmps_["123_aa_Lvo"]("R,a,i") -= r2["abab"]("R,a,e,i,n") * reused_["11_bb_ov"]("n,e");
    tmps_["123_aa_Lvo"]("R,a,i") += f["aa_ov"]("j,d") * r2_1["aaaa"]("R,d,a,i,j");
    tmps_["123_aa_Lvo"]("R,a,i") += t1_1["aa"]("a,i") * tmps_["43_L"]("R");
    tmps_["123_aa_Lvo"]("R,a,i") -= t1_1["aa"]("a,j") * tmps_["121_aa_Loo"]("R,i,j");
    tmps_["123_aa_Lvo"]("R,a,i") -= t2_1["abab"]("a,b,i,k") * tmps_["68_bb_Lov"]("R,k,b");
    tmps_["123_aa_Lvo"]("R,a,i") += scalars_["1"] * tmps_["122_aa_Lvo"]("R,a,i");
    tmps_["123_aa_Lvo"]("R,a,i") -= 2.00 * t1_1["aa"]("a,j") * tmps_["34_aa_Loo"]("R,i,j");
    tmps_["123_aa_Lvo"]("R,a,i") -= t0_1 * tmps_["120_aa_Lvo"]("R,a,i");
    tmps_["123_aa_Lvo"]("R,a,i") += scalars_["2"] * tmps_["122_aa_Lvo"]("R,a,i");
    tmps_["123_aa_Lvo"]("R,a,i") -= t2["abab"]("a,b,i,k") * tmps_["119_bb_Lov"]("R,k,b");
    tmps_["123_aa_Lvo"]("R,a,i") -= t1_1["aa"]("a,j") * tmps_["37_aa_Loo"]("R,j,i");
    tmps_["123_aa_Lvo"]("R,a,i") += 0.50 * t1_1["aa"]("a,j") * tmps_["36_aa_Loo"]("R,i,j");
    tmps_["119_bb_Lov"].~TArrayD();

    // sigmar1_1_aa  = -1.00 <k,j||c,b>_abab r2_aaaa(c,a,i,k) t1_1_bb(b,j)
    //              += +1.00 d-_aa(j,b) r2_1_aaaa(b,a,i,j) t0_1
    //              += -1.00 f_aa(j,b) r1_aa(a,j) t1_1_aa(b,i)
    //              += +1.00 r1_aa(a,i) t0_1 w0
    //              += +0.25 <k,j||b,c>_aaaa r1_aa(a,i) t2_1_aaaa(b,c,k,j)
    //              += +0.25 <k,j||b,c>_abab r1_aa(a,i) t2_1_abab(b,c,k,j)
    //              += +0.25 <j,k||b,c>_abab r1_aa(a,i) t2_1_abab(b,c,j,k)
    //              += +0.25 <k,j||c,b>_abab r1_aa(a,i) t2_1_abab(c,b,k,j)
    //              += +0.25 <j,k||c,b>_abab r1_aa(a,i) t2_1_abab(c,b,j,k)
    //              += -1.00 d-_aa(j,b) r1_aa(a,i) t0_1 t1_1_aa(b,j)
    //              += +1.00 f_bb(j,b) r1_aa(a,i) t1_1_bb(b,j)
    //              += -1.00 d-_bb(j,b) r1_aa(a,i) t0_1 t1_1_bb(b,j)
    //              += +0.25 <k,j||b,c>_bbbb r1_aa(a,i) t2_1_bbbb(b,c,k,j)
    //              += +1.00 f_aa(j,b) r1_aa(a,i) t1_1_aa(b,j)
    //              += +0.50 <a,j||b,c>_abab r2_1_abab(b,c,i,j)
    //              += +0.50 <a,j||c,b>_abab r2_1_abab(c,b,i,j)
    //              += -2.00 d-_aa(j,b) r1_1_aa(a,i) t1_1_aa(b,j)
    //              += +1.00 <j,a||b,c>_aaaa r1_aa(c,i) t1_1_aa(b,j)
    //              += -0.50 <k,j||b,c>_aaaa r1_aa(c,i) t2_1_aaaa(b,a,k,j)
    //              += -0.50 <k,j||b,c>_aaaa r1_1_aa(c,i) t2_aaaa(b,a,k,j)
    //              += -2.00 d-_bb(j,b) r1_1_aa(a,i) t1_1_bb(b,j)
    //              += +1.00 f_bb(j,b) r2_1_abab(a,b,i,j)
    //              += -1.00 d-_aa(j,j) r1_1_aa(a,i) t0_1
    //              += +1.00 f_aa(k,k) l2_bbbb(i,j,a,b)
    //              += -0.50 <l,k||l,k>_abab l2_bbbb(i,j,a,b)
    //              += -0.50 <k,l||k,l>_abab l2_bbbb(i,j,a,b)
    //              += +0.25 <l,k||c,d>_aaaa t2_aaaa(c,d,l,k) l2_bbbb(i,j,a,b)
    //              += -0.50 <l,k||l,k>_bbbb l2_bbbb(i,j,a,b)
    //              += -0.50 <l,k||l,k>_aaaa l2_bbbb(i,j,a,b)
    //              += +0.25 <l,k||c,d>_bbbb t2_bbbb(c,d,l,k) l2_bbbb(i,j,a,b)
    //              += +1.00 f_bb(k,k) l2_bbbb(i,j,a,b)
    //              += -1.00 d-_bb(k,k) t0_1 l2_bbbb(i,j,a,b)
    //              += +0.25 <l,k||c,d>_abab t2_abab(c,d,l,k) l2_bbbb(i,j,a,b)
    //              += +0.25 <k,l||c,d>_abab t2_abab(c,d,k,l) l2_bbbb(i,j,a,b)
    //              += +0.25 <l,k||d,c>_abab t2_abab(d,c,l,k) l2_bbbb(i,j,a,b)
    //              += +0.25 <k,l||d,c>_abab t2_abab(d,c,k,l) l2_bbbb(i,j,a,b)
    //              += +1.00 <k,j||b,c>_aaaa r1_1_aa(c,k) t2_aaaa(b,a,i,j)
    //              += +1.00 <a,j||c,b>_abab r1_aa(c,i) t1_1_bb(b,j)
    //              += -0.50 <k,j||c,b>_abab r1_aa(c,i) t2_1_abab(a,b,k,j)
    //              += -0.50 <j,k||c,b>_abab r1_aa(c,i) t2_1_abab(a,b,j,k)
    //              += +1.00 f_aa(a,i) r0_1
    //              += +1.00 <a,j||i,b>_abab r1_1_bb(b,j)
    //              += +1.00 f_aa(a,b) r1_1_aa(b,i)
    //              += -1.00 d+_bb(j,j) r1_aa(a,i)
    //              += -1.00 <j,k||b,c>_abab r1_1_bb(c,k) t2_aaaa(b,a,i,j)
    //              += +1.00 <k,j||b,c>_aaaa r1_aa(c,k) t2_1_aaaa(b,a,i,j)
    //              += -0.50 <k,j||b,c>_abab r2_abab(a,c,k,j) t1_1_aa(b,i)
    //              += -0.50 <j,k||b,c>_abab r2_abab(a,c,j,k) t1_1_aa(b,i)
    //              += -0.50 <k,j||c,b>_abab r1_1_aa(c,i) t2_abab(a,b,k,j)
    //              += -0.50 <j,k||c,b>_abab r1_1_aa(c,i) t2_abab(a,b,j,k)
    //              += -1.00 f_aa(j,i) r1_1_aa(a,j)
    //              += -1.00 d-_bb(j,b) r2_1_abab(a,b,i,j) t0_1
    //              += +2.00 d-_aa(j,i) r0_1 t1_1_aa(a,j)
    //              += +2.00 d-_aa(j,b) r0_1 t2_1_aaaa(b,a,i,j)
    //              += -2.00 d-_aa(a,b) r0_1 t1_1_aa(b,i)
    //              += -2.00 d-_bb(j,b) r0_1 t2_1_abab(a,b,i,j)
    //              += -0.50 <j,a||b,c>_aaaa r2_1_aaaa(b,c,i,j)
    //              += +1.00 <j,k||b,c>_abab r2_abab(a,c,i,k) t1_1_aa(b,j)
    //              += -1.00 <j,a||b,c>_aaaa r1_aa(c,j) t1_1_aa(b,i)
    //              += +1.00 d-_aa(j,i) r1_1_aa(a,j) t0_1
    //              += -1.00 d-_aa(a,b) r1_1_aa(b,i) t0_1
    //              += -0.50 <k,j||b,c>_aaaa r1_1_aa(a,k) t2_aaaa(b,c,i,j)
    //              += +1.00 <k,j||b,i>_aaaa r1_aa(a,k) t1_1_aa(b,j)
    //              += -0.50 <k,j||b,c>_aaaa r1_aa(a,k) t2_1_aaaa(b,c,i,j)
    //              += -1.00 d+_aa(j,j) r1_aa(a,i)
    //              += +1.00 r1_1_aa(a,i) w0
    //              += +2.00 d-_aa(j,b) r1_1_aa(a,j) t1_1_aa(b,i)
    //              += -0.50 <k,j||b,i>_aaaa r2_1_aaaa(b,a,k,j)
    //              += +1.00 <a,j||b,c>_abab r1_bb(c,j) t1_1_aa(b,i)
    //              += -0.50 <k,j||b,c>_abab r1_1_aa(a,k) t2_abab(b,c,i,j)
    //              += -0.50 <k,j||c,b>_abab r1_1_aa(a,k) t2_abab(c,b,i,j)
    //              += +1.00 <j,a||b,i>_aaaa r1_1_aa(b,j)
    //              += +1.00 <k,j||b,c>_aaaa r2_aaaa(c,a,i,k) t1_1_aa(b,j)
    //              += -0.50 <k,j||i,b>_abab r2_1_abab(a,b,k,j)
    //              += -0.50 <j,k||i,b>_abab r2_1_abab(a,b,j,k)
    //              += +0.50 <k,j||b,c>_aaaa r2_aaaa(c,a,k,j) t1_1_aa(b,i)
    //              += -1.00 <k,j||b,c>_bbbb r1_1_bb(c,k) t2_abab(a,b,i,j)
    //              += -1.00 <k,j||b,c>_bbbb r1_bb(c,k) t2_1_abab(a,b,i,j)
    //              += -0.50 <j,a||b,c>_aaaa r0_1 t2_aaaa(b,c,i,j)
    //              += -0.50 <k,j||i,b>_abab r0_1 t2_abab(a,b,k,j)
    //              += -0.50 <j,k||i,b>_abab r0_1 t2_abab(a,b,j,k)
    //              += -0.50 <k,j||b,i>_aaaa r0_1 t2_aaaa(b,a,k,j)
    //              += -1.00 f_aa(j,b) r0_1 t2_aaaa(b,a,i,j)
    //              += +0.50 <a,j||b,c>_abab r0_1 t2_abab(b,c,i,j)
    //              += +0.50 <a,j||c,b>_abab r0_1 t2_abab(c,b,i,j)
    //              += -1.00 d-_aa(a,i) r0_1 t0_1
    //              += +1.00 f_bb(j,b) r0_1 t2_abab(a,b,i,j)
    //              += -1.00 <k,j||i,b>_abab r1_aa(a,k) t1_1_bb(b,j)
    //              += -0.50 <k,j||b,c>_abab r1_aa(a,k) t2_1_abab(b,c,i,j)
    //              += -0.50 <k,j||c,b>_abab r1_aa(a,k) t2_1_abab(c,b,i,j)
    //              += -1.00 <j,k||b,c>_abab r1_bb(c,k) t2_1_aaaa(b,a,i,j)
    //              += +1.00 d-_aa(j,b) r0_1 t2_aaaa(b,a,i,j) t0_1
    //              += -1.00 d-_bb(j,b) r0_1 t2_abab(a,b,i,j) t0_1
    //              += -1.00 <k,j||b,c>_bbbb r2_abab(a,c,i,k) t1_1_bb(b,j)
    //              += -1.00 f_aa(j,b) r2_1_aaaa(b,a,i,j)
    //              += -1.00 d-_bb(j,b) r1_1_bb(b,j) t1_1_aa(a,i)
    //              += -1.00 d-_aa(j,b) r1_1_aa(b,j) t1_1_aa(a,i)
    //              += +1.00 d-_aa(j,b) r1_aa(b,i) t0_1 t1_1_aa(a,j)
    //              += +1.00 <k,j||c,b>_abab r1_aa(c,k) t2_1_abab(a,b,i,j)
    //              += -1.00 d-_aa(j,j) r0_1 t1_1_aa(a,i)
    //              += +2.00 d-_aa(j,b) r1_1_aa(b,i) t1_1_aa(a,j)
    //              += +1.00 d-_aa(j,b) r1_aa(a,j) t0_1 t1_1_aa(b,i)
    //              += -1.00 d-_bb(j,j) r0_1 t1_1_aa(a,i)
    //              += +1.00 <k,j||c,b>_abab r1_1_aa(c,k) t2_abab(a,b,i,j)
    //              += -1.00 <j,k||i,b>_abab r1_bb(b,k) t1_1_aa(a,j)
    //              += -1.00 <k,j||b,i>_aaaa r1_aa(b,k) t1_1_aa(a,j)
    //              += +0.50 <k,j||b,c>_aaaa r2_aaaa(b,c,i,k) t1_1_aa(a,j)
    //              += -1.00 f_aa(j,b) r1_aa(b,i) t1_1_aa(a,j)
    //              += -0.50 <j,k||b,c>_abab r2_abab(b,c,i,k) t1_1_aa(a,j)
    //              += -0.50 <j,k||c,b>_abab r2_abab(c,b,i,k) t1_1_aa(a,j)
    sigmar1_1_aa("R,a,i")  = -1.00 * tmps_["123_aa_Lvo"]("R,a,i");
    tmps_["123_aa_Lvo"].~TArrayD();

    // flops: o1v1L1  = o2v1L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["124_aa_Lvo"]("R,a,i")  = 0.50 * t1_1["aa"]("a,j") * tmps_["36_aa_Loo"]("R,i,j");
    tmps_["124_aa_Lvo"].~TArrayD();

    // flops: o2v2L1  = o3v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["125_aaaa_Lvvoo"]("R,a,b,i,j")  = 0.50 * t2_1["aaaa"]("a,b,i,k") * tmps_["36_aa_Loo"]("R,j,k");
    tmps_["125_aaaa_Lvvoo"].~TArrayD();

    // flops: o2v2L1  = o2v2L1 o3v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1
    tmps_["126_aaaa_Lovvo"]("R,i,a,b,j")  = r0_1("R") * reused_["75_aaaa_ovvo"]("i,a,b,j");
    tmps_["126_aaaa_Lovvo"]("R,i,a,b,j") += r2["aaaa"]("R,a,b,j,k") * reused_["120_aa_oo"]("k,i");

    // sigmar2_aaaa += +1.00 P(i,j) d-_aa(k,j) r0_1 t2_aaaa(a,b,i,k)
    //              += +1.00 P(i,j) d-_aa(k,c) r2_aaaa(a,b,i,k) t1_1_aa(c,j)
    sigmar2_aaaa("R,a,b,i,j") += tmps_["126_aaaa_Lovvo"]("R,j,a,b,i");
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["126_aaaa_Lovvo"]("R,i,a,b,j");

    // flops: o2v0L1  = o3v1L1 o3v1L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
    tmps_["129_aa_Loo"]("R,i,j")  = eri["aaaa_oovo"]("i,k,a,j") * r1_1["aa"]("R,a,k");
    tmps_["129_aa_Loo"]("R,i,j") += eri["abba_oovo"]("i,l,b,j") * r1_1["bb"]("R,b,l");

    // flops: o2v0L1  = o2v1L1 o3v2L1 o2v0L1 o3v2L1 o2v0L1 o2v1L1 o2v0L1 o3v1L1 o2v0L1 o2v1L1 o2v0L1 o2v1L1 o2v0L1 o3v1L1 o2v0L1
    //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1 o2v0L1
    tmps_["128_aa_Loo"]("R,i,j")  = -1.00 * r1["aa"]("R,a,i") * reused_["128_aa_ov"]("j,a");
    tmps_["128_aa_Loo"]("R,i,j") += eri["abab_oovv"]("j,l,b,d") * r2_1["abab"]("R,b,d,i,l");
    tmps_["128_aa_Loo"]("R,i,j") += 0.50 * eri["aaaa_oovv"]("j,k,b,a") * r2_1["aaaa"]("R,b,a,i,k");
    tmps_["128_aa_Loo"]("R,i,j") -= r1_1["aa"]("R,b,i") * reused_["13_aa_ov"]("j,b");
    tmps_["128_aa_Loo"]("R,i,j") += r1["bb"]("R,d,l") * reused_["116_abba_oovo"]("j,l,d,i");
    tmps_["128_aa_Loo"]("R,i,j") += f["aa_ov"]("j,b") * r1_1["aa"]("R,b,i");
    tmps_["128_aa_Loo"]("R,i,j") += r1["aa"]("R,a,i") * reused_["9_aa_ov"]("j,a");
    tmps_["128_aa_Loo"]("R,i,j") += r1["aa"]("R,a,k") * reused_["119_aaaa_ooov"]("i,j,k,a");

    // flops: o2v2L1  = o3v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["127_aaaa_Lvvoo"]("R,a,b,i,j")  = 2.00 * r2_1["aaaa"]("R,a,b,i,k") * reused_["17_aa_oo"]("k,j");

    // flops: o2v2L1  = o3v2L1 o3v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o4v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j")  = -1.00 * t2_1["aaaa"]("a,b,i,k") * tmps_["121_aa_Loo"]("R,j,k");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") -= 2.00 * t2_1["aaaa"]("a,b,i,k") * tmps_["34_aa_Loo"]("R,j,k");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") += t2["aaaa"]("a,b,i,k") * tmps_["128_aa_Loo"]("R,j,k");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") += t2_1["aaaa"]("a,b,j,k") * tmps_["37_aa_Loo"]("R,k,i");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") += 0.50 * r0_1("R") * reused_["96_aaaa_ovvo"]("i,a,b,j");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") -= 0.50 * r1["aa"]("R,c,j") * reused_["127_aaaa_vovv"]("c,i,a,b");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") -= f["aa_oo"]("k,i") * r2_1["aaaa"]("R,a,b,j,k");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") += eri["aaaa_vvvo"]("a,b,c,i") * r1_1["aa"]("R,c,j");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") += 2.00 * r2_1["aaaa"]("R,a,b,j,k") * reused_["120_aa_oo"]("k,i");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") -= 0.50 * r2["aaaa"]("R,a,b,l,k") * reused_["60_aaaa_oooo"]("k,l,i,j");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") -= r2_1["aaaa"]("R,a,b,j,l") * reused_["16_aa_oo"]("l,i");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") -= reused_["126_aaaa_vvvo"]("a,b,e,i") * r1["aa"]("R,e,j");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") -= r2["aaaa"]("R,a,b,j,l") * reused_["124_aa_oo"]("l,i");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") += 0.50 * r2_1["aaaa"]("R,a,b,j,l") * reused_["15_aa_oo"]("i,l");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") += dp["aa_oo"]("k,i") * r2["aaaa"]("R,a,b,j,k");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") -= 0.50 * r1_1["aa"]("R,c,j") * reused_["14_aaaa_vovv"]("c,i,a,b");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") -= r2["aaaa"]("R,a,b,j,k") * reused_["121_aa_oo"]("i,k");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") += r2["aaaa"]("R,a,b,j,l") * reused_["125_aa_oo"]("l,i");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") += 0.50 * tmps_["127_aaaa_Lvvoo"]("R,a,b,j,i");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") += t0_1 * tmps_["126_aaaa_Lovvo"]("R,i,a,b,j");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") += t2["aaaa"]("a,b,j,k") * tmps_["129_aa_Loo"]("R,k,i");
    tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j") += 0.50 * t2_1["aaaa"]("a,b,i,k") * tmps_["36_aa_Loo"]("R,j,k");
    tmps_["127_aaaa_Lvvoo"].~TArrayD();
    tmps_["126_aaaa_Lovvo"].~TArrayD();

    // sigmar2_1_aaaa += +0.50 P(i,j) <k,l||c,d>_abab r2_1_abab(c,d,i,l) t2_aaaa(a,b,j,k)
    //                += +0.50 P(i,j) <k,l||d,c>_abab r2_1_abab(d,c,i,l) t2_aaaa(a,b,j,k)
    //                += +1.00 P(i,j) <l,k||c,d>_aaaa r1_aa(d,i) t2_aaaa(a,b,j,k) t1_1_aa(c,l)
    //                += -0.50 P(i,j) <l,k||c,d>_aaaa r2_1_aaaa(c,d,i,l) t2_aaaa(a,b,j,k)
    //                += -1.00 P(i,j) d-_aa(k,c) r1_1_aa(c,i) t2_aaaa(a,b,j,k) t0_1
    //                += +1.00 P(i,j) <k,l||c,d>_abab r1_bb(d,l) t2_aaaa(a,b,j,k) t1_1_aa(c,i)
    //                += +1.00 P(i,j) f_aa(k,c) r1_1_aa(c,i) t2_aaaa(a,b,j,k)
    //                += +1.00 P(i,j) <k,l||d,c>_abab r1_aa(d,i) t2_aaaa(a,b,j,k) t1_1_bb(c,l)
    //                += -1.00 P(i,j) <l,k||c,d>_aaaa r1_aa(d,l) t2_aaaa(a,b,j,k) t1_1_aa(c,i)
    //                += -2.00 P(i,j) d-_aa(k,c) r1_1_aa(c,i) t2_1_aaaa(a,b,j,k)
    //                += -1.00 P(i,j) <k,l||j,c>_abab r1_bb(c,l) t2_1_aaaa(a,b,i,k)
    //                += -1.00 P(i,j) <l,k||c,j>_aaaa r1_aa(c,l) t2_1_aaaa(a,b,i,k)
    //                += -1.00 P(i,j) d-_aa(k,c) r1_aa(c,i) t0_1 t2_1_aaaa(a,b,j,k)
    //                += -0.50 P(i,j) <l,k||c,d>_aaaa r0_1 t2_aaaa(a,b,i,l) t2_aaaa(c,d,j,k)
    //                += +1.00 P(i,j) <l,k||c,d>_bbbb r0_1 t2_abab(a,c,j,k) t2_abab(b,d,i,l)
    //                += -0.50 P(i,j) <l,k||c,d>_abab r0_1 t2_aaaa(a,b,i,l) t2_abab(c,d,j,k)
    //                += -0.50 P(i,j) <l,k||d,c>_abab r0_1 t2_aaaa(a,b,i,l) t2_abab(d,c,j,k)
    //                += -2.00 P(i,j) d-_aa(k,c) r0_1 t2_aaaa(a,b,j,k) t1_1_aa(c,i)
    //                += +2.00 P(i,j) d-_aa(k,j) r0_1 t2_1_aaaa(a,b,i,k)
    //                += -1.00 P(i,j) f_aa(k,j) r0_1 t2_aaaa(a,b,i,k)
    //                += +1.00 P(i,j) <l,k||c,d>_aaaa r0_1 t2_aaaa(c,a,j,k) t2_aaaa(d,b,i,l)
    //                += +0.50 P(i,j) <l,k||c,j>_aaaa r1_aa(c,i) t2_1_aaaa(a,b,l,k)
    //                += -1.00 P(i,j) f_aa(k,j) r2_1_aaaa(a,b,i,k)
    //                += +1.00 P(i,j) <a,b||c,j>_aaaa r1_1_aa(c,i)
    //                += +2.00 P(i,j) d-_aa(k,c) r2_1_aaaa(a,b,i,k) t1_1_aa(c,j)
    //                += +0.50 P(i,j) <l,k||c,j>_aaaa r2_aaaa(a,b,l,k) t1_1_aa(c,i)
    //                += -0.50 P(i,j) <l,k||c,d>_abab r2_1_aaaa(a,b,i,l) t2_abab(c,d,j,k)
    //                += -0.50 P(i,j) <l,k||d,c>_abab r2_1_aaaa(a,b,i,l) t2_abab(d,c,j,k)
    //                += -1.00 P(i,j) <a,b||c,d>_aaaa r1_aa(d,i) t1_1_aa(c,j)
    //                += -0.50 P(i,j) <l,k||c,d>_aaaa r1_aa(d,i) t2_aaaa(a,b,l,k) t1_1_aa(c,j)
    //                += +1.00 P(i,j) <l,k||c,j>_aaaa r2_aaaa(a,b,i,l) t1_1_aa(c,k)
    //                += -0.50 P(i,j) <l,k||c,d>_aaaa r2_aaaa(a,b,i,l) t2_1_aaaa(c,d,j,k)
    //                += -0.50 P(i,j) <l,k||c,d>_aaaa r2_1_aaaa(a,b,i,l) t2_aaaa(c,d,j,k)
    //                += +1.00 P(i,j) d+_aa(k,j) r2_aaaa(a,b,i,k)
    //                += +0.50 P(i,j) <l,k||c,j>_aaaa r1_1_aa(c,i) t2_aaaa(a,b,l,k)
    //                += -1.00 P(i,j) f_aa(k,c) r2_aaaa(a,b,i,k) t1_1_aa(c,j)
    //                += -1.00 P(i,j) <l,k||j,c>_abab r2_aaaa(a,b,i,l) t1_1_bb(c,k)
    //                += -0.50 P(i,j) <l,k||c,d>_abab r2_aaaa(a,b,i,l) t2_1_abab(c,d,j,k)
    //                += -0.50 P(i,j) <l,k||d,c>_abab r2_aaaa(a,b,i,l) t2_1_abab(d,c,j,k)
    //                += +1.00 P(i,j) d-_aa(k,j) r2_1_aaaa(a,b,i,k) t0_1
    //                += +1.00 P(i,j) d-_aa(k,j) r0_1 t2_aaaa(a,b,i,k) t0_1
    //                += +1.00 P(i,j) d-_aa(k,c) r2_aaaa(a,b,i,k) t0_1 t1_1_aa(c,j)
    //                += -1.00 P(i,j) <l,k||c,j>_aaaa r1_1_aa(c,l) t2_aaaa(a,b,i,k)
    //                += -1.00 P(i,j) <k,l||j,c>_abab r1_1_bb(c,l) t2_aaaa(a,b,i,k)
    //                += -0.50 P(i,j) <l,k||c,d>_aaaa r2_aaaa(c,d,i,l) t2_1_aaaa(a,b,j,k)
    //                += +1.00 P(i,j) f_aa(k,c) r1_aa(c,i) t2_1_aaaa(a,b,j,k)
    //                += +0.50 P(i,j) <k,l||c,d>_abab r2_abab(c,d,i,l) t2_1_aaaa(a,b,j,k)
    //                += +0.50 P(i,j) <k,l||d,c>_abab r2_abab(d,c,i,l) t2_1_aaaa(a,b,j,k)
    sigmar2_1_aaaa("R,a,b,i,j") += tmps_["130_aaaa_Lvvoo"]("R,a,b,j,i");
    sigmar2_1_aaaa("R,a,b,i,j") -= tmps_["130_aaaa_Lvvoo"]("R,a,b,i,j");
    tmps_["130_aaaa_Lvvoo"].~TArrayD();

    // flops: o2v2L1  = o3v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["131_aaaa_Lvoov"]("R,a,i,j,b")  = reused_["82_aaaa_vooo"]("a,i,j,k") * r1["aa"]("R,b,k");

    // sigmar2_1_aaaa += -1.00 P(a,b) d+_aa(k,c) r1_aa(b,k) t2_aaaa(c,a,i,j)
    sigmar2_1_aaaa("R,a,b,i,j") -= tmps_["131_aaaa_Lvoov"]("R,a,i,j,b");
    sigmar2_1_aaaa("R,a,b,i,j") += tmps_["131_aaaa_Lvoov"]("R,b,i,j,a");

    // flops: o0v2L1  = o2v3L1 o2v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["135_aa_Lvv"]("R,a,b")  = r2["abab"]("R,a,d,i,k") * reused_["129_abab_oovv"]("i,k,b,d");
    tmps_["135_aa_Lvv"]("R,a,b") += eri["aaaa_oovv"]("j,i,b,c") * r2["aaaa"]("R,c,a,i,j");

    // flops: o0v2L1  = o1v3L1 o1v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["134_aa_Lvv"]("R,a,b")  = eri["aaaa_vovv"]("a,i,b,c") * r1["aa"]("R,c,i");
    tmps_["134_aa_Lvv"]("R,a,b") += eri["abab_vovv"]("a,j,b,d") * r1["bb"]("R,d,j");

    // flops: o2v2L1  = o2v2L1
    //  mems: o2v2L1  = o2v2L1
    tmps_["132_aaaa_Lvoov"]("R,a,i,j,b")  = -2.00 * t0_1 * tmps_["131_aaaa_Lvoov"]("R,a,i,j,b");
    tmps_["132_aaaa_Lvoov"].~TArrayD();

    // flops: o3v1L1  = o3v2L1
    //  mems: o3v1L1  = o3v1L1
    tmps_["133_aaaa_Lvooo"]("R,a,i,j,k")  = r2["aaaa"]("R,b,a,i,j") * dp["aa_ov"]("k,b");

    // flops: o2v2L1  = o2v2L1 o3v2L1 o2v3L1 o2v2L1 o3v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["136_aaaa_Lvoov"]("R,a,i,j,b")  = -2.00 * t0_1 * tmps_["131_aaaa_Lvoov"]("R,a,i,j,b");
    tmps_["136_aaaa_Lvoov"]("R,a,i,j,b") += 2.00 * t1_1["aa"]("a,k") * tmps_["133_aaaa_Lvooo"]("R,b,i,j,k");
    tmps_["136_aaaa_Lvoov"]("R,a,i,j,b") += t2["aaaa"]("c,a,i,j") * tmps_["135_aa_Lvv"]("R,b,c");
    tmps_["136_aaaa_Lvoov"]("R,a,i,j,b") -= 2.00 * eri["aaaa_vooo"]("a,k,i,j") * r1["aa"]("R,b,k");
    tmps_["136_aaaa_Lvoov"]("R,a,i,j,b") -= 2.00 * reused_["57_aa_vv"]("d,a") * r2["aaaa"]("R,d,b,i,j");
    tmps_["136_aaaa_Lvoov"]("R,a,i,j,b") -= 2.00 * reused_["59_aa_vv"]("a,c") * r2["aaaa"]("R,c,b,i,j");
    tmps_["136_aaaa_Lvoov"]("R,a,i,j,b") += 2.00 * f["aa_vv"]("a,c") * r2["aaaa"]("R,c,b,i,j");
    tmps_["136_aaaa_Lvoov"]("R,a,i,j,b") += 2.00 * reused_["130_aaaa_vooo"]("a,i,j,k") * r1["aa"]("R,b,k");
    tmps_["136_aaaa_Lvoov"]("R,a,i,j,b") += reused_["23_aa_vv"]("a,d") * r2["aaaa"]("R,d,b,i,j");
    tmps_["136_aaaa_Lvoov"]("R,a,i,j,b") -= 2.00 * dp["aa_vv"]("a,c") * r2_1["aaaa"]("R,c,b,i,j");
    tmps_["136_aaaa_Lvoov"]("R,a,i,j,b") += 2.00 * tmps_["134_aa_Lvv"]("R,a,c") * t2["aaaa"]("c,b,i,j");
    tmps_["133_aaaa_Lvooo"].~TArrayD();
    tmps_["131_aaaa_Lvoov"].~TArrayD();

    // sigmar2_aaaa += -0.50 P(a,b) <l,k||c,d>_aaaa r2_aaaa(d,b,l,k) t2_aaaa(c,a,i,j)
    //              += +0.50 P(a,b) <l,k||c,d>_abab r2_abab(b,d,l,k) t2_aaaa(c,a,i,j)
    //              += +0.50 P(a,b) <k,l||c,d>_abab r2_abab(b,d,k,l) t2_aaaa(c,a,i,j)
    //              += +1.00 P(a,b) d-_aa(k,c) r2_aaaa(c,b,i,j) t1_1_aa(a,k)
    //              += +1.00 P(a,b) <k,a||i,j>_aaaa r1_aa(b,k)
    //              += -0.50 P(a,b) <l,k||d,c>_abab r2_aaaa(d,b,i,j) t2_abab(a,c,l,k)
    //              += -0.50 P(a,b) <k,l||d,c>_abab r2_aaaa(d,b,i,j) t2_abab(a,c,k,l)
    //              += -1.00 P(a,b) d-_aa(a,c) r2_aaaa(c,b,i,j) t0_1
    //              += +1.00 P(a,b) f_aa(a,c) r2_aaaa(c,b,i,j)
    //              += +1.00 P(a,b) f_aa(k,c) r1_aa(b,k) t2_aaaa(c,a,i,j)
    //              += +0.50 P(a,b) <k,a||c,d>_aaaa r1_aa(b,k) t2_aaaa(c,d,i,j)
    //              += -0.50 P(a,b) <l,k||c,d>_aaaa r2_aaaa(d,b,i,j) t2_aaaa(c,a,l,k)
    //              += -1.00 P(a,b) d-_aa(a,c) r2_1_aaaa(c,b,i,j)
    //              += -1.00 P(a,b) <k,a||c,d>_aaaa r1_aa(d,k) t2_aaaa(c,b,i,j)
    //              += +1.00 P(a,b) <a,k||c,d>_abab r1_bb(d,k) t2_aaaa(c,b,i,j)
    //              += -1.00 P(a,b) d-_aa(k,c) r1_aa(b,k) t2_aaaa(c,a,i,j) t0_1
    sigmar2_aaaa("R,a,b,i,j") += 0.50 * tmps_["136_aaaa_Lvoov"]("R,a,i,j,b");
    sigmar2_aaaa("R,a,b,i,j") -= 0.50 * tmps_["136_aaaa_Lvoov"]("R,b,i,j,a");
    tmps_["136_aaaa_Lvoov"].~TArrayD();

    // flops: o0v2L1  = o2v3L1 o2v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["140_aa_Lvv"]("L,a,b")  = l2["abab"]("L,i,j,a,c") * t2["abab"]("b,c,i,j");
    tmps_["140_aa_Lvv"]("L,a,b") += l2_1["abab"]("L,i,j,a,c") * t2_1["abab"]("b,c,i,j");

    // flops: o0v2L1  = o2v3L1 o2v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["139_aa_Lvv"]("L,a,b")  = l2["aaaa"]("L,i,j,c,a") * t2["aaaa"]("b,c,i,j");
    tmps_["139_aa_Lvv"]("L,a,b") += l2_1["aaaa"]("L,i,j,c,a") * t2_1["aaaa"]("b,c,i,j");

    // flops: o3v1L1  = o3v2L1
    //  mems: o3v1L1  = o3v1L1
    tmps_["138_aaaa_Loovo"]("L,i,j,a,k")  = l2_1["aaaa"]("L,i,j,b,a") * t1_1["aa"]("b,k");

    // flops: o3v1L1  = o3v2L1
    //  mems: o3v1L1  = o3v1L1
    tmps_["137_aaaa_Loovo"]("L,i,j,a,k")  = l2["aaaa"]("L,i,j,b,a") * t1_1["aa"]("b,k");

    // flops: o2v2L1  = o3v2L1 o2v3L1 o2v3L1 o3v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["141_aaaa_Lvoov"]("L,a,i,j,b")  = -1.00 * f["aa_ov"]("k,a") * tmps_["138_aaaa_Loovo"]("L,i,j,b,k");
    tmps_["141_aaaa_Lvoov"]("L,a,i,j,b") -= eri["aaaa_oovv"]("i,j,a,c") * tmps_["140_aa_Lvv"]("L,b,c");
    tmps_["141_aaaa_Lvoov"]("L,a,i,j,b") += 0.50 * eri["aaaa_oovv"]("i,j,a,c") * tmps_["139_aa_Lvv"]("L,b,c");
    tmps_["141_aaaa_Lvoov"]("L,a,i,j,b") += dp["aa_ov"]("k,a") * tmps_["137_aaaa_Loovo"]("L,i,j,b,k");
    tmps_["141_aaaa_Lvoov"]("L,a,i,j,b") += reused_["13_aa_ov"]("k,a") * tmps_["138_aaaa_Loovo"]("L,i,j,b,k");
    tmps_["141_aaaa_Lvoov"]("L,a,i,j,b") -= eri["aaaa_oovo"]("i,j,a,k") * l1["aa"]("L,k,b");
    tmps_["141_aaaa_Lvoov"]("L,a,i,j,b") -= l1_1["aa"]("L,k,b") * reused_["131_aaaa_oovo"]("i,j,a,k");
    tmps_["141_aaaa_Lvoov"]("L,a,i,j,b") += reused_["132_aa_vv"]("d,a") * l2_1["aaaa"]("L,i,j,d,b");
    tmps_["141_aaaa_Lvoov"]("L,a,i,j,b") -= dp["aa_vv"]("d,a") * l2_1["aaaa"]("L,i,j,d,b");
    tmps_["141_aaaa_Lvoov"]("L,a,i,j,b") -= reused_["59_aa_vv"]("d,a") * l2["aaaa"]("L,i,j,d,b");
    tmps_["141_aaaa_Lvoov"]("L,a,i,j,b") += f["aa_vv"]("d,a") * l2["aaaa"]("L,i,j,d,b");
    tmps_["141_aaaa_Lvoov"]("L,a,i,j,b") += reused_["123_aa_vv"]("d,a") * l2_1["aaaa"]("L,i,j,d,b");
    tmps_["141_aaaa_Lvoov"]("L,a,i,j,b") -= reused_["57_aa_vv"]("a,d") * l2["aaaa"]("L,i,j,d,b");
    tmps_["141_aaaa_Lvoov"]("L,a,i,j,b") -= 0.50 * reused_["18_aa_vv"]("a,d") * l2["aaaa"]("L,i,j,d,b");
    tmps_["137_aaaa_Loovo"].~TArrayD();

    // sigmal2_aaaa += +1.00 P(a,b) d-_aa(k,a) t1_1_aa(c,k) l2_aaaa(i,j,c,b)
    //              += +1.00 P(a,b) d-_aa(k,a) t0_1 t1_1_aa(c,k) l2_1_aaaa(i,j,c,b)
    //              += -0.50 P(a,b) <i,j||d,a>_aaaa t2_aaaa(d,c,k,l) l2_aaaa(k,l,c,b)
    //              += -0.50 P(a,b) <i,j||d,a>_aaaa t2_1_aaaa(d,c,k,l) l2_1_aaaa(k,l,c,b)
    //              += +0.50 P(a,b) <i,j||d,a>_aaaa t2_abab(d,c,k,l) l2_abab(k,l,b,c)
    //              += +0.50 P(a,b) <i,j||d,a>_aaaa t2_abab(d,c,l,k) l2_abab(l,k,b,c)
    //              += +0.50 P(a,b) <i,j||d,a>_aaaa t2_1_abab(d,c,k,l) l2_1_abab(k,l,b,c)
    //              += +0.50 P(a,b) <i,j||d,a>_aaaa t2_1_abab(d,c,l,k) l2_1_abab(l,k,b,c)
    //              += -1.00 P(a,b) f_aa(k,a) t1_1_aa(c,k) l2_1_aaaa(i,j,c,b)
    //              += -1.00 P(a,b) <i,j||a,k>_aaaa l1_aa(k,b)
    //              += +1.00 P(a,b) <i,j||c,a>_aaaa t1_1_aa(c,k) l1_1_aa(k,b)
    //              += +1.00 P(a,b) <k,c||d,a>_aaaa t1_1_aa(d,k) l2_1_aaaa(i,j,c,b)
    //              += -0.50 P(a,b) <l,k||d,a>_aaaa t2_1_aaaa(d,c,l,k) l2_1_aaaa(i,j,c,b)
    //              += -1.00 P(a,b) d+_aa(c,a) l2_1_aaaa(i,j,c,b)
    //              += -1.00 P(a,b) d-_aa(c,a) t0_1 l2_aaaa(i,j,c,b)
    //              += +1.00 P(a,b) f_aa(c,a) l2_aaaa(i,j,c,b)
    //              += +1.00 P(a,b) <c,k||a,d>_abab t1_1_bb(d,k) l2_1_aaaa(i,j,c,b)
    //              += -0.50 P(a,b) <l,k||a,d>_abab t2_1_abab(c,d,l,k) l2_1_aaaa(i,j,c,b)
    //              += -0.50 P(a,b) <k,l||a,d>_abab t2_1_abab(c,d,k,l) l2_1_aaaa(i,j,c,b)
    //              += -0.50 P(a,b) <l,k||a,d>_abab t2_abab(c,d,l,k) l2_aaaa(i,j,c,b)
    //              += -0.50 P(a,b) <k,l||a,d>_abab t2_abab(c,d,k,l) l2_aaaa(i,j,c,b)
    //              += -0.50 P(a,b) <l,k||d,a>_aaaa t2_aaaa(d,c,l,k) l2_aaaa(i,j,c,b)
    sigmal2_aaaa("L,a,b,i,j") += tmps_["141_aaaa_Lvoov"]("L,a,i,j,b");
    sigmal2_aaaa("L,a,b,i,j") -= tmps_["141_aaaa_Lvoov"]("L,b,i,j,a");
    tmps_["141_aaaa_Lvoov"].~TArrayD();

    // flops: o2v2L1  = o2v2L1 o3v2L1 o3v2L1 o2v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["142_aaaa_Lvoov"]("R,a,i,j,b")  = r0_1("R") * reused_["76_aaaa_vvoo"]("a,b,i,j");
    tmps_["142_aaaa_Lvoov"]("R,a,i,j,b") += r1["aa"]("R,b,k") * reused_["134_aaaa_vooo"]("a,i,j,k");
    tmps_["142_aaaa_Lvoov"]("R,a,i,j,b") += r1_1["aa"]("R,b,k") * reused_["82_aaaa_vooo"]("a,i,j,k");

    // sigmar2_aaaa += -1.00 P(a,b) d-_aa(k,c) r1_aa(b,k) t2_1_aaaa(c,a,i,j)
    //              += -1.00 P(a,b) d-_aa(k,c) r1_1_aa(b,k) t2_aaaa(c,a,i,j)
    //              += -1.00 P(a,b) d-_aa(a,c) r0_1 t2_aaaa(c,b,i,j)
    sigmar2_aaaa("R,a,b,i,j") -= tmps_["142_aaaa_Lvoov"]("R,a,i,j,b");
    sigmar2_aaaa("R,a,b,i,j") += tmps_["142_aaaa_Lvoov"]("R,b,i,j,a");

    // flops: o0v2L1  = o1v3L1 o1v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["147_aa_Lvv"]("R,a,b")  = eri["aaaa_vovv"]("a,j,b,d") * r1_1["aa"]("R,d,j");
    tmps_["147_aa_Lvv"]("R,a,b") += eri["abab_vovv"]("a,i,b,c") * r1_1["bb"]("R,c,i");

    // flops: o0v2L1  = o2v3L1 o2v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["146_aa_Lvv"]("R,a,b")  = r2_1["abab"]("R,a,c,i,j") * eri["abab_oovv"]("i,j,b,c");
    tmps_["146_aa_Lvv"]("R,a,b") += 0.50 * eri["aaaa_oovv"]("k,i,b,d") * r2_1["aaaa"]("R,d,a,i,k");

    // flops: o1v1L1  = o2v2L1 o2v2L1 o1v1L1
    //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1
    tmps_["144_aa_Lov"]("R,i,a")  = eri["aaaa_oovv"]("i,k,a,c") * r1["aa"]("R,c,k");
    tmps_["144_aa_Lov"]("R,i,a") += eri["abab_oovv"]("i,j,a,b") * r1["bb"]("R,b,j");

    // flops: o3v1L1  = o3v2L1 o3v3L1 o3v1L1
    //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1
    tmps_["145_aaaa_Lvooo"]("R,a,i,j,k")  = t2["aaaa"]("b,a,i,j") * tmps_["144_aa_Lov"]("R,k,b");
    tmps_["145_aaaa_Lvooo"]("R,a,i,j,k") -= 0.50 * eri["aaaa_vovv"]("a,k,b,c") * r2["aaaa"]("R,b,c,i,j");

    // flops: o3v1L1  = o3v2L1 o3v2L1 o3v1L1 o3v2L1 o3v1L1
    //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1
    tmps_["143_aaaa_Lvooo"]("R,a,i,j,k")  = -0.50 * r2["aaaa"]("R,b,a,i,j") * f["aa_ov"]("k,b");
    tmps_["143_aaaa_Lvooo"]("R,a,i,j,k") += r2_1["aaaa"]("R,b,a,i,j") * dp["aa_ov"]("k,b");
    tmps_["143_aaaa_Lvooo"]("R,a,i,j,k") += 0.50 * r2["aaaa"]("R,b,a,i,j") * reused_["13_aa_ov"]("k,b");

    // flops: o2v2L1  = o2v3L1 o2v3L1 o2v3L1 o2v3L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v3L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b")  = -2.00 * dp["aa_vv"]("a,c") * r2["aaaa"]("R,c,b,i,j");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") -= 2.00 * reused_["57_aa_vv"]("d,a") * r2_1["aaaa"]("R,d,b,i,j");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") += 2.00 * f["aa_vv"]("a,c") * r2_1["aaaa"]("R,c,b,i,j");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") += r2_1["aaaa"]("R,d,b,i,j") * reused_["23_aa_vv"]("a,d");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") += 2.00 * r1["aa"]("R,b,n") * reused_["133_aaaa_ooov"]("n,i,j,a");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") += 2.00 * r1_1["aa"]("R,b,k") * reused_["130_aaaa_vooo"]("a,i,j,k");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") -= 4.00 * r1_1["aa"]("R,b,k") * reused_["134_aaaa_vooo"]("a,i,j,k");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") -= 2.00 * eri["aaaa_vooo"]("a,k,i,j") * r1_1["aa"]("R,b,k");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") -= 2.00 * r2["aaaa"]("R,d,b,i,j") * reused_["122_aa_vv"]("a,d");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") += 2.00 * r1["aa"]("R,b,n") * reused_["136_aaaa_vooo"]("a,i,j,n");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") -= 4.00 * r0_1("R") * reused_["97_aaaa_vvoo"]("b,a,i,j");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") += 2.00 * r2["aaaa"]("R,d,b,i,j") * reused_["123_aa_vv"]("a,d");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") -= reused_["135_aaaa_vooo"]("a,k,i,j") * r1["aa"]("R,b,k");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") -= 2.00 * reused_["59_aa_vv"]("a,c") * r2_1["aaaa"]("R,c,b,i,j");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") += 2.00 * tmps_["134_aa_Lvv"]("R,a,c") * t2_1["aaaa"]("c,b,i,j");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") -= 2.00 * t0_1 * tmps_["142_aaaa_Lvoov"]("R,a,i,j,b");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") += 2.00 * t2["aaaa"]("c,a,i,j") * tmps_["146_aa_Lvv"]("R,b,c");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") += t2_1["aaaa"]("c,a,i,j") * tmps_["135_aa_Lvv"]("R,b,c");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") += 4.00 * t1_1["aa"]("a,k") * tmps_["143_aaaa_Lvooo"]("R,b,i,j,k");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") += 2.00 * tmps_["147_aa_Lvv"]("R,a,c") * t2["aaaa"]("c,b,i,j");
    tmps_["148_aaaa_Lvoov"]("R,a,i,j,b") += 2.00 * tmps_["145_aaaa_Lvooo"]("R,a,i,j,k") * t1_1["aa"]("b,k");
    tmps_["145_aaaa_Lvooo"].~TArrayD();
    tmps_["143_aaaa_Lvooo"].~TArrayD();
    tmps_["142_aaaa_Lvoov"].~TArrayD();

    // sigmar2_1_aaaa += -0.50 P(a,b) <l,k||c,d>_aaaa r2_aaaa(d,b,l,k) t2_1_aaaa(c,a,i,j)
    //                += +0.50 P(a,b) <l,k||c,d>_abab r2_abab(b,d,l,k) t2_1_aaaa(c,a,i,j)
    //                += +0.50 P(a,b) <k,l||c,d>_abab r2_abab(b,d,k,l) t2_1_aaaa(c,a,i,j)
    //                += +0.50 P(a,b) <l,k||c,d>_abab r2_1_abab(b,d,l,k) t2_aaaa(c,a,i,j)
    //                += +0.50 P(a,b) <k,l||c,d>_abab r2_1_abab(b,d,k,l) t2_aaaa(c,a,i,j)
    //                += -0.50 P(a,b) <l,k||c,d>_aaaa r2_1_aaaa(d,b,l,k) t2_aaaa(c,a,i,j)
    //                += -1.00 P(a,b) d-_aa(k,c) r1_aa(b,k) t0_1 t2_1_aaaa(c,a,i,j)
    //                += -1.00 P(a,b) d-_aa(k,c) r1_1_aa(b,k) t2_aaaa(c,a,i,j) t0_1
    //                += -1.00 P(a,b) d-_aa(a,c) r0_1 t2_aaaa(c,b,i,j) t0_1
    //                += +2.00 P(a,b) d-_aa(k,c) r2_1_aaaa(c,b,i,j) t1_1_aa(a,k)
    //                += -1.00 P(a,b) f_aa(k,c) r2_aaaa(c,b,i,j) t1_1_aa(a,k)
    //                += +1.00 P(a,b) d-_aa(k,c) r2_aaaa(c,b,i,j) t0_1 t1_1_aa(a,k)
    //                += +1.00 P(a,b) <a,k||c,d>_abab r1_1_bb(d,k) t2_aaaa(c,b,i,j)
    //                += -1.00 P(a,b) <k,a||c,d>_aaaa r1_1_aa(d,k) t2_aaaa(c,b,i,j)
    //                += -1.00 P(a,b) <k,a||c,d>_aaaa r1_aa(d,k) t2_1_aaaa(c,b,i,j)
    //                += +1.00 P(a,b) <a,k||c,d>_abab r1_bb(d,k) t2_1_aaaa(c,b,i,j)
    //                += -0.50 P(a,b) <l,k||c,d>_aaaa r2_1_aaaa(d,b,i,j) t2_aaaa(c,a,l,k)
    //                += +1.00 P(a,b) f_aa(a,c) r2_1_aaaa(c,b,i,j)
    //                += -0.50 P(a,b) <l,k||d,c>_abab r2_1_aaaa(d,b,i,j) t2_abab(a,c,l,k)
    //                += -0.50 P(a,b) <k,l||d,c>_abab r2_1_aaaa(d,b,i,j) t2_abab(a,c,k,l)
    //                += -1.00 P(a,b) <l,k||i,j>_aaaa r1_aa(b,l) t1_1_aa(a,k)
    //                += +1.00 P(a,b) <l,k||c,d>_aaaa r1_aa(b,l) t2_aaaa(c,a,i,j) t1_1_aa(d,k)
    //                += -0.50 P(a,b) <l,k||c,d>_aaaa r1_aa(b,l) t2_aaaa(c,d,i,j) t1_1_aa(a,k)
    //                += +1.00 P(a,b) f_aa(k,c) r1_1_aa(b,k) t2_aaaa(c,a,i,j)
    //                += +0.50 P(a,b) <k,a||c,d>_aaaa r1_1_aa(b,k) t2_aaaa(c,d,i,j)
    //                += -1.00 P(a,b) d+_aa(a,c) r2_aaaa(c,b,i,j)
    //                += -2.00 P(a,b) d-_aa(k,c) r1_1_aa(b,k) t2_1_aaaa(c,a,i,j)
    //                += +1.00 P(a,b) <k,a||i,j>_aaaa r1_1_aa(b,k)
    //                += +1.00 P(a,b) <k,a||c,d>_aaaa r2_aaaa(d,b,i,j) t1_1_aa(c,k)
    //                += -0.50 P(a,b) <l,k||c,d>_aaaa r2_aaaa(d,b,i,j) t2_1_aaaa(c,a,l,k)
    //                += +1.00 P(a,b) <l,k||c,d>_abab r1_aa(b,l) t2_aaaa(c,a,i,j) t1_1_bb(d,k)
    //                += -2.00 P(a,b) d-_aa(k,c) r0_1 t2_aaaa(c,a,i,j) t1_1_aa(b,k)
    //                += -2.00 P(a,b) d-_aa(a,c) r0_1 t2_1_aaaa(c,b,i,j)
    //                += +1.00 P(a,b) f_aa(a,c) r0_1 t2_aaaa(c,b,i,j)
    //                += +1.00 P(a,b) <a,k||d,c>_abab r2_aaaa(d,b,i,j) t1_1_bb(c,k)
    //                += -0.50 P(a,b) <l,k||d,c>_abab r2_aaaa(d,b,i,j) t2_1_abab(a,c,l,k)
    //                += -0.50 P(a,b) <k,l||d,c>_abab r2_aaaa(d,b,i,j) t2_1_abab(a,c,k,l)
    //                += +0.50 P(a,b) <k,a||c,d>_aaaa r1_aa(b,k) t2_1_aaaa(c,d,i,j)
    //                += +1.00 P(a,b) f_aa(k,c) r1_aa(b,k) t2_1_aaaa(c,a,i,j)
    //                += -1.00 P(a,b) d-_aa(a,c) r2_1_aaaa(c,b,i,j) t0_1
    //                += +1.00 P(a,b) <k,l||c,d>_abab r1_bb(d,l) t2_aaaa(c,a,i,j) t1_1_aa(b,k)
    //                += -1.00 P(a,b) <l,k||c,d>_aaaa r1_aa(d,l) t2_aaaa(c,a,i,j) t1_1_aa(b,k)
    //                += +0.50 P(a,b) <k,a||c,d>_aaaa r2_aaaa(c,d,i,j) t1_1_aa(b,k)
    sigmar2_1_aaaa("R,a,b,i,j") += 0.50 * tmps_["148_aaaa_Lvoov"]("R,a,i,j,b");
    sigmar2_1_aaaa("R,a,b,i,j") -= 0.50 * tmps_["148_aaaa_Lvoov"]("R,b,i,j,a");
    tmps_["148_aaaa_Lvoov"].~TArrayD();

    // flops: o2v2L1  = o3v2L1 o3v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1
    tmps_["149_bbbb_Lvoov"]("R,a,i,j,b")  = reused_["51_bbbb_vooo"]("a,i,j,k") * r1_1["bb"]("R,b,k");
    tmps_["149_bbbb_Lvoov"]("R,a,i,j,b") += reused_["138_bbbb_vooo"]("a,i,j,k") * r1["bb"]("R,b,k");

    // sigmar2_bbbb += -1.00 P(a,b) d-_bb(k,c) r1_1_bb(b,k) t2_bbbb(c,a,i,j)
    //              += -1.00 P(a,b) d-_bb(k,c) r1_bb(b,k) t2_1_bbbb(c,a,i,j)
    sigmar2_bbbb("R,a,b,i,j") -= tmps_["149_bbbb_Lvoov"]("R,a,i,j,b");
    sigmar2_bbbb("R,a,b,i,j") += tmps_["149_bbbb_Lvoov"]("R,b,i,j,a");

    // flops: o0v2L1  = o2v3L1 o2v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["154_bb_Lvv"]("R,a,b")  = r2_1["abab"]("R,c,a,i,j") * eri["abab_oovv"]("i,j,c,b");
    tmps_["154_bb_Lvv"]("R,a,b") += 0.50 * eri["bbbb_oovv"]("j,k,b,d") * r2_1["bbbb"]("R,d,a,k,j");

    // flops: o0v2L1  = o1v3L1 o1v3L1 o0v2L1
    //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
    tmps_["153_bb_Lvv"]("R,a,b")  = -1.00 * eri["baab_vovv"]("a,j,d,b") * r1_1["aa"]("R,d,j");
    tmps_["153_bb_Lvv"]("R,a,b") += eri["bbbb_vovv"]("a,i,b,c") * r1_1["bb"]("R,c,i");

    // flops: o1v1L1  = o2v2L1
    //  mems: o1v1L1  = o1v1L1
    tmps_["151_bb_Lov"]("R,i,a")  = eri["bbbb_oovv"]("i,j,a,b") * r1["bb"]("R,b,j");

    // flops: o3v1L1  = o1v1L1 o3v2L1 o3v3L1 o3v1L1
    //  mems: o3v1L1  = o1v1L1 o3v1L1 o3v1L1 o3v1L1
    tmps_["152_bbbb_Lovoo"]("R,i,a,j,k")  = (tmps_["68_bb_Lov"]("R,i,b") + tmps_["151_bb_Lov"]("R,i,b")) * t2["bbbb"]("b,a,j,k");
    tmps_["152_bbbb_Lovoo"]("R,i,a,j,k") -= 0.50 * eri["bbbb_vovv"]("a,i,b,c") * r2["bbbb"]("R,b,c,j,k");

    // flops: o3v1L1  = o1v1 o3v2L1 o3v2L1 o3v1L1
    //  mems: o3v1L1  = o1v1 o3v1L1 o3v1L1 o3v1L1
    tmps_["150_bbbb_Lovoo"]("R,i,a,j,k")  = (reused_["12_bb_ov"]("i,b") + -1.00 * f["bb_ov"]("i,b")) * r2["bbbb"]("R,b,a,j,k");
    tmps_["150_bbbb_Lovoo"]("R,i,a,j,k") += 2.00 * dp["bb_ov"]("i,b") * r2_1["bbbb"]("R,b,a,j,k");

    // flops: o2v2L1  = o2v3L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o0v0L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v3L1 o2v2L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1 o2v3L1 o2v2L1 o2v2L1 o3v2L1 o2v2L1
    //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o0v0L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b")  = tmps_["153_bb_Lvv"]("R,a,c") * t2["bbbb"]("c,b,i,j");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") -= t0_1 * tmps_["149_bbbb_Lvoov"]("R,a,i,j,b");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") -= 2.00 * r0_1("R") * reused_["100_bbbb_voov"]("a,i,j,b");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") -= 2.00 * reused_["138_bbbb_vooo"]("a,i,j,k") * r1_1["bb"]("R,b,k");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") -= dp["bb_vv"]("a,c") * r2["bbbb"]("R,c,b,i,j");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") += reused_["53_bbbb_vooo"]("a,i,j,k") * r1_1["bb"]("R,b,k");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") += f["bb_vv"]("a,c") * r2_1["bbbb"]("R,c,b,i,j");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") -= reused_["34_bb_vv"]("a,c") * r2_1["bbbb"]("R,c,b,i,j");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") += 0.50 * reused_["107_bb_vv"]("a,d") * r2["bbbb"]("R,d,b,i,j");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") -= t0_1 * r0_1("R") * reused_["52_bbbb_vvoo"]("a,b,i,j");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") -= eri["bbbb_vooo"]("a,k,i,j") * r1_1["bb"]("R,b,k");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") += reused_["140_bbbb_vooo"]("a,i,j,m") * r1["bb"]("R,b,m");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") -= reused_["108_bb_vv"]("a,d") * r2["bbbb"]("R,d,b,i,j");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") -= reused_["3_bb_vv"]("d,a") * r2_1["bbbb"]("R,d,b,i,j");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") -= reused_["137_bbbb_ovoo"]("m,a,i,j") * r1["bb"]("R,b,m");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") += 0.50 * reused_["29_bb_vv"]("a,d") * r2_1["bbbb"]("R,d,b,i,j");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") += reused_["139_bbbb_vooo"]("a,i,j,k") * r1["bb"]("R,b,k");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") += 0.50 * t2_1["bbbb"]("c,a,i,j") * tmps_["64_bb_Lvv"]("R,b,c");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") += t2["bbbb"]("c,a,i,j") * tmps_["154_bb_Lvv"]("R,b,c");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") += t1_1["bb"]("a,k") * tmps_["150_bbbb_Lovoo"]("R,k,b,i,j");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") -= tmps_["65_bb_Lvv"]("R,a,c") * t2_1["bbbb"]("c,b,i,j");
    tmps_["155_bbbb_Lvoov"]("R,a,i,j,b") += tmps_["152_bbbb_Lovoo"]("R,k,a,i,j") * t1_1["bb"]("b,k");
    tmps_["152_bbbb_Lovoo"].~TArrayD();
    tmps_["150_bbbb_Lovoo"].~TArrayD();
    tmps_["149_bbbb_Lvoov"].~TArrayD();

    // sigmar2_1_bbbb += +1.00 P(a,b) <l,k||d,c>_abab r1_aa(d,l) t2_bbbb(c,a,i,j) t1_1_bb(b,k)
    //                += -1.00 P(a,b) <l,k||c,d>_bbbb r1_bb(d,l) t2_bbbb(c,a,i,j) t1_1_bb(b,k)
    //                += +0.50 P(a,b) <k,a||c,d>_bbbb r2_bbbb(c,d,i,j) t1_1_bb(b,k)
    //                += +0.50 P(a,b) <l,k||d,c>_abab r2_1_abab(d,b,l,k) t2_bbbb(c,a,i,j)
    //                += +0.50 P(a,b) <k,l||d,c>_abab r2_1_abab(d,b,k,l) t2_bbbb(c,a,i,j)
    //                += -0.50 P(a,b) <l,k||c,d>_bbbb r2_1_bbbb(d,b,l,k) t2_bbbb(c,a,i,j)
    //                += -0.50 P(a,b) <l,k||c,d>_bbbb r2_bbbb(d,b,l,k) t2_1_bbbb(c,a,i,j)
    //                += +0.50 P(a,b) <l,k||d,c>_abab r2_abab(d,b,l,k) t2_1_bbbb(c,a,i,j)
    //                += +0.50 P(a,b) <k,l||d,c>_abab r2_abab(d,b,k,l) t2_1_bbbb(c,a,i,j)
    //                += -2.00 P(a,b) d-_bb(k,c) r0_1 t2_bbbb(c,a,i,j) t1_1_bb(b,k)
    //                += -2.00 P(a,b) d-_bb(a,c) r0_1 t2_1_bbbb(c,b,i,j)
    //                += +1.00 P(a,b) f_bb(a,c) r0_1 t2_bbbb(c,b,i,j)
    //                += -2.00 P(a,b) d-_bb(k,c) r1_1_bb(b,k) t2_1_bbbb(c,a,i,j)
    //                += -1.00 P(a,b) d+_bb(a,c) r2_bbbb(c,b,i,j)
    //                += +1.00 P(a,b) f_bb(k,c) r1_1_bb(b,k) t2_bbbb(c,a,i,j)
    //                += +0.50 P(a,b) <k,a||c,d>_bbbb r1_1_bb(b,k) t2_bbbb(c,d,i,j)
    //                += +1.00 P(a,b) f_bb(a,c) r2_1_bbbb(c,b,i,j)
    //                += -1.00 P(a,b) d-_bb(a,c) r2_1_bbbb(c,b,i,j) t0_1
    //                += -0.50 P(a,b) <l,k||c,d>_bbbb r2_bbbb(d,b,i,j) t2_1_bbbb(c,a,l,k)
    //                += +1.00 P(a,b) <k,a||c,d>_bbbb r2_bbbb(d,b,i,j) t1_1_bb(c,k)
    //                += -1.00 P(a,b) d-_bb(a,c) r0_1 t2_bbbb(c,b,i,j) t0_1
    //                += +1.00 P(a,b) <k,a||i,j>_bbbb r1_1_bb(b,k)
    //                += +1.00 P(a,b) <k,l||d,c>_abab r1_bb(b,l) t2_bbbb(c,a,i,j) t1_1_aa(d,k)
    //                += +1.00 P(a,b) <k,a||c,d>_abab r2_bbbb(d,b,i,j) t1_1_aa(c,k)
    //                += -0.50 P(a,b) <l,k||c,d>_abab r2_bbbb(d,b,i,j) t2_1_abab(c,a,l,k)
    //                += -0.50 P(a,b) <k,l||c,d>_abab r2_bbbb(d,b,i,j) t2_1_abab(c,a,k,l)
    //                += -0.50 P(a,b) <l,k||c,d>_abab r2_1_bbbb(d,b,i,j) t2_abab(c,a,l,k)
    //                += -0.50 P(a,b) <k,l||c,d>_abab r2_1_bbbb(d,b,i,j) t2_abab(c,a,k,l)
    //                += +1.00 P(a,b) <l,k||c,d>_bbbb r1_bb(b,l) t2_bbbb(c,a,i,j) t1_1_bb(d,k)
    //                += -1.00 P(a,b) <l,k||i,j>_bbbb r1_bb(b,l) t1_1_bb(a,k)
    //                += -0.50 P(a,b) <l,k||c,d>_bbbb r1_bb(b,l) t2_bbbb(c,d,i,j) t1_1_bb(a,k)
    //                += -0.50 P(a,b) <l,k||c,d>_bbbb r2_1_bbbb(d,b,i,j) t2_bbbb(c,a,l,k)
    //                += +1.00 P(a,b) f_bb(k,c) r1_bb(b,k) t2_1_bbbb(c,a,i,j)
    //                += +0.50 P(a,b) <k,a||c,d>_bbbb r1_bb(b,k) t2_1_bbbb(c,d,i,j)
    //                += -1.00 P(a,b) d-_bb(k,c) r1_1_bb(b,k) t2_bbbb(c,a,i,j) t0_1
    //                += -1.00 P(a,b) d-_bb(k,c) r1_bb(b,k) t0_1 t2_1_bbbb(c,a,i,j)
    //                += +1.00 P(a,b) d-_bb(k,c) r2_bbbb(c,b,i,j) t0_1 t1_1_bb(a,k)
    //                += -1.00 P(a,b) f_bb(k,c) r2_bbbb(c,b,i,j) t1_1_bb(a,k)
    //                += +2.00 P(a,b) d-_bb(k,c) r2_1_bbbb(c,b,i,j) t1_1_bb(a,k)
    //                += +1.00 P(a,b) <k,a||d,c>_abab r1_aa(d,k) t2_1_bbbb(c,b,i,j)
    //                += -1.00 P(a,b) <k,a||c,d>_bbbb r1_bb(d,k) t2_1_bbbb(c,b,i,j)
    //                += -1.00 P(a,b) <k,a||c,d>_bbbb r1_1_bb(d,k) t2_bbbb(c,b,i,j)
    //                += +1.00 P(a,b) <k,a||d,c>_abab r1_1_aa(d,k) t2_bbbb(c,b,i,j)
    sigmar2_1_bbbb("R,a,b,i,j") += tmps_["155_bbbb_Lvoov"]("R,a,i,j,b");
    sigmar2_1_bbbb("R,a,b,i,j") -= tmps_["155_bbbb_Lvoov"]("R,b,i,j,a");
    tmps_["155_bbbb_Lvoov"].~TArrayD();
}
