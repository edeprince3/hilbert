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

#include "cc_cavity/include/qed_ccsd_21/eom_ea_qed_ccsd_21.h"
#include "cc_cavity/include/qed_ccsd_21/qed_ccsd_21.h"

void hilbert::EOM_EA_QED_CCSD_21::sigma_ea_21_3() {

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

    /// reference right operators
    TArrayMap r1 = {
            {"a", evec_blks_["r1_a"]},
            {"b", evec_blks_["r1_b"]}
    };
    TArrayMap r2 = {
            {"aaa", evec_blks_["r2_aaa"]},
            {"abb", evec_blks_["r2_abb"]},
    };
    TArrayMap r1_1 = {
            {"a", evec_blks_["r1_1_a"]},
            {"b", evec_blks_["r1_1_b"]}
    };
    TArrayMap r2_1 = {
            {"aaa", evec_blks_["r2_1_aaa"]},
            {"abb", evec_blks_["r2_1_abb"]},
    };

    /// reference left operators

    TArrayMap l1 = {
            {"a", evec_blks_["l1_a"]},
            {"b", evec_blks_["l1_b"]}
    };
    TArrayMap l2 = {
            {"aaa", evec_blks_["l2_aaa"]},
            {"bab", evec_blks_["l2_bab"]},
    };
    TArrayMap l1_1 = {
            {"a", evec_blks_["l1_1_a"]},
            {"b", evec_blks_["l1_1_b"]}
    };
    TArrayMap l2_1 = {
            {"aaa", evec_blks_["l2_1_aaa"]},
            {"bab", evec_blks_["l2_1_bab"]},
    };

    /// reference sigma vectors

    // electronic

    auto &sigmar1_a = sigvec_blks_["sigmar1_a"];
    auto &sigmal1_a = sigvec_blks_["sigmal1_a"];

    auto &sigmar2_aaa = sigvec_blks_["sigmar2_aaa"];
    auto &sigmar2_abb = sigvec_blks_["sigmar2_abb"];
    auto &sigmal2_aaa = sigvec_blks_["sigmal2_aaa"];
    auto &sigmal2_abb = sigvec_blks_["sigmal2_abb"];

    // photonic

    auto &sigmar1_1_a = sigvec_blks_["sigmar1_1_a"];
    auto &sigmal1_1_a = sigvec_blks_["sigmal1_1_a"];

    auto &sigmar2_1_aaa = sigvec_blks_["sigmar2_1_aaa"];
    auto &sigmar2_1_abb = sigvec_blks_["sigmar2_1_abb"];
    auto &sigmal2_1_aaa = sigvec_blks_["sigmal2_1_aaa"];
    auto &sigmal2_1_abb = sigvec_blks_["sigmal2_1_abb"];


    /// extract coherent state basis operators

    // electronic

    auto &csigmar1_a = tmps_["csigmar1_a"];
    auto &csigmal1_a = tmps_["csigmal1_a"];

    auto &csigmar2_aaa = tmps_["csigmar2_aaa"];
    auto &csigmar2_abb = tmps_["csigmar2_abb"];
    auto &csigmal2_aaa = tmps_["csigmal2_aaa"];
    auto &csigmal2_abb = tmps_["csigmal2_abb"];

    // photonic

    auto &csigmar1_1_a = tmps_["csigmar1_1_a"];
    auto &csigmal1_1_a = tmps_["csigmal1_1_a"];

    auto &csigmar2_1_aaa = tmps_["csigmar2_1_aaa"];
    auto &csigmar2_1_abb = tmps_["csigmar2_1_abb"];
    auto &csigmal2_1_aaa = tmps_["csigmal2_1_aaa"];
    auto &csigmal2_1_abb = tmps_["csigmal2_1_abb"];


    // flops: o1v2L1  = o1v2L1 o1v2L1 o1v3L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o2v3L1 o1v3L1 o1v3L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["75_aaa_Lvov"]("R,a,i,b")  = -1.00 * r1["a"]("R,a") * reused_["34_aa_vo"]("b,i");
    tmps_["75_aaa_Lvov"]("R,a,i,b") -= r1["a"]("R,a") * reused_["36_aa_vo"]("b,i");
    tmps_["75_aaa_Lvov"]("R,a,i,b") -= 0.50 * r2["aaa"]("R,e,a,i") * reused_["1_aa_vv"]("b,e");
    tmps_["75_aaa_Lvov"]("R,a,i,b") -= r1["a"]("R,e") * reused_["26_aaaa_vvvo"]("b,e,a,i");
    tmps_["75_aaa_Lvov"]("R,a,i,b") -= r1_1["a"]("R,a") * dp["aa_vo"]("b,i");
    tmps_["75_aaa_Lvov"]("R,a,i,b") -= r1["a"]("R,e") * reused_["64_aaaa_vvvo"]("b,e,a,i");
    tmps_["75_aaa_Lvov"]("R,a,i,b") += r1["a"]("R,a") * f["aa_vo"]("b,i");
    tmps_["75_aaa_Lvov"]("R,a,i,b") += r2["aaa"]("R,c,a,l") * eri["aaaa_vovo"]("b,l,c,i");
    tmps_["75_aaa_Lvov"]("R,a,i,b") += r2["aaa"]("R,c,a,i") * reused_["3_aa_vv"]("b,c");
    tmps_["75_aaa_Lvov"]("R,a,i,b") += dp["aa_vv"]("b,c") * r2_1["aaa"]("R,c,a,i");
    tmps_["75_aaa_Lvov"]("R,a,i,b") += r2["aaa"]("R,e,a,k") * reused_["44_aaaa_voov"]("b,i,k,e");
    tmps_["75_aaa_Lvov"]("R,a,i,b") += r2["aaa"]("R,e,a,i") * reused_["2_aa_vv"]("b,e");
    tmps_["75_aaa_Lvov"]("R,a,i,b") -= r2["aaa"]("R,c,a,i") * f["aa_vv"]("b,c");
    tmps_["75_aaa_Lvov"]("R,a,i,b") -= r2["abb"]("R,a,d,j") * eri["abba_vovo"]("b,j,d,i");
    tmps_["75_aaa_Lvov"]("R,a,i,b") -= r2["abb"]("R,a,f,m") * reused_["47_aabb_voov"]("b,i,m,f");
    tmps_["75_aaa_Lvov"]("R,a,i,b") += r2["abb"]("R,a,f,m") * reused_["52_aabb_voov"]("b,i,m,f");
    tmps_["75_aaa_Lvov"]("R,a,i,b") += tmps_["71_aaa_Lvvo"]("R,a,b,i");
    tmps_["75_aaa_Lvov"]("R,a,i,b") -= t1_1["aa"]("b,l") * tmps_["70_aaa_Lvoo"]("R,a,i,l");
    tmps_["75_aaa_Lvov"]("R,a,i,b") += t1_1["aa"]("b,i") * tmps_["74_a_Lv"]("R,a");
    tmps_["75_aaa_Lvov"]("R,a,i,b") -= t2["abab"]("b,d,i,j") * tmps_["72_abb_Lvov"]("R,a,j,d");
    tmps_["75_aaa_Lvov"]("R,a,i,b") += t1_1["aa"]("a,i") * tmps_["73_a_Lv"]("R,b");
    tmps_["71_aaa_Lvvo"].~TArrayD();
    tmps_["70_aaa_Lvoo"].~TArrayD();

    // sigmar2_aaa += -1.00 P(a,b) d-_aa(a,c) r2_aaa(c,b,i) t0_1
    //             += -1.00 P(a,b) d-_aa(a,c) r2_1_aaa(c,b,i)
    //             += +1.00 P(a,b) <j,a||c,i>_aaaa r2_aaa(c,b,j)
    //             += +1.00 P(a,b) <k,j||c,d>_aaaa r2_aaa(d,b,k) t2_aaaa(c,a,i,j)
    //             += -1.00 P(a,b) f_aa(a,i) r1_a(b)
    //             += -1.00 P(a,b) <j,a||c,d>_aaaa r1_a(d) t2_aaaa(c,b,i,j)
    //             += +1.00 P(a,b) d-_aa(a,i) r1_1_a(b)
    //             += -0.50 P(a,b) <k,j||d,c>_abab r2_aaa(d,b,i) t2_abab(a,c,k,j)
    //             += -0.50 P(a,b) <j,k||d,c>_abab r2_aaa(d,b,i) t2_abab(a,c,j,k)
    //             += +1.00 P(a,b) <a,j||d,c>_abab r1_a(d) t2_abab(b,c,i,j)
    //             += -0.50 P(a,b) <k,j||c,d>_aaaa r2_aaa(d,b,i) t2_aaaa(c,a,k,j)
    //             += +1.00 P(a,b) d-_aa(a,i) r1_a(b) t0_1
    //             += +0.50 P(a,b) <k,j||c,i>_aaaa r1_a(b) t2_aaaa(c,a,k,j)
    //             += +0.50 P(a,b) <k,j||i,c>_abab r1_a(b) t2_abab(a,c,k,j)
    //             += +0.50 P(a,b) <j,k||i,c>_abab r1_a(b) t2_abab(a,c,j,k)
    //             += +1.00 P(a,b) f_aa(j,c) r1_a(b) t2_aaaa(c,a,i,j)
    //             += -1.00 P(a,b) f_bb(j,c) r1_a(b) t2_abab(a,c,i,j)
    //             += -0.50 P(a,b) <a,j||c,d>_abab r1_a(b) t2_abab(c,d,i,j)
    //             += -0.50 P(a,b) <a,j||d,c>_abab r1_a(b) t2_abab(d,c,i,j)
    //             += +1.00 P(a,b) d-_bb(j,j) r1_a(b) t1_1_aa(a,i)
    //             += +1.00 P(a,b) d-_aa(j,j) r1_a(b) t1_1_aa(a,i)
    //             += +0.50 P(a,b) <j,a||c,d>_aaaa r1_a(b) t2_aaaa(c,d,i,j)
    //             += +1.00 P(a,b) d-_bb(j,c) r1_a(b) t2_abab(a,c,i,j) t0_1
    //             += -1.00 P(a,b) d-_aa(j,c) r1_a(b) t2_aaaa(c,a,i,j) t0_1
    //             += +1.00 P(a,b) f_aa(a,c) r2_aaa(c,b,i)
    //             += +1.00 P(a,b) d-_aa(a,c) r1_a(b) t1_1_aa(c,i)
    //             += +1.00 P(a,b) d-_bb(j,c) r1_a(b) t2_1_abab(a,c,i,j)
    //             += -1.00 P(a,b) d-_aa(j,c) r1_a(b) t2_1_aaaa(c,a,i,j)
    //             += -1.00 P(a,b) d-_aa(j,i) r1_a(b) t1_1_aa(a,j)
    //             += -1.00 P(a,b) <a,j||i,c>_abab r2_abb(b,c,j)
    //             += +1.00 P(a,b) <j,k||c,d>_abab r2_abb(b,d,k) t2_aaaa(c,a,i,j)
    //             += +1.00 P(a,b) <k,j||c,d>_bbbb r2_abb(b,d,k) t2_abab(a,c,i,j)
    //             += +1.00 P(a,b) d-_bb(j,c) r1_1_a(b) t2_abab(a,c,i,j)
    //             += -1.00 P(a,b) d-_aa(j,c) r1_1_a(b) t2_aaaa(c,a,i,j)
    //             += +1.00 P(a,b) d-_aa(j,c) r2_aaa(c,b,i) t1_1_aa(a,j)
    //             += -1.00 P(a,b) d-_aa(j,c) r2_aaa(c,b,j) t1_1_aa(a,i)
    //             += +1.00 P(a,b) d-_bb(j,c) r2_abb(b,c,j) t1_1_aa(a,i)
    //             += +1.00 P(a,b) <k,j||d,c>_abab r2_aaa(d,b,k) t2_abab(a,c,i,j)
    //             += -1.00 P(a,b) d-_aa(a,c) r1_a(c) t1_1_aa(b,i)
    sigmar2_aaa("R,a,b,i") -= tmps_["75_aaa_Lvov"]("R,b,i,a");
    sigmar2_aaa("R,a,b,i") += tmps_["75_aaa_Lvov"]("R,a,i,b");
    tmps_["75_aaa_Lvov"].~TArrayD();

    // flops: o1v2L1  = o1v1L1 o2v2L1
    //  mems: o1v2L1  = o1v0L1 o1v2L1
    tmps_["77_abb_Lvvo"]("R,a,b,i")  = dp["aa_ov"]("j,c") * r1["a"]("R,c") * t2["abab"]("a,b,j,i");

    // sigmar2_1_abb  = +1.00 d+_aa(j,c) r1_a(c) t2_abab(a,b,j,i)
    sigmar2_1_abb("R,a,b,i")  = tmps_["77_abb_Lvvo"]("R,a,b,i");

    // flops: o3v0L1  = o3v2L1
    //  mems: o3v0L1  = o3v0L1
    tmps_["81_abb_Looo"]("R,i,j,k")  = eri["abab_oovv"]("i,j,a,b") * r2["abb"]("R,a,b,k");

    // flops: o1v2L1  = o2v3L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["80_bab_Lvov"]("R,a,i,b")  = r2["abb"]("R,c,a,j") * eri["abab_oovv"]("i,j,c,b");

    // flops: o2v1L1  = o2v2L1
    //  mems: o2v1L1  = o2v1L1
    tmps_["79_abb_Lvoo"]("R,a,i,j")  = r2["abb"]("R,a,b,i") * dp["bb_ov"]("j,b");

    // flops: o2v1L1  = o2v2L1
    //  mems: o2v1L1  = o2v1L1
    tmps_["78_bba_Lvoo"]("R,a,i,j")  = r2["abb"]("R,b,a,i") * dp["aa_ov"]("j,b");

    // flops: o1v0L1  = o1v1L1
    //  mems: o1v0L1  = o1v0L1
    tmps_["76_a_Lo"]("R,i")  = dp["aa_ov"]("i,a") * r1["a"]("R,a");

    // flops: o1v2L1  = o2v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o3v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v3L1 o2v2L1 o1v3L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o2v3L1 o1v2L1 o1v4L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["82_abb_Lvvo"]("R,a,b,i")  = -0.50 * t2["abab"]("a,b,j,i") * tmps_["26_a_Lo"]("R,j");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += t1_1["bb"]("b,l") * tmps_["79_abb_Lvoo"]("R,a,i,l");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += t1_1["bb"]("b,i") * tmps_["74_a_Lv"]("R,a");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += t2_1["abab"]("a,b,j,i") * tmps_["76_a_Lo"]("R,j");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += t2["abab"]("a,b,m,l") * tmps_["81_abb_Looo"]("R,m,l,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += t1_1["aa"]("a,j") * tmps_["78_bba_Lvoo"]("R,b,i,j");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += t2["abab"]("a,b,j,i") * tmps_["24_a_Lo"]("R,j");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += t2["abab"]("a,d,j,i") * tmps_["80_bab_Lvov"]("R,b,j,d");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += t0_1 * tmps_["77_abb_Lvvo"]("R,a,b,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= dp["aa_vv"]("a,c") * r2_1["abb"]("R,c,b,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= r1["a"]("R,f") * reused_["28_aabb_vvvo"]("a,f,b,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= r2["abb"]("R,a,b,l") * f["bb_oo"]("l,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= reused_["3_aa_vv"]("a,c") * r2["abb"]("R,c,b,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= r2["abb"]("R,a,d,i") * reused_["66_bb_vv"]("b,d");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= scalars_["3"] * r2["abb"]("R,a,b,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += r1["a"]("R,c") * reused_["30_abab_vovv"]("c,i,a,b");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += r1_1["a"]("R,a") * reused_["33_bb_vo"]("b,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += r2["aaa"]("R,f,a,m") * reused_["50_bbaa_voov"]("b,i,m,f");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += r2["abb"]("R,a,e,k") * reused_["51_bbbb_voov"]("b,i,k,e");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= scalars_["4"] * r2["abb"]("R,a,b,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += 0.50 * r2["abb"]("R,f,b,i") * reused_["1_aa_vv"]("a,f");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += r2["abb"]("R,a,b,l") * reused_["67_bb_oo"]("l,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= eri["abab_vovo"]("a,l,c,i") * r2["abb"]("R,c,b,l");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += eri["abab_vvvv"]("a,b,c,e") * r2["abb"]("R,c,e,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= scalars_["5"] * r2["abb"]("R,a,b,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += r1["a"]("R,a") * f["bb_vo"]("b,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= scalars_["1"] * r2_1["abb"]("R,a,b,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += r2["aaa"]("R,c,a,j") * eri["baab_vovo"]("b,j,c,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += 0.50 * r2["abb"]("R,a,e,i") * reused_["56_bb_vv"]("b,e");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += f["aa_vv"]("a,c") * r2["abb"]("R,c,b,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += f["bb_vv"]("b,d") * r2["abb"]("R,a,d,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += dp["bb_oo"]("l,i") * r2_1["abb"]("R,a,b,l");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= r2["aaa"]("R,f,a,m") * reused_["48_bbaa_voov"]("b,i,m,f");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= scalars_["2"] * r2_1["abb"]("R,a,b,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= r2["abb"]("R,f,b,i") * reused_["2_aa_vv"]("a,f");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += r1["a"]("R,a") * reused_["37_bb_vo"]("b,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += r1["a"]("R,f") * reused_["29_baab_vvvo"]("b,f,a,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += eri["abab_vvvo"]("a,b,c,i") * r1["a"]("R,c");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= dp["bb_vo"]("b,i") * r1_1["a"]("R,a");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= r2["abb"]("R,a,b,k") * reused_["57_bb_oo"]("k,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= r2["abb"]("R,a,e,i") * reused_["55_bb_vv"]("b,e");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= dp["bb_vv"]("b,d") * r2_1["abb"]("R,a,d,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += r1["a"]("R,a") * reused_["35_bb_ov"]("i,b");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= r1["a"]("R,f") * reused_["65_aabb_vvvo"]("a,f,b,i");
    tmps_["82_abb_Lvvo"]("R,a,b,i") += 0.50 * r2["abb"]("R,a,b,k") * reused_["58_bb_oo"]("i,k");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= eri["bbbb_vovo"]("b,l,d,i") * r2["abb"]("R,a,d,l");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= r2["abb"]("R,a,e,k") * reused_["53_bbbb_voov"]("b,i,k,e");
    tmps_["82_abb_Lvvo"]("R,a,b,i") -= t1_1["bb"]("b,i") * tmps_["73_a_Lv"]("R,a");
    tmps_["79_abb_Lvoo"].~TArrayD();
    tmps_["78_bba_Lvoo"].~TArrayD();
    tmps_["77_abb_Lvvo"].~TArrayD();
    tmps_["74_a_Lv"].~TArrayD();
    tmps_["73_a_Lv"].~TArrayD();

    // sigmar2_abb  = +1.00 d-_aa(j,c) r1_a(c) t2_1_abab(a,b,j,i)
    //             += +1.00 d-_aa(j,c) r2_aaa(c,a,j) t1_1_bb(b,i)
    //             += -1.00 d-_bb(j,c) r2_abb(a,c,j) t1_1_bb(b,i)
    //             += +1.00 d-_bb(j,c) r2_abb(a,c,i) t1_1_bb(b,j)
    //             += +0.50 <k,j||c,d>_aaaa r2_aaa(c,d,k) t2_abab(a,b,j,i)
    //             += -0.50 <j,k||c,d>_abab r2_abb(c,d,k) t2_abab(a,b,j,i)
    //             += -0.50 <j,k||d,c>_abab r2_abb(d,c,k) t2_abab(a,b,j,i)
    //             += -1.00 f_aa(j,c) r1_a(c) t2_abab(a,b,j,i)
    //             += +0.25 <k,j||c,d>_abab r2_abb(c,d,i) t2_abab(a,b,k,j)
    //             += +0.25 <j,k||c,d>_abab r2_abb(c,d,i) t2_abab(a,b,j,k)
    //             += +0.25 <k,j||d,c>_abab r2_abb(d,c,i) t2_abab(a,b,k,j)
    //             += +0.25 <j,k||d,c>_abab r2_abb(d,c,i) t2_abab(a,b,j,k)
    //             += +1.00 d-_aa(j,c) r2_abb(c,b,i) t1_1_aa(a,j)
    //             += +1.00 d-_aa(j,c) r1_1_a(c) t2_abab(a,b,j,i)
    //             += +1.00 <j,k||d,c>_abab r2_abb(d,b,k) t2_abab(a,c,j,i)
    //             += +1.00 d-_aa(j,c) r1_a(c) t2_abab(a,b,j,i) t0_1
    //             += +1.00 <k,j||d,c>_abab r2_aaa(d,a,k) t2_bbbb(c,b,i,j)
    //             += +1.00 d-_bb(j,c) r1_1_a(a) t2_bbbb(c,b,i,j)
    //             += -1.00 d-_aa(j,c) r1_1_a(a) t2_abab(c,b,j,i)
    //             += +0.50 <k,j||c,i>_abab r1_a(c) t2_abab(a,b,k,j)
    //             += +0.50 <j,k||c,i>_abab r1_a(c) t2_abab(a,b,j,k)
    //             += -1.00 d-_bb(j,c) r2_abb(a,b,i) t1_1_bb(c,j)
    //             += +1.00 <j,k||c,d>_abab r2_abb(a,d,k) t2_abab(c,b,j,i)
    //             += -1.00 d-_bb(b,c) r2_abb(a,c,i) t0_1
    //             += -1.00 d-_aa(a,c) r2_abb(c,b,i) t0_1
    //             += -1.00 d-_aa(j,c) r2_abb(a,b,i) t1_1_aa(c,j)
    //             += -0.50 <k,j||c,d>_aaaa r2_abb(d,b,i) t2_aaaa(c,a,k,j)
    //             += +1.00 d-_bb(j,i) r2_abb(a,b,j) t0_1
    //             += -1.00 <a,j||c,i>_abab r2_abb(c,b,j)
    //             += +0.50 <a,b||c,d>_abab r2_abb(c,d,i)
    //             += +0.50 <a,b||d,c>_abab r2_abb(d,c,i)
    //             += -0.50 <k,j||k,j>_abab r2_abb(a,b,i)
    //             += -0.50 <j,k||j,k>_abab r2_abb(a,b,i)
    //             += +0.25 <k,j||c,d>_bbbb r2_abb(a,b,i) t2_bbbb(c,d,k,j)
    //             += +1.00 f_bb(j,j) r2_abb(a,b,i)
    //             += +0.25 <k,j||c,d>_aaaa r2_abb(a,b,i) t2_aaaa(c,d,k,j)
    //             += -1.00 d-_aa(j,j) r2_abb(a,b,i) t0_1
    //             += +1.00 f_aa(j,j) r2_abb(a,b,i)
    //             += -0.50 <k,j||k,j>_bbbb r2_abb(a,b,i)
    //             += -0.50 <k,j||k,j>_aaaa r2_abb(a,b,i)
    //             += +0.25 <k,j||c,d>_abab r2_abb(a,b,i) t2_abab(c,d,k,j)
    //             += +0.25 <j,k||c,d>_abab r2_abb(a,b,i) t2_abab(c,d,j,k)
    //             += +0.25 <k,j||d,c>_abab r2_abb(a,b,i) t2_abab(d,c,k,j)
    //             += +0.25 <j,k||d,c>_abab r2_abb(a,b,i) t2_abab(d,c,j,k)
    //             += -1.00 d-_bb(j,j) r2_abb(a,b,i) t0_1
    //             += +1.00 f_bb(b,i) r1_a(a)
    //             += -1.00 d-_bb(j,j) r2_1_abb(a,b,i)
    //             += -1.00 f_bb(j,i) r2_abb(a,b,j)
    //             += -1.00 <j,b||c,i>_abab r2_aaa(c,a,j)
    //             += -1.00 <a,j||d,c>_abab r1_a(d) t2_bbbb(c,b,i,j)
    //             += -0.50 <k,j||c,d>_bbbb r2_abb(a,d,i) t2_bbbb(c,b,k,j)
    //             += +1.00 f_aa(a,c) r2_abb(c,b,i)
    //             += -1.00 d-_aa(a,c) r2_1_abb(c,b,i)
    //             += +1.00 f_bb(b,c) r2_abb(a,c,i)
    //             += +1.00 d-_bb(j,i) r2_1_abb(a,b,j)
    //             += +1.00 <k,j||c,d>_aaaa r2_aaa(d,a,k) t2_abab(c,b,j,i)
    //             += -1.00 d-_aa(j,j) r2_1_abb(a,b,i)
    //             += -0.50 <k,j||d,c>_abab r2_abb(d,b,i) t2_abab(a,c,k,j)
    //             += -0.50 <j,k||d,c>_abab r2_abb(d,b,i) t2_abab(a,c,j,k)
    //             += +1.00 d-_bb(j,c) r1_a(a) t2_bbbb(c,b,i,j) t0_1
    //             += -1.00 d-_aa(j,c) r1_a(a) t2_abab(c,b,j,i) t0_1
    //             += -1.00 d-_bb(b,i) r1_a(a) t0_1
    //             += -1.00 d-_aa(j,j) r1_a(a) t1_1_bb(b,i)
    //             += -1.00 d-_bb(j,j) r1_a(a) t1_1_bb(b,i)
    //             += -0.50 <k,j||c,i>_abab r1_a(a) t2_abab(c,b,k,j)
    //             += -0.50 <j,k||c,i>_abab r1_a(a) t2_abab(c,b,j,k)
    //             += -0.50 <j,b||c,d>_bbbb r1_a(a) t2_bbbb(c,d,i,j)
    //             += +1.00 f_aa(j,c) r1_a(a) t2_abab(c,b,j,i)
    //             += +0.50 <j,b||c,d>_abab r1_a(a) t2_abab(c,d,j,i)
    //             += +0.50 <j,b||d,c>_abab r1_a(a) t2_abab(d,c,j,i)
    //             += -0.50 <k,j||c,i>_bbbb r1_a(a) t2_bbbb(c,b,k,j)
    //             += -1.00 f_bb(j,c) r1_a(a) t2_bbbb(c,b,i,j)
    //             += -1.00 <j,b||d,c>_abab r1_a(d) t2_abab(a,c,j,i)
    //             += +1.00 <a,b||c,i>_abab r1_a(c)
    //             += -1.00 d-_bb(b,i) r1_1_a(a)
    //             += -0.50 <j,k||c,d>_abab r2_abb(a,b,k) t2_abab(c,d,j,i)
    //             += -0.50 <j,k||d,c>_abab r2_abb(a,b,k) t2_abab(d,c,j,i)
    //             += -0.50 <k,j||c,d>_abab r2_abb(a,d,i) t2_abab(c,b,k,j)
    //             += -0.50 <j,k||c,d>_abab r2_abb(a,d,i) t2_abab(c,b,j,k)
    //             += -1.00 d-_bb(b,c) r2_1_abb(a,c,i)
    //             += +1.00 d-_bb(j,i) r1_a(a) t1_1_bb(b,j)
    //             += -1.00 d-_aa(j,c) r1_a(a) t2_1_abab(c,b,j,i)
    //             += +1.00 d-_bb(j,c) r1_a(a) t2_1_bbbb(c,b,i,j)
    //             += -1.00 d-_bb(b,c) r1_a(a) t1_1_bb(c,i)
    //             += +1.00 <j,a||c,d>_aaaa r1_a(d) t2_abab(c,b,j,i)
    //             += -0.50 <k,j||c,d>_bbbb r2_abb(a,b,k) t2_bbbb(c,d,i,j)
    //             += +1.00 <j,b||c,i>_bbbb r2_abb(a,c,j)
    //             += +1.00 <k,j||c,d>_bbbb r2_abb(a,d,k) t2_bbbb(c,b,i,j)
    //             += -1.00 d-_aa(a,c) r1_a(c) t1_1_bb(b,i)
    sigmar2_abb("R,a,b,i")  = tmps_["82_abb_Lvvo"]("R,a,b,i");
    tmps_["82_abb_Lvvo"].~TArrayD();

    // flops: o0v1L1  = o1v2L1
    //  mems: o0v1L1  = o0v1L1
    tmps_["85_a_Lv"]("L,a")  = l2["aaa"]("L,i,b,a") * t1_1["aa"]("b,i");

    // csigmal1_a += -1.00 t1_1_aa(b,i) l2_aaa(i,b,a)
    csigmal1_a("L,a") -= tmps_["85_a_Lv"]("L,a");

    // flops: o1v2L1  = o2v3L1 o2v3L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1
    tmps_["84_abb_Lvvo"]("L,a,b,i")  = l2["aaa"]("L,j,c,a") * t2["abab"]("c,b,j,i");
    tmps_["84_abb_Lvvo"]("L,a,b,i") += l2_1["aaa"]("L,j,c,a") * t2_1["abab"]("c,b,j,i");

    // flops: o2v1L1  = o2v2L1
    //  mems: o2v1L1  = o2v1L1
    tmps_["83_aaa_Lovo"]("L,i,a,j")  = l2["aaa"]("L,i,b,a") * t1_1["aa"]("b,j");

    // flops: o1v2L1  = o1v2L1 o3v2L1 o2v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o1v3L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o3v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["86_aaa_Lvov"]("L,a,i,b")  = -1.00 * tmps_["85_a_Lv"]("L,b") * dp["aa_ov"]("i,a");
    tmps_["86_aaa_Lvov"]("L,a,i,b") -= tmps_["61_aaa_Lovo"]("L,k,b,j") * eri["aaaa_oovo"]("i,j,a,k");
    tmps_["86_aaa_Lvov"]("L,a,i,b") += dp["aa_ov"]("k,a") * tmps_["83_aaa_Lovo"]("L,i,b,k");
    tmps_["86_aaa_Lvov"]("L,a,i,b") -= f["aa_ov"]("k,a") * tmps_["61_aaa_Lovo"]("L,i,b,k");
    tmps_["86_aaa_Lvov"]("L,a,i,b") += reused_["4_aa_ov"]("k,a") * tmps_["61_aaa_Lovo"]("L,i,b,k");
    tmps_["86_aaa_Lvov"]("L,a,i,b") -= l2["aaa"]("L,i,c,b") * reused_["2_aa_vv"]("c,a");
    tmps_["86_aaa_Lvov"]("L,a,i,b") -= l2["bab"]("L,m,b,d") * reused_["71_bbaa_voov"]("d,m,i,a");
    tmps_["86_aaa_Lvov"]("L,a,i,b") += l2_1["bab"]("L,m,b,d") * reused_["73_baab_vovo"]("d,i,a,m");
    tmps_["86_aaa_Lvov"]("L,a,i,b") += l2_1["bab"]("L,m,b,d") * reused_["72_bbaa_voov"]("d,m,i,a");
    tmps_["86_aaa_Lvov"]("L,a,i,b") += l2["bab"]("L,m,b,d") * reused_["50_bbaa_voov"]("d,m,i,a");
    tmps_["86_aaa_Lvov"]("L,a,i,b") -= l1["a"]("L,b") * f["aa_ov"]("i,a");
    tmps_["86_aaa_Lvov"]("L,a,i,b") -= l2["aaa"]("L,k,c,b") * reused_["70_aaaa_voov"]("c,k,i,a");
    tmps_["86_aaa_Lvov"]("L,a,i,b") += l2_1["aaa"]("L,i,c,b") * reused_["61_aa_vv"]("c,a");
    tmps_["86_aaa_Lvov"]("L,a,i,b") += l1["a"]("L,b") * reused_["4_aa_ov"]("i,a");
    tmps_["86_aaa_Lvov"]("L,a,i,b") -= l2["aaa"]("L,i,c,b") * reused_["3_aa_vv"]("c,a");
    tmps_["86_aaa_Lvov"]("L,a,i,b") -= l1_1["a"]("L,b") * reused_["11_aa_ov"]("i,a");
    tmps_["86_aaa_Lvov"]("L,a,i,b") += l1_1["a"]("L,b") * dp["aa_ov"]("i,a");
    tmps_["86_aaa_Lvov"]("L,a,i,b") -= 0.50 * l2["aaa"]("L,i,c,b") * reused_["31_aa_vv"]("a,c");
    tmps_["86_aaa_Lvov"]("L,a,i,b") -= dp["aa_vv"]("c,a") * l2_1["aaa"]("L,i,c,b");
    tmps_["86_aaa_Lvov"]("L,a,i,b") -= l2_1["aaa"]("L,i,c,b") * reused_["13_aa_vv"]("c,a");
    tmps_["86_aaa_Lvov"]("L,a,i,b") += f["aa_vv"]("c,a") * l2["aaa"]("L,i,c,b");
    tmps_["86_aaa_Lvov"]("L,a,i,b") -= l1_1["a"]("L,b") * reused_["54_aa_ov"]("i,a");
    tmps_["86_aaa_Lvov"]("L,a,i,b") -= l2["aaa"]("L,k,c,b") * eri["aaaa_vovo"]("c,i,a,k");
    tmps_["86_aaa_Lvov"]("L,a,i,b") -= l2_1["bab"]("L,m,b,d") * reused_["68_bbaa_voov"]("d,m,i,a");
    tmps_["86_aaa_Lvov"]("L,a,i,b") += l2["bab"]("L,m,b,d") * eri["baab_vovo"]("d,i,a,m");
    tmps_["86_aaa_Lvov"]("L,a,i,b") -= l2_1["aaa"]("L,k,c,b") * reused_["69_aaaa_vovo"]("c,i,a,k");
    tmps_["86_aaa_Lvov"]("L,a,i,b") += tmps_["84_abb_Lvvo"]("L,b,e,l") * eri["abab_oovv"]("i,l,a,e");
    tmps_["86_aaa_Lvov"]("L,a,i,b") += eri["abab_oovo"]("i,l,a,m") * tmps_["65_bab_Lovo"]("L,m,b,l");
    tmps_["86_aaa_Lvov"]("L,a,i,b") += dp["aa_ov"]("i,a") * tmps_["67_a_Lv"]("L,b");
    tmps_["84_abb_Lvvo"].~TArrayD();
    tmps_["83_aaa_Lovo"].~TArrayD();

    // sigmal2_aaa += +1.00 P(a,b) d-_aa(j,a) t1_1_aa(c,j) l2_aaa(i,c,b)
    //             += -1.00 P(a,b) f_aa(j,a) t1_1_aa(c,j) l2_1_aaa(i,c,b)
    //             += +1.00 P(a,b) d-_aa(j,a) t0_1 t1_1_aa(c,j) l2_1_aaa(i,c,b)
    //             += -0.50 P(a,b) <k,j||a,d>_abab t2_abab(c,d,k,j) l2_aaa(i,c,b)
    //             += -0.50 P(a,b) <j,k||a,d>_abab t2_abab(c,d,j,k) l2_aaa(i,c,b)
    //             += +1.00 P(a,b) <i,k||d,a>_aaaa t2_abab(d,c,k,j) l2_bab(j,b,c)
    //             += -1.00 P(a,b) <i,c||a,d>_abab t1_1_bb(d,j) l2_1_bab(j,b,c)
    //             += +1.00 P(a,b) <i,k||a,d>_abab t2_1_bbbb(d,c,j,k) l2_1_bab(j,b,c)
    //             += +1.00 P(a,b) <i,k||a,d>_abab t2_bbbb(d,c,j,k) l2_bab(j,b,c)
    //             += -1.00 P(a,b) f_aa(i,a) l1_a(b)
    //             += +1.00 P(a,b) <i,k||d,a>_aaaa t2_aaaa(d,c,j,k) l2_aaa(j,c,b)
    //             += +1.00 P(a,b) <j,c||d,a>_aaaa t1_1_aa(d,j) l2_1_aaa(i,c,b)
    //             += -0.50 P(a,b) <k,j||d,a>_aaaa t2_1_aaaa(d,c,k,j) l2_1_aaa(i,c,b)
    //             += +1.00 P(a,b) d-_aa(i,a) t0_1 l1_a(b)
    //             += -1.00 P(a,b) d-_aa(c,a) t0_1 l2_aaa(i,c,b)
    //             += -1.00 P(a,b) <i,j||a,c>_abab t1_1_bb(c,j) l1_1_a(b)
    //             += +1.00 P(a,b) d+_aa(i,a) l1_1_a(b)
    //             += -0.50 P(a,b) <k,j||d,a>_aaaa t2_aaaa(d,c,k,j) l2_aaa(i,c,b)
    //             += -1.00 P(a,b) d+_aa(c,a) l2_1_aaa(i,c,b)
    //             += -0.50 P(a,b) <k,j||a,d>_abab t2_1_abab(c,d,k,j) l2_1_aaa(i,c,b)
    //             += -0.50 P(a,b) <j,k||a,d>_abab t2_1_abab(c,d,j,k) l2_1_aaa(i,c,b)
    //             += +1.00 P(a,b) <c,j||a,d>_abab t1_1_bb(d,j) l2_1_aaa(i,c,b)
    //             += +1.00 P(a,b) f_aa(c,a) l2_aaa(i,c,b)
    //             += +1.00 P(a,b) <i,j||c,a>_aaaa t1_1_aa(c,j) l1_1_a(b)
    //             += +1.00 P(a,b) <i,c||a,j>_aaaa l2_aaa(j,c,b)
    //             += +1.00 P(a,b) <i,k||d,a>_aaaa t2_1_abab(d,c,k,j) l2_1_bab(j,b,c)
    //             += -1.00 P(a,b) <i,c||a,j>_abab l2_bab(j,b,c)
    //             += -1.00 P(a,b) <i,c||d,a>_aaaa t1_1_aa(d,j) l2_1_aaa(j,c,b)
    //             += +1.00 P(a,b) <i,k||d,a>_aaaa t2_1_aaaa(d,c,j,k) l2_1_aaa(j,c,b)
    //             += +1.00 P(a,b) <i,k||a,d>_abab t2_1_abab(c,d,j,k) l2_1_aaa(j,c,b)
    //             += +1.00 P(a,b) <i,k||a,d>_abab t2_abab(c,d,j,k) l2_aaa(j,c,b)
    //             += +1.00 P(a,b) <i,k||a,j>_abab t1_1_bb(c,k) l2_1_bab(j,b,c)
    //             += +1.00 P(a,b) d-_aa(i,a) t1_1_bb(c,j) l2_bab(j,b,c)
    //             += -1.00 P(a,b) <i,k||a,j>_aaaa t1_1_aa(c,k) l2_1_aaa(j,c,b)
    //             += -1.00 P(a,b) d-_aa(i,a) t1_1_aa(c,j) l2_aaa(j,c,b)
    sigmal2_aaa("L,a,b,i") += tmps_["86_aaa_Lvov"]("L,a,i,b");
    sigmal2_aaa("L,a,b,i") -= tmps_["86_aaa_Lvov"]("L,b,i,a");
    tmps_["86_aaa_Lvov"].~TArrayD();

    // flops: o1v2L1  = o2v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["87_abb_Lvvo"]("L,a,b,i")  = l2_1["bab"]("L,j,a,b") * reused_["59_bb_oo"]("i,j");

    // sigmal2_1_abb += +2.00 d-_bb(i,c) t1_1_bb(c,j) l2_1_bab(j,a,b)
    sigmal2_1_abb("L,a,b,i") += 2.00 * tmps_["87_abb_Lvvo"]("L,a,b,i");

    // flops: o1v2L1  = o2v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["92_bab_Lovv"]("L,i,a,b")  = -0.50 * eri["abab_oovv"]("j,i,a,b") * tmps_["38_a_Lo"]("L,j");
    tmps_["38_a_Lo"].~TArrayD();

    // flops: o1v2L1  = o2v3L1 o2v3L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1
    tmps_["91_bba_Lvvo"]("L,a,b,i")  = l2["bab"]("L,j,c,a") * t2["abab"]("c,b,i,j");
    tmps_["91_bba_Lvvo"]("L,a,b,i") += l2_1["bab"]("L,j,c,a") * t2_1["abab"]("c,b,i,j");

    // flops: o3v0L1  = o3v2L1 o3v2L1 o3v0L1
    //  mems: o3v0L1  = o3v0L1 o3v0L1 o3v0L1
    tmps_["90_abb_Looo"]("L,i,j,k")  = t2["abab"]("a,b,i,j") * l2["bab"]("L,k,a,b");
    tmps_["90_abb_Looo"]("L,i,j,k") += t2_1["abab"]("a,b,i,j") * l2_1["bab"]("L,k,a,b");

    // flops: o2v1L1  = o2v2L1
    //  mems: o2v1L1  = o2v1L1
    tmps_["89_bab_Lovo"]("L,i,a,j")  = l2["bab"]("L,i,a,b") * t1_1["bb"]("b,j");

    // flops: o2v1L1  = o2v2L1
    //  mems: o2v1L1  = o2v1L1
    tmps_["88_abb_Loov"]("L,i,j,a")  = t1_1["aa"]("b,i") * l2["bab"]("L,j,b,a");

    // flops: o1v2L1  = o2v3L1 o2v2L1 o2v2L1 o2v3L1 o2v3L1 o3v2L1 o1v4L1 o2v2L1 o1v3L1 o2v3L1 o2v3L1 o2v2L1 o2v3L1 o2v3L1 o1v3L1 o1v3L1 o1v3L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o2v2L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o3v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o3v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o3v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["93_bba_Lovv"]("L,i,a,b")  = -1.00 * eri["baab_vovv"]("f,m,b,a") * tmps_["62_abb_Loov"]("L,m,i,f");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= dp["aa_ov"]("m,b") * tmps_["88_abb_Loov"]("L,m,i,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= 0.50 * eri["abab_oovv"]("l,i,b,a") * tmps_["39_a_Lo"]("L,l");
    tmps_["93_bba_Lovv"]("L,i,a,b") += eri["abab_vovv"]("e,k,b,a") * tmps_["65_bab_Lovo"]("L,i,e,k");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= eri["abab_oovv"]("l,i,b,f") * tmps_["91_bba_Lvvo"]("L,a,f,l");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= eri["abab_oovv"]("l,k,b,a") * tmps_["90_abb_Looo"]("L,l,k,i");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= eri["abab_vvvv"]("e,c,b,a") * l2["bab"]("L,i,e,c");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= l2["bab"]("L,k,b,a") * reused_["67_bb_oo"]("i,k");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= reused_["78_bb_vv"]("c,a") * l2_1["bab"]("L,i,b,c");
    tmps_["93_bba_Lovv"]("L,i,a,b") += eri["bbbb_vovo"]("c,i,a,k") * l2["bab"]("L,k,b,c");
    tmps_["93_bba_Lovv"]("L,i,a,b") += l2_1["bab"]("L,k,b,c") * reused_["76_bbbb_voov"]("c,k,i,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= dp["bb_oo"]("i,k") * l2_1["bab"]("L,k,b,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") += l2_1["bab"]("L,k,d,a") * reused_["83_abab_vovo"]("d,i,b,k");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= eri["abba_vovo"]("d,i,a,m") * l2["aaa"]("L,m,d,b");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= f["aa_vv"]("d,b") * l2["bab"]("L,i,d,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") += reused_["66_bb_vv"]("c,a") * l2["bab"]("L,i,b,c");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= eri["abab_vovv"]("d,i,b,a") * l1["a"]("L,d");
    tmps_["93_bba_Lovv"]("L,i,a,b") += reused_["88_bb_vv"]("c,a") * l2_1["bab"]("L,i,b,c");
    tmps_["93_bba_Lovv"]("L,i,a,b") += l2["bab"]("L,k,b,c") * reused_["82_bbbb_voov"]("c,k,i,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") += 0.50 * reused_["31_aa_vv"]("b,d") * l2["bab"]("L,i,d,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= l2["aaa"]("L,m,d,b") * reused_["47_aabb_voov"]("d,m,i,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= l2_1["bab"]("L,k,b,c") * reused_["80_bbbb_voov"]("c,k,i,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") += reused_["2_aa_vv"]("d,b") * l2["bab"]("L,i,d,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") += l2_1["aaa"]("L,m,d,b") * reused_["74_aabb_voov"]("d,m,i,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") += scalars_["2"] * l2_1["bab"]("L,i,b,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= l2_1["aaa"]("L,m,d,b") * reused_["79_aabb_voov"]("d,m,i,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") += scalars_["5"] * l2["bab"]("L,i,b,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") += dp["bb_ov"]("i,a") * l1_1["a"]("L,b");
    tmps_["93_bba_Lovv"]("L,i,a,b") += l2["aaa"]("L,m,d,b") * reused_["81_aabb_voov"]("d,m,i,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") += scalars_["3"] * l2["bab"]("L,i,b,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") += dp["aa_vv"]("d,b") * l2_1["bab"]("L,i,d,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") += l2_1["bab"]("L,k,b,a") * reused_["87_bb_oo"]("k,i");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= l2["bab"]("L,k,b,c") * reused_["51_bbbb_voov"]("c,k,i,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") += scalars_["4"] * l2["bab"]("L,i,b,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") += l2_1["aaa"]("L,m,d,b") * reused_["84_abba_vovo"]("d,i,a,m");
    tmps_["93_bba_Lovv"]("L,i,a,b") += scalars_["1"] * l2_1["bab"]("L,i,b,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= f["bb_vv"]("c,a") * l2["bab"]("L,i,b,c");
    tmps_["93_bba_Lovv"]("L,i,a,b") += l2["bab"]("L,k,b,a") * f["bb_oo"]("i,k");
    tmps_["93_bba_Lovv"]("L,i,a,b") += l1["a"]("L,b") * reused_["5_bb_ov"]("i,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") += eri["abab_vovo"]("d,i,b,k") * l2["bab"]("L,k,d,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") += l2_1["bab"]("L,k,b,a") * reused_["89_bb_oo"]("i,k");
    tmps_["93_bba_Lovv"]("L,i,a,b") += l2["bab"]("L,i,d,a") * reused_["3_aa_vv"]("d,b");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= f["bb_ov"]("i,a") * l1["a"]("L,b");
    tmps_["93_bba_Lovv"]("L,i,a,b") += l2["bab"]("L,i,b,c") * reused_["55_bb_vv"]("c,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= l1_1["a"]("L,b") * reused_["12_bb_ov"]("i,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= l1_1["a"]("L,b") * reused_["75_bb_ov"]("i,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= scalars_["6"] * l2_1["bab"]("L,i,b,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") += dp["bb_vv"]("c,a") * l2_1["bab"]("L,i,b,c");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= l2_1["bab"]("L,i,d,a") * reused_["61_aa_vv"]("d,b");
    tmps_["93_bba_Lovv"]("L,i,a,b") += 0.50 * l2["bab"]("L,i,b,c") * reused_["85_bb_vv"]("a,c");
    tmps_["93_bba_Lovv"]("L,i,a,b") += l2["bab"]("L,k,b,a") * reused_["57_bb_oo"]("i,k");
    tmps_["93_bba_Lovv"]("L,i,a,b") += 0.50 * l2["bab"]("L,k,b,a") * reused_["86_bb_oo"]("i,k");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= l2_1["bab"]("L,k,b,a") * reused_["77_bb_oo"]("i,k");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= l2["bab"]("L,k,b,a") * reused_["59_bb_oo"]("i,k");
    tmps_["93_bba_Lovv"]("L,i,a,b") += l2_1["bab"]("L,i,d,a") * reused_["13_aa_vv"]("d,b");
    tmps_["93_bba_Lovv"]("L,i,a,b") += eri["abab_oovv"]("l,i,b,a") * tmps_["41_a_Lo"]("L,l");
    tmps_["93_bba_Lovv"]("L,i,a,b") += f["bb_ov"]("k,a") * tmps_["65_bab_Lovo"]("L,i,b,k");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= reused_["5_bb_ov"]("k,a") * tmps_["65_bab_Lovo"]("L,i,b,k");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= t0_1 * tmps_["87_abb_Lvvo"]("L,b,a,i");
    tmps_["93_bba_Lovv"]("L,i,a,b") += eri["abab_oovv"]("l,i,b,a") * tmps_["33_a_Lo"]("L,l");
    tmps_["93_bba_Lovv"]("L,i,a,b") += f["aa_ov"]("m,b") * tmps_["62_abb_Loov"]("L,m,i,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= reused_["4_aa_ov"]("m,b") * tmps_["62_abb_Loov"]("L,m,i,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= eri["abab_oovo"]("l,i,b,k") * tmps_["62_abb_Loov"]("L,l,k,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= dp["bb_ov"]("k,a") * tmps_["89_bab_Lovo"]("L,i,b,k");
    tmps_["93_bba_Lovv"]("L,i,a,b") += eri["abab_oovv"]("m,i,b,a") * tmps_["42_a_Lo"]("L,m");
    tmps_["93_bba_Lovv"]("L,i,a,b") += tmps_["92_bab_Lovv"]("L,i,b,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") += tmps_["61_aaa_Lovo"]("L,m,b,l") * eri["abba_oovo"]("l,i,a,m");
    tmps_["93_bba_Lovv"]("L,i,a,b") += tmps_["67_a_Lv"]("L,b") * dp["bb_ov"]("i,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") -= tmps_["85_a_Lv"]("L,b") * dp["bb_ov"]("i,a");
    tmps_["93_bba_Lovv"]("L,i,a,b") += tmps_["65_bab_Lovo"]("L,k,b,j") * eri["bbbb_oovo"]("i,j,a,k");
    tmps_["92_bab_Lovv"].~TArrayD();
    tmps_["91_bba_Lvvo"].~TArrayD();
    tmps_["90_abb_Looo"].~TArrayD();
    tmps_["89_bab_Lovo"].~TArrayD();
    tmps_["88_abb_Loov"].~TArrayD();
    tmps_["87_abb_Lvvo"].~TArrayD();
    tmps_["85_a_Lv"].~TArrayD();
    tmps_["67_a_Lv"].~TArrayD();
    tmps_["42_a_Lo"].~TArrayD();
    tmps_["41_a_Lo"].~TArrayD();
    tmps_["39_a_Lo"].~TArrayD();
    tmps_["33_a_Lo"].~TArrayD();

    // sigmal2_abb  = +1.00 <k,i||j,b>_abab t1_1_aa(c,k) l2_1_aaa(j,c,a)
    //             += -0.50 <k,i||a,b>_abab t2_abab(d,c,k,j) l2_bab(j,d,c)
    //             += -0.50 <k,i||a,b>_abab t2_abab(c,d,k,j) l2_bab(j,c,d)
    //             += -0.50 <k,j||d,b>_abab t2_1_abab(d,c,k,j) l2_1_bab(i,a,c)
    //             += -0.50 <j,k||d,b>_abab t2_1_abab(d,c,j,k) l2_1_bab(i,a,c)
    //             += +1.00 <j,c||d,b>_abab t1_1_aa(d,j) l2_1_bab(i,a,c)
    //             += +1.00 <c,i||a,b>_abab l1_a(c)
    //             += -1.00 d-_bb(c,b) t0_1 l2_bab(i,a,c)
    //             += +1.00 f_aa(c,a) l2_bab(i,c,b)
    //             += +1.00 <i,k||d,b>_bbbb t2_bbbb(d,c,j,k) l2_bab(j,a,c)
    //             += -0.50 <k,j||d,a>_aaaa t2_aaaa(d,c,k,j) l2_bab(i,c,b)
    //             += +1.00 <k,i||d,b>_abab t2_aaaa(d,c,j,k) l2_aaa(j,c,a)
    //             += +1.00 <k,i||d,b>_abab t2_1_abab(d,c,k,j) l2_1_bab(j,a,c)
    //             += -0.50 <k,j||a,d>_abab t2_abab(c,d,k,j) l2_bab(i,c,b)
    //             += -0.50 <j,k||a,d>_abab t2_abab(c,d,j,k) l2_bab(i,c,b)
    //             += +1.00 <i,k||d,b>_bbbb t2_1_abab(c,d,j,k) l2_1_aaa(j,c,a)
    //             += -1.00 d+_aa(j,j) l2_1_bab(i,a,b)
    //             += +1.00 <k,i||d,b>_abab t2_1_aaaa(d,c,j,k) l2_1_aaa(j,c,a)
    //             += -0.50 <k,j||k,j>_abab l2_bab(i,a,b)
    //             += -0.50 <j,k||j,k>_abab l2_bab(i,a,b)
    //             += +0.25 <k,j||c,d>_bbbb t2_bbbb(c,d,k,j) l2_bab(i,a,b)
    //             += +1.00 f_bb(j,j) l2_bab(i,a,b)
    //             += +0.25 <k,j||c,d>_aaaa t2_aaaa(c,d,k,j) l2_bab(i,a,b)
    //             += -1.00 d-_aa(j,j) t0_1 l2_bab(i,a,b)
    //             += +1.00 f_aa(j,j) l2_bab(i,a,b)
    //             += -0.50 <k,j||k,j>_bbbb l2_bab(i,a,b)
    //             += -0.50 <k,j||k,j>_aaaa l2_bab(i,a,b)
    //             += +0.25 <k,j||c,d>_abab t2_abab(c,d,k,j) l2_bab(i,a,b)
    //             += +0.25 <j,k||c,d>_abab t2_abab(c,d,j,k) l2_bab(i,a,b)
    //             += +0.25 <k,j||d,c>_abab t2_abab(d,c,k,j) l2_bab(i,a,b)
    //             += +0.25 <j,k||d,c>_abab t2_abab(d,c,j,k) l2_bab(i,a,b)
    //             += -1.00 d-_bb(j,j) t0_1 l2_bab(i,a,b)
    //             += -1.00 <c,i||j,b>_abab l2_aaa(j,c,a)
    //             += -1.00 d+_bb(i,b) l1_1_a(a)
    //             += -1.00 <c,i||a,d>_abab t1_1_bb(d,j) l2_1_bab(j,c,b)
    //             += +1.00 d+_bb(i,j) l2_1_bab(j,a,b)
    //             += +1.00 <i,k||d,b>_bbbb t2_1_bbbb(d,c,j,k) l2_1_bab(j,a,c)
    //             += -1.00 <i,c||d,b>_bbbb t1_1_bb(d,j) l2_1_bab(j,a,c)
    //             += +1.00 <i,c||b,j>_bbbb l2_bab(j,a,c)
    //             += +1.00 <j,c||d,b>_bbbb t1_1_bb(d,j) l2_1_bab(i,a,c)
    //             += -0.50 <k,j||d,b>_bbbb t2_1_bbbb(d,c,k,j) l2_1_bab(i,a,c)
    //             += +1.00 d-_bb(i,j) t0_1 l2_bab(j,a,b)
    //             += +1.00 <i,k||d,b>_bbbb t2_abab(c,d,j,k) l2_aaa(j,c,a)
    //             += -1.00 d-_bb(j,c) t1_1_bb(c,j) l2_bab(i,a,b)
    //             += -1.00 d+_aa(c,a) l2_1_bab(i,c,b)
    //             += -1.00 f_bb(i,c) t1_1_bb(c,j) l2_1_bab(j,a,b)
    //             += +1.00 <k,i||d,b>_abab t2_abab(d,c,k,j) l2_bab(j,a,c)
    //             += -1.00 d-_aa(j,c) t1_1_aa(c,j) l2_bab(i,a,b)
    //             += -1.00 <c,i||d,b>_abab t1_1_aa(d,j) l2_1_aaa(j,c,a)
    //             += -1.00 d+_bb(j,j) l2_1_bab(i,a,b)
    //             += +1.00 f_bb(c,b) l2_bab(i,a,c)
    //             += -1.00 f_bb(i,j) l2_bab(j,a,b)
    //             += -1.00 d-_bb(i,b) t0_1 l1_a(a)
    //             += +0.50 <d,c||a,b>_abab l2_bab(i,d,c)
    //             += +0.50 <c,d||a,b>_abab l2_bab(i,c,d)
    //             += -1.00 <c,i||a,j>_abab l2_bab(j,c,b)
    //             += -1.00 <k,i||c,j>_abab t1_1_aa(c,k) l2_1_bab(j,a,b)
    //             += -0.50 <k,i||c,d>_abab t2_1_abab(c,d,k,j) l2_1_bab(j,a,b)
    //             += -0.50 <k,i||d,c>_abab t2_1_abab(d,c,k,j) l2_1_bab(j,a,b)
    //             += -1.00 d-_aa(c,a) t0_1 l2_bab(i,c,b)
    //             += +1.00 f_bb(i,b) l1_a(a)
    //             += -0.50 <k,j||d,b>_abab t2_abab(d,c,k,j) l2_bab(i,a,c)
    //             += -0.50 <j,k||d,b>_abab t2_abab(d,c,j,k) l2_bab(i,a,c)
    //             += +1.00 <j,i||c,b>_abab t1_1_aa(c,j) l1_1_a(a)
    //             += -1.00 <i,j||c,b>_bbbb t1_1_bb(c,j) l1_1_a(a)
    //             += +1.00 t0_1 l2_1_bab(i,a,b) w0
    //             += +0.25 <k,j||c,d>_abab t2_1_abab(c,d,k,j) l2_1_bab(i,a,b)
    //             += +0.25 <j,k||c,d>_abab t2_1_abab(c,d,j,k) l2_1_bab(i,a,b)
    //             += +0.25 <k,j||d,c>_abab t2_1_abab(d,c,k,j) l2_1_bab(i,a,b)
    //             += +0.25 <j,k||d,c>_abab t2_1_abab(d,c,j,k) l2_1_bab(i,a,b)
    //             += -1.00 d-_bb(j,c) t0_1 t1_1_bb(c,j) l2_1_bab(i,a,b)
    //             += +0.25 <k,j||c,d>_bbbb t2_1_bbbb(c,d,k,j) l2_1_bab(i,a,b)
    //             += -1.00 d-_aa(j,c) t0_1 t1_1_aa(c,j) l2_1_bab(i,a,b)
    //             += +1.00 f_bb(j,c) t1_1_bb(c,j) l2_1_bab(i,a,b)
    //             += +0.25 <k,j||c,d>_aaaa t2_1_aaaa(c,d,k,j) l2_1_bab(i,a,b)
    //             += +1.00 f_aa(j,c) t1_1_aa(c,j) l2_1_bab(i,a,b)
    //             += -1.00 d+_bb(c,b) l2_1_bab(i,a,c)
    //             += +1.00 <j,c||d,a>_aaaa t1_1_aa(d,j) l2_1_bab(i,c,b)
    //             += -0.50 <k,j||d,a>_aaaa t2_1_aaaa(d,c,k,j) l2_1_bab(i,c,b)
    //             += -0.50 <k,j||d,b>_bbbb t2_bbbb(d,c,k,j) l2_bab(i,a,c)
    //             += -0.50 <k,i||c,d>_abab t2_abab(c,d,k,j) l2_bab(j,a,b)
    //             += -0.50 <k,i||d,c>_abab t2_abab(d,c,k,j) l2_bab(j,a,b)
    //             += -0.50 <i,k||c,d>_bbbb t2_bbbb(c,d,j,k) l2_bab(j,a,b)
    //             += +1.00 <i,k||c,j>_bbbb t1_1_bb(c,k) l2_1_bab(j,a,b)
    //             += -0.50 <i,k||c,d>_bbbb t2_1_bbbb(c,d,j,k) l2_1_bab(j,a,b)
    //             += +1.00 d-_bb(i,c) t1_1_bb(c,j) l2_bab(j,a,b)
    //             += -0.50 <k,j||a,d>_abab t2_1_abab(c,d,k,j) l2_1_bab(i,c,b)
    //             += -0.50 <j,k||a,d>_abab t2_1_abab(c,d,j,k) l2_1_bab(i,c,b)
    //             += +1.00 <c,j||a,d>_abab t1_1_bb(d,j) l2_1_bab(i,c,b)
    //             += +0.25 <k,j||a,b>_abab t2_1_abab(d,c,k,j) l2_1_bab(i,d,c)
    //             += +0.25 <j,k||a,b>_abab t2_1_abab(d,c,j,k) l2_1_bab(i,d,c)
    //             += +0.25 <k,j||a,b>_abab t2_1_abab(c,d,k,j) l2_1_bab(i,c,d)
    //             += +0.25 <j,k||a,b>_abab t2_1_abab(c,d,j,k) l2_1_bab(i,c,d)
    //             += +0.25 <k,j||a,b>_abab t2_abab(d,c,k,j) l2_bab(i,d,c)
    //             += +0.25 <j,k||a,b>_abab t2_abab(d,c,j,k) l2_bab(i,d,c)
    //             += +0.25 <k,j||a,b>_abab t2_abab(c,d,k,j) l2_bab(i,c,d)
    //             += +0.25 <j,k||a,b>_abab t2_abab(c,d,j,k) l2_bab(i,c,d)
    //             += -1.00 f_bb(j,b) t1_1_bb(c,j) l2_1_bab(i,a,c)
    //             += +1.00 d-_bb(j,b) t0_1 t1_1_bb(c,j) l2_1_bab(i,a,c)
    //             += +1.00 d-_bb(i,c) t0_1 t1_1_bb(c,j) l2_1_bab(j,a,b)
    //             += +1.00 <k,i||a,d>_abab t2_abab(c,d,k,j) l2_bab(j,c,b)
    //             += +1.00 <k,i||a,d>_abab t2_1_abab(c,d,k,j) l2_1_bab(j,c,b)
    //             += -1.00 <d,j||a,b>_abab t1_1_bb(c,j) l2_1_bab(i,d,c)
    //             += +0.50 <k,i||a,b>_abab t2_aaaa(d,c,j,k) l2_aaa(j,d,c)
    //             += +1.00 d-_aa(j,a) t1_1_aa(c,j) l2_bab(i,c,b)
    //             += -0.50 <k,i||a,b>_abab t2_1_abab(d,c,k,j) l2_1_bab(j,d,c)
    //             += -0.50 <k,i||a,b>_abab t2_1_abab(c,d,k,j) l2_1_bab(j,c,d)
    //             += -1.00 f_aa(j,a) t1_1_aa(c,j) l2_1_bab(i,c,b)
    //             += +1.00 d-_aa(j,a) t0_1 t1_1_aa(c,j) l2_1_bab(i,c,b)
    //             += +1.00 <k,i||a,j>_abab t1_1_aa(c,k) l2_1_bab(j,c,b)
    //             += +1.00 d-_bb(j,b) t1_1_bb(c,j) l2_bab(i,a,c)
    //             += -1.00 <j,d||a,b>_abab t1_1_aa(c,j) l2_1_bab(i,c,d)
    //             += -1.00 <j,i||a,b>_abab t1_1_aa(c,j) l1_1_a(c)
    //             += +0.50 <k,i||a,b>_abab t2_1_aaaa(d,c,j,k) l2_1_aaa(j,d,c)
    //             += -1.00 d-_bb(i,b) t1_1_bb(c,j) l2_bab(j,a,c)
    //             += +1.00 d-_bb(i,b) t1_1_aa(c,j) l2_aaa(j,c,a)
    //             += -1.00 <i,k||b,j>_bbbb t1_1_bb(c,k) l2_1_bab(j,a,c)
    sigmal2_abb("L,a,b,i")  = -1.00 * tmps_["93_bba_Lovv"]("L,i,b,a");
    tmps_["93_bba_Lovv"].~TArrayD();

    // flops: o0v1L1  = o1v2L1
    //  mems: o0v1L1  = o0v1L1
    tmps_["97_a_Lv"]("L,a")  = l2_1["aaa"]("L,i,b,a") * t1_1["aa"]("b,i");

    // csigmal1_1_a += -1.00 t1_1_aa(b,i) l2_1_aaa(i,b,a)
    csigmal1_1_a("L,a") -= tmps_["97_a_Lv"]("L,a");

    // flops: o1v2L1  = o2v2L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["96_bba_Lvov"]("L,a,i,b")  = -2.00 * dp["bb_ov"]("j,a") * tmps_["65_bab_Lovo"]("L,i,b,j");
    tmps_["65_bab_Lovo"].~TArrayD();

    // flops: o1v2L1  = o2v3L1
    //  mems: o1v2L1  = o1v2L1
    tmps_["95_bba_Lvvo"]("L,a,b,i")  = l2_1["bab"]("L,j,c,a") * t2["abab"]("c,b,i,j");

    // flops: o3v0L1  = o3v2L1
    //  mems: o3v0L1  = o3v0L1
    tmps_["94_abb_Looo"]("L,i,j,k")  = t2["abab"]("a,b,i,j") * l2_1["bab"]("L,k,a,b");

    // flops: o1v2L1  = o1v2L1 o2v3L1 o1v2L1 o2v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v3L1 o1v2L1 o2v2L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o1v2L1 o2v3L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v4L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o2v2L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v3L1 o1v2L1 o1v2L1 o3v2L1 o1v2L1 o2v2L1 o1v2L1 o2v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    //  mems: o1v2L1  = o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1 o1v2L1
    tmps_["98_bba_Lovv"]("L,i,a,b")  = -1.00 * tmps_["53_a_Lv"]("L,b") * dp["bb_ov"]("i,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") += tmps_["95_bba_Lvvo"]("L,a,c,j") * eri["abab_oovv"]("j,i,b,c");
    tmps_["98_bba_Lovv"]("L,i,a,b") += tmps_["97_a_Lv"]("L,b") * dp["bb_ov"]("i,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= eri["abab_oovv"]("j,i,b,a") * tmps_["20_a_Lo"]("L,j");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= l2_1["bab"]("L,k,b,a") * f["bb_oo"]("i,k");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= 2.00 * scalars_["4"] * l2_1["bab"]("L,i,b,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= scalars_["2"] * l2["bab"]("L,i,b,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") += l2_1["aaa"]("L,l,d,b") * eri["abba_vovo"]("d,i,a,l");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= l2_1["aaa"]("L,l,d,b") * reused_["81_aabb_voov"]("d,l,i,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= l1_1["a"]("L,b") * reused_["5_bb_ov"]("i,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") += l2_1["aaa"]("L,l,d,b") * reused_["47_aabb_voov"]("d,l,i,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= l2_1["bab"]("L,k,b,e") * reused_["82_bbbb_voov"]("e,k,i,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= 0.50 * l2_1["bab"]("L,i,d,a") * reused_["31_aa_vv"]("b,d");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= dp["bb_vv"]("e,a") * l2["bab"]("L,i,b,e");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= l2_1["bab"]("L,i,b,e") * reused_["66_bb_vv"]("e,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= l2_1["bab"]("L,i,d,a") * reused_["2_aa_vv"]("d,b");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= l2_1["bab"]("L,i,b,e") * reused_["55_bb_vv"]("e,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= 2.00 * scalars_["3"] * l2_1["bab"]("L,i,b,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") += l1_1["a"]("L,b") * f["bb_ov"]("i,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") += l2_1["bab"]("L,k,b,e") * reused_["51_bbbb_voov"]("e,k,i,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= l2_1["bab"]("L,i,d,a") * reused_["3_aa_vv"]("d,b");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= 0.50 * l2_1["bab"]("L,k,b,a") * reused_["86_bb_oo"]("i,k");
    tmps_["98_bba_Lovv"]("L,i,a,b") += f["bb_vv"]("e,a") * l2_1["bab"]("L,i,b,e");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= scalars_["1"] * l2["bab"]("L,i,b,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= l2_1["bab"]("L,k,b,e") * eri["bbbb_vovo"]("e,i,a,k");
    tmps_["98_bba_Lovv"]("L,i,a,b") += w0 * l2_1["bab"]("L,i,b,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= l1["a"]("L,b") * dp["bb_ov"]("i,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") += eri["abab_vvvv"]("f,e,b,a") * l2_1["bab"]("L,i,f,e");
    tmps_["98_bba_Lovv"]("L,i,a,b") += l2["bab"]("L,k,b,a") * dp["bb_oo"]("i,k");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= l2_1["bab"]("L,k,b,a") * reused_["57_bb_oo"]("i,k");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= scalars_["5"] * l2_1["bab"]("L,i,b,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") += l2_1["bab"]("L,k,b,a") * reused_["67_bb_oo"]("i,k");
    tmps_["98_bba_Lovv"]("L,i,a,b") += f["aa_vv"]("d,b") * l2_1["bab"]("L,i,d,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") += eri["abab_vovv"]("d,i,b,a") * l1_1["a"]("L,d");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= dp["aa_vv"]("d,b") * l2["bab"]("L,i,d,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= 0.50 * l2_1["bab"]("L,i,b,e") * reused_["85_bb_vv"]("a,e");
    tmps_["98_bba_Lovv"]("L,i,a,b") += eri["abab_oovv"]("j,k,b,a") * tmps_["94_abb_Looo"]("L,j,k,i");
    tmps_["98_bba_Lovv"]("L,i,a,b") += 2.00 * dp["aa_ov"]("l,b") * tmps_["62_abb_Loov"]("L,l,i,a");
    tmps_["98_bba_Lovv"]("L,i,a,b") += 0.50 * eri["abab_oovv"]("j,i,b,a") * tmps_["19_a_Lo"]("L,j");
    tmps_["98_bba_Lovv"]("L,i,a,b") -= tmps_["96_bba_Lvov"]("L,a,i,b");
    tmps_["96_bba_Lvov"].~TArrayD();
    tmps_["95_bba_Lvvo"].~TArrayD();
    tmps_["94_abb_Looo"].~TArrayD();
    tmps_["62_abb_Loov"].~TArrayD();
    tmps_["20_a_Lo"].~TArrayD();
    tmps_["19_a_Lo"].~TArrayD();
}