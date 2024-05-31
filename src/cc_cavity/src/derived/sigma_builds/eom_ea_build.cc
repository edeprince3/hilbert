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

#include "cc_cavity/include/derived/eom_ea_driver.h"

void hilbert::EOM_EA_CCSD::build_common_ops() {
    
    if (!cc_wfn_->has_t1_integrals_) cc_wfn_->transform_integrals(true);

    TArrayMap &eri = cc_wfn_->V_blks_;
    TArrayMap &f = cc_wfn_->F_blks_;

    TArrayMap t2 = {
        {"aaaa_vvoo", cc_wfn_->amplitudes_["t2_aaaa"]},
        {"abab_vvoo", cc_wfn_->amplitudes_["t2_abab"]},
        {"bbbb_vvoo", cc_wfn_->amplitudes_["t2_bbbb"]}
    };

    {
        reuse_tmps_["1_bbbb_vvoo"]("a,f,i,m")  = eri["abab_oovv"]("j,m,b,f") * t2["abab_vvoo"]("b,a,j,i");
        reuse_tmps_["2_aaaa_vvoo"]("a,e,i,m")  = eri["abab_oovv"]("m,j,e,b") * t2["abab_vvoo"]("a,b,i,j");
        reuse_tmps_["3_bbbb_vvoo"]("f,b,m,j")  = eri["bbbb_oovv"]("j,i,a,b") * t2["bbbb_vvoo"]("a,f,m,i");
        reuse_tmps_["4_aaaa_vvoo"]("f,b,m,j")  = eri["aaaa_oovv"]("j,i,a,b") * t2["aaaa_vvoo"]("a,f,m,i");
        reuse_tmps_["5_abab_vvoo"]("f,a,m,i")  = eri["abab_oovv"]("m,j,f,b") * t2["bbbb_vvoo"]("b,a,i,j");
        reuse_tmps_["6_abab_vvoo"]("e,a,m,i")  = eri["aaaa_oovv"]("m,j,b,e") * t2["abab_vvoo"]("b,a,j,i");
        reuse_tmps_["7_abab_vvoo"]("a,f,i,m")  = eri["abab_oovv"]("j,m,b,f") * t2["aaaa_vvoo"]("b,a,i,j");
        reuse_tmps_["8_abab_vvoo"]("a,f,i,m")  = eri["bbbb_oovv"]("m,j,b,f") * t2["abab_vvoo"]("a,b,i,j");
        reuse_tmps_["9_aabb_vvoo"]("e,b,m,j")  = eri["abab_oovv"]("i,j,b,a") * t2["abab_vvoo"]("e,a,i,m");
        reuse_tmps_["10_bbaa_vvoo"]("f,b,m,j")  = eri["abab_oovv"]("j,i,a,b") * t2["abab_vvoo"]("a,f,m,i");
        reuse_tmps_["11_abba_vvvo"]("b,a,e,i")  = eri["abab_vovv"]("b,j,c,e") * t2["abab_vvoo"]("c,a,i,j");
        reuse_tmps_["12_bbbb_vvvo"]("a,e,b,i")  = eri["baab_vovv"]("b,j,c,e") * t2["abab_vvoo"]("c,a,j,i");
        reuse_tmps_["13_bbbb_vvvo"]("a,e,b,i")  = eri["bbbb_vovv"]("b,j,c,e") * t2["bbbb_vvoo"]("c,a,i,j");
        reuse_tmps_["14_aaaa_vvvo"]("a,e,b,i")  = eri["aaaa_vovv"]("b,j,c,e") * t2["aaaa_vvoo"]("c,a,i,j");
        reuse_tmps_["15_aaaa_vvvo"]("e,b,f,m")  = eri["abab_vovv"]("f,i,b,a") * t2["abab_vvoo"]("e,a,m,i");
        reuse_tmps_["16_aabb_vvvo"]("a,e,b,i")  = eri["baab_vovv"]("b,j,e,c") * t2["abab_vvoo"]("a,c,j,i");
        reuse_tmps_["17_abba_vvvo"]("a,e,b,i")  = eri["baab_vovv"]("b,j,c,e") * t2["aaaa_vvoo"]("c,a,i,j");
        reuse_tmps_["18_aabb_vvvo"]("e,b,a,i")  = eri["abab_vovv"]("b,j,e,c") * t2["bbbb_vvoo"]("c,a,i,j");
        reuse_tmps_["19_aabb_vvvo"]("e,b,a,i")  = eri["aaaa_vovv"]("b,j,c,e") * t2["abab_vvoo"]("c,a,j,i");
        reuse_tmps_["20_abba_vvvo"]("e,b,f,m")  = eri["bbbb_vovv"]("f,i,a,b") * t2["abab_vvoo"]("e,a,m,i");
        reuse_tmps_["21_aa_vv"]("a,e")  = eri["abab_oovv"]("j,i,e,b") * t2["abab_vvoo"]("a,b,j,i");
        reuse_tmps_["22_bb_vv"]("a,e")  = eri["abab_oovv"]("j,i,b,e") * t2["abab_vvoo"]("b,a,j,i");
        reuse_tmps_["23_bb_oo"]("m,j")  = eri["abab_oovv"]("i,j,a,b") * t2["abab_vvoo"]("a,b,i,m");
        reuse_tmps_["24_aa_oo"]("m,j")  = eri["abab_oovv"]("j,i,a,b") * t2["abab_vvoo"]("a,b,m,i");
        reuse_tmps_["25_bb_vv"]("e,b")  = 0.50 * eri["bbbb_oovv"]("j,i,a,b") * t2["bbbb_vvoo"]("a,e,j,i");
        reuse_tmps_["26_aa_vv"]("e,b")  = 0.50 * eri["aaaa_oovv"]("j,i,a,b") * t2["aaaa_vvoo"]("a,e,j,i");
        reuse_tmps_["27_aa_oo"]("m,j")  = 0.50 * eri["aaaa_oovv"]("j,i,a,b") * t2["aaaa_vvoo"]("a,b,m,i");
        reuse_tmps_["28_bb_oo"]("m,j")  = 0.50 * eri["bbbb_oovv"]("j,i,a,b") * t2["bbbb_vvoo"]("a,b,m,i");
        reuse_tmps_["29_abba_vvvo"]("b,a,e,i")  = eri["abba_oovo"]("k,j,e,i") * t2["abab_vvoo"]("b,a,k,j");
        reuse_tmps_["30_bbbb_vvvo"]("a,b,e,i")  = eri["bbbb_oovo"]("k,j,e,i") * t2["bbbb_vvoo"]("b,a,k,j");
        reuse_tmps_["31_aabb_vvvo"]("b,e,a,i")  = eri["abab_oovo"]("k,j,e,i") * t2["abab_vvoo"]("b,a,k,j");
        reuse_tmps_["32_aaaa_vvvo"]("a,b,e,i")  = eri["aaaa_oovo"]("k,j,e,i") * t2["aaaa_vvoo"]("b,a,k,j");
        reuse_tmps_["33_bb_vo"]("a,i")  = f["aa_ov"]("j,b") * t2["abab_vvoo"]("b,a,j,i");
        reuse_tmps_["34_bb_vo"]("f,m")  = f["bb_ov"]("i,a") * t2["bbbb_vvoo"]("a,f,m,i");
        reuse_tmps_["35_aa_vo"]("f,m")  = f["aa_ov"]("i,a") * t2["aaaa_vvoo"]("a,f,m,i");
        reuse_tmps_["36_aa_vo"]("f,m")  = f["bb_ov"]("i,a") * t2["abab_vvoo"]("f,a,m,i");
        reuse_tmps_["37_bb_vo"]("f,m")  = 0.50 * eri["bbbb_vovv"]("f,i,a,b") * t2["bbbb_vvoo"]("a,b,m,i");
        reuse_tmps_["38_bb_vo"]("f,m")  = eri["baab_vovv"]("f,i,a,b") * t2["abab_vvoo"]("a,b,i,m");
        reuse_tmps_["39_aa_vo"]("e,m")  = 0.50 * eri["aaaa_vovv"]("e,i,a,b") * t2["aaaa_vvoo"]("a,b,m,i");
        reuse_tmps_["40_aa_vo"]("f,m")  = eri["abab_vovv"]("f,i,a,b") * t2["abab_vvoo"]("a,b,m,i");
        reuse_tmps_["41_bb_vo"]("f,m")  = eri["abab_oovo"]("j,i,a,m") * t2["abab_vvoo"]("a,f,j,i");
        reuse_tmps_["42_bb_vo"]("f,m")  = 0.50 * eri["bbbb_oovo"]("j,i,a,m") * t2["bbbb_vvoo"]("a,f,j,i");
        reuse_tmps_["43_aa_vo"]("f,m")  = 0.50 * eri["aaaa_oovo"]("j,i,a,m") * t2["aaaa_vvoo"]("a,f,j,i");
        reuse_tmps_["44_aa_vo"]("e,m")  = eri["abba_oovo"]("j,i,a,m") * t2["abab_vvoo"]("e,a,j,i");

    }

    world_.gop.fence();
}

void hilbert::EOM_EA_CCSD::build_Hc_cH(size_t L) {

    // ensure that the integrals are transformed with t1 amplitudes
    if (!cc_wfn_->has_t1_integrals_)
        cc_wfn_->transform_integrals(true);

    /// build sigma vectors for this trial
    sigvec_blks_.clear();

    // reduce sigma into spins
    sigvec_blks_["sigmar1_a"] = HelperD::makeTensor(world_, {L, va_}, true);
    sigvec_blks_["sigmar1_b"] = HelperD::makeTensor(world_, {L, vb_}, true);
    sigvec_blks_["sigmar2_aaa"] = HelperD::makeTensor(world_, {L, va_, va_, oa_}, true);
    sigvec_blks_["sigmar2_abb"] = HelperD::makeTensor(world_, {L, va_, vb_, ob_}, true);
    sigvec_blks_["sigmar2_aba"] = HelperD::makeTensor(world_, {L, va_, vb_, oa_}, true);
    sigvec_blks_["sigmar2_bbb"] = HelperD::makeTensor(world_, {L, vb_, vb_, ob_}, true);

    sigvec_blks_["sigmal1_a"] = HelperD::makeTensor(world_, {L, va_}, true);
    sigvec_blks_["sigmal1_b"] = HelperD::makeTensor(world_, {L, vb_}, true);
    sigvec_blks_["sigmal2_aaa"] = HelperD::makeTensor(world_, {L, va_, va_, oa_}, true);
    sigvec_blks_["sigmal2_abb"] = HelperD::makeTensor(world_, {L, va_, vb_, ob_}, true);
    sigvec_blks_["sigmal2_aba"] = HelperD::makeTensor(world_, {L, va_, vb_, oa_}, true);
    sigvec_blks_["sigmal2_bbb"] = HelperD::makeTensor(world_, {L, vb_, vb_, ob_}, true);

    // get reference to electronic integrals
    TArrayMap &V_blks_ = cc_wfn_->V_blks_;
    TArrayMap &F_blks_ = cc_wfn_->F_blks_;
    TArrayMap &Id_blks_ = cc_wfn_->Id_blks_;
    TArrayMap &eri = cc_wfn_->V_blks_;
    TArrayMap &f = cc_wfn_->F_blks_;
    TArrayMap &Id = cc_wfn_->Id_blks_;

    TArrayMap t2 = {
            {"aaaa_vvoo", cc_wfn_->amplitudes_["t2_aaaa"]},
            {"abab_vvoo", cc_wfn_->amplitudes_["t2_abab"]},
            {"bbbb_vvoo", cc_wfn_->amplitudes_["t2_bbbb"]}
    };

    /// reference right operators
    TArrayMap r1 = {
            {"a_Lv", evec_blks_["r1_a"]},
            {"b_Lv", evec_blks_["r1_b"]}
    };

    TArrayMap r2 = {
            {"aaa_Lvvo", evec_blks_["r2_aaa"]},
            {"abb_Lvvo", evec_blks_["r2_abb"]},
            {"aba_Lvvo", evec_blks_["r2_aba"]},
            {"bbb_Lvvo", evec_blks_["r2_bbb"]}
    };



    /// reference left operators

    TArrayMap l1 = {
            {"a_Lv", evec_blks_["l1_a"]},
            {"b_Lv", evec_blks_["l1_b"]}
    };
    TArrayMap l2 = {
            {"aaa_Lovv", evec_blks_["l2_aaa"]},
            {"bab_Lovv", evec_blks_["l2_bab"]},
            {"aab_Lovv", evec_blks_["l2_aab"]},
            {"bbb_Lovv", evec_blks_["l2_bbb"]}
    };

    /// reference sigma vectors
    auto &sigmar1_a = sigvec_blks_["sigmar1_a"];
    auto &sigmar1_b = sigvec_blks_["sigmar1_b"];
    auto &sigmal1_a = sigvec_blks_["sigmal1_a"];
    auto &sigmal1_b = sigvec_blks_["sigmal1_b"];

    auto &sigmar2_aaa = sigvec_blks_["sigmar2_aaa"];
    auto &sigmar2_abb = sigvec_blks_["sigmar2_abb"];
    auto &sigmar2_aba = sigvec_blks_["sigmar2_aba"];
    auto &sigmar2_bbb = sigvec_blks_["sigmar2_bbb"];
    auto &sigmal2_aaa = sigvec_blks_["sigmal2_aaa"];
    auto &sigmal2_abb = sigvec_blks_["sigmal2_abb"];
    auto &sigmal2_aba = sigvec_blks_["sigmal2_aba"];
    auto &sigmal2_bbb = sigvec_blks_["sigmal2_bbb"];

    {
        TArrayMap tmps_, perm_tmps;
        std::map<std::string, double> scalars_;

        // sigmar1_a = -0.50 <j,i||a,b>_aaaa t2_aaaa(a,e,j,i) r1_a(b)  // flops: o0v1L1 = o0v2L1 | mem: o0v1L1 = o0v1L1
        sigmar1_a("X,e")  = -1.00 * reuse_tmps_["26_aa_vv"]("e,b") * r1["a_Lv"]("X,b");

        // sigmar1_a += -0.50 <j,i||b,a>_abab t2_abab(e,a,j,i) r1_a(b) @ -0.50 <i,j||b,a>_abab t2_abab(e,a,i,j) r1_a(b)  // flops: o0v1L1 += o0v2L1 | mem: o0v1L1 += o0v1L1
        sigmar1_a("X,e") -= reuse_tmps_["21_aa_vv"]("e,b") * r1["a_Lv"]("X,b");

        // sigmar1_a += +1.00 f_aa(e,a) r1_a(a)  // flops: o0v1L1 += o0v2L1 | mem: o0v1L1 += o0v1L1
        sigmar1_a("X,e") += f["aa_vv"]("e,a") * r1["a_Lv"]("X,a");

        // sigmar1_b = -0.50 <j,i||a,b>_bbbb t2_bbbb(a,e,j,i) r1_b(b)  // flops: o0v1L1 = o0v2L1 | mem: o0v1L1 = o0v1L1
        sigmar1_b("X,e")  = -1.00 * reuse_tmps_["25_bb_vv"]("e,b") * r1["b_Lv"]("X,b");

        // sigmar1_b += -0.50 <j,i||a,b>_abab t2_abab(a,e,j,i) r1_b(b) @ -0.50 <i,j||a,b>_abab t2_abab(a,e,i,j) r1_b(b)  // flops: o0v1L1 += o0v2L1 | mem: o0v1L1 += o0v1L1
        sigmar1_b("X,e") -= reuse_tmps_["22_bb_vv"]("e,b") * r1["b_Lv"]("X,b");

        // sigmar1_b += +1.00 f_bb(e,a) r1_b(a)  // flops: o0v1L1 += o0v2L1 | mem: o0v1L1 += o0v1L1
        sigmar1_b("X,e") += f["bb_vv"]("e,a") * r1["b_Lv"]("X,a");

        // sigmal1_a = +1.00 f_aa(a,e) l1_a(a)  // flops: o0v1L1 = o0v2L1 | mem: o0v1L1 = o0v1L1
        sigmal1_a("X,e")  = f["aa_vv"]("a,e") * l1["a_Lv"]("X,a");

        // sigmal1_a += -0.50 <j,i||e,b>_abab l1_a(a) t2_abab(a,b,j,i) @ -0.50 <i,j||e,b>_abab l1_a(a) t2_abab(a,b,i,j)  // flops: o0v1L1 += o0v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") -= reuse_tmps_["21_aa_vv"]("a,e") * l1["a_Lv"]("X,a");

        // sigmal1_a += -0.50 <j,i||b,e>_aaaa l1_a(a) t2_aaaa(b,a,j,i)  // flops: o0v1L1 += o0v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") -= reuse_tmps_["26_aa_vv"]("a,e") * l1["a_Lv"]("X,a");

        // sigmal1_b = +1.00 f_bb(a,e) l1_b(a)  // flops: o0v1L1 = o0v2L1 | mem: o0v1L1 = o0v1L1
        sigmal1_b("X,e")  = f["bb_vv"]("a,e") * l1["b_Lv"]("X,a");

        // sigmal1_b += -0.50 <j,i||b,e>_bbbb l1_b(a) t2_bbbb(b,a,j,i)  // flops: o0v1L1 += o0v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") -= reuse_tmps_["25_bb_vv"]("a,e") * l1["b_Lv"]("X,a");

        // sigmal1_b += -0.50 <j,i||b,e>_abab l1_b(a) t2_abab(b,a,j,i) @ -0.50 <i,j||b,e>_abab l1_b(a) t2_abab(b,a,i,j)  // flops: o0v1L1 += o0v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") -= reuse_tmps_["22_bb_vv"]("a,e") * l1["b_Lv"]("X,a");

        // sigmar1_a += -1.00 f_aa(i,a) r2_aaa(a,e,i)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmar1_a("X,e") -= f["aa_ov"]("i,a") * r2["aaa_Lvvo"]("X,a,e,i");

        // sigmar1_a += +1.00 f_bb(i,a) r2_abb(e,a,i)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmar1_a("X,e") += f["bb_ov"]("i,a") * r2["abb_Lvvo"]("X,e,a,i");

        // sigmar1_b += -1.00 f_aa(i,a) r2_aba(a,e,i)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmar1_b("X,e") -= f["aa_ov"]("i,a") * r2["aba_Lvvo"]("X,a,e,i");

        // sigmar1_b += -1.00 f_bb(i,a) r2_bbb(a,e,i)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmar1_b("X,e") -= f["bb_ov"]("i,a") * r2["bbb_Lvvo"]("X,a,e,i");

        // sigmar2_aaa = -0.50 P(e,f) <e,i||a,b>_abab t2_abab(a,b,m,i) r1_a(f) @ -0.50 P(e,f) <e,i||b,a>_abab t2_abab(b,a,m,i) r1_a(f)  // flops: o1v2L1 = o1v2L1 | mem: o1v2L1 = o1v2L1
        sigmar2_aaa("X,e,f,m")  = -1.00 * reuse_tmps_["40_aa_vo"]("e,m") * r1["a_Lv"]("X,f");
        sigmar2_aaa("X,e,f,m") -= reuse_tmps_["35_aa_vo"]("f,m") * r1["a_Lv"]("X,e");
        sigmar2_aaa("X,e,f,m") += reuse_tmps_["39_aa_vo"]("f,m") * r1["a_Lv"]("X,e");

        // sigmar2_aaa += +0.50 P(e,f) <j,i||m,a>_abab t2_abab(e,a,j,i) r1_a(f) @ +0.50 P(e,f) <i,j||m,a>_abab t2_abab(e,a,i,j) r1_a(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") -= reuse_tmps_["44_aa_vo"]("e,m") * r1["a_Lv"]("X,f");
        sigmar2_aaa("X,e,f,m") += reuse_tmps_["36_aa_vo"]("f,m") * r1["a_Lv"]("X,e");
        sigmar2_aaa("X,e,f,m") += reuse_tmps_["44_aa_vo"]("f,m") * r1["a_Lv"]("X,e");
        sigmar2_aaa("X,e,f,m") += f["aa_vo"]("f,m") * r1["a_Lv"]("X,e");

        // sigmar2_aaa += -1.00 P(e,f) f_aa(e,m) r1_a(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") -= f["aa_vo"]("e,m") * r1["a_Lv"]("X,f");

        // sigmar2_aaa += +0.50 P(e,f) <j,i||a,m>_aaaa t2_aaaa(a,e,j,i) r1_a(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") += reuse_tmps_["43_aa_vo"]("e,m") * r1["a_Lv"]("X,f");

        // sigmar2_aaa += +1.00 P(e,f) f_aa(i,a) t2_aaaa(a,e,m,i) r1_a(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") += reuse_tmps_["35_aa_vo"]("e,m") * r1["a_Lv"]("X,f");
        sigmar2_aaa("X,e,f,m") += reuse_tmps_["40_aa_vo"]("f,m") * r1["a_Lv"]("X,e");

        // sigmar2_aaa += -1.00 P(e,f) f_bb(i,a) t2_abab(e,a,m,i) r1_a(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") -= reuse_tmps_["36_aa_vo"]("e,m") * r1["a_Lv"]("X,f");
        sigmar2_aaa("X,e,f,m") -= reuse_tmps_["43_aa_vo"]("f,m") * r1["a_Lv"]("X,e");

        // sigmar2_aaa += +0.50 P(e,f) <i,e||a,b>_aaaa t2_aaaa(a,b,m,i) r1_a(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") -= reuse_tmps_["39_aa_vo"]("e,m") * r1["a_Lv"]("X,f");

        // sigmar2_aba = -0.50 <e,i||a,b>_abab t2_abab(a,b,m,i) r1_b(f) @ -0.50 <e,i||b,a>_abab t2_abab(b,a,m,i) r1_b(f)  // flops: o1v2L1 = o1v2L1 | mem: o1v2L1 = o1v2L1
        sigmar2_aba("X,e,f,m")  = -1.00 * reuse_tmps_["40_aa_vo"]("e,m") * r1["b_Lv"]("X,f");

        // sigmar2_aba += -1.00 f_bb(i,a) t2_abab(e,a,m,i) r1_b(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= reuse_tmps_["36_aa_vo"]("e,m") * r1["b_Lv"]("X,f");

        // sigmar2_aba += +1.00 f_aa(i,a) t2_aaaa(a,e,m,i) r1_b(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") += reuse_tmps_["35_aa_vo"]("e,m") * r1["b_Lv"]("X,f");

        // sigmar2_aba += +0.50 <j,i||m,a>_abab t2_abab(e,a,j,i) r1_b(f) @ +0.50 <i,j||m,a>_abab t2_abab(e,a,i,j) r1_b(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= reuse_tmps_["44_aa_vo"]("e,m") * r1["b_Lv"]("X,f");

        // sigmar2_aba += -1.00 f_aa(e,m) r1_b(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= f["aa_vo"]("e,m") * r1["b_Lv"]("X,f");

        // sigmar2_aba += +0.50 <i,e||a,b>_aaaa t2_aaaa(a,b,m,i) r1_b(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= reuse_tmps_["39_aa_vo"]("e,m") * r1["b_Lv"]("X,f");

        // sigmar2_aba += +0.50 <j,i||a,m>_aaaa t2_aaaa(a,e,j,i) r1_b(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") += reuse_tmps_["43_aa_vo"]("e,m") * r1["b_Lv"]("X,f");

        // sigmar2_abb = +0.50 <i,f||a,b>_abab t2_abab(a,b,i,m) r1_a(e) @ +0.50 <i,f||b,a>_abab t2_abab(b,a,i,m) r1_a(e)  // flops: o1v2L1 = o1v2L1 | mem: o1v2L1 = o1v2L1
        sigmar2_abb("X,e,f,m")  = -1.00 * reuse_tmps_["38_bb_vo"]("f,m") * r1["a_Lv"]("X,e");

        // sigmar2_abb += +1.00 f_aa(i,a) t2_abab(a,f,i,m) r1_a(e)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") += reuse_tmps_["33_bb_vo"]("f,m") * r1["a_Lv"]("X,e");

        // sigmar2_abb += -0.50 <j,i||a,m>_bbbb t2_bbbb(a,f,j,i) r1_a(e)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") -= reuse_tmps_["42_bb_vo"]("f,m") * r1["a_Lv"]("X,e");

        // sigmar2_abb += -0.50 <j,i||a,m>_abab t2_abab(a,f,j,i) r1_a(e) @ -0.50 <i,j||a,m>_abab t2_abab(a,f,i,j) r1_a(e)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") -= reuse_tmps_["41_bb_vo"]("f,m") * r1["a_Lv"]("X,e");

        // sigmar2_abb += -1.00 f_bb(i,a) t2_bbbb(a,f,m,i) r1_a(e)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") -= reuse_tmps_["34_bb_vo"]("f,m") * r1["a_Lv"]("X,e");

        // sigmar2_abb += -0.50 <i,f||a,b>_bbbb t2_bbbb(a,b,m,i) r1_a(e)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") += reuse_tmps_["37_bb_vo"]("f,m") * r1["a_Lv"]("X,e");

        // sigmar2_abb += +1.00 f_bb(f,m) r1_a(e)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") += f["bb_vo"]("f,m") * r1["a_Lv"]("X,e");
        sigmar2_bbb("X,e,f,m")  = -1.00 * reuse_tmps_["41_bb_vo"]("f,m") * r1["b_Lv"]("X,e");
        sigmar2_bbb("X,e,f,m") -= reuse_tmps_["34_bb_vo"]("f,m") * r1["b_Lv"]("X,e");

        // sigmar2_bbb += +1.00 P(e,f) f_bb(i,a) t2_bbbb(a,e,m,i) r1_b(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") += reuse_tmps_["34_bb_vo"]("e,m") * r1["b_Lv"]("X,f");
        sigmar2_bbb("X,e,f,m") -= reuse_tmps_["42_bb_vo"]("f,m") * r1["b_Lv"]("X,e");
        sigmar2_bbb("X,e,f,m") -= reuse_tmps_["38_bb_vo"]("f,m") * r1["b_Lv"]("X,e");
        sigmar2_bbb("X,e,f,m") += reuse_tmps_["37_bb_vo"]("f,m") * r1["b_Lv"]("X,e");

        // sigmar2_bbb += +0.50 P(e,f) <j,i||a,m>_abab t2_abab(a,e,j,i) r1_b(f) @ +0.50 P(e,f) <i,j||a,m>_abab t2_abab(a,e,i,j) r1_b(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") += reuse_tmps_["41_bb_vo"]("e,m") * r1["b_Lv"]("X,f");
        sigmar2_bbb("X,e,f,m") += f["bb_vo"]("f,m") * r1["b_Lv"]("X,e");
        sigmar2_bbb("X,e,f,m") += reuse_tmps_["33_bb_vo"]("f,m") * r1["b_Lv"]("X,e");

        // sigmar2_bbb += -0.50 P(e,f) <i,e||a,b>_abab t2_abab(a,b,i,m) r1_b(f) @ -0.50 P(e,f) <i,e||b,a>_abab t2_abab(b,a,i,m) r1_b(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") += reuse_tmps_["38_bb_vo"]("e,m") * r1["b_Lv"]("X,f");

        // sigmar2_bbb += -1.00 P(e,f) f_bb(e,m) r1_b(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") -= f["bb_vo"]("e,m") * r1["b_Lv"]("X,f");

        // sigmar2_bbb += -1.00 P(e,f) f_aa(i,a) t2_abab(a,e,i,m) r1_b(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") -= reuse_tmps_["33_bb_vo"]("e,m") * r1["b_Lv"]("X,f");

        // sigmar2_bbb += +0.50 P(e,f) <i,e||a,b>_bbbb t2_bbbb(a,b,m,i) r1_b(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") -= reuse_tmps_["37_bb_vo"]("e,m") * r1["b_Lv"]("X,f");

        // sigmar2_bbb += +0.50 P(e,f) <j,i||a,m>_bbbb t2_bbbb(a,e,j,i) r1_b(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") += reuse_tmps_["42_bb_vo"]("e,m") * r1["b_Lv"]("X,f");

        // sigmal1_a += -0.50 <j,a||b,c>_bbbb l2_bab(i,e,a) t2_bbbb(b,c,i,j)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") += reuse_tmps_["37_bb_vo"]("a,i") * l2["bab_Lovv"]("X,i,e,a");

        // sigmal1_a += +0.50 <j,a||b,c>_abab l2_bab(i,e,a) t2_abab(b,c,j,i) @ +0.50 <j,a||c,b>_abab l2_bab(i,e,a) t2_abab(c,b,j,i)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") -= reuse_tmps_["38_bb_vo"]("a,i") * l2["bab_Lovv"]("X,i,e,a");

        // sigmal1_a += +0.50 <k,j||i,b>_abab l2_aaa(i,a,e) t2_abab(a,b,k,j) @ +0.50 <j,k||i,b>_abab l2_aaa(i,a,e) t2_abab(a,b,j,k)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") -= reuse_tmps_["44_aa_vo"]("a,i") * l2["aaa_Lovv"]("X,i,a,e");

        // sigmal1_a += +0.50 <k,j||b,i>_aaaa l2_aaa(i,a,e) t2_aaaa(b,a,k,j)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") += reuse_tmps_["43_aa_vo"]("a,i") * l2["aaa_Lovv"]("X,i,a,e");

        // sigmal1_a += -1.00 f_bb(j,b) l2_bab(i,e,a) t2_bbbb(b,a,i,j)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") -= reuse_tmps_["34_bb_vo"]("a,i") * l2["bab_Lovv"]("X,i,e,a");

        // sigmal1_a += +0.50 <j,a||b,c>_aaaa l2_aaa(i,a,e) t2_aaaa(b,c,i,j)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") -= reuse_tmps_["39_aa_vo"]("a,i") * l2["aaa_Lovv"]("X,i,a,e");

        // sigmal1_a += -0.50 <k,j||b,i>_abab l2_bab(i,e,a) t2_abab(b,a,k,j) @ -0.50 <j,k||b,i>_abab l2_bab(i,e,a) t2_abab(b,a,j,k)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") -= reuse_tmps_["41_bb_vo"]("a,i") * l2["bab_Lovv"]("X,i,e,a");

        // sigmal1_a += -1.00 f_aa(a,i) l2_aaa(i,a,e)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") -= f["aa_vo"]("a,i") * l2["aaa_Lovv"]("X,i,a,e");

        // sigmal1_a += -1.00 f_bb(j,b) l2_aaa(i,a,e) t2_abab(a,b,i,j)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") -= reuse_tmps_["36_aa_vo"]("a,i") * l2["aaa_Lovv"]("X,i,a,e");

        // sigmal1_a += +1.00 f_aa(j,b) l2_aaa(i,a,e) t2_aaaa(b,a,i,j)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") += reuse_tmps_["35_aa_vo"]("a,i") * l2["aaa_Lovv"]("X,i,a,e");

        // sigmal1_a += +1.00 f_bb(a,i) l2_bab(i,e,a)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") += f["bb_vo"]("a,i") * l2["bab_Lovv"]("X,i,e,a");

        // sigmal1_a += -0.50 <k,j||b,i>_bbbb l2_bab(i,e,a) t2_bbbb(b,a,k,j)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") -= reuse_tmps_["42_bb_vo"]("a,i") * l2["bab_Lovv"]("X,i,e,a");

        // sigmal1_a += -0.50 <a,j||b,c>_abab l2_aaa(i,a,e) t2_abab(b,c,i,j) @ -0.50 <a,j||c,b>_abab l2_aaa(i,a,e) t2_abab(c,b,i,j)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") -= reuse_tmps_["40_aa_vo"]("a,i") * l2["aaa_Lovv"]("X,i,a,e");

        // sigmal1_a += +1.00 f_aa(j,b) l2_bab(i,e,a) t2_abab(b,a,j,i)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") += reuse_tmps_["33_bb_vo"]("a,i") * l2["bab_Lovv"]("X,i,e,a");

        // sigmal1_b += +1.00 f_aa(j,b) l2_aab(i,a,e) t2_aaaa(b,a,i,j)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") += reuse_tmps_["35_aa_vo"]("a,i") * l2["aab_Lovv"]("X,i,a,e");

        // sigmal1_b += +0.50 <k,j||i,b>_abab l2_aab(i,a,e) t2_abab(a,b,k,j) @ +0.50 <j,k||i,b>_abab l2_aab(i,a,e) t2_abab(a,b,j,k)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") -= reuse_tmps_["44_aa_vo"]("a,i") * l2["aab_Lovv"]("X,i,a,e");

        // sigmal1_b += +0.50 <k,j||b,i>_aaaa l2_aab(i,a,e) t2_aaaa(b,a,k,j)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") += reuse_tmps_["43_aa_vo"]("a,i") * l2["aab_Lovv"]("X,i,a,e");

        // sigmal1_b += +1.00 f_bb(j,b) l2_bbb(i,a,e) t2_bbbb(b,a,i,j)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") += reuse_tmps_["34_bb_vo"]("a,i") * l2["bbb_Lovv"]("X,i,a,e");

        // sigmal1_b += +0.50 <j,a||b,c>_bbbb l2_bbb(i,a,e) t2_bbbb(b,c,i,j)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") -= reuse_tmps_["37_bb_vo"]("a,i") * l2["bbb_Lovv"]("X,i,a,e");

        // sigmal1_b += -1.00 f_bb(a,i) l2_bbb(i,a,e)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") -= f["bb_vo"]("a,i") * l2["bbb_Lovv"]("X,i,a,e");

        // sigmal1_b += -1.00 f_aa(j,b) l2_bbb(i,a,e) t2_abab(b,a,j,i)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") -= reuse_tmps_["33_bb_vo"]("a,i") * l2["bbb_Lovv"]("X,i,a,e");

        // sigmal1_b += -0.50 <a,j||b,c>_abab l2_aab(i,a,e) t2_abab(b,c,i,j) @ -0.50 <a,j||c,b>_abab l2_aab(i,a,e) t2_abab(c,b,i,j)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") -= reuse_tmps_["40_aa_vo"]("a,i") * l2["aab_Lovv"]("X,i,a,e");

        // sigmal1_b += -1.00 f_aa(a,i) l2_aab(i,a,e)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") -= f["aa_vo"]("a,i") * l2["aab_Lovv"]("X,i,a,e");

        // sigmal1_b += +0.50 <k,j||b,i>_bbbb l2_bbb(i,a,e) t2_bbbb(b,a,k,j)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") += reuse_tmps_["42_bb_vo"]("a,i") * l2["bbb_Lovv"]("X,i,a,e");

        // sigmal1_b += -1.00 f_bb(j,b) l2_aab(i,a,e) t2_abab(a,b,i,j)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") -= reuse_tmps_["36_aa_vo"]("a,i") * l2["aab_Lovv"]("X,i,a,e");

        // sigmal1_b += +0.50 <k,j||b,i>_abab l2_bbb(i,a,e) t2_abab(b,a,k,j) @ +0.50 <j,k||b,i>_abab l2_bbb(i,a,e) t2_abab(b,a,j,k)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") += reuse_tmps_["41_bb_vo"]("a,i") * l2["bbb_Lovv"]("X,i,a,e");

        // sigmal1_b += -0.50 <j,a||b,c>_abab l2_bbb(i,a,e) t2_abab(b,c,j,i) @ -0.50 <j,a||c,b>_abab l2_bbb(i,a,e) t2_abab(c,b,j,i)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") += reuse_tmps_["38_bb_vo"]("a,i") * l2["bbb_Lovv"]("X,i,a,e");

        // sigmal1_b += +0.50 <j,a||b,c>_aaaa l2_aab(i,a,e) t2_aaaa(b,c,i,j)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") -= reuse_tmps_["39_aa_vo"]("a,i") * l2["aab_Lovv"]("X,i,a,e");
        sigmal2_aaa("X,e,f,m")  = f["aa_ov"]("m,f") * l1["a_Lv"]("X,e");

        // sigmal2_aaa += -1.00 P(e,f) f_aa(m,e) l1_a(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aaa("X,e,f,m") -= f["aa_ov"]("m,e") * l1["a_Lv"]("X,f");

        // sigmal2_aba = -1.00 f_aa(m,e) l1_b(f)  // flops: o1v2L1 = o1v2L1 | mem: o1v2L1 = o1v2L1
        sigmal2_aba("X,e,f,m")  = -1.00 * f["aa_ov"]("m,e") * l1["b_Lv"]("X,f");

        // sigmal2_abb = +1.00 f_bb(m,f) l1_a(e)  // flops: o1v2L1 = o1v2L1 | mem: o1v2L1 = o1v2L1
        sigmal2_abb("X,e,f,m")  = f["bb_ov"]("m,f") * l1["a_Lv"]("X,e");
        sigmal2_bbb("X,e,f,m")  = f["bb_ov"]("m,f") * l1["b_Lv"]("X,e");

        // sigmal2_bbb += -1.00 P(e,f) f_bb(m,e) l1_b(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_bbb("X,e,f,m") -= f["bb_ov"]("m,e") * l1["b_Lv"]("X,f");

        // sigmar2_aaa += -0.50 <j,i||a,b>_abab t2_abab(a,b,m,i) r2_aaa(e,f,j) @ -0.50 <j,i||b,a>_abab t2_abab(b,a,m,i) r2_aaa(e,f,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") -= reuse_tmps_["24_aa_oo"]("m,j") * r2["aaa_Lvvo"]("X,e,f,j");

        // sigmar2_aaa += -0.50 <j,i||a,b>_aaaa t2_aaaa(a,b,m,i) r2_aaa(e,f,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") -= reuse_tmps_["27_aa_oo"]("m,j") * r2["aaa_Lvvo"]("X,e,f,j");

        // sigmar2_aaa += -1.00 f_aa(i,m) r2_aaa(e,f,i)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") -= f["aa_oo"]("i,m") * r2["aaa_Lvvo"]("X,e,f,i");

        // sigmar2_aba += -1.00 f_aa(i,m) r2_aba(e,f,i)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= f["aa_oo"]("i,m") * r2["aba_Lvvo"]("X,e,f,i");

        // sigmar2_aba += -0.50 <j,i||a,b>_aaaa t2_aaaa(a,b,m,i) r2_aba(e,f,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= reuse_tmps_["27_aa_oo"]("m,j") * r2["aba_Lvvo"]("X,e,f,j");

        // sigmar2_aba += -0.50 <j,i||a,b>_abab t2_abab(a,b,m,i) r2_aba(e,f,j) @ -0.50 <j,i||b,a>_abab t2_abab(b,a,m,i) r2_aba(e,f,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= reuse_tmps_["24_aa_oo"]("m,j") * r2["aba_Lvvo"]("X,e,f,j");

        // sigmar2_abb += -1.00 f_bb(i,m) r2_abb(e,f,i)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") -= f["bb_oo"]("i,m") * r2["abb_Lvvo"]("X,e,f,i");

        // sigmar2_abb += -0.50 <j,i||a,b>_bbbb t2_bbbb(a,b,m,i) r2_abb(e,f,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") -= reuse_tmps_["28_bb_oo"]("m,j") * r2["abb_Lvvo"]("X,e,f,j");

        // sigmar2_abb += -0.50 <i,j||a,b>_abab t2_abab(a,b,i,m) r2_abb(e,f,j) @ -0.50 <i,j||b,a>_abab t2_abab(b,a,i,m) r2_abb(e,f,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") -= reuse_tmps_["23_bb_oo"]("m,j") * r2["abb_Lvvo"]("X,e,f,j");

        // sigmar2_bbb += -1.00 f_bb(i,m) r2_bbb(e,f,i)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") -= f["bb_oo"]("i,m") * r2["bbb_Lvvo"]("X,e,f,i");

        // sigmar2_bbb += -0.50 <i,j||a,b>_abab t2_abab(a,b,i,m) r2_bbb(e,f,j) @ -0.50 <i,j||b,a>_abab t2_abab(b,a,i,m) r2_bbb(e,f,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") -= reuse_tmps_["23_bb_oo"]("m,j") * r2["bbb_Lvvo"]("X,e,f,j");

        // sigmar2_bbb += -0.50 <j,i||a,b>_bbbb t2_bbbb(a,b,m,i) r2_bbb(e,f,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") -= reuse_tmps_["28_bb_oo"]("m,j") * r2["bbb_Lvvo"]("X,e,f,j");

        // sigmal2_aaa += -1.00 f_aa(m,i) l2_aaa(i,e,f)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aaa("X,e,f,m") -= f["aa_oo"]("m,i") * l2["aaa_Lovv"]("X,i,e,f");

        // sigmal2_aaa += -0.50 <m,j||a,b>_aaaa l2_aaa(i,e,f) t2_aaaa(a,b,i,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aaa("X,e,f,m") -= reuse_tmps_["27_aa_oo"]("i,m") * l2["aaa_Lovv"]("X,i,e,f");

        // sigmal2_aaa += -0.50 <m,j||a,b>_abab l2_aaa(i,e,f) t2_abab(a,b,i,j) @ -0.50 <m,j||b,a>_abab l2_aaa(i,e,f) t2_abab(b,a,i,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aaa("X,e,f,m") -= reuse_tmps_["24_aa_oo"]("i,m") * l2["aaa_Lovv"]("X,i,e,f");

        // sigmal2_aba += -0.50 <m,j||a,b>_abab l2_aab(i,e,f) t2_abab(a,b,i,j) @ -0.50 <m,j||b,a>_abab l2_aab(i,e,f) t2_abab(b,a,i,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") -= reuse_tmps_["24_aa_oo"]("i,m") * l2["aab_Lovv"]("X,i,e,f");

        // sigmal2_aba += -1.00 f_aa(m,i) l2_aab(i,e,f)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") -= f["aa_oo"]("m,i") * l2["aab_Lovv"]("X,i,e,f");

        // sigmal2_aba += -0.50 <m,j||a,b>_aaaa l2_aab(i,e,f) t2_aaaa(a,b,i,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") -= reuse_tmps_["27_aa_oo"]("i,m") * l2["aab_Lovv"]("X,i,e,f");

        // sigmal2_abb += -0.50 <j,m||a,b>_abab l2_bab(i,e,f) t2_abab(a,b,j,i) @ -0.50 <j,m||b,a>_abab l2_bab(i,e,f) t2_abab(b,a,j,i)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") -= reuse_tmps_["23_bb_oo"]("i,m") * l2["bab_Lovv"]("X,i,e,f");

        // sigmal2_abb += -0.50 <m,j||a,b>_bbbb l2_bab(i,e,f) t2_bbbb(a,b,i,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") -= reuse_tmps_["28_bb_oo"]("i,m") * l2["bab_Lovv"]("X,i,e,f");

        // sigmal2_abb += -1.00 f_bb(m,i) l2_bab(i,e,f)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") -= f["bb_oo"]("m,i") * l2["bab_Lovv"]("X,i,e,f");

        // sigmal2_bbb += -1.00 f_bb(m,i) l2_bbb(i,e,f)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_bbb("X,e,f,m") -= f["bb_oo"]("m,i") * l2["bbb_Lovv"]("X,i,e,f");

        // sigmal2_bbb += -0.50 <j,m||a,b>_abab l2_bbb(i,e,f) t2_abab(a,b,j,i) @ -0.50 <j,m||b,a>_abab l2_bbb(i,e,f) t2_abab(b,a,j,i)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_bbb("X,e,f,m") -= reuse_tmps_["23_bb_oo"]("i,m") * l2["bbb_Lovv"]("X,i,e,f");

        // sigmal2_bbb += -0.50 <m,j||a,b>_bbbb l2_bbb(i,e,f) t2_bbbb(a,b,i,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_bbb("X,e,f,m") -= reuse_tmps_["28_bb_oo"]("i,m") * l2["bbb_Lovv"]("X,i,e,f");

        // sigmar1_a += -0.50 <i,e||a,b>_aaaa r2_aaa(a,b,i)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmar1_a("X,e") += 0.50 * eri["aaaa_vovv"]("e,i,a,b") * r2["aaa_Lvvo"]("X,a,b,i");

        // sigmar1_a += +0.50 <e,i||a,b>_abab r2_abb(a,b,i) @ +0.50 <e,i||b,a>_abab r2_abb(b,a,i)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmar1_a("X,e") += eri["abab_vovv"]("e,i,a,b") * r2["abb_Lvvo"]("X,a,b,i");

        // sigmar1_b += -0.50 <i,e||a,b>_abab r2_aba(a,b,i) @ -0.50 <i,e||b,a>_abab r2_aba(b,a,i)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmar1_b("X,e") += eri["baab_vovv"]("e,i,a,b") * r2["aba_Lvvo"]("X,a,b,i");

        // sigmar1_b += -0.50 <i,e||a,b>_bbbb r2_bbb(a,b,i)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmar1_b("X,e") += 0.50 * eri["bbbb_vovv"]("e,i,a,b") * r2["bbb_Lvvo"]("X,a,b,i");

        // sigmar2_aaa += +1.00 <e,f||a,m>_aaaa r1_a(a)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") += eri["aaaa_vvvo"]("e,f,a,m") * r1["a_Lv"]("X,a");

        // sigmar2_aaa += +0.50 <j,i||a,m>_aaaa t2_aaaa(e,f,j,i) r1_a(a)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") += 0.50 * reuse_tmps_["32_aaaa_vvvo"]("f,e,a,m") * r1["a_Lv"]("X,a");

        // sigmar2_aba += +1.00 <e,i||a,b>_abab t2_abab(a,f,m,i) r1_b(b)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") += reuse_tmps_["11_abba_vvvo"]("e,f,b,m") * r1["b_Lv"]("X,b");

        // sigmar2_aba += -1.00 <i,f||a,b>_bbbb t2_abab(e,a,m,i) r1_b(b)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") += reuse_tmps_["20_abba_vvvo"]("e,b,f,m") * r1["b_Lv"]("X,b");

        // sigmar2_aba += -0.50 <j,i||a,b>_bbbb t2_bbbb(a,f,j,i) r2_aba(e,b,m)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= reuse_tmps_["25_bb_vv"]("f,b") * r2["aba_Lvvo"]("X,e,b,m");

        // sigmar2_aba += +1.00 f_bb(f,a) r2_aba(e,a,m)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") += f["bb_vv"]("f,a") * r2["aba_Lvvo"]("X,e,a,m");

        // sigmar2_aba += -0.50 <j,i||m,a>_abab t2_abab(e,f,j,i) r1_b(a) @ -0.50 <i,j||m,a>_abab t2_abab(e,f,i,j) r1_b(a)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") += reuse_tmps_["29_abba_vvvo"]("e,f,a,m") * r1["b_Lv"]("X,a");

        // sigmar2_aba += -0.50 <j,i||a,b>_aaaa t2_aaaa(a,e,j,i) r2_aba(b,f,m)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= reuse_tmps_["26_aa_vv"]("e,b") * r2["aba_Lvvo"]("X,b,f,m");

        // sigmar2_aba += -1.00 <e,f||m,a>_abab r1_b(a)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") += eri["abba_vvvo"]("e,f,a,m") * r1["b_Lv"]("X,a");

        // sigmar2_aba += +1.00 <i,f||a,b>_abab t2_aaaa(a,e,m,i) r1_b(b)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= reuse_tmps_["17_abba_vvvo"]("e,b,f,m") * r1["b_Lv"]("X,b");

        // sigmar2_aba += +1.00 f_aa(e,a) r2_aba(a,f,m)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") += f["aa_vv"]("e,a") * r2["aba_Lvvo"]("X,a,f,m");

        // sigmar2_aba += -0.50 <j,i||b,a>_abab t2_abab(e,a,j,i) r2_aba(b,f,m) @ -0.50 <i,j||b,a>_abab t2_abab(e,a,i,j) r2_aba(b,f,m)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= reuse_tmps_["21_aa_vv"]("e,b") * r2["aba_Lvvo"]("X,b,f,m");

        // sigmar2_aba += -0.50 <j,i||a,b>_abab t2_abab(a,f,j,i) r2_aba(e,b,m) @ -0.50 <i,j||a,b>_abab t2_abab(a,f,i,j) r2_aba(e,b,m)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= reuse_tmps_["22_bb_vv"]("f,b") * r2["aba_Lvvo"]("X,e,b,m");

        // sigmar2_abb += -1.00 <i,f||b,a>_abab t2_abab(e,a,i,m) r1_a(b)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") += reuse_tmps_["16_aabb_vvvo"]("e,b,f,m") * r1["a_Lv"]("X,b");

        // sigmar2_abb += +1.00 <i,e||a,b>_aaaa t2_abab(a,f,i,m) r1_a(b)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") -= reuse_tmps_["19_aabb_vvvo"]("b,e,f,m") * r1["a_Lv"]("X,b");

        // sigmar2_abb += -0.50 <j,i||a,b>_bbbb t2_bbbb(a,f,j,i) r2_abb(e,b,m)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") -= reuse_tmps_["25_bb_vv"]("f,b") * r2["abb_Lvvo"]("X,e,b,m");

        // sigmar2_abb += +0.50 <j,i||a,m>_abab t2_abab(e,f,j,i) r1_a(a) @ +0.50 <i,j||a,m>_abab t2_abab(e,f,i,j) r1_a(a)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") += reuse_tmps_["31_aabb_vvvo"]("e,a,f,m") * r1["a_Lv"]("X,a");

        // sigmar2_abb += +1.00 f_aa(e,a) r2_abb(a,f,m)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") += f["aa_vv"]("e,a") * r2["abb_Lvvo"]("X,a,f,m");

        // sigmar2_abb += -0.50 <j,i||a,b>_aaaa t2_aaaa(a,e,j,i) r2_abb(b,f,m)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") -= reuse_tmps_["26_aa_vv"]("e,b") * r2["abb_Lvvo"]("X,b,f,m");

        // sigmar2_abb += -1.00 <e,i||b,a>_abab t2_bbbb(a,f,m,i) r1_a(b)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") -= reuse_tmps_["18_aabb_vvvo"]("b,e,f,m") * r1["a_Lv"]("X,b");

        // sigmar2_abb += +1.00 <e,f||a,m>_abab r1_a(a)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") += eri["abab_vvvo"]("e,f,a,m") * r1["a_Lv"]("X,a");

        // sigmar2_abb += +1.00 f_bb(f,a) r2_abb(e,a,m)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") += f["bb_vv"]("f,a") * r2["abb_Lvvo"]("X,e,a,m");

        // sigmar2_abb += -0.50 <j,i||b,a>_abab t2_abab(e,a,j,i) r2_abb(b,f,m) @ -0.50 <i,j||b,a>_abab t2_abab(e,a,i,j) r2_abb(b,f,m)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") -= reuse_tmps_["21_aa_vv"]("e,b") * r2["abb_Lvvo"]("X,b,f,m");

        // sigmar2_abb += -0.50 <j,i||a,b>_abab t2_abab(a,f,j,i) r2_abb(e,b,m) @ -0.50 <i,j||a,b>_abab t2_abab(a,f,i,j) r2_abb(e,b,m)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") -= reuse_tmps_["22_bb_vv"]("f,b") * r2["abb_Lvvo"]("X,e,b,m");

        // sigmar2_bbb += +0.50 <j,i||a,m>_bbbb t2_bbbb(e,f,j,i) r1_b(a)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") += 0.50 * reuse_tmps_["30_bbbb_vvvo"]("f,e,a,m") * r1["b_Lv"]("X,a");

        // sigmar2_bbb += +1.00 <e,f||a,m>_bbbb r1_b(a)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") += eri["bbbb_vvvo"]("e,f,a,m") * r1["b_Lv"]("X,a");

        // sigmal1_a += -1.00 <j,b||e,c>_abab l2_bab(i,a,b) t2_abab(a,c,j,i)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") += reuse_tmps_["16_aabb_vvvo"]("a,e,b,i") * l2["bab_Lovv"]("X,i,a,b");

        // sigmal1_a += +0.25 <k,j||e,i>_aaaa l2_aaa(i,b,a) t2_aaaa(b,a,k,j)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") += 0.25 * reuse_tmps_["32_aaaa_vvvo"]("a,b,e,i") * l2["aaa_Lovv"]("X,i,b,a");

        // sigmal1_a += -1.00 <j,b||c,e>_aaaa l2_aaa(i,b,a) t2_aaaa(c,a,i,j)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") += reuse_tmps_["14_aaaa_vvvo"]("a,e,b,i") * l2["aaa_Lovv"]("X,i,b,a");

        // sigmal1_a += +0.25 <k,j||e,i>_abab l2_bab(i,b,a) t2_abab(b,a,k,j) @ +0.25 <j,k||e,i>_abab l2_bab(i,b,a) t2_abab(b,a,j,k) @ +0.25 <k,j||e,i>_abab l2_bab(i,a,b) t2_abab(a,b,k,j) @ +0.25 <j,k||e,i>_abab l2_bab(i,a,b) t2_abab(a,b,j,k)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") += reuse_tmps_["31_aabb_vvvo"]("b,e,a,i") * l2["bab_Lovv"]("X,i,b,a");

        // sigmal1_a += +1.00 <b,j||e,c>_abab l2_aaa(i,b,a) t2_abab(a,c,i,j)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") += reuse_tmps_["15_aaaa_vvvo"]("a,e,b,i") * l2["aaa_Lovv"]("X,i,b,a");

        // sigmal1_a += +0.50 <b,a||e,i>_aaaa l2_aaa(i,b,a)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") += 0.50 * eri["aaaa_vvvo"]("b,a,e,i") * l2["aaa_Lovv"]("X,i,b,a");

        // sigmal1_a += +1.00 <j,b||c,e>_aaaa l2_bab(i,b,a) t2_abab(c,a,j,i)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") -= reuse_tmps_["19_aabb_vvvo"]("e,b,a,i") * l2["bab_Lovv"]("X,i,b,a");

        // sigmal1_a += -1.00 <b,j||e,c>_abab l2_bab(i,b,a) t2_bbbb(c,a,i,j)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") -= reuse_tmps_["18_aabb_vvvo"]("e,b,a,i") * l2["bab_Lovv"]("X,i,b,a");

        // sigmal1_a += +0.50 <b,a||e,i>_abab l2_bab(i,b,a) @ +0.50 <a,b||e,i>_abab l2_bab(i,a,b)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") += eri["abab_vvvo"]("b,a,e,i") * l2["bab_Lovv"]("X,i,b,a");

        // sigmal1_b += -0.25 <k,j||i,e>_abab l2_aab(i,b,a) t2_abab(b,a,k,j) @ -0.25 <j,k||i,e>_abab l2_aab(i,b,a) t2_abab(b,a,j,k) @ -0.25 <k,j||i,e>_abab l2_aab(i,a,b) t2_abab(a,b,k,j) @ -0.25 <j,k||i,e>_abab l2_aab(i,a,b) t2_abab(a,b,j,k)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") += reuse_tmps_["29_abba_vvvo"]("b,a,e,i") * l2["aab_Lovv"]("X,i,b,a");

        // sigmal1_b += -1.00 <j,b||c,e>_bbbb l2_bbb(i,b,a) t2_bbbb(c,a,i,j)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") += reuse_tmps_["13_bbbb_vvvo"]("a,e,b,i") * l2["bbb_Lovv"]("X,i,b,a");

        // sigmal1_b += +1.00 <b,j||c,e>_abab l2_aab(i,b,a) t2_abab(c,a,i,j)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") += reuse_tmps_["11_abba_vvvo"]("b,a,e,i") * l2["aab_Lovv"]("X,i,b,a");

        // sigmal1_b += -1.00 <j,b||c,e>_bbbb l2_aab(i,a,b) t2_abab(a,c,i,j)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") += reuse_tmps_["20_abba_vvvo"]("a,e,b,i") * l2["aab_Lovv"]("X,i,a,b");

        // sigmal1_b += +0.25 <k,j||e,i>_bbbb l2_bbb(i,b,a) t2_bbbb(b,a,k,j)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") += 0.25 * reuse_tmps_["30_bbbb_vvvo"]("a,b,e,i") * l2["bbb_Lovv"]("X,i,b,a");

        // sigmal1_b += +1.00 <j,b||c,e>_abab l2_bbb(i,b,a) t2_abab(c,a,j,i)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") -= reuse_tmps_["12_bbbb_vvvo"]("a,e,b,i") * l2["bbb_Lovv"]("X,i,b,a");

        // sigmal1_b += -0.50 <b,a||i,e>_abab l2_aab(i,b,a) @ -0.50 <a,b||i,e>_abab l2_aab(i,a,b)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") += eri["abba_vvvo"]("b,a,e,i") * l2["aab_Lovv"]("X,i,b,a");

        // sigmal1_b += +1.00 <j,b||c,e>_abab l2_aab(i,a,b) t2_aaaa(c,a,i,j)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") -= reuse_tmps_["17_abba_vvvo"]("a,e,b,i") * l2["aab_Lovv"]("X,i,a,b");

        // sigmal1_b += +0.50 <b,a||e,i>_bbbb l2_bbb(i,b,a)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") += 0.50 * eri["bbbb_vvvo"]("b,a,e,i") * l2["bbb_Lovv"]("X,i,b,a");

        // sigmal2_aaa += -1.00 <m,a||e,f>_aaaa l1_a(a)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aaa("X,e,f,m") += eri["aaaa_vovv"]("a,m,e,f") * l1["a_Lv"]("X,a");

        // sigmal2_aba += +1.00 f_bb(a,f) l2_aab(m,e,a)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") += f["bb_vv"]("a,f") * l2["aab_Lovv"]("X,m,e,a");

        // sigmal2_aba += -1.00 <m,a||e,f>_abab l1_b(a)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") += eri["baab_vovv"]("a,m,e,f") * l1["b_Lv"]("X,a");

        // sigmal2_aba += +1.00 f_aa(a,e) l2_aab(m,a,f)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") += f["aa_vv"]("a,e") * l2["aab_Lovv"]("X,m,a,f");

        // sigmal2_aba += -0.50 <j,i||b,f>_abab l2_aab(m,e,a) t2_abab(b,a,j,i) @ -0.50 <i,j||b,f>_abab l2_aab(m,e,a) t2_abab(b,a,i,j)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") -= reuse_tmps_["22_bb_vv"]("a,f") * l2["aab_Lovv"]("X,m,e,a");

        // sigmal2_aba += -0.50 <j,i||e,b>_abab l2_aab(m,a,f) t2_abab(a,b,j,i) @ -0.50 <i,j||e,b>_abab l2_aab(m,a,f) t2_abab(a,b,i,j)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") -= reuse_tmps_["21_aa_vv"]("a,e") * l2["aab_Lovv"]("X,m,a,f");

        // sigmal2_aba += -0.50 <j,i||b,e>_aaaa l2_aab(m,a,f) t2_aaaa(b,a,j,i)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") -= reuse_tmps_["26_aa_vv"]("a,e") * l2["aab_Lovv"]("X,m,a,f");

        // sigmal2_aba += -0.50 <j,i||b,f>_bbbb l2_aab(m,e,a) t2_bbbb(b,a,j,i)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") -= reuse_tmps_["25_bb_vv"]("a,f") * l2["aab_Lovv"]("X,m,e,a");

        // sigmal2_abb += +1.00 f_bb(a,f) l2_bab(m,e,a)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") += f["bb_vv"]("a,f") * l2["bab_Lovv"]("X,m,e,a");

        // sigmal2_abb += +1.00 f_aa(a,e) l2_bab(m,a,f)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") += f["aa_vv"]("a,e") * l2["bab_Lovv"]("X,m,a,f");

        // sigmal2_abb += +1.00 <a,m||e,f>_abab l1_a(a)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") += eri["abab_vovv"]("a,m,e,f") * l1["a_Lv"]("X,a");

        // sigmal2_abb += -0.50 <j,i||b,e>_aaaa l2_bab(m,a,f) t2_aaaa(b,a,j,i)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") -= reuse_tmps_["26_aa_vv"]("a,e") * l2["bab_Lovv"]("X,m,a,f");

        // sigmal2_abb += -0.50 <j,i||e,b>_abab l2_bab(m,a,f) t2_abab(a,b,j,i) @ -0.50 <i,j||e,b>_abab l2_bab(m,a,f) t2_abab(a,b,i,j)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") -= reuse_tmps_["21_aa_vv"]("a,e") * l2["bab_Lovv"]("X,m,a,f");

        // sigmal2_abb += -0.50 <j,i||b,f>_abab l2_bab(m,e,a) t2_abab(b,a,j,i) @ -0.50 <i,j||b,f>_abab l2_bab(m,e,a) t2_abab(b,a,i,j)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") -= reuse_tmps_["22_bb_vv"]("a,f") * l2["bab_Lovv"]("X,m,e,a");

        // sigmal2_abb += -0.50 <j,i||b,f>_bbbb l2_bab(m,e,a) t2_bbbb(b,a,j,i)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") -= reuse_tmps_["25_bb_vv"]("a,f") * l2["bab_Lovv"]("X,m,e,a");

        // sigmal2_bbb += -1.00 <m,a||e,f>_bbbb l1_b(a)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_bbb("X,e,f,m") += eri["bbbb_vovv"]("a,m,e,f") * l1["b_Lv"]("X,a");

        // sigmar2_aba += -1.00 <i,f||m,a>_abab r2_aba(e,a,i)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= eri["baba_vovo"]("f,i,a,m") * r2["aba_Lvvo"]("X,e,a,i");

        // sigmar2_aba += +1.00 <e,i||m,a>_abab r2_bbb(a,f,i)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= eri["abba_vovo"]("e,i,a,m") * r2["bbb_Lvvo"]("X,a,f,i");

        // sigmar2_aba += +1.00 <i,e||a,m>_aaaa r2_aba(a,f,i)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= eri["aaaa_vovo"]("e,i,a,m") * r2["aba_Lvvo"]("X,a,f,i");

        // sigmar2_aba += -1.00 <j,i||a,b>_bbbb t2_abab(e,a,m,i) r2_bbb(b,f,j)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= reuse_tmps_["8_abab_vvoo"]("e,b,m,j") * r2["bbb_Lvvo"]("X,b,f,j");

        // sigmar2_aba += +1.00 <j,i||a,b>_aaaa t2_aaaa(a,e,m,i) r2_aba(b,f,j)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") += reuse_tmps_["4_aaaa_vvoo"]("e,b,m,j") * r2["aba_Lvvo"]("X,b,f,j");

        // sigmar2_aba += -1.00 <i,j||a,b>_abab t2_aaaa(a,e,m,i) r2_bbb(b,f,j)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= reuse_tmps_["7_abab_vvoo"]("e,b,m,j") * r2["bbb_Lvvo"]("X,b,f,j");

        // sigmar2_aba += +1.00 <j,i||a,b>_abab t2_abab(a,f,m,i) r2_aba(e,b,j)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") += reuse_tmps_["10_bbaa_vvoo"]("f,b,m,j") * r2["aba_Lvvo"]("X,e,b,j");

        // sigmar2_aba += +1.00 <j,i||b,a>_abab t2_abab(e,a,m,i) r2_aba(b,f,j)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") += reuse_tmps_["2_aaaa_vvoo"]("e,b,m,j") * r2["aba_Lvvo"]("X,b,f,j");

        // sigmar2_abb += +1.00 <i,f||a,m>_bbbb r2_abb(e,a,i)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") -= eri["bbbb_vovo"]("f,i,a,m") * r2["abb_Lvvo"]("X,e,a,i");

        // sigmar2_abb += +1.00 <j,i||a,b>_bbbb t2_bbbb(a,f,m,i) r2_abb(e,b,j)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") += reuse_tmps_["3_bbbb_vvoo"]("f,b,m,j") * r2["abb_Lvvo"]("X,e,b,j");

        // sigmar2_abb += -1.00 <i,f||a,m>_abab r2_aaa(a,e,i)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") += eri["baab_vovo"]("f,i,a,m") * r2["aaa_Lvvo"]("X,a,e,i");

        // sigmar2_abb += +1.00 <j,i||b,a>_abab t2_bbbb(a,f,m,i) r2_aaa(b,e,j)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") += reuse_tmps_["5_abab_vvoo"]("b,f,j,m") * r2["aaa_Lvvo"]("X,b,e,j");

        // sigmar2_abb += -1.00 <e,i||a,m>_abab r2_abb(a,f,i)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") -= eri["abab_vovo"]("e,i,a,m") * r2["abb_Lvvo"]("X,a,f,i");

        // sigmar2_abb += +1.00 <j,i||a,b>_aaaa t2_abab(a,f,i,m) r2_aaa(b,e,j)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") += reuse_tmps_["6_abab_vvoo"]("b,f,j,m") * r2["aaa_Lvvo"]("X,b,e,j");

        // sigmar2_abb += +1.00 <i,j||a,b>_abab t2_abab(a,f,i,m) r2_abb(e,b,j)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") += reuse_tmps_["1_bbbb_vvoo"]("f,b,m,j") * r2["abb_Lvvo"]("X,e,b,j");

        // sigmar2_abb += +1.00 <i,j||b,a>_abab t2_abab(e,a,i,m) r2_abb(b,f,j)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") += reuse_tmps_["9_aabb_vvoo"]("e,b,m,j") * r2["abb_Lvvo"]("X,b,f,j");

        // sigmal2_aba += -1.00 <m,a||i,f>_abab l2_aab(i,e,a)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") -= eri["baba_vovo"]("a,m,f,i") * l2["aab_Lovv"]("X,i,e,a");

        // sigmal2_aba += +1.00 <m,a||e,i>_abab l2_bbb(i,a,f)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") -= eri["baab_vovo"]("a,m,e,i") * l2["bbb_Lovv"]("X,i,a,f");

        // sigmal2_aba += +1.00 <m,j||b,f>_abab l2_aab(i,e,a) t2_abab(b,a,i,j)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") += reuse_tmps_["10_bbaa_vvoo"]("a,f,i,m") * l2["aab_Lovv"]("X,i,e,a");

        // sigmal2_aba += +1.00 <m,j||e,b>_abab l2_aab(i,a,f) t2_abab(a,b,i,j)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") += reuse_tmps_["2_aaaa_vvoo"]("a,e,i,m") * l2["aab_Lovv"]("X,i,a,f");

        // sigmal2_aba += +1.00 <m,j||b,e>_aaaa l2_aab(i,a,f) t2_aaaa(b,a,i,j)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") += reuse_tmps_["4_aaaa_vvoo"]("a,e,i,m") * l2["aab_Lovv"]("X,i,a,f");

        // sigmal2_aba += +1.00 <m,a||e,i>_aaaa l2_aab(i,a,f)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") -= eri["aaaa_vovo"]("a,m,e,i") * l2["aab_Lovv"]("X,i,a,f");

        // sigmal2_aba += -1.00 <m,j||e,b>_abab l2_bbb(i,a,f) t2_bbbb(b,a,i,j)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") -= reuse_tmps_["5_abab_vvoo"]("e,a,m,i") * l2["bbb_Lovv"]("X,i,a,f");

        // sigmal2_aba += -1.00 <m,j||b,e>_aaaa l2_bbb(i,a,f) t2_abab(b,a,j,i)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") -= reuse_tmps_["6_abab_vvoo"]("e,a,m,i") * l2["bbb_Lovv"]("X,i,a,f");

        // sigmal2_abb += +1.00 <m,a||f,i>_bbbb l2_bab(i,e,a)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") -= eri["bbbb_vovo"]("a,m,f,i") * l2["bab_Lovv"]("X,i,e,a");

        // sigmal2_abb += -1.00 <a,m||e,i>_abab l2_bab(i,a,f)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") -= eri["abab_vovo"]("a,m,e,i") * l2["bab_Lovv"]("X,i,a,f");

        // sigmal2_abb += +1.00 <j,m||e,b>_abab l2_bab(i,a,f) t2_abab(a,b,j,i)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") += reuse_tmps_["9_aabb_vvoo"]("a,e,i,m") * l2["bab_Lovv"]("X,i,a,f");

        // sigmal2_abb += +1.00 <j,m||b,f>_abab l2_bab(i,e,a) t2_abab(b,a,j,i)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") += reuse_tmps_["1_bbbb_vvoo"]("a,f,i,m") * l2["bab_Lovv"]("X,i,e,a");

        // sigmal2_abb += +1.00 <j,m||b,f>_abab l2_aaa(i,a,e) t2_aaaa(b,a,i,j)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") += reuse_tmps_["7_abab_vvoo"]("a,f,i,m") * l2["aaa_Lovv"]("X,i,a,e");

        // sigmal2_abb += +1.00 <m,j||b,f>_bbbb l2_bab(i,e,a) t2_bbbb(b,a,i,j)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") += reuse_tmps_["3_bbbb_vvoo"]("a,f,i,m") * l2["bab_Lovv"]("X,i,e,a");

        // sigmal2_abb += -1.00 <a,m||i,f>_abab l2_aaa(i,a,e)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") += eri["abba_vovo"]("a,m,f,i") * l2["aaa_Lovv"]("X,i,a,e");

        // sigmal2_abb += +1.00 <m,j||b,f>_bbbb l2_aaa(i,a,e) t2_abab(a,b,i,j)  // flops: o1v2L1 += o2v3L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") += reuse_tmps_["8_abab_vvoo"]("a,f,i,m") * l2["aaa_Lovv"]("X,i,a,e");

        // sigmar2_aaa += +0.50 <e,f||a,b>_aaaa r2_aaa(a,b,m)  // flops: o1v2L1 += o1v4L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") += 0.50 * eri["aaaa_vvvv"]("e,f,a,b") * r2["aaa_Lvvo"]("X,a,b,m");

        // sigmar2_aba += +0.50 <e,f||a,b>_abab r2_aba(a,b,m) @ +0.50 <e,f||b,a>_abab r2_aba(b,a,m)  // flops: o1v2L1 += o1v4L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") += eri["abab_vvvv"]("e,f,a,b") * r2["aba_Lvvo"]("X,a,b,m");

        // sigmar2_abb += +0.50 <e,f||a,b>_abab r2_abb(a,b,m) @ +0.50 <e,f||b,a>_abab r2_abb(b,a,m)  // flops: o1v2L1 += o1v4L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") += eri["abab_vvvv"]("e,f,a,b") * r2["abb_Lvvo"]("X,a,b,m");

        // sigmar2_bbb += +0.50 <e,f||a,b>_bbbb r2_bbb(a,b,m)  // flops: o1v2L1 += o1v4L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") += 0.50 * eri["bbbb_vvvv"]("e,f,a,b") * r2["bbb_Lvvo"]("X,a,b,m");

        // sigmal2_aaa += +0.50 <b,a||e,f>_aaaa l2_aaa(m,b,a)  // flops: o1v2L1 += o1v4L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aaa("X,e,f,m") += 0.50 * eri["aaaa_vvvv"]("b,a,e,f") * l2["aaa_Lovv"]("X,m,b,a");

        // sigmal2_aba += +0.50 <b,a||e,f>_abab l2_aab(m,b,a) @ +0.50 <a,b||e,f>_abab l2_aab(m,a,b)  // flops: o1v2L1 += o1v4L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") += eri["abab_vvvv"]("b,a,e,f") * l2["aab_Lovv"]("X,m,b,a");

        // sigmal2_abb += +0.50 <b,a||e,f>_abab l2_bab(m,b,a) @ +0.50 <a,b||e,f>_abab l2_bab(m,a,b)  // flops: o1v2L1 += o1v4L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") += eri["abab_vvvv"]("b,a,e,f") * l2["bab_Lovv"]("X,m,b,a");

        // sigmal2_bbb += +0.50 <b,a||e,f>_bbbb l2_bbb(m,b,a)  // flops: o1v2L1 += o1v4L1 | mem: o1v2L1 += o1v2L1
        sigmal2_bbb("X,e,f,m") += 0.50 * eri["bbbb_vvvv"]("b,a,e,f") * l2["bbb_Lovv"]("X,m,b,a");

        // sigmar2_aaa += +0.25 <j,i||a,b>_aaaa t2_aaaa(e,f,j,i) r2_aaa(a,b,m)  // flops: o1v2L1 += o3v2L1 o3v2L1 | mem: o1v2L1 += o3v0L1 o1v2L1
        sigmar2_aaa("X,e,f,m") += 0.25 * eri["aaaa_oovv"]("j,i,a,b") * r2["aaa_Lvvo"]("X,a,b,m") * t2["aaaa_vvoo"]("e,f,j,i");

        // sigmar2_aba += +0.25 <j,i||a,b>_abab t2_abab(e,f,j,i) r2_aba(a,b,m) @ +0.25 <j,i||b,a>_abab t2_abab(e,f,j,i) r2_aba(b,a,m) @ +0.25 <i,j||a,b>_abab t2_abab(e,f,i,j) r2_aba(a,b,m) @ +0.25 <i,j||b,a>_abab t2_abab(e,f,i,j) r2_aba(b,a,m)  // flops: o1v2L1 += o3v2L1 o3v2L1 | mem: o1v2L1 += o3v0L1 o1v2L1
        sigmar2_aba("X,e,f,m") += eri["abab_oovv"]("j,i,a,b") * r2["aba_Lvvo"]("X,a,b,m") * t2["abab_vvoo"]("e,f,j,i");

        // sigmar2_abb += +0.25 <j,i||a,b>_abab t2_abab(e,f,j,i) r2_abb(a,b,m) @ +0.25 <j,i||b,a>_abab t2_abab(e,f,j,i) r2_abb(b,a,m) @ +0.25 <i,j||a,b>_abab t2_abab(e,f,i,j) r2_abb(a,b,m) @ +0.25 <i,j||b,a>_abab t2_abab(e,f,i,j) r2_abb(b,a,m)  // flops: o1v2L1 += o3v2L1 o3v2L1 | mem: o1v2L1 += o3v0L1 o1v2L1
        sigmar2_abb("X,e,f,m") += eri["abab_oovv"]("j,i,a,b") * r2["abb_Lvvo"]("X,a,b,m") * t2["abab_vvoo"]("e,f,j,i");

        // sigmar2_bbb += +0.25 <j,i||a,b>_bbbb t2_bbbb(e,f,j,i) r2_bbb(a,b,m)  // flops: o1v2L1 += o3v2L1 o3v2L1 | mem: o1v2L1 += o3v0L1 o1v2L1
        sigmar2_bbb("X,e,f,m") += 0.25 * eri["bbbb_oovv"]("j,i,a,b") * r2["bbb_Lvvo"]("X,a,b,m") * t2["bbbb_vvoo"]("e,f,j,i");

        // sigmal2_aaa += +0.25 <j,i||e,f>_aaaa l2_aaa(m,b,a) t2_aaaa(b,a,j,i)  // flops: o1v2L1 += o3v2L1 o3v2L1 | mem: o1v2L1 += o3v0L1 o1v2L1
        sigmal2_aaa("X,e,f,m") += 0.25 * l2["aaa_Lovv"]("X,m,b,a") * t2["aaaa_vvoo"]("b,a,j,i") * eri["aaaa_oovv"]("j,i,e,f");

        // sigmal2_aba += +0.25 <j,i||e,f>_abab l2_aab(m,b,a) t2_abab(b,a,j,i) @ +0.25 <i,j||e,f>_abab l2_aab(m,b,a) t2_abab(b,a,i,j) @ +0.25 <j,i||e,f>_abab l2_aab(m,a,b) t2_abab(a,b,j,i) @ +0.25 <i,j||e,f>_abab l2_aab(m,a,b) t2_abab(a,b,i,j)  // flops: o1v2L1 += o3v2L1 o3v2L1 | mem: o1v2L1 += o3v0L1 o1v2L1
        sigmal2_aba("X,e,f,m") += l2["aab_Lovv"]("X,m,b,a") * t2["abab_vvoo"]("b,a,j,i") * eri["abab_oovv"]("j,i,e,f");

        // sigmal2_abb += +0.25 <j,i||e,f>_abab l2_bab(m,b,a) t2_abab(b,a,j,i) @ +0.25 <i,j||e,f>_abab l2_bab(m,b,a) t2_abab(b,a,i,j) @ +0.25 <j,i||e,f>_abab l2_bab(m,a,b) t2_abab(a,b,j,i) @ +0.25 <i,j||e,f>_abab l2_bab(m,a,b) t2_abab(a,b,i,j)  // flops: o1v2L1 += o3v2L1 o3v2L1 | mem: o1v2L1 += o3v0L1 o1v2L1
        sigmal2_abb("X,e,f,m") += l2["bab_Lovv"]("X,m,b,a") * t2["abab_vvoo"]("b,a,j,i") * eri["abab_oovv"]("j,i,e,f");

        // sigmal2_bbb += +0.25 <j,i||e,f>_bbbb l2_bbb(m,b,a) t2_bbbb(b,a,j,i)  // flops: o1v2L1 += o3v2L1 o3v2L1 | mem: o1v2L1 += o3v0L1 o1v2L1
        sigmal2_bbb("X,e,f,m") += 0.25 * l2["bbb_Lovv"]("X,m,b,a") * t2["bbbb_vvoo"]("b,a,j,i") * eri["bbbb_oovv"]("j,i,e,f");

        // tmps_[1_bbb_Lvvo](X,e,f,m) = 1.00 l2[aab_Lovv](X,i,a,f) * eri[abab_oovv](j,m,b,e) * t2[aaaa_vvoo](b,a,i,j) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["1_bbb_Lvvo"]("X,e,f,m")  = l2["aab_Lovv"]("X,i,a,f") * reuse_tmps_["7_abab_vvoo"]("a,e,i,m");

        // sigmal2_bbb += -1.00 P(e,f) <j,m||b,e>_abab l2_aab(i,a,f) t2_aaaa(b,a,i,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_bbb("X,e,f,m") -= tmps_["1_bbb_Lvvo"]("X,e,f,m");
        sigmal2_bbb("X,e,f,m") += tmps_["1_bbb_Lvvo"]("X,f,e,m");
        tmps_["1_bbb_Lvvo"].~TArrayD();

        // tmps_[2_bbb_Lvvo](X,f,e,m) = 1.00 l2[bbb_Lovv](X,i,a,e) * eri[bbbb_oovv](m,j,b,f) * t2[bbbb_vvoo](b,a,i,j) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["2_bbb_Lvvo"]("X,f,e,m")  = l2["bbb_Lovv"]("X,i,a,e") * reuse_tmps_["3_bbbb_vvoo"]("a,f,i,m");
        sigmal2_bbb("X,e,f,m") -= tmps_["2_bbb_Lvvo"]("X,f,e,m");

        // sigmal2_bbb += +1.00 P(e,f) <m,j||b,e>_bbbb l2_bbb(i,a,f) t2_bbbb(b,a,i,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_bbb("X,e,f,m") += tmps_["2_bbb_Lvvo"]("X,e,f,m");
        tmps_["2_bbb_Lvvo"].~TArrayD();

        // tmps_[3_bbb_Lvvo](X,e,f,m) = 1.00 eri[bbbb_vovo](a,m,f,i) * l2[bbb_Lovv](X,i,a,e) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["3_bbb_Lvvo"]("X,e,f,m")  = eri["bbbb_vovo"]("a,m,f,i") * l2["bbb_Lovv"]("X,i,a,e");

        // sigmal2_bbb += +1.00 P(e,f) <m,a||e,i>_bbbb l2_bbb(i,a,f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_bbb("X,e,f,m") -= tmps_["3_bbb_Lvvo"]("X,f,e,m");
        sigmal2_bbb("X,e,f,m") += tmps_["3_bbb_Lvvo"]("X,e,f,m");
        tmps_["3_bbb_Lvvo"].~TArrayD();

        // tmps_[4_bbb_Lvvo](X,f,e,m) = 1.00 l2[bbb_Lovv](X,i,a,e) * eri[abab_oovv](j,m,b,f) * t2[abab_vvoo](b,a,j,i) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["4_bbb_Lvvo"]("X,f,e,m")  = l2["bbb_Lovv"]("X,i,a,e") * reuse_tmps_["1_bbbb_vvoo"]("a,f,i,m");
        sigmal2_bbb("X,e,f,m") -= tmps_["4_bbb_Lvvo"]("X,f,e,m");

        // sigmal2_bbb += +1.00 P(e,f) <j,m||b,e>_abab l2_bbb(i,a,f) t2_abab(b,a,j,i)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_bbb("X,e,f,m") += tmps_["4_bbb_Lvvo"]("X,e,f,m");
        tmps_["4_bbb_Lvvo"].~TArrayD();

        // tmps_[5_bbb_Lvvo](X,e,f,m) = 1.00 l2[aab_Lovv](X,i,a,f) * eri[bbbb_oovv](m,j,b,e) * t2[abab_vvoo](a,b,i,j) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["5_bbb_Lvvo"]("X,e,f,m")  = l2["aab_Lovv"]("X,i,a,f") * reuse_tmps_["8_abab_vvoo"]("a,e,i,m");

        // sigmal2_bbb += -1.00 P(e,f) <m,j||b,e>_bbbb l2_aab(i,a,f) t2_abab(a,b,i,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_bbb("X,e,f,m") -= tmps_["5_bbb_Lvvo"]("X,e,f,m");
        sigmal2_bbb("X,e,f,m") += tmps_["5_bbb_Lvvo"]("X,f,e,m");
        tmps_["5_bbb_Lvvo"].~TArrayD();

        // tmps_[6_bbb_Lvvo](X,e,f,m) = 1.00 eri[abba_vovo](a,m,f,i) * l2[aab_Lovv](X,i,a,e) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["6_bbb_Lvvo"]("X,e,f,m")  = eri["abba_vovo"]("a,m,f,i") * l2["aab_Lovv"]("X,i,a,e");

        // sigmal2_bbb += +1.00 P(e,f) <a,m||i,e>_abab l2_aab(i,a,f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_bbb("X,e,f,m") -= tmps_["6_bbb_Lvvo"]("X,f,e,m");
        sigmal2_bbb("X,e,f,m") += tmps_["6_bbb_Lvvo"]("X,e,f,m");
        tmps_["6_bbb_Lvvo"].~TArrayD();

        // tmps_[7_aaa_Lvvo](X,f,e,m) = 1.00 eri[abab_oovv](m,j,e,b) * t2[abab_vvoo](a,b,i,j) * l2[aaa_Lovv](X,i,a,f) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["7_aaa_Lvvo"]("X,f,e,m")  = reuse_tmps_["2_aaaa_vvoo"]("a,e,i,m") * l2["aaa_Lovv"]("X,i,a,f");
        sigmal2_aaa("X,e,f,m") -= tmps_["7_aaa_Lvvo"]("X,e,f,m");

        // sigmal2_aaa += +1.00 P(e,f) <m,j||e,b>_abab l2_aaa(i,a,f) t2_abab(a,b,i,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aaa("X,e,f,m") += tmps_["7_aaa_Lvvo"]("X,f,e,m");
        tmps_["7_aaa_Lvvo"].~TArrayD();

        // tmps_[8_aaa_Lvvo](X,f,e,m) = 1.00 l2[aaa_Lovv](X,i,a,e) * eri[aaaa_oovv](m,j,b,f) * t2[aaaa_vvoo](b,a,i,j) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["8_aaa_Lvvo"]("X,f,e,m")  = l2["aaa_Lovv"]("X,i,a,e") * reuse_tmps_["4_aaaa_vvoo"]("a,f,i,m");
        sigmal2_aaa("X,e,f,m") -= tmps_["8_aaa_Lvvo"]("X,f,e,m");

        // sigmal2_aaa += +1.00 P(e,f) <m,j||b,e>_aaaa l2_aaa(i,a,f) t2_aaaa(b,a,i,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aaa("X,e,f,m") += tmps_["8_aaa_Lvvo"]("X,e,f,m");
        tmps_["8_aaa_Lvvo"].~TArrayD();

        // tmps_[9_aaa_Lvvo](X,e,f,m) = 1.00 eri[baab_vovo](a,m,f,i) * l2[bab_Lovv](X,i,e,a) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["9_aaa_Lvvo"]("X,e,f,m")  = eri["baab_vovo"]("a,m,f,i") * l2["bab_Lovv"]("X,i,e,a");
        sigmal2_aaa("X,e,f,m") -= tmps_["9_aaa_Lvvo"]("X,e,f,m");

        // sigmal2_aaa += -1.00 P(e,f) <m,a||e,i>_abab l2_bab(i,f,a)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aaa("X,e,f,m") += tmps_["9_aaa_Lvvo"]("X,f,e,m");
        tmps_["9_aaa_Lvvo"].~TArrayD();

        // tmps_[10_aaa_Lvvo](X,f,e,m) = 1.00 eri[aaaa_vovo](a,m,e,i) * l2[aaa_Lovv](X,i,a,f) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["10_aaa_Lvvo"]("X,f,e,m")  = eri["aaaa_vovo"]("a,m,e,i") * l2["aaa_Lovv"]("X,i,a,f");

        // sigmal2_aaa += +1.00 P(e,f) <m,a||e,i>_aaaa l2_aaa(i,a,f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aaa("X,e,f,m") -= tmps_["10_aaa_Lvvo"]("X,f,e,m");
        sigmal2_aaa("X,e,f,m") += tmps_["10_aaa_Lvvo"]("X,e,f,m");
        tmps_["10_aaa_Lvvo"].~TArrayD();

        // tmps_[11_aaa_Lvvo](X,e,f,m) = 1.00 l2[bab_Lovv](X,i,f,a) * eri[aaaa_oovv](m,j,b,e) * t2[abab_vvoo](b,a,j,i) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["11_aaa_Lvvo"]("X,e,f,m")  = l2["bab_Lovv"]("X,i,f,a") * reuse_tmps_["6_abab_vvoo"]("e,a,m,i");

        // sigmal2_aaa += +1.00 P(e,f) <m,j||b,e>_aaaa l2_bab(i,f,a) t2_abab(b,a,j,i)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aaa("X,e,f,m") += tmps_["11_aaa_Lvvo"]("X,e,f,m");
        sigmal2_aaa("X,e,f,m") -= tmps_["11_aaa_Lvvo"]("X,f,e,m");
        tmps_["11_aaa_Lvvo"].~TArrayD();

        // tmps_[12_aaa_Lvvo](X,e,f,m) = 1.00 l2[bab_Lovv](X,i,f,a) * eri[abab_oovv](m,j,e,b) * t2[bbbb_vvoo](b,a,i,j) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["12_aaa_Lvvo"]("X,e,f,m")  = l2["bab_Lovv"]("X,i,f,a") * reuse_tmps_["5_abab_vvoo"]("e,a,m,i");

        // sigmal2_aaa += +1.00 P(e,f) <m,j||e,b>_abab l2_bab(i,f,a) t2_bbbb(b,a,i,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aaa("X,e,f,m") += tmps_["12_aaa_Lvvo"]("X,e,f,m");
        sigmal2_aaa("X,e,f,m") -= tmps_["12_aaa_Lvvo"]("X,f,e,m");
        tmps_["12_aaa_Lvvo"].~TArrayD();

        // tmps_[13_bbb_Lvvo](X,e,f,m) = 1.00 r2[bbb_Lvvo](X,b,f,j) * eri[bbbb_oovv](j,i,a,b) * t2[bbbb_vvoo](a,e,m,i) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["13_bbb_Lvvo"]("X,e,f,m")  = r2["bbb_Lvvo"]("X,b,f,j") * reuse_tmps_["3_bbbb_vvoo"]("e,b,m,j");

        // sigmar2_bbb += +1.00 P(e,f) <j,i||a,b>_bbbb t2_bbbb(a,e,m,i) r2_bbb(b,f,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") += tmps_["13_bbb_Lvvo"]("X,e,f,m");
        sigmar2_bbb("X,e,f,m") -= tmps_["13_bbb_Lvvo"]("X,f,e,m");
        tmps_["13_bbb_Lvvo"].~TArrayD();

        // tmps_[14_bbb_Lvvo](X,f,e,m) = 1.00 r2[aba_Lvvo](X,b,e,j) * eri[aaaa_oovv](j,i,a,b) * t2[abab_vvoo](a,f,i,m) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["14_bbb_Lvvo"]("X,f,e,m")  = r2["aba_Lvvo"]("X,b,e,j") * reuse_tmps_["6_abab_vvoo"]("b,f,j,m");

        // sigmar2_bbb += -1.00 P(e,f) <j,i||a,b>_aaaa t2_abab(a,e,i,m) r2_aba(b,f,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") -= tmps_["14_bbb_Lvvo"]("X,e,f,m");
        sigmar2_bbb("X,e,f,m") += tmps_["14_bbb_Lvvo"]("X,f,e,m");
        tmps_["14_bbb_Lvvo"].~TArrayD();

        // tmps_[15_bbb_Lvvo](X,e,f,m) = 1.00 eri[bbbb_vovo](f,i,a,m) * r2[bbb_Lvvo](X,a,e,i) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["15_bbb_Lvvo"]("X,e,f,m")  = eri["bbbb_vovo"]("f,i,a,m") * r2["bbb_Lvvo"]("X,a,e,i");

        // sigmar2_bbb += +1.00 P(e,f) <i,e||a,m>_bbbb r2_bbb(a,f,i)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") -= tmps_["15_bbb_Lvvo"]("X,f,e,m");
        sigmar2_bbb("X,e,f,m") += tmps_["15_bbb_Lvvo"]("X,e,f,m");
        tmps_["15_bbb_Lvvo"].~TArrayD();

        // tmps_[16_bbb_Lvvo](X,e,f,m) = 1.00 r2[aba_Lvvo](X,b,f,j) * eri[abab_oovv](j,i,b,a) * t2[bbbb_vvoo](a,e,m,i) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["16_bbb_Lvvo"]("X,e,f,m")  = r2["aba_Lvvo"]("X,b,f,j") * reuse_tmps_["5_abab_vvoo"]("b,e,j,m");

        // sigmar2_bbb += -1.00 P(e,f) <j,i||b,a>_abab t2_bbbb(a,e,m,i) r2_aba(b,f,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") -= tmps_["16_bbb_Lvvo"]("X,e,f,m");
        sigmar2_bbb("X,e,f,m") += tmps_["16_bbb_Lvvo"]("X,f,e,m");
        tmps_["16_bbb_Lvvo"].~TArrayD();

        // tmps_[17_bbb_Lvvo](X,f,e,m) = 1.00 eri[baab_vovo](e,i,a,m) * r2[aba_Lvvo](X,a,f,i) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["17_bbb_Lvvo"]("X,f,e,m")  = eri["baab_vovo"]("e,i,a,m") * r2["aba_Lvvo"]("X,a,f,i");
        sigmar2_bbb("X,e,f,m") += tmps_["17_bbb_Lvvo"]("X,e,f,m");

        // sigmar2_bbb += +1.00 P(e,f) <i,e||a,m>_abab r2_aba(a,f,i)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") -= tmps_["17_bbb_Lvvo"]("X,f,e,m");
        tmps_["17_bbb_Lvvo"].~TArrayD();

        // tmps_[18_bbb_Lvvo](X,f,e,m) = 1.00 r2[bbb_Lvvo](X,b,e,j) * eri[abab_oovv](i,j,a,b) * t2[abab_vvoo](a,f,i,m) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["18_bbb_Lvvo"]("X,f,e,m")  = r2["bbb_Lvvo"]("X,b,e,j") * reuse_tmps_["1_bbbb_vvoo"]("f,b,m,j");
        sigmar2_bbb("X,e,f,m") -= tmps_["18_bbb_Lvvo"]("X,f,e,m");

        // sigmar2_bbb += +1.00 P(e,f) <i,j||a,b>_abab t2_abab(a,e,i,m) r2_bbb(b,f,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") += tmps_["18_bbb_Lvvo"]("X,e,f,m");
        tmps_["18_bbb_Lvvo"].~TArrayD();

        // tmps_[19_aaa_Lvvo](X,e,f,m) = 1.00 r2[aaa_Lvvo](X,b,f,j) * eri[abab_oovv](j,i,b,a) * t2[abab_vvoo](e,a,m,i) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["19_aaa_Lvvo"]("X,e,f,m")  = r2["aaa_Lvvo"]("X,b,f,j") * reuse_tmps_["2_aaaa_vvoo"]("e,b,m,j");
        sigmar2_aaa("X,e,f,m") -= tmps_["19_aaa_Lvvo"]("X,f,e,m");

        // sigmar2_aaa += +1.00 P(e,f) <j,i||b,a>_abab t2_abab(e,a,m,i) r2_aaa(b,f,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") += tmps_["19_aaa_Lvvo"]("X,e,f,m");
        tmps_["19_aaa_Lvvo"].~TArrayD();

        // tmps_[20_aaa_Lvvo](X,f,e,m) = 1.00 eri[aaaa_vovo](e,i,a,m) * r2[aaa_Lvvo](X,a,f,i) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["20_aaa_Lvvo"]("X,f,e,m")  = eri["aaaa_vovo"]("e,i,a,m") * r2["aaa_Lvvo"]("X,a,f,i");

        // sigmar2_aaa += +1.00 P(e,f) <i,e||a,m>_aaaa r2_aaa(a,f,i)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") -= tmps_["20_aaa_Lvvo"]("X,f,e,m");
        sigmar2_aaa("X,e,f,m") += tmps_["20_aaa_Lvvo"]("X,e,f,m");
        tmps_["20_aaa_Lvvo"].~TArrayD();

        // tmps_[21_aaa_Lvvo](X,e,f,m) = 1.00 r2[aaa_Lvvo](X,b,f,j) * eri[aaaa_oovv](j,i,a,b) * t2[aaaa_vvoo](a,e,m,i) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["21_aaa_Lvvo"]("X,e,f,m")  = r2["aaa_Lvvo"]("X,b,f,j") * reuse_tmps_["4_aaaa_vvoo"]("e,b,m,j");

        // sigmar2_aaa += +1.00 P(e,f) <j,i||a,b>_aaaa t2_aaaa(a,e,m,i) r2_aaa(b,f,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") += tmps_["21_aaa_Lvvo"]("X,e,f,m");
        sigmar2_aaa("X,e,f,m") -= tmps_["21_aaa_Lvvo"]("X,f,e,m");
        tmps_["21_aaa_Lvvo"].~TArrayD();

        // tmps_[22_aaa_Lvvo](X,f,e,m) = 1.00 r2[abb_Lvvo](X,e,b,j) * eri[abab_oovv](i,j,a,b) * t2[aaaa_vvoo](a,f,m,i) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["22_aaa_Lvvo"]("X,f,e,m")  = r2["abb_Lvvo"]("X,e,b,j") * reuse_tmps_["7_abab_vvoo"]("f,b,m,j");

        // sigmar2_aaa += +1.00 P(e,f) <i,j||a,b>_abab t2_aaaa(a,e,m,i) r2_abb(f,b,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") += tmps_["22_aaa_Lvvo"]("X,e,f,m");
        sigmar2_aaa("X,e,f,m") -= tmps_["22_aaa_Lvvo"]("X,f,e,m");
        tmps_["22_aaa_Lvvo"].~TArrayD();

        // tmps_[23_aaa_Lvvo](X,f,e,m) = 1.00 eri[abba_vovo](e,i,a,m) * r2[abb_Lvvo](X,f,a,i) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["23_aaa_Lvvo"]("X,f,e,m")  = eri["abba_vovo"]("e,i,a,m") * r2["abb_Lvvo"]("X,f,a,i");

        // sigmar2_aaa += -1.00 P(e,f) <e,i||m,a>_abab r2_abb(f,a,i)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") += tmps_["23_aaa_Lvvo"]("X,f,e,m");
        sigmar2_aaa("X,e,f,m") -= tmps_["23_aaa_Lvvo"]("X,e,f,m");
        tmps_["23_aaa_Lvvo"].~TArrayD();

        // tmps_[24_aaa_Lvvo](X,e,f,m) = 1.00 r2[abb_Lvvo](X,f,b,j) * eri[bbbb_oovv](j,i,a,b) * t2[abab_vvoo](e,a,m,i) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["24_aaa_Lvvo"]("X,e,f,m")  = r2["abb_Lvvo"]("X,f,b,j") * reuse_tmps_["8_abab_vvoo"]("e,b,m,j");

        // sigmar2_aaa += +1.00 P(e,f) <j,i||a,b>_bbbb t2_abab(e,a,m,i) r2_abb(f,b,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") += tmps_["24_aaa_Lvvo"]("X,e,f,m");
        sigmar2_aaa("X,e,f,m") -= tmps_["24_aaa_Lvvo"]("X,f,e,m");
        tmps_["24_aaa_Lvvo"].~TArrayD();

        // tmps_[25_bbb_Lvvo](X,f,e,m) = 1.00 l2[bbb_Lovv](X,m,a,e) * eri[bbbb_oovv](j,i,b,f) * t2[bbbb_vvoo](b,a,j,i) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["25_bbb_Lvvo"]("X,f,e,m")  = l2["bbb_Lovv"]("X,m,a,e") * reuse_tmps_["25_bb_vv"]("a,f");
        sigmal2_bbb("X,e,f,m") += tmps_["25_bbb_Lvvo"]("X,f,e,m");

        // sigmal2_bbb += -0.50 P(e,f) <j,i||b,e>_bbbb l2_bbb(m,a,f) t2_bbbb(b,a,j,i)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_bbb("X,e,f,m") -= tmps_["25_bbb_Lvvo"]("X,e,f,m");
        tmps_["25_bbb_Lvvo"].~TArrayD();

        // tmps_[26_bbb_Lvvo](X,f,e,m) = 1.00 f[bb_vv](a,e) * l2[bbb_Lovv](X,m,a,f) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["26_bbb_Lvvo"]("X,f,e,m")  = f["bb_vv"]("a,e") * l2["bbb_Lovv"]("X,m,a,f");

        // sigmal2_bbb += +1.00 P(e,f) f_bb(a,e) l2_bbb(m,a,f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_bbb("X,e,f,m") += tmps_["26_bbb_Lvvo"]("X,f,e,m");
        sigmal2_bbb("X,e,f,m") -= tmps_["26_bbb_Lvvo"]("X,e,f,m");
        tmps_["26_bbb_Lvvo"].~TArrayD();

        // tmps_[27_bbb_Lvvo](X,f,e,m) = 1.00 l2[bbb_Lovv](X,m,a,e) * eri[abab_oovv](j,i,b,f) * t2[abab_vvoo](b,a,j,i) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["27_bbb_Lvvo"]("X,f,e,m")  = l2["bbb_Lovv"]("X,m,a,e") * reuse_tmps_["22_bb_vv"]("a,f");
        sigmal2_bbb("X,e,f,m") += tmps_["27_bbb_Lvvo"]("X,f,e,m");

        // sigmal2_bbb += -0.50 P(e,f) <j,i||b,e>_abab l2_bbb(m,a,f) t2_abab(b,a,j,i) @ -0.50 P(e,f) <i,j||b,e>_abab l2_bbb(m,a,f) t2_abab(b,a,i,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_bbb("X,e,f,m") -= tmps_["27_bbb_Lvvo"]("X,e,f,m");
        tmps_["27_bbb_Lvvo"].~TArrayD();

        // tmps_[28_aaa_Lvvo](X,f,e,m) = 1.00 l2[aaa_Lovv](X,m,a,e) * eri[aaaa_oovv](j,i,b,f) * t2[aaaa_vvoo](b,a,j,i) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["28_aaa_Lvvo"]("X,f,e,m")  = l2["aaa_Lovv"]("X,m,a,e") * reuse_tmps_["26_aa_vv"]("a,f");
        sigmal2_aaa("X,e,f,m") += tmps_["28_aaa_Lvvo"]("X,f,e,m");

        // sigmal2_aaa += -0.50 P(e,f) <j,i||b,e>_aaaa l2_aaa(m,a,f) t2_aaaa(b,a,j,i)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aaa("X,e,f,m") -= tmps_["28_aaa_Lvvo"]("X,e,f,m");
        tmps_["28_aaa_Lvvo"].~TArrayD();

        // tmps_[29_aaa_Lvvo](X,f,e,m) = 1.00 eri[abab_oovv](j,i,e,b) * t2[abab_vvoo](a,b,j,i) * l2[aaa_Lovv](X,m,a,f) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["29_aaa_Lvvo"]("X,f,e,m")  = reuse_tmps_["21_aa_vv"]("a,e") * l2["aaa_Lovv"]("X,m,a,f");
        sigmal2_aaa("X,e,f,m") += tmps_["29_aaa_Lvvo"]("X,e,f,m");

        // sigmal2_aaa += -0.50 P(e,f) <j,i||e,b>_abab l2_aaa(m,a,f) t2_abab(a,b,j,i) @ -0.50 P(e,f) <i,j||e,b>_abab l2_aaa(m,a,f) t2_abab(a,b,i,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aaa("X,e,f,m") -= tmps_["29_aaa_Lvvo"]("X,f,e,m");
        tmps_["29_aaa_Lvvo"].~TArrayD();

        // tmps_[30_aaa_Lvvo](X,f,e,m) = 1.00 f[aa_vv](a,e) * l2[aaa_Lovv](X,m,a,f) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["30_aaa_Lvvo"]("X,f,e,m")  = f["aa_vv"]("a,e") * l2["aaa_Lovv"]("X,m,a,f");

        // sigmal2_aaa += +1.00 P(e,f) f_aa(a,e) l2_aaa(m,a,f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aaa("X,e,f,m") += tmps_["30_aaa_Lvvo"]("X,f,e,m");
        sigmal2_aaa("X,e,f,m") -= tmps_["30_aaa_Lvvo"]("X,e,f,m");
        tmps_["30_aaa_Lvvo"].~TArrayD();

        // tmps_[31_bbb_Lvvo](X,f,e,m) = 1.00 r2[bbb_Lvvo](X,b,e,m) * eri[abab_oovv](j,i,a,b) * t2[abab_vvoo](a,f,j,i) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["31_bbb_Lvvo"]("X,f,e,m")  = r2["bbb_Lvvo"]("X,b,e,m") * reuse_tmps_["22_bb_vv"]("f,b");
        sigmar2_bbb("X,e,f,m") += tmps_["31_bbb_Lvvo"]("X,f,e,m");

        // sigmar2_bbb += -0.50 P(e,f) <j,i||a,b>_abab t2_abab(a,e,j,i) r2_bbb(b,f,m) @ -0.50 P(e,f) <i,j||a,b>_abab t2_abab(a,e,i,j) r2_bbb(b,f,m)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") -= tmps_["31_bbb_Lvvo"]("X,e,f,m");
        tmps_["31_bbb_Lvvo"].~TArrayD();

        // tmps_[32_bbb_Lvvo](X,f,e,m) = 1.00 r1[b_Lv](X,b) * eri[baab_vovv](f,i,a,b) * t2[abab_vvoo](a,e,i,m) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["32_bbb_Lvvo"]("X,f,e,m")  = r1["b_Lv"]("X,b") * reuse_tmps_["12_bbbb_vvvo"]("e,b,f,m");
        sigmar2_bbb("X,e,f,m") += tmps_["32_bbb_Lvvo"]("X,f,e,m");

        // sigmar2_bbb += +1.00 P(e,f) <i,e||a,b>_abab t2_abab(a,f,i,m) r1_b(b)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") -= tmps_["32_bbb_Lvvo"]("X,e,f,m");
        tmps_["32_bbb_Lvvo"].~TArrayD();

        // tmps_[33_bbb_Lvvo](X,e,f,m) = 1.00 r2[bbb_Lvvo](X,b,f,m) * eri[bbbb_oovv](j,i,a,b) * t2[bbbb_vvoo](a,e,j,i) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["33_bbb_Lvvo"]("X,e,f,m")  = r2["bbb_Lvvo"]("X,b,f,m") * reuse_tmps_["25_bb_vv"]("e,b");

        // sigmar2_bbb += -0.50 P(e,f) <j,i||a,b>_bbbb t2_bbbb(a,e,j,i) r2_bbb(b,f,m)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") -= tmps_["33_bbb_Lvvo"]("X,e,f,m");
        sigmar2_bbb("X,e,f,m") += tmps_["33_bbb_Lvvo"]("X,f,e,m");
        tmps_["33_bbb_Lvvo"].~TArrayD();

        // tmps_[34_bbb_Lvvo](X,e,f,m) = 1.00 r1[b_Lv](X,b) * eri[bbbb_vovv](e,i,a,b) * t2[bbbb_vvoo](a,f,m,i) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["34_bbb_Lvvo"]("X,e,f,m")  = r1["b_Lv"]("X,b") * reuse_tmps_["13_bbbb_vvvo"]("f,b,e,m");

        // sigmar2_bbb += -1.00 P(e,f) <i,e||a,b>_bbbb t2_bbbb(a,f,m,i) r1_b(b)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") += tmps_["34_bbb_Lvvo"]("X,e,f,m");
        sigmar2_bbb("X,e,f,m") -= tmps_["34_bbb_Lvvo"]("X,f,e,m");
        tmps_["34_bbb_Lvvo"].~TArrayD();

        // tmps_[35_bbb_Lvvo](X,e,f,m) = 1.00 f[bb_vv](f,a) * r2[bbb_Lvvo](X,a,e,m) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["35_bbb_Lvvo"]("X,e,f,m")  = f["bb_vv"]("f,a") * r2["bbb_Lvvo"]("X,a,e,m");
        sigmar2_bbb("X,e,f,m") -= tmps_["35_bbb_Lvvo"]("X,e,f,m");

        // sigmar2_bbb += +1.00 P(e,f) f_bb(e,a) r2_bbb(a,f,m)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") += tmps_["35_bbb_Lvvo"]("X,f,e,m");
        tmps_["35_bbb_Lvvo"].~TArrayD();

        // tmps_[36_aaa_Lvvo](X,e,f,m) = 1.00 r2[aaa_Lvvo](X,b,f,m) * eri[abab_oovv](j,i,b,a) * t2[abab_vvoo](e,a,j,i) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["36_aaa_Lvvo"]("X,e,f,m")  = r2["aaa_Lvvo"]("X,b,f,m") * reuse_tmps_["21_aa_vv"]("e,b");
        sigmar2_aaa("X,e,f,m") += tmps_["36_aaa_Lvvo"]("X,f,e,m");

        // sigmar2_aaa += -0.50 P(e,f) <j,i||b,a>_abab t2_abab(e,a,j,i) r2_aaa(b,f,m) @ -0.50 P(e,f) <i,j||b,a>_abab t2_abab(e,a,i,j) r2_aaa(b,f,m)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") -= tmps_["36_aaa_Lvvo"]("X,e,f,m");
        tmps_["36_aaa_Lvvo"].~TArrayD();

        // tmps_[37_aaa_Lvvo](X,e,f,m) = 1.00 r2[aaa_Lvvo](X,b,f,m) * eri[aaaa_oovv](j,i,a,b) * t2[aaaa_vvoo](a,e,j,i) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["37_aaa_Lvvo"]("X,e,f,m")  = r2["aaa_Lvvo"]("X,b,f,m") * reuse_tmps_["26_aa_vv"]("e,b");

        // sigmar2_aaa += -0.50 P(e,f) <j,i||a,b>_aaaa t2_aaaa(a,e,j,i) r2_aaa(b,f,m)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") -= tmps_["37_aaa_Lvvo"]("X,e,f,m");
        sigmar2_aaa("X,e,f,m") += tmps_["37_aaa_Lvvo"]("X,f,e,m");
        tmps_["37_aaa_Lvvo"].~TArrayD();

        // tmps_[38_aaa_Lvvo](X,e,f,m) = 1.00 f[aa_vv](f,a) * r2[aaa_Lvvo](X,a,e,m) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["38_aaa_Lvvo"]("X,e,f,m")  = f["aa_vv"]("f,a") * r2["aaa_Lvvo"]("X,a,e,m");

        // sigmar2_aaa += +1.00 P(e,f) f_aa(e,a) r2_aaa(a,f,m)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") += tmps_["38_aaa_Lvvo"]("X,f,e,m");
        sigmar2_aaa("X,e,f,m") -= tmps_["38_aaa_Lvvo"]("X,e,f,m");
        tmps_["38_aaa_Lvvo"].~TArrayD();

        // tmps_[39_aaa_Lvvo](X,e,f,m) = 1.00 r1[a_Lv](X,b) * eri[abab_vovv](e,i,b,a) * t2[abab_vvoo](f,a,m,i) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["39_aaa_Lvvo"]("X,e,f,m")  = r1["a_Lv"]("X,b") * reuse_tmps_["15_aaaa_vvvo"]("f,b,e,m");

        // sigmar2_aaa += +1.00 P(e,f) <e,i||b,a>_abab t2_abab(f,a,m,i) r1_a(b)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") += tmps_["39_aaa_Lvvo"]("X,e,f,m");
        sigmar2_aaa("X,e,f,m") -= tmps_["39_aaa_Lvvo"]("X,f,e,m");
        tmps_["39_aaa_Lvvo"].~TArrayD();

        // tmps_[40_aaa_Lvvo](X,f,e,m) = 1.00 r1[a_Lv](X,b) * eri[aaaa_vovv](f,i,a,b) * t2[aaaa_vvoo](a,e,m,i) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
        tmps_["40_aaa_Lvvo"]("X,f,e,m")  = r1["a_Lv"]("X,b") * reuse_tmps_["14_aaaa_vvvo"]("e,b,f,m");
        sigmar2_aaa("X,e,f,m") -= tmps_["40_aaa_Lvvo"]("X,f,e,m");

        // sigmar2_aaa += -1.00 P(e,f) <i,e||a,b>_aaaa t2_aaaa(a,f,m,i) r1_a(b)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") += tmps_["40_aaa_Lvvo"]("X,e,f,m");
        tmps_["40_aaa_Lvvo"].~TArrayD();

        // tmps_[41_b_Lo](X,j) = 1.00 l2[aab_Lovv](X,i,b,a) * t2[abab_vvoo](b,a,i,j) // flops: o1v0L1 = o2v2L1 | mem: o1v0L1 = o1v0L1
        tmps_["41_b_Lo"]("X,j")  = l2["aab_Lovv"]("X,i,b,a") * t2["abab_vvoo"]("b,a,i,j");

        // sigmal1_b += +0.50 f_bb(j,e) l2_aab(i,b,a) t2_abab(b,a,i,j) @ +0.50 f_bb(j,e) l2_aab(i,a,b) t2_abab(a,b,i,j)  // flops: o0v1L1 += o1v1L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") += tmps_["41_b_Lo"]("X,j") * f["bb_ov"]("j,e");

        // sigmal2_aba += -0.50 <m,j||e,f>_abab l2_aab(i,b,a) t2_abab(b,a,i,j) @ -0.50 <m,j||e,f>_abab l2_aab(i,a,b) t2_abab(a,b,i,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") -= tmps_["41_b_Lo"]("X,j") * eri["abab_oovv"]("m,j,e,f");

        // sigmal2_bbb += -0.50 <m,j||e,f>_bbbb l2_aab(i,b,a) t2_abab(b,a,i,j) @ -0.50 <m,j||e,f>_bbbb l2_aab(i,a,b) t2_abab(a,b,i,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_bbb("X,e,f,m") -= tmps_["41_b_Lo"]("X,j") * eri["bbbb_oovv"]("m,j,e,f");
        tmps_["41_b_Lo"].~TArrayD();

        // tmps_[42_a_Lo](X,j) = 0.50 l2[aaa_Lovv](X,i,b,a) * t2[aaaa_vvoo](b,a,i,j) // flops: o1v0L1 = o2v2L1 | mem: o1v0L1 = o1v0L1
        tmps_["42_a_Lo"]("X,j")  = 0.50 * l2["aaa_Lovv"]("X,i,b,a") * t2["aaaa_vvoo"]("b,a,i,j");

        // sigmal1_a += +0.50 f_aa(j,e) l2_aaa(i,b,a) t2_aaaa(b,a,i,j)  // flops: o0v1L1 += o1v1L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") += tmps_["42_a_Lo"]("X,j") * f["aa_ov"]("j,e");

        // sigmal2_aaa += -0.50 <m,j||e,f>_aaaa l2_aaa(i,b,a) t2_aaaa(b,a,i,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aaa("X,e,f,m") -= tmps_["42_a_Lo"]("X,j") * eri["aaaa_oovv"]("m,j,e,f");

        // sigmal2_abb += +0.50 <j,m||e,f>_abab l2_aaa(i,b,a) t2_aaaa(b,a,i,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") += tmps_["42_a_Lo"]("X,j") * eri["abab_oovv"]("j,m,e,f");
        tmps_["42_a_Lo"].~TArrayD();

        // tmps_[43_a_Lo](X,j) = 1.00 l2[bab_Lovv](X,i,b,a) * t2[abab_vvoo](b,a,j,i) // flops: o1v0L1 = o2v2L1 | mem: o1v0L1 = o1v0L1
        tmps_["43_a_Lo"]("X,j")  = l2["bab_Lovv"]("X,i,b,a") * t2["abab_vvoo"]("b,a,j,i");

        // sigmal1_a += -0.50 f_aa(j,e) l2_bab(i,b,a) t2_abab(b,a,j,i) @ -0.50 f_aa(j,e) l2_bab(i,a,b) t2_abab(a,b,j,i)  // flops: o0v1L1 += o1v1L1 | mem: o0v1L1 += o0v1L1
        sigmal1_a("X,e") -= tmps_["43_a_Lo"]("X,j") * f["aa_ov"]("j,e");

        // sigmal2_aaa += +0.50 <m,j||e,f>_aaaa l2_bab(i,b,a) t2_abab(b,a,j,i) @ +0.50 <m,j||e,f>_aaaa l2_bab(i,a,b) t2_abab(a,b,j,i)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aaa("X,e,f,m") += tmps_["43_a_Lo"]("X,j") * eri["aaaa_oovv"]("m,j,e,f");

        // sigmal2_abb += -0.50 <j,m||e,f>_abab l2_bab(i,b,a) t2_abab(b,a,j,i) @ -0.50 <j,m||e,f>_abab l2_bab(i,a,b) t2_abab(a,b,j,i)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_abb("X,e,f,m") -= tmps_["43_a_Lo"]("X,j") * eri["abab_oovv"]("j,m,e,f");
        tmps_["43_a_Lo"].~TArrayD();

        // tmps_[44_b_Lo](X,j) = 0.50 l2[bbb_Lovv](X,i,b,a) * t2[bbbb_vvoo](b,a,i,j) // flops: o1v0L1 = o2v2L1 | mem: o1v0L1 = o1v0L1
        tmps_["44_b_Lo"]("X,j")  = 0.50 * l2["bbb_Lovv"]("X,i,b,a") * t2["bbbb_vvoo"]("b,a,i,j");

        // sigmal1_b += +0.50 f_bb(j,e) l2_bbb(i,b,a) t2_bbbb(b,a,i,j)  // flops: o0v1L1 += o1v1L1 | mem: o0v1L1 += o0v1L1
        sigmal1_b("X,e") += tmps_["44_b_Lo"]("X,j") * f["bb_ov"]("j,e");

        // sigmal2_aba += -0.50 <m,j||e,f>_abab l2_bbb(i,b,a) t2_bbbb(b,a,i,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_aba("X,e,f,m") -= tmps_["44_b_Lo"]("X,j") * eri["abab_oovv"]("m,j,e,f");

        // sigmal2_bbb += -0.50 <m,j||e,f>_bbbb l2_bbb(i,b,a) t2_bbbb(b,a,i,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmal2_bbb("X,e,f,m") -= tmps_["44_b_Lo"]("X,j") * eri["bbbb_oovv"]("m,j,e,f");
        tmps_["44_b_Lo"].~TArrayD();

        // tmps_[45_b_Lo](X,i) = 0.50 eri[bbbb_oovv](j,i,a,b) * r2[bbb_Lvvo](X,a,b,j) // flops: o1v0L1 = o2v2L1 | mem: o1v0L1 = o1v0L1
        tmps_["45_b_Lo"]("X,i")  = 0.50 * eri["bbbb_oovv"]("j,i,a,b") * r2["bbb_Lvvo"]("X,a,b,j");

        // sigmar2_aba += -0.50 <j,i||a,b>_bbbb t2_abab(e,f,m,i) r2_bbb(a,b,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= tmps_["45_b_Lo"]("X,i") * t2["abab_vvoo"]("e,f,m,i");

        // sigmar2_bbb += -0.50 <j,i||a,b>_bbbb t2_bbbb(e,f,m,i) r2_bbb(a,b,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") -= tmps_["45_b_Lo"]("X,i") * t2["bbbb_vvoo"]("e,f,m,i");
        tmps_["45_b_Lo"].~TArrayD();

        // tmps_[46_b_Lo](X,i) = 1.00 eri[abab_oovv](j,i,a,b) * r2[aba_Lvvo](X,a,b,j) // flops: o1v0L1 = o2v2L1 | mem: o1v0L1 = o1v0L1
        tmps_["46_b_Lo"]("X,i")  = eri["abab_oovv"]("j,i,a,b") * r2["aba_Lvvo"]("X,a,b,j");

        // sigmar2_aba += -0.50 <j,i||a,b>_abab t2_abab(e,f,m,i) r2_aba(a,b,j) @ -0.50 <j,i||b,a>_abab t2_abab(e,f,m,i) r2_aba(b,a,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") -= tmps_["46_b_Lo"]("X,i") * t2["abab_vvoo"]("e,f,m,i");

        // sigmar2_bbb += -0.50 <j,i||a,b>_abab t2_bbbb(e,f,m,i) r2_aba(a,b,j) @ -0.50 <j,i||b,a>_abab t2_bbbb(e,f,m,i) r2_aba(b,a,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") -= tmps_["46_b_Lo"]("X,i") * t2["bbbb_vvoo"]("e,f,m,i");
        tmps_["46_b_Lo"].~TArrayD();

        // tmps_[47_a_Lo](X,i) = 1.00 eri[abab_oovv](i,j,a,b) * r2[abb_Lvvo](X,a,b,j) // flops: o1v0L1 = o2v2L1 | mem: o1v0L1 = o1v0L1
        tmps_["47_a_Lo"]("X,i")  = eri["abab_oovv"]("i,j,a,b") * r2["abb_Lvvo"]("X,a,b,j");

        // sigmar2_aaa += +0.50 <i,j||a,b>_abab t2_aaaa(e,f,m,i) r2_abb(a,b,j) @ +0.50 <i,j||b,a>_abab t2_aaaa(e,f,m,i) r2_abb(b,a,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") += tmps_["47_a_Lo"]("X,i") * t2["aaaa_vvoo"]("e,f,m,i");

        // sigmar2_abb += -0.50 <i,j||a,b>_abab t2_abab(e,f,i,m) r2_abb(a,b,j) @ -0.50 <i,j||b,a>_abab t2_abab(e,f,i,m) r2_abb(b,a,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") -= tmps_["47_a_Lo"]("X,i") * t2["abab_vvoo"]("e,f,i,m");
        tmps_["47_a_Lo"].~TArrayD();

        // tmps_[48_a_Lo](X,i) = 0.50 eri[aaaa_oovv](j,i,a,b) * r2[aaa_Lvvo](X,a,b,j) // flops: o1v0L1 = o2v2L1 | mem: o1v0L1 = o1v0L1
        tmps_["48_a_Lo"]("X,i")  = 0.50 * eri["aaaa_oovv"]("j,i,a,b") * r2["aaa_Lvvo"]("X,a,b,j");

        // sigmar2_aaa += -0.50 <j,i||a,b>_aaaa t2_aaaa(e,f,m,i) r2_aaa(a,b,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") -= tmps_["48_a_Lo"]("X,i") * t2["aaaa_vvoo"]("e,f,m,i");

        // sigmar2_abb += +0.50 <j,i||a,b>_aaaa t2_abab(e,f,i,m) r2_aaa(a,b,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") += tmps_["48_a_Lo"]("X,i") * t2["abab_vvoo"]("e,f,i,m");
        tmps_["48_a_Lo"].~TArrayD();

        // tmps_[49_b_Lo](X,i) = 1.00 f[bb_ov](i,a) * r1[b_Lv](X,a) // flops: o1v0L1 = o1v1L1 | mem: o1v0L1 = o1v0L1
        tmps_["49_b_Lo"]("X,i")  = f["bb_ov"]("i,a") * r1["b_Lv"]("X,a");

        // sigmar2_aba += +1.00 f_bb(i,a) t2_abab(e,f,m,i) r1_b(a)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aba("X,e,f,m") += tmps_["49_b_Lo"]("X,i") * t2["abab_vvoo"]("e,f,m,i");

        // sigmar2_bbb += +1.00 f_bb(i,a) t2_bbbb(e,f,m,i) r1_b(a)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_bbb("X,e,f,m") += tmps_["49_b_Lo"]("X,i") * t2["bbbb_vvoo"]("e,f,m,i");
        tmps_["49_b_Lo"].~TArrayD();

        // tmps_[50_a_Lo](X,i) = 1.00 f[aa_ov](i,a) * r1[a_Lv](X,a) // flops: o1v0L1 = o1v1L1 | mem: o1v0L1 = o1v0L1
        tmps_["50_a_Lo"]("X,i")  = f["aa_ov"]("i,a") * r1["a_Lv"]("X,a");

        // sigmar2_aaa += +1.00 f_aa(i,a) t2_aaaa(e,f,m,i) r1_a(a)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_aaa("X,e,f,m") += tmps_["50_a_Lo"]("X,i") * t2["aaaa_vvoo"]("e,f,m,i");

        // sigmar2_abb += -1.00 f_aa(i,a) t2_abab(e,f,i,m) r1_a(a)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
        sigmar2_abb("X,e,f,m") -= tmps_["50_a_Lo"]("X,i") * t2["abab_vvoo"]("e,f,i,m");
        tmps_["50_a_Lo"].~TArrayD();

    }

    double cc_energy = cc_wfn_->cc_energy_;

    // add ground state energy
    sigmar1_a("I,e") += r1["a_Lv"]("I,e") * cc_energy;
    sigmar1_b("I,e") += r1["b_Lv"]("I,e") * cc_energy;

    sigmal1_a("I,e") += l1["a_Lv"]("I,e") * cc_energy;
    sigmal1_b("I,e") += l1["b_Lv"]("I,e") * cc_energy;

    sigmar2_aaa("I,e,f,n") += r2["aaa_Lvvo"]("I,e,f,n") * cc_energy;
    sigmar2_abb("I,e,f,n") += r2["abb_Lvvo"]("I,e,f,n") * cc_energy;
    sigmar2_aba("I,e,f,n") += r2["aba_Lvvo"]("I,e,f,n") * cc_energy;
    sigmar2_bbb("I,e,f,n") += r2["bbb_Lvvo"]("I,e,f,n") * cc_energy;

    sigmal2_aaa("I,e,f,n") += l2["aaa_Lovv"]("I,n,e,f") * cc_energy;
    sigmal2_abb("I,e,f,n") += l2["bab_Lovv"]("I,n,e,f") * cc_energy;
    sigmal2_aba("I,e,f,n") += l2["aab_Lovv"]("I,n,e,f") * cc_energy;
    sigmal2_bbb("I,e,f,n") += l2["bbb_Lovv"]("I,n,e,f") * cc_energy;

    world_.gop.fence();
}