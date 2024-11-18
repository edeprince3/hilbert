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

#include "cc_cavity/include/ccsd/eom_ea_ccsd.h"

double* hilbert::EOM_EA_CCSD::build_ss_diagonal() {
    throw PsiException("EOM_EE_CCSD::build_ss_diagonal() is not implemented", __FILE__, __LINE__);
    return nullptr;
}

void hilbert::EOM_EA_CCSD::build_hamiltonian() {
    // ensure that the integrals are transformed with t1 amplitudes
    if (!cc_wfn_->has_t1_integrals_) cc_wfn_->transform_integrals(true);

    // get reference to electronic integrals
    TArrayMap &V_blks_ = cc_wfn_->V_blks_;
    TArrayMap &F_blks_ = cc_wfn_->F_blks_;
    TArrayMap &Id_blks_ = cc_wfn_->Id_blks_;
    TArrayMap &eri = cc_wfn_->V_blks_;
    TArrayMap &F = cc_wfn_->F_blks_;
    TArrayMap &Id = cc_wfn_->Id_blks_;

    TArrayMap t1, t2;

    t1["aa_vo"] =     cc_wfn_->amplitudes_["t1_aa"];
    t1["bb_vo"] =     cc_wfn_->amplitudes_["t1_bb"];
    t2["aaaa_vvoo"] = cc_wfn_->amplitudes_["t2_aaaa"];
    t2["abab_vvoo"] = cc_wfn_->amplitudes_["t2_abab"];
    t2["bbbb_vvoo"] = cc_wfn_->amplitudes_["t2_bbbb"];


    TArrayMap tmps_, perm_tmps;
    std::map<std::string, double> scalars_;

    TArrayD H11_a_a;
    TArrayD H11_b_b;
    TArrayD H12_a_aaa;
    TArrayD H12_a_abb;
    TArrayD H12_b_aba;
    TArrayD H12_b_bbb;
    TArrayD H21_aaa_a;
    TArrayD H21_aba_b;
    TArrayD H21_abb_a;
    TArrayD H21_bbb_b;
    TArrayD H22_aaa_aaa;
    TArrayD H22_aaa_abb;
    TArrayD H22_aba_aba;
    TArrayD H22_aba_bbb;
    TArrayD H22_abb_aaa;
    TArrayD H22_abb_abb;
    TArrayD H22_bbb_aba;
    TArrayD H22_bbb_bbb;

    {
        // H11_a_a = +1.00 f_aa(a,e)  // flops: o0v2 = o0v2 | mem: o0v2 = o0v2
        H11_a_a("a,e")  = F["aa_vv"]("a,e");

        // H11_b_b = +1.00 f_bb(a,e)  // flops: o0v2 = o0v2 | mem: o0v2 = o0v2
        H11_b_b("a,e")  = F["bb_vv"]("a,e");

        // H12_a_aaa = -1.00 <m,a||e,f>_aaaa  // flops: o1v3 = o1v3 | mem: o1v3 = o1v3
        H12_a_aaa("a,e,f,m")  = eri["aaaa_vovv"]("a,m,e,f");

        // H12_a_abb = +1.00 <a,m||e,f>_abab  // flops: o1v3 = o1v3 | mem: o1v3 = o1v3
        H12_a_abb("a,e,f,m")  = eri["abab_vovv"]("a,m,e,f");

        // H12_b_aba = -1.00 <m,a||e,f>_abab  // flops: o1v3 = o1v3 | mem: o1v3 = o1v3
        H12_b_aba("a,e,f,m")  = eri["baab_vovv"]("a,m,e,f");

        // H12_b_bbb = -1.00 <m,a||e,f>_bbbb  // flops: o1v3 = o1v3 | mem: o1v3 = o1v3
        H12_b_bbb("a,e,f,m")  = eri["bbbb_vovv"]("a,m,e,f");

        // H21_aaa_a = +1.00 <a,b||e,i>_aaaa  // flops: o1v3 = o1v3 | mem: o1v3 = o1v3
        H21_aaa_a("a,b,i,e")  = eri["aaaa_vvvo"]("a,b,e,i");

        // H21_aba_b = -1.00 <a,b||i,e>_abab  // flops: o1v3 = o1v3 | mem: o1v3 = o1v3
        H21_aba_b("a,b,i,e")  = eri["abba_vvvo"]("a,b,e,i");

        // H21_abb_a = +1.00 <a,b||e,i>_abab  // flops: o1v3 = o1v3 | mem: o1v3 = o1v3
        H21_abb_a("a,b,i,e")  = eri["abab_vvvo"]("a,b,e,i");

        // H21_bbb_b = +1.00 <a,b||e,i>_bbbb  // flops: o1v3 = o1v3 | mem: o1v3 = o1v3
        H21_bbb_b("a,b,i,e")  = eri["bbbb_vvvo"]("a,b,e,i");

        // H12_a_aaa += +1.00 d_aa(a,e) f_aa(m,f)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H12_a_aaa("a,e,f,m") += Id["aa_vv"]("a,e") * F["aa_ov"]("m,f");

        // H12_a_aaa += -1.00 d_aa(a,f) f_aa(m,e)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H12_a_aaa("a,e,f,m") -= Id["aa_vv"]("a,f") * F["aa_ov"]("m,e");

        // H12_a_abb += +1.00 d_aa(a,e) f_bb(m,f)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H12_a_abb("a,e,f,m") += Id["aa_vv"]("a,e") * F["bb_ov"]("m,f");

        // H12_b_aba += -1.00 d_bb(a,f) f_aa(m,e)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H12_b_aba("a,e,f,m") -= Id["bb_vv"]("a,f") * F["aa_ov"]("m,e");

        // H12_b_bbb += +1.00 d_bb(a,e) f_bb(m,f)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H12_b_bbb("a,e,f,m") += Id["bb_vv"]("a,e") * F["bb_ov"]("m,f");

        // H12_b_bbb += -1.00 d_bb(a,f) f_bb(m,e)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H12_b_bbb("a,e,f,m") -= Id["bb_vv"]("a,f") * F["bb_ov"]("m,e");

        // H21_aaa_a += +1.00 d_aa(a,e) f_aa(b,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aaa_a("a,b,i,e") += Id["aa_vv"]("a,e") * F["aa_vo"]("b,i");

        // H21_aaa_a += -1.00 d_aa(b,e) f_aa(a,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aaa_a("a,b,i,e") -= Id["aa_vv"]("b,e") * F["aa_vo"]("a,i");

        // H21_aba_b += -1.00 d_bb(b,e) f_aa(a,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aba_b("a,b,i,e") -= Id["bb_vv"]("b,e") * F["aa_vo"]("a,i");

        // H21_abb_a += +1.00 d_aa(a,e) f_bb(b,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_abb_a("a,b,i,e") += Id["aa_vv"]("a,e") * F["bb_vo"]("b,i");

        // H21_bbb_b += -1.00 d_bb(b,e) f_bb(a,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_bbb_b("a,b,i,e") -= Id["bb_vv"]("b,e") * F["bb_vo"]("a,i");

        // H21_bbb_b += +1.00 d_bb(a,e) f_bb(b,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_bbb_b("a,b,i,e") += Id["bb_vv"]("a,e") * F["bb_vo"]("b,i");

        // H21_aaa_a += +1.00 f_aa(j,e) t2_aaaa(a,b,i,j)  // flops: o1v3 += o2v3 | mem: o1v3 += o1v3
        H21_aaa_a("a,b,i,e") += F["aa_ov"]("j,e") * t2["aaaa_vvoo"]("a,b,i,j");

        // H21_aba_b += +1.00 f_bb(j,e) t2_abab(a,b,i,j)  // flops: o1v3 += o2v3 | mem: o1v3 += o1v3
        H21_aba_b("a,b,i,e") += F["bb_ov"]("j,e") * t2["abab_vvoo"]("a,b,i,j");

        // H21_abb_a += -1.00 f_aa(j,e) t2_abab(a,b,j,i)  // flops: o1v3 += o2v3 | mem: o1v3 += o1v3
        H21_abb_a("a,b,i,e") -= F["aa_ov"]("j,e") * t2["abab_vvoo"]("a,b,j,i");

        // H21_bbb_b += +1.00 f_bb(j,e) t2_bbbb(a,b,i,j)  // flops: o1v3 += o2v3 | mem: o1v3 += o1v3
        H21_bbb_b("a,b,i,e") += F["bb_ov"]("j,e") * t2["bbbb_vvoo"]("a,b,i,j");

        // H21_aaa_a += +0.50 <k,j||e,i>_aaaa t2_aaaa(a,b,k,j)  // flops: o1v3 += o3v3 | mem: o1v3 += o1v3
        H21_aaa_a("a,b,i,e") += 0.50 * eri["aaaa_oovo"]("k,j,e,i") * t2["aaaa_vvoo"]("a,b,k,j");

        // H21_aba_b += -0.50 <k,j||i,e>_abab t2_abab(a,b,k,j) @ -0.50 <j,k||i,e>_abab t2_abab(a,b,j,k)  // flops: o1v3 += o3v3 | mem: o1v3 += o1v3
        H21_aba_b("a,b,i,e") += eri["abba_oovo"]("k,j,e,i") * t2["abab_vvoo"]("a,b,k,j");

        // H21_abb_a += +0.50 <k,j||e,i>_abab t2_abab(a,b,k,j) @ +0.50 <j,k||e,i>_abab t2_abab(a,b,j,k)  // flops: o1v3 += o3v3 | mem: o1v3 += o1v3
        H21_abb_a("a,b,i,e") += eri["abab_oovo"]("k,j,e,i") * t2["abab_vvoo"]("a,b,k,j");

        // H21_bbb_b += +0.50 <k,j||e,i>_bbbb t2_bbbb(a,b,k,j)  // flops: o1v3 += o3v3 | mem: o1v3 += o1v3
        H21_bbb_b("a,b,i,e") += 0.50 * eri["bbbb_oovo"]("k,j,e,i") * t2["bbbb_vvoo"]("a,b,k,j");

        // H21_aba_b += -1.00 <j,b||c,e>_bbbb t2_abab(a,c,i,j)  // flops: o1v3 += o2v4 | mem: o1v3 += o1v3
        H21_aba_b("a,b,i,e") += eri["bbbb_vovv"]("b,j,c,e") * t2["abab_vvoo"]("a,c,i,j");

        // H21_aba_b += +1.00 <a,j||c,e>_abab t2_abab(c,b,i,j)  // flops: o1v3 += o2v4 | mem: o1v3 += o1v3
        H21_aba_b("a,b,i,e") += eri["abab_vovv"]("a,j,c,e") * t2["abab_vvoo"]("c,b,i,j");

        // H21_aba_b += +1.00 <j,b||c,e>_abab t2_aaaa(c,a,i,j)  // flops: o1v3 += o2v4 | mem: o1v3 += o1v3
        H21_aba_b("a,b,i,e") -= eri["baab_vovv"]("b,j,c,e") * t2["aaaa_vvoo"]("c,a,i,j");

        // H21_abb_a += -1.00 <j,b||e,c>_abab t2_abab(a,c,j,i)  // flops: o1v3 += o2v4 | mem: o1v3 += o1v3
        H21_abb_a("a,b,i,e") += eri["baab_vovv"]("b,j,e,c") * t2["abab_vvoo"]("a,c,j,i");

        // H21_abb_a += -1.00 <a,j||e,c>_abab t2_bbbb(c,b,i,j)  // flops: o1v3 += o2v4 | mem: o1v3 += o1v3
        H21_abb_a("a,b,i,e") -= eri["abab_vovv"]("a,j,e,c") * t2["bbbb_vvoo"]("c,b,i,j");

        // H21_abb_a += +1.00 <j,a||c,e>_aaaa t2_abab(c,b,j,i)  // flops: o1v3 += o2v4 | mem: o1v3 += o1v3
        H21_abb_a("a,b,i,e") -= eri["aaaa_vovv"]("a,j,c,e") * t2["abab_vvoo"]("c,b,j,i");

        // H22_aaa_aaa = +1.00 d_aa(a,e) <m,b||f,i>_aaaa  // flops: o2v4 = o2v4 | mem: o2v4 = o2v4
        H22_aaa_aaa("a,b,i,e,f,m")  = -1.00 * Id["aa_vv"]("a,e") * eri["aaaa_vovo"]("b,m,f,i");

        // H22_aaa_aaa += -1.00 d_aa(b,e) <m,a||f,i>_aaaa  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") += Id["aa_vv"]("b,e") * eri["aaaa_vovo"]("a,m,f,i");

        // H22_aaa_aaa += -1.00 d_aa(a,f) <m,b||e,i>_aaaa  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") += Id["aa_vv"]("a,f") * eri["aaaa_vovo"]("b,m,e,i");

        // H22_aaa_aaa += +1.00 d_aa(i,m) <a,b||e,f>_aaaa  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") += Id["aa_oo"]("i,m") * eri["aaaa_vvvv"]("a,b,e,f");

        // H22_aaa_aaa += +1.00 d_aa(b,f) <m,a||e,i>_aaaa  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") -= Id["aa_vv"]("b,f") * eri["aaaa_vovo"]("a,m,e,i");

        // H22_aaa_abb = +1.00 d_aa(a,e) <b,m||i,f>_abab  // flops: o2v4 = o2v4 | mem: o2v4 = o2v4
        H22_aaa_abb("a,b,i,e,f,m")  = -1.00 * Id["aa_vv"]("a,e") * eri["abba_vovo"]("b,m,f,i");

        // H22_aaa_abb += -1.00 d_aa(b,e) <a,m||i,f>_abab  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_abb("a,b,i,e,f,m") += Id["aa_vv"]("b,e") * eri["abba_vovo"]("a,m,f,i");

        // H22_aba_aba = +1.00 d_aa(i,m) <a,b||e,f>_abab  // flops: o2v4 = o2v4 | mem: o2v4 = o2v4
        H22_aba_aba("a,b,i,e,f,m")  = Id["aa_oo"]("i,m") * eri["abab_vvvv"]("a,b,e,f");

        // H22_aba_aba += -1.00 d_aa(a,e) <m,b||i,f>_abab  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_aba("a,b,i,e,f,m") -= Id["aa_vv"]("a,e") * eri["baba_vovo"]("b,m,f,i");

        // H22_aba_aba += +1.00 d_bb(b,f) <m,a||e,i>_aaaa  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_aba("a,b,i,e,f,m") -= Id["bb_vv"]("b,f") * eri["aaaa_vovo"]("a,m,e,i");

        // H22_aba_bbb = -1.00 d_bb(b,e) <a,m||i,f>_abab  // flops: o2v4 = o2v4 | mem: o2v4 = o2v4
        H22_aba_bbb("a,b,i,e,f,m")  = Id["bb_vv"]("b,e") * eri["abba_vovo"]("a,m,f,i");

        // H22_aba_bbb += +1.00 d_bb(b,f) <a,m||i,e>_abab  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_bbb("a,b,i,e,f,m") -= Id["bb_vv"]("b,f") * eri["abba_vovo"]("a,m,e,i");

        // H22_abb_aaa = +1.00 d_aa(a,e) <m,b||f,i>_abab  // flops: o2v4 = o2v4 | mem: o2v4 = o2v4
        H22_abb_aaa("a,b,i,e,f,m")  = -1.00 * Id["aa_vv"]("a,e") * eri["baab_vovo"]("b,m,f,i");

        // H22_abb_aaa += -1.00 d_aa(a,f) <m,b||e,i>_abab  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_aaa("a,b,i,e,f,m") += Id["aa_vv"]("a,f") * eri["baab_vovo"]("b,m,e,i");

        // H22_abb_abb = +1.00 d_aa(a,e) <m,b||f,i>_bbbb  // flops: o2v4 = o2v4 | mem: o2v4 = o2v4
        H22_abb_abb("a,b,i,e,f,m")  = -1.00 * Id["aa_vv"]("a,e") * eri["bbbb_vovo"]("b,m,f,i");

        // H22_abb_abb += +1.00 d_bb(i,m) <a,b||e,f>_abab  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_abb("a,b,i,e,f,m") += Id["bb_oo"]("i,m") * eri["abab_vvvv"]("a,b,e,f");

        // H22_abb_abb += -1.00 d_bb(b,f) <a,m||e,i>_abab  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_abb("a,b,i,e,f,m") -= Id["bb_vv"]("b,f") * eri["abab_vovo"]("a,m,e,i");

        // H22_bbb_aba = -1.00 d_bb(a,f) <m,b||e,i>_abab  // flops: o2v4 = o2v4 | mem: o2v4 = o2v4
        H22_bbb_aba("a,b,i,e,f,m")  = Id["bb_vv"]("a,f") * eri["baab_vovo"]("b,m,e,i");

        // H22_bbb_aba += +1.00 d_bb(b,f) <m,a||e,i>_abab  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_aba("a,b,i,e,f,m") -= Id["bb_vv"]("b,f") * eri["baab_vovo"]("a,m,e,i");

        // H22_bbb_bbb = +1.00 d_bb(a,e) <m,b||f,i>_bbbb  // flops: o2v4 = o2v4 | mem: o2v4 = o2v4
        H22_bbb_bbb("a,b,i,e,f,m")  = -1.00 * Id["bb_vv"]("a,e") * eri["bbbb_vovo"]("b,m,f,i");

        // H22_bbb_bbb += +1.00 d_bb(i,m) <a,b||e,f>_bbbb  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") += Id["bb_oo"]("i,m") * eri["bbbb_vvvv"]("a,b,e,f");

        // H22_bbb_bbb += -1.00 d_bb(a,f) <m,b||e,i>_bbbb  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") += Id["bb_vv"]("a,f") * eri["bbbb_vovo"]("b,m,e,i");

        // H22_bbb_bbb += -1.00 d_bb(b,e) <m,a||f,i>_bbbb  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") += Id["bb_vv"]("b,e") * eri["bbbb_vovo"]("a,m,f,i");

        // H22_bbb_bbb += +1.00 d_bb(b,f) <m,a||e,i>_bbbb  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") -= Id["bb_vv"]("b,f") * eri["bbbb_vovo"]("a,m,e,i");

        // H22_aaa_aaa += -1.00 <m,j||e,f>_aaaa t2_aaaa(a,b,i,j)  // flops: o2v4 += o3v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") -= eri["aaaa_oovv"]("m,j,e,f") * t2["aaaa_vvoo"]("a,b,i,j");

        // H22_aaa_abb += +1.00 <j,m||e,f>_abab t2_aaaa(a,b,i,j)  // flops: o2v4 += o3v4 | mem: o2v4 += o2v4
        H22_aaa_abb("a,b,i,e,f,m") += eri["abab_oovv"]("j,m,e,f") * t2["aaaa_vvoo"]("a,b,i,j");

        // H22_aba_aba += -1.00 <m,j||e,f>_abab t2_abab(a,b,i,j)  // flops: o2v4 += o3v4 | mem: o2v4 += o2v4
        H22_aba_aba("a,b,i,e,f,m") -= eri["abab_oovv"]("m,j,e,f") * t2["abab_vvoo"]("a,b,i,j");

        // H22_aba_bbb += -1.00 <m,j||e,f>_bbbb t2_abab(a,b,i,j)  // flops: o2v4 += o3v4 | mem: o2v4 += o2v4
        H22_aba_bbb("a,b,i,e,f,m") -= eri["bbbb_oovv"]("m,j,e,f") * t2["abab_vvoo"]("a,b,i,j");

        // H22_abb_aaa += +1.00 <m,j||e,f>_aaaa t2_abab(a,b,j,i)  // flops: o2v4 += o3v4 | mem: o2v4 += o2v4
        H22_abb_aaa("a,b,i,e,f,m") += eri["aaaa_oovv"]("m,j,e,f") * t2["abab_vvoo"]("a,b,j,i");

        // H22_abb_abb += -1.00 <j,m||e,f>_abab t2_abab(a,b,j,i)  // flops: o2v4 += o3v4 | mem: o2v4 += o2v4
        H22_abb_abb("a,b,i,e,f,m") -= eri["abab_oovv"]("j,m,e,f") * t2["abab_vvoo"]("a,b,j,i");

        // H22_bbb_aba += -1.00 <m,j||e,f>_abab t2_bbbb(a,b,i,j)  // flops: o2v4 += o3v4 | mem: o2v4 += o2v4
        H22_bbb_aba("a,b,i,e,f,m") -= eri["abab_oovv"]("m,j,e,f") * t2["bbbb_vvoo"]("a,b,i,j");

        // H22_bbb_bbb += -1.00 <m,j||e,f>_bbbb t2_bbbb(a,b,i,j)  // flops: o2v4 += o3v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") -= eri["bbbb_oovv"]("m,j,e,f") * t2["bbbb_vvoo"]("a,b,i,j");

        // H22_aaa_aaa += +0.50 d_aa(i,m) <k,j||e,f>_aaaa t2_aaaa(a,b,k,j)  // flops: o2v4 += o2v4 o2v4 | mem: o2v4 += o0v4 o2v4
        H22_aaa_aaa("a,b,i,e,f,m") += 0.50 * eri["aaaa_oovv"]("k,j,e,f") * t2["aaaa_vvoo"]("a,b,k,j") * Id["aa_oo"]("i,m");

        // H22_aba_aba += +1.00 d_aa(a,e) <m,j||c,f>_abab t2_abab(c,b,i,j)  // flops: o2v4 += o3v3 o2v4 | mem: o2v4 += o2v2 o2v4
        H22_aba_aba("a,b,i,e,f,m") += eri["abab_oovv"]("m,j,c,f") * t2["abab_vvoo"]("c,b,i,j") * Id["aa_vv"]("a,e");

        // H22_abb_abb += +1.00 d_bb(b,f) <j,m||e,c>_abab t2_abab(a,c,j,i)  // flops: o2v4 += o3v3 o2v4 | mem: o2v4 += o2v2 o2v4
        H22_abb_abb("a,b,i,e,f,m") += eri["abab_oovv"]("j,m,e,c") * t2["abab_vvoo"]("a,c,j,i") * Id["bb_vv"]("b,f");

        // H22_bbb_bbb += +0.50 d_bb(i,m) <k,j||e,f>_bbbb t2_bbbb(a,b,k,j)  // flops: o2v4 += o2v4 o2v4 | mem: o2v4 += o0v4 o2v4
        H22_bbb_bbb("a,b,i,e,f,m") += 0.50 * eri["bbbb_oovv"]("k,j,e,f") * t2["bbbb_vvoo"]("a,b,k,j") * Id["bb_oo"]("i,m");

        // tmps_[1_bbbb_vvvo](b,e,a,i) = 1.00 eri[bbbb_vovv](a,j,c,e) * t2[bbbb_vvoo](c,b,i,j) // flops: o1v3 = o2v4 | mem: o1v3 = o1v3
        tmps_["1_bbbb_vvvo"]("b,e,a,i")  = eri["bbbb_vovv"]("a,j,c,e") * t2["bbbb_vvoo"]("c,b,i,j");
        H21_bbb_b("a,b,i,e") -= tmps_["1_bbbb_vvvo"]("a,e,b,i");

        // H21_bbb_b += -1.00 P(a,b) <j,a||c,e>_bbbb t2_bbbb(c,b,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_bbb_b("a,b,i,e") += tmps_["1_bbbb_vvvo"]("b,e,a,i");
        tmps_["1_bbbb_vvvo"].~TArrayD();

        // tmps_[2_bbbb_vvvo](a,e,b,i) = 1.00 eri[baab_vovv](b,j,c,e) * t2[abab_vvoo](c,a,j,i) // flops: o1v3 = o2v4 | mem: o1v3 = o1v3
        tmps_["2_bbbb_vvvo"]("a,e,b,i")  = eri["baab_vovv"]("b,j,c,e") * t2["abab_vvoo"]("c,a,j,i");

        // H21_bbb_b += +1.00 P(a,b) <j,a||c,e>_abab t2_abab(c,b,j,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_bbb_b("a,b,i,e") -= tmps_["2_bbbb_vvvo"]("b,e,a,i");
        H21_bbb_b("a,b,i,e") += tmps_["2_bbbb_vvvo"]("a,e,b,i");
        tmps_["2_bbbb_vvvo"].~TArrayD();

        // tmps_[3_aaaa_vvvo](a,e,b,i) = 1.00 eri[aaaa_vovv](b,j,c,e) * t2[aaaa_vvoo](c,a,i,j) // flops: o1v3 = o2v4 | mem: o1v3 = o1v3
        tmps_["3_aaaa_vvvo"]("a,e,b,i")  = eri["aaaa_vovv"]("b,j,c,e") * t2["aaaa_vvoo"]("c,a,i,j");
        H21_aaa_a("a,b,i,e") -= tmps_["3_aaaa_vvvo"]("a,e,b,i");

        // H21_aaa_a += -1.00 P(a,b) <j,a||c,e>_aaaa t2_aaaa(c,b,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aaa_a("a,b,i,e") += tmps_["3_aaaa_vvvo"]("b,e,a,i");
        tmps_["3_aaaa_vvvo"].~TArrayD();

        // tmps_[4_aaaa_vvvo](b,e,a,i) = 1.00 eri[abab_vovv](a,j,e,c) * t2[abab_vvoo](b,c,i,j) // flops: o1v3 = o2v4 | mem: o1v3 = o1v3
        tmps_["4_aaaa_vvvo"]("b,e,a,i")  = eri["abab_vovv"]("a,j,e,c") * t2["abab_vvoo"]("b,c,i,j");
        H21_aaa_a("a,b,i,e") -= tmps_["4_aaaa_vvvo"]("a,e,b,i");

        // H21_aaa_a += +1.00 P(a,b) <a,j||e,c>_abab t2_abab(b,c,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aaa_a("a,b,i,e") += tmps_["4_aaaa_vvvo"]("b,e,a,i");
        tmps_["4_aaaa_vvvo"].~TArrayD();

        // tmps_[5_aabb_vvvv](a,e,b,f) = 1.00 eri[abab_oovv](k,j,e,f) * t2[abab_vvoo](a,b,k,j) // flops: o0v4 = o2v4 | mem: o0v4 = o0v4
        tmps_["5_aabb_vvvv"]("a,e,b,f")  = eri["abab_oovv"]("k,j,e,f") * t2["abab_vvoo"]("a,b,k,j");

        // H22_aba_aba += +0.50 d_aa(i,m) <k,j||e,f>_abab t2_abab(a,b,k,j) @ +0.50 d_aa(i,m) <j,k||e,f>_abab t2_abab(a,b,j,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_aba("a,b,i,e,f,m") += tmps_["5_aabb_vvvv"]("a,e,b,f") * Id["aa_oo"]("i,m");

        // H22_abb_abb += +0.50 d_bb(i,m) <k,j||e,f>_abab t2_abab(a,b,k,j) @ +0.50 d_bb(i,m) <j,k||e,f>_abab t2_abab(a,b,j,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_abb("a,b,i,e,f,m") += tmps_["5_aabb_vvvv"]("a,e,b,f") * Id["bb_oo"]("i,m");
        tmps_["5_aabb_vvvv"].~TArrayD();

        // tmps_[6_bbbb_vvoo](a,e,i,m) = 1.00 eri[abab_oovv](j,m,c,e) * t2[abab_vvoo](c,a,j,i) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["6_bbbb_vvoo"]("a,e,i,m")  = eri["abab_oovv"]("j,m,c,e") * t2["abab_vvoo"]("c,a,j,i");

        // H22_abb_abb += +1.00 d_aa(a,e) <j,m||c,f>_abab t2_abab(c,b,j,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_abb("a,b,i,e,f,m") += tmps_["6_bbbb_vvoo"]("b,f,i,m") * Id["aa_vv"]("a,e");

        // H22_bbb_bbb += -1.00 d_bb(b,e) <j,m||c,f>_abab t2_abab(c,a,j,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") -= tmps_["6_bbbb_vvoo"]("a,f,i,m") * Id["bb_vv"]("b,e");

        // H22_bbb_bbb += +1.00 d_bb(b,f) <j,m||c,e>_abab t2_abab(c,a,j,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") += tmps_["6_bbbb_vvoo"]("a,e,i,m") * Id["bb_vv"]("b,f");

        // H22_bbb_bbb += +1.00 d_bb(a,e) <j,m||c,f>_abab t2_abab(c,b,j,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") += tmps_["6_bbbb_vvoo"]("b,f,i,m") * Id["bb_vv"]("a,e");

        // H22_bbb_bbb += -1.00 d_bb(a,f) <j,m||c,e>_abab t2_abab(c,b,j,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") -= tmps_["6_bbbb_vvoo"]("b,e,i,m") * Id["bb_vv"]("a,f");
        tmps_["6_bbbb_vvoo"].~TArrayD();

        // tmps_[7_bbbb_vvoo](b,f,i,m) = 1.00 eri[bbbb_oovv](m,j,c,f) * t2[bbbb_vvoo](c,b,i,j) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["7_bbbb_vvoo"]("b,f,i,m")  = eri["bbbb_oovv"]("m,j,c,f") * t2["bbbb_vvoo"]("c,b,i,j");

        // H22_abb_abb += +1.00 d_aa(a,e) <m,j||c,f>_bbbb t2_bbbb(c,b,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_abb("a,b,i,e,f,m") += tmps_["7_bbbb_vvoo"]("b,f,i,m") * Id["aa_vv"]("a,e");

        // H22_bbb_bbb += -1.00 d_bb(b,e) <m,j||c,f>_bbbb t2_bbbb(c,a,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") -= tmps_["7_bbbb_vvoo"]("a,f,i,m") * Id["bb_vv"]("b,e");

        // H22_bbb_bbb += -1.00 d_bb(a,f) <m,j||c,e>_bbbb t2_bbbb(c,b,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") -= tmps_["7_bbbb_vvoo"]("b,e,i,m") * Id["bb_vv"]("a,f");

        // H22_bbb_bbb += +1.00 d_bb(a,e) <m,j||c,f>_bbbb t2_bbbb(c,b,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") += tmps_["7_bbbb_vvoo"]("b,f,i,m") * Id["bb_vv"]("a,e");

        // H22_bbb_bbb += +1.00 d_bb(b,f) <m,j||c,e>_bbbb t2_bbbb(c,a,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") += tmps_["7_bbbb_vvoo"]("a,e,i,m") * Id["bb_vv"]("b,f");
        tmps_["7_bbbb_vvoo"].~TArrayD();

        // tmps_[8_aaaa_vvoo](a,f,i,m) = 1.00 eri[abab_oovv](m,j,f,c) * t2[abab_vvoo](a,c,i,j) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["8_aaaa_vvoo"]("a,f,i,m")  = eri["abab_oovv"]("m,j,f,c") * t2["abab_vvoo"]("a,c,i,j");

        // H22_aaa_aaa += +1.00 d_aa(a,e) <m,j||f,c>_abab t2_abab(b,c,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") += tmps_["8_aaaa_vvoo"]("b,f,i,m") * Id["aa_vv"]("a,e");

        // H22_aaa_aaa += -1.00 d_aa(a,f) <m,j||e,c>_abab t2_abab(b,c,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") -= tmps_["8_aaaa_vvoo"]("b,e,i,m") * Id["aa_vv"]("a,f");

        // H22_aaa_aaa += -1.00 d_aa(b,e) <m,j||f,c>_abab t2_abab(a,c,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") -= tmps_["8_aaaa_vvoo"]("a,f,i,m") * Id["aa_vv"]("b,e");

        // H22_aaa_aaa += +1.00 d_aa(b,f) <m,j||e,c>_abab t2_abab(a,c,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") += tmps_["8_aaaa_vvoo"]("a,e,i,m") * Id["aa_vv"]("b,f");

        // H22_aba_aba += +1.00 d_bb(b,f) <m,j||e,c>_abab t2_abab(a,c,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_aba("a,b,i,e,f,m") += tmps_["8_aaaa_vvoo"]("a,e,i,m") * Id["bb_vv"]("b,f");
        tmps_["8_aaaa_vvoo"].~TArrayD();

        // tmps_[9_aaaa_vvoo](b,f,i,m) = 1.00 eri[aaaa_oovv](m,j,c,f) * t2[aaaa_vvoo](c,b,i,j) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["9_aaaa_vvoo"]("b,f,i,m")  = eri["aaaa_oovv"]("m,j,c,f") * t2["aaaa_vvoo"]("c,b,i,j");

        // H22_aaa_aaa += -1.00 d_aa(a,f) <m,j||c,e>_aaaa t2_aaaa(c,b,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") -= tmps_["9_aaaa_vvoo"]("b,e,i,m") * Id["aa_vv"]("a,f");

        // H22_aaa_aaa += -1.00 d_aa(b,e) <m,j||c,f>_aaaa t2_aaaa(c,a,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") -= tmps_["9_aaaa_vvoo"]("a,f,i,m") * Id["aa_vv"]("b,e");

        // H22_aaa_aaa += +1.00 d_aa(a,e) <m,j||c,f>_aaaa t2_aaaa(c,b,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") += tmps_["9_aaaa_vvoo"]("b,f,i,m") * Id["aa_vv"]("a,e");

        // H22_aaa_aaa += +1.00 d_aa(b,f) <m,j||c,e>_aaaa t2_aaaa(c,a,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") += tmps_["9_aaaa_vvoo"]("a,e,i,m") * Id["aa_vv"]("b,f");

        // H22_aba_aba += +1.00 d_bb(b,f) <m,j||c,e>_aaaa t2_aaaa(c,a,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_aba("a,b,i,e,f,m") += tmps_["9_aaaa_vvoo"]("a,e,i,m") * Id["bb_vv"]("b,f");
        tmps_["9_aaaa_vvoo"].~TArrayD();

        // tmps_[10_abab_vvoo](e,a,m,i) = 1.00 eri[aaaa_oovv](m,j,c,e) * t2[abab_vvoo](c,a,j,i) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["10_abab_vvoo"]("e,a,m,i")  = eri["aaaa_oovv"]("m,j,c,e") * t2["abab_vvoo"]("c,a,j,i");

        // H22_abb_aaa += -1.00 d_aa(a,e) <m,j||c,f>_aaaa t2_abab(c,b,j,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_aaa("a,b,i,e,f,m") -= tmps_["10_abab_vvoo"]("f,b,m,i") * Id["aa_vv"]("a,e");

        // H22_abb_aaa += +1.00 d_aa(a,f) <m,j||c,e>_aaaa t2_abab(c,b,j,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_aaa("a,b,i,e,f,m") += tmps_["10_abab_vvoo"]("e,b,m,i") * Id["aa_vv"]("a,f");

        // H22_bbb_aba += +1.00 d_bb(a,f) <m,j||c,e>_aaaa t2_abab(c,b,j,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_aba("a,b,i,e,f,m") += tmps_["10_abab_vvoo"]("e,b,m,i") * Id["bb_vv"]("a,f");

        // H22_bbb_aba += -1.00 d_bb(b,f) <m,j||c,e>_aaaa t2_abab(c,a,j,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_aba("a,b,i,e,f,m") -= tmps_["10_abab_vvoo"]("e,a,m,i") * Id["bb_vv"]("b,f");
        tmps_["10_abab_vvoo"].~TArrayD();

        // tmps_[11_abab_vvoo](e,b,m,i) = 1.00 eri[abab_oovv](m,j,e,c) * t2[bbbb_vvoo](c,b,i,j) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["11_abab_vvoo"]("e,b,m,i")  = eri["abab_oovv"]("m,j,e,c") * t2["bbbb_vvoo"]("c,b,i,j");

        // H22_abb_aaa += -1.00 d_aa(a,e) <m,j||f,c>_abab t2_bbbb(c,b,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_aaa("a,b,i,e,f,m") -= tmps_["11_abab_vvoo"]("f,b,m,i") * Id["aa_vv"]("a,e");

        // H22_abb_aaa += +1.00 d_aa(a,f) <m,j||e,c>_abab t2_bbbb(c,b,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_aaa("a,b,i,e,f,m") += tmps_["11_abab_vvoo"]("e,b,m,i") * Id["aa_vv"]("a,f");

        // H22_bbb_aba += -1.00 d_bb(b,f) <m,j||e,c>_abab t2_bbbb(c,a,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_aba("a,b,i,e,f,m") -= tmps_["11_abab_vvoo"]("e,a,m,i") * Id["bb_vv"]("b,f");

        // H22_bbb_aba += +1.00 d_bb(a,f) <m,j||e,c>_abab t2_bbbb(c,b,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_aba("a,b,i,e,f,m") += tmps_["11_abab_vvoo"]("e,b,m,i") * Id["bb_vv"]("a,f");
        tmps_["11_abab_vvoo"].~TArrayD();

        // tmps_[12_abab_vvoo](b,f,i,m) = 1.00 eri[bbbb_oovv](m,j,c,f) * t2[abab_vvoo](b,c,i,j) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["12_abab_vvoo"]("b,f,i,m")  = eri["bbbb_oovv"]("m,j,c,f") * t2["abab_vvoo"]("b,c,i,j");

        // H22_aaa_abb += -1.00 d_aa(a,e) <m,j||c,f>_bbbb t2_abab(b,c,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_abb("a,b,i,e,f,m") -= tmps_["12_abab_vvoo"]("b,f,i,m") * Id["aa_vv"]("a,e");

        // H22_aaa_abb += +1.00 d_aa(b,e) <m,j||c,f>_bbbb t2_abab(a,c,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_abb("a,b,i,e,f,m") += tmps_["12_abab_vvoo"]("a,f,i,m") * Id["aa_vv"]("b,e");

        // H22_aba_bbb += +1.00 d_bb(b,e) <m,j||c,f>_bbbb t2_abab(a,c,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_bbb("a,b,i,e,f,m") += tmps_["12_abab_vvoo"]("a,f,i,m") * Id["bb_vv"]("b,e");

        // H22_aba_bbb += -1.00 d_bb(b,f) <m,j||c,e>_bbbb t2_abab(a,c,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_bbb("a,b,i,e,f,m") -= tmps_["12_abab_vvoo"]("a,e,i,m") * Id["bb_vv"]("b,f");
        tmps_["12_abab_vvoo"].~TArrayD();

        // tmps_[13_abab_vvoo](b,f,i,m) = 1.00 eri[abab_oovv](j,m,c,f) * t2[aaaa_vvoo](c,b,i,j) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["13_abab_vvoo"]("b,f,i,m")  = eri["abab_oovv"]("j,m,c,f") * t2["aaaa_vvoo"]("c,b,i,j");

        // H22_aaa_abb += -1.00 d_aa(a,e) <j,m||c,f>_abab t2_aaaa(c,b,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_abb("a,b,i,e,f,m") -= tmps_["13_abab_vvoo"]("b,f,i,m") * Id["aa_vv"]("a,e");

        // H22_aaa_abb += +1.00 d_aa(b,e) <j,m||c,f>_abab t2_aaaa(c,a,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_abb("a,b,i,e,f,m") += tmps_["13_abab_vvoo"]("a,f,i,m") * Id["aa_vv"]("b,e");

        // H22_aba_bbb += -1.00 d_bb(b,f) <j,m||c,e>_abab t2_aaaa(c,a,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_bbb("a,b,i,e,f,m") -= tmps_["13_abab_vvoo"]("a,e,i,m") * Id["bb_vv"]("b,f");

        // H22_aba_bbb += +1.00 d_bb(b,e) <j,m||c,f>_abab t2_aaaa(c,a,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_bbb("a,b,i,e,f,m") += tmps_["13_abab_vvoo"]("a,f,i,m") * Id["bb_vv"]("b,e");
        tmps_["13_abab_vvoo"].~TArrayD();

        // tmps_[14_bb_vv](a,e) = 0.50 eri[bbbb_oovv](k,j,c,e) * t2[bbbb_vvoo](c,a,k,j) // flops: o0v2 = o2v3 | mem: o0v2 = o0v2
        tmps_["14_bb_vv"]("a,e")  = 0.50 * eri["bbbb_oovv"]("k,j,c,e") * t2["bbbb_vvoo"]("c,a,k,j");

        // H11_b_b += -0.50 <j,i||b,e>_bbbb t2_bbbb(b,a,j,i)  // flops: o0v2 += o0v2 | mem: o0v2 += o0v2
        H11_b_b("a,e") -= tmps_["14_bb_vv"]("a,e");

        // tmps_[30_bbbb_vvoo](f,b,m,i) = 1.00 Id[bb_oo](i,m) * Id[bb_vv](b,f) // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
        tmps_["30_bbbb_vvoo"]("f,b,m,i")  = Id["bb_oo"]("i,m") * Id["bb_vv"]("b,f");

        // H22_bbb_bbb += +0.50 d_bb(b,e) d_bb(i,m) <k,j||c,f>_bbbb t2_bbbb(c,a,k,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") += tmps_["14_bb_vv"]("a,f") * tmps_["30_bbbb_vvoo"]("e,b,m,i");

        // H22_bbb_bbb += -0.50 d_bb(a,e) d_bb(i,m) <k,j||c,f>_bbbb t2_bbbb(c,b,k,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") -= tmps_["14_bb_vv"]("b,f") * tmps_["30_bbbb_vvoo"]("e,a,m,i");

        // H22_bbb_bbb += +0.50 d_bb(a,f) d_bb(i,m) <k,j||c,e>_bbbb t2_bbbb(c,b,k,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") += tmps_["14_bb_vv"]("b,e") * tmps_["30_bbbb_vvoo"]("f,a,m,i");

        // H22_bbb_bbb += -0.50 d_bb(b,f) d_bb(i,m) <k,j||c,e>_bbbb t2_bbbb(c,a,k,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") -= tmps_["14_bb_vv"]("a,e") * tmps_["30_bbbb_vvoo"]("f,b,m,i");

        // tmps_[31_aaaa_vvoo](e,b,m,i) = 1.00 Id[aa_oo](i,m) * Id[aa_vv](b,e) // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
        tmps_["31_aaaa_vvoo"]("e,b,m,i")  = Id["aa_oo"]("i,m") * Id["aa_vv"]("b,e");

        // H22_aba_aba += -0.50 d_aa(a,e) d_aa(i,m) <k,j||c,f>_bbbb t2_bbbb(c,b,k,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_aba("a,b,i,e,f,m") -= tmps_["14_bb_vv"]("b,f") * tmps_["31_aaaa_vvoo"]("e,a,m,i");

        // tmps_[37_aabb_vvoo](e,a,m,i) = 1.00 Id[aa_vv](a,e) * Id[bb_oo](i,m) // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
        tmps_["37_aabb_vvoo"]("e,a,m,i")  = Id["aa_vv"]("a,e") * Id["bb_oo"]("i,m");

        // H22_abb_abb += -0.50 d_aa(a,e) d_bb(i,m) <k,j||c,f>_bbbb t2_bbbb(c,b,k,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_abb("a,b,i,e,f,m") -= tmps_["14_bb_vv"]("b,f") * tmps_["37_aabb_vvoo"]("e,a,m,i");
        tmps_["14_bb_vv"].~TArrayD();

        // tmps_[15_bb_vv](a,f) = 1.00 eri[abab_oovv](k,j,c,f) * t2[abab_vvoo](c,a,k,j) // flops: o0v2 = o2v3 | mem: o0v2 = o0v2
        tmps_["15_bb_vv"]("a,f")  = eri["abab_oovv"]("k,j,c,f") * t2["abab_vvoo"]("c,a,k,j");

        // H11_b_b += -0.50 <j,i||b,e>_abab t2_abab(b,a,j,i) @ -0.50 <i,j||b,e>_abab t2_abab(b,a,i,j)  // flops: o0v2 += o0v2 | mem: o0v2 += o0v2
        H11_b_b("a,e") -= tmps_["15_bb_vv"]("a,e");

        // H22_bbb_bbb += -0.50 d_bb(a,e) d_bb(i,m) <k,j||c,f>_abab t2_abab(c,b,k,j) @ -0.50 d_bb(a,e) d_bb(i,m) <j,k||c,f>_abab t2_abab(c,b,j,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") -= tmps_["15_bb_vv"]("b,f") * tmps_["30_bbbb_vvoo"]("e,a,m,i");

        // H22_bbb_bbb += +0.50 d_bb(a,f) d_bb(i,m) <k,j||c,e>_abab t2_abab(c,b,k,j) @ +0.50 d_bb(a,f) d_bb(i,m) <j,k||c,e>_abab t2_abab(c,b,j,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") += tmps_["15_bb_vv"]("b,e") * tmps_["30_bbbb_vvoo"]("f,a,m,i");

        // H22_bbb_bbb += -0.50 d_bb(b,f) d_bb(i,m) <k,j||c,e>_abab t2_abab(c,a,k,j) @ -0.50 d_bb(b,f) d_bb(i,m) <j,k||c,e>_abab t2_abab(c,a,j,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") -= tmps_["15_bb_vv"]("a,e") * tmps_["30_bbbb_vvoo"]("f,b,m,i");

        // H22_bbb_bbb += +0.50 d_bb(b,e) d_bb(i,m) <k,j||c,f>_abab t2_abab(c,a,k,j) @ +0.50 d_bb(b,e) d_bb(i,m) <j,k||c,f>_abab t2_abab(c,a,j,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") += tmps_["15_bb_vv"]("a,f") * tmps_["30_bbbb_vvoo"]("e,b,m,i");

        // H22_aba_aba += -0.50 d_aa(a,e) d_aa(i,m) <k,j||c,f>_abab t2_abab(c,b,k,j) @ -0.50 d_aa(a,e) d_aa(i,m) <j,k||c,f>_abab t2_abab(c,b,j,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_aba("a,b,i,e,f,m") -= tmps_["15_bb_vv"]("b,f") * tmps_["31_aaaa_vvoo"]("e,a,m,i");

        // H22_abb_abb += -0.50 d_aa(a,e) d_bb(i,m) <k,j||c,f>_abab t2_abab(c,b,k,j) @ -0.50 d_aa(a,e) d_bb(i,m) <j,k||c,f>_abab t2_abab(c,b,j,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_abb("a,b,i,e,f,m") -= tmps_["15_bb_vv"]("b,f") * tmps_["37_aabb_vvoo"]("e,a,m,i");
        tmps_["15_bb_vv"].~TArrayD();

        // tmps_[16_aa_vv](b,e) = 1.00 eri[abab_oovv](k,j,e,c) * t2[abab_vvoo](b,c,k,j) // flops: o0v2 = o2v3 | mem: o0v2 = o0v2
        tmps_["16_aa_vv"]("b,e")  = eri["abab_oovv"]("k,j,e,c") * t2["abab_vvoo"]("b,c,k,j");

        // H11_a_a += -0.50 <j,i||e,b>_abab t2_abab(a,b,j,i) @ -0.50 <i,j||e,b>_abab t2_abab(a,b,i,j)  // flops: o0v2 += o0v2 | mem: o0v2 += o0v2
        H11_a_a("a,e") -= tmps_["16_aa_vv"]("a,e");

        // H22_abb_abb += -0.50 d_bb(b,f) d_bb(i,m) <k,j||e,c>_abab t2_abab(a,c,k,j) @ -0.50 d_bb(b,f) d_bb(i,m) <j,k||e,c>_abab t2_abab(a,c,j,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_abb("a,b,i,e,f,m") -= tmps_["16_aa_vv"]("a,e") * tmps_["30_bbbb_vvoo"]("f,b,m,i");

        // H22_aaa_aaa += -0.50 d_aa(a,e) d_aa(i,m) <k,j||f,c>_abab t2_abab(b,c,k,j) @ -0.50 d_aa(a,e) d_aa(i,m) <j,k||f,c>_abab t2_abab(b,c,j,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") -= tmps_["16_aa_vv"]("b,f") * tmps_["31_aaaa_vvoo"]("e,a,m,i");

        // H22_aaa_aaa += +0.50 d_aa(a,f) d_aa(i,m) <k,j||e,c>_abab t2_abab(b,c,k,j) @ +0.50 d_aa(a,f) d_aa(i,m) <j,k||e,c>_abab t2_abab(b,c,j,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") += tmps_["16_aa_vv"]("b,e") * tmps_["31_aaaa_vvoo"]("f,a,m,i");

        // H22_aaa_aaa += +0.50 d_aa(b,e) d_aa(i,m) <k,j||f,c>_abab t2_abab(a,c,k,j) @ +0.50 d_aa(b,e) d_aa(i,m) <j,k||f,c>_abab t2_abab(a,c,j,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") += tmps_["16_aa_vv"]("a,f") * tmps_["31_aaaa_vvoo"]("e,b,m,i");

        // H22_aaa_aaa += -0.50 d_aa(b,f) d_aa(i,m) <k,j||e,c>_abab t2_abab(a,c,k,j) @ -0.50 d_aa(b,f) d_aa(i,m) <j,k||e,c>_abab t2_abab(a,c,j,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") -= tmps_["16_aa_vv"]("a,e") * tmps_["31_aaaa_vvoo"]("f,b,m,i");

        // tmps_[38_bbaa_vvoo](f,b,m,i) = 1.00 Id[aa_oo](i,m) * Id[bb_vv](b,f) // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
        tmps_["38_bbaa_vvoo"]("f,b,m,i")  = Id["aa_oo"]("i,m") * Id["bb_vv"]("b,f");

        // H22_aba_aba += -0.50 d_bb(b,f) d_aa(i,m) <k,j||e,c>_abab t2_abab(a,c,k,j) @ -0.50 d_bb(b,f) d_aa(i,m) <j,k||e,c>_abab t2_abab(a,c,j,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_aba("a,b,i,e,f,m") -= tmps_["16_aa_vv"]("a,e") * tmps_["38_bbaa_vvoo"]("f,b,m,i");
        tmps_["16_aa_vv"].~TArrayD();

        // tmps_[17_aa_vv](b,e) = 0.50 eri[aaaa_oovv](k,j,c,e) * t2[aaaa_vvoo](c,b,k,j) // flops: o0v2 = o2v3 | mem: o0v2 = o0v2
        tmps_["17_aa_vv"]("b,e")  = 0.50 * eri["aaaa_oovv"]("k,j,c,e") * t2["aaaa_vvoo"]("c,b,k,j");

        // H11_a_a += -0.50 <j,i||b,e>_aaaa t2_aaaa(b,a,j,i)  // flops: o0v2 += o0v2 | mem: o0v2 += o0v2
        H11_a_a("a,e") -= tmps_["17_aa_vv"]("a,e");

        // H22_abb_abb += -0.50 d_bb(b,f) d_bb(i,m) <k,j||c,e>_aaaa t2_aaaa(c,a,k,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_abb("a,b,i,e,f,m") -= tmps_["17_aa_vv"]("a,e") * tmps_["30_bbbb_vvoo"]("f,b,m,i");

        // H22_aaa_aaa += -0.50 d_aa(b,f) d_aa(i,m) <k,j||c,e>_aaaa t2_aaaa(c,a,k,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") -= tmps_["17_aa_vv"]("a,e") * tmps_["31_aaaa_vvoo"]("f,b,m,i");

        // H22_aaa_aaa += +0.50 d_aa(b,e) d_aa(i,m) <k,j||c,f>_aaaa t2_aaaa(c,a,k,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") += tmps_["17_aa_vv"]("a,f") * tmps_["31_aaaa_vvoo"]("e,b,m,i");

        // H22_aaa_aaa += +0.50 d_aa(a,f) d_aa(i,m) <k,j||c,e>_aaaa t2_aaaa(c,b,k,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") += tmps_["17_aa_vv"]("b,e") * tmps_["31_aaaa_vvoo"]("f,a,m,i");

        // H22_aaa_aaa += -0.50 d_aa(a,e) d_aa(i,m) <k,j||c,f>_aaaa t2_aaaa(c,b,k,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") -= tmps_["17_aa_vv"]("b,f") * tmps_["31_aaaa_vvoo"]("e,a,m,i");

        // H22_aba_aba += -0.50 d_bb(b,f) d_aa(i,m) <k,j||c,e>_aaaa t2_aaaa(c,a,k,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_aba("a,b,i,e,f,m") -= tmps_["17_aa_vv"]("a,e") * tmps_["38_bbaa_vvoo"]("f,b,m,i");
        tmps_["17_aa_vv"].~TArrayD();

        // tmps_[18_bb_vo](b,i) = 0.50 eri[bbbb_vovv](b,j,c,d) * t2[bbbb_vvoo](c,d,i,j) // flops: o1v1 = o2v3 | mem: o1v1 = o1v1
        tmps_["18_bb_vo"]("b,i")  = 0.50 * eri["bbbb_vovv"]("b,j,c,d") * t2["bbbb_vvoo"]("c,d,i,j");

        // H21_abb_a += -0.50 d_aa(a,e) <j,b||c,d>_bbbb t2_bbbb(c,d,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_abb_a("a,b,i,e") += tmps_["18_bb_vo"]("b,i") * Id["aa_vv"]("a,e");

        // H21_bbb_b += +0.50 d_bb(b,e) <j,a||c,d>_bbbb t2_bbbb(c,d,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_bbb_b("a,b,i,e") -= tmps_["18_bb_vo"]("a,i") * Id["bb_vv"]("b,e");

        // H21_bbb_b += -0.50 d_bb(a,e) <j,b||c,d>_bbbb t2_bbbb(c,d,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_bbb_b("a,b,i,e") += tmps_["18_bb_vo"]("b,i") * Id["bb_vv"]("a,e");
        tmps_["18_bb_vo"].~TArrayD();

        // tmps_[19_bb_vo](b,i) = 1.00 eri[baab_vovv](b,j,c,d) * t2[abab_vvoo](c,d,j,i) // flops: o1v1 = o2v3 | mem: o1v1 = o1v1
        tmps_["19_bb_vo"]("b,i")  = eri["baab_vovv"]("b,j,c,d") * t2["abab_vvoo"]("c,d,j,i");

        // H21_abb_a += +0.50 d_aa(a,e) <j,b||c,d>_abab t2_abab(c,d,j,i) @ +0.50 d_aa(a,e) <j,b||d,c>_abab t2_abab(d,c,j,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_abb_a("a,b,i,e") -= tmps_["19_bb_vo"]("b,i") * Id["aa_vv"]("a,e");

        // H21_bbb_b += +0.50 d_bb(a,e) <j,b||c,d>_abab t2_abab(c,d,j,i) @ +0.50 d_bb(a,e) <j,b||d,c>_abab t2_abab(d,c,j,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_bbb_b("a,b,i,e") -= tmps_["19_bb_vo"]("b,i") * Id["bb_vv"]("a,e");

        // H21_bbb_b += -0.50 d_bb(b,e) <j,a||c,d>_abab t2_abab(c,d,j,i) @ -0.50 d_bb(b,e) <j,a||d,c>_abab t2_abab(d,c,j,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_bbb_b("a,b,i,e") += tmps_["19_bb_vo"]("a,i") * Id["bb_vv"]("b,e");
        tmps_["19_bb_vo"].~TArrayD();

        // tmps_[20_aa_vo](a,i) = 1.00 eri[abab_vovv](a,j,c,d) * t2[abab_vvoo](c,d,i,j) // flops: o1v1 = o2v3 | mem: o1v1 = o1v1
        tmps_["20_aa_vo"]("a,i")  = eri["abab_vovv"]("a,j,c,d") * t2["abab_vvoo"]("c,d,i,j");

        // H21_aaa_a += +0.50 d_aa(a,e) <b,j||c,d>_abab t2_abab(c,d,i,j) @ +0.50 d_aa(a,e) <b,j||d,c>_abab t2_abab(d,c,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aaa_a("a,b,i,e") += tmps_["20_aa_vo"]("b,i") * Id["aa_vv"]("a,e");

        // H21_aaa_a += -0.50 d_aa(b,e) <a,j||c,d>_abab t2_abab(c,d,i,j) @ -0.50 d_aa(b,e) <a,j||d,c>_abab t2_abab(d,c,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aaa_a("a,b,i,e") -= tmps_["20_aa_vo"]("a,i") * Id["aa_vv"]("b,e");

        // H21_aba_b += -0.50 d_bb(b,e) <a,j||c,d>_abab t2_abab(c,d,i,j) @ -0.50 d_bb(b,e) <a,j||d,c>_abab t2_abab(d,c,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aba_b("a,b,i,e") -= tmps_["20_aa_vo"]("a,i") * Id["bb_vv"]("b,e");
        tmps_["20_aa_vo"].~TArrayD();

        // tmps_[21_aa_vo](b,i) = 0.50 eri[aaaa_vovv](b,j,c,d) * t2[aaaa_vvoo](c,d,i,j) // flops: o1v1 = o2v3 | mem: o1v1 = o1v1
        tmps_["21_aa_vo"]("b,i")  = 0.50 * eri["aaaa_vovv"]("b,j,c,d") * t2["aaaa_vvoo"]("c,d,i,j");

        // H21_aaa_a += -0.50 d_aa(a,e) <j,b||c,d>_aaaa t2_aaaa(c,d,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aaa_a("a,b,i,e") += tmps_["21_aa_vo"]("b,i") * Id["aa_vv"]("a,e");

        // H21_aaa_a += +0.50 d_aa(b,e) <j,a||c,d>_aaaa t2_aaaa(c,d,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aaa_a("a,b,i,e") -= tmps_["21_aa_vo"]("a,i") * Id["aa_vv"]("b,e");

        // H21_aba_b += +0.50 d_bb(b,e) <j,a||c,d>_aaaa t2_aaaa(c,d,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aba_b("a,b,i,e") -= tmps_["21_aa_vo"]("a,i") * Id["bb_vv"]("b,e");
        tmps_["21_aa_vo"].~TArrayD();

        // tmps_[22_bb_oo](i,m) = 0.50 eri[bbbb_oovv](m,j,c,d) * t2[bbbb_vvoo](c,d,i,j) // flops: o2v0 = o3v2 | mem: o2v0 = o2v0
        tmps_["22_bb_oo"]("i,m")  = 0.50 * eri["bbbb_oovv"]("m,j,c,d") * t2["bbbb_vvoo"]("c,d,i,j");

        // tmps_[41_bbbb_vvoo](f,b,m,i) = 1.00 Id[bb_vv](b,f) * eri[bbbb_oovv](m,j,c,d) * t2[bbbb_vvoo](c,d,i,j) // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
        tmps_["41_bbbb_vvoo"]("f,b,m,i")  = Id["bb_vv"]("b,f") * tmps_["22_bb_oo"]("i,m");
        tmps_["22_bb_oo"].~TArrayD();

        // H22_abb_abb += -0.50 d_aa(a,e) d_bb(b,f) <m,j||c,d>_bbbb t2_bbbb(c,d,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_abb("a,b,i,e,f,m") -= Id["aa_vv"]("a,e") * tmps_["41_bbbb_vvoo"]("f,b,m,i");

        // H22_bbb_bbb += +0.50 d_bb(b,e) d_bb(a,f) <m,j||c,d>_bbbb t2_bbbb(c,d,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") += tmps_["41_bbbb_vvoo"]("e,b,m,i") * Id["bb_vv"]("a,f");

        // H22_bbb_bbb += -0.50 d_bb(a,e) d_bb(b,f) <m,j||c,d>_bbbb t2_bbbb(c,d,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") -= tmps_["41_bbbb_vvoo"]("e,a,m,i") * Id["bb_vv"]("b,f");
        tmps_["41_bbbb_vvoo"].~TArrayD();

        // tmps_[23_aa_oo](i,m) = 1.00 eri[abab_oovv](m,j,c,d) * t2[abab_vvoo](c,d,i,j) // flops: o2v0 = o3v2 | mem: o2v0 = o2v0
        tmps_["23_aa_oo"]("i,m")  = eri["abab_oovv"]("m,j,c,d") * t2["abab_vvoo"]("c,d,i,j");

        // tmps_[43_aaaa_vvoo](e,a,m,i) = 1.00 Id[aa_vv](a,e) * eri[abab_oovv](m,j,c,d) * t2[abab_vvoo](c,d,i,j) // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
        tmps_["43_aaaa_vvoo"]("e,a,m,i")  = Id["aa_vv"]("a,e") * tmps_["23_aa_oo"]("i,m");
        tmps_["23_aa_oo"].~TArrayD();

        // H22_aaa_aaa += -0.50 d_aa(a,e) d_aa(b,f) <m,j||c,d>_abab t2_abab(c,d,i,j) @ -0.50 d_aa(a,e) d_aa(b,f) <m,j||d,c>_abab t2_abab(d,c,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") -= tmps_["43_aaaa_vvoo"]("e,a,m,i") * Id["aa_vv"]("b,f");

        // H22_aaa_aaa += +0.50 d_aa(b,e) d_aa(a,f) <m,j||c,d>_abab t2_abab(c,d,i,j) @ +0.50 d_aa(b,e) d_aa(a,f) <m,j||d,c>_abab t2_abab(d,c,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") += tmps_["43_aaaa_vvoo"]("e,b,m,i") * Id["aa_vv"]("a,f");

        // H22_aba_aba += -0.50 d_aa(a,e) d_bb(b,f) <m,j||c,d>_abab t2_abab(c,d,i,j) @ -0.50 d_aa(a,e) d_bb(b,f) <m,j||d,c>_abab t2_abab(d,c,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_aba("a,b,i,e,f,m") -= tmps_["43_aaaa_vvoo"]("e,a,m,i") * Id["bb_vv"]("b,f");
        tmps_["43_aaaa_vvoo"].~TArrayD();

        // tmps_[24_aa_oo](i,m) = 0.50 eri[aaaa_oovv](m,j,c,d) * t2[aaaa_vvoo](c,d,i,j) // flops: o2v0 = o3v2 | mem: o2v0 = o2v0
        tmps_["24_aa_oo"]("i,m")  = 0.50 * eri["aaaa_oovv"]("m,j,c,d") * t2["aaaa_vvoo"]("c,d,i,j");

        // tmps_[42_aaaa_vvoo](e,b,m,i) = 1.00 Id[aa_vv](b,e) * eri[aaaa_oovv](m,j,c,d) * t2[aaaa_vvoo](c,d,i,j) // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
        tmps_["42_aaaa_vvoo"]("e,b,m,i")  = Id["aa_vv"]("b,e") * tmps_["24_aa_oo"]("i,m");
        tmps_["24_aa_oo"].~TArrayD();

        // H22_aaa_aaa += +0.50 d_aa(b,e) d_aa(a,f) <m,j||c,d>_aaaa t2_aaaa(c,d,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") += tmps_["42_aaaa_vvoo"]("e,b,m,i") * Id["aa_vv"]("a,f");

        // H22_aaa_aaa += -0.50 d_aa(a,e) d_aa(b,f) <m,j||c,d>_aaaa t2_aaaa(c,d,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") -= tmps_["42_aaaa_vvoo"]("e,a,m,i") * Id["aa_vv"]("b,f");

        // H22_aba_aba += -0.50 d_aa(a,e) d_bb(b,f) <m,j||c,d>_aaaa t2_aaaa(c,d,i,j)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_aba("a,b,i,e,f,m") -= tmps_["42_aaaa_vvoo"]("e,a,m,i") * Id["bb_vv"]("b,f");
        tmps_["42_aaaa_vvoo"].~TArrayD();

        // tmps_[25_bb_oo](i,m) = 1.00 eri[abab_oovv](j,m,c,d) * t2[abab_vvoo](c,d,j,i) // flops: o2v0 = o3v2 | mem: o2v0 = o2v0
        tmps_["25_bb_oo"]("i,m")  = eri["abab_oovv"]("j,m,c,d") * t2["abab_vvoo"]("c,d,j,i");

        // tmps_[40_bbbb_vvoo](f,b,m,i) = 1.00 Id[bb_vv](b,f) * eri[abab_oovv](j,m,c,d) * t2[abab_vvoo](c,d,j,i) // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
        tmps_["40_bbbb_vvoo"]("f,b,m,i")  = Id["bb_vv"]("b,f") * tmps_["25_bb_oo"]("i,m");
        tmps_["25_bb_oo"].~TArrayD();

        // H22_abb_abb += -0.50 d_aa(a,e) d_bb(b,f) <j,m||c,d>_abab t2_abab(c,d,j,i) @ -0.50 d_aa(a,e) d_bb(b,f) <j,m||d,c>_abab t2_abab(d,c,j,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_abb("a,b,i,e,f,m") -= Id["aa_vv"]("a,e") * tmps_["40_bbbb_vvoo"]("f,b,m,i");

        // H22_bbb_bbb += +0.50 d_bb(b,e) d_bb(a,f) <j,m||c,d>_abab t2_abab(c,d,j,i) @ +0.50 d_bb(b,e) d_bb(a,f) <j,m||d,c>_abab t2_abab(d,c,j,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") += tmps_["40_bbbb_vvoo"]("e,b,m,i") * Id["bb_vv"]("a,f");

        // H22_bbb_bbb += -0.50 d_bb(a,e) d_bb(b,f) <j,m||c,d>_abab t2_abab(c,d,j,i) @ -0.50 d_bb(a,e) d_bb(b,f) <j,m||d,c>_abab t2_abab(d,c,j,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") -= tmps_["40_bbbb_vvoo"]("e,a,m,i") * Id["bb_vv"]("b,f");
        tmps_["40_bbbb_vvoo"].~TArrayD();

        // tmps_[26_bb_vo](a,i) = 1.00 eri[abab_oovo](k,j,c,i) * t2[abab_vvoo](c,a,k,j) // flops: o1v1 = o3v2 | mem: o1v1 = o1v1
        tmps_["26_bb_vo"]("a,i")  = eri["abab_oovo"]("k,j,c,i") * t2["abab_vvoo"]("c,a,k,j");

        // H21_abb_a += -0.50 d_aa(a,e) <k,j||c,i>_abab t2_abab(c,b,k,j) @ -0.50 d_aa(a,e) <j,k||c,i>_abab t2_abab(c,b,j,k)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_abb_a("a,b,i,e") -= tmps_["26_bb_vo"]("b,i") * Id["aa_vv"]("a,e");

        // H21_bbb_b += +0.50 d_bb(b,e) <k,j||c,i>_abab t2_abab(c,a,k,j) @ +0.50 d_bb(b,e) <j,k||c,i>_abab t2_abab(c,a,j,k)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_bbb_b("a,b,i,e") += tmps_["26_bb_vo"]("a,i") * Id["bb_vv"]("b,e");

        // H21_bbb_b += -0.50 d_bb(a,e) <k,j||c,i>_abab t2_abab(c,b,k,j) @ -0.50 d_bb(a,e) <j,k||c,i>_abab t2_abab(c,b,j,k)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_bbb_b("a,b,i,e") -= tmps_["26_bb_vo"]("b,i") * Id["bb_vv"]("a,e");
        tmps_["26_bb_vo"].~TArrayD();

        // tmps_[27_bb_vo](b,i) = 0.50 eri[bbbb_oovo](k,j,c,i) * t2[bbbb_vvoo](c,b,k,j) // flops: o1v1 = o3v2 | mem: o1v1 = o1v1
        tmps_["27_bb_vo"]("b,i")  = 0.50 * eri["bbbb_oovo"]("k,j,c,i") * t2["bbbb_vvoo"]("c,b,k,j");

        // H21_abb_a += -0.50 d_aa(a,e) <k,j||c,i>_bbbb t2_bbbb(c,b,k,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_abb_a("a,b,i,e") -= tmps_["27_bb_vo"]("b,i") * Id["aa_vv"]("a,e");

        // H21_bbb_b += +0.50 d_bb(b,e) <k,j||c,i>_bbbb t2_bbbb(c,a,k,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_bbb_b("a,b,i,e") += tmps_["27_bb_vo"]("a,i") * Id["bb_vv"]("b,e");

        // H21_bbb_b += -0.50 d_bb(a,e) <k,j||c,i>_bbbb t2_bbbb(c,b,k,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_bbb_b("a,b,i,e") -= tmps_["27_bb_vo"]("b,i") * Id["bb_vv"]("a,e");
        tmps_["27_bb_vo"].~TArrayD();

        // tmps_[28_aa_vo](a,i) = 0.50 eri[aaaa_oovo](k,j,c,i) * t2[aaaa_vvoo](c,a,k,j) // flops: o1v1 = o3v2 | mem: o1v1 = o1v1
        tmps_["28_aa_vo"]("a,i")  = 0.50 * eri["aaaa_oovo"]("k,j,c,i") * t2["aaaa_vvoo"]("c,a,k,j");

        // H21_aaa_a += +0.50 d_aa(b,e) <k,j||c,i>_aaaa t2_aaaa(c,a,k,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aaa_a("a,b,i,e") += tmps_["28_aa_vo"]("a,i") * Id["aa_vv"]("b,e");

        // H21_aaa_a += -0.50 d_aa(a,e) <k,j||c,i>_aaaa t2_aaaa(c,b,k,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aaa_a("a,b,i,e") -= tmps_["28_aa_vo"]("b,i") * Id["aa_vv"]("a,e");

        // H21_aba_b += +0.50 d_bb(b,e) <k,j||c,i>_aaaa t2_aaaa(c,a,k,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aba_b("a,b,i,e") += tmps_["28_aa_vo"]("a,i") * Id["bb_vv"]("b,e");
        tmps_["28_aa_vo"].~TArrayD();

        // tmps_[29_aa_vo](b,i) = 1.00 eri[abba_oovo](k,j,c,i) * t2[abab_vvoo](b,c,k,j) // flops: o1v1 = o3v2 | mem: o1v1 = o1v1
        tmps_["29_aa_vo"]("b,i")  = eri["abba_oovo"]("k,j,c,i") * t2["abab_vvoo"]("b,c,k,j");

        // H21_aaa_a += -0.50 d_aa(a,e) <k,j||i,c>_abab t2_abab(b,c,k,j) @ -0.50 d_aa(a,e) <j,k||i,c>_abab t2_abab(b,c,j,k)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aaa_a("a,b,i,e") += tmps_["29_aa_vo"]("b,i") * Id["aa_vv"]("a,e");

        // H21_aaa_a += +0.50 d_aa(b,e) <k,j||i,c>_abab t2_abab(a,c,k,j) @ +0.50 d_aa(b,e) <j,k||i,c>_abab t2_abab(a,c,j,k)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aaa_a("a,b,i,e") -= tmps_["29_aa_vo"]("a,i") * Id["aa_vv"]("b,e");

        // H21_aba_b += +0.50 d_bb(b,e) <k,j||i,c>_abab t2_abab(a,c,k,j) @ +0.50 d_bb(b,e) <j,k||i,c>_abab t2_abab(a,c,j,k)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aba_b("a,b,i,e") -= tmps_["29_aa_vo"]("a,i") * Id["bb_vv"]("b,e");
        tmps_["29_aa_vo"].~TArrayD();

        // H22_abb_abb += +1.00 d_bb(b,f) d_bb(i,m) f_aa(a,e)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_abb("a,b,i,e,f,m") += tmps_["30_bbbb_vvoo"]("f,b,m,i") * F["aa_vv"]("a,e");

        // H22_bbb_bbb += +1.00 d_bb(a,e) d_bb(i,m) f_bb(b,f)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") += tmps_["30_bbbb_vvoo"]("e,a,m,i") * F["bb_vv"]("b,f");

        // H22_bbb_bbb += -1.00 d_bb(a,f) d_bb(i,m) f_bb(b,e)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") -= tmps_["30_bbbb_vvoo"]("f,a,m,i") * F["bb_vv"]("b,e");

        // H22_bbb_bbb += +1.00 d_bb(b,f) d_bb(i,m) f_bb(a,e)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") += tmps_["30_bbbb_vvoo"]("f,b,m,i") * F["bb_vv"]("a,e");

        // H22_bbb_bbb += -1.00 d_bb(b,e) d_bb(i,m) f_bb(a,f)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") -= tmps_["30_bbbb_vvoo"]("e,b,m,i") * F["bb_vv"]("a,f");
        tmps_["30_bbbb_vvoo"].~TArrayD();

        // H22_aaa_aaa += -1.00 d_aa(a,f) d_aa(i,m) f_aa(b,e)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") -= tmps_["31_aaaa_vvoo"]("f,a,m,i") * F["aa_vv"]("b,e");

        // H22_aaa_aaa += +1.00 d_aa(b,f) d_aa(i,m) f_aa(a,e)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") += tmps_["31_aaaa_vvoo"]("f,b,m,i") * F["aa_vv"]("a,e");

        // H22_aaa_aaa += -1.00 d_aa(b,e) d_aa(i,m) f_aa(a,f)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") -= tmps_["31_aaaa_vvoo"]("e,b,m,i") * F["aa_vv"]("a,f");

        // H22_aaa_aaa += +1.00 d_aa(a,e) d_aa(i,m) f_aa(b,f)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") += tmps_["31_aaaa_vvoo"]("e,a,m,i") * F["aa_vv"]("b,f");

        // H22_aba_aba += +1.00 d_aa(a,e) d_aa(i,m) f_bb(b,f)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_aba("a,b,i,e,f,m") += tmps_["31_aaaa_vvoo"]("e,a,m,i") * F["bb_vv"]("b,f");
        tmps_["31_aaaa_vvoo"].~TArrayD();

        // tmps_[32_bb_vo](a,i) = 1.00 F[aa_ov](j,c) * t2[abab_vvoo](c,a,j,i) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
        tmps_["32_bb_vo"]("a,i")  = F["aa_ov"]("j,c") * t2["abab_vvoo"]("c,a,j,i");

        // H21_abb_a += +1.00 d_aa(a,e) f_aa(j,c) t2_abab(c,b,j,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_abb_a("a,b,i,e") += tmps_["32_bb_vo"]("b,i") * Id["aa_vv"]("a,e");

        // H21_bbb_b += +1.00 d_bb(a,e) f_aa(j,c) t2_abab(c,b,j,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_bbb_b("a,b,i,e") += tmps_["32_bb_vo"]("b,i") * Id["bb_vv"]("a,e");

        // H21_bbb_b += -1.00 d_bb(b,e) f_aa(j,c) t2_abab(c,a,j,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_bbb_b("a,b,i,e") -= tmps_["32_bb_vo"]("a,i") * Id["bb_vv"]("b,e");
        tmps_["32_bb_vo"].~TArrayD();

        // tmps_[33_bb_vo](b,i) = 1.00 F[bb_ov](j,c) * t2[bbbb_vvoo](c,b,i,j) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
        tmps_["33_bb_vo"]("b,i")  = F["bb_ov"]("j,c") * t2["bbbb_vvoo"]("c,b,i,j");

        // H21_abb_a += -1.00 d_aa(a,e) f_bb(j,c) t2_bbbb(c,b,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_abb_a("a,b,i,e") -= tmps_["33_bb_vo"]("b,i") * Id["aa_vv"]("a,e");

        // H21_bbb_b += -1.00 d_bb(a,e) f_bb(j,c) t2_bbbb(c,b,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_bbb_b("a,b,i,e") -= tmps_["33_bb_vo"]("b,i") * Id["bb_vv"]("a,e");

        // H21_bbb_b += +1.00 d_bb(b,e) f_bb(j,c) t2_bbbb(c,a,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_bbb_b("a,b,i,e") += tmps_["33_bb_vo"]("a,i") * Id["bb_vv"]("b,e");
        tmps_["33_bb_vo"].~TArrayD();

        // tmps_[34_aa_vo](a,i) = 1.00 F[bb_ov](j,c) * t2[abab_vvoo](a,c,i,j) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
        tmps_["34_aa_vo"]("a,i")  = F["bb_ov"]("j,c") * t2["abab_vvoo"]("a,c,i,j");

        // H21_aaa_a += -1.00 d_aa(b,e) f_bb(j,c) t2_abab(a,c,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aaa_a("a,b,i,e") -= tmps_["34_aa_vo"]("a,i") * Id["aa_vv"]("b,e");

        // H21_aaa_a += +1.00 d_aa(a,e) f_bb(j,c) t2_abab(b,c,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aaa_a("a,b,i,e") += tmps_["34_aa_vo"]("b,i") * Id["aa_vv"]("a,e");

        // H21_aba_b += -1.00 d_bb(b,e) f_bb(j,c) t2_abab(a,c,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aba_b("a,b,i,e") -= tmps_["34_aa_vo"]("a,i") * Id["bb_vv"]("b,e");
        tmps_["34_aa_vo"].~TArrayD();

        // tmps_[35_aa_vo](b,i) = 1.00 F[aa_ov](j,c) * t2[aaaa_vvoo](c,b,i,j) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
        tmps_["35_aa_vo"]("b,i")  = F["aa_ov"]("j,c") * t2["aaaa_vvoo"]("c,b,i,j");

        // H21_aaa_a += -1.00 d_aa(a,e) f_aa(j,c) t2_aaaa(c,b,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aaa_a("a,b,i,e") -= tmps_["35_aa_vo"]("b,i") * Id["aa_vv"]("a,e");

        // H21_aaa_a += +1.00 d_aa(b,e) f_aa(j,c) t2_aaaa(c,a,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aaa_a("a,b,i,e") += tmps_["35_aa_vo"]("a,i") * Id["aa_vv"]("b,e");

        // H21_aba_b += +1.00 d_bb(b,e) f_aa(j,c) t2_aaaa(c,a,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
        H21_aba_b("a,b,i,e") += tmps_["35_aa_vo"]("a,i") * Id["bb_vv"]("b,e");
        tmps_["35_aa_vo"].~TArrayD();

        // tmps_[36_bbbb_vvoo](f,b,i,m) = 1.00 Id[bb_vv](b,f) * F[bb_oo](m,i) // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
        tmps_["36_bbbb_vvoo"]("f,b,i,m")  = Id["bb_vv"]("b,f") * F["bb_oo"]("m,i");

        // H22_abb_abb += -1.00 d_aa(a,e) d_bb(b,f) f_bb(m,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_abb("a,b,i,e,f,m") -= Id["aa_vv"]("a,e") * tmps_["36_bbbb_vvoo"]("f,b,i,m");

        // H22_bbb_bbb += -1.00 d_bb(a,e) d_bb(b,f) f_bb(m,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") -= tmps_["36_bbbb_vvoo"]("e,a,i,m") * Id["bb_vv"]("b,f");

        // H22_bbb_bbb += +1.00 d_bb(b,e) d_bb(a,f) f_bb(m,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_bbb_bbb("a,b,i,e,f,m") += tmps_["36_bbbb_vvoo"]("e,b,i,m") * Id["bb_vv"]("a,f");
        tmps_["36_bbbb_vvoo"].~TArrayD();

        // H22_abb_abb += +1.00 d_aa(a,e) d_bb(i,m) f_bb(b,f)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_abb_abb("a,b,i,e,f,m") += tmps_["37_aabb_vvoo"]("e,a,m,i") * F["bb_vv"]("b,f");
        tmps_["37_aabb_vvoo"].~TArrayD();

        // H22_aba_aba += +1.00 d_bb(b,f) d_aa(i,m) f_aa(a,e)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_aba("a,b,i,e,f,m") += tmps_["38_bbaa_vvoo"]("f,b,m,i") * F["aa_vv"]("a,e");
        tmps_["38_bbaa_vvoo"].~TArrayD();

        // tmps_[39_aaaa_vvoo](e,a,i,m) = 1.00 Id[aa_vv](a,e) * F[aa_oo](m,i) // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
        tmps_["39_aaaa_vvoo"]("e,a,i,m")  = Id["aa_vv"]("a,e") * F["aa_oo"]("m,i");

        // H22_aaa_aaa += -1.00 d_aa(a,e) d_aa(b,f) f_aa(m,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") -= tmps_["39_aaaa_vvoo"]("e,a,i,m") * Id["aa_vv"]("b,f");

        // H22_aaa_aaa += +1.00 d_aa(b,e) d_aa(a,f) f_aa(m,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aaa_aaa("a,b,i,e,f,m") += tmps_["39_aaaa_vvoo"]("e,b,i,m") * Id["aa_vv"]("a,f");

        // H22_aba_aba += -1.00 d_aa(a,e) d_bb(b,f) f_aa(m,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
        H22_aba_aba("a,b,i,e,f,m") -= tmps_["39_aaaa_vvoo"]("e,a,i,m") * Id["bb_vv"]("b,f");
        tmps_["39_aaaa_vvoo"].~TArrayD();

    }
    world_.gop.fence();

    // build the full Hamiltonian

    // grab the dimensions
    size_t va = va_, vb = vb_, oa = oa_, ob = ob_;

    // get reduced dimensions
    size_t vaa = (va*va-va)/2, vbb = (vb*vb-vb)/2, oaa = (oa*oa-oa)/2, obb = (ob*ob-ob)/2;
    if (va == 1) vaa = 1;
    if (vb == 1) vbb = 1;
    if (oa == 1) oaa = 1;
    if (ob == 1) obb = 1;


    //  a, b, aaa, abb, aba, bbb
    size_t dim = va + vb + vaa*oa + va*vb*ob + va*vb*oa + vbb*ob;

    // initialize Hamiltonian
    auto pH = (double**) malloc(dim*sizeof(double*));
    pH[0] = (double*) malloc(dim*dim*sizeof(double));
    memset(pH[0], 0, dim*dim*sizeof(double));
    for (size_t i = 1; i < dim; i++)
        pH[i] = pH[i-1] + dim;

    // setup offsets; start with left indices
    size_t lid = 0;
    size_t rid = 0;

    // left a block
    {
        // right a block
        foreach_inplace(H11_a_a, [pH, lid, rid](auto &tile){
            for (auto &x : tile.range()) {
                size_t a = x[0], e = x[1];
                pH[a + lid][e + rid] = tile[x];
            }
        }); rid += va;

        rid += vb; // skip right b block

        // right aaa block
        foreach_inplace(H12_a_aaa, [pH, lid, rid, va, oa](auto &tile){
            for (auto &x : tile.range()) {
                size_t a = x[0], e = x[1], f = x[2], m = x[3];
                if (e < f) {
                    size_t ef = EOM_Driver::sqr_2_tri_idx(e, f, va);
                    size_t efm = ef * oa + m;
                    pH[a + lid][efm + rid] = tile[x];
                }
            }
        }); rid += vaa*oa;

        // right abb block
        foreach_inplace(H12_a_abb, [pH, lid, rid, vb, ob](auto &tile){
            for (auto &x : tile.range()) {
                size_t a = x[0], e = x[1], f = x[2], m = x[3];
                size_t efm = e * vb * ob + f * ob + m;
                pH[a + lid][efm + rid] = tile[x];
            }
        }); rid += va*vb*ob;

        // other right blocks are empty
    } lid += va;

    // left b block
    rid = 0;
    {
        rid += va; // skip right a block

        // right b block
        foreach_inplace(H11_b_b, [pH, lid, rid](auto &tile){
            for (auto &x : tile.range()) {
                size_t a = x[0], e = x[1];
                pH[a + lid][e + rid] = tile[x];
            }
        }); rid += vb;

        rid += vaa*oa;   // skip right aaa block
        rid += va*vb*ob; // skip right abb block

        // right aba block
        foreach_inplace(H12_b_aba, [pH, lid, rid, va, vb, oa](auto &tile){
            for (auto &x : tile.range()) {
                size_t a = x[0], e = x[1], f = x[2], m = x[3];
                size_t efm = e * vb * oa + f * oa + m;
                pH[a + lid][efm + rid] = tile[x];
            }
        }); rid += va*vb*oa;

        // right bbb block
        foreach_inplace(H12_b_bbb, [pH, lid, rid, vb, ob](auto &tile){
            for (auto &x : tile.range()) {
                size_t a = x[0], e = x[1], f = x[2], m = x[3];
                if (e < f) {
                    size_t ef = EOM_Driver::sqr_2_tri_idx(e, f, vb);
                    size_t efm = ef * ob + m;
                    pH[a + lid][efm + rid] = tile[x];
                }
            }
        }); rid += vbb*ob;

    } lid += vb;

    // left aaa block
    rid = 0;
    {
        // right a block
        foreach_inplace(H21_aaa_a, [pH, lid, rid, va, oa](auto &tile){
            for (auto &x : tile.range()) {
                size_t a = x[0], b = x[1], i = x[2], e = x[3];
                if (a < b) {
                    size_t ab = EOM_Driver::sqr_2_tri_idx(a, b, va);
                    size_t abi = ab * oa + i;
                    pH[abi + lid][e + rid] = tile[x];
                }
            }
         }); rid += va;

         rid += vb; // skip right b block

        // right aaa block
        foreach_inplace(H22_aaa_aaa, [pH, lid, rid, va, oa](auto &tile){
            for (auto &x : tile.range()) {
                size_t a = x[0], b = x[1], i = x[2], e = x[3], f = x[4], m = x[5];
                if (a < b && e < f) {
                    size_t ab = EOM_Driver::sqr_2_tri_idx(a, b, va);
                    size_t ef = EOM_Driver::sqr_2_tri_idx(e, f, va);
                    size_t abi = ab * oa + i;
                    size_t efm = ef * oa + m;
                    pH[abi + lid][efm + rid] = tile[x];
                }
            }
        }); rid += vaa*oa;

        // right abb block
        foreach_inplace(H22_aaa_abb, [pH, lid, rid, va, vb, oa, ob](auto &tile){
            for (auto &x : tile.range()) {
                size_t a = x[0], b = x[1], i = x[2], e = x[3], f = x[4], m = x[5];
                if (a < b) {
                    size_t ab = EOM_Driver::sqr_2_tri_idx(a, b, va);
                    size_t abi = ab * oa + i;
                    size_t efm = e * vb * ob + f * ob + m;
                    pH[abi + lid][efm + rid] = tile[x];
                }
            }
        }); rid += va*vb*ob;

        // other right blocks are empty
    } lid += vaa*oa;

    // left abb block
    rid = 0;
    {
        // right a block
        foreach_inplace(H21_abb_a, [pH, lid, rid, va, vb, ob](auto &tile){
            for (auto &x : tile.range()) {
                size_t a = x[0], b = x[1], i = x[2], e = x[3];
                size_t abi = a * vb * ob + b * ob + i;
                pH[abi + lid][e + rid] = tile[x];
            }
        }); rid += va;

        rid += vb; // skip right b block

        // right aaa block
        foreach_inplace(H22_abb_aaa, [pH, lid, rid, va, vb, oa, ob](auto &tile){
            for (auto &x : tile.range()) {
                size_t a = x[0], b = x[1], i = x[2], e = x[3], f = x[4], m = x[5];
                if (e < f) {
                    size_t abi = a * vb * ob + b * ob + i;
                    size_t ef = EOM_Driver::sqr_2_tri_idx(e, f, va);
                    size_t efm = ef * oa + m;
                    pH[abi + lid][efm + rid] = tile[x];
                }
            }
        }); rid += vaa*oa;

        // right abb block
        foreach_inplace(H22_abb_abb, [pH, lid, rid, va, vb, ob](auto &tile){
            for (auto &x : tile.range()) {
                size_t a = x[0], b = x[1], i = x[2], e = x[3], f = x[4], m = x[5];
                size_t abi = a * vb * ob + b * ob + i;
                size_t efm = e * vb * ob + f * ob + m;
                pH[abi + lid][efm + rid] = tile[x];
            }
        }); rid += va*vb*ob;

        // other right blocks are empty
    } lid += va*vb*ob;

    // left aba block
    rid = 0;
    {
        rid += va; // skip right a block

        // right b block
        foreach_inplace(H21_aba_b, [pH, lid, rid, va, vb, oa](auto &tile){
            for (auto &x : tile.range()) {
                size_t a = x[0], b = x[1], i = x[2], e = x[3];
                size_t abi = a * vb * oa + b * oa + i;
                pH[abi + lid][e + rid] = tile[x];
            }
        }); rid += vb;

        rid += vaa*oa;   // skip right aaa block
        rid += va*vb*ob; // skip right abb block

        // right aba block
        foreach_inplace(H22_aba_aba, [pH, lid, rid, va, vb, oa](auto &tile){
            for (auto &x : tile.range()) {
                size_t a = x[0], b = x[1], i = x[2], e = x[3], f = x[4], m = x[5];
                size_t abi = a * vb * oa + b * oa + i;
                size_t efm = e * vb * oa + f * oa + m;
                pH[abi + lid][efm + rid] = tile[x];
            }
        }); rid += va*vb*oa;

        // right bbb block
        foreach_inplace(H22_aba_bbb, [pH, lid, rid, va, vb, oa, ob](auto &tile){
            for (auto &x : tile.range()) {
                size_t a = x[0], b = x[1], i = x[2], e = x[3], f = x[4], m = x[5];
                if (e < f) {
                    size_t abi = a * vb * oa + b * oa + i;
                    size_t ef = EOM_Driver::sqr_2_tri_idx(e, f, vb);
                    size_t efm = ef * ob + m;
                    pH[abi + lid][efm + rid] = tile[x];
                }
            }
        }); rid += vbb*ob;
    } lid += va*vb*oa;

    // left bbb block
    rid = 0;
    {
        rid += va; // skip right a block

        // right b block
        foreach_inplace(H21_bbb_b, [pH, lid, rid, vb, ob](auto &tile){
            for (auto &x : tile.range()) {
                size_t a = x[0], b = x[1], i = x[2], e = x[3];
                if (a < b) {
                    size_t ab = EOM_Driver::sqr_2_tri_idx(a, b, vb);
                    size_t abi = ab * ob + i;
                    pH[abi + lid][e + rid] = tile[x];
                }
            }
        }); rid += vb;

        rid += vaa*oa;   // skip right aaa block
        rid += va*vb*ob; // skip right abb block

        // right aba block
        foreach_inplace(H22_bbb_aba, [pH, lid, rid, va, vb, oa, ob](auto &tile){
            for (auto &x : tile.range()) {
                size_t a = x[0], b = x[1], i = x[2], e = x[3], f = x[4], m = x[5];
                if (a < b) {
                    size_t ab = EOM_Driver::sqr_2_tri_idx(a, b, vb);
                    size_t abi = ab * ob + i;
                    size_t efm = e * vb * oa + f * oa + m;
                    pH[abi + lid][efm + rid] = tile[x];
                }
            }
        }); rid += va*vb*oa;

        // right bbb block
        foreach_inplace(H22_bbb_bbb, [pH, lid, rid, vb, ob](auto &tile){
            for (auto &x : tile.range()) {
                size_t a = x[0], b = x[1], i = x[2], e = x[3], f = x[4], m = x[5];
                if (a < b && e < f) {
                    size_t ab = EOM_Driver::sqr_2_tri_idx(a, b, vb);
                    size_t abi = ab * ob + i;
                    size_t ef = EOM_Driver::sqr_2_tri_idx(e, f, vb);
                    size_t efm = ef * ob + m;
                    pH[abi + lid][efm + rid] = tile[x];
                }
            }
        }); rid += vbb*ob;
    } lid += vbb * ob;

    // DONE Building!
    world_.gop.fence();

    Printf("dim = %d; lid = %d; rid = %d\n", dim, lid, rid);

    // add energy to diagonal
    for (size_t i = 0; i < dim; i++) {
        pH[i][i] += cc_wfn_->cc_energy_;
    }

    // initialize containers for eigenvalues and eigenvectors
    double *eval_r = (double*) malloc(dim*sizeof(double));
    double *eval_i = (double*) malloc(dim*sizeof(double));
    double *levecs = (double*) malloc(dim*dim*sizeof(double));
    double *revecs = (double*) malloc(dim*dim*sizeof(double));

    long int info;
    char vl = 'V';
    char vr = 'V';
    long int n = dim;
    long int lwork = 4 * n;
    auto *work = (double *) malloc(lwork * sizeof(double));

    // diagonalize the Hamiltonian
    psi::fnocc::DGEEV(vl, vr, n, *pH, n, eval_r, eval_i, levecs, n, revecs, n, work, lwork, info);

    // perform argument sort by pairing eigenvalues with index
    auto *eig_pair = (std::pair<std::complex<double>, size_t> *) malloc(
            n * sizeof(std::pair<std::complex<double>, size_t>));

    for (size_t i = 0; i < n; i++) {
        eig_pair[i].first = std::complex<double>(eval_r[i], eval_i[i]);
        eig_pair[i].second = i;
    }

    // sort by eigenvalue
    bool ascending = true;
    std::sort(eig_pair, eig_pair + n,
              [ascending](std::pair<std::complex<double>, size_t> &l, std::pair<std::complex<double>, size_t> &r) {
                  if (ascending) {
                      bool moveOrder = l.first.real() < r.first.real();
                      if (!moveOrder) {
                          if (fabs(l.first.real() - r.first.real()) <= 1e-16) {
                              return l.first.imag() > r.first.imag();
                          }
                      }
                      return moveOrder;
                  } else {
                      bool moveOrder = l.first.real() > r.first.real();
                      if (!moveOrder) {
                          if (fabs(l.first.real() - r.first.real()) <= 1e-16) {
                              return l.first.imag() < r.first.imag();
                          }
                      }
                      return moveOrder;
                  }
              }
    );

    double** revec_p = revec_->pointer();
    double** levec_p = levec_->pointer();
    double* eigvals_p = eigvals_->pointer();


    // apply sorted indices to reorder eigensystem
    for (size_t i = 0; i < dim; ++i) {

        // get ordered index
        size_t min_i = eig_pair[i].second;

        // reorder eigenvalue
        eigvals_p[i] = eval_r[min_i];

        // reorder eigenvector
        for (size_t j = 0; j < n; j++) {
            revec_p[i][j] = revecs[min_i*dim + j];
            levec_p[i][j] = levecs[min_i*dim + j];
        }

        if (fabs(eval_i[min_i]) >= 1e-16) {

            // get conjugate pair index
            size_t min_i1 = eig_pair[i+1].second;
            eigvals_p[i+1] = eval_r[min_i1];

            // reorder eigenvector
            for (size_t j = 0; j < dim; j++) {
                revec_p[i+1][j] = revecs[min_i1*dim + j];
                levec_p[i+1][j] = levecs[min_i1*dim + j];
            }

            // skip next index
            ++i;
        }
    }

    // clean up
    free(pH[0]); free(pH);
    free(eval_r); free(eval_i);
    free(levecs); free(revecs);
    free(work);
}

void hilbert::EOM_EA_CCSD::build_common_ops() {

    if (!cc_wfn_->has_t1_integrals_)
        cc_wfn_->transform_integrals(true);

    // get reference to electronic integrals
    TArrayMap &V_blks_ = cc_wfn_->V_blks_;
    TArrayMap &F_blks_ = cc_wfn_->F_blks_;
    TArrayMap &Id_blks_ = cc_wfn_->Id_blks_;
    TArrayMap &eri = cc_wfn_->V_blks_;
    TArrayMap &f = cc_wfn_->F_blks_;
    TArrayMap &Id = cc_wfn_->Id_blks_;

    TArrayMap t2 = {
        {"aaaa", cc_wfn_->amplitudes_["t2_aaaa"]},
        {"abab", cc_wfn_->amplitudes_["t2_abab"]},
        {"bbbb", cc_wfn_->amplitudes_["t2_bbbb"]}
    };

    {
        scalars_["1"]  = -2.00 * dot(Id["aa_oo"]("i,j"), f["aa_oo"]("i,j"));
        scalars_["1"] -= 2.00 * dot(Id["bb_oo"]("l,m"), f["bb_oo"]("l,m"));
        scalars_["1"] -= 2.00 * dot(eri["abab_oovv"]("k,l,a,c"), t2["abab"]("a,c,k,l"));
        scalars_["1"] += 0.50 * dot(eri["bbbb_oovv"]("l,n,d,c"), t2["bbbb"]("d,c,n,l"));
        scalars_["1"] += 1.00 * dot(Id["aa_oo"]("i,j") * eri["aaaa_oooo"]("i,k,j,I"), Id["aa_oo"]("k,I"));
        scalars_["1"] += 1.00 * dot(Id["bb_oo"]("l,m") * eri["bbbb_oooo"]("l,n,m,o"), Id["bb_oo"]("n,o"));
        scalars_["1"] += 2.00 * dot(Id["bb_oo"]("l,m") * eri["abab_oooo"]("k,l,I,m"), Id["aa_oo"]("k,I"));
        scalars_["1"] += 0.50 * dot(eri["aaaa_oovv"]("i,k,a,b"), t2["aaaa"]("a,b,k,i"));
    }

    {
        reused_["1_bb_vv"]("a,b")  = eri["abab_oovv"]("i,j,c,a") * t2["abab"]("c,b,i,j");
        reused_["2_bb_vv"]("a,b")  = t2["bbbb"]("c,a,i,j") * eri["bbbb_oovv"]("j,i,c,b");
        reused_["3_aa_vv"]("a,b")  = t2["aaaa"]("c,a,i,j") * eri["aaaa_oovv"]("j,i,c,b");
        reused_["4_aa_vv"]("a,b")  = t2["abab"]("a,c,i,j") * eri["abab_oovv"]("i,j,b,c");
        reused_["5_aaaa_vovv"]("a,i,b,c")  = t2["aaaa"]("d,a,i,j") * eri["aaaa_vovv"]("b,j,c,d");
        reused_["6_aaaa_vvvo"]("a,b,c,i")  = eri["abab_vovv"]("a,j,b,d") * t2["abab"]("c,d,i,j");
        reused_["7_aabb_vvvo"]("a,b,c,i")  = eri["aaaa_vovv"]("a,j,b,d") * t2["abab"]("d,c,j,i");
        reused_["8_aabb_vvvo"]("a,b,c,i")  = eri["abab_vovv"]("a,j,b,d") * t2["bbbb"]("d,c,i,j");
        reused_["9_baab_vvvo"]("a,b,c,i")  = eri["baab_vovv"]("a,j,b,d") * t2["abab"]("c,d,j,i");
        reused_["10_aaaa_vovv"]("a,i,b,c")  = eri["aaaa_oovo"]("j,k,a,i") * t2["aaaa"]("b,c,k,j");
        reused_["11_abab_vovv"]("a,i,b,c")  = eri["abab_oovo"]("j,k,a,i") * t2["abab"]("b,c,j,k");
        reused_["12_aa_vv"]("a,b")  = eri["aaaa_oovv"]("i,j,a,c") * t2["aaaa"]("c,b,j,i");
        reused_["13_aa_ov"]("i,a")  = -2.00 * f["aa_ov"]("l,c") * t2["aaaa"]("c,a,i,l");
        reused_["13_aa_ov"]("i,a") += 2.00 * eri["abab_vovv"]("a,j,c,e") * t2["abab"]("c,e,i,j");
        reused_["13_aa_ov"]("i,a") += eri["aaaa_oovo"]("l,k,c,i") * t2["aaaa"]("c,a,k,l");
        reused_["13_aa_ov"]("i,a") += eri["aaaa_vovv"]("a,l,c,d") * t2["aaaa"]("c,d,i,l");
        reused_["13_aa_ov"]("i,a") += 2.00 * eri["abba_oovo"]("k,j,b,i") * t2["abab"]("a,b,k,j");
        reused_["13_aa_ov"]("i,a") += 2.00 * f["bb_ov"]("j,b") * t2["abab"]("a,b,i,j");
        reused_["14_bb_vo"]("a,i")  = -1.00 * eri["abab_oovo"]("j,k,b,i") * t2["abab"]("b,a,j,k");
        reused_["14_bb_vo"]("a,i") -= eri["baab_vovv"]("a,m,b,d") * t2["abab"]("b,d,m,i");
        reused_["14_bb_vo"]("a,i") += 0.50 * eri["bbbb_vovv"]("a,k,c,d") * t2["bbbb"]("c,d,i,k");
        reused_["14_bb_vo"]("a,i") += f["aa_ov"]("m,b") * t2["abab"]("b,a,m,i");
        reused_["14_bb_vo"]("a,i") -= f["bb_ov"]("k,c") * t2["bbbb"]("c,a,i,k");
        reused_["14_bb_vo"]("a,i") += 0.50 * eri["bbbb_oovo"]("k,l,c,i") * t2["bbbb"]("c,a,l,k");
        reused_["15_aabb_vvvo"]("a,b,c,i")  = eri["aaaa_vovv"]("a,j,d,b") * t2["abab"]("d,c,j,i");
        reused_["16_bbaa_voov"]("a,i,j,b")  = t2["abab"]("c,a,k,i") * eri["aaaa_oovv"]("k,j,c,b");
        reused_["17_bbbb_voov"]("a,i,j,b")  = t2["bbbb"]("c,a,i,k") * eri["bbbb_oovv"]("k,j,c,b");
        reused_["18_bb_oo"]("i,j")  = t2["bbbb"]("a,b,i,k") * eri["bbbb_oovv"]("k,j,a,b");
        reused_["19_bbbb_voov"]("a,i,j,b")  = t2["abab"]("c,a,k,i") * eri["abab_oovv"]("k,j,c,b");
        reused_["20_bbaa_voov"]("a,i,j,b")  = t2["bbbb"]("c,a,i,k") * eri["abab_oovv"]("j,k,b,c");
        reused_["21_bb_oo"]("i,j")  = t2["abab"]("a,b,k,i") * eri["abab_oovv"]("k,j,a,b");
        reused_["22_bb_vv"]("a,b")  = eri["bbbb_oovv"]("i,j,a,c") * t2["bbbb"]("c,b,j,i");
        reused_["23_abba_vvvo"]("a,b,c,i")  = eri["abab_vovv"]("a,j,d,b") * t2["abab"]("d,c,i,j");
        reused_["23_abba_vvvo"]("a,b,c,i") += eri["abba_oovo"]("k,j,b,i") * t2["abab"]("a,c,k,j");
        reused_["24_aabb_vovv"]("a,i,b,c")  = t2["aaaa"]("d,a,i,j") * eri["baab_vovv"]("b,j,d,c");
        reused_["24_aabb_vovv"]("a,i,b,c") += t2["abab"]("a,e,i,k") * eri["bbbb_vovv"]("b,k,c,e");
        reused_["25_bbbb_vovv"]("a,i,b,c")  = t2["abab"]("e,a,k,i") * eri["baab_vovv"]("b,k,e,c");
        reused_["25_bbbb_vovv"]("a,i,b,c") += 0.25 * eri["bbbb_oovo"]("j,l,c,i") * t2["bbbb"]("b,a,l,j");
        reused_["25_bbbb_vovv"]("a,i,b,c") += t2["bbbb"]("d,a,i,j") * eri["bbbb_vovv"]("b,j,c,d");
        reused_["26_aaaa_vvvo"]("a,b,c,i")  = eri["aaaa_vovv"]("a,j,d,b") * t2["aaaa"]("d,c,i,j");
        reused_["27_aaaa_voov"]("a,i,j,b")  = t2["aaaa"]("c,a,i,k") * eri["aaaa_oovv"]("k,j,c,b");
        reused_["28_aabb_voov"]("a,i,j,b")  = t2["abab"]("a,c,i,k") * eri["bbbb_oovv"]("k,j,c,b");
        reused_["29_aabb_voov"]("a,i,j,b")  = t2["aaaa"]("c,a,i,k") * eri["abab_oovv"]("k,j,c,b");
        reused_["30_aa_oo"]("i,j")  = eri["aaaa_oovv"]("i,k,a,b") * t2["aaaa"]("a,b,j,k");
        reused_["31_aa_oo"]("i,j")  = t2["abab"]("a,b,i,k") * eri["abab_oovv"]("j,k,a,b");
        reused_["32_aa_oo"]("i,j")  = t2["aaaa"]("a,b,i,k") * eri["aaaa_oovv"]("k,j,a,b");
        reused_["33_aabb_voov"]("a,i,j,b")  = t2["abab"]("a,c,i,k") * eri["bbbb_oovv"]("j,k,b,c");
        reused_["34_bbbb_voov"]("a,i,j,b")  = t2["bbbb"]("c,a,i,k") * eri["bbbb_oovv"]("j,k,b,c");
        reused_["35_bb_oo"]("i,j")  = eri["bbbb_oovv"]("i,k,a,b") * t2["bbbb"]("a,b,j,k");
        reused_["36_aaaa_voov"]("a,i,j,b")  = t2["aaaa"]("c,a,i,k") * eri["aaaa_oovv"]("j,k,b,c");
        reused_["37_bbaa_voov"]("a,i,j,b")  = t2["abab"]("c,a,k,i") * eri["aaaa_oovv"]("j,k,b,c");
    }

    world_.gop.fence();
}

