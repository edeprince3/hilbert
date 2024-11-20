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
 *  You should have renergyeived a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 *  @END LICENSE
 */

#include "cc_cavity/include/qed_ccsd_21/qed_ccsd_21.h"


using namespace std;
using namespace TA;
using namespace hilbert;

void QED_CCSD_21::resid_21_2() {

    // unpack integrals
    TArrayMap &Id         = Id_blks_;
    TArrayMap &f          = F_blks_;
    TArrayMap &eri        = V_blks_;

    // get effective dipole integrals
    TArrayMap dp = effective_dipole();

    // extract scalars
    double &energy  = scalars_["energy"];
    double &cenergy = scalars_["cenergy"];
    double  w0      = scalars_["w0"];

    // extract residuals
    TArrayD &rt1_aa     = residuals_["t1_aa"];
    TArrayD &rt1_bb     = residuals_["t1_bb"];
    TArrayD &rt1_1_aa   = residuals_["t1_1_aa"];
    TArrayD &rt1_1_bb   = residuals_["t1_1_bb"];
    TArrayD &rt2_aaaa   = residuals_["t2_aaaa"];
    TArrayD &rt2_abab   = residuals_["t2_abab"];
    TArrayD &rt2_bbbb   = residuals_["t2_bbbb"];
    TArrayD &rt2_1_aaaa = residuals_["t2_1_aaaa"];
    TArrayD &rt2_1_abab = residuals_["t2_1_abab"];
    TArrayD &rt2_1_bbbb = residuals_["t2_1_bbbb"];
    double  &rt0_1 = scalars_["rt0_1"];

    // extract coherents
    TArrayD &crt1_aa     = tmps_["ct1_aa"];
    TArrayD &crt1_bb     = tmps_["ct1_bb"];
    TArrayD &crt1_1_aa   = tmps_["ct1_1_aa"];
    TArrayD &crt1_1_bb   = tmps_["ct1_1_bb"];
    TArrayD &crt2_aaaa   = tmps_["ct2_aaaa"];
    TArrayD &crt2_abab   = tmps_["ct2_abab"];
    TArrayD &crt2_bbbb   = tmps_["ct2_bbbb"];
    TArrayD &crt2_1_aaaa = tmps_["ct2_1_aaaa"];
    TArrayD &crt2_1_abab = tmps_["ct2_1_abab"];
    TArrayD &crt2_1_bbbb = tmps_["ct2_1_bbbb"];
    double  &crt0_1 = scalars_["crt0_1"];

    // extract amplitudes
    double t0_1 = scalars_["t0_1"];

    // extract 1-body amplitudes
    std::map<std::string, TA::TArrayD> t1 {
            {"aa", amplitudes_["t1_aa"]},
            {"bb", amplitudes_["t1_bb"]}
    };
    std::map<std::string, TA::TArrayD> t1_1 {
            {"aa", amplitudes_["t1_1_aa"]},
            {"bb", amplitudes_["t1_1_bb"]}
    };

    // extract 2-body amplitudes
    std::map<std::string, TA::TArrayD> t2 {
            {"aaaa", amplitudes_["t2_aaaa"]},
            {"abab", amplitudes_["t2_abab"]},
            {"bbbb", amplitudes_["t2_bbbb"]}
    };
    std::map<std::string, TA::TArrayD> t2_1 {
            {"aaaa", amplitudes_["t2_1_aaaa"]},
            {"abab", amplitudes_["t2_1_abab"]},
            {"bbbb", amplitudes_["t2_1_bbbb"]}
    };


    // flops: o1v1  = o1v1 o3v2 o2v3 o2v3 o2v2 o1v1 o1v1 o1v1 o3v2 o1v1 o2v2 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1
    //  mems: o1v1  = o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1
    tmps_["31_aa_vo"]("a,i")  = -1.00 * t0_1 * tmps_["25_aa_vo"]("a,i");
    tmps_["31_aa_vo"]("a,i") -= 0.50 * eri["aaaa_oovo"]("j,k,b,i") * t2["aaaa"]("b,a,k,j");
    tmps_["31_aa_vo"]("a,i") -= 0.50 * eri["aaaa_vovv"]("a,j,b,c") * t2["aaaa"]("b,c,i,j");
    tmps_["31_aa_vo"]("a,i") -= eri["abab_vovv"]("a,l,b,e") * t2["abab"]("b,e,i,l");
    tmps_["31_aa_vo"]("a,i") += f["aa_ov"]("j,b") * t2["aaaa"]("b,a,i,j");
    tmps_["31_aa_vo"]("a,i") += scalars_["3"] * t1_1["aa"]("a,i");
    tmps_["31_aa_vo"]("a,i") -= eri["abba_oovo"]("k,l,d,i") * t2["abab"]("a,d,k,l");
    tmps_["31_aa_vo"]("a,i") -= f["bb_ov"]("l,d") * t2["abab"]("a,d,i,l");
    tmps_["31_aa_vo"]("a,i") -= f["aa_vo"]("a,i");
    tmps_["31_aa_vo"]("a,i") += t0_1 * dp["aa_vo"]("a,i");

    // rt1_aa  = -1.00 f_aa(j,b) t2_aaaa(b,a,i,j)
    //        += -1.00 d-_bb(j,j) t1_1_aa(a,i)
    //        += -1.00 d-_aa(j,j) t1_1_aa(a,i)
    //        += +0.50 <a,j||b,c>_abab t2_abab(b,c,i,j)
    //        += +0.50 <a,j||c,b>_abab t2_abab(c,b,i,j)
    //        += -0.50 <k,j||i,b>_abab t2_abab(a,b,k,j)
    //        += -0.50 <j,k||i,b>_abab t2_abab(a,b,j,k)
    //        += +1.00 f_bb(j,b) t2_abab(a,b,i,j)
    //        += +1.00 f_aa(a,i)
    //        += -1.00 d-_aa(a,i) t0_1
    //        += -0.50 <j,a||b,c>_aaaa t2_aaaa(b,c,i,j)
    //        += -0.50 <k,j||b,i>_aaaa t2_aaaa(b,a,k,j)
    //        += +1.00 d-_aa(j,b) t2_aaaa(b,a,i,j) t0_1
    //        += -1.00 d-_bb(j,b) t2_abab(a,b,i,j) t0_1
    rt1_aa("a,i")  = -1.00 * tmps_["31_aa_vo"]("a,i");
    tmps_["31_aa_vo"].~TArrayD();

    // flops: o1v1  = o1v1
    //  mems: o1v1  = o1v1
    tmps_["32_aa_vo"]("a,i")  = -1.00 * t0_1 * tmps_["25_aa_vo"]("a,i");
    tmps_["32_aa_vo"].~TArrayD();

    // flops: o2v2  = o3v3 o2v2 o3v3 o2v2 o2v2 o2v2 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["33_bbbb_vovo"]("a,i,b,j")  = -1.00 * eri["bbbb_vovo"]("a,l,d,i") * t2["bbbb"]("d,b,j,l");
    tmps_["33_bbbb_vovo"]("a,i,b,j") += dp["bb_vo"]("a,i") * t1_1["bb"]("b,j");
    tmps_["33_bbbb_vovo"]("a,i,b,j") += eri["baab_vovo"]("a,k,c,i") * t2["abab"]("c,b,k,j");
    tmps_["33_bbbb_vovo"]("a,i,b,j") += t1_1["bb"]("b,j") * tmps_["29_bb_vo"]("a,i");

    // rt2_bbbb += +1.00 P(i,j) P(a,b) d-_bb(a,j) t1_1_bb(b,i)
    //          += -1.00 P(i,j) P(a,b) <k,a||c,j>_abab t2_abab(c,b,k,i)
    //          += +1.00 P(i,j) P(a,b) <k,a||c,j>_bbbb t2_bbbb(c,b,i,k)
    //          += +1.00 P(i,j) P(a,b) d-_aa(k,c) t2_abab(c,a,k,j) t1_1_bb(b,i)
    //          += -1.00 P(i,j) P(a,b) d-_bb(k,c) t2_bbbb(c,a,j,k) t1_1_bb(b,i)
    rt2_bbbb("a,b,i,j") += tmps_["33_bbbb_vovo"]("a,j,b,i");
    rt2_bbbb("a,b,i,j") -= tmps_["33_bbbb_vovo"]("a,i,b,j");
    rt2_bbbb("a,b,i,j") -= tmps_["33_bbbb_vovo"]("b,j,a,i");
    rt2_bbbb("a,b,i,j") += tmps_["33_bbbb_vovo"]("b,i,a,j");
    tmps_["33_bbbb_vovo"].~TArrayD();

    // flops: o0v2  = o2v3
    //  mems: o0v2  = o0v2
    tmps_["38_aa_vv"]("a,b")  = t2["aaaa"]("c,a,i,j") * eri["aaaa_oovv"]("j,i,c,b");

    // flops: o4v0  = o4v2
    //  mems: o4v0  = o4v0
    tmps_["37_aaaa_oooo"]("i,j,k,l")  = t2["aaaa"]("a,b,i,j") * eri["aaaa_oovv"]("k,l,a,b");

    // flops: o2v2  = o2v2
    //  mems: o2v2  = o2v2
    tmps_["36_aaaa_vvoo"]("a,b,i,j")  = 2.00 * scalars_["3"] * t2_1["aaaa"]("a,b,i,j");

    // flops: o1v1  = o1v1
    //  mems: o1v1  = o1v1
    tmps_["34_bb_vo"]("a,i")  = -1.00 * t0_1 * tmps_["29_bb_vo"]("a,i");
    tmps_["34_bb_vo"].~TArrayD();

    // flops: o0v2  = o2v3
    //  mems: o0v2  = o0v2
    tmps_["35_aa_vv"]("a,b")  = t2["aaaa"]("c,a,i,j") * eri["aaaa_oovv"]("j,i,b,c");

    // flops: o2v2  = o2v4 o4v2 o2v2 o2v2 o2v2 o2v3 o2v2 o4v2 o2v2 o2v3 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["39_aaaa_oovv"]("i,j,a,b")  = -1.00 * eri["aaaa_vvvv"]("a,b,c,d") * t2["aaaa"]("c,d,i,j");
    tmps_["39_aaaa_oovv"]("i,j,a,b") += eri["aaaa_oooo"]("l,k,i,j") * t2["aaaa"]("a,b,k,l");
    tmps_["39_aaaa_oovv"]("i,j,a,b") -= 2.00 * eri["aaaa_vvoo"]("a,b,i,j");
    tmps_["39_aaaa_oovv"]("i,j,a,b") += tmps_["36_aaaa_vvoo"]("a,b,i,j");
    tmps_["39_aaaa_oovv"]("i,j,a,b") -= t2["aaaa"]("d,b,i,j") * tmps_["38_aa_vv"]("a,d");
    tmps_["39_aaaa_oovv"]("i,j,a,b") += 0.50 * t2["aaaa"]("a,b,k,l") * tmps_["37_aaaa_oooo"]("i,j,l,k");
    tmps_["39_aaaa_oovv"]("i,j,a,b") -= t2["aaaa"]("c,a,i,j") * tmps_["35_aa_vv"]("b,c");
    tmps_["36_aaaa_vvoo"].~TArrayD();
    tmps_["35_aa_vv"].~TArrayD();

    // rt2_aaaa += +0.50 <l,k||i,j>_aaaa t2_aaaa(a,b,l,k)
    //          += +0.50 <a,b||c,d>_aaaa t2_aaaa(c,d,i,j)
    //          += +1.00 <a,b||i,j>_aaaa
    //          += -1.00 d-_bb(k,k) t2_1_aaaa(a,b,i,j)
    //          += -1.00 d-_aa(k,k) t2_1_aaaa(a,b,i,j)
    //          += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,a,l,k) t2_aaaa(d,b,i,j)
    //          += +0.25 <l,k||c,d>_aaaa t2_aaaa(a,b,l,k) t2_aaaa(c,d,i,j)
    //          += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,a,i,j) t2_aaaa(d,b,l,k)
    rt2_aaaa("a,b,i,j") -= 0.50 * tmps_["39_aaaa_oovv"]("i,j,a,b");
    tmps_["39_aaaa_oovv"].~TArrayD();

    // flops: o1v1  = o2v2 o2v2 o1v1
    //  mems: o1v1  = o1v1 o1v1 o1v1
    tmps_["46_bb_vo"]("a,i")  = -1.00 * dp["aa_ov"]("k,c") * t2_1["abab"]("c,a,k,i");
    tmps_["46_bb_vo"]("a,i") += dp["bb_ov"]("j,b") * t2_1["bbbb"]("b,a,i,j");

    // rt1_bb += +1.00 d-_bb(j,b) t2_1_bbbb(b,a,i,j)
    //        += -1.00 d-_aa(j,b) t2_1_abab(b,a,j,i)
    rt1_bb("a,i") += tmps_["46_bb_vo"]("a,i");

    // flops: o1v1  = o1v2 o2v1 o1v1
    //  mems: o1v1  = o1v1 o1v1 o1v1
    tmps_["47_bb_ov"]("i,a")  = -1.00 * dp["bb_vv"]("a,b") * t1_1["bb"]("b,i");
    tmps_["47_bb_ov"]("i,a") += dp["bb_oo"]("j,i") * t1_1["bb"]("a,j");

    // rt1_bb += +1.00 d-_bb(j,i) t1_1_bb(a,j)
    //        += -1.00 d-_bb(a,b) t1_1_bb(b,i)
    rt1_bb("a,i") += tmps_["47_bb_ov"]("i,a");

    // flops: o2v0  = o2v1
    //  mems: o2v0  = o2v0
    tmps_["45_bb_oo"]("i,j")  = t1_1["bb"]("a,i") * dp["bb_ov"]("j,a");

    // flops: o1v1  = o2v2
    //  mems: o1v1  = o1v1
    tmps_["44_bb_ov"]("i,a")  = eri["abab_oovv"]("j,i,b,a") * t1_1["aa"]("b,j");

    // flops: o1v1  = o2v2
    //  mems: o1v1  = o1v1
    tmps_["40_bb_vo"]("a,i")  = -1.00 * dp["aa_ov"]("j,b") * t2_1["abab"]("b,a,j,i");
    tmps_["40_bb_vo"].~TArrayD();

    // flops: o1v1  = o1v2
    //  mems: o1v1  = o1v1
    tmps_["41_bb_vo"]("a,i")  = -1.00 * dp["bb_vv"]("a,b") * t1_1["bb"]("b,i");
    tmps_["41_bb_vo"].~TArrayD();

    // flops: o1v1  = o2v2
    //  mems: o1v1  = o1v1
    tmps_["42_bb_vo"]("a,i")  = 2.00 * t1_1["bb"]("b,j") * tmps_["4_bbbb_voov"]("a,i,j,b");
    tmps_["42_bb_vo"].~TArrayD();

    // flops: o0v2  = o2v3
    //  mems: o0v2  = o0v2
    tmps_["43_bb_vv"]("a,b")  = t2["bbbb"]("c,a,i,j") * eri["bbbb_oovv"]("j,i,c,b");

    // flops: o1v1  = o2v3 o3v2 o1v1 o3v2 o1v1 o2v2 o1v1 o2v2 o1v1 o1v1 o2v1 o1v1 o1v1 o1v1 o2v2 o1v1 o1v1 o1v1 o1v2 o1v1 o2v2 o1v1 o2v3 o1v1 o1v1 o2v1 o1v1 o2v2 o1v1 o1v2 o1v1 o1v1 o1v2 o1v1 o2v2 o1v1 o2v1 o1v1 o1v1 o1v1 o2v2 o1v1 o2v2 o1v1
    //  mems: o1v1  = o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1
    tmps_["48_bb_vo"]("a,i")  = -2.00 * eri["baab_vovv"]("a,l,d,b") * t2_1["abab"]("d,b,l,i");
    tmps_["48_bb_vo"]("a,i") += eri["bbbb_oovo"]("k,j,c,i") * t2_1["bbbb"]("c,a,j,k");
    tmps_["48_bb_vo"]("a,i") -= 2.00 * eri["abab_oovo"]("m,k,d,i") * t2_1["abab"]("d,a,m,k");
    tmps_["48_bb_vo"]("a,i") -= 2.00 * eri["bbbb_vovo"]("a,k,c,i") * t1_1["bb"]("c,k");
    tmps_["48_bb_vo"]("a,i") -= 2.00 * f["bb_ov"]("k,c") * t2_1["bbbb"]("c,a,i,k");
    tmps_["48_bb_vo"]("a,i") -= 2.00 * dp["bb_vo"]("a,i");
    tmps_["48_bb_vo"]("a,i") -= 2.00 * f["bb_oo"]("k,i") * t1_1["bb"]("a,k");
    tmps_["48_bb_vo"]("a,i") += 2.00 * w0 * t1_1["bb"]("a,i");
    tmps_["48_bb_vo"]("a,i") -= 2.00 * eri["baab_vovo"]("a,l,d,i") * t1_1["aa"]("d,l");
    tmps_["48_bb_vo"]("a,i") -= 2.00 * scalars_["4"] * t1_1["bb"]("a,i");
    tmps_["48_bb_vo"]("a,i") += 2.00 * f["bb_vv"]("a,c") * t1_1["bb"]("c,i");
    tmps_["48_bb_vo"]("a,i") += 2.00 * f["aa_ov"]("l,d") * t2_1["abab"]("d,a,l,i");
    tmps_["48_bb_vo"]("a,i") += eri["bbbb_vovv"]("a,k,c,b") * t2_1["bbbb"]("c,b,i,k");
    tmps_["48_bb_vo"]("a,i") += 2.00 * t0_1 * tmps_["46_bb_vo"]("a,i");
    tmps_["48_bb_vo"]("a,i") += t1_1["bb"]("a,j") * tmps_["16_bb_oo"]("i,j");
    tmps_["48_bb_vo"]("a,i") += 2.00 * t1_1["aa"]("e,m") * tmps_["14_bbaa_voov"]("a,i,m,e");
    tmps_["48_bb_vo"]("a,i") += t1_1["bb"]("b,i") * tmps_["43_bb_vv"]("a,b");
    tmps_["48_bb_vo"]("a,i") -= 2.00 * t1_1["bb"]("b,i") * tmps_["6_bb_vv"]("a,b");
    tmps_["48_bb_vo"]("a,i") -= 2.00 * t1_1["bb"]("b,j") * tmps_["15_bbbb_voov"]("a,i,j,b");
    tmps_["48_bb_vo"]("a,i") += 4.00 * t1_1["bb"]("a,k") * tmps_["45_bb_oo"]("i,k");
    tmps_["48_bb_vo"]("a,i") += 2.00 * t0_1 * tmps_["47_bb_ov"]("i,a");
    tmps_["48_bb_vo"]("a,i") -= 2.00 * t2["bbbb"]("c,a,i,k") * tmps_["44_bb_ov"]("k,c");
    tmps_["48_bb_vo"]("a,i") += 2.00 * t1_1["bb"]("b,j") * tmps_["4_bbbb_voov"]("a,i,j,b");

    // rt1_1_bb += -0.50 <k,j||b,c>_bbbb t2_bbbb(b,c,i,j) t1_1_bb(a,k)
    //          += -0.50 <j,k||b,c>_abab t2_abab(b,c,j,i) t1_1_bb(a,k)
    //          += -0.50 <j,k||c,b>_abab t2_abab(c,b,j,i) t1_1_bb(a,k)
    //          += +1.00 d-_bb(j,b) t0_1 t2_1_bbbb(b,a,i,j)
    //          += -1.00 d-_aa(j,b) t0_1 t2_1_abab(b,a,j,i)
    //          += -1.00 <k,j||b,c>_aaaa t2_abab(b,a,j,i) t1_1_aa(c,k)
    //          += -0.50 <k,j||b,c>_bbbb t2_bbbb(b,a,k,j) t1_1_bb(c,i)
    //          += +0.50 <j,a||b,c>_abab t2_1_abab(b,c,j,i)
    //          += +0.50 <j,a||c,b>_abab t2_1_abab(c,b,j,i)
    //          += -0.50 <k,j||b,i>_bbbb t2_1_bbbb(b,a,k,j)
    //          += -0.50 <k,j||b,i>_abab t2_1_abab(b,a,k,j)
    //          += -0.50 <j,k||b,i>_abab t2_1_abab(b,a,j,k)
    //          += +1.00 <j,a||b,i>_bbbb t1_1_bb(b,j)
    //          += -1.00 f_bb(j,b) t2_1_bbbb(b,a,i,j)
    //          += -1.00 d+_bb(a,i)
    //          += -1.00 f_bb(j,i) t1_1_bb(a,j)
    //          += +1.00 t1_1_bb(a,i) w0
    //          += +1.00 <j,a||b,i>_abab t1_1_aa(b,j)
    //          += -1.00 d-_aa(j,b) t1_1_bb(a,i) t1_1_aa(b,j)
    //          += -1.00 d-_bb(j,b) t1_1_bb(a,i) t1_1_bb(b,j)
    //          += +1.00 f_bb(a,b) t1_1_bb(b,i)
    //          += +1.00 f_aa(j,b) t2_1_abab(b,a,j,i)
    //          += -0.50 <j,a||b,c>_bbbb t2_1_bbbb(b,c,i,j)
    //          += -0.50 <k,j||b,c>_abab t2_abab(b,a,k,j) t1_1_bb(c,i)
    //          += -0.50 <j,k||b,c>_abab t2_abab(b,a,j,k) t1_1_bb(c,i)
    //          += +1.00 <k,j||b,c>_bbbb t2_bbbb(b,a,i,j) t1_1_bb(c,k)
    //          += +2.00 d-_bb(j,b) t1_1_bb(a,j) t1_1_bb(b,i)
    //          += +1.00 d-_bb(j,i) t0_1 t1_1_bb(a,j)
    //          += -1.00 d-_bb(a,b) t0_1 t1_1_bb(b,i)
    //          += -1.00 <k,j||c,b>_abab t2_bbbb(b,a,i,j) t1_1_aa(c,k)
    //          += +1.00 <j,k||b,c>_abab t2_abab(b,a,j,i) t1_1_bb(c,k)
    rt1_1_bb("a,i") += 0.50 * tmps_["48_bb_vo"]("a,i");
    tmps_["48_bb_vo"].~TArrayD();

    // flops: o1v1  = o2v2 o2v2 o1v1
    //  mems: o1v1  = o1v1 o1v1 o1v1
    tmps_["53_aa_vo"]("a,i")  = -1.00 * dp["bb_ov"]("k,c") * t2_1["abab"]("a,c,i,k");
    tmps_["53_aa_vo"]("a,i") += dp["aa_ov"]("j,b") * t2_1["aaaa"]("b,a,i,j");

    // rt1_aa += +1.00 d-_aa(j,b) t2_1_aaaa(b,a,i,j)
    //        += -1.00 d-_bb(j,b) t2_1_abab(a,b,i,j)
    rt1_aa("a,i") += tmps_["53_aa_vo"]("a,i");

    // flops: o1v1  = o1v2 o2v1 o1v1
    //  mems: o1v1  = o1v1 o1v1 o1v1
    tmps_["54_aa_vo"]("a,i")  = dp["aa_vv"]("a,b") * t1_1["aa"]("b,i");
    tmps_["54_aa_vo"]("a,i") -= dp["aa_oo"]("j,i") * t1_1["aa"]("a,j");

    // rt1_aa += -1.00 d-_aa(a,b) t1_1_aa(b,i)
    //        += +1.00 d-_aa(j,i) t1_1_aa(a,j)
    rt1_aa("a,i") -= tmps_["54_aa_vo"]("a,i");

    // flops: o2v0  = o2v1
    //  mems: o2v0  = o2v0
    tmps_["52_aa_oo"]("i,j")  = t1_1["aa"]("a,i") * dp["aa_ov"]("j,a");

    // flops: o0v2  = o2v3
    //  mems: o0v2  = o0v2
    tmps_["51_aa_vv"]("a,b")  = t2["abab"]("a,c,i,j") * eri["abab_oovv"]("i,j,b,c");

    // flops: o1v1  = o2v2
    //  mems: o1v1  = o1v1
    tmps_["49_aa_vo"]("a,i")  = -1.00 * dp["bb_ov"]("j,b") * t2_1["abab"]("a,b,i,j");
    tmps_["49_aa_vo"].~TArrayD();

    // flops: o1v1  = o2v2
    //  mems: o1v1  = o1v1
    tmps_["50_aa_vo"]("a,i")  = -1.00 * f["bb_ov"]("j,b") * t2_1["abab"]("a,b,i,j");

    // flops: o1v1  = o2v3 o2v3 o2v2 o1v1 o1v1 o1v1 o1v2 o1v1 o2v2 o1v1 o1v1 o1v1 o1v1 o1v1 o2v1 o1v1 o2v2 o1v1 o3v2 o1v1 o3v2 o1v1 o1v1 o2v1 o2v2 o2v2 o1v1 o1v1 o1v1 o1v2 o1v1 o2v2 o1v1 o2v1 o1v1 o1v1 o1v1 o2v2 o1v1 o1v1 o1v1 o1v2 o1v1
    //  mems: o1v1  = o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1
    tmps_["55_aa_vo"]("a,i")  = -0.50 * eri["aaaa_vovv"]("a,k,c,b") * t2_1["aaaa"]("c,b,i,k");
    tmps_["55_aa_vo"]("a,i") -= eri["abab_vovv"]("a,l,c,e") * t2_1["abab"]("c,e,i,l");
    tmps_["55_aa_vo"]("a,i") += eri["aaaa_vovo"]("a,k,c,i") * t1_1["aa"]("c,k");
    tmps_["55_aa_vo"]("a,i") -= w0 * t1_1["aa"]("a,i");
    tmps_["55_aa_vo"]("a,i") -= f["aa_vv"]("a,c") * t1_1["aa"]("c,i");
    tmps_["55_aa_vo"]("a,i") += f["aa_ov"]("k,c") * t2_1["aaaa"]("c,a,i,k");
    tmps_["55_aa_vo"]("a,i") += scalars_["4"] * t1_1["aa"]("a,i");
    tmps_["55_aa_vo"]("a,i") += dp["aa_vo"]("a,i");
    tmps_["55_aa_vo"]("a,i") += t1_1["aa"]("a,k") * f["aa_oo"]("k,i");
    tmps_["55_aa_vo"]("a,i") += eri["abba_vovo"]("a,l,d,i") * t1_1["bb"]("d,l");
    tmps_["55_aa_vo"]("a,i") -= t2_1["abab"]("a,d,j,l") * eri["abba_oovo"]("j,l,d,i");
    tmps_["55_aa_vo"]("a,i") -= 0.50 * t2_1["aaaa"]("c,a,j,k") * eri["aaaa_oovo"]("k,j,c,i");
    tmps_["55_aa_vo"]("a,i") += tmps_["50_aa_vo"]("a,i");
    tmps_["55_aa_vo"]("a,i") -= 0.50 * t1_1["aa"]("a,j") * tmps_["11_aa_oo"]("i,j");
    tmps_["55_aa_vo"]("a,i") -= t1_1["bb"]("e,m") * tmps_["10_aabb_voov"]("a,i,m,e");
    tmps_["55_aa_vo"]("a,i") += t1_1["aa"]("b,j") * tmps_["9_aaaa_voov"]("a,i,j,b");
    tmps_["55_aa_vo"]("a,i") += t1_1["aa"]("b,i") * tmps_["51_aa_vv"]("a,b");
    tmps_["55_aa_vo"]("a,i") += t1_1["bb"]("e,m") * tmps_["2_aabb_voov"]("a,i,m,e");
    tmps_["55_aa_vo"]("a,i") -= 2.00 * t1_1["aa"]("a,k") * tmps_["52_aa_oo"]("i,k");
    tmps_["55_aa_vo"]("a,i") += t0_1 * tmps_["54_aa_vo"]("a,i");
    tmps_["55_aa_vo"]("a,i") -= t2["abab"]("a,d,i,l") * tmps_["44_bb_ov"]("l,d");
    tmps_["55_aa_vo"]("a,i") -= t0_1 * tmps_["53_aa_vo"]("a,i");
    tmps_["55_aa_vo"]("a,i") -= 0.50 * t1_1["aa"]("b,i") * tmps_["38_aa_vv"]("a,b");
    tmps_["50_aa_vo"].~TArrayD();
    tmps_["44_bb_ov"].~TArrayD();

    // rt1_1_aa += +1.00 <k,j||b,c>_aaaa t2_aaaa(b,a,i,j) t1_1_aa(c,k)
    //          += -1.00 <k,j||b,c>_bbbb t2_abab(a,b,i,j) t1_1_bb(c,k)
    //          += -0.50 <k,j||b,c>_aaaa t2_aaaa(b,c,i,j) t1_1_aa(a,k)
    //          += -0.50 <k,j||b,c>_abab t2_abab(b,c,i,j) t1_1_aa(a,k)
    //          += -0.50 <k,j||c,b>_abab t2_abab(c,b,i,j) t1_1_aa(a,k)
    //          += +1.00 <j,a||b,i>_aaaa t1_1_aa(b,j)
    //          += +0.50 <a,j||b,c>_abab t2_1_abab(b,c,i,j)
    //          += +0.50 <a,j||c,b>_abab t2_1_abab(c,b,i,j)
    //          += +1.00 t1_1_aa(a,i) w0
    //          += +1.00 f_aa(a,b) t1_1_aa(b,i)
    //          += -1.00 f_aa(j,b) t2_1_aaaa(b,a,i,j)
    //          += -1.00 d-_aa(j,b) t1_1_aa(a,i) t1_1_aa(b,j)
    //          += -1.00 d-_bb(j,b) t1_1_aa(a,i) t1_1_bb(b,j)
    //          += -0.50 <j,a||b,c>_aaaa t2_1_aaaa(b,c,i,j)
    //          += -1.00 d+_aa(a,i)
    //          += -1.00 f_aa(j,i) t1_1_aa(a,j)
    //          += +1.00 <a,j||i,b>_abab t1_1_bb(b,j)
    //          += -0.50 <k,j||i,b>_abab t2_1_abab(a,b,k,j)
    //          += -0.50 <j,k||i,b>_abab t2_1_abab(a,b,j,k)
    //          += -0.50 <k,j||b,i>_aaaa t2_1_aaaa(b,a,k,j)
    //          += +1.00 f_bb(j,b) t2_1_abab(a,b,i,j)
    //          += -0.50 <k,j||c,b>_abab t2_abab(a,b,k,j) t1_1_aa(c,i)
    //          += -0.50 <j,k||c,b>_abab t2_abab(a,b,j,k) t1_1_aa(c,i)
    //          += -1.00 <j,k||b,c>_abab t2_aaaa(b,a,i,j) t1_1_bb(c,k)
    //          += +2.00 d-_aa(j,b) t1_1_aa(a,j) t1_1_aa(b,i)
    //          += -1.00 d-_aa(a,b) t0_1 t1_1_aa(b,i)
    //          += +1.00 d-_aa(j,i) t0_1 t1_1_aa(a,j)
    //          += +1.00 <k,j||c,b>_abab t2_abab(a,b,i,j) t1_1_aa(c,k)
    //          += +1.00 d-_aa(j,b) t0_1 t2_1_aaaa(b,a,i,j)
    //          += -1.00 d-_bb(j,b) t0_1 t2_1_abab(a,b,i,j)
    //          += -0.50 <k,j||b,c>_aaaa t2_aaaa(b,a,k,j) t1_1_aa(c,i)
    rt1_1_aa("a,i") -= tmps_["55_aa_vo"]("a,i");
    tmps_["55_aa_vo"].~TArrayD();

    // flops: o2v2  = o3v2 o2v1 o3v2 o2v2
    //  mems: o2v2  = o2v2 o2v0 o2v2 o2v2
    tmps_["57_bbbb_ovvo"]("i,a,b,j")  = dp["bb_oo"]("k,i") * t2_1["bbbb"]("a,b,j,k");
    tmps_["57_bbbb_ovvo"]("i,a,b,j") -= dp["bb_ov"]("k,c") * t1_1["bb"]("c,j") * t2["bbbb"]("a,b,i,k");

    // rt2_bbbb += +1.00 P(i,j) d-_bb(k,j) t2_1_bbbb(a,b,i,k)
    //          += -1.00 P(i,j) d-_bb(k,c) t2_bbbb(a,b,j,k) t1_1_bb(c,i)
    rt2_bbbb("a,b,i,j") += tmps_["57_bbbb_ovvo"]("j,a,b,i");
    rt2_bbbb("a,b,i,j") -= tmps_["57_bbbb_ovvo"]("i,a,b,j");

    // flops: o2v0  = o3v1 o3v1 o2v0
    //  mems: o2v0  = o2v0 o2v0 o2v0
    tmps_["59_bb_oo"]("i,j")  = -1.00 * eri["abab_oovo"]("l,i,b,j") * t1_1["aa"]("b,l");
    tmps_["59_bb_oo"]("i,j") += eri["bbbb_oovo"]("i,k,a,j") * t1_1["bb"]("a,k");

    // flops: o2v0  = o3v2 o2v1 o3v2 o2v0 o2v0
    //  mems: o2v0  = o2v0 o2v0 o2v0 o2v0 o2v0
    tmps_["58_bb_oo"]("i,j")  = 0.50 * t2_1["bbbb"]("c,b,i,l") * eri["bbbb_oovv"]("j,l,c,b");
    tmps_["58_bb_oo"]("i,j") += t1_1["bb"]("c,i") * f["bb_ov"]("j,c");
    tmps_["58_bb_oo"]("i,j") += t2_1["abab"]("a,b,k,i") * eri["abab_oovv"]("k,j,a,b");

    // flops: o4v0  = o4v1
    //  mems: o4v0  = o4v0
    tmps_["56_bbbb_oooo"]("i,j,k,l")  = eri["bbbb_oovo"]("i,j,a,k") * t1_1["bb"]("a,l");

    // flops: o2v2  = o3v2 o3v2 o2v2 o3v2 o2v2 o2v2 o2v2 o2v3 o3v2 o2v2 o2v2 o3v2 o2v2 o4v2 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["60_bbbb_vvoo"]("a,b,i,j")  = 0.25 * t2_1["bbbb"]("a,b,i,k") * tmps_["16_bb_oo"]("j,k");
    tmps_["60_bbbb_vvoo"]("a,b,i,j") += t2_1["bbbb"]("a,b,i,l") * tmps_["45_bb_oo"]("j,l");
    tmps_["60_bbbb_vvoo"]("a,b,i,j") += 0.50 * t2["bbbb"]("a,b,i,l") * tmps_["59_bb_oo"]("l,j");
    tmps_["60_bbbb_vvoo"]("a,b,i,j") += 0.50 * t0_1 * tmps_["57_bbbb_ovvo"]("j,a,b,i");
    tmps_["60_bbbb_vvoo"]("a,b,i,j") += 0.50 * eri["bbbb_vvvo"]("a,b,c,j") * t1_1["bb"]("c,i");
    tmps_["60_bbbb_vvoo"]("a,b,i,j") -= 0.50 * f["bb_oo"]("l,j") * t2_1["bbbb"]("a,b,i,l");
    tmps_["60_bbbb_vvoo"]("a,b,i,j") += 0.50 * t2["bbbb"]("a,b,j,l") * tmps_["58_bb_oo"]("i,l");
    tmps_["60_bbbb_vvoo"]("a,b,i,j") -= 0.25 * t2["bbbb"]("a,b,k,l") * tmps_["56_bbbb_oooo"]("l,k,j,i");
    tmps_["57_bbbb_ovvo"].~TArrayD();
    tmps_["56_bbbb_oooo"].~TArrayD();

    // rt2_1_bbbb += +2.00 P(i,j) d-_bb(k,c) t2_1_bbbb(a,b,i,k) t1_1_bb(c,j)
    //            += -0.50 P(i,j) <l,k||c,d>_bbbb t2_bbbb(c,d,j,k) t2_1_bbbb(a,b,i,l)
    //            += -0.50 P(i,j) <k,l||c,d>_abab t2_abab(c,d,k,j) t2_1_bbbb(a,b,i,l)
    //            += -0.50 P(i,j) <k,l||d,c>_abab t2_abab(d,c,k,j) t2_1_bbbb(a,b,i,l)
    //            += -1.00 P(i,j) <l,k||c,j>_bbbb t2_bbbb(a,b,i,k) t1_1_bb(c,l)
    //            += -1.00 P(i,j) <l,k||c,j>_abab t2_bbbb(a,b,i,k) t1_1_aa(c,l)
    //            += +1.00 P(i,j) d-_bb(k,j) t0_1 t2_1_bbbb(a,b,i,k)
    //            += -1.00 P(i,j) d-_bb(k,c) t2_bbbb(a,b,j,k) t0_1 t1_1_bb(c,i)
    //            += +1.00 P(i,j) <a,b||c,j>_bbbb t1_1_bb(c,i)
    //            += -1.00 P(i,j) f_bb(k,j) t2_1_bbbb(a,b,i,k)
    //            += +0.50 P(i,j) <l,k||c,d>_abab t2_bbbb(a,b,j,k) t2_1_abab(c,d,l,i)
    //            += +0.50 P(i,j) <l,k||d,c>_abab t2_bbbb(a,b,j,k) t2_1_abab(d,c,l,i)
    //            += +1.00 P(i,j) f_bb(k,c) t2_bbbb(a,b,j,k) t1_1_bb(c,i)
    //            += -0.50 P(i,j) <l,k||c,d>_bbbb t2_bbbb(a,b,j,k) t2_1_bbbb(c,d,i,l)
    //            += +0.50 P(i,j) <l,k||c,j>_bbbb t2_bbbb(a,b,l,k) t1_1_bb(c,i)
    rt2_1_bbbb("a,b,i,j") += 2.00 * tmps_["60_bbbb_vvoo"]("a,b,i,j");
    rt2_1_bbbb("a,b,i,j") -= 2.00 * tmps_["60_bbbb_vvoo"]("a,b,j,i");
    tmps_["60_bbbb_vvoo"].~TArrayD();

    // flops: o2v2  = o3v2 o3v2 o2v3 o2v2
    //  mems: o2v2  = o3v1 o2v2 o2v2 o2v2
    tmps_["62_bbbb_voov"]("a,i,j,b")  = dp["bb_ov"]("k,c") * t2["bbbb"]("c,a,i,j") * t1_1["bb"]("b,k");
    tmps_["62_bbbb_voov"]("a,i,j,b") += dp["bb_vv"]("a,c") * t2_1["bbbb"]("c,b,i,j");

    // rt2_bbbb += -1.00 P(a,b) d-_bb(k,c) t2_bbbb(c,a,i,j) t1_1_bb(b,k)
    //          += -1.00 P(a,b) d-_bb(a,c) t2_1_bbbb(c,b,i,j)
    rt2_bbbb("a,b,i,j") -= tmps_["62_bbbb_voov"]("a,i,j,b");
    rt2_bbbb("a,b,i,j") += tmps_["62_bbbb_voov"]("b,i,j,a");

    // flops: o0v2  = o1v3 o1v3 o0v2
    //  mems: o0v2  = o0v2 o0v2 o0v2
    tmps_["65_bb_vv"]("a,b")  = -1.00 * eri["baab_vovv"]("a,j,d,b") * t1_1["aa"]("d,j");
    tmps_["65_bb_vv"]("a,b") += eri["bbbb_vovv"]("a,i,b,c") * t1_1["bb"]("c,i");

    // flops: o0v2  = o2v3 o2v3 o0v2
    //  mems: o0v2  = o0v2 o0v2 o0v2
    tmps_["64_bb_vv"]("a,b")  = 2.00 * t2_1["abab"]("d,a,k,j") * eri["abab_oovv"]("k,j,d,b");
    tmps_["64_bb_vv"]("a,b") += t2_1["bbbb"]("c,a,i,j") * eri["bbbb_oovv"]("j,i,b,c");
}
