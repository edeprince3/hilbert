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

void QED_CCSD::resid_21_4() {

    // unpack integrals
    TArrayMap &Id         = Id_blks_;
    TArrayMap &f          = F_blks_;
    TArrayMap &eri        = V_blks_;

    // get effective dipole integrals
    TArrayMap dp = effective_dipole();

    // extract scalars_
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


    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["95_bbbb_voov"]("a,i,j,b")  = t2_1["abab"]("c,a,k,i") * eri["abab_oovv"]("k,j,c,b");

    // flops: o2v0  = o2v1
    //  mems: o2v0  = o2v0
    tmps_["92_bb_oo"]("i,j")  = dp["bb_ov"]("i,a") * t1_1["bb"]("a,j");
    tmps_["92_bb_oo"].~TArrayD();

    // flops: o2v0  = o2v1
    //  mems: o2v0  = o2v0
    tmps_["93_aa_oo"]("i,j")  = dp["aa_ov"]("i,a") * t1_1["aa"]("a,j");
    tmps_["93_aa_oo"].~TArrayD();

    // flops: o3v1  = o4v2 o4v2 o3v1
    //  mems: o3v1  = o3v1 o3v1 o3v1
    tmps_["94_bbbb_oovo"]("i,j,a,k")  = -1.00 * eri["abab_oovo"]("m,i,c,j") * t2["abab"]("c,a,m,k");
    tmps_["94_bbbb_oovo"]("i,j,a,k") += eri["bbbb_oovo"]("l,i,b,j") * t2["bbbb"]("b,a,k,l");

    // flops: o2v2  = o3v3 o3v3 o3v3 o2v2 o3v3 o3v3 o2v2 o2v2 o2v2 o2v2 o3v2 o2v2 o2v2 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["98_bbbb_vovo"]("a,i,b,j")  = -1.00 * t2["bbbb"]("d,a,i,l") * tmps_["97_bbbb_vovo"]("b,l,d,j");
    tmps_["98_bbbb_vovo"]("a,i,b,j") -= tmps_["4_bbbb_voov"]("b,i,n,e") * t2_1["bbbb"]("e,a,j,n");
    tmps_["98_bbbb_vovo"]("a,i,b,j") += tmps_["14_bbaa_voov"]("b,i,m,f") * t2_1["abab"]("f,a,m,j");
    tmps_["98_bbbb_vovo"]("a,i,b,j") += eri["bbbb_vovo"]("b,l,d,i") * t2_1["bbbb"]("d,a,j,l");
    tmps_["98_bbbb_vovo"]("a,i,b,j") -= eri["baab_vovo"]("b,k,c,i") * t2_1["abab"]("c,a,k,j");
    tmps_["98_bbbb_vovo"]("a,i,b,j") += t1_1["bb"]("a,j") * tmps_["47_bb_ov"]("i,b");
    tmps_["98_bbbb_vovo"]("a,i,b,j") -= t1_1["bb"]("a,n") * tmps_["94_bbbb_oovo"]("n,i,b,j");
    tmps_["98_bbbb_vovo"]("a,i,b,j") += t1_1["bb"]("b,i") * tmps_["46_bb_vo"]("a,j");
    tmps_["98_bbbb_vovo"]("a,i,b,j") += tmps_["15_bbbb_voov"]("b,i,n,e") * t2_1["bbbb"]("e,a,j,n");
    tmps_["98_bbbb_vovo"]("a,i,b,j") -= t2["bbbb"]("d,b,i,l") * tmps_["95_bbbb_voov"]("a,j,l,d");
    tmps_["98_bbbb_vovo"]("a,i,b,j") += t2["abab"]("c,a,k,i") * tmps_["96_baab_vovo"]("b,k,c,j");
    tmps_["94_bbbb_oovo"].~TArrayD();

    // rt2_1_bbbb += +1.00 P(i,j) P(a,b) <k,a||c,d>_abab t2_abab(c,b,k,j) t1_1_bb(d,i)
    //            += +1.00 P(i,j) P(a,b) <k,l||c,d>_abab t2_abab(c,a,k,j) t2_1_bbbb(d,b,i,l)
    //            += +1.00 P(i,j) P(a,b) <l,k||c,d>_aaaa t2_abab(c,a,k,j) t2_1_abab(d,b,l,i)
    //            += +1.00 P(i,j) P(a,b) <k,a||c,j>_bbbb t2_1_bbbb(c,b,i,k)
    //            += -1.00 P(i,j) P(a,b) <k,a||c,j>_abab t2_1_abab(c,b,k,i)
    //            += -1.00 P(i,j) P(a,b) d-_bb(k,j) t1_1_bb(a,k) t1_1_bb(b,i)
    //            += +1.00 P(i,j) P(a,b) d-_bb(a,c) t1_1_bb(b,i) t1_1_bb(c,j)
    //            += -1.00 P(i,j) P(a,b) <l,k||c,j>_bbbb t2_bbbb(c,a,i,k) t1_1_bb(b,l)
    //            += -1.00 P(i,j) P(a,b) <k,l||c,j>_abab t2_abab(c,a,k,i) t1_1_bb(b,l)
    //            += -1.00 P(i,j) P(a,b) d-_bb(k,c) t1_1_bb(a,j) t2_1_bbbb(c,b,i,k)
    //            += +1.00 P(i,j) P(a,b) d-_aa(k,c) t1_1_bb(a,j) t2_1_abab(c,b,k,i)
    //            += +1.00 P(i,j) P(a,b) <l,k||c,d>_bbbb t2_bbbb(c,a,j,k) t2_1_bbbb(d,b,i,l)
    //            += +1.00 P(i,j) P(a,b) <l,k||d,c>_abab t2_bbbb(c,a,j,k) t2_1_abab(d,b,l,i)
    //            += -1.00 P(i,j) P(a,b) <k,a||c,d>_bbbb t2_bbbb(c,b,j,k) t1_1_bb(d,i)
    rt2_1_bbbb("a,b,i,j") -= tmps_["98_bbbb_vovo"]("b,j,a,i");
    rt2_1_bbbb("a,b,i,j") += tmps_["98_bbbb_vovo"]("b,i,a,j");
    rt2_1_bbbb("a,b,i,j") += tmps_["98_bbbb_vovo"]("a,j,b,i");
    rt2_1_bbbb("a,b,i,j") -= tmps_["98_bbbb_vovo"]("a,i,b,j");
    tmps_["98_bbbb_vovo"].~TArrayD();

    // flops: o2v2  = o2v3 o2v3 o3v2 o3v2 o2v1 o3v2 o2v2 o2v2 o2v2 o3v2 o2v2 o2v1 o3v2 o2v2 o3v2 o3v2 o2v2 o3v2 o2v2
    //  mems: o2v2  = o2v2 o2v2 o3v1 o2v2 o2v0 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v0 o2v2 o2v2 o3v1 o2v2 o2v2 o2v2 o2v2
    tmps_["99_baba_ovvo"]("i,a,b,j")  = -1.00 * dp["aa_vv"]("a,d") * t2_1["abab"]("d,b,j,i");
    tmps_["99_baba_ovvo"]("i,a,b,j") -= dp["bb_vv"]("b,c") * t2_1["abab"]("a,c,j,i");
    tmps_["99_baba_ovvo"]("i,a,b,j") += dp["aa_ov"]("l,d") * t2["abab"]("d,b,j,i") * t1_1["aa"]("a,l");
    tmps_["99_baba_ovvo"]("i,a,b,j") += dp["bb_ov"]("k,c") * t1_1["bb"]("c,i") * t2["abab"]("a,b,j,k");
    tmps_["99_baba_ovvo"]("i,a,b,j") += dp["aa_oo"]("l,j") * t2_1["abab"]("a,b,l,i");
    tmps_["99_baba_ovvo"]("i,a,b,j") += dp["aa_ov"]("l,d") * t1_1["aa"]("d,j") * t2["abab"]("a,b,l,i");
    tmps_["99_baba_ovvo"]("i,a,b,j") += dp["bb_ov"]("k,c") * t2["abab"]("a,c,j,i") * t1_1["bb"]("b,k");
    tmps_["99_baba_ovvo"]("i,a,b,j") += dp["bb_oo"]("k,i") * t2_1["abab"]("a,b,j,k");

    // rt2_abab += +1.00 d-_bb(k,c) t2_abab(a,b,i,k) t1_1_bb(c,j)
    //          += +1.00 d-_aa(k,c) t2_abab(c,b,i,j) t1_1_aa(a,k)
    //          += -1.00 d-_bb(b,c) t2_1_abab(a,c,i,j)
    //          += -1.00 d-_aa(a,c) t2_1_abab(c,b,i,j)
    //          += +1.00 d-_aa(k,i) t2_1_abab(a,b,k,j)
    //          += +1.00 d-_aa(k,c) t2_abab(a,b,k,j) t1_1_aa(c,i)
    //          += +1.00 d-_bb(k,c) t2_abab(a,c,i,j) t1_1_bb(b,k)
    //          += +1.00 d-_bb(k,j) t2_1_abab(a,b,i,k)
    rt2_abab("a,b,i,j") += tmps_["99_baba_ovvo"]("j,a,b,i");

    // flops: o2v2  = o2v3
    //  mems: o2v2  = o2v2
    tmps_["108_abba_vovo"]("a,i,b,j")  = eri["abab_vovv"]("a,i,c,b") * t1_1["aa"]("c,j");

    // flops: o2v2  = o2v3
    //  mems: o2v2  = o2v2
    tmps_["107_aaaa_vovo"]("a,i,b,j")  = eri["aaaa_vovv"]("a,i,b,c") * t1_1["aa"]("c,j");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["106_aabb_voov"]("a,i,j,b")  = t2_1["aaaa"]("c,a,i,k") * eri["abab_oovv"]("k,j,c,b");

    // flops: o2v2  = o3v3 o2v3 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2
    tmps_["105_baab_voov"]("a,i,j,b")  = t2_1["abab"]("c,a,i,k") * eri["abab_oovv"]("j,k,c,b");
    tmps_["105_baab_voov"]("a,i,j,b") += eri["baab_vovv"]("a,j,c,b") * t1_1["aa"]("c,i");

    // flops: o3v1  = o3v2 o3v3 o3v1 o3v2 o3v1
    //  mems: o3v1  = o3v1 o3v1 o3v1 o3v1 o3v1
    tmps_["104_abab_vooo"]("a,i,j,k")  = -2.00 * t2_1["abab"]("a,b,j,k") * dp["bb_ov"]("i,b");
    tmps_["104_abab_vooo"]("a,i,j,k") += eri["abab_vovv"]("a,i,c,d") * t2["abab"]("c,d,j,k");
    tmps_["104_abab_vooo"]("a,i,j,k") += t2["abab"]("a,b,j,k") * f["bb_ov"]("i,b");

    // flops: o3v1  = o4v2 o4v2 o3v1 o4v2 o3v1
    //  mems: o3v1  = o3v1 o3v1 o3v1 o3v1 o3v1
    tmps_["103_bbaa_oovo"]("i,j,a,k")  = eri["abba_oovo"]("m,i,b,k") * t2["abab"]("a,b,m,j");
    tmps_["103_bbaa_oovo"]("i,j,a,k") += eri["bbbb_oovo"]("l,i,b,j") * t2["abab"]("a,b,k,l");
    tmps_["103_bbaa_oovo"]("i,j,a,k") -= eri["abab_oovo"]("m,i,c,j") * t2["aaaa"]("c,a,k,m");

    // flops: o3v1  = o4v2 o4v2 o3v1 o4v2 o3v1
    //  mems: o3v1  = o3v1 o3v1 o3v1 o3v1 o3v1
    tmps_["102_aabb_oovo"]("i,j,a,k")  = eri["aaaa_oovo"]("l,i,b,j") * t2["abab"]("b,a,l,k");
    tmps_["102_aabb_oovo"]("i,j,a,k") += eri["abba_oovo"]("i,m,c,j") * t2["bbbb"]("c,a,k,m");
    tmps_["102_aabb_oovo"]("i,j,a,k") -= eri["abab_oovo"]("i,m,b,k") * t2["abab"]("b,a,j,m");

    // flops: o3v1  = o3v2 o3v2 o3v3 o3v1 o3v1
    //  mems: o3v1  = o3v1 o3v1 o3v1 o3v1 o3v1
    tmps_["101_baba_vooo"]("a,i,j,k")  = -2.00 * t2_1["abab"]("b,a,i,j") * dp["aa_ov"]("k,b");
    tmps_["101_baba_vooo"]("a,i,j,k") += t2["abab"]("b,a,i,j") * f["aa_ov"]("k,b");
    tmps_["101_baba_vooo"]("a,i,j,k") -= eri["baab_vovv"]("a,k,b,c") * t2["abab"]("b,c,i,j");

    // flops: o4v0  = o4v1 o4v1 o4v2 o4v0 o4v0
    //  mems: o4v0  = o4v0 o4v0 o4v0 o4v0 o4v0
    tmps_["100_abba_oooo"]("i,j,k,l")  = -1.00 * eri["abba_oovo"]("i,j,c,l") * t1_1["bb"]("c,k");
    tmps_["100_abba_oooo"]("i,j,k,l") += eri["abab_oovo"]("i,j,a,k") * t1_1["aa"]("a,l");
    tmps_["100_abba_oooo"]("i,j,k,l") += eri["abab_oovv"]("i,j,a,b") * t2_1["abab"]("a,b,l,k");

    // flops: o2v2  = o2v2 o3v3 o2v2 o2v3 o2v2 o4v2 o2v2 o3v3 o2v2 o3v3 o2v2 o3v2 o2v2 o2v2 o2v2 o2v3 o2v2 o2v4 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o2v3 o2v2 o3v3 o2v2 o3v3 o2v2 o2v3 o3v3 o2v2 o3v3 o2v2 o2v3 o2v2 o3v2 o3v2 o2v2 o3v2 o3v2 o3v3 o3v2 o3v2 o3v2 o3v2 o3v3 o3v3 o3v3 o2v2 o2v2 o2v2 o2v2 o3v3 o2v2 o3v3 o2v2 o2v2 o2v2 o2v2 o2v3 o2v2 o2v3 o2v2 o4v2 o2v2 o3v2 o2v2 o3v2 o2v2 o2v2 o3v3 o2v2 o2v2 o2v2 o2v3 o2v2 o3v3 o2v2 o3v3 o2v2 o4v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v3 o2v2 o2v3 o2v2 o3v3 o2v2 o2v3 o2v2 o2v3 o2v2 o3v3 o2v2 o2v2 o2v2 o2v2 o3v3 o2v2 o2v2 o2v2 o3v2 o2v2 o3v2 o2v2 o2v2 o2v2 o2v3 o2v2 o3v3 o2v2 o3v3 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["109_abba_vovo"]("a,i,b,j")  = -1.00 * scalars_["4"] * t2_1["abab"]("a,b,j,i");
    tmps_["109_abba_vovo"]("a,i,b,j") -= eri["abab_vovo"]("a,l,c,i") * t2_1["abab"]("c,b,j,l");
    tmps_["109_abba_vovo"]("a,i,b,j") += eri["abab_vvvo"]("a,b,c,i") * t1_1["aa"]("c,j");
    tmps_["109_abba_vovo"]("a,i,b,j") += eri["abab_oooo"]("m,l,j,i") * t2_1["abab"]("a,b,m,l");
    tmps_["109_abba_vovo"]("a,i,b,j") -= eri["aaaa_vovo"]("a,k,c,j") * t2_1["abab"]("c,b,k,i");
    tmps_["109_abba_vovo"]("a,i,b,j") += eri["abba_vovo"]("a,l,d,j") * t2_1["bbbb"]("d,b,i,l");
    tmps_["109_abba_vovo"]("a,i,b,j") -= f["aa_oo"]("k,j") * t2_1["abab"]("a,b,k,i");
    tmps_["109_abba_vovo"]("a,i,b,j") += w0 * t2_1["abab"]("a,b,j,i");
    tmps_["109_abba_vovo"]("a,i,b,j") -= eri["abba_vvvo"]("a,b,d,j") * t1_1["bb"]("d,i");
    tmps_["109_abba_vovo"]("a,i,b,j") += eri["abab_vvvv"]("a,b,c,f") * t2_1["abab"]("c,f,j,i");
    tmps_["109_abba_vovo"]("a,i,b,j") -= f["bb_oo"]("l,i") * t2_1["abab"]("a,b,j,l");
    tmps_["109_abba_vovo"]("a,i,b,j") -= eri["abab_vooo"]("a,l,j,i") * t1_1["bb"]("b,l");
    tmps_["109_abba_vovo"]("a,i,b,j") += eri["baab_vooo"]("b,k,j,i") * t1_1["aa"]("a,k");
    tmps_["109_abba_vovo"]("a,i,b,j") += f["aa_vv"]("a,c") * t2_1["abab"]("c,b,j,i");
    tmps_["109_abba_vovo"]("a,i,b,j") += eri["baab_vovo"]("b,k,c,i") * t2_1["aaaa"]("c,a,j,k");
    tmps_["109_abba_vovo"]("a,i,b,j") -= eri["baba_vovo"]("b,k,d,j") * t2_1["abab"]("a,d,k,i");
    tmps_["109_abba_vovo"]("a,i,b,j") -= eri["abab_vovv"]("a,l,c,f") * t1_1["bb"]("f,i") * t2["abab"]("c,b,j,l");
    tmps_["109_abba_vovo"]("a,i,b,j") -= eri["bbbb_vovo"]("b,l,d,i") * t2_1["abab"]("a,d,j,l");
    tmps_["109_abba_vovo"]("a,i,b,j") += f["bb_vv"]("b,d") * t2_1["abab"]("a,d,j,i");
    tmps_["109_abba_vovo"]("a,i,b,j") -= t1_1["aa"]("a,k") * tmps_["101_baba_vooo"]("b,j,i,k");
    tmps_["109_abba_vovo"]("a,i,b,j") -= t1_1["aa"]("a,m") * tmps_["102_aabb_oovo"]("m,j,b,i");
    tmps_["109_abba_vovo"]("a,i,b,j") += t1_1["aa"]("a,j") * tmps_["46_bb_vo"]("b,i");
    tmps_["109_abba_vovo"]("a,i,b,j") -= t1_1["bb"]("b,l") * tmps_["104_abab_vooo"]("a,l,j,i");
    tmps_["109_abba_vovo"]("a,i,b,j") -= t2["abab"]("a,b,k,i") * tmps_["89_aa_oo"]("j,k");
    tmps_["109_abba_vovo"]("a,i,b,j") -= t2_1["abab"]("a,f,j,n") * tmps_["15_bbbb_voov"]("b,i,n,f");
    tmps_["109_abba_vovo"]("a,i,b,j") += 0.50 * t2_1["abab"]("a,b,j,n") * tmps_["16_bb_oo"]("i,n");
    tmps_["109_abba_vovo"]("a,i,b,j") += 2.00 * t2_1["abab"]("a,b,j,l") * tmps_["45_bb_oo"]("i,l");
    tmps_["109_abba_vovo"]("a,i,b,j") += 0.50 * t2_1["abab"]("a,b,m,i") * tmps_["11_aa_oo"]("j,m");
    tmps_["109_abba_vovo"]("a,i,b,j") += t2["abab"]("a,b,j,l") * tmps_["59_bb_oo"]("l,i");
    tmps_["109_abba_vovo"]("a,i,b,j") -= t2_1["bbbb"]("f,b,i,n") * tmps_["10_aabb_voov"]("a,j,n,f");
    tmps_["109_abba_vovo"]("a,i,b,j") += t2["bbbb"]("d,b,i,l") * tmps_["106_aabb_voov"]("a,j,l,d");
    tmps_["109_abba_vovo"]("a,i,b,j") += t2_1["abab"]("a,f,m,i") * tmps_["75_baab_voov"]("b,j,m,f");
    tmps_["109_abba_vovo"]("a,i,b,j") += t2["abab"]("a,d,j,l") * tmps_["95_bbbb_voov"]("b,i,l,d");
    tmps_["109_abba_vovo"]("a,i,b,j") += t2["abab"]("a,d,k,i") * tmps_["105_baab_voov"]("b,j,k,d");
    tmps_["109_abba_vovo"]("a,i,b,j") += t1_1["aa"]("a,j") * tmps_["47_bb_ov"]("i,b");
    tmps_["109_abba_vovo"]("a,i,b,j") += t2["abab"]("c,b,j,i") * tmps_["82_aa_vv"]("a,c");
    tmps_["109_abba_vovo"]("a,i,b,j") += t2["abab"]("a,d,j,i") * tmps_["65_bb_vv"]("b,d");
    tmps_["109_abba_vovo"]("a,i,b,j") += t2_1["abab"]("a,b,m,l") * tmps_["76_abab_oooo"]("j,i,m,l");
    tmps_["109_abba_vovo"]("a,i,b,j") += 2.00 * t2_1["abab"]("a,b,k,i") * tmps_["52_aa_oo"]("j,k");
    tmps_["109_abba_vovo"]("a,i,b,j") += t2["abab"]("a,b,k,i") * tmps_["90_aa_oo"]("k,j");
    tmps_["109_abba_vovo"]("a,i,b,j") += t2_1["abab"]("a,f,j,n") * tmps_["4_bbbb_voov"]("b,i,n,f");
    tmps_["109_abba_vovo"]("a,i,b,j") -= t2_1["abab"]("e,b,j,i") * tmps_["51_aa_vv"]("a,e");
    tmps_["109_abba_vovo"]("a,i,b,j") -= t2_1["abab"]("e,b,m,i") * tmps_["9_aaaa_voov"]("a,j,m,e");
    tmps_["109_abba_vovo"]("a,i,b,j") += t2_1["bbbb"]("f,b,i,n") * tmps_["2_aabb_voov"]("a,j,n,f");
    tmps_["109_abba_vovo"]("a,i,b,j") += t2["abab"]("a,b,m,l") * tmps_["100_abba_oooo"]("m,l,i,j");
    tmps_["109_abba_vovo"]("a,i,b,j") -= t2_1["abab"]("a,f,j,i") * tmps_["6_bb_vv"]("b,f");
    tmps_["109_abba_vovo"]("a,i,b,j") -= 0.50 * t2["abab"]("a,d,j,i") * tmps_["64_bb_vv"]("b,d");
    tmps_["109_abba_vovo"]("a,i,b,j") += t2["aaaa"]("c,a,j,k") * tmps_["96_baab_vovo"]("b,k,c,i");
    tmps_["109_abba_vovo"]("a,i,b,j") -= tmps_["81_aa_vv"]("a,c") * t2["abab"]("c,b,j,i");
    tmps_["109_abba_vovo"]("a,i,b,j") += 0.50 * t2_1["abab"]("a,f,j,i") * tmps_["43_bb_vv"]("b,f");
    tmps_["109_abba_vovo"]("a,i,b,j") -= t2["abab"]("a,d,j,l") * tmps_["97_bbbb_vovo"]("b,l,d,i");
    tmps_["109_abba_vovo"]("a,i,b,j") += t0_1 * tmps_["99_baba_ovvo"]("i,a,b,j");
    tmps_["109_abba_vovo"]("a,i,b,j") -= t2_1["aaaa"]("e,a,j,m") * tmps_["14_bbaa_voov"]("b,i,m,e");
    tmps_["109_abba_vovo"]("a,i,b,j") -= t1_1["bb"]("b,i") * tmps_["54_aa_vo"]("a,j");
    tmps_["109_abba_vovo"]("a,i,b,j") -= t1_1["bb"]("b,n") * tmps_["103_bbaa_oovo"]("n,i,a,j");
    tmps_["109_abba_vovo"]("a,i,b,j") -= t2["abab"]("a,b,j,l") * tmps_["58_bb_oo"]("i,l");
    tmps_["109_abba_vovo"]("a,i,b,j") += t1_1["bb"]("b,i") * tmps_["53_aa_vo"]("a,j");
    tmps_["109_abba_vovo"]("a,i,b,j") += 0.50 * t2_1["abab"]("e,b,j,i") * tmps_["38_aa_vv"]("a,e");
    tmps_["109_abba_vovo"]("a,i,b,j") -= t2["bbbb"]("d,b,i,l") * tmps_["108_abba_vovo"]("a,l,d,j");
    tmps_["109_abba_vovo"]("a,i,b,j") -= t2["abab"]("c,b,k,i") * tmps_["107_aaaa_vovo"]("a,k,c,j");
    tmps_["105_baab_voov"].~TArrayD();
    tmps_["104_abab_vooo"].~TArrayD();
    tmps_["103_bbaa_oovo"].~TArrayD();
    tmps_["102_aabb_oovo"].~TArrayD();
    tmps_["101_baba_vooo"].~TArrayD();
    tmps_["100_abba_oooo"].~TArrayD();
    tmps_["99_baba_ovvo"].~TArrayD();
    tmps_["97_bbbb_vovo"].~TArrayD();
    tmps_["96_baab_vovo"].~TArrayD();
    tmps_["95_bbbb_voov"].~TArrayD();
    tmps_["90_aa_oo"].~TArrayD();
    tmps_["89_aa_oo"].~TArrayD();
    tmps_["82_aa_vv"].~TArrayD();
    tmps_["81_aa_vv"].~TArrayD();
    tmps_["76_abab_oooo"].~TArrayD();
    tmps_["75_baab_voov"].~TArrayD();
    tmps_["65_bb_vv"].~TArrayD();
    tmps_["64_bb_vv"].~TArrayD();
    tmps_["59_bb_oo"].~TArrayD();
    tmps_["58_bb_oo"].~TArrayD();
    tmps_["52_aa_oo"].~TArrayD();
    tmps_["51_aa_vv"].~TArrayD();
    tmps_["47_bb_ov"].~TArrayD();
    tmps_["46_bb_vo"].~TArrayD();
    tmps_["45_bb_oo"].~TArrayD();
    tmps_["43_bb_vv"].~TArrayD();
    tmps_["38_aa_vv"].~TArrayD();
    tmps_["16_bb_oo"].~TArrayD();
    tmps_["15_bbbb_voov"].~TArrayD();
    tmps_["14_bbaa_voov"].~TArrayD();
    tmps_["11_aa_oo"].~TArrayD();
    tmps_["6_bb_vv"].~TArrayD();
    tmps_["4_bbbb_voov"].~TArrayD();

    // rt2_1_abab += +1.00 <l,k||c,d>_abab t2_abab(c,b,i,k) t2_1_abab(a,d,l,j)
    //            += +1.00 <l,k||d,c>_abab t2_bbbb(c,b,j,k) t2_1_aaaa(d,a,i,l)
    //            += +1.00 <l,k||c,d>_bbbb t2_abab(a,c,i,k) t2_1_bbbb(d,b,j,l)
    //            += -1.00 <l,k||c,j>_bbbb t2_abab(a,b,i,k) t1_1_bb(c,l)
    //            += -1.00 <l,k||c,j>_abab t2_abab(a,b,i,k) t1_1_aa(c,l)
    //            += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,d,i,k) t2_1_abab(a,b,l,j)
    //            += -0.50 <l,k||c,d>_abab t2_abab(c,d,i,k) t2_1_abab(a,b,l,j)
    //            += -0.50 <l,k||d,c>_abab t2_abab(d,c,i,k) t2_1_abab(a,b,l,j)
    //            += +1.00 <l,k||d,c>_abab t2_abab(a,c,i,k) t2_1_abab(d,b,l,j)
    //            += +1.00 <k,l||d,c>_abab t2_abab(a,c,k,j) t2_1_abab(d,b,i,l)
    //            += -1.00 <k,b||d,c>_abab t2_abab(a,c,k,j) t1_1_aa(d,i)
    //            += +1.00 d-_bb(k,j) t1_1_aa(a,i) t1_1_bb(b,k)
    //            += -1.00 d-_bb(b,c) t1_1_aa(a,i) t1_1_bb(c,j)
    //            += +2.00 d-_bb(k,c) t2_1_abab(a,b,i,k) t1_1_bb(c,j)
    //            += -1.00 <k,a||c,d>_aaaa t2_abab(c,b,i,j) t1_1_aa(d,k)
    //            += +1.00 <a,k||c,d>_abab t2_abab(c,b,i,j) t1_1_bb(d,k)
    //            += -1.00 <k,b||c,d>_bbbb t2_abab(a,c,i,j) t1_1_bb(d,k)
    //            += +1.00 <k,b||d,c>_abab t2_abab(a,c,i,j) t1_1_aa(d,k)
    //            += +0.25 <l,k||c,d>_abab t2_abab(c,d,i,j) t2_1_abab(a,b,l,k)
    //            += +0.25 <k,l||c,d>_abab t2_abab(c,d,i,j) t2_1_abab(a,b,k,l)
    //            += +0.25 <l,k||d,c>_abab t2_abab(d,c,i,j) t2_1_abab(a,b,l,k)
    //            += +0.25 <k,l||d,c>_abab t2_abab(d,c,i,j) t2_1_abab(a,b,k,l)
    //            += +2.00 d-_aa(k,c) t2_1_abab(a,b,k,j) t1_1_aa(c,i)
    //            += -1.00 <k,l||i,c>_abab t2_abab(a,b,k,j) t1_1_bb(c,l)
    //            += -1.00 <l,k||c,i>_aaaa t2_abab(a,b,k,j) t1_1_aa(c,l)
    //            += -0.50 <l,k||c,d>_bbbb t2_bbbb(c,d,j,k) t2_1_abab(a,b,i,l)
    //            += -0.50 <k,l||c,d>_abab t2_abab(c,d,k,j) t2_1_abab(a,b,i,l)
    //            += -0.50 <k,l||d,c>_abab t2_abab(d,c,k,j) t2_1_abab(a,b,i,l)
    //            += +1.00 <k,l||c,d>_abab t2_abab(c,b,k,j) t2_1_abab(a,d,i,l)
    //            += +1.00 <l,k||c,d>_bbbb t2_bbbb(c,b,j,k) t2_1_abab(a,d,i,l)
    //            += -0.50 <k,l||c,d>_abab t2_abab(a,b,k,j) t2_1_abab(c,d,i,l)
    //            += -0.50 <k,l||d,c>_abab t2_abab(a,b,k,j) t2_1_abab(d,c,i,l)
    //            += -1.00 f_aa(k,c) t2_abab(a,b,k,j) t1_1_aa(c,i)
    //            += +0.50 <l,k||c,d>_aaaa t2_abab(a,b,k,j) t2_1_aaaa(c,d,i,l)
    //            += -0.50 <l,k||d,c>_abab t2_abab(a,c,l,k) t2_1_abab(d,b,i,j)
    //            += -0.50 <k,l||d,c>_abab t2_abab(a,c,k,l) t2_1_abab(d,b,i,j)
    //            += +1.00 <l,k||c,d>_aaaa t2_aaaa(c,a,i,k) t2_1_abab(d,b,l,j)
    //            += +1.00 <k,l||c,d>_abab t2_aaaa(c,a,i,k) t2_1_bbbb(d,b,j,l)
    //            += +0.50 <l,k||c,j>_abab t2_abab(a,b,l,k) t1_1_aa(c,i)
    //            += +0.50 <k,l||c,j>_abab t2_abab(a,b,k,l) t1_1_aa(c,i)
    //            += +0.25 <l,k||c,d>_abab t2_abab(a,b,l,k) t2_1_abab(c,d,i,j)
    //            += +0.25 <l,k||d,c>_abab t2_abab(a,b,l,k) t2_1_abab(d,c,i,j)
    //            += +0.25 <k,l||c,d>_abab t2_abab(a,b,k,l) t2_1_abab(c,d,i,j)
    //            += +0.25 <k,l||d,c>_abab t2_abab(a,b,k,l) t2_1_abab(d,c,i,j)
    //            += +0.50 <l,k||i,c>_abab t2_abab(a,b,l,k) t1_1_bb(c,j)
    //            += +0.50 <k,l||i,c>_abab t2_abab(a,b,k,l) t1_1_bb(c,j)
    //            += -0.50 <a,k||c,d>_abab t2_abab(c,d,i,j) t1_1_bb(b,k)
    //            += -0.50 <a,k||d,c>_abab t2_abab(d,c,i,j) t1_1_bb(b,k)
    //            += +2.00 d-_bb(k,c) t2_1_abab(a,c,i,j) t1_1_bb(b,k)
    //            += -1.00 f_bb(k,c) t2_abab(a,c,i,j) t1_1_bb(b,k)
    //            += +1.00 d-_bb(k,c) t1_1_aa(a,i) t2_1_bbbb(c,b,j,k)
    //            += -1.00 d-_aa(k,c) t1_1_aa(a,i) t2_1_abab(c,b,k,j)
    //            += +1.00 <l,k||c,i>_aaaa t2_abab(c,b,k,j) t1_1_aa(a,l)
    //            += +1.00 <l,k||i,c>_abab t2_bbbb(c,b,j,k) t1_1_aa(a,l)
    //            += +1.00 <l,k||c,j>_abab t2_abab(c,b,i,k) t1_1_aa(a,l)
    //            += -1.00 f_aa(k,c) t2_abab(c,b,i,j) t1_1_aa(a,k)
    //            += -0.50 <k,b||c,d>_abab t2_abab(c,d,i,j) t1_1_aa(a,k)
    //            += -0.50 <k,b||d,c>_abab t2_abab(d,c,i,j) t1_1_aa(a,k)
    //            += +2.00 d-_aa(k,c) t1_1_aa(a,k) t2_1_abab(c,b,i,j)
    //            += -0.50 <l,k||c,d>_abab t2_abab(c,b,l,k) t2_1_abab(a,d,i,j)
    //            += -0.50 <k,l||c,d>_abab t2_abab(c,b,k,l) t2_1_abab(a,d,i,j)
    //            += +0.50 <l,k||c,d>_bbbb t2_abab(a,c,i,j) t2_1_bbbb(d,b,l,k)
    //            += -0.50 <l,k||d,c>_abab t2_abab(a,c,i,j) t2_1_abab(d,b,l,k)
    //            += -0.50 <k,l||d,c>_abab t2_abab(a,c,i,j) t2_1_abab(d,b,k,l)
    //            += -1.00 <k,b||c,d>_abab t2_aaaa(c,a,i,k) t1_1_bb(d,j)
    //            += -0.50 <l,k||c,d>_abab t2_abab(c,b,i,j) t2_1_abab(a,d,l,k)
    //            += -0.50 <k,l||c,d>_abab t2_abab(c,b,i,j) t2_1_abab(a,d,k,l)
    //            += +0.50 <l,k||c,d>_aaaa t2_abab(c,b,i,j) t2_1_aaaa(d,a,l,k)
    //            += -0.50 <l,k||c,d>_bbbb t2_bbbb(c,b,l,k) t2_1_abab(a,d,i,j)
    //            += +1.00 <k,b||c,d>_bbbb t2_abab(a,c,i,k) t1_1_bb(d,j)
    //            += -1.00 d-_aa(k,c) t2_1_abab(a,b,i,j) t1_1_aa(c,k)
    //            += -1.00 d-_bb(k,c) t2_1_abab(a,b,i,j) t1_1_bb(c,k)
    //            += -1.00 <a,k||c,j>_abab t2_1_abab(c,b,i,k)
    //            += +1.00 <a,b||c,j>_abab t1_1_aa(c,i)
    //            += +0.50 <l,k||i,j>_abab t2_1_abab(a,b,l,k)
    //            += +0.50 <k,l||i,j>_abab t2_1_abab(a,b,k,l)
    //            += +1.00 <k,a||c,i>_aaaa t2_1_abab(c,b,k,j)
    //            += -1.00 <a,k||i,c>_abab t2_1_bbbb(c,b,j,k)
    //            += -1.00 f_aa(k,i) t2_1_abab(a,b,k,j)
    //            += +1.00 t2_1_abab(a,b,i,j) w0
    //            += +1.00 <a,b||i,c>_abab t1_1_bb(c,j)
    //            += +0.50 <a,b||c,d>_abab t2_1_abab(c,d,i,j)
    //            += +0.50 <a,b||d,c>_abab t2_1_abab(d,c,i,j)
    //            += -1.00 f_bb(k,j) t2_1_abab(a,b,i,k)
    //            += -1.00 <a,k||i,j>_abab t1_1_bb(b,k)
    //            += -1.00 <k,b||i,j>_abab t1_1_aa(a,k)
    //            += +1.00 f_aa(a,c) t2_1_abab(c,b,i,j)
    //            += -1.00 <k,b||c,j>_abab t2_1_aaaa(c,a,i,k)
    //            += -1.00 <k,b||i,c>_abab t2_1_abab(a,c,k,j)
    //            += -1.00 <a,k||c,d>_abab t2_abab(c,b,i,k) t1_1_bb(d,j)
    //            += +1.00 <k,b||c,j>_bbbb t2_1_abab(a,c,i,k)
    //            += +1.00 f_bb(b,c) t2_1_abab(a,c,i,j)
    //            += +1.00 d-_bb(k,c) t2_abab(a,b,i,k) t0_1 t1_1_bb(c,j)
    //            += +1.00 d-_aa(k,c) t2_abab(c,b,i,j) t0_1 t1_1_aa(a,k)
    //            += -1.00 d-_bb(b,c) t0_1 t2_1_abab(a,c,i,j)
    //            += -1.00 d-_aa(a,c) t0_1 t2_1_abab(c,b,i,j)
    //            += +1.00 d-_aa(k,i) t0_1 t2_1_abab(a,b,k,j)
    //            += +1.00 d-_aa(k,c) t2_abab(a,b,k,j) t0_1 t1_1_aa(c,i)
    //            += +1.00 d-_bb(k,c) t2_abab(a,c,i,j) t0_1 t1_1_bb(b,k)
    //            += +1.00 d-_bb(k,j) t0_1 t2_1_abab(a,b,i,k)
    //            += +1.00 <l,k||c,d>_aaaa t2_abab(c,b,k,j) t2_1_aaaa(d,a,i,l)
    //            += -1.00 d-_aa(a,c) t1_1_bb(b,j) t1_1_aa(c,i)
    //            += +1.00 d-_aa(k,i) t1_1_aa(a,k) t1_1_bb(b,j)
    //            += +1.00 <l,k||c,j>_bbbb t2_abab(a,c,i,k) t1_1_bb(b,l)
    //            += +1.00 <k,l||i,c>_abab t2_abab(a,c,k,j) t1_1_bb(b,l)
    //            += +1.00 <k,l||c,j>_abab t2_aaaa(c,a,i,k) t1_1_bb(b,l)
    //            += -0.50 <l,k||c,d>_abab t2_abab(a,b,i,k) t2_1_abab(c,d,l,j)
    //            += -0.50 <l,k||d,c>_abab t2_abab(a,b,i,k) t2_1_abab(d,c,l,j)
    //            += -1.00 f_bb(k,c) t2_abab(a,b,i,k) t1_1_bb(c,j)
    //            += +0.50 <l,k||c,d>_bbbb t2_abab(a,b,i,k) t2_1_bbbb(c,d,j,l)
    //            += +1.00 d-_aa(k,c) t2_1_aaaa(c,a,i,k) t1_1_bb(b,j)
    //            += -1.00 d-_bb(k,c) t2_1_abab(a,c,i,k) t1_1_bb(b,j)
    //            += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,a,l,k) t2_1_abab(d,b,i,j)
    //            += -1.00 <a,k||d,c>_abab t2_bbbb(c,b,j,k) t1_1_aa(d,i)
    //            += +1.00 <k,a||c,d>_aaaa t2_abab(c,b,k,j) t1_1_aa(d,i)
    rt2_1_abab("a,b,i,j") += tmps_["109_abba_vovo"]("a,j,b,i");
    tmps_["109_abba_vovo"].~TArrayD();

    // flops: o3v1  = o4v2 o4v2 o3v1
    //  mems: o3v1  = o3v1 o3v1 o3v1
    tmps_["110_aaaa_oovo"]("i,j,a,k")  = eri["aaaa_oovo"]("m,i,c,j") * t2["aaaa"]("c,a,k,m");
    tmps_["110_aaaa_oovo"]("i,j,a,k") += eri["abba_oovo"]("i,l,b,j") * t2["abab"]("a,b,k,l");

    // flops: o2v2  = o3v3 o2v2 o3v3 o3v3 o2v2 o3v3 o3v2 o2v2 o3v3 o2v2 o3v3 o2v2 o2v2 o2v2 o3v3 o2v2 o2v2 o2v2 o3v3 o2v2 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["111_aaaa_vovo"]("a,i,b,j")  = t2["aaaa"]("c,a,i,k") * tmps_["107_aaaa_vovo"]("b,k,c,j");
    tmps_["111_aaaa_vovo"]("a,i,b,j") -= tmps_["53_aa_vo"]("a,j") * t1_1["aa"]("b,i");
    tmps_["111_aaaa_vovo"]("a,i,b,j") -= eri["aaaa_vovo"]("b,k,c,i") * t2_1["aaaa"]("c,a,j,k");
    tmps_["111_aaaa_vovo"]("a,i,b,j") += eri["abba_vovo"]("b,l,d,i") * t2_1["abab"]("a,d,j,l");
    tmps_["111_aaaa_vovo"]("a,i,b,j") -= t2_1["aaaa"]("e,a,j,m") * tmps_["9_aaaa_voov"]("b,i,m,e");
    tmps_["111_aaaa_vovo"]("a,i,b,j") += t1_1["aa"]("a,m") * tmps_["110_aaaa_oovo"]("m,i,b,j");
    tmps_["111_aaaa_vovo"]("a,i,b,j") += t2_1["abab"]("a,f,j,n") * tmps_["2_aabb_voov"]("b,i,n,f");
    tmps_["111_aaaa_vovo"]("a,i,b,j") -= t2_1["abab"]("a,f,j,n") * tmps_["10_aabb_voov"]("b,i,n,f");
    tmps_["111_aaaa_vovo"]("a,i,b,j") += t1_1["aa"]("a,j") * tmps_["54_aa_vo"]("b,i");
    tmps_["111_aaaa_vovo"]("a,i,b,j") += tmps_["106_aabb_voov"]("a,j,l,d") * t2["abab"]("b,d,i,l");
    tmps_["111_aaaa_vovo"]("a,i,b,j") += t2["abab"]("a,d,i,l") * tmps_["108_abba_vovo"]("b,l,d,j");
    tmps_["110_aaaa_oovo"].~TArrayD();
    tmps_["108_abba_vovo"].~TArrayD();
    tmps_["107_aaaa_vovo"].~TArrayD();
    tmps_["106_aabb_voov"].~TArrayD();
    tmps_["54_aa_vo"].~TArrayD();
    tmps_["53_aa_vo"].~TArrayD();
    tmps_["10_aabb_voov"].~TArrayD();
    tmps_["9_aaaa_voov"].~TArrayD();
    tmps_["2_aabb_voov"].~TArrayD();

    // rt2_1_aaaa += +1.00 P(i,j) P(a,b) <a,k||d,c>_abab t2_abab(b,c,j,k) t1_1_aa(d,i)
    //            += -1.00 P(i,j) P(a,b) <l,k||j,c>_abab t2_abab(a,c,i,k) t1_1_aa(b,l)
    //            += -1.00 P(i,j) P(a,b) <l,k||c,j>_aaaa t2_aaaa(c,a,i,k) t1_1_aa(b,l)
    //            += +1.00 P(i,j) P(a,b) <l,k||c,d>_aaaa t2_aaaa(c,a,j,k) t2_1_aaaa(d,b,i,l)
    //            += +1.00 P(i,j) P(a,b) <k,l||c,d>_abab t2_aaaa(c,a,j,k) t2_1_abab(b,d,i,l)
    //            += +1.00 P(i,j) P(a,b) <l,k||c,d>_bbbb t2_abab(a,c,j,k) t2_1_abab(b,d,i,l)
    //            += +1.00 P(i,j) P(a,b) d-_aa(a,c) t1_1_aa(b,i) t1_1_aa(c,j)
    //            += -1.00 P(i,j) P(a,b) d-_aa(k,j) t1_1_aa(a,k) t1_1_aa(b,i)
    //            += +1.00 P(i,j) P(a,b) <l,k||d,c>_abab t2_abab(a,c,j,k) t2_1_aaaa(d,b,i,l)
    //            += +1.00 P(i,j) P(a,b) <k,a||c,j>_aaaa t2_1_aaaa(c,b,i,k)
    //            += -1.00 P(i,j) P(a,b) <a,k||j,c>_abab t2_1_abab(b,c,i,k)
    //            += -1.00 P(i,j) P(a,b) d-_aa(k,c) t1_1_aa(a,j) t2_1_aaaa(c,b,i,k)
    //            += +1.00 P(i,j) P(a,b) d-_bb(k,c) t1_1_aa(a,j) t2_1_abab(b,c,i,k)
    //            += -1.00 P(i,j) P(a,b) <k,a||c,d>_aaaa t2_aaaa(c,b,j,k) t1_1_aa(d,i)
    rt2_1_aaaa("a,b,i,j") += tmps_["111_aaaa_vovo"]("b,j,a,i");
    rt2_1_aaaa("a,b,i,j") -= tmps_["111_aaaa_vovo"]("b,i,a,j");
    rt2_1_aaaa("a,b,i,j") -= tmps_["111_aaaa_vovo"]("a,j,b,i");
    rt2_1_aaaa("a,b,i,j") += tmps_["111_aaaa_vovo"]("a,i,b,j");
    tmps_["111_aaaa_vovo"].~TArrayD();
}
