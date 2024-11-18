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

void QED_CCSD::resid_21_1() {
    
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

    ///// Scalars /////

    scalars_["1"]  = -1.00 * dot(eri["abab_oovv"]("j,k,a,d"), t2_1["abab"]("a,d,j,k"));
    scalars_["1"] += 1.00 * dot(Id["aa_oo"]("i,n"), dp["aa_oo"]("i,n"));
    scalars_["1"] += 1.00 * dot(dp["bb_ov"]("k,c"), t1_1["bb"]("c,k")) * t0_1;
    scalars_["1"] += 0.25 * dot(eri["bbbb_oovv"]("k,m,c,d"), t2_1["bbbb"]("c,d,m,k"));
    scalars_["1"] -= 1.00 * dot(f["aa_ov"]("i,a"), t1_1["aa"]("a,i"));
    scalars_["1"] += 1.00 * dot(Id["bb_oo"]("k,l"), dp["bb_oo"]("k,l"));
    scalars_["1"] -= t0_1 * w0;
    scalars_["1"] += 1.00 * dot(dp["aa_ov"]("i,a"), t1_1["aa"]("a,i")) * t0_1;
    scalars_["1"] -= 1.00 * dot(f["bb_ov"]("k,c"), t1_1["bb"]("c,k"));
    scalars_["1"] += 0.25 * dot(eri["aaaa_oovv"]("i,j,a,b"), t2_1["aaaa"]("a,b,j,i"));
    scalars_["2"]  = -1.00 * dot(Id["aa_oo"]("i,j"), f["aa_oo"]("i,j"));
    scalars_["2"] -= 1.00 * dot(Id["bb_oo"]("m,I"), f["bb_oo"]("m,I"));
    scalars_["2"] += 1.00 * dot(Id["bb_oo"]("m,I") * eri["abab_oooo"]("n,m,o,I"), Id["aa_oo"]("n,o"));
    scalars_["2"] += 0.25 * dot(eri["bbbb_oovv"]("m,k,a,d"), t2["bbbb"]("a,d,k,m"));
    scalars_["2"] += 1.00 * dot(dp["aa_ov"]("i,b"), t1_1["aa"]("b,i"));
    scalars_["2"] -= 1.00 * dot(eri["abab_oovv"]("n,m,b,d"), t2["abab"]("b,d,n,m"));
    scalars_["2"] += 0.25 * dot(eri["aaaa_oovv"]("i,n,b,c"), t2["aaaa"]("b,c,n,i"));
    scalars_["2"] += 1.00 * dot(Id["bb_oo"]("m,I"), dp["bb_oo"]("m,I")) * t0_1;
    scalars_["2"] += 1.00 * dot(dp["bb_ov"]("m,a"), t1_1["bb"]("a,m"));
    scalars_["2"] += 0.50 * dot(Id["aa_oo"]("i,j") * eri["aaaa_oooo"]("i,n,j,o"), Id["aa_oo"]("n,o"));
    scalars_["2"] += 0.50 * dot(Id["bb_oo"]("m,I") * eri["bbbb_oooo"]("m,k,I,l"), Id["bb_oo"]("k,l"));
    scalars_["2"] += 1.00 * dot(Id["aa_oo"]("i,j"), dp["aa_oo"]("i,j")) * t0_1;
    scalars_["3"]  = 1.00 * dot(Id["aa_oo"]("k,l"), dp["aa_oo"]("k,l"));
    scalars_["3"] += 1.00 * dot(Id["bb_oo"]("i,j"), dp["bb_oo"]("i,j"));
    scalars_["4"]  = 1.00 * dot(dp["aa_ov"]("i,a"), t1_1["aa"]("a,i"));
    scalars_["4"] += 1.00 * dot(dp["bb_ov"]("j,b"), t1_1["bb"]("b,j"));

///// End of Scalars /////

/////////////////// Evaluate Equations ///////////////////


    // crt0_1  = +1.00
    crt0_1  = 1.00;

    // crt1_aa  = +1.00 t1_1_aa(a,i)
    crt1_aa("a,i")  = t1_1["aa"]("a,i");

    // crt1_bb  = +1.00 t1_1_bb(a,i)
    crt1_bb("a,i")  = t1_1["bb"]("a,i");

    // crt2_aaaa  = +1.00 t2_1_aaaa(a,b,i,j)
    crt2_aaaa("a,b,i,j")  = t2_1["aaaa"]("a,b,i,j");

    // crt2_abab  = +1.00 t2_1_abab(a,b,i,j)
    crt2_abab("a,b,i,j")  = t2_1["abab"]("a,b,i,j");

    // crt2_bbbb  = +1.00 t2_1_bbbb(a,b,i,j)
    crt2_bbbb("a,b,i,j")  = t2_1["bbbb"]("a,b,i,j");

    // energy  = -0.50 <j,i||j,i>_abab
    //        += -0.50 <i,j||i,j>_abab
    //        += +1.00 f_bb(i,i)
    //        += +0.25 <j,i||a,b>_bbbb t2_bbbb(a,b,j,i)
    //        += -1.00 d-_aa(i,a) t1_1_aa(a,i)
    //        += +0.25 <j,i||a,b>_abab t2_abab(a,b,j,i)
    //        += +0.25 <i,j||a,b>_abab t2_abab(a,b,i,j)
    //        += +0.25 <j,i||b,a>_abab t2_abab(b,a,j,i)
    //        += +0.25 <i,j||b,a>_abab t2_abab(b,a,i,j)
    //        += +0.25 <j,i||a,b>_aaaa t2_aaaa(a,b,j,i)
    //        += -1.00 d-_bb(i,i) t0_1
    //        += -1.00 d-_bb(i,a) t1_1_bb(a,i)
    //        += -0.50 <j,i||j,i>_aaaa
    //        += -0.50 <j,i||j,i>_bbbb
    //        += +1.00 f_aa(i,i)
    //        += -1.00 d-_aa(i,i) t0_1
    energy  = -1.00 * scalars_["2"];

    // rt0_1  = -1.00 d+_aa(i,i)
    //       += +0.25 <j,i||a,b>_abab t2_1_abab(a,b,j,i)
    //       += +0.25 <i,j||a,b>_abab t2_1_abab(a,b,i,j)
    //       += +0.25 <j,i||b,a>_abab t2_1_abab(b,a,j,i)
    //       += +0.25 <i,j||b,a>_abab t2_1_abab(b,a,i,j)
    //       += -1.00 d-_bb(i,a) t0_1 t1_1_bb(a,i)
    //       += +0.25 <j,i||a,b>_bbbb t2_1_bbbb(a,b,j,i)
    //       += +1.00 f_aa(i,a) t1_1_aa(a,i)
    //       += -1.00 d+_bb(i,i)
    //       += +1.00 t0_1 w0
    //       += -1.00 d-_aa(i,a) t0_1 t1_1_aa(a,i)
    //       += +1.00 f_bb(i,a) t1_1_bb(a,i)
    //       += +0.25 <j,i||a,b>_aaaa t2_1_aaaa(a,b,j,i)
    rt0_1  = -1.00 * scalars_["1"];

    // cenergy  = +1.00 t0_1
    cenergy  = t0_1;

    // flops: o2v2  = o2v3 o2v3
    //  mems: o2v2  = o0v2 o2v2
    tmps_["1_aaaa_vvoo"]("a,b,i,j")  = eri["abab_oovv"]("k,l,c,d") * t2["abab"]("a,d,k,l") * t2["aaaa"]("c,b,i,j");

    // rt2_aaaa  = +0.50 <l,k||c,d>_abab t2_aaaa(c,a,i,j) t2_abab(b,d,l,k)
    //          += +0.50 <k,l||c,d>_abab t2_aaaa(c,a,i,j) t2_abab(b,d,k,l)
    rt2_aaaa("a,b,i,j")  = tmps_["1_aaaa_vvoo"]("b,a,i,j");

    // rt2_aaaa += -0.50 <l,k||d,c>_abab t2_abab(a,c,l,k) t2_aaaa(d,b,i,j)
    //          += -0.50 <k,l||d,c>_abab t2_abab(a,c,k,l) t2_aaaa(d,b,i,j)
    rt2_aaaa("a,b,i,j") -= tmps_["1_aaaa_vvoo"]("a,b,i,j");
    tmps_["1_aaaa_vvoo"].~TArrayD();

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["2_aabb_voov"]("a,i,j,b")  = t2["aaaa"]("c,a,i,k") * eri["abab_oovv"]("k,j,c,b");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["3_aaaa_vovo"]("a,i,b,j")  = t2["abab"]("a,c,i,k") * tmps_["2_aabb_voov"]("b,j,k,c");

    // rt2_aaaa += +1.00 P(i,j) <l,k||d,c>_abab t2_abab(a,c,j,k) t2_aaaa(d,b,i,l)
    rt2_aaaa("a,b,i,j") += tmps_["3_aaaa_vovo"]("a,j,b,i");
    rt2_aaaa("a,b,i,j") -= tmps_["3_aaaa_vovo"]("a,i,b,j");

    // rt2_aaaa += +1.00 P(i,j) <k,l||c,d>_abab t2_aaaa(c,a,j,k) t2_abab(b,d,i,l)
    rt2_aaaa("a,b,i,j") += tmps_["3_aaaa_vovo"]("b,i,a,j");
    rt2_aaaa("a,b,i,j") -= tmps_["3_aaaa_vovo"]("b,j,a,i");
    tmps_["3_aaaa_vovo"].~TArrayD();

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["4_bbbb_voov"]("a,i,j,b")  = t2["abab"]("c,a,k,i") * eri["abab_oovv"]("k,j,c,b");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["5_bbbb_vovo"]("a,i,b,j")  = tmps_["4_bbbb_voov"]("a,i,k,c") * t2["bbbb"]("c,b,j,k");

    // rt2_bbbb  = +1.00 P(i,j) <l,k||d,c>_abab t2_bbbb(c,a,j,k) t2_abab(d,b,l,i)
    rt2_bbbb("a,b,i,j")  = tmps_["5_bbbb_vovo"]("b,i,a,j");
    rt2_bbbb("a,b,i,j") -= tmps_["5_bbbb_vovo"]("b,j,a,i");

    // rt2_bbbb += +1.00 P(i,j) <k,l||c,d>_abab t2_abab(c,a,k,j) t2_bbbb(d,b,i,l)
    rt2_bbbb("a,b,i,j") += tmps_["5_bbbb_vovo"]("a,j,b,i");
    rt2_bbbb("a,b,i,j") -= tmps_["5_bbbb_vovo"]("a,i,b,j");
    tmps_["5_bbbb_vovo"].~TArrayD();

    // flops: o0v2  = o2v3
    //  mems: o0v2  = o0v2
    tmps_["6_bb_vv"]("a,b")  = t2["abab"]("c,a,i,j") * eri["abab_oovv"]("i,j,c,b");

    // flops: o2v2  = o2v3
    //  mems: o2v2  = o2v2
    tmps_["7_bbbb_voov"]("a,i,j,b")  = t2["bbbb"]("c,a,i,j") * tmps_["6_bb_vv"]("b,c");

    // rt2_bbbb += +0.50 <l,k||d,c>_abab t2_bbbb(c,a,i,j) t2_abab(d,b,l,k)
    //          += +0.50 <k,l||d,c>_abab t2_bbbb(c,a,i,j) t2_abab(d,b,k,l)
    rt2_bbbb("a,b,i,j") += tmps_["7_bbbb_voov"]("a,i,j,b");

    // rt2_bbbb += -0.50 <l,k||c,d>_abab t2_abab(c,a,l,k) t2_bbbb(d,b,i,j)
    //          += -0.50 <k,l||c,d>_abab t2_abab(c,a,k,l) t2_bbbb(d,b,i,j)
    rt2_bbbb("a,b,i,j") -= tmps_["7_bbbb_voov"]("b,i,j,a");
    tmps_["7_bbbb_voov"].~TArrayD();

    // flops: o2v2  = o3v2
    //  mems: o2v2  = o2v2
    tmps_["8_aaaa_ovvo"]("i,a,b,j")  = dp["aa_oo"]("k,i") * t2["aaaa"]("a,b,j,k");

    // rt2_1_aaaa  = +1.00 P(i,j) d+_aa(k,j) t2_aaaa(a,b,i,k)
    rt2_1_aaaa("a,b,i,j")  = tmps_["8_aaaa_ovvo"]("j,a,b,i");
    rt2_1_aaaa("a,b,i,j") -= tmps_["8_aaaa_ovvo"]("i,a,b,j");

    // flops: o2v0  = o3v2 o3v2 o2v0
    //  mems: o2v0  = o2v0 o2v0 o2v0
    tmps_["11_aa_oo"]("i,j")  = -2.00 * t2["abab"]("a,c,i,l") * eri["abab_oovv"]("j,l,a,c");
    tmps_["11_aa_oo"]("i,j") += t2["aaaa"]("a,b,i,k") * eri["aaaa_oovv"]("k,j,a,b");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["10_aabb_voov"]("a,i,j,b")  = t2["abab"]("a,c,i,k") * eri["bbbb_oovv"]("k,j,c,b");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["9_aaaa_voov"]("a,i,j,b")  = t2["aaaa"]("c,a,i,k") * eri["aaaa_oovv"]("k,j,c,b");

    // flops: o2v2  = o3v2 o2v2 o2v2 o3v3 o2v2 o3v3 o2v2 o3v2 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["12_aaaa_ovvo"]("i,a,b,j")  = -1.00 * t2["aaaa"]("a,b,j,m") * f["aa_oo"]("m,i");
    tmps_["12_aaaa_ovvo"]("i,a,b,j") += t0_1 * tmps_["8_aaaa_ovvo"]("i,a,b,j");
    tmps_["12_aaaa_ovvo"]("i,a,b,j") -= t2["abab"]("b,d,j,l") * tmps_["10_aabb_voov"]("a,i,l,d");
    tmps_["12_aaaa_ovvo"]("i,a,b,j") -= t2["aaaa"]("c,b,j,k") * tmps_["9_aaaa_voov"]("a,i,k,c");
    tmps_["12_aaaa_ovvo"]("i,a,b,j") += 0.50 * t2["aaaa"]("a,b,j,k") * tmps_["11_aa_oo"]("i,k");
    tmps_["8_aaaa_ovvo"].~TArrayD();

    // rt2_aaaa += +1.00 P(i,j) d-_aa(k,j) t2_aaaa(a,b,i,k) t0_1
    //          += -1.00 P(i,j) f_aa(k,j) t2_aaaa(a,b,i,k)
    //          += +1.00 P(i,j) <l,k||c,d>_bbbb t2_abab(a,c,j,k) t2_abab(b,d,i,l)
    //          += +1.00 P(i,j) <l,k||c,d>_aaaa t2_aaaa(c,a,j,k) t2_aaaa(d,b,i,l)
    //          += -0.50 P(i,j) <l,k||c,d>_aaaa t2_aaaa(a,b,i,l) t2_aaaa(c,d,j,k)
    //          += -0.50 P(i,j) <l,k||c,d>_abab t2_aaaa(a,b,i,l) t2_abab(c,d,j,k)
    //          += -0.50 P(i,j) <l,k||d,c>_abab t2_aaaa(a,b,i,l) t2_abab(d,c,j,k)
    rt2_aaaa("a,b,i,j") += tmps_["12_aaaa_ovvo"]("j,a,b,i");
    rt2_aaaa("a,b,i,j") -= tmps_["12_aaaa_ovvo"]("i,a,b,j");
    tmps_["12_aaaa_ovvo"].~TArrayD();

    // flops: o2v2  = o3v2
    //  mems: o2v2  = o2v2
    tmps_["13_bbbb_ovvo"]("i,a,b,j")  = dp["bb_oo"]("k,i") * t2["bbbb"]("a,b,j,k");

    // rt2_1_bbbb  = +1.00 P(i,j) d+_bb(k,j) t2_bbbb(a,b,i,k)
    rt2_1_bbbb("a,b,i,j")  = tmps_["13_bbbb_ovvo"]("j,a,b,i");
    rt2_1_bbbb("a,b,i,j") -= tmps_["13_bbbb_ovvo"]("i,a,b,j");

    // flops: o2v0  = o3v2 o3v2 o2v0
    //  mems: o2v0  = o2v0 o2v0 o2v0
    tmps_["16_bb_oo"]("i,j")  = -2.00 * t2["abab"]("c,b,l,i") * eri["abab_oovv"]("l,j,c,b");
    tmps_["16_bb_oo"]("i,j") += t2["bbbb"]("a,b,i,k") * eri["bbbb_oovv"]("k,j,a,b");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["15_bbbb_voov"]("a,i,j,b")  = t2["bbbb"]("c,a,i,k") * eri["bbbb_oovv"]("k,j,c,b");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["14_bbaa_voov"]("a,i,j,b")  = t2["abab"]("c,a,k,i") * eri["aaaa_oovv"]("k,j,c,b");

    // flops: o2v2  = o3v2 o2v2 o3v2 o2v2 o2v2 o3v3 o2v2 o3v3 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["17_bbbb_ovvo"]("i,a,b,j")  = -0.50 * t2["bbbb"]("a,b,j,k") * tmps_["16_bb_oo"]("i,k");
    tmps_["17_bbbb_ovvo"]("i,a,b,j") -= t0_1 * tmps_["13_bbbb_ovvo"]("i,a,b,j");
    tmps_["17_bbbb_ovvo"]("i,a,b,j") += f["bb_oo"]("m,i") * t2["bbbb"]("a,b,j,m");
    tmps_["17_bbbb_ovvo"]("i,a,b,j") += t2["abab"]("d,b,l,j") * tmps_["14_bbaa_voov"]("a,i,l,d");
    tmps_["17_bbbb_ovvo"]("i,a,b,j") += t2["bbbb"]("c,b,j,k") * tmps_["15_bbbb_voov"]("a,i,k,c");
    tmps_["13_bbbb_ovvo"].~TArrayD();

    // rt2_bbbb += -1.00 P(i,j) f_bb(k,j) t2_bbbb(a,b,i,k)
    //          += +1.00 P(i,j) d-_bb(k,j) t2_bbbb(a,b,i,k) t0_1
    //          += -0.50 P(i,j) <l,k||c,d>_bbbb t2_bbbb(a,b,i,l) t2_bbbb(c,d,j,k)
    //          += -0.50 P(i,j) <k,l||c,d>_abab t2_bbbb(a,b,i,l) t2_abab(c,d,k,j)
    //          += -0.50 P(i,j) <k,l||d,c>_abab t2_bbbb(a,b,i,l) t2_abab(d,c,k,j)
    //          += +1.00 P(i,j) <l,k||c,d>_aaaa t2_abab(c,a,k,j) t2_abab(d,b,l,i)
    //          += +1.00 P(i,j) <l,k||c,d>_bbbb t2_bbbb(c,a,j,k) t2_bbbb(d,b,i,l)
    rt2_bbbb("a,b,i,j") -= tmps_["17_bbbb_ovvo"]("j,a,b,i");
    rt2_bbbb("a,b,i,j") += tmps_["17_bbbb_ovvo"]("i,a,b,j");
    tmps_["17_bbbb_ovvo"].~TArrayD();

    // flops: o2v2  = o2v3
    //  mems: o2v2  = o2v2
    tmps_["18_bbbb_vvoo"]("a,b,i,j")  = dp["bb_vv"]("a,c") * t2["bbbb"]("c,b,i,j");

    // rt2_1_bbbb += -1.00 P(a,b) d+_bb(a,c) t2_bbbb(c,b,i,j)
    rt2_1_bbbb("a,b,i,j") -= tmps_["18_bbbb_vvoo"]("a,b,i,j");
    rt2_1_bbbb("a,b,i,j") += tmps_["18_bbbb_vvoo"]("b,a,i,j");

    // flops: o2v2  = o2v2 o2v3 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2
    tmps_["19_bbbb_vvoo"]("a,b,i,j")  = -1.00 * t0_1 * tmps_["18_bbbb_vvoo"]("a,b,i,j");
    tmps_["19_bbbb_vvoo"]("a,b,i,j") += f["bb_vv"]("a,c") * t2["bbbb"]("c,b,i,j");

    // rt2_bbbb += +1.00 P(a,b) f_bb(a,c) t2_bbbb(c,b,i,j)
    //          += -1.00 P(a,b) d-_bb(a,c) t2_bbbb(c,b,i,j) t0_1
    rt2_bbbb("a,b,i,j") += tmps_["19_bbbb_vvoo"]("a,b,i,j");
    rt2_bbbb("a,b,i,j") -= tmps_["19_bbbb_vvoo"]("b,a,i,j");
    tmps_["19_bbbb_vvoo"].~TArrayD();

    // flops: o2v2  = o2v2
    //  mems: o2v2  = o2v2
    tmps_["20_bbbb_vvoo"]("a,b,i,j")  = -1.00 * t0_1 * tmps_["18_bbbb_vvoo"]("a,b,i,j");
    tmps_["20_bbbb_vvoo"].~TArrayD();
    tmps_["18_bbbb_vvoo"].~TArrayD();

    // flops: o2v2  = o2v3
    //  mems: o2v2  = o2v2
    tmps_["21_aaaa_vvoo"]("a,b,i,j")  = dp["aa_vv"]("a,c") * t2["aaaa"]("c,b,i,j");

    // rt2_1_aaaa += -1.00 P(a,b) d+_aa(a,c) t2_aaaa(c,b,i,j)
    rt2_1_aaaa("a,b,i,j") -= tmps_["21_aaaa_vvoo"]("a,b,i,j");
    rt2_1_aaaa("a,b,i,j") += tmps_["21_aaaa_vvoo"]("b,a,i,j");

    // flops: o2v2  = o2v2 o2v3 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2
    tmps_["22_aaaa_vvoo"]("a,b,i,j")  = -1.00 * t0_1 * tmps_["21_aaaa_vvoo"]("a,b,i,j");
    tmps_["22_aaaa_vvoo"]("a,b,i,j") += f["aa_vv"]("a,c") * t2["aaaa"]("c,b,i,j");

    // rt2_aaaa += +1.00 P(a,b) f_aa(a,c) t2_aaaa(c,b,i,j)
    //          += -1.00 P(a,b) d-_aa(a,c) t2_aaaa(c,b,i,j) t0_1
    rt2_aaaa("a,b,i,j") += tmps_["22_aaaa_vvoo"]("a,b,i,j");
    rt2_aaaa("a,b,i,j") -= tmps_["22_aaaa_vvoo"]("b,a,i,j");
    tmps_["22_aaaa_vvoo"].~TArrayD();

    // flops: o2v2  = o2v2
    //  mems: o2v2  = o2v2
    tmps_["23_aaaa_vvoo"]("a,b,i,j")  = -1.00 * t0_1 * tmps_["21_aaaa_vvoo"]("a,b,i,j");
    tmps_["23_aaaa_vvoo"].~TArrayD();
    tmps_["21_aaaa_vvoo"].~TArrayD();

    // flops: o1v1  = o2v2
    //  mems: o1v1  = o1v1
    tmps_["24_aa_vo"]("a,i")  = -1.00 * dp["bb_ov"]("j,b") * t2["abab"]("a,b,i,j");
    tmps_["24_aa_vo"].~TArrayD();

    // flops: o1v1  = o2v2 o2v2 o1v1
    //  mems: o1v1  = o1v1 o1v1 o1v1
    tmps_["25_aa_vo"]("a,i")  = -1.00 * dp["bb_ov"]("k,c") * t2["abab"]("a,c,i,k");
    tmps_["25_aa_vo"]("a,i") += dp["aa_ov"]("j,b") * t2["aaaa"]("b,a,i,j");

    // rt1_1_aa  = +1.00 d+_aa(j,b) t2_aaaa(b,a,i,j)
    //          += -1.00 d+_bb(j,b) t2_abab(a,b,i,j)
    rt1_1_aa("a,i")  = tmps_["25_aa_vo"]("a,i");

    // flops: o2v2  = o3v3 o2v2 o3v3 o2v2 o2v2 o2v2 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["26_aaaa_vovo"]("a,i,b,j")  = -1.00 * eri["aaaa_vovo"]("a,l,d,i") * t2["aaaa"]("d,b,j,l");
    tmps_["26_aaaa_vovo"]("a,i,b,j") += dp["aa_vo"]("a,i") * t1_1["aa"]("b,j");
    tmps_["26_aaaa_vovo"]("a,i,b,j") += eri["abba_vovo"]("a,k,c,i") * t2["abab"]("b,c,j,k");
    tmps_["26_aaaa_vovo"]("a,i,b,j") -= t1_1["aa"]("b,j") * tmps_["25_aa_vo"]("a,i");

    // rt2_aaaa += -1.00 P(i,j) P(a,b) <a,k||j,c>_abab t2_abab(b,c,i,k)
    //          += +1.00 P(i,j) P(a,b) d-_aa(a,j) t1_1_aa(b,i)
    //          += +1.00 P(i,j) P(a,b) <k,a||c,j>_aaaa t2_aaaa(c,b,i,k)
    //          += -1.00 P(i,j) P(a,b) d-_aa(k,c) t2_aaaa(c,a,j,k) t1_1_aa(b,i)
    //          += +1.00 P(i,j) P(a,b) d-_bb(k,c) t2_abab(a,c,j,k) t1_1_aa(b,i)
    rt2_aaaa("a,b,i,j") += tmps_["26_aaaa_vovo"]("a,j,b,i");
    rt2_aaaa("a,b,i,j") -= tmps_["26_aaaa_vovo"]("a,i,b,j");
    rt2_aaaa("a,b,i,j") -= tmps_["26_aaaa_vovo"]("b,j,a,i");
    rt2_aaaa("a,b,i,j") += tmps_["26_aaaa_vovo"]("b,i,a,j");
    tmps_["26_aaaa_vovo"].~TArrayD();

    // flops: o1v1  = o2v2 o2v2 o1v1
    //  mems: o1v1  = o1v1 o1v1 o1v1
    tmps_["29_bb_vo"]("a,i")  = -1.00 * dp["bb_ov"]("k,c") * t2["bbbb"]("c,a,i,k");
    tmps_["29_bb_vo"]("a,i") += dp["aa_ov"]("j,b") * t2["abab"]("b,a,j,i");

    // rt1_1_bb  = -1.00 d+_aa(j,b) t2_abab(b,a,j,i)
    //          += +1.00 d+_bb(j,b) t2_bbbb(b,a,i,j)
    rt1_1_bb("a,i")  = -1.00 * tmps_["29_bb_vo"]("a,i");

    // flops: o1v1  = o2v2
    //  mems: o1v1  = o1v1
    tmps_["27_bb_vo"]("a,i")  = -1.00 * dp["bb_ov"]("j,b") * t2["bbbb"]("b,a,i,j");
    tmps_["27_bb_vo"].~TArrayD();

    // flops: o1v1  = o2v2
    //  mems: o1v1  = o1v1
    tmps_["28_bb_vo"]("a,i")  = -1.00 * f["bb_ov"]("j,b") * t2["bbbb"]("b,a,i,j");

    // flops: o1v1  = o2v3 o1v1 o2v3 o3v2 o1v1 o2v2 o1v1 o1v1 o1v1 o3v2 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1
    //  mems: o1v1  = o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1
    tmps_["30_bb_vo"]("a,i")  = -1.00 * eri["baab_vovv"]("a,k,c,d") * t2["abab"]("c,d,k,i");
    tmps_["30_bb_vo"]("a,i") -= scalars_["3"] * t1_1["bb"]("a,i");
    tmps_["30_bb_vo"]("a,i") += 0.50 * eri["bbbb_vovv"]("a,j,b,d") * t2["bbbb"]("b,d,i,j");
    tmps_["30_bb_vo"]("a,i") += 0.50 * t2["bbbb"]("b,a,o,j") * eri["bbbb_oovo"]("j,o,b,i");
    tmps_["30_bb_vo"]("a,i") += f["bb_vo"]("a,i");
    tmps_["30_bb_vo"]("a,i") += f["aa_ov"]("k,c") * t2["abab"]("c,a,k,i");
    tmps_["30_bb_vo"]("a,i") -= t2["abab"]("c,a,l,j") * eri["abab_oovo"]("l,j,c,i");
    tmps_["30_bb_vo"]("a,i") -= t0_1 * dp["bb_vo"]("a,i");
    tmps_["30_bb_vo"]("a,i") += tmps_["28_bb_vo"]("a,i");
    tmps_["30_bb_vo"]("a,i") -= t0_1 * tmps_["29_bb_vo"]("a,i");
    tmps_["28_bb_vo"].~TArrayD();

    // rt1_bb  = +1.00 f_bb(a,i)
    //        += -0.50 <k,j||b,i>_bbbb t2_bbbb(b,a,k,j)
    //        += +1.00 f_aa(j,b) t2_abab(b,a,j,i)
    //        += -0.50 <j,a||b,c>_bbbb t2_bbbb(b,c,i,j)
    //        += -1.00 d-_bb(j,j) t1_1_bb(a,i)
    //        += -1.00 d-_aa(j,j) t1_1_bb(a,i)
    //        += -0.50 <k,j||b,i>_abab t2_abab(b,a,k,j)
    //        += -0.50 <j,k||b,i>_abab t2_abab(b,a,j,k)
    //        += +0.50 <j,a||b,c>_abab t2_abab(b,c,j,i)
    //        += +0.50 <j,a||c,b>_abab t2_abab(c,b,j,i)
    //        += -1.00 d-_bb(a,i) t0_1
    //        += -1.00 f_bb(j,b) t2_bbbb(b,a,i,j)
    //        += -1.00 d-_aa(j,b) t2_abab(b,a,j,i) t0_1
    //        += +1.00 d-_bb(j,b) t2_bbbb(b,a,i,j) t0_1
    rt1_bb("a,i")  = tmps_["30_bb_vo"]("a,i");
    tmps_["30_bb_vo"].~TArrayD();
}
