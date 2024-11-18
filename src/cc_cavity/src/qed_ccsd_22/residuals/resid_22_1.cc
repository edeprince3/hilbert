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

#include "cc_cavity/include/qed_ccsd_22/qed_ccsd_22.h"


using namespace std;
using namespace TA;
using namespace hilbert;

#if MAX_PHOTON_LEVEL >= 2
void QED_CCSD_22::resid_22_1() {

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
    TArrayD &rt1_2_aa   = residuals_["t1_2_aa"];
    TArrayD &rt1_2_bb   = residuals_["t1_2_bb"];
    TArrayD &rt2_aaaa   = residuals_["t2_aaaa"];
    TArrayD &rt2_abab   = residuals_["t2_abab"];
    TArrayD &rt2_bbbb   = residuals_["t2_bbbb"];
    TArrayD &rt2_1_aaaa = residuals_["t2_1_aaaa"];
    TArrayD &rt2_1_abab = residuals_["t2_1_abab"];
    TArrayD &rt2_1_bbbb = residuals_["t2_1_bbbb"];
    TArrayD &rt2_2_aaaa = residuals_["t2_2_aaaa"];
    TArrayD &rt2_2_abab = residuals_["t2_2_abab"];
    TArrayD &rt2_2_bbbb = residuals_["t2_2_bbbb"];
    double  &rt0_1 = scalars_["rt0_1"];
    double  &rt0_2 = scalars_["rt0_2"];

    // extract coherents
    TArrayD &crt1_aa     = tmps_["ct1_aa"];
    TArrayD &crt1_bb     = tmps_["ct1_bb"];
    TArrayD &crt1_1_aa   = tmps_["ct1_1_aa"];
    TArrayD &crt1_1_bb   = tmps_["ct1_1_bb"];
    TArrayD &crt1_2_aa   = tmps_["ct1_2_aa"];
    TArrayD &crt1_2_bb   = tmps_["ct1_2_bb"];
    TArrayD &crt2_aaaa   = tmps_["ct2_aaaa"];
    TArrayD &crt2_abab   = tmps_["ct2_abab"];
    TArrayD &crt2_bbbb   = tmps_["ct2_bbbb"];
    TArrayD &crt2_1_aaaa = tmps_["ct2_1_aaaa"];
    TArrayD &crt2_1_abab = tmps_["ct2_1_abab"];
    TArrayD &crt2_1_bbbb = tmps_["ct2_1_bbbb"];
    TArrayD &crt2_2_aaaa = tmps_["ct2_2_aaaa"];
    TArrayD &crt2_2_abab = tmps_["ct2_2_abab"];
    TArrayD &crt2_2_bbbb = tmps_["ct2_2_bbbb"];
    double  &crt0_1 = scalars_["crt0_1"];
    double  &crt0_2 = scalars_["crt0_2"];

    // extract amplitudes
    double t0_1 = scalars_["t0_1"];
    double t0_2 = scalars_["t0_2"];

    // extract 1-body amplitudes
    std::map<std::string, TA::TArrayD> t1 {
            {"aa", amplitudes_["t1_aa"]},
            {"bb", amplitudes_["t1_bb"]}
    };
    std::map<std::string, TA::TArrayD> t1_1 {
            {"aa", amplitudes_["t1_1_aa"]},
            {"bb", amplitudes_["t1_1_bb"]}
    };
    std::map<std::string, TA::TArrayD> t1_2 {
            {"aa", amplitudes_["t1_2_aa"]},
            {"bb", amplitudes_["t1_2_bb"]}
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
    std::map<std::string, TA::TArrayD> t2_2 {
            {"aaaa", amplitudes_["t2_2_aaaa"]},
            {"abab", amplitudes_["t2_2_abab"]},
            {"bbbb", amplitudes_["t2_2_bbbb"]}
    };

    scalars_["1"]  = 1.00 * dot(dp["aa_ov"]("i,a"), t1_1["aa"]("a,i"));
    scalars_["1"] += 1.00 * dot(dp["bb_ov"]("j,b"), t1_1["bb"]("b,j"));
    scalars_["1"]  = 1.00 * dot(dp["aa_ov"]("i,a"), t1_1["aa"]("a,i"));
    scalars_["1"] += 1.00 * dot(dp["bb_ov"]("j,b"), t1_1["bb"]("b,j"));
    scalars_["1"]  = 1.00 * dot(dp["aa_ov"]("i,a"), t1_1["aa"]("a,i"));
    scalars_["1"] += 1.00 * dot(dp["bb_ov"]("j,b"), t1_1["bb"]("b,j"));
    scalars_["2"]  = 2.00 * dot(Id["aa_oo"]("i,m"), dp["aa_oo"]("i,m")) * t0_1;
    scalars_["2"] += 1.00 * dot(Id["aa_oo"]("j,I") * eri["aaaa_oooo"]("i,j,m,I"), Id["aa_oo"]("i,m"));
    scalars_["2"] += 2.00 * dot(Id["bb_oo"]("k,l") * eri["abab_oooo"]("j,k,I,l"), Id["aa_oo"]("j,I"));
    scalars_["2"] -= 2.00 * dot(Id["aa_oo"]("i,m"), f["aa_oo"]("i,m"));
    scalars_["2"] += 1.00 * dot(Id["bb_oo"]("n,o") * eri["bbbb_oooo"]("k,n,l,o"), Id["bb_oo"]("k,l"));
    scalars_["2"] += 2.00 * dot(Id["bb_oo"]("k,l"), dp["bb_oo"]("k,l")) * t0_1;
    scalars_["2"] += 0.50 * dot(eri["bbbb_oovv"]("k,n,d,c"), t2["bbbb"]("d,c,n,k"));
    scalars_["2"] -= 2.00 * dot(eri["abab_oovv"]("j,k,a,c"), t2["abab"]("a,c,j,k"));
    scalars_["2"] -= 2.00 * dot(Id["bb_oo"]("k,l"), f["bb_oo"]("k,l"));
    scalars_["2"] += 0.50 * dot(eri["aaaa_oovv"]("i,j,a,b"), t2["aaaa"]("a,b,j,i"));
    scalars_["3"]  = -1.00 * dot(Id["aa_oo"]("i,n"), dp["aa_oo"]("i,n"));
    scalars_["3"] += 1.00 * dot(eri["abab_oovv"]("m,j,a,c"), t2_1["abab"]("a,c,m,j"));
    scalars_["3"] -= 0.25 * dot(eri["aaaa_oovv"]("i,m,a,d"), t2_1["aaaa"]("a,d,m,i"));
    scalars_["3"] -= 2.00 * dot(Id["aa_oo"]("i,n"), dp["aa_oo"]("i,n")) * t0_2;
    scalars_["3"] -= 2.00 * dot(dp["bb_ov"]("j,b"), t1_2["bb"]("b,j"));
    scalars_["3"] += t0_1 * w0;
    scalars_["3"] += 1.00 * dot(f["bb_ov"]("j,b"), t1_1["bb"]("b,j"));
    scalars_["3"] -= 1.00 * dot(dp["bb_ov"]("j,b"), t1_1["bb"]("b,j")) * t0_1;
    scalars_["3"] -= 2.00 * dot(dp["aa_ov"]("i,a"), t1_2["aa"]("a,i"));
    scalars_["3"] -= 2.00 * dot(Id["bb_oo"]("j,k"), dp["bb_oo"]("j,k")) * t0_2;
    scalars_["3"] -= 0.25 * dot(eri["bbbb_oovv"]("j,l,b,c"), t2_1["bbbb"]("b,c,l,j"));
    scalars_["3"] += 1.00 * dot(f["aa_ov"]("i,a"), t1_1["aa"]("a,i"));
    scalars_["3"] -= 1.00 * dot(Id["bb_oo"]("j,k"), dp["bb_oo"]("j,k"));
    scalars_["3"] -= 1.00 * dot(dp["aa_ov"]("i,a"), t1_1["aa"]("a,i")) * t0_1;
    scalars_["4"]  = -2.00 * dot(dp["aa_ov"]("j,b"), t1_1["aa"]("b,j")) * t0_2;
    scalars_["4"] -= 1.00 * dot(dp["bb_ov"]("i,a"), t1_2["bb"]("a,i")) * t0_1;
    scalars_["4"] += 0.50 * dot(eri["aaaa_oovv"]("j,l,b,d") * t1_1["aa"]("d,l"), t1_1["aa"]("b,j"));
    scalars_["4"] += 1.00 * dot(f["bb_ov"]("i,a"), t1_2["bb"]("a,i"));
    scalars_["4"] += 1.00 * dot(eri["abab_oovv"]("j,k,b,c") * t1_1["aa"]("b,j"), t1_1["bb"]("c,k"));
    scalars_["4"] -= 0.25 * dot(eri["bbbb_oovv"]("i,k,a,c"), t2_2["bbbb"]("a,c,k,i"));
    scalars_["4"] += 2.00 * t0_2 * w0;
    scalars_["4"] -= 2.00 * dot(dp["bb_ov"]("i,a"), t1_1["bb"]("a,i")) * t0_2;
    scalars_["4"] += 1.00 * dot(eri["abab_oovv"]("l,i,b,c"), t2_2["abab"]("b,c,l,i"));
    scalars_["4"] -= 0.25 * dot(eri["aaaa_oovv"]("j,l,b,d"), t2_2["aaaa"]("b,d,l,j"));
    scalars_["4"] -= 1.00 * dot(dp["aa_ov"]("j,b"), t1_2["aa"]("b,j")) * t0_1;
    scalars_["4"] += 1.00 * dot(f["aa_ov"]("j,b"), t1_2["aa"]("b,j"));
    scalars_["4"] += 0.50 * dot(eri["bbbb_oovv"]("i,k,a,c") * t1_1["bb"]("c,k"), t1_1["bb"]("a,i"));
    scalars_["5"]  = 1.00 * dot(Id["aa_oo"]("i,j"), dp["aa_oo"]("i,j"));
    scalars_["5"] += 1.00 * dot(Id["bb_oo"]("k,l"), dp["bb_oo"]("k,l"));
    scalars_["6"]  = 1.00 * dot(dp["aa_ov"]("i,a"), t1_2["aa"]("a,i"));
    scalars_["6"] += 1.00 * dot(dp["bb_ov"]("j,b"), t1_2["bb"]("b,j"));

    // crt0_1  = +1.00
    crt0_1  = 1.00;

    // cenergy  = +1.00 t0_1
    cenergy  = t0_1;

    // rt0_2  = -2.00 d+_bb(i,a) t1_1_bb(a,i)
    //       += -2.00 d+_aa(i,a) t1_1_aa(a,i)
    rt0_2  = -2.00 * scalars_["1"];

    // rt0_1  = +0.25 <j,i||a,b>_abab t2_1_abab(a,b,j,i)
    //       += +0.25 <i,j||a,b>_abab t2_1_abab(a,b,i,j)
    //       += +0.25 <j,i||b,a>_abab t2_1_abab(b,a,j,i)
    //       += +0.25 <i,j||b,a>_abab t2_1_abab(b,a,i,j)
    //       += -1.00 d+_aa(i,i)
    //       += +0.25 <j,i||a,b>_aaaa t2_1_aaaa(a,b,j,i)
    //       += -2.00 d-_aa(i,i) t0_2
    //       += -2.00 d-_bb(i,a) t1_2_bb(a,i)
    //       += +1.00 t0_1 w0
    //       += +1.00 f_bb(i,a) t1_1_bb(a,i)
    //       += -1.00 d-_bb(i,a) t0_1 t1_1_bb(a,i)
    //       += -2.00 d-_aa(i,a) t1_2_aa(a,i)
    //       += -2.00 d-_bb(i,i) t0_2
    //       += +0.25 <j,i||a,b>_bbbb t2_1_bbbb(a,b,j,i)
    //       += +1.00 f_aa(i,a) t1_1_aa(a,i)
    //       += -1.00 d+_bb(i,i)
    //       += -1.00 d-_aa(i,a) t0_1 t1_1_aa(a,i)
    rt0_1  = scalars_["3"];

    // energy  = -1.00 d-_bb(i,a) t1_1_bb(a,i)
    //        += -1.00 d-_aa(i,a) t1_1_aa(a,i)
    energy  = -1.00 * scalars_["1"];

    // crt2_bbbb  = +1.00 t2_1_bbbb(a,b,i,j)
    crt2_bbbb("a,b,i,j")  = t2_1["bbbb"]("a,b,i,j");

    // crt2_abab  = +1.00 t2_1_abab(a,b,i,j)
    crt2_abab("a,b,i,j")  = t2_1["abab"]("a,b,i,j");

    // crt2_aaaa  = +1.00 t2_1_aaaa(a,b,i,j)
    crt2_aaaa("a,b,i,j")  = t2_1["aaaa"]("a,b,i,j");

    // crt2_1_bbbb  = +2.00 t2_2_bbbb(a,b,i,j)
    crt2_1_bbbb("a,b,i,j")  = 2.00 * t2_2["bbbb"]("a,b,i,j");

    // crt2_1_abab  = +2.00 t2_2_abab(a,b,i,j)
    crt2_1_abab("a,b,i,j")  = 2.00 * t2_2["abab"]("a,b,i,j");

    // crt2_1_aaaa  = +2.00 t2_2_aaaa(a,b,i,j)
    crt2_1_aaaa("a,b,i,j")  = 2.00 * t2_2["aaaa"]("a,b,i,j");

    // crt1_bb  = +1.00 t1_1_bb(a,i)
    crt1_bb("a,i")  = t1_1["bb"]("a,i");

    // crt1_aa  = +1.00 t1_1_aa(a,i)
    crt1_aa("a,i")  = t1_1["aa"]("a,i");

    // crt1_1_bb  = +2.00 t1_2_bb(a,i)
    crt1_1_bb("a,i")  = 2.00 * t1_2["bb"]("a,i");

    // crt1_1_aa  = +2.00 t1_2_aa(a,i)
    crt1_1_aa("a,i")  = 2.00 * t1_2["aa"]("a,i");

    // crt0_1 += +2.00 t0_2
    crt0_1 += 2.00 * t0_2;

    // energy += -0.50 <j,i||j,i>_aaaa
    //        += -0.50 <j,i||j,i>_abab
    //        += -0.50 <i,j||i,j>_abab
    //        += -0.50 <j,i||j,i>_bbbb
    //        += -1.00 d-_bb(i,i) t0_1
    //        += +0.25 <j,i||a,b>_bbbb t2_bbbb(a,b,j,i)
    //        += +0.25 <j,i||a,b>_abab t2_abab(a,b,j,i)
    //        += +0.25 <i,j||a,b>_abab t2_abab(a,b,i,j)
    //        += +0.25 <j,i||b,a>_abab t2_abab(b,a,j,i)
    //        += +0.25 <i,j||b,a>_abab t2_abab(b,a,i,j)
    //        += +1.00 f_aa(i,i)
    //        += +1.00 f_bb(i,i)
    //        += -1.00 d-_aa(i,i) t0_1
    //        += +0.25 <j,i||a,b>_aaaa t2_aaaa(a,b,j,i)
    energy -= 0.50 * scalars_["2"];

    // rt0_2 += +2.00 f_bb(i,a) t1_2_bb(a,i)
    //       += -1.00 <j,i||a,b>_aaaa t1_1_aa(a,i) t1_1_aa(b,j)
    //       += +1.00 <i,j||a,b>_abab t1_1_aa(a,i) t1_1_bb(b,j)
    //       += +1.00 <j,i||b,a>_abab t1_1_bb(a,i) t1_1_aa(b,j)
    //       += -2.00 d-_bb(i,a) t0_1 t1_2_bb(a,i)
    //       += -4.00 d-_aa(i,a) t1_1_aa(a,i) t0_2
    //       += +0.50 <j,i||a,b>_bbbb t2_2_bbbb(a,b,j,i)
    //       += +4.00 t0_2 w0
    //       += +0.50 <j,i||a,b>_abab t2_2_abab(a,b,j,i)
    //       += +0.50 <i,j||a,b>_abab t2_2_abab(a,b,i,j)
    //       += +0.50 <j,i||b,a>_abab t2_2_abab(b,a,j,i)
    //       += +0.50 <i,j||b,a>_abab t2_2_abab(b,a,i,j)
    //       += +0.50 <j,i||a,b>_aaaa t2_2_aaaa(a,b,j,i)
    //       += -4.00 d-_bb(i,a) t1_1_bb(a,i) t0_2
    //       += +2.00 f_aa(i,a) t1_2_aa(a,i)
    //       += -1.00 <j,i||a,b>_bbbb t1_1_bb(a,i) t1_1_bb(b,j)
    //       += -2.00 d-_aa(i,a) t0_1 t1_2_aa(a,i)
    rt0_2 += 2.00 * scalars_["4"];

    // flops: o1v1  = o1v1
    //  mems: o1v1  = o1v1
    tmps_["15_bb_ov"]("i,a")  = t0_2 * dp["bb_ov"]("i,a");

    // flops: o1v1  = o1v1
    //  mems: o1v1  = o1v1
    tmps_["14_bb_ov"]("i,a")  = t0_1 * dp["bb_ov"]("i,a");

    // flops: o1v1  = o1v1
    //  mems: o1v1  = o1v1
    tmps_["13_aa_ov"]("i,a")  = t0_2 * dp["aa_ov"]("i,a");

    // flops: o1v1  = o1v1
    //  mems: o1v1  = o1v1
    tmps_["12_aa_vo"]("a,i")  = t0_1 * t1_1["aa"]("a,i");

    // flops: o1v1  = o1v1
    //  mems: o1v1  = o1v1
    tmps_["11_aa_ov"]("i,a")  = t0_1 * dp["aa_ov"]("i,a");

    // flops: o0v2  = o0v2
    //  mems: o0v2  = o0v2
    tmps_["10_aa_vv"]("a,b")  = t0_1 * dp["aa_vv"]("a,b");

    // flops: o2v0  = o2v1
    //  mems: o2v0  = o2v0
    tmps_["9_aa_oo"]("i,j")  = dp["aa_ov"]("i,a") * t1_1["aa"]("a,j");

    // flops: o1v1  = o2v2
    //  mems: o1v1  = o1v1
    tmps_["8_bb_ov"]("i,a")  = eri["abab_oovv"]("j,i,b,a") * t1_1["aa"]("b,j");

    // flops: o2v0  = o3v2
    //  mems: o2v0  = o2v0
    tmps_["7_aa_oo"]("i,j")  = eri["abab_oovv"]("i,k,a,b") * t2["abab"]("a,b,j,k");

    // flops: o2v0  = o3v2
    //  mems: o2v0  = o2v0
    tmps_["6_aa_oo"]("i,j")  = eri["aaaa_oovv"]("k,i,a,b") * t2["aaaa"]("a,b,j,k");

    // flops: o0v2  = o2v3
    //  mems: o0v2  = o0v2
    tmps_["5_aa_vv"]("a,b")  = eri["abab_oovv"]("i,j,a,c") * t2["abab"]("b,c,i,j");

    // flops: o0v2  = o2v3
    //  mems: o0v2  = o0v2
    tmps_["4_aa_vv"]("a,b")  = eri["aaaa_oovv"]("i,j,c,a") * t2["aaaa"]("c,b,j,i");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["3_bbaa_ovvo"]("i,a,b,j")  = eri["bbbb_oovv"]("k,i,c,a") * t2["abab"]("b,c,j,k");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["2_bbaa_ovvo"]("i,a,b,j")  = eri["abab_oovv"]("k,i,c,a") * t2["aaaa"]("c,b,j,k");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["1_aaaa_ovvo"]("i,a,b,j")  = eri["aaaa_oovv"]("k,i,c,a") * t2["aaaa"]("c,b,j,k");

    // flops: o1v1  = o2v2 o2v3 o1v1 o2v2 o1v1 o2v1 o1v1 o1v1 o1v1 o1v1 o2v2 o1v1 o1v1 o1v1 o2v3 o1v1 o1v1 o1v1 o3v2 o1v1 o1v1 o1v1 o1v2 o1v1 o3v2 o1v1 o2v2 o1v1 o2v2 o2v2 o1v1 o2v2 o2v1 o1v1 o1v1 o1v2 o1v1 o2v2 o2v2 o1v1 o1v1 o2v2 o1v1 o1v1 o1v2 o1v1 o2v2 o1v1 o1v2 o1v1 o2v2 o2v2 o1v1 o1v1 o2v1 o1v1 o1v2 o1v1 o2v1 o1v1 o2v2 o1v1 o2v2 o1v1 o2v1 o1v1 o2v2 o1v1 o2v1 o1v1
    //  mems: o1v1  = o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1
    tmps_["16_aa_vo"]("a,i")  = -0.50 * eri["aaaa_vovo"]("a,j,b,i") * t1_1["aa"]("b,j");
    tmps_["16_aa_vo"]("a,i") += 0.25 * eri["aaaa_vovv"]("a,j,b,e") * t2_1["aaaa"]("b,e,i,j");
    tmps_["16_aa_vo"]("a,i") -= 0.50 * f["aa_ov"]("j,b") * t2_1["aaaa"]("b,a,i,j");
    tmps_["16_aa_vo"]("a,i") -= 0.50 * f["aa_oo"]("j,i") * t1_1["aa"]("a,j");
    tmps_["16_aa_vo"]("a,i") -= 0.50 * scalars_["1"] * t1_1["aa"]("a,i");
    tmps_["16_aa_vo"]("a,i") -= 0.50 * dp["aa_vo"]("a,i");
    tmps_["16_aa_vo"]("a,i") += 0.50 * f["bb_ov"]("k,c") * t2_1["abab"]("a,c,i,k");
    tmps_["16_aa_vo"]("a,i") += 0.50 * w0 * t1_1["aa"]("a,i");
    tmps_["16_aa_vo"]("a,i") += 0.50 * eri["abab_vovv"]("a,k,b,d") * t2_1["abab"]("b,d,i,k");
    tmps_["16_aa_vo"]("a,i") -= scalars_["5"] * t1_2["aa"]("a,i");
    tmps_["16_aa_vo"]("a,i") += 0.50 * eri["abba_oovo"]("l,k,c,i") * t2_1["abab"]("a,c,l,k");
    tmps_["16_aa_vo"]("a,i") -= t0_2 * dp["aa_vo"]("a,i");
    tmps_["16_aa_vo"]("a,i") += 0.50 * f["aa_vv"]("a,b") * t1_1["aa"]("b,i");
    tmps_["16_aa_vo"]("a,i") += 0.25 * eri["aaaa_oovo"]("j,l,b,i") * t2_1["aaaa"]("b,a,l,j");
    tmps_["16_aa_vo"]("a,i") -= 0.50 * eri["abba_vovo"]("a,k,c,i") * t1_1["bb"]("c,k");
    tmps_["16_aa_vo"]("a,i") += 0.50 * dp["aa_ov"]("j,b") * t2["aaaa"]("b,a,i,j");
    tmps_["16_aa_vo"]("a,i") -= 0.50 * dp["bb_ov"]("k,c") * t2["abab"]("a,c,i,k");
    tmps_["16_aa_vo"]("a,i") -= 0.50 * t1_1["aa"]("e,l") * tmps_["1_aaaa_ovvo"]("l,e,a,i");
    tmps_["16_aa_vo"]("a,i") += 0.50 * dp["aa_oo"]("j,i") * tmps_["12_aa_vo"]("a,j");
    tmps_["16_aa_vo"]("a,i") += 0.25 * t1_1["aa"]("e,i") * tmps_["4_aa_vv"]("e,a");
    tmps_["16_aa_vo"]("a,i") -= dp["bb_ov"]("k,c") * t2_2["abab"]("a,c,i,k");
    tmps_["16_aa_vo"]("a,i") += dp["aa_ov"]("j,b") * t2_2["aaaa"]("b,a,i,j");
    tmps_["16_aa_vo"]("a,i") += 0.50 * t1_1["bb"]("d,m") * tmps_["3_bbaa_ovvo"]("m,d,a,i");
    tmps_["16_aa_vo"]("a,i") -= 0.50 * t1_1["aa"]("b,i") * tmps_["10_aa_vv"]("a,b");
    tmps_["16_aa_vo"]("a,i") -= 0.50 * t1_1["bb"]("d,m") * tmps_["2_bbaa_ovvo"]("m,d,a,i");
    tmps_["16_aa_vo"]("a,i") -= 0.50 * t1_1["aa"]("e,i") * tmps_["5_aa_vv"]("e,a");
    tmps_["16_aa_vo"]("a,i") += 0.50 * tmps_["8_bb_ov"]("k,c") * t2["abab"]("a,c,i,k");
    tmps_["16_aa_vo"]("a,i") += 0.50 * tmps_["11_aa_ov"]("j,b") * t2_1["aaaa"]("b,a,i,j");
    tmps_["16_aa_vo"]("a,i") += 0.25 * tmps_["6_aa_oo"]("l,i") * t1_1["aa"]("a,l");
    tmps_["16_aa_vo"]("a,i") -= dp["aa_vv"]("a,b") * t1_2["aa"]("b,i");
    tmps_["16_aa_vo"]("a,i") -= 0.50 * tmps_["7_aa_oo"]("l,i") * t1_1["aa"]("a,l");
    tmps_["16_aa_vo"]("a,i") -= tmps_["15_bb_ov"]("k,c") * t2["abab"]("a,c,i,k");
    tmps_["16_aa_vo"]("a,i") -= 0.50 * tmps_["14_bb_ov"]("k,c") * t2_1["abab"]("a,c,i,k");
    tmps_["16_aa_vo"]("a,i") += dp["aa_oo"]("j,i") * t1_2["aa"]("a,j");
    tmps_["16_aa_vo"]("a,i") += tmps_["13_aa_ov"]("j,b") * t2["aaaa"]("b,a,i,j");
    tmps_["16_aa_vo"]("a,i") += tmps_["9_aa_oo"]("j,i") * t1_1["aa"]("a,j");

    // rt1_1_aa  = +2.00 d-_aa(j,b) t1_1_aa(a,j) t1_1_aa(b,i)
    //          += +1.00 d-_aa(j,b) t0_1 t2_1_aaaa(b,a,i,j)
    //          += +1.00 <k,j||c,b>_abab t2_abab(a,b,i,j) t1_1_aa(c,k)
    //          += -0.50 <k,j||b,c>_aaaa t2_aaaa(b,a,k,j) t1_1_aa(c,i)
    //          += +1.00 d-_aa(j,i) t0_1 t1_1_aa(a,j)
    //          += +1.00 <k,j||b,c>_aaaa t2_aaaa(b,a,i,j) t1_1_aa(c,k)
    //          += -1.00 d+_bb(j,b) t2_abab(a,b,i,j)
    //          += +1.00 d+_aa(j,b) t2_aaaa(b,a,i,j)
    //          += +2.00 d-_aa(j,b) t2_2_aaaa(b,a,i,j)
    //          += -2.00 d-_bb(j,b) t2_2_abab(a,b,i,j)
    //          += -1.00 <k,j||b,c>_bbbb t2_abab(a,b,i,j) t1_1_bb(c,k)
    //          += -0.50 <j,a||b,c>_aaaa t2_1_aaaa(b,c,i,j)
    //          += +1.00 <j,a||b,i>_aaaa t1_1_aa(b,j)
    //          += -1.00 f_aa(j,b) t2_1_aaaa(b,a,i,j)
    //          += -1.00 f_aa(j,i) t1_1_aa(a,j)
    //          += -1.00 d-_bb(j,b) t1_1_aa(a,i) t1_1_bb(b,j)
    //          += -1.00 d-_aa(j,b) t1_1_aa(a,i) t1_1_aa(b,j)
    //          += +1.00 f_bb(j,b) t2_1_abab(a,b,i,j)
    //          += -1.00 d+_aa(a,i)
    //          += +1.00 t1_1_aa(a,i) w0
    //          += +0.50 <a,j||b,c>_abab t2_1_abab(b,c,i,j)
    //          += +0.50 <a,j||c,b>_abab t2_1_abab(c,b,i,j)
    //          += -2.00 d-_aa(j,j) t1_2_aa(a,i)
    //          += -2.00 d-_bb(j,j) t1_2_aa(a,i)
    //          += -0.50 <k,j||i,b>_abab t2_1_abab(a,b,k,j)
    //          += -0.50 <j,k||i,b>_abab t2_1_abab(a,b,j,k)
    //          += -2.00 d-_aa(a,i) t0_2
    //          += +1.00 f_aa(a,b) t1_1_aa(b,i)
    //          += -0.50 <k,j||b,i>_aaaa t2_1_aaaa(b,a,k,j)
    //          += +1.00 <a,j||i,b>_abab t1_1_bb(b,j)
    //          += -1.00 d-_aa(a,b) t0_1 t1_1_aa(b,i)
    //          += -1.00 <j,k||b,c>_abab t2_aaaa(b,a,i,j) t1_1_bb(c,k)
    //          += -0.50 <k,j||c,b>_abab t2_abab(a,b,k,j) t1_1_aa(c,i)
    //          += -0.50 <j,k||c,b>_abab t2_abab(a,b,j,k) t1_1_aa(c,i)
    //          += -0.50 <k,j||b,c>_aaaa t2_aaaa(b,c,i,j) t1_1_aa(a,k)
    //          += +2.00 d-_aa(j,i) t1_2_aa(a,j)
    //          += -2.00 d-_aa(a,b) t1_2_aa(b,i)
    //          += -0.50 <k,j||b,c>_abab t2_abab(b,c,i,j) t1_1_aa(a,k)
    //          += -0.50 <k,j||c,b>_abab t2_abab(c,b,i,j) t1_1_aa(a,k)
    //          += -2.00 d-_bb(j,b) t2_abab(a,b,i,j) t0_2
    //          += -1.00 d-_bb(j,b) t0_1 t2_1_abab(a,b,i,j)
    //          += +2.00 d-_aa(j,b) t2_aaaa(b,a,i,j) t0_2
    rt1_1_aa("a,i")  = 2.00 * tmps_["16_aa_vo"]("a,i");
    tmps_["16_aa_vo"].~TArrayD();

    // flops: o1v1  = o1v1
    //  mems: o1v1  = o1v1
    tmps_["26_bb_vo"]("a,i")  = t0_1 * t1_1["bb"]("a,i");

    // flops: o0v2  = o0v2
    //  mems: o0v2  = o0v2
    tmps_["25_bb_vv"]("a,b")  = t0_1 * dp["bb_vv"]("a,b");

    // flops: o2v0  = o2v1
    //  mems: o2v0  = o2v0
    tmps_["24_bb_oo"]("i,j")  = dp["bb_ov"]("i,a") * t1_1["bb"]("a,j");

    // flops: o2v0  = o3v2
    //  mems: o2v0  = o2v0
    tmps_["23_bb_oo"]("i,j")  = eri["bbbb_oovv"]("k,i,a,b") * t2["bbbb"]("a,b,j,k");

    // flops: o2v0  = o3v2
    //  mems: o2v0  = o2v0
    tmps_["22_bb_oo"]("i,j")  = eri["abab_oovv"]("k,i,a,b") * t2["abab"]("a,b,k,j");

    // flops: o0v2  = o2v3
    //  mems: o0v2  = o0v2
    tmps_["21_bb_vv"]("a,b")  = eri["bbbb_oovv"]("i,j,c,a") * t2["bbbb"]("c,b,j,i");

    // flops: o0v2  = o2v3
    //  mems: o0v2  = o0v2
    tmps_["20_bb_vv"]("a,b")  = eri["abab_oovv"]("i,j,c,a") * t2["abab"]("c,b,i,j");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["19_bbbb_ovvo"]("i,a,b,j")  = eri["bbbb_oovv"]("k,i,c,a") * t2["bbbb"]("c,b,j,k");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["18_bbbb_ovvo"]("i,a,b,j")  = eri["abab_oovv"]("k,i,c,a") * t2["abab"]("c,b,k,j");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["17_aabb_ovvo"]("i,a,b,j")  = eri["aaaa_oovv"]("k,i,c,a") * t2["abab"]("c,b,k,j");

    // flops: o1v1  = o2v2 o1v1 o1v1 o1v1 o1v1 o3v2 o1v1 o2v2 o1v1 o1v1 o2v1 o1v1 o1v2 o1v1 o1v1 o1v1 o2v2 o1v1 o2v3 o1v1 o3v2 o1v1 o1v1 o1v1 o2v2 o1v1 o2v3 o1v1 o2v2 o1v2 o1v1 o1v2 o1v1 o1v1 o2v1 o1v1 o2v2 o2v2 o1v1 o1v1 o2v2 o1v1 o2v2 o2v2 o1v1 o1v1 o1v2 o1v1 o2v2 o1v1 o2v2 o2v1 o1v1 o2v2 o1v1 o2v2 o1v1 o2v1 o1v1 o2v2 o1v1 o2v2 o1v1 o1v1 o2v1 o1v1 o1v2 o1v1 o2v1 o1v1
    //  mems: o1v1  = o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1
    tmps_["27_bb_vo"]("a,i")  = 0.50 * f["aa_ov"]("l,d") * t2_1["abab"]("d,a,l,i");
    tmps_["27_bb_vo"]("a,i") -= t0_2 * dp["bb_vo"]("a,i");
    tmps_["27_bb_vo"]("a,i") -= scalars_["5"] * t1_2["bb"]("a,i");
    tmps_["27_bb_vo"]("a,i") += 0.25 * eri["bbbb_oovo"]("j,k,b,i") * t2_1["bbbb"]("b,a,k,j");
    tmps_["27_bb_vo"]("a,i") -= 0.50 * eri["baab_vovo"]("a,l,d,i") * t1_1["aa"]("d,l");
    tmps_["27_bb_vo"]("a,i") -= 0.50 * dp["bb_vo"]("a,i");
    tmps_["27_bb_vo"]("a,i") -= 0.50 * f["bb_oo"]("j,i") * t1_1["bb"]("a,j");
    tmps_["27_bb_vo"]("a,i") += 0.50 * f["bb_vv"]("a,b") * t1_1["bb"]("b,i");
    tmps_["27_bb_vo"]("a,i") += 0.50 * w0 * t1_1["bb"]("a,i");
    tmps_["27_bb_vo"]("a,i") -= 0.50 * f["bb_ov"]("j,b") * t2_1["bbbb"]("b,a,i,j");
    tmps_["27_bb_vo"]("a,i") += 0.25 * eri["bbbb_vovv"]("a,j,b,c") * t2_1["bbbb"]("b,c,i,j");
    tmps_["27_bb_vo"]("a,i") -= 0.50 * eri["abab_oovo"]("m,j,d,i") * t2_1["abab"]("d,a,m,j");
    tmps_["27_bb_vo"]("a,i") -= 0.50 * scalars_["1"] * t1_1["bb"]("a,i");
    tmps_["27_bb_vo"]("a,i") -= 0.50 * eri["bbbb_vovo"]("a,j,b,i") * t1_1["bb"]("b,j");
    tmps_["27_bb_vo"]("a,i") -= 0.50 * eri["baab_vovv"]("a,l,d,c") * t2_1["abab"]("d,c,l,i");
    tmps_["27_bb_vo"]("a,i") += 0.50 * t1_1["bb"]("c,k") * tmps_["18_bbbb_ovvo"]("k,c,a,i");
    tmps_["27_bb_vo"]("a,i") -= 0.50 * t1_1["bb"]("c,i") * tmps_["20_bb_vv"]("c,a");
    tmps_["27_bb_vo"]("a,i") -= 0.50 * t1_1["bb"]("b,i") * tmps_["25_bb_vv"]("a,b");
    tmps_["27_bb_vo"]("a,i") += 0.50 * dp["bb_oo"]("j,i") * tmps_["26_bb_vo"]("a,j");
    tmps_["27_bb_vo"]("a,i") += 0.50 * dp["bb_ov"]("j,b") * t2["bbbb"]("b,a,i,j");
    tmps_["27_bb_vo"]("a,i") -= 0.50 * dp["aa_ov"]("l,d") * t2["abab"]("d,a,l,i");
    tmps_["27_bb_vo"]("a,i") += 0.50 * t1_1["aa"]("e,m") * tmps_["17_aabb_ovvo"]("m,e,a,i");
    tmps_["27_bb_vo"]("a,i") -= dp["aa_ov"]("l,d") * t2_2["abab"]("d,a,l,i");
    tmps_["27_bb_vo"]("a,i") += dp["bb_ov"]("j,b") * t2_2["bbbb"]("b,a,i,j");
    tmps_["27_bb_vo"]("a,i") += 0.25 * t1_1["bb"]("c,i") * tmps_["21_bb_vv"]("c,a");
    tmps_["27_bb_vo"]("a,i") -= 0.50 * t1_1["bb"]("c,k") * tmps_["19_bbbb_ovvo"]("k,c,a,i");
    tmps_["27_bb_vo"]("a,i") += 0.50 * tmps_["14_bb_ov"]("j,b") * t2_1["bbbb"]("b,a,i,j");
    tmps_["27_bb_vo"]("a,i") -= 0.50 * tmps_["22_bb_oo"]("k,i") * t1_1["bb"]("a,k");
    tmps_["27_bb_vo"]("a,i") += tmps_["15_bb_ov"]("j,b") * t2["bbbb"]("b,a,i,j");
    tmps_["27_bb_vo"]("a,i") -= 0.50 * tmps_["8_bb_ov"]("j,b") * t2["bbbb"]("b,a,i,j");
    tmps_["27_bb_vo"]("a,i") += 0.25 * tmps_["23_bb_oo"]("k,i") * t1_1["bb"]("a,k");
    tmps_["27_bb_vo"]("a,i") -= 0.50 * tmps_["11_aa_ov"]("l,d") * t2_1["abab"]("d,a,l,i");
    tmps_["27_bb_vo"]("a,i") -= tmps_["13_aa_ov"]("l,d") * t2["abab"]("d,a,l,i");
    tmps_["27_bb_vo"]("a,i") += dp["bb_oo"]("j,i") * t1_2["bb"]("a,j");
    tmps_["27_bb_vo"]("a,i") -= dp["bb_vv"]("a,b") * t1_2["bb"]("b,i");
    tmps_["27_bb_vo"]("a,i") += tmps_["24_bb_oo"]("j,i") * t1_1["bb"]("a,j");

    // rt1_1_bb  = +2.00 d-_bb(j,b) t1_1_bb(a,j) t1_1_bb(b,i)
    //          += -0.50 <j,k||b,c>_abab t2_abab(b,c,j,i) t1_1_bb(a,k)
    //          += -0.50 <j,k||c,b>_abab t2_abab(c,b,j,i) t1_1_bb(a,k)
    //          += +1.00 d-_bb(j,b) t0_1 t2_1_bbbb(b,a,i,j)
    //          += +2.00 d-_bb(j,b) t2_bbbb(b,a,i,j) t0_2
    //          += -1.00 <k,j||c,b>_abab t2_bbbb(b,a,i,j) t1_1_aa(c,k)
    //          += -1.00 d-_aa(j,b) t0_1 t2_1_abab(b,a,j,i)
    //          += -0.50 <k,j||b,c>_bbbb t2_bbbb(b,c,i,j) t1_1_bb(a,k)
    //          += -2.00 d-_aa(j,b) t2_abab(b,a,j,i) t0_2
    //          += -0.50 <k,j||b,c>_abab t2_abab(b,a,k,j) t1_1_bb(c,i)
    //          += -0.50 <j,k||b,c>_abab t2_abab(b,a,j,k) t1_1_bb(c,i)
    //          += +1.00 <j,k||b,c>_abab t2_abab(b,a,j,i) t1_1_bb(c,k)
    //          += -1.00 d-_bb(a,b) t0_1 t1_1_bb(b,i)
    //          += -2.00 d-_bb(a,i) t0_2
    //          += +1.00 f_aa(j,b) t2_1_abab(b,a,j,i)
    //          += -2.00 d-_aa(j,j) t1_2_bb(a,i)
    //          += -2.00 d-_bb(j,j) t1_2_bb(a,i)
    //          += +1.00 <j,a||b,i>_abab t1_1_aa(b,j)
    //          += -0.50 <k,j||b,i>_bbbb t2_1_bbbb(b,a,k,j)
    //          += -1.00 d+_bb(a,i)
    //          += -1.00 f_bb(j,i) t1_1_bb(a,j)
    //          += +1.00 f_bb(a,b) t1_1_bb(b,i)
    //          += +1.00 t1_1_bb(a,i) w0
    //          += -1.00 f_bb(j,b) t2_1_bbbb(b,a,i,j)
    //          += -0.50 <k,j||b,i>_abab t2_1_abab(b,a,k,j)
    //          += -0.50 <j,k||b,i>_abab t2_1_abab(b,a,j,k)
    //          += -0.50 <j,a||b,c>_bbbb t2_1_bbbb(b,c,i,j)
    //          += -1.00 d-_bb(j,b) t1_1_bb(a,i) t1_1_bb(b,j)
    //          += -1.00 d-_aa(j,b) t1_1_bb(a,i) t1_1_aa(b,j)
    //          += +1.00 <j,a||b,i>_bbbb t1_1_bb(b,j)
    //          += +0.50 <j,a||b,c>_abab t2_1_abab(b,c,j,i)
    //          += +0.50 <j,a||c,b>_abab t2_1_abab(c,b,j,i)
    //          += +1.00 <k,j||b,c>_bbbb t2_bbbb(b,a,i,j) t1_1_bb(c,k)
    //          += +1.00 d-_bb(j,i) t0_1 t1_1_bb(a,j)
    //          += -1.00 d+_aa(j,b) t2_abab(b,a,j,i)
    //          += +1.00 d+_bb(j,b) t2_bbbb(b,a,i,j)
    //          += -1.00 <k,j||b,c>_aaaa t2_abab(b,a,j,i) t1_1_aa(c,k)
    //          += +2.00 d-_bb(j,b) t2_2_bbbb(b,a,i,j)
    //          += -2.00 d-_aa(j,b) t2_2_abab(b,a,j,i)
    //          += -0.50 <k,j||b,c>_bbbb t2_bbbb(b,a,k,j) t1_1_bb(c,i)
    //          += -2.00 d-_bb(a,b) t1_2_bb(b,i)
    //          += +2.00 d-_bb(j,i) t1_2_bb(a,j)
    rt1_1_bb("a,i")  = 2.00 * tmps_["27_bb_vo"]("a,i");
    tmps_["27_bb_vo"].~TArrayD();

    // flops: o1v1  = o2v2
    //  mems: o1v1  = o1v1
    tmps_["44_aa_ov"]("i,a")  = eri["aaaa_oovv"]("j,i,b,a") * t1_1["aa"]("b,j");

    // flops: o1v1  = o2v2
    //  mems: o1v1  = o1v1
    tmps_["43_bb_ov"]("i,a")  = eri["bbbb_oovv"]("j,i,b,a") * t1_1["bb"]("b,j");

    // flops: o2v0  = o2v0
    //  mems: o2v0  = o2v0
    tmps_["42_aa_oo"]("i,j")  = t0_2 * dp["aa_oo"]("i,j");

    // flops: o2v0  = o2v0
    //  mems: o2v0  = o2v0
    tmps_["41_aa_oo"]("i,j")  = t0_1 * dp["aa_oo"]("i,j");

    // flops: o0v2  = o0v2
    //  mems: o0v2  = o0v2
    tmps_["40_aa_vv"]("a,b")  = t0_2 * dp["aa_vv"]("a,b");

    // flops: o2v0  = o2v1
    //  mems: o2v0  = o2v0
    tmps_["39_aa_oo"]("i,j")  = f["aa_ov"]("i,a") * t1_1["aa"]("a,j");

    // flops: o2v0  = o2v1
    //  mems: o2v0  = o2v0
    tmps_["38_aa_oo"]("i,j")  = dp["aa_ov"]("i,a") * t1_2["aa"]("a,j");

    // flops: o2v0  = o3v1
    //  mems: o2v0  = o2v0
    tmps_["37_aa_oo"]("i,j")  = eri["abba_oovo"]("i,k,a,j") * t1_1["bb"]("a,k");

    // flops: o2v0  = o3v1
    //  mems: o2v0  = o2v0
    tmps_["36_aa_oo"]("i,j")  = eri["aaaa_oovo"]("k,i,a,j") * t1_1["aa"]("a,k");

    // flops: o1v1  = o2v2
    //  mems: o1v1  = o1v1
    tmps_["35_bb_ov"]("i,a")  = eri["abab_oovv"]("j,i,b,a") * t1_2["aa"]("b,j");

    // flops: o2v0  = o3v2
    //  mems: o2v0  = o2v0
    tmps_["34_aa_oo"]("i,j")  = eri["abab_oovv"]("i,k,a,b") * t2_1["abab"]("a,b,j,k");

    // flops: o2v0  = o3v2
    //  mems: o2v0  = o2v0
    tmps_["33_aa_oo"]("i,j")  = eri["aaaa_oovv"]("i,k,a,b") * t2_1["aaaa"]("a,b,j,k");

    // flops: o0v2  = o2v3
    //  mems: o0v2  = o0v2
    tmps_["32_aa_vv"]("a,b")  = eri["abab_oovv"]("i,j,a,c") * t2_1["abab"]("b,c,i,j");

    // flops: o2v2  = o2v3
    //  mems: o2v2  = o2v2
    tmps_["31_abba_vovo"]("a,i,b,j")  = eri["abab_vovv"]("a,i,c,b") * t1_1["aa"]("c,j");

    // flops: o2v2  = o2v3
    //  mems: o2v2  = o2v2
    tmps_["30_aaaa_vovo"]("a,i,b,j")  = eri["aaaa_vovv"]("a,i,b,c") * t1_1["aa"]("c,j");

    // flops: o0v2  = o2v3
    //  mems: o0v2  = o0v2
    tmps_["29_aa_vv"]("a,b")  = eri["aaaa_oovv"]("i,j,a,c") * t2_1["aaaa"]("c,b,j,i");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["28_bbaa_ovvo"]("i,a,b,j")  = eri["abab_oovv"]("k,i,c,a") * t2_1["aaaa"]("c,b,j,k");

    // flops: o1v1  = o2v1 o2v2 o1v1 o2v1 o1v1 o2v1 o1v1 o1v2 o1v1 o1v1 o2v2 o1v1 o2v3 o1v1 o2v2 o1v1 o2v3 o1v1 o3v2 o1v1 o2v2 o1v1 o1v1 o1v1 o3v2 o1v1 o1v1 o1v1 o2v2 o1v1 o2v1 o1v1 o2v1 o2v2 o1v1 o1v2 o1v1 o1v2 o1v1 o1v2 o1v1 o1v1 o2v1 o1v1 o1v2 o1v1 o2v2 o1v1 o2v2 o1v1 o1v2 o1v1 o2v1 o1v1 o2v2 o1v1 o2v2 o1v1 o2v2 o1v1 o1v2 o1v1 o2v2 o1v1 o1v1 o2v2 o1v1 o2v2 o1v1 o2v2 o1v1 o2v1 o1v1 o2v2 o1v1 o2v1 o1v1 o2v2 o1v1 o2v1 o1v1 o2v1 o1v1 o2v1 o1v1 o2v2 o1v1 o2v1 o1v1
    //  mems: o1v1  = o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1
    tmps_["45_aa_vo"]("a,i")  = -0.17 * tmps_["33_aa_oo"]("j,i") * t1_1["aa"]("a,j");
    tmps_["45_aa_vo"]("a,i") += 0.33 * t2_1["abab"]("a,e,i,m") * tmps_["8_bb_ov"]("m,e");
    tmps_["45_aa_vo"]("a,i") += tmps_["38_aa_oo"]("j,i") * t1_1["aa"]("a,j");
    tmps_["45_aa_vo"]("a,i") -= 0.33 * tmps_["34_aa_oo"]("j,i") * t1_1["aa"]("a,j");
    tmps_["45_aa_vo"]("a,i") += 0.33 * f["aa_vv"]("a,c") * t1_2["aa"]("c,i");
    tmps_["45_aa_vo"]("a,i") -= 0.33 * scalars_["6"] * t1_1["aa"]("a,i");
    tmps_["45_aa_vo"]("a,i") -= 0.33 * eri["abba_vovo"]("a,k,b,i") * t1_2["bb"]("b,k");
    tmps_["45_aa_vo"]("a,i") += 0.33 * eri["abab_vovv"]("a,k,c,e") * t2_2["abab"]("c,e,i,k");
    tmps_["45_aa_vo"]("a,i") += 0.33 * f["bb_ov"]("k,b") * t2_2["abab"]("a,b,i,k");
    tmps_["45_aa_vo"]("a,i") += 0.1650 * eri["aaaa_vovv"]("a,j,c,d") * t2_2["aaaa"]("c,d,i,j");
    tmps_["45_aa_vo"]("a,i") += 0.33 * eri["abba_oovo"]("l,k,b,i") * t2_2["abab"]("a,b,l,k");
    tmps_["45_aa_vo"]("a,i") -= 0.33 * eri["aaaa_vovo"]("a,j,c,i") * t1_2["aa"]("c,j");
    tmps_["45_aa_vo"]("a,i") += 0.66 * w0 * t1_2["aa"]("a,i");
    tmps_["45_aa_vo"]("a,i") += 0.1650 * eri["aaaa_oovo"]("j,l,c,i") * t2_2["aaaa"]("c,a,l,j");
    tmps_["45_aa_vo"]("a,i") -= 0.66 * scalars_["1"] * t1_2["aa"]("a,i");
    tmps_["45_aa_vo"]("a,i") -= 0.33 * f["aa_ov"]("j,c") * t2_2["aaaa"]("c,a,i,j");
    tmps_["45_aa_vo"]("a,i") -= 0.33 * f["aa_oo"]("j,i") * t1_2["aa"]("a,j");
    tmps_["45_aa_vo"]("a,i") += 0.33 * tmps_["12_aa_vo"]("a,j") * tmps_["9_aa_oo"]("j,i");
    tmps_["45_aa_vo"]("a,i") -= 0.33 * t1_2["aa"]("d,l") * tmps_["1_aaaa_ovvo"]("l,d,a,i");
    tmps_["45_aa_vo"]("a,i") -= 0.33 * t1_2["aa"]("c,i") * tmps_["10_aa_vv"]("a,c");
    tmps_["45_aa_vo"]("a,i") -= 0.33 * t1_2["aa"]("d,i") * tmps_["5_aa_vv"]("d,a");
    tmps_["45_aa_vo"]("a,i") -= 0.33 * t1_1["aa"]("c,i") * tmps_["32_aa_vv"]("c,a");
    tmps_["45_aa_vo"]("a,i") += 0.66 * t1_1["aa"]("a,j") * tmps_["42_aa_oo"]("j,i");
    tmps_["45_aa_vo"]("a,i") -= 0.1650 * t1_1["aa"]("c,i") * tmps_["29_aa_vv"]("c,a");
    tmps_["45_aa_vo"]("a,i") -= 0.33 * t1_2["bb"]("e,m") * tmps_["2_bbaa_ovvo"]("m,e,a,i");
    tmps_["45_aa_vo"]("a,i") += 0.33 * t1_1["bb"]("b,k") * tmps_["31_abba_vovo"]("a,k,b,i");
    tmps_["45_aa_vo"]("a,i") += 0.1650 * t1_2["aa"]("d,i") * tmps_["4_aa_vv"]("d,a");
    tmps_["45_aa_vo"]("a,i") += 0.33 * t1_2["aa"]("a,j") * tmps_["41_aa_oo"]("j,i");
    tmps_["45_aa_vo"]("a,i") -= 0.33 * tmps_["44_aa_ov"]("l,d") * t2_1["aaaa"]("d,a,i,l");
    tmps_["45_aa_vo"]("a,i") -= 0.33 * t1_1["aa"]("c,j") * tmps_["30_aaaa_vovo"]("a,j,c,i");
    tmps_["45_aa_vo"]("a,i") += 0.33 * t1_2["bb"]("e,m") * tmps_["3_bbaa_ovvo"]("m,e,a,i");
    tmps_["45_aa_vo"]("a,i") -= 0.66 * t1_1["aa"]("c,i") * tmps_["40_aa_vv"]("a,c");
    tmps_["45_aa_vo"]("a,i") -= 0.33 * t1_1["bb"]("b,k") * tmps_["28_bbaa_ovvo"]("k,b,a,i");
    tmps_["45_aa_vo"]("a,i") += 0.67 * t2_1["aaaa"]("c,a,i,j") * tmps_["13_aa_ov"]("j,c");
    tmps_["45_aa_vo"]("a,i") -= 0.33 * t2_2["abab"]("a,b,i,k") * tmps_["14_bb_ov"]("k,b");
    tmps_["45_aa_vo"]("a,i") += 0.33 * t2_1["abab"]("a,e,i,m") * tmps_["43_bb_ov"]("m,e");
    tmps_["45_aa_vo"]("a,i") -= 0.33 * tmps_["36_aa_oo"]("l,i") * t1_1["aa"]("a,l");
    tmps_["45_aa_vo"]("a,i") += 0.33 * t2_2["aaaa"]("c,a,i,j") * tmps_["11_aa_ov"]("j,c");
    tmps_["45_aa_vo"]("a,i") += tmps_["9_aa_oo"]("j,i") * t1_2["aa"]("a,j");
    tmps_["45_aa_vo"]("a,i") -= 0.67 * t2_1["abab"]("a,b,i,k") * tmps_["15_bb_ov"]("k,b");
    tmps_["45_aa_vo"]("a,i") += 0.33 * tmps_["37_aa_oo"]("l,i") * t1_1["aa"]("a,l");
    tmps_["45_aa_vo"]("a,i") -= 0.33 * tmps_["7_aa_oo"]("l,i") * t1_2["aa"]("a,l");
    tmps_["45_aa_vo"]("a,i") += 0.17 * tmps_["6_aa_oo"]("l,i") * t1_2["aa"]("a,l");
    tmps_["45_aa_vo"]("a,i") += 0.33 * t2["abab"]("a,b,i,k") * tmps_["35_bb_ov"]("k,b");
    tmps_["45_aa_vo"]("a,i") -= 0.33 * tmps_["39_aa_oo"]("j,i") * t1_1["aa"]("a,j");
    tmps_["44_aa_ov"].~TArrayD();
    tmps_["43_bb_ov"].~TArrayD();
}
#endif
