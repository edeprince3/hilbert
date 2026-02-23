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

void QED_CCSD_22::resid_22_2() {

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


        // rt1_2_aa  = +6.00 d-_aa(j,b) t1_1_aa(a,j) t1_2_aa(b,i)
        //          += +2.00 <j,k||b,c>_abab t2_1_abab(a,c,i,k) t1_1_aa(b,j)
        //          += +1.00 <k,j||b,c>_aaaa t1_1_aa(a,j) t2_1_aaaa(b,c,i,k)
        //          += -1.00 <j,k||b,c>_abab t1_1_aa(a,j) t2_1_abab(b,c,i,k)
        //          += -1.00 <j,k||c,b>_abab t1_1_aa(a,j) t2_1_abab(c,b,i,k)
        //          += +2.00 <k,j||b,c>_aaaa t2_aaaa(b,a,i,j) t1_2_aa(c,k)
        //          += +2.00 d-_aa(j,b) t0_1 t1_1_aa(a,j) t1_1_aa(b,i)
        //          += -2.00 d-_aa(a,b) t0_1 t1_2_aa(b,i)
        //          += -1.00 <k,j||c,b>_abab t2_abab(a,b,k,j) t1_2_aa(c,i)
        //          += -1.00 <j,k||c,b>_abab t2_abab(a,b,j,k) t1_2_aa(c,i)
        //          += -1.00 <k,j||b,c>_abab t2_1_abab(a,c,k,j) t1_1_aa(b,i)
        //          += -1.00 <j,k||b,c>_abab t2_1_abab(a,c,j,k) t1_1_aa(b,i)
        //          += -2.00 d-_aa(j,b) t1_1_aa(a,i) t1_2_aa(b,j)
        //          += -2.00 d-_bb(j,b) t1_1_aa(a,i) t1_2_bb(b,j)
        //          += +2.00 f_aa(a,b) t1_2_aa(b,i)
        //          += +2.00 <a,j||i,b>_abab t1_2_bb(b,j)
        //          += +2.00 f_bb(j,b) t2_2_abab(a,b,i,j)
        //          += +1.00 <a,j||b,c>_abab t2_2_abab(b,c,i,j)
        //          += +1.00 <a,j||c,b>_abab t2_2_abab(c,b,i,j)
        //          += +2.00 <j,a||b,i>_aaaa t1_2_aa(b,j)
        //          += -1.00 <j,a||b,c>_aaaa t2_2_aaaa(b,c,i,j)
        //          += -1.00 <k,j||i,b>_abab t2_2_abab(a,b,k,j)
        //          += -1.00 <j,k||i,b>_abab t2_2_abab(a,b,j,k)
        //          += +4.00 t1_2_aa(a,i) w0
        //          += -4.00 d-_bb(j,b) t1_1_bb(b,j) t1_2_aa(a,i)
        //          += -4.00 d-_aa(j,b) t1_1_aa(b,j) t1_2_aa(a,i)
        //          += -1.00 <k,j||b,i>_aaaa t2_2_aaaa(b,a,k,j)
        //          += -2.00 f_aa(j,b) t2_2_aaaa(b,a,i,j)
        //          += -2.00 f_aa(j,i) t1_2_aa(a,j)
        //          += +4.00 d-_aa(j,i) t1_1_aa(a,j) t0_2
        //          += +1.00 <k,j||b,c>_aaaa t2_1_aaaa(c,a,k,j) t1_1_aa(b,i)
        //          += -2.00 <j,k||b,c>_abab t2_aaaa(b,a,i,j) t1_2_bb(c,k)
        //          += +2.00 <a,j||c,b>_abab t1_1_bb(b,j) t1_1_aa(c,i)
        //          += +2.00 d-_aa(j,i) t0_1 t1_2_aa(a,j)
        //          += -1.00 <k,j||b,c>_aaaa t2_aaaa(b,a,k,j) t1_2_aa(c,i)
        //          += +2.00 <k,j||b,c>_aaaa t2_1_aaaa(c,a,i,k) t1_1_aa(b,j)
        //          += +2.00 <j,a||b,c>_aaaa t1_1_aa(b,j) t1_1_aa(c,i)
        //          += -2.00 <k,j||b,c>_bbbb t2_abab(a,b,i,j) t1_2_bb(c,k)
        //          += -4.00 d-_aa(a,b) t1_1_aa(b,i) t0_2
        //          += -2.00 <k,j||c,b>_abab t2_1_aaaa(c,a,i,k) t1_1_bb(b,j)
        //          += +4.00 d-_aa(j,b) t2_1_aaaa(b,a,i,j) t0_2
        //          += -2.00 d-_bb(j,b) t0_1 t2_2_abab(a,b,i,j)
        //          += -2.00 <k,j||b,c>_bbbb t2_1_abab(a,c,i,k) t1_1_bb(b,j)
        //          += +2.00 <k,j||b,i>_aaaa t1_1_aa(a,k) t1_1_aa(b,j)
        //          += +2.00 d-_aa(j,b) t0_1 t2_2_aaaa(b,a,i,j)
        //          += +6.00 d-_aa(j,b) t1_1_aa(b,i) t1_2_aa(a,j)
        //          += -4.00 d-_bb(j,b) t2_1_abab(a,b,i,j) t0_2
        //          += -2.00 <k,j||i,b>_abab t1_1_aa(a,k) t1_1_bb(b,j)
        //          += -1.00 <k,j||b,c>_abab t2_abab(b,c,i,j) t1_2_aa(a,k)
        //          += -1.00 <k,j||c,b>_abab t2_abab(c,b,i,j) t1_2_aa(a,k)
        //          += -1.00 <k,j||b,c>_aaaa t2_aaaa(b,c,i,j) t1_2_aa(a,k)
        //          += +2.00 <k,j||c,b>_abab t2_abab(a,b,i,j) t1_2_aa(c,k)
        //          += -2.00 f_aa(j,b) t1_1_aa(a,j) t1_1_aa(b,i)
        rt1_2_aa("a,i")  = 6.00 * tmps_["45_aa_vo"]("a,i");
        tmps_["45_aa_vo"].~TArrayD();

        // flops: o1v1  = o2v1 o1v2 o1v1 o2v2 o2v2 o1v1 o1v1
        //  mems: o1v1  = o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1
        tmps_["46_aa_ov"]("i,a")  = -1.00 * dp["aa_oo"]("k,i") * t1_1["aa"]("a,k");
        tmps_["46_aa_ov"]("i,a") += dp["aa_vv"]("a,b") * t1_1["aa"]("b,i");
        tmps_["46_aa_ov"]("i,a") -= dp["aa_ov"]("k,b") * t2_1["aaaa"]("b,a,i,k");
        tmps_["46_aa_ov"]("i,a") += dp["bb_ov"]("j,c") * t2_1["abab"]("a,c,i,j");

        // rt1_aa  = -1.00 d-_bb(j,b) t2_1_abab(a,b,i,j)
        //        += +1.00 d-_aa(j,b) t2_1_aaaa(b,a,i,j)
        //        += -1.00 d-_aa(a,b) t1_1_aa(b,i)
        //        += +1.00 d-_aa(j,i) t1_1_aa(a,j)
        rt1_aa("a,i")  = -1.00 * tmps_["46_aa_ov"]("i,a");

        // rt1_2_aa += -2.00 d+_bb(j,b) t2_1_abab(a,b,i,j)
        //          += +2.00 d+_aa(j,b) t2_1_aaaa(b,a,i,j)
        //          += -2.00 d+_aa(a,b) t1_1_aa(b,i)
        //          += +2.00 d+_aa(j,i) t1_1_aa(a,j)
        rt1_2_aa("a,i") -= 2.00 * tmps_["46_aa_ov"]("i,a");
        tmps_["46_aa_ov"].~TArrayD();

        // flops: o2v0  = o2v0
        //  mems: o2v0  = o2v0
        tmps_["62_bb_oo"]("i,j")  = t0_2 * dp["bb_oo"]("i,j");

        // flops: o2v0  = o2v0
        //  mems: o2v0  = o2v0
        tmps_["61_bb_oo"]("i,j")  = t0_1 * dp["bb_oo"]("i,j");

        // flops: o0v2  = o0v2
        //  mems: o0v2  = o0v2
        tmps_["60_bb_vv"]("a,b")  = t0_2 * dp["bb_vv"]("a,b");

        // flops: o2v0  = o2v1
        //  mems: o2v0  = o2v0
        tmps_["59_bb_oo"]("i,j")  = f["bb_ov"]("i,a") * t1_1["bb"]("a,j");

        // flops: o2v0  = o2v1
        //  mems: o2v0  = o2v0
        tmps_["58_bb_oo"]("i,j")  = dp["bb_ov"]("i,a") * t1_2["bb"]("a,j");

        // flops: o2v0  = o3v1
        //  mems: o2v0  = o2v0
        tmps_["57_bb_oo"]("i,j")  = eri["bbbb_oovo"]("k,i,a,j") * t1_1["bb"]("a,k");

        // flops: o2v0  = o3v1
        //  mems: o2v0  = o2v0
        tmps_["56_bb_oo"]("i,j")  = eri["abab_oovo"]("k,i,a,j") * t1_1["aa"]("a,k");

        // flops: o2v0  = o3v2
        //  mems: o2v0  = o2v0
        tmps_["55_bb_oo"]("i,j")  = eri["bbbb_oovv"]("i,k,a,b") * t2_1["bbbb"]("a,b,j,k");

        // flops: o2v0  = o3v2
        //  mems: o2v0  = o2v0
        tmps_["54_bb_oo"]("i,j")  = eri["abab_oovv"]("k,i,a,b") * t2_1["abab"]("a,b,k,j");

        // flops: o2v2  = o2v3
        //  mems: o2v2  = o2v2
        tmps_["53_bbbb_vovo"]("a,i,b,j")  = eri["bbbb_vovv"]("a,i,b,c") * t1_1["bb"]("c,j");

        // flops: o0v2  = o2v3
        //  mems: o0v2  = o0v2
        tmps_["52_bb_vv"]("a,b")  = eri["bbbb_oovv"]("i,j,a,c") * t2_1["bbbb"]("c,b,j,i");

        // flops: o2v2  = o2v3
        //  mems: o2v2  = o2v2
        tmps_["51_baab_vovo"]("a,i,b,j")  = eri["baab_vovv"]("a,i,b,c") * t1_1["bb"]("c,j");

        // flops: o0v2  = o2v3
        //  mems: o0v2  = o0v2
        tmps_["50_bb_vv"]("a,b")  = eri["abab_oovv"]("i,j,c,a") * t2_1["abab"]("c,b,i,j");

        // flops: o2v2  = o3v3
        //  mems: o2v2  = o2v2
        tmps_["49_bbbb_ovvo"]("i,a,b,j")  = eri["bbbb_oovv"]("i,k,a,c") * t2_1["bbbb"]("c,b,j,k");

        // flops: o2v2  = o3v3
        //  mems: o2v2  = o2v2
        tmps_["48_bbbb_ovvo"]("i,a,b,j")  = eri["abab_oovv"]("k,i,c,a") * t2_1["abab"]("c,b,k,j");

        // flops: o2v2  = o3v3
        //  mems: o2v2  = o2v2
        tmps_["47_aabb_ovvo"]("i,a,b,j")  = eri["aaaa_oovv"]("i,k,a,c") * t2_1["abab"]("c,b,k,j");

        // flops: o1v1  = o1v2 o2v1 o1v1 o1v1 o1v2 o1v1 o3v2 o1v1 o2v2 o1v1 o2v3 o1v1 o2v1 o1v1 o2v2 o1v1 o3v2 o1v1 o2v2 o1v1 o2v2 o1v1 o1v1 o2v3 o1v1 o1v1 o1v1 o1v1 o1v1 o1v2 o2v2 o2v2 o1v1 o1v2 o1v1 o2v2 o1v1 o2v1 o1v1 o2v2 o1v1 o1v1 o2v1 o1v1 o2v2 o1v1 o2v1 o1v1 o2v2 o1v1 o1v2 o1v1 o2v2 o1v1 o1v2 o1v1 o1v2 o1v1 o2v2 o1v1 o2v2 o2v2 o1v1 o1v1 o2v1 o1v1 o2v1 o1v1 o2v2 o1v1 o2v2 o1v1 o2v1 o1v1 o2v1 o1v1 o2v1 o1v1 o2v1 o1v1 o2v2 o1v1 o2v1 o1v1 o2v2 o1v1 o2v1 o1v1
        //  mems: o1v1  = o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o2v2 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1
        tmps_["63_bb_vo"]("a,i")  = -0.50 * t1_1["bb"]("c,i") * tmps_["52_bb_vv"]("c,a");
        tmps_["63_bb_vo"]("a,i") -= 0.50 * tmps_["55_bb_oo"]("j,i") * t1_1["bb"]("a,j");
        tmps_["63_bb_vo"]("a,i") -= 2.00 * scalars_["1"] * t1_2["bb"]("a,i");
        tmps_["63_bb_vo"]("a,i") += f["bb_vv"]("a,c") * t1_2["bb"]("c,i");
        tmps_["63_bb_vo"]("a,i") += 0.50 * eri["bbbb_oovo"]("j,l,c,i") * t2_2["bbbb"]("c,a,l,j");
        tmps_["63_bb_vo"]("a,i") -= eri["bbbb_vovo"]("a,j,c,i") * t1_2["bb"]("c,j");
        tmps_["63_bb_vo"]("a,i") += 0.50 * eri["bbbb_vovv"]("a,j,c,d") * t2_2["bbbb"]("c,d,i,j");
        tmps_["63_bb_vo"]("a,i") -= f["bb_oo"]("j,i") * t1_2["bb"]("a,j");
        tmps_["63_bb_vo"]("a,i") -= f["bb_ov"]("j,c") * t2_2["bbbb"]("c,a,i,j");
        tmps_["63_bb_vo"]("a,i") -= eri["abab_oovo"]("m,j,b,i") * t2_2["abab"]("b,a,m,j");
        tmps_["63_bb_vo"]("a,i") -= eri["baab_vovo"]("a,k,b,i") * t1_2["aa"]("b,k");
        tmps_["63_bb_vo"]("a,i") += f["aa_ov"]("k,b") * t2_2["abab"]("b,a,k,i");
        tmps_["63_bb_vo"]("a,i") -= scalars_["6"] * t1_1["bb"]("a,i");
        tmps_["63_bb_vo"]("a,i") -= eri["baab_vovv"]("a,k,b,d") * t2_2["abab"]("b,d,k,i");
        tmps_["63_bb_vo"]("a,i") += 2.00 * w0 * t1_2["bb"]("a,i");
        tmps_["63_bb_vo"]("a,i") -= t1_1["bb"]("c,i") * tmps_["50_bb_vv"]("c,a");
        tmps_["63_bb_vo"]("a,i") += t1_1["bb"]("c,j") * (-1.00 * tmps_["53_bbbb_vovo"]("a,j,c,i") + tmps_["48_bbbb_ovvo"]("j,c,a,i"));
        tmps_["63_bb_vo"]("a,i") -= t1_2["bb"]("c,i") * tmps_["25_bb_vv"]("a,c");
        tmps_["63_bb_vo"]("a,i") += t1_2["bb"]("d,l") * tmps_["18_bbbb_ovvo"]("l,d,a,i");
        tmps_["63_bb_vo"]("a,i") += tmps_["24_bb_oo"]("j,i") * tmps_["26_bb_vo"]("a,j");
        tmps_["63_bb_vo"]("a,i") += t1_1["aa"]("b,k") * tmps_["47_aabb_ovvo"]("k,b,a,i");
        tmps_["63_bb_vo"]("a,i") += t1_2["bb"]("a,j") * tmps_["61_bb_oo"]("j,i");
        tmps_["63_bb_vo"]("a,i") -= t1_1["bb"]("c,j") * tmps_["49_bbbb_ovvo"]("j,c,a,i");
        tmps_["63_bb_vo"]("a,i") += 2.00 * tmps_["62_bb_oo"]("j,i") * t1_1["bb"]("a,j");
        tmps_["63_bb_vo"]("a,i") += t1_2["aa"]("e,m") * tmps_["17_aabb_ovvo"]("m,e,a,i");
        tmps_["63_bb_vo"]("a,i") -= 2.00 * t1_1["bb"]("c,i") * tmps_["60_bb_vv"]("a,c");
        tmps_["63_bb_vo"]("a,i") -= t1_2["bb"]("d,l") * tmps_["19_bbbb_ovvo"]("l,d,a,i");
        tmps_["63_bb_vo"]("a,i") += 0.50 * t1_2["bb"]("d,i") * tmps_["21_bb_vv"]("d,a");
        tmps_["63_bb_vo"]("a,i") -= t1_2["bb"]("d,i") * tmps_["20_bb_vv"]("d,a");
        tmps_["63_bb_vo"]("a,i") -= t1_1["aa"]("b,k") * tmps_["51_baab_vovo"]("a,k,b,i");
        tmps_["63_bb_vo"]("a,i") -= 2.00 * t2_1["abab"]("b,a,k,i") * tmps_["13_aa_ov"]("k,b");
        tmps_["63_bb_vo"]("a,i") += tmps_["14_bb_ov"]("j,c") * t2_2["bbbb"]("c,a,i,j");
        tmps_["63_bb_vo"]("a,i") += 3.00 * tmps_["24_bb_oo"]("j,i") * t1_2["bb"]("a,j");
        tmps_["63_bb_vo"]("a,i") -= tmps_["59_bb_oo"]("j,i") * t1_1["bb"]("a,j");
        tmps_["63_bb_vo"]("a,i") += 2.00 * t2_1["bbbb"]("c,a,i,j") * tmps_["15_bb_ov"]("j,c");
        tmps_["63_bb_vo"]("a,i") -= t2["bbbb"]("c,a,i,j") * tmps_["35_bb_ov"]("j,c");
        tmps_["63_bb_vo"]("a,i") -= tmps_["57_bb_oo"]("l,i") * t1_1["bb"]("a,l");
        tmps_["63_bb_vo"]("a,i") -= tmps_["22_bb_oo"]("l,i") * t1_2["bb"]("a,l");
        tmps_["63_bb_vo"]("a,i") += 0.50 * tmps_["23_bb_oo"]("l,i") * t1_2["bb"]("a,l");
        tmps_["63_bb_vo"]("a,i") -= tmps_["56_bb_oo"]("l,i") * t1_1["bb"]("a,l");
        tmps_["63_bb_vo"]("a,i") -= t2_1["bbbb"]("d,a,i,l") * tmps_["8_bb_ov"]("l,d");
        tmps_["63_bb_vo"]("a,i") -= tmps_["54_bb_oo"]("j,i") * t1_1["bb"]("a,j");
        tmps_["63_bb_vo"]("a,i") -= t2_2["abab"]("b,a,k,i") * tmps_["11_aa_ov"]("k,b");
        tmps_["63_bb_vo"]("a,i") += 3.00 * tmps_["58_bb_oo"]("j,i") * t1_1["bb"]("a,j");
        tmps_["35_bb_ov"].~TArrayD();
        tmps_["15_bb_ov"].~TArrayD();
        tmps_["13_aa_ov"].~TArrayD();

        // rt1_2_bb  = +2.00 d-_bb(j,b) t0_1 t2_2_bbbb(b,a,i,j)
        //          += -4.00 d-_aa(j,b) t2_1_abab(b,a,j,i) t0_2
        //          += +2.00 <k,j||c,b>_abab t2_1_abab(c,a,k,i) t1_1_bb(b,j)
        //          += +2.00 <j,a||b,c>_bbbb t1_1_bb(b,j) t1_1_bb(c,i)
        //          += -1.00 <k,j||c,b>_abab t2_1_abab(c,a,k,j) t1_1_bb(b,i)
        //          += -1.00 <j,k||c,b>_abab t2_1_abab(c,a,j,k) t1_1_bb(b,i)
        //          += +2.00 <j,k||b,c>_abab t2_abab(b,a,j,i) t1_2_bb(c,k)
        //          += +2.00 d-_bb(j,b) t0_1 t1_1_bb(a,j) t1_1_bb(b,i)
        //          += -2.00 <k,j||b,c>_aaaa t2_1_abab(c,a,k,i) t1_1_aa(b,j)
        //          += -2.00 d-_bb(a,b) t0_1 t1_2_bb(b,i)
        //          += +1.00 <j,a||b,c>_abab t2_2_abab(b,c,j,i)
        //          += +1.00 <j,a||c,b>_abab t2_2_abab(c,b,j,i)
        //          += -2.00 d-_aa(j,b) t1_1_bb(a,i) t1_2_aa(b,j)
        //          += -2.00 d-_bb(j,b) t1_1_bb(a,i) t1_2_bb(b,j)
        //          += -4.00 d-_bb(j,b) t1_1_bb(b,j) t1_2_bb(a,i)
        //          += -4.00 d-_aa(j,b) t1_1_aa(b,j) t1_2_bb(a,i)
        //          += +2.00 f_bb(a,b) t1_2_bb(b,i)
        //          += -1.00 <k,j||b,i>_bbbb t2_2_bbbb(b,a,k,j)
        //          += +2.00 <j,a||b,i>_bbbb t1_2_bb(b,j)
        //          += -1.00 <j,a||b,c>_bbbb t2_2_bbbb(b,c,i,j)
        //          += -2.00 f_bb(j,i) t1_2_bb(a,j)
        //          += -2.00 f_bb(j,b) t2_2_bbbb(b,a,i,j)
        //          += -1.00 <k,j||b,i>_abab t2_2_abab(b,a,k,j)
        //          += -1.00 <j,k||b,i>_abab t2_2_abab(b,a,j,k)
        //          += +2.00 <j,a||b,i>_abab t1_2_aa(b,j)
        //          += +2.00 f_aa(j,b) t2_2_abab(b,a,j,i)
        //          += +4.00 t1_2_bb(a,i) w0
        //          += +1.00 <k,j||b,c>_bbbb t2_1_bbbb(c,a,k,j) t1_1_bb(b,i)
        //          += +2.00 d-_bb(j,i) t0_1 t1_2_bb(a,j)
        //          += +2.00 <k,j||b,c>_bbbb t2_1_bbbb(c,a,i,k) t1_1_bb(b,j)
        //          += +4.00 d-_bb(j,i) t1_1_bb(a,j) t0_2
        //          += -2.00 <k,j||b,c>_aaaa t2_abab(b,a,j,i) t1_2_aa(c,k)
        //          += -4.00 d-_bb(a,b) t1_1_bb(b,i) t0_2
        //          += +2.00 <k,j||b,c>_bbbb t2_bbbb(b,a,i,j) t1_2_bb(c,k)
        //          += -1.00 <k,j||b,c>_bbbb t2_bbbb(b,a,k,j) t1_2_bb(c,i)
        //          += -1.00 <k,j||b,c>_abab t2_abab(b,a,k,j) t1_2_bb(c,i)
        //          += -1.00 <j,k||b,c>_abab t2_abab(b,a,j,k) t1_2_bb(c,i)
        //          += +2.00 <j,a||b,c>_abab t1_1_aa(b,j) t1_1_bb(c,i)
        //          += +6.00 d-_bb(j,b) t1_1_bb(b,i) t1_2_bb(a,j)
        //          += -2.00 f_bb(j,b) t1_1_bb(a,j) t1_1_bb(b,i)
        //          += +4.00 d-_bb(j,b) t2_1_bbbb(b,a,i,j) t0_2
        //          += +1.00 <k,j||b,c>_bbbb t1_1_bb(a,j) t2_1_bbbb(b,c,i,k)
        //          += -2.00 <k,j||c,b>_abab t2_bbbb(b,a,i,j) t1_2_aa(c,k)
        //          += +2.00 <k,j||b,i>_bbbb t1_1_bb(a,k) t1_1_bb(b,j)
        //          += -1.00 <j,k||b,c>_abab t2_abab(b,c,j,i) t1_2_bb(a,k)
        //          += -1.00 <j,k||c,b>_abab t2_abab(c,b,j,i) t1_2_bb(a,k)
        //          += -1.00 <k,j||b,c>_bbbb t2_bbbb(b,c,i,j) t1_2_bb(a,k)
        //          += -2.00 <j,k||b,i>_abab t1_1_bb(a,k) t1_1_aa(b,j)
        //          += -2.00 <j,k||b,c>_abab t2_1_bbbb(c,a,i,k) t1_1_aa(b,j)
        //          += -1.00 <k,j||b,c>_abab t1_1_bb(a,j) t2_1_abab(b,c,k,i)
        //          += -1.00 <k,j||c,b>_abab t1_1_bb(a,j) t2_1_abab(c,b,k,i)
        //          += -2.00 d-_aa(j,b) t0_1 t2_2_abab(b,a,j,i)
        //          += +6.00 d-_bb(j,b) t1_1_bb(a,j) t1_2_bb(b,i)
        rt1_2_bb("a,i")  = 2.00 * tmps_["63_bb_vo"]("a,i");
        tmps_["63_bb_vo"].~TArrayD();

        // flops: o1v1  = o2v2 o2v2 o1v1 o1v2 o2v1 o1v1 o1v1
        //  mems: o1v1  = o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1
        tmps_["64_bb_vo"]("a,i")  = -1.00 * dp["aa_ov"]("k,c") * t2_1["abab"]("c,a,k,i");
        tmps_["64_bb_vo"]("a,i") += dp["bb_ov"]("j,b") * t2_1["bbbb"]("b,a,i,j");
        tmps_["64_bb_vo"]("a,i") -= dp["bb_vv"]("a,b") * t1_1["bb"]("b,i");
        tmps_["64_bb_vo"]("a,i") += dp["bb_oo"]("j,i") * t1_1["bb"]("a,j");

        // rt1_bb  = +1.00 d-_bb(j,i) t1_1_bb(a,j)
        //        += -1.00 d-_bb(a,b) t1_1_bb(b,i)
        //        += +1.00 d-_bb(j,b) t2_1_bbbb(b,a,i,j)
        //        += -1.00 d-_aa(j,b) t2_1_abab(b,a,j,i)
        rt1_bb("a,i")  = tmps_["64_bb_vo"]("a,i");

        // rt1_2_bb += +2.00 d+_bb(j,i) t1_1_bb(a,j)
        //          += -2.00 d+_bb(a,b) t1_1_bb(b,i)
        //          += +2.00 d+_bb(j,b) t2_1_bbbb(b,a,i,j)
        //          += -2.00 d+_aa(j,b) t2_1_abab(b,a,j,i)
        rt1_2_bb("a,i") += 2.00 * tmps_["64_bb_vo"]("a,i");
        tmps_["64_bb_vo"].~TArrayD();

        // flops: o1v1  = o2v2 o2v3 o1v1 o3v2 o1v1 o1v1 o1v1 o1v1 o2v2 o1v1 o3v2 o1v1 o2v3 o1v1 o1v1 o1v1 o2v2 o1v1 o2v2 o1v1
        //  mems: o1v1  = o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1
        tmps_["65_aa_vo"]("a,i")  = f["aa_ov"]("k,c") * t2["aaaa"]("c,a,i,k");
        tmps_["65_aa_vo"]("a,i") -= 0.50 * eri["aaaa_vovv"]("a,k,c,e") * t2["aaaa"]("c,e,i,k");
        tmps_["65_aa_vo"]("a,i") -= eri["abba_oovo"]("l,j,b,i") * t2["abab"]("a,b,l,j");
        tmps_["65_aa_vo"]("a,i") -= f["aa_vo"]("a,i");
        tmps_["65_aa_vo"]("a,i") += t0_1 * dp["aa_vo"]("a,i");
        tmps_["65_aa_vo"]("a,i") -= f["bb_ov"]("j,b") * t2["abab"]("a,b,i,j");
        tmps_["65_aa_vo"]("a,i") -= 0.50 * eri["aaaa_oovo"]("k,l,c,i") * t2["aaaa"]("c,a,l,k");
        tmps_["65_aa_vo"]("a,i") -= eri["abab_vovv"]("a,j,c,d") * t2["abab"]("c,d,i,j");
        tmps_["65_aa_vo"]("a,i") += scalars_["5"] * t1_1["aa"]("a,i");
        tmps_["65_aa_vo"]("a,i") += tmps_["14_bb_ov"]("j,b") * t2["abab"]("a,b,i,j");
        tmps_["65_aa_vo"]("a,i") -= t2["aaaa"]("c,a,i,k") * tmps_["11_aa_ov"]("k,c");

        // rt1_aa += -1.00 d-_bb(j,b) t2_abab(a,b,i,j) t0_1
        //        += -0.50 <j,a||b,c>_aaaa t2_aaaa(b,c,i,j)
        //        += -1.00 f_aa(j,b) t2_aaaa(b,a,i,j)
        //        += -0.50 <k,j||i,b>_abab t2_abab(a,b,k,j)
        //        += -0.50 <j,k||i,b>_abab t2_abab(a,b,j,k)
        //        += +1.00 f_aa(a,i)
        //        += -1.00 d-_aa(a,i) t0_1
        //        += +1.00 f_bb(j,b) t2_abab(a,b,i,j)
        //        += -0.50 <k,j||b,i>_aaaa t2_aaaa(b,a,k,j)
        //        += +0.50 <a,j||b,c>_abab t2_abab(b,c,i,j)
        //        += +0.50 <a,j||c,b>_abab t2_abab(c,b,i,j)
        //        += -1.00 d-_aa(j,j) t1_1_aa(a,i)
        //        += -1.00 d-_bb(j,j) t1_1_aa(a,i)
        //        += +1.00 d-_aa(j,b) t2_aaaa(b,a,i,j) t0_1
        rt1_aa("a,i") -= tmps_["65_aa_vo"]("a,i");
        tmps_["65_aa_vo"].~TArrayD();

        // flops: o1v1  = o1v1 o3v2 o1v1 o2v3 o1v1 o2v2 o1v1 o2v2 o1v1 o2v3 o1v1 o3v2 o1v1 o1v1 o1v1 o1v1 o2v2 o1v1 o2v2 o1v1
        //  mems: o1v1  = o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1
        tmps_["66_bb_vo"]("a,i")  = -1.00 * scalars_["5"] * t1_1["bb"]("a,i");
        tmps_["66_bb_vo"]("a,i") -= eri["abab_oovo"]("m,j,c,i") * t2["abab"]("c,a,m,j");
        tmps_["66_bb_vo"]("a,i") -= eri["baab_vovv"]("a,k,c,d") * t2["abab"]("c,d,k,i");
        tmps_["66_bb_vo"]("a,i") -= f["bb_ov"]("j,b") * t2["bbbb"]("b,a,i,j");
        tmps_["66_bb_vo"]("a,i") += f["aa_ov"]("k,c") * t2["abab"]("c,a,k,i");
        tmps_["66_bb_vo"]("a,i") += 0.50 * eri["bbbb_vovv"]("a,j,b,d") * t2["bbbb"]("b,d,i,j");
        tmps_["66_bb_vo"]("a,i") += 0.50 * eri["bbbb_oovo"]("j,l,b,i") * t2["bbbb"]("b,a,l,j");
        tmps_["66_bb_vo"]("a,i") += f["bb_vo"]("a,i");
        tmps_["66_bb_vo"]("a,i") += -1.00 * t0_1 * dp["bb_vo"]("a,i");
        tmps_["66_bb_vo"]("a,i") -= t2["abab"]("c,a,k,i") * tmps_["11_aa_ov"]("k,c");
        tmps_["66_bb_vo"]("a,i") += tmps_["14_bb_ov"]("j,b") * t2["bbbb"]("b,a,i,j");
        tmps_["14_bb_ov"].~TArrayD();
        tmps_["11_aa_ov"].~TArrayD();

        // rt1_bb += +1.00 d-_bb(j,b) t2_bbbb(b,a,i,j) t0_1
        //        += +1.00 f_aa(j,b) t2_abab(b,a,j,i)
        //        += -1.00 d-_aa(j,j) t1_1_bb(a,i)
        //        += -1.00 d-_bb(j,j) t1_1_bb(a,i)
        //        += -1.00 d-_bb(a,i) t0_1
        //        += +0.50 <j,a||b,c>_abab t2_abab(b,c,j,i)
        //        += +0.50 <j,a||c,b>_abab t2_abab(c,b,j,i)
        //        += -1.00 f_bb(j,b) t2_bbbb(b,a,i,j)
        //        += -0.50 <j,a||b,c>_bbbb t2_bbbb(b,c,i,j)
        //        += -0.50 <k,j||b,i>_abab t2_abab(b,a,k,j)
        //        += -0.50 <j,k||b,i>_abab t2_abab(b,a,j,k)
        //        += -0.50 <k,j||b,i>_bbbb t2_bbbb(b,a,k,j)
        //        += +1.00 f_bb(a,i)
        //        += -1.00 d-_aa(j,b) t2_abab(b,a,j,i) t0_1
        rt1_bb("a,i") += tmps_["66_bb_vo"]("a,i");
        tmps_["66_bb_vo"].~TArrayD();

        // flops: o0v2  = o1v3
        //  mems: o0v2  = o0v2
        tmps_["72_aa_vv"]("a,b")  = eri["abab_vovv"]("a,i,b,c") * t1_1["bb"]("c,i");

        // flops: o0v2  = o1v3
        //  mems: o0v2  = o0v2
        tmps_["71_aa_vv"]("a,b")  = eri["aaaa_vovv"]("a,i,b,c") * t1_1["aa"]("c,i");

        // flops: o3v1  = o3v2
        //  mems: o3v1  = o3v1
        tmps_["70_aaaa_ovoo"]("i,a,j,k")  = f["aa_ov"]("i,b") * t2["aaaa"]("b,a,j,k");

        // flops: o3v1  = o3v2
        //  mems: o3v1  = o3v1
        tmps_["69_aaaa_ovoo"]("i,a,j,k")  = dp["aa_ov"]("i,b") * t2_1["aaaa"]("b,a,j,k");

        // flops: o3v1  = o3v2
        //  mems: o3v1  = o3v1
        tmps_["68_aaaa_ovoo"]("i,a,j,k")  = dp["aa_ov"]("i,b") * t2["aaaa"]("b,a,j,k");

        // flops: o3v1  = o3v3
        //  mems: o3v1  = o3v1
        tmps_["67_aaaa_vooo"]("a,i,j,k")  = eri["aaaa_vovv"]("a,i,b,c") * t2["aaaa"]("b,c,j,k");

        // flops: o2v2  = o2v3 o2v3 o2v2 o2v3 o2v2 o3v2 o2v2 o3v2 o2v2 o2v3 o2v2 o2v3 o2v2 o3v2 o2v3 o2v2 o3v2 o2v2 o2v3 o2v2 o2v3 o2v2 o3v2 o2v2 o2v2 o3v2 o2v2 o2v3 o2v2 o2v3 o2v2 o2v3 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["73_aaaa_vvoo"]("a,b,i,j")  = -0.50 * dp["aa_vv"]("a,c") * t2["aaaa"]("c,b,i,j");
        tmps_["73_aaaa_vvoo"]("a,b,i,j") += 0.50 * f["aa_vv"]("a,c") * t2_1["aaaa"]("c,b,i,j");
        tmps_["73_aaaa_vvoo"]("a,b,i,j") -= dp["aa_vv"]("a,c") * t2_2["aaaa"]("c,b,i,j");
        tmps_["73_aaaa_vvoo"]("a,b,i,j") -= 0.50 * eri["aaaa_vooo"]("a,k,i,j") * t1_1["aa"]("b,k");
        tmps_["73_aaaa_vvoo"]("a,b,i,j") -= 0.50 * tmps_["12_aa_vo"]("b,k") * tmps_["68_aaaa_ovoo"]("k,a,i,j");
        tmps_["73_aaaa_vvoo"]("a,b,i,j") -= 0.50 * t2_1["aaaa"]("c,b,i,j") * tmps_["10_aa_vv"]("a,c");
        tmps_["73_aaaa_vvoo"]("a,b,i,j") -= t2["aaaa"]("c,b,i,j") * tmps_["40_aa_vv"]("a,c");
        tmps_["73_aaaa_vvoo"]("a,b,i,j") -= 0.25 * t1_1["aa"]("b,k") * tmps_["67_aaaa_vooo"]("a,k,i,j");
        tmps_["73_aaaa_vvoo"]("a,b,i,j") -= 0.50 * t2_1["aaaa"]("e,b,i,j") * tmps_["5_aa_vv"]("e,a");
        tmps_["73_aaaa_vvoo"]("a,b,i,j") += tmps_["69_aaaa_ovoo"]("k,b,i,j") * t1_1["aa"]("a,k");
        tmps_["73_aaaa_vvoo"]("a,b,i,j") += 0.25 * tmps_["29_aa_vv"]("c,b") * t2["aaaa"]("c,a,i,j");
        tmps_["73_aaaa_vvoo"]("a,b,i,j") += 0.25 * t2_1["aaaa"]("e,b,i,j") * tmps_["4_aa_vv"]("e,a");
        tmps_["73_aaaa_vvoo"]("a,b,i,j") += 0.50 * t1_1["aa"]("b,k") * tmps_["70_aaaa_ovoo"]("k,a,i,j");
        tmps_["73_aaaa_vvoo"]("a,b,i,j") -= t1_2["aa"]("b,k") * tmps_["68_aaaa_ovoo"]("k,a,i,j");
        tmps_["73_aaaa_vvoo"]("a,b,i,j") += 0.50 * tmps_["32_aa_vv"]("c,b") * t2["aaaa"]("c,a,i,j");
        tmps_["73_aaaa_vvoo"]("a,b,i,j") += 0.50 * t2["aaaa"]("c,b,i,j") * tmps_["71_aa_vv"]("a,c");
        tmps_["73_aaaa_vvoo"]("a,b,i,j") += 0.50 * t2["aaaa"]("c,b,i,j") * tmps_["72_aa_vv"]("a,c");

        // rt2_1_aaaa  = +2.00 P(a,b) d-_aa(k,c) t1_1_aa(a,k) t2_1_aaaa(c,b,i,j)
        //            += -0.50 P(a,b) <l,k||d,c>_abab t2_abab(a,c,l,k) t2_1_aaaa(d,b,i,j)
        //            += -0.50 P(a,b) <k,l||d,c>_abab t2_abab(a,c,k,l) t2_1_aaaa(d,b,i,j)
        //            += +0.50 P(a,b) <k,a||c,d>_aaaa t2_aaaa(c,d,i,j) t1_1_aa(b,k)
        //            += -0.50 P(a,b) <l,k||c,d>_aaaa t2_aaaa(c,a,i,j) t2_1_aaaa(d,b,l,k)
        //            += -0.50 P(a,b) <l,k||c,d>_aaaa t2_aaaa(c,a,l,k) t2_1_aaaa(d,b,i,j)
        //            += +1.00 P(a,b) f_aa(k,c) t2_aaaa(c,a,i,j) t1_1_aa(b,k)
        //            += +1.00 P(a,b) f_aa(a,c) t2_1_aaaa(c,b,i,j)
        //            += -1.00 P(a,b) d+_aa(a,c) t2_aaaa(c,b,i,j)
        //            += -2.00 P(a,b) d-_aa(a,c) t2_2_aaaa(c,b,i,j)
        //            += +1.00 P(a,b) <k,a||i,j>_aaaa t1_1_aa(b,k)
        //            += -1.00 P(a,b) d-_aa(k,c) t2_aaaa(c,a,i,j) t0_1 t1_1_aa(b,k)
        //            += -1.00 P(a,b) d-_aa(a,c) t0_1 t2_1_aaaa(c,b,i,j)
        //            += -2.00 P(a,b) d-_aa(a,c) t2_aaaa(c,b,i,j) t0_2
        //            += +0.50 P(a,b) <l,k||c,d>_abab t2_aaaa(c,a,i,j) t2_1_abab(b,d,l,k)
        //            += +0.50 P(a,b) <k,l||c,d>_abab t2_aaaa(c,a,i,j) t2_1_abab(b,d,k,l)
        //            += -1.00 P(a,b) <k,a||c,d>_aaaa t2_aaaa(c,b,i,j) t1_1_aa(d,k)
        //            += +1.00 P(a,b) <a,k||c,d>_abab t2_aaaa(c,b,i,j) t1_1_bb(d,k)
        //            += -2.00 P(a,b) d-_aa(k,c) t2_aaaa(c,a,i,j) t1_2_aa(b,k)
        rt2_1_aaaa("a,b,i,j")  = 2.00 * tmps_["73_aaaa_vvoo"]("a,b,i,j");
        rt2_1_aaaa("a,b,i,j") -= 2.00 * tmps_["73_aaaa_vvoo"]("b,a,i,j");
        tmps_["73_aaaa_vvoo"].~TArrayD();

        // flops: o3v1  = o4v2 o4v2 o3v1
        //  mems: o3v1  = o3v1 o3v1 o3v1
        tmps_["77_aaaa_oovo"]("i,j,a,k")  = eri["aaaa_oovo"]("l,i,b,j") * t2["aaaa"]("b,a,k,l");
        tmps_["77_aaaa_oovo"]("i,j,a,k") += eri["abba_oovo"]("i,m,c,j") * t2["abab"]("a,c,k,m");

        // flops: o1v1  = o2v1 o1v2 o1v1
        //  mems: o1v1  = o1v1 o1v1 o1v1
        tmps_["76_aa_ov"]("i,a")  = -1.00 * dp["aa_oo"]("j,i") * t1_1["aa"]("a,j");
        tmps_["76_aa_ov"]("i,a") += dp["aa_vv"]("a,b") * t1_1["aa"]("b,i");

        // flops: o1v1  = o2v2 o2v2 o1v1
        //  mems: o1v1  = o1v1 o1v1 o1v1
        tmps_["75_aa_vo"]("a,i")  = -1.00 * dp["aa_ov"]("k,c") * t2["aaaa"]("c,a,i,k");
        tmps_["75_aa_vo"]("a,i") += dp["bb_ov"]("j,b") * t2["abab"]("a,b,i,j");

        // flops: o1v1  = o2v2 o2v2 o1v1
        //  mems: o1v1  = o1v1 o1v1 o1v1
        tmps_["74_aa_vo"]("a,i")  = -1.00 * dp["aa_ov"]("k,c") * t2_1["aaaa"]("c,a,i,k");
        tmps_["74_aa_vo"]("a,i") += dp["bb_ov"]("j,b") * t2_1["abab"]("a,b,i,j");

        // flops: o2v2  = o3v3 o3v3 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o3v2 o2v2 o3v3 o3v3 o2v2 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["78_aaaa_vovo"]("a,i,b,j")  = eri["aaaa_vovo"]("a,l,d,i") * t2_1["aaaa"]("d,b,j,l");
        tmps_["78_aaaa_vovo"]("a,i,b,j") -= eri["abba_vovo"]("a,k,c,i") * t2_1["abab"]("b,c,j,k");
        tmps_["78_aaaa_vovo"]("a,i,b,j") -= 2.00 * dp["aa_vo"]("a,i") * t1_2["aa"]("b,j");
        tmps_["78_aaaa_vovo"]("a,i,b,j") -= t1_1["aa"]("a,i") * tmps_["74_aa_vo"]("b,j");
        tmps_["78_aaaa_vovo"]("a,i,b,j") -= t1_1["aa"]("b,j") * tmps_["76_aa_ov"]("i,a");
        tmps_["78_aaaa_vovo"]("a,i,b,j") -= 2.00 * t1_2["aa"]("b,j") * tmps_["75_aa_vo"]("a,i");
        tmps_["78_aaaa_vovo"]("a,i,b,j") -= t1_1["aa"]("b,n") * tmps_["77_aaaa_oovo"]("n,i,a,j");
        tmps_["78_aaaa_vovo"]("a,i,b,j") += tmps_["1_aaaa_ovvo"]("n,e,a,i") * t2_1["aaaa"]("e,b,j,n");
        tmps_["78_aaaa_vovo"]("a,i,b,j") += tmps_["3_bbaa_ovvo"]("m,f,a,i") * t2_1["abab"]("b,f,j,m");
        tmps_["78_aaaa_vovo"]("a,i,b,j") -= tmps_["28_bbaa_ovvo"]("k,c,b,j") * t2["abab"]("a,c,i,k");
        tmps_["78_aaaa_vovo"]("a,i,b,j") -= t2_1["abab"]("b,f,j,m") * tmps_["2_bbaa_ovvo"]("m,f,a,i");
        tmps_["78_aaaa_vovo"]("a,i,b,j") -= t2["aaaa"]("d,b,i,l") * tmps_["30_aaaa_vovo"]("a,l,d,j");
        tmps_["78_aaaa_vovo"]("a,i,b,j") -= t2["abab"]("b,c,i,k") * tmps_["31_abba_vovo"]("a,k,c,j");

        // rt2_1_aaaa += +1.00 P(i,j) P(a,b) <l,k||c,d>_bbbb t2_abab(a,c,j,k) t2_1_abab(b,d,i,l)
        //            += +1.00 P(i,j) P(a,b) <l,k||c,d>_aaaa t2_aaaa(c,a,j,k) t2_1_aaaa(d,b,i,l)
        //            += +1.00 P(i,j) P(a,b) d-_bb(k,c) t1_1_aa(a,j) t2_1_abab(b,c,i,k)
        //            += -1.00 P(i,j) P(a,b) d-_aa(k,c) t1_1_aa(a,j) t2_1_aaaa(c,b,i,k)
        //            += -1.00 P(i,j) P(a,b) <a,k||j,c>_abab t2_1_abab(b,c,i,k)
        //            += +1.00 P(i,j) P(a,b) <k,a||c,j>_aaaa t2_1_aaaa(c,b,i,k)
        //            += +2.00 P(i,j) P(a,b) d-_aa(a,j) t1_2_aa(b,i)
        //            += +1.00 P(i,j) P(a,b) d-_aa(a,c) t1_1_aa(b,i) t1_1_aa(c,j)
        //            += -1.00 P(i,j) P(a,b) d-_aa(k,j) t1_1_aa(a,k) t1_1_aa(b,i)
        //            += +2.00 P(i,j) P(a,b) d-_bb(k,c) t2_abab(a,c,j,k) t1_2_aa(b,i)
        //            += -2.00 P(i,j) P(a,b) d-_aa(k,c) t2_aaaa(c,a,j,k) t1_2_aa(b,i)
        //            += -1.00 P(i,j) P(a,b) <l,k||c,j>_aaaa t2_aaaa(c,a,i,k) t1_1_aa(b,l)
        //            += -1.00 P(i,j) P(a,b) <l,k||j,c>_abab t2_abab(a,c,i,k) t1_1_aa(b,l)
        //            += +1.00 P(i,j) P(a,b) <l,k||d,c>_abab t2_abab(a,c,j,k) t2_1_aaaa(d,b,i,l)
        //            += +1.00 P(i,j) P(a,b) <k,l||c,d>_abab t2_aaaa(c,a,j,k) t2_1_abab(b,d,i,l)
        //            += -1.00 P(i,j) P(a,b) <k,a||c,d>_aaaa t2_aaaa(c,b,j,k) t1_1_aa(d,i)
        //            += +1.00 P(i,j) P(a,b) <a,k||d,c>_abab t2_abab(b,c,j,k) t1_1_aa(d,i)
        rt2_1_aaaa("a,b,i,j") -= tmps_["78_aaaa_vovo"]("a,j,b,i");
        rt2_1_aaaa("a,b,i,j") += tmps_["78_aaaa_vovo"]("a,i,b,j");
        rt2_1_aaaa("a,b,i,j") += tmps_["78_aaaa_vovo"]("b,j,a,i");
        rt2_1_aaaa("a,b,i,j") -= tmps_["78_aaaa_vovo"]("b,i,a,j");
        tmps_["78_aaaa_vovo"].~TArrayD();

        // flops: o2v0  = o2v0
        //  mems: o2v0  = o2v0
        tmps_["81_aa_oo"]("i,j")  = t0_1 * tmps_["9_aa_oo"]("i,j");

        // flops: o2v0  = o3v1
        //  mems: o2v0  = o2v0
        tmps_["80_aa_oo"]("i,j")  = eri["aaaa_oovo"]("i,k,a,j") * t1_1["aa"]("a,k");

        // flops: o4v0  = o4v1
        //  mems: o4v0  = o4v0
        tmps_["79_aaaa_oooo"]("i,j,k,l")  = eri["aaaa_oovo"]("i,j,a,k") * t1_1["aa"]("a,l");

        // flops: o2v2  = o4v2 o3v2 o2v3 o2v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o3v2 o2v2 o2v2 o3v2 o2v2 o2v2 o3v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o2v2 o3v2 o2v2 o3v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["82_aaaa_ovvo"]("i,a,b,j")  = -0.25 * tmps_["79_aaaa_oooo"]("k,l,i,j") * t2["aaaa"]("a,b,l,k");
        tmps_["82_aaaa_ovvo"]("i,a,b,j") += 0.50 * (dp["aa_oo"]("k,i") * t2["aaaa"]("a,b,j,k") + eri["aaaa_vvvo"]("a,b,c,i") * t1_1["aa"]("c,j"));
        tmps_["82_aaaa_ovvo"]("i,a,b,j") -= 0.50 * f["aa_oo"]("k,i") * t2_1["aaaa"]("a,b,j,k");
        tmps_["82_aaaa_ovvo"]("i,a,b,j") += dp["aa_oo"]("k,i") * t2_2["aaaa"]("a,b,j,k");
        tmps_["82_aaaa_ovvo"]("i,a,b,j") += t2["aaaa"]("a,b,j,k") * tmps_["42_aa_oo"]("k,i");
        tmps_["82_aaaa_ovvo"]("i,a,b,j") += 0.50 * t2_1["aaaa"]("a,b,j,k") * tmps_["41_aa_oo"]("k,i");
        tmps_["82_aaaa_ovvo"]("i,a,b,j") -= 0.50 * t2["aaaa"]("a,b,i,k") * tmps_["81_aa_oo"]("k,j");
        tmps_["82_aaaa_ovvo"]("i,a,b,j") -= tmps_["38_aa_oo"]("k,j") * t2["aaaa"]("a,b,i,k");
        tmps_["82_aaaa_ovvo"]("i,a,b,j") += 0.25 * tmps_["6_aa_oo"]("l,i") * t2_1["aaaa"]("a,b,j,l");
        tmps_["82_aaaa_ovvo"]("i,a,b,j") -= 0.50 * tmps_["7_aa_oo"]("l,i") * t2_1["aaaa"]("a,b,j,l");
        tmps_["82_aaaa_ovvo"]("i,a,b,j") += tmps_["9_aa_oo"]("k,i") * t2_1["aaaa"]("a,b,j,k");
        tmps_["82_aaaa_ovvo"]("i,a,b,j") += 0.50 * tmps_["37_aa_oo"]("k,i") * t2["aaaa"]("a,b,j,k");
        tmps_["82_aaaa_ovvo"]("i,a,b,j") += 0.50 * tmps_["80_aa_oo"]("k,i") * t2["aaaa"]("a,b,j,k");
        tmps_["82_aaaa_ovvo"]("i,a,b,j") += 0.50 * tmps_["34_aa_oo"]("k,j") * t2["aaaa"]("a,b,i,k");
        tmps_["82_aaaa_ovvo"]("i,a,b,j") += 0.25 * tmps_["33_aa_oo"]("k,j") * t2["aaaa"]("a,b,i,k");
        tmps_["82_aaaa_ovvo"]("i,a,b,j") += 0.50 * tmps_["39_aa_oo"]("k,j") * t2["aaaa"]("a,b,i,k");

        // rt2_1_aaaa += +2.00 P(i,j) d-_aa(k,c) t2_1_aaaa(a,b,i,k) t1_1_aa(c,j)
        //            += -0.50 P(i,j) <l,k||c,d>_aaaa t2_aaaa(c,d,j,k) t2_1_aaaa(a,b,i,l)
        //            += -1.00 P(i,j) <k,l||j,c>_abab t2_aaaa(a,b,i,k) t1_1_bb(c,l)
        //            += -2.00 P(i,j) d-_aa(k,c) t2_aaaa(a,b,j,k) t1_2_aa(c,i)
        //            += -0.50 P(i,j) <l,k||c,d>_abab t2_abab(c,d,j,k) t2_1_aaaa(a,b,i,l)
        //            += -0.50 P(i,j) <l,k||d,c>_abab t2_abab(d,c,j,k) t2_1_aaaa(a,b,i,l)
        //            += -1.00 P(i,j) <l,k||c,j>_aaaa t2_aaaa(a,b,i,k) t1_1_aa(c,l)
        //            += +0.50 P(i,j) <k,l||c,d>_abab t2_aaaa(a,b,j,k) t2_1_abab(c,d,i,l)
        //            += +0.50 P(i,j) <k,l||d,c>_abab t2_aaaa(a,b,j,k) t2_1_abab(d,c,i,l)
        //            += -1.00 P(i,j) d-_aa(k,c) t2_aaaa(a,b,j,k) t0_1 t1_1_aa(c,i)
        //            += +1.00 P(i,j) d-_aa(k,j) t0_1 t2_1_aaaa(a,b,i,k)
        //            += +2.00 P(i,j) d-_aa(k,j) t2_aaaa(a,b,i,k) t0_2
        //            += +2.00 P(i,j) d-_aa(k,j) t2_2_aaaa(a,b,i,k)
        //            += +1.00 P(i,j) <a,b||c,j>_aaaa t1_1_aa(c,i)
        //            += +1.00 P(i,j) d+_aa(k,j) t2_aaaa(a,b,i,k)
        //            += -1.00 P(i,j) f_aa(k,j) t2_1_aaaa(a,b,i,k)
        //            += +0.50 P(i,j) <l,k||c,j>_aaaa t2_aaaa(a,b,l,k) t1_1_aa(c,i)
        //            += -0.50 P(i,j) <l,k||c,d>_aaaa t2_aaaa(a,b,j,k) t2_1_aaaa(c,d,i,l)
        //            += +1.00 P(i,j) f_aa(k,c) t2_aaaa(a,b,j,k) t1_1_aa(c,i)
        rt2_1_aaaa("a,b,i,j") += 2.00 * tmps_["82_aaaa_ovvo"]("j,a,b,i");
        rt2_1_aaaa("a,b,i,j") -= 2.00 * tmps_["82_aaaa_ovvo"]("i,a,b,j");
        tmps_["82_aaaa_ovvo"].~TArrayD();

        // flops: o4v0  = o4v2
        //  mems: o4v0  = o4v0
        tmps_["84_aaaa_oooo"]("i,j,k,l")  = eri["aaaa_oovv"]("i,j,a,b") * t2_1["aaaa"]("a,b,k,l");

        // flops: o4v0  = o4v2
        //  mems: o4v0  = o4v0
        tmps_["83_aaaa_oooo"]("i,j,k,l")  = eri["aaaa_oovv"]("i,j,a,b") * t2["aaaa"]("a,b,k,l");

        // flops: o2v2  = o2v4 o2v2 o2v2 o4v2 o2v2 o2v2 o2v2 o2v2 o2v2 o4v2 o2v2 o4v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["85_aaaa_vvoo"]("a,b,i,j")  = -0.50 * eri["aaaa_vvvv"]("a,b,c,d") * t2_1["aaaa"]("c,d,i,j");
        tmps_["85_aaaa_vvoo"]("a,b,i,j") += scalars_["1"] * t2_1["aaaa"]("a,b,i,j");
        tmps_["85_aaaa_vvoo"]("a,b,i,j") += 0.50 * eri["aaaa_oooo"]("l,k,i,j") * t2_1["aaaa"]("a,b,k,l");
        tmps_["85_aaaa_vvoo"]("a,b,i,j") += 2.00 * scalars_["5"] * t2_2["aaaa"]("a,b,i,j");
        tmps_["85_aaaa_vvoo"]("a,b,i,j") -= w0 * t2_1["aaaa"]("a,b,i,j");
        tmps_["85_aaaa_vvoo"]("a,b,i,j") += 0.25 * tmps_["84_aaaa_oooo"]("l,k,i,j") * t2["aaaa"]("a,b,k,l");
        tmps_["85_aaaa_vvoo"]("a,b,i,j") += 0.25 * tmps_["83_aaaa_oooo"]("l,k,i,j") * t2_1["aaaa"]("a,b,k,l");

        // rt2_1_aaaa += -1.00 d-_bb(k,c) t2_1_aaaa(a,b,i,j) t1_1_bb(c,k)
        //            += -1.00 d-_aa(k,c) t2_1_aaaa(a,b,i,j) t1_1_aa(c,k)
        //            += +0.50 <a,b||c,d>_aaaa t2_1_aaaa(c,d,i,j)
        //            += +0.50 <l,k||i,j>_aaaa t2_1_aaaa(a,b,l,k)
        //            += -2.00 d-_aa(k,k) t2_2_aaaa(a,b,i,j)
        //            += -2.00 d-_bb(k,k) t2_2_aaaa(a,b,i,j)
        //            += +1.00 t2_1_aaaa(a,b,i,j) w0
        //            += +0.25 <l,k||c,d>_aaaa t2_aaaa(a,b,l,k) t2_1_aaaa(c,d,i,j)
        //            += +0.25 <l,k||c,d>_aaaa t2_aaaa(c,d,i,j) t2_1_aaaa(a,b,l,k)
        rt2_1_aaaa("a,b,i,j") -= tmps_["85_aaaa_vvoo"]("a,b,i,j");
        tmps_["85_aaaa_vvoo"].~TArrayD();

        // flops: o2v0  = o2v0
        //  mems: o2v0  = o2v0
        tmps_["106_bb_oo"]("i,j")  = t0_1 * tmps_["24_bb_oo"]("i,j");

        // flops: o4v0  = o4v1 o4v1 o4v0 o4v2 o4v0
        //  mems: o4v0  = o4v0 o4v0 o4v0 o4v0 o4v0
        tmps_["105_abba_oooo"]("i,j,k,l")  = -1.00 * eri["abab_oovo"]("i,j,b,k") * t1_1["aa"]("b,l");
        tmps_["105_abba_oooo"]("i,j,k,l") += eri["abba_oovo"]("i,j,a,l") * t1_1["bb"]("a,k");
        tmps_["105_abba_oooo"]("i,j,k,l") -= eri["abab_oovv"]("i,j,b,c") * t2_1["abab"]("b,c,l,k");
}