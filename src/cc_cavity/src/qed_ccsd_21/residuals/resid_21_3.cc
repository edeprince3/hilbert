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

void QED_CCSD::resid_21_3() {
    
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


    // flops: o3v1  = o3v3 o3v2 o3v1
    //  mems: o3v1  = o3v1 o3v1 o3v1
    tmps_["63_bbbb_vooo"]("a,i,j,k")  = eri["bbbb_vovv"]("a,i,b,c") * t2["bbbb"]("b,c,j,k");
    tmps_["63_bbbb_vooo"]("a,i,j,k") -= 2.00 * t2["bbbb"]("b,a,j,k") * f["bb_ov"]("i,b");

    // flops: o3v1  = o3v2
    //  mems: o3v1  = o3v1
    tmps_["61_bbbb_vooo"]("a,i,j,k")  = t2_1["bbbb"]("b,a,i,j") * dp["bb_ov"]("k,b");

    // flops: o2v2  = o2v3 o2v2 o2v3 o2v2 o2v2 o3v2 o2v2 o2v3 o2v2 o3v2 o2v3 o2v2 o2v2 o3v2 o2v2 o2v3 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["66_bbbb_voov"]("a,i,j,b")  = -2.00 * tmps_["6_bb_vv"]("a,d") * t2_1["bbbb"]("d,b,i,j");
    tmps_["66_bbbb_voov"]("a,i,j,b") -= 2.00 * t0_1 * tmps_["62_bbbb_voov"]("a,i,j,b");
    tmps_["66_bbbb_voov"]("a,i,j,b") += t2["bbbb"]("c,a,i,j") * tmps_["64_bb_vv"]("b,c");
    tmps_["66_bbbb_voov"]("a,i,j,b") -= tmps_["63_bbbb_vooo"]("a,k,i,j") * t1_1["bb"]("b,k");
    tmps_["66_bbbb_voov"]("a,i,j,b") += tmps_["43_bb_vv"]("a,d") * t2_1["bbbb"]("d,b,i,j");
    tmps_["66_bbbb_voov"]("a,i,j,b") -= 2.00 * eri["bbbb_vooo"]("a,k,i,j") * t1_1["bb"]("b,k");
    tmps_["66_bbbb_voov"]("a,i,j,b") += 2.00 * f["bb_vv"]("a,c") * t2_1["bbbb"]("c,b,i,j");
    tmps_["66_bbbb_voov"]("a,i,j,b") += 4.00 * t1_1["bb"]("a,k") * tmps_["61_bbbb_vooo"]("b,i,j,k");
    tmps_["66_bbbb_voov"]("a,i,j,b") += 2.00 * tmps_["65_bb_vv"]("a,c") * t2["bbbb"]("c,b,i,j");
    tmps_["63_bbbb_vooo"].~TArrayD();
    tmps_["62_bbbb_voov"].~TArrayD();
    tmps_["61_bbbb_vooo"].~TArrayD();

    // rt2_1_bbbb += -0.50 P(a,b) <l,k||c,d>_bbbb t2_bbbb(c,a,i,j) t2_1_bbbb(d,b,l,k)
    //            += +0.50 P(a,b) <l,k||d,c>_abab t2_bbbb(c,a,i,j) t2_1_abab(d,b,l,k)
    //            += +0.50 P(a,b) <k,l||d,c>_abab t2_bbbb(c,a,i,j) t2_1_abab(d,b,k,l)
    //            += -1.00 P(a,b) d-_bb(k,c) t2_bbbb(c,a,i,j) t0_1 t1_1_bb(b,k)
    //            += -1.00 P(a,b) d-_bb(a,c) t0_1 t2_1_bbbb(c,b,i,j)
    //            += -0.50 P(a,b) <l,k||c,d>_abab t2_abab(c,a,l,k) t2_1_bbbb(d,b,i,j)
    //            += -0.50 P(a,b) <k,l||c,d>_abab t2_abab(c,a,k,l) t2_1_bbbb(d,b,i,j)
    //            += +0.50 P(a,b) <k,a||c,d>_bbbb t2_bbbb(c,d,i,j) t1_1_bb(b,k)
    //            += +1.00 P(a,b) f_bb(k,c) t2_bbbb(c,a,i,j) t1_1_bb(b,k)
    //            += -0.50 P(a,b) <l,k||c,d>_bbbb t2_bbbb(c,a,l,k) t2_1_bbbb(d,b,i,j)
    //            += +1.00 P(a,b) <k,a||i,j>_bbbb t1_1_bb(b,k)
    //            += +1.00 P(a,b) f_bb(a,c) t2_1_bbbb(c,b,i,j)
    //            += +2.00 P(a,b) d-_bb(k,c) t1_1_bb(a,k) t2_1_bbbb(c,b,i,j)
    //            += -1.00 P(a,b) <k,a||c,d>_bbbb t2_bbbb(c,b,i,j) t1_1_bb(d,k)
    //            += +1.00 P(a,b) <k,a||d,c>_abab t2_bbbb(c,b,i,j) t1_1_aa(d,k)
    rt2_1_bbbb("a,b,i,j") += 0.50 * tmps_["66_bbbb_voov"]("a,i,j,b");
    rt2_1_bbbb("a,b,i,j") -= 0.50 * tmps_["66_bbbb_voov"]("b,i,j,a");
    tmps_["66_bbbb_voov"].~TArrayD();

    // flops: o4v0  = o4v2
    //  mems: o4v0  = o4v0
    tmps_["69_bbbb_oooo"]("i,j,k,l")  = t2["bbbb"]("a,b,i,j") * eri["bbbb_oovv"]("k,l,a,b");

    // flops: o2v2  = o2v2
    //  mems: o2v2  = o2v2
    tmps_["68_bbbb_vvoo"]("a,b,i,j")  = -2.00 * scalars_["4"] * t2_1["bbbb"]("a,b,i,j");

    // flops: o4v0  = o4v2
    //  mems: o4v0  = o4v0
    tmps_["67_bbbb_oooo"]("i,j,k,l")  = t2_1["bbbb"]("a,b,i,j") * eri["bbbb_oovv"]("k,l,a,b");

    // flops: o2v2  = o2v4 o4v2 o2v2 o2v2 o2v2 o2v2 o4v2 o4v2 o2v2 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["70_bbbb_vvoo"]("a,b,i,j")  = -2.00 * eri["bbbb_vvvv"]("a,b,c,d") * t2_1["bbbb"]("c,d,i,j");
    tmps_["70_bbbb_vvoo"]("a,b,i,j") += 2.00 * eri["bbbb_oooo"]("l,k,i,j") * t2_1["bbbb"]("a,b,k,l");
    tmps_["70_bbbb_vvoo"]("a,b,i,j") -= 4.00 * w0 * t2_1["bbbb"]("a,b,i,j");
    tmps_["70_bbbb_vvoo"]("a,b,i,j") -= 2.00 * tmps_["68_bbbb_vvoo"]("a,b,i,j");
    tmps_["70_bbbb_vvoo"]("a,b,i,j") += t2["bbbb"]("a,b,k,l") * tmps_["67_bbbb_oooo"]("i,j,l,k");
    tmps_["70_bbbb_vvoo"]("a,b,i,j") += t2_1["bbbb"]("a,b,k,l") * tmps_["69_bbbb_oooo"]("i,j,l,k");
    tmps_["68_bbbb_vvoo"].~TArrayD();
    tmps_["67_bbbb_oooo"].~TArrayD();

    // rt2_1_bbbb += +0.25 <l,k||c,d>_bbbb t2_bbbb(c,d,i,j) t2_1_bbbb(a,b,l,k)
    //            += +0.25 <l,k||c,d>_bbbb t2_bbbb(a,b,l,k) t2_1_bbbb(c,d,i,j)
    //            += +0.50 <a,b||c,d>_bbbb t2_1_bbbb(c,d,i,j)
    //            += +0.50 <l,k||i,j>_bbbb t2_1_bbbb(a,b,l,k)
    //            += +1.00 t2_1_bbbb(a,b,i,j) w0
    //            += -1.00 d-_aa(k,c) t2_1_bbbb(a,b,i,j) t1_1_aa(c,k)
    //            += -1.00 d-_bb(k,c) t2_1_bbbb(a,b,i,j) t1_1_bb(c,k)
    rt2_1_bbbb("a,b,i,j") -= 0.25 * tmps_["70_bbbb_vvoo"]("a,b,i,j");
    tmps_["70_bbbb_vvoo"].~TArrayD();

    // flops: o0v2  = o2v3
    //  mems: o0v2  = o0v2
    tmps_["72_bb_vv"]("a,b")  = t2["bbbb"]("c,a,i,j") * eri["bbbb_oovv"]("j,i,b,c");
    tmps_["71_bbbb_vvoo"]("a,b,i,j")  = -2.00 * eri["bbbb_vvoo"]("a,b,i,j");

    // flops: o2v2  = o2v3 o4v2 o2v2 o4v2 o2v4 o2v2 o2v2 o2v2 o2v2 o2v2 o2v3 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["73_bbbb_vvoo"]("a,b,i,j")  = -2.00 * tmps_["43_bb_vv"]("a,d") * t2["bbbb"]("d,b,i,j");
    tmps_["73_bbbb_vvoo"]("a,b,i,j") += t2["bbbb"]("a,b,k,l") * tmps_["69_bbbb_oooo"]("i,j,l,k");
    tmps_["73_bbbb_vvoo"]("a,b,i,j") += 2.00 * eri["bbbb_oooo"]("l,k,i,j") * t2["bbbb"]("a,b,k,l");
    tmps_["73_bbbb_vvoo"]("a,b,i,j") -= 2.00 * eri["bbbb_vvvv"]("a,b,c,d") * t2["bbbb"]("c,d,i,j");
    tmps_["73_bbbb_vvoo"]("a,b,i,j") += 4.00 * scalars_["3"] * t2_1["bbbb"]("a,b,i,j");
    tmps_["73_bbbb_vvoo"]("a,b,i,j") += 2.00 * tmps_["71_bbbb_vvoo"]("a,b,i,j");
    tmps_["73_bbbb_vvoo"]("a,b,i,j") -= 2.00 * t2["bbbb"]("c,a,i,j") * tmps_["72_bb_vv"]("b,c");
    tmps_["71_bbbb_vvoo"].~TArrayD();
    tmps_["69_bbbb_oooo"].~TArrayD();

    // rt2_bbbb += +0.25 <l,k||c,d>_bbbb t2_bbbb(a,b,l,k) t2_bbbb(c,d,i,j)
    //          += -0.50 <l,k||c,d>_bbbb t2_bbbb(c,a,l,k) t2_bbbb(d,b,i,j)
    //          += +0.50 <l,k||i,j>_bbbb t2_bbbb(a,b,l,k)
    //          += +0.50 <a,b||c,d>_bbbb t2_bbbb(c,d,i,j)
    //          += -1.00 d-_bb(k,k) t2_1_bbbb(a,b,i,j)
    //          += -1.00 d-_aa(k,k) t2_1_bbbb(a,b,i,j)
    //          += +1.00 <a,b||i,j>_bbbb
    //          += -0.50 <l,k||c,d>_bbbb t2_bbbb(c,a,i,j) t2_bbbb(d,b,l,k)
    rt2_bbbb("a,b,i,j") -= 0.25 * tmps_["73_bbbb_vvoo"]("a,b,i,j");
    tmps_["73_bbbb_vvoo"].~TArrayD();

    // flops: o2v2  = o2v3 o2v3 o3v2 o3v2 o2v2 o2v2 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["74_aabb_ovvo"]("i,a,b,j")  = -1.00 * dp["aa_vv"]("a,c") * t2["abab"]("c,b,i,j");
    tmps_["74_aabb_ovvo"]("i,a,b,j") -= dp["bb_vv"]("b,d") * t2["abab"]("a,d,i,j");
    tmps_["74_aabb_ovvo"]("i,a,b,j") += dp["aa_oo"]("k,i") * t2["abab"]("a,b,k,j");
    tmps_["74_aabb_ovvo"]("i,a,b,j") += dp["bb_oo"]("l,j") * t2["abab"]("a,b,i,l");

    // rt2_1_abab  = +1.00 d+_aa(k,i) t2_abab(a,b,k,j)
    //            += +1.00 d+_bb(k,j) t2_abab(a,b,i,k)
    //            += -1.00 d+_bb(b,c) t2_abab(a,c,i,j)
    //            += -1.00 d+_aa(a,c) t2_abab(c,b,i,j)
    rt2_1_abab("a,b,i,j")  = tmps_["74_aabb_ovvo"]("i,a,b,j");

    // flops: o4v0  = o4v2
    //  mems: o4v0  = o4v0
    tmps_["76_abab_oooo"]("i,j,k,l")  = t2["abab"]("a,b,i,j") * eri["abab_oovv"]("k,l,a,b");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["75_baab_voov"]("a,i,j,b")  = t2["abab"]("c,a,i,k") * eri["abab_oovv"]("j,k,c,b");

    // flops: o2v2  = o2v3 o2v3 o2v3 o2v3 o2v2 o3v2 o2v2 o3v3 o2v2 o2v3 o2v2 o4v2 o2v2 o3v3 o2v2 o2v2 o2v2 o3v2 o2v2 o3v3 o2v2 o3v3 o2v2 o2v2 o2v2 o2v3 o2v2 o2v4 o2v2 o3v3 o2v2 o2v2 o2v2 o3v3 o2v2 o3v2 o2v2 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o2v2 o2v2 o2v2 o2v2 o3v2 o2v2 o3v3 o2v2 o3v3 o2v2 o4v2 o2v2 o2v2 o2v2 o2v2 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["77_abba_vvoo"]("a,b,i,j")  = -2.00 * tmps_["51_aa_vv"]("a,c") * t2["abab"]("c,b,j,i");
    tmps_["77_abba_vvoo"]("a,b,i,j") -= t2["abab"]("a,d,j,i") * tmps_["72_bb_vv"]("b,d");
    tmps_["77_abba_vvoo"]("a,b,i,j") -= 2.00 * t2["abab"]("a,d,j,i") * tmps_["6_bb_vv"]("b,d");
    tmps_["77_abba_vvoo"]("a,b,i,j") += tmps_["38_aa_vv"]("a,c") * t2["abab"]("c,b,j,i");
    tmps_["77_abba_vvoo"]("a,b,i,j") -= 2.00 * t1_1["aa"]("a,j") * tmps_["29_bb_vo"]("b,i");
    tmps_["77_abba_vvoo"]("a,b,i,j") += 2.00 * eri["abab_vvoo"]("a,b,j,i");
    tmps_["77_abba_vvoo"]("a,b,i,j") -= 2.00 * f["aa_oo"]("n,j") * t2["abab"]("a,b,n,i");
    tmps_["77_abba_vvoo"]("a,b,i,j") -= 2.00 * eri["abab_vovo"]("a,l,f,i") * t2["abab"]("f,b,j,l");
    tmps_["77_abba_vvoo"]("a,b,i,j") += 2.00 * f["bb_vv"]("b,d") * t2["abab"]("a,d,j,i");
    tmps_["77_abba_vvoo"]("a,b,i,j") += 2.00 * eri["abab_oooo"]("k,l,j,i") * t2["abab"]("a,b,k,l");
    tmps_["77_abba_vvoo"]("a,b,i,j") += 2.00 * eri["abba_vovo"]("a,l,d,j") * t2["bbbb"]("d,b,i,l");
    tmps_["77_abba_vvoo"]("a,b,i,j") -= 2.00 * dp["bb_vo"]("b,i") * t1_1["aa"]("a,j");
    tmps_["77_abba_vvoo"]("a,b,i,j") -= 2.00 * f["bb_oo"]("l,i") * t2["abab"]("a,b,j,l");
    tmps_["77_abba_vvoo"]("a,b,i,j") += 2.00 * eri["baab_vovo"]("b,n,f,i") * t2["aaaa"]("f,a,j,n");
    tmps_["77_abba_vvoo"]("a,b,i,j") -= 2.00 * eri["baba_vovo"]("b,n,d,j") * t2["abab"]("a,d,n,i");
    tmps_["77_abba_vvoo"]("a,b,i,j") -= 2.00 * dp["aa_vo"]("a,j") * t1_1["bb"]("b,i");
    tmps_["77_abba_vvoo"]("a,b,i,j") += 2.00 * f["aa_vv"]("a,f") * t2["abab"]("f,b,j,i");
    tmps_["77_abba_vvoo"]("a,b,i,j") += 2.00 * eri["abab_vvvv"]("a,b,f,e") * t2["abab"]("f,e,j,i");
    tmps_["77_abba_vvoo"]("a,b,i,j") -= 2.00 * eri["aaaa_vovo"]("a,n,f,j") * t2["abab"]("f,b,n,i");
    tmps_["77_abba_vvoo"]("a,b,i,j") -= 2.00 * scalars_["3"] * t2_1["abab"]("a,b,j,i");
    tmps_["77_abba_vvoo"]("a,b,i,j") -= 2.00 * eri["bbbb_vovo"]("b,l,d,i") * t2["abab"]("a,d,j,l");
    tmps_["77_abba_vvoo"]("a,b,i,j") += t2["abab"]("a,b,k,i") * tmps_["11_aa_oo"]("j,k");
    tmps_["77_abba_vvoo"]("a,b,i,j") += 2.00 * t2["bbbb"]("e,b,i,m") * tmps_["2_aabb_voov"]("a,j,m,e");
    tmps_["77_abba_vvoo"]("a,b,i,j") -= 2.00 * t2["bbbb"]("e,b,i,m") * tmps_["10_aabb_voov"]("a,j,m,e");
    tmps_["77_abba_vvoo"]("a,b,i,j") += 2.00 * t2["abab"]("a,d,j,l") * tmps_["4_bbbb_voov"]("b,i,l,d");
    tmps_["77_abba_vvoo"]("a,b,i,j") += 2.00 * t1_1["bb"]("b,i") * tmps_["25_aa_vo"]("a,j");
    tmps_["77_abba_vvoo"]("a,b,i,j") += t2["abab"]("a,b,j,m") * tmps_["16_bb_oo"]("i,m");
    tmps_["77_abba_vvoo"]("a,b,i,j") += 2.00 * t2["abab"]("a,d,n,i") * tmps_["75_baab_voov"]("b,j,n,d");
    tmps_["77_abba_vvoo"]("a,b,i,j") -= 2.00 * tmps_["9_aaaa_voov"]("a,j,k,c") * t2["abab"]("c,b,k,i");
    tmps_["77_abba_vvoo"]("a,b,i,j") += 2.00 * t2["abab"]("a,b,k,l") * tmps_["76_abab_oooo"]("j,i,k,l");
    tmps_["77_abba_vvoo"]("a,b,i,j") += 2.00 * t0_1 * tmps_["74_aabb_ovvo"]("j,a,b,i");
    tmps_["74_aabb_ovvo"].~TArrayD();
    tmps_["72_bb_vv"].~TArrayD();
    tmps_["29_bb_vo"].~TArrayD();
    tmps_["25_aa_vo"].~TArrayD();

    // rt2_abab  = -0.50 <l,k||c,d>_aaaa t2_abab(a,b,l,j) t2_aaaa(c,d,i,k)
    //          += -0.50 <l,k||c,d>_abab t2_abab(a,b,l,j) t2_abab(c,d,i,k)
    //          += -0.50 <l,k||d,c>_abab t2_abab(a,b,l,j) t2_abab(d,c,i,k)
    //          += +1.00 <a,b||i,j>_abab
    //          += -1.00 f_aa(k,i) t2_abab(a,b,k,j)
    //          += -1.00 <a,k||c,j>_abab t2_abab(c,b,i,k)
    //          += +1.00 f_bb(b,c) t2_abab(a,c,i,j)
    //          += +0.50 <l,k||i,j>_abab t2_abab(a,b,l,k)
    //          += +0.50 <k,l||i,j>_abab t2_abab(a,b,k,l)
    //          += -1.00 <a,k||i,c>_abab t2_bbbb(c,b,j,k)
    //          += -1.00 d-_bb(b,j) t1_1_aa(a,i)
    //          += -1.00 f_bb(k,j) t2_abab(a,b,i,k)
    //          += -1.00 <k,b||c,j>_abab t2_aaaa(c,a,i,k)
    //          += -1.00 <k,b||i,c>_abab t2_abab(a,c,k,j)
    //          += -1.00 d-_aa(a,i) t1_1_bb(b,j)
    //          += +1.00 f_aa(a,c) t2_abab(c,b,i,j)
    //          += +0.50 <a,b||c,d>_abab t2_abab(c,d,i,j)
    //          += +0.50 <a,b||d,c>_abab t2_abab(d,c,i,j)
    //          += +1.00 <k,a||c,i>_aaaa t2_abab(c,b,k,j)
    //          += -1.00 d-_bb(k,k) t2_1_abab(a,b,i,j)
    //          += -1.00 d-_aa(k,k) t2_1_abab(a,b,i,j)
    //          += +1.00 <k,b||c,j>_bbbb t2_abab(a,c,i,k)
    //          += -1.00 d-_aa(k,c) t2_abab(c,b,k,j) t1_1_aa(a,i)
    //          += +1.00 d-_bb(k,c) t2_bbbb(c,b,j,k) t1_1_aa(a,i)
    //          += +1.00 <k,l||c,d>_abab t2_aaaa(c,a,i,k) t2_bbbb(d,b,j,l)
    //          += +1.00 <l,k||c,d>_bbbb t2_abab(a,c,i,k) t2_bbbb(d,b,j,l)
    //          += +1.00 <l,k||d,c>_abab t2_abab(a,c,i,k) t2_abab(d,b,l,j)
    //          += +1.00 d-_aa(k,c) t2_aaaa(c,a,i,k) t1_1_bb(b,j)
    //          += -1.00 d-_bb(k,c) t2_abab(a,c,i,k) t1_1_bb(b,j)
    //          += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,a,l,k) t2_abab(d,b,i,j)
    //          += -0.50 <l,k||d,c>_abab t2_abab(a,c,i,j) t2_abab(d,b,l,k)
    //          += -0.50 <k,l||d,c>_abab t2_abab(a,c,i,j) t2_abab(d,b,k,l)
    //          += -0.50 <l,k||c,d>_bbbb t2_abab(a,b,i,l) t2_bbbb(c,d,j,k)
    //          += -0.50 <k,l||c,d>_abab t2_abab(a,b,i,l) t2_abab(c,d,k,j)
    //          += -0.50 <k,l||d,c>_abab t2_abab(a,b,i,l) t2_abab(d,c,k,j)
    //          += +1.00 <k,l||d,c>_abab t2_abab(a,c,k,j) t2_abab(d,b,i,l)
    //          += +1.00 <l,k||c,d>_aaaa t2_aaaa(c,a,i,k) t2_abab(d,b,l,j)
    //          += +0.25 <l,k||c,d>_abab t2_abab(a,b,l,k) t2_abab(c,d,i,j)
    //          += +0.25 <l,k||d,c>_abab t2_abab(a,b,l,k) t2_abab(d,c,i,j)
    //          += +0.25 <k,l||c,d>_abab t2_abab(a,b,k,l) t2_abab(c,d,i,j)
    //          += +0.25 <k,l||d,c>_abab t2_abab(a,b,k,l) t2_abab(d,c,i,j)
    //          += +1.00 d-_aa(k,i) t2_abab(a,b,k,j) t0_1
    //          += +1.00 d-_bb(k,j) t2_abab(a,b,i,k) t0_1
    //          += -1.00 d-_bb(b,c) t2_abab(a,c,i,j) t0_1
    //          += -1.00 d-_aa(a,c) t2_abab(c,b,i,j) t0_1
    //          += +0.50 <l,k||c,d>_bbbb t2_abab(a,c,i,j) t2_bbbb(d,b,l,k)
    //          += -0.50 <l,k||d,c>_abab t2_abab(a,c,l,k) t2_abab(d,b,i,j)
    //          += -0.50 <k,l||d,c>_abab t2_abab(a,c,k,l) t2_abab(d,b,i,j)
    rt2_abab("a,b,i,j")  = 0.50 * tmps_["77_abba_vvoo"]("a,b,j,i");
    tmps_["77_abba_vvoo"].~TArrayD();

    // flops: o2v2  = o2v3 o3v2 o3v2 o2v2
    //  mems: o2v2  = o2v2 o3v1 o2v2 o2v2
    tmps_["79_aaaa_vvoo"]("a,b,i,j")  = dp["aa_vv"]("a,c") * t2_1["aaaa"]("c,b,i,j");
    tmps_["79_aaaa_vvoo"]("a,b,i,j") += dp["aa_ov"]("k,c") * t2["aaaa"]("c,a,i,j") * t1_1["aa"]("b,k");

    // rt2_aaaa += -1.00 P(a,b) d-_aa(a,c) t2_1_aaaa(c,b,i,j)
    //          += -1.00 P(a,b) d-_aa(k,c) t2_aaaa(c,a,i,j) t1_1_aa(b,k)
    rt2_aaaa("a,b,i,j") -= tmps_["79_aaaa_vvoo"]("a,b,i,j");
    rt2_aaaa("a,b,i,j") += tmps_["79_aaaa_vvoo"]("b,a,i,j");

    // flops: o0v2  = o1v3 o1v3 o0v2
    //  mems: o0v2  = o0v2 o0v2 o0v2
    tmps_["82_aa_vv"]("a,b")  = eri["aaaa_vovv"]("a,i,b,c") * t1_1["aa"]("c,i");
    tmps_["82_aa_vv"]("a,b") += eri["abab_vovv"]("a,j,b,d") * t1_1["bb"]("d,j");

    // flops: o0v2  = o2v3 o2v3 o0v2
    //  mems: o0v2  = o0v2 o0v2 o0v2
    tmps_["81_aa_vv"]("a,b")  = 0.50 * t2_1["aaaa"]("d,a,i,k") * eri["aaaa_oovv"]("k,i,b,d");
    tmps_["81_aa_vv"]("a,b") += t2_1["abab"]("a,c,i,j") * eri["abab_oovv"]("i,j,b,c");

    // flops: o3v1  = o3v3 o3v2 o3v1
    //  mems: o3v1  = o3v1 o3v1 o3v1
    tmps_["80_aaaa_vooo"]("a,i,j,k")  = eri["aaaa_vovv"]("a,i,b,c") * t2["aaaa"]("b,c,j,k");
    tmps_["80_aaaa_vooo"]("a,i,j,k") -= 2.00 * t2["aaaa"]("b,a,j,k") * f["aa_ov"]("i,b");

    // flops: o3v1  = o3v2
    //  mems: o3v1  = o3v1
    tmps_["78_aaaa_vooo"]("a,i,j,k")  = t2_1["aaaa"]("b,a,i,j") * dp["aa_ov"]("k,b");

    // flops: o2v2  = o2v2 o2v3 o3v2 o2v3 o2v2 o2v2 o2v2 o2v3 o2v2 o3v2 o2v3 o2v2 o2v2 o2v3 o2v2 o3v2 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["83_aaaa_vvoo"]("a,b,i,j")  = -2.00 * t0_1 * tmps_["79_aaaa_vvoo"]("a,b,i,j");
    tmps_["83_aaaa_vvoo"]("a,b,i,j") += 2.00 * t2["aaaa"]("d,a,i,j") * tmps_["81_aa_vv"]("b,d");
    tmps_["83_aaaa_vvoo"]("a,b,i,j") -= t1_1["aa"]("b,k") * tmps_["80_aaaa_vooo"]("a,k,i,j");
    tmps_["83_aaaa_vvoo"]("a,b,i,j") += t2_1["aaaa"]("c,b,i,j") * tmps_["38_aa_vv"]("a,c");
    tmps_["83_aaaa_vvoo"]("a,b,i,j") += 2.00 * t2["aaaa"]("d,b,i,j") * tmps_["82_aa_vv"]("a,d");
    tmps_["83_aaaa_vvoo"]("a,b,i,j") -= 2.00 * eri["aaaa_vooo"]("a,k,i,j") * t1_1["aa"]("b,k");
    tmps_["83_aaaa_vvoo"]("a,b,i,j") += 2.00 * f["aa_vv"]("a,d") * t2_1["aaaa"]("d,b,i,j");
    tmps_["83_aaaa_vvoo"]("a,b,i,j") -= 2.00 * tmps_["51_aa_vv"]("a,c") * t2_1["aaaa"]("c,b,i,j");
    tmps_["83_aaaa_vvoo"]("a,b,i,j") += 4.00 * t1_1["aa"]("a,k") * tmps_["78_aaaa_vooo"]("b,i,j,k");
    tmps_["80_aaaa_vooo"].~TArrayD();
    tmps_["79_aaaa_vvoo"].~TArrayD();
    tmps_["78_aaaa_vooo"].~TArrayD();

    // rt2_1_aaaa += -0.50 P(a,b) <l,k||c,d>_aaaa t2_aaaa(c,a,l,k) t2_1_aaaa(d,b,i,j)
    //            += +0.50 P(a,b) <k,a||c,d>_aaaa t2_aaaa(c,d,i,j) t1_1_aa(b,k)
    //            += +1.00 P(a,b) f_aa(k,c) t2_aaaa(c,a,i,j) t1_1_aa(b,k)
    //            += +0.50 P(a,b) <l,k||c,d>_abab t2_aaaa(c,a,i,j) t2_1_abab(b,d,l,k)
    //            += +0.50 P(a,b) <k,l||c,d>_abab t2_aaaa(c,a,i,j) t2_1_abab(b,d,k,l)
    //            += -0.50 P(a,b) <l,k||c,d>_aaaa t2_aaaa(c,a,i,j) t2_1_aaaa(d,b,l,k)
    //            += -1.00 P(a,b) d-_aa(a,c) t0_1 t2_1_aaaa(c,b,i,j)
    //            += -1.00 P(a,b) d-_aa(k,c) t2_aaaa(c,a,i,j) t0_1 t1_1_aa(b,k)
    //            += -1.00 P(a,b) <k,a||c,d>_aaaa t2_aaaa(c,b,i,j) t1_1_aa(d,k)
    //            += +1.00 P(a,b) <a,k||c,d>_abab t2_aaaa(c,b,i,j) t1_1_bb(d,k)
    //            += +1.00 P(a,b) <k,a||i,j>_aaaa t1_1_aa(b,k)
    //            += +1.00 P(a,b) f_aa(a,c) t2_1_aaaa(c,b,i,j)
    //            += -0.50 P(a,b) <l,k||d,c>_abab t2_abab(a,c,l,k) t2_1_aaaa(d,b,i,j)
    //            += -0.50 P(a,b) <k,l||d,c>_abab t2_abab(a,c,k,l) t2_1_aaaa(d,b,i,j)
    //            += +2.00 P(a,b) d-_aa(k,c) t1_1_aa(a,k) t2_1_aaaa(c,b,i,j)
    rt2_1_aaaa("a,b,i,j") += 0.50 * tmps_["83_aaaa_vvoo"]("a,b,i,j");
    rt2_1_aaaa("a,b,i,j") -= 0.50 * tmps_["83_aaaa_vvoo"]("b,a,i,j");
    tmps_["83_aaaa_vvoo"].~TArrayD();

    // flops: o2v2  = o2v2
    //  mems: o2v2  = o2v2
    tmps_["85_aaaa_vvoo"]("a,b,i,j")  = -2.00 * w0 * t2_1["aaaa"]("a,b,i,j");

    // flops: o4v0  = o4v2
    //  mems: o4v0  = o4v0
    tmps_["84_aaaa_oooo"]("i,j,k,l")  = t2_1["aaaa"]("a,b,i,j") * eri["aaaa_oovv"]("k,l,a,b");

    // flops: o2v2  = o4v2 o4v2 o2v2 o2v2 o2v4 o2v2 o2v2 o2v2 o4v2 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["86_aaaa_vvoo"]("a,b,i,j")  = t2["aaaa"]("a,b,k,l") * tmps_["84_aaaa_oooo"]("i,j,l,k");
    tmps_["86_aaaa_vvoo"]("a,b,i,j") += 2.00 * eri["aaaa_oooo"]("l,k,i,j") * t2_1["aaaa"]("a,b,k,l");
    tmps_["86_aaaa_vvoo"]("a,b,i,j") += 4.00 * scalars_["4"] * t2_1["aaaa"]("a,b,i,j");
    tmps_["86_aaaa_vvoo"]("a,b,i,j") -= 2.00 * eri["aaaa_vvvv"]("a,b,c,d") * t2_1["aaaa"]("c,d,i,j");
    tmps_["86_aaaa_vvoo"]("a,b,i,j") += 2.00 * tmps_["85_aaaa_vvoo"]("a,b,i,j");
    tmps_["86_aaaa_vvoo"]("a,b,i,j") += t2_1["aaaa"]("a,b,k,l") * tmps_["37_aaaa_oooo"]("i,j,l,k");
    tmps_["85_aaaa_vvoo"].~TArrayD();
    tmps_["84_aaaa_oooo"].~TArrayD();
    tmps_["37_aaaa_oooo"].~TArrayD();

    // rt2_1_aaaa += +0.25 <l,k||c,d>_aaaa t2_aaaa(a,b,l,k) t2_1_aaaa(c,d,i,j)
    //            += +0.50 <l,k||i,j>_aaaa t2_1_aaaa(a,b,l,k)
    //            += -1.00 d-_aa(k,c) t2_1_aaaa(a,b,i,j) t1_1_aa(c,k)
    //            += -1.00 d-_bb(k,c) t2_1_aaaa(a,b,i,j) t1_1_bb(c,k)
    //            += +0.50 <a,b||c,d>_aaaa t2_1_aaaa(c,d,i,j)
    //            += +1.00 t2_1_aaaa(a,b,i,j) w0
    //            += +0.25 <l,k||c,d>_aaaa t2_aaaa(c,d,i,j) t2_1_aaaa(a,b,l,k)
    rt2_1_aaaa("a,b,i,j") -= 0.25 * tmps_["86_aaaa_vvoo"]("a,b,i,j");
    tmps_["86_aaaa_vvoo"].~TArrayD();

    // flops: o2v2  = o2v1 o3v2 o3v2 o2v2
    //  mems: o2v2  = o2v0 o2v2 o2v2 o2v2
    tmps_["88_aaaa_ovvo"]("i,a,b,j")  = dp["aa_ov"]("k,c") * t1_1["aa"]("c,i") * t2["aaaa"]("a,b,j,k");
    tmps_["88_aaaa_ovvo"]("i,a,b,j") -= dp["aa_oo"]("k,j") * t2_1["aaaa"]("a,b,i,k");

    // rt2_aaaa += -1.00 P(i,j) d-_aa(k,c) t2_aaaa(a,b,j,k) t1_1_aa(c,i)
    //          += +1.00 P(i,j) d-_aa(k,j) t2_1_aaaa(a,b,i,k)
    rt2_aaaa("a,b,i,j") -= tmps_["88_aaaa_ovvo"]("i,a,b,j");
    rt2_aaaa("a,b,i,j") += tmps_["88_aaaa_ovvo"]("j,a,b,i");

    // flops: o2v0  = o3v1 o3v1 o2v0
    //  mems: o2v0  = o2v0 o2v0 o2v0
    tmps_["90_aa_oo"]("i,j")  = eri["aaaa_oovo"]("i,l,b,j") * t1_1["aa"]("b,l");
    tmps_["90_aa_oo"]("i,j") += eri["abba_oovo"]("i,k,a,j") * t1_1["bb"]("a,k");

    // flops: o2v0  = o3v2 o2v1 o3v2 o2v0 o2v0
    //  mems: o2v0  = o2v0 o2v0 o2v0 o2v0 o2v0
    tmps_["89_aa_oo"]("i,j")  = 0.50 * t2_1["aaaa"]("a,c,i,l") * eri["aaaa_oovv"]("j,l,a,c");
    tmps_["89_aa_oo"]("i,j") += t1_1["aa"]("a,i") * f["aa_ov"]("j,a");
    tmps_["89_aa_oo"]("i,j") += t2_1["abab"]("a,b,i,k") * eri["abab_oovv"]("j,k,a,b");

    // flops: o4v0  = o4v1
    //  mems: o4v0  = o4v0
    tmps_["87_aaaa_oooo"]("i,j,k,l")  = eri["aaaa_oovo"]("i,j,a,k") * t1_1["aa"]("a,l");

    // flops: o2v2  = o3v2 o2v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o4v2 o2v2 o2v3 o3v2 o2v2 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["91_aaaa_ovvo"]("i,a,b,j")  = -0.50 * t2_1["aaaa"]("a,b,i,k") * tmps_["11_aa_oo"]("j,k");
    tmps_["91_aaaa_ovvo"]("i,a,b,j") += t0_1 * tmps_["88_aaaa_ovvo"]("i,a,b,j");
    tmps_["91_aaaa_ovvo"]("i,a,b,j") -= 2.00 * t2_1["aaaa"]("a,b,i,l") * tmps_["52_aa_oo"]("j,l");
    tmps_["91_aaaa_ovvo"]("i,a,b,j") -= t2["aaaa"]("a,b,i,l") * tmps_["90_aa_oo"]("l,j");
    tmps_["91_aaaa_ovvo"]("i,a,b,j") -= t2["aaaa"]("a,b,j,l") * tmps_["89_aa_oo"]("i,l");
    tmps_["91_aaaa_ovvo"]("i,a,b,j") += 0.50 * t2["aaaa"]("a,b,k,l") * tmps_["87_aaaa_oooo"]("l,k,j,i");
    tmps_["91_aaaa_ovvo"]("i,a,b,j") -= eri["aaaa_vvvo"]("a,b,c,j") * t1_1["aa"]("c,i");
    tmps_["91_aaaa_ovvo"]("i,a,b,j") += f["aa_oo"]("l,j") * t2_1["aaaa"]("a,b,i,l");
    tmps_["88_aaaa_ovvo"].~TArrayD();
    tmps_["87_aaaa_oooo"].~TArrayD();

    // rt2_1_aaaa += -1.00 P(i,j) d-_aa(k,c) t2_aaaa(a,b,j,k) t0_1 t1_1_aa(c,i)
    //            += +1.00 P(i,j) d-_aa(k,j) t0_1 t2_1_aaaa(a,b,i,k)
    //            += -0.50 P(i,j) <l,k||c,d>_aaaa t2_aaaa(c,d,j,k) t2_1_aaaa(a,b,i,l)
    //            += -0.50 P(i,j) <l,k||c,d>_abab t2_abab(c,d,j,k) t2_1_aaaa(a,b,i,l)
    //            += -0.50 P(i,j) <l,k||d,c>_abab t2_abab(d,c,j,k) t2_1_aaaa(a,b,i,l)
    //            += +2.00 P(i,j) d-_aa(k,c) t2_1_aaaa(a,b,i,k) t1_1_aa(c,j)
    //            += -1.00 P(i,j) <k,l||j,c>_abab t2_aaaa(a,b,i,k) t1_1_bb(c,l)
    //            += -1.00 P(i,j) <l,k||c,j>_aaaa t2_aaaa(a,b,i,k) t1_1_aa(c,l)
    //            += +0.50 P(i,j) <k,l||c,d>_abab t2_aaaa(a,b,j,k) t2_1_abab(c,d,i,l)
    //            += +0.50 P(i,j) <k,l||d,c>_abab t2_aaaa(a,b,j,k) t2_1_abab(d,c,i,l)
    //            += +1.00 P(i,j) f_aa(k,c) t2_aaaa(a,b,j,k) t1_1_aa(c,i)
    //            += -0.50 P(i,j) <l,k||c,d>_aaaa t2_aaaa(a,b,j,k) t2_1_aaaa(c,d,i,l)
    //            += +0.50 P(i,j) <l,k||c,j>_aaaa t2_aaaa(a,b,l,k) t1_1_aa(c,i)
    //            += +1.00 P(i,j) <a,b||c,j>_aaaa t1_1_aa(c,i)
    //            += -1.00 P(i,j) f_aa(k,j) t2_1_aaaa(a,b,i,k)
    rt2_1_aaaa("a,b,i,j") -= tmps_["91_aaaa_ovvo"]("i,a,b,j");
    rt2_1_aaaa("a,b,i,j") += tmps_["91_aaaa_ovvo"]("j,a,b,i");
    tmps_["91_aaaa_ovvo"].~TArrayD();

    // flops: o2v2  = o2v3
    //  mems: o2v2  = o2v2
    tmps_["97_bbbb_vovo"]("a,i,b,j")  = eri["bbbb_vovv"]("a,i,b,c") * t1_1["bb"]("c,j");

    // flops: o2v2  = o2v3
    //  mems: o2v2  = o2v2
    tmps_["96_baab_vovo"]("a,i,b,j")  = eri["baab_vovv"]("a,i,b,c") * t1_1["bb"]("c,j");
}
