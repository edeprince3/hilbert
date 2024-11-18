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

void QED_CCSD_22::resid_22_5() {

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


        // flops: o2v2  = o3v2 o3v3 o2v2 o3v2 o2v2 o2v3 o2v2 o3v3 o2v2 o2v3 o2v2 o3v3 o2v2 o3v2 o2v2 o3v2 o2v2 o3v3 o2v2 o2v3 o2v2 o3v2 o2v2 o3v2 o2v2 o3v3 o2v2 o2v2 o2v2 o3v3 o2v2 o2v2 o2v2 o2v3 o2v2 o3v2 o2v2 o2v3 o2v2 o3v2 o2v2 o3v3 o2v2 o2v3 o2v2 o3v3 o2v2 o2v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v2 o2v2 o4v2 o2v2 o2v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v2 o2v2 o3v2 o2v2 o2v3 o2v2 o3v2 o2v2 o3v2 o2v2 o3v3 o2v2 o2v2 o2v2 o3v2 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v2 o2v2 o3v3 o2v2 o3v2 o2v2 o3v2 o2v2 o3v3 o2v2 o3v2 o2v2 o3v3 o2v2 o2v3 o2v2 o2v3 o3v3 o2v2 o2v3 o3v3 o2v2 o2v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o4v2 o2v2 o3v3 o2v2 o1v4 o2v3 o2v2 o2v4 o2v2 o2v3 o2v2 o3v3 o2v2 o2v3 o2v2 o2v3 o2v2 o3v3 o2v2 o2v2 o2v2 o3v3 o2v2 o3v2 o2v2 o2v2 o2v2 o3v3 o2v2 o2v2 o2v2 o2v2 o3v3 o2v2 o2v2 o2v2 o3v2 o2v2 o3v2 o2v2 o2v3 o2v2 o3v2 o2v2 o2v3 o2v2 o2v2 o2v2 o3v2 o2v2 o3v2 o2v2 o2v3 o2v2 o3v2 o2v2 o2v3 o2v2 o3v3 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o2v2 o3v2 o2v2 o2v3 o2v2 o3v2 o2v2 o4v2 o2v2 o4v2 o2v2 o3v3 o2v2 o2v2 o2v2 o2v2 o2v2 o2v3 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o2v3 o2v2 o2v3 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o3v2 o2v2 o3v2 o2v2 o3v2 o3v2 o2v2 o2v2 o2v2 o3v2 o2v2 o2v2 o3v2 o2v2 o2v3 o2v2 o2v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v2 o2v2 o2v2 o2v2 o3v2 o2v2 o3v2 o3v2 o2v2 o4v2 o2v2 o2v2 o2v2 o3v2 o2v2 o2v2 o2v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o1v3 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["191_abab_vvoo"]("a,b,i,j")  = 3.00 * tmps_["94_abab_ovoo"]("l,b,i,j") * t1_2["aa"]("a,l");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["163_baab_vovo"]("b,l,c,j") * t2["aaaa"]("c,a,i,l");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 3.00 * tmps_["9_aa_oo"]("l,i") * t2_2["abab"]("a,b,l,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["72_aa_vv"]("a,f") * t2_1["abab"]("f,b,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["17_aabb_ovvo"]("n,f,b,j") * t2_2["aaaa"]("f,a,i,n");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["32_aa_vv"]("f,a") * t2_1["abab"]("f,b,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["47_aabb_ovvo"]("l,c,b,j") * t2_1["aaaa"]("c,a,i,l");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["59_bb_oo"]("m,j") * t2_1["abab"]("a,b,i,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 3.00 * tmps_["38_aa_oo"]("l,i") * t2_1["abab"]("a,b,l,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["121_aaaa_vovo"]("a,l,c,i") * t2["abab"]("c,b,l,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 0.50 * tmps_["4_aa_vv"]("f,a") * t2_2["abab"]("f,b,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["36_aa_oo"]("n,i") * t2_1["abab"]("a,b,n,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["34_aa_oo"]("n,i") * t2_1["abab"]("a,b,n,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["89_abab_vovo"]("a,m,f,j") * t2_1["abab"]("f,b,i,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["170_bb_vo"]("b,j") * t1_1["aa"]("a,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["2_bbaa_ovvo"]("k,e,a,i") * t2_2["bbbb"]("e,b,j,k");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["171_bb_ov"]("j,b") * t1_1["aa"]("a,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["95_bb_vv"]("b,e") * t2_1["abab"]("a,e,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 0.50 * tmps_["6_aa_oo"]("n,i") * t2_2["abab"]("a,b,n,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["50_bb_vv"]("d,b") * t2_1["abab"]("a,d,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 0.50 * tmps_["23_bb_oo"]("k,j") * t2_2["abab"]("a,b,i,k");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["18_bbbb_ovvo"]("k,e,b,j") * t2_2["abab"]("a,e,i,k");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["167_bb_vv"]("b,e") * t2_1["abab"]("a,e,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["48_bbbb_ovvo"]("m,d,b,j") * t2_1["abab"]("a,d,i,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= 0.50 * tmps_["52_bb_vv"]("d,b") * t2_1["abab"]("a,d,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["87_abba_ovvo"]("l,d,b,i") * t2_1["abab"]("a,d,l,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["31_abba_vovo"]("a,m,e,i") * t2_1["bbbb"]("e,b,j,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["54_bb_oo"]("k,j") * t2_1["abab"]("a,b,i,k");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["88_abab_oooo"]("n,m,i,j") * t2_2["abab"]("a,b,n,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 0.50 * tmps_["144_aa_vv"]("f,a") * t2_1["abab"]("f,b,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["28_bbaa_ovvo"]("k,e,a,i") * t2_1["bbbb"]("e,b,j,k");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["49_bbbb_ovvo"]("m,d,b,j") * t2_1["abab"]("a,d,i,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 3.00 * tmps_["24_bb_oo"]("m,j") * t2_2["abab"]("a,b,i,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["57_bb_oo"]("k,j") * t2_1["abab"]("a,b,i,k");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["153_aa_vv"]("a,f") * t2_1["abab"]("f,b,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 0.50 * tmps_["134_aa_oo"]("n,i") * t2_1["abab"]("a,b,n,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["37_aa_oo"]("n,i") * t2_1["abab"]("a,b,n,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["120_bbaa_ovvo"]("m,d,a,i") * t2["bbbb"]("d,b,j,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["125_aa_vo"]("a,i") * t1_1["bb"]("b,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["39_aa_oo"]("l,i") * t2_1["abab"]("a,b,l,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["122_aaaa_vovo"]("a,l,f,i") * t2_1["abab"]("f,b,l,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["123_abba_vovo"]("a,m,d,i") * t2["bbbb"]("d,b,j,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["165_bbbb_vovo"]("b,m,d,j") * t2["abab"]("a,d,i,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 3.00 * tmps_["58_bb_oo"]("m,j") * t2_1["abab"]("a,b,i,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["162_bbbb_ovvo"]("m,d,b,j") * t2["abab"]("a,d,i,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 0.50 * tmps_["166_bb_oo"]("k,j") * t2_1["abab"]("a,b,i,k");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["7_aa_oo"]("n,i") * t2_2["abab"]("a,b,n,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["90_baba_vovo"]("b,l,e,i") * t2_1["abab"]("a,e,l,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["56_bb_oo"]("k,j") * t2_1["abab"]("a,b,i,k");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["19_bbbb_ovvo"]("k,e,b,j") * t2_2["abab"]("a,e,i,k");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 0.50 * tmps_["21_bb_vv"]("e,b") * t2_2["abab"]("a,e,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += eri["abab_vvvo"]("a,b,c,j") * t1_2["aa"]("c,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= eri["aaaa_vovo"]("a,l,c,i") * t2_2["abab"]("c,b,l,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= eri["abab_vovv"]("a,m,c,e") * t1_2["bb"]("e,j") * t2["abab"]("c,b,i,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 2.00 * w0 * t2_2["abab"]("a,b,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= f["aa_oo"]("l,i") * t2_2["abab"]("a,b,l,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += eri["baab_vooo"]("b,l,i,j") * t1_2["aa"]("a,l");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= f["bb_oo"]("m,j") * t2_2["abab"]("a,b,i,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") += eri["abab_oooo"]("n,m,i,j") * t2_2["abab"]("a,b,n,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= eri["abab_vovo"]("a,m,c,j") * t2_2["abab"]("c,b,i,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") += eri["abab_vvvv"]("a,b,f,d") * t1_1["aa"]("f,i") * t1_1["bb"]("d,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += eri["abab_vvvv"]("a,b,c,e") * t2_2["abab"]("c,e,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += f["bb_vv"]("b,d") * t2_2["abab"]("a,d,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += eri["baab_vovo"]("b,l,c,j") * t2_2["aaaa"]("c,a,i,l");
        tmps_["191_abab_vvoo"]("a,b,i,j") += f["aa_vv"]("a,c") * t2_2["abab"]("c,b,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= eri["abba_vvvo"]("a,b,d,i") * t1_2["bb"]("d,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += eri["abba_vovo"]("a,m,d,i") * t2_2["bbbb"]("d,b,j,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= 2.00 * scalars_["1"] * t2_2["abab"]("a,b,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= eri["baba_vovo"]("b,l,d,i") * t2_2["abab"]("a,d,l,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= eri["abab_vooo"]("a,m,i,j") * t1_2["bb"]("b,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= scalars_["6"] * t2_1["abab"]("a,b,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= eri["bbbb_vovo"]("b,m,d,j") * t2_2["abab"]("a,d,i,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") += t1_1["bb"]("b,j") * tmps_["128_aa_vo"]("a,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") += t1_1["aa"]("a,i") * tmps_["174_bb_vo"]("b,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += t2["abab"]("a,d,l,j") * tmps_["181_abba_ovvo"]("l,d,b,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 2.00 * t1_2["aa"]("a,i") * tmps_["103_bb_vo"]("b,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 2.00 * tmps_["91_baab_ovoo"]("m,a,i,j") * tmps_["177_bb_vo"]("b,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 2.00 * t2_1["abab"]("a,b,i,m") * tmps_["62_bb_oo"]("m,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= t2["abab"]("c,b,i,j") * tmps_["155_aa_vv"]("c,a");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= t2["abab"]("a,b,i,m") * tmps_["176_bb_oo"]("m,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= t2_2["abab"]("c,b,i,j") * tmps_["10_aa_vv"]("a,c");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= 2.00 * t1_2["bb"]("b,j") * tmps_["76_aa_ov"]("i,a");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= t2["abab"]("a,b,i,m") * tmps_["175_bb_oo"]("m,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= t1_2["bb"]("b,m") * tmps_["101_abab_vooo"]("a,m,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= 2.00 * t2_1["abab"]("c,b,i,j") * tmps_["40_aa_vv"]("a,c");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= t1_1["bb"]("b,m") * tmps_["180_abab_vooo"]("a,m,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= 2.00 * t2_1["abab"]("a,d,i,j") * tmps_["60_bb_vv"]("b,d");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= t2_2["abab"]("f,b,n,j") * tmps_["1_aaaa_ovvo"]("n,f,a,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") += t2_2["abab"]("a,b,l,j") * tmps_["41_aa_oo"]("l,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= t1_2["aa"]("a,n") * tmps_["100_abba_oovo"]("n,j,b,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= t1_1["aa"]("a,l") * tmps_["179_aabb_ooov"]("l,i,j,b");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 2.00 * t2_1["abab"]("a,b,l,j") * tmps_["42_aa_oo"]("l,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") += t2["abab"]("a,b,l,j") * tmps_["140_aa_oo"]("l,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") += t2["abab"]("c,b,i,j") * tmps_["156_aa_vv"]("a,c");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["12_aa_vo"]("a,l") * tmps_["94_abab_ovoo"]("l,b,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= t2_1["abab"]("a,b,n,m") * tmps_["105_abba_oooo"]("n,m,j,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") += t2["abab"]("a,b,n,m") * tmps_["178_abab_oooo"]("n,m,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= t2_2["bbbb"]("e,b,j,k") * tmps_["3_bbaa_ovvo"]("k,e,a,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= 2.00 * t1_2["bb"]("b,j") * tmps_["74_aa_vo"]("a,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 2.00 * t1_2["aa"]("a,i") * tmps_["104_bb_vo"]("b,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= t2["abab"]("a,d,i,j") * tmps_["173_bb_vv"]("b,d");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["93_abab_ovoo"]("l,b,i,j") * tmps_["157_aa_vo"]("a,l");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= t2["abab"]("a,b,l,j") * tmps_["139_aa_oo"]("l,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= t1_2["bb"]("b,k") * tmps_["98_bbaa_oovo"]("k,j,a,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["26_bb_vo"]("b,m") * tmps_["92_baab_ovoo"]("m,a,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= t2["abab"]("a,d,i,j") * tmps_["172_bb_vv"]("d,b");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= t2_2["abab"]("a,d,i,j") * tmps_["25_bb_vv"]("b,d");
        tmps_["191_abab_vvoo"]("a,b,i,j") += t2_2["abab"]("a,b,i,m") * tmps_["61_bb_oo"]("m,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += t1_2["aa"]("a,l") * tmps_["99_abab_ovoo"]("l,b,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += t1_1["aa"]("a,n") * tmps_["183_baba_vooo"]("b,i,j,n");
        tmps_["191_abab_vvoo"]("a,b,i,j") += t2_1["abab"]("a,b,l,j") * tmps_["81_aa_oo"]("l,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") += t2_1["abab"]("a,b,i,m") * tmps_["106_bb_oo"]("m,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += t2["abab"]("a,b,i,m") * tmps_["182_bb_oo"]("m,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += t1_1["bb"]("b,k") * tmps_["184_bbaa_oovo"]("j,k,a,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 2.00 * t2["abab"]("a,b,l,j") * tmps_["141_aa_oo"]("l,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["22_bb_oo"]("k,j") * t2_2["abab"]("a,b,i,k");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["20_bb_vv"]("e,b") * t2_2["abab"]("a,e,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["5_aa_vv"]("f,a") * t2_2["abab"]("f,b,i,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["164_bbbb_vovo"]("b,m,e,j") * t2_1["abab"]("a,e,i,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["86_abba_ovvo"]("n,e,b,i") * t2_2["abab"]("a,e,n,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["51_baab_vovo"]("b,l,f,j") * t2_1["aaaa"]("f,a,i,l");
        tmps_["191_abab_vvoo"]("a,b,i,j") += 3.00 * tmps_["92_baab_ovoo"]("m,a,i,j") * t1_2["bb"]("b,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["126_aa_ov"]("i,a") * t1_1["bb"]("b,j");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= t1_1["bb"]("b,k") * tmps_["190_abba_vooo"]("a,j,k,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") -= tmps_["186_bb_oo"]("j,m") * t2["abab"]("a,b,i,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["189_bbaa_vooo"]("b,j,n,i") * t1_1["aa"]("a,n");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["187_baba_oooo"]("j,n,m,i") * t2["abab"]("a,b,n,m");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["185_bb_vo"]("b,j") * t1_1["aa"]("a,i");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["188_baab_vooo"]("b,l,i,j") * t1_1["aa"]("a,l");
        tmps_["191_abab_vvoo"]("a,b,i,j") += tmps_["131_aa_vo"]("a,i") * t1_1["bb"]("b,j");
        tmps_["190_abba_vooo"].~TArrayD();
        tmps_["189_bbaa_vooo"].~TArrayD();
        tmps_["188_baab_vooo"].~TArrayD();
        tmps_["187_baba_oooo"].~TArrayD();
        tmps_["184_bbaa_oovo"].~TArrayD();
        tmps_["183_baba_vooo"].~TArrayD();
        tmps_["181_abba_ovvo"].~TArrayD();
        tmps_["180_abab_vooo"].~TArrayD();
        tmps_["179_aabb_ooov"].~TArrayD();
        tmps_["178_abab_oooo"].~TArrayD();
        tmps_["157_aa_vo"].~TArrayD();
        tmps_["156_aa_vv"].~TArrayD();
        tmps_["155_aa_vv"].~TArrayD();
        tmps_["153_aa_vv"].~TArrayD();
        tmps_["144_aa_vv"].~TArrayD();
        tmps_["141_aa_oo"].~TArrayD();
        tmps_["140_aa_oo"].~TArrayD();
        tmps_["139_aa_oo"].~TArrayD();
        tmps_["134_aa_oo"].~TArrayD();
        tmps_["131_aa_vo"].~TArrayD();
        tmps_["128_aa_vo"].~TArrayD();
        tmps_["126_aa_ov"].~TArrayD();
        tmps_["125_aa_vo"].~TArrayD();
        tmps_["123_abba_vovo"].~TArrayD();
        tmps_["122_aaaa_vovo"].~TArrayD();
        tmps_["121_aaaa_vovo"].~TArrayD();
        tmps_["120_bbaa_ovvo"].~TArrayD();
        tmps_["105_abba_oooo"].~TArrayD();
        tmps_["101_abab_vooo"].~TArrayD();
        tmps_["100_abba_oovo"].~TArrayD();
        tmps_["99_abab_ovoo"].~TArrayD();
        tmps_["98_bbaa_oovo"].~TArrayD();
        tmps_["94_abab_ovoo"].~TArrayD();
        tmps_["92_baab_ovoo"].~TArrayD();
        tmps_["90_baba_vovo"].~TArrayD();
        tmps_["89_abab_vovo"].~TArrayD();
        tmps_["87_abba_ovvo"].~TArrayD();
        tmps_["81_aa_oo"].~TArrayD();
        tmps_["76_aa_ov"].~TArrayD();
        tmps_["74_aa_vo"].~TArrayD();
        tmps_["72_aa_vv"].~TArrayD();
        tmps_["42_aa_oo"].~TArrayD();
        tmps_["40_aa_vv"].~TArrayD();
        tmps_["39_aa_oo"].~TArrayD();
        tmps_["38_aa_oo"].~TArrayD();
        tmps_["37_aa_oo"].~TArrayD();
        tmps_["36_aa_oo"].~TArrayD();
        tmps_["34_aa_oo"].~TArrayD();
        tmps_["32_aa_vv"].~TArrayD();
        tmps_["31_abba_vovo"].~TArrayD();
        tmps_["28_bbaa_ovvo"].~TArrayD();
        tmps_["12_aa_vo"].~TArrayD();

        // rt2_2_abab  = +1.00 <l,k||d,c>_abab t2_abab(a,b,l,k) t1_1_bb(c,j) t1_1_aa(d,i)
        //            += +1.00 <k,l||d,c>_abab t2_abab(a,b,k,l) t1_1_bb(c,j) t1_1_aa(d,i)
        //            += +2.00 <l,k||d,c>_abab t2_bbbb(c,b,j,k) t1_1_aa(a,l) t1_1_aa(d,i)
        //            += +2.00 d-_aa(k,c) t1_1_aa(a,k) t1_1_bb(b,j) t1_1_aa(c,i)
        //            += +2.00 d-_bb(k,c) t1_1_aa(a,i) t1_1_bb(b,k) t1_1_bb(c,j)
        //            += +1.00 <k,l||c,d>_abab t2_abab(c,d,i,j) t1_1_aa(a,k) t1_1_bb(b,l)
        //            += +1.00 <k,l||d,c>_abab t2_abab(d,c,i,j) t1_1_aa(a,k) t1_1_bb(b,l)
        //            += -2.00 <l,k||c,d>_abab t2_abab(a,b,i,k) t1_1_aa(c,l) t1_1_bb(d,j)
        //            += -2.00 <k,l||d,c>_abab t2_abab(a,c,i,j) t1_1_bb(b,l) t1_1_aa(d,k)
        //            += -2.00 <l,k||c,d>_bbbb t2_abab(a,c,i,j) t1_1_bb(b,l) t1_1_bb(d,k)
        //            += +2.00 <k,l||d,c>_abab t2_abab(a,c,k,j) t1_1_bb(b,l) t1_1_aa(d,i)
        //            += -2.00 <k,b||c,d>_abab t2_aaaa(c,a,i,k) t1_2_bb(d,j)
        //            += +6.00 d-_aa(k,c) t2_1_abab(c,b,i,j) t1_2_aa(a,k)
        //            += +6.00 d-_aa(k,c) t1_1_aa(c,i) t2_2_abab(a,b,k,j)
        //            += +2.00 <a,k||d,c>_abab t2_1_abab(d,b,i,j) t1_1_bb(c,k)
        //            += +6.00 d-_aa(k,c) t2_1_abab(a,b,k,j) t1_2_aa(c,i)
        //            += +2.00 <l,k||c,d>_aaaa t2_abab(c,b,k,j) t2_2_aaaa(d,a,i,l)
        //            += -1.00 <l,k||d,c>_abab t2_1_abab(a,c,l,k) t2_1_abab(d,b,i,j)
        //            += -1.00 <k,l||d,c>_abab t2_1_abab(a,c,k,l) t2_1_abab(d,b,i,j)
        //            += +2.00 <l,k||c,d>_aaaa t2_1_aaaa(c,a,i,k) t2_1_abab(d,b,l,j)
        //            += -2.00 f_bb(k,c) t2_1_abab(a,b,i,k) t1_1_bb(c,j)
        //            += +2.00 <k,l||c,d>_abab t2_aaaa(c,a,i,k) t2_2_bbbb(d,b,j,l)
        //            += +2.00 <k,a||c,d>_aaaa t2_abab(c,b,k,j) t1_2_aa(d,i)
        //            += -1.00 <l,k||c,d>_aaaa t2_aaaa(c,a,l,k) t2_2_abab(d,b,i,j)
        //            += +2.00 <l,k||c,i>_aaaa t2_1_abab(a,b,l,j) t1_1_aa(c,k)
        //            += -1.00 <l,k||c,d>_abab t2_1_abab(a,b,l,j) t2_1_abab(c,d,i,k)
        //            += -1.00 <l,k||d,c>_abab t2_1_abab(a,b,l,j) t2_1_abab(d,c,i,k)
        //            += -2.00 <a,k||d,c>_abab t2_1_abab(d,b,i,k) t1_1_bb(c,j)
        //            += -2.00 d-_bb(b,c) t1_1_aa(a,i) t1_2_bb(c,j)
        //            += +2.00 d-_bb(k,j) t1_1_aa(a,i) t1_2_bb(b,k)
        //            += +2.00 <k,l||c,d>_abab t2_abab(c,b,k,j) t2_2_abab(a,d,i,l)
        //            += +2.00 <k,b||c,d>_abab t2_1_abab(a,d,i,j) t1_1_aa(c,k)
        //            += -1.00 <l,k||c,d>_aaaa t2_aaaa(c,d,i,k) t2_2_abab(a,b,l,j)
        //            += -1.00 <l,k||d,c>_abab t2_1_abab(a,c,i,j) t2_1_abab(d,b,l,k)
        //            += -1.00 <k,l||d,c>_abab t2_1_abab(a,c,i,j) t2_1_abab(d,b,k,l)
        //            += -1.00 <l,k||c,d>_bbbb t2_bbbb(c,d,j,k) t2_2_abab(a,b,i,l)
        //            += +2.00 <l,k||d,c>_abab t2_1_abab(a,c,i,k) t2_1_abab(d,b,l,j)
        //            += +2.00 <k,b||c,d>_bbbb t2_1_abab(a,d,i,j) t1_1_bb(c,k)
        //            += +2.00 <k,l||d,c>_abab t2_1_abab(a,c,k,j) t2_1_abab(d,b,i,l)
        //            += +1.00 <l,k||c,d>_bbbb t2_1_abab(a,c,i,j) t2_1_bbbb(d,b,l,k)
        //            += +0.50 <l,k||c,d>_abab t2_abab(c,d,i,j) t2_2_abab(a,b,l,k)
        //            += +0.50 <k,l||c,d>_abab t2_abab(c,d,i,j) t2_2_abab(a,b,k,l)
        //            += +0.50 <l,k||d,c>_abab t2_abab(d,c,i,j) t2_2_abab(a,b,l,k)
        //            += +0.50 <k,l||d,c>_abab t2_abab(d,c,i,j) t2_2_abab(a,b,k,l)
        //            += -2.00 <a,k||c,d>_abab t2_1_bbbb(d,b,j,k) t1_1_aa(c,i)
        //            += -1.00 <k,l||c,d>_abab t2_1_abab(a,b,i,l) t2_1_abab(c,d,k,j)
        //            += -1.00 <k,l||d,c>_abab t2_1_abab(a,b,i,l) t2_1_abab(d,c,k,j)
        //            += +2.00 <k,l||c,d>_abab t2_1_aaaa(c,a,i,k) t2_1_bbbb(d,b,j,l)
        //            += -1.00 <l,k||c,d>_aaaa t2_1_aaaa(c,a,l,k) t2_1_abab(d,b,i,j)
        //            += +6.00 d-_bb(k,c) t1_1_bb(c,j) t2_2_abab(a,b,i,k)
        //            += +2.00 <l,k||c,d>_bbbb t2_1_abab(a,c,i,k) t2_1_bbbb(d,b,j,l)
        //            += -2.00 <l,k||i,c>_abab t2_1_abab(a,b,l,j) t1_1_bb(c,k)
        //            += +2.00 <l,k||c,j>_bbbb t2_1_abab(a,b,i,l) t1_1_bb(c,k)
        //            += +2.00 <k,a||c,d>_aaaa t2_1_abab(d,b,i,j) t1_1_aa(c,k)
        //            += -1.00 <l,k||c,d>_aaaa t2_1_abab(a,b,l,j) t2_1_aaaa(c,d,i,k)
        //            += +2.00 <l,k||d,c>_abab t2_bbbb(c,b,j,k) t2_2_aaaa(d,a,i,l)
        //            += -2.00 <k,a||c,d>_aaaa t2_1_abab(d,b,k,j) t1_1_aa(c,i)
        //            += -2.00 d-_aa(a,c) t1_1_bb(b,j) t1_2_aa(c,i)
        //            += -2.00 f_aa(k,c) t2_1_abab(a,b,k,j) t1_1_aa(c,i)
        //            += +6.00 d-_bb(k,c) t2_1_abab(a,b,i,k) t1_2_bb(c,j)
        //            += -2.00 <a,k||d,c>_abab t2_bbbb(c,b,j,k) t1_2_aa(d,i)
        //            += +2.00 <k,b||c,d>_bbbb t2_abab(a,c,i,k) t1_2_bb(d,j)
        //            += +2.00 <l,k||d,c>_abab t2_abab(a,c,i,k) t2_2_abab(d,b,l,j)
        //            += -2.00 <k,b||c,d>_abab t2_1_abab(a,d,k,j) t1_1_aa(c,i)
        //            += -1.00 <l,k||c,d>_bbbb t2_1_abab(a,b,i,l) t2_1_bbbb(c,d,j,k)
        //            += -1.00 <l,k||c,d>_abab t2_abab(c,d,i,k) t2_2_abab(a,b,l,j)
        //            += -1.00 <l,k||d,c>_abab t2_abab(d,c,i,k) t2_2_abab(a,b,l,j)
        //            += -2.00 <k,b||c,d>_bbbb t2_1_abab(a,d,i,k) t1_1_bb(c,j)
        //            += -2.00 <k,l||c,j>_abab t2_1_abab(a,b,i,l) t1_1_aa(c,k)
        //            += +2.00 <l,k||c,d>_bbbb t2_bbbb(c,b,j,k) t2_2_abab(a,d,i,l)
        //            += -1.00 <l,k||c,d>_bbbb t2_bbbb(c,b,l,k) t2_2_abab(a,d,i,j)
        //            += +2.00 <k,l||c,d>_abab t2_aaaa(c,a,i,k) t1_1_bb(b,l) t1_1_bb(d,j)
        //            += +2.00 <l,k||c,d>_bbbb t2_abab(a,c,i,k) t1_1_bb(b,l) t1_1_bb(d,j)
        //            += -2.00 <l,k||c,d>_bbbb t2_abab(a,b,i,k) t1_1_bb(c,l) t1_1_bb(d,j)
        //            += +2.00 d-_bb(k,c) t2_abab(a,b,i,k) t0_1 t1_2_bb(c,j)
        //            += +4.00 d-_bb(k,c) t2_abab(a,b,i,k) t1_1_bb(c,j) t0_2
        //            += +2.00 <l,k||c,d>_abab t2_abab(c,b,i,k) t1_1_aa(a,l) t1_1_bb(d,j)
        //            += -2.00 <l,k||c,d>_aaaa t2_abab(c,b,i,j) t1_1_aa(a,l) t1_1_aa(d,k)
        //            += -2.00 <l,k||c,d>_abab t2_abab(c,b,i,j) t1_1_aa(a,l) t1_1_bb(d,k)
        //            += +2.00 <l,k||c,d>_aaaa t2_abab(c,b,k,j) t1_1_aa(a,l) t1_1_aa(d,i)
        //            += +2.00 d-_aa(k,c) t0_1 t2_1_abab(a,b,k,j) t1_1_aa(c,i)
        //            += +2.00 d-_bb(k,c) t0_1 t2_1_abab(a,b,i,k) t1_1_bb(c,j)
        //            += +2.00 d-_bb(k,c) t1_1_aa(a,i) t2_2_bbbb(c,b,j,k)
        //            += -2.00 d-_aa(k,c) t1_1_aa(a,i) t2_2_abab(c,b,k,j)
        //            += +2.00 d-_aa(k,c) t1_1_bb(b,j) t2_2_aaaa(c,a,i,k)
        //            += -2.00 d-_bb(k,c) t1_1_bb(b,j) t2_2_abab(a,c,i,k)
        //            += -2.00 <k,b||d,c>_abab t2_abab(a,c,k,j) t1_2_aa(d,i)
        //            += +2.00 <k,l||d,c>_abab t2_abab(a,c,k,j) t2_2_abab(d,b,i,l)
        //            += -2.00 <k,l||i,c>_abab t2_abab(a,b,k,j) t1_2_bb(c,l)
        //            += -2.00 <l,k||c,i>_aaaa t2_abab(a,b,k,j) t1_2_aa(c,l)
        //            += +4.00 d-_bb(k,c) t2_1_bbbb(c,b,j,k) t1_2_aa(a,i)
        //            += -4.00 d-_aa(k,c) t2_1_abab(c,b,k,j) t1_2_aa(a,i)
        //            += +4.00 d-_bb(k,c) t2_abab(a,c,i,j) t1_1_bb(b,k) t0_2
        //            += +2.00 d-_bb(k,c) t2_abab(a,c,i,j) t0_1 t1_2_bb(b,k)
        //            += +4.00 d-_bb(k,j) t2_1_abab(a,b,i,k) t0_2
        //            += -1.00 <l,k||c,d>_abab t2_abab(c,b,i,j) t2_2_abab(a,d,l,k)
        //            += -1.00 <k,l||c,d>_abab t2_abab(c,b,i,j) t2_2_abab(a,d,k,l)
        //            += +1.00 <l,k||c,d>_aaaa t2_abab(c,b,i,j) t2_2_aaaa(d,a,l,k)
        //            += -2.00 f_bb(k,c) t2_abab(a,b,i,k) t1_2_bb(c,j)
        //            += +1.00 <l,k||c,d>_bbbb t2_abab(a,b,i,k) t2_2_bbbb(c,d,j,l)
        //            += -1.00 <l,k||c,d>_abab t2_abab(a,b,i,k) t2_2_abab(c,d,l,j)
        //            += -1.00 <l,k||d,c>_abab t2_abab(a,b,i,k) t2_2_abab(d,c,l,j)
        //            += -2.00 d-_aa(a,c) t0_1 t2_2_abab(c,b,i,j)
        //            += -4.00 d-_aa(a,c) t1_1_aa(c,i) t1_2_bb(b,j)
        //            += +4.00 d-_aa(k,i) t1_1_aa(a,k) t1_2_bb(b,j)
        //            += -2.00 <l,k||c,j>_abab t2_abab(a,b,i,k) t1_2_aa(c,l)
        //            += -2.00 <l,k||c,j>_bbbb t2_abab(a,b,i,k) t1_2_bb(c,l)
        //            += -2.00 f_bb(k,c) t2_abab(a,c,i,j) t1_2_bb(b,k)
        //            += -1.00 <a,k||c,d>_abab t2_abab(c,d,i,j) t1_2_bb(b,k)
        //            += -1.00 <a,k||d,c>_abab t2_abab(d,c,i,j) t1_2_bb(b,k)
        //            += -4.00 d-_aa(a,c) t2_1_abab(c,b,i,j) t0_2
        //            += +2.00 <l,k||i,c>_abab t2_1_abab(a,c,l,j) t1_1_bb(b,k)
        //            += -2.00 <a,k||i,c>_abab t1_1_bb(b,k) t1_1_bb(c,j)
        //            += -1.00 <a,k||c,d>_abab t1_1_bb(b,k) t2_1_abab(c,d,i,j)
        //            += -1.00 <a,k||d,c>_abab t1_1_bb(b,k) t2_1_abab(d,c,i,j)
        //            += +6.00 d-_bb(k,c) t1_1_bb(b,k) t2_2_abab(a,c,i,j)
        //            += -2.00 <l,k||c,j>_bbbb t2_1_abab(a,c,i,l) t1_1_bb(b,k)
        //            += +2.00 <l,k||c,j>_abab t2_1_aaaa(c,a,i,l) t1_1_bb(b,k)
        //            += -2.00 <a,k||c,j>_abab t1_1_bb(b,k) t1_1_aa(c,i)
        //            += -2.00 f_bb(k,c) t2_1_abab(a,c,i,j) t1_1_bb(b,k)
        //            += -4.00 d-_bb(b,c) t2_1_abab(a,c,i,j) t0_2
        //            += +2.00 <l,k||c,d>_aaaa t2_aaaa(c,a,i,k) t2_2_abab(d,b,l,j)
        //            += +2.00 d-_aa(k,i) t0_1 t2_2_abab(a,b,k,j)
        //            += +2.00 <l,k||i,c>_abab t2_bbbb(c,b,j,k) t1_2_aa(a,l)
        //            += +2.00 <l,k||c,j>_abab t2_abab(c,b,i,k) t1_2_aa(a,l)
        //            += +2.00 <l,k||c,i>_aaaa t2_abab(c,b,k,j) t1_2_aa(a,l)
        //            += -2.00 <k,b||i,c>_abab t1_1_aa(a,k) t1_1_bb(c,j)
        //            += -2.00 f_aa(k,c) t1_1_aa(a,k) t2_1_abab(c,b,i,j)
        //            += +2.00 <k,l||c,j>_abab t1_1_aa(a,k) t2_1_abab(c,b,i,l)
        //            += +2.00 <k,l||i,j>_abab t1_1_aa(a,k) t1_1_bb(b,l)
        //            += +6.00 d-_aa(k,c) t1_1_aa(a,k) t2_2_abab(c,b,i,j)
        //            += +2.00 <k,l||i,c>_abab t1_1_aa(a,k) t2_1_bbbb(c,b,j,l)
        //            += -2.00 <l,k||c,i>_aaaa t1_1_aa(a,k) t2_1_abab(c,b,l,j)
        //            += -2.00 <k,b||c,j>_abab t1_1_aa(a,k) t1_1_aa(c,i)
        //            += -1.00 <k,b||c,d>_abab t1_1_aa(a,k) t2_1_abab(c,d,i,j)
        //            += -1.00 <k,b||d,c>_abab t1_1_aa(a,k) t2_1_abab(d,c,i,j)
        //            += +4.00 d-_aa(k,i) t2_1_abab(a,b,k,j) t0_2
        //            += +2.00 <k,a||c,i>_aaaa t2_2_abab(c,b,k,j)
        //            += +2.00 <a,b||c,j>_abab t1_2_aa(c,i)
        //            += -2.00 <a,k||c,d>_abab t2_abab(c,b,i,k) t1_2_bb(d,j)
        //            += +4.00 t2_2_abab(a,b,i,j) w0
        //            += -2.00 f_aa(k,i) t2_2_abab(a,b,k,j)
        //            += -2.00 <k,b||i,j>_abab t1_2_aa(a,k)
        //            += -2.00 f_bb(k,j) t2_2_abab(a,b,i,k)
        //            += -2.00 <a,k||c,j>_abab t2_2_abab(c,b,i,k)
        //            += +1.00 <l,k||i,j>_abab t2_2_abab(a,b,l,k)
        //            += +1.00 <k,l||i,j>_abab t2_2_abab(a,b,k,l)
        //            += +2.00 <a,b||d,c>_abab t1_1_bb(c,j) t1_1_aa(d,i)
        //            += +1.00 <a,b||c,d>_abab t2_2_abab(c,d,i,j)
        //            += +1.00 <a,b||d,c>_abab t2_2_abab(d,c,i,j)
        //            += +2.00 f_bb(b,c) t2_2_abab(a,c,i,j)
        //            += -2.00 <k,b||c,j>_abab t2_2_aaaa(c,a,i,k)
        //            += +2.00 f_aa(a,c) t2_2_abab(c,b,i,j)
        //            += +2.00 <a,b||i,c>_abab t1_2_bb(c,j)
        //            += -2.00 <a,k||i,c>_abab t2_2_bbbb(c,b,j,k)
        //            += -4.00 d-_bb(k,c) t1_1_bb(c,k) t2_2_abab(a,b,i,j)
        //            += -4.00 d-_aa(k,c) t1_1_aa(c,k) t2_2_abab(a,b,i,j)
        //            += -2.00 <k,b||i,c>_abab t2_2_abab(a,c,k,j)
        //            += -2.00 <a,k||i,j>_abab t1_2_bb(b,k)
        //            += -2.00 d-_aa(k,c) t2_1_abab(a,b,i,j) t1_2_aa(c,k)
        //            += -2.00 d-_bb(k,c) t2_1_abab(a,b,i,j) t1_2_bb(c,k)
        //            += +2.00 <k,b||c,j>_bbbb t2_2_abab(a,c,i,k)
        //            += +2.00 <a,k||c,d>_abab t2_abab(c,b,i,j) t1_2_bb(d,k)
        //            += -2.00 <k,a||c,d>_aaaa t2_abab(c,b,i,j) t1_2_aa(d,k)
        //            += +2.00 d-_aa(k,c) t0_1 t1_1_aa(a,k) t2_1_abab(c,b,i,j)
        //            += +1.00 <l,k||i,c>_abab t2_1_abab(a,b,l,k) t1_1_bb(c,j)
        //            += +1.00 <k,l||i,c>_abab t2_1_abab(a,b,k,l) t1_1_bb(c,j)
        //            += +1.00 <l,k||c,j>_abab t2_1_abab(a,b,l,k) t1_1_aa(c,i)
        //            += +1.00 <k,l||c,j>_abab t2_1_abab(a,b,k,l) t1_1_aa(c,i)
        //            += +0.50 <l,k||c,d>_abab t2_1_abab(a,b,l,k) t2_1_abab(c,d,i,j)
        //            += +0.50 <l,k||d,c>_abab t2_1_abab(a,b,l,k) t2_1_abab(d,c,i,j)
        //            += +0.50 <k,l||c,d>_abab t2_1_abab(a,b,k,l) t2_1_abab(c,d,i,j)
        //            += +0.50 <k,l||d,c>_abab t2_1_abab(a,b,k,l) t2_1_abab(d,c,i,j)
        //            += +1.00 <l,k||c,j>_abab t2_abab(a,b,l,k) t1_2_aa(c,i)
        //            += +1.00 <k,l||c,j>_abab t2_abab(a,b,k,l) t1_2_aa(c,i)
        //            += +0.50 <l,k||c,d>_abab t2_abab(a,b,l,k) t2_2_abab(c,d,i,j)
        //            += +0.50 <l,k||d,c>_abab t2_abab(a,b,l,k) t2_2_abab(d,c,i,j)
        //            += +0.50 <k,l||c,d>_abab t2_abab(a,b,k,l) t2_2_abab(c,d,i,j)
        //            += +0.50 <k,l||d,c>_abab t2_abab(a,b,k,l) t2_2_abab(d,c,i,j)
        //            += +1.00 <l,k||i,c>_abab t2_abab(a,b,l,k) t1_2_bb(c,j)
        //            += +1.00 <k,l||i,c>_abab t2_abab(a,b,k,l) t1_2_bb(c,j)
        //            += +2.00 <l,k||c,d>_bbbb t2_abab(a,c,i,k) t2_2_bbbb(d,b,j,l)
        //            += -4.00 d-_bb(k,c) t2_1_abab(a,c,i,k) t1_2_bb(b,j)
        //            += +4.00 d-_aa(k,c) t2_1_aaaa(c,a,i,k) t1_2_bb(b,j)
        //            += +4.00 d-_bb(k,j) t1_1_bb(b,k) t1_2_aa(a,i)
        //            += -4.00 d-_bb(b,c) t1_1_bb(c,j) t1_2_aa(a,i)
        //            += +2.00 <k,b||d,c>_abab t2_abab(a,c,i,j) t1_2_aa(d,k)
        //            += -2.00 <k,b||c,d>_bbbb t2_abab(a,c,i,j) t1_2_bb(d,k)
        //            += +2.00 d-_aa(k,c) t2_abab(c,b,i,j) t0_1 t1_2_aa(a,k)
        //            += +4.00 d-_aa(k,c) t2_abab(c,b,i,j) t1_1_aa(a,k) t0_2
        //            += -2.00 f_aa(k,c) t2_abab(a,b,k,j) t1_2_aa(c,i)
        //            += -1.00 <k,l||c,d>_abab t2_abab(a,b,k,j) t2_2_abab(c,d,i,l)
        //            += -1.00 <k,l||d,c>_abab t2_abab(a,b,k,j) t2_2_abab(d,c,i,l)
        //            += +1.00 <l,k||c,d>_aaaa t2_abab(a,b,k,j) t2_2_aaaa(c,d,i,l)
        //            += +2.00 <l,k||c,j>_bbbb t2_abab(a,c,i,k) t1_2_bb(b,l)
        //            += +2.00 <k,l||i,c>_abab t2_abab(a,c,k,j) t1_2_bb(b,l)
        //            += +2.00 <k,l||c,j>_abab t2_aaaa(c,a,i,k) t1_2_bb(b,l)
        //            += +2.00 d-_bb(k,c) t0_1 t2_1_abab(a,c,i,j) t1_1_bb(b,k)
        //            += -1.00 <l,k||d,c>_abab t2_abab(a,c,i,j) t2_2_abab(d,b,l,k)
        //            += -1.00 <k,l||d,c>_abab t2_abab(a,c,i,j) t2_2_abab(d,b,k,l)
        //            += +1.00 <l,k||c,d>_bbbb t2_abab(a,c,i,j) t2_2_bbbb(d,b,l,k)
        //            += -2.00 d-_bb(b,c) t0_1 t2_2_abab(a,c,i,j)
        //            += +2.00 d-_bb(k,j) t0_1 t2_2_abab(a,b,i,k)
        //            += -1.00 <k,b||c,d>_abab t2_abab(c,d,i,j) t1_2_aa(a,k)
        //            += -1.00 <k,b||d,c>_abab t2_abab(d,c,i,j) t1_2_aa(a,k)
        //            += -2.00 f_aa(k,c) t2_abab(c,b,i,j) t1_2_aa(a,k)
        //            += +4.00 d-_aa(k,c) t2_abab(a,b,k,j) t1_1_aa(c,i) t0_2
        //            += -2.00 <k,l||d,c>_abab t2_abab(a,b,k,j) t1_1_bb(c,l) t1_1_aa(d,i)
        //            += -2.00 <l,k||c,d>_aaaa t2_abab(a,b,k,j) t1_1_aa(c,l) t1_1_aa(d,i)
        //            += +2.00 d-_aa(k,c) t2_abab(a,b,k,j) t0_1 t1_2_aa(c,i)
        //            += -1.00 <k,l||c,d>_abab t2_abab(c,d,k,j) t2_2_abab(a,b,i,l)
        //            += -1.00 <k,l||d,c>_abab t2_abab(d,c,k,j) t2_2_abab(a,b,i,l)
        //            += -1.00 <l,k||c,d>_abab t2_abab(c,b,l,k) t2_2_abab(a,d,i,j)
        //            += -1.00 <k,l||c,d>_abab t2_abab(c,b,k,l) t2_2_abab(a,d,i,j)
        //            += -1.00 <l,k||d,c>_abab t2_abab(a,c,l,k) t2_2_abab(d,b,i,j)
        //            += -1.00 <k,l||d,c>_abab t2_abab(a,c,k,l) t2_2_abab(d,b,i,j)
        //            += +2.00 <l,k||c,d>_abab t2_abab(c,b,i,k) t2_2_abab(a,d,l,j)
        //            += -2.00 <k,b||d,c>_abab t2_1_aaaa(d,a,i,k) t1_1_bb(c,j)
        //            += +6.00 d-_bb(k,c) t2_1_abab(a,c,i,j) t1_2_bb(b,k)
        //            += +2.00 d-_aa(k,i) t1_1_bb(b,j) t1_2_aa(a,k)
        rt2_2_abab("a,b,i,j")  = 2.00 * tmps_["191_abab_vvoo"]("a,b,i,j");
        tmps_["191_abab_vvoo"].~TArrayD();

        // flops: o2v2  = o2v3 o3v2 o2v2 o2v3 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["192_abab_vvoo"]("a,b,i,j")  = -1.00 * dp["aa_vv"]("a,c") * t2_1["abab"]("c,b,i,j");
        tmps_["192_abab_vvoo"]("a,b,i,j") += dp["bb_oo"]("l,j") * t2_1["abab"]("a,b,i,l");
        tmps_["192_abab_vvoo"]("a,b,i,j") -= dp["bb_vv"]("b,d") * t2_1["abab"]("a,d,i,j");
        tmps_["192_abab_vvoo"]("a,b,i,j") += dp["aa_oo"]("k,i") * t2_1["abab"]("a,b,k,j");
        tmps_["192_abab_vvoo"]("a,b,i,j") += tmps_["9_aa_oo"]("k,i") * t2["abab"]("a,b,k,j");
        tmps_["192_abab_vvoo"]("a,b,i,j") += tmps_["24_bb_oo"]("l,j") * t2["abab"]("a,b,i,l");
        tmps_["192_abab_vvoo"]("a,b,i,j") += tmps_["91_baab_ovoo"]("l,a,i,j") * t1_1["bb"]("b,l");
        tmps_["192_abab_vvoo"]("a,b,i,j") += tmps_["93_abab_ovoo"]("k,b,i,j") * t1_1["aa"]("a,k");
        tmps_["93_abab_ovoo"].~TArrayD();
        tmps_["91_baab_ovoo"].~TArrayD();
        tmps_["9_aa_oo"].~TArrayD();

        // rt2_abab  = +1.00 d-_bb(k,c) t2_abab(a,c,i,j) t1_1_bb(b,k)
        //          += +1.00 d-_aa(k,c) t2_abab(a,b,k,j) t1_1_aa(c,i)
        //          += +1.00 d-_bb(k,j) t2_1_abab(a,b,i,k)
        //          += -1.00 d-_aa(a,c) t2_1_abab(c,b,i,j)
        //          += -1.00 d-_bb(b,c) t2_1_abab(a,c,i,j)
        //          += +1.00 d-_aa(k,i) t2_1_abab(a,b,k,j)
        //          += +1.00 d-_bb(k,c) t2_abab(a,b,i,k) t1_1_bb(c,j)
        //          += +1.00 d-_aa(k,c) t2_abab(c,b,i,j) t1_1_aa(a,k)
        rt2_abab("a,b,i,j")  = tmps_["192_abab_vvoo"]("a,b,i,j");

        // rt2_2_abab += +2.00 d+_bb(k,c) t2_abab(a,c,i,j) t1_1_bb(b,k)
        //            += +2.00 d+_aa(k,c) t2_abab(a,b,k,j) t1_1_aa(c,i)
        //            += +2.00 d+_bb(k,j) t2_1_abab(a,b,i,k)
        //            += -2.00 d+_aa(a,c) t2_1_abab(c,b,i,j)
        //            += -2.00 d+_bb(b,c) t2_1_abab(a,c,i,j)
        //            += +2.00 d+_aa(k,i) t2_1_abab(a,b,k,j)
        //            += +2.00 d+_bb(k,c) t2_abab(a,b,i,k) t1_1_bb(c,j)
        //            += +2.00 d+_aa(k,c) t2_abab(c,b,i,j) t1_1_aa(a,k)
        rt2_2_abab("a,b,i,j") += 2.00 * tmps_["192_abab_vvoo"]("a,b,i,j");
        tmps_["192_abab_vvoo"].~TArrayD();

        // flops: o4v0  = o4v1
        //  mems: o4v0  = o4v0
        tmps_["197_bbbb_oooo"]("i,j,k,l")  = t1_1["bb"]("a,i") * tmps_["168_bbbb_oovo"]("j,k,a,l");
        tmps_["168_bbbb_oovo"].~TArrayD();

        // flops: o3v1  = o4v1
        //  mems: o3v1  = o3v1
        tmps_["196_bbbb_vooo"]("a,i,j,k")  = tmps_["117_bbbb_oooo"]("i,l,j,k") * t1_1["bb"]("a,l");

        // flops: o3v1  = o4v1
        //  mems: o3v1  = o3v1
        tmps_["195_bbbb_ooov"]("i,j,k,a")  = eri["bbbb_oooo"]("l,i,j,k") * t1_1["bb"]("a,l");

        // flops: o4v0  = o4v2
        //  mems: o4v0  = o4v0
        tmps_["194_bbbb_oooo"]("i,j,k,l")  = eri["bbbb_oovv"]("i,j,a,b") * t2_2["bbbb"]("a,b,k,l");

        // flops: o0v2  = o2v3
        //  mems: o0v2  = o0v2
        tmps_["193_bb_vv"]("a,b")  = eri["bbbb_oovv"]("i,j,c,a") * t2_1["bbbb"]("c,b,j,i");

        // flops: o2v2  = o4v2 o2v3 o2v2 o4v2 o2v2 o2v3 o2v2 o4v2 o2v2 o1v4 o2v3 o2v2 o2v2 o4v2 o2v2 o2v2 o2v2 o2v4 o2v2 o2v2 o2v2 o3v2 o2v2 o3v2 o2v2 o2v2 o4v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o1v3 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["198_bbbb_vvoo"]("a,b,i,j")  = -0.50 * tmps_["194_bbbb_oooo"]("l,k,i,j") * t2["bbbb"]("a,b,k,l");
        tmps_["198_bbbb_vvoo"]("a,b,i,j") += tmps_["193_bb_vv"]("d,a") * t2_1["bbbb"]("d,b,i,j");
        tmps_["198_bbbb_vvoo"]("a,b,i,j") -= 0.50 * tmps_["118_bbbb_oooo"]("l,k,i,j") * t2_1["bbbb"]("a,b,k,l");
        tmps_["198_bbbb_vvoo"]("a,b,i,j") += tmps_["52_bb_vv"]("c,b") * t2_1["bbbb"]("c,a,i,j");
        tmps_["198_bbbb_vvoo"]("a,b,i,j") -= 0.50 * tmps_["117_bbbb_oooo"]("l,k,i,j") * t2_2["bbbb"]("a,b,k,l");
        tmps_["198_bbbb_vvoo"]("a,b,i,j") -= 2.00 * eri["bbbb_vvvv"]("a,b,c,d") * t1_1["bb"]("d,i") * t1_1["bb"]("c,j");
        tmps_["198_bbbb_vvoo"]("a,b,i,j") -= 2.00 * scalars_["6"] * t2_1["bbbb"]("a,b,i,j");
        tmps_["198_bbbb_vvoo"]("a,b,i,j") -= eri["bbbb_oooo"]("l,k,i,j") * t2_2["bbbb"]("a,b,k,l");
        tmps_["198_bbbb_vvoo"]("a,b,i,j") += 4.00 * w0 * t2_2["bbbb"]("a,b,i,j");
        tmps_["198_bbbb_vvoo"]("a,b,i,j") += eri["bbbb_vvvv"]("a,b,c,d") * t2_2["bbbb"]("c,d,i,j");
        tmps_["198_bbbb_vvoo"]("a,b,i,j") -= 4.00 * scalars_["1"] * t2_2["bbbb"]("a,b,i,j");
        tmps_["198_bbbb_vvoo"]("a,b,i,j") += 2.00 * t1_1["bb"]("b,k") * tmps_["195_bbbb_ooov"]("k,i,j,a");
        tmps_["198_bbbb_vvoo"]("a,b,i,j") += t1_1["bb"]("a,l") * tmps_["196_bbbb_vooo"]("b,l,i,j");
        tmps_["198_bbbb_vvoo"]("a,b,i,j") += tmps_["197_bbbb_oooo"]("j,l,k,i") * t2["bbbb"]("a,b,k,l");
        tmps_["197_bbbb_oooo"].~TArrayD();
        tmps_["196_bbbb_vooo"].~TArrayD();
        tmps_["195_bbbb_ooov"].~TArrayD();
        tmps_["194_bbbb_oooo"].~TArrayD();
        tmps_["193_bb_vv"].~TArrayD();
        tmps_["118_bbbb_oooo"].~TArrayD();
        tmps_["52_bb_vv"].~TArrayD();

        // rt2_2_bbbb  = -1.00 <l,k||c,d>_bbbb t2_bbbb(a,b,l,k) t1_1_bb(c,j) t1_1_bb(d,i)
        //            += -1.00 <l,k||c,d>_bbbb t2_1_bbbb(c,a,l,k) t2_1_bbbb(d,b,i,j)
        //            += +0.50 <l,k||c,d>_bbbb t2_bbbb(a,b,l,k) t2_2_bbbb(c,d,i,j)
        //            += +0.50 <l,k||c,d>_bbbb t2_1_bbbb(a,b,l,k) t2_1_bbbb(c,d,i,j)
        //            += -1.00 <l,k||c,d>_bbbb t2_1_bbbb(c,a,i,j) t2_1_bbbb(d,b,l,k)
        //            += +0.50 <l,k||c,d>_bbbb t2_bbbb(c,d,i,j) t2_2_bbbb(a,b,l,k)
        //            += -2.00 <l,k||i,j>_bbbb t1_1_bb(a,k) t1_1_bb(b,l)
        //            += -2.00 d-_aa(k,c) t2_1_bbbb(a,b,i,j) t1_2_aa(c,k)
        //            += -2.00 d-_bb(k,c) t2_1_bbbb(a,b,i,j) t1_2_bb(c,k)
        //            += -2.00 <a,b||c,d>_bbbb t1_1_bb(c,j) t1_1_bb(d,i)
        //            += +4.00 t2_2_bbbb(a,b,i,j) w0
        //            += +1.00 <l,k||i,j>_bbbb t2_2_bbbb(a,b,l,k)
        //            += -4.00 d-_bb(k,c) t1_1_bb(c,k) t2_2_bbbb(a,b,i,j)
        //            += -4.00 d-_aa(k,c) t1_1_aa(c,k) t2_2_bbbb(a,b,i,j)
        //            += +1.00 <a,b||c,d>_bbbb t2_2_bbbb(c,d,i,j)
        //            += -1.00 <l,k||c,d>_bbbb t2_bbbb(c,d,i,j) t1_1_bb(a,k) t1_1_bb(b,l)
        rt2_2_bbbb("a,b,i,j")  = tmps_["198_bbbb_vvoo"]("a,b,i,j");
        tmps_["198_bbbb_vvoo"].~TArrayD();

        // flops: o3v1  = o2v2 o3v2
        //  mems: o3v1  = o2v2 o3v1
        tmps_["201_bbbb_oovo"]("i,j,a,k")  = t1_1["bb"]("b,i") * (-1.00 * tmps_["19_bbbb_ovvo"]("j,b,a,k") + tmps_["18_bbbb_ovvo"]("j,b,a,k"));

        // flops: o3v1  = o4v2 o4v2 o3v1
        //  mems: o3v1  = o3v1 o3v1 o3v1
        tmps_["200_bbbb_oovo"]("i,j,a,k")  = eri["abab_oovo"]("l,i,b,j") * t2_1["abab"]("b,a,l,k");
        tmps_["200_bbbb_oovo"]("i,j,a,k") += eri["bbbb_oovo"]("i,m,c,j") * t2_1["bbbb"]("c,a,k,m");

        // flops: o3v1  = o3v2
        //  mems: o3v1  = o3v1
        tmps_["199_bbbb_vooo"]("a,i,j,k")  = eri["bbbb_vovo"]("a,i,b,j") * t1_1["bb"]("b,k");

        // flops: o2v2  = o2v2 o3v3 o3v3 o2v2 o2v2 o2v2 o3v2 o2v2 o2v2 o2v2 o3v2 o2v2 o2v2 o2v2 o3v2 o2v2 o3v3 o2v2 o3v3 o2v2 o2v2 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v2 o2v2 o2v2 o2v2 o3v3 o2v2 o3v3 o2v2 o2v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["202_bbbb_vovo"]("a,i,b,j")  = -1.00 * -1.00 * (tmps_["174_bb_vo"]("a,i") * t1_1["bb"]("b,j") + -1.00 * eri["baab_vovo"]("b,m,e,j") * t2_2["abab"]("e,a,m,i") + eri["bbbb_vovo"]("b,k,c,j") * t2_2["bbbb"]("c,a,i,k"));
        tmps_["202_bbbb_vovo"]("a,i,b,j") -= t1_1["bb"]("b,k") * tmps_["200_bbbb_oovo"]("k,j,a,i");
        tmps_["202_bbbb_vovo"]("a,i,b,j") += 2.00 * t1_2["bb"]("a,i") * tmps_["104_bb_vo"]("b,j");
        tmps_["202_bbbb_vovo"]("a,i,b,j") -= t1_2["bb"]("a,l") * tmps_["110_bbbb_oovo"]("l,j,b,i");
        tmps_["202_bbbb_vovo"]("a,i,b,j") += 2.00 * t1_2["bb"]("a,i") * tmps_["103_bb_vo"]("b,j");
        tmps_["202_bbbb_vovo"]("a,i,b,j") -= t1_1["bb"]("a,l") * tmps_["201_bbbb_oovo"]("i,l,b,j");
        tmps_["202_bbbb_vovo"]("a,i,b,j") += tmps_["17_aabb_ovvo"]("n,f,b,j") * t2_2["abab"]("f,a,n,i");
        tmps_["202_bbbb_vovo"]("a,i,b,j") -= tmps_["162_bbbb_ovvo"]("k,c,a,i") * t2["bbbb"]("c,b,j,k");
        tmps_["202_bbbb_vovo"]("a,i,b,j") += tmps_["170_bb_vo"]("b,i") * t1_1["bb"]("a,j");
        tmps_["202_bbbb_vovo"]("a,i,b,j") -= tmps_["51_baab_vovo"]("b,m,f,j") * t2_1["abab"]("f,a,m,i");
        tmps_["202_bbbb_vovo"]("a,i,b,j") += tmps_["163_baab_vovo"]("b,m,e,i") * t2["abab"]("e,a,m,j");
        tmps_["202_bbbb_vovo"]("a,i,b,j") -= tmps_["18_bbbb_ovvo"]("l,d,b,j") * t2_2["bbbb"]("d,a,i,l");
        tmps_["202_bbbb_vovo"]("a,i,b,j") -= tmps_["164_bbbb_vovo"]("b,k,d,j") * t2_1["bbbb"]("d,a,i,k");
        tmps_["202_bbbb_vovo"]("a,i,b,j") += tmps_["199_bbbb_vooo"]("b,k,j,i") * t1_1["bb"]("a,k");
        tmps_["202_bbbb_vovo"]("a,i,b,j") -= tmps_["171_bb_ov"]("j,a") * t1_1["bb"]("b,i");
        tmps_["202_bbbb_vovo"]("a,i,b,j") -= tmps_["165_bbbb_vovo"]("b,k,c,i") * t2["bbbb"]("c,a,j,k");
        tmps_["202_bbbb_vovo"]("a,i,b,j") += tmps_["19_bbbb_ovvo"]("l,d,b,j") * t2_2["bbbb"]("d,a,i,l");
        tmps_["202_bbbb_vovo"]("a,i,b,j") += tmps_["185_bb_vo"]("b,j") * t1_1["bb"]("a,i");
        tmps_["201_bbbb_oovo"].~TArrayD();
        tmps_["200_bbbb_oovo"].~TArrayD();
        tmps_["199_bbbb_vooo"].~TArrayD();
        tmps_["185_bb_vo"].~TArrayD();
        tmps_["174_bb_vo"].~TArrayD();
        tmps_["171_bb_ov"].~TArrayD();
        tmps_["170_bb_vo"].~TArrayD();
        tmps_["165_bbbb_vovo"].~TArrayD();
        tmps_["164_bbbb_vovo"].~TArrayD();
        tmps_["163_baab_vovo"].~TArrayD();
        tmps_["162_bbbb_ovvo"].~TArrayD();
        tmps_["110_bbbb_oovo"].~TArrayD();
        tmps_["104_bb_vo"].~TArrayD();
        tmps_["103_bb_vo"].~TArrayD();
        tmps_["51_baab_vovo"].~TArrayD();
}