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
void QED_CCSD_22::resid_22_3() {
    
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


        // flops: o1v1  = o1v2 o2v1 o1v1
        //  mems: o1v1  = o1v1 o1v1 o1v1
        tmps_["104_bb_vo"]("a,i")  = -1.00 * dp["bb_vv"]("a,b") * t1_1["bb"]("b,i");
        tmps_["104_bb_vo"]("a,i") += dp["bb_oo"]("j,i") * t1_1["bb"]("a,j");

        // flops: o1v1  = o2v2 o2v2 o1v1
        //  mems: o1v1  = o1v1 o1v1 o1v1
        tmps_["103_bb_vo"]("a,i")  = -1.00 * dp["aa_ov"]("k,c") * t2_1["abab"]("c,a,k,i");
        tmps_["103_bb_vo"]("a,i") += dp["bb_ov"]("j,b") * t2_1["bbbb"]("b,a,i,j");

        // flops: o1v1  = o2v2 o2v2 o1v1
        //  mems: o1v1  = o1v1 o1v1 o1v1
        tmps_["102_bb_vo"]("a,i")  = -1.00 * dp["bb_ov"]("k,c") * t2["bbbb"]("c,a,i,k");
        tmps_["102_bb_vo"]("a,i") += dp["aa_ov"]("j,b") * t2["abab"]("b,a,j,i");

        // flops: o3v1  = o3v3 o3v2 o3v1
        //  mems: o3v1  = o3v1 o3v1 o3v1
        tmps_["101_abab_vooo"]("a,i,j,k")  = eri["abab_vovv"]("a,i,b,c") * t2["abab"]("b,c,j,k");
        tmps_["101_abab_vooo"]("a,i,j,k") += f["bb_ov"]("i,d") * t2["abab"]("a,d,j,k");

        // flops: o3v1  = o4v2 o4v2 o3v1 o4v2 o3v1
        //  mems: o3v1  = o3v1 o3v1 o3v1 o3v1 o3v1
        tmps_["100_abba_oovo"]("i,j,a,k")  = -1.00 * eri["abab_oovo"]("i,m,b,j") * t2["abab"]("b,a,k,m");
        tmps_["100_abba_oovo"]("i,j,a,k") += eri["abba_oovo"]("i,m,c,k") * t2["bbbb"]("c,a,j,m");
        tmps_["100_abba_oovo"]("i,j,a,k") += eri["aaaa_oovo"]("l,i,b,k") * t2["abab"]("b,a,l,j");

        // flops: o3v1  = o3v2 o3v3 o3v1
        //  mems: o3v1  = o3v1 o3v1 o3v1
        tmps_["99_abab_ovoo"]("i,a,j,k")  = -1.00 * f["aa_ov"]("i,b") * t2["abab"]("b,a,j,k");
        tmps_["99_abab_ovoo"]("i,a,j,k") += eri["baab_vovv"]("a,i,b,c") * t2["abab"]("b,c,j,k");

        // flops: o3v1  = o4v2 o4v2 o4v2 o3v1 o3v1
        //  mems: o3v1  = o3v1 o3v1 o3v1 o3v1 o3v1
        tmps_["98_bbaa_oovo"]("i,j,a,k")  = -1.00 * eri["abab_oovo"]("l,i,b,j") * t2["aaaa"]("b,a,k,l");
        tmps_["98_bbaa_oovo"]("i,j,a,k") += eri["abba_oovo"]("l,i,c,k") * t2["abab"]("a,c,l,j");
        tmps_["98_bbaa_oovo"]("i,j,a,k") += eri["bbbb_oovo"]("m,i,c,j") * t2["abab"]("a,c,k,m");

        // flops: o2v0  = o3v1
        //  mems: o2v0  = o2v0
        tmps_["97_bb_oo"]("i,j")  = eri["bbbb_oovo"]("i,k,a,j") * t1_1["bb"]("a,k");

        // flops: o0v2  = o1v3
        //  mems: o0v2  = o0v2
        tmps_["96_bb_vv"]("a,b")  = eri["bbbb_vovv"]("a,i,b,c") * t1_1["bb"]("c,i");

        // flops: o0v2  = o1v3
        //  mems: o0v2  = o0v2
        tmps_["95_bb_vv"]("a,b")  = eri["baab_vovv"]("a,i,c,b") * t1_1["aa"]("c,i");

        // flops: o3v1  = o3v2
        //  mems: o3v1  = o3v1
        tmps_["94_abab_ovoo"]("i,a,j,k")  = dp["aa_ov"]("i,b") * t2_1["abab"]("b,a,j,k");

        // flops: o3v1  = o3v2
        //  mems: o3v1  = o3v1
        tmps_["93_abab_ovoo"]("i,a,j,k")  = dp["aa_ov"]("i,b") * t2["abab"]("b,a,j,k");

        // flops: o3v1  = o3v2
        //  mems: o3v1  = o3v1
        tmps_["92_baab_ovoo"]("i,a,j,k")  = dp["bb_ov"]("i,b") * t2_1["abab"]("a,b,j,k");

        // flops: o3v1  = o3v2
        //  mems: o3v1  = o3v1
        tmps_["91_baab_ovoo"]("i,a,j,k")  = dp["bb_ov"]("i,b") * t2["abab"]("a,b,j,k");

        // flops: o2v2  = o2v3
        //  mems: o2v2  = o2v2
        tmps_["90_baba_vovo"]("a,i,b,j")  = eri["baab_vovv"]("a,i,c,b") * t1_1["aa"]("c,j");

        // flops: o2v2  = o2v3
        //  mems: o2v2  = o2v2
        tmps_["89_abab_vovo"]("a,i,b,j")  = eri["abab_vovv"]("a,i,b,c") * t1_1["bb"]("c,j");

        // flops: o4v0  = o4v2
        //  mems: o4v0  = o4v0
        tmps_["88_abab_oooo"]("i,j,k,l")  = eri["abab_oovv"]("i,j,a,b") * t2["abab"]("a,b,k,l");

        // flops: o2v2  = o3v3
        //  mems: o2v2  = o2v2
        tmps_["87_abba_ovvo"]("i,a,b,j")  = eri["abab_oovv"]("i,k,c,a") * t2_1["abab"]("c,b,j,k");

        // flops: o2v2  = o3v3
        //  mems: o2v2  = o2v2
        tmps_["86_abba_ovvo"]("i,a,b,j")  = eri["abab_oovv"]("i,k,c,a") * t2["abab"]("c,b,j,k");

        // flops: o2v2  = o3v2 o3v3 o2v2 o2v2 o3v2 o2v2 o4v2 o2v2 o3v3 o2v2 o2v4 o2v2 o2v3 o2v2 o2v3 o2v2 o3v2 o2v2 o3v3 o2v2 o2v3 o2v2 o3v3 o2v2 o3v2 o2v2 o3v2 o2v2 o2v3 o2v2 o3v2 o2v2 o2v2 o2v2 o3v3 o2v2 o3v3 o2v2 o2v2 o2v2 o2v3 o2v2 o3v2 o2v2 o2v3 o2v2 o2v2 o2v2 o3v2 o2v2 o2v3 o2v2 o2v3 o2v2 o3v2 o2v2 o2v2 o2v2 o3v2 o3v2 o2v2 o3v3 o2v2 o2v3 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o2v2 o2v2 o4v2 o2v2 o2v3 o2v2 o2v2 o2v2 o3v2 o2v2 o2v3 o2v2 o3v3 o2v2 o3v2 o2v2 o2v2 o2v2 o3v2 o2v2 o2v2 o2v3 o2v2 o2v2 o2v2 o2v2 o2v2 o3v2 o2v2 o2v2 o2v2 o3v2 o3v2 o2v2 o2v2 o2v2 o2v3 o3v2 o2v2 o2v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v2 o2v2 o4v2 o2v2 o3v2 o2v2 o2v3 o2v2 o2v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o2v3 o2v2 o3v2 o2v2 o3v3 o2v2 o2v3 o2v2 o3v2 o2v2 o3v2 o2v2 o3v3 o2v2 o3v2 o2v2 o3v3 o2v2 o3v2 o2v2 o3v2 o2v2 o2v3 o2v2 o3v2 o2v2 o3v3 o2v2 o3v3 o2v2 o2v3 o2v2 o3v2 o2v2 o3v2 o2v2 o3v3 o2v2 o3v2 o2v2 o2v3 o2v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v3 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o2v3 o2v2 o3v2 o2v2 o2v3 o2v2 o2v3 o2v2 o3v2 o2v2 o3v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["107_bbaa_vovo"]("a,i,b,j")  = -0.50 * tmps_["23_bb_oo"]("n,i") * t2_1["abab"]("b,a,j,n");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += eri["bbbb_vovo"]("a,k,c,i") * t2_1["abab"]("b,c,j,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += 2.00 * dp["aa_vo"]("b,j") * t1_2["bb"]("a,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += f["aa_oo"]("l,j") * t2_1["abab"]("b,a,l,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= eri["abab_oooo"]("m,k,j,i") * t2_1["abab"]("b,a,m,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += eri["abab_vovo"]("b,k,e,i") * t2_1["abab"]("e,a,j,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= eri["abab_vvvv"]("b,a,e,f") * t2_1["abab"]("e,f,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += 2.00 * dp["bb_vv"]("a,c") * t2_2["abab"]("b,c,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= eri["abab_vvvo"]("b,a,e,i") * t1_1["aa"]("e,j");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= 2.00 * dp["bb_oo"]("k,i") * t2_2["abab"]("b,a,j,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += eri["aaaa_vovo"]("b,l,e,j") * t2_1["abab"]("e,a,l,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += 2.00 * dp["aa_vv"]("b,e") * t2_2["abab"]("e,a,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= eri["abba_vovo"]("b,k,c,j") * t2_1["bbbb"]("c,a,i,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= dp["aa_oo"]("l,j") * t2["abab"]("b,a,l,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= 2.00 * dp["aa_oo"]("l,j") * t2_2["abab"]("b,a,l,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += eri["abba_vvvo"]("b,a,c,j") * t1_1["bb"]("c,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= dp["bb_oo"]("k,i") * t2["abab"]("b,a,j,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += 2.00 * dp["bb_vo"]("a,i") * t1_2["aa"]("b,j");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= eri["baab_vovo"]("a,l,e,i") * t2_1["aaaa"]("e,b,j,l");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += eri["baba_vovo"]("a,l,c,j") * t2_1["abab"]("b,c,l,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= w0 * t2_1["abab"]("b,a,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += dp["bb_vv"]("a,c") * t2["abab"]("b,c,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += eri["abab_vooo"]("b,k,j,i") * t1_1["bb"]("a,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= f["bb_vv"]("a,c") * t2_1["abab"]("b,c,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += scalars_["1"] * t2_1["abab"]("b,a,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= eri["baab_vooo"]("a,l,j,i") * t1_1["aa"]("b,l");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += dp["aa_vv"]("b,e") * t2["abab"]("e,a,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= f["aa_vv"]("b,e") * t2_1["abab"]("e,a,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += f["bb_oo"]("k,i") * t2_1["abab"]("b,a,j,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += 2.00 * scalars_["5"] * t2_2["abab"]("b,a,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += t1_1["bb"]("a,k") * tmps_["101_abab_vooo"]("b,k,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= 2.00 * t2["abab"]("b,a,l,i") * tmps_["42_aa_oo"]("l,j");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += t2_1["bbbb"]("f,a,i,n") * tmps_["3_bbaa_ovvo"]("n,f,b,j");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += 2.00 * t2["abab"]("b,c,j,i") * tmps_["60_bb_vv"]("a,c");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= t1_1["aa"]("b,l") * tmps_["99_abab_ovoo"]("l,a,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= tmps_["12_aa_vo"]("b,l") * tmps_["93_abab_ovoo"]("l,a,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= t2_1["abab"]("b,a,l,i") * tmps_["41_aa_oo"]("l,j");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= 2.00 * t2["abab"]("b,a,j,k") * tmps_["62_bb_oo"]("k,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= t1_1["aa"]("b,j") * tmps_["103_bb_vo"]("a,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += t2["abab"]("b,a,m,k") * tmps_["105_abba_oooo"]("m,k,i,j");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += t2_1["abab"]("e,a,j,i") * tmps_["10_aa_vv"]("b,e");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= t1_1["aa"]("b,j") * tmps_["104_bb_vo"]("a,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += t1_1["bb"]("a,n") * tmps_["98_bbaa_oovo"]("n,i,b,j");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += 2.00 * t2["abab"]("e,a,j,i") * tmps_["40_aa_vv"]("b,e");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += t2_1["abab"]("d,a,m,i") * tmps_["1_aaaa_ovvo"]("m,d,b,j");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += t1_1["aa"]("b,m") * tmps_["100_abba_oovo"]("m,i,a,j");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += 2.00 * t1_2["bb"]("a,i") * tmps_["75_aa_vo"]("b,j");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= t2_1["abab"]("b,a,j,k") * tmps_["61_bb_oo"]("k,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += t2_1["abab"]("b,c,j,i") * tmps_["25_bb_vv"]("a,c");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += t1_1["bb"]("a,i") * tmps_["76_aa_ov"]("j,b");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += t1_1["bb"]("a,i") * tmps_["74_aa_vo"]("b,j");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= tmps_["26_bb_vo"]("a,k") * tmps_["91_baab_ovoo"]("k,b,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += 2.00 * t1_2["aa"]("b,j") * tmps_["102_bb_vo"]("a,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= t2["abab"]("b,a,j,k") * tmps_["106_bb_oo"]("k,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= t2["abab"]("b,a,l,i") * tmps_["81_aa_oo"]("l,j");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= 0.50 * t2_1["abab"]("d,a,j,i") * tmps_["4_aa_vv"]("d,b");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= 0.50 * tmps_["6_aa_oo"]("m,j") * t2_1["abab"]("b,a,m,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= 0.50 * tmps_["21_bb_vv"]("f,a") * t2_1["abab"]("b,f,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= t2["bbbb"]("c,a,i,k") * tmps_["28_bbaa_ovvo"]("k,c,b,j");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= tmps_["51_baab_vovo"]("a,l,e,i") * t2["aaaa"]("e,b,j,l");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= tmps_["87_abba_ovvo"]("l,c,a,j") * t2["abab"]("b,c,l,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= tmps_["37_aa_oo"]("l,j") * t2["abab"]("b,a,l,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= tmps_["88_abab_oooo"]("m,k,j,i") * t2_1["abab"]("b,a,m,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= 2.00 * tmps_["94_abab_ovoo"]("l,a,j,i") * t1_1["aa"]("b,l");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += tmps_["50_bb_vv"]("c,a") * t2["abab"]("b,c,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += tmps_["20_bb_vv"]("f,a") * t2_1["abab"]("b,f,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += tmps_["31_abba_vovo"]("b,k,c,j") * t2["bbbb"]("c,a,i,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += tmps_["30_aaaa_vovo"]("b,l,e,j") * t2["abab"]("e,a,l,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= t2_1["bbbb"]("f,a,i,n") * tmps_["2_bbaa_ovvo"]("n,f,b,j");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= tmps_["18_bbbb_ovvo"]("n,f,a,i") * t2_1["abab"]("b,f,j,n");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= t2["abab"]("e,a,j,i") * tmps_["71_aa_vv"]("b,e");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= 2.00 * tmps_["93_abab_ovoo"]("l,a,j,i") * t1_2["aa"]("b,l");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= tmps_["86_abba_ovvo"]("m,f,a,j") * t2_1["abab"]("b,f,m,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= t2["abab"]("e,a,j,i") * tmps_["72_aa_vv"]("b,e");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= tmps_["97_bb_oo"]("k,i") * t2["abab"]("b,a,j,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += tmps_["59_bb_oo"]("k,i") * t2["abab"]("b,a,j,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= tmps_["48_bbbb_ovvo"]("k,c,a,i") * t2["abab"]("b,c,j,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= 2.00 * tmps_["58_bb_oo"]("k,i") * t2["abab"]("b,a,j,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += tmps_["89_abab_vovo"]("b,k,e,i") * t2["abab"]("e,a,j,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += tmps_["7_aa_oo"]("m,j") * t2_1["abab"]("b,a,m,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += 0.50 * tmps_["55_bb_oo"]("k,i") * t2["abab"]("b,a,j,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += tmps_["5_aa_vv"]("d,b") * t2_1["abab"]("d,a,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += tmps_["34_aa_oo"]("l,j") * t2["abab"]("b,a,l,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += tmps_["19_bbbb_ovvo"]("n,f,a,i") * t2_1["abab"]("b,f,j,n");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += tmps_["53_bbbb_vovo"]("a,k,c,i") * t2["abab"]("b,c,j,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += tmps_["32_aa_vv"]("e,b") * t2["abab"]("e,a,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= 2.00 * tmps_["38_aa_oo"]("l,j") * t2["abab"]("b,a,l,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += 0.50 * tmps_["33_aa_oo"]("l,j") * t2["abab"]("b,a,l,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += tmps_["17_aabb_ovvo"]("m,d,a,i") * t2_1["aaaa"]("d,b,j,m");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= 2.00 * t1_1["bb"]("a,k") * tmps_["92_baab_ovoo"]("k,b,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += 0.50 * tmps_["52_bb_vv"]("c,a") * t2["abab"]("b,c,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += tmps_["22_bb_oo"]("n,i") * t2_1["abab"]("b,a,j,n");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= 2.00 * tmps_["9_aa_oo"]("l,j") * t2_1["abab"]("b,a,l,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= tmps_["90_baba_vovo"]("a,l,c,j") * t2["abab"]("b,c,l,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += tmps_["56_bb_oo"]("k,i") * t2["abab"]("b,a,j,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += tmps_["54_bb_oo"]("k,i") * t2["abab"]("b,a,j,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= tmps_["80_aa_oo"]("l,j") * t2["abab"]("b,a,l,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += 0.50 * t2["abab"]("e,a,j,i") * tmps_["29_aa_vv"]("e,b");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += tmps_["39_aa_oo"]("l,j") * t2["abab"]("b,a,l,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= tmps_["96_bb_vv"]("a,c") * t2["abab"]("b,c,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") += tmps_["95_bb_vv"]("a,c") * t2["abab"]("b,c,j,i");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= 2.00 * tmps_["24_bb_oo"]("k,i") * t2_1["abab"]("b,a,j,k");
        tmps_["107_bbaa_vovo"]("a,i,b,j") -= 2.00 * t1_2["bb"]("a,k") * tmps_["91_baab_ovoo"]("k,b,j,i");
        tmps_["80_aa_oo"].~TArrayD();
        tmps_["71_aa_vv"].~TArrayD();
        tmps_["33_aa_oo"].~TArrayD();
        tmps_["30_aaaa_vovo"].~TArrayD();

        // rt2_1_abab  = -0.50 <l,k||d,c>_abab t2_abab(a,c,i,j) t2_1_abab(d,b,l,k)
        //            += -0.50 <k,l||d,c>_abab t2_abab(a,c,i,j) t2_1_abab(d,b,k,l)
        //            += +2.00 d-_aa(k,c) t1_1_aa(a,k) t2_1_abab(c,b,i,j)
        //            += +1.00 <l,k||d,c>_abab t2_bbbb(c,b,j,k) t2_1_aaaa(d,a,i,l)
        //            += -0.50 <l,k||c,d>_abab t2_abab(c,b,l,k) t2_1_abab(a,d,i,j)
        //            += -0.50 <k,l||c,d>_abab t2_abab(c,b,k,l) t2_1_abab(a,d,i,j)
        //            += -1.00 <a,k||d,c>_abab t2_bbbb(c,b,j,k) t1_1_aa(d,i)
        //            += +1.00 <k,a||c,d>_aaaa t2_abab(c,b,k,j) t1_1_aa(d,i)
        //            += +1.00 <k,l||d,c>_abab t2_abab(a,c,k,j) t2_1_abab(d,b,i,l)
        //            += -1.00 <k,l||i,c>_abab t2_abab(a,b,k,j) t1_1_bb(c,l)
        //            += +2.00 d-_aa(k,c) t2_abab(c,b,i,j) t1_2_aa(a,k)
        //            += -1.00 f_bb(k,c) t2_abab(a,b,i,k) t1_1_bb(c,j)
        //            += +2.00 d-_bb(k,c) t2_abab(a,b,i,k) t1_2_bb(c,j)
        //            += -1.00 <a,k||c,d>_abab t2_abab(c,b,i,k) t1_1_bb(d,j)
        //            += -1.00 <k,b||c,d>_abab t2_aaaa(c,a,i,k) t1_1_bb(d,j)
        //            += -0.50 <l,k||c,d>_abab t2_abab(c,d,i,k) t2_1_abab(a,b,l,j)
        //            += -0.50 <l,k||d,c>_abab t2_abab(d,c,i,k) t2_1_abab(a,b,l,j)
        //            += -0.50 <l,k||c,d>_bbbb t2_bbbb(c,b,l,k) t2_1_abab(a,d,i,j)
        //            += +0.50 <l,k||c,d>_bbbb t2_abab(a,b,i,k) t2_1_bbbb(c,d,j,l)
        //            += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,a,l,k) t2_1_abab(d,b,i,j)
        //            += -0.50 <l,k||d,c>_abab t2_abab(a,c,l,k) t2_1_abab(d,b,i,j)
        //            += -0.50 <k,l||d,c>_abab t2_abab(a,c,k,l) t2_1_abab(d,b,i,j)
        //            += -0.50 <k,l||c,d>_abab t2_abab(a,b,k,j) t2_1_abab(c,d,i,l)
        //            += -0.50 <k,l||d,c>_abab t2_abab(a,b,k,j) t2_1_abab(d,c,i,l)
        //            += +0.25 <l,k||c,d>_abab t2_abab(c,d,i,j) t2_1_abab(a,b,l,k)
        //            += +0.25 <k,l||c,d>_abab t2_abab(c,d,i,j) t2_1_abab(a,b,k,l)
        //            += +0.25 <l,k||d,c>_abab t2_abab(d,c,i,j) t2_1_abab(a,b,l,k)
        //            += +0.25 <k,l||d,c>_abab t2_abab(d,c,i,j) t2_1_abab(a,b,k,l)
        //            += +1.00 <l,k||c,d>_bbbb t2_bbbb(c,b,j,k) t2_1_abab(a,d,i,l)
        //            += +1.00 <k,b||c,d>_bbbb t2_abab(a,c,i,k) t1_1_bb(d,j)
        //            += -0.50 <l,k||c,d>_abab t2_abab(c,b,i,j) t2_1_abab(a,d,l,k)
        //            += -0.50 <k,l||c,d>_abab t2_abab(c,b,i,j) t2_1_abab(a,d,k,l)
        //            += +2.00 d-_aa(k,c) t2_abab(a,b,k,j) t1_2_aa(c,i)
        //            += +0.50 <l,k||c,d>_aaaa t2_abab(a,b,k,j) t2_1_aaaa(c,d,i,l)
        //            += +1.00 <k,l||c,d>_abab t2_aaaa(c,a,i,k) t2_1_bbbb(d,b,j,l)
        //            += +1.00 <k,l||c,d>_abab t2_abab(c,b,k,j) t2_1_abab(a,d,i,l)
        //            += +1.00 <l,k||c,d>_aaaa t2_abab(c,b,k,j) t2_1_aaaa(d,a,i,l)
        //            += -1.00 <k,a||c,d>_aaaa t2_abab(c,b,i,j) t1_1_aa(d,k)
        //            += +1.00 <l,k||c,d>_abab t2_abab(c,b,i,k) t2_1_abab(a,d,l,j)
        //            += +2.00 d-_bb(k,c) t2_1_abab(a,c,i,j) t1_1_bb(b,k)
        //            += +1.00 <a,k||c,d>_abab t2_abab(c,b,i,j) t1_1_bb(d,k)
        //            += -1.00 <l,k||c,j>_bbbb t2_abab(a,b,i,k) t1_1_bb(c,l)
        //            += +0.50 <l,k||c,d>_bbbb t2_abab(a,c,i,j) t2_1_bbbb(d,b,l,k)
        //            += +1.00 <l,k||d,c>_abab t2_abab(a,c,i,k) t2_1_abab(d,b,l,j)
        //            += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,d,i,k) t2_1_abab(a,b,l,j)
        //            += +1.00 d-_aa(k,c) t2_abab(a,b,k,j) t0_1 t1_1_aa(c,i)
        //            += +1.00 d-_bb(k,c) t2_abab(a,b,i,k) t0_1 t1_1_bb(c,j)
        //            += +2.00 d-_aa(k,i) t2_abab(a,b,k,j) t0_2
        //            += -1.00 f_bb(k,c) t2_abab(a,c,i,j) t1_1_bb(b,k)
        //            += -0.50 <a,k||c,d>_abab t2_abab(c,d,i,j) t1_1_bb(b,k)
        //            += -0.50 <a,k||d,c>_abab t2_abab(d,c,i,j) t1_1_bb(b,k)
        //            += -2.00 d-_bb(b,c) t2_abab(a,c,i,j) t0_2
        //            += +1.00 <l,k||c,d>_bbbb t2_abab(a,c,i,k) t2_1_bbbb(d,b,j,l)
        //            += -0.50 <k,b||c,d>_abab t2_abab(c,d,i,j) t1_1_aa(a,k)
        //            += -0.50 <k,b||d,c>_abab t2_abab(d,c,i,j) t1_1_aa(a,k)
        //            += -1.00 f_aa(k,c) t2_abab(c,b,i,j) t1_1_aa(a,k)
        //            += +1.00 d-_aa(k,c) t2_abab(c,b,i,j) t0_1 t1_1_aa(a,k)
        //            += +1.00 d-_aa(k,i) t0_1 t2_1_abab(a,b,k,j)
        //            += +2.00 d-_bb(k,j) t2_abab(a,b,i,k) t0_2
        //            += +1.00 d-_bb(k,c) t1_1_aa(a,i) t2_1_bbbb(c,b,j,k)
        //            += -1.00 d-_aa(k,c) t1_1_aa(a,i) t2_1_abab(c,b,k,j)
        //            += +0.50 <l,k||i,c>_abab t2_abab(a,b,l,k) t1_1_bb(c,j)
        //            += +0.50 <k,l||i,c>_abab t2_abab(a,b,k,l) t1_1_bb(c,j)
        //            += +0.50 <l,k||c,j>_abab t2_abab(a,b,l,k) t1_1_aa(c,i)
        //            += +0.50 <k,l||c,j>_abab t2_abab(a,b,k,l) t1_1_aa(c,i)
        //            += +0.25 <l,k||c,d>_abab t2_abab(a,b,l,k) t2_1_abab(c,d,i,j)
        //            += +0.25 <l,k||d,c>_abab t2_abab(a,b,l,k) t2_1_abab(d,c,i,j)
        //            += +0.25 <k,l||c,d>_abab t2_abab(a,b,k,l) t2_1_abab(c,d,i,j)
        //            += +0.25 <k,l||d,c>_abab t2_abab(a,b,k,l) t2_1_abab(d,c,i,j)
        //            += +1.00 d-_bb(k,j) t1_1_aa(a,i) t1_1_bb(b,k)
        //            += -1.00 d-_bb(b,c) t1_1_aa(a,i) t1_1_bb(c,j)
        //            += -1.00 d-_aa(a,c) t0_1 t2_1_abab(c,b,i,j)
        //            += -2.00 d-_aa(a,c) t2_abab(c,b,i,j) t0_2
        //            += +1.00 <l,k||c,j>_bbbb t2_abab(a,c,i,k) t1_1_bb(b,l)
        //            += +1.00 <k,l||i,c>_abab t2_abab(a,c,k,j) t1_1_bb(b,l)
        //            += +1.00 <k,l||c,j>_abab t2_aaaa(c,a,i,k) t1_1_bb(b,l)
        //            += -2.00 d-_bb(k,c) t2_abab(a,c,i,k) t1_2_bb(b,j)
        //            += +2.00 d-_aa(k,c) t2_aaaa(c,a,i,k) t1_2_bb(b,j)
        //            += +1.00 <l,k||c,d>_aaaa t2_aaaa(c,a,i,k) t2_1_abab(d,b,l,j)
        //            += +1.00 <l,k||i,c>_abab t2_bbbb(c,b,j,k) t1_1_aa(a,l)
        //            += +1.00 <l,k||c,j>_abab t2_abab(c,b,i,k) t1_1_aa(a,l)
        //            += +1.00 <l,k||c,i>_aaaa t2_abab(c,b,k,j) t1_1_aa(a,l)
        //            += +1.00 d-_bb(k,j) t0_1 t2_1_abab(a,b,i,k)
        //            += -2.00 d-_aa(a,i) t1_2_bb(b,j)
        //            += +1.00 <k,b||c,j>_bbbb t2_1_abab(a,c,i,k)
        //            += -1.00 f_aa(k,i) t2_1_abab(a,b,k,j)
        //            += -2.00 d-_bb(b,c) t2_2_abab(a,c,i,j)
        //            += +0.50 <l,k||i,j>_abab t2_1_abab(a,b,l,k)
        //            += +0.50 <k,l||i,j>_abab t2_1_abab(a,b,k,l)
        //            += -1.00 <a,k||c,j>_abab t2_1_abab(c,b,i,k)
        //            += +0.50 <a,b||c,d>_abab t2_1_abab(c,d,i,j)
        //            += +0.50 <a,b||d,c>_abab t2_1_abab(d,c,i,j)
        //            += -2.00 d-_aa(a,c) t2_2_abab(c,b,i,j)
        //            += +1.00 <a,b||c,j>_abab t1_1_aa(c,i)
        //            += +2.00 d-_bb(k,j) t2_2_abab(a,b,i,k)
        //            += +1.00 <k,a||c,i>_aaaa t2_1_abab(c,b,k,j)
        //            += -2.00 d-_bb(b,j) t1_2_aa(a,i)
        //            += -1.00 <a,k||i,c>_abab t2_1_bbbb(c,b,j,k)
        //            += +1.00 d+_aa(k,i) t2_abab(a,b,k,j)
        //            += +2.00 d-_aa(k,i) t2_2_abab(a,b,k,j)
        //            += +1.00 <a,b||i,c>_abab t1_1_bb(c,j)
        //            += +1.00 d+_bb(k,j) t2_abab(a,b,i,k)
        //            += +1.00 t2_1_abab(a,b,i,j) w0
        //            += -1.00 <k,b||c,j>_abab t2_1_aaaa(c,a,i,k)
        //            += -1.00 <k,b||i,c>_abab t2_1_abab(a,c,k,j)
        //            += -1.00 d-_bb(k,c) t2_1_abab(a,b,i,j) t1_1_bb(c,k)
        //            += -1.00 d-_aa(k,c) t2_1_abab(a,b,i,j) t1_1_aa(c,k)
        //            += -1.00 d+_bb(b,c) t2_abab(a,c,i,j)
        //            += -1.00 <a,k||i,j>_abab t1_1_bb(b,k)
        //            += +1.00 f_bb(b,c) t2_1_abab(a,c,i,j)
        //            += -1.00 f_bb(k,j) t2_1_abab(a,b,i,k)
        //            += -1.00 <k,b||i,j>_abab t1_1_aa(a,k)
        //            += -1.00 d+_aa(a,c) t2_abab(c,b,i,j)
        //            += +1.00 f_aa(a,c) t2_1_abab(c,b,i,j)
        //            += -2.00 d-_aa(k,k) t2_2_abab(a,b,i,j)
        //            += -2.00 d-_bb(k,k) t2_2_abab(a,b,i,j)
        //            += +1.00 d-_bb(k,c) t2_abab(a,c,i,j) t0_1 t1_1_bb(b,k)
        //            += -1.00 d-_bb(b,c) t0_1 t2_1_abab(a,c,i,j)
        //            += -1.00 d-_aa(a,c) t1_1_bb(b,j) t1_1_aa(c,i)
        //            += +1.00 d-_aa(k,i) t1_1_aa(a,k) t1_1_bb(b,j)
        //            += -1.00 d-_bb(k,c) t2_1_abab(a,c,i,k) t1_1_bb(b,j)
        //            += +1.00 d-_aa(k,c) t2_1_aaaa(c,a,i,k) t1_1_bb(b,j)
        //            += -2.00 d-_aa(k,c) t2_abab(c,b,k,j) t1_2_aa(a,i)
        //            += +2.00 d-_bb(k,c) t2_bbbb(c,b,j,k) t1_2_aa(a,i)
        //            += -0.50 <k,l||c,d>_abab t2_abab(c,d,k,j) t2_1_abab(a,b,i,l)
        //            += -0.50 <k,l||d,c>_abab t2_abab(d,c,k,j) t2_1_abab(a,b,i,l)
        //            += +2.00 d-_aa(k,c) t2_1_abab(a,b,k,j) t1_1_aa(c,i)
        //            += -1.00 <k,b||d,c>_abab t2_abab(a,c,k,j) t1_1_aa(d,i)
        //            += -1.00 <l,k||c,j>_abab t2_abab(a,b,i,k) t1_1_aa(c,l)
        //            += -0.50 <l,k||c,d>_bbbb t2_bbbb(c,d,j,k) t2_1_abab(a,b,i,l)
        //            += -0.50 <l,k||c,d>_abab t2_abab(a,b,i,k) t2_1_abab(c,d,l,j)
        //            += -0.50 <l,k||d,c>_abab t2_abab(a,b,i,k) t2_1_abab(d,c,l,j)
        //            += -1.00 <l,k||c,i>_aaaa t2_abab(a,b,k,j) t1_1_aa(c,l)
        //            += +0.50 <l,k||c,d>_aaaa t2_abab(c,b,i,j) t2_1_aaaa(d,a,l,k)
        //            += -1.00 f_aa(k,c) t2_abab(a,b,k,j) t1_1_aa(c,i)
        //            += -1.00 <k,b||c,d>_bbbb t2_abab(a,c,i,j) t1_1_bb(d,k)
        //            += +1.00 <k,b||d,c>_abab t2_abab(a,c,i,j) t1_1_aa(d,k)
        //            += +2.00 d-_bb(k,c) t2_1_abab(a,b,i,k) t1_1_bb(c,j)
        //            += +2.00 d-_bb(k,c) t2_abab(a,c,i,j) t1_2_bb(b,k)
        rt2_1_abab("a,b,i,j")  = -1.00 * tmps_["107_bbaa_vovo"]("b,j,a,i");
        tmps_["107_bbaa_vovo"].~TArrayD();

        // flops: o4v0  = o4v1
        //  mems: o4v0  = o4v0
        tmps_["108_bbbb_oooo"]("i,j,k,l")  = eri["bbbb_oovo"]("i,j,a,k") * t1_1["bb"]("a,l");

        // flops: o2v2  = o3v2 o3v2 o2v2 o3v2 o3v2 o2v2 o2v3 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o2v2 o3v2 o2v2 o3v2 o3v2 o2v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o4v2 o2v2 o3v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["109_bbbb_ovvo"]("i,a,b,j")  = -0.50 * tmps_["55_bb_oo"]("m,j") * t2["bbbb"]("a,b,i,m");
        tmps_["109_bbbb_ovvo"]("i,a,b,j") -= 0.50 * tmps_["23_bb_oo"]("k,i") * t2_1["bbbb"]("a,b,j,k");
        tmps_["109_bbbb_ovvo"]("i,a,b,j") += f["bb_oo"]("m,i") * t2_1["bbbb"]("a,b,j,m");
        tmps_["109_bbbb_ovvo"]("i,a,b,j") -= dp["bb_oo"]("m,i") * t2["bbbb"]("a,b,j,m");
        tmps_["109_bbbb_ovvo"]("i,a,b,j") -= eri["bbbb_vvvo"]("a,b,e,i") * t1_1["bb"]("e,j");
        tmps_["109_bbbb_ovvo"]("i,a,b,j") -= 2.00 * dp["bb_oo"]("m,i") * t2_2["bbbb"]("a,b,j,m");
        tmps_["109_bbbb_ovvo"]("i,a,b,j") -= t2_1["bbbb"]("a,b,j,m") * tmps_["61_bb_oo"]("m,i");
        tmps_["109_bbbb_ovvo"]("i,a,b,j") -= 2.00 * t2["bbbb"]("a,b,j,m") * tmps_["62_bb_oo"]("m,i");
        tmps_["109_bbbb_ovvo"]("i,a,b,j") += t2["bbbb"]("a,b,i,m") * tmps_["106_bb_oo"]("m,j");
        tmps_["109_bbbb_ovvo"]("i,a,b,j") -= 2.00 * tmps_["24_bb_oo"]("m,i") * t2_1["bbbb"]("a,b,j,m");
        tmps_["109_bbbb_ovvo"]("i,a,b,j") += tmps_["56_bb_oo"]("m,i") * t2["bbbb"]("a,b,j,m");
        tmps_["109_bbbb_ovvo"]("i,a,b,j") -= tmps_["59_bb_oo"]("m,j") * t2["bbbb"]("a,b,i,m");
        tmps_["109_bbbb_ovvo"]("i,a,b,j") -= tmps_["97_bb_oo"]("m,i") * t2["bbbb"]("a,b,j,m");
        tmps_["109_bbbb_ovvo"]("i,a,b,j") -= tmps_["54_bb_oo"]("m,j") * t2["bbbb"]("a,b,i,m");
        tmps_["109_bbbb_ovvo"]("i,a,b,j") += 2.00 * tmps_["58_bb_oo"]("m,j") * t2["bbbb"]("a,b,i,m");
        tmps_["109_bbbb_ovvo"]("i,a,b,j") += 0.50 * tmps_["108_bbbb_oooo"]("m,k,i,j") * t2["bbbb"]("a,b,k,m");
        tmps_["109_bbbb_ovvo"]("i,a,b,j") += tmps_["22_bb_oo"]("k,i") * t2_1["bbbb"]("a,b,j,k");
        tmps_["97_bb_oo"].~TArrayD();
        tmps_["55_bb_oo"].~TArrayD();

        // rt2_1_bbbb  = -1.00 P(i,j) <l,k||c,j>_abab t2_bbbb(a,b,i,k) t1_1_aa(c,l)
        //            += +2.00 P(i,j) d-_bb(k,c) t2_1_bbbb(a,b,i,k) t1_1_bb(c,j)
        //            += -1.00 P(i,j) d-_bb(k,c) t2_bbbb(a,b,j,k) t0_1 t1_1_bb(c,i)
        //            += +1.00 P(i,j) d+_bb(k,j) t2_bbbb(a,b,i,k)
        //            += -1.00 P(i,j) f_bb(k,j) t2_1_bbbb(a,b,i,k)
        //            += +1.00 P(i,j) <a,b||c,j>_bbbb t1_1_bb(c,i)
        //            += +2.00 P(i,j) d-_bb(k,j) t2_2_bbbb(a,b,i,k)
        //            += +1.00 P(i,j) d-_bb(k,j) t0_1 t2_1_bbbb(a,b,i,k)
        //            += +2.00 P(i,j) d-_bb(k,j) t2_bbbb(a,b,i,k) t0_2
        //            += +1.00 P(i,j) f_bb(k,c) t2_bbbb(a,b,j,k) t1_1_bb(c,i)
        //            += -0.50 P(i,j) <l,k||c,d>_bbbb t2_bbbb(c,d,j,k) t2_1_bbbb(a,b,i,l)
        //            += -0.50 P(i,j) <l,k||c,d>_bbbb t2_bbbb(a,b,j,k) t2_1_bbbb(c,d,i,l)
        //            += -1.00 P(i,j) <l,k||c,j>_bbbb t2_bbbb(a,b,i,k) t1_1_bb(c,l)
        //            += +0.50 P(i,j) <l,k||c,d>_abab t2_bbbb(a,b,j,k) t2_1_abab(c,d,l,i)
        //            += +0.50 P(i,j) <l,k||d,c>_abab t2_bbbb(a,b,j,k) t2_1_abab(d,c,l,i)
        //            += -2.00 P(i,j) d-_bb(k,c) t2_bbbb(a,b,j,k) t1_2_bb(c,i)
        //            += +0.50 P(i,j) <l,k||c,j>_bbbb t2_bbbb(a,b,l,k) t1_1_bb(c,i)
        //            += -0.50 P(i,j) <k,l||c,d>_abab t2_abab(c,d,k,j) t2_1_bbbb(a,b,i,l)
        //            += -0.50 P(i,j) <k,l||d,c>_abab t2_abab(d,c,k,j) t2_1_bbbb(a,b,i,l)
        rt2_1_bbbb("a,b,i,j")  = -1.00 * tmps_["109_bbbb_ovvo"]("j,a,b,i");
        rt2_1_bbbb("a,b,i,j") += tmps_["109_bbbb_ovvo"]("i,a,b,j");
        tmps_["109_bbbb_ovvo"].~TArrayD();

        // flops: o3v1  = o4v2 o4v2 o3v1
        //  mems: o3v1  = o3v1 o3v1 o3v1
        tmps_["110_bbbb_oovo"]("i,j,a,k")  = -1.00 * eri["abab_oovo"]("m,i,c,j") * t2["abab"]("c,a,m,k");
        tmps_["110_bbbb_oovo"]("i,j,a,k") += eri["bbbb_oovo"]("l,i,b,j") * t2["bbbb"]("b,a,k,l");

        // flops: o2v2  = o2v2 o2v2 o2v2 o3v2 o2v2 o3v3 o2v2 o2v2 o3v3 o2v2 o2v2 o2v2 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["111_bbbb_vovo"]("a,i,b,j")  = t1_1["bb"]("b,j") * tmps_["103_bb_vo"]("a,i");
        tmps_["111_bbbb_vovo"]("a,i,b,j") -= 2.00 * t1_2["bb"]("a,i") * tmps_["102_bb_vo"]("b,j");
        tmps_["111_bbbb_vovo"]("a,i,b,j") -= t1_1["bb"]("a,n") * tmps_["110_bbbb_oovo"]("n,j,b,i");
        tmps_["111_bbbb_vovo"]("a,i,b,j") += eri["bbbb_vovo"]("b,l,d,j") * t2_1["bbbb"]("d,a,i,l");
        tmps_["111_bbbb_vovo"]("a,i,b,j") -= 2.00 * dp["bb_vo"]("b,j") * t1_2["bb"]("a,i");
        tmps_["111_bbbb_vovo"]("a,i,b,j") -= eri["baab_vovo"]("b,m,e,j") * t2_1["abab"]("e,a,m,i");
        tmps_["111_bbbb_vovo"]("a,i,b,j") += t1_1["bb"]("a,i") * tmps_["104_bb_vo"]("b,j");
        tmps_["111_bbbb_vovo"]("a,i,b,j") += tmps_["19_bbbb_ovvo"]("n,f,b,j") * t2_1["bbbb"]("f,a,i,n");
        tmps_["111_bbbb_vovo"]("a,i,b,j") -= tmps_["18_bbbb_ovvo"]("n,f,b,j") * t2_1["bbbb"]("f,a,i,n");
        tmps_["111_bbbb_vovo"]("a,i,b,j") -= t2["bbbb"]("d,b,j,l") * tmps_["48_bbbb_ovvo"]("l,d,a,i");
        tmps_["111_bbbb_vovo"]("a,i,b,j") += tmps_["51_baab_vovo"]("b,m,e,i") * t2["abab"]("e,a,m,j");
        tmps_["111_bbbb_vovo"]("a,i,b,j") += tmps_["17_aabb_ovvo"]("k,c,b,j") * t2_1["abab"]("c,a,k,i");
        tmps_["111_bbbb_vovo"]("a,i,b,j") -= tmps_["53_bbbb_vovo"]("b,l,d,i") * t2["bbbb"]("d,a,j,l");
        tmps_["53_bbbb_vovo"].~TArrayD();

        // rt2_1_bbbb += +1.00 P(i,j) P(a,b) <l,k||c,d>_bbbb t2_bbbb(c,a,j,k) t2_1_bbbb(d,b,i,l)
        //            += -1.00 P(i,j) P(a,b) <l,k||c,j>_bbbb t2_bbbb(c,a,i,k) t1_1_bb(b,l)
        //            += -1.00 P(i,j) P(a,b) <k,l||c,j>_abab t2_abab(c,a,k,i) t1_1_bb(b,l)
        //            += +2.00 P(i,j) P(a,b) d-_aa(k,c) t2_abab(c,a,k,j) t1_2_bb(b,i)
        //            += -2.00 P(i,j) P(a,b) d-_bb(k,c) t2_bbbb(c,a,j,k) t1_2_bb(b,i)
        //            += -1.00 P(i,j) P(a,b) d-_bb(k,c) t1_1_bb(a,j) t2_1_bbbb(c,b,i,k)
        //            += +1.00 P(i,j) P(a,b) d-_aa(k,c) t1_1_bb(a,j) t2_1_abab(c,b,k,i)
        //            += +2.00 P(i,j) P(a,b) d-_bb(a,j) t1_2_bb(b,i)
        //            += +1.00 P(i,j) P(a,b) <k,a||c,j>_bbbb t2_1_bbbb(c,b,i,k)
        //            += -1.00 P(i,j) P(a,b) <k,a||c,j>_abab t2_1_abab(c,b,k,i)
        //            += -1.00 P(i,j) P(a,b) d-_bb(k,j) t1_1_bb(a,k) t1_1_bb(b,i)
        //            += +1.00 P(i,j) P(a,b) d-_bb(a,c) t1_1_bb(b,i) t1_1_bb(c,j)
        //            += +1.00 P(i,j) P(a,b) <k,l||c,d>_abab t2_abab(c,a,k,j) t2_1_bbbb(d,b,i,l)
        //            += +1.00 P(i,j) P(a,b) <l,k||d,c>_abab t2_bbbb(c,a,j,k) t2_1_abab(d,b,l,i)
        //            += +1.00 P(i,j) P(a,b) <k,a||c,d>_abab t2_abab(c,b,k,j) t1_1_bb(d,i)
        //            += +1.00 P(i,j) P(a,b) <l,k||c,d>_aaaa t2_abab(c,a,k,j) t2_1_abab(d,b,l,i)
        //            += -1.00 P(i,j) P(a,b) <k,a||c,d>_bbbb t2_bbbb(c,b,j,k) t1_1_bb(d,i)
        rt2_1_bbbb("a,b,i,j") -= tmps_["111_bbbb_vovo"]("b,i,a,j");
        rt2_1_bbbb("a,b,i,j") += tmps_["111_bbbb_vovo"]("b,j,a,i");
        rt2_1_bbbb("a,b,i,j") += tmps_["111_bbbb_vovo"]("a,i,b,j");
        rt2_1_bbbb("a,b,i,j") -= tmps_["111_bbbb_vovo"]("a,j,b,i");
        tmps_["111_bbbb_vovo"].~TArrayD();

        // flops: o3v1  = o3v2
        //  mems: o3v1  = o3v1
        tmps_["115_bbbb_ovoo"]("i,a,j,k")  = f["bb_ov"]("i,b") * t2["bbbb"]("b,a,j,k");

        // flops: o3v1  = o3v2
        //  mems: o3v1  = o3v1
        tmps_["114_bbbb_ovoo"]("i,a,j,k")  = dp["bb_ov"]("i,b") * t2_1["bbbb"]("b,a,j,k");

        // flops: o3v1  = o3v2
        //  mems: o3v1  = o3v1
        tmps_["113_bbbb_ovoo"]("i,a,j,k")  = dp["bb_ov"]("i,b") * t2["bbbb"]("b,a,j,k");

        // flops: o3v1  = o3v3
        //  mems: o3v1  = o3v1
        tmps_["112_bbbb_vooo"]("a,i,j,k")  = eri["bbbb_vovv"]("a,i,b,c") * t2["bbbb"]("b,c,j,k");

        // flops: o2v2  = o2v3 o3v2 o2v3 o2v2 o2v3 o2v2 o2v3 o2v2 o2v3 o2v3 o2v2 o2v2 o3v2 o2v2 o2v2 o3v2 o2v2 o2v3 o2v2 o2v3 o2v2 o2v3 o2v2 o2v3 o2v2 o3v2 o2v2 o2v3 o2v2 o3v2 o2v2 o3v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["116_bbbb_voov"]("a,i,j,b")  = -1.00 * t2["bbbb"]("c,b,i,j") * tmps_["95_bb_vv"]("a,c");
        tmps_["116_bbbb_voov"]("a,i,j,b") -= eri["bbbb_vooo"]("a,k,i,j") * t1_1["bb"]("b,k");
        tmps_["116_bbbb_voov"]("a,i,j,b") -= 2.00 * dp["bb_vv"]("a,c") * t2_2["bbbb"]("c,b,i,j");
        tmps_["116_bbbb_voov"]("a,i,j,b") -= dp["bb_vv"]("a,c") * t2["bbbb"]("c,b,i,j");
        tmps_["116_bbbb_voov"]("a,i,j,b") += f["bb_vv"]("a,c") * t2_1["bbbb"]("c,b,i,j");
        tmps_["116_bbbb_voov"]("a,i,j,b") -= t2_1["bbbb"]("c,b,i,j") * tmps_["25_bb_vv"]("a,c");
        tmps_["116_bbbb_voov"]("a,i,j,b") -= 2.00 * t2["bbbb"]("c,b,i,j") * tmps_["60_bb_vv"]("a,c");
        tmps_["116_bbbb_voov"]("a,i,j,b") -= tmps_["113_bbbb_ovoo"]("k,a,i,j") * tmps_["26_bb_vo"]("b,k");
        tmps_["116_bbbb_voov"]("a,i,j,b") -= 0.50 * t1_1["bb"]("b,k") * tmps_["112_bbbb_vooo"]("a,k,i,j");
        tmps_["116_bbbb_voov"]("a,i,j,b") -= t2_1["bbbb"]("d,b,i,j") * tmps_["20_bb_vv"]("d,a");
        tmps_["116_bbbb_voov"]("a,i,j,b") += tmps_["50_bb_vv"]("c,b") * t2["bbbb"]("c,a,i,j");
        tmps_["116_bbbb_voov"]("a,i,j,b") += tmps_["96_bb_vv"]("a,c") * t2["bbbb"]("c,b,i,j");
        tmps_["116_bbbb_voov"]("a,i,j,b") += 0.50 * tmps_["52_bb_vv"]("c,b") * t2["bbbb"]("c,a,i,j");
        tmps_["116_bbbb_voov"]("a,i,j,b") += 2.00 * tmps_["114_bbbb_ovoo"]("k,b,i,j") * t1_1["bb"]("a,k");
        tmps_["116_bbbb_voov"]("a,i,j,b") += 0.50 * t2_1["bbbb"]("d,b,i,j") * tmps_["21_bb_vv"]("d,a");
        tmps_["116_bbbb_voov"]("a,i,j,b") += tmps_["115_bbbb_ovoo"]("k,a,i,j") * t1_1["bb"]("b,k");
        tmps_["116_bbbb_voov"]("a,i,j,b") -= 2.00 * t1_2["bb"]("b,k") * tmps_["113_bbbb_ovoo"]("k,a,i,j");
        tmps_["96_bb_vv"].~TArrayD();

        // rt2_1_bbbb += +0.50 P(a,b) <l,k||d,c>_abab t2_bbbb(c,a,i,j) t2_1_abab(d,b,l,k)
        //            += +0.50 P(a,b) <k,l||d,c>_abab t2_bbbb(c,a,i,j) t2_1_abab(d,b,k,l)
        //            += -2.00 P(a,b) d-_bb(a,c) t2_bbbb(c,b,i,j) t0_2
        //            += -1.00 P(a,b) d-_bb(a,c) t0_1 t2_1_bbbb(c,b,i,j)
        //            += +1.00 P(a,b) f_bb(a,c) t2_1_bbbb(c,b,i,j)
        //            += -2.00 P(a,b) d-_bb(a,c) t2_2_bbbb(c,b,i,j)
        //            += +1.00 P(a,b) <k,a||i,j>_bbbb t1_1_bb(b,k)
        //            += -1.00 P(a,b) d+_bb(a,c) t2_bbbb(c,b,i,j)
        //            += -1.00 P(a,b) d-_bb(k,c) t2_bbbb(c,a,i,j) t0_1 t1_1_bb(b,k)
        //            += -1.00 P(a,b) <k,a||c,d>_bbbb t2_bbbb(c,b,i,j) t1_1_bb(d,k)
        //            += -0.50 P(a,b) <l,k||c,d>_bbbb t2_bbbb(c,a,i,j) t2_1_bbbb(d,b,l,k)
        //            += +2.00 P(a,b) d-_bb(k,c) t1_1_bb(a,k) t2_1_bbbb(c,b,i,j)
        //            += +1.00 P(a,b) <k,a||d,c>_abab t2_bbbb(c,b,i,j) t1_1_aa(d,k)
        //            += +0.50 P(a,b) <k,a||c,d>_bbbb t2_bbbb(c,d,i,j) t1_1_bb(b,k)
        //            += -0.50 P(a,b) <l,k||c,d>_abab t2_abab(c,a,l,k) t2_1_bbbb(d,b,i,j)
        //            += -0.50 P(a,b) <k,l||c,d>_abab t2_abab(c,a,k,l) t2_1_bbbb(d,b,i,j)
        //            += -0.50 P(a,b) <l,k||c,d>_bbbb t2_bbbb(c,a,l,k) t2_1_bbbb(d,b,i,j)
        //            += +1.00 P(a,b) f_bb(k,c) t2_bbbb(c,a,i,j) t1_1_bb(b,k)
        //            += -2.00 P(a,b) d-_bb(k,c) t2_bbbb(c,a,i,j) t1_2_bb(b,k)
        rt2_1_bbbb("a,b,i,j") += tmps_["116_bbbb_voov"]("a,i,j,b");
        rt2_1_bbbb("a,b,i,j") -= tmps_["116_bbbb_voov"]("b,i,j,a");
        tmps_["116_bbbb_voov"].~TArrayD();

        // flops: o4v0  = o4v2
        //  mems: o4v0  = o4v0
        tmps_["118_bbbb_oooo"]("i,j,k,l")  = eri["bbbb_oovv"]("i,j,a,b") * t2_1["bbbb"]("a,b,k,l");

        // flops: o4v0  = o4v2
        //  mems: o4v0  = o4v0
        tmps_["117_bbbb_oooo"]("i,j,k,l")  = eri["bbbb_oovv"]("i,j,a,b") * t2["bbbb"]("a,b,k,l");

        // flops: o2v2  = 0 o2v2 o2v2 o2v2 o2v4 o2v2 o4v2 o2v2 o4v2 o4v2 o2v2 o2v2
        //  mems: o2v2  = 0 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["119_bbbb_vvoo"]("a,b,i,j")  = 4.00 * (-1.00 * w0 + scalars_["1"]) * t2_1["bbbb"]("a,b,i,j");
        tmps_["119_bbbb_vvoo"]("a,b,i,j") += 8.00 * scalars_["5"] * t2_2["bbbb"]("a,b,i,j");
        tmps_["119_bbbb_vvoo"]("a,b,i,j") -= 2.00 * eri["bbbb_vvvv"]("a,b,c,d") * t2_1["bbbb"]("c,d,i,j");
        tmps_["119_bbbb_vvoo"]("a,b,i,j") += 2.00 * eri["bbbb_oooo"]("l,k,i,j") * t2_1["bbbb"]("a,b,k,l");
        tmps_["119_bbbb_vvoo"]("a,b,i,j") += tmps_["118_bbbb_oooo"]("l,k,i,j") * t2["bbbb"]("a,b,k,l");
        tmps_["119_bbbb_vvoo"]("a,b,i,j") += tmps_["117_bbbb_oooo"]("l,k,i,j") * t2_1["bbbb"]("a,b,k,l");

        // rt2_1_bbbb += +0.25 <l,k||c,d>_bbbb t2_bbbb(a,b,l,k) t2_1_bbbb(c,d,i,j)
        //            += +0.25 <l,k||c,d>_bbbb t2_bbbb(c,d,i,j) t2_1_bbbb(a,b,l,k)
        //            += -1.00 d-_bb(k,c) t2_1_bbbb(a,b,i,j) t1_1_bb(c,k)
        //            += -1.00 d-_aa(k,c) t2_1_bbbb(a,b,i,j) t1_1_aa(c,k)
        //            += +1.00 t2_1_bbbb(a,b,i,j) w0
        //            += -2.00 d-_aa(k,k) t2_2_bbbb(a,b,i,j)
        //            += -2.00 d-_bb(k,k) t2_2_bbbb(a,b,i,j)
        //            += +0.50 <a,b||c,d>_bbbb t2_1_bbbb(c,d,i,j)
        //            += +0.50 <l,k||i,j>_bbbb t2_1_bbbb(a,b,l,k)
        rt2_1_bbbb("a,b,i,j") -= 0.25 * tmps_["119_bbbb_vvoo"]("a,b,i,j");
        tmps_["119_bbbb_vvoo"].~TArrayD();

        // flops: o3v1  = o3v2
        //  mems: o3v1  = o3v1
        tmps_["124_abba_oovo"]("i,j,a,k")  = eri["abab_oovv"]("i,j,b,a") * t1_1["aa"]("b,k");

        // flops: o3v1  = o4v2
        //  mems: o3v1  = o3v1
        tmps_["132_aaaa_vooo"]("a,i,j,k")  = t2["abab"]("a,b,i,l") * tmps_["124_abba_oovo"]("j,l,b,k");

        // flops: o1v1  = o2v1
        //  mems: o1v1  = o1v1
        tmps_["131_aa_vo"]("a,i")  = t1_1["aa"]("a,j") * tmps_["9_aa_oo"]("j,i");

        // flops: o3v1  = o3v2
        //  mems: o3v1  = o3v1
        tmps_["130_aaaa_ovoo"]("i,a,j,k")  = t1_1["aa"]("b,k") * tmps_["1_aaaa_ovvo"]("i,b,a,j");

        // flops: o3v1  = o4v2 o4v2 o3v1
        //  mems: o3v1  = o3v1 o3v1 o3v1
        tmps_["129_aaaa_oovo"]("i,j,a,k")  = -1.00 * eri["abba_oovo"]("i,m,c,j") * t2_1["abab"]("a,c,k,m");
        tmps_["129_aaaa_oovo"]("i,j,a,k") += eri["aaaa_oovo"]("i,l,b,j") * t2_1["aaaa"]("b,a,k,l");

        // flops: o1v1  = o2v2 o2v2 o1v1
        //  mems: o1v1  = o1v1 o1v1 o1v1
        tmps_["128_aa_vo"]("a,i")  = -1.00 * dp["bb_ov"]("k,c") * t2_2["abab"]("a,c,i,k");
        tmps_["128_aa_vo"]("a,i") += dp["aa_ov"]("j,b") * t2_2["aaaa"]("b,a,i,j");

        // flops: o3v1  = o3v2
        //  mems: o3v1  = o3v1
        tmps_["127_aaaa_vooo"]("a,i,j,k")  = eri["aaaa_vovo"]("a,i,b,j") * t1_1["aa"]("b,k");
}
#endif
