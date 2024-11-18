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
void QED_CCSD_22::resid_22_4() {
    
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


        // flops: o1v1  = o2v1
        //  mems: o1v1  = o1v1
        tmps_["126_aa_ov"]("i,a")  = dp["aa_oo"]("j,i") * t1_2["aa"]("a,j");

        // flops: o1v1  = o1v2
        //  mems: o1v1  = o1v1
        tmps_["125_aa_vo"]("a,i")  = dp["aa_vv"]("a,b") * t1_2["aa"]("b,i");

        // flops: o2v2  = o2v3
        //  mems: o2v2  = o2v2
        tmps_["123_abba_vovo"]("a,i,b,j")  = eri["abab_vovv"]("a,i,c,b") * t1_2["aa"]("c,j");

        // flops: o2v2  = o2v3
        //  mems: o2v2  = o2v2
        tmps_["122_aaaa_vovo"]("a,i,b,j")  = eri["aaaa_vovv"]("a,i,c,b") * t1_1["aa"]("c,j");

        // flops: o2v2  = o2v3
        //  mems: o2v2  = o2v2
        tmps_["121_aaaa_vovo"]("a,i,b,j")  = eri["aaaa_vovv"]("a,i,b,c") * t1_2["aa"]("c,j");

        // flops: o2v2  = o3v3
        //  mems: o2v2  = o2v2
        tmps_["120_bbaa_ovvo"]("i,a,b,j")  = eri["abab_oovv"]("k,i,c,a") * t2_2["aaaa"]("c,b,j,k");

        // flops: o2v2  = o3v2 o2v2 o2v2 o3v3 o2v2 o2v2 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o2v2 o3v3 o3v3 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o3v2 o2v2 o3v2 o2v2 o2v2 o3v2 o2v2 o2v2 o2v2 o3v3 o2v2 o3v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["133_aaaa_vvoo"]("a,b,i,j")  = -1.00 * tmps_["130_aaaa_ovoo"]("k,b,i,j") * t1_1["aa"]("a,k");
        tmps_["133_aaaa_vvoo"]("a,b,i,j") -= tmps_["125_aa_vo"]("b,j") * t1_1["aa"]("a,i");
        tmps_["133_aaaa_vvoo"]("a,b,i,j") += tmps_["123_abba_vovo"]("b,l,c,j") * t2["abab"]("a,c,i,l");
        tmps_["133_aaaa_vvoo"]("a,b,i,j") += tmps_["126_aa_ov"]("i,a") * t1_1["aa"]("b,j");
        tmps_["133_aaaa_vvoo"]("a,b,i,j") -= tmps_["3_bbaa_ovvo"]("n,f,b,i") * t2_2["abab"]("a,f,j,n");
        tmps_["133_aaaa_vvoo"]("a,b,i,j") += tmps_["120_bbaa_ovvo"]("l,c,a,j") * t2["abab"]("b,c,i,l");
        tmps_["133_aaaa_vvoo"]("a,b,i,j") += tmps_["122_aaaa_vovo"]("b,m,e,i") * t2_1["aaaa"]("e,a,j,m");
        tmps_["133_aaaa_vvoo"]("a,b,i,j") -= tmps_["31_abba_vovo"]("b,l,f,i") * t2_1["abab"]("a,f,j,l");
        tmps_["133_aaaa_vvoo"]("a,b,i,j") -= tmps_["1_aaaa_ovvo"]("k,e,b,i") * t2_2["aaaa"]("e,a,j,k");
        tmps_["133_aaaa_vvoo"]("a,b,i,j") += tmps_["2_bbaa_ovvo"]("n,f,b,i") * t2_2["abab"]("a,f,j,n");
        tmps_["133_aaaa_vvoo"]("a,b,i,j") += 2.00 * -0.50 * (tmps_["128_aa_vo"]("a,j") * t1_1["aa"]("b,i") + -1.00 * eri["abba_vovo"]("b,l,c,i") * t2_2["abab"]("a,c,j,l") + eri["aaaa_vovo"]("b,m,d,i") * t2_2["aaaa"]("d,a,j,m"));
        tmps_["133_aaaa_vvoo"]("a,b,i,j") += 2.00 * t1_2["aa"]("a,j") * tmps_["76_aa_ov"]("i,b");
        tmps_["133_aaaa_vvoo"]("a,b,i,j") += 2.00 * t1_2["aa"]("a,j") * tmps_["74_aa_vo"]("b,i");
        tmps_["133_aaaa_vvoo"]("a,b,i,j") += t1_2["aa"]("a,k") * tmps_["77_aaaa_oovo"]("k,i,b,j");
        tmps_["133_aaaa_vvoo"]("a,b,i,j") += t1_1["aa"]("b,m") * tmps_["129_aaaa_oovo"]("m,i,a,j");
        tmps_["133_aaaa_vvoo"]("a,b,i,j") -= tmps_["127_aaaa_vooo"]("b,m,i,j") * t1_1["aa"]("a,m");
        tmps_["133_aaaa_vvoo"]("a,b,i,j") -= tmps_["131_aa_vo"]("b,i") * t1_1["aa"]("a,j");
        tmps_["133_aaaa_vvoo"]("a,b,i,j") += tmps_["121_aaaa_vovo"]("b,m,d,j") * t2["aaaa"]("d,a,i,m");
        tmps_["133_aaaa_vvoo"]("a,b,i,j") += tmps_["132_aaaa_vooo"]("b,i,k,j") * t1_1["aa"]("a,k");
        tmps_["132_aaaa_vooo"].~TArrayD();
        tmps_["130_aaaa_ovoo"].~TArrayD();
        tmps_["129_aaaa_oovo"].~TArrayD();
        tmps_["127_aaaa_vooo"].~TArrayD();
        tmps_["77_aaaa_oovo"].~TArrayD();

        // rt2_2_aaaa  = +2.00 P(i,j) P(a,b) <l,k||d,c>_abab t2_abab(a,c,j,k) t1_1_aa(b,l) t1_1_aa(d,i)
        //            += +2.00 P(i,j) P(a,b) <a,k||d,c>_abab t2_abab(b,c,j,k) t1_2_aa(d,i)
        //            += +2.00 P(i,j) P(a,b) <l,k||c,d>_aaaa t2_aaaa(c,a,j,k) t1_1_aa(b,l) t1_1_aa(d,i)
        //            += -2.00 P(i,j) P(a,b) d-_aa(a,c) t1_1_aa(b,j) t1_2_aa(c,i)
        //            += +2.00 P(i,j) P(a,b) d-_aa(k,j) t1_1_aa(a,i) t1_2_aa(b,k)
        //            += +2.00 P(i,j) P(a,b) <l,k||c,d>_bbbb t2_abab(a,c,j,k) t2_2_abab(b,d,i,l)
        //            += +2.00 P(i,j) P(a,b) <l,k||d,c>_abab t2_abab(a,c,j,k) t2_2_aaaa(d,b,i,l)
        //            += -2.00 P(i,j) P(a,b) <k,a||c,d>_aaaa t2_1_aaaa(d,b,i,k) t1_1_aa(c,j)
        //            += -2.00 P(i,j) P(a,b) <a,k||c,d>_abab t2_1_abab(b,d,i,k) t1_1_aa(c,j)
        //            += +2.00 P(i,j) P(a,b) <l,k||c,d>_aaaa t2_aaaa(c,a,j,k) t2_2_aaaa(d,b,i,l)
        //            += +2.00 P(i,j) P(a,b) <k,l||c,d>_abab t2_aaaa(c,a,j,k) t2_2_abab(b,d,i,l)
        //            += +4.00 P(i,j) P(a,b) d-_aa(a,c) t1_1_aa(c,j) t1_2_aa(b,i)
        //            += -4.00 P(i,j) P(a,b) d-_aa(k,j) t1_1_aa(a,k) t1_2_aa(b,i)
        //            += -2.00 P(i,j) P(a,b) d-_aa(k,c) t1_1_aa(a,j) t2_2_aaaa(c,b,i,k)
        //            += +2.00 P(i,j) P(a,b) d-_bb(k,c) t1_1_aa(a,j) t2_2_abab(b,c,i,k)
        //            += +2.00 P(i,j) P(a,b) <k,a||c,j>_aaaa t2_2_aaaa(c,b,i,k)
        //            += -2.00 P(i,j) P(a,b) <a,k||j,c>_abab t2_2_abab(b,c,i,k)
        //            += +4.00 P(i,j) P(a,b) d-_bb(k,c) t2_1_abab(a,c,j,k) t1_2_aa(b,i)
        //            += -4.00 P(i,j) P(a,b) d-_aa(k,c) t2_1_aaaa(c,a,j,k) t1_2_aa(b,i)
        //            += -2.00 P(i,j) P(a,b) <l,k||c,j>_aaaa t2_aaaa(c,a,i,k) t1_2_aa(b,l)
        //            += -2.00 P(i,j) P(a,b) <l,k||j,c>_abab t2_abab(a,c,i,k) t1_2_aa(b,l)
        //            += -2.00 P(i,j) P(a,b) <l,k||c,j>_aaaa t1_1_aa(a,k) t2_1_aaaa(c,b,i,l)
        //            += +2.00 P(i,j) P(a,b) <k,l||j,c>_abab t1_1_aa(a,k) t2_1_abab(b,c,i,l)
        //            += -2.00 P(i,j) P(a,b) <k,a||c,d>_aaaa t2_aaaa(c,b,j,k) t1_2_aa(d,i)
        //            += +2.00 P(i,j) P(a,b) <k,a||c,j>_aaaa t1_1_aa(b,k) t1_1_aa(c,i)
        //            += -2.00 P(i,j) P(a,b) d-_aa(k,c) t1_1_aa(a,k) t1_1_aa(b,i) t1_1_aa(c,j)
        rt2_2_aaaa("a,b,i,j")  = 2.00 * tmps_["133_aaaa_vvoo"]("b,a,j,i");
        rt2_2_aaaa("a,b,i,j") -= 2.00 * tmps_["133_aaaa_vvoo"]("b,a,i,j");
        rt2_2_aaaa("a,b,i,j") -= 2.00 * tmps_["133_aaaa_vvoo"]("a,b,j,i");
        rt2_2_aaaa("a,b,i,j") += 2.00 * tmps_["133_aaaa_vvoo"]("a,b,i,j");
        tmps_["133_aaaa_vvoo"].~TArrayD();

        // flops: o3v1  = o3v2
        //  mems: o3v1  = o3v1
        tmps_["135_aaaa_oovo"]("i,j,a,k")  = eri["aaaa_oovv"]("i,j,a,b") * t1_1["aa"]("b,k");

        // flops: o2v0  = o3v1 o2v0 o2v0 o3v1 o2v0 o2v0 o2v0
        //  mems: o2v0  = o2v0 o2v0 o2v0 o2v0 o2v0 o2v0 o2v0
        tmps_["141_aa_oo"]("i,j")  = -0.50 * t1_1["bb"]("b,l") * tmps_["124_abba_oovo"]("i,l,b,j");
        tmps_["141_aa_oo"]("i,j") += t0_2 * tmps_["9_aa_oo"]("i,j");
        tmps_["141_aa_oo"]("i,j") += 0.50 * t1_1["aa"]("a,k") * tmps_["135_aaaa_oovo"]("i,k,a,j");
        tmps_["141_aa_oo"]("i,j") += 0.50 * t0_1 * tmps_["38_aa_oo"]("i,j");

        // flops: o2v0  = o3v1 o3v1 o2v0
        //  mems: o2v0  = o2v0 o2v0 o2v0
        tmps_["140_aa_oo"]("i,j")  = eri["aaaa_oovo"]("i,k,a,j") * t1_2["aa"]("a,k");
        tmps_["140_aa_oo"]("i,j") += eri["abba_oovo"]("i,l,b,j") * t1_2["bb"]("b,l");

        // flops: o2v0  = o3v2 o3v2 o2v1 o2v0 o2v0
        //  mems: o2v0  = o2v0 o2v0 o2v0 o2v0 o2v0
        tmps_["139_aa_oo"]("i,j")  = 0.50 * eri["aaaa_oovv"]("i,k,a,b") * t2_2["aaaa"]("a,b,j,k");
        tmps_["139_aa_oo"]("i,j") += eri["abab_oovv"]("i,l,a,c") * t2_2["abab"]("a,c,j,l");
        tmps_["139_aa_oo"]("i,j") += f["aa_ov"]("i,a") * t1_2["aa"]("a,j");

        // flops: o2v2  = o3v3
        //  mems: o2v2  = o2v2
        tmps_["138_aaaa_ovvo"]("i,a,b,j")  = eri["aaaa_oovv"]("k,i,c,a") * t2_1["aaaa"]("c,b,j,k");

        // flops: o2v2  = o3v3
        //  mems: o2v2  = o2v2
        tmps_["137_bbaa_ovvo"]("i,a,b,j")  = eri["bbbb_oovv"]("k,i,c,a") * t2_1["abab"]("b,c,j,k");

        // flops: o4v0  = o4v1
        //  mems: o4v0  = o4v0
        tmps_["136_aaaa_oooo"]("i,j,k,l")  = eri["aaaa_oovo"]("i,j,a,k") * t1_2["aa"]("a,l");

        // flops: o2v0  = o3v2
        //  mems: o2v0  = o2v0
        tmps_["134_aa_oo"]("i,j")  = eri["aaaa_oovv"]("k,i,a,b") * t2_1["aaaa"]("a,b,j,k");

        // flops: o2v2  = o4v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o4v2 o2v2 o3v3 o3v2 o2v2 o3v2 o2v2 o3v3 o2v2 o3v2 o2v2 o2v3 o3v2 o2v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["142_aaaa_vovo"]("a,i,b,j")  = -1.00 * tmps_["136_aaaa_oooo"]("l,k,i,j") * t2["aaaa"]("a,b,k,l");
        tmps_["142_aaaa_vovo"]("a,i,b,j") += tmps_["134_aa_oo"]("k,i") * t2_1["aaaa"]("a,b,j,k");
        tmps_["142_aaaa_vovo"]("a,i,b,j") -= 6.00 * tmps_["38_aa_oo"]("l,j") * t2_1["aaaa"]("a,b,i,l");
        tmps_["142_aaaa_vovo"]("a,i,b,j") += 2.00 * tmps_["37_aa_oo"]("k,i") * t2_1["aaaa"]("a,b,j,k");
        tmps_["142_aaaa_vovo"]("a,i,b,j") -= 2.00 * tmps_["7_aa_oo"]("k,i") * t2_2["aaaa"]("a,b,j,k");
        tmps_["142_aaaa_vovo"]("a,i,b,j") -= 2.00 * tmps_["39_aa_oo"]("l,i") * t2_1["aaaa"]("a,b,j,l");
        tmps_["142_aaaa_vovo"]("a,i,b,j") -= tmps_["79_aaaa_oooo"]("l,k,i,j") * t2_1["aaaa"]("a,b,k,l");
        tmps_["142_aaaa_vovo"]("a,i,b,j") -= 2.00 * t2_1["abab"]("b,e,j,n") * tmps_["137_bbaa_ovvo"]("n,e,a,i");
        tmps_["142_aaaa_vovo"]("a,i,b,j") += 2.00 * t2_2["aaaa"]("a,b,j,l") * tmps_["41_aa_oo"]("l,i");
        tmps_["142_aaaa_vovo"]("a,i,b,j") += 2.00 * t2["aaaa"]("a,b,j,l") * tmps_["140_aa_oo"]("l,i");
        tmps_["142_aaaa_vovo"]("a,i,b,j") -= 2.00 * t2_1["aaaa"]("d,b,j,k") * tmps_["138_aaaa_ovvo"]("k,d,a,i");
        tmps_["142_aaaa_vovo"]("a,i,b,j") += 4.00 * t2_1["aaaa"]("a,b,j,l") * tmps_["42_aa_oo"]("l,i");
        tmps_["142_aaaa_vovo"]("a,i,b,j") += 2.00 * eri["aaaa_vvvo"]("a,b,c,i") * t1_2["aa"]("c,j");
        tmps_["142_aaaa_vovo"]("a,i,b,j") -= 2.00 * f["aa_oo"]("l,i") * t2_2["aaaa"]("a,b,j,l");
        tmps_["142_aaaa_vovo"]("a,i,b,j") += 2.00 * t2["aaaa"]("a,b,i,l") * tmps_["139_aa_oo"]("l,j");
        tmps_["142_aaaa_vovo"]("a,i,b,j") -= 4.00 * t2["aaaa"]("a,b,i,l") * tmps_["141_aa_oo"]("l,j");
        tmps_["142_aaaa_vovo"]("a,i,b,j") += 2.00 * t2_1["aaaa"]("a,b,j,l") * tmps_["81_aa_oo"]("l,i");
        tmps_["142_aaaa_vovo"]("a,i,b,j") -= 2.00 * tmps_["34_aa_oo"]("k,i") * t2_1["aaaa"]("a,b,j,k");
        tmps_["142_aaaa_vovo"]("a,i,b,j") -= 2.00 * tmps_["36_aa_oo"]("k,i") * t2_1["aaaa"]("a,b,j,k");
        tmps_["142_aaaa_vovo"]("a,i,b,j") += 6.00 * tmps_["9_aa_oo"]("l,i") * t2_2["aaaa"]("a,b,j,l");
        tmps_["142_aaaa_vovo"]("a,i,b,j") += tmps_["6_aa_oo"]("k,i") * t2_2["aaaa"]("a,b,j,k");
        tmps_["138_aaaa_ovvo"].~TArrayD();
        tmps_["137_bbaa_ovvo"].~TArrayD();
        tmps_["136_aaaa_oooo"].~TArrayD();
        tmps_["79_aaaa_oooo"].~TArrayD();

        // rt2_2_aaaa += -1.00 P(i,j) <l,k||c,d>_aaaa t2_1_aaaa(a,b,i,l) t2_1_aaaa(c,d,j,k)
        //            += +1.00 P(i,j) <l,k||c,j>_aaaa t2_aaaa(a,b,l,k) t1_2_aa(c,i)
        //            += -6.00 P(i,j) d-_aa(k,c) t2_1_aaaa(a,b,j,k) t1_2_aa(c,i)
        //            += -2.00 P(i,j) <l,k||j,c>_abab t2_1_aaaa(a,b,i,l) t1_1_bb(c,k)
        //            += -1.00 P(i,j) <l,k||c,d>_abab t2_abab(c,d,j,k) t2_2_aaaa(a,b,i,l)
        //            += -1.00 P(i,j) <l,k||d,c>_abab t2_abab(d,c,j,k) t2_2_aaaa(a,b,i,l)
        //            += -2.00 P(i,j) f_aa(k,c) t2_1_aaaa(a,b,i,k) t1_1_aa(c,j)
        //            += -4.00 P(i,j) d-_aa(k,c) t2_aaaa(a,b,j,k) t1_1_aa(c,i) t0_2
        //            += +2.00 P(i,j) <k,l||d,c>_abab t2_aaaa(a,b,j,k) t1_1_bb(c,l) t1_1_aa(d,i)
        //            += +2.00 P(i,j) <l,k||c,d>_aaaa t2_aaaa(a,b,j,k) t1_1_aa(c,l) t1_1_aa(d,i)
        //            += -2.00 P(i,j) d-_aa(k,c) t2_aaaa(a,b,j,k) t0_1 t1_2_aa(c,i)
        //            += -2.00 P(i,j) <k,l||j,c>_abab t2_aaaa(a,b,i,k) t1_2_bb(c,l)
        //            += -2.00 P(i,j) <l,k||c,j>_aaaa t2_aaaa(a,b,i,k) t1_2_aa(c,l)
        //            += +2.00 P(i,j) d-_aa(k,j) t0_1 t2_2_aaaa(a,b,i,k)
        //            += +2.00 P(i,j) <l,k||c,d>_bbbb t2_1_abab(a,c,j,k) t2_1_abab(b,d,i,l)
        //            += +2.00 P(i,j) <l,k||c,d>_aaaa t2_1_aaaa(c,a,j,k) t2_1_aaaa(d,b,i,l)
        //            += +4.00 P(i,j) d-_aa(k,j) t2_1_aaaa(a,b,i,k) t0_2
        //            += +2.00 P(i,j) f_aa(k,c) t2_aaaa(a,b,j,k) t1_2_aa(c,i)
        //            += +1.00 P(i,j) <k,l||c,d>_abab t2_aaaa(a,b,j,k) t2_2_abab(c,d,i,l)
        //            += +1.00 P(i,j) <k,l||d,c>_abab t2_aaaa(a,b,j,k) t2_2_abab(d,c,i,l)
        //            += -1.00 P(i,j) <l,k||c,d>_aaaa t2_aaaa(a,b,j,k) t2_2_aaaa(c,d,i,l)
        //            += -2.00 P(i,j) f_aa(k,j) t2_2_aaaa(a,b,i,k)
        //            += +2.00 P(i,j) <a,b||c,j>_aaaa t1_2_aa(c,i)
        //            += +2.00 P(i,j) d-_aa(k,c) t0_1 t2_1_aaaa(a,b,i,k) t1_1_aa(c,j)
        //            += +1.00 P(i,j) <l,k||c,j>_aaaa t2_1_aaaa(a,b,l,k) t1_1_aa(c,i)
        //            += -1.00 P(i,j) <l,k||c,d>_abab t2_1_aaaa(a,b,i,l) t2_1_abab(c,d,j,k)
        //            += -1.00 P(i,j) <l,k||d,c>_abab t2_1_aaaa(a,b,i,l) t2_1_abab(d,c,j,k)
        //            += +2.00 P(i,j) <l,k||c,j>_aaaa t2_1_aaaa(a,b,i,l) t1_1_aa(c,k)
        //            += +6.00 P(i,j) d-_aa(k,c) t1_1_aa(c,j) t2_2_aaaa(a,b,i,k)
        //            += -1.00 P(i,j) <l,k||c,d>_aaaa t2_aaaa(c,d,j,k) t2_2_aaaa(a,b,i,l)
        rt2_2_aaaa("a,b,i,j") += tmps_["142_aaaa_vovo"]("a,j,b,i");
        rt2_2_aaaa("a,b,i,j") -= tmps_["142_aaaa_vovo"]("a,i,b,j");
        tmps_["142_aaaa_vovo"].~TArrayD();

        // flops: o2v2  = o3v2 o2v3 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2
        tmps_["143_aaaa_vvoo"]("a,b,i,j")  = tmps_["68_aaaa_ovoo"]("k,b,i,j") * t1_1["aa"]("a,k");
        tmps_["143_aaaa_vvoo"]("a,b,i,j") += dp["aa_vv"]("b,c") * t2_1["aaaa"]("c,a,i,j");

        // rt2_aaaa  = -1.00 P(a,b) d-_aa(a,c) t2_1_aaaa(c,b,i,j)
        //          += -1.00 P(a,b) d-_aa(k,c) t2_aaaa(c,a,i,j) t1_1_aa(b,k)
        rt2_aaaa("a,b,i,j")  = -1.00 * tmps_["143_aaaa_vvoo"]("b,a,i,j");
        rt2_aaaa("a,b,i,j") += tmps_["143_aaaa_vvoo"]("a,b,i,j");

        // rt2_2_aaaa += -2.00 P(a,b) d+_aa(a,c) t2_1_aaaa(c,b,i,j)
        //            += -2.00 P(a,b) d+_aa(k,c) t2_aaaa(c,a,i,j) t1_1_aa(b,k)
        rt2_2_aaaa("a,b,i,j") -= 2.00 * tmps_["143_aaaa_vvoo"]("b,a,i,j");
        rt2_2_aaaa("a,b,i,j") += 2.00 * tmps_["143_aaaa_vvoo"]("a,b,i,j");
        tmps_["143_aaaa_vvoo"].~TArrayD();

        // flops: o4v0  = o4v1
        //  mems: o4v0  = o4v0
        tmps_["148_aaaa_oooo"]("i,j,k,l")  = t1_1["aa"]("a,i") * tmps_["135_aaaa_oovo"]("j,k,a,l");
        tmps_["135_aaaa_oovo"].~TArrayD();

        // flops: o3v1  = o4v1
        //  mems: o3v1  = o3v1
        tmps_["147_aaaa_vooo"]("a,i,j,k")  = tmps_["83_aaaa_oooo"]("i,l,j,k") * t1_1["aa"]("a,l");

        // flops: o3v1  = o4v1
        //  mems: o3v1  = o3v1
        tmps_["146_aaaa_ooov"]("i,j,k,a")  = eri["aaaa_oooo"]("l,i,j,k") * t1_1["aa"]("a,l");

        // flops: o4v0  = o4v2
        //  mems: o4v0  = o4v0
        tmps_["145_aaaa_oooo"]("i,j,k,l")  = eri["aaaa_oovv"]("i,j,a,b") * t2_2["aaaa"]("a,b,k,l");

        // flops: o0v2  = o2v3
        //  mems: o0v2  = o0v2
        tmps_["144_aa_vv"]("a,b")  = eri["aaaa_oovv"]("i,j,c,a") * t2_1["aaaa"]("c,b,j,i");

        // flops: o2v2  = o4v2 o2v2 o2v2 o1v4 o2v3 o2v2 o2v4 o2v2 o2v2 o2v2 o2v2 o2v2 o3v2 o2v2 o3v2 o2v2 o2v3 o2v2 o4v2 o2v2 o4v2 o2v2 o2v3 o2v2 o4v2 o2v2 o4v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o1v3 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["149_aaaa_oovv"]("i,j,a,b")  = -1.00 * eri["aaaa_oooo"]("l,k,i,j") * t2_2["aaaa"]("a,b,k,l");
        tmps_["149_aaaa_oovv"]("i,j,a,b") += 4.00 * w0 * t2_2["aaaa"]("a,b,i,j");
        tmps_["149_aaaa_oovv"]("i,j,a,b") -= 2.00 * eri["aaaa_vvvv"]("a,b,c,d") * t1_1["aa"]("d,i") * t1_1["aa"]("c,j");
        tmps_["149_aaaa_oovv"]("i,j,a,b") += eri["aaaa_vvvv"]("a,b,c,d") * t2_2["aaaa"]("c,d,i,j");
        tmps_["149_aaaa_oovv"]("i,j,a,b") -= 4.00 * scalars_["1"] * t2_2["aaaa"]("a,b,i,j");
        tmps_["149_aaaa_oovv"]("i,j,a,b") -= 2.00 * scalars_["6"] * t2_1["aaaa"]("a,b,i,j");
        tmps_["149_aaaa_oovv"]("i,j,a,b") += 2.00 * t1_1["aa"]("b,k") * tmps_["146_aaaa_ooov"]("k,i,j,a");
        tmps_["149_aaaa_oovv"]("i,j,a,b") += t1_1["aa"]("a,l") * tmps_["147_aaaa_vooo"]("b,l,i,j");
        tmps_["149_aaaa_oovv"]("i,j,a,b") += tmps_["29_aa_vv"]("c,b") * t2_1["aaaa"]("c,a,i,j");
        tmps_["149_aaaa_oovv"]("i,j,a,b") -= 0.50 * tmps_["84_aaaa_oooo"]("l,k,i,j") * t2_1["aaaa"]("a,b,k,l");
        tmps_["149_aaaa_oovv"]("i,j,a,b") -= 0.50 * tmps_["145_aaaa_oooo"]("l,k,i,j") * t2["aaaa"]("a,b,k,l");
        tmps_["149_aaaa_oovv"]("i,j,a,b") += tmps_["144_aa_vv"]("d,a") * t2_1["aaaa"]("d,b,i,j");
        tmps_["149_aaaa_oovv"]("i,j,a,b") -= 0.50 * tmps_["83_aaaa_oooo"]("l,k,i,j") * t2_2["aaaa"]("a,b,k,l");
        tmps_["149_aaaa_oovv"]("i,j,a,b") += tmps_["148_aaaa_oooo"]("j,l,k,i") * t2["aaaa"]("a,b,k,l");
        tmps_["148_aaaa_oooo"].~TArrayD();
        tmps_["147_aaaa_vooo"].~TArrayD();
        tmps_["146_aaaa_ooov"].~TArrayD();
        tmps_["145_aaaa_oooo"].~TArrayD();
        tmps_["84_aaaa_oooo"].~TArrayD();
        tmps_["29_aa_vv"].~TArrayD();

        // rt2_2_aaaa += -1.00 <l,k||c,d>_aaaa t2_aaaa(a,b,l,k) t1_1_aa(c,j) t1_1_aa(d,i)
        //            += +0.50 <l,k||c,d>_aaaa t2_1_aaaa(a,b,l,k) t2_1_aaaa(c,d,i,j)
        //            += +4.00 t2_2_aaaa(a,b,i,j) w0
        //            += +1.00 <l,k||i,j>_aaaa t2_2_aaaa(a,b,l,k)
        //            += -2.00 <a,b||c,d>_aaaa t1_1_aa(c,j) t1_1_aa(d,i)
        //            += +1.00 <a,b||c,d>_aaaa t2_2_aaaa(c,d,i,j)
        //            += -4.00 d-_bb(k,c) t1_1_bb(c,k) t2_2_aaaa(a,b,i,j)
        //            += -4.00 d-_aa(k,c) t1_1_aa(c,k) t2_2_aaaa(a,b,i,j)
        //            += -2.00 d-_aa(k,c) t2_1_aaaa(a,b,i,j) t1_2_aa(c,k)
        //            += -2.00 d-_bb(k,c) t2_1_aaaa(a,b,i,j) t1_2_bb(c,k)
        //            += -2.00 <l,k||i,j>_aaaa t1_1_aa(a,k) t1_1_aa(b,l)
        //            += -1.00 <l,k||c,d>_aaaa t2_aaaa(c,d,i,j) t1_1_aa(a,k) t1_1_aa(b,l)
        //            += -1.00 <l,k||c,d>_aaaa t2_1_aaaa(c,a,i,j) t2_1_aaaa(d,b,l,k)
        //            += +0.50 <l,k||c,d>_aaaa t2_aaaa(a,b,l,k) t2_2_aaaa(c,d,i,j)
        //            += +0.50 <l,k||c,d>_aaaa t2_aaaa(c,d,i,j) t2_2_aaaa(a,b,l,k)
        //            += -1.00 <l,k||c,d>_aaaa t2_1_aaaa(c,a,l,k) t2_1_aaaa(d,b,i,j)
        rt2_2_aaaa("a,b,i,j") += tmps_["149_aaaa_oovv"]("i,j,a,b");
        tmps_["149_aaaa_oovv"].~TArrayD();

        // flops: o2v2  = o3v3
        //  mems: o2v2  = o2v2
        tmps_["150_aaaa_vovo"]("a,i,b,j")  = t2_1["abab"]("a,c,i,k") * tmps_["28_bbaa_ovvo"]("k,c,b,j");

        // rt2_2_aaaa += +2.00 P(i,j) <k,l||c,d>_abab t2_1_aaaa(c,a,j,k) t2_1_abab(b,d,i,l)
        rt2_2_aaaa("a,b,i,j") += 2.00 * tmps_["150_aaaa_vovo"]("b,i,a,j");
        rt2_2_aaaa("a,b,i,j") -= 2.00 * tmps_["150_aaaa_vovo"]("b,j,a,i");

        // rt2_2_aaaa += +2.00 P(i,j) <l,k||d,c>_abab t2_1_abab(a,c,j,k) t2_1_aaaa(d,b,i,l)
        rt2_2_aaaa("a,b,i,j") += 2.00 * tmps_["150_aaaa_vovo"]("a,j,b,i");
        rt2_2_aaaa("a,b,i,j") -= 2.00 * tmps_["150_aaaa_vovo"]("a,i,b,j");
        tmps_["150_aaaa_vovo"].~TArrayD();

        // flops: o2v2  = o3v2 o3v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2
        tmps_["151_aaaa_vvoo"]("a,b,i,j")  = -1.00 * tmps_["9_aa_oo"]("k,j") * t2["aaaa"]("a,b,i,k");
        tmps_["151_aaaa_vvoo"]("a,b,i,j") += dp["aa_oo"]("k,i") * t2_1["aaaa"]("a,b,j,k");

        // rt2_2_aaaa += +2.00 P(i,j) d+_aa(k,j) t2_1_aaaa(a,b,i,k)
        //            += -2.00 P(i,j) d+_aa(k,c) t2_aaaa(a,b,j,k) t1_1_aa(c,i)
        rt2_2_aaaa("a,b,i,j") += 2.00 * tmps_["151_aaaa_vvoo"]("a,b,j,i");
        rt2_2_aaaa("a,b,i,j") -= 2.00 * tmps_["151_aaaa_vvoo"]("a,b,i,j");

        // rt2_aaaa += +1.00 P(i,j) d-_aa(k,j) t2_1_aaaa(a,b,i,k)
        //          += -1.00 P(i,j) d-_aa(k,c) t2_aaaa(a,b,j,k) t1_1_aa(c,i)
        rt2_aaaa("a,b,i,j") += tmps_["151_aaaa_vvoo"]("a,b,j,i");
        rt2_aaaa("a,b,i,j") -= tmps_["151_aaaa_vvoo"]("a,b,i,j");
        tmps_["151_aaaa_vvoo"].~TArrayD();

        // flops: o2v2  = o2v3
        //  mems: o2v2  = o2v2
        tmps_["152_aaaa_voov"]("a,i,j,b")  = tmps_["32_aa_vv"]("c,b") * t2_1["aaaa"]("c,a,i,j");

        // rt2_2_aaaa += -1.00 <l,k||d,c>_abab t2_1_abab(a,c,l,k) t2_1_aaaa(d,b,i,j)
        //            += -1.00 <k,l||d,c>_abab t2_1_abab(a,c,k,l) t2_1_aaaa(d,b,i,j)
        rt2_2_aaaa("a,b,i,j") -= 2.00 * tmps_["152_aaaa_voov"]("b,i,j,a");

        // rt2_2_aaaa += +1.00 <l,k||c,d>_abab t2_1_aaaa(c,a,i,j) t2_1_abab(b,d,l,k)
        //            += +1.00 <k,l||c,d>_abab t2_1_aaaa(c,a,i,j) t2_1_abab(b,d,k,l)
        rt2_2_aaaa("a,b,i,j") += 2.00 * tmps_["152_aaaa_voov"]("a,i,j,b");
        tmps_["152_aaaa_voov"].~TArrayD();

        // flops: o1v1  = o2v2 o2v2 o1v1
        //  mems: o1v1  = o1v1 o1v1 o1v1
        tmps_["159_aa_ov"]("i,a")  = -1.00 * eri["abab_oovv"]("i,k,a,c") * t1_1["bb"]("c,k");
        tmps_["159_aa_ov"]("i,a") += eri["aaaa_oovv"]("j,i,a,b") * t1_1["aa"]("b,j");

        // flops: o3v1  = o3v2
        //  mems: o3v1  = o3v1
        tmps_["160_aaaa_vooo"]("a,i,j,k")  = tmps_["159_aa_ov"]("k,b") * t2["aaaa"]("b,a,i,j");

        // flops: o3v1  = o3v2 o3v2 o3v1
        //  mems: o3v1  = o3v1 o3v1 o3v1
        tmps_["158_aaaa_ovoo"]("i,a,j,k")  = -0.33 * f["aa_ov"]("i,b") * t2_1["aaaa"]("b,a,j,k");
        tmps_["158_aaaa_ovoo"]("i,a,j,k") += dp["aa_ov"]("i,b") * t2_2["aaaa"]("b,a,j,k");

        // flops: o1v1  = o1v1 o1v1 o1v1
        //  mems: o1v1  = o1v1 o1v1 o1v1
        tmps_["157_aa_vo"]("a,i")  = 2.00 * t0_2 * t1_1["aa"]("a,i");
        tmps_["157_aa_vo"]("a,i") += t0_1 * t1_2["aa"]("a,i");

        // flops: o0v2  = o1v3 o1v3 o0v2
        //  mems: o0v2  = o0v2 o0v2 o0v2
        tmps_["156_aa_vv"]("a,b")  = eri["aaaa_vovv"]("a,i,b,c") * t1_2["aa"]("c,i");
        tmps_["156_aa_vv"]("a,b") += eri["abab_vovv"]("a,j,b,d") * t1_2["bb"]("d,j");

        // flops: o0v2  = o2v3 o2v3 o0v2
        //  mems: o0v2  = o0v2 o0v2 o0v2
        tmps_["155_aa_vv"]("a,b")  = 0.50 * eri["aaaa_oovv"]("k,i,a,d") * t2_2["aaaa"]("d,b,i,k");
        tmps_["155_aa_vv"]("a,b") += eri["abab_oovv"]("i,j,a,c") * t2_2["abab"]("b,c,i,j");

        // flops: o3v1  = o3v3
        //  mems: o3v1  = o3v1
        tmps_["154_aaaa_vooo"]("a,i,j,k")  = eri["aaaa_vovv"]("a,i,b,c") * t2_1["aaaa"]("b,c,j,k");

        // flops: o0v2  = o1v3
        //  mems: o0v2  = o0v2
        tmps_["153_aa_vv"]("a,b")  = eri["aaaa_vovv"]("a,i,c,b") * t1_1["aa"]("c,i");

        // flops: o2v2  = o3v2 o2v3 o3v2 o2v3 o2v2 o2v2 o2v3 o2v2 o3v2 o2v2 o2v3 o2v2 o2v3 o3v2 o2v2 o2v2 o3v2 o2v2 o2v2 o2v3 o2v2 o2v3 o3v2 o2v2 o3v2 o2v2 o2v2 o3v2 o2v2 o3v2 o2v2 o2v3 o2v2 o2v3 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["161_aaaa_vvoo"]("a,b,i,j")  = -2.00 * tmps_["70_aaaa_ovoo"]("m,b,i,j") * t1_2["aa"]("a,m");
        tmps_["161_aaaa_vvoo"]("a,b,i,j") -= 2.00 * t2["aaaa"]("e,b,i,j") * tmps_["155_aa_vv"]("e,a");
        tmps_["161_aaaa_vvoo"]("a,b,i,j") += 2.00 * tmps_["68_aaaa_ovoo"]("m,b,i,j") * tmps_["157_aa_vo"]("a,m");
        tmps_["161_aaaa_vvoo"]("a,b,i,j") += 2.00 * t2_2["aaaa"]("e,a,i,j") * tmps_["10_aa_vv"]("b,e");
        tmps_["161_aaaa_vvoo"]("a,b,i,j") -= 2.00 * t2["aaaa"]("e,a,i,j") * tmps_["156_aa_vv"]("b,e");
        tmps_["161_aaaa_vvoo"]("a,b,i,j") -= 2.00 * tmps_["12_aa_vo"]("b,m") * tmps_["69_aaaa_ovoo"]("m,a,i,j");
        tmps_["161_aaaa_vvoo"]("a,b,i,j") += 4.00 * t2_1["aaaa"]("e,a,i,j") * tmps_["40_aa_vv"]("b,e");
        tmps_["161_aaaa_vvoo"]("a,b,i,j") -= 2.00 * f["aa_vv"]("b,e") * t2_2["aaaa"]("e,a,i,j");
        tmps_["161_aaaa_vvoo"]("a,b,i,j") += 2.00 * eri["aaaa_vooo"]("b,m,i,j") * t1_2["aa"]("a,m");
        tmps_["161_aaaa_vvoo"]("a,b,i,j") -= 6.00 * t1_1["aa"]("b,m") * tmps_["158_aaaa_ovoo"]("m,a,i,j");
        tmps_["161_aaaa_vvoo"]("a,b,i,j") -= tmps_["4_aa_vv"]("c,b") * t2_2["aaaa"]("c,a,i,j");
        tmps_["161_aaaa_vvoo"]("a,b,i,j") -= 2.00 * tmps_["72_aa_vv"]("b,c") * t2_1["aaaa"]("c,a,i,j");
        tmps_["161_aaaa_vvoo"]("a,b,i,j") += tmps_["154_aaaa_vooo"]("b,m,i,j") * t1_1["aa"]("a,m");
        tmps_["161_aaaa_vvoo"]("a,b,i,j") += tmps_["67_aaaa_vooo"]("b,m,i,j") * t1_2["aa"]("a,m");
        tmps_["161_aaaa_vvoo"]("a,b,i,j") += 2.00 * tmps_["160_aaaa_vooo"]("b,i,j,k") * t1_1["aa"]("a,k");
        tmps_["161_aaaa_vvoo"]("a,b,i,j") += 6.00 * tmps_["69_aaaa_ovoo"]("m,b,i,j") * t1_2["aa"]("a,m");
        tmps_["161_aaaa_vvoo"]("a,b,i,j") += 2.00 * tmps_["153_aa_vv"]("b,c") * t2_1["aaaa"]("c,a,i,j");
        tmps_["161_aaaa_vvoo"]("a,b,i,j") += 2.00 * tmps_["5_aa_vv"]("c,b") * t2_2["aaaa"]("c,a,i,j");
        tmps_["160_aaaa_vooo"].~TArrayD();
        tmps_["158_aaaa_ovoo"].~TArrayD();
        tmps_["154_aaaa_vooo"].~TArrayD();
        tmps_["70_aaaa_ovoo"].~TArrayD();
        tmps_["69_aaaa_ovoo"].~TArrayD();
        tmps_["68_aaaa_ovoo"].~TArrayD();
        tmps_["67_aaaa_vooo"].~TArrayD();

        // rt2_2_aaaa += +1.00 P(a,b) <k,a||c,d>_aaaa t1_1_aa(b,k) t2_1_aaaa(c,d,i,j)
        //            += +2.00 P(a,b) <a,k||d,c>_abab t2_1_aaaa(d,b,i,j) t1_1_bb(c,k)
        //            += +1.00 P(a,b) <k,a||c,d>_aaaa t2_aaaa(c,d,i,j) t1_2_aa(b,k)
        //            += -2.00 P(a,b) d-_aa(a,c) t0_1 t2_2_aaaa(c,b,i,j)
        //            += -2.00 P(a,b) d-_aa(k,c) t2_aaaa(c,a,i,j) t0_1 t1_2_aa(b,k)
        //            += -4.00 P(a,b) d-_aa(k,c) t2_aaaa(c,a,i,j) t1_1_aa(b,k) t0_2
        //            += +1.00 P(a,b) <l,k||c,d>_abab t2_aaaa(c,a,i,j) t2_2_abab(b,d,l,k)
        //            += +1.00 P(a,b) <k,l||c,d>_abab t2_aaaa(c,a,i,j) t2_2_abab(b,d,k,l)
        //            += -1.00 P(a,b) <l,k||c,d>_aaaa t2_aaaa(c,a,i,j) t2_2_aaaa(d,b,l,k)
        //            += +2.00 P(a,b) <a,k||c,d>_abab t2_aaaa(c,b,i,j) t1_2_bb(d,k)
        //            += -2.00 P(a,b) <k,a||c,d>_aaaa t2_aaaa(c,b,i,j) t1_2_aa(d,k)
        //            += -4.00 P(a,b) d-_aa(a,c) t2_1_aaaa(c,b,i,j) t0_2
        //            += +2.00 P(a,b) d-_aa(k,c) t0_1 t1_1_aa(a,k) t2_1_aaaa(c,b,i,j)
        //            += +6.00 P(a,b) d-_aa(k,c) t1_1_aa(a,k) t2_2_aaaa(c,b,i,j)
        //            += -2.00 P(a,b) f_aa(k,c) t1_1_aa(a,k) t2_1_aaaa(c,b,i,j)
        //            += +2.00 P(a,b) <k,a||i,j>_aaaa t1_2_aa(b,k)
        //            += +2.00 P(a,b) f_aa(a,c) t2_2_aaaa(c,b,i,j)
        //            += +2.00 P(a,b) <l,k||c,d>_aaaa t2_aaaa(c,a,i,j) t1_1_aa(b,l) t1_1_aa(d,k)
        //            += +2.00 P(a,b) <l,k||c,d>_abab t2_aaaa(c,a,i,j) t1_1_aa(b,l) t1_1_bb(d,k)
        //            += -6.00 P(a,b) d-_aa(k,c) t2_1_aaaa(c,a,i,j) t1_2_aa(b,k)
        //            += +2.00 P(a,b) f_aa(k,c) t2_aaaa(c,a,i,j) t1_2_aa(b,k)
        //            += -1.00 P(a,b) <l,k||c,d>_aaaa t2_aaaa(c,a,l,k) t2_2_aaaa(d,b,i,j)
        //            += +2.00 P(a,b) <k,a||c,d>_aaaa t2_1_aaaa(d,b,i,j) t1_1_aa(c,k)
        //            += -1.00 P(a,b) <l,k||d,c>_abab t2_abab(a,c,l,k) t2_2_aaaa(d,b,i,j)
        //            += -1.00 P(a,b) <k,l||d,c>_abab t2_abab(a,c,k,l) t2_2_aaaa(d,b,i,j)
        rt2_2_aaaa("a,b,i,j") -= tmps_["161_aaaa_vvoo"]("b,a,i,j");
        rt2_2_aaaa("a,b,i,j") += tmps_["161_aaaa_vvoo"]("a,b,i,j");
        tmps_["161_aaaa_vvoo"].~TArrayD();

        // flops: o1v1  = o2v2
        //  mems: o1v1  = o1v1
        tmps_["169_bb_ov"]("i,a")  = eri["bbbb_oovv"]("j,i,a,b") * t1_1["bb"]("b,j");

        // flops: o3v1  = o4v2 o1v1 o3v2 o3v1
        //  mems: o3v1  = o3v1 o1v1 o3v1 o3v1
        tmps_["190_abba_vooo"]("a,i,j,k")  = -1.00 * tmps_["124_abba_oovo"]("l,j,b,k") * t2["abab"]("a,b,l,i");
        tmps_["190_abba_vooo"]("a,i,j,k") += t2["abab"]("a,b,k,i") * (-1.00 * tmps_["169_bb_ov"]("j,b") + tmps_["8_bb_ov"]("j,b"));

        // flops: o3v1  = o4v2
        //  mems: o3v1  = o3v1
        tmps_["189_bbaa_vooo"]("a,i,j,k")  = t2["bbbb"]("b,a,i,l") * tmps_["124_abba_oovo"]("j,l,b,k");

        // flops: o3v1  = o4v1
        //  mems: o3v1  = o3v1
        tmps_["188_baab_vooo"]("a,i,j,k")  = t1_1["bb"]("a,l") * tmps_["88_abab_oooo"]("i,l,j,k");

        // flops: o4v0  = o4v1
        //  mems: o4v0  = o4v0
        tmps_["187_baba_oooo"]("i,j,k,l")  = t1_1["bb"]("a,i") * tmps_["124_abba_oovo"]("j,k,a,l");
        tmps_["124_abba_oovo"].~TArrayD();

        // flops: o2v0  = o2v1
        //  mems: o2v0  = o2v0
        tmps_["186_bb_oo"]("i,j")  = t1_1["bb"]("a,i") * tmps_["8_bb_ov"]("j,a");

        // flops: o1v1  = o2v1
        //  mems: o1v1  = o1v1
        tmps_["185_bb_vo"]("a,i")  = t1_1["bb"]("a,j") * tmps_["24_bb_oo"]("j,i");

        // flops: o3v1  = o2v2 o3v2
        //  mems: o3v1  = o2v2 o3v1
        tmps_["184_bbaa_oovo"]("i,j,a,k")  = t1_1["bb"]("b,i") * (-1.00 * tmps_["3_bbaa_ovvo"]("j,b,a,k") + tmps_["2_bbaa_ovvo"]("j,b,a,k"));

        // flops: o3v1  = o3v2 o3v2 o3v1 o3v2 o3v1
        //  mems: o3v1  = o3v1 o3v1 o3v1 o3v1 o3v1
        tmps_["183_baba_vooo"]("a,i,j,k")  = tmps_["159_aa_ov"]("k,c") * t2["abab"]("c,a,i,j");
        tmps_["183_baba_vooo"]("a,i,j,k") -= tmps_["17_aabb_ovvo"]("k,d,a,j") * t1_1["aa"]("d,i");
        tmps_["183_baba_vooo"]("a,i,j,k") += t1_1["bb"]("b,j") * tmps_["86_abba_ovvo"]("k,b,a,i");
        tmps_["159_aa_ov"].~TArrayD();

        // flops: o3v1  = o3v2
        //  mems: o3v1  = o3v1
        tmps_["168_bbbb_oovo"]("i,j,a,k")  = eri["bbbb_oovv"]("i,j,a,b") * t1_1["bb"]("b,k");

        // flops: o2v0  = o2v0 o2v0 o2v0 o3v1 o2v0
        //  mems: o2v0  = o2v0 o2v0 o2v0 o2v0 o2v0
        tmps_["182_bb_oo"]("i,j")  = 2.00 * t0_2 * tmps_["24_bb_oo"]("i,j");
        tmps_["182_bb_oo"]("i,j") += t0_1 * tmps_["58_bb_oo"]("i,j");
        tmps_["182_bb_oo"]("i,j") += t1_1["bb"]("a,k") * tmps_["168_bbbb_oovo"]("i,k,a,j");

        // flops: o2v2  = o3v3 o2v3 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2
        tmps_["181_abba_ovvo"]("i,a,b,j")  = eri["abab_oovv"]("i,k,c,a") * t2_2["abab"]("c,b,j,k");
        tmps_["181_abba_ovvo"]("i,a,b,j") += eri["baab_vovv"]("b,i,c,a") * t1_2["aa"]("c,j");

        // flops: o3v1  = o3v2 o4v2 o3v1 o4v2 o3v1 o3v3 o3v1 o3v2 o3v1 o4v2 o3v1 o3v2 o3v1 o3v2 o3v1
        //  mems: o3v1  = o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1
        tmps_["180_abab_vooo"]("a,i,j,k")  = -1.00 * eri["abba_vovo"]("a,i,b,j") * t1_1["bb"]("b,k");
        tmps_["180_abab_vooo"]("a,i,j,k") += eri["abba_oovo"]("m,i,b,j") * t2_1["abab"]("a,b,m,k");
        tmps_["180_abab_vooo"]("a,i,j,k") -= eri["abab_oovo"]("m,i,c,k") * t2_1["aaaa"]("c,a,j,m");
        tmps_["180_abab_vooo"]("a,i,j,k") += eri["abab_vovv"]("a,i,c,d") * t2_1["abab"]("c,d,j,k");
        tmps_["180_abab_vooo"]("a,i,j,k") -= 3.00 * dp["bb_ov"]("i,b") * t2_2["abab"]("a,b,j,k");
        tmps_["180_abab_vooo"]("a,i,j,k") -= eri["bbbb_oovo"]("i,l,b,k") * t2_1["abab"]("a,b,j,l");
        tmps_["180_abab_vooo"]("a,i,j,k") += eri["abab_vovo"]("a,i,c,k") * t1_1["aa"]("c,j");
        tmps_["180_abab_vooo"]("a,i,j,k") += f["bb_ov"]("i,b") * t2_1["abab"]("a,b,j,k");

        // flops: o3v1  = o3v2 o3v2 o3v2 o3v1 o3v1 o4v2 o3v1 o3v2 o3v1 o3v3 o3v1 o4v1 o4v2 o3v1 o3v1 o3v1 o4v2 o3v1
        //  mems: o3v1  = o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1 o3v1
        tmps_["179_aabb_ooov"]("i,j,k,a")  = -3.00 * dp["aa_ov"]("i,b") * t2_2["abab"]("b,a,j,k");
        tmps_["179_aabb_ooov"]("i,j,k,a") += eri["baba_vovo"]("a,i,d,j") * t1_1["bb"]("d,k");
        tmps_["179_aabb_ooov"]("i,j,k,a") += f["aa_ov"]("i,b") * t2_1["abab"]("b,a,j,k");
        tmps_["179_aabb_ooov"]("i,j,k,a") -= eri["aaaa_oovo"]("i,m,b,j") * t2_1["abab"]("b,a,m,k");
        tmps_["179_aabb_ooov"]("i,j,k,a") -= eri["baab_vovo"]("a,i,b,k") * t1_1["aa"]("b,j");
        tmps_["179_aabb_ooov"]("i,j,k,a") -= eri["baab_vovv"]("a,i,b,c") * t2_1["abab"]("b,c,j,k");
        tmps_["179_aabb_ooov"]("i,j,k,a") += -1.00 * (eri["abab_oooo"]("i,l,j,k") * t1_1["bb"]("a,l") + eri["abab_oovo"]("i,l,b,k") * t2_1["abab"]("b,a,j,l"));
        tmps_["179_aabb_ooov"]("i,j,k,a") += eri["abba_oovo"]("i,l,d,j") * t2_1["bbbb"]("d,a,k,l");

        // flops: o4v0  = o4v1 o4v1 o4v2 o4v0 o4v0
        //  mems: o4v0  = o4v0 o4v0 o4v0 o4v0 o4v0
        tmps_["178_abab_oooo"]("i,j,k,l")  = -1.00 * eri["abba_oovo"]("i,j,a,k") * t1_2["bb"]("a,l");
        tmps_["178_abab_oooo"]("i,j,k,l") += eri["abab_oovo"]("i,j,b,l") * t1_2["aa"]("b,k");
        tmps_["178_abab_oooo"]("i,j,k,l") += eri["abab_oovv"]("i,j,b,c") * t2_2["abab"]("b,c,k,l");

        // flops: o1v1  = o1v1 o1v1 o1v1
        //  mems: o1v1  = o1v1 o1v1 o1v1
        tmps_["177_bb_vo"]("a,i")  = 0.50 * t0_1 * t1_2["bb"]("a,i");
        tmps_["177_bb_vo"]("a,i") += t0_2 * t1_1["bb"]("a,i");

        // flops: o2v0  = o3v2 o2v1 o2v0 o3v2 o2v0
        //  mems: o2v0  = o2v0 o2v0 o2v0 o2v0 o2v0
        tmps_["176_bb_oo"]("i,j")  = 0.50 * eri["bbbb_oovv"]("i,l,c,b") * t2_2["bbbb"]("c,b,j,l");
        tmps_["176_bb_oo"]("i,j") += f["bb_ov"]("i,c") * t1_2["bb"]("c,j");
        tmps_["176_bb_oo"]("i,j") += eri["abab_oovv"]("k,i,a,b") * t2_2["abab"]("a,b,k,j");

        // flops: o2v0  = o3v1 o3v1 o2v0
        //  mems: o2v0  = o2v0 o2v0 o2v0
        tmps_["175_bb_oo"]("i,j")  = -1.00 * eri["bbbb_oovo"]("i,l,b,j") * t1_2["bb"]("b,l");
        tmps_["175_bb_oo"]("i,j") += eri["abab_oovo"]("k,i,a,j") * t1_2["aa"]("a,k");

        // flops: o1v1  = o2v2 o2v2 o1v1
        //  mems: o1v1  = o1v1 o1v1 o1v1
        tmps_["174_bb_vo"]("a,i")  = -1.00 * dp["aa_ov"]("k,c") * t2_2["abab"]("c,a,k,i");
        tmps_["174_bb_vo"]("a,i") += dp["bb_ov"]("j,b") * t2_2["bbbb"]("b,a,i,j");

        // flops: o0v2  = o1v3 o1v3 o0v2
        //  mems: o0v2  = o0v2 o0v2 o0v2
        tmps_["173_bb_vv"]("a,b")  = -1.00 * eri["bbbb_vovv"]("a,j,b,d") * t1_2["bb"]("d,j");
        tmps_["173_bb_vv"]("a,b") += eri["baab_vovv"]("a,i,c,b") * t1_2["aa"]("c,i");

        // flops: o0v2  = o2v3 o2v3 o0v2
        //  mems: o0v2  = o0v2 o0v2 o0v2
        tmps_["172_bb_vv"]("a,b")  = 0.50 * eri["bbbb_oovv"]("j,k,a,d") * t2_2["bbbb"]("d,b,k,j");
        tmps_["172_bb_vv"]("a,b") += eri["abab_oovv"]("i,j,c,a") * t2_2["abab"]("c,b,i,j");

        // flops: o1v1  = o2v1
        //  mems: o1v1  = o1v1
        tmps_["171_bb_ov"]("i,a")  = dp["bb_oo"]("j,i") * t1_2["bb"]("a,j");

        // flops: o1v1  = o1v2
        //  mems: o1v1  = o1v1
        tmps_["170_bb_vo"]("a,i")  = dp["bb_vv"]("a,b") * t1_2["bb"]("b,i");

        // flops: o0v2  = o1v3
        //  mems: o0v2  = o0v2
        tmps_["167_bb_vv"]("a,b")  = eri["bbbb_vovv"]("a,i,c,b") * t1_1["bb"]("c,i");

        // flops: o2v0  = o3v2
        //  mems: o2v0  = o2v0
        tmps_["166_bb_oo"]("i,j")  = eri["bbbb_oovv"]("k,i,a,b") * t2_1["bbbb"]("a,b,j,k");

        // flops: o2v2  = o2v3
        //  mems: o2v2  = o2v2
        tmps_["165_bbbb_vovo"]("a,i,b,j")  = eri["bbbb_vovv"]("a,i,b,c") * t1_2["bb"]("c,j");

        // flops: o2v2  = o2v3
        //  mems: o2v2  = o2v2
        tmps_["164_bbbb_vovo"]("a,i,b,j")  = eri["bbbb_vovv"]("a,i,c,b") * t1_1["bb"]("c,j");

        // flops: o2v2  = o2v3
        //  mems: o2v2  = o2v2
        tmps_["163_baab_vovo"]("a,i,b,j")  = eri["baab_vovv"]("a,i,b,c") * t1_2["bb"]("c,j");

        // flops: o2v2  = o3v3
        //  mems: o2v2  = o2v2
        tmps_["162_bbbb_ovvo"]("i,a,b,j")  = eri["abab_oovv"]("k,i,c,a") * t2_2["abab"]("c,b,k,j");
}
#endif
