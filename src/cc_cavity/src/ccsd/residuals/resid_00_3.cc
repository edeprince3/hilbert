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

#include "cc_cavity/include/ccsd/ccsd.h"


using namespace std;
using namespace TA;
using namespace hilbert;

void CCSD::resid_00_3() {

    // unpack integrals
    TArrayMap &Id         = Id_blks_;
    TArrayMap &f          = F_blks_;
    TArrayMap &eri        = V_blks_;

    // extract energy
    double &energy  = scalars_["energy"];

    // extract residuals
    TArrayD &rt1_aa     = residuals_["t1_aa"];
    TArrayD &rt1_bb     = residuals_["t1_bb"];
    TArrayD &rt2_aaaa   = residuals_["t2_aaaa"];
    TArrayD &rt2_abab   = residuals_["t2_abab"];
    TArrayD &rt2_bbbb   = residuals_["t2_bbbb"];

    // extract 1-body amplitudes
    std::map<std::string, TA::TArrayD> t1 {
            {"aa_vo", amplitudes_["t1_aa"]},
            {"bb_vo", amplitudes_["t1_bb"]}
    };

    // extract 2-body amplitudes
    std::map<std::string, TA::TArrayD> t2 {
            {"aaaa_vvoo", amplitudes_["t2_aaaa"]},
            {"abab_vvoo", amplitudes_["t2_abab"]},
            {"bbbb_vvoo", amplitudes_["t2_bbbb"]}
    };

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["15_aaaa_vvoo"]("a,b,j,i")  = t2["abab_vvoo"]("b,d,i,l") * tmps_["1_abab_vvoo"]("a,d,j,l");
    tmps_["1_abab_vvoo"].~TArrayD();

    // rt2_aaaa += +1.00 P(i,j) <l,k||d,c>_abab t2_abab(a,c,j,k) t2_aaaa(d,b,i,l)
    rt2_aaaa("a,b,i,j") += tmps_["15_aaaa_vvoo"]("b,a,i,j");
    rt2_aaaa("a,b,i,j") -= tmps_["15_aaaa_vvoo"]("b,a,j,i");

    // rt2_aaaa += +1.00 P(i,j) <k,l||c,d>_abab t2_aaaa(c,a,j,k) t2_abab(b,d,i,l)
    rt2_aaaa("a,b,i,j") += tmps_["15_aaaa_vvoo"]("a,b,j,i");
    rt2_aaaa("a,b,i,j") -= tmps_["15_aaaa_vvoo"]("a,b,i,j");
    tmps_["15_aaaa_vvoo"].~TArrayD();

    // flops: o2v2  = o3v3 o3v3 o2v2 o3v2 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["21_bbbb_vvoo"]("b,a,i,j")  = t2["abab_vvoo"]("c,a,k,j") * tmps_["2_abab_vvoo"]("c,b,k,i");
    tmps_["21_bbbb_vvoo"]("b,a,i,j") += t2["bbbb_vvoo"]("c,a,j,k") * tmps_["4_bbbb_vvoo"]("b,c,i,k");
    tmps_["21_bbbb_vvoo"]("b,a,i,j") -= 0.50 * t2["bbbb_vvoo"]("a,b,i,l") * tmps_["10_bb_oo"]("l,j");
    tmps_["10_bb_oo"].~TArrayD();
    tmps_["4_bbbb_vvoo"].~TArrayD();
    tmps_["2_abab_vvoo"].~TArrayD();

    // rt2_bbbb += +1.00 P(i,j) <l,k||c,d>_aaaa t2_abab(c,a,k,j) t2_abab(d,b,l,i)
    //          += +1.00 P(i,j) <l,k||c,d>_bbbb t2_bbbb(c,a,j,k) t2_bbbb(d,b,i,l)
    //          += -0.50 P(i,j) <l,k||c,d>_bbbb t2_bbbb(a,b,i,l) t2_bbbb(c,d,j,k)
    //          += -0.50 P(i,j) <k,l||c,d>_abab t2_bbbb(a,b,i,l) t2_abab(c,d,k,j)
    //          += -0.50 P(i,j) <k,l||d,c>_abab t2_bbbb(a,b,i,l) t2_abab(d,c,k,j)
    rt2_bbbb("a,b,i,j") -= tmps_["21_bbbb_vvoo"]("b,a,i,j");
    rt2_bbbb("a,b,i,j") += tmps_["21_bbbb_vvoo"]("b,a,j,i");
    tmps_["21_bbbb_vvoo"].~TArrayD();

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["16_bbbb_vvoo"]("a,b,j,i")  = t2["bbbb_vvoo"]("d,b,i,l") * tmps_["3_bbbb_vvoo"]("a,d,j,l");
    tmps_["3_bbbb_vvoo"].~TArrayD();

    // rt2_bbbb += +1.00 P(i,j) <l,k||d,c>_abab t2_bbbb(c,a,j,k) t2_abab(d,b,l,i)
    rt2_bbbb("a,b,i,j") += tmps_["16_bbbb_vvoo"]("b,a,i,j");
    rt2_bbbb("a,b,i,j") -= tmps_["16_bbbb_vvoo"]("b,a,j,i");

    // rt2_bbbb += +1.00 P(i,j) <k,l||c,d>_abab t2_abab(c,a,k,j) t2_bbbb(d,b,i,l)
    rt2_bbbb("a,b,i,j") += tmps_["16_bbbb_vvoo"]("a,b,j,i");
    rt2_bbbb("a,b,i,j") -= tmps_["16_bbbb_vvoo"]("a,b,i,j");
    tmps_["16_bbbb_vvoo"].~TArrayD();

    // flops: o2v2  = o2v3 o2v3
    //  mems: o2v2  = o0v2 o2v2
    tmps_["5_aaaa_vvoo"]("b,a,j,i")  = eri["abab_oovv"]("l,k,d,c") * t2["abab_vvoo"]("a,c,l,k") * t2["aaaa_vvoo"]("d,b,i,j");

    // rt2_aaaa += +0.50 <l,k||c,d>_abab t2_aaaa(c,a,i,j) t2_abab(b,d,l,k)
    //          += +0.50 <k,l||c,d>_abab t2_aaaa(c,a,i,j) t2_abab(b,d,k,l)
    rt2_aaaa("a,b,i,j") += tmps_["5_aaaa_vvoo"]("a,b,j,i");

    // rt2_aaaa += -0.50 <l,k||d,c>_abab t2_abab(a,c,l,k) t2_aaaa(d,b,i,j)
    //          += -0.50 <k,l||d,c>_abab t2_abab(a,c,k,l) t2_aaaa(d,b,i,j)
    rt2_aaaa("a,b,i,j") -= tmps_["5_aaaa_vvoo"]("b,a,j,i");
    tmps_["5_aaaa_vvoo"].~TArrayD();

    // flops: o2v2  = o2v3 o2v4 o4v2 o2v2 o2v2 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["19_aaaa_vvoo"]("b,a,i,j")  = t2["aaaa_vvoo"]("d,b,i,j") * tmps_["6_aa_vv"]("a,d");
    tmps_["19_aaaa_vvoo"]("b,a,i,j") += eri["aaaa_vvvv"]("a,b,c,d") * t2["aaaa_vvoo"]("c,d,i,j");
    tmps_["19_aaaa_vvoo"]("b,a,i,j") += 2.00 * eri["aaaa_vvoo"]("a,b,i,j");
    tmps_["19_aaaa_vvoo"]("b,a,i,j") -= eri["aaaa_oooo"]("k,l,i,j") * t2["aaaa_vvoo"]("a,b,l,k");
    tmps_["6_aa_vv"].~TArrayD();

    // rt2_aaaa += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,a,l,k) t2_aaaa(d,b,i,j)
    //          += +0.50 <a,b||c,d>_aaaa t2_aaaa(c,d,i,j)
    //          += +1.00 <a,b||i,j>_aaaa
    //          += +0.50 <l,k||i,j>_aaaa t2_aaaa(a,b,l,k)
    rt2_aaaa("a,b,i,j") += 0.50 * tmps_["19_aaaa_vvoo"]("b,a,i,j");
    tmps_["19_aaaa_vvoo"].~TArrayD();

    // flops: o2v2  = o2v3
    //  mems: o2v2  = o2v2
    tmps_["17_bbbb_vvoo"]("b,a,j,i")  = t2["bbbb_vvoo"]("c,a,i,j") * tmps_["7_bb_vv"]("b,c");
    tmps_["7_bb_vv"].~TArrayD();

    // rt2_bbbb += +0.50 <l,k||d,c>_abab t2_bbbb(c,a,i,j) t2_abab(d,b,l,k)
    //          += +0.50 <k,l||d,c>_abab t2_bbbb(c,a,i,j) t2_abab(d,b,k,l)
    rt2_bbbb("a,b,i,j") += tmps_["17_bbbb_vvoo"]("b,a,j,i");

    // rt2_bbbb += -0.50 <l,k||c,d>_abab t2_abab(c,a,l,k) t2_bbbb(d,b,i,j)
    //          += -0.50 <k,l||c,d>_abab t2_abab(c,a,k,l) t2_bbbb(d,b,i,j)
    rt2_bbbb("a,b,i,j") -= tmps_["17_bbbb_vvoo"]("a,b,j,i");
    tmps_["17_bbbb_vvoo"].~TArrayD();

    // flops: o2v2  = o4v2 o2v4 o2v2 o2v2 o2v3 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["18_bbbb_vvoo"]("b,a,j,i")  = eri["bbbb_oooo"]("k,l,i,j") * t2["bbbb_vvoo"]("a,b,l,k");
    tmps_["18_bbbb_vvoo"]("b,a,j,i") -= eri["bbbb_vvvv"]("a,b,c,d") * t2["bbbb_vvoo"]("c,d,i,j");
    tmps_["18_bbbb_vvoo"]("b,a,j,i") -= 2.00 * eri["bbbb_vvoo"]("a,b,i,j");
    tmps_["18_bbbb_vvoo"]("b,a,j,i") -= t2["bbbb_vvoo"]("c,a,i,j") * tmps_["8_bb_vv"]("b,c");
    tmps_["8_bb_vv"].~TArrayD();

    // rt2_bbbb += +0.50 <l,k||i,j>_bbbb t2_bbbb(a,b,l,k)
    //          += +0.50 <a,b||c,d>_bbbb t2_bbbb(c,d,i,j)
    //          += +1.00 <a,b||i,j>_bbbb
    //          += -0.50 <l,k||c,d>_bbbb t2_bbbb(c,a,i,j) t2_bbbb(d,b,l,k)
    rt2_bbbb("a,b,i,j") -= 0.50 * tmps_["18_bbbb_vvoo"]("b,a,j,i");
    tmps_["18_bbbb_vvoo"].~TArrayD();

}
