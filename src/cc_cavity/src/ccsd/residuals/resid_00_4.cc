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

void CCSD::resid_00_4() {

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

    // flops: o2v2  = o3v2 o3v2 o3v3 o3v3 o2v2 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["22_aaaa_vvoo"]("a,b,i,j")  = t2["aaaa_vvoo"]("a,b,i,l") * tmps_["9_aa_oo"]("l,j");
    tmps_["22_aaaa_vvoo"]("a,b,i,j") += f["aa_oo"]("k,j") * t2["aaaa_vvoo"]("a,b,i,k");
    tmps_["22_aaaa_vvoo"]("a,b,i,j") += eri["aaaa_oovv"]("k,l,c,d") * t2["aaaa_vvoo"]("c,a,j,k") * t2["aaaa_vvoo"]("d,b,i,l");
    tmps_["9_aa_oo"].~TArrayD();

    // rt2_aaaa += -0.50 P(i,j) <l,k||c,d>_abab t2_aaaa(a,b,i,l) t2_abab(c,d,j,k)
    //          += -0.50 P(i,j) <l,k||d,c>_abab t2_aaaa(a,b,i,l) t2_abab(d,c,j,k)
    //          += -0.50 P(i,j) <l,k||c,d>_aaaa t2_aaaa(a,b,i,l) t2_aaaa(c,d,j,k)
    //          += -1.00 P(i,j) f_aa(k,j) t2_aaaa(a,b,i,k)
    //          += +1.00 P(i,j) <l,k||c,d>_aaaa t2_aaaa(c,a,j,k) t2_aaaa(d,b,i,l)
    rt2_aaaa("a,b,i,j") -= tmps_["22_aaaa_vvoo"]("a,b,i,j");
    rt2_aaaa("a,b,i,j") += tmps_["22_aaaa_vvoo"]("a,b,j,i");
    tmps_["22_aaaa_vvoo"].~TArrayD();

    // flops: o2v2  = o3v3 o3v3 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2
    tmps_["11_bbbb_vvoo"]("a,b,j,i")  = eri["baab_vovo"]("a,k,c,j") * t2["abab_vvoo"]("c,b,k,i");
    tmps_["11_bbbb_vvoo"]("a,b,j,i") -= eri["bbbb_vovo"]("a,k,c,j") * t2["bbbb_vvoo"]("c,b,i,k");

    // rt2_bbbb += -1.00 P(i,j) P(a,b) <k,a||c,j>_abab t2_abab(c,b,k,i)
    //          += +1.00 P(i,j) P(a,b) <k,a||c,j>_bbbb t2_bbbb(c,b,i,k)
    rt2_bbbb("a,b,i,j") += tmps_["11_bbbb_vvoo"]("a,b,j,i");
    rt2_bbbb("a,b,i,j") -= tmps_["11_bbbb_vvoo"]("a,b,i,j");
    rt2_bbbb("a,b,i,j") -= tmps_["11_bbbb_vvoo"]("b,a,j,i");
    rt2_bbbb("a,b,i,j") += tmps_["11_bbbb_vvoo"]("b,a,i,j");
    tmps_["11_bbbb_vvoo"].~TArrayD();

    // flops: o2v2  = o3v3 o3v3 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2
    tmps_["12_aaaa_vvoo"]("a,b,j,i")  = eri["aaaa_vovo"]("a,k,c,j") * t2["aaaa_vvoo"]("c,b,i,k");
    tmps_["12_aaaa_vvoo"]("a,b,j,i") -= eri["abba_vovo"]("a,k,c,j") * t2["abab_vvoo"]("b,c,i,k");

    // rt2_aaaa += +1.00 P(i,j) P(a,b) <k,a||c,j>_aaaa t2_aaaa(c,b,i,k)
    //          += -1.00 P(i,j) P(a,b) <a,k||j,c>_abab t2_abab(b,c,i,k)
    rt2_aaaa("a,b,i,j") -= tmps_["12_aaaa_vvoo"]("a,b,j,i");
    rt2_aaaa("a,b,i,j") += tmps_["12_aaaa_vvoo"]("a,b,i,j");
    rt2_aaaa("a,b,i,j") += tmps_["12_aaaa_vvoo"]("b,a,j,i");
    rt2_aaaa("a,b,i,j") -= tmps_["12_aaaa_vvoo"]("b,a,i,j");
    tmps_["12_aaaa_vvoo"].~TArrayD();

    // flops: o1v1  = o2v2 o2v2 o1v1 o2v3 o1v1 o2v3 o1v1 o3v2 o1v1 o3v2 o1v1 o1v1
    //  mems: o1v1  = o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1
    tmps_["13_aa_vo"]("a,i")  = f["aa_ov"]("j,b") * t2["aaaa_vvoo"]("b,a,i,j");
    tmps_["13_aa_vo"]("a,i") -= f["bb_ov"]("j,b") * t2["abab_vvoo"]("a,b,i,j");
    tmps_["13_aa_vo"]("a,i") -= f["aa_vo"]("a,i");
    tmps_["13_aa_vo"]("a,i") -= 0.50 * eri["aaaa_vovv"]("a,j,b,c") * t2["aaaa_vvoo"]("b,c,i,j");
    tmps_["13_aa_vo"]("a,i") -= eri["abab_vovv"]("a,j,b,c") * t2["abab_vvoo"]("b,c,i,j");
    tmps_["13_aa_vo"]("a,i") -= 0.50 * eri["aaaa_oovo"]("j,k,b,i") * t2["aaaa_vvoo"]("b,a,k,j");
    tmps_["13_aa_vo"]("a,i") -= eri["abba_oovo"]("k,j,b,i") * t2["abab_vvoo"]("a,b,k,j");

    // rt1_aa  = -1.00 f_aa(j,b) t2_aaaa(b,a,i,j)
    //        += +1.00 f_bb(j,b) t2_abab(a,b,i,j)
    //        += +1.00 f_aa(a,i)
    //        += -0.50 <j,a||b,c>_aaaa t2_aaaa(b,c,i,j)
    //        += +0.50 <a,j||b,c>_abab t2_abab(b,c,i,j)
    //        += +0.50 <a,j||c,b>_abab t2_abab(c,b,i,j)
    //        += -0.50 <k,j||b,i>_aaaa t2_aaaa(b,a,k,j)
    //        += -0.50 <k,j||i,b>_abab t2_abab(a,b,k,j)
    //        += -0.50 <j,k||i,b>_abab t2_abab(a,b,j,k)
    rt1_aa("a,i")  = -1.00 * tmps_["13_aa_vo"]("a,i");
    tmps_["13_aa_vo"].~TArrayD();

    // flops: o1v1  = o2v2 o2v3 o2v3 o1v1 o3v2 o1v1 o3v2 o1v1 o2v2 o1v1 o1v1 o1v1
    //  mems: o1v1  = o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1 o1v1
    tmps_["14_bb_vo"]("a,i")  = f["aa_ov"]("j,b") * t2["abab_vvoo"]("b,a,j,i");
    tmps_["14_bb_vo"]("a,i") += 0.50 * eri["bbbb_vovv"]("a,j,b,c") * t2["bbbb_vvoo"]("b,c,i,j");
    tmps_["14_bb_vo"]("a,i") -= eri["baab_vovv"]("a,j,b,c") * t2["abab_vvoo"]("b,c,j,i");
    tmps_["14_bb_vo"]("a,i") += f["bb_vo"]("a,i");
    tmps_["14_bb_vo"]("a,i") -= eri["abab_oovo"]("k,j,b,i") * t2["abab_vvoo"]("b,a,k,j");
    tmps_["14_bb_vo"]("a,i") += 0.50 * eri["bbbb_oovo"]("j,k,b,i") * t2["bbbb_vvoo"]("b,a,k,j");
    tmps_["14_bb_vo"]("a,i") -= f["bb_ov"]("j,b") * t2["bbbb_vvoo"]("b,a,i,j");

    // rt1_bb  = +1.00 f_aa(j,b) t2_abab(b,a,j,i)
    //        += -0.50 <j,a||b,c>_bbbb t2_bbbb(b,c,i,j)
    //        += +0.50 <j,a||b,c>_abab t2_abab(b,c,j,i)
    //        += +0.50 <j,a||c,b>_abab t2_abab(c,b,j,i)
    //        += +1.00 f_bb(a,i)
    //        += -0.50 <k,j||b,i>_abab t2_abab(b,a,k,j)
    //        += -0.50 <j,k||b,i>_abab t2_abab(b,a,j,k)
    //        += -0.50 <k,j||b,i>_bbbb t2_bbbb(b,a,k,j)
    //        += -1.00 f_bb(j,b) t2_bbbb(b,a,i,j)
    rt1_bb("a,i")  = tmps_["14_bb_vo"]("a,i");
    tmps_["14_bb_vo"].~TArrayD();

}
