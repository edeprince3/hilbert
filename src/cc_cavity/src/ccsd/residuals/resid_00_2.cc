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

void CCSD::resid_00_2() {

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

    // rt2_bbbb += +0.25 <l,k||c,d>_bbbb t2_bbbb(a,b,l,k) t2_bbbb(c,d,i,j)
    // flops: o2v2 += o4v2 o4v2
    //  mems: o2v2 += o4v0 o2v2
    rt2_bbbb("a,b,i,j") -= 0.25 * eri["bbbb_oovv"]("k,l,c,d") * t2["bbbb_vvoo"]("c,d,i,j") * t2["bbbb_vvoo"]("a,b,l,k");

    // flops: o2v0  = o3v2 o3v2 o2v0
    //  mems: o2v0  = o2v0 o2v0 o2v0
    tmps_["10_bb_oo"]("l,j")  = eri["bbbb_oovv"]("k,l,c,d") * t2["bbbb_vvoo"]("c,d,j,k");
    tmps_["10_bb_oo"]("l,j") -= 2.00 * eri["abab_oovv"]("k,l,c,d") * t2["abab_vvoo"]("c,d,k,j");

    // flops: o2v0  = o3v2 o3v2 o2v0
    //  mems: o2v0  = o2v0 o2v0 o2v0
    tmps_["9_aa_oo"]("l,j")  = eri["abab_oovv"]("l,k,c,d") * t2["abab_vvoo"]("c,d,j,k");
    tmps_["9_aa_oo"]("l,j") -= 0.50 * eri["aaaa_oovv"]("k,l,c,d") * t2["aaaa_vvoo"]("c,d,j,k");

    // flops: o0v2  = o2v3
    //  mems: o0v2  = o0v2
    tmps_["8_bb_vv"]("b,c")  = eri["bbbb_oovv"]("k,l,c,d") * t2["bbbb_vvoo"]("d,b,l,k");

    // flops: o0v2  = o2v3
    //  mems: o0v2  = o0v2
    tmps_["7_bb_vv"]("b,c")  = eri["abab_oovv"]("l,k,d,c") * t2["abab_vvoo"]("d,b,l,k");

    // flops: o0v2  = o2v3
    //  mems: o0v2  = o0v2
    tmps_["6_aa_vv"]("a,d")  = eri["aaaa_oovv"]("k,l,c,d") * t2["aaaa_vvoo"]("c,a,l,k");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["4_bbbb_vvoo"]("b,c,j,k")  = eri["bbbb_oovv"]("k,l,c,d") * t2["bbbb_vvoo"]("d,b,j,l");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["3_bbbb_vvoo"]("b,c,j,k")  = eri["abab_oovv"]("l,k,d,c") * t2["abab_vvoo"]("d,b,l,j");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["2_abab_vvoo"]("c,b,k,j")  = eri["aaaa_oovv"]("k,l,c,d") * t2["abab_vvoo"]("d,b,l,j");

    // flops: o2v2  = o3v3
    //  mems: o2v2  = o2v2
    tmps_["1_abab_vvoo"]("b,c,i,k")  = eri["abab_oovv"]("l,k,d,c") * t2["aaaa_vvoo"]("d,b,i,l");

    // flops: o2v2  = o3v2 o3v3 o2v3 o2v2 o2v3 o2v2 o3v3 o3v3 o2v2 o2v2 o3v3 o2v3 o2v2 o2v4 o2v3 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o3v2 o2v2 o2v2 o2v3 o2v2 o4v2 o2v2 o3v3 o2v2 o3v2 o2v2 o2v2 o2v2 o2v2 o2v2 o3v2 o2v2
    //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
    tmps_["20_abab_vvoo"]("a,b,i,j")  = t2["abab_vvoo"]("a,b,i,l") * tmps_["10_bb_oo"]("l,j");
    tmps_["20_abab_vvoo"]("a,b,i,j") += 2.00 * t2["abab_vvoo"]("a,c,i,k") * tmps_["3_bbbb_vvoo"]("b,c,j,k");
    tmps_["20_abab_vvoo"]("a,b,i,j") += t2["abab_vvoo"]("d,b,i,j") * tmps_["6_aa_vv"]("a,d");
    tmps_["20_abab_vvoo"]("a,b,i,j") -= t2["abab_vvoo"]("a,c,i,j") * tmps_["8_bb_vv"]("b,c");
    tmps_["20_abab_vvoo"]("a,b,i,j") -= 2.00 * t2["abab_vvoo"]("a,c,i,k") * tmps_["4_bbbb_vvoo"]("b,c,j,k");
    tmps_["20_abab_vvoo"]("a,b,i,j") += 2.00 * t2["bbbb_vvoo"]("d,b,j,l") * tmps_["1_abab_vvoo"]("a,d,i,l");
    tmps_["20_abab_vvoo"]("a,b,i,j") -= 2.00 * t2["aaaa_vvoo"]("c,a,i,k") * tmps_["2_abab_vvoo"]("c,b,k,j");
    tmps_["20_abab_vvoo"]("a,b,i,j") -= 2.00 * t2["abab_vvoo"]("a,c,i,j") * tmps_["7_bb_vv"]("b,c");
    tmps_["20_abab_vvoo"]("a,b,i,j") += 2.00 * eri["abab_vvvv"]("a,b,c,d") * t2["abab_vvoo"]("c,d,i,j");
    tmps_["20_abab_vvoo"]("a,b,i,j") += 2.00 * f["bb_vv"]("b,c") * t2["abab_vvoo"]("a,c,i,j");
    tmps_["20_abab_vvoo"]("a,b,i,j") += 2.00 * eri["baab_vovo"]("b,k,c,j") * t2["aaaa_vvoo"]("c,a,i,k");
    tmps_["20_abab_vvoo"]("a,b,i,j") -= 2.00 * eri["baba_vovo"]("b,k,c,i") * t2["abab_vvoo"]("a,c,k,j");
    tmps_["20_abab_vvoo"]("a,b,i,j") -= 2.00 * eri["aaaa_vovo"]("a,k,c,i") * t2["abab_vvoo"]("c,b,k,j");
    tmps_["20_abab_vvoo"]("a,b,i,j") -= 2.00 * eri["bbbb_vovo"]("b,k,c,j") * t2["abab_vvoo"]("a,c,i,k");
    tmps_["20_abab_vvoo"]("a,b,i,j") -= 2.00 * eri["abab_vovo"]("a,k,c,j") * t2["abab_vvoo"]("c,b,i,k");
    tmps_["20_abab_vvoo"]("a,b,i,j") -= 2.00 * f["bb_oo"]("k,j") * t2["abab_vvoo"]("a,b,i,k");
    tmps_["20_abab_vvoo"]("a,b,i,j") += 2.00 * eri["abab_vvoo"]("a,b,i,j");
    tmps_["20_abab_vvoo"]("a,b,i,j") += 2.00 * f["aa_vv"]("a,c") * t2["abab_vvoo"]("c,b,i,j");
    tmps_["20_abab_vvoo"]("a,b,i,j") += 2.00 * eri["abab_oooo"]("l,k,i,j") * t2["abab_vvoo"]("a,b,l,k");
    tmps_["20_abab_vvoo"]("a,b,i,j") += 2.00 * eri["abba_vovo"]("a,k,c,i") * t2["bbbb_vvoo"]("c,b,j,k");
    tmps_["20_abab_vvoo"]("a,b,i,j") -= 2.00 * f["aa_oo"]("k,i") * t2["abab_vvoo"]("a,b,k,j");
    tmps_["20_abab_vvoo"]("a,b,i,j") -= 2.00 * t2["abab_vvoo"]("a,b,l,j") * tmps_["9_aa_oo"]("l,i");

    // rt2_abab += -0.50 <l,k||c,d>_bbbb t2_abab(a,b,i,l) t2_bbbb(c,d,j,k)
    //          += -0.50 <k,l||c,d>_abab t2_abab(a,b,i,l) t2_abab(c,d,k,j)
    //          += -0.50 <k,l||d,c>_abab t2_abab(a,b,i,l) t2_abab(d,c,k,j)
    //          += +1.00 <l,k||d,c>_abab t2_abab(a,c,i,k) t2_abab(d,b,l,j)
    //          += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,a,l,k) t2_abab(d,b,i,j)
    //          += +0.50 <l,k||c,d>_bbbb t2_abab(a,c,i,j) t2_bbbb(d,b,l,k)
    //          += +1.00 <l,k||c,d>_bbbb t2_abab(a,c,i,k) t2_bbbb(d,b,j,l)
    //          += +1.00 <k,l||c,d>_abab t2_aaaa(c,a,i,k) t2_bbbb(d,b,j,l)
    //          += +1.00 <l,k||c,d>_aaaa t2_aaaa(c,a,i,k) t2_abab(d,b,l,j)
    //          += -0.50 <l,k||d,c>_abab t2_abab(a,c,i,j) t2_abab(d,b,l,k)
    //          += -0.50 <k,l||d,c>_abab t2_abab(a,c,i,j) t2_abab(d,b,k,l)
    //          += +0.50 <a,b||c,d>_abab t2_abab(c,d,i,j)
    //          += +0.50 <a,b||d,c>_abab t2_abab(d,c,i,j)
    //          += +1.00 f_bb(b,c) t2_abab(a,c,i,j)
    //          += -1.00 <k,b||c,j>_abab t2_aaaa(c,a,i,k)
    //          += -1.00 <k,b||i,c>_abab t2_abab(a,c,k,j)
    //          += +1.00 <k,a||c,i>_aaaa t2_abab(c,b,k,j)
    //          += +1.00 <k,b||c,j>_bbbb t2_abab(a,c,i,k)
    //          += -1.00 <a,k||c,j>_abab t2_abab(c,b,i,k)
    //          += -1.00 f_bb(k,j) t2_abab(a,b,i,k)
    //          += +1.00 <a,b||i,j>_abab
    //          += +1.00 f_aa(a,c) t2_abab(c,b,i,j)
    //          += +0.50 <l,k||i,j>_abab t2_abab(a,b,l,k)
    //          += +0.50 <k,l||i,j>_abab t2_abab(a,b,k,l)
    //          += -1.00 <a,k||i,c>_abab t2_bbbb(c,b,j,k)
    //          += -1.00 f_aa(k,i) t2_abab(a,b,k,j)
    //          += -0.50 <l,k||c,d>_abab t2_abab(a,b,l,j) t2_abab(c,d,i,k)
    //          += -0.50 <l,k||d,c>_abab t2_abab(a,b,l,j) t2_abab(d,c,i,k)
    //          += -0.50 <l,k||c,d>_aaaa t2_abab(a,b,l,j) t2_aaaa(c,d,i,k)
    rt2_abab("a,b,i,j") += 0.50 * tmps_["20_abab_vvoo"]("a,b,i,j");
    tmps_["20_abab_vvoo"].~TArrayD();

}
