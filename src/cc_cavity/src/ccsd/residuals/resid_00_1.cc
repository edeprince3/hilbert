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

void CCSD::resid_00_1() {

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

    scalars_["1"]  = 1.00 * dot(eri["aaaa_oovv"]("i,j,a,b"), t2["aaaa_vvoo"]("a,b,j,i"));
    scalars_["2"]  = 1.00 * dot(eri["aaaa_oovv"]("i,j,a,b"), t2["aaaa_vvoo"]("a,b,j,i"));
    scalars_["2"] += 2.00 * dot(Id["bb_oo"]("i,i1") * eri["bbbb_oooo"]("i,j,i1,j1"), Id["bb_oo"]("j,j1"));
    scalars_["2"] += 4.00 * dot(Id["bb_oo"]("i,i1") * eri["abab_oooo"]("j,i,j1,i1"), Id["aa_oo"]("j,j1"));
    scalars_["2"] += 2.00 * dot(Id["aa_oo"]("i,i1") * eri["aaaa_oooo"]("i,j,i1,j1"), Id["aa_oo"]("j,j1"));

    // energy  = +0.25 <j,i||a,b>_aaaa t2_aaaa(a,b,j,i)
    //        += -0.50 <j,i||j,i>_bbbb
    //        += -0.50 <j,i||j,i>_abab
    //        += -0.50 <i,j||i,j>_abab
    //        += -0.50 <j,i||j,i>_aaaa
    energy  = -0.25 * scalars_["2"];

    // energy += +1.00 f_aa(i,i)
    // flops: 0 += o2v0
    //  mems: 0 += 0
    energy += 1.00 * dot(Id["aa_oo"]("i,i1"), f["aa_oo"]("i,i1"));

    // energy += +1.00 f_bb(i,i)
    // flops: 0 += o2v0
    //  mems: 0 += 0
    energy += 1.00 * dot(Id["bb_oo"]("i,i1"), f["bb_oo"]("i,i1"));

    // energy += +0.25 <j,i||a,b>_bbbb t2_bbbb(a,b,j,i)
    // flops: 0 += o2v2
    //  mems: 0 += 0
    energy -= 0.25 * dot(eri["bbbb_oovv"]("i,j,a,b"), t2["bbbb_vvoo"]("a,b,j,i"));

    // energy += +0.25 <j,i||a,b>_abab t2_abab(a,b,j,i)
    //        += +0.25 <i,j||a,b>_abab t2_abab(a,b,i,j)
    //        += +0.25 <j,i||b,a>_abab t2_abab(b,a,j,i)
    //        += +0.25 <i,j||b,a>_abab t2_abab(b,a,i,j)
    // flops: 0 += o2v2
    //  mems: 0 += 0
    energy += 1.00 * dot(eri["abab_oovv"]("j,i,a,b"), t2["abab_vvoo"]("a,b,j,i"));

    // rt2_aaaa  = +1.00 P(a,b) f_aa(a,c) t2_aaaa(c,b,i,j)
    // flops: o2v2  = o2v3
    //  mems: o2v2  = o2v2
    tmps_["perm_aaaa_vvoo"]("a,b,i,j")  = f["aa_vv"]("a,c") * t2["aaaa_vvoo"]("c,b,i,j");
    rt2_aaaa("a,b,i,j")  = tmps_["perm_aaaa_vvoo"]("a,b,i,j");
    rt2_aaaa("a,b,i,j") -= tmps_["perm_aaaa_vvoo"]("b,a,i,j");

    // rt2_aaaa += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,a,i,j) t2_aaaa(d,b,l,k)
    // flops: o2v2 += o2v3 o2v3
    //  mems: o2v2 += o0v2 o2v2
    rt2_aaaa("a,b,i,j") += 0.50 * eri["aaaa_oovv"]("k,l,c,d") * t2["aaaa_vvoo"]("d,b,l,k") * t2["aaaa_vvoo"]("c,a,i,j");

    // rt2_aaaa += +0.25 <l,k||c,d>_aaaa t2_aaaa(a,b,l,k) t2_aaaa(c,d,i,j)
    // flops: o2v2 += o4v2 o4v2
    //  mems: o2v2 += o4v0 o2v2
    rt2_aaaa("a,b,i,j") -= 0.25 * eri["aaaa_oovv"]("k,l,c,d") * t2["aaaa_vvoo"]("c,d,i,j") * t2["aaaa_vvoo"]("a,b,l,k");

    // rt2_aaaa += +1.00 P(i,j) <l,k||c,d>_bbbb t2_abab(a,c,j,k) t2_abab(b,d,i,l)
    // flops: o2v2 += o3v3 o3v3
    //  mems: o2v2 += o2v2 o2v2
    tmps_["perm_aaaa_vvoo"]("a,b,i,j")  = eri["bbbb_oovv"]("k,l,c,d") * t2["abab_vvoo"]("a,c,j,k") * t2["abab_vvoo"]("b,d,i,l");
    rt2_aaaa("a,b,i,j") -= tmps_["perm_aaaa_vvoo"]("a,b,i,j");
    rt2_aaaa("a,b,i,j") += tmps_["perm_aaaa_vvoo"]("a,b,j,i");

    // rt2_abab  = -0.50 <l,k||d,c>_abab t2_abab(a,c,l,k) t2_abab(d,b,i,j)
    //          += -0.50 <k,l||d,c>_abab t2_abab(a,c,k,l) t2_abab(d,b,i,j)
    // flops: o2v2  = o2v3 o2v3
    //  mems: o2v2  = o0v2 o2v2
    rt2_abab("a,b,i,j")  = -1.00 * eri["abab_oovv"]("l,k,d,c") * t2["abab_vvoo"]("a,c,l,k") * t2["abab_vvoo"]("d,b,i,j");

    // rt2_abab += +0.25 <l,k||c,d>_abab t2_abab(a,b,l,k) t2_abab(c,d,i,j)
    //          += +0.25 <l,k||d,c>_abab t2_abab(a,b,l,k) t2_abab(d,c,i,j)
    //          += +0.25 <k,l||c,d>_abab t2_abab(a,b,k,l) t2_abab(c,d,i,j)
    //          += +0.25 <k,l||d,c>_abab t2_abab(a,b,k,l) t2_abab(d,c,i,j)
    // flops: o2v2 += o4v2 o4v2
    //  mems: o2v2 += o4v0 o2v2
    rt2_abab("a,b,i,j") += eri["abab_oovv"]("l,k,c,d") * t2["abab_vvoo"]("c,d,i,j") * t2["abab_vvoo"]("a,b,l,k");

    // rt2_abab += +1.00 <k,l||d,c>_abab t2_abab(a,c,k,j) t2_abab(d,b,i,l)
    // flops: o2v2 += o3v3 o3v3
    //  mems: o2v2 += o2v2 o2v2
    rt2_abab("a,b,i,j") += eri["abab_oovv"]("k,l,d,c") * t2["abab_vvoo"]("d,b,i,l") * t2["abab_vvoo"]("a,c,k,j");

    // rt2_bbbb  = -1.00 P(i,j) f_bb(k,j) t2_bbbb(a,b,i,k)
    // flops: o2v2  = o3v2
    //  mems: o2v2  = o2v2
    tmps_["perm_bbbb_vvoo"]("a,b,i,j")  = f["bb_oo"]("k,j") * t2["bbbb_vvoo"]("a,b,i,k");
    rt2_bbbb("a,b,i,j")  = -1.00 * tmps_["perm_bbbb_vvoo"]("a,b,i,j");
    rt2_bbbb("a,b,i,j") += tmps_["perm_bbbb_vvoo"]("a,b,j,i");

    // rt2_bbbb += +1.00 P(a,b) f_bb(a,c) t2_bbbb(c,b,i,j)
    // flops: o2v2 += o2v3
    //  mems: o2v2 += o2v2
    tmps_["perm_bbbb_vvoo"]("a,b,i,j")  = f["bb_vv"]("a,c") * t2["bbbb_vvoo"]("c,b,i,j");
    rt2_bbbb("a,b,i,j") += tmps_["perm_bbbb_vvoo"]("a,b,i,j");
    rt2_bbbb("a,b,i,j") -= tmps_["perm_bbbb_vvoo"]("b,a,i,j");

    // rt2_bbbb += -0.50 <l,k||c,d>_bbbb t2_bbbb(c,a,l,k) t2_bbbb(d,b,i,j)
    // flops: o2v2 += o2v3 o2v3
    //  mems: o2v2 += o0v2 o2v2
    rt2_bbbb("a,b,i,j") += 0.50 * eri["bbbb_oovv"]("k,l,c,d") * t2["bbbb_vvoo"]("c,a,l,k") * t2["bbbb_vvoo"]("d,b,i,j");

}
