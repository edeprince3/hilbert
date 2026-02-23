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

void QED_CCSD_22::resid_22_6() {

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


        // rt2_2_bbbb += -2.00 P(i,j) P(a,b) d-_bb(k,c) t1_1_bb(a,k) t1_1_bb(b,i) t1_1_bb(c,j)
        //            += -2.00 P(i,j) P(a,b) <l,k||c,j>_bbbb t1_1_bb(a,k) t2_1_bbbb(c,b,i,l)
        //            += +2.00 P(i,j) P(a,b) <l,k||c,j>_abab t1_1_bb(a,k) t2_1_abab(c,b,l,i)
        //            += -2.00 P(i,j) P(a,b) d-_bb(k,c) t1_1_bb(a,j) t2_2_bbbb(c,b,i,k)
        //            += +2.00 P(i,j) P(a,b) d-_aa(k,c) t1_1_bb(a,j) t2_2_abab(c,b,k,i)
        //            += +2.00 P(i,j) P(a,b) <k,a||c,j>_bbbb t2_2_bbbb(c,b,i,k)
        //            += -2.00 P(i,j) P(a,b) <k,a||c,j>_abab t2_2_abab(c,b,k,i)
        //            += -4.00 P(i,j) P(a,b) d-_bb(k,j) t1_1_bb(a,k) t1_2_bb(b,i)
        //            += +4.00 P(i,j) P(a,b) d-_bb(a,c) t1_1_bb(c,j) t1_2_bb(b,i)
        //            += -2.00 P(i,j) P(a,b) <l,k||c,j>_bbbb t2_bbbb(c,a,i,k) t1_2_bb(b,l)
        //            += -2.00 P(i,j) P(a,b) <k,l||c,j>_abab t2_abab(c,a,k,i) t1_2_bb(b,l)
        //            += -4.00 P(i,j) P(a,b) d-_bb(k,c) t2_1_bbbb(c,a,j,k) t1_2_bb(b,i)
        //            += +4.00 P(i,j) P(a,b) d-_aa(k,c) t2_1_abab(c,a,k,j) t1_2_bb(b,i)
        //            += +2.00 P(i,j) P(a,b) <k,l||c,d>_abab t2_abab(c,a,k,j) t1_1_bb(b,l) t1_1_bb(d,i)
        //            += +2.00 P(i,j) P(a,b) <l,k||c,d>_bbbb t2_bbbb(c,a,j,k) t1_1_bb(b,l) t1_1_bb(d,i)
        //            += +2.00 P(i,j) P(a,b) <l,k||c,d>_aaaa t2_abab(c,a,k,j) t2_2_abab(d,b,l,i)
        //            += +2.00 P(i,j) P(a,b) <l,k||d,c>_abab t2_bbbb(c,a,j,k) t2_2_abab(d,b,l,i)
        //            += -2.00 P(i,j) P(a,b) d-_bb(a,c) t1_1_bb(b,j) t1_2_bb(c,i)
        //            += -2.00 P(i,j) P(a,b) <k,a||d,c>_abab t2_1_abab(d,b,k,i) t1_1_bb(c,j)
        //            += +2.00 P(i,j) P(a,b) <k,a||c,d>_abab t2_abab(c,b,k,j) t1_2_bb(d,i)
        //            += +2.00 P(i,j) P(a,b) <k,l||c,d>_abab t2_abab(c,a,k,j) t2_2_bbbb(d,b,i,l)
        //            += -2.00 P(i,j) P(a,b) <k,a||c,d>_bbbb t2_1_bbbb(d,b,i,k) t1_1_bb(c,j)
        //            += +2.00 P(i,j) P(a,b) <k,a||c,j>_bbbb t1_1_bb(b,k) t1_1_bb(c,i)
        //            += +2.00 P(i,j) P(a,b) d-_bb(k,j) t1_1_bb(a,i) t1_2_bb(b,k)
        //            += -2.00 P(i,j) P(a,b) <k,a||c,d>_bbbb t2_bbbb(c,b,j,k) t1_2_bb(d,i)
        //            += +2.00 P(i,j) P(a,b) <l,k||c,d>_bbbb t2_bbbb(c,a,j,k) t2_2_bbbb(d,b,i,l)
        rt2_2_bbbb("a,b,i,j") -= 2.00 * tmps_["202_bbbb_vovo"]("b,i,a,j");
        rt2_2_bbbb("a,b,i,j") += 2.00 * tmps_["202_bbbb_vovo"]("b,j,a,i");
        rt2_2_bbbb("a,b,i,j") += 2.00 * tmps_["202_bbbb_vovo"]("a,i,b,j");
        rt2_2_bbbb("a,b,i,j") -= 2.00 * tmps_["202_bbbb_vovo"]("a,j,b,i");
        tmps_["202_bbbb_vovo"].~TArrayD();

        // flops: o2v2  = o2v3 o3v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2
        tmps_["203_bbbb_vvoo"]("a,b,i,j")  = dp["bb_vv"]("a,c") * t2_1["bbbb"]("c,b,i,j");
        tmps_["203_bbbb_vvoo"]("a,b,i,j") += tmps_["113_bbbb_ovoo"]("k,a,i,j") * t1_1["bb"]("b,k");

        // rt2_bbbb  = -1.00 P(a,b) d-_bb(k,c) t2_bbbb(c,a,i,j) t1_1_bb(b,k)
        //          += -1.00 P(a,b) d-_bb(a,c) t2_1_bbbb(c,b,i,j)
        rt2_bbbb("a,b,i,j")  = -1.00 * tmps_["203_bbbb_vvoo"]("a,b,i,j");
        rt2_bbbb("a,b,i,j") += tmps_["203_bbbb_vvoo"]("b,a,i,j");

        // rt2_2_bbbb += -2.00 P(a,b) d+_bb(k,c) t2_bbbb(c,a,i,j) t1_1_bb(b,k)
        //            += -2.00 P(a,b) d+_bb(a,c) t2_1_bbbb(c,b,i,j)
        rt2_2_bbbb("a,b,i,j") -= 2.00 * tmps_["203_bbbb_vvoo"]("a,b,i,j");
        rt2_2_bbbb("a,b,i,j") += 2.00 * tmps_["203_bbbb_vvoo"]("b,a,i,j");
        tmps_["203_bbbb_vvoo"].~TArrayD();

        // flops: o2v2  = o2v3
        //  mems: o2v2  = o2v2
        tmps_["204_bbbb_voov"]("a,i,j,b")  = t2_1["bbbb"]("c,a,i,j") * tmps_["50_bb_vv"]("c,b");
        tmps_["50_bb_vv"].~TArrayD();

        // rt2_2_bbbb += +1.00 <l,k||d,c>_abab t2_1_bbbb(c,a,i,j) t2_1_abab(d,b,l,k)
        //            += +1.00 <k,l||d,c>_abab t2_1_bbbb(c,a,i,j) t2_1_abab(d,b,k,l)
        rt2_2_bbbb("a,b,i,j") += 2.00 * tmps_["204_bbbb_voov"]("a,i,j,b");

        // rt2_2_bbbb += -1.00 <l,k||c,d>_abab t2_1_abab(c,a,l,k) t2_1_bbbb(d,b,i,j)
        //            += -1.00 <k,l||c,d>_abab t2_1_abab(c,a,k,l) t2_1_bbbb(d,b,i,j)
        rt2_2_bbbb("a,b,i,j") -= 2.00 * tmps_["204_bbbb_voov"]("b,i,j,a");
        tmps_["204_bbbb_voov"].~TArrayD();

        // flops: o4v0  = o4v1
        //  mems: o4v0  = o4v0
        tmps_["205_bbbb_oooo"]("i,j,k,l")  = eri["bbbb_oovo"]("i,j,a,k") * t1_2["bb"]("a,l");

        // flops: o2v2  = o3v2 o4v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v3 o2v2 o4v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v3 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o2v2 o3v2 o3v2 o2v2 o2v2 o2v2 o3v2 o2v2 o3v3 o2v2 o3v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["206_bbbb_vvoo"]("a,b,i,j")  = 3.00 * tmps_["24_bb_oo"]("k,j") * t2_2["bbbb"]("a,b,i,k");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") -= 0.50 * tmps_["108_bbbb_oooo"]("k,n,j,i") * t2_1["bbbb"]("a,b,n,k");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") -= tmps_["54_bb_oo"]("n,j") * t2_1["bbbb"]("a,b,i,n");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") -= 3.00 * tmps_["58_bb_oo"]("k,i") * t2_1["bbbb"]("a,b,j,k");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") -= tmps_["59_bb_oo"]("k,j") * t2_1["bbbb"]("a,b,i,k");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") += 0.50 * tmps_["23_bb_oo"]("n,j") * t2_2["bbbb"]("a,b,i,n");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") -= tmps_["49_bbbb_ovvo"]("k,f,b,i") * t2_1["bbbb"]("f,a,j,k");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") -= 0.50 * tmps_["205_bbbb_oooo"]("k,n,j,i") * t2["bbbb"]("a,b,n,k");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") -= tmps_["57_bb_oo"]("n,j") * t2_1["bbbb"]("a,b,i,n");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") -= tmps_["22_bb_oo"]("n,j") * t2_2["bbbb"]("a,b,i,n");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") -= tmps_["56_bb_oo"]("n,j") * t2_1["bbbb"]("a,b,i,n");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") -= f["bb_oo"]("k,j") * t2_2["bbbb"]("a,b,i,k");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") += eri["bbbb_vvvo"]("a,b,f,j") * t1_2["bb"]("f,i");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") -= t2["bbbb"]("a,b,i,k") * tmps_["175_bb_oo"]("k,j");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") += t2_2["bbbb"]("a,b,i,k") * tmps_["61_bb_oo"]("k,j");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") += 2.00 * t2_1["bbbb"]("a,b,i,k") * tmps_["62_bb_oo"]("k,j");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") += t2["bbbb"]("a,b,j,k") * tmps_["176_bb_oo"]("k,i");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") -= t2["bbbb"]("a,b,j,k") * tmps_["182_bb_oo"]("k,i");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") += t2_1["bbbb"]("a,b,i,k") * tmps_["106_bb_oo"]("k,j");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") += 0.50 * tmps_["166_bb_oo"]("n,j") * t2_1["bbbb"]("a,b,i,n");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") -= tmps_["47_aabb_ovvo"]("m,d,b,i") * t2_1["abab"]("d,a,m,j");
        tmps_["206_bbbb_vvoo"]("a,b,i,j") += tmps_["186_bb_oo"]("i,k") * t2["bbbb"]("a,b,j,k");
        tmps_["205_bbbb_oooo"].~TArrayD();
        tmps_["186_bb_oo"].~TArrayD();
        tmps_["182_bb_oo"].~TArrayD();
        tmps_["176_bb_oo"].~TArrayD();
        tmps_["175_bb_oo"].~TArrayD();
        tmps_["166_bb_oo"].~TArrayD();
        tmps_["108_bbbb_oooo"].~TArrayD();
        tmps_["106_bb_oo"].~TArrayD();
        tmps_["62_bb_oo"].~TArrayD();
        tmps_["59_bb_oo"].~TArrayD();
        tmps_["58_bb_oo"].~TArrayD();
        tmps_["57_bb_oo"].~TArrayD();
        tmps_["56_bb_oo"].~TArrayD();
        tmps_["54_bb_oo"].~TArrayD();
        tmps_["49_bbbb_ovvo"].~TArrayD();
        tmps_["47_aabb_ovvo"].~TArrayD();

        // rt2_2_bbbb += +2.00 P(i,j) <l,k||c,d>_abab t2_bbbb(a,b,j,k) t1_1_aa(c,l) t1_1_bb(d,i)
        //            += +1.00 P(i,j) <l,k||c,j>_bbbb t2_1_bbbb(a,b,l,k) t1_1_bb(c,i)
        //            += +6.00 P(i,j) d-_bb(k,c) t1_1_bb(c,j) t2_2_bbbb(a,b,i,k)
        //            += -1.00 P(i,j) <k,l||c,d>_abab t2_1_bbbb(a,b,i,l) t2_1_abab(c,d,k,j)
        //            += -1.00 P(i,j) <k,l||d,c>_abab t2_1_bbbb(a,b,i,l) t2_1_abab(d,c,k,j)
        //            += -6.00 P(i,j) d-_bb(k,c) t2_1_bbbb(a,b,j,k) t1_2_bb(c,i)
        //            += -2.00 P(i,j) f_bb(k,c) t2_1_bbbb(a,b,i,k) t1_1_bb(c,j)
        //            += +2.00 P(i,j) <l,k||c,d>_bbbb t2_1_bbbb(c,a,j,k) t2_1_bbbb(d,b,i,l)
        //            += -1.00 P(i,j) <l,k||c,d>_bbbb t2_bbbb(c,d,j,k) t2_2_bbbb(a,b,i,l)
        //            += +1.00 P(i,j) <l,k||c,j>_bbbb t2_bbbb(a,b,l,k) t1_2_bb(c,i)
        //            += +2.00 P(i,j) <l,k||c,j>_bbbb t2_1_bbbb(a,b,i,l) t1_1_bb(c,k)
        //            += -1.00 P(i,j) <k,l||c,d>_abab t2_abab(c,d,k,j) t2_2_bbbb(a,b,i,l)
        //            += -1.00 P(i,j) <k,l||d,c>_abab t2_abab(d,c,k,j) t2_2_bbbb(a,b,i,l)
        //            += -2.00 P(i,j) <k,l||c,j>_abab t2_1_bbbb(a,b,i,l) t1_1_aa(c,k)
        //            += +2.00 P(i,j) <l,k||c,d>_aaaa t2_1_abab(c,a,k,j) t2_1_abab(d,b,l,i)
        //            += +2.00 P(i,j) d-_bb(k,c) t0_1 t2_1_bbbb(a,b,i,k) t1_1_bb(c,j)
        //            += +2.00 P(i,j) <l,k||c,d>_bbbb t2_bbbb(a,b,j,k) t1_1_bb(c,l) t1_1_bb(d,i)
        //            += -2.00 P(i,j) d-_bb(k,c) t2_bbbb(a,b,j,k) t0_1 t1_2_bb(c,i)
        //            += -4.00 P(i,j) d-_bb(k,c) t2_bbbb(a,b,j,k) t1_1_bb(c,i) t0_2
        //            += -2.00 P(i,j) <l,k||c,j>_abab t2_bbbb(a,b,i,k) t1_2_aa(c,l)
        //            += -2.00 P(i,j) <l,k||c,j>_bbbb t2_bbbb(a,b,i,k) t1_2_bb(c,l)
        //            += +2.00 P(i,j) <a,b||c,j>_bbbb t1_2_bb(c,i)
        //            += -2.00 P(i,j) f_bb(k,j) t2_2_bbbb(a,b,i,k)
        //            += +2.00 P(i,j) d-_bb(k,j) t0_1 t2_2_bbbb(a,b,i,k)
        //            += +4.00 P(i,j) d-_bb(k,j) t2_1_bbbb(a,b,i,k) t0_2
        //            += +2.00 P(i,j) f_bb(k,c) t2_bbbb(a,b,j,k) t1_2_bb(c,i)
        //            += -1.00 P(i,j) <l,k||c,d>_bbbb t2_bbbb(a,b,j,k) t2_2_bbbb(c,d,i,l)
        //            += +1.00 P(i,j) <l,k||c,d>_abab t2_bbbb(a,b,j,k) t2_2_abab(c,d,l,i)
        //            += +1.00 P(i,j) <l,k||d,c>_abab t2_bbbb(a,b,j,k) t2_2_abab(d,c,l,i)
        //            += -1.00 P(i,j) <l,k||c,d>_bbbb t2_1_bbbb(a,b,i,l) t2_1_bbbb(c,d,j,k)
        rt2_2_bbbb("a,b,i,j") += 2.00 * tmps_["206_bbbb_vvoo"]("a,b,i,j");
        rt2_2_bbbb("a,b,i,j") -= 2.00 * tmps_["206_bbbb_vvoo"]("a,b,j,i");
        tmps_["206_bbbb_vvoo"].~TArrayD();

        // flops: o2v2  = o3v3
        //  mems: o2v2  = o2v2
        tmps_["207_bbbb_vovo"]("a,i,b,j")  = t2_1["bbbb"]("c,a,i,k") * tmps_["48_bbbb_ovvo"]("k,c,b,j");
        tmps_["48_bbbb_ovvo"].~TArrayD();

        // rt2_2_bbbb += +2.00 P(i,j) <l,k||d,c>_abab t2_1_bbbb(c,a,j,k) t2_1_abab(d,b,l,i)
        rt2_2_bbbb("a,b,i,j") += 2.00 * tmps_["207_bbbb_vovo"]("a,j,b,i");
        rt2_2_bbbb("a,b,i,j") -= 2.00 * tmps_["207_bbbb_vovo"]("a,i,b,j");

        // rt2_2_bbbb += +2.00 P(i,j) <k,l||c,d>_abab t2_1_abab(c,a,k,j) t2_1_bbbb(d,b,i,l)
        rt2_2_bbbb("a,b,i,j") += 2.00 * tmps_["207_bbbb_vovo"]("b,i,a,j");
        rt2_2_bbbb("a,b,i,j") -= 2.00 * tmps_["207_bbbb_vovo"]("b,j,a,i");
        tmps_["207_bbbb_vovo"].~TArrayD();

        // flops: o3v1  = o1v1 o3v2
        //  mems: o3v1  = o1v1 o3v1
        tmps_["210_bbbb_vooo"]("a,i,j,k")  = (-1.00 * tmps_["169_bb_ov"]("k,b") + tmps_["8_bb_ov"]("k,b")) * t2["bbbb"]("b,a,i,j");
        tmps_["169_bb_ov"].~TArrayD();
        tmps_["8_bb_ov"].~TArrayD();

        // flops: o3v1  = o3v2 o3v2 o3v1
        //  mems: o3v1  = o3v1 o3v1 o3v1
        tmps_["209_bbbb_ovoo"]("i,a,j,k")  = -3.00 * dp["bb_ov"]("i,b") * t2_2["bbbb"]("b,a,j,k");
        tmps_["209_bbbb_ovoo"]("i,a,j,k") += f["bb_ov"]("i,b") * t2_1["bbbb"]("b,a,j,k");

        // flops: o3v1  = o3v3
        //  mems: o3v1  = o3v1
        tmps_["208_bbbb_vooo"]("a,i,j,k")  = eri["bbbb_vovv"]("a,i,b,c") * t2_1["bbbb"]("b,c,j,k");

        // flops: o2v2  = o3v2 o3v2 o2v2 o3v2 o2v3 o2v2 o3v2 o2v2 o3v2 o2v2 o2v3 o2v2 o2v3 o2v2 o3v2 o2v2 o2v3 o2v2 o2v3 o2v2 o2v2 o2v3 o2v2 o2v3 o2v2 o2v3 o2v2 o3v2 o2v3 o2v2 o2v2 o3v2 o2v2 o3v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["211_bbbb_voov"]("a,i,j,b")  = -1.00 * t1_2["bb"]("b,l") * tmps_["112_bbbb_vooo"]("a,l,i,j");
        tmps_["211_bbbb_voov"]("a,i,j,b") -= t1_1["bb"]("b,l") * tmps_["208_bbbb_vooo"]("a,l,i,j");
        tmps_["211_bbbb_voov"]("a,i,j,b") -= 2.00 * eri["bbbb_vooo"]("a,l,i,j") * t1_2["bb"]("b,l");
        tmps_["211_bbbb_voov"]("a,i,j,b") += 2.00 * f["bb_vv"]("a,c") * t2_2["bbbb"]("c,b,i,j");
        tmps_["211_bbbb_voov"]("a,i,j,b") -= 4.00 * tmps_["113_bbbb_ovoo"]("l,a,i,j") * tmps_["177_bb_vo"]("b,l");
        tmps_["211_bbbb_voov"]("a,i,j,b") -= 2.00 * t1_1["bb"]("a,l") * tmps_["209_bbbb_ovoo"]("l,b,i,j");
        tmps_["211_bbbb_voov"]("a,i,j,b") -= 4.00 * t2_1["bbbb"]("c,b,i,j") * tmps_["60_bb_vv"]("a,c");
        tmps_["211_bbbb_voov"]("a,i,j,b") += 2.00 * t2["bbbb"]("c,a,i,j") * tmps_["172_bb_vv"]("c,b");
        tmps_["211_bbbb_voov"]("a,i,j,b") += 2.00 * tmps_["114_bbbb_ovoo"]("l,b,i,j") * tmps_["26_bb_vo"]("a,l");
        tmps_["211_bbbb_voov"]("a,i,j,b") -= 2.00 * t2_2["bbbb"]("c,b,i,j") * tmps_["25_bb_vv"]("a,c");
        tmps_["211_bbbb_voov"]("a,i,j,b") -= 2.00 * t2["bbbb"]("c,b,i,j") * tmps_["173_bb_vv"]("a,c");
        tmps_["211_bbbb_voov"]("a,i,j,b") -= 2.00 * t2_2["bbbb"]("d,b,i,j") * tmps_["20_bb_vv"]("d,a");
        tmps_["211_bbbb_voov"]("a,i,j,b") -= 2.00 * t2_1["bbbb"]("d,b,i,j") * tmps_["95_bb_vv"]("a,d");
        tmps_["211_bbbb_voov"]("a,i,j,b") -= 2.00 * t2_1["bbbb"]("d,b,i,j") * tmps_["167_bb_vv"]("a,d");
        tmps_["211_bbbb_voov"]("a,i,j,b") += 2.00 * t1_2["bb"]("b,l") * tmps_["115_bbbb_ovoo"]("l,a,i,j");
        tmps_["211_bbbb_voov"]("a,i,j,b") += tmps_["21_bb_vv"]("d,a") * t2_2["bbbb"]("d,b,i,j");
        tmps_["211_bbbb_voov"]("a,i,j,b") -= 6.00 * t1_2["bb"]("b,l") * tmps_["114_bbbb_ovoo"]("l,a,i,j");
        tmps_["211_bbbb_voov"]("a,i,j,b") += 2.00 * t1_1["bb"]("b,k") * tmps_["210_bbbb_vooo"]("a,i,j,k");
        tmps_["210_bbbb_vooo"].~TArrayD();
        tmps_["209_bbbb_ovoo"].~TArrayD();
        tmps_["208_bbbb_vooo"].~TArrayD();
        tmps_["177_bb_vo"].~TArrayD();
        tmps_["173_bb_vv"].~TArrayD();
        tmps_["172_bb_vv"].~TArrayD();
        tmps_["167_bb_vv"].~TArrayD();
        tmps_["115_bbbb_ovoo"].~TArrayD();
        tmps_["114_bbbb_ovoo"].~TArrayD();
        tmps_["113_bbbb_ovoo"].~TArrayD();
        tmps_["112_bbbb_vooo"].~TArrayD();
        tmps_["95_bb_vv"].~TArrayD();
        tmps_["60_bb_vv"].~TArrayD();
        tmps_["26_bb_vo"].~TArrayD();

        // rt2_2_bbbb += -1.00 P(a,b) <l,k||c,d>_bbbb t2_bbbb(c,a,l,k) t2_2_bbbb(d,b,i,j)
        //            += +2.00 P(a,b) f_bb(k,c) t2_bbbb(c,a,i,j) t1_2_bb(b,k)
        //            += -4.00 P(a,b) d-_bb(k,c) t2_bbbb(c,a,i,j) t1_1_bb(b,k) t0_2
        //            += -2.00 P(a,b) d-_bb(k,c) t2_bbbb(c,a,i,j) t0_1 t1_2_bb(b,k)
        //            += +2.00 P(a,b) f_bb(a,c) t2_2_bbbb(c,b,i,j)
        //            += +2.00 P(a,b) <k,a||i,j>_bbbb t1_2_bb(b,k)
        //            += -2.00 P(a,b) f_bb(k,c) t1_1_bb(a,k) t2_1_bbbb(c,b,i,j)
        //            += +6.00 P(a,b) d-_bb(k,c) t1_1_bb(a,k) t2_2_bbbb(c,b,i,j)
        //            += -4.00 P(a,b) d-_bb(a,c) t2_1_bbbb(c,b,i,j) t0_2
        //            += +1.00 P(a,b) <l,k||d,c>_abab t2_bbbb(c,a,i,j) t2_2_abab(d,b,l,k)
        //            += +1.00 P(a,b) <k,l||d,c>_abab t2_bbbb(c,a,i,j) t2_2_abab(d,b,k,l)
        //            += -1.00 P(a,b) <l,k||c,d>_bbbb t2_bbbb(c,a,i,j) t2_2_bbbb(d,b,l,k)
        //            += +2.00 P(a,b) d-_bb(k,c) t0_1 t1_1_bb(a,k) t2_1_bbbb(c,b,i,j)
        //            += -2.00 P(a,b) d-_bb(a,c) t0_1 t2_2_bbbb(c,b,i,j)
        //            += +2.00 P(a,b) <k,a||d,c>_abab t2_bbbb(c,b,i,j) t1_2_aa(d,k)
        //            += -2.00 P(a,b) <k,a||c,d>_bbbb t2_bbbb(c,b,i,j) t1_2_bb(d,k)
        //            += +1.00 P(a,b) <k,a||c,d>_bbbb t1_1_bb(b,k) t2_1_bbbb(c,d,i,j)
        //            += -1.00 P(a,b) <l,k||c,d>_abab t2_abab(c,a,l,k) t2_2_bbbb(d,b,i,j)
        //            += -1.00 P(a,b) <k,l||c,d>_abab t2_abab(c,a,k,l) t2_2_bbbb(d,b,i,j)
        //            += +2.00 P(a,b) <k,a||c,d>_abab t2_1_bbbb(d,b,i,j) t1_1_aa(c,k)
        //            += +1.00 P(a,b) <k,a||c,d>_bbbb t2_bbbb(c,d,i,j) t1_2_bb(b,k)
        //            += +2.00 P(a,b) <k,a||c,d>_bbbb t2_1_bbbb(d,b,i,j) t1_1_bb(c,k)
        //            += -6.00 P(a,b) d-_bb(k,c) t2_1_bbbb(c,a,i,j) t1_2_bb(b,k)
        //            += +2.00 P(a,b) <k,l||d,c>_abab t2_bbbb(c,a,i,j) t1_1_bb(b,l) t1_1_aa(d,k)
        //            += +2.00 P(a,b) <l,k||c,d>_bbbb t2_bbbb(c,a,i,j) t1_1_bb(b,l) t1_1_bb(d,k)
        rt2_2_bbbb("a,b,i,j") += tmps_["211_bbbb_voov"]("a,i,j,b");
        rt2_2_bbbb("a,b,i,j") -= tmps_["211_bbbb_voov"]("b,i,j,a");
        tmps_["211_bbbb_voov"].~TArrayD();

        // flops: o2v2  = o3v2 o3v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2
        tmps_["212_bbbb_ovvo"]("i,a,b,j")  = -1.00 * dp["bb_oo"]("k,i") * t2_1["bbbb"]("a,b,j,k");
        tmps_["212_bbbb_ovvo"]("i,a,b,j") += tmps_["24_bb_oo"]("k,j") * t2["bbbb"]("a,b,i,k");
        tmps_["24_bb_oo"].~TArrayD();

        // rt2_2_bbbb += -2.00 P(i,j) d+_bb(k,c) t2_bbbb(a,b,j,k) t1_1_bb(c,i)
        //            += +2.00 P(i,j) d+_bb(k,j) t2_1_bbbb(a,b,i,k)
        rt2_2_bbbb("a,b,i,j") -= 2.00 * tmps_["212_bbbb_ovvo"]("j,a,b,i");
        rt2_2_bbbb("a,b,i,j") += 2.00 * tmps_["212_bbbb_ovvo"]("i,a,b,j");

        // rt2_bbbb += -1.00 P(i,j) d-_bb(k,c) t2_bbbb(a,b,j,k) t1_1_bb(c,i)
        //          += +1.00 P(i,j) d-_bb(k,j) t2_1_bbbb(a,b,i,k)
        rt2_bbbb("a,b,i,j") -= tmps_["212_bbbb_ovvo"]("j,a,b,i");
        rt2_bbbb("a,b,i,j") += tmps_["212_bbbb_ovvo"]("i,a,b,j");
        tmps_["212_bbbb_ovvo"].~TArrayD();

        // flops: o2v2  = o3v2 o3v2 o2v2 o3v3 o2v2 o3v2 o2v2 o3v3 o2v2 o3v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["213_aaaa_ovvo"]("i,a,b,j")  = -0.50 * tmps_["6_aa_oo"]("k,i") * t2["aaaa"]("a,b,j,k");
        tmps_["213_aaaa_ovvo"]("i,a,b,j") -= tmps_["41_aa_oo"]("n,i") * t2["aaaa"]("a,b,j,n");
        tmps_["213_aaaa_ovvo"]("i,a,b,j") += t2["aaaa"]("f,b,j,k") * tmps_["1_aaaa_ovvo"]("k,f,a,i");
        tmps_["213_aaaa_ovvo"]("i,a,b,j") += f["aa_oo"]("n,i") * t2["aaaa"]("a,b,j,n");
        tmps_["213_aaaa_ovvo"]("i,a,b,j") += t2["abab"]("b,d,j,m") * tmps_["3_bbaa_ovvo"]("m,d,a,i");
        tmps_["213_aaaa_ovvo"]("i,a,b,j") += tmps_["7_aa_oo"]("k,i") * t2["aaaa"]("a,b,j,k");

        // rt2_aaaa += -0.50 P(i,j) <l,k||c,d>_abab t2_aaaa(a,b,i,l) t2_abab(c,d,j,k)
        //          += -0.50 P(i,j) <l,k||d,c>_abab t2_aaaa(a,b,i,l) t2_abab(d,c,j,k)
        //          += +1.00 P(i,j) <l,k||c,d>_aaaa t2_aaaa(c,a,j,k) t2_aaaa(d,b,i,l)
        //          += -1.00 P(i,j) f_aa(k,j) t2_aaaa(a,b,i,k)
        //          += +1.00 P(i,j) d-_aa(k,j) t2_aaaa(a,b,i,k) t0_1
        //          += +1.00 P(i,j) <l,k||c,d>_bbbb t2_abab(a,c,j,k) t2_abab(b,d,i,l)
        //          += -0.50 P(i,j) <l,k||c,d>_aaaa t2_aaaa(a,b,i,l) t2_aaaa(c,d,j,k)
        rt2_aaaa("a,b,i,j") -= tmps_["213_aaaa_ovvo"]("j,a,b,i");
        rt2_aaaa("a,b,i,j") += tmps_["213_aaaa_ovvo"]("i,a,b,j");
        tmps_["213_aaaa_ovvo"].~TArrayD();

        // flops: o2v2  = o3v3
        //  mems: o2v2  = o2v2
        tmps_["214_aaaa_vovo"]("a,i,b,j")  = t2["abab"]("a,c,i,k") * tmps_["2_bbaa_ovvo"]("k,c,b,j");

        // rt2_aaaa += +1.00 P(i,j) <k,l||c,d>_abab t2_aaaa(c,a,j,k) t2_abab(b,d,i,l)
        rt2_aaaa("a,b,i,j") += tmps_["214_aaaa_vovo"]("b,i,a,j");
        rt2_aaaa("a,b,i,j") -= tmps_["214_aaaa_vovo"]("b,j,a,i");

        // rt2_aaaa += +1.00 P(i,j) <l,k||d,c>_abab t2_abab(a,c,j,k) t2_aaaa(d,b,i,l)
        rt2_aaaa("a,b,i,j") += tmps_["214_aaaa_vovo"]("a,j,b,i");
        rt2_aaaa("a,b,i,j") -= tmps_["214_aaaa_vovo"]("a,i,b,j");
        tmps_["214_aaaa_vovo"].~TArrayD();

        // flops: o2v2  = o3v3 o2v2 o3v3 o2v2 o2v2 o2v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["215_aaaa_vovo"]("a,i,b,j")  = -1.00 * eri["aaaa_vovo"]("a,k,c,i") * t2["aaaa"]("c,b,j,k");
        tmps_["215_aaaa_vovo"]("a,i,b,j") += dp["aa_vo"]("a,i") * t1_1["aa"]("b,j");
        tmps_["215_aaaa_vovo"]("a,i,b,j") += eri["abba_vovo"]("a,l,d,i") * t2["abab"]("b,d,j,l");
        tmps_["215_aaaa_vovo"]("a,i,b,j") += t1_1["aa"]("b,j") * tmps_["75_aa_vo"]("a,i");

        // rt2_aaaa += +1.00 P(i,j) P(a,b) d-_bb(k,c) t2_abab(a,c,j,k) t1_1_aa(b,i)
        //          += -1.00 P(i,j) P(a,b) d-_aa(k,c) t2_aaaa(c,a,j,k) t1_1_aa(b,i)
        //          += -1.00 P(i,j) P(a,b) <a,k||j,c>_abab t2_abab(b,c,i,k)
        //          += +1.00 P(i,j) P(a,b) d-_aa(a,j) t1_1_aa(b,i)
        //          += +1.00 P(i,j) P(a,b) <k,a||c,j>_aaaa t2_aaaa(c,b,i,k)
        rt2_aaaa("a,b,i,j") += tmps_["215_aaaa_vovo"]("a,j,b,i");
        rt2_aaaa("a,b,i,j") -= tmps_["215_aaaa_vovo"]("a,i,b,j");
        rt2_aaaa("a,b,i,j") -= tmps_["215_aaaa_vovo"]("b,j,a,i");
        rt2_aaaa("a,b,i,j") += tmps_["215_aaaa_vovo"]("b,i,a,j");
        tmps_["215_aaaa_vovo"].~TArrayD();

        // flops: o0v2  = o2v3
        //  mems: o0v2  = o0v2
        tmps_["216_aa_vv"]("a,b")  = eri["aaaa_oovv"]("i,j,a,c") * t2["aaaa"]("c,b,j,i");

        // flops: o2v2  = o4v2 o4v2 o2v2 o2v4 o2v2 o2v2 o2v2 o2v2 o2v3 o2v3 o2v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["217_aaaa_oovv"]("i,j,a,b")  = -0.50 * tmps_["83_aaaa_oooo"]("l,k,i,j") * t2["aaaa"]("a,b,k,l");
        tmps_["217_aaaa_oovv"]("i,j,a,b") -= eri["aaaa_oooo"]("l,k,i,j") * t2["aaaa"]("a,b,k,l");
        tmps_["217_aaaa_oovv"]("i,j,a,b") += 2.00 * eri["aaaa_vvoo"]("a,b,i,j");
        tmps_["217_aaaa_oovv"]("i,j,a,b") += eri["aaaa_vvvv"]("a,b,c,d") * t2["aaaa"]("c,d,i,j");
        tmps_["217_aaaa_oovv"]("i,j,a,b") -= 2.00 * scalars_["5"] * t2_1["aaaa"]("a,b,i,j");
        tmps_["217_aaaa_oovv"]("i,j,a,b") += tmps_["216_aa_vv"]("c,b") * t2["aaaa"]("c,a,i,j");
        tmps_["217_aaaa_oovv"]("i,j,a,b") += tmps_["4_aa_vv"]("d,a") * t2["aaaa"]("d,b,i,j");
        tmps_["216_aa_vv"].~TArrayD();
        tmps_["83_aaaa_oooo"].~TArrayD();

        // rt2_aaaa += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,a,l,k) t2_aaaa(d,b,i,j)
        //          += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,a,i,j) t2_aaaa(d,b,l,k)
        //          += +1.00 <a,b||i,j>_aaaa
        //          += +0.50 <l,k||i,j>_aaaa t2_aaaa(a,b,l,k)
        //          += +0.50 <a,b||c,d>_aaaa t2_aaaa(c,d,i,j)
        //          += -1.00 d-_aa(k,k) t2_1_aaaa(a,b,i,j)
        //          += -1.00 d-_bb(k,k) t2_1_aaaa(a,b,i,j)
        //          += +0.25 <l,k||c,d>_aaaa t2_aaaa(a,b,l,k) t2_aaaa(c,d,i,j)
        rt2_aaaa("a,b,i,j") += 0.50 * tmps_["217_aaaa_oovv"]("i,j,a,b");
        tmps_["217_aaaa_oovv"].~TArrayD();

        // flops: o2v2  = o0v2 o2v3
        //  mems: o2v2  = o0v2 o2v2
        tmps_["218_aaaa_vvoo"]("a,b,i,j")  = (-1.00 * tmps_["10_aa_vv"]("a,c") + f["aa_vv"]("a,c")) * t2["aaaa"]("c,b,i,j");

        // rt2_aaaa += +1.00 P(a,b) f_aa(a,c) t2_aaaa(c,b,i,j)
        //          += -1.00 P(a,b) d-_aa(a,c) t2_aaaa(c,b,i,j) t0_1
        rt2_aaaa("a,b,i,j") += tmps_["218_aaaa_vvoo"]("a,b,i,j");
        rt2_aaaa("a,b,i,j") -= tmps_["218_aaaa_vvoo"]("b,a,i,j");
        tmps_["218_aaaa_vvoo"].~TArrayD();

        // flops: o2v2  = o2v3
        //  mems: o2v2  = o2v2
        tmps_["219_aaaa_voov"]("a,i,j,b")  = tmps_["5_aa_vv"]("c,b") * t2["aaaa"]("c,a,i,j");

        // rt2_aaaa += -0.50 <l,k||d,c>_abab t2_abab(a,c,l,k) t2_aaaa(d,b,i,j)
        //          += -0.50 <k,l||d,c>_abab t2_abab(a,c,k,l) t2_aaaa(d,b,i,j)
        rt2_aaaa("a,b,i,j") -= tmps_["219_aaaa_voov"]("b,i,j,a");

        // rt2_aaaa += +0.50 <l,k||c,d>_abab t2_aaaa(c,a,i,j) t2_abab(b,d,l,k)
        //          += +0.50 <k,l||c,d>_abab t2_aaaa(c,a,i,j) t2_abab(b,d,k,l)
        rt2_aaaa("a,b,i,j") += tmps_["219_aaaa_voov"]("a,i,j,b");
        tmps_["219_aaaa_voov"].~TArrayD();

        // flops: o0v2  = o2v3
        //  mems: o0v2  = o0v2
        tmps_["220_bb_vv"]("a,b")  = eri["bbbb_oovv"]("i,j,a,c") * t2["bbbb"]("c,b,j,i");

        // flops: o2v2  = o3v3 o3v3 o3v2 o2v2 o2v2 o4v2 o2v2 o3v2 o2v2 o2v4 o2v2 o3v3 o2v2 o2v2 o2v2 o2v2 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o2v2 o2v2 o2v2 o2v3 o2v2 o2v3 o2v2 o3v3 o2v2 o2v3 o2v2 o2v2 o2v2 o3v3 o2v2 o3v2 o2v2 o3v2 o2v2 o2v2 o2v2 o2v3 o2v2 o2v3 o2v3 o2v2 o3v2 o2v2 o2v3 o2v2 o3v2 o2v2 o3v2 o2v2 o3v3 o2v2 o3v3 o2v2 o3v3 o2v2 o4v2 o2v2 o2v2 o3v2 o2v2 o2v3 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["221_aabb_vovo"]("a,i,b,j")  = eri["abba_vovo"]("a,l,c,i") * t2["bbbb"]("c,b,j,l");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= eri["abab_vovo"]("a,l,e,j") * t2["abab"]("e,b,i,l");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= f["aa_oo"]("m,i") * t2["abab"]("a,b,m,j");
        tmps_["221_aabb_vovo"]("a,i,b,j") += eri["abab_oooo"]("k,l,i,j") * t2["abab"]("a,b,k,l");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= f["bb_oo"]("l,j") * t2["abab"]("a,b,i,l");
        tmps_["221_aabb_vovo"]("a,i,b,j") += eri["abab_vvvv"]("a,b,e,f") * t2["abab"]("e,f,i,j");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= eri["aaaa_vovo"]("a,m,e,i") * t2["abab"]("e,b,m,j");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= scalars_["5"] * t2_1["abab"]("a,b,i,j");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= dp["bb_vo"]("b,j") * t1_1["aa"]("a,i");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= eri["bbbb_vovo"]("b,l,c,j") * t2["abab"]("a,c,i,l");
        tmps_["221_aabb_vovo"]("a,i,b,j") += eri["baab_vovo"]("b,m,e,j") * t2["aaaa"]("e,a,i,m");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= eri["baba_vovo"]("b,m,c,i") * t2["abab"]("a,c,m,j");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= dp["aa_vo"]("a,i") * t1_1["bb"]("b,j");
        tmps_["221_aabb_vovo"]("a,i,b,j") += eri["abab_vvoo"]("a,b,i,j");
        tmps_["221_aabb_vovo"]("a,i,b,j") += f["aa_vv"]("a,e") * t2["abab"]("e,b,i,j");
        tmps_["221_aabb_vovo"]("a,i,b,j") += f["bb_vv"]("b,c") * t2["abab"]("a,c,i,j");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= t2["abab"]("d,b,k,j") * tmps_["1_aaaa_ovvo"]("k,d,a,i");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= t2["abab"]("e,b,i,j") * tmps_["10_aa_vv"]("a,e");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= t1_1["bb"]("b,j") * tmps_["75_aa_vo"]("a,i");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= t2["bbbb"]("f,b,j,n") * tmps_["3_bbaa_ovvo"]("n,f,a,i");
        tmps_["221_aabb_vovo"]("a,i,b,j") += t2["abab"]("a,b,i,l") * tmps_["61_bb_oo"]("l,j");
        tmps_["221_aabb_vovo"]("a,i,b,j") += t2["abab"]("a,b,m,j") * tmps_["41_aa_oo"]("m,i");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= t1_1["aa"]("a,i") * tmps_["102_bb_vo"]("b,j");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= t2["abab"]("a,c,i,j") * tmps_["25_bb_vv"]("b,c");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= t2["abab"]("d,b,i,j") * tmps_["5_aa_vv"]("d,a");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= 0.50 * tmps_["220_bb_vv"]("c,b") * t2["abab"]("a,c,i,j");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= tmps_["22_bb_oo"]("n,j") * t2["abab"]("a,b,i,n");
        tmps_["221_aabb_vovo"]("a,i,b,j") += 0.50 * t2["abab"]("d,b,i,j") * tmps_["4_aa_vv"]("d,a");
        tmps_["221_aabb_vovo"]("a,i,b,j") += 0.50 * tmps_["6_aa_oo"]("k,i") * t2["abab"]("a,b,k,j");
        tmps_["221_aabb_vovo"]("a,i,b,j") += 0.50 * tmps_["23_bb_oo"]("n,j") * t2["abab"]("a,b,i,n");
        tmps_["221_aabb_vovo"]("a,i,b,j") += tmps_["2_bbaa_ovvo"]("n,f,a,i") * t2["bbbb"]("f,b,j,n");
        tmps_["221_aabb_vovo"]("a,i,b,j") += tmps_["18_bbbb_ovvo"]("l,c,b,j") * t2["abab"]("a,c,i,l");
        tmps_["221_aabb_vovo"]("a,i,b,j") += tmps_["86_abba_ovvo"]("m,c,b,i") * t2["abab"]("a,c,m,j");
        tmps_["221_aabb_vovo"]("a,i,b,j") += tmps_["88_abab_oooo"]("k,l,i,j") * t2["abab"]("a,b,k,l");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= tmps_["7_aa_oo"]("k,i") * t2["abab"]("a,b,k,j");
        tmps_["221_aabb_vovo"]("a,i,b,j") -= tmps_["20_bb_vv"]("c,b") * t2["abab"]("a,c,i,j");
        tmps_["88_abab_oooo"].~TArrayD();
        tmps_["86_abba_ovvo"].~TArrayD();
        tmps_["75_aa_vo"].~TArrayD();
        tmps_["41_aa_oo"].~TArrayD();
        tmps_["10_aa_vv"].~TArrayD();
        tmps_["7_aa_oo"].~TArrayD();
        tmps_["6_aa_oo"].~TArrayD();
        tmps_["5_aa_vv"].~TArrayD();
        tmps_["4_aa_vv"].~TArrayD();
        tmps_["3_bbaa_ovvo"].~TArrayD();
        tmps_["2_bbaa_ovvo"].~TArrayD();
        tmps_["1_aaaa_ovvo"].~TArrayD();

        // rt2_abab += +1.00 <k,l||c,d>_abab t2_aaaa(c,a,i,k) t2_bbbb(d,b,j,l)
        //          += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,a,l,k) t2_abab(d,b,i,j)
        //          += +1.00 <l,k||d,c>_abab t2_abab(a,c,i,k) t2_abab(d,b,l,j)
        //          += +0.50 <l,k||c,d>_bbbb t2_abab(a,c,i,j) t2_bbbb(d,b,l,k)
        //          += -0.50 <l,k||c,d>_bbbb t2_abab(a,b,i,l) t2_bbbb(c,d,j,k)
        //          += +1.00 <k,l||d,c>_abab t2_abab(a,c,k,j) t2_abab(d,b,i,l)
        //          += +0.25 <l,k||c,d>_abab t2_abab(a,b,l,k) t2_abab(c,d,i,j)
        //          += +0.25 <l,k||d,c>_abab t2_abab(a,b,l,k) t2_abab(d,c,i,j)
        //          += +0.25 <k,l||c,d>_abab t2_abab(a,b,k,l) t2_abab(c,d,i,j)
        //          += +0.25 <k,l||d,c>_abab t2_abab(a,b,k,l) t2_abab(d,c,i,j)
        //          += -0.50 <k,l||c,d>_abab t2_abab(a,b,i,l) t2_abab(c,d,k,j)
        //          += -0.50 <k,l||d,c>_abab t2_abab(a,b,i,l) t2_abab(d,c,k,j)
        //          += -0.50 <l,k||c,d>_aaaa t2_abab(a,b,l,j) t2_aaaa(c,d,i,k)
        //          += -0.50 <l,k||d,c>_abab t2_abab(a,c,l,k) t2_abab(d,b,i,j)
        //          += -0.50 <k,l||d,c>_abab t2_abab(a,c,k,l) t2_abab(d,b,i,j)
        //          += -1.00 f_aa(k,i) t2_abab(a,b,k,j)
        //          += -1.00 f_bb(k,j) t2_abab(a,b,i,k)
        //          += -1.00 <a,k||c,j>_abab t2_abab(c,b,i,k)
        //          += -1.00 <a,k||i,c>_abab t2_bbbb(c,b,j,k)
        //          += +0.50 <l,k||i,j>_abab t2_abab(a,b,l,k)
        //          += +0.50 <k,l||i,j>_abab t2_abab(a,b,k,l)
        //          += +1.00 <k,a||c,i>_aaaa t2_abab(c,b,k,j)
        //          += +0.50 <a,b||c,d>_abab t2_abab(c,d,i,j)
        //          += +0.50 <a,b||d,c>_abab t2_abab(d,c,i,j)
        //          += -1.00 d-_aa(k,k) t2_1_abab(a,b,i,j)
        //          += -1.00 d-_bb(k,k) t2_1_abab(a,b,i,j)
        //          += -1.00 d-_bb(b,j) t1_1_aa(a,i)
        //          += +1.00 <k,b||c,j>_bbbb t2_abab(a,c,i,k)
        //          += -1.00 <k,b||c,j>_abab t2_aaaa(c,a,i,k)
        //          += -1.00 <k,b||i,c>_abab t2_abab(a,c,k,j)
        //          += -1.00 d-_aa(a,i) t1_1_bb(b,j)
        //          += +1.00 f_aa(a,c) t2_abab(c,b,i,j)
        //          += +1.00 <a,b||i,j>_abab
        //          += +1.00 f_bb(b,c) t2_abab(a,c,i,j)
        //          += +1.00 <l,k||c,d>_aaaa t2_aaaa(c,a,i,k) t2_abab(d,b,l,j)
        //          += -1.00 d-_aa(a,c) t2_abab(c,b,i,j) t0_1
        //          += -1.00 d-_bb(k,c) t2_abab(a,c,i,k) t1_1_bb(b,j)
        //          += +1.00 d-_aa(k,c) t2_aaaa(c,a,i,k) t1_1_bb(b,j)
        //          += +1.00 <l,k||c,d>_bbbb t2_abab(a,c,i,k) t2_bbbb(d,b,j,l)
        //          += +1.00 d-_bb(k,j) t2_abab(a,b,i,k) t0_1
        //          += +1.00 d-_aa(k,i) t2_abab(a,b,k,j) t0_1
        //          += -1.00 d-_aa(k,c) t2_abab(c,b,k,j) t1_1_aa(a,i)
        //          += +1.00 d-_bb(k,c) t2_bbbb(c,b,j,k) t1_1_aa(a,i)
        //          += -1.00 d-_bb(b,c) t2_abab(a,c,i,j) t0_1
        //          += -0.50 <l,k||c,d>_abab t2_abab(a,b,l,j) t2_abab(c,d,i,k)
        //          += -0.50 <l,k||d,c>_abab t2_abab(a,b,l,j) t2_abab(d,c,i,k)
        //          += -0.50 <l,k||d,c>_abab t2_abab(a,c,i,j) t2_abab(d,b,l,k)
        //          += -0.50 <k,l||d,c>_abab t2_abab(a,c,i,j) t2_abab(d,b,k,l)
        rt2_abab("a,b,i,j") += tmps_["221_aabb_vovo"]("a,i,b,j");
        tmps_["221_aabb_vovo"].~TArrayD();

        // flops: o2v2  = o3v2 o3v2 o2v2 o3v3 o2v2 o3v3 o2v2 o3v2 o2v2 o3v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["222_bbbb_ovvo"]("i,a,b,j")  = -0.50 * tmps_["23_bb_oo"]("k,i") * t2["bbbb"]("a,b,j,k");
        tmps_["222_bbbb_ovvo"]("i,a,b,j") -= tmps_["61_bb_oo"]("m,i") * t2["bbbb"]("a,b,j,m");
        tmps_["222_bbbb_ovvo"]("i,a,b,j") += t2["bbbb"]("d,b,j,k") * tmps_["19_bbbb_ovvo"]("k,d,a,i");
        tmps_["222_bbbb_ovvo"]("i,a,b,j") += t2["abab"]("e,b,n,j") * tmps_["17_aabb_ovvo"]("n,e,a,i");
        tmps_["222_bbbb_ovvo"]("i,a,b,j") += f["bb_oo"]("m,i") * t2["bbbb"]("a,b,j,m");
        tmps_["222_bbbb_ovvo"]("i,a,b,j") += tmps_["22_bb_oo"]("k,i") * t2["bbbb"]("a,b,j,k");
        tmps_["61_bb_oo"].~TArrayD();
        tmps_["23_bb_oo"].~TArrayD();
        tmps_["22_bb_oo"].~TArrayD();
        tmps_["19_bbbb_ovvo"].~TArrayD();
        tmps_["17_aabb_ovvo"].~TArrayD();

        // rt2_bbbb += +1.00 P(i,j) <l,k||c,d>_bbbb t2_bbbb(c,a,j,k) t2_bbbb(d,b,i,l)
        //          += +1.00 P(i,j) d-_bb(k,j) t2_bbbb(a,b,i,k) t0_1
        //          += +1.00 P(i,j) <l,k||c,d>_aaaa t2_abab(c,a,k,j) t2_abab(d,b,l,i)
        //          += -1.00 P(i,j) f_bb(k,j) t2_bbbb(a,b,i,k)
        //          += -0.50 P(i,j) <k,l||c,d>_abab t2_bbbb(a,b,i,l) t2_abab(c,d,k,j)
        //          += -0.50 P(i,j) <k,l||d,c>_abab t2_bbbb(a,b,i,l) t2_abab(d,c,k,j)
        //          += -0.50 P(i,j) <l,k||c,d>_bbbb t2_bbbb(a,b,i,l) t2_bbbb(c,d,j,k)
        rt2_bbbb("a,b,i,j") -= tmps_["222_bbbb_ovvo"]("j,a,b,i");
        rt2_bbbb("a,b,i,j") += tmps_["222_bbbb_ovvo"]("i,a,b,j");
        tmps_["222_bbbb_ovvo"].~TArrayD();

        // flops: o2v2  = o2v3
        //  mems: o2v2  = o2v2
        tmps_["223_bbbb_voov"]("a,i,j,b")  = t2["bbbb"]("c,a,i,j") * tmps_["20_bb_vv"]("c,b");
        tmps_["20_bb_vv"].~TArrayD();

        // rt2_bbbb += -0.50 <l,k||c,d>_abab t2_abab(c,a,l,k) t2_bbbb(d,b,i,j)
        //          += -0.50 <k,l||c,d>_abab t2_abab(c,a,k,l) t2_bbbb(d,b,i,j)
        rt2_bbbb("a,b,i,j") -= tmps_["223_bbbb_voov"]("b,i,j,a");

        // rt2_bbbb += +0.50 <l,k||d,c>_abab t2_bbbb(c,a,i,j) t2_abab(d,b,l,k)
        //          += +0.50 <k,l||d,c>_abab t2_bbbb(c,a,i,j) t2_abab(d,b,k,l)
        rt2_bbbb("a,b,i,j") += tmps_["223_bbbb_voov"]("a,i,j,b");
        tmps_["223_bbbb_voov"].~TArrayD();

        // flops: o2v2  = o3v3
        //  mems: o2v2  = o2v2
        tmps_["224_bbbb_vovo"]("a,i,b,j")  = t2["bbbb"]("c,a,i,k") * tmps_["18_bbbb_ovvo"]("k,c,b,j");
        tmps_["18_bbbb_ovvo"].~TArrayD();

        // rt2_bbbb += +1.00 P(i,j) <l,k||d,c>_abab t2_bbbb(c,a,j,k) t2_abab(d,b,l,i)
        rt2_bbbb("a,b,i,j") += tmps_["224_bbbb_vovo"]("a,j,b,i");
        rt2_bbbb("a,b,i,j") -= tmps_["224_bbbb_vovo"]("a,i,b,j");

        // rt2_bbbb += +1.00 P(i,j) <k,l||c,d>_abab t2_abab(c,a,k,j) t2_bbbb(d,b,i,l)
        rt2_bbbb("a,b,i,j") += tmps_["224_bbbb_vovo"]("b,i,a,j");
        rt2_bbbb("a,b,i,j") -= tmps_["224_bbbb_vovo"]("b,j,a,i");
        tmps_["224_bbbb_vovo"].~TArrayD();

        // flops: o2v2  = o3v3 o2v2 o2v2 o3v3 o2v2 o2v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["225_bbbb_vovo"]("a,i,b,j")  = -1.00 * eri["bbbb_vovo"]("a,l,d,i") * t2["bbbb"]("d,b,j,l");
        tmps_["225_bbbb_vovo"]("a,i,b,j") += dp["bb_vo"]("a,i") * t1_1["bb"]("b,j");
        tmps_["225_bbbb_vovo"]("a,i,b,j") += eri["baab_vovo"]("a,k,c,i") * t2["abab"]("c,b,k,j");
        tmps_["225_bbbb_vovo"]("a,i,b,j") += t1_1["bb"]("b,j") * tmps_["102_bb_vo"]("a,i");
        tmps_["102_bb_vo"].~TArrayD();

        // rt2_bbbb += +1.00 P(i,j) P(a,b) d-_aa(k,c) t2_abab(c,a,k,j) t1_1_bb(b,i)
        //          += -1.00 P(i,j) P(a,b) d-_bb(k,c) t2_bbbb(c,a,j,k) t1_1_bb(b,i)
        //          += +1.00 P(i,j) P(a,b) d-_bb(a,j) t1_1_bb(b,i)
        //          += +1.00 P(i,j) P(a,b) <k,a||c,j>_bbbb t2_bbbb(c,b,i,k)
        //          += -1.00 P(i,j) P(a,b) <k,a||c,j>_abab t2_abab(c,b,k,i)
        rt2_bbbb("a,b,i,j") += tmps_["225_bbbb_vovo"]("a,j,b,i");
        rt2_bbbb("a,b,i,j") -= tmps_["225_bbbb_vovo"]("a,i,b,j");
        rt2_bbbb("a,b,i,j") -= tmps_["225_bbbb_vovo"]("b,j,a,i");
        rt2_bbbb("a,b,i,j") += tmps_["225_bbbb_vovo"]("b,i,a,j");
        tmps_["225_bbbb_vovo"].~TArrayD();

        // flops: o2v2  = o4v2 o4v2 o2v2 o2v4 o2v2 o2v2 o2v2 o2v2 o2v3 o2v2 o2v3 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["226_bbbb_oovv"]("i,j,a,b")  = -0.50 * tmps_["117_bbbb_oooo"]("k,l,i,j") * t2["bbbb"]("a,b,l,k");
        tmps_["226_bbbb_oovv"]("i,j,a,b") -= eri["bbbb_oooo"]("k,l,i,j") * t2["bbbb"]("a,b,l,k");
        tmps_["226_bbbb_oovv"]("i,j,a,b") += eri["bbbb_vvvv"]("a,b,d,c") * t2["bbbb"]("d,c,i,j");
        tmps_["226_bbbb_oovv"]("i,j,a,b") += 2.00 * eri["bbbb_vvoo"]("a,b,i,j");
        tmps_["226_bbbb_oovv"]("i,j,a,b") += -2.00 * scalars_["5"] * t2_1["bbbb"]("a,b,i,j");
        tmps_["226_bbbb_oovv"]("i,j,a,b") += tmps_["220_bb_vv"]("d,b") * t2["bbbb"]("d,a,i,j");
        tmps_["226_bbbb_oovv"]("i,j,a,b") += tmps_["21_bb_vv"]("c,a") * t2["bbbb"]("c,b,i,j");
        tmps_["220_bb_vv"].~TArrayD();
        tmps_["117_bbbb_oooo"].~TArrayD();
        tmps_["21_bb_vv"].~TArrayD();

        // rt2_bbbb += -0.50 <l,k||c,d>_bbbb t2_bbbb(c,a,i,j) t2_bbbb(d,b,l,k)
        //          += +0.50 <a,b||c,d>_bbbb t2_bbbb(c,d,i,j)
        //          += +0.50 <l,k||i,j>_bbbb t2_bbbb(a,b,l,k)
        //          += +1.00 <a,b||i,j>_bbbb
        //          += -1.00 d-_aa(k,k) t2_1_bbbb(a,b,i,j)
        //          += -1.00 d-_bb(k,k) t2_1_bbbb(a,b,i,j)
        //          += -0.50 <l,k||c,d>_bbbb t2_bbbb(c,a,l,k) t2_bbbb(d,b,i,j)
        //          += +0.25 <l,k||c,d>_bbbb t2_bbbb(a,b,l,k) t2_bbbb(c,d,i,j)
        rt2_bbbb("a,b,i,j") += 0.50 * tmps_["226_bbbb_oovv"]("i,j,a,b");
        tmps_["226_bbbb_oovv"].~TArrayD();

        // flops: o2v2  = o0v2 o2v3
        //  mems: o2v2  = o0v2 o2v2
        tmps_["227_bbbb_vvoo"]("a,b,i,j")  = (-1.00 * tmps_["25_bb_vv"]("a,c") + f["bb_vv"]("a,c")) * t2["bbbb"]("c,b,i,j");
        tmps_["25_bb_vv"].~TArrayD();

        // rt2_bbbb += +1.00 P(a,b) f_bb(a,c) t2_bbbb(c,b,i,j)
        //          += -1.00 P(a,b) d-_bb(a,c) t2_bbbb(c,b,i,j) t0_1
        rt2_bbbb("a,b,i,j") += tmps_["227_bbbb_vvoo"]("a,b,i,j");
        rt2_bbbb("a,b,i,j") -= tmps_["227_bbbb_vvoo"]("b,a,i,j");
        tmps_["227_bbbb_vvoo"].~TArrayD();

}