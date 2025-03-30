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
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 *  @END LICENSE
 */

#include "cc_cavity/include/ccsd/eom_ee_ccsd.h"

double* hilbert::EOM_EE_CCSD::build_ss_diagonal() {

    if (!cc_wfn_->has_t1_integrals_) cc_wfn_->transform_integrals(true);

    // extract 2-body amplitudes
    std::map<std::string, TA::TArrayD> t2 {
            {"aaaa", cc_wfn_->amplitudes_["t2_aaaa"]},
            {"abab", cc_wfn_->amplitudes_["t2_abab"]},
            {"bbbb", cc_wfn_->amplitudes_["t2_bbbb"]}
    };

    world_.gop.fence();

    TArrayMap &eri = cc_wfn_->V_blks_;
    TArrayMap &f = cc_wfn_->F_blks_;
    TArrayMap &Id = cc_wfn_->Id_blks_;

    // initialize diagonal blocks
    TArrayD    H_ss_aaaa = makeTensor(world_, {oa_,va_, oa_,va_}, true);
    TArrayD    H_ss_bbbb = makeTensor(world_, {ob_,vb_, ob_,vb_}, true);

    {

        scalars_["1"]  = 1.00 * dot(Id["aa_oo"]("j,j1"), f["aa_oo"]("j,j1"));
        scalars_["2"]  = 1.00 * dot(eri["abab_oovv"]("k,j,b,c"), t2["abab"]("b,c,k,j"));
        scalars_["3"]  = 1.00 * dot(eri["bbbb_oovv"]("j,k,b,c"), t2["bbbb"]("b,c,k,j"));
        scalars_["4"]  = 1.00 * dot(eri["aaaa_oovv"]("j,k,b,c"), t2["aaaa"]("b,c,k,j"));
        scalars_["5"]  = 1.00 * dot(Id["bb_oo"]("j,j1") * eri["abab_oooo"]("k,j,k1,j1"), Id["aa_oo"]("k,k1"));
        scalars_["6"]  = 1.00 * dot(Id["bb_oo"]("j,j1") * eri["bbbb_oooo"]("j,k,j1,k1"), Id["bb_oo"]("k,k1"));
        scalars_["7"]  = 1.00 * dot(Id["aa_oo"]("j,j1") * eri["aaaa_oooo"]("j,k,j1,k1"), Id["aa_oo"]("k,k1"));
        scalars_["8"]  = 1.00 * dot(Id["bb_oo"]("j,j1"), f["bb_oo"]("j,j1"));
        scalars_["9"]  = scalars_["1"];
        scalars_["9"] -= 0.25 * scalars_["3"];
        scalars_["9"] -= 0.25 * scalars_["4"];
        scalars_["9"] -= 0.50 * scalars_["6"];
        scalars_["9"] -= scalars_["5"];
        scalars_["9"] += scalars_["2"];
        scalars_["9"] -= 0.50 * scalars_["7"];
        scalars_["9"] += scalars_["8"];


        // flops: o0v2  = o2v3 o2v3 o0v2
        //  mems: o0v2  = o0v2 o0v2 o0v2
        tmps_["2_aa_vv"]("a,e")  = eri["abab_oovv"]("k,j,a,b") * t2["abab"]("e,b,k,j");
        tmps_["2_aa_vv"]("a,e") += 0.50 * eri["aaaa_oovv"]("j,k,a,b") * t2["aaaa"]("b,e,k,j");

        // flops: o2v0  = o2v0 o3v2 o2v0 o3v2 o2v0
        //  mems: o2v0  = o2v0 o2v0 o2v0 o2v0 o2v0
        tmps_["1_aa_oo"]("m,i")  = scalars_["9"] * Id["aa_oo"]("m,i");
        tmps_["1_aa_oo"]("m,i") -= eri["abab_oovv"]("i,j,b,c") * t2["abab"]("b,c,m,j");
        tmps_["1_aa_oo"]("m,i") -= 0.50 * eri["aaaa_oovv"]("i,j,b,c") * t2["aaaa"]("b,c,m,j");

        // flops: o2v2  = o2v2 o3v3 o3v3 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["3_aaaa_vvoo"]("a,e,m,i")  = Id["aa_oo"]("m,i") * tmps_["2_aa_vv"]("a,e");
        tmps_["3_aaaa_vvoo"]("a,e,m,i") += t2["aaaa"]("b,e,m,j") * eri["aaaa_oovv"]("i,j,a,b");
        tmps_["3_aaaa_vvoo"]("a,e,m,i") -= t2["abab"]("e,b,m,j") * eri["abab_oovv"]("i,j,a,b");
        tmps_["3_aaaa_vvoo"]("a,e,m,i") -= f["aa_vv"]("e,a") * Id["aa_oo"]("m,i");
        tmps_["3_aaaa_vvoo"]("a,e,m,i") += Id["aa_vv"]("e,a") * f["aa_oo"]("i,m");
        tmps_["3_aaaa_vvoo"]("a,e,m,i") += eri["aaaa_vovo"]("e,i,a,m");
        tmps_["3_aaaa_vvoo"]("a,e,m,i") -= Id["aa_vv"]("e,a") * tmps_["1_aa_oo"]("m,i");
        tmps_["2_aa_vv"].~TArrayD();
        tmps_["1_aa_oo"].~TArrayD();

        // H_ss_aaaa  = -0.50 d_aa(m,i) <k,j||a,b>_abab t2_abab(e,b,k,j)
        //           += -0.50 d_aa(m,i) <j,k||a,b>_abab t2_abab(e,b,j,k)
        //           += -0.50 d_aa(m,i) <k,j||b,a>_aaaa t2_aaaa(b,e,k,j)
        //           += +1.00 <i,j||b,a>_aaaa t2_aaaa(b,e,m,j)
        //           += +1.00 <i,j||a,b>_abab t2_abab(e,b,m,j)
        //           += +1.00 d_aa(m,i) f_aa(e,a)
        //           += -1.00 d_aa(e,a) f_aa(i,m)
        //           += +1.00 <i,e||a,m>_aaaa
        //           += +1.00 d_aa(e,a) d_aa(m,i) f_aa(j,j)
        //           += +0.25 d_aa(e,a) d_aa(m,i) <k,j||b,c>_bbbb t2_bbbb(b,c,k,j)
        //           += +0.25 d_aa(e,a) d_aa(m,i) <k,j||b,c>_aaaa t2_aaaa(b,c,k,j)
        //           += -0.50 d_aa(e,a) d_aa(m,i) <k,j||k,j>_bbbb
        //           += -0.50 d_aa(e,a) d_aa(m,i) <k,j||k,j>_abab
        //           += -0.50 d_aa(e,a) d_aa(m,i) <j,k||j,k>_abab
        //           += +0.25 d_aa(e,a) d_aa(m,i) <k,j||b,c>_abab t2_abab(b,c,k,j)
        //           += +0.25 d_aa(e,a) d_aa(m,i) <j,k||b,c>_abab t2_abab(b,c,j,k)
        //           += +0.25 d_aa(e,a) d_aa(m,i) <k,j||c,b>_abab t2_abab(c,b,k,j)
        //           += +0.25 d_aa(e,a) d_aa(m,i) <j,k||c,b>_abab t2_abab(c,b,j,k)
        //           += -0.50 d_aa(e,a) d_aa(m,i) <k,j||k,j>_aaaa
        //           += +1.00 d_aa(e,a) d_aa(m,i) f_bb(j,j)
        //           += -0.50 d_aa(e,a) <i,j||b,c>_abab t2_abab(b,c,m,j)
        //           += -0.50 d_aa(e,a) <i,j||c,b>_abab t2_abab(c,b,m,j)
        //           += -0.50 d_aa(e,a) <i,j||b,c>_aaaa t2_aaaa(b,c,m,j)
        H_ss_aaaa("a,e,i,m")  = -1.00 * tmps_["3_aaaa_vvoo"]("a,e,m,i");
        tmps_["3_aaaa_vvoo"].~TArrayD();

        // flops: o0v2  = o2v3 o2v3 o0v2
        //  mems: o0v2  = o0v2 o0v2 o0v2
        tmps_["5_bb_vv"]("a,e")  = eri["bbbb_oovv"]("j,k,a,b") * t2["bbbb"]("b,e,k,j");
        tmps_["5_bb_vv"]("a,e") += 2.00 * eri["abab_oovv"]("k,j,b,a") * t2["abab"]("b,e,k,j");

        // flops: o2v0  = o2v0 o3v2 o2v0 o3v2 o2v0
        //  mems: o2v0  = o2v0 o2v0 o2v0 o2v0 o2v0
        tmps_["4_bb_oo"]("m,i")  = scalars_["9"] * Id["bb_oo"]("m,i");
        tmps_["4_bb_oo"]("m,i") -= eri["abab_oovv"]("j,i,b,c") * t2["abab"]("b,c,j,m");
        tmps_["4_bb_oo"]("m,i") -= 0.50 * eri["bbbb_oovv"]("i,j,b,c") * t2["bbbb"]("b,c,m,j");

        // flops: o2v2  = o2v2 o2v2 o2v2 o3v3 o2v2 o3v3 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["6_bbbb_vvoo"]("e,a,m,i")  = Id["bb_vv"]("e,a") * tmps_["4_bb_oo"]("m,i");
        tmps_["6_bbbb_vvoo"]("e,a,m,i") -= 0.50 * Id["bb_oo"]("m,i") * tmps_["5_bb_vv"]("a,e");
        tmps_["6_bbbb_vvoo"]("e,a,m,i") += eri["abab_oovv"]("j,i,b,a") * t2["abab"]("b,e,j,m");
        tmps_["6_bbbb_vvoo"]("e,a,m,i") -= eri["bbbb_vovo"]("e,i,a,m");
        tmps_["6_bbbb_vvoo"]("e,a,m,i") -= eri["bbbb_oovv"]("i,j,a,b") * t2["bbbb"]("b,e,m,j");
        tmps_["6_bbbb_vvoo"]("e,a,m,i") -= Id["bb_vv"]("e,a") * f["bb_oo"]("i,m");
        tmps_["6_bbbb_vvoo"]("e,a,m,i") += Id["bb_oo"]("m,i") * f["bb_vv"]("e,a");
        tmps_["5_bb_vv"].~TArrayD();
        tmps_["4_bb_oo"].~TArrayD();

        // H_ss_bbbb  = +1.00 d_bb(e,a) d_bb(m,i) f_aa(j,j)
        //           += +0.25 d_bb(e,a) d_bb(m,i) <k,j||b,c>_bbbb t2_bbbb(b,c,k,j)
        //           += +0.25 d_bb(e,a) d_bb(m,i) <k,j||b,c>_aaaa t2_aaaa(b,c,k,j)
        //           += -0.50 d_bb(e,a) d_bb(m,i) <k,j||k,j>_bbbb
        //           += -0.50 d_bb(e,a) d_bb(m,i) <k,j||k,j>_abab
        //           += -0.50 d_bb(e,a) d_bb(m,i) <j,k||j,k>_abab
        //           += +0.25 d_bb(e,a) d_bb(m,i) <k,j||b,c>_abab t2_abab(b,c,k,j)
        //           += +0.25 d_bb(e,a) d_bb(m,i) <j,k||b,c>_abab t2_abab(b,c,j,k)
        //           += +0.25 d_bb(e,a) d_bb(m,i) <k,j||c,b>_abab t2_abab(c,b,k,j)
        //           += +0.25 d_bb(e,a) d_bb(m,i) <j,k||c,b>_abab t2_abab(c,b,j,k)
        //           += -0.50 d_bb(e,a) d_bb(m,i) <k,j||k,j>_aaaa
        //           += +1.00 d_bb(e,a) d_bb(m,i) f_bb(j,j)
        //           += -0.50 d_bb(e,a) <j,i||b,c>_abab t2_abab(b,c,j,m)
        //           += -0.50 d_bb(e,a) <j,i||c,b>_abab t2_abab(c,b,j,m)
        //           += -0.50 d_bb(e,a) <i,j||b,c>_bbbb t2_bbbb(b,c,m,j)
        //           += -0.50 d_bb(m,i) <k,j||b,a>_bbbb t2_bbbb(b,e,k,j)
        //           += -0.50 d_bb(m,i) <k,j||b,a>_abab t2_abab(b,e,k,j)
        //           += -0.50 d_bb(m,i) <j,k||b,a>_abab t2_abab(b,e,j,k)
        //           += +1.00 <j,i||b,a>_abab t2_abab(b,e,j,m)
        //           += +1.00 <i,e||a,m>_bbbb
        //           += +1.00 <i,j||b,a>_bbbb t2_bbbb(b,e,m,j)
        //           += -1.00 d_bb(e,a) f_bb(i,m)
        //           += +1.00 d_bb(m,i) f_bb(e,a)
        H_ss_bbbb("a,e,i,m")  = tmps_["6_bbbb_vvoo"]("e,a,m,i");
        tmps_["6_bbbb_vvoo"].~TArrayD();
    }

    // pluck out diagonals
    auto * ss_diag = (double*) calloc(singleDim_, sizeof(double));

    size_t id = 0;
    // extract singles
    size_t oa = oa_, ob = ob_, va = va_, vb = vb_;
    foreach_inplace(H_ss_aaaa, [id,oa,va,ob,vb,ss_diag](auto &tile){
        for (auto &x : tile.range()) {
            size_t m = x[0], e = x[1], i = x[2], a = x[3];
            if (m == i && e == a) {
                size_t index = e * oa + i + id;
                ss_diag[index] = tile[x];
            }
        }
    }); id += oa_*va_;

    foreach_inplace(H_ss_bbbb, [id,oa,va,ob,vb,ss_diag](auto &tile){
        for (auto &x : tile.range()) {
            size_t m = x[0], e = x[1], i = x[2], a = x[3];
            if (m == i && e == a) {
                size_t index = e * ob + i + id;
                ss_diag[index] = tile[x];
            }
        }
    }); id += ob_*vb_;

    tmps_.clear();
    world_.gop.fence(); // TODO: broadcast to all nodes
    return ss_diag;

}

void hilbert::EOM_EE_CCSD::build_common_ops() {

    if (!cc_wfn_->has_t1_integrals_) cc_wfn_->transform_integrals(true);

    TArrayMap &eri = cc_wfn_->V_blks_;
    TArrayMap &f = cc_wfn_->F_blks_;
    TArrayMap &Id = cc_wfn_->Id_blks_;

    TArrayMap t2 = {
            {"aaaa", cc_wfn_->amplitudes_["t2_aaaa"]},
            {"abab", cc_wfn_->amplitudes_["t2_abab"]},
            {"bbbb", cc_wfn_->amplitudes_["t2_bbbb"]}
    };

    {
        scalars_["1"]  = 1.00 * dot(Id["bb_oo"]("j,j1") * eri["bbbb_oooo"]("j,k,j1,k1"), Id["bb_oo"]("k,k1"));
        scalars_["2"]  = 1.00 * dot(Id["aa_oo"]("j,j1") * eri["aaaa_oooo"]("j,k,j1,k1"), Id["aa_oo"]("k,k1"));
        scalars_["3"]  = 1.00 * dot(Id["aa_oo"]("j,j1"), f["aa_oo"]("j,j1"));
        scalars_["4"]  = 1.00 * dot(eri["abab_oovv"]("k,j,b,c"), t2["abab"]("b,c,k,j"));
        scalars_["5"]  = 1.00 * dot(eri["bbbb_oovv"]("j,k,b,c"), t2["bbbb"]("b,c,k,j"));
        scalars_["6"]  = 1.00 * dot(Id["bb_oo"]("j,j1"), f["bb_oo"]("j,j1"));
        scalars_["7"]  = 1.00 * dot(eri["aaaa_oovv"]("j,k,b,c"), t2["aaaa"]("b,c,k,j"));
        scalars_["8"]  = 1.00 * dot(Id["bb_oo"]("j,j1") * eri["abab_oooo"]("k,j,k1,j1"), Id["aa_oo"]("k,k1"));
        scalars_["9"]  = scalars_["5"];
        scalars_["9"] -= 4.00 * scalars_["3"];
        scalars_["9"] -= 4.00 * scalars_["6"];
        scalars_["9"] += 2.00 * scalars_["2"];
        scalars_["9"] += 2.00 * scalars_["1"];
        scalars_["9"] += scalars_["7"];
        scalars_["9"] += 4.00 * scalars_["8"];
        scalars_["9"] -= 4.00 * scalars_["4"];
    }

    {
        reused_["1_bbbb_vvvo"]("a,d,b,j")  = eri["baab_vovv"]("a,k,c,d") * t2["abab"]("c,b,k,j");
        reused_["2_bbbb_vvvo"]("a,d,b,j")  = eri["bbbb_vovv"]("a,k,c,d") * t2["bbbb"]("c,b,j,k");
        reused_["3_bbbb_vvvo"]("a,d,b,j")  = reused_["2_bbbb_vvvo"]("a,d,b,j");
        reused_["3_bbbb_vvvo"]("a,d,b,j") -= reused_["1_bbbb_vvvo"]("a,d,b,j");
        reused_["4_bbbb_voov"]("c,k,j,b")  = eri["abab_oovv"]("l,j,d,b") * t2["abab"]("d,c,l,k");
        reused_["5_bbaa_voov"]("c,k,j,a")  = eri["aaaa_oovv"]("j,l,a,d") * t2["abab"]("d,c,l,k");
        reused_["6_bbaa_voov"]("b,j,l,d")  = eri["aaaa_oovv"]("k,l,c,d") * t2["abab"]("c,b,k,j");
        reused_["7_aabb_voov"]("c,k,j,b")  = eri["abab_oovv"]("l,j,d,b") * t2["aaaa"]("d,c,k,l");
        reused_["8_abab_vovv"]("a,j,c,b")  = eri["abab_oovo"]("l,k,a,j") * t2["abab"]("c,b,l,k");
        reused_["9_abab_vooo"]("b,i,j,k")  = eri["abab_vovv"]("b,i,c,d") * t2["abab"]("c,d,j,k");
        reused_["10_aaaa_vovv"]("a,j,c,b")  = eri["aaaa_oovo"]("k,l,a,j") * t2["aaaa"]("c,b,l,k");
        reused_["11_aaaa_vooo"]("b,i,j,k")  = eri["aaaa_vovv"]("b,i,c,d") * t2["aaaa"]("c,d,j,k");
        reused_["12_aaaa_voov"]("a,j,l,d")  = eri["aaaa_oovv"]("k,l,c,d") * t2["aaaa"]("c,a,j,k");
        reused_["13_aaaa_voov"]("c,k,j,a")  = eri["aaaa_oovv"]("j,l,a,d") * t2["aaaa"]("d,c,k,l");
        reused_["14_bbbb_vovv"]("b,j,c,a")  = eri["bbbb_vovv"]("c,k,a,d") * t2["bbbb"]("d,b,j,k");
        reused_["15_bbaa_voov"]("c,k,j,a")  = eri["abab_oovv"]("j,l,a,d") * t2["bbbb"]("d,c,k,l");
        reused_["16_aabb_vovv"]("b,j,c,a")  = eri["bbbb_vovv"]("c,k,a,d") * t2["abab"]("b,d,j,k");
        reused_["17_abba_vvvo"]("c,a,b,j")  = eri["abab_vovv"]("c,k,d,a") * t2["abab"]("d,b,j,k");
        reused_["18_bbaa_vvvo"]("b,d,a,i")  = t2["aaaa"]("c,a,i,k") * eri["baab_vovv"]("b,k,c,d");
        reused_["19_bbaa_vvvo"]("b,d,a,i")  = t2["abab"]("a,c,i,k") * eri["bbbb_vovv"]("b,k,c,d");
        reused_["20_bbaa_vvvo"]("b,d,a,i")  = reused_["18_bbaa_vvvo"]("b,d,a,i");
        reused_["20_bbaa_vvvo"]("b,d,a,i") -= reused_["19_bbaa_vvvo"]("b,d,a,i");
        reused_["20_bbaa_vvvo"]("b,d,a,i") -= reused_["17_abba_vvvo"]("a,d,b,i");
        reused_["21_aabb_vvvo"]("a,d,b,j")  = eri["aaaa_vovv"]("a,k,c,d") * t2["abab"]("c,b,k,j");
        reused_["22_aabb_vvvo"]("c,a,b,j")  = eri["abab_vovv"]("c,k,a,d") * t2["bbbb"]("d,b,j,k");
        reused_["23_baab_vvvo"]("b,d,a,j")  = t2["abab"]("a,c,k,j") * eri["baab_vovv"]("b,k,d,c");
        reused_["24_baab_vvvo"]("b,d,a,j")  = reused_["23_baab_vvvo"]("b,d,a,j");
        reused_["24_baab_vvvo"]("b,d,a,j") -= reused_["21_aabb_vvvo"]("a,d,b,j");
        reused_["24_baab_vvvo"]("b,d,a,j") -= reused_["22_aabb_vvvo"]("a,d,b,j");
        reused_["25_aabb_vvvo"]("c,a,b,j")  = eri["aaaa_vovv"]("c,k,a,d") * t2["abab"]("d,b,k,j");
        reused_["26_aaaa_oooo"]("i,j,k,l")  = eri["aaaa_oovv"]("k,l,c,d") * t2["aaaa"]("c,d,i,j");
        reused_["27_bbaa_vooo"]("b,j,i,k")  = eri["aaaa_oovo"]("i,l,c,k") * t2["abab"]("c,b,l,j");
        reused_["28_aabb_oovo"]("l,i,b,j")  = t2["abab"]("c,b,k,j") * eri["aaaa_oovo"]("k,l,c,i");
        reused_["29_aaaa_oovo"]("l,j,a,i")  = t2["abab"]("a,c,i,k") * eri["abba_oovo"]("l,k,c,j");
        reused_["30_abab_oooo"]("i,j,k,l")  = eri["abab_oovv"]("i,j,c,d") * t2["abab"]("c,d,k,l");
        reused_["31_baab_oovo"]("l,i,a,j")  = t2["abab"]("a,c,k,j") * eri["abba_oovo"]("k,l,c,i");
        reused_["32_abba_oovo"]("l,j,b,i")  = t2["abab"]("c,b,i,k") * eri["abab_oovo"]("l,k,c,j");
        reused_["33_aabb_vooo"]("b,j,i,k")  = eri["bbbb_oovo"]("i,l,c,k") * t2["abab"]("b,c,j,l");
        reused_["34_bbaa_oovo"]("l,j,a,i")  = t2["abab"]("a,c,i,k") * eri["bbbb_oovo"]("k,l,c,j");
        reused_["35_bbbb_oovo"]("l,j,a,i")  = t2["abab"]("c,a,k,i") * eri["abab_oovo"]("k,l,c,j");
        reused_["36_aaaa_vooo"]("b,j,i,k")  = eri["aaaa_oovo"]("i,l,c,k") * t2["aaaa"]("c,b,j,l");
        reused_["37_aaaa_oovo"]("l,j,a,i")  = t2["aaaa"]("c,a,i,k") * eri["aaaa_oovo"]("k,l,c,j");
        reused_["38_bbaa_oovo"]("l,j,a,i")  = t2["aaaa"]("c,a,i,k") * eri["abab_oovo"]("k,l,c,j");
        reused_["39_bbbb_vovv"]("a,j,c,b")  = eri["bbbb_oovo"]("k,l,a,j") * t2["bbbb"]("c,b,l,k");
        reused_["40_bbbb_vooo"]("b,i,j,k")  = eri["bbbb_vovv"]("b,i,c,d") * t2["bbbb"]("c,d,j,k");
        reused_["41_bbbb_voov"]("c,k,j,b")  = eri["bbbb_oovv"]("j,l,b,d") * t2["bbbb"]("d,c,k,l");
        reused_["42_bbbb_voov"]("b,j,l,d")  = eri["bbbb_oovv"]("k,l,c,d") * t2["bbbb"]("c,b,j,k");
        reused_["43_baab_vooo"]("b,i,j,k")  = eri["baab_vovv"]("b,i,c,d") * t2["abab"]("c,d,j,k");
        reused_["44_baab_voov"]("c,k,i,b")  = eri["abab_oovv"]("i,l,d,b") * t2["abab"]("d,c,k,l");
        reused_["45_baab_vovv"]("c,i,a,b")  = t2["abab"]("a,b,l,k") * eri["abba_oovo"]("l,k,c,i");
        reused_["46_aabb_voov"]("c,k,j,b")  = eri["bbbb_oovv"]("j,l,b,d") * t2["abab"]("c,d,k,l");
        reused_["47_aabb_voov"]("a,j,l,d")  = eri["bbbb_oovv"]("k,l,c,d") * t2["abab"]("a,c,j,k");
        reused_["48_aa_vo"]("b,j")  = eri["abab_vovv"]("b,k,c,d") * t2["abab"]("c,d,j,k");
        reused_["49_bb_ov"]("j,b")  = t2["abab"]("c,b,l,k") * eri["abab_oovo"]("l,k,c,j");
        reused_["50_bb_oo"]("i,k")  = eri["abab_oovv"]("j,k,b,c") * t2["abab"]("b,c,j,i");
        reused_["51_aaaa_vooo"]("b,j,k,i")  = f["aa_ov"]("i,c") * t2["aaaa"]("c,b,j,k");
        reused_["52_aa_oo"]("i,j")  = eri["aaaa_oovv"]("i,k,b,c") * t2["aaaa"]("b,c,j,k");
        reused_["53_aa_oo"]("i,k")  = eri["aaaa_oovv"]("j,k,b,c") * t2["aaaa"]("b,c,i,j");
        reused_["54_aa_ov"]("j,a")  = t2["aaaa"]("c,a,l,k") * eri["aaaa_oovo"]("k,l,c,j");
        reused_["55_bb_vv"]("a,b")  = eri["bbbb_oovv"]("j,k,a,c") * t2["bbbb"]("c,b,k,j");
        reused_["56_bb_vv"]("a,c")  = eri["bbbb_oovv"]("j,k,b,c") * t2["bbbb"]("b,a,k,j");
        reused_["57_bb_vo"]("b,j")  = eri["bbbb_vovv"]("b,k,c,d") * t2["bbbb"]("c,d,j,k");
        reused_["58_bb_vo"]("b,j")  = eri["baab_vovv"]("b,k,c,d") * t2["abab"]("c,d,k,j");
        reused_["59_bb_vv"]("a,b")  = eri["abab_oovv"]("k,j,c,a") * t2["abab"]("c,b,k,j");
        reused_["60_aabb_vooo"]("b,j,k,i")  = f["bb_ov"]("i,c") * t2["abab"]("b,c,j,k");
        reused_["61_aa_vv"]("a,c")  = eri["abab_oovv"]("k,j,c,b") * t2["abab"]("a,b,k,j");
        reused_["62_aa_vv"]("a,b")  = eri["aaaa_oovv"]("j,k,a,c") * t2["aaaa"]("c,b,k,j");
        reused_["63_aa_vv"]("a,c")  = eri["aaaa_oovv"]("j,k,b,c") * t2["aaaa"]("b,a,k,j");
        reused_["64_aa_vo"]("b,j")  = eri["aaaa_vovv"]("b,k,c,d") * t2["aaaa"]("c,d,j,k");
        reused_["65_aabb_oovo"]("l,i,b,j")  = t2["bbbb"]("c,b,j,k") * eri["abba_oovo"]("l,k,c,i");
        reused_["66_bbbb_vooo"]("b,j,i,k")  = eri["bbbb_oovo"]("i,l,c,k") * t2["bbbb"]("c,b,j,l");
        reused_["67_bbbb_oovo"]("l,j,a,i")  = t2["bbbb"]("c,a,i,k") * eri["bbbb_oovo"]("k,l,c,j");
        reused_["68_bbbb_oooo"]("i,j,k,l")  = eri["bbbb_oovv"]("k,l,c,d") * t2["bbbb"]("c,d,i,j");
        reused_["69_aaaa_vvvo"]("a,d,b,j")  = eri["abab_vovv"]("a,k,d,c") * t2["abab"]("b,c,j,k");
        reused_["70_aaaa_vovv"]("b,j,c,a")  = eri["aaaa_vovv"]("c,k,a,d") * t2["aaaa"]("d,b,j,k");
        reused_["71_bbbb_vooo"]("b,j,k,i")  = f["bb_ov"]("i,c") * t2["bbbb"]("c,b,j,k");
        reused_["72_bbbb_vooo"]("b,i,j,k")  = reused_["40_bbbb_vooo"]("b,i,j,k");
        reused_["72_bbbb_vooo"]("b,i,j,k") -= 2.00 * reused_["71_bbbb_vooo"]("b,j,k,i");
        reused_["73_aaaa_vvvo"]("a,d,b,j")  = eri["aaaa_vovv"]("a,k,c,d") * t2["aaaa"]("c,b,j,k");
        reused_["74_aaaa_vvvo"]("a,d,b,j")  = reused_["69_aaaa_vvvo"]("a,d,b,j");
        reused_["74_aaaa_vvvo"]("a,d,b,j") += reused_["73_aaaa_vvvo"]("a,d,b,j");
        reused_["75_bb_ov"]("j,b")  = t2["bbbb"]("c,b,l,k") * eri["bbbb_oovo"]("k,l,c,j");
        reused_["76_bb_vo"]("b,j")  = f["aa_ov"]("k,c") * t2["abab"]("c,b,k,j");
        reused_["77_bb_vo"]("b,j")  = f["bb_ov"]("k,c") * t2["bbbb"]("c,b,j,k");
        reused_["78_bb_vo"]("b,j")  = reused_["77_bb_vo"]("b,j");
        reused_["78_bb_vo"]("b,j") -= 0.50 * reused_["75_bb_ov"]("j,b");
        reused_["78_bb_vo"]("b,j") += reused_["49_bb_ov"]("j,b");
        reused_["78_bb_vo"]("b,j") -= reused_["76_bb_vo"]("b,j");
        reused_["78_bb_vo"]("b,j") += reused_["58_bb_vo"]("b,j");
        reused_["78_bb_vo"]("b,j") -= 0.50 * reused_["57_bb_vo"]("b,j");
        reused_["79_aa_ov"]("j,a")  = t2["abab"]("a,c,l,k") * eri["abba_oovo"]("l,k,c,j");
        reused_["80_aa_vo"]("b,j")  = f["aa_ov"]("k,c") * t2["aaaa"]("c,b,j,k");
        reused_["81_aa_vo"]("b,j")  = f["bb_ov"]("k,c") * t2["abab"]("b,c,j,k");
        reused_["82_aa_vo"]("b,j")  = reused_["81_aa_vo"]("b,j");
        reused_["82_aa_vo"]("b,j") += reused_["79_aa_ov"]("j,b");
        reused_["82_aa_vo"]("b,j") += reused_["48_aa_vo"]("b,j");
        reused_["82_aa_vo"]("b,j") -= reused_["80_aa_vo"]("b,j");
        reused_["82_aa_vo"]("b,j") += 0.50 * reused_["54_aa_ov"]("j,b");
        reused_["82_aa_vo"]("b,j") += 0.50 * reused_["64_aa_vo"]("b,j");
        reused_["83_aabb_oovo"]("l,i,b,j")  = reused_["65_aabb_oovo"]("l,i,b,j");
        reused_["83_aabb_oovo"]("l,i,b,j") += reused_["28_aabb_oovo"]("l,i,b,j");
        reused_["83_aabb_oovo"]("l,i,b,j") -= reused_["32_abba_oovo"]("l,j,b,i");
        reused_["84_bbaa_oovo"]("l,j,a,i")  = reused_["34_bbaa_oovo"]("l,j,a,i");
        reused_["84_bbaa_oovo"]("l,j,a,i") -= reused_["38_bbaa_oovo"]("l,j,a,i");
        reused_["84_bbaa_oovo"]("l,j,a,i") += reused_["31_baab_oovo"]("l,i,a,j");
        reused_["85_bbbb_oovo"]("l,j,a,i")  = reused_["35_bbbb_oovo"]("l,j,a,i");
        reused_["85_bbbb_oovo"]("l,j,a,i") -= reused_["67_bbbb_oovo"]("l,j,a,i");
        reused_["86_aaaa_oovo"]("l,j,a,i")  = reused_["37_aaaa_oovo"]("l,j,a,i");
        reused_["86_aaaa_oovo"]("l,j,a,i") += reused_["29_aaaa_oovo"]("l,j,a,i");
        reused_["87_bb_oo"]("i,j")  = eri["bbbb_oovv"]("i,k,b,c") * t2["bbbb"]("b,c,j,k");
        reused_["88_bb_oo"]("i,k")  = eri["bbbb_oovv"]("j,k,b,c") * t2["bbbb"]("b,c,i,j");
        reused_["89_aaaa_vooo"]("b,j,k,i")  = reused_["51_aaaa_vooo"]("b,j,k,i");
        reused_["89_aaaa_vooo"]("b,j,k,i") -= 0.50 * reused_["11_aaaa_vooo"]("b,i,j,k");
        reused_["90_baba_vooo"]("b,j,k,i")  = f["aa_ov"]("i,c") * t2["abab"]("c,b,j,k");
        reused_["91_baab_vooo"]("b,i,j,k")  = reused_["43_baab_vooo"]("b,i,j,k");
        reused_["91_baab_vooo"]("b,i,j,k") -= reused_["90_baba_vooo"]("b,j,k,i");
        reused_["92_aa_oo"]("i,k")  = eri["abab_oovv"]("k,j,b,c") * t2["abab"]("b,c,i,j");
        reused_["93_abab_vooo"]("b,i,j,k")  = reused_["9_abab_vooo"]("b,i,j,k");
        reused_["93_abab_vooo"]("b,i,j,k") += reused_["60_aabb_vooo"]("b,j,k,i");
    }

    world_.gop.fence();
}