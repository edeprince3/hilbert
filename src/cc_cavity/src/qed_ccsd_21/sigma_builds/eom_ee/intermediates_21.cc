
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

#include "cc_cavity/include/qed_ccsd_21/eom_ee_qed_ccsd_21.h"
#include "cc_cavity/include/qed_ccsd_21/qed_ccsd_21.h"

double* hilbert::EOM_EE_QED_CCSD::build_ss_diagonal() {

    if (!cc_wfn_->has_t1_integrals_) cc_wfn_->transform_integrals(true);

    // Get cavity information
    double w0 = cc_wfn_->cavity_frequency_[2];
    double coupling_factor_z = w0 * cc_wfn_->cavity_coupling_strength_[2];

    // get dipole integrals
    std::map<std::string, TA::TArrayD> dp;

    dp["aa_oo"]("i, j") = coupling_factor_z * cc_wfn_->Dip_blks_["dz_aa_oo"]("i, j");
    dp["aa_ov"]("i, a") = coupling_factor_z * cc_wfn_->Dip_blks_["dz_aa_ov"]("i, a");
    dp["aa_vo"]("a, i") = coupling_factor_z * cc_wfn_->Dip_blks_["dz_aa_vo"]("a, i");
    dp["aa_vv"]("a, b") = coupling_factor_z * cc_wfn_->Dip_blks_["dz_aa_vv"]("a, b");

    dp["bb_oo"]("i, j") = coupling_factor_z * cc_wfn_->Dip_blks_["dz_bb_oo"]("i, j");
    dp["bb_ov"]("i, a") = coupling_factor_z * cc_wfn_->Dip_blks_["dz_bb_ov"]("i, a");
    dp["bb_vo"]("a, i") = coupling_factor_z * cc_wfn_->Dip_blks_["dz_bb_vo"]("a, i");
    dp["bb_vv"]("a, b") = coupling_factor_z * cc_wfn_->Dip_blks_["dz_bb_vv"]("a, b");

    // extract 0-body amplitudes
    double t0_1;
    foreach_inplace(cc_wfn_->amplitudes_["t0_1"], [&t0_1](auto &tile){
        for(auto &x : tile.range())
            t0_1 = tile[x];
    });
    cc_wfn_->scalars_["t0_1"] = t0_1;
    scalars_["t0_1"] = t0_1;

    // extract 1-body amplitudes
    std::map<std::string, TA::TArrayD> t1 {
            {"aa", cc_wfn_->amplitudes_["t1_aa"]},
            {"bb", cc_wfn_->amplitudes_["t1_bb"]}
    };
    std::map<std::string, TA::TArrayD> t1_1 {
            {"aa", cc_wfn_->amplitudes_["t1_1_aa"]},
            {"bb", cc_wfn_->amplitudes_["t1_1_bb"]}
    };

    // extract 2-body amplitudes
    std::map<std::string, TA::TArrayD> t2 {
            {"aaaa", cc_wfn_->amplitudes_["t2_aaaa"]},
            {"abab", cc_wfn_->amplitudes_["t2_abab"]},
            {"bbbb", cc_wfn_->amplitudes_["t2_bbbb"]}
    };
    std::map<std::string, TA::TArrayD> t2_1 {
            {"aaaa", cc_wfn_->amplitudes_["t2_1_aaaa"]},
            {"abab", cc_wfn_->amplitudes_["t2_1_abab"]},
            {"bbbb", cc_wfn_->amplitudes_["t2_1_bbbb"]}
    };

    world_.gop.fence();

    TArrayMap &eri = cc_wfn_->V_blks_;
    TArrayMap &f = cc_wfn_->F_blks_;
    TArrayMap &Id = cc_wfn_->Id_blks_;

    // initialize diagonal blocks
    double  H_hw0 = 0;
    double cH_hw0 = 0;
    TArrayD    H_ss_aaaa = makeTensor(world_, {oa_,va_, oa_,va_}, true);
    TArrayD    H_ss_bbbb = makeTensor(world_, {ob_,vb_, ob_,vb_}, true);
    TArrayD   cH_ss_aaaa = makeTensor(world_, {oa_,va_, oa_,va_}, true);
    TArrayD   cH_ss_bbbb = makeTensor(world_, {ob_,vb_, ob_,vb_}, true);

    TArrayD  H_hw1_aaaa = makeTensor(world_, {oa_,va_, oa_,va_}, true);
    TArrayD  H_hw1_bbbb = makeTensor(world_, {ob_,vb_, ob_,vb_}, true);
    TArrayD cH_hw1_aaaa = makeTensor(world_, {oa_,va_, oa_,va_}, true);
    TArrayD cH_hw1_bbbb = makeTensor(world_, {ob_,vb_, ob_,vb_}, true);

    {
        scalars_["1"]  = 1.00 * dot(eri["aaaa_oovv"]("i,j,a,b"), t2["aaaa"]("a,b,j,i"));
        scalars_["2"]  = 1.00 * dot(eri["abab_oovv"]("j,i,a,b"), t2["abab"]("a,b,j,i"));
        scalars_["3"]  = 1.00 * dot(eri["bbbb_oovv"]("i,j,a,b"), t2["bbbb"]("a,b,j,i"));
        scalars_["4"]  = 1.00 * dot(Id["bb_oo"]("i,i1") * eri["bbbb_oooo"]("i,j,i1,j1"), Id["bb_oo"]("j,j1"));
        scalars_["5"]  = 1.00 * dot(Id["bb_oo"]("i,i1") * eri["abab_oooo"]("j,i,j1,i1"), Id["aa_oo"]("j,j1"));
        scalars_["6"]  = 1.00 * dot(Id["aa_oo"]("i,i1") * eri["aaaa_oooo"]("i,j,i1,j1"), Id["aa_oo"]("j,j1"));
        scalars_["7"]  = 1.00 * dot(dp["aa_ov"]("i,a"), t1_1["aa"]("a,i"));
        scalars_["8"]  = 1.00 * dot(dp["bb_ov"]("i,a"), t1_1["bb"]("a,i"));
        scalars_["9"]  = 1.00 * dot(Id["bb_oo"]("i,i1"), dp["bb_oo"]("i,i1")) * t0_1;
        scalars_["10"]  = 1.00 * dot(Id["bb_oo"]("i,i1"), f["bb_oo"]("i,i1"));
        scalars_["11"]  = 1.00 * dot(Id["aa_oo"]("i,i1"), dp["aa_oo"]("i,i1")) * t0_1;
        scalars_["12"]  = 1.00 * dot(Id["aa_oo"]("i,i1"), f["aa_oo"]("i,i1"));
        scalars_["13"]  = scalars_["12"];
        scalars_["13"] -= 0.25 * scalars_["1"];
        scalars_["13"] -= scalars_["5"];
        scalars_["13"] -= 2.00 * scalars_["7"];
        scalars_["13"] += w0;
        scalars_["13"] -= 0.25 * scalars_["3"];
        scalars_["13"] -= 0.50 * scalars_["4"];
        scalars_["13"] -= 2.00 * scalars_["8"];
        scalars_["13"] -= scalars_["9"];
        scalars_["13"] -= scalars_["11"];
        scalars_["13"] += scalars_["2"];
        scalars_["13"] -= 0.50 * scalars_["6"];
        scalars_["13"] += scalars_["10"];
        scalars_["13"]  = scalars_["12"];
        scalars_["13"] -= 0.25 * scalars_["1"];
        scalars_["13"] -= scalars_["5"];
        scalars_["13"] -= 2.00 * scalars_["7"];
        scalars_["13"] += w0;
        scalars_["13"] -= 0.25 * scalars_["3"];
        scalars_["13"] -= 0.50 * scalars_["4"];
        scalars_["13"] -= 2.00 * scalars_["8"];
        scalars_["13"] -= scalars_["9"];
        scalars_["13"] -= scalars_["11"];
        scalars_["13"] += scalars_["2"];
        scalars_["13"] -= 0.50 * scalars_["6"];
        scalars_["13"] += scalars_["10"];
        scalars_["14"]  = 1.00 * dot(Id["aa_oo"]("j,j1"), dp["aa_oo"]("j,j1"));
        scalars_["15"]  = 1.00 * dot(Id["bb_oo"]("j,j1"), dp["bb_oo"]("j,j1"));
        scalars_["16"]  = scalars_["12"];
        scalars_["16"] -= scalars_["5"];
        scalars_["16"] -= 0.50 * scalars_["4"];
        scalars_["16"] += scalars_["10"];
        scalars_["16"] += scalars_["2"];
        scalars_["16"] -= 0.50 * scalars_["6"];
        scalars_["16"] -= 0.25 * scalars_["3"];
        scalars_["16"] -= scalars_["9"];
        scalars_["16"] -= scalars_["11"];
        scalars_["16"] -= 0.25 * scalars_["1"];
        scalars_["17"]  = scalars_["8"];
        scalars_["17"] += scalars_["7"];

        // H_hw0  = +1.00 f_aa(i,i)
        //       += +0.25 <j,i||a,b>_aaaa t2_aaaa(a,b,j,i)
        //       += -0.50 <j,i||j,i>_abab
        //       += -0.50 <i,j||i,j>_abab
        //       += -2.00 d-_aa(i,a) t1_1_aa(a,i)
        //       += +1.00 w0
        //       += +0.25 <j,i||a,b>_bbbb t2_bbbb(a,b,j,i)
        //       += -0.50 <j,i||j,i>_bbbb
        //       += -2.00 d-_bb(i,a) t1_1_bb(a,i)
        //       += -1.00 d-_bb(i,i) t0_1
        //       += -1.00 d-_aa(i,i) t0_1
        //       += +0.25 <j,i||a,b>_abab t2_abab(a,b,j,i)
        //       += +0.25 <i,j||a,b>_abab t2_abab(a,b,i,j)
        //       += +0.25 <j,i||b,a>_abab t2_abab(b,a,j,i)
        //       += +0.25 <i,j||b,a>_abab t2_abab(b,a,i,j)
        //       += -0.50 <j,i||j,i>_aaaa
        //       += +1.00 f_bb(i,i)
        H_hw0  = scalars_["13"];

        // cH_hw0  = +1.00 t0_1
        cH_hw0  = t0_1;

        // flops: o2v0  = o2v0
        //  mems: o2v0  = o2v0
        tmps_["3_aa_oo"]("m,i")  = scalars_["16"] * Id["aa_oo"]("m,i");

        // flops: o2v0  = o2v0
        //  mems: o2v0  = o2v0
        tmps_["2_aa_oo"]("i,m")  = t0_1 * dp["aa_oo"]("i,m");

        // flops: o2v0  = o2v0
        //  mems: o2v0  = o2v0
        tmps_["1_aa_oo"]("m,i")  = t0_1 * Id["aa_oo"]("m,i");

        // flops: o2v2  = o2v2 o3v3 o2v2 o2v2 o2v2 o2v2 o3v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o3v2 o2v2 o2v2 o3v3 o2v2 o2v2 o2v2 o2v3 o2v2 o2v2 o2v3 o2v2 o2v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v0 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v0 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o0v2 o2v2 o2v2 o0v2 o2v2 o2v2 o2v2
        tmps_["4_aaaa_vvoo"]("e,a,m,i")  = dp["aa_vv"]("e,a") * tmps_["1_aa_oo"]("m,i");
        tmps_["4_aaaa_vvoo"]("e,a,m,i") += eri["aaaa_oovv"]("i,j,a,b") * t2["aaaa"]("b,e,m,j");
        tmps_["4_aaaa_vvoo"]("e,a,m,i") -= tmps_["2_aa_oo"]("i,m") * Id["aa_vv"]("e,a");
        tmps_["4_aaaa_vvoo"]("e,a,m,i") += Id["aa_vv"]("e,a") * f["aa_oo"]("i,m");
        tmps_["4_aaaa_vvoo"]("e,a,m,i") += eri["abab_oovv"]("i,j,b,c") * t2["abab"]("b,c,m,j") * Id["aa_vv"]("e,a");
        tmps_["4_aaaa_vvoo"]("e,a,m,i") -= Id["aa_oo"]("m,i") * f["aa_vv"]("e,a");
        tmps_["4_aaaa_vvoo"]("e,a,m,i") += dp["aa_ov"]("i,a") * t1_1["aa"]("e,m");
        tmps_["4_aaaa_vvoo"]("e,a,m,i") += eri["aaaa_vovo"]("e,i,a,m");
        tmps_["4_aaaa_vvoo"]("e,a,m,i") += 0.50 * eri["aaaa_oovv"]("i,j,b,c") * t2["aaaa"]("b,c,m,j") * Id["aa_vv"]("e,a");
        tmps_["4_aaaa_vvoo"]("e,a,m,i") -= eri["abab_oovv"]("i,j,a,b") * t2["abab"]("e,b,m,j");
        tmps_["4_aaaa_vvoo"]("e,a,m,i") -= tmps_["3_aa_oo"]("m,i") * Id["aa_vv"]("e,a");
        tmps_["4_aaaa_vvoo"]("e,a,m,i") += 0.50 * eri["aaaa_oovv"]("j,k,a,b") * t2["aaaa"]("b,e,k,j") * Id["aa_oo"]("m,i");
        tmps_["4_aaaa_vvoo"]("e,a,m,i") += eri["abab_oovv"]("k,j,a,b") * t2["abab"]("e,b,k,j") * Id["aa_oo"]("m,i");
        tmps_["3_aa_oo"].~TArrayD();
        tmps_["2_aa_oo"].~TArrayD();

        // H_hw1_aaaa  = -1.00 d_aa(m,i) d-_aa(e,a) t0_1
        //            += +1.00 <i,j||b,a>_aaaa t2_aaaa(b,e,m,j)
        //            += +1.00 d_aa(e,a) d-_aa(i,m) t0_1
        //            += -1.00 d_aa(e,a) f_aa(i,m)
        //            += -0.50 d_aa(e,a) <i,j||b,c>_abab t2_abab(b,c,m,j)
        //            += -0.50 d_aa(e,a) <i,j||c,b>_abab t2_abab(c,b,m,j)
        //            += +1.00 d_aa(m,i) f_aa(e,a)
        //            += -1.00 d-_aa(i,a) t1_1_aa(e,m)
        //            += +1.00 <i,e||a,m>_aaaa
        //            += -0.50 d_aa(e,a) <i,j||b,c>_aaaa t2_aaaa(b,c,m,j)
        //            += +1.00 <i,j||a,b>_abab t2_abab(e,b,m,j)
        //            += +1.00 d_aa(e,a) d_aa(m,i) f_aa(j,j)
        //            += -0.50 d_aa(e,a) d_aa(m,i) <k,j||k,j>_abab
        //            += -0.50 d_aa(e,a) d_aa(m,i) <j,k||j,k>_abab
        //            += -0.50 d_aa(e,a) d_aa(m,i) <k,j||k,j>_bbbb
        //            += +1.00 d_aa(e,a) d_aa(m,i) f_bb(j,j)
        //            += +0.25 d_aa(e,a) d_aa(m,i) <k,j||b,c>_abab t2_abab(b,c,k,j)
        //            += +0.25 d_aa(e,a) d_aa(m,i) <j,k||b,c>_abab t2_abab(b,c,j,k)
        //            += +0.25 d_aa(e,a) d_aa(m,i) <k,j||c,b>_abab t2_abab(c,b,k,j)
        //            += +0.25 d_aa(e,a) d_aa(m,i) <j,k||c,b>_abab t2_abab(c,b,j,k)
        //            += -0.50 d_aa(e,a) d_aa(m,i) <k,j||k,j>_aaaa
        //            += +0.25 d_aa(e,a) d_aa(m,i) <k,j||b,c>_bbbb t2_bbbb(b,c,k,j)
        //            += -1.00 d_aa(e,a) d_aa(m,i) d-_bb(j,j) t0_1
        //            += -1.00 d_aa(e,a) d_aa(m,i) d-_aa(j,j) t0_1
        //            += +0.25 d_aa(e,a) d_aa(m,i) <k,j||b,c>_aaaa t2_aaaa(b,c,k,j)
        //            += -0.50 d_aa(m,i) <k,j||b,a>_aaaa t2_aaaa(b,e,k,j)
        //            += -0.50 d_aa(m,i) <k,j||a,b>_abab t2_abab(e,b,k,j)
        //            += -0.50 d_aa(m,i) <j,k||a,b>_abab t2_abab(e,b,j,k)
        H_hw1_aaaa("a,e,i,m")  = -1.00 * tmps_["4_aaaa_vvoo"]("e,a,m,i");

        // H_ss_aaaa  = -1.00 d_aa(m,i) d-_aa(e,a) t0_1
        //           += +1.00 <i,j||b,a>_aaaa t2_aaaa(b,e,m,j)
        //           += +1.00 d_aa(e,a) d-_aa(i,m) t0_1
        //           += -1.00 d_aa(e,a) f_aa(i,m)
        //           += -0.50 d_aa(e,a) <i,j||b,c>_abab t2_abab(b,c,m,j)
        //           += -0.50 d_aa(e,a) <i,j||c,b>_abab t2_abab(c,b,m,j)
        //           += +1.00 d_aa(m,i) f_aa(e,a)
        //           += -1.00 d-_aa(i,a) t1_1_aa(e,m)
        //           += +1.00 <i,e||a,m>_aaaa
        //           += -0.50 d_aa(e,a) <i,j||b,c>_aaaa t2_aaaa(b,c,m,j)
        //           += +1.00 <i,j||a,b>_abab t2_abab(e,b,m,j)
        //           += +1.00 d_aa(e,a) d_aa(m,i) f_aa(j,j)
        //           += -0.50 d_aa(e,a) d_aa(m,i) <k,j||k,j>_abab
        //           += -0.50 d_aa(e,a) d_aa(m,i) <j,k||j,k>_abab
        //           += -0.50 d_aa(e,a) d_aa(m,i) <k,j||k,j>_bbbb
        //           += +1.00 d_aa(e,a) d_aa(m,i) f_bb(j,j)
        //           += +0.25 d_aa(e,a) d_aa(m,i) <k,j||b,c>_abab t2_abab(b,c,k,j)
        //           += +0.25 d_aa(e,a) d_aa(m,i) <j,k||b,c>_abab t2_abab(b,c,j,k)
        //           += +0.25 d_aa(e,a) d_aa(m,i) <k,j||c,b>_abab t2_abab(c,b,k,j)
        //           += +0.25 d_aa(e,a) d_aa(m,i) <j,k||c,b>_abab t2_abab(c,b,j,k)
        //           += -0.50 d_aa(e,a) d_aa(m,i) <k,j||k,j>_aaaa
        //           += +0.25 d_aa(e,a) d_aa(m,i) <k,j||b,c>_bbbb t2_bbbb(b,c,k,j)
        //           += -1.00 d_aa(e,a) d_aa(m,i) d-_bb(j,j) t0_1
        //           += -1.00 d_aa(e,a) d_aa(m,i) d-_aa(j,j) t0_1
        //           += +0.25 d_aa(e,a) d_aa(m,i) <k,j||b,c>_aaaa t2_aaaa(b,c,k,j)
        //           += -0.50 d_aa(m,i) <k,j||b,a>_aaaa t2_aaaa(b,e,k,j)
        //           += -0.50 d_aa(m,i) <k,j||a,b>_abab t2_abab(e,b,k,j)
        //           += -0.50 d_aa(m,i) <j,k||a,b>_abab t2_abab(e,b,j,k)
        H_ss_aaaa("a,e,i,m")  = -1.00 * tmps_["4_aaaa_vvoo"]("e,a,m,i");
        tmps_["4_aaaa_vvoo"].~TArrayD();

        // flops: o2v0  = o2v0
        //  mems: o2v0  = o2v0
        tmps_["5_aa_oo"]("m,i")  = w0 * Id["aa_oo"]("m,i");

        // flops: o2v2  = o2v2
        //  mems: o2v2  = o2v2
        tmps_["6_aaaa_vvoo"]("e,a,m,i")  = Id["aa_vv"]("e,a") * tmps_["5_aa_oo"]("m,i");
        tmps_["5_aa_oo"].~TArrayD();

        // H_hw1_aaaa += +1.00 d_aa(e,a) d_aa(m,i) w0
        H_hw1_aaaa("a,e,i,m") += tmps_["6_aaaa_vvoo"]("e,a,m,i");
        tmps_["6_aaaa_vvoo"].~TArrayD();

        // flops: o0v2  = o1v2
        //  mems: o0v2  = o0v2
        tmps_["9_aa_vv"]("a,e")  = dp["aa_ov"]("j,a") * t1_1["aa"]("e,j");

        // flops: o2v0  = o2v1
        //  mems: o2v0  = o2v0
        tmps_["8_aa_oo"]("i,m")  = dp["aa_ov"]("i,b") * t1_1["aa"]("b,m");

        // flops: o2v0  = o2v0
        //  mems: o2v0  = o2v0
        tmps_["7_aa_oo"]("m,i")  = scalars_["17"] * Id["aa_oo"]("m,i");

        // flops: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["10_aaaa_vvoo"]("a,e,m,i")  = Id["aa_oo"]("m,i") * tmps_["9_aa_vv"]("a,e");
        tmps_["10_aaaa_vvoo"]("a,e,m,i") -= Id["aa_vv"]("e,a") * tmps_["7_aa_oo"]("m,i");
        tmps_["10_aaaa_vvoo"]("a,e,m,i") += Id["aa_vv"]("e,a") * tmps_["8_aa_oo"]("i,m");
        tmps_["9_aa_vv"].~TArrayD();
        tmps_["8_aa_oo"].~TArrayD();
        tmps_["7_aa_oo"].~TArrayD();

        // H_hw1_aaaa += +2.00 d_aa(m,i) d-_aa(j,a) t1_1_aa(e,j)
        //            += -2.00 d_aa(e,a) d_aa(m,i) d-_bb(j,b) t1_1_bb(b,j)
        //            += -2.00 d_aa(e,a) d_aa(m,i) d-_aa(j,b) t1_1_aa(b,j)
        //            += +2.00 d_aa(e,a) d-_aa(i,b) t1_1_aa(b,m)
        H_hw1_aaaa("a,e,i,m") += 2.00 * tmps_["10_aaaa_vvoo"]("a,e,m,i");

        // H_ss_aaaa += +1.00 d_aa(m,i) d-_aa(j,a) t1_1_aa(e,j)
        //           += -1.00 d_aa(e,a) d_aa(m,i) d-_bb(j,b) t1_1_bb(b,j)
        //           += -1.00 d_aa(e,a) d_aa(m,i) d-_aa(j,b) t1_1_aa(b,j)
        //           += +1.00 d_aa(e,a) d-_aa(i,b) t1_1_aa(b,m)
        H_ss_aaaa("a,e,i,m") += tmps_["10_aaaa_vvoo"]("a,e,m,i");
        tmps_["10_aaaa_vvoo"].~TArrayD();
        tmps_["14_bbbb_vovo"]("e,i,a,m")  = 2.00 * eri["bbbb_vovo"]("e,i,a,m");

        // flops: o2v0  = o2v0
        //  mems: o2v0  = o2v0
        tmps_["13_bb_oo"]("i,m")  = t0_1 * dp["bb_oo"]("i,m");

        // flops: o2v0  = o2v0
        //  mems: o2v0  = o2v0
        tmps_["12_bb_oo"]("m,i")  = scalars_["16"] * Id["bb_oo"]("m,i");

        // flops: o2v0  = o2v0
        //  mems: o2v0  = o2v0
        tmps_["11_bb_oo"]("m,i")  = t0_1 * Id["bb_oo"]("m,i");

        // flops: o2v2  = o2v2 o2v3 o2v2 o3v2 o2v2 o2v2 o3v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o3v3 o2v2 o2v2 o2v2 o2v2 o2v2 o3v3 o2v2 o2v2 o2v2 o2v3 o2v2 o2v2 o2v2 o2v2
        //  mems: o2v2  = o2v2 o0v2 o2v2 o2v0 o2v2 o2v2 o2v0 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o2v2 o0v2 o2v2 o2v2 o2v2 o2v2
        tmps_["15_bbbb_vvoo"]("e,a,m,i")  = dp["bb_vv"]("e,a") * tmps_["11_bb_oo"]("m,i");
        tmps_["15_bbbb_vvoo"]("e,a,m,i") += 0.50 * eri["bbbb_oovv"]("j,k,a,b") * t2["bbbb"]("b,e,k,j") * Id["bb_oo"]("m,i");
        tmps_["15_bbbb_vvoo"]("e,a,m,i") += eri["abab_oovv"]("j,i,b,c") * t2["abab"]("b,c,j,m") * Id["bb_vv"]("e,a");
        tmps_["15_bbbb_vvoo"]("e,a,m,i") += 0.50 * eri["bbbb_oovv"]("i,j,b,c") * t2["bbbb"]("b,c,m,j") * Id["bb_vv"]("e,a");
        tmps_["15_bbbb_vvoo"]("e,a,m,i") -= tmps_["12_bb_oo"]("m,i") * Id["bb_vv"]("e,a");
        tmps_["15_bbbb_vvoo"]("e,a,m,i") -= tmps_["13_bb_oo"]("i,m") * Id["bb_vv"]("e,a");
        tmps_["15_bbbb_vvoo"]("e,a,m,i") -= eri["abab_oovv"]("j,i,b,a") * t2["abab"]("b,e,j,m");
        tmps_["15_bbbb_vvoo"]("e,a,m,i") += Id["bb_vv"]("e,a") * f["bb_oo"]("i,m");
        tmps_["15_bbbb_vvoo"]("e,a,m,i") += dp["bb_ov"]("i,a") * t1_1["bb"]("e,m");
        tmps_["15_bbbb_vvoo"]("e,a,m,i") += eri["bbbb_oovv"]("i,j,a,b") * t2["bbbb"]("b,e,m,j");
        tmps_["15_bbbb_vvoo"]("e,a,m,i") -= Id["bb_oo"]("m,i") * f["bb_vv"]("e,a");
        tmps_["15_bbbb_vvoo"]("e,a,m,i") += eri["abab_oovv"]("k,j,b,a") * t2["abab"]("b,e,k,j") * Id["bb_oo"]("m,i");
        tmps_["15_bbbb_vvoo"]("e,a,m,i") += 0.50 * tmps_["14_bbbb_vovo"]("e,i,a,m");
        tmps_["14_bbbb_vovo"].~TArrayD();
        tmps_["13_bb_oo"].~TArrayD();
        tmps_["12_bb_oo"].~TArrayD();

        // H_hw1_bbbb  = -1.00 d_bb(m,i) d-_bb(e,a) t0_1
        //            += -0.50 d_bb(m,i) <k,j||b,a>_bbbb t2_bbbb(b,e,k,j)
        //            += -0.50 d_bb(e,a) <j,i||b,c>_abab t2_abab(b,c,j,m)
        //            += -0.50 d_bb(e,a) <j,i||c,b>_abab t2_abab(c,b,j,m)
        //            += -0.50 d_bb(e,a) <i,j||b,c>_bbbb t2_bbbb(b,c,m,j)
        //            += +1.00 d_bb(e,a) d_bb(m,i) f_aa(j,j)
        //            += -0.50 d_bb(e,a) d_bb(m,i) <k,j||k,j>_abab
        //            += -0.50 d_bb(e,a) d_bb(m,i) <j,k||j,k>_abab
        //            += -0.50 d_bb(e,a) d_bb(m,i) <k,j||k,j>_bbbb
        //            += +1.00 d_bb(e,a) d_bb(m,i) f_bb(j,j)
        //            += +0.25 d_bb(e,a) d_bb(m,i) <k,j||b,c>_abab t2_abab(b,c,k,j)
        //            += +0.25 d_bb(e,a) d_bb(m,i) <j,k||b,c>_abab t2_abab(b,c,j,k)
        //            += +0.25 d_bb(e,a) d_bb(m,i) <k,j||c,b>_abab t2_abab(c,b,k,j)
        //            += +0.25 d_bb(e,a) d_bb(m,i) <j,k||c,b>_abab t2_abab(c,b,j,k)
        //            += -0.50 d_bb(e,a) d_bb(m,i) <k,j||k,j>_aaaa
        //            += +0.25 d_bb(e,a) d_bb(m,i) <k,j||b,c>_bbbb t2_bbbb(b,c,k,j)
        //            += -1.00 d_bb(e,a) d_bb(m,i) d-_bb(j,j) t0_1
        //            += -1.00 d_bb(e,a) d_bb(m,i) d-_aa(j,j) t0_1
        //            += +0.25 d_bb(e,a) d_bb(m,i) <k,j||b,c>_aaaa t2_aaaa(b,c,k,j)
        //            += +1.00 d_bb(e,a) d-_bb(i,m) t0_1
        //            += +1.00 <j,i||b,a>_abab t2_abab(b,e,j,m)
        //            += -1.00 d_bb(e,a) f_bb(i,m)
        //            += -1.00 d-_bb(i,a) t1_1_bb(e,m)
        //            += +1.00 <i,j||b,a>_bbbb t2_bbbb(b,e,m,j)
        //            += +1.00 d_bb(m,i) f_bb(e,a)
        //            += -0.50 d_bb(m,i) <k,j||b,a>_abab t2_abab(b,e,k,j)
        //            += -0.50 d_bb(m,i) <j,k||b,a>_abab t2_abab(b,e,j,k)
        //            += +1.00 <i,e||a,m>_bbbb
        H_hw1_bbbb("a,e,i,m")  = -1.00 * tmps_["15_bbbb_vvoo"]("e,a,m,i");

        // H_ss_bbbb  = -1.00 d_bb(m,i) d-_bb(e,a) t0_1
        //           += -0.50 d_bb(m,i) <k,j||b,a>_bbbb t2_bbbb(b,e,k,j)
        //           += -0.50 d_bb(e,a) <j,i||b,c>_abab t2_abab(b,c,j,m)
        //           += -0.50 d_bb(e,a) <j,i||c,b>_abab t2_abab(c,b,j,m)
        //           += -0.50 d_bb(e,a) <i,j||b,c>_bbbb t2_bbbb(b,c,m,j)
        //           += +1.00 d_bb(e,a) d_bb(m,i) f_aa(j,j)
        //           += -0.50 d_bb(e,a) d_bb(m,i) <k,j||k,j>_abab
        //           += -0.50 d_bb(e,a) d_bb(m,i) <j,k||j,k>_abab
        //           += -0.50 d_bb(e,a) d_bb(m,i) <k,j||k,j>_bbbb
        //           += +1.00 d_bb(e,a) d_bb(m,i) f_bb(j,j)
        //           += +0.25 d_bb(e,a) d_bb(m,i) <k,j||b,c>_abab t2_abab(b,c,k,j)
        //           += +0.25 d_bb(e,a) d_bb(m,i) <j,k||b,c>_abab t2_abab(b,c,j,k)
        //           += +0.25 d_bb(e,a) d_bb(m,i) <k,j||c,b>_abab t2_abab(c,b,k,j)
        //           += +0.25 d_bb(e,a) d_bb(m,i) <j,k||c,b>_abab t2_abab(c,b,j,k)
        //           += -0.50 d_bb(e,a) d_bb(m,i) <k,j||k,j>_aaaa
        //           += +0.25 d_bb(e,a) d_bb(m,i) <k,j||b,c>_bbbb t2_bbbb(b,c,k,j)
        //           += -1.00 d_bb(e,a) d_bb(m,i) d-_bb(j,j) t0_1
        //           += -1.00 d_bb(e,a) d_bb(m,i) d-_aa(j,j) t0_1
        //           += +0.25 d_bb(e,a) d_bb(m,i) <k,j||b,c>_aaaa t2_aaaa(b,c,k,j)
        //           += +1.00 d_bb(e,a) d-_bb(i,m) t0_1
        //           += +1.00 <j,i||b,a>_abab t2_abab(b,e,j,m)
        //           += -1.00 d_bb(e,a) f_bb(i,m)
        //           += -1.00 d-_bb(i,a) t1_1_bb(e,m)
        //           += +1.00 <i,j||b,a>_bbbb t2_bbbb(b,e,m,j)
        //           += +1.00 d_bb(m,i) f_bb(e,a)
        //           += -0.50 d_bb(m,i) <k,j||b,a>_abab t2_abab(b,e,k,j)
        //           += -0.50 d_bb(m,i) <j,k||b,a>_abab t2_abab(b,e,j,k)
        //           += +1.00 <i,e||a,m>_bbbb
        H_ss_bbbb("a,e,i,m")  = -1.00 * tmps_["15_bbbb_vvoo"]("e,a,m,i");
        tmps_["15_bbbb_vvoo"].~TArrayD();

        // flops: o2v0  = o2v0
        //  mems: o2v0  = o2v0
        tmps_["16_bb_oo"]("m,i")  = w0 * Id["bb_oo"]("m,i");

        // flops: o2v2  = o2v2
        //  mems: o2v2  = o2v2
        tmps_["17_bbbb_vvoo"]("e,a,m,i")  = Id["bb_vv"]("e,a") * tmps_["16_bb_oo"]("m,i");
        tmps_["16_bb_oo"].~TArrayD();

        // H_hw1_bbbb += +1.00 d_bb(e,a) d_bb(m,i) w0
        H_hw1_bbbb("a,e,i,m") += tmps_["17_bbbb_vvoo"]("e,a,m,i");
        tmps_["17_bbbb_vvoo"].~TArrayD();

        // flops: o0v2  = o1v2
        //  mems: o0v2  = o0v2
        tmps_["20_bb_vv"]("a,e")  = dp["bb_ov"]("j,a") * t1_1["bb"]("e,j");

        // flops: o2v0  = o2v1
        //  mems: o2v0  = o2v0
        tmps_["19_bb_oo"]("i,m")  = dp["bb_ov"]("i,b") * t1_1["bb"]("b,m");

        // flops: o2v0  = o2v0
        //  mems: o2v0  = o2v0
        tmps_["18_bb_oo"]("m,i")  = scalars_["17"] * Id["bb_oo"]("m,i");

        // flops: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2
        //  mems: o2v2  = o2v2 o2v2 o2v2 o2v2 o2v2
        tmps_["21_bbbb_vvoo"]("a,e,m,i")  = Id["bb_oo"]("m,i") * tmps_["20_bb_vv"]("a,e");
        tmps_["21_bbbb_vvoo"]("a,e,m,i") += Id["bb_vv"]("e,a") * tmps_["19_bb_oo"]("i,m");
        tmps_["21_bbbb_vvoo"]("a,e,m,i") -= Id["bb_vv"]("e,a") * tmps_["18_bb_oo"]("m,i");
        tmps_["20_bb_vv"].~TArrayD();
        tmps_["19_bb_oo"].~TArrayD();
        tmps_["18_bb_oo"].~TArrayD();

        // H_hw1_bbbb += +2.00 d_bb(m,i) d-_bb(j,a) t1_1_bb(e,j)
        //            += +2.00 d_bb(e,a) d-_bb(i,b) t1_1_bb(b,m)
        //            += -2.00 d_bb(e,a) d_bb(m,i) d-_bb(j,b) t1_1_bb(b,j)
        //            += -2.00 d_bb(e,a) d_bb(m,i) d-_aa(j,b) t1_1_aa(b,j)
        H_hw1_bbbb("a,e,i,m") += 2.00 * tmps_["21_bbbb_vvoo"]("a,e,m,i");

        // H_ss_bbbb += +1.00 d_bb(m,i) d-_bb(j,a) t1_1_bb(e,j)
        //           += +1.00 d_bb(e,a) d-_bb(i,b) t1_1_bb(b,m)
        //           += -1.00 d_bb(e,a) d_bb(m,i) d-_bb(j,b) t1_1_bb(b,j)
        //           += -1.00 d_bb(e,a) d_bb(m,i) d-_aa(j,b) t1_1_aa(b,j)
        H_ss_bbbb("a,e,i,m") += tmps_["21_bbbb_vvoo"]("a,e,m,i");
        tmps_["21_bbbb_vvoo"].~TArrayD();

        // flops: o2v2  = o2v2
        //  mems: o2v2  = o2v2
        tmps_["22_aaaa_vvoo"]("e,a,m,i")  = Id["aa_vv"]("e,a") * tmps_["1_aa_oo"]("m,i");
        tmps_["1_aa_oo"].~TArrayD();

        // cH_hw1_aaaa  = +1.00 d_aa(e,a) d_aa(m,i) t0_1
        cH_hw1_aaaa("a,e,i,m")  = tmps_["22_aaaa_vvoo"]("e,a,m,i");

        // cH_ss_aaaa  = +1.00 d_aa(e,a) d_aa(m,i) t0_1
        cH_ss_aaaa("a,e,i,m")  = tmps_["22_aaaa_vvoo"]("e,a,m,i");
        tmps_["22_aaaa_vvoo"].~TArrayD();

        // flops: o2v2  = o2v2
        //  mems: o2v2  = o2v2
        tmps_["23_bbbb_vvoo"]("e,a,m,i")  = Id["bb_vv"]("e,a") * tmps_["11_bb_oo"]("m,i");
        tmps_["11_bb_oo"].~TArrayD();

        // cH_hw1_bbbb  = +1.00 d_bb(e,a) d_bb(m,i) t0_1
        cH_hw1_bbbb("a,e,i,m")  = tmps_["23_bbbb_vvoo"]("e,a,m,i");

        // cH_ss_bbbb  = +1.00 d_bb(e,a) d_bb(m,i) t0_1
        cH_ss_bbbb("a,e,i,m")  = tmps_["23_bbbb_vvoo"]("e,a,m,i");
        tmps_["23_bbbb_vvoo"].~TArrayD();

    }

    // add coherent state basis terms
    double coherent_scalar;
    double e_dip_z_ = cc_wfn_->e_dip_z_;
    double nuc_dip_z_ = cc_wfn_->nuc_dip_z_;
    if ( options_.get_bool("QED_USE_RELAXED_ORBITALS")) {
        coherent_scalar = coupling_factor_z * e_dip_z_;
    } else {
        coherent_scalar = -coupling_factor_z * nuc_dip_z_;
    }

    H_hw0 += coherent_scalar * cH_hw0;
    H_ss_aaaa("m,e,i,a") += coherent_scalar * cH_ss_aaaa("m,e,i,a");
    H_ss_bbbb("m,e,i,a") += coherent_scalar * cH_ss_bbbb("m,e,i,a");
    H_hw1_aaaa("m,e,i,a") += coherent_scalar * cH_hw1_aaaa("m,e,i,a");
    H_hw1_bbbb("m,e,i,a") += coherent_scalar * cH_hw1_bbbb("m,e,i,a");

    // pluck out diagonals
    size_t dim_tot = singleDim_;
    if (includes_["t0_1"]) ++dim_tot;
    if (includes_["t1_1"]) dim_tot += singleDim_;
    auto * ss_diag = (double*) calloc(dim_tot, sizeof(double));

    size_t id = 0;
    ss_diag[id++] = H_hw0;

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

    // extract singles w/ hw
    foreach_inplace(H_hw1_aaaa, [id,oa,va,ob,vb,ss_diag](auto &tile){
        for (auto &x : tile.range()) {
            size_t m = x[0], e = x[1], i = x[2], a = x[3];
            if (m == i && e == a) {
                size_t index = e * oa + i + id;
                ss_diag[index] = tile[x];
            }
        }
    }); id += oa_*va_;

    foreach_inplace(H_hw1_bbbb, [id,oa,va,ob,vb,ss_diag](auto &tile){
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

void hilbert::EOM_EE_QED_CCSD::build_common_ops() {
    
    if (!cc_wfn_->has_t1_integrals_) cc_wfn_->transform_integrals(true);

    // Get cavity information
    double w0 = cc_wfn_->cavity_frequency_[2];
    double coupling_factor_z = w0 * cc_wfn_->cavity_coupling_strength_[2];

    // get dipole integrals
    TArrayMap dp = reinterpret_pointer_cast<QED_CCSD>(cc_wfn_)->effective_dipole();

    // extract 0-body amplitudes
    // extract 0-body amplitudes
    double t0_1;
    foreach_inplace(cc_wfn_->amplitudes_["t0_1"], [&t0_1](auto &tile){
        for(auto &x : tile.range())
            t0_1 = tile[x];
    });
    cc_wfn_->scalars_["t0_1"] = t0_1;
    scalars_["t0_1"] = t0_1;

    // extract 1-body amplitudes
    std::map<std::string, TA::TArrayD> t1 {
            {"aa", cc_wfn_->amplitudes_["t1_aa"]},
            {"bb", cc_wfn_->amplitudes_["t1_bb"]}
    };
    std::map<std::string, TA::TArrayD> t1_1 {
            {"aa", cc_wfn_->amplitudes_["t1_1_aa"]},
            {"bb", cc_wfn_->amplitudes_["t1_1_bb"]}
    };

    // extract 2-body amplitudes
    std::map<std::string, TA::TArrayD> t2 {
            {"aaaa", cc_wfn_->amplitudes_["t2_aaaa"]},
            {"abab", cc_wfn_->amplitudes_["t2_abab"]},
            {"bbbb", cc_wfn_->amplitudes_["t2_bbbb"]}
    };
    std::map<std::string, TA::TArrayD> t2_1 {
            {"aaaa", cc_wfn_->amplitudes_["t2_1_aaaa"]},
            {"abab", cc_wfn_->amplitudes_["t2_1_abab"]},
            {"bbbb", cc_wfn_->amplitudes_["t2_1_bbbb"]}
    };

    world_.gop.fence();

    TArrayMap &eri = cc_wfn_->V_blks_;
    TArrayMap &f = cc_wfn_->F_blks_;
    TArrayMap &Id = cc_wfn_->Id_blks_;

    {

        scalars_["1"]  = 1.00 * dot(Id["aa_oo"]("i,i1"), dp["aa_oo"]("i,i1"));
        scalars_["2"]  = 1.00 * dot(Id["bb_oo"]("k,k1"), dp["bb_oo"]("k,k1"));
        scalars_["3"]  = 1.00 * dot(dp["bb_ov"]("i,a"), t1_1["bb"]("a,i"));
        scalars_["4"]  = 1.00 * dot(dp["aa_ov"]("i,a"), t1_1["aa"]("a,i"));
        scalars_["5"]  = -1.00 * dot(Id["aa_oo"]("o,I"), f["aa_oo"]("o,I"));
        scalars_["5"] += scalars_["1"] * t0_1;
        scalars_["5"] += 1.00 * dot(Id["bb_oo"]("j,k") * eri["abab_oooo"]("i,j,m,k"), Id["aa_oo"]("i,m"));
        scalars_["5"] += 0.25 * dot(eri["aaaa_oovv"]("o,i,a,d"), t2["aaaa"]("a,d,i,o"));
        scalars_["5"] += 0.50 * dot(Id["bb_oo"]("j,k") * eri["bbbb_oooo"]("j,l,k,n"), Id["bb_oo"]("l,n"));
        scalars_["5"] += 0.50 * dot(Id["aa_oo"]("o,I") * eri["aaaa_oooo"]("o,i,I,m"), Id["aa_oo"]("i,m"));
        scalars_["5"] += 0.25 * dot(eri["bbbb_oovv"]("j,l,c,b"), t2["bbbb"]("c,b,l,j"));
        scalars_["5"] -= 1.00 * dot(Id["bb_oo"]("j,k"), f["bb_oo"]("j,k"));
        scalars_["5"] += scalars_["2"] * t0_1;
        scalars_["5"] -= 1.00 * dot(eri["abab_oovv"]("i,j,a,b"), t2["abab"]("a,b,i,j"));
        scalars_["6"]  = -0.25 * dot(eri["aaaa_oovv"]("i,l,a,d"), t2_1["aaaa"]("a,d,l,i"));
        scalars_["6"] += t0_1 * w0;
        scalars_["6"] += 1.00 * dot(eri["abab_oovv"]("l,j,a,c"), t2_1["abab"]("a,c,l,j"));
        scalars_["6"] -= 1.00 * dot(dp["aa_ov"]("i,a"), t1_1["aa"]("a,i")) * t0_1;
        scalars_["6"] += 1.00 * dot(f["bb_ov"]("j,b"), t1_1["bb"]("b,j"));
        scalars_["6"] -= 1.00 * dot(dp["bb_ov"]("j,b"), t1_1["bb"]("b,j")) * t0_1;
        scalars_["6"] -= 0.25 * dot(eri["bbbb_oovv"]("j,k,b,c"), t2_1["bbbb"]("b,c,k,j"));
        scalars_["6"] += 1.00 * dot(f["aa_ov"]("i,a"), t1_1["aa"]("a,i"));

    }

    {
        reused_["1_bbaa_voov"]("a,i,j,b")  = t2["bbbb"]("c,a,i,k") * eri["abab_oovv"]("j,k,b,c");
        reused_["2_bbbb_vovo"]("a,i,b,j")  = t2["abab"]("c,a,k,i") * reused_["1_bbaa_voov"]("b,j,k,c");
        reused_["3_bb_vv"]("a,b")  = eri["abab_oovv"]("i,j,c,a") * t2["abab"]("c,b,i,j");
        reused_["4_bbbb_voov"]("a,i,j,b")  = t2["bbbb"]("c,a,i,j") * reused_["3_bb_vv"]("c,b");
        reused_["5_aaaa_vvoo"]("a,b,i,j")  = eri["abab_oovv"]("k,l,c,d") * t2["abab"]("a,d,k,l") * t2["aaaa"]("c,b,i,j");
        reused_["6_aabb_voov"]("a,i,j,b")  = t2["aaaa"]("c,a,i,k") * eri["abab_oovv"]("k,j,c,b");
        reused_["7_aaaa_vovo"]("a,i,b,j")  = t2["abab"]("a,c,i,k") * reused_["6_aabb_voov"]("b,j,k,c");
        reused_["8_bb_ov"]("i,a")  = eri["abab_oovv"]("j,i,b,a") * t1_1["aa"]("b,j");
        reused_["9_aa_ov"]("i,a")  = eri["abab_oovv"]("i,j,a,b") * t1_1["bb"]("b,j");
        reused_["10_aa_ov"]("i,a")  = eri["aaaa_oovv"]("j,i,b,a") * t1_1["aa"]("b,j");
        reused_["11_bb_ov"]("i,a")  = eri["bbbb_oovv"]("j,i,b,a") * t1_1["bb"]("b,j");
        reused_["12_bb_ov"]("i,a")  = t0_1 * dp["bb_ov"]("i,a");
        reused_["13_aa_ov"]("i,a")  = t0_1 * dp["aa_ov"]("i,a");
        reused_["14_aaaa_vovv"]("a,i,b,c")  = eri["aaaa_oovo"]("j,k,a,i") * t2["aaaa"]("b,c,k,j");
        reused_["15_aa_oo"]("i,j")  = t2["aaaa"]("a,b,i,k") * eri["aaaa_oovv"]("k,j,a,b");
        reused_["16_aa_oo"]("i,j")  = eri["abab_oovv"]("i,k,a,b") * t2["abab"]("a,b,j,k");
        reused_["17_aa_oo"]("i,j")  = t0_1 * dp["aa_oo"]("i,j");
        reused_["18_aa_vv"]("a,b")  = eri["aaaa_oovv"]("i,j,a,c") * t2["aaaa"]("c,b,j,i");
        reused_["19_aaaa_voov"]("a,i,j,b")  = reused_["18_aa_vv"]("c,b") * t2["aaaa"]("c,a,i,j");
        reused_["20_aaaa_vvoo"]("a,b,i,j")  = eri["aaaa_vvvv"]("a,b,c,d") * t2["aaaa"]("c,d,i,j");
        reused_["21_aaaa_oooo"]("i,j,k,l")  = eri["aaaa_oovv"]("i,j,a,b") * t2["aaaa"]("a,b,k,l");
        reused_["22_aaaa_oooo"]("i,j,k,l")  = t2_1["aaaa"]("a,b,i,j") * eri["aaaa_oovv"]("k,l,a,b");
        reused_["23_aa_vv"]("a,b")  = t2["aaaa"]("c,a,i,j") * eri["aaaa_oovv"]("j,i,c,b");
        reused_["24_aaaa_vvoo"]("a,b,i,j")  = reused_["23_aa_vv"]("a,c") * t2["aaaa"]("c,b,i,j");
        reused_["25_aaaa_oovv"]("i,j,a,b")  = (reused_["21_aaaa_oooo"]("k,l,i,j") + 2.00 * eri["aaaa_oooo"]("k,l,i,j")) * t2["aaaa"]("a,b,l,k");
        reused_["26_bbbb_voov"]("a,i,j,b")  = t2["abab"]("c,a,k,i") * eri["abab_oovv"]("k,j,c,b");
        reused_["27_bbaa_voov"]("a,i,j,b")  = t2["abab"]("c,a,k,i") * eri["aaaa_oovv"]("k,j,c,b");
        reused_["28_bbbb_voov"]("a,i,j,b")  = t2["bbbb"]("c,a,i,k") * eri["bbbb_oovv"]("k,j,c,b");
        reused_["29_bb_vv"]("a,b")  = t2["bbbb"]("c,a,i,j") * eri["bbbb_oovv"]("j,i,c,b");
        reused_["30_bb_oo"]("i,j")  = eri["abab_oovv"]("k,i,a,b") * t2["abab"]("a,b,k,j");
        reused_["31_bb_oo"]("i,j")  = t2["bbbb"]("a,b,i,k") * eri["bbbb_oovv"]("k,j,a,b");
        reused_["32_bb_vo"]("a,i")  = dp["aa_ov"]("j,b") * t2["abab"]("b,a,j,i");
        reused_["33_bb_vo"]("a,i")  = dp["bb_ov"]("j,b") * t2["bbbb"]("b,a,i,j");
        reused_["34_bb_vv"]("a,b")  = t0_1 * dp["bb_vv"]("a,b");
        reused_["35_bb_oo"]("i,j")  = t0_1 * dp["bb_oo"]("i,j");
        reused_["36_bbbb_oooo"]("i,j,k,l")  = eri["bbbb_oovv"]("i,j,a,b") * t2["bbbb"]("a,b,k,l");
        reused_["37_bbbb_vvvo"]("a,b,c,i")  = (eri["bbbb_vvvv"]("a,b,d,c") + -0.50 * t2["bbbb"]("a,b,j,k") * eri["bbbb_oovv"]("k,j,d,c")) * t1_1["bb"]("d,i");
        reused_["38_bb_oo"]("i,j")  = dp["bb_ov"]("i,a") * t1_1["bb"]("a,j");
        reused_["39_bbbb_ovvo"]("i,a,b,j")  = dp["bb_oo"]("k,i") * t2["bbbb"]("a,b,j,k");
        reused_["40_bb_oo"]("i,j")  = f["bb_ov"]("i,a") * t1_1["bb"]("a,j");
        reused_["41_bbbb_vovv"]("a,i,b,c")  = eri["bbbb_oovo"]("j,k,a,i") * t2["bbbb"]("b,c,k,j");
        reused_["42_bbbb_vovv"]("a,i,b,c")  = eri["bbbb_oovo"]("j,k,a,i") * t2_1["bbbb"]("b,c,k,j");
        reused_["43_abab_oovo"]("i,j,a,k")  = eri["abab_oovv"]("i,j,a,b") * t1_1["bb"]("b,k");
        reused_["44_bbbb_ooov"]("i,j,k,a")  = t1_1["bb"]("b,i") * eri["bbbb_oovv"]("j,k,b,a");
        reused_["45_bbbb_oooo"]("i,j,k,l")  = eri["bbbb_oovo"]("i,j,a,k") * t1_1["bb"]("a,l");
        reused_["46_bb_ov"]("i,a")  = eri["bbbb_oovv"]("i,j,b,a") * t1_1["bb"]("b,j");
        reused_["47_bb_oo"]("i,j")  = t2_1["bbbb"]("a,b,i,k") * eri["bbbb_oovv"]("k,j,a,b");
        reused_["47_bb_oo"]("i,j") -= 2.00 * eri["bbbb_oovo"]("k,j,a,i") * t1_1["bb"]("a,k");
        reused_["48_bb_oo"]("i,j")  = eri["abab_oovo"]("k,i,a,j") * t1_1["aa"]("a,k");
        reused_["48_bb_oo"]("i,j") += eri["abab_oovv"]("k,i,a,b") * t2_1["abab"]("a,b,k,j");
        reused_["49_bbbb_ovvo"]("i,a,b,j")  = -2.00 * reused_["28_bbbb_voov"]("a,i,k,c") * t2["bbbb"]("c,b,j,k");
        reused_["49_bbbb_ovvo"]("i,a,b,j") -= 2.00 * f["bb_oo"]("n,i") * t2["bbbb"]("a,b,j,n");
        reused_["49_bbbb_ovvo"]("i,a,b,j") -= 4.00 * dp["bb_ov"]("n,f") * t1_1["bb"]("f,j") * t2["bbbb"]("a,b,i,n");
        reused_["49_bbbb_ovvo"]("i,a,b,j") += 4.00 * dp["bb_oo"]("n,i") * t2_1["bbbb"]("a,b,j,n");
        reused_["49_bbbb_ovvo"]("i,a,b,j") += reused_["31_bb_oo"]("i,k") * t2["bbbb"]("a,b,j,k");
        reused_["49_bbbb_ovvo"]("i,a,b,j") -= 2.00 * reused_["27_bbaa_voov"]("a,i,m,e") * t2["abab"]("e,b,m,j");
        reused_["49_bbbb_ovvo"]("i,a,b,j") -= 2.00 * t2["bbbb"]("a,b,j,k") * reused_["30_bb_oo"]("k,i");
        reused_["50_bbbb_oooo"]("i,j,k,l")  = t2_1["bbbb"]("a,b,i,j") * eri["bbbb_oovv"]("k,l,a,b");
        reused_["51_bbbb_vooo"]("a,i,j,k")  = t2["bbbb"]("b,a,i,j") * dp["bb_ov"]("k,b");
        reused_["52_bbbb_vvoo"]("a,b,i,j")  = dp["bb_vv"]("a,c") * t2["bbbb"]("c,b,i,j");
        reused_["53_bbbb_vooo"]("a,i,j,k")  = t2["bbbb"]("b,a,i,j") * f["bb_ov"]("k,b");
        reused_["53_bbbb_vooo"]("a,i,j,k") -= 0.50 * eri["bbbb_vovv"]("a,k,b,c") * t2["bbbb"]("b,c,i,j");
        reused_["54_aa_vo"]("a,i")  = dp["aa_ov"]("j,b") * t2["aaaa"]("b,a,i,j");
        reused_["55_aaaa_voov"]("a,i,j,b")  = t2["aaaa"]("c,a,i,k") * eri["aaaa_oovv"]("k,j,c,b");
        reused_["56_aabb_voov"]("a,i,j,b")  = t2["abab"]("a,c,i,k") * eri["bbbb_oovv"]("k,j,c,b");
        reused_["57_aa_vv"]("a,b")  = eri["abab_oovv"]("i,j,a,c") * t2["abab"]("b,c,i,j");
        reused_["58_aa_vo"]("a,i")  = dp["bb_ov"]("j,b") * t2["abab"]("a,b,i,j");
        reused_["59_aa_vv"]("a,b")  = t0_1 * dp["aa_vv"]("a,b");
        reused_["60_aaaa_oooo"]("i,j,k,l")  = eri["aaaa_oovo"]("i,j,a,k") * t1_1["aa"]("a,l");
        reused_["61_aaaa_ovvo"]("i,a,b,j")  = Id["aa_oo"]("l,m") * eri["aaaa_oooo"]("l,k,m,i") * t2["aaaa"]("a,b,j,k");
        reused_["62_aaaa_ovvo"]("i,a,b,j")  = Id["aa_oo"]("l,m") * eri["aaaa_oooo"]("k,l,i,m") * t2["aaaa"]("a,b,j,k");
        reused_["63_baba_ovvo"]("i,a,b,j")  = Id["bb_oo"]("l,m") * eri["bbbb_oooo"]("k,l,i,m") * t2["abab"]("a,b,j,k");
        reused_["64_baba_ovvo"]("i,a,b,j")  = Id["bb_oo"]("l,m") * eri["bbbb_oooo"]("l,k,m,i") * t2["abab"]("a,b,j,k");
        reused_["65_aabb_ovvo"]("i,a,b,j")  = Id["aa_oo"]("l,m") * eri["aaaa_oooo"]("l,k,m,i") * t2["abab"]("a,b,k,j");
        reused_["66_aabb_ovvo"]("i,a,b,j")  = Id["aa_oo"]("l,m") * eri["aaaa_oooo"]("k,l,i,m") * t2["abab"]("a,b,k,j");
        reused_["67_bbbb_ovvo"]("i,a,b,j")  = Id["bb_oo"]("l,m") * eri["bbbb_oooo"]("l,k,m,i") * t2["bbbb"]("a,b,j,k");
        reused_["68_bbbb_ovvo"]("i,a,b,j")  = Id["bb_oo"]("l,m") * eri["bbbb_oooo"]("k,l,i,m") * t2["bbbb"]("a,b,j,k");
        reused_["69_aaaa_voov"]("a,i,j,b")  = t2["aaaa"]("c,a,i,k") * eri["aaaa_oovv"]("j,k,b,c");
        reused_["70_aabb_voov"]("a,i,j,b")  = t2["abab"]("a,c,i,k") * eri["bbbb_oovv"]("j,k,b,c");
        reused_["71_aabb_vovo"]("a,i,b,j")  = reused_["69_aaaa_voov"]("a,i,l,d") * t2["abab"]("d,b,l,j");
        reused_["71_aabb_vovo"]("a,i,b,j") -= eri["bbbb_oovv"]("k,m,c,e") * t2["bbbb"]("c,b,m,k") * t2["abab"]("a,e,i,j");
        reused_["71_aabb_vovo"]("a,i,b,j") += reused_["70_aabb_voov"]("a,i,k,c") * t2["bbbb"]("c,b,j,k");
        reused_["72_aabb_ovvo"]("i,a,b,j")  = dp["aa_oo"]("k,i") * t2["abab"]("a,b,k,j");
        reused_["73_abab_vvoo"]("a,b,i,j")  = -1.00 * t2["abab"]("a,b,i,k") * dp["bb_oo"]("k,j");
        reused_["73_abab_vvoo"]("a,b,i,j") += dp["aa_vv"]("a,c") * t2["abab"]("c,b,i,j");
        reused_["74_baab_vvoo"]("a,b,i,j")  = dp["bb_vv"]("a,c") * t2["abab"]("b,c,i,j");
        reused_["75_aaaa_ovvo"]("i,a,b,j")  = dp["aa_oo"]("k,i") * t2["aaaa"]("a,b,j,k");
        reused_["76_aaaa_vvoo"]("a,b,i,j")  = dp["aa_vv"]("a,c") * t2["aaaa"]("c,b,i,j");
        reused_["77_aabb_vooo"]("a,i,j,k")  = t2["abab"]("a,b,i,j") * dp["bb_ov"]("k,b");
        reused_["78_baab_voov"]("a,i,j,b")  = t2["abab"]("c,a,i,k") * eri["abab_oovv"]("j,k,c,b");
        reused_["79_bbbb_vvoo"]("a,b,i,j")  = eri["bbbb_vvvv"]("a,b,c,d") * t2["bbbb"]("c,d,i,j");
        reused_["80_abab_oooo"]("i,j,k,l")  = eri["abab_oovv"]("i,j,a,b") * t2["abab"]("a,b,k,l");
        reused_["81_abba_vovo"]("a,i,b,j")  = t2["abab"]("a,c,k,i") * eri["baba_vovo"]("b,k,c,j");
        reused_["82_aaaa_vooo"]("a,i,j,k")  = t2["aaaa"]("b,a,i,j") * dp["aa_ov"]("k,b");
        reused_["83_baba_vooo"]("a,i,j,k")  = t2["abab"]("b,a,i,j") * dp["aa_ov"]("k,b");
        reused_["84_bb_vo"]("a,i")  = -0.50 * eri["bbbb_vovv"]("a,k,d,c") * t2["bbbb"]("d,c,i,k");
        reused_["84_bb_vo"]("a,i") -= 0.50 * t2["bbbb"]("d,a,m,k") * eri["bbbb_oovo"]("k,m,d,i");
        reused_["84_bb_vo"]("a,i") += f["bb_ov"]("k,d") * t2["bbbb"]("d,a,i,k");
        reused_["84_bb_vo"]("a,i") += t0_1 * dp["bb_vo"]("a,i");
        reused_["84_bb_vo"]("a,i") += eri["abab_oovo"]("l,k,b,i") * t2["abab"]("b,a,l,k");
        reused_["84_bb_vo"]("a,i") -= f["aa_ov"]("j,b") * t2["abab"]("b,a,j,i");
        reused_["84_bb_vo"]("a,i") += eri["baab_vovv"]("a,j,b,c") * t2["abab"]("b,c,j,i");
        reused_["85_aaaa_vovo"]("a,i,b,j")  = -1.00 * eri["aaaa_vovo"]("a,l,d,i") * t2["aaaa"]("d,b,j,l");
        reused_["85_aaaa_vovo"]("a,i,b,j") += eri["abba_vovo"]("a,k,c,i") * t2["abab"]("b,c,j,k");
        reused_["86_bbbb_vovo"]("a,i,b,j")  = -1.00 * eri["baab_vovo"]("a,l,d,i") * t2["abab"]("d,b,l,j");
        reused_["86_bbbb_vovo"]("a,i,b,j") += eri["bbbb_vovo"]("a,k,c,i") * t2["bbbb"]("c,b,j,k");
        reused_["87_aa_vo"]("a,i")  = -2.00 * f["aa_ov"]("k,c") * t2["aaaa"]("c,a,i,k");
        reused_["87_aa_vo"]("a,i") += 2.00 * t2["abab"]("a,b,l,j") * eri["abba_oovo"]("l,j,b,i");
        reused_["87_aa_vo"]("a,i") += eri["aaaa_vovv"]("a,k,c,e") * t2["aaaa"]("c,e,i,k");
        reused_["87_aa_vo"]("a,i") += t2["aaaa"]("c,a,l,k") * eri["aaaa_oovo"]("k,l,c,i");
        reused_["87_aa_vo"]("a,i") += 2.00 * eri["abab_vovv"]("a,j,c,d") * t2["abab"]("c,d,i,j");
        reused_["87_aa_vo"]("a,i") -= 2.00 * t0_1 * dp["aa_vo"]("a,i");
        reused_["87_aa_vo"]("a,i") += 2.00 * f["bb_ov"]("j,b") * t2["abab"]("a,b,i,j");
        reused_["88_bb_vo"]("a,i")  = -1.00 * dp["bb_ov"]("k,c") * t2_1["bbbb"]("c,a,i,k");
        reused_["88_bb_vo"]("a,i") -= t1_1["bb"]("a,k") * dp["bb_oo"]("k,i");
        reused_["88_bb_vo"]("a,i") += dp["bb_vv"]("a,c") * t1_1["bb"]("c,i");
        reused_["88_bb_vo"]("a,i") += dp["aa_ov"]("j,b") * t2_1["abab"]("b,a,j,i");
        reused_["89_aa_ov"]("i,a")  = -1.00 * dp["aa_vv"]("a,c") * t1_1["aa"]("c,i");
        reused_["89_aa_ov"]("i,a") += dp["aa_oo"]("k,i") * t1_1["aa"]("a,k");
        reused_["89_aa_ov"]("i,a") += dp["aa_ov"]("k,c") * t2_1["aaaa"]("c,a,i,k");
        reused_["89_aa_ov"]("i,a") -= dp["bb_ov"]("j,b") * t2_1["abab"]("a,b,i,j");
        reused_["90_aa_vo"]("a,i")  = (scalars_["2"] + scalars_["1"]) * t1_1["aa"]("a,i");
        reused_["91_bb_vo"]("a,i")  = (scalars_["1"] + scalars_["2"]) * t1_1["bb"]("a,i");
        reused_["92_bbbb_voov"]("a,i,j,b")  = t2["bbbb"]("c,a,i,j") * reused_["29_bb_vv"]("b,c");
        reused_["93_aabb_vovo"]("a,i,b,j")  = reused_["55_aaaa_voov"]("a,i,l,d") * t2["abab"]("d,b,l,j");
        reused_["93_aabb_vovo"]("a,i,b,j") += reused_["56_aabb_voov"]("a,i,k,c") * t2["bbbb"]("c,b,j,k");
        reused_["94_bbbb_oovv"]("i,j,a,b")  = (reused_["36_bbbb_oooo"]("k,l,i,j") + 2.00 * eri["bbbb_oooo"]("k,l,i,j")) * t2["bbbb"]("a,b,l,k");
        reused_["95_aabb_voov"]("a,i,j,b")  = t2["abab"]("a,c,i,j") * reused_["3_bb_vv"]("c,b");
        reused_["95_aabb_voov"]("a,i,j,b") -= eri["baab_vovo"]("b,k,e,j") * t2["aaaa"]("e,a,i,k");
        reused_["95_aabb_voov"]("a,i,j,b") += eri["bbbb_vovo"]("b,m,f,j") * t2["abab"]("a,f,i,m");
        reused_["95_aabb_voov"]("a,i,j,b") -= t2["abab"]("a,f,i,j") * f["bb_vv"]("b,f");
        reused_["95_aabb_voov"]("a,i,j,b") += 2.00 * dp["bb_vv"]("b,f") * t2_1["abab"]("a,f,i,j");
        reused_["95_aabb_voov"]("a,i,j,b") -= reused_["1_bbaa_voov"]("b,j,l,d") * t2["aaaa"]("d,a,i,l");
        reused_["95_aabb_voov"]("a,i,j,b") -= 2.00 * reused_["83_baba_vooo"]("b,i,j,k") * t1_1["aa"]("a,k");
        reused_["96_aaaa_ovvo"]("i,a,b,j")  = -2.00 * eri["bbbb_oovv"]("n,m,f,e") * t2["abab"]("a,f,i,n") * t2["abab"]("b,e,j,m");
        reused_["96_aaaa_ovvo"]("i,a,b,j") -= 2.00 * eri["abab_oovv"]("k,n,d,e") * t2["abab"]("d,e,i,n") * t2["aaaa"]("a,b,j,k");
        reused_["96_aaaa_ovvo"]("i,a,b,j") -= 4.00 * dp["aa_ov"]("l,d") * t1_1["aa"]("d,j") * t2["aaaa"]("a,b,i,l");
        reused_["96_aaaa_ovvo"]("i,a,b,j") += 4.00 * dp["aa_oo"]("l,i") * t2_1["aaaa"]("a,b,j,l");
        reused_["96_aaaa_ovvo"]("i,a,b,j") -= 2.00 * f["aa_oo"]("l,i") * t2["aaaa"]("a,b,j,l");
        reused_["96_aaaa_ovvo"]("i,a,b,j") += reused_["15_aa_oo"]("i,k") * t2["aaaa"]("a,b,j,k");
        reused_["96_aaaa_ovvo"]("i,a,b,j") -= 2.00 * reused_["55_aaaa_voov"]("a,i,k,c") * t2["aaaa"]("c,b,j,k");
        reused_["97_aaaa_vvoo"]("a,b,i,j")  = t1_1["aa"]("a,k") * reused_["82_aaaa_vooo"]("b,i,j,k");
        reused_["97_aaaa_vvoo"]("a,b,i,j") -= 0.50 * f["aa_vv"]("b,c") * t2["aaaa"]("c,a,i,j");
        reused_["97_aaaa_vvoo"]("a,b,i,j") += dp["aa_vv"]("b,c") * t2_1["aaaa"]("c,a,i,j");
        reused_["98_abba_vvoo"]("a,b,i,j")  = t2["abab"]("a,b,k,i") * reused_["16_aa_oo"]("k,j");
        reused_["98_abba_vvoo"]("a,b,i,j") -= eri["abba_vovo"]("a,l,c,j") * t2["bbbb"]("c,b,i,l");
        reused_["98_abba_vvoo"]("a,b,i,j") -= 2.00 * dp["bb_ov"]("l,c") * t1_1["bb"]("c,i") * t2["abab"]("a,b,j,l");
        reused_["98_abba_vvoo"]("a,b,i,j") -= 2.00 * t2_1["abab"]("a,b,m,i") * dp["aa_oo"]("m,j");
        reused_["98_abba_vvoo"]("a,b,i,j") += eri["aaaa_vovo"]("a,m,d,j") * t2["abab"]("d,b,m,i");
        reused_["98_abba_vvoo"]("a,b,i,j") += f["aa_oo"]("m,j") * t2["abab"]("a,b,m,i");
        reused_["98_abba_vvoo"]("a,b,i,j") -= reused_["26_bbbb_voov"]("b,i,l,c") * t2["abab"]("a,c,j,l");
        reused_["98_abba_vvoo"]("a,b,i,j") -= 0.50 * reused_["15_aa_oo"]("j,k") * t2["abab"]("a,b,k,i");
        reused_["99_abab_vvoo"]("a,b,i,j")  = -0.50 * reused_["23_aa_vv"]("a,c") * t2["abab"]("c,b,i,j");
        reused_["99_abab_vvoo"]("a,b,i,j") -= reused_["80_abab_oooo"]("k,l,i,j") * t2["abab"]("a,b,k,l");
        reused_["99_abab_vvoo"]("a,b,i,j") -= 2.00 * reused_["77_aabb_vooo"]("a,i,j,l") * t1_1["bb"]("b,l");
        reused_["99_abab_vvoo"]("a,b,i,j") -= t2["abab"]("a,f,m,j") * reused_["78_baab_voov"]("b,i,m,f");
        reused_["99_abab_vvoo"]("a,b,i,j") -= eri["abab_oooo"]("k,l,i,j") * t2["abab"]("a,b,k,l");
        reused_["99_abab_vvoo"]("a,b,i,j") -= eri["abab_vvvv"]("a,b,d,e") * t2["abab"]("d,e,i,j");
        reused_["99_abab_vvoo"]("a,b,i,j") -= f["aa_vv"]("a,d") * t2["abab"]("d,b,i,j");
        reused_["99_abab_vvoo"]("a,b,i,j") -= 2.00 * dp["bb_oo"]("l,j") * t2_1["abab"]("a,b,i,l");
        reused_["99_abab_vvoo"]("a,b,i,j") += f["bb_oo"]("l,j") * t2["abab"]("a,b,i,l");
        reused_["99_abab_vvoo"]("a,b,i,j") -= 2.00 * dp["aa_ov"]("m,d") * t1_1["aa"]("d,i") * t2["abab"]("a,b,m,j");
        reused_["99_abab_vvoo"]("a,b,i,j") -= 0.50 * eri["bbbb_oovv"]("l,n,f,e") * t2["bbbb"]("f,e,j,l") * t2["abab"]("a,b,i,n");
        reused_["99_abab_vvoo"]("a,b,i,j") += 2.00 * dp["aa_vv"]("a,d") * t2_1["abab"]("d,b,i,j");
        reused_["99_abab_vvoo"]("a,b,i,j") += eri["abab_vovo"]("a,l,d,j") * t2["abab"]("d,b,i,l");
        reused_["99_abab_vvoo"]("a,b,i,j") += reused_["30_bb_oo"]("n,j") * t2["abab"]("a,b,i,n");
        reused_["100_bbbb_voov"]("a,i,j,b")  = -0.50 * f["bb_vv"]("a,c") * t2["bbbb"]("c,b,i,j");
        reused_["100_bbbb_voov"]("a,i,j,b") += dp["bb_vv"]("a,c") * t2_1["bbbb"]("c,b,i,j");
        reused_["100_bbbb_voov"]("a,i,j,b") += reused_["51_bbbb_vooo"]("a,i,j,k") * t1_1["bb"]("b,k");
        reused_["101_bbaa_voov"]("a,i,j,b")  = t2_1["abab"]("c,a,k,i") * eri["aaaa_oovv"]("k,j,c,b");
        reused_["102_bbbb_voov"]("a,i,j,b")  = t2_1["abab"]("c,a,k,i") * eri["abab_oovv"]("k,j,c,b");
        reused_["103_bbaa_voov"]("a,i,j,b")  = t2_1["bbbb"]("c,a,i,k") * eri["abab_oovv"]("j,k,b,c");
        reused_["104_bbbb_voov"]("a,i,j,b")  = t2_1["bbbb"]("c,a,i,k") * eri["bbbb_oovv"]("k,j,c,b");
        reused_["105_baab_vovo"]("a,i,b,j")  = eri["baab_vovv"]("a,i,b,c") * t1_1["bb"]("c,j");
        reused_["106_bbbb_vovo"]("a,i,b,j")  = eri["bbbb_vovv"]("a,i,c,b") * t1_1["bb"]("c,j");
        reused_["107_bb_vv"]("a,b")  = -2.00 * eri["bbbb_vovv"]("a,j,c,b") * t1_1["bb"]("c,j");
        reused_["107_bb_vv"]("a,b") += eri["bbbb_oovv"]("j,i,c,b") * t2_1["bbbb"]("c,a,i,j");
        reused_["108_bb_vv"]("a,b")  = eri["baab_vovv"]("a,i,c,b") * t1_1["aa"]("c,i");
        reused_["108_bb_vv"]("a,b") += t2_1["abab"]("c,a,j,k") * eri["abab_oovv"]("j,k,c,b");
        reused_["109_bb_vv"]("a,b")  = eri["bbbb_oovv"]("i,j,a,c") * t2["bbbb"]("c,b,j,i");
        reused_["110_bbbb_voov"]("a,i,j,b")  = reused_["109_bb_vv"]("c,b") * t2["bbbb"]("c,a,i,j");
        reused_["111_bbbb_oovo"]("i,j,a,k")  = eri["bbbb_oovv"]("i,j,a,b") * t1_1["bb"]("b,k");
        reused_["112_bb_vv"]("a,b")  = eri["bbbb_vovv"]("a,i,b,c") * t1_1["bb"]("c,i");
        reused_["112_bb_vv"]("a,b") -= 0.50 * eri["bbbb_oovv"]("i,j,b,c") * t2_1["bbbb"]("c,a,j,i");
        reused_["113_aaaa_voov"]("a,i,j,b")  = t2_1["aaaa"]("c,a,i,k") * eri["aaaa_oovv"]("k,j,c,b");
        reused_["114_aabb_voov"]("a,i,j,b")  = t2_1["aaaa"]("c,a,i,k") * eri["abab_oovv"]("k,j,c,b");
        reused_["115_aabb_voov"]("a,i,j,b")  = t2_1["abab"]("a,c,i,k") * eri["bbbb_oovv"]("k,j,c,b");
        reused_["116_abba_oovo"]("i,j,a,k")  = eri["abab_oovv"]("i,j,b,a") * t1_1["aa"]("b,k");
        reused_["117_aaaa_vovo"]("a,i,b,j")  = eri["aaaa_vovv"]("a,i,c,b") * t1_1["aa"]("c,j");
        reused_["118_abba_vovo"]("a,i,b,j")  = eri["abab_vovv"]("a,i,c,b") * t1_1["aa"]("c,j");
        reused_["119_aaaa_ooov"]("i,j,k,a")  = t1_1["aa"]("b,i") * eri["aaaa_oovv"]("j,k,b,a");
        reused_["120_aa_oo"]("i,j")  = dp["aa_ov"]("i,a") * t1_1["aa"]("a,j");
        reused_["121_aa_oo"]("i,j")  = t1_1["aa"]("a,i") * f["aa_ov"]("j,a");
        reused_["122_aa_vv"]("a,b")  = -0.50 * t2_1["aaaa"]("c,a,j,i") * eri["aaaa_oovv"]("i,j,c,b");
        reused_["122_aa_vv"]("a,b") += eri["aaaa_vovv"]("a,i,c,b") * t1_1["aa"]("c,i");
        reused_["123_aa_vv"]("a,b")  = -1.00 * t2_1["abab"]("a,c,j,i") * eri["abab_oovv"]("j,i,b,c");
        reused_["123_aa_vv"]("a,b") += eri["abab_vovv"]("a,i,b,c") * t1_1["bb"]("c,i");
        reused_["124_aa_oo"]("i,j")  = -0.50 * eri["aaaa_oovv"]("k,i,a,b") * t2_1["aaaa"]("a,b,j,k");
        reused_["124_aa_oo"]("i,j") += eri["aaaa_oovo"]("k,i,a,j") * t1_1["aa"]("a,k");
        reused_["125_aa_oo"]("i,j")  = -1.00 * eri["abab_oovv"]("i,k,b,c") * t2_1["abab"]("b,c,j,k");
        reused_["125_aa_oo"]("i,j") += eri["abba_oovo"]("i,k,a,j") * t1_1["bb"]("a,k");
        reused_["126_aaaa_vvvo"]("a,b,c,i")  = (eri["aaaa_vvvv"]("a,b,d,c") + -0.50 * t2["aaaa"]("a,b,j,k") * eri["aaaa_oovv"]("k,j,d,c")) * t1_1["aa"]("d,i");
        reused_["127_aaaa_vovv"]("a,i,b,c")  = eri["aaaa_oovo"]("j,k,a,i") * t2_1["aaaa"]("b,c,k,j");
        reused_["128_aa_ov"]("i,a")  = eri["aaaa_oovv"]("i,j,b,a") * t1_1["aa"]("b,j");
        reused_["129_abab_oovv"]("i,j,a,b")  = 2.00 * eri["abab_oovv"]("i,j,a,b");
        reused_["130_aaaa_vooo"]("a,i,j,k")  = t2["aaaa"]("b,a,i,j") * f["aa_ov"]("k,b");
        reused_["130_aaaa_vooo"]("a,i,j,k") -= 0.50 * eri["aaaa_vovv"]("a,k,b,c") * t2["aaaa"]("b,c,i,j");
        reused_["131_aaaa_oovo"]("i,j,a,k")  = eri["aaaa_oovv"]("i,j,a,b") * t1_1["aa"]("b,k");
        reused_["132_aa_vv"]("a,b")  = eri["aaaa_vovv"]("a,i,b,c") * t1_1["aa"]("c,i");
        reused_["132_aa_vv"]("a,b") -= 0.50 * eri["aaaa_oovv"]("i,j,b,c") * t2_1["aaaa"]("c,a,j,i");
        reused_["133_aaaa_ooov"]("i,j,k,a")  = -1.00 * eri["aaaa_oovv"]("l,i,b,c") * t1_1["aa"]("c,l") * t2["aaaa"]("b,a,j,k");
        reused_["133_aaaa_ooov"]("i,j,k,a") += eri["aaaa_oooo"]("l,i,j,k") * t1_1["aa"]("a,l");
        reused_["133_aaaa_ooov"]("i,j,k,a") += 0.50 * reused_["21_aaaa_oooo"]("l,i,j,k") * t1_1["aa"]("a,l");
        reused_["134_aaaa_vooo"]("a,i,j,k")  = t2_1["aaaa"]("b,a,i,j") * dp["aa_ov"]("k,b");
        reused_["135_aaaa_vooo"]("a,i,j,k")  = eri["aaaa_vovv"]("a,i,b,c") * t2_1["aaaa"]("b,c,j,k");
        reused_["135_aaaa_vooo"]("a,i,j,k") -= 2.00 * t2_1["aaaa"]("b,a,j,k") * f["aa_ov"]("i,b");
        reused_["136_aaaa_vooo"]("a,i,j,k")  = t2["aaaa"]("b,a,i,j") * reused_["9_aa_ov"]("k,b");
        reused_["137_bbbb_ovoo"]("i,a,j,k")  = eri["bbbb_oovv"]("l,i,b,c") * t1_1["bb"]("c,l") * t2["bbbb"]("b,a,j,k");
        reused_["137_bbbb_ovoo"]("i,a,j,k") -= eri["bbbb_oooo"]("l,i,j,k") * t1_1["bb"]("a,l");
        reused_["137_bbbb_ovoo"]("i,a,j,k") -= 0.50 * reused_["36_bbbb_oooo"]("l,i,j,k") * t1_1["bb"]("a,l");
        reused_["138_bbbb_vooo"]("a,i,j,k")  = t2_1["bbbb"]("b,a,i,j") * dp["bb_ov"]("k,b");
        reused_["139_bbbb_vooo"]("a,i,j,k")  = t2_1["bbbb"]("b,a,i,j") * f["bb_ov"]("k,b");
        reused_["139_bbbb_vooo"]("a,i,j,k") -= 0.50 * eri["bbbb_vovv"]("a,k,b,c") * t2_1["bbbb"]("b,c,i,j");
        reused_["140_bbbb_vooo"]("a,i,j,k")  = t2["bbbb"]("b,a,i,j") * reused_["8_bb_ov"]("k,b");
        reused_["141_aaaa_vvvo"]("a,b,c,i")  = eri["aaaa_vovv"]("a,j,d,b") * t2_1["aaaa"]("d,c,i,j");
        reused_["142_aaaa_oovo"]("i,j,a,k")  = -1.00 * reused_["55_aaaa_voov"]("a,j,i,c") * t1_1["aa"]("c,k");
        reused_["142_aaaa_oovo"]("i,j,a,k") += eri["aaaa_oovo"]("l,i,b,j") * t2_1["aaaa"]("b,a,k,l");
        reused_["143_aaaa_voov"]("a,i,j,b")  = t2["abab"]("a,c,i,k") * eri["abab_oovv"]("j,k,b,c");
        reused_["144_aaaa_vvvo"]("a,b,c,i")  = eri["aaaa_vovv"]("a,j,d,b") * t2["aaaa"]("d,c,i,j");
        reused_["145_aaaa_vvvo"]("a,b,c,i")  = eri["abab_vovv"]("a,j,b,d") * t2["abab"]("c,d,i,j");
        reused_["146_aaaa_vvvo"]("a,b,c,i")  = eri["abab_vovv"]("a,j,b,d") * t2_1["abab"]("c,d,i,j");
        reused_["147_aa_ov"]("i,a")  = eri["aaaa_oovv"]("i,j,a,b") * t1_1["aa"]("b,j");
        reused_["148_bb_ov"]("i,a")  = eri["bbbb_oovv"]("i,j,a,b") * t1_1["bb"]("b,j");
        reused_["149_aaaa_oovo"]("i,j,a,k")  = eri["aaaa_oovo"]("l,i,b,j") * t2["aaaa"]("b,a,k,l");
        reused_["150_aaaa_oovo"]("i,j,a,k")  = eri["abba_oovo"]("i,l,b,j") * t2["abab"]("a,b,k,l");
        reused_["151_aaaa_vooo"]("a,i,j,k")  = eri["aaaa_vovo"]("a,i,b,j") * t1_1["aa"]("b,k");
        reused_["152_aa_vo"]("a,i")  = t0_1 * reused_["58_aa_vo"]("a,i");
        reused_["153_aa_vo"]("a,i")  = -0.50 * t2["aaaa"]("b,a,j,k") * reused_["131_aaaa_oovo"]("k,j,b,i");
        reused_["153_aa_vo"]("a,i") -= 2.00 * t1_1["aa"]("a,k") * reused_["120_aa_oo"]("k,i");
        reused_["153_aa_vo"]("a,i") += reused_["6_aabb_voov"]("a,i,m,d") * t1_1["bb"]("d,m");
        reused_["153_aa_vo"]("a,i") -= t0_1 * reused_["89_aa_ov"]("i,a");
        reused_["153_aa_vo"]("a,i") += reused_["116_abba_oovo"]("j,l,e,i") * t2["abab"]("a,e,j,l");
        reused_["153_aa_vo"]("a,i") += t1_1["aa"]("a,j") * reused_["16_aa_oo"]("j,i");
        reused_["153_aa_vo"]("a,i") -= w0 * t1_1["aa"]("a,i");
        reused_["153_aa_vo"]("a,i") += f["aa_ov"]("k,b") * t2_1["aaaa"]("b,a,i,k");
        reused_["153_aa_vo"]("a,i") += scalars_["3"] * t1_1["aa"]("a,i");
        reused_["153_aa_vo"]("a,i") -= t2_1["abab"]("a,e,j,l") * eri["abba_oovo"]("j,l,e,i");
        reused_["153_aa_vo"]("a,i") += scalars_["4"] * t1_1["aa"]("a,i");
        reused_["153_aa_vo"]("a,i") += eri["aaaa_vovo"]("a,k,b,i") * t1_1["aa"]("b,k");
        reused_["153_aa_vo"]("a,i") -= f["bb_ov"]("l,e") * t2_1["abab"]("a,e,i,l");
        reused_["153_aa_vo"]("a,i") += f["aa_oo"]("k,i") * t1_1["aa"]("a,k");
        reused_["153_aa_vo"]("a,i") -= 0.50 * eri["aaaa_vovv"]("a,k,b,c") * t2_1["aaaa"]("b,c,i,k");
        reused_["153_aa_vo"]("a,i") += eri["abba_vovo"]("a,l,e,i") * t1_1["bb"]("e,l");
        reused_["153_aa_vo"]("a,i") -= 0.50 * eri["aaaa_oovo"]("k,j,b,i") * t2_1["aaaa"]("b,a,j,k");
        reused_["153_aa_vo"]("a,i") -= eri["abab_vovv"]("a,l,b,d") * t2_1["abab"]("b,d,i,l");
        reused_["153_aa_vo"]("a,i") -= f["aa_vv"]("a,b") * t1_1["aa"]("b,i");
        reused_["153_aa_vo"]("a,i") -= 0.50 * reused_["15_aa_oo"]("i,j") * t1_1["aa"]("a,j");
        reused_["154_aaaa_vooo"]("a,i,j,k")  = reused_["143_aaaa_voov"]("a,i,j,b") * t1_1["aa"]("b,k");
        reused_["154_aaaa_vooo"]("a,i,j,k") += eri["abba_oovo"]("j,l,c,i") * t2_1["abab"]("a,c,k,l");
        reused_["155_aa_vo"]("a,i")  = reused_["147_aa_ov"]("j,b") * t2["aaaa"]("b,a,i,j");
        reused_["156_aa_vo"]("a,i")  = t0_1 * reused_["54_aa_vo"]("a,i");
        reused_["157_aa_vo"]("a,i")  = (reused_["8_bb_ov"]("j,b") + reused_["148_bb_ov"]("j,b")) * t2["abab"]("a,b,i,j");
        reused_["158_bbaa_voov"]("a,i,j,b")  = t2["abab"]("c,a,k,i") * eri["aaaa_oovv"]("j,k,b,c");
        reused_["159_bbaa_voov"]("a,i,j,b")  = t2_1["abab"]("c,a,k,i") * eri["aaaa_oovv"]("j,k,b,c");
        reused_["160_aaaa_vovo"]("a,i,b,j")  = eri["aaaa_vovv"]("a,i,b,c") * t1_1["aa"]("c,j");
        reused_["160_aaaa_vovo"]("a,i,b,j") += t2_1["aaaa"]("c,a,j,k") * eri["aaaa_oovv"]("i,k,b,c");
        reused_["161_bb_oo"]("i,j")  = eri["bbbb_oovv"]("i,k,a,b") * t2["bbbb"]("a,b,j,k");
        reused_["162_aabb_voov"]("a,i,j,b")  = t2_1["abab"]("a,c,i,k") * eri["bbbb_oovv"]("j,k,b,c");
        reused_["163_bbbb_voov"]("a,i,j,b")  = t2["bbbb"]("c,a,i,k") * eri["bbbb_oovv"]("j,k,b,c");
        reused_["164_bbbb_voov"]("a,i,j,b")  = t2_1["bbbb"]("c,a,i,k") * eri["bbbb_oovv"]("j,k,b,c");
        reused_["164_bbbb_voov"]("a,i,j,b") += eri["bbbb_vovv"]("a,j,b,c") * t1_1["bb"]("c,i");
        reused_["165_aaaa_vovv"]("a,i,b,c")  = t2["aaaa"]("d,a,i,j") * eri["aaaa_vovv"]("b,j,c,d");
        reused_["166_aabb_vvvo"]("a,b,c,i")  = eri["aaaa_vovv"]("a,j,b,d") * t2["abab"]("d,c,j,i");
        reused_["167_aabb_vvvo"]("a,b,c,i")  = eri["abab_vovv"]("a,j,b,d") * t2["bbbb"]("d,c,i,j");
        reused_["168_abba_vovv"]("a,i,b,c")  = t2["abab"]("a,d,j,i") * eri["baab_vovv"]("b,j,c,d");
        reused_["169_abab_vovv"]("a,i,b,c")  = eri["abab_oovo"]("j,k,a,i") * t2["abab"]("b,c,j,k");
        reused_["170_aaaa_vooo"]("a,i,j,k")  = t2["aaaa"]("b,a,i,l") * eri["aaaa_oovo"]("j,l,b,k");
        reused_["171_abba_oovo"]("i,j,a,k")  = eri["abab_oovo"]("i,l,b,j") * t2["abab"]("b,a,k,l");
        reused_["172_bbaa_vooo"]("a,i,j,k")  = t2["abab"]("b,a,l,i") * eri["aaaa_oovo"]("j,l,b,k");
        reused_["173_bbaa_vooo"]("a,i,j,k")  = t2["bbbb"]("b,a,i,l") * eri["abba_oovo"]("j,l,b,k");
        reused_["174_aa_oo"]("i,j")  = eri["aaaa_oovv"]("i,k,a,b") * t2["aaaa"]("a,b,j,k");
        reused_["175_baab_vooo"]("a,i,j,k")  = eri["baab_vovv"]("a,i,b,c") * t2["abab"]("b,c,j,k");
        reused_["175_baab_vooo"]("a,i,j,k") -= t2["abab"]("b,a,j,k") * f["aa_ov"]("i,b");
        reused_["176_bbbb_vvvo"]("a,b,c,i")  = eri["baab_vovv"]("a,j,d,b") * t2["abab"]("d,c,j,i");
        reused_["177_bbbb_vvvo"]("a,b,c,i")  = eri["bbbb_vovv"]("a,j,d,b") * t2["bbbb"]("d,c,i,j");
        reused_["178_bbbb_oovo"]("i,j,a,k")  = eri["abab_oovo"]("l,i,b,j") * t2["abab"]("b,a,l,k");
        reused_["179_bbbb_oovo"]("i,j,a,k")  = eri["bbbb_oovo"]("l,i,b,j") * t2["bbbb"]("b,a,k,l");
        reused_["180_bb_vo"]("a,i")  = t0_1 * (reused_["32_bb_vo"]("a,i") + -1.00 * reused_["33_bb_vo"]("a,i"));
        reused_["181_bb_oo"]("i,j")  = -0.50 * eri["bbbb_oovv"]("i,k,a,b") * t2_1["bbbb"]("a,b,j,k");
        reused_["181_bb_oo"]("i,j") += eri["bbbb_oovo"]("i,k,a,j") * t1_1["bb"]("a,k");
        reused_["182_aa_oo"]("i,j")  = -0.50 * eri["aaaa_oovv"]("i,k,a,b") * t2_1["aaaa"]("a,b,j,k");
        reused_["182_aa_oo"]("i,j") += eri["aaaa_oovo"]("i,k,a,j") * t1_1["aa"]("a,k");
        reused_["183_aaaa_voov"]("a,i,j,b")  = t2["aaaa"]("c,a,i,k") * eri["aaaa_oovv"]("k,j,b,c");
        reused_["184_bbaa_voov"]("a,i,j,b")  = t2["abab"]("c,a,k,i") * eri["aaaa_oovv"]("k,j,b,c");
        reused_["185_aabb_vvvo"]("a,b,c,i")  = eri["aaaa_vovv"]("a,j,b,d") * t2_1["abab"]("d,c,j,i");
        reused_["186_bbaa_vooo"]("a,i,j,k")  = t2_1["abab"]("b,a,l,i") * eri["aaaa_oovo"]("j,l,b,k");
        reused_["187_aaaa_vovv"]("a,i,b,c")  = 0.50 * eri["aaaa_vvvv"]("a,b,c,d") * t1_1["aa"]("d,i");
        reused_["187_aaaa_vovv"]("a,i,b,c") += 0.25 * eri["aaaa_oovv"]("j,k,c,d") * t1_1["aa"]("d,i") * t2["aaaa"]("b,a,k,j");
        reused_["187_aaaa_vovv"]("a,i,b,c") += eri["aaaa_vovv"]("b,j,c,d") * t2_1["aaaa"]("d,a,i,j");
        reused_["188_aabb_oovo"]("i,j,a,k")  = -1.00 * t2["abab"]("b,a,j,k") * reused_["147_aa_ov"]("i,b");
        reused_["188_aabb_oovo"]("i,j,a,k") += eri["aaaa_oovv"]("i,l,b,c") * t1_1["aa"]("c,j") * t2["abab"]("b,a,l,k");
        reused_["189_aaaa_vooo"]("a,i,j,k")  = t2["aaaa"]("b,a,i,j") * reused_["147_aa_ov"]("k,b");
        reused_["189_aaaa_vooo"]("a,i,j,k") -= 2.00 * t2_1["aaaa"]("b,a,i,l") * eri["aaaa_oovo"]("k,l,b,j");
        reused_["189_aaaa_vooo"]("a,i,j,k") -= t1_1["aa"]("a,l") * eri["aaaa_oooo"]("k,l,i,j");
        reused_["189_aaaa_vooo"]("a,i,j,k") += 2.00 * eri["aaaa_oovv"]("k,l,b,c") * t1_1["aa"]("c,i") * t2["aaaa"]("b,a,j,l");
        reused_["189_aaaa_vooo"]("a,i,j,k") -= 0.50 * t1_1["aa"]("a,l") * reused_["21_aaaa_oooo"]("k,l,i,j");
        reused_["190_baba_vooo"]("a,i,j,k")  = t2_1["abab"]("b,a,i,j") * dp["aa_ov"]("k,b");
        reused_["191_abba_voov"]("a,i,j,b")  = t2["abab"]("a,c,k,i") * eri["abab_oovv"]("k,j,b,c");
        reused_["192_abba_vovv"]("a,i,b,c")  = t2_1["abab"]("a,d,j,i") * eri["baab_vovv"]("b,j,c,d");
        reused_["193_abab_vovv"]("a,i,b,c")  = eri["abab_oovo"]("j,k,a,i") * t2_1["abab"]("b,c,j,k");
        reused_["194_baab_vooo"]("a,i,j,k")  = eri["baba_vovo"]("a,i,b,j") * t1_1["bb"]("b,k");
        reused_["195_baba_vooo"]("a,i,j,k")  = t2_1["abab"]("b,a,i,j") * f["aa_ov"]("k,b");
        reused_["195_baba_vooo"]("a,i,j,k") -= eri["baab_vovo"]("a,k,b,j") * t1_1["aa"]("b,i");
        reused_["195_baba_vooo"]("a,i,j,k") -= eri["baab_vovv"]("a,k,b,c") * t2_1["abab"]("b,c,i,j");
        reused_["196_aabb_vvvo"]("a,b,c,i")  = eri["abab_vovv"]("a,j,b,d") * t2_1["bbbb"]("d,c,i,j");
        reused_["196_aabb_vvvo"]("a,b,c,i") -= eri["abab_oovv"]("k,j,b,d") * t1_1["bb"]("d,i") * t2["abab"]("a,c,k,j");
        reused_["196_aabb_vvvo"]("a,b,c,i") -= eri["abab_vvvv"]("a,c,b,d") * t1_1["bb"]("d,i");
        reused_["197_baba_vooo"]("a,i,j,k")  = -1.00 * t2_1["abab"]("c,a,i,l") * eri["abab_oovo"]("k,l,c,j");
        reused_["197_baba_vooo"]("a,i,j,k") -= eri["abab_oooo"]("k,l,i,j") * t1_1["bb"]("a,l");
        reused_["197_baba_vooo"]("a,i,j,k") += reused_["9_aa_ov"]("k,c") * t2["abab"]("c,a,i,j");
        reused_["197_baba_vooo"]("a,i,j,k") -= reused_["80_abab_oooo"]("k,l,i,j") * t1_1["bb"]("a,l");
        reused_["197_baba_vooo"]("a,i,j,k") -= reused_["1_bbaa_voov"]("a,j,k,b") * t1_1["aa"]("b,i");
        reused_["198_baab_vooo"]("a,i,j,k")  = reused_["78_baab_voov"]("a,i,j,b") * t1_1["bb"]("b,k");
        reused_["198_baab_vooo"]("a,i,j,k") -= eri["abba_oovo"]("j,l,d,i") * t2_1["bbbb"]("d,a,k,l");
        reused_["199_bb_ov"]("i,a")  = -2.00 * reused_["1_bbaa_voov"]("a,i,m,e") * t1_1["aa"]("e,m");
        reused_["199_bb_ov"]("i,a") += reused_["31_bb_oo"]("i,j") * t1_1["bb"]("a,j");
        reused_["199_bb_ov"]("i,a") -= 2.00 * t2["abab"]("d,a,m,k") * reused_["43_abab_oovo"]("m,k,d,i");
        reused_["199_bb_ov"]("i,a") += 2.00 * f["aa_ov"]("l,d") * t2_1["abab"]("d,a,l,i");
        reused_["199_bb_ov"]("i,a") -= 2.00 * t1_1["bb"]("a,k") * f["bb_oo"]("k,i");
        reused_["199_bb_ov"]("i,a") += eri["bbbb_vovv"]("a,k,b,c") * t2_1["bbbb"]("b,c,i,k");
        reused_["199_bb_ov"]("i,a") -= 2.00 * eri["baab_vovo"]("a,l,d,i") * t1_1["aa"]("d,l");
        reused_["199_bb_ov"]("i,a") -= 2.00 * f["bb_ov"]("k,b") * t2_1["bbbb"]("b,a,i,k");
        reused_["199_bb_ov"]("i,a") += 2.00 * f["bb_vv"]("a,b") * t1_1["bb"]("b,i");
        reused_["199_bb_ov"]("i,a") -= 2.00 * eri["baab_vovv"]("a,l,d,c") * t2_1["abab"]("d,c,l,i");
        reused_["199_bb_ov"]("i,a") -= 2.00 * eri["bbbb_vovo"]("a,k,b,i") * t1_1["bb"]("b,k");
        reused_["199_bb_ov"]("i,a") += eri["bbbb_oovo"]("k,j,b,i") * t2_1["bbbb"]("b,a,j,k");
        reused_["199_bb_ov"]("i,a") -= 2.00 * scalars_["3"] * t1_1["bb"]("a,i");
        reused_["199_bb_ov"]("i,a") -= 2.00 * scalars_["4"] * t1_1["bb"]("a,i");
        reused_["199_bb_ov"]("i,a") += 2.00 * w0 * t1_1["bb"]("a,i");
        reused_["199_bb_ov"]("i,a") -= 2.00 * t2_1["abab"]("d,a,m,k") * eri["abab_oovo"]("m,k,d,i");
        reused_["199_bb_ov"]("i,a") -= 2.00 * t1_1["bb"]("a,j") * reused_["30_bb_oo"]("j,i");
        reused_["199_bb_ov"]("i,a") += 4.00 * t1_1["bb"]("a,k") * reused_["38_bb_oo"]("k,i");
        reused_["199_bb_ov"]("i,a") -= 2.00 * t0_1 * reused_["88_bb_vo"]("a,i");
        reused_["199_bb_ov"]("i,a") += t2["bbbb"]("b,a,j,k") * reused_["111_bbbb_oovo"]("k,j,b,i");
        reused_["199_bb_ov"]("i,a") += 2.00 * reused_["26_bbbb_voov"]("a,i,j,c") * t1_1["bb"]("c,j");
        reused_["200_bb_vo"]("a,i")  = -1.00 * reused_["147_aa_ov"]("l,d") * t2["abab"]("d,a,l,i");
        reused_["200_bb_vo"]("a,i") += reused_["148_bb_ov"]("j,b") * t2["bbbb"]("b,a,i,j");
        reused_["201_aabb_voov"]("a,i,j,b")  = t2["abab"]("a,c,i,k") * eri["bbbb_oovv"]("k,j,b,c");
        reused_["202_bbbb_voov"]("a,i,j,b")  = t2["bbbb"]("c,a,i,k") * eri["bbbb_oovv"]("k,j,b,c");
        reused_["203_aabb_vovv"]("a,i,b,c")  = t2_1["abab"]("a,d,i,j") * eri["bbbb_vovv"]("b,j,c,d");
        reused_["204_bbaa_oovo"]("i,j,a,k")  = eri["bbbb_oovv"]("i,l,b,c") * t1_1["bb"]("c,j") * t2["abab"]("a,b,k,l");
        reused_["205_bbbb_vovv"]("a,i,b,c")  = eri["bbbb_oovv"]("j,k,a,d") * t1_1["bb"]("d,i") * t2["bbbb"]("b,c,k,j");
        reused_["205_bbbb_vovv"]("a,i,b,c") += 4.00 * eri["bbbb_vovv"]("b,j,a,d") * t2_1["bbbb"]("d,c,i,j");
        reused_["205_bbbb_vovv"]("a,i,b,c") += 2.00 * eri["bbbb_vvvv"]("c,b,a,d") * t1_1["bb"]("d,i");
        reused_["206_aabb_vooo"]("a,i,j,k")  = t2["abab"]("a,b,i,j") * reused_["148_bb_ov"]("k,b");
        reused_["206_aabb_vooo"]("a,i,j,k") -= t2_1["abab"]("a,b,i,l") * eri["bbbb_oovo"]("k,l,b,j");
        reused_["207_bbbb_vooo"]("a,i,j,k")  = t2["bbbb"]("b,a,i,j") * reused_["148_bb_ov"]("k,b");
        reused_["207_bbbb_vooo"]("a,i,j,k") += 2.00 * eri["bbbb_oovv"]("k,l,b,c") * t1_1["bb"]("c,i") * t2["bbbb"]("b,a,j,l");
        reused_["207_bbbb_vooo"]("a,i,j,k") -= 2.00 * t2_1["bbbb"]("b,a,i,l") * eri["bbbb_oovo"]("k,l,b,j");
        reused_["207_bbbb_vooo"]("a,i,j,k") -= t1_1["bb"]("a,l") * eri["bbbb_oooo"]("k,l,i,j");
        reused_["207_bbbb_vooo"]("a,i,j,k") -= 0.50 * t1_1["bb"]("a,l") * reused_["36_bbbb_oooo"]("k,l,i,j");
        reused_["208_aabb_vooo"]("a,i,j,k")  = t2_1["abab"]("a,b,i,j") * dp["bb_ov"]("k,b");
        reused_["209_abba_vvvo"]("a,b,c,i")  = eri["abab_vovv"]("a,j,d,b") * t2["abab"]("d,c,i,j");
        reused_["210_bbaa_vvvo"]("a,b,c,i")  = eri["baab_vovv"]("a,j,d,b") * t2["aaaa"]("d,c,i,j");
        reused_["211_bbaa_vvvo"]("a,b,c,i")  = eri["baab_vovv"]("a,j,d,b") * t2_1["aaaa"]("d,c,i,j");
        reused_["212_bbbb_vvvo"]("a,b,c,i")  = eri["baab_vovv"]("a,j,d,b") * t2_1["abab"]("d,c,j,i");
        reused_["213_aabb_vovv"]("a,i,b,c")  = t2["abab"]("a,d,i,j") * eri["bbbb_vovv"]("b,j,c,d");
        reused_["214_bbbb_vovv"]("a,i,b,c")  = t2["bbbb"]("d,a,i,j") * eri["bbbb_vovv"]("b,j,c,d");
        reused_["215_baab_vovv"]("a,i,b,c")  = eri["abba_oovo"]("j,k,a,i") * t2["abab"]("b,c,j,k");
        reused_["216_baab_vovv"]("a,i,b,c")  = eri["abba_oovo"]("j,k,a,i") * t2_1["abab"]("b,c,j,k");
        reused_["217_bbaa_oovo"]("i,j,a,k")  = eri["abab_oovo"]("l,i,b,j") * t2["aaaa"]("b,a,k,l");
        reused_["218_aabb_vooo"]("a,i,j,k")  = t2["abab"]("a,b,i,l") * eri["bbbb_oovo"]("j,l,b,k");
        reused_["219_abba_vooo"]("a,i,j,k")  = t2["abab"]("a,b,l,i") * eri["abba_oovo"]("l,j,b,k");
        reused_["220_bbbb_vooo"]("a,i,j,k")  = t2["bbbb"]("b,a,i,l") * eri["bbbb_oovo"]("j,l,b,k");
        reused_["221_abab_vooo"]("a,i,j,k")  = eri["abba_vovo"]("a,i,b,j") * t1_1["bb"]("b,k");
        reused_["222_bbbb_vooo"]("a,i,j,k")  = eri["bbbb_vovo"]("a,i,b,j") * t1_1["bb"]("b,k");
        reused_["223_abba_vvvo"]("a,b,c,i")  = eri["abab_vovv"]("a,j,d,b") * t2_1["abab"]("d,c,i,j");
        reused_["223_abba_vvvo"]("a,b,c,i") -= eri["abab_vvvv"]("a,c,d,b") * t1_1["aa"]("d,i");
        reused_["223_abba_vvvo"]("a,b,c,i") -= eri["abab_oovv"]("k,j,d,b") * t1_1["aa"]("d,i") * t2["abab"]("a,c,k,j");
        reused_["224_abba_vooo"]("a,i,j,k")  = t2_1["abab"]("a,d,k,j") * f["bb_ov"]("i,d");
        reused_["224_abba_vooo"]("a,i,j,k") += eri["abab_vovo"]("a,i,b,j") * t1_1["aa"]("b,k");
        reused_["224_abba_vooo"]("a,i,j,k") += eri["abab_vovv"]("a,i,b,c") * t2_1["abab"]("b,c,k,j");
        reused_["225_aabb_vooo"]("a,i,j,k")  = t2["abab"]("a,b,i,j") * f["bb_ov"]("k,b");
        reused_["225_aabb_vooo"]("a,i,j,k") += eri["abab_vovv"]("a,k,c,d") * t2["abab"]("c,d,i,j");
        reused_["226_baab_oovo"]("i,j,a,k")  = -1.00 * reused_["6_aabb_voov"]("a,j,i,c") * t1_1["bb"]("c,k");
        reused_["226_baab_oovo"]("i,j,a,k") += eri["abba_oovo"]("l,i,b,j") * t2_1["abab"]("a,b,l,k");
        reused_["227_aabb_vooo"]("a,i,j,k")  = t2["abab"]("a,b,i,j") * reused_["8_bb_ov"]("k,b");
        reused_["227_aabb_vooo"]("a,i,j,k") -= eri["abab_oovo"]("l,k,d,j") * t2_1["aaaa"]("d,a,i,l");
        reused_["227_aabb_vooo"]("a,i,j,k") -= reused_["191_abba_voov"]("a,j,k,c") * t1_1["aa"]("c,i");
        reused_["228_bbbb_vooo"]("a,i,j,k")  = reused_["26_bbbb_voov"]("a,i,j,b") * t1_1["bb"]("b,k");
        reused_["228_bbbb_vooo"]("a,i,j,k") -= eri["abab_oovo"]("l,j,c,i") * t2_1["abab"]("c,a,l,k");
        reused_["229_aa_vo"]("a,i")  = -1.00 * reused_["156_aa_vo"]("a,i");
        reused_["230_aabb_vvvo"]("a,b,c,i")  = eri["aaaa_vovv"]("a,j,d,b") * t2["abab"]("d,c,j,i");
        reused_["231_bbaa_vvvo"]("a,b,c,i")  = eri["bbbb_vovv"]("a,j,d,b") * t2["abab"]("c,d,i,j");
        reused_["232_bbaa_oovo"]("i,j,a,k")  = eri["bbbb_oovo"]("l,i,b,j") * t2["abab"]("a,b,k,l");
        reused_["233_aabb_oovo"]("i,j,a,k")  = eri["aaaa_oovo"]("l,i,b,j") * t2["abab"]("b,a,l,k");
        reused_["234_bb_vo"]("a,i")  = -1.00 * t1_1["bb"]("a,i");
        reused_["235_aabb_vvvo"]("a,b,c,i")  = eri["aaaa_vovv"]("a,j,d,b") * t2_1["abab"]("d,c,j,i");
        reused_["236_bbaa_vvvo"]("a,b,c,i")  = eri["bbbb_vovv"]("a,j,d,b") * t2_1["abab"]("c,d,i,j");
        reused_["237_baab_vvoo"]("a,b,i,j")  = eri["bbbb_oovv"]("k,l,c,d") * t2["bbbb"]("d,a,l,k") * t2["abab"]("b,c,i,j");
        reused_["238_aabb_oovo"]("i,j,a,k")  = -1.00 * eri["aaaa_oovv"]("l,i,c,b") * t1_1["aa"]("b,l") * t2["abab"]("c,a,j,k");
        reused_["238_aabb_oovo"]("i,j,a,k") += eri["aaaa_oovo"]("l,i,c,j") * t2_1["abab"]("c,a,l,k");
        reused_["238_aabb_oovo"]("i,j,a,k") += reused_["27_bbaa_voov"]("a,k,i,b") * t1_1["aa"]("b,j");
        reused_["239_aabb_vooo"]("a,i,j,k")  = reused_["56_aabb_voov"]("a,i,j,b") * t1_1["bb"]("b,k");
        reused_["239_aabb_vooo"]("a,i,j,k") -= eri["bbbb_oovv"]("l,j,c,b") * t1_1["bb"]("b,l") * t2["abab"]("a,c,i,k");
        reused_["239_aabb_vooo"]("a,i,j,k") += eri["bbbb_oovo"]("l,j,c,k") * t2_1["abab"]("a,c,i,l");
        reused_["240_abab_vvoo"]("a,b,i,j")  = reused_["57_aa_vv"]("c,a") * t2["abab"]("c,b,i,j");
        reused_["241_baab_voov"]("a,i,j,b")  = t2_1["abab"]("c,a,i,k") * eri["abab_oovv"]("j,k,c,b");
        reused_["242_abab_vovo"]("a,i,b,j")  = eri["abab_vovv"]("a,i,b,c") * t1_1["bb"]("c,j");
        reused_["243_baba_vovo"]("a,i,b,j")  = eri["baab_vovv"]("a,i,c,b") * t1_1["aa"]("c,j");
        reused_["244_abab_oooo"]("i,j,k,l")  = eri["abba_oovo"]("i,j,a,k") * t1_1["bb"]("a,l");
        reused_["245_abab_oooo"]("i,j,k,l")  = t2_1["abab"]("a,b,i,j") * eri["abab_oovv"]("k,l,a,b");
        reused_["245_abab_oooo"]("i,j,k,l") += eri["abab_oovo"]("k,l,a,j") * t1_1["aa"]("a,i");
        reused_["246_abab_oovv"]("i,j,a,b")  = -0.50 * eri["abab_oovv"]("i,j,a,b");
        reused_["247_bbbb_vvvo"]("a,b,c,i")  = eri["bbbb_vovv"]("a,j,d,b") * t2_1["bbbb"]("d,c,i,j");
        reused_["248_bbbb_vooo"]("a,i,j,k")  = reused_["28_bbbb_voov"]("a,i,j,b") * t1_1["bb"]("b,k");
        reused_["248_bbbb_vooo"]("a,i,j,k") -= eri["bbbb_oovo"]("l,j,c,i") * t2_1["bbbb"]("c,a,k,l");
    }

    world_.gop.fence();
}