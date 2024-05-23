
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

#include "../../../include/derived/eom_ee_qed_ccsd.h"

void hilbert::EOM_EE_QED_CCSD::build_common_ops() {
    
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

        // extract 1-body amplitudes
    std::map<std::string, TA::TArrayD> t1 {
            {"aa_vo", cc_wfn_->amplitudes_["t1_aa"]},
            {"bb_vo", cc_wfn_->amplitudes_["t1_bb"]}
    };
    std::map<std::string, TA::TArrayD> u1 {
            {"aa_vo", cc_wfn_->amplitudes_["u1_aa"]},
            {"bb_vo", cc_wfn_->amplitudes_["u1_bb"]}
    };

    // extract 2-body amplitudes
    std::map<std::string, TA::TArrayD> t2 {
            {"aaaa_vvoo", cc_wfn_->amplitudes_["t2_aaaa"]},
            {"abab_vvoo", cc_wfn_->amplitudes_["t2_abab"]},
            {"bbbb_vvoo", cc_wfn_->amplitudes_["t2_bbbb"]}
    };
    std::map<std::string, TA::TArrayD> u2 {
            {"aaaa_vvoo", cc_wfn_->amplitudes_["u2_aaaa"]},
            {"abab_vvoo", cc_wfn_->amplitudes_["u2_abab"]},
            {"bbbb_vvoo", cc_wfn_->amplitudes_["u2_bbbb"]}
    };

    world_.gop.fence();

    TArrayMap &eri = cc_wfn_->V_blks_;
    TArrayMap &f = cc_wfn_->F_blks_;
    TArrayMap &Id = cc_wfn_->Id_blks_;

    map<string, bool> includes_ = {
        {"u0", include_u0_},
        {"u1", include_u1_},
        {"u2", include_u2_}
    };

    {
        scalars_["0"]  = 1.00 * dot(Id["aa_oo"]("i,i1"), f["aa_oo"]("i,i1"));
        scalars_["1"]  = 1.00 * dot(Id["bb_oo"]("i,i1"), f["bb_oo"]("i,i1"));

        if (includes_["u1"]) {
            scalars_["2"]  = 1.00 * dot(f["aa_ov"]("i,a"), u1["aa_vo"]("a,i"));
            scalars_["3"]  = 1.00 * dot(f["bb_ov"]("i,a"), u1["bb_vo"]("a,i"));
        }
        scalars_["4"]  = 1.00 * dot(Id["aa_oo"]("j,j1") * eri["aaaa_oooo"]("j,i,j1,i1"), Id["aa_oo"]("i,i1"));
        scalars_["5"]  = 1.00 * dot(Id["bb_oo"]("i,i1") * eri["abab_oooo"]("j,i,j1,i1"), Id["aa_oo"]("j,j1"));
        scalars_["6"]  = 1.00 * dot(Id["bb_oo"]("j,j1") * eri["bbbb_oooo"]("j,i,j1,i1"), Id["bb_oo"]("i,i1"));
        scalars_["7"]  = 1.00 * dot(eri["aaaa_oovv"]("j,i,a,b"), t2["aaaa_vvoo"]("a,b,j,i"));
        scalars_["8"]  = 1.00 * dot(eri["abab_oovv"]("j,i,a,b"), t2["abab_vvoo"]("a,b,j,i"));
        scalars_["9"]  = 1.00 * dot(eri["bbbb_oovv"]("j,i,a,b"), t2["bbbb_vvoo"]("a,b,j,i"));

        if (includes_["u2"]) {
            scalars_["10"]  = 1.00 * dot(eri["aaaa_oovv"]("j,i,a,b"), u2["aaaa_vvoo"]("a,b,j,i"));
            scalars_["11"]  = 1.00 * dot(eri["abab_oovv"]("j,i,a,b"), u2["abab_vvoo"]("a,b,j,i"));
            scalars_["12"]  = 1.00 * dot(eri["bbbb_oovv"]("j,i,a,b"), u2["bbbb_vvoo"]("a,b,j,i"));
        }
        scalars_["13"]  = 1.00 * dot(Id["aa_oo"]("i,i1"), dp["aa_oo"]("i,i1"));
        scalars_["14"]  = 1.00 * dot(Id["bb_oo"]("i,i1"), dp["bb_oo"]("i,i1"));

        if (includes_["u1"]) {
            scalars_["15"]  = 1.00 * dot(dp["aa_ov"]("i,a"), u1["aa_vo"]("a,i"));
            scalars_["16"]  = 1.00 * dot(dp["bb_ov"]("i,a"), u1["bb_vo"]("a,i"));
        }
    }

    {
        reuse_tmps_["1_aabb_vvvv"]("b,e,a,f")  = eri["abab_oovv"]("j,i,e,f") * t2["abab_vvoo"]("b,a,j,i");

        if (includes_["u1"]) {
            reuse_tmps_["262_abba_vvvo"]("e,b,f,m")  = u1["aa_vo"]("a,m") * reuse_tmps_["1_aabb_vvvv"]("e,a,f,b");
            reuse_tmps_["263_aabb_vvvo"]("e,b,a,i")  = u1["bb_vo"]("c,i") * reuse_tmps_["1_aabb_vvvv"]("b,e,a,c");
        }
        reuse_tmps_["2_bbaa_vvoo"]("f,a,m,i")  = eri["abab_oovv"]("i,j,b,a") * t2["abab_vvoo"]("b,f,m,j");

        if (includes_["u1"]) {
            reuse_tmps_["269_baab_vooo"]("f,j,m,n")  = reuse_tmps_["2_bbaa_vvoo"]("f,b,m,j") * u1["bb_vo"]("b,n");
        }
        reuse_tmps_["311_abab_vvoo"]("a,b,j,i")  = reuse_tmps_["2_bbaa_vvoo"]("b,d,j,l") * t2["abab_vvoo"]("a,d,l,i");
        reuse_tmps_["3_bbbb_vvoo"]("f,a,n,i")  = eri["abab_oovv"]("j,i,b,a") * t2["abab_vvoo"]("b,f,j,n");

        if (includes_["u1"]) {
            reuse_tmps_["264_bbbb_vooo"]("e,m,j,n")  = reuse_tmps_["3_bbbb_vvoo"]("e,b,n,j") * u1["bb_vo"]("b,m");
            reuse_tmps_["273_bb_vo"]("a,i")  = u1["bb_vo"]("c,k") * reuse_tmps_["3_bbbb_vvoo"]("a,c,i,k");
        }
        reuse_tmps_["303_bbbb_vvoo"]("f,e,m,n")  = t2["bbbb_vvoo"]("a,e,n,i") * reuse_tmps_["3_bbbb_vvoo"]("f,a,m,i");
        reuse_tmps_["304_abab_vvoo"]("e,f,m,n")  = t2["abab_vvoo"]("e,a,m,i") * reuse_tmps_["3_bbbb_vvoo"]("f,a,n,i");
        reuse_tmps_["4_aaaa_vvoo"]("f,b,m,j")  = eri["abab_oovv"]("j,i,b,a") * t2["abab_vvoo"]("f,a,m,i");

        if (includes_["u1"]) {
            reuse_tmps_["266_aaaa_vooo"]("e,m,j,n")  = reuse_tmps_["4_aaaa_vvoo"]("e,b,n,j") * u1["aa_vo"]("b,m");
            reuse_tmps_["272_aa_vo"]("f,m")  = u1["aa_vo"]("b,j") * reuse_tmps_["4_aaaa_vvoo"]("f,b,m,j");
        }
        reuse_tmps_["302_aaaa_vvoo"]("f,e,n,m")  = t2["aaaa_vvoo"]("a,e,m,i") * reuse_tmps_["4_aaaa_vvoo"]("f,a,n,i");

        if (includes_["u2"]) {
            reuse_tmps_["5_aaaa_vvoo"]("f,b,m,j")  = eri["abab_oovv"]("j,i,b,a") * u2["abab_vvoo"]("f,a,m,i");
        }
        reuse_tmps_["6_bbbb_vvvv"]("a,b,f,e")  = 0.50 * eri["bbbb_oovv"]("j,i,e,f") * t2["bbbb_vvoo"]("b,a,j,i");

        if (includes_["u1"]) {
            reuse_tmps_["260_bbbb_vvvo"]("e,b,a,i")  = u1["bb_vo"]("c,i") * reuse_tmps_["6_bbbb_vvvv"]("a,b,e,c");
        }
        reuse_tmps_["7_aaaa_vvvv"]("a,b,f,e")  = 0.50 * eri["aaaa_oovv"]("j,i,e,f") * t2["aaaa_vvoo"]("b,a,j,i");

        if (includes_["u1"]) {
            reuse_tmps_["261_aaaa_vvvo"]("e,b,a,i")  = u1["aa_vo"]("c,i") * reuse_tmps_["7_aaaa_vvvv"]("a,b,e,c");
        }
        reuse_tmps_["8_aaaa_vvoo"]("f,a,m,i")  = eri["aaaa_oovv"]("j,i,a,b") * t2["aaaa_vvoo"]("b,f,m,j");
        reuse_tmps_["306_aaaa_vvoo"]("f,e,n,m")  = t2["aaaa_vvoo"]("a,e,m,i") * reuse_tmps_["8_aaaa_vvoo"]("f,a,n,i");
        reuse_tmps_["313_abab_vvoo"]("a,b,i,j")  = t2["abab_vvoo"]("c,b,k,j") * reuse_tmps_["8_aaaa_vvoo"]("a,c,i,k");
        reuse_tmps_["9_bbbb_vvoo"]("f,a,n,i")  = eri["bbbb_oovv"]("j,i,a,b") * t2["bbbb_vvoo"]("b,f,n,j");
        reuse_tmps_["305_bbbb_vvoo"]("f,e,n,m")  = t2["bbbb_vvoo"]("a,e,m,i") * reuse_tmps_["9_bbbb_vvoo"]("f,a,n,i");
        reuse_tmps_["308_abab_vvoo"]("e,f,m,n")  = t2["abab_vvoo"]("e,a,m,i") * reuse_tmps_["9_bbbb_vvoo"]("f,a,n,i");
        reuse_tmps_["10_aaaa_vvoo"]("a,e,i,m")  = eri["aaaa_oovv"]("m,j,b,e") * t2["aaaa_vvoo"]("b,a,i,j");

        if (includes_["u1"]) {
            reuse_tmps_["268_aaaa_vooo"]("f,m,j,n")  = reuse_tmps_["10_aaaa_vvoo"]("f,b,n,j") * u1["aa_vo"]("b,m");
            reuse_tmps_["280_aa_vo"]("f,m")  = u1["aa_vo"]("b,j") * reuse_tmps_["10_aaaa_vvoo"]("f,b,m,j");
        }
        reuse_tmps_["307_abab_vvoo"]("e,f,m,n")  = t2["abab_vvoo"]("b,f,j,n") * reuse_tmps_["10_aaaa_vvoo"]("e,b,m,j");
        reuse_tmps_["11_bbbb_vvoo"]("b,e,i,k")  = eri["bbbb_oovv"]("k,j,c,e") * t2["bbbb_vvoo"]("c,b,i,j");

        if (includes_["u1"]) {
            reuse_tmps_["267_bbbb_vooo"]("f,m,j,n")  = reuse_tmps_["11_bbbb_vvoo"]("f,b,n,j") * u1["bb_vo"]("b,m");
            reuse_tmps_["281_bb_vo"]("f,n")  = u1["bb_vo"]("b,j") * reuse_tmps_["11_bbbb_vvoo"]("f,b,n,j");
        }
        reuse_tmps_["314_abab_vvoo"]("a,b,i,j")  = reuse_tmps_["11_bbbb_vvoo"]("b,d,j,l") * t2["abab_vvoo"]("a,d,i,l");
        reuse_tmps_["12_abab_vvoo"]("f,a,m,i")  = eri["abab_oovv"]("j,i,b,a") * t2["aaaa_vvoo"]("b,f,m,j");

        if (includes_["u1"]) {
            reuse_tmps_["289_aabb_vooo"]("e,m,j,n")  = u1["bb_vo"]("b,n") * reuse_tmps_["12_abab_vvoo"]("e,b,m,j");
            reuse_tmps_["299_aa_vo"]("f,m")  = u1["bb_vo"]("b,j") * reuse_tmps_["12_abab_vvoo"]("f,b,m,j");
        }
        reuse_tmps_["312_abab_vvoo"]("e,f,m,n")  = t2["bbbb_vvoo"]("b,f,n,j") * reuse_tmps_["12_abab_vvoo"]("e,b,m,j");
        reuse_tmps_["13_abab_vvoo"]("f,a,n,i")  = eri["abab_oovv"]("n,j,f,b") * t2["bbbb_vvoo"]("b,a,i,j");

        if (includes_["u1"]) {
            reuse_tmps_["290_baab_vooo"]("f,m,j,n")  = reuse_tmps_["13_abab_vvoo"]("b,f,j,n") * u1["aa_vo"]("b,m");
            reuse_tmps_["301_bb_vo"]("a,i")  = u1["aa_vo"]("c,k") * reuse_tmps_["13_abab_vvoo"]("c,a,k,i");
        }
        reuse_tmps_["14_abab_vvoo"]("a,f,i,n")  = eri["aaaa_oovv"]("j,i,a,b") * t2["abab_vvoo"]("b,f,j,n");
        reuse_tmps_["309_bbbb_vvoo"]("e,f,m,n")  = reuse_tmps_["14_abab_vvoo"]("a,f,i,n") * t2["abab_vvoo"]("a,e,i,m");
        reuse_tmps_["15_abab_vvoo"]("f,a,m,i")  = eri["bbbb_oovv"]("j,i,a,b") * t2["abab_vvoo"]("f,b,m,j");
        reuse_tmps_["310_aaaa_vvoo"]("a,b,i,j")  = t2["abab_vvoo"]("b,c,j,k") * reuse_tmps_["15_abab_vvoo"]("a,c,i,k");
        reuse_tmps_["16_abab_vvoo"]("f,b,m,j")  = eri["bbbb_oovv"]("j,i,a,b") * t2["abab_vvoo"]("f,a,m,i");

        if (includes_["u1"]) {
            reuse_tmps_["288_aabb_vooo"]("e,m,j,n")  = u1["bb_vo"]("b,n") * reuse_tmps_["16_abab_vvoo"]("e,b,m,j");
            reuse_tmps_["298_aa_vo"]("f,m")  = u1["bb_vo"]("b,j") * reuse_tmps_["16_abab_vvoo"]("f,b,m,j");
        }
        reuse_tmps_["17_abab_vvoo"]("b,e,j,m")  = eri["aaaa_oovv"]("j,i,a,b") * t2["abab_vvoo"]("a,e,i,m");

        if (includes_["u1"]) {
            reuse_tmps_["291_baab_vooo"]("f,m,j,n")  = reuse_tmps_["17_abab_vvoo"]("b,f,j,n") * u1["aa_vo"]("b,m");
            reuse_tmps_["300_bb_vo"]("a,i")  = reuse_tmps_["17_abab_vvoo"]("c,a,k,i") * u1["aa_vo"]("c,k");
        }

        if (includes_["u2"]) {
            reuse_tmps_["18_bbbb_vvoo"]("e,b,m,j")  = eri["abab_oovv"]("i,j,a,b") * u2["abab_vvoo"]("a,e,i,m");
        }

        if (includes_["u1"]) {
            reuse_tmps_["19_bbbb_vvoo"]("e,a,i,m")  = eri["bbbb_vovv"]("a,m,b,e") * u1["bb_vo"]("b,i");
            reuse_tmps_["20_aaaa_vvoo"]("e,a,i,m")  = eri["aaaa_vovv"]("a,m,b,e") * u1["aa_vo"]("b,i");
        }

        if (includes_["u2"]) {
            reuse_tmps_["21_bbbb_vvoo"]("a,e,i,m")  = eri["bbbb_oovv"]("m,j,b,e") * u2["bbbb_vvoo"]("b,a,i,j");
            reuse_tmps_["22_aaaa_vvoo"]("a,e,i,m")  = eri["aaaa_oovv"]("m,j,b,e") * u2["aaaa_vvoo"]("b,a,i,j");
            reuse_tmps_["23_abab_vvoo"]("f,b,m,j")  = eri["bbbb_oovv"]("j,i,a,b") * u2["abab_vvoo"]("f,a,m,i");
            reuse_tmps_["24_abab_vvoo"]("b,e,j,m")  = eri["aaaa_oovv"]("j,i,a,b") * u2["abab_vvoo"]("a,e,i,m");
            reuse_tmps_["25_abab_vvoo"]("b,e,j,m")  = eri["abab_oovv"]("j,i,b,a") * u2["bbbb_vvoo"]("a,e,m,i");
            reuse_tmps_["26_abab_vvoo"]("e,b,m,j")  = eri["abab_oovv"]("i,j,a,b") * u2["aaaa_vvoo"]("a,e,m,i");
        }

        if (includes_["u1"]) {
            reuse_tmps_["27_abab_vvoo"]("e,b,m,i")  = eri["abab_vovv"]("e,i,a,b") * u1["aa_vo"]("a,m");
            reuse_tmps_["28_abab_vvoo"]("f,a,n,i")  = eri["baab_vovv"]("a,n,f,b") * u1["bb_vo"]("b,i");
        }
        reuse_tmps_["29_aabb_vvoo"]("a,c,i,k")  = eri["abab_oovv"]("l,k,c,d") * t2["abab_vvoo"]("a,d,l,i");

        if (includes_["u1"]) {
            reuse_tmps_["265_aabb_vooo"]("e,m,j,n")  = u1["aa_vo"]("b,m") * reuse_tmps_["29_aabb_vvoo"]("e,b,n,j");
        }

        if (includes_["u2"]) {
            reuse_tmps_["30_aabb_vvoo"]("a,e,i,n")  = eri["abab_oovv"]("j,n,e,b") * u2["abab_vvoo"]("a,b,j,i");
        }

        if (includes_["u1"]) {
            reuse_tmps_["31_aabb_vvoo"]("e,a,i,n")  = eri["abab_vovv"]("a,n,e,b") * u1["bb_vo"]("b,i");
            reuse_tmps_["32_bbaa_vvoo"]("f,a,i,m")  = eri["baab_vovv"]("a,m,b,f") * u1["aa_vo"]("b,i");
            reuse_tmps_["33_bb_vv"]("e,a")  = eri["bbbb_vovv"]("a,i,b,e") * u1["bb_vo"]("b,i");
            reuse_tmps_["34_aa_vv"]("e,a")  = eri["aaaa_vovv"]("a,i,b,e") * u1["aa_vo"]("b,i");
        }

        if (includes_["u2"]) {
            reuse_tmps_["35_bbaa_vvoo"]("a,f,i,m")  = eri["abab_oovv"]("m,j,b,f") * u2["abab_vvoo"]("b,a,i,j");
        }

        if (includes_["u1"]) {
            reuse_tmps_["36_bb_vv"]("e,a")  = eri["baab_vovv"]("a,i,b,e") * u1["aa_vo"]("b,i");
            reuse_tmps_["37_aa_vv"]("e,a")  = eri["abab_vovv"]("a,i,e,b") * u1["bb_vo"]("b,i");
            reuse_tmps_["38_bbbb_oooo"]("i,j,n,m")  = 0.50 * eri["bbbb_oovo"]("m,n,a,j") * u1["bb_vo"]("a,i");
            reuse_tmps_["39_aaaa_oooo"]("i,j,n,m")  = 0.50 * eri["aaaa_oovo"]("m,n,a,j") * u1["aa_vo"]("a,i");
        }
        reuse_tmps_["40_aabb_oooo"]("m,i,n,j")  = eri["abab_oovv"]("i,j,a,b") * t2["abab_vvoo"]("a,b,m,n");

        if (includes_["u1"]) {
            reuse_tmps_["270_aabb_vooo"]("e,m,j,n")  = u1["aa_vo"]("e,i") * reuse_tmps_["40_aabb_oooo"]("m,i,n,j");
            reuse_tmps_["271_baab_vooo"]("a,m,i,j")  = reuse_tmps_["40_aabb_oooo"]("i,m,j,k") * u1["bb_vo"]("a,k");
        }
        reuse_tmps_["315_abab_vvoo"]("b,a,i,j")  = t2["abab_vvoo"]("b,a,l,k") * reuse_tmps_["40_aabb_oooo"]("i,l,j,k");
        reuse_tmps_["41_aaaa_oooo"]("j,i,k,l")  = 0.25 * eri["aaaa_oovv"]("l,k,c,d") * t2["aaaa_vvoo"]("c,d,i,j");

        if (includes_["u1"]) {
            reuse_tmps_["279_aaaa_vooo"]("a,m,i,j")  = u1["aa_vo"]("a,k") * reuse_tmps_["41_aaaa_oooo"]("j,i,k,m");
        }
        reuse_tmps_["316_aaaa_vvoo"]("a,b,i,j")  = t2["aaaa_vvoo"]("b,a,l,k") * reuse_tmps_["41_aaaa_oooo"]("j,i,k,l");
        reuse_tmps_["42_bbbb_oooo"]("j,i,k,l")  = 0.25 * eri["bbbb_oovv"]("l,k,c,d") * t2["bbbb_vvoo"]("c,d,i,j");

        if (includes_["u1"]) {
            reuse_tmps_["278_bbbb_vooo"]("a,m,i,j")  = reuse_tmps_["42_bbbb_oooo"]("j,i,k,m") * u1["bb_vo"]("a,k");
        }
        reuse_tmps_["317_bbbb_vvoo"]("a,b,i,j")  = t2["bbbb_vvoo"]("b,a,l,k") * reuse_tmps_["42_bbbb_oooo"]("j,i,k,l");

        if (includes_["u1"]) {
            reuse_tmps_["43_bb_oo"]("i,m")  = eri["bbbb_oovo"]("m,j,a,i") * u1["bb_vo"]("a,j");
            reuse_tmps_["44_aa_oo"]("i,m")  = eri["aaaa_oovo"]("m,j,a,i") * u1["aa_vo"]("a,j");
        }
        reuse_tmps_["45_bb_vv"]("f,a")  = eri["abab_oovv"]("j,i,b,a") * t2["abab_vvoo"]("b,f,j,i");

        if (includes_["u1"]) {
            reuse_tmps_["274_bb_vo"]("a,i")  = reuse_tmps_["45_bb_vv"]("a,c") * u1["bb_vo"]("c,i");
        }
        reuse_tmps_["322_abab_vvoo"]("e,f,m,n")  = t2["abab_vvoo"]("e,a,m,n") * reuse_tmps_["45_bb_vv"]("f,a");
        reuse_tmps_["323_bbbb_vvoo"]("e,f,n,m")  = t2["bbbb_vvoo"]("b,f,m,n") * reuse_tmps_["45_bb_vv"]("e,b");
        reuse_tmps_["46_aa_vv"]("f,b")  = eri["abab_oovv"]("j,i,b,a") * t2["abab_vvoo"]("f,a,j,i");

        if (includes_["u1"]) {
            reuse_tmps_["275_aa_vo"]("a,i")  = u1["aa_vo"]("c,i") * reuse_tmps_["46_aa_vv"]("a,c");
        }
        reuse_tmps_["326_aaaa_vvoo"]("a,b,j,i")  = reuse_tmps_["46_aa_vv"]("b,d") * t2["aaaa_vvoo"]("d,a,i,j");
        reuse_tmps_["327_abab_vvoo"]("b,a,i,j")  = reuse_tmps_["46_aa_vv"]("b,d") * t2["abab_vvoo"]("d,a,i,j");
        reuse_tmps_["47_bb_vv"]("f,a")  = 0.50 * eri["bbbb_oovv"]("j,i,a,b") * t2["bbbb_vvoo"]("b,f,j,i");
        reuse_tmps_["324_bbbb_vvoo"]("f,e,n,m")  = t2["bbbb_vvoo"]("a,e,m,n") * reuse_tmps_["47_bb_vv"]("f,a");
        reuse_tmps_["334_abab_vvoo"]("e,f,m,n")  = t2["abab_vvoo"]("e,a,m,n") * reuse_tmps_["47_bb_vv"]("f,a");
        reuse_tmps_["48_aa_vv"]("f,a")  = eri["aaaa_oovv"]("j,i,a,b") * t2["aaaa_vvoo"]("b,f,j,i");
        reuse_tmps_["333_aaaa_vvoo"]("e,f,n,m")  = reuse_tmps_["48_aa_vv"]("f,a") * t2["aaaa_vvoo"]("a,e,m,n");
        reuse_tmps_["49_aa_vv"]("a,e")  = 0.50 * eri["aaaa_oovv"]("j,i,b,e") * t2["aaaa_vvoo"]("b,a,j,i");

        if (includes_["u1"]) {
            reuse_tmps_["282_aa_vo"]("a,i")  = u1["aa_vo"]("c,i") * reuse_tmps_["49_aa_vv"]("a,c");
        }
        reuse_tmps_["325_aaaa_vvoo"]("f,e,n,m")  = reuse_tmps_["49_aa_vv"]("e,b") * t2["aaaa_vvoo"]("b,f,m,n");
        reuse_tmps_["335_abab_vvoo"]("b,a,i,j")  = t2["abab_vvoo"]("d,a,i,j") * reuse_tmps_["49_aa_vv"]("b,d");
        reuse_tmps_["50_bb_vv"]("a,c")  = 0.50 * eri["bbbb_oovv"]("k,j,b,c") * t2["bbbb_vvoo"]("b,a,k,j");

        if (includes_["u1"]) {
            reuse_tmps_["283_bb_vo"]("a,i")  = u1["bb_vo"]("c,i") * reuse_tmps_["50_bb_vv"]("a,c");
        }
        reuse_tmps_["332_bbbb_vvoo"]("e,f,n,m")  = t2["bbbb_vvoo"]("b,f,m,n") * reuse_tmps_["50_bb_vv"]("e,b");
        reuse_tmps_["345_abab_vvoo"]("a,b,i,j")  = t2["abab_vvoo"]("a,d,i,j") * reuse_tmps_["50_bb_vv"]("b,d");
        reuse_tmps_["51_aa_oo"]("i,m")  = eri["abab_oovv"]("m,j,a,b") * t2["abab_vvoo"]("a,b,i,j");

        if (includes_["u1"]) {
            reuse_tmps_["277_aa_vo"]("a,i")  = u1["aa_vo"]("a,k") * reuse_tmps_["51_aa_oo"]("i,k");
        }
        reuse_tmps_["340_aaaa_vvoo"]("a,b,j,i")  = t2["aaaa_vvoo"]("b,a,i,l") * reuse_tmps_["51_aa_oo"]("j,l");
        reuse_tmps_["347_abab_vvoo"]("e,f,m,n")  = reuse_tmps_["51_aa_oo"]("m,j") * t2["abab_vvoo"]("e,f,j,n");
        reuse_tmps_["52_bb_oo"]("i,k")  = eri["abab_oovv"]("j,k,b,c") * t2["abab_vvoo"]("b,c,j,i");

        if (includes_["u1"]) {
            reuse_tmps_["276_bb_vo"]("a,i")  = u1["bb_vo"]("a,k") * reuse_tmps_["52_bb_oo"]("i,k");
        }
        reuse_tmps_["336_abab_vvoo"]("e,f,m,n")  = reuse_tmps_["52_bb_oo"]("n,j") * t2["abab_vvoo"]("e,f,m,j");
        reuse_tmps_["343_bbbb_vvoo"]("a,b,j,i")  = t2["bbbb_vvoo"]("b,a,i,l") * reuse_tmps_["52_bb_oo"]("j,l");
        reuse_tmps_["53_aa_oo"]("i,m")  = 0.50 * eri["aaaa_oovv"]("m,j,a,b") * t2["aaaa_vvoo"]("a,b,i,j");

        if (includes_["u1"]) {
            reuse_tmps_["285_aa_vo"]("a,i")  = u1["aa_vo"]("a,k") * reuse_tmps_["53_aa_oo"]("i,k");
        }
        reuse_tmps_["342_aaaa_vvoo"]("a,b,i,j")  = reuse_tmps_["53_aa_oo"]("j,l") * t2["aaaa_vvoo"]("b,a,i,l");
        reuse_tmps_["346_abab_vvoo"]("e,f,m,n")  = reuse_tmps_["53_aa_oo"]("m,j") * t2["abab_vvoo"]("e,f,j,n");
        reuse_tmps_["54_bb_oo"]("i,k")  = 0.50 * eri["bbbb_oovv"]("k,j,b,c") * t2["bbbb_vvoo"]("b,c,i,j");

        if (includes_["u1"]) {
            reuse_tmps_["284_bb_vo"]("a,i")  = u1["bb_vo"]("a,k") * reuse_tmps_["54_bb_oo"]("i,k");
        }
        reuse_tmps_["341_bbbb_vvoo"]("a,b,j,i")  = t2["bbbb_vvoo"]("b,a,i,l") * reuse_tmps_["54_bb_oo"]("j,l");
        reuse_tmps_["350_abab_vvoo"]("b,a,i,j")  = reuse_tmps_["54_bb_oo"]("j,l") * t2["abab_vvoo"]("b,a,i,l");

        if (includes_["u2"]) {
            reuse_tmps_["55_aabb_oooo"]("m,j,n,i")  = eri["abab_oovv"]("j,i,a,b") * u2["abab_vvoo"]("a,b,m,n");
            reuse_tmps_["56_bbbb_oooo"]("j,i,n,m")  = 0.25 * eri["bbbb_oovv"]("m,n,a,b") * u2["bbbb_vvoo"]("a,b,i,j");
            reuse_tmps_["57_aaaa_oooo"]("j,i,n,m")  = 0.25 * eri["aaaa_oovv"]("m,n,a,b") * u2["aaaa_vvoo"]("a,b,i,j");
            reuse_tmps_["58_aa_vv"]("f,b")  = eri["abab_oovv"]("j,i,b,a") * u2["abab_vvoo"]("f,a,j,i");
            reuse_tmps_["59_bb_vv"]("e,b")  = eri["abab_oovv"]("j,i,a,b") * u2["abab_vvoo"]("a,e,j,i");
            reuse_tmps_["60_bb_vv"]("a,e")  = 0.50 * eri["bbbb_oovv"]("j,i,b,e") * u2["bbbb_vvoo"]("b,a,j,i");
            reuse_tmps_["61_aa_vv"]("a,e")  = 0.50 * eri["aaaa_oovv"]("j,i,b,e") * u2["aaaa_vvoo"]("b,a,j,i");
            reuse_tmps_["62_bb_oo"]("i,m")  = eri["abab_oovv"]("j,m,a,b") * u2["abab_vvoo"]("a,b,j,i");
            reuse_tmps_["63_bb_oo"]("i,m")  = 0.50 * eri["bbbb_oovv"]("m,j,a,b") * u2["bbbb_vvoo"]("a,b,i,j");
            reuse_tmps_["64_aa_oo"]("i,m")  = 0.50 * eri["aaaa_oovv"]("m,j,a,b") * u2["aaaa_vvoo"]("a,b,i,j");
            reuse_tmps_["65_aa_oo"]("i,m")  = eri["abab_oovv"]("m,j,a,b") * u2["abab_vvoo"]("a,b,i,j");
        }

        if (includes_["u1"]) {
            reuse_tmps_["66_aabb_oooo"]("m,j,n,i")  = eri["abba_oovo"]("j,i,a,m") * u1["bb_vo"]("a,n");
            reuse_tmps_["67_aabb_oooo"]("m,j,n,i")  = eri["abab_oovo"]("j,i,a,n") * u1["aa_vo"]("a,m");
            reuse_tmps_["68_aa_oo"]("m,j")  = eri["abba_oovo"]("j,i,a,m") * u1["bb_vo"]("a,i");
            reuse_tmps_["69_bb_oo"]("i,m")  = eri["abab_oovo"]("j,m,a,i") * u1["aa_vo"]("a,j");
        }
        reuse_tmps_["70_bbbb_vvvo"]("e,b,f,n")  = eri["baab_vovv"]("f,i,a,b") * t2["abab_vvoo"]("a,e,i,n");
        reuse_tmps_["71_abba_vvvo"]("e,f,b,m")  = eri["abab_vovv"]("e,i,a,b") * t2["abab_vvoo"]("a,f,m,i");
        reuse_tmps_["72_bbbb_vvvo"]("e,b,f,m")  = eri["bbbb_vovv"]("f,i,a,b") * t2["bbbb_vvoo"]("a,e,m,i");
        reuse_tmps_["73_aaaa_vvvo"]("f,b,e,m")  = eri["aaaa_vovv"]("e,i,a,b") * t2["aaaa_vvoo"]("a,f,m,i");
        reuse_tmps_["74_aaaa_vvvo"]("f,b,e,m")  = eri["abab_vovv"]("e,i,b,a") * t2["abab_vvoo"]("f,a,m,i");

        if (includes_["u2"]) {
            reuse_tmps_["75_abba_vvvo"]("e,f,b,m")  = eri["abab_vovv"]("e,i,a,b") * u2["abab_vvoo"]("a,f,m,i");
            reuse_tmps_["76_aabb_vvvo"]("e,b,f,n")  = eri["baab_vovv"]("f,i,b,a") * u2["abab_vvoo"]("e,a,i,n");
            reuse_tmps_["77_bbbb_vvvo"]("e,b,f,n")  = eri["bbbb_vovv"]("f,i,a,b") * u2["bbbb_vvoo"]("a,e,n,i");
            reuse_tmps_["78_aaaa_vvvo"]("e,b,f,n")  = eri["aaaa_vovv"]("f,i,a,b") * u2["aaaa_vvoo"]("a,e,n,i");
        }
        reuse_tmps_["79_aabb_vvvo"]("e,b,f,n")  = eri["baab_vovv"]("f,i,b,a") * t2["abab_vvoo"]("e,a,i,n");
        reuse_tmps_["80_aabb_vvvo"]("b,e,f,n")  = eri["aaaa_vovv"]("e,i,a,b") * t2["abab_vvoo"]("a,f,i,n");
        reuse_tmps_["81_abba_vvvo"]("a,e,b,i")  = eri["bbbb_vovv"]("b,j,c,e") * t2["abab_vvoo"]("a,c,i,j");

        if (includes_["u2"]) {
            reuse_tmps_["82_aaaa_vvvo"]("f,b,e,m")  = eri["abab_vovv"]("e,i,b,a") * u2["abab_vvoo"]("f,a,m,i");
            reuse_tmps_["83_bbbb_vvvo"]("e,b,f,m")  = eri["baab_vovv"]("f,i,a,b") * u2["abab_vvoo"]("a,e,i,m");
        }
        reuse_tmps_["84_bbbb_vvvo"]("f,e,a,m")  = 0.50 * eri["bbbb_oovo"]("j,i,a,m") * t2["bbbb_vvoo"]("e,f,j,i");
        reuse_tmps_["85_aaaa_vvvo"]("f,e,a,n")  = 0.50 * eri["aaaa_oovo"]("j,i,a,n") * t2["aaaa_vvoo"]("e,f,j,i");
        reuse_tmps_["86_bbbb_vooo"]("a,j,i,m")  = 0.50 * eri["bbbb_vovv"]("a,m,b,c") * t2["bbbb_vvoo"]("b,c,i,j");
        reuse_tmps_["87_aaaa_vooo"]("a,j,i,m")  = 0.50 * eri["aaaa_vovv"]("a,m,b,c") * t2["aaaa_vvoo"]("b,c,i,j");
        reuse_tmps_["88_bb_vo"]("a,i")  = 0.50 * eri["bbbb_vovv"]("a,j,b,c") * t2["bbbb_vvoo"]("b,c,i,j");
        reuse_tmps_["89_aa_vo"]("a,i")  = 0.50 * eri["aaaa_vovv"]("a,j,b,c") * t2["aaaa_vvoo"]("b,c,i,j");
        reuse_tmps_["90_abba_vvvo"]("e,b,f,m")  = eri["baab_vovv"]("f,i,a,b") * t2["aaaa_vvoo"]("a,e,m,i");
        reuse_tmps_["91_aabb_vvvo"]("b,e,f,n")  = eri["abab_vovv"]("e,i,b,a") * t2["bbbb_vvoo"]("a,f,n,i");

        if (includes_["u2"]) {
            reuse_tmps_["92_abba_vvvo"]("a,e,b,i")  = eri["bbbb_vovv"]("b,j,c,e") * u2["abab_vvoo"]("a,c,i,j");
            reuse_tmps_["93_aabb_vvvo"]("b,e,f,n")  = eri["aaaa_vovv"]("e,i,a,b") * u2["abab_vvoo"]("a,f,i,n");
        }
        reuse_tmps_["94_abba_vvvo"]("e,f,a,m")  = eri["abba_oovo"]("j,i,a,m") * t2["abab_vvoo"]("e,f,j,i");
        reuse_tmps_["95_aabb_vvvo"]("e,a,f,n")  = eri["abab_oovo"]("j,i,a,n") * t2["abab_vvoo"]("e,f,j,i");
        reuse_tmps_["96_aabb_vooo"]("e,m,n,i")  = eri["abab_vovv"]("e,i,a,b") * t2["abab_vvoo"]("a,b,m,n");
        reuse_tmps_["97_baab_vooo"]("f,m,i,n")  = eri["baab_vovv"]("f,i,a,b") * t2["abab_vvoo"]("a,b,m,n");
        reuse_tmps_["98_aa_vo"]("a,i")  = eri["abab_vovv"]("a,j,b,c") * t2["abab_vvoo"]("b,c,i,j");
        reuse_tmps_["99_bb_vo"]("a,i")  = eri["baab_vovv"]("a,j,b,c") * t2["abab_vvoo"]("b,c,j,i");

        if (includes_["u2"]) {
            reuse_tmps_["100_bbbb_vvvo"]("f,e,a,n")  = 0.50 * eri["bbbb_oovo"]("j,i,a,n") * u2["bbbb_vvoo"]("e,f,j,i");
            reuse_tmps_["101_aaaa_vvvo"]("f,e,a,n")  = 0.50 * eri["aaaa_oovo"]("j,i,a,n") * u2["aaaa_vvoo"]("e,f,j,i");
            reuse_tmps_["102_bbbb_vooo"]("a,j,i,m")  = 0.50 * eri["bbbb_vovv"]("a,m,b,c") * u2["bbbb_vvoo"]("b,c,i,j");
            reuse_tmps_["103_aaaa_vooo"]("a,j,i,m")  = 0.50 * eri["aaaa_vovv"]("a,m,b,c") * u2["aaaa_vvoo"]("b,c,i,j");
            reuse_tmps_["104_aa_vo"]("a,i")  = 0.50 * eri["aaaa_vovv"]("a,j,b,c") * u2["aaaa_vvoo"]("b,c,i,j");
            reuse_tmps_["105_bb_vo"]("a,i")  = 0.50 * eri["bbbb_vovv"]("a,j,b,c") * u2["bbbb_vvoo"]("b,c,i,j");
        }

        if (includes_["u1"]) {
            reuse_tmps_["106_bbbb_vvvo"]("b,f,e,m")  = eri["bbbb_vvvv"]("e,f,a,b") * u1["bb_vo"]("a,m");
            reuse_tmps_["107_aaaa_vvvo"]("b,f,e,m")  = eri["aaaa_vvvv"]("e,f,a,b") * u1["aa_vo"]("a,m");
        }

        if (includes_["u2"]) {
            reuse_tmps_["108_abba_vvvo"]("e,b,f,m")  = eri["baab_vovv"]("f,i,a,b") * u2["aaaa_vvoo"]("a,e,m,i");
            reuse_tmps_["109_aabb_vvvo"]("b,e,f,n")  = eri["abab_vovv"]("e,i,b,a") * u2["bbbb_vvoo"]("a,f,n,i");
            reuse_tmps_["110_aabb_vvvo"]("e,a,f,n")  = eri["abab_oovo"]("j,i,a,n") * u2["abab_vvoo"]("e,f,j,i");
            reuse_tmps_["111_abba_vvvo"]("e,f,a,m")  = eri["abba_oovo"]("j,i,a,m") * u2["abab_vvoo"]("e,f,j,i");
            reuse_tmps_["112_aabb_vooo"]("e,m,n,i")  = eri["abab_vovv"]("e,i,a,b") * u2["abab_vvoo"]("a,b,m,n");
            reuse_tmps_["113_baab_vooo"]("f,m,i,n")  = eri["baab_vovv"]("f,i,a,b") * u2["abab_vvoo"]("a,b,m,n");
            reuse_tmps_["114_bb_vo"]("a,i")  = eri["baab_vovv"]("a,j,b,c") * u2["abab_vvoo"]("b,c,j,i");
            reuse_tmps_["115_aa_vo"]("a,i")  = eri["abab_vovv"]("a,j,b,c") * u2["abab_vvoo"]("b,c,i,j");
        }

        if (includes_["u1"]) {
            reuse_tmps_["116_abba_vvvo"]("b,e,a,i")  = eri["abab_vvvv"]("b,a,c,e") * u1["aa_vo"]("c,i");
            reuse_tmps_["117_aabb_vvvo"]("e,b,a,i")  = eri["abab_vvvv"]("b,a,e,c") * u1["bb_vo"]("c,i");
        }
        reuse_tmps_["118_abab_vvoo"]("e,f,m,n")  = eri["abab_vvvv"]("e,f,a,b") * t2["abab_vvoo"]("a,b,m,n");
        reuse_tmps_["119_aaaa_vvoo"]("f,e,n,m")  = eri["aaaa_vvvv"]("e,f,a,b") * t2["aaaa_vvoo"]("a,b,m,n");
        reuse_tmps_["120_bbbb_vvoo"]("a,b,j,i")  = eri["bbbb_vvvv"]("b,a,c,d") * t2["bbbb_vvoo"]("c,d,i,j");

        if (includes_["u1"]) {
            reuse_tmps_["121_aa_oo"]("i,m")  = dp["aa_ov"]("m,a") * u1["aa_vo"]("a,i");
            reuse_tmps_["318_aa_vo"]("a,i")  = u1["aa_vo"]("a,j") * reuse_tmps_["121_aa_oo"]("i,j");
            reuse_tmps_["344_aaaa_vvoo"]("a,b,i,j")  = t2["aaaa_vvoo"]("b,a,j,k") * reuse_tmps_["121_aa_oo"]("i,k");
            reuse_tmps_["349_abab_vvoo"]("b,a,i,j")  = t2["abab_vvoo"]("b,a,k,j") * reuse_tmps_["121_aa_oo"]("i,k");
            reuse_tmps_["122_bb_oo"]("i,j")  = dp["bb_ov"]("j,b") * u1["bb_vo"]("b,i");
            reuse_tmps_["319_bb_vo"]("a,i")  = u1["bb_vo"]("a,j") * reuse_tmps_["122_bb_oo"]("i,j");
            reuse_tmps_["339_bbbb_vvoo"]("a,b,j,i")  = reuse_tmps_["122_bb_oo"]("i,k") * t2["bbbb_vvoo"]("b,a,j,k");
            reuse_tmps_["348_abab_vvoo"]("b,a,j,i")  = t2["abab_vvoo"]("b,a,j,k") * reuse_tmps_["122_bb_oo"]("i,k");
            reuse_tmps_["123_bbbb_vooo"]("e,i,n,m")  = eri["bbbb_oovv"]("m,n,a,e") * u1["bb_vo"]("a,i");
            reuse_tmps_["124_aaaa_vooo"]("e,i,n,m")  = eri["aaaa_oovv"]("m,n,a,e") * u1["aa_vo"]("a,i");
            reuse_tmps_["125_aabb_vooo"]("b,j,m,i")  = eri["abab_oovv"]("j,i,b,a") * u1["bb_vo"]("a,m");
            reuse_tmps_["126_baab_vooo"]("b,m,i,j")  = eri["abab_oovv"]("i,j,a,b") * u1["aa_vo"]("a,m");
            reuse_tmps_["127_bb_vo"]("b,i")  = eri["bbbb_oovv"]("j,i,a,b") * u1["bb_vo"]("a,j");
            reuse_tmps_["128_aa_vo"]("b,i")  = eri["aaaa_oovv"]("j,i,a,b") * u1["aa_vo"]("a,j");
        }
        reuse_tmps_["129_bbbb_vooo"]("b,j,i,k")  = dp["bb_ov"]("k,c") * t2["bbbb_vvoo"]("c,b,i,j");

        if (includes_["u1"]) {
            reuse_tmps_["320_bbbb_vvoo"]("a,b,i,j")  = reuse_tmps_["129_bbbb_vooo"]("b,j,i,k") * u1["bb_vo"]("a,k");
        }
        reuse_tmps_["130_aaaa_vooo"]("a,j,i,m")  = dp["aa_ov"]("m,b") * t2["aaaa_vvoo"]("b,a,i,j");

        if (includes_["u1"]) {
            reuse_tmps_["321_aaaa_vvoo"]("b,a,i,j")  = u1["aa_vo"]("a,k") * reuse_tmps_["130_aaaa_vooo"]("b,j,i,k");
            reuse_tmps_["131_aa_vo"]("b,j")  = eri["abab_oovv"]("j,i,b,a") * u1["bb_vo"]("a,i");
            reuse_tmps_["293_aaaa_vooo"]("a,m,j,i")  = t2["aaaa_vvoo"]("b,a,i,j") * reuse_tmps_["131_aa_vo"]("b,m");
            reuse_tmps_["294_baab_vooo"]("f,m,j,n")  = reuse_tmps_["131_aa_vo"]("a,j") * t2["abab_vvoo"]("a,f,m,n");
            reuse_tmps_["132_bb_vo"]("e,m")  = eri["abab_oovv"]("i,m,a,e") * u1["aa_vo"]("a,i");
            reuse_tmps_["292_bbbb_vooo"]("a,m,j,i")  = t2["bbbb_vvoo"]("b,a,i,j") * reuse_tmps_["132_bb_vo"]("b,m");
            reuse_tmps_["295_aabb_vooo"]("e,m,n,j")  = reuse_tmps_["132_bb_vo"]("a,j") * t2["abab_vvoo"]("e,a,m,n");
        }
        reuse_tmps_["133_aabb_vooo"]("e,m,n,j")  = eri["abba_oovo"]("i,j,a,m") * t2["abab_vvoo"]("e,a,i,n");
        reuse_tmps_["134_baab_vooo"]("f,m,j,n")  = eri["abab_oovo"]("j,i,a,n") * t2["abab_vvoo"]("a,f,m,i");

        if (includes_["u2"]) {
            reuse_tmps_["135_bbbb_vooo"]("a,j,i,m")  = dp["bb_ov"]("m,b") * u2["bbbb_vvoo"]("b,a,i,j");
            reuse_tmps_["136_aaaa_vooo"]("a,j,i,m")  = dp["aa_ov"]("m,b") * u2["aaaa_vvoo"]("b,a,i,j");
        }

        if (includes_["u1"]) {
            reuse_tmps_["137_bb_vo"]("e,m")  = eri["bbbb_oovv"]("m,i,a,e") * u1["bb_vo"]("a,i");
            reuse_tmps_["138_aa_vo"]("e,m")  = eri["aaaa_oovv"]("m,i,a,e") * u1["aa_vo"]("a,i");
            reuse_tmps_["139_bb_vo"]("b,m")  = eri["bbbb_oovv"]("m,k,b,c") * u1["bb_vo"]("c,k");
            reuse_tmps_["286_bbbb_vooo"]("a,m,j,i")  = t2["bbbb_vvoo"]("b,a,i,j") * reuse_tmps_["139_bb_vo"]("b,m");
            reuse_tmps_["296_aabb_vooo"]("e,m,n,j")  = reuse_tmps_["139_bb_vo"]("a,j") * t2["abab_vvoo"]("e,a,m,n");
            reuse_tmps_["140_aa_vo"]("b,m")  = eri["aaaa_oovv"]("m,k,b,c") * u1["aa_vo"]("c,k");
            reuse_tmps_["287_aaaa_vooo"]("a,m,j,i")  = t2["aaaa_vvoo"]("b,a,i,j") * reuse_tmps_["140_aa_vo"]("b,m");
            reuse_tmps_["297_baab_vooo"]("f,m,j,n")  = reuse_tmps_["140_aa_vo"]("a,j") * t2["abab_vvoo"]("a,f,m,n");
            reuse_tmps_["141_bb_oo"]("i,m")  = f["bb_ov"]("m,a") * u1["bb_vo"]("a,i");
            reuse_tmps_["142_aa_oo"]("i,m")  = f["aa_ov"]("m,a") * u1["aa_vo"]("a,i");
        }
        reuse_tmps_["143_bbbb_vooo"]("a,i,j,m")  = eri["bbbb_oovo"]("m,k,b,j") * t2["bbbb_vvoo"]("b,a,i,k");
        reuse_tmps_["144_aaaa_vooo"]("a,i,j,m")  = eri["aaaa_oovo"]("m,k,b,j") * t2["aaaa_vvoo"]("b,a,i,k");
        reuse_tmps_["145_baab_vooo"]("f,m,i,n")  = dp["aa_ov"]("i,a") * t2["abab_vvoo"]("a,f,m,n");

        if (includes_["u1"]) {
            reuse_tmps_["337_abab_vvoo"]("e,f,m,n")  = reuse_tmps_["145_baab_vooo"]("f,m,i,n") * u1["aa_vo"]("e,i");
        }
        reuse_tmps_["146_aabb_vooo"]("a,i,j,m")  = dp["bb_ov"]("m,b") * t2["abab_vvoo"]("a,b,i,j");

        if (includes_["u1"]) {
            reuse_tmps_["338_abab_vvoo"]("b,a,i,j")  = u1["bb_vo"]("a,k") * reuse_tmps_["146_aabb_vooo"]("b,i,j,k");
            reuse_tmps_["147_bbbb_vooo"]("a,i,j,m")  = eri["bbbb_vovo"]("a,m,b,j") * u1["bb_vo"]("b,i");
            reuse_tmps_["148_aaaa_vooo"]("a,i,j,m")  = eri["aaaa_vovo"]("a,m,b,j") * u1["aa_vo"]("b,i");
        }
        reuse_tmps_["149_aaaa_vooo"]("e,n,m,j")  = eri["abba_oovo"]("j,i,a,m") * t2["abab_vvoo"]("e,a,n,i");
        reuse_tmps_["150_bbbb_vooo"]("a,i,j,m")  = eri["abab_oovo"]("k,m,b,j") * t2["abab_vvoo"]("b,a,k,i");
        reuse_tmps_["151_bb_vo"]("a,i")  = dp["bb_ov"]("j,b") * t2["bbbb_vvoo"]("b,a,i,j");

        if (includes_["u1"]) {
            reuse_tmps_["329_bbbb_vvoo"]("b,a,j,i")  = u1["bb_vo"]("a,i") * reuse_tmps_["151_bb_vo"]("b,j");
            reuse_tmps_["352_abab_vvoo"]("a,b,i,j")  = u1["aa_vo"]("a,i") * reuse_tmps_["151_bb_vo"]("b,j");
        }
        reuse_tmps_["152_aa_vo"]("a,i")  = dp["aa_ov"]("j,b") * t2["aaaa_vvoo"]("b,a,i,j");

        if (includes_["u1"]) {
            reuse_tmps_["331_aaaa_vvoo"]("b,a,j,i")  = u1["aa_vo"]("a,i") * reuse_tmps_["152_aa_vo"]("b,j");
            reuse_tmps_["353_abab_vvoo"]("b,a,j,i")  = u1["bb_vo"]("a,i") * reuse_tmps_["152_aa_vo"]("b,j");
        }

        if (includes_["u2"]) {
            reuse_tmps_["153_aabb_vooo"]("a,i,j,m")  = dp["bb_ov"]("m,b") * u2["abab_vvoo"]("a,b,i,j");
            reuse_tmps_["154_baab_vooo"]("a,i,m,j")  = dp["aa_ov"]("m,b") * u2["abab_vvoo"]("b,a,i,j");
        }
        reuse_tmps_["155_bbbb_vooo"]("a,j,i,m")  = f["bb_ov"]("m,b") * t2["bbbb_vvoo"]("b,a,i,j");
        reuse_tmps_["156_aaaa_vooo"]("a,j,i,m")  = f["aa_ov"]("m,b") * t2["aaaa_vvoo"]("b,a,i,j");

        if (includes_["u1"]) {
            reuse_tmps_["157_aa_vo"]("a,i")  = eri["aaaa_vovo"]("a,j,b,i") * u1["aa_vo"]("b,j");
            reuse_tmps_["158_bb_vo"]("a,i")  = eri["bbbb_vovo"]("a,j,b,i") * u1["bb_vo"]("b,j");
        }

        if (includes_["u2"]) {
            reuse_tmps_["159_baab_vooo"]("f,m,j,n")  = eri["abab_oovo"]("j,i,a,n") * u2["abab_vvoo"]("a,f,m,i");
            reuse_tmps_["160_aaaa_vooo"]("f,n,m,j")  = eri["abba_oovo"]("j,i,a,m") * u2["abab_vvoo"]("f,a,n,i");
            reuse_tmps_["161_bb_vo"]("a,i")  = dp["bb_ov"]("j,b") * u2["bbbb_vvoo"]("b,a,i,j");
            reuse_tmps_["162_aa_vo"]("a,i")  = dp["aa_ov"]("j,b") * u2["aaaa_vvoo"]("b,a,i,j");
            reuse_tmps_["163_bbbb_vooo"]("a,i,j,m")  = eri["bbbb_oovo"]("m,k,b,j") * u2["bbbb_vvoo"]("b,a,i,k");
            reuse_tmps_["164_aaaa_vooo"]("a,i,j,m")  = eri["aaaa_oovo"]("m,k,b,j") * u2["aaaa_vvoo"]("b,a,i,k");
        }
        reuse_tmps_["165_aabb_vooo"]("e,m,n,j")  = eri["bbbb_oovo"]("j,i,a,n") * t2["abab_vvoo"]("e,a,m,i");
        reuse_tmps_["166_baab_vooo"]("a,j,m,i")  = eri["aaaa_oovo"]("m,k,b,j") * t2["abab_vvoo"]("b,a,k,i");
        reuse_tmps_["167_aa_vo"]("f,m")  = dp["bb_ov"]("i,a") * t2["abab_vvoo"]("f,a,m,i");

        if (includes_["u1"]) {
            reuse_tmps_["328_aaaa_vvoo"]("b,a,j,i")  = u1["aa_vo"]("a,i") * reuse_tmps_["167_aa_vo"]("b,j");
            reuse_tmps_["351_abab_vvoo"]("b,a,j,i")  = u1["bb_vo"]("a,i") * reuse_tmps_["167_aa_vo"]("b,j");
        }
        reuse_tmps_["168_bb_vo"]("e,m")  = dp["aa_ov"]("i,a") * t2["abab_vvoo"]("a,e,i,m");

        if (includes_["u1"]) {
            reuse_tmps_["330_bbbb_vvoo"]("b,a,j,i")  = u1["bb_vo"]("a,i") * reuse_tmps_["168_bb_vo"]("b,j");
            reuse_tmps_["354_abab_vvoo"]("a,b,i,j")  = u1["aa_vo"]("a,i") * reuse_tmps_["168_bb_vo"]("b,j");
        }
        reuse_tmps_["169_baab_vooo"]("f,m,i,n")  = f["aa_ov"]("i,a") * t2["abab_vvoo"]("a,f,m,n");
        reuse_tmps_["170_aabb_vooo"]("a,i,j,m")  = f["bb_ov"]("m,b") * t2["abab_vvoo"]("a,b,i,j");

        if (includes_["u2"]) {
            reuse_tmps_["171_bbbb_vooo"]("a,i,j,m")  = eri["abab_oovo"]("k,m,b,j") * u2["abab_vvoo"]("b,a,k,i");
            reuse_tmps_["172_aa_vo"]("f,m")  = dp["bb_ov"]("i,a") * u2["abab_vvoo"]("f,a,m,i");
            reuse_tmps_["173_bb_vo"]("a,i")  = dp["aa_ov"]("j,b") * u2["abab_vvoo"]("b,a,j,i");
        }
        reuse_tmps_["174_aa_vo"]("a,i")  = 0.50 * eri["aaaa_oovo"]("k,j,b,i") * t2["aaaa_vvoo"]("b,a,k,j");
        reuse_tmps_["175_bb_vo"]("a,i")  = 0.50 * eri["bbbb_oovo"]("k,j,b,i") * t2["bbbb_vvoo"]("b,a,k,j");
        reuse_tmps_["176_aa_vo"]("a,i")  = f["aa_ov"]("j,b") * t2["aaaa_vvoo"]("b,a,i,j");
        reuse_tmps_["177_bb_vo"]("a,i")  = f["bb_ov"]("j,b") * t2["bbbb_vvoo"]("b,a,i,j");

        if (includes_["u2"]) {
            reuse_tmps_["178_bbbb_vooo"]("a,j,i,m")  = f["bb_ov"]("m,b") * u2["bbbb_vvoo"]("b,a,i,j");
            reuse_tmps_["179_aaaa_vooo"]("a,j,i,m")  = f["aa_ov"]("m,b") * u2["aaaa_vvoo"]("b,a,i,j");
        }

        if (includes_["u1"]) {
            reuse_tmps_["180_baab_vooo"]("f,m,i,n")  = eri["baab_vovo"]("f,i,a,n") * u1["aa_vo"]("a,m");
            reuse_tmps_["181_aabb_vooo"]("a,j,i,m")  = eri["abba_vovo"]("a,m,b,j") * u1["bb_vo"]("b,i");
        }
        reuse_tmps_["182_aabb_vooo"]("e,m,n,j")  = eri["abab_oovo"]("i,j,a,n") * t2["aaaa_vvoo"]("a,e,m,i");
        reuse_tmps_["183_baab_vooo"]("a,j,m,i")  = eri["abba_oovo"]("m,k,b,j") * t2["bbbb_vvoo"]("b,a,i,k");

        if (includes_["u2"]) {
            reuse_tmps_["184_aabb_vooo"]("e,m,n,j")  = eri["abba_oovo"]("i,j,a,m") * u2["abab_vvoo"]("e,a,i,n");
            reuse_tmps_["185_aabb_vooo"]("e,m,n,j")  = eri["bbbb_oovo"]("j,i,a,n") * u2["abab_vvoo"]("e,a,m,i");
            reuse_tmps_["186_baab_vooo"]("a,j,m,i")  = eri["aaaa_oovo"]("m,k,b,j") * u2["abab_vvoo"]("b,a,k,i");
        }
        reuse_tmps_["187_aa_vo"]("f,n")  = eri["abba_oovo"]("j,i,a,n") * t2["abab_vvoo"]("f,a,j,i");
        reuse_tmps_["188_bb_vo"]("a,i")  = eri["abab_oovo"]("k,j,b,i") * t2["abab_vvoo"]("b,a,k,j");
        reuse_tmps_["189_aa_vo"]("f,m")  = f["bb_ov"]("i,a") * t2["abab_vvoo"]("f,a,m,i");
        reuse_tmps_["190_bb_vo"]("a,i")  = f["aa_ov"]("j,b") * t2["abab_vvoo"]("b,a,j,i");

        if (includes_["u2"]) {
            reuse_tmps_["191_baab_vooo"]("f,m,i,n")  = f["aa_ov"]("i,a") * u2["abab_vvoo"]("a,f,m,n");
            reuse_tmps_["192_aabb_vooo"]("a,i,j,m")  = f["bb_ov"]("m,b") * u2["abab_vvoo"]("a,b,i,j");
        }

        if (includes_["u1"]) {
            reuse_tmps_["193_aa_vo"]("a,i")  = eri["abba_vovo"]("a,j,b,i") * u1["bb_vo"]("b,j");
            reuse_tmps_["194_bb_vo"]("a,i")  = eri["baab_vovo"]("a,j,b,i") * u1["aa_vo"]("b,j");
        }

        if (includes_["u2"]) {
            reuse_tmps_["195_aa_vo"]("a,i")  = 0.50 * eri["aaaa_oovo"]("k,j,b,i") * u2["aaaa_vvoo"]("b,a,k,j");
            reuse_tmps_["196_bb_vo"]("a,i")  = 0.50 * eri["bbbb_oovo"]("k,j,b,i") * u2["bbbb_vvoo"]("b,a,k,j");
            reuse_tmps_["197_aa_vo"]("a,i")  = f["aa_ov"]("j,b") * u2["aaaa_vvoo"]("b,a,i,j");
            reuse_tmps_["198_bb_vo"]("a,i")  = f["bb_ov"]("j,b") * u2["bbbb_vvoo"]("b,a,i,j");
        }

        if (includes_["u1"]) {
            reuse_tmps_["199_bbbb_vooo"]("a,j,i,m")  = eri["bbbb_oooo"]("m,k,i,j") * u1["bb_vo"]("a,k");
            reuse_tmps_["200_aaaa_vooo"]("a,j,i,m")  = eri["aaaa_oooo"]("m,k,i,j") * u1["aa_vo"]("a,k");
        }

        if (includes_["u2"]) {
            reuse_tmps_["201_aabb_vooo"]("e,m,n,j")  = eri["abab_oovo"]("i,j,a,n") * u2["aaaa_vvoo"]("a,e,m,i");
            reuse_tmps_["202_baab_vooo"]("a,j,m,i")  = eri["abba_oovo"]("m,k,b,j") * u2["bbbb_vvoo"]("b,a,i,k");
            reuse_tmps_["203_aa_vo"]("f,n")  = eri["abba_oovo"]("j,i,a,n") * u2["abab_vvoo"]("f,a,j,i");
            reuse_tmps_["204_bb_vo"]("a,i")  = eri["abab_oovo"]("k,j,b,i") * u2["abab_vvoo"]("b,a,k,j");
            reuse_tmps_["205_bb_vo"]("f,n")  = f["aa_ov"]("i,a") * u2["abab_vvoo"]("a,f,i,n");
            reuse_tmps_["206_aa_vo"]("f,n")  = f["bb_ov"]("i,a") * u2["abab_vvoo"]("f,a,n,i");
        }

        if (includes_["u1"]) {
            reuse_tmps_["207_baab_vooo"]("f,m,i,n")  = eri["baba_vovo"]("f,i,a,m") * u1["bb_vo"]("a,n");
            reuse_tmps_["208_aabb_vooo"]("a,i,j,m")  = eri["abab_vovo"]("a,m,b,j") * u1["aa_vo"]("b,i");
            reuse_tmps_["209_baab_vooo"]("f,m,j,n")  = eri["abab_oooo"]("j,i,m,n") * u1["bb_vo"]("f,i");
            reuse_tmps_["210_aabb_vooo"]("a,i,j,m")  = eri["abab_oooo"]("k,m,i,j") * u1["aa_vo"]("a,k");
        }
        reuse_tmps_["211_bbbb_vvoo"]("e,f,n,m")  = eri["bbbb_vovo"]("f,i,a,m") * t2["bbbb_vvoo"]("a,e,n,i");
        reuse_tmps_["212_aaaa_vvoo"]("f,e,m,n")  = eri["abba_vovo"]("e,i,a,n") * t2["abab_vvoo"]("f,a,m,i");
        reuse_tmps_["213_bbbb_vvoo"]("a,b,i,j")  = eri["baab_vovo"]("b,k,c,j") * t2["abab_vvoo"]("c,a,k,i");
        reuse_tmps_["214_aaaa_vvoo"]("a,b,i,j")  = eri["aaaa_vovo"]("b,k,c,j") * t2["aaaa_vvoo"]("c,a,i,k");
        reuse_tmps_["215_abab_vvoo"]("e,f,m,n")  = eri["abba_vovo"]("e,i,a,m") * t2["bbbb_vvoo"]("a,f,n,i");
        reuse_tmps_["216_abab_vvoo"]("e,f,m,n")  = eri["baba_vovo"]("f,i,a,m") * t2["abab_vvoo"]("e,a,i,n");
        reuse_tmps_["217_abab_vvoo"]("e,f,m,n")  = eri["abab_vovo"]("e,i,a,n") * t2["abab_vvoo"]("a,f,m,i");
        reuse_tmps_["218_abab_vvoo"]("e,f,m,n")  = eri["baab_vovo"]("f,i,a,n") * t2["aaaa_vvoo"]("a,e,m,i");
        reuse_tmps_["219_abab_vvoo"]("a,b,i,j")  = eri["bbbb_vovo"]("b,k,c,j") * t2["abab_vvoo"]("a,c,i,k");
        reuse_tmps_["220_abab_vvoo"]("b,a,j,i")  = eri["aaaa_vovo"]("b,k,c,j") * t2["abab_vvoo"]("c,a,k,i");
        reuse_tmps_["221_bbbb_vvoo"]("a,b,j,i")  = eri["bbbb_oooo"]("l,k,i,j") * t2["bbbb_vvoo"]("b,a,l,k");
        reuse_tmps_["222_aaaa_vvoo"]("a,b,j,i")  = eri["aaaa_oooo"]("l,k,i,j") * t2["aaaa_vvoo"]("b,a,l,k");
        reuse_tmps_["223_abab_vvoo"]("b,a,i,j")  = eri["abab_oooo"]("l,k,i,j") * t2["abab_vvoo"]("b,a,l,k");

        if (includes_["u1"]) {
            reuse_tmps_["224_aa_vo"]("a,i")  = dp["aa_vv"]("a,b") * u1["aa_vo"]("b,i");
            reuse_tmps_["225_bb_vo"]("a,i")  = dp["bb_vv"]("a,b") * u1["bb_vo"]("b,i");
            reuse_tmps_["226_aa_vo"]("a,i")  = dp["aa_oo"]("j,i") * u1["aa_vo"]("a,j");
            reuse_tmps_["227_bb_vo"]("a,i")  = dp["bb_oo"]("j,i") * u1["bb_vo"]("a,j");
            reuse_tmps_["228_bb_vo"]("a,i")  = f["bb_vv"]("a,b") * u1["bb_vo"]("b,i");
            reuse_tmps_["229_aa_vo"]("a,i")  = f["aa_vv"]("a,b") * u1["aa_vo"]("b,i");
            reuse_tmps_["230_bb_vo"]("a,i")  = f["bb_oo"]("j,i") * u1["bb_vo"]("a,j");
            reuse_tmps_["231_aa_vo"]("a,i")  = f["aa_oo"]("j,i") * u1["aa_vo"]("a,j");
        }
        reuse_tmps_["232_bbbb_vvoo"]("f,e,n,m")  = dp["bb_vv"]("e,a") * t2["bbbb_vvoo"]("a,f,m,n");
        reuse_tmps_["233_aaaa_vvoo"]("e,f,n,m")  = dp["aa_vv"]("f,a") * t2["aaaa_vvoo"]("a,e,m,n");
        reuse_tmps_["234_aaaa_vvoo"]("a,b,i,j")  = dp["aa_oo"]("k,j") * t2["aaaa_vvoo"]("b,a,i,k");
        reuse_tmps_["235_bbbb_vvoo"]("a,b,i,j")  = dp["bb_oo"]("k,j") * t2["bbbb_vvoo"]("b,a,i,k");

        if (includes_["u1"]) {
            reuse_tmps_["236_aaaa_vvoo"]("a,b,i,j")  = dp["aa_vo"]("b,j") * u1["aa_vo"]("a,i");
            reuse_tmps_["237_bbbb_vvoo"]("a,b,i,j")  = dp["bb_vo"]("b,j") * u1["bb_vo"]("a,i");
        }
        reuse_tmps_["238_abab_vvoo"]("e,f,m,n")  = dp["aa_vv"]("e,a") * t2["abab_vvoo"]("a,f,m,n");
        reuse_tmps_["239_abab_vvoo"]("a,b,i,j")  = dp["bb_vv"]("b,c") * t2["abab_vvoo"]("a,c,i,j");
        reuse_tmps_["240_abab_vvoo"]("e,f,m,n")  = dp["aa_oo"]("i,m") * t2["abab_vvoo"]("e,f,i,n");
        reuse_tmps_["241_abab_vvoo"]("b,a,i,j")  = dp["bb_oo"]("k,j") * t2["abab_vvoo"]("b,a,i,k");
        reuse_tmps_["242_bbbb_vvoo"]("f,e,n,m")  = f["bb_vv"]("e,a") * t2["bbbb_vvoo"]("a,f,m,n");

        if (includes_["u2"]) {
            reuse_tmps_["243_bbbb_vvoo"]("e,f,n,m")  = dp["bb_vv"]("f,a") * u2["bbbb_vvoo"]("a,e,m,n");
            reuse_tmps_["244_aaaa_vvoo"]("e,f,n,m")  = dp["aa_vv"]("f,a") * u2["aaaa_vvoo"]("a,e,m,n");
        }
        reuse_tmps_["245_aaaa_vvoo"]("f,e,n,m")  = f["aa_vv"]("e,a") * t2["aaaa_vvoo"]("a,f,m,n");
        reuse_tmps_["246_aaaa_vvoo"]("a,b,i,j")  = f["aa_oo"]("k,j") * t2["aaaa_vvoo"]("b,a,i,k");

        if (includes_["u2"]) {
            reuse_tmps_["247_aaaa_vvoo"]("a,b,i,j")  = dp["aa_oo"]("k,j") * u2["aaaa_vvoo"]("b,a,i,k");
            reuse_tmps_["248_bbbb_vvoo"]("a,b,i,j")  = dp["bb_oo"]("k,j") * u2["bbbb_vvoo"]("b,a,i,k");
        }
        reuse_tmps_["249_bbbb_vvoo"]("a,b,i,j")  = f["bb_oo"]("k,j") * t2["bbbb_vvoo"]("b,a,i,k");

        if (includes_["u2"]) {
            reuse_tmps_["250_abab_vvoo"]("e,f,m,n")  = dp["aa_vv"]("e,a") * u2["abab_vvoo"]("a,f,m,n");
        }
        reuse_tmps_["251_abab_vvoo"]("e,f,m,n")  = f["aa_vv"]("e,a") * t2["abab_vvoo"]("a,f,m,n");
        reuse_tmps_["252_abab_vvoo"]("a,b,i,j")  = f["bb_vv"]("b,c") * t2["abab_vvoo"]("a,c,i,j");

        if (includes_["u2"]) {
            reuse_tmps_["253_abab_vvoo"]("a,b,i,j")  = dp["bb_vv"]("b,c") * u2["abab_vvoo"]("a,c,i,j");
            reuse_tmps_["254_abab_vvoo"]("e,f,m,n")  = dp["aa_oo"]("i,m") * u2["abab_vvoo"]("e,f,i,n");
        }
        reuse_tmps_["255_abab_vvoo"]("e,f,m,n")  = f["aa_oo"]("i,m") * t2["abab_vvoo"]("e,f,i,n");
        reuse_tmps_["256_abab_vvoo"]("b,a,i,j")  = f["bb_oo"]("k,j") * t2["abab_vvoo"]("b,a,i,k");

        if (includes_["u2"]) {
            reuse_tmps_["257_abab_vvoo"]("b,a,i,j")  = dp["bb_oo"]("k,j") * u2["abab_vvoo"]("b,a,i,k");
        }

        if (includes_["u1"]) {
            reuse_tmps_["258_abab_vvoo"]("a,b,i,j")  = dp["bb_vo"]("b,j") * u1["aa_vo"]("a,i");
            reuse_tmps_["259_abab_vvoo"]("b,a,j,i")  = dp["aa_vo"]("b,j") * u1["bb_vo"]("a,i");
        }

    }

    world_.gop.fence();
}