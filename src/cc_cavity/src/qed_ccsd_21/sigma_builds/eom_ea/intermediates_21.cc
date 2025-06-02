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

#include "cc_cavity/include/qed_ccsd_21/eom_ea_qed_ccsd_21.h"
#include "cc_cavity/include/qed_ccsd_21/qed_ccsd_21.h"

void hilbert::EOM_EA_QED_CCSD_21::build_hamiltonian() {
    // TODO: Implement build_hamiltonian
    throw PsiException("EOM_EA_CCSD::build_hamiltonian not implemented", __FILE__, __LINE__);
}

void hilbert::EOM_EA_QED_CCSD_21::build_common_ops() {

    if (!cc_wfn_->has_t1_integrals_) cc_wfn_->transform_integrals(true);

    // Get cavity information
    double w0 = cc_wfn_->cavity_frequency_[2];
    double coupling_factor_z = w0 * cc_wfn_->cavity_coupling_strength_[2];

    // get dipole integrals
    TArrayMap dp = reinterpret_pointer_cast<QED_CCSD_21>(cc_wfn_)->effective_dipole();

    // extract 0-body amplitudes
    double t0_1;
    foreach_inplace(cc_wfn_->amplitudes_["t0_1"], [&t0_1](auto &tile){
        for(auto &x : tile.range())
            t0_1 = tile[x];
    });

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

        scalars_["1"]  = 1.00 * dot(Id["bb_oo"]("i,i1"), dp["bb_oo"]("i,i1"));
        scalars_["2"]  = 1.00 * dot(Id["aa_oo"]("i,i1"), dp["aa_oo"]("i,i1"));
        scalars_["3"]  = 1.00 * dot(dp["bb_ov"]("i,b"), t1_1["bb"]("b,i"));
        scalars_["4"]  = 1.00 * dot(dp["aa_ov"]("i,b"), t1_1["aa"]("b,i"));
        scalars_["5"]  = -1.00 * dot(Id["aa_oo"]("o,I"), f["aa_oo"]("o,I"));
        scalars_["5"] -= 1.00 * dot(Id["bb_oo"]("i,j"), f["bb_oo"]("i,j"));
        scalars_["5"] += 0.25 * dot(eri["bbbb_oovv"]("i,m,d,b"), t2["bbbb"]("d,b,m,i"));
        scalars_["5"] += 1.00 * dot(Id["bb_oo"]("i,j") * eri["abab_oooo"]("k,i,l,j"), Id["aa_oo"]("k,l"));
        scalars_["5"] += 0.25 * dot(eri["aaaa_oovv"]("o,k,a,c"), t2["aaaa"]("a,c,k,o"));
        scalars_["5"] += scalars_["2"] * t0_1;
        scalars_["5"] += 0.50 * dot(Id["bb_oo"]("i,j") * eri["bbbb_oooo"]("i,m,j,n"), Id["bb_oo"]("m,n"));
        scalars_["5"] += 0.50 * dot(Id["aa_oo"]("o,I") * eri["aaaa_oooo"]("o,k,I,l"), Id["aa_oo"]("k,l"));
        scalars_["5"] -= 1.00 * dot(eri["abab_oovv"]("k,i,a,b"), t2["abab"]("a,b,k,i"));
        scalars_["5"] += scalars_["1"] * t0_1;
        scalars_["6"]  = -0.25 * dot(eri["aaaa_oovv"]("i,k,a,c"), t2_1["aaaa"]("a,c,k,i"));
        scalars_["6"] -= 0.25 * dot(eri["bbbb_oovv"]("j,l,b,d"), t2_1["bbbb"]("b,d,l,j"));
        scalars_["6"] -= 1.00 * dot(dp["bb_ov"]("j,b"), t1_1["bb"]("b,j")) * t0_1;
        scalars_["6"] += 1.00 * dot(eri["abab_oovv"]("k,j,a,d"), t2_1["abab"]("a,d,k,j"));
        scalars_["6"] += t0_1 * w0;
        scalars_["6"] -= 1.00 * dot(dp["aa_ov"]("i,a"), t1_1["aa"]("a,i")) * t0_1;
        scalars_["6"] += 1.00 * dot(f["bb_ov"]("j,b"), t1_1["bb"]("b,j"));
        scalars_["6"] += 1.00 * dot(f["aa_ov"]("i,a"), t1_1["aa"]("a,i"));
    }

    {
        reused_["1_aa_vv"]("a,b")  = eri["aaaa_oovv"]("j,i,c,b") * t2["aaaa"]("c,a,i,j");
        reused_["2_aa_vv"]("a,b")  = eri["abab_oovv"]("i,j,b,c") * t2["abab"]("a,c,i,j");
        reused_["3_aa_vv"]("a,b")  = t0_1 * dp["aa_vv"]("a,b");
        reused_["4_aa_ov"]("i,a")  = t0_1 * dp["aa_ov"]("i,a");
        reused_["5_bb_ov"]("i,a")  = t0_1 * dp["bb_ov"]("i,a");
        reused_["6_aa_oo"]("i,j")  = eri["abab_oovv"]("i,k,a,b") * t2["abab"]("a,b,j,k");
        reused_["7_aa_oo"]("i,j")  = eri["aaaa_oovv"]("i,k,a,b") * t2["aaaa"]("a,b,j,k");
        reused_["8_aa_oo"]("i,j")  = t0_1 * dp["aa_oo"]("i,j");
        reused_["9_aa_ov"]("i,a")  = eri["aaaa_oovv"]("j,i,b,a") * t1_1["aa"]("b,j");
        reused_["10_bb_ov"]("i,a")  = eri["bbbb_oovv"]("j,i,b,a") * t1_1["bb"]("b,j");
        reused_["11_aa_ov"]("i,a")  = eri["abab_oovv"]("i,j,a,b") * t1_1["bb"]("b,j");
        reused_["12_bb_ov"]("i,a")  = eri["abab_oovv"]("j,i,b,a") * t1_1["aa"]("b,j");
        reused_["13_aa_vv"]("a,b")  = -1.00 * eri["abab_vovv"]("a,j,b,c") * t1_1["bb"]("c,j");
        reused_["13_aa_vv"]("a,b") += eri["abab_oovv"]("i,j,b,c") * t2_1["abab"]("a,c,i,j");
        reused_["14_aa_vv"]("a,b")  = eri["aaaa_vovv"]("a,i,c,b") * t1_1["aa"]("c,i");
        reused_["14_aa_vv"]("a,b") -= 0.50 * eri["aaaa_oovv"]("i,j,c,b") * t2_1["aaaa"]("c,a,j,i");
        reused_["15_aaaa_vovv"]("a,i,b,c")  = eri["aaaa_oovo"]("j,k,a,i") * t2["aaaa"]("b,c,k,j");
        reused_["16_aa_oo"]("i,j")  = eri["aaaa_oovv"]("k,j,a,b") * t2["aaaa"]("a,b,i,k");
        reused_["17_aa_oo"]("i,j")  = -0.50 * eri["aaaa_oovv"]("i,k,a,b") * t2_1["aaaa"]("a,b,j,k");
        reused_["17_aa_oo"]("i,j") += eri["aaaa_oovo"]("i,k,a,j") * t1_1["aa"]("a,k");
        reused_["18_aa_oo"]("i,j")  = dp["aa_ov"]("i,a") * t1_1["aa"]("a,j");
        reused_["19_aa_oo"]("i,j")  = f["aa_ov"]("i,a") * t1_1["aa"]("a,j");
        reused_["20_aa_oo"]("i,j")  = -1.00 * eri["abab_oovv"]("i,k,b,c") * t2_1["abab"]("b,c,j,k");
        reused_["20_aa_oo"]("i,j") += eri["abba_oovo"]("i,k,a,j") * t1_1["bb"]("a,k");
        reused_["21_aaaa_vvvo"]("a,b,c,i")  = (eri["aaaa_vvvv"]("a,b,d,c") + -0.50 * t2["aaaa"]("a,b,j,k") * eri["aaaa_oovv"]("k,j,d,c")) * t1_1["aa"]("d,i");
        reused_["22_aa_oo"]("i,j")  = -2.00 * eri["aaaa_oovo"]("k,j,a,i") * t1_1["aa"]("a,k");
        reused_["22_aa_oo"]("i,j") += eri["aaaa_oovv"]("k,j,a,b") * t2_1["aaaa"]("a,b,i,k");
        reused_["23_aaaa_vovv"]("a,i,b,c")  = eri["aaaa_oovo"]("j,k,a,i") * t2_1["aaaa"]("b,c,k,j");
        reused_["24_aa_ov"]("i,a")  = eri["aaaa_oovv"]("i,j,b,a") * t1_1["aa"]("b,j");
        reused_["25_aaaa_vovv"]("a,i,b,c")  = eri["aaaa_vovv"]("b,j,c,d") * t2["aaaa"]("d,a,i,j");
        reused_["26_aaaa_vvvo"]("a,b,c,i")  = eri["abab_vovv"]("a,j,b,d") * t2["abab"]("c,d,i,j");
        reused_["27_aabb_vvvo"]("a,b,c,i")  = eri["aaaa_vovv"]("a,j,b,d") * t2["abab"]("d,c,j,i");
        reused_["28_aabb_vvvo"]("a,b,c,i")  = eri["abab_vovv"]("a,j,b,d") * t2["bbbb"]("d,c,i,j");
        reused_["29_baab_vvvo"]("a,b,c,i")  = t2["abab"]("c,d,j,i") * eri["baab_vovv"]("a,j,b,d");
        reused_["30_abab_vovv"]("a,i,b,c")  = eri["abab_oovo"]("j,k,a,i") * t2["abab"]("b,c,j,k");
        reused_["31_aa_vv"]("a,b")  = eri["aaaa_oovv"]("i,j,a,c") * t2["aaaa"]("c,b,j,i");
        reused_["32_aa_vo"]("a,i")  = -1.00 * dp["aa_ov"]("k,c") * t2["aaaa"]("c,a,i,k");
        reused_["32_aa_vo"]("a,i") += dp["bb_ov"]("j,b") * t2["abab"]("a,b,i,j");
        reused_["33_bb_vo"]("a,i")  = -1.00 * dp["aa_ov"]("k,c") * t2["abab"]("c,a,k,i");
        reused_["33_bb_vo"]("a,i") += dp["bb_ov"]("j,b") * t2["bbbb"]("b,a,i,j");
        reused_["34_aa_vo"]("a,i")  = -1.00 * dp["aa_ov"]("j,b") * t2_1["aaaa"]("b,a,i,j");
        reused_["34_aa_vo"]("a,i") += dp["aa_vv"]("a,b") * t1_1["aa"]("b,i");
        reused_["34_aa_vo"]("a,i") += dp["bb_ov"]("k,c") * t2_1["abab"]("a,c,i,k");
        reused_["34_aa_vo"]("a,i") -= t1_1["aa"]("a,j") * dp["aa_oo"]("j,i");
        reused_["35_bb_ov"]("i,a")  = -1.00 * dp["aa_ov"]("k,c") * t2_1["abab"]("c,a,k,i");
        reused_["35_bb_ov"]("i,a") += dp["bb_oo"]("j,i") * t1_1["bb"]("a,j");
        reused_["35_bb_ov"]("i,a") += dp["bb_ov"]("j,b") * t2_1["bbbb"]("b,a,i,j");
        reused_["35_bb_ov"]("i,a") -= dp["bb_vv"]("a,b") * t1_1["bb"]("b,i");
        reused_["36_aa_vo"]("a,i")  = -0.50 * eri["aaaa_vovv"]("a,j,b,c") * t2["aaaa"]("b,c,i,j");
        reused_["36_aa_vo"]("a,i") -= 0.50 * t2["aaaa"]("b,a,n,j") * eri["aaaa_oovo"]("j,n,b,i");
        reused_["36_aa_vo"]("a,i") += t0_1 * dp["aa_vo"]("a,i");
        reused_["36_aa_vo"]("a,i") -= t2["abab"]("a,d,n,k") * eri["abba_oovo"]("n,k,d,i");
        reused_["36_aa_vo"]("a,i") += f["aa_ov"]("j,b") * t2["aaaa"]("b,a,i,j");
        reused_["36_aa_vo"]("a,i") -= f["bb_ov"]("k,d") * t2["abab"]("a,d,i,k");
        reused_["36_aa_vo"]("a,i") -= eri["abab_vovv"]("a,k,b,e") * t2["abab"]("b,e,i,k");
        reused_["36_aa_vo"]("a,i") += scalars_["1"] * t1_1["aa"]("a,i");
        reused_["36_aa_vo"]("a,i") += scalars_["2"] * t1_1["aa"]("a,i");
        reused_["36_aa_vo"]("a,i") += t0_1 * reused_["32_aa_vo"]("a,i");
        reused_["37_bb_vo"]("a,i")  = -1.00 * t0_1 * dp["bb_vo"]("a,i");
        reused_["37_bb_vo"]("a,i") -= scalars_["2"] * t1_1["bb"]("a,i");
        reused_["37_bb_vo"]("a,i") -= scalars_["1"] * t1_1["bb"]("a,i");
        reused_["37_bb_vo"]("a,i") -= eri["abab_oovo"]("m,j,c,i") * t2["abab"]("c,a,m,j");
        reused_["37_bb_vo"]("a,i") += 0.50 * eri["bbbb_vovv"]("a,j,b,d") * t2["bbbb"]("b,d,i,j");
        reused_["37_bb_vo"]("a,i") += f["aa_ov"]("k,c") * t2["abab"]("c,a,k,i");
        reused_["37_bb_vo"]("a,i") -= eri["baab_vovv"]("a,k,c,d") * t2["abab"]("c,d,k,i");
        reused_["37_bb_vo"]("a,i") += 0.50 * eri["bbbb_oovo"]("j,l,b,i") * t2["bbbb"]("b,a,l,j");
        reused_["37_bb_vo"]("a,i") -= f["bb_ov"]("j,b") * t2["bbbb"]("b,a,i,j");
        reused_["37_bb_vo"]("a,i") += t0_1 * reused_["33_bb_vo"]("a,i");
        reused_["38_aabb_vvvo"]("a,b,c,i")  = eri["aaaa_vovv"]("a,j,b,d") * t2_1["abab"]("d,c,j,i");
        reused_["39_aaaa_voov"]("a,i,j,b")  = eri["aaaa_oovv"]("k,j,b,c") * t2["aaaa"]("c,a,i,k");
        reused_["40_bbaa_voov"]("a,i,j,b")  = eri["aaaa_oovv"]("k,j,b,c") * t2["abab"]("c,a,k,i");
        reused_["41_aaaa_vvvo"]("a,b,c,i")  = 0.50 * eri["aaaa_oovv"]("j,k,c,d") * t1_1["aa"]("d,i") * t2["aaaa"]("b,a,k,j");
        reused_["41_aaaa_vvvo"]("a,b,c,i") += 2.00 * eri["aaaa_vovv"]("b,j,c,d") * t2_1["aaaa"]("d,a,i,j");
        reused_["41_aaaa_vvvo"]("a,b,c,i") += eri["aaaa_vvvv"]("a,b,c,d") * t1_1["aa"]("d,i");
        reused_["42_aaaa_vvvo"]("a,b,c,i")  = eri["abab_vovv"]("a,j,b,d") * t2_1["abab"]("c,d,i,j");
        reused_["43_abba_vovv"]("a,i,b,c")  = eri["baab_vovv"]("b,j,c,d") * t2_1["abab"]("a,d,j,i");
        reused_["44_aaaa_voov"]("a,i,j,b")  = eri["aaaa_oovv"]("k,j,c,b") * t2["aaaa"]("c,a,i,k");
        reused_["45_abba_voov"]("a,i,j,b")  = eri["abab_oovv"]("k,j,b,c") * t2["abab"]("a,c,k,i");
        reused_["46_abab_vovv"]("a,i,b,c")  = eri["abab_oovo"]("j,k,a,i") * t2_1["abab"]("b,c,j,k");
        reused_["47_aabb_voov"]("a,i,j,b")  = eri["abab_oovv"]("k,j,c,b") * t2["aaaa"]("c,a,i,k");
        reused_["48_bbaa_voov"]("a,i,j,b")  = eri["aaaa_oovv"]("k,j,c,b") * t2["abab"]("c,a,k,i");
        reused_["49_aaaa_voov"]("a,i,j,b")  = eri["abab_oovv"]("j,k,b,c") * t2["abab"]("a,c,i,k");
        reused_["50_bbaa_voov"]("a,i,j,b")  = eri["abab_oovv"]("j,k,b,c") * t2["bbbb"]("c,a,i,k");
        reused_["51_bbbb_voov"]("a,i,j,b")  = eri["abab_oovv"]("k,j,c,b") * t2["abab"]("c,a,k,i");
        reused_["52_aabb_voov"]("a,i,j,b")  = eri["bbbb_oovv"]("k,j,c,b") * t2["abab"]("a,c,i,k");
        reused_["53_bbbb_voov"]("a,i,j,b")  = eri["bbbb_oovv"]("k,j,c,b") * t2["bbbb"]("c,a,i,k");
        reused_["54_aa_ov"]("i,a")  = eri["aaaa_oovv"]("i,j,a,b") * t1_1["aa"]("b,j");
        reused_["55_bb_vv"]("a,b")  = eri["abab_oovv"]("i,j,c,b") * t2["abab"]("c,a,i,j");
        reused_["56_bb_vv"]("a,b")  = eri["bbbb_oovv"]("j,i,c,b") * t2["bbbb"]("c,a,i,j");
        reused_["57_bb_oo"]("i,j")  = eri["abab_oovv"]("k,i,a,b") * t2["abab"]("a,b,k,j");
        reused_["58_bb_oo"]("i,j")  = eri["bbbb_oovv"]("k,j,a,b") * t2["bbbb"]("a,b,i,k");
        reused_["59_bb_oo"]("i,j")  = dp["bb_ov"]("i,a") * t1_1["bb"]("a,j");
        reused_["60_abab_vovv"]("a,i,b,c")  = -1.00 * eri["abab_vovv"]("b,k,a,d") * t2_1["bbbb"]("d,c,i,k");
        reused_["60_abab_vovv"]("a,i,b,c") += eri["abab_oovv"]("j,k,a,d") * t1_1["bb"]("d,i") * t2["abab"]("b,c,j,k");
        reused_["60_abab_vovv"]("a,i,b,c") += eri["abab_vvvv"]("b,c,a,d") * t1_1["bb"]("d,i");
        reused_["61_aa_vv"]("a,b")  = -0.50 * t2_1["aaaa"]("c,a,j,i") * eri["aaaa_oovv"]("i,j,b,c");
        reused_["61_aa_vv"]("a,b") += eri["aaaa_vovv"]("a,i,b,c") * t1_1["aa"]("c,i");
        reused_["62_aa_vo"]("a,i")  = -0.50 * eri["aaaa_vovv"]("a,l,d,c") * t2_1["aaaa"]("d,c,i,l");
        reused_["62_aa_vo"]("a,i") -= 0.50 * t2_1["aaaa"]("d,a,k,l") * eri["aaaa_oovo"]("l,k,d,i");
        reused_["62_aa_vo"]("a,i") -= eri["abab_vovv"]("a,j,d,e") * t2_1["abab"]("d,e,i,j");
        reused_["62_aa_vo"]("a,i") -= f["aa_vv"]("a,d") * t1_1["aa"]("d,i");
        reused_["62_aa_vo"]("a,i") += f["aa_oo"]("l,i") * t1_1["aa"]("a,l");
        reused_["62_aa_vo"]("a,i") += f["aa_ov"]("l,d") * t2_1["aaaa"]("d,a,i,l");
        reused_["62_aa_vo"]("a,i") += eri["abba_vovo"]("a,j,b,i") * t1_1["bb"]("b,j");
        reused_["62_aa_vo"]("a,i") += eri["aaaa_vovo"]("a,l,d,i") * t1_1["aa"]("d,l");
        reused_["62_aa_vo"]("a,i") -= f["bb_ov"]("j,b") * t2_1["abab"]("a,b,i,j");
        reused_["62_aa_vo"]("a,i") -= eri["abba_oovo"]("k,j,b,i") * t2_1["abab"]("a,b,k,j");
        reused_["62_aa_vo"]("a,i") -= w0 * t1_1["aa"]("a,i");
        reused_["62_aa_vo"]("a,i") += scalars_["4"] * t1_1["aa"]("a,i");
        reused_["62_aa_vo"]("a,i") += scalars_["3"] * t1_1["aa"]("a,i");
        reused_["62_aa_vo"]("a,i") -= 0.50 * reused_["1_aa_vv"]("a,c") * t1_1["aa"]("c,i");
        reused_["62_aa_vo"]("a,i") -= 0.50 * t1_1["aa"]("a,k") * reused_["16_aa_oo"]("i,k");
        reused_["62_aa_vo"]("a,i") -= reused_["52_aabb_voov"]("a,i,m,e") * t1_1["bb"]("e,m");
        reused_["62_aa_vo"]("a,i") += reused_["6_aa_oo"]("k,i") * t1_1["aa"]("a,k");
        reused_["62_aa_vo"]("a,i") += reused_["47_aabb_voov"]("a,i,m,e") * t1_1["bb"]("e,m");
        reused_["62_aa_vo"]("a,i") += reused_["2_aa_vv"]("a,c") * t1_1["aa"]("c,i");
        reused_["62_aa_vo"]("a,i") += reused_["44_aaaa_voov"]("a,i,k,c") * t1_1["aa"]("c,k");
        reused_["62_aa_vo"]("a,i") -= 2.00 * reused_["18_aa_oo"]("l,i") * t1_1["aa"]("a,l");
        reused_["62_aa_vo"]("a,i") += t0_1 * reused_["34_aa_vo"]("a,i");
        reused_["63_bb_ov"]("i,a")  = -2.00 * reused_["50_bbaa_voov"]("a,i,k,e") * t1_1["aa"]("e,k");
        reused_["63_bb_ov"]("i,a") -= 2.00 * reused_["53_bbbb_voov"]("a,i,j,b") * t1_1["bb"]("b,j");
        reused_["63_bb_ov"]("i,a") += 2.00 * w0 * t1_1["bb"]("a,i");
        reused_["63_bb_ov"]("i,a") -= 2.00 * eri["abab_oovo"]("k,l,c,i") * t2_1["abab"]("c,a,k,l");
        reused_["63_bb_ov"]("i,a") -= 2.00 * scalars_["3"] * t1_1["bb"]("a,i");
        reused_["63_bb_ov"]("i,a") += 2.00 * f["aa_ov"]("m,c") * t2_1["abab"]("c,a,m,i");
        reused_["63_bb_ov"]("i,a") -= 2.00 * f["bb_oo"]("l,i") * t1_1["bb"]("a,l");
        reused_["63_bb_ov"]("i,a") -= 2.00 * eri["bbbb_vovo"]("a,l,d,i") * t1_1["bb"]("d,l");
        reused_["63_bb_ov"]("i,a") -= 2.00 * scalars_["4"] * t1_1["bb"]("a,i");
        reused_["63_bb_ov"]("i,a") -= 2.00 * f["bb_ov"]("l,d") * t2_1["bbbb"]("d,a,i,l");
        reused_["63_bb_ov"]("i,a") += eri["bbbb_vovv"]("a,l,d,b") * t2_1["bbbb"]("d,b,i,l");
        reused_["63_bb_ov"]("i,a") += 2.00 * f["bb_vv"]("a,d") * t1_1["bb"]("d,i");
        reused_["63_bb_ov"]("i,a") += eri["bbbb_oovo"]("l,j,d,i") * t2_1["bbbb"]("d,a,j,l");
        reused_["63_bb_ov"]("i,a") -= 2.00 * eri["baab_vovv"]("a,m,c,b") * t2_1["abab"]("c,b,m,i");
        reused_["63_bb_ov"]("i,a") -= 2.00 * eri["baab_vovo"]("a,m,c,i") * t1_1["aa"]("c,m");
        reused_["63_bb_ov"]("i,a") += reused_["58_bb_oo"]("i,j") * t1_1["bb"]("a,j");
        reused_["63_bb_ov"]("i,a") += reused_["56_bb_vv"]("a,b") * t1_1["bb"]("b,i");
        reused_["63_bb_ov"]("i,a") += 2.00 * reused_["48_bbaa_voov"]("a,i,k,e") * t1_1["aa"]("e,k");
        reused_["63_bb_ov"]("i,a") += 4.00 * reused_["59_bb_oo"]("l,i") * t1_1["bb"]("a,l");
        reused_["63_bb_ov"]("i,a") -= 2.00 * reused_["57_bb_oo"]("j,i") * t1_1["bb"]("a,j");
        reused_["63_bb_ov"]("i,a") += 2.00 * reused_["51_bbbb_voov"]("a,i,j,b") * t1_1["bb"]("b,j");
        reused_["63_bb_ov"]("i,a") -= 2.00 * reused_["55_bb_vv"]("a,b") * t1_1["bb"]("b,i");
        reused_["63_bb_ov"]("i,a") += 2.00 * t0_1 * reused_["35_bb_ov"]("i,a");
        reused_["64_aaaa_vvvo"]("a,b,c,i")  = eri["aaaa_vovv"]("a,j,d,b") * t2["aaaa"]("d,c,i,j");
        reused_["65_aabb_vvvo"]("a,b,c,i")  = eri["aaaa_vovv"]("a,j,d,b") * t2["abab"]("d,c,j,i");
        reused_["66_bb_vv"]("a,b")  = t0_1 * dp["bb_vv"]("a,b");
        reused_["67_bb_oo"]("i,j")  = t0_1 * dp["bb_oo"]("i,j");
        reused_["68_bbaa_voov"]("a,i,j,b")  = eri["aaaa_oovv"]("j,k,b,c") * t2_1["abab"]("c,a,k,i");
        reused_["69_aaaa_vovo"]("a,i,b,j")  = eri["aaaa_oovv"]("i,k,b,c") * t2_1["aaaa"]("c,a,j,k");
        reused_["69_aaaa_vovo"]("a,i,b,j") += eri["aaaa_vovv"]("a,i,b,c") * t1_1["aa"]("c,j");
        reused_["70_aaaa_voov"]("a,i,j,b")  = eri["aaaa_oovv"]("j,k,b,c") * t2["aaaa"]("c,a,i,k");
        reused_["71_bbaa_voov"]("a,i,j,b")  = eri["aaaa_oovv"]("j,k,b,c") * t2["abab"]("c,a,k,i");
        reused_["72_bbaa_voov"]("a,i,j,b")  = eri["abab_oovv"]("j,k,b,c") * t2_1["bbbb"]("c,a,i,k");
        reused_["73_baab_vovo"]("a,i,b,j")  = eri["baab_vovv"]("a,i,b,c") * t1_1["bb"]("c,j");
        reused_["74_aabb_voov"]("a,i,j,b")  = eri["bbbb_oovv"]("j,k,b,c") * t2_1["abab"]("a,c,i,k");
        reused_["75_bb_ov"]("i,a")  = eri["bbbb_oovv"]("i,j,a,b") * t1_1["bb"]("b,j");
        reused_["76_bbbb_voov"]("a,i,j,b")  = eri["bbbb_vovv"]("a,j,b,c") * t1_1["bb"]("c,i");
        reused_["76_bbbb_voov"]("a,i,j,b") += eri["bbbb_oovv"]("j,k,b,c") * t2_1["bbbb"]("c,a,i,k");
        reused_["77_bb_oo"]("i,j")  = -0.50 * eri["bbbb_oovv"]("i,k,a,b") * t2_1["bbbb"]("a,b,j,k");
        reused_["77_bb_oo"]("i,j") += eri["bbbb_oovo"]("i,k,a,j") * t1_1["bb"]("a,k");
        reused_["78_bb_vv"]("a,b")  = -0.50 * t2_1["bbbb"]("c,a,j,i") * eri["bbbb_oovv"]("i,j,b,c");
        reused_["78_bb_vv"]("a,b") += eri["bbbb_vovv"]("a,i,b,c") * t1_1["bb"]("c,i");
        reused_["79_aabb_voov"]("a,i,j,b")  = eri["abab_oovv"]("k,j,c,b") * t2_1["aaaa"]("c,a,i,k");
        reused_["80_bbbb_voov"]("a,i,j,b")  = eri["abab_oovv"]("k,j,c,b") * t2_1["abab"]("c,a,k,i");
        reused_["81_aabb_voov"]("a,i,j,b")  = eri["bbbb_oovv"]("j,k,b,c") * t2["abab"]("a,c,i,k");
        reused_["82_bbbb_voov"]("a,i,j,b")  = eri["bbbb_oovv"]("j,k,b,c") * t2["bbbb"]("c,a,i,k");
        reused_["83_abab_vovo"]("a,i,b,j")  = eri["abab_vovv"]("a,i,b,c") * t1_1["bb"]("c,j");
        reused_["84_abba_vovo"]("a,i,b,j")  = eri["abab_vovv"]("a,i,c,b") * t1_1["aa"]("c,j");
        reused_["85_bb_vv"]("a,b")  = eri["bbbb_oovv"]("i,j,a,c") * t2["bbbb"]("c,b,j,i");
        reused_["86_bb_oo"]("i,j")  = eri["bbbb_oovv"]("i,k,a,b") * t2["bbbb"]("a,b,j,k");
        reused_["87_bb_oo"]("i,j")  = f["bb_ov"]("j,a") * t1_1["bb"]("a,i");
        reused_["88_bb_vv"]("a,b")  = eri["baab_vovv"]("a,k,c,b") * t1_1["aa"]("c,k");
        reused_["88_bb_vv"]("a,b") += eri["abab_oovv"]("i,j,c,b") * t2_1["abab"]("c,a,i,j");
        reused_["89_bb_oo"]("i,j")  = eri["abab_oovo"]("k,i,a,j") * t1_1["aa"]("a,k");
        reused_["89_bb_oo"]("i,j") += eri["abab_oovv"]("k,i,a,b") * t2_1["abab"]("a,b,k,j");
        reused_["90_aabb_vvvo"]("a,b,c,i")  = eri["aaaa_vovv"]("a,j,d,b") * t2_1["abab"]("d,c,j,i");
        reused_["91_bbaa_voov"]("a,i,j,b")  = eri["aaaa_oovv"]("k,j,c,b") * t2_1["abab"]("c,a,k,i");
        reused_["92_bbbb_voov"]("a,i,j,b")  = eri["bbbb_oovv"]("k,j,c,b") * t2_1["bbbb"]("c,a,i,k");
        reused_["93_bbbb_vovo"]("a,i,b,j")  = eri["bbbb_vovv"]("a,i,c,b") * t1_1["bb"]("c,j");
        reused_["94_bb_vv"]("a,b")  = eri["bbbb_vovv"]("a,i,c,b") * t1_1["bb"]("c,i");
        reused_["94_bb_vv"]("a,b") -= 0.50 * eri["bbbb_oovv"]("i,j,c,b") * t2_1["bbbb"]("c,a,j,i");
        reused_["95_bb_oo"]("i,j")  = -0.50 * eri["bbbb_oovv"]("k,i,a,b") * t2_1["bbbb"]("a,b,j,k");
        reused_["95_bb_oo"]("i,j") += eri["bbbb_oovo"]("k,i,a,j") * t1_1["bb"]("a,k");
        reused_["96_aaaa_vvvo"]("a,b,c,i")  = eri["aaaa_vovv"]("a,j,d,b") * t2_1["aaaa"]("d,c,i,j");
        reused_["97_aaaa_voov"]("a,i,j,b")  = eri["aaaa_oovv"]("k,j,c,b") * t2_1["aaaa"]("c,a,i,k");
        reused_["98_aabb_voov"]("a,i,j,b")  = eri["bbbb_oovv"]("k,j,c,b") * t2_1["abab"]("a,c,i,k");
        reused_["99_aaaa_vovo"]("a,i,b,j")  = eri["aaaa_vovv"]("a,i,c,b") * t1_1["aa"]("c,j");
        reused_["100_aa_vo"]("a,i")  = reused_["12_bb_ov"]("j,b") * t2["abab"]("a,b,i,j");
    }

    world_.gop.fence();
}

