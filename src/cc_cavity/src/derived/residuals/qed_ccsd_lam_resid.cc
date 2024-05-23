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

#include "cc_cavity/include/lambda_driver.h"


using namespace std;
using namespace TA;
using namespace hilbert;

void LambdaDriver::build_intermediates() {

    if (!cc_wfn_->has_t1_integrals_) cc_wfn_->transform_integrals(true);

    // Get cavity information
    double w0 = cc_wfn_->cavity_frequency_[2];
    double coupling_factor_z = w0 * cc_wfn_->cavity_coupling_strength_[2];

    TA::TArrayD dp_aa_oo, dp_bb_oo, dp_aa_ov, dp_bb_ov, dp_aa_vo, dp_bb_vo, dp_aa_vv, dp_bb_vv;

    // get dipole integrals
    dp_aa_oo("i, j") = coupling_factor_z * cc_wfn_->Dip_blks_["dz_aa_oo"]("i, j");
    dp_aa_ov("i, a") = coupling_factor_z * cc_wfn_->Dip_blks_["dz_aa_ov"]("i, a");
    dp_aa_vo("a, i") = coupling_factor_z * cc_wfn_->Dip_blks_["dz_aa_vo"]("a, i");
    dp_aa_vv("a, b") = coupling_factor_z * cc_wfn_->Dip_blks_["dz_aa_vv"]("a, b");

    dp_bb_oo("i, j") = coupling_factor_z * cc_wfn_->Dip_blks_["dz_bb_oo"]("i, j");
    dp_bb_ov("i, a") = coupling_factor_z * cc_wfn_->Dip_blks_["dz_bb_ov"]("i, a");
    dp_bb_vo("a, i") = coupling_factor_z * cc_wfn_->Dip_blks_["dz_bb_vo"]("a, i");
    dp_bb_vv("a, b") = coupling_factor_z * cc_wfn_->Dip_blks_["dz_bb_vv"]("a, b");

    TA::TArrayD &t1_aa_vo =     cc_wfn_->amplitudes_["t1_aa"];
    TA::TArrayD &t1_bb_vo =     cc_wfn_->amplitudes_["t1_bb"];
    TA::TArrayD &u1_aa_vo =     cc_wfn_->amplitudes_["u1_aa"];
    TA::TArrayD &u1_bb_vo =     cc_wfn_->amplitudes_["u1_bb"];
    TA::TArrayD &t2_aaaa_vvoo = cc_wfn_->amplitudes_["t2_aaaa"];
    TA::TArrayD &t2_abab_vvoo = cc_wfn_->amplitudes_["t2_abab"];
    TA::TArrayD &t2_bbbb_vvoo = cc_wfn_->amplitudes_["t2_bbbb"];
    TA::TArrayD &u2_aaaa_vvoo = cc_wfn_->amplitudes_["u2_aaaa"];
    TA::TArrayD &u2_abab_vvoo = cc_wfn_->amplitudes_["u2_abab"];
    TA::TArrayD &u2_bbbb_vvoo = cc_wfn_->amplitudes_["u2_bbbb"];

    TArrayMap &V_blks_ = cc_wfn_->V_blks_;
    TArrayMap &F_blks_ = cc_wfn_->F_blks_;

    sharedOps = TArrayMap();

    {

        // sharedOps["abab_vvoo_0"] += 1.000000 V_blks_["abab_oovv"]("l,k,d,c") t2_aaaa_vvoo("d,a,i,l")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_0"]("a,c,i,k") = V_blks_["abab_oovv"]("l,k,d,c") * t2_aaaa_vvoo("d,a,i,l");

        // sharedOps["bbbb_vvoo_1"] += 1.000000 V_blks_["abab_oovv"]("l,k,d,c") t2_abab_vvoo("d,a,l,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_1"]("c,a,k,i") = V_blks_["abab_oovv"]("l,k,d,c") * t2_abab_vvoo("d,a,l,i");

        // sharedOps["bbbb_vvoo_2"] += 1.000000 t2_bbbb_vvoo("d,a,i,l") V_blks_["bbbb_oovv"]("l,k,c,d")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_2"]("a,c,i,k") = t2_bbbb_vvoo("d,a,i,l") * V_blks_["bbbb_oovv"]("l,k,c,d");

        // sharedOps["abab_vvoo_3"] += 1.000000 t2_abab_vvoo("a,d,i,l") V_blks_["bbbb_oovv"]("l,k,c,d")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_3"]("a,c,i,k") = t2_abab_vvoo("a,d,i,l") * V_blks_["bbbb_oovv"]("l,k,c,d");

        // sharedOps["abab_vvoo_4"] += 1.000000 V_blks_["aaaa_oovv"]("l,k,c,d") t2_abab_vvoo("d,a,l,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_4"]("c,a,k,i") = V_blks_["aaaa_oovv"]("l,k,c,d") * t2_abab_vvoo("d,a,l,i");

        // sharedOps["aaaa_vvoo_5"] += 1.000000 V_blks_["aaaa_oovv"]("l,k,c,d") t2_aaaa_vvoo("d,a,i,l")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_5"]("c,a,k,i") = V_blks_["aaaa_oovv"]("l,k,c,d") * t2_aaaa_vvoo("d,a,i,l");

        // sharedOps["aaaa_vvoo_6"] += 1.000000 t2_abab_vvoo("a,d,i,l") V_blks_["abab_oovv"]("k,l,c,d")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_6"]("a,c,i,k") = t2_abab_vvoo("a,d,i,l") * V_blks_["abab_oovv"]("k,l,c,d");

        // sharedOps["abab_vvoo_7"] += 1.000000 t2_bbbb_vvoo("d,a,i,l") V_blks_["abab_oovv"]("k,l,c,d")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_7"]("c,a,k,i") = t2_bbbb_vvoo("d,a,i,l") * V_blks_["abab_oovv"]("k,l,c,d");

        // sharedOps["bbaa_vvoo_8"] += 1.000000 V_blks_["abab_oovv"]("k,l,d,c") t2_abab_vvoo("d,a,i,l")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbaa_vvoo_8"]("c,a,k,i") = V_blks_["abab_oovv"]("k,l,d,c") * t2_abab_vvoo("d,a,i,l");

        // sharedOps["aabb_vvoo_9"] += 1.000000 V_blks_["abab_oovv"]("l,k,c,d") t2_abab_vvoo("a,d,l,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aabb_vvoo_9"]("c,a,k,i") = V_blks_["abab_oovv"]("l,k,c,d") * t2_abab_vvoo("a,d,l,i");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_10"] += 1.000000 V_blks_["abab_vvvv"]("b,a,c,d") u2_abab_vvoo("c,d,i,j")
        // flops: o2v4: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_10"]("b,a,i,j") = V_blks_["abab_vvvv"]("b,a,c,d") * u2_abab_vvoo("c,d,i,j");
    }

    {

        // sharedOps["abab_vvoo_11"] += 1.000000 V_blks_["abab_vvvv"]("b,a,c,d") t2_abab_vvoo("c,d,i,j")
        // flops: o2v4: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_11"]("b,a,i,j") = V_blks_["abab_vvvv"]("b,a,c,d") * t2_abab_vvoo("c,d,i,j");

        // sharedOps["aabb_vvvo_12"] += 1.000000 t2_bbbb_vvoo("c,a,i,j") V_blks_["abab_vovv"]("b,j,e,c")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["aabb_vvvo_12"]("b,e,a,i") = t2_bbbb_vvoo("c,a,i,j") * V_blks_["abab_vovv"]("b,j,e,c");

        // sharedOps["aabb_vvvo_13"] += 1.000000 V_blks_["aaaa_vovv"]("b,j,e,c") t2_abab_vvoo("c,a,j,i")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["aabb_vvvo_13"]("b,e,a,i") = V_blks_["aaaa_vovv"]("b,j,e,c") * t2_abab_vvoo("c,a,j,i");

        // sharedOps["bbbb_vvvo_14"] += 1.000000 V_blks_["baab_vovv"]("b,j,c,e") t2_abab_vvoo("c,a,j,i")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["bbbb_vvvo_14"]("b,e,a,i") = V_blks_["baab_vovv"]("b,j,c,e") * t2_abab_vvoo("c,a,j,i");

        // sharedOps["abba_vvvo_15"] += 1.000000 V_blks_["abab_vovv"]("b,j,c,e") t2_abab_vvoo("c,a,i,j")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["abba_vvvo_15"]("b,e,a,i") = V_blks_["abab_vovv"]("b,j,c,e") * t2_abab_vvoo("c,a,i,j");

        // sharedOps["abba_vvvo_16"] += 1.000000 V_blks_["baab_vovv"]("b,j,c,e") t2_aaaa_vvoo("c,a,i,j")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["abba_vvvo_16"]("a,b,e,i") = V_blks_["baab_vovv"]("b,j,c,e") * t2_aaaa_vvoo("c,a,i,j");

        // sharedOps["aaaa_vvvo_17"] += 1.000000 V_blks_["abab_vovv"]("b,j,e,c") t2_abab_vvoo("a,c,i,j")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["aaaa_vvvo_17"]("b,e,a,i") = V_blks_["abab_vovv"]("b,j,e,c") * t2_abab_vvoo("a,c,i,j");

        // sharedOps["aabb_vvvo_18"] += 1.000000 V_blks_["baab_vovv"]("b,j,e,c") t2_abab_vvoo("a,c,j,i")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["aabb_vvvo_18"]("e,a,b,i") = V_blks_["baab_vovv"]("b,j,e,c") * t2_abab_vvoo("a,c,j,i");

        // sharedOps["aaaa_vvvo_19"] += 1.000000 t2_aaaa_vvoo("c,a,i,j") V_blks_["aaaa_vovv"]("b,j,e,c")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["aaaa_vvvo_19"]("a,b,e,i") = t2_aaaa_vvoo("c,a,i,j") * V_blks_["aaaa_vovv"]("b,j,e,c");

        // sharedOps["bbbb_vvvo_20"] += 1.000000 t2_bbbb_vvoo("c,a,i,j") V_blks_["bbbb_vovv"]("b,j,e,c")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["bbbb_vvvo_20"]("a,b,e,i") = t2_bbbb_vvoo("c,a,i,j") * V_blks_["bbbb_vovv"]("b,j,e,c");

        // sharedOps["abba_vvvo_21"] += 1.000000 t2_abab_vvoo("a,c,i,j") V_blks_["bbbb_vovv"]("b,j,e,c")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["abba_vvvo_21"]("a,b,e,i") = t2_abab_vvoo("a,c,i,j") * V_blks_["bbbb_vovv"]("b,j,e,c");
    }

    if (include_u2_) {

        // sharedOps["aabb_vvvo_22"] += 1.000000 V_blks_["abab_vovv"]("b,j,e,c") u2_bbbb_vvoo("c,a,i,j")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["aabb_vvvo_22"]("b,e,a,i") = V_blks_["abab_vovv"]("b,j,e,c") * u2_bbbb_vvoo("c,a,i,j");

        // sharedOps["aaaa_vvvo_23"] += 1.000000 u2_abab_vvoo("a,c,i,j") V_blks_["abab_vovv"]("b,j,e,c")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["aaaa_vvvo_23"]("a,b,e,i") = u2_abab_vvoo("a,c,i,j") * V_blks_["abab_vovv"]("b,j,e,c");

        // sharedOps["aaaa_vvoo_24"] += 1.000000 u2_aaaa_vvoo("c,d,i,j") V_blks_["aaaa_vvvv"]("b,a,c,d")
        // flops: o2v4: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_24"]("b,a,i,j") = u2_aaaa_vvoo("c,d,i,j") * V_blks_["aaaa_vvvv"]("b,a,c,d");

        // sharedOps["abba_vvvo_25"] += 1.000000 u2_abab_vvoo("a,c,i,j") V_blks_["bbbb_vovv"]("b,j,e,c")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["abba_vvvo_25"]("a,b,e,i") = u2_abab_vvoo("a,c,i,j") * V_blks_["bbbb_vovv"]("b,j,e,c");
    }

    {

        // sharedOps["aaaa_vvoo_26"] += 1.000000 t2_aaaa_vvoo("c,d,i,j") V_blks_["aaaa_vvvv"]("b,a,c,d")
        // flops: o2v4: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_26"]("b,a,i,j") = t2_aaaa_vvoo("c,d,i,j") * V_blks_["aaaa_vvvv"]("b,a,c,d");
    }

    if (include_u2_) {

        // sharedOps["bbbb_vvoo_27"] += 1.000000 V_blks_["bbbb_vvvv"]("b,a,c,d") u2_bbbb_vvoo("c,d,i,j")
        // flops: o2v4: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_27"]("b,a,i,j") = V_blks_["bbbb_vvvv"]("b,a,c,d") * u2_bbbb_vvoo("c,d,i,j");
    }

    {

        // sharedOps["bbbb_vvoo_28"] += 1.000000 V_blks_["bbbb_vvvv"]("b,a,c,d") t2_bbbb_vvoo("c,d,i,j")
        // flops: o2v4: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_28"]("b,a,i,j") = V_blks_["bbbb_vvvv"]("b,a,c,d") * t2_bbbb_vvoo("c,d,i,j");
    }

    if (include_u2_) {

        // sharedOps["aaaa_vvvo_29"] += 1.000000 V_blks_["aaaa_vovv"]("b,j,e,c") u2_aaaa_vvoo("c,a,i,j")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["aaaa_vvvo_29"]("b,e,a,i") = V_blks_["aaaa_vovv"]("b,j,e,c") * u2_aaaa_vvoo("c,a,i,j");

        // sharedOps["abba_vvvo_30"] += 1.000000 u2_aaaa_vvoo("c,a,i,j") V_blks_["baab_vovv"]("b,j,c,e")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["abba_vvvo_30"]("a,b,e,i") = u2_aaaa_vvoo("c,a,i,j") * V_blks_["baab_vovv"]("b,j,c,e");

        // sharedOps["aabb_vvvo_31"] += 1.000000 V_blks_["aaaa_vovv"]("b,j,e,c") u2_abab_vvoo("c,a,j,i")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["aabb_vvvo_31"]("b,e,a,i") = V_blks_["aaaa_vovv"]("b,j,e,c") * u2_abab_vvoo("c,a,j,i");

        // sharedOps["aabb_vvvo_32"] += 1.000000 V_blks_["baab_vovv"]("b,j,e,c") u2_abab_vvoo("a,c,j,i")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["aabb_vvvo_32"]("e,a,b,i") = V_blks_["baab_vovv"]("b,j,e,c") * u2_abab_vvoo("a,c,j,i");

        // sharedOps["bbbb_vvvo_33"] += 1.000000 V_blks_["bbbb_vovv"]("b,j,e,c") u2_bbbb_vvoo("c,a,i,j")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["bbbb_vvvo_33"]("b,e,a,i") = V_blks_["bbbb_vovv"]("b,j,e,c") * u2_bbbb_vvoo("c,a,i,j");

        // sharedOps["abba_vvvo_34"] += 1.000000 u2_abab_vvoo("c,a,i,j") V_blks_["abab_vovv"]("b,j,c,e")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["abba_vvvo_34"]("b,a,e,i") = u2_abab_vvoo("c,a,i,j") * V_blks_["abab_vovv"]("b,j,c,e");

        // sharedOps["bbbb_vvvo_35"] += 1.000000 u2_abab_vvoo("c,a,j,i") V_blks_["baab_vovv"]("b,j,c,e")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["bbbb_vvvo_35"]("a,b,e,i") = u2_abab_vvoo("c,a,j,i") * V_blks_["baab_vovv"]("b,j,c,e");
    }

    if (include_u1_) {

        // sharedOps["bb_vo_36"] += 1.000000 V_blks_["abab_oovv"]("j,k,b,c") u1_aa_vo("b,j")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_36"]("c,k") = V_blks_["abab_oovv"]("j,k,b,c") * u1_aa_vo("b,j");

        // sharedOps["aa_vo_37"] += 1.000000 u1_bb_vo("b,j") V_blks_["abab_oovv"]("k,j,c,b")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_37"]("c,k") = u1_bb_vo("b,j") * V_blks_["abab_oovv"]("k,j,c,b");

        // sharedOps["aa_vo_38"] += 1.000000 V_blks_["aaaa_oovv"]("k,j,e,c") u1_aa_vo("c,j")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_38"]("e,k") = V_blks_["aaaa_oovv"]("k,j,e,c") * u1_aa_vo("c,j");

        // sharedOps["bb_vo_39"] += 1.000000 u1_bb_vo("c,j") V_blks_["bbbb_oovv"]("k,j,e,c")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_39"]("e,k") = u1_bb_vo("c,j") * V_blks_["bbbb_oovv"]("k,j,e,c");
    }

    {

        // sharedOps["baab_vooo_40"] += 1.000000 V_blks_["baab_vovv"]("b,k,c,d") t2_abab_vvoo("c,d,j,i")
        // flops: o3v3: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_40"]("b,k,j,i") = V_blks_["baab_vovv"]("b,k,c,d") * t2_abab_vvoo("c,d,j,i");

        // sharedOps["aabb_vooo_41"] += 1.000000 V_blks_["abab_vovv"]("b,k,c,d") t2_abab_vvoo("c,d,i,j")
        // flops: o3v3: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_41"]("b,i,k,j") = V_blks_["abab_vovv"]("b,k,c,d") * t2_abab_vvoo("c,d,i,j");

        // sharedOps["bbbb_vooo_42"] += 1.000000 t2_bbbb_vvoo("c,d,i,j") V_blks_["bbbb_vovv"]("b,k,c,d")
        // flops: o3v3: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["bbbb_vooo_42"]("b,i,j,k") = t2_bbbb_vvoo("c,d,i,j") * V_blks_["bbbb_vovv"]("b,k,c,d");

        // sharedOps["aaaa_vooo_43"] += 1.000000 V_blks_["aaaa_vovv"]("b,k,c,d") t2_aaaa_vvoo("c,d,i,j")
        // flops: o3v3: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aaaa_vooo_43"]("b,k,i,j") = V_blks_["aaaa_vovv"]("b,k,c,d") * t2_aaaa_vvoo("c,d,i,j");

        // sharedOps["bbbb_vvvo_44"] += 1.000000 t2_bbbb_vvoo("b,a,k,j") V_blks_["bbbb_oovo"]("k,j,e,i")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["bbbb_vvvo_44"]("b,a,e,i") = t2_bbbb_vvoo("b,a,k,j") * V_blks_["bbbb_oovo"]("k,j,e,i");
    }

    if (include_u2_) {

        // sharedOps["bbbb_vvoo_45"] += 1.000000 V_blks_["bbbb_oovv"]("j,n,f,b") u2_bbbb_vvoo("b,a,i,j")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_45"]("f,a,n,i") = V_blks_["bbbb_oovv"]("j,n,f,b") * u2_bbbb_vvoo("b,a,i,j");
    }

    {

        // sharedOps["abba_vvvo_46"] += 1.000000 V_blks_["abba_oovo"]("j,k,e,i") t2_abab_vvoo("b,a,j,k")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["abba_vvvo_46"]("b,e,a,i") = V_blks_["abba_oovo"]("j,k,e,i") * t2_abab_vvoo("b,a,j,k");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_47"] += 1.000000 u2_abab_vvoo("b,a,j,i") V_blks_["aaaa_oovv"]("j,n,e,b")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_47"]("e,a,n,i") = u2_abab_vvoo("b,a,j,i") * V_blks_["aaaa_oovv"]("j,n,e,b");

        // sharedOps["abab_vvoo_48"] += 1.000000 u2_abab_vvoo("a,b,i,j") V_blks_["bbbb_oovv"]("j,n,f,b")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_48"]("a,f,i,n") = u2_abab_vvoo("a,b,i,j") * V_blks_["bbbb_oovv"]("j,n,f,b");

        // sharedOps["bbbb_vvoo_49"] += 1.000000 u2_abab_vvoo("c,b,k,j") V_blks_["abab_oovv"]("k,l,c,d")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_49"]("b,d,j,l") = u2_abab_vvoo("c,b,k,j") * V_blks_["abab_oovv"]("k,l,c,d");

        // sharedOps["abab_vvoo_50"] += 1.000000 u2_aaaa_vvoo("c,b,j,k") V_blks_["abab_oovv"]("k,l,c,d")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_50"]("b,d,j,l") = u2_aaaa_vvoo("c,b,j,k") * V_blks_["abab_oovv"]("k,l,c,d");

        // sharedOps["aaaa_vvoo_51"] += 1.000000 u2_aaaa_vvoo("b,a,i,j") V_blks_["aaaa_oovv"]("j,n,e,b")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_51"]("a,e,i,n") = u2_aaaa_vvoo("b,a,i,j") * V_blks_["aaaa_oovv"]("j,n,e,b");
    }

    {

        // sharedOps["aabb_vvvo_53"] += 1.000000 V_blks_["abab_oovo"]("k,j,e,i") t2_abab_vvoo("b,a,k,j")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["aabb_vvvo_53"]("e,b,a,i") = V_blks_["abab_oovo"]("k,j,e,i") * t2_abab_vvoo("b,a,k,j");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_54"] += 1.000000 u2_bbbb_vvoo("c,b,j,k") V_blks_["abab_oovv"]("l,k,d,c")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_54"]("d,b,l,j") = u2_bbbb_vvoo("c,b,j,k") * V_blks_["abab_oovv"]("l,k,d,c");

        // sharedOps["aaaa_vvoo_55"] += 1.000000 u2_abab_vvoo("b,c,j,k") V_blks_["abab_oovv"]("l,k,d,c")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_55"]("b,d,j,l") = u2_abab_vvoo("b,c,j,k") * V_blks_["abab_oovv"]("l,k,d,c");
    }

    {

        // sharedOps["aaaa_vvvo_56"] += 1.000000 V_blks_["aaaa_oovo"]("k,j,e,i") t2_aaaa_vvoo("b,a,k,j")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["aaaa_vvvo_56"]("e,b,a,i") = V_blks_["aaaa_oovo"]("k,j,e,i") * t2_aaaa_vvoo("b,a,k,j");

        // sharedOps["aabb_vvvv_58"] += 1.000000 V_blks_["abab_oovv"]("k,j,e,c") t2_abab_vvoo("a,b,k,j")
        // flops: o2v4: 1, o0v4: 1 | mem: o0v4: 2,
        sharedOps["aabb_vvvv_58"]("e,a,c,b") = V_blks_["abab_oovv"]("k,j,e,c") * t2_abab_vvoo("a,b,k,j");
    }

    if (include_u1_) {

        // sharedOps["bbbb_vvoo_59"] += 1.000000 V_blks_["baab_vovv"]("b,k,d,c") u1_bb_vo("c,j") t2_abab_vvoo("d,a,k,i")
        // flops: o3v3: 1, o2v3: 1, o2v2: 1 | mem: o2v2: 3,
        sharedOps["bbbb_vvoo_59"]("b,a,j,i") = V_blks_["baab_vovv"]("b,k,d,c") * u1_bb_vo("c,j") * t2_abab_vvoo("d,a,k,i");

        // sharedOps["abab_vvoo_60"] += 1.000000 V_blks_["aaaa_vovv"]("b,k,c,d") u1_aa_vo("c,j") t2_abab_vvoo("d,a,k,i")
        // flops: o3v3: 1, o2v3: 1, o2v2: 1 | mem: o2v2: 3,
        sharedOps["abab_vvoo_60"]("b,a,j,i") = V_blks_["aaaa_vovv"]("b,k,c,d") * u1_aa_vo("c,j") * t2_abab_vvoo("d,a,k,i");

        // sharedOps["abab_vvoo_61"] += 1.000000 V_blks_["abab_vovv"]("b,k,d,c") u1_bb_vo("c,j") t2_abab_vvoo("d,a,i,k")
        // flops: o3v3: 1, o2v3: 1, o2v2: 1 | mem: o2v2: 3,
        sharedOps["abab_vvoo_61"]("b,a,i,j") = V_blks_["abab_vovv"]("b,k,d,c") * u1_bb_vo("c,j") * t2_abab_vvoo("d,a,i,k");

        // sharedOps["aaaa_vvoo_62"] += 1.000000 V_blks_["aaaa_vovv"]("b,k,c,d") u1_aa_vo("c,j") t2_aaaa_vvoo("d,a,i,k")
        // flops: o3v3: 1, o2v3: 1, o2v2: 1 | mem: o2v2: 3,
        sharedOps["aaaa_vvoo_62"]("b,a,j,i") = V_blks_["aaaa_vovv"]("b,k,c,d") * u1_aa_vo("c,j") * t2_aaaa_vvoo("d,a,i,k");

        // sharedOps["abab_vvoo_63"] += 1.000000 V_blks_["baab_vovv"]("b,k,c,d") u1_aa_vo("c,j") t2_abab_vvoo("a,d,k,i")
        // flops: o3v3: 1, o2v3: 1, o2v2: 1 | mem: o2v2: 3,
        sharedOps["abab_vvoo_63"]("a,b,j,i") = V_blks_["baab_vovv"]("b,k,c,d") * u1_aa_vo("c,j") * t2_abab_vvoo("a,d,k,i");

        // sharedOps["abab_vvoo_64"] += 1.000000 V_blks_["bbbb_vovv"]("b,k,c,d") u1_bb_vo("c,j") t2_abab_vvoo("a,d,i,k")
        // flops: o3v3: 1, o2v3: 1, o2v2: 1 | mem: o2v2: 3,
        sharedOps["abab_vvoo_64"]("a,b,i,j") = V_blks_["bbbb_vovv"]("b,k,c,d") * u1_bb_vo("c,j") * t2_abab_vvoo("a,d,i,k");

        // sharedOps["abab_vvoo_65"] += 1.000000 V_blks_["abab_vovv"]("b,k,c,d") u1_aa_vo("c,j") t2_bbbb_vvoo("d,a,i,k")
        // flops: o3v3: 1, o2v3: 1, o2v2: 1 | mem: o2v2: 3,
        sharedOps["abab_vvoo_65"]("b,a,j,i") = V_blks_["abab_vovv"]("b,k,c,d") * u1_aa_vo("c,j") * t2_bbbb_vvoo("d,a,i,k");

        // sharedOps["abab_vvoo_66"] += 1.000000 V_blks_["baab_vovv"]("b,k,d,c") u1_bb_vo("c,j") t2_aaaa_vvoo("d,a,i,k")
        // flops: o3v3: 1, o2v3: 1, o2v2: 1 | mem: o2v2: 3,
        sharedOps["abab_vvoo_66"]("a,b,i,j") = V_blks_["baab_vovv"]("b,k,d,c") * u1_bb_vo("c,j") * t2_aaaa_vvoo("d,a,i,k");

        // sharedOps["aaaa_vvoo_67"] += 1.000000 V_blks_["abab_vovv"]("b,k,c,d") u1_aa_vo("c,j") t2_abab_vvoo("a,d,i,k")
        // flops: o3v3: 1, o2v3: 1, o2v2: 1 | mem: o2v2: 3,
        sharedOps["aaaa_vvoo_67"]("b,a,j,i") = V_blks_["abab_vovv"]("b,k,c,d") * u1_aa_vo("c,j") * t2_abab_vvoo("a,d,i,k");

        // sharedOps["bbbb_vvoo_68"] += 1.000000 V_blks_["bbbb_vovv"]("b,k,c,d") u1_bb_vo("c,j") t2_bbbb_vvoo("d,a,i,k")
        // flops: o3v3: 1, o2v3: 1, o2v2: 1 | mem: o2v2: 3,
        sharedOps["bbbb_vvoo_68"]("b,a,j,i") = V_blks_["bbbb_vovv"]("b,k,c,d") * u1_bb_vo("c,j") * t2_bbbb_vvoo("d,a,i,k");

        // sharedOps["bbbb_vvvo_69"] += 1.000000 V_blks_["bbbb_oovv"]("k,j,e,c") u1_bb_vo("c,i") t2_bbbb_vvoo("b,a,k,j")
        // flops: o3v3: 1, o3v2: 1, o1v3: 1 | mem: o1v3: 2, o3v1: 1,
        sharedOps["bbbb_vvvo_69"]("e,b,a,i") = V_blks_["bbbb_oovv"]("k,j,e,c") * u1_bb_vo("c,i") * t2_bbbb_vvoo("b,a,k,j");

        // sharedOps["aaaa_vvvo_70"] += 1.000000 V_blks_["aaaa_oovv"]("k,j,e,c") u1_aa_vo("c,i") t2_aaaa_vvoo("b,a,k,j")
        // flops: o3v3: 1, o3v2: 1, o1v3: 1 | mem: o1v3: 2, o3v1: 1,
        sharedOps["aaaa_vvvo_70"]("e,b,a,i") = V_blks_["aaaa_oovv"]("k,j,e,c") * u1_aa_vo("c,i") * t2_aaaa_vvoo("b,a,k,j");
    }

    if (include_u2_) {

        // sharedOps["abba_vvvo_71"] += 1.000000 u2_abab_vvoo("a,b,k,j") V_blks_["abba_oovo"]("k,j,e,i")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["abba_vvvo_71"]("a,b,e,i") = u2_abab_vvoo("a,b,k,j") * V_blks_["abba_oovo"]("k,j,e,i");
    }

    {

        // sharedOps["bbbb_vvoo_72"] += 1.000000 t2_bbbb_vvoo("c,a,i,k") V_blks_["bbbb_vovo"]("b,k,c,j")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_72"]("a,b,i,j") = t2_bbbb_vvoo("c,a,i,k") * V_blks_["bbbb_vovo"]("b,k,c,j");

        // sharedOps["bbbb_vvoo_73"] += 1.000000 t2_abab_vvoo("c,a,k,i") V_blks_["baab_vovo"]("b,k,c,j")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_73"]("a,b,i,j") = t2_abab_vvoo("c,a,k,i") * V_blks_["baab_vovo"]("b,k,c,j");
    }

    if (include_u2_) {

        // sharedOps["bbbb_vvoo_74"] += 1.000000 V_blks_["bbbb_vovo"]("b,k,c,j") u2_bbbb_vvoo("c,a,i,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_74"]("b,a,j,i") = V_blks_["bbbb_vovo"]("b,k,c,j") * u2_bbbb_vvoo("c,a,i,k");

        // sharedOps["abab_vvoo_75"] += 1.000000 u2_abab_vvoo("c,a,i,k") V_blks_["abab_vovo"]("b,k,c,j")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_75"]("b,a,i,j") = u2_abab_vvoo("c,a,i,k") * V_blks_["abab_vovo"]("b,k,c,j");

        // sharedOps["abab_vvoo_76"] += 1.000000 u2_aaaa_vvoo("c,a,i,k") V_blks_["baab_vovo"]("b,k,c,j")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_76"]("a,b,i,j") = u2_aaaa_vvoo("c,a,i,k") * V_blks_["baab_vovo"]("b,k,c,j");
    }

    {

        // sharedOps["abab_vvoo_77"] += 1.000000 V_blks_["baba_vovo"]("b,k,c,j") t2_abab_vvoo("a,c,k,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_77"]("a,b,j,i") = V_blks_["baba_vovo"]("b,k,c,j") * t2_abab_vvoo("a,c,k,i");

        // sharedOps["abab_vvoo_78"] += 1.000000 t2_abab_vvoo("c,a,k,i") V_blks_["aaaa_vovo"]("b,k,c,j")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_78"]("b,a,j,i") = t2_abab_vvoo("c,a,k,i") * V_blks_["aaaa_vovo"]("b,k,c,j");

        // sharedOps["abab_vvoo_79"] += 1.000000 V_blks_["abba_vovo"]("b,k,c,j") t2_bbbb_vvoo("c,a,i,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_79"]("b,a,j,i") = V_blks_["abba_vovo"]("b,k,c,j") * t2_bbbb_vvoo("c,a,i,k");

        // sharedOps["abab_vvoo_80"] += 1.000000 V_blks_["bbbb_vovo"]("b,k,c,j") t2_abab_vvoo("a,c,i,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_80"]("a,b,i,j") = V_blks_["bbbb_vovo"]("b,k,c,j") * t2_abab_vvoo("a,c,i,k");
    }

    if (include_u2_) {

        // sharedOps["aabb_vooo_81"] += 1.000000 u2_abab_vvoo("b,c,i,j") V_blks_["abab_vovv"]("a,m,b,c")
        // flops: o3v3: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_81"]("a,i,j,m") = u2_abab_vvoo("b,c,i,j") * V_blks_["abab_vovv"]("a,m,b,c");

        // sharedOps["bbbb_vvoo_82"] += 1.000000 u2_abab_vvoo("c,a,k,i") V_blks_["baab_vovo"]("b,k,c,j")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_82"]("a,b,i,j") = u2_abab_vvoo("c,a,k,i") * V_blks_["baab_vovo"]("b,k,c,j");

        // sharedOps["abab_vvoo_83"] += 1.000000 u2_abab_vvoo("a,c,k,i") V_blks_["baba_vovo"]("b,k,c,j")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_83"]("a,b,j,i") = u2_abab_vvoo("a,c,k,i") * V_blks_["baba_vovo"]("b,k,c,j");
    }

    {

        // sharedOps["aaaa_vvoo_84"] += 1.000000 V_blks_["aaaa_vovo"]("b,k,c,j") t2_aaaa_vvoo("c,a,i,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_84"]("b,a,j,i") = V_blks_["aaaa_vovo"]("b,k,c,j") * t2_aaaa_vvoo("c,a,i,k");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_85"] += 1.000000 V_blks_["abba_vovo"]("b,k,c,j") u2_bbbb_vvoo("c,a,i,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_85"]("b,a,j,i") = V_blks_["abba_vovo"]("b,k,c,j") * u2_bbbb_vvoo("c,a,i,k");

        // sharedOps["abab_vvoo_86"] += 1.000000 V_blks_["aaaa_vovo"]("b,k,c,j") u2_abab_vvoo("c,a,k,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_86"]("b,a,j,i") = V_blks_["aaaa_vovo"]("b,k,c,j") * u2_abab_vvoo("c,a,k,i");
    }

    {

        // sharedOps["abab_vvoo_87"] += 1.000000 V_blks_["abab_vovo"]("b,k,c,j") t2_abab_vvoo("c,a,i,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_87"]("b,a,i,j") = V_blks_["abab_vovo"]("b,k,c,j") * t2_abab_vvoo("c,a,i,k");

        // sharedOps["aaaa_vvoo_88"] += 1.000000 t2_abab_vvoo("a,c,i,k") V_blks_["abba_vovo"]("b,k,c,j")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_88"]("a,b,i,j") = t2_abab_vvoo("a,c,i,k") * V_blks_["abba_vovo"]("b,k,c,j");

        // sharedOps["abab_vvoo_89"] += 1.000000 V_blks_["baab_vovo"]("b,k,c,j") t2_aaaa_vvoo("c,a,i,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_89"]("a,b,i,j") = V_blks_["baab_vovo"]("b,k,c,j") * t2_aaaa_vvoo("c,a,i,k");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_90"] += 1.000000 u2_abab_vvoo("a,c,i,k") V_blks_["bbbb_vovo"]("b,k,c,j")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_90"]("a,b,i,j") = u2_abab_vvoo("a,c,i,k") * V_blks_["bbbb_vovo"]("b,k,c,j");

        // sharedOps["aaaa_vvoo_91"] += 1.000000 V_blks_["aaaa_vovo"]("b,k,c,j") u2_aaaa_vvoo("c,a,i,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_91"]("b,a,j,i") = V_blks_["aaaa_vovo"]("b,k,c,j") * u2_aaaa_vvoo("c,a,i,k");

        // sharedOps["bbbb_vvvo_92"] += 1.000000 V_blks_["bbbb_oovo"]("k,j,e,i") u2_bbbb_vvoo("b,a,k,j")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["bbbb_vvvo_92"]("e,b,a,i") = V_blks_["bbbb_oovo"]("k,j,e,i") * u2_bbbb_vvoo("b,a,k,j");

        // sharedOps["aabb_vvoo_93"] += 1.000000 u2_abab_vvoo("b,c,k,j") V_blks_["abab_oovv"]("k,l,d,c")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aabb_vvoo_93"]("b,d,j,l") = u2_abab_vvoo("b,c,k,j") * V_blks_["abab_oovv"]("k,l,d,c");

        // sharedOps["aaaa_vvvo_94"] += 1.000000 u2_aaaa_vvoo("b,a,k,j") V_blks_["aaaa_oovo"]("k,j,e,i")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["aaaa_vvvo_94"]("b,a,e,i") = u2_aaaa_vvoo("b,a,k,j") * V_blks_["aaaa_oovo"]("k,j,e,i");

        // sharedOps["bbaa_vvoo_95"] += 1.000000 V_blks_["abab_oovv"]("l,k,c,d") u2_abab_vvoo("c,b,j,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbaa_vvoo_95"]("d,b,l,j") = V_blks_["abab_oovv"]("l,k,c,d") * u2_abab_vvoo("c,b,j,k");

        // sharedOps["aaaa_vooo_96"] += 1.000000 u2_aaaa_vvoo("b,c,i,j") V_blks_["aaaa_vovv"]("a,m,b,c")
        // flops: o3v3: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aaaa_vooo_96"]("a,i,j,m") = u2_aaaa_vvoo("b,c,i,j") * V_blks_["aaaa_vovv"]("a,m,b,c");

        // sharedOps["aaaa_vvoo_97"] += 1.000000 V_blks_["abba_vovo"]("b,k,c,j") u2_abab_vvoo("a,c,i,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_97"]("b,a,j,i") = V_blks_["abba_vovo"]("b,k,c,j") * u2_abab_vvoo("a,c,i,k");

        // sharedOps["aabb_vvvo_98"] += 1.000000 u2_abab_vvoo("b,a,j,k") V_blks_["abab_oovo"]("j,k,e,i")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["aabb_vvvo_98"]("b,e,a,i") = u2_abab_vvoo("b,a,j,k") * V_blks_["abab_oovo"]("j,k,e,i");

        // sharedOps["baab_vooo_99"] += 1.000000 u2_abab_vvoo("b,c,i,j") V_blks_["baab_vovv"]("a,m,b,c")
        // flops: o3v3: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_99"]("a,i,m,j") = u2_abab_vvoo("b,c,i,j") * V_blks_["baab_vovv"]("a,m,b,c");

        // sharedOps["bbbb_vooo_100"] += 1.000000 V_blks_["bbbb_vovv"]("a,m,b,c") u2_bbbb_vvoo("b,c,i,j")
        // flops: o3v3: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["bbbb_vooo_100"]("a,m,i,j") = V_blks_["bbbb_vovv"]("a,m,b,c") * u2_bbbb_vvoo("b,c,i,j");
    }

    {

        // sharedOps["aabb_oooo_103"] += 1.000000 V_blks_["abab_oovv"]("k,l,c,d") t2_abab_vvoo("c,d,j,i")
        // flops: o4v2: 1, o4v0: 1 | mem: o4v0: 2,
        sharedOps["aabb_oooo_103"]("k,j,l,i") = V_blks_["abab_oovv"]("k,l,c,d") * t2_abab_vvoo("c,d,j,i");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_104"] += 1.000000 V_blks_["abab_oovv"]("k,l,d,c") u2_abab_vvoo("d,c,j,i") t2_abab_vvoo("b,a,k,l")
        // flops: o4v2: 2, o2v2: 1 | mem: o2v2: 2, o4v0: 1,
        sharedOps["abab_vvoo_104"]("b,a,j,i") = V_blks_["abab_oovv"]("k,l,d,c") * u2_abab_vvoo("d,c,j,i") * t2_abab_vvoo("b,a,k,l");
    }

    {

        // sharedOps["abab_vvoo_105"] += 1.000000 t2_abab_vvoo("b,a,l,k") V_blks_["abab_oooo"]("l,k,i,j")
        // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_105"]("b,a,i,j") = t2_abab_vvoo("b,a,l,k") * V_blks_["abab_oooo"]("l,k,i,j");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_106"] += 1.000000 V_blks_["abab_oooo"]("l,k,i,j") u2_abab_vvoo("b,a,l,k")
        // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_106"]("b,a,i,j") = V_blks_["abab_oooo"]("l,k,i,j") * u2_abab_vvoo("b,a,l,k");
    }

    {

        // sharedOps["bbbb_oooo_107"] += 1.000000 t2_bbbb_vvoo("c,d,i,j") V_blks_["bbbb_oovv"]("l,k,c,d")
        // flops: o4v2: 1, o4v0: 1 | mem: o4v0: 2,
        sharedOps["bbbb_oooo_107"]("i,j,l,k") = t2_bbbb_vvoo("c,d,i,j") * V_blks_["bbbb_oovv"]("l,k,c,d");

        // sharedOps["aaaa_oooo_108"] += 1.000000 t2_aaaa_vvoo("c,d,i,j") V_blks_["aaaa_oovv"]("l,k,c,d")
        // flops: o4v2: 1, o4v0: 1 | mem: o4v0: 2,
        sharedOps["aaaa_oooo_108"]("i,j,l,k") = t2_aaaa_vvoo("c,d,i,j") * V_blks_["aaaa_oovv"]("l,k,c,d");
    }

    if (include_u1_) {

        // sharedOps["abab_vvoo_109"] += 1.000000 V_blks_["abab_oovo"]("l,k,c,j") u1_aa_vo("c,i") t2_abab_vvoo("a,b,l,k")
        // flops: o4v2: 1, o4v1: 1, o2v2: 1 | mem: o2v2: 2, o4v0: 1,
        sharedOps["abab_vvoo_109"]("a,b,i,j") = V_blks_["abab_oovo"]("l,k,c,j") * u1_aa_vo("c,i") * t2_abab_vvoo("a,b,l,k");

        // sharedOps["abab_vvoo_110"] += 1.000000 V_blks_["abba_oovo"]("k,l,c,j") u1_bb_vo("c,i") t2_abab_vvoo("b,a,k,l")
        // flops: o4v2: 1, o4v1: 1, o2v2: 1 | mem: o2v2: 2, o4v0: 1,
        sharedOps["abab_vvoo_110"]("b,a,j,i") = V_blks_["abba_oovo"]("k,l,c,j") * u1_bb_vo("c,i") * t2_abab_vvoo("b,a,k,l");
    }

    {

        // sharedOps["baab_vooo_111"] += 1.000000 t2_bbbb_vvoo("c,a,i,l") V_blks_["abba_oovo"]("k,l,c,j")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_111"]("a,k,j,i") = t2_bbbb_vvoo("c,a,i,l") * V_blks_["abba_oovo"]("k,l,c,j");

        // sharedOps["bbbb_vooo_112"] += 1.000000 t2_bbbb_vvoo("c,a,i,l") V_blks_["bbbb_oovo"]("l,k,c,j")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["bbbb_vooo_112"]("a,i,k,j") = t2_bbbb_vvoo("c,a,i,l") * V_blks_["bbbb_oovo"]("l,k,c,j");

        // sharedOps["aabb_vooo_113"] += 1.000000 V_blks_["abba_oovo"]("l,k,c,j") t2_abab_vvoo("a,c,l,i")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_113"]("a,j,k,i") = V_blks_["abba_oovo"]("l,k,c,j") * t2_abab_vvoo("a,c,l,i");

        // sharedOps["aabb_vooo_114"] += 1.000000 V_blks_["bbbb_oovo"]("l,k,c,j") t2_abab_vvoo("a,c,i,l")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_114"]("a,i,k,j") = V_blks_["bbbb_oovo"]("l,k,c,j") * t2_abab_vvoo("a,c,i,l");

        // sharedOps["aaaa_vooo_115"] += 1.000000 V_blks_["abba_oovo"]("k,l,c,j") t2_abab_vvoo("a,c,i,l")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aaaa_vooo_115"]("a,k,j,i") = V_blks_["abba_oovo"]("k,l,c,j") * t2_abab_vvoo("a,c,i,l");

        // sharedOps["bbbb_vooo_116"] += 1.000000 t2_abab_vvoo("c,a,l,i") V_blks_["abab_oovo"]("l,k,c,j")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["bbbb_vooo_116"]("a,i,k,j") = t2_abab_vvoo("c,a,l,i") * V_blks_["abab_oovo"]("l,k,c,j");

        // sharedOps["baab_vooo_117"] += 1.000000 V_blks_["aaaa_oovo"]("l,k,c,j") t2_abab_vvoo("c,a,l,i")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_117"]("a,k,j,i") = V_blks_["aaaa_oovo"]("l,k,c,j") * t2_abab_vvoo("c,a,l,i");

        // sharedOps["baab_vooo_118"] += 1.000000 V_blks_["abab_oovo"]("k,l,c,j") t2_abab_vvoo("c,a,i,l")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_118"]("a,k,i,j") = V_blks_["abab_oovo"]("k,l,c,j") * t2_abab_vvoo("c,a,i,l");

        // sharedOps["aaaa_vooo_119"] += 1.000000 V_blks_["aaaa_oovo"]("l,k,c,j") t2_aaaa_vvoo("c,a,i,l")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aaaa_vooo_119"]("a,k,j,i") = V_blks_["aaaa_oovo"]("l,k,c,j") * t2_aaaa_vvoo("c,a,i,l");

        // sharedOps["aabb_vooo_120"] += 1.000000 V_blks_["abab_oovo"]("l,k,c,j") t2_aaaa_vvoo("c,a,i,l")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_120"]("a,i,k,j") = V_blks_["abab_oovo"]("l,k,c,j") * t2_aaaa_vvoo("c,a,i,l");
    }

    if (include_u2_) {

        // sharedOps["aaaa_oooo_121"] += 1.000000 V_blks_["aaaa_oovv"]("l,k,c,d") u2_aaaa_vvoo("c,d,i,j")
        // flops: o4v2: 1, o4v0: 1 | mem: o4v0: 2,
        sharedOps["aaaa_oooo_121"]("l,k,i,j") = V_blks_["aaaa_oovv"]("l,k,c,d") * u2_aaaa_vvoo("c,d,i,j");

        // sharedOps["bbbb_vvoo_122"] += 1.000000 V_blks_["bbbb_oovv"]("l,k,c,d") u2_bbbb_vvoo("c,d,i,j") t2_bbbb_vvoo("b,a,l,k")
        // flops: o4v2: 2, o2v2: 1 | mem: o2v2: 2, o4v0: 1,
        sharedOps["bbbb_vvoo_122"]("b,a,i,j") = V_blks_["bbbb_oovv"]("l,k,c,d") * u2_bbbb_vvoo("c,d,i,j") * t2_bbbb_vvoo("b,a,l,k");
    }

    if (include_u1_) {

        // sharedOps["aaaa_vvoo_123"] += 1.000000 V_blks_["aaaa_oovo"]("l,k,c,j") u1_aa_vo("c,i") t2_aaaa_vvoo("b,a,l,k")
        // flops: o4v2: 1, o4v1: 1, o2v2: 1 | mem: o2v2: 2, o4v0: 1,
        sharedOps["aaaa_vvoo_123"]("b,a,j,i") = V_blks_["aaaa_oovo"]("l,k,c,j") * u1_aa_vo("c,i") * t2_aaaa_vvoo("b,a,l,k");

        // sharedOps["bbbb_vvoo_124"] += 1.000000 V_blks_["bbbb_oovo"]("l,k,c,j") u1_bb_vo("c,i") t2_bbbb_vvoo("b,a,l,k")
        // flops: o4v2: 1, o4v1: 1, o2v2: 1 | mem: o2v2: 2, o4v0: 1,
        sharedOps["bbbb_vvoo_124"]("b,a,j,i") = V_blks_["bbbb_oovo"]("l,k,c,j") * u1_bb_vo("c,i") * t2_bbbb_vvoo("b,a,l,k");
    }

    if (include_u2_) {

        // sharedOps["aabb_oooo_125"] += 1.000000 u2_abab_vvoo("d,c,j,i") V_blks_["abab_oovv"]("k,l,d,c")
        // flops: o4v2: 1, o4v0: 1 | mem: o4v0: 2,
        sharedOps["aabb_oooo_125"]("j,k,i,l") = u2_abab_vvoo("d,c,j,i") * V_blks_["abab_oovv"]("k,l,d,c");

        // sharedOps["baab_vooo_126"] += 1.000000 V_blks_["aaaa_oovo"]("k,m,b,j") u2_abab_vvoo("b,a,k,i")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_126"]("a,m,j,i") = V_blks_["aaaa_oovo"]("k,m,b,j") * u2_abab_vvoo("b,a,k,i");
    }

    {

        // sharedOps["bbbb_vvoo_127"] += 1.000000 V_blks_["bbbb_oooo"]("l,k,i,j") t2_bbbb_vvoo("b,a,l,k")
        // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_127"]("b,a,i,j") = V_blks_["bbbb_oooo"]("l,k,i,j") * t2_bbbb_vvoo("b,a,l,k");
    }

    if (include_u2_) {

        // sharedOps["bbbb_vvoo_128"] += 1.000000 V_blks_["bbbb_oooo"]("l,k,i,j") u2_bbbb_vvoo("b,a,l,k")
        // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_128"]("b,a,i,j") = V_blks_["bbbb_oooo"]("l,k,i,j") * u2_bbbb_vvoo("b,a,l,k");

        // sharedOps["baab_vooo_129"] += 1.000000 u2_bbbb_vvoo("b,a,i,k") V_blks_["abba_oovo"]("m,k,b,j")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_129"]("a,m,j,i") = u2_bbbb_vvoo("b,a,i,k") * V_blks_["abba_oovo"]("m,k,b,j");

        // sharedOps["aabb_vooo_130"] += 1.000000 u2_aaaa_vvoo("b,a,i,k") V_blks_["abab_oovo"]("k,m,b,j")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_130"]("a,i,m,j") = u2_aaaa_vvoo("b,a,i,k") * V_blks_["abab_oovo"]("k,m,b,j");

        // sharedOps["aaaa_vvoo_131"] += 1.000000 V_blks_["aaaa_oooo"]("l,k,i,j") u2_aaaa_vvoo("b,a,l,k")
        // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_131"]("b,a,i,j") = V_blks_["aaaa_oooo"]("l,k,i,j") * u2_aaaa_vvoo("b,a,l,k");

        // sharedOps["baab_vooo_132"] += 1.000000 u2_abab_vvoo("b,a,i,k") V_blks_["abab_oovo"]("m,k,b,j")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_132"]("a,i,m,j") = u2_abab_vvoo("b,a,i,k") * V_blks_["abab_oovo"]("m,k,b,j");

        // sharedOps["bbbb_vooo_133"] += 1.000000 V_blks_["abab_oovo"]("k,m,b,j") u2_abab_vvoo("b,a,k,i")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["bbbb_vooo_133"]("a,m,j,i") = V_blks_["abab_oovo"]("k,m,b,j") * u2_abab_vvoo("b,a,k,i");

        // sharedOps["bbbb_vooo_134"] += 1.000000 V_blks_["bbbb_oovo"]("k,m,b,j") u2_bbbb_vvoo("b,a,i,k")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["bbbb_vooo_134"]("a,m,j,i") = V_blks_["bbbb_oovo"]("k,m,b,j") * u2_bbbb_vvoo("b,a,i,k");

        // sharedOps["bbbb_oooo_135"] += 1.000000 u2_bbbb_vvoo("c,d,i,j") V_blks_["bbbb_oovv"]("l,k,c,d")
        // flops: o4v2: 1, o4v0: 1 | mem: o4v0: 2,
        sharedOps["bbbb_oooo_135"]("i,j,l,k") = u2_bbbb_vvoo("c,d,i,j") * V_blks_["bbbb_oovv"]("l,k,c,d");

        // sharedOps["aabb_vooo_136"] += 1.000000 V_blks_["bbbb_oovo"]("k,m,b,j") u2_abab_vvoo("a,b,i,k")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_136"]("a,i,m,j") = V_blks_["bbbb_oovo"]("k,m,b,j") * u2_abab_vvoo("a,b,i,k");

        // sharedOps["aaaa_vooo_137"] += 1.000000 V_blks_["abba_oovo"]("m,k,b,j") u2_abab_vvoo("a,b,i,k")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aaaa_vooo_137"]("a,m,j,i") = V_blks_["abba_oovo"]("m,k,b,j") * u2_abab_vvoo("a,b,i,k");

        // sharedOps["aaaa_vooo_138"] += 1.000000 u2_aaaa_vvoo("b,a,i,k") V_blks_["aaaa_oovo"]("k,m,b,j")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aaaa_vooo_138"]("a,i,m,j") = u2_aaaa_vvoo("b,a,i,k") * V_blks_["aaaa_oovo"]("k,m,b,j");
    }

    {

        // sharedOps["aaaa_vvoo_139"] += 1.000000 t2_aaaa_vvoo("b,a,l,k") V_blks_["aaaa_oooo"]("l,k,i,j")
        // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_139"]("b,a,i,j") = t2_aaaa_vvoo("b,a,l,k") * V_blks_["aaaa_oooo"]("l,k,i,j");
    }

    if (include_u2_) {

        // sharedOps["aabb_vooo_140"] += 1.000000 V_blks_["abba_oovo"]("k,m,b,j") u2_abab_vvoo("a,b,k,i")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_140"]("a,j,m,i") = V_blks_["abba_oovo"]("k,m,b,j") * u2_abab_vvoo("a,b,k,i");
    }

    if (include_u1_) {

        // sharedOps["aabb_vvvo_141"] += 1.000000 V_blks_["abab_vvvv"]("b,a,e,c") u1_bb_vo("c,i")
        // flops: o1v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["aabb_vvvo_141"]("b,e,a,i") = V_blks_["abab_vvvv"]("b,a,e,c") * u1_bb_vo("c,i");

        // sharedOps["abba_vvvo_142"] += 1.000000 V_blks_["abab_vvvv"]("b,a,c,e") u1_aa_vo("c,i")
        // flops: o1v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["abba_vvvo_142"]("b,a,e,i") = V_blks_["abab_vvvv"]("b,a,c,e") * u1_aa_vo("c,i");

        // sharedOps["aaaa_vvvo_143"] += 1.000000 u1_aa_vo("c,i") V_blks_["aaaa_vvvv"]("b,a,e,c")
        // flops: o1v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["aaaa_vvvo_143"]("b,a,e,i") = u1_aa_vo("c,i") * V_blks_["aaaa_vvvv"]("b,a,e,c");

        // sharedOps["bbbb_vvvo_144"] += 1.000000 u1_bb_vo("c,i") V_blks_["bbbb_vvvv"]("b,a,e,c")
        // flops: o1v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["bbbb_vvvo_144"]("b,a,e,i") = u1_bb_vo("c,i") * V_blks_["bbbb_vvvv"]("b,a,e,c");

        // sharedOps["aabb_vooo_145"] += 1.000000 u1_bb_vo("b,i") V_blks_["abab_oovv"]("k,j,c,b")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_145"]("c,k,i,j") = u1_bb_vo("b,i") * V_blks_["abab_oovv"]("k,j,c,b");

        // sharedOps["baab_vooo_146"] += 1.000000 u1_aa_vo("b,i") V_blks_["abab_oovv"]("k,j,b,c")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_146"]("c,i,k,j") = u1_aa_vo("b,i") * V_blks_["abab_oovv"]("k,j,b,c");

        // sharedOps["aaaa_vooo_147"] += 1.000000 u1_aa_vo("c,j") V_blks_["aaaa_oovv"]("k,m,e,c")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aaaa_vooo_147"]("e,j,k,m") = u1_aa_vo("c,j") * V_blks_["aaaa_oovv"]("k,m,e,c");

        // sharedOps["bbbb_vooo_148"] += 1.000000 u1_bb_vo("c,j") V_blks_["bbbb_oovv"]("k,m,e,c")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["bbbb_vooo_148"]("e,j,k,m") = u1_bb_vo("c,j") * V_blks_["bbbb_oovv"]("k,m,e,c");
    }

    {

        // sharedOps["bb_vv_149"] += 1.000000 V_blks_["abab_oovv"]("k,l,d,c") t2_abab_vvoo("d,a,k,l")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sharedOps["bb_vv_149"]("c,a") = V_blks_["abab_oovv"]("k,l,d,c") * t2_abab_vvoo("d,a,k,l");

        // sharedOps["aa_vv_150"] += 1.000000 t2_abab_vvoo("a,d,l,k") V_blks_["abab_oovv"]("l,k,c,d")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sharedOps["aa_vv_150"]("a,c") = t2_abab_vvoo("a,d,l,k") * V_blks_["abab_oovv"]("l,k,c,d");

        // sharedOps["aa_vv_151"] += 1.000000 t2_aaaa_vvoo("d,a,l,k") V_blks_["aaaa_oovv"]("l,k,c,d")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sharedOps["aa_vv_151"]("a,c") = t2_aaaa_vvoo("d,a,l,k") * V_blks_["aaaa_oovv"]("l,k,c,d");

        // sharedOps["bb_vv_152"] += 1.000000 V_blks_["bbbb_oovv"]("l,k,c,d") t2_bbbb_vvoo("d,a,l,k")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sharedOps["bb_vv_152"]("c,a") = V_blks_["bbbb_oovv"]("l,k,c,d") * t2_bbbb_vvoo("d,a,l,k");
    }

    if (include_u2_) {

        // sharedOps["aa_vv_153"] += 1.000000 u2_abab_vvoo("b,c,k,l") V_blks_["abab_oovv"]("k,l,d,c")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sharedOps["aa_vv_153"]("b,d") = u2_abab_vvoo("b,c,k,l") * V_blks_["abab_oovv"]("k,l,d,c");

        // sharedOps["bb_vv_154"] += 1.000000 u2_abab_vvoo("c,b,l,k") V_blks_["abab_oovv"]("l,k,c,d")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sharedOps["bb_vv_154"]("b,d") = u2_abab_vvoo("c,b,l,k") * V_blks_["abab_oovv"]("l,k,c,d");
    }

    {

        // sharedOps["abab_vvoo_155"] += 1.000000 dp_aa_vv("b,c") t2_abab_vvoo("c,a,i,j")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_155"]("b,a,i,j") = dp_aa_vv("b,c") * t2_abab_vvoo("c,a,i,j");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_156"] += 1.000000 u2_abab_vvoo("c,b,j,i") dp_bb_vv("a,b")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_156"]("c,a,j,i") = u2_abab_vvoo("c,b,j,i") * dp_bb_vv("a,b");

        // sharedOps["abab_vvoo_157"] += 1.000000 dp_aa_vv("a,b") u2_abab_vvoo("b,c,j,i")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_157"]("a,c,j,i") = dp_aa_vv("a,b") * u2_abab_vvoo("b,c,j,i");
    }

    {

        // sharedOps["abab_vvoo_158"] += 1.000000 dp_bb_vv("b,c") t2_abab_vvoo("a,c,i,j")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_158"]("a,b,i,j") = dp_bb_vv("b,c") * t2_abab_vvoo("a,c,i,j");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_159"] += 1.000000 V_blks_["aaaa_oovv"]("l,k,c,d") u2_aaaa_vvoo("c,b,l,k") t2_abab_vvoo("d,a,j,i")
        // flops: o2v3: 2, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
        sharedOps["abab_vvoo_159"]("b,a,j,i") = V_blks_["aaaa_oovv"]("l,k,c,d") * u2_aaaa_vvoo("c,b,l,k") * t2_abab_vvoo("d,a,j,i");

        // sharedOps["abab_vvoo_160"] += 1.000000 V_blks_["bbbb_oovv"]("l,k,c,d") u2_bbbb_vvoo("c,b,l,k") t2_abab_vvoo("a,d,i,j")
        // flops: o2v3: 2, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
        sharedOps["abab_vvoo_160"]("a,b,i,j") = V_blks_["bbbb_oovv"]("l,k,c,d") * u2_bbbb_vvoo("c,b,l,k") * t2_abab_vvoo("a,d,i,j");
    }

    {

        // sharedOps["abab_vvoo_161"] += 1.000000 V_blks_["aaaa_oovv"]("l,k,c,d") t2_aaaa_vvoo("c,b,l,k") t2_abab_vvoo("d,a,i,j")
        // flops: o2v3: 2, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
        sharedOps["abab_vvoo_161"]("b,a,i,j") = V_blks_["aaaa_oovv"]("l,k,c,d") * t2_aaaa_vvoo("c,b,l,k") * t2_abab_vvoo("d,a,i,j");

        // sharedOps["abab_vvoo_162"] += 1.000000 V_blks_["bbbb_oovv"]("l,k,c,d") t2_bbbb_vvoo("c,b,l,k") t2_abab_vvoo("a,d,i,j")
        // flops: o2v3: 2, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
        sharedOps["abab_vvoo_162"]("a,b,i,j") = V_blks_["bbbb_oovv"]("l,k,c,d") * t2_bbbb_vvoo("c,b,l,k") * t2_abab_vvoo("a,d,i,j");
    }

    if (include_u1_) {

        // sharedOps["abab_vvoo_163"] += 1.000000 V_blks_["baab_vovv"]("b,k,d,c") u1_bb_vo("c,j")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_163"]("d,b,k,j") = V_blks_["baab_vovv"]("b,k,d,c") * u1_bb_vo("c,j");
    }

    if (include_u2_) {

        // sharedOps["aaaa_vvoo_164"] += 1.000000 dp_aa_vv("a,b") u2_aaaa_vvoo("b,c,j,i")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_164"]("a,c,j,i") = dp_aa_vv("a,b") * u2_aaaa_vvoo("b,c,j,i");
    }

    if (include_u1_) {

        // sharedOps["aaaa_vvoo_165"] += 1.000000 V_blks_["aaaa_vovv"]("a,m,e,b") u1_aa_vo("b,i")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_165"]("a,e,m,i") = V_blks_["aaaa_vovv"]("a,m,e,b") * u1_aa_vo("b,i");

        // sharedOps["bbbb_vvoo_166"] += 1.000000 u1_bb_vo("b,i") V_blks_["bbbb_vovv"]("a,m,e,b")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_166"]("a,e,i,m") = u1_bb_vo("b,i") * V_blks_["bbbb_vovv"]("a,m,e,b");

        // sharedOps["abab_vvoo_167"] += 1.000000 V_blks_["abab_vovv"]("b,k,c,d") u1_aa_vo("c,j")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_167"]("b,d,j,k") = V_blks_["abab_vovv"]("b,k,c,d") * u1_aa_vo("c,j");
    }

    {

        // sharedOps["bbbb_vvoo_168"] += 1.000000 dp_bb_vv("b,c") t2_bbbb_vvoo("c,a,i,j")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_168"]("b,a,i,j") = dp_bb_vv("b,c") * t2_bbbb_vvoo("c,a,i,j");
    }

    if (include_u2_) {

        // sharedOps["bb_vv_169"] += 1.000000 u2_bbbb_vvoo("b,a,j,i") V_blks_["bbbb_oovv"]("j,i,e,b")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sharedOps["bb_vv_169"]("a,e") = u2_bbbb_vvoo("b,a,j,i") * V_blks_["bbbb_oovv"]("j,i,e,b");

        // sharedOps["bbbb_vvoo_170"] += 1.000000 dp_bb_vv("a,b") u2_bbbb_vvoo("b,c,j,i")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_170"]("a,c,j,i") = dp_bb_vv("a,b") * u2_bbbb_vvoo("b,c,j,i");

        // sharedOps["aa_vv_171"] += 1.000000 u2_aaaa_vvoo("b,a,j,i") V_blks_["aaaa_oovv"]("j,i,e,b")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sharedOps["aa_vv_171"]("a,e") = u2_aaaa_vvoo("b,a,j,i") * V_blks_["aaaa_oovv"]("j,i,e,b");
    }

    {

        // sharedOps["aaaa_vvoo_172"] += 1.000000 dp_aa_vv("b,c") t2_aaaa_vvoo("c,a,i,j")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_172"]("b,a,i,j") = dp_aa_vv("b,c") * t2_aaaa_vvoo("c,a,i,j");
    }

    if (include_u1_) {

        // sharedOps["abab_vvoo_173"] += 1.000000 V_blks_["bbbb_vovv"]("b,k,c,d") u1_bb_vo("c,k") t2_abab_vvoo("a,d,i,j")
        // flops: o2v3: 1, o1v3: 1, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
        sharedOps["abab_vvoo_173"]("a,b,i,j") = V_blks_["bbbb_vovv"]("b,k,c,d") * u1_bb_vo("c,k") * t2_abab_vvoo("a,d,i,j");

        // sharedOps["abab_vvoo_174"] += 1.000000 V_blks_["baab_vovv"]("b,k,c,d") u1_aa_vo("c,k") t2_abab_vvoo("a,d,i,j")
        // flops: o2v3: 1, o1v3: 1, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
        sharedOps["abab_vvoo_174"]("a,b,i,j") = V_blks_["baab_vovv"]("b,k,c,d") * u1_aa_vo("c,k") * t2_abab_vvoo("a,d,i,j");

        // sharedOps["abab_vvoo_175"] += 1.000000 V_blks_["aaaa_vovv"]("b,k,c,d") u1_aa_vo("c,k") t2_abab_vvoo("d,a,i,j")
        // flops: o2v3: 1, o1v3: 1, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
        sharedOps["abab_vvoo_175"]("b,a,i,j") = V_blks_["aaaa_vovv"]("b,k,c,d") * u1_aa_vo("c,k") * t2_abab_vvoo("d,a,i,j");

        // sharedOps["abab_vvoo_176"] += 1.000000 V_blks_["abab_vovv"]("b,k,d,c") u1_bb_vo("c,k") t2_abab_vvoo("d,a,j,i")
        // flops: o2v3: 1, o1v3: 1, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
        sharedOps["abab_vvoo_176"]("b,a,j,i") = V_blks_["abab_vovv"]("b,k,d,c") * u1_bb_vo("c,k") * t2_abab_vvoo("d,a,j,i");
    }

    if (include_u2_) {

        // sharedOps["aa_vo_177"] += 1.000000 V_blks_["abab_vovv"]("a,j,b,c") u2_abab_vvoo("b,c,i,j")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_177"]("a,i") = V_blks_["abab_vovv"]("a,j,b,c") * u2_abab_vvoo("b,c,i,j");
    }

    {

        // sharedOps["bbbb_vvoo_178"] += 1.000000 V_blks_["bbbb_oovv"]("l,k,c,d") t2_bbbb_vvoo("c,b,l,k") t2_bbbb_vvoo("d,a,i,j")
        // flops: o2v3: 2, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
        sharedOps["bbbb_vvoo_178"]("b,a,i,j") = V_blks_["bbbb_oovv"]("l,k,c,d") * t2_bbbb_vvoo("c,b,l,k") * t2_bbbb_vvoo("d,a,i,j");
    }

    if (include_u2_) {

        // sharedOps["bbbb_vvoo_179"] += 1.000000 V_blks_["bbbb_oovv"]("l,k,c,d") u2_bbbb_vvoo("c,b,l,k") t2_bbbb_vvoo("d,a,i,j")
        // flops: o2v3: 2, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
        sharedOps["bbbb_vvoo_179"]("b,a,i,j") = V_blks_["bbbb_oovv"]("l,k,c,d") * u2_bbbb_vvoo("c,b,l,k") * t2_bbbb_vvoo("d,a,i,j");

        // sharedOps["abab_vvoo_180"] += 1.000000 u2_abab_vvoo("a,c,i,j") F_blks_["bb_vv"]("b,c")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_180"]("a,b,i,j") = u2_abab_vvoo("a,c,i,j") * F_blks_["bb_vv"]("b,c");
    }

    {

        // sharedOps["abab_vvoo_181"] += 1.000000 F_blks_["aa_vv"]("b,c") t2_abab_vvoo("c,a,i,j")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_181"]("b,a,i,j") = F_blks_["aa_vv"]("b,c") * t2_abab_vvoo("c,a,i,j");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_182"] += 1.000000 F_blks_["aa_vv"]("b,c") u2_abab_vvoo("c,a,i,j")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_182"]("b,a,i,j") = F_blks_["aa_vv"]("b,c") * u2_abab_vvoo("c,a,i,j");
    }

    {

        // sharedOps["bb_vo_183"] += 1.000000 V_blks_["baab_vovv"]("a,j,b,c") t2_abab_vvoo("b,c,j,i")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_183"]("a,i") = V_blks_["baab_vovv"]("a,j,b,c") * t2_abab_vvoo("b,c,j,i");

        // sharedOps["aa_vo_184"] += 1.000000 V_blks_["abab_vovv"]("a,j,b,c") t2_abab_vvoo("b,c,i,j")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_184"]("a,i") = V_blks_["abab_vovv"]("a,j,b,c") * t2_abab_vvoo("b,c,i,j");
    }

    if (include_u1_) {

        // sharedOps["abab_vvoo_185"] += 1.000000 u1_aa_vo("c,i") V_blks_["abab_vvvo"]("b,a,c,j")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_185"]("b,a,i,j") = u1_aa_vo("c,i") * V_blks_["abab_vvvo"]("b,a,c,j");
    }

    {

        // sharedOps["aaaa_vvoo_186"] += 1.000000 V_blks_["aaaa_oovv"]("l,k,c,d") t2_aaaa_vvoo("c,b,l,k") t2_aaaa_vvoo("d,a,i,j")
        // flops: o2v3: 2, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
        sharedOps["aaaa_vvoo_186"]("b,a,i,j") = V_blks_["aaaa_oovv"]("l,k,c,d") * t2_aaaa_vvoo("c,b,l,k") * t2_aaaa_vvoo("d,a,i,j");
    }

    if (include_u2_) {

        // sharedOps["aaaa_vvoo_187"] += 1.000000 V_blks_["aaaa_oovv"]("l,k,c,d") u2_aaaa_vvoo("c,b,l,k") t2_aaaa_vvoo("d,a,i,j")
        // flops: o2v3: 2, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
        sharedOps["aaaa_vvoo_187"]("b,a,i,j") = V_blks_["aaaa_oovv"]("l,k,c,d") * u2_aaaa_vvoo("c,b,l,k") * t2_aaaa_vvoo("d,a,i,j");

        // sharedOps["bb_vo_188"] += 1.000000 V_blks_["baab_vovv"]("a,j,b,c") u2_abab_vvoo("b,c,j,i")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_188"]("a,i") = V_blks_["baab_vovv"]("a,j,b,c") * u2_abab_vvoo("b,c,j,i");
    }

    {

        // sharedOps["abab_vvoo_189"] += 1.000000 t2_abab_vvoo("a,c,i,j") F_blks_["bb_vv"]("b,c")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_189"]("a,b,i,j") = t2_abab_vvoo("a,c,i,j") * F_blks_["bb_vv"]("b,c");
    }

    if (include_u1_) {

        // sharedOps["abab_vvoo_190"] += 1.000000 V_blks_["abba_vvvo"]("b,a,c,j") u1_bb_vo("c,i")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_190"]("b,a,j,i") = V_blks_["abba_vvvo"]("b,a,c,j") * u1_bb_vo("c,i");

        // sharedOps["aaaa_vvoo_191"] += 1.000000 V_blks_["abab_vovv"]("b,k,d,c") u1_bb_vo("c,k") t2_aaaa_vvoo("d,a,i,j")
        // flops: o2v3: 1, o1v3: 1, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
        sharedOps["aaaa_vvoo_191"]("b,a,i,j") = V_blks_["abab_vovv"]("b,k,d,c") * u1_bb_vo("c,k") * t2_aaaa_vvoo("d,a,i,j");

        // sharedOps["aaaa_vvoo_192"] += 1.000000 V_blks_["aaaa_vovv"]("b,k,c,d") u1_aa_vo("c,k") t2_aaaa_vvoo("d,a,i,j")
        // flops: o2v3: 1, o1v3: 1, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
        sharedOps["aaaa_vvoo_192"]("b,a,i,j") = V_blks_["aaaa_vovv"]("b,k,c,d") * u1_aa_vo("c,k") * t2_aaaa_vvoo("d,a,i,j");

        // sharedOps["bbbb_vvoo_193"] += 1.000000 V_blks_["bbbb_vovv"]("b,k,c,d") u1_bb_vo("c,k") t2_bbbb_vvoo("d,a,i,j")
        // flops: o2v3: 1, o1v3: 1, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
        sharedOps["bbbb_vvoo_193"]("b,a,i,j") = V_blks_["bbbb_vovv"]("b,k,c,d") * u1_bb_vo("c,k") * t2_bbbb_vvoo("d,a,i,j");

        // sharedOps["bbbb_vvoo_194"] += 1.000000 V_blks_["baab_vovv"]("b,k,c,d") u1_aa_vo("c,k") t2_bbbb_vvoo("d,a,i,j")
        // flops: o2v3: 1, o1v3: 1, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
        sharedOps["bbbb_vvoo_194"]("b,a,i,j") = V_blks_["baab_vovv"]("b,k,c,d") * u1_aa_vo("c,k") * t2_bbbb_vvoo("d,a,i,j");

        // sharedOps["aaaa_vvoo_195"] += 1.000000 V_blks_["aaaa_vvvo"]("b,a,c,j") u1_aa_vo("c,i")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_195"]("b,a,j,i") = V_blks_["aaaa_vvvo"]("b,a,c,j") * u1_aa_vo("c,i");
    }

    if (include_u2_) {

        // sharedOps["bbbb_vvoo_196"] += 1.000000 F_blks_["bb_vv"]("b,c") u2_bbbb_vvoo("c,a,i,j")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_196"]("b,a,i,j") = F_blks_["bb_vv"]("b,c") * u2_bbbb_vvoo("c,a,i,j");
    }

    {

        // sharedOps["bbbb_vvoo_197"] += 1.000000 t2_bbbb_vvoo("c,a,i,j") F_blks_["bb_vv"]("b,c")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_197"]("a,b,i,j") = t2_bbbb_vvoo("c,a,i,j") * F_blks_["bb_vv"]("b,c");

        // sharedOps["aa_vo_198"] += 1.000000 t2_aaaa_vvoo("b,c,i,j") V_blks_["aaaa_vovv"]("a,j,b,c")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_198"]("a,i") = t2_aaaa_vvoo("b,c,i,j") * V_blks_["aaaa_vovv"]("a,j,b,c");
    }

    if (include_u2_) {

        // sharedOps["aa_vo_199"] += 1.000000 u2_aaaa_vvoo("b,c,i,j") V_blks_["aaaa_vovv"]("a,j,b,c")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_199"]("a,i") = u2_aaaa_vvoo("b,c,i,j") * V_blks_["aaaa_vovv"]("a,j,b,c");

        // sharedOps["aaaa_vvoo_200"] += 1.000000 F_blks_["aa_vv"]("b,c") u2_aaaa_vvoo("c,a,i,j")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_200"]("b,a,i,j") = F_blks_["aa_vv"]("b,c") * u2_aaaa_vvoo("c,a,i,j");
    }

    if (include_u1_) {

        // sharedOps["bbaa_vvoo_201"] += 1.000000 u1_aa_vo("c,j") V_blks_["baab_vovv"]("b,k,c,d")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbaa_vvoo_201"]("b,d,j,k") = u1_aa_vo("c,j") * V_blks_["baab_vovv"]("b,k,c,d");
    }

    if (include_u2_) {

        // sharedOps["bb_vo_202"] += 1.000000 V_blks_["bbbb_vovv"]("a,j,b,c") u2_bbbb_vvoo("b,c,i,j")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_202"]("a,i") = V_blks_["bbbb_vovv"]("a,j,b,c") * u2_bbbb_vvoo("b,c,i,j");
    }

    if (include_u1_) {

        // sharedOps["bbbb_vvoo_203"] += 1.000000 V_blks_["bbbb_vvvo"]("b,a,c,j") u1_bb_vo("c,i")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_203"]("b,a,j,i") = V_blks_["bbbb_vvvo"]("b,a,c,j") * u1_bb_vo("c,i");
    }

    {

        // sharedOps["bb_vo_204"] += 1.000000 t2_bbbb_vvoo("b,c,i,j") V_blks_["bbbb_vovv"]("a,j,b,c")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_204"]("a,i") = t2_bbbb_vvoo("b,c,i,j") * V_blks_["bbbb_vovv"]("a,j,b,c");
    }

    if (include_u1_) {

        // sharedOps["aabb_vvoo_205"] += 1.000000 u1_bb_vo("c,j") V_blks_["abab_vovv"]("b,k,d,c")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aabb_vvoo_205"]("b,d,j,k") = u1_bb_vo("c,j") * V_blks_["abab_vovv"]("b,k,d,c");
    }

    {

        // sharedOps["aaaa_vvoo_206"] += 1.000000 t2_aaaa_vvoo("c,a,i,j") F_blks_["aa_vv"]("b,c")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_206"]("a,b,i,j") = t2_aaaa_vvoo("c,a,i,j") * F_blks_["aa_vv"]("b,c");

        // sharedOps["aa_oo_207"] += 1.000000 V_blks_["abab_oovv"]("j,k,b,c") t2_abab_vvoo("b,c,i,k")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sharedOps["aa_oo_207"]("j,i") = V_blks_["abab_oovv"]("j,k,b,c") * t2_abab_vvoo("b,c,i,k");

        // sharedOps["bb_oo_208"] += 1.000000 V_blks_["abab_oovv"]("k,j,b,c") t2_abab_vvoo("b,c,k,i")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sharedOps["bb_oo_208"]("j,i") = V_blks_["abab_oovv"]("k,j,b,c") * t2_abab_vvoo("b,c,k,i");
    }

    if (include_u1_) {

        // sharedOps["abab_vvoo_209"] += 1.000000 dp_bb_ov("k,c") t2_abab_vvoo("a,c,i,j") u1_bb_vo("b,k")
        // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
        sharedOps["abab_vvoo_209"]("a,b,i,j") = dp_bb_ov("k,c") * t2_abab_vvoo("a,c,i,j") * u1_bb_vo("b,k");

        // sharedOps["abab_vvoo_210"] += 1.000000 dp_aa_ov("k,c") t2_abab_vvoo("c,a,i,j") u1_aa_vo("b,k")
        // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
        sharedOps["abab_vvoo_210"]("b,a,i,j") = dp_aa_ov("k,c") * t2_abab_vvoo("c,a,i,j") * u1_aa_vo("b,k");
    }

    {

        // sharedOps["aa_oo_211"] += 1.000000 t2_aaaa_vvoo("b,c,i,k") V_blks_["aaaa_oovv"]("k,j,b,c")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sharedOps["aa_oo_211"]("i,j") = t2_aaaa_vvoo("b,c,i,k") * V_blks_["aaaa_oovv"]("k,j,b,c");

        // sharedOps["bb_oo_212"] += 1.000000 V_blks_["bbbb_oovv"]("k,j,b,c") t2_bbbb_vvoo("b,c,i,k")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sharedOps["bb_oo_212"]("j,i") = V_blks_["bbbb_oovv"]("k,j,b,c") * t2_bbbb_vvoo("b,c,i,k");
    }

    if (include_u2_) {

        // sharedOps["aa_oo_213"] += 1.000000 V_blks_["abab_oovv"]("l,k,c,d") u2_abab_vvoo("c,d,j,k")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sharedOps["aa_oo_213"]("l,j") = V_blks_["abab_oovv"]("l,k,c,d") * u2_abab_vvoo("c,d,j,k");

        // sharedOps["bb_oo_214"] += 1.000000 u2_abab_vvoo("d,c,k,j") V_blks_["abab_oovv"]("k,l,d,c")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sharedOps["bb_oo_214"]("j,l") = u2_abab_vvoo("d,c,k,j") * V_blks_["abab_oovv"]("k,l,d,c");
    }

    {

        // sharedOps["aaaa_vooo_215"] += 1.000000 dp_aa_ov("k,c") t2_aaaa_vvoo("c,a,i,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aaaa_vooo_215"]("a,k,i,j") = dp_aa_ov("k,c") * t2_aaaa_vvoo("c,a,i,j");

        // sharedOps["bbbb_vooo_216"] += 1.000000 t2_bbbb_vvoo("c,a,i,j") dp_bb_ov("k,c")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["bbbb_vooo_216"]("a,i,j,k") = t2_bbbb_vvoo("c,a,i,j") * dp_bb_ov("k,c");
    }

    if (include_u1_) {

        // sharedOps["abab_vvoo_217"] += 1.000000 dp_aa_ov("k,c") u1_aa_vo("c,j") t2_abab_vvoo("b,a,k,i")
        // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["abab_vvoo_217"]("b,a,j,i") = dp_aa_ov("k,c") * u1_aa_vo("c,j") * t2_abab_vvoo("b,a,k,i");

        // sharedOps["abab_vvoo_218"] += 1.000000 dp_bb_ov("k,c") u1_bb_vo("c,j") t2_abab_vvoo("b,a,i,k")
        // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["abab_vvoo_218"]("b,a,i,j") = dp_bb_ov("k,c") * u1_bb_vo("c,j") * t2_abab_vvoo("b,a,i,k");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_219"] += 1.000000 dp_bb_oo("j,i") u2_abab_vvoo("a,b,k,j")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_219"]("a,b,k,i") = dp_bb_oo("j,i") * u2_abab_vvoo("a,b,k,j");
    }

    {

        // sharedOps["abab_vvoo_220"] += 1.000000 dp_aa_oo("k,j") t2_abab_vvoo("b,a,k,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_220"]("b,a,j,i") = dp_aa_oo("k,j") * t2_abab_vvoo("b,a,k,i");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_221"] += 1.000000 dp_aa_oo("j,i") u2_abab_vvoo("a,b,j,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_221"]("a,b,i,k") = dp_aa_oo("j,i") * u2_abab_vvoo("a,b,j,k");
    }

    {

        // sharedOps["abab_vvoo_222"] += 1.000000 t2_abab_vvoo("b,a,i,k") dp_bb_oo("k,j")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_222"]("b,a,i,j") = t2_abab_vvoo("b,a,i,k") * dp_bb_oo("k,j");
    }

    if (include_u2_) {

        // sharedOps["baab_vooo_223"] += 1.000000 u2_abab_vvoo("c,a,i,j") dp_aa_ov("k,c")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_223"]("a,i,k,j") = u2_abab_vvoo("c,a,i,j") * dp_aa_ov("k,c");

        // sharedOps["aabb_vooo_224"] += 1.000000 u2_abab_vvoo("a,c,i,j") dp_bb_ov("k,c")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_224"]("a,i,j,k") = u2_abab_vvoo("a,c,i,j") * dp_bb_ov("k,c");
    }

    if (include_u1_) {

        // sharedOps["abab_vvoo_225"] += 1.000000 F_blks_["aa_ov"]("k,c") t2_abab_vvoo("c,a,i,j") u1_aa_vo("b,k")
        // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
        sharedOps["abab_vvoo_225"]("b,a,i,j") = F_blks_["aa_ov"]("k,c") * t2_abab_vvoo("c,a,i,j") * u1_aa_vo("b,k");
    }

    {

        // sharedOps["baab_vooo_226"] += 1.000000 t2_abab_vvoo("c,a,i,j") dp_aa_ov("k,c")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_226"]("a,i,k,j") = t2_abab_vvoo("c,a,i,j") * dp_aa_ov("k,c");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_227"] += 1.000000 V_blks_["aaaa_oovv"]("l,k,c,d") u2_aaaa_vvoo("c,d,j,k") t2_abab_vvoo("a,b,l,i")
        // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["abab_vvoo_227"]("a,b,j,i") = V_blks_["aaaa_oovv"]("l,k,c,d") * u2_aaaa_vvoo("c,d,j,k") * t2_abab_vvoo("a,b,l,i");
    }

    {

        // sharedOps["aabb_vooo_228"] += 1.000000 t2_abab_vvoo("a,c,i,j") dp_bb_ov("k,c")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_228"]("a,i,j,k") = t2_abab_vvoo("a,c,i,j") * dp_bb_ov("k,c");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_229"] += 1.000000 V_blks_["bbbb_oovv"]("l,k,c,d") u2_bbbb_vvoo("c,d,j,k") t2_abab_vvoo("a,b,i,l")
        // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["abab_vvoo_229"]("a,b,i,j") = V_blks_["bbbb_oovv"]("l,k,c,d") * u2_bbbb_vvoo("c,d,j,k") * t2_abab_vvoo("a,b,i,l");
    }

    if (include_u1_) {

        // sharedOps["abab_vvoo_230"] += 1.000000 F_blks_["bb_ov"]("k,c") t2_abab_vvoo("a,c,i,j") u1_bb_vo("b,k")
        // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
        sharedOps["abab_vvoo_230"]("a,b,i,j") = F_blks_["bb_ov"]("k,c") * t2_abab_vvoo("a,c,i,j") * u1_bb_vo("b,k");
    }

    {

        // sharedOps["abab_vvoo_231"] += 1.000000 V_blks_["aaaa_oovv"]("l,k,c,d") t2_aaaa_vvoo("c,d,j,k") t2_abab_vvoo("b,a,l,i")
        // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["abab_vvoo_231"]("b,a,j,i") = V_blks_["aaaa_oovv"]("l,k,c,d") * t2_aaaa_vvoo("c,d,j,k") * t2_abab_vvoo("b,a,l,i");
    }

    if (include_u2_) {

        // sharedOps["bbbb_vooo_232"] += 1.000000 dp_bb_ov("k,c") u2_bbbb_vvoo("c,a,i,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["bbbb_vooo_232"]("a,k,i,j") = dp_bb_ov("k,c") * u2_bbbb_vvoo("c,a,i,j");
    }

    {

        // sharedOps["abab_vvoo_233"] += 1.000000 V_blks_["bbbb_oovv"]("l,k,c,d") t2_bbbb_vvoo("c,d,j,k") t2_abab_vvoo("b,a,i,l")
        // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["abab_vvoo_233"]("b,a,i,j") = V_blks_["bbbb_oovv"]("l,k,c,d") * t2_bbbb_vvoo("c,d,j,k") * t2_abab_vvoo("b,a,i,l");
    }

    if (include_u2_) {

        // sharedOps["aaaa_vooo_234"] += 1.000000 u2_aaaa_vvoo("c,a,i,j") dp_aa_ov("k,c")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aaaa_vooo_234"]("a,i,j,k") = u2_aaaa_vvoo("c,a,i,j") * dp_aa_ov("k,c");
    }

    if (include_u1_) {

        // sharedOps["aaaa_vvoo_235"] += 1.000000 dp_aa_ov("k,c") u1_aa_vo("c,j") t2_aaaa_vvoo("b,a,i,k")
        // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["aaaa_vvoo_235"]("b,a,j,i") = dp_aa_ov("k,c") * u1_aa_vo("c,j") * t2_aaaa_vvoo("b,a,i,k");

        // sharedOps["bbbb_vvoo_236"] += 1.000000 dp_bb_ov("k,c") u1_bb_vo("c,j") t2_bbbb_vvoo("b,a,i,k")
        // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["bbbb_vvoo_236"]("b,a,j,i") = dp_bb_ov("k,c") * u1_bb_vo("c,j") * t2_bbbb_vvoo("b,a,i,k");
    }

    if (include_u2_) {

        // sharedOps["aaaa_vvoo_237"] += 1.000000 dp_aa_oo("j,i") u2_aaaa_vvoo("a,b,k,j")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_237"]("a,b,i,k") = dp_aa_oo("j,i") * u2_aaaa_vvoo("a,b,k,j");
    }

    {

        // sharedOps["bbbb_vvoo_238"] += 1.000000 dp_bb_oo("k,j") t2_bbbb_vvoo("b,a,i,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_238"]("b,a,j,i") = dp_bb_oo("k,j") * t2_bbbb_vvoo("b,a,i,k");

        // sharedOps["aaaa_vvoo_239"] += 1.000000 t2_aaaa_vvoo("b,a,i,k") dp_aa_oo("k,j")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_239"]("b,a,i,j") = t2_aaaa_vvoo("b,a,i,k") * dp_aa_oo("k,j");

        // sharedOps["aaaa_vooo_240"] += 1.000000 F_blks_["aa_ov"]("k,c") t2_aaaa_vvoo("c,a,i,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aaaa_vooo_240"]("a,k,i,j") = F_blks_["aa_ov"]("k,c") * t2_aaaa_vvoo("c,a,i,j");
    }

    if (include_u2_) {

        // sharedOps["bbbb_vvoo_241"] += 1.000000 u2_bbbb_vvoo("a,b,k,j") dp_bb_oo("j,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_241"]("a,b,k,i") = u2_bbbb_vvoo("a,b,k,j") * dp_bb_oo("j,i");

        // sharedOps["bb_oo_242"] += 1.000000 u2_bbbb_vvoo("a,b,i,j") V_blks_["bbbb_oovv"]("j,m,a,b")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sharedOps["bb_oo_242"]("i,m") = u2_bbbb_vvoo("a,b,i,j") * V_blks_["bbbb_oovv"]("j,m,a,b");
    }

    {

        // sharedOps["bbbb_vooo_243"] += 1.000000 F_blks_["bb_ov"]("k,c") t2_bbbb_vvoo("c,a,i,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["bbbb_vooo_243"]("a,k,i,j") = F_blks_["bb_ov"]("k,c") * t2_bbbb_vvoo("c,a,i,j");
    }

    if (include_u2_) {

        // sharedOps["aa_oo_244"] += 1.000000 u2_aaaa_vvoo("a,b,i,j") V_blks_["aaaa_oovv"]("j,m,a,b")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sharedOps["aa_oo_244"]("i,m") = u2_aaaa_vvoo("a,b,i,j") * V_blks_["aaaa_oovv"]("j,m,a,b");
    }

    if (include_u1_) {

        // sharedOps["abab_vvoo_245"] += 1.000000 V_blks_["abab_oovo"]("k,l,c,j") u1_aa_vo("c,k") t2_abab_vvoo("b,a,i,l")
        // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["abab_vvoo_245"]("b,a,i,j") = V_blks_["abab_oovo"]("k,l,c,j") * u1_aa_vo("c,k") * t2_abab_vvoo("b,a,i,l");

        // sharedOps["abab_vvoo_246"] += 1.000000 V_blks_["bbbb_oovo"]("l,k,c,j") u1_bb_vo("c,k") t2_abab_vvoo("b,a,i,l")
        // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["abab_vvoo_246"]("b,a,i,j") = V_blks_["bbbb_oovo"]("l,k,c,j") * u1_bb_vo("c,k") * t2_abab_vvoo("b,a,i,l");

        // sharedOps["abab_vvoo_247"] += 1.000000 V_blks_["aaaa_oovo"]("l,k,c,j") u1_aa_vo("c,k") t2_abab_vvoo("b,a,l,i")
        // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["abab_vvoo_247"]("b,a,j,i") = V_blks_["aaaa_oovo"]("l,k,c,j") * u1_aa_vo("c,k") * t2_abab_vvoo("b,a,l,i");

        // sharedOps["abab_vvoo_248"] += 1.000000 V_blks_["abba_oovo"]("l,k,c,j") u1_bb_vo("c,k") t2_abab_vvoo("b,a,l,i")
        // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["abab_vvoo_248"]("b,a,j,i") = V_blks_["abba_oovo"]("l,k,c,j") * u1_bb_vo("c,k") * t2_abab_vvoo("b,a,l,i");
    }

    if (include_u1_ && include_u2_) {

        // sharedOps["abab_vvoo_249"] += 1.000000 dp_aa_ov("k,c") u1_aa_vo("c,j") u2_abab_vvoo("b,a,k,i")
        // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["abab_vvoo_249"]("b,a,j,i") = dp_aa_ov("k,c") * u1_aa_vo("c,j") * u2_abab_vvoo("b,a,k,i");
    }

    if (include_u1_) {

        // sharedOps["abab_vvoo_250"] += 1.000000 F_blks_["bb_ov"]("k,c") u1_bb_vo("c,j") t2_abab_vvoo("b,a,i,k")
        // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["abab_vvoo_250"]("b,a,i,j") = F_blks_["bb_ov"]("k,c") * u1_bb_vo("c,j") * t2_abab_vvoo("b,a,i,k");

        // sharedOps["abab_vvoo_251"] += 1.000000 F_blks_["aa_ov"]("k,c") u1_aa_vo("c,j") t2_abab_vvoo("a,b,k,i")
        // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["abab_vvoo_251"]("a,b,j,i") = F_blks_["aa_ov"]("k,c") * u1_aa_vo("c,j") * t2_abab_vvoo("a,b,k,i");
    }

    if (include_u1_ && include_u2_) {

        // sharedOps["abab_vvoo_252"] += 1.000000 dp_bb_ov("k,c") u1_bb_vo("c,j") u2_abab_vvoo("b,a,i,k")
        // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["abab_vvoo_252"]("b,a,i,j") = dp_bb_ov("k,c") * u1_bb_vo("c,j") * u2_abab_vvoo("b,a,i,k");
    }

    if (include_u1_) {

        // sharedOps["aa_vo_253"] += 1.000000 V_blks_["aaaa_oovv"]("k,j,b,c") u1_aa_vo("b,i") t2_aaaa_vvoo("c,a,k,j")
        // flops: o3v2: 2, o1v1: 1 | mem: o3v1: 1, o1v1: 2,
        sharedOps["aa_vo_253"]("a,i") = V_blks_["aaaa_oovv"]("k,j,b,c") * u1_aa_vo("b,i") * t2_aaaa_vvoo("c,a,k,j");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_254"] += 1.000000 F_blks_["bb_oo"]("k,j") u2_abab_vvoo("b,a,i,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_254"]("b,a,i,j") = F_blks_["bb_oo"]("k,j") * u2_abab_vvoo("b,a,i,k");

        // sharedOps["aa_vo_255"] += 1.000000 u2_abab_vvoo("a,b,k,j") V_blks_["abba_oovo"]("k,j,b,i")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_255"]("a,i") = u2_abab_vvoo("a,b,k,j") * V_blks_["abba_oovo"]("k,j,b,i");
    }

    {

        // sharedOps["bbbb_vvoo_256"] += 1.000000 V_blks_["bbbb_oovv"]("l,k,c,d") t2_bbbb_vvoo("c,d,j,k") t2_bbbb_vvoo("b,a,i,l")
        // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["bbbb_vvoo_256"]("b,a,j,i") = V_blks_["bbbb_oovv"]("l,k,c,d") * t2_bbbb_vvoo("c,d,j,k") * t2_bbbb_vvoo("b,a,i,l");

        // sharedOps["aaaa_vvoo_257"] += 1.000000 V_blks_["aaaa_oovv"]("l,k,c,d") t2_aaaa_vvoo("c,d,j,k") t2_aaaa_vvoo("b,a,i,l")
        // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["aaaa_vvoo_257"]("b,a,j,i") = V_blks_["aaaa_oovv"]("l,k,c,d") * t2_aaaa_vvoo("c,d,j,k") * t2_aaaa_vvoo("b,a,i,l");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_258"] += 1.000000 F_blks_["aa_oo"]("k,j") u2_abab_vvoo("b,a,k,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_258"]("b,a,j,i") = F_blks_["aa_oo"]("k,j") * u2_abab_vvoo("b,a,k,i");
    }

    {

        // sharedOps["baab_vooo_259"] += 1.000000 F_blks_["aa_ov"]("k,c") t2_abab_vvoo("c,a,i,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_259"]("a,k,i,j") = F_blks_["aa_ov"]("k,c") * t2_abab_vvoo("c,a,i,j");
    }

    if (include_u1_) {

        // sharedOps["abab_vvoo_260"] += 1.000000 V_blks_["baab_vooo"]("b,k,i,j") u1_aa_vo("a,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_260"]("a,b,i,j") = V_blks_["baab_vooo"]("b,k,i,j") * u1_aa_vo("a,k");

        // sharedOps["bb_vo_261"] += 1.000000 V_blks_["bbbb_oovv"]("k,j,b,c") u1_bb_vo("b,i") t2_bbbb_vvoo("c,a,k,j")
        // flops: o3v2: 2, o1v1: 1 | mem: o3v1: 1, o1v1: 2,
        sharedOps["bb_vo_261"]("a,i") = V_blks_["bbbb_oovv"]("k,j,b,c") * u1_bb_vo("b,i") * t2_bbbb_vvoo("c,a,k,j");
    }

    {

        // sharedOps["aabb_vooo_262"] += 1.000000 t2_abab_vvoo("a,c,i,j") F_blks_["bb_ov"]("k,c")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_262"]("a,i,j,k") = t2_abab_vvoo("a,c,i,j") * F_blks_["bb_ov"]("k,c");
    }

    if (include_u1_) {

        // sharedOps["abab_vvoo_263"] += 1.000000 u1_bb_vo("a,k") V_blks_["abab_vooo"]("b,k,i,j")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_263"]("b,a,i,j") = u1_bb_vo("a,k") * V_blks_["abab_vooo"]("b,k,i,j");
    }

    if (include_u2_) {

        // sharedOps["bb_vo_264"] += 1.000000 u2_abab_vvoo("b,a,k,j") V_blks_["abab_oovo"]("k,j,b,i")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_264"]("a,i") = u2_abab_vvoo("b,a,k,j") * V_blks_["abab_oovo"]("k,j,b,i");
    }

    {

        // sharedOps["aa_vo_265"] += 1.000000 V_blks_["abba_oovo"]("k,j,b,i") t2_abab_vvoo("a,b,k,j")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_265"]("a,i") = V_blks_["abba_oovo"]("k,j,b,i") * t2_abab_vvoo("a,b,k,j");

        // sharedOps["bb_vo_266"] += 1.000000 t2_abab_vvoo("b,a,k,j") V_blks_["abab_oovo"]("k,j,b,i")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_266"]("a,i") = t2_abab_vvoo("b,a,k,j") * V_blks_["abab_oovo"]("k,j,b,i");
    }

    if (include_u2_) {

        // sharedOps["bbbb_vvoo_267"] += 1.000000 V_blks_["bbbb_oovv"]("l,k,c,d") u2_bbbb_vvoo("c,d,j,k") t2_bbbb_vvoo("b,a,i,l")
        // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["bbbb_vvoo_267"]("b,a,j,i") = V_blks_["bbbb_oovv"]("l,k,c,d") * u2_bbbb_vvoo("c,d,j,k") * t2_bbbb_vvoo("b,a,i,l");
    }

    {

        // sharedOps["abab_vvoo_268"] += 1.000000 t2_abab_vvoo("b,a,k,i") F_blks_["aa_oo"]("k,j")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_268"]("b,a,j,i") = t2_abab_vvoo("b,a,k,i") * F_blks_["aa_oo"]("k,j");

        // sharedOps["abab_vvoo_269"] += 1.000000 F_blks_["bb_oo"]("k,j") t2_abab_vvoo("b,a,i,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_269"]("b,a,i,j") = F_blks_["bb_oo"]("k,j") * t2_abab_vvoo("b,a,i,k");
    }

    if (include_u2_) {

        // sharedOps["aaaa_vvoo_270"] += 1.000000 V_blks_["aaaa_oovv"]("l,k,c,d") u2_aaaa_vvoo("c,d,j,k") t2_aaaa_vvoo("b,a,i,l")
        // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["aaaa_vvoo_270"]("b,a,j,i") = V_blks_["aaaa_oovv"]("l,k,c,d") * u2_aaaa_vvoo("c,d,j,k") * t2_aaaa_vvoo("b,a,i,l");
    }

    if (include_u1_) {

        // sharedOps["baab_vooo_271"] += 1.000000 V_blks_["aaaa_oovv"]("k,m,b,c") u1_aa_vo("b,k") t2_abab_vvoo("c,a,j,i")
        // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o3v1: 2, o1v1: 1,
        sharedOps["baab_vooo_271"]("a,m,j,i") = V_blks_["aaaa_oovv"]("k,m,b,c") * u1_aa_vo("b,k") * t2_abab_vvoo("c,a,j,i");

        // sharedOps["aabb_vooo_272"] += 1.000000 V_blks_["bbbb_oovv"]("k,m,b,c") u1_bb_vo("b,k") t2_abab_vvoo("a,c,i,j")
        // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o3v1: 2, o1v1: 1,
        sharedOps["aabb_vooo_272"]("a,i,m,j") = V_blks_["bbbb_oovv"]("k,m,b,c") * u1_bb_vo("b,k") * t2_abab_vvoo("a,c,i,j");

        // sharedOps["aaaa_vooo_273"] += 1.000000 V_blks_["aaaa_oovv"]("k,m,b,c") u1_aa_vo("b,k") t2_aaaa_vvoo("c,a,i,j")
        // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o3v1: 2, o1v1: 1,
        sharedOps["aaaa_vooo_273"]("a,m,i,j") = V_blks_["aaaa_oovv"]("k,m,b,c") * u1_aa_vo("b,k") * t2_aaaa_vvoo("c,a,i,j");

        // sharedOps["bbbb_vooo_274"] += 1.000000 V_blks_["bbbb_oovv"]("k,m,b,c") u1_bb_vo("b,k") t2_bbbb_vvoo("c,a,i,j")
        // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o3v1: 2, o1v1: 1,
        sharedOps["bbbb_vooo_274"]("a,m,i,j") = V_blks_["bbbb_oovv"]("k,m,b,c") * u1_bb_vo("b,k") * t2_bbbb_vvoo("c,a,i,j");

        // sharedOps["aaaa_vvoo_275"] += 1.000000 V_blks_["abba_oovo"]("l,k,c,j") u1_bb_vo("c,k") t2_aaaa_vvoo("b,a,i,l")
        // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["aaaa_vvoo_275"]("b,a,j,i") = V_blks_["abba_oovo"]("l,k,c,j") * u1_bb_vo("c,k") * t2_aaaa_vvoo("b,a,i,l");

        // sharedOps["bbbb_vvoo_276"] += 1.000000 V_blks_["bbbb_oovo"]("l,k,c,j") u1_bb_vo("c,k") t2_bbbb_vvoo("b,a,i,l")
        // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["bbbb_vvoo_276"]("b,a,j,i") = V_blks_["bbbb_oovo"]("l,k,c,j") * u1_bb_vo("c,k") * t2_bbbb_vvoo("b,a,i,l");

        // sharedOps["aaaa_vvoo_277"] += 1.000000 V_blks_["aaaa_oovo"]("l,k,c,j") u1_aa_vo("c,k") t2_aaaa_vvoo("b,a,i,l")
        // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["aaaa_vvoo_277"]("b,a,j,i") = V_blks_["aaaa_oovo"]("l,k,c,j") * u1_aa_vo("c,k") * t2_aaaa_vvoo("b,a,i,l");

        // sharedOps["bbbb_vvoo_278"] += 1.000000 V_blks_["abab_oovo"]("k,l,c,j") u1_aa_vo("c,k") t2_bbbb_vvoo("b,a,i,l")
        // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["bbbb_vvoo_278"]("b,a,j,i") = V_blks_["abab_oovo"]("k,l,c,j") * u1_aa_vo("c,k") * t2_bbbb_vvoo("b,a,i,l");
    }

    if (include_u1_ && include_u2_) {

        // sharedOps["bbbb_vvoo_279"] += 1.000000 dp_bb_ov("k,c") u1_bb_vo("c,j") u2_bbbb_vvoo("b,a,i,k")
        // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["bbbb_vvoo_279"]("b,a,j,i") = dp_bb_ov("k,c") * u1_bb_vo("c,j") * u2_bbbb_vvoo("b,a,i,k");

        // sharedOps["aaaa_vvoo_280"] += 1.000000 dp_aa_ov("k,c") u1_aa_vo("c,j") u2_aaaa_vvoo("b,a,i,k")
        // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["aaaa_vvoo_280"]("b,a,j,i") = dp_aa_ov("k,c") * u1_aa_vo("c,j") * u2_aaaa_vvoo("b,a,i,k");
    }

    if (include_u1_) {

        // sharedOps["bbbb_vvoo_281"] += 1.000000 F_blks_["bb_ov"]("k,c") u1_bb_vo("c,j") t2_bbbb_vvoo("b,a,i,k")
        // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["bbbb_vvoo_281"]("b,a,j,i") = F_blks_["bb_ov"]("k,c") * u1_bb_vo("c,j") * t2_bbbb_vvoo("b,a,i,k");

        // sharedOps["aaaa_vvoo_282"] += 1.000000 F_blks_["aa_ov"]("k,c") u1_aa_vo("c,j") t2_aaaa_vvoo("b,a,i,k")
        // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
        sharedOps["aaaa_vvoo_282"]("b,a,j,i") = F_blks_["aa_ov"]("k,c") * u1_aa_vo("c,j") * t2_aaaa_vvoo("b,a,i,k");
    }

    {

        // sharedOps["aaaa_vvoo_283"] += 1.000000 t2_aaaa_vvoo("b,a,i,k") F_blks_["aa_oo"]("k,j")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_283"]("b,a,i,j") = t2_aaaa_vvoo("b,a,i,k") * F_blks_["aa_oo"]("k,j");
    }

    if (include_u1_) {

        // sharedOps["baab_vooo_284"] += 1.000000 V_blks_["baab_vovo"]("a,m,b,j") u1_aa_vo("b,i")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_284"]("a,m,i,j") = V_blks_["baab_vovo"]("a,m,b,j") * u1_aa_vo("b,i");
    }

    if (include_u2_) {

        // sharedOps["aaaa_vvoo_285"] += 1.000000 u2_aaaa_vvoo("b,a,i,k") F_blks_["aa_oo"]("k,j")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_285"]("b,a,i,j") = u2_aaaa_vvoo("b,a,i,k") * F_blks_["aa_oo"]("k,j");
    }

    {

        // sharedOps["bb_vo_286"] += 1.000000 t2_bbbb_vvoo("b,a,k,j") V_blks_["bbbb_oovo"]("k,j,b,i")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_286"]("a,i") = t2_bbbb_vvoo("b,a,k,j") * V_blks_["bbbb_oovo"]("k,j,b,i");
    }

    if (include_u2_) {

        // sharedOps["aaaa_vooo_287"] += 1.000000 F_blks_["aa_ov"]("m,b") u2_aaaa_vvoo("b,a,i,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aaaa_vooo_287"]("a,m,i,j") = F_blks_["aa_ov"]("m,b") * u2_aaaa_vvoo("b,a,i,j");
    }

    if (include_u1_) {

        // sharedOps["aaaa_vvoo_288"] += 1.000000 V_blks_["aaaa_vooo"]("b,k,i,j") u1_aa_vo("a,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_288"]("b,a,i,j") = V_blks_["aaaa_vooo"]("b,k,i,j") * u1_aa_vo("a,k");
    }

    if (include_u2_) {

        // sharedOps["bbbb_vvoo_289"] += 1.000000 F_blks_["bb_oo"]("k,j") u2_bbbb_vvoo("b,a,i,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_289"]("b,a,j,i") = F_blks_["bb_oo"]("k,j") * u2_bbbb_vvoo("b,a,i,k");

        // sharedOps["aa_vo_290"] += 1.000000 V_blks_["aaaa_oovo"]("k,j,b,i") u2_aaaa_vvoo("b,a,k,j")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_290"]("a,i") = V_blks_["aaaa_oovo"]("k,j,b,i") * u2_aaaa_vvoo("b,a,k,j");
    }

    {

        // sharedOps["aa_vo_291"] += 1.000000 t2_aaaa_vvoo("b,a,k,j") V_blks_["aaaa_oovo"]("k,j,b,i")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_291"]("a,i") = t2_aaaa_vvoo("b,a,k,j") * V_blks_["aaaa_oovo"]("k,j,b,i");

        // sharedOps["bbbb_vvoo_292"] += 1.000000 t2_bbbb_vvoo("b,a,i,k") F_blks_["bb_oo"]("k,j")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_292"]("b,a,i,j") = t2_bbbb_vvoo("b,a,i,k") * F_blks_["bb_oo"]("k,j");
    }

    if (include_u1_) {

        // sharedOps["bbbb_vvoo_293"] += 1.000000 V_blks_["bbbb_vooo"]("b,k,i,j") u1_bb_vo("a,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_293"]("b,a,i,j") = V_blks_["bbbb_vooo"]("b,k,i,j") * u1_bb_vo("a,k");

        // sharedOps["aabb_vooo_294"] += 1.000000 u1_aa_vo("b,i") V_blks_["abab_vovo"]("a,m,b,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_294"]("a,i,m,j") = u1_aa_vo("b,i") * V_blks_["abab_vovo"]("a,m,b,j");

        // sharedOps["aabb_vooo_295"] += 1.000000 V_blks_["abba_vovo"]("a,m,b,j") u1_bb_vo("b,i")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_295"]("a,j,m,i") = V_blks_["abba_vovo"]("a,m,b,j") * u1_bb_vo("b,i");

        // sharedOps["baab_vooo_296"] += 1.000000 u1_bb_vo("b,i") V_blks_["baba_vovo"]("a,m,b,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_296"]("a,m,j,i") = u1_bb_vo("b,i") * V_blks_["baba_vovo"]("a,m,b,j");
    }

    if (include_u2_) {

        // sharedOps["bb_vo_297"] += 1.000000 u2_bbbb_vvoo("b,a,k,j") V_blks_["bbbb_oovo"]("k,j,b,i")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_297"]("a,i") = u2_bbbb_vvoo("b,a,k,j") * V_blks_["bbbb_oovo"]("k,j,b,i");

        // sharedOps["aabb_vooo_298"] += 1.000000 F_blks_["bb_ov"]("m,b") u2_abab_vvoo("a,b,i,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_298"]("a,i,m,j") = F_blks_["bb_ov"]("m,b") * u2_abab_vvoo("a,b,i,j");
    }

    if (include_u1_) {

        // sharedOps["aaaa_vooo_299"] += 1.000000 V_blks_["aaaa_vovo"]("a,m,b,j") u1_aa_vo("b,i")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aaaa_vooo_299"]("a,m,j,i") = V_blks_["aaaa_vovo"]("a,m,b,j") * u1_aa_vo("b,i");
    }

    if (include_u2_) {

        // sharedOps["baab_vooo_300"] += 1.000000 u2_abab_vvoo("b,a,i,j") F_blks_["aa_ov"]("m,b")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_300"]("a,i,m,j") = u2_abab_vvoo("b,a,i,j") * F_blks_["aa_ov"]("m,b");
    }

    if (include_u1_) {

        // sharedOps["bbbb_vooo_301"] += 1.000000 u1_bb_vo("b,i") V_blks_["bbbb_vovo"]("a,m,b,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["bbbb_vooo_301"]("a,i,m,j") = u1_bb_vo("b,i") * V_blks_["bbbb_vovo"]("a,m,b,j");
    }

    if (include_u2_) {

        // sharedOps["bbbb_vooo_302"] += 1.000000 u2_bbbb_vvoo("b,a,i,j") F_blks_["bb_ov"]("m,b")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["bbbb_vooo_302"]("a,i,j,m") = u2_bbbb_vvoo("b,a,i,j") * F_blks_["bb_ov"]("m,b");
    }

    if (include_u1_) {

        // sharedOps["aaaa_vooo_304"] += 1.000000 V_blks_["aaaa_oooo"]("k,m,i,j") u1_aa_vo("a,k")
        // flops: o4v1: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aaaa_vooo_304"]("a,m,i,j") = V_blks_["aaaa_oooo"]("k,m,i,j") * u1_aa_vo("a,k");

        // sharedOps["baab_vooo_305"] += 1.000000 V_blks_["abab_oooo"]("m,k,i,j") u1_bb_vo("a,k")
        // flops: o4v1: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_305"]("a,m,i,j") = V_blks_["abab_oooo"]("m,k,i,j") * u1_bb_vo("a,k");

        // sharedOps["aabb_oooo_306"] += 1.000000 V_blks_["abab_oovo"]("l,k,c,j") u1_aa_vo("c,i")
        // flops: o4v1: 1, o4v0: 1 | mem: o4v0: 2,
        sharedOps["aabb_oooo_306"]("l,i,k,j") = V_blks_["abab_oovo"]("l,k,c,j") * u1_aa_vo("c,i");

        // sharedOps["bbbb_oooo_307"] += 1.000000 u1_bb_vo("c,i") V_blks_["bbbb_oovo"]("l,k,c,j")
        // flops: o4v1: 1, o4v0: 1 | mem: o4v0: 2,
        sharedOps["bbbb_oooo_307"]("i,l,k,j") = u1_bb_vo("c,i") * V_blks_["bbbb_oovo"]("l,k,c,j");

        // sharedOps["aaaa_oooo_308"] += 1.000000 V_blks_["aaaa_oovo"]("l,k,c,j") u1_aa_vo("c,i")
        // flops: o4v1: 1, o4v0: 1 | mem: o4v0: 2,
        sharedOps["aaaa_oooo_308"]("l,k,j,i") = V_blks_["aaaa_oovo"]("l,k,c,j") * u1_aa_vo("c,i");

        // sharedOps["aabb_oooo_309"] += 1.000000 V_blks_["abba_oovo"]("k,l,c,j") u1_bb_vo("c,i")
        // flops: o4v1: 1, o4v0: 1 | mem: o4v0: 2,
        sharedOps["aabb_oooo_309"]("k,j,l,i") = V_blks_["abba_oovo"]("k,l,c,j") * u1_bb_vo("c,i");

        // sharedOps["bbbb_vooo_310"] += 1.000000 V_blks_["bbbb_oooo"]("k,m,i,j") u1_bb_vo("a,k")
        // flops: o4v1: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["bbbb_vooo_310"]("a,m,i,j") = V_blks_["bbbb_oooo"]("k,m,i,j") * u1_bb_vo("a,k");

        // sharedOps["aabb_vooo_311"] += 1.000000 V_blks_["abab_oooo"]("k,m,i,j") u1_aa_vo("a,k")
        // flops: o4v1: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_311"]("a,i,m,j") = V_blks_["abab_oooo"]("k,m,i,j") * u1_aa_vo("a,k");

        // sharedOps["aa_oo_312"] += 1.000000 u1_aa_vo("c,j") dp_aa_ov("k,c")
        // flops: o2v1: 1, o2v0: 1 | mem: o2v0: 2,
        sharedOps["aa_oo_312"]("j,k") = u1_aa_vo("c,j") * dp_aa_ov("k,c");

        // sharedOps["bb_oo_313"] += 1.000000 dp_bb_ov("k,c") u1_bb_vo("c,j")
        // flops: o2v1: 1, o2v0: 1 | mem: o2v0: 2,
        sharedOps["bb_oo_313"]("k,j") = dp_bb_ov("k,c") * u1_bb_vo("c,j");

        // sharedOps["bb_vv_314"] += 1.000000 u1_bb_vo("b,i") V_blks_["bbbb_vovv"]("a,i,e,b")
        // flops: o1v3: 1, o0v2: 1 | mem: o0v2: 2,
        sharedOps["bb_vv_314"]("a,e") = u1_bb_vo("b,i") * V_blks_["bbbb_vovv"]("a,i,e,b");

        // sharedOps["bb_vv_315"] += 1.000000 u1_aa_vo("c,k") V_blks_["baab_vovv"]("b,k,c,d")
        // flops: o1v3: 1, o0v2: 1 | mem: o0v2: 2,
        sharedOps["bb_vv_315"]("b,d") = u1_aa_vo("c,k") * V_blks_["baab_vovv"]("b,k,c,d");

        // sharedOps["aa_vv_316"] += 1.000000 u1_bb_vo("c,k") V_blks_["abab_vovv"]("b,k,d,c")
        // flops: o1v3: 1, o0v2: 1 | mem: o0v2: 2,
        sharedOps["aa_vv_316"]("b,d") = u1_bb_vo("c,k") * V_blks_["abab_vovv"]("b,k,d,c");

        // sharedOps["aa_vv_317"] += 1.000000 u1_aa_vo("b,i") V_blks_["aaaa_vovv"]("a,i,e,b")
        // flops: o1v3: 1, o0v2: 1 | mem: o0v2: 2,
        sharedOps["aa_vv_317"]("a,e") = u1_aa_vo("b,i") * V_blks_["aaaa_vovv"]("a,i,e,b");
    }

    {

        // sharedOps["aa_vo_318"] += 1.000000 dp_bb_ov("j,b") t2_abab_vvoo("a,b,i,j")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_318"]("a,i") = dp_bb_ov("j,b") * t2_abab_vvoo("a,b,i,j");
    }

    if (include_u2_) {

        // sharedOps["bb_vo_319"] += 1.000000 u2_abab_vvoo("c,a,k,i") dp_aa_ov("k,c")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_319"]("a,i") = u2_abab_vvoo("c,a,k,i") * dp_aa_ov("k,c");
    }

    {

        // sharedOps["aa_vo_320"] += 1.000000 dp_aa_ov("j,b") t2_aaaa_vvoo("b,a,i,j")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_320"]("a,i") = dp_aa_ov("j,b") * t2_aaaa_vvoo("b,a,i,j");

        // sharedOps["bb_vo_321"] += 1.000000 t2_abab_vvoo("b,a,j,i") dp_aa_ov("j,b")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_321"]("a,i") = t2_abab_vvoo("b,a,j,i") * dp_aa_ov("j,b");
    }

    if (include_u2_) {

        // sharedOps["bb_vo_322"] += 1.000000 u2_bbbb_vvoo("c,a,i,k") dp_bb_ov("k,c")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_322"]("a,i") = u2_bbbb_vvoo("c,a,i,k") * dp_bb_ov("k,c");
    }

    {

        // sharedOps["bb_vo_323"] += 1.000000 dp_bb_ov("j,b") t2_bbbb_vvoo("b,a,i,j")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_323"]("a,i") = dp_bb_ov("j,b") * t2_bbbb_vvoo("b,a,i,j");
    }

    if (include_u2_) {

        // sharedOps["aa_vo_324"] += 1.000000 dp_aa_ov("k,c") u2_aaaa_vvoo("c,a,i,k")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_324"]("a,i") = dp_aa_ov("k,c") * u2_aaaa_vvoo("c,a,i,k");

        // sharedOps["aa_vo_325"] += 1.000000 u2_abab_vvoo("a,c,i,k") dp_bb_ov("k,c")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_325"]("a,i") = u2_abab_vvoo("a,c,i,k") * dp_bb_ov("k,c");
    }

    if (include_u1_) {

        // sharedOps["bb_vo_326"] += 1.000000 V_blks_["bbbb_oovv"]("i,m,e,a") u1_bb_vo("a,i")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_326"]("e,m") = V_blks_["bbbb_oovv"]("i,m,e,a") * u1_bb_vo("a,i");

        // sharedOps["aa_vo_327"] += 1.000000 V_blks_["aaaa_oovv"]("i,m,e,a") u1_aa_vo("a,i")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_327"]("e,m") = V_blks_["aaaa_oovv"]("i,m,e,a") * u1_aa_vo("a,i");

        // sharedOps["bb_vo_328"] += 1.000000 V_blks_["bbbb_oovv"]("k,j,b,c") u1_bb_vo("b,j") t2_bbbb_vvoo("c,a,i,k")
        // flops: o2v2: 2, o1v1: 1 | mem: o1v1: 3,
        sharedOps["bb_vo_328"]("a,i") = V_blks_["bbbb_oovv"]("k,j,b,c") * u1_bb_vo("b,j") * t2_bbbb_vvoo("c,a,i,k");

        // sharedOps["aa_vo_329"] += 1.000000 V_blks_["bbbb_oovv"]("k,j,b,c") u1_bb_vo("b,j") t2_abab_vvoo("a,c,i,k")
        // flops: o2v2: 2, o1v1: 1 | mem: o1v1: 3,
        sharedOps["aa_vo_329"]("a,i") = V_blks_["bbbb_oovv"]("k,j,b,c") * u1_bb_vo("b,j") * t2_abab_vvoo("a,c,i,k");

        // sharedOps["bb_vo_330"] += 1.000000 V_blks_["aaaa_oovv"]("k,j,b,c") u1_aa_vo("b,j") t2_abab_vvoo("c,a,k,i")
        // flops: o2v2: 2, o1v1: 1 | mem: o1v1: 3,
        sharedOps["bb_vo_330"]("a,i") = V_blks_["aaaa_oovv"]("k,j,b,c") * u1_aa_vo("b,j") * t2_abab_vvoo("c,a,k,i");

        // sharedOps["aa_vo_331"] += 1.000000 V_blks_["aaaa_oovv"]("k,j,b,c") u1_aa_vo("b,j") t2_aaaa_vvoo("c,a,i,k")
        // flops: o2v2: 2, o1v1: 1 | mem: o1v1: 3,
        sharedOps["aa_vo_331"]("a,i") = V_blks_["aaaa_oovv"]("k,j,b,c") * u1_aa_vo("b,j") * t2_aaaa_vvoo("c,a,i,k");
    }

    {

        // sharedOps["bb_vo_332"] += 1.000000 t2_abab_vvoo("b,a,j,i") F_blks_["aa_ov"]("j,b")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_332"]("a,i") = t2_abab_vvoo("b,a,j,i") * F_blks_["aa_ov"]("j,b");

        // sharedOps["bb_vo_333"] += 1.000000 t2_bbbb_vvoo("b,a,i,j") F_blks_["bb_ov"]("j,b")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_333"]("a,i") = t2_bbbb_vvoo("b,a,i,j") * F_blks_["bb_ov"]("j,b");
    }

    if (include_u1_) {

        // sharedOps["aa_vo_334"] += 1.000000 V_blks_["aaaa_vovo"]("a,j,b,i") u1_aa_vo("b,j")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_334"]("a,i") = V_blks_["aaaa_vovo"]("a,j,b,i") * u1_aa_vo("b,j");
    }

    if (include_u2_) {

        // sharedOps["bb_vo_335"] += 1.000000 u2_bbbb_vvoo("b,a,i,j") F_blks_["bb_ov"]("j,b")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_335"]("a,i") = u2_bbbb_vvoo("b,a,i,j") * F_blks_["bb_ov"]("j,b");
    }

    {

        // sharedOps["aa_vo_336"] += 1.000000 F_blks_["bb_ov"]("j,b") t2_abab_vvoo("a,b,i,j")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_336"]("a,i") = F_blks_["bb_ov"]("j,b") * t2_abab_vvoo("a,b,i,j");

        // sharedOps["aa_vo_337"] += 1.000000 F_blks_["aa_ov"]("j,b") t2_aaaa_vvoo("b,a,i,j")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_337"]("a,i") = F_blks_["aa_ov"]("j,b") * t2_aaaa_vvoo("b,a,i,j");
    }

    if (include_u2_) {

        // sharedOps["bb_vo_338"] += 1.000000 F_blks_["aa_ov"]("j,b") u2_abab_vvoo("b,a,j,i")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_338"]("a,i") = F_blks_["aa_ov"]("j,b") * u2_abab_vvoo("b,a,j,i");
    }

    if (include_u1_) {

        // sharedOps["aa_vo_339"] += 1.000000 u1_bb_vo("b,j") V_blks_["abba_vovo"]("a,j,b,i")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_339"]("a,i") = u1_bb_vo("b,j") * V_blks_["abba_vovo"]("a,j,b,i");

        // sharedOps["bb_vo_340"] += 1.000000 V_blks_["baab_vovo"]("a,j,b,i") u1_aa_vo("b,j")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_340"]("a,i") = V_blks_["baab_vovo"]("a,j,b,i") * u1_aa_vo("b,j");
    }

    if (include_u2_) {

        // sharedOps["aa_vo_341"] += 1.000000 u2_aaaa_vvoo("b,a,i,j") F_blks_["aa_ov"]("j,b")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_341"]("a,i") = u2_aaaa_vvoo("b,a,i,j") * F_blks_["aa_ov"]("j,b");

        // sharedOps["aa_vo_342"] += 1.000000 F_blks_["bb_ov"]("j,b") u2_abab_vvoo("a,b,i,j")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_342"]("a,i") = F_blks_["bb_ov"]("j,b") * u2_abab_vvoo("a,b,i,j");
    }

    if (include_u1_) {

        // sharedOps["bb_vo_343"] += 1.000000 u1_bb_vo("b,j") V_blks_["bbbb_vovo"]("a,j,b,i")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_343"]("a,i") = u1_bb_vo("b,j") * V_blks_["bbbb_vovo"]("a,j,b,i");

        // sharedOps["bb_oo_344"] += 1.000000 u1_bb_vo("a,j") V_blks_["bbbb_oovo"]("j,m,a,i")
        // flops: o3v1: 1, o2v0: 1 | mem: o2v0: 2,
        sharedOps["bb_oo_344"]("m,i") = u1_bb_vo("a,j") * V_blks_["bbbb_oovo"]("j,m,a,i");

        // sharedOps["aa_oo_345"] += 1.000000 V_blks_["abba_oovo"]("l,k,c,j") u1_bb_vo("c,k")
        // flops: o3v1: 1, o2v0: 1 | mem: o2v0: 2,
        sharedOps["aa_oo_345"]("l,j") = V_blks_["abba_oovo"]("l,k,c,j") * u1_bb_vo("c,k");

        // sharedOps["bb_oo_346"] += 1.000000 V_blks_["abab_oovo"]("k,l,c,j") u1_aa_vo("c,k")
        // flops: o3v1: 1, o2v0: 1 | mem: o2v0: 2,
        sharedOps["bb_oo_346"]("l,j") = V_blks_["abab_oovo"]("k,l,c,j") * u1_aa_vo("c,k");

        // sharedOps["aa_oo_347"] += 1.000000 u1_aa_vo("a,j") V_blks_["aaaa_oovo"]("j,m,a,i")
        // flops: o3v1: 1, o2v0: 1 | mem: o2v0: 2,
        sharedOps["aa_oo_347"]("m,i") = u1_aa_vo("a,j") * V_blks_["aaaa_oovo"]("j,m,a,i");

        // sharedOps["bbbb_vvoo_348"] += 1.000000 dp_bb_vv("b,c") u1_bb_vo("c,j") u1_bb_vo("a,i")
        // flops: o2v2: 2, o1v2: 1 | mem: o2v2: 2, o1v1: 1,
        sharedOps["bbbb_vvoo_348"]("b,a,j,i") = dp_bb_vv("b,c") * u1_bb_vo("c,j") * u1_bb_vo("a,i");

        // sharedOps["abab_vvoo_349"] += 1.000000 dp_aa_vv("b,c") u1_aa_vo("c,j") u1_bb_vo("a,i")
        // flops: o2v2: 2, o1v2: 1 | mem: o2v2: 2, o1v1: 1,
        sharedOps["abab_vvoo_349"]("b,a,j,i") = dp_aa_vv("b,c") * u1_aa_vo("c,j") * u1_bb_vo("a,i");

        // sharedOps["aaaa_vvoo_350"] += 1.000000 dp_aa_vv("b,c") u1_aa_vo("c,j") u1_aa_vo("a,i")
        // flops: o2v2: 2, o1v2: 1 | mem: o2v2: 2, o1v1: 1,
        sharedOps["aaaa_vvoo_350"]("b,a,j,i") = dp_aa_vv("b,c") * u1_aa_vo("c,j") * u1_aa_vo("a,i");

        // sharedOps["abab_vvoo_351"] += 1.000000 dp_bb_vv("b,c") u1_bb_vo("c,j") u1_aa_vo("a,i")
        // flops: o2v2: 2, o1v2: 1 | mem: o2v2: 2, o1v1: 1,
        sharedOps["abab_vvoo_351"]("a,b,i,j") = dp_bb_vv("b,c") * u1_bb_vo("c,j") * u1_aa_vo("a,i");

        // sharedOps["bbbb_vvoo_352"] += 1.000000 dp_bb_oo("k,j") u1_bb_vo("b,k") u1_bb_vo("a,i")
        // flops: o2v2: 2, o2v1: 1 | mem: o2v2: 2, o1v1: 1,
        sharedOps["bbbb_vvoo_352"]("b,a,j,i") = dp_bb_oo("k,j") * u1_bb_vo("b,k") * u1_bb_vo("a,i");

        // sharedOps["aaaa_vvoo_353"] += 1.000000 dp_aa_oo("k,j") u1_aa_vo("a,i") u1_aa_vo("b,k")
        // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o2v2: 2, o3v1: 1,
        sharedOps["aaaa_vvoo_353"]("a,b,j,i") = dp_aa_oo("k,j") * u1_aa_vo("a,i") * u1_aa_vo("b,k");

        // sharedOps["abab_vvoo_354"] += 1.000000 dp_bb_oo("k,j") u1_bb_vo("b,k") u1_aa_vo("a,i")
        // flops: o2v2: 2, o2v1: 1 | mem: o2v2: 2, o1v1: 1,
        sharedOps["abab_vvoo_354"]("a,b,i,j") = dp_bb_oo("k,j") * u1_bb_vo("b,k") * u1_aa_vo("a,i");

        // sharedOps["abab_vvoo_355"] += 1.000000 dp_aa_oo("k,j") u1_aa_vo("b,k") u1_bb_vo("a,i")
        // flops: o2v2: 2, o2v1: 1 | mem: o2v2: 2, o1v1: 1,
        sharedOps["abab_vvoo_355"]("b,a,j,i") = dp_aa_oo("k,j") * u1_aa_vo("b,k") * u1_bb_vo("a,i");

        // sharedOps["aaaa_vvoo_356"] += 1.000000 dp_aa_vo("a,i") u1_aa_vo("b,j")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_356"]("a,b,i,j") = dp_aa_vo("a,i") * u1_aa_vo("b,j");

        // sharedOps["abab_vvoo_357"] += 1.000000 dp_bb_vo("a,i") u1_aa_vo("b,j")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["abab_vvoo_357"]("b,a,j,i") = dp_bb_vo("a,i") * u1_aa_vo("b,j");

        // sharedOps["bbbb_vvoo_358"] += 1.000000 u1_bb_vo("b,j") dp_bb_vo("a,i")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_358"]("b,a,j,i") = u1_bb_vo("b,j") * dp_bb_vo("a,i");

        // sharedOps["abab_vvoo_359"] += 1.000000 dp_aa_vo("a,i") u1_bb_vo("b,j")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["abab_vvoo_359"]("a,b,i,j") = dp_aa_vo("a,i") * u1_bb_vo("b,j");

        // sharedOps["aa_vo_360"] += 1.000000 dp_aa_vv("a,b") u1_aa_vo("b,i")
        // flops: o1v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_360"]("a,i") = dp_aa_vv("a,b") * u1_aa_vo("b,i");

        // sharedOps["bb_vo_361"] += 1.000000 dp_bb_vv("b,c") u1_bb_vo("c,j")
        // flops: o1v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_361"]("b,j") = dp_bb_vv("b,c") * u1_bb_vo("c,j");

        // sharedOps["aa_vo_362"] += 1.000000 F_blks_["aa_vv"]("a,b") u1_aa_vo("b,i")
        // flops: o1v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_362"]("a,i") = F_blks_["aa_vv"]("a,b") * u1_aa_vo("b,i");

        // sharedOps["bb_vo_363"] += 1.000000 F_blks_["bb_vv"]("a,b") u1_bb_vo("b,i")
        // flops: o1v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_363"]("a,i") = F_blks_["bb_vv"]("a,b") * u1_bb_vo("b,i");

        // sharedOps["aa_vo_364"] += 1.000000 dp_aa_oo("j,i") u1_aa_vo("a,j")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_364"]("a,i") = dp_aa_oo("j,i") * u1_aa_vo("a,j");

        // sharedOps["bb_vo_365"] += 1.000000 u1_bb_vo("a,j") dp_bb_oo("j,i")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_365"]("a,i") = u1_bb_vo("a,j") * dp_bb_oo("j,i");

        // sharedOps["bb_oo_366"] += 1.000000 u1_bb_vo("c,j") F_blks_["bb_ov"]("k,c")
        // flops: o2v1: 1, o2v0: 1 | mem: o2v0: 2,
        sharedOps["bb_oo_366"]("j,k") = u1_bb_vo("c,j") * F_blks_["bb_ov"]("k,c");

        // sharedOps["aa_oo_367"] += 1.000000 u1_aa_vo("c,j") F_blks_["aa_ov"]("k,c")
        // flops: o2v1: 1, o2v0: 1 | mem: o2v0: 2,
        sharedOps["aa_oo_367"]("j,k") = u1_aa_vo("c,j") * F_blks_["aa_ov"]("k,c");

        // sharedOps["aa_vo_368"] += 1.000000 F_blks_["aa_oo"]("j,i") u1_aa_vo("a,j")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_368"]("a,i") = F_blks_["aa_oo"]("j,i") * u1_aa_vo("a,j");

        // sharedOps["bb_vo_369"] += 1.000000 F_blks_["bb_oo"]("j,i") u1_bb_vo("a,j")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_369"]("a,i") = F_blks_["bb_oo"]("j,i") * u1_bb_vo("a,j");
    }

    {

        // sharedOps["abab_vvoo_52"] += 1.000000 tempOp1_abab_vvoo("a,c,i,k") t2_bbbb_vvoo("c,b,j,k") 1 V_blks_["abab_oovv"]("l,k,d,c") t2_aaaa_vvoo("d,a,i,l")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_52"]("a,b,i,j") = sharedOps["abab_vvoo_0"]("a,c,i,k") * t2_bbbb_vvoo("c,b,j,k");

        // sharedOps["aaaa_vvoo_57"] += 1.000000 t2_abab_vvoo("b,c,j,k") tempOp1_abab_vvoo("a,c,i,k") 1 V_blks_["abab_oovv"]("l,k,d,c") t2_aaaa_vvoo("d,a,i,l")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_57"]("b,a,j,i") = t2_abab_vvoo("b,c,j,k") * sharedOps["abab_vvoo_0"]("a,c,i,k");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_101"] += 1.000000 u2_bbbb_vvoo("c,b,j,k") tempOp1_abab_vvoo("a,c,i,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_101"]("a,b,i,j") = u2_bbbb_vvoo("c,b,j,k") * sharedOps["abab_vvoo_0"]("a,c,i,k");

        // sharedOps["aaaa_vvoo_102"] += 1.000000 u2_abab_vvoo("b,c,j,k") tempOp1_abab_vvoo("a,c,i,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_102"]("b,a,j,i") = u2_abab_vvoo("b,c,j,k") * sharedOps["abab_vvoo_0"]("a,c,i,k");
    }

    if (include_u1_) {

        // sharedOps["aabb_vooo_303"] += 1.000000 u1_bb_vo("b,j") tempOp1_abab_vvoo("a,b,i,m")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_303"]("a,i,j,m") = u1_bb_vo("b,j") * sharedOps["abab_vvoo_0"]("a,b,i,m");
    }

    {

        // sharedOps["abab_vvoo_371"] += 1.000000 t2_abab_vvoo("a,d,i,l") tempOp2_bbbb_vvoo("d,b,l,j") 1 V_blks_["abab_oovv"]("l,k,d,c") t2_abab_vvoo("d,a,l,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_371"]("a,b,i,j") = t2_abab_vvoo("a,d,i,l") * sharedOps["bbbb_vvoo_1"]("d,b,l,j");

        // sharedOps["bbbb_vvoo_372"] += 1.000000 tempOp2_bbbb_vvoo("c,a,k,i") t2_bbbb_vvoo("c,b,j,k") 1 V_blks_["abab_oovv"]("l,k,d,c") t2_abab_vvoo("d,a,l,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_372"]("a,b,i,j") = sharedOps["bbbb_vvoo_1"]("c,a,k,i") * t2_bbbb_vvoo("c,b,j,k");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_374"] += 1.000000 tempOp2_bbbb_vvoo("c,a,k,i") u2_abab_vvoo("b,c,j,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_374"]("b,a,j,i") = sharedOps["bbbb_vvoo_1"]("c,a,k,i") * u2_abab_vvoo("b,c,j,k");

        // sharedOps["bbbb_vvoo_383"] += 1.000000 u2_bbbb_vvoo("c,b,j,k") tempOp2_bbbb_vvoo("c,a,k,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_383"]("b,a,j,i") = u2_bbbb_vvoo("c,b,j,k") * sharedOps["bbbb_vvoo_1"]("c,a,k,i");
    }

    if (include_u1_) {

        // sharedOps["bbbb_vooo_471"] += 1.000000 u1_bb_vo("b,j") tempOp2_bbbb_vvoo("b,a,m,i")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["bbbb_vooo_471"]("a,j,m,i") = u1_bb_vo("b,j") * sharedOps["bbbb_vvoo_1"]("b,a,m,i");
    }

    {

        // sharedOps["abab_vvoo_373"] += 1.000000 tempOp3_bbbb_vvoo("a,c,i,k") t2_abab_vvoo("b,c,j,k") 1 t2_bbbb_vvoo("d,a,i,l") V_blks_["bbbb_oovv"]("l,k,c,d")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_373"]("b,a,j,i") = sharedOps["bbbb_vvoo_2"]("a,c,i,k") * t2_abab_vvoo("b,c,j,k");

        // sharedOps["bbbb_vvoo_378"] += 1.000000 tempOp3_bbbb_vvoo("a,c,i,k") t2_bbbb_vvoo("c,b,j,k") 1 t2_bbbb_vvoo("d,a,i,l") V_blks_["bbbb_oovv"]("l,k,c,d")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_378"]("a,b,i,j") = sharedOps["bbbb_vvoo_2"]("a,c,i,k") * t2_bbbb_vvoo("c,b,j,k");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_380"] += 1.000000 u2_abab_vvoo("b,c,j,k") tempOp3_bbbb_vvoo("a,c,i,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_380"]("b,a,j,i") = u2_abab_vvoo("b,c,j,k") * sharedOps["bbbb_vvoo_2"]("a,c,i,k");

        // sharedOps["bbbb_vvoo_388"] += 1.000000 tempOp3_bbbb_vvoo("a,c,i,k") u2_bbbb_vvoo("c,b,j,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_388"]("a,b,i,j") = sharedOps["bbbb_vvoo_2"]("a,c,i,k") * u2_bbbb_vvoo("c,b,j,k");
    }

    if (include_u1_) {

        // sharedOps["bbbb_vooo_450"] += 1.000000 tempOp3_bbbb_vvoo("a,b,i,m") u1_bb_vo("b,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["bbbb_vooo_450"]("a,i,m,j") = sharedOps["bbbb_vvoo_2"]("a,b,i,m") * u1_bb_vo("b,j");
    }

    {

        // sharedOps["abab_vvoo_375"] += 1.000000 tempOp4_abab_vvoo("a,c,i,k") t2_bbbb_vvoo("c,b,j,k") 1 t2_abab_vvoo("a,d,i,l") V_blks_["bbbb_oovv"]("l,k,c,d")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_375"]("a,b,i,j") = sharedOps["abab_vvoo_3"]("a,c,i,k") * t2_bbbb_vvoo("c,b,j,k");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_382"] += 1.000000 u2_bbbb_vvoo("c,b,j,k") tempOp4_abab_vvoo("a,c,i,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_382"]("a,b,i,j") = u2_bbbb_vvoo("c,b,j,k") * sharedOps["abab_vvoo_3"]("a,c,i,k");

        // sharedOps["aaaa_vvoo_385"] += 1.000000 u2_abab_vvoo("b,c,j,k") tempOp4_abab_vvoo("a,c,i,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_385"]("b,a,j,i") = u2_abab_vvoo("b,c,j,k") * sharedOps["abab_vvoo_3"]("a,c,i,k");
    }

    {

        // sharedOps["aaaa_vvoo_391"] += 1.000000 tempOp4_abab_vvoo("a,c,i,k") t2_abab_vvoo("b,c,j,k") 1 t2_abab_vvoo("a,d,i,l") V_blks_["bbbb_oovv"]("l,k,c,d")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_391"]("a,b,i,j") = sharedOps["abab_vvoo_3"]("a,c,i,k") * t2_abab_vvoo("b,c,j,k");
    }

    if (include_u1_) {

        // sharedOps["aabb_vooo_461"] += 1.000000 u1_bb_vo("b,j") tempOp4_abab_vvoo("a,b,i,m")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_461"]("a,i,j,m") = u1_bb_vo("b,j") * sharedOps["abab_vvoo_3"]("a,b,i,m");
    }

    if (include_u2_) {

        // sharedOps["bbbb_vvoo_379"] += 1.000000 u2_abab_vvoo("c,b,k,j") tempOp5_abab_vvoo("c,a,k,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_379"]("b,a,j,i") = u2_abab_vvoo("c,b,k,j") * sharedOps["abab_vvoo_4"]("c,a,k,i");
    }

    {

        // sharedOps["bbbb_vvoo_381"] += 1.000000 tempOp5_abab_vvoo("c,a,k,i") t2_abab_vvoo("c,b,k,j") 1 V_blks_["aaaa_oovv"]("l,k,c,d") t2_abab_vvoo("d,a,l,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_381"]("a,b,i,j") = sharedOps["abab_vvoo_4"]("c,a,k,i") * t2_abab_vvoo("c,b,k,j");

        // sharedOps["abab_vvoo_394"] += 1.000000 t2_aaaa_vvoo("c,b,j,k") tempOp5_abab_vvoo("c,a,k,i") 1 V_blks_["aaaa_oovv"]("l,k,c,d") t2_abab_vvoo("d,a,l,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_394"]("b,a,j,i") = t2_aaaa_vvoo("c,b,j,k") * sharedOps["abab_vvoo_4"]("c,a,k,i");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_395"] += 1.000000 tempOp5_abab_vvoo("c,a,k,i") u2_aaaa_vvoo("c,b,j,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_395"]("b,a,j,i") = sharedOps["abab_vvoo_4"]("c,a,k,i") * u2_aaaa_vvoo("c,b,j,k");
    }

    if (include_u1_) {

        // sharedOps["baab_vooo_468"] += 1.000000 tempOp5_abab_vvoo("b,a,m,i") u1_aa_vo("b,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_468"]("a,m,j,i") = sharedOps["abab_vvoo_4"]("b,a,m,i") * u1_aa_vo("b,j");
    }

    if (include_u2_) {

        // sharedOps["aaaa_vvoo_387"] += 1.000000 tempOp6_aaaa_vvoo("c,a,k,i") u2_aaaa_vvoo("c,b,j,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_387"]("a,b,i,j") = sharedOps["aaaa_vvoo_5"]("c,a,k,i") * u2_aaaa_vvoo("c,b,j,k");
    }

    {

        // sharedOps["aaaa_vvoo_392"] += 1.000000 t2_aaaa_vvoo("c,b,j,k") tempOp6_aaaa_vvoo("c,a,k,i") 1 V_blks_["aaaa_oovv"]("l,k,c,d") t2_aaaa_vvoo("d,a,i,l")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_392"]("b,a,j,i") = t2_aaaa_vvoo("c,b,j,k") * sharedOps["aaaa_vvoo_5"]("c,a,k,i");

        // sharedOps["abab_vvoo_393"] += 1.000000 tempOp6_aaaa_vvoo("c,a,k,i") t2_abab_vvoo("c,b,k,j") 1 V_blks_["aaaa_oovv"]("l,k,c,d") t2_aaaa_vvoo("d,a,i,l")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_393"]("a,b,i,j") = sharedOps["aaaa_vvoo_5"]("c,a,k,i") * t2_abab_vvoo("c,b,k,j");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_396"] += 1.000000 tempOp6_aaaa_vvoo("c,a,k,i") u2_abab_vvoo("c,b,k,j")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_396"]("a,b,i,j") = sharedOps["aaaa_vvoo_5"]("c,a,k,i") * u2_abab_vvoo("c,b,k,j");
    }

    if (include_u1_) {

        // sharedOps["aaaa_vooo_465"] += 1.000000 tempOp6_aaaa_vvoo("b,a,m,i") u1_aa_vo("b,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aaaa_vooo_465"]("a,m,i,j") = sharedOps["aaaa_vvoo_5"]("b,a,m,i") * u1_aa_vo("b,j");
    }

    if (include_u2_) {

        // sharedOps["aaaa_vvoo_376"] += 1.000000 u2_aaaa_vvoo("c,b,j,k") tempOp7_aaaa_vvoo("a,c,i,k")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_376"]("b,a,j,i") = u2_aaaa_vvoo("c,b,j,k") * sharedOps["aaaa_vvoo_6"]("a,c,i,k");

        // sharedOps["abab_vvoo_384"] += 1.000000 tempOp7_aaaa_vvoo("a,c,i,k") u2_abab_vvoo("c,b,k,j")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_384"]("a,b,i,j") = sharedOps["aaaa_vvoo_6"]("a,c,i,k") * u2_abab_vvoo("c,b,k,j");
    }

    if (include_u1_) {

        // sharedOps["aaaa_vooo_453"] += 1.000000 tempOp7_aaaa_vvoo("a,b,i,m") u1_aa_vo("b,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aaaa_vooo_453"]("a,i,m,j") = sharedOps["aaaa_vvoo_6"]("a,b,i,m") * u1_aa_vo("b,j");
    }

    if (include_u2_) {

        // sharedOps["bbbb_vvoo_377"] += 1.000000 u2_abab_vvoo("c,b,k,j") tempOp8_abab_vvoo("c,a,k,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_377"]("b,a,j,i") = u2_abab_vvoo("c,b,k,j") * sharedOps["abab_vvoo_7"]("c,a,k,i");

        // sharedOps["abab_vvoo_386"] += 1.000000 u2_aaaa_vvoo("c,b,j,k") tempOp8_abab_vvoo("c,a,k,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_386"]("b,a,j,i") = u2_aaaa_vvoo("c,b,j,k") * sharedOps["abab_vvoo_7"]("c,a,k,i");
    }

    if (include_u1_) {

        // sharedOps["baab_vooo_467"] += 1.000000 u1_aa_vo("b,j") tempOp8_abab_vvoo("b,a,m,i")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_467"]("a,j,m,i") = u1_aa_vo("b,j") * sharedOps["abab_vvoo_7"]("b,a,m,i");
    }

    {

        // sharedOps["abab_vvoo_370"] += 1.000000 t2_abab_vvoo("b,c,k,j") tempOp9_bbaa_vvoo("c,a,k,i") 1 V_blks_["abab_oovv"]("k,l,d,c") t2_abab_vvoo("d,a,i,l")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_370"]("b,a,i,j") = t2_abab_vvoo("b,c,k,j") * sharedOps["bbaa_vvoo_8"]("c,a,k,i");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_390"] += 1.000000 u2_abab_vvoo("b,c,k,j") tempOp9_bbaa_vvoo("c,a,k,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_390"]("b,a,i,j") = u2_abab_vvoo("b,c,k,j") * sharedOps["bbaa_vvoo_8"]("c,a,k,i");
    }

    if (include_u1_) {

        // sharedOps["baab_vooo_466"] += 1.000000 tempOp9_bbaa_vvoo("b,a,m,i") u1_bb_vo("b,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_466"]("a,m,i,j") = sharedOps["bbaa_vvoo_8"]("b,a,m,i") * u1_bb_vo("b,j");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_389"] += 1.000000 u2_abab_vvoo("c,b,j,k") tempOp10_aabb_vvoo("c,a,k,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_389"]("a,b,j,i") = u2_abab_vvoo("c,b,j,k") * sharedOps["aabb_vvoo_9"]("c,a,k,i");
    }

    if (include_u1_) {

        // sharedOps["aabb_vooo_473"] += 1.000000 tempOp10_aabb_vvoo("b,a,m,i") u1_aa_vo("b,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_473"]("a,j,m,i") = sharedOps["aabb_vvoo_9"]("b,a,m,i") * u1_aa_vo("b,j");
    }

    {

        // sharedOps["bbbb_vooo_472"] += 1.000000 tempOp37_bb_vo("c,m") t2_bbbb_vvoo("c,a,i,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["bbbb_vooo_472"]("a,m,i,j") = sharedOps["bb_vo_36"]("c,m") * t2_bbbb_vvoo("c,a,i,j");

        // sharedOps["aabb_vooo_474"] += 1.000000 tempOp37_bb_vo("c,m") t2_abab_vvoo("a,c,i,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_474"]("a,i,m,j") = sharedOps["bb_vo_36"]("c,m") * t2_abab_vvoo("a,c,i,j");

        // sharedOps["bb_vo_480"] += 1.000000 t2_bbbb_vvoo("c,a,i,k") tempOp37_bb_vo("c,k")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_480"]("a,i") = t2_bbbb_vvoo("c,a,i,k") * sharedOps["bb_vo_36"]("c,k");

        // sharedOps["aa_vo_482"] += 1.000000 tempOp37_bb_vo("c,k") t2_abab_vvoo("a,c,i,k")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_482"]("a,i") = sharedOps["bb_vo_36"]("c,k") * t2_abab_vvoo("a,c,i,k");

        // sharedOps["baab_vooo_469"] += 1.000000 tempOp38_aa_vo("c,m") t2_abab_vvoo("c,a,i,j")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_469"]("a,m,i,j") = sharedOps["aa_vo_37"]("c,m") * t2_abab_vvoo("c,a,i,j");

        // sharedOps["aaaa_vooo_470"] += 1.000000 t2_aaaa_vvoo("c,a,i,j") tempOp38_aa_vo("c,m")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aaaa_vooo_470"]("a,i,j,m") = t2_aaaa_vvoo("c,a,i,j") * sharedOps["aa_vo_37"]("c,m");

        // sharedOps["aa_vo_479"] += 1.000000 tempOp38_aa_vo("c,k") t2_aaaa_vvoo("c,a,i,k")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_479"]("a,i") = sharedOps["aa_vo_37"]("c,k") * t2_aaaa_vvoo("c,a,i,k");

        // sharedOps["bb_vo_481"] += 1.000000 tempOp38_aa_vo("c,k") t2_abab_vvoo("c,a,k,i")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_481"]("a,i") = sharedOps["aa_vo_37"]("c,k") * t2_abab_vvoo("c,a,k,i");
    }

    if (include_u1_) {

        // sharedOps["abab_vvoo_427"] += 1.000000 tempOp41_baab_vooo("b,k,i,j") u1_aa_vo("a,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_427"]("a,b,i,j") = sharedOps["baab_vooo_40"]("b,k,i,j") * u1_aa_vo("a,k");

        // sharedOps["abab_vvoo_423"] += 1.000000 tempOp42_aabb_vooo("b,j,k,i") u1_bb_vo("a,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_423"]("b,a,j,i") = sharedOps["aabb_vooo_41"]("b,j,k,i") * u1_bb_vo("a,k");

        // sharedOps["bbbb_vvoo_447"] += 1.000000 tempOp43_bbbb_vooo("b,i,j,k") u1_bb_vo("a,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_447"]("b,a,i,j") = sharedOps["bbbb_vooo_42"]("b,i,j,k") * u1_bb_vo("a,k");

        // sharedOps["aaaa_vvoo_445"] += 1.000000 tempOp44_aaaa_vooo("b,k,i,j") u1_aa_vo("a,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_445"]("b,a,i,j") = sharedOps["aaaa_vooo_43"]("b,k,i,j") * u1_aa_vo("a,k");

        // sharedOps["aabb_vvvo_404"] += 1.000000 tempOp59_aabb_vvvv("e,b,c,a") u1_bb_vo("c,i")
        // flops: o1v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["aabb_vvvo_404"]("e,b,a,i") = sharedOps["aabb_vvvv_58"]("e,b,c,a") * u1_bb_vo("c,i");

        // sharedOps["abba_vvvo_405"] += 1.000000 u1_aa_vo("c,i") tempOp59_aabb_vvvv("c,b,e,a")
        // flops: o1v4: 1, o1v3: 1 | mem: o1v3: 2,
        sharedOps["abba_vvvo_405"]("b,e,a,i") = u1_aa_vo("c,i") * sharedOps["aabb_vvvv_58"]("c,b,e,a");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_397"] += 1.000000 u2_abab_vvoo("a,b,l,k") tempOp104_aabb_oooo("l,j,k,i")
        // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_397"]("a,b,j,i") = u2_abab_vvoo("a,b,l,k") * sharedOps["aabb_oooo_103"]("l,j,k,i");
    }

    {

        // sharedOps["abab_vvoo_398"] += 1.000000 tempOp104_aabb_oooo("l,j,k,i") t2_abab_vvoo("a,b,l,k") 1 V_blks_["abab_oovv"]("k,l,c,d") t2_abab_vvoo("c,d,j,i")
        // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_398"]("a,b,j,i") = sharedOps["aabb_oooo_103"]("l,j,k,i") * t2_abab_vvoo("a,b,l,k");
    }

    if (include_u1_) {

        // sharedOps["baab_vooo_476"] += 1.000000 tempOp104_aabb_oooo("m,i,k,j") u1_bb_vo("a,k")
        // flops: o4v1: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["baab_vooo_476"]("a,m,i,j") = sharedOps["aabb_oooo_103"]("m,i,k,j") * u1_bb_vo("a,k");

        // sharedOps["aabb_vooo_478"] += 1.000000 tempOp104_aabb_oooo("k,i,m,j") u1_aa_vo("a,k")
        // flops: o4v1: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aabb_vooo_478"]("a,i,m,j") = sharedOps["aabb_oooo_103"]("k,i,m,j") * u1_aa_vo("a,k");
    }

    if (include_u2_) {

        // sharedOps["bbbb_vvoo_399"] += 1.000000 tempOp108_bbbb_oooo("i,j,l,k") u2_bbbb_vvoo("b,a,l,k")
        // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_399"]("b,a,i,j") = sharedOps["bbbb_oooo_107"]("i,j,l,k") * u2_bbbb_vvoo("b,a,l,k");
    }

    {

        // sharedOps["bbbb_vvoo_400"] += 1.000000 tempOp108_bbbb_oooo("i,j,l,k") t2_bbbb_vvoo("b,a,l,k") 1 t2_bbbb_vvoo("c,d,i,j") V_blks_["bbbb_oovv"]("l,k,c,d")
        // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_400"]("b,a,i,j") = sharedOps["bbbb_oooo_107"]("i,j,l,k") * t2_bbbb_vvoo("b,a,l,k");
    }

    if (include_u1_) {

        // sharedOps["bbbb_vooo_477"] += 1.000000 u1_bb_vo("a,k") tempOp108_bbbb_oooo("i,j,k,m")
        // flops: o4v1: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["bbbb_vooo_477"]("a,i,j,m") = u1_bb_vo("a,k") * sharedOps["bbbb_oooo_107"]("i,j,k,m");
    }

    {

        // sharedOps["aaaa_vvoo_401"] += 1.000000 t2_aaaa_vvoo("b,a,l,k") tempOp109_aaaa_oooo("i,j,l,k") 1 t2_aaaa_vvoo("c,d,i,j") V_blks_["aaaa_oovv"]("l,k,c,d")
        // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_401"]("b,a,i,j") = t2_aaaa_vvoo("b,a,l,k") * sharedOps["aaaa_oooo_108"]("i,j,l,k");
    }

    if (include_u2_) {

        // sharedOps["aaaa_vvoo_403"] += 1.000000 tempOp109_aaaa_oooo("i,j,l,k") u2_aaaa_vvoo("b,a,l,k")
        // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_403"]("b,a,i,j") = sharedOps["aaaa_oooo_108"]("i,j,l,k") * u2_aaaa_vvoo("b,a,l,k");
    }

    if (include_u1_) {

        // sharedOps["aaaa_vooo_475"] += 1.000000 tempOp109_aaaa_oooo("i,j,k,m") u1_aa_vo("a,k")
        // flops: o4v1: 1, o3v1: 1 | mem: o3v1: 2,
        sharedOps["aaaa_vooo_475"]("a,i,j,m") = sharedOps["aaaa_oooo_108"]("i,j,k,m") * u1_aa_vo("a,k");

        // sharedOps["abab_vvoo_448"] += 1.000000 u1_aa_vo("b,k") tempOp112_baab_vooo("a,k,j,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_448"]("b,a,j,i") = u1_aa_vo("b,k") * sharedOps["baab_vooo_111"]("a,k,j,i");

        // sharedOps["bbbb_vvoo_456"] += 1.000000 u1_bb_vo("b,k") tempOp113_bbbb_vooo("a,i,k,j")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_456"]("b,a,i,j") = u1_bb_vo("b,k") * sharedOps["bbbb_vooo_112"]("a,i,k,j");

        // sharedOps["abab_vvoo_463"] += 1.000000 u1_bb_vo("b,k") tempOp114_aabb_vooo("a,j,k,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_463"]("a,b,j,i") = u1_bb_vo("b,k") * sharedOps["aabb_vooo_113"]("a,j,k,i");

        // sharedOps["abab_vvoo_464"] += 1.000000 tempOp115_aabb_vooo("a,i,k,j") u1_bb_vo("b,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_464"]("a,b,i,j") = sharedOps["aabb_vooo_114"]("a,i,k,j") * u1_bb_vo("b,k");

        // sharedOps["aaaa_vvoo_459"] += 1.000000 u1_aa_vo("b,k") tempOp116_aaaa_vooo("a,k,j,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_459"]("b,a,j,i") = u1_aa_vo("b,k") * sharedOps["aaaa_vooo_115"]("a,k,j,i");

        // sharedOps["bbbb_vvoo_452"] += 1.000000 u1_bb_vo("b,k") tempOp117_bbbb_vooo("a,i,k,j")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_452"]("b,a,i,j") = u1_bb_vo("b,k") * sharedOps["bbbb_vooo_116"]("a,i,k,j");

        // sharedOps["abab_vvoo_458"] += 1.000000 tempOp118_baab_vooo("a,k,j,i") u1_aa_vo("b,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_458"]("b,a,j,i") = sharedOps["baab_vooo_117"]("a,k,j,i") * u1_aa_vo("b,k");

        // sharedOps["abab_vvoo_454"] += 1.000000 u1_aa_vo("b,k") tempOp119_baab_vooo("a,k,i,j")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_454"]("b,a,i,j") = u1_aa_vo("b,k") * sharedOps["baab_vooo_118"]("a,k,i,j");

        // sharedOps["aaaa_vvoo_462"] += 1.000000 u1_aa_vo("b,k") tempOp120_aaaa_vooo("a,k,j,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_462"]("b,a,j,i") = u1_aa_vo("b,k") * sharedOps["aaaa_vooo_119"]("a,k,j,i");

        // sharedOps["abab_vvoo_457"] += 1.000000 u1_bb_vo("b,k") tempOp121_aabb_vooo("a,i,k,j")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_457"]("a,b,i,j") = u1_bb_vo("b,k") * sharedOps["aabb_vooo_120"]("a,i,k,j");
    }

    {

        // sharedOps["aaaa_vvoo_402"] += 1.000000 t2_aaaa_vvoo("b,a,l,k") tempOp122_aaaa_oooo("l,k,i,j")
        // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_402"]("b,a,i,j") = t2_aaaa_vvoo("b,a,l,k") * sharedOps["aaaa_oooo_121"]("l,k,i,j");

        // sharedOps["bb_vo_435"] += 1.000000 t2_abab_vvoo("c,a,k,j") tempOp146_aabb_vooo("c,k,i,j")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_435"]("a,i") = t2_abab_vvoo("c,a,k,j") * sharedOps["aabb_vooo_145"]("c,k,i,j");

        // sharedOps["aa_vo_436"] += 1.000000 tempOp147_baab_vooo("c,i,k,j") t2_abab_vvoo("a,c,k,j")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_436"]("a,i") = sharedOps["baab_vooo_146"]("c,i,k,j") * t2_abab_vvoo("a,c,k,j");

        // sharedOps["abab_vvoo_406"] += 1.000000 t2_abab_vvoo("a,d,j,i") tempOp150_bb_vv("d,b") 1 V_blks_["abab_oovv"]("k,l,d,c") t2_abab_vvoo("d,a,k,l")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_406"]("a,b,j,i") = t2_abab_vvoo("a,d,j,i") * sharedOps["bb_vv_149"]("d,b");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_409"] += 1.000000 tempOp150_bb_vv("c,a") u2_abab_vvoo("b,c,j,i")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_409"]("b,a,j,i") = sharedOps["bb_vv_149"]("c,a") * u2_abab_vvoo("b,c,j,i");

        // sharedOps["bbbb_vvoo_415"] += 1.000000 u2_bbbb_vvoo("c,b,i,j") tempOp150_bb_vv("c,a")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_415"]("b,a,i,j") = u2_bbbb_vvoo("c,b,i,j") * sharedOps["bb_vv_149"]("c,a");
    }

    {

        // sharedOps["bbbb_vvoo_417"] += 1.000000 tempOp150_bb_vv("d,b") t2_bbbb_vvoo("d,a,i,j") 1 V_blks_["abab_oovv"]("k,l,d,c") t2_abab_vvoo("d,a,k,l")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_417"]("b,a,i,j") = sharedOps["bb_vv_149"]("d,b") * t2_bbbb_vvoo("d,a,i,j");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_407"] += 1.000000 u2_abab_vvoo("c,b,j,i") tempOp151_aa_vv("a,c")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_407"]("a,b,j,i") = u2_abab_vvoo("c,b,j,i") * sharedOps["aa_vv_150"]("a,c");
    }

    {

        // sharedOps["abab_vvoo_410"] += 1.000000 tempOp151_aa_vv("b,d") t2_abab_vvoo("d,a,i,j") 1 t2_abab_vvoo("a,d,l,k") V_blks_["abab_oovv"]("l,k,c,d")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_410"]("b,a,i,j") = sharedOps["aa_vv_150"]("b,d") * t2_abab_vvoo("d,a,i,j");
    }

    if (include_u2_) {

        // sharedOps["aaaa_vvoo_412"] += 1.000000 u2_aaaa_vvoo("c,b,i,j") tempOp151_aa_vv("a,c")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_412"]("b,a,i,j") = u2_aaaa_vvoo("c,b,i,j") * sharedOps["aa_vv_150"]("a,c");
    }

    {

        // sharedOps["aaaa_vvoo_418"] += 1.000000 tempOp151_aa_vv("b,d") t2_aaaa_vvoo("d,a,i,j") 1 t2_abab_vvoo("a,d,l,k") V_blks_["abab_oovv"]("l,k,c,d")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_418"]("b,a,i,j") = sharedOps["aa_vv_150"]("b,d") * t2_aaaa_vvoo("d,a,i,j");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_414"] += 1.000000 u2_abab_vvoo("c,b,i,j") tempOp152_aa_vv("a,c")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_414"]("a,b,i,j") = u2_abab_vvoo("c,b,i,j") * sharedOps["aa_vv_151"]("a,c");

        // sharedOps["aaaa_vvoo_421"] += 1.000000 tempOp152_aa_vv("a,c") u2_aaaa_vvoo("c,b,i,j")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_421"]("a,b,i,j") = sharedOps["aa_vv_151"]("a,c") * u2_aaaa_vvoo("c,b,i,j");

        // sharedOps["abab_vvoo_416"] += 1.000000 tempOp153_bb_vv("c,a") u2_abab_vvoo("b,c,j,i")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_416"]("b,a,j,i") = sharedOps["bb_vv_152"]("c,a") * u2_abab_vvoo("b,c,j,i");

        // sharedOps["bbbb_vvoo_420"] += 1.000000 u2_bbbb_vvoo("c,b,i,j") tempOp153_bb_vv("c,a")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_420"]("b,a,i,j") = u2_bbbb_vvoo("c,b,i,j") * sharedOps["bb_vv_152"]("c,a");
    }

    {

        // sharedOps["abab_vvoo_408"] += 1.000000 t2_abab_vvoo("d,a,i,j") tempOp154_aa_vv("b,d")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_408"]("b,a,i,j") = t2_abab_vvoo("d,a,i,j") * sharedOps["aa_vv_153"]("b,d");

        // sharedOps["aaaa_vvoo_413"] += 1.000000 t2_aaaa_vvoo("d,a,i,j") tempOp154_aa_vv("b,d")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_413"]("a,b,i,j") = t2_aaaa_vvoo("d,a,i,j") * sharedOps["aa_vv_153"]("b,d");

        // sharedOps["abab_vvoo_411"] += 1.000000 t2_abab_vvoo("a,d,j,i") tempOp155_bb_vv("b,d")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_411"]("a,b,j,i") = t2_abab_vvoo("a,d,j,i") * sharedOps["bb_vv_154"]("b,d");

        // sharedOps["bbbb_vvoo_419"] += 1.000000 tempOp155_bb_vv("b,d") t2_bbbb_vvoo("d,a,i,j")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_419"]("b,a,i,j") = sharedOps["bb_vv_154"]("b,d") * t2_bbbb_vvoo("d,a,i,j");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_422"] += 1.000000 u2_abab_vvoo("b,a,k,j") tempOp208_aa_oo("k,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_422"]("b,a,i,j") = u2_abab_vvoo("b,a,k,j") * sharedOps["aa_oo_207"]("k,i");
    }

    {

        // sharedOps["abab_vvoo_426"] += 1.000000 tempOp208_aa_oo("l,j") t2_abab_vvoo("b,a,l,i") 1 V_blks_["abab_oovv"]("j,k,b,c") t2_abab_vvoo("b,c,i,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_426"]("b,a,j,i") = sharedOps["aa_oo_207"]("l,j") * t2_abab_vvoo("b,a,l,i");
    }

    if (include_u2_) {

        // sharedOps["aaaa_vvoo_437"] += 1.000000 tempOp208_aa_oo("k,i") u2_aaaa_vvoo("b,a,j,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_437"]("b,a,i,j") = sharedOps["aa_oo_207"]("k,i") * u2_aaaa_vvoo("b,a,j,k");
    }

    {

        // sharedOps["aaaa_vvoo_439"] += 1.000000 tempOp208_aa_oo("l,j") t2_aaaa_vvoo("b,a,i,l") 1 V_blks_["abab_oovv"]("j,k,b,c") t2_abab_vvoo("b,c,i,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_439"]("b,a,j,i") = sharedOps["aa_oo_207"]("l,j") * t2_aaaa_vvoo("b,a,i,l");
    }

    if (include_u1_) {

        // sharedOps["aa_vo_500"] += 1.000000 tempOp208_aa_oo("j,i") u1_aa_vo("a,j")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_500"]("a,i") = sharedOps["aa_oo_207"]("j,i") * u1_aa_vo("a,j");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_424"] += 1.000000 u2_abab_vvoo("a,b,j,k") tempOp209_bb_oo("k,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_424"]("a,b,j,i") = u2_abab_vvoo("a,b,j,k") * sharedOps["bb_oo_208"]("k,i");
    }

    {

        // sharedOps["abab_vvoo_425"] += 1.000000 tempOp209_bb_oo("l,j") t2_abab_vvoo("b,a,i,l") 1 V_blks_["abab_oovv"]("k,j,b,c") t2_abab_vvoo("b,c,k,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_425"]("b,a,i,j") = sharedOps["bb_oo_208"]("l,j") * t2_abab_vvoo("b,a,i,l");
    }

    if (include_u2_) {

        // sharedOps["bbbb_vvoo_440"] += 1.000000 u2_bbbb_vvoo("b,a,j,k") tempOp209_bb_oo("k,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_440"]("b,a,j,i") = u2_bbbb_vvoo("b,a,j,k") * sharedOps["bb_oo_208"]("k,i");
    }

    {

        // sharedOps["bbbb_vvoo_442"] += 1.000000 t2_bbbb_vvoo("b,a,i,l") tempOp209_bb_oo("l,j") 1 V_blks_["abab_oovv"]("k,j,b,c") t2_abab_vvoo("b,c,k,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_442"]("b,a,i,j") = t2_bbbb_vvoo("b,a,i,l") * sharedOps["bb_oo_208"]("l,j");
    }

    if (include_u1_) {

        // sharedOps["bb_vo_499"] += 1.000000 tempOp209_bb_oo("j,i") u1_bb_vo("a,j")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_499"]("a,i") = sharedOps["bb_oo_208"]("j,i") * u1_bb_vo("a,j");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_434"] += 1.000000 tempOp212_aa_oo("i,k") u2_abab_vvoo("b,a,k,j")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_434"]("b,a,i,j") = sharedOps["aa_oo_211"]("i,k") * u2_abab_vvoo("b,a,k,j");

        // sharedOps["aaaa_vvoo_460"] += 1.000000 tempOp212_aa_oo("i,k") u2_aaaa_vvoo("b,a,j,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_460"]("b,a,i,j") = sharedOps["aa_oo_211"]("i,k") * u2_aaaa_vvoo("b,a,j,k");
    }

    if (include_u1_) {

        // sharedOps["aa_vo_503"] += 1.000000 tempOp212_aa_oo("i,j") u1_aa_vo("a,j")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_503"]("a,i") = sharedOps["aa_oo_211"]("i,j") * u1_aa_vo("a,j");
    }

    if (include_u2_) {

        // sharedOps["abab_vvoo_433"] += 1.000000 tempOp213_bb_oo("k,i") u2_abab_vvoo("b,a,j,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_433"]("b,a,j,i") = sharedOps["bb_oo_212"]("k,i") * u2_abab_vvoo("b,a,j,k");

        // sharedOps["bbbb_vvoo_455"] += 1.000000 u2_bbbb_vvoo("b,a,j,k") tempOp213_bb_oo("k,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_455"]("b,a,j,i") = u2_bbbb_vvoo("b,a,j,k") * sharedOps["bb_oo_212"]("k,i");
    }

    if (include_u1_) {

        // sharedOps["bb_vo_504"] += 1.000000 tempOp213_bb_oo("j,i") u1_bb_vo("a,j")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_504"]("a,i") = sharedOps["bb_oo_212"]("j,i") * u1_bb_vo("a,j");
    }

    {

        // sharedOps["abab_vvoo_428"] += 1.000000 tempOp214_aa_oo("l,j") t2_abab_vvoo("a,b,l,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_428"]("a,b,j,i") = sharedOps["aa_oo_213"]("l,j") * t2_abab_vvoo("a,b,l,i");

        // sharedOps["aaaa_vvoo_443"] += 1.000000 tempOp214_aa_oo("l,j") t2_aaaa_vvoo("b,a,i,l")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_443"]("b,a,j,i") = sharedOps["aa_oo_213"]("l,j") * t2_aaaa_vvoo("b,a,i,l");

        // sharedOps["abab_vvoo_429"] += 1.000000 t2_abab_vvoo("b,a,i,l") tempOp215_bb_oo("j,l")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_429"]("b,a,i,j") = t2_abab_vvoo("b,a,i,l") * sharedOps["bb_oo_214"]("j,l");

        // sharedOps["bbbb_vvoo_438"] += 1.000000 tempOp215_bb_oo("j,l") t2_bbbb_vvoo("b,a,i,l")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_438"]("b,a,j,i") = sharedOps["bb_oo_214"]("j,l") * t2_bbbb_vvoo("b,a,i,l");
    }

    if (include_u1_) {

        // sharedOps["aaaa_vvoo_431"] += 1.000000 u1_aa_vo("b,k") tempOp216_aaaa_vooo("a,k,i,j")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_431"]("b,a,i,j") = u1_aa_vo("b,k") * sharedOps["aaaa_vooo_215"]("a,k,i,j");

        // sharedOps["bbbb_vvoo_430"] += 1.000000 tempOp217_bbbb_vooo("a,i,j,k") u1_bb_vo("b,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_430"]("a,b,i,j") = sharedOps["bbbb_vooo_216"]("a,i,j,k") * u1_bb_vo("b,k");

        // sharedOps["abab_vvoo_441"] += 1.000000 tempOp224_baab_vooo("a,i,k,j") u1_aa_vo("b,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_441"]("b,a,i,j") = sharedOps["baab_vooo_223"]("a,i,k,j") * u1_aa_vo("b,k");

        // sharedOps["abab_vvoo_432"] += 1.000000 u1_bb_vo("b,k") tempOp225_aabb_vooo("a,j,i,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["abab_vvoo_432"]("a,b,j,i") = u1_bb_vo("b,k") * sharedOps["aabb_vooo_224"]("a,j,i,k");

        // sharedOps["bbbb_vvoo_446"] += 1.000000 u1_bb_vo("b,k") tempOp233_bbbb_vooo("a,k,i,j")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_446"]("b,a,i,j") = u1_bb_vo("b,k") * sharedOps["bbbb_vooo_232"]("a,k,i,j");

        // sharedOps["aaaa_vvoo_444"] += 1.000000 u1_aa_vo("b,k") tempOp235_aaaa_vooo("a,i,j,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_444"]("b,a,i,j") = u1_aa_vo("b,k") * sharedOps["aaaa_vooo_234"]("a,i,j,k");

        // sharedOps["aaaa_vvoo_451"] += 1.000000 u1_aa_vo("b,k") tempOp241_aaaa_vooo("a,k,i,j")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_451"]("b,a,i,j") = u1_aa_vo("b,k") * sharedOps["aaaa_vooo_240"]("a,k,i,j");

        // sharedOps["bbbb_vvoo_449"] += 1.000000 tempOp244_bbbb_vooo("a,k,i,j") u1_bb_vo("b,k")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_449"]("a,b,i,j") = sharedOps["bbbb_vooo_243"]("a,k,i,j") * u1_bb_vo("b,k");

        // sharedOps["aa_vo_502"] += 1.000000 tempOp313_aa_oo("i,j") u1_aa_vo("a,j")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["aa_vo_502"]("a,i") = sharedOps["aa_oo_312"]("i,j") * u1_aa_vo("a,j");

        // sharedOps["bb_vo_501"] += 1.000000 u1_bb_vo("a,j") tempOp314_bb_oo("j,i")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sharedOps["bb_vo_501"]("a,i") = u1_bb_vo("a,j") * sharedOps["bb_oo_313"]("j,i");

        // sharedOps["aaaa_vvoo_494"] += 1.000000 tempOp319_aa_vo("a,i") u1_aa_vo("b,j")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_494"]("a,b,i,j") = sharedOps["aa_vo_318"]("a,i") * u1_aa_vo("b,j");

        // sharedOps["abab_vvoo_496"] += 1.000000 u1_bb_vo("b,j") tempOp319_aa_vo("a,i")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["abab_vvoo_496"]("a,b,i,j") = u1_bb_vo("b,j") * sharedOps["aa_vo_318"]("a,i");

        // sharedOps["abab_vvoo_495"] += 1.000000 u1_aa_vo("b,j") tempOp320_bb_vo("a,i")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["abab_vvoo_495"]("b,a,j,i") = u1_aa_vo("b,j") * sharedOps["bb_vo_319"]("a,i");

        // sharedOps["bbbb_vvoo_497"] += 1.000000 tempOp320_bb_vo("a,i") u1_bb_vo("b,j")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_497"]("a,b,i,j") = sharedOps["bb_vo_319"]("a,i") * u1_bb_vo("b,j");

        // sharedOps["abab_vvoo_484"] += 1.000000 tempOp321_aa_vo("a,i") u1_bb_vo("b,j")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["abab_vvoo_484"]("a,b,i,j") = sharedOps["aa_vo_320"]("a,i") * u1_bb_vo("b,j");

        // sharedOps["aaaa_vvoo_493"] += 1.000000 tempOp321_aa_vo("a,i") u1_aa_vo("b,j")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_493"]("a,b,i,j") = sharedOps["aa_vo_320"]("a,i") * u1_aa_vo("b,j");

        // sharedOps["abab_vvoo_491"] += 1.000000 tempOp322_bb_vo("a,i") u1_aa_vo("b,j")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["abab_vvoo_491"]("b,a,j,i") = sharedOps["bb_vo_321"]("a,i") * u1_aa_vo("b,j");

        // sharedOps["bbbb_vvoo_492"] += 1.000000 tempOp322_bb_vo("a,i") u1_bb_vo("b,j")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_492"]("a,b,i,j") = sharedOps["bb_vo_321"]("a,i") * u1_bb_vo("b,j");

        // sharedOps["abab_vvoo_483"] += 1.000000 tempOp323_bb_vo("a,i") u1_aa_vo("b,j")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["abab_vvoo_483"]("b,a,j,i") = sharedOps["bb_vo_322"]("a,i") * u1_aa_vo("b,j");

        // sharedOps["bbbb_vvoo_490"] += 1.000000 tempOp323_bb_vo("a,i") u1_bb_vo("b,j")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_490"]("a,b,i,j") = sharedOps["bb_vo_322"]("a,i") * u1_bb_vo("b,j");

        // sharedOps["abab_vvoo_489"] += 1.000000 u1_aa_vo("b,j") tempOp324_bb_vo("a,i")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["abab_vvoo_489"]("b,a,j,i") = u1_aa_vo("b,j") * sharedOps["bb_vo_323"]("a,i");

        // sharedOps["bbbb_vvoo_498"] += 1.000000 tempOp324_bb_vo("a,i") u1_bb_vo("b,j")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["bbbb_vvoo_498"]("a,b,i,j") = sharedOps["bb_vo_323"]("a,i") * u1_bb_vo("b,j");

        // sharedOps["aaaa_vvoo_486"] += 1.000000 tempOp325_aa_vo("a,i") u1_aa_vo("b,j")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_486"]("a,b,i,j") = sharedOps["aa_vo_324"]("a,i") * u1_aa_vo("b,j");

        // sharedOps["abab_vvoo_487"] += 1.000000 u1_bb_vo("b,j") tempOp325_aa_vo("a,i")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["abab_vvoo_487"]("a,b,i,j") = u1_bb_vo("b,j") * sharedOps["aa_vo_324"]("a,i");

        // sharedOps["abab_vvoo_485"] += 1.000000 u1_bb_vo("b,j") tempOp326_aa_vo("a,i")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["abab_vvoo_485"]("a,b,i,j") = u1_bb_vo("b,j") * sharedOps["aa_vo_325"]("a,i");

        // sharedOps["aaaa_vvoo_488"] += 1.000000 tempOp326_aa_vo("a,i") u1_aa_vo("b,j")
        // flops: o2v2: 2 | mem: o2v2: 2,
        sharedOps["aaaa_vvoo_488"]("a,b,i,j") = sharedOps["aa_vo_325"]("a,i") * u1_aa_vo("b,j");
    }

}


double LambdaDriver::build_residuals() {
    // Get cavity information
    double w0 = cc_wfn_->cavity_frequency_[2];
    double coupling_factor_z = w0 * cc_wfn_->cavity_coupling_strength_[2];

    TA::TArrayD dp_aa_oo, dp_bb_oo, dp_aa_ov, dp_bb_ov, dp_aa_vo, dp_bb_vo, dp_aa_vv, dp_bb_vv;

    TArrayMap Dip_blks_ = cc_wfn_->Dip_blks_;

    // get dipole integrals
    dp_aa_oo("i, j") = coupling_factor_z * Dip_blks_["dz_aa_oo"]("i, j");
    dp_aa_ov("i, a") = coupling_factor_z * Dip_blks_["dz_aa_ov"]("i, a");
    dp_aa_vo("a, i") = coupling_factor_z * Dip_blks_["dz_aa_vo"]("a, i");
    dp_aa_vv("a, b") = coupling_factor_z * Dip_blks_["dz_aa_vv"]("a, b");

    dp_bb_oo("i, j") = coupling_factor_z * Dip_blks_["dz_bb_oo"]("i, j");
    dp_bb_ov("i, a") = coupling_factor_z * Dip_blks_["dz_bb_ov"]("i, a");
    dp_bb_vo("a, i") = coupling_factor_z * Dip_blks_["dz_bb_vo"]("a, i");
    dp_bb_vv("a, b") = coupling_factor_z * Dip_blks_["dz_bb_vv"]("a, b");

    // get residual tensors
    TArrayD &rl1_aa_vo = L_residuals_["l1_aa"];
    TArrayD &rl1_bb_vo = L_residuals_["l1_bb"];
    TArrayD &rm1_aa_vo = L_residuals_["m1_aa"];
    TArrayD &rm1_bb_vo = L_residuals_["m1_bb"];
    TArrayD &rl2_aaaa_vvoo = L_residuals_["l2_aaaa"];
    TArrayD &rl2_abab_vvoo = L_residuals_["l2_abab"];
    TArrayD &rl2_bbbb_vvoo = L_residuals_["l2_bbbb"];
    TArrayD &rm2_aaaa_vvoo = L_residuals_["m2_aaaa"];
    TArrayD &rm2_abab_vvoo = L_residuals_["m2_abab"];
    TArrayD &rm2_bbbb_vvoo = L_residuals_["m2_bbbb"];
    double rm0 = 0.0;

    // initialize to zero
    double energy = 0;
    zero_tiles(rl1_aa_vo);
    zero_tiles(rl1_bb_vo);
    zero_tiles(rl2_aaaa_vvoo);
    zero_tiles(rl2_abab_vvoo);
    zero_tiles(rl2_bbbb_vvoo);

    if ( include_u1_){
        zero_tiles(rm1_aa_vo);
        zero_tiles(rm1_bb_vo);
    }
    if ( include_u2_){
        zero_tiles(rm2_aaaa_vvoo);
        zero_tiles(rm2_abab_vvoo);
        zero_tiles(rm2_bbbb_vvoo);
    }

    // initialize coherent state amplitudes to zero
    TA::TArrayD crl1_aa_vo     = HelperD::makeTensor(world_, {va_, oa_}, true);
    TA::TArrayD crl1_bb_vo     = HelperD::makeTensor(world_, {vb_, ob_}, true);
    TA::TArrayD crl2_aaaa_vvoo = HelperD::makeTensor(world_, {va_, va_, oa_, oa_}, true);
    TA::TArrayD crl2_abab_vvoo = HelperD::makeTensor(world_, {va_, vb_, oa_, ob_}, true);
    TA::TArrayD crl2_bbbb_vvoo = HelperD::makeTensor(world_, {vb_, vb_, ob_, ob_}, true);

    double cenergy = 0.0;
    double crm0 = 0.0;
    TA::TArrayD crm1_aa_vo;
    TA::TArrayD crm1_bb_vo;
    TA::TArrayD crm2_aaaa_vvoo;
    TA::TArrayD crm2_abab_vvoo;
    TA::TArrayD crm2_bbbb_vvoo;

    if (include_u1_){
        crm1_aa_vo = HelperD::makeTensor(world_, {va_, oa_}, true);
        crm1_bb_vo = HelperD::makeTensor(world_, {vb_, ob_}, true);
    }
    if (include_u2_){
        crm2_aaaa_vvoo = HelperD::makeTensor(world_, {va_, va_, oa_, oa_}, true);
        crm2_abab_vvoo = HelperD::makeTensor(world_, {va_, vb_, oa_, ob_}, true);
        crm2_bbbb_vvoo = HelperD::makeTensor(world_, {vb_, vb_, ob_, ob_}, true);
    }

    /// get reference amplitudes

    // extract scalar amplitude
    double& u0 = cc_wfn_->scalar_amps_["u0"];

    // extract 1-body cluster amplitudes
    TA::TArrayD &t1_aa_vo = T_amplitudes_["t1_aa"];
    TA::TArrayD &t1_bb_vo = T_amplitudes_["t1_bb"];
    TA::TArrayD &u1_aa_vo = T_amplitudes_["u1_aa"];
    TA::TArrayD &u1_bb_vo = T_amplitudes_["u1_bb"];

    // extract 2-body cluster amplitudes
    TA::TArrayD &t2_aaaa_vvoo = T_amplitudes_["t2_aaaa"];
    TA::TArrayD &t2_abab_vvoo = T_amplitudes_["t2_abab"];
    TA::TArrayD &t2_bbbb_vvoo = T_amplitudes_["t2_bbbb"];
    TA::TArrayD &u2_aaaa_vvoo = T_amplitudes_["u2_aaaa"];
    TA::TArrayD &u2_abab_vvoo = T_amplitudes_["u2_abab"];
    TA::TArrayD &u2_bbbb_vvoo = T_amplitudes_["u2_bbbb"];

    // extract lambda amplitudes

    // extract scalar amplitude
    double& m0 = L_scalar_amps_["m0"];
    if(include_u0_) {
        world_.gop.fence();
        HelperD::forall(L_amplitudes_["m0"], [&m0](auto &tile, auto &x){
            m0 = tile[x];
        });
        world_.gop.fence();
    }

    // extract 1-body amplitudes
    TA::TArrayD &l1_aa_ov = L_amplitudes_["l1_aa"];
    TA::TArrayD &l1_bb_ov = L_amplitudes_["l1_bb"];
    TA::TArrayD &m1_aa_ov = L_amplitudes_["m1_aa"];
    TA::TArrayD &m1_bb_ov = L_amplitudes_["m1_bb"];

    // extract 2-body amplitudes
    TA::TArrayD &l2_aaaa_oovv = L_amplitudes_["l2_aaaa"];
    TA::TArrayD &l2_abab_oovv = L_amplitudes_["l2_abab"];
    TA::TArrayD &l2_bbbb_oovv = L_amplitudes_["l2_bbbb"];
    TA::TArrayD &m2_aaaa_oovv = L_amplitudes_["m2_aaaa"];
    TA::TArrayD &m2_abab_oovv = L_amplitudes_["m2_abab"];
    TA::TArrayD &m2_bbbb_oovv = L_amplitudes_["m2_bbbb"];

    // grab integrals
    TArrayMap &V_blks_ = cc_wfn_->V_blks_;
    TArrayMap &F_blks_ = cc_wfn_->F_blks_;
    TArrayMap &Id_blks_ = cc_wfn_->F_blks_;

    bool include_m0_ = include_u0_;
    bool include_m1_ = include_u1_;
    bool include_m2_ = include_u2_;

    world_.gop.fence();

    // compute energy and residuals
    {

        TA::TArrayD tempPerm_aaaa_vvoo = HelperD::makeTensor(world_,{va_,va_,oa_,oa_}, true);
        TA::TArrayD tempPerm_bbbb_vvoo = HelperD::makeTensor(world_,{vb_,vb_,ob_,ob_}, true);

// %%%% Prepare Temp Scalars %%%%
        double scalar_tmp_100, scalar_tmp_101, scalar_tmp_102, scalar_tmp_103, scalar_tmp_104, scalar_tmp_105, scalar_tmp_106,
                scalar_tmp_107, scalar_tmp_108, scalar_tmp_109, scalar_tmp_110, scalar_tmp_111, scalar_tmp_112, scalar_tmp_113, scalar_tmp_114,
                scalar_tmp_115, scalar_tmp_116, scalar_tmp_117, scalar_tmp_118, scalar_tmp_119, scalar_tmp_120, scalar_tmp_121, scalar_tmp_122,
                scalar_tmp_123, scalar_tmp_124, scalar_tmp_125, scalar_tmp_126, scalar_tmp_127, scalar_tmp_128, scalar_tmp_129, scalar_tmp_130,
                scalar_tmp_131, scalar_tmp_132, scalar_tmp_133, scalar_tmp_134, scalar_tmp_135, scalar_tmp_136, scalar_tmp_137, scalar_tmp_138,
                scalar_tmp_139, scalar_tmp_140, scalar_tmp_141, scalar_tmp_142, scalar_tmp_143, scalar_tmp_144, scalar_tmp_145, scalar_tmp_146,
                scalar_tmp_147, scalar_tmp_148, scalar_tmp_149, scalar_tmp_150, scalar_tmp_151, scalar_tmp_152, scalar_tmp_153, scalar_tmp_154,
                scalar_tmp_155, scalar_tmp_156, scalar_tmp_157, scalar_tmp_158, scalar_tmp_159, scalar_tmp_160, scalar_tmp_161, scalar_tmp_162,
                scalar_tmp_163, scalar_tmp_164, scalar_tmp_165, scalar_tmp_166, scalar_tmp_167, scalar_tmp_168, scalar_tmp_169, scalar_tmp_170,
                scalar_tmp_171, scalar_tmp_172, scalar_tmp_173, scalar_tmp_174, scalar_tmp_175, scalar_tmp_176, scalar_tmp_177, scalar_tmp_178,
                scalar_tmp_201, scalar_tmp_202, scalar_tmp_203, scalar_tmp_204, scalar_tmp_205, scalar_tmp_206, scalar_tmp_207, scalar_tmp_208,
                scalar_tmp_209, scalar_tmp_210, scalar_tmp_211, scalar_tmp_212, scalar_tmp_213, scalar_tmp_214, scalar_tmp_215, scalar_tmp_216,
                scalar_tmp_217, scalar_tmp_218, scalar_tmp_219, scalar_tmp_220, scalar_tmp_221, scalar_tmp_222, scalar_tmp_223, scalar_tmp_224,
                scalar_tmp_225, scalar_tmp_226, scalar_tmp_227, scalar_tmp_228, scalar_tmp_229, scalar_tmp_230, scalar_tmp_46, scalar_tmp_47,
                scalar_tmp_48, scalar_tmp_49, scalar_tmp_50, scalar_tmp_51, scalar_tmp_52, scalar_tmp_53, scalar_tmp_54, scalar_tmp_59,
                scalar_tmp_60, scalar_tmp_61, scalar_tmp_62, scalar_tmp_67, scalar_tmp_68, scalar_tmp_69, scalar_tmp_70, scalar_tmp_71,
                scalar_tmp_72, scalar_tmp_73, scalar_tmp_74, scalar_tmp_75, scalar_tmp_76, scalar_tmp_77, scalar_tmp_78, scalar_tmp_79,
                scalar_tmp_80, scalar_tmp_81, scalar_tmp_82, scalar_tmp_83, scalar_tmp_84, scalar_tmp_85, scalar_tmp_86, scalar_tmp_87,
                scalar_tmp_88, scalar_tmp_89, scalar_tmp_90, scalar_tmp_91, scalar_tmp_92, scalar_tmp_93, scalar_tmp_94, scalar_tmp_95;

// %%%% Prepare Scalars %%%%
        double scalar_0, scalar_1, scalar_2, scalar_3, scalar_4, scalar_5, scalar_6,
                scalar_7, scalar_8, scalar_9, scalar_10, scalar_11, scalar_12, scalar_13, scalar_14,
                scalar_15, scalar_16, scalar_17, scalar_18, scalar_19, scalar_20, scalar_21, scalar_22,
                scalar_23, scalar_24, scalar_25, scalar_26, scalar_27, scalar_28, scalar_29, scalar_30,
                scalar_31, scalar_32, scalar_33, scalar_34, scalar_35, scalar_36, scalar_37, scalar_38,
                scalar_39, scalar_40, scalar_41, scalar_42, scalar_43, scalar_44, scalar_45, scalar_46,
                scalar_47, scalar_48, scalar_49, scalar_50, scalar_51, scalar_52, scalar_53, scalar_54,
                scalar_55, scalar_56, scalar_57, scalar_58, scalar_59, scalar_60, scalar_61, scalar_62,
                scalar_63, scalar_64, scalar_65, scalar_66, scalar_67, scalar_68, scalar_69, scalar_70,
                scalar_71, scalar_72, scalar_73, scalar_74, scalar_75, scalar_76, scalar_77, scalar_78,
                scalar_79, scalar_80, scalar_81, scalar_82, scalar_83, scalar_84, scalar_85, scalar_86,
                scalar_87, scalar_88, scalar_89, scalar_90, scalar_91, scalar_92, scalar_93, scalar_94,
                scalar_95, scalar_96, scalar_97, scalar_98, scalar_99, scalar_100, scalar_101, scalar_102,
                scalar_103, scalar_104, scalar_105, scalar_106, scalar_107, scalar_108, scalar_109, scalar_110,
                scalar_111, scalar_112, scalar_113, scalar_114, scalar_115, scalar_116, scalar_117, scalar_118,
                scalar_119, scalar_120, scalar_121, scalar_122, scalar_123, scalar_124, scalar_125, scalar_126,
                scalar_127, scalar_128, scalar_129, scalar_130, scalar_131, scalar_132, scalar_133, scalar_134,
                scalar_135, scalar_136, scalar_137, scalar_138, scalar_139, scalar_140, scalar_141, scalar_142,
                scalar_143, scalar_144, scalar_145, scalar_146, scalar_147, scalar_148, scalar_149, scalar_150,
                scalar_151, scalar_152, scalar_153, scalar_154, scalar_155, scalar_156, scalar_157, scalar_158,
                scalar_159, scalar_160, scalar_161, scalar_162, scalar_163, scalar_164, scalar_165, scalar_166,
                scalar_167, scalar_168, scalar_169, scalar_170, scalar_171, scalar_172, scalar_173, scalar_174,
                scalar_175, scalar_176, scalar_177, scalar_178, scalar_179, scalar_180, scalar_181, scalar_182,
                scalar_183, scalar_184, scalar_185, scalar_186, scalar_187, scalar_188, scalar_189, scalar_190,
                scalar_191, scalar_192, scalar_193, scalar_194, scalar_195, scalar_196, scalar_197, scalar_198,
                scalar_199, scalar_200, scalar_201, scalar_202, scalar_203, scalar_204, scalar_205, scalar_206,
                scalar_207, scalar_208, scalar_209, scalar_210, scalar_211, scalar_212, scalar_213, scalar_214,
                scalar_215, scalar_216, scalar_217, scalar_218, scalar_219, scalar_220, scalar_221, scalar_222,
                scalar_223, scalar_224, scalar_225, scalar_226, scalar_227, scalar_228, scalar_229, scalar_230,
                scalar_231, scalar_232, scalar_233, scalar_234, scalar_235, scalar_236, scalar_237, scalar_238,
                scalar_239, scalar_240, scalar_241, scalar_242, scalar_243, scalar_244, scalar_245, scalar_246,
                scalar_247, scalar_248, scalar_249, scalar_250, scalar_251, scalar_252, scalar_253, scalar_254,
                scalar_255, scalar_256, scalar_257, scalar_258, scalar_259, scalar_260, scalar_261, scalar_262,
                scalar_263, scalar_264, scalar_265, scalar_266, scalar_267, scalar_268, scalar_269, scalar_270,
                scalar_271, scalar_272, scalar_273, scalar_274, scalar_275, scalar_276, scalar_277, scalar_278,
                scalar_279, scalar_280, scalar_281, scalar_282, scalar_283, scalar_284, scalar_285, scalar_286,
                scalar_287, scalar_288, scalar_289, scalar_290, scalar_291, scalar_292, scalar_293, scalar_294,
                scalar_295, scalar_296, scalar_297, scalar_298, scalar_299, scalar_300, scalar_301, scalar_302,
                scalar_303, scalar_304, scalar_305, scalar_306, scalar_307, scalar_308, scalar_309, scalar_310,
                scalar_311, scalar_312, scalar_313, scalar_314, scalar_315, scalar_316, scalar_317, scalar_318,
                scalar_319, scalar_320, scalar_321, scalar_322, scalar_323, scalar_324, scalar_325, scalar_326,
                scalar_327, scalar_328, scalar_329, scalar_330, scalar_331, scalar_332, scalar_333, scalar_334,
                scalar_335, scalar_336, scalar_337, scalar_338, scalar_339, scalar_340, scalar_341, scalar_342,
                scalar_343, scalar_344, scalar_345, scalar_346, scalar_347, scalar_348, scalar_349, scalar_350,
                scalar_351, scalar_352, scalar_353, scalar_354, scalar_355, scalar_356, scalar_357, scalar_358,
                scalar_359, scalar_360, scalar_361, scalar_362, scalar_363, scalar_364, scalar_365, scalar_366,
                scalar_367, scalar_368, scalar_369, scalar_370, scalar_371, scalar_372, scalar_373, scalar_374,
                scalar_375, scalar_376, scalar_377, scalar_378, scalar_379, scalar_380, scalar_381, scalar_382,
                scalar_383, scalar_384, scalar_385, scalar_386, scalar_387, scalar_388, scalar_389, scalar_390,
                scalar_391;

// %%%% Assign Scalars %%%%
        scalar_0 = dot(F_blks_["aa_oo"]("o0,o1"), Id_blks_["aa_oo"]("o0,o1"));
        scalar_1 = dot(F_blks_["bb_oo"]("o0,o1"), Id_blks_["bb_oo"]("o0,o1"));
        if (include_u1_)
            scalar_2 = dot(F_blks_["aa_ov"]("i,a"), u1_aa_vo("a,i"));
        if (include_u1_)
            scalar_3 = dot(F_blks_["bb_ov"]("i,a"), u1_bb_vo("a,i"));
        scalar_4 = dot(dp_aa_oo("o0,o1"), Id_blks_["aa_oo"]("o0,o1"));
        scalar_5 = dot(dp_bb_oo("o0,o1"), Id_blks_["bb_oo"]("o0,o1"));
        if (include_u1_)
            scalar_6 = dot(dp_aa_ov("i,a"), u1_aa_vo("a,i"));
        if (include_u1_)
            scalar_7 = dot(dp_bb_ov("i,a"), u1_bb_vo("a,i"));
        scalar_8 = dot(V_blks_["aaaa_oooo"]("o0,o1,o2,o3"), Id_blks_["aaaa_oooo"]("o0,o1,o2,o3"));
        scalar_9 = dot(V_blks_["abab_oooo"]("o0,o1,o2,o3"), Id_blks_["abab_oooo"]("o0,o1,o2,o3"));
        scalar_10 = dot(V_blks_["bbbb_oooo"]("o0,o1,o2,o3"), Id_blks_["bbbb_oooo"]("o0,o1,o2,o3"));
        if (include_u2_)
            scalar_11 = dot(V_blks_["aaaa_oovv"]("j,i,a,b"), u2_aaaa_vvoo("a,b,j,i"));
        if (include_u2_)
            scalar_12 = dot(V_blks_["abab_oovv"]("j,i,a,b"), u2_abab_vvoo("a,b,j,i"));
        if (include_u2_)
            scalar_13 = dot(V_blks_["bbbb_oovv"]("j,i,a,b"), u2_bbbb_vvoo("a,b,j,i"));
        scalar_14 = dot(V_blks_["aaaa_oovv"]("j,i,a,b"), t2_aaaa_vvoo("a,b,j,i"));
        scalar_15 = dot(V_blks_["abab_oovv"]("j,i,a,b"), t2_abab_vvoo("a,b,j,i"));
        scalar_16 = dot(V_blks_["bbbb_oovv"]("j,i,a,b"), t2_bbbb_vvoo("a,b,j,i"));
        if (include_m2_ && include_u2_)
            scalar_17 = dot(sharedOps["abab_vvoo_397"]("b,a,i,j"), m2_abab_oovv("i,j,b,a"));
        scalar_18 = dot(sharedOps["abab_vvoo_398"]("b,a,i,j"), l2_abab_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_19 = dot(sharedOps["abab_vvoo_104"]("a,b,j,i"), m2_abab_oovv("j,i,a,b"));
        if (include_m2_ && include_u2_)
            scalar_20 = dot(u2_abab_vvoo("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_u2_)
            scalar_21 = dot(l2_abab_oovv("k,j,a,b"), u2_abab_vvoo("a,b,k,j"));
        if (include_m2_ && include_u2_)
            scalar_22 = dot(sharedOps["abab_vvoo_10"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        scalar_23 = dot(l2_abab_oovv("j,i,b,a"), sharedOps["abab_vvoo_105"]("b,a,j,i"));
        if (include_m2_ && include_u2_)
            scalar_24 = dot(m2_abab_oovv("i,j,b,a"), sharedOps["abab_vvoo_106"]("b,a,i,j"));
        scalar_25 = dot(l2_abab_oovv("i,j,a,b"), sharedOps["abab_vvoo_11"]("a,b,i,j"));
        if (include_m2_ && include_u2_)
            scalar_26 = dot(u2_aaaa_vvoo("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        if (include_u2_)
            scalar_27 = dot(u2_bbbb_vvoo("a,b,k,j"), l2_bbbb_oovv("k,j,a,b"));
        if (include_u2_)
            scalar_28 = dot(l2_aaaa_oovv("k,j,a,b"), u2_aaaa_vvoo("a,b,k,j"));
        if (include_m2_ && include_u2_)
            scalar_29 = dot(u2_bbbb_vvoo("b,a,i,j"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_30 = dot(sharedOps["abab_vvoo_411"]("a,b,j,i"), m2_abab_oovv("j,i,a,b"));
        if (include_m2_ && include_u2_)
            scalar_31 = dot(m2_abab_oovv("i,j,a,b"), sharedOps["abab_vvoo_407"]("a,b,i,j"));
        if (include_m2_ && include_u0_ && include_u2_)
            scalar_32 = dot(sharedOps["abab_vvoo_156"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_u0_)
            scalar_33 = dot(l2_abab_oovv("j,i,b,a"), sharedOps["abab_vvoo_155"]("b,a,j,i"));
        scalar_34 = dot(sharedOps["abab_vvoo_406"]("a,b,j,i"), l2_abab_oovv("j,i,a,b"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_35 = dot(m2_abab_oovv("j,i,b,a"), sharedOps["abab_vvoo_423"]("b,a,j,i"));
        if (include_u0_)
            scalar_36 = dot(l2_abab_oovv("j,i,a,b"), sharedOps["abab_vvoo_220"]("a,b,j,i"));
        if (include_m2_ && include_u2_)
            scalar_37 = dot(sharedOps["abab_vvoo_429"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u0_ && include_u2_)
            scalar_38 = dot(m2_abab_oovv("j,i,b,a"), sharedOps["abab_vvoo_221"]("b,a,j,i"));
        scalar_39 = dot(sharedOps["abab_vvoo_410"]("b,a,i,j"), l2_abab_oovv("i,j,b,a"));
        scalar_40 = dot(sharedOps["abab_vvoo_425"]("a,b,i,j"), l2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u2_)
            scalar_41 = dot(m2_abab_oovv("i,j,a,b"), sharedOps["abab_vvoo_422"]("a,b,i,j"));
        if (include_m2_ && include_u0_ && include_u1_ && include_u2_)
            scalar_42 = dot(sharedOps["abab_vvoo_218"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_43 = dot(m2_abab_oovv("j,i,a,b"), sharedOps["abab_vvoo_110"]("a,b,j,i"));
        scalar_44 = dot(sharedOps["abab_vvoo_426"]("b,a,j,i"), l2_abab_oovv("j,i,b,a"));
        if (include_u0_)
            scalar_45 = dot(l2_abab_oovv("i,j,a,b"), sharedOps["abab_vvoo_158"]("a,b,i,j"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_46 = dot(m2_abab_oovv("i,j,a,b"), sharedOps["abab_vvoo_109"]("a,b,i,j"));
        if (include_m2_ && include_u0_ && include_u2_)
            scalar_47 = dot(sharedOps["abab_vvoo_157"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        if (include_m2_ && include_u0_ && include_u1_ && include_u2_)
            scalar_48 = dot(sharedOps["abab_vvoo_210"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        if (include_m2_ && include_u0_ && include_u1_ && include_u2_)
            scalar_49 = dot(m2_abab_oovv("j,i,a,b"), sharedOps["abab_vvoo_209"]("a,b,j,i"));
        if (include_m2_ && include_u2_)
            scalar_50 = dot(m2_abab_oovv("j,i,a,b"), sharedOps["abab_vvoo_424"]("a,b,j,i"));
        if (include_u0_)
            scalar_51 = dot(sharedOps["abab_vvoo_222"]("a,b,i,j"), l2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u0_ && include_u1_ && include_u2_)
            scalar_52 = dot(m2_abab_oovv("j,i,b,a"), sharedOps["abab_vvoo_217"]("b,a,j,i"));
        if (include_m2_ && include_u0_ && include_u2_)
            scalar_53 = dot(m2_abab_oovv("i,j,b,a"), sharedOps["abab_vvoo_219"]("b,a,i,j"));
        scalar_54 = dot(V_blks_["abab_vvoo"]("a,b,j,i"), l2_abab_oovv("j,i,a,b"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_55 = dot(m2_abab_oovv("i,j,a,b"), sharedOps["abab_vvoo_427"]("a,b,i,j"));
        if (include_m2_ && include_u2_)
            scalar_56 = dot(m2_abab_oovv("j,i,b,a"), sharedOps["abab_vvoo_428"]("b,a,j,i"));
        if (include_m2_ && include_u2_)
            scalar_57 = dot(m2_abab_oovv("i,j,b,a"), sharedOps["abab_vvoo_408"]("b,a,i,j"));
        if (include_m2_ && include_u2_)
            scalar_58 = dot(m2_abab_oovv("j,i,b,a"), sharedOps["abab_vvoo_409"]("b,a,j,i"));
        if (include_u1_)
            scalar_59 = dot(sharedOps["abab_vvoo_209"]("a,b,i,j"), l2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u2_)
            scalar_60 = dot(m2_abab_oovv("j,i,a,b"), sharedOps["abab_vvoo_158"]("a,b,j,i"));
        if (include_m2_ && include_u2_)
            scalar_61 = dot(sharedOps["abab_vvoo_254"]("b,a,i,j"), m2_abab_oovv("i,j,b,a"));
        if (include_u2_)
            scalar_62 = dot(sharedOps["abab_vvoo_156"]("c,a,i,j"), l2_abab_oovv("i,j,c,a"));
        if (include_u0_)
            scalar_63 = dot(l2_bbbb_oovv("i,j,b,a"), sharedOps["bbbb_vvoo_238"]("b,a,j,i"));
        if (include_m2_ && include_u2_)
            scalar_64 = dot(m2_abab_oovv("j,i,b,a"), sharedOps["abab_vvoo_155"]("b,a,j,i"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_65 = dot(sharedOps["abab_vvoo_175"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        scalar_66 = dot(l2_abab_oovv("i,j,b,a"), sharedOps["abab_vvoo_370"]("b,a,i,j"));
        if (include_m2_ && include_u2_)
            scalar_67 = dot(sharedOps["abab_vvoo_220"]("a,b,j,i"), m2_abab_oovv("j,i,a,b"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_68 = dot(m2_abab_oovv("j,i,a,b"), sharedOps["abab_vvoo_249"]("a,b,j,i"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_69 = dot(sharedOps["abab_vvoo_247"]("a,b,j,i"), m2_abab_oovv("j,i,a,b"));
        if (include_m2_ && include_u2_)
            scalar_70 = dot(sharedOps["aaaa_vvoo_412"]("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_71 = dot(m2_abab_oovv("j,i,a,b"), sharedOps["abab_vvoo_180"]("a,b,j,i"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_72 = dot(sharedOps["abab_vvoo_246"]("b,a,i,j"), m2_abab_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_73 = dot(sharedOps["abab_vvoo_222"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_74 = dot(m2_abab_oovv("i,j,a,b"), sharedOps["abab_vvoo_245"]("a,b,i,j"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_75 = dot(sharedOps["abab_vvoo_190"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        if (include_m2_ && include_u2_)
            scalar_76 = dot(sharedOps["bbbb_vvoo_440"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        scalar_77 = dot(sharedOps["abab_vvoo_189"]("a,b,i,j"), l2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u0_ && include_u2_)
            scalar_78 = dot(m2_bbbb_oovv("i,j,b,a"), sharedOps["bbbb_vvoo_170"]("b,a,i,j"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_79 = dot(m2_abab_oovv("i,j,b,a"), sharedOps["abab_vvoo_252"]("b,a,i,j"));
        if (include_u1_)
            scalar_80 = dot(l2_abab_oovv("i,j,b,a"), sharedOps["abab_vvoo_218"]("b,a,i,j"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_81 = dot(sharedOps["abab_vvoo_176"]("b,a,i,j"), m2_abab_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_82 = dot(sharedOps["abab_vvoo_174"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_83 = dot(sharedOps["abab_vvoo_260"]("a,b,j,i"), m2_abab_oovv("j,i,a,b"));
        scalar_84 = dot(sharedOps["abab_vvoo_371"]("b,a,j,i"), l2_abab_oovv("j,i,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_85 = dot(sharedOps["abab_vvoo_173"]("a,b,j,i"), m2_abab_oovv("j,i,a,b"));
        scalar_86 = dot(l2_abab_oovv("j,i,a,b"), sharedOps["abab_vvoo_162"]("a,b,j,i"));
        if (include_m2_ && include_u2_)
            scalar_87 = dot(sharedOps["abab_vvoo_258"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        if (include_m2_ && include_u0_ && include_u2_)
            scalar_88 = dot(sharedOps["aaaa_vvoo_237"]("b,a,j,i"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_89 = dot(m2_abab_oovv("j,i,b,a"), sharedOps["abab_vvoo_182"]("b,a,j,i"));
        if (include_u0_)
            scalar_90 = dot(sharedOps["aaaa_vvoo_172"]("b,a,i,j"), l2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u0_ && include_u1_ && include_u2_)
            scalar_91 = dot(sharedOps["bbbb_vvoo_236"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u0_ && include_u2_)
            scalar_92 = dot(sharedOps["bbbb_vvoo_241"]("b,a,i,j"), m2_bbbb_oovv("i,j,b,a"));
        if (include_u1_)
            scalar_93 = dot(l2_abab_oovv("i,j,b,a"), sharedOps["abab_vvoo_210"]("b,a,i,j"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_94 = dot(sharedOps["abab_vvoo_248"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        if (include_u0_)
            scalar_95 = dot(sharedOps["aaaa_vvoo_239"]("b,a,i,j"), l2_aaaa_oovv("i,j,b,a"));
        scalar_96 = dot(l2_abab_oovv("i,j,a,b"), sharedOps["abab_vvoo_233"]("a,b,i,j"));
        scalar_97 = dot(l2_abab_oovv("i,j,b,a"), sharedOps["abab_vvoo_269"]("b,a,i,j"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_98 = dot(m2_abab_oovv("i,j,b,a"), sharedOps["abab_vvoo_441"]("b,a,i,j"));
        if (include_m2_ && include_u2_)
            scalar_99 = dot(m2_aaaa_oovv("i,j,b,a"), sharedOps["aaaa_vvoo_443"]("b,a,j,i"));
        if (include_m2_ && include_u2_)
            scalar_100 = dot(m2_abab_oovv("i,j,b,a"), sharedOps["abab_vvoo_229"]("b,a,i,j"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_101 = dot(sharedOps["abab_vvoo_230"]("a,b,j,i"), m2_abab_oovv("j,i,a,b"));
        if (include_m2_ && include_u2_)
            scalar_102 = dot(m2_abab_oovv("j,i,b,a"), sharedOps["abab_vvoo_159"]("b,a,j,i"));
        if (include_m2_ && include_u0_ && include_u2_)
            scalar_103 = dot(m2_aaaa_oovv("i,j,b,a"), sharedOps["aaaa_vvoo_164"]("b,a,i,j"));
        scalar_104 = dot(sharedOps["abab_vvoo_231"]("b,a,j,i"), l2_abab_oovv("j,i,b,a"));
        scalar_105 = dot(sharedOps["abab_vvoo_181"]("b,a,i,j"), l2_abab_oovv("i,j,b,a"));
        if (include_u2_)
            scalar_106 = dot(l2_abab_oovv("k,i,a,b"), sharedOps["abab_vvoo_219"]("a,b,k,i"));
        if (include_u2_)
            scalar_107 = dot(sharedOps["abab_vvoo_157"]("a,c,i,j"), l2_abab_oovv("i,j,a,c"));
        if (include_u1_)
            scalar_108 = dot(sharedOps["abab_vvoo_217"]("a,b,j,i"), l2_abab_oovv("j,i,a,b"));
        if (include_u0_)
            scalar_109 = dot(l2_bbbb_oovv("i,j,b,a"), sharedOps["bbbb_vvoo_168"]("b,a,i,j"));
        scalar_110 = dot(sharedOps["abab_vvoo_268"]("a,b,j,i"), l2_abab_oovv("j,i,a,b"));
        if (include_u2_)
            scalar_111 = dot(l2_abab_oovv("i,k,a,b"), sharedOps["abab_vvoo_221"]("a,b,i,k"));
        if (include_m2_ && include_u2_)
            scalar_112 = dot(m2_abab_oovv("j,i,a,b"), sharedOps["abab_vvoo_414"]("a,b,j,i"));
        if (include_m2_ && include_u2_)
            scalar_113 = dot(sharedOps["abab_vvoo_160"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_114 = dot(sharedOps["abab_vvoo_250"]("b,a,i,j"), m2_abab_oovv("i,j,b,a"));
        if (include_m2_ && include_u0_ && include_u1_ && include_u2_)
            scalar_115 = dot(m2_aaaa_oovv("i,j,b,a"), sharedOps["aaaa_vvoo_235"]("b,a,j,i"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_116 = dot(m2_abab_oovv("j,i,b,a"), sharedOps["abab_vvoo_251"]("b,a,j,i"));
        scalar_117 = dot(sharedOps["bbbb_vvoo_442"]("b,a,i,j"), l2_bbbb_oovv("i,j,b,a"));
        scalar_118 = dot(l2_abab_oovv("j,i,b,a"), sharedOps["abab_vvoo_161"]("b,a,j,i"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_119 = dot(sharedOps["abab_vvoo_225"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_120 = dot(sharedOps["abab_vvoo_263"]("b,a,i,j"), m2_abab_oovv("i,j,b,a"));
        if (include_m2_ && include_u0_ && include_u1_ && include_u2_)
            scalar_121 = dot(m2_bbbb_oovv("i,j,b,a"), sharedOps["bbbb_vvoo_430"]("a,b,i,j"));
        if (include_m2_ && include_u2_)
            scalar_122 = dot(sharedOps["abab_vvoo_434"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_123 = dot(sharedOps["abab_vvoo_185"]("b,a,i,j"), m2_abab_oovv("i,j,b,a"));
        if (include_m2_ && include_u0_ && include_u1_ && include_u2_)
            scalar_124 = dot(sharedOps["aaaa_vvoo_431"]("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_125 = dot(m2_abab_oovv("j,i,a,b"), sharedOps["abab_vvoo_433"]("a,b,j,i"));
        if (include_m2_ && include_u2_)
            scalar_126 = dot(m2_bbbb_oovv("i,j,b,a"), sharedOps["bbbb_vvoo_419"]("b,a,i,j"));
        if (include_m2_ && include_u2_)
            scalar_127 = dot(m2_aaaa_oovv("i,j,b,a"), sharedOps["aaaa_vvoo_413"]("a,b,i,j"));
        scalar_128 = dot(sharedOps["abab_vvoo_52"]("b,a,j,i"), l2_abab_oovv("j,i,b,a"));
        if (include_m2_ && include_u2_)
            scalar_129 = dot(m2_bbbb_oovv("i,j,b,a"), sharedOps["bbbb_vvoo_415"]("b,a,i,j"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_130 = dot(m2_abab_oovv("i,j,a,b"), sharedOps["abab_vvoo_432"]("a,b,i,j"));
        if (include_m2_ && include_u2_)
            scalar_131 = dot(m2_abab_oovv("j,i,b,a"), sharedOps["abab_vvoo_416"]("b,a,j,i"));
        if (include_m2_ && include_u2_)
            scalar_132 = dot(m2_abab_oovv("j,i,a,b"), sharedOps["abab_vvoo_227"]("a,b,j,i"));
        if (include_m2_ && include_u2_)
            scalar_133 = dot(m2_bbbb_oovv("i,j,b,a"), sharedOps["bbbb_vvoo_438"]("b,a,j,i"));
        if (include_m2_ && include_u2_)
            scalar_134 = dot(m2_aaaa_oovv("i,j,b,a"), sharedOps["aaaa_vvoo_437"]("b,a,i,j"));
        scalar_135 = dot(l2_aaaa_oovv("i,j,b,a"), sharedOps["aaaa_vvoo_439"]("b,a,j,i"));
        scalar_136 = dot(l2_bbbb_oovv("i,j,b,a"), sharedOps["bbbb_vvoo_417"]("b,a,i,j"));
        scalar_137 = dot(l2_aaaa_oovv("i,j,b,a"), sharedOps["aaaa_vvoo_418"]("b,a,i,j"));
        if (include_m1_ && include_u1_)
            scalar_138 = dot(u1_bb_vo("a,i"), m1_bb_ov("i,a"));
        if (include_u1_)
            scalar_139 = dot(l1_bb_ov("j,a"), u1_bb_vo("a,j"));
        if (include_u1_)
            scalar_140 = dot(l1_aa_ov("j,a"), u1_aa_vo("a,j"));
        if (include_m1_ && include_u1_)
            scalar_141 = dot(m1_aa_ov("i,a"), u1_aa_vo("a,i"));
        if (include_m1_ && include_u0_ && include_u1_)
            scalar_142 = dot(sharedOps["aa_vo_360"]("a,i"), m1_aa_ov("i,a"));
        scalar_143 = dot(sharedOps["bb_vo_266"]("a,i"), l1_bb_ov("i,a"));
        scalar_144 = dot(sharedOps["aa_vo_265"]("a,i"), l1_aa_ov("i,a"));
        scalar_145 = dot(sharedOps["aa_vo_184"]("a,i"), l1_aa_ov("i,a"));
        if (include_u0_)
            scalar_146 = dot(l1_aa_ov("i,a"), sharedOps["aa_vo_318"]("a,i"));
        if (include_u0_)
            scalar_147 = dot(l1_aa_ov("i,a"), sharedOps["aa_vo_320"]("a,i"));
        if (include_m1_ && include_u1_ && include_u2_)
            scalar_148 = dot(sharedOps["bb_vo_188"]("a,i"), m1_bb_ov("i,a"));
        if (include_m1_ && include_u0_ && include_u1_)
            scalar_149 = dot(sharedOps["bb_vo_361"]("a,i"), m1_bb_ov("i,a"));
        if (include_u0_)
            scalar_150 = dot(sharedOps["bb_vo_321"]("a,i"), l1_bb_ov("i,a"));
        if (include_m1_ && include_u0_ && include_u1_ && include_u2_)
            scalar_151 = dot(sharedOps["bb_vo_319"]("a,i"), m1_bb_ov("i,a"));
        if (include_m1_ && include_u0_ && include_u1_ && include_u2_)
            scalar_152 = dot(m1_aa_ov("i,a"), sharedOps["aa_vo_324"]("a,i"));
        if (include_m1_ && include_u1_)
            scalar_153 = dot(sharedOps["bb_vo_435"]("a,i"), m1_bb_ov("i,a"));
        scalar_154 = dot(l1_bb_ov("i,a"), sharedOps["bb_vo_183"]("a,i"));
        scalar_155 = dot(dp_aa_vo("a,i"), l1_aa_ov("i,a"));
        if (include_m1_ && include_u1_)
            scalar_156 = dot(sharedOps["bb_vo_499"]("a,i"), m1_bb_ov("i,a"));
        if (include_m1_ && include_u1_ && include_u2_)
            scalar_157 = dot(m1_aa_ov("i,a"), sharedOps["aa_vo_177"]("a,i"));
        if (include_m1_ && include_u0_ && include_u1_ && include_u2_)
            scalar_158 = dot(sharedOps["bb_vo_322"]("a,i"), m1_bb_ov("i,a"));
        if (include_m1_ && include_u1_ && include_u2_)
            scalar_159 = dot(sharedOps["bb_vo_264"]("a,i"), m1_bb_ov("i,a"));
        if (include_m1_ && include_u0_ && include_u1_ && include_u2_)
            scalar_160 = dot(sharedOps["aa_vo_325"]("a,i"), m1_aa_ov("i,a"));
        if (include_m1_ && include_u0_ && include_u1_)
            scalar_161 = dot(m1_bb_ov("i,a"), sharedOps["bb_vo_365"]("a,i"));
        scalar_162 = dot(l1_bb_ov("i,a"), dp_bb_vo("a,i"));
        if (include_m1_ && include_u1_ && include_u2_)
            scalar_163 = dot(sharedOps["aa_vo_255"]("a,i"), m1_aa_ov("i,a"));
        if (include_m1_ && include_u0_ && include_u1_)
            scalar_164 = dot(sharedOps["aa_vo_364"]("a,i"), m1_aa_ov("i,a"));
        if (include_u0_)
            scalar_165 = dot(sharedOps["bb_vo_323"]("a,i"), l1_bb_ov("i,a"));
        if (include_m1_ && include_u1_)
            scalar_166 = dot(m1_aa_ov("i,a"), sharedOps["aa_vo_436"]("a,i"));
        if (include_m1_ && include_u1_)
            scalar_167 = dot(m1_aa_ov("i,a"), sharedOps["aa_vo_500"]("a,i"));
        if (include_m2_ && include_u2_)
            scalar_168 = dot(sharedOps["aaaa_vvoo_131"]("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_169 = dot(sharedOps["bbbb_vvoo_279"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_170 = dot(sharedOps["bbbb_vvoo_278"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_171 = dot(sharedOps["aaaa_vvoo_277"]("b,a,j,i"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m1_ && include_u1_)
            scalar_172 = dot(sharedOps["bb_vo_328"]("a,i"), m1_bb_ov("i,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_173 = dot(sharedOps["bbbb_vvoo_276"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m1_ && include_u1_)
            scalar_174 = dot(sharedOps["bb_vo_321"]("a,i"), m1_bb_ov("i,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_175 = dot(sharedOps["aaaa_vvoo_275"]("b,a,j,i"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m1_ && include_u1_)
            scalar_176 = dot(sharedOps["bb_vo_330"]("a,i"), m1_bb_ov("i,a"));
        if (include_m2_ && include_u2_)
            scalar_177 = dot(sharedOps["aaaa_vvoo_270"]("b,a,j,i"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m1_ && include_u1_)
            scalar_178 = dot(sharedOps["aa_vo_331"]("a,i"), m1_aa_ov("i,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_179 = dot(sharedOps["bbbb_vvoo_281"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        scalar_180 = dot(sharedOps["bb_vo_332"]("a,i"), l1_bb_ov("i,a"));
        scalar_181 = dot(sharedOps["bb_vo_333"]("a,i"), l1_bb_ov("i,a"));
        if (include_m1_ && include_u1_)
            scalar_182 = dot(sharedOps["aa_vo_334"]("a,i"), m1_aa_ov("i,a"));
        if (include_m1_ && include_u1_ && include_u2_)
            scalar_183 = dot(sharedOps["bb_vo_335"]("a,i"), m1_bb_ov("i,a"));
        if (include_m1_ && include_u1_ && include_u2_)
            scalar_184 = dot(sharedOps["bb_vo_338"]("a,i"), m1_bb_ov("i,a"));
        scalar_185 = dot(sharedOps["aa_vo_336"]("a,i"), l1_aa_ov("i,a"));
        if (include_m2_ && include_u2_)
            scalar_186 = dot(sharedOps["bbbb_vvoo_267"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m1_ && include_u1_)
            scalar_187 = dot(sharedOps["aa_vo_253"]("a,i"), m1_aa_ov("i,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_188 = dot(sharedOps["aaaa_vvoo_280"]("b,a,j,i"), m2_aaaa_oovv("i,j,b,a"));
        if (include_u1_)
            scalar_189 = dot(sharedOps["bbbb_vvoo_236"]("b,a,j,i"), l2_bbbb_oovv("i,j,b,a"));
        if (include_u1_)
            scalar_190 = dot(sharedOps["aaaa_vvoo_235"]("b,a,j,i"), l2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_191 = dot(sharedOps["aaaa_vvoo_239"]("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        scalar_192 = dot(sharedOps["aaaa_vvoo_206"]("a,b,i,j"), l2_aaaa_oovv("i,j,b,a"));
        scalar_193 = dot(sharedOps["bbbb_vvoo_256"]("b,a,j,i"), l2_bbbb_oovv("i,j,b,a"));
        scalar_194 = dot(sharedOps["aaaa_vvoo_257"]("b,a,j,i"), l2_aaaa_oovv("i,j,b,a"));
        scalar_195 = dot(sharedOps["bb_vo_204"]("a,i"), l1_bb_ov("i,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_196 = dot(sharedOps["bbbb_vvoo_203"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m1_ && include_u1_ && include_u2_)
            scalar_197 = dot(sharedOps["bb_vo_202"]("a,i"), m1_bb_ov("i,a"));
        if (include_m1_ && include_u1_)
            scalar_198 = dot(sharedOps["aa_vo_329"]("a,i"), m1_aa_ov("i,a"));
        scalar_199 = dot(sharedOps["aa_vo_337"]("a,i"), l1_aa_ov("i,a"));
        if (include_m2_ && include_u2_)
            scalar_200 = dot(sharedOps["aaaa_vvoo_200"]("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m1_ && include_u1_ && include_u2_)
            scalar_201 = dot(sharedOps["aa_vo_199"]("a,i"), m1_aa_ov("i,a"));
        if (include_u2_)
            scalar_202 = dot(sharedOps["aaaa_vvoo_237"]("a,b,i,k"), l2_aaaa_oovv("k,i,a,b"));
        if (include_m2_ && include_u2_)
            scalar_203 = dot(sharedOps["bbbb_vvoo_238"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        if (include_u2_)
            scalar_204 = dot(sharedOps["bbbb_vvoo_241"]("a,b,k,i"), l2_bbbb_oovv("k,i,a,b"));
        scalar_205 = dot(sharedOps["aaaa_vvoo_139"]("b,a,i,j"), l2_aaaa_oovv("i,j,b,a"));
        if (include_m1_ && include_u1_)
            scalar_206 = dot(sharedOps["bb_vo_323"]("a,i"), m1_bb_ov("i,a"));
        if (include_u2_)
            scalar_207 = dot(sharedOps["aa_vo_324"]("b,j"), l1_aa_ov("j,b"));
        scalar_208 = dot(sharedOps["aa_vo_198"]("a,i"), l1_aa_ov("i,a"));
        scalar_209 = dot(sharedOps["bbbb_vvoo_197"]("a,b,i,j"), l2_bbbb_oovv("i,j,b,a"));
        if (include_u2_)
            scalar_210 = dot(sharedOps["aa_vo_325"]("b,j"), l1_aa_ov("j,b"));
        if (include_u2_)
            scalar_211 = dot(sharedOps["bb_vo_322"]("b,j"), l1_bb_ov("j,b"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_212 = dot(sharedOps["abab_vvoo_355"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_213 = dot(sharedOps["abab_vvoo_354"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u2_)
            scalar_214 = dot(sharedOps["bbbb_vvoo_379"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_215 = dot(sharedOps["bbbb_vvoo_352"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_216 = dot(sharedOps["abab_vvoo_351"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_217 = dot(sharedOps["aaaa_vvoo_350"]("b,a,j,i"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_218 = dot(sharedOps["abab_vvoo_349"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_219 = dot(sharedOps["bbbb_vvoo_348"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m1_ && include_u1_)
            scalar_220 = dot(sharedOps["bb_vo_343"]("a,i"), m1_bb_ov("i,a"));
        if (include_m1_ && include_u1_ && include_u2_)
            scalar_221 = dot(sharedOps["aa_vo_342"]("a,i"), m1_aa_ov("i,a"));
        if (include_m1_ && include_u1_ && include_u2_)
            scalar_222 = dot(sharedOps["aa_vo_341"]("a,i"), m1_aa_ov("i,a"));
        if (include_m1_ && include_u1_)
            scalar_223 = dot(sharedOps["bb_vo_340"]("a,i"), m1_bb_ov("i,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_224 = dot(sharedOps["aaaa_vvoo_353"]("a,b,j,i"), m2_aaaa_oovv("i,j,b,a"));
        scalar_225 = dot(sharedOps["bbbb_vvoo_378"]("a,b,i,j"), l2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_226 = dot(sharedOps["bbbb_vvoo_377"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        scalar_227 = dot(sharedOps["bbbb_vvoo_292"]("b,a,i,j"), l2_bbbb_oovv("i,j,b,a"));
        scalar_228 = dot(sharedOps["abab_vvoo_375"]("a,b,i,j"), l2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u2_)
            scalar_229 = dot(sharedOps["abab_vvoo_374"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        scalar_230 = dot(sharedOps["abab_vvoo_373"]("b,a,j,i"), l2_abab_oovv("j,i,b,a"));
        scalar_231 = dot(sharedOps["bbbb_vvoo_372"]("b,a,j,i"), l2_bbbb_oovv("i,j,b,a"));
        scalar_232 = dot(sharedOps["bbbb_vvoo_372"]("a,b,i,j"), l2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_233 = dot(sharedOps["aaaa_vvoo_187"]("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        if (include_u2_)
            scalar_234 = dot(sharedOps["bbbb_vvoo_170"]("a,c,j,i"), l2_bbbb_oovv("j,i,a,c"));
        if (include_u2_)
            scalar_235 = dot(sharedOps["aaaa_vvoo_164"]("a,c,j,i"), l2_aaaa_oovv("j,i,a,c"));
        if (include_m2_ && include_u2_)
            scalar_236 = dot(sharedOps["bbbb_vvoo_168"]("b,a,i,j"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m1_ && include_u1_)
            scalar_237 = dot(sharedOps["bb_vo_369"]("a,i"), m1_bb_ov("i,a"));
        if (include_m1_ && include_u1_)
            scalar_238 = dot(sharedOps["aa_vo_368"]("a,i"), m1_aa_ov("i,a"));
        if (include_u1_)
            scalar_239 = dot(sharedOps["bb_vo_365"]("a,i"), l1_bb_ov("i,a"));
        if (include_u1_)
            scalar_240 = dot(sharedOps["aa_vo_364"]("a,i"), l1_aa_ov("i,a"));
        if (include_m2_ && include_u2_)
            scalar_241 = dot(sharedOps["aaaa_vvoo_172"]("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m1_ && include_u1_)
            scalar_242 = dot(sharedOps["bb_vo_363"]("a,i"), m1_bb_ov("i,a"));
        if (include_m2_ && include_u2_)
            scalar_243 = dot(sharedOps["aaaa_vvoo_376"]("b,a,j,i"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_244 = dot(sharedOps["bbbb_vvoo_122"]("b,a,i,j"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_245 = dot(sharedOps["aaaa_vvoo_123"]("b,a,j,i"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_246 = dot(sharedOps["bbbb_vvoo_124"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        scalar_247 = dot(sharedOps["bbbb_vvoo_127"]("b,a,i,j"), l2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_248 = dot(sharedOps["bbbb_vvoo_128"]("b,a,i,j"), m2_bbbb_oovv("i,j,b,a"));
        scalar_249 = dot(sharedOps["bbbb_vvoo_178"]("b,a,i,j"), l2_bbbb_oovv("i,j,b,a"));
        if (include_m1_ && include_u1_)
            scalar_250 = dot(sharedOps["bb_vo_261"]("a,i"), m1_bb_ov("i,a"));
        if (include_m2_ && include_u2_)
            scalar_251 = dot(sharedOps["bbbb_vvoo_179"]("b,a,i,j"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_252 = dot(sharedOps["abab_vvoo_396"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u2_)
            scalar_253 = dot(sharedOps["abab_vvoo_395"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        scalar_254 = dot(sharedOps["abab_vvoo_394"]("b,a,j,i"), l2_abab_oovv("j,i,b,a"));
        scalar_255 = dot(sharedOps["abab_vvoo_393"]("a,b,i,j"), l2_abab_oovv("i,j,a,b"));
        scalar_256 = dot(sharedOps["aaaa_vvoo_392"]("b,a,j,i"), l2_aaaa_oovv("i,j,b,a"));
        scalar_257 = dot(sharedOps["aaaa_vvoo_391"]("a,b,i,j"), l2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_258 = dot(sharedOps["abab_vvoo_390"]("b,a,i,j"), m2_abab_oovv("i,j,b,a"));
        if (include_u1_)
            scalar_259 = dot(sharedOps["aaaa_vvoo_356"]("a,b,i,j"), l2_aaaa_oovv("j,i,a,b"));
        if (include_u1_)
            scalar_260 = dot(sharedOps["abab_vvoo_357"]("b,a,j,i"), l2_abab_oovv("j,i,b,a"));
        if (include_u1_)
            scalar_261 = dot(sharedOps["bbbb_vvoo_358"]("b,a,j,i"), l2_bbbb_oovv("j,i,a,b"));
        if (include_u1_)
            scalar_262 = dot(sharedOps["abab_vvoo_359"]("a,b,i,j"), l2_abab_oovv("i,j,a,b"));
        if (include_u1_)
            scalar_263 = dot(sharedOps["aa_vo_360"]("a,i"), l1_aa_ov("i,a"));
        scalar_264 = dot(sharedOps["aaaa_vvoo_186"]("b,a,i,j"), l2_aaaa_oovv("i,j,b,a"));
        if (include_u1_)
            scalar_265 = dot(sharedOps["bb_vo_361"]("a,i"), l1_bb_ov("i,a"));
        if (include_m2_ && include_u2_)
            scalar_266 = dot(sharedOps["abab_vvoo_389"]("a,b,j,i"), m2_abab_oovv("j,i,a,b"));
        if (include_m1_ && include_u1_)
            scalar_267 = dot(sharedOps["aa_vo_339"]("a,i"), m1_aa_ov("i,a"));
        scalar_268 = dot(sharedOps["bbbb_vvoo_381"]("a,b,i,j"), l2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_269 = dot(sharedOps["abab_vvoo_382"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u2_)
            scalar_270 = dot(sharedOps["bbbb_vvoo_383"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_271 = dot(sharedOps["abab_vvoo_384"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u2_)
            scalar_272 = dot(sharedOps["aaaa_vvoo_385"]("b,a,j,i"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_273 = dot(sharedOps["abab_vvoo_386"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        if (include_m2_ && include_u2_)
            scalar_274 = dot(sharedOps["aaaa_vvoo_387"]("a,b,i,j"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_275 = dot(sharedOps["bbbb_vvoo_388"]("a,b,i,j"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_276 = dot(sharedOps["abab_vvoo_380"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        scalar_277 = dot(sharedOps["abab_vvoo_78"]("b,a,j,i"), l2_abab_oovv("j,i,b,a"));
        if (include_m1_)
            scalar_278 = dot(dp_bb_vo("a,i"), m1_bb_ov("i,a"));
        if (include_m1_)
            scalar_279 = dot(dp_aa_vo("a,i"), m1_aa_ov("i,a"));
        scalar_280 = dot(sharedOps["bbbb_vvoo_400"]("b,a,i,j"), l2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_281 = dot(sharedOps["bbbb_vvoo_399"]("b,a,i,j"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m1_ && include_u1_)
            scalar_282 = dot(sharedOps["aa_vo_362"]("a,i"), m1_aa_ov("i,a"));
        if (include_m2_ && include_u2_)
            scalar_283 = dot(sharedOps["aaaa_vvoo_402"]("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_284 = dot(sharedOps["aaaa_vvoo_403"]("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        scalar_285 = dot(sharedOps["abab_vvoo_77"]("a,b,j,i"), l2_abab_oovv("j,i,a,b"));
        if (include_m2_ && include_u2_)
            scalar_286 = dot(sharedOps["abab_vvoo_76"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u2_)
            scalar_287 = dot(sharedOps["abab_vvoo_75"]("b,a,i,j"), m2_abab_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_288 = dot(sharedOps["bbbb_vvoo_74"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        scalar_289 = dot(sharedOps["bbbb_vvoo_73"]("a,b,i,j"), l2_bbbb_oovv("i,j,b,a"));
        scalar_290 = dot(sharedOps["bbbb_vvoo_72"]("a,b,i,j"), l2_bbbb_oovv("i,j,b,a"));
        scalar_291 = dot(sharedOps["aaaa_vvoo_401"]("b,a,i,j"), l2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_292 = dot(sharedOps["aaaa_vvoo_460"]("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_293 = dot(sharedOps["aaaa_vvoo_462"]("b,a,j,i"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_294 = dot(sharedOps["abab_vvoo_464"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_m1_ && include_u1_)
            scalar_295 = dot(sharedOps["aa_vo_479"]("a,i"), m1_aa_ov("i,a"));
        if (include_m1_ && include_u1_)
            scalar_296 = dot(sharedOps["bb_vo_480"]("a,i"), m1_bb_ov("i,a"));
        if (include_m1_ && include_u1_)
            scalar_297 = dot(sharedOps["bb_vo_481"]("a,i"), m1_bb_ov("i,a"));
        if (include_m1_ && include_u1_)
            scalar_298 = dot(sharedOps["aa_vo_482"]("a,i"), m1_aa_ov("i,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_299 = dot(m2_abab_oovv("j,i,b,a"), sharedOps["abab_vvoo_483"]("b,a,j,i"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_300 = dot(m2_abab_oovv("i,j,a,b"), sharedOps["abab_vvoo_485"]("a,b,i,j"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_301 = dot(m2_aaaa_oovv("i,j,b,a"), sharedOps["aaaa_vvoo_486"]("a,b,i,j"));
        scalar_302 = dot(V_blks_["bbbb_vvoo"]("b,a,i,j"), l2_bbbb_oovv("i,j,b,a"));
        if (include_u1_)
            scalar_303 = dot(l2_abab_oovv("i,j,a,b"), sharedOps["abab_vvoo_484"]("a,b,i,j"));
        if (include_m2_ && include_u2_)
            scalar_304 = dot(sharedOps["abab_vvoo_83"]("a,b,j,i"), m2_abab_oovv("j,i,a,b"));
        if (include_m2_ && include_u2_)
            scalar_305 = dot(sharedOps["bbbb_vvoo_82"]("a,b,i,j"), m2_bbbb_oovv("i,j,b,a"));
        scalar_306 = dot(sharedOps["abab_vvoo_80"]("a,b,i,j"), l2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u2_)
            scalar_307 = dot(sharedOps["bbbb_vvoo_420"]("b,a,i,j"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_308 = dot(sharedOps["aaaa_vvoo_421"]("a,b,i,j"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_309 = dot(sharedOps["abab_vvoo_90"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_310 = dot(sharedOps["aaaa_vvoo_67"]("b,a,j,i"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_311 = dot(sharedOps["abab_vvoo_65"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_312 = dot(sharedOps["abab_vvoo_66"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_313 = dot(sharedOps["abab_vvoo_63"]("a,b,j,i"), m2_abab_oovv("j,i,a,b"));
        scalar_314 = dot(sharedOps["aaaa_vvoo_57"]("a,b,i,j"), l2_aaaa_oovv("i,j,b,a"));
        scalar_315 = dot(sharedOps["aaaa_vvoo_57"]("b,a,j,i"), l2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_316 = dot(sharedOps["bbbb_vvoo_59"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_317 = dot(sharedOps["abab_vvoo_60"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_318 = dot(sharedOps["abab_vvoo_61"]("b,a,i,j"), m2_abab_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_319 = dot(sharedOps["aaaa_vvoo_62"]("b,a,j,i"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_320 = dot(sharedOps["abab_vvoo_64"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_321 = dot(sharedOps["bbbb_vvoo_68"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        scalar_322 = dot(F_blks_["bb_vo"]("a,i"), l1_bb_ov("i,a"));
        scalar_323 = dot(F_blks_["aa_vo"]("a,i"), l1_aa_ov("i,a"));
        scalar_324 = dot(sharedOps["aaaa_vvoo_84"]("b,a,j,i"), l2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_325 = dot(sharedOps["abab_vvoo_101"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u2_)
            scalar_326 = dot(sharedOps["aaaa_vvoo_97"]("b,a,j,i"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_327 = dot(sharedOps["aaaa_vvoo_91"]("b,a,j,i"), m2_aaaa_oovv("i,j,b,a"));
        scalar_328 = dot(sharedOps["abab_vvoo_89"]("a,b,i,j"), l2_abab_oovv("i,j,a,b"));
        scalar_329 = dot(sharedOps["aaaa_vvoo_88"]("a,b,i,j"), l2_aaaa_oovv("i,j,b,a"));
        scalar_330 = dot(sharedOps["abab_vvoo_87"]("b,a,i,j"), l2_abab_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_331 = dot(sharedOps["abab_vvoo_86"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        if (include_m2_ && include_u2_)
            scalar_332 = dot(sharedOps["abab_vvoo_85"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        if (include_m2_ && include_u2_)
            scalar_333 = dot(sharedOps["aaaa_vvoo_102"]("b,a,j,i"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_334 = dot(sharedOps["bbbb_vvoo_196"]("b,a,i,j"), m2_bbbb_oovv("i,j,b,a"));
        scalar_335 = dot(sharedOps["bbbb_vvoo_28"]("b,a,i,j"), l2_bbbb_oovv("i,j,b,a"));
        scalar_336 = dot(V_blks_["aaaa_vvoo"]("b,a,i,j"), l2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_337 = dot(sharedOps["aaaa_vvoo_459"]("b,a,j,i"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_338 = dot(sharedOps["aaaa_vvoo_445"]("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_339 = dot(sharedOps["bbbb_vvoo_446"]("b,a,i,j"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_340 = dot(sharedOps["bbbb_vvoo_447"]("b,a,i,j"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_341 = dot(sharedOps["abab_vvoo_448"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_342 = dot(sharedOps["bbbb_vvoo_449"]("a,b,i,j"), m2_bbbb_oovv("i,j,b,a"));
        if (include_u1_)
            scalar_343 = dot(sharedOps["aaaa_vvoo_431"]("b,a,i,j"), l2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_344 = dot(sharedOps["aaaa_vvoo_444"]("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_345 = dot(sharedOps["aaaa_vvoo_282"]("b,a,j,i"), m2_aaaa_oovv("i,j,b,a"));
        scalar_346 = dot(sharedOps["aaaa_vvoo_283"]("b,a,i,j"), l2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_347 = dot(sharedOps["aaaa_vvoo_285"]("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        scalar_348 = dot(sharedOps["bb_vo_286"]("a,i"), l1_bb_ov("i,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_349 = dot(sharedOps["aaaa_vvoo_288"]("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_350 = dot(sharedOps["bbbb_vvoo_289"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m1_ && include_u1_ && include_u2_)
            scalar_351 = dot(sharedOps["aa_vo_290"]("a,i"), m1_aa_ov("i,a"));
        scalar_352 = dot(sharedOps["aa_vo_291"]("a,i"), l1_aa_ov("i,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_353 = dot(sharedOps["aaaa_vvoo_191"]("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m1_ && include_u1_ && include_u2_)
            scalar_354 = dot(sharedOps["bb_vo_297"]("a,i"), m1_bb_ov("i,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_355 = dot(sharedOps["aaaa_vvoo_195"]("b,a,j,i"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m1_ && include_u1_)
            scalar_356 = dot(sharedOps["aa_vo_318"]("a,i"), m1_aa_ov("i,a"));
        if (include_u2_)
            scalar_357 = dot(sharedOps["bb_vo_319"]("b,j"), l1_bb_ov("j,b"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_358 = dot(sharedOps["bbbb_vvoo_194"]("b,a,i,j"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m1_ && include_u1_)
            scalar_359 = dot(sharedOps["aa_vo_320"]("a,i"), m1_aa_ov("i,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_360 = dot(sharedOps["bbbb_vvoo_193"]("b,a,i,j"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_361 = dot(sharedOps["aaaa_vvoo_192"]("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_362 = dot(sharedOps["bbbb_vvoo_293"]("b,a,i,j"), m2_bbbb_oovv("i,j,b,a"));
        if (include_u1_)
            scalar_363 = dot(l2_aaaa_oovv("i,j,b,a"), sharedOps["aaaa_vvoo_494"]("a,b,i,j"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_364 = dot(m2_abab_oovv("j,i,b,a"), sharedOps["abab_vvoo_495"]("b,a,j,i"));
        if (include_u1_)
            scalar_365 = dot(l2_abab_oovv("i,j,a,b"), sharedOps["abab_vvoo_496"]("a,b,i,j"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_366 = dot(m2_abab_oovv("i,j,a,b"), sharedOps["abab_vvoo_487"]("a,b,i,j"));
        if (include_u1_)
            scalar_367 = dot(l2_bbbb_oovv("i,j,b,a"), sharedOps["bbbb_vvoo_498"]("a,b,i,j"));
        scalar_368 = dot(sharedOps["abab_vvoo_79"]("b,a,j,i"), l2_abab_oovv("j,i,b,a"));
        if (include_m1_ && include_u1_)
            scalar_369 = dot(sharedOps["bb_vo_501"]("a,i"), m1_bb_ov("i,a"));
        if (include_m1_ && include_u1_)
            scalar_370 = dot(sharedOps["aa_vo_502"]("a,i"), m1_aa_ov("i,a"));
        if (include_m1_ && include_u1_)
            scalar_371 = dot(sharedOps["bb_vo_504"]("a,i"), m1_bb_ov("i,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_372 = dot(sharedOps["abab_vvoo_463"]("a,b,j,i"), m2_abab_oovv("j,i,a,b"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_373 = dot(sharedOps["bbbb_vvoo_452"]("b,a,i,j"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_374 = dot(sharedOps["abab_vvoo_454"]("b,a,i,j"), m2_abab_oovv("i,j,b,a"));
        if (include_m2_ && include_u2_)
            scalar_375 = dot(sharedOps["bbbb_vvoo_455"]("b,a,j,i"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_376 = dot(sharedOps["bbbb_vvoo_456"]("b,a,i,j"), m2_bbbb_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_377 = dot(sharedOps["abab_vvoo_457"]("a,b,i,j"), m2_abab_oovv("i,j,a,b"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_378 = dot(sharedOps["abab_vvoo_458"]("b,a,j,i"), m2_abab_oovv("j,i,b,a"));
        if (include_m1_ && include_u1_)
            scalar_379 = dot(sharedOps["aa_vo_503"]("a,i"), m1_aa_ov("i,a"));
        if (include_m2_ && include_u2_)
            scalar_380 = dot(sharedOps["bbbb_vvoo_27"]("b,a,i,j"), m2_bbbb_oovv("i,j,b,a"));
        if (include_u1_)
            scalar_381 = dot(sharedOps["bbbb_vvoo_430"]("a,b,i,j"), l2_bbbb_oovv("i,j,b,a"));
        scalar_382 = dot(sharedOps["aaaa_vvoo_26"]("b,a,i,j"), l2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_383 = dot(sharedOps["aaaa_vvoo_451"]("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        if (include_u1_)
            scalar_384 = dot(l2_aaaa_oovv("i,j,b,a"), sharedOps["aaaa_vvoo_493"]("a,b,i,j"));
        if (include_m2_ && include_u2_)
            scalar_385 = dot(sharedOps["aaaa_vvoo_24"]("b,a,i,j"), m2_aaaa_oovv("i,j,b,a"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_386 = dot(m2_bbbb_oovv("i,j,b,a"), sharedOps["bbbb_vvoo_497"]("a,b,i,j"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_387 = dot(m2_aaaa_oovv("i,j,b,a"), sharedOps["aaaa_vvoo_488"]("a,b,i,j"));
        if (include_u1_)
            scalar_388 = dot(l2_abab_oovv("j,i,b,a"), sharedOps["abab_vvoo_489"]("b,a,j,i"));
        if (include_m2_ && include_u1_ && include_u2_)
            scalar_389 = dot(m2_bbbb_oovv("i,j,b,a"), sharedOps["bbbb_vvoo_490"]("a,b,i,j"));
        if (include_u1_)
            scalar_390 = dot(l2_abab_oovv("j,i,b,a"), sharedOps["abab_vvoo_491"]("b,a,j,i"));
        if (include_u1_)
            scalar_391 = dot(l2_bbbb_oovv("i,j,b,a"), sharedOps["bbbb_vvoo_492"]("a,b,i,j"));
        auto tempOps = map<std::string, TA::TArrayD>();



/// ****** p†q ****** ///



        if (include_u2_) {

            // tempOps["aa_vv_0"] += 1.000000 m2_xabab_Loovv("I,i,j,b,a") t2_abab_vvoo("c,a,i,j")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempOps["aa_vv_0"]("b,c") = m2_abab_oovv("i,j,b,a") * t2_abab_vvoo("c,a,i,j");

            // rm2_abab_vvoo += -0.500000 <m,n||b,f>_abab t2_abab(b,a,i,j) m2_abab(i,j,e,a)
            // +                -0.500000 <m,n||b,f>_abab t2_abab(b,a,j,i) m2_abab(j,i,e,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= tempOps["aa_vv_0"]("e,b") * V_blks_["abab_oovv"]("m,n,b,f");

            // tempPerm_aaaa_vvoo += -0.500000 P(e,f) <m,n||e,b>_aaaa t2_abab(b,a,i,j) m2_abab(i,j,f,a)
            // +                -0.500000 P(e,f) <m,n||e,b>_aaaa t2_abab(b,a,j,i) m2_abab(j,i,f,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempOps["aa_vv_0"]("f,b") * V_blks_["aaaa_oovv"]("m,n,e,b");
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += 0.500000 <b,m||c,e>_abab t2_abab(c,a,i,j) m2_abab(i,j,b,a)
            // +            0.500000 <b,m||c,e>_abab t2_abab(c,a,j,i) m2_abab(j,i,b,a)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += tempOps["aa_vv_0"]("b,c") * V_blks_["abab_vovv"]("b,m,c,e");

            // rm1_aa_vo += 0.500000 <m,b||e,c>_aaaa t2_abab(c,a,i,j) m2_abab(i,j,b,a)
            // +            0.500000 <m,b||e,c>_aaaa t2_abab(c,a,j,i) m2_abab(j,i,b,a)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= tempOps["aa_vv_0"]("b,c") * V_blks_["aaaa_vovv"]("b,m,e,c");
        }
        tempOps["aa_vv_0"] = TArrayD();

        if (include_u2_) {

            // tempOps["bb_vv_1"] += 1.000000 m2_xbbbb_Loovv("I,i,j,b,a") t2_bbbb_vvoo("c,a,i,j")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempOps["bb_vv_1"]("b,c") = 0.500000 * m2_bbbb_oovv("i,j,b,a") * t2_bbbb_vvoo("c,a,i,j");

            // tempPerm_bbbb_vvoo += -0.500000 P(e,f) <m,n||e,b>_bbbb t2_bbbb(b,a,i,j) m2_bbbb(i,j,f,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempOps["bb_vv_1"]("f,b") * V_blks_["bbbb_oovv"]("m,n,e,b");

            // rm2_abab_vvoo += -0.500000 <m,n||e,b>_abab t2_bbbb(b,a,i,j) m2_bbbb(i,j,f,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= tempOps["bb_vv_1"]("f,b") * V_blks_["abab_oovv"]("m,n,e,b");
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += 0.500000 <m,b||e,c>_bbbb t2_bbbb(c,a,i,j) m2_bbbb(i,j,b,a)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= tempOps["bb_vv_1"]("b,c") * V_blks_["bbbb_vovv"]("b,m,e,c");

            // rm1_aa_vo += 0.500000 <m,b||e,c>_abab t2_bbbb(c,a,i,j) m2_bbbb(i,j,b,a)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= tempOps["bb_vv_1"]("b,c") * V_blks_["baab_vovv"]("b,m,e,c");
        }
        tempOps["bb_vv_1"] = TArrayD();

        if (include_u2_) {

            // tempOps["bb_vv_2"] += 1.000000 t2_abab_vvoo("a,c,i,j") m2_xabab_Loovv("I,i,j,a,b")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempOps["bb_vv_2"]("c,b") = t2_abab_vvoo("a,c,i,j") * m2_abab_oovv("i,j,a,b");

            // tempPerm_bbbb_vvoo += -0.500000 P(e,f) <m,n||e,b>_bbbb t2_abab(a,b,i,j) m2_abab(i,j,a,f)
            // +                -0.500000 P(e,f) <m,n||e,b>_bbbb t2_abab(a,b,j,i) m2_abab(j,i,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= tempOps["bb_vv_2"]("b,f") * V_blks_["bbbb_oovv"]("m,n,e,b");

            // rm2_abab_vvoo += -0.500000 <m,n||e,b>_abab t2_abab(a,b,i,j) m2_abab(i,j,a,f)
            // +                -0.500000 <m,n||e,b>_abab t2_abab(a,b,j,i) m2_abab(j,i,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= tempOps["bb_vv_2"]("b,f") * V_blks_["abab_oovv"]("m,n,e,b");
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += 0.500000 <m,b||e,c>_bbbb t2_abab(a,c,i,j) m2_abab(i,j,a,b)
            // +            0.500000 <m,b||e,c>_bbbb t2_abab(a,c,j,i) m2_abab(j,i,a,b)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= tempOps["bb_vv_2"]("c,b") * V_blks_["bbbb_vovv"]("b,m,e,c");

            // rm1_aa_vo += 0.500000 <m,b||e,c>_abab t2_abab(a,c,i,j) m2_abab(i,j,a,b)
            // +            0.500000 <m,b||e,c>_abab t2_abab(a,c,j,i) m2_abab(j,i,a,b)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= tempOps["bb_vv_2"]("c,b") * V_blks_["baab_vovv"]("b,m,e,c");
        }
        tempOps["bb_vv_2"] = TArrayD();

        if (include_u2_) {

            // tempOps["aa_vv_3"] += 1.000000 m2_xaaaa_Loovv("I,i,j,b,a") t2_aaaa_vvoo("c,a,i,j")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempOps["aa_vv_3"]("b,c") = 0.500000 * m2_aaaa_oovv("i,j,b,a") * t2_aaaa_vvoo("c,a,i,j");

            // rm2_abab_vvoo += -0.500000 <m,n||b,f>_abab t2_aaaa(b,a,i,j) m2_aaaa(i,j,e,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= tempOps["aa_vv_3"]("e,b") * V_blks_["abab_oovv"]("m,n,b,f");

            // tempPerm_aaaa_vvoo += -0.500000 P(e,f) <m,n||e,b>_aaaa t2_aaaa(b,a,i,j) m2_aaaa(i,j,f,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= tempOps["aa_vv_3"]("f,b") * V_blks_["aaaa_oovv"]("m,n,e,b");
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += 0.500000 <b,m||c,e>_abab t2_aaaa(c,a,i,j) m2_aaaa(i,j,b,a)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += tempOps["aa_vv_3"]("b,c") * V_blks_["abab_vovv"]("b,m,c,e");

            // rm1_aa_vo += 0.500000 <m,b||e,c>_aaaa t2_aaaa(c,a,i,j) m2_aaaa(i,j,b,a)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= tempOps["aa_vv_3"]("b,c") * V_blks_["aaaa_vovv"]("b,m,e,c");
        }
        tempOps["aa_vv_3"] = TArrayD();

        if (include_u2_) {

            // tempOps["bb_vv_4"] += 1.000000 m2_xbbbb_Loovv("I,i,j,b,a") u2_bbbb_vvoo("c,a,i,j")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempOps["bb_vv_4"]("b,c") = 0.500000 * m2_bbbb_oovv("i,j,b,a") * u2_bbbb_vvoo("c,a,i,j");

            // rm2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rm2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += -0.500000 P(e,f) <m,n||e,b>_bbbb u2_bbbb(b,a,i,j) m2_bbbb(i,j,f,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempOps["bb_vv_4"]("f,b") * V_blks_["bbbb_oovv"]("m,n,e,b");

            // rl2_abab_vvoo += -0.500000 <m,n||e,b>_abab u2_bbbb(b,a,i,j) m2_bbbb(i,j,f,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["bb_vv_4"]("f,b") * V_blks_["abab_oovv"]("m,n,e,b");

            // rl1_bb_vo += 0.500000 <m,b||e,c>_bbbb u2_bbbb(c,a,i,j) m2_bbbb(i,j,b,a)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bb_vv_4"]("b,c") * V_blks_["bbbb_vovv"]("b,m,e,c");

            // rl1_aa_vo += 0.500000 <m,b||e,c>_abab u2_bbbb(c,a,i,j) m2_bbbb(i,j,b,a)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["bb_vv_4"]("b,c") * V_blks_["baab_vovv"]("b,m,e,c");
        }
        tempOps["bb_vv_4"] = TArrayD();

        {

            // tempOps["aa_vv_5"] += 1.000000 l2_xabab_Loovv("I,i,j,b,a") t2_abab_vvoo("c,a,i,j")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempOps["aa_vv_5"]("b,c") = l2_abab_oovv("i,j,b,a") * t2_abab_vvoo("c,a,i,j");

            // rl2_abab_vvoo += -0.500000 <m,n||b,f>_abab l2_abab(i,j,e,a) t2_abab(b,a,i,j)
            // +                -0.500000 <m,n||b,f>_abab l2_abab(j,i,e,a) t2_abab(b,a,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["aa_vv_5"]("e,b") * V_blks_["abab_oovv"]("m,n,b,f");
        }

        if (include_u2_) {

            // rm2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rm2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");

            // tempPerm_aaaa_vvoo += -0.500000 P(e,f) <m,n||e,b>_aaaa l2_abab(i,j,f,a) t2_abab(b,a,i,j)
            // +                -0.500000 P(e,f) <m,n||e,b>_aaaa l2_abab(j,i,f,a) t2_abab(b,a,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempOps["aa_vv_5"]("f,b") * V_blks_["aaaa_oovv"]("m,n,e,b");

            // rl1_bb_vo += 0.500000 <b,m||c,e>_abab l2_abab(i,j,b,a) t2_abab(c,a,i,j)
            // +            0.500000 <b,m||c,e>_abab l2_abab(j,i,b,a) t2_abab(c,a,j,i)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["aa_vv_5"]("b,c") * V_blks_["abab_vovv"]("b,m,c,e");

            // rl1_aa_vo += 0.500000 <m,b||e,c>_aaaa l2_abab(i,j,b,a) t2_abab(c,a,i,j)
            // +            0.500000 <m,b||e,c>_aaaa l2_abab(j,i,b,a) t2_abab(c,a,j,i)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aa_vv_5"]("b,c") * V_blks_["aaaa_vovv"]("b,m,e,c");
            tempOps["aa_vv_5"] = TArrayD();
        }

        if (include_u2_) {

            // tempOps["bb_vv_6"] += 1.000000 u2_abab_vvoo("a,c,i,j") m2_xabab_Loovv("I,i,j,a,b")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempOps["bb_vv_6"]("c,b") = u2_abab_vvoo("a,c,i,j") * m2_abab_oovv("i,j,a,b");

            // tempPerm_bbbb_vvoo += -0.500000 P(e,f) <m,n||e,b>_bbbb u2_abab(a,b,i,j) m2_abab(i,j,a,f)
            // +                -0.500000 P(e,f) <m,n||e,b>_bbbb u2_abab(a,b,j,i) m2_abab(j,i,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= tempOps["bb_vv_6"]("b,f") * V_blks_["bbbb_oovv"]("m,n,e,b");

            // rl2_abab_vvoo += -0.500000 <m,n||e,b>_abab u2_abab(a,b,i,j) m2_abab(i,j,a,f)
            // +                -0.500000 <m,n||e,b>_abab u2_abab(a,b,j,i) m2_abab(j,i,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["bb_vv_6"]("b,f") * V_blks_["abab_oovv"]("m,n,e,b");

            // rl1_bb_vo += 0.500000 <m,b||e,c>_bbbb u2_abab(a,c,i,j) m2_abab(i,j,a,b)
            // +            0.500000 <m,b||e,c>_bbbb u2_abab(a,c,j,i) m2_abab(j,i,a,b)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bb_vv_6"]("c,b") * V_blks_["bbbb_vovv"]("b,m,e,c");

            // rl1_aa_vo += 0.500000 <m,b||e,c>_abab u2_abab(a,c,i,j) m2_abab(i,j,a,b)
            // +            0.500000 <m,b||e,c>_abab u2_abab(a,c,j,i) m2_abab(j,i,a,b)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["bb_vv_6"]("c,b") * V_blks_["baab_vovv"]("b,m,e,c");
        }
        tempOps["bb_vv_6"] = TArrayD();

        {

            // tempOps["aa_vv_7"] += 1.000000 l2_xaaaa_Loovv("I,i,j,b,a") t2_aaaa_vvoo("c,a,i,j")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempOps["aa_vv_7"]("b,c") = 0.500000 * l2_aaaa_oovv("i,j,b,a") * t2_aaaa_vvoo("c,a,i,j");

            // rl2_abab_vvoo += -0.500000 <m,n||b,f>_abab l2_aaaa(i,j,e,a) t2_aaaa(b,a,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["aa_vv_7"]("e,b") * V_blks_["abab_oovv"]("m,n,b,f");

            // tempPerm_aaaa_vvoo += -0.500000 P(e,f) <m,n||e,b>_aaaa l2_aaaa(i,j,f,a) t2_aaaa(b,a,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= tempOps["aa_vv_7"]("f,b") * V_blks_["aaaa_oovv"]("m,n,e,b");

            // rl1_bb_vo += 0.500000 <b,m||c,e>_abab l2_aaaa(i,j,b,a) t2_aaaa(c,a,i,j)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["aa_vv_7"]("b,c") * V_blks_["abab_vovv"]("b,m,c,e");

            // rl1_aa_vo += 0.500000 <m,b||e,c>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(c,a,i,j)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aa_vv_7"]("b,c") * V_blks_["aaaa_vovv"]("b,m,e,c");
            tempOps["aa_vv_7"] = TArrayD();

            // tempOps["bb_vv_8"] += 1.000000 l2_xbbbb_Loovv("I,i,j,b,a") t2_bbbb_vvoo("c,a,i,j")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempOps["bb_vv_8"]("b,c") = 0.500000 * l2_bbbb_oovv("i,j,b,a") * t2_bbbb_vvoo("c,a,i,j");

            // tempPerm_bbbb_vvoo += -0.500000 P(e,f) <m,n||e,b>_bbbb l2_bbbb(i,j,f,a) t2_bbbb(b,a,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= tempOps["bb_vv_8"]("f,b") * V_blks_["bbbb_oovv"]("m,n,e,b");

            // rl2_abab_vvoo += -0.500000 <m,n||e,b>_abab l2_bbbb(i,j,f,a) t2_bbbb(b,a,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["bb_vv_8"]("f,b") * V_blks_["abab_oovv"]("m,n,e,b");

            // rl1_bb_vo += 0.500000 <m,b||e,c>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(c,a,i,j)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bb_vv_8"]("b,c") * V_blks_["bbbb_vovv"]("b,m,e,c");

            // rl1_aa_vo += 0.500000 <m,b||e,c>_abab l2_bbbb(i,j,b,a) t2_bbbb(c,a,i,j)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["bb_vv_8"]("b,c") * V_blks_["baab_vovv"]("b,m,e,c");
            tempOps["bb_vv_8"] = TArrayD();

            // tempOps["bb_vv_9"] += 1.000000 l2_xabab_Loovv("I,i,j,a,b") t2_abab_vvoo("a,c,i,j")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempOps["bb_vv_9"]("b,c") = l2_abab_oovv("i,j,a,b") * t2_abab_vvoo("a,c,i,j");

            // tempPerm_bbbb_vvoo += -0.500000 P(e,f) <m,n||e,b>_bbbb l2_abab(i,j,a,f) t2_abab(a,b,i,j)
            // +                -0.500000 P(e,f) <m,n||e,b>_bbbb l2_abab(j,i,a,f) t2_abab(a,b,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= tempOps["bb_vv_9"]("f,b") * V_blks_["bbbb_oovv"]("m,n,e,b");

            // rl2_abab_vvoo += -0.500000 <m,n||e,b>_abab l2_abab(i,j,a,f) t2_abab(a,b,i,j)
            // +                -0.500000 <m,n||e,b>_abab l2_abab(j,i,a,f) t2_abab(a,b,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["bb_vv_9"]("f,b") * V_blks_["abab_oovv"]("m,n,e,b");

            // rl1_bb_vo += 0.500000 <m,b||e,c>_bbbb l2_abab(i,j,a,b) t2_abab(a,c,i,j)
            // +            0.500000 <m,b||e,c>_bbbb l2_abab(j,i,a,b) t2_abab(a,c,j,i)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bb_vv_9"]("b,c") * V_blks_["bbbb_vovv"]("b,m,e,c");

            // rl1_aa_vo += 0.500000 <m,b||e,c>_abab l2_abab(i,j,a,b) t2_abab(a,c,i,j)
            // +            0.500000 <m,b||e,c>_abab l2_abab(j,i,a,b) t2_abab(a,c,j,i)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["bb_vv_9"]("b,c") * V_blks_["baab_vovv"]("b,m,e,c");
            tempOps["bb_vv_9"] = TArrayD();
        }

        if (include_u2_) {

            // tempOps["aa_vv_10"] += 1.000000 m2_xabab_Loovv("I,i,j,b,a") u2_abab_vvoo("c,a,i,j")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempOps["aa_vv_10"]("b,c") = m2_abab_oovv("i,j,b,a") * u2_abab_vvoo("c,a,i,j");

            // rl2_abab_vvoo += -0.500000 <m,n||b,f>_abab u2_abab(b,a,i,j) m2_abab(i,j,e,a)
            // +                -0.500000 <m,n||b,f>_abab u2_abab(b,a,j,i) m2_abab(j,i,e,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["aa_vv_10"]("e,b") * V_blks_["abab_oovv"]("m,n,b,f");

            // tempPerm_aaaa_vvoo += -0.500000 P(e,f) <m,n||e,b>_aaaa u2_abab(b,a,i,j) m2_abab(i,j,f,a)
            // +                -0.500000 P(e,f) <m,n||e,b>_aaaa u2_abab(b,a,j,i) m2_abab(j,i,f,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= tempOps["aa_vv_10"]("f,b") * V_blks_["aaaa_oovv"]("m,n,e,b");

            // rl1_bb_vo += 0.500000 <b,m||c,e>_abab u2_abab(c,a,i,j) m2_abab(i,j,b,a)
            // +            0.500000 <b,m||c,e>_abab u2_abab(c,a,j,i) m2_abab(j,i,b,a)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["aa_vv_10"]("b,c") * V_blks_["abab_vovv"]("b,m,c,e");

            // rl1_aa_vo += 0.500000 <m,b||e,c>_aaaa u2_abab(c,a,i,j) m2_abab(i,j,b,a)
            // +            0.500000 <m,b||e,c>_aaaa u2_abab(c,a,j,i) m2_abab(j,i,b,a)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aa_vv_10"]("b,c") * V_blks_["aaaa_vovv"]("b,m,e,c");

            // tempOps["aa_vv_11"] += 1.000000 m2_xaaaa_Loovv("I,i,j,b,a") u2_aaaa_vvoo("c,a,i,j")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempOps["aa_vv_11"]("b,c") = 0.500000 * m2_aaaa_oovv("i,j,b,a") * u2_aaaa_vvoo("c,a,i,j");

            // rl2_abab_vvoo += -0.500000 <m,n||b,f>_abab u2_aaaa(b,a,i,j) m2_aaaa(i,j,e,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["aa_vv_11"]("e,b") * V_blks_["abab_oovv"]("m,n,b,f");

            // tempPerm_aaaa_vvoo += -0.500000 P(e,f) <m,n||e,b>_aaaa u2_aaaa(b,a,i,j) m2_aaaa(i,j,f,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= tempOps["aa_vv_11"]("f,b") * V_blks_["aaaa_oovv"]("m,n,e,b");

            // rl1_bb_vo += 0.500000 <b,m||c,e>_abab u2_aaaa(c,a,i,j) m2_aaaa(i,j,b,a)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["aa_vv_11"]("b,c") * V_blks_["abab_vovv"]("b,m,c,e");

            // rl1_aa_vo += 0.500000 <m,b||e,c>_aaaa u2_aaaa(c,a,i,j) m2_aaaa(i,j,b,a)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aa_vv_11"]("b,c") * V_blks_["aaaa_vovv"]("b,m,e,c");
        }
        tempOps["aa_vv_10"] = TArrayD();
        tempOps["aa_vv_11"] = TArrayD();

        {

            // tempOps["bbbb_vvoo_12"] += 1.000000 l2_xbbbb_Loovv("I,m,n,a,f") dp_bb_vv("a,e")
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempOps["bbbb_vvoo_12"]("f,e,m,n") = l2_bbbb_oovv("m,n,a,f") * dp_bb_vv("a,e");

            // rl2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");

            // tempPerm_bbbb_vvoo += -1.000000 P(e,f) dp_bb(a,e) l2_bbbb(m,n,a,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempOps["bbbb_vvoo_12"]("f,e,m,n");
        }

        if (include_u2_) {

            // rm2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rm2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u0_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(e,f) dp_bb(a,e) l2_bbbb(m,n,a,f) u0
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempOps["bbbb_vvoo_12"]("f,e,m,n") * u0;
        }
        tempOps["bbbb_vvoo_12"] = TArrayD();

        if (include_u2_) {

            // tempOps["bbbb_vvoo_13"] += 1.000000 m2_xbbbb_Loovv("I,m,n,a,f") dp_bb_vv("a,e")
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempOps["bbbb_vvoo_13"]("f,e,m,n") = m2_bbbb_oovv("m,n,a,f") * dp_bb_vv("a,e");
        }

        {

            // rl2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u0_ && include_u2_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(e,f) dp_bb(a,e) u0 m2_bbbb(m,n,a,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempOps["bbbb_vvoo_13"]("f,e,m,n") * u0;
        }

        if (include_u2_) {

            // rm2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rm2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(e,f) dp_bb(a,e) m2_bbbb(m,n,a,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempOps["bbbb_vvoo_13"]("f,e,m,n");
        }
        tempOps["bbbb_vvoo_13"] = TArrayD();

        {

            // tempOps["aaaa_vvoo_14"] += 1.000000 l2_xaaaa_Loovv("I,m,n,a,f") dp_aa_vv("a,e")
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempOps["aaaa_vvoo_14"]("f,e,m,n") = l2_aaaa_oovv("m,n,a,f") * dp_aa_vv("a,e");

            // rl2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");

            // tempPerm_aaaa_vvoo += -1.000000 P(e,f) dp_aa(a,e) l2_aaaa(m,n,a,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempOps["aaaa_vvoo_14"]("f,e,m,n");
        }

        if (include_u2_) {

            // rm2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rm2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u0_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(e,f) dp_aa(a,e) l2_aaaa(m,n,a,f) u0
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempOps["aaaa_vvoo_14"]("f,e,m,n") * u0;
        }
        tempOps["aaaa_vvoo_14"] = TArrayD();

        if (include_u2_) {

            // tempOps["aaaa_vvoo_15"] += 1.000000 dp_aa_vv("a,e") m2_xaaaa_Loovv("I,m,n,a,f")
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempOps["aaaa_vvoo_15"]("e,f,m,n") = dp_aa_vv("a,e") * m2_aaaa_oovv("m,n,a,f");
        }

        {

            // rl2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u0_ && include_u2_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(e,f) dp_aa(a,e) u0 m2_aaaa(m,n,a,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempOps["aaaa_vvoo_15"]("e,f,m,n") * u0;
        }

        if (include_u2_) {

            // rm2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rm2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(e,f) dp_aa(a,e) m2_aaaa(m,n,a,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempOps["aaaa_vvoo_15"]("e,f,m,n");

            // tempOps["abab_vvoo_16"] += 1.000000 m2_xabab_Loovv("I,m,n,e,a") dp_bb_vv("a,f")
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempOps["abab_vvoo_16"]("e,f,m,n") = m2_abab_oovv("m,n,e,a") * dp_bb_vv("a,f");
        }
        tempOps["aaaa_vvoo_15"] = TArrayD();

        if (include_u0_ && include_u2_) {

            // rm2_abab_vvoo += -1.000000 dp_bb(a,f) u0 m2_abab(m,n,e,a)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_abab_vvoo("e,f,m,n") -= tempOps["abab_vvoo_16"]("e,f,m,n") * u0;
        }

        if (include_u2_) {

            // rl2_abab_vvoo += -1.000000 dp_bb(a,f) m2_abab(m,n,e,a)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["abab_vvoo_16"]("e,f,m,n");
        }
        tempOps["abab_vvoo_16"] = TArrayD();

        {

            // tempOps["abab_vvoo_17"] += 1.000000 dp_bb_vv("a,f") l2_xabab_Loovv("I,m,n,e,a")
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempOps["abab_vvoo_17"]("e,f,m,n") = dp_bb_vv("a,f") * l2_abab_oovv("m,n,e,a");
        }

        if (include_u2_) {

            // rm2_abab_vvoo += -1.000000 dp_bb(a,f) l2_abab(m,n,e,a)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_abab_vvoo("e,f,m,n") -= tempOps["abab_vvoo_17"]("e,f,m,n");
        }

        if (include_u0_) {

            // rl2_abab_vvoo += -1.000000 dp_bb(a,f) l2_abab(m,n,e,a) u0
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["abab_vvoo_17"]("e,f,m,n") * u0;
        }
        tempOps["abab_vvoo_17"] = TArrayD();

        if (include_u2_) {

            // tempOps["abab_vvoo_18"] += 1.000000 dp_aa_vv("a,e") m2_xabab_Loovv("I,m,n,a,f")
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempOps["abab_vvoo_18"]("e,f,m,n") = dp_aa_vv("a,e") * m2_abab_oovv("m,n,a,f");
        }

        if (include_u0_ && include_u2_) {

            // rm2_abab_vvoo += -1.000000 dp_aa(a,e) u0 m2_abab(m,n,a,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_abab_vvoo("e,f,m,n") -= tempOps["abab_vvoo_18"]("e,f,m,n") * u0;
        }

        if (include_u2_) {

            // rl2_abab_vvoo += -1.000000 dp_aa(a,e) m2_abab(m,n,a,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["abab_vvoo_18"]("e,f,m,n");
        }
        tempOps["abab_vvoo_18"] = TArrayD();

        {

            // tempOps["abab_vvoo_19"] += 1.000000 l2_xabab_Loovv("I,m,n,a,f") dp_aa_vv("a,e")
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempOps["abab_vvoo_19"]("e,f,m,n") = l2_abab_oovv("m,n,a,f") * dp_aa_vv("a,e");
        }

        if (include_u2_) {

            // rm2_abab_vvoo += -1.000000 dp_aa(a,e) l2_abab(m,n,a,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_abab_vvoo("e,f,m,n") -= tempOps["abab_vvoo_19"]("e,f,m,n");
        }

        if (include_u0_) {

            // rl2_abab_vvoo += -1.000000 dp_aa(a,e) l2_abab(m,n,a,f) u0
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["abab_vvoo_19"]("e,f,m,n") * u0;
        }
        tempOps["abab_vvoo_19"] = TArrayD();

        if (include_u1_ && include_u2_) {

            // tempOps["baab_vooo_20"] += 1.000000 u1_aa_vo("b,k") m2_xabab_Loovv("I,i,j,b,a")
            // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
            tempOps["baab_vooo_20"]("a,k,i,j") = u1_aa_vo("b,k") * m2_abab_oovv("i,j,b,a");

            // rm2_abab_vvoo += 2.000000 dp_aa(i,e) u1_aa(a,i) m2_abab(m,n,a,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += 2.000000 * tempOps["baab_vooo_20"]("f,i,m,n") * dp_aa_ov("i,e");
        }

        {

            // rl2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <j,n||i,e>_abab u1_aa(a,j) m2_abab(i,m,a,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempOps["baab_vooo_20"]("f,j,i,m") * V_blks_["abba_oovo"]("j,n,e,i");

            // rl2_abab_vvoo += 1.000000 <j,m||e,i>_aaaa u1_aa(a,j) m2_abab(i,n,a,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += tempOps["baab_vooo_20"]("f,j,i,n") * V_blks_["aaaa_oovo"]("j,m,e,i");

            // rl2_abab_vvoo += -1.000000 f_aa(i,e) u1_aa(a,i) m2_abab(m,n,a,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["baab_vooo_20"]("f,i,m,n") * F_blks_["aa_ov"]("i,e");

            // rl2_abab_vvoo += 1.000000 <j,n||e,i>_abab u1_aa(a,j) m2_abab(m,i,a,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += tempOps["baab_vooo_20"]("f,j,m,i") * V_blks_["abab_oovo"]("j,n,e,i");

            // rl2_abab_vvoo += -1.000000 <i,b||e,f>_abab u1_aa(a,i) m2_abab(m,n,a,b)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += tempOps["baab_vooo_20"]("b,i,m,n") * V_blks_["baab_vovv"]("b,i,e,f");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // rl2_abab_vvoo += 1.000000 dp_aa(i,e) u0 u1_aa(a,i) m2_abab(m,n,a,f)
            // flops: o3v2L1: 1, o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 2, o0v0: 1,
            rl2_abab_vvoo("e,f,m,n") += u0 * dp_aa_ov("i,e") * tempOps["baab_vooo_20"]("f,i,m,n");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += 1.000000 <j,k||c,e>_abab t2_abab(c,a,i,k) u1_aa(b,j) m2_abab(i,m,b,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["baab_vooo_20"]("a,j,i,m") * sharedOps["bbaa_vvoo_8"]("e,a,j,i");

            // rl1_bb_vo += -1.000000 <j,b||i,e>_abab u1_aa(a,j) m2_abab(i,m,a,b)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["baab_vooo_20"]("b,j,i,m") * V_blks_["baba_vovo"]("b,j,e,i");

            // rl1_bb_vo += -0.500000 <k,m||c,e>_abab t2_abab(c,a,i,j) u1_aa(b,k) m2_abab(i,j,b,a)
            // +            -0.500000 <k,m||c,e>_abab t2_abab(c,a,j,i) u1_aa(b,k) m2_abab(j,i,b,a)
            // flops: o3v2L1: 1, o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 3,
            rl1_bb_vo("e,m") -= t2_abab_vvoo("c,a,i,j") * tempOps["baab_vooo_20"]("a,k,i,j") * V_blks_["abab_oovv"]("k,m,c,e");

            // rl1_aa_vo += 1.000000 <k,j||e,c>_aaaa t2_abab(c,a,k,i) u1_aa(b,j) m2_abab(m,i,b,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["baab_vooo_20"]("a,j,m,i") * sharedOps["abab_vvoo_4"]("e,a,j,i");

            // rl1_aa_vo += 0.500000 <k,m||e,c>_aaaa t2_abab(c,a,i,j) u1_aa(b,k) m2_abab(i,j,b,a)
            // +            0.500000 <k,m||e,c>_aaaa t2_abab(c,a,j,i) u1_aa(b,k) m2_abab(j,i,b,a)
            // flops: o3v2L1: 1, o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 3,
            rl1_aa_vo("e,m") += t2_abab_vvoo("c,a,i,j") * tempOps["baab_vooo_20"]("a,k,i,j") * V_blks_["aaaa_oovv"]("k,m,e,c");

            // rl1_aa_vo += 1.000000 <j,k||e,c>_abab t2_bbbb(c,a,i,k) u1_aa(b,j) m2_abab(m,i,b,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["baab_vooo_20"]("a,j,m,i") * sharedOps["abab_vvoo_7"]("e,a,j,i");

            // rl1_aa_vo += -1.000000 <j,b||e,i>_abab u1_aa(a,j) m2_abab(m,i,a,b)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["baab_vooo_20"]("b,j,m,i") * V_blks_["baab_vovo"]("b,j,e,i");

            // tempOps["aabb_vooo_21"] += 1.000000 m2_xabab_Loovv("I,i,j,a,b") u1_bb_vo("b,k")
            // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
            tempOps["aabb_vooo_21"]("a,i,j,k") = m2_abab_oovv("i,j,a,b") * u1_bb_vo("b,k");

            // rm2_abab_vvoo += 2.000000 dp_bb(i,f) u1_bb(a,i) m2_abab(m,n,e,a)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += 2.000000 * tempOps["aabb_vooo_21"]("e,m,n,i") * dp_bb_ov("i,f");

            // rl2_abab_vvoo += -1.000000 f_bb(i,f) u1_bb(a,i) m2_abab(m,n,e,a)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["aabb_vooo_21"]("e,m,n,i") * F_blks_["bb_ov"]("i,f");

            // rl2_abab_vvoo += 1.000000 <m,j||i,f>_abab u1_bb(a,j) m2_abab(i,n,e,a)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["aabb_vooo_21"]("e,i,n,j") * V_blks_["abba_oovo"]("m,j,f,i");

            // rl2_abab_vvoo += -1.000000 <b,i||e,f>_abab u1_bb(a,i) m2_abab(m,n,b,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["aabb_vooo_21"]("b,m,n,i") * V_blks_["abab_vovv"]("b,i,e,f");
        }
        tempOps["baab_vooo_20"] = TArrayD();

        if (include_u0_ && include_u1_ && include_u2_) {

            // rl2_abab_vvoo += 1.000000 dp_bb(i,f) u0 u1_bb(a,i) m2_abab(m,n,e,a)
            // flops: o3v2L1: 1, o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 2, o0v0: 1,
            rl2_abab_vvoo("e,f,m,n") += u0 * dp_bb_ov("i,f") * tempOps["aabb_vooo_21"]("e,m,n,i");
        }

        if (include_u1_ && include_u2_) {

            // rl2_abab_vvoo += 1.000000 <j,n||f,i>_bbbb u1_bb(a,j) m2_abab(m,i,e,a)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += tempOps["aabb_vooo_21"]("e,m,i,j") * V_blks_["bbbb_oovo"]("j,n,f,i");
        }

        {

            // rl2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <n,j||e,i>_abab u1_bb(a,j) m2_abab(m,i,f,a)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempOps["aabb_vooo_21"]("f,m,i,j") * V_blks_["abab_oovo"]("n,j,e,i");

            // rl1_bb_vo += -1.000000 <b,j||i,e>_abab u1_bb(a,j) m2_abab(i,m,b,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["aabb_vooo_21"]("b,i,m,j") * V_blks_["abba_vovo"]("b,j,e,i");

            // rl1_bb_vo += 1.000000 <k,j||c,e>_abab t2_aaaa(c,a,i,k) u1_bb(b,j) m2_abab(i,m,a,b)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["aabb_vooo_21"]("a,i,m,j") * sharedOps["abab_vvoo_0"]("a,e,i,j");

            // rl1_bb_vo += 0.500000 <k,m||e,c>_bbbb t2_abab(a,c,i,j) u1_bb(b,k) m2_abab(i,j,a,b)
            // +            0.500000 <k,m||e,c>_bbbb t2_abab(a,c,j,i) u1_bb(b,k) m2_abab(j,i,a,b)
            // flops: o3v2L1: 1, o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 3,
            rl1_bb_vo("e,m") += t2_abab_vvoo("a,c,i,j") * tempOps["aabb_vooo_21"]("a,i,j,k") * V_blks_["bbbb_oovv"]("k,m,e,c");

            // rl1_bb_vo += 1.000000 <k,j||e,c>_bbbb t2_abab(a,c,i,k) u1_bb(b,j) m2_abab(i,m,a,b)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["aabb_vooo_21"]("a,i,m,j") * sharedOps["abab_vvoo_3"]("a,e,i,j");

            // rl1_aa_vo += -1.000000 <b,j||e,i>_abab u1_bb(a,j) m2_abab(m,i,b,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aabb_vooo_21"]("b,m,i,j") * V_blks_["abab_vovo"]("b,j,e,i");

            // rl1_aa_vo += 1.000000 <k,j||e,c>_abab t2_abab(a,c,k,i) u1_bb(b,j) m2_abab(m,i,a,b)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["aabb_vooo_21"]("a,m,i,j") * sharedOps["aabb_vvoo_9"]("e,a,j,i");

            // rl1_aa_vo += -0.500000 <m,k||e,c>_abab t2_abab(a,c,i,j) u1_bb(b,k) m2_abab(i,j,a,b)
            // +            -0.500000 <m,k||e,c>_abab t2_abab(a,c,j,i) u1_bb(b,k) m2_abab(j,i,a,b)
            // flops: o3v2L1: 1, o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 3,
            rl1_aa_vo("e,m") -= t2_abab_vvoo("a,c,i,j") * tempOps["aabb_vooo_21"]("a,i,j,k") * V_blks_["abab_oovv"]("m,k,e,c");
        }
        tempOps["aabb_vooo_21"] = TArrayD();

        if (include_u2_) {

            // tempOps["aa_oo_22"] += 1.000000 t2_aaaa_vvoo("b,a,i,j") m2_xaaaa_Loovv("I,i,m,b,a")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempOps["aa_oo_22"]("j,m") = 0.500000 * t2_aaaa_vvoo("b,a,i,j") * m2_aaaa_oovv("i,m,b,a");

            // rm2_abab_vvoo += -0.500000 <j,n||e,f>_abab t2_aaaa(b,a,i,j) m2_aaaa(i,m,b,a)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= tempOps["aa_oo_22"]("j,m") * V_blks_["abab_oovv"]("j,n,e,f");
        }

        {

            // rl2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("f,e,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += -0.500000 P(m,n) <j,n||e,f>_aaaa t2_aaaa(b,a,i,j) m2_aaaa(i,m,b,a)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempOps["aa_oo_22"]("j,m") * V_blks_["aaaa_oovv"]("j,n,e,f");
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += -0.500000 <k,m||j,e>_abab t2_aaaa(b,a,i,k) m2_aaaa(i,j,b,a)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += tempOps["aa_oo_22"]("k,j") * V_blks_["abba_oovo"]("k,m,e,j");

            // rm1_aa_vo += -0.500000 f_aa(j,e) t2_aaaa(b,a,i,j) m2_aaaa(i,m,b,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= tempOps["aa_oo_22"]("j,m") * F_blks_["aa_ov"]("j,e");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // rm1_aa_vo += 0.500000 dp_aa(j,e) t2_aaaa(b,a,i,j) u0 m2_aaaa(i,m,b,a)
            // flops: o2v1L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rm1_aa_vo("e,m") += tempOps["aa_oo_22"]("j,m") * dp_aa_ov("j,e") * u0;
        }

        if (include_u1_ && include_u2_) {

            // rm1_aa_vo += 0.500000 <k,m||e,j>_aaaa t2_aaaa(b,a,i,k) m2_aaaa(i,j,b,a)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += tempOps["aa_oo_22"]("k,j") * V_blks_["aaaa_oovo"]("k,m,e,j");

            // rl1_bb_vo += -0.500000 <k,m||c,e>_abab t2_aaaa(b,a,i,k) u1_aa(c,j) m2_aaaa(i,j,b,a)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["aa_oo_22"]("k,j") * sharedOps["baab_vooo_146"]("e,j,k,m");
        }

        if (include_u2_) {

            // rl1_aa_vo += 0.500000 dp_aa(j,e) t2_aaaa(b,a,i,j) m2_aaaa(i,m,b,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["aa_oo_22"]("j,m") * dp_aa_ov("j,e");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += -0.500000 <k,j||e,c>_abab t2_aaaa(b,a,i,k) u1_bb(c,j) m2_aaaa(i,m,b,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aa_oo_22"]("k,m") * sharedOps["aa_vo_37"]("e,k");

            // rl1_aa_vo += 0.500000 <k,m||e,c>_aaaa t2_aaaa(b,a,i,k) u1_aa(c,j) m2_aaaa(i,j,b,a)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["aa_oo_22"]("k,j") * sharedOps["aaaa_vooo_147"]("e,j,k,m");

            // rl1_aa_vo += -0.500000 <k,j||e,c>_aaaa t2_aaaa(b,a,i,k) u1_aa(c,j) m2_aaaa(i,m,b,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aa_oo_22"]("k,m") * sharedOps["aa_vo_38"]("e,k");
        }
        tempOps["aa_oo_22"] = TArrayD();

        if (include_u2_) {

            // tempOps["aa_oo_23"] += 1.000000 m2_xabab_Loovv("I,m,i,b,a") t2_abab_vvoo("b,a,j,i")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempOps["aa_oo_23"]("m,j") = m2_abab_oovv("m,i,b,a") * t2_abab_vvoo("b,a,j,i");

            // rm2_abab_vvoo += -0.500000 <j,n||e,f>_abab t2_abab(b,a,j,i) m2_abab(m,i,b,a)
            // +                -0.500000 <j,n||e,f>_abab t2_abab(a,b,j,i) m2_abab(m,i,a,b)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= tempOps["aa_oo_23"]("m,j") * V_blks_["abab_oovv"]("j,n,e,f");

            // tempPerm_aaaa_vvoo += -0.500000 P(m,n) <j,n||e,f>_aaaa t2_abab(b,a,j,i) m2_abab(m,i,b,a)
            // +                -0.500000 P(m,n) <j,n||e,f>_aaaa t2_abab(a,b,j,i) m2_abab(m,i,a,b)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= tempOps["aa_oo_23"]("m,j") * V_blks_["aaaa_oovv"]("j,n,e,f");
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += -0.500000 <k,m||j,e>_abab t2_abab(b,a,k,i) m2_abab(j,i,b,a)
            // +            -0.500000 <k,m||j,e>_abab t2_abab(a,b,k,i) m2_abab(j,i,a,b)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += tempOps["aa_oo_23"]("j,k") * V_blks_["abba_oovo"]("k,m,e,j");

            // rm1_aa_vo += -0.500000 f_aa(j,e) t2_abab(b,a,j,i) m2_abab(m,i,b,a)
            // +            -0.500000 f_aa(j,e) t2_abab(a,b,j,i) m2_abab(m,i,a,b)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= tempOps["aa_oo_23"]("m,j") * F_blks_["aa_ov"]("j,e");

            // rm1_aa_vo += 0.500000 <k,m||e,j>_aaaa t2_abab(b,a,k,i) m2_abab(j,i,b,a)
            // +            0.500000 <k,m||e,j>_aaaa t2_abab(a,b,k,i) m2_abab(j,i,a,b)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += tempOps["aa_oo_23"]("j,k") * V_blks_["aaaa_oovo"]("k,m,e,j");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // rm1_aa_vo += 0.500000 dp_aa(j,e) t2_abab(b,a,j,i) u0 m2_abab(m,i,b,a)
            // +            0.500000 dp_aa(j,e) t2_abab(a,b,j,i) u0 m2_abab(m,i,a,b)
            // flops: o2v1L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rm1_aa_vo("e,m") += tempOps["aa_oo_23"]("m,j") * dp_aa_ov("j,e") * u0;
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += -0.500000 <k,m||c,e>_abab t2_abab(b,a,k,i) u1_aa(c,j) m2_abab(j,i,b,a)
            // +            -0.500000 <k,m||c,e>_abab t2_abab(a,b,k,i) u1_aa(c,j) m2_abab(j,i,a,b)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["aa_oo_23"]("j,k") * sharedOps["baab_vooo_146"]("e,j,k,m");

            // rl1_aa_vo += -0.500000 <k,j||e,c>_abab t2_abab(b,a,k,i) u1_bb(c,j) m2_abab(m,i,b,a)
            // +            -0.500000 <k,j||e,c>_abab t2_abab(a,b,k,i) u1_bb(c,j) m2_abab(m,i,a,b)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aa_oo_23"]("m,k") * sharedOps["aa_vo_37"]("e,k");
        }

        if (include_u2_) {

            // rl1_aa_vo += 0.500000 dp_aa(j,e) t2_abab(b,a,j,i) m2_abab(m,i,b,a)
            // +            0.500000 dp_aa(j,e) t2_abab(a,b,j,i) m2_abab(m,i,a,b)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["aa_oo_23"]("m,j") * dp_aa_ov("j,e");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += 0.500000 <k,m||e,c>_aaaa t2_abab(b,a,k,i) u1_aa(c,j) m2_abab(j,i,b,a)
            // +            0.500000 <k,m||e,c>_aaaa t2_abab(a,b,k,i) u1_aa(c,j) m2_abab(j,i,a,b)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["aa_oo_23"]("j,k") * sharedOps["aaaa_vooo_147"]("e,j,k,m");

            // rl1_aa_vo += -0.500000 <k,j||e,c>_aaaa t2_abab(b,a,k,i) u1_aa(c,j) m2_abab(m,i,b,a)
            // +            -0.500000 <k,j||e,c>_aaaa t2_abab(a,b,k,i) u1_aa(c,j) m2_abab(m,i,a,b)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aa_oo_23"]("m,k") * sharedOps["aa_vo_38"]("e,k");
        }
        tempOps["aa_oo_23"] = TArrayD();

        if (include_u2_) {

            // tempOps["bb_oo_24"] += 1.000000 t2_abab_vvoo("b,a,i,k") m2_xabab_Loovv("I,i,j,b,a")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempOps["bb_oo_24"]("k,j") = t2_abab_vvoo("b,a,i,k") * m2_abab_oovv("i,j,b,a");
        }

        {

            // rl2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("f,e,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += -0.500000 P(m,n) <j,n||e,f>_bbbb t2_abab(b,a,i,j) m2_abab(i,m,b,a)
            // +                -0.500000 P(m,n) <j,n||e,f>_bbbb t2_abab(a,b,i,j) m2_abab(i,m,a,b)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempOps["bb_oo_24"]("j,m") * V_blks_["bbbb_oovv"]("j,n,e,f");

            // rm2_abab_vvoo += -0.500000 <m,j||e,f>_abab t2_abab(b,a,i,j) m2_abab(i,n,b,a)
            // +                -0.500000 <m,j||e,f>_abab t2_abab(a,b,i,j) m2_abab(i,n,a,b)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= tempOps["bb_oo_24"]("j,n") * V_blks_["abab_oovv"]("m,j,e,f");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // rm1_bb_vo += 0.500000 dp_bb(j,e) t2_abab(b,a,i,j) u0 m2_abab(i,m,b,a)
            // +            0.500000 dp_bb(j,e) t2_abab(a,b,i,j) u0 m2_abab(i,m,a,b)
            // flops: o2v1L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rm1_bb_vo("e,m") += tempOps["bb_oo_24"]("j,m") * dp_bb_ov("j,e") * u0;
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += 0.500000 <k,m||e,j>_bbbb t2_abab(b,a,i,k) m2_abab(i,j,b,a)
            // +            0.500000 <k,m||e,j>_bbbb t2_abab(a,b,i,k) m2_abab(i,j,a,b)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += tempOps["bb_oo_24"]("k,j") * V_blks_["bbbb_oovo"]("k,m,e,j");

            // rm1_bb_vo += -0.500000 f_bb(j,e) t2_abab(b,a,i,j) m2_abab(i,m,b,a)
            // +            -0.500000 f_bb(j,e) t2_abab(a,b,i,j) m2_abab(i,m,a,b)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= tempOps["bb_oo_24"]("j,m") * F_blks_["bb_ov"]("j,e");

            // rm1_aa_vo += -0.500000 <m,k||e,j>_abab t2_abab(b,a,i,k) m2_abab(i,j,b,a)
            // +            -0.500000 <m,k||e,j>_abab t2_abab(a,b,i,k) m2_abab(i,j,a,b)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= tempOps["bb_oo_24"]("k,j") * V_blks_["abab_oovo"]("m,k,e,j");
        }

        if (include_u2_) {

            // rl1_bb_vo += 0.500000 dp_bb(j,e) t2_abab(b,a,i,j) m2_abab(i,m,b,a)
            // +            0.500000 dp_bb(j,e) t2_abab(a,b,i,j) m2_abab(i,m,a,b)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["bb_oo_24"]("j,m") * dp_bb_ov("j,e");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += -0.500000 <j,k||c,e>_abab t2_abab(b,a,i,k) u1_aa(c,j) m2_abab(i,m,b,a)
            // +            -0.500000 <j,k||c,e>_abab t2_abab(a,b,i,k) u1_aa(c,j) m2_abab(i,m,a,b)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bb_oo_24"]("k,m") * sharedOps["bb_vo_36"]("e,k");

            // rl1_bb_vo += -0.500000 <k,j||e,c>_bbbb t2_abab(b,a,i,k) u1_bb(c,j) m2_abab(i,m,b,a)
            // +            -0.500000 <k,j||e,c>_bbbb t2_abab(a,b,i,k) u1_bb(c,j) m2_abab(i,m,a,b)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bb_oo_24"]("k,m") * sharedOps["bb_vo_39"]("e,k");

            // rl1_bb_vo += 0.500000 <k,m||e,c>_bbbb t2_abab(b,a,i,k) u1_bb(c,j) m2_abab(i,j,b,a)
            // +            0.500000 <k,m||e,c>_bbbb t2_abab(a,b,i,k) u1_bb(c,j) m2_abab(i,j,a,b)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["bb_oo_24"]("k,j") * sharedOps["bbbb_vooo_148"]("e,j,k,m");

            // rl1_aa_vo += -0.500000 <m,k||e,c>_abab t2_abab(b,a,i,k) u1_bb(c,j) m2_abab(i,j,b,a)
            // +            -0.500000 <m,k||e,c>_abab t2_abab(a,b,i,k) u1_bb(c,j) m2_abab(i,j,a,b)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["bb_oo_24"]("k,j") * sharedOps["aabb_vooo_145"]("e,m,j,k");
        }
        tempOps["bb_oo_24"] = TArrayD();

        if (include_u2_) {

            // tempOps["bb_oo_25"] += 1.000000 t2_bbbb_vvoo("b,a,i,k") m2_xbbbb_Loovv("I,i,j,b,a")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempOps["bb_oo_25"]("k,j") = 0.500000 * t2_bbbb_vvoo("b,a,i,k") * m2_bbbb_oovv("i,j,b,a");

            // tempPerm_bbbb_vvoo += -0.500000 P(m,n) <j,n||e,f>_bbbb t2_bbbb(b,a,i,j) m2_bbbb(i,m,b,a)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= tempOps["bb_oo_25"]("j,m") * V_blks_["bbbb_oovv"]("j,n,e,f");

            // rm2_abab_vvoo += -0.500000 <m,j||e,f>_abab t2_bbbb(b,a,i,j) m2_bbbb(i,n,b,a)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= tempOps["bb_oo_25"]("j,n") * V_blks_["abab_oovv"]("m,j,e,f");
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += -0.500000 f_bb(j,e) t2_bbbb(b,a,i,j) m2_bbbb(i,m,b,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= tempOps["bb_oo_25"]("j,m") * F_blks_["bb_ov"]("j,e");

            // rm1_bb_vo += 0.500000 <k,m||e,j>_bbbb t2_bbbb(b,a,i,k) m2_bbbb(i,j,b,a)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += tempOps["bb_oo_25"]("k,j") * V_blks_["bbbb_oovo"]("k,m,e,j");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // rm1_bb_vo += 0.500000 dp_bb(j,e) t2_bbbb(b,a,i,j) u0 m2_bbbb(i,m,b,a)
            // flops: o2v1L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rm1_bb_vo("e,m") += tempOps["bb_oo_25"]("j,m") * dp_bb_ov("j,e") * u0;
        }

        if (include_u1_ && include_u2_) {

            // rm1_aa_vo += -0.500000 <m,k||e,j>_abab t2_bbbb(b,a,i,k) m2_bbbb(i,j,b,a)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= tempOps["bb_oo_25"]("k,j") * V_blks_["abab_oovo"]("m,k,e,j");

            // rl1_bb_vo += 0.500000 <k,m||e,c>_bbbb t2_bbbb(b,a,i,k) u1_bb(c,j) m2_bbbb(i,j,b,a)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["bb_oo_25"]("k,j") * sharedOps["bbbb_vooo_148"]("e,j,k,m");

            // rl1_bb_vo += -0.500000 <k,j||e,c>_bbbb t2_bbbb(b,a,i,k) u1_bb(c,j) m2_bbbb(i,m,b,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bb_oo_25"]("k,m") * sharedOps["bb_vo_39"]("e,k");

            // rl1_bb_vo += -0.500000 <j,k||c,e>_abab t2_bbbb(b,a,i,k) u1_aa(c,j) m2_bbbb(i,m,b,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bb_oo_25"]("k,m") * sharedOps["bb_vo_36"]("e,k");
        }

        if (include_u2_) {

            // rl1_bb_vo += 0.500000 dp_bb(j,e) t2_bbbb(b,a,i,j) m2_bbbb(i,m,b,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["bb_oo_25"]("j,m") * dp_bb_ov("j,e");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += -0.500000 <m,k||e,c>_abab t2_bbbb(b,a,i,k) u1_bb(c,j) m2_bbbb(i,j,b,a)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["bb_oo_25"]("k,j") * sharedOps["aabb_vooo_145"]("e,m,j,k");
        }
        tempOps["bb_oo_25"] = TArrayD();

        if (include_u2_) {

            // tempOps["aa_oo_26"] += 1.000000 u2_abab_vvoo("b,a,j,i") m2_xabab_Loovv("I,m,i,b,a")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempOps["aa_oo_26"]("j,m") = u2_abab_vvoo("b,a,j,i") * m2_abab_oovv("m,i,b,a");
        }

        if (include_u1_ && include_u2_) {

            // rm1_aa_vo += 1.000000 dp_aa(j,e) u2_abab(b,a,j,i) m2_abab(m,i,b,a)
            // +            1.000000 dp_aa(j,e) u2_abab(a,b,j,i) m2_abab(m,i,a,b)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += 2.000000 * tempOps["aa_oo_26"]("j,m") * dp_aa_ov("j,e");
        }

        if (include_u2_) {

            // rl2_abab_vvoo += -0.500000 <j,n||e,f>_abab u2_abab(b,a,j,i) m2_abab(m,i,b,a)
            // +                -0.500000 <j,n||e,f>_abab u2_abab(a,b,j,i) m2_abab(m,i,a,b)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["aa_oo_26"]("j,m") * V_blks_["abab_oovv"]("j,n,e,f");

            // rm2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rm2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += -0.500000 P(m,n) <j,n||e,f>_aaaa u2_abab(b,a,j,i) m2_abab(m,i,b,a)
            // +                -0.500000 P(m,n) <j,n||e,f>_aaaa u2_abab(a,b,j,i) m2_abab(m,i,a,b)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempOps["aa_oo_26"]("j,m") * V_blks_["aaaa_oovv"]("j,n,e,f");

            // rl1_bb_vo += -0.500000 <k,m||j,e>_abab u2_abab(b,a,k,i) m2_abab(j,i,b,a)
            // +            -0.500000 <k,m||j,e>_abab u2_abab(a,b,k,i) m2_abab(j,i,a,b)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["aa_oo_26"]("k,j") * V_blks_["abba_oovo"]("k,m,e,j");

            // rl1_aa_vo += 0.500000 <k,m||e,j>_aaaa u2_abab(b,a,k,i) m2_abab(j,i,b,a)
            // +            0.500000 <k,m||e,j>_aaaa u2_abab(a,b,k,i) m2_abab(j,i,a,b)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["aa_oo_26"]("k,j") * V_blks_["aaaa_oovo"]("k,m,e,j");
        }

        if (include_u0_ && include_u2_) {

            // rl1_aa_vo += 0.500000 dp_aa(j,e) u0 u2_abab(b,a,j,i) m2_abab(m,i,b,a)
            // +            0.500000 dp_aa(j,e) u0 u2_abab(a,b,j,i) m2_abab(m,i,a,b)
            // flops: o2v1L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_aa_vo("e,m") += u0 * dp_aa_ov("j,e") * tempOps["aa_oo_26"]("j,m");
        }

        if (include_u2_) {

            // rl1_aa_vo += -0.500000 f_aa(j,e) u2_abab(b,a,j,i) m2_abab(m,i,b,a)
            // +            -0.500000 f_aa(j,e) u2_abab(a,b,j,i) m2_abab(m,i,a,b)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aa_oo_26"]("j,m") * F_blks_["aa_ov"]("j,e");

            // tempOps["bb_oo_27"] += 1.000000 m2_xabab_Loovv("I,i,j,b,a") u2_abab_vvoo("b,a,i,k")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempOps["bb_oo_27"]("j,k") = m2_abab_oovv("i,j,b,a") * u2_abab_vvoo("b,a,i,k");
        }
        tempOps["aa_oo_26"] = TArrayD();

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += 1.000000 dp_bb(j,e) u2_abab(b,a,i,j) m2_abab(i,m,b,a)
            // +            1.000000 dp_bb(j,e) u2_abab(a,b,i,j) m2_abab(i,m,a,b)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += 2.000000 * tempOps["bb_oo_27"]("m,j") * dp_bb_ov("j,e");
        }

        if (include_u2_) {

            // rm2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rm2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += -0.500000 P(m,n) <j,n||e,f>_bbbb u2_abab(b,a,i,j) m2_abab(i,m,b,a)
            // +                -0.500000 P(m,n) <j,n||e,f>_bbbb u2_abab(a,b,i,j) m2_abab(i,m,a,b)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempOps["bb_oo_27"]("m,j") * V_blks_["bbbb_oovv"]("j,n,e,f");

            // rl2_abab_vvoo += -0.500000 <m,j||e,f>_abab u2_abab(b,a,i,j) m2_abab(i,n,b,a)
            // +                -0.500000 <m,j||e,f>_abab u2_abab(a,b,i,j) m2_abab(i,n,a,b)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["bb_oo_27"]("n,j") * V_blks_["abab_oovv"]("m,j,e,f");
        }

        if (include_u0_ && include_u2_) {

            // rl1_bb_vo += 0.500000 dp_bb(j,e) u0 u2_abab(b,a,i,j) m2_abab(i,m,b,a)
            // +            0.500000 dp_bb(j,e) u0 u2_abab(a,b,i,j) m2_abab(i,m,a,b)
            // flops: o2v1L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_bb_vo("e,m") += u0 * dp_bb_ov("j,e") * tempOps["bb_oo_27"]("m,j");
        }

        if (include_u2_) {

            // rl1_bb_vo += 0.500000 <k,m||e,j>_bbbb u2_abab(b,a,i,k) m2_abab(i,j,b,a)
            // +            0.500000 <k,m||e,j>_bbbb u2_abab(a,b,i,k) m2_abab(i,j,a,b)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["bb_oo_27"]("j,k") * V_blks_["bbbb_oovo"]("k,m,e,j");

            // rl1_bb_vo += -0.500000 f_bb(j,e) u2_abab(b,a,i,j) m2_abab(i,m,b,a)
            // +            -0.500000 f_bb(j,e) u2_abab(a,b,i,j) m2_abab(i,m,a,b)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bb_oo_27"]("m,j") * F_blks_["bb_ov"]("j,e");

            // rl1_aa_vo += -0.500000 <m,k||e,j>_abab u2_abab(b,a,i,k) m2_abab(i,j,b,a)
            // +            -0.500000 <m,k||e,j>_abab u2_abab(a,b,i,k) m2_abab(i,j,a,b)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["bb_oo_27"]("j,k") * V_blks_["abab_oovo"]("m,k,e,j");

            // tempOps["aa_oo_28"] += 1.000000 m2_xaaaa_Loovv("I,i,m,b,a") u2_aaaa_vvoo("b,a,i,j")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempOps["aa_oo_28"]("m,j") = 0.500000 * m2_aaaa_oovv("i,m,b,a") * u2_aaaa_vvoo("b,a,i,j");
        }
        tempOps["bb_oo_27"] = TArrayD();

        if (include_u1_ && include_u2_) {

            // rm1_aa_vo += 1.000000 dp_aa(j,e) u2_aaaa(b,a,i,j) m2_aaaa(i,m,b,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += 2.000000 * tempOps["aa_oo_28"]("m,j") * dp_aa_ov("j,e");
        }

        if (include_u2_) {

            // rl2_abab_vvoo += -0.500000 <j,n||e,f>_abab u2_aaaa(b,a,i,j) m2_aaaa(i,m,b,a)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["aa_oo_28"]("m,j") * V_blks_["abab_oovv"]("j,n,e,f");

            // tempPerm_aaaa_vvoo += -0.500000 P(m,n) <j,n||e,f>_aaaa u2_aaaa(b,a,i,j) m2_aaaa(i,m,b,a)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= tempOps["aa_oo_28"]("m,j") * V_blks_["aaaa_oovv"]("j,n,e,f");

            // rl1_bb_vo += -0.500000 <k,m||j,e>_abab u2_aaaa(b,a,i,k) m2_aaaa(i,j,b,a)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["aa_oo_28"]("j,k") * V_blks_["abba_oovo"]("k,m,e,j");

            // rl1_aa_vo += 0.500000 <k,m||e,j>_aaaa u2_aaaa(b,a,i,k) m2_aaaa(i,j,b,a)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["aa_oo_28"]("j,k") * V_blks_["aaaa_oovo"]("k,m,e,j");
        }

        if (include_u0_ && include_u2_) {

            // rl1_aa_vo += 0.500000 dp_aa(j,e) u0 u2_aaaa(b,a,i,j) m2_aaaa(i,m,b,a)
            // flops: o2v1L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_aa_vo("e,m") += u0 * dp_aa_ov("j,e") * tempOps["aa_oo_28"]("m,j");
        }

        if (include_u2_) {

            // rl1_aa_vo += -0.500000 f_aa(j,e) u2_aaaa(b,a,i,j) m2_aaaa(i,m,b,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aa_oo_28"]("m,j") * F_blks_["aa_ov"]("j,e");
        }
        tempOps["aa_oo_28"] = TArrayD();

        {

            // tempOps["bb_oo_29"] += 1.000000 l2_xbbbb_Loovv("I,i,j,b,a") t2_bbbb_vvoo("b,a,i,k")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempOps["bb_oo_29"]("j,k") = 0.500000 * l2_bbbb_oovv("i,j,b,a") * t2_bbbb_vvoo("b,a,i,k");
        }

        if (include_u1_) {

            // rm1_bb_vo += 0.500000 dp_bb(j,e) l2_bbbb(i,m,b,a) t2_bbbb(b,a,i,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += tempOps["bb_oo_29"]("m,j") * dp_bb_ov("j,e");
        }

        {

            // tempPerm_bbbb_vvoo += -0.500000 P(m,n) <j,n||e,f>_bbbb l2_bbbb(i,m,b,a) t2_bbbb(b,a,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= tempOps["bb_oo_29"]("m,j") * V_blks_["bbbb_oovv"]("j,n,e,f");

            // rl2_abab_vvoo += -0.500000 <m,j||e,f>_abab l2_bbbb(i,n,b,a) t2_bbbb(b,a,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["bb_oo_29"]("n,j") * V_blks_["abab_oovv"]("m,j,e,f");

            // rl1_bb_vo += -0.500000 f_bb(j,e) l2_bbbb(i,m,b,a) t2_bbbb(b,a,i,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bb_oo_29"]("m,j") * F_blks_["bb_ov"]("j,e");

            // rl1_bb_vo += 0.500000 <k,m||e,j>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(b,a,i,k)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["bb_oo_29"]("j,k") * V_blks_["bbbb_oovo"]("k,m,e,j");
        }

        if (include_u0_) {

            // rl1_bb_vo += 0.500000 dp_bb(j,e) l2_bbbb(i,m,b,a) t2_bbbb(b,a,i,j) u0
            // flops: o2v1L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_bb_vo("e,m") += tempOps["bb_oo_29"]("m,j") * u0 * dp_bb_ov("j,e");
        }

        {

            // rl1_aa_vo += -0.500000 <m,k||e,j>_abab l2_bbbb(i,j,b,a) t2_bbbb(b,a,i,k)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["bb_oo_29"]("j,k") * V_blks_["abab_oovo"]("m,k,e,j");
            tempOps["bb_oo_29"] = TArrayD();

            // tempOps["bb_oo_30"] += 1.000000 l2_xabab_Loovv("I,i,j,b,a") t2_abab_vvoo("b,a,i,k")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempOps["bb_oo_30"]("j,k") = l2_abab_oovv("i,j,b,a") * t2_abab_vvoo("b,a,i,k");
        }

        if (include_u1_) {

            // rm1_bb_vo += 0.500000 dp_bb(j,e) l2_abab(i,m,b,a) t2_abab(b,a,i,j)
            // +            0.500000 dp_bb(j,e) l2_abab(i,m,a,b) t2_abab(a,b,i,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += tempOps["bb_oo_30"]("m,j") * dp_bb_ov("j,e");
        }

        {

            // tempPerm_bbbb_vvoo += -0.500000 P(m,n) <j,n||e,f>_bbbb l2_abab(i,m,b,a) t2_abab(b,a,i,j)
            // +                -0.500000 P(m,n) <j,n||e,f>_bbbb l2_abab(i,m,a,b) t2_abab(a,b,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= tempOps["bb_oo_30"]("m,j") * V_blks_["bbbb_oovv"]("j,n,e,f");

            // rl2_abab_vvoo += -0.500000 <m,j||e,f>_abab l2_abab(i,n,b,a) t2_abab(b,a,i,j)
            // +                -0.500000 <m,j||e,f>_abab l2_abab(i,n,a,b) t2_abab(a,b,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["bb_oo_30"]("n,j") * V_blks_["abab_oovv"]("m,j,e,f");
        }

        if (include_u0_) {

            // rl1_bb_vo += 0.500000 dp_bb(j,e) l2_abab(i,m,b,a) t2_abab(b,a,i,j) u0
            // +            0.500000 dp_bb(j,e) l2_abab(i,m,a,b) t2_abab(a,b,i,j) u0
            // flops: o2v1L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_bb_vo("e,m") += tempOps["bb_oo_30"]("m,j") * u0 * dp_bb_ov("j,e");
        }

        {

            // rl1_bb_vo += 0.500000 <k,m||e,j>_bbbb l2_abab(i,j,b,a) t2_abab(b,a,i,k)
            // +            0.500000 <k,m||e,j>_bbbb l2_abab(i,j,a,b) t2_abab(a,b,i,k)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["bb_oo_30"]("j,k") * V_blks_["bbbb_oovo"]("k,m,e,j");

            // rl1_bb_vo += -0.500000 f_bb(j,e) l2_abab(i,m,b,a) t2_abab(b,a,i,j)
            // +            -0.500000 f_bb(j,e) l2_abab(i,m,a,b) t2_abab(a,b,i,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bb_oo_30"]("m,j") * F_blks_["bb_ov"]("j,e");

            // rl1_aa_vo += -0.500000 <m,k||e,j>_abab l2_abab(i,j,b,a) t2_abab(b,a,i,k)
            // +            -0.500000 <m,k||e,j>_abab l2_abab(i,j,a,b) t2_abab(a,b,i,k)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["bb_oo_30"]("j,k") * V_blks_["abab_oovo"]("m,k,e,j");
            tempOps["bb_oo_30"] = TArrayD();

            // tempOps["aa_oo_31"] += 1.000000 l2_xaaaa_Loovv("I,i,m,b,a") t2_aaaa_vvoo("b,a,i,j")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempOps["aa_oo_31"]("m,j") = 0.500000 * l2_aaaa_oovv("i,m,b,a") * t2_aaaa_vvoo("b,a,i,j");
        }

        if (include_u1_) {

            // rm1_aa_vo += 0.500000 dp_aa(j,e) l2_aaaa(i,m,b,a) t2_aaaa(b,a,i,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += tempOps["aa_oo_31"]("m,j") * dp_aa_ov("j,e");
        }

        {

            // rl2_abab_vvoo += -0.500000 <j,n||e,f>_abab l2_aaaa(i,m,b,a) t2_aaaa(b,a,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["aa_oo_31"]("m,j") * V_blks_["abab_oovv"]("j,n,e,f");

            // tempPerm_aaaa_vvoo += -0.500000 P(m,n) <j,n||e,f>_aaaa l2_aaaa(i,m,b,a) t2_aaaa(b,a,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= tempOps["aa_oo_31"]("m,j") * V_blks_["aaaa_oovv"]("j,n,e,f");

            // rl1_bb_vo += -0.500000 <k,m||j,e>_abab l2_aaaa(i,j,b,a) t2_aaaa(b,a,i,k)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["aa_oo_31"]("j,k") * V_blks_["abba_oovo"]("k,m,e,j");

            // rl1_aa_vo += -0.500000 f_aa(j,e) l2_aaaa(i,m,b,a) t2_aaaa(b,a,i,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aa_oo_31"]("m,j") * F_blks_["aa_ov"]("j,e");
        }

        if (include_u0_) {

            // rl1_aa_vo += 0.500000 dp_aa(j,e) l2_aaaa(i,m,b,a) t2_aaaa(b,a,i,j) u0
            // flops: o2v1L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_aa_vo("e,m") += tempOps["aa_oo_31"]("m,j") * u0 * dp_aa_ov("j,e");
        }

        {

            // rl1_aa_vo += 0.500000 <k,m||e,j>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(b,a,i,k)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["aa_oo_31"]("j,k") * V_blks_["aaaa_oovo"]("k,m,e,j");
            tempOps["aa_oo_31"] = TArrayD();

            // tempOps["aa_oo_32"] += 1.000000 t2_abab_vvoo("b,a,j,i") l2_xabab_Loovv("I,m,i,b,a")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempOps["aa_oo_32"]("j,m") = t2_abab_vvoo("b,a,j,i") * l2_abab_oovv("m,i,b,a");
        }

        if (include_u1_) {

            // rm1_aa_vo += 0.500000 dp_aa(j,e) l2_abab(m,i,b,a) t2_abab(b,a,j,i)
            // +            0.500000 dp_aa(j,e) l2_abab(m,i,a,b) t2_abab(a,b,j,i)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += tempOps["aa_oo_32"]("j,m") * dp_aa_ov("j,e");
        }

        {

            // rl2_abab_vvoo += -0.500000 <j,n||e,f>_abab l2_abab(m,i,b,a) t2_abab(b,a,j,i)
            // +                -0.500000 <j,n||e,f>_abab l2_abab(m,i,a,b) t2_abab(a,b,j,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["aa_oo_32"]("j,m") * V_blks_["abab_oovv"]("j,n,e,f");

            // tempPerm_aaaa_vvoo += -0.500000 P(m,n) <j,n||e,f>_aaaa l2_abab(m,i,b,a) t2_abab(b,a,j,i)
            // +                -0.500000 P(m,n) <j,n||e,f>_aaaa l2_abab(m,i,a,b) t2_abab(a,b,j,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= tempOps["aa_oo_32"]("j,m") * V_blks_["aaaa_oovv"]("j,n,e,f");

            // rl1_bb_vo += -0.500000 <k,m||j,e>_abab l2_abab(j,i,b,a) t2_abab(b,a,k,i)
            // +            -0.500000 <k,m||j,e>_abab l2_abab(j,i,a,b) t2_abab(a,b,k,i)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["aa_oo_32"]("k,j") * V_blks_["abba_oovo"]("k,m,e,j");

            // rl1_aa_vo += 0.500000 <k,m||e,j>_aaaa l2_abab(j,i,b,a) t2_abab(b,a,k,i)
            // +            0.500000 <k,m||e,j>_aaaa l2_abab(j,i,a,b) t2_abab(a,b,k,i)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["aa_oo_32"]("k,j") * V_blks_["aaaa_oovo"]("k,m,e,j");

            // rl1_aa_vo += -0.500000 f_aa(j,e) l2_abab(m,i,b,a) t2_abab(b,a,j,i)
            // +            -0.500000 f_aa(j,e) l2_abab(m,i,a,b) t2_abab(a,b,j,i)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aa_oo_32"]("j,m") * F_blks_["aa_ov"]("j,e");
        }

        if (include_u0_) {

            // rl1_aa_vo += 0.500000 dp_aa(j,e) l2_abab(m,i,b,a) t2_abab(b,a,j,i) u0
            // +            0.500000 dp_aa(j,e) l2_abab(m,i,a,b) t2_abab(a,b,j,i) u0
            // flops: o2v1L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_aa_vo("e,m") += tempOps["aa_oo_32"]("j,m") * u0 * dp_aa_ov("j,e");
        }
        tempOps["aa_oo_32"] = TArrayD();

        if (include_u2_) {

            // tempOps["bb_oo_33"] += 1.000000 m2_xbbbb_Loovv("I,i,j,b,a") u2_bbbb_vvoo("b,a,i,k")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempOps["bb_oo_33"]("j,k") = 0.500000 * m2_bbbb_oovv("i,j,b,a") * u2_bbbb_vvoo("b,a,i,k");
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += 1.000000 dp_bb(j,e) u2_bbbb(b,a,i,j) m2_bbbb(i,m,b,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += 2.000000 * tempOps["bb_oo_33"]("m,j") * dp_bb_ov("j,e");
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += -0.500000 P(m,n) <j,n||e,f>_bbbb u2_bbbb(b,a,i,j) m2_bbbb(i,m,b,a)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= tempOps["bb_oo_33"]("m,j") * V_blks_["bbbb_oovv"]("j,n,e,f");

            // rl2_abab_vvoo += -0.500000 <m,j||e,f>_abab u2_bbbb(b,a,i,j) m2_bbbb(i,n,b,a)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["bb_oo_33"]("n,j") * V_blks_["abab_oovv"]("m,j,e,f");

            // rl1_bb_vo += -0.500000 f_bb(j,e) u2_bbbb(b,a,i,j) m2_bbbb(i,m,b,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bb_oo_33"]("m,j") * F_blks_["bb_ov"]("j,e");

            // rl1_bb_vo += 0.500000 <k,m||e,j>_bbbb u2_bbbb(b,a,i,k) m2_bbbb(i,j,b,a)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["bb_oo_33"]("j,k") * V_blks_["bbbb_oovo"]("k,m,e,j");
        }

        if (include_u0_ && include_u2_) {

            // rl1_bb_vo += 0.500000 dp_bb(j,e) u0 u2_bbbb(b,a,i,j) m2_bbbb(i,m,b,a)
            // flops: o2v1L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_bb_vo("e,m") += u0 * dp_bb_ov("j,e") * tempOps["bb_oo_33"]("m,j");
        }

        if (include_u2_) {

            // rl1_aa_vo += -0.500000 <m,k||e,j>_abab u2_bbbb(b,a,i,k) m2_bbbb(i,j,b,a)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["bb_oo_33"]("j,k") * V_blks_["abab_oovo"]("m,k,e,j");
        }
        tempOps["bb_oo_33"] = TArrayD();

        if (include_u1_ && include_u2_) {

            // tempOps["bbbb_vooo_34"] += 1.000000 m2_xbbbb_Loovv("I,m,i,b,a") u1_bb_vo("a,j")
            // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
            tempOps["bbbb_vooo_34"]("b,m,i,j") = m2_bbbb_oovv("m,i,b,a") * u1_bb_vo("a,j");
        }

        {

            // rl2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_bbbb_vvoo += -2.000000 P(e,f) dp_bb(i,e) u1_bb(a,i) m2_bbbb(m,n,f,a)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -2.000000 * tempOps["bbbb_vooo_34"]("f,m,n,i") * dp_bb_ov("i,e");
        }

        if (include_u2_) {

            // rm2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rm2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) <j,n||e,i>_bbbb u1_bb(a,j) m2_bbbb(m,i,f,a)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempOps["bbbb_vooo_34"]("f,m,i,j") * V_blks_["bbbb_oovo"]("j,n,e,i");
        }

        {

            // rl2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("f,e,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(e,f) f_bb(i,e) u1_bb(a,i) m2_bbbb(m,n,f,a)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempOps["bbbb_vooo_34"]("f,m,n,i") * F_blks_["bb_ov"]("i,e");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(e,f) dp_bb(i,e) u0 u1_bb(a,i) m2_bbbb(m,n,f,a)
            // flops: o3v2L1: 1, o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 2, o0v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") -= u0 * tempOps["bbbb_vooo_34"]("f,m,n,i") * dp_bb_ov("i,e");
        }

        if (include_u1_ && include_u2_) {

            // rl2_bbbb_vvoo += 1.000000 <i,b||e,f>_bbbb u1_bb(a,i) m2_bbbb(m,n,b,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_bbbb_vvoo("e,f,m,n") -= tempOps["bbbb_vooo_34"]("b,m,n,i") * V_blks_["bbbb_vovv"]("b,i,e,f");

            // rl2_abab_vvoo += -1.000000 <m,j||e,i>_abab u1_bb(a,j) m2_bbbb(n,i,f,a)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["bbbb_vooo_34"]("f,n,i,j") * V_blks_["abab_oovo"]("m,j,e,i");

            // rl1_bb_vo += 1.000000 <j,b||e,i>_bbbb u1_bb(a,j) m2_bbbb(m,i,b,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bbbb_vooo_34"]("b,m,i,j") * V_blks_["bbbb_vovo"]("b,j,e,i");

            // tempOps["aaaa_vooo_35"] += 1.000000 u1_aa_vo("a,j") m2_xaaaa_Loovv("I,m,i,b,a")
            // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
            tempOps["aaaa_vooo_35"]("b,j,m,i") = u1_aa_vo("a,j") * m2_aaaa_oovv("m,i,b,a");
        }
        tempOps["bbbb_vooo_34"] = TArrayD();

        {

            // rl2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += -2.000000 P(e,f) dp_aa(i,e) u1_aa(a,i) m2_aaaa(m,n,f,a)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -2.000000 * tempOps["aaaa_vooo_35"]("f,i,m,n") * dp_aa_ov("i,e");

            // rl2_abab_vvoo += -1.000000 <j,n||i,f>_abab u1_aa(a,j) m2_aaaa(m,i,e,a)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += tempOps["aaaa_vooo_35"]("e,j,m,i") * V_blks_["abba_oovo"]("j,n,f,i");
        }

        if (include_u2_) {

            // rm2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rm2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) <j,n||e,i>_aaaa u1_aa(a,j) m2_aaaa(m,i,f,a)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempOps["aaaa_vooo_35"]("f,j,m,i") * V_blks_["aaaa_oovo"]("j,n,e,i");
        }

        {

            // rl2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("f,e,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(e,f) dp_aa(i,e) u0 u1_aa(a,i) m2_aaaa(m,n,f,a)
            // flops: o3v2L1: 1, o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 2, o0v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * u0 * tempOps["aaaa_vooo_35"]("f,i,m,n") * dp_aa_ov("i,e");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(e,f) f_aa(i,e) u1_aa(a,i) m2_aaaa(m,n,f,a)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += tempOps["aaaa_vooo_35"]("f,i,m,n") * F_blks_["aa_ov"]("i,e");

            // rl2_aaaa_vvoo += 1.000000 <i,b||e,f>_aaaa u1_aa(a,i) m2_aaaa(m,n,b,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_aaaa_vvoo("e,f,m,n") -= tempOps["aaaa_vooo_35"]("b,i,m,n") * V_blks_["aaaa_vovv"]("b,i,e,f");

            // rl1_aa_vo += 1.000000 <j,b||e,i>_aaaa u1_aa(a,j) m2_aaaa(m,i,b,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aaaa_vooo_35"]("b,j,m,i") * V_blks_["aaaa_vovo"]("b,j,e,i");

            // tempOps["aaaa_vooo_36"] += 1.000000 m2_xaaaa_Loovv("I,i,j,b,a") u1_aa_vo("b,k")
            // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
            tempOps["aaaa_vooo_36"]("a,i,j,k") = m2_aaaa_oovv("i,j,b,a") * u1_aa_vo("b,k");

            // rl1_bb_vo += -0.500000 <k,m||c,e>_abab t2_aaaa(c,a,i,j) u1_aa(b,k) m2_aaaa(i,j,b,a)
            // flops: o3v2L1: 1, o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 3,
            rl1_bb_vo("e,m") -= 0.500000 * t2_aaaa_vvoo("c,a,i,j") * tempOps["aaaa_vooo_36"]("a,i,j,k") * V_blks_["abab_oovv"]("k,m,c,e");

            // rl1_aa_vo += 0.500000 <k,m||e,c>_aaaa t2_aaaa(c,a,i,j) u1_aa(b,k) m2_aaaa(i,j,b,a)
            // flops: o3v2L1: 1, o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 3,
            rl1_aa_vo("e,m") += 0.500000 * t2_aaaa_vvoo("c,a,i,j") * tempOps["aaaa_vooo_36"]("a,i,j,k") * V_blks_["aaaa_oovv"]("k,m,e,c");

            // rl1_aa_vo += 1.000000 <j,k||e,c>_abab t2_abab(a,c,i,k) u1_aa(b,j) m2_aaaa(i,m,b,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["aaaa_vooo_36"]("a,i,m,j") * sharedOps["aaaa_vvoo_6"]("a,e,i,j");

            // rl1_aa_vo += 1.000000 <k,j||e,c>_aaaa t2_aaaa(c,a,i,k) u1_aa(b,j) m2_aaaa(i,m,b,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["aaaa_vooo_36"]("a,i,m,j") * sharedOps["aaaa_vvoo_5"]("e,a,j,i");

            // tempOps["bbbb_vooo_37"] += 1.000000 m2_xbbbb_Loovv("I,i,j,b,a") u1_bb_vo("b,k")
            // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
            tempOps["bbbb_vooo_37"]("a,i,j,k") = m2_bbbb_oovv("i,j,b,a") * u1_bb_vo("b,k");

            // rl1_bb_vo += 1.000000 <k,j||c,e>_abab t2_abab(c,a,k,i) u1_bb(b,j) m2_bbbb(i,m,b,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["bbbb_vooo_37"]("a,i,m,j") * sharedOps["bbbb_vvoo_1"]("e,a,j,i");

            // rl1_bb_vo += 0.500000 <k,m||e,c>_bbbb t2_bbbb(c,a,i,j) u1_bb(b,k) m2_bbbb(i,j,b,a)
            // flops: o3v2L1: 1, o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 3,
            rl1_bb_vo("e,m") += 0.500000 * t2_bbbb_vvoo("c,a,i,j") * tempOps["bbbb_vooo_37"]("a,i,j,k") * V_blks_["bbbb_oovv"]("k,m,e,c");

            // rl1_bb_vo += 1.000000 <k,j||e,c>_bbbb t2_bbbb(c,a,i,k) u1_bb(b,j) m2_bbbb(i,m,b,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["bbbb_vooo_37"]("a,i,m,j") * sharedOps["bbbb_vvoo_2"]("a,e,i,j");

            // rl1_aa_vo += -0.500000 <m,k||e,c>_abab t2_bbbb(c,a,i,j) u1_bb(b,k) m2_bbbb(i,j,b,a)
            // flops: o3v2L1: 1, o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 3,
            rl1_aa_vo("e,m") -= 0.500000 * t2_bbbb_vvoo("c,a,i,j") * tempOps["bbbb_vooo_37"]("a,i,j,k") * V_blks_["abab_oovv"]("m,k,e,c");
        }
        tempOps["aaaa_vooo_35"] = TArrayD();
        tempOps["aaaa_vooo_36"] = TArrayD();
        tempOps["bbbb_vooo_37"] = TArrayD();

        {

            // tempOps["bbbb_vvoo_38"] += 1.000000 l2_xbbbb_Loovv("I,m,i,e,f") dp_bb_oo("n,i")
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempOps["bbbb_vvoo_38"]("e,f,m,n") = l2_bbbb_oovv("m,i,e,f") * dp_bb_oo("n,i");

            // rl2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) dp_bb(n,i) l2_bbbb(m,i,e,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempOps["bbbb_vvoo_38"]("e,f,m,n");
        }

        if (include_u2_) {

            // rm2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rm2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u0_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) dp_bb(n,i) l2_bbbb(m,i,e,f) u0
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempOps["bbbb_vvoo_38"]("e,f,m,n") * u0;
        }
        tempOps["bbbb_vvoo_38"] = TArrayD();

        if (include_u2_) {

            // tempOps["bbbb_vvoo_39"] += 1.000000 m2_xbbbb_Loovv("I,m,i,e,f") dp_bb_oo("n,i")
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempOps["bbbb_vvoo_39"]("e,f,m,n") = m2_bbbb_oovv("m,i,e,f") * dp_bb_oo("n,i");
        }

        {

            // rl2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u0_ && include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) dp_bb(n,i) u0 m2_bbbb(m,i,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempOps["bbbb_vvoo_39"]("e,f,m,n") * u0;
        }

        if (include_u2_) {

            // rm2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rm2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) dp_bb(n,i) m2_bbbb(m,i,e,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempOps["bbbb_vvoo_39"]("e,f,m,n");

            // tempOps["aaaa_vvoo_40"] += 1.000000 m2_xaaaa_Loovv("I,m,i,e,f") dp_aa_oo("n,i")
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempOps["aaaa_vvoo_40"]("e,f,m,n") = m2_aaaa_oovv("m,i,e,f") * dp_aa_oo("n,i");
        }
        tempOps["bbbb_vvoo_39"] = TArrayD();

        {

            // rl2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u0_ && include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) dp_aa(n,i) u0 m2_aaaa(m,i,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempOps["aaaa_vvoo_40"]("e,f,m,n") * u0;
        }

        if (include_u2_) {

            // rm2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rm2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) dp_aa(n,i) m2_aaaa(m,i,e,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempOps["aaaa_vvoo_40"]("e,f,m,n");
        }
        tempOps["aaaa_vvoo_40"] = TArrayD();

        {

            // tempOps["aaaa_vvoo_41"] += 1.000000 l2_xaaaa_Loovv("I,m,i,e,f") dp_aa_oo("n,i")
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempOps["aaaa_vvoo_41"]("e,f,m,n") = l2_aaaa_oovv("m,i,e,f") * dp_aa_oo("n,i");

            // rl2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) dp_aa(n,i) l2_aaaa(m,i,e,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempOps["aaaa_vvoo_41"]("e,f,m,n");
        }

        if (include_u2_) {

            // rm2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rm2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u0_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) dp_aa(n,i) l2_aaaa(m,i,e,f) u0
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempOps["aaaa_vvoo_41"]("e,f,m,n") * u0;
        }
        tempOps["aaaa_vvoo_41"] = TArrayD();

        {

            // tempOps["abab_vvoo_42"] += 1.000000 l2_xabab_Loovv("I,i,n,e,f") dp_aa_oo("m,i")
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempOps["abab_vvoo_42"]("e,f,m,n") = l2_abab_oovv("i,n,e,f") * dp_aa_oo("m,i");
        }

        if (include_u2_) {

            // rm2_abab_vvoo += 1.000000 dp_aa(m,i) l2_abab(i,n,e,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_abab_vvoo("e,f,m,n") += tempOps["abab_vvoo_42"]("e,f,m,n");
        }

        if (include_u0_) {

            // rl2_abab_vvoo += 1.000000 dp_aa(m,i) l2_abab(i,n,e,f) u0
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rl2_abab_vvoo("e,f,m,n") += tempOps["abab_vvoo_42"]("e,f,m,n") * u0;
        }
        tempOps["abab_vvoo_42"] = TArrayD();

        {

            // tempOps["abab_vvoo_43"] += 1.000000 dp_bb_oo("n,i") l2_xabab_Loovv("I,m,i,e,f")
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempOps["abab_vvoo_43"]("e,f,m,n") = dp_bb_oo("n,i") * l2_abab_oovv("m,i,e,f");
        }

        if (include_u2_) {

            // rm2_abab_vvoo += 1.000000 dp_bb(n,i) l2_abab(m,i,e,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_abab_vvoo("e,f,m,n") += tempOps["abab_vvoo_43"]("e,f,m,n");
        }

        if (include_u0_) {

            // rl2_abab_vvoo += 1.000000 dp_bb(n,i) l2_abab(m,i,e,f) u0
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rl2_abab_vvoo("e,f,m,n") += tempOps["abab_vvoo_43"]("e,f,m,n") * u0;
        }
        tempOps["abab_vvoo_43"] = TArrayD();

        if (include_u2_) {

            // tempOps["abab_vvoo_44"] += 1.000000 m2_xabab_Loovv("I,m,i,e,f") dp_bb_oo("n,i")
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempOps["abab_vvoo_44"]("e,f,m,n") = m2_abab_oovv("m,i,e,f") * dp_bb_oo("n,i");
        }

        if (include_u0_ && include_u2_) {

            // rm2_abab_vvoo += 1.000000 dp_bb(n,i) u0 m2_abab(m,i,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_abab_vvoo("e,f,m,n") += tempOps["abab_vvoo_44"]("e,f,m,n") * u0;
        }

        if (include_u2_) {

            // rl2_abab_vvoo += 1.000000 dp_bb(n,i) m2_abab(m,i,e,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_abab_vvoo("e,f,m,n") += tempOps["abab_vvoo_44"]("e,f,m,n");

            // tempOps["abab_vvoo_45"] += 1.000000 dp_aa_oo("m,i") m2_xabab_Loovv("I,i,n,e,f")
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempOps["abab_vvoo_45"]("e,f,m,n") = dp_aa_oo("m,i") * m2_abab_oovv("i,n,e,f");
        }
        tempOps["abab_vvoo_44"] = TArrayD();

        if (include_u0_ && include_u2_) {

            // rm2_abab_vvoo += 1.000000 dp_aa(m,i) u0 m2_abab(i,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_abab_vvoo("e,f,m,n") += tempOps["abab_vvoo_45"]("e,f,m,n") * u0;
        }

        if (include_u2_) {

            // rl2_abab_vvoo += 1.000000 dp_aa(m,i) m2_abab(i,n,e,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_abab_vvoo("e,f,m,n") += tempOps["abab_vvoo_45"]("e,f,m,n");

            // scalar_tmp_46 += 1.000000 sharedOps["abab_vvoo_397"]("b,a,i,j") m2_xabab_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_46 = 0.062500 * scalar_17;

            // energy += 0.062500 <k,l||c,d>_abab t2_abab(c,d,i,j) u2_abab(b,a,k,l) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_46;

            // energy += 0.062500 <l,k||c,d>_abab t2_abab(c,d,i,j) u2_abab(b,a,l,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_46;

            // energy += 0.062500 <k,l||d,c>_abab t2_abab(d,c,j,i) u2_abab(a,b,k,l) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_46;

            // energy += 0.062500 <l,k||d,c>_abab t2_abab(d,c,i,j) u2_abab(a,b,l,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_46;

            // energy += 0.062500 <l,k||d,c>_abab t2_abab(d,c,j,i) u2_abab(b,a,l,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_46;

            // energy += 0.062500 <l,k||c,d>_abab t2_abab(c,d,j,i) u2_abab(a,b,l,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_46;

            // energy += 0.062500 <k,l||c,d>_abab t2_abab(c,d,j,i) u2_abab(a,b,k,l) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_46;

            // energy += 0.062500 <k,l||c,d>_abab t2_abab(c,d,i,j) u2_abab(a,b,k,l) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_46;

            // energy += 0.062500 <l,k||c,d>_abab t2_abab(c,d,j,i) u2_abab(b,a,l,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_46;

            // energy += 0.062500 <k,l||c,d>_abab t2_abab(c,d,j,i) u2_abab(b,a,k,l) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_46;

            // energy += 0.062500 <k,l||d,c>_abab t2_abab(d,c,j,i) u2_abab(b,a,k,l) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_46;

            // energy += 0.062500 <k,l||d,c>_abab t2_abab(d,c,i,j) u2_abab(b,a,k,l) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_46;

            // energy += 0.062500 <k,l||d,c>_abab t2_abab(d,c,i,j) u2_abab(a,b,k,l) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_46;

            // energy += 0.062500 <l,k||d,c>_abab t2_abab(d,c,i,j) u2_abab(b,a,l,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_46;

            // energy += 0.062500 <l,k||d,c>_abab t2_abab(d,c,j,i) u2_abab(a,b,l,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_46;

            // energy += 0.062500 <l,k||c,d>_abab t2_abab(c,d,i,j) u2_abab(a,b,l,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_46;
        }
        tempOps["abab_vvoo_45"] = TArrayD();

        {

            // scalar_tmp_47 += 1.000000 sharedOps["abab_vvoo_398"]("b,a,i,j") l2_xabab_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_47 = 0.062500 * scalar_18;

            // energy += 0.062500 <l,k||c,d>_abab l2_abab(i,j,b,a) t2_abab(c,d,i,j) t2_abab(b,a,l,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_47;

            // energy += 0.062500 <l,k||d,c>_abab l2_abab(i,j,a,b) t2_abab(d,c,i,j) t2_abab(a,b,l,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_47;

            // energy += 0.062500 <l,k||c,d>_abab l2_abab(i,j,a,b) t2_abab(c,d,i,j) t2_abab(a,b,l,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_47;

            // energy += 0.062500 <l,k||d,c>_abab l2_abab(j,i,b,a) t2_abab(d,c,j,i) t2_abab(b,a,l,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_47;

            // energy += 0.062500 <l,k||d,c>_abab l2_abab(j,i,a,b) t2_abab(d,c,j,i) t2_abab(a,b,l,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_47;

            // energy += 0.062500 <k,l||d,c>_abab l2_abab(i,j,b,a) t2_abab(d,c,i,j) t2_abab(b,a,k,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_47;

            // energy += 0.062500 <l,k||d,c>_abab l2_abab(i,j,b,a) t2_abab(d,c,i,j) t2_abab(b,a,l,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_47;

            // energy += 0.062500 <k,l||c,d>_abab l2_abab(i,j,b,a) t2_abab(c,d,i,j) t2_abab(b,a,k,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_47;

            // energy += 0.062500 <k,l||d,c>_abab l2_abab(j,i,a,b) t2_abab(d,c,j,i) t2_abab(a,b,k,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_47;

            // energy += 0.062500 <k,l||c,d>_abab l2_abab(j,i,a,b) t2_abab(c,d,j,i) t2_abab(a,b,k,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_47;

            // energy += 0.062500 <k,l||d,c>_abab l2_abab(j,i,b,a) t2_abab(d,c,j,i) t2_abab(b,a,k,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_47;

            // energy += 0.062500 <l,k||c,d>_abab l2_abab(j,i,a,b) t2_abab(c,d,j,i) t2_abab(a,b,l,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_47;

            // energy += 0.062500 <k,l||c,d>_abab l2_abab(i,j,a,b) t2_abab(c,d,i,j) t2_abab(a,b,k,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_47;

            // energy += 0.062500 <k,l||d,c>_abab l2_abab(i,j,a,b) t2_abab(d,c,i,j) t2_abab(a,b,k,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_47;

            // energy += 0.062500 <l,k||c,d>_abab l2_abab(j,i,b,a) t2_abab(c,d,j,i) t2_abab(b,a,l,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_47;

            // energy += 0.062500 <k,l||c,d>_abab l2_abab(j,i,b,a) t2_abab(c,d,j,i) t2_abab(b,a,k,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_47;
        }

        if (include_u2_) {

            // scalar_tmp_48 += 1.000000 sharedOps["abab_vvoo_104"]("a,b,j,i") m2_xabab_Loovv("I,j,i,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_48 = 0.062500 * scalar_19;

            // energy += 0.062500 <k,l||d,c>_abab t2_abab(a,b,k,l) u2_abab(d,c,j,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_48;

            // energy += 0.062500 <l,k||d,c>_abab t2_abab(a,b,l,k) u2_abab(d,c,j,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_48;

            // energy += 0.062500 <l,k||d,c>_abab t2_abab(a,b,l,k) u2_abab(d,c,i,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_48;

            // energy += 0.062500 <k,l||c,d>_abab t2_abab(a,b,k,l) u2_abab(c,d,i,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_48;

            // energy += 0.062500 <l,k||d,c>_abab t2_abab(b,a,l,k) u2_abab(d,c,i,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_48;

            // energy += 0.062500 <k,l||d,c>_abab t2_abab(b,a,k,l) u2_abab(d,c,j,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_48;

            // energy += 0.062500 <l,k||c,d>_abab t2_abab(b,a,l,k) u2_abab(c,d,j,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_48;

            // energy += 0.062500 <k,l||d,c>_abab t2_abab(b,a,k,l) u2_abab(d,c,i,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_48;

            // energy += 0.062500 <k,l||c,d>_abab t2_abab(b,a,k,l) u2_abab(c,d,j,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_48;

            // energy += 0.062500 <k,l||c,d>_abab t2_abab(a,b,k,l) u2_abab(c,d,j,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_48;

            // energy += 0.062500 <k,l||c,d>_abab t2_abab(b,a,k,l) u2_abab(c,d,i,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_48;

            // energy += 0.062500 <l,k||c,d>_abab t2_abab(a,b,l,k) u2_abab(c,d,i,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_48;

            // energy += 0.062500 <l,k||d,c>_abab t2_abab(b,a,l,k) u2_abab(d,c,j,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_48;

            // energy += 0.062500 <l,k||c,d>_abab t2_abab(b,a,l,k) u2_abab(c,d,i,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_48;

            // energy += 0.062500 <l,k||c,d>_abab t2_abab(a,b,l,k) u2_abab(c,d,j,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_48;

            // energy += 0.062500 <k,l||d,c>_abab t2_abab(a,b,k,l) u2_abab(d,c,i,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_48;

            // scalar_tmp_49 += 1.000000 u2_abab_vvoo("a,b,i,j") m2_xabab_Loovv("I,i,j,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_49 = 0.250000 * scalar_20;
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += -0.250000 dp_bb(m,e) u2_abab(b,a,i,j) m2_abab(i,j,b,a)
            // +            -0.250000 dp_bb(m,e) u2_abab(b,a,j,i) m2_abab(j,i,b,a)
            // +            -0.250000 dp_bb(m,e) u2_abab(a,b,i,j) m2_abab(i,j,a,b)
            // +            -0.250000 dp_bb(m,e) u2_abab(a,b,j,i) m2_abab(j,i,a,b)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= 4.000000 * scalar_tmp_49 * dp_bb_ov("m,e");

            // rm1_aa_vo += -0.250000 dp_aa(m,e) u2_abab(b,a,i,j) m2_abab(i,j,b,a)
            // +            -0.250000 dp_aa(m,e) u2_abab(b,a,j,i) m2_abab(j,i,b,a)
            // +            -0.250000 dp_aa(m,e) u2_abab(a,b,i,j) m2_abab(i,j,a,b)
            // +            -0.250000 dp_aa(m,e) u2_abab(a,b,j,i) m2_abab(j,i,a,b)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= 4.000000 * scalar_tmp_49 * dp_aa_ov("m,e");
        }

        if (include_u2_) {

            // energy += 0.250000 u2_abab(a,b,i,j) m2_abab(i,j,a,b) w0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_49 * w0;
        }

        if (include_u1_ && include_u2_) {

            // energy += -0.250000 dp_bb(k,c) u1_bb(c,k) u2_abab(a,b,i,j) m2_abab(i,j,a,b)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_7 * scalar_tmp_49;
        }

        if (include_u2_) {

            // energy += 0.250000 u2_abab(b,a,j,i) m2_abab(j,i,b,a) w0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_49 * w0;
        }

        if (include_u1_ && include_u2_) {

            // energy += -0.250000 dp_bb(k,c) u1_bb(c,k) u2_abab(a,b,j,i) m2_abab(j,i,a,b)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_7 * scalar_tmp_49;
        }

        if (include_u2_) {

            // energy += 0.250000 u2_abab(a,b,j,i) m2_abab(j,i,a,b) w0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_49 * w0;

            // energy += 0.250000 u2_abab(b,a,i,j) m2_abab(i,j,b,a) w0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_49 * w0;
        }

        if (include_u1_ && include_u2_) {

            // energy += -0.250000 dp_aa(k,c) u1_aa(c,k) u2_abab(b,a,i,j) m2_abab(i,j,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_6 * scalar_tmp_49;

            // energy += -0.250000 dp_bb(k,c) u1_bb(c,k) u2_abab(b,a,j,i) m2_abab(j,i,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_7 * scalar_tmp_49;

            // energy += -0.250000 dp_aa(k,c) u1_aa(c,k) u2_abab(a,b,j,i) m2_abab(j,i,a,b)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_6 * scalar_tmp_49;

            // energy += -0.250000 dp_aa(k,c) u1_aa(c,k) u2_abab(b,a,j,i) m2_abab(j,i,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_6 * scalar_tmp_49;

            // energy += -0.250000 dp_aa(k,c) u1_aa(c,k) u2_abab(a,b,i,j) m2_abab(i,j,a,b)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_6 * scalar_tmp_49;

            // energy += -0.250000 dp_bb(k,c) u1_bb(c,k) u2_abab(b,a,i,j) m2_abab(i,j,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_7 * scalar_tmp_49;
        }

        if (include_u2_) {

            // scalar_tmp_50 += 1.000000 l2_xabab_Loovv("I,k,j,a,b") u2_abab_vvoo("a,b,k,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_50 = 0.250000 * scalar_21;

            // rl1_bb_vo += -0.250000 dp_bb(m,e) l2_abab(j,i,a,b) u2_abab(a,b,j,i)
            // +            -0.250000 dp_bb(m,e) l2_abab(j,i,b,a) u2_abab(b,a,j,i)
            // +            -0.250000 dp_bb(m,e) l2_abab(i,j,a,b) u2_abab(a,b,i,j)
            // +            -0.250000 dp_bb(m,e) l2_abab(i,j,b,a) u2_abab(b,a,i,j)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= 4.000000 * scalar_tmp_50 * dp_bb_ov("m,e");

            // rl1_aa_vo += -0.250000 dp_aa(m,e) l2_abab(j,i,a,b) u2_abab(a,b,j,i)
            // +            -0.250000 dp_aa(m,e) l2_abab(j,i,b,a) u2_abab(b,a,j,i)
            // +            -0.250000 dp_aa(m,e) l2_abab(i,j,a,b) u2_abab(a,b,i,j)
            // +            -0.250000 dp_aa(m,e) l2_abab(i,j,b,a) u2_abab(b,a,i,j)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= 4.000000 * scalar_tmp_50 * dp_aa_ov("m,e");

            // cenergy += 0.250000 l2_abab(i,j,b,a) u2_abab(b,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            cenergy += scalar_tmp_50;

            // cenergy += 0.250000 l2_abab(i,j,a,b) u2_abab(a,b,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            cenergy += scalar_tmp_50;

            // cenergy += 0.250000 l2_abab(j,i,b,a) u2_abab(b,a,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            cenergy += scalar_tmp_50;

            // cenergy += 0.250000 l2_abab(j,i,a,b) u2_abab(a,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            cenergy += scalar_tmp_50;

            // energy += -0.250000 dp_bb(i,i) l2_abab(j,k,a,b) u2_abab(a,b,j,k)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_5 * scalar_tmp_50;

            // energy += -0.250000 dp_aa(i,i) l2_abab(j,k,a,b) u2_abab(a,b,j,k)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_4 * scalar_tmp_50;

            // energy += -0.250000 dp_aa(i,i) l2_abab(k,j,a,b) u2_abab(a,b,k,j)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_4 * scalar_tmp_50;

            // energy += -0.250000 dp_bb(i,i) l2_abab(k,j,a,b) u2_abab(a,b,k,j)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_5 * scalar_tmp_50;

            // energy += -0.250000 dp_aa(i,i) l2_abab(k,j,b,a) u2_abab(b,a,k,j)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_4 * scalar_tmp_50;

            // energy += -0.250000 dp_bb(i,i) l2_abab(j,k,b,a) u2_abab(b,a,j,k)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_5 * scalar_tmp_50;

            // energy += -0.250000 dp_aa(i,i) l2_abab(j,k,b,a) u2_abab(b,a,j,k)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_4 * scalar_tmp_50;

            // energy += -0.250000 dp_bb(i,i) l2_abab(k,j,b,a) u2_abab(b,a,k,j)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_5 * scalar_tmp_50;

            // scalar_tmp_51 += 1.000000 sharedOps["abab_vvoo_10"]("b,a,j,i") m2_xabab_Loovv("I,j,i,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_51 = 0.125000 * scalar_22;

            // energy += 0.125000 <a,b||c,d>_abab u2_abab(c,d,j,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_51;

            // energy += 0.125000 <b,a||d,c>_abab u2_abab(d,c,i,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_51;

            // energy += 0.125000 <a,b||d,c>_abab u2_abab(d,c,i,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_51;

            // energy += 0.125000 <b,a||d,c>_abab u2_abab(d,c,j,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_51;

            // energy += 0.125000 <a,b||d,c>_abab u2_abab(d,c,j,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_51;

            // energy += 0.125000 <a,b||c,d>_abab u2_abab(c,d,i,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_51;

            // energy += 0.125000 <b,a||c,d>_abab u2_abab(c,d,j,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_51;

            // energy += 0.125000 <b,a||c,d>_abab u2_abab(c,d,i,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_51;
        }

        {

            // scalar_tmp_52 += 1.000000 l2_xabab_Loovv("I,j,i,b,a") sharedOps["abab_vvoo_105"]("b,a,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_52 = 0.125000 * scalar_23;

            // energy += 0.125000 <k,l||j,i>_abab l2_abab(j,i,a,b) t2_abab(a,b,k,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_52;

            // energy += 0.125000 <k,l||i,j>_abab l2_abab(i,j,a,b) t2_abab(a,b,k,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_52;

            // energy += 0.125000 <k,l||i,j>_abab l2_abab(i,j,b,a) t2_abab(b,a,k,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_52;

            // energy += 0.125000 <k,l||j,i>_abab l2_abab(j,i,b,a) t2_abab(b,a,k,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_52;

            // energy += 0.125000 <l,k||j,i>_abab l2_abab(j,i,b,a) t2_abab(b,a,l,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_52;

            // energy += 0.125000 <l,k||j,i>_abab l2_abab(j,i,a,b) t2_abab(a,b,l,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_52;

            // energy += 0.125000 <l,k||i,j>_abab l2_abab(i,j,b,a) t2_abab(b,a,l,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_52;

            // energy += 0.125000 <l,k||i,j>_abab l2_abab(i,j,a,b) t2_abab(a,b,l,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_52;
        }

        if (include_u2_) {

            // scalar_tmp_53 += 1.000000 m2_xabab_Loovv("I,i,j,b,a") sharedOps["abab_vvoo_106"]("b,a,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_53 = 0.125000 * scalar_24;

            // energy += 0.125000 <k,l||i,j>_abab u2_abab(b,a,k,l) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_53;

            // energy += 0.125000 <l,k||i,j>_abab u2_abab(a,b,l,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_53;

            // energy += 0.125000 <k,l||j,i>_abab u2_abab(a,b,k,l) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_53;

            // energy += 0.125000 <k,l||i,j>_abab u2_abab(a,b,k,l) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_53;

            // energy += 0.125000 <l,k||j,i>_abab u2_abab(a,b,l,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_53;

            // energy += 0.125000 <k,l||j,i>_abab u2_abab(b,a,k,l) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_53;

            // energy += 0.125000 <l,k||j,i>_abab u2_abab(b,a,l,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_53;

            // energy += 0.125000 <l,k||i,j>_abab u2_abab(b,a,l,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_53;
        }

        {

            // scalar_tmp_54 += 1.000000 l2_xabab_Loovv("I,i,j,a,b") sharedOps["abab_vvoo_11"]("a,b,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_54 = 0.125000 * scalar_25;

            // energy += 0.125000 <a,b||c,d>_abab l2_abab(i,j,a,b) t2_abab(c,d,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_54;

            // energy += 0.125000 <b,a||d,c>_abab l2_abab(j,i,b,a) t2_abab(d,c,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_54;

            // energy += 0.125000 <b,a||c,d>_abab l2_abab(j,i,b,a) t2_abab(c,d,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_54;

            // energy += 0.125000 <a,b||d,c>_abab l2_abab(j,i,a,b) t2_abab(d,c,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_54;

            // energy += 0.125000 <a,b||d,c>_abab l2_abab(i,j,a,b) t2_abab(d,c,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_54;

            // energy += 0.125000 <b,a||d,c>_abab l2_abab(i,j,b,a) t2_abab(d,c,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_54;

            // energy += 0.125000 <b,a||c,d>_abab l2_abab(i,j,b,a) t2_abab(c,d,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_54;

            // energy += 0.125000 <a,b||c,d>_abab l2_abab(j,i,a,b) t2_abab(c,d,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_54;
        }

        if (include_u1_ && include_u2_) {

            // tempOps["bb_vo_55"] += 1.000000 u1_bb_vo("a,i") m2_xbbbb_Loovv("I,i,j,e,a")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["bb_vo_55"]("e,j") = u1_bb_vo("a,i") * m2_bbbb_oovv("i,j,e,a");
        }

        {

            // rl2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) dp_bb(n,e) u1_bb(a,i) m2_bbbb(i,m,f,a)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempOps["bb_vo_55"]("f,m") * dp_bb_ov("n,e");

            // rm2_abab_vvoo += 1.000000 dp_aa(m,e) u1_bb(a,i) m2_bbbb(i,n,f,a)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += tempOps["bb_vo_55"]("f,n") * dp_aa_ov("m,e");

            // rm1_bb_vo += 1.000000 dp_bb(b,e) u1_bb(a,i) m2_bbbb(i,m,b,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += tempOps["bb_vo_55"]("b,m") * dp_bb_vv("b,e");

            // rm1_bb_vo += -1.000000 dp_bb(m,j) u1_bb(a,i) m2_bbbb(i,j,e,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= tempOps["bb_vo_55"]("e,j") * dp_bb_oo("m,j");

            // rl1_bb_vo += -1.000000 dp_bb(j,e) u1_bb(b,j) u1_bb(a,i) m2_bbbb(i,m,b,a)
            // flops: o2v1L1: 2, o1v1L1: 1 | mem: o1v1L1: 2, o2v0L1: 1,
            rl1_bb_vo("e,m") -= u1_bb_vo("b,j") * tempOps["bb_vo_55"]("b,m") * dp_bb_ov("j,e");

            // rl1_bb_vo += -1.000000 dp_bb(m,b) u1_bb(b,j) u1_bb(a,i) m2_bbbb(i,j,e,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bb_vo_55"]("e,j") * sharedOps["bb_oo_313"]("m,j");

            // tempOps["bb_vo_56"] += 1.000000 u1_aa_vo("a,i") m2_xabab_Loovv("I,i,j,a,e")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["bb_vo_56"]("e,j") = u1_aa_vo("a,i") * m2_abab_oovv("i,j,a,e");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) dp_bb(n,e) u1_aa(a,i) m2_abab(i,m,a,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += tempOps["bb_vo_56"]("f,m") * dp_bb_ov("n,e");

            // rm2_abab_vvoo += -1.000000 dp_aa(m,e) u1_aa(a,i) m2_abab(i,n,a,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= tempOps["bb_vo_56"]("f,n") * dp_aa_ov("m,e");

            // rm1_bb_vo += -1.000000 dp_bb(b,e) u1_aa(a,i) m2_abab(i,m,a,b)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= tempOps["bb_vo_56"]("b,m") * dp_bb_vv("b,e");

            // rm1_bb_vo += 1.000000 dp_bb(m,j) u1_aa(a,i) m2_abab(i,j,a,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += tempOps["bb_vo_56"]("e,j") * dp_bb_oo("m,j");

            // rl1_bb_vo += 1.000000 dp_bb(j,e) u1_bb(b,j) u1_aa(a,i) m2_abab(i,m,a,b)
            // flops: o2v1L1: 2, o1v1L1: 1 | mem: o1v1L1: 2, o2v0L1: 1,
            rl1_bb_vo("e,m") += u1_bb_vo("b,j") * tempOps["bb_vo_56"]("b,m") * dp_bb_ov("j,e");

            // rl1_bb_vo += 1.000000 dp_bb(m,b) u1_bb(b,j) u1_aa(a,i) m2_abab(i,j,a,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["bb_vo_56"]("e,j") * sharedOps["bb_oo_313"]("m,j");

            // tempOps["aa_vo_57"] += 1.000000 m2_xabab_Loovv("I,j,i,e,a") u1_bb_vo("a,i")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["aa_vo_57"]("e,j") = m2_abab_oovv("j,i,e,a") * u1_bb_vo("a,i");

            // rm2_abab_vvoo += -1.000000 dp_bb(n,f) u1_bb(a,i) m2_abab(m,i,e,a)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= tempOps["aa_vo_57"]("e,m") * dp_bb_ov("n,f");
        }
        tempOps["bb_vo_55"] = TArrayD();
        tempOps["bb_vo_56"] = TArrayD();

        {

            // rl2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) dp_aa(n,e) u1_bb(a,i) m2_abab(m,i,f,a)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempOps["aa_vo_57"]("f,m") * dp_aa_ov("n,e");

            // rm1_aa_vo += -1.000000 dp_aa(b,e) u1_bb(a,i) m2_abab(m,i,b,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= tempOps["aa_vo_57"]("b,m") * dp_aa_vv("b,e");

            // rm1_aa_vo += 1.000000 dp_aa(m,j) u1_bb(a,i) m2_abab(j,i,e,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += tempOps["aa_vo_57"]("e,j") * dp_aa_oo("m,j");

            // rl1_aa_vo += 1.000000 dp_aa(m,b) u1_aa(b,j) u1_bb(a,i) m2_abab(j,i,e,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["aa_vo_57"]("e,j") * sharedOps["aa_oo_312"]("j,m");

            // rl1_aa_vo += 1.000000 dp_aa(j,e) u1_aa(b,j) u1_bb(a,i) m2_abab(m,i,b,a)
            // flops: o2v1L1: 2, o1v1L1: 1 | mem: o1v1L1: 2, o2v0L1: 1,
            rl1_aa_vo("e,m") += u1_aa_vo("b,j") * tempOps["aa_vo_57"]("b,m") * dp_aa_ov("j,e");

            // tempOps["aa_vo_58"] += 1.000000 u1_aa_vo("a,i") m2_xaaaa_Loovv("I,i,j,e,a")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["aa_vo_58"]("e,j") = u1_aa_vo("a,i") * m2_aaaa_oovv("i,j,e,a");

            // rm2_abab_vvoo += 1.000000 dp_bb(n,f) u1_aa(a,i) m2_aaaa(i,m,e,a)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += tempOps["aa_vo_58"]("e,m") * dp_bb_ov("n,f");

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) dp_aa(n,e) u1_aa(a,i) m2_aaaa(i,m,f,a)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= tempOps["aa_vo_58"]("f,m") * dp_aa_ov("n,e");

            // rm1_aa_vo += 1.000000 dp_aa(b,e) u1_aa(a,i) m2_aaaa(i,m,b,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += tempOps["aa_vo_58"]("b,m") * dp_aa_vv("b,e");

            // rm1_aa_vo += -1.000000 dp_aa(m,j) u1_aa(a,i) m2_aaaa(i,j,e,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= tempOps["aa_vo_58"]("e,j") * dp_aa_oo("m,j");

            // rl1_aa_vo += -1.000000 dp_aa(j,e) u1_aa(b,j) u1_aa(a,i) m2_aaaa(i,m,b,a)
            // flops: o2v1L1: 2, o1v1L1: 1 | mem: o1v1L1: 2, o2v0L1: 1,
            rl1_aa_vo("e,m") -= u1_aa_vo("b,j") * tempOps["aa_vo_58"]("b,m") * dp_aa_ov("j,e");

            // rl1_aa_vo += -1.000000 dp_aa(m,b) u1_aa(b,j) u1_aa(a,i) m2_aaaa(i,j,e,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aa_vo_58"]("e,j") * sharedOps["aa_oo_312"]("j,m");
        }
        tempOps["aa_vo_57"] = TArrayD();
        tempOps["aa_vo_58"] = TArrayD();

        if (include_u2_) {

            // scalar_tmp_59 += 1.000000 u2_aaaa_vvoo("b,a,i,j") m2_xaaaa_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_59 = 0.250000 * scalar_26;
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += -0.250000 dp_bb(m,e) u2_aaaa(b,a,i,j) m2_aaaa(i,j,b,a)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= scalar_tmp_59 * dp_bb_ov("m,e");

            // rm1_aa_vo += -0.250000 dp_aa(m,e) u2_aaaa(b,a,i,j) m2_aaaa(i,j,b,a)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= scalar_tmp_59 * dp_aa_ov("m,e");

            // energy += -0.250000 dp_aa(k,c) u1_aa(c,k) u2_aaaa(b,a,i,j) m2_aaaa(i,j,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_6 * scalar_tmp_59;
        }

        if (include_u2_) {

            // energy += 0.250000 u2_aaaa(b,a,i,j) m2_aaaa(i,j,b,a) w0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_59 * w0;
        }

        if (include_u1_ && include_u2_) {

            // energy += -0.250000 dp_bb(k,c) u1_bb(c,k) u2_aaaa(b,a,i,j) m2_aaaa(i,j,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_7 * scalar_tmp_59;
        }

        if (include_u2_) {

            // scalar_tmp_60 += 1.000000 u2_bbbb_vvoo("a,b,k,j") l2_xbbbb_Loovv("I,k,j,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_60 = 0.250000 * scalar_27;

            // rl1_bb_vo += -0.250000 dp_bb(m,e) l2_bbbb(j,i,a,b) u2_bbbb(a,b,j,i)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= scalar_tmp_60 * dp_bb_ov("m,e");

            // rl1_aa_vo += -0.250000 dp_aa(m,e) l2_bbbb(j,i,a,b) u2_bbbb(a,b,j,i)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= scalar_tmp_60 * dp_aa_ov("m,e");

            // cenergy += 0.250000 l2_bbbb(j,i,a,b) u2_bbbb(a,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            cenergy += scalar_tmp_60;

            // energy += -0.250000 dp_aa(i,i) l2_bbbb(k,j,a,b) u2_bbbb(a,b,k,j)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_4 * scalar_tmp_60;

            // energy += -0.250000 dp_bb(i,i) l2_bbbb(k,j,a,b) u2_bbbb(a,b,k,j)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_5 * scalar_tmp_60;

            // scalar_tmp_61 += 1.000000 l2_xaaaa_Loovv("I,k,j,a,b") u2_aaaa_vvoo("a,b,k,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_61 = 0.250000 * scalar_28;

            // rl1_bb_vo += -0.250000 dp_bb(m,e) l2_aaaa(j,i,a,b) u2_aaaa(a,b,j,i)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= scalar_tmp_61 * dp_bb_ov("m,e");

            // rl1_aa_vo += -0.250000 dp_aa(m,e) l2_aaaa(j,i,a,b) u2_aaaa(a,b,j,i)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= scalar_tmp_61 * dp_aa_ov("m,e");

            // cenergy += 0.250000 l2_aaaa(j,i,a,b) u2_aaaa(a,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            cenergy += scalar_tmp_61;

            // energy += -0.250000 dp_aa(i,i) l2_aaaa(k,j,a,b) u2_aaaa(a,b,k,j)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_4 * scalar_tmp_61;

            // energy += -0.250000 dp_bb(i,i) l2_aaaa(k,j,a,b) u2_aaaa(a,b,k,j)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_5 * scalar_tmp_61;

            // scalar_tmp_62 += 1.000000 u2_bbbb_vvoo("b,a,i,j") m2_xbbbb_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_62 = 0.250000 * scalar_29;
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += -0.250000 dp_bb(m,e) u2_bbbb(b,a,i,j) m2_bbbb(i,j,b,a)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= scalar_tmp_62 * dp_bb_ov("m,e");

            // rm1_aa_vo += -0.250000 dp_aa(m,e) u2_bbbb(b,a,i,j) m2_bbbb(i,j,b,a)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= scalar_tmp_62 * dp_aa_ov("m,e");
        }

        if (include_u2_) {

            // energy += 0.250000 u2_bbbb(b,a,i,j) m2_bbbb(i,j,b,a) w0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_62 * w0;
        }

        if (include_u1_ && include_u2_) {

            // energy += -0.250000 dp_aa(k,c) u1_aa(c,k) u2_bbbb(b,a,i,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_6 * scalar_tmp_62;

            // energy += -0.250000 dp_bb(k,c) u1_bb(c,k) u2_bbbb(b,a,i,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_7 * scalar_tmp_62;
        }

        if (include_u1_) {

            // tempOps["bbbb_vvoo_63"] += 1.000000 dp_bb_ov("n,e") m1_xbb_Lov("I,m,f")
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempOps["bbbb_vvoo_63"]("e,f,n,m") = dp_bb_ov("n,e") * m1_bb_ov("m,f");
        }

        if (include_u0_ && include_u1_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) dp_bb(n,e) u0 m1_bb(m,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") += tempOps["bbbb_vvoo_63"]("e,f,n,m") * u0;
        }

        if (include_u2_) {

            // rm2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rm2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
            rm2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
            rm2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("f,e,n,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) dp_bb(n,e) m1_bb(m,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempOps["bbbb_vvoo_63"]("e,f,n,m");
        }
        tempOps["bbbb_vvoo_63"] = TArrayD();

        {

            // tempOps["bbbb_vvoo_64"] += 1.000000 dp_bb_ov("n,e") l1_xbb_Lov("I,m,f")
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempOps["bbbb_vvoo_64"]("e,f,n,m") = dp_bb_ov("n,e") * l1_bb_ov("m,f");

            // rl2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("f,e,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) dp_bb(n,e) l1_bb(m,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempOps["bbbb_vvoo_64"]("e,f,n,m");
        }

        if (include_u2_) {

            // rm2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rm2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
            rm2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
            rm2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("f,e,n,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u0_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) dp_bb(n,e) l1_bb(m,f) u0
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempOps["bbbb_vvoo_64"]("e,f,n,m") * u0;
        }
        tempOps["bbbb_vvoo_64"] = TArrayD();

        if (include_u1_) {

            // tempOps["aaaa_vvoo_65"] += 1.000000 dp_aa_ov("n,e") m1_xaa_Lov("I,m,f")
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempOps["aaaa_vvoo_65"]("e,f,n,m") = dp_aa_ov("n,e") * m1_aa_ov("m,f");
        }

        if (include_u0_ && include_u1_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) dp_aa(n,e) u0 m1_aa(m,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") += tempOps["aaaa_vvoo_65"]("e,f,n,m") * u0;
        }

        if (include_u2_) {

            // rm2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rm2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
            rm2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
            rm2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("f,e,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) dp_aa(n,e) m1_aa(m,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempOps["aaaa_vvoo_65"]("e,f,n,m");
        }
        tempOps["aaaa_vvoo_65"] = TArrayD();

        {

            // tempOps["aaaa_vvoo_66"] += 1.000000 l1_xaa_Lov("I,m,f") dp_aa_ov("n,e")
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempOps["aaaa_vvoo_66"]("f,e,m,n") = l1_aa_ov("m,f") * dp_aa_ov("n,e");

            // rl2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("f,e,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) dp_aa(n,e) l1_aa(m,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempOps["aaaa_vvoo_66"]("f,e,m,n");
        }

        if (include_u2_) {

            // rm2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rm2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
            rm2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
            rm2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("f,e,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u0_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) dp_aa(n,e) l1_aa(m,f) u0
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempOps["aaaa_vvoo_66"]("f,e,m,n") * u0;
        }
        tempOps["aaaa_vvoo_66"] = TArrayD();

        if (include_u2_) {

            // scalar_tmp_67 += 1.000000 sharedOps["abab_vvoo_411"]("a,b,j,i") m2_xabab_Loovv("I,j,i,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_67 = 0.250000 * scalar_30;

            // energy += -0.250000 <k,l||c,d>_abab t2_abab(a,d,j,i) u2_abab(c,b,k,l) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_67;

            // energy += -0.250000 <l,k||c,d>_abab t2_abab(a,d,i,j) u2_abab(c,b,l,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_67;

            // energy += -0.250000 <l,k||c,d>_abab t2_abab(a,d,j,i) u2_abab(c,b,l,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_67;

            // energy += -0.250000 <k,l||c,d>_abab t2_abab(a,d,i,j) u2_abab(c,b,k,l) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_67;

            // scalar_tmp_68 += 1.000000 m2_xabab_Loovv("I,i,j,a,b") sharedOps["abab_vvoo_407"]("a,b,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_68 = 0.250000 * scalar_31;

            // energy += -0.250000 <l,k||c,d>_abab t2_abab(a,d,l,k) u2_abab(c,b,i,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_68;

            // energy += -0.250000 <k,l||c,d>_abab t2_abab(a,d,k,l) u2_abab(c,b,i,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_68;

            // energy += -0.250000 <l,k||c,d>_abab t2_abab(a,d,l,k) u2_abab(c,b,j,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_68;

            // energy += -0.250000 <k,l||c,d>_abab t2_abab(a,d,k,l) u2_abab(c,b,j,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_68;
        }

        if (include_u0_ && include_u2_) {

            // scalar_tmp_69 += 1.000000 sharedOps["abab_vvoo_156"]("a,b,i,j") m2_xabab_Loovv("I,i,j,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_69 = 0.500000 * scalar_32;

            // rm0 += -0.500000 dp_bb(b,c) u2_abab(a,c,j,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_69;

            // rm0 += -0.500000 dp_bb(b,c) u2_abab(a,c,i,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_69;

            // energy += -0.500000 dp_bb(b,c) u0 u2_abab(a,c,j,i) m2_abab(j,i,a,b)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_69 * u0;

            // energy += -0.500000 dp_bb(b,c) u0 u2_abab(a,c,i,j) m2_abab(i,j,a,b)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_69 * u0;
        }

        if (include_u0_) {

            // scalar_tmp_70 += 1.000000 l2_xabab_Loovv("I,j,i,b,a") sharedOps["abab_vvoo_155"]("b,a,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_70 = 0.500000 * scalar_33;

            // rm0 += -0.500000 dp_aa(b,c) l2_abab(j,i,b,a) t2_abab(c,a,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_70;

            // rm0 += -0.500000 dp_aa(b,c) l2_abab(i,j,b,a) t2_abab(c,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_70;

            // energy += -0.500000 dp_aa(b,c) l2_abab(j,i,b,a) t2_abab(c,a,j,i) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_70 * u0;

            // energy += -0.500000 dp_aa(b,c) l2_abab(i,j,b,a) t2_abab(c,a,i,j) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_70 * u0;
        }

        {

            // scalar_tmp_71 += 1.000000 sharedOps["abab_vvoo_406"]("a,b,j,i") l2_xabab_Loovv("I,j,i,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_71 = 0.250000 * scalar_34;

            // energy += -0.250000 <l,k||c,d>_abab l2_abab(i,j,a,b) t2_abab(c,b,l,k) t2_abab(a,d,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_71;

            // energy += -0.250000 <k,l||c,d>_abab l2_abab(i,j,a,b) t2_abab(c,b,k,l) t2_abab(a,d,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_71;

            // energy += -0.250000 <k,l||c,d>_abab l2_abab(j,i,a,b) t2_abab(c,b,k,l) t2_abab(a,d,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_71;

            // energy += -0.250000 <l,k||c,d>_abab l2_abab(j,i,a,b) t2_abab(c,b,l,k) t2_abab(a,d,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_71;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_72 += 1.000000 m2_xabab_Loovv("I,j,i,b,a") sharedOps["abab_vvoo_423"]("b,a,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_72 = 0.250000 * scalar_35;

            // energy += -0.250000 <b,k||d,c>_abab t2_abab(d,c,i,j) u1_bb(a,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_72;

            // energy += -0.250000 <b,k||c,d>_abab t2_abab(c,d,i,j) u1_bb(a,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_72;

            // energy += -0.250000 <b,k||c,d>_abab t2_abab(c,d,j,i) u1_bb(a,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_72;

            // energy += -0.250000 <b,k||d,c>_abab t2_abab(d,c,j,i) u1_bb(a,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_72;
        }

        if (include_u0_) {

            // scalar_tmp_73 += 1.000000 l2_xabab_Loovv("I,j,i,a,b") sharedOps["abab_vvoo_220"]("a,b,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_73 = 0.500000 * scalar_36;

            // rm0 += 0.500000 dp_aa(k,j) l2_abab(j,i,b,a) t2_abab(b,a,k,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_73;

            // rm0 += 0.500000 dp_aa(k,j) l2_abab(j,i,a,b) t2_abab(a,b,k,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_73;

            // energy += 0.500000 dp_aa(k,j) l2_abab(j,i,a,b) t2_abab(a,b,k,i) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_73 * u0;

            // energy += 0.500000 dp_aa(k,j) l2_abab(j,i,b,a) t2_abab(b,a,k,i) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_73 * u0;
        }

        if (include_u2_) {

            // scalar_tmp_74 += 1.000000 sharedOps["abab_vvoo_429"]("a,b,i,j") m2_xabab_Loovv("I,i,j,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_74 = 0.250000 * scalar_37;

            // energy += -0.250000 <k,l||c,d>_abab t2_abab(a,b,i,l) u2_abab(c,d,k,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_74;

            // energy += -0.250000 <k,l||d,c>_abab t2_abab(a,b,i,l) u2_abab(d,c,k,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_74;

            // energy += -0.250000 <k,l||d,c>_abab t2_abab(b,a,i,l) u2_abab(d,c,k,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_74;

            // energy += -0.250000 <k,l||c,d>_abab t2_abab(b,a,i,l) u2_abab(c,d,k,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_74;
        }

        if (include_u0_ && include_u2_) {

            // scalar_tmp_75 += 1.000000 m2_xabab_Loovv("I,j,i,b,a") sharedOps["abab_vvoo_221"]("b,a,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_75 = 0.500000 * scalar_38;

            // rm0 += 0.500000 dp_aa(k,j) u2_abab(b,a,k,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_75;

            // rm0 += 0.500000 dp_aa(k,j) u2_abab(a,b,k,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_75;

            // energy += 0.500000 dp_aa(k,j) u0 u2_abab(b,a,k,i) m2_abab(j,i,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_75 * u0;

            // energy += 0.500000 dp_aa(k,j) u0 u2_abab(a,b,k,i) m2_abab(j,i,a,b)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_75 * u0;
        }

        {

            // scalar_tmp_76 += 1.000000 sharedOps["abab_vvoo_410"]("b,a,i,j") l2_xabab_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_76 = 0.250000 * scalar_39;

            // energy += -0.250000 <l,k||d,c>_abab l2_abab(i,j,b,a) t2_abab(b,c,l,k) t2_abab(d,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_76;

            // energy += -0.250000 <k,l||d,c>_abab l2_abab(i,j,b,a) t2_abab(b,c,k,l) t2_abab(d,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_76;

            // energy += -0.250000 <k,l||d,c>_abab l2_abab(j,i,b,a) t2_abab(b,c,k,l) t2_abab(d,a,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_76;

            // energy += -0.250000 <l,k||d,c>_abab l2_abab(j,i,b,a) t2_abab(b,c,l,k) t2_abab(d,a,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_76;

            // scalar_tmp_77 += 1.000000 sharedOps["abab_vvoo_425"]("a,b,i,j") l2_xabab_Loovv("I,i,j,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_77 = 0.250000 * scalar_40;

            // energy += -0.250000 <k,l||d,c>_abab l2_abab(i,j,a,b) t2_abab(d,c,k,j) t2_abab(a,b,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_77;

            // energy += -0.250000 <k,l||d,c>_abab l2_abab(i,j,b,a) t2_abab(d,c,k,j) t2_abab(b,a,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_77;

            // energy += -0.250000 <k,l||c,d>_abab l2_abab(i,j,b,a) t2_abab(c,d,k,j) t2_abab(b,a,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_77;

            // energy += -0.250000 <k,l||c,d>_abab l2_abab(i,j,a,b) t2_abab(c,d,k,j) t2_abab(a,b,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_77;
        }

        if (include_u2_) {

            // scalar_tmp_78 += 1.000000 m2_xabab_Loovv("I,i,j,a,b") sharedOps["abab_vvoo_422"]("a,b,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_78 = 0.250000 * scalar_41;

            // energy += -0.250000 <k,l||d,c>_abab t2_abab(d,c,i,l) u2_abab(a,b,k,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_78;

            // energy += -0.250000 <k,l||d,c>_abab t2_abab(d,c,i,l) u2_abab(b,a,k,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_78;

            // energy += -0.250000 <k,l||c,d>_abab t2_abab(c,d,i,l) u2_abab(a,b,k,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_78;

            // energy += -0.250000 <k,l||c,d>_abab t2_abab(c,d,i,l) u2_abab(b,a,k,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_78;
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // scalar_tmp_79 += 1.000000 sharedOps["abab_vvoo_218"]("a,b,i,j") m2_xabab_Loovv("I,i,j,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_79 = 0.500000 * scalar_42;

            // rm0 += 0.500000 dp_bb(k,c) t2_abab(b,a,i,k) u1_bb(c,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_79;

            // rm0 += 0.500000 dp_bb(k,c) t2_abab(a,b,i,k) u1_bb(c,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_79;

            // energy += 0.500000 dp_bb(k,c) t2_abab(a,b,i,k) u0 u1_bb(c,j) m2_abab(i,j,a,b)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_79 * u0;

            // energy += 0.500000 dp_bb(k,c) t2_abab(b,a,i,k) u0 u1_bb(c,j) m2_abab(i,j,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_79 * u0;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_80 += 1.000000 m2_xabab_Loovv("I,j,i,a,b") sharedOps["abab_vvoo_110"]("a,b,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_80 = 0.250000 * scalar_43;

            // energy += 0.250000 <k,l||j,c>_abab t2_abab(b,a,k,l) u1_bb(c,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_80;

            // energy += 0.250000 <l,k||j,c>_abab t2_abab(b,a,l,k) u1_bb(c,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_80;

            // energy += 0.250000 <k,l||j,c>_abab t2_abab(a,b,k,l) u1_bb(c,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_80;

            // energy += 0.250000 <l,k||j,c>_abab t2_abab(a,b,l,k) u1_bb(c,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_80;
        }

        {

            // scalar_tmp_81 += 1.000000 sharedOps["abab_vvoo_426"]("b,a,j,i") l2_xabab_Loovv("I,j,i,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_81 = 0.250000 * scalar_44;

            // energy += -0.250000 <l,k||c,d>_abab l2_abab(j,i,b,a) t2_abab(c,d,j,k) t2_abab(b,a,l,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_81;

            // energy += -0.250000 <l,k||d,c>_abab l2_abab(j,i,a,b) t2_abab(d,c,j,k) t2_abab(a,b,l,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_81;

            // energy += -0.250000 <l,k||c,d>_abab l2_abab(j,i,a,b) t2_abab(c,d,j,k) t2_abab(a,b,l,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_81;

            // energy += -0.250000 <l,k||d,c>_abab l2_abab(j,i,b,a) t2_abab(d,c,j,k) t2_abab(b,a,l,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_81;
        }

        if (include_u0_) {

            // scalar_tmp_82 += 1.000000 l2_xabab_Loovv("I,i,j,a,b") sharedOps["abab_vvoo_158"]("a,b,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_82 = 0.500000 * scalar_45;

            // rm0 += -0.500000 dp_bb(b,c) l2_abab(i,j,a,b) t2_abab(a,c,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_82;

            // rm0 += -0.500000 dp_bb(b,c) l2_abab(j,i,a,b) t2_abab(a,c,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_82;

            // energy += -0.500000 dp_bb(b,c) l2_abab(j,i,a,b) t2_abab(a,c,j,i) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_82 * u0;

            // energy += -0.500000 dp_bb(b,c) l2_abab(i,j,a,b) t2_abab(a,c,i,j) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_82 * u0;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_83 += 1.000000 m2_xabab_Loovv("I,i,j,a,b") sharedOps["abab_vvoo_109"]("a,b,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_83 = 0.250000 * scalar_46;

            // energy += 0.250000 <k,l||c,j>_abab t2_abab(b,a,k,l) u1_aa(c,i) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_83;

            // energy += 0.250000 <k,l||c,j>_abab t2_abab(a,b,k,l) u1_aa(c,i) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_83;

            // energy += 0.250000 <l,k||c,j>_abab t2_abab(b,a,l,k) u1_aa(c,i) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_83;

            // energy += 0.250000 <l,k||c,j>_abab t2_abab(a,b,l,k) u1_aa(c,i) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_83;
        }

        if (include_u0_ && include_u2_) {

            // scalar_tmp_84 += 1.000000 sharedOps["abab_vvoo_157"]("b,a,j,i") m2_xabab_Loovv("I,j,i,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_84 = 0.500000 * scalar_47;

            // rm0 += -0.500000 dp_aa(b,c) u2_abab(c,a,i,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_84;

            // rm0 += -0.500000 dp_aa(b,c) u2_abab(c,a,j,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_84;

            // energy += -0.500000 dp_aa(b,c) u0 u2_abab(c,a,j,i) m2_abab(j,i,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_84 * u0;

            // energy += -0.500000 dp_aa(b,c) u0 u2_abab(c,a,i,j) m2_abab(i,j,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_84 * u0;
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // scalar_tmp_85 += 1.000000 sharedOps["abab_vvoo_210"]("b,a,j,i") m2_xabab_Loovv("I,j,i,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_85 = 0.500000 * scalar_48;

            // rm0 += 0.500000 dp_aa(k,c) t2_abab(c,a,j,i) u1_aa(b,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_85;

            // rm0 += 0.500000 dp_aa(k,c) t2_abab(c,a,i,j) u1_aa(b,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_85;

            // energy += 0.500000 dp_aa(k,c) t2_abab(c,a,j,i) u0 u1_aa(b,k) m2_abab(j,i,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_85 * u0;

            // energy += 0.500000 dp_aa(k,c) t2_abab(c,a,i,j) u0 u1_aa(b,k) m2_abab(i,j,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_85 * u0;

            // scalar_tmp_86 += 1.000000 m2_xabab_Loovv("I,j,i,a,b") sharedOps["abab_vvoo_209"]("a,b,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_86 = 0.500000 * scalar_49;

            // rm0 += 0.500000 dp_bb(k,c) t2_abab(a,c,i,j) u1_bb(b,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_86;

            // rm0 += 0.500000 dp_bb(k,c) t2_abab(a,c,j,i) u1_bb(b,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_86;

            // energy += 0.500000 dp_bb(k,c) t2_abab(a,c,j,i) u0 u1_bb(b,k) m2_abab(j,i,a,b)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_86 * u0;

            // energy += 0.500000 dp_bb(k,c) t2_abab(a,c,i,j) u0 u1_bb(b,k) m2_abab(i,j,a,b)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_86 * u0;
        }

        if (include_u2_) {

            // scalar_tmp_87 += 1.000000 m2_xabab_Loovv("I,j,i,a,b") sharedOps["abab_vvoo_424"]("a,b,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_87 = 0.250000 * scalar_50;

            // energy += -0.250000 <l,k||c,d>_abab t2_abab(c,d,l,i) u2_abab(a,b,j,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_87;

            // energy += -0.250000 <l,k||d,c>_abab t2_abab(d,c,l,i) u2_abab(a,b,j,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_87;

            // energy += -0.250000 <l,k||c,d>_abab t2_abab(c,d,l,i) u2_abab(b,a,j,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_87;

            // energy += -0.250000 <l,k||d,c>_abab t2_abab(d,c,l,i) u2_abab(b,a,j,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_87;
        }

        if (include_u0_) {

            // scalar_tmp_88 += 1.000000 sharedOps["abab_vvoo_222"]("a,b,i,j") l2_xabab_Loovv("I,i,j,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_88 = 0.500000 * scalar_51;

            // rm0 += 0.500000 dp_bb(k,j) l2_abab(i,j,b,a) t2_abab(b,a,i,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_88;

            // rm0 += 0.500000 dp_bb(k,j) l2_abab(i,j,a,b) t2_abab(a,b,i,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_88;

            // energy += 0.500000 dp_bb(k,j) l2_abab(i,j,b,a) t2_abab(b,a,i,k) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_88 * u0;

            // energy += 0.500000 dp_bb(k,j) l2_abab(i,j,a,b) t2_abab(a,b,i,k) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_88 * u0;
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // scalar_tmp_89 += 1.000000 m2_xabab_Loovv("I,j,i,b,a") sharedOps["abab_vvoo_217"]("b,a,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_89 = 0.500000 * scalar_52;

            // rm0 += 0.500000 dp_aa(k,c) t2_abab(b,a,k,i) u1_aa(c,j) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_89;

            // rm0 += 0.500000 dp_aa(k,c) t2_abab(a,b,k,i) u1_aa(c,j) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_89;

            // energy += 0.500000 dp_aa(k,c) t2_abab(a,b,k,i) u0 u1_aa(c,j) m2_abab(j,i,a,b)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_89 * u0;

            // energy += 0.500000 dp_aa(k,c) t2_abab(b,a,k,i) u0 u1_aa(c,j) m2_abab(j,i,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_89 * u0;
        }

        if (include_u0_ && include_u2_) {

            // scalar_tmp_90 += 1.000000 m2_xabab_Loovv("I,i,j,b,a") sharedOps["abab_vvoo_219"]("b,a,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_90 = 0.500000 * scalar_53;

            // rm0 += 0.500000 dp_bb(k,j) u2_abab(b,a,i,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_90;

            // rm0 += 0.500000 dp_bb(k,j) u2_abab(a,b,i,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_90;

            // energy += 0.500000 dp_bb(k,j) u0 u2_abab(a,b,i,k) m2_abab(i,j,a,b)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_90 * u0;

            // energy += 0.500000 dp_bb(k,j) u0 u2_abab(b,a,i,k) m2_abab(i,j,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_90 * u0;
        }

        {

            // scalar_tmp_91 += 1.000000 V_blks_["abab_vvoo"]("a,b,j,i") l2_xabab_Loovv("I,j,i,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_91 = 0.250000 * scalar_54;

            // energy += 0.250000 <b,a||j,i>_abab l2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_91;

            // energy += 0.250000 <a,b||i,j>_abab l2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_91;

            // energy += 0.250000 <b,a||i,j>_abab l2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_91;

            // energy += 0.250000 <a,b||j,i>_abab l2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_91;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_92 += 1.000000 m2_xabab_Loovv("I,i,j,a,b") sharedOps["abab_vvoo_427"]("a,b,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_92 = 0.250000 * scalar_55;

            // energy += -0.250000 <k,b||c,d>_abab t2_abab(c,d,i,j) u1_aa(a,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_92;

            // energy += -0.250000 <k,b||c,d>_abab t2_abab(c,d,j,i) u1_aa(a,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_92;

            // energy += -0.250000 <k,b||d,c>_abab t2_abab(d,c,j,i) u1_aa(a,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_92;

            // energy += -0.250000 <k,b||d,c>_abab t2_abab(d,c,i,j) u1_aa(a,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_92;
        }

        if (include_u2_) {

            // scalar_tmp_93 += 1.000000 m2_xabab_Loovv("I,j,i,b,a") sharedOps["abab_vvoo_428"]("b,a,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_93 = 0.250000 * scalar_56;

            // energy += -0.250000 <l,k||d,c>_abab t2_abab(a,b,l,i) u2_abab(d,c,j,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_93;

            // energy += -0.250000 <l,k||c,d>_abab t2_abab(b,a,l,i) u2_abab(c,d,j,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_93;

            // energy += -0.250000 <l,k||d,c>_abab t2_abab(b,a,l,i) u2_abab(d,c,j,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_93;

            // energy += -0.250000 <l,k||c,d>_abab t2_abab(a,b,l,i) u2_abab(c,d,j,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_93;

            // scalar_tmp_94 += 1.000000 m2_xabab_Loovv("I,i,j,b,a") sharedOps["abab_vvoo_408"]("b,a,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_94 = 0.250000 * scalar_57;

            // energy += -0.250000 <l,k||d,c>_abab t2_abab(d,a,i,j) u2_abab(b,c,l,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_94;

            // energy += -0.250000 <l,k||d,c>_abab t2_abab(d,a,j,i) u2_abab(b,c,l,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_94;

            // energy += -0.250000 <k,l||d,c>_abab t2_abab(d,a,j,i) u2_abab(b,c,k,l) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_94;

            // energy += -0.250000 <k,l||d,c>_abab t2_abab(d,a,i,j) u2_abab(b,c,k,l) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_94;

            // scalar_tmp_95 += 1.000000 m2_xabab_Loovv("I,j,i,b,a") sharedOps["abab_vvoo_409"]("b,a,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_95 = 0.250000 * scalar_58;

            // energy += -0.250000 <k,l||d,c>_abab t2_abab(d,a,k,l) u2_abab(b,c,j,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_95;

            // energy += -0.250000 <l,k||d,c>_abab t2_abab(d,a,l,k) u2_abab(b,c,j,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_95;

            // energy += -0.250000 <k,l||d,c>_abab t2_abab(d,a,k,l) u2_abab(b,c,i,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_95;

            // energy += -0.250000 <l,k||d,c>_abab t2_abab(d,a,l,k) u2_abab(b,c,i,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_95;
        }

        if (include_u1_) {

            // tempOps["bb_vo_96"] += 1.000000 l2_xbbbb_Loovv("I,j,i,e,a") u1_bb_vo("a,j")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["bb_vo_96"]("e,i") = l2_bbbb_oovv("j,i,e,a") * u1_bb_vo("a,j");

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) dp_bb(n,e) l2_bbbb(i,m,f,a) u1_bb(a,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= tempOps["bb_vo_96"]("f,m") * dp_bb_ov("n,e");

            // rl2_abab_vvoo += 1.000000 dp_aa(m,e) l2_bbbb(i,n,f,a) u1_bb(a,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += tempOps["bb_vo_96"]("f,n") * dp_aa_ov("m,e");

            // rl1_bb_vo += -1.000000 dp_bb(m,i) l2_bbbb(j,i,e,a) u1_bb(a,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bb_vo_96"]("e,i") * dp_bb_oo("m,i");

            // rl1_bb_vo += 1.000000 dp_bb(a,e) l2_bbbb(i,m,a,b) u1_bb(b,i)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["bb_vo_96"]("a,m") * dp_bb_vv("a,e");

            // tempOps["bb_vo_97"] += 1.000000 u1_aa_vo("a,j") l2_xabab_Loovv("I,j,i,a,e")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["bb_vo_97"]("e,i") = u1_aa_vo("a,j") * l2_abab_oovv("j,i,a,e");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) dp_bb(n,e) l2_abab(i,m,a,f) u1_aa(a,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += tempOps["bb_vo_97"]("f,m") * dp_bb_ov("n,e");

            // rl2_abab_vvoo += -1.000000 dp_aa(m,e) l2_abab(i,n,a,f) u1_aa(a,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["bb_vo_97"]("f,n") * dp_aa_ov("m,e");

            // rl1_bb_vo += -1.000000 dp_bb(a,e) l2_abab(i,m,b,a) u1_aa(b,i)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bb_vo_97"]("a,m") * dp_bb_vv("a,e");

            // rl1_bb_vo += 1.000000 dp_bb(m,i) l2_abab(j,i,a,e) u1_aa(a,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["bb_vo_97"]("e,i") * dp_bb_oo("m,i");

            // tempOps["aa_vo_98"] += 1.000000 l2_xabab_Loovv("I,i,j,e,a") u1_bb_vo("a,j")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["aa_vo_98"]("e,i") = l2_abab_oovv("i,j,e,a") * u1_bb_vo("a,j");

            // rl2_abab_vvoo += -1.000000 dp_bb(n,f) l2_abab(m,i,e,a) u1_bb(a,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["aa_vo_98"]("e,m") * dp_bb_ov("n,f");

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) dp_aa(n,e) l2_abab(m,i,f,a) u1_bb(a,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += tempOps["aa_vo_98"]("f,m") * dp_aa_ov("n,e");

            // rl1_aa_vo += 1.000000 dp_aa(m,i) l2_abab(i,j,e,a) u1_bb(a,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["aa_vo_98"]("e,i") * dp_aa_oo("m,i");

            // rl1_aa_vo += -1.000000 dp_aa(a,e) l2_abab(m,i,a,b) u1_bb(b,i)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aa_vo_98"]("a,m") * dp_aa_vv("a,e");

            // tempOps["aa_vo_99"] += 1.000000 l2_xaaaa_Loovv("I,i,m,a,b") u1_aa_vo("b,i")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["aa_vo_99"]("a,m") = l2_aaaa_oovv("i,m,a,b") * u1_aa_vo("b,i");

            // rl2_abab_vvoo += 1.000000 dp_bb(n,f) l2_aaaa(i,m,e,a) u1_aa(a,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += tempOps["aa_vo_99"]("e,m") * dp_bb_ov("n,f");

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) dp_aa(n,e) l2_aaaa(i,m,f,a) u1_aa(a,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= tempOps["aa_vo_99"]("f,m") * dp_aa_ov("n,e");

            // rl1_aa_vo += 1.000000 dp_aa(a,e) l2_aaaa(i,m,a,b) u1_aa(b,i)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["aa_vo_99"]("a,m") * dp_aa_vv("a,e");

            // rl1_aa_vo += -1.000000 dp_aa(m,i) l2_aaaa(j,i,e,a) u1_aa(a,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aa_vo_99"]("e,i") * dp_aa_oo("m,i");

            // scalar_tmp_100 += 1.000000 sharedOps["abab_vvoo_209"]("a,b,i,j") l2_xabab_Loovv("I,i,j,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_100 = 0.500000 * scalar_59;

            // energy += 0.500000 dp_bb(k,c) l2_abab(i,j,a,b) t2_abab(a,c,i,j) u1_bb(b,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_100;

            // energy += 0.500000 dp_bb(k,c) l2_abab(j,i,a,b) t2_abab(a,c,j,i) u1_bb(b,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_100;
        }
        tempOps["bb_vo_96"] = TArrayD();
        tempOps["bb_vo_97"] = TArrayD();
        tempOps["aa_vo_98"] = TArrayD();
        tempOps["aa_vo_99"] = TArrayD();

        if (include_u2_) {

            // scalar_tmp_101 += 1.000000 m2_xabab_Loovv("I,j,i,a,b") sharedOps["abab_vvoo_158"]("a,b,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_101 = 0.500000 * scalar_60;

            // energy += -0.500000 dp_bb(b,c) t2_abab(a,c,j,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_101;

            // energy += -0.500000 dp_bb(b,c) t2_abab(a,c,i,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_101;

            // scalar_tmp_102 += 1.000000 sharedOps["abab_vvoo_254"]("b,a,i,j") m2_xabab_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_102 = 0.500000 * scalar_61;

            // energy += -0.500000 f_bb(k,j) u2_abab(b,a,i,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_102;

            // energy += -0.500000 f_bb(k,j) u2_abab(a,b,i,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_102;

            // scalar_tmp_103 += 1.000000 sharedOps["abab_vvoo_156"]("c,a,i,j") l2_xabab_Loovv("I,i,j,c,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_103 = 0.500000 * scalar_62;

            // energy += -0.500000 dp_bb(a,b) l2_abab(j,i,c,a) u2_abab(c,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_103;

            // energy += -0.500000 dp_bb(a,b) l2_abab(i,j,c,a) u2_abab(c,b,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_103;
        }

        if (include_u0_) {

            // scalar_tmp_104 += 1.000000 l2_xbbbb_Loovv("I,i,j,b,a") sharedOps["bbbb_vvoo_238"]("b,a,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_104 = 0.500000 * scalar_63;

            // rm0 += 0.500000 dp_bb(k,j) l2_bbbb(i,j,b,a) t2_bbbb(b,a,i,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_104;

            // energy += 0.500000 dp_bb(k,j) l2_bbbb(i,j,b,a) t2_bbbb(b,a,i,k) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_104 * u0;
        }

        if (include_u2_) {

            // scalar_tmp_105 += 1.000000 m2_xabab_Loovv("I,j,i,b,a") sharedOps["abab_vvoo_155"]("b,a,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_105 = 0.500000 * scalar_64;

            // energy += -0.500000 dp_aa(b,c) t2_abab(c,a,i,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_105;

            // energy += -0.500000 dp_aa(b,c) t2_abab(c,a,j,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_105;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_106 += 1.000000 sharedOps["abab_vvoo_175"]("b,a,j,i") m2_xabab_Loovv("I,j,i,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_106 = 0.500000 * scalar_65;

            // energy += 0.500000 <k,b||c,d>_aaaa t2_abab(d,a,j,i) u1_aa(c,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_106;

            // energy += 0.500000 <k,b||c,d>_aaaa t2_abab(d,a,i,j) u1_aa(c,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_106;
        }

        {

            // scalar_tmp_107 += 1.000000 l2_xabab_Loovv("I,i,j,b,a") sharedOps["abab_vvoo_370"]("b,a,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_107 = 0.500000 * scalar_66;

            // energy += 0.500000 <k,l||d,c>_abab l2_abab(i,j,b,a) t2_abab(b,c,k,j) t2_abab(d,a,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_107;

            // energy += 0.500000 <l,k||c,d>_abab l2_abab(j,i,a,b) t2_abab(c,b,j,k) t2_abab(a,d,l,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_107;
        }

        if (include_u2_) {

            // scalar_tmp_108 += 1.000000 sharedOps["abab_vvoo_220"]("a,b,j,i") m2_xabab_Loovv("I,j,i,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_108 = 0.500000 * scalar_67;

            // energy += 0.500000 dp_aa(k,j) t2_abab(b,a,k,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_108;

            // energy += 0.500000 dp_aa(k,j) t2_abab(a,b,k,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_108;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_109 += 1.000000 m2_xabab_Loovv("I,j,i,a,b") sharedOps["abab_vvoo_249"]("a,b,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_109 = scalar_68;

            // energy += 1.000000 dp_aa(k,c) u1_aa(c,j) u2_abab(b,a,k,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_109;

            // energy += 1.000000 dp_aa(k,c) u1_aa(c,j) u2_abab(a,b,k,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_109;

            // scalar_tmp_110 += 1.000000 sharedOps["abab_vvoo_247"]("a,b,j,i") m2_xabab_Loovv("I,j,i,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_110 = 0.500000 * scalar_69;

            // energy += 0.500000 <l,k||c,j>_aaaa t2_abab(b,a,l,i) u1_aa(c,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_110;

            // energy += 0.500000 <l,k||c,j>_aaaa t2_abab(a,b,l,i) u1_aa(c,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_110;
        }

        if (include_u2_) {

            // scalar_tmp_111 += 1.000000 sharedOps["aaaa_vvoo_412"]("b,a,i,j") m2_xaaaa_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_111 = 0.250000 * scalar_70;

            // energy += 0.250000 <k,l||c,d>_abab t2_abab(a,d,k,l) u2_aaaa(c,b,i,j) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_111;

            // energy += 0.250000 <l,k||c,d>_abab t2_abab(a,d,l,k) u2_aaaa(c,b,i,j) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_111;

            // scalar_tmp_112 += 1.000000 m2_xabab_Loovv("I,j,i,a,b") sharedOps["abab_vvoo_180"]("a,b,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_112 = 0.500000 * scalar_71;

            // energy += 0.500000 f_bb(b,c) u2_abab(a,c,i,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_112;

            // energy += 0.500000 f_bb(b,c) u2_abab(a,c,j,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_112;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_113 += 1.000000 sharedOps["abab_vvoo_246"]("b,a,i,j") m2_xabab_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_113 = 0.500000 * scalar_72;

            // energy += 0.500000 <l,k||c,j>_bbbb t2_abab(a,b,i,l) u1_bb(c,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_113;

            // energy += 0.500000 <l,k||c,j>_bbbb t2_abab(b,a,i,l) u1_bb(c,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_113;
        }

        if (include_u2_) {

            // scalar_tmp_114 += 1.000000 sharedOps["abab_vvoo_222"]("a,b,i,j") m2_xabab_Loovv("I,i,j,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_114 = 0.500000 * scalar_73;

            // energy += 0.500000 dp_bb(k,j) t2_abab(b,a,i,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_114;

            // energy += 0.500000 dp_bb(k,j) t2_abab(a,b,i,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_114;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_115 += 1.000000 m2_xabab_Loovv("I,i,j,a,b") sharedOps["abab_vvoo_245"]("a,b,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_115 = 0.500000 * scalar_74;

            // energy += -0.500000 <k,l||c,j>_abab t2_abab(b,a,i,l) u1_aa(c,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_115;

            // energy += -0.500000 <k,l||c,j>_abab t2_abab(a,b,i,l) u1_aa(c,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_115;

            // scalar_tmp_116 += 1.000000 sharedOps["abab_vvoo_190"]("b,a,j,i") m2_xabab_Loovv("I,j,i,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_116 = 0.500000 * scalar_75;

            // energy += 0.500000 <a,b||j,c>_abab u1_bb(c,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_116;

            // energy += 0.500000 <b,a||j,c>_abab u1_bb(c,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_116;
        }

        if (include_u2_) {

            // scalar_tmp_117 += 1.000000 sharedOps["bbbb_vvoo_440"]("b,a,j,i") m2_xbbbb_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_117 = 0.250000 * scalar_76;

            // energy += 0.250000 <l,k||d,c>_abab t2_abab(d,c,l,i) u2_bbbb(b,a,j,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_117;

            // energy += 0.250000 <l,k||c,d>_abab t2_abab(c,d,l,i) u2_bbbb(b,a,j,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_117;
        }

        {

            // scalar_tmp_118 += 1.000000 sharedOps["abab_vvoo_189"]("a,b,i,j") l2_xabab_Loovv("I,i,j,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_118 = 0.500000 * scalar_77;

            // energy += 0.500000 f_bb(b,c) l2_abab(j,i,a,b) t2_abab(a,c,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_118;

            // energy += 0.500000 f_bb(b,c) l2_abab(i,j,a,b) t2_abab(a,c,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_118;
        }

        if (include_u0_ && include_u2_) {

            // scalar_tmp_119 += 1.000000 m2_xbbbb_Loovv("I,i,j,b,a") sharedOps["bbbb_vvoo_170"]("b,a,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_119 = 0.500000 * scalar_78;

            // rm0 += -0.500000 dp_bb(b,c) u2_bbbb(c,a,i,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_119;

            // energy += -0.500000 dp_bb(b,c) u0 u2_bbbb(c,a,i,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_119 * u0;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_120 += 1.000000 m2_xabab_Loovv("I,i,j,b,a") sharedOps["abab_vvoo_252"]("b,a,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_120 = scalar_79;

            // energy += 1.000000 dp_bb(k,c) u1_bb(c,j) u2_abab(b,a,i,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_120;

            // energy += 1.000000 dp_bb(k,c) u1_bb(c,j) u2_abab(a,b,i,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_120;
        }

        if (include_u1_) {

            // scalar_tmp_121 += 1.000000 l2_xabab_Loovv("I,i,j,b,a") sharedOps["abab_vvoo_218"]("b,a,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_121 = 0.500000 * scalar_80;

            // energy += 0.500000 dp_bb(k,c) l2_abab(i,j,b,a) t2_abab(b,a,i,k) u1_bb(c,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_121;

            // energy += 0.500000 dp_bb(k,c) l2_abab(i,j,a,b) t2_abab(a,b,i,k) u1_bb(c,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_121;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_122 += 1.000000 sharedOps["abab_vvoo_176"]("b,a,i,j") m2_xabab_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_122 = 0.500000 * scalar_81;

            // energy += 0.500000 <b,k||d,c>_abab t2_abab(d,a,j,i) u1_bb(c,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_122;

            // energy += 0.500000 <b,k||d,c>_abab t2_abab(d,a,i,j) u1_bb(c,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_122;

            // scalar_tmp_123 += 1.000000 sharedOps["abab_vvoo_174"]("a,b,i,j") m2_xabab_Loovv("I,i,j,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_123 = 0.500000 * scalar_82;

            // energy += 0.500000 <k,b||c,d>_abab t2_abab(a,d,j,i) u1_aa(c,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_123;

            // energy += 0.500000 <k,b||c,d>_abab t2_abab(a,d,i,j) u1_aa(c,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_123;

            // scalar_tmp_124 += 1.000000 sharedOps["abab_vvoo_260"]("a,b,j,i") m2_xabab_Loovv("I,j,i,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_124 = 0.500000 * scalar_83;

            // energy += -0.500000 <k,b||i,j>_abab u1_aa(a,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_124;

            // energy += -0.500000 <k,b||j,i>_abab u1_aa(a,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_124;
        }

        {

            // scalar_tmp_125 += 1.000000 sharedOps["abab_vvoo_371"]("b,a,j,i") l2_xabab_Loovv("I,j,i,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_125 = 0.500000 * scalar_84;

            // energy += 0.500000 <l,k||d,c>_abab l2_abab(j,i,b,a) t2_abab(b,c,j,k) t2_abab(d,a,l,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_125;

            // energy += 0.500000 <k,l||c,d>_abab l2_abab(i,j,a,b) t2_abab(c,b,k,j) t2_abab(a,d,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_125;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_126 += 1.000000 sharedOps["abab_vvoo_173"]("a,b,j,i") m2_xabab_Loovv("I,j,i,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_126 = 0.500000 * scalar_85;

            // energy += 0.500000 <k,b||c,d>_bbbb t2_abab(a,d,i,j) u1_bb(c,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_126;

            // energy += 0.500000 <k,b||c,d>_bbbb t2_abab(a,d,j,i) u1_bb(c,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_126;
        }

        {

            // scalar_tmp_127 += 1.000000 l2_xabab_Loovv("I,j,i,a,b") sharedOps["abab_vvoo_162"]("a,b,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_127 = 0.250000 * scalar_86;

            // energy += -0.250000 <l,k||c,d>_bbbb l2_abab(i,j,a,b) t2_bbbb(c,b,l,k) t2_abab(a,d,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_127;

            // energy += -0.250000 <l,k||c,d>_bbbb l2_abab(j,i,a,b) t2_bbbb(c,b,l,k) t2_abab(a,d,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_127;
        }

        if (include_u2_) {

            // scalar_tmp_128 += 1.000000 sharedOps["abab_vvoo_258"]("b,a,j,i") m2_xabab_Loovv("I,j,i,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_128 = 0.500000 * scalar_87;

            // energy += -0.500000 f_aa(k,j) u2_abab(a,b,k,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_128;

            // energy += -0.500000 f_aa(k,j) u2_abab(b,a,k,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_128;
        }

        if (include_u0_ && include_u2_) {

            // scalar_tmp_129 += 1.000000 sharedOps["aaaa_vvoo_237"]("b,a,j,i") m2_xaaaa_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_129 = 0.500000 * scalar_88;

            // rm0 += 0.500000 dp_aa(k,j) u2_aaaa(b,a,i,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_129;

            // energy += 0.500000 dp_aa(k,j) u0 u2_aaaa(b,a,i,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_129 * u0;
        }

        if (include_u2_) {

            // scalar_tmp_130 += 1.000000 m2_xabab_Loovv("I,j,i,b,a") sharedOps["abab_vvoo_182"]("b,a,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_130 = 0.500000 * scalar_89;

            // energy += 0.500000 f_aa(b,c) u2_abab(c,a,j,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_130;

            // energy += 0.500000 f_aa(b,c) u2_abab(c,a,i,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_130;
        }

        if (include_u0_) {

            // scalar_tmp_131 += 1.000000 sharedOps["aaaa_vvoo_172"]("b,a,i,j") l2_xaaaa_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_131 = 0.500000 * scalar_90;

            // rm0 += -0.500000 dp_aa(b,c) l2_aaaa(i,j,b,a) t2_aaaa(c,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_131;

            // energy += -0.500000 dp_aa(b,c) l2_aaaa(i,j,b,a) t2_aaaa(c,a,i,j) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_131 * u0;
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // scalar_tmp_132 += 1.000000 sharedOps["bbbb_vvoo_236"]("b,a,j,i") m2_xbbbb_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_132 = 0.500000 * scalar_91;

            // rm0 += 0.500000 dp_bb(k,c) t2_bbbb(b,a,i,k) u1_bb(c,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_132;

            // energy += 0.500000 dp_bb(k,c) t2_bbbb(b,a,i,k) u0 u1_bb(c,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_132 * u0;
        }

        if (include_u0_ && include_u2_) {

            // scalar_tmp_133 += 1.000000 sharedOps["bbbb_vvoo_241"]("b,a,i,j") m2_xbbbb_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_133 = 0.500000 * scalar_92;

            // rm0 += 0.500000 dp_bb(k,j) u2_bbbb(b,a,i,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_133;

            // energy += 0.500000 dp_bb(k,j) u0 u2_bbbb(b,a,i,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_133 * u0;
        }

        if (include_u1_) {

            // scalar_tmp_134 += 1.000000 l2_xabab_Loovv("I,i,j,b,a") sharedOps["abab_vvoo_210"]("b,a,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_134 = 0.500000 * scalar_93;

            // energy += 0.500000 dp_aa(k,c) l2_abab(i,j,b,a) t2_abab(c,a,i,j) u1_aa(b,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_134;

            // energy += 0.500000 dp_aa(k,c) l2_abab(j,i,b,a) t2_abab(c,a,j,i) u1_aa(b,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_134;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_135 += 1.000000 sharedOps["abab_vvoo_248"]("b,a,j,i") m2_xabab_Loovv("I,j,i,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_135 = 0.500000 * scalar_94;

            // energy += -0.500000 <l,k||j,c>_abab t2_abab(a,b,l,i) u1_bb(c,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_135;

            // energy += -0.500000 <l,k||j,c>_abab t2_abab(b,a,l,i) u1_bb(c,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_135;
        }

        if (include_u0_) {

            // scalar_tmp_136 += 1.000000 sharedOps["aaaa_vvoo_239"]("b,a,i,j") l2_xaaaa_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_136 = 0.500000 * scalar_95;

            // rm0 += 0.500000 dp_aa(k,j) l2_aaaa(i,j,b,a) t2_aaaa(b,a,i,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_136;

            // energy += 0.500000 dp_aa(k,j) l2_aaaa(i,j,b,a) t2_aaaa(b,a,i,k) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_136 * u0;
        }

        {

            // scalar_tmp_137 += 1.000000 l2_xabab_Loovv("I,i,j,a,b") sharedOps["abab_vvoo_233"]("a,b,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_137 = 0.250000 * scalar_96;

            // energy += -0.250000 <l,k||c,d>_bbbb l2_abab(i,j,b,a) t2_bbbb(c,d,j,k) t2_abab(b,a,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_137;

            // energy += -0.250000 <l,k||c,d>_bbbb l2_abab(i,j,a,b) t2_bbbb(c,d,j,k) t2_abab(a,b,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_137;

            // scalar_tmp_138 += 1.000000 l2_xabab_Loovv("I,i,j,b,a") sharedOps["abab_vvoo_269"]("b,a,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_138 = 0.500000 * scalar_97;

            // energy += -0.500000 f_bb(k,j) l2_abab(i,j,b,a) t2_abab(b,a,i,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_138;

            // energy += -0.500000 f_bb(k,j) l2_abab(i,j,a,b) t2_abab(a,b,i,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_138;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_139 += 1.000000 m2_xabab_Loovv("I,i,j,b,a") sharedOps["abab_vvoo_441"]("b,a,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_139 = scalar_98;

            // energy += 1.000000 dp_aa(k,c) u1_aa(b,k) u2_abab(c,a,i,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_139;

            // energy += 1.000000 dp_aa(k,c) u1_aa(b,k) u2_abab(c,a,j,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_139;
        }

        if (include_u2_) {

            // scalar_tmp_140 += 1.000000 m2_xaaaa_Loovv("I,i,j,b,a") sharedOps["aaaa_vvoo_443"]("b,a,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_140 = 0.250000 * scalar_99;

            // energy += -0.250000 <l,k||d,c>_abab t2_aaaa(b,a,i,l) u2_abab(d,c,j,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_140;

            // energy += -0.250000 <l,k||c,d>_abab t2_aaaa(b,a,i,l) u2_abab(c,d,j,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_140;

            // scalar_tmp_141 += 1.000000 m2_xabab_Loovv("I,i,j,b,a") sharedOps["abab_vvoo_229"]("b,a,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_141 = 0.250000 * scalar_100;

            // energy += -0.250000 <l,k||c,d>_bbbb t2_abab(b,a,i,l) u2_bbbb(c,d,j,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_141;

            // energy += -0.250000 <l,k||c,d>_bbbb t2_abab(a,b,i,l) u2_bbbb(c,d,j,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_141;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_142 += 1.000000 sharedOps["abab_vvoo_230"]("a,b,j,i") m2_xabab_Loovv("I,j,i,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_142 = 0.500000 * scalar_101;

            // energy += -0.500000 f_bb(k,c) t2_abab(a,c,i,j) u1_bb(b,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_142;

            // energy += -0.500000 f_bb(k,c) t2_abab(a,c,j,i) u1_bb(b,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_142;
        }

        if (include_u2_) {

            // scalar_tmp_143 += 1.000000 m2_xabab_Loovv("I,j,i,b,a") sharedOps["abab_vvoo_159"]("b,a,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_143 = 0.250000 * scalar_102;

            // energy += -0.250000 <l,k||c,d>_aaaa t2_abab(d,a,j,i) u2_aaaa(c,b,l,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_143;

            // energy += -0.250000 <l,k||c,d>_aaaa t2_abab(d,a,i,j) u2_aaaa(c,b,l,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_143;
        }

        if (include_u0_ && include_u2_) {

            // scalar_tmp_144 += 1.000000 m2_xaaaa_Loovv("I,i,j,b,a") sharedOps["aaaa_vvoo_164"]("b,a,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_144 = 0.500000 * scalar_103;

            // rm0 += -0.500000 dp_aa(b,c) u2_aaaa(c,a,i,j) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_144;

            // energy += -0.500000 dp_aa(b,c) u0 u2_aaaa(c,a,i,j) m2_aaaa(i,j,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_144 * u0;
        }

        {

            // scalar_tmp_145 += 1.000000 sharedOps["abab_vvoo_231"]("b,a,j,i") l2_xabab_Loovv("I,j,i,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_145 = 0.250000 * scalar_104;

            // energy += -0.250000 <l,k||c,d>_aaaa l2_abab(j,i,a,b) t2_aaaa(c,d,j,k) t2_abab(a,b,l,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_145;

            // energy += -0.250000 <l,k||c,d>_aaaa l2_abab(j,i,b,a) t2_aaaa(c,d,j,k) t2_abab(b,a,l,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_145;

            // scalar_tmp_146 += 1.000000 sharedOps["abab_vvoo_181"]("b,a,i,j") l2_xabab_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_146 = 0.500000 * scalar_105;

            // energy += 0.500000 f_aa(b,c) l2_abab(j,i,b,a) t2_abab(c,a,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_146;

            // energy += 0.500000 f_aa(b,c) l2_abab(i,j,b,a) t2_abab(c,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_146;
        }

        if (include_u2_) {

            // scalar_tmp_147 += 1.000000 l2_xabab_Loovv("I,k,i,a,b") sharedOps["abab_vvoo_219"]("a,b,k,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_147 = 0.500000 * scalar_106;

            // energy += 0.500000 dp_bb(j,i) l2_abab(k,i,a,b) u2_abab(a,b,k,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_147;

            // energy += 0.500000 dp_bb(j,i) l2_abab(k,i,b,a) u2_abab(b,a,k,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_147;

            // scalar_tmp_148 += 1.000000 sharedOps["abab_vvoo_157"]("a,c,i,j") l2_xabab_Loovv("I,i,j,a,c")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_148 = 0.500000 * scalar_107;

            // energy += -0.500000 dp_aa(a,b) l2_abab(j,i,a,c) u2_abab(b,c,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_148;

            // energy += -0.500000 dp_aa(a,b) l2_abab(i,j,a,c) u2_abab(b,c,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_148;
        }

        if (include_u1_) {

            // scalar_tmp_149 += 1.000000 sharedOps["abab_vvoo_217"]("a,b,j,i") l2_xabab_Loovv("I,j,i,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_149 = 0.500000 * scalar_108;

            // energy += 0.500000 dp_aa(k,c) l2_abab(j,i,b,a) t2_abab(b,a,k,i) u1_aa(c,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_149;

            // energy += 0.500000 dp_aa(k,c) l2_abab(j,i,a,b) t2_abab(a,b,k,i) u1_aa(c,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_149;
        }

        if (include_u0_) {

            // scalar_tmp_150 += 1.000000 l2_xbbbb_Loovv("I,i,j,b,a") sharedOps["bbbb_vvoo_168"]("b,a,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_150 = 0.500000 * scalar_109;

            // rm0 += -0.500000 dp_bb(b,c) l2_bbbb(i,j,b,a) t2_bbbb(c,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_150;

            // energy += -0.500000 dp_bb(b,c) l2_bbbb(i,j,b,a) t2_bbbb(c,a,i,j) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_150 * u0;
        }

        {

            // scalar_tmp_151 += 1.000000 sharedOps["abab_vvoo_268"]("a,b,j,i") l2_xabab_Loovv("I,j,i,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_151 = 0.500000 * scalar_110;

            // energy += -0.500000 f_aa(k,j) l2_abab(j,i,b,a) t2_abab(b,a,k,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_151;

            // energy += -0.500000 f_aa(k,j) l2_abab(j,i,a,b) t2_abab(a,b,k,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_151;
        }

        if (include_u2_) {

            // scalar_tmp_152 += 1.000000 l2_xabab_Loovv("I,i,k,a,b") sharedOps["abab_vvoo_221"]("a,b,i,k")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_152 = 0.500000 * scalar_111;

            // energy += 0.500000 dp_aa(j,i) l2_abab(i,k,b,a) u2_abab(b,a,j,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_152;

            // energy += 0.500000 dp_aa(j,i) l2_abab(i,k,a,b) u2_abab(a,b,j,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_152;

            // scalar_tmp_153 += 1.000000 m2_xabab_Loovv("I,j,i,a,b") sharedOps["abab_vvoo_414"]("a,b,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_153 = 0.250000 * scalar_112;

            // energy += 0.250000 <l,k||c,d>_aaaa t2_aaaa(d,a,l,k) u2_abab(c,b,i,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_153;

            // energy += 0.250000 <l,k||c,d>_aaaa t2_aaaa(d,a,l,k) u2_abab(c,b,j,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_153;

            // scalar_tmp_154 += 1.000000 sharedOps["abab_vvoo_160"]("a,b,i,j") m2_xabab_Loovv("I,i,j,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_154 = 0.250000 * scalar_113;

            // energy += -0.250000 <l,k||c,d>_bbbb t2_abab(a,d,i,j) u2_bbbb(c,b,l,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_154;

            // energy += -0.250000 <l,k||c,d>_bbbb t2_abab(a,d,j,i) u2_bbbb(c,b,l,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_154;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_155 += 1.000000 sharedOps["abab_vvoo_250"]("b,a,i,j") m2_xabab_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_155 = 0.500000 * scalar_114;

            // energy += -0.500000 f_bb(k,c) t2_abab(a,b,i,k) u1_bb(c,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_155;

            // energy += -0.500000 f_bb(k,c) t2_abab(b,a,i,k) u1_bb(c,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_155;
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // scalar_tmp_156 += 1.000000 m2_xaaaa_Loovv("I,i,j,b,a") sharedOps["aaaa_vvoo_235"]("b,a,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_156 = 0.500000 * scalar_115;

            // rm0 += 0.500000 dp_aa(k,c) t2_aaaa(b,a,i,k) u1_aa(c,j) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_156;

            // energy += 0.500000 dp_aa(k,c) t2_aaaa(b,a,i,k) u0 u1_aa(c,j) m2_aaaa(i,j,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_156 * u0;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_157 += 1.000000 m2_xabab_Loovv("I,j,i,b,a") sharedOps["abab_vvoo_251"]("b,a,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_157 = 0.500000 * scalar_116;

            // energy += -0.500000 f_aa(k,c) t2_abab(b,a,k,i) u1_aa(c,j) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_157;

            // energy += -0.500000 f_aa(k,c) t2_abab(a,b,k,i) u1_aa(c,j) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_157;
        }

        {

            // scalar_tmp_158 += 1.000000 sharedOps["bbbb_vvoo_442"]("b,a,i,j") l2_xbbbb_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_158 = 0.250000 * scalar_117;

            // energy += -0.250000 <k,l||d,c>_abab l2_bbbb(i,j,b,a) t2_abab(d,c,k,j) t2_bbbb(b,a,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_158;

            // energy += -0.250000 <k,l||c,d>_abab l2_bbbb(i,j,b,a) t2_abab(c,d,k,j) t2_bbbb(b,a,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_158;

            // scalar_tmp_159 += 1.000000 l2_xabab_Loovv("I,j,i,b,a") sharedOps["abab_vvoo_161"]("b,a,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_159 = 0.250000 * scalar_118;

            // energy += -0.250000 <l,k||c,d>_aaaa l2_abab(j,i,b,a) t2_aaaa(c,b,l,k) t2_abab(d,a,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_159;

            // energy += -0.250000 <l,k||c,d>_aaaa l2_abab(i,j,b,a) t2_aaaa(c,b,l,k) t2_abab(d,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_159;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_160 += 1.000000 sharedOps["abab_vvoo_225"]("b,a,j,i") m2_xabab_Loovv("I,j,i,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_160 = 0.500000 * scalar_119;

            // energy += -0.500000 f_aa(k,c) t2_abab(c,a,j,i) u1_aa(b,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_160;

            // energy += -0.500000 f_aa(k,c) t2_abab(c,a,i,j) u1_aa(b,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_160;

            // scalar_tmp_161 += 1.000000 sharedOps["abab_vvoo_263"]("b,a,i,j") m2_xabab_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_161 = 0.500000 * scalar_120;

            // energy += -0.500000 <b,k||i,j>_abab u1_bb(a,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_161;

            // energy += -0.500000 <b,k||j,i>_abab u1_bb(a,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_161;
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // scalar_tmp_162 += 1.000000 m2_xbbbb_Loovv("I,i,j,b,a") sharedOps["bbbb_vvoo_430"]("a,b,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_162 = 0.500000 * scalar_121;

            // rm0 += 0.500000 dp_bb(k,c) t2_bbbb(c,a,i,j) u1_bb(b,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_162;

            // energy += 0.500000 dp_bb(k,c) t2_bbbb(c,a,i,j) u0 u1_bb(b,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_162 * u0;
        }

        if (include_u2_) {

            // scalar_tmp_163 += 1.000000 sharedOps["abab_vvoo_434"]("a,b,i,j") m2_xabab_Loovv("I,i,j,a,b")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_163 = 0.250000 * scalar_122;

            // energy += 0.250000 <l,k||c,d>_aaaa t2_aaaa(c,d,i,l) u2_abab(a,b,k,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_163;

            // energy += 0.250000 <l,k||c,d>_aaaa t2_aaaa(c,d,i,l) u2_abab(b,a,k,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_163;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_164 += 1.000000 sharedOps["abab_vvoo_185"]("b,a,i,j") m2_xabab_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_164 = 0.500000 * scalar_123;

            // energy += 0.500000 <a,b||c,j>_abab u1_aa(c,i) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_164;

            // energy += 0.500000 <b,a||c,j>_abab u1_aa(c,i) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_164;
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // scalar_tmp_165 += 1.000000 sharedOps["aaaa_vvoo_431"]("b,a,i,j") m2_xaaaa_Loovv("I,i,j,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_165 = 0.500000 * scalar_124;

            // rm0 += 0.500000 dp_aa(k,c) t2_aaaa(c,a,i,j) u1_aa(b,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_165;

            // energy += 0.500000 dp_aa(k,c) t2_aaaa(c,a,i,j) u0 u1_aa(b,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_165 * u0;
        }

        if (include_u2_) {

            // scalar_tmp_166 += 1.000000 m2_xabab_Loovv("I,j,i,a,b") sharedOps["abab_vvoo_433"]("a,b,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_166 = 0.250000 * scalar_125;

            // energy += 0.250000 <l,k||c,d>_bbbb t2_bbbb(c,d,i,l) u2_abab(a,b,j,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_166;

            // energy += 0.250000 <l,k||c,d>_bbbb t2_bbbb(c,d,i,l) u2_abab(b,a,j,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_166;

            // scalar_tmp_167 += 1.000000 m2_xbbbb_Loovv("I,i,j,b,a") sharedOps["bbbb_vvoo_419"]("b,a,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_167 = 0.250000 * scalar_126;

            // energy += -0.250000 <k,l||c,d>_abab t2_bbbb(d,a,i,j) u2_abab(c,b,k,l) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_167;

            // energy += -0.250000 <l,k||c,d>_abab t2_bbbb(d,a,i,j) u2_abab(c,b,l,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_167;

            // scalar_tmp_168 += 1.000000 m2_xaaaa_Loovv("I,i,j,b,a") sharedOps["aaaa_vvoo_413"]("a,b,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_168 = 0.250000 * scalar_127;

            // energy += -0.250000 <k,l||d,c>_abab t2_aaaa(d,a,i,j) u2_abab(b,c,k,l) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_168;

            // energy += -0.250000 <l,k||d,c>_abab t2_aaaa(d,a,i,j) u2_abab(b,c,l,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_168;
        }

        {

            // scalar_tmp_169 += 1.000000 sharedOps["abab_vvoo_52"]("b,a,j,i") l2_xabab_Loovv("I,j,i,b,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_169 = 0.500000 * scalar_128;

            // energy += 0.500000 <k,l||c,d>_abab l2_abab(j,i,b,a) t2_aaaa(c,b,j,k) t2_bbbb(d,a,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_169;

            // energy += 0.500000 <l,k||d,c>_abab l2_abab(i,j,a,b) t2_bbbb(c,b,j,k) t2_aaaa(d,a,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_169;
        }

        if (include_u2_) {

            // scalar_tmp_170 += 1.000000 m2_xbbbb_Loovv("I,i,j,b,a") sharedOps["bbbb_vvoo_415"]("b,a,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_170 = 0.250000 * scalar_129;

            // energy += 0.250000 <l,k||d,c>_abab t2_abab(d,a,l,k) u2_bbbb(c,b,i,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_170;

            // energy += 0.250000 <k,l||d,c>_abab t2_abab(d,a,k,l) u2_bbbb(c,b,i,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_170;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_171 += 1.000000 m2_xabab_Loovv("I,i,j,a,b") sharedOps["abab_vvoo_432"]("a,b,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_171 = scalar_130;

            // energy += 1.000000 dp_bb(k,c) u1_bb(b,k) u2_abab(a,c,j,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_171;

            // energy += 1.000000 dp_bb(k,c) u1_bb(b,k) u2_abab(a,c,i,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_171;
        }

        if (include_u2_) {

            // scalar_tmp_172 += 1.000000 m2_xabab_Loovv("I,j,i,b,a") sharedOps["abab_vvoo_416"]("b,a,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_172 = 0.250000 * scalar_131;

            // energy += 0.250000 <l,k||c,d>_bbbb t2_bbbb(d,a,l,k) u2_abab(b,c,j,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_172;

            // energy += 0.250000 <l,k||c,d>_bbbb t2_bbbb(d,a,l,k) u2_abab(b,c,i,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_172;

            // scalar_tmp_173 += 1.000000 m2_xabab_Loovv("I,j,i,a,b") sharedOps["abab_vvoo_227"]("a,b,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_173 = 0.250000 * scalar_132;

            // energy += -0.250000 <l,k||c,d>_aaaa t2_abab(b,a,l,i) u2_aaaa(c,d,j,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_173;

            // energy += -0.250000 <l,k||c,d>_aaaa t2_abab(a,b,l,i) u2_aaaa(c,d,j,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_173;

            // scalar_tmp_174 += 1.000000 m2_xbbbb_Loovv("I,i,j,b,a") sharedOps["bbbb_vvoo_438"]("b,a,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_174 = 0.250000 * scalar_133;

            // energy += -0.250000 <k,l||d,c>_abab t2_bbbb(b,a,i,l) u2_abab(d,c,k,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_174;

            // energy += -0.250000 <k,l||c,d>_abab t2_bbbb(b,a,i,l) u2_abab(c,d,k,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_174;

            // scalar_tmp_175 += 1.000000 m2_xaaaa_Loovv("I,i,j,b,a") sharedOps["aaaa_vvoo_437"]("b,a,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_175 = 0.250000 * scalar_134;

            // energy += 0.250000 <k,l||d,c>_abab t2_abab(d,c,i,l) u2_aaaa(b,a,j,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_175;

            // energy += 0.250000 <k,l||c,d>_abab t2_abab(c,d,i,l) u2_aaaa(b,a,j,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_175;
        }

        {

            // scalar_tmp_176 += 1.000000 l2_xaaaa_Loovv("I,i,j,b,a") sharedOps["aaaa_vvoo_439"]("b,a,j,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_176 = 0.250000 * scalar_135;

            // energy += -0.250000 <l,k||d,c>_abab l2_aaaa(i,j,b,a) t2_abab(d,c,j,k) t2_aaaa(b,a,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_176;

            // energy += -0.250000 <l,k||c,d>_abab l2_aaaa(i,j,b,a) t2_abab(c,d,j,k) t2_aaaa(b,a,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_176;

            // scalar_tmp_177 += 1.000000 l2_xbbbb_Loovv("I,i,j,b,a") sharedOps["bbbb_vvoo_417"]("b,a,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_177 = 0.250000 * scalar_136;

            // energy += -0.250000 <l,k||c,d>_abab l2_bbbb(i,j,b,a) t2_abab(c,b,l,k) t2_bbbb(d,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_177;

            // energy += -0.250000 <k,l||c,d>_abab l2_bbbb(i,j,b,a) t2_abab(c,b,k,l) t2_bbbb(d,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_177;

            // scalar_tmp_178 += 1.000000 l2_xaaaa_Loovv("I,i,j,b,a") sharedOps["aaaa_vvoo_418"]("b,a,i,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_178 = 0.250000 * scalar_137;

            // energy += -0.250000 <k,l||d,c>_abab l2_aaaa(i,j,b,a) t2_abab(b,c,k,l) t2_aaaa(d,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_178;

            // energy += -0.250000 <l,k||d,c>_abab l2_aaaa(i,j,b,a) t2_abab(b,c,l,k) t2_aaaa(d,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_178;
        }

        if (include_u1_) {

            // tempOps["bb_vo_179"] += 1.000000 m1_xaa_Lov("I,i,a") t2_abab_vvoo("a,b,i,j")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["bb_vo_179"]("b,j") = m1_aa_ov("i,a") * t2_abab_vvoo("a,b,i,j");

            // rm1_bb_vo += -1.000000 <j,m||e,b>_bbbb t2_abab(a,b,i,j) m1_aa(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= tempOps["bb_vo_179"]("b,j") * V_blks_["bbbb_oovv"]("j,m,e,b");

            // rm1_aa_vo += 1.000000 <m,j||e,b>_abab t2_abab(a,b,i,j) m1_aa(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += tempOps["bb_vo_179"]("b,j") * V_blks_["abab_oovv"]("m,j,e,b");

            // tempOps["bb_vo_180"] += 1.000000 t2_bbbb_vvoo("b,a,i,j") m1_xbb_Lov("I,i,a")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["bb_vo_180"]("b,j") = t2_bbbb_vvoo("b,a,i,j") * m1_bb_ov("i,a");

            // rm1_bb_vo += 1.000000 <j,m||e,b>_bbbb t2_bbbb(b,a,i,j) m1_bb(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += tempOps["bb_vo_180"]("b,j") * V_blks_["bbbb_oovv"]("j,m,e,b");

            // rm1_aa_vo += -1.000000 <m,j||e,b>_abab t2_bbbb(b,a,i,j) m1_bb(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= tempOps["bb_vo_180"]("b,j") * V_blks_["abab_oovv"]("m,j,e,b");

            // tempOps["aa_vo_181"] += 1.000000 t2_aaaa_vvoo("b,a,i,j") m1_xaa_Lov("I,i,a")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["aa_vo_181"]("b,j") = t2_aaaa_vvoo("b,a,i,j") * m1_aa_ov("i,a");

            // rm1_bb_vo += -1.000000 <j,m||b,e>_abab t2_aaaa(b,a,i,j) m1_aa(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= tempOps["aa_vo_181"]("b,j") * V_blks_["abab_oovv"]("j,m,b,e");

            // rm1_aa_vo += 1.000000 <j,m||e,b>_aaaa t2_aaaa(b,a,i,j) m1_aa(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += tempOps["aa_vo_181"]("b,j") * V_blks_["aaaa_oovv"]("j,m,e,b");

            // tempOps["aa_vo_182"] += 1.000000 m1_xbb_Lov("I,i,a") t2_abab_vvoo("b,a,j,i")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["aa_vo_182"]("b,j") = m1_bb_ov("i,a") * t2_abab_vvoo("b,a,j,i");

            // rm1_bb_vo += 1.000000 <j,m||b,e>_abab t2_abab(b,a,j,i) m1_bb(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += tempOps["aa_vo_182"]("b,j") * V_blks_["abab_oovv"]("j,m,b,e");

            // rm1_aa_vo += -1.000000 <j,m||e,b>_aaaa t2_abab(b,a,j,i) m1_bb(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= tempOps["aa_vo_182"]("b,j") * V_blks_["aaaa_oovv"]("j,m,e,b");
        }
        tempOps["bb_vo_179"] = TArrayD();
        tempOps["bb_vo_180"] = TArrayD();
        tempOps["aa_vo_181"] = TArrayD();
        tempOps["aa_vo_182"] = TArrayD();

        {

            // tempOps["aa_vo_183"] += 1.000000 t2_aaaa_vvoo("b,a,i,j") l1_xaa_Lov("I,i,a")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["aa_vo_183"]("b,j") = t2_aaaa_vvoo("b,a,i,j") * l1_aa_ov("i,a");

            // rl1_bb_vo += -1.000000 <j,m||b,e>_abab l1_aa(i,a) t2_aaaa(b,a,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["aa_vo_183"]("b,j") * V_blks_["abab_oovv"]("j,m,b,e");

            // rl1_aa_vo += 1.000000 <j,m||e,b>_aaaa l1_aa(i,a) t2_aaaa(b,a,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["aa_vo_183"]("b,j") * V_blks_["aaaa_oovv"]("j,m,e,b");
            tempOps["aa_vo_183"] = TArrayD();
        }

        if (include_u1_ && include_u2_) {

            // tempOps["bb_vo_184"] += 1.000000 u2_abab_vvoo("a,b,i,j") m1_xaa_Lov("I,i,a")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["bb_vo_184"]("b,j") = u2_abab_vvoo("a,b,i,j") * m1_aa_ov("i,a");

            // rl1_bb_vo += -1.000000 <j,m||e,b>_bbbb u2_abab(a,b,i,j) m1_aa(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bb_vo_184"]("b,j") * V_blks_["bbbb_oovv"]("j,m,e,b");

            // rl1_aa_vo += 1.000000 <m,j||e,b>_abab u2_abab(a,b,i,j) m1_aa(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["bb_vo_184"]("b,j") * V_blks_["abab_oovv"]("m,j,e,b");

            // tempOps["aa_vo_185"] += 1.000000 u2_aaaa_vvoo("b,a,i,j") m1_xaa_Lov("I,i,a")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["aa_vo_185"]("b,j") = u2_aaaa_vvoo("b,a,i,j") * m1_aa_ov("i,a");

            // rl1_bb_vo += -1.000000 <j,m||b,e>_abab u2_aaaa(b,a,i,j) m1_aa(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["aa_vo_185"]("b,j") * V_blks_["abab_oovv"]("j,m,b,e");

            // rl1_aa_vo += 1.000000 <j,m||e,b>_aaaa u2_aaaa(b,a,i,j) m1_aa(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["aa_vo_185"]("b,j") * V_blks_["aaaa_oovv"]("j,m,e,b");

            // tempOps["bb_vo_186"] += 1.000000 u2_bbbb_vvoo("b,a,i,j") m1_xbb_Lov("I,i,a")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["bb_vo_186"]("b,j") = u2_bbbb_vvoo("b,a,i,j") * m1_bb_ov("i,a");

            // rl1_bb_vo += 1.000000 <j,m||e,b>_bbbb u2_bbbb(b,a,i,j) m1_bb(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["bb_vo_186"]("b,j") * V_blks_["bbbb_oovv"]("j,m,e,b");

            // rl1_aa_vo += -1.000000 <m,j||e,b>_abab u2_bbbb(b,a,i,j) m1_bb(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["bb_vo_186"]("b,j") * V_blks_["abab_oovv"]("m,j,e,b");
        }
        tempOps["bb_vo_184"] = TArrayD();
        tempOps["aa_vo_185"] = TArrayD();
        tempOps["bb_vo_186"] = TArrayD();

        {

            // tempOps["bb_vo_187"] += 1.000000 l1_xaa_Lov("I,i,a") t2_abab_vvoo("a,b,i,j")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["bb_vo_187"]("b,j") = l1_aa_ov("i,a") * t2_abab_vvoo("a,b,i,j");

            // rl1_bb_vo += -1.000000 <j,m||e,b>_bbbb l1_aa(i,a) t2_abab(a,b,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bb_vo_187"]("b,j") * V_blks_["bbbb_oovv"]("j,m,e,b");

            // rl1_aa_vo += 1.000000 <m,j||e,b>_abab l1_aa(i,a) t2_abab(a,b,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["bb_vo_187"]("b,j") * V_blks_["abab_oovv"]("m,j,e,b");
            tempOps["bb_vo_187"] = TArrayD();

            // tempOps["aa_vo_188"] += 1.000000 t2_abab_vvoo("b,a,j,i") l1_xbb_Lov("I,i,a")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["aa_vo_188"]("b,j") = t2_abab_vvoo("b,a,j,i") * l1_bb_ov("i,a");

            // rl1_bb_vo += 1.000000 <j,m||b,e>_abab l1_bb(i,a) t2_abab(b,a,j,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["aa_vo_188"]("b,j") * V_blks_["abab_oovv"]("j,m,b,e");

            // rl1_aa_vo += -1.000000 <j,m||e,b>_aaaa l1_bb(i,a) t2_abab(b,a,j,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aa_vo_188"]("b,j") * V_blks_["aaaa_oovv"]("j,m,e,b");
            tempOps["aa_vo_188"] = TArrayD();

            // tempOps["bb_vo_189"] += 1.000000 l1_xbb_Lov("I,i,a") t2_bbbb_vvoo("b,a,i,j")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["bb_vo_189"]("b,j") = l1_bb_ov("i,a") * t2_bbbb_vvoo("b,a,i,j");

            // rl1_bb_vo += 1.000000 <j,m||e,b>_bbbb l1_bb(i,a) t2_bbbb(b,a,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["bb_vo_189"]("b,j") * V_blks_["bbbb_oovv"]("j,m,e,b");

            // rl1_aa_vo += -1.000000 <m,j||e,b>_abab l1_bb(i,a) t2_bbbb(b,a,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["bb_vo_189"]("b,j") * V_blks_["abab_oovv"]("m,j,e,b");
            tempOps["bb_vo_189"] = TArrayD();
        }

        if (include_u1_ && include_u2_) {

            // tempOps["aa_vo_190"] += 1.000000 m1_xbb_Lov("I,i,a") u2_abab_vvoo("b,a,j,i")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["aa_vo_190"]("b,j") = m1_bb_ov("i,a") * u2_abab_vvoo("b,a,j,i");

            // rl1_bb_vo += 1.000000 <j,m||b,e>_abab u2_abab(b,a,j,i) m1_bb(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["aa_vo_190"]("b,j") * V_blks_["abab_oovv"]("j,m,b,e");

            // rl1_aa_vo += -1.000000 <j,m||e,b>_aaaa u2_abab(b,a,j,i) m1_bb(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aa_vo_190"]("b,j") * V_blks_["aaaa_oovv"]("j,m,e,b");
        }
        tempOps["aa_vo_190"] = TArrayD();

        if (include_u1_) {

            // tempOps["bb_vo_191"] += 1.000000 dp_bb_vv("a,e") m1_xbb_Lov("I,m,a")
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["bb_vo_191"]("e,m") = dp_bb_vv("a,e") * m1_bb_ov("m,a");
        }

        if (include_u0_ && include_u1_) {

            // rm1_bb_vo += -1.000000 dp_bb(a,e) u0 m1_bb(m,a)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rm1_bb_vo("e,m") -= tempOps["bb_vo_191"]("e,m") * u0;
        }

        if (include_u1_) {

            // rl1_bb_vo += -1.000000 dp_bb(a,e) m1_bb(m,a)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rl1_bb_vo("e,m") -= tempOps["bb_vo_191"]("e,m");
        }
        tempOps["bb_vo_191"] = TArrayD();

        {

            // tempOps["bb_vo_192"] += 1.000000 l1_xbb_Lov("I,m,a") dp_bb_vv("a,e")
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["bb_vo_192"]("e,m") = l1_bb_ov("m,a") * dp_bb_vv("a,e");
        }

        if (include_u1_) {

            // rm1_bb_vo += -1.000000 dp_bb(a,e) l1_bb(m,a)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rm1_bb_vo("e,m") -= tempOps["bb_vo_192"]("e,m");
        }

        if (include_u0_) {

            // rl1_bb_vo += -1.000000 dp_bb(a,e) l1_bb(m,a) u0
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rl1_bb_vo("e,m") -= tempOps["bb_vo_192"]("e,m") * u0;
        }
        tempOps["bb_vo_192"] = TArrayD();

        {

            // tempOps["aa_vo_193"] += 1.000000 dp_aa_vv("a,e") l1_xaa_Lov("I,m,a")
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["aa_vo_193"]("e,m") = dp_aa_vv("a,e") * l1_aa_ov("m,a");
        }

        if (include_u1_) {

            // rm1_aa_vo += -1.000000 dp_aa(a,e) l1_aa(m,a)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rm1_aa_vo("e,m") -= tempOps["aa_vo_193"]("e,m");
        }

        if (include_u0_) {

            // rl1_aa_vo += -1.000000 dp_aa(a,e) l1_aa(m,a) u0
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rl1_aa_vo("e,m") -= tempOps["aa_vo_193"]("e,m") * u0;
        }
        tempOps["aa_vo_193"] = TArrayD();

        if (include_u1_) {

            // tempOps["aa_vo_194"] += 1.000000 dp_aa_vv("a,e") m1_xaa_Lov("I,m,a")
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["aa_vo_194"]("e,m") = dp_aa_vv("a,e") * m1_aa_ov("m,a");
        }

        if (include_u0_ && include_u1_) {

            // rm1_aa_vo += -1.000000 dp_aa(a,e) u0 m1_aa(m,a)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rm1_aa_vo("e,m") -= tempOps["aa_vo_194"]("e,m") * u0;
        }

        if (include_u1_) {

            // rl1_aa_vo += -1.000000 dp_aa(a,e) m1_aa(m,a)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rl1_aa_vo("e,m") -= tempOps["aa_vo_194"]("e,m");

            // tempOps["aa_oo_195"] += 1.000000 u1_aa_vo("a,j") m1_xaa_Lov("I,i,a")
            // flops: o2v1L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempOps["aa_oo_195"]("j,i") = u1_aa_vo("a,j") * m1_aa_ov("i,a");

            // rm1_aa_vo += 2.000000 dp_aa(i,e) u1_aa(a,i) m1_aa(m,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += 2.000000 * tempOps["aa_oo_195"]("i,m") * dp_aa_ov("i,e");

            // rl2_abab_vvoo += -1.000000 <i,n||e,f>_abab u1_aa(a,i) m1_aa(m,a)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["aa_oo_195"]("i,m") * V_blks_["abab_oovv"]("i,n,e,f");
        }
        tempOps["aa_vo_194"] = TArrayD();

        {

            // rl2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("f,e,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) <i,n||e,f>_aaaa u1_aa(a,i) m1_aa(m,a)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempOps["aa_oo_195"]("i,m") * V_blks_["aaaa_oovv"]("i,n,e,f");

            // rl1_bb_vo += -1.000000 <j,m||i,e>_abab u1_aa(a,j) m1_aa(i,a)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["aa_oo_195"]("j,i") * V_blks_["abba_oovo"]("j,m,e,i");

            // rl1_aa_vo += 1.000000 <j,m||e,i>_aaaa u1_aa(a,j) m1_aa(i,a)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += tempOps["aa_oo_195"]("j,i") * V_blks_["aaaa_oovo"]("j,m,e,i");
        }

        if (include_u0_ && include_u1_) {

            // rl1_aa_vo += 1.000000 dp_aa(i,e) u0 u1_aa(a,i) m1_aa(m,a)
            // flops: o2v1L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_aa_vo("e,m") += u0 * dp_aa_ov("i,e") * tempOps["aa_oo_195"]("i,m");
        }

        if (include_u1_) {

            // rl1_aa_vo += -1.000000 f_aa(i,e) u1_aa(a,i) m1_aa(m,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["aa_oo_195"]("i,m") * F_blks_["aa_ov"]("i,e");

            // tempOps["bb_oo_196"] += 1.000000 m1_xbb_Lov("I,i,a") u1_bb_vo("a,j")
            // flops: o2v1L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempOps["bb_oo_196"]("i,j") = m1_bb_ov("i,a") * u1_bb_vo("a,j");

            // rm1_bb_vo += 2.000000 dp_bb(i,e) u1_bb(a,i) m1_bb(m,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += 2.000000 * tempOps["bb_oo_196"]("m,i") * dp_bb_ov("i,e");
        }
        tempOps["aa_oo_195"] = TArrayD();

        {

            // rl2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("f,e,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) <i,n||e,f>_bbbb u1_bb(a,i) m1_bb(m,a)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempOps["bb_oo_196"]("m,i") * V_blks_["bbbb_oovv"]("i,n,e,f");

            // rl2_abab_vvoo += -1.000000 <m,i||e,f>_abab u1_bb(a,i) m1_bb(n,a)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= tempOps["bb_oo_196"]("n,i") * V_blks_["abab_oovv"]("m,i,e,f");
        }

        if (include_u0_ && include_u1_) {

            // rl1_bb_vo += 1.000000 dp_bb(i,e) u0 u1_bb(a,i) m1_bb(m,a)
            // flops: o2v1L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_bb_vo("e,m") += u0 * dp_bb_ov("i,e") * tempOps["bb_oo_196"]("m,i");
        }

        if (include_u1_) {

            // rl1_bb_vo += -1.000000 f_bb(i,e) u1_bb(a,i) m1_bb(m,a)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= tempOps["bb_oo_196"]("m,i") * F_blks_["bb_ov"]("i,e");

            // rl1_bb_vo += 1.000000 <j,m||e,i>_bbbb u1_bb(a,j) m1_bb(i,a)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += tempOps["bb_oo_196"]("i,j") * V_blks_["bbbb_oovo"]("j,m,e,i");

            // rl1_aa_vo += -1.000000 <m,j||e,i>_abab u1_bb(a,j) m1_bb(i,a)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= tempOps["bb_oo_196"]("i,j") * V_blks_["abab_oovo"]("m,j,e,i");
        }
        tempOps["bb_oo_196"] = TArrayD();

        {

            // tempOps["bb_vo_197"] += 1.000000 l1_xbb_Lov("I,i,e") dp_bb_oo("m,i")
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["bb_vo_197"]("e,m") = l1_bb_ov("i,e") * dp_bb_oo("m,i");
        }

        if (include_u1_) {

            // rm1_bb_vo += 1.000000 dp_bb(m,i) l1_bb(i,e)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rm1_bb_vo("e,m") += tempOps["bb_vo_197"]("e,m");
        }

        if (include_u0_) {

            // rl1_bb_vo += 1.000000 dp_bb(m,i) l1_bb(i,e) u0
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rl1_bb_vo("e,m") += tempOps["bb_vo_197"]("e,m") * u0;
        }
        tempOps["bb_vo_197"] = TArrayD();

        if (include_u1_) {

            // tempOps["bb_vo_198"] += 1.000000 m1_xbb_Lov("I,i,e") dp_bb_oo("m,i")
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["bb_vo_198"]("e,m") = m1_bb_ov("i,e") * dp_bb_oo("m,i");
        }

        if (include_u0_ && include_u1_) {

            // rm1_bb_vo += 1.000000 dp_bb(m,i) u0 m1_bb(i,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rm1_bb_vo("e,m") += tempOps["bb_vo_198"]("e,m") * u0;
        }

        if (include_u1_) {

            // rl1_bb_vo += 1.000000 dp_bb(m,i) m1_bb(i,e)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rl1_bb_vo("e,m") += tempOps["bb_vo_198"]("e,m");
        }
        tempOps["bb_vo_198"] = TArrayD();

        {

            // tempOps["aa_vo_199"] += 1.000000 l1_xaa_Lov("I,i,e") dp_aa_oo("m,i")
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["aa_vo_199"]("e,m") = l1_aa_ov("i,e") * dp_aa_oo("m,i");
        }

        if (include_u1_) {

            // rm1_aa_vo += 1.000000 dp_aa(m,i) l1_aa(i,e)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rm1_aa_vo("e,m") += tempOps["aa_vo_199"]("e,m");
        }

        if (include_u0_) {

            // rl1_aa_vo += 1.000000 dp_aa(m,i) l1_aa(i,e) u0
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rl1_aa_vo("e,m") += tempOps["aa_vo_199"]("e,m") * u0;
        }
        tempOps["aa_vo_199"] = TArrayD();

        if (include_u1_) {

            // tempOps["aa_vo_200"] += 1.000000 dp_aa_oo("m,i") m1_xaa_Lov("I,i,e")
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempOps["aa_vo_200"]("e,m") = dp_aa_oo("m,i") * m1_aa_ov("i,e");
        }

        if (include_u0_ && include_u1_) {

            // rm1_aa_vo += 1.000000 dp_aa(m,i) u0 m1_aa(i,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rm1_aa_vo("e,m") += tempOps["aa_vo_200"]("e,m") * u0;
        }

        if (include_u1_) {

            // rl1_aa_vo += 1.000000 dp_aa(m,i) m1_aa(i,e)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rl1_aa_vo("e,m") += tempOps["aa_vo_200"]("e,m");

            // scalar_tmp_201 += 1.000000 u1_bb_vo("a,i") m1_xbb_Lov("I,i,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_201 = scalar_138;

            // rm1_bb_vo += -1.000000 dp_bb(m,e) u1_bb(a,i) m1_bb(i,a)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= scalar_tmp_201 * dp_bb_ov("m,e");

            // rm1_aa_vo += -1.000000 dp_aa(m,e) u1_bb(a,i) m1_bb(i,a)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= scalar_tmp_201 * dp_aa_ov("m,e");

            // energy += -1.000000 dp_aa(j,b) u1_aa(b,j) u1_bb(a,i) m1_bb(i,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_6 * scalar_tmp_201;

            // energy += -1.000000 dp_bb(j,b) u1_bb(b,j) u1_bb(a,i) m1_bb(i,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_7 * scalar_tmp_201;

            // energy += 1.000000 u1_bb(a,i) m1_bb(i,a) w0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_201 * w0;

            // scalar_tmp_202 += 1.000000 l1_xbb_Lov("I,j,a") u1_bb_vo("a,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_202 = scalar_139;

            // rl1_bb_vo += -1.000000 dp_bb(m,e) l1_bb(i,a) u1_bb(a,i)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= scalar_tmp_202 * dp_bb_ov("m,e");

            // rl1_aa_vo += -1.000000 dp_aa(m,e) l1_bb(i,a) u1_bb(a,i)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= scalar_tmp_202 * dp_aa_ov("m,e");

            // cenergy += 1.000000 l1_bb(i,a) u1_bb(a,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            cenergy += scalar_tmp_202;

            // energy += -1.000000 dp_bb(i,i) l1_bb(j,a) u1_bb(a,j)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_5 * scalar_tmp_202;

            // energy += -1.000000 dp_aa(i,i) l1_bb(j,a) u1_bb(a,j)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_4 * scalar_tmp_202;

            // scalar_tmp_203 += 1.000000 l1_xaa_Lov("I,j,a") u1_aa_vo("a,j")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_203 = scalar_140;

            // rl1_bb_vo += -1.000000 dp_bb(m,e) l1_aa(i,a) u1_aa(a,i)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= scalar_tmp_203 * dp_bb_ov("m,e");

            // rl1_aa_vo += -1.000000 dp_aa(m,e) l1_aa(i,a) u1_aa(a,i)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= scalar_tmp_203 * dp_aa_ov("m,e");

            // cenergy += 1.000000 l1_aa(i,a) u1_aa(a,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            cenergy += scalar_tmp_203;

            // energy += -1.000000 dp_aa(i,i) l1_aa(j,a) u1_aa(a,j)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_4 * scalar_tmp_203;

            // energy += -1.000000 dp_bb(i,i) l1_aa(j,a) u1_aa(a,j)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_5 * scalar_tmp_203;

            // scalar_tmp_204 += 1.000000 m1_xaa_Lov("I,i,a") u1_aa_vo("a,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_204 = scalar_141;

            // rm1_bb_vo += -1.000000 dp_bb(m,e) u1_aa(a,i) m1_aa(i,a)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= scalar_tmp_204 * dp_bb_ov("m,e");

            // rm1_aa_vo += -1.000000 dp_aa(m,e) u1_aa(a,i) m1_aa(i,a)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= scalar_tmp_204 * dp_aa_ov("m,e");

            // energy += -1.000000 dp_bb(j,b) u1_bb(b,j) u1_aa(a,i) m1_aa(i,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_7 * scalar_tmp_204;

            // energy += 1.000000 u1_aa(a,i) m1_aa(i,a) w0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_204 * w0;

            // energy += -1.000000 dp_aa(j,b) u1_aa(b,j) u1_aa(a,i) m1_aa(i,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_6 * scalar_tmp_204;
        }
        tempOps["aa_vo_200"] = TArrayD();

        if (include_u0_ && include_u1_) {

            // scalar_tmp_205 += 1.000000 sharedOps["aa_vo_360"]("a,i") m1_xaa_Lov("I,i,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_205 = scalar_142;

            // rm0 += -1.000000 dp_aa(a,b) u1_aa(b,i) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_205;

            // energy += -1.000000 dp_aa(a,b) u0 u1_aa(b,i) m1_aa(i,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_205 * u0;
        }

        {

            // scalar_tmp_206 += 1.000000 sharedOps["bb_vo_266"]("a,i") l1_xbb_Lov("I,i,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_206 = 0.500000 * scalar_143;

            // energy += -0.500000 <k,j||b,i>_abab l1_bb(i,a) t2_abab(b,a,k,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_206;

            // energy += -0.500000 <j,k||b,i>_abab l1_bb(i,a) t2_abab(b,a,j,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_206;

            // scalar_tmp_207 += 1.000000 sharedOps["aa_vo_265"]("a,i") l1_xaa_Lov("I,i,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_207 = 0.500000 * scalar_144;

            // energy += -0.500000 <k,j||i,b>_abab l1_aa(i,a) t2_abab(a,b,k,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_207;

            // energy += -0.500000 <j,k||i,b>_abab l1_aa(i,a) t2_abab(a,b,j,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_207;

            // scalar_tmp_208 += 1.000000 sharedOps["aa_vo_184"]("a,i") l1_xaa_Lov("I,i,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_208 = 0.500000 * scalar_145;

            // energy += 0.500000 <a,j||b,c>_abab l1_aa(i,a) t2_abab(b,c,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_208;

            // energy += 0.500000 <a,j||c,b>_abab l1_aa(i,a) t2_abab(c,b,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_208;
        }

        if (include_u0_) {

            // scalar_tmp_209 += 1.000000 l1_xaa_Lov("I,i,a") sharedOps["aa_vo_318"]("a,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_209 = scalar_146;

            // rm0 += -1.000000 dp_bb(j,b) l1_aa(i,a) t2_abab(a,b,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_209;

            // energy += -1.000000 dp_bb(j,b) l1_aa(i,a) t2_abab(a,b,i,j) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_209 * u0;

            // scalar_tmp_210 += 1.000000 l1_xaa_Lov("I,i,a") sharedOps["aa_vo_320"]("a,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_210 = scalar_147;

            // rm0 += 1.000000 dp_aa(j,b) l1_aa(i,a) t2_aaaa(b,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_210;

            // energy += 1.000000 dp_aa(j,b) l1_aa(i,a) t2_aaaa(b,a,i,j) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_210 * u0;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_211 += 1.000000 sharedOps["bb_vo_188"]("a,i") m1_xbb_Lov("I,i,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_211 = 0.500000 * scalar_148;

            // energy += 0.500000 <j,a||c,b>_abab u2_abab(c,b,j,i) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_211;

            // energy += 0.500000 <j,a||b,c>_abab u2_abab(b,c,j,i) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_211;
        }

        if (include_u0_ && include_u1_) {

            // scalar_tmp_212 += 1.000000 sharedOps["bb_vo_361"]("a,i") m1_xbb_Lov("I,i,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_212 = scalar_149;

            // rm0 += -1.000000 dp_bb(a,b) u1_bb(b,i) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_212;

            // energy += -1.000000 dp_bb(a,b) u0 u1_bb(b,i) m1_bb(i,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_212 * u0;
        }

        if (include_u0_) {

            // scalar_tmp_213 += 1.000000 sharedOps["bb_vo_321"]("a,i") l1_xbb_Lov("I,i,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_213 = scalar_150;

            // rm0 += -1.000000 dp_aa(j,b) l1_bb(i,a) t2_abab(b,a,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_213;

            // energy += -1.000000 dp_aa(j,b) l1_bb(i,a) t2_abab(b,a,j,i) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_213 * u0;
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // scalar_tmp_214 += 1.000000 sharedOps["bb_vo_319"]("a,i") m1_xbb_Lov("I,i,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_214 = scalar_151;

            // rm0 += -1.000000 dp_aa(j,b) u2_abab(b,a,j,i) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_214;

            // energy += -1.000000 dp_aa(j,b) u0 u2_abab(b,a,j,i) m1_bb(i,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_214 * u0;

            // scalar_tmp_215 += 1.000000 m1_xaa_Lov("I,i,a") sharedOps["aa_vo_324"]("a,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_215 = scalar_152;

            // rm0 += 1.000000 dp_aa(j,b) u2_aaaa(b,a,i,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_215;

            // energy += 1.000000 dp_aa(j,b) u0 u2_aaaa(b,a,i,j) m1_aa(i,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_215 * u0;
        }

        if (include_u1_) {

            // scalar_tmp_216 += 1.000000 sharedOps["bb_vo_435"]("a,i") m1_xbb_Lov("I,i,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_216 = 0.500000 * scalar_153;

            // energy += -0.500000 <k,j||c,b>_abab t2_abab(c,a,k,j) u1_bb(b,i) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_216;

            // energy += -0.500000 <j,k||c,b>_abab t2_abab(c,a,j,k) u1_bb(b,i) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_216;
        }

        {

            // scalar_tmp_217 += 1.000000 l1_xbb_Lov("I,i,a") sharedOps["bb_vo_183"]("a,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_217 = 0.500000 * scalar_154;

            // energy += 0.500000 <j,a||c,b>_abab l1_bb(i,a) t2_abab(c,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_217;

            // energy += 0.500000 <j,a||b,c>_abab l1_bb(i,a) t2_abab(b,c,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_217;

            // scalar_tmp_218 += 1.000000 dp_aa_vo("a,i") l1_xaa_Lov("I,i,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_218 = scalar_155;
        }

        if (include_u0_) {

            // rm0 += -1.000000 dp_aa(a,i) l1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_218;

            // energy += -1.000000 dp_aa(a,i) l1_aa(i,a) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_218 * u0;
        }

        if (include_u1_) {

            // scalar_tmp_219 += 1.000000 sharedOps["bb_vo_499"]("a,i") m1_xbb_Lov("I,i,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_219 = 0.500000 * scalar_156;

            // energy += -0.500000 <k,j||b,c>_abab t2_abab(b,c,k,i) u1_bb(a,j) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_219;

            // energy += -0.500000 <k,j||c,b>_abab t2_abab(c,b,k,i) u1_bb(a,j) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_219;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_220 += 1.000000 m1_xaa_Lov("I,i,a") sharedOps["aa_vo_177"]("a,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_220 = 0.500000 * scalar_157;

            // energy += 0.500000 <a,j||b,c>_abab u2_abab(b,c,i,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_220;

            // energy += 0.500000 <a,j||c,b>_abab u2_abab(c,b,i,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_220;
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // scalar_tmp_221 += 1.000000 sharedOps["bb_vo_322"]("a,i") m1_xbb_Lov("I,i,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_221 = scalar_158;

            // rm0 += 1.000000 dp_bb(j,b) u2_bbbb(b,a,i,j) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_221;

            // energy += 1.000000 dp_bb(j,b) u0 u2_bbbb(b,a,i,j) m1_bb(i,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_221 * u0;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_222 += 1.000000 sharedOps["bb_vo_264"]("a,i") m1_xbb_Lov("I,i,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_222 = 0.500000 * scalar_159;

            // energy += -0.500000 <j,k||b,i>_abab u2_abab(b,a,j,k) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_222;

            // energy += -0.500000 <k,j||b,i>_abab u2_abab(b,a,k,j) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_222;
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // scalar_tmp_223 += 1.000000 sharedOps["aa_vo_325"]("a,i") m1_xaa_Lov("I,i,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_223 = scalar_160;

            // rm0 += -1.000000 dp_bb(j,b) u2_abab(a,b,i,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_223;

            // energy += -1.000000 dp_bb(j,b) u0 u2_abab(a,b,i,j) m1_aa(i,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_223 * u0;
        }

        if (include_u0_ && include_u1_) {

            // scalar_tmp_224 += 1.000000 m1_xbb_Lov("I,i,a") sharedOps["bb_vo_365"]("a,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_224 = scalar_161;

            // rm0 += 1.000000 dp_bb(j,i) u1_bb(a,j) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_224;

            // energy += 1.000000 dp_bb(j,i) u0 u1_bb(a,j) m1_bb(i,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_224 * u0;
        }

        {

            // scalar_tmp_225 += 1.000000 l1_xbb_Lov("I,i,a") dp_bb_vo("a,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_225 = scalar_162;
        }

        if (include_u0_) {

            // rm0 += -1.000000 dp_bb(a,i) l1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_tmp_225;

            // energy += -1.000000 dp_bb(a,i) l1_bb(i,a) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_tmp_225 * u0;
        }

        if (include_u1_ && include_u2_) {

            // scalar_tmp_226 += 1.000000 sharedOps["aa_vo_255"]("a,i") m1_xaa_Lov("I,i,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_226 = 0.500000 * scalar_163;

            // energy += -0.500000 <k,j||i,b>_abab u2_abab(a,b,k,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_226;

            // energy += -0.500000 <j,k||i,b>_abab u2_abab(a,b,j,k) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_tmp_226;
        }

        if (include_u0_ && include_u1_) {

            // scalar_tmp_227 += 1.000000 sharedOps["aa_vo_364"]("a,i") m1_xaa_Lov("I,i,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_227 = scalar_164;

            // rm0 += 1.000000 dp_aa(j,i) u1_aa(a,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_227;

            // energy += 1.000000 dp_aa(j,i) u0 u1_aa(a,j) m1_aa(i,a)
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_227 * u0;
        }

        if (include_u0_) {

            // scalar_tmp_228 += 1.000000 sharedOps["bb_vo_323"]("a,i") l1_xbb_Lov("I,i,a")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_228 = scalar_165;

            // rm0 += 1.000000 dp_bb(j,b) l1_bb(i,a) t2_bbbb(b,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 += scalar_tmp_228;

            // energy += 1.000000 dp_bb(j,b) l1_bb(i,a) t2_bbbb(b,a,i,j) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_tmp_228 * u0;
        }

        if (include_u1_) {

            // scalar_tmp_229 += 1.000000 m1_xaa_Lov("I,i,a") sharedOps["aa_vo_436"]("a,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_229 = 0.500000 * scalar_166;

            // energy += -0.500000 <j,k||b,c>_abab t2_abab(a,c,j,k) u1_aa(b,i) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_229;

            // energy += -0.500000 <k,j||b,c>_abab t2_abab(a,c,k,j) u1_aa(b,i) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_229;

            // scalar_tmp_230 += 1.000000 m1_xaa_Lov("I,i,a") sharedOps["aa_vo_500"]("a,i")
            // flops: o0v0: 1 | mem: o0v0: 1,
            scalar_tmp_230 = 0.500000 * scalar_167;

            // energy += -0.500000 <j,k||b,c>_abab t2_abab(b,c,i,k) u1_aa(a,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_230;

            // energy += -0.500000 <j,k||c,b>_abab t2_abab(c,b,i,k) u1_aa(a,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_tmp_230;
        }

        if (include_u2_) {

            // crm2_bbbb_vvoo += 1.000000 l2_bbbb(m,n,e,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            crm2_bbbb_vvoo("e,f,m,n") += l2_bbbb_oovv("m,n,e,f");

            // crm2_abab_vvoo += 1.000000 l2_abab(m,n,e,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            crm2_abab_vvoo("e,f,m,n") += l2_abab_oovv("m,n,e,f");

            // crm2_aaaa_vvoo += 1.000000 l2_aaaa(m,n,e,f)
            // flops: o2v2: 1 | mem: o2v2: 1,
            crm2_aaaa_vvoo("e,f,m,n") += l2_aaaa_oovv("m,n,e,f");
        }

        {

            // rl2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <j,n||b,e>_abab t2_aaaa(b,a,i,j) m2_abab(i,m,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = sharedOps["abab_vvoo_0"]("a,e,i,n") * m2_abab_oovv("i,m,a,f");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <j,n||b,e>_abab t2_abab(b,a,j,i) m2_bbbb(i,m,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["bbbb_vvoo_1"]("e,a,n,i") * m2_bbbb_oovv("i,m,f,a");

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) <a,n||i,e>_abab m2_abab(i,m,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += V_blks_["abba_vovo"]("a,n,e,i") * m2_abab_oovv("i,m,a,f");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <j,n||e,b>_bbbb t2_bbbb(b,a,i,j) m2_bbbb(i,m,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["bbbb_vvoo_2"]("a,e,i,n") * m2_bbbb_oovv("i,m,f,a");
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) f_bb(n,e) m1_bb(m,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= F_blks_["bb_ov"]("n,e") * m1_bb_ov("m,f");
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <n,a||e,i>_bbbb m2_bbbb(m,i,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= V_blks_["bbbb_vovo"]("a,n,e,i") * m2_bbbb_oovv("m,i,a,f");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <j,n||e,b>_bbbb t2_abab(a,b,i,j) m2_abab(i,m,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["abab_vvoo_3"]("a,e,i,n") * m2_abab_oovv("i,m,a,f");

            // rm2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rm2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
            rm2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
            rm2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("f,e,n,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += 0.500000 P(m,n) <j,n||a,b>_abab t2_abab(a,b,j,i) m2_bbbb(i,m,e,f)
            // +                0.500000 P(m,n) <j,n||b,a>_abab t2_abab(b,a,j,i) m2_bbbb(i,m,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = sharedOps["bb_oo_208"]("n,i") * m2_bbbb_oovv("i,m,e,f");

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) f_bb(n,i) m2_bbbb(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= F_blks_["bb_oo"]("n,i") * m2_bbbb_oovv("m,i,e,f");

            // tempPerm_bbbb_vvoo += -0.500000 P(m,n) <j,n||a,b>_bbbb t2_bbbb(a,b,i,j) m2_bbbb(i,m,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= 0.500000 * sharedOps["bb_oo_212"]("n,i") * m2_bbbb_oovv("i,m,e,f");
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) <n,a||e,f>_bbbb m1_bb(m,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += V_blks_["bbbb_vovv"]("a,n,e,f") * m1_bb_ov("m,a");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_bbbb_vvoo += -2.000000 P(m,n) dp_bb(n,a) u1_bb(a,i) m2_bbbb(i,m,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= 2.000000 * sharedOps["bb_oo_313"]("n,i") * m2_bbbb_oovv("i,m,e,f");
        }

        if (include_u2_) {

            // rm2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rm2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(e,f) <m,n||e,i>_bbbb m1_bb(i,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * V_blks_["bbbb_oovo"]("m,n,e,i") * m1_bb_ov("i,f");
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(e,f) f_bb(a,e) m2_bbbb(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += F_blks_["bb_vv"]("a,e") * m2_bbbb_oovv("m,n,a,f");

            // tempPerm_bbbb_vvoo += 0.500000 P(e,f) <j,i||b,e>_abab t2_abab(b,a,j,i) m2_bbbb(m,n,f,a)
            // +                0.500000 P(e,f) <i,j||b,e>_abab t2_abab(b,a,i,j) m2_bbbb(m,n,f,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["bb_vv_149"]("e,a") * m2_bbbb_oovv("m,n,f,a");

            // tempPerm_bbbb_vvoo += -0.500000 P(e,f) <j,i||e,b>_bbbb t2_bbbb(b,a,j,i) m2_bbbb(m,n,f,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= 0.500000 * sharedOps["bb_vv_152"]("e,a") * m2_bbbb_oovv("m,n,f,a");
        }

        if (include_u0_ && include_u2_) {

            // rm2_bbbb_vvoo += 1.000000 <m,n||e,f>_bbbb m0
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rm2_bbbb_vvoo("e,f,m,n") += V_blks_["bbbb_oovv"]("m,n,e,f") * m0;
        }

        if (include_u2_) {

            // rm2_bbbb_vvoo += -1.000000 dp_aa(i,i) l2_bbbb(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_bbbb_vvoo("e,f,m,n") -= scalar_4 * l2_bbbb_oovv("m,n,e,f");

            // rm2_bbbb_vvoo += 0.500000 <b,a||e,f>_bbbb m2_bbbb(m,n,b,a)
            // flops: o2v4L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_bbbb_vvoo("e,f,m,n") += 0.500000 * V_blks_["bbbb_vvvv"]("b,a,e,f") * m2_bbbb_oovv("m,n,b,a");

            // rm2_bbbb_vvoo += -1.000000 dp_bb(i,i) l2_bbbb(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_bbbb_vvoo("e,f,m,n") -= scalar_5 * l2_bbbb_oovv("m,n,e,f");
        }

        if (include_u1_ && include_u2_) {

            // rm2_bbbb_vvoo += -1.000000 dp_aa(i,a) u1_aa(a,i) m2_bbbb(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_bbbb_vvoo("e,f,m,n") -= scalar_6 * m2_bbbb_oovv("m,n,e,f");
        }

        if (include_u2_) {

            // rm2_bbbb_vvoo += 1.000000 m2_bbbb(m,n,e,f) w0
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_bbbb_vvoo("e,f,m,n") += m2_bbbb_oovv("m,n,e,f") * w0;

            // rm2_bbbb_vvoo += 0.250000 <j,i||e,f>_bbbb t2_bbbb(b,a,j,i) m2_bbbb(m,n,b,a)
            // flops: o4v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o4v0L1: 1,
            rm2_bbbb_vvoo("e,f,m,n") += 0.250000 * t2_bbbb_vvoo("b,a,j,i") * m2_bbbb_oovv("m,n,b,a") * V_blks_["bbbb_oovv"]("j,i,e,f");
        }

        if (include_u1_ && include_u2_) {

            // rm2_bbbb_vvoo += -1.000000 dp_bb(i,a) u1_bb(a,i) m2_bbbb(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_bbbb_vvoo("e,f,m,n") -= scalar_7 * m2_bbbb_oovv("m,n,e,f");
        }

        if (include_u2_) {

            // rm2_bbbb_vvoo += 0.250000 <m,n||a,b>_bbbb t2_bbbb(a,b,i,j) m2_bbbb(i,j,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_bbbb_vvoo("e,f,m,n") += 0.250000 * sharedOps["bbbb_oooo_107"]("i,j,m,n") * m2_bbbb_oovv("i,j,e,f");

            // rm2_bbbb_vvoo += 0.500000 <m,n||i,j>_bbbb m2_bbbb(i,j,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_bbbb_vvoo("e,f,m,n") += 0.500000 * V_blks_["bbbb_oooo"]("m,n,i,j") * m2_bbbb_oovv("i,j,e,f");

            // rm2_abab_vvoo += 0.500000 <j,i||e,b>_aaaa t2_aaaa(b,a,j,i) m2_abab(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += 0.500000 * sharedOps["aa_vv_151"]("a,e") * m2_abab_oovv("m,n,a,f");

            // rm2_abab_vvoo += 0.250000 <j,i||e,f>_abab t2_abab(b,a,j,i) m2_abab(m,n,b,a)
            // +                0.250000 <i,j||e,f>_abab t2_abab(b,a,i,j) m2_abab(m,n,b,a)
            // +                0.250000 <j,i||e,f>_abab t2_abab(a,b,j,i) m2_abab(m,n,a,b)
            // +                0.250000 <i,j||e,f>_abab t2_abab(a,b,i,j) m2_abab(m,n,a,b)
            // flops: o4v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o4v0L1: 1,
            rm2_abab_vvoo("e,f,m,n") += t2_abab_vvoo("b,a,j,i") * m2_abab_oovv("m,n,b,a") * V_blks_["abab_oovv"]("j,i,e,f");
        }

        if (include_u1_ && include_u2_) {

            // rm2_abab_vvoo += 2.000000 dp_bb(n,a) u1_bb(a,i) m2_abab(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += 2.000000 * sharedOps["bb_oo_313"]("n,i") * m2_abab_oovv("m,i,e,f");

            // rm2_abab_vvoo += -1.000000 <m,n||e,i>_abab m1_bb(i,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= V_blks_["abab_oovo"]("m,n,e,i") * m1_bb_ov("i,f");
        }

        if (include_u2_) {

            // rm2_abab_vvoo += 1.000000 <m,j||e,b>_abab t2_bbbb(b,a,i,j) m2_bbbb(i,n,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += sharedOps["abab_vvoo_7"]("e,a,m,i") * m2_bbbb_oovv("i,n,f,a");

            // rm2_abab_vvoo += 1.000000 m2_abab(m,n,e,f) w0
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_abab_vvoo("e,f,m,n") += m2_abab_oovv("m,n,e,f") * w0;

            // rm2_abab_vvoo += 0.250000 <m,n||a,b>_abab t2_abab(a,b,i,j) m2_abab(i,j,e,f)
            // +                0.250000 <m,n||a,b>_abab t2_abab(a,b,j,i) m2_abab(j,i,e,f)
            // +                0.250000 <m,n||b,a>_abab t2_abab(b,a,i,j) m2_abab(i,j,e,f)
            // +                0.250000 <m,n||b,a>_abab t2_abab(b,a,j,i) m2_abab(j,i,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += sharedOps["aabb_oooo_103"]("m,i,n,j") * m2_abab_oovv("i,j,e,f");
        }

        if (include_u1_ && include_u2_) {

            // rm2_abab_vvoo += 1.000000 <a,n||e,f>_abab m1_aa(m,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += V_blks_["abab_vovv"]("a,n,e,f") * m1_aa_ov("m,a");
        }

        if (include_u2_) {

            // rm2_abab_vvoo += -1.000000 <m,a||e,i>_abab m2_bbbb(n,i,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += V_blks_["baab_vovo"]("a,m,e,i") * m2_bbbb_oovv("n,i,a,f");

            // rm2_abab_vvoo += -1.000000 <a,n||e,i>_abab m2_abab(m,i,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= V_blks_["abab_vovo"]("a,n,e,i") * m2_abab_oovv("m,i,a,f");

            // rm2_abab_vvoo += -1.000000 dp_aa(i,i) l2_abab(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_abab_vvoo("e,f,m,n") -= scalar_4 * l2_abab_oovv("m,n,e,f");

            // rm2_abab_vvoo += -1.000000 <a,n||i,f>_abab m2_aaaa(m,i,a,e)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += V_blks_["abba_vovo"]("a,n,f,i") * m2_aaaa_oovv("m,i,a,e");

            // rm2_abab_vvoo += 1.000000 <n,a||f,i>_bbbb m2_abab(m,i,e,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= V_blks_["bbbb_vovo"]("a,n,f,i") * m2_abab_oovv("m,i,e,a");

            // rm2_abab_vvoo += -0.500000 <m,j||a,b>_abab t2_abab(a,b,i,j) m2_abab(i,n,e,f)
            // +                -0.500000 <m,j||b,a>_abab t2_abab(b,a,i,j) m2_abab(i,n,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= sharedOps["aa_oo_207"]("m,i") * m2_abab_oovv("i,n,e,f");

            // rm2_abab_vvoo += -0.500000 <j,n||a,b>_abab t2_abab(a,b,j,i) m2_abab(m,i,e,f)
            // +                -0.500000 <j,n||b,a>_abab t2_abab(b,a,j,i) m2_abab(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= sharedOps["bb_oo_208"]("n,i") * m2_abab_oovv("m,i,e,f");
        }

        if (include_u1_ && include_u2_) {

            // rm2_abab_vvoo += -1.000000 <m,n||i,f>_abab m1_aa(i,e)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += V_blks_["abba_oovo"]("m,n,f,i") * m1_aa_ov("i,e");
        }

        if (include_u0_ && include_u2_) {

            // rm2_abab_vvoo += 1.000000 <m,n||e,f>_abab m0
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += V_blks_["abab_oovv"]("m,n,e,f") * m0;
        }

        if (include_u2_) {

            // rm2_abab_vvoo += -0.500000 <j,i||b,f>_abab t2_abab(b,a,j,i) m2_abab(m,n,e,a)
            // +                -0.500000 <i,j||b,f>_abab t2_abab(b,a,i,j) m2_abab(m,n,e,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= sharedOps["bb_vv_149"]("f,a") * m2_abab_oovv("m,n,e,a");

            // rm2_abab_vvoo += -1.000000 dp_bb(i,i) l2_abab(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_abab_vvoo("e,f,m,n") -= scalar_5 * l2_abab_oovv("m,n,e,f");

            // rm2_abab_vvoo += 1.000000 <m,a||e,i>_aaaa m2_abab(i,n,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= V_blks_["aaaa_vovo"]("a,m,e,i") * m2_abab_oovv("i,n,a,f");

            // rm2_abab_vvoo += 0.500000 <b,a||e,f>_abab m2_abab(m,n,b,a)
            // +                0.500000 <a,b||e,f>_abab m2_abab(m,n,a,b)
            // flops: o2v4L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += V_blks_["abab_vvvv"]("b,a,e,f") * m2_abab_oovv("m,n,b,a");

            // rm2_abab_vvoo += 1.000000 <j,n||b,f>_abab t2_abab(b,a,j,i) m2_abab(m,i,e,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += sharedOps["bbbb_vvoo_1"]("f,a,n,i") * m2_abab_oovv("m,i,e,a");

            // rm2_abab_vvoo += 1.000000 <m,j||b,f>_abab t2_abab(b,a,i,j) m2_abab(i,n,e,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += sharedOps["bbaa_vvoo_8"]("f,a,m,i") * m2_abab_oovv("i,n,e,a");

            // rm2_abab_vvoo += 0.500000 <j,i||f,b>_bbbb t2_bbbb(b,a,j,i) m2_abab(m,n,e,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += 0.500000 * sharedOps["bb_vv_152"]("f,a") * m2_abab_oovv("m,n,e,a");

            // rm2_abab_vvoo += -1.000000 <m,a||i,f>_abab m2_abab(i,n,e,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= V_blks_["baba_vovo"]("a,m,f,i") * m2_abab_oovv("i,n,e,a");

            // rm2_abab_vvoo += 1.000000 <m,j||e,b>_abab t2_abab(a,b,i,j) m2_abab(i,n,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += sharedOps["aaaa_vvoo_6"]("a,e,i,m") * m2_abab_oovv("i,n,a,f");
        }

        if (include_u1_ && include_u2_) {

            // rm2_abab_vvoo += -1.000000 dp_aa(i,a) u1_aa(a,i) m2_abab(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_abab_vvoo("e,f,m,n") -= scalar_6 * m2_abab_oovv("m,n,e,f");
        }

        if (include_u2_) {

            // rm2_abab_vvoo += 1.000000 <j,n||f,b>_bbbb t2_bbbb(b,a,i,j) m2_abab(m,i,e,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += sharedOps["bbbb_vvoo_2"]("a,f,i,n") * m2_abab_oovv("m,i,e,a");

            // rm2_abab_vvoo += -1.000000 f_bb(n,i) m2_abab(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= F_blks_["bb_oo"]("n,i") * m2_abab_oovv("m,i,e,f");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // rm2_abab_vvoo += -1.000000 dp_bb(n,f) u0 m1_aa(m,e)
            // flops: o2v2L1: 2, o0v0: 1 | mem: o2v2L1: 2, o0v0: 1,
            rm2_abab_vvoo("e,f,m,n") -= dp_bb_ov("n,f") * u0 * m1_aa_ov("m,e");
        }

        if (include_u2_) {

            // rm2_abab_vvoo += 1.000000 f_aa(a,e) m2_abab(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += F_blks_["aa_vv"]("a,e") * m2_abab_oovv("m,n,a,f");
        }

        if (include_u1_ && include_u2_) {

            // rm2_abab_vvoo += -1.000000 dp_bb(i,a) u1_bb(a,i) m2_abab(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_abab_vvoo("e,f,m,n") -= scalar_7 * m2_abab_oovv("m,n,e,f");
        }

        if (include_u2_) {

            // rm2_abab_vvoo += 1.000000 <j,n||e,b>_abab t2_abab(a,b,j,i) m2_abab(m,i,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += sharedOps["aabb_vvoo_9"]("e,a,n,i") * m2_abab_oovv("m,i,a,f");
        }

        if (include_u1_ && include_u2_) {

            // rm2_abab_vvoo += 1.000000 <m,a||e,f>_abab m1_bb(n,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= V_blks_["baab_vovv"]("a,m,e,f") * m1_bb_ov("n,a");

            // rm2_abab_vvoo += 1.000000 f_bb(n,f) m1_aa(m,e)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += F_blks_["bb_ov"]("n,f") * m1_aa_ov("m,e");
        }

        if (include_u2_) {

            // rm2_abab_vvoo += 0.500000 <m,n||i,j>_abab m2_abab(i,j,e,f)
            // +                0.500000 <m,n||j,i>_abab m2_abab(j,i,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += V_blks_["abab_oooo"]("m,n,i,j") * m2_abab_oovv("i,j,e,f");
        }

        if (include_u1_ && include_u2_) {

            // rm2_abab_vvoo += 1.000000 f_aa(m,e) m1_bb(n,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += F_blks_["aa_ov"]("m,e") * m1_bb_ov("n,f");
        }

        if (include_u2_) {

            // rm2_abab_vvoo += 1.000000 <j,m||e,b>_aaaa t2_abab(b,a,j,i) m2_bbbb(i,n,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += sharedOps["abab_vvoo_4"]("e,a,m,i") * m2_bbbb_oovv("i,n,f,a");
        }

        if (include_u1_ && include_u2_) {

            // rm2_abab_vvoo += 2.000000 dp_aa(m,a) u1_aa(a,i) m2_abab(i,n,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += 2.000000 * sharedOps["aa_oo_312"]("i,m") * m2_abab_oovv("i,n,e,f");
        }

        if (include_u2_) {

            // rm2_abab_vvoo += -1.000000 dp_bb(n,f) l1_aa(m,e)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= dp_bb_ov("n,f") * l1_aa_ov("m,e");

            // rm2_abab_vvoo += 1.000000 <j,m||e,b>_aaaa t2_aaaa(b,a,i,j) m2_abab(i,n,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += sharedOps["aaaa_vvoo_5"]("e,a,m,i") * m2_abab_oovv("i,n,a,f");

            // rm2_abab_vvoo += 0.500000 <j,m||a,b>_aaaa t2_aaaa(a,b,i,j) m2_abab(i,n,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += 0.500000 * sharedOps["aa_oo_211"]("i,m") * m2_abab_oovv("i,n,e,f");

            // rm2_abab_vvoo += 1.000000 f_bb(a,f) m2_abab(m,n,e,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += F_blks_["bb_vv"]("a,f") * m2_abab_oovv("m,n,e,a");

            // rm2_abab_vvoo += 1.000000 <j,n||f,b>_bbbb t2_abab(a,b,i,j) m2_aaaa(i,m,e,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += sharedOps["abab_vvoo_3"]("a,f,i,n") * m2_aaaa_oovv("i,m,e,a");

            // rm2_abab_vvoo += -1.000000 f_aa(m,i) m2_abab(i,n,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= F_blks_["aa_oo"]("m,i") * m2_abab_oovv("i,n,e,f");

            // rm2_abab_vvoo += -0.500000 <j,i||e,b>_abab t2_abab(a,b,j,i) m2_abab(m,n,a,f)
            // +                -0.500000 <i,j||e,b>_abab t2_abab(a,b,i,j) m2_abab(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= sharedOps["aa_vv_150"]("a,e") * m2_abab_oovv("m,n,a,f");

            // rm2_abab_vvoo += 1.000000 <j,n||b,f>_abab t2_aaaa(b,a,i,j) m2_aaaa(i,m,e,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += sharedOps["abab_vvoo_0"]("a,f,i,n") * m2_aaaa_oovv("i,m,e,a");

            // rm2_abab_vvoo += 0.500000 <j,n||a,b>_bbbb t2_bbbb(a,b,i,j) m2_abab(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") += 0.500000 * sharedOps["bb_oo_212"]("n,i") * m2_abab_oovv("m,i,e,f");

            // rm2_abab_vvoo += -1.000000 dp_aa(m,e) l1_bb(n,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rm2_abab_vvoo("e,f,m,n") -= dp_aa_ov("m,e") * l1_bb_ov("n,f");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // rm2_abab_vvoo += -1.000000 dp_aa(m,e) u0 m1_bb(n,f)
            // flops: o2v2L1: 2, o0v0: 1 | mem: o2v2L1: 2, o0v0: 1,
            rm2_abab_vvoo("e,f,m,n") -= dp_aa_ov("m,e") * u0 * m1_bb_ov("n,f");
        }

        {

            // rl2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <n,j||e,b>_abab t2_bbbb(b,a,i,j) m2_abab(m,i,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = sharedOps["abab_vvoo_7"]("e,a,n,i") * m2_abab_oovv("m,i,f,a");

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) <n,a||e,i>_abab m2_abab(m,i,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += V_blks_["baab_vovo"]("a,n,e,i") * m2_abab_oovv("m,i,f,a");

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <j,n||e,b>_aaaa t2_aaaa(b,a,i,j) m2_aaaa(i,m,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["aaaa_vvoo_5"]("e,a,n,i") * m2_aaaa_oovv("i,m,f,a");

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <n,j||e,b>_abab t2_abab(a,b,i,j) m2_aaaa(i,m,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["aaaa_vvoo_6"]("a,e,i,n") * m2_aaaa_oovv("i,m,f,a");

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <j,n||e,b>_aaaa t2_abab(b,a,j,i) m2_abab(m,i,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["abab_vvoo_4"]("e,a,n,i") * m2_abab_oovv("m,i,f,a");
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) f_aa(n,e) m1_aa(m,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= F_blks_["aa_ov"]("n,e") * m1_aa_ov("m,f");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <n,a||e,i>_aaaa m2_aaaa(m,i,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= V_blks_["aaaa_vovo"]("a,n,e,i") * m2_aaaa_oovv("m,i,a,f");

            // rm2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rm2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
            rm2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
            rm2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("f,e,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += -0.500000 P(m,n) <j,n||a,b>_aaaa t2_aaaa(a,b,i,j) m2_aaaa(i,m,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -0.500000 * sharedOps["aa_oo_211"]("i,n") * m2_aaaa_oovv("i,m,e,f");

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) f_aa(n,i) m2_aaaa(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= F_blks_["aa_oo"]("n,i") * m2_aaaa_oovv("m,i,e,f");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += -2.000000 P(m,n) dp_aa(n,a) u1_aa(a,i) m2_aaaa(i,m,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= 2.000000 * sharedOps["aa_oo_312"]("i,n") * m2_aaaa_oovv("i,m,e,f");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += 0.500000 P(m,n) <n,j||a,b>_abab t2_abab(a,b,i,j) m2_aaaa(i,m,e,f)
            // +                0.500000 P(m,n) <n,j||b,a>_abab t2_abab(b,a,i,j) m2_aaaa(i,m,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["aa_oo_207"]("n,i") * m2_aaaa_oovv("i,m,e,f");
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) <n,a||e,f>_aaaa m1_aa(m,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += V_blks_["aaaa_vovv"]("a,n,e,f") * m1_aa_ov("m,a");
        }

        if (include_u2_) {

            // rm2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rm2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(e,f) <m,n||e,i>_aaaa m1_aa(i,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * V_blks_["aaaa_oovo"]("m,n,e,i") * m1_aa_ov("i,f");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += -0.500000 P(e,f) <j,i||e,b>_aaaa t2_aaaa(b,a,j,i) m2_aaaa(m,n,f,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= 0.500000 * sharedOps["aa_vv_151"]("a,e") * m2_aaaa_oovv("m,n,f,a");

            // tempPerm_aaaa_vvoo += 1.000000 P(e,f) f_aa(a,e) m2_aaaa(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += F_blks_["aa_vv"]("a,e") * m2_aaaa_oovv("m,n,a,f");

            // tempPerm_aaaa_vvoo += 0.500000 P(e,f) <j,i||e,b>_abab t2_abab(a,b,j,i) m2_aaaa(m,n,f,a)
            // +                0.500000 P(e,f) <i,j||e,b>_abab t2_abab(a,b,i,j) m2_aaaa(m,n,f,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["aa_vv_150"]("a,e") * m2_aaaa_oovv("m,n,f,a");

            // rm2_aaaa_vvoo += 0.250000 <m,n||a,b>_aaaa t2_aaaa(a,b,i,j) m2_aaaa(i,j,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_aaaa_vvoo("e,f,m,n") += 0.250000 * sharedOps["aaaa_oooo_108"]("i,j,m,n") * m2_aaaa_oovv("i,j,e,f");

            // rm2_aaaa_vvoo += 0.250000 <j,i||e,f>_aaaa t2_aaaa(b,a,j,i) m2_aaaa(m,n,b,a)
            // flops: o4v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o4v0L1: 1,
            rm2_aaaa_vvoo("e,f,m,n") += 0.250000 * t2_aaaa_vvoo("b,a,j,i") * m2_aaaa_oovv("m,n,b,a") * V_blks_["aaaa_oovv"]("j,i,e,f");
        }

        if (include_u0_ && include_u2_) {

            // rm2_aaaa_vvoo += 1.000000 <m,n||e,f>_aaaa m0
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rm2_aaaa_vvoo("e,f,m,n") += V_blks_["aaaa_oovv"]("m,n,e,f") * m0;
        }

        if (include_u2_) {

            // rm2_aaaa_vvoo += 1.000000 m2_aaaa(m,n,e,f) w0
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_aaaa_vvoo("e,f,m,n") += m2_aaaa_oovv("m,n,e,f") * w0;

            // rm2_aaaa_vvoo += -1.000000 dp_aa(i,i) l2_aaaa(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_aaaa_vvoo("e,f,m,n") -= scalar_4 * l2_aaaa_oovv("m,n,e,f");
        }

        if (include_u1_ && include_u2_) {

            // rm2_aaaa_vvoo += -1.000000 dp_aa(i,a) u1_aa(a,i) m2_aaaa(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_aaaa_vvoo("e,f,m,n") -= scalar_6 * m2_aaaa_oovv("m,n,e,f");
        }

        if (include_u2_) {

            // rm2_aaaa_vvoo += 0.500000 <m,n||i,j>_aaaa m2_aaaa(i,j,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_aaaa_vvoo("e,f,m,n") += 0.500000 * V_blks_["aaaa_oooo"]("m,n,i,j") * m2_aaaa_oovv("i,j,e,f");

            // rm2_aaaa_vvoo += -1.000000 dp_bb(i,i) l2_aaaa(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_aaaa_vvoo("e,f,m,n") -= scalar_5 * l2_aaaa_oovv("m,n,e,f");

            // rm2_aaaa_vvoo += 0.500000 <b,a||e,f>_aaaa m2_aaaa(m,n,b,a)
            // flops: o2v4L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rm2_aaaa_vvoo("e,f,m,n") += 0.500000 * V_blks_["aaaa_vvvv"]("b,a,e,f") * m2_aaaa_oovv("m,n,b,a");
        }

        if (include_u1_ && include_u2_) {

            // rm2_aaaa_vvoo += -1.000000 dp_bb(i,a) u1_bb(a,i) m2_aaaa(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            rm2_aaaa_vvoo("e,f,m,n") -= scalar_7 * m2_aaaa_oovv("m,n,e,f");
        }

        if (include_u1_) {

            // crm1_bb_vo += 1.000000 l1_bb(m,e)
            // flops: o1v1: 1 | mem: o1v1: 1,
            crm1_bb_vo("e,m") += l1_bb_ov("m,e");

            // crm1_aa_vo += 1.000000 l1_aa(m,e)
            // flops: o1v1: 1 | mem: o1v1: 1,
            crm1_aa_vo("e,m") += l1_aa_ov("m,e");

            // rm1_bb_vo += 1.000000 <a,m||i,e>_abab m1_aa(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= V_blks_["abba_vovo"]("a,m,e,i") * m1_aa_ov("i,a");
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += 1.000000 <k,m||b,j>_abab t2_aaaa(b,a,i,k) m2_abab(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += sharedOps["aabb_vooo_120"]("a,i,m,j") * m2_abab_oovv("i,j,a,e");
        }

        if (include_u1_) {

            // rm1_bb_vo += -1.000000 dp_bb(m,e)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rm1_bb_vo("e,m") -= dp_bb_ov("m,e");

            // rm1_bb_vo += 0.500000 dp_bb(m,b) l2_bbbb(i,j,e,a) t2_bbbb(b,a,i,j)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += 0.500000 * sharedOps["bbbb_vooo_216"]("a,i,j,m") * l2_bbbb_oovv("i,j,e,a");

            // rm1_bb_vo += 0.500000 <j,i||e,b>_bbbb t2_bbbb(b,a,j,i) m1_bb(m,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += 0.500000 * sharedOps["bb_vv_152"]("e,a") * m1_bb_ov("m,a");
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += 1.000000 <k,m||b,j>_abab t2_abab(b,a,k,i) m2_bbbb(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += sharedOps["bbbb_vooo_116"]("a,i,m,j") * m2_bbbb_oovv("i,j,e,a");

            // rm1_bb_vo += -1.000000 <j,b||e,c>_bbbb t2_bbbb(c,a,i,j) m2_bbbb(i,m,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += sharedOps["bbbb_vvvo_20"]("a,b,e,i") * m2_bbbb_oovv("i,m,b,a");

            // rm1_bb_vo += -0.500000 f_bb(m,b) t2_bbbb(b,a,i,j) m2_bbbb(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= 0.500000 * sharedOps["bbbb_vooo_243"]("a,m,i,j") * m2_bbbb_oovv("i,j,e,a");
        }

        if (include_u1_) {

            // rm1_bb_vo += -1.000000 dp_bb(j,b) l2_bbbb(i,m,e,a) t2_bbbb(b,a,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= sharedOps["bb_vo_323"]("a,i") * l2_bbbb_oovv("i,m,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += 1.000000 dp_aa(j,b) u2_aaaa(b,a,i,j) m2_abab(i,m,a,e)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += sharedOps["aa_vo_324"]("a,i") * m2_abab_oovv("i,m,a,e");

            // rm1_bb_vo += 0.250000 <k,j||e,i>_bbbb t2_bbbb(b,a,k,j) m2_bbbb(m,i,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += 0.250000 * sharedOps["bbbb_vvvo_44"]("b,a,e,i") * m2_bbbb_oovv("m,i,b,a");
        }

        if (include_u1_) {

            // rm1_bb_vo += 1.000000 dp_aa(j,b) l2_bbbb(i,m,e,a) t2_abab(b,a,j,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += sharedOps["bb_vo_321"]("a,i") * l2_bbbb_oovv("i,m,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += -1.000000 dp_aa(a,b) u1_aa(b,i) m2_abab(i,m,a,e)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= sharedOps["aa_vo_360"]("a,i") * m2_abab_oovv("i,m,a,e");

            // rm1_bb_vo += -1.000000 dp_bb(a,b) u1_bb(b,i) m2_bbbb(i,m,a,e)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= sharedOps["bb_vo_361"]("a,i") * m2_bbbb_oovv("i,m,a,e");
        }

        if (include_u1_) {

            // rm1_bb_vo += -1.000000 dp_aa(i,i) l1_bb(m,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rm1_bb_vo("e,m") -= scalar_4 * l1_bb_ov("m,e");

            // rm1_bb_vo += -0.500000 <j,i||b,e>_abab t2_abab(b,a,j,i) m1_bb(m,a)
            // +            -0.500000 <i,j||b,e>_abab t2_abab(b,a,i,j) m1_bb(m,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= sharedOps["bb_vv_149"]("e,a") * m1_bb_ov("m,a");
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += -1.000000 <j,b||e,c>_bbbb t2_abab(a,c,i,j) m2_abab(i,m,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += sharedOps["abba_vvvo_21"]("a,b,e,i") * m2_abab_oovv("i,m,a,b");

            // rm1_bb_vo += 0.250000 <k,j||i,e>_abab t2_abab(b,a,k,j) m2_abab(i,m,b,a)
            // +            0.250000 <j,k||i,e>_abab t2_abab(b,a,j,k) m2_abab(i,m,b,a)
            // +            0.250000 <k,j||i,e>_abab t2_abab(a,b,k,j) m2_abab(i,m,a,b)
            // +            0.250000 <j,k||i,e>_abab t2_abab(a,b,j,k) m2_abab(i,m,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= sharedOps["abba_vvvo_46"]("b,e,a,i") * m2_abab_oovv("i,m,b,a");
        }

        if (include_u1_) {

            // rm1_bb_vo += -0.500000 <j,m||a,b>_abab t2_abab(a,b,j,i) m1_bb(i,e)
            // +            -0.500000 <j,m||b,a>_abab t2_abab(b,a,j,i) m1_bb(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= sharedOps["bb_oo_208"]("m,i") * m1_bb_ov("i,e");

            // rm1_bb_vo += 2.000000 dp_bb(m,a) u1_bb(a,i) m1_bb(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += 2.000000 * sharedOps["bb_oo_313"]("m,i") * m1_bb_ov("i,e");
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += 1.000000 dp_aa(j,i) u1_aa(a,j) m2_abab(i,m,a,e)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += sharedOps["aa_vo_364"]("a,i") * m2_abab_oovv("i,m,a,e");

            // rm1_bb_vo += -0.500000 f_bb(m,b) t2_abab(a,b,i,j) m2_abab(i,j,a,e)
            // +            -0.500000 f_bb(m,b) t2_abab(a,b,j,i) m2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= sharedOps["aabb_vooo_262"]("a,i,j,m") * m2_abab_oovv("i,j,a,e");
        }

        if (include_u1_) {

            // rm1_bb_vo += 1.000000 <m,a||e,i>_bbbb m1_bb(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= V_blks_["bbbb_vovo"]("a,m,e,i") * m1_bb_ov("i,a");
        }

        if (include_u0_ && include_u1_) {

            // rm1_bb_vo += 1.000000 f_bb(m,e) m0
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += F_blks_["bb_ov"]("m,e") * m0;
        }

        if (include_u1_) {

            // rm1_bb_vo += 1.000000 dp_aa(j,b) l2_abab(i,m,a,e) t2_aaaa(b,a,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += sharedOps["aa_vo_320"]("a,i") * l2_abab_oovv("i,m,a,e");
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += 1.000000 dp_bb(j,i) u1_bb(a,j) m2_bbbb(m,i,e,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += sharedOps["bb_vo_365"]("a,i") * m2_bbbb_oovv("m,i,e,a");

            // rm1_bb_vo += -1.000000 <k,m||b,j>_bbbb t2_bbbb(b,a,i,k) m2_bbbb(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= sharedOps["bbbb_vooo_112"]("a,i,m,j") * m2_bbbb_oovv("i,j,e,a");

            // rm1_bb_vo += 0.500000 <m,a||i,j>_bbbb m2_bbbb(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= 0.500000 * V_blks_["bbbb_vooo"]("a,m,i,j") * m2_bbbb_oovv("i,j,a,e");

            // rm1_bb_vo += 0.500000 <b,a||e,i>_bbbb m2_bbbb(m,i,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += 0.500000 * V_blks_["bbbb_vvvo"]("b,a,e,i") * m2_bbbb_oovv("m,i,b,a");

            // rm1_bb_vo += 1.000000 <k,m||j,b>_abab t2_abab(a,b,k,i) m2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= sharedOps["aabb_vooo_113"]("a,j,m,i") * m2_abab_oovv("j,i,a,e");

            // rm1_bb_vo += -1.000000 dp_bb(j,b) u2_bbbb(b,a,i,j) m2_bbbb(i,m,e,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= sharedOps["bb_vo_322"]("a,i") * m2_bbbb_oovv("i,m,e,a");
        }

        if (include_u1_) {

            // rm1_bb_vo += 1.000000 f_bb(a,e) m1_bb(m,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += F_blks_["bb_vv"]("a,e") * m1_bb_ov("m,a");
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += -1.000000 <j,b||c,e>_abab t2_aaaa(c,a,i,j) m2_abab(i,m,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += sharedOps["abba_vvvo_16"]("a,b,e,i") * m2_abab_oovv("i,m,a,b");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // rm1_bb_vo += 0.500000 dp_bb(m,b) t2_bbbb(b,a,i,j) u0 m2_bbbb(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rm1_bb_vo("e,m") += 0.500000 * sharedOps["bbbb_vooo_216"]("a,i,j,m") * m2_bbbb_oovv("i,j,e,a") * u0;
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += -1.000000 dp_bb(j,b) u2_abab(a,b,i,j) m2_abab(i,m,a,e)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= sharedOps["aa_vo_325"]("a,i") * m2_abab_oovv("i,m,a,e");

            // rm1_bb_vo += -1.000000 <k,m||b,j>_bbbb t2_abab(a,b,i,k) m2_abab(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= sharedOps["aabb_vooo_114"]("a,i,m,j") * m2_abab_oovv("i,j,a,e");

            // rm1_bb_vo += 0.250000 <m,a||b,c>_bbbb t2_bbbb(b,c,i,j) m2_bbbb(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= 0.250000 * sharedOps["bbbb_vooo_42"]("a,i,j,m") * m2_bbbb_oovv("i,j,a,e");

            // rm1_bb_vo += -1.000000 <j,b||c,e>_abab t2_abab(c,a,j,i) m2_bbbb(i,m,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += sharedOps["bbbb_vvvo_14"]("b,e,a,i") * m2_bbbb_oovv("i,m,b,a");

            // rm1_bb_vo += -1.000000 <b,j||c,e>_abab t2_abab(c,a,i,j) m2_abab(i,m,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= sharedOps["abba_vvvo_15"]("b,e,a,i") * m2_abab_oovv("i,m,b,a");

            // rm1_bb_vo += 1.000000 dp_bb(m,b) u2_abab(a,b,i,j) m2_abab(i,j,a,e)
            // +            1.000000 dp_bb(m,b) u2_abab(a,b,j,i) m2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += 2.000000 * sharedOps["aabb_vooo_224"]("a,i,j,m") * m2_abab_oovv("i,j,a,e");
        }

        if (include_u1_) {

            // rm1_bb_vo += -1.000000 dp_bb(j,b) l2_abab(i,m,a,e) t2_abab(a,b,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= sharedOps["aa_vo_318"]("a,i") * l2_abab_oovv("i,m,a,e");

            // rm1_bb_vo += 1.000000 dp_bb(a,i) l2_bbbb(m,i,a,e)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += dp_bb_vo("a,i") * l2_bbbb_oovv("m,i,a,e");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // rm1_bb_vo += 0.500000 dp_bb(m,b) t2_abab(a,b,i,j) u0 m2_abab(i,j,a,e)
            // +            0.500000 dp_bb(m,b) t2_abab(a,b,j,i) u0 m2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rm1_bb_vo("e,m") += sharedOps["aabb_vooo_228"]("a,i,j,m") * m2_abab_oovv("i,j,a,e") * u0;
        }

        if (include_u1_) {

            // rm1_bb_vo += 0.500000 dp_bb(m,b) l2_abab(i,j,a,e) t2_abab(a,b,i,j)
            // +            0.500000 dp_bb(m,b) l2_abab(j,i,a,e) t2_abab(a,b,j,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += sharedOps["aabb_vooo_228"]("a,i,j,m") * l2_abab_oovv("i,j,a,e");

            // rm1_bb_vo += -1.000000 dp_bb(i,i) l1_bb(m,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rm1_bb_vo("e,m") -= scalar_5 * l1_bb_ov("m,e");

            // rm1_bb_vo += -1.000000 dp_aa(i,a) u1_aa(a,i) m1_bb(m,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rm1_bb_vo("e,m") -= scalar_6 * m1_bb_ov("m,e");

            // rm1_bb_vo += 1.000000 m1_bb(m,e) w0
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rm1_bb_vo("e,m") += m1_bb_ov("m,e") * w0;
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += 0.500000 <b,a||i,e>_abab m2_abab(i,m,b,a)
            // +            0.500000 <a,b||i,e>_abab m2_abab(i,m,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= V_blks_["abba_vvvo"]("b,a,e,i") * m2_abab_oovv("i,m,b,a");

            // rm1_bb_vo += -0.500000 <a,m||i,j>_abab m2_abab(i,j,a,e)
            // +            -0.500000 <a,m||j,i>_abab m2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= V_blks_["abab_vooo"]("a,m,i,j") * m2_abab_oovv("i,j,a,e");
        }

        if (include_u1_) {

            // rm1_bb_vo += -1.000000 f_bb(m,i) m1_bb(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= F_blks_["bb_oo"]("m,i") * m1_bb_ov("i,e");

            // rm1_bb_vo += -1.000000 dp_bb(i,a) u1_bb(a,i) m1_bb(m,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rm1_bb_vo("e,m") -= scalar_7 * m1_bb_ov("m,e");
        }

        if (include_u1_ && include_u2_) {

            // rm1_bb_vo += -0.250000 <a,m||b,c>_abab t2_abab(b,c,i,j) m2_abab(i,j,a,e)
            // +            -0.250000 <a,m||b,c>_abab t2_abab(b,c,j,i) m2_abab(j,i,a,e)
            // +            -0.250000 <a,m||c,b>_abab t2_abab(c,b,i,j) m2_abab(i,j,a,e)
            // +            -0.250000 <a,m||c,b>_abab t2_abab(c,b,j,i) m2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= sharedOps["aabb_vooo_41"]("a,i,m,j") * m2_abab_oovv("i,j,a,e");

            // rm1_bb_vo += 1.000000 dp_aa(j,b) u2_abab(b,a,j,i) m2_bbbb(i,m,e,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += sharedOps["bb_vo_319"]("a,i") * m2_bbbb_oovv("i,m,e,a");

            // rm1_bb_vo += 1.000000 dp_bb(m,b) u2_bbbb(b,a,i,j) m2_bbbb(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += sharedOps["bbbb_vooo_232"]("a,m,i,j") * m2_bbbb_oovv("i,j,e,a");
        }

        if (include_u1_) {

            // rm1_bb_vo += 0.500000 <j,m||a,b>_bbbb t2_bbbb(a,b,i,j) m1_bb(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") += 0.500000 * sharedOps["bb_oo_212"]("m,i") * m1_bb_ov("i,e");

            // rm1_bb_vo += -1.000000 dp_aa(a,i) l2_abab(i,m,a,e)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_bb_vo("e,m") -= dp_aa_vo("a,i") * l2_abab_oovv("i,m,a,e");
        }

        if (include_u0_ && include_u1_) {

            // rm1_bb_vo += -1.000000 dp_bb(m,e) u0 m0
            // flops: o1v1L1: 2, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rm1_bb_vo("e,m") -= dp_bb_ov("m,e") * u0 * m0;
        }

        if (include_u1_) {

            // rm1_aa_vo += 1.000000 dp_aa(a,i) l2_aaaa(m,i,a,e)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += dp_aa_vo("a,i") * l2_aaaa_oovv("m,i,a,e");
        }

        if (include_u1_ && include_u2_) {

            // rm1_aa_vo += 1.000000 dp_aa(m,b) u2_abab(b,a,i,j) m2_abab(i,j,e,a)
            // +            1.000000 dp_aa(m,b) u2_abab(b,a,j,i) m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += 2.000000 * sharedOps["baab_vooo_223"]("a,i,m,j") * m2_abab_oovv("i,j,e,a");

            // rm1_aa_vo += -1.000000 dp_aa(a,b) u1_aa(b,i) m2_aaaa(i,m,a,e)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= sharedOps["aa_vo_360"]("a,i") * m2_aaaa_oovv("i,m,a,e");

            // rm1_aa_vo += 1.000000 <m,k||b,j>_abab t2_abab(b,a,i,k) m2_abab(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += sharedOps["baab_vooo_118"]("a,m,i,j") * m2_abab_oovv("i,j,e,a");
        }

        if (include_u1_) {

            // rm1_aa_vo += -0.500000 <m,j||a,b>_abab t2_abab(a,b,i,j) m1_aa(i,e)
            // +            -0.500000 <m,j||b,a>_abab t2_abab(b,a,i,j) m1_aa(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= sharedOps["aa_oo_207"]("m,i") * m1_aa_ov("i,e");
        }

        if (include_u1_ && include_u2_) {

            // rm1_aa_vo += 0.500000 <b,a||e,i>_aaaa m2_aaaa(m,i,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += 0.500000 * V_blks_["aaaa_vvvo"]("b,a,e,i") * m2_aaaa_oovv("m,i,b,a");

            // rm1_aa_vo += -1.000000 <b,j||e,c>_abab t2_bbbb(c,a,i,j) m2_abab(m,i,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= sharedOps["aabb_vvvo_12"]("b,e,a,i") * m2_abab_oovv("m,i,b,a");
        }

        if (include_u1_) {

            // rm1_aa_vo += 0.500000 dp_aa(m,b) l2_aaaa(i,j,e,a) t2_aaaa(b,a,i,j)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += 0.500000 * sharedOps["aaaa_vooo_215"]("a,m,i,j") * l2_aaaa_oovv("i,j,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rm1_aa_vo += -1.000000 <j,b||e,c>_abab t2_abab(a,c,j,i) m2_abab(m,i,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += sharedOps["aabb_vvvo_18"]("e,a,b,i") * m2_abab_oovv("m,i,a,b");
        }

        if (include_u1_) {

            // rm1_aa_vo += 0.500000 <j,m||a,b>_aaaa t2_aaaa(a,b,i,j) m1_aa(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += 0.500000 * sharedOps["aa_oo_211"]("i,m") * m1_aa_ov("i,e");

            // rm1_aa_vo += 2.000000 dp_aa(m,a) u1_aa(a,i) m1_aa(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += 2.000000 * sharedOps["aa_oo_312"]("i,m") * m1_aa_ov("i,e");

            // rm1_aa_vo += -1.000000 dp_bb(i,i) l1_aa(m,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rm1_aa_vo("e,m") -= scalar_5 * l1_aa_ov("m,e");

            // rm1_aa_vo += 1.000000 dp_bb(j,b) l2_aaaa(i,m,e,a) t2_abab(a,b,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += sharedOps["aa_vo_318"]("a,i") * l2_aaaa_oovv("i,m,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rm1_aa_vo += 1.000000 dp_bb(j,b) u2_abab(a,b,i,j) m2_aaaa(i,m,e,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += sharedOps["aa_vo_325"]("a,i") * m2_aaaa_oovv("i,m,e,a");
        }

        if (include_u1_) {

            // rm1_aa_vo += -1.000000 dp_aa(i,a) u1_aa(a,i) m1_aa(m,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rm1_aa_vo("e,m") -= scalar_6 * m1_aa_ov("m,e");
        }

        if (include_u1_ && include_u2_) {

            // rm1_aa_vo += 1.000000 <m,k||j,b>_abab t2_abab(a,b,i,k) m2_aaaa(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= sharedOps["aaaa_vooo_115"]("a,m,j,i") * m2_aaaa_oovv("i,j,e,a");

            // rm1_aa_vo += 1.000000 dp_bb(j,i) u1_bb(a,j) m2_abab(m,i,e,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += sharedOps["bb_vo_365"]("a,i") * m2_abab_oovv("m,i,e,a");

            // rm1_aa_vo += -1.000000 dp_aa(j,b) u2_aaaa(b,a,i,j) m2_aaaa(i,m,e,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= sharedOps["aa_vo_324"]("a,i") * m2_aaaa_oovv("i,m,e,a");

            // rm1_aa_vo += -1.000000 dp_bb(a,b) u1_bb(b,i) m2_abab(m,i,e,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= sharedOps["bb_vo_361"]("a,i") * m2_abab_oovv("m,i,e,a");

            // rm1_aa_vo += -1.000000 dp_aa(j,b) u2_abab(b,a,j,i) m2_abab(m,i,e,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= sharedOps["bb_vo_319"]("a,i") * m2_abab_oovv("m,i,e,a");

            // rm1_aa_vo += -1.000000 <j,b||e,c>_aaaa t2_aaaa(c,a,i,j) m2_aaaa(i,m,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += sharedOps["aaaa_vvvo_19"]("a,b,e,i") * m2_aaaa_oovv("i,m,b,a");
        }

        if (include_u1_) {

            // rm1_aa_vo += -1.000000 dp_bb(i,a) u1_bb(a,i) m1_aa(m,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rm1_aa_vo("e,m") -= scalar_7 * m1_aa_ov("m,e");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // rm1_aa_vo += 0.500000 dp_aa(m,b) t2_aaaa(b,a,i,j) u0 m2_aaaa(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rm1_aa_vo("e,m") += 0.500000 * sharedOps["aaaa_vooo_215"]("a,m,i,j") * m2_aaaa_oovv("i,j,e,a") * u0;
        }

        if (include_u1_ && include_u2_) {

            // rm1_aa_vo += -0.500000 <m,a||i,j>_abab m2_abab(i,j,e,a)
            // +            -0.500000 <m,a||j,i>_abab m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += V_blks_["baab_vooo"]("a,m,i,j") * m2_abab_oovv("i,j,e,a");
        }

        if (include_u1_) {

            // rm1_aa_vo += -1.000000 dp_aa(i,i) l1_aa(m,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rm1_aa_vo("e,m") -= scalar_4 * l1_aa_ov("m,e");
        }

        if (include_u1_ && include_u2_) {

            // rm1_aa_vo += 1.000000 dp_aa(m,b) u2_aaaa(b,a,i,j) m2_aaaa(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += sharedOps["aaaa_vooo_234"]("a,i,j,m") * m2_aaaa_oovv("i,j,e,a");

            // rm1_aa_vo += -1.000000 <k,m||b,j>_aaaa t2_abab(b,a,k,i) m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= sharedOps["baab_vooo_117"]("a,m,j,i") * m2_abab_oovv("j,i,e,a");
        }

        if (include_u1_) {

            // rm1_aa_vo += 1.000000 m1_aa(m,e) w0
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rm1_aa_vo("e,m") += m1_aa_ov("m,e") * w0;
        }

        if (include_u1_ && include_u2_) {

            // rm1_aa_vo += 0.250000 <k,j||e,i>_abab t2_abab(b,a,k,j) m2_abab(m,i,b,a)
            // +            0.250000 <j,k||e,i>_abab t2_abab(b,a,j,k) m2_abab(m,i,b,a)
            // +            0.250000 <k,j||e,i>_abab t2_abab(a,b,k,j) m2_abab(m,i,a,b)
            // +            0.250000 <j,k||e,i>_abab t2_abab(a,b,j,k) m2_abab(m,i,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += sharedOps["aabb_vvvo_53"]("e,b,a,i") * m2_abab_oovv("m,i,b,a");
        }

        if (include_u0_ && include_u1_) {

            // rm1_aa_vo += 1.000000 f_aa(m,e) m0
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += F_blks_["aa_ov"]("m,e") * m0;
        }

        if (include_u1_ && include_u2_) {

            // rm1_aa_vo += 1.000000 dp_aa(j,i) u1_aa(a,j) m2_aaaa(m,i,e,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += sharedOps["aa_vo_364"]("a,i") * m2_aaaa_oovv("m,i,e,a");
        }

        if (include_u1_) {

            // rm1_aa_vo += 1.000000 <m,a||e,i>_aaaa m1_aa(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= V_blks_["aaaa_vovo"]("a,m,e,i") * m1_aa_ov("i,a");

            // rm1_aa_vo += 1.000000 f_aa(a,e) m1_aa(m,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += F_blks_["aa_vv"]("a,e") * m1_aa_ov("m,a");
        }

        if (include_u0_ && include_u1_) {

            // rm1_aa_vo += -1.000000 dp_aa(m,e) u0 m0
            // flops: o1v1L1: 2, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rm1_aa_vo("e,m") -= dp_aa_ov("m,e") * u0 * m0;
        }

        if (include_u1_) {

            // rm1_aa_vo += -0.500000 <j,i||e,b>_abab t2_abab(a,b,j,i) m1_aa(m,a)
            // +            -0.500000 <i,j||e,b>_abab t2_abab(a,b,i,j) m1_aa(m,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= sharedOps["aa_vv_150"]("a,e") * m1_aa_ov("m,a");

            // rm1_aa_vo += 0.500000 <j,i||e,b>_aaaa t2_aaaa(b,a,j,i) m1_aa(m,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += 0.500000 * sharedOps["aa_vv_151"]("a,e") * m1_aa_ov("m,a");
        }

        if (include_u1_ && include_u2_) {

            // rm1_aa_vo += 0.500000 <m,a||i,j>_aaaa m2_aaaa(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= 0.500000 * V_blks_["aaaa_vooo"]("a,m,i,j") * m2_aaaa_oovv("i,j,a,e");
        }

        if (include_u1_) {

            // rm1_aa_vo += -1.000000 f_aa(m,i) m1_aa(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= F_blks_["aa_oo"]("m,i") * m1_aa_ov("i,e");

            // rm1_aa_vo += 0.500000 dp_aa(m,b) l2_abab(i,j,e,a) t2_abab(b,a,i,j)
            // +            0.500000 dp_aa(m,b) l2_abab(j,i,e,a) t2_abab(b,a,j,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += sharedOps["baab_vooo_226"]("a,i,m,j") * l2_abab_oovv("i,j,e,a");

            // rm1_aa_vo += -1.000000 dp_aa(j,b) l2_abab(m,i,e,a) t2_abab(b,a,j,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= sharedOps["bb_vo_321"]("a,i") * l2_abab_oovv("m,i,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rm1_aa_vo += -1.000000 <j,b||e,c>_aaaa t2_abab(c,a,j,i) m2_abab(m,i,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += sharedOps["aabb_vvvo_13"]("b,e,a,i") * m2_abab_oovv("m,i,b,a");
        }

        if (include_u1_) {

            // rm1_aa_vo += -1.000000 dp_aa(m,e)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rm1_aa_vo("e,m") -= dp_aa_ov("m,e");

            // rm1_aa_vo += 1.000000 dp_bb(j,b) l2_abab(m,i,e,a) t2_bbbb(b,a,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += sharedOps["bb_vo_323"]("a,i") * l2_abab_oovv("m,i,e,a");

            // rm1_aa_vo += -1.000000 dp_aa(j,b) l2_aaaa(i,m,e,a) t2_aaaa(b,a,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= sharedOps["aa_vo_320"]("a,i") * l2_aaaa_oovv("i,m,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rm1_aa_vo += -0.500000 f_aa(m,b) t2_abab(b,a,i,j) m2_abab(i,j,e,a)
            // +            -0.500000 f_aa(m,b) t2_abab(b,a,j,i) m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= sharedOps["baab_vooo_259"]("a,m,i,j") * m2_abab_oovv("i,j,e,a");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // rm1_aa_vo += 0.500000 dp_aa(m,b) t2_abab(b,a,i,j) u0 m2_abab(i,j,e,a)
            // +            0.500000 dp_aa(m,b) t2_abab(b,a,j,i) u0 m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rm1_aa_vo("e,m") += sharedOps["baab_vooo_226"]("a,i,m,j") * m2_abab_oovv("i,j,e,a") * u0;
        }

        if (include_u1_) {

            // rm1_aa_vo += -1.000000 dp_bb(a,i) l2_abab(m,i,e,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= dp_bb_vo("a,i") * l2_abab_oovv("m,i,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rm1_aa_vo += 1.000000 dp_bb(j,b) u2_bbbb(b,a,i,j) m2_abab(m,i,e,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += sharedOps["bb_vo_322"]("a,i") * m2_abab_oovv("m,i,e,a");
        }

        if (include_u1_) {

            // rm1_aa_vo += 1.000000 <m,a||e,i>_abab m1_bb(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= V_blks_["baab_vovo"]("a,m,e,i") * m1_bb_ov("i,a");
        }

        if (include_u1_ && include_u2_) {

            // rm1_aa_vo += 0.250000 <k,j||e,i>_aaaa t2_aaaa(b,a,k,j) m2_aaaa(m,i,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += 0.250000 * sharedOps["aaaa_vvvo_56"]("e,b,a,i") * m2_aaaa_oovv("m,i,b,a");

            // rm1_aa_vo += -1.000000 <b,j||e,c>_abab t2_abab(a,c,i,j) m2_aaaa(i,m,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= sharedOps["aaaa_vvvo_17"]("b,e,a,i") * m2_aaaa_oovv("i,m,b,a");

            // rm1_aa_vo += 0.250000 <m,a||b,c>_aaaa t2_aaaa(b,c,i,j) m2_aaaa(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= 0.250000 * sharedOps["aaaa_vooo_43"]("a,m,i,j") * m2_aaaa_oovv("i,j,a,e");

            // rm1_aa_vo += -0.250000 <m,a||b,c>_abab t2_abab(b,c,i,j) m2_abab(i,j,e,a)
            // +            -0.250000 <m,a||b,c>_abab t2_abab(b,c,j,i) m2_abab(j,i,e,a)
            // +            -0.250000 <m,a||c,b>_abab t2_abab(c,b,i,j) m2_abab(i,j,e,a)
            // +            -0.250000 <m,a||c,b>_abab t2_abab(c,b,j,i) m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += sharedOps["baab_vooo_40"]("a,m,i,j") * m2_abab_oovv("i,j,e,a");

            // rm1_aa_vo += 0.500000 <b,a||e,i>_abab m2_abab(m,i,b,a)
            // +            0.500000 <a,b||e,i>_abab m2_abab(m,i,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") += V_blks_["abab_vvvo"]("b,a,e,i") * m2_abab_oovv("m,i,b,a");

            // rm1_aa_vo += -1.000000 <k,m||b,j>_aaaa t2_aaaa(b,a,i,k) m2_aaaa(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= sharedOps["aaaa_vooo_119"]("a,m,j,i") * m2_aaaa_oovv("i,j,e,a");

            // rm1_aa_vo += -0.500000 f_aa(m,b) t2_aaaa(b,a,i,j) m2_aaaa(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= 0.500000 * sharedOps["aaaa_vooo_240"]("a,m,i,j") * m2_aaaa_oovv("i,j,e,a");

            // rm1_aa_vo += 1.000000 <m,k||j,b>_abab t2_bbbb(b,a,i,k) m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rm1_aa_vo("e,m") -= sharedOps["baab_vooo_111"]("a,m,j,i") * m2_abab_oovv("j,i,e,a");
        }

        if (include_u0_) {

            // crm0 += 1.000000
            // flops: o0v0: 1 | mem: o0v0: 1,
            crm0 += 1.000000;
        }

        if (include_u0_ && include_u1_) {

            // rm0 += -1.000000 dp_bb(i,a) u1_bb(a,i) m0
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            rm0 -= scalar_7 * m0;
        }

        if (include_u0_) {

            // rm0 += -1.000000 dp_bb(i,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_5;

            // rm0 += -1.000000 dp_aa(i,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            rm0 -= scalar_4;
        }

        if (include_u0_ && include_u1_) {

            // rm0 += -1.000000 dp_aa(i,a) u1_aa(a,i) m0
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            rm0 -= scalar_6 * m0;
        }

        if (include_u0_) {

            // rm0 += 1.000000 m0 w0
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            rm0 += m0 * w0;
        }

        if (include_u2_) {

            // rm2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rm2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <j,n||b,e>_abab u2_abab(b,a,j,i) m2_bbbb(i,m,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = sharedOps["bbbb_vvoo_49"]("a,e,i,n") * m2_bbbb_oovv("i,m,f,a");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <j,n||b,e>_abab u2_aaaa(b,a,i,j) m2_abab(i,m,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["abab_vvoo_50"]("a,e,i,n") * m2_abab_oovv("i,m,a,f");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <j,n||e,b>_bbbb u2_abab(a,b,i,j) m2_abab(i,m,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["abab_vvoo_48"]("a,e,i,n") * m2_abab_oovv("i,m,a,f");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) <n,a||e,b>_bbbb u1_bb(b,i) m2_bbbb(i,m,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["bbbb_vvoo_166"]("a,e,i,n") * m2_bbbb_oovv("i,m,a,f");
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <i,n||e,a>_bbbb u1_bb(a,i) m1_bb(m,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["bb_vo_326"]("e,n") * m1_bb_ov("m,f");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <n,a||e,i>_bbbb l2_bbbb(m,i,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= V_blks_["bbbb_vovo"]("a,n,e,i") * l2_bbbb_oovv("m,i,a,f");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <j,n||b,e>_abab l2_abab(i,m,a,f) t2_aaaa(b,a,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["abab_vvoo_0"]("a,e,i,n") * l2_abab_oovv("i,m,a,f");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <j,n||e,b>_bbbb l2_bbbb(i,m,f,a) t2_bbbb(b,a,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["bbbb_vvoo_2"]("a,e,i,n") * l2_bbbb_oovv("i,m,f,a");
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <j,n||e,b>_bbbb u2_bbbb(b,a,i,j) m2_bbbb(i,m,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["bbbb_vvoo_45"]("e,a,n,i") * m2_bbbb_oovv("i,m,f,a");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) <a,n||b,e>_abab u1_aa(b,i) m2_abab(i,m,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= sharedOps["abab_vvoo_167"]("a,e,i,n") * m2_abab_oovv("i,m,a,f");
        }

        {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) f_bb(n,e) l1_bb(m,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= F_blks_["bb_ov"]("n,e") * l1_bb_ov("m,f");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <j,n||e,b>_bbbb l2_abab(i,m,a,f) t2_abab(a,b,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["abab_vvoo_3"]("a,e,i,n") * l2_abab_oovv("i,m,a,f");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <j,n||b,e>_abab l2_bbbb(i,m,f,a) t2_abab(b,a,j,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["bbbb_vvoo_1"]("e,a,n,i") * l2_bbbb_oovv("i,m,f,a");
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) <i,n||a,e>_abab u1_aa(a,i) m1_bb(m,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= sharedOps["bb_vo_36"]("e,n") * m1_bb_ov("m,f");
        }

        {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) <a,n||i,e>_abab l2_abab(i,m,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += V_blks_["abba_vovo"]("a,n,e,i") * l2_abab_oovv("i,m,a,f");

            // rl2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("f,e,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) f_bb(n,i) l2_bbbb(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * F_blks_["bb_oo"]("n,i") * l2_bbbb_oovv("m,i,e,f");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) <j,n||a,i>_abab u1_aa(a,j) m2_bbbb(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= sharedOps["bb_oo_346"]("n,i") * m2_bbbb_oovv("m,i,e,f");
        }

        {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) <n,a||e,f>_bbbb l1_bb(m,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += V_blks_["bbbb_vovv"]("a,n,e,f") * l1_bb_ov("m,a");
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) dp_bb(n,a) l2_bbbb(i,m,e,f) u1_bb(a,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= sharedOps["bb_oo_313"]("n,i") * l2_bbbb_oovv("i,m,e,f");
        }

        {

            // tempPerm_bbbb_vvoo += -0.500000 P(m,n) <j,n||a,b>_bbbb l2_bbbb(i,m,e,f) t2_bbbb(a,b,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= 0.500000 * sharedOps["bb_oo_212"]("n,i") * l2_bbbb_oovv("i,m,e,f");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) f_bb(n,a) u1_bb(a,i) m2_bbbb(i,m,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["bb_oo_366"]("i,n") * m2_bbbb_oovv("i,m,e,f");
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += -0.500000 P(m,n) <j,n||a,b>_bbbb u2_bbbb(a,b,i,j) m2_bbbb(i,m,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= 0.500000 * sharedOps["bb_oo_242"]("i,n") * m2_bbbb_oovv("i,m,e,f");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) dp_bb(n,a) u0 u1_bb(a,i) m2_bbbb(i,m,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 2, o0v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") -= sharedOps["bb_oo_313"]("n,i") * u0 * m2_bbbb_oovv("i,m,e,f");
        }

        {

            // tempPerm_bbbb_vvoo += 0.500000 P(m,n) <j,n||a,b>_abab l2_bbbb(i,m,e,f) t2_abab(a,b,j,i)
            // +                0.500000 P(m,n) <j,n||b,a>_abab l2_bbbb(i,m,e,f) t2_abab(b,a,j,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["bb_oo_208"]("n,i") * l2_bbbb_oovv("i,m,e,f");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) <j,n||a,i>_bbbb u1_bb(a,j) m2_bbbb(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= sharedOps["bb_oo_344"]("n,i") * m2_bbbb_oovv("m,i,e,f");
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += 0.500000 P(m,n) <j,n||a,b>_abab u2_abab(a,b,j,i) m2_bbbb(i,m,e,f)
            // +                0.500000 P(m,n) <j,n||b,a>_abab u2_abab(b,a,j,i) m2_bbbb(i,m,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["bb_oo_214"]("i,n") * m2_bbbb_oovv("i,m,e,f");
        }

        {

            // rl2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");

            // tempPerm_bbbb_vvoo += -1.000000 P(e,f) <m,n||e,i>_bbbb l1_bb(i,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * V_blks_["bbbb_oovo"]("m,n,e,i") * l1_bb_ov("i,f");

            // tempPerm_bbbb_vvoo += 0.500000 P(e,f) <j,i||b,e>_abab l2_bbbb(m,n,f,a) t2_abab(b,a,j,i)
            // +                0.500000 P(e,f) <i,j||b,e>_abab l2_bbbb(m,n,f,a) t2_abab(b,a,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["bb_vv_149"]("e,a") * l2_bbbb_oovv("m,n,f,a");

            // tempPerm_bbbb_vvoo += 1.000000 P(e,f) f_bb(a,e) l2_bbbb(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += F_blks_["bb_vv"]("a,e") * l2_bbbb_oovv("m,n,a,f");

            // tempPerm_bbbb_vvoo += -0.500000 P(e,f) <j,i||e,b>_bbbb l2_bbbb(m,n,f,a) t2_bbbb(b,a,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= 0.500000 * sharedOps["bb_vv_152"]("e,a") * l2_bbbb_oovv("m,n,f,a");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(e,f) <i,a||b,e>_abab u1_aa(b,i) m2_bbbb(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= sharedOps["bb_vv_315"]("a,e") * m2_bbbb_oovv("m,n,a,f");

            // tempPerm_bbbb_vvoo += -1.000000 P(e,f) <i,a||e,b>_bbbb u1_bb(b,i) m2_bbbb(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["bb_vv_314"]("a,e") * m2_bbbb_oovv("m,n,a,f");
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(e,f) <m,n||e,a>_bbbb u1_bb(a,i) m1_bb(i,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= sharedOps["bbbb_vooo_148"]("e,i,m,n") * m1_bb_ov("i,f");
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += 0.500000 P(e,f) <j,i||b,e>_abab u2_abab(b,a,j,i) m2_bbbb(m,n,f,a)
            // +                0.500000 P(e,f) <i,j||b,e>_abab u2_abab(b,a,i,j) m2_bbbb(m,n,f,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += sharedOps["bb_vv_154"]("a,e") * m2_bbbb_oovv("m,n,f,a");

            // tempPerm_bbbb_vvoo += -0.500000 P(e,f) <j,i||e,b>_bbbb u2_bbbb(b,a,j,i) m2_bbbb(m,n,f,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= 0.500000 * sharedOps["bb_vv_169"]("a,e") * m2_bbbb_oovv("m,n,f,a");
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(e,f) dp_bb(i,e) l2_bbbb(m,n,f,a) u1_bb(a,i)
            // flops: o3v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o3v1L1: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") -= u1_bb_vo("a,i") * l2_bbbb_oovv("m,n,f,a") * dp_bb_ov("i,e");
        }

        {

            // rl2_bbbb_vvoo += 1.000000 <m,n||e,f>_bbbb
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_bbbb_vvoo("e,f,m,n") += V_blks_["bbbb_oovv"]("m,n,e,f");

            // rl2_bbbb_vvoo += 0.250000 <j,i||e,f>_bbbb l2_bbbb(m,n,b,a) t2_bbbb(b,a,j,i)
            // flops: o4v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o4v0L1: 1,
            rl2_bbbb_vvoo("e,f,m,n") += 0.250000 * l2_bbbb_oovv("m,n,b,a") * t2_bbbb_vvoo("b,a,j,i") * V_blks_["bbbb_oovv"]("j,i,e,f");

            // rl2_bbbb_vvoo += 0.250000 <m,n||a,b>_bbbb l2_bbbb(i,j,e,f) t2_bbbb(a,b,i,j)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_bbbb_vvoo("e,f,m,n") += 0.250000 * sharedOps["bbbb_oooo_107"]("i,j,m,n") * l2_bbbb_oovv("i,j,e,f");

            // rl2_bbbb_vvoo += 0.500000 <b,a||e,f>_bbbb l2_bbbb(m,n,b,a)
            // flops: o2v4L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_bbbb_vvoo("e,f,m,n") += 0.500000 * V_blks_["bbbb_vvvv"]("b,a,e,f") * l2_bbbb_oovv("m,n,b,a");
        }

        if (include_u1_ && include_u2_) {

            // rl2_bbbb_vvoo += 1.000000 <m,n||a,j>_bbbb u1_bb(a,i) m2_bbbb(i,j,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_bbbb_vvoo("e,f,m,n") += sharedOps["bbbb_oooo_307"]("i,m,n,j") * m2_bbbb_oovv("i,j,e,f");
        }

        if (include_u2_) {

            // rl2_bbbb_vvoo += 0.250000 <m,n||a,b>_bbbb u2_bbbb(a,b,i,j) m2_bbbb(i,j,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_bbbb_vvoo("e,f,m,n") += 0.250000 * sharedOps["bbbb_oooo_135"]("i,j,m,n") * m2_bbbb_oovv("i,j,e,f");

            // rl2_bbbb_vvoo += 0.250000 <j,i||e,f>_bbbb u2_bbbb(b,a,j,i) m2_bbbb(m,n,b,a)
            // flops: o4v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o4v0L1: 1,
            rl2_bbbb_vvoo("e,f,m,n") += 0.250000 * u2_bbbb_vvoo("b,a,j,i") * m2_bbbb_oovv("m,n,b,a") * V_blks_["bbbb_oovv"]("j,i,e,f");
        }

        {

            // rl2_bbbb_vvoo += 0.500000 <m,n||i,j>_bbbb l2_bbbb(i,j,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_bbbb_vvoo("e,f,m,n") += 0.500000 * V_blks_["bbbb_oooo"]("m,n,i,j") * l2_bbbb_oovv("i,j,e,f");

            // rl2_abab_vvoo += 1.000000 <m,n||e,f>_abab
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_abab_vvoo("e,f,m,n") += V_blks_["abab_oovv"]("m,n,e,f");

            // rl2_abab_vvoo += 1.000000 <j,n||b,f>_abab l2_abab(m,i,e,a) t2_abab(b,a,j,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["bbbb_vvoo_1"]("f,a,n,i") * l2_abab_oovv("m,i,e,a");
        }

        if (include_u2_) {

            // rl2_abab_vvoo += -0.500000 <m,j||a,b>_abab u2_abab(a,b,i,j) m2_abab(i,n,e,f)
            // +                -0.500000 <m,j||b,a>_abab u2_abab(b,a,i,j) m2_abab(i,n,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["aa_oo_213"]("m,i") * m2_abab_oovv("i,n,e,f");

            // rl2_abab_vvoo += 1.000000 <m,j||b,f>_abab u2_abab(b,a,i,j) m2_abab(i,n,e,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["bbaa_vvoo_95"]("f,a,m,i") * m2_abab_oovv("i,n,e,a");

            // rl2_abab_vvoo += 0.500000 <j,i||f,b>_bbbb u2_bbbb(b,a,j,i) m2_abab(m,n,e,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += 0.500000 * sharedOps["bb_vv_169"]("a,f") * m2_abab_oovv("m,n,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rl2_abab_vvoo += -1.000000 f_bb(n,a) u1_bb(a,i) m2_abab(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["bb_oo_366"]("i,n") * m2_abab_oovv("m,i,e,f");
        }

        if (include_u1_) {

            // rl2_abab_vvoo += -1.000000 dp_aa(m,e) m1_bb(n,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= dp_aa_ov("m,e") * m1_bb_ov("n,f");
        }

        {

            // rl2_abab_vvoo += 1.000000 f_aa(a,e) l2_abab(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += F_blks_["aa_vv"]("a,e") * l2_abab_oovv("m,n,a,f");

            // rl2_abab_vvoo += -1.000000 <m,n||i,f>_abab l1_aa(i,e)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += V_blks_["abba_oovo"]("m,n,f,i") * l1_aa_ov("i,e");
        }

        if (include_u1_ && include_u2_) {

            // rl2_abab_vvoo += -1.000000 <j,n||a,i>_bbbb u1_bb(a,j) m2_abab(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["bb_oo_344"]("n,i") * m2_abab_oovv("m,i,e,f");
        }

        {

            // rl2_abab_vvoo += 0.500000 <j,i||f,b>_bbbb l2_abab(m,n,e,a) t2_bbbb(b,a,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += 0.500000 * sharedOps["bb_vv_152"]("f,a") * l2_abab_oovv("m,n,e,a");

            // rl2_abab_vvoo += -0.500000 <m,j||a,b>_abab l2_abab(i,n,e,f) t2_abab(a,b,i,j)
            // +                -0.500000 <m,j||b,a>_abab l2_abab(i,n,e,f) t2_abab(b,a,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["aa_oo_207"]("m,i") * l2_abab_oovv("i,n,e,f");

            // rl2_abab_vvoo += 1.000000 <n,a||f,i>_bbbb l2_abab(m,i,e,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= V_blks_["bbbb_vovo"]("a,n,f,i") * l2_abab_oovv("m,i,e,a");

            // rl2_abab_vvoo += 1.000000 <j,n||f,b>_bbbb l2_aaaa(i,m,e,a) t2_abab(a,b,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["abab_vvoo_3"]("a,f,i,n") * l2_aaaa_oovv("i,m,e,a");
        }

        if (include_u2_) {

            // rl2_abab_vvoo += -0.500000 <j,i||b,f>_abab u2_abab(b,a,j,i) m2_abab(m,n,e,a)
            // +                -0.500000 <i,j||b,f>_abab u2_abab(b,a,i,j) m2_abab(m,n,e,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["bb_vv_154"]("a,f") * m2_abab_oovv("m,n,e,a");
        }

        if (include_u1_) {

            // rl2_abab_vvoo += 1.000000 <m,i||e,a>_abab u1_bb(a,i) m1_bb(n,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["aa_vo_37"]("e,m") * m1_bb_ov("n,f");
        }

        if (include_u2_) {

            // rl2_abab_vvoo += 1.000000 <j,n||e,b>_abab u2_abab(a,b,j,i) m2_abab(m,i,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["aabb_vvoo_93"]("a,e,i,n") * m2_abab_oovv("m,i,a,f");
        }

        {

            // rl2_abab_vvoo += 1.000000 <m,j||e,b>_abab l2_abab(i,n,a,f) t2_abab(a,b,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["aaaa_vvoo_6"]("a,e,i,m") * l2_abab_oovv("i,n,a,f");

            // rl2_abab_vvoo += 0.500000 <j,m||a,b>_aaaa l2_abab(i,n,e,f) t2_aaaa(a,b,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += 0.500000 * sharedOps["aa_oo_211"]("i,m") * l2_abab_oovv("i,n,e,f");
        }

        if (include_u2_) {

            // rl2_abab_vvoo += 1.000000 <j,n||b,f>_abab u2_abab(b,a,j,i) m2_abab(m,i,e,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["bbbb_vvoo_49"]("a,f,i,n") * m2_abab_oovv("m,i,e,a");
        }

        if (include_u1_) {

            // rl2_abab_vvoo += -1.000000 <m,n||e,a>_abab u1_bb(a,i) m1_bb(i,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["aabb_vooo_145"]("e,m,i,n") * m1_bb_ov("i,f");

            // rl2_abab_vvoo += 1.000000 dp_bb(i,f) l2_abab(m,n,e,a) u1_bb(a,i)
            // flops: o3v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o3v1L1: 1,
            rl2_abab_vvoo("e,f,m,n") += l2_abab_oovv("m,n,e,a") * u1_bb_vo("a,i") * dp_bb_ov("i,f");

            // rl2_abab_vvoo += -1.000000 <i,m||e,a>_aaaa u1_aa(a,i) m1_bb(n,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["aa_vo_327"]("e,m") * m1_bb_ov("n,f");
        }

        {

            // rl2_abab_vvoo += 1.000000 <m,j||e,b>_abab l2_bbbb(i,n,f,a) t2_bbbb(b,a,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["abab_vvoo_7"]("e,a,m,i") * l2_bbbb_oovv("i,n,f,a");
        }

        if (include_u0_) {

            // rl2_abab_vvoo += -1.000000 dp_aa(m,e) l1_bb(n,f) u0
            // flops: o2v2L1: 2, o0v0: 1 | mem: o2v2L1: 2, o0v0: 1,
            rl2_abab_vvoo("e,f,m,n") -= dp_aa_ov("m,e") * l1_bb_ov("n,f") * u0;
        }

        if (include_u2_) {

            // rl2_abab_vvoo += 1.000000 <j,n||b,f>_abab u2_aaaa(b,a,i,j) m2_aaaa(i,m,e,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["abab_vvoo_50"]("a,f,i,n") * m2_aaaa_oovv("i,m,e,a");
        }

        {

            // rl2_abab_vvoo += 1.000000 f_aa(m,e) l1_bb(n,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += F_blks_["aa_ov"]("m,e") * l1_bb_ov("n,f");
        }

        if (include_u1_ && include_u2_) {

            // rl2_abab_vvoo += 1.000000 <m,n||j,a>_abab u1_bb(a,i) m2_abab(j,i,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["aabb_oooo_309"]("m,j,n,i") * m2_abab_oovv("j,i,e,f");

            // rl2_abab_vvoo += -1.000000 <m,j||i,a>_abab u1_bb(a,j) m2_abab(i,n,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["aa_oo_345"]("m,i") * m2_abab_oovv("i,n,e,f");
        }

        if (include_u1_) {

            // rl2_abab_vvoo += 1.000000 <i,n||a,f>_abab u1_aa(a,i) m1_aa(m,e)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["bb_vo_36"]("f,n") * m1_aa_ov("m,e");
        }

        {

            // rl2_abab_vvoo += 1.000000 <j,n||b,f>_abab l2_aaaa(i,m,e,a) t2_aaaa(b,a,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["abab_vvoo_0"]("a,f,i,n") * l2_aaaa_oovv("i,m,e,a");

            // rl2_abab_vvoo += 0.250000 <j,i||e,f>_abab l2_abab(m,n,b,a) t2_abab(b,a,j,i)
            // +                0.250000 <i,j||e,f>_abab l2_abab(m,n,b,a) t2_abab(b,a,i,j)
            // +                0.250000 <j,i||e,f>_abab l2_abab(m,n,a,b) t2_abab(a,b,j,i)
            // +                0.250000 <i,j||e,f>_abab l2_abab(m,n,a,b) t2_abab(a,b,i,j)
            // flops: o4v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o4v0L1: 1,
            rl2_abab_vvoo("e,f,m,n") += l2_abab_oovv("m,n,b,a") * t2_abab_vvoo("b,a,j,i") * V_blks_["abab_oovv"]("j,i,e,f");

            // rl2_abab_vvoo += 0.250000 <m,n||a,b>_abab l2_abab(i,j,e,f) t2_abab(a,b,i,j)
            // +                0.250000 <m,n||b,a>_abab l2_abab(i,j,e,f) t2_abab(b,a,i,j)
            // +                0.250000 <m,n||a,b>_abab l2_abab(j,i,e,f) t2_abab(a,b,j,i)
            // +                0.250000 <m,n||b,a>_abab l2_abab(j,i,e,f) t2_abab(b,a,j,i)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["aabb_oooo_103"]("m,i,n,j") * l2_abab_oovv("i,j,e,f");

            // rl2_abab_vvoo += 1.000000 <m,a||e,i>_aaaa l2_abab(i,n,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= V_blks_["aaaa_vovo"]("a,m,e,i") * l2_abab_oovv("i,n,a,f");
        }

        if (include_u2_) {

            // rl2_abab_vvoo += 1.000000 <m,j||e,b>_abab u2_abab(a,b,i,j) m2_abab(i,n,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["aaaa_vvoo_55"]("a,e,i,m") * m2_abab_oovv("i,n,a,f");
        }

        {

            // rl2_abab_vvoo += -1.000000 <m,a||i,f>_abab l2_abab(i,n,e,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= V_blks_["baba_vovo"]("a,m,f,i") * l2_abab_oovv("i,n,e,a");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // rl2_abab_vvoo += 1.000000 dp_aa(m,a) u0 u1_aa(a,i) m2_abab(i,n,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 2, o0v0: 1,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["aa_oo_312"]("i,m") * m2_abab_oovv("i,n,e,f") * u0;
        }

        if (include_u1_ && include_u2_) {

            // rl2_abab_vvoo += 1.000000 <m,a||e,b>_abab u1_bb(b,i) m2_bbbb(i,n,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["abab_vvoo_163"]("e,a,m,i") * m2_bbbb_oovv("i,n,a,f");
        }

        if (include_u1_) {

            // rl2_abab_vvoo += -1.000000 <m,n||a,f>_abab u1_aa(a,i) m1_aa(i,e)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["baab_vooo_146"]("f,i,m,n") * m1_aa_ov("i,e");
        }

        {

            // rl2_abab_vvoo += 1.000000 <j,m||e,b>_aaaa l2_abab(i,n,a,f) t2_aaaa(b,a,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["aaaa_vvoo_5"]("e,a,m,i") * l2_abab_oovv("i,n,a,f");
        }

        if (include_u1_ && include_u2_) {

            // rl2_abab_vvoo += -1.000000 <m,a||b,f>_abab u1_aa(b,i) m2_abab(i,n,e,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["bbaa_vvoo_201"]("a,f,i,m") * m2_abab_oovv("i,n,e,a");
        }

        {

            // rl2_abab_vvoo += -0.500000 <j,i||b,f>_abab l2_abab(m,n,e,a) t2_abab(b,a,j,i)
            // +                -0.500000 <i,j||b,f>_abab l2_abab(m,n,e,a) t2_abab(b,a,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["bb_vv_149"]("f,a") * l2_abab_oovv("m,n,e,a");
        }

        if (include_u2_) {

            // rl2_abab_vvoo += 0.250000 <m,n||a,b>_abab u2_abab(a,b,i,j) m2_abab(i,j,e,f)
            // +                0.250000 <m,n||a,b>_abab u2_abab(a,b,j,i) m2_abab(j,i,e,f)
            // +                0.250000 <m,n||b,a>_abab u2_abab(b,a,i,j) m2_abab(i,j,e,f)
            // +                0.250000 <m,n||b,a>_abab u2_abab(b,a,j,i) m2_abab(j,i,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["aabb_oooo_125"]("i,m,j,n") * m2_abab_oovv("i,j,e,f");
        }

        {

            // rl2_abab_vvoo += 1.000000 <m,j||b,f>_abab l2_abab(i,n,e,a) t2_abab(b,a,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["bbaa_vvoo_8"]("f,a,m,i") * l2_abab_oovv("i,n,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rl2_abab_vvoo += -1.000000 <j,m||a,i>_aaaa u1_aa(a,j) m2_abab(i,n,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["aa_oo_347"]("m,i") * m2_abab_oovv("i,n,e,f");
        }

        {

            // rl2_abab_vvoo += 1.000000 f_bb(a,f) l2_abab(m,n,e,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += F_blks_["bb_vv"]("a,f") * l2_abab_oovv("m,n,e,a");
        }

        if (include_u2_) {

            // rl2_abab_vvoo += 1.000000 <j,n||f,b>_bbbb u2_abab(a,b,i,j) m2_aaaa(i,m,e,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["abab_vvoo_48"]("a,f,i,n") * m2_aaaa_oovv("i,m,e,a");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // rl2_abab_vvoo += 1.000000 dp_bb(n,a) u0 u1_bb(a,i) m2_abab(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 2, o0v0: 1,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["bb_oo_313"]("n,i") * m2_abab_oovv("m,i,e,f") * u0;
        }

        if (include_u2_) {

            // rl2_abab_vvoo += 0.250000 <j,i||e,f>_abab u2_abab(b,a,j,i) m2_abab(m,n,b,a)
            // +                0.250000 <i,j||e,f>_abab u2_abab(b,a,i,j) m2_abab(m,n,b,a)
            // +                0.250000 <j,i||e,f>_abab u2_abab(a,b,j,i) m2_abab(m,n,a,b)
            // +                0.250000 <i,j||e,f>_abab u2_abab(a,b,i,j) m2_abab(m,n,a,b)
            // flops: o4v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o4v0L1: 1,
            rl2_abab_vvoo("e,f,m,n") += u2_abab_vvoo("b,a,j,i") * m2_abab_oovv("m,n,b,a") * V_blks_["abab_oovv"]("j,i,e,f");

            // rl2_abab_vvoo += 0.500000 <j,m||a,b>_aaaa u2_aaaa(a,b,i,j) m2_abab(i,n,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += 0.500000 * sharedOps["aa_oo_244"]("i,m") * m2_abab_oovv("i,n,e,f");
        }

        if (include_u1_ && include_u2_) {

            // rl2_abab_vvoo += 1.000000 <a,n||b,f>_abab u1_aa(b,i) m2_aaaa(i,m,a,e)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["abab_vvoo_167"]("a,f,i,n") * m2_aaaa_oovv("i,m,a,e");
        }

        if (include_u1_) {

            // rl2_abab_vvoo += -1.000000 dp_bb(n,f) m1_aa(m,e)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= dp_bb_ov("n,f") * m1_aa_ov("m,e");
        }

        if (include_u0_) {

            // rl2_abab_vvoo += -1.000000 dp_bb(n,f) l1_aa(m,e) u0
            // flops: o2v2L1: 2, o0v0: 1 | mem: o2v2L1: 2, o0v0: 1,
            rl2_abab_vvoo("e,f,m,n") -= dp_bb_ov("n,f") * l1_aa_ov("m,e") * u0;
        }

        {

            // rl2_abab_vvoo += -1.000000 <a,n||i,f>_abab l2_aaaa(m,i,a,e)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += V_blks_["abba_vovo"]("a,n,f,i") * l2_aaaa_oovv("m,i,a,e");
        }

        if (include_u2_) {

            // rl2_abab_vvoo += -0.500000 <j,i||e,b>_abab u2_abab(a,b,j,i) m2_abab(m,n,a,f)
            // +                -0.500000 <i,j||e,b>_abab u2_abab(a,b,i,j) m2_abab(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["aa_vv_153"]("a,e") * m2_abab_oovv("m,n,a,f");
        }

        if (include_u1_ && include_u2_) {

            // rl2_abab_vvoo += 1.000000 <n,a||f,b>_bbbb u1_bb(b,i) m2_abab(m,i,e,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["bbbb_vvoo_166"]("a,f,i,n") * m2_abab_oovv("m,i,e,a");

            // rl2_abab_vvoo += -1.000000 f_aa(m,a) u1_aa(a,i) m2_abab(i,n,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["aa_oo_367"]("i,m") * m2_abab_oovv("i,n,e,f");
        }

        {

            // rl2_abab_vvoo += -1.000000 f_bb(n,i) l2_abab(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= F_blks_["bb_oo"]("n,i") * l2_abab_oovv("m,i,e,f");

            // rl2_abab_vvoo += -1.000000 f_aa(m,i) l2_abab(i,n,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= F_blks_["aa_oo"]("m,i") * l2_abab_oovv("i,n,e,f");
        }

        if (include_u2_) {

            // rl2_abab_vvoo += 0.500000 <j,n||a,b>_bbbb u2_bbbb(a,b,i,j) m2_abab(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += 0.500000 * sharedOps["bb_oo_242"]("i,n") * m2_abab_oovv("m,i,e,f");

            // rl2_abab_vvoo += 1.000000 <m,j||e,b>_abab u2_bbbb(b,a,i,j) m2_bbbb(i,n,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["abab_vvoo_54"]("e,a,m,i") * m2_bbbb_oovv("i,n,f,a");

            // rl2_abab_vvoo += -0.500000 <j,n||a,b>_abab u2_abab(a,b,j,i) m2_abab(m,i,e,f)
            // +                -0.500000 <j,n||b,a>_abab u2_abab(b,a,j,i) m2_abab(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["bb_oo_214"]("i,n") * m2_abab_oovv("m,i,e,f");

            // rl2_abab_vvoo += 0.500000 <j,i||e,b>_aaaa u2_aaaa(b,a,j,i) m2_abab(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += 0.500000 * sharedOps["aa_vv_171"]("a,e") * m2_abab_oovv("m,n,a,f");
        }

        if (include_u1_) {

            // rl2_abab_vvoo += 1.000000 dp_aa(i,e) l2_abab(m,n,a,f) u1_aa(a,i)
            // flops: o3v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o3v1L1: 1,
            rl2_abab_vvoo("e,f,m,n") += l2_abab_oovv("m,n,a,f") * u1_aa_vo("a,i") * dp_aa_ov("i,e");
        }

        {

            // rl2_abab_vvoo += 1.000000 <j,n||e,b>_abab l2_abab(m,i,a,f) t2_abab(a,b,j,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["aabb_vvoo_9"]("e,a,n,i") * l2_abab_oovv("m,i,a,f");
        }

        if (include_u1_ && include_u2_) {

            // rl2_abab_vvoo += -1.000000 <i,a||f,b>_bbbb u1_bb(b,i) m2_abab(m,n,e,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["bb_vv_314"]("a,f") * m2_abab_oovv("m,n,e,a");
        }

        if (include_u2_) {

            // rl2_abab_vvoo += 1.000000 <j,m||e,b>_aaaa u2_abab(b,a,j,i) m2_bbbb(i,n,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["abab_vvoo_47"]("e,a,m,i") * m2_bbbb_oovv("i,n,f,a");
        }

        {

            // rl2_abab_vvoo += -0.500000 <j,n||a,b>_abab l2_abab(m,i,e,f) t2_abab(a,b,j,i)
            // +                -0.500000 <j,n||b,a>_abab l2_abab(m,i,e,f) t2_abab(b,a,j,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["bb_oo_208"]("n,i") * l2_abab_oovv("m,i,e,f");

            // rl2_abab_vvoo += -1.000000 <a,n||e,i>_abab l2_abab(m,i,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= V_blks_["abab_vovo"]("a,n,e,i") * l2_abab_oovv("m,i,a,f");

            // rl2_abab_vvoo += -0.500000 <j,i||e,b>_abab l2_abab(m,n,a,f) t2_abab(a,b,j,i)
            // +                -0.500000 <i,j||e,b>_abab l2_abab(m,n,a,f) t2_abab(a,b,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["aa_vv_150"]("a,e") * l2_abab_oovv("m,n,a,f");
        }

        if (include_u1_) {

            // rl2_abab_vvoo += 1.000000 dp_bb(n,a) l2_abab(m,i,e,f) u1_bb(a,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["bb_oo_313"]("n,i") * l2_abab_oovv("m,i,e,f");
        }

        if (include_u2_) {

            // rl2_abab_vvoo += 1.000000 <j,m||e,b>_aaaa u2_aaaa(b,a,i,j) m2_abab(i,n,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["aaaa_vvoo_51"]("a,e,i,m") * m2_abab_oovv("i,n,a,f");
        }

        if (include_u1_ && include_u2_) {

            // rl2_abab_vvoo += -1.000000 <a,n||e,b>_abab u1_bb(b,i) m2_abab(m,i,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["aabb_vvoo_205"]("a,e,i,n") * m2_abab_oovv("m,i,a,f");

            // rl2_abab_vvoo += -1.000000 <i,a||e,b>_aaaa u1_aa(b,i) m2_abab(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["aa_vv_317"]("a,e") * m2_abab_oovv("m,n,a,f");
        }

        {

            // rl2_abab_vvoo += 1.000000 <j,n||f,b>_bbbb l2_abab(m,i,e,a) t2_bbbb(b,a,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["bbbb_vvoo_2"]("a,f,i,n") * l2_abab_oovv("m,i,e,a");

            // rl2_abab_vvoo += 1.000000 f_bb(n,f) l1_aa(m,e)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += F_blks_["bb_ov"]("n,f") * l1_aa_ov("m,e");

            // rl2_abab_vvoo += 0.500000 <j,i||e,b>_aaaa l2_abab(m,n,a,f) t2_aaaa(b,a,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += 0.500000 * sharedOps["aa_vv_151"]("a,e") * l2_abab_oovv("m,n,a,f");
        }

        if (include_u2_) {

            // rl2_abab_vvoo += 1.000000 <j,n||f,b>_bbbb u2_bbbb(b,a,i,j) m2_abab(m,i,e,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["bbbb_vvoo_45"]("f,a,n,i") * m2_abab_oovv("m,i,e,a");
        }

        {

            // rl2_abab_vvoo += 1.000000 <j,m||e,b>_aaaa l2_bbbb(i,n,f,a) t2_abab(b,a,j,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["abab_vvoo_4"]("e,a,m,i") * l2_bbbb_oovv("i,n,f,a");

            // rl2_abab_vvoo += 0.500000 <m,n||i,j>_abab l2_abab(i,j,e,f)
            // +                0.500000 <m,n||j,i>_abab l2_abab(j,i,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += V_blks_["abab_oooo"]("m,n,i,j") * l2_abab_oovv("i,j,e,f");

            // rl2_abab_vvoo += -1.000000 <m,n||e,i>_abab l1_bb(i,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= V_blks_["abab_oovo"]("m,n,e,i") * l1_bb_ov("i,f");

            // rl2_abab_vvoo += -1.000000 <m,a||e,i>_abab l2_bbbb(n,i,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += V_blks_["baab_vovo"]("a,m,e,i") * l2_bbbb_oovv("n,i,a,f");
        }

        if (include_u1_) {

            // rl2_abab_vvoo += 1.000000 dp_aa(m,a) l2_abab(i,n,e,f) u1_aa(a,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["aa_oo_312"]("i,m") * l2_abab_oovv("i,n,e,f");
        }

        if (include_u1_ && include_u2_) {

            // rl2_abab_vvoo += 1.000000 <a,i||e,b>_abab u1_bb(b,i) m2_abab(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["aa_vv_316"]("a,e") * m2_abab_oovv("m,n,a,f");

            // rl2_abab_vvoo += 1.000000 <m,a||e,b>_aaaa u1_aa(b,i) m2_abab(i,n,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["aaaa_vvoo_165"]("a,e,m,i") * m2_abab_oovv("i,n,a,f");
        }

        {

            // rl2_abab_vvoo += 1.000000 <m,a||e,f>_abab l1_bb(n,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= V_blks_["baab_vovv"]("a,m,e,f") * l1_bb_ov("n,a");
        }

        if (include_u1_ && include_u2_) {

            // rl2_abab_vvoo += -1.000000 <j,n||a,i>_abab u1_aa(a,j) m2_abab(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["bb_oo_346"]("n,i") * m2_abab_oovv("m,i,e,f");
        }

        {

            // rl2_abab_vvoo += 1.000000 <a,n||e,f>_abab l1_aa(m,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += V_blks_["abab_vovv"]("a,n,e,f") * l1_aa_ov("m,a");
        }

        if (include_u1_ && include_u2_) {

            // rl2_abab_vvoo += 1.000000 <m,n||a,j>_abab u1_aa(a,i) m2_abab(i,j,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += sharedOps["aabb_oooo_306"]("m,i,n,j") * m2_abab_oovv("i,j,e,f");
        }

        {

            // rl2_abab_vvoo += 0.500000 <b,a||e,f>_abab l2_abab(m,n,b,a)
            // +                0.500000 <a,b||e,f>_abab l2_abab(m,n,a,b)
            // flops: o2v4L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += V_blks_["abab_vvvv"]("b,a,e,f") * l2_abab_oovv("m,n,b,a");
        }

        if (include_u1_ && include_u2_) {

            // rl2_abab_vvoo += 1.000000 <i,a||b,f>_abab u1_aa(b,i) m2_abab(m,n,e,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["bb_vv_315"]("a,f") * m2_abab_oovv("m,n,e,a");
        }

        if (include_u1_) {

            // rl2_abab_vvoo += -1.000000 <i,n||f,a>_bbbb u1_bb(a,i) m1_aa(m,e)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") -= sharedOps["bb_vo_326"]("f,n") * m1_aa_ov("m,e");
        }

        {

            // rl2_abab_vvoo += 0.500000 <j,n||a,b>_bbbb l2_abab(m,i,e,f) t2_bbbb(a,b,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_abab_vvoo("e,f,m,n") += 0.500000 * sharedOps["bb_oo_212"]("n,i") * l2_abab_oovv("m,i,e,f");
        }

        if (include_u2_) {

            // rm2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rm2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rm2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <j,n||e,b>_aaaa l2_abab(m,i,f,a) t2_abab(b,a,j,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = sharedOps["abab_vvoo_4"]("e,a,n,i") * l2_abab_oovv("m,i,f,a");
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <i,n||e,a>_aaaa u1_aa(a,i) m1_aa(m,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["aa_vo_327"]("e,n") * m1_aa_ov("m,f");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <j,n||e,b>_aaaa u2_aaaa(b,a,i,j) m2_aaaa(i,m,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["aaaa_vvoo_51"]("a,e,i,n") * m2_aaaa_oovv("i,m,f,a");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <n,j||e,b>_abab l2_abab(m,i,f,a) t2_bbbb(b,a,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["abab_vvoo_7"]("e,a,n,i") * l2_abab_oovv("m,i,f,a");

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <n,j||e,b>_abab l2_aaaa(i,m,f,a) t2_abab(a,b,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["aaaa_vvoo_6"]("a,e,i,n") * l2_aaaa_oovv("i,m,f,a");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) <n,a||e,b>_aaaa u1_aa(b,i) m2_aaaa(i,m,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["aaaa_vvoo_165"]("a,e,n,i") * m2_aaaa_oovv("i,m,a,f");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <n,j||e,b>_abab u2_abab(a,b,i,j) m2_aaaa(i,m,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["aaaa_vvoo_55"]("a,e,i,n") * m2_aaaa_oovv("i,m,f,a");
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) <n,i||e,a>_abab u1_bb(a,i) m1_aa(m,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= sharedOps["aa_vo_37"]("e,n") * m1_aa_ov("m,f");
        }

        {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) f_aa(n,e) l1_aa(m,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= F_blks_["aa_ov"]("n,e") * l1_aa_ov("m,f");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <n,j||e,b>_abab u2_bbbb(b,a,i,j) m2_abab(m,i,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["abab_vvoo_54"]("e,a,n,i") * m2_abab_oovv("m,i,f,a");
        }

        {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) <n,a||e,i>_abab l2_abab(m,i,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += V_blks_["baab_vovo"]("a,n,e,i") * l2_abab_oovv("m,i,f,a");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) <n,a||e,b>_abab u1_bb(b,i) m2_abab(m,i,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["abab_vvoo_163"]("e,a,n,i") * m2_abab_oovv("m,i,f,a");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <n,a||e,i>_aaaa l2_aaaa(m,i,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= V_blks_["aaaa_vovo"]("a,n,e,i") * l2_aaaa_oovv("m,i,a,f");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <j,n||e,b>_aaaa u2_abab(b,a,j,i) m2_abab(m,i,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["abab_vvoo_47"]("e,a,n,i") * m2_abab_oovv("m,i,f,a");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <j,n||e,b>_aaaa l2_aaaa(i,m,f,a) t2_aaaa(b,a,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["aaaa_vvoo_5"]("e,a,n,i") * l2_aaaa_oovv("i,m,f,a");

            // rl2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("f,e,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) <j,n||a,i>_aaaa u1_aa(a,j) m2_aaaa(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * sharedOps["aa_oo_347"]("n,i") * m2_aaaa_oovv("m,i,e,f");
        }

        {

            // tempPerm_aaaa_vvoo += -0.500000 P(m,n) <j,n||a,b>_aaaa l2_aaaa(i,m,e,f) t2_aaaa(a,b,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= 0.500000 * sharedOps["aa_oo_211"]("i,n") * l2_aaaa_oovv("i,m,e,f");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += -0.500000 P(m,n) <j,n||a,b>_aaaa u2_aaaa(a,b,i,j) m2_aaaa(i,m,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= 0.500000 * sharedOps["aa_oo_244"]("i,n") * m2_aaaa_oovv("i,m,e,f");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) <n,j||i,a>_abab u1_bb(a,j) m2_aaaa(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["aa_oo_345"]("n,i") * m2_aaaa_oovv("m,i,e,f");
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) dp_aa(n,a) l2_aaaa(i,m,e,f) u1_aa(a,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= sharedOps["aa_oo_312"]("i,n") * l2_aaaa_oovv("i,m,e,f");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) dp_aa(n,a) u0 u1_aa(a,i) m2_aaaa(i,m,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 2, o0v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") -= sharedOps["aa_oo_312"]("i,n") * u0 * m2_aaaa_oovv("i,m,e,f");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) f_aa(n,a) u1_aa(a,i) m2_aaaa(i,m,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["aa_oo_367"]("i,n") * m2_aaaa_oovv("i,m,e,f");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += 0.500000 P(m,n) <n,j||a,b>_abab u2_abab(a,b,i,j) m2_aaaa(i,m,e,f)
            // +                0.500000 P(m,n) <n,j||b,a>_abab u2_abab(b,a,i,j) m2_aaaa(i,m,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["aa_oo_213"]("n,i") * m2_aaaa_oovv("i,m,e,f");
        }

        {

            // tempPerm_aaaa_vvoo += 0.500000 P(m,n) <n,j||a,b>_abab l2_aaaa(i,m,e,f) t2_abab(a,b,i,j)
            // +                0.500000 P(m,n) <n,j||b,a>_abab l2_aaaa(i,m,e,f) t2_abab(b,a,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["aa_oo_207"]("n,i") * l2_aaaa_oovv("i,m,e,f");

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) <n,a||e,f>_aaaa l1_aa(m,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += V_blks_["aaaa_vovv"]("a,n,e,f") * l1_aa_ov("m,a");

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) f_aa(n,i) l2_aaaa(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= F_blks_["aa_oo"]("n,i") * l2_aaaa_oovv("m,i,e,f");

            // rl2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(e,f) dp_aa(i,e) l2_aaaa(m,n,f,a) u1_aa(a,i)
            // flops: o3v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o3v1L1: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * u1_aa_vo("a,i") * l2_aaaa_oovv("m,n,f,a") * dp_aa_ov("i,e");

            // tempPerm_aaaa_vvoo += -1.000000 P(e,f) <m,n||e,a>_aaaa u1_aa(a,i) m1_aa(i,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= sharedOps["aaaa_vooo_147"]("e,i,m,n") * m1_aa_ov("i,f");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(e,f) <i,a||e,b>_aaaa u1_aa(b,i) m2_aaaa(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["aa_vv_317"]("a,e") * m2_aaaa_oovv("m,n,a,f");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 P(e,f) f_aa(a,e) l2_aaaa(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += F_blks_["aa_vv"]("a,e") * l2_aaaa_oovv("m,n,a,f");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(e,f) <a,i||e,b>_abab u1_bb(b,i) m2_aaaa(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["aa_vv_316"]("a,e") * m2_aaaa_oovv("m,n,a,f");
        }

        {

            // tempPerm_aaaa_vvoo += -1.000000 P(e,f) <m,n||e,i>_aaaa l1_aa(i,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= V_blks_["aaaa_oovo"]("m,n,e,i") * l1_aa_ov("i,f");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += -0.500000 P(e,f) <j,i||e,b>_aaaa u2_aaaa(b,a,j,i) m2_aaaa(m,n,f,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= 0.500000 * sharedOps["aa_vv_171"]("a,e") * m2_aaaa_oovv("m,n,f,a");

            // tempPerm_aaaa_vvoo += 0.500000 P(e,f) <j,i||e,b>_abab u2_abab(a,b,j,i) m2_aaaa(m,n,f,a)
            // +                0.500000 P(e,f) <i,j||e,b>_abab u2_abab(a,b,i,j) m2_aaaa(m,n,f,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["aa_vv_153"]("a,e") * m2_aaaa_oovv("m,n,f,a");
        }

        {

            // tempPerm_aaaa_vvoo += -0.500000 P(e,f) <j,i||e,b>_aaaa l2_aaaa(m,n,f,a) t2_aaaa(b,a,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= 0.500000 * sharedOps["aa_vv_151"]("a,e") * l2_aaaa_oovv("m,n,f,a");

            // tempPerm_aaaa_vvoo += 0.500000 P(e,f) <j,i||e,b>_abab l2_aaaa(m,n,f,a) t2_abab(a,b,j,i)
            // +                0.500000 P(e,f) <i,j||e,b>_abab l2_aaaa(m,n,f,a) t2_abab(a,b,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += sharedOps["aa_vv_150"]("a,e") * l2_aaaa_oovv("m,n,f,a");

            // rl2_aaaa_vvoo += 1.000000 <m,n||e,f>_aaaa
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_aaaa_vvoo("e,f,m,n") += V_blks_["aaaa_oovv"]("m,n,e,f");
        }

        if (include_u2_) {

            // rl2_aaaa_vvoo += 0.250000 <j,i||e,f>_aaaa u2_aaaa(b,a,j,i) m2_aaaa(m,n,b,a)
            // flops: o4v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o4v0L1: 1,
            rl2_aaaa_vvoo("e,f,m,n") += 0.250000 * u2_aaaa_vvoo("b,a,j,i") * m2_aaaa_oovv("m,n,b,a") * V_blks_["aaaa_oovv"]("j,i,e,f");

            // rl2_aaaa_vvoo += 0.250000 <m,n||a,b>_aaaa u2_aaaa(a,b,i,j) m2_aaaa(i,j,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_aaaa_vvoo("e,f,m,n") += 0.250000 * sharedOps["aaaa_oooo_121"]("m,n,i,j") * m2_aaaa_oovv("i,j,e,f");
        }

        if (include_u1_ && include_u2_) {

            // rl2_aaaa_vvoo += 1.000000 <m,n||a,j>_aaaa u1_aa(a,i) m2_aaaa(i,j,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_aaaa_vvoo("e,f,m,n") += sharedOps["aaaa_oooo_308"]("m,n,j,i") * m2_aaaa_oovv("i,j,e,f");
        }

        {

            // rl2_aaaa_vvoo += 0.250000 <m,n||a,b>_aaaa l2_aaaa(i,j,e,f) t2_aaaa(a,b,i,j)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_aaaa_vvoo("e,f,m,n") += 0.250000 * sharedOps["aaaa_oooo_108"]("i,j,m,n") * l2_aaaa_oovv("i,j,e,f");

            // rl2_aaaa_vvoo += 0.250000 <j,i||e,f>_aaaa l2_aaaa(m,n,b,a) t2_aaaa(b,a,j,i)
            // flops: o4v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o4v0L1: 1,
            rl2_aaaa_vvoo("e,f,m,n") += 0.250000 * l2_aaaa_oovv("m,n,b,a") * t2_aaaa_vvoo("b,a,j,i") * V_blks_["aaaa_oovv"]("j,i,e,f");

            // rl2_aaaa_vvoo += 0.500000 <b,a||e,f>_aaaa l2_aaaa(m,n,b,a)
            // flops: o2v4L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_aaaa_vvoo("e,f,m,n") += 0.500000 * V_blks_["aaaa_vvvv"]("b,a,e,f") * l2_aaaa_oovv("m,n,b,a");

            // rl2_aaaa_vvoo += 0.500000 <m,n||i,j>_aaaa l2_aaaa(i,j,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            rl2_aaaa_vvoo("e,f,m,n") += 0.500000 * V_blks_["aaaa_oooo"]("m,n,i,j") * l2_aaaa_oovv("i,j,e,f");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += 0.250000 <k,m||b,c>_abab t2_abab(b,c,i,j) u1_aa(a,k) m2_abab(i,j,a,e)
            // +            0.250000 <k,m||b,c>_abab t2_abab(b,c,j,i) u1_aa(a,k) m2_abab(j,i,a,e)
            // +            0.250000 <k,m||c,b>_abab t2_abab(c,b,i,j) u1_aa(a,k) m2_abab(i,j,a,e)
            // +            0.250000 <k,m||c,b>_abab t2_abab(c,b,j,i) u1_aa(a,k) m2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["aabb_vooo_478"]("a,i,m,j") * m2_abab_oovv("i,j,a,e");
        }

        if (include_u2_) {

            // rl1_bb_vo += -0.500000 f_bb(m,b) u2_bbbb(b,a,i,j) m2_bbbb(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= 0.500000 * sharedOps["bbbb_vooo_302"]("a,i,j,m") * m2_bbbb_oovv("i,j,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += 0.500000 <j,m||a,b>_bbbb u2_bbbb(a,b,i,j) m1_bb(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += 0.500000 * sharedOps["bb_oo_242"]("i,m") * m1_bb_ov("i,e");
        }

        if (include_u2_) {

            // rl1_bb_vo += -1.000000 <k,m||b,j>_bbbb u2_bbbb(b,a,i,k) m2_bbbb(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["bbbb_vooo_134"]("a,m,j,i") * m2_bbbb_oovv("i,j,e,a");
        }

        {

            // rl1_bb_vo += 1.000000 <k,m||b,j>_abab l2_abab(i,j,a,e) t2_aaaa(b,a,i,k)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["aabb_vooo_120"]("a,i,m,j") * l2_abab_oovv("i,j,a,e");
        }

        if (include_u0_) {

            // rl1_bb_vo += -1.000000 dp_bb(m,e) u0
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rl1_bb_vo("e,m") -= dp_bb_ov("m,e") * u0;
        }

        if (include_u2_) {

            // rl1_bb_vo += 0.500000 dp_bb(i,e) l2_bbbb(j,m,a,b) u2_bbbb(a,b,j,i)
            // flops: o3v2L1: 1, o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2, o2v0L1: 1,
            rl1_bb_vo("e,m") += 0.500000 * l2_bbbb_oovv("j,m,a,b") * u2_bbbb_vvoo("a,b,j,i") * dp_bb_ov("i,e");

            // rl1_bb_vo += -1.000000 <j,b||c,e>_abab u2_abab(c,a,j,i) m2_bbbb(i,m,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["bbbb_vvvo_35"]("a,b,e,i") * m2_bbbb_oovv("i,m,b,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += -0.500000 <j,i||b,e>_abab u2_abab(b,a,j,i) m1_bb(m,a)
            // +            -0.500000 <i,j||b,e>_abab u2_abab(b,a,i,j) m1_bb(m,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["bb_vv_154"]("a,e") * m1_bb_ov("m,a");
        }

        if (include_u0_ && include_u1_) {

            // rl1_bb_vo += -1.000000 <i,m||e,a>_bbbb u1_bb(a,i) m0
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["bb_vo_326"]("e,m") * m0;
        }

        {

            // rl1_bb_vo += 0.500000 <j,m||a,b>_bbbb l1_bb(i,e) t2_bbbb(a,b,i,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += 0.500000 * sharedOps["bb_oo_212"]("m,i") * l1_bb_ov("i,e");

            // rl1_bb_vo += 0.500000 <b,a||e,i>_bbbb l2_bbbb(m,i,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += 0.500000 * V_blks_["bbbb_vvvo"]("b,a,e,i") * l2_bbbb_oovv("m,i,b,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += -0.250000 <k,m||b,c>_bbbb t2_bbbb(b,c,i,j) u1_bb(a,k) m2_bbbb(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= 0.250000 * sharedOps["bbbb_vooo_477"]("a,i,j,m") * m2_bbbb_oovv("i,j,e,a");
        }

        {

            // rl1_bb_vo += 0.500000 <j,i||e,b>_bbbb l1_bb(m,a) t2_bbbb(b,a,j,i)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += 0.500000 * sharedOps["bb_vv_152"]("e,a") * l1_bb_ov("m,a");

            // rl1_bb_vo += -0.500000 <a,m||i,j>_abab l2_abab(i,j,a,e)
            // +            -0.500000 <a,m||j,i>_abab l2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= V_blks_["abab_vooo"]("a,m,i,j") * l2_abab_oovv("i,j,a,e");
        }

        if (include_u0_ && include_u1_) {

            // rl1_bb_vo += 1.000000 <i,m||a,e>_abab u1_aa(a,i) m0
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["bb_vo_36"]("e,m") * m0;
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += 1.000000 <k,m||c,b>_abab t2_abab(c,a,k,i) u1_bb(b,j) m2_bbbb(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["bbbb_vooo_471"]("a,j,m,i") * m2_bbbb_oovv("i,j,e,a");
        }

        {

            // rl1_bb_vo += 0.500000 <m,a||i,j>_bbbb l2_bbbb(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= 0.500000 * V_blks_["bbbb_vooo"]("a,m,i,j") * l2_bbbb_oovv("i,j,a,e");
        }

        if (include_u2_) {

            // rl1_bb_vo += -0.500000 f_bb(m,b) u2_abab(a,b,i,j) m2_abab(i,j,a,e)
            // +            -0.500000 f_bb(m,b) u2_abab(a,b,j,i) m2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["aabb_vooo_298"]("a,i,m,j") * m2_abab_oovv("i,j,a,e");
        }

        {

            // rl1_bb_vo += 0.500000 <b,a||i,e>_abab l2_abab(i,m,b,a)
            // +            0.500000 <a,b||i,e>_abab l2_abab(i,m,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= V_blks_["abba_vvvo"]("b,a,e,i") * l2_abab_oovv("i,m,b,a");
        }

        if (include_u2_) {

            // rl1_bb_vo += 1.000000 <k,m||j,b>_abab u2_abab(a,b,k,i) m2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["aabb_vooo_140"]("a,j,m,i") * m2_abab_oovv("j,i,a,e");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += 0.500000 <k,m||i,j>_abab u1_aa(a,k) m2_abab(i,j,a,e)
            // +            0.500000 <k,m||j,i>_abab u1_aa(a,k) m2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["aabb_vooo_311"]("a,i,m,j") * m2_abab_oovv("i,j,a,e");
        }

        if (include_u0_ && include_u2_) {

            // rl1_bb_vo += 0.500000 dp_bb(m,b) u0 u2_bbbb(b,a,i,j) m2_bbbb(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_bb_vo("e,m") += 0.500000 * sharedOps["bbbb_vooo_232"]("a,m,i,j") * m2_bbbb_oovv("i,j,e,a") * u0;
        }

        {

            // rl1_bb_vo += -0.500000 <j,m||a,b>_abab l1_bb(i,e) t2_abab(a,b,j,i)
            // +            -0.500000 <j,m||b,a>_abab l1_bb(i,e) t2_abab(b,a,j,i)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["bb_oo_208"]("m,i") * l1_bb_ov("i,e");
        }

        if (include_u2_) {

            // rl1_bb_vo += -1.000000 <j,b||e,c>_bbbb u2_abab(a,c,i,j) m2_abab(i,m,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["abba_vvvo_25"]("a,b,e,i") * m2_abab_oovv("i,m,a,b");
        }

        {

            // rl1_bb_vo += -1.000000 <j,b||e,c>_bbbb l2_abab(i,m,a,b) t2_abab(a,c,i,j)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["abba_vvvo_21"]("a,b,e,i") * l2_abab_oovv("i,m,a,b");

            // rl1_bb_vo += -1.000000 <j,b||c,e>_abab l2_bbbb(i,m,b,a) t2_abab(c,a,j,i)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["bbbb_vvvo_14"]("b,e,a,i") * l2_bbbb_oovv("i,m,b,a");
        }

        if (include_u0_) {

            // rl1_bb_vo += 0.500000 dp_bb(m,b) l2_bbbb(i,j,e,a) t2_bbbb(b,a,i,j) u0
            // flops: o3v2L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_bb_vo("e,m") += 0.500000 * sharedOps["bbbb_vooo_216"]("a,i,j,m") * u0 * l2_bbbb_oovv("i,j,e,a");
        }

        if (include_u1_) {

            // rl1_bb_vo += 1.000000 <a,m||b,e>_abab u1_aa(b,i) m1_aa(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["abab_vvoo_167"]("a,e,i,m") * m1_aa_ov("i,a");
        }

        if (include_u2_) {

            // rl1_bb_vo += -1.000000 <b,j||c,e>_abab u2_abab(c,a,i,j) m2_abab(i,m,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["abba_vvvo_34"]("b,a,e,i") * m2_abab_oovv("i,m,b,a");
        }

        {

            // rl1_bb_vo += 1.000000 <a,m||i,e>_abab l1_aa(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= V_blks_["abba_vovo"]("a,m,e,i") * l1_aa_ov("i,a");

            // rl1_bb_vo += -1.000000 <b,j||c,e>_abab l2_abab(i,m,b,a) t2_abab(c,a,i,j)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["abba_vvvo_15"]("b,e,a,i") * l2_abab_oovv("i,m,b,a");

            // rl1_bb_vo += -0.500000 f_bb(m,b) l2_bbbb(i,j,e,a) t2_bbbb(b,a,i,j)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= 0.500000 * sharedOps["bbbb_vooo_243"]("a,m,i,j") * l2_bbbb_oovv("i,j,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += 1.000000 <k,m||b,c>_bbbb t2_abab(a,c,i,k) u1_bb(b,j) m2_abab(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["aabb_vooo_461"]("a,i,j,m") * m2_abab_oovv("i,j,a,e");
        }

        if (include_u2_) {

            // rl1_bb_vo += -0.250000 <a,m||b,c>_abab u2_abab(b,c,i,j) m2_abab(i,j,a,e)
            // +            -0.250000 <a,m||b,c>_abab u2_abab(b,c,j,i) m2_abab(j,i,a,e)
            // +            -0.250000 <a,m||c,b>_abab u2_abab(c,b,i,j) m2_abab(i,j,a,e)
            // +            -0.250000 <a,m||c,b>_abab u2_abab(c,b,j,i) m2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["aabb_vooo_81"]("a,i,j,m") * m2_abab_oovv("i,j,a,e");
        }

        {

            // rl1_bb_vo += -1.000000 <k,m||b,j>_bbbb l2_abab(i,j,a,e) t2_abab(a,b,i,k)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["aabb_vooo_114"]("a,i,m,j") * l2_abab_oovv("i,j,a,e");
        }

        if (include_u2_) {

            // rl1_bb_vo += 0.500000 dp_bb(m,b) t2_bbbb(b,a,i,j) m2_bbbb(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += 0.500000 * sharedOps["bbbb_vooo_216"]("a,i,j,m") * m2_bbbb_oovv("i,j,e,a");
        }

        if (include_u1_) {

            // rl1_bb_vo += 1.000000 <m,a||e,b>_bbbb u1_bb(b,i) m1_bb(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["bbbb_vvoo_166"]("a,e,i,m") * m1_bb_ov("i,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += -1.000000 <a,m||j,b>_abab u1_bb(b,i) m2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["aabb_vooo_295"]("a,j,m,i") * m2_abab_oovv("j,i,a,e");

            // rl1_bb_vo += -0.500000 <k,m||b,c>_abab t2_abab(a,c,i,j) u1_aa(b,k) m2_abab(i,j,a,e)
            // +            -0.500000 <k,m||b,c>_abab t2_abab(a,c,j,i) u1_aa(b,k) m2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["aabb_vooo_474"]("a,i,m,j") * m2_abab_oovv("i,j,a,e");
        }

        {

            // rl1_bb_vo += 0.250000 <k,j||i,e>_abab l2_abab(i,m,b,a) t2_abab(b,a,k,j)
            // +            0.250000 <j,k||i,e>_abab l2_abab(i,m,b,a) t2_abab(b,a,j,k)
            // +            0.250000 <k,j||i,e>_abab l2_abab(i,m,a,b) t2_abab(a,b,k,j)
            // +            0.250000 <j,k||i,e>_abab l2_abab(i,m,a,b) t2_abab(a,b,j,k)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["abba_vvvo_46"]("b,e,a,i") * l2_abab_oovv("i,m,b,a");
        }

        if (include_u0_ && include_u1_) {

            // rl1_bb_vo += 1.000000 dp_bb(m,a) u0 u1_bb(a,i) m1_bb(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_bb_vo("e,m") += sharedOps["bb_oo_313"]("m,i") * m1_bb_ov("i,e") * u0;
        }

        {

            // rl1_bb_vo += -0.500000 <j,i||b,e>_abab l1_bb(m,a) t2_abab(b,a,j,i)
            // +            -0.500000 <i,j||b,e>_abab l1_bb(m,a) t2_abab(b,a,i,j)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["bb_vv_149"]("e,a") * l1_bb_ov("m,a");
        }

        if (include_u2_) {

            // rl1_bb_vo += -1.000000 <j,b||c,e>_abab u2_aaaa(c,a,i,j) m2_abab(i,m,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["abba_vvvo_30"]("a,b,e,i") * m2_abab_oovv("i,m,a,b");

            // rl1_bb_vo += -1.000000 <k,m||b,j>_bbbb u2_abab(a,b,i,k) m2_abab(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["aabb_vooo_136"]("a,i,m,j") * m2_abab_oovv("i,j,a,e");

            // rl1_bb_vo += 0.250000 <k,j||i,e>_abab u2_abab(b,a,k,j) m2_abab(i,m,b,a)
            // +            0.250000 <j,k||i,e>_abab u2_abab(b,a,j,k) m2_abab(i,m,b,a)
            // +            0.250000 <k,j||i,e>_abab u2_abab(a,b,k,j) m2_abab(i,m,a,b)
            // +            0.250000 <j,k||i,e>_abab u2_abab(a,b,j,k) m2_abab(i,m,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["abba_vvvo_71"]("b,a,e,i") * m2_abab_oovv("i,m,b,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += 1.000000 <k,m||b,c>_abab t2_abab(a,c,k,i) u1_aa(b,j) m2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["aabb_vooo_473"]("a,j,m,i") * m2_abab_oovv("j,i,a,e");
        }

        if (include_u0_) {

            // rl1_bb_vo += -1.000000 dp_bb(m,e) m0
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= dp_bb_ov("m,e") * m0;
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += 1.000000 <m,a||b,j>_bbbb u1_bb(b,i) m2_bbbb(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["bbbb_vooo_301"]("a,i,m,j") * m2_bbbb_oovv("i,j,a,e");
        }

        if (include_u2_) {

            // rl1_bb_vo += 0.250000 <m,a||b,c>_bbbb u2_bbbb(b,c,i,j) m2_bbbb(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= 0.250000 * sharedOps["bbbb_vooo_100"]("a,m,i,j") * m2_bbbb_oovv("i,j,a,e");

            // rl1_bb_vo += 0.250000 <k,j||e,i>_bbbb u2_bbbb(b,a,k,j) m2_bbbb(m,i,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += 0.250000 * sharedOps["bbbb_vvvo_92"]("e,b,a,i") * m2_bbbb_oovv("m,i,b,a");
        }

        {

            // rl1_bb_vo += -0.500000 f_bb(m,b) l2_abab(i,j,a,e) t2_abab(a,b,i,j)
            // +            -0.500000 f_bb(m,b) l2_abab(j,i,a,e) t2_abab(a,b,j,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["aabb_vooo_262"]("a,i,j,m") * l2_abab_oovv("i,j,a,e");
        }

        if (include_u1_) {

            // rl1_bb_vo += -1.000000 f_bb(m,a) u1_bb(a,i) m1_bb(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["bb_oo_366"]("i,m") * m1_bb_ov("i,e");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += -0.500000 <k,m||b,c>_abab t2_bbbb(c,a,i,j) u1_aa(b,k) m2_bbbb(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= 0.500000 * sharedOps["bbbb_vooo_472"]("a,m,i,j") * m2_bbbb_oovv("i,j,e,a");

            // rl1_bb_vo += 1.000000 <k,m||c,b>_abab t2_aaaa(c,a,i,k) u1_bb(b,j) m2_abab(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["aabb_vooo_303"]("a,i,j,m") * m2_abab_oovv("i,j,a,e");
        }

        if (include_u2_) {

            // rl1_bb_vo += 0.500000 dp_bb(m,a) l2_abab(j,i,b,e) u2_abab(b,a,j,i)
            // +            0.500000 dp_bb(m,a) l2_abab(i,j,b,e) u2_abab(b,a,i,j)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["aabb_vooo_224"]("b,j,i,m") * l2_abab_oovv("j,i,b,e");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += -0.250000 <k,j||e,c>_bbbb t2_bbbb(b,a,k,j) u1_bb(c,i) m2_bbbb(i,m,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= 0.250000 * sharedOps["bbbb_vvvo_69"]("e,b,a,i") * m2_bbbb_oovv("i,m,b,a");
        }

        if (include_u1_) {

            // rl1_bb_vo += -1.000000 <j,m||a,i>_abab u1_aa(a,j) m1_bb(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["bb_oo_346"]("m,i") * m1_bb_ov("i,e");
        }

        if (include_u0_ && include_u2_) {

            // rl1_bb_vo += 0.500000 dp_bb(m,b) u0 u2_abab(a,b,i,j) m2_abab(i,j,a,e)
            // +            0.500000 dp_bb(m,b) u0 u2_abab(a,b,j,i) m2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_bb_vo("e,m") += sharedOps["aabb_vooo_224"]("a,i,j,m") * m2_abab_oovv("i,j,a,e") * u0;
        }

        {

            // rl1_bb_vo += 1.000000 f_bb(m,e)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rl1_bb_vo("e,m") += F_blks_["bb_ov"]("m,e");

            // rl1_bb_vo += 0.250000 <k,j||e,i>_bbbb l2_bbbb(m,i,b,a) t2_bbbb(b,a,k,j)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += 0.250000 * sharedOps["bbbb_vvvo_44"]("b,a,e,i") * l2_bbbb_oovv("m,i,b,a");

            // rl1_bb_vo += 1.000000 <k,m||b,j>_abab l2_bbbb(i,j,e,a) t2_abab(b,a,k,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["bbbb_vooo_116"]("a,i,m,j") * l2_bbbb_oovv("i,j,e,a");
        }

        if (include_u1_) {

            // rl1_bb_vo += -1.000000 <j,m||a,i>_bbbb u1_bb(a,j) m1_bb(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["bb_oo_344"]("m,i") * m1_bb_ov("i,e");
        }

        if (include_u2_) {

            // rl1_bb_vo += 1.000000 <k,m||b,j>_abab u2_abab(b,a,k,i) m2_bbbb(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["bbbb_vooo_133"]("a,m,j,i") * m2_bbbb_oovv("i,j,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += 0.250000 <k,j||c,e>_abab t2_abab(b,a,k,j) u1_aa(c,i) m2_abab(i,m,b,a)
            // +            0.250000 <j,k||c,e>_abab t2_abab(b,a,j,k) u1_aa(c,i) m2_abab(i,m,b,a)
            // +            0.250000 <k,j||c,e>_abab t2_abab(a,b,k,j) u1_aa(c,i) m2_abab(i,m,a,b)
            // +            0.250000 <j,k||c,e>_abab t2_abab(a,b,j,k) u1_aa(c,i) m2_abab(i,m,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["abba_vvvo_405"]("b,e,a,i") * m2_abab_oovv("i,m,b,a");

            // rl1_bb_vo += 0.500000 <b,a||c,e>_abab u1_aa(c,i) m2_abab(i,m,b,a)
            // +            0.500000 <a,b||c,e>_abab u1_aa(c,i) m2_abab(i,m,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["abba_vvvo_142"]("b,a,e,i") * m2_abab_oovv("i,m,b,a");
        }

        if (include_u1_) {

            // rl1_bb_vo += 1.000000 <i,a||b,e>_abab u1_aa(b,i) m1_bb(m,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["bb_vv_315"]("a,e") * m1_bb_ov("m,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += -0.500000 <b,a||e,c>_bbbb u1_bb(c,i) m2_bbbb(i,m,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= 0.500000 * sharedOps["bbbb_vvvo_144"]("b,a,e,i") * m2_bbbb_oovv("i,m,b,a");
        }

        {

            // rl1_bb_vo += 0.250000 <m,a||b,c>_bbbb l2_bbbb(i,j,a,e) t2_bbbb(b,c,i,j)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= 0.250000 * sharedOps["bbbb_vooo_42"]("a,i,j,m") * l2_bbbb_oovv("i,j,a,e");
        }

        if (include_u1_) {

            // rl1_bb_vo += -1.000000 <i,a||e,b>_bbbb u1_bb(b,i) m1_bb(m,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["bb_vv_314"]("a,e") * m1_bb_ov("m,a");
        }

        {

            // rl1_bb_vo += -1.000000 <j,b||c,e>_abab l2_abab(i,m,a,b) t2_aaaa(c,a,i,j)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["abba_vvvo_16"]("a,b,e,i") * l2_abab_oovv("i,m,a,b");
        }

        if (include_u2_) {

            // rl1_bb_vo += 1.000000 <k,m||b,j>_abab u2_aaaa(b,a,i,k) m2_abab(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["aabb_vooo_130"]("a,i,m,j") * m2_abab_oovv("i,j,a,e");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += 1.000000 <k,m||b,c>_bbbb t2_bbbb(c,a,i,k) u1_bb(b,j) m2_bbbb(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["bbbb_vooo_450"]("a,i,m,j") * m2_bbbb_oovv("i,j,e,a");

            // rl1_bb_vo += -1.000000 <a,m||b,j>_abab u1_aa(b,i) m2_abab(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["aabb_vooo_294"]("a,i,m,j") * m2_abab_oovv("i,j,a,e");
        }

        if (include_u2_) {

            // rl1_bb_vo += 0.500000 dp_bb(m,b) t2_abab(a,b,i,j) m2_abab(i,j,a,e)
            // +            0.500000 dp_bb(m,b) t2_abab(a,b,j,i) m2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["aabb_vooo_228"]("a,i,j,m") * m2_abab_oovv("i,j,a,e");
        }

        if (include_u0_) {

            // rl1_bb_vo += 0.500000 dp_bb(m,b) l2_abab(i,j,a,e) t2_abab(a,b,i,j) u0
            // +            0.500000 dp_bb(m,b) l2_abab(j,i,a,e) t2_abab(a,b,j,i) u0
            // flops: o3v2L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_bb_vo("e,m") += sharedOps["aabb_vooo_228"]("a,i,j,m") * u0 * l2_abab_oovv("i,j,a,e");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += -0.500000 <k,m||b,c>_bbbb t2_abab(a,c,i,j) u1_bb(b,k) m2_abab(i,j,a,e)
            // +            -0.500000 <k,m||b,c>_bbbb t2_abab(a,c,j,i) u1_bb(b,k) m2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["aabb_vooo_272"]("a,i,m,j") * m2_abab_oovv("i,j,a,e");
        }

        if (include_u1_) {

            // rl1_bb_vo += 1.000000 dp_bb(i,e) l1_bb(m,a) u1_bb(a,i)
            // flops: o2v1L1: 2, o1v1L1: 1 | mem: o1v1L1: 2, o2v0L1: 1,
            rl1_bb_vo("e,m") += l1_bb_ov("m,a") * u1_bb_vo("a,i") * dp_bb_ov("i,e");
        }

        if (include_u2_) {

            // rl1_bb_vo += 0.500000 dp_bb(i,e) l2_abab(j,m,a,b) u2_abab(a,b,j,i)
            // +            0.500000 dp_bb(i,e) l2_abab(j,m,b,a) u2_abab(b,a,j,i)
            // flops: o3v2L1: 1, o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2, o2v0L1: 1,
            rl1_bb_vo("e,m") += l2_abab_oovv("j,m,a,b") * u2_abab_vvoo("a,b,j,i") * dp_bb_ov("i,e");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += -0.500000 <j,m||a,b>_abab u2_abab(a,b,j,i) m1_bb(i,e)
            // +            -0.500000 <j,m||b,a>_abab u2_abab(b,a,j,i) m1_bb(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["bb_oo_214"]("i,m") * m1_bb_ov("i,e");
        }

        if (include_u1_) {

            // rl1_bb_vo += 1.000000 dp_bb(m,a) l1_bb(i,e) u1_bb(a,i)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["bb_oo_313"]("m,i") * l1_bb_ov("i,e");
        }

        {

            // rl1_bb_vo += -1.000000 <k,m||b,j>_bbbb l2_bbbb(i,j,e,a) t2_bbbb(b,a,i,k)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["bbbb_vooo_112"]("a,i,m,j") * l2_bbbb_oovv("i,j,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += -0.500000 <k,m||b,c>_bbbb t2_bbbb(c,a,i,j) u1_bb(b,k) m2_bbbb(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= 0.500000 * sharedOps["bbbb_vooo_274"]("a,m,i,j") * m2_bbbb_oovv("i,j,e,a");

            // rl1_bb_vo += 0.500000 <j,i||e,b>_bbbb u2_bbbb(b,a,j,i) m1_bb(m,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += 0.500000 * sharedOps["bb_vv_169"]("a,e") * m1_bb_ov("m,a");
        }

        {

            // rl1_bb_vo += 1.000000 f_bb(a,e) l1_bb(m,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += F_blks_["bb_vv"]("a,e") * l1_bb_ov("m,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_bb_vo += -0.500000 <k,m||i,j>_bbbb u1_bb(a,k) m2_bbbb(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= 0.500000 * sharedOps["bbbb_vooo_310"]("a,m,i,j") * m2_bbbb_oovv("i,j,e,a");
        }

        {

            // rl1_bb_vo += -1.000000 <j,b||e,c>_bbbb l2_bbbb(i,m,b,a) t2_bbbb(c,a,i,j)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["bbbb_vvvo_20"]("a,b,e,i") * l2_bbbb_oovv("i,m,b,a");

            // rl1_bb_vo += 1.000000 <k,m||j,b>_abab l2_abab(j,i,a,e) t2_abab(a,b,k,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["aabb_vooo_113"]("a,j,m,i") * l2_abab_oovv("j,i,a,e");

            // rl1_bb_vo += 1.000000 <m,a||e,i>_bbbb l1_bb(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= V_blks_["bbbb_vovo"]("a,m,e,i") * l1_bb_ov("i,a");
        }

        if (include_u2_) {

            // rl1_bb_vo += 0.500000 dp_bb(m,a) l2_bbbb(j,i,e,b) u2_bbbb(a,b,j,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += 0.500000 * sharedOps["bbbb_vooo_232"]("b,m,j,i") * l2_bbbb_oovv("j,i,e,b");
        }

        {

            // rl1_bb_vo += -1.000000 f_bb(m,i) l1_bb(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= F_blks_["bb_oo"]("m,i") * l1_bb_ov("i,e");
        }

        if (include_u2_) {

            // rl1_bb_vo += -1.000000 <j,b||e,c>_bbbb u2_bbbb(c,a,i,j) m2_bbbb(i,m,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") += sharedOps["bbbb_vvvo_33"]("b,e,a,i") * m2_bbbb_oovv("i,m,b,a");
        }

        {

            // rl1_bb_vo += -0.250000 <a,m||b,c>_abab l2_abab(i,j,a,e) t2_abab(b,c,i,j)
            // +            -0.250000 <a,m||c,b>_abab l2_abab(i,j,a,e) t2_abab(c,b,i,j)
            // +            -0.250000 <a,m||b,c>_abab l2_abab(j,i,a,e) t2_abab(b,c,j,i)
            // +            -0.250000 <a,m||c,b>_abab l2_abab(j,i,a,e) t2_abab(c,b,j,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_bb_vo("e,m") -= sharedOps["aabb_vooo_41"]("a,i,m,j") * l2_abab_oovv("i,j,a,e");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += -0.500000 <m,k||c,b>_abab t2_aaaa(c,a,i,j) u1_bb(b,k) m2_aaaa(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= 0.500000 * sharedOps["aaaa_vooo_470"]("a,i,j,m") * m2_aaaa_oovv("i,j,e,a");
        }

        if (include_u2_) {

            // rl1_aa_vo += -0.250000 <m,a||b,c>_abab u2_abab(b,c,i,j) m2_abab(i,j,e,a)
            // +            -0.250000 <m,a||b,c>_abab u2_abab(b,c,j,i) m2_abab(j,i,e,a)
            // +            -0.250000 <m,a||c,b>_abab u2_abab(c,b,i,j) m2_abab(i,j,e,a)
            // +            -0.250000 <m,a||c,b>_abab u2_abab(c,b,j,i) m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["baab_vooo_99"]("a,i,m,j") * m2_abab_oovv("i,j,e,a");
        }

        {

            // rl1_aa_vo += 0.500000 <j,m||a,b>_aaaa l1_aa(i,e) t2_aaaa(a,b,i,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += 0.500000 * sharedOps["aa_oo_211"]("i,m") * l1_aa_ov("i,e");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += 1.000000 <k,m||b,c>_aaaa t2_aaaa(c,a,i,k) u1_aa(b,j) m2_aaaa(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["aaaa_vooo_465"]("a,m,i,j") * m2_aaaa_oovv("i,j,e,a");
        }

        if (include_u1_) {

            // rl1_aa_vo += -1.000000 <i,a||e,b>_aaaa u1_aa(b,i) m1_aa(m,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["aa_vv_317"]("a,e") * m1_aa_ov("m,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += 0.500000 <b,a||e,c>_abab u1_bb(c,i) m2_abab(m,i,b,a)
            // +            0.500000 <a,b||e,c>_abab u1_bb(c,i) m2_abab(m,i,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["aabb_vvvo_141"]("b,e,a,i") * m2_abab_oovv("m,i,b,a");
        }

        {

            // rl1_aa_vo += 1.000000 <m,k||j,b>_abab l2_aaaa(i,j,e,a) t2_abab(a,b,i,k)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["aaaa_vooo_115"]("a,m,j,i") * l2_aaaa_oovv("i,j,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += 1.000000 <m,a||b,j>_aaaa u1_aa(b,i) m2_aaaa(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["aaaa_vooo_299"]("a,m,j,i") * m2_aaaa_oovv("i,j,a,e");

            // rl1_aa_vo += -0.500000 <k,m||b,c>_aaaa t2_aaaa(c,a,i,j) u1_aa(b,k) m2_aaaa(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= 0.500000 * sharedOps["aaaa_vooo_273"]("a,m,i,j") * m2_aaaa_oovv("i,j,e,a");
        }

        if (include_u1_) {

            // rl1_aa_vo += -1.000000 <j,m||a,i>_aaaa u1_aa(a,j) m1_aa(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["aa_oo_347"]("m,i") * m1_aa_ov("i,e");
        }

        {

            // rl1_aa_vo += -0.500000 <m,a||i,j>_abab l2_abab(i,j,e,a)
            // +            -0.500000 <m,a||j,i>_abab l2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += V_blks_["baab_vooo"]("a,m,i,j") * l2_abab_oovv("i,j,e,a");

            // rl1_aa_vo += 1.000000 <m,a||e,i>_aaaa l1_aa(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= V_blks_["aaaa_vovo"]("a,m,e,i") * l1_aa_ov("i,a");
        }

        if (include_u2_) {

            // rl1_aa_vo += -0.500000 f_aa(m,b) u2_abab(b,a,i,j) m2_abab(i,j,e,a)
            // +            -0.500000 f_aa(m,b) u2_abab(b,a,j,i) m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["baab_vooo_300"]("a,i,m,j") * m2_abab_oovv("i,j,e,a");
        }

        {

            // rl1_aa_vo += -0.250000 <m,a||b,c>_abab l2_abab(i,j,e,a) t2_abab(b,c,i,j)
            // +            -0.250000 <m,a||c,b>_abab l2_abab(i,j,e,a) t2_abab(c,b,i,j)
            // +            -0.250000 <m,a||b,c>_abab l2_abab(j,i,e,a) t2_abab(b,c,j,i)
            // +            -0.250000 <m,a||c,b>_abab l2_abab(j,i,e,a) t2_abab(c,b,j,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["baab_vooo_40"]("a,m,i,j") * l2_abab_oovv("i,j,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += 0.250000 <m,k||b,c>_abab t2_abab(b,c,i,j) u1_bb(a,k) m2_abab(i,j,e,a)
            // +            0.250000 <m,k||b,c>_abab t2_abab(b,c,j,i) u1_bb(a,k) m2_abab(j,i,e,a)
            // +            0.250000 <m,k||c,b>_abab t2_abab(c,b,i,j) u1_bb(a,k) m2_abab(i,j,e,a)
            // +            0.250000 <m,k||c,b>_abab t2_abab(c,b,j,i) u1_bb(a,k) m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["baab_vooo_476"]("a,m,i,j") * m2_abab_oovv("i,j,e,a");

            // rl1_aa_vo += 0.250000 <k,j||e,c>_abab t2_abab(b,a,k,j) u1_bb(c,i) m2_abab(m,i,b,a)
            // +            0.250000 <j,k||e,c>_abab t2_abab(b,a,j,k) u1_bb(c,i) m2_abab(m,i,b,a)
            // +            0.250000 <k,j||e,c>_abab t2_abab(a,b,k,j) u1_bb(c,i) m2_abab(m,i,a,b)
            // +            0.250000 <j,k||e,c>_abab t2_abab(a,b,j,k) u1_bb(c,i) m2_abab(m,i,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["aabb_vvvo_404"]("e,b,a,i") * m2_abab_oovv("m,i,b,a");
        }

        {

            // rl1_aa_vo += -1.000000 <k,m||b,j>_aaaa l2_abab(j,i,e,a) t2_abab(b,a,k,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["baab_vooo_117"]("a,m,j,i") * l2_abab_oovv("j,i,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += 1.000000 <m,k||b,c>_abab t2_abab(a,c,i,k) u1_aa(b,j) m2_aaaa(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["aaaa_vooo_453"]("a,i,m,j") * m2_aaaa_oovv("i,j,e,a");

            // rl1_aa_vo += -0.500000 <b,a||e,c>_aaaa u1_aa(c,i) m2_aaaa(i,m,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= 0.500000 * sharedOps["aaaa_vvvo_143"]("b,a,e,i") * m2_aaaa_oovv("i,m,b,a");
        }

        {

            // rl1_aa_vo += 0.500000 <m,a||i,j>_aaaa l2_aaaa(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= 0.500000 * V_blks_["aaaa_vooo"]("a,m,i,j") * l2_aaaa_oovv("i,j,a,e");
        }

        if (include_u2_) {

            // rl1_aa_vo += 1.000000 <m,k||b,j>_abab u2_abab(b,a,i,k) m2_abab(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["baab_vooo_132"]("a,i,m,j") * m2_abab_oovv("i,j,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += -0.500000 <k,m||i,j>_aaaa u1_aa(a,k) m2_aaaa(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= 0.500000 * sharedOps["aaaa_vooo_304"]("a,m,i,j") * m2_aaaa_oovv("i,j,e,a");
        }

        if (include_u2_) {

            // rl1_aa_vo += 0.500000 dp_aa(m,a) l2_aaaa(j,i,e,b) u2_aaaa(a,b,j,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += 0.500000 * sharedOps["aaaa_vooo_234"]("b,j,i,m") * l2_aaaa_oovv("j,i,e,b");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += -0.500000 <j,i||e,b>_abab u2_abab(a,b,j,i) m1_aa(m,a)
            // +            -0.500000 <i,j||e,b>_abab u2_abab(a,b,i,j) m1_aa(m,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["aa_vv_153"]("a,e") * m1_aa_ov("m,a");
        }

        if (include_u2_) {

            // rl1_aa_vo += 0.500000 dp_aa(m,b) t2_abab(b,a,i,j) m2_abab(i,j,e,a)
            // +            0.500000 dp_aa(m,b) t2_abab(b,a,j,i) m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["baab_vooo_226"]("a,i,m,j") * m2_abab_oovv("i,j,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += 1.000000 <m,k||b,c>_abab t2_bbbb(c,a,i,k) u1_aa(b,j) m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["baab_vooo_467"]("a,j,m,i") * m2_abab_oovv("j,i,e,a");

            // rl1_aa_vo += -0.500000 <m,k||c,b>_abab t2_abab(c,a,i,j) u1_bb(b,k) m2_abab(i,j,e,a)
            // +            -0.500000 <m,k||c,b>_abab t2_abab(c,a,j,i) u1_bb(b,k) m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["baab_vooo_469"]("a,m,i,j") * m2_abab_oovv("i,j,e,a");
        }

        if (include_u2_) {

            // rl1_aa_vo += 0.250000 <k,j||e,i>_abab u2_abab(b,a,k,j) m2_abab(m,i,b,a)
            // +            0.250000 <j,k||e,i>_abab u2_abab(b,a,j,k) m2_abab(m,i,b,a)
            // +            0.250000 <k,j||e,i>_abab u2_abab(a,b,k,j) m2_abab(m,i,a,b)
            // +            0.250000 <j,k||e,i>_abab u2_abab(a,b,j,k) m2_abab(m,i,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["aabb_vvvo_98"]("b,e,a,i") * m2_abab_oovv("m,i,b,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += -0.500000 <m,j||a,b>_abab u2_abab(a,b,i,j) m1_aa(i,e)
            // +            -0.500000 <m,j||b,a>_abab u2_abab(b,a,i,j) m1_aa(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["aa_oo_213"]("m,i") * m1_aa_ov("i,e");
        }

        if (include_u0_ && include_u2_) {

            // rl1_aa_vo += 0.500000 dp_aa(m,b) u0 u2_aaaa(b,a,i,j) m2_aaaa(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_aa_vo("e,m") += 0.500000 * sharedOps["aaaa_vooo_234"]("a,i,j,m") * m2_aaaa_oovv("i,j,e,a") * u0;
        }

        {

            // rl1_aa_vo += -0.500000 f_aa(m,b) l2_abab(i,j,e,a) t2_abab(b,a,i,j)
            // +            -0.500000 f_aa(m,b) l2_abab(j,i,e,a) t2_abab(b,a,j,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["baab_vooo_259"]("a,m,i,j") * l2_abab_oovv("i,j,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += 0.500000 <m,k||i,j>_abab u1_bb(a,k) m2_abab(i,j,e,a)
            // +            0.500000 <m,k||j,i>_abab u1_bb(a,k) m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["baab_vooo_305"]("a,m,i,j") * m2_abab_oovv("i,j,e,a");
        }

        {

            // rl1_aa_vo += 1.000000 <m,a||e,i>_abab l1_bb(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= V_blks_["baab_vovo"]("a,m,e,i") * l1_bb_ov("i,a");

            // rl1_aa_vo += -1.000000 <k,m||b,j>_aaaa l2_aaaa(i,j,e,a) t2_aaaa(b,a,i,k)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["aaaa_vooo_119"]("a,m,j,i") * l2_aaaa_oovv("i,j,e,a");

            // rl1_aa_vo += -1.000000 <b,j||e,c>_abab l2_aaaa(i,m,b,a) t2_abab(a,c,i,j)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["aaaa_vvvo_17"]("b,e,a,i") * l2_aaaa_oovv("i,m,b,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += -0.250000 <k,j||e,c>_aaaa t2_aaaa(b,a,k,j) u1_aa(c,i) m2_aaaa(i,m,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= 0.250000 * sharedOps["aaaa_vvvo_70"]("e,b,a,i") * m2_aaaa_oovv("i,m,b,a");
        }

        {

            // rl1_aa_vo += 1.000000 <m,k||j,b>_abab l2_abab(j,i,e,a) t2_bbbb(b,a,i,k)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["baab_vooo_111"]("a,m,j,i") * l2_abab_oovv("j,i,e,a");
        }

        if (include_u2_) {

            // rl1_aa_vo += 0.250000 <m,a||b,c>_aaaa u2_aaaa(b,c,i,j) m2_aaaa(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= 0.250000 * sharedOps["aaaa_vooo_96"]("a,i,j,m") * m2_aaaa_oovv("i,j,a,e");

            // rl1_aa_vo += 0.500000 dp_aa(m,b) t2_aaaa(b,a,i,j) m2_aaaa(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += 0.500000 * sharedOps["aaaa_vooo_215"]("a,m,i,j") * m2_aaaa_oovv("i,j,e,a");

            // rl1_aa_vo += 0.500000 dp_aa(i,e) l2_abab(m,j,a,b) u2_abab(a,b,i,j)
            // +            0.500000 dp_aa(i,e) l2_abab(m,j,b,a) u2_abab(b,a,i,j)
            // flops: o3v2L1: 1, o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2, o2v0L1: 1,
            rl1_aa_vo("e,m") += l2_abab_oovv("m,j,a,b") * u2_abab_vvoo("a,b,i,j") * dp_aa_ov("i,e");
        }

        if (include_u1_) {

            // rl1_aa_vo += 1.000000 <m,a||e,b>_aaaa u1_aa(b,i) m1_aa(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["aaaa_vvoo_165"]("a,e,m,i") * m1_aa_ov("i,a");
        }

        {

            // rl1_aa_vo += 0.500000 <b,a||e,i>_aaaa l2_aaaa(m,i,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += 0.500000 * V_blks_["aaaa_vvvo"]("b,a,e,i") * l2_aaaa_oovv("m,i,b,a");
        }

        if (include_u0_) {

            // rl1_aa_vo += 0.500000 dp_aa(m,b) l2_aaaa(i,j,e,a) t2_aaaa(b,a,i,j) u0
            // flops: o3v2L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_aa_vo("e,m") += 0.500000 * sharedOps["aaaa_vooo_215"]("a,m,i,j") * u0 * l2_aaaa_oovv("i,j,e,a");
        }

        {

            // rl1_aa_vo += 0.250000 <k,j||e,i>_aaaa l2_aaaa(m,i,b,a) t2_aaaa(b,a,k,j)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += 0.250000 * sharedOps["aaaa_vvvo_56"]("e,b,a,i") * l2_aaaa_oovv("m,i,b,a");
        }

        if (include_u1_) {

            // rl1_aa_vo += -1.000000 <m,j||i,a>_abab u1_bb(a,j) m1_aa(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["aa_oo_345"]("m,i") * m1_aa_ov("i,e");
        }

        {

            // rl1_aa_vo += 0.250000 <m,a||b,c>_aaaa l2_aaaa(i,j,a,e) t2_aaaa(b,c,i,j)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= 0.250000 * sharedOps["aaaa_vooo_43"]("a,m,i,j") * l2_aaaa_oovv("i,j,a,e");

            // rl1_aa_vo += -1.000000 <j,b||e,c>_abab l2_abab(m,i,a,b) t2_abab(a,c,j,i)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["aabb_vvvo_18"]("e,a,b,i") * l2_abab_oovv("m,i,a,b");
        }

        if (include_u2_) {

            // rl1_aa_vo += 0.500000 dp_aa(i,e) l2_aaaa(j,m,a,b) u2_aaaa(a,b,j,i)
            // flops: o3v2L1: 1, o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2, o2v0L1: 1,
            rl1_aa_vo("e,m") += 0.500000 * l2_aaaa_oovv("j,m,a,b") * u2_aaaa_vvoo("a,b,j,i") * dp_aa_ov("i,e");
        }

        if (include_u1_) {

            // rl1_aa_vo += 1.000000 <a,i||e,b>_abab u1_bb(b,i) m1_aa(m,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["aa_vv_316"]("a,e") * m1_aa_ov("m,a");
        }

        if (include_u2_) {

            // rl1_aa_vo += 0.500000 dp_aa(m,a) l2_abab(j,i,e,b) u2_abab(a,b,j,i)
            // +            0.500000 dp_aa(m,a) l2_abab(i,j,e,b) u2_abab(a,b,i,j)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["baab_vooo_223"]("b,j,m,i") * l2_abab_oovv("j,i,e,b");

            // rl1_aa_vo += -1.000000 <b,j||e,c>_abab u2_abab(a,c,i,j) m2_aaaa(i,m,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["aaaa_vvvo_23"]("a,b,e,i") * m2_aaaa_oovv("i,m,b,a");
        }

        if (include_u0_) {

            // rl1_aa_vo += 0.500000 dp_aa(m,b) l2_abab(i,j,e,a) t2_abab(b,a,i,j) u0
            // +            0.500000 dp_aa(m,b) l2_abab(j,i,e,a) t2_abab(b,a,j,i) u0
            // flops: o3v2L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_aa_vo("e,m") += sharedOps["baab_vooo_226"]("a,i,m,j") * u0 * l2_abab_oovv("i,j,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += -1.000000 <m,a||b,j>_abab u1_aa(b,i) m2_abab(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["baab_vooo_284"]("a,m,i,j") * m2_abab_oovv("i,j,e,a");
        }

        if (include_u1_) {

            // rl1_aa_vo += -1.000000 f_aa(m,a) u1_aa(a,i) m1_aa(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["aa_oo_367"]("i,m") * m1_aa_ov("i,e");
        }

        if (include_u2_) {

            // rl1_aa_vo += -1.000000 <k,m||b,j>_aaaa u2_aaaa(b,a,i,k) m2_aaaa(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["aaaa_vooo_138"]("a,i,m,j") * m2_aaaa_oovv("i,j,e,a");
        }

        {

            // rl1_aa_vo += -0.500000 f_aa(m,b) l2_aaaa(i,j,e,a) t2_aaaa(b,a,i,j)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= 0.500000 * sharedOps["aaaa_vooo_240"]("a,m,i,j") * l2_aaaa_oovv("i,j,e,a");

            // rl1_aa_vo += 1.000000 <m,k||b,j>_abab l2_abab(i,j,e,a) t2_abab(b,a,i,k)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["baab_vooo_118"]("a,m,i,j") * l2_abab_oovv("i,j,e,a");
        }

        if (include_u1_) {

            // rl1_aa_vo += 1.000000 dp_aa(m,a) l1_aa(i,e) u1_aa(a,i)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["aa_oo_312"]("i,m") * l1_aa_ov("i,e");
        }

        {

            // rl1_aa_vo += 0.250000 <k,j||e,i>_abab l2_abab(m,i,b,a) t2_abab(b,a,k,j)
            // +            0.250000 <j,k||e,i>_abab l2_abab(m,i,b,a) t2_abab(b,a,j,k)
            // +            0.250000 <k,j||e,i>_abab l2_abab(m,i,a,b) t2_abab(a,b,k,j)
            // +            0.250000 <j,k||e,i>_abab l2_abab(m,i,a,b) t2_abab(a,b,j,k)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["aabb_vvvo_53"]("e,b,a,i") * l2_abab_oovv("m,i,b,a");

            // rl1_aa_vo += -0.500000 <j,i||e,b>_abab l1_aa(m,a) t2_abab(a,b,j,i)
            // +            -0.500000 <i,j||e,b>_abab l1_aa(m,a) t2_abab(a,b,i,j)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["aa_vv_150"]("a,e") * l1_aa_ov("m,a");
        }

        if (include_u0_ && include_u1_) {

            // rl1_aa_vo += 1.000000 <m,i||e,a>_abab u1_bb(a,i) m0
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["aa_vo_37"]("e,m") * m0;
        }

        {

            // rl1_aa_vo += 1.000000 f_aa(a,e) l1_aa(m,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += F_blks_["aa_vv"]("a,e") * l1_aa_ov("m,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += 1.000000 <m,k||c,b>_abab t2_abab(c,a,i,k) u1_bb(b,j) m2_abab(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["baab_vooo_466"]("a,m,i,j") * m2_abab_oovv("i,j,e,a");
        }

        if (include_u2_) {

            // rl1_aa_vo += 1.000000 <m,k||j,b>_abab u2_abab(a,b,i,k) m2_aaaa(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["aaaa_vooo_137"]("a,m,j,i") * m2_aaaa_oovv("i,j,e,a");
        }

        if (include_u0_ && include_u1_) {

            // rl1_aa_vo += 1.000000 dp_aa(m,a) u0 u1_aa(a,i) m1_aa(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_aa_vo("e,m") += sharedOps["aa_oo_312"]("i,m") * m1_aa_ov("i,e") * u0;
        }

        {

            // rl1_aa_vo += 1.000000 f_aa(m,e)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rl1_aa_vo("e,m") += F_blks_["aa_ov"]("m,e");
        }

        if (include_u2_) {

            // rl1_aa_vo += -1.000000 <j,b||e,c>_aaaa u2_aaaa(c,a,i,j) m2_aaaa(i,m,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["aaaa_vvvo_29"]("b,e,a,i") * m2_aaaa_oovv("i,m,b,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += -0.250000 <k,m||b,c>_aaaa t2_aaaa(b,c,i,j) u1_aa(a,k) m2_aaaa(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= 0.250000 * sharedOps["aaaa_vooo_475"]("a,i,j,m") * m2_aaaa_oovv("i,j,e,a");
        }

        if (include_u1_) {

            // rl1_aa_vo += 1.000000 <m,a||e,b>_abab u1_bb(b,i) m1_bb(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["abab_vvoo_163"]("e,a,m,i") * m1_bb_ov("i,a");
        }

        {

            // rl1_aa_vo += -1.000000 <j,b||e,c>_aaaa l2_abab(m,i,b,a) t2_abab(c,a,j,i)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["aabb_vvvo_13"]("b,e,a,i") * l2_abab_oovv("m,i,b,a");
        }

        if (include_u0_ && include_u1_) {

            // rl1_aa_vo += -1.000000 <i,m||e,a>_aaaa u1_aa(a,i) m0
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["aa_vo_327"]("e,m") * m0;
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += -0.500000 <k,m||b,c>_aaaa t2_abab(c,a,i,j) u1_aa(b,k) m2_abab(i,j,e,a)
            // +            -0.500000 <k,m||b,c>_aaaa t2_abab(c,a,j,i) u1_aa(b,k) m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["baab_vooo_271"]("a,m,i,j") * m2_abab_oovv("i,j,e,a");

            // rl1_aa_vo += -1.000000 <m,a||j,b>_abab u1_bb(b,i) m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["baab_vooo_296"]("a,m,j,i") * m2_abab_oovv("j,i,e,a");
        }

        if (include_u1_) {

            // rl1_aa_vo += 1.000000 dp_aa(i,e) l1_aa(m,a) u1_aa(a,i)
            // flops: o2v1L1: 2, o1v1L1: 1 | mem: o1v1L1: 2, o2v0L1: 1,
            rl1_aa_vo("e,m") += l1_aa_ov("m,a") * u1_aa_vo("a,i") * dp_aa_ov("i,e");
        }

        {

            // rl1_aa_vo += -1.000000 f_aa(m,i) l1_aa(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= F_blks_["aa_oo"]("m,i") * l1_aa_ov("i,e");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += 0.500000 <j,i||e,b>_aaaa u2_aaaa(b,a,j,i) m1_aa(m,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += 0.500000 * sharedOps["aa_vv_171"]("a,e") * m1_aa_ov("m,a");
        }

        if (include_u2_) {

            // rl1_aa_vo += 1.000000 <m,k||j,b>_abab u2_bbbb(b,a,i,k) m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["baab_vooo_129"]("a,m,j,i") * m2_abab_oovv("j,i,e,a");
        }

        if (include_u0_ && include_u2_) {

            // rl1_aa_vo += 0.500000 dp_aa(m,b) u0 u2_abab(b,a,i,j) m2_abab(i,j,e,a)
            // +            0.500000 dp_aa(m,b) u0 u2_abab(b,a,j,i) m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 2, o0v0: 1,
            rl1_aa_vo("e,m") += sharedOps["baab_vooo_223"]("a,i,m,j") * m2_abab_oovv("i,j,e,a") * u0;
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += 1.000000 <k,m||b,c>_aaaa t2_abab(c,a,k,i) u1_aa(b,j) m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["baab_vooo_468"]("a,m,j,i") * m2_abab_oovv("j,i,e,a");
        }

        {

            // rl1_aa_vo += 0.500000 <b,a||e,i>_abab l2_abab(m,i,b,a)
            // +            0.500000 <a,b||e,i>_abab l2_abab(m,i,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += V_blks_["abab_vvvo"]("b,a,e,i") * l2_abab_oovv("m,i,b,a");
        }

        if (include_u0_) {

            // rl1_aa_vo += -1.000000 dp_aa(m,e) u0
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            rl1_aa_vo("e,m") -= dp_aa_ov("m,e") * u0;
        }

        if (include_u2_) {

            // rl1_aa_vo += -1.000000 <j,b||e,c>_abab u2_abab(a,c,j,i) m2_abab(m,i,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["aabb_vvvo_32"]("e,a,b,i") * m2_abab_oovv("m,i,a,b");

            // rl1_aa_vo += 0.250000 <k,j||e,i>_aaaa u2_aaaa(b,a,k,j) m2_aaaa(m,i,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += 0.250000 * sharedOps["aaaa_vvvo_94"]("b,a,e,i") * m2_aaaa_oovv("m,i,b,a");
        }

        {

            // rl1_aa_vo += -1.000000 <j,b||e,c>_aaaa l2_aaaa(i,m,b,a) t2_aaaa(c,a,i,j)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["aaaa_vvvo_19"]("a,b,e,i") * l2_aaaa_oovv("i,m,b,a");
        }

        if (include_u2_) {

            // rl1_aa_vo += -1.000000 <j,b||e,c>_aaaa u2_abab(c,a,j,i) m2_abab(m,i,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += sharedOps["aabb_vvvo_31"]("b,e,a,i") * m2_abab_oovv("m,i,b,a");
        }

        {

            // rl1_aa_vo += 0.500000 <j,i||e,b>_aaaa l1_aa(m,a) t2_aaaa(b,a,j,i)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += 0.500000 * sharedOps["aa_vv_151"]("a,e") * l1_aa_ov("m,a");

            // rl1_aa_vo += -1.000000 <b,j||e,c>_abab l2_abab(m,i,b,a) t2_bbbb(c,a,i,j)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["aabb_vvvo_12"]("b,e,a,i") * l2_abab_oovv("m,i,b,a");
        }

        if (include_u2_) {

            // rl1_aa_vo += -0.500000 f_aa(m,b) u2_aaaa(b,a,i,j) m2_aaaa(i,j,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= 0.500000 * sharedOps["aaaa_vooo_287"]("a,m,i,j") * m2_aaaa_oovv("i,j,e,a");

            // rl1_aa_vo += -1.000000 <b,j||e,c>_abab u2_bbbb(c,a,i,j) m2_abab(m,i,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["aabb_vvvo_22"]("b,e,a,i") * m2_abab_oovv("m,i,b,a");

            // rl1_aa_vo += -1.000000 <k,m||b,j>_aaaa u2_abab(b,a,k,i) m2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["baab_vooo_126"]("a,m,j,i") * m2_abab_oovv("j,i,e,a");
        }

        if (include_u1_ && include_u2_) {

            // rl1_aa_vo += 0.500000 <j,m||a,b>_aaaa u2_aaaa(a,b,i,j) m1_aa(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") += 0.500000 * sharedOps["aa_oo_244"]("i,m") * m1_aa_ov("i,e");
        }

        {

            // rl1_aa_vo += -0.500000 <m,j||a,b>_abab l1_aa(i,e) t2_abab(a,b,i,j)
            // +            -0.500000 <m,j||b,a>_abab l1_aa(i,e) t2_abab(b,a,i,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= sharedOps["aa_oo_207"]("m,i") * l1_aa_ov("i,e");
        }

        if (include_u0_) {

            // rl1_aa_vo += -1.000000 dp_aa(m,e) m0
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            rl1_aa_vo("e,m") -= dp_aa_ov("m,e") * m0;

            // cenergy += 1.000000 u0
            // flops: o0v0: 1 | mem: o0v0: 1,
            cenergy += u0;

            // cenergy += 1.000000 m0
            // flops: o0v0: 1 | mem: o0v0: 1,
            cenergy += m0;
        }

        if (include_u2_) {

            // energy += 0.125000 <l,k||i,j>_aaaa u2_aaaa(b,a,l,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.125000 * scalar_168;
        }

        if (include_u1_ && include_u2_) {

            // energy += 1.000000 dp_bb(k,c) u1_bb(c,j) u2_bbbb(b,a,i,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_169;

            // energy += -0.500000 <k,l||c,j>_abab t2_bbbb(b,a,i,l) u1_aa(c,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_170;

            // energy += 0.500000 <l,k||c,j>_aaaa t2_aaaa(b,a,i,l) u1_aa(c,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_171;
        }

        if (include_u1_) {

            // energy += 1.000000 <k,j||b,c>_bbbb t2_bbbb(c,a,i,k) u1_bb(b,j) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_172;
        }

        if (include_u1_ && include_u2_) {

            // energy += 0.500000 <l,k||c,j>_bbbb t2_bbbb(b,a,i,l) u1_bb(c,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_173;
        }

        if (include_u1_) {

            // energy += -1.000000 dp_aa(j,b) t2_abab(b,a,j,i) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_174;
        }

        if (include_u1_ && include_u2_) {

            // energy += -0.500000 <l,k||j,c>_abab t2_aaaa(b,a,i,l) u1_bb(c,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_175;
        }

        if (include_u1_) {

            // energy += -1.000000 <k,j||b,c>_aaaa t2_abab(c,a,k,i) u1_aa(b,j) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_176;
        }

        if (include_u2_) {

            // energy += -0.250000 <l,k||c,d>_aaaa t2_aaaa(b,a,i,l) u2_aaaa(c,d,j,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.250000 * scalar_177;
        }

        if (include_u1_) {

            // energy += 1.000000 <k,j||b,c>_aaaa t2_aaaa(c,a,i,k) u1_aa(b,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_178;
        }

        if (include_u1_ && include_u2_) {

            // energy += -0.500000 f_bb(k,c) t2_bbbb(b,a,i,k) u1_bb(c,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_179;
        }

        {

            // energy += 1.000000 f_aa(j,b) l1_bb(i,a) t2_abab(b,a,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_180;

            // energy += -1.000000 f_bb(j,b) l1_bb(i,a) t2_bbbb(b,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_181;
        }

        if (include_u1_) {

            // energy += 1.000000 <j,a||b,i>_aaaa u1_aa(b,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_182;
        }

        if (include_u1_ && include_u2_) {

            // energy += -1.000000 f_bb(j,b) u2_bbbb(b,a,i,j) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_183;

            // energy += 1.000000 f_aa(j,b) u2_abab(b,a,j,i) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_184;
        }

        {

            // energy += 1.000000 f_bb(j,b) l1_aa(i,a) t2_abab(a,b,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_185;
        }

        if (include_u2_) {

            // energy += -0.250000 <l,k||c,d>_bbbb t2_bbbb(b,a,i,l) u2_bbbb(c,d,j,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.250000 * scalar_186;
        }

        if (include_u1_) {

            // energy += 0.500000 <k,j||b,c>_aaaa t2_aaaa(c,a,k,j) u1_aa(b,i) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_187;
        }

        if (include_u1_ && include_u2_) {

            // energy += 1.000000 dp_aa(k,c) u1_aa(c,j) u2_aaaa(b,a,i,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_188;
        }

        if (include_u1_) {

            // energy += 0.500000 dp_bb(k,c) l2_bbbb(i,j,b,a) t2_bbbb(b,a,i,k) u1_bb(c,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_189;

            // energy += 0.500000 dp_aa(k,c) l2_aaaa(i,j,b,a) t2_aaaa(b,a,i,k) u1_aa(c,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_190;
        }

        if (include_u2_) {

            // energy += 0.500000 dp_aa(k,j) t2_aaaa(b,a,i,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_191;
        }

        {

            // energy += 0.500000 f_aa(b,c) l2_aaaa(i,j,b,a) t2_aaaa(c,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_192;

            // energy += -0.250000 <l,k||c,d>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(c,d,j,k) t2_bbbb(b,a,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.250000 * scalar_193;

            // energy += -0.250000 <l,k||c,d>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(c,d,j,k) t2_aaaa(b,a,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.250000 * scalar_194;

            // energy += -0.500000 <j,a||b,c>_bbbb l1_bb(i,a) t2_bbbb(b,c,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_195;
        }

        if (include_u1_ && include_u2_) {

            // energy += 0.500000 <b,a||c,j>_bbbb u1_bb(c,i) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_196;

            // energy += -0.500000 <j,a||b,c>_bbbb u2_bbbb(b,c,i,j) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_197;
        }

        if (include_u1_) {

            // energy += -1.000000 <k,j||b,c>_bbbb t2_abab(a,c,i,k) u1_bb(b,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_198;
        }

        {

            // energy += -1.000000 f_aa(j,b) l1_aa(i,a) t2_aaaa(b,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_199;
        }

        if (include_u2_) {

            // energy += 0.500000 f_aa(b,c) u2_aaaa(c,a,i,j) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_200;
        }

        if (include_u1_ && include_u2_) {

            // energy += -0.500000 <j,a||b,c>_aaaa u2_aaaa(b,c,i,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_201;
        }

        if (include_u2_) {

            // energy += 0.500000 dp_aa(j,i) l2_aaaa(k,i,a,b) u2_aaaa(a,b,k,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_202;

            // energy += 0.500000 dp_bb(k,j) t2_bbbb(b,a,i,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_203;

            // energy += 0.500000 dp_bb(j,i) l2_bbbb(k,i,a,b) u2_bbbb(a,b,k,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_204;
        }

        {

            // energy += 0.125000 <l,k||i,j>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(b,a,l,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.125000 * scalar_205;
        }

        if (include_u1_) {

            // energy += 1.000000 dp_bb(j,b) t2_bbbb(b,a,i,j) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_206;
        }

        if (include_u2_) {

            // energy += 1.000000 dp_aa(i,a) l1_aa(j,b) u2_aaaa(a,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_207;
        }

        {

            // energy += -0.500000 <j,a||b,c>_aaaa l1_aa(i,a) t2_aaaa(b,c,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_208;

            // energy += 0.500000 f_bb(b,c) l2_bbbb(i,j,b,a) t2_bbbb(c,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_209;
        }

        if (include_u2_) {

            // energy += -1.000000 dp_bb(i,a) l1_aa(j,b) u2_abab(b,a,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_210;

            // energy += 1.000000 dp_bb(i,a) l1_bb(j,b) u2_bbbb(a,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_211;
        }

        if (include_u1_ && include_u2_) {

            // energy += 1.000000 dp_aa(k,j) u1_aa(b,k) u1_bb(a,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_212;

            // energy += 1.000000 dp_bb(k,j) u1_bb(b,k) u1_aa(a,i) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_213;
        }

        if (include_u2_) {

            // energy += 1.000000 <l,k||c,d>_aaaa t2_abab(d,a,l,i) u2_abab(c,b,k,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_214;
        }

        if (include_u1_ && include_u2_) {

            // energy += -1.000000 dp_bb(k,j) u1_bb(b,k) u1_bb(a,i) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_215;

            // energy += -1.000000 dp_bb(b,c) u1_bb(c,j) u1_aa(a,i) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_216;

            // energy += 1.000000 dp_aa(b,c) u1_aa(c,j) u1_aa(a,i) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_217;

            // energy += -1.000000 dp_aa(b,c) u1_aa(c,j) u1_bb(a,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_218;

            // energy += 1.000000 dp_bb(b,c) u1_bb(c,j) u1_bb(a,i) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_219;
        }

        if (include_u1_) {

            // energy += 1.000000 <j,a||b,i>_bbbb u1_bb(b,j) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_220;
        }

        if (include_u1_ && include_u2_) {

            // energy += 1.000000 f_bb(j,b) u2_abab(a,b,i,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_221;

            // energy += -1.000000 f_aa(j,b) u2_aaaa(b,a,i,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_222;
        }

        if (include_u1_) {

            // energy += 1.000000 <j,a||b,i>_abab u1_aa(b,j) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_223;
        }

        if (include_u1_ && include_u2_) {

            // energy += -1.000000 dp_aa(k,j) u1_aa(b,k) u1_aa(a,i) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_224;
        }

        {

            // energy += 0.500000 <l,k||c,d>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(c,b,j,k) t2_bbbb(d,a,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_225;
        }

        if (include_u2_) {

            // energy += 1.000000 <k,l||c,d>_abab t2_bbbb(d,a,i,l) u2_abab(c,b,k,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_226;
        }

        {

            // energy += -0.500000 f_bb(k,j) l2_bbbb(i,j,b,a) t2_bbbb(b,a,i,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_227;

            // energy += 0.500000 <l,k||c,d>_bbbb l2_abab(i,j,a,b) t2_bbbb(c,b,j,k) t2_abab(a,d,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_228;
        }

        if (include_u2_) {

            // energy += 1.000000 <l,k||d,c>_abab t2_abab(d,a,l,i) u2_abab(b,c,j,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_229;
        }

        {

            // energy += 0.500000 <l,k||c,d>_bbbb l2_abab(j,i,b,a) t2_abab(b,c,j,k) t2_bbbb(d,a,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_230;

            // energy += 0.500000 <k,l||c,d>_abab l2_bbbb(i,j,b,a) t2_abab(c,b,k,j) t2_bbbb(d,a,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_231;

            // energy += 0.500000 <l,k||d,c>_abab l2_bbbb(i,j,b,a) t2_bbbb(c,b,j,k) t2_abab(d,a,l,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_232;
        }

        if (include_u2_) {

            // energy += -0.250000 <l,k||c,d>_aaaa t2_aaaa(d,a,i,j) u2_aaaa(c,b,l,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.250000 * scalar_233;

            // energy += -0.500000 dp_bb(a,b) l2_bbbb(j,i,a,c) u2_bbbb(b,c,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_234;

            // energy += -0.500000 dp_aa(a,b) l2_aaaa(j,i,a,c) u2_aaaa(b,c,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_235;

            // energy += -0.500000 dp_bb(b,c) t2_bbbb(c,a,i,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_236;
        }

        if (include_u1_) {

            // energy += -1.000000 f_bb(j,i) u1_bb(a,j) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_237;

            // energy += -1.000000 f_aa(j,i) u1_aa(a,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_238;

            // energy += 1.000000 dp_bb(j,i) l1_bb(i,a) u1_bb(a,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_239;

            // energy += 1.000000 dp_aa(j,i) l1_aa(i,a) u1_aa(a,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_240;
        }

        if (include_u2_) {

            // energy += -0.500000 dp_aa(b,c) t2_aaaa(c,a,i,j) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_241;
        }

        if (include_u1_) {

            // energy += 1.000000 f_bb(a,b) u1_bb(b,i) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_242;
        }

        if (include_u2_) {

            // energy += 1.000000 <k,l||c,d>_abab t2_abab(a,d,i,l) u2_aaaa(c,b,j,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_243;

            // energy += 0.062500 <l,k||c,d>_bbbb t2_bbbb(b,a,l,k) u2_bbbb(c,d,i,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.062500 * scalar_244;
        }

        if (include_u1_ && include_u2_) {

            // energy += 0.250000 <l,k||c,j>_aaaa t2_aaaa(b,a,l,k) u1_aa(c,i) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.250000 * scalar_245;

            // energy += 0.250000 <l,k||c,j>_bbbb t2_bbbb(b,a,l,k) u1_bb(c,i) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.250000 * scalar_246;
        }

        {

            // energy += 0.125000 <l,k||i,j>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(b,a,l,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.125000 * scalar_247;
        }

        if (include_u2_) {

            // energy += 0.125000 <l,k||i,j>_bbbb u2_bbbb(b,a,l,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.125000 * scalar_248;
        }

        {

            // energy += -0.250000 <l,k||c,d>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(c,b,l,k) t2_bbbb(d,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.250000 * scalar_249;
        }

        if (include_u1_) {

            // energy += 0.500000 <k,j||b,c>_bbbb t2_bbbb(c,a,k,j) u1_bb(b,i) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_250;
        }

        if (include_u2_) {

            // energy += -0.250000 <l,k||c,d>_bbbb t2_bbbb(d,a,i,j) u2_bbbb(c,b,l,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.250000 * scalar_251;

            // energy += 1.000000 <l,k||c,d>_aaaa t2_aaaa(d,a,i,l) u2_abab(c,b,k,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_252;

            // energy += 1.000000 <l,k||c,d>_aaaa t2_abab(d,a,l,i) u2_aaaa(c,b,j,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_253;
        }

        {

            // energy += 0.500000 <l,k||c,d>_aaaa l2_abab(j,i,b,a) t2_aaaa(c,b,j,k) t2_abab(d,a,l,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_254;

            // energy += 0.500000 <l,k||c,d>_aaaa l2_abab(i,j,a,b) t2_abab(c,b,k,j) t2_aaaa(d,a,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_255;

            // energy += 0.500000 <l,k||c,d>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(c,b,j,k) t2_aaaa(d,a,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_256;

            // energy += 0.500000 <l,k||c,d>_bbbb l2_aaaa(i,j,b,a) t2_abab(b,c,j,k) t2_abab(a,d,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_257;
        }

        if (include_u2_) {

            // energy += 1.000000 <k,l||d,c>_abab t2_abab(d,a,i,l) u2_abab(b,c,k,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_258;
        }

        if (include_u1_) {

            // energy += 1.000000 dp_aa(a,i) l2_aaaa(j,i,a,b) u1_aa(b,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_259;

            // energy += -1.000000 dp_bb(a,i) l2_abab(j,i,b,a) u1_aa(b,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_260;

            // energy += 1.000000 dp_bb(a,i) l2_bbbb(j,i,a,b) u1_bb(b,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_261;

            // energy += -1.000000 dp_aa(a,i) l2_abab(i,j,a,b) u1_bb(b,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_262;

            // energy += -1.000000 dp_aa(a,b) l1_aa(i,a) u1_aa(b,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_263;
        }

        {

            // energy += -0.250000 <l,k||c,d>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(c,b,l,k) t2_aaaa(d,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.250000 * scalar_264;
        }

        if (include_u1_) {

            // energy += -1.000000 dp_bb(a,b) l1_bb(i,a) u1_bb(b,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_265;
        }

        if (include_u2_) {

            // energy += 1.000000 <l,k||c,d>_abab t2_abab(a,d,l,i) u2_abab(c,b,j,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_266;
        }

        if (include_u1_) {

            // energy += 1.000000 <a,j||i,b>_abab u1_bb(b,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_267;
        }

        {

            // energy += 0.500000 <l,k||c,d>_aaaa l2_bbbb(i,j,b,a) t2_abab(c,b,k,j) t2_abab(d,a,l,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_268;
        }

        if (include_u2_) {

            // energy += 1.000000 <l,k||c,d>_bbbb t2_abab(a,d,i,l) u2_bbbb(c,b,j,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_269;

            // energy += 1.000000 <l,k||d,c>_abab t2_abab(d,a,l,i) u2_bbbb(c,b,j,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_270;

            // energy += 1.000000 <k,l||c,d>_abab t2_abab(a,d,i,l) u2_abab(c,b,k,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_271;

            // energy += 1.000000 <l,k||c,d>_bbbb t2_abab(a,d,i,l) u2_abab(b,c,j,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_272;

            // energy += 1.000000 <k,l||c,d>_abab t2_bbbb(d,a,i,l) u2_aaaa(c,b,j,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_273;

            // energy += 1.000000 <l,k||c,d>_aaaa t2_aaaa(d,a,i,l) u2_aaaa(c,b,j,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_274;

            // energy += 1.000000 <l,k||c,d>_bbbb t2_bbbb(d,a,i,l) u2_bbbb(c,b,j,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_275;

            // energy += 1.000000 <l,k||c,d>_bbbb t2_bbbb(d,a,i,l) u2_abab(b,c,j,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_276;
        }

        if (include_u0_ && include_u2_) {

            // energy += 0.250000 <j,i||a,b>_aaaa u2_aaaa(a,b,j,i) m0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += 0.250000 * scalar_11 * m0;
        }

        {

            // energy += 1.000000 <k,b||c,j>_aaaa l2_abab(j,i,b,a) t2_abab(c,a,k,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_277;
        }

        if (include_u1_) {

            // energy += -1.000000 dp_bb(a,i) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_278;

            // energy += -1.000000 dp_aa(a,i) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_279;
        }

        if (include_u0_) {

            // energy += -1.000000 dp_bb(i,i) m0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_5 * m0;
        }

        {

            // energy += 0.062500 <l,k||c,d>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(c,d,i,j) t2_bbbb(b,a,l,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.062500 * scalar_280;
        }

        if (include_u2_) {

            // energy += 0.062500 <l,k||c,d>_bbbb t2_bbbb(c,d,i,j) u2_bbbb(b,a,l,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.062500 * scalar_281;
        }

        if (include_u0_) {

            // energy += -1.000000 dp_aa(i,i) m0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_4 * m0;
        }

        if (include_u1_) {

            // energy += 1.000000 f_aa(a,b) u1_aa(b,i) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_282;
        }

        if (include_u2_) {

            // energy += 0.062500 <l,k||c,d>_aaaa t2_aaaa(b,a,l,k) u2_aaaa(c,d,i,j) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.062500 * scalar_283;

            // energy += 0.062500 <l,k||c,d>_aaaa t2_aaaa(c,d,i,j) u2_aaaa(b,a,l,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.062500 * scalar_284;
        }

        {

            // energy += -1.000000 <k,b||j,c>_abab l2_abab(j,i,a,b) t2_abab(a,c,k,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_285;
        }

        if (include_u2_) {

            // energy += -1.000000 <k,b||c,j>_abab u2_aaaa(c,a,i,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_286;

            // energy += -1.000000 <b,k||c,j>_abab u2_abab(c,a,i,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_287;

            // energy += 1.000000 <k,b||c,j>_bbbb u2_bbbb(c,a,i,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_288;
        }

        {

            // energy += -1.000000 <k,b||c,j>_abab l2_bbbb(i,j,b,a) t2_abab(c,a,k,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_289;

            // energy += 1.000000 <k,b||c,j>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(c,a,i,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_290;

            // energy += 0.062500 <l,k||c,d>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(c,d,i,j) t2_aaaa(b,a,l,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.062500 * scalar_291;
        }

        if (include_u2_) {

            // energy += -0.250000 <l,k||c,d>_aaaa t2_aaaa(c,d,i,l) u2_aaaa(b,a,j,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.250000 * scalar_292;
        }

        if (include_u1_ && include_u2_) {

            // energy += -1.000000 <l,k||c,j>_aaaa t2_aaaa(c,a,i,l) u1_aa(b,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_293;

            // energy += -1.000000 <l,k||c,j>_bbbb t2_abab(a,c,i,l) u1_bb(b,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_294;
        }

        if (include_u1_) {

            // energy += -1.000000 <k,j||c,b>_abab t2_aaaa(c,a,i,k) u1_bb(b,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_295;

            // energy += -1.000000 <j,k||b,c>_abab t2_bbbb(c,a,i,k) u1_aa(b,j) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_296;

            // energy += 1.000000 <k,j||c,b>_abab t2_abab(c,a,k,i) u1_bb(b,j) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_297;

            // energy += 1.000000 <j,k||b,c>_abab t2_abab(a,c,i,k) u1_aa(b,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_298;
        }

        if (include_u1_ && include_u2_) {

            // energy += 1.000000 dp_bb(k,c) u1_aa(b,j) u2_bbbb(c,a,i,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_299;
        }

        {

            // energy += 0.250000 <j,i||b,a>_abab t2_abab(b,a,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.250000 * scalar_15;
        }

        if (include_u1_ && include_u2_) {

            // energy += -1.000000 dp_bb(k,c) u1_bb(b,j) u2_abab(a,c,i,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_300;

            // energy += -1.000000 dp_aa(k,c) u1_aa(b,j) u2_aaaa(c,a,i,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_301;
        }

        if (include_u0_) {

            // energy += 1.000000 u0 m0 w0
            // flops: o0v0: 3 | mem: o0v0: 3,
            energy += u0 * m0 * w0;
        }

        {

            // energy += 0.250000 <j,i||a,b>_abab t2_abab(a,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.250000 * scalar_15;

            // energy += 0.250000 <j,i||a,b>_aaaa t2_aaaa(a,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.250000 * scalar_14;

            // energy += 0.250000 <b,a||i,j>_bbbb l2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.250000 * scalar_302;
        }

        if (include_u0_ && include_u2_) {

            // energy += 0.250000 <j,i||b,a>_abab u2_abab(b,a,j,i) m0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += 0.250000 * scalar_12 * m0;

            // energy += 0.250000 <i,j||a,b>_abab u2_abab(a,b,i,j) m0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += 0.250000 * scalar_12 * m0;
        }

        if (include_u1_) {

            // energy += 1.000000 dp_aa(k,c) l2_abab(i,j,a,b) t2_aaaa(c,a,i,k) u1_bb(b,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_303;
        }

        if (include_u2_) {

            // energy += -1.000000 <k,b||j,c>_abab u2_abab(a,c,k,i) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_304;

            // energy += -1.000000 <k,b||c,j>_abab u2_abab(c,a,k,i) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_305;
        }

        {

            // energy += 1.000000 <k,b||c,j>_bbbb l2_abab(i,j,a,b) t2_abab(a,c,i,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_306;
        }

        if (include_u2_) {

            // energy += -0.250000 <l,k||c,d>_bbbb t2_bbbb(d,a,l,k) u2_bbbb(c,b,i,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.250000 * scalar_307;

            // energy += -0.250000 <l,k||c,d>_aaaa t2_aaaa(d,a,l,k) u2_aaaa(c,b,i,j) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.250000 * scalar_308;

            // energy += 1.000000 <k,b||c,j>_bbbb u2_abab(a,c,i,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_309;
        }

        if (include_u1_ && include_u2_) {

            // energy += -1.000000 <b,k||c,d>_abab t2_abab(a,d,i,k) u1_aa(c,j) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_310;

            // energy += -1.000000 <b,k||c,d>_abab t2_bbbb(d,a,i,k) u1_aa(c,j) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_311;
        }

        if (include_u0_ && include_u1_) {

            // energy += 1.000000 f_bb(i,a) u1_bb(a,i) m0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_3 * m0;
        }

        if (include_u1_ && include_u2_) {

            // energy += -1.000000 <k,b||d,c>_abab t2_aaaa(d,a,i,k) u1_bb(c,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_312;

            // energy += -1.000000 <k,b||c,d>_abab t2_abab(a,d,k,i) u1_aa(c,j) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_313;
        }

        {

            // energy += 0.500000 <k,l||c,d>_abab l2_aaaa(i,j,b,a) t2_aaaa(c,b,j,k) t2_abab(a,d,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_314;

            // energy += 0.500000 <l,k||d,c>_abab l2_aaaa(i,j,b,a) t2_abab(b,c,j,k) t2_aaaa(d,a,i,l)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_315;
        }

        if (include_u1_ && include_u2_) {

            // energy += -1.000000 <k,b||d,c>_abab t2_abab(d,a,k,i) u1_bb(c,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_316;

            // energy += -1.000000 <k,b||c,d>_aaaa t2_abab(d,a,k,i) u1_aa(c,j) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_317;

            // energy += -1.000000 <b,k||d,c>_abab t2_abab(d,a,i,k) u1_bb(c,j) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_318;

            // energy += -1.000000 <k,b||c,d>_aaaa t2_aaaa(d,a,i,k) u1_aa(c,j) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_319;

            // energy += -1.000000 <k,b||c,d>_bbbb t2_abab(a,d,i,k) u1_bb(c,j) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_320;
        }

        if (include_u0_ && include_u1_) {

            // energy += 1.000000 f_aa(i,a) u1_aa(a,i) m0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += scalar_2 * m0;
        }

        {

            // energy += 1.000000 f_bb(i,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_1;
        }

        if (include_u0_ && include_u1_) {

            // energy += -1.000000 dp_bb(i,a) u0 u1_bb(a,i) m0
            // flops: o0v0: 3 | mem: o0v0: 3,
            energy -= scalar_7 * u0 * m0;
        }

        if (include_u0_ && include_u2_) {

            // energy += 0.250000 <j,i||a,b>_abab u2_abab(a,b,j,i) m0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += 0.250000 * scalar_12 * m0;
        }

        if (include_u1_ && include_u2_) {

            // energy += -1.000000 <k,b||c,d>_bbbb t2_bbbb(d,a,i,k) u1_bb(c,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_321;
        }

        {

            // energy += 1.000000 f_aa(i,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_0;

            // energy += 1.000000 f_bb(a,i) l1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_322;

            // energy += 1.000000 f_aa(a,i) l1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_323;

            // energy += 1.000000 <k,b||c,j>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(c,a,i,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_324;
        }

        if (include_u2_) {

            // energy += 1.000000 <l,k||d,c>_abab t2_aaaa(d,a,i,l) u2_bbbb(c,b,j,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_325;

            // energy += -1.000000 <b,k||j,c>_abab u2_abab(a,c,i,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_326;

            // energy += 1.000000 <k,b||c,j>_aaaa u2_aaaa(c,a,i,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_327;
        }

        {

            // energy += -1.000000 <k,b||c,j>_abab l2_abab(i,j,a,b) t2_aaaa(c,a,i,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_328;

            // energy += -1.000000 <b,k||j,c>_abab l2_aaaa(i,j,b,a) t2_abab(a,c,i,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_329;

            // energy += -1.000000 <b,k||c,j>_abab l2_abab(i,j,b,a) t2_abab(c,a,i,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_330;
        }

        if (include_u2_) {

            // energy += 1.000000 <k,b||c,j>_aaaa u2_abab(c,a,k,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_331;

            // energy += -1.000000 <b,k||j,c>_abab u2_bbbb(c,a,i,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_332;

            // energy += 1.000000 <l,k||d,c>_abab t2_aaaa(d,a,i,l) u2_abab(b,c,j,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_333;

            // energy += 0.500000 f_bb(b,c) u2_bbbb(c,a,i,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_334;
        }

        if (include_u0_ && include_u1_) {

            // energy += -1.000000 dp_aa(i,a) u0 u1_aa(a,i) m0
            // flops: o0v0: 3 | mem: o0v0: 3,
            energy -= scalar_6 * u0 * m0;
        }

        {

            // energy += 0.125000 <b,a||c,d>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(c,d,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.125000 * scalar_335;

            // energy += 0.250000 <b,a||i,j>_aaaa l2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.250000 * scalar_336;

            // energy += -0.500000 <j,i||j,i>_bbbb
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_10;

            // energy += -0.500000 <i,j||i,j>_abab
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_9;

            // energy += -0.500000 <j,i||j,i>_abab
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_9;

            // energy += -0.500000 <j,i||j,i>_aaaa
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_8;
        }

        if (include_u1_ && include_u2_) {

            // energy += 1.000000 <k,l||j,c>_abab t2_abab(a,c,i,l) u1_aa(b,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_337;

            // energy += 0.250000 <k,b||c,d>_aaaa t2_aaaa(c,d,i,j) u1_aa(a,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.250000 * scalar_338;

            // energy += 1.000000 dp_bb(k,c) u1_bb(b,k) u2_bbbb(c,a,i,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_339;

            // energy += 0.250000 <k,b||c,d>_bbbb t2_bbbb(c,d,i,j) u1_bb(a,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.250000 * scalar_340;

            // energy += 1.000000 <k,l||j,c>_abab t2_bbbb(c,a,i,l) u1_aa(b,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_341;

            // energy += -0.500000 f_bb(k,c) t2_bbbb(c,a,i,j) u1_bb(b,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_342;
        }

        if (include_u1_) {

            // energy += 0.500000 dp_aa(k,c) l2_aaaa(i,j,b,a) t2_aaaa(c,a,i,j) u1_aa(b,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_343;
        }

        if (include_u0_ && include_u2_) {

            // energy += 0.250000 <j,i||a,b>_bbbb u2_bbbb(a,b,j,i) m0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += 0.250000 * scalar_13 * m0;

            // energy += 0.250000 <i,j||b,a>_abab u2_abab(b,a,i,j) m0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy += 0.250000 * scalar_12 * m0;
        }

        if (include_u1_ && include_u2_) {

            // energy += 1.000000 dp_aa(k,c) u1_aa(b,k) u2_aaaa(c,a,i,j) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_344;

            // energy += -0.500000 f_aa(k,c) t2_aaaa(b,a,i,k) u1_aa(c,j) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_345;
        }

        {

            // energy += -0.500000 f_aa(k,j) l2_aaaa(i,j,b,a) t2_aaaa(b,a,i,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_346;
        }

        if (include_u2_) {

            // energy += -0.500000 f_aa(k,j) u2_aaaa(b,a,i,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_347;
        }

        {

            // energy += -0.500000 <k,j||b,i>_bbbb l1_bb(i,a) t2_bbbb(b,a,k,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_348;
        }

        if (include_u1_ && include_u2_) {

            // energy += 0.500000 <k,b||i,j>_aaaa u1_aa(a,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_349;
        }

        if (include_u2_) {

            // energy += -0.500000 f_bb(k,j) u2_bbbb(b,a,i,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_350;
        }

        if (include_u1_ && include_u2_) {

            // energy += -0.500000 <k,j||b,i>_aaaa u2_aaaa(b,a,k,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_351;
        }

        {

            // energy += -0.500000 <k,j||b,i>_aaaa l1_aa(i,a) t2_aaaa(b,a,k,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_352;
        }

        if (include_u1_ && include_u2_) {

            // energy += 0.500000 <b,k||d,c>_abab t2_aaaa(d,a,i,j) u1_bb(c,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_353;

            // energy += -0.500000 <k,j||b,i>_bbbb u2_bbbb(b,a,k,j) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_354;

            // energy += 0.500000 <b,a||c,j>_aaaa u1_aa(c,i) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_355;
        }

        if (include_u1_) {

            // energy += -1.000000 dp_bb(j,b) t2_abab(a,b,i,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_356;
        }

        if (include_u2_) {

            // energy += -1.000000 dp_aa(i,a) l1_bb(j,b) u2_abab(a,b,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_357;
        }

        if (include_u1_ && include_u2_) {

            // energy += 0.500000 <k,b||c,d>_abab t2_bbbb(d,a,i,j) u1_aa(c,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_358;
        }

        if (include_u1_) {

            // energy += 1.000000 dp_aa(j,b) t2_aaaa(b,a,i,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_359;
        }

        if (include_u1_ && include_u2_) {

            // energy += 0.500000 <k,b||c,d>_bbbb t2_bbbb(d,a,i,j) u1_bb(c,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_360;

            // energy += 0.500000 <k,b||c,d>_aaaa t2_aaaa(d,a,i,j) u1_aa(c,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_361;

            // energy += 0.500000 <k,b||i,j>_bbbb u1_bb(a,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_362;
        }

        if (include_u1_) {

            // energy += 1.000000 dp_bb(k,c) l2_aaaa(i,j,b,a) t2_abab(a,c,i,k) u1_aa(b,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_363;
        }

        if (include_u1_ && include_u2_) {

            // energy += -1.000000 dp_aa(k,c) u1_aa(b,j) u2_abab(c,a,k,i) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_364;
        }

        if (include_u1_) {

            // energy += -1.000000 dp_bb(k,c) l2_abab(i,j,a,b) t2_abab(a,c,i,k) u1_bb(b,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_365;
        }

        if (include_u1_ && include_u2_) {

            // energy += 1.000000 dp_aa(k,c) u1_bb(b,j) u2_aaaa(c,a,i,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_366;
        }

        if (include_u1_) {

            // energy += -1.000000 dp_bb(k,c) l2_bbbb(i,j,b,a) t2_bbbb(c,a,i,k) u1_bb(b,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_367;
        }

        {

            // energy += -1.000000 <b,k||j,c>_abab l2_abab(j,i,b,a) t2_bbbb(c,a,i,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_368;
        }

        if (include_u1_) {

            // energy += 2.000000 dp_bb(j,b) u1_bb(b,i) u1_bb(a,j) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 2.000000 * scalar_369;

            // energy += 2.000000 dp_aa(j,b) u1_aa(b,i) u1_aa(a,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 2.000000 * scalar_370;
        }

        {

            // energy += 0.250000 <j,i||a,b>_bbbb t2_bbbb(a,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.250000 * scalar_16;
        }

        if (include_u1_) {

            // energy += 0.500000 <k,j||b,c>_bbbb t2_bbbb(b,c,i,k) u1_bb(a,j) m1_bb(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_371;
        }

        if (include_u1_ && include_u2_) {

            // energy += 1.000000 <l,k||j,c>_abab t2_abab(a,c,l,i) u1_bb(b,k) m2_abab(j,i,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_372;

            // energy += 1.000000 <l,k||c,j>_abab t2_abab(c,a,l,i) u1_bb(b,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_373;

            // energy += 1.000000 <k,l||c,j>_abab t2_abab(c,a,i,l) u1_aa(b,k) m2_abab(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_374;
        }

        if (include_u2_) {

            // energy += -0.250000 <l,k||c,d>_bbbb t2_bbbb(c,d,i,l) u2_bbbb(b,a,j,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.250000 * scalar_375;
        }

        if (include_u1_ && include_u2_) {

            // energy += -1.000000 <l,k||c,j>_bbbb t2_bbbb(c,a,i,l) u1_bb(b,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_376;

            // energy += 1.000000 <l,k||c,j>_abab t2_aaaa(c,a,i,l) u1_bb(b,k) m2_abab(i,j,a,b)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_377;

            // energy += -1.000000 <l,k||c,j>_aaaa t2_abab(c,a,l,i) u1_aa(b,k) m2_abab(j,i,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_378;
        }

        if (include_u1_) {

            // energy += 0.500000 <k,j||b,c>_aaaa t2_aaaa(b,c,i,k) u1_aa(a,j) m1_aa(i,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_379;

            // energy += -1.000000 dp_bb(i,a) u1_bb(a,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_7;

            // energy += -1.000000 dp_aa(i,a) u1_aa(a,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_6;
        }

        if (include_u2_) {

            // energy += 0.125000 <b,a||c,d>_bbbb u2_bbbb(c,d,i,j) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.125000 * scalar_380;
        }

        {

            // energy += 0.250000 <i,j||b,a>_abab t2_abab(b,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.250000 * scalar_15;
        }

        if (include_u1_) {

            // energy += 0.500000 dp_bb(k,c) l2_bbbb(i,j,b,a) t2_bbbb(c,a,i,j) u1_bb(b,k)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.500000 * scalar_381;
        }

        {

            // energy += 0.125000 <b,a||c,d>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(c,d,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.125000 * scalar_382;
        }

        if (include_u1_ && include_u2_) {

            // energy += -0.500000 f_aa(k,c) t2_aaaa(c,a,i,j) u1_aa(b,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar_383;
        }

        {

            // energy += 0.250000 <i,j||a,b>_abab t2_abab(a,b,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.250000 * scalar_15;
        }

        if (include_u1_) {

            // energy += -1.000000 dp_aa(k,c) l2_aaaa(i,j,b,a) t2_aaaa(c,a,i,k) u1_aa(b,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_384;
        }

        if (include_u0_) {

            // energy += -1.000000 dp_bb(i,i) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_5 * u0;
        }

        if (include_u2_) {

            // energy += 0.125000 <b,a||c,d>_aaaa u2_aaaa(c,d,i,j) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.125000 * scalar_385;
        }

        if (include_u1_ && include_u2_) {

            // energy += 1.000000 dp_aa(k,c) u1_bb(b,j) u2_abab(c,a,k,i) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_386;

            // energy += 1.000000 dp_bb(k,c) u1_aa(b,j) u2_abab(a,c,i,k) m2_aaaa(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_387;
        }

        if (include_u1_) {

            // energy += 1.000000 dp_bb(k,c) l2_abab(j,i,b,a) t2_bbbb(c,a,i,k) u1_aa(b,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_388;
        }

        if (include_u1_ && include_u2_) {

            // energy += -1.000000 dp_bb(k,c) u1_bb(b,j) u2_bbbb(c,a,i,k) m2_bbbb(i,j,b,a)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_389;
        }

        if (include_u1_) {

            // energy += -1.000000 dp_aa(k,c) l2_abab(j,i,b,a) t2_abab(c,a,k,i) u1_aa(b,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= scalar_390;

            // energy += 1.000000 dp_aa(k,c) l2_bbbb(i,j,b,a) t2_abab(c,a,k,i) u1_bb(b,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar_391;
        }

        if (include_u0_) {

            // energy += -1.000000 dp_aa(i,i) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy -= scalar_4 * u0;
        }

        {

            // rl2_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo tempPerm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rl2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_xaaaa_Lvvoo tempPerm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = 0.000000 * tempPerm_aaaa_vvoo("e,f,m,n");

            // rl2_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo tempPerm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rl2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rl2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_xbbbb_Lvvoo tempPerm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 0.000000 * tempPerm_bbbb_vvoo("e,f,m,n");
        }

        /*
            Total Number of Terms: 1994
            Cost Metric: (last) 6285 -> (new) 4644

            Total FLOP scaling:
            ------------------
               Scaling :  initial |  reorder | optimize |     diff
              -------- : -------- | -------- | -------- | --------
                o4v4L1 :      300 |        0 |        0 |        0
              -------- : -------- | -------- | -------- | --------
                o3v4L1 :       56 |        0 |        0 |        0
                o4v3L1 :       56 |        0 |        0 |        0
              -------- : -------- | -------- | -------- | --------
                o2v4L1 :       34 |        8 |        6 |       -2
                o3v3L1 :      490 |       88 |       88 |        0
                  o3v4 :      132 |        0 |        0 |        0
                o4v2L1 :       64 |       76 |       47 |      -29
                  o4v3 :      132 |        0 |        0 |        0
              -------- : -------- | -------- | -------- | --------
                o2v3L1 :      254 |      282 |      150 |     -132
                  o2v4 :      110 |       50 |       27 |      -23
                o3v2L1 :      254 |      562 |      292 |     -270
                  o3v3 :      230 |      246 |       99 |     -147
                  o4v2 :      110 |      216 |       47 |     -169
              -------- : -------- | -------- | -------- | --------
                o1v3L1 :        0 |       36 |       24 |      -12
                  o1v4 :       20 |        6 |        6 |        0
                o2v2L1 :     1501 |     1601 |      978 |     -623
                  o2v3 :      246 |      274 |       92 |     -182
                o3v1L1 :       12 |       52 |       42 |      -10
                  o3v2 :      246 |      470 |      171 |     -299
                  o4v1 :       20 |       30 |       16 |      -14
              -------- : -------- | -------- | -------- | --------
                  o0v4 :        0 |        0 |        1 |        1
                o1v2L1 :       36 |       50 |       32 |      -18
                  o1v3 :       20 |       24 |       48 |       24
                o2v1L1 :       36 |      154 |      112 |      -42
                  o2v2 :       88 |      114 |      357 |      243
                  o3v1 :       20 |       24 |       93 |       69
                  o4v0 :        0 |        3 |       13 |       10
              -------- : -------- | -------- | -------- | --------
                o0v2L1 :        0 |        0 |       12 |       12
                o1v1L1 :      804 |      880 |      662 |     -218
                  o1v2 :       48 |       16 |        8 |       -8
                o2v0L1 :        2 |        0 |       14 |       14
                  o2v1 :       52 |       80 |       29 |      -51
              -------- : -------- | -------- | -------- | --------
                  o0v2 :        0 |        0 |       12 |       12
                  o1v1 :       34 |        4 |       72 |       68
                  o2v0 :        0 |        4 |       20 |       16
              -------- : -------- | -------- | -------- | --------
                o0v0L1 :       56 |       56 |      207 |      151
              -------- : -------- | -------- | -------- | --------
                  o0v0 :      851 |      879 |      867 |      -12

            Total MEM scaling:
            ------------------
               Scaling :  initial |     last |      new |     diff
              -------- : -------- | -------- | -------- | --------
                o4v4L1 :      150 |        0 |        0 |        0
              -------- : -------- | -------- | -------- | --------
                o3v3L1 :      188 |        0 |        0 |        0
                o4v2L1 :       18 |        0 |        0 |        0
              -------- : -------- | -------- | -------- | --------
                  o2v4 :       66 |        0 |        0 |        0
                  o3v3 :      108 |        0 |        0 |        0
                  o4v2 :       66 |        0 |        0 |        0
              -------- : -------- | -------- | -------- | --------
                o2v2L1 :     1153 |      959 |      927 |      -32
                o3v1L1 :        8 |       60 |       16 |      -44
                o4v0L1 :        0 |       18 |        9 |       -9
              -------- : -------- | -------- | -------- | --------
                  o0v4 :       44 |        0 |        2 |        2
                  o1v3 :      170 |       86 |       72 |      -14
                  o2v2 :      490 |      614 |      574 |      -40
                  o3v1 :      170 |      244 |      169 |      -75
                  o4v0 :       44 |       96 |       26 |      -70
              -------- : -------- | -------- | -------- | --------
                o0v2L1 :        0 |       76 |       24 |      -52
                o1v1L1 :     1378 |     1466 |     1070 |     -396
                o2v0L1 :       18 |      176 |       38 |     -138
              -------- : -------- | -------- | -------- | --------
                  o0v2 :      100 |      132 |       40 |      -92
                  o1v1 :      102 |      178 |      151 |      -27
                  o2v0 :      102 |      194 |       60 |     -134
              -------- : -------- | -------- | -------- | --------
                o0v0L1 :      718 |      766 |      582 |     -184
              -------- : -------- | -------- | -------- | --------
                  o0v0 :      897 |      896 |      884 |      -12

         * */

        double coherent_scalar;
        if ( options_.get_bool("QED_USE_RELAXED_ORBITALS")) {
            coherent_scalar = coupling_factor_z * cc_wfn_->e_dip_z_;
        } else {
            coherent_scalar = -coupling_factor_z * cc_wfn_->nuc_dip_z_;
        }

        energy += coherent_scalar * cenergy;
        rl1_aa_vo(idx_map_[2])     += coherent_scalar * crl1_aa_vo(idx_map_[2]);
        rl1_bb_vo(idx_map_[2])     += coherent_scalar * crl1_bb_vo(idx_map_[2]);
        rl2_aaaa_vvoo(idx_map_[4]) += coherent_scalar * crl2_aaaa_vvoo(idx_map_[4]);
        rl2_abab_vvoo(idx_map_[4]) += coherent_scalar * crl2_abab_vvoo(idx_map_[4]);
        rl2_bbbb_vvoo(idx_map_[4]) += coherent_scalar * crl2_bbbb_vvoo(idx_map_[4]);
        if (include_u0_) rm0 = rm0 + coherent_scalar * crm0;
        if (include_u1_) {
            rm1_aa_vo(idx_map_[2]) += coherent_scalar * crm1_aa_vo(idx_map_[2]);
            rm1_bb_vo(idx_map_[2]) += coherent_scalar * crm1_bb_vo(idx_map_[2]);
        }
        if (include_u2_) {
            rm2_aaaa_vvoo(idx_map_[4]) += coherent_scalar * crm2_aaaa_vvoo(idx_map_[4]);
            rm2_abab_vvoo(idx_map_[4]) += coherent_scalar * crm2_abab_vvoo(idx_map_[4]);
            rm2_bbbb_vvoo(idx_map_[4]) += coherent_scalar * crm2_bbbb_vvoo(idx_map_[4]);
        }

        L_scalar_resids_["m0"] = rm0;

    }
    world_.gop.fence();
    return energy + cc_wfn_->enuc_ + cc_wfn_->average_electric_dipole_self_energy_;
}
