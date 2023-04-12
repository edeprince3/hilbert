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

#include "../../include/derived/eom_ea_driver.h"

void hilbert::EOM_EA_Driver::build_common_ops() {
    
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

    TA::TArrayD &t1_aa_vo =     cc_wfn_->amplitudes_["t1_aa"],
                &t1_bb_vo =     cc_wfn_->amplitudes_["t1_bb"],
                &u1_aa_vo =     cc_wfn_->amplitudes_["u1_aa"],
                &u1_bb_vo =     cc_wfn_->amplitudes_["u1_bb"],
                &t2_aaaa_vvoo = cc_wfn_->amplitudes_["t2_aaaa"],
                &t2_abab_vvoo = cc_wfn_->amplitudes_["t2_abab"],
                &t2_bbbb_vvoo = cc_wfn_->amplitudes_["t2_bbbb"],
                &u2_aaaa_vvoo = cc_wfn_->amplitudes_["u2_aaaa"],
                &u2_abab_vvoo = cc_wfn_->amplitudes_["u2_abab"],
                &u2_bbbb_vvoo = cc_wfn_->amplitudes_["u2_bbbb"];
    
    TArrayMap &V_blks_ = cc_wfn_->V_blks_;
    TArrayMap &F_blks_ = cc_wfn_->F_blks_;

    sigmaOps = std::vector<TA::TArrayD>(159);

    {

        // sigmaOps[0] += 1.000000 V_blks_["baab_vovv"]("e,i,a,b") t2_abab_vvoo("a,f,i,n")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[0]("e,b,f,n") = V_blks_["baab_vovv"]("e,i,a,b") * t2_abab_vvoo("a,f,i,n");

        // sigmaOps[1] += 1.000000 t2_abab_vvoo("a,e,n,i") V_blks_["abab_vovv"]("f,i,a,b")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[1]("f,e,b,n") = t2_abab_vvoo("a,e,n,i") * V_blks_["abab_vovv"]("f,i,a,b");

        // sigmaOps[2] += 1.000000 V_blks_["bbbb_vovv"]("e,i,a,b") t2_bbbb_vvoo("a,f,n,i")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[2]("e,b,f,n") = V_blks_["bbbb_vovv"]("e,i,a,b") * t2_bbbb_vvoo("a,f,n,i");

        // sigmaOps[3] += 1.000000 V_blks_["aaaa_vovv"]("e,i,a,b") t2_aaaa_vvoo("a,f,n,i")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[3]("e,b,f,n") = V_blks_["aaaa_vovv"]("e,i,a,b") * t2_aaaa_vvoo("a,f,n,i");

        // sigmaOps[4] += 1.000000 V_blks_["abab_vovv"]("e,i,b,a") t2_abab_vvoo("f,a,n,i")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[4]("e,b,f,n") = V_blks_["abab_vovv"]("e,i,b,a") * t2_abab_vvoo("f,a,n,i");

        // sigmaOps[5] += 1.000000 V_blks_["baab_vovv"]("f,i,b,a") t2_abab_vvoo("e,a,i,n")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[5]("b,e,f,n") = V_blks_["baab_vovv"]("f,i,b,a") * t2_abab_vvoo("e,a,i,n");

        // sigmaOps[6] += 1.000000 V_blks_["abab_vovv"]("e,i,b,a") t2_bbbb_vvoo("a,f,n,i")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[6]("e,b,f,n") = V_blks_["abab_vovv"]("e,i,b,a") * t2_bbbb_vvoo("a,f,n,i");
    }

    if (include_u2_) {

        // sigmaOps[7] += 1.000000 u2_abab_vvoo("a,c,j,i") V_blks_["baab_vovv"]("b,j,e,c")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[7]("a,e,b,i") = u2_abab_vvoo("a,c,j,i") * V_blks_["baab_vovv"]("b,j,e,c");

        // sigmaOps[8] += 1.000000 V_blks_["baab_vovv"]("b,j,c,e") u2_abab_vvoo("c,a,j,i")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[8]("b,e,a,i") = V_blks_["baab_vovv"]("b,j,c,e") * u2_abab_vvoo("c,a,j,i");

        // sigmaOps[9] += 1.000000 u2_abab_vvoo("c,a,i,j") V_blks_["abab_vovv"]("b,j,c,e")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[9]("b,a,e,i") = u2_abab_vvoo("c,a,i,j") * V_blks_["abab_vovv"]("b,j,c,e");

        // sigmaOps[10] += 1.000000 V_blks_["bbbb_vovv"]("b,j,c,e") u2_bbbb_vvoo("c,a,i,j")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[10]("b,e,a,i") = V_blks_["bbbb_vovv"]("b,j,c,e") * u2_bbbb_vvoo("c,a,i,j");

        // sigmaOps[11] += 1.000000 u2_abab_vvoo("a,c,i,j") V_blks_["abab_vovv"]("b,j,e,c")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[11]("a,b,e,i") = u2_abab_vvoo("a,c,i,j") * V_blks_["abab_vovv"]("b,j,e,c");

        // sigmaOps[12] += 1.000000 V_blks_["aaaa_vovv"]("b,j,c,e") u2_aaaa_vvoo("c,a,i,j")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[12]("b,e,a,i") = V_blks_["aaaa_vovv"]("b,j,c,e") * u2_aaaa_vvoo("c,a,i,j");
    }

    {

        // sigmaOps[13] += 1.000000 V_blks_["aaaa_vovv"]("e,i,a,b") t2_abab_vvoo("a,f,i,n")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[13]("e,b,f,n") = V_blks_["aaaa_vovv"]("e,i,a,b") * t2_abab_vvoo("a,f,i,n");

        // sigmaOps[14] += 1.000000 V_blks_["bbbb_vovv"]("e,i,a,b") t2_abab_vvoo("f,a,n,i")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[14]("f,e,b,n") = V_blks_["bbbb_vovv"]("e,i,a,b") * t2_abab_vvoo("f,a,n,i");

        // sigmaOps[15] += 1.000000 t2_aaaa_vvoo("a,f,n,i") V_blks_["baab_vovv"]("e,i,a,b")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[15]("f,e,b,n") = t2_aaaa_vvoo("a,f,n,i") * V_blks_["baab_vovv"]("e,i,a,b");
    }

    if (include_u2_) {

        // sigmaOps[16] += 1.000000 V_blks_["baab_vovv"]("e,i,a,b") u2_aaaa_vvoo("a,f,n,i")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[16]("f,e,b,n") = V_blks_["baab_vovv"]("e,i,a,b") * u2_aaaa_vvoo("a,f,n,i");

        // sigmaOps[17] += 1.000000 V_blks_["abab_vovv"]("e,i,b,a") u2_bbbb_vvoo("a,f,n,i")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[17]("e,b,f,n") = V_blks_["abab_vovv"]("e,i,b,a") * u2_bbbb_vvoo("a,f,n,i");

        // sigmaOps[18] += 1.000000 u2_abab_vvoo("a,f,i,n") V_blks_["aaaa_vovv"]("e,i,a,b")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[18]("e,b,f,n") = u2_abab_vvoo("a,f,i,n") * V_blks_["aaaa_vovv"]("e,i,a,b");

        // sigmaOps[19] += 1.000000 u2_abab_vvoo("f,a,n,i") V_blks_["bbbb_vovv"]("e,i,a,b")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[19]("f,e,b,n") = u2_abab_vvoo("f,a,n,i") * V_blks_["bbbb_vovv"]("e,i,a,b");
    }

    {

        // sigmaOps[20] += 1.000000 t2_abab_vvoo("e,a,n,i") V_blks_["abab_oovv"]("j,i,b,a")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[20]("e,b,n,j") = t2_abab_vvoo("e,a,n,i") * V_blks_["abab_oovv"]("j,i,b,a");

        // sigmaOps[21] += 1.000000 V_blks_["aaaa_oovv"]("j,i,a,b") t2_aaaa_vvoo("a,e,n,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[21]("b,e,j,n") = V_blks_["aaaa_oovv"]("j,i,a,b") * t2_aaaa_vvoo("a,e,n,i");

        // sigmaOps[22] += 1.000000 V_blks_["aaaa_oovv"]("j,i,a,b") t2_abab_vvoo("a,f,i,n")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[22]("b,f,j,n") = V_blks_["aaaa_oovv"]("j,i,a,b") * t2_abab_vvoo("a,f,i,n");

        // sigmaOps[23] += 1.000000 t2_abab_vvoo("f,a,n,i") V_blks_["bbbb_oovv"]("j,i,a,b")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[23]("f,b,n,j") = t2_abab_vvoo("f,a,n,i") * V_blks_["bbbb_oovv"]("j,i,a,b");

        // sigmaOps[24] += 1.000000 t2_aaaa_vvoo("a,f,n,i") V_blks_["abab_oovv"]("i,j,a,b")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[24]("f,b,n,j") = t2_aaaa_vvoo("a,f,n,i") * V_blks_["abab_oovv"]("i,j,a,b");

        // sigmaOps[25] += 1.000000 V_blks_["bbbb_oovv"]("j,i,a,b") t2_bbbb_vvoo("a,e,n,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[25]("b,e,j,n") = V_blks_["bbbb_oovv"]("j,i,a,b") * t2_bbbb_vvoo("a,e,n,i");

        // sigmaOps[26] += 1.000000 V_blks_["abab_oovv"]("j,i,b,a") t2_bbbb_vvoo("a,f,n,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[26]("b,f,j,n") = V_blks_["abab_oovv"]("j,i,b,a") * t2_bbbb_vvoo("a,f,n,i");

        // sigmaOps[27] += 1.000000 t2_abab_vvoo("a,e,i,n") V_blks_["abab_oovv"]("i,j,a,b")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[27]("e,b,n,j") = t2_abab_vvoo("a,e,i,n") * V_blks_["abab_oovv"]("i,j,a,b");

        // sigmaOps[28] += 1.000000 V_blks_["bbbb_oovo"]("j,i,a,n") t2_bbbb_vvoo("e,f,j,i")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[28]("a,e,f,n") = V_blks_["bbbb_oovo"]("j,i,a,n") * t2_bbbb_vvoo("e,f,j,i");

        // sigmaOps[29] += 1.000000 V_blks_["aaaa_oovo"]("j,i,a,n") t2_aaaa_vvoo("e,f,j,i")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[29]("a,e,f,n") = V_blks_["aaaa_oovo"]("j,i,a,n") * t2_aaaa_vvoo("e,f,j,i");

        // sigmaOps[30] += 1.000000 t2_abab_vvoo("f,e,i,j") V_blks_["abba_oovo"]("i,j,a,n")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[30]("f,e,a,n") = t2_abab_vvoo("f,e,i,j") * V_blks_["abba_oovo"]("i,j,a,n");

        // sigmaOps[31] += 1.000000 t2_abab_vvoo("e,f,j,i") V_blks_["abab_oovo"]("j,i,a,n")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[31]("e,a,f,n") = t2_abab_vvoo("e,f,j,i") * V_blks_["abab_oovo"]("j,i,a,n");
    }

    if (include_u2_) {

        // sigmaOps[32] += 1.000000 V_blks_["bbbb_oovv"]("n,j,b,f") u2_bbbb_vvoo("b,a,i,j")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[32]("f,a,n,i") = V_blks_["bbbb_oovv"]("n,j,b,f") * u2_bbbb_vvoo("b,a,i,j");

        // sigmaOps[33] += 1.000000 u2_abab_vvoo("a,b,i,j") V_blks_["bbbb_oovv"]("n,j,b,f")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[33]("a,f,i,n") = u2_abab_vvoo("a,b,i,j") * V_blks_["bbbb_oovv"]("n,j,b,f");

        // sigmaOps[34] += 1.000000 u2_abab_vvoo("b,a,j,i") V_blks_["abab_oovv"]("j,n,b,f")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[34]("a,f,i,n") = u2_abab_vvoo("b,a,j,i") * V_blks_["abab_oovv"]("j,n,b,f");

        // sigmaOps[35] += 1.000000 u2_aaaa_vvoo("b,a,i,j") V_blks_["aaaa_oovv"]("n,j,b,e")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[35]("a,e,i,n") = u2_aaaa_vvoo("b,a,i,j") * V_blks_["aaaa_oovv"]("n,j,b,e");

        // sigmaOps[36] += 1.000000 u2_abab_vvoo("b,a,j,i") V_blks_["aaaa_oovv"]("n,j,b,e")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[36]("e,a,n,i") = u2_abab_vvoo("b,a,j,i") * V_blks_["aaaa_oovv"]("n,j,b,e");

        // sigmaOps[37] += 1.000000 u2_bbbb_vvoo("b,a,i,j") V_blks_["abab_oovv"]("n,j,e,b")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[37]("e,a,n,i") = u2_bbbb_vvoo("b,a,i,j") * V_blks_["abab_oovv"]("n,j,e,b");

        // sigmaOps[38] += 1.000000 u2_abab_vvoo("a,b,i,j") V_blks_["abab_oovv"]("n,j,e,b")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[38]("a,e,i,n") = u2_abab_vvoo("a,b,i,j") * V_blks_["abab_oovv"]("n,j,e,b");

        // sigmaOps[39] += 1.000000 V_blks_["abab_oovv"]("j,n,b,f") u2_aaaa_vvoo("b,a,i,j")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[39]("a,f,i,n") = V_blks_["abab_oovv"]("j,n,b,f") * u2_aaaa_vvoo("b,a,i,j");
    }

    if (include_u1_) {

        // sigmaOps[40] += 1.000000 V_blks_["abab_oovv"]("k,j,c,e") u1_aa_vo("c,i") t2_abab_vvoo("b,a,k,j")
        // flops: o3v3: 1, o3v2: 1, o1v3: 1 | mem: o1v3: 2, o3v1: 1,
        sigmaOps[40]("b,e,a,i") = V_blks_["abab_oovv"]("k,j,c,e") * u1_aa_vo("c,i") * t2_abab_vvoo("b,a,k,j");

        // sigmaOps[41] += 1.000000 V_blks_["bbbb_oovv"]("k,j,c,e") u1_bb_vo("c,i") t2_bbbb_vvoo("b,a,k,j")
        // flops: o3v3: 1, o3v2: 1, o1v3: 1 | mem: o1v3: 2, o3v1: 1,
        sigmaOps[41]("e,b,a,i") = V_blks_["bbbb_oovv"]("k,j,c,e") * u1_bb_vo("c,i") * t2_bbbb_vvoo("b,a,k,j");

        // sigmaOps[42] += 1.000000 V_blks_["aaaa_oovv"]("k,j,c,e") u1_aa_vo("c,i") t2_aaaa_vvoo("b,a,k,j")
        // flops: o3v3: 1, o3v2: 1, o1v3: 1 | mem: o1v3: 2, o3v1: 1,
        sigmaOps[42]("e,b,a,i") = V_blks_["aaaa_oovv"]("k,j,c,e") * u1_aa_vo("c,i") * t2_aaaa_vvoo("b,a,k,j");

        // sigmaOps[43] += 1.000000 V_blks_["abab_oovv"]("k,j,e,c") u1_bb_vo("c,i") t2_abab_vvoo("a,b,k,j")
        // flops: o3v3: 1, o3v2: 1, o1v3: 1 | mem: o1v3: 2, o3v1: 1,
        sigmaOps[43]("e,a,b,i") = V_blks_["abab_oovv"]("k,j,e,c") * u1_bb_vo("c,i") * t2_abab_vvoo("a,b,k,j");
    }

    if (include_u2_) {

        // sigmaOps[44] += 1.000000 u2_bbbb_vvoo("b,a,k,j") V_blks_["bbbb_oovo"]("k,j,e,i")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[44]("b,a,e,i") = u2_bbbb_vvoo("b,a,k,j") * V_blks_["bbbb_oovo"]("k,j,e,i");

        // sigmaOps[45] += 1.000000 V_blks_["abab_oovo"]("j,k,e,i") u2_abab_vvoo("a,b,j,k")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[45]("e,a,b,i") = V_blks_["abab_oovo"]("j,k,e,i") * u2_abab_vvoo("a,b,j,k");

        // sigmaOps[46] += 1.000000 V_blks_["abba_oovo"]("k,j,e,i") u2_abab_vvoo("b,a,k,j")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[46]("b,e,a,i") = V_blks_["abba_oovo"]("k,j,e,i") * u2_abab_vvoo("b,a,k,j");

        // sigmaOps[47] += 1.000000 V_blks_["aaaa_oovo"]("k,j,e,i") u2_aaaa_vvoo("b,a,k,j")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[47]("e,b,a,i") = V_blks_["aaaa_oovo"]("k,j,e,i") * u2_aaaa_vvoo("b,a,k,j");
    }

    {

        // sigmaOps[48] += 1.000000 t2_abab_vvoo("a,e,n,i") V_blks_["abab_oovv"]("j,i,a,b")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[48]("e,b,n,j") = t2_abab_vvoo("a,e,n,i") * V_blks_["abab_oovv"]("j,i,a,b");

        // sigmaOps[49] += 1.000000 t2_abab_vvoo("e,a,i,n") V_blks_["abab_oovv"]("i,j,b,a")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[49]("e,b,n,j") = t2_abab_vvoo("e,a,i,n") * V_blks_["abab_oovv"]("i,j,b,a");
    }

    if (include_u1_) {

        // sigmaOps[50] += 1.000000 V_blks_["bbbb_vvvv"]("b,a,c,e") u1_bb_vo("c,i")
        // flops: o1v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[50]("b,a,e,i") = V_blks_["bbbb_vvvv"]("b,a,c,e") * u1_bb_vo("c,i");

        // sigmaOps[51] += 1.000000 V_blks_["aaaa_vvvv"]("b,a,c,e") u1_aa_vo("c,i")
        // flops: o1v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[51]("b,a,e,i") = V_blks_["aaaa_vvvv"]("b,a,c,e") * u1_aa_vo("c,i");

        // sigmaOps[52] += 1.000000 V_blks_["abab_vvvv"]("b,a,c,e") u1_aa_vo("c,i")
        // flops: o1v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[52]("b,a,e,i") = V_blks_["abab_vvvv"]("b,a,c,e") * u1_aa_vo("c,i");

        // sigmaOps[53] += 1.000000 V_blks_["abab_vvvv"]("a,b,e,c") u1_bb_vo("c,i")
        // flops: o1v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[53]("a,e,b,i") = V_blks_["abab_vvvv"]("a,b,e,c") * u1_bb_vo("c,i");

        // sigmaOps[54] += 1.000000 u1_bb_vo("c,k") V_blks_["abab_oovv"]("j,k,b,c")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[54]("b,j") = u1_bb_vo("c,k") * V_blks_["abab_oovv"]("j,k,b,c");

        // sigmaOps[55] += 1.000000 V_blks_["abab_oovv"]("k,j,c,b") u1_aa_vo("c,k")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[55]("b,j") = V_blks_["abab_oovv"]("k,j,c,b") * u1_aa_vo("c,k");

        // sigmaOps[56] += 1.000000 u1_aa_vo("c,k") V_blks_["aaaa_oovv"]("k,j,c,e")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[56]("e,j") = u1_aa_vo("c,k") * V_blks_["aaaa_oovv"]("k,j,c,e");

        // sigmaOps[57] += 1.000000 u1_bb_vo("c,k") V_blks_["bbbb_oovv"]("k,j,c,e")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[57]("e,j") = u1_bb_vo("c,k") * V_blks_["bbbb_oovv"]("k,j,c,e");
    }

    {

        // sigmaOps[58] += 1.000000 V_blks_["abab_oovv"]("j,i,a,b") t2_abab_vvoo("a,e,j,i")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sigmaOps[58]("b,e") = V_blks_["abab_oovv"]("j,i,a,b") * t2_abab_vvoo("a,e,j,i");

        // sigmaOps[59] += 1.000000 V_blks_["abab_oovv"]("j,i,b,a") t2_abab_vvoo("e,a,j,i")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sigmaOps[59]("b,e") = V_blks_["abab_oovv"]("j,i,b,a") * t2_abab_vvoo("e,a,j,i");

        // sigmaOps[60] += 1.000000 V_blks_["aaaa_oovv"]("j,i,a,b") t2_aaaa_vvoo("a,e,j,i")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sigmaOps[60]("b,e") = V_blks_["aaaa_oovv"]("j,i,a,b") * t2_aaaa_vvoo("a,e,j,i");

        // sigmaOps[61] += 1.000000 t2_bbbb_vvoo("a,e,j,i") V_blks_["bbbb_oovv"]("j,i,a,b")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sigmaOps[61]("e,b") = t2_bbbb_vvoo("a,e,j,i") * V_blks_["bbbb_oovv"]("j,i,a,b");

        // sigmaOps[62] += 1.000000 t2_abab_vvoo("a,b,n,i") V_blks_["abab_vovv"]("e,i,a,b")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[62]("e,n") = t2_abab_vvoo("a,b,n,i") * V_blks_["abab_vovv"]("e,i,a,b");

        // sigmaOps[63] += 1.000000 V_blks_["baab_vovv"]("f,i,a,b") t2_abab_vvoo("a,b,i,n")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[63]("f,n") = V_blks_["baab_vovv"]("f,i,a,b") * t2_abab_vvoo("a,b,i,n");

        // sigmaOps[64] += 1.000000 V_blks_["bbbb_vovv"]("f,i,a,b") t2_bbbb_vvoo("a,b,n,i")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[64]("f,n") = V_blks_["bbbb_vovv"]("f,i,a,b") * t2_bbbb_vvoo("a,b,n,i");

        // sigmaOps[65] += 1.000000 t2_aaaa_vvoo("a,b,n,i") V_blks_["aaaa_vovv"]("e,i,a,b")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[65]("e,n") = t2_aaaa_vvoo("a,b,n,i") * V_blks_["aaaa_vovv"]("e,i,a,b");
    }

    if (include_u2_) {

        // sigmaOps[66] += 1.000000 u2_abab_vvoo("a,b,j,i") V_blks_["abab_oovv"]("j,i,e,b")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sigmaOps[66]("a,e") = u2_abab_vvoo("a,b,j,i") * V_blks_["abab_oovv"]("j,i,e,b");

        // sigmaOps[67] += 1.000000 V_blks_["bbbb_oovv"]("j,i,b,e") u2_bbbb_vvoo("b,a,j,i")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sigmaOps[67]("e,a") = V_blks_["bbbb_oovv"]("j,i,b,e") * u2_bbbb_vvoo("b,a,j,i");

        // sigmaOps[68] += 1.000000 V_blks_["abab_oovv"]("i,j,b,e") u2_abab_vvoo("b,a,i,j")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sigmaOps[68]("e,a") = V_blks_["abab_oovv"]("i,j,b,e") * u2_abab_vvoo("b,a,i,j");

        // sigmaOps[69] += 1.000000 u2_aaaa_vvoo("b,a,j,i") V_blks_["aaaa_oovv"]("j,i,b,e")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sigmaOps[69]("a,e") = u2_aaaa_vvoo("b,a,j,i") * V_blks_["aaaa_oovv"]("j,i,b,e");

        // sigmaOps[70] += 1.000000 V_blks_["abab_vovv"]("a,j,b,c") u2_abab_vvoo("b,c,i,j")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[70]("a,i") = V_blks_["abab_vovv"]("a,j,b,c") * u2_abab_vvoo("b,c,i,j");

        // sigmaOps[71] += 1.000000 V_blks_["baab_vovv"]("a,j,b,c") u2_abab_vvoo("b,c,j,i")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[71]("a,i") = V_blks_["baab_vovv"]("a,j,b,c") * u2_abab_vvoo("b,c,j,i");

        // sigmaOps[72] += 1.000000 u2_aaaa_vvoo("b,c,i,j") V_blks_["aaaa_vovv"]("a,j,b,c")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[72]("a,i") = u2_aaaa_vvoo("b,c,i,j") * V_blks_["aaaa_vovv"]("a,j,b,c");

        // sigmaOps[73] += 1.000000 u2_bbbb_vvoo("b,c,i,j") V_blks_["bbbb_vovv"]("a,j,b,c")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[73]("a,i") = u2_bbbb_vvoo("b,c,i,j") * V_blks_["bbbb_vovv"]("a,j,b,c");
    }

    if (include_u1_) {

        // sigmaOps[74] += 1.000000 u1_aa_vo("b,i") V_blks_["abab_vovv"]("a,n,b,f")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[74]("a,f,i,n") = u1_aa_vo("b,i") * V_blks_["abab_vovv"]("a,n,b,f");

        // sigmaOps[75] += 1.000000 u1_bb_vo("b,i") V_blks_["bbbb_vovv"]("a,n,b,f")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[75]("a,f,i,n") = u1_bb_vo("b,i") * V_blks_["bbbb_vovv"]("a,n,b,f");

        // sigmaOps[76] += 1.000000 u1_aa_vo("b,i") V_blks_["aaaa_vovv"]("a,n,b,e")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[76]("a,e,i,n") = u1_aa_vo("b,i") * V_blks_["aaaa_vovv"]("a,n,b,e");

        // sigmaOps[77] += 1.000000 u1_bb_vo("b,i") V_blks_["baab_vovv"]("a,n,e,b")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[77]("e,a,n,i") = u1_bb_vo("b,i") * V_blks_["baab_vovv"]("a,n,e,b");
    }

    {

        // sigmaOps[78] += 1.000000 t2_aaaa_vvoo("a,e,j,i") V_blks_["aaaa_oovo"]("j,i,a,n")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[78]("e,n") = t2_aaaa_vvoo("a,e,j,i") * V_blks_["aaaa_oovo"]("j,i,a,n");

        // sigmaOps[79] += 1.000000 V_blks_["abab_oovv"]("i,j,b,a") t2_abab_vvoo("b,a,i,n")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[79]("j,n") = V_blks_["abab_oovv"]("i,j,b,a") * t2_abab_vvoo("b,a,i,n");
    }

    if (include_u1_) {

        // sigmaOps[80] += 1.000000 V_blks_["aaaa_oovv"]("k,j,b,c") u1_aa_vo("c,i") t2_aaaa_vvoo("b,a,k,j")
        // flops: o3v2: 2, o1v1: 1 | mem: o3v1: 1, o1v1: 2,
        sigmaOps[80]("a,i") = V_blks_["aaaa_oovv"]("k,j,b,c") * u1_aa_vo("c,i") * t2_aaaa_vvoo("b,a,k,j");
    }

    {

        // sigmaOps[81] += 1.000000 V_blks_["bbbb_oovv"]("j,i,a,b") t2_bbbb_vvoo("a,b,n,i")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[81]("j,n") = V_blks_["bbbb_oovv"]("j,i,a,b") * t2_bbbb_vvoo("a,b,n,i");

        // sigmaOps[82] += 1.000000 V_blks_["aaaa_oovv"]("j,i,a,b") t2_aaaa_vvoo("a,b,n,i")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[82]("j,n") = V_blks_["aaaa_oovv"]("j,i,a,b") * t2_aaaa_vvoo("a,b,n,i");

        // sigmaOps[83] += 1.000000 t2_abab_vvoo("e,a,j,i") V_blks_["abba_oovo"]("j,i,a,n")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[83]("e,n") = t2_abab_vvoo("e,a,j,i") * V_blks_["abba_oovo"]("j,i,a,n");

        // sigmaOps[84] += 1.000000 t2_abab_vvoo("a,b,n,i") V_blks_["abab_oovv"]("j,i,a,b")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[84]("n,j") = t2_abab_vvoo("a,b,n,i") * V_blks_["abab_oovv"]("j,i,a,b");
    }

    if (include_u1_) {

        // sigmaOps[85] += 1.000000 V_blks_["abab_oovv"]("k,j,b,c") u1_bb_vo("c,i") t2_abab_vvoo("b,a,k,j")
        // flops: o3v2: 2, o1v1: 1 | mem: o3v1: 1, o1v1: 2,
        sigmaOps[85]("a,i") = V_blks_["abab_oovv"]("k,j,b,c") * u1_bb_vo("c,i") * t2_abab_vvoo("b,a,k,j");
    }

    {

        // sigmaOps[86] += 1.000000 V_blks_["bbbb_oovo"]("j,i,a,n") t2_bbbb_vvoo("a,f,j,i")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[86]("f,n") = V_blks_["bbbb_oovo"]("j,i,a,n") * t2_bbbb_vvoo("a,f,j,i");
    }

    if (include_u1_) {

        // sigmaOps[87] += 1.000000 V_blks_["abab_oovv"]("k,j,c,b") u1_aa_vo("c,i") t2_abab_vvoo("a,b,k,j")
        // flops: o3v2: 2, o1v1: 1 | mem: o3v1: 1, o1v1: 2,
        sigmaOps[87]("a,i") = V_blks_["abab_oovv"]("k,j,c,b") * u1_aa_vo("c,i") * t2_abab_vvoo("a,b,k,j");
    }

    {

        // sigmaOps[88] += 1.000000 V_blks_["abab_oovo"]("j,i,a,n") t2_abab_vvoo("a,f,j,i")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[88]("f,n") = V_blks_["abab_oovo"]("j,i,a,n") * t2_abab_vvoo("a,f,j,i");
    }

    if (include_u1_) {

        // sigmaOps[89] += 1.000000 V_blks_["bbbb_oovv"]("k,j,b,c") u1_bb_vo("c,i") t2_bbbb_vvoo("b,a,k,j")
        // flops: o3v2: 2, o1v1: 1 | mem: o3v1: 1, o1v1: 2,
        sigmaOps[89]("a,i") = V_blks_["bbbb_oovv"]("k,j,b,c") * u1_bb_vo("c,i") * t2_bbbb_vvoo("b,a,k,j");
    }

    if (include_u2_) {

        // sigmaOps[90] += 1.000000 V_blks_["abba_oovo"]("j,k,b,i") u2_abab_vvoo("a,b,j,k")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[90]("a,i") = V_blks_["abba_oovo"]("j,k,b,i") * u2_abab_vvoo("a,b,j,k");

        // sigmaOps[91] += 1.000000 u2_aaaa_vvoo("b,a,k,j") V_blks_["aaaa_oovo"]("k,j,b,i")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[91]("a,i") = u2_aaaa_vvoo("b,a,k,j") * V_blks_["aaaa_oovo"]("k,j,b,i");

        // sigmaOps[92] += 1.000000 u2_abab_vvoo("b,a,k,j") V_blks_["abab_oovo"]("k,j,b,i")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[92]("a,i") = u2_abab_vvoo("b,a,k,j") * V_blks_["abab_oovo"]("k,j,b,i");

        // sigmaOps[93] += 1.000000 V_blks_["bbbb_oovo"]("k,j,b,i") u2_bbbb_vvoo("b,a,k,j")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[93]("a,i") = V_blks_["bbbb_oovo"]("k,j,b,i") * u2_bbbb_vvoo("b,a,k,j");

        // sigmaOps[94] += 1.000000 u2_abab_vvoo("a,b,j,i") V_blks_["abab_oovv"]("j,n,a,b")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[94]("i,n") = u2_abab_vvoo("a,b,j,i") * V_blks_["abab_oovv"]("j,n,a,b");

        // sigmaOps[95] += 1.000000 V_blks_["abab_oovv"]("n,j,a,b") u2_abab_vvoo("a,b,i,j")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[95]("n,i") = V_blks_["abab_oovv"]("n,j,a,b") * u2_abab_vvoo("a,b,i,j");

        // sigmaOps[96] += 1.000000 V_blks_["aaaa_oovv"]("n,j,a,b") u2_aaaa_vvoo("a,b,i,j")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[96]("n,i") = V_blks_["aaaa_oovv"]("n,j,a,b") * u2_aaaa_vvoo("a,b,i,j");

        // sigmaOps[97] += 1.000000 V_blks_["bbbb_oovv"]("n,j,a,b") u2_bbbb_vvoo("a,b,i,j")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[97]("n,i") = V_blks_["bbbb_oovv"]("n,j,a,b") * u2_bbbb_vvoo("a,b,i,j");
    }

    {

        // sigmaOps[98] += 1.000000 dp_aa_ov("i,a") t2_abab_vvoo("a,f,i,n")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[98]("f,n") = dp_aa_ov("i,a") * t2_abab_vvoo("a,f,i,n");

        // sigmaOps[99] += 1.000000 dp_aa_ov("i,a") t2_aaaa_vvoo("a,e,n,i")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[99]("e,n") = dp_aa_ov("i,a") * t2_aaaa_vvoo("a,e,n,i");

        // sigmaOps[100] += 1.000000 dp_bb_ov("i,a") t2_bbbb_vvoo("a,f,n,i")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[100]("f,n") = dp_bb_ov("i,a") * t2_bbbb_vvoo("a,f,n,i");

        // sigmaOps[101] += 1.000000 dp_bb_ov("i,a") t2_abab_vvoo("e,a,n,i")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[101]("e,n") = dp_bb_ov("i,a") * t2_abab_vvoo("e,a,n,i");
    }

    if (include_u2_) {

        // sigmaOps[102] += 1.000000 dp_bb_ov("i,a") u2_bbbb_vvoo("a,f,n,i")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[102]("f,n") = dp_bb_ov("i,a") * u2_bbbb_vvoo("a,f,n,i");

        // sigmaOps[103] += 1.000000 dp_aa_ov("i,a") u2_abab_vvoo("a,f,i,n")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[103]("f,n") = dp_aa_ov("i,a") * u2_abab_vvoo("a,f,i,n");

        // sigmaOps[104] += 1.000000 u2_abab_vvoo("e,a,n,i") dp_bb_ov("i,a")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[104]("e,n") = u2_abab_vvoo("e,a,n,i") * dp_bb_ov("i,a");

        // sigmaOps[105] += 1.000000 dp_aa_ov("i,a") u2_aaaa_vvoo("a,e,n,i")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[105]("e,n") = dp_aa_ov("i,a") * u2_aaaa_vvoo("a,e,n,i");
    }

    if (include_u1_) {

        // sigmaOps[106] += 1.000000 u1_bb_vo("b,i") V_blks_["bbbb_vovv"]("a,i,b,e")
        // flops: o1v3: 1, o0v2: 1 | mem: o0v2: 2,
        sigmaOps[106]("a,e") = u1_bb_vo("b,i") * V_blks_["bbbb_vovv"]("a,i,b,e");

        // sigmaOps[107] += 1.000000 u1_aa_vo("b,i") V_blks_["baab_vovv"]("a,i,b,e")
        // flops: o1v3: 1, o0v2: 1 | mem: o0v2: 2,
        sigmaOps[107]("a,e") = u1_aa_vo("b,i") * V_blks_["baab_vovv"]("a,i,b,e");

        // sigmaOps[108] += 1.000000 V_blks_["aaaa_vovv"]("a,i,b,e") u1_aa_vo("b,i")
        // flops: o1v3: 1, o0v2: 1 | mem: o0v2: 2,
        sigmaOps[108]("a,e") = V_blks_["aaaa_vovv"]("a,i,b,e") * u1_aa_vo("b,i");

        // sigmaOps[109] += 1.000000 u1_bb_vo("b,i") V_blks_["abab_vovv"]("a,i,e,b")
        // flops: o1v3: 1, o0v2: 1 | mem: o0v2: 2,
        sigmaOps[109]("a,e") = u1_bb_vo("b,i") * V_blks_["abab_vovv"]("a,i,e,b");

        // sigmaOps[110] += 1.000000 dp_aa_ov("i,a") u1_aa_vo("a,n")
        // flops: o2v1: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[110]("i,n") = dp_aa_ov("i,a") * u1_aa_vo("a,n");

        // sigmaOps[111] += 1.000000 u1_bb_vo("a,n") dp_bb_ov("i,a")
        // flops: o2v1: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[111]("n,i") = u1_bb_vo("a,n") * dp_bb_ov("i,a");

        // sigmaOps[112] += 1.000000 dp_aa_vv("e,a") u1_aa_vo("a,n")
        // flops: o1v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[112]("e,n") = dp_aa_vv("e,a") * u1_aa_vo("a,n");

        // sigmaOps[113] += 1.000000 dp_bb_vv("f,a") u1_bb_vo("a,n")
        // flops: o1v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[113]("f,n") = dp_bb_vv("f,a") * u1_bb_vo("a,n");
    }

    {

        // sigmaOps[114] += 1.000000 F_blks_["bb_ov"]("i,a") t2_abab_vvoo("e,a,n,i")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[114]("e,n") = F_blks_["bb_ov"]("i,a") * t2_abab_vvoo("e,a,n,i");
    }

    if (include_u1_) {

        // sigmaOps[115] += 1.000000 V_blks_["bbbb_oovv"]("k,j,b,c") u1_bb_vo("c,k") t2_bbbb_vvoo("b,a,i,j")
        // flops: o2v2: 2, o1v1: 1 | mem: o1v1: 3,
        sigmaOps[115]("a,i") = V_blks_["bbbb_oovv"]("k,j,b,c") * u1_bb_vo("c,k") * t2_bbbb_vvoo("b,a,i,j");

        // sigmaOps[116] += 1.000000 u1_aa_vo("c,k") V_blks_["aaaa_oovv"]("k,j,b,c")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[116]("b,j") = u1_aa_vo("c,k") * V_blks_["aaaa_oovv"]("k,j,b,c");
    }

    {

        // sigmaOps[117] += 1.000000 t2_abab_vvoo("a,f,i,n") F_blks_["aa_ov"]("i,a")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[117]("f,n") = t2_abab_vvoo("a,f,i,n") * F_blks_["aa_ov"]("i,a");

        // sigmaOps[118] += 1.000000 t2_bbbb_vvoo("a,f,n,i") F_blks_["bb_ov"]("i,a")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[118]("f,n") = t2_bbbb_vvoo("a,f,n,i") * F_blks_["bb_ov"]("i,a");

        // sigmaOps[119] += 1.000000 F_blks_["aa_ov"]("i,a") t2_aaaa_vvoo("a,e,n,i")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[119]("e,n") = F_blks_["aa_ov"]("i,a") * t2_aaaa_vvoo("a,e,n,i");
    }

    if (include_u1_) {

        // sigmaOps[120] += 1.000000 V_blks_["bbbb_oovv"]("k,j,b,c") u1_bb_vo("c,k") t2_abab_vvoo("a,b,i,j")
        // flops: o2v2: 2, o1v1: 1 | mem: o1v1: 3,
        sigmaOps[120]("a,i") = V_blks_["bbbb_oovv"]("k,j,b,c") * u1_bb_vo("c,k") * t2_abab_vvoo("a,b,i,j");
    }

    if (include_u2_) {

        // sigmaOps[121] += 1.000000 u2_aaaa_vvoo("b,a,i,j") F_blks_["aa_ov"]("j,b")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[121]("a,i") = u2_aaaa_vvoo("b,a,i,j") * F_blks_["aa_ov"]("j,b");
    }

    {

        // sigmaOps[122] += 1.000000 t2_abab_vvoo("a,b,i,j") dp_bb_ov("j,b")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[122]("a,i") = t2_abab_vvoo("a,b,i,j") * dp_bb_ov("j,b");
    }

    if (include_u1_) {

        // sigmaOps[123] += 1.000000 V_blks_["abba_vovo"]("a,j,b,i") u1_bb_vo("b,j")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[123]("a,i") = V_blks_["abba_vovo"]("a,j,b,i") * u1_bb_vo("b,j");

        // sigmaOps[124] += 1.000000 V_blks_["baab_vovo"]("a,j,b,i") u1_aa_vo("b,j")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[124]("a,i") = V_blks_["baab_vovo"]("a,j,b,i") * u1_aa_vo("b,j");

        // sigmaOps[125] += 1.000000 V_blks_["bbbb_vovo"]("a,j,b,i") u1_bb_vo("b,j")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[125]("a,i") = V_blks_["bbbb_vovo"]("a,j,b,i") * u1_bb_vo("b,j");

        // sigmaOps[126] += 1.000000 V_blks_["aaaa_vovo"]("a,j,b,i") u1_aa_vo("b,j")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[126]("a,i") = V_blks_["aaaa_vovo"]("a,j,b,i") * u1_aa_vo("b,j");
    }

    if (include_u2_) {

        // sigmaOps[127] += 1.000000 F_blks_["bb_ov"]("j,b") u2_abab_vvoo("a,b,i,j")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[127]("a,i") = F_blks_["bb_ov"]("j,b") * u2_abab_vvoo("a,b,i,j");

        // sigmaOps[128] += 1.000000 F_blks_["bb_ov"]("j,b") u2_bbbb_vvoo("b,a,i,j")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[128]("a,i") = F_blks_["bb_ov"]("j,b") * u2_bbbb_vvoo("b,a,i,j");
    }

    {

        // sigmaOps[129] += 1.000000 dp_aa_ov("j,b") t2_aaaa_vvoo("b,a,i,j")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[129]("a,i") = dp_aa_ov("j,b") * t2_aaaa_vvoo("b,a,i,j");

        // sigmaOps[130] += 1.000000 dp_aa_ov("j,b") t2_abab_vvoo("b,a,j,i")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[130]("a,i") = dp_aa_ov("j,b") * t2_abab_vvoo("b,a,j,i");

        // sigmaOps[131] += 1.000000 dp_bb_ov("j,b") t2_bbbb_vvoo("b,a,i,j")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[131]("a,i") = dp_bb_ov("j,b") * t2_bbbb_vvoo("b,a,i,j");
    }

    if (include_u2_) {

        // sigmaOps[132] += 1.000000 u2_abab_vvoo("b,a,j,i") F_blks_["aa_ov"]("j,b")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[132]("a,i") = u2_abab_vvoo("b,a,j,i") * F_blks_["aa_ov"]("j,b");
    }

    if (include_u1_) {

        // sigmaOps[133] += 1.000000 u1_aa_vo("a,i") V_blks_["aaaa_oovv"]("n,i,a,e")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[133]("e,n") = u1_aa_vo("a,i") * V_blks_["aaaa_oovv"]("n,i,a,e");

        // sigmaOps[134] += 1.000000 u1_bb_vo("a,i") V_blks_["bbbb_oovv"]("n,i,a,f")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[134]("f,n") = u1_bb_vo("a,i") * V_blks_["bbbb_oovv"]("n,i,a,f");

        // sigmaOps[135] += 1.000000 u1_bb_vo("f,i") dp_bb_oo("i,n")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[135]("f,n") = u1_bb_vo("f,i") * dp_bb_oo("i,n");

        // sigmaOps[136] += 1.000000 dp_aa_oo("i,n") u1_aa_vo("e,i")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[136]("e,n") = dp_aa_oo("i,n") * u1_aa_vo("e,i");

        // sigmaOps[137] += 1.000000 u1_aa_vo("a,j") V_blks_["abab_oovo"]("j,n,a,i")
        // flops: o3v1: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[137]("n,i") = u1_aa_vo("a,j") * V_blks_["abab_oovo"]("j,n,a,i");

        // sigmaOps[138] += 1.000000 u1_bb_vo("a,j") V_blks_["bbbb_oovo"]("n,j,a,i")
        // flops: o3v1: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[138]("n,i") = u1_bb_vo("a,j") * V_blks_["bbbb_oovo"]("n,j,a,i");

        // sigmaOps[139] += 1.000000 V_blks_["aaaa_oovo"]("n,j,a,i") u1_aa_vo("a,j")
        // flops: o3v1: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[139]("n,i") = V_blks_["aaaa_oovo"]("n,j,a,i") * u1_aa_vo("a,j");

        // sigmaOps[140] += 1.000000 u1_bb_vo("a,j") V_blks_["abba_oovo"]("n,j,a,i")
        // flops: o3v1: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[140]("n,i") = u1_bb_vo("a,j") * V_blks_["abba_oovo"]("n,j,a,i");

        // sigmaOps[141] += 1.000000 F_blks_["bb_vv"]("a,b") u1_bb_vo("b,i")
        // flops: o1v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[141]("a,i") = F_blks_["bb_vv"]("a,b") * u1_bb_vo("b,i");

        // sigmaOps[142] += 1.000000 u1_aa_vo("b,i") F_blks_["aa_vv"]("a,b")
        // flops: o1v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[142]("a,i") = u1_aa_vo("b,i") * F_blks_["aa_vv"]("a,b");

        // sigmaOps[143] += 1.000000 F_blks_["aa_oo"]("j,i") u1_aa_vo("a,j")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[143]("a,i") = F_blks_["aa_oo"]("j,i") * u1_aa_vo("a,j");

        // sigmaOps[144] += 1.000000 u1_bb_vo("a,j") F_blks_["bb_oo"]("j,i")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[144]("a,i") = u1_bb_vo("a,j") * F_blks_["bb_oo"]("j,i");

        // sigmaOps[145] += 1.000000 u1_bb_vo("a,i") F_blks_["bb_ov"]("n,a")
        // flops: o2v1: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[145]("i,n") = u1_bb_vo("a,i") * F_blks_["bb_ov"]("n,a");

        // sigmaOps[146] += 1.000000 u1_aa_vo("a,i") F_blks_["aa_ov"]("n,a")
        // flops: o2v1: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[146]("i,n") = u1_aa_vo("a,i") * F_blks_["aa_ov"]("n,a");

        // sigmaOps[151] += 1.000000 t2_aaaa_vvoo("b,a,i,j") tempOp55_aa_vo("b,j") 1 u1_bb_vo("c,k") V_blks_["abab_oovv"]("j,k,b,c")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[151]("a,i") = t2_aaaa_vvoo("b,a,i,j") * sigmaOps[54]("b,j");

        // sigmaOps[152] += 1.000000 t2_abab_vvoo("b,a,j,i") tempOp55_aa_vo("b,j") 1 u1_bb_vo("c,k") V_blks_["abab_oovv"]("j,k,b,c")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[152]("a,i") = t2_abab_vvoo("b,a,j,i") * sigmaOps[54]("b,j");

        // sigmaOps[147] += 1.000000 tempOp56_bb_vo("b,j") t2_abab_vvoo("a,b,i,j") 1 V_blks_["abab_oovv"]("k,j,c,b") u1_aa_vo("c,k")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[147]("a,i") = sigmaOps[55]("b,j") * t2_abab_vvoo("a,b,i,j");

        // sigmaOps[150] += 1.000000 tempOp56_bb_vo("b,j") t2_bbbb_vvoo("b,a,i,j") 1 V_blks_["abab_oovv"]("k,j,c,b") u1_aa_vo("c,k")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[150]("a,i") = sigmaOps[55]("b,j") * t2_bbbb_vvoo("b,a,i,j");

        // sigmaOps[153] += 1.000000 tempOp80_bb_oo("k,i") u1_bb_vo("a,k") 1 V_blks_["abab_oovv"]("i,j,b,a") t2_abab_vvoo("b,a,i,n")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[153]("a,i") = sigmaOps[79]("k,i") * u1_bb_vo("a,k");

        // sigmaOps[158] += 1.000000 tempOp82_bb_oo("k,i") u1_bb_vo("a,k") 1 V_blks_["bbbb_oovv"]("j,i,a,b") t2_bbbb_vvoo("a,b,n,i")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[158]("a,i") = sigmaOps[81]("k,i") * u1_bb_vo("a,k");

        // sigmaOps[156] += 1.000000 tempOp83_aa_oo("k,i") u1_aa_vo("a,k") 1 V_blks_["aaaa_oovv"]("j,i,a,b") t2_aaaa_vvoo("a,b,n,i")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[156]("a,i") = sigmaOps[82]("k,i") * u1_aa_vo("a,k");

        // sigmaOps[157] += 1.000000 u1_aa_vo("a,k") tempOp85_aa_oo("i,k") 1 t2_abab_vvoo("a,b,n,i") V_blks_["abab_oovv"]("j,i,a,b")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[157]("a,i") = u1_aa_vo("a,k") * sigmaOps[84]("i,k");

        // sigmaOps[154] += 1.000000 u1_aa_vo("a,j") tempOp111_aa_oo("j,i") 1 dp_aa_ov("i,a") u1_aa_vo("a,n")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[154]("a,i") = u1_aa_vo("a,j") * sigmaOps[110]("j,i");

        // sigmaOps[155] += 1.000000 u1_bb_vo("a,j") tempOp112_bb_oo("i,j") 1 u1_bb_vo("a,n") dp_bb_ov("i,a")
        // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[155]("a,i") = u1_bb_vo("a,j") * sigmaOps[111]("i,j");

        // sigmaOps[148] += 1.000000 t2_abab_vvoo("b,a,j,i") tempOp117_aa_vo("b,j") 1 u1_aa_vo("c,k") V_blks_["aaaa_oovv"]("k,j,b,c")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[148]("a,i") = t2_abab_vvoo("b,a,j,i") * sigmaOps[116]("b,j");

        // sigmaOps[149] += 1.000000 t2_aaaa_vvoo("b,a,i,j") tempOp117_aa_vo("b,j") 1 u1_aa_vo("c,k") V_blks_["aaaa_oovv"]("k,j,b,c")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[149]("a,i") = t2_aaaa_vvoo("b,a,i,j") * sigmaOps[116]("b,j");
    }

    world_.gop.fence();
}

void hilbert::EOM_EA_Driver::build_Hc_cH(size_t L) {

    // ensure that the integrals are transformed with t1 amplitudes
    if (!cc_wfn_->has_t1_integrals_) cc_wfn_->transform_integrals(true);

    /// build sigma vectors for this trial
    sigvec_blks_.clear();

    // reduce sigma into spins
    sigvec_blks_["sigmar1_a"] = HelperD::makeTensor(world_, {L, va_}, true);
    sigvec_blks_["sigmar1_b"] = HelperD::makeTensor(world_, {L, vb_}, true);
    sigvec_blks_["sigmal1_a"] = HelperD::makeTensor(world_, {L, va_}, true);
    sigvec_blks_["sigmal1_b"] = HelperD::makeTensor(world_, {L, vb_}, true);

    sigvec_blks_["sigmar2_aaa"] = HelperD::makeTensor(world_, {L, va_, va_, oa_}, true);
    sigvec_blks_["sigmar2_abb"] = HelperD::makeTensor(world_, {L, va_, vb_, ob_}, true);
    sigvec_blks_["sigmar2_baa"] = HelperD::makeTensor(world_, {L, vb_, va_, oa_}, true);
    sigvec_blks_["sigmar2_bbb"] = HelperD::makeTensor(world_, {L, vb_, vb_, ob_}, true);
    sigvec_blks_["sigmal2_aaa"] = HelperD::makeTensor(world_, {L, va_, va_, oa_}, true);
    sigvec_blks_["sigmal2_abb"] = HelperD::makeTensor(world_, {L, va_, vb_, ob_}, true);
    sigvec_blks_["sigmal2_baa"] = HelperD::makeTensor(world_, {L, vb_, va_, oa_}, true);
    sigvec_blks_["sigmal2_bbb"] = HelperD::makeTensor(world_, {L, vb_, vb_, ob_}, true);

    if (include_u1_) {
        sigvec_blks_["sigmas1_a"] = HelperD::makeTensor(world_, {L, va_}, true);
        sigvec_blks_["sigmas1_b"] = HelperD::makeTensor(world_, {L, vb_}, true);
        sigvec_blks_["sigmam1_a"] = HelperD::makeTensor(world_, {L, va_}, true);
        sigvec_blks_["sigmam1_b"] = HelperD::makeTensor(world_, {L, vb_}, true);
    }
    if (include_u2_) {
        sigvec_blks_["sigmas2_aaa"] = HelperD::makeTensor(world_, {L, va_, va_, oa_}, true);
        sigvec_blks_["sigmas2_abb"] = HelperD::makeTensor(world_, {L, va_, vb_, ob_}, true);
        sigvec_blks_["sigmas2_baa"] = HelperD::makeTensor(world_, {L, vb_, va_, oa_}, true);
        sigvec_blks_["sigmas2_bbb"] = HelperD::makeTensor(world_, {L, vb_, vb_, ob_}, true);
        sigvec_blks_["sigmam2_aaa"] = HelperD::makeTensor(world_, {L, va_, va_, oa_}, true);
        sigvec_blks_["sigmam2_abb"] = HelperD::makeTensor(world_, {L, va_, vb_, ob_}, true);
        sigvec_blks_["sigmam2_baa"] = HelperD::makeTensor(world_, {L, vb_, va_, oa_}, true);
        sigvec_blks_["sigmam2_bbb"] = HelperD::makeTensor(world_, {L, vb_, vb_, ob_}, true);
    }
    
    // get cavity information

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

    // get reference to electronic integrals
    TArrayMap &V_blks_ = cc_wfn_->V_blks_;
    TArrayMap &F_blks_ = cc_wfn_->F_blks_;
    TArrayMap &Id_blks_ = cc_wfn_->Id_blks_;

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
    
    double u0 = cc_wfn_->scalar_amps_["u0"];

    bool include_s0_ = include_u0_;
    bool include_m0_ = include_u0_;

    bool include_s1_ = include_u1_;
    bool include_m1_ = include_u1_;

    bool include_s2_ = include_u2_;
    bool include_m2_ = include_u2_;
    
    /// reference right operators
    auto &r1_xa_Lv = evec_blks_["r1_a"];
    auto &r1_xb_Lv = evec_blks_["r1_b"];
    auto &r2_xaaa_Lvvo = evec_blks_["r2_aaa"];
    auto &r2_xabb_Lvvo = evec_blks_["r2_abb"];
    auto &r2_xbaa_Lvvo = evec_blks_["r2_baa"];
    auto &r2_xbbb_Lvvo = evec_blks_["r2_bbb"];
            
    auto &s1_xa_Lv = evec_blks_["s1_a"];
    auto &s1_xb_Lv = evec_blks_["s1_b"];
    auto &s2_xaaa_Lvvo = evec_blks_["s2_aaa"];
    auto &s2_xabb_Lvvo = evec_blks_["s2_abb"];
    auto &s2_xbaa_Lvvo = evec_blks_["s2_baa"];
    auto &s2_xbbb_Lvvo = evec_blks_["s2_bbb"];

    /// reference left operators
    auto &l1_xa_Lv = evec_blks_["l1_a"];
    auto &l1_xb_Lv = evec_blks_["l1_b"];
    auto &l2_xaaa_Lovv = evec_blks_["l2_aaa"];
    auto &l2_xbba_Lovv = evec_blks_["l2_bba"];
    auto &l2_xaab_Lovv = evec_blks_["l2_aab"];
    auto &l2_xbbb_Lovv = evec_blks_["l2_bbb"];

    auto &m1_xa_Lv = evec_blks_["m1_a"];
    auto &m1_xb_Lv = evec_blks_["m1_b"];
    auto &m2_xaaa_Lovv = evec_blks_["m2_aaa"];
    auto &m2_xbba_Lovv = evec_blks_["m2_bba"];
    auto &m2_xaab_Lovv = evec_blks_["m2_aab"];
    auto &m2_xbbb_Lovv = evec_blks_["m2_bbb"];

    /// reference sigma vectors
    auto &sigmar1_xa_Lv = sigvec_blks_["sigmar1_a"];
    auto &sigmar1_xb_Lv = sigvec_blks_["sigmar1_b"];
    auto &sigmal1_xa_Lv = sigvec_blks_["sigmal1_a"];
    auto &sigmal1_xb_Lv = sigvec_blks_["sigmal1_b"];

    auto &sigmar2_xaaa_Lvvo = sigvec_blks_["sigmar2_aaa"];
    auto &sigmar2_xabb_Lvvo = sigvec_blks_["sigmar2_abb"];
    auto &sigmar2_xbaa_Lvvo = sigvec_blks_["sigmar2_baa"];
    auto &sigmar2_xbbb_Lvvo = sigvec_blks_["sigmar2_bbb"];
    auto &sigmal2_xaaa_Lvvo = sigvec_blks_["sigmal2_aaa"];
    auto &sigmal2_xabb_Lvvo = sigvec_blks_["sigmal2_abb"];
    auto &sigmal2_xbaa_Lvvo = sigvec_blks_["sigmal2_baa"];
    auto &sigmal2_xbbb_Lvvo = sigvec_blks_["sigmal2_bbb"];

    auto &sigmas1_xa_Lv = sigvec_blks_["sigmas1_a"];
    auto &sigmas1_xb_Lv = sigvec_blks_["sigmas1_b"];
    auto &sigmam1_xa_Lv = sigvec_blks_["sigmam1_a"];
    auto &sigmam1_xb_Lv = sigvec_blks_["sigmam1_b"];

    auto &sigmas2_xaaa_Lvvo = sigvec_blks_["sigmas2_aaa"];
    auto &sigmas2_xabb_Lvvo = sigvec_blks_["sigmas2_abb"];
    auto &sigmas2_xbaa_Lvvo = sigvec_blks_["sigmas2_baa"];
    auto &sigmas2_xbbb_Lvvo = sigvec_blks_["sigmas2_bbb"];
    auto &sigmam2_xaaa_Lvvo = sigvec_blks_["sigmam2_aaa"];
    auto &sigmam2_xabb_Lvvo = sigvec_blks_["sigmam2_abb"];
    auto &sigmam2_xbaa_Lvvo = sigvec_blks_["sigmam2_baa"];
    auto &sigmam2_xbbb_Lvvo = sigvec_blks_["sigmam2_bbb"];

    // coherent state basis
    auto csigmar1_xa_Lv = sigvec_blks_["sigmar1_a"].clone();
    auto csigmar1_xb_Lv = sigvec_blks_["sigmar1_b"].clone();
    auto csigmal1_xa_Lv = sigvec_blks_["sigmal1_a"].clone();
    auto csigmal1_xb_Lv = sigvec_blks_["sigmal1_b"].clone();

    auto csigmar2_xaaa_Lvvo = sigvec_blks_["sigmar2_aaa"].clone();
    auto csigmar2_xabb_Lvvo = sigvec_blks_["sigmar2_abb"].clone();
    auto csigmar2_xbaa_Lvvo = sigvec_blks_["sigmar2_baa"].clone();
    auto csigmar2_xbbb_Lvvo = sigvec_blks_["sigmar2_bbb"].clone();
    auto csigmal2_xaaa_Lvvo = sigvec_blks_["sigmal2_aaa"].clone();
    auto csigmal2_xabb_Lvvo = sigvec_blks_["sigmal2_abb"].clone();
    auto csigmal2_xbaa_Lvvo = sigvec_blks_["sigmal2_baa"].clone();
    auto csigmal2_xbbb_Lvvo = sigvec_blks_["sigmal2_bbb"].clone();

    auto csigmas1_xa_Lv = sigvec_blks_["sigmas1_a"].clone();
    auto csigmas1_xb_Lv = sigvec_blks_["sigmas1_b"].clone();
    auto csigmam1_xa_Lv = sigvec_blks_["sigmam1_a"].clone();
    auto csigmam1_xb_Lv = sigvec_blks_["sigmam1_b"].clone();

    auto csigmas2_xaaa_Lvvo = sigvec_blks_["sigmas2_aaa"].clone();
    auto csigmas2_xabb_Lvvo = sigvec_blks_["sigmas2_abb"].clone();
    auto csigmas2_xbaa_Lvvo = sigvec_blks_["sigmas2_baa"].clone();
    auto csigmas2_xbbb_Lvvo = sigvec_blks_["sigmas2_bbb"].clone();
    auto csigmam2_xaaa_Lvvo = sigvec_blks_["sigmam2_aaa"].clone();
    auto csigmam2_xabb_Lvvo = sigvec_blks_["sigmam2_abb"].clone();
    auto csigmam2_xbaa_Lvvo = sigvec_blks_["sigmam2_baa"].clone();
    auto csigmam2_xbbb_Lvvo = sigvec_blks_["sigmam2_bbb"].clone();

    {
        TA::TArrayD tempPerm_xaaa_Lvvo;
        TA::TArrayD tempPerm_xbbb_Lvvo;
        double scalar0;
        double scalar1;
        double scalar2;
        double scalar3;
        double scalar4;
        double scalar5;
        double scalar6;
        double scalar7;
        double scalar8;
        double scalar9;
        double scalar10;
        double scalar11;
        double scalar12;
        double scalar13;
        double scalar14;
        double scalar15;
        double scalar16;
        double scalar17;
        double scalar18;
        scalar0 = dot(dp_aa_ov("i,a"), u1_aa_vo("a,i"));
        scalar1 = dot(dp_bb_ov("i,a"), u1_bb_vo("a,i"));
        scalar2 = dot(V_blks_["aaaa_oovv"]("j,i,a,b"), t2_aaaa_vvoo("a,b,j,i"));
        scalar3 = dot(V_blks_["abab_oovv"]("j,i,a,b"), t2_abab_vvoo("a,b,j,i"));
        scalar4 = dot(V_blks_["bbbb_oovv"]("j,i,a,b"), t2_bbbb_vvoo("a,b,j,i"));
        scalar5 = dot(F_blks_["aa_oo"]("o0,o1"), Id_blks_["aa_oo"]("o0,o1"));
        scalar6 = dot(F_blks_["bb_oo"]("o0,o1"), Id_blks_["bb_oo"]("o0,o1"));
        scalar7 = dot(dp_aa_oo("o0,o1"), Id_blks_["aa_oo"]("o0,o1"));
        scalar8 = dot(dp_bb_oo("o0,o1"), Id_blks_["bb_oo"]("o0,o1"));
        scalar9 = dot(V_blks_["aaaa_oooo"]("o0,o1,o2,o3"), Id_blks_["aaaa_oooo"]("o0,o1,o2,o3"));
        scalar10 = dot(V_blks_["abab_oooo"]("o0,o1,o2,o3"), Id_blks_["abab_oooo"]("o0,o1,o2,o3"));
        scalar11 = dot(V_blks_["bbbb_oooo"]("o0,o1,o2,o3"), Id_blks_["bbbb_oooo"]("o0,o1,o2,o3"));
        scalar12 = dot(F_blks_["aa_ov"]("i,a"), u1_aa_vo("a,i"));
        scalar13 = dot(F_blks_["bb_ov"]("i,a"), u1_bb_vo("a,i"));
        scalar14 = dot(V_blks_["aaaa_oovv"]("j,i,a,b"), u2_aaaa_vvoo("a,b,j,i"));
        scalar15 = dot(V_blks_["abab_oovv"]("j,i,a,b"), u2_abab_vvoo("a,b,j,i"));
        scalar16 = dot(V_blks_["bbbb_oovv"]("j,i,a,b"), u2_bbbb_vvoo("a,b,j,i"));
        scalar17 = dot(dp_aa_oo("o0,o1"), Id_blks_["aa_oo"]("o0,o1"));
        scalar18 = dot(dp_bb_oo("o0,o1"), Id_blks_["bb_oo"]("o0,o1"));


        /// ****** sigmar1_a(e) ****** ///

        {

            // sigmar1_xa_Lv += 1.000000 f_aa(i,i) r1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xa_Lv("I,e") += scalar5 * r1_xa_Lv("I,e");

            // sigmar1_xa_Lv += 1.000000 f_bb(i,i) r1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xa_Lv("I,e") += scalar6 * r1_xa_Lv("I,e");

            // sigmar1_xa_Lv += 1.000000 f_aa(e,a) r1_a(a) 
            // flops: o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmar1_xa_Lv("I,e") += F_blks_["aa_vv"]("e,a") * r1_xa_Lv("I,a");

            // sigmar1_xa_Lv += -1.000000 f_aa(i,a) r2_aaa(a,e,i) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmar1_xa_Lv("I,e") -= F_blks_["aa_ov"]("i,a") * r2_xaaa_Lvvo("I,a,e,i");
        }

        if (include_u1_) {

            // sigmar1_xa_Lv += -1.000000 d-_aa(i,i) s1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xa_Lv("I,e") -= scalar7 * s1_xa_Lv("I,e");

            // sigmar1_xa_Lv += -1.000000 d-_bb(i,i) s1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xa_Lv("I,e") -= scalar8 * s1_xa_Lv("I,e");

            // sigmar1_xa_Lv += -1.000000 d-_aa(e,a) s1_a(a) 
            // flops: o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmar1_xa_Lv("I,e") -= dp_aa_vv("e,a") * s1_xa_Lv("I,a");
        }

        if (include_u2_) {

            // sigmar1_xa_Lv += 1.000000 d-_aa(i,a) s2_aaa(a,e,i) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmar1_xa_Lv("I,e") += dp_aa_ov("i,a") * s2_xaaa_Lvvo("I,a,e,i");
        }

        if (include_u0_) {

            // sigmar1_xa_Lv += -1.000000 d-_aa(i,i) r1_a(e) u0 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmar1_xa_Lv("I,e") -= scalar7 * u0 * r1_xa_Lv("I,e");

            // sigmar1_xa_Lv += -1.000000 d-_bb(i,i) r1_a(e) u0 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmar1_xa_Lv("I,e") -= scalar8 * u0 * r1_xa_Lv("I,e");

            // sigmar1_xa_Lv += -1.000000 d-_aa(e,a) r1_a(a) u0 
            // flops: o0v2L1: 1, o0v1L1: 2 | mem: o0v1L1: 2, o0v1: 1, 
            sigmar1_xa_Lv("I,e") -= dp_aa_vv("e,a") * r1_xa_Lv("I,a") * u0;

            // sigmar1_xa_Lv += 1.000000 d-_aa(i,a) r2_aaa(a,e,i) u0 
            // flops: o1v2L1: 1, o0v1L1: 2 | mem: o0v1L1: 2, o0v1: 1, 
            sigmar1_xa_Lv("I,e") += dp_aa_ov("i,a") * r2_xaaa_Lvvo("I,a,e,i") * u0;
        }

        if (include_u1_) {

            // sigmar1_xa_Lv += -1.000000 d-_aa(i,a) r1_a(e) u1_aa(a,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xa_Lv("I,e") -= scalar0 * r1_xa_Lv("I,e");

            // sigmar1_xa_Lv += -1.000000 d-_bb(i,a) r1_a(e) u1_bb(a,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xa_Lv("I,e") -= scalar1 * r1_xa_Lv("I,e");

            // sigmar1_xa_Lv += 1.000000 d-_aa(i,a) r1_a(a) u1_aa(e,i) 
            // flops: o1v1L1: 2, o0v1L1: 1 | mem: o0v1L1: 2, o1v0L1: 1, 
            sigmar1_xa_Lv("I,e") += dp_aa_ov("i,a") * r1_xa_Lv("I,a") * u1_aa_vo("e,i");
        }

        {

            // sigmar1_xa_Lv += -0.500000 <j,i||j,i>_aaaa r1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xa_Lv("I,e") -= 0.500000 * scalar9 * r1_xa_Lv("I,e");

            // sigmar1_xa_Lv += -0.500000 <j,i||j,i>_abab r1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xa_Lv("I,e") -= 0.500000 * scalar10 * r1_xa_Lv("I,e");

            // sigmar1_xa_Lv += -0.500000 <i,j||i,j>_abab r1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xa_Lv("I,e") -= 0.500000 * scalar10 * r1_xa_Lv("I,e");

            // sigmar1_xa_Lv += -0.500000 <j,i||j,i>_bbbb r1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xa_Lv("I,e") -= 0.500000 * scalar11 * r1_xa_Lv("I,e");

            // sigmar1_xa_Lv += -0.500000 <i,e||a,b>_aaaa r2_aaa(a,b,i) 
            // flops: o1v3L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmar1_xa_Lv("I,e") += 0.500000 * V_blks_["aaaa_vovv"]("e,i,a,b") * r2_xaaa_Lvvo("I,a,b,i");

            // sigmar1_xa_Lv += 0.250000 <j,i||a,b>_aaaa r1_a(e) t2_aaaa(a,b,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xa_Lv("I,e") += 0.250000 * scalar2 * r1_xa_Lv("I,e");

            // sigmar1_xa_Lv += 0.250000 <j,i||a,b>_abab r1_a(e) t2_abab(a,b,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xa_Lv("I,e") += 0.250000 * scalar3 * r1_xa_Lv("I,e");

            // sigmar1_xa_Lv += 0.250000 <i,j||a,b>_abab r1_a(e) t2_abab(a,b,i,j) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xa_Lv("I,e") += 0.250000 * scalar3 * r1_xa_Lv("I,e");

            // sigmar1_xa_Lv += 0.250000 <j,i||b,a>_abab r1_a(e) t2_abab(b,a,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xa_Lv("I,e") += 0.250000 * scalar3 * r1_xa_Lv("I,e");

            // sigmar1_xa_Lv += 0.250000 <i,j||b,a>_abab r1_a(e) t2_abab(b,a,i,j) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xa_Lv("I,e") += 0.250000 * scalar3 * r1_xa_Lv("I,e");

            // sigmar1_xa_Lv += 0.250000 <j,i||a,b>_bbbb r1_a(e) t2_bbbb(a,b,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xa_Lv("I,e") += 0.250000 * scalar4 * r1_xa_Lv("I,e");

            // sigmar1_xa_Lv += -0.500000 <j,i||a,b>_aaaa r1_a(b) t2_aaaa(a,e,j,i) 
            // flops: o2v2L1: 2, o0v1L1: 1 | mem: o2v1L1: 1, o0v1L1: 2, 
            sigmar1_xa_Lv("I,e") -= 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r1_xa_Lv("I,b") * t2_aaaa_vvoo("a,e,j,i");

            // sigmar1_xa_Lv += -0.500000 <j,i||b,a>_abab r1_a(b) t2_abab(e,a,j,i) 
            // flops: o2v2L1: 2, o0v1L1: 1 | mem: o2v1L1: 1, o0v1L1: 2, 
            sigmar1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("j,i,b,a") * r1_xa_Lv("I,b") * t2_abab_vvoo("e,a,j,i");

            // sigmar1_xa_Lv += -0.500000 <i,j||b,a>_abab r1_a(b) t2_abab(e,a,i,j) 
            // flops: o2v2L1: 2, o0v1L1: 1 | mem: o2v1L1: 1, o0v1L1: 2, 
            sigmar1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("i,j,b,a") * r1_xa_Lv("I,b") * t2_abab_vvoo("e,a,i,j");
        }




/// ****** sigmar1_b(e) ****** ///



        {

            // sigmar1_xb_Lv += 1.000000 f_aa(i,i) r1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xb_Lv("I,e") += scalar5 * r1_xb_Lv("I,e");

            // sigmar1_xb_Lv += 1.000000 f_bb(i,i) r1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xb_Lv("I,e") += scalar6 * r1_xb_Lv("I,e");

            // sigmar1_xb_Lv += 1.000000 f_bb(e,a) r1_b(a) 
            // flops: o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmar1_xb_Lv("I,e") += F_blks_["bb_vv"]("e,a") * r1_xb_Lv("I,a");

            // sigmar1_xb_Lv += -1.000000 f_bb(i,a) r2_bbb(a,e,i) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmar1_xb_Lv("I,e") -= F_blks_["bb_ov"]("i,a") * r2_xbbb_Lvvo("I,a,e,i");
        }

        if (include_u1_) {

            // sigmar1_xb_Lv += -1.000000 d-_aa(i,i) s1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xb_Lv("I,e") -= scalar7 * s1_xb_Lv("I,e");

            // sigmar1_xb_Lv += -1.000000 d-_bb(i,i) s1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xb_Lv("I,e") -= scalar8 * s1_xb_Lv("I,e");

            // sigmar1_xb_Lv += -1.000000 d-_bb(e,a) s1_b(a) 
            // flops: o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmar1_xb_Lv("I,e") -= dp_bb_vv("e,a") * s1_xb_Lv("I,a");
        }

        if (include_u2_) {

            // sigmar1_xb_Lv += 1.000000 d-_bb(i,a) s2_bbb(a,e,i) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmar1_xb_Lv("I,e") += dp_bb_ov("i,a") * s2_xbbb_Lvvo("I,a,e,i");
        }

        if (include_u0_) {

            // sigmar1_xb_Lv += -1.000000 d-_aa(i,i) r1_b(e) u0 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmar1_xb_Lv("I,e") -= scalar7 * u0 * r1_xb_Lv("I,e");

            // sigmar1_xb_Lv += -1.000000 d-_bb(i,i) r1_b(e) u0 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmar1_xb_Lv("I,e") -= scalar8 * u0 * r1_xb_Lv("I,e");

            // sigmar1_xb_Lv += -1.000000 d-_bb(e,a) r1_b(a) u0 
            // flops: o0v2L1: 1, o0v1L1: 2 | mem: o0v1L1: 2, o0v1: 1, 
            sigmar1_xb_Lv("I,e") -= dp_bb_vv("e,a") * r1_xb_Lv("I,a") * u0;

            // sigmar1_xb_Lv += 1.000000 d-_bb(i,a) r2_bbb(a,e,i) u0 
            // flops: o1v2L1: 1, o0v1L1: 2 | mem: o0v1L1: 2, o0v1: 1, 
            sigmar1_xb_Lv("I,e") += dp_bb_ov("i,a") * r2_xbbb_Lvvo("I,a,e,i") * u0;
        }

        if (include_u1_) {

            // sigmar1_xb_Lv += -1.000000 d-_aa(i,a) r1_b(e) u1_aa(a,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xb_Lv("I,e") -= scalar0 * r1_xb_Lv("I,e");

            // sigmar1_xb_Lv += -1.000000 d-_bb(i,a) r1_b(e) u1_bb(a,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xb_Lv("I,e") -= scalar1 * r1_xb_Lv("I,e");

            // sigmar1_xb_Lv += 1.000000 d-_bb(i,a) r1_b(a) u1_bb(e,i) 
            // flops: o1v1L1: 2, o0v1L1: 1 | mem: o0v1L1: 2, o1v0L1: 1, 
            sigmar1_xb_Lv("I,e") += dp_bb_ov("i,a") * r1_xb_Lv("I,a") * u1_bb_vo("e,i");
        }

        {

            // sigmar1_xb_Lv += -0.500000 <j,i||j,i>_aaaa r1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xb_Lv("I,e") -= 0.500000 * scalar9 * r1_xb_Lv("I,e");

            // sigmar1_xb_Lv += -0.500000 <j,i||j,i>_abab r1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xb_Lv("I,e") -= 0.500000 * scalar10 * r1_xb_Lv("I,e");

            // sigmar1_xb_Lv += -0.500000 <i,j||i,j>_abab r1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xb_Lv("I,e") -= 0.500000 * scalar10 * r1_xb_Lv("I,e");

            // sigmar1_xb_Lv += -0.500000 <j,i||j,i>_bbbb r1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xb_Lv("I,e") -= 0.500000 * scalar11 * r1_xb_Lv("I,e");

            // sigmar1_xb_Lv += -0.500000 <i,e||a,b>_bbbb r2_bbb(a,b,i) 
            // flops: o1v3L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmar1_xb_Lv("I,e") += 0.500000 * V_blks_["bbbb_vovv"]("e,i,a,b") * r2_xbbb_Lvvo("I,a,b,i");

            // sigmar1_xb_Lv += 0.250000 <j,i||a,b>_aaaa r1_b(e) t2_aaaa(a,b,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xb_Lv("I,e") += 0.250000 * scalar2 * r1_xb_Lv("I,e");

            // sigmar1_xb_Lv += 0.250000 <j,i||a,b>_abab r1_b(e) t2_abab(a,b,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xb_Lv("I,e") += 0.250000 * scalar3 * r1_xb_Lv("I,e");

            // sigmar1_xb_Lv += 0.250000 <i,j||a,b>_abab r1_b(e) t2_abab(a,b,i,j) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xb_Lv("I,e") += 0.250000 * scalar3 * r1_xb_Lv("I,e");

            // sigmar1_xb_Lv += 0.250000 <j,i||b,a>_abab r1_b(e) t2_abab(b,a,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xb_Lv("I,e") += 0.250000 * scalar3 * r1_xb_Lv("I,e");

            // sigmar1_xb_Lv += 0.250000 <i,j||b,a>_abab r1_b(e) t2_abab(b,a,i,j) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xb_Lv("I,e") += 0.250000 * scalar3 * r1_xb_Lv("I,e");

            // sigmar1_xb_Lv += 0.250000 <j,i||a,b>_bbbb r1_b(e) t2_bbbb(a,b,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmar1_xb_Lv("I,e") += 0.250000 * scalar4 * r1_xb_Lv("I,e");

            // sigmar1_xb_Lv += -0.500000 <j,i||a,b>_abab r1_b(b) t2_abab(a,e,j,i) 
            // flops: o2v2L1: 2, o0v1L1: 1 | mem: o2v1L1: 1, o0v1L1: 2, 
            sigmar1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("j,i,a,b") * r1_xb_Lv("I,b") * t2_abab_vvoo("a,e,j,i");

            // sigmar1_xb_Lv += -0.500000 <i,j||a,b>_abab r1_b(b) t2_abab(a,e,i,j) 
            // flops: o2v2L1: 2, o0v1L1: 1 | mem: o2v1L1: 1, o0v1L1: 2, 
            sigmar1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("i,j,a,b") * r1_xb_Lv("I,b") * t2_abab_vvoo("a,e,i,j");

            // sigmar1_xb_Lv += -0.500000 <j,i||a,b>_bbbb r1_b(b) t2_bbbb(a,e,j,i) 
            // flops: o2v2L1: 2, o0v1L1: 1 | mem: o2v1L1: 1, o0v1L1: 2, 
            sigmar1_xb_Lv("I,e") -= 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r1_xb_Lv("I,b") * t2_bbbb_vvoo("a,e,j,i");
        }




/// ****** sigmar2_aaa(e,f,n) ****** ///



        {

            // sigmar2_xaaa_Lvvo += -1.000000 P(e,f) f_aa(e,n) r1_a(f) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = F_blks_["aa_vo"]("e,n") * r1_xa_Lv("I,f");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += 1.000000 f_aa(i,i) r2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") += scalar5 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmar2_xaaa_Lvvo += 1.000000 f_bb(i,i) r2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") += scalar6 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmar2_xaaa_Lvvo += -1.000000 f_aa(i,n) r2_aaa(e,f,i) 
            // flops: o2v2L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmar2_xaaa_Lvvo("I,e,f,n") -= F_blks_["aa_oo"]("i,n") * r2_xaaa_Lvvo("I,e,f,i");

            // sigmar2_xaaa_Lvvo += 1.000000 P(e,f) f_aa(e,a) r2_aaa(a,f,n) 
            // flops: o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = F_blks_["aa_vv"]("e,a") * r2_xaaa_Lvvo("I,a,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += 1.000000 P(e,f) f_aa(i,a) r1_a(f) t2_aaaa(a,e,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = F_blks_["aa_ov"]("i,a") * r1_xa_Lv("I,f") * t2_aaaa_vvoo("a,e,n,i");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += -1.000000 P(e,f) f_bb(i,a) r1_a(f) t2_abab(e,a,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = F_blks_["bb_ov"]("i,a") * r1_xa_Lv("I,f") * t2_abab_vvoo("e,a,n,i");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += 1.000000 f_aa(i,a) r1_a(a) t2_aaaa(e,f,n,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o1v1L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") += F_blks_["aa_ov"]("i,a") * r1_xa_Lv("I,a") * t2_aaaa_vvoo("e,f,n,i");
        }

        if (include_u1_) {

            // sigmar2_xaaa_Lvvo += 1.000000 P(e,f) d-_aa(e,n) s1_a(f) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_vo("e,n") * s1_xa_Lv("I,f");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmar2_xaaa_Lvvo += -1.000000 d-_aa(i,i) s2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") -= scalar7 * s2_xaaa_Lvvo("I,e,f,n");

            // sigmar2_xaaa_Lvvo += -1.000000 d-_bb(i,i) s2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") -= scalar8 * s2_xaaa_Lvvo("I,e,f,n");

            // sigmar2_xaaa_Lvvo += 1.000000 d-_aa(i,n) s2_aaa(e,f,i) 
            // flops: o2v2L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmar2_xaaa_Lvvo("I,e,f,n") += dp_aa_oo("i,n") * s2_xaaa_Lvvo("I,e,f,i");

            // sigmar2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(e,a) s2_aaa(a,f,n) 
            // flops: o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_vv("e,a") * s2_xaaa_Lvvo("I,a,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u1_) {

            // sigmar2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(i,a) t2_aaaa(a,e,n,i) s1_a(f) 
            // flops: o1v2L1: 3, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("i,a") * t2_aaaa_vvoo("a,e,n,i") * s1_xa_Lv("I,f");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += 1.000000 P(e,f) d-_bb(i,a) t2_abab(e,a,n,i) s1_a(f) 
            // flops: o1v2L1: 3, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_bb_ov("i,a") * t2_abab_vvoo("e,a,n,i") * s1_xa_Lv("I,f");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += -1.000000 d-_aa(i,a) t2_aaaa(e,f,n,i) s1_a(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") -= dp_aa_ov("i,a") * t2_aaaa_vvoo("e,f,n,i") * s1_xa_Lv("I,a");
        }

        if (include_u0_) {

            // sigmar2_xaaa_Lvvo += 1.000000 P(e,f) d-_aa(e,n) r1_a(f) u0 
            // flops: o1v2L1: 4 | mem: o1v2L1: 2, o1v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_vo("e,n") * r1_xa_Lv("I,f") * u0;
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += -1.000000 d-_aa(i,i) r2_aaa(e,f,n) u0 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") -= scalar7 * u0 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmar2_xaaa_Lvvo += -1.000000 d-_bb(i,i) r2_aaa(e,f,n) u0 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") -= scalar8 * u0 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmar2_xaaa_Lvvo += 1.000000 d-_aa(i,n) r2_aaa(e,f,i) u0 
            // flops: o2v2L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v2: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") += dp_aa_oo("i,n") * r2_xaaa_Lvvo("I,e,f,i") * u0;

            // sigmar2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(e,a) r2_aaa(a,f,n) u0 
            // flops: o1v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_vv("e,a") * r2_xaaa_Lvvo("I,a,f,n") * u0;
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u1_) {

            // sigmar2_xaaa_Lvvo += 1.000000 P(e,f) d-_aa(i,i) r1_a(f) u1_aa(e,n) 
            // flops: o1v2L1: 2, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = r1_xa_Lv("I,f") * scalar7 * u1_aa_vo("e,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += 1.000000 P(e,f) d-_bb(i,i) r1_a(f) u1_aa(e,n) 
            // flops: o1v2L1: 2, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = r1_xa_Lv("I,f") * scalar8 * u1_aa_vo("e,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(i,n) r1_a(f) u1_aa(e,i) 
            // flops: o2v2L1: 1, o1v2L1: 2, o2v1L1: 1 | mem: o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_oo("i,n") * r1_xa_Lv("I,f") * u1_aa_vo("e,i");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += 1.000000 P(e,f) d-_aa(e,a) r1_a(f) u1_aa(a,n) 
            // flops: o1v3L1: 1, o0v3L1: 1, o1v2L1: 2 | mem: o0v3L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_vv("e,a") * r1_xa_Lv("I,f") * u1_aa_vo("a,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(e,a) r1_a(a) u1_aa(f,n) 
            // flops: o1v2L1: 3, o0v2L1: 1 | mem: o1v2L1: 2, o0v1L1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_vv("e,a") * r1_xa_Lv("I,a") * u1_aa_vo("f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += -1.000000 d-_aa(i,a) r2_aaa(e,f,n) u1_aa(a,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") -= scalar0 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmar2_xaaa_Lvvo += -1.000000 d-_bb(i,a) r2_aaa(e,f,n) u1_bb(a,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") -= scalar1 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmar2_xaaa_Lvvo += 1.000000 d-_aa(i,a) r2_aaa(e,f,i) u1_aa(a,n) 
            // flops: o1v3L1: 2, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, 
            sigmar2_xaaa_Lvvo("I,e,f,n") += dp_aa_ov("i,a") * r2_xaaa_Lvvo("I,e,f,i") * u1_aa_vo("a,n");

            // sigmar2_xaaa_Lvvo += 1.000000 P(e,f) d-_aa(i,a) r2_aaa(a,f,n) u1_aa(e,i) 
            // flops: o2v2L1: 2, o1v2L1: 2 | mem: o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("i,a") * r2_xaaa_Lvvo("I,a,f,n") * u1_aa_vo("e,i");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(i,a) r2_aaa(a,f,i) u1_aa(e,n) 
            // flops: o1v2L1: 4 | mem: o1v2L1: 2, o0v1L1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("i,a") * r2_xaaa_Lvvo("I,a,f,i") * u1_aa_vo("e,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmar2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(i,a) r1_a(f) u2_aaaa(a,e,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("i,a") * r1_xa_Lv("I,f") * u2_aaaa_vvoo("a,e,n,i");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += 1.000000 P(e,f) d-_bb(i,a) r1_a(f) u2_abab(e,a,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_bb_ov("i,a") * r1_xa_Lv("I,f") * u2_abab_vvoo("e,a,n,i");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += -1.000000 d-_aa(i,a) r1_a(a) u2_aaaa(e,f,n,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o1v1L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") -= dp_aa_ov("i,a") * r1_xa_Lv("I,a") * u2_aaaa_vvoo("e,f,n,i");
        }

        if (include_u0_) {

            // sigmar2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(i,a) r1_a(f) t2_aaaa(a,e,n,i) u0 
            // flops: o2v3L1: 1, o1v2L1: 4 | mem: o1v2L1: 3, o1v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("i,a") * r1_xa_Lv("I,f") * t2_aaaa_vvoo("a,e,n,i") * u0;
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += 1.000000 P(e,f) d-_bb(i,a) r1_a(f) t2_abab(e,a,n,i) u0 
            // flops: o2v3L1: 1, o1v2L1: 4 | mem: o1v2L1: 3, o1v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_bb_ov("i,a") * r1_xa_Lv("I,f") * t2_abab_vvoo("e,a,n,i") * u0;
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += -1.000000 d-_aa(i,a) r1_a(a) t2_aaaa(e,f,n,i) u0 
            // flops: o2v2L1: 1, o1v2L1: 2, o1v1L1: 1 | mem: o1v2L1: 2, o1v2: 1, o1v0L1: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") -= dp_aa_ov("i,a") * r1_xa_Lv("I,a") * t2_aaaa_vvoo("e,f,n,i") * u0;
        }

        {

            // sigmar2_xaaa_Lvvo += -0.500000 <j,i||j,i>_aaaa r2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * scalar9 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmar2_xaaa_Lvvo += -0.500000 <j,i||j,i>_abab r2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * scalar10 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmar2_xaaa_Lvvo += -0.500000 <i,j||i,j>_abab r2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * scalar10 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmar2_xaaa_Lvvo += -0.500000 <j,i||j,i>_bbbb r2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * scalar11 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmar2_xaaa_Lvvo += 1.000000 <e,f||a,n>_aaaa r1_a(a) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmar2_xaaa_Lvvo("I,e,f,n") += V_blks_["aaaa_vvvo"]("e,f,a,n") * r1_xa_Lv("I,a");

            // sigmar2_xaaa_Lvvo += 1.000000 P(e,f) <i,e||a,n>_aaaa r2_aaa(a,f,i) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_vovo"]("e,i,a,n") * r2_xaaa_Lvvo("I,a,f,i");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += 0.500000 <e,f||a,b>_aaaa r2_aaa(a,b,n) 
            // flops: o1v4L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmar2_xaaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["aaaa_vvvv"]("e,f,a,b") * r2_xaaa_Lvvo("I,a,b,n");

            // sigmar2_xaaa_Lvvo += 0.500000 P(e,f) <j,i||a,n>_aaaa r1_a(f) t2_aaaa(a,e,j,i) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["aaaa_oovo"]("j,i,a,n") * r1_xa_Lv("I,f") * t2_aaaa_vvoo("a,e,j,i");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += 0.500000 P(e,f) <j,i||n,a>_abab r1_a(f) t2_abab(e,a,j,i) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abba_oovo"]("j,i,a,n") * r1_xa_Lv("I,f") * t2_abab_vvoo("e,a,j,i");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += 0.500000 P(e,f) <i,j||n,a>_abab r1_a(f) t2_abab(e,a,i,j) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abba_oovo"]("i,j,a,n") * r1_xa_Lv("I,f") * t2_abab_vvoo("e,a,i,j");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += 0.500000 <j,i||a,n>_aaaa r1_a(a) t2_aaaa(e,f,j,i) 
            // flops: o3v2L1: 1, o3v1L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o3v0L1: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["aaaa_oovo"]("j,i,a,n") * r1_xa_Lv("I,a") * t2_aaaa_vvoo("e,f,j,i");

            // sigmar2_xaaa_Lvvo += 0.500000 P(e,f) <i,e||a,b>_aaaa r1_a(f) t2_aaaa(a,b,n,i) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 2 | mem: o1v4L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["aaaa_vovv"]("e,i,a,b") * r1_xa_Lv("I,f") * t2_aaaa_vvoo("a,b,n,i");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += -0.500000 P(e,f) <e,i||a,b>_abab r1_a(f) t2_abab(a,b,n,i) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 2 | mem: o1v4L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_vovv"]("e,i,a,b") * r1_xa_Lv("I,f") * t2_abab_vvoo("a,b,n,i");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += -0.500000 P(e,f) <e,i||b,a>_abab r1_a(f) t2_abab(b,a,n,i) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 2 | mem: o1v4L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_vovv"]("e,i,b,a") * r1_xa_Lv("I,f") * t2_abab_vvoo("b,a,n,i");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += -1.000000 P(e,f) <i,e||a,b>_aaaa r1_a(b) t2_aaaa(a,f,n,i) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_vovv"]("e,i,a,b") * r1_xa_Lv("I,b") * t2_aaaa_vvoo("a,f,n,i");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += 1.000000 P(e,f) <e,i||b,a>_abab r1_a(b) t2_abab(f,a,n,i) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_vovv"]("e,i,b,a") * r1_xa_Lv("I,b") * t2_abab_vvoo("f,a,n,i");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += 0.250000 <j,i||a,b>_aaaa r2_aaa(e,f,n) t2_aaaa(a,b,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar2 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmar2_xaaa_Lvvo += 0.250000 <j,i||a,b>_abab r2_aaa(e,f,n) t2_abab(a,b,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar3 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmar2_xaaa_Lvvo += 0.250000 <i,j||a,b>_abab r2_aaa(e,f,n) t2_abab(a,b,i,j) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar3 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmar2_xaaa_Lvvo += 0.250000 <j,i||b,a>_abab r2_aaa(e,f,n) t2_abab(b,a,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar3 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmar2_xaaa_Lvvo += 0.250000 <i,j||b,a>_abab r2_aaa(e,f,n) t2_abab(b,a,i,j) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar3 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmar2_xaaa_Lvvo += 0.250000 <j,i||a,b>_bbbb r2_aaa(e,f,n) t2_bbbb(a,b,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar4 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmar2_xaaa_Lvvo += -0.500000 <j,i||a,b>_aaaa r2_aaa(e,f,j) t2_aaaa(a,b,n,i) 
            // flops: o2v4L1: 2, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmar2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaa_Lvvo("I,e,f,j") * t2_aaaa_vvoo("a,b,n,i");

            // sigmar2_xaaa_Lvvo += -0.500000 <j,i||a,b>_abab r2_aaa(e,f,j) t2_abab(a,b,n,i) 
            // flops: o2v4L1: 2, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmar2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("j,i,a,b") * r2_xaaa_Lvvo("I,e,f,j") * t2_abab_vvoo("a,b,n,i");

            // sigmar2_xaaa_Lvvo += -0.500000 <j,i||b,a>_abab r2_aaa(e,f,j) t2_abab(b,a,n,i) 
            // flops: o2v4L1: 2, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmar2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("j,i,b,a") * r2_xaaa_Lvvo("I,e,f,j") * t2_abab_vvoo("b,a,n,i");

            // sigmar2_xaaa_Lvvo += -0.500000 P(e,f) <j,i||a,b>_aaaa r2_aaa(b,f,n) t2_aaaa(a,e,j,i) 
            // flops: o3v3L1: 2, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaa_Lvvo("I,b,f,n") * t2_aaaa_vvoo("a,e,j,i");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += -0.500000 P(e,f) <j,i||b,a>_abab r2_aaa(b,f,n) t2_abab(e,a,j,i) 
            // flops: o3v3L1: 2, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("j,i,b,a") * r2_xaaa_Lvvo("I,b,f,n") * t2_abab_vvoo("e,a,j,i");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += -0.500000 P(e,f) <i,j||b,a>_abab r2_aaa(b,f,n) t2_abab(e,a,i,j) 
            // flops: o3v3L1: 2, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("i,j,b,a") * r2_xaaa_Lvvo("I,b,f,n") * t2_abab_vvoo("e,a,i,j");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += 1.000000 P(e,f) <j,i||a,b>_aaaa r2_aaa(b,f,j) t2_aaaa(a,e,n,i) 
            // flops: o2v3L1: 2, o1v2L1: 2 | mem: o1v2L1: 3, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaa_Lvvo("I,b,f,j") * t2_aaaa_vvoo("a,e,n,i");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += 1.000000 P(e,f) <j,i||b,a>_abab r2_aaa(b,f,j) t2_abab(e,a,n,i) 
            // flops: o2v3L1: 2, o1v2L1: 2 | mem: o1v2L1: 3, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("j,i,b,a") * r2_xaaa_Lvvo("I,b,f,j") * t2_abab_vvoo("e,a,n,i");
            sigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmar2_xaaa_Lvvo += 0.250000 <j,i||a,b>_aaaa r2_aaa(a,b,n) t2_aaaa(e,f,j,i) 
            // flops: o3v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o3v0L1: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") += 0.250000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaa_Lvvo("I,a,b,n") * t2_aaaa_vvoo("e,f,j,i");

            // sigmar2_xaaa_Lvvo += -0.500000 <j,i||a,b>_aaaa r2_aaa(a,b,j) t2_aaaa(e,f,n,i) 
            // flops: o2v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmar2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaa_Lvvo("I,a,b,j") * t2_aaaa_vvoo("e,f,n,i");
        }




/// ****** sigmar2_abb(e,f,n) ****** ///



        {

            // sigmar2_xabb_Lvvo += 1.000000 f_bb(f,n) r1_a(e) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            sigmar2_xabb_Lvvo("I,e,f,n") += F_blks_["bb_vo"]("f,n") * r1_xa_Lv("I,e");

            // sigmar2_xabb_Lvvo += 1.000000 f_aa(i,a) r1_a(e) t2_abab(a,f,i,n) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            sigmar2_xabb_Lvvo("I,e,f,n") += F_blks_["aa_ov"]("i,a") * r1_xa_Lv("I,e") * t2_abab_vvoo("a,f,i,n");

            // sigmar2_xabb_Lvvo += -1.000000 f_bb(i,a) r1_a(e) t2_bbbb(a,f,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            sigmar2_xabb_Lvvo("I,e,f,n") -= F_blks_["bb_ov"]("i,a") * r1_xa_Lv("I,e") * t2_bbbb_vvoo("a,f,n,i");

            // sigmar2_xabb_Lvvo += -1.000000 f_aa(i,a) r1_a(a) t2_abab(e,f,i,n) 
            // flops: o2v2L1: 1, o1v2L1: 1, o1v1L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmar2_xabb_Lvvo("I,e,f,n") -= F_blks_["aa_ov"]("i,a") * r1_xa_Lv("I,a") * t2_abab_vvoo("e,f,i,n");
        }

        if (include_u1_) {

            // sigmar2_xabb_Lvvo += -1.000000 d-_bb(f,n) s1_a(e) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            sigmar2_xabb_Lvvo("I,e,f,n") -= dp_bb_vo("f,n") * s1_xa_Lv("I,e");

            // sigmar2_xabb_Lvvo += -1.000000 d-_aa(i,a) t2_abab(a,f,i,n) s1_a(e) 
            // flops: o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmar2_xabb_Lvvo("I,e,f,n") -= dp_aa_ov("i,a") * t2_abab_vvoo("a,f,i,n") * s1_xa_Lv("I,e");

            // sigmar2_xabb_Lvvo += 1.000000 d-_bb(i,a) t2_bbbb(a,f,n,i) s1_a(e) 
            // flops: o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmar2_xabb_Lvvo("I,e,f,n") += dp_bb_ov("i,a") * t2_bbbb_vvoo("a,f,n,i") * s1_xa_Lv("I,e");

            // sigmar2_xabb_Lvvo += 1.000000 d-_aa(i,a) t2_abab(e,f,i,n) s1_a(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmar2_xabb_Lvvo("I,e,f,n") += dp_aa_ov("i,a") * t2_abab_vvoo("e,f,i,n") * s1_xa_Lv("I,a");
        }

        if (include_u0_) {

            // sigmar2_xabb_Lvvo += -1.000000 d-_bb(f,n) r1_a(e) u0 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, o1v2: 1, 
            sigmar2_xabb_Lvvo("I,e,f,n") -= dp_bb_vo("f,n") * r1_xa_Lv("I,e") * u0;
        }

        if (include_u1_) {

            // sigmar2_xabb_Lvvo += -1.000000 d-_aa(i,i) r1_a(e) u1_bb(f,n) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            sigmar2_xabb_Lvvo("I,e,f,n") -= r1_xa_Lv("I,e") * scalar7 * u1_bb_vo("f,n");

            // sigmar2_xabb_Lvvo += -1.000000 d-_bb(i,i) r1_a(e) u1_bb(f,n) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            sigmar2_xabb_Lvvo("I,e,f,n") -= r1_xa_Lv("I,e") * scalar8 * u1_bb_vo("f,n");

            // sigmar2_xabb_Lvvo += 1.000000 d-_bb(i,n) r1_a(e) u1_bb(f,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o2v1L1: 1 | mem: o1v2L1: 2, o2v1L1: 1, 
            sigmar2_xabb_Lvvo("I,e,f,n") += dp_bb_oo("i,n") * r1_xa_Lv("I,e") * u1_bb_vo("f,i");

            // sigmar2_xabb_Lvvo += -1.000000 d-_bb(f,a) r1_a(e) u1_bb(a,n) 
            // flops: o1v3L1: 1, o0v3L1: 1, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, 
            sigmar2_xabb_Lvvo("I,e,f,n") -= dp_bb_vv("f,a") * r1_xa_Lv("I,e") * u1_bb_vo("a,n");

            // sigmar2_xabb_Lvvo += -1.000000 d-_aa(e,a) r1_a(a) u1_bb(f,n) 
            // flops: o1v2L1: 2, o0v2L1: 1 | mem: o1v2L1: 2, o0v1L1: 1, 
            sigmar2_xabb_Lvvo("I,e,f,n") -= dp_aa_vv("e,a") * r1_xa_Lv("I,a") * u1_bb_vo("f,n");

            // sigmar2_xabb_Lvvo += 1.000000 d-_aa(i,a) r2_aaa(a,e,i) u1_bb(f,n) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, o0v1L1: 1, 
            sigmar2_xabb_Lvvo("I,e,f,n") += dp_aa_ov("i,a") * r2_xaaa_Lvvo("I,a,e,i") * u1_bb_vo("f,n");
        }

        if (include_u2_) {

            // sigmar2_xabb_Lvvo += -1.000000 d-_aa(i,a) r1_a(e) u2_abab(a,f,i,n) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            sigmar2_xabb_Lvvo("I,e,f,n") -= dp_aa_ov("i,a") * r1_xa_Lv("I,e") * u2_abab_vvoo("a,f,i,n");

            // sigmar2_xabb_Lvvo += 1.000000 d-_bb(i,a) r1_a(e) u2_bbbb(a,f,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            sigmar2_xabb_Lvvo("I,e,f,n") += dp_bb_ov("i,a") * r1_xa_Lv("I,e") * u2_bbbb_vvoo("a,f,n,i");

            // sigmar2_xabb_Lvvo += 1.000000 d-_aa(i,a) r1_a(a) u2_abab(e,f,i,n) 
            // flops: o2v2L1: 1, o1v2L1: 1, o1v1L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmar2_xabb_Lvvo("I,e,f,n") += dp_aa_ov("i,a") * r1_xa_Lv("I,a") * u2_abab_vvoo("e,f,i,n");
        }

        if (include_u0_) {

            // sigmar2_xabb_Lvvo += -1.000000 d-_aa(i,a) r1_a(e) t2_abab(a,f,i,n) u0 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, o1v2: 1, 
            sigmar2_xabb_Lvvo("I,e,f,n") -= dp_aa_ov("i,a") * r1_xa_Lv("I,e") * t2_abab_vvoo("a,f,i,n") * u0;

            // sigmar2_xabb_Lvvo += 1.000000 d-_bb(i,a) r1_a(e) t2_bbbb(a,f,n,i) u0 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, o1v2: 1, 
            sigmar2_xabb_Lvvo("I,e,f,n") += dp_bb_ov("i,a") * r1_xa_Lv("I,e") * t2_bbbb_vvoo("a,f,n,i") * u0;

            // sigmar2_xabb_Lvvo += 1.000000 d-_aa(i,a) r1_a(a) t2_abab(e,f,i,n) u0 
            // flops: o2v2L1: 1, o1v2L1: 2, o1v1L1: 1 | mem: o1v2L1: 2, o1v2: 1, o1v0L1: 1, 
            sigmar2_xabb_Lvvo("I,e,f,n") += dp_aa_ov("i,a") * r1_xa_Lv("I,a") * t2_abab_vvoo("e,f,i,n") * u0;
        }

        {

            // sigmar2_xabb_Lvvo += 1.000000 <e,f||a,n>_abab r1_a(a) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmar2_xabb_Lvvo("I,e,f,n") += V_blks_["abab_vvvo"]("e,f,a,n") * r1_xa_Lv("I,a");

            // sigmar2_xabb_Lvvo += -1.000000 <i,f||a,n>_abab r2_aaa(a,e,i) 
            // flops: o2v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmar2_xabb_Lvvo("I,e,f,n") += V_blks_["baab_vovo"]("f,i,a,n") * r2_xaaa_Lvvo("I,a,e,i");

            // sigmar2_xabb_Lvvo += -0.500000 <j,i||a,n>_abab r1_a(e) t2_abab(a,f,j,i) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 1 | mem: o3v2L1: 1, o1v2L1: 2, 
            sigmar2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovo"]("j,i,a,n") * r1_xa_Lv("I,e") * t2_abab_vvoo("a,f,j,i");

            // sigmar2_xabb_Lvvo += -0.500000 <i,j||a,n>_abab r1_a(e) t2_abab(a,f,i,j) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 1 | mem: o3v2L1: 1, o1v2L1: 2, 
            sigmar2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovo"]("i,j,a,n") * r1_xa_Lv("I,e") * t2_abab_vvoo("a,f,i,j");

            // sigmar2_xabb_Lvvo += -0.500000 <j,i||a,n>_bbbb r1_a(e) t2_bbbb(a,f,j,i) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 1 | mem: o3v2L1: 1, o1v2L1: 2, 
            sigmar2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovo"]("j,i,a,n") * r1_xa_Lv("I,e") * t2_bbbb_vvoo("a,f,j,i");

            // sigmar2_xabb_Lvvo += 0.500000 <j,i||a,n>_abab r1_a(a) t2_abab(e,f,j,i) 
            // flops: o3v2L1: 1, o3v1L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o3v0L1: 1, 
            sigmar2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovo"]("j,i,a,n") * r1_xa_Lv("I,a") * t2_abab_vvoo("e,f,j,i");

            // sigmar2_xabb_Lvvo += 0.500000 <i,j||a,n>_abab r1_a(a) t2_abab(e,f,i,j) 
            // flops: o3v2L1: 1, o3v1L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o3v0L1: 1, 
            sigmar2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovo"]("i,j,a,n") * r1_xa_Lv("I,a") * t2_abab_vvoo("e,f,i,j");

            // sigmar2_xabb_Lvvo += 0.500000 <i,f||a,b>_abab r1_a(e) t2_abab(a,b,i,n) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmar2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["baab_vovv"]("f,i,a,b") * r1_xa_Lv("I,e") * t2_abab_vvoo("a,b,i,n");

            // sigmar2_xabb_Lvvo += 0.500000 <i,f||b,a>_abab r1_a(e) t2_abab(b,a,i,n) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmar2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["baab_vovv"]("f,i,b,a") * r1_xa_Lv("I,e") * t2_abab_vvoo("b,a,i,n");

            // sigmar2_xabb_Lvvo += -0.500000 <i,f||a,b>_bbbb r1_a(e) t2_bbbb(a,b,n,i) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmar2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["bbbb_vovv"]("f,i,a,b") * r1_xa_Lv("I,e") * t2_bbbb_vvoo("a,b,n,i");

            // sigmar2_xabb_Lvvo += 1.000000 <i,e||a,b>_aaaa r1_a(b) t2_abab(a,f,i,n) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmar2_xabb_Lvvo("I,e,f,n") -= V_blks_["aaaa_vovv"]("e,i,a,b") * r1_xa_Lv("I,b") * t2_abab_vvoo("a,f,i,n");

            // sigmar2_xabb_Lvvo += -1.000000 <e,i||b,a>_abab r1_a(b) t2_bbbb(a,f,n,i) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmar2_xabb_Lvvo("I,e,f,n") -= V_blks_["abab_vovv"]("e,i,b,a") * r1_xa_Lv("I,b") * t2_bbbb_vvoo("a,f,n,i");

            // sigmar2_xabb_Lvvo += -1.000000 <i,f||b,a>_abab r1_a(b) t2_abab(e,a,i,n) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmar2_xabb_Lvvo("I,e,f,n") += V_blks_["baab_vovv"]("f,i,b,a") * r1_xa_Lv("I,b") * t2_abab_vvoo("e,a,i,n");

            // sigmar2_xabb_Lvvo += 1.000000 <j,i||a,b>_aaaa r2_aaa(b,e,j) t2_abab(a,f,i,n) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmar2_xabb_Lvvo("I,e,f,n") += V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaa_Lvvo("I,b,e,j") * t2_abab_vvoo("a,f,i,n");

            // sigmar2_xabb_Lvvo += 1.000000 <j,i||b,a>_abab r2_aaa(b,e,j) t2_bbbb(a,f,n,i) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmar2_xabb_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("j,i,b,a") * r2_xaaa_Lvvo("I,b,e,j") * t2_bbbb_vvoo("a,f,n,i");

            // sigmar2_xabb_Lvvo += 0.500000 <j,i||a,b>_aaaa r2_aaa(a,b,j) t2_abab(e,f,i,n) 
            // flops: o2v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmar2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaa_Lvvo("I,a,b,j") * t2_abab_vvoo("e,f,i,n");
        }




/// ****** sigmar2_baa(e,f,n) ****** ///



        {

            // sigmar2_xbaa_Lvvo += 1.000000 f_aa(f,n) r1_b(e) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += F_blks_["aa_vo"]("f,n") * r1_xb_Lv("I,e");

            // sigmar2_xbaa_Lvvo += -1.000000 f_aa(i,a) r1_b(e) t2_aaaa(a,f,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            sigmar2_xbaa_Lvvo("I,e,f,n") -= F_blks_["aa_ov"]("i,a") * r1_xb_Lv("I,e") * t2_aaaa_vvoo("a,f,n,i");

            // sigmar2_xbaa_Lvvo += 1.000000 f_bb(i,a) r1_b(e) t2_abab(f,a,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += F_blks_["bb_ov"]("i,a") * r1_xb_Lv("I,e") * t2_abab_vvoo("f,a,n,i");

            // sigmar2_xbaa_Lvvo += -1.000000 f_bb(i,a) r1_b(a) t2_abab(f,e,n,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o1v1L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmar2_xbaa_Lvvo("I,e,f,n") -= F_blks_["bb_ov"]("i,a") * r1_xb_Lv("I,a") * t2_abab_vvoo("f,e,n,i");
        }

        if (include_u1_) {

            // sigmar2_xbaa_Lvvo += -1.000000 d-_aa(f,n) s1_b(e) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            sigmar2_xbaa_Lvvo("I,e,f,n") -= dp_aa_vo("f,n") * s1_xb_Lv("I,e");

            // sigmar2_xbaa_Lvvo += 1.000000 d-_aa(i,a) t2_aaaa(a,f,n,i) s1_b(e) 
            // flops: o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += dp_aa_ov("i,a") * t2_aaaa_vvoo("a,f,n,i") * s1_xb_Lv("I,e");

            // sigmar2_xbaa_Lvvo += -1.000000 d-_bb(i,a) t2_abab(f,a,n,i) s1_b(e) 
            // flops: o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmar2_xbaa_Lvvo("I,e,f,n") -= dp_bb_ov("i,a") * t2_abab_vvoo("f,a,n,i") * s1_xb_Lv("I,e");

            // sigmar2_xbaa_Lvvo += 1.000000 d-_bb(i,a) t2_abab(f,e,n,i) s1_b(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += dp_bb_ov("i,a") * t2_abab_vvoo("f,e,n,i") * s1_xb_Lv("I,a");
        }

        if (include_u0_) {

            // sigmar2_xbaa_Lvvo += -1.000000 d-_aa(f,n) r1_b(e) u0 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, o1v2: 1, 
            sigmar2_xbaa_Lvvo("I,e,f,n") -= dp_aa_vo("f,n") * r1_xb_Lv("I,e") * u0;
        }

        if (include_u1_) {

            // sigmar2_xbaa_Lvvo += -1.000000 d-_aa(i,i) r1_b(e) u1_aa(f,n) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            sigmar2_xbaa_Lvvo("I,e,f,n") -= r1_xb_Lv("I,e") * scalar7 * u1_aa_vo("f,n");

            // sigmar2_xbaa_Lvvo += -1.000000 d-_bb(i,i) r1_b(e) u1_aa(f,n) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            sigmar2_xbaa_Lvvo("I,e,f,n") -= r1_xb_Lv("I,e") * scalar8 * u1_aa_vo("f,n");

            // sigmar2_xbaa_Lvvo += 1.000000 d-_aa(i,n) r1_b(e) u1_aa(f,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o2v1L1: 1 | mem: o1v2L1: 2, o2v1L1: 1, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += dp_aa_oo("i,n") * r1_xb_Lv("I,e") * u1_aa_vo("f,i");

            // sigmar2_xbaa_Lvvo += -1.000000 d-_aa(f,a) r1_b(e) u1_aa(a,n) 
            // flops: o1v3L1: 1, o0v3L1: 1, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, 
            sigmar2_xbaa_Lvvo("I,e,f,n") -= dp_aa_vv("f,a") * r1_xb_Lv("I,e") * u1_aa_vo("a,n");

            // sigmar2_xbaa_Lvvo += -1.000000 d-_bb(e,a) r1_b(a) u1_aa(f,n) 
            // flops: o1v2L1: 2, o0v2L1: 1 | mem: o1v2L1: 2, o0v1L1: 1, 
            sigmar2_xbaa_Lvvo("I,e,f,n") -= dp_bb_vv("e,a") * r1_xb_Lv("I,a") * u1_aa_vo("f,n");

            // sigmar2_xbaa_Lvvo += 1.000000 d-_bb(i,a) r2_bbb(a,e,i) u1_aa(f,n) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, o0v1L1: 1, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += dp_bb_ov("i,a") * r2_xbbb_Lvvo("I,a,e,i") * u1_aa_vo("f,n");
        }

        if (include_u2_) {

            // sigmar2_xbaa_Lvvo += 1.000000 d-_aa(i,a) r1_b(e) u2_aaaa(a,f,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += dp_aa_ov("i,a") * r1_xb_Lv("I,e") * u2_aaaa_vvoo("a,f,n,i");

            // sigmar2_xbaa_Lvvo += -1.000000 d-_bb(i,a) r1_b(e) u2_abab(f,a,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            sigmar2_xbaa_Lvvo("I,e,f,n") -= dp_bb_ov("i,a") * r1_xb_Lv("I,e") * u2_abab_vvoo("f,a,n,i");

            // sigmar2_xbaa_Lvvo += 1.000000 d-_bb(i,a) r1_b(a) u2_abab(f,e,n,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o1v1L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += dp_bb_ov("i,a") * r1_xb_Lv("I,a") * u2_abab_vvoo("f,e,n,i");
        }

        if (include_u0_) {

            // sigmar2_xbaa_Lvvo += 1.000000 d-_aa(i,a) r1_b(e) t2_aaaa(a,f,n,i) u0 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, o1v2: 1, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += dp_aa_ov("i,a") * r1_xb_Lv("I,e") * t2_aaaa_vvoo("a,f,n,i") * u0;

            // sigmar2_xbaa_Lvvo += -1.000000 d-_bb(i,a) r1_b(e) t2_abab(f,a,n,i) u0 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, o1v2: 1, 
            sigmar2_xbaa_Lvvo("I,e,f,n") -= dp_bb_ov("i,a") * r1_xb_Lv("I,e") * t2_abab_vvoo("f,a,n,i") * u0;

            // sigmar2_xbaa_Lvvo += 1.000000 d-_bb(i,a) r1_b(a) t2_abab(f,e,n,i) u0 
            // flops: o2v2L1: 1, o1v2L1: 2, o1v1L1: 1 | mem: o1v2L1: 2, o1v2: 1, o1v0L1: 1, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += dp_bb_ov("i,a") * r1_xb_Lv("I,a") * t2_abab_vvoo("f,e,n,i") * u0;
        }

        {

            // sigmar2_xbaa_Lvvo += 1.000000 <f,e||n,a>_abab r1_b(a) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmar2_xbaa_Lvvo("I,e,f,n") -= V_blks_["abba_vvvo"]("f,e,a,n") * r1_xb_Lv("I,a");

            // sigmar2_xbaa_Lvvo += -1.000000 <f,i||n,a>_abab r2_bbb(a,e,i) 
            // flops: o2v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += V_blks_["abba_vovo"]("f,i,a,n") * r2_xbbb_Lvvo("I,a,e,i");

            // sigmar2_xbaa_Lvvo += -0.500000 <j,i||a,n>_aaaa r1_b(e) t2_aaaa(a,f,j,i) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 1 | mem: o3v2L1: 1, o1v2L1: 2, 
            sigmar2_xbaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovo"]("j,i,a,n") * r1_xb_Lv("I,e") * t2_aaaa_vvoo("a,f,j,i");

            // sigmar2_xbaa_Lvvo += -0.500000 <j,i||n,a>_abab r1_b(e) t2_abab(f,a,j,i) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 1 | mem: o3v2L1: 1, o1v2L1: 2, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abba_oovo"]("j,i,a,n") * r1_xb_Lv("I,e") * t2_abab_vvoo("f,a,j,i");

            // sigmar2_xbaa_Lvvo += -0.500000 <i,j||n,a>_abab r1_b(e) t2_abab(f,a,i,j) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 1 | mem: o3v2L1: 1, o1v2L1: 2, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abba_oovo"]("i,j,a,n") * r1_xb_Lv("I,e") * t2_abab_vvoo("f,a,i,j");

            // sigmar2_xbaa_Lvvo += 0.500000 <j,i||n,a>_abab r1_b(a) t2_abab(f,e,j,i) 
            // flops: o3v2L1: 1, o3v1L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o3v0L1: 1, 
            sigmar2_xbaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abba_oovo"]("j,i,a,n") * r1_xb_Lv("I,a") * t2_abab_vvoo("f,e,j,i");

            // sigmar2_xbaa_Lvvo += 0.500000 <i,j||n,a>_abab r1_b(a) t2_abab(f,e,i,j) 
            // flops: o3v2L1: 1, o3v1L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o3v0L1: 1, 
            sigmar2_xbaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abba_oovo"]("i,j,a,n") * r1_xb_Lv("I,a") * t2_abab_vvoo("f,e,i,j");

            // sigmar2_xbaa_Lvvo += -0.500000 <i,f||a,b>_aaaa r1_b(e) t2_aaaa(a,b,n,i) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["aaaa_vovv"]("f,i,a,b") * r1_xb_Lv("I,e") * t2_aaaa_vvoo("a,b,n,i");

            // sigmar2_xbaa_Lvvo += 0.500000 <f,i||a,b>_abab r1_b(e) t2_abab(a,b,n,i) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_vovv"]("f,i,a,b") * r1_xb_Lv("I,e") * t2_abab_vvoo("a,b,n,i");

            // sigmar2_xbaa_Lvvo += 0.500000 <f,i||b,a>_abab r1_b(e) t2_abab(b,a,n,i) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_vovv"]("f,i,b,a") * r1_xb_Lv("I,e") * t2_abab_vvoo("b,a,n,i");

            // sigmar2_xbaa_Lvvo += -1.000000 <i,e||a,b>_abab r1_b(b) t2_aaaa(a,f,n,i) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += V_blks_["baab_vovv"]("e,i,a,b") * r1_xb_Lv("I,b") * t2_aaaa_vvoo("a,f,n,i");

            // sigmar2_xbaa_Lvvo += 1.000000 <i,e||a,b>_bbbb r1_b(b) t2_abab(f,a,n,i) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmar2_xbaa_Lvvo("I,e,f,n") -= V_blks_["bbbb_vovv"]("e,i,a,b") * r1_xb_Lv("I,b") * t2_abab_vvoo("f,a,n,i");

            // sigmar2_xbaa_Lvvo += -1.000000 <f,i||a,b>_abab r1_b(b) t2_abab(a,e,n,i) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmar2_xbaa_Lvvo("I,e,f,n") -= V_blks_["abab_vovv"]("f,i,a,b") * r1_xb_Lv("I,b") * t2_abab_vvoo("a,e,n,i");

            // sigmar2_xbaa_Lvvo += 1.000000 <i,j||a,b>_abab r2_bbb(b,e,j) t2_aaaa(a,f,n,i) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("i,j,a,b") * r2_xbbb_Lvvo("I,b,e,j") * t2_aaaa_vvoo("a,f,n,i");

            // sigmar2_xbaa_Lvvo += 1.000000 <j,i||a,b>_bbbb r2_bbb(b,e,j) t2_abab(f,a,n,i) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += V_blks_["bbbb_oovv"]("j,i,a,b") * r2_xbbb_Lvvo("I,b,e,j") * t2_abab_vvoo("f,a,n,i");

            // sigmar2_xbaa_Lvvo += 0.500000 <j,i||a,b>_bbbb r2_bbb(a,b,j) t2_abab(f,e,n,i) 
            // flops: o2v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmar2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r2_xbbb_Lvvo("I,a,b,j") * t2_abab_vvoo("f,e,n,i");
        }




/// ****** sigmar2_bbb(e,f,n) ****** ///



        {

            // sigmar2_xbbb_Lvvo += -1.000000 P(e,f) f_bb(e,n) r1_b(f) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = F_blks_["bb_vo"]("e,n") * r1_xb_Lv("I,f");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += 1.000000 f_aa(i,i) r2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") += scalar5 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmar2_xbbb_Lvvo += 1.000000 f_bb(i,i) r2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") += scalar6 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmar2_xbbb_Lvvo += -1.000000 f_bb(i,n) r2_bbb(e,f,i) 
            // flops: o2v2L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmar2_xbbb_Lvvo("I,e,f,n") -= F_blks_["bb_oo"]("i,n") * r2_xbbb_Lvvo("I,e,f,i");

            // sigmar2_xbbb_Lvvo += 1.000000 P(e,f) f_bb(e,a) r2_bbb(a,f,n) 
            // flops: o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = F_blks_["bb_vv"]("e,a") * r2_xbbb_Lvvo("I,a,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += -1.000000 P(e,f) f_aa(i,a) r1_b(f) t2_abab(a,e,i,n) 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = F_blks_["aa_ov"]("i,a") * r1_xb_Lv("I,f") * t2_abab_vvoo("a,e,i,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += 1.000000 P(e,f) f_bb(i,a) r1_b(f) t2_bbbb(a,e,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = F_blks_["bb_ov"]("i,a") * r1_xb_Lv("I,f") * t2_bbbb_vvoo("a,e,n,i");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += 1.000000 f_bb(i,a) r1_b(a) t2_bbbb(e,f,n,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o1v1L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") += F_blks_["bb_ov"]("i,a") * r1_xb_Lv("I,a") * t2_bbbb_vvoo("e,f,n,i");
        }

        if (include_u1_) {

            // sigmar2_xbbb_Lvvo += 1.000000 P(e,f) d-_bb(e,n) s1_b(f) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_vo("e,n") * s1_xb_Lv("I,f");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmar2_xbbb_Lvvo += -1.000000 d-_aa(i,i) s2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") -= scalar7 * s2_xbbb_Lvvo("I,e,f,n");

            // sigmar2_xbbb_Lvvo += -1.000000 d-_bb(i,i) s2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") -= scalar8 * s2_xbbb_Lvvo("I,e,f,n");

            // sigmar2_xbbb_Lvvo += 1.000000 d-_bb(i,n) s2_bbb(e,f,i) 
            // flops: o2v2L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmar2_xbbb_Lvvo("I,e,f,n") += dp_bb_oo("i,n") * s2_xbbb_Lvvo("I,e,f,i");

            // sigmar2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(e,a) s2_bbb(a,f,n) 
            // flops: o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_vv("e,a") * s2_xbbb_Lvvo("I,a,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u1_) {

            // sigmar2_xbbb_Lvvo += 1.000000 P(e,f) d-_aa(i,a) t2_abab(a,e,i,n) s1_b(f) 
            // flops: o1v2L1: 3, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_aa_ov("i,a") * t2_abab_vvoo("a,e,i,n") * s1_xb_Lv("I,f");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(i,a) t2_bbbb(a,e,n,i) s1_b(f) 
            // flops: o1v2L1: 3, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("i,a") * t2_bbbb_vvoo("a,e,n,i") * s1_xb_Lv("I,f");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += -1.000000 d-_bb(i,a) t2_bbbb(e,f,n,i) s1_b(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") -= dp_bb_ov("i,a") * t2_bbbb_vvoo("e,f,n,i") * s1_xb_Lv("I,a");
        }

        if (include_u0_) {

            // sigmar2_xbbb_Lvvo += 1.000000 P(e,f) d-_bb(e,n) r1_b(f) u0 
            // flops: o1v2L1: 4 | mem: o1v2L1: 2, o1v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_vo("e,n") * r1_xb_Lv("I,f") * u0;
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += -1.000000 d-_aa(i,i) r2_bbb(e,f,n) u0 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") -= scalar7 * u0 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmar2_xbbb_Lvvo += -1.000000 d-_bb(i,i) r2_bbb(e,f,n) u0 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") -= scalar8 * u0 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmar2_xbbb_Lvvo += 1.000000 d-_bb(i,n) r2_bbb(e,f,i) u0 
            // flops: o2v2L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v2: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") += dp_bb_oo("i,n") * r2_xbbb_Lvvo("I,e,f,i") * u0;

            // sigmar2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(e,a) r2_bbb(a,f,n) u0 
            // flops: o1v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_vv("e,a") * r2_xbbb_Lvvo("I,a,f,n") * u0;
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u1_) {

            // sigmar2_xbbb_Lvvo += 1.000000 P(e,f) d-_aa(i,i) r1_b(f) u1_bb(e,n) 
            // flops: o1v2L1: 2, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = r1_xb_Lv("I,f") * scalar7 * u1_bb_vo("e,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += 1.000000 P(e,f) d-_bb(i,i) r1_b(f) u1_bb(e,n) 
            // flops: o1v2L1: 2, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = r1_xb_Lv("I,f") * scalar8 * u1_bb_vo("e,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(i,n) r1_b(f) u1_bb(e,i) 
            // flops: o2v2L1: 1, o1v2L1: 2, o2v1L1: 1 | mem: o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_oo("i,n") * r1_xb_Lv("I,f") * u1_bb_vo("e,i");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += 1.000000 P(e,f) d-_bb(e,a) r1_b(f) u1_bb(a,n) 
            // flops: o1v3L1: 1, o0v3L1: 1, o1v2L1: 2 | mem: o0v3L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_vv("e,a") * r1_xb_Lv("I,f") * u1_bb_vo("a,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(e,a) r1_b(a) u1_bb(f,n) 
            // flops: o1v2L1: 3, o0v2L1: 1 | mem: o1v2L1: 2, o0v1L1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_vv("e,a") * r1_xb_Lv("I,a") * u1_bb_vo("f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += -1.000000 d-_aa(i,a) r2_bbb(e,f,n) u1_aa(a,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") -= scalar0 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmar2_xbbb_Lvvo += -1.000000 d-_bb(i,a) r2_bbb(e,f,n) u1_bb(a,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") -= scalar1 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmar2_xbbb_Lvvo += 1.000000 d-_bb(i,a) r2_bbb(e,f,i) u1_bb(a,n) 
            // flops: o1v3L1: 2, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, 
            sigmar2_xbbb_Lvvo("I,e,f,n") += dp_bb_ov("i,a") * r2_xbbb_Lvvo("I,e,f,i") * u1_bb_vo("a,n");

            // sigmar2_xbbb_Lvvo += 1.000000 P(e,f) d-_bb(i,a) r2_bbb(a,f,n) u1_bb(e,i) 
            // flops: o2v2L1: 2, o1v2L1: 2 | mem: o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("i,a") * r2_xbbb_Lvvo("I,a,f,n") * u1_bb_vo("e,i");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(i,a) r2_bbb(a,f,i) u1_bb(e,n) 
            // flops: o1v2L1: 4 | mem: o1v2L1: 2, o0v1L1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("i,a") * r2_xbbb_Lvvo("I,a,f,i") * u1_bb_vo("e,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmar2_xbbb_Lvvo += 1.000000 P(e,f) d-_aa(i,a) r1_b(f) u2_abab(a,e,i,n) 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_aa_ov("i,a") * r1_xb_Lv("I,f") * u2_abab_vvoo("a,e,i,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(i,a) r1_b(f) u2_bbbb(a,e,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("i,a") * r1_xb_Lv("I,f") * u2_bbbb_vvoo("a,e,n,i");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += -1.000000 d-_bb(i,a) r1_b(a) u2_bbbb(e,f,n,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o1v1L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") -= dp_bb_ov("i,a") * r1_xb_Lv("I,a") * u2_bbbb_vvoo("e,f,n,i");
        }

        if (include_u0_) {

            // sigmar2_xbbb_Lvvo += 1.000000 P(e,f) d-_aa(i,a) r1_b(f) t2_abab(a,e,i,n) u0 
            // flops: o2v3L1: 1, o1v2L1: 4 | mem: o1v2L1: 3, o1v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_aa_ov("i,a") * r1_xb_Lv("I,f") * t2_abab_vvoo("a,e,i,n") * u0;
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(i,a) r1_b(f) t2_bbbb(a,e,n,i) u0 
            // flops: o2v3L1: 1, o1v2L1: 4 | mem: o1v2L1: 3, o1v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("i,a") * r1_xb_Lv("I,f") * t2_bbbb_vvoo("a,e,n,i") * u0;
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += -1.000000 d-_bb(i,a) r1_b(a) t2_bbbb(e,f,n,i) u0 
            // flops: o2v2L1: 1, o1v2L1: 2, o1v1L1: 1 | mem: o1v2L1: 2, o1v2: 1, o1v0L1: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") -= dp_bb_ov("i,a") * r1_xb_Lv("I,a") * t2_bbbb_vvoo("e,f,n,i") * u0;
        }

        {

            // sigmar2_xbbb_Lvvo += -0.500000 <j,i||j,i>_aaaa r2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * scalar9 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmar2_xbbb_Lvvo += -0.500000 <j,i||j,i>_abab r2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * scalar10 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmar2_xbbb_Lvvo += -0.500000 <i,j||i,j>_abab r2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * scalar10 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmar2_xbbb_Lvvo += -0.500000 <j,i||j,i>_bbbb r2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * scalar11 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmar2_xbbb_Lvvo += 1.000000 <e,f||a,n>_bbbb r1_b(a) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmar2_xbbb_Lvvo("I,e,f,n") += V_blks_["bbbb_vvvo"]("e,f,a,n") * r1_xb_Lv("I,a");

            // sigmar2_xbbb_Lvvo += 1.000000 P(e,f) <i,e||a,n>_bbbb r2_bbb(a,f,i) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_vovo"]("e,i,a,n") * r2_xbbb_Lvvo("I,a,f,i");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += 0.500000 <e,f||a,b>_bbbb r2_bbb(a,b,n) 
            // flops: o1v4L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmar2_xbbb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["bbbb_vvvv"]("e,f,a,b") * r2_xbbb_Lvvo("I,a,b,n");

            // sigmar2_xbbb_Lvvo += 0.500000 P(e,f) <j,i||a,n>_abab r1_b(f) t2_abab(a,e,j,i) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovo"]("j,i,a,n") * r1_xb_Lv("I,f") * t2_abab_vvoo("a,e,j,i");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += 0.500000 P(e,f) <i,j||a,n>_abab r1_b(f) t2_abab(a,e,i,j) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovo"]("i,j,a,n") * r1_xb_Lv("I,f") * t2_abab_vvoo("a,e,i,j");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += 0.500000 P(e,f) <j,i||a,n>_bbbb r1_b(f) t2_bbbb(a,e,j,i) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["bbbb_oovo"]("j,i,a,n") * r1_xb_Lv("I,f") * t2_bbbb_vvoo("a,e,j,i");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += 0.500000 <j,i||a,n>_bbbb r1_b(a) t2_bbbb(e,f,j,i) 
            // flops: o3v2L1: 1, o3v1L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o3v0L1: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["bbbb_oovo"]("j,i,a,n") * r1_xb_Lv("I,a") * t2_bbbb_vvoo("e,f,j,i");

            // sigmar2_xbbb_Lvvo += -0.500000 P(e,f) <i,e||a,b>_abab r1_b(f) t2_abab(a,b,i,n) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 2 | mem: o1v4L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["baab_vovv"]("e,i,a,b") * r1_xb_Lv("I,f") * t2_abab_vvoo("a,b,i,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += -0.500000 P(e,f) <i,e||b,a>_abab r1_b(f) t2_abab(b,a,i,n) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 2 | mem: o1v4L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["baab_vovv"]("e,i,b,a") * r1_xb_Lv("I,f") * t2_abab_vvoo("b,a,i,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += 0.500000 P(e,f) <i,e||a,b>_bbbb r1_b(f) t2_bbbb(a,b,n,i) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 2 | mem: o1v4L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["bbbb_vovv"]("e,i,a,b") * r1_xb_Lv("I,f") * t2_bbbb_vvoo("a,b,n,i");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += 1.000000 P(e,f) <i,e||a,b>_abab r1_b(b) t2_abab(a,f,i,n) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["baab_vovv"]("e,i,a,b") * r1_xb_Lv("I,b") * t2_abab_vvoo("a,f,i,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += -1.000000 P(e,f) <i,e||a,b>_bbbb r1_b(b) t2_bbbb(a,f,n,i) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_vovv"]("e,i,a,b") * r1_xb_Lv("I,b") * t2_bbbb_vvoo("a,f,n,i");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += 0.250000 <j,i||a,b>_aaaa r2_bbb(e,f,n) t2_aaaa(a,b,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar2 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmar2_xbbb_Lvvo += 0.250000 <j,i||a,b>_abab r2_bbb(e,f,n) t2_abab(a,b,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar3 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmar2_xbbb_Lvvo += 0.250000 <i,j||a,b>_abab r2_bbb(e,f,n) t2_abab(a,b,i,j) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar3 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmar2_xbbb_Lvvo += 0.250000 <j,i||b,a>_abab r2_bbb(e,f,n) t2_abab(b,a,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar3 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmar2_xbbb_Lvvo += 0.250000 <i,j||b,a>_abab r2_bbb(e,f,n) t2_abab(b,a,i,j) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar3 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmar2_xbbb_Lvvo += 0.250000 <j,i||a,b>_bbbb r2_bbb(e,f,n) t2_bbbb(a,b,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar4 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmar2_xbbb_Lvvo += -0.500000 <i,j||a,b>_abab r2_bbb(e,f,j) t2_abab(a,b,i,n) 
            // flops: o2v4L1: 2, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmar2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("i,j,a,b") * r2_xbbb_Lvvo("I,e,f,j") * t2_abab_vvoo("a,b,i,n");

            // sigmar2_xbbb_Lvvo += -0.500000 <i,j||b,a>_abab r2_bbb(e,f,j) t2_abab(b,a,i,n) 
            // flops: o2v4L1: 2, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmar2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("i,j,b,a") * r2_xbbb_Lvvo("I,e,f,j") * t2_abab_vvoo("b,a,i,n");

            // sigmar2_xbbb_Lvvo += -0.500000 <j,i||a,b>_bbbb r2_bbb(e,f,j) t2_bbbb(a,b,n,i) 
            // flops: o2v4L1: 2, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmar2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r2_xbbb_Lvvo("I,e,f,j") * t2_bbbb_vvoo("a,b,n,i");

            // sigmar2_xbbb_Lvvo += -0.500000 P(e,f) <j,i||a,b>_abab r2_bbb(b,f,n) t2_abab(a,e,j,i) 
            // flops: o3v3L1: 2, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("j,i,a,b") * r2_xbbb_Lvvo("I,b,f,n") * t2_abab_vvoo("a,e,j,i");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += -0.500000 P(e,f) <i,j||a,b>_abab r2_bbb(b,f,n) t2_abab(a,e,i,j) 
            // flops: o3v3L1: 2, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("i,j,a,b") * r2_xbbb_Lvvo("I,b,f,n") * t2_abab_vvoo("a,e,i,j");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += -0.500000 P(e,f) <j,i||a,b>_bbbb r2_bbb(b,f,n) t2_bbbb(a,e,j,i) 
            // flops: o3v3L1: 2, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r2_xbbb_Lvvo("I,b,f,n") * t2_bbbb_vvoo("a,e,j,i");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += 1.000000 P(e,f) <i,j||a,b>_abab r2_bbb(b,f,j) t2_abab(a,e,i,n) 
            // flops: o2v3L1: 2, o1v2L1: 2 | mem: o1v2L1: 3, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("i,j,a,b") * r2_xbbb_Lvvo("I,b,f,j") * t2_abab_vvoo("a,e,i,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += 1.000000 P(e,f) <j,i||a,b>_bbbb r2_bbb(b,f,j) t2_bbbb(a,e,n,i) 
            // flops: o2v3L1: 2, o1v2L1: 2 | mem: o1v2L1: 3, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_oovv"]("j,i,a,b") * r2_xbbb_Lvvo("I,b,f,j") * t2_bbbb_vvoo("a,e,n,i");
            sigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmar2_xbbb_Lvvo += 0.250000 <j,i||a,b>_bbbb r2_bbb(a,b,n) t2_bbbb(e,f,j,i) 
            // flops: o3v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o3v0L1: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") += 0.250000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r2_xbbb_Lvvo("I,a,b,n") * t2_bbbb_vvoo("e,f,j,i");

            // sigmar2_xbbb_Lvvo += -0.500000 <j,i||a,b>_bbbb r2_bbb(a,b,j) t2_bbbb(e,f,n,i) 
            // flops: o2v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmar2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r2_xbbb_Lvvo("I,a,b,j") * t2_bbbb_vvoo("e,f,n,i");
        }




/// ****** sigmal1_a(e) ****** ///



        if (include_u0_ && include_u1_) {

            // sigmal1_xa_Lv += 1.000000 u0 m1_a(e) w0 
            // flops: o0v1L1: 2, o0v1: 1 | mem: o0v1L1: 1, o0v1: 2, 
            sigmal1_xa_Lv("I,e") += u0 * m1_xa_Lv("I,e") * w0;
        }

        if (include_u1_ && include_u2_) {

            // sigmal1_xa_Lv += -1.000000 u1_aa(a,i) m2_aaa(i,a,e) w0 
            // flops: o1v2L1: 1, o0v1L1: 2 | mem: o0v1L1: 2, o0v1: 1, 
            sigmal1_xa_Lv("I,e") -= u1_aa_vo("a,i") * m2_xaaa_Lovv("I,i,a,e") * w0;

            // sigmal1_xa_Lv += -1.000000 u1_bb(a,i) m2_bba(i,a,e) w0 
            // flops: o1v2L1: 1, o0v1L1: 2 | mem: o0v1L1: 2, o0v1: 1, 
            sigmal1_xa_Lv("I,e") -= u1_bb_vo("a,i") * m2_xbba_Lovv("I,i,a,e") * w0;
        }

        {

            // sigmal1_xa_Lv += 1.000000 f_aa(i,i) l1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += scalar5 * l1_xa_Lv("I,e");

            // sigmal1_xa_Lv += 1.000000 f_bb(i,i) l1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += scalar6 * l1_xa_Lv("I,e");

            // sigmal1_xa_Lv += 1.000000 f_aa(a,e) l1_a(a) 
            // flops: o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += F_blks_["aa_vv"]("a,e") * l1_xa_Lv("I,a");

            // sigmal1_xa_Lv += -1.000000 f_aa(a,i) l2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= F_blks_["aa_vo"]("a,i") * l2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -1.000000 f_bb(a,i) l2_bba(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= F_blks_["bb_vo"]("a,i") * l2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 f_aa(j,b) l2_aaa(i,a,e) t2_aaaa(b,a,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += F_blks_["aa_ov"]("j,b") * l2_xaaa_Lovv("I,i,a,e") * t2_aaaa_vvoo("b,a,i,j");

            // sigmal1_xa_Lv += -1.000000 f_bb(j,b) l2_aaa(i,a,e) t2_abab(a,b,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= F_blks_["bb_ov"]("j,b") * l2_xaaa_Lovv("I,i,a,e") * t2_abab_vvoo("a,b,i,j");

            // sigmal1_xa_Lv += -1.000000 f_aa(j,b) l2_bba(i,a,e) t2_abab(b,a,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= F_blks_["aa_ov"]("j,b") * l2_xbba_Lovv("I,i,a,e") * t2_abab_vvoo("b,a,j,i");

            // sigmal1_xa_Lv += 1.000000 f_bb(j,b) l2_bba(i,a,e) t2_bbbb(b,a,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += F_blks_["bb_ov"]("j,b") * l2_xbba_Lovv("I,i,a,e") * t2_bbbb_vvoo("b,a,i,j");

            // sigmal1_xa_Lv += 0.500000 f_aa(j,e) l2_aaa(i,b,a) t2_aaaa(b,a,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += 0.500000 * F_blks_["aa_ov"]("j,e") * l2_xaaa_Lovv("I,i,b,a") * t2_aaaa_vvoo("b,a,i,j");

            // sigmal1_xa_Lv += 0.500000 f_aa(j,e) l2_bba(i,b,a) t2_abab(a,b,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += 0.500000 * F_blks_["aa_ov"]("j,e") * l2_xbba_Lovv("I,i,b,a") * t2_abab_vvoo("a,b,j,i");
        }

        if (include_u1_) {

            // sigmal1_xa_Lv += 1.000000 f_aa(i,a) u1_aa(a,i) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += scalar12 * m1_xa_Lv("I,e");

            // sigmal1_xa_Lv += 1.000000 f_bb(i,a) u1_bb(a,i) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += scalar13 * m1_xa_Lv("I,e");

            // sigmal1_xa_Lv += -1.000000 f_aa(i,e) u1_aa(a,i) m1_a(a) 
            // flops: o0v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmal1_xa_Lv("I,e") -= F_blks_["aa_ov"]("i,e") * u1_aa_vo("a,i") * m1_xa_Lv("I,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmal1_xa_Lv += 1.000000 f_aa(j,i) u1_aa(a,j) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += F_blks_["aa_oo"]("j,i") * u1_aa_vo("a,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 f_bb(j,i) u1_bb(a,j) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += F_blks_["bb_oo"]("j,i") * u1_bb_vo("a,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -1.000000 f_aa(a,b) u1_aa(b,i) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") -= F_blks_["aa_vv"]("a,b") * u1_aa_vo("b,i") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -1.000000 f_bb(a,b) u1_bb(b,i) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") -= F_blks_["bb_vv"]("a,b") * u1_bb_vo("b,i") * m2_xbba_Lovv("I,i,a,e");
        }

        if (include_u2_) {

            // sigmal1_xa_Lv += 1.000000 f_aa(j,b) u2_aaaa(b,a,i,j) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += F_blks_["aa_ov"]("j,b") * u2_aaaa_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -1.000000 f_aa(j,b) u2_abab(b,a,j,i) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") -= F_blks_["aa_ov"]("j,b") * u2_abab_vvoo("b,a,j,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -1.000000 f_bb(j,b) u2_abab(a,b,i,j) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") -= F_blks_["bb_ov"]("j,b") * u2_abab_vvoo("a,b,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 f_bb(j,b) u2_bbbb(b,a,i,j) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += F_blks_["bb_ov"]("j,b") * u2_bbbb_vvoo("b,a,i,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.500000 f_aa(j,e) u2_aaaa(b,a,i,j) m2_aaa(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += 0.500000 * F_blks_["aa_ov"]("j,e") * u2_aaaa_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += 0.500000 f_aa(j,e) u2_abab(a,b,j,i) m2_bba(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += 0.500000 * F_blks_["aa_ov"]("j,e") * u2_abab_vvoo("a,b,j,i") * m2_xbba_Lovv("I,i,b,a");
        }

        if (include_u1_) {

            // sigmal1_xa_Lv += -1.000000 d+_aa(i,i) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") -= scalar17 * m1_xa_Lv("I,e");

            // sigmal1_xa_Lv += -1.000000 d+_bb(i,i) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") -= scalar18 * m1_xa_Lv("I,e");

            // sigmal1_xa_Lv += -1.000000 d+_aa(a,e) m1_a(a) 
            // flops: o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= dp_aa_vv("a,e") * m1_xa_Lv("I,a");
        }

        if (include_u2_) {

            // sigmal1_xa_Lv += 1.000000 d+_aa(a,i) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += dp_aa_vo("a,i") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 d+_bb(a,i) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += dp_bb_vo("a,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -1.000000 d+_aa(j,b) t2_aaaa(b,a,i,j) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") -= dp_aa_ov("j,b") * t2_aaaa_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 d+_aa(j,b) t2_abab(b,a,j,i) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += dp_aa_ov("j,b") * t2_abab_vvoo("b,a,j,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 d+_bb(j,b) t2_abab(a,b,i,j) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += dp_bb_ov("j,b") * t2_abab_vvoo("a,b,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -1.000000 d+_bb(j,b) t2_bbbb(b,a,i,j) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") -= dp_bb_ov("j,b") * t2_bbbb_vvoo("b,a,i,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -0.500000 d+_aa(j,e) t2_aaaa(b,a,i,j) m2_aaa(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * dp_aa_ov("j,e") * t2_aaaa_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += -0.500000 d+_aa(j,e) t2_abab(a,b,j,i) m2_bba(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * dp_aa_ov("j,e") * t2_abab_vvoo("a,b,j,i") * m2_xbba_Lovv("I,i,b,a");
        }

        if (include_u0_) {

            // sigmal1_xa_Lv += -1.000000 d-_aa(i,i) l1_a(e) u0 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmal1_xa_Lv("I,e") -= scalar7 * u0 * l1_xa_Lv("I,e");

            // sigmal1_xa_Lv += -1.000000 d-_bb(i,i) l1_a(e) u0 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmal1_xa_Lv("I,e") -= scalar8 * u0 * l1_xa_Lv("I,e");

            // sigmal1_xa_Lv += -1.000000 d-_aa(a,e) l1_a(a) u0 
            // flops: o0v2L1: 1, o0v1L1: 2 | mem: o0v1L1: 2, o0v1: 1, 
            sigmal1_xa_Lv("I,e") -= dp_aa_vv("a,e") * l1_xa_Lv("I,a") * u0;

            // sigmal1_xa_Lv += 1.000000 d-_aa(a,i) l2_aaa(i,a,e) u0 
            // flops: o1v2L1: 1, o0v1L1: 2 | mem: o0v1L1: 2, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += dp_aa_vo("a,i") * l2_xaaa_Lovv("I,i,a,e") * u0;

            // sigmal1_xa_Lv += 1.000000 d-_bb(a,i) l2_bba(i,a,e) u0 
            // flops: o1v2L1: 1, o0v1L1: 2 | mem: o0v1L1: 2, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += dp_bb_vo("a,i") * l2_xbba_Lovv("I,i,a,e") * u0;
        }

        if (include_u1_) {

            // sigmal1_xa_Lv += -1.000000 d-_aa(i,a) l1_a(e) u1_aa(a,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") -= scalar0 * l1_xa_Lv("I,e");

            // sigmal1_xa_Lv += -1.000000 d-_bb(i,a) l1_a(e) u1_bb(a,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") -= scalar1 * l1_xa_Lv("I,e");

            // sigmal1_xa_Lv += 1.000000 d-_aa(i,e) l1_a(a) u1_aa(a,i) 
            // flops: o1v2L1: 2, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += dp_aa_ov("i,e") * l1_xa_Lv("I,a") * u1_aa_vo("a,i");

            // sigmal1_xa_Lv += 1.000000 d-_aa(i,i) l2_aaa(j,a,e) u1_aa(a,j) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += scalar7 * u1_aa_vo("a,j") * l2_xaaa_Lovv("I,j,a,e");

            // sigmal1_xa_Lv += 1.000000 d-_bb(i,i) l2_aaa(j,a,e) u1_aa(a,j) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += scalar8 * u1_aa_vo("a,j") * l2_xaaa_Lovv("I,j,a,e");

            // sigmal1_xa_Lv += 1.000000 d-_aa(i,i) l2_bba(j,a,e) u1_bb(a,j) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += scalar7 * u1_bb_vo("a,j") * l2_xbba_Lovv("I,j,a,e");

            // sigmal1_xa_Lv += 1.000000 d-_bb(i,i) l2_bba(j,a,e) u1_bb(a,j) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += scalar8 * u1_bb_vo("a,j") * l2_xbba_Lovv("I,j,a,e");

            // sigmal1_xa_Lv += -1.000000 d-_aa(j,i) l2_aaa(i,a,e) u1_aa(a,j) 
            // flops: o2v2L1: 1, o1v2L1: 1, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= dp_aa_oo("j,i") * l2_xaaa_Lovv("I,i,a,e") * u1_aa_vo("a,j");

            // sigmal1_xa_Lv += -1.000000 d-_bb(j,i) l2_bba(i,a,e) u1_bb(a,j) 
            // flops: o2v2L1: 1, o1v2L1: 1, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= dp_bb_oo("j,i") * l2_xbba_Lovv("I,i,a,e") * u1_bb_vo("a,j");

            // sigmal1_xa_Lv += 1.000000 d-_aa(a,b) l2_aaa(i,a,e) u1_aa(b,i) 
            // flops: o1v3L1: 1, o1v2L1: 1, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += dp_aa_vv("a,b") * l2_xaaa_Lovv("I,i,a,e") * u1_aa_vo("b,i");

            // sigmal1_xa_Lv += 1.000000 d-_bb(a,b) l2_bba(i,a,e) u1_bb(b,i) 
            // flops: o1v3L1: 1, o1v2L1: 1, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += dp_bb_vv("a,b") * l2_xbba_Lovv("I,i,a,e") * u1_bb_vo("b,i");

            // sigmal1_xa_Lv += -1.000000 d-_aa(a,e) l2_aaa(i,a,b) u1_aa(b,i) 
            // flops: o1v3L1: 1, o1v2L1: 1, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= dp_aa_vv("a,e") * l2_xaaa_Lovv("I,i,a,b") * u1_aa_vo("b,i");
        }

        if (include_u2_) {

            // sigmal1_xa_Lv += -1.000000 d-_aa(i,a) l2_aaa(j,b,e) u2_aaaa(a,b,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= dp_aa_ov("i,a") * l2_xaaa_Lovv("I,j,b,e") * u2_aaaa_vvoo("a,b,j,i");

            // sigmal1_xa_Lv += 1.000000 d-_bb(i,a) l2_aaa(j,b,e) u2_abab(b,a,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += dp_bb_ov("i,a") * l2_xaaa_Lovv("I,j,b,e") * u2_abab_vvoo("b,a,j,i");

            // sigmal1_xa_Lv += 1.000000 d-_aa(i,a) l2_bba(j,b,e) u2_abab(a,b,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += dp_aa_ov("i,a") * l2_xbba_Lovv("I,j,b,e") * u2_abab_vvoo("a,b,i,j");

            // sigmal1_xa_Lv += -1.000000 d-_bb(i,a) l2_bba(j,b,e) u2_bbbb(a,b,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= dp_bb_ov("i,a") * l2_xbba_Lovv("I,j,b,e") * u2_bbbb_vvoo("a,b,j,i");

            // sigmal1_xa_Lv += -0.500000 d-_aa(i,e) l2_aaa(j,a,b) u2_aaaa(a,b,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * dp_aa_ov("i,e") * l2_xaaa_Lovv("I,j,a,b") * u2_aaaa_vvoo("a,b,j,i");

            // sigmal1_xa_Lv += -0.500000 d-_aa(i,e) l2_bba(j,a,b) u2_abab(b,a,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * dp_aa_ov("i,e") * l2_xbba_Lovv("I,j,a,b") * u2_abab_vvoo("b,a,i,j");
        }

        if (include_u0_) {

            // sigmal1_xa_Lv += -1.000000 d-_aa(j,b) l2_aaa(i,a,e) t2_aaaa(b,a,i,j) u0 
            // flops: o2v3L1: 2, o0v1L1: 2 | mem: o2v3L1: 1, o0v1L1: 2, o0v1: 1, 
            sigmal1_xa_Lv("I,e") -= dp_aa_ov("j,b") * l2_xaaa_Lovv("I,i,a,e") * t2_aaaa_vvoo("b,a,i,j") * u0;

            // sigmal1_xa_Lv += 1.000000 d-_bb(j,b) l2_aaa(i,a,e) t2_abab(a,b,i,j) u0 
            // flops: o2v3L1: 2, o0v1L1: 2 | mem: o2v3L1: 1, o0v1L1: 2, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += dp_bb_ov("j,b") * l2_xaaa_Lovv("I,i,a,e") * t2_abab_vvoo("a,b,i,j") * u0;

            // sigmal1_xa_Lv += 1.000000 d-_aa(j,b) l2_bba(i,a,e) t2_abab(b,a,j,i) u0 
            // flops: o2v3L1: 2, o0v1L1: 2 | mem: o2v3L1: 1, o0v1L1: 2, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += dp_aa_ov("j,b") * l2_xbba_Lovv("I,i,a,e") * t2_abab_vvoo("b,a,j,i") * u0;

            // sigmal1_xa_Lv += -1.000000 d-_bb(j,b) l2_bba(i,a,e) t2_bbbb(b,a,i,j) u0 
            // flops: o2v3L1: 2, o0v1L1: 2 | mem: o2v3L1: 1, o0v1L1: 2, o0v1: 1, 
            sigmal1_xa_Lv("I,e") -= dp_bb_ov("j,b") * l2_xbba_Lovv("I,i,a,e") * t2_bbbb_vvoo("b,a,i,j") * u0;

            // sigmal1_xa_Lv += -0.500000 d-_aa(j,e) l2_aaa(i,b,a) t2_aaaa(b,a,i,j) u0 
            // flops: o2v3L1: 2, o0v1L1: 2 | mem: o2v3L1: 1, o0v1L1: 2, o0v1: 1, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * dp_aa_ov("j,e") * l2_xaaa_Lovv("I,i,b,a") * t2_aaaa_vvoo("b,a,i,j") * u0;

            // sigmal1_xa_Lv += -0.500000 d-_aa(j,e) l2_bba(i,b,a) t2_abab(a,b,j,i) u0 
            // flops: o2v3L1: 2, o0v1L1: 2 | mem: o2v3L1: 1, o0v1L1: 2, o0v1: 1, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * dp_aa_ov("j,e") * l2_xbba_Lovv("I,i,b,a") * t2_abab_vvoo("a,b,j,i") * u0;
        }

        if (include_u0_ && include_u1_) {

            // sigmal1_xa_Lv += -1.000000 d-_aa(i,a) u0 u1_aa(a,i) m1_a(e) 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmal1_xa_Lv("I,e") -= scalar0 * u0 * m1_xa_Lv("I,e");

            // sigmal1_xa_Lv += -1.000000 d-_bb(i,a) u0 u1_bb(a,i) m1_a(e) 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmal1_xa_Lv("I,e") -= scalar1 * u0 * m1_xa_Lv("I,e");

            // sigmal1_xa_Lv += 1.000000 d-_aa(i,e) u0 u1_aa(a,i) m1_a(a) 
            // flops: o0v2L1: 1, o1v2: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o0v2: 1, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += dp_aa_ov("i,e") * u0 * u1_aa_vo("a,i") * m1_xa_Lv("I,a");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmal1_xa_Lv += -1.000000 d-_aa(j,i) u0 u1_aa(a,j) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o2v1: 1, o0v1L1: 1, o2v0: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xa_Lv("I,e") -= dp_aa_oo("j,i") * u0 * u1_aa_vo("a,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -1.000000 d-_bb(j,i) u0 u1_bb(a,j) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o2v1: 1, o0v1L1: 1, o2v0: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xa_Lv("I,e") -= dp_bb_oo("j,i") * u0 * u1_bb_vo("a,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 d-_aa(a,b) u0 u1_aa(b,i) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1, o0v2: 1 | mem: o0v1L1: 2, o0v2: 1, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += dp_aa_vv("a,b") * u0 * u1_aa_vo("b,i") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 d-_bb(a,b) u0 u1_bb(b,i) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1, o0v2: 1 | mem: o0v1L1: 2, o0v2: 1, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += dp_bb_vv("a,b") * u0 * u1_bb_vo("b,i") * m2_xbba_Lovv("I,i,a,e");
        }

        if (include_u0_ && include_u2_) {

            // sigmal1_xa_Lv += -1.000000 d-_aa(j,b) u0 u2_aaaa(b,a,i,j) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 2, 
            sigmal1_xa_Lv("I,e") -= dp_aa_ov("j,b") * u0 * u2_aaaa_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 d-_aa(j,b) u0 u2_abab(b,a,j,i) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 2, 
            sigmal1_xa_Lv("I,e") += dp_aa_ov("j,b") * u0 * u2_abab_vvoo("b,a,j,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 d-_bb(j,b) u0 u2_abab(a,b,i,j) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 2, 
            sigmal1_xa_Lv("I,e") += dp_bb_ov("j,b") * u0 * u2_abab_vvoo("a,b,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -1.000000 d-_bb(j,b) u0 u2_bbbb(b,a,i,j) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 2, 
            sigmal1_xa_Lv("I,e") -= dp_bb_ov("j,b") * u0 * u2_bbbb_vvoo("b,a,i,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -0.500000 d-_aa(j,e) u0 u2_aaaa(b,a,i,j) m2_aaa(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1, o1v1: 1 | mem: o1v3: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * dp_aa_ov("j,e") * u0 * u2_aaaa_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += -0.500000 d-_aa(j,e) u0 u2_abab(a,b,j,i) m2_bba(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1, o1v1: 1 | mem: o1v3: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * dp_aa_ov("j,e") * u0 * u2_abab_vvoo("a,b,j,i") * m2_xbba_Lovv("I,i,b,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmal1_xa_Lv += 1.000000 d-_aa(j,b) u1_aa(b,j) u1_aa(a,i) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += scalar0 * u1_aa_vo("a,i") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 d-_aa(j,b) u1_aa(b,j) u1_bb(a,i) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += scalar0 * u1_bb_vo("a,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 d-_bb(j,b) u1_bb(b,j) u1_aa(a,i) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += scalar1 * u1_aa_vo("a,i") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 d-_bb(j,b) u1_bb(b,j) u1_bb(a,i) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += scalar1 * u1_bb_vo("a,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -2.000000 d-_aa(j,b) u1_aa(b,i) u1_aa(a,j) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o2v1: 2, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xa_Lv("I,e") -= 2.000000 * dp_aa_ov("j,b") * u1_aa_vo("b,i") * u1_aa_vo("a,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -2.000000 d-_bb(j,b) u1_bb(b,i) u1_bb(a,j) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o2v1: 2, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xa_Lv("I,e") -= 2.000000 * dp_bb_ov("j,b") * u1_bb_vo("b,i") * u1_bb_vo("a,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 d-_aa(j,e) u1_aa(b,j) u1_aa(a,i) m2_aaa(i,b,a) 
            // flops: o1v3L1: 1, o1v3: 1, o1v2: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, o0v2: 1, 
            sigmal1_xa_Lv("I,e") += dp_aa_ov("j,e") * u1_aa_vo("b,j") * u1_aa_vo("a,i") * m2_xaaa_Lovv("I,i,b,a");
        }

        {

            // sigmal1_xa_Lv += -0.500000 <j,i||j,i>_aaaa l1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * scalar9 * l1_xa_Lv("I,e");

            // sigmal1_xa_Lv += -0.500000 <j,i||j,i>_abab l1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * scalar10 * l1_xa_Lv("I,e");

            // sigmal1_xa_Lv += -0.500000 <i,j||i,j>_abab l1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * scalar10 * l1_xa_Lv("I,e");

            // sigmal1_xa_Lv += -0.500000 <j,i||j,i>_bbbb l1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * scalar11 * l1_xa_Lv("I,e");

            // sigmal1_xa_Lv += 0.500000 <b,a||e,i>_aaaa l2_aaa(i,b,a) 
            // flops: o1v3L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["aaaa_vvvo"]("b,a,e,i") * l2_xaaa_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += -0.500000 <a,b||e,i>_abab l2_bba(i,b,a) 
            // flops: o1v3L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_vvvo"]("a,b,e,i") * l2_xbba_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += 0.250000 <j,i||a,b>_aaaa l1_a(e) t2_aaaa(a,b,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.250000 * scalar2 * l1_xa_Lv("I,e");

            // sigmal1_xa_Lv += 0.250000 <j,i||a,b>_abab l1_a(e) t2_abab(a,b,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.250000 * scalar3 * l1_xa_Lv("I,e");

            // sigmal1_xa_Lv += 0.250000 <i,j||a,b>_abab l1_a(e) t2_abab(a,b,i,j) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.250000 * scalar3 * l1_xa_Lv("I,e");

            // sigmal1_xa_Lv += 0.250000 <j,i||b,a>_abab l1_a(e) t2_abab(b,a,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.250000 * scalar3 * l1_xa_Lv("I,e");

            // sigmal1_xa_Lv += 0.250000 <i,j||b,a>_abab l1_a(e) t2_abab(b,a,i,j) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.250000 * scalar3 * l1_xa_Lv("I,e");

            // sigmal1_xa_Lv += 0.250000 <j,i||a,b>_bbbb l1_a(e) t2_bbbb(a,b,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.250000 * scalar4 * l1_xa_Lv("I,e");

            // sigmal1_xa_Lv += -0.500000 <j,i||b,e>_aaaa l1_a(a) t2_aaaa(b,a,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["aaaa_oovv"]("j,i,b,e") * l1_xa_Lv("I,a") * t2_aaaa_vvoo("b,a,j,i");

            // sigmal1_xa_Lv += -0.500000 <j,i||e,b>_abab l1_a(a) t2_abab(a,b,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("j,i,e,b") * l1_xa_Lv("I,a") * t2_abab_vvoo("a,b,j,i");

            // sigmal1_xa_Lv += -0.500000 <i,j||e,b>_abab l1_a(a) t2_abab(a,b,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("i,j,e,b") * l1_xa_Lv("I,a") * t2_abab_vvoo("a,b,i,j");

            // sigmal1_xa_Lv += 0.500000 <k,j||b,i>_aaaa l2_aaa(i,a,e) t2_aaaa(b,a,k,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["aaaa_oovo"]("k,j,b,i") * l2_xaaa_Lovv("I,i,a,e") * t2_aaaa_vvoo("b,a,k,j");

            // sigmal1_xa_Lv += 0.500000 <k,j||i,b>_abab l2_aaa(i,a,e) t2_abab(a,b,k,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["abba_oovo"]("k,j,b,i") * l2_xaaa_Lovv("I,i,a,e") * t2_abab_vvoo("a,b,k,j");

            // sigmal1_xa_Lv += 0.500000 <j,k||i,b>_abab l2_aaa(i,a,e) t2_abab(a,b,j,k) 
            // flops: o3v3L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["abba_oovo"]("j,k,b,i") * l2_xaaa_Lovv("I,i,a,e") * t2_abab_vvoo("a,b,j,k");

            // sigmal1_xa_Lv += 0.500000 <k,j||b,i>_abab l2_bba(i,a,e) t2_abab(b,a,k,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["abab_oovo"]("k,j,b,i") * l2_xbba_Lovv("I,i,a,e") * t2_abab_vvoo("b,a,k,j");

            // sigmal1_xa_Lv += 0.500000 <j,k||b,i>_abab l2_bba(i,a,e) t2_abab(b,a,j,k) 
            // flops: o3v3L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["abab_oovo"]("j,k,b,i") * l2_xbba_Lovv("I,i,a,e") * t2_abab_vvoo("b,a,j,k");

            // sigmal1_xa_Lv += 0.500000 <k,j||b,i>_bbbb l2_bba(i,a,e) t2_bbbb(b,a,k,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["bbbb_oovo"]("k,j,b,i") * l2_xbba_Lovv("I,i,a,e") * t2_bbbb_vvoo("b,a,k,j");

            // sigmal1_xa_Lv += 0.250000 <k,j||e,i>_aaaa l2_aaa(i,b,a) t2_aaaa(b,a,k,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += 0.250000 * V_blks_["aaaa_oovo"]("k,j,e,i") * l2_xaaa_Lovv("I,i,b,a") * t2_aaaa_vvoo("b,a,k,j");

            // sigmal1_xa_Lv += -0.250000 <k,j||e,i>_abab l2_bba(i,b,a) t2_abab(a,b,k,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.250000 * V_blks_["abab_oovo"]("k,j,e,i") * l2_xbba_Lovv("I,i,b,a") * t2_abab_vvoo("a,b,k,j");

            // sigmal1_xa_Lv += -0.250000 <j,k||e,i>_abab l2_bba(i,b,a) t2_abab(a,b,j,k) 
            // flops: o3v3L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.250000 * V_blks_["abab_oovo"]("j,k,e,i") * l2_xbba_Lovv("I,i,b,a") * t2_abab_vvoo("a,b,j,k");

            // sigmal1_xa_Lv += 0.500000 <j,a||b,c>_aaaa l2_aaa(i,a,e) t2_aaaa(b,c,i,j) 
            // flops: o2v4L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["aaaa_vovv"]("a,j,b,c") * l2_xaaa_Lovv("I,i,a,e") * t2_aaaa_vvoo("b,c,i,j");

            // sigmal1_xa_Lv += -0.500000 <a,j||b,c>_abab l2_aaa(i,a,e) t2_abab(b,c,i,j) 
            // flops: o2v4L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_vovv"]("a,j,b,c") * l2_xaaa_Lovv("I,i,a,e") * t2_abab_vvoo("b,c,i,j");

            // sigmal1_xa_Lv += -0.500000 <a,j||c,b>_abab l2_aaa(i,a,e) t2_abab(c,b,i,j) 
            // flops: o2v4L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_vovv"]("a,j,c,b") * l2_xaaa_Lovv("I,i,a,e") * t2_abab_vvoo("c,b,i,j");

            // sigmal1_xa_Lv += -0.500000 <j,a||b,c>_abab l2_bba(i,a,e) t2_abab(b,c,j,i) 
            // flops: o2v4L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["baab_vovv"]("a,j,b,c") * l2_xbba_Lovv("I,i,a,e") * t2_abab_vvoo("b,c,j,i");

            // sigmal1_xa_Lv += -0.500000 <j,a||c,b>_abab l2_bba(i,a,e) t2_abab(c,b,j,i) 
            // flops: o2v4L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["baab_vovv"]("a,j,c,b") * l2_xbba_Lovv("I,i,a,e") * t2_abab_vvoo("c,b,j,i");

            // sigmal1_xa_Lv += 0.500000 <j,a||b,c>_bbbb l2_bba(i,a,e) t2_bbbb(b,c,i,j) 
            // flops: o2v4L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["bbbb_vovv"]("a,j,b,c") * l2_xbba_Lovv("I,i,a,e") * t2_bbbb_vvoo("b,c,i,j");

            // sigmal1_xa_Lv += -1.000000 <j,b||c,e>_aaaa l2_aaa(i,b,a) t2_aaaa(c,a,i,j) 
            // flops: o2v4L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += V_blks_["aaaa_vovv"]("b,j,c,e") * l2_xaaa_Lovv("I,i,b,a") * t2_aaaa_vvoo("c,a,i,j");

            // sigmal1_xa_Lv += 1.000000 <b,j||e,c>_abab l2_aaa(i,b,a) t2_abab(a,c,i,j) 
            // flops: o2v4L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += V_blks_["abab_vovv"]("b,j,e,c") * l2_xaaa_Lovv("I,i,b,a") * t2_abab_vvoo("a,c,i,j");

            // sigmal1_xa_Lv += 1.000000 <j,b||e,c>_abab l2_bba(i,b,a) t2_abab(a,c,j,i) 
            // flops: o2v4L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= V_blks_["baab_vovv"]("b,j,e,c") * l2_xbba_Lovv("I,i,b,a") * t2_abab_vvoo("a,c,j,i");
        }

        if (include_u1_) {

            // sigmal1_xa_Lv += 1.000000 <i,a||b,e>_aaaa u1_aa(b,i) m1_a(a) 
            // flops: o1v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmal1_xa_Lv("I,e") -= V_blks_["aaaa_vovv"]("a,i,b,e") * u1_aa_vo("b,i") * m1_xa_Lv("I,a");

            // sigmal1_xa_Lv += 1.000000 <a,i||e,b>_abab u1_bb(b,i) m1_a(a) 
            // flops: o1v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmal1_xa_Lv("I,e") += V_blks_["abab_vovv"]("a,i,e,b") * u1_bb_vo("b,i") * m1_xa_Lv("I,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmal1_xa_Lv += -1.000000 <j,a||b,i>_aaaa u1_aa(b,j) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += V_blks_["aaaa_vovo"]("a,j,b,i") * u1_aa_vo("b,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -1.000000 <j,a||b,i>_abab u1_aa(b,j) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += V_blks_["baab_vovo"]("a,j,b,i") * u1_aa_vo("b,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -1.000000 <a,j||i,b>_abab u1_bb(b,j) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += V_blks_["abba_vovo"]("a,j,b,i") * u1_bb_vo("b,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -1.000000 <j,a||b,i>_bbbb u1_bb(b,j) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += V_blks_["bbbb_vovo"]("a,j,b,i") * u1_bb_vo("b,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 <j,b||e,i>_aaaa u1_aa(a,j) m2_aaa(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= V_blks_["aaaa_vovo"]("b,j,e,i") * u1_aa_vo("a,j") * m2_xaaa_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += 1.000000 <j,b||e,i>_abab u1_aa(a,j) m2_bba(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= V_blks_["baab_vovo"]("b,j,e,i") * u1_aa_vo("a,j") * m2_xbba_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += -0.500000 <b,a||c,e>_aaaa u1_aa(c,i) m2_aaa(i,b,a) 
            // flops: o1v3L1: 1, o1v4: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["aaaa_vvvv"]("b,a,c,e") * u1_aa_vo("c,i") * m2_xaaa_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += -0.500000 <a,b||e,c>_abab u1_bb(c,i) m2_bba(i,b,a) 
            // flops: o1v3L1: 1, o1v4: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_vvvv"]("a,b,e,c") * u1_bb_vo("c,i") * m2_xbba_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += 0.250000 <j,i||a,b>_aaaa u2_aaaa(a,b,j,i) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.250000 * scalar14 * m1_xa_Lv("I,e");

            // sigmal1_xa_Lv += 0.250000 <j,i||a,b>_abab u2_abab(a,b,j,i) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.250000 * scalar15 * m1_xa_Lv("I,e");

            // sigmal1_xa_Lv += 0.250000 <i,j||a,b>_abab u2_abab(a,b,i,j) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.250000 * scalar15 * m1_xa_Lv("I,e");

            // sigmal1_xa_Lv += 0.250000 <j,i||b,a>_abab u2_abab(b,a,j,i) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.250000 * scalar15 * m1_xa_Lv("I,e");

            // sigmal1_xa_Lv += 0.250000 <i,j||b,a>_abab u2_abab(b,a,i,j) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.250000 * scalar15 * m1_xa_Lv("I,e");

            // sigmal1_xa_Lv += 0.250000 <j,i||a,b>_bbbb u2_bbbb(a,b,j,i) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.250000 * scalar16 * m1_xa_Lv("I,e");

            // sigmal1_xa_Lv += -0.500000 <j,i||b,e>_aaaa u2_aaaa(b,a,j,i) m1_a(a) 
            // flops: o2v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["aaaa_oovv"]("j,i,b,e") * u2_aaaa_vvoo("b,a,j,i") * m1_xa_Lv("I,a");

            // sigmal1_xa_Lv += -0.500000 <j,i||e,b>_abab u2_abab(a,b,j,i) m1_a(a) 
            // flops: o2v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("j,i,e,b") * u2_abab_vvoo("a,b,j,i") * m1_xa_Lv("I,a");

            // sigmal1_xa_Lv += -0.500000 <i,j||e,b>_abab u2_abab(a,b,i,j) m1_a(a) 
            // flops: o2v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("i,j,e,b") * u2_abab_vvoo("a,b,i,j") * m1_xa_Lv("I,a");
        }

        if (include_u2_) {

            // sigmal1_xa_Lv += 0.500000 <k,j||b,i>_aaaa u2_aaaa(b,a,k,j) m2_aaa(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["aaaa_oovo"]("k,j,b,i") * u2_aaaa_vvoo("b,a,k,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.500000 <k,j||b,i>_abab u2_abab(b,a,k,j) m2_bba(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["abab_oovo"]("k,j,b,i") * u2_abab_vvoo("b,a,k,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.500000 <j,k||b,i>_abab u2_abab(b,a,j,k) m2_bba(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["abab_oovo"]("j,k,b,i") * u2_abab_vvoo("b,a,j,k") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.500000 <k,j||i,b>_abab u2_abab(a,b,k,j) m2_aaa(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["abba_oovo"]("k,j,b,i") * u2_abab_vvoo("a,b,k,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.500000 <j,k||i,b>_abab u2_abab(a,b,j,k) m2_aaa(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["abba_oovo"]("j,k,b,i") * u2_abab_vvoo("a,b,j,k") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.500000 <k,j||b,i>_bbbb u2_bbbb(b,a,k,j) m2_bba(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["bbbb_oovo"]("k,j,b,i") * u2_bbbb_vvoo("b,a,k,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.250000 <k,j||e,i>_aaaa u2_aaaa(b,a,k,j) m2_aaa(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += 0.250000 * V_blks_["aaaa_oovo"]("k,j,e,i") * u2_aaaa_vvoo("b,a,k,j") * m2_xaaa_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += -0.250000 <k,j||e,i>_abab u2_abab(a,b,k,j) m2_bba(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.250000 * V_blks_["abab_oovo"]("k,j,e,i") * u2_abab_vvoo("a,b,k,j") * m2_xbba_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += -0.250000 <j,k||e,i>_abab u2_abab(a,b,j,k) m2_bba(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.250000 * V_blks_["abab_oovo"]("j,k,e,i") * u2_abab_vvoo("a,b,j,k") * m2_xbba_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += 0.500000 <j,a||b,c>_aaaa u2_aaaa(b,c,i,j) m2_aaa(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["aaaa_vovv"]("a,j,b,c") * u2_aaaa_vvoo("b,c,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -0.500000 <a,j||b,c>_abab u2_abab(b,c,i,j) m2_aaa(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_vovv"]("a,j,b,c") * u2_abab_vvoo("b,c,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -0.500000 <j,a||b,c>_abab u2_abab(b,c,j,i) m2_bba(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["baab_vovv"]("a,j,b,c") * u2_abab_vvoo("b,c,j,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -0.500000 <a,j||c,b>_abab u2_abab(c,b,i,j) m2_aaa(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_vovv"]("a,j,c,b") * u2_abab_vvoo("c,b,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -0.500000 <j,a||c,b>_abab u2_abab(c,b,j,i) m2_bba(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["baab_vovv"]("a,j,c,b") * u2_abab_vvoo("c,b,j,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.500000 <j,a||b,c>_bbbb u2_bbbb(b,c,i,j) m2_bba(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") -= 0.500000 * V_blks_["bbbb_vovv"]("a,j,b,c") * u2_bbbb_vvoo("b,c,i,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -1.000000 <j,b||c,e>_aaaa u2_aaaa(c,a,i,j) m2_aaa(i,b,a) 
            // flops: o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += V_blks_["aaaa_vovv"]("b,j,c,e") * u2_aaaa_vvoo("c,a,i,j") * m2_xaaa_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += 1.000000 <b,j||e,c>_abab u2_abab(a,c,i,j) m2_aaa(i,b,a) 
            // flops: o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += V_blks_["abab_vovv"]("b,j,e,c") * u2_abab_vvoo("a,c,i,j") * m2_xaaa_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += 1.000000 <j,b||e,c>_abab u2_abab(a,c,j,i) m2_bba(i,b,a) 
            // flops: o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= V_blks_["baab_vovv"]("b,j,e,c") * u2_abab_vvoo("a,c,j,i") * m2_xbba_Lovv("I,i,b,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmal1_xa_Lv += 0.500000 <k,j||b,c>_aaaa t2_aaaa(b,c,i,j) u1_aa(a,k) m2_aaa(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["aaaa_oovv"]("k,j,b,c") * t2_aaaa_vvoo("b,c,i,j") * u1_aa_vo("a,k") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.500000 <k,j||b,c>_abab t2_abab(b,c,i,j) u1_aa(a,k) m2_aaa(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("k,j,b,c") * t2_abab_vvoo("b,c,i,j") * u1_aa_vo("a,k") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.500000 <j,k||b,c>_abab t2_abab(b,c,j,i) u1_bb(a,k) m2_bba(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("j,k,b,c") * t2_abab_vvoo("b,c,j,i") * u1_bb_vo("a,k") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.500000 <k,j||c,b>_abab t2_abab(c,b,i,j) u1_aa(a,k) m2_aaa(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("k,j,c,b") * t2_abab_vvoo("c,b,i,j") * u1_aa_vo("a,k") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.500000 <j,k||c,b>_abab t2_abab(c,b,j,i) u1_bb(a,k) m2_bba(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("j,k,c,b") * t2_abab_vvoo("c,b,j,i") * u1_bb_vo("a,k") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.500000 <k,j||b,c>_bbbb t2_bbbb(b,c,i,j) u1_bb(a,k) m2_bba(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["bbbb_oovv"]("k,j,b,c") * t2_bbbb_vvoo("b,c,i,j") * u1_bb_vo("a,k") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.500000 <k,j||b,c>_aaaa t2_aaaa(b,a,k,j) u1_aa(c,i) m2_aaa(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["aaaa_oovv"]("k,j,b,c") * t2_aaaa_vvoo("b,a,k,j") * u1_aa_vo("c,i") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.500000 <k,j||b,c>_abab t2_abab(b,a,k,j) u1_bb(c,i) m2_bba(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("k,j,b,c") * t2_abab_vvoo("b,a,k,j") * u1_bb_vo("c,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.500000 <j,k||b,c>_abab t2_abab(b,a,j,k) u1_bb(c,i) m2_bba(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("j,k,b,c") * t2_abab_vvoo("b,a,j,k") * u1_bb_vo("c,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.500000 <k,j||c,b>_abab t2_abab(a,b,k,j) u1_aa(c,i) m2_aaa(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("k,j,c,b") * t2_abab_vvoo("a,b,k,j") * u1_aa_vo("c,i") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.500000 <j,k||c,b>_abab t2_abab(a,b,j,k) u1_aa(c,i) m2_aaa(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("j,k,c,b") * t2_abab_vvoo("a,b,j,k") * u1_aa_vo("c,i") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 0.500000 <k,j||b,c>_bbbb t2_bbbb(b,a,k,j) u1_bb(c,i) m2_bba(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["bbbb_oovv"]("k,j,b,c") * t2_bbbb_vvoo("b,a,k,j") * u1_bb_vo("c,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -1.000000 <k,j||b,c>_aaaa t2_aaaa(b,a,i,j) u1_aa(c,k) m2_aaa(i,a,e) 
            // flops: o3v3: 1, o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o2v2: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") -= V_blks_["aaaa_oovv"]("k,j,b,c") * t2_aaaa_vvoo("b,a,i,j") * u1_aa_vo("c,k") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 <j,k||b,c>_abab t2_aaaa(b,a,i,j) u1_bb(c,k) m2_aaa(i,a,e) 
            // flops: o3v3: 1, o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o2v2: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += V_blks_["abab_oovv"]("j,k,b,c") * t2_aaaa_vvoo("b,a,i,j") * u1_bb_vo("c,k") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 <k,j||b,c>_aaaa t2_abab(b,a,j,i) u1_aa(c,k) m2_bba(i,a,e) 
            // flops: o3v3: 1, o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o2v2: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += V_blks_["aaaa_oovv"]("k,j,b,c") * t2_abab_vvoo("b,a,j,i") * u1_aa_vo("c,k") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -1.000000 <j,k||b,c>_abab t2_abab(b,a,j,i) u1_bb(c,k) m2_bba(i,a,e) 
            // flops: o3v3: 1, o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o2v2: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") -= V_blks_["abab_oovv"]("j,k,b,c") * t2_abab_vvoo("b,a,j,i") * u1_bb_vo("c,k") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -1.000000 <k,j||c,b>_abab t2_abab(a,b,i,j) u1_aa(c,k) m2_aaa(i,a,e) 
            // flops: o3v3: 1, o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o2v2: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") -= V_blks_["abab_oovv"]("k,j,c,b") * t2_abab_vvoo("a,b,i,j") * u1_aa_vo("c,k") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 <k,j||b,c>_bbbb t2_abab(a,b,i,j) u1_bb(c,k) m2_aaa(i,a,e) 
            // flops: o3v3: 1, o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o2v2: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += V_blks_["bbbb_oovv"]("k,j,b,c") * t2_abab_vvoo("a,b,i,j") * u1_bb_vo("c,k") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 <k,j||c,b>_abab t2_bbbb(b,a,i,j) u1_aa(c,k) m2_bba(i,a,e) 
            // flops: o3v3: 1, o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o2v2: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") += V_blks_["abab_oovv"]("k,j,c,b") * t2_bbbb_vvoo("b,a,i,j") * u1_aa_vo("c,k") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += -1.000000 <k,j||b,c>_bbbb t2_bbbb(b,a,i,j) u1_bb(c,k) m2_bba(i,a,e) 
            // flops: o3v3: 1, o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o2v2: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xa_Lv("I,e") -= V_blks_["bbbb_oovv"]("k,j,b,c") * t2_bbbb_vvoo("b,a,i,j") * u1_bb_vo("c,k") * m2_xbba_Lovv("I,i,a,e");

            // sigmal1_xa_Lv += 1.000000 <k,j||c,e>_aaaa t2_aaaa(c,b,i,j) u1_aa(a,k) m2_aaa(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o2v2: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += V_blks_["aaaa_oovv"]("k,j,c,e") * t2_aaaa_vvoo("c,b,i,j") * u1_aa_vo("a,k") * m2_xaaa_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += -1.000000 <k,j||c,e>_aaaa t2_abab(c,b,j,i) u1_aa(a,k) m2_bba(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o2v2: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= V_blks_["aaaa_oovv"]("k,j,c,e") * t2_abab_vvoo("c,b,j,i") * u1_aa_vo("a,k") * m2_xbba_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += 1.000000 <k,j||e,c>_abab t2_abab(b,c,i,j) u1_aa(a,k) m2_aaa(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o2v2: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += V_blks_["abab_oovv"]("k,j,e,c") * t2_abab_vvoo("b,c,i,j") * u1_aa_vo("a,k") * m2_xaaa_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += -1.000000 <k,j||e,c>_abab t2_bbbb(c,b,i,j) u1_aa(a,k) m2_bba(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o2v2: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= V_blks_["abab_oovv"]("k,j,e,c") * t2_bbbb_vvoo("c,b,i,j") * u1_aa_vo("a,k") * m2_xbba_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += -0.250000 <k,j||c,e>_aaaa t2_aaaa(b,a,k,j) u1_aa(c,i) m2_aaa(i,b,a) 
            // flops: o2v4: 1, o1v3L1: 1, o1v4: 1, o0v1L1: 1 | mem: o0v4: 1, o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.250000 * V_blks_["aaaa_oovv"]("k,j,c,e") * t2_aaaa_vvoo("b,a,k,j") * u1_aa_vo("c,i") * m2_xaaa_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += -0.250000 <k,j||e,c>_abab t2_abab(a,b,k,j) u1_bb(c,i) m2_bba(i,b,a) 
            // flops: o2v4: 1, o1v3L1: 1, o1v4: 1, o0v1L1: 1 | mem: o0v4: 1, o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.250000 * V_blks_["abab_oovv"]("k,j,e,c") * t2_abab_vvoo("a,b,k,j") * u1_bb_vo("c,i") * m2_xbba_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += -0.250000 <j,k||e,c>_abab t2_abab(a,b,j,k) u1_bb(c,i) m2_bba(i,b,a) 
            // flops: o2v4: 1, o1v3L1: 1, o1v4: 1, o0v1L1: 1 | mem: o0v4: 1, o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") -= 0.250000 * V_blks_["abab_oovv"]("j,k,e,c") * t2_abab_vvoo("a,b,j,k") * u1_bb_vo("c,i") * m2_xbba_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += 0.500000 <k,j||c,e>_aaaa t2_aaaa(b,a,i,j) u1_aa(c,k) m2_aaa(i,b,a) 
            // flops: o3v4: 1, o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o2v4: 1, o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["aaaa_oovv"]("k,j,c,e") * t2_aaaa_vvoo("b,a,i,j") * u1_aa_vo("c,k") * m2_xaaa_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += 0.500000 <j,k||e,c>_abab t2_aaaa(b,a,i,j) u1_bb(c,k) m2_aaa(i,b,a) 
            // flops: o3v4: 1, o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o2v4: 1, o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("j,k,e,c") * t2_aaaa_vvoo("b,a,i,j") * u1_bb_vo("c,k") * m2_xaaa_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += 0.500000 <k,j||c,e>_aaaa t2_abab(a,b,j,i) u1_aa(c,k) m2_bba(i,b,a) 
            // flops: o3v4: 1, o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o2v4: 1, o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["aaaa_oovv"]("k,j,c,e") * t2_abab_vvoo("a,b,j,i") * u1_aa_vo("c,k") * m2_xbba_Lovv("I,i,b,a");

            // sigmal1_xa_Lv += 0.500000 <j,k||e,c>_abab t2_abab(a,b,j,i) u1_bb(c,k) m2_bba(i,b,a) 
            // flops: o3v4: 1, o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o2v4: 1, o1v3: 1, o0v1L1: 2, 
            sigmal1_xa_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("j,k,e,c") * t2_abab_vvoo("a,b,j,i") * u1_bb_vo("c,k") * m2_xbba_Lovv("I,i,b,a");
        }




/// ****** sigmal1_b(e) ****** ///



        if (include_u0_ && include_u1_) {

            // sigmal1_xb_Lv += 1.000000 u0 m1_b(e) w0 
            // flops: o0v1L1: 2, o0v1: 1 | mem: o0v1L1: 1, o0v1: 2, 
            sigmal1_xb_Lv("I,e") += u0 * m1_xb_Lv("I,e") * w0;
        }

        if (include_u1_ && include_u2_) {

            // sigmal1_xb_Lv += -1.000000 u1_aa(a,i) m2_aab(i,a,e) w0 
            // flops: o1v2L1: 1, o0v1L1: 2 | mem: o0v1L1: 2, o0v1: 1, 
            sigmal1_xb_Lv("I,e") -= u1_aa_vo("a,i") * m2_xaab_Lovv("I,i,a,e") * w0;

            // sigmal1_xb_Lv += -1.000000 u1_bb(a,i) m2_bbb(i,a,e) w0 
            // flops: o1v2L1: 1, o0v1L1: 2 | mem: o0v1L1: 2, o0v1: 1, 
            sigmal1_xb_Lv("I,e") -= u1_bb_vo("a,i") * m2_xbbb_Lovv("I,i,a,e") * w0;
        }

        {

            // sigmal1_xb_Lv += 1.000000 f_aa(i,i) l1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += scalar5 * l1_xb_Lv("I,e");

            // sigmal1_xb_Lv += 1.000000 f_bb(i,i) l1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += scalar6 * l1_xb_Lv("I,e");

            // sigmal1_xb_Lv += 1.000000 f_bb(a,e) l1_b(a) 
            // flops: o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += F_blks_["bb_vv"]("a,e") * l1_xb_Lv("I,a");

            // sigmal1_xb_Lv += -1.000000 f_aa(a,i) l2_aab(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= F_blks_["aa_vo"]("a,i") * l2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -1.000000 f_bb(a,i) l2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= F_blks_["bb_vo"]("a,i") * l2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 f_aa(j,b) l2_aab(i,a,e) t2_aaaa(b,a,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += F_blks_["aa_ov"]("j,b") * l2_xaab_Lovv("I,i,a,e") * t2_aaaa_vvoo("b,a,i,j");

            // sigmal1_xb_Lv += -1.000000 f_bb(j,b) l2_aab(i,a,e) t2_abab(a,b,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= F_blks_["bb_ov"]("j,b") * l2_xaab_Lovv("I,i,a,e") * t2_abab_vvoo("a,b,i,j");

            // sigmal1_xb_Lv += -1.000000 f_aa(j,b) l2_bbb(i,a,e) t2_abab(b,a,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= F_blks_["aa_ov"]("j,b") * l2_xbbb_Lovv("I,i,a,e") * t2_abab_vvoo("b,a,j,i");

            // sigmal1_xb_Lv += 1.000000 f_bb(j,b) l2_bbb(i,a,e) t2_bbbb(b,a,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += F_blks_["bb_ov"]("j,b") * l2_xbbb_Lovv("I,i,a,e") * t2_bbbb_vvoo("b,a,i,j");

            // sigmal1_xb_Lv += 0.500000 f_bb(j,e) l2_aab(i,b,a) t2_abab(b,a,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.500000 * F_blks_["bb_ov"]("j,e") * l2_xaab_Lovv("I,i,b,a") * t2_abab_vvoo("b,a,i,j");

            // sigmal1_xb_Lv += 0.500000 f_bb(j,e) l2_bbb(i,b,a) t2_bbbb(b,a,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.500000 * F_blks_["bb_ov"]("j,e") * l2_xbbb_Lovv("I,i,b,a") * t2_bbbb_vvoo("b,a,i,j");
        }

        if (include_u1_) {

            // sigmal1_xb_Lv += 1.000000 f_aa(i,a) u1_aa(a,i) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += scalar12 * m1_xb_Lv("I,e");

            // sigmal1_xb_Lv += 1.000000 f_bb(i,a) u1_bb(a,i) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += scalar13 * m1_xb_Lv("I,e");

            // sigmal1_xb_Lv += -1.000000 f_bb(i,e) u1_bb(a,i) m1_b(a) 
            // flops: o0v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmal1_xb_Lv("I,e") -= F_blks_["bb_ov"]("i,e") * u1_bb_vo("a,i") * m1_xb_Lv("I,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmal1_xb_Lv += 1.000000 f_aa(j,i) u1_aa(a,j) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += F_blks_["aa_oo"]("j,i") * u1_aa_vo("a,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 f_bb(j,i) u1_bb(a,j) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += F_blks_["bb_oo"]("j,i") * u1_bb_vo("a,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -1.000000 f_aa(a,b) u1_aa(b,i) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") -= F_blks_["aa_vv"]("a,b") * u1_aa_vo("b,i") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -1.000000 f_bb(a,b) u1_bb(b,i) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") -= F_blks_["bb_vv"]("a,b") * u1_bb_vo("b,i") * m2_xbbb_Lovv("I,i,a,e");
        }

        if (include_u2_) {

            // sigmal1_xb_Lv += 1.000000 f_aa(j,b) u2_aaaa(b,a,i,j) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += F_blks_["aa_ov"]("j,b") * u2_aaaa_vvoo("b,a,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -1.000000 f_aa(j,b) u2_abab(b,a,j,i) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") -= F_blks_["aa_ov"]("j,b") * u2_abab_vvoo("b,a,j,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -1.000000 f_bb(j,b) u2_abab(a,b,i,j) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") -= F_blks_["bb_ov"]("j,b") * u2_abab_vvoo("a,b,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 f_bb(j,b) u2_bbbb(b,a,i,j) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += F_blks_["bb_ov"]("j,b") * u2_bbbb_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 0.500000 f_bb(j,e) u2_abab(b,a,i,j) m2_aab(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.500000 * F_blks_["bb_ov"]("j,e") * u2_abab_vvoo("b,a,i,j") * m2_xaab_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += 0.500000 f_bb(j,e) u2_bbbb(b,a,i,j) m2_bbb(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.500000 * F_blks_["bb_ov"]("j,e") * u2_bbbb_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,i,b,a");
        }

        if (include_u1_) {

            // sigmal1_xb_Lv += -1.000000 d+_aa(i,i) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") -= scalar17 * m1_xb_Lv("I,e");

            // sigmal1_xb_Lv += -1.000000 d+_bb(i,i) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") -= scalar18 * m1_xb_Lv("I,e");

            // sigmal1_xb_Lv += -1.000000 d+_bb(a,e) m1_b(a) 
            // flops: o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= dp_bb_vv("a,e") * m1_xb_Lv("I,a");
        }

        if (include_u2_) {

            // sigmal1_xb_Lv += 1.000000 d+_aa(a,i) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += dp_aa_vo("a,i") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 d+_bb(a,i) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += dp_bb_vo("a,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -1.000000 d+_aa(j,b) t2_aaaa(b,a,i,j) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") -= dp_aa_ov("j,b") * t2_aaaa_vvoo("b,a,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 d+_aa(j,b) t2_abab(b,a,j,i) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += dp_aa_ov("j,b") * t2_abab_vvoo("b,a,j,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 d+_bb(j,b) t2_abab(a,b,i,j) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += dp_bb_ov("j,b") * t2_abab_vvoo("a,b,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -1.000000 d+_bb(j,b) t2_bbbb(b,a,i,j) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") -= dp_bb_ov("j,b") * t2_bbbb_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -0.500000 d+_bb(j,e) t2_abab(b,a,i,j) m2_aab(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * dp_bb_ov("j,e") * t2_abab_vvoo("b,a,i,j") * m2_xaab_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += -0.500000 d+_bb(j,e) t2_bbbb(b,a,i,j) m2_bbb(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * dp_bb_ov("j,e") * t2_bbbb_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,i,b,a");
        }

        if (include_u0_) {

            // sigmal1_xb_Lv += -1.000000 d-_aa(i,i) l1_b(e) u0 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmal1_xb_Lv("I,e") -= scalar7 * u0 * l1_xb_Lv("I,e");

            // sigmal1_xb_Lv += -1.000000 d-_bb(i,i) l1_b(e) u0 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmal1_xb_Lv("I,e") -= scalar8 * u0 * l1_xb_Lv("I,e");

            // sigmal1_xb_Lv += -1.000000 d-_bb(a,e) l1_b(a) u0 
            // flops: o0v2L1: 1, o0v1L1: 2 | mem: o0v1L1: 2, o0v1: 1, 
            sigmal1_xb_Lv("I,e") -= dp_bb_vv("a,e") * l1_xb_Lv("I,a") * u0;

            // sigmal1_xb_Lv += 1.000000 d-_aa(a,i) l2_aab(i,a,e) u0 
            // flops: o1v2L1: 1, o0v1L1: 2 | mem: o0v1L1: 2, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += dp_aa_vo("a,i") * l2_xaab_Lovv("I,i,a,e") * u0;

            // sigmal1_xb_Lv += 1.000000 d-_bb(a,i) l2_bbb(i,a,e) u0 
            // flops: o1v2L1: 1, o0v1L1: 2 | mem: o0v1L1: 2, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += dp_bb_vo("a,i") * l2_xbbb_Lovv("I,i,a,e") * u0;
        }

        if (include_u1_) {

            // sigmal1_xb_Lv += -1.000000 d-_aa(i,a) l1_b(e) u1_aa(a,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") -= scalar0 * l1_xb_Lv("I,e");

            // sigmal1_xb_Lv += -1.000000 d-_bb(i,a) l1_b(e) u1_bb(a,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") -= scalar1 * l1_xb_Lv("I,e");

            // sigmal1_xb_Lv += 1.000000 d-_bb(i,e) l1_b(a) u1_bb(a,i) 
            // flops: o1v2L1: 2, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += dp_bb_ov("i,e") * l1_xb_Lv("I,a") * u1_bb_vo("a,i");

            // sigmal1_xb_Lv += 1.000000 d-_aa(i,i) l2_aab(j,a,e) u1_aa(a,j) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += scalar7 * u1_aa_vo("a,j") * l2_xaab_Lovv("I,j,a,e");

            // sigmal1_xb_Lv += 1.000000 d-_bb(i,i) l2_aab(j,a,e) u1_aa(a,j) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += scalar8 * u1_aa_vo("a,j") * l2_xaab_Lovv("I,j,a,e");

            // sigmal1_xb_Lv += 1.000000 d-_aa(i,i) l2_bbb(j,a,e) u1_bb(a,j) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += scalar7 * u1_bb_vo("a,j") * l2_xbbb_Lovv("I,j,a,e");

            // sigmal1_xb_Lv += 1.000000 d-_bb(i,i) l2_bbb(j,a,e) u1_bb(a,j) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += scalar8 * u1_bb_vo("a,j") * l2_xbbb_Lovv("I,j,a,e");

            // sigmal1_xb_Lv += -1.000000 d-_aa(j,i) l2_aab(i,a,e) u1_aa(a,j) 
            // flops: o2v2L1: 1, o1v2L1: 1, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= dp_aa_oo("j,i") * l2_xaab_Lovv("I,i,a,e") * u1_aa_vo("a,j");

            // sigmal1_xb_Lv += -1.000000 d-_bb(j,i) l2_bbb(i,a,e) u1_bb(a,j) 
            // flops: o2v2L1: 1, o1v2L1: 1, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= dp_bb_oo("j,i") * l2_xbbb_Lovv("I,i,a,e") * u1_bb_vo("a,j");

            // sigmal1_xb_Lv += 1.000000 d-_aa(a,b) l2_aab(i,a,e) u1_aa(b,i) 
            // flops: o1v3L1: 1, o1v2L1: 1, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += dp_aa_vv("a,b") * l2_xaab_Lovv("I,i,a,e") * u1_aa_vo("b,i");

            // sigmal1_xb_Lv += 1.000000 d-_bb(a,b) l2_bbb(i,a,e) u1_bb(b,i) 
            // flops: o1v3L1: 1, o1v2L1: 1, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += dp_bb_vv("a,b") * l2_xbbb_Lovv("I,i,a,e") * u1_bb_vo("b,i");

            // sigmal1_xb_Lv += -1.000000 d-_bb(a,e) l2_bbb(i,a,b) u1_bb(b,i) 
            // flops: o1v3L1: 1, o1v2L1: 1, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= dp_bb_vv("a,e") * l2_xbbb_Lovv("I,i,a,b") * u1_bb_vo("b,i");
        }

        if (include_u2_) {

            // sigmal1_xb_Lv += -1.000000 d-_aa(i,a) l2_aab(j,b,e) u2_aaaa(a,b,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= dp_aa_ov("i,a") * l2_xaab_Lovv("I,j,b,e") * u2_aaaa_vvoo("a,b,j,i");

            // sigmal1_xb_Lv += 1.000000 d-_bb(i,a) l2_aab(j,b,e) u2_abab(b,a,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += dp_bb_ov("i,a") * l2_xaab_Lovv("I,j,b,e") * u2_abab_vvoo("b,a,j,i");

            // sigmal1_xb_Lv += 1.000000 d-_aa(i,a) l2_bbb(j,b,e) u2_abab(a,b,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += dp_aa_ov("i,a") * l2_xbbb_Lovv("I,j,b,e") * u2_abab_vvoo("a,b,i,j");

            // sigmal1_xb_Lv += -1.000000 d-_bb(i,a) l2_bbb(j,b,e) u2_bbbb(a,b,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= dp_bb_ov("i,a") * l2_xbbb_Lovv("I,j,b,e") * u2_bbbb_vvoo("a,b,j,i");

            // sigmal1_xb_Lv += -0.500000 d-_bb(i,e) l2_aab(j,a,b) u2_abab(a,b,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * dp_bb_ov("i,e") * l2_xaab_Lovv("I,j,a,b") * u2_abab_vvoo("a,b,j,i");

            // sigmal1_xb_Lv += -0.500000 d-_bb(i,e) l2_bbb(j,a,b) u2_bbbb(a,b,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * dp_bb_ov("i,e") * l2_xbbb_Lovv("I,j,a,b") * u2_bbbb_vvoo("a,b,j,i");
        }

        if (include_u0_) {

            // sigmal1_xb_Lv += -1.000000 d-_aa(j,b) l2_aab(i,a,e) t2_aaaa(b,a,i,j) u0 
            // flops: o2v3L1: 2, o0v1L1: 2 | mem: o2v3L1: 1, o0v1L1: 2, o0v1: 1, 
            sigmal1_xb_Lv("I,e") -= dp_aa_ov("j,b") * l2_xaab_Lovv("I,i,a,e") * t2_aaaa_vvoo("b,a,i,j") * u0;

            // sigmal1_xb_Lv += 1.000000 d-_bb(j,b) l2_aab(i,a,e) t2_abab(a,b,i,j) u0 
            // flops: o2v3L1: 2, o0v1L1: 2 | mem: o2v3L1: 1, o0v1L1: 2, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += dp_bb_ov("j,b") * l2_xaab_Lovv("I,i,a,e") * t2_abab_vvoo("a,b,i,j") * u0;

            // sigmal1_xb_Lv += 1.000000 d-_aa(j,b) l2_bbb(i,a,e) t2_abab(b,a,j,i) u0 
            // flops: o2v3L1: 2, o0v1L1: 2 | mem: o2v3L1: 1, o0v1L1: 2, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += dp_aa_ov("j,b") * l2_xbbb_Lovv("I,i,a,e") * t2_abab_vvoo("b,a,j,i") * u0;

            // sigmal1_xb_Lv += -1.000000 d-_bb(j,b) l2_bbb(i,a,e) t2_bbbb(b,a,i,j) u0 
            // flops: o2v3L1: 2, o0v1L1: 2 | mem: o2v3L1: 1, o0v1L1: 2, o0v1: 1, 
            sigmal1_xb_Lv("I,e") -= dp_bb_ov("j,b") * l2_xbbb_Lovv("I,i,a,e") * t2_bbbb_vvoo("b,a,i,j") * u0;

            // sigmal1_xb_Lv += -0.500000 d-_bb(j,e) l2_aab(i,b,a) t2_abab(b,a,i,j) u0 
            // flops: o2v3L1: 2, o0v1L1: 2 | mem: o2v3L1: 1, o0v1L1: 2, o0v1: 1, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * dp_bb_ov("j,e") * l2_xaab_Lovv("I,i,b,a") * t2_abab_vvoo("b,a,i,j") * u0;

            // sigmal1_xb_Lv += -0.500000 d-_bb(j,e) l2_bbb(i,b,a) t2_bbbb(b,a,i,j) u0 
            // flops: o2v3L1: 2, o0v1L1: 2 | mem: o2v3L1: 1, o0v1L1: 2, o0v1: 1, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * dp_bb_ov("j,e") * l2_xbbb_Lovv("I,i,b,a") * t2_bbbb_vvoo("b,a,i,j") * u0;
        }

        if (include_u0_ && include_u1_) {

            // sigmal1_xb_Lv += -1.000000 d-_aa(i,a) u0 u1_aa(a,i) m1_b(e) 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmal1_xb_Lv("I,e") -= scalar0 * u0 * m1_xb_Lv("I,e");

            // sigmal1_xb_Lv += -1.000000 d-_bb(i,a) u0 u1_bb(a,i) m1_b(e) 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmal1_xb_Lv("I,e") -= scalar1 * u0 * m1_xb_Lv("I,e");

            // sigmal1_xb_Lv += 1.000000 d-_bb(i,e) u0 u1_bb(a,i) m1_b(a) 
            // flops: o0v2L1: 1, o1v2: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o0v2: 1, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += dp_bb_ov("i,e") * u0 * u1_bb_vo("a,i") * m1_xb_Lv("I,a");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmal1_xb_Lv += -1.000000 d-_aa(j,i) u0 u1_aa(a,j) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o2v1: 1, o0v1L1: 1, o2v0: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xb_Lv("I,e") -= dp_aa_oo("j,i") * u0 * u1_aa_vo("a,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -1.000000 d-_bb(j,i) u0 u1_bb(a,j) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o2v1: 1, o0v1L1: 1, o2v0: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xb_Lv("I,e") -= dp_bb_oo("j,i") * u0 * u1_bb_vo("a,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 d-_aa(a,b) u0 u1_aa(b,i) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1, o0v2: 1 | mem: o0v1L1: 2, o0v2: 1, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += dp_aa_vv("a,b") * u0 * u1_aa_vo("b,i") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 d-_bb(a,b) u0 u1_bb(b,i) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1, o0v2: 1 | mem: o0v1L1: 2, o0v2: 1, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += dp_bb_vv("a,b") * u0 * u1_bb_vo("b,i") * m2_xbbb_Lovv("I,i,a,e");
        }

        if (include_u0_ && include_u2_) {

            // sigmal1_xb_Lv += -1.000000 d-_aa(j,b) u0 u2_aaaa(b,a,i,j) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 2, 
            sigmal1_xb_Lv("I,e") -= dp_aa_ov("j,b") * u0 * u2_aaaa_vvoo("b,a,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 d-_aa(j,b) u0 u2_abab(b,a,j,i) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 2, 
            sigmal1_xb_Lv("I,e") += dp_aa_ov("j,b") * u0 * u2_abab_vvoo("b,a,j,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 d-_bb(j,b) u0 u2_abab(a,b,i,j) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 2, 
            sigmal1_xb_Lv("I,e") += dp_bb_ov("j,b") * u0 * u2_abab_vvoo("a,b,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -1.000000 d-_bb(j,b) u0 u2_bbbb(b,a,i,j) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 2, 
            sigmal1_xb_Lv("I,e") -= dp_bb_ov("j,b") * u0 * u2_bbbb_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -0.500000 d-_bb(j,e) u0 u2_abab(b,a,i,j) m2_aab(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1, o1v1: 1 | mem: o1v3: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * dp_bb_ov("j,e") * u0 * u2_abab_vvoo("b,a,i,j") * m2_xaab_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += -0.500000 d-_bb(j,e) u0 u2_bbbb(b,a,i,j) m2_bbb(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1, o1v1: 1 | mem: o1v3: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * dp_bb_ov("j,e") * u0 * u2_bbbb_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,i,b,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmal1_xb_Lv += 1.000000 d-_aa(j,b) u1_aa(b,j) u1_aa(a,i) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += scalar0 * u1_aa_vo("a,i") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 d-_aa(j,b) u1_aa(b,j) u1_bb(a,i) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += scalar0 * u1_bb_vo("a,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 d-_bb(j,b) u1_bb(b,j) u1_aa(a,i) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += scalar1 * u1_aa_vo("a,i") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 d-_bb(j,b) u1_bb(b,j) u1_bb(a,i) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += scalar1 * u1_bb_vo("a,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -2.000000 d-_aa(j,b) u1_aa(b,i) u1_aa(a,j) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o2v1: 2, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xb_Lv("I,e") -= 2.000000 * dp_aa_ov("j,b") * u1_aa_vo("b,i") * u1_aa_vo("a,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -2.000000 d-_bb(j,b) u1_bb(b,i) u1_bb(a,j) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o2v1: 2, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xb_Lv("I,e") -= 2.000000 * dp_bb_ov("j,b") * u1_bb_vo("b,i") * u1_bb_vo("a,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 d-_bb(j,e) u1_bb(b,j) u1_bb(a,i) m2_bbb(i,b,a) 
            // flops: o1v3L1: 1, o1v3: 1, o1v2: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, o0v2: 1, 
            sigmal1_xb_Lv("I,e") += dp_bb_ov("j,e") * u1_bb_vo("b,j") * u1_bb_vo("a,i") * m2_xbbb_Lovv("I,i,b,a");
        }

        {

            // sigmal1_xb_Lv += -0.500000 <j,i||j,i>_aaaa l1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * scalar9 * l1_xb_Lv("I,e");

            // sigmal1_xb_Lv += -0.500000 <j,i||j,i>_abab l1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * scalar10 * l1_xb_Lv("I,e");

            // sigmal1_xb_Lv += -0.500000 <i,j||i,j>_abab l1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * scalar10 * l1_xb_Lv("I,e");

            // sigmal1_xb_Lv += -0.500000 <j,i||j,i>_bbbb l1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * scalar11 * l1_xb_Lv("I,e");

            // sigmal1_xb_Lv += -0.500000 <b,a||i,e>_abab l2_aab(i,b,a) 
            // flops: o1v3L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["abba_vvvo"]("b,a,e,i") * l2_xaab_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += 0.500000 <b,a||e,i>_bbbb l2_bbb(i,b,a) 
            // flops: o1v3L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["bbbb_vvvo"]("b,a,e,i") * l2_xbbb_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += 0.250000 <j,i||a,b>_aaaa l1_b(e) t2_aaaa(a,b,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.250000 * scalar2 * l1_xb_Lv("I,e");

            // sigmal1_xb_Lv += 0.250000 <j,i||a,b>_abab l1_b(e) t2_abab(a,b,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.250000 * scalar3 * l1_xb_Lv("I,e");

            // sigmal1_xb_Lv += 0.250000 <i,j||a,b>_abab l1_b(e) t2_abab(a,b,i,j) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.250000 * scalar3 * l1_xb_Lv("I,e");

            // sigmal1_xb_Lv += 0.250000 <j,i||b,a>_abab l1_b(e) t2_abab(b,a,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.250000 * scalar3 * l1_xb_Lv("I,e");

            // sigmal1_xb_Lv += 0.250000 <i,j||b,a>_abab l1_b(e) t2_abab(b,a,i,j) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.250000 * scalar3 * l1_xb_Lv("I,e");

            // sigmal1_xb_Lv += 0.250000 <j,i||a,b>_bbbb l1_b(e) t2_bbbb(a,b,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.250000 * scalar4 * l1_xb_Lv("I,e");

            // sigmal1_xb_Lv += -0.500000 <j,i||b,e>_abab l1_b(a) t2_abab(b,a,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("j,i,b,e") * l1_xb_Lv("I,a") * t2_abab_vvoo("b,a,j,i");

            // sigmal1_xb_Lv += -0.500000 <i,j||b,e>_abab l1_b(a) t2_abab(b,a,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("i,j,b,e") * l1_xb_Lv("I,a") * t2_abab_vvoo("b,a,i,j");

            // sigmal1_xb_Lv += -0.500000 <j,i||b,e>_bbbb l1_b(a) t2_bbbb(b,a,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["bbbb_oovv"]("j,i,b,e") * l1_xb_Lv("I,a") * t2_bbbb_vvoo("b,a,j,i");

            // sigmal1_xb_Lv += 0.500000 <k,j||b,i>_aaaa l2_aab(i,a,e) t2_aaaa(b,a,k,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["aaaa_oovo"]("k,j,b,i") * l2_xaab_Lovv("I,i,a,e") * t2_aaaa_vvoo("b,a,k,j");

            // sigmal1_xb_Lv += 0.500000 <k,j||i,b>_abab l2_aab(i,a,e) t2_abab(a,b,k,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["abba_oovo"]("k,j,b,i") * l2_xaab_Lovv("I,i,a,e") * t2_abab_vvoo("a,b,k,j");

            // sigmal1_xb_Lv += 0.500000 <j,k||i,b>_abab l2_aab(i,a,e) t2_abab(a,b,j,k) 
            // flops: o3v3L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["abba_oovo"]("j,k,b,i") * l2_xaab_Lovv("I,i,a,e") * t2_abab_vvoo("a,b,j,k");

            // sigmal1_xb_Lv += 0.500000 <k,j||b,i>_abab l2_bbb(i,a,e) t2_abab(b,a,k,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["abab_oovo"]("k,j,b,i") * l2_xbbb_Lovv("I,i,a,e") * t2_abab_vvoo("b,a,k,j");

            // sigmal1_xb_Lv += 0.500000 <j,k||b,i>_abab l2_bbb(i,a,e) t2_abab(b,a,j,k) 
            // flops: o3v3L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["abab_oovo"]("j,k,b,i") * l2_xbbb_Lovv("I,i,a,e") * t2_abab_vvoo("b,a,j,k");

            // sigmal1_xb_Lv += 0.500000 <k,j||b,i>_bbbb l2_bbb(i,a,e) t2_bbbb(b,a,k,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["bbbb_oovo"]("k,j,b,i") * l2_xbbb_Lovv("I,i,a,e") * t2_bbbb_vvoo("b,a,k,j");

            // sigmal1_xb_Lv += -0.250000 <k,j||i,e>_abab l2_aab(i,b,a) t2_abab(b,a,k,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.250000 * V_blks_["abba_oovo"]("k,j,e,i") * l2_xaab_Lovv("I,i,b,a") * t2_abab_vvoo("b,a,k,j");

            // sigmal1_xb_Lv += -0.250000 <j,k||i,e>_abab l2_aab(i,b,a) t2_abab(b,a,j,k) 
            // flops: o3v3L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.250000 * V_blks_["abba_oovo"]("j,k,e,i") * l2_xaab_Lovv("I,i,b,a") * t2_abab_vvoo("b,a,j,k");

            // sigmal1_xb_Lv += 0.250000 <k,j||e,i>_bbbb l2_bbb(i,b,a) t2_bbbb(b,a,k,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.250000 * V_blks_["bbbb_oovo"]("k,j,e,i") * l2_xbbb_Lovv("I,i,b,a") * t2_bbbb_vvoo("b,a,k,j");

            // sigmal1_xb_Lv += 0.500000 <j,a||b,c>_aaaa l2_aab(i,a,e) t2_aaaa(b,c,i,j) 
            // flops: o2v4L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["aaaa_vovv"]("a,j,b,c") * l2_xaab_Lovv("I,i,a,e") * t2_aaaa_vvoo("b,c,i,j");

            // sigmal1_xb_Lv += -0.500000 <a,j||b,c>_abab l2_aab(i,a,e) t2_abab(b,c,i,j) 
            // flops: o2v4L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_vovv"]("a,j,b,c") * l2_xaab_Lovv("I,i,a,e") * t2_abab_vvoo("b,c,i,j");

            // sigmal1_xb_Lv += -0.500000 <a,j||c,b>_abab l2_aab(i,a,e) t2_abab(c,b,i,j) 
            // flops: o2v4L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_vovv"]("a,j,c,b") * l2_xaab_Lovv("I,i,a,e") * t2_abab_vvoo("c,b,i,j");

            // sigmal1_xb_Lv += -0.500000 <j,a||b,c>_abab l2_bbb(i,a,e) t2_abab(b,c,j,i) 
            // flops: o2v4L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["baab_vovv"]("a,j,b,c") * l2_xbbb_Lovv("I,i,a,e") * t2_abab_vvoo("b,c,j,i");

            // sigmal1_xb_Lv += -0.500000 <j,a||c,b>_abab l2_bbb(i,a,e) t2_abab(c,b,j,i) 
            // flops: o2v4L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["baab_vovv"]("a,j,c,b") * l2_xbbb_Lovv("I,i,a,e") * t2_abab_vvoo("c,b,j,i");

            // sigmal1_xb_Lv += 0.500000 <j,a||b,c>_bbbb l2_bbb(i,a,e) t2_bbbb(b,c,i,j) 
            // flops: o2v4L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["bbbb_vovv"]("a,j,b,c") * l2_xbbb_Lovv("I,i,a,e") * t2_bbbb_vvoo("b,c,i,j");

            // sigmal1_xb_Lv += 1.000000 <b,j||c,e>_abab l2_aab(i,b,a) t2_abab(c,a,i,j) 
            // flops: o2v4L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += V_blks_["abab_vovv"]("b,j,c,e") * l2_xaab_Lovv("I,i,b,a") * t2_abab_vvoo("c,a,i,j");

            // sigmal1_xb_Lv += 1.000000 <j,b||c,e>_abab l2_bbb(i,b,a) t2_abab(c,a,j,i) 
            // flops: o2v4L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= V_blks_["baab_vovv"]("b,j,c,e") * l2_xbbb_Lovv("I,i,b,a") * t2_abab_vvoo("c,a,j,i");

            // sigmal1_xb_Lv += -1.000000 <j,b||c,e>_bbbb l2_bbb(i,b,a) t2_bbbb(c,a,i,j) 
            // flops: o2v4L1: 1, o2v3L1: 1, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += V_blks_["bbbb_vovv"]("b,j,c,e") * l2_xbbb_Lovv("I,i,b,a") * t2_bbbb_vvoo("c,a,i,j");
        }

        if (include_u1_) {

            // sigmal1_xb_Lv += 1.000000 <i,a||b,e>_abab u1_aa(b,i) m1_b(a) 
            // flops: o1v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmal1_xb_Lv("I,e") -= V_blks_["baab_vovv"]("a,i,b,e") * u1_aa_vo("b,i") * m1_xb_Lv("I,a");

            // sigmal1_xb_Lv += 1.000000 <i,a||b,e>_bbbb u1_bb(b,i) m1_b(a) 
            // flops: o1v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmal1_xb_Lv("I,e") -= V_blks_["bbbb_vovv"]("a,i,b,e") * u1_bb_vo("b,i") * m1_xb_Lv("I,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmal1_xb_Lv += -1.000000 <j,a||b,i>_aaaa u1_aa(b,j) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += V_blks_["aaaa_vovo"]("a,j,b,i") * u1_aa_vo("b,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -1.000000 <j,a||b,i>_abab u1_aa(b,j) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += V_blks_["baab_vovo"]("a,j,b,i") * u1_aa_vo("b,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -1.000000 <a,j||i,b>_abab u1_bb(b,j) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += V_blks_["abba_vovo"]("a,j,b,i") * u1_bb_vo("b,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -1.000000 <j,a||b,i>_bbbb u1_bb(b,j) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += V_blks_["bbbb_vovo"]("a,j,b,i") * u1_bb_vo("b,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 <b,j||i,e>_abab u1_bb(a,j) m2_aab(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= V_blks_["abba_vovo"]("b,j,e,i") * u1_bb_vo("a,j") * m2_xaab_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += 1.000000 <j,b||e,i>_bbbb u1_bb(a,j) m2_bbb(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= V_blks_["bbbb_vovo"]("b,j,e,i") * u1_bb_vo("a,j") * m2_xbbb_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += -0.500000 <b,a||c,e>_abab u1_aa(c,i) m2_aab(i,b,a) 
            // flops: o1v3L1: 1, o1v4: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_vvvv"]("b,a,c,e") * u1_aa_vo("c,i") * m2_xaab_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += -0.500000 <b,a||c,e>_bbbb u1_bb(c,i) m2_bbb(i,b,a) 
            // flops: o1v3L1: 1, o1v4: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["bbbb_vvvv"]("b,a,c,e") * u1_bb_vo("c,i") * m2_xbbb_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += 0.250000 <j,i||a,b>_aaaa u2_aaaa(a,b,j,i) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.250000 * scalar14 * m1_xb_Lv("I,e");

            // sigmal1_xb_Lv += 0.250000 <j,i||a,b>_abab u2_abab(a,b,j,i) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.250000 * scalar15 * m1_xb_Lv("I,e");

            // sigmal1_xb_Lv += 0.250000 <i,j||a,b>_abab u2_abab(a,b,i,j) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.250000 * scalar15 * m1_xb_Lv("I,e");

            // sigmal1_xb_Lv += 0.250000 <j,i||b,a>_abab u2_abab(b,a,j,i) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.250000 * scalar15 * m1_xb_Lv("I,e");

            // sigmal1_xb_Lv += 0.250000 <i,j||b,a>_abab u2_abab(b,a,i,j) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.250000 * scalar15 * m1_xb_Lv("I,e");

            // sigmal1_xb_Lv += 0.250000 <j,i||a,b>_bbbb u2_bbbb(a,b,j,i) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.250000 * scalar16 * m1_xb_Lv("I,e");

            // sigmal1_xb_Lv += -0.500000 <j,i||b,e>_abab u2_abab(b,a,j,i) m1_b(a) 
            // flops: o2v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("j,i,b,e") * u2_abab_vvoo("b,a,j,i") * m1_xb_Lv("I,a");

            // sigmal1_xb_Lv += -0.500000 <i,j||b,e>_abab u2_abab(b,a,i,j) m1_b(a) 
            // flops: o2v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("i,j,b,e") * u2_abab_vvoo("b,a,i,j") * m1_xb_Lv("I,a");

            // sigmal1_xb_Lv += -0.500000 <j,i||b,e>_bbbb u2_bbbb(b,a,j,i) m1_b(a) 
            // flops: o2v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["bbbb_oovv"]("j,i,b,e") * u2_bbbb_vvoo("b,a,j,i") * m1_xb_Lv("I,a");
        }

        if (include_u2_) {

            // sigmal1_xb_Lv += 0.500000 <k,j||b,i>_aaaa u2_aaaa(b,a,k,j) m2_aab(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["aaaa_oovo"]("k,j,b,i") * u2_aaaa_vvoo("b,a,k,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 0.500000 <k,j||b,i>_abab u2_abab(b,a,k,j) m2_bbb(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["abab_oovo"]("k,j,b,i") * u2_abab_vvoo("b,a,k,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 0.500000 <j,k||b,i>_abab u2_abab(b,a,j,k) m2_bbb(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["abab_oovo"]("j,k,b,i") * u2_abab_vvoo("b,a,j,k") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 0.500000 <k,j||i,b>_abab u2_abab(a,b,k,j) m2_aab(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["abba_oovo"]("k,j,b,i") * u2_abab_vvoo("a,b,k,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 0.500000 <j,k||i,b>_abab u2_abab(a,b,j,k) m2_aab(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["abba_oovo"]("j,k,b,i") * u2_abab_vvoo("a,b,j,k") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 0.500000 <k,j||b,i>_bbbb u2_bbbb(b,a,k,j) m2_bbb(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["bbbb_oovo"]("k,j,b,i") * u2_bbbb_vvoo("b,a,k,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -0.250000 <k,j||i,e>_abab u2_abab(b,a,k,j) m2_aab(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.250000 * V_blks_["abba_oovo"]("k,j,e,i") * u2_abab_vvoo("b,a,k,j") * m2_xaab_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += -0.250000 <j,k||i,e>_abab u2_abab(b,a,j,k) m2_aab(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.250000 * V_blks_["abba_oovo"]("j,k,e,i") * u2_abab_vvoo("b,a,j,k") * m2_xaab_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += 0.250000 <k,j||e,i>_bbbb u2_bbbb(b,a,k,j) m2_bbb(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.250000 * V_blks_["bbbb_oovo"]("k,j,e,i") * u2_bbbb_vvoo("b,a,k,j") * m2_xbbb_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += 0.500000 <j,a||b,c>_aaaa u2_aaaa(b,c,i,j) m2_aab(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["aaaa_vovv"]("a,j,b,c") * u2_aaaa_vvoo("b,c,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -0.500000 <a,j||b,c>_abab u2_abab(b,c,i,j) m2_aab(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_vovv"]("a,j,b,c") * u2_abab_vvoo("b,c,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -0.500000 <j,a||b,c>_abab u2_abab(b,c,j,i) m2_bbb(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["baab_vovv"]("a,j,b,c") * u2_abab_vvoo("b,c,j,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -0.500000 <a,j||c,b>_abab u2_abab(c,b,i,j) m2_aab(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_vovv"]("a,j,c,b") * u2_abab_vvoo("c,b,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -0.500000 <j,a||c,b>_abab u2_abab(c,b,j,i) m2_bbb(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["baab_vovv"]("a,j,c,b") * u2_abab_vvoo("c,b,j,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 0.500000 <j,a||b,c>_bbbb u2_bbbb(b,c,i,j) m2_bbb(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") -= 0.500000 * V_blks_["bbbb_vovv"]("a,j,b,c") * u2_bbbb_vvoo("b,c,i,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 <b,j||c,e>_abab u2_abab(c,a,i,j) m2_aab(i,b,a) 
            // flops: o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += V_blks_["abab_vovv"]("b,j,c,e") * u2_abab_vvoo("c,a,i,j") * m2_xaab_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += 1.000000 <j,b||c,e>_abab u2_abab(c,a,j,i) m2_bbb(i,b,a) 
            // flops: o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= V_blks_["baab_vovv"]("b,j,c,e") * u2_abab_vvoo("c,a,j,i") * m2_xbbb_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += -1.000000 <j,b||c,e>_bbbb u2_bbbb(c,a,i,j) m2_bbb(i,b,a) 
            // flops: o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += V_blks_["bbbb_vovv"]("b,j,c,e") * u2_bbbb_vvoo("c,a,i,j") * m2_xbbb_Lovv("I,i,b,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmal1_xb_Lv += 0.500000 <k,j||b,c>_aaaa t2_aaaa(b,c,i,j) u1_aa(a,k) m2_aab(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["aaaa_oovv"]("k,j,b,c") * t2_aaaa_vvoo("b,c,i,j") * u1_aa_vo("a,k") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 0.500000 <k,j||b,c>_abab t2_abab(b,c,i,j) u1_aa(a,k) m2_aab(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("k,j,b,c") * t2_abab_vvoo("b,c,i,j") * u1_aa_vo("a,k") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 0.500000 <j,k||b,c>_abab t2_abab(b,c,j,i) u1_bb(a,k) m2_bbb(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("j,k,b,c") * t2_abab_vvoo("b,c,j,i") * u1_bb_vo("a,k") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 0.500000 <k,j||c,b>_abab t2_abab(c,b,i,j) u1_aa(a,k) m2_aab(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("k,j,c,b") * t2_abab_vvoo("c,b,i,j") * u1_aa_vo("a,k") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 0.500000 <j,k||c,b>_abab t2_abab(c,b,j,i) u1_bb(a,k) m2_bbb(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("j,k,c,b") * t2_abab_vvoo("c,b,j,i") * u1_bb_vo("a,k") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 0.500000 <k,j||b,c>_bbbb t2_bbbb(b,c,i,j) u1_bb(a,k) m2_bbb(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, o2v0: 1, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["bbbb_oovv"]("k,j,b,c") * t2_bbbb_vvoo("b,c,i,j") * u1_bb_vo("a,k") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 0.500000 <k,j||b,c>_aaaa t2_aaaa(b,a,k,j) u1_aa(c,i) m2_aab(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["aaaa_oovv"]("k,j,b,c") * t2_aaaa_vvoo("b,a,k,j") * u1_aa_vo("c,i") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 0.500000 <k,j||b,c>_abab t2_abab(b,a,k,j) u1_bb(c,i) m2_bbb(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("k,j,b,c") * t2_abab_vvoo("b,a,k,j") * u1_bb_vo("c,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 0.500000 <j,k||b,c>_abab t2_abab(b,a,j,k) u1_bb(c,i) m2_bbb(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("j,k,b,c") * t2_abab_vvoo("b,a,j,k") * u1_bb_vo("c,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 0.500000 <k,j||c,b>_abab t2_abab(a,b,k,j) u1_aa(c,i) m2_aab(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("k,j,c,b") * t2_abab_vvoo("a,b,k,j") * u1_aa_vo("c,i") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 0.500000 <j,k||c,b>_abab t2_abab(a,b,j,k) u1_aa(c,i) m2_aab(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("j,k,c,b") * t2_abab_vvoo("a,b,j,k") * u1_aa_vo("c,i") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 0.500000 <k,j||b,c>_bbbb t2_bbbb(b,a,k,j) u1_bb(c,i) m2_bbb(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["bbbb_oovv"]("k,j,b,c") * t2_bbbb_vvoo("b,a,k,j") * u1_bb_vo("c,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -1.000000 <k,j||b,c>_aaaa t2_aaaa(b,a,i,j) u1_aa(c,k) m2_aab(i,a,e) 
            // flops: o3v3: 1, o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o2v2: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") -= V_blks_["aaaa_oovv"]("k,j,b,c") * t2_aaaa_vvoo("b,a,i,j") * u1_aa_vo("c,k") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 <j,k||b,c>_abab t2_aaaa(b,a,i,j) u1_bb(c,k) m2_aab(i,a,e) 
            // flops: o3v3: 1, o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o2v2: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += V_blks_["abab_oovv"]("j,k,b,c") * t2_aaaa_vvoo("b,a,i,j") * u1_bb_vo("c,k") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 <k,j||b,c>_aaaa t2_abab(b,a,j,i) u1_aa(c,k) m2_bbb(i,a,e) 
            // flops: o3v3: 1, o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o2v2: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += V_blks_["aaaa_oovv"]("k,j,b,c") * t2_abab_vvoo("b,a,j,i") * u1_aa_vo("c,k") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -1.000000 <j,k||b,c>_abab t2_abab(b,a,j,i) u1_bb(c,k) m2_bbb(i,a,e) 
            // flops: o3v3: 1, o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o2v2: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") -= V_blks_["abab_oovv"]("j,k,b,c") * t2_abab_vvoo("b,a,j,i") * u1_bb_vo("c,k") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -1.000000 <k,j||c,b>_abab t2_abab(a,b,i,j) u1_aa(c,k) m2_aab(i,a,e) 
            // flops: o3v3: 1, o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o2v2: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") -= V_blks_["abab_oovv"]("k,j,c,b") * t2_abab_vvoo("a,b,i,j") * u1_aa_vo("c,k") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 <k,j||b,c>_bbbb t2_abab(a,b,i,j) u1_bb(c,k) m2_aab(i,a,e) 
            // flops: o3v3: 1, o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o2v2: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += V_blks_["bbbb_oovv"]("k,j,b,c") * t2_abab_vvoo("a,b,i,j") * u1_bb_vo("c,k") * m2_xaab_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += 1.000000 <k,j||c,b>_abab t2_bbbb(b,a,i,j) u1_aa(c,k) m2_bbb(i,a,e) 
            // flops: o3v3: 1, o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o2v2: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") += V_blks_["abab_oovv"]("k,j,c,b") * t2_bbbb_vvoo("b,a,i,j") * u1_aa_vo("c,k") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -1.000000 <k,j||b,c>_bbbb t2_bbbb(b,a,i,j) u1_bb(c,k) m2_bbb(i,a,e) 
            // flops: o3v3: 1, o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o2v2: 1, o0v1L1: 2, o1v1: 1, 
            sigmal1_xb_Lv("I,e") -= V_blks_["bbbb_oovv"]("k,j,b,c") * t2_bbbb_vvoo("b,a,i,j") * u1_bb_vo("c,k") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal1_xb_Lv += -1.000000 <j,k||c,e>_abab t2_aaaa(c,b,i,j) u1_bb(a,k) m2_aab(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o2v2: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= V_blks_["abab_oovv"]("j,k,c,e") * t2_aaaa_vvoo("c,b,i,j") * u1_bb_vo("a,k") * m2_xaab_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += 1.000000 <j,k||c,e>_abab t2_abab(c,b,j,i) u1_bb(a,k) m2_bbb(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o2v2: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += V_blks_["abab_oovv"]("j,k,c,e") * t2_abab_vvoo("c,b,j,i") * u1_bb_vo("a,k") * m2_xbbb_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += -1.000000 <k,j||c,e>_bbbb t2_abab(b,c,i,j) u1_bb(a,k) m2_aab(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o2v2: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= V_blks_["bbbb_oovv"]("k,j,c,e") * t2_abab_vvoo("b,c,i,j") * u1_bb_vo("a,k") * m2_xaab_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += 1.000000 <k,j||c,e>_bbbb t2_bbbb(c,b,i,j) u1_bb(a,k) m2_bbb(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o2v2: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += V_blks_["bbbb_oovv"]("k,j,c,e") * t2_bbbb_vvoo("c,b,i,j") * u1_bb_vo("a,k") * m2_xbbb_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += -0.250000 <k,j||c,e>_abab t2_abab(b,a,k,j) u1_aa(c,i) m2_aab(i,b,a) 
            // flops: o2v4: 1, o1v3L1: 1, o1v4: 1, o0v1L1: 1 | mem: o0v4: 1, o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= 0.250000 * V_blks_["abab_oovv"]("k,j,c,e") * t2_abab_vvoo("b,a,k,j") * u1_aa_vo("c,i") * m2_xaab_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += -0.250000 <j,k||c,e>_abab t2_abab(b,a,j,k) u1_aa(c,i) m2_aab(i,b,a) 
            // flops: o2v4: 1, o1v3L1: 1, o1v4: 1, o0v1L1: 1 | mem: o0v4: 1, o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= 0.250000 * V_blks_["abab_oovv"]("j,k,c,e") * t2_abab_vvoo("b,a,j,k") * u1_aa_vo("c,i") * m2_xaab_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += -0.250000 <k,j||c,e>_bbbb t2_bbbb(b,a,k,j) u1_bb(c,i) m2_bbb(i,b,a) 
            // flops: o2v4: 1, o1v3L1: 1, o1v4: 1, o0v1L1: 1 | mem: o0v4: 1, o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") -= 0.250000 * V_blks_["bbbb_oovv"]("k,j,c,e") * t2_bbbb_vvoo("b,a,k,j") * u1_bb_vo("c,i") * m2_xbbb_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += 0.500000 <k,j||c,e>_abab t2_abab(b,a,i,j) u1_aa(c,k) m2_aab(i,b,a) 
            // flops: o3v4: 1, o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o2v4: 1, o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("k,j,c,e") * t2_abab_vvoo("b,a,i,j") * u1_aa_vo("c,k") * m2_xaab_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += 0.500000 <k,j||c,e>_bbbb t2_abab(b,a,i,j) u1_bb(c,k) m2_aab(i,b,a) 
            // flops: o3v4: 1, o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o2v4: 1, o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["bbbb_oovv"]("k,j,c,e") * t2_abab_vvoo("b,a,i,j") * u1_bb_vo("c,k") * m2_xaab_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += 0.500000 <k,j||c,e>_abab t2_bbbb(b,a,i,j) u1_aa(c,k) m2_bbb(i,b,a) 
            // flops: o3v4: 1, o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o2v4: 1, o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["abab_oovv"]("k,j,c,e") * t2_bbbb_vvoo("b,a,i,j") * u1_aa_vo("c,k") * m2_xbbb_Lovv("I,i,b,a");

            // sigmal1_xb_Lv += 0.500000 <k,j||c,e>_bbbb t2_bbbb(b,a,i,j) u1_bb(c,k) m2_bbb(i,b,a) 
            // flops: o3v4: 1, o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o2v4: 1, o1v3: 1, o0v1L1: 2, 
            sigmal1_xb_Lv("I,e") += 0.500000 * V_blks_["bbbb_oovv"]("k,j,c,e") * t2_bbbb_vvoo("b,a,i,j") * u1_bb_vo("c,k") * m2_xbbb_Lovv("I,i,b,a");
        }




/// ****** sigmal2_aaa(e,f,n) ****** ///



        if (include_u0_ && include_u2_) {

            // sigmal2_xaaa_Lvvo += 1.000000 u0 m2_aaa(n,e,f) w0 
            // flops: o1v2L1: 2, o1v2: 1 | mem: o1v2L1: 1, o1v2: 2, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += u0 * m2_xaaa_Lovv("I,n,e,f") * w0;
        }

        {

            // sigmal2_xaaa_Lvvo += -1.000000 P(e,f) f_aa(n,e) l1_a(f) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = F_blks_["aa_ov"]("n,e") * l1_xa_Lv("I,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += 1.000000 f_aa(i,i) l2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += scalar5 * l2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += 1.000000 f_bb(i,i) l2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += scalar6 * l2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += -1.000000 f_aa(n,i) l2_aaa(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= F_blks_["aa_oo"]("n,i") * l2_xaaa_Lovv("I,i,e,f");

            // sigmal2_xaaa_Lvvo += 1.000000 P(e,f) f_aa(a,e) l2_aaa(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = F_blks_["aa_vv"]("a,e") * l2_xaaa_Lovv("I,n,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u1_ && include_u2_) {

            // sigmal2_xaaa_Lvvo += 1.000000 f_aa(i,a) u1_aa(a,i) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += scalar12 * m2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += 1.000000 f_bb(i,a) u1_bb(a,i) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += scalar13 * m2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += -1.000000 f_aa(n,a) u1_aa(a,i) m2_aaa(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1, o2v1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= F_blks_["aa_ov"]("n,a") * u1_aa_vo("a,i") * m2_xaaa_Lovv("I,i,e,f");

            // sigmal2_xaaa_Lvvo += -1.000000 P(e,f) f_aa(i,e) u1_aa(a,i) m2_aaa(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2, o1v2: 1 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = F_blks_["aa_ov"]("i,e") * u1_aa_vo("a,i") * m2_xaaa_Lovv("I,n,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u1_) {

            // sigmal2_xaaa_Lvvo += 1.000000 P(e,f) d+_aa(n,e) m1_a(f) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("n,e") * m1_xa_Lv("I,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmal2_xaaa_Lvvo += -1.000000 d+_aa(i,i) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= scalar17 * m2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += -1.000000 d+_bb(i,i) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= scalar18 * m2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += 1.000000 d+_aa(n,i) m2_aaa(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += dp_aa_oo("n,i") * m2_xaaa_Lovv("I,i,e,f");

            // sigmal2_xaaa_Lvvo += -1.000000 P(e,f) d+_aa(a,e) m2_aaa(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_vv("a,e") * m2_xaaa_Lovv("I,n,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u0_) {

            // sigmal2_xaaa_Lvvo += 1.000000 P(e,f) d-_aa(n,e) l1_a(f) u0 
            // flops: o1v2L1: 4 | mem: o1v2L1: 2, o1v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("n,e") * l1_xa_Lv("I,f") * u0;
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += -1.000000 d-_aa(i,i) l2_aaa(n,e,f) u0 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= scalar7 * u0 * l2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += -1.000000 d-_bb(i,i) l2_aaa(n,e,f) u0 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= scalar8 * u0 * l2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += 1.000000 d-_aa(n,i) l2_aaa(i,e,f) u0 
            // flops: o2v2L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += dp_aa_oo("n,i") * l2_xaaa_Lovv("I,i,e,f") * u0;

            // sigmal2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(a,e) l2_aaa(n,a,f) u0 
            // flops: o1v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_vv("a,e") * l2_xaaa_Lovv("I,n,a,f") * u0;
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u1_) {

            // sigmal2_xaaa_Lvvo += -1.000000 d-_aa(i,a) l2_aaa(n,e,f) u1_aa(a,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= scalar0 * l2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += -1.000000 d-_bb(i,a) l2_aaa(n,e,f) u1_bb(a,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= scalar1 * l2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += 1.000000 d-_aa(n,a) l2_aaa(i,e,f) u1_aa(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += dp_aa_ov("n,a") * l2_xaaa_Lovv("I,i,e,f") * u1_aa_vo("a,i");

            // sigmal2_xaaa_Lvvo += 1.000000 P(e,f) d-_aa(i,e) l2_aaa(n,a,f) u1_aa(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("i,e") * l2_xaaa_Lovv("I,n,a,f") * u1_aa_vo("a,i");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(n,e) l2_aaa(i,a,f) u1_aa(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("n,e") * l2_xaaa_Lovv("I,i,a,f") * u1_aa_vo("a,i");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(n,e) l2_bba(i,a,f) u1_bb(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("n,e") * l2_xbba_Lovv("I,i,a,f") * u1_bb_vo("a,i");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmal2_xaaa_Lvvo += -1.000000 d-_aa(i,a) u0 u1_aa(a,i) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= scalar0 * u0 * m2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += -1.000000 d-_bb(i,a) u0 u1_bb(a,i) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= scalar1 * u0 * m2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += 1.000000 d-_aa(n,a) u0 u1_aa(a,i) m2_aaa(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1, o2v1: 1, o1v1: 1 | mem: o1v2L1: 2, o1v1: 1, o2v0: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += dp_aa_ov("n,a") * u0 * u1_aa_vo("a,i") * m2_xaaa_Lovv("I,i,e,f");

            // sigmal2_xaaa_Lvvo += 1.000000 P(e,f) d-_aa(i,e) u0 u1_aa(a,i) m2_aaa(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2, o1v2: 1, o1v1: 1 | mem: o1v2L1: 2, o0v2: 1, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("i,e") * u0 * u1_aa_vo("a,i") * m2_xaaa_Lovv("I,n,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        {

            // sigmal2_xaaa_Lvvo += -0.500000 <j,i||j,i>_aaaa l2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * scalar9 * l2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += -0.500000 <j,i||j,i>_abab l2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * scalar10 * l2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += -0.500000 <i,j||i,j>_abab l2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * scalar10 * l2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += -0.500000 <j,i||j,i>_bbbb l2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * scalar11 * l2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += -1.000000 <n,a||e,f>_aaaa l1_a(a) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += V_blks_["aaaa_vovv"]("a,n,e,f") * l1_xa_Lv("I,a");

            // sigmal2_xaaa_Lvvo += 1.000000 P(e,f) <n,a||e,i>_aaaa l2_aaa(i,a,f) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_vovo"]("a,n,e,i") * l2_xaaa_Lovv("I,i,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += 1.000000 P(e,f) <n,a||e,i>_abab l2_bba(i,a,f) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["baab_vovo"]("a,n,e,i") * l2_xbba_Lovv("I,i,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += 0.500000 <b,a||e,f>_aaaa l2_aaa(n,b,a) 
            // flops: o1v4L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["aaaa_vvvv"]("b,a,e,f") * l2_xaaa_Lovv("I,n,b,a");

            // sigmal2_xaaa_Lvvo += 0.250000 <j,i||a,b>_aaaa l2_aaa(n,e,f) t2_aaaa(a,b,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar2 * l2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += 0.250000 <j,i||a,b>_abab l2_aaa(n,e,f) t2_abab(a,b,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar3 * l2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += 0.250000 <i,j||a,b>_abab l2_aaa(n,e,f) t2_abab(a,b,i,j) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar3 * l2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += 0.250000 <j,i||b,a>_abab l2_aaa(n,e,f) t2_abab(b,a,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar3 * l2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += 0.250000 <i,j||b,a>_abab l2_aaa(n,e,f) t2_abab(b,a,i,j) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar3 * l2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += 0.250000 <j,i||a,b>_bbbb l2_aaa(n,e,f) t2_bbbb(a,b,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar4 * l2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += -0.500000 <n,j||a,b>_aaaa l2_aaa(i,e,f) t2_aaaa(a,b,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovv"]("n,j,a,b") * l2_xaaa_Lovv("I,i,e,f") * t2_aaaa_vvoo("a,b,i,j");

            // sigmal2_xaaa_Lvvo += -0.500000 <n,j||a,b>_abab l2_aaa(i,e,f) t2_abab(a,b,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("n,j,a,b") * l2_xaaa_Lovv("I,i,e,f") * t2_abab_vvoo("a,b,i,j");

            // sigmal2_xaaa_Lvvo += -0.500000 <n,j||b,a>_abab l2_aaa(i,e,f) t2_abab(b,a,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("n,j,b,a") * l2_xaaa_Lovv("I,i,e,f") * t2_abab_vvoo("b,a,i,j");

            // sigmal2_xaaa_Lvvo += -0.500000 P(e,f) <j,i||b,e>_aaaa l2_aaa(n,a,f) t2_aaaa(b,a,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 2 | mem: o3v4L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["aaaa_oovv"]("j,i,b,e") * l2_xaaa_Lovv("I,n,a,f") * t2_aaaa_vvoo("b,a,j,i");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += -0.500000 P(e,f) <j,i||e,b>_abab l2_aaa(n,a,f) t2_abab(a,b,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 2 | mem: o3v4L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("j,i,e,b") * l2_xaaa_Lovv("I,n,a,f") * t2_abab_vvoo("a,b,j,i");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += -0.500000 P(e,f) <i,j||e,b>_abab l2_aaa(n,a,f) t2_abab(a,b,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 2 | mem: o3v4L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("i,j,e,b") * l2_xaaa_Lovv("I,n,a,f") * t2_abab_vvoo("a,b,i,j");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += 1.000000 P(e,f) <n,j||b,e>_aaaa l2_aaa(i,a,f) t2_aaaa(b,a,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 2 | mem: o3v4L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_oovv"]("n,j,b,e") * l2_xaaa_Lovv("I,i,a,f") * t2_aaaa_vvoo("b,a,i,j");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += 1.000000 P(e,f) <n,j||e,b>_abab l2_aaa(i,a,f) t2_abab(a,b,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 2 | mem: o3v4L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("n,j,e,b") * l2_xaaa_Lovv("I,i,a,f") * t2_abab_vvoo("a,b,i,j");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += -1.000000 P(e,f) <n,j||b,e>_aaaa l2_bba(i,a,f) t2_abab(b,a,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 2 | mem: o3v4L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_oovv"]("n,j,b,e") * l2_xbba_Lovv("I,i,a,f") * t2_abab_vvoo("b,a,j,i");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += -1.000000 P(e,f) <n,j||e,b>_abab l2_bba(i,a,f) t2_bbbb(b,a,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 2 | mem: o3v4L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("n,j,e,b") * l2_xbba_Lovv("I,i,a,f") * t2_bbbb_vvoo("b,a,i,j");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += 0.250000 <j,i||e,f>_aaaa l2_aaa(n,b,a) t2_aaaa(b,a,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += 0.250000 * V_blks_["aaaa_oovv"]("j,i,e,f") * l2_xaaa_Lovv("I,n,b,a") * t2_aaaa_vvoo("b,a,j,i");

            // sigmal2_xaaa_Lvvo += -0.500000 <n,j||e,f>_aaaa l2_aaa(i,b,a) t2_aaaa(b,a,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovv"]("n,j,e,f") * l2_xaaa_Lovv("I,i,b,a") * t2_aaaa_vvoo("b,a,i,j");

            // sigmal2_xaaa_Lvvo += -0.500000 <n,j||e,f>_aaaa l2_bba(i,b,a) t2_abab(a,b,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovv"]("n,j,e,f") * l2_xbba_Lovv("I,i,b,a") * t2_abab_vvoo("a,b,j,i");
        }

        if (include_u1_) {

            // sigmal2_xaaa_Lvvo += 1.000000 P(e,f) <n,i||a,e>_aaaa u1_aa(a,i) m1_a(f) 
            // flops: o1v2L1: 3, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_oovv"]("n,i,a,e") * u1_aa_vo("a,i") * m1_xa_Lv("I,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += -1.000000 P(e,f) <n,i||e,a>_abab u1_bb(a,i) m1_a(f) 
            // flops: o1v2L1: 3, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("n,i,e,a") * u1_bb_vo("a,i") * m1_xa_Lv("I,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += 1.000000 <n,i||e,f>_aaaa u1_aa(a,i) m1_a(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += V_blks_["aaaa_oovv"]("n,i,e,f") * u1_aa_vo("a,i") * m1_xa_Lv("I,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmal2_xaaa_Lvvo += 1.000000 <n,j||a,i>_aaaa u1_aa(a,j) m2_aaa(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1, o3v1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += V_blks_["aaaa_oovo"]("n,j,a,i") * u1_aa_vo("a,j") * m2_xaaa_Lovv("I,i,e,f");

            // sigmal2_xaaa_Lvvo += -1.000000 <n,j||i,a>_abab u1_bb(a,j) m2_aaa(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1, o3v1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += V_blks_["abba_oovo"]("n,j,a,i") * u1_bb_vo("a,j") * m2_xaaa_Lovv("I,i,e,f");

            // sigmal2_xaaa_Lvvo += -1.000000 P(e,f) <n,j||e,i>_aaaa u1_aa(a,j) m2_aaa(i,a,f) 
            // flops: o2v3L1: 1, o3v2: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_oovo"]("n,j,e,i") * u1_aa_vo("a,j") * m2_xaaa_Lovv("I,i,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += -1.000000 P(e,f) <n,j||e,i>_abab u1_bb(a,j) m2_bba(i,a,f) 
            // flops: o2v3L1: 1, o3v2: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_oovo"]("n,j,e,i") * u1_bb_vo("a,j") * m2_xbba_Lovv("I,i,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += 1.000000 P(e,f) <i,a||b,e>_aaaa u1_aa(b,i) m2_aaa(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2, o1v3: 1 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_vovv"]("a,i,b,e") * u1_aa_vo("b,i") * m2_xaaa_Lovv("I,n,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += 1.000000 P(e,f) <a,i||e,b>_abab u1_bb(b,i) m2_aaa(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2, o1v3: 1 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_vovv"]("a,i,e,b") * u1_bb_vo("b,i") * m2_xaaa_Lovv("I,n,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += -1.000000 P(e,f) <n,a||b,e>_aaaa u1_aa(b,i) m2_aaa(i,a,f) 
            // flops: o2v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_vovv"]("a,n,b,e") * u1_aa_vo("b,i") * m2_xaaa_Lovv("I,i,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += 1.000000 P(e,f) <n,a||e,b>_abab u1_bb(b,i) m2_bba(i,a,f) 
            // flops: o2v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["baab_vovv"]("a,n,e,b") * u1_bb_vo("b,i") * m2_xbba_Lovv("I,i,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += 1.000000 <i,b||e,f>_aaaa u1_aa(a,i) m2_aaa(n,b,a) 
            // flops: o1v4L1: 1, o1v4: 1, o1v2L1: 1 | mem: o0v4: 1, o1v2L1: 2, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= V_blks_["aaaa_vovv"]("b,i,e,f") * u1_aa_vo("a,i") * m2_xaaa_Lovv("I,n,b,a");
        }

        if (include_u2_) {

            // sigmal2_xaaa_Lvvo += 0.250000 <j,i||a,b>_aaaa u2_aaaa(a,b,j,i) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar14 * m2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += 0.250000 <j,i||a,b>_abab u2_abab(a,b,j,i) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar15 * m2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += 0.250000 <i,j||a,b>_abab u2_abab(a,b,i,j) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar15 * m2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += 0.250000 <j,i||b,a>_abab u2_abab(b,a,j,i) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar15 * m2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += 0.250000 <i,j||b,a>_abab u2_abab(b,a,i,j) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar15 * m2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += 0.250000 <j,i||a,b>_bbbb u2_bbbb(a,b,j,i) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar16 * m2_xaaa_Lovv("I,n,e,f");

            // sigmal2_xaaa_Lvvo += -0.500000 <n,j||a,b>_aaaa u2_aaaa(a,b,i,j) m2_aaa(i,e,f) 
            // flops: o2v2L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovv"]("n,j,a,b") * u2_aaaa_vvoo("a,b,i,j") * m2_xaaa_Lovv("I,i,e,f");

            // sigmal2_xaaa_Lvvo += -0.500000 <n,j||a,b>_abab u2_abab(a,b,i,j) m2_aaa(i,e,f) 
            // flops: o2v2L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("n,j,a,b") * u2_abab_vvoo("a,b,i,j") * m2_xaaa_Lovv("I,i,e,f");

            // sigmal2_xaaa_Lvvo += -0.500000 <n,j||b,a>_abab u2_abab(b,a,i,j) m2_aaa(i,e,f) 
            // flops: o2v2L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("n,j,b,a") * u2_abab_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,e,f");

            // sigmal2_xaaa_Lvvo += -0.500000 P(e,f) <j,i||b,e>_aaaa u2_aaaa(b,a,j,i) m2_aaa(n,a,f) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["aaaa_oovv"]("j,i,b,e") * u2_aaaa_vvoo("b,a,j,i") * m2_xaaa_Lovv("I,n,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += -0.500000 P(e,f) <j,i||e,b>_abab u2_abab(a,b,j,i) m2_aaa(n,a,f) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("j,i,e,b") * u2_abab_vvoo("a,b,j,i") * m2_xaaa_Lovv("I,n,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += -0.500000 P(e,f) <i,j||e,b>_abab u2_abab(a,b,i,j) m2_aaa(n,a,f) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("i,j,e,b") * u2_abab_vvoo("a,b,i,j") * m2_xaaa_Lovv("I,n,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += 1.000000 P(e,f) <n,j||b,e>_aaaa u2_aaaa(b,a,i,j) m2_aaa(i,a,f) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_oovv"]("n,j,b,e") * u2_aaaa_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += -1.000000 P(e,f) <n,j||b,e>_aaaa u2_abab(b,a,j,i) m2_bba(i,a,f) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_oovv"]("n,j,b,e") * u2_abab_vvoo("b,a,j,i") * m2_xbba_Lovv("I,i,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += 1.000000 P(e,f) <n,j||e,b>_abab u2_abab(a,b,i,j) m2_aaa(i,a,f) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("n,j,e,b") * u2_abab_vvoo("a,b,i,j") * m2_xaaa_Lovv("I,i,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += -1.000000 P(e,f) <n,j||e,b>_abab u2_bbbb(b,a,i,j) m2_bba(i,a,f) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("n,j,e,b") * u2_bbbb_vvoo("b,a,i,j") * m2_xbba_Lovv("I,i,a,f");
            sigmal2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmal2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmal2_xaaa_Lvvo += 0.250000 <j,i||e,f>_aaaa u2_aaaa(b,a,j,i) m2_aaa(n,b,a) 
            // flops: o1v4L1: 1, o2v4: 1, o1v2L1: 1 | mem: o0v4: 1, o1v2L1: 2, 
            sigmal2_xaaa_Lvvo("I,e,f,n") += 0.250000 * V_blks_["aaaa_oovv"]("j,i,e,f") * u2_aaaa_vvoo("b,a,j,i") * m2_xaaa_Lovv("I,n,b,a");

            // sigmal2_xaaa_Lvvo += -0.500000 <n,j||e,f>_aaaa u2_aaaa(b,a,i,j) m2_aaa(i,b,a) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovv"]("n,j,e,f") * u2_aaaa_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,b,a");

            // sigmal2_xaaa_Lvvo += -0.500000 <n,j||e,f>_aaaa u2_abab(a,b,j,i) m2_bba(i,b,a) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmal2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovv"]("n,j,e,f") * u2_abab_vvoo("a,b,j,i") * m2_xbba_Lovv("I,i,b,a");
        }




/// ****** sigmal2_abb(e,f,n) ****** ///



        {

            // sigmal2_xabb_Lvvo += 1.000000 f_bb(n,f) l1_a(e) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") += F_blks_["bb_ov"]("n,f") * l1_xa_Lv("I,e");

            // sigmal2_xabb_Lvvo += -1.000000 f_bb(a,f) l2_bba(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= F_blks_["bb_vv"]("a,f") * l2_xbba_Lovv("I,n,a,e");
        }

        if (include_u1_ && include_u2_) {

            // sigmal2_xabb_Lvvo += 1.000000 f_bb(i,f) u1_bb(a,i) m2_bba(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1, o1v2: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") += F_blks_["bb_ov"]("i,f") * u1_bb_vo("a,i") * m2_xbba_Lovv("I,n,a,e");
        }

        if (include_u1_) {

            // sigmal2_xabb_Lvvo += -1.000000 d+_bb(n,f) m1_a(e) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= dp_bb_ov("n,f") * m1_xa_Lv("I,e");
        }

        if (include_u2_) {

            // sigmal2_xabb_Lvvo += 1.000000 d+_bb(a,f) m2_bba(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") += dp_bb_vv("a,f") * m2_xbba_Lovv("I,n,a,e");
        }

        if (include_u0_) {

            // sigmal2_xabb_Lvvo += -1.000000 d-_bb(n,f) l1_a(e) u0 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, o1v2: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= dp_bb_ov("n,f") * l1_xa_Lv("I,e") * u0;

            // sigmal2_xabb_Lvvo += 1.000000 d-_bb(a,f) l2_bba(n,a,e) u0 
            // flops: o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v2: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") += dp_bb_vv("a,f") * l2_xbba_Lovv("I,n,a,e") * u0;
        }

        if (include_u1_) {

            // sigmal2_xabb_Lvvo += -1.000000 d-_bb(i,f) l2_bba(n,a,e) u1_bb(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= dp_bb_ov("i,f") * l2_xbba_Lovv("I,n,a,e") * u1_bb_vo("a,i");

            // sigmal2_xabb_Lvvo += 1.000000 d-_bb(n,f) l2_aaa(i,a,e) u1_aa(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") += dp_bb_ov("n,f") * l2_xaaa_Lovv("I,i,a,e") * u1_aa_vo("a,i");

            // sigmal2_xabb_Lvvo += 1.000000 d-_bb(n,f) l2_bba(i,a,e) u1_bb(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") += dp_bb_ov("n,f") * l2_xbba_Lovv("I,i,a,e") * u1_bb_vo("a,i");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmal2_xabb_Lvvo += -1.000000 d-_bb(i,f) u0 u1_bb(a,i) m2_bba(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1, o1v2: 1, o1v1: 1 | mem: o1v2L1: 2, o0v2: 1, o1v1: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= dp_bb_ov("i,f") * u0 * u1_bb_vo("a,i") * m2_xbba_Lovv("I,n,a,e");
        }

        {

            // sigmal2_xabb_Lvvo += 1.000000 <a,n||e,f>_abab l1_a(a) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") += V_blks_["abab_vovv"]("a,n,e,f") * l1_xa_Lv("I,a");

            // sigmal2_xabb_Lvvo += -1.000000 <a,n||i,f>_abab l2_aaa(i,a,e) 
            // flops: o2v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") += V_blks_["abba_vovo"]("a,n,f,i") * l2_xaaa_Lovv("I,i,a,e");

            // sigmal2_xabb_Lvvo += -1.000000 <n,a||f,i>_bbbb l2_bba(i,a,e) 
            // flops: o2v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") += V_blks_["bbbb_vovo"]("a,n,f,i") * l2_xbba_Lovv("I,i,a,e");

            // sigmal2_xabb_Lvvo += -0.500000 <a,b||e,f>_abab l2_bba(n,b,a) 
            // flops: o1v4L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_vvvv"]("a,b,e,f") * l2_xbba_Lovv("I,n,b,a");

            // sigmal2_xabb_Lvvo += 0.500000 <j,i||b,f>_abab l2_bba(n,a,e) t2_abab(b,a,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("j,i,b,f") * l2_xbba_Lovv("I,n,a,e") * t2_abab_vvoo("b,a,j,i");

            // sigmal2_xabb_Lvvo += 0.500000 <i,j||b,f>_abab l2_bba(n,a,e) t2_abab(b,a,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("i,j,b,f") * l2_xbba_Lovv("I,n,a,e") * t2_abab_vvoo("b,a,i,j");

            // sigmal2_xabb_Lvvo += 0.500000 <j,i||b,f>_bbbb l2_bba(n,a,e) t2_bbbb(b,a,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["bbbb_oovv"]("j,i,b,f") * l2_xbba_Lovv("I,n,a,e") * t2_bbbb_vvoo("b,a,j,i");

            // sigmal2_xabb_Lvvo += 1.000000 <j,n||b,f>_abab l2_aaa(i,a,e) t2_aaaa(b,a,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("j,n,b,f") * l2_xaaa_Lovv("I,i,a,e") * t2_aaaa_vvoo("b,a,i,j");

            // sigmal2_xabb_Lvvo += 1.000000 <n,j||b,f>_bbbb l2_aaa(i,a,e) t2_abab(a,b,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") += V_blks_["bbbb_oovv"]("n,j,b,f") * l2_xaaa_Lovv("I,i,a,e") * t2_abab_vvoo("a,b,i,j");

            // sigmal2_xabb_Lvvo += -1.000000 <j,n||b,f>_abab l2_bba(i,a,e) t2_abab(b,a,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= V_blks_["abab_oovv"]("j,n,b,f") * l2_xbba_Lovv("I,i,a,e") * t2_abab_vvoo("b,a,j,i");

            // sigmal2_xabb_Lvvo += -1.000000 <n,j||b,f>_bbbb l2_bba(i,a,e) t2_bbbb(b,a,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= V_blks_["bbbb_oovv"]("n,j,b,f") * l2_xbba_Lovv("I,i,a,e") * t2_bbbb_vvoo("b,a,i,j");

            // sigmal2_xabb_Lvvo += -0.250000 <j,i||e,f>_abab l2_bba(n,b,a) t2_abab(a,b,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= 0.250000 * V_blks_["abab_oovv"]("j,i,e,f") * l2_xbba_Lovv("I,n,b,a") * t2_abab_vvoo("a,b,j,i");

            // sigmal2_xabb_Lvvo += -0.250000 <i,j||e,f>_abab l2_bba(n,b,a) t2_abab(a,b,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= 0.250000 * V_blks_["abab_oovv"]("i,j,e,f") * l2_xbba_Lovv("I,n,b,a") * t2_abab_vvoo("a,b,i,j");

            // sigmal2_xabb_Lvvo += 0.500000 <j,n||e,f>_abab l2_aaa(i,b,a) t2_aaaa(b,a,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("j,n,e,f") * l2_xaaa_Lovv("I,i,b,a") * t2_aaaa_vvoo("b,a,i,j");

            // sigmal2_xabb_Lvvo += 0.500000 <j,n||e,f>_abab l2_bba(i,b,a) t2_abab(a,b,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("j,n,e,f") * l2_xbba_Lovv("I,i,b,a") * t2_abab_vvoo("a,b,j,i");
        }

        if (include_u1_) {

            // sigmal2_xabb_Lvvo += 1.000000 <i,n||a,f>_abab u1_aa(a,i) m1_a(e) 
            // flops: o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("i,n,a,f") * u1_aa_vo("a,i") * m1_xa_Lv("I,e");

            // sigmal2_xabb_Lvvo += -1.000000 <n,i||a,f>_bbbb u1_bb(a,i) m1_a(e) 
            // flops: o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= V_blks_["bbbb_oovv"]("n,i,a,f") * u1_bb_vo("a,i") * m1_xa_Lv("I,e");

            // sigmal2_xabb_Lvvo += -1.000000 <i,n||e,f>_abab u1_aa(a,i) m1_a(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= V_blks_["abab_oovv"]("i,n,e,f") * u1_aa_vo("a,i") * m1_xa_Lv("I,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmal2_xabb_Lvvo += 1.000000 <j,n||i,f>_abab u1_aa(a,j) m2_aaa(i,a,e) 
            // flops: o2v3L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= V_blks_["abba_oovo"]("j,n,f,i") * u1_aa_vo("a,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal2_xabb_Lvvo += 1.000000 <n,j||f,i>_bbbb u1_bb(a,j) m2_bba(i,a,e) 
            // flops: o2v3L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") += V_blks_["bbbb_oovo"]("n,j,f,i") * u1_bb_vo("a,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmal2_xabb_Lvvo += -1.000000 <i,a||b,f>_abab u1_aa(b,i) m2_bba(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1, o1v3: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") += V_blks_["baab_vovv"]("a,i,b,f") * u1_aa_vo("b,i") * m2_xbba_Lovv("I,n,a,e");

            // sigmal2_xabb_Lvvo += -1.000000 <i,a||b,f>_bbbb u1_bb(b,i) m2_bba(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1, o1v3: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") += V_blks_["bbbb_vovv"]("a,i,b,f") * u1_bb_vo("b,i") * m2_xbba_Lovv("I,n,a,e");

            // sigmal2_xabb_Lvvo += -1.000000 <a,n||b,f>_abab u1_aa(b,i) m2_aaa(i,a,e) 
            // flops: o2v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= V_blks_["abab_vovv"]("a,n,b,f") * u1_aa_vo("b,i") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal2_xabb_Lvvo += 1.000000 <n,a||b,f>_bbbb u1_bb(b,i) m2_bba(i,a,e) 
            // flops: o2v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= V_blks_["bbbb_vovv"]("a,n,b,f") * u1_bb_vo("b,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmal2_xabb_Lvvo += 1.000000 <i,b||e,f>_abab u1_aa(a,i) m2_bba(n,b,a) 
            // flops: o1v4L1: 1, o1v4: 1, o1v2L1: 1 | mem: o0v4: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= V_blks_["baab_vovv"]("b,i,e,f") * u1_aa_vo("a,i") * m2_xbba_Lovv("I,n,b,a");
        }

        if (include_u2_) {

            // sigmal2_xabb_Lvvo += 0.500000 <j,i||b,f>_abab u2_abab(b,a,j,i) m2_bba(n,a,e) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("j,i,b,f") * u2_abab_vvoo("b,a,j,i") * m2_xbba_Lovv("I,n,a,e");

            // sigmal2_xabb_Lvvo += 0.500000 <i,j||b,f>_abab u2_abab(b,a,i,j) m2_bba(n,a,e) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("i,j,b,f") * u2_abab_vvoo("b,a,i,j") * m2_xbba_Lovv("I,n,a,e");

            // sigmal2_xabb_Lvvo += 0.500000 <j,i||b,f>_bbbb u2_bbbb(b,a,j,i) m2_bba(n,a,e) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["bbbb_oovv"]("j,i,b,f") * u2_bbbb_vvoo("b,a,j,i") * m2_xbba_Lovv("I,n,a,e");

            // sigmal2_xabb_Lvvo += 1.000000 <j,n||b,f>_abab u2_aaaa(b,a,i,j) m2_aaa(i,a,e) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("j,n,b,f") * u2_aaaa_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal2_xabb_Lvvo += -1.000000 <j,n||b,f>_abab u2_abab(b,a,j,i) m2_bba(i,a,e) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= V_blks_["abab_oovv"]("j,n,b,f") * u2_abab_vvoo("b,a,j,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmal2_xabb_Lvvo += 1.000000 <n,j||b,f>_bbbb u2_abab(a,b,i,j) m2_aaa(i,a,e) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") += V_blks_["bbbb_oovv"]("n,j,b,f") * u2_abab_vvoo("a,b,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmal2_xabb_Lvvo += -1.000000 <n,j||b,f>_bbbb u2_bbbb(b,a,i,j) m2_bba(i,a,e) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= V_blks_["bbbb_oovv"]("n,j,b,f") * u2_bbbb_vvoo("b,a,i,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmal2_xabb_Lvvo += -0.250000 <j,i||e,f>_abab u2_abab(a,b,j,i) m2_bba(n,b,a) 
            // flops: o1v4L1: 1, o2v4: 1, o1v2L1: 1 | mem: o0v4: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= 0.250000 * V_blks_["abab_oovv"]("j,i,e,f") * u2_abab_vvoo("a,b,j,i") * m2_xbba_Lovv("I,n,b,a");

            // sigmal2_xabb_Lvvo += -0.250000 <i,j||e,f>_abab u2_abab(a,b,i,j) m2_bba(n,b,a) 
            // flops: o1v4L1: 1, o2v4: 1, o1v2L1: 1 | mem: o0v4: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") -= 0.250000 * V_blks_["abab_oovv"]("i,j,e,f") * u2_abab_vvoo("a,b,i,j") * m2_xbba_Lovv("I,n,b,a");

            // sigmal2_xabb_Lvvo += 0.500000 <j,n||e,f>_abab u2_aaaa(b,a,i,j) m2_aaa(i,b,a) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("j,n,e,f") * u2_aaaa_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,b,a");

            // sigmal2_xabb_Lvvo += 0.500000 <j,n||e,f>_abab u2_abab(a,b,j,i) m2_bba(i,b,a) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmal2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("j,n,e,f") * u2_abab_vvoo("a,b,j,i") * m2_xbba_Lovv("I,i,b,a");
        }




/// ****** sigmal2_baa(e,f,n) ****** ///



        {

            // sigmal2_xbaa_Lvvo += 1.000000 f_aa(n,f) l1_b(e) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += F_blks_["aa_ov"]("n,f") * l1_xb_Lv("I,e");

            // sigmal2_xbaa_Lvvo += -1.000000 f_aa(a,f) l2_aab(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= F_blks_["aa_vv"]("a,f") * l2_xaab_Lovv("I,n,a,e");
        }

        if (include_u1_ && include_u2_) {

            // sigmal2_xbaa_Lvvo += 1.000000 f_aa(i,f) u1_aa(a,i) m2_aab(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1, o1v2: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += F_blks_["aa_ov"]("i,f") * u1_aa_vo("a,i") * m2_xaab_Lovv("I,n,a,e");
        }

        if (include_u1_) {

            // sigmal2_xbaa_Lvvo += -1.000000 d+_aa(n,f) m1_b(e) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= dp_aa_ov("n,f") * m1_xb_Lv("I,e");
        }

        if (include_u2_) {

            // sigmal2_xbaa_Lvvo += 1.000000 d+_aa(a,f) m2_aab(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += dp_aa_vv("a,f") * m2_xaab_Lovv("I,n,a,e");
        }

        if (include_u0_) {

            // sigmal2_xbaa_Lvvo += -1.000000 d-_aa(n,f) l1_b(e) u0 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, o1v2: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= dp_aa_ov("n,f") * l1_xb_Lv("I,e") * u0;

            // sigmal2_xbaa_Lvvo += 1.000000 d-_aa(a,f) l2_aab(n,a,e) u0 
            // flops: o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v2: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += dp_aa_vv("a,f") * l2_xaab_Lovv("I,n,a,e") * u0;
        }

        if (include_u1_) {

            // sigmal2_xbaa_Lvvo += -1.000000 d-_aa(i,f) l2_aab(n,a,e) u1_aa(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= dp_aa_ov("i,f") * l2_xaab_Lovv("I,n,a,e") * u1_aa_vo("a,i");

            // sigmal2_xbaa_Lvvo += 1.000000 d-_aa(n,f) l2_aab(i,a,e) u1_aa(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += dp_aa_ov("n,f") * l2_xaab_Lovv("I,i,a,e") * u1_aa_vo("a,i");

            // sigmal2_xbaa_Lvvo += 1.000000 d-_aa(n,f) l2_bbb(i,a,e) u1_bb(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += dp_aa_ov("n,f") * l2_xbbb_Lovv("I,i,a,e") * u1_bb_vo("a,i");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmal2_xbaa_Lvvo += -1.000000 d-_aa(i,f) u0 u1_aa(a,i) m2_aab(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1, o1v2: 1, o1v1: 1 | mem: o1v2L1: 2, o0v2: 1, o1v1: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= dp_aa_ov("i,f") * u0 * u1_aa_vo("a,i") * m2_xaab_Lovv("I,n,a,e");
        }

        {

            // sigmal2_xbaa_Lvvo += 1.000000 <n,a||f,e>_abab l1_b(a) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= V_blks_["baab_vovv"]("a,n,f,e") * l1_xb_Lv("I,a");

            // sigmal2_xbaa_Lvvo += -1.000000 <n,a||f,i>_aaaa l2_aab(i,a,e) 
            // flops: o2v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += V_blks_["aaaa_vovo"]("a,n,f,i") * l2_xaab_Lovv("I,i,a,e");

            // sigmal2_xbaa_Lvvo += -1.000000 <n,a||f,i>_abab l2_bbb(i,a,e) 
            // flops: o2v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += V_blks_["baab_vovo"]("a,n,f,i") * l2_xbbb_Lovv("I,i,a,e");

            // sigmal2_xbaa_Lvvo += -0.500000 <b,a||f,e>_abab l2_aab(n,b,a) 
            // flops: o1v4L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_vvvv"]("b,a,f,e") * l2_xaab_Lovv("I,n,b,a");

            // sigmal2_xbaa_Lvvo += 0.500000 <j,i||b,f>_aaaa l2_aab(n,a,e) t2_aaaa(b,a,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["aaaa_oovv"]("j,i,b,f") * l2_xaab_Lovv("I,n,a,e") * t2_aaaa_vvoo("b,a,j,i");

            // sigmal2_xbaa_Lvvo += 0.500000 <j,i||f,b>_abab l2_aab(n,a,e) t2_abab(a,b,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("j,i,f,b") * l2_xaab_Lovv("I,n,a,e") * t2_abab_vvoo("a,b,j,i");

            // sigmal2_xbaa_Lvvo += 0.500000 <i,j||f,b>_abab l2_aab(n,a,e) t2_abab(a,b,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("i,j,f,b") * l2_xaab_Lovv("I,n,a,e") * t2_abab_vvoo("a,b,i,j");

            // sigmal2_xbaa_Lvvo += -1.000000 <n,j||b,f>_aaaa l2_aab(i,a,e) t2_aaaa(b,a,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= V_blks_["aaaa_oovv"]("n,j,b,f") * l2_xaab_Lovv("I,i,a,e") * t2_aaaa_vvoo("b,a,i,j");

            // sigmal2_xbaa_Lvvo += -1.000000 <n,j||f,b>_abab l2_aab(i,a,e) t2_abab(a,b,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= V_blks_["abab_oovv"]("n,j,f,b") * l2_xaab_Lovv("I,i,a,e") * t2_abab_vvoo("a,b,i,j");

            // sigmal2_xbaa_Lvvo += 1.000000 <n,j||b,f>_aaaa l2_bbb(i,a,e) t2_abab(b,a,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += V_blks_["aaaa_oovv"]("n,j,b,f") * l2_xbbb_Lovv("I,i,a,e") * t2_abab_vvoo("b,a,j,i");

            // sigmal2_xbaa_Lvvo += 1.000000 <n,j||f,b>_abab l2_bbb(i,a,e) t2_bbbb(b,a,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("n,j,f,b") * l2_xbbb_Lovv("I,i,a,e") * t2_bbbb_vvoo("b,a,i,j");

            // sigmal2_xbaa_Lvvo += -0.250000 <j,i||f,e>_abab l2_aab(n,b,a) t2_abab(b,a,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= 0.250000 * V_blks_["abab_oovv"]("j,i,f,e") * l2_xaab_Lovv("I,n,b,a") * t2_abab_vvoo("b,a,j,i");

            // sigmal2_xbaa_Lvvo += -0.250000 <i,j||f,e>_abab l2_aab(n,b,a) t2_abab(b,a,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= 0.250000 * V_blks_["abab_oovv"]("i,j,f,e") * l2_xaab_Lovv("I,n,b,a") * t2_abab_vvoo("b,a,i,j");

            // sigmal2_xbaa_Lvvo += 0.500000 <n,j||f,e>_abab l2_aab(i,b,a) t2_abab(b,a,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("n,j,f,e") * l2_xaab_Lovv("I,i,b,a") * t2_abab_vvoo("b,a,i,j");

            // sigmal2_xbaa_Lvvo += 0.500000 <n,j||f,e>_abab l2_bbb(i,b,a) t2_bbbb(b,a,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("n,j,f,e") * l2_xbbb_Lovv("I,i,b,a") * t2_bbbb_vvoo("b,a,i,j");
        }

        if (include_u1_) {

            // sigmal2_xbaa_Lvvo += -1.000000 <n,i||a,f>_aaaa u1_aa(a,i) m1_b(e) 
            // flops: o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= V_blks_["aaaa_oovv"]("n,i,a,f") * u1_aa_vo("a,i") * m1_xb_Lv("I,e");

            // sigmal2_xbaa_Lvvo += 1.000000 <n,i||f,a>_abab u1_bb(a,i) m1_b(e) 
            // flops: o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("n,i,f,a") * u1_bb_vo("a,i") * m1_xb_Lv("I,e");

            // sigmal2_xbaa_Lvvo += -1.000000 <n,i||f,e>_abab u1_bb(a,i) m1_b(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= V_blks_["abab_oovv"]("n,i,f,e") * u1_bb_vo("a,i") * m1_xb_Lv("I,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmal2_xbaa_Lvvo += 1.000000 <n,j||f,i>_aaaa u1_aa(a,j) m2_aab(i,a,e) 
            // flops: o2v3L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += V_blks_["aaaa_oovo"]("n,j,f,i") * u1_aa_vo("a,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal2_xbaa_Lvvo += 1.000000 <n,j||f,i>_abab u1_bb(a,j) m2_bbb(i,a,e) 
            // flops: o2v3L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += V_blks_["abab_oovo"]("n,j,f,i") * u1_bb_vo("a,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal2_xbaa_Lvvo += -1.000000 <i,a||b,f>_aaaa u1_aa(b,i) m2_aab(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1, o1v3: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += V_blks_["aaaa_vovv"]("a,i,b,f") * u1_aa_vo("b,i") * m2_xaab_Lovv("I,n,a,e");

            // sigmal2_xbaa_Lvvo += -1.000000 <a,i||f,b>_abab u1_bb(b,i) m2_aab(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1, o1v3: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= V_blks_["abab_vovv"]("a,i,f,b") * u1_bb_vo("b,i") * m2_xaab_Lovv("I,n,a,e");

            // sigmal2_xbaa_Lvvo += 1.000000 <n,a||b,f>_aaaa u1_aa(b,i) m2_aab(i,a,e) 
            // flops: o2v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= V_blks_["aaaa_vovv"]("a,n,b,f") * u1_aa_vo("b,i") * m2_xaab_Lovv("I,i,a,e");

            // sigmal2_xbaa_Lvvo += -1.000000 <n,a||f,b>_abab u1_bb(b,i) m2_bbb(i,a,e) 
            // flops: o2v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += V_blks_["baab_vovv"]("a,n,f,b") * u1_bb_vo("b,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal2_xbaa_Lvvo += 1.000000 <b,i||f,e>_abab u1_bb(a,i) m2_aab(n,b,a) 
            // flops: o1v4L1: 1, o1v4: 1, o1v2L1: 1 | mem: o0v4: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += V_blks_["abab_vovv"]("b,i,f,e") * u1_bb_vo("a,i") * m2_xaab_Lovv("I,n,b,a");
        }

        if (include_u2_) {

            // sigmal2_xbaa_Lvvo += 0.500000 <j,i||b,f>_aaaa u2_aaaa(b,a,j,i) m2_aab(n,a,e) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["aaaa_oovv"]("j,i,b,f") * u2_aaaa_vvoo("b,a,j,i") * m2_xaab_Lovv("I,n,a,e");

            // sigmal2_xbaa_Lvvo += 0.500000 <j,i||f,b>_abab u2_abab(a,b,j,i) m2_aab(n,a,e) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("j,i,f,b") * u2_abab_vvoo("a,b,j,i") * m2_xaab_Lovv("I,n,a,e");

            // sigmal2_xbaa_Lvvo += 0.500000 <i,j||f,b>_abab u2_abab(a,b,i,j) m2_aab(n,a,e) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("i,j,f,b") * u2_abab_vvoo("a,b,i,j") * m2_xaab_Lovv("I,n,a,e");

            // sigmal2_xbaa_Lvvo += -1.000000 <n,j||b,f>_aaaa u2_aaaa(b,a,i,j) m2_aab(i,a,e) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= V_blks_["aaaa_oovv"]("n,j,b,f") * u2_aaaa_vvoo("b,a,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal2_xbaa_Lvvo += 1.000000 <n,j||b,f>_aaaa u2_abab(b,a,j,i) m2_bbb(i,a,e) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += V_blks_["aaaa_oovv"]("n,j,b,f") * u2_abab_vvoo("b,a,j,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal2_xbaa_Lvvo += -1.000000 <n,j||f,b>_abab u2_abab(a,b,i,j) m2_aab(i,a,e) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= V_blks_["abab_oovv"]("n,j,f,b") * u2_abab_vvoo("a,b,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmal2_xbaa_Lvvo += 1.000000 <n,j||f,b>_abab u2_bbbb(b,a,i,j) m2_bbb(i,a,e) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("n,j,f,b") * u2_bbbb_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmal2_xbaa_Lvvo += -0.250000 <j,i||f,e>_abab u2_abab(b,a,j,i) m2_aab(n,b,a) 
            // flops: o1v4L1: 1, o2v4: 1, o1v2L1: 1 | mem: o0v4: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= 0.250000 * V_blks_["abab_oovv"]("j,i,f,e") * u2_abab_vvoo("b,a,j,i") * m2_xaab_Lovv("I,n,b,a");

            // sigmal2_xbaa_Lvvo += -0.250000 <i,j||f,e>_abab u2_abab(b,a,i,j) m2_aab(n,b,a) 
            // flops: o1v4L1: 1, o2v4: 1, o1v2L1: 1 | mem: o0v4: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") -= 0.250000 * V_blks_["abab_oovv"]("i,j,f,e") * u2_abab_vvoo("b,a,i,j") * m2_xaab_Lovv("I,n,b,a");

            // sigmal2_xbaa_Lvvo += 0.500000 <n,j||f,e>_abab u2_abab(b,a,i,j) m2_aab(i,b,a) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("n,j,f,e") * u2_abab_vvoo("b,a,i,j") * m2_xaab_Lovv("I,i,b,a");

            // sigmal2_xbaa_Lvvo += 0.500000 <n,j||f,e>_abab u2_bbbb(b,a,i,j) m2_bbb(i,b,a) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmal2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("n,j,f,e") * u2_bbbb_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,i,b,a");
        }




/// ****** sigmal2_bbb(e,f,n) ****** ///



        if (include_u0_ && include_u2_) {

            // sigmal2_xbbb_Lvvo += 1.000000 u0 m2_bbb(n,e,f) w0 
            // flops: o1v2L1: 2, o1v2: 1 | mem: o1v2L1: 1, o1v2: 2, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += u0 * m2_xbbb_Lovv("I,n,e,f") * w0;
        }

        {

            // sigmal2_xbbb_Lvvo += -1.000000 P(e,f) f_bb(n,e) l1_b(f) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = F_blks_["bb_ov"]("n,e") * l1_xb_Lv("I,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += 1.000000 f_aa(i,i) l2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += scalar5 * l2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += 1.000000 f_bb(i,i) l2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += scalar6 * l2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += -1.000000 f_bb(n,i) l2_bbb(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= F_blks_["bb_oo"]("n,i") * l2_xbbb_Lovv("I,i,e,f");

            // sigmal2_xbbb_Lvvo += 1.000000 P(e,f) f_bb(a,e) l2_bbb(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = F_blks_["bb_vv"]("a,e") * l2_xbbb_Lovv("I,n,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u1_ && include_u2_) {

            // sigmal2_xbbb_Lvvo += 1.000000 f_aa(i,a) u1_aa(a,i) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += scalar12 * m2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += 1.000000 f_bb(i,a) u1_bb(a,i) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += scalar13 * m2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += -1.000000 f_bb(n,a) u1_bb(a,i) m2_bbb(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1, o2v1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= F_blks_["bb_ov"]("n,a") * u1_bb_vo("a,i") * m2_xbbb_Lovv("I,i,e,f");

            // sigmal2_xbbb_Lvvo += -1.000000 P(e,f) f_bb(i,e) u1_bb(a,i) m2_bbb(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2, o1v2: 1 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = F_blks_["bb_ov"]("i,e") * u1_bb_vo("a,i") * m2_xbbb_Lovv("I,n,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u1_) {

            // sigmal2_xbbb_Lvvo += 1.000000 P(e,f) d+_bb(n,e) m1_b(f) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("n,e") * m1_xb_Lv("I,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmal2_xbbb_Lvvo += -1.000000 d+_aa(i,i) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= scalar17 * m2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += -1.000000 d+_bb(i,i) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= scalar18 * m2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += 1.000000 d+_bb(n,i) m2_bbb(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += dp_bb_oo("n,i") * m2_xbbb_Lovv("I,i,e,f");

            // sigmal2_xbbb_Lvvo += -1.000000 P(e,f) d+_bb(a,e) m2_bbb(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_vv("a,e") * m2_xbbb_Lovv("I,n,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u0_) {

            // sigmal2_xbbb_Lvvo += 1.000000 P(e,f) d-_bb(n,e) l1_b(f) u0 
            // flops: o1v2L1: 4 | mem: o1v2L1: 2, o1v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("n,e") * l1_xb_Lv("I,f") * u0;
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += -1.000000 d-_aa(i,i) l2_bbb(n,e,f) u0 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= scalar7 * u0 * l2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += -1.000000 d-_bb(i,i) l2_bbb(n,e,f) u0 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= scalar8 * u0 * l2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += 1.000000 d-_bb(n,i) l2_bbb(i,e,f) u0 
            // flops: o2v2L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += dp_bb_oo("n,i") * l2_xbbb_Lovv("I,i,e,f") * u0;

            // sigmal2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(a,e) l2_bbb(n,a,f) u0 
            // flops: o1v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_vv("a,e") * l2_xbbb_Lovv("I,n,a,f") * u0;
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u1_) {

            // sigmal2_xbbb_Lvvo += -1.000000 d-_aa(i,a) l2_bbb(n,e,f) u1_aa(a,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= scalar0 * l2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += -1.000000 d-_bb(i,a) l2_bbb(n,e,f) u1_bb(a,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= scalar1 * l2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += 1.000000 d-_bb(n,a) l2_bbb(i,e,f) u1_bb(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += dp_bb_ov("n,a") * l2_xbbb_Lovv("I,i,e,f") * u1_bb_vo("a,i");

            // sigmal2_xbbb_Lvvo += 1.000000 P(e,f) d-_bb(i,e) l2_bbb(n,a,f) u1_bb(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("i,e") * l2_xbbb_Lovv("I,n,a,f") * u1_bb_vo("a,i");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(n,e) l2_aab(i,a,f) u1_aa(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("n,e") * l2_xaab_Lovv("I,i,a,f") * u1_aa_vo("a,i");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(n,e) l2_bbb(i,a,f) u1_bb(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("n,e") * l2_xbbb_Lovv("I,i,a,f") * u1_bb_vo("a,i");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmal2_xbbb_Lvvo += -1.000000 d-_aa(i,a) u0 u1_aa(a,i) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= scalar0 * u0 * m2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += -1.000000 d-_bb(i,a) u0 u1_bb(a,i) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= scalar1 * u0 * m2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += 1.000000 d-_bb(n,a) u0 u1_bb(a,i) m2_bbb(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1, o2v1: 1, o1v1: 1 | mem: o1v2L1: 2, o1v1: 1, o2v0: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += dp_bb_ov("n,a") * u0 * u1_bb_vo("a,i") * m2_xbbb_Lovv("I,i,e,f");

            // sigmal2_xbbb_Lvvo += 1.000000 P(e,f) d-_bb(i,e) u0 u1_bb(a,i) m2_bbb(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2, o1v2: 1, o1v1: 1 | mem: o1v2L1: 2, o0v2: 1, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("i,e") * u0 * u1_bb_vo("a,i") * m2_xbbb_Lovv("I,n,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        {

            // sigmal2_xbbb_Lvvo += -0.500000 <j,i||j,i>_aaaa l2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * scalar9 * l2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += -0.500000 <j,i||j,i>_abab l2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * scalar10 * l2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += -0.500000 <i,j||i,j>_abab l2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * scalar10 * l2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += -0.500000 <j,i||j,i>_bbbb l2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * scalar11 * l2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += -1.000000 <n,a||e,f>_bbbb l1_b(a) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += V_blks_["bbbb_vovv"]("a,n,e,f") * l1_xb_Lv("I,a");

            // sigmal2_xbbb_Lvvo += 1.000000 P(e,f) <a,n||i,e>_abab l2_aab(i,a,f) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["abba_vovo"]("a,n,e,i") * l2_xaab_Lovv("I,i,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += 1.000000 P(e,f) <n,a||e,i>_bbbb l2_bbb(i,a,f) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_vovo"]("a,n,e,i") * l2_xbbb_Lovv("I,i,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += 0.500000 <b,a||e,f>_bbbb l2_bbb(n,b,a) 
            // flops: o1v4L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["bbbb_vvvv"]("b,a,e,f") * l2_xbbb_Lovv("I,n,b,a");

            // sigmal2_xbbb_Lvvo += 0.250000 <j,i||a,b>_aaaa l2_bbb(n,e,f) t2_aaaa(a,b,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar2 * l2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += 0.250000 <j,i||a,b>_abab l2_bbb(n,e,f) t2_abab(a,b,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar3 * l2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += 0.250000 <i,j||a,b>_abab l2_bbb(n,e,f) t2_abab(a,b,i,j) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar3 * l2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += 0.250000 <j,i||b,a>_abab l2_bbb(n,e,f) t2_abab(b,a,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar3 * l2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += 0.250000 <i,j||b,a>_abab l2_bbb(n,e,f) t2_abab(b,a,i,j) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar3 * l2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += 0.250000 <j,i||a,b>_bbbb l2_bbb(n,e,f) t2_bbbb(a,b,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar4 * l2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += -0.500000 <j,n||a,b>_abab l2_bbb(i,e,f) t2_abab(a,b,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("j,n,a,b") * l2_xbbb_Lovv("I,i,e,f") * t2_abab_vvoo("a,b,j,i");

            // sigmal2_xbbb_Lvvo += -0.500000 <j,n||b,a>_abab l2_bbb(i,e,f) t2_abab(b,a,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("j,n,b,a") * l2_xbbb_Lovv("I,i,e,f") * t2_abab_vvoo("b,a,j,i");

            // sigmal2_xbbb_Lvvo += -0.500000 <n,j||a,b>_bbbb l2_bbb(i,e,f) t2_bbbb(a,b,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovv"]("n,j,a,b") * l2_xbbb_Lovv("I,i,e,f") * t2_bbbb_vvoo("a,b,i,j");

            // sigmal2_xbbb_Lvvo += -0.500000 P(e,f) <j,i||b,e>_abab l2_bbb(n,a,f) t2_abab(b,a,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 2 | mem: o3v4L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("j,i,b,e") * l2_xbbb_Lovv("I,n,a,f") * t2_abab_vvoo("b,a,j,i");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += -0.500000 P(e,f) <i,j||b,e>_abab l2_bbb(n,a,f) t2_abab(b,a,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 2 | mem: o3v4L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("i,j,b,e") * l2_xbbb_Lovv("I,n,a,f") * t2_abab_vvoo("b,a,i,j");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += -0.500000 P(e,f) <j,i||b,e>_bbbb l2_bbb(n,a,f) t2_bbbb(b,a,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 2 | mem: o3v4L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["bbbb_oovv"]("j,i,b,e") * l2_xbbb_Lovv("I,n,a,f") * t2_bbbb_vvoo("b,a,j,i");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += -1.000000 P(e,f) <j,n||b,e>_abab l2_aab(i,a,f) t2_aaaa(b,a,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 2 | mem: o3v4L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("j,n,b,e") * l2_xaab_Lovv("I,i,a,f") * t2_aaaa_vvoo("b,a,i,j");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += -1.000000 P(e,f) <n,j||b,e>_bbbb l2_aab(i,a,f) t2_abab(a,b,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 2 | mem: o3v4L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_oovv"]("n,j,b,e") * l2_xaab_Lovv("I,i,a,f") * t2_abab_vvoo("a,b,i,j");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += 1.000000 P(e,f) <j,n||b,e>_abab l2_bbb(i,a,f) t2_abab(b,a,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 2 | mem: o3v4L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("j,n,b,e") * l2_xbbb_Lovv("I,i,a,f") * t2_abab_vvoo("b,a,j,i");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += 1.000000 P(e,f) <n,j||b,e>_bbbb l2_bbb(i,a,f) t2_bbbb(b,a,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 2 | mem: o3v4L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_oovv"]("n,j,b,e") * l2_xbbb_Lovv("I,i,a,f") * t2_bbbb_vvoo("b,a,i,j");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += 0.250000 <j,i||e,f>_bbbb l2_bbb(n,b,a) t2_bbbb(b,a,j,i) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += 0.250000 * V_blks_["bbbb_oovv"]("j,i,e,f") * l2_xbbb_Lovv("I,n,b,a") * t2_bbbb_vvoo("b,a,j,i");

            // sigmal2_xbbb_Lvvo += -0.500000 <n,j||e,f>_bbbb l2_aab(i,b,a) t2_abab(b,a,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovv"]("n,j,e,f") * l2_xaab_Lovv("I,i,b,a") * t2_abab_vvoo("b,a,i,j");

            // sigmal2_xbbb_Lvvo += -0.500000 <n,j||e,f>_bbbb l2_bbb(i,b,a) t2_bbbb(b,a,i,j) 
            // flops: o3v4L1: 2, o1v2L1: 1 | mem: o3v4L1: 1, o1v2L1: 2, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovv"]("n,j,e,f") * l2_xbbb_Lovv("I,i,b,a") * t2_bbbb_vvoo("b,a,i,j");
        }

        if (include_u1_) {

            // sigmal2_xbbb_Lvvo += -1.000000 P(e,f) <i,n||a,e>_abab u1_aa(a,i) m1_b(f) 
            // flops: o1v2L1: 3, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("i,n,a,e") * u1_aa_vo("a,i") * m1_xb_Lv("I,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += 1.000000 P(e,f) <n,i||a,e>_bbbb u1_bb(a,i) m1_b(f) 
            // flops: o1v2L1: 3, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_oovv"]("n,i,a,e") * u1_bb_vo("a,i") * m1_xb_Lv("I,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += 1.000000 <n,i||e,f>_bbbb u1_bb(a,i) m1_b(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += V_blks_["bbbb_oovv"]("n,i,e,f") * u1_bb_vo("a,i") * m1_xb_Lv("I,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmal2_xbbb_Lvvo += -1.000000 <j,n||a,i>_abab u1_aa(a,j) m2_bbb(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1, o3v1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= V_blks_["abab_oovo"]("j,n,a,i") * u1_aa_vo("a,j") * m2_xbbb_Lovv("I,i,e,f");

            // sigmal2_xbbb_Lvvo += 1.000000 <n,j||a,i>_bbbb u1_bb(a,j) m2_bbb(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1, o3v1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += V_blks_["bbbb_oovo"]("n,j,a,i") * u1_bb_vo("a,j") * m2_xbbb_Lovv("I,i,e,f");

            // sigmal2_xbbb_Lvvo += -1.000000 P(e,f) <j,n||i,e>_abab u1_aa(a,j) m2_aab(i,a,f) 
            // flops: o2v3L1: 1, o3v2: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["abba_oovo"]("j,n,e,i") * u1_aa_vo("a,j") * m2_xaab_Lovv("I,i,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += -1.000000 P(e,f) <n,j||e,i>_bbbb u1_bb(a,j) m2_bbb(i,a,f) 
            // flops: o2v3L1: 1, o3v2: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_oovo"]("n,j,e,i") * u1_bb_vo("a,j") * m2_xbbb_Lovv("I,i,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += 1.000000 P(e,f) <i,a||b,e>_abab u1_aa(b,i) m2_bbb(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2, o1v3: 1 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["baab_vovv"]("a,i,b,e") * u1_aa_vo("b,i") * m2_xbbb_Lovv("I,n,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += 1.000000 P(e,f) <i,a||b,e>_bbbb u1_bb(b,i) m2_bbb(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2, o1v3: 1 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_vovv"]("a,i,b,e") * u1_bb_vo("b,i") * m2_xbbb_Lovv("I,n,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += 1.000000 P(e,f) <a,n||b,e>_abab u1_aa(b,i) m2_aab(i,a,f) 
            // flops: o2v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["abab_vovv"]("a,n,b,e") * u1_aa_vo("b,i") * m2_xaab_Lovv("I,i,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += -1.000000 P(e,f) <n,a||b,e>_bbbb u1_bb(b,i) m2_bbb(i,a,f) 
            // flops: o2v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_vovv"]("a,n,b,e") * u1_bb_vo("b,i") * m2_xbbb_Lovv("I,i,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += 1.000000 <i,b||e,f>_bbbb u1_bb(a,i) m2_bbb(n,b,a) 
            // flops: o1v4L1: 1, o1v4: 1, o1v2L1: 1 | mem: o0v4: 1, o1v2L1: 2, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= V_blks_["bbbb_vovv"]("b,i,e,f") * u1_bb_vo("a,i") * m2_xbbb_Lovv("I,n,b,a");
        }

        if (include_u2_) {

            // sigmal2_xbbb_Lvvo += 0.250000 <j,i||a,b>_aaaa u2_aaaa(a,b,j,i) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar14 * m2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += 0.250000 <j,i||a,b>_abab u2_abab(a,b,j,i) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar15 * m2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += 0.250000 <i,j||a,b>_abab u2_abab(a,b,i,j) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar15 * m2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += 0.250000 <j,i||b,a>_abab u2_abab(b,a,j,i) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar15 * m2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += 0.250000 <i,j||b,a>_abab u2_abab(b,a,i,j) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar15 * m2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += 0.250000 <j,i||a,b>_bbbb u2_bbbb(a,b,j,i) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar16 * m2_xbbb_Lovv("I,n,e,f");

            // sigmal2_xbbb_Lvvo += -0.500000 <j,n||a,b>_abab u2_abab(a,b,j,i) m2_bbb(i,e,f) 
            // flops: o2v2L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("j,n,a,b") * u2_abab_vvoo("a,b,j,i") * m2_xbbb_Lovv("I,i,e,f");

            // sigmal2_xbbb_Lvvo += -0.500000 <j,n||b,a>_abab u2_abab(b,a,j,i) m2_bbb(i,e,f) 
            // flops: o2v2L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("j,n,b,a") * u2_abab_vvoo("b,a,j,i") * m2_xbbb_Lovv("I,i,e,f");

            // sigmal2_xbbb_Lvvo += -0.500000 <n,j||a,b>_bbbb u2_bbbb(a,b,i,j) m2_bbb(i,e,f) 
            // flops: o2v2L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovv"]("n,j,a,b") * u2_bbbb_vvoo("a,b,i,j") * m2_xbbb_Lovv("I,i,e,f");

            // sigmal2_xbbb_Lvvo += -0.500000 P(e,f) <j,i||b,e>_abab u2_abab(b,a,j,i) m2_bbb(n,a,f) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("j,i,b,e") * u2_abab_vvoo("b,a,j,i") * m2_xbbb_Lovv("I,n,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += -0.500000 P(e,f) <i,j||b,e>_abab u2_abab(b,a,i,j) m2_bbb(n,a,f) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("i,j,b,e") * u2_abab_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,n,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += -0.500000 P(e,f) <j,i||b,e>_bbbb u2_bbbb(b,a,j,i) m2_bbb(n,a,f) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["bbbb_oovv"]("j,i,b,e") * u2_bbbb_vvoo("b,a,j,i") * m2_xbbb_Lovv("I,n,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += -1.000000 P(e,f) <j,n||b,e>_abab u2_aaaa(b,a,i,j) m2_aab(i,a,f) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("j,n,b,e") * u2_aaaa_vvoo("b,a,i,j") * m2_xaab_Lovv("I,i,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += 1.000000 P(e,f) <j,n||b,e>_abab u2_abab(b,a,j,i) m2_bbb(i,a,f) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("j,n,b,e") * u2_abab_vvoo("b,a,j,i") * m2_xbbb_Lovv("I,i,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += -1.000000 P(e,f) <n,j||b,e>_bbbb u2_abab(a,b,i,j) m2_aab(i,a,f) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_oovv"]("n,j,b,e") * u2_abab_vvoo("a,b,i,j") * m2_xaab_Lovv("I,i,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += 1.000000 P(e,f) <n,j||b,e>_bbbb u2_bbbb(b,a,i,j) m2_bbb(i,a,f) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_oovv"]("n,j,b,e") * u2_bbbb_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,i,a,f");
            sigmal2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmal2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmal2_xbbb_Lvvo += 0.250000 <j,i||e,f>_bbbb u2_bbbb(b,a,j,i) m2_bbb(n,b,a) 
            // flops: o1v4L1: 1, o2v4: 1, o1v2L1: 1 | mem: o0v4: 1, o1v2L1: 2, 
            sigmal2_xbbb_Lvvo("I,e,f,n") += 0.250000 * V_blks_["bbbb_oovv"]("j,i,e,f") * u2_bbbb_vvoo("b,a,j,i") * m2_xbbb_Lovv("I,n,b,a");

            // sigmal2_xbbb_Lvvo += -0.500000 <n,j||e,f>_bbbb u2_abab(b,a,i,j) m2_aab(i,b,a) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovv"]("n,j,e,f") * u2_abab_vvoo("b,a,i,j") * m2_xaab_Lovv("I,i,b,a");

            // sigmal2_xbbb_Lvvo += -0.500000 <n,j||e,f>_bbbb u2_bbbb(b,a,i,j) m2_bbb(i,b,a) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmal2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovv"]("n,j,e,f") * u2_bbbb_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,i,b,a");
        }




/// ****** sigmas1_a(e) ****** ///



        if (include_u1_) {

            // sigmas1_xa_Lv += 1.000000 s1_a(e) w0 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") += s1_xa_Lv("I,e") * w0;
        }

        if (include_u0_ && include_u1_) {

            // sigmas1_xa_Lv += 1.000000 r1_a(e) u0 w0 
            // flops: o0v1L1: 2, o0v1: 1 | mem: o0v1L1: 1, o0v1: 2, 
            sigmas1_xa_Lv("I,e") += r1_xa_Lv("I,e") * u0 * w0;
        }

        if (include_u1_) {

            // sigmas1_xa_Lv += 1.000000 f_aa(i,i) s1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") += scalar5 * s1_xa_Lv("I,e");

            // sigmas1_xa_Lv += 1.000000 f_bb(i,i) s1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") += scalar6 * s1_xa_Lv("I,e");

            // sigmas1_xa_Lv += 1.000000 f_aa(e,a) s1_a(a) 
            // flops: o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmas1_xa_Lv("I,e") += F_blks_["aa_vv"]("e,a") * s1_xa_Lv("I,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmas1_xa_Lv += -1.000000 f_aa(i,a) s2_aaa(a,e,i) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmas1_xa_Lv("I,e") -= F_blks_["aa_ov"]("i,a") * s2_xaaa_Lvvo("I,a,e,i");
        }

        if (include_u1_) {

            // sigmas1_xa_Lv += 1.000000 f_aa(i,a) r1_a(e) u1_aa(a,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") += scalar12 * r1_xa_Lv("I,e");

            // sigmas1_xa_Lv += 1.000000 f_bb(i,a) r1_a(e) u1_bb(a,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") += scalar13 * r1_xa_Lv("I,e");

            // sigmas1_xa_Lv += -1.000000 f_aa(i,a) r1_a(a) u1_aa(e,i) 
            // flops: o1v1L1: 2, o0v1L1: 1 | mem: o0v1L1: 2, o1v0L1: 1, 
            sigmas1_xa_Lv("I,e") -= F_blks_["aa_ov"]("i,a") * r1_xa_Lv("I,a") * u1_aa_vo("e,i");

            // sigmas1_xa_Lv += -1.000000 d+_aa(i,i) r1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") -= scalar17 * r1_xa_Lv("I,e");

            // sigmas1_xa_Lv += -1.000000 d+_bb(i,i) r1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") -= scalar18 * r1_xa_Lv("I,e");

            // sigmas1_xa_Lv += -1.000000 d+_aa(e,a) r1_a(a) 
            // flops: o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmas1_xa_Lv("I,e") -= dp_aa_vv("e,a") * r1_xa_Lv("I,a");

            // sigmas1_xa_Lv += 1.000000 d+_aa(i,a) r2_aaa(a,e,i) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmas1_xa_Lv("I,e") += dp_aa_ov("i,a") * r2_xaaa_Lvvo("I,a,e,i");
        }

        if (include_u0_ && include_u1_) {

            // sigmas1_xa_Lv += -1.000000 d-_aa(i,i) u0 s1_a(e) 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmas1_xa_Lv("I,e") -= u0 * scalar7 * s1_xa_Lv("I,e");

            // sigmas1_xa_Lv += -1.000000 d-_bb(i,i) u0 s1_a(e) 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmas1_xa_Lv("I,e") -= u0 * scalar8 * s1_xa_Lv("I,e");

            // sigmas1_xa_Lv += -1.000000 d-_aa(e,a) u0 s1_a(a) 
            // flops: o0v2L1: 1, o0v1L1: 1, o0v2: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmas1_xa_Lv("I,e") -= dp_aa_vv("e,a") * u0 * s1_xa_Lv("I,a");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmas1_xa_Lv += 1.000000 d-_aa(i,a) u0 s2_aaa(a,e,i) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmas1_xa_Lv("I,e") += dp_aa_ov("i,a") * u0 * s2_xaaa_Lvvo("I,a,e,i");
        }

        if (include_u1_) {

            // sigmas1_xa_Lv += -2.000000 d-_aa(i,a) u1_aa(a,i) s1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") -= 2.000000 * scalar0 * s1_xa_Lv("I,e");

            // sigmas1_xa_Lv += -2.000000 d-_bb(i,a) u1_bb(a,i) s1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") -= 2.000000 * scalar1 * s1_xa_Lv("I,e");

            // sigmas1_xa_Lv += 2.000000 d-_aa(i,a) u1_aa(e,i) s1_a(a) 
            // flops: o0v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmas1_xa_Lv("I,e") += 2.000000 * dp_aa_ov("i,a") * u1_aa_vo("e,i") * s1_xa_Lv("I,a");
        }

        if (include_u0_ && include_u1_) {

            // sigmas1_xa_Lv += -1.000000 d-_aa(i,a) r1_a(e) u0 u1_aa(a,i) 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmas1_xa_Lv("I,e") -= u0 * scalar0 * r1_xa_Lv("I,e");

            // sigmas1_xa_Lv += -1.000000 d-_bb(i,a) r1_a(e) u0 u1_bb(a,i) 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmas1_xa_Lv("I,e") -= u0 * scalar1 * r1_xa_Lv("I,e");

            // sigmas1_xa_Lv += 1.000000 d-_aa(i,a) r1_a(a) u0 u1_aa(e,i) 
            // flops: o1v1L1: 1, o0v1L1: 1, o1v0L1: 1, o1v1: 1 | mem: o0v1L1: 1, o1v0L1: 1, o0v1: 1, o1v0: 1, 
            sigmas1_xa_Lv("I,e") += dp_aa_ov("i,a") * r1_xa_Lv("I,a") * u0 * u1_aa_vo("e,i");
        }

        if (include_u1_) {

            // sigmas1_xa_Lv += -0.500000 <j,i||j,i>_aaaa s1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") -= 0.500000 * scalar9 * s1_xa_Lv("I,e");

            // sigmas1_xa_Lv += -0.500000 <j,i||j,i>_abab s1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") -= 0.500000 * scalar10 * s1_xa_Lv("I,e");

            // sigmas1_xa_Lv += -0.500000 <i,j||i,j>_abab s1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") -= 0.500000 * scalar10 * s1_xa_Lv("I,e");

            // sigmas1_xa_Lv += -0.500000 <j,i||j,i>_bbbb s1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") -= 0.500000 * scalar11 * s1_xa_Lv("I,e");
        }

        if (include_u1_ && include_u2_) {

            // sigmas1_xa_Lv += -0.500000 <i,e||a,b>_aaaa s2_aaa(a,b,i) 
            // flops: o1v3L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmas1_xa_Lv("I,e") += 0.500000 * V_blks_["aaaa_vovv"]("e,i,a,b") * s2_xaaa_Lvvo("I,a,b,i");
        }

        if (include_u1_) {

            // sigmas1_xa_Lv += 0.250000 <j,i||a,b>_aaaa t2_aaaa(a,b,j,i) s1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") += 0.250000 * scalar2 * s1_xa_Lv("I,e");

            // sigmas1_xa_Lv += 0.250000 <j,i||a,b>_abab t2_abab(a,b,j,i) s1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") += 0.250000 * scalar3 * s1_xa_Lv("I,e");

            // sigmas1_xa_Lv += 0.250000 <i,j||a,b>_abab t2_abab(a,b,i,j) s1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") += 0.250000 * scalar3 * s1_xa_Lv("I,e");

            // sigmas1_xa_Lv += 0.250000 <j,i||b,a>_abab t2_abab(b,a,j,i) s1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") += 0.250000 * scalar3 * s1_xa_Lv("I,e");

            // sigmas1_xa_Lv += 0.250000 <i,j||b,a>_abab t2_abab(b,a,i,j) s1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") += 0.250000 * scalar3 * s1_xa_Lv("I,e");

            // sigmas1_xa_Lv += 0.250000 <j,i||a,b>_bbbb t2_bbbb(a,b,j,i) s1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") += 0.250000 * scalar4 * s1_xa_Lv("I,e");

            // sigmas1_xa_Lv += -0.500000 <j,i||a,b>_aaaa t2_aaaa(a,e,j,i) s1_a(b) 
            // flops: o2v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmas1_xa_Lv("I,e") -= 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * t2_aaaa_vvoo("a,e,j,i") * s1_xa_Lv("I,b");

            // sigmas1_xa_Lv += -0.500000 <j,i||b,a>_abab t2_abab(e,a,j,i) s1_a(b) 
            // flops: o2v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmas1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("j,i,b,a") * t2_abab_vvoo("e,a,j,i") * s1_xa_Lv("I,b");

            // sigmas1_xa_Lv += -0.500000 <i,j||b,a>_abab t2_abab(e,a,i,j) s1_a(b) 
            // flops: o2v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmas1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("i,j,b,a") * t2_abab_vvoo("e,a,i,j") * s1_xa_Lv("I,b");

            // sigmas1_xa_Lv += 1.000000 <i,e||a,b>_aaaa r1_a(b) u1_aa(a,i) 
            // flops: o1v3L1: 1, o1v2L1: 1, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmas1_xa_Lv("I,e") -= V_blks_["aaaa_vovv"]("e,i,a,b") * r1_xa_Lv("I,b") * u1_aa_vo("a,i");

            // sigmas1_xa_Lv += 1.000000 <e,i||b,a>_abab r1_a(b) u1_bb(a,i) 
            // flops: o1v3L1: 1, o1v2L1: 1, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmas1_xa_Lv("I,e") += V_blks_["abab_vovv"]("e,i,b,a") * r1_xa_Lv("I,b") * u1_bb_vo("a,i");

            // sigmas1_xa_Lv += 1.000000 <j,i||a,b>_aaaa r2_aaa(b,e,j) u1_aa(a,i) 
            // flops: o2v3L1: 1, o1v2L1: 1, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmas1_xa_Lv("I,e") += V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaa_Lvvo("I,b,e,j") * u1_aa_vo("a,i");

            // sigmas1_xa_Lv += -1.000000 <j,i||b,a>_abab r2_aaa(b,e,j) u1_bb(a,i) 
            // flops: o2v3L1: 1, o1v2L1: 1, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmas1_xa_Lv("I,e") -= V_blks_["abab_oovv"]("j,i,b,a") * r2_xaaa_Lvvo("I,b,e,j") * u1_bb_vo("a,i");

            // sigmas1_xa_Lv += 0.500000 <j,i||a,b>_aaaa r2_aaa(a,b,j) u1_aa(e,i) 
            // flops: o2v2L1: 1, o1v1L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v0L1: 1, 
            sigmas1_xa_Lv("I,e") += 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaa_Lvvo("I,a,b,j") * u1_aa_vo("e,i");
        }

        if (include_u1_ && include_u2_) {

            // sigmas1_xa_Lv += 0.250000 <j,i||a,b>_aaaa r1_a(e) u2_aaaa(a,b,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") += 0.250000 * scalar14 * r1_xa_Lv("I,e");

            // sigmas1_xa_Lv += 0.250000 <j,i||a,b>_abab r1_a(e) u2_abab(a,b,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") += 0.250000 * scalar15 * r1_xa_Lv("I,e");

            // sigmas1_xa_Lv += 0.250000 <i,j||a,b>_abab r1_a(e) u2_abab(a,b,i,j) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") += 0.250000 * scalar15 * r1_xa_Lv("I,e");

            // sigmas1_xa_Lv += 0.250000 <j,i||b,a>_abab r1_a(e) u2_abab(b,a,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") += 0.250000 * scalar15 * r1_xa_Lv("I,e");

            // sigmas1_xa_Lv += 0.250000 <i,j||b,a>_abab r1_a(e) u2_abab(b,a,i,j) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") += 0.250000 * scalar15 * r1_xa_Lv("I,e");

            // sigmas1_xa_Lv += 0.250000 <j,i||a,b>_bbbb r1_a(e) u2_bbbb(a,b,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xa_Lv("I,e") += 0.250000 * scalar16 * r1_xa_Lv("I,e");

            // sigmas1_xa_Lv += -0.500000 <j,i||a,b>_aaaa r1_a(b) u2_aaaa(a,e,j,i) 
            // flops: o2v2L1: 2, o0v1L1: 1 | mem: o2v1L1: 1, o0v1L1: 2, 
            sigmas1_xa_Lv("I,e") -= 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r1_xa_Lv("I,b") * u2_aaaa_vvoo("a,e,j,i");

            // sigmas1_xa_Lv += -0.500000 <j,i||b,a>_abab r1_a(b) u2_abab(e,a,j,i) 
            // flops: o2v2L1: 2, o0v1L1: 1 | mem: o2v1L1: 1, o0v1L1: 2, 
            sigmas1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("j,i,b,a") * r1_xa_Lv("I,b") * u2_abab_vvoo("e,a,j,i");

            // sigmas1_xa_Lv += -0.500000 <i,j||b,a>_abab r1_a(b) u2_abab(e,a,i,j) 
            // flops: o2v2L1: 2, o0v1L1: 1 | mem: o2v1L1: 1, o0v1L1: 2, 
            sigmas1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("i,j,b,a") * r1_xa_Lv("I,b") * u2_abab_vvoo("e,a,i,j");
        }




/// ****** sigmas1_b(e) ****** ///



        if (include_u1_) {

            // sigmas1_xb_Lv += 1.000000 s1_b(e) w0 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") += s1_xb_Lv("I,e") * w0;
        }

        if (include_u0_ && include_u1_) {

            // sigmas1_xb_Lv += 1.000000 r1_b(e) u0 w0 
            // flops: o0v1L1: 2, o0v1: 1 | mem: o0v1L1: 1, o0v1: 2, 
            sigmas1_xb_Lv("I,e") += r1_xb_Lv("I,e") * u0 * w0;
        }

        if (include_u1_) {

            // sigmas1_xb_Lv += 1.000000 f_aa(i,i) s1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") += scalar5 * s1_xb_Lv("I,e");

            // sigmas1_xb_Lv += 1.000000 f_bb(i,i) s1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") += scalar6 * s1_xb_Lv("I,e");

            // sigmas1_xb_Lv += 1.000000 f_bb(e,a) s1_b(a) 
            // flops: o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmas1_xb_Lv("I,e") += F_blks_["bb_vv"]("e,a") * s1_xb_Lv("I,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmas1_xb_Lv += -1.000000 f_bb(i,a) s2_bbb(a,e,i) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmas1_xb_Lv("I,e") -= F_blks_["bb_ov"]("i,a") * s2_xbbb_Lvvo("I,a,e,i");
        }

        if (include_u1_) {

            // sigmas1_xb_Lv += 1.000000 f_aa(i,a) r1_b(e) u1_aa(a,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") += scalar12 * r1_xb_Lv("I,e");

            // sigmas1_xb_Lv += 1.000000 f_bb(i,a) r1_b(e) u1_bb(a,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") += scalar13 * r1_xb_Lv("I,e");

            // sigmas1_xb_Lv += -1.000000 f_bb(i,a) r1_b(a) u1_bb(e,i) 
            // flops: o1v1L1: 2, o0v1L1: 1 | mem: o0v1L1: 2, o1v0L1: 1, 
            sigmas1_xb_Lv("I,e") -= F_blks_["bb_ov"]("i,a") * r1_xb_Lv("I,a") * u1_bb_vo("e,i");

            // sigmas1_xb_Lv += -1.000000 d+_aa(i,i) r1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") -= scalar17 * r1_xb_Lv("I,e");

            // sigmas1_xb_Lv += -1.000000 d+_bb(i,i) r1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") -= scalar18 * r1_xb_Lv("I,e");

            // sigmas1_xb_Lv += -1.000000 d+_bb(e,a) r1_b(a) 
            // flops: o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmas1_xb_Lv("I,e") -= dp_bb_vv("e,a") * r1_xb_Lv("I,a");

            // sigmas1_xb_Lv += 1.000000 d+_bb(i,a) r2_bbb(a,e,i) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmas1_xb_Lv("I,e") += dp_bb_ov("i,a") * r2_xbbb_Lvvo("I,a,e,i");
        }

        if (include_u0_ && include_u1_) {

            // sigmas1_xb_Lv += -1.000000 d-_aa(i,i) u0 s1_b(e) 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmas1_xb_Lv("I,e") -= u0 * scalar7 * s1_xb_Lv("I,e");

            // sigmas1_xb_Lv += -1.000000 d-_bb(i,i) u0 s1_b(e) 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmas1_xb_Lv("I,e") -= u0 * scalar8 * s1_xb_Lv("I,e");

            // sigmas1_xb_Lv += -1.000000 d-_bb(e,a) u0 s1_b(a) 
            // flops: o0v2L1: 1, o0v1L1: 1, o0v2: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmas1_xb_Lv("I,e") -= dp_bb_vv("e,a") * u0 * s1_xb_Lv("I,a");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmas1_xb_Lv += 1.000000 d-_bb(i,a) u0 s2_bbb(a,e,i) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmas1_xb_Lv("I,e") += dp_bb_ov("i,a") * u0 * s2_xbbb_Lvvo("I,a,e,i");
        }

        if (include_u1_) {

            // sigmas1_xb_Lv += -2.000000 d-_aa(i,a) u1_aa(a,i) s1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") -= 2.000000 * scalar0 * s1_xb_Lv("I,e");

            // sigmas1_xb_Lv += -2.000000 d-_bb(i,a) u1_bb(a,i) s1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") -= 2.000000 * scalar1 * s1_xb_Lv("I,e");

            // sigmas1_xb_Lv += 2.000000 d-_bb(i,a) u1_bb(e,i) s1_b(a) 
            // flops: o0v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmas1_xb_Lv("I,e") += 2.000000 * dp_bb_ov("i,a") * u1_bb_vo("e,i") * s1_xb_Lv("I,a");
        }

        if (include_u0_ && include_u1_) {

            // sigmas1_xb_Lv += -1.000000 d-_aa(i,a) r1_b(e) u0 u1_aa(a,i) 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmas1_xb_Lv("I,e") -= u0 * scalar0 * r1_xb_Lv("I,e");

            // sigmas1_xb_Lv += -1.000000 d-_bb(i,a) r1_b(e) u0 u1_bb(a,i) 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmas1_xb_Lv("I,e") -= u0 * scalar1 * r1_xb_Lv("I,e");

            // sigmas1_xb_Lv += 1.000000 d-_bb(i,a) r1_b(a) u0 u1_bb(e,i) 
            // flops: o1v1L1: 1, o0v1L1: 1, o1v0L1: 1, o1v1: 1 | mem: o0v1L1: 1, o1v0L1: 1, o0v1: 1, o1v0: 1, 
            sigmas1_xb_Lv("I,e") += dp_bb_ov("i,a") * r1_xb_Lv("I,a") * u0 * u1_bb_vo("e,i");
        }

        if (include_u1_) {

            // sigmas1_xb_Lv += -0.500000 <j,i||j,i>_aaaa s1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") -= 0.500000 * scalar9 * s1_xb_Lv("I,e");

            // sigmas1_xb_Lv += -0.500000 <j,i||j,i>_abab s1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") -= 0.500000 * scalar10 * s1_xb_Lv("I,e");

            // sigmas1_xb_Lv += -0.500000 <i,j||i,j>_abab s1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") -= 0.500000 * scalar10 * s1_xb_Lv("I,e");

            // sigmas1_xb_Lv += -0.500000 <j,i||j,i>_bbbb s1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") -= 0.500000 * scalar11 * s1_xb_Lv("I,e");
        }

        if (include_u1_ && include_u2_) {

            // sigmas1_xb_Lv += -0.500000 <i,e||a,b>_bbbb s2_bbb(a,b,i) 
            // flops: o1v3L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmas1_xb_Lv("I,e") += 0.500000 * V_blks_["bbbb_vovv"]("e,i,a,b") * s2_xbbb_Lvvo("I,a,b,i");
        }

        if (include_u1_) {

            // sigmas1_xb_Lv += 0.250000 <j,i||a,b>_aaaa t2_aaaa(a,b,j,i) s1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") += 0.250000 * scalar2 * s1_xb_Lv("I,e");

            // sigmas1_xb_Lv += 0.250000 <j,i||a,b>_abab t2_abab(a,b,j,i) s1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") += 0.250000 * scalar3 * s1_xb_Lv("I,e");

            // sigmas1_xb_Lv += 0.250000 <i,j||a,b>_abab t2_abab(a,b,i,j) s1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") += 0.250000 * scalar3 * s1_xb_Lv("I,e");

            // sigmas1_xb_Lv += 0.250000 <j,i||b,a>_abab t2_abab(b,a,j,i) s1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") += 0.250000 * scalar3 * s1_xb_Lv("I,e");

            // sigmas1_xb_Lv += 0.250000 <i,j||b,a>_abab t2_abab(b,a,i,j) s1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") += 0.250000 * scalar3 * s1_xb_Lv("I,e");

            // sigmas1_xb_Lv += 0.250000 <j,i||a,b>_bbbb t2_bbbb(a,b,j,i) s1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") += 0.250000 * scalar4 * s1_xb_Lv("I,e");

            // sigmas1_xb_Lv += -0.500000 <j,i||a,b>_abab t2_abab(a,e,j,i) s1_b(b) 
            // flops: o2v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmas1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("j,i,a,b") * t2_abab_vvoo("a,e,j,i") * s1_xb_Lv("I,b");

            // sigmas1_xb_Lv += -0.500000 <i,j||a,b>_abab t2_abab(a,e,i,j) s1_b(b) 
            // flops: o2v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmas1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("i,j,a,b") * t2_abab_vvoo("a,e,i,j") * s1_xb_Lv("I,b");

            // sigmas1_xb_Lv += -0.500000 <j,i||a,b>_bbbb t2_bbbb(a,e,j,i) s1_b(b) 
            // flops: o2v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmas1_xb_Lv("I,e") -= 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * t2_bbbb_vvoo("a,e,j,i") * s1_xb_Lv("I,b");

            // sigmas1_xb_Lv += 1.000000 <i,e||a,b>_abab r1_b(b) u1_aa(a,i) 
            // flops: o1v3L1: 1, o1v2L1: 1, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmas1_xb_Lv("I,e") -= V_blks_["baab_vovv"]("e,i,a,b") * r1_xb_Lv("I,b") * u1_aa_vo("a,i");

            // sigmas1_xb_Lv += 1.000000 <i,e||a,b>_bbbb r1_b(b) u1_bb(a,i) 
            // flops: o1v3L1: 1, o1v2L1: 1, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmas1_xb_Lv("I,e") -= V_blks_["bbbb_vovv"]("e,i,a,b") * r1_xb_Lv("I,b") * u1_bb_vo("a,i");

            // sigmas1_xb_Lv += -1.000000 <i,j||a,b>_abab r2_bbb(b,e,j) u1_aa(a,i) 
            // flops: o2v3L1: 1, o1v2L1: 1, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmas1_xb_Lv("I,e") -= V_blks_["abab_oovv"]("i,j,a,b") * r2_xbbb_Lvvo("I,b,e,j") * u1_aa_vo("a,i");

            // sigmas1_xb_Lv += 1.000000 <j,i||a,b>_bbbb r2_bbb(b,e,j) u1_bb(a,i) 
            // flops: o2v3L1: 1, o1v2L1: 1, o0v1L1: 1 | mem: o1v2L1: 1, o0v1L1: 2, 
            sigmas1_xb_Lv("I,e") += V_blks_["bbbb_oovv"]("j,i,a,b") * r2_xbbb_Lvvo("I,b,e,j") * u1_bb_vo("a,i");

            // sigmas1_xb_Lv += 0.500000 <j,i||a,b>_bbbb r2_bbb(a,b,j) u1_bb(e,i) 
            // flops: o2v2L1: 1, o1v1L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v0L1: 1, 
            sigmas1_xb_Lv("I,e") += 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r2_xbbb_Lvvo("I,a,b,j") * u1_bb_vo("e,i");
        }

        if (include_u1_ && include_u2_) {

            // sigmas1_xb_Lv += 0.250000 <j,i||a,b>_aaaa r1_b(e) u2_aaaa(a,b,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") += 0.250000 * scalar14 * r1_xb_Lv("I,e");

            // sigmas1_xb_Lv += 0.250000 <j,i||a,b>_abab r1_b(e) u2_abab(a,b,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") += 0.250000 * scalar15 * r1_xb_Lv("I,e");

            // sigmas1_xb_Lv += 0.250000 <i,j||a,b>_abab r1_b(e) u2_abab(a,b,i,j) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") += 0.250000 * scalar15 * r1_xb_Lv("I,e");

            // sigmas1_xb_Lv += 0.250000 <j,i||b,a>_abab r1_b(e) u2_abab(b,a,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") += 0.250000 * scalar15 * r1_xb_Lv("I,e");

            // sigmas1_xb_Lv += 0.250000 <i,j||b,a>_abab r1_b(e) u2_abab(b,a,i,j) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") += 0.250000 * scalar15 * r1_xb_Lv("I,e");

            // sigmas1_xb_Lv += 0.250000 <j,i||a,b>_bbbb r1_b(e) u2_bbbb(a,b,j,i) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmas1_xb_Lv("I,e") += 0.250000 * scalar16 * r1_xb_Lv("I,e");

            // sigmas1_xb_Lv += -0.500000 <j,i||a,b>_abab r1_b(b) u2_abab(a,e,j,i) 
            // flops: o2v2L1: 2, o0v1L1: 1 | mem: o2v1L1: 1, o0v1L1: 2, 
            sigmas1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("j,i,a,b") * r1_xb_Lv("I,b") * u2_abab_vvoo("a,e,j,i");

            // sigmas1_xb_Lv += -0.500000 <i,j||a,b>_abab r1_b(b) u2_abab(a,e,i,j) 
            // flops: o2v2L1: 2, o0v1L1: 1 | mem: o2v1L1: 1, o0v1L1: 2, 
            sigmas1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("i,j,a,b") * r1_xb_Lv("I,b") * u2_abab_vvoo("a,e,i,j");

            // sigmas1_xb_Lv += -0.500000 <j,i||a,b>_bbbb r1_b(b) u2_bbbb(a,e,j,i) 
            // flops: o2v2L1: 2, o0v1L1: 1 | mem: o2v1L1: 1, o0v1L1: 2, 
            sigmas1_xb_Lv("I,e") -= 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r1_xb_Lv("I,b") * u2_bbbb_vvoo("a,e,j,i");
        }




/// ****** sigmas2_aaa(e,f,n) ****** ///



        if (include_u2_) {

            // sigmas2_xaaa_Lvvo += 1.000000 s2_aaa(e,f,n) w0 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += s2_xaaa_Lvvo("I,e,f,n") * w0;
        }

        if (include_u0_ && include_u2_) {

            // sigmas2_xaaa_Lvvo += 1.000000 r2_aaa(e,f,n) u0 w0 
            // flops: o1v2L1: 2, o1v2: 1 | mem: o1v2L1: 1, o1v2: 2, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += r2_xaaa_Lvvo("I,e,f,n") * u0 * w0;
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) r1_a(f) u1_aa(e,n) w0 
            // flops: o1v2L1: 4 | mem: o1v2L1: 2, o1v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = r1_xa_Lv("I,f") * u1_aa_vo("e,n") * w0;
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) f_aa(e,n) s1_a(f) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = F_blks_["aa_vo"]("e,n") * s1_xa_Lv("I,f");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmas2_xaaa_Lvvo += 1.000000 f_aa(i,i) s2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += scalar5 * s2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += 1.000000 f_bb(i,i) s2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += scalar6 * s2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += -1.000000 f_aa(i,n) s2_aaa(e,f,i) 
            // flops: o2v2L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= F_blks_["aa_oo"]("i,n") * s2_xaaa_Lvvo("I,e,f,i");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) f_aa(e,a) s2_aaa(a,f,n) 
            // flops: o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = F_blks_["aa_vv"]("e,a") * s2_xaaa_Lvvo("I,a,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) f_aa(i,a) t2_aaaa(a,e,n,i) s1_a(f) 
            // flops: o1v2L1: 3, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = F_blks_["aa_ov"]("i,a") * t2_aaaa_vvoo("a,e,n,i") * s1_xa_Lv("I,f");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) f_bb(i,a) t2_abab(e,a,n,i) s1_a(f) 
            // flops: o1v2L1: 3, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = F_blks_["bb_ov"]("i,a") * t2_abab_vvoo("e,a,n,i") * s1_xa_Lv("I,f");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 f_aa(i,a) t2_aaaa(e,f,n,i) s1_a(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += F_blks_["aa_ov"]("i,a") * t2_aaaa_vvoo("e,f,n,i") * s1_xa_Lv("I,a");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) f_aa(i,n) r1_a(f) u1_aa(e,i) 
            // flops: o2v2L1: 1, o1v2L1: 2, o2v1L1: 1 | mem: o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = F_blks_["aa_oo"]("i,n") * r1_xa_Lv("I,f") * u1_aa_vo("e,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) f_aa(e,a) r1_a(f) u1_aa(a,n) 
            // flops: o1v3L1: 1, o0v3L1: 1, o1v2L1: 2 | mem: o0v3L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = F_blks_["aa_vv"]("e,a") * r1_xa_Lv("I,f") * u1_aa_vo("a,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 f_aa(i,a) r2_aaa(e,f,n) u1_aa(a,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += scalar12 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += 1.000000 f_bb(i,a) r2_aaa(e,f,n) u1_bb(a,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += scalar13 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += -1.000000 f_aa(i,a) r2_aaa(e,f,i) u1_aa(a,n) 
            // flops: o1v3L1: 2, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= F_blks_["aa_ov"]("i,a") * r2_xaaa_Lvvo("I,e,f,i") * u1_aa_vo("a,n");

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) f_aa(i,a) r2_aaa(a,f,n) u1_aa(e,i) 
            // flops: o2v2L1: 2, o1v2L1: 2 | mem: o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = F_blks_["aa_ov"]("i,a") * r2_xaaa_Lvvo("I,a,f,n") * u1_aa_vo("e,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) f_aa(i,a) r1_a(f) u2_aaaa(a,e,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = F_blks_["aa_ov"]("i,a") * r1_xa_Lv("I,f") * u2_aaaa_vvoo("a,e,n,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) f_bb(i,a) r1_a(f) u2_abab(e,a,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = F_blks_["bb_ov"]("i,a") * r1_xa_Lv("I,f") * u2_abab_vvoo("e,a,n,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 f_aa(i,a) r1_a(a) u2_aaaa(e,f,n,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o1v1L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += F_blks_["aa_ov"]("i,a") * r1_xa_Lv("I,a") * u2_aaaa_vvoo("e,f,n,i");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) d+_aa(e,n) r1_a(f) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_vo("e,n") * r1_xa_Lv("I,f");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 d+_aa(i,i) r2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= scalar17 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += -1.000000 d+_bb(i,i) r2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= scalar18 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += 1.000000 d+_aa(i,n) r2_aaa(e,f,i) 
            // flops: o2v2L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += dp_aa_oo("i,n") * r2_xaaa_Lvvo("I,e,f,i");

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) d+_aa(e,a) r2_aaa(a,f,n) 
            // flops: o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_vv("e,a") * r2_xaaa_Lvvo("I,a,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) d+_aa(i,a) r1_a(f) t2_aaaa(a,e,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("i,a") * r1_xa_Lv("I,f") * t2_aaaa_vvoo("a,e,n,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) d+_bb(i,a) r1_a(f) t2_abab(e,a,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_bb_ov("i,a") * r1_xa_Lv("I,f") * t2_abab_vvoo("e,a,n,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 d+_aa(i,a) r1_a(a) t2_aaaa(e,f,n,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o1v1L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= dp_aa_ov("i,a") * r1_xa_Lv("I,a") * t2_aaaa_vvoo("e,f,n,i");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) d-_aa(e,n) u0 s1_a(f) 
            // flops: o1v2L1: 3, o1v1: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_vo("e,n") * u0 * s1_xa_Lv("I,f");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u0_ && include_u2_) {

            // sigmas2_xaaa_Lvvo += -1.000000 d-_aa(i,i) u0 s2_aaa(e,f,n) 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= u0 * scalar7 * s2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += -1.000000 d-_bb(i,i) u0 s2_aaa(e,f,n) 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= u0 * scalar8 * s2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += 1.000000 d-_aa(i,n) u0 s2_aaa(e,f,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o2v0: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += dp_aa_oo("i,n") * u0 * s2_xaaa_Lvvo("I,e,f,i");

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(e,a) u0 s2_aaa(a,f,n) 
            // flops: o1v3L1: 1, o1v2L1: 2, o0v2: 1 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_vv("e,a") * u0 * s2_xaaa_Lvvo("I,a,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) d-_aa(i,i) u1_aa(e,n) s1_a(f) 
            // flops: o1v2L1: 2, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = scalar7 * s1_xa_Lv("I,f") * u1_aa_vo("e,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) d-_bb(i,i) u1_aa(e,n) s1_a(f) 
            // flops: o1v2L1: 2, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = scalar8 * s1_xa_Lv("I,f") * u1_aa_vo("e,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -2.000000 P(e,f) d-_aa(i,n) u1_aa(e,i) s1_a(f) 
            // flops: o1v2L1: 3, o2v1: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 2.000000 * dp_aa_oo("i,n") * u1_aa_vo("e,i") * s1_xa_Lv("I,f");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 2.000000 P(e,f) d-_aa(e,a) u1_aa(a,n) s1_a(f) 
            // flops: o1v2L1: 3, o1v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 2.000000 * dp_aa_vv("e,a") * u1_aa_vo("a,n") * s1_xa_Lv("I,f");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(e,a) u1_aa(f,n) s1_a(a) 
            // flops: o1v3L1: 1, o1v2L1: 2, o1v3: 1 | mem: o1v2L1: 2, o1v3: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_vv("e,a") * u1_aa_vo("f,n") * s1_xa_Lv("I,a");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -2.000000 d-_aa(i,a) u1_aa(a,i) s2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= 2.000000 * scalar0 * s2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += -2.000000 d-_bb(i,a) u1_bb(a,i) s2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= 2.000000 * scalar1 * s2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += 2.000000 d-_aa(i,a) u1_aa(a,n) s2_aaa(e,f,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o2v1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += 2.000000 * dp_aa_ov("i,a") * u1_aa_vo("a,n") * s2_xaaa_Lvvo("I,e,f,i");

            // sigmas2_xaaa_Lvvo += 2.000000 P(e,f) d-_aa(i,a) u1_aa(e,i) s2_aaa(a,f,n) 
            // flops: o1v3L1: 1, o1v2L1: 2, o1v2: 1 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 2.000000 * dp_aa_ov("i,a") * u1_aa_vo("e,i") * s2_xaaa_Lvvo("I,a,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(i,a) u1_aa(e,n) s2_aaa(a,f,i) 
            // flops: o2v3L1: 1, o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("i,a") * u1_aa_vo("e,n") * s2_xaaa_Lvvo("I,a,f,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -2.000000 P(e,f) d-_aa(i,a) u2_aaaa(a,e,n,i) s1_a(f) 
            // flops: o1v2L1: 3, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 2.000000 * dp_aa_ov("i,a") * u2_aaaa_vvoo("a,e,n,i") * s1_xa_Lv("I,f");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 2.000000 P(e,f) d-_bb(i,a) u2_abab(e,a,n,i) s1_a(f) 
            // flops: o1v2L1: 3, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 2.000000 * dp_bb_ov("i,a") * u2_abab_vvoo("e,a,n,i") * s1_xa_Lv("I,f");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -2.000000 d-_aa(i,a) u2_aaaa(e,f,n,i) s1_a(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= 2.000000 * dp_aa_ov("i,a") * u2_aaaa_vvoo("e,f,n,i") * s1_xa_Lv("I,a");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(i,a) t2_aaaa(a,e,n,i) u0 s1_a(f) 
            // flops: o1v2L1: 3, o2v2: 1, o1v1: 1 | mem: o1v2L1: 2, o1v1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("i,a") * t2_aaaa_vvoo("a,e,n,i") * u0 * s1_xa_Lv("I,f");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) d-_bb(i,a) t2_abab(e,a,n,i) u0 s1_a(f) 
            // flops: o1v2L1: 3, o2v2: 1, o1v1: 1 | mem: o1v2L1: 2, o1v1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_bb_ov("i,a") * t2_abab_vvoo("e,a,n,i") * u0 * s1_xa_Lv("I,f");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 d-_aa(i,a) t2_aaaa(e,f,n,i) u0 s1_a(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1, o1v3: 1 | mem: o1v2L1: 2, o1v3: 2, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= dp_aa_ov("i,a") * t2_aaaa_vvoo("e,f,n,i") * u0 * s1_xa_Lv("I,a");

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(i,n) r1_a(f) u0 u1_aa(e,i) 
            // flops: o1v2L1: 2, o2v1L1: 2, o2v2: 1 | mem: o1v2L1: 1, o2v1L1: 1, o1v2: 1, o2v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_oo("i,n") * r1_xa_Lv("I,f") * u0 * u1_aa_vo("e,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) d-_aa(e,a) r1_a(f) u0 u1_aa(a,n) 
            // flops: o0v3L1: 2, o1v2L1: 2, o1v3: 1 | mem: o0v3L1: 1, o1v2L1: 1, o0v3: 1, o1v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_vv("e,a") * r1_xa_Lv("I,f") * u0 * u1_aa_vo("a,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 d-_aa(i,a) r2_aaa(e,f,n) u0 u1_aa(a,i) 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= u0 * scalar0 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += -1.000000 d-_bb(i,a) r2_aaa(e,f,n) u0 u1_bb(a,i) 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= u0 * scalar1 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += 1.000000 d-_aa(i,a) r2_aaa(e,f,i) u0 u1_aa(a,n) 
            // flops: o1v3L1: 1, o0v3L1: 1, o1v2L1: 1, o1v3: 1 | mem: o0v3L1: 1, o1v2L1: 1, o0v3: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += dp_aa_ov("i,a") * r2_xaaa_Lvvo("I,e,f,i") * u0 * u1_aa_vo("a,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) d-_aa(i,a) r2_aaa(a,f,n) u0 u1_aa(e,i) 
            // flops: o2v2L1: 1, o1v2L1: 2, o2v1L1: 1, o2v2: 1 | mem: o1v2L1: 1, o2v1L1: 1, o1v2: 1, o2v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("i,a") * r2_xaaa_Lvvo("I,a,f,n") * u0 * u1_aa_vo("e,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u0_ && include_u2_) {

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(i,a) r1_a(f) u0 u2_aaaa(a,e,n,i) 
            // flops: o2v3: 1, o1v2L1: 4 | mem: o1v2L1: 2, o1v2: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("i,a") * r1_xa_Lv("I,f") * u0 * u2_aaaa_vvoo("a,e,n,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) d-_bb(i,a) r1_a(f) u0 u2_abab(e,a,n,i) 
            // flops: o2v3: 1, o1v2L1: 4 | mem: o1v2L1: 2, o1v2: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_bb_ov("i,a") * r1_xa_Lv("I,f") * u0 * u2_abab_vvoo("e,a,n,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 d-_aa(i,a) r1_a(a) u0 u2_aaaa(e,f,n,i) 
            // flops: o1v2L1: 1, o2v2: 1, o1v1L1: 1, o1v0L1: 1 | mem: o1v2L1: 1, o1v2: 1, o1v0L1: 1, o1v0: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= dp_aa_ov("i,a") * r1_xa_Lv("I,a") * u0 * u2_aaaa_vvoo("e,f,n,i");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) d-_aa(i,a) r1_a(f) u1_aa(a,i) u1_aa(e,n) 
            // flops: o1v2L1: 2, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = scalar0 * r1_xa_Lv("I,f") * u1_aa_vo("e,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) d-_bb(i,a) r1_a(f) u1_bb(a,i) u1_aa(e,n) 
            // flops: o1v2L1: 2, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = scalar1 * r1_xa_Lv("I,f") * u1_aa_vo("e,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -2.000000 P(e,f) d-_aa(i,a) r1_a(f) u1_aa(a,n) u1_aa(e,i) 
            // flops: o2v2L1: 2, o1v2L1: 3 | mem: o1v2L1: 3, o2v1L1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 2.000000 * dp_aa_ov("i,a") * r1_xa_Lv("I,f") * u1_aa_vo("a,n") * u1_aa_vo("e,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) d-_aa(i,a) r1_a(a) u1_aa(e,i) u1_aa(f,n) 
            // flops: o1v2L1: 3, o1v1L1: 2 | mem: o1v2L1: 2, o0v1L1: 1, o1v0L1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("i,a") * r1_xa_Lv("I,a") * u1_aa_vo("e,i") * u1_aa_vo("f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmas2_xaaa_Lvvo += -0.500000 <j,i||j,i>_aaaa s2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * scalar9 * s2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += -0.500000 <j,i||j,i>_abab s2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * scalar10 * s2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += -0.500000 <i,j||i,j>_abab s2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * scalar10 * s2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += -0.500000 <j,i||j,i>_bbbb s2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * scalar11 * s2_xaaa_Lvvo("I,e,f,n");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xaaa_Lvvo += 1.000000 <e,f||a,n>_aaaa s1_a(a) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += V_blks_["aaaa_vvvo"]("e,f,a,n") * s1_xa_Lv("I,a");
        }

        if (include_u2_) {

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) <i,e||a,n>_aaaa s2_aaa(a,f,i) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_vovo"]("e,i,a,n") * s2_xaaa_Lvvo("I,a,f,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 0.500000 <e,f||a,b>_aaaa s2_aaa(a,b,n) 
            // flops: o1v4L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["aaaa_vvvv"]("e,f,a,b") * s2_xaaa_Lvvo("I,a,b,n");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xaaa_Lvvo += 0.500000 P(e,f) <j,i||a,n>_aaaa t2_aaaa(a,e,j,i) s1_a(f) 
            // flops: o3v2: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["aaaa_oovo"]("j,i,a,n") * t2_aaaa_vvoo("a,e,j,i") * s1_xa_Lv("I,f");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 0.500000 P(e,f) <j,i||n,a>_abab t2_abab(e,a,j,i) s1_a(f) 
            // flops: o3v2: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abba_oovo"]("j,i,a,n") * t2_abab_vvoo("e,a,j,i") * s1_xa_Lv("I,f");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 0.500000 P(e,f) <i,j||n,a>_abab t2_abab(e,a,i,j) s1_a(f) 
            // flops: o3v2: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abba_oovo"]("i,j,a,n") * t2_abab_vvoo("e,a,i,j") * s1_xa_Lv("I,f");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 0.500000 <j,i||a,n>_aaaa t2_aaaa(e,f,j,i) s1_a(a) 
            // flops: o3v3: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["aaaa_oovo"]("j,i,a,n") * t2_aaaa_vvoo("e,f,j,i") * s1_xa_Lv("I,a");

            // sigmas2_xaaa_Lvvo += 0.500000 P(e,f) <i,e||a,b>_aaaa t2_aaaa(a,b,n,i) s1_a(f) 
            // flops: o2v3: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["aaaa_vovv"]("e,i,a,b") * t2_aaaa_vvoo("a,b,n,i") * s1_xa_Lv("I,f");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -0.500000 P(e,f) <e,i||a,b>_abab t2_abab(a,b,n,i) s1_a(f) 
            // flops: o2v3: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_vovv"]("e,i,a,b") * t2_abab_vvoo("a,b,n,i") * s1_xa_Lv("I,f");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -0.500000 P(e,f) <e,i||b,a>_abab t2_abab(b,a,n,i) s1_a(f) 
            // flops: o2v3: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_vovv"]("e,i,b,a") * t2_abab_vvoo("b,a,n,i") * s1_xa_Lv("I,f");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) <i,e||a,b>_aaaa t2_aaaa(a,f,n,i) s1_a(b) 
            // flops: o2v4: 1, o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v3: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_vovv"]("e,i,a,b") * t2_aaaa_vvoo("a,f,n,i") * s1_xa_Lv("I,b");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) <e,i||b,a>_abab t2_abab(f,a,n,i) s1_a(b) 
            // flops: o2v4: 1, o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v3: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_vovv"]("e,i,b,a") * t2_abab_vvoo("f,a,n,i") * s1_xa_Lv("I,b");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmas2_xaaa_Lvvo += 0.250000 <j,i||a,b>_aaaa t2_aaaa(a,b,j,i) s2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar2 * s2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += 0.250000 <j,i||a,b>_abab t2_abab(a,b,j,i) s2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar3 * s2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += 0.250000 <i,j||a,b>_abab t2_abab(a,b,i,j) s2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar3 * s2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += 0.250000 <j,i||b,a>_abab t2_abab(b,a,j,i) s2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar3 * s2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += 0.250000 <i,j||b,a>_abab t2_abab(b,a,i,j) s2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar3 * s2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += 0.250000 <j,i||a,b>_bbbb t2_bbbb(a,b,j,i) s2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar4 * s2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += -0.500000 <j,i||a,b>_aaaa t2_aaaa(a,b,n,i) s2_aaa(e,f,j) 
            // flops: o2v2L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * t2_aaaa_vvoo("a,b,n,i") * s2_xaaa_Lvvo("I,e,f,j");

            // sigmas2_xaaa_Lvvo += -0.500000 <j,i||a,b>_abab t2_abab(a,b,n,i) s2_aaa(e,f,j) 
            // flops: o2v2L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("j,i,a,b") * t2_abab_vvoo("a,b,n,i") * s2_xaaa_Lvvo("I,e,f,j");

            // sigmas2_xaaa_Lvvo += -0.500000 <j,i||b,a>_abab t2_abab(b,a,n,i) s2_aaa(e,f,j) 
            // flops: o2v2L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("j,i,b,a") * t2_abab_vvoo("b,a,n,i") * s2_xaaa_Lvvo("I,e,f,j");

            // sigmas2_xaaa_Lvvo += -0.500000 P(e,f) <j,i||a,b>_aaaa t2_aaaa(a,e,j,i) s2_aaa(b,f,n) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * t2_aaaa_vvoo("a,e,j,i") * s2_xaaa_Lvvo("I,b,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -0.500000 P(e,f) <j,i||b,a>_abab t2_abab(e,a,j,i) s2_aaa(b,f,n) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("j,i,b,a") * t2_abab_vvoo("e,a,j,i") * s2_xaaa_Lvvo("I,b,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -0.500000 P(e,f) <i,j||b,a>_abab t2_abab(e,a,i,j) s2_aaa(b,f,n) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("i,j,b,a") * t2_abab_vvoo("e,a,i,j") * s2_xaaa_Lvvo("I,b,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) <j,i||a,b>_aaaa t2_aaaa(a,e,n,i) s2_aaa(b,f,j) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_oovv"]("j,i,a,b") * t2_aaaa_vvoo("a,e,n,i") * s2_xaaa_Lvvo("I,b,f,j");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) <j,i||b,a>_abab t2_abab(e,a,n,i) s2_aaa(b,f,j) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("j,i,b,a") * t2_abab_vvoo("e,a,n,i") * s2_xaaa_Lvvo("I,b,f,j");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 0.250000 <j,i||a,b>_aaaa t2_aaaa(e,f,j,i) s2_aaa(a,b,n) 
            // flops: o1v4L1: 1, o2v4: 1, o1v2L1: 1 | mem: o0v4: 1, o1v2L1: 2, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += 0.250000 * V_blks_["aaaa_oovv"]("j,i,a,b") * t2_aaaa_vvoo("e,f,j,i") * s2_xaaa_Lvvo("I,a,b,n");

            // sigmas2_xaaa_Lvvo += -0.500000 <j,i||a,b>_aaaa t2_aaaa(e,f,n,i) s2_aaa(a,b,j) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * t2_aaaa_vvoo("e,f,n,i") * s2_xaaa_Lvvo("I,a,b,j");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) <i,e||a,n>_aaaa r1_a(f) u1_aa(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_vovo"]("e,i,a,n") * r1_xa_Lv("I,f") * u1_aa_vo("a,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) <e,i||n,a>_abab r1_a(f) u1_bb(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abba_vovo"]("e,i,a,n") * r1_xa_Lv("I,f") * u1_bb_vo("a,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) <i,e||a,n>_aaaa r1_a(a) u1_aa(f,i) 
            // flops: o2v2L1: 2, o1v2L1: 2 | mem: o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_vovo"]("e,i,a,n") * r1_xa_Lv("I,a") * u1_aa_vo("f,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 <e,f||a,b>_aaaa r1_a(b) u1_aa(a,n) 
            // flops: o0v4L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= V_blks_["aaaa_vvvv"]("e,f,a,b") * r1_xa_Lv("I,b") * u1_aa_vo("a,n");

            // sigmas2_xaaa_Lvvo += 1.000000 <j,i||a,n>_aaaa r2_aaa(e,f,j) u1_aa(a,i) 
            // flops: o3v3L1: 1, o2v3L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += V_blks_["aaaa_oovo"]("j,i,a,n") * r2_xaaa_Lvvo("I,e,f,j") * u1_aa_vo("a,i");

            // sigmas2_xaaa_Lvvo += -1.000000 <j,i||n,a>_abab r2_aaa(e,f,j) u1_bb(a,i) 
            // flops: o3v3L1: 1, o2v3L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += V_blks_["abba_oovo"]("j,i,a,n") * r2_xaaa_Lvvo("I,e,f,j") * u1_bb_vo("a,i");

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) <j,i||a,n>_aaaa r2_aaa(a,f,j) u1_aa(e,i) 
            // flops: o3v2L1: 1, o2v2L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_oovo"]("j,i,a,n") * r2_xaaa_Lvvo("I,a,f,j") * u1_aa_vo("e,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) <i,e||a,b>_aaaa r2_aaa(b,f,n) u1_aa(a,i) 
            // flops: o2v4L1: 1, o2v3L1: 1, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_vovv"]("e,i,a,b") * r2_xaaa_Lvvo("I,b,f,n") * u1_aa_vo("a,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) <e,i||b,a>_abab r2_aaa(b,f,n) u1_bb(a,i) 
            // flops: o2v4L1: 1, o2v3L1: 1, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_vovv"]("e,i,b,a") * r2_xaaa_Lvvo("I,b,f,n") * u1_bb_vo("a,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) <i,e||a,b>_aaaa r2_aaa(b,f,i) u1_aa(a,n) 
            // flops: o1v4L1: 1, o1v3L1: 1, o1v2L1: 2 | mem: o0v3L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_vovv"]("e,i,a,b") * r2_xaaa_Lvvo("I,b,f,i") * u1_aa_vo("a,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 0.500000 P(e,f) <i,e||a,b>_aaaa r2_aaa(a,b,n) u1_aa(f,i) 
            // flops: o2v3L1: 1, o2v2L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["aaaa_vovv"]("e,i,a,b") * r2_xaaa_Lvvo("I,a,b,n") * u1_aa_vo("f,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmas2_xaaa_Lvvo += 0.500000 P(e,f) <j,i||a,n>_aaaa r1_a(f) u2_aaaa(a,e,j,i) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["aaaa_oovo"]("j,i,a,n") * r1_xa_Lv("I,f") * u2_aaaa_vvoo("a,e,j,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 0.500000 P(e,f) <j,i||n,a>_abab r1_a(f) u2_abab(e,a,j,i) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abba_oovo"]("j,i,a,n") * r1_xa_Lv("I,f") * u2_abab_vvoo("e,a,j,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 0.500000 P(e,f) <i,j||n,a>_abab r1_a(f) u2_abab(e,a,i,j) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abba_oovo"]("i,j,a,n") * r1_xa_Lv("I,f") * u2_abab_vvoo("e,a,i,j");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 0.500000 <j,i||a,n>_aaaa r1_a(a) u2_aaaa(e,f,j,i) 
            // flops: o3v2L1: 1, o3v1L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o3v0L1: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["aaaa_oovo"]("j,i,a,n") * r1_xa_Lv("I,a") * u2_aaaa_vvoo("e,f,j,i");

            // sigmas2_xaaa_Lvvo += 0.500000 P(e,f) <i,e||a,b>_aaaa r1_a(f) u2_aaaa(a,b,n,i) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 2 | mem: o1v4L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["aaaa_vovv"]("e,i,a,b") * r1_xa_Lv("I,f") * u2_aaaa_vvoo("a,b,n,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -0.500000 P(e,f) <e,i||a,b>_abab r1_a(f) u2_abab(a,b,n,i) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 2 | mem: o1v4L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_vovv"]("e,i,a,b") * r1_xa_Lv("I,f") * u2_abab_vvoo("a,b,n,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -0.500000 P(e,f) <e,i||b,a>_abab r1_a(f) u2_abab(b,a,n,i) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 2 | mem: o1v4L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_vovv"]("e,i,b,a") * r1_xa_Lv("I,f") * u2_abab_vvoo("b,a,n,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) <i,e||a,b>_aaaa r1_a(b) u2_aaaa(a,f,n,i) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_vovv"]("e,i,a,b") * r1_xa_Lv("I,b") * u2_aaaa_vvoo("a,f,n,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) <e,i||b,a>_abab r1_a(b) u2_abab(f,a,n,i) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_vovv"]("e,i,b,a") * r1_xa_Lv("I,b") * u2_abab_vvoo("f,a,n,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 0.250000 <j,i||a,b>_aaaa r2_aaa(e,f,n) u2_aaaa(a,b,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar14 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += 0.250000 <j,i||a,b>_abab r2_aaa(e,f,n) u2_abab(a,b,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar15 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += 0.250000 <i,j||a,b>_abab r2_aaa(e,f,n) u2_abab(a,b,i,j) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar15 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += 0.250000 <j,i||b,a>_abab r2_aaa(e,f,n) u2_abab(b,a,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar15 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += 0.250000 <i,j||b,a>_abab r2_aaa(e,f,n) u2_abab(b,a,i,j) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar15 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += 0.250000 <j,i||a,b>_bbbb r2_aaa(e,f,n) u2_bbbb(a,b,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar16 * r2_xaaa_Lvvo("I,e,f,n");

            // sigmas2_xaaa_Lvvo += -0.500000 <j,i||a,b>_aaaa r2_aaa(e,f,j) u2_aaaa(a,b,n,i) 
            // flops: o2v4L1: 2, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaa_Lvvo("I,e,f,j") * u2_aaaa_vvoo("a,b,n,i");

            // sigmas2_xaaa_Lvvo += -0.500000 <j,i||a,b>_abab r2_aaa(e,f,j) u2_abab(a,b,n,i) 
            // flops: o2v4L1: 2, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("j,i,a,b") * r2_xaaa_Lvvo("I,e,f,j") * u2_abab_vvoo("a,b,n,i");

            // sigmas2_xaaa_Lvvo += -0.500000 <j,i||b,a>_abab r2_aaa(e,f,j) u2_abab(b,a,n,i) 
            // flops: o2v4L1: 2, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("j,i,b,a") * r2_xaaa_Lvvo("I,e,f,j") * u2_abab_vvoo("b,a,n,i");

            // sigmas2_xaaa_Lvvo += -0.500000 P(e,f) <j,i||a,b>_aaaa r2_aaa(b,f,n) u2_aaaa(a,e,j,i) 
            // flops: o3v3L1: 2, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaa_Lvvo("I,b,f,n") * u2_aaaa_vvoo("a,e,j,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -0.500000 P(e,f) <j,i||b,a>_abab r2_aaa(b,f,n) u2_abab(e,a,j,i) 
            // flops: o3v3L1: 2, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("j,i,b,a") * r2_xaaa_Lvvo("I,b,f,n") * u2_abab_vvoo("e,a,j,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -0.500000 P(e,f) <i,j||b,a>_abab r2_aaa(b,f,n) u2_abab(e,a,i,j) 
            // flops: o3v3L1: 2, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("i,j,b,a") * r2_xaaa_Lvvo("I,b,f,n") * u2_abab_vvoo("e,a,i,j");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) <j,i||a,b>_aaaa r2_aaa(b,f,j) u2_aaaa(a,e,n,i) 
            // flops: o2v3L1: 2, o1v2L1: 2 | mem: o1v2L1: 3, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaa_Lvvo("I,b,f,j") * u2_aaaa_vvoo("a,e,n,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) <j,i||b,a>_abab r2_aaa(b,f,j) u2_abab(e,a,n,i) 
            // flops: o2v3L1: 2, o1v2L1: 2 | mem: o1v2L1: 3, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("j,i,b,a") * r2_xaaa_Lvvo("I,b,f,j") * u2_abab_vvoo("e,a,n,i");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 0.250000 <j,i||a,b>_aaaa r2_aaa(a,b,n) u2_aaaa(e,f,j,i) 
            // flops: o3v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o3v0L1: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += 0.250000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaa_Lvvo("I,a,b,n") * u2_aaaa_vvoo("e,f,j,i");

            // sigmas2_xaaa_Lvvo += -0.500000 <j,i||a,b>_aaaa r2_aaa(a,b,j) u2_aaaa(e,f,n,i) 
            // flops: o2v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaa_Lvvo("I,a,b,j") * u2_aaaa_vvoo("e,f,n,i");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xaaa_Lvvo += 0.500000 P(e,f) <j,i||a,b>_aaaa r1_a(f) t2_aaaa(a,b,n,i) u1_aa(e,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r1_xa_Lv("I,f") * t2_aaaa_vvoo("a,b,n,i") * u1_aa_vo("e,j");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 0.500000 P(e,f) <j,i||a,b>_abab r1_a(f) t2_abab(a,b,n,i) u1_aa(e,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("j,i,a,b") * r1_xa_Lv("I,f") * t2_abab_vvoo("a,b,n,i") * u1_aa_vo("e,j");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 0.500000 P(e,f) <j,i||b,a>_abab r1_a(f) t2_abab(b,a,n,i) u1_aa(e,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("j,i,b,a") * r1_xa_Lv("I,f") * t2_abab_vvoo("b,a,n,i") * u1_aa_vo("e,j");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 0.500000 P(e,f) <j,i||a,b>_aaaa r1_a(f) t2_aaaa(a,e,j,i) u1_aa(b,n) 
            // flops: o2v4L1: 1, o2v3L1: 1, o1v3L1: 1, o1v2L1: 2 | mem: o2v3L1: 1, o0v3L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r1_xa_Lv("I,f") * t2_aaaa_vvoo("a,e,j,i") * u1_aa_vo("b,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 0.500000 P(e,f) <j,i||b,a>_abab r1_a(f) t2_abab(e,a,j,i) u1_aa(b,n) 
            // flops: o2v4L1: 1, o2v3L1: 1, o1v3L1: 1, o1v2L1: 2 | mem: o2v3L1: 1, o0v3L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("j,i,b,a") * r1_xa_Lv("I,f") * t2_abab_vvoo("e,a,j,i") * u1_aa_vo("b,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 0.500000 P(e,f) <i,j||b,a>_abab r1_a(f) t2_abab(e,a,i,j) u1_aa(b,n) 
            // flops: o2v4L1: 1, o2v3L1: 1, o1v3L1: 1, o1v2L1: 2 | mem: o2v3L1: 1, o0v3L1: 1, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("i,j,b,a") * r1_xa_Lv("I,f") * t2_abab_vvoo("e,a,i,j") * u1_aa_vo("b,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) <j,i||a,b>_aaaa r1_a(f) t2_aaaa(a,e,n,i) u1_aa(b,j) 
            // flops: o3v4L1: 1, o2v3L1: 2, o1v2L1: 2 | mem: o2v3L1: 2, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_oovv"]("j,i,a,b") * r1_xa_Lv("I,f") * t2_aaaa_vvoo("a,e,n,i") * u1_aa_vo("b,j");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) <i,j||a,b>_abab r1_a(f) t2_aaaa(a,e,n,i) u1_bb(b,j) 
            // flops: o3v4L1: 1, o2v3L1: 2, o1v2L1: 2 | mem: o2v3L1: 2, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("i,j,a,b") * r1_xa_Lv("I,f") * t2_aaaa_vvoo("a,e,n,i") * u1_bb_vo("b,j");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -1.000000 P(e,f) <j,i||b,a>_abab r1_a(f) t2_abab(e,a,n,i) u1_aa(b,j) 
            // flops: o3v4L1: 1, o2v3L1: 2, o1v2L1: 2 | mem: o2v3L1: 2, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("j,i,b,a") * r1_xa_Lv("I,f") * t2_abab_vvoo("e,a,n,i") * u1_aa_vo("b,j");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) <j,i||a,b>_bbbb r1_a(f) t2_abab(e,a,n,i) u1_bb(b,j) 
            // flops: o3v4L1: 1, o2v3L1: 2, o1v2L1: 2 | mem: o2v3L1: 2, o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["bbbb_oovv"]("j,i,a,b") * r1_xa_Lv("I,f") * t2_abab_vvoo("e,a,n,i") * u1_bb_vo("b,j");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) <j,i||a,b>_aaaa r1_a(b) t2_aaaa(a,e,n,i) u1_aa(f,j) 
            // flops: o3v2L1: 1, o2v2L1: 2, o1v2L1: 2 | mem: o1v2L1: 2, o2v1L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_oovv"]("j,i,a,b") * r1_xa_Lv("I,b") * t2_aaaa_vvoo("a,e,n,i") * u1_aa_vo("f,j");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += 1.000000 P(e,f) <j,i||b,a>_abab r1_a(b) t2_abab(e,a,n,i) u1_aa(f,j) 
            // flops: o3v2L1: 1, o2v2L1: 2, o1v2L1: 2 | mem: o1v2L1: 2, o2v1L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("j,i,b,a") * r1_xa_Lv("I,b") * t2_abab_vvoo("e,a,n,i") * u1_aa_vo("f,j");
            sigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmas2_xaaa_Lvvo += -0.500000 <j,i||a,b>_aaaa r1_a(b) t2_aaaa(e,f,j,i) u1_aa(a,n) 
            // flops: o2v3L1: 1, o1v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r1_xa_Lv("I,b") * t2_aaaa_vvoo("e,f,j,i") * u1_aa_vo("a,n");

            // sigmas2_xaaa_Lvvo += 1.000000 <j,i||a,b>_aaaa r1_a(b) t2_aaaa(e,f,n,i) u1_aa(a,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += V_blks_["aaaa_oovv"]("j,i,a,b") * r1_xa_Lv("I,b") * t2_aaaa_vvoo("e,f,n,i") * u1_aa_vo("a,j");

            // sigmas2_xaaa_Lvvo += 1.000000 <i,j||b,a>_abab r1_a(b) t2_aaaa(e,f,n,i) u1_bb(a,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xaaa_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("i,j,b,a") * r1_xa_Lv("I,b") * t2_aaaa_vvoo("e,f,n,i") * u1_bb_vo("a,j");
        }




/// ****** sigmas2_abb(e,f,n) ****** ///



        if (include_u1_ && include_u2_) {

            // sigmas2_xabb_Lvvo += 1.000000 r1_a(e) u1_bb(f,n) w0 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, o1v2: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += r1_xa_Lv("I,e") * u1_bb_vo("f,n") * w0;

            // sigmas2_xabb_Lvvo += 1.000000 f_bb(f,n) s1_a(e) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") += F_blks_["bb_vo"]("f,n") * s1_xa_Lv("I,e");

            // sigmas2_xabb_Lvvo += 1.000000 f_aa(i,a) t2_abab(a,f,i,n) s1_a(e) 
            // flops: o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += F_blks_["aa_ov"]("i,a") * t2_abab_vvoo("a,f,i,n") * s1_xa_Lv("I,e");

            // sigmas2_xabb_Lvvo += -1.000000 f_bb(i,a) t2_bbbb(a,f,n,i) s1_a(e) 
            // flops: o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= F_blks_["bb_ov"]("i,a") * t2_bbbb_vvoo("a,f,n,i") * s1_xa_Lv("I,e");

            // sigmas2_xabb_Lvvo += -1.000000 f_aa(i,a) t2_abab(e,f,i,n) s1_a(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= F_blks_["aa_ov"]("i,a") * t2_abab_vvoo("e,f,i,n") * s1_xa_Lv("I,a");

            // sigmas2_xabb_Lvvo += -1.000000 f_bb(i,n) r1_a(e) u1_bb(f,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o2v1L1: 1 | mem: o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= F_blks_["bb_oo"]("i,n") * r1_xa_Lv("I,e") * u1_bb_vo("f,i");

            // sigmas2_xabb_Lvvo += 1.000000 f_bb(f,a) r1_a(e) u1_bb(a,n) 
            // flops: o1v3L1: 1, o0v3L1: 1, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") += F_blks_["bb_vv"]("f,a") * r1_xa_Lv("I,e") * u1_bb_vo("a,n");
        }

        if (include_u2_) {

            // sigmas2_xabb_Lvvo += 1.000000 f_aa(i,a) r1_a(e) u2_abab(a,f,i,n) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            sigmas2_xabb_Lvvo("I,e,f,n") += F_blks_["aa_ov"]("i,a") * r1_xa_Lv("I,e") * u2_abab_vvoo("a,f,i,n");

            // sigmas2_xabb_Lvvo += -1.000000 f_bb(i,a) r1_a(e) u2_bbbb(a,f,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= F_blks_["bb_ov"]("i,a") * r1_xa_Lv("I,e") * u2_bbbb_vvoo("a,f,n,i");

            // sigmas2_xabb_Lvvo += -1.000000 f_aa(i,a) r1_a(a) u2_abab(e,f,i,n) 
            // flops: o2v2L1: 1, o1v2L1: 1, o1v1L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= F_blks_["aa_ov"]("i,a") * r1_xa_Lv("I,a") * u2_abab_vvoo("e,f,i,n");

            // sigmas2_xabb_Lvvo += -1.000000 d+_bb(f,n) r1_a(e) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= dp_bb_vo("f,n") * r1_xa_Lv("I,e");

            // sigmas2_xabb_Lvvo += -1.000000 d+_aa(i,a) r1_a(e) t2_abab(a,f,i,n) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= dp_aa_ov("i,a") * r1_xa_Lv("I,e") * t2_abab_vvoo("a,f,i,n");

            // sigmas2_xabb_Lvvo += 1.000000 d+_bb(i,a) r1_a(e) t2_bbbb(a,f,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            sigmas2_xabb_Lvvo("I,e,f,n") += dp_bb_ov("i,a") * r1_xa_Lv("I,e") * t2_bbbb_vvoo("a,f,n,i");

            // sigmas2_xabb_Lvvo += 1.000000 d+_aa(i,a) r1_a(a) t2_abab(e,f,i,n) 
            // flops: o2v2L1: 1, o1v2L1: 1, o1v1L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += dp_aa_ov("i,a") * r1_xa_Lv("I,a") * t2_abab_vvoo("e,f,i,n");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmas2_xabb_Lvvo += -1.000000 d-_bb(f,n) u0 s1_a(e) 
            // flops: o1v2L1: 2, o1v1: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= dp_bb_vo("f,n") * u0 * s1_xa_Lv("I,e");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xabb_Lvvo += -1.000000 d-_aa(i,i) u1_bb(f,n) s1_a(e) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= scalar7 * s1_xa_Lv("I,e") * u1_bb_vo("f,n");

            // sigmas2_xabb_Lvvo += -1.000000 d-_bb(i,i) u1_bb(f,n) s1_a(e) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= scalar8 * s1_xa_Lv("I,e") * u1_bb_vo("f,n");

            // sigmas2_xabb_Lvvo += 2.000000 d-_bb(i,n) u1_bb(f,i) s1_a(e) 
            // flops: o1v2L1: 2, o2v1: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += 2.000000 * dp_bb_oo("i,n") * u1_bb_vo("f,i") * s1_xa_Lv("I,e");

            // sigmas2_xabb_Lvvo += -2.000000 d-_bb(f,a) u1_bb(a,n) s1_a(e) 
            // flops: o1v2L1: 2, o1v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= 2.000000 * dp_bb_vv("f,a") * u1_bb_vo("a,n") * s1_xa_Lv("I,e");

            // sigmas2_xabb_Lvvo += -1.000000 d-_aa(e,a) u1_bb(f,n) s1_a(a) 
            // flops: o1v3L1: 1, o1v2L1: 1, o1v3: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= dp_aa_vv("e,a") * u1_bb_vo("f,n") * s1_xa_Lv("I,a");

            // sigmas2_xabb_Lvvo += 1.000000 d-_aa(i,a) u1_bb(f,n) s2_aaa(a,e,i) 
            // flops: o2v3L1: 1, o1v2L1: 1, o2v2: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += dp_aa_ov("i,a") * u1_bb_vo("f,n") * s2_xaaa_Lvvo("I,a,e,i");

            // sigmas2_xabb_Lvvo += -2.000000 d-_aa(i,a) u2_abab(a,f,i,n) s1_a(e) 
            // flops: o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= 2.000000 * dp_aa_ov("i,a") * u2_abab_vvoo("a,f,i,n") * s1_xa_Lv("I,e");

            // sigmas2_xabb_Lvvo += 2.000000 d-_bb(i,a) u2_bbbb(a,f,n,i) s1_a(e) 
            // flops: o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += 2.000000 * dp_bb_ov("i,a") * u2_bbbb_vvoo("a,f,n,i") * s1_xa_Lv("I,e");

            // sigmas2_xabb_Lvvo += 2.000000 d-_aa(i,a) u2_abab(e,f,i,n) s1_a(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += 2.000000 * dp_aa_ov("i,a") * u2_abab_vvoo("e,f,i,n") * s1_xa_Lv("I,a");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmas2_xabb_Lvvo += -1.000000 d-_aa(i,a) t2_abab(a,f,i,n) u0 s1_a(e) 
            // flops: o1v2L1: 2, o2v2: 1, o1v1: 1 | mem: o1v2L1: 2, o1v1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= dp_aa_ov("i,a") * t2_abab_vvoo("a,f,i,n") * u0 * s1_xa_Lv("I,e");

            // sigmas2_xabb_Lvvo += 1.000000 d-_bb(i,a) t2_bbbb(a,f,n,i) u0 s1_a(e) 
            // flops: o1v2L1: 2, o2v2: 1, o1v1: 1 | mem: o1v2L1: 2, o1v1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") += dp_bb_ov("i,a") * t2_bbbb_vvoo("a,f,n,i") * u0 * s1_xa_Lv("I,e");

            // sigmas2_xabb_Lvvo += 1.000000 d-_aa(i,a) t2_abab(e,f,i,n) u0 s1_a(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1, o1v3: 1 | mem: o1v2L1: 2, o1v3: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") += dp_aa_ov("i,a") * t2_abab_vvoo("e,f,i,n") * u0 * s1_xa_Lv("I,a");

            // sigmas2_xabb_Lvvo += 1.000000 d-_bb(i,n) r1_a(e) u0 u1_bb(f,i) 
            // flops: o1v2L1: 1, o2v1L1: 2, o2v2: 1 | mem: o1v2L1: 1, o2v1L1: 1, o1v2: 1, o2v1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += dp_bb_oo("i,n") * r1_xa_Lv("I,e") * u0 * u1_bb_vo("f,i");

            // sigmas2_xabb_Lvvo += -1.000000 d-_bb(f,a) r1_a(e) u0 u1_bb(a,n) 
            // flops: o0v3L1: 2, o1v2L1: 1, o1v3: 1 | mem: o0v3L1: 1, o1v2L1: 1, o0v3: 1, o1v2: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= dp_bb_vv("f,a") * r1_xa_Lv("I,e") * u0 * u1_bb_vo("a,n");
        }

        if (include_u0_ && include_u2_) {

            // sigmas2_xabb_Lvvo += -1.000000 d-_aa(i,a) r1_a(e) u0 u2_abab(a,f,i,n) 
            // flops: o2v3: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v2: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= dp_aa_ov("i,a") * r1_xa_Lv("I,e") * u0 * u2_abab_vvoo("a,f,i,n");

            // sigmas2_xabb_Lvvo += 1.000000 d-_bb(i,a) r1_a(e) u0 u2_bbbb(a,f,n,i) 
            // flops: o2v3: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v2: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") += dp_bb_ov("i,a") * r1_xa_Lv("I,e") * u0 * u2_bbbb_vvoo("a,f,n,i");

            // sigmas2_xabb_Lvvo += 1.000000 d-_aa(i,a) r1_a(a) u0 u2_abab(e,f,i,n) 
            // flops: o1v2L1: 1, o2v2: 1, o1v1L1: 1, o1v0L1: 1 | mem: o1v2L1: 1, o1v2: 1, o1v0L1: 1, o1v0: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += dp_aa_ov("i,a") * r1_xa_Lv("I,a") * u0 * u2_abab_vvoo("e,f,i,n");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xabb_Lvvo += -1.000000 d-_aa(i,a) r1_a(e) u1_aa(a,i) u1_bb(f,n) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= scalar0 * r1_xa_Lv("I,e") * u1_bb_vo("f,n");

            // sigmas2_xabb_Lvvo += -1.000000 d-_bb(i,a) r1_a(e) u1_bb(a,i) u1_bb(f,n) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= scalar1 * r1_xa_Lv("I,e") * u1_bb_vo("f,n");

            // sigmas2_xabb_Lvvo += 2.000000 d-_bb(i,a) r1_a(e) u1_bb(a,n) u1_bb(f,i) 
            // flops: o2v2L1: 2, o1v2L1: 2 | mem: o1v2L1: 3, o2v1L1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += 2.000000 * dp_bb_ov("i,a") * r1_xa_Lv("I,e") * u1_bb_vo("a,n") * u1_bb_vo("f,i");

            // sigmas2_xabb_Lvvo += 1.000000 d-_aa(i,a) r1_a(a) u1_aa(e,i) u1_bb(f,n) 
            // flops: o1v2L1: 2, o1v1L1: 2 | mem: o1v2L1: 2, o0v1L1: 1, o1v0L1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += dp_aa_ov("i,a") * r1_xa_Lv("I,a") * u1_aa_vo("e,i") * u1_bb_vo("f,n");

            // sigmas2_xabb_Lvvo += 1.000000 <e,f||a,n>_abab s1_a(a) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") += V_blks_["abab_vvvo"]("e,f,a,n") * s1_xa_Lv("I,a");
        }

        if (include_u2_) {

            // sigmas2_xabb_Lvvo += -1.000000 <i,f||a,n>_abab s2_aaa(a,e,i) 
            // flops: o2v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") += V_blks_["baab_vovo"]("f,i,a,n") * s2_xaaa_Lvvo("I,a,e,i");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xabb_Lvvo += -0.500000 <j,i||a,n>_abab t2_abab(a,f,j,i) s1_a(e) 
            // flops: o3v2: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovo"]("j,i,a,n") * t2_abab_vvoo("a,f,j,i") * s1_xa_Lv("I,e");

            // sigmas2_xabb_Lvvo += -0.500000 <i,j||a,n>_abab t2_abab(a,f,i,j) s1_a(e) 
            // flops: o3v2: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovo"]("i,j,a,n") * t2_abab_vvoo("a,f,i,j") * s1_xa_Lv("I,e");

            // sigmas2_xabb_Lvvo += -0.500000 <j,i||a,n>_bbbb t2_bbbb(a,f,j,i) s1_a(e) 
            // flops: o3v2: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovo"]("j,i,a,n") * t2_bbbb_vvoo("a,f,j,i") * s1_xa_Lv("I,e");

            // sigmas2_xabb_Lvvo += 0.500000 <j,i||a,n>_abab t2_abab(e,f,j,i) s1_a(a) 
            // flops: o3v3: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovo"]("j,i,a,n") * t2_abab_vvoo("e,f,j,i") * s1_xa_Lv("I,a");

            // sigmas2_xabb_Lvvo += 0.500000 <i,j||a,n>_abab t2_abab(e,f,i,j) s1_a(a) 
            // flops: o3v3: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovo"]("i,j,a,n") * t2_abab_vvoo("e,f,i,j") * s1_xa_Lv("I,a");

            // sigmas2_xabb_Lvvo += 0.500000 <i,f||a,b>_abab t2_abab(a,b,i,n) s1_a(e) 
            // flops: o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["baab_vovv"]("f,i,a,b") * t2_abab_vvoo("a,b,i,n") * s1_xa_Lv("I,e");

            // sigmas2_xabb_Lvvo += 0.500000 <i,f||b,a>_abab t2_abab(b,a,i,n) s1_a(e) 
            // flops: o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["baab_vovv"]("f,i,b,a") * t2_abab_vvoo("b,a,i,n") * s1_xa_Lv("I,e");

            // sigmas2_xabb_Lvvo += -0.500000 <i,f||a,b>_bbbb t2_bbbb(a,b,n,i) s1_a(e) 
            // flops: o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["bbbb_vovv"]("f,i,a,b") * t2_bbbb_vvoo("a,b,n,i") * s1_xa_Lv("I,e");

            // sigmas2_xabb_Lvvo += 1.000000 <i,e||a,b>_aaaa t2_abab(a,f,i,n) s1_a(b) 
            // flops: o2v4: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= V_blks_["aaaa_vovv"]("e,i,a,b") * t2_abab_vvoo("a,f,i,n") * s1_xa_Lv("I,b");

            // sigmas2_xabb_Lvvo += -1.000000 <e,i||b,a>_abab t2_bbbb(a,f,n,i) s1_a(b) 
            // flops: o2v4: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= V_blks_["abab_vovv"]("e,i,b,a") * t2_bbbb_vvoo("a,f,n,i") * s1_xa_Lv("I,b");

            // sigmas2_xabb_Lvvo += -1.000000 <i,f||b,a>_abab t2_abab(e,a,i,n) s1_a(b) 
            // flops: o2v4: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += V_blks_["baab_vovv"]("f,i,b,a") * t2_abab_vvoo("e,a,i,n") * s1_xa_Lv("I,b");
        }

        if (include_u2_) {

            // sigmas2_xabb_Lvvo += 1.000000 <j,i||a,b>_aaaa t2_abab(a,f,i,n) s2_aaa(b,e,j) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += V_blks_["aaaa_oovv"]("j,i,a,b") * t2_abab_vvoo("a,f,i,n") * s2_xaaa_Lvvo("I,b,e,j");

            // sigmas2_xabb_Lvvo += 1.000000 <j,i||b,a>_abab t2_bbbb(a,f,n,i) s2_aaa(b,e,j) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("j,i,b,a") * t2_bbbb_vvoo("a,f,n,i") * s2_xaaa_Lvvo("I,b,e,j");

            // sigmas2_xabb_Lvvo += 0.500000 <j,i||a,b>_aaaa t2_abab(e,f,i,n) s2_aaa(a,b,j) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * t2_abab_vvoo("e,f,i,n") * s2_xaaa_Lvvo("I,a,b,j");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xabb_Lvvo += 1.000000 <i,f||a,n>_abab r1_a(e) u1_aa(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= V_blks_["baab_vovo"]("f,i,a,n") * r1_xa_Lv("I,e") * u1_aa_vo("a,i");

            // sigmas2_xabb_Lvvo += 1.000000 <i,f||a,n>_bbbb r1_a(e) u1_bb(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= V_blks_["bbbb_vovo"]("f,i,a,n") * r1_xa_Lv("I,e") * u1_bb_vo("a,i");

            // sigmas2_xabb_Lvvo += -1.000000 <e,i||a,n>_abab r1_a(a) u1_bb(f,i) 
            // flops: o2v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= V_blks_["abab_vovo"]("e,i,a,n") * r1_xa_Lv("I,a") * u1_bb_vo("f,i");

            // sigmas2_xabb_Lvvo += -1.000000 <i,f||a,n>_abab r1_a(a) u1_aa(e,i) 
            // flops: o2v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += V_blks_["baab_vovo"]("f,i,a,n") * r1_xa_Lv("I,a") * u1_aa_vo("e,i");

            // sigmas2_xabb_Lvvo += 1.000000 <e,f||b,a>_abab r1_a(b) u1_bb(a,n) 
            // flops: o0v4L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") += V_blks_["abab_vvvv"]("e,f,b,a") * r1_xa_Lv("I,b") * u1_bb_vo("a,n");

            // sigmas2_xabb_Lvvo += 1.000000 <j,i||a,n>_abab r2_aaa(a,e,j) u1_bb(f,i) 
            // flops: o3v2L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += V_blks_["abab_oovo"]("j,i,a,n") * r2_xaaa_Lvvo("I,a,e,j") * u1_bb_vo("f,i");

            // sigmas2_xabb_Lvvo += -1.000000 <i,f||b,a>_abab r2_aaa(b,e,i) u1_bb(a,n) 
            // flops: o1v4L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") += V_blks_["baab_vovv"]("f,i,b,a") * r2_xaaa_Lvvo("I,b,e,i") * u1_bb_vo("a,n");
        }

        if (include_u2_) {

            // sigmas2_xabb_Lvvo += -0.500000 <j,i||a,n>_abab r1_a(e) u2_abab(a,f,j,i) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 1 | mem: o3v2L1: 1, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovo"]("j,i,a,n") * r1_xa_Lv("I,e") * u2_abab_vvoo("a,f,j,i");

            // sigmas2_xabb_Lvvo += -0.500000 <i,j||a,n>_abab r1_a(e) u2_abab(a,f,i,j) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 1 | mem: o3v2L1: 1, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovo"]("i,j,a,n") * r1_xa_Lv("I,e") * u2_abab_vvoo("a,f,i,j");

            // sigmas2_xabb_Lvvo += -0.500000 <j,i||a,n>_bbbb r1_a(e) u2_bbbb(a,f,j,i) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 1 | mem: o3v2L1: 1, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovo"]("j,i,a,n") * r1_xa_Lv("I,e") * u2_bbbb_vvoo("a,f,j,i");

            // sigmas2_xabb_Lvvo += 0.500000 <j,i||a,n>_abab r1_a(a) u2_abab(e,f,j,i) 
            // flops: o3v2L1: 1, o3v1L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o3v0L1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovo"]("j,i,a,n") * r1_xa_Lv("I,a") * u2_abab_vvoo("e,f,j,i");

            // sigmas2_xabb_Lvvo += 0.500000 <i,j||a,n>_abab r1_a(a) u2_abab(e,f,i,j) 
            // flops: o3v2L1: 1, o3v1L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o3v0L1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovo"]("i,j,a,n") * r1_xa_Lv("I,a") * u2_abab_vvoo("e,f,i,j");

            // sigmas2_xabb_Lvvo += 0.500000 <i,f||a,b>_abab r1_a(e) u2_abab(a,b,i,n) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["baab_vovv"]("f,i,a,b") * r1_xa_Lv("I,e") * u2_abab_vvoo("a,b,i,n");

            // sigmas2_xabb_Lvvo += 0.500000 <i,f||b,a>_abab r1_a(e) u2_abab(b,a,i,n) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["baab_vovv"]("f,i,b,a") * r1_xa_Lv("I,e") * u2_abab_vvoo("b,a,i,n");

            // sigmas2_xabb_Lvvo += -0.500000 <i,f||a,b>_bbbb r1_a(e) u2_bbbb(a,b,n,i) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["bbbb_vovv"]("f,i,a,b") * r1_xa_Lv("I,e") * u2_bbbb_vvoo("a,b,n,i");

            // sigmas2_xabb_Lvvo += 1.000000 <i,e||a,b>_aaaa r1_a(b) u2_abab(a,f,i,n) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= V_blks_["aaaa_vovv"]("e,i,a,b") * r1_xa_Lv("I,b") * u2_abab_vvoo("a,f,i,n");

            // sigmas2_xabb_Lvvo += -1.000000 <e,i||b,a>_abab r1_a(b) u2_bbbb(a,f,n,i) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= V_blks_["abab_vovv"]("e,i,b,a") * r1_xa_Lv("I,b") * u2_bbbb_vvoo("a,f,n,i");

            // sigmas2_xabb_Lvvo += -1.000000 <i,f||b,a>_abab r1_a(b) u2_abab(e,a,i,n) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmas2_xabb_Lvvo("I,e,f,n") += V_blks_["baab_vovv"]("f,i,b,a") * r1_xa_Lv("I,b") * u2_abab_vvoo("e,a,i,n");

            // sigmas2_xabb_Lvvo += 1.000000 <j,i||a,b>_aaaa r2_aaa(b,e,j) u2_abab(a,f,i,n) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmas2_xabb_Lvvo("I,e,f,n") += V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaa_Lvvo("I,b,e,j") * u2_abab_vvoo("a,f,i,n");

            // sigmas2_xabb_Lvvo += 1.000000 <j,i||b,a>_abab r2_aaa(b,e,j) u2_bbbb(a,f,n,i) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmas2_xabb_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("j,i,b,a") * r2_xaaa_Lvvo("I,b,e,j") * u2_bbbb_vvoo("a,f,n,i");

            // sigmas2_xabb_Lvvo += 0.500000 <j,i||a,b>_aaaa r2_aaa(a,b,j) u2_abab(e,f,i,n) 
            // flops: o2v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaa_Lvvo("I,a,b,j") * u2_abab_vvoo("e,f,i,n");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xabb_Lvvo += -0.500000 <i,j||a,b>_abab r1_a(e) t2_abab(a,b,i,n) u1_bb(f,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("i,j,a,b") * r1_xa_Lv("I,e") * t2_abab_vvoo("a,b,i,n") * u1_bb_vo("f,j");

            // sigmas2_xabb_Lvvo += -0.500000 <i,j||b,a>_abab r1_a(e) t2_abab(b,a,i,n) u1_bb(f,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("i,j,b,a") * r1_xa_Lv("I,e") * t2_abab_vvoo("b,a,i,n") * u1_bb_vo("f,j");

            // sigmas2_xabb_Lvvo += -0.500000 <j,i||a,b>_bbbb r1_a(e) t2_bbbb(a,b,n,i) u1_bb(f,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r1_xa_Lv("I,e") * t2_bbbb_vvoo("a,b,n,i") * u1_bb_vo("f,j");

            // sigmas2_xabb_Lvvo += -0.500000 <j,i||a,b>_abab r1_a(e) t2_abab(a,f,j,i) u1_bb(b,n) 
            // flops: o2v4L1: 1, o2v3L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o0v3L1: 1, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("j,i,a,b") * r1_xa_Lv("I,e") * t2_abab_vvoo("a,f,j,i") * u1_bb_vo("b,n");

            // sigmas2_xabb_Lvvo += -0.500000 <i,j||a,b>_abab r1_a(e) t2_abab(a,f,i,j) u1_bb(b,n) 
            // flops: o2v4L1: 1, o2v3L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o0v3L1: 1, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("i,j,a,b") * r1_xa_Lv("I,e") * t2_abab_vvoo("a,f,i,j") * u1_bb_vo("b,n");

            // sigmas2_xabb_Lvvo += -0.500000 <j,i||a,b>_bbbb r1_a(e) t2_bbbb(a,f,j,i) u1_bb(b,n) 
            // flops: o2v4L1: 1, o2v3L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o0v3L1: 1, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r1_xa_Lv("I,e") * t2_bbbb_vvoo("a,f,j,i") * u1_bb_vo("b,n");

            // sigmas2_xabb_Lvvo += -1.000000 <j,i||a,b>_aaaa r1_a(e) t2_abab(a,f,i,n) u1_aa(b,j) 
            // flops: o3v4L1: 1, o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 2, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= V_blks_["aaaa_oovv"]("j,i,a,b") * r1_xa_Lv("I,e") * t2_abab_vvoo("a,f,i,n") * u1_aa_vo("b,j");

            // sigmas2_xabb_Lvvo += 1.000000 <i,j||a,b>_abab r1_a(e) t2_abab(a,f,i,n) u1_bb(b,j) 
            // flops: o3v4L1: 1, o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 2, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("i,j,a,b") * r1_xa_Lv("I,e") * t2_abab_vvoo("a,f,i,n") * u1_bb_vo("b,j");

            // sigmas2_xabb_Lvvo += -1.000000 <j,i||b,a>_abab r1_a(e) t2_bbbb(a,f,n,i) u1_aa(b,j) 
            // flops: o3v4L1: 1, o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 2, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= V_blks_["abab_oovv"]("j,i,b,a") * r1_xa_Lv("I,e") * t2_bbbb_vvoo("a,f,n,i") * u1_aa_vo("b,j");

            // sigmas2_xabb_Lvvo += 1.000000 <j,i||a,b>_bbbb r1_a(e) t2_bbbb(a,f,n,i) u1_bb(b,j) 
            // flops: o3v4L1: 1, o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 2, o1v2L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") += V_blks_["bbbb_oovv"]("j,i,a,b") * r1_xa_Lv("I,e") * t2_bbbb_vvoo("a,f,n,i") * u1_bb_vo("b,j");

            // sigmas2_xabb_Lvvo += 1.000000 <i,j||b,a>_abab r1_a(b) t2_abab(e,a,i,n) u1_bb(f,j) 
            // flops: o3v2L1: 1, o2v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o2v1L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("i,j,b,a") * r1_xa_Lv("I,b") * t2_abab_vvoo("e,a,i,n") * u1_bb_vo("f,j");

            // sigmas2_xabb_Lvvo += 1.000000 <j,i||a,b>_aaaa r1_a(b) t2_abab(a,f,i,n) u1_aa(e,j) 
            // flops: o3v2L1: 1, o2v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o2v1L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") += V_blks_["aaaa_oovv"]("j,i,a,b") * r1_xa_Lv("I,b") * t2_abab_vvoo("a,f,i,n") * u1_aa_vo("e,j");

            // sigmas2_xabb_Lvvo += 1.000000 <j,i||b,a>_abab r1_a(b) t2_bbbb(a,f,n,i) u1_aa(e,j) 
            // flops: o3v2L1: 1, o2v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o2v1L1: 2, 
            sigmas2_xabb_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("j,i,b,a") * r1_xa_Lv("I,b") * t2_bbbb_vvoo("a,f,n,i") * u1_aa_vo("e,j");

            // sigmas2_xabb_Lvvo += 0.500000 <j,i||b,a>_abab r1_a(b) t2_abab(e,f,j,i) u1_bb(a,n) 
            // flops: o2v3L1: 1, o1v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("j,i,b,a") * r1_xa_Lv("I,b") * t2_abab_vvoo("e,f,j,i") * u1_bb_vo("a,n");

            // sigmas2_xabb_Lvvo += 0.500000 <i,j||b,a>_abab r1_a(b) t2_abab(e,f,i,j) u1_bb(a,n) 
            // flops: o2v3L1: 1, o1v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("i,j,b,a") * r1_xa_Lv("I,b") * t2_abab_vvoo("e,f,i,j") * u1_bb_vo("a,n");

            // sigmas2_xabb_Lvvo += -1.000000 <j,i||a,b>_aaaa r1_a(b) t2_abab(e,f,i,n) u1_aa(a,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= V_blks_["aaaa_oovv"]("j,i,a,b") * r1_xa_Lv("I,b") * t2_abab_vvoo("e,f,i,n") * u1_aa_vo("a,j");

            // sigmas2_xabb_Lvvo += -1.000000 <i,j||b,a>_abab r1_a(b) t2_abab(e,f,i,n) u1_bb(a,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xabb_Lvvo("I,e,f,n") -= V_blks_["abab_oovv"]("i,j,b,a") * r1_xa_Lv("I,b") * t2_abab_vvoo("e,f,i,n") * u1_bb_vo("a,j");
        }




/// ****** sigmas2_baa(e,f,n) ****** ///



        if (include_u1_ && include_u2_) {

            // sigmas2_xbaa_Lvvo += 1.000000 r1_b(e) u1_aa(f,n) w0 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, o1v2: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += r1_xb_Lv("I,e") * u1_aa_vo("f,n") * w0;

            // sigmas2_xbaa_Lvvo += 1.000000 f_aa(f,n) s1_b(e) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += F_blks_["aa_vo"]("f,n") * s1_xb_Lv("I,e");

            // sigmas2_xbaa_Lvvo += -1.000000 f_aa(i,a) t2_aaaa(a,f,n,i) s1_b(e) 
            // flops: o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= F_blks_["aa_ov"]("i,a") * t2_aaaa_vvoo("a,f,n,i") * s1_xb_Lv("I,e");

            // sigmas2_xbaa_Lvvo += 1.000000 f_bb(i,a) t2_abab(f,a,n,i) s1_b(e) 
            // flops: o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += F_blks_["bb_ov"]("i,a") * t2_abab_vvoo("f,a,n,i") * s1_xb_Lv("I,e");

            // sigmas2_xbaa_Lvvo += -1.000000 f_bb(i,a) t2_abab(f,e,n,i) s1_b(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= F_blks_["bb_ov"]("i,a") * t2_abab_vvoo("f,e,n,i") * s1_xb_Lv("I,a");

            // sigmas2_xbaa_Lvvo += -1.000000 f_aa(i,n) r1_b(e) u1_aa(f,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o2v1L1: 1 | mem: o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= F_blks_["aa_oo"]("i,n") * r1_xb_Lv("I,e") * u1_aa_vo("f,i");

            // sigmas2_xbaa_Lvvo += 1.000000 f_aa(f,a) r1_b(e) u1_aa(a,n) 
            // flops: o1v3L1: 1, o0v3L1: 1, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += F_blks_["aa_vv"]("f,a") * r1_xb_Lv("I,e") * u1_aa_vo("a,n");
        }

        if (include_u2_) {

            // sigmas2_xbaa_Lvvo += -1.000000 f_aa(i,a) r1_b(e) u2_aaaa(a,f,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= F_blks_["aa_ov"]("i,a") * r1_xb_Lv("I,e") * u2_aaaa_vvoo("a,f,n,i");

            // sigmas2_xbaa_Lvvo += 1.000000 f_bb(i,a) r1_b(e) u2_abab(f,a,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += F_blks_["bb_ov"]("i,a") * r1_xb_Lv("I,e") * u2_abab_vvoo("f,a,n,i");

            // sigmas2_xbaa_Lvvo += -1.000000 f_bb(i,a) r1_b(a) u2_abab(f,e,n,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o1v1L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= F_blks_["bb_ov"]("i,a") * r1_xb_Lv("I,a") * u2_abab_vvoo("f,e,n,i");

            // sigmas2_xbaa_Lvvo += -1.000000 d+_aa(f,n) r1_b(e) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= dp_aa_vo("f,n") * r1_xb_Lv("I,e");

            // sigmas2_xbaa_Lvvo += 1.000000 d+_aa(i,a) r1_b(e) t2_aaaa(a,f,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += dp_aa_ov("i,a") * r1_xb_Lv("I,e") * t2_aaaa_vvoo("a,f,n,i");

            // sigmas2_xbaa_Lvvo += -1.000000 d+_bb(i,a) r1_b(e) t2_abab(f,a,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= dp_bb_ov("i,a") * r1_xb_Lv("I,e") * t2_abab_vvoo("f,a,n,i");

            // sigmas2_xbaa_Lvvo += 1.000000 d+_bb(i,a) r1_b(a) t2_abab(f,e,n,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o1v1L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += dp_bb_ov("i,a") * r1_xb_Lv("I,a") * t2_abab_vvoo("f,e,n,i");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmas2_xbaa_Lvvo += -1.000000 d-_aa(f,n) u0 s1_b(e) 
            // flops: o1v2L1: 2, o1v1: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= dp_aa_vo("f,n") * u0 * s1_xb_Lv("I,e");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xbaa_Lvvo += -1.000000 d-_aa(i,i) u1_aa(f,n) s1_b(e) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= scalar7 * s1_xb_Lv("I,e") * u1_aa_vo("f,n");

            // sigmas2_xbaa_Lvvo += -1.000000 d-_bb(i,i) u1_aa(f,n) s1_b(e) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= scalar8 * s1_xb_Lv("I,e") * u1_aa_vo("f,n");

            // sigmas2_xbaa_Lvvo += 2.000000 d-_aa(i,n) u1_aa(f,i) s1_b(e) 
            // flops: o1v2L1: 2, o2v1: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += 2.000000 * dp_aa_oo("i,n") * u1_aa_vo("f,i") * s1_xb_Lv("I,e");

            // sigmas2_xbaa_Lvvo += -2.000000 d-_aa(f,a) u1_aa(a,n) s1_b(e) 
            // flops: o1v2L1: 2, o1v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= 2.000000 * dp_aa_vv("f,a") * u1_aa_vo("a,n") * s1_xb_Lv("I,e");

            // sigmas2_xbaa_Lvvo += -1.000000 d-_bb(e,a) u1_aa(f,n) s1_b(a) 
            // flops: o1v3L1: 1, o1v2L1: 1, o1v3: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= dp_bb_vv("e,a") * u1_aa_vo("f,n") * s1_xb_Lv("I,a");

            // sigmas2_xbaa_Lvvo += 1.000000 d-_bb(i,a) u1_aa(f,n) s2_bbb(a,e,i) 
            // flops: o2v3L1: 1, o1v2L1: 1, o2v2: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += dp_bb_ov("i,a") * u1_aa_vo("f,n") * s2_xbbb_Lvvo("I,a,e,i");

            // sigmas2_xbaa_Lvvo += 2.000000 d-_aa(i,a) u2_aaaa(a,f,n,i) s1_b(e) 
            // flops: o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += 2.000000 * dp_aa_ov("i,a") * u2_aaaa_vvoo("a,f,n,i") * s1_xb_Lv("I,e");

            // sigmas2_xbaa_Lvvo += -2.000000 d-_bb(i,a) u2_abab(f,a,n,i) s1_b(e) 
            // flops: o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= 2.000000 * dp_bb_ov("i,a") * u2_abab_vvoo("f,a,n,i") * s1_xb_Lv("I,e");

            // sigmas2_xbaa_Lvvo += 2.000000 d-_bb(i,a) u2_abab(f,e,n,i) s1_b(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += 2.000000 * dp_bb_ov("i,a") * u2_abab_vvoo("f,e,n,i") * s1_xb_Lv("I,a");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmas2_xbaa_Lvvo += 1.000000 d-_aa(i,a) t2_aaaa(a,f,n,i) u0 s1_b(e) 
            // flops: o1v2L1: 2, o2v2: 1, o1v1: 1 | mem: o1v2L1: 2, o1v1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += dp_aa_ov("i,a") * t2_aaaa_vvoo("a,f,n,i") * u0 * s1_xb_Lv("I,e");

            // sigmas2_xbaa_Lvvo += -1.000000 d-_bb(i,a) t2_abab(f,a,n,i) u0 s1_b(e) 
            // flops: o1v2L1: 2, o2v2: 1, o1v1: 1 | mem: o1v2L1: 2, o1v1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= dp_bb_ov("i,a") * t2_abab_vvoo("f,a,n,i") * u0 * s1_xb_Lv("I,e");

            // sigmas2_xbaa_Lvvo += 1.000000 d-_bb(i,a) t2_abab(f,e,n,i) u0 s1_b(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1, o1v3: 1 | mem: o1v2L1: 2, o1v3: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += dp_bb_ov("i,a") * t2_abab_vvoo("f,e,n,i") * u0 * s1_xb_Lv("I,a");

            // sigmas2_xbaa_Lvvo += 1.000000 d-_aa(i,n) r1_b(e) u0 u1_aa(f,i) 
            // flops: o1v2L1: 1, o2v1L1: 2, o2v2: 1 | mem: o1v2L1: 1, o2v1L1: 1, o1v2: 1, o2v1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += dp_aa_oo("i,n") * r1_xb_Lv("I,e") * u0 * u1_aa_vo("f,i");

            // sigmas2_xbaa_Lvvo += -1.000000 d-_aa(f,a) r1_b(e) u0 u1_aa(a,n) 
            // flops: o0v3L1: 2, o1v2L1: 1, o1v3: 1 | mem: o0v3L1: 1, o1v2L1: 1, o0v3: 1, o1v2: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= dp_aa_vv("f,a") * r1_xb_Lv("I,e") * u0 * u1_aa_vo("a,n");
        }

        if (include_u0_ && include_u2_) {

            // sigmas2_xbaa_Lvvo += 1.000000 d-_aa(i,a) r1_b(e) u0 u2_aaaa(a,f,n,i) 
            // flops: o2v3: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v2: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += dp_aa_ov("i,a") * r1_xb_Lv("I,e") * u0 * u2_aaaa_vvoo("a,f,n,i");

            // sigmas2_xbaa_Lvvo += -1.000000 d-_bb(i,a) r1_b(e) u0 u2_abab(f,a,n,i) 
            // flops: o2v3: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v2: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= dp_bb_ov("i,a") * r1_xb_Lv("I,e") * u0 * u2_abab_vvoo("f,a,n,i");

            // sigmas2_xbaa_Lvvo += 1.000000 d-_bb(i,a) r1_b(a) u0 u2_abab(f,e,n,i) 
            // flops: o1v2L1: 1, o2v2: 1, o1v1L1: 1, o1v0L1: 1 | mem: o1v2L1: 1, o1v2: 1, o1v0L1: 1, o1v0: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += dp_bb_ov("i,a") * r1_xb_Lv("I,a") * u0 * u2_abab_vvoo("f,e,n,i");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xbaa_Lvvo += -1.000000 d-_aa(i,a) r1_b(e) u1_aa(a,i) u1_aa(f,n) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= scalar0 * r1_xb_Lv("I,e") * u1_aa_vo("f,n");

            // sigmas2_xbaa_Lvvo += -1.000000 d-_bb(i,a) r1_b(e) u1_bb(a,i) u1_aa(f,n) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= scalar1 * r1_xb_Lv("I,e") * u1_aa_vo("f,n");

            // sigmas2_xbaa_Lvvo += 2.000000 d-_aa(i,a) r1_b(e) u1_aa(a,n) u1_aa(f,i) 
            // flops: o2v2L1: 2, o1v2L1: 2 | mem: o1v2L1: 3, o2v1L1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += 2.000000 * dp_aa_ov("i,a") * r1_xb_Lv("I,e") * u1_aa_vo("a,n") * u1_aa_vo("f,i");

            // sigmas2_xbaa_Lvvo += 1.000000 d-_bb(i,a) r1_b(a) u1_bb(e,i) u1_aa(f,n) 
            // flops: o1v2L1: 2, o1v1L1: 2 | mem: o1v2L1: 2, o0v1L1: 1, o1v0L1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += dp_bb_ov("i,a") * r1_xb_Lv("I,a") * u1_bb_vo("e,i") * u1_aa_vo("f,n");

            // sigmas2_xbaa_Lvvo += 1.000000 <f,e||n,a>_abab s1_b(a) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= V_blks_["abba_vvvo"]("f,e,a,n") * s1_xb_Lv("I,a");
        }

        if (include_u2_) {

            // sigmas2_xbaa_Lvvo += -1.000000 <f,i||n,a>_abab s2_bbb(a,e,i) 
            // flops: o2v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += V_blks_["abba_vovo"]("f,i,a,n") * s2_xbbb_Lvvo("I,a,e,i");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xbaa_Lvvo += -0.500000 <j,i||a,n>_aaaa t2_aaaa(a,f,j,i) s1_b(e) 
            // flops: o3v2: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovo"]("j,i,a,n") * t2_aaaa_vvoo("a,f,j,i") * s1_xb_Lv("I,e");

            // sigmas2_xbaa_Lvvo += -0.500000 <j,i||n,a>_abab t2_abab(f,a,j,i) s1_b(e) 
            // flops: o3v2: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abba_oovo"]("j,i,a,n") * t2_abab_vvoo("f,a,j,i") * s1_xb_Lv("I,e");

            // sigmas2_xbaa_Lvvo += -0.500000 <i,j||n,a>_abab t2_abab(f,a,i,j) s1_b(e) 
            // flops: o3v2: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abba_oovo"]("i,j,a,n") * t2_abab_vvoo("f,a,i,j") * s1_xb_Lv("I,e");

            // sigmas2_xbaa_Lvvo += 0.500000 <j,i||n,a>_abab t2_abab(f,e,j,i) s1_b(a) 
            // flops: o3v3: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abba_oovo"]("j,i,a,n") * t2_abab_vvoo("f,e,j,i") * s1_xb_Lv("I,a");

            // sigmas2_xbaa_Lvvo += 0.500000 <i,j||n,a>_abab t2_abab(f,e,i,j) s1_b(a) 
            // flops: o3v3: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abba_oovo"]("i,j,a,n") * t2_abab_vvoo("f,e,i,j") * s1_xb_Lv("I,a");

            // sigmas2_xbaa_Lvvo += -0.500000 <i,f||a,b>_aaaa t2_aaaa(a,b,n,i) s1_b(e) 
            // flops: o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["aaaa_vovv"]("f,i,a,b") * t2_aaaa_vvoo("a,b,n,i") * s1_xb_Lv("I,e");

            // sigmas2_xbaa_Lvvo += 0.500000 <f,i||a,b>_abab t2_abab(a,b,n,i) s1_b(e) 
            // flops: o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_vovv"]("f,i,a,b") * t2_abab_vvoo("a,b,n,i") * s1_xb_Lv("I,e");

            // sigmas2_xbaa_Lvvo += 0.500000 <f,i||b,a>_abab t2_abab(b,a,n,i) s1_b(e) 
            // flops: o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_vovv"]("f,i,b,a") * t2_abab_vvoo("b,a,n,i") * s1_xb_Lv("I,e");

            // sigmas2_xbaa_Lvvo += -1.000000 <i,e||a,b>_abab t2_aaaa(a,f,n,i) s1_b(b) 
            // flops: o2v4: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += V_blks_["baab_vovv"]("e,i,a,b") * t2_aaaa_vvoo("a,f,n,i") * s1_xb_Lv("I,b");

            // sigmas2_xbaa_Lvvo += 1.000000 <i,e||a,b>_bbbb t2_abab(f,a,n,i) s1_b(b) 
            // flops: o2v4: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= V_blks_["bbbb_vovv"]("e,i,a,b") * t2_abab_vvoo("f,a,n,i") * s1_xb_Lv("I,b");

            // sigmas2_xbaa_Lvvo += -1.000000 <f,i||a,b>_abab t2_abab(a,e,n,i) s1_b(b) 
            // flops: o2v4: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= V_blks_["abab_vovv"]("f,i,a,b") * t2_abab_vvoo("a,e,n,i") * s1_xb_Lv("I,b");
        }

        if (include_u2_) {

            // sigmas2_xbaa_Lvvo += 1.000000 <i,j||a,b>_abab t2_aaaa(a,f,n,i) s2_bbb(b,e,j) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("i,j,a,b") * t2_aaaa_vvoo("a,f,n,i") * s2_xbbb_Lvvo("I,b,e,j");

            // sigmas2_xbaa_Lvvo += 1.000000 <j,i||a,b>_bbbb t2_abab(f,a,n,i) s2_bbb(b,e,j) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += V_blks_["bbbb_oovv"]("j,i,a,b") * t2_abab_vvoo("f,a,n,i") * s2_xbbb_Lvvo("I,b,e,j");

            // sigmas2_xbaa_Lvvo += 0.500000 <j,i||a,b>_bbbb t2_abab(f,e,n,i) s2_bbb(a,b,j) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * t2_abab_vvoo("f,e,n,i") * s2_xbbb_Lvvo("I,a,b,j");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xbaa_Lvvo += 1.000000 <i,f||a,n>_aaaa r1_b(e) u1_aa(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= V_blks_["aaaa_vovo"]("f,i,a,n") * r1_xb_Lv("I,e") * u1_aa_vo("a,i");

            // sigmas2_xbaa_Lvvo += 1.000000 <f,i||n,a>_abab r1_b(e) u1_bb(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= V_blks_["abba_vovo"]("f,i,a,n") * r1_xb_Lv("I,e") * u1_bb_vo("a,i");

            // sigmas2_xbaa_Lvvo += -1.000000 <i,e||n,a>_abab r1_b(a) u1_aa(f,i) 
            // flops: o2v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= V_blks_["baba_vovo"]("e,i,a,n") * r1_xb_Lv("I,a") * u1_aa_vo("f,i");

            // sigmas2_xbaa_Lvvo += -1.000000 <f,i||n,a>_abab r1_b(a) u1_bb(e,i) 
            // flops: o2v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += V_blks_["abba_vovo"]("f,i,a,n") * r1_xb_Lv("I,a") * u1_bb_vo("e,i");

            // sigmas2_xbaa_Lvvo += 1.000000 <f,e||a,b>_abab r1_b(b) u1_aa(a,n) 
            // flops: o0v4L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += V_blks_["abab_vvvv"]("f,e,a,b") * r1_xb_Lv("I,b") * u1_aa_vo("a,n");

            // sigmas2_xbaa_Lvvo += 1.000000 <i,j||n,a>_abab r2_bbb(a,e,j) u1_aa(f,i) 
            // flops: o3v2L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= V_blks_["abba_oovo"]("i,j,a,n") * r2_xbbb_Lvvo("I,a,e,j") * u1_aa_vo("f,i");

            // sigmas2_xbaa_Lvvo += -1.000000 <f,i||a,b>_abab r2_bbb(b,e,i) u1_aa(a,n) 
            // flops: o1v4L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= V_blks_["abab_vovv"]("f,i,a,b") * r2_xbbb_Lvvo("I,b,e,i") * u1_aa_vo("a,n");
        }

        if (include_u2_) {

            // sigmas2_xbaa_Lvvo += -0.500000 <j,i||a,n>_aaaa r1_b(e) u2_aaaa(a,f,j,i) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 1 | mem: o3v2L1: 1, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovo"]("j,i,a,n") * r1_xb_Lv("I,e") * u2_aaaa_vvoo("a,f,j,i");

            // sigmas2_xbaa_Lvvo += -0.500000 <j,i||n,a>_abab r1_b(e) u2_abab(f,a,j,i) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 1 | mem: o3v2L1: 1, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abba_oovo"]("j,i,a,n") * r1_xb_Lv("I,e") * u2_abab_vvoo("f,a,j,i");

            // sigmas2_xbaa_Lvvo += -0.500000 <i,j||n,a>_abab r1_b(e) u2_abab(f,a,i,j) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 1 | mem: o3v2L1: 1, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abba_oovo"]("i,j,a,n") * r1_xb_Lv("I,e") * u2_abab_vvoo("f,a,i,j");

            // sigmas2_xbaa_Lvvo += 0.500000 <j,i||n,a>_abab r1_b(a) u2_abab(f,e,j,i) 
            // flops: o3v2L1: 1, o3v1L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o3v0L1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abba_oovo"]("j,i,a,n") * r1_xb_Lv("I,a") * u2_abab_vvoo("f,e,j,i");

            // sigmas2_xbaa_Lvvo += 0.500000 <i,j||n,a>_abab r1_b(a) u2_abab(f,e,i,j) 
            // flops: o3v2L1: 1, o3v1L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o3v0L1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abba_oovo"]("i,j,a,n") * r1_xb_Lv("I,a") * u2_abab_vvoo("f,e,i,j");

            // sigmas2_xbaa_Lvvo += -0.500000 <i,f||a,b>_aaaa r1_b(e) u2_aaaa(a,b,n,i) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["aaaa_vovv"]("f,i,a,b") * r1_xb_Lv("I,e") * u2_aaaa_vvoo("a,b,n,i");

            // sigmas2_xbaa_Lvvo += 0.500000 <f,i||a,b>_abab r1_b(e) u2_abab(a,b,n,i) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_vovv"]("f,i,a,b") * r1_xb_Lv("I,e") * u2_abab_vvoo("a,b,n,i");

            // sigmas2_xbaa_Lvvo += 0.500000 <f,i||b,a>_abab r1_b(e) u2_abab(b,a,n,i) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_vovv"]("f,i,b,a") * r1_xb_Lv("I,e") * u2_abab_vvoo("b,a,n,i");

            // sigmas2_xbaa_Lvvo += -1.000000 <i,e||a,b>_abab r1_b(b) u2_aaaa(a,f,n,i) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += V_blks_["baab_vovv"]("e,i,a,b") * r1_xb_Lv("I,b") * u2_aaaa_vvoo("a,f,n,i");

            // sigmas2_xbaa_Lvvo += 1.000000 <i,e||a,b>_bbbb r1_b(b) u2_abab(f,a,n,i) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= V_blks_["bbbb_vovv"]("e,i,a,b") * r1_xb_Lv("I,b") * u2_abab_vvoo("f,a,n,i");

            // sigmas2_xbaa_Lvvo += -1.000000 <f,i||a,b>_abab r1_b(b) u2_abab(a,e,n,i) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= V_blks_["abab_vovv"]("f,i,a,b") * r1_xb_Lv("I,b") * u2_abab_vvoo("a,e,n,i");

            // sigmas2_xbaa_Lvvo += 1.000000 <i,j||a,b>_abab r2_bbb(b,e,j) u2_aaaa(a,f,n,i) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("i,j,a,b") * r2_xbbb_Lvvo("I,b,e,j") * u2_aaaa_vvoo("a,f,n,i");

            // sigmas2_xbaa_Lvvo += 1.000000 <j,i||a,b>_bbbb r2_bbb(b,e,j) u2_abab(f,a,n,i) 
            // flops: o2v3L1: 2, o1v2L1: 1 | mem: o1v2L1: 3, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += V_blks_["bbbb_oovv"]("j,i,a,b") * r2_xbbb_Lvvo("I,b,e,j") * u2_abab_vvoo("f,a,n,i");

            // sigmas2_xbaa_Lvvo += 0.500000 <j,i||a,b>_bbbb r2_bbb(a,b,j) u2_abab(f,e,n,i) 
            // flops: o2v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r2_xbbb_Lvvo("I,a,b,j") * u2_abab_vvoo("f,e,n,i");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xbaa_Lvvo += -0.500000 <j,i||a,b>_aaaa r1_b(e) t2_aaaa(a,b,n,i) u1_aa(f,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r1_xb_Lv("I,e") * t2_aaaa_vvoo("a,b,n,i") * u1_aa_vo("f,j");

            // sigmas2_xbaa_Lvvo += -0.500000 <j,i||a,b>_abab r1_b(e) t2_abab(a,b,n,i) u1_aa(f,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("j,i,a,b") * r1_xb_Lv("I,e") * t2_abab_vvoo("a,b,n,i") * u1_aa_vo("f,j");

            // sigmas2_xbaa_Lvvo += -0.500000 <j,i||b,a>_abab r1_b(e) t2_abab(b,a,n,i) u1_aa(f,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("j,i,b,a") * r1_xb_Lv("I,e") * t2_abab_vvoo("b,a,n,i") * u1_aa_vo("f,j");

            // sigmas2_xbaa_Lvvo += -0.500000 <j,i||a,b>_aaaa r1_b(e) t2_aaaa(a,f,j,i) u1_aa(b,n) 
            // flops: o2v4L1: 1, o2v3L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o0v3L1: 1, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r1_xb_Lv("I,e") * t2_aaaa_vvoo("a,f,j,i") * u1_aa_vo("b,n");

            // sigmas2_xbaa_Lvvo += -0.500000 <j,i||b,a>_abab r1_b(e) t2_abab(f,a,j,i) u1_aa(b,n) 
            // flops: o2v4L1: 1, o2v3L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o0v3L1: 1, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("j,i,b,a") * r1_xb_Lv("I,e") * t2_abab_vvoo("f,a,j,i") * u1_aa_vo("b,n");

            // sigmas2_xbaa_Lvvo += -0.500000 <i,j||b,a>_abab r1_b(e) t2_abab(f,a,i,j) u1_aa(b,n) 
            // flops: o2v4L1: 1, o2v3L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o0v3L1: 1, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("i,j,b,a") * r1_xb_Lv("I,e") * t2_abab_vvoo("f,a,i,j") * u1_aa_vo("b,n");

            // sigmas2_xbaa_Lvvo += 1.000000 <j,i||a,b>_aaaa r1_b(e) t2_aaaa(a,f,n,i) u1_aa(b,j) 
            // flops: o3v4L1: 1, o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 2, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += V_blks_["aaaa_oovv"]("j,i,a,b") * r1_xb_Lv("I,e") * t2_aaaa_vvoo("a,f,n,i") * u1_aa_vo("b,j");

            // sigmas2_xbaa_Lvvo += -1.000000 <i,j||a,b>_abab r1_b(e) t2_aaaa(a,f,n,i) u1_bb(b,j) 
            // flops: o3v4L1: 1, o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 2, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= V_blks_["abab_oovv"]("i,j,a,b") * r1_xb_Lv("I,e") * t2_aaaa_vvoo("a,f,n,i") * u1_bb_vo("b,j");

            // sigmas2_xbaa_Lvvo += 1.000000 <j,i||b,a>_abab r1_b(e) t2_abab(f,a,n,i) u1_aa(b,j) 
            // flops: o3v4L1: 1, o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 2, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("j,i,b,a") * r1_xb_Lv("I,e") * t2_abab_vvoo("f,a,n,i") * u1_aa_vo("b,j");

            // sigmas2_xbaa_Lvvo += -1.000000 <j,i||a,b>_bbbb r1_b(e) t2_abab(f,a,n,i) u1_bb(b,j) 
            // flops: o3v4L1: 1, o2v3L1: 2, o1v2L1: 1 | mem: o2v3L1: 2, o1v2L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= V_blks_["bbbb_oovv"]("j,i,a,b") * r1_xb_Lv("I,e") * t2_abab_vvoo("f,a,n,i") * u1_bb_vo("b,j");

            // sigmas2_xbaa_Lvvo += 1.000000 <j,i||a,b>_abab r1_b(b) t2_abab(a,e,n,i) u1_aa(f,j) 
            // flops: o3v2L1: 1, o2v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o2v1L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("j,i,a,b") * r1_xb_Lv("I,b") * t2_abab_vvoo("a,e,n,i") * u1_aa_vo("f,j");

            // sigmas2_xbaa_Lvvo += 1.000000 <i,j||a,b>_abab r1_b(b) t2_aaaa(a,f,n,i) u1_bb(e,j) 
            // flops: o3v2L1: 1, o2v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o2v1L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("i,j,a,b") * r1_xb_Lv("I,b") * t2_aaaa_vvoo("a,f,n,i") * u1_bb_vo("e,j");

            // sigmas2_xbaa_Lvvo += 1.000000 <j,i||a,b>_bbbb r1_b(b) t2_abab(f,a,n,i) u1_bb(e,j) 
            // flops: o3v2L1: 1, o2v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o2v1L1: 2, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += V_blks_["bbbb_oovv"]("j,i,a,b") * r1_xb_Lv("I,b") * t2_abab_vvoo("f,a,n,i") * u1_bb_vo("e,j");

            // sigmas2_xbaa_Lvvo += 0.500000 <j,i||a,b>_abab r1_b(b) t2_abab(f,e,j,i) u1_aa(a,n) 
            // flops: o2v3L1: 1, o1v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("j,i,a,b") * r1_xb_Lv("I,b") * t2_abab_vvoo("f,e,j,i") * u1_aa_vo("a,n");

            // sigmas2_xbaa_Lvvo += 0.500000 <i,j||a,b>_abab r1_b(b) t2_abab(f,e,i,j) u1_aa(a,n) 
            // flops: o2v3L1: 1, o1v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("i,j,a,b") * r1_xb_Lv("I,b") * t2_abab_vvoo("f,e,i,j") * u1_aa_vo("a,n");

            // sigmas2_xbaa_Lvvo += -1.000000 <j,i||a,b>_abab r1_b(b) t2_abab(f,e,n,i) u1_aa(a,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= V_blks_["abab_oovv"]("j,i,a,b") * r1_xb_Lv("I,b") * t2_abab_vvoo("f,e,n,i") * u1_aa_vo("a,j");

            // sigmas2_xbaa_Lvvo += -1.000000 <j,i||a,b>_bbbb r1_b(b) t2_abab(f,e,n,i) u1_bb(a,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xbaa_Lvvo("I,e,f,n") -= V_blks_["bbbb_oovv"]("j,i,a,b") * r1_xb_Lv("I,b") * t2_abab_vvoo("f,e,n,i") * u1_bb_vo("a,j");
        }




/// ****** sigmas2_bbb(e,f,n) ****** ///



        if (include_u2_) {

            // sigmas2_xbbb_Lvvo += 1.000000 s2_bbb(e,f,n) w0 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += s2_xbbb_Lvvo("I,e,f,n") * w0;
        }

        if (include_u0_ && include_u2_) {

            // sigmas2_xbbb_Lvvo += 1.000000 r2_bbb(e,f,n) u0 w0 
            // flops: o1v2L1: 2, o1v2: 1 | mem: o1v2L1: 1, o1v2: 2, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += r2_xbbb_Lvvo("I,e,f,n") * u0 * w0;
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) r1_b(f) u1_bb(e,n) w0 
            // flops: o1v2L1: 4 | mem: o1v2L1: 2, o1v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = r1_xb_Lv("I,f") * u1_bb_vo("e,n") * w0;
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) f_bb(e,n) s1_b(f) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = F_blks_["bb_vo"]("e,n") * s1_xb_Lv("I,f");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmas2_xbbb_Lvvo += 1.000000 f_aa(i,i) s2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += scalar5 * s2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += 1.000000 f_bb(i,i) s2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += scalar6 * s2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += -1.000000 f_bb(i,n) s2_bbb(e,f,i) 
            // flops: o2v2L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= F_blks_["bb_oo"]("i,n") * s2_xbbb_Lvvo("I,e,f,i");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) f_bb(e,a) s2_bbb(a,f,n) 
            // flops: o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = F_blks_["bb_vv"]("e,a") * s2_xbbb_Lvvo("I,a,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) f_aa(i,a) t2_abab(a,e,i,n) s1_b(f) 
            // flops: o1v2L1: 3, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = F_blks_["aa_ov"]("i,a") * t2_abab_vvoo("a,e,i,n") * s1_xb_Lv("I,f");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) f_bb(i,a) t2_bbbb(a,e,n,i) s1_b(f) 
            // flops: o1v2L1: 3, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = F_blks_["bb_ov"]("i,a") * t2_bbbb_vvoo("a,e,n,i") * s1_xb_Lv("I,f");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 f_bb(i,a) t2_bbbb(e,f,n,i) s1_b(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += F_blks_["bb_ov"]("i,a") * t2_bbbb_vvoo("e,f,n,i") * s1_xb_Lv("I,a");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) f_bb(i,n) r1_b(f) u1_bb(e,i) 
            // flops: o2v2L1: 1, o1v2L1: 2, o2v1L1: 1 | mem: o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = F_blks_["bb_oo"]("i,n") * r1_xb_Lv("I,f") * u1_bb_vo("e,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) f_bb(e,a) r1_b(f) u1_bb(a,n) 
            // flops: o1v3L1: 1, o0v3L1: 1, o1v2L1: 2 | mem: o0v3L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = F_blks_["bb_vv"]("e,a") * r1_xb_Lv("I,f") * u1_bb_vo("a,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 f_aa(i,a) r2_bbb(e,f,n) u1_aa(a,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += scalar12 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += 1.000000 f_bb(i,a) r2_bbb(e,f,n) u1_bb(a,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += scalar13 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += -1.000000 f_bb(i,a) r2_bbb(e,f,i) u1_bb(a,n) 
            // flops: o1v3L1: 2, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= F_blks_["bb_ov"]("i,a") * r2_xbbb_Lvvo("I,e,f,i") * u1_bb_vo("a,n");

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) f_bb(i,a) r2_bbb(a,f,n) u1_bb(e,i) 
            // flops: o2v2L1: 2, o1v2L1: 2 | mem: o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = F_blks_["bb_ov"]("i,a") * r2_xbbb_Lvvo("I,a,f,n") * u1_bb_vo("e,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) f_aa(i,a) r1_b(f) u2_abab(a,e,i,n) 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = F_blks_["aa_ov"]("i,a") * r1_xb_Lv("I,f") * u2_abab_vvoo("a,e,i,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) f_bb(i,a) r1_b(f) u2_bbbb(a,e,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = F_blks_["bb_ov"]("i,a") * r1_xb_Lv("I,f") * u2_bbbb_vvoo("a,e,n,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 f_bb(i,a) r1_b(a) u2_bbbb(e,f,n,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o1v1L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += F_blks_["bb_ov"]("i,a") * r1_xb_Lv("I,a") * u2_bbbb_vvoo("e,f,n,i");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) d+_bb(e,n) r1_b(f) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_vo("e,n") * r1_xb_Lv("I,f");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 d+_aa(i,i) r2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= scalar17 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += -1.000000 d+_bb(i,i) r2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= scalar18 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += 1.000000 d+_bb(i,n) r2_bbb(e,f,i) 
            // flops: o2v2L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += dp_bb_oo("i,n") * r2_xbbb_Lvvo("I,e,f,i");

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) d+_bb(e,a) r2_bbb(a,f,n) 
            // flops: o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_vv("e,a") * r2_xbbb_Lvvo("I,a,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) d+_aa(i,a) r1_b(f) t2_abab(a,e,i,n) 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_aa_ov("i,a") * r1_xb_Lv("I,f") * t2_abab_vvoo("a,e,i,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) d+_bb(i,a) r1_b(f) t2_bbbb(a,e,n,i) 
            // flops: o2v3L1: 1, o1v2L1: 3 | mem: o1v2L1: 3, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("i,a") * r1_xb_Lv("I,f") * t2_bbbb_vvoo("a,e,n,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 d+_bb(i,a) r1_b(a) t2_bbbb(e,f,n,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o1v1L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= dp_bb_ov("i,a") * r1_xb_Lv("I,a") * t2_bbbb_vvoo("e,f,n,i");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) d-_bb(e,n) u0 s1_b(f) 
            // flops: o1v2L1: 3, o1v1: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_vo("e,n") * u0 * s1_xb_Lv("I,f");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u0_ && include_u2_) {

            // sigmas2_xbbb_Lvvo += -1.000000 d-_aa(i,i) u0 s2_bbb(e,f,n) 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= u0 * scalar7 * s2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += -1.000000 d-_bb(i,i) u0 s2_bbb(e,f,n) 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= u0 * scalar8 * s2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += 1.000000 d-_bb(i,n) u0 s2_bbb(e,f,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o2v0: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += dp_bb_oo("i,n") * u0 * s2_xbbb_Lvvo("I,e,f,i");

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(e,a) u0 s2_bbb(a,f,n) 
            // flops: o1v3L1: 1, o1v2L1: 2, o0v2: 1 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_vv("e,a") * u0 * s2_xbbb_Lvvo("I,a,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) d-_aa(i,i) u1_bb(e,n) s1_b(f) 
            // flops: o1v2L1: 2, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = scalar7 * s1_xb_Lv("I,f") * u1_bb_vo("e,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) d-_bb(i,i) u1_bb(e,n) s1_b(f) 
            // flops: o1v2L1: 2, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = scalar8 * s1_xb_Lv("I,f") * u1_bb_vo("e,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -2.000000 P(e,f) d-_bb(i,n) u1_bb(e,i) s1_b(f) 
            // flops: o1v2L1: 3, o2v1: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 2.000000 * dp_bb_oo("i,n") * u1_bb_vo("e,i") * s1_xb_Lv("I,f");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 2.000000 P(e,f) d-_bb(e,a) u1_bb(a,n) s1_b(f) 
            // flops: o1v2L1: 3, o1v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 2.000000 * dp_bb_vv("e,a") * u1_bb_vo("a,n") * s1_xb_Lv("I,f");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(e,a) u1_bb(f,n) s1_b(a) 
            // flops: o1v3L1: 1, o1v2L1: 2, o1v3: 1 | mem: o1v2L1: 2, o1v3: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_vv("e,a") * u1_bb_vo("f,n") * s1_xb_Lv("I,a");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -2.000000 d-_aa(i,a) u1_aa(a,i) s2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= 2.000000 * scalar0 * s2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += -2.000000 d-_bb(i,a) u1_bb(a,i) s2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= 2.000000 * scalar1 * s2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += 2.000000 d-_bb(i,a) u1_bb(a,n) s2_bbb(e,f,i) 
            // flops: o2v2L1: 1, o1v2L1: 1, o2v1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += 2.000000 * dp_bb_ov("i,a") * u1_bb_vo("a,n") * s2_xbbb_Lvvo("I,e,f,i");

            // sigmas2_xbbb_Lvvo += 2.000000 P(e,f) d-_bb(i,a) u1_bb(e,i) s2_bbb(a,f,n) 
            // flops: o1v3L1: 1, o1v2L1: 2, o1v2: 1 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 2.000000 * dp_bb_ov("i,a") * u1_bb_vo("e,i") * s2_xbbb_Lvvo("I,a,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(i,a) u1_bb(e,n) s2_bbb(a,f,i) 
            // flops: o2v3L1: 1, o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("i,a") * u1_bb_vo("e,n") * s2_xbbb_Lvvo("I,a,f,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 2.000000 P(e,f) d-_aa(i,a) u2_abab(a,e,i,n) s1_b(f) 
            // flops: o1v2L1: 3, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 2.000000 * dp_aa_ov("i,a") * u2_abab_vvoo("a,e,i,n") * s1_xb_Lv("I,f");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -2.000000 P(e,f) d-_bb(i,a) u2_bbbb(a,e,n,i) s1_b(f) 
            // flops: o1v2L1: 3, o2v2: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 2.000000 * dp_bb_ov("i,a") * u2_bbbb_vvoo("a,e,n,i") * s1_xb_Lv("I,f");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -2.000000 d-_bb(i,a) u2_bbbb(e,f,n,i) s1_b(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= 2.000000 * dp_bb_ov("i,a") * u2_bbbb_vvoo("e,f,n,i") * s1_xb_Lv("I,a");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) d-_aa(i,a) t2_abab(a,e,i,n) u0 s1_b(f) 
            // flops: o1v2L1: 3, o2v2: 1, o1v1: 1 | mem: o1v2L1: 2, o1v1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_aa_ov("i,a") * t2_abab_vvoo("a,e,i,n") * u0 * s1_xb_Lv("I,f");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(i,a) t2_bbbb(a,e,n,i) u0 s1_b(f) 
            // flops: o1v2L1: 3, o2v2: 1, o1v1: 1 | mem: o1v2L1: 2, o1v1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("i,a") * t2_bbbb_vvoo("a,e,n,i") * u0 * s1_xb_Lv("I,f");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 d-_bb(i,a) t2_bbbb(e,f,n,i) u0 s1_b(a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1, o1v3: 1 | mem: o1v2L1: 2, o1v3: 2, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= dp_bb_ov("i,a") * t2_bbbb_vvoo("e,f,n,i") * u0 * s1_xb_Lv("I,a");

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(i,n) r1_b(f) u0 u1_bb(e,i) 
            // flops: o1v2L1: 2, o2v1L1: 2, o2v2: 1 | mem: o1v2L1: 1, o2v1L1: 1, o1v2: 1, o2v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_oo("i,n") * r1_xb_Lv("I,f") * u0 * u1_bb_vo("e,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) d-_bb(e,a) r1_b(f) u0 u1_bb(a,n) 
            // flops: o0v3L1: 2, o1v2L1: 2, o1v3: 1 | mem: o0v3L1: 1, o1v2L1: 1, o0v3: 1, o1v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_vv("e,a") * r1_xb_Lv("I,f") * u0 * u1_bb_vo("a,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 d-_aa(i,a) r2_bbb(e,f,n) u0 u1_aa(a,i) 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= u0 * scalar0 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += -1.000000 d-_bb(i,a) r2_bbb(e,f,n) u0 u1_bb(a,i) 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= u0 * scalar1 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += 1.000000 d-_bb(i,a) r2_bbb(e,f,i) u0 u1_bb(a,n) 
            // flops: o1v3L1: 1, o0v3L1: 1, o1v2L1: 1, o1v3: 1 | mem: o0v3L1: 1, o1v2L1: 1, o0v3: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += dp_bb_ov("i,a") * r2_xbbb_Lvvo("I,e,f,i") * u0 * u1_bb_vo("a,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) d-_bb(i,a) r2_bbb(a,f,n) u0 u1_bb(e,i) 
            // flops: o2v2L1: 1, o1v2L1: 2, o2v1L1: 1, o2v2: 1 | mem: o1v2L1: 1, o2v1L1: 1, o1v2: 1, o2v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("i,a") * r2_xbbb_Lvvo("I,a,f,n") * u0 * u1_bb_vo("e,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u0_ && include_u2_) {

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) d-_aa(i,a) r1_b(f) u0 u2_abab(a,e,i,n) 
            // flops: o2v3: 1, o1v2L1: 4 | mem: o1v2L1: 2, o1v2: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_aa_ov("i,a") * r1_xb_Lv("I,f") * u0 * u2_abab_vvoo("a,e,i,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(i,a) r1_b(f) u0 u2_bbbb(a,e,n,i) 
            // flops: o2v3: 1, o1v2L1: 4 | mem: o1v2L1: 2, o1v2: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("i,a") * r1_xb_Lv("I,f") * u0 * u2_bbbb_vvoo("a,e,n,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 d-_bb(i,a) r1_b(a) u0 u2_bbbb(e,f,n,i) 
            // flops: o1v2L1: 1, o2v2: 1, o1v1L1: 1, o1v0L1: 1 | mem: o1v2L1: 1, o1v2: 1, o1v0L1: 1, o1v0: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= dp_bb_ov("i,a") * r1_xb_Lv("I,a") * u0 * u2_bbbb_vvoo("e,f,n,i");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) d-_aa(i,a) r1_b(f) u1_aa(a,i) u1_bb(e,n) 
            // flops: o1v2L1: 2, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = scalar0 * r1_xb_Lv("I,f") * u1_bb_vo("e,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) d-_bb(i,a) r1_b(f) u1_bb(a,i) u1_bb(e,n) 
            // flops: o1v2L1: 2, o1v2: 1, o0v1L1: 1 | mem: o1v2L1: 1, o1v2: 1, o0v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = scalar1 * r1_xb_Lv("I,f") * u1_bb_vo("e,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -2.000000 P(e,f) d-_bb(i,a) r1_b(f) u1_bb(a,n) u1_bb(e,i) 
            // flops: o2v2L1: 2, o1v2L1: 3 | mem: o1v2L1: 3, o2v1L1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 2.000000 * dp_bb_ov("i,a") * r1_xb_Lv("I,f") * u1_bb_vo("a,n") * u1_bb_vo("e,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) d-_bb(i,a) r1_b(a) u1_bb(e,i) u1_bb(f,n) 
            // flops: o1v2L1: 3, o1v1L1: 2 | mem: o1v2L1: 2, o0v1L1: 1, o1v0L1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("i,a") * r1_xb_Lv("I,a") * u1_bb_vo("e,i") * u1_bb_vo("f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmas2_xbbb_Lvvo += -0.500000 <j,i||j,i>_aaaa s2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * scalar9 * s2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += -0.500000 <j,i||j,i>_abab s2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * scalar10 * s2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += -0.500000 <i,j||i,j>_abab s2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * scalar10 * s2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += -0.500000 <j,i||j,i>_bbbb s2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * scalar11 * s2_xbbb_Lvvo("I,e,f,n");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xbbb_Lvvo += 1.000000 <e,f||a,n>_bbbb s1_b(a) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += V_blks_["bbbb_vvvo"]("e,f,a,n") * s1_xb_Lv("I,a");
        }

        if (include_u2_) {

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) <i,e||a,n>_bbbb s2_bbb(a,f,i) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_vovo"]("e,i,a,n") * s2_xbbb_Lvvo("I,a,f,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 0.500000 <e,f||a,b>_bbbb s2_bbb(a,b,n) 
            // flops: o1v4L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["bbbb_vvvv"]("e,f,a,b") * s2_xbbb_Lvvo("I,a,b,n");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xbbb_Lvvo += 0.500000 P(e,f) <j,i||a,n>_abab t2_abab(a,e,j,i) s1_b(f) 
            // flops: o3v2: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovo"]("j,i,a,n") * t2_abab_vvoo("a,e,j,i") * s1_xb_Lv("I,f");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 0.500000 P(e,f) <i,j||a,n>_abab t2_abab(a,e,i,j) s1_b(f) 
            // flops: o3v2: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovo"]("i,j,a,n") * t2_abab_vvoo("a,e,i,j") * s1_xb_Lv("I,f");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 0.500000 P(e,f) <j,i||a,n>_bbbb t2_bbbb(a,e,j,i) s1_b(f) 
            // flops: o3v2: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["bbbb_oovo"]("j,i,a,n") * t2_bbbb_vvoo("a,e,j,i") * s1_xb_Lv("I,f");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 0.500000 <j,i||a,n>_bbbb t2_bbbb(e,f,j,i) s1_b(a) 
            // flops: o3v3: 1, o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o1v3: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["bbbb_oovo"]("j,i,a,n") * t2_bbbb_vvoo("e,f,j,i") * s1_xb_Lv("I,a");

            // sigmas2_xbbb_Lvvo += -0.500000 P(e,f) <i,e||a,b>_abab t2_abab(a,b,i,n) s1_b(f) 
            // flops: o2v3: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["baab_vovv"]("e,i,a,b") * t2_abab_vvoo("a,b,i,n") * s1_xb_Lv("I,f");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -0.500000 P(e,f) <i,e||b,a>_abab t2_abab(b,a,i,n) s1_b(f) 
            // flops: o2v3: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["baab_vovv"]("e,i,b,a") * t2_abab_vvoo("b,a,i,n") * s1_xb_Lv("I,f");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 0.500000 P(e,f) <i,e||a,b>_bbbb t2_bbbb(a,b,n,i) s1_b(f) 
            // flops: o2v3: 1, o1v2L1: 3 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["bbbb_vovv"]("e,i,a,b") * t2_bbbb_vvoo("a,b,n,i") * s1_xb_Lv("I,f");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) <i,e||a,b>_abab t2_abab(a,f,i,n) s1_b(b) 
            // flops: o2v4: 1, o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v3: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["baab_vovv"]("e,i,a,b") * t2_abab_vvoo("a,f,i,n") * s1_xb_Lv("I,b");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) <i,e||a,b>_bbbb t2_bbbb(a,f,n,i) s1_b(b) 
            // flops: o2v4: 1, o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, o1v3: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_vovv"]("e,i,a,b") * t2_bbbb_vvoo("a,f,n,i") * s1_xb_Lv("I,b");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmas2_xbbb_Lvvo += 0.250000 <j,i||a,b>_aaaa t2_aaaa(a,b,j,i) s2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar2 * s2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += 0.250000 <j,i||a,b>_abab t2_abab(a,b,j,i) s2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar3 * s2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += 0.250000 <i,j||a,b>_abab t2_abab(a,b,i,j) s2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar3 * s2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += 0.250000 <j,i||b,a>_abab t2_abab(b,a,j,i) s2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar3 * s2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += 0.250000 <i,j||b,a>_abab t2_abab(b,a,i,j) s2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar3 * s2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += 0.250000 <j,i||a,b>_bbbb t2_bbbb(a,b,j,i) s2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar4 * s2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += -0.500000 <i,j||a,b>_abab t2_abab(a,b,i,n) s2_bbb(e,f,j) 
            // flops: o2v2L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("i,j,a,b") * t2_abab_vvoo("a,b,i,n") * s2_xbbb_Lvvo("I,e,f,j");

            // sigmas2_xbbb_Lvvo += -0.500000 <i,j||b,a>_abab t2_abab(b,a,i,n) s2_bbb(e,f,j) 
            // flops: o2v2L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("i,j,b,a") * t2_abab_vvoo("b,a,i,n") * s2_xbbb_Lvvo("I,e,f,j");

            // sigmas2_xbbb_Lvvo += -0.500000 <j,i||a,b>_bbbb t2_bbbb(a,b,n,i) s2_bbb(e,f,j) 
            // flops: o2v2L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * t2_bbbb_vvoo("a,b,n,i") * s2_xbbb_Lvvo("I,e,f,j");

            // sigmas2_xbbb_Lvvo += -0.500000 P(e,f) <j,i||a,b>_abab t2_abab(a,e,j,i) s2_bbb(b,f,n) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("j,i,a,b") * t2_abab_vvoo("a,e,j,i") * s2_xbbb_Lvvo("I,b,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -0.500000 P(e,f) <i,j||a,b>_abab t2_abab(a,e,i,j) s2_bbb(b,f,n) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("i,j,a,b") * t2_abab_vvoo("a,e,i,j") * s2_xbbb_Lvvo("I,b,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -0.500000 P(e,f) <j,i||a,b>_bbbb t2_bbbb(a,e,j,i) s2_bbb(b,f,n) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * t2_bbbb_vvoo("a,e,j,i") * s2_xbbb_Lvvo("I,b,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) <i,j||a,b>_abab t2_abab(a,e,i,n) s2_bbb(b,f,j) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("i,j,a,b") * t2_abab_vvoo("a,e,i,n") * s2_xbbb_Lvvo("I,b,f,j");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) <j,i||a,b>_bbbb t2_bbbb(a,e,n,i) s2_bbb(b,f,j) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_oovv"]("j,i,a,b") * t2_bbbb_vvoo("a,e,n,i") * s2_xbbb_Lvvo("I,b,f,j");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 0.250000 <j,i||a,b>_bbbb t2_bbbb(e,f,j,i) s2_bbb(a,b,n) 
            // flops: o1v4L1: 1, o2v4: 1, o1v2L1: 1 | mem: o0v4: 1, o1v2L1: 2, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += 0.250000 * V_blks_["bbbb_oovv"]("j,i,a,b") * t2_bbbb_vvoo("e,f,j,i") * s2_xbbb_Lvvo("I,a,b,n");

            // sigmas2_xbbb_Lvvo += -0.500000 <j,i||a,b>_bbbb t2_bbbb(e,f,n,i) s2_bbb(a,b,j) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * t2_bbbb_vvoo("e,f,n,i") * s2_xbbb_Lvvo("I,a,b,j");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) <i,e||a,n>_abab r1_b(f) u1_aa(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["baab_vovo"]("e,i,a,n") * r1_xb_Lv("I,f") * u1_aa_vo("a,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) <i,e||a,n>_bbbb r1_b(f) u1_bb(a,i) 
            // flops: o2v3L1: 2, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_vovo"]("e,i,a,n") * r1_xb_Lv("I,f") * u1_bb_vo("a,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) <i,e||a,n>_bbbb r1_b(a) u1_bb(f,i) 
            // flops: o2v2L1: 2, o1v2L1: 2 | mem: o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_vovo"]("e,i,a,n") * r1_xb_Lv("I,a") * u1_bb_vo("f,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 <e,f||a,b>_bbbb r1_b(b) u1_bb(a,n) 
            // flops: o0v4L1: 1, o1v3L1: 1, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= V_blks_["bbbb_vvvv"]("e,f,a,b") * r1_xb_Lv("I,b") * u1_bb_vo("a,n");

            // sigmas2_xbbb_Lvvo += -1.000000 <i,j||a,n>_abab r2_bbb(e,f,j) u1_aa(a,i) 
            // flops: o3v3L1: 1, o2v3L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= V_blks_["abab_oovo"]("i,j,a,n") * r2_xbbb_Lvvo("I,e,f,j") * u1_aa_vo("a,i");

            // sigmas2_xbbb_Lvvo += 1.000000 <j,i||a,n>_bbbb r2_bbb(e,f,j) u1_bb(a,i) 
            // flops: o3v3L1: 1, o2v3L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += V_blks_["bbbb_oovo"]("j,i,a,n") * r2_xbbb_Lvvo("I,e,f,j") * u1_bb_vo("a,i");

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) <j,i||a,n>_bbbb r2_bbb(a,f,j) u1_bb(e,i) 
            // flops: o3v2L1: 1, o2v2L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_oovo"]("j,i,a,n") * r2_xbbb_Lvvo("I,a,f,j") * u1_bb_vo("e,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) <i,e||a,b>_abab r2_bbb(b,f,n) u1_aa(a,i) 
            // flops: o2v4L1: 1, o2v3L1: 1, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["baab_vovv"]("e,i,a,b") * r2_xbbb_Lvvo("I,b,f,n") * u1_aa_vo("a,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) <i,e||a,b>_bbbb r2_bbb(b,f,n) u1_bb(a,i) 
            // flops: o2v4L1: 1, o2v3L1: 1, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_vovv"]("e,i,a,b") * r2_xbbb_Lvvo("I,b,f,n") * u1_bb_vo("a,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) <i,e||a,b>_bbbb r2_bbb(b,f,i) u1_bb(a,n) 
            // flops: o1v4L1: 1, o1v3L1: 1, o1v2L1: 2 | mem: o0v3L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_vovv"]("e,i,a,b") * r2_xbbb_Lvvo("I,b,f,i") * u1_bb_vo("a,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 0.500000 P(e,f) <i,e||a,b>_bbbb r2_bbb(a,b,n) u1_bb(f,i) 
            // flops: o2v3L1: 1, o2v2L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["bbbb_vovv"]("e,i,a,b") * r2_xbbb_Lvvo("I,a,b,n") * u1_bb_vo("f,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmas2_xbbb_Lvvo += 0.500000 P(e,f) <j,i||a,n>_abab r1_b(f) u2_abab(a,e,j,i) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovo"]("j,i,a,n") * r1_xb_Lv("I,f") * u2_abab_vvoo("a,e,j,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 0.500000 P(e,f) <i,j||a,n>_abab r1_b(f) u2_abab(a,e,i,j) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovo"]("i,j,a,n") * r1_xb_Lv("I,f") * u2_abab_vvoo("a,e,i,j");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 0.500000 P(e,f) <j,i||a,n>_bbbb r1_b(f) u2_bbbb(a,e,j,i) 
            // flops: o3v3L1: 1, o3v2L1: 1, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["bbbb_oovo"]("j,i,a,n") * r1_xb_Lv("I,f") * u2_bbbb_vvoo("a,e,j,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 0.500000 <j,i||a,n>_bbbb r1_b(a) u2_bbbb(e,f,j,i) 
            // flops: o3v2L1: 1, o3v1L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, o3v0L1: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["bbbb_oovo"]("j,i,a,n") * r1_xb_Lv("I,a") * u2_bbbb_vvoo("e,f,j,i");

            // sigmas2_xbbb_Lvvo += -0.500000 P(e,f) <i,e||a,b>_abab r1_b(f) u2_abab(a,b,i,n) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 2 | mem: o1v4L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["baab_vovv"]("e,i,a,b") * r1_xb_Lv("I,f") * u2_abab_vvoo("a,b,i,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -0.500000 P(e,f) <i,e||b,a>_abab r1_b(f) u2_abab(b,a,i,n) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 2 | mem: o1v4L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["baab_vovv"]("e,i,b,a") * r1_xb_Lv("I,f") * u2_abab_vvoo("b,a,i,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 0.500000 P(e,f) <i,e||a,b>_bbbb r1_b(f) u2_bbbb(a,b,n,i) 
            // flops: o2v4L1: 1, o1v4L1: 1, o1v2L1: 2 | mem: o1v4L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["bbbb_vovv"]("e,i,a,b") * r1_xb_Lv("I,f") * u2_bbbb_vvoo("a,b,n,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) <i,e||a,b>_abab r1_b(b) u2_abab(a,f,i,n) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["baab_vovv"]("e,i,a,b") * r1_xb_Lv("I,b") * u2_abab_vvoo("a,f,i,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) <i,e||a,b>_bbbb r1_b(b) u2_bbbb(a,f,n,i) 
            // flops: o2v3L1: 1, o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 3, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_vovv"]("e,i,a,b") * r1_xb_Lv("I,b") * u2_bbbb_vvoo("a,f,n,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 0.250000 <j,i||a,b>_aaaa r2_bbb(e,f,n) u2_aaaa(a,b,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar14 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += 0.250000 <j,i||a,b>_abab r2_bbb(e,f,n) u2_abab(a,b,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar15 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += 0.250000 <i,j||a,b>_abab r2_bbb(e,f,n) u2_abab(a,b,i,j) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar15 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += 0.250000 <j,i||b,a>_abab r2_bbb(e,f,n) u2_abab(b,a,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar15 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += 0.250000 <i,j||b,a>_abab r2_bbb(e,f,n) u2_abab(b,a,i,j) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar15 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += 0.250000 <j,i||a,b>_bbbb r2_bbb(e,f,n) u2_bbbb(a,b,j,i) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar16 * r2_xbbb_Lvvo("I,e,f,n");

            // sigmas2_xbbb_Lvvo += -0.500000 <i,j||a,b>_abab r2_bbb(e,f,j) u2_abab(a,b,i,n) 
            // flops: o2v4L1: 2, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("i,j,a,b") * r2_xbbb_Lvvo("I,e,f,j") * u2_abab_vvoo("a,b,i,n");

            // sigmas2_xbbb_Lvvo += -0.500000 <i,j||b,a>_abab r2_bbb(e,f,j) u2_abab(b,a,i,n) 
            // flops: o2v4L1: 2, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("i,j,b,a") * r2_xbbb_Lvvo("I,e,f,j") * u2_abab_vvoo("b,a,i,n");

            // sigmas2_xbbb_Lvvo += -0.500000 <j,i||a,b>_bbbb r2_bbb(e,f,j) u2_bbbb(a,b,n,i) 
            // flops: o2v4L1: 2, o1v2L1: 1 | mem: o1v4L1: 1, o1v2L1: 2, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r2_xbbb_Lvvo("I,e,f,j") * u2_bbbb_vvoo("a,b,n,i");

            // sigmas2_xbbb_Lvvo += -0.500000 P(e,f) <j,i||a,b>_abab r2_bbb(b,f,n) u2_abab(a,e,j,i) 
            // flops: o3v3L1: 2, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("j,i,a,b") * r2_xbbb_Lvvo("I,b,f,n") * u2_abab_vvoo("a,e,j,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -0.500000 P(e,f) <i,j||a,b>_abab r2_bbb(b,f,n) u2_abab(a,e,i,j) 
            // flops: o3v3L1: 2, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("i,j,a,b") * r2_xbbb_Lvvo("I,b,f,n") * u2_abab_vvoo("a,e,i,j");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -0.500000 P(e,f) <j,i||a,b>_bbbb r2_bbb(b,f,n) u2_bbbb(a,e,j,i) 
            // flops: o3v3L1: 2, o1v2L1: 2 | mem: o3v2L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r2_xbbb_Lvvo("I,b,f,n") * u2_bbbb_vvoo("a,e,j,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) <i,j||a,b>_abab r2_bbb(b,f,j) u2_abab(a,e,i,n) 
            // flops: o2v3L1: 2, o1v2L1: 2 | mem: o1v2L1: 3, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("i,j,a,b") * r2_xbbb_Lvvo("I,b,f,j") * u2_abab_vvoo("a,e,i,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) <j,i||a,b>_bbbb r2_bbb(b,f,j) u2_bbbb(a,e,n,i) 
            // flops: o2v3L1: 2, o1v2L1: 2 | mem: o1v2L1: 3, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_oovv"]("j,i,a,b") * r2_xbbb_Lvvo("I,b,f,j") * u2_bbbb_vvoo("a,e,n,i");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 0.250000 <j,i||a,b>_bbbb r2_bbb(a,b,n) u2_bbbb(e,f,j,i) 
            // flops: o3v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o3v0L1: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += 0.250000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r2_xbbb_Lvvo("I,a,b,n") * u2_bbbb_vvoo("e,f,j,i");

            // sigmas2_xbbb_Lvvo += -0.500000 <j,i||a,b>_bbbb r2_bbb(a,b,j) u2_bbbb(e,f,n,i) 
            // flops: o2v2L1: 2, o1v2L1: 1 | mem: o1v2L1: 2, o1v0L1: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r2_xbbb_Lvvo("I,a,b,j") * u2_bbbb_vvoo("e,f,n,i");
        }

        if (include_u1_ && include_u2_) {

            // sigmas2_xbbb_Lvvo += 0.500000 P(e,f) <i,j||a,b>_abab r1_b(f) t2_abab(a,b,i,n) u1_bb(e,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("i,j,a,b") * r1_xb_Lv("I,f") * t2_abab_vvoo("a,b,i,n") * u1_bb_vo("e,j");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 0.500000 P(e,f) <i,j||b,a>_abab r1_b(f) t2_abab(b,a,i,n) u1_bb(e,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("i,j,b,a") * r1_xb_Lv("I,f") * t2_abab_vvoo("b,a,i,n") * u1_bb_vo("e,j");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 0.500000 P(e,f) <j,i||a,b>_bbbb r1_b(f) t2_bbbb(a,b,n,i) u1_bb(e,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 2 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r1_xb_Lv("I,f") * t2_bbbb_vvoo("a,b,n,i") * u1_bb_vo("e,j");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 0.500000 P(e,f) <j,i||a,b>_abab r1_b(f) t2_abab(a,e,j,i) u1_bb(b,n) 
            // flops: o2v4L1: 1, o2v3L1: 1, o1v3L1: 1, o1v2L1: 2 | mem: o2v3L1: 1, o0v3L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("j,i,a,b") * r1_xb_Lv("I,f") * t2_abab_vvoo("a,e,j,i") * u1_bb_vo("b,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 0.500000 P(e,f) <i,j||a,b>_abab r1_b(f) t2_abab(a,e,i,j) u1_bb(b,n) 
            // flops: o2v4L1: 1, o2v3L1: 1, o1v3L1: 1, o1v2L1: 2 | mem: o2v3L1: 1, o0v3L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("i,j,a,b") * r1_xb_Lv("I,f") * t2_abab_vvoo("a,e,i,j") * u1_bb_vo("b,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 0.500000 P(e,f) <j,i||a,b>_bbbb r1_b(f) t2_bbbb(a,e,j,i) u1_bb(b,n) 
            // flops: o2v4L1: 1, o2v3L1: 1, o1v3L1: 1, o1v2L1: 2 | mem: o2v3L1: 1, o0v3L1: 1, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r1_xb_Lv("I,f") * t2_bbbb_vvoo("a,e,j,i") * u1_bb_vo("b,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) <j,i||a,b>_aaaa r1_b(f) t2_abab(a,e,i,n) u1_aa(b,j) 
            // flops: o3v4L1: 1, o2v3L1: 2, o1v2L1: 2 | mem: o2v3L1: 2, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["aaaa_oovv"]("j,i,a,b") * r1_xb_Lv("I,f") * t2_abab_vvoo("a,e,i,n") * u1_aa_vo("b,j");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) <i,j||a,b>_abab r1_b(f) t2_abab(a,e,i,n) u1_bb(b,j) 
            // flops: o3v4L1: 1, o2v3L1: 2, o1v2L1: 2 | mem: o2v3L1: 2, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("i,j,a,b") * r1_xb_Lv("I,f") * t2_abab_vvoo("a,e,i,n") * u1_bb_vo("b,j");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) <j,i||b,a>_abab r1_b(f) t2_bbbb(a,e,n,i) u1_aa(b,j) 
            // flops: o3v4L1: 1, o2v3L1: 2, o1v2L1: 2 | mem: o2v3L1: 2, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("j,i,b,a") * r1_xb_Lv("I,f") * t2_bbbb_vvoo("a,e,n,i") * u1_aa_vo("b,j");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -1.000000 P(e,f) <j,i||a,b>_bbbb r1_b(f) t2_bbbb(a,e,n,i) u1_bb(b,j) 
            // flops: o3v4L1: 1, o2v3L1: 2, o1v2L1: 2 | mem: o2v3L1: 2, o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_oovv"]("j,i,a,b") * r1_xb_Lv("I,f") * t2_bbbb_vvoo("a,e,n,i") * u1_bb_vo("b,j");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) <i,j||a,b>_abab r1_b(b) t2_abab(a,e,i,n) u1_bb(f,j) 
            // flops: o3v2L1: 1, o2v2L1: 2, o1v2L1: 2 | mem: o1v2L1: 2, o2v1L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("i,j,a,b") * r1_xb_Lv("I,b") * t2_abab_vvoo("a,e,i,n") * u1_bb_vo("f,j");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += 1.000000 P(e,f) <j,i||a,b>_bbbb r1_b(b) t2_bbbb(a,e,n,i) u1_bb(f,j) 
            // flops: o3v2L1: 1, o2v2L1: 2, o1v2L1: 2 | mem: o1v2L1: 2, o2v1L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_oovv"]("j,i,a,b") * r1_xb_Lv("I,b") * t2_bbbb_vvoo("a,e,n,i") * u1_bb_vo("f,j");
            sigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmas2_xbbb_Lvvo += -0.500000 <j,i||a,b>_bbbb r1_b(b) t2_bbbb(e,f,j,i) u1_bb(a,n) 
            // flops: o2v3L1: 1, o1v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o0v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r1_xb_Lv("I,b") * t2_bbbb_vvoo("e,f,j,i") * u1_bb_vo("a,n");

            // sigmas2_xbbb_Lvvo += 1.000000 <j,i||a,b>_abab r1_b(b) t2_bbbb(e,f,n,i) u1_aa(a,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("j,i,a,b") * r1_xb_Lv("I,b") * t2_bbbb_vvoo("e,f,n,i") * u1_aa_vo("a,j");

            // sigmas2_xbbb_Lvvo += 1.000000 <j,i||a,b>_bbbb r1_b(b) t2_bbbb(e,f,n,i) u1_bb(a,j) 
            // flops: o3v3L1: 1, o2v3L1: 1, o2v2L1: 1, o1v2L1: 1 | mem: o2v3L1: 1, o1v2L1: 2, o2v1L1: 1, 
            sigmas2_xbbb_Lvvo("I,e,f,n") += V_blks_["bbbb_oovv"]("j,i,a,b") * r1_xb_Lv("I,b") * t2_bbbb_vvoo("e,f,n,i") * u1_bb_vo("a,j");
        }




/// ****** sigmam1_a(e) ****** ///



        if (include_u1_) {

            // sigmam1_xa_Lv += 1.000000 m1_a(e) w0 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xa_Lv("I,e") += m1_xa_Lv("I,e") * w0;

            // sigmam1_xa_Lv += 1.000000 f_aa(i,i) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xa_Lv("I,e") += scalar5 * m1_xa_Lv("I,e");

            // sigmam1_xa_Lv += 1.000000 f_bb(i,i) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xa_Lv("I,e") += scalar6 * m1_xa_Lv("I,e");

            // sigmam1_xa_Lv += 1.000000 f_aa(a,e) m1_a(a) 
            // flops: o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") += F_blks_["aa_vv"]("a,e") * m1_xa_Lv("I,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmam1_xa_Lv += -1.000000 f_aa(a,i) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") -= F_blks_["aa_vo"]("a,i") * m2_xaaa_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += -1.000000 f_bb(a,i) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") -= F_blks_["bb_vo"]("a,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += 1.000000 f_aa(j,b) t2_aaaa(b,a,i,j) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") += F_blks_["aa_ov"]("j,b") * t2_aaaa_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += -1.000000 f_aa(j,b) t2_abab(b,a,j,i) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") -= F_blks_["aa_ov"]("j,b") * t2_abab_vvoo("b,a,j,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += -1.000000 f_bb(j,b) t2_abab(a,b,i,j) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") -= F_blks_["bb_ov"]("j,b") * t2_abab_vvoo("a,b,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += 1.000000 f_bb(j,b) t2_bbbb(b,a,i,j) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") += F_blks_["bb_ov"]("j,b") * t2_bbbb_vvoo("b,a,i,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += 0.500000 f_aa(j,e) t2_aaaa(b,a,i,j) m2_aaa(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") += 0.500000 * F_blks_["aa_ov"]("j,e") * t2_aaaa_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,b,a");

            // sigmam1_xa_Lv += 0.500000 f_aa(j,e) t2_abab(a,b,j,i) m2_bba(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") += 0.500000 * F_blks_["aa_ov"]("j,e") * t2_abab_vvoo("a,b,j,i") * m2_xbba_Lovv("I,i,b,a");
        }

        if (include_u1_) {

            // sigmam1_xa_Lv += -1.000000 d-_aa(i,i) l1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xa_Lv("I,e") -= scalar7 * l1_xa_Lv("I,e");

            // sigmam1_xa_Lv += -1.000000 d-_bb(i,i) l1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xa_Lv("I,e") -= scalar8 * l1_xa_Lv("I,e");

            // sigmam1_xa_Lv += -1.000000 d-_aa(a,e) l1_a(a) 
            // flops: o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") -= dp_aa_vv("a,e") * l1_xa_Lv("I,a");

            // sigmam1_xa_Lv += 1.000000 d-_aa(a,i) l2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") += dp_aa_vo("a,i") * l2_xaaa_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += 1.000000 d-_bb(a,i) l2_bba(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") += dp_bb_vo("a,i") * l2_xbba_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += -1.000000 d-_aa(j,b) l2_aaa(i,a,e) t2_aaaa(b,a,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") -= dp_aa_ov("j,b") * l2_xaaa_Lovv("I,i,a,e") * t2_aaaa_vvoo("b,a,i,j");

            // sigmam1_xa_Lv += 1.000000 d-_bb(j,b) l2_aaa(i,a,e) t2_abab(a,b,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") += dp_bb_ov("j,b") * l2_xaaa_Lovv("I,i,a,e") * t2_abab_vvoo("a,b,i,j");

            // sigmam1_xa_Lv += 1.000000 d-_aa(j,b) l2_bba(i,a,e) t2_abab(b,a,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") += dp_aa_ov("j,b") * l2_xbba_Lovv("I,i,a,e") * t2_abab_vvoo("b,a,j,i");

            // sigmam1_xa_Lv += -1.000000 d-_bb(j,b) l2_bba(i,a,e) t2_bbbb(b,a,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") -= dp_bb_ov("j,b") * l2_xbba_Lovv("I,i,a,e") * t2_bbbb_vvoo("b,a,i,j");

            // sigmam1_xa_Lv += -0.500000 d-_aa(j,e) l2_aaa(i,b,a) t2_aaaa(b,a,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") -= 0.500000 * dp_aa_ov("j,e") * l2_xaaa_Lovv("I,i,b,a") * t2_aaaa_vvoo("b,a,i,j");

            // sigmam1_xa_Lv += -0.500000 d-_aa(j,e) l2_bba(i,b,a) t2_abab(a,b,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") -= 0.500000 * dp_aa_ov("j,e") * l2_xbba_Lovv("I,i,b,a") * t2_abab_vvoo("a,b,j,i");
        }

        if (include_u0_ && include_u1_) {

            // sigmam1_xa_Lv += -1.000000 d-_aa(i,i) u0 m1_a(e) 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmam1_xa_Lv("I,e") -= u0 * scalar7 * m1_xa_Lv("I,e");

            // sigmam1_xa_Lv += -1.000000 d-_bb(i,i) u0 m1_a(e) 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmam1_xa_Lv("I,e") -= u0 * scalar8 * m1_xa_Lv("I,e");

            // sigmam1_xa_Lv += -1.000000 d-_aa(a,e) u0 m1_a(a) 
            // flops: o0v2L1: 1, o0v1L1: 1, o0v2: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmam1_xa_Lv("I,e") -= dp_aa_vv("a,e") * u0 * m1_xa_Lv("I,a");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmam1_xa_Lv += 1.000000 d-_aa(a,i) u0 m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") += dp_aa_vo("a,i") * u0 * m2_xaaa_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += 1.000000 d-_bb(a,i) u0 m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") += dp_bb_vo("a,i") * u0 * m2_xbba_Lovv("I,i,a,e");
        }

        if (include_u1_) {

            // sigmam1_xa_Lv += -2.000000 d-_aa(i,a) u1_aa(a,i) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xa_Lv("I,e") -= 2.000000 * scalar0 * m1_xa_Lv("I,e");

            // sigmam1_xa_Lv += -2.000000 d-_bb(i,a) u1_bb(a,i) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xa_Lv("I,e") -= 2.000000 * scalar1 * m1_xa_Lv("I,e");

            // sigmam1_xa_Lv += 2.000000 d-_aa(i,e) u1_aa(a,i) m1_a(a) 
            // flops: o0v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmam1_xa_Lv("I,e") += 2.000000 * dp_aa_ov("i,e") * u1_aa_vo("a,i") * m1_xa_Lv("I,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmam1_xa_Lv += 1.000000 d-_aa(i,i) u1_aa(a,j) m2_aaa(j,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") += u1_aa_vo("a,j") * scalar7 * m2_xaaa_Lovv("I,j,a,e");

            // sigmam1_xa_Lv += 1.000000 d-_bb(i,i) u1_aa(a,j) m2_aaa(j,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") += u1_aa_vo("a,j") * scalar8 * m2_xaaa_Lovv("I,j,a,e");

            // sigmam1_xa_Lv += 1.000000 d-_aa(i,i) u1_bb(a,j) m2_bba(j,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") += u1_bb_vo("a,j") * scalar7 * m2_xbba_Lovv("I,j,a,e");

            // sigmam1_xa_Lv += 1.000000 d-_bb(i,i) u1_bb(a,j) m2_bba(j,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") += u1_bb_vo("a,j") * scalar8 * m2_xbba_Lovv("I,j,a,e");

            // sigmam1_xa_Lv += -2.000000 d-_aa(j,i) u1_aa(a,j) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") -= 2.000000 * dp_aa_oo("j,i") * u1_aa_vo("a,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += -2.000000 d-_bb(j,i) u1_bb(a,j) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") -= 2.000000 * dp_bb_oo("j,i") * u1_bb_vo("a,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += 2.000000 d-_aa(a,b) u1_aa(b,i) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") += 2.000000 * dp_aa_vv("a,b") * u1_aa_vo("b,i") * m2_xaaa_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += 2.000000 d-_bb(a,b) u1_bb(b,i) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") += 2.000000 * dp_bb_vv("a,b") * u1_bb_vo("b,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += -1.000000 d-_aa(b,e) u1_aa(a,i) m2_aaa(i,b,a) 
            // flops: o1v3L1: 1, o1v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") -= dp_aa_vv("b,e") * u1_aa_vo("a,i") * m2_xaaa_Lovv("I,i,b,a");

            // sigmam1_xa_Lv += -2.000000 d-_aa(j,b) u2_aaaa(b,a,i,j) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") -= 2.000000 * dp_aa_ov("j,b") * u2_aaaa_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += 2.000000 d-_aa(j,b) u2_abab(b,a,j,i) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") += 2.000000 * dp_aa_ov("j,b") * u2_abab_vvoo("b,a,j,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += 2.000000 d-_bb(j,b) u2_abab(a,b,i,j) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") += 2.000000 * dp_bb_ov("j,b") * u2_abab_vvoo("a,b,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += -2.000000 d-_bb(j,b) u2_bbbb(b,a,i,j) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") -= 2.000000 * dp_bb_ov("j,b") * u2_bbbb_vvoo("b,a,i,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += -1.000000 d-_aa(j,e) u2_aaaa(b,a,i,j) m2_aaa(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") -= dp_aa_ov("j,e") * u2_aaaa_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,b,a");

            // sigmam1_xa_Lv += -1.000000 d-_aa(j,e) u2_abab(a,b,j,i) m2_bba(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") -= dp_aa_ov("j,e") * u2_abab_vvoo("a,b,j,i") * m2_xbba_Lovv("I,i,b,a");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmam1_xa_Lv += -1.000000 d-_aa(j,b) t2_aaaa(b,a,i,j) u0 m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 2, 
            sigmam1_xa_Lv("I,e") -= dp_aa_ov("j,b") * t2_aaaa_vvoo("b,a,i,j") * u0 * m2_xaaa_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += 1.000000 d-_aa(j,b) t2_abab(b,a,j,i) u0 m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 2, 
            sigmam1_xa_Lv("I,e") += dp_aa_ov("j,b") * t2_abab_vvoo("b,a,j,i") * u0 * m2_xbba_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += 1.000000 d-_bb(j,b) t2_abab(a,b,i,j) u0 m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 2, 
            sigmam1_xa_Lv("I,e") += dp_bb_ov("j,b") * t2_abab_vvoo("a,b,i,j") * u0 * m2_xaaa_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += -1.000000 d-_bb(j,b) t2_bbbb(b,a,i,j) u0 m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 2, 
            sigmam1_xa_Lv("I,e") -= dp_bb_ov("j,b") * t2_bbbb_vvoo("b,a,i,j") * u0 * m2_xbba_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += -0.500000 d-_aa(j,e) t2_aaaa(b,a,i,j) u0 m2_aaa(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v3: 1, o0v1L1: 1 | mem: o1v3: 2, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") -= 0.500000 * dp_aa_ov("j,e") * t2_aaaa_vvoo("b,a,i,j") * u0 * m2_xaaa_Lovv("I,i,b,a");

            // sigmam1_xa_Lv += -0.500000 d-_aa(j,e) t2_abab(a,b,j,i) u0 m2_bba(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v3: 1, o0v1L1: 1 | mem: o1v3: 2, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") -= 0.500000 * dp_aa_ov("j,e") * t2_abab_vvoo("a,b,j,i") * u0 * m2_xbba_Lovv("I,i,b,a");
        }

        if (include_u1_) {

            // sigmam1_xa_Lv += -0.500000 <j,i||j,i>_aaaa m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xa_Lv("I,e") -= 0.500000 * scalar9 * m1_xa_Lv("I,e");

            // sigmam1_xa_Lv += -0.500000 <j,i||j,i>_abab m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xa_Lv("I,e") -= 0.500000 * scalar10 * m1_xa_Lv("I,e");

            // sigmam1_xa_Lv += -0.500000 <i,j||i,j>_abab m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xa_Lv("I,e") -= 0.500000 * scalar10 * m1_xa_Lv("I,e");

            // sigmam1_xa_Lv += -0.500000 <j,i||j,i>_bbbb m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xa_Lv("I,e") -= 0.500000 * scalar11 * m1_xa_Lv("I,e");
        }

        if (include_u1_ && include_u2_) {

            // sigmam1_xa_Lv += 0.500000 <b,a||e,i>_aaaa m2_aaa(i,b,a) 
            // flops: o1v3L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") += 0.500000 * V_blks_["aaaa_vvvo"]("b,a,e,i") * m2_xaaa_Lovv("I,i,b,a");

            // sigmam1_xa_Lv += -0.500000 <a,b||e,i>_abab m2_bba(i,b,a) 
            // flops: o1v3L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_vvvo"]("a,b,e,i") * m2_xbba_Lovv("I,i,b,a");
        }

        if (include_u1_) {

            // sigmam1_xa_Lv += 0.250000 <j,i||a,b>_aaaa t2_aaaa(a,b,j,i) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xa_Lv("I,e") += 0.250000 * scalar2 * m1_xa_Lv("I,e");

            // sigmam1_xa_Lv += 0.250000 <j,i||a,b>_abab t2_abab(a,b,j,i) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xa_Lv("I,e") += 0.250000 * scalar3 * m1_xa_Lv("I,e");

            // sigmam1_xa_Lv += 0.250000 <i,j||a,b>_abab t2_abab(a,b,i,j) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xa_Lv("I,e") += 0.250000 * scalar3 * m1_xa_Lv("I,e");

            // sigmam1_xa_Lv += 0.250000 <j,i||b,a>_abab t2_abab(b,a,j,i) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xa_Lv("I,e") += 0.250000 * scalar3 * m1_xa_Lv("I,e");

            // sigmam1_xa_Lv += 0.250000 <i,j||b,a>_abab t2_abab(b,a,i,j) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xa_Lv("I,e") += 0.250000 * scalar3 * m1_xa_Lv("I,e");

            // sigmam1_xa_Lv += 0.250000 <j,i||a,b>_bbbb t2_bbbb(a,b,j,i) m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xa_Lv("I,e") += 0.250000 * scalar4 * m1_xa_Lv("I,e");

            // sigmam1_xa_Lv += -0.500000 <j,i||b,e>_aaaa t2_aaaa(b,a,j,i) m1_a(a) 
            // flops: o2v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmam1_xa_Lv("I,e") -= 0.500000 * V_blks_["aaaa_oovv"]("j,i,b,e") * t2_aaaa_vvoo("b,a,j,i") * m1_xa_Lv("I,a");

            // sigmam1_xa_Lv += -0.500000 <j,i||e,b>_abab t2_abab(a,b,j,i) m1_a(a) 
            // flops: o2v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmam1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("j,i,e,b") * t2_abab_vvoo("a,b,j,i") * m1_xa_Lv("I,a");

            // sigmam1_xa_Lv += -0.500000 <i,j||e,b>_abab t2_abab(a,b,i,j) m1_a(a) 
            // flops: o2v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmam1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("i,j,e,b") * t2_abab_vvoo("a,b,i,j") * m1_xa_Lv("I,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmam1_xa_Lv += 0.500000 <k,j||b,i>_aaaa t2_aaaa(b,a,k,j) m2_aaa(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") += 0.500000 * V_blks_["aaaa_oovo"]("k,j,b,i") * t2_aaaa_vvoo("b,a,k,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += 0.500000 <k,j||b,i>_abab t2_abab(b,a,k,j) m2_bba(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") += 0.500000 * V_blks_["abab_oovo"]("k,j,b,i") * t2_abab_vvoo("b,a,k,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += 0.500000 <j,k||b,i>_abab t2_abab(b,a,j,k) m2_bba(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") += 0.500000 * V_blks_["abab_oovo"]("j,k,b,i") * t2_abab_vvoo("b,a,j,k") * m2_xbba_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += 0.500000 <k,j||i,b>_abab t2_abab(a,b,k,j) m2_aaa(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") -= 0.500000 * V_blks_["abba_oovo"]("k,j,b,i") * t2_abab_vvoo("a,b,k,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += 0.500000 <j,k||i,b>_abab t2_abab(a,b,j,k) m2_aaa(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") -= 0.500000 * V_blks_["abba_oovo"]("j,k,b,i") * t2_abab_vvoo("a,b,j,k") * m2_xaaa_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += 0.500000 <k,j||b,i>_bbbb t2_bbbb(b,a,k,j) m2_bba(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") += 0.500000 * V_blks_["bbbb_oovo"]("k,j,b,i") * t2_bbbb_vvoo("b,a,k,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += 0.250000 <k,j||e,i>_aaaa t2_aaaa(b,a,k,j) m2_aaa(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") += 0.250000 * V_blks_["aaaa_oovo"]("k,j,e,i") * t2_aaaa_vvoo("b,a,k,j") * m2_xaaa_Lovv("I,i,b,a");

            // sigmam1_xa_Lv += -0.250000 <k,j||e,i>_abab t2_abab(a,b,k,j) m2_bba(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") -= 0.250000 * V_blks_["abab_oovo"]("k,j,e,i") * t2_abab_vvoo("a,b,k,j") * m2_xbba_Lovv("I,i,b,a");

            // sigmam1_xa_Lv += -0.250000 <j,k||e,i>_abab t2_abab(a,b,j,k) m2_bba(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") -= 0.250000 * V_blks_["abab_oovo"]("j,k,e,i") * t2_abab_vvoo("a,b,j,k") * m2_xbba_Lovv("I,i,b,a");

            // sigmam1_xa_Lv += 0.500000 <j,a||b,c>_aaaa t2_aaaa(b,c,i,j) m2_aaa(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") -= 0.500000 * V_blks_["aaaa_vovv"]("a,j,b,c") * t2_aaaa_vvoo("b,c,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += -0.500000 <a,j||b,c>_abab t2_abab(b,c,i,j) m2_aaa(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_vovv"]("a,j,b,c") * t2_abab_vvoo("b,c,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += -0.500000 <j,a||b,c>_abab t2_abab(b,c,j,i) m2_bba(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") += 0.500000 * V_blks_["baab_vovv"]("a,j,b,c") * t2_abab_vvoo("b,c,j,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += -0.500000 <a,j||c,b>_abab t2_abab(c,b,i,j) m2_aaa(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") -= 0.500000 * V_blks_["abab_vovv"]("a,j,c,b") * t2_abab_vvoo("c,b,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += -0.500000 <j,a||c,b>_abab t2_abab(c,b,j,i) m2_bba(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") += 0.500000 * V_blks_["baab_vovv"]("a,j,c,b") * t2_abab_vvoo("c,b,j,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += 0.500000 <j,a||b,c>_bbbb t2_bbbb(b,c,i,j) m2_bba(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xa_Lv("I,e") -= 0.500000 * V_blks_["bbbb_vovv"]("a,j,b,c") * t2_bbbb_vvoo("b,c,i,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmam1_xa_Lv += -1.000000 <j,b||c,e>_aaaa t2_aaaa(c,a,i,j) m2_aaa(i,b,a) 
            // flops: o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") += V_blks_["aaaa_vovv"]("b,j,c,e") * t2_aaaa_vvoo("c,a,i,j") * m2_xaaa_Lovv("I,i,b,a");

            // sigmam1_xa_Lv += 1.000000 <b,j||e,c>_abab t2_abab(a,c,i,j) m2_aaa(i,b,a) 
            // flops: o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") += V_blks_["abab_vovv"]("b,j,e,c") * t2_abab_vvoo("a,c,i,j") * m2_xaaa_Lovv("I,i,b,a");

            // sigmam1_xa_Lv += 1.000000 <j,b||e,c>_abab t2_abab(a,c,j,i) m2_bba(i,b,a) 
            // flops: o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xa_Lv("I,e") -= V_blks_["baab_vovv"]("b,j,e,c") * t2_abab_vvoo("a,c,j,i") * m2_xbba_Lovv("I,i,b,a");
        }




/// ****** sigmam1_b(e) ****** ///



        if (include_u1_) {

            // sigmam1_xb_Lv += 1.000000 m1_b(e) w0 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xb_Lv("I,e") += m1_xb_Lv("I,e") * w0;

            // sigmam1_xb_Lv += 1.000000 f_aa(i,i) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xb_Lv("I,e") += scalar5 * m1_xb_Lv("I,e");

            // sigmam1_xb_Lv += 1.000000 f_bb(i,i) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xb_Lv("I,e") += scalar6 * m1_xb_Lv("I,e");

            // sigmam1_xb_Lv += 1.000000 f_bb(a,e) m1_b(a) 
            // flops: o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") += F_blks_["bb_vv"]("a,e") * m1_xb_Lv("I,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmam1_xb_Lv += -1.000000 f_aa(a,i) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") -= F_blks_["aa_vo"]("a,i") * m2_xaab_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += -1.000000 f_bb(a,i) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") -= F_blks_["bb_vo"]("a,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += 1.000000 f_aa(j,b) t2_aaaa(b,a,i,j) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") += F_blks_["aa_ov"]("j,b") * t2_aaaa_vvoo("b,a,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += -1.000000 f_aa(j,b) t2_abab(b,a,j,i) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") -= F_blks_["aa_ov"]("j,b") * t2_abab_vvoo("b,a,j,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += -1.000000 f_bb(j,b) t2_abab(a,b,i,j) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") -= F_blks_["bb_ov"]("j,b") * t2_abab_vvoo("a,b,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += 1.000000 f_bb(j,b) t2_bbbb(b,a,i,j) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") += F_blks_["bb_ov"]("j,b") * t2_bbbb_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += 0.500000 f_bb(j,e) t2_abab(b,a,i,j) m2_aab(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") += 0.500000 * F_blks_["bb_ov"]("j,e") * t2_abab_vvoo("b,a,i,j") * m2_xaab_Lovv("I,i,b,a");

            // sigmam1_xb_Lv += 0.500000 f_bb(j,e) t2_bbbb(b,a,i,j) m2_bbb(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") += 0.500000 * F_blks_["bb_ov"]("j,e") * t2_bbbb_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,i,b,a");
        }

        if (include_u1_) {

            // sigmam1_xb_Lv += -1.000000 d-_aa(i,i) l1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xb_Lv("I,e") -= scalar7 * l1_xb_Lv("I,e");

            // sigmam1_xb_Lv += -1.000000 d-_bb(i,i) l1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xb_Lv("I,e") -= scalar8 * l1_xb_Lv("I,e");

            // sigmam1_xb_Lv += -1.000000 d-_bb(a,e) l1_b(a) 
            // flops: o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") -= dp_bb_vv("a,e") * l1_xb_Lv("I,a");

            // sigmam1_xb_Lv += 1.000000 d-_aa(a,i) l2_aab(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") += dp_aa_vo("a,i") * l2_xaab_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += 1.000000 d-_bb(a,i) l2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") += dp_bb_vo("a,i") * l2_xbbb_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += -1.000000 d-_aa(j,b) l2_aab(i,a,e) t2_aaaa(b,a,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") -= dp_aa_ov("j,b") * l2_xaab_Lovv("I,i,a,e") * t2_aaaa_vvoo("b,a,i,j");

            // sigmam1_xb_Lv += 1.000000 d-_bb(j,b) l2_aab(i,a,e) t2_abab(a,b,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") += dp_bb_ov("j,b") * l2_xaab_Lovv("I,i,a,e") * t2_abab_vvoo("a,b,i,j");

            // sigmam1_xb_Lv += 1.000000 d-_aa(j,b) l2_bbb(i,a,e) t2_abab(b,a,j,i) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") += dp_aa_ov("j,b") * l2_xbbb_Lovv("I,i,a,e") * t2_abab_vvoo("b,a,j,i");

            // sigmam1_xb_Lv += -1.000000 d-_bb(j,b) l2_bbb(i,a,e) t2_bbbb(b,a,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") -= dp_bb_ov("j,b") * l2_xbbb_Lovv("I,i,a,e") * t2_bbbb_vvoo("b,a,i,j");

            // sigmam1_xb_Lv += -0.500000 d-_bb(j,e) l2_aab(i,b,a) t2_abab(b,a,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") -= 0.500000 * dp_bb_ov("j,e") * l2_xaab_Lovv("I,i,b,a") * t2_abab_vvoo("b,a,i,j");

            // sigmam1_xb_Lv += -0.500000 d-_bb(j,e) l2_bbb(i,b,a) t2_bbbb(b,a,i,j) 
            // flops: o2v3L1: 2, o0v1L1: 1 | mem: o2v3L1: 1, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") -= 0.500000 * dp_bb_ov("j,e") * l2_xbbb_Lovv("I,i,b,a") * t2_bbbb_vvoo("b,a,i,j");
        }

        if (include_u0_ && include_u1_) {

            // sigmam1_xb_Lv += -1.000000 d-_aa(i,i) u0 m1_b(e) 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmam1_xb_Lv("I,e") -= u0 * scalar7 * m1_xb_Lv("I,e");

            // sigmam1_xb_Lv += -1.000000 d-_bb(i,i) u0 m1_b(e) 
            // flops: o0v1L1: 2, o0v0: 1 | mem: o0v1L1: 1, o0v1: 1, o0v0: 1, 
            sigmam1_xb_Lv("I,e") -= u0 * scalar8 * m1_xb_Lv("I,e");

            // sigmam1_xb_Lv += -1.000000 d-_bb(a,e) u0 m1_b(a) 
            // flops: o0v2L1: 1, o0v1L1: 1, o0v2: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmam1_xb_Lv("I,e") -= dp_bb_vv("a,e") * u0 * m1_xb_Lv("I,a");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmam1_xb_Lv += 1.000000 d-_aa(a,i) u0 m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") += dp_aa_vo("a,i") * u0 * m2_xaab_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += 1.000000 d-_bb(a,i) u0 m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") += dp_bb_vo("a,i") * u0 * m2_xbbb_Lovv("I,i,a,e");
        }

        if (include_u1_) {

            // sigmam1_xb_Lv += -2.000000 d-_aa(i,a) u1_aa(a,i) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xb_Lv("I,e") -= 2.000000 * scalar0 * m1_xb_Lv("I,e");

            // sigmam1_xb_Lv += -2.000000 d-_bb(i,a) u1_bb(a,i) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xb_Lv("I,e") -= 2.000000 * scalar1 * m1_xb_Lv("I,e");

            // sigmam1_xb_Lv += 2.000000 d-_bb(i,e) u1_bb(a,i) m1_b(a) 
            // flops: o0v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmam1_xb_Lv("I,e") += 2.000000 * dp_bb_ov("i,e") * u1_bb_vo("a,i") * m1_xb_Lv("I,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmam1_xb_Lv += 1.000000 d-_aa(i,i) u1_aa(a,j) m2_aab(j,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") += u1_aa_vo("a,j") * scalar7 * m2_xaab_Lovv("I,j,a,e");

            // sigmam1_xb_Lv += 1.000000 d-_bb(i,i) u1_aa(a,j) m2_aab(j,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") += u1_aa_vo("a,j") * scalar8 * m2_xaab_Lovv("I,j,a,e");

            // sigmam1_xb_Lv += 1.000000 d-_aa(i,i) u1_bb(a,j) m2_bbb(j,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") += u1_bb_vo("a,j") * scalar7 * m2_xbbb_Lovv("I,j,a,e");

            // sigmam1_xb_Lv += 1.000000 d-_bb(i,i) u1_bb(a,j) m2_bbb(j,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") += u1_bb_vo("a,j") * scalar8 * m2_xbbb_Lovv("I,j,a,e");

            // sigmam1_xb_Lv += -2.000000 d-_aa(j,i) u1_aa(a,j) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") -= 2.000000 * dp_aa_oo("j,i") * u1_aa_vo("a,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += -2.000000 d-_bb(j,i) u1_bb(a,j) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o2v1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") -= 2.000000 * dp_bb_oo("j,i") * u1_bb_vo("a,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += 2.000000 d-_aa(a,b) u1_aa(b,i) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") += 2.000000 * dp_aa_vv("a,b") * u1_aa_vo("b,i") * m2_xaab_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += 2.000000 d-_bb(a,b) u1_bb(b,i) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o1v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") += 2.000000 * dp_bb_vv("a,b") * u1_bb_vo("b,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += -1.000000 d-_bb(b,e) u1_bb(a,i) m2_bbb(i,b,a) 
            // flops: o1v3L1: 1, o1v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") -= dp_bb_vv("b,e") * u1_bb_vo("a,i") * m2_xbbb_Lovv("I,i,b,a");

            // sigmam1_xb_Lv += -2.000000 d-_aa(j,b) u2_aaaa(b,a,i,j) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") -= 2.000000 * dp_aa_ov("j,b") * u2_aaaa_vvoo("b,a,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += 2.000000 d-_aa(j,b) u2_abab(b,a,j,i) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") += 2.000000 * dp_aa_ov("j,b") * u2_abab_vvoo("b,a,j,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += 2.000000 d-_bb(j,b) u2_abab(a,b,i,j) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") += 2.000000 * dp_bb_ov("j,b") * u2_abab_vvoo("a,b,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += -2.000000 d-_bb(j,b) u2_bbbb(b,a,i,j) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") -= 2.000000 * dp_bb_ov("j,b") * u2_bbbb_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += -1.000000 d-_bb(j,e) u2_abab(b,a,i,j) m2_aab(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") -= dp_bb_ov("j,e") * u2_abab_vvoo("b,a,i,j") * m2_xaab_Lovv("I,i,b,a");

            // sigmam1_xb_Lv += -1.000000 d-_bb(j,e) u2_bbbb(b,a,i,j) m2_bbb(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") -= dp_bb_ov("j,e") * u2_bbbb_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,i,b,a");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmam1_xb_Lv += -1.000000 d-_aa(j,b) t2_aaaa(b,a,i,j) u0 m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 2, 
            sigmam1_xb_Lv("I,e") -= dp_aa_ov("j,b") * t2_aaaa_vvoo("b,a,i,j") * u0 * m2_xaab_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += 1.000000 d-_aa(j,b) t2_abab(b,a,j,i) u0 m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 2, 
            sigmam1_xb_Lv("I,e") += dp_aa_ov("j,b") * t2_abab_vvoo("b,a,j,i") * u0 * m2_xbbb_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += 1.000000 d-_bb(j,b) t2_abab(a,b,i,j) u0 m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 2, 
            sigmam1_xb_Lv("I,e") += dp_bb_ov("j,b") * t2_abab_vvoo("a,b,i,j") * u0 * m2_xaab_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += -1.000000 d-_bb(j,b) t2_bbbb(b,a,i,j) u0 m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o2v2: 1, o0v1L1: 1, o1v1: 1 | mem: o0v1L1: 2, o1v1: 2, 
            sigmam1_xb_Lv("I,e") -= dp_bb_ov("j,b") * t2_bbbb_vvoo("b,a,i,j") * u0 * m2_xbbb_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += -0.500000 d-_bb(j,e) t2_abab(b,a,i,j) u0 m2_aab(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v3: 1, o0v1L1: 1 | mem: o1v3: 2, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") -= 0.500000 * dp_bb_ov("j,e") * t2_abab_vvoo("b,a,i,j") * u0 * m2_xaab_Lovv("I,i,b,a");

            // sigmam1_xb_Lv += -0.500000 d-_bb(j,e) t2_bbbb(b,a,i,j) u0 m2_bbb(i,b,a) 
            // flops: o1v3L1: 1, o2v3: 1, o1v3: 1, o0v1L1: 1 | mem: o1v3: 2, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") -= 0.500000 * dp_bb_ov("j,e") * t2_bbbb_vvoo("b,a,i,j") * u0 * m2_xbbb_Lovv("I,i,b,a");
        }

        if (include_u1_) {

            // sigmam1_xb_Lv += -0.500000 <j,i||j,i>_aaaa m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xb_Lv("I,e") -= 0.500000 * scalar9 * m1_xb_Lv("I,e");

            // sigmam1_xb_Lv += -0.500000 <j,i||j,i>_abab m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xb_Lv("I,e") -= 0.500000 * scalar10 * m1_xb_Lv("I,e");

            // sigmam1_xb_Lv += -0.500000 <i,j||i,j>_abab m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xb_Lv("I,e") -= 0.500000 * scalar10 * m1_xb_Lv("I,e");

            // sigmam1_xb_Lv += -0.500000 <j,i||j,i>_bbbb m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xb_Lv("I,e") -= 0.500000 * scalar11 * m1_xb_Lv("I,e");
        }

        if (include_u1_ && include_u2_) {

            // sigmam1_xb_Lv += -0.500000 <b,a||i,e>_abab m2_aab(i,b,a) 
            // flops: o1v3L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") += 0.500000 * V_blks_["abba_vvvo"]("b,a,e,i") * m2_xaab_Lovv("I,i,b,a");

            // sigmam1_xb_Lv += 0.500000 <b,a||e,i>_bbbb m2_bbb(i,b,a) 
            // flops: o1v3L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") += 0.500000 * V_blks_["bbbb_vvvo"]("b,a,e,i") * m2_xbbb_Lovv("I,i,b,a");
        }

        if (include_u1_) {

            // sigmam1_xb_Lv += 0.250000 <j,i||a,b>_aaaa t2_aaaa(a,b,j,i) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xb_Lv("I,e") += 0.250000 * scalar2 * m1_xb_Lv("I,e");

            // sigmam1_xb_Lv += 0.250000 <j,i||a,b>_abab t2_abab(a,b,j,i) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xb_Lv("I,e") += 0.250000 * scalar3 * m1_xb_Lv("I,e");

            // sigmam1_xb_Lv += 0.250000 <i,j||a,b>_abab t2_abab(a,b,i,j) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xb_Lv("I,e") += 0.250000 * scalar3 * m1_xb_Lv("I,e");

            // sigmam1_xb_Lv += 0.250000 <j,i||b,a>_abab t2_abab(b,a,j,i) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xb_Lv("I,e") += 0.250000 * scalar3 * m1_xb_Lv("I,e");

            // sigmam1_xb_Lv += 0.250000 <i,j||b,a>_abab t2_abab(b,a,i,j) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xb_Lv("I,e") += 0.250000 * scalar3 * m1_xb_Lv("I,e");

            // sigmam1_xb_Lv += 0.250000 <j,i||a,b>_bbbb t2_bbbb(a,b,j,i) m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            sigmam1_xb_Lv("I,e") += 0.250000 * scalar4 * m1_xb_Lv("I,e");

            // sigmam1_xb_Lv += -0.500000 <j,i||b,e>_abab t2_abab(b,a,j,i) m1_b(a) 
            // flops: o2v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmam1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("j,i,b,e") * t2_abab_vvoo("b,a,j,i") * m1_xb_Lv("I,a");

            // sigmam1_xb_Lv += -0.500000 <i,j||b,e>_abab t2_abab(b,a,i,j) m1_b(a) 
            // flops: o2v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmam1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_oovv"]("i,j,b,e") * t2_abab_vvoo("b,a,i,j") * m1_xb_Lv("I,a");

            // sigmam1_xb_Lv += -0.500000 <j,i||b,e>_bbbb t2_bbbb(b,a,j,i) m1_b(a) 
            // flops: o2v3: 1, o0v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o0v2: 1, 
            sigmam1_xb_Lv("I,e") -= 0.500000 * V_blks_["bbbb_oovv"]("j,i,b,e") * t2_bbbb_vvoo("b,a,j,i") * m1_xb_Lv("I,a");
        }

        if (include_u1_ && include_u2_) {

            // sigmam1_xb_Lv += 0.500000 <k,j||b,i>_aaaa t2_aaaa(b,a,k,j) m2_aab(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") += 0.500000 * V_blks_["aaaa_oovo"]("k,j,b,i") * t2_aaaa_vvoo("b,a,k,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += 0.500000 <k,j||b,i>_abab t2_abab(b,a,k,j) m2_bbb(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") += 0.500000 * V_blks_["abab_oovo"]("k,j,b,i") * t2_abab_vvoo("b,a,k,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += 0.500000 <j,k||b,i>_abab t2_abab(b,a,j,k) m2_bbb(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") += 0.500000 * V_blks_["abab_oovo"]("j,k,b,i") * t2_abab_vvoo("b,a,j,k") * m2_xbbb_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += 0.500000 <k,j||i,b>_abab t2_abab(a,b,k,j) m2_aab(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") -= 0.500000 * V_blks_["abba_oovo"]("k,j,b,i") * t2_abab_vvoo("a,b,k,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += 0.500000 <j,k||i,b>_abab t2_abab(a,b,j,k) m2_aab(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") -= 0.500000 * V_blks_["abba_oovo"]("j,k,b,i") * t2_abab_vvoo("a,b,j,k") * m2_xaab_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += 0.500000 <k,j||b,i>_bbbb t2_bbbb(b,a,k,j) m2_bbb(i,a,e) 
            // flops: o3v2: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") += 0.500000 * V_blks_["bbbb_oovo"]("k,j,b,i") * t2_bbbb_vvoo("b,a,k,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += -0.250000 <k,j||i,e>_abab t2_abab(b,a,k,j) m2_aab(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") += 0.250000 * V_blks_["abba_oovo"]("k,j,e,i") * t2_abab_vvoo("b,a,k,j") * m2_xaab_Lovv("I,i,b,a");

            // sigmam1_xb_Lv += -0.250000 <j,k||i,e>_abab t2_abab(b,a,j,k) m2_aab(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") += 0.250000 * V_blks_["abba_oovo"]("j,k,e,i") * t2_abab_vvoo("b,a,j,k") * m2_xaab_Lovv("I,i,b,a");

            // sigmam1_xb_Lv += 0.250000 <k,j||e,i>_bbbb t2_bbbb(b,a,k,j) m2_bbb(i,b,a) 
            // flops: o3v3: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") += 0.250000 * V_blks_["bbbb_oovo"]("k,j,e,i") * t2_bbbb_vvoo("b,a,k,j") * m2_xbbb_Lovv("I,i,b,a");

            // sigmam1_xb_Lv += 0.500000 <j,a||b,c>_aaaa t2_aaaa(b,c,i,j) m2_aab(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") -= 0.500000 * V_blks_["aaaa_vovv"]("a,j,b,c") * t2_aaaa_vvoo("b,c,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += -0.500000 <a,j||b,c>_abab t2_abab(b,c,i,j) m2_aab(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_vovv"]("a,j,b,c") * t2_abab_vvoo("b,c,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += -0.500000 <j,a||b,c>_abab t2_abab(b,c,j,i) m2_bbb(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") += 0.500000 * V_blks_["baab_vovv"]("a,j,b,c") * t2_abab_vvoo("b,c,j,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += -0.500000 <a,j||c,b>_abab t2_abab(c,b,i,j) m2_aab(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") -= 0.500000 * V_blks_["abab_vovv"]("a,j,c,b") * t2_abab_vvoo("c,b,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += -0.500000 <j,a||c,b>_abab t2_abab(c,b,j,i) m2_bbb(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") += 0.500000 * V_blks_["baab_vovv"]("a,j,c,b") * t2_abab_vvoo("c,b,j,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += 0.500000 <j,a||b,c>_bbbb t2_bbbb(b,c,i,j) m2_bbb(i,a,e) 
            // flops: o2v3: 1, o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, o1v1: 1, 
            sigmam1_xb_Lv("I,e") -= 0.500000 * V_blks_["bbbb_vovv"]("a,j,b,c") * t2_bbbb_vvoo("b,c,i,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmam1_xb_Lv += 1.000000 <b,j||c,e>_abab t2_abab(c,a,i,j) m2_aab(i,b,a) 
            // flops: o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") += V_blks_["abab_vovv"]("b,j,c,e") * t2_abab_vvoo("c,a,i,j") * m2_xaab_Lovv("I,i,b,a");

            // sigmam1_xb_Lv += 1.000000 <j,b||c,e>_abab t2_abab(c,a,j,i) m2_bbb(i,b,a) 
            // flops: o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") -= V_blks_["baab_vovv"]("b,j,c,e") * t2_abab_vvoo("c,a,j,i") * m2_xbbb_Lovv("I,i,b,a");

            // sigmam1_xb_Lv += -1.000000 <j,b||c,e>_bbbb t2_bbbb(c,a,i,j) m2_bbb(i,b,a) 
            // flops: o2v4: 1, o1v3L1: 1, o0v1L1: 1 | mem: o1v3: 1, o0v1L1: 2, 
            sigmam1_xb_Lv("I,e") += V_blks_["bbbb_vovv"]("b,j,c,e") * t2_bbbb_vvoo("c,a,i,j") * m2_xbbb_Lovv("I,i,b,a");
        }




/// ****** sigmam2_aaa(e,f,n) ****** ///



        if (include_u2_) {

            // sigmam2_xaaa_Lvvo += 1.000000 m2_aaa(n,e,f) w0 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") += m2_xaaa_Lovv("I,n,e,f") * w0;
        }

        if (include_u1_ && include_u2_) {

            // sigmam2_xaaa_Lvvo += -1.000000 P(e,f) f_aa(n,e) m1_a(f) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = F_blks_["aa_ov"]("n,e") * m1_xa_Lv("I,f");
            sigmam2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmam2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmam2_xaaa_Lvvo += 1.000000 f_aa(i,i) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") += scalar5 * m2_xaaa_Lovv("I,n,e,f");

            // sigmam2_xaaa_Lvvo += 1.000000 f_bb(i,i) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") += scalar6 * m2_xaaa_Lovv("I,n,e,f");

            // sigmam2_xaaa_Lvvo += -1.000000 f_aa(n,i) m2_aaa(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xaaa_Lvvo("I,e,f,n") -= F_blks_["aa_oo"]("n,i") * m2_xaaa_Lovv("I,i,e,f");

            // sigmam2_xaaa_Lvvo += 1.000000 P(e,f) f_aa(a,e) m2_aaa(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = F_blks_["aa_vv"]("a,e") * m2_xaaa_Lovv("I,n,a,f");
            sigmam2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmam2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmam2_xaaa_Lvvo += 1.000000 P(e,f) d-_aa(n,e) l1_a(f) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("n,e") * l1_xa_Lv("I,f");
            sigmam2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmam2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmam2_xaaa_Lvvo += -1.000000 d-_aa(i,i) l2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") -= scalar7 * l2_xaaa_Lovv("I,n,e,f");

            // sigmam2_xaaa_Lvvo += -1.000000 d-_bb(i,i) l2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") -= scalar8 * l2_xaaa_Lovv("I,n,e,f");

            // sigmam2_xaaa_Lvvo += 1.000000 d-_aa(n,i) l2_aaa(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xaaa_Lvvo("I,e,f,n") += dp_aa_oo("n,i") * l2_xaaa_Lovv("I,i,e,f");

            // sigmam2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(a,e) l2_aaa(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_vv("a,e") * l2_xaaa_Lovv("I,n,a,f");
            sigmam2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmam2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmam2_xaaa_Lvvo += 1.000000 P(e,f) d-_aa(n,e) u0 m1_a(f) 
            // flops: o1v2L1: 3, o1v1: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("n,e") * u0 * m1_xa_Lv("I,f");
            sigmam2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmam2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u0_ && include_u2_) {

            // sigmam2_xaaa_Lvvo += -1.000000 d-_aa(i,i) u0 m2_aaa(n,e,f) 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") -= u0 * scalar7 * m2_xaaa_Lovv("I,n,e,f");

            // sigmam2_xaaa_Lvvo += -1.000000 d-_bb(i,i) u0 m2_aaa(n,e,f) 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") -= u0 * scalar8 * m2_xaaa_Lovv("I,n,e,f");

            // sigmam2_xaaa_Lvvo += 1.000000 d-_aa(n,i) u0 m2_aaa(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1, o2v0: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") += dp_aa_oo("n,i") * u0 * m2_xaaa_Lovv("I,i,e,f");

            // sigmam2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(a,e) u0 m2_aaa(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2, o0v2: 1 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_vv("a,e") * u0 * m2_xaaa_Lovv("I,n,a,f");
            sigmam2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmam2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u1_ && include_u2_) {

            // sigmam2_xaaa_Lvvo += -2.000000 d-_aa(i,a) u1_aa(a,i) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") -= 2.000000 * scalar0 * m2_xaaa_Lovv("I,n,e,f");

            // sigmam2_xaaa_Lvvo += -2.000000 d-_bb(i,a) u1_bb(a,i) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") -= 2.000000 * scalar1 * m2_xaaa_Lovv("I,n,e,f");

            // sigmam2_xaaa_Lvvo += 2.000000 d-_aa(n,a) u1_aa(a,i) m2_aaa(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1, o2v1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") += 2.000000 * dp_aa_ov("n,a") * u1_aa_vo("a,i") * m2_xaaa_Lovv("I,i,e,f");

            // sigmam2_xaaa_Lvvo += 2.000000 P(e,f) d-_aa(i,e) u1_aa(a,i) m2_aaa(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2, o1v2: 1 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 2.000000 * dp_aa_ov("i,e") * u1_aa_vo("a,i") * m2_xaaa_Lovv("I,n,a,f");
            sigmam2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmam2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmam2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(n,e) u1_aa(a,i) m2_aaa(i,a,f) 
            // flops: o2v3L1: 1, o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("n,e") * u1_aa_vo("a,i") * m2_xaaa_Lovv("I,i,a,f");
            sigmam2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmam2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmam2_xaaa_Lvvo += -1.000000 P(e,f) d-_aa(n,e) u1_bb(a,i) m2_bba(i,a,f) 
            // flops: o2v3L1: 1, o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = dp_aa_ov("n,e") * u1_bb_vo("a,i") * m2_xbba_Lovv("I,i,a,f");
            sigmam2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmam2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmam2_xaaa_Lvvo += -0.500000 <j,i||j,i>_aaaa m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * scalar9 * m2_xaaa_Lovv("I,n,e,f");

            // sigmam2_xaaa_Lvvo += -0.500000 <j,i||j,i>_abab m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * scalar10 * m2_xaaa_Lovv("I,n,e,f");

            // sigmam2_xaaa_Lvvo += -0.500000 <i,j||i,j>_abab m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * scalar10 * m2_xaaa_Lovv("I,n,e,f");

            // sigmam2_xaaa_Lvvo += -0.500000 <j,i||j,i>_bbbb m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * scalar11 * m2_xaaa_Lovv("I,n,e,f");
        }

        if (include_u1_ && include_u2_) {

            // sigmam2_xaaa_Lvvo += -1.000000 <n,a||e,f>_aaaa m1_a(a) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xaaa_Lvvo("I,e,f,n") += V_blks_["aaaa_vovv"]("a,n,e,f") * m1_xa_Lv("I,a");
        }

        if (include_u2_) {

            // sigmam2_xaaa_Lvvo += 1.000000 P(e,f) <n,a||e,i>_aaaa m2_aaa(i,a,f) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_vovo"]("a,n,e,i") * m2_xaaa_Lovv("I,i,a,f");
            sigmam2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmam2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmam2_xaaa_Lvvo += 1.000000 P(e,f) <n,a||e,i>_abab m2_bba(i,a,f) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["baab_vovo"]("a,n,e,i") * m2_xbba_Lovv("I,i,a,f");
            sigmam2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmam2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmam2_xaaa_Lvvo += 0.500000 <b,a||e,f>_aaaa m2_aaa(n,b,a) 
            // flops: o1v4L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xaaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["aaaa_vvvv"]("b,a,e,f") * m2_xaaa_Lovv("I,n,b,a");

            // sigmam2_xaaa_Lvvo += 0.250000 <j,i||a,b>_aaaa t2_aaaa(a,b,j,i) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar2 * m2_xaaa_Lovv("I,n,e,f");

            // sigmam2_xaaa_Lvvo += 0.250000 <j,i||a,b>_abab t2_abab(a,b,j,i) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar3 * m2_xaaa_Lovv("I,n,e,f");

            // sigmam2_xaaa_Lvvo += 0.250000 <i,j||a,b>_abab t2_abab(a,b,i,j) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar3 * m2_xaaa_Lovv("I,n,e,f");

            // sigmam2_xaaa_Lvvo += 0.250000 <j,i||b,a>_abab t2_abab(b,a,j,i) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar3 * m2_xaaa_Lovv("I,n,e,f");

            // sigmam2_xaaa_Lvvo += 0.250000 <i,j||b,a>_abab t2_abab(b,a,i,j) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar3 * m2_xaaa_Lovv("I,n,e,f");

            // sigmam2_xaaa_Lvvo += 0.250000 <j,i||a,b>_bbbb t2_bbbb(a,b,j,i) m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") += 0.250000 * scalar4 * m2_xaaa_Lovv("I,n,e,f");

            // sigmam2_xaaa_Lvvo += -0.500000 <n,j||a,b>_aaaa t2_aaaa(a,b,i,j) m2_aaa(i,e,f) 
            // flops: o2v2L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovv"]("n,j,a,b") * t2_aaaa_vvoo("a,b,i,j") * m2_xaaa_Lovv("I,i,e,f");

            // sigmam2_xaaa_Lvvo += -0.500000 <n,j||a,b>_abab t2_abab(a,b,i,j) m2_aaa(i,e,f) 
            // flops: o2v2L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("n,j,a,b") * t2_abab_vvoo("a,b,i,j") * m2_xaaa_Lovv("I,i,e,f");

            // sigmam2_xaaa_Lvvo += -0.500000 <n,j||b,a>_abab t2_abab(b,a,i,j) m2_aaa(i,e,f) 
            // flops: o2v2L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmam2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("n,j,b,a") * t2_abab_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,e,f");

            // sigmam2_xaaa_Lvvo += -0.500000 P(e,f) <j,i||b,e>_aaaa t2_aaaa(b,a,j,i) m2_aaa(n,a,f) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["aaaa_oovv"]("j,i,b,e") * t2_aaaa_vvoo("b,a,j,i") * m2_xaaa_Lovv("I,n,a,f");
            sigmam2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmam2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmam2_xaaa_Lvvo += -0.500000 P(e,f) <j,i||e,b>_abab t2_abab(a,b,j,i) m2_aaa(n,a,f) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("j,i,e,b") * t2_abab_vvoo("a,b,j,i") * m2_xaaa_Lovv("I,n,a,f");
            sigmam2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmam2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmam2_xaaa_Lvvo += -0.500000 P(e,f) <i,j||e,b>_abab t2_abab(a,b,i,j) m2_aaa(n,a,f) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("i,j,e,b") * t2_abab_vvoo("a,b,i,j") * m2_xaaa_Lovv("I,n,a,f");
            sigmam2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmam2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmam2_xaaa_Lvvo += 1.000000 P(e,f) <n,j||b,e>_aaaa t2_aaaa(b,a,i,j) m2_aaa(i,a,f) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_oovv"]("n,j,b,e") * t2_aaaa_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,a,f");
            sigmam2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmam2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmam2_xaaa_Lvvo += -1.000000 P(e,f) <n,j||b,e>_aaaa t2_abab(b,a,j,i) m2_bba(i,a,f) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["aaaa_oovv"]("n,j,b,e") * t2_abab_vvoo("b,a,j,i") * m2_xbba_Lovv("I,i,a,f");
            sigmam2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmam2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmam2_xaaa_Lvvo += 1.000000 P(e,f) <n,j||e,b>_abab t2_abab(a,b,i,j) m2_aaa(i,a,f) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("n,j,e,b") * t2_abab_vvoo("a,b,i,j") * m2_xaaa_Lovv("I,i,a,f");
            sigmam2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmam2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmam2_xaaa_Lvvo += -1.000000 P(e,f) <n,j||e,b>_abab t2_bbbb(b,a,i,j) m2_bba(i,a,f) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("n,j,e,b") * t2_bbbb_vvoo("b,a,i,j") * m2_xbba_Lovv("I,i,a,f");
            sigmam2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            sigmam2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");

            // sigmam2_xaaa_Lvvo += 0.250000 <j,i||e,f>_aaaa t2_aaaa(b,a,j,i) m2_aaa(n,b,a) 
            // flops: o1v4L1: 1, o2v4: 1, o1v2L1: 1 | mem: o0v4: 1, o1v2L1: 2, 
            sigmam2_xaaa_Lvvo("I,e,f,n") += 0.250000 * V_blks_["aaaa_oovv"]("j,i,e,f") * t2_aaaa_vvoo("b,a,j,i") * m2_xaaa_Lovv("I,n,b,a");

            // sigmam2_xaaa_Lvvo += -0.500000 <n,j||e,f>_aaaa t2_aaaa(b,a,i,j) m2_aaa(i,b,a) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmam2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovv"]("n,j,e,f") * t2_aaaa_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,b,a");

            // sigmam2_xaaa_Lvvo += -0.500000 <n,j||e,f>_aaaa t2_abab(a,b,j,i) m2_bba(i,b,a) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmam2_xaaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["aaaa_oovv"]("n,j,e,f") * t2_abab_vvoo("a,b,j,i") * m2_xbba_Lovv("I,i,b,a");
        }




/// ****** sigmam2_abb(e,f,n) ****** ///



        if (include_u1_ && include_u2_) {

            // sigmam2_xabb_Lvvo += 1.000000 f_bb(n,f) m1_a(e) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            sigmam2_xabb_Lvvo("I,e,f,n") += F_blks_["bb_ov"]("n,f") * m1_xa_Lv("I,e");
        }

        if (include_u2_) {

            // sigmam2_xabb_Lvvo += -1.000000 f_bb(a,f) m2_bba(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xabb_Lvvo("I,e,f,n") -= F_blks_["bb_vv"]("a,f") * m2_xbba_Lovv("I,n,a,e");

            // sigmam2_xabb_Lvvo += -1.000000 d-_bb(n,f) l1_a(e) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            sigmam2_xabb_Lvvo("I,e,f,n") -= dp_bb_ov("n,f") * l1_xa_Lv("I,e");

            // sigmam2_xabb_Lvvo += 1.000000 d-_bb(a,f) l2_bba(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xabb_Lvvo("I,e,f,n") += dp_bb_vv("a,f") * l2_xbba_Lovv("I,n,a,e");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmam2_xabb_Lvvo += -1.000000 d-_bb(n,f) u0 m1_a(e) 
            // flops: o1v2L1: 2, o1v1: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmam2_xabb_Lvvo("I,e,f,n") -= dp_bb_ov("n,f") * u0 * m1_xa_Lv("I,e");
        }

        if (include_u0_ && include_u2_) {

            // sigmam2_xabb_Lvvo += 1.000000 d-_bb(a,f) u0 m2_bba(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1, o0v2: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmam2_xabb_Lvvo("I,e,f,n") += dp_bb_vv("a,f") * u0 * m2_xbba_Lovv("I,n,a,e");
        }

        if (include_u1_ && include_u2_) {

            // sigmam2_xabb_Lvvo += -2.000000 d-_bb(i,f) u1_bb(a,i) m2_bba(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1, o1v2: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmam2_xabb_Lvvo("I,e,f,n") -= 2.000000 * dp_bb_ov("i,f") * u1_bb_vo("a,i") * m2_xbba_Lovv("I,n,a,e");

            // sigmam2_xabb_Lvvo += 1.000000 d-_bb(n,f) u1_aa(a,i) m2_aaa(i,a,e) 
            // flops: o2v3L1: 1, o1v2L1: 1, o2v2: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmam2_xabb_Lvvo("I,e,f,n") += dp_bb_ov("n,f") * u1_aa_vo("a,i") * m2_xaaa_Lovv("I,i,a,e");

            // sigmam2_xabb_Lvvo += 1.000000 d-_bb(n,f) u1_bb(a,i) m2_bba(i,a,e) 
            // flops: o2v3L1: 1, o1v2L1: 1, o2v2: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmam2_xabb_Lvvo("I,e,f,n") += dp_bb_ov("n,f") * u1_bb_vo("a,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmam2_xabb_Lvvo += 1.000000 <a,n||e,f>_abab m1_a(a) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xabb_Lvvo("I,e,f,n") += V_blks_["abab_vovv"]("a,n,e,f") * m1_xa_Lv("I,a");
        }

        if (include_u2_) {

            // sigmam2_xabb_Lvvo += -1.000000 <a,n||i,f>_abab m2_aaa(i,a,e) 
            // flops: o2v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xabb_Lvvo("I,e,f,n") += V_blks_["abba_vovo"]("a,n,f,i") * m2_xaaa_Lovv("I,i,a,e");

            // sigmam2_xabb_Lvvo += -1.000000 <n,a||f,i>_bbbb m2_bba(i,a,e) 
            // flops: o2v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xabb_Lvvo("I,e,f,n") += V_blks_["bbbb_vovo"]("a,n,f,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmam2_xabb_Lvvo += -0.500000 <a,b||e,f>_abab m2_bba(n,b,a) 
            // flops: o1v4L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xabb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_vvvv"]("a,b,e,f") * m2_xbba_Lovv("I,n,b,a");

            // sigmam2_xabb_Lvvo += 0.500000 <j,i||b,f>_abab t2_abab(b,a,j,i) m2_bba(n,a,e) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmam2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("j,i,b,f") * t2_abab_vvoo("b,a,j,i") * m2_xbba_Lovv("I,n,a,e");

            // sigmam2_xabb_Lvvo += 0.500000 <i,j||b,f>_abab t2_abab(b,a,i,j) m2_bba(n,a,e) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmam2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("i,j,b,f") * t2_abab_vvoo("b,a,i,j") * m2_xbba_Lovv("I,n,a,e");

            // sigmam2_xabb_Lvvo += 0.500000 <j,i||b,f>_bbbb t2_bbbb(b,a,j,i) m2_bba(n,a,e) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmam2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["bbbb_oovv"]("j,i,b,f") * t2_bbbb_vvoo("b,a,j,i") * m2_xbba_Lovv("I,n,a,e");

            // sigmam2_xabb_Lvvo += 1.000000 <j,n||b,f>_abab t2_aaaa(b,a,i,j) m2_aaa(i,a,e) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmam2_xabb_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("j,n,b,f") * t2_aaaa_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmam2_xabb_Lvvo += -1.000000 <j,n||b,f>_abab t2_abab(b,a,j,i) m2_bba(i,a,e) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmam2_xabb_Lvvo("I,e,f,n") -= V_blks_["abab_oovv"]("j,n,b,f") * t2_abab_vvoo("b,a,j,i") * m2_xbba_Lovv("I,i,a,e");

            // sigmam2_xabb_Lvvo += 1.000000 <n,j||b,f>_bbbb t2_abab(a,b,i,j) m2_aaa(i,a,e) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmam2_xabb_Lvvo("I,e,f,n") += V_blks_["bbbb_oovv"]("n,j,b,f") * t2_abab_vvoo("a,b,i,j") * m2_xaaa_Lovv("I,i,a,e");

            // sigmam2_xabb_Lvvo += -1.000000 <n,j||b,f>_bbbb t2_bbbb(b,a,i,j) m2_bba(i,a,e) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmam2_xabb_Lvvo("I,e,f,n") -= V_blks_["bbbb_oovv"]("n,j,b,f") * t2_bbbb_vvoo("b,a,i,j") * m2_xbba_Lovv("I,i,a,e");

            // sigmam2_xabb_Lvvo += -0.250000 <j,i||e,f>_abab t2_abab(a,b,j,i) m2_bba(n,b,a) 
            // flops: o1v4L1: 1, o2v4: 1, o1v2L1: 1 | mem: o0v4: 1, o1v2L1: 2, 
            sigmam2_xabb_Lvvo("I,e,f,n") -= 0.250000 * V_blks_["abab_oovv"]("j,i,e,f") * t2_abab_vvoo("a,b,j,i") * m2_xbba_Lovv("I,n,b,a");

            // sigmam2_xabb_Lvvo += -0.250000 <i,j||e,f>_abab t2_abab(a,b,i,j) m2_bba(n,b,a) 
            // flops: o1v4L1: 1, o2v4: 1, o1v2L1: 1 | mem: o0v4: 1, o1v2L1: 2, 
            sigmam2_xabb_Lvvo("I,e,f,n") -= 0.250000 * V_blks_["abab_oovv"]("i,j,e,f") * t2_abab_vvoo("a,b,i,j") * m2_xbba_Lovv("I,n,b,a");

            // sigmam2_xabb_Lvvo += 0.500000 <j,n||e,f>_abab t2_aaaa(b,a,i,j) m2_aaa(i,b,a) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmam2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("j,n,e,f") * t2_aaaa_vvoo("b,a,i,j") * m2_xaaa_Lovv("I,i,b,a");

            // sigmam2_xabb_Lvvo += 0.500000 <j,n||e,f>_abab t2_abab(a,b,j,i) m2_bba(i,b,a) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmam2_xabb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("j,n,e,f") * t2_abab_vvoo("a,b,j,i") * m2_xbba_Lovv("I,i,b,a");
        }




/// ****** sigmam2_baa(e,f,n) ****** ///



        if (include_u1_ && include_u2_) {

            // sigmam2_xbaa_Lvvo += 1.000000 f_aa(n,f) m1_b(e) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            sigmam2_xbaa_Lvvo("I,e,f,n") += F_blks_["aa_ov"]("n,f") * m1_xb_Lv("I,e");
        }

        if (include_u2_) {

            // sigmam2_xbaa_Lvvo += -1.000000 f_aa(a,f) m2_aab(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xbaa_Lvvo("I,e,f,n") -= F_blks_["aa_vv"]("a,f") * m2_xaab_Lovv("I,n,a,e");

            // sigmam2_xbaa_Lvvo += -1.000000 d-_aa(n,f) l1_b(e) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            sigmam2_xbaa_Lvvo("I,e,f,n") -= dp_aa_ov("n,f") * l1_xb_Lv("I,e");

            // sigmam2_xbaa_Lvvo += 1.000000 d-_aa(a,f) l2_aab(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xbaa_Lvvo("I,e,f,n") += dp_aa_vv("a,f") * l2_xaab_Lovv("I,n,a,e");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmam2_xbaa_Lvvo += -1.000000 d-_aa(n,f) u0 m1_b(e) 
            // flops: o1v2L1: 2, o1v1: 1 | mem: o1v2L1: 2, o1v1: 1, 
            sigmam2_xbaa_Lvvo("I,e,f,n") -= dp_aa_ov("n,f") * u0 * m1_xb_Lv("I,e");
        }

        if (include_u0_ && include_u2_) {

            // sigmam2_xbaa_Lvvo += 1.000000 d-_aa(a,f) u0 m2_aab(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1, o0v2: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmam2_xbaa_Lvvo("I,e,f,n") += dp_aa_vv("a,f") * u0 * m2_xaab_Lovv("I,n,a,e");
        }

        if (include_u1_ && include_u2_) {

            // sigmam2_xbaa_Lvvo += -2.000000 d-_aa(i,f) u1_aa(a,i) m2_aab(n,a,e) 
            // flops: o1v3L1: 1, o1v2L1: 1, o1v2: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmam2_xbaa_Lvvo("I,e,f,n") -= 2.000000 * dp_aa_ov("i,f") * u1_aa_vo("a,i") * m2_xaab_Lovv("I,n,a,e");

            // sigmam2_xbaa_Lvvo += 1.000000 d-_aa(n,f) u1_aa(a,i) m2_aab(i,a,e) 
            // flops: o2v3L1: 1, o1v2L1: 1, o2v2: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmam2_xbaa_Lvvo("I,e,f,n") += dp_aa_ov("n,f") * u1_aa_vo("a,i") * m2_xaab_Lovv("I,i,a,e");

            // sigmam2_xbaa_Lvvo += 1.000000 d-_aa(n,f) u1_bb(a,i) m2_bbb(i,a,e) 
            // flops: o2v3L1: 1, o1v2L1: 1, o2v2: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmam2_xbaa_Lvvo("I,e,f,n") += dp_aa_ov("n,f") * u1_bb_vo("a,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmam2_xbaa_Lvvo += 1.000000 <n,a||f,e>_abab m1_b(a) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xbaa_Lvvo("I,e,f,n") -= V_blks_["baab_vovv"]("a,n,f,e") * m1_xb_Lv("I,a");
        }

        if (include_u2_) {

            // sigmam2_xbaa_Lvvo += -1.000000 <n,a||f,i>_aaaa m2_aab(i,a,e) 
            // flops: o2v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xbaa_Lvvo("I,e,f,n") += V_blks_["aaaa_vovo"]("a,n,f,i") * m2_xaab_Lovv("I,i,a,e");

            // sigmam2_xbaa_Lvvo += -1.000000 <n,a||f,i>_abab m2_bbb(i,a,e) 
            // flops: o2v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xbaa_Lvvo("I,e,f,n") += V_blks_["baab_vovo"]("a,n,f,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmam2_xbaa_Lvvo += -0.500000 <b,a||f,e>_abab m2_aab(n,b,a) 
            // flops: o1v4L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xbaa_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_vvvv"]("b,a,f,e") * m2_xaab_Lovv("I,n,b,a");

            // sigmam2_xbaa_Lvvo += 0.500000 <j,i||b,f>_aaaa t2_aaaa(b,a,j,i) m2_aab(n,a,e) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmam2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["aaaa_oovv"]("j,i,b,f") * t2_aaaa_vvoo("b,a,j,i") * m2_xaab_Lovv("I,n,a,e");

            // sigmam2_xbaa_Lvvo += 0.500000 <j,i||f,b>_abab t2_abab(a,b,j,i) m2_aab(n,a,e) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmam2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("j,i,f,b") * t2_abab_vvoo("a,b,j,i") * m2_xaab_Lovv("I,n,a,e");

            // sigmam2_xbaa_Lvvo += 0.500000 <i,j||f,b>_abab t2_abab(a,b,i,j) m2_aab(n,a,e) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o0v2: 1, 
            sigmam2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("i,j,f,b") * t2_abab_vvoo("a,b,i,j") * m2_xaab_Lovv("I,n,a,e");

            // sigmam2_xbaa_Lvvo += -1.000000 <n,j||b,f>_aaaa t2_aaaa(b,a,i,j) m2_aab(i,a,e) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmam2_xbaa_Lvvo("I,e,f,n") -= V_blks_["aaaa_oovv"]("n,j,b,f") * t2_aaaa_vvoo("b,a,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmam2_xbaa_Lvvo += 1.000000 <n,j||b,f>_aaaa t2_abab(b,a,j,i) m2_bbb(i,a,e) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmam2_xbaa_Lvvo("I,e,f,n") += V_blks_["aaaa_oovv"]("n,j,b,f") * t2_abab_vvoo("b,a,j,i") * m2_xbbb_Lovv("I,i,a,e");

            // sigmam2_xbaa_Lvvo += -1.000000 <n,j||f,b>_abab t2_abab(a,b,i,j) m2_aab(i,a,e) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmam2_xbaa_Lvvo("I,e,f,n") -= V_blks_["abab_oovv"]("n,j,f,b") * t2_abab_vvoo("a,b,i,j") * m2_xaab_Lovv("I,i,a,e");

            // sigmam2_xbaa_Lvvo += 1.000000 <n,j||f,b>_abab t2_bbbb(b,a,i,j) m2_bbb(i,a,e) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v2: 1, 
            sigmam2_xbaa_Lvvo("I,e,f,n") += V_blks_["abab_oovv"]("n,j,f,b") * t2_bbbb_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,i,a,e");

            // sigmam2_xbaa_Lvvo += -0.250000 <j,i||f,e>_abab t2_abab(b,a,j,i) m2_aab(n,b,a) 
            // flops: o1v4L1: 1, o2v4: 1, o1v2L1: 1 | mem: o0v4: 1, o1v2L1: 2, 
            sigmam2_xbaa_Lvvo("I,e,f,n") -= 0.250000 * V_blks_["abab_oovv"]("j,i,f,e") * t2_abab_vvoo("b,a,j,i") * m2_xaab_Lovv("I,n,b,a");

            // sigmam2_xbaa_Lvvo += -0.250000 <i,j||f,e>_abab t2_abab(b,a,i,j) m2_aab(n,b,a) 
            // flops: o1v4L1: 1, o2v4: 1, o1v2L1: 1 | mem: o0v4: 1, o1v2L1: 2, 
            sigmam2_xbaa_Lvvo("I,e,f,n") -= 0.250000 * V_blks_["abab_oovv"]("i,j,f,e") * t2_abab_vvoo("b,a,i,j") * m2_xaab_Lovv("I,n,b,a");

            // sigmam2_xbaa_Lvvo += 0.500000 <n,j||f,e>_abab t2_abab(b,a,i,j) m2_aab(i,b,a) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmam2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("n,j,f,e") * t2_abab_vvoo("b,a,i,j") * m2_xaab_Lovv("I,i,b,a");

            // sigmam2_xbaa_Lvvo += 0.500000 <n,j||f,e>_abab t2_bbbb(b,a,i,j) m2_bbb(i,b,a) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmam2_xbaa_Lvvo("I,e,f,n") += 0.500000 * V_blks_["abab_oovv"]("n,j,f,e") * t2_bbbb_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,i,b,a");
        }




/// ****** sigmam2_bbb(e,f,n) ****** ///



        if (include_u2_) {

            // sigmam2_xbbb_Lvvo += 1.000000 m2_bbb(n,e,f) w0 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") += m2_xbbb_Lovv("I,n,e,f") * w0;
        }

        if (include_u1_ && include_u2_) {

            // sigmam2_xbbb_Lvvo += -1.000000 P(e,f) f_bb(n,e) m1_b(f) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = F_blks_["bb_ov"]("n,e") * m1_xb_Lv("I,f");
            sigmam2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmam2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmam2_xbbb_Lvvo += 1.000000 f_aa(i,i) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") += scalar5 * m2_xbbb_Lovv("I,n,e,f");

            // sigmam2_xbbb_Lvvo += 1.000000 f_bb(i,i) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") += scalar6 * m2_xbbb_Lovv("I,n,e,f");

            // sigmam2_xbbb_Lvvo += -1.000000 f_bb(n,i) m2_bbb(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xbbb_Lvvo("I,e,f,n") -= F_blks_["bb_oo"]("n,i") * m2_xbbb_Lovv("I,i,e,f");

            // sigmam2_xbbb_Lvvo += 1.000000 P(e,f) f_bb(a,e) m2_bbb(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = F_blks_["bb_vv"]("a,e") * m2_xbbb_Lovv("I,n,a,f");
            sigmam2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmam2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmam2_xbbb_Lvvo += 1.000000 P(e,f) d-_bb(n,e) l1_b(f) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("n,e") * l1_xb_Lv("I,f");
            sigmam2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmam2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmam2_xbbb_Lvvo += -1.000000 d-_aa(i,i) l2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") -= scalar7 * l2_xbbb_Lovv("I,n,e,f");

            // sigmam2_xbbb_Lvvo += -1.000000 d-_bb(i,i) l2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") -= scalar8 * l2_xbbb_Lovv("I,n,e,f");

            // sigmam2_xbbb_Lvvo += 1.000000 d-_bb(n,i) l2_bbb(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xbbb_Lvvo("I,e,f,n") += dp_bb_oo("n,i") * l2_xbbb_Lovv("I,i,e,f");

            // sigmam2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(a,e) l2_bbb(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_vv("a,e") * l2_xbbb_Lovv("I,n,a,f");
            sigmam2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmam2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // sigmam2_xbbb_Lvvo += 1.000000 P(e,f) d-_bb(n,e) u0 m1_b(f) 
            // flops: o1v2L1: 3, o1v1: 1 | mem: o1v2L1: 2, o1v1: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("n,e") * u0 * m1_xb_Lv("I,f");
            sigmam2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmam2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u0_ && include_u2_) {

            // sigmam2_xbbb_Lvvo += -1.000000 d-_aa(i,i) u0 m2_bbb(n,e,f) 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") -= u0 * scalar7 * m2_xbbb_Lovv("I,n,e,f");

            // sigmam2_xbbb_Lvvo += -1.000000 d-_bb(i,i) u0 m2_bbb(n,e,f) 
            // flops: o1v2L1: 2, o0v0: 1 | mem: o1v2L1: 1, o1v2: 1, o0v0: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") -= u0 * scalar8 * m2_xbbb_Lovv("I,n,e,f");

            // sigmam2_xbbb_Lvvo += 1.000000 d-_bb(n,i) u0 m2_bbb(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1, o2v0: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") += dp_bb_oo("n,i") * u0 * m2_xbbb_Lovv("I,i,e,f");

            // sigmam2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(a,e) u0 m2_bbb(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2, o0v2: 1 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_vv("a,e") * u0 * m2_xbbb_Lovv("I,n,a,f");
            sigmam2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmam2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u1_ && include_u2_) {

            // sigmam2_xbbb_Lvvo += -2.000000 d-_aa(i,a) u1_aa(a,i) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") -= 2.000000 * scalar0 * m2_xbbb_Lovv("I,n,e,f");

            // sigmam2_xbbb_Lvvo += -2.000000 d-_bb(i,a) u1_bb(a,i) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") -= 2.000000 * scalar1 * m2_xbbb_Lovv("I,n,e,f");

            // sigmam2_xbbb_Lvvo += 2.000000 d-_bb(n,a) u1_bb(a,i) m2_bbb(i,e,f) 
            // flops: o2v2L1: 1, o1v2L1: 1, o2v1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") += 2.000000 * dp_bb_ov("n,a") * u1_bb_vo("a,i") * m2_xbbb_Lovv("I,i,e,f");

            // sigmam2_xbbb_Lvvo += 2.000000 P(e,f) d-_bb(i,e) u1_bb(a,i) m2_bbb(n,a,f) 
            // flops: o1v3L1: 1, o1v2L1: 2, o1v2: 1 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 2.000000 * dp_bb_ov("i,e") * u1_bb_vo("a,i") * m2_xbbb_Lovv("I,n,a,f");
            sigmam2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmam2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmam2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(n,e) u1_aa(a,i) m2_aab(i,a,f) 
            // flops: o2v3L1: 1, o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("n,e") * u1_aa_vo("a,i") * m2_xaab_Lovv("I,i,a,f");
            sigmam2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmam2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmam2_xbbb_Lvvo += -1.000000 P(e,f) d-_bb(n,e) u1_bb(a,i) m2_bbb(i,a,f) 
            // flops: o2v3L1: 1, o1v2L1: 2, o2v2: 1 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = dp_bb_ov("n,e") * u1_bb_vo("a,i") * m2_xbbb_Lovv("I,i,a,f");
            sigmam2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmam2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");
        }

        if (include_u2_) {

            // sigmam2_xbbb_Lvvo += -0.500000 <j,i||j,i>_aaaa m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * scalar9 * m2_xbbb_Lovv("I,n,e,f");

            // sigmam2_xbbb_Lvvo += -0.500000 <j,i||j,i>_abab m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * scalar10 * m2_xbbb_Lovv("I,n,e,f");

            // sigmam2_xbbb_Lvvo += -0.500000 <i,j||i,j>_abab m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * scalar10 * m2_xbbb_Lovv("I,n,e,f");

            // sigmam2_xbbb_Lvvo += -0.500000 <j,i||j,i>_bbbb m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * scalar11 * m2_xbbb_Lovv("I,n,e,f");
        }

        if (include_u1_ && include_u2_) {

            // sigmam2_xbbb_Lvvo += -1.000000 <n,a||e,f>_bbbb m1_b(a) 
            // flops: o1v3L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xbbb_Lvvo("I,e,f,n") += V_blks_["bbbb_vovv"]("a,n,e,f") * m1_xb_Lv("I,a");
        }

        if (include_u2_) {

            // sigmam2_xbbb_Lvvo += 1.000000 P(e,f) <a,n||i,e>_abab m2_aab(i,a,f) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["abba_vovo"]("a,n,e,i") * m2_xaab_Lovv("I,i,a,f");
            sigmam2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmam2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmam2_xbbb_Lvvo += 1.000000 P(e,f) <n,a||e,i>_bbbb m2_bbb(i,a,f) 
            // flops: o2v3L1: 1, o1v2L1: 2 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_vovo"]("a,n,e,i") * m2_xbbb_Lovv("I,i,a,f");
            sigmam2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmam2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmam2_xbbb_Lvvo += 0.500000 <b,a||e,f>_bbbb m2_bbb(n,b,a) 
            // flops: o1v4L1: 1, o1v2L1: 1 | mem: o1v2L1: 2, 
            sigmam2_xbbb_Lvvo("I,e,f,n") += 0.500000 * V_blks_["bbbb_vvvv"]("b,a,e,f") * m2_xbbb_Lovv("I,n,b,a");

            // sigmam2_xbbb_Lvvo += 0.250000 <j,i||a,b>_aaaa t2_aaaa(a,b,j,i) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar2 * m2_xbbb_Lovv("I,n,e,f");

            // sigmam2_xbbb_Lvvo += 0.250000 <j,i||a,b>_abab t2_abab(a,b,j,i) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar3 * m2_xbbb_Lovv("I,n,e,f");

            // sigmam2_xbbb_Lvvo += 0.250000 <i,j||a,b>_abab t2_abab(a,b,i,j) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar3 * m2_xbbb_Lovv("I,n,e,f");

            // sigmam2_xbbb_Lvvo += 0.250000 <j,i||b,a>_abab t2_abab(b,a,j,i) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar3 * m2_xbbb_Lovv("I,n,e,f");

            // sigmam2_xbbb_Lvvo += 0.250000 <i,j||b,a>_abab t2_abab(b,a,i,j) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar3 * m2_xbbb_Lovv("I,n,e,f");

            // sigmam2_xbbb_Lvvo += 0.250000 <j,i||a,b>_bbbb t2_bbbb(a,b,j,i) m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") += 0.250000 * scalar4 * m2_xbbb_Lovv("I,n,e,f");

            // sigmam2_xbbb_Lvvo += -0.500000 <j,n||a,b>_abab t2_abab(a,b,j,i) m2_bbb(i,e,f) 
            // flops: o2v2L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("j,n,a,b") * t2_abab_vvoo("a,b,j,i") * m2_xbbb_Lovv("I,i,e,f");

            // sigmam2_xbbb_Lvvo += -0.500000 <j,n||b,a>_abab t2_abab(b,a,j,i) m2_bbb(i,e,f) 
            // flops: o2v2L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["abab_oovv"]("j,n,b,a") * t2_abab_vvoo("b,a,j,i") * m2_xbbb_Lovv("I,i,e,f");

            // sigmam2_xbbb_Lvvo += -0.500000 <n,j||a,b>_bbbb t2_bbbb(a,b,i,j) m2_bbb(i,e,f) 
            // flops: o2v2L1: 1, o3v2: 1, o1v2L1: 1 | mem: o1v2L1: 2, o2v0: 1, 
            sigmam2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovv"]("n,j,a,b") * t2_bbbb_vvoo("a,b,i,j") * m2_xbbb_Lovv("I,i,e,f");

            // sigmam2_xbbb_Lvvo += -0.500000 P(e,f) <j,i||b,e>_abab t2_abab(b,a,j,i) m2_bbb(n,a,f) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("j,i,b,e") * t2_abab_vvoo("b,a,j,i") * m2_xbbb_Lovv("I,n,a,f");
            sigmam2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmam2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmam2_xbbb_Lvvo += -0.500000 P(e,f) <i,j||b,e>_abab t2_abab(b,a,i,j) m2_bbb(n,a,f) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["abab_oovv"]("i,j,b,e") * t2_abab_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,n,a,f");
            sigmam2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmam2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmam2_xbbb_Lvvo += -0.500000 P(e,f) <j,i||b,e>_bbbb t2_bbbb(b,a,j,i) m2_bbb(n,a,f) 
            // flops: o1v3L1: 1, o2v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o0v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = 0.500000 * V_blks_["bbbb_oovv"]("j,i,b,e") * t2_bbbb_vvoo("b,a,j,i") * m2_xbbb_Lovv("I,n,a,f");
            sigmam2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmam2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmam2_xbbb_Lvvo += -1.000000 P(e,f) <j,n||b,e>_abab t2_aaaa(b,a,i,j) m2_aab(i,a,f) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("j,n,b,e") * t2_aaaa_vvoo("b,a,i,j") * m2_xaab_Lovv("I,i,a,f");
            sigmam2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmam2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmam2_xbbb_Lvvo += 1.000000 P(e,f) <j,n||b,e>_abab t2_abab(b,a,j,i) m2_bbb(i,a,f) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["abab_oovv"]("j,n,b,e") * t2_abab_vvoo("b,a,j,i") * m2_xbbb_Lovv("I,i,a,f");
            sigmam2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmam2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmam2_xbbb_Lvvo += -1.000000 P(e,f) <n,j||b,e>_bbbb t2_abab(a,b,i,j) m2_aab(i,a,f) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_oovv"]("n,j,b,e") * t2_abab_vvoo("a,b,i,j") * m2_xaab_Lovv("I,i,a,f");
            sigmam2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmam2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmam2_xbbb_Lvvo += 1.000000 P(e,f) <n,j||b,e>_bbbb t2_bbbb(b,a,i,j) m2_bbb(i,a,f) 
            // flops: o2v3L1: 1, o3v3: 1, o1v2L1: 2 | mem: o1v2L1: 2, o2v2: 1, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = V_blks_["bbbb_oovv"]("n,j,b,e") * t2_bbbb_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,i,a,f");
            sigmam2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,e,f,n");
            sigmam2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,f,e,n");

            // sigmam2_xbbb_Lvvo += 0.250000 <j,i||e,f>_bbbb t2_bbbb(b,a,j,i) m2_bbb(n,b,a) 
            // flops: o1v4L1: 1, o2v4: 1, o1v2L1: 1 | mem: o0v4: 1, o1v2L1: 2, 
            sigmam2_xbbb_Lvvo("I,e,f,n") += 0.250000 * V_blks_["bbbb_oovv"]("j,i,e,f") * t2_bbbb_vvoo("b,a,j,i") * m2_xbbb_Lovv("I,n,b,a");

            // sigmam2_xbbb_Lvvo += -0.500000 <n,j||e,f>_bbbb t2_abab(b,a,i,j) m2_aab(i,b,a) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmam2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovv"]("n,j,e,f") * t2_abab_vvoo("b,a,i,j") * m2_xaab_Lovv("I,i,b,a");

            // sigmam2_xbbb_Lvvo += -0.500000 <n,j||e,f>_bbbb t2_bbbb(b,a,i,j) m2_bbb(i,b,a) 
            // flops: o2v4L1: 1, o3v4: 1, o1v2L1: 1 | mem: o2v4: 1, o1v2L1: 2, 
            sigmam2_xbbb_Lvvo("I,e,f,n") -= 0.500000 * V_blks_["bbbb_oovv"]("n,j,e,f") * t2_bbbb_vvoo("b,a,i,j") * m2_xbbb_Lovv("I,i,b,a");
        }




/// ****** csigmar1_a(e) ****** ///



        if (include_u1_) {

            // csigmar1_xa_Lv += 1.000000 s1_a(e) 
            // flops: o0v1L1: 1 | mem: o0v1L1: 1, 
            csigmar1_xa_Lv("I,e") += s1_xa_Lv("I,e");
        }

        if (include_u0_) {

            // csigmar1_xa_Lv += 1.000000 r1_a(e) u0 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            csigmar1_xa_Lv("I,e") += r1_xa_Lv("I,e") * u0;
        }




/// ****** csigmar1_b(e) ****** ///



        if (include_u1_) {

            // csigmar1_xb_Lv += 1.000000 s1_b(e) 
            // flops: o0v1L1: 1 | mem: o0v1L1: 1, 
            csigmar1_xb_Lv("I,e") += s1_xb_Lv("I,e");
        }

        if (include_u0_) {

            // csigmar1_xb_Lv += 1.000000 r1_b(e) u0 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            csigmar1_xb_Lv("I,e") += r1_xb_Lv("I,e") * u0;
        }




/// ****** csigmar2_aaa(e,f,n) ****** ///



        if (include_u2_) {

            // csigmar2_xaaa_Lvvo += 1.000000 s2_aaa(e,f,n) 
            // flops: o1v2L1: 1 | mem: o1v2L1: 1, 
            csigmar2_xaaa_Lvvo("I,e,f,n") += s2_xaaa_Lvvo("I,e,f,n");
        }

        if (include_u0_) {

            // csigmar2_xaaa_Lvvo += 1.000000 r2_aaa(e,f,n) u0 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            csigmar2_xaaa_Lvvo("I,e,f,n") += r2_xaaa_Lvvo("I,e,f,n") * u0;
        }

        if (include_u1_) {

            // csigmar2_xaaa_Lvvo += -1.000000 P(e,f) r1_a(f) u1_aa(e,n) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = r1_xa_Lv("I,f") * u1_aa_vo("e,n");
            csigmar2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            csigmar2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");
        }




/// ****** csigmar2_abb(e,f,n) ****** ///



        if (include_u1_) {

            // csigmar2_xabb_Lvvo += 1.000000 r1_a(e) u1_bb(f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            csigmar2_xabb_Lvvo("I,e,f,n") += r1_xa_Lv("I,e") * u1_bb_vo("f,n");
        }




/// ****** csigmar2_baa(e,f,n) ****** ///



        if (include_u1_) {

            // csigmar2_xbaa_Lvvo += 1.000000 r1_b(e) u1_aa(f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            csigmar2_xbaa_Lvvo("I,e,f,n") += r1_xb_Lv("I,e") * u1_aa_vo("f,n");
        }




/// ****** csigmar2_bbb(e,f,n) ****** ///



        if (include_u2_) {

            // csigmar2_xbbb_Lvvo += 1.000000 s2_bbb(e,f,n) 
            // flops: o1v2L1: 1 | mem: o1v2L1: 1, 
            csigmar2_xbbb_Lvvo("I,e,f,n") += s2_xbbb_Lvvo("I,e,f,n");
        }

        if (include_u0_) {

            // csigmar2_xbbb_Lvvo += 1.000000 r2_bbb(e,f,n) u0 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            csigmar2_xbbb_Lvvo("I,e,f,n") += r2_xbbb_Lvvo("I,e,f,n") * u0;
        }

        if (include_u1_) {

            // csigmar2_xbbb_Lvvo += -1.000000 P(e,f) r1_b(f) u1_bb(e,n) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = r1_xb_Lv("I,f") * u1_bb_vo("e,n");
            csigmar2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            csigmar2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");
        }




/// ****** csigmal1_a(e) ****** ///



        if (include_u1_) {

            // csigmal1_xa_Lv += 1.000000 m1_a(e) 
            // flops: o0v1L1: 1 | mem: o0v1L1: 1, 
            csigmal1_xa_Lv("I,e") += m1_xa_Lv("I,e");
        }

        if (include_u0_) {

            // csigmal1_xa_Lv += 1.000000 l1_a(e) u0 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            csigmal1_xa_Lv("I,e") += l1_xa_Lv("I,e") * u0;
        }

        if (include_u1_) {

            // csigmal1_xa_Lv += -1.000000 l2_aaa(i,a,e) u1_aa(a,i) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            csigmal1_xa_Lv("I,e") -= l2_xaaa_Lovv("I,i,a,e") * u1_aa_vo("a,i");

            // csigmal1_xa_Lv += -1.000000 l2_bba(i,a,e) u1_bb(a,i) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            csigmal1_xa_Lv("I,e") -= l2_xbba_Lovv("I,i,a,e") * u1_bb_vo("a,i");
        }




/// ****** csigmal1_b(e) ****** ///



        if (include_u1_) {

            // csigmal1_xb_Lv += 1.000000 m1_b(e) 
            // flops: o0v1L1: 1 | mem: o0v1L1: 1, 
            csigmal1_xb_Lv("I,e") += m1_xb_Lv("I,e");
        }

        if (include_u0_) {

            // csigmal1_xb_Lv += 1.000000 l1_b(e) u0 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            csigmal1_xb_Lv("I,e") += l1_xb_Lv("I,e") * u0;
        }

        if (include_u1_) {

            // csigmal1_xb_Lv += -1.000000 l2_aab(i,a,e) u1_aa(a,i) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            csigmal1_xb_Lv("I,e") -= l2_xaab_Lovv("I,i,a,e") * u1_aa_vo("a,i");

            // csigmal1_xb_Lv += -1.000000 l2_bbb(i,a,e) u1_bb(a,i) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            csigmal1_xb_Lv("I,e") -= l2_xbbb_Lovv("I,i,a,e") * u1_bb_vo("a,i");
        }




/// ****** csigmal2_aaa(e,f,n) ****** ///



        if (include_u2_) {

            // csigmal2_xaaa_Lvvo += 1.000000 m2_aaa(n,e,f) 
            // flops: o1v2L1: 1 | mem: o1v2L1: 1, 
            csigmal2_xaaa_Lvvo("I,e,f,n") += m2_xaaa_Lovv("I,n,e,f");
        }

        if (include_u0_) {

            // csigmal2_xaaa_Lvvo += 1.000000 l2_aaa(n,e,f) u0 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            csigmal2_xaaa_Lvvo("I,e,f,n") += l2_xaaa_Lovv("I,n,e,f") * u0;
        }




/// ****** csigmal2_bbb(e,f,n) ****** ///



        if (include_u2_) {

            // csigmal2_xbbb_Lvvo += 1.000000 m2_bbb(n,e,f) 
            // flops: o1v2L1: 1 | mem: o1v2L1: 1, 
            csigmal2_xbbb_Lvvo("I,e,f,n") += m2_xbbb_Lovv("I,n,e,f");
        }

        if (include_u0_) {

            // csigmal2_xbbb_Lvvo += 1.000000 l2_bbb(n,e,f) u0 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            csigmal2_xbbb_Lvvo("I,e,f,n") += l2_xbbb_Lovv("I,n,e,f") * u0;
        }




/// ****** csigmas1_a(e) ****** ///



        if (include_u1_) {

            // csigmas1_xa_Lv += 1.000000 r1_a(e) 
            // flops: o0v1L1: 1 | mem: o0v1L1: 1, 
            csigmas1_xa_Lv("I,e") += r1_xa_Lv("I,e");
        }

        if (include_u0_ && include_u1_) {

            // csigmas1_xa_Lv += 1.000000 u0 s1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            csigmas1_xa_Lv("I,e") += u0 * s1_xa_Lv("I,e");
        }




/// ****** csigmas1_b(e) ****** ///



        if (include_u1_) {

            // csigmas1_xb_Lv += 1.000000 r1_b(e) 
            // flops: o0v1L1: 1 | mem: o0v1L1: 1, 
            csigmas1_xb_Lv("I,e") += r1_xb_Lv("I,e");
        }

        if (include_u0_ && include_u1_) {

            // csigmas1_xb_Lv += 1.000000 u0 s1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            csigmas1_xb_Lv("I,e") += u0 * s1_xb_Lv("I,e");
        }




/// ****** csigmas2_aaa(e,f,n) ****** ///



        if (include_u2_) {

            // csigmas2_xaaa_Lvvo += 1.000000 r2_aaa(e,f,n) 
            // flops: o1v2L1: 1 | mem: o1v2L1: 1, 
            csigmas2_xaaa_Lvvo("I,e,f,n") += r2_xaaa_Lvvo("I,e,f,n");
        }

        if (include_u0_ && include_u2_) {

            // csigmas2_xaaa_Lvvo += 1.000000 u0 s2_aaa(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            csigmas2_xaaa_Lvvo("I,e,f,n") += u0 * s2_xaaa_Lvvo("I,e,f,n");
        }

        if (include_u1_ && include_u2_) {

            // csigmas2_xaaa_Lvvo += -1.000000 P(e,f) u1_aa(e,n) s1_a(f) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xaaa_Lvvo("I,e,f,n") = u1_aa_vo("e,n") * s1_xa_Lv("I,f");
            csigmas2_xaaa_Lvvo("I,e,f,n") -= tempPerm_xaaa_Lvvo("I,e,f,n");
            csigmas2_xaaa_Lvvo("I,e,f,n") += tempPerm_xaaa_Lvvo("I,f,e,n");
        }




/// ****** csigmas2_abb(e,f,n) ****** ///



        if (include_u1_ && include_u2_) {

            // csigmas2_xabb_Lvvo += 1.000000 u1_bb(f,n) s1_a(e) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            csigmas2_xabb_Lvvo("I,e,f,n") += u1_bb_vo("f,n") * s1_xa_Lv("I,e");
        }




/// ****** csigmas2_baa(e,f,n) ****** ///



        if (include_u1_ && include_u2_) {

            // csigmas2_xbaa_Lvvo += 1.000000 u1_aa(f,n) s1_b(e) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 2, 
            csigmas2_xbaa_Lvvo("I,e,f,n") += u1_aa_vo("f,n") * s1_xb_Lv("I,e");
        }




/// ****** csigmas2_bbb(e,f,n) ****** ///



        if (include_u2_) {

            // csigmas2_xbbb_Lvvo += 1.000000 r2_bbb(e,f,n) 
            // flops: o1v2L1: 1 | mem: o1v2L1: 1, 
            csigmas2_xbbb_Lvvo("I,e,f,n") += r2_xbbb_Lvvo("I,e,f,n");
        }

        if (include_u0_ && include_u2_) {

            // csigmas2_xbbb_Lvvo += 1.000000 u0 s2_bbb(e,f,n) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            csigmas2_xbbb_Lvvo("I,e,f,n") += u0 * s2_xbbb_Lvvo("I,e,f,n");
        }

        if (include_u1_ && include_u2_) {

            // csigmas2_xbbb_Lvvo += -1.000000 P(e,f) u1_bb(e,n) s1_b(f) 
            // flops: o1v2L1: 3 | mem: o1v2L1: 2, 
            tempPerm_xbbb_Lvvo("I,e,f,n") = u1_bb_vo("e,n") * s1_xb_Lv("I,f");
            csigmas2_xbbb_Lvvo("I,e,f,n") -= tempPerm_xbbb_Lvvo("I,e,f,n");
            csigmas2_xbbb_Lvvo("I,e,f,n") += tempPerm_xbbb_Lvvo("I,f,e,n");
        }




/// ****** csigmam1_a(e) ****** ///



        if (include_u1_) {

            // csigmam1_xa_Lv += 1.000000 l1_a(e) 
            // flops: o0v1L1: 1 | mem: o0v1L1: 1, 
            csigmam1_xa_Lv("I,e") += l1_xa_Lv("I,e");
        }

        if (include_u0_ && include_u1_) {

            // csigmam1_xa_Lv += 1.000000 u0 m1_a(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            csigmam1_xa_Lv("I,e") += u0 * m1_xa_Lv("I,e");
        }

        if (include_u1_ && include_u2_) {

            // csigmam1_xa_Lv += -1.000000 u1_aa(a,i) m2_aaa(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            csigmam1_xa_Lv("I,e") -= u1_aa_vo("a,i") * m2_xaaa_Lovv("I,i,a,e");

            // csigmam1_xa_Lv += -1.000000 u1_bb(a,i) m2_bba(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            csigmam1_xa_Lv("I,e") -= u1_bb_vo("a,i") * m2_xbba_Lovv("I,i,a,e");
        }




/// ****** csigmam1_b(e) ****** ///



        if (include_u1_) {

            // csigmam1_xb_Lv += 1.000000 l1_b(e) 
            // flops: o0v1L1: 1 | mem: o0v1L1: 1, 
            csigmam1_xb_Lv("I,e") += l1_xb_Lv("I,e");
        }

        if (include_u0_ && include_u1_) {

            // csigmam1_xb_Lv += 1.000000 u0 m1_b(e) 
            // flops: o0v1L1: 2 | mem: o0v1L1: 1, o0v1: 1, 
            csigmam1_xb_Lv("I,e") += u0 * m1_xb_Lv("I,e");
        }

        if (include_u1_ && include_u2_) {

            // csigmam1_xb_Lv += -1.000000 u1_aa(a,i) m2_aab(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            csigmam1_xb_Lv("I,e") -= u1_aa_vo("a,i") * m2_xaab_Lovv("I,i,a,e");

            // csigmam1_xb_Lv += -1.000000 u1_bb(a,i) m2_bbb(i,a,e) 
            // flops: o1v2L1: 1, o0v1L1: 1 | mem: o0v1L1: 2, 
            csigmam1_xb_Lv("I,e") -= u1_bb_vo("a,i") * m2_xbbb_Lovv("I,i,a,e");
        }




/// ****** csigmam2_aaa(e,f,n) ****** ///



        if (include_u2_) {

            // csigmam2_xaaa_Lvvo += 1.000000 l2_aaa(n,e,f) 
            // flops: o1v2L1: 1 | mem: o1v2L1: 1, 
            csigmam2_xaaa_Lvvo("I,e,f,n") += l2_xaaa_Lovv("I,n,e,f");
        }

        if (include_u0_ && include_u2_) {

            // csigmam2_xaaa_Lvvo += 1.000000 u0 m2_aaa(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            csigmam2_xaaa_Lvvo("I,e,f,n") += u0 * m2_xaaa_Lovv("I,n,e,f");
        }
        
        /// ****** csigmam2_bbb(e,f,n) ****** ///

        if (include_u2_) {

            // csigmam2_xbbb_Lvvo += 1.000000 l2_bbb(n,e,f) 
            // flops: o1v2L1: 1 | mem: o1v2L1: 1, 
            csigmam2_xbbb_Lvvo("I,e,f,n") += l2_xbbb_Lovv("I,n,e,f");
        }

        if (include_u0_ && include_u2_) {

            // csigmam2_xbbb_Lvvo += 1.000000 u0 m2_bbb(n,e,f) 
            // flops: o1v2L1: 2 | mem: o1v2L1: 1, o1v2: 1, 
            csigmam2_xbbb_Lvvo("I,e,f,n") += u0 * m2_xbbb_Lovv("I,n,e,f");
        }
    }

    double coherent_scalar;
    double e_dip_z_ = cc_wfn_->e_dip_z_;
    double nuc_dip_z_ = cc_wfn_->nuc_dip_z_;
    if ( options_.get_bool("QED_USE_RELAXED_ORBITALS")) {
        coherent_scalar = coupling_factor_z * e_dip_z_;
    } else {
        coherent_scalar = -coupling_factor_z * nuc_dip_z_;
    }

    sigmar1_xa_Lv("I,e") += coherent_scalar * csigmar1_xa_Lv("I,e");
    sigmar1_xb_Lv("I,e") += coherent_scalar * csigmar1_xb_Lv("I,e");

    sigmar2_xaaa_Lvvo("I,e,f,n") += coherent_scalar * csigmar2_xaaa_Lvvo("I,e,f,n");
    sigmar2_xabb_Lvvo("I,e,f,n") += coherent_scalar * csigmar2_xabb_Lvvo("I,e,f,n");
    sigmar2_xbaa_Lvvo("I,e,f,n") += coherent_scalar * csigmar2_xbaa_Lvvo("I,e,f,n");
    sigmar2_xbbb_Lvvo("I,e,f,n") += coherent_scalar * csigmar2_xbbb_Lvvo("I,e,f,n");

    if (include_u1_) {
        sigmas1_xa_Lv("I,e") += coherent_scalar * csigmas1_xa_Lv("I,e");
        sigmas1_xb_Lv("I,e") += coherent_scalar * csigmas1_xb_Lv("I,e");
    }

    if (include_u2_) {
        sigmas2_xaaa_Lvvo("I,e,f,n") += coherent_scalar * csigmas2_xaaa_Lvvo("I,e,f,n");
        sigmas2_xabb_Lvvo("I,e,f,n") += coherent_scalar * csigmas2_xabb_Lvvo("I,e,f,n");
        sigmas2_xbaa_Lvvo("I,e,f,n") += coherent_scalar * csigmas2_xbaa_Lvvo("I,e,f,n");
        sigmas2_xbbb_Lvvo("I,e,f,n") += coherent_scalar * csigmas2_xbbb_Lvvo("I,e,f,n");
    }

    sigmal1_xa_Lv("I,e") += coherent_scalar * csigmal1_xa_Lv("I,e");
    sigmal1_xb_Lv("I,e") += coherent_scalar * csigmal1_xb_Lv("I,e");

    sigmal2_xaaa_Lvvo("I,e,f,n") += coherent_scalar * csigmal2_xaaa_Lvvo("I,e,f,n");
    sigmal2_xabb_Lvvo("I,e,f,n") += coherent_scalar * csigmal2_xabb_Lvvo("I,e,f,n");
    sigmal2_xbaa_Lvvo("I,e,f,n") += coherent_scalar * csigmal2_xbaa_Lvvo("I,e,f,n");
    sigmal2_xbbb_Lvvo("I,e,f,n") += coherent_scalar * csigmal2_xbbb_Lvvo("I,e,f,n");

    if (include_u1_) {
        sigmam1_xa_Lv("I,e") += coherent_scalar * csigmam1_xa_Lv("I,e");
        sigmam1_xb_Lv("I,e") += coherent_scalar * csigmam1_xb_Lv("I,e");
    }

    if (include_u2_) {
        sigmam2_xaaa_Lvvo("I,e,f,n") += coherent_scalar * csigmam2_xaaa_Lvvo("I,e,f,n");
        sigmam2_xabb_Lvvo("I,e,f,n") += coherent_scalar * csigmam2_xabb_Lvvo("I,e,f,n");
        sigmam2_xbaa_Lvvo("I,e,f,n") += coherent_scalar * csigmam2_xbaa_Lvvo("I,e,f,n");
        sigmam2_xbbb_Lvvo("I,e,f,n") += coherent_scalar * csigmam2_xbbb_Lvvo("I,e,f,n");
    }

    double average_electric_dipole_self_energy_ = cc_wfn_->average_electric_dipole_self_energy_;
    double enuc_ = cc_wfn_->enuc_;

    sigmar1_xa_Lv("I,e") += r1_xa_Lv("I,e") * (average_electric_dipole_self_energy_ + enuc_);
    sigmar1_xb_Lv("I,e") += r1_xb_Lv("I,e") * (average_electric_dipole_self_energy_ + enuc_);

    sigmal1_xa_Lv("I,e") += l1_xa_Lv("I,e") * (average_electric_dipole_self_energy_ + enuc_);
    sigmal1_xb_Lv("I,e") += l1_xb_Lv("I,e") * (average_electric_dipole_self_energy_ + enuc_);

    sigmar2_xaaa_Lvvo("I,e,f,n") += r2_xaaa_Lvvo("I,e,f,n") * (average_electric_dipole_self_energy_ + enuc_);
    sigmar2_xabb_Lvvo("I,e,f,n") += r2_xabb_Lvvo("I,e,f,n") * (average_electric_dipole_self_energy_ + enuc_);
    sigmar2_xbaa_Lvvo("I,e,f,n") += r2_xbaa_Lvvo("I,e,f,n") * (average_electric_dipole_self_energy_ + enuc_);
    sigmar2_xbbb_Lvvo("I,e,f,n") += r2_xbbb_Lvvo("I,e,f,n") * (average_electric_dipole_self_energy_ + enuc_);

    sigmal2_xaaa_Lvvo("I,e,f,n") += l2_xaaa_Lovv("I,n,f,e") * (average_electric_dipole_self_energy_ + enuc_);
    sigmal2_xabb_Lvvo("I,e,f,n") += l2_xbba_Lovv("I,n,f,e") * (average_electric_dipole_self_energy_ + enuc_);
    sigmal2_xbaa_Lvvo("I,e,f,n") += l2_xaab_Lovv("I,n,f,e") * (average_electric_dipole_self_energy_ + enuc_);
    sigmal2_xbbb_Lvvo("I,e,f,n") += l2_xbbb_Lovv("I,n,f,e") * (average_electric_dipole_self_energy_ + enuc_);

    if ( include_u1_ ) {
        sigmas1_xa_Lv("I,e") += s1_xa_Lv("I,e") * (average_electric_dipole_self_energy_ + enuc_);
        sigmas1_xb_Lv("I,e") += s1_xb_Lv("I,e") * (average_electric_dipole_self_energy_ + enuc_);

        sigmam1_xa_Lv("I,e") += m1_xa_Lv("I,e") * (average_electric_dipole_self_energy_ + enuc_);
        sigmam1_xb_Lv("I,e") += m1_xb_Lv("I,e") * (average_electric_dipole_self_energy_ + enuc_);
    }

    if ( include_u2_ ) {
        sigmas2_xaaa_Lvvo("I,e,f,n") += s2_xaaa_Lvvo("I,e,f,n") * (average_electric_dipole_self_energy_ + enuc_);
        sigmas2_xabb_Lvvo("I,e,f,n") += s2_xabb_Lvvo("I,e,f,n") * (average_electric_dipole_self_energy_ + enuc_);
        sigmas2_xbaa_Lvvo("I,e,f,n") += s2_xbaa_Lvvo("I,e,f,n") * (average_electric_dipole_self_energy_ + enuc_);
        sigmas2_xbbb_Lvvo("I,e,f,n") += s2_xbbb_Lvvo("I,e,f,n") * (average_electric_dipole_self_energy_ + enuc_);

        sigmam2_xaaa_Lvvo("I,e,f,n") += m2_xaaa_Lovv("I,n,f,e") * (average_electric_dipole_self_energy_ + enuc_);
        sigmam2_xabb_Lvvo("I,e,f,n") += m2_xbba_Lovv("I,n,f,e") * (average_electric_dipole_self_energy_ + enuc_);
        sigmam2_xbaa_Lvvo("I,e,f,n") += m2_xaab_Lovv("I,n,f,e") * (average_electric_dipole_self_energy_ + enuc_);
        sigmam2_xbbb_Lvvo("I,e,f,n") += m2_xbbb_Lovv("I,n,f,e") * (average_electric_dipole_self_energy_ + enuc_);
    }

    world_.gop.fence();
}