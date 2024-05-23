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

#include "../../../include/derived/eom_ee_ccsd.h"

double* hilbert::EOM_EE_CCSD::build_ss_diagonal() {

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
    
    // initialize diagonal blocks
    TArrayD    H_ss_aaaa_ovov = HelperD::makeTensor(world_, {oa_, va_,oa_, va_}, true);
    TArrayD    H_ss_bbbb_ovov = HelperD::makeTensor(world_, {ob_, vb_,ob_, vb_}, true);

    {
        double scalar_0;
        double scalar_1;
        double scalar_2;
        double scalar_3;
        double scalar_4;
        double scalar_5;
        double scalar_6;
        double scalar_7;
        scalar_0 = dot(V_blks_["abab_oooo"]("o0,o1,o2,o3"), Id_blks_["abab_oooo"]("o0,o1,o2,o3"));
        scalar_1 = dot(V_blks_["bbbb_oooo"]("o0,o1,o2,o3"), Id_blks_["bbbb_oooo"]("o0,o1,o2,o3"));
        scalar_2 = dot(F_blks_["aa_oo"]("o0,o1"), Id_blks_["aa_oo"]("o0,o1"));
        scalar_3 = dot(V_blks_["abab_oovv"]("k,j,b,c"), t2_abab_vvoo("b,c,k,j"));
        scalar_4 = dot(V_blks_["aaaa_oovv"]("k,j,b,c"), t2_aaaa_vvoo("b,c,k,j"));
        scalar_5 = dot(F_blks_["bb_oo"]("o0,o1"), Id_blks_["bb_oo"]("o0,o1"));
        scalar_6 = dot(V_blks_["aaaa_oooo"]("o0,o1,o2,o3"), Id_blks_["aaaa_oooo"]("o0,o1,o2,o3"));
        scalar_7 = dot(V_blks_["bbbb_oovv"]("k,j,b,c"), t2_bbbb_vvoo("b,c,k,j"));
        auto tempArray = vector<TA::TArrayD>(2);

        /// ****** p†q ****** ///

        {

            // tempArray[0] += 1.000000 Id_blks_["bb_vv"]("e,a") Id_blks_["bb_oo"]("m,i") 
            // flops: o2v2: 2 | mem: o2v2: 2, 
            tempArray[0]("e,a,m,i") = Id_blks_["bb_vv"]("e,a") * Id_blks_["bb_oo"]("m,i");

            // H_ss_bbbb_ovov += -0.500000 d_bb(e,a) d_bb(m,i) <k,j||k,j>_bbbb 
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1, 
            H_ss_bbbb_ovov("m,e,i,a") -= 0.500000 * scalar_1 * tempArray[0]("e,a,m,i");

            // H_ss_bbbb_ovov += -0.500000 d_bb(e,a) d_bb(m,i) <k,j||k,j>_abab 
            // +                 -0.500000 d_bb(e,a) d_bb(m,i) <j,k||j,k>_abab 
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1, 
            H_ss_bbbb_ovov("m,e,i,a") -= scalar_0 * tempArray[0]("e,a,m,i");

            // H_ss_bbbb_ovov += 1.000000 d_bb(e,a) d_bb(m,i) f_aa(j,j) 
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1, 
            H_ss_bbbb_ovov("m,e,i,a") += scalar_2 * tempArray[0]("e,a,m,i");

            // H_ss_bbbb_ovov += 0.250000 d_bb(e,a) d_bb(m,i) <k,j||b,c>_bbbb t2_bbbb(b,c,k,j) 
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1, 
            H_ss_bbbb_ovov("m,e,i,a") += 0.250000 * scalar_7 * tempArray[0]("e,a,m,i");

            // H_ss_bbbb_ovov += -0.500000 d_bb(e,a) d_bb(m,i) <k,j||k,j>_aaaa 
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1, 
            H_ss_bbbb_ovov("m,e,i,a") -= 0.500000 * scalar_6 * tempArray[0]("e,a,m,i");

            // H_ss_bbbb_ovov += 0.250000 d_bb(e,a) d_bb(m,i) <k,j||b,c>_aaaa t2_aaaa(b,c,k,j) 
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1, 
            H_ss_bbbb_ovov("m,e,i,a") += 0.250000 * scalar_4 * tempArray[0]("e,a,m,i");

            // H_ss_bbbb_ovov += 1.000000 d_bb(e,a) d_bb(m,i) f_bb(j,j) 
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1, 
            H_ss_bbbb_ovov("m,e,i,a") += scalar_5 * tempArray[0]("e,a,m,i");

            // H_ss_bbbb_ovov += 0.250000 d_bb(e,a) d_bb(m,i) <k,j||b,c>_abab t2_abab(b,c,k,j) 
            // +                 0.250000 d_bb(e,a) d_bb(m,i) <j,k||b,c>_abab t2_abab(b,c,j,k) 
            // +                 0.250000 d_bb(e,a) d_bb(m,i) <k,j||c,b>_abab t2_abab(c,b,k,j) 
            // +                 0.250000 d_bb(e,a) d_bb(m,i) <j,k||c,b>_abab t2_abab(c,b,j,k) 
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1, 
            H_ss_bbbb_ovov("m,e,i,a") += scalar_3 * tempArray[0]("e,a,m,i");

            // tempArray[1] += 1.000000 Id_blks_["aa_oo"]("m,i") Id_blks_["aa_vv"]("e,a") 
            // flops: o2v2: 2 | mem: o2v2: 2, 
            tempArray[1]("e,a,m,i") = Id_blks_["aa_oo"]("m,i") * Id_blks_["aa_vv"]("e,a");

            // H_ss_aaaa_ovov += 1.000000 d_aa(e,a) d_aa(m,i) f_aa(j,j) 
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1, 
            H_ss_aaaa_ovov("m,e,i,a") += scalar_2 * tempArray[1]("e,a,m,i");

            // H_ss_aaaa_ovov += 0.250000 d_aa(e,a) d_aa(m,i) <k,j||b,c>_abab t2_abab(b,c,k,j) 
            // +                 0.250000 d_aa(e,a) d_aa(m,i) <j,k||b,c>_abab t2_abab(b,c,j,k) 
            // +                 0.250000 d_aa(e,a) d_aa(m,i) <k,j||c,b>_abab t2_abab(c,b,k,j) 
            // +                 0.250000 d_aa(e,a) d_aa(m,i) <j,k||c,b>_abab t2_abab(c,b,j,k) 
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1, 
            H_ss_aaaa_ovov("m,e,i,a") += scalar_3 * tempArray[1]("e,a,m,i");

            // H_ss_aaaa_ovov += 0.250000 d_aa(e,a) d_aa(m,i) <k,j||b,c>_aaaa t2_aaaa(b,c,k,j) 
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1, 
            H_ss_aaaa_ovov("m,e,i,a") += 0.250000 * scalar_4 * tempArray[1]("e,a,m,i");

            // H_ss_aaaa_ovov += -0.500000 d_aa(e,a) d_aa(m,i) <k,j||k,j>_aaaa 
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1, 
            H_ss_aaaa_ovov("m,e,i,a") -= 0.500000 * scalar_6 * tempArray[1]("e,a,m,i");

            // H_ss_aaaa_ovov += 0.250000 d_aa(e,a) d_aa(m,i) <k,j||b,c>_bbbb t2_bbbb(b,c,k,j) 
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1, 
            H_ss_aaaa_ovov("m,e,i,a") += 0.250000 * scalar_7 * tempArray[1]("e,a,m,i");

            // H_ss_aaaa_ovov += 1.000000 d_aa(e,a) d_aa(m,i) f_bb(j,j) 
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1, 
            H_ss_aaaa_ovov("m,e,i,a") += scalar_5 * tempArray[1]("e,a,m,i");

            // H_ss_aaaa_ovov += -0.500000 d_aa(e,a) d_aa(m,i) <k,j||k,j>_abab 
            // +                 -0.500000 d_aa(e,a) d_aa(m,i) <j,k||j,k>_abab 
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1, 
            H_ss_aaaa_ovov("m,e,i,a") -= scalar_0 * tempArray[1]("e,a,m,i");

            // H_ss_aaaa_ovov += -0.500000 d_aa(e,a) d_aa(m,i) <k,j||k,j>_bbbb 
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1, 
            H_ss_aaaa_ovov("m,e,i,a") -= 0.500000 * scalar_1 * tempArray[1]("e,a,m,i");

            // H_ss_aaaa_ovov += 1.000000 d_aa(m,i) f_aa(e,a) 
            // flops: o2v2: 2 | mem: o2v2: 2, 
            H_ss_aaaa_ovov("m,e,i,a") += Id_blks_["aa_oo"]("m,i") * F_blks_["aa_vv"]("e,a");

            // H_ss_aaaa_ovov += -0.500000 d_aa(e,a) <i,j||b,c>_abab t2_abab(b,c,m,j) 
            // +                 -0.500000 d_aa(e,a) <i,j||c,b>_abab t2_abab(c,b,m,j) 
            // flops: o3v2: 1, o2v2: 2 | mem: o2v2: 2, o2v0: 1, 
            H_ss_aaaa_ovov("m,e,i,a") -= V_blks_["abab_oovv"]("i,j,b,c") * t2_abab_vvoo("b,c,m,j") * Id_blks_["aa_vv"]("e,a");

            // H_ss_aaaa_ovov += -0.500000 d_aa(e,a) <i,j||b,c>_aaaa t2_aaaa(b,c,m,j) 
            // flops: o3v2: 1, o2v2: 2 | mem: o2v2: 2, o2v0: 1, 
            H_ss_aaaa_ovov("m,e,i,a") -= 0.500000 * V_blks_["aaaa_oovv"]("i,j,b,c") * t2_aaaa_vvoo("b,c,m,j") * Id_blks_["aa_vv"]("e,a");

            // H_ss_aaaa_ovov += -0.500000 d_aa(m,i) <k,j||a,b>_abab t2_abab(e,b,k,j) 
            // +                 -0.500000 d_aa(m,i) <j,k||a,b>_abab t2_abab(e,b,j,k) 
            // flops: o2v3: 1, o2v2: 2 | mem: o2v2: 2, o0v2: 1, 
            H_ss_aaaa_ovov("m,e,i,a") -= V_blks_["abab_oovv"]("k,j,a,b") * t2_abab_vvoo("e,b,k,j") * Id_blks_["aa_oo"]("m,i");

            // H_ss_aaaa_ovov += 1.000000 <i,j||b,a>_aaaa t2_aaaa(b,e,m,j) 
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2, 
            H_ss_aaaa_ovov("m,e,i,a") += V_blks_["aaaa_oovv"]("i,j,b,a") * t2_aaaa_vvoo("b,e,m,j");

            // H_ss_aaaa_ovov += -1.000000 d_aa(e,a) f_aa(i,m) 
            // flops: o2v2: 2 | mem: o2v2: 2, 
            H_ss_aaaa_ovov("m,e,i,a") -= Id_blks_["aa_vv"]("e,a") * F_blks_["aa_oo"]("i,m");

            // H_ss_aaaa_ovov += 1.000000 <i,j||a,b>_abab t2_abab(e,b,m,j) 
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2, 
            H_ss_aaaa_ovov("m,e,i,a") += V_blks_["abab_oovv"]("i,j,a,b") * t2_abab_vvoo("e,b,m,j");

            // H_ss_aaaa_ovov += -0.500000 d_aa(m,i) <k,j||b,a>_aaaa t2_aaaa(b,e,k,j) 
            // flops: o2v3: 1, o2v2: 2 | mem: o2v2: 2, o0v2: 1, 
            H_ss_aaaa_ovov("m,e,i,a") -= 0.500000 * V_blks_["aaaa_oovv"]("k,j,b,a") * t2_aaaa_vvoo("b,e,k,j") * Id_blks_["aa_oo"]("m,i");

            // H_ss_aaaa_ovov += 1.000000 <i,e||a,m>_aaaa 
            // flops: o2v2: 1 | mem: o2v2: 1, 
            H_ss_aaaa_ovov("m,e,i,a") -= V_blks_["aaaa_vovo"]("e,i,a,m");

            // H_ss_bbbb_ovov += -0.500000 d_bb(e,a) <i,j||b,c>_bbbb t2_bbbb(b,c,m,j) 
            // flops: o3v2: 1, o2v2: 2 | mem: o2v2: 2, o2v0: 1, 
            H_ss_bbbb_ovov("m,e,i,a") -= 0.500000 * V_blks_["bbbb_oovv"]("i,j,b,c") * t2_bbbb_vvoo("b,c,m,j") * Id_blks_["bb_vv"]("e,a");

            // H_ss_bbbb_ovov += 1.000000 <j,i||b,a>_abab t2_abab(b,e,j,m) 
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2, 
            H_ss_bbbb_ovov("m,e,i,a") += V_blks_["abab_oovv"]("j,i,b,a") * t2_abab_vvoo("b,e,j,m");

            // H_ss_bbbb_ovov += -1.000000 d_bb(e,a) f_bb(i,m) 
            // flops: o2v2: 2 | mem: o2v2: 2, 
            H_ss_bbbb_ovov("m,e,i,a") -= Id_blks_["bb_vv"]("e,a") * F_blks_["bb_oo"]("i,m");

            // H_ss_bbbb_ovov += -0.500000 d_bb(m,i) <k,j||b,a>_bbbb t2_bbbb(b,e,k,j) 
            // flops: o2v3: 1, o2v2: 2 | mem: o2v2: 2, o0v2: 1, 
            H_ss_bbbb_ovov("m,e,i,a") -= 0.500000 * V_blks_["bbbb_oovv"]("k,j,b,a") * t2_bbbb_vvoo("b,e,k,j") * Id_blks_["bb_oo"]("m,i");

            // H_ss_bbbb_ovov += -0.500000 d_bb(m,i) <k,j||b,a>_abab t2_abab(b,e,k,j) 
            // +                 -0.500000 d_bb(m,i) <j,k||b,a>_abab t2_abab(b,e,j,k) 
            // flops: o2v3: 1, o2v2: 2 | mem: o2v2: 2, o0v2: 1, 
            H_ss_bbbb_ovov("m,e,i,a") -= V_blks_["abab_oovv"]("k,j,b,a") * t2_abab_vvoo("b,e,k,j") * Id_blks_["bb_oo"]("m,i");

            // H_ss_bbbb_ovov += 1.000000 <i,j||b,a>_bbbb t2_bbbb(b,e,m,j) 
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2, 
            H_ss_bbbb_ovov("m,e,i,a") += V_blks_["bbbb_oovv"]("i,j,b,a") * t2_bbbb_vvoo("b,e,m,j");

            // H_ss_bbbb_ovov += 1.000000 <i,e||a,m>_bbbb 
            // flops: o2v2: 1 | mem: o2v2: 1, 
            H_ss_bbbb_ovov("m,e,i,a") -= V_blks_["bbbb_vovo"]("e,i,a,m");

            // H_ss_bbbb_ovov += 1.000000 d_bb(m,i) f_bb(e,a) 
            // flops: o2v2: 2 | mem: o2v2: 2, 
            H_ss_bbbb_ovov("m,e,i,a") += Id_blks_["bb_oo"]("m,i") * F_blks_["bb_vv"]("e,a");

            // H_ss_bbbb_ovov += -0.500000 d_bb(e,a) <j,i||b,c>_abab t2_abab(b,c,j,m) 
            // +                 -0.500000 d_bb(e,a) <j,i||c,b>_abab t2_abab(c,b,j,m) 
            // flops: o3v2: 1, o2v2: 2 | mem: o2v2: 2, o2v0: 1, 
            H_ss_bbbb_ovov("m,e,i,a") -= V_blks_["abab_oovv"]("j,i,b,c") * t2_abab_vvoo("b,c,j,m") * Id_blks_["bb_vv"]("e,a");
        }
        tempArray[0].~TArrayD();
        tempArray[1].~TArrayD();
        
    }
        // pluck out diagonals
    auto * ss_diag = (double*) calloc(singleDim_, sizeof(double));

    size_t id = 0;
    // extract singles
    size_t oa = oa_, ob = ob_, va = va_, vb = vb_;
    HelperD::forall(H_ss_aaaa_ovov, [id,oa,va,ob,vb,ss_diag](auto &tile, auto &x){
        size_t m = x[0], e = x[1], i = x[2], a = x[3];
        if (m==i && e==a){
            size_t index = e*oa + i + id;
            ss_diag[index] = tile[x];
        }
    }); id += oa_*va_;

    HelperD::forall(H_ss_bbbb_ovov, [id,oa,va,ob,vb,ss_diag](auto &tile, auto &x){
        size_t m = x[0], e = x[1], i = x[2], a = x[3];
        if (m==i && e==a){
            size_t index = e*ob + i + id;
            ss_diag[index] = tile[x];
        }
    }); id += ob_*vb_;

    world_.gop.fence(); // TODO: broadcast to all nodes
    return ss_diag;

}

void hilbert::EOM_EE_CCSD::build_common_ops() {

    if (!cc_wfn_->has_t1_integrals_) cc_wfn_->transform_integrals(true);

    TArrayMap &eri = cc_wfn_->V_blks_;
    TArrayMap &f = cc_wfn_->F_blks_;

    TArrayMap t2 = {
            {"aaaa_vvoo", cc_wfn_->amplitudes_["t2_aaaa"]},
            {"abab_vvoo", cc_wfn_->amplitudes_["t2_abab"]},
            {"bbbb_vvoo", cc_wfn_->amplitudes_["t2_bbbb"]}
    };

    {
        reuse_tmps_["1_bbbb_vvoo"]("a,e,i,m")  = eri["abab_oovv"]("j,m,b,e") * t2["abab_vvoo"]("b,a,j,i");
        reuse_tmps_["2_aaaa_vvoo"]("a,e,i,m")  = eri["abab_oovv"]("m,j,e,b") * t2["abab_vvoo"]("a,b,i,j");
        reuse_tmps_["3_bbbb_vvoo"]("e,b,m,j")  = eri["bbbb_oovv"]("j,i,a,b") * t2["bbbb_vvoo"]("a,e,m,i");
        reuse_tmps_["4_aaaa_vvoo"]("e,b,m,j")  = eri["aaaa_oovv"]("j,i,a,b") * t2["aaaa_vvoo"]("a,e,m,i");
        reuse_tmps_["5_abab_vvoo"]("f,a,m,i")  = eri["abab_oovv"]("m,j,f,b") * t2["bbbb_vvoo"]("b,a,i,j");
        reuse_tmps_["6_abab_vvoo"]("e,a,m,i")  = eri["aaaa_oovv"]("m,j,b,e") * t2["abab_vvoo"]("b,a,j,i");
        reuse_tmps_["7_abab_vvoo"]("e,b,n,j")  = eri["abab_oovv"]("i,j,a,b") * t2["aaaa_vvoo"]("a,e,n,i");
        reuse_tmps_["8_abab_vvoo"]("e,b,n,j")  = eri["bbbb_oovv"]("j,i,a,b") * t2["abab_vvoo"]("e,a,n,i");
        reuse_tmps_["9_aabb_vvoo"]("e,b,n,j")  = eri["abab_oovv"]("i,j,b,a") * t2["abab_vvoo"]("e,a,i,n");
        reuse_tmps_["10_bbaa_vvoo"]("f,b,m,j")  = eri["abab_oovv"]("j,i,a,b") * t2["abab_vvoo"]("a,f,m,i");
        reuse_tmps_["11_aabb_oooo"]("i,m,j,n")  = eri["abab_oovv"]("m,n,a,b") * t2["abab_vvoo"]("a,b,i,j");
        reuse_tmps_["12_bbbb_oooo"]("n,m,i,j")  = 0.25 * eri["bbbb_oovv"]("j,i,a,b") * t2["bbbb_vvoo"]("a,b,m,n");
        reuse_tmps_["13_aaaa_oooo"]("n,m,i,j")  = 0.25 * eri["aaaa_oovv"]("j,i,a,b") * t2["aaaa_vvoo"]("a,b,m,n");
        reuse_tmps_["14_aa_vv"]("a,e")  = eri["abab_oovv"]("j,i,e,b") * t2["abab_vvoo"]("a,b,j,i");
        reuse_tmps_["15_bb_vv"]("a,e")  = eri["abab_oovv"]("j,i,b,e") * t2["abab_vvoo"]("b,a,j,i");
        reuse_tmps_["16_bb_vv"]("e,b")  = 0.50 * eri["bbbb_oovv"]("j,i,a,b") * t2["bbbb_vvoo"]("a,e,j,i");
        reuse_tmps_["17_aa_vv"]("e,b")  = 0.50 * eri["aaaa_oovv"]("j,i,a,b") * t2["aaaa_vvoo"]("a,e,j,i");
        reuse_tmps_["18_bb_oo"]("m,j")  = eri["abab_oovv"]("i,j,a,b") * t2["abab_vvoo"]("a,b,i,m");
        reuse_tmps_["19_bb_oo"]("m,j")  = 0.50 * eri["bbbb_oovv"]("j,i,a,b") * t2["bbbb_vvoo"]("a,b,m,i");
        reuse_tmps_["20_aa_oo"]("m,j")  = eri["abab_oovv"]("j,i,a,b") * t2["abab_vvoo"]("a,b,m,i");
        reuse_tmps_["21_aa_oo"]("m,j")  = 0.50 * eri["aaaa_oovv"]("j,i,a,b") * t2["aaaa_vvoo"]("a,b,m,i");
        reuse_tmps_["22_bbbb_vvvo"]("a,e,b,i")  = eri["baab_vovv"]("b,j,c,e") * t2["abab_vvoo"]("c,a,j,i");
        reuse_tmps_["23_abba_vvvo"]("b,a,e,i")  = eri["abab_vovv"]("b,j,c,e") * t2["abab_vvoo"]("c,a,i,j");
        reuse_tmps_["24_bbbb_vvvo"]("a,e,b,i")  = eri["bbbb_vovv"]("b,j,c,e") * t2["bbbb_vvoo"]("c,a,i,j");
        reuse_tmps_["25_aaaa_vvvo"]("a,e,b,i")  = eri["aaaa_vovv"]("b,j,c,e") * t2["aaaa_vvoo"]("c,a,i,j");
        reuse_tmps_["26_aaaa_vvvo"]("f,b,e,m")  = eri["abab_vovv"]("e,i,b,a") * t2["abab_vvoo"]("f,a,m,i");
        reuse_tmps_["27_aabb_vvvo"]("a,e,b,i")  = eri["baab_vovv"]("b,j,e,c") * t2["abab_vvoo"]("a,c,j,i");
        reuse_tmps_["28_aabb_vvvo"]("e,b,a,i")  = eri["aaaa_vovv"]("b,j,c,e") * t2["abab_vvoo"]("c,a,j,i");
        reuse_tmps_["29_abba_vvvo"]("e,b,f,m")  = eri["bbbb_vovv"]("f,i,a,b") * t2["abab_vvoo"]("e,a,m,i");
        reuse_tmps_["30_bbbb_vvvo"]("a,b,e,i")  = 0.50 * eri["bbbb_oovo"]("k,j,e,i") * t2["bbbb_vvoo"]("b,a,k,j");
        reuse_tmps_["31_aaaa_vvvo"]("a,b,e,i")  = 0.50 * eri["aaaa_oovo"]("k,j,e,i") * t2["aaaa_vvoo"]("b,a,k,j");
        reuse_tmps_["32_bbbb_vooo"]("e,n,m,i")  = 0.50 * eri["bbbb_vovv"]("e,i,a,b") * t2["bbbb_vvoo"]("a,b,m,n");
        reuse_tmps_["33_aaaa_vooo"]("e,n,m,i")  = 0.50 * eri["aaaa_vovv"]("e,i,a,b") * t2["aaaa_vvoo"]("a,b,m,n");
        reuse_tmps_["34_bb_vo"]("f,n")  = 0.50 * eri["bbbb_vovv"]("f,i,a,b") * t2["bbbb_vvoo"]("a,b,n,i");
        reuse_tmps_["35_aa_vo"]("e,m")  = 0.50 * eri["aaaa_vovv"]("e,i,a,b") * t2["aaaa_vvoo"]("a,b,m,i");
        reuse_tmps_["36_abba_vvvo"]("a,e,b,i")  = eri["baab_vovv"]("b,j,c,e") * t2["aaaa_vvoo"]("c,a,i,j");
        reuse_tmps_["37_aabb_vvvo"]("e,b,a,i")  = eri["abab_vovv"]("b,j,e,c") * t2["bbbb_vvoo"]("c,a,i,j");
        reuse_tmps_["38_abba_vvvo"]("b,a,e,i")  = eri["abba_oovo"]("k,j,e,i") * t2["abab_vvoo"]("b,a,k,j");
        reuse_tmps_["39_aabb_vvvo"]("b,e,a,i")  = eri["abab_oovo"]("k,j,e,i") * t2["abab_vvoo"]("b,a,k,j");
        reuse_tmps_["40_aabb_vooo"]("a,i,j,m")  = eri["abab_vovv"]("a,m,b,c") * t2["abab_vvoo"]("b,c,i,j");
        reuse_tmps_["41_baab_vooo"]("a,i,m,j")  = eri["baab_vovv"]("a,m,b,c") * t2["abab_vvoo"]("b,c,i,j");
        reuse_tmps_["42_bb_vo"]("f,n")  = eri["baab_vovv"]("f,i,a,b") * t2["abab_vvoo"]("a,b,i,n");
        reuse_tmps_["43_aa_vo"]("e,m")  = eri["abab_vovv"]("e,i,a,b") * t2["abab_vvoo"]("a,b,m,i");
        reuse_tmps_["44_aabb_vooo"]("a,j,i,m")  = eri["abba_oovo"]("k,m,b,j") * t2["abab_vvoo"]("a,b,k,i");
        reuse_tmps_["45_baab_vooo"]("a,i,m,j")  = eri["abab_oovo"]("m,k,b,j") * t2["abab_vvoo"]("b,a,i,k");
        reuse_tmps_["46_bbbb_vooo"]("f,m,n,j")  = eri["bbbb_oovo"]("j,i,a,n") * t2["bbbb_vvoo"]("a,f,m,i");
        reuse_tmps_["47_aaaa_vooo"]("f,m,n,j")  = eri["aaaa_oovo"]("j,i,a,n") * t2["aaaa_vvoo"]("a,f,m,i");
        reuse_tmps_["48_aaaa_vooo"]("a,i,j,m")  = eri["abba_oovo"]("m,k,b,j") * t2["abab_vvoo"]("a,b,i,k");
        reuse_tmps_["49_bbbb_vooo"]("f,n,m,j")  = eri["abab_oovo"]("i,j,a,m") * t2["abab_vvoo"]("a,f,i,n");
        reuse_tmps_["50_bbbb_vooo"]("e,n,m,i")  = f["bb_ov"]("i,a") * t2["bbbb_vvoo"]("a,e,m,n");
        reuse_tmps_["51_aaaa_vooo"]("e,n,m,i")  = f["aa_ov"]("i,a") * t2["aaaa_vvoo"]("a,e,m,n");
        reuse_tmps_["52_aabb_vooo"]("a,i,j,m")  = eri["bbbb_oovo"]("m,k,b,j") * t2["abab_vvoo"]("a,b,i,k");
        reuse_tmps_["53_baab_vooo"]("f,m,j,n")  = eri["aaaa_oovo"]("j,i,a,m") * t2["abab_vvoo"]("a,f,i,n");
        reuse_tmps_["54_baab_vooo"]("a,i,m,j")  = f["aa_ov"]("m,b") * t2["abab_vvoo"]("b,a,i,j");
        reuse_tmps_["55_aabb_vooo"]("e,m,n,i")  = f["bb_ov"]("i,a") * t2["abab_vvoo"]("e,a,m,n");
        reuse_tmps_["56_bb_vo"]("f,n")  = 0.50 * eri["bbbb_oovo"]("j,i,a,n") * t2["bbbb_vvoo"]("a,f,j,i");
        reuse_tmps_["57_aa_vo"]("f,m")  = 0.50 * eri["aaaa_oovo"]("j,i,a,m") * t2["aaaa_vvoo"]("a,f,j,i");
        reuse_tmps_["58_bb_vo"]("f,n")  = f["bb_ov"]("i,a") * t2["bbbb_vvoo"]("a,f,n,i");
        reuse_tmps_["59_aa_vo"]("e,m")  = f["aa_ov"]("i,a") * t2["aaaa_vvoo"]("a,e,m,i");
        reuse_tmps_["60_aabb_vooo"]("a,i,j,m")  = eri["abab_oovo"]("k,m,b,j") * t2["aaaa_vvoo"]("b,a,i,k");
        reuse_tmps_["61_baab_vooo"]("f,m,j,n")  = eri["abba_oovo"]("j,i,a,m") * t2["bbbb_vvoo"]("a,f,n,i");
        reuse_tmps_["62_bb_vo"]("f,n")  = eri["abab_oovo"]("j,i,a,n") * t2["abab_vvoo"]("a,f,j,i");
        reuse_tmps_["63_aa_vo"]("f,n")  = eri["abba_oovo"]("j,i,a,n") * t2["abab_vvoo"]("f,a,j,i");
        reuse_tmps_["64_bb_vo"]("a,i")  = f["aa_ov"]("j,b") * t2["abab_vvoo"]("b,a,j,i");
        reuse_tmps_["65_aa_vo"]("f,m")  = f["bb_ov"]("i,a") * t2["abab_vvoo"]("f,a,m,i");
    }

    world_.gop.fence();
}

void hilbert::EOM_EE_CCSD::build_Hc_cH(size_t L) {

    // ensure that the integrals are transformed with t1 amplitudes
    if (!cc_wfn_->has_t1_integrals_) cc_wfn_->transform_integrals(true);

    /// build sigma vectors for this trial
    sigvec_blks_.clear();

    // reduce sigma into spins
    sigvec_blks_["sigmar0"]   = HelperD::makeTensor(world_, {L}, true);
    sigvec_blks_["sigmar1_aa"]   = HelperD::makeTensor(world_, {L, va_, oa_}, true);
    sigvec_blks_["sigmar1_bb"]   = HelperD::makeTensor(world_, {L, vb_, ob_}, true);
    sigvec_blks_["sigmar2_aaaa"] = HelperD::makeTensor(world_, {L, va_, va_, oa_, oa_}, true);
    sigvec_blks_["sigmar2_abab"] = HelperD::makeTensor(world_, {L, va_, vb_, oa_, ob_}, true);
    sigvec_blks_["sigmar2_bbbb"] = HelperD::makeTensor(world_, {L, vb_, vb_, ob_, ob_}, true);

    sigvec_blks_["sigmal0"]      = HelperD::makeTensor(world_, {L}, true);
    sigvec_blks_["sigmal1_aa"]   = HelperD::makeTensor(world_, {L, va_, oa_}, true);
    sigvec_blks_["sigmal1_bb"]   = HelperD::makeTensor(world_, {L, vb_, ob_}, true);
    sigvec_blks_["sigmal2_aaaa"] = HelperD::makeTensor(world_, {L, va_, va_, oa_, oa_}, true);
    sigvec_blks_["sigmal2_abab"] = HelperD::makeTensor(world_, {L, va_, vb_, oa_, ob_}, true);
    sigvec_blks_["sigmal2_bbbb"] = HelperD::makeTensor(world_, {L, vb_, vb_, ob_, ob_}, true);

    // get reference to electronic integrals
    TArrayMap &V_blks_ = cc_wfn_->V_blks_;
    TArrayMap &F_blks_ = cc_wfn_->F_blks_;
    TArrayMap &Id_blks_ = cc_wfn_->Id_blks_;
    TArrayMap &eri = cc_wfn_->V_blks_;
    TArrayMap &f = cc_wfn_->F_blks_;
    TArrayMap &Id = cc_wfn_->Id_blks_;

    TArrayMap t2 = {
            {"aaaa_vvoo", cc_wfn_->amplitudes_["t2_aaaa"]},
            {"abab_vvoo", cc_wfn_->amplitudes_["t2_abab"]},
            {"bbbb_vvoo", cc_wfn_->amplitudes_["t2_bbbb"]}
    };

    /// reference right operators
    TArrayD r0 = evec_blks_["r0"];
    TArrayMap r1 = {
            {"aa_Lvo", evec_blks_["r1_aa"]},
            {"bb_Lvo", evec_blks_["r1_bb"]}
    };

    TArrayMap r2 = {
            {"aaaa_Lvvoo", evec_blks_["r2_aaaa"]},
            {"abab_Lvvoo", evec_blks_["r2_abab"]},
            {"bbbb_Lvvoo", evec_blks_["r2_bbbb"]}
    };

    /// reference left operators
    TArrayD l0 = evec_blks_["l0"];
    TArrayMap l1 = {
            {"aa_Lov", evec_blks_["l1_aa"]},
            {"bb_Lov", evec_blks_["l1_bb"]}
    };
    TArrayMap l2 = {
            {"aaaa_Loovv", evec_blks_["l2_aaaa"]},
            {"abab_Loovv", evec_blks_["l2_abab"]},
            {"bbbb_Loovv", evec_blks_["l2_bbbb"]}
    };

    /// reference sigma vectors
    auto &sigmar0 = sigvec_blks_["sigmar0"];
    auto &sigmal0 = sigvec_blks_["sigmal0"];

    auto &sigmar1_aa = sigvec_blks_["sigmar1_aa"];
    auto &sigmar1_bb = sigvec_blks_["sigmar1_bb"];
    auto &sigmal1_aa = sigvec_blks_["sigmal1_aa"];
    auto &sigmal1_bb = sigvec_blks_["sigmal1_bb"];

    auto &sigmar2_aaaa = sigvec_blks_["sigmar2_aaaa"];
    auto &sigmar2_abab = sigvec_blks_["sigmar2_abab"];
    auto &sigmar2_bbbb = sigvec_blks_["sigmar2_bbbb"];
    auto &sigmal2_aaaa = sigvec_blks_["sigmal2_aaaa"];
    auto &sigmal2_abab = sigvec_blks_["sigmal2_abab"];
    auto &sigmal2_bbbb = sigvec_blks_["sigmal2_bbbb"];

    {
        TArrayMap tmps_, perm_tmps;
        std::map<std::string, double> scalars_;

        // sigmar1_aa = -0.50 <j,i||a,b>_aaaa t2_aaaa(a,b,m,i) r1_aa(e,j)  // flops: o1v1L1 = o2v1L1 | mem: o1v1L1 = o1v1L1
        sigmar1_aa("X,e,m")  = -1.00 * reuse_tmps_["21_aa_oo"]("m,j") * r1["aa_Lvo"]("X,e,j");

        // sigmar1_aa += -1.00 f_aa(i,m) r1_aa(e,i)  // flops: o1v1L1 += o2v1L1 | mem: o1v1L1 += o1v1L1
        sigmar1_aa("X,e,m") -= f["aa_oo"]("i,m") * r1["aa_Lvo"]("X,e,i");

        // sigmar1_aa += -0.50 <j,i||a,b>_abab t2_abab(a,b,m,i) r1_aa(e,j) @ -0.50 <j,i||b,a>_abab t2_abab(b,a,m,i) r1_aa(e,j)  // flops: o1v1L1 += o2v1L1 | mem: o1v1L1 += o1v1L1
        sigmar1_aa("X,e,m") -= reuse_tmps_["20_aa_oo"]("m,j") * r1["aa_Lvo"]("X,e,j");

        // sigmar1_bb = -1.00 f_bb(i,m) r1_bb(e,i)  // flops: o1v1L1 = o2v1L1 | mem: o1v1L1 = o1v1L1
        sigmar1_bb("X,e,m")  = -1.00 * f["bb_oo"]("i,m") * r1["bb_Lvo"]("X,e,i");

        // sigmar1_bb += -0.50 <i,j||a,b>_abab t2_abab(a,b,i,m) r1_bb(e,j) @ -0.50 <i,j||b,a>_abab t2_abab(b,a,i,m) r1_bb(e,j)  // flops: o1v1L1 += o2v1L1 | mem: o1v1L1 += o1v1L1
        sigmar1_bb("X,e,m") -= reuse_tmps_["18_bb_oo"]("m,j") * r1["bb_Lvo"]("X,e,j");

        // sigmar1_bb += -0.50 <j,i||a,b>_bbbb t2_bbbb(a,b,m,i) r1_bb(e,j)  // flops: o1v1L1 += o2v1L1 | mem: o1v1L1 += o1v1L1
        sigmar1_bb("X,e,m") -= reuse_tmps_["19_bb_oo"]("m,j") * r1["bb_Lvo"]("X,e,j");

        // sigmal1_aa = -0.50 <m,j||a,b>_abab l1_aa(i,e) t2_abab(a,b,i,j) @ -0.50 <m,j||b,a>_abab l1_aa(i,e) t2_abab(b,a,i,j)  // flops: o1v1L1 = o2v1L1 | mem: o1v1L1 = o1v1L1
        sigmal1_aa("X,e,m")  = -1.00 * reuse_tmps_["20_aa_oo"]("i,m") * l1["aa_Lov"]("X,i,e");

        // sigmal1_aa += -0.50 <m,j||a,b>_aaaa l1_aa(i,e) t2_aaaa(a,b,i,j)  // flops: o1v1L1 += o2v1L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= reuse_tmps_["21_aa_oo"]("i,m") * l1["aa_Lov"]("X,i,e");

        // sigmal1_aa += -1.00 f_aa(m,i) l1_aa(i,e)  // flops: o1v1L1 += o2v1L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= f["aa_oo"]("m,i") * l1["aa_Lov"]("X,i,e");

        // sigmal1_bb = -0.50 <m,j||a,b>_bbbb l1_bb(i,e) t2_bbbb(a,b,i,j)  // flops: o1v1L1 = o2v1L1 | mem: o1v1L1 = o1v1L1
        sigmal1_bb("X,e,m")  = -1.00 * reuse_tmps_["19_bb_oo"]("i,m") * l1["bb_Lov"]("X,i,e");

        // sigmal1_bb += -0.50 <j,m||a,b>_abab l1_bb(i,e) t2_abab(a,b,j,i) @ -0.50 <j,m||b,a>_abab l1_bb(i,e) t2_abab(b,a,j,i)  // flops: o1v1L1 += o2v1L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= reuse_tmps_["18_bb_oo"]("i,m") * l1["bb_Lov"]("X,i,e");

        // sigmal1_bb += -1.00 f_bb(m,i) l1_bb(i,e)  // flops: o1v1L1 += o2v1L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= f["bb_oo"]("m,i") * l1["bb_Lov"]("X,i,e");

        // sigmar1_aa += -0.50 <j,i||a,b>_aaaa t2_aaaa(a,e,j,i) r1_aa(b,m)  // flops: o1v1L1 += o1v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_aa("X,e,m") -= reuse_tmps_["17_aa_vv"]("e,b") * r1["aa_Lvo"]("X,b,m");

        // sigmar1_aa += +1.00 f_aa(e,a) r1_aa(a,m)  // flops: o1v1L1 += o1v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_aa("X,e,m") += f["aa_vv"]("e,a") * r1["aa_Lvo"]("X,a,m");

        // sigmar1_aa += -0.50 <j,i||b,a>_abab t2_abab(e,a,j,i) r1_aa(b,m) @ -0.50 <i,j||b,a>_abab t2_abab(e,a,i,j) r1_aa(b,m)  // flops: o1v1L1 += o1v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_aa("X,e,m") -= reuse_tmps_["14_aa_vv"]("e,b") * r1["aa_Lvo"]("X,b,m");

        // sigmar1_bb += +1.00 f_bb(e,a) r1_bb(a,m)  // flops: o1v1L1 += o1v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_bb("X,e,m") += f["bb_vv"]("e,a") * r1["bb_Lvo"]("X,a,m");

        // sigmar1_bb += -0.50 <j,i||a,b>_abab t2_abab(a,e,j,i) r1_bb(b,m) @ -0.50 <i,j||a,b>_abab t2_abab(a,e,i,j) r1_bb(b,m)  // flops: o1v1L1 += o1v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_bb("X,e,m") -= reuse_tmps_["15_bb_vv"]("e,b") * r1["bb_Lvo"]("X,b,m");

        // sigmar1_bb += -0.50 <j,i||a,b>_bbbb t2_bbbb(a,e,j,i) r1_bb(b,m)  // flops: o1v1L1 += o1v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_bb("X,e,m") -= reuse_tmps_["16_bb_vv"]("e,b") * r1["bb_Lvo"]("X,b,m");

        // sigmal1_aa += -0.50 <j,i||e,b>_abab l1_aa(m,a) t2_abab(a,b,j,i) @ -0.50 <i,j||e,b>_abab l1_aa(m,a) t2_abab(a,b,i,j)  // flops: o1v1L1 += o1v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= reuse_tmps_["14_aa_vv"]("a,e") * l1["aa_Lov"]("X,m,a");

        // sigmal1_aa += -0.50 <j,i||b,e>_aaaa l1_aa(m,a) t2_aaaa(b,a,j,i)  // flops: o1v1L1 += o1v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= reuse_tmps_["17_aa_vv"]("a,e") * l1["aa_Lov"]("X,m,a");

        // sigmal1_aa += +1.00 f_aa(a,e) l1_aa(m,a)  // flops: o1v1L1 += o1v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += f["aa_vv"]("a,e") * l1["aa_Lov"]("X,m,a");

        // sigmal1_bb += -0.50 <j,i||b,e>_bbbb l1_bb(m,a) t2_bbbb(b,a,j,i)  // flops: o1v1L1 += o1v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= reuse_tmps_["16_bb_vv"]("a,e") * l1["bb_Lov"]("X,m,a");

        // sigmal1_bb += -0.50 <j,i||b,e>_abab l1_bb(m,a) t2_abab(b,a,j,i) @ -0.50 <i,j||b,e>_abab l1_bb(m,a) t2_abab(b,a,i,j)  // flops: o1v1L1 += o1v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= reuse_tmps_["15_bb_vv"]("a,e") * l1["bb_Lov"]("X,m,a");

        // sigmal1_bb += +1.00 f_bb(a,e) l1_bb(m,a)  // flops: o1v1L1 += o1v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += f["bb_vv"]("a,e") * l1["bb_Lov"]("X,m,a");

        // sigmar1_aa += +1.00 <j,i||a,b>_aaaa t2_aaaa(a,e,m,i) r1_aa(b,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_aa("X,e,m") += reuse_tmps_["4_aaaa_vvoo"]("e,b,m,j") * r1["aa_Lvo"]("X,b,j");

        // sigmar1_aa += -1.00 <i,j||a,b>_abab t2_aaaa(a,e,m,i) r1_bb(b,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_aa("X,e,m") -= reuse_tmps_["7_abab_vvoo"]("e,b,m,j") * r1["bb_Lvo"]("X,b,j");

        // sigmar1_aa += -1.00 f_aa(i,a) r2_aaaa(a,e,m,i)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_aa("X,e,m") -= f["aa_ov"]("i,a") * r2["aaaa_Lvvoo"]("X,a,e,m,i");

        // sigmar1_aa += +1.00 <e,i||m,a>_abab r1_bb(a,i)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_aa("X,e,m") -= eri["abba_vovo"]("e,i,a,m") * r1["bb_Lvo"]("X,a,i");

        // sigmar1_aa += +1.00 <i,e||a,m>_aaaa r1_aa(a,i)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_aa("X,e,m") -= eri["aaaa_vovo"]("e,i,a,m") * r1["aa_Lvo"]("X,a,i");

        // sigmar1_aa += +1.00 f_bb(i,a) r2_abab(e,a,m,i)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_aa("X,e,m") += f["bb_ov"]("i,a") * r2["abab_Lvvoo"]("X,e,a,m,i");

        // sigmar1_aa += -1.00 <j,i||a,b>_bbbb t2_abab(e,a,m,i) r1_bb(b,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_aa("X,e,m") -= reuse_tmps_["8_abab_vvoo"]("e,b,m,j") * r1["bb_Lvo"]("X,b,j");

        // sigmar1_aa += +1.00 <j,i||b,a>_abab t2_abab(e,a,m,i) r1_aa(b,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_aa("X,e,m") += reuse_tmps_["2_aaaa_vvoo"]("e,b,m,j") * r1["aa_Lvo"]("X,b,j");

        // sigmar1_bb += -1.00 f_bb(i,a) r2_bbbb(a,e,m,i)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_bb("X,e,m") -= f["bb_ov"]("i,a") * r2["bbbb_Lvvoo"]("X,a,e,m,i");

        // sigmar1_bb += -1.00 <j,i||b,a>_abab t2_bbbb(a,e,m,i) r1_aa(b,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_bb("X,e,m") -= reuse_tmps_["5_abab_vvoo"]("b,e,j,m") * r1["aa_Lvo"]("X,b,j");

        // sigmar1_bb += -1.00 <j,i||a,b>_aaaa t2_abab(a,e,i,m) r1_aa(b,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_bb("X,e,m") -= reuse_tmps_["6_abab_vvoo"]("b,e,j,m") * r1["aa_Lvo"]("X,b,j");

        // sigmar1_bb += +1.00 <i,e||a,m>_abab r1_aa(a,i)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_bb("X,e,m") -= eri["baab_vovo"]("e,i,a,m") * r1["aa_Lvo"]("X,a,i");

        // sigmar1_bb += +1.00 <i,j||a,b>_abab t2_abab(a,e,i,m) r1_bb(b,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_bb("X,e,m") += reuse_tmps_["1_bbbb_vvoo"]("e,b,m,j") * r1["bb_Lvo"]("X,b,j");

        // sigmar1_bb += +1.00 <i,e||a,m>_bbbb r1_bb(a,i)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_bb("X,e,m") -= eri["bbbb_vovo"]("e,i,a,m") * r1["bb_Lvo"]("X,a,i");

        // sigmar1_bb += +1.00 <j,i||a,b>_bbbb t2_bbbb(a,e,m,i) r1_bb(b,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_bb("X,e,m") += reuse_tmps_["3_bbbb_vvoo"]("e,b,m,j") * r1["bb_Lvo"]("X,b,j");

        // sigmar1_bb += +1.00 f_aa(i,a) r2_abab(a,e,i,m)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_bb("X,e,m") += f["aa_ov"]("i,a") * r2["abab_Lvvoo"]("X,a,e,i,m");
        sigmar2_aaaa("X,e,f,m,n")  = reuse_tmps_["59_aa_vo"]("f,m") * r1["aa_Lvo"]("X,e,n");
        sigmar2_aaaa("X,e,f,m,n") -= reuse_tmps_["35_aa_vo"]("f,m") * r1["aa_Lvo"]("X,e,n");

        // sigmar2_aaaa += +0.50 P(m,n) P(e,f) <j,i||a,n>_aaaa t2_aaaa(a,e,j,i) r1_aa(f,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += reuse_tmps_["57_aa_vo"]("e,n") * r1["aa_Lvo"]("X,f,m");
        sigmar2_aaaa("X,e,f,m,n") += reuse_tmps_["65_aa_vo"]("f,n") * r1["aa_Lvo"]("X,e,m");

        // sigmar2_aaaa += -1.00 P(m,n) P(e,f) f_aa(e,n) r1_aa(f,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") -= f["aa_vo"]("e,n") * r1["aa_Lvo"]("X,f,m");
        sigmar2_aaaa("X,e,f,m,n") += reuse_tmps_["65_aa_vo"]("e,m") * r1["aa_Lvo"]("X,f,n");

        // sigmar2_aaaa += +0.50 P(m,n) P(e,f) <i,e||a,b>_aaaa t2_aaaa(a,b,n,i) r1_aa(f,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") -= reuse_tmps_["35_aa_vo"]("e,n") * r1["aa_Lvo"]("X,f,m");
        sigmar2_aaaa("X,e,f,m,n") += f["aa_vo"]("e,m") * r1["aa_Lvo"]("X,f,n");

        // sigmar2_aaaa += -1.00 P(m,n) P(e,f) f_bb(i,a) t2_abab(e,a,n,i) r1_aa(f,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") -= reuse_tmps_["65_aa_vo"]("e,n") * r1["aa_Lvo"]("X,f,m");
        sigmar2_aaaa("X,e,f,m,n") += reuse_tmps_["63_aa_vo"]("f,n") * r1["aa_Lvo"]("X,e,m");
        sigmar2_aaaa("X,e,f,m,n") -= reuse_tmps_["65_aa_vo"]("f,m") * r1["aa_Lvo"]("X,e,n");

        // sigmar2_aaaa += +0.50 P(m,n) P(e,f) <j,i||n,a>_abab t2_abab(e,a,j,i) r1_aa(f,m) @ +0.50 P(m,n) P(e,f) <i,j||n,a>_abab t2_abab(e,a,i,j) r1_aa(f,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") -= reuse_tmps_["63_aa_vo"]("e,n") * r1["aa_Lvo"]("X,f,m");

        // sigmar2_aaaa += +1.00 P(m,n) P(e,f) f_aa(i,a) t2_aaaa(a,e,n,i) r1_aa(f,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += reuse_tmps_["59_aa_vo"]("e,n") * r1["aa_Lvo"]("X,f,m");

        // sigmar2_aaaa += -0.50 P(m,n) P(e,f) <e,i||a,b>_abab t2_abab(a,b,n,i) r1_aa(f,m) @ -0.50 P(m,n) P(e,f) <e,i||b,a>_abab t2_abab(b,a,n,i) r1_aa(f,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") -= reuse_tmps_["43_aa_vo"]("e,n") * r1["aa_Lvo"]("X,f,m");
        sigmar2_aaaa("X,e,f,m,n") -= reuse_tmps_["57_aa_vo"]("f,n") * r1["aa_Lvo"]("X,e,m");
        sigmar2_aaaa("X,e,f,m,n") += f["aa_vo"]("f,n") * r1["aa_Lvo"]("X,e,m");
        sigmar2_aaaa("X,e,f,m,n") -= reuse_tmps_["59_aa_vo"]("e,m") * r1["aa_Lvo"]("X,f,n");
        sigmar2_aaaa("X,e,f,m,n") -= reuse_tmps_["63_aa_vo"]("f,m") * r1["aa_Lvo"]("X,e,n");
        sigmar2_aaaa("X,e,f,m,n") += reuse_tmps_["35_aa_vo"]("f,n") * r1["aa_Lvo"]("X,e,m");
        sigmar2_aaaa("X,e,f,m,n") -= reuse_tmps_["57_aa_vo"]("e,m") * r1["aa_Lvo"]("X,f,n");
        sigmar2_aaaa("X,e,f,m,n") -= reuse_tmps_["59_aa_vo"]("f,n") * r1["aa_Lvo"]("X,e,m");
        sigmar2_aaaa("X,e,f,m,n") += reuse_tmps_["43_aa_vo"]("f,n") * r1["aa_Lvo"]("X,e,m");
        sigmar2_aaaa("X,e,f,m,n") -= f["aa_vo"]("f,m") * r1["aa_Lvo"]("X,e,n");
        sigmar2_aaaa("X,e,f,m,n") += reuse_tmps_["63_aa_vo"]("e,m") * r1["aa_Lvo"]("X,f,n");
        sigmar2_aaaa("X,e,f,m,n") -= reuse_tmps_["43_aa_vo"]("f,m") * r1["aa_Lvo"]("X,e,n");
        sigmar2_aaaa("X,e,f,m,n") += reuse_tmps_["57_aa_vo"]("f,m") * r1["aa_Lvo"]("X,e,n");
        sigmar2_aaaa("X,e,f,m,n") += reuse_tmps_["35_aa_vo"]("e,m") * r1["aa_Lvo"]("X,f,n");
        sigmar2_aaaa("X,e,f,m,n") += reuse_tmps_["43_aa_vo"]("e,m") * r1["aa_Lvo"]("X,f,n");

        // sigmar2_abab = -0.50 <i,e||a,b>_aaaa t2_aaaa(a,b,m,i) r1_bb(f,n)  // flops: o2v2L1 = o2v2L1 | mem: o2v2L1 = o2v2L1
        sigmar2_abab("X,e,f,m,n")  = reuse_tmps_["35_aa_vo"]("e,m") * r1["bb_Lvo"]("X,f,n");

        // sigmar2_abab += -1.00 f_bb(i,a) t2_bbbb(a,f,n,i) r1_aa(e,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["58_bb_vo"]("f,n") * r1["aa_Lvo"]("X,e,m");

        // sigmar2_abab += -0.50 <j,i||a,m>_aaaa t2_aaaa(a,e,j,i) r1_bb(f,n)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["57_aa_vo"]("e,m") * r1["bb_Lvo"]("X,f,n");

        // sigmar2_abab += -0.50 <j,i||m,a>_abab t2_abab(e,a,j,i) r1_bb(f,n) @ -0.50 <i,j||m,a>_abab t2_abab(e,a,i,j) r1_bb(f,n)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["63_aa_vo"]("e,m") * r1["bb_Lvo"]("X,f,n");

        // sigmar2_abab += -0.50 <i,f||a,b>_bbbb t2_bbbb(a,b,n,i) r1_aa(e,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["34_bb_vo"]("f,n") * r1["aa_Lvo"]("X,e,m");

        // sigmar2_abab += +0.50 <i,f||a,b>_abab t2_abab(a,b,i,n) r1_aa(e,m) @ +0.50 <i,f||b,a>_abab t2_abab(b,a,i,n) r1_aa(e,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["42_bb_vo"]("f,n") * r1["aa_Lvo"]("X,e,m");

        // sigmar2_abab += +1.00 f_bb(i,a) t2_abab(e,a,m,i) r1_bb(f,n)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["65_aa_vo"]("e,m") * r1["bb_Lvo"]("X,f,n");

        // sigmar2_abab += +1.00 f_bb(f,n) r1_aa(e,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += f["bb_vo"]("f,n") * r1["aa_Lvo"]("X,e,m");

        // sigmar2_abab += +0.50 <e,i||a,b>_abab t2_abab(a,b,m,i) r1_bb(f,n) @ +0.50 <e,i||b,a>_abab t2_abab(b,a,m,i) r1_bb(f,n)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["43_aa_vo"]("e,m") * r1["bb_Lvo"]("X,f,n");

        // sigmar2_abab += -0.50 <j,i||a,n>_bbbb t2_bbbb(a,f,j,i) r1_aa(e,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["56_bb_vo"]("f,n") * r1["aa_Lvo"]("X,e,m");

        // sigmar2_abab += -1.00 f_aa(i,a) t2_aaaa(a,e,m,i) r1_bb(f,n)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["59_aa_vo"]("e,m") * r1["bb_Lvo"]("X,f,n");

        // sigmar2_abab += +1.00 f_aa(i,a) t2_abab(a,f,i,n) r1_aa(e,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["64_bb_vo"]("f,n") * r1["aa_Lvo"]("X,e,m");

        // sigmar2_abab += -0.50 <j,i||a,n>_abab t2_abab(a,f,j,i) r1_aa(e,m) @ -0.50 <i,j||a,n>_abab t2_abab(a,f,i,j) r1_aa(e,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["62_bb_vo"]("f,n") * r1["aa_Lvo"]("X,e,m");

        // sigmar2_abab += +1.00 f_aa(e,m) r1_bb(f,n)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += f["aa_vo"]("e,m") * r1["bb_Lvo"]("X,f,n");
        sigmar2_bbbb("X,e,f,m,n")  = f["bb_vo"]("f,n") * r1["bb_Lvo"]("X,e,m");
        sigmar2_bbbb("X,e,f,m,n") -= reuse_tmps_["56_bb_vo"]("e,m") * r1["bb_Lvo"]("X,f,n");
        sigmar2_bbbb("X,e,f,m,n") -= reuse_tmps_["34_bb_vo"]("f,m") * r1["bb_Lvo"]("X,e,n");
        sigmar2_bbbb("X,e,f,m,n") += reuse_tmps_["64_bb_vo"]("f,n") * r1["bb_Lvo"]("X,e,m");

        // sigmar2_bbbb += +0.50 P(m,n) P(e,f) <j,i||a,n>_abab t2_abab(a,e,j,i) r1_bb(f,m) @ +0.50 P(m,n) P(e,f) <i,j||a,n>_abab t2_abab(a,e,i,j) r1_bb(f,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += reuse_tmps_["62_bb_vo"]("e,n") * r1["bb_Lvo"]("X,f,m");
        sigmar2_bbbb("X,e,f,m,n") += reuse_tmps_["62_bb_vo"]("f,m") * r1["bb_Lvo"]("X,e,n");
        sigmar2_bbbb("X,e,f,m,n") -= reuse_tmps_["56_bb_vo"]("f,n") * r1["bb_Lvo"]("X,e,m");
        sigmar2_bbbb("X,e,f,m,n") += reuse_tmps_["34_bb_vo"]("f,n") * r1["bb_Lvo"]("X,e,m");
        sigmar2_bbbb("X,e,f,m,n") -= f["bb_vo"]("f,m") * r1["bb_Lvo"]("X,e,n");
        sigmar2_bbbb("X,e,f,m,n") += reuse_tmps_["64_bb_vo"]("e,m") * r1["bb_Lvo"]("X,f,n");
        sigmar2_bbbb("X,e,f,m,n") -= reuse_tmps_["58_bb_vo"]("f,n") * r1["bb_Lvo"]("X,e,m");

        // sigmar2_bbbb += -1.00 P(m,n) P(e,f) f_aa(i,a) t2_abab(a,e,i,n) r1_bb(f,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= reuse_tmps_["64_bb_vo"]("e,n") * r1["bb_Lvo"]("X,f,m");

        // sigmar2_bbbb += +1.00 P(m,n) P(e,f) f_bb(i,a) t2_bbbb(a,e,n,i) r1_bb(f,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += reuse_tmps_["58_bb_vo"]("e,n") * r1["bb_Lvo"]("X,f,m");
        sigmar2_bbbb("X,e,f,m,n") -= reuse_tmps_["42_bb_vo"]("e,m") * r1["bb_Lvo"]("X,f,n");

        // sigmar2_bbbb += +0.50 P(m,n) P(e,f) <j,i||a,n>_bbbb t2_bbbb(a,e,j,i) r1_bb(f,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += reuse_tmps_["56_bb_vo"]("e,n") * r1["bb_Lvo"]("X,f,m");
        sigmar2_bbbb("X,e,f,m,n") -= reuse_tmps_["64_bb_vo"]("f,m") * r1["bb_Lvo"]("X,e,n");
        sigmar2_bbbb("X,e,f,m,n") += f["bb_vo"]("e,m") * r1["bb_Lvo"]("X,f,n");
        sigmar2_bbbb("X,e,f,m,n") += reuse_tmps_["42_bb_vo"]("f,m") * r1["bb_Lvo"]("X,e,n");
        sigmar2_bbbb("X,e,f,m,n") -= reuse_tmps_["42_bb_vo"]("f,n") * r1["bb_Lvo"]("X,e,m");
        sigmar2_bbbb("X,e,f,m,n") += reuse_tmps_["34_bb_vo"]("e,m") * r1["bb_Lvo"]("X,f,n");
        sigmar2_bbbb("X,e,f,m,n") -= reuse_tmps_["62_bb_vo"]("e,m") * r1["bb_Lvo"]("X,f,n");
        sigmar2_bbbb("X,e,f,m,n") -= reuse_tmps_["58_bb_vo"]("e,m") * r1["bb_Lvo"]("X,f,n");
        sigmar2_bbbb("X,e,f,m,n") -= reuse_tmps_["62_bb_vo"]("f,n") * r1["bb_Lvo"]("X,e,m");
        sigmar2_bbbb("X,e,f,m,n") += reuse_tmps_["56_bb_vo"]("f,m") * r1["bb_Lvo"]("X,e,n");
        sigmar2_bbbb("X,e,f,m,n") += reuse_tmps_["58_bb_vo"]("f,m") * r1["bb_Lvo"]("X,e,n");

        // sigmar2_bbbb += -0.50 P(m,n) P(e,f) <i,e||a,b>_abab t2_abab(a,b,i,n) r1_bb(f,m) @ -0.50 P(m,n) P(e,f) <i,e||b,a>_abab t2_abab(b,a,i,n) r1_bb(f,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += reuse_tmps_["42_bb_vo"]("e,n") * r1["bb_Lvo"]("X,f,m");

        // sigmar2_bbbb += +0.50 P(m,n) P(e,f) <i,e||a,b>_bbbb t2_bbbb(a,b,n,i) r1_bb(f,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= reuse_tmps_["34_bb_vo"]("e,n") * r1["bb_Lvo"]("X,f,m");

        // sigmar2_bbbb += -1.00 P(m,n) P(e,f) f_bb(e,n) r1_bb(f,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= f["bb_vo"]("e,n") * r1["bb_Lvo"]("X,f,m");

        // sigmal1_aa += +1.00 <m,a||e,i>_abab l1_bb(i,a)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= eri["baab_vovo"]("a,m,e,i") * l1["bb_Lov"]("X,i,a");

        // sigmal1_aa += -1.00 f_aa(a,i) l2_aaaa(m,i,a,e)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= f["aa_vo"]("a,i") * l2["aaaa_Loovv"]("X,m,i,a,e");

        // sigmal1_aa += +0.50 <k,j||b,i>_aaaa l2_aaaa(m,i,a,e) t2_aaaa(b,a,k,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += reuse_tmps_["57_aa_vo"]("a,i") * l2["aaaa_Loovv"]("X,m,i,a,e");

        // sigmal1_aa += -1.00 <m,j||e,b>_abab l1_bb(i,a) t2_bbbb(b,a,i,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= l1["bb_Lov"]("X,i,a") * reuse_tmps_["5_abab_vvoo"]("e,a,m,i");

        // sigmal1_aa += -0.50 <j,a||b,c>_bbbb l2_abab(m,i,e,a) t2_bbbb(b,c,i,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += reuse_tmps_["34_bb_vo"]("a,i") * l2["abab_Loovv"]("X,m,i,e,a");

        // sigmal1_aa += -0.50 <a,j||b,c>_abab l2_aaaa(m,i,a,e) t2_abab(b,c,i,j) @ -0.50 <a,j||c,b>_abab l2_aaaa(m,i,a,e) t2_abab(c,b,i,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= reuse_tmps_["43_aa_vo"]("a,i") * l2["aaaa_Loovv"]("X,m,i,a,e");

        // sigmal1_aa += +1.00 f_aa(j,b) l2_aaaa(m,i,a,e) t2_aaaa(b,a,i,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += reuse_tmps_["59_aa_vo"]("a,i") * l2["aaaa_Loovv"]("X,m,i,a,e");

        // sigmal1_aa += +0.50 <j,a||b,c>_abab l2_abab(m,i,e,a) t2_abab(b,c,j,i) @ +0.50 <j,a||c,b>_abab l2_abab(m,i,e,a) t2_abab(c,b,j,i)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= reuse_tmps_["42_bb_vo"]("a,i") * l2["abab_Loovv"]("X,m,i,e,a");

        // sigmal1_aa += -0.50 <k,j||b,i>_bbbb l2_abab(m,i,e,a) t2_bbbb(b,a,k,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= reuse_tmps_["56_bb_vo"]("a,i") * l2["abab_Loovv"]("X,m,i,e,a");

        // sigmal1_aa += +1.00 f_aa(j,b) l2_abab(m,i,e,a) t2_abab(b,a,j,i)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += reuse_tmps_["64_bb_vo"]("a,i") * l2["abab_Loovv"]("X,m,i,e,a");

        // sigmal1_aa += -1.00 f_bb(j,b) l2_aaaa(m,i,a,e) t2_abab(a,b,i,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= reuse_tmps_["65_aa_vo"]("a,i") * l2["aaaa_Loovv"]("X,m,i,a,e");

        // sigmal1_aa += +0.50 <k,j||i,b>_abab l2_aaaa(m,i,a,e) t2_abab(a,b,k,j) @ +0.50 <j,k||i,b>_abab l2_aaaa(m,i,a,e) t2_abab(a,b,j,k)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= reuse_tmps_["63_aa_vo"]("a,i") * l2["aaaa_Loovv"]("X,m,i,a,e");

        // sigmal1_aa += +0.50 <j,a||b,c>_aaaa l2_aaaa(m,i,a,e) t2_aaaa(b,c,i,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= reuse_tmps_["35_aa_vo"]("a,i") * l2["aaaa_Loovv"]("X,m,i,a,e");

        // sigmal1_aa += -1.00 <m,j||b,e>_aaaa l1_bb(i,a) t2_abab(b,a,j,i)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= l1["bb_Lov"]("X,i,a") * reuse_tmps_["6_abab_vvoo"]("e,a,m,i");

        // sigmal1_aa += +1.00 f_bb(a,i) l2_abab(m,i,e,a)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += f["bb_vo"]("a,i") * l2["abab_Loovv"]("X,m,i,e,a");

        // sigmal1_aa += +1.00 <m,a||e,i>_aaaa l1_aa(i,a)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= eri["aaaa_vovo"]("a,m,e,i") * l1["aa_Lov"]("X,i,a");

        // sigmal1_aa += -1.00 f_bb(j,b) l2_abab(m,i,e,a) t2_bbbb(b,a,i,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= reuse_tmps_["58_bb_vo"]("a,i") * l2["abab_Loovv"]("X,m,i,e,a");

        // sigmal1_aa += +1.00 <m,j||b,e>_aaaa l1_aa(i,a) t2_aaaa(b,a,i,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += l1["aa_Lov"]("X,i,a") * reuse_tmps_["4_aaaa_vvoo"]("a,e,i,m");

        // sigmal1_aa += +1.00 <m,j||e,b>_abab l1_aa(i,a) t2_abab(a,b,i,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += l1["aa_Lov"]("X,i,a") * reuse_tmps_["2_aaaa_vvoo"]("a,e,i,m");

        // sigmal1_aa += -0.50 <k,j||b,i>_abab l2_abab(m,i,e,a) t2_abab(b,a,k,j) @ -0.50 <j,k||b,i>_abab l2_abab(m,i,e,a) t2_abab(b,a,j,k)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= reuse_tmps_["62_bb_vo"]("a,i") * l2["abab_Loovv"]("X,m,i,e,a");

        // sigmal1_bb += -0.50 <j,a||b,c>_aaaa l2_abab(i,m,a,e) t2_aaaa(b,c,i,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += reuse_tmps_["35_aa_vo"]("a,i") * l2["abab_Loovv"]("X,i,m,a,e");

        // sigmal1_bb += +0.50 <a,j||b,c>_abab l2_abab(i,m,a,e) t2_abab(b,c,i,j) @ +0.50 <a,j||c,b>_abab l2_abab(i,m,a,e) t2_abab(c,b,i,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += reuse_tmps_["43_aa_vo"]("a,i") * l2["abab_Loovv"]("X,i,m,a,e");

        // sigmal1_bb += +0.50 <k,j||b,i>_bbbb l2_bbbb(m,i,a,e) t2_bbbb(b,a,k,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += reuse_tmps_["56_bb_vo"]("a,i") * l2["bbbb_Loovv"]("X,m,i,a,e");

        // sigmal1_bb += -0.50 <k,j||b,i>_aaaa l2_abab(i,m,a,e) t2_aaaa(b,a,k,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= reuse_tmps_["57_aa_vo"]("a,i") * l2["abab_Loovv"]("X,i,m,a,e");

        // sigmal1_bb += +1.00 f_aa(a,i) l2_abab(i,m,a,e)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += f["aa_vo"]("a,i") * l2["abab_Loovv"]("X,i,m,a,e");

        // sigmal1_bb += +0.50 <k,j||b,i>_abab l2_bbbb(m,i,a,e) t2_abab(b,a,k,j) @ +0.50 <j,k||b,i>_abab l2_bbbb(m,i,a,e) t2_abab(b,a,j,k)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += reuse_tmps_["62_bb_vo"]("a,i") * l2["bbbb_Loovv"]("X,m,i,a,e");

        // sigmal1_bb += +1.00 <a,m||i,e>_abab l1_aa(i,a)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= eri["abba_vovo"]("a,m,e,i") * l1["aa_Lov"]("X,i,a");

        // sigmal1_bb += +1.00 f_bb(j,b) l2_bbbb(m,i,a,e) t2_bbbb(b,a,i,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += reuse_tmps_["58_bb_vo"]("a,i") * l2["bbbb_Loovv"]("X,m,i,a,e");

        // sigmal1_bb += +1.00 <m,j||b,e>_bbbb l1_bb(i,a) t2_bbbb(b,a,i,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += l1["bb_Lov"]("X,i,a") * reuse_tmps_["3_bbbb_vvoo"]("a,e,i,m");

        // sigmal1_bb += +1.00 f_bb(j,b) l2_abab(i,m,a,e) t2_abab(a,b,i,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += reuse_tmps_["65_aa_vo"]("a,i") * l2["abab_Loovv"]("X,i,m,a,e");

        // sigmal1_bb += -1.00 f_aa(j,b) l2_abab(i,m,a,e) t2_aaaa(b,a,i,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= reuse_tmps_["59_aa_vo"]("a,i") * l2["abab_Loovv"]("X,i,m,a,e");

        // sigmal1_bb += -1.00 f_bb(a,i) l2_bbbb(m,i,a,e)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= f["bb_vo"]("a,i") * l2["bbbb_Loovv"]("X,m,i,a,e");

        // sigmal1_bb += -1.00 <m,j||b,e>_bbbb l1_aa(i,a) t2_abab(a,b,i,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= l1["aa_Lov"]("X,i,a") * reuse_tmps_["8_abab_vvoo"]("a,e,i,m");

        // sigmal1_bb += +0.50 <j,a||b,c>_bbbb l2_bbbb(m,i,a,e) t2_bbbb(b,c,i,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= reuse_tmps_["34_bb_vo"]("a,i") * l2["bbbb_Loovv"]("X,m,i,a,e");

        // sigmal1_bb += +1.00 <j,m||b,e>_abab l1_bb(i,a) t2_abab(b,a,j,i)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += l1["bb_Lov"]("X,i,a") * reuse_tmps_["1_bbbb_vvoo"]("a,e,i,m");

        // sigmal1_bb += -1.00 f_aa(j,b) l2_bbbb(m,i,a,e) t2_abab(b,a,j,i)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= reuse_tmps_["64_bb_vo"]("a,i") * l2["bbbb_Loovv"]("X,m,i,a,e");

        // sigmal1_bb += +1.00 <m,a||e,i>_bbbb l1_bb(i,a)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= eri["bbbb_vovo"]("a,m,e,i") * l1["bb_Lov"]("X,i,a");

        // sigmal1_bb += -1.00 <j,m||b,e>_abab l1_aa(i,a) t2_aaaa(b,a,i,j)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= l1["aa_Lov"]("X,i,a") * reuse_tmps_["7_abab_vvoo"]("a,e,i,m");

        // sigmal1_bb += -0.50 <k,j||i,b>_abab l2_abab(i,m,a,e) t2_abab(a,b,k,j) @ -0.50 <j,k||i,b>_abab l2_abab(i,m,a,e) t2_abab(a,b,j,k)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += reuse_tmps_["63_aa_vo"]("a,i") * l2["abab_Loovv"]("X,i,m,a,e");

        // sigmal1_bb += -0.50 <j,a||b,c>_abab l2_bbbb(m,i,a,e) t2_abab(b,c,j,i) @ -0.50 <j,a||c,b>_abab l2_bbbb(m,i,a,e) t2_abab(c,b,j,i)  // flops: o1v1L1 += o2v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += reuse_tmps_["42_bb_vo"]("a,i") * l2["bbbb_Loovv"]("X,m,i,a,e");
        sigmal2_aaaa("X,e,f,m,n")  = -1.00 * f["aa_ov"]("m,f") * l1["aa_Lov"]("X,n,e");
        sigmal2_aaaa("X,e,f,m,n") += f["aa_ov"]("n,f") * l1["aa_Lov"]("X,m,e");

        // sigmal2_aaaa += -1.00 P(m,n) P(e,f) f_aa(n,e) l1_aa(m,f)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") -= f["aa_ov"]("n,e") * l1["aa_Lov"]("X,m,f");
        sigmal2_aaaa("X,e,f,m,n") += f["aa_ov"]("m,e") * l1["aa_Lov"]("X,n,f");

        // sigmal2_abab = +1.00 f_aa(m,e) l1_bb(n,f)  // flops: o2v2L1 = o2v2L1 | mem: o2v2L1 = o2v2L1
        sigmal2_abab("X,e,f,m,n")  = f["aa_ov"]("m,e") * l1["bb_Lov"]("X,n,f");

        // sigmal2_abab += +1.00 f_bb(n,f) l1_aa(m,e)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += f["bb_ov"]("n,f") * l1["aa_Lov"]("X,m,e");

        // sigmal2_bbbb = -1.00 P(m,n) P(e,f) f_bb(n,e) l1_bb(m,f)  // flops: o2v2L1 = o2v2L1 | mem: o2v2L1 = o2v2L1
        sigmal2_bbbb("X,e,f,m,n")  = -1.00 * f["bb_ov"]("n,e") * l1["bb_Lov"]("X,m,f");
        sigmal2_bbbb("X,e,f,m,n") += f["bb_ov"]("n,f") * l1["bb_Lov"]("X,m,e");
        sigmal2_bbbb("X,e,f,m,n") -= f["bb_ov"]("m,f") * l1["bb_Lov"]("X,n,e");
        sigmal2_bbbb("X,e,f,m,n") += f["bb_ov"]("m,e") * l1["bb_Lov"]("X,n,f");

        // sigmar1_aa += -0.50 <j,i||a,m>_aaaa r2_aaaa(a,e,j,i)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_aa("X,e,m") -= 0.50 * eri["aaaa_oovo"]("j,i,a,m") * r2["aaaa_Lvvoo"]("X,a,e,j,i");

        // sigmar1_aa += -0.50 <j,i||m,a>_abab r2_abab(e,a,j,i) @ -0.50 <i,j||m,a>_abab r2_abab(e,a,i,j)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_aa("X,e,m") += eri["abba_oovo"]("j,i,a,m") * r2["abab_Lvvoo"]("X,e,a,j,i");

        // sigmar1_bb += -0.50 <j,i||a,m>_bbbb r2_bbbb(a,e,j,i)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_bb("X,e,m") -= 0.50 * eri["bbbb_oovo"]("j,i,a,m") * r2["bbbb_Lvvoo"]("X,a,e,j,i");

        // sigmar1_bb += -0.50 <j,i||a,m>_abab r2_abab(a,e,j,i) @ -0.50 <i,j||a,m>_abab r2_abab(a,e,i,j)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmar1_bb("X,e,m") -= eri["abab_oovo"]("j,i,a,m") * r2["abab_Lvvoo"]("X,a,e,j,i");

        // sigmar2_abab += +1.00 <j,i||m,a>_abab t2_bbbb(a,f,n,i) r1_aa(e,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["61_baab_vooo"]("f,m,j,n") * r1["aa_Lvo"]("X,e,j");

        // sigmar2_abab += -1.00 f_aa(i,m) r2_abab(e,f,i,n)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= f["aa_oo"]("i,m") * r2["abab_Lvvoo"]("X,e,f,i,n");

        // sigmar2_abab += -0.50 <j,i||a,b>_aaaa t2_aaaa(a,b,m,i) r2_abab(e,f,j,n)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["21_aa_oo"]("m,j") * r2["abab_Lvvoo"]("X,e,f,j,n");

        // sigmar2_abab += -1.00 <i,f||m,n>_abab r1_aa(e,i)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += eri["baab_vooo"]("f,i,m,n") * r1["aa_Lvo"]("X,e,i");

        // sigmar2_abab += -1.00 f_bb(i,n) r2_abab(e,f,m,i)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= f["bb_oo"]("i,n") * r2["abab_Lvvoo"]("X,e,f,m,i");

        // sigmar2_abab += +1.00 <i,j||a,n>_abab t2_aaaa(a,e,m,i) r1_bb(f,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["60_aabb_vooo"]("e,m,n,j") * r1["bb_Lvo"]("X,f,j");

        // sigmar2_abab += +1.00 <i,j||m,a>_abab t2_abab(e,a,i,n) r1_bb(f,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["44_aabb_vooo"]("e,m,n,j") * r1["bb_Lvo"]("X,f,j");

        // sigmar2_abab += -0.50 <i,f||a,b>_abab t2_abab(a,b,m,n) r1_aa(e,i) @ -0.50 <i,f||b,a>_abab t2_abab(b,a,m,n) r1_aa(e,i)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["41_baab_vooo"]("f,m,i,n") * r1["aa_Lvo"]("X,e,i");

        // sigmar2_abab += -1.00 f_bb(i,a) t2_abab(e,a,m,n) r1_bb(f,i)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["55_aabb_vooo"]("e,m,n,i") * r1["bb_Lvo"]("X,f,i");

        // sigmar2_abab += -0.50 <e,i||a,b>_abab t2_abab(a,b,m,n) r1_bb(f,i) @ -0.50 <e,i||b,a>_abab t2_abab(b,a,m,n) r1_bb(f,i)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["40_aabb_vooo"]("e,m,n,i") * r1["bb_Lvo"]("X,f,i");

        // sigmar2_abab += -0.50 <i,j||a,b>_abab t2_abab(a,b,i,n) r2_abab(e,f,m,j) @ -0.50 <i,j||b,a>_abab t2_abab(b,a,i,n) r2_abab(e,f,m,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["18_bb_oo"]("n,j") * r2["abab_Lvvoo"]("X,e,f,m,j");

        // sigmar2_abab += -0.50 <j,i||a,b>_abab t2_abab(a,b,m,i) r2_abab(e,f,j,n) @ -0.50 <j,i||b,a>_abab t2_abab(b,a,m,i) r2_abab(e,f,j,n)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["20_aa_oo"]("m,j") * r2["abab_Lvvoo"]("X,e,f,j,n");

        // sigmar2_abab += -1.00 f_aa(i,a) t2_abab(a,f,m,n) r1_aa(e,i)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["54_baab_vooo"]("f,m,i,n") * r1["aa_Lvo"]("X,e,i");

        // sigmar2_abab += +1.00 <j,i||a,m>_aaaa t2_abab(a,f,i,n) r1_aa(e,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["53_baab_vooo"]("f,m,j,n") * r1["aa_Lvo"]("X,e,j");

        // sigmar2_abab += -0.50 <j,i||a,b>_bbbb t2_bbbb(a,b,n,i) r2_abab(e,f,m,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["19_bb_oo"]("n,j") * r2["abab_Lvvoo"]("X,e,f,m,j");

        // sigmar2_abab += +1.00 <j,i||a,n>_bbbb t2_abab(e,a,m,i) r1_bb(f,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["52_aabb_vooo"]("e,m,n,j") * r1["bb_Lvo"]("X,f,j");

        // sigmar2_abab += +1.00 <j,i||a,n>_abab t2_abab(a,f,m,i) r1_aa(e,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["45_baab_vooo"]("f,m,j,n") * r1["aa_Lvo"]("X,e,j");

        // sigmar2_abab += -1.00 <e,i||m,n>_abab r1_bb(f,i)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= eri["abab_vooo"]("e,i,m,n") * r1["bb_Lvo"]("X,f,i");

        // sigmal1_aa += +1.00 <m,k||j,b>_abab l2_abab(j,i,e,a) t2_bbbb(b,a,i,k)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= reuse_tmps_["61_baab_vooo"]("a,j,m,i") * l2["abab_Loovv"]("X,j,i,e,a");

        // sigmal1_aa += -1.00 <m,k||b,j>_aaaa l2_aaaa(i,j,a,e) t2_aaaa(b,a,i,k)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= reuse_tmps_["47_aaaa_vooo"]("a,i,j,m") * l2["aaaa_Loovv"]("X,i,j,a,e");

        // sigmal1_aa += +0.50 f_aa(m,b) l2_aaaa(i,j,a,e) t2_aaaa(b,a,i,j)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += 0.50 * reuse_tmps_["51_aaaa_vooo"]("a,j,i,m") * l2["aaaa_Loovv"]("X,i,j,a,e");

        // sigmal1_aa += -0.25 <m,a||b,c>_abab l2_abab(i,j,e,a) t2_abab(b,c,i,j) @ -0.25 <m,a||c,b>_abab l2_abab(i,j,e,a) t2_abab(c,b,i,j) @ -0.25 <m,a||b,c>_abab l2_abab(j,i,e,a) t2_abab(b,c,j,i) @ -0.25 <m,a||c,b>_abab l2_abab(j,i,e,a) t2_abab(c,b,j,i)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += reuse_tmps_["41_baab_vooo"]("a,i,m,j") * l2["abab_Loovv"]("X,i,j,e,a");

        // sigmal1_aa += -0.50 f_aa(m,b) l2_abab(i,j,e,a) t2_abab(b,a,i,j) @ -0.50 f_aa(m,b) l2_abab(j,i,e,a) t2_abab(b,a,j,i)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= reuse_tmps_["54_baab_vooo"]("a,i,m,j") * l2["abab_Loovv"]("X,i,j,e,a");

        // sigmal1_aa += -1.00 <m,k||j,b>_abab l2_aaaa(i,j,a,e) t2_abab(a,b,i,k)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += reuse_tmps_["48_aaaa_vooo"]("a,i,j,m") * l2["aaaa_Loovv"]("X,i,j,a,e");

        // sigmal1_aa += +0.25 <m,a||b,c>_aaaa l2_aaaa(i,j,a,e) t2_aaaa(b,c,i,j)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= 0.50 * reuse_tmps_["33_aaaa_vooo"]("a,j,i,m") * l2["aaaa_Loovv"]("X,i,j,a,e");

        // sigmal1_aa += +0.50 <m,a||i,j>_aaaa l2_aaaa(i,j,a,e)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= 0.50 * eri["aaaa_vooo"]("a,m,i,j") * l2["aaaa_Loovv"]("X,i,j,a,e");

        // sigmal1_aa += +1.00 <m,k||b,j>_aaaa l2_abab(j,i,e,a) t2_abab(b,a,k,i)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += reuse_tmps_["53_baab_vooo"]("a,j,m,i") * l2["abab_Loovv"]("X,j,i,e,a");

        // sigmal1_aa += -0.50 <m,a||i,j>_abab l2_abab(i,j,e,a) @ -0.50 <m,a||j,i>_abab l2_abab(j,i,e,a)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += eri["baab_vooo"]("a,m,i,j") * l2["abab_Loovv"]("X,i,j,e,a");

        // sigmal1_aa += +1.00 <m,k||b,j>_abab l2_abab(i,j,e,a) t2_abab(b,a,i,k)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += reuse_tmps_["45_baab_vooo"]("a,i,m,j") * l2["abab_Loovv"]("X,i,j,e,a");

        // sigmal1_bb += -0.25 <a,m||b,c>_abab l2_abab(i,j,a,e) t2_abab(b,c,i,j) @ -0.25 <a,m||c,b>_abab l2_abab(i,j,a,e) t2_abab(c,b,i,j) @ -0.25 <a,m||b,c>_abab l2_abab(j,i,a,e) t2_abab(b,c,j,i) @ -0.25 <a,m||c,b>_abab l2_abab(j,i,a,e) t2_abab(c,b,j,i)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= reuse_tmps_["40_aabb_vooo"]("a,i,j,m") * l2["abab_Loovv"]("X,i,j,a,e");

        // sigmal1_bb += -1.00 <m,k||b,j>_bbbb l2_bbbb(i,j,a,e) t2_bbbb(b,a,i,k)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= reuse_tmps_["46_bbbb_vooo"]("a,i,j,m") * l2["bbbb_Loovv"]("X,i,j,a,e");

        // sigmal1_bb += -1.00 <k,m||b,j>_abab l2_bbbb(i,j,a,e) t2_abab(b,a,k,i)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= reuse_tmps_["49_bbbb_vooo"]("a,i,j,m") * l2["bbbb_Loovv"]("X,i,j,a,e");

        // sigmal1_bb += +0.50 f_bb(m,b) l2_bbbb(i,j,a,e) t2_bbbb(b,a,i,j)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += 0.50 * reuse_tmps_["50_bbbb_vooo"]("a,j,i,m") * l2["bbbb_Loovv"]("X,i,j,a,e");

        // sigmal1_bb += -0.50 f_bb(m,b) l2_abab(i,j,a,e) t2_abab(a,b,i,j) @ -0.50 f_bb(m,b) l2_abab(j,i,a,e) t2_abab(a,b,j,i)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= reuse_tmps_["55_aabb_vooo"]("a,i,j,m") * l2["abab_Loovv"]("X,i,j,a,e");

        // sigmal1_bb += +1.00 <k,m||b,j>_abab l2_abab(i,j,a,e) t2_aaaa(b,a,i,k)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += reuse_tmps_["60_aabb_vooo"]("a,i,j,m") * l2["abab_Loovv"]("X,i,j,a,e");

        // sigmal1_bb += +0.50 <m,a||i,j>_bbbb l2_bbbb(i,j,a,e)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= 0.50 * eri["bbbb_vooo"]("a,m,i,j") * l2["bbbb_Loovv"]("X,i,j,a,e");

        // sigmal1_bb += +0.25 <m,a||b,c>_bbbb l2_bbbb(i,j,a,e) t2_bbbb(b,c,i,j)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= 0.50 * reuse_tmps_["32_bbbb_vooo"]("a,j,i,m") * l2["bbbb_Loovv"]("X,i,j,a,e");

        // sigmal1_bb += -0.50 <a,m||i,j>_abab l2_abab(i,j,a,e) @ -0.50 <a,m||j,i>_abab l2_abab(j,i,a,e)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= eri["abab_vooo"]("a,m,i,j") * l2["abab_Loovv"]("X,i,j,a,e");

        // sigmal1_bb += +1.00 <m,k||b,j>_bbbb l2_abab(i,j,a,e) t2_abab(a,b,i,k)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += reuse_tmps_["52_aabb_vooo"]("a,i,j,m") * l2["abab_Loovv"]("X,i,j,a,e");

        // sigmal1_bb += +1.00 <k,m||j,b>_abab l2_abab(j,i,a,e) t2_abab(a,b,k,i)  // flops: o1v1L1 += o3v2L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= reuse_tmps_["44_aabb_vooo"]("a,j,i,m") * l2["abab_Loovv"]("X,j,i,a,e");

        // sigmal2_abab += -0.50 <n,j||a,b>_bbbb l2_abab(m,i,e,f) t2_bbbb(a,b,i,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= reuse_tmps_["19_bb_oo"]("i,n") * l2["abab_Loovv"]("X,m,i,e,f");

        // sigmal2_abab += -0.50 <j,n||a,b>_abab l2_abab(m,i,e,f) t2_abab(a,b,j,i) @ -0.50 <j,n||b,a>_abab l2_abab(m,i,e,f) t2_abab(b,a,j,i)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= reuse_tmps_["18_bb_oo"]("i,n") * l2["abab_Loovv"]("X,m,i,e,f");

        // sigmal2_abab += -0.50 <m,j||a,b>_abab l2_abab(i,n,e,f) t2_abab(a,b,i,j) @ -0.50 <m,j||b,a>_abab l2_abab(i,n,e,f) t2_abab(b,a,i,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= reuse_tmps_["20_aa_oo"]("i,m") * l2["abab_Loovv"]("X,i,n,e,f");

        // sigmal2_abab += -1.00 <m,n||i,f>_abab l1_aa(i,e)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += eri["abba_oovo"]("m,n,f,i") * l1["aa_Lov"]("X,i,e");

        // sigmal2_abab += -1.00 <m,n||e,i>_abab l1_bb(i,f)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= eri["abab_oovo"]("m,n,e,i") * l1["bb_Lov"]("X,i,f");

        // sigmal2_abab += -0.50 <m,j||a,b>_aaaa l2_abab(i,n,e,f) t2_aaaa(a,b,i,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= reuse_tmps_["21_aa_oo"]("i,m") * l2["abab_Loovv"]("X,i,n,e,f");

        // sigmal2_abab += -1.00 f_bb(n,i) l2_abab(m,i,e,f)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= f["bb_oo"]("n,i") * l2["abab_Loovv"]("X,m,i,e,f");

        // sigmal2_abab += -1.00 f_aa(m,i) l2_abab(i,n,e,f)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= f["aa_oo"]("m,i") * l2["abab_Loovv"]("X,i,n,e,f");

        // sigmar1_aa += -0.50 <i,e||a,b>_aaaa r2_aaaa(a,b,m,i)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmar1_aa("X,e,m") += 0.50 * eri["aaaa_vovv"]("e,i,a,b") * r2["aaaa_Lvvoo"]("X,a,b,m,i");

        // sigmar1_aa += +0.50 <e,i||a,b>_abab r2_abab(a,b,m,i) @ +0.50 <e,i||b,a>_abab r2_abab(b,a,m,i)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmar1_aa("X,e,m") += eri["abab_vovv"]("e,i,a,b") * r2["abab_Lvvoo"]("X,a,b,m,i");

        // sigmar1_bb += +0.50 <i,e||a,b>_abab r2_abab(a,b,i,m) @ +0.50 <i,e||b,a>_abab r2_abab(b,a,i,m)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmar1_bb("X,e,m") -= eri["baab_vovv"]("e,i,a,b") * r2["abab_Lvvoo"]("X,a,b,i,m");

        // sigmar1_bb += -0.50 <i,e||a,b>_bbbb r2_bbbb(a,b,m,i)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmar1_bb("X,e,m") += 0.50 * eri["bbbb_vovv"]("e,i,a,b") * r2["bbbb_Lvvoo"]("X,a,b,m,i");

        // sigmar2_abab += +1.00 <e,f||m,a>_abab r1_bb(a,n)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= eri["abba_vvvo"]("e,f,a,m") * r1["bb_Lvo"]("X,a,n");

        // sigmar2_abab += -0.50 <j,i||a,b>_aaaa t2_aaaa(a,e,j,i) r2_abab(b,f,m,n)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["17_aa_vv"]("e,b") * r2["abab_Lvvoo"]("X,b,f,m,n");

        // sigmar2_abab += -1.00 <i,f||a,b>_abab t2_aaaa(a,e,m,i) r1_bb(b,n)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["36_abba_vvvo"]("e,b,f,m") * r1["bb_Lvo"]("X,b,n");

        // sigmar2_abab += -1.00 <e,i||b,a>_abab t2_bbbb(a,f,n,i) r1_aa(b,m)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["37_aabb_vvvo"]("b,e,f,n") * r1["aa_Lvo"]("X,b,m");

        // sigmar2_abab += +0.50 <j,i||m,a>_abab t2_abab(e,f,j,i) r1_bb(a,n) @ +0.50 <i,j||m,a>_abab t2_abab(e,f,i,j) r1_bb(a,n)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["38_abba_vvvo"]("e,f,a,m") * r1["bb_Lvo"]("X,a,n");

        // sigmar2_abab += -1.00 <i,f||b,a>_abab t2_abab(e,a,i,n) r1_aa(b,m)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["27_aabb_vvvo"]("e,b,f,n") * r1["aa_Lvo"]("X,b,m");

        // sigmar2_abab += +1.00 <i,e||a,b>_aaaa t2_abab(a,f,i,n) r1_aa(b,m)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["28_aabb_vvvo"]("b,e,f,n") * r1["aa_Lvo"]("X,b,m");

        // sigmar2_abab += -1.00 <e,i||a,b>_abab t2_abab(a,f,m,i) r1_bb(b,n)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["23_abba_vvvo"]("e,f,b,m") * r1["bb_Lvo"]("X,b,n");

        // sigmar2_abab += -0.50 <j,i||a,b>_abab t2_abab(a,f,j,i) r2_abab(e,b,m,n) @ -0.50 <i,j||a,b>_abab t2_abab(a,f,i,j) r2_abab(e,b,m,n)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["15_bb_vv"]("f,b") * r2["abab_Lvvoo"]("X,e,b,m,n");

        // sigmar2_abab += -0.50 <j,i||b,a>_abab t2_abab(e,a,j,i) r2_abab(b,f,m,n) @ -0.50 <i,j||b,a>_abab t2_abab(e,a,i,j) r2_abab(b,f,m,n)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["14_aa_vv"]("e,b") * r2["abab_Lvvoo"]("X,b,f,m,n");

        // sigmar2_abab += +1.00 f_aa(e,a) r2_abab(a,f,m,n)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += f["aa_vv"]("e,a") * r2["abab_Lvvoo"]("X,a,f,m,n");

        // sigmar2_abab += -0.50 <j,i||a,b>_bbbb t2_bbbb(a,f,j,i) r2_abab(e,b,m,n)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["16_bb_vv"]("f,b") * r2["abab_Lvvoo"]("X,e,b,m,n");

        // sigmar2_abab += +1.00 <e,f||a,n>_abab r1_aa(a,m)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += eri["abab_vvvo"]("e,f,a,n") * r1["aa_Lvo"]("X,a,m");

        // sigmar2_abab += +0.50 <j,i||a,n>_abab t2_abab(e,f,j,i) r1_aa(a,m) @ +0.50 <i,j||a,n>_abab t2_abab(e,f,i,j) r1_aa(a,m)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["39_aabb_vvvo"]("e,a,f,n") * r1["aa_Lvo"]("X,a,m");

        // sigmar2_abab += +1.00 <i,f||a,b>_bbbb t2_abab(e,a,m,i) r1_bb(b,n)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= reuse_tmps_["29_abba_vvvo"]("e,b,f,m") * r1["bb_Lvo"]("X,b,n");

        // sigmar2_abab += +1.00 f_bb(f,a) r2_abab(e,a,m,n)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += f["bb_vv"]("f,a") * r2["abab_Lvvoo"]("X,e,a,m,n");

        // sigmal1_aa += +1.00 <j,b||c,e>_aaaa l2_abab(m,i,b,a) t2_abab(c,a,j,i)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= reuse_tmps_["28_aabb_vvvo"]("e,b,a,i") * l2["abab_Loovv"]("X,m,i,b,a");

        // sigmal1_aa += +0.25 <k,j||e,i>_aaaa l2_aaaa(m,i,b,a) t2_aaaa(b,a,k,j)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += 0.50 * reuse_tmps_["31_aaaa_vvvo"]("a,b,e,i") * l2["aaaa_Loovv"]("X,m,i,b,a");

        // sigmal1_aa += +1.00 <b,j||e,c>_abab l2_aaaa(m,i,b,a) t2_abab(a,c,i,j)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += reuse_tmps_["26_aaaa_vvvo"]("a,e,b,i") * l2["aaaa_Loovv"]("X,m,i,b,a");

        // sigmal1_aa += -1.00 <j,b||e,c>_abab l2_abab(m,i,a,b) t2_abab(a,c,j,i)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += reuse_tmps_["27_aabb_vvvo"]("a,e,b,i") * l2["abab_Loovv"]("X,m,i,a,b");

        // sigmal1_aa += -1.00 <j,b||c,e>_aaaa l2_aaaa(m,i,b,a) t2_aaaa(c,a,i,j)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += reuse_tmps_["25_aaaa_vvvo"]("a,e,b,i") * l2["aaaa_Loovv"]("X,m,i,b,a");

        // sigmal1_aa += +0.50 <b,a||e,i>_abab l2_abab(m,i,b,a) @ +0.50 <a,b||e,i>_abab l2_abab(m,i,a,b)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += eri["abab_vvvo"]("b,a,e,i") * l2["abab_Loovv"]("X,m,i,b,a");

        // sigmal1_aa += +0.50 <b,a||e,i>_aaaa l2_aaaa(m,i,b,a)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += 0.50 * eri["aaaa_vvvo"]("b,a,e,i") * l2["aaaa_Loovv"]("X,m,i,b,a");

        // sigmal1_aa += -1.00 <b,j||e,c>_abab l2_abab(m,i,b,a) t2_bbbb(c,a,i,j)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= reuse_tmps_["37_aabb_vvvo"]("e,b,a,i") * l2["abab_Loovv"]("X,m,i,b,a");

        // sigmal1_aa += +0.25 <k,j||e,i>_abab l2_abab(m,i,b,a) t2_abab(b,a,k,j) @ +0.25 <j,k||e,i>_abab l2_abab(m,i,b,a) t2_abab(b,a,j,k) @ +0.25 <k,j||e,i>_abab l2_abab(m,i,a,b) t2_abab(a,b,k,j) @ +0.25 <j,k||e,i>_abab l2_abab(m,i,a,b) t2_abab(a,b,j,k)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += reuse_tmps_["39_aabb_vvvo"]("b,e,a,i") * l2["abab_Loovv"]("X,m,i,b,a");

        // sigmal1_bb += -1.00 <b,j||c,e>_abab l2_abab(i,m,b,a) t2_abab(c,a,i,j)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= reuse_tmps_["23_abba_vvvo"]("b,a,e,i") * l2["abab_Loovv"]("X,i,m,b,a");

        // sigmal1_bb += +0.25 <k,j||e,i>_bbbb l2_bbbb(m,i,b,a) t2_bbbb(b,a,k,j)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += 0.50 * reuse_tmps_["30_bbbb_vvvo"]("a,b,e,i") * l2["bbbb_Loovv"]("X,m,i,b,a");

        // sigmal1_bb += -1.00 <j,b||c,e>_abab l2_abab(i,m,a,b) t2_aaaa(c,a,i,j)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += reuse_tmps_["36_abba_vvvo"]("a,e,b,i") * l2["abab_Loovv"]("X,i,m,a,b");

        // sigmal1_bb += +1.00 <j,b||c,e>_bbbb l2_abab(i,m,a,b) t2_abab(a,c,i,j)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= reuse_tmps_["29_abba_vvvo"]("a,e,b,i") * l2["abab_Loovv"]("X,i,m,a,b");

        // sigmal1_bb += -1.00 <j,b||c,e>_bbbb l2_bbbb(m,i,b,a) t2_bbbb(c,a,i,j)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += reuse_tmps_["24_bbbb_vvvo"]("a,e,b,i") * l2["bbbb_Loovv"]("X,m,i,b,a");

        // sigmal1_bb += +0.50 <b,a||i,e>_abab l2_abab(i,m,b,a) @ +0.50 <a,b||i,e>_abab l2_abab(i,m,a,b)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= eri["abba_vvvo"]("b,a,e,i") * l2["abab_Loovv"]("X,i,m,b,a");

        // sigmal1_bb += +0.50 <b,a||e,i>_bbbb l2_bbbb(m,i,b,a)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += 0.50 * eri["bbbb_vvvo"]("b,a,e,i") * l2["bbbb_Loovv"]("X,m,i,b,a");

        // sigmal1_bb += +0.25 <k,j||i,e>_abab l2_abab(i,m,b,a) t2_abab(b,a,k,j) @ +0.25 <j,k||i,e>_abab l2_abab(i,m,b,a) t2_abab(b,a,j,k) @ +0.25 <k,j||i,e>_abab l2_abab(i,m,a,b) t2_abab(a,b,k,j) @ +0.25 <j,k||i,e>_abab l2_abab(i,m,a,b) t2_abab(a,b,j,k)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= reuse_tmps_["38_abba_vvvo"]("b,a,e,i") * l2["abab_Loovv"]("X,i,m,b,a");

        // sigmal1_bb += +1.00 <j,b||c,e>_abab l2_bbbb(m,i,b,a) t2_abab(c,a,j,i)  // flops: o1v1L1 += o2v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= reuse_tmps_["22_bbbb_vvvo"]("a,e,b,i") * l2["bbbb_Loovv"]("X,m,i,b,a");

        // sigmal2_abab += -0.50 <j,i||b,f>_bbbb l2_abab(m,n,e,a) t2_bbbb(b,a,j,i)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= reuse_tmps_["16_bb_vv"]("a,f") * l2["abab_Loovv"]("X,m,n,e,a");

        // sigmal2_abab += +1.00 f_bb(a,f) l2_abab(m,n,e,a)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += f["bb_vv"]("a,f") * l2["abab_Loovv"]("X,m,n,e,a");

        // sigmal2_abab += -0.50 <j,i||e,b>_abab l2_abab(m,n,a,f) t2_abab(a,b,j,i) @ -0.50 <i,j||e,b>_abab l2_abab(m,n,a,f) t2_abab(a,b,i,j)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= reuse_tmps_["14_aa_vv"]("a,e") * l2["abab_Loovv"]("X,m,n,a,f");

        // sigmal2_abab += -0.50 <j,i||b,f>_abab l2_abab(m,n,e,a) t2_abab(b,a,j,i) @ -0.50 <i,j||b,f>_abab l2_abab(m,n,e,a) t2_abab(b,a,i,j)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= reuse_tmps_["15_bb_vv"]("a,f") * l2["abab_Loovv"]("X,m,n,e,a");

        // sigmal2_abab += +1.00 <a,n||e,f>_abab l1_aa(m,a)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += eri["abab_vovv"]("a,n,e,f") * l1["aa_Lov"]("X,m,a");

        // sigmal2_abab += -0.50 <j,i||b,e>_aaaa l2_abab(m,n,a,f) t2_aaaa(b,a,j,i)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= reuse_tmps_["17_aa_vv"]("a,e") * l2["abab_Loovv"]("X,m,n,a,f");

        // sigmal2_abab += +1.00 f_aa(a,e) l2_abab(m,n,a,f)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += f["aa_vv"]("a,e") * l2["abab_Loovv"]("X,m,n,a,f");

        // sigmal2_abab += +1.00 <m,a||e,f>_abab l1_bb(n,a)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= eri["baab_vovv"]("a,m,e,f") * l1["bb_Lov"]("X,n,a");

        // sigmar2_aaaa += +0.25 <j,i||a,b>_aaaa t2_aaaa(a,b,m,n) r2_aaaa(e,f,j,i)  // flops: o2v2L1 += o4v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += reuse_tmps_["13_aaaa_oooo"]("n,m,i,j") * r2["aaaa_Lvvoo"]("X,e,f,j,i");

        // sigmar2_aaaa += +0.50 <j,i||m,n>_aaaa r2_aaaa(e,f,j,i)  // flops: o2v2L1 += o4v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += 0.50 * eri["aaaa_oooo"]("j,i,m,n") * r2["aaaa_Lvvoo"]("X,e,f,j,i");

        // sigmar2_abab += +0.25 <j,i||a,b>_abab t2_abab(a,b,m,n) r2_abab(e,f,j,i) @ +0.25 <i,j||a,b>_abab t2_abab(a,b,m,n) r2_abab(e,f,i,j) @ +0.25 <j,i||b,a>_abab t2_abab(b,a,m,n) r2_abab(e,f,j,i) @ +0.25 <i,j||b,a>_abab t2_abab(b,a,m,n) r2_abab(e,f,i,j)  // flops: o2v2L1 += o4v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["11_aabb_oooo"]("m,j,n,i") * r2["abab_Lvvoo"]("X,e,f,j,i");

        // sigmar2_abab += +0.50 <j,i||m,n>_abab r2_abab(e,f,j,i) @ +0.50 <i,j||m,n>_abab r2_abab(e,f,i,j)  // flops: o2v2L1 += o4v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += eri["abab_oooo"]("j,i,m,n") * r2["abab_Lvvoo"]("X,e,f,j,i");

        // sigmar2_bbbb += +0.50 <j,i||m,n>_bbbb r2_bbbb(e,f,j,i)  // flops: o2v2L1 += o4v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += 0.50 * eri["bbbb_oooo"]("j,i,m,n") * r2["bbbb_Lvvoo"]("X,e,f,j,i");

        // sigmar2_bbbb += +0.25 <j,i||a,b>_bbbb t2_bbbb(a,b,m,n) r2_bbbb(e,f,j,i)  // flops: o2v2L1 += o4v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += reuse_tmps_["12_bbbb_oooo"]("n,m,i,j") * r2["bbbb_Lvvoo"]("X,e,f,j,i");

        // sigmal2_aaaa += +0.25 <m,n||a,b>_aaaa l2_aaaa(i,j,e,f) t2_aaaa(a,b,i,j)  // flops: o2v2L1 += o4v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") += reuse_tmps_["13_aaaa_oooo"]("j,i,n,m") * l2["aaaa_Loovv"]("X,i,j,e,f");

        // sigmal2_aaaa += +0.50 <m,n||i,j>_aaaa l2_aaaa(i,j,e,f)  // flops: o2v2L1 += o4v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") += 0.50 * eri["aaaa_oooo"]("m,n,i,j") * l2["aaaa_Loovv"]("X,i,j,e,f");

        // sigmal2_abab += +0.25 <m,n||a,b>_abab l2_abab(i,j,e,f) t2_abab(a,b,i,j) @ +0.25 <m,n||b,a>_abab l2_abab(i,j,e,f) t2_abab(b,a,i,j) @ +0.25 <m,n||a,b>_abab l2_abab(j,i,e,f) t2_abab(a,b,j,i) @ +0.25 <m,n||b,a>_abab l2_abab(j,i,e,f) t2_abab(b,a,j,i)  // flops: o2v2L1 += o4v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += reuse_tmps_["11_aabb_oooo"]("i,m,j,n") * l2["abab_Loovv"]("X,i,j,e,f");

        // sigmal2_abab += +0.50 <m,n||i,j>_abab l2_abab(i,j,e,f) @ +0.50 <m,n||j,i>_abab l2_abab(j,i,e,f)  // flops: o2v2L1 += o4v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += eri["abab_oooo"]("m,n,i,j") * l2["abab_Loovv"]("X,i,j,e,f");

        // sigmal2_bbbb += +0.50 <m,n||i,j>_bbbb l2_bbbb(i,j,e,f)  // flops: o2v2L1 += o4v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") += 0.50 * eri["bbbb_oooo"]("m,n,i,j") * l2["bbbb_Loovv"]("X,i,j,e,f");

        // sigmal2_bbbb += +0.25 <m,n||a,b>_bbbb l2_bbbb(i,j,e,f) t2_bbbb(a,b,i,j)  // flops: o2v2L1 += o4v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") += reuse_tmps_["12_bbbb_oooo"]("j,i,n,m") * l2["bbbb_Loovv"]("X,i,j,e,f");

        // sigmar2_abab += +1.00 <i,e||a,m>_aaaa r2_abab(a,f,i,n)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= eri["aaaa_vovo"]("e,i,a,m") * r2["abab_Lvvoo"]("X,a,f,i,n");

        // sigmar2_abab += +1.00 <j,i||b,a>_abab t2_bbbb(a,f,n,i) r2_aaaa(b,e,m,j)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["5_abab_vvoo"]("b,f,j,n") * r2["aaaa_Lvvoo"]("X,b,e,m,j");

        // sigmar2_abab += +1.00 <j,i||a,b>_aaaa t2_aaaa(a,e,m,i) r2_abab(b,f,j,n)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["4_aaaa_vvoo"]("e,b,m,j") * r2["abab_Lvvoo"]("X,b,f,j,n");

        // sigmar2_abab += -1.00 <i,f||m,a>_abab r2_abab(e,a,i,n)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= eri["baba_vovo"]("f,i,a,m") * r2["abab_Lvvoo"]("X,e,a,i,n");

        // sigmar2_abab += +1.00 <i,j||a,b>_abab t2_aaaa(a,e,m,i) r2_bbbb(b,f,n,j)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["7_abab_vvoo"]("e,b,m,j") * r2["bbbb_Lvvoo"]("X,b,f,n,j");

        // sigmar2_abab += -1.00 <e,i||a,n>_abab r2_abab(a,f,m,i)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= eri["abab_vovo"]("e,i,a,n") * r2["abab_Lvvoo"]("X,a,f,m,i");

        // sigmar2_abab += +1.00 <i,f||a,n>_bbbb r2_abab(e,a,m,i)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= eri["bbbb_vovo"]("f,i,a,n") * r2["abab_Lvvoo"]("X,e,a,m,i");

        // sigmar2_abab += -1.00 <e,i||m,a>_abab r2_bbbb(a,f,n,i)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += eri["abba_vovo"]("e,i,a,m") * r2["bbbb_Lvvoo"]("X,a,f,n,i");

        // sigmar2_abab += +1.00 <i,j||b,a>_abab t2_abab(e,a,i,n) r2_abab(b,f,m,j)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["9_aabb_vvoo"]("e,b,n,j") * r2["abab_Lvvoo"]("X,b,f,m,j");

        // sigmar2_abab += +1.00 <j,i||a,b>_abab t2_abab(a,f,m,i) r2_abab(e,b,j,n)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["10_bbaa_vvoo"]("f,b,m,j") * r2["abab_Lvvoo"]("X,e,b,j,n");

        // sigmar2_abab += +1.00 <i,j||a,b>_abab t2_abab(a,f,i,n) r2_abab(e,b,m,j)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["1_bbbb_vvoo"]("f,b,n,j") * r2["abab_Lvvoo"]("X,e,b,m,j");

        // sigmar2_abab += +1.00 <j,i||b,a>_abab t2_abab(e,a,m,i) r2_abab(b,f,j,n)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["2_aaaa_vvoo"]("e,b,m,j") * r2["abab_Lvvoo"]("X,b,f,j,n");

        // sigmar2_abab += +1.00 <j,i||a,b>_bbbb t2_abab(e,a,m,i) r2_bbbb(b,f,n,j)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["8_abab_vvoo"]("e,b,m,j") * r2["bbbb_Lvvoo"]("X,b,f,n,j");

        // sigmar2_abab += +1.00 <j,i||a,b>_bbbb t2_bbbb(a,f,n,i) r2_abab(e,b,m,j)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["3_bbbb_vvoo"]("f,b,n,j") * r2["abab_Lvvoo"]("X,e,b,m,j");

        // sigmar2_abab += +1.00 <j,i||a,b>_aaaa t2_abab(a,f,i,n) r2_aaaa(b,e,m,j)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += reuse_tmps_["6_abab_vvoo"]("b,f,j,n") * r2["aaaa_Lvvoo"]("X,b,e,m,j");

        // sigmar2_abab += -1.00 <i,f||a,n>_abab r2_aaaa(a,e,m,i)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += eri["baab_vovo"]("f,i,a,n") * r2["aaaa_Lvvoo"]("X,a,e,m,i");

        // sigmal2_abab += +1.00 <n,j||b,f>_bbbb l2_abab(m,i,e,a) t2_bbbb(b,a,i,j)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += reuse_tmps_["3_bbbb_vvoo"]("a,f,i,n") * l2["abab_Loovv"]("X,m,i,e,a");

        // sigmal2_abab += +1.00 <m,j||b,e>_aaaa l2_bbbb(n,i,a,f) t2_abab(b,a,j,i)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += reuse_tmps_["6_abab_vvoo"]("e,a,m,i") * l2["bbbb_Loovv"]("X,n,i,a,f");

        // sigmal2_abab += +1.00 <m,j||e,b>_abab l2_abab(i,n,a,f) t2_abab(a,b,i,j)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += reuse_tmps_["2_aaaa_vvoo"]("a,e,i,m") * l2["abab_Loovv"]("X,i,n,a,f");

        // sigmal2_abab += +1.00 <j,n||b,f>_abab l2_abab(m,i,e,a) t2_abab(b,a,j,i)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += reuse_tmps_["1_bbbb_vvoo"]("a,f,i,n") * l2["abab_Loovv"]("X,m,i,e,a");

        // sigmal2_abab += +1.00 <j,n||e,b>_abab l2_abab(m,i,a,f) t2_abab(a,b,j,i)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += reuse_tmps_["9_aabb_vvoo"]("a,e,i,n") * l2["abab_Loovv"]("X,m,i,a,f");

        // sigmal2_abab += +1.00 <m,j||b,f>_abab l2_abab(i,n,e,a) t2_abab(b,a,i,j)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += reuse_tmps_["10_bbaa_vvoo"]("a,f,i,m") * l2["abab_Loovv"]("X,i,n,e,a");

        // sigmal2_abab += -1.00 <m,a||i,f>_abab l2_abab(i,n,e,a)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= eri["baba_vovo"]("a,m,f,i") * l2["abab_Loovv"]("X,i,n,e,a");

        // sigmal2_abab += +1.00 <n,a||f,i>_bbbb l2_abab(m,i,e,a)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= eri["bbbb_vovo"]("a,n,f,i") * l2["abab_Loovv"]("X,m,i,e,a");

        // sigmal2_abab += +1.00 <m,j||b,e>_aaaa l2_abab(i,n,a,f) t2_aaaa(b,a,i,j)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += reuse_tmps_["4_aaaa_vvoo"]("a,e,i,m") * l2["abab_Loovv"]("X,i,n,a,f");

        // sigmal2_abab += +1.00 <n,j||b,f>_bbbb l2_aaaa(m,i,a,e) t2_abab(a,b,i,j)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += reuse_tmps_["8_abab_vvoo"]("a,f,i,n") * l2["aaaa_Loovv"]("X,m,i,a,e");

        // sigmal2_abab += +1.00 <j,n||b,f>_abab l2_aaaa(m,i,a,e) t2_aaaa(b,a,i,j)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += reuse_tmps_["7_abab_vvoo"]("a,f,i,n") * l2["aaaa_Loovv"]("X,m,i,a,e");

        // sigmal2_abab += +1.00 <m,j||e,b>_abab l2_bbbb(n,i,a,f) t2_bbbb(b,a,i,j)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += reuse_tmps_["5_abab_vvoo"]("e,a,m,i") * l2["bbbb_Loovv"]("X,n,i,a,f");

        // sigmal2_abab += -1.00 <a,n||i,f>_abab l2_aaaa(m,i,a,e)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += eri["abba_vovo"]("a,n,f,i") * l2["aaaa_Loovv"]("X,m,i,a,e");

        // sigmal2_abab += -1.00 <a,n||e,i>_abab l2_abab(m,i,a,f)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= eri["abab_vovo"]("a,n,e,i") * l2["abab_Loovv"]("X,m,i,a,f");

        // sigmal2_abab += -1.00 <m,a||e,i>_abab l2_bbbb(n,i,a,f)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += eri["baab_vovo"]("a,m,e,i") * l2["bbbb_Loovv"]("X,n,i,a,f");

        // sigmal2_abab += +1.00 <m,a||e,i>_aaaa l2_abab(i,n,a,f)  // flops: o2v2L1 += o3v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= eri["aaaa_vovo"]("a,m,e,i") * l2["abab_Loovv"]("X,i,n,a,f");

        // sigmar2_aaaa += +0.50 <e,f||a,b>_aaaa r2_aaaa(a,b,m,n)  // flops: o2v2L1 += o2v4L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += 0.50 * eri["aaaa_vvvv"]("e,f,a,b") * r2["aaaa_Lvvoo"]("X,a,b,m,n");

        // sigmar2_abab += +0.50 <e,f||a,b>_abab r2_abab(a,b,m,n) @ +0.50 <e,f||b,a>_abab r2_abab(b,a,m,n)  // flops: o2v2L1 += o2v4L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += eri["abab_vvvv"]("e,f,a,b") * r2["abab_Lvvoo"]("X,a,b,m,n");

        // sigmar2_bbbb += +0.50 <e,f||a,b>_bbbb r2_bbbb(a,b,m,n)  // flops: o2v2L1 += o2v4L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += 0.50 * eri["bbbb_vvvv"]("e,f,a,b") * r2["bbbb_Lvvoo"]("X,a,b,m,n");

        // sigmal2_aaaa += +0.50 <b,a||e,f>_aaaa l2_aaaa(m,n,b,a)  // flops: o2v2L1 += o2v4L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") += 0.50 * eri["aaaa_vvvv"]("b,a,e,f") * l2["aaaa_Loovv"]("X,m,n,b,a");

        // sigmal2_abab += +0.50 <b,a||e,f>_abab l2_abab(m,n,b,a) @ +0.50 <a,b||e,f>_abab l2_abab(m,n,a,b)  // flops: o2v2L1 += o2v4L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += eri["abab_vvvv"]("b,a,e,f") * l2["abab_Loovv"]("X,m,n,b,a");

        // sigmal2_bbbb += +0.50 <b,a||e,f>_bbbb l2_bbbb(m,n,b,a)  // flops: o2v2L1 += o2v4L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") += 0.50 * eri["bbbb_vvvv"]("b,a,e,f") * l2["bbbb_Loovv"]("X,m,n,b,a");

        // sigmar2_aaaa += +0.25 <j,i||a,b>_aaaa t2_aaaa(e,f,j,i) r2_aaaa(a,b,m,n)  // flops: o2v2L1 += o4v2L1 o4v2L1 | mem: o2v2L1 += o4v0L1 o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += 0.25 * eri["aaaa_oovv"]("j,i,a,b") * r2["aaaa_Lvvoo"]("X,a,b,m,n") * t2["aaaa_vvoo"]("e,f,j,i");

        // sigmar2_abab += +0.25 <j,i||a,b>_abab t2_abab(e,f,j,i) r2_abab(a,b,m,n) @ +0.25 <j,i||b,a>_abab t2_abab(e,f,j,i) r2_abab(b,a,m,n) @ +0.25 <i,j||a,b>_abab t2_abab(e,f,i,j) r2_abab(a,b,m,n) @ +0.25 <i,j||b,a>_abab t2_abab(e,f,i,j) r2_abab(b,a,m,n)  // flops: o2v2L1 += o4v2L1 o4v2L1 | mem: o2v2L1 += o4v0L1 o2v2L1
        sigmar2_abab("X,e,f,m,n") += eri["abab_oovv"]("j,i,a,b") * r2["abab_Lvvoo"]("X,a,b,m,n") * t2["abab_vvoo"]("e,f,j,i");

        // sigmar2_bbbb += +0.25 <j,i||a,b>_bbbb t2_bbbb(e,f,j,i) r2_bbbb(a,b,m,n)  // flops: o2v2L1 += o4v2L1 o4v2L1 | mem: o2v2L1 += o4v0L1 o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += 0.25 * eri["bbbb_oovv"]("j,i,a,b") * r2["bbbb_Lvvoo"]("X,a,b,m,n") * t2["bbbb_vvoo"]("e,f,j,i");

        // sigmal2_aaaa += +0.25 <j,i||e,f>_aaaa l2_aaaa(m,n,b,a) t2_aaaa(b,a,j,i)  // flops: o2v2L1 += o4v2L1 o4v2L1 | mem: o2v2L1 += o4v0L1 o2v2L1
        sigmal2_aaaa("X,e,f,m,n") += 0.25 * l2["aaaa_Loovv"]("X,m,n,b,a") * t2["aaaa_vvoo"]("b,a,j,i") * eri["aaaa_oovv"]("j,i,e,f");

        // sigmal2_abab += +0.25 <j,i||e,f>_abab l2_abab(m,n,b,a) t2_abab(b,a,j,i) @ +0.25 <i,j||e,f>_abab l2_abab(m,n,b,a) t2_abab(b,a,i,j) @ +0.25 <j,i||e,f>_abab l2_abab(m,n,a,b) t2_abab(a,b,j,i) @ +0.25 <i,j||e,f>_abab l2_abab(m,n,a,b) t2_abab(a,b,i,j)  // flops: o2v2L1 += o4v2L1 o4v2L1 | mem: o2v2L1 += o4v0L1 o2v2L1
        sigmal2_abab("X,e,f,m,n") += l2["abab_Loovv"]("X,m,n,b,a") * t2["abab_vvoo"]("b,a,j,i") * eri["abab_oovv"]("j,i,e,f");

        // sigmal2_bbbb += +0.25 <j,i||e,f>_bbbb l2_bbbb(m,n,b,a) t2_bbbb(b,a,j,i)  // flops: o2v2L1 += o4v2L1 o4v2L1 | mem: o2v2L1 += o4v0L1 o2v2L1
        sigmal2_bbbb("X,e,f,m,n") += 0.25 * l2["bbbb_Loovv"]("X,m,n,b,a") * t2["bbbb_vvoo"]("b,a,j,i") * eri["bbbb_oovv"]("j,i,e,f");

        // tmps_[1_bbbb_Lvvoo](X,f,e,m,n) = 1.00 l2[abab_Loovv](X,i,n,a,e) * eri[abab_oovv](j,m,b,f) * t2[aaaa_vvoo](b,a,i,j) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["1_bbbb_Lvvoo"]("X,f,e,m,n")  = l2["abab_Loovv"]("X,i,n,a,e") * reuse_tmps_["7_abab_vvoo"]("a,f,i,m");
        sigmal2_bbbb("X,e,f,m,n") += tmps_["1_bbbb_Lvvoo"]("X,f,e,m,n");
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["1_bbbb_Lvvoo"]("X,e,f,m,n");
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["1_bbbb_Lvvoo"]("X,f,e,n,m");

        // sigmal2_bbbb += +1.00 P(m,n) P(e,f) <j,n||b,e>_abab l2_abab(i,m,a,f) t2_aaaa(b,a,i,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") += tmps_["1_bbbb_Lvvoo"]("X,e,f,n,m");
        tmps_["1_bbbb_Lvvoo"].~TArrayD();

        // tmps_[2_bbbb_Lvvoo](X,e,f,n,m) = 1.00 eri[bbbb_oovv](m,j,b,f) * t2[bbbb_vvoo](b,a,i,j) * l2[bbbb_Loovv](X,n,i,a,e) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["2_bbbb_Lvvoo"]("X,e,f,n,m")  = reuse_tmps_["3_bbbb_vvoo"]("a,f,i,m") * l2["bbbb_Loovv"]("X,n,i,a,e");

        // sigmal2_bbbb += +1.00 P(m,n) P(e,f) <n,j||b,e>_bbbb l2_bbbb(m,i,a,f) t2_bbbb(b,a,i,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") += tmps_["2_bbbb_Lvvoo"]("X,f,e,m,n");
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["2_bbbb_Lvvoo"]("X,f,e,n,m");
        sigmal2_bbbb("X,e,f,m,n") += tmps_["2_bbbb_Lvvoo"]("X,e,f,n,m");
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["2_bbbb_Lvvoo"]("X,e,f,m,n");
        tmps_["2_bbbb_Lvvoo"].~TArrayD();

        // tmps_[3_bbbb_Lvvoo](X,f,e,n,m) = 1.00 eri[bbbb_vovo](a,m,e,i) * l2[bbbb_Loovv](X,n,i,a,f) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["3_bbbb_Lvvoo"]("X,f,e,n,m")  = eri["bbbb_vovo"]("a,m,e,i") * l2["bbbb_Loovv"]("X,n,i,a,f");

        // sigmal2_bbbb += +1.00 P(m,n) P(e,f) <n,a||e,i>_bbbb l2_bbbb(m,i,a,f)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["3_bbbb_Lvvoo"]("X,f,e,m,n");
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["3_bbbb_Lvvoo"]("X,e,f,n,m");
        sigmal2_bbbb("X,e,f,m,n") += tmps_["3_bbbb_Lvvoo"]("X,f,e,n,m");
        sigmal2_bbbb("X,e,f,m,n") += tmps_["3_bbbb_Lvvoo"]("X,e,f,m,n");
        tmps_["3_bbbb_Lvvoo"].~TArrayD();

        // tmps_[4_bbbb_Lvvoo](X,f,e,m,n) = 1.00 l2[abab_Loovv](X,i,n,a,e) * eri[bbbb_oovv](m,j,b,f) * t2[abab_vvoo](a,b,i,j) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["4_bbbb_Lvvoo"]("X,f,e,m,n")  = l2["abab_Loovv"]("X,i,n,a,e") * reuse_tmps_["8_abab_vvoo"]("a,f,i,m");
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["4_bbbb_Lvvoo"]("X,e,f,m,n");
        sigmal2_bbbb("X,e,f,m,n") += tmps_["4_bbbb_Lvvoo"]("X,f,e,m,n");
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["4_bbbb_Lvvoo"]("X,f,e,n,m");

        // sigmal2_bbbb += +1.00 P(m,n) P(e,f) <n,j||b,e>_bbbb l2_abab(i,m,a,f) t2_abab(a,b,i,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") += tmps_["4_bbbb_Lvvoo"]("X,e,f,n,m");
        tmps_["4_bbbb_Lvvoo"].~TArrayD();

        // tmps_[5_bbbb_Lvvoo](X,f,e,m,n) = 1.00 eri[abba_vovo](a,n,e,i) * l2[abab_Loovv](X,i,m,a,f) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["5_bbbb_Lvvoo"]("X,f,e,m,n")  = eri["abba_vovo"]("a,n,e,i") * l2["abab_Loovv"]("X,i,m,a,f");
        sigmal2_bbbb("X,e,f,m,n") += tmps_["5_bbbb_Lvvoo"]("X,e,f,n,m");

        // sigmal2_bbbb += -1.00 P(m,n) P(e,f) <a,n||i,e>_abab l2_abab(i,m,a,f)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") += tmps_["5_bbbb_Lvvoo"]("X,f,e,m,n");
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["5_bbbb_Lvvoo"]("X,f,e,n,m");
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["5_bbbb_Lvvoo"]("X,e,f,m,n");
        tmps_["5_bbbb_Lvvoo"].~TArrayD();

        // tmps_[6_bbbb_Lvvoo](X,e,f,n,m) = 1.00 eri[abab_oovv](j,m,b,f) * t2[abab_vvoo](b,a,j,i) * l2[bbbb_Loovv](X,n,i,a,e) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["6_bbbb_Lvvoo"]("X,e,f,n,m")  = reuse_tmps_["1_bbbb_vvoo"]("a,f,i,m") * l2["bbbb_Loovv"]("X,n,i,a,e");
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["6_bbbb_Lvvoo"]("X,e,f,m,n");
        sigmal2_bbbb("X,e,f,m,n") += tmps_["6_bbbb_Lvvoo"]("X,e,f,n,m");

        // sigmal2_bbbb += +1.00 P(m,n) P(e,f) <j,n||b,e>_abab l2_bbbb(m,i,a,f) t2_abab(b,a,j,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") += tmps_["6_bbbb_Lvvoo"]("X,f,e,m,n");
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["6_bbbb_Lvvoo"]("X,f,e,n,m");
        tmps_["6_bbbb_Lvvoo"].~TArrayD();

        // tmps_[7_aaaa_Lvvoo](X,e,f,n,m) = 1.00 eri[aaaa_vovo](a,m,f,i) * l2[aaaa_Loovv](X,n,i,a,e) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["7_aaaa_Lvvoo"]("X,e,f,n,m")  = eri["aaaa_vovo"]("a,m,f,i") * l2["aaaa_Loovv"]("X,n,i,a,e");
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["7_aaaa_Lvvoo"]("X,e,f,n,m");
        sigmal2_aaaa("X,e,f,m,n") += tmps_["7_aaaa_Lvvoo"]("X,e,f,m,n");
        sigmal2_aaaa("X,e,f,m,n") += tmps_["7_aaaa_Lvvoo"]("X,f,e,n,m");

        // sigmal2_aaaa += +1.00 P(m,n) P(e,f) <n,a||e,i>_aaaa l2_aaaa(m,i,a,f)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["7_aaaa_Lvvoo"]("X,f,e,m,n");
        tmps_["7_aaaa_Lvvoo"].~TArrayD();

        // tmps_[8_aaaa_Lvvoo](X,f,e,n,m) = 1.00 eri[baab_vovo](a,m,e,i) * l2[abab_Loovv](X,n,i,f,a) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["8_aaaa_Lvvoo"]("X,f,e,n,m")  = eri["baab_vovo"]("a,m,e,i") * l2["abab_Loovv"]("X,n,i,f,a");
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["8_aaaa_Lvvoo"]("X,f,e,n,m");

        // sigmal2_aaaa += -1.00 P(m,n) P(e,f) <n,a||e,i>_abab l2_abab(m,i,f,a)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") += tmps_["8_aaaa_Lvvoo"]("X,f,e,m,n");
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["8_aaaa_Lvvoo"]("X,e,f,m,n");
        sigmal2_aaaa("X,e,f,m,n") += tmps_["8_aaaa_Lvvoo"]("X,e,f,n,m");
        tmps_["8_aaaa_Lvvoo"].~TArrayD();

        // tmps_[9_aaaa_Lvvoo](X,e,f,m,n) = 1.00 l2[aaaa_Loovv](X,n,i,a,f) * eri[aaaa_oovv](m,j,b,e) * t2[aaaa_vvoo](b,a,i,j) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["9_aaaa_Lvvoo"]("X,e,f,m,n")  = l2["aaaa_Loovv"]("X,n,i,a,f") * reuse_tmps_["4_aaaa_vvoo"]("a,e,i,m");
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["9_aaaa_Lvvoo"]("X,f,e,n,m");
        sigmal2_aaaa("X,e,f,m,n") += tmps_["9_aaaa_Lvvoo"]("X,f,e,m,n");

        // sigmal2_aaaa += +1.00 P(m,n) P(e,f) <n,j||b,e>_aaaa l2_aaaa(m,i,a,f) t2_aaaa(b,a,i,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") += tmps_["9_aaaa_Lvvoo"]("X,e,f,n,m");
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["9_aaaa_Lvvoo"]("X,e,f,m,n");
        tmps_["9_aaaa_Lvvoo"].~TArrayD();

        // tmps_[10_aaaa_Lvvoo](X,f,e,m,n) = 1.00 l2[abab_Loovv](X,n,i,e,a) * eri[abab_oovv](m,j,f,b) * t2[bbbb_vvoo](b,a,i,j) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["10_aaaa_Lvvoo"]("X,f,e,m,n")  = l2["abab_Loovv"]("X,n,i,e,a") * reuse_tmps_["5_abab_vvoo"]("f,a,m,i");
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["10_aaaa_Lvvoo"]("X,e,f,m,n");
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["10_aaaa_Lvvoo"]("X,f,e,n,m");

        // sigmal2_aaaa += +1.00 P(m,n) P(e,f) <n,j||e,b>_abab l2_abab(m,i,f,a) t2_bbbb(b,a,i,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") += tmps_["10_aaaa_Lvvoo"]("X,e,f,n,m");
        sigmal2_aaaa("X,e,f,m,n") += tmps_["10_aaaa_Lvvoo"]("X,f,e,m,n");
        tmps_["10_aaaa_Lvvoo"].~TArrayD();

        // tmps_[11_aaaa_Lvvoo](X,f,e,m,n) = 1.00 l2[aaaa_Loovv](X,n,i,a,e) * eri[abab_oovv](m,j,f,b) * t2[abab_vvoo](a,b,i,j) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["11_aaaa_Lvvoo"]("X,f,e,m,n")  = l2["aaaa_Loovv"]("X,n,i,a,e") * reuse_tmps_["2_aaaa_vvoo"]("a,f,i,m");
        sigmal2_aaaa("X,e,f,m,n") += tmps_["11_aaaa_Lvvoo"]("X,f,e,m,n");
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["11_aaaa_Lvvoo"]("X,e,f,m,n");
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["11_aaaa_Lvvoo"]("X,f,e,n,m");

        // sigmal2_aaaa += +1.00 P(m,n) P(e,f) <n,j||e,b>_abab l2_aaaa(m,i,a,f) t2_abab(a,b,i,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") += tmps_["11_aaaa_Lvvoo"]("X,e,f,n,m");
        tmps_["11_aaaa_Lvvoo"].~TArrayD();

        // tmps_[12_aaaa_Lvvoo](X,f,e,m,n) = 1.00 l2[abab_Loovv](X,n,i,e,a) * eri[aaaa_oovv](m,j,b,f) * t2[abab_vvoo](b,a,j,i) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["12_aaaa_Lvvoo"]("X,f,e,m,n")  = l2["abab_Loovv"]("X,n,i,e,a") * reuse_tmps_["6_abab_vvoo"]("f,a,m,i");

        // sigmal2_aaaa += +1.00 P(m,n) P(e,f) <n,j||b,e>_aaaa l2_abab(m,i,f,a) t2_abab(b,a,j,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") += tmps_["12_aaaa_Lvvoo"]("X,e,f,n,m");
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["12_aaaa_Lvvoo"]("X,e,f,m,n");
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["12_aaaa_Lvvoo"]("X,f,e,n,m");
        sigmal2_aaaa("X,e,f,m,n") += tmps_["12_aaaa_Lvvoo"]("X,f,e,m,n");
        tmps_["12_aaaa_Lvvoo"].~TArrayD();

        // tmps_[13_bbbb_Lvvoo](X,f,e,m,n) = 1.00 r2[bbbb_Lvvoo](X,b,e,n,j) * eri[abab_oovv](i,j,a,b) * t2[abab_vvoo](a,f,i,m) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["13_bbbb_Lvvoo"]("X,f,e,m,n")  = r2["bbbb_Lvvoo"]("X,b,e,n,j") * reuse_tmps_["1_bbbb_vvoo"]("f,b,m,j");

        // sigmar2_bbbb += +1.00 P(m,n) P(e,f) <i,j||a,b>_abab t2_abab(a,e,i,n) r2_bbbb(b,f,m,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += tmps_["13_bbbb_Lvvoo"]("X,e,f,n,m");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["13_bbbb_Lvvoo"]("X,f,e,n,m");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["13_bbbb_Lvvoo"]("X,e,f,m,n");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["13_bbbb_Lvvoo"]("X,f,e,m,n");
        tmps_["13_bbbb_Lvvoo"].~TArrayD();

        // tmps_[14_bbbb_Lvvoo](X,f,e,m,n) = 1.00 eri[bbbb_vovo](e,i,a,n) * r2[bbbb_Lvvoo](X,a,f,m,i) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["14_bbbb_Lvvoo"]("X,f,e,m,n")  = eri["bbbb_vovo"]("e,i,a,n") * r2["bbbb_Lvvoo"]("X,a,f,m,i");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["14_bbbb_Lvvoo"]("X,e,f,n,m");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["14_bbbb_Lvvoo"]("X,f,e,n,m");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["14_bbbb_Lvvoo"]("X,e,f,m,n");

        // sigmar2_bbbb += +1.00 P(m,n) P(e,f) <i,e||a,n>_bbbb r2_bbbb(a,f,m,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["14_bbbb_Lvvoo"]("X,f,e,m,n");
        tmps_["14_bbbb_Lvvoo"].~TArrayD();

        // tmps_[15_bbbb_Lvvoo](X,e,f,m,n) = 1.00 r2[bbbb_Lvvoo](X,b,f,n,j) * eri[bbbb_oovv](j,i,a,b) * t2[bbbb_vvoo](a,e,m,i) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["15_bbbb_Lvvoo"]("X,e,f,m,n")  = r2["bbbb_Lvvoo"]("X,b,f,n,j") * reuse_tmps_["3_bbbb_vvoo"]("e,b,m,j");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["15_bbbb_Lvvoo"]("X,f,e,n,m");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["15_bbbb_Lvvoo"]("X,e,f,m,n");

        // sigmar2_bbbb += +1.00 P(m,n) P(e,f) <j,i||a,b>_bbbb t2_bbbb(a,e,n,i) r2_bbbb(b,f,m,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += tmps_["15_bbbb_Lvvoo"]("X,e,f,n,m");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["15_bbbb_Lvvoo"]("X,f,e,m,n");
        tmps_["15_bbbb_Lvvoo"].~TArrayD();

        // tmps_[16_bbbb_Lvvoo](X,e,f,n,m) = 1.00 eri[baab_vovo](f,i,a,m) * r2[abab_Lvvoo](X,a,e,i,n) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["16_bbbb_Lvvoo"]("X,e,f,n,m")  = eri["baab_vovo"]("f,i,a,m") * r2["abab_Lvvoo"]("X,a,e,i,n");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["16_bbbb_Lvvoo"]("X,e,f,n,m");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["16_bbbb_Lvvoo"]("X,e,f,m,n");

        // sigmar2_bbbb += -1.00 P(m,n) P(e,f) <i,e||a,n>_abab r2_abab(a,f,i,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += tmps_["16_bbbb_Lvvoo"]("X,f,e,m,n");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["16_bbbb_Lvvoo"]("X,f,e,n,m");
        tmps_["16_bbbb_Lvvoo"].~TArrayD();

        // tmps_[17_bbbb_Lvvoo](X,e,f,n,m) = 1.00 r2[abab_Lvvoo](X,b,f,j,m) * eri[abab_oovv](j,i,b,a) * t2[bbbb_vvoo](a,e,n,i) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["17_bbbb_Lvvoo"]("X,e,f,n,m")  = r2["abab_Lvvoo"]("X,b,f,j,m") * reuse_tmps_["5_abab_vvoo"]("b,e,j,n");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["17_bbbb_Lvvoo"]("X,f,e,m,n");

        // sigmar2_bbbb += +1.00 P(m,n) P(e,f) <j,i||b,a>_abab t2_bbbb(a,e,n,i) r2_abab(b,f,j,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += tmps_["17_bbbb_Lvvoo"]("X,e,f,n,m");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["17_bbbb_Lvvoo"]("X,f,e,n,m");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["17_bbbb_Lvvoo"]("X,e,f,m,n");
        tmps_["17_bbbb_Lvvoo"].~TArrayD();

        // tmps_[18_bbbb_Lvvoo](X,e,f,n,m) = 1.00 r2[abab_Lvvoo](X,b,f,j,m) * eri[aaaa_oovv](j,i,a,b) * t2[abab_vvoo](a,e,i,n) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["18_bbbb_Lvvoo"]("X,e,f,n,m")  = r2["abab_Lvvoo"]("X,b,f,j,m") * reuse_tmps_["6_abab_vvoo"]("b,e,j,n");

        // sigmar2_bbbb += +1.00 P(m,n) P(e,f) <j,i||a,b>_aaaa t2_abab(a,e,i,n) r2_abab(b,f,j,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += tmps_["18_bbbb_Lvvoo"]("X,e,f,n,m");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["18_bbbb_Lvvoo"]("X,f,e,m,n");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["18_bbbb_Lvvoo"]("X,e,f,m,n");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["18_bbbb_Lvvoo"]("X,f,e,n,m");
        tmps_["18_bbbb_Lvvoo"].~TArrayD();

        // tmps_[19_aaaa_Lvvoo](X,e,f,m,n) = 1.00 eri[abba_vovo](f,i,a,n) * r2[abab_Lvvoo](X,e,a,m,i) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["19_aaaa_Lvvoo"]("X,e,f,m,n")  = eri["abba_vovo"]("f,i,a,n") * r2["abab_Lvvoo"]("X,e,a,m,i");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["19_aaaa_Lvvoo"]("X,e,f,m,n");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["19_aaaa_Lvvoo"]("X,f,e,n,m");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["19_aaaa_Lvvoo"]("X,e,f,n,m");

        // sigmar2_aaaa += -1.00 P(m,n) P(e,f) <e,i||n,a>_abab r2_abab(f,a,m,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += tmps_["19_aaaa_Lvvoo"]("X,f,e,m,n");
        tmps_["19_aaaa_Lvvoo"].~TArrayD();

        // tmps_[20_aaaa_Lvvoo](X,e,f,n,m) = 1.00 r2[aaaa_Lvvoo](X,b,f,m,j) * eri[aaaa_oovv](j,i,a,b) * t2[aaaa_vvoo](a,e,n,i) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["20_aaaa_Lvvoo"]("X,e,f,n,m")  = r2["aaaa_Lvvoo"]("X,b,f,m,j") * reuse_tmps_["4_aaaa_vvoo"]("e,b,n,j");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["20_aaaa_Lvvoo"]("X,f,e,m,n");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["20_aaaa_Lvvoo"]("X,e,f,m,n");

        // sigmar2_aaaa += +1.00 P(m,n) P(e,f) <j,i||a,b>_aaaa t2_aaaa(a,e,n,i) r2_aaaa(b,f,m,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += tmps_["20_aaaa_Lvvoo"]("X,e,f,n,m");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["20_aaaa_Lvvoo"]("X,f,e,n,m");
        tmps_["20_aaaa_Lvvoo"].~TArrayD();

        // tmps_[21_aaaa_Lvvoo](X,f,e,m,n) = 1.00 r2[abab_Lvvoo](X,e,b,n,j) * eri[bbbb_oovv](j,i,a,b) * t2[abab_vvoo](f,a,m,i) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["21_aaaa_Lvvoo"]("X,f,e,m,n")  = r2["abab_Lvvoo"]("X,e,b,n,j") * reuse_tmps_["8_abab_vvoo"]("f,b,m,j");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["21_aaaa_Lvvoo"]("X,e,f,m,n");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["21_aaaa_Lvvoo"]("X,f,e,n,m");

        // sigmar2_aaaa += +1.00 P(m,n) P(e,f) <j,i||a,b>_bbbb t2_abab(e,a,n,i) r2_abab(f,b,m,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += tmps_["21_aaaa_Lvvoo"]("X,e,f,n,m");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["21_aaaa_Lvvoo"]("X,f,e,m,n");
        tmps_["21_aaaa_Lvvoo"].~TArrayD();

        // tmps_[22_aaaa_Lvvoo](X,f,e,m,n) = 1.00 r2[aaaa_Lvvoo](X,b,e,n,j) * eri[abab_oovv](j,i,b,a) * t2[abab_vvoo](f,a,m,i) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["22_aaaa_Lvvoo"]("X,f,e,m,n")  = r2["aaaa_Lvvoo"]("X,b,e,n,j") * reuse_tmps_["2_aaaa_vvoo"]("f,b,m,j");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["22_aaaa_Lvvoo"]("X,f,e,n,m");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["22_aaaa_Lvvoo"]("X,f,e,m,n");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["22_aaaa_Lvvoo"]("X,e,f,m,n");

        // sigmar2_aaaa += +1.00 P(m,n) P(e,f) <j,i||b,a>_abab t2_abab(e,a,n,i) r2_aaaa(b,f,m,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += tmps_["22_aaaa_Lvvoo"]("X,e,f,n,m");
        tmps_["22_aaaa_Lvvoo"].~TArrayD();

        // tmps_[23_aaaa_Lvvoo](X,e,f,m,n) = 1.00 eri[aaaa_vovo](f,i,a,n) * r2[aaaa_Lvvoo](X,a,e,m,i) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["23_aaaa_Lvvoo"]("X,e,f,m,n")  = eri["aaaa_vovo"]("f,i,a,n") * r2["aaaa_Lvvoo"]("X,a,e,m,i");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["23_aaaa_Lvvoo"]("X,e,f,n,m");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["23_aaaa_Lvvoo"]("X,e,f,m,n");

        // sigmar2_aaaa += +1.00 P(m,n) P(e,f) <i,e||a,n>_aaaa r2_aaaa(a,f,m,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["23_aaaa_Lvvoo"]("X,f,e,m,n");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["23_aaaa_Lvvoo"]("X,f,e,n,m");
        tmps_["23_aaaa_Lvvoo"].~TArrayD();

        // tmps_[24_aaaa_Lvvoo](X,f,e,m,n) = 1.00 r2[abab_Lvvoo](X,e,b,n,j) * eri[abab_oovv](i,j,a,b) * t2[aaaa_vvoo](a,f,m,i) // flops: o2v2L1 = o3v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["24_aaaa_Lvvoo"]("X,f,e,m,n")  = r2["abab_Lvvoo"]("X,e,b,n,j") * reuse_tmps_["7_abab_vvoo"]("f,b,m,j");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["24_aaaa_Lvvoo"]("X,e,f,m,n");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["24_aaaa_Lvvoo"]("X,f,e,n,m");

        // sigmar2_aaaa += +1.00 P(m,n) P(e,f) <i,j||a,b>_abab t2_aaaa(a,e,n,i) r2_abab(f,b,m,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += tmps_["24_aaaa_Lvvoo"]("X,e,f,n,m");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["24_aaaa_Lvvoo"]("X,f,e,m,n");
        tmps_["24_aaaa_Lvvoo"].~TArrayD();

        // tmps_[25_aa_Lvv](X,b,f) = 1.00 l2[abab_Loovv](X,i,j,f,a) * t2[abab_vvoo](b,a,i,j) // flops: o0v2L1 = o2v3L1 | mem: o0v2L1 = o0v2L1
        tmps_["25_aa_Lvv"]("X,b,f")  = l2["abab_Loovv"]("X,i,j,f,a") * t2["abab_vvoo"]("b,a,i,j");

        // sigmal1_aa += -0.50 <m,b||c,e>_aaaa l2_abab(i,j,b,a) t2_abab(c,a,i,j) @ -0.50 <m,b||c,e>_aaaa l2_abab(j,i,b,a) t2_abab(c,a,j,i)  // flops: o1v1L1 += o1v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += tmps_["25_aa_Lvv"]("X,c,b") * eri["aaaa_vovv"]("b,m,c,e");

        // sigmal1_bb += +0.50 <b,m||c,e>_abab l2_abab(i,j,b,a) t2_abab(c,a,i,j) @ +0.50 <b,m||c,e>_abab l2_abab(j,i,b,a) t2_abab(c,a,j,i)  // flops: o1v1L1 += o1v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += tmps_["25_aa_Lvv"]("X,c,b") * eri["abab_vovv"]("b,m,c,e");

        // sigmal2_abab += -0.50 <m,n||b,f>_abab l2_abab(i,j,e,a) t2_abab(b,a,i,j) @ -0.50 <m,n||b,f>_abab l2_abab(j,i,e,a) t2_abab(b,a,j,i)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= tmps_["25_aa_Lvv"]("X,b,e") * eri["abab_oovv"]("m,n,b,f");

        // tmps_[104_aaaa_Lvvoo](X,f,e,n,m) = 1.00 l2[abab_Loovv](X,i,j,e,a) * t2[abab_vvoo](b,a,i,j) * eri[aaaa_oovv](m,n,b,f) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["104_aaaa_Lvvoo"]("X,f,e,n,m")  = tmps_["25_aa_Lvv"]("X,b,e") * eri["aaaa_oovv"]("m,n,b,f");
        tmps_["25_aa_Lvv"].~TArrayD();

        // sigmal2_aaaa += +0.50 P(e,f) <m,n||b,e>_aaaa l2_abab(i,j,f,a) t2_abab(b,a,i,j) @ +0.50 P(e,f) <m,n||b,e>_aaaa l2_abab(j,i,f,a) t2_abab(b,a,j,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") += tmps_["104_aaaa_Lvvoo"]("X,e,f,n,m");
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["104_aaaa_Lvvoo"]("X,f,e,n,m");
        tmps_["104_aaaa_Lvvoo"].~TArrayD();

        // tmps_[26_bb_Lvv](X,c,b) = 1.00 l2[abab_Loovv](X,i,j,a,b) * t2[abab_vvoo](a,c,i,j) // flops: o0v2L1 = o2v3L1 | mem: o0v2L1 = o0v2L1
        tmps_["26_bb_Lvv"]("X,c,b")  = l2["abab_Loovv"]("X,i,j,a,b") * t2["abab_vvoo"]("a,c,i,j");

        // sigmal1_aa += +0.50 <m,b||e,c>_abab l2_abab(i,j,a,b) t2_abab(a,c,i,j) @ +0.50 <m,b||e,c>_abab l2_abab(j,i,a,b) t2_abab(a,c,j,i)  // flops: o1v1L1 += o1v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= tmps_["26_bb_Lvv"]("X,c,b") * eri["baab_vovv"]("b,m,e,c");

        // sigmal1_bb += -0.50 <m,b||c,e>_bbbb l2_abab(i,j,a,b) t2_abab(a,c,i,j) @ -0.50 <m,b||c,e>_bbbb l2_abab(j,i,a,b) t2_abab(a,c,j,i)  // flops: o1v1L1 += o1v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += tmps_["26_bb_Lvv"]("X,c,b") * eri["bbbb_vovv"]("b,m,c,e");

        // sigmal2_abab += -0.50 <m,n||e,b>_abab l2_abab(i,j,a,f) t2_abab(a,b,i,j) @ -0.50 <m,n||e,b>_abab l2_abab(j,i,a,f) t2_abab(a,b,j,i)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= tmps_["26_bb_Lvv"]("X,b,f") * eri["abab_oovv"]("m,n,e,b");

        // tmps_[101_bbbb_Lvvoo](X,f,e,n,m) = 1.00 l2[abab_Loovv](X,i,j,a,e) * t2[abab_vvoo](a,b,i,j) * eri[bbbb_oovv](m,n,b,f) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["101_bbbb_Lvvoo"]("X,f,e,n,m")  = tmps_["26_bb_Lvv"]("X,b,e") * eri["bbbb_oovv"]("m,n,b,f");
        tmps_["26_bb_Lvv"].~TArrayD();
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["101_bbbb_Lvvoo"]("X,f,e,n,m");

        // sigmal2_bbbb += +0.50 P(e,f) <m,n||b,e>_bbbb l2_abab(i,j,a,f) t2_abab(a,b,i,j) @ +0.50 P(e,f) <m,n||b,e>_bbbb l2_abab(j,i,a,f) t2_abab(a,b,j,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") += tmps_["101_bbbb_Lvvoo"]("X,e,f,n,m");
        tmps_["101_bbbb_Lvvoo"].~TArrayD();

        // tmps_[27_bbbb_Lvvoo](X,f,e,m,n) = 1.00 r1[bb_Lvo](X,b,n) * eri[baab_vovv](f,i,a,b) * t2[abab_vvoo](a,e,i,m) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["27_bbbb_Lvvoo"]("X,f,e,m,n")  = r1["bb_Lvo"]("X,b,n") * reuse_tmps_["22_bbbb_vvvo"]("e,b,f,m");

        // sigmar2_bbbb += +1.00 P(m,n) P(e,f) <i,e||a,b>_abab t2_abab(a,f,i,n) r1_bb(b,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["27_bbbb_Lvvoo"]("X,e,f,n,m");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["27_bbbb_Lvvoo"]("X,e,f,m,n");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["27_bbbb_Lvvoo"]("X,f,e,m,n");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["27_bbbb_Lvvoo"]("X,f,e,n,m");
        tmps_["27_bbbb_Lvvoo"].~TArrayD();

        // tmps_[28_bbbb_Lvvoo](X,f,e,n,m) = 1.00 r1[bb_Lvo](X,b,m) * eri[bbbb_vovv](f,i,a,b) * t2[bbbb_vvoo](a,e,n,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["28_bbbb_Lvvoo"]("X,f,e,n,m")  = r1["bb_Lvo"]("X,b,m") * reuse_tmps_["24_bbbb_vvvo"]("e,b,f,n");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["28_bbbb_Lvvoo"]("X,f,e,m,n");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["28_bbbb_Lvvoo"]("X,e,f,m,n");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["28_bbbb_Lvvoo"]("X,f,e,n,m");

        // sigmar2_bbbb += -1.00 P(m,n) P(e,f) <i,e||a,b>_bbbb t2_bbbb(a,f,n,i) r1_bb(b,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += tmps_["28_bbbb_Lvvoo"]("X,e,f,n,m");
        tmps_["28_bbbb_Lvvoo"].~TArrayD();

        // tmps_[29_aaaa_Lvvoo](X,f,e,n,m) = 1.00 r1[aa_Lvo](X,b,m) * eri[aaaa_vovv](f,i,a,b) * t2[aaaa_vvoo](a,e,n,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["29_aaaa_Lvvoo"]("X,f,e,n,m")  = r1["aa_Lvo"]("X,b,m") * reuse_tmps_["25_aaaa_vvvo"]("e,b,f,n");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["29_aaaa_Lvvoo"]("X,f,e,n,m");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["29_aaaa_Lvvoo"]("X,f,e,m,n");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["29_aaaa_Lvvoo"]("X,e,f,m,n");

        // sigmar2_aaaa += -1.00 P(m,n) P(e,f) <i,e||a,b>_aaaa t2_aaaa(a,f,n,i) r1_aa(b,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += tmps_["29_aaaa_Lvvoo"]("X,e,f,n,m");
        tmps_["29_aaaa_Lvvoo"].~TArrayD();

        // tmps_[30_aaaa_Lvvoo](X,f,e,n,m) = 1.00 r1[aa_Lvo](X,b,m) * eri[abab_vovv](f,i,b,a) * t2[abab_vvoo](e,a,n,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["30_aaaa_Lvvoo"]("X,f,e,n,m")  = r1["aa_Lvo"]("X,b,m") * reuse_tmps_["26_aaaa_vvvo"]("e,b,f,n");

        // sigmar2_aaaa += +1.00 P(m,n) P(e,f) <e,i||b,a>_abab t2_abab(f,a,n,i) r1_aa(b,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += tmps_["30_aaaa_Lvvoo"]("X,e,f,n,m");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["30_aaaa_Lvvoo"]("X,f,e,n,m");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["30_aaaa_Lvvoo"]("X,e,f,m,n");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["30_aaaa_Lvvoo"]("X,f,e,m,n");
        tmps_["30_aaaa_Lvvoo"].~TArrayD();

        // tmps_[31_bb_Lvv](X,b,f) = 0.50 l2[bbbb_Loovv](X,i,j,a,f) * t2[bbbb_vvoo](b,a,i,j) // flops: o0v2L1 = o2v3L1 | mem: o0v2L1 = o0v2L1
        tmps_["31_bb_Lvv"]("X,b,f")  = 0.50 * l2["bbbb_Loovv"]("X,i,j,a,f") * t2["bbbb_vvoo"]("b,a,i,j");

        // sigmal2_abab += +0.50 <m,n||e,b>_abab l2_bbbb(i,j,a,f) t2_bbbb(b,a,i,j)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += tmps_["31_bb_Lvv"]("X,b,f") * eri["abab_oovv"]("m,n,e,b");

        // tmps_[102_bbbb_Lvvoo](X,e,f,n,m) = 1.00 eri[bbbb_oovv](m,n,b,f) * l2[bbbb_Loovv](X,i,j,a,e) * t2[bbbb_vvoo](b,a,i,j) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["102_bbbb_Lvvoo"]("X,e,f,n,m")  = eri["bbbb_oovv"]("m,n,b,f") * tmps_["31_bb_Lvv"]("X,b,e");
        tmps_["31_bb_Lvv"].~TArrayD();

        // sigmal2_bbbb += -0.50 P(e,f) <m,n||b,e>_bbbb l2_bbbb(i,j,a,f) t2_bbbb(b,a,i,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["102_bbbb_Lvvoo"]("X,f,e,n,m");
        sigmal2_bbbb("X,e,f,m,n") += tmps_["102_bbbb_Lvvoo"]("X,e,f,n,m");
        tmps_["102_bbbb_Lvvoo"].~TArrayD();

        // tmps_[32_aa_Lvv](X,b,f) = 0.50 l2[aaaa_Loovv](X,i,j,a,f) * t2[aaaa_vvoo](b,a,i,j) // flops: o0v2L1 = o2v3L1 | mem: o0v2L1 = o0v2L1
        tmps_["32_aa_Lvv"]("X,b,f")  = 0.50 * l2["aaaa_Loovv"]("X,i,j,a,f") * t2["aaaa_vvoo"]("b,a,i,j");

        // sigmal2_abab += +0.50 <m,n||b,f>_abab l2_aaaa(i,j,a,e) t2_aaaa(b,a,i,j)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += tmps_["32_aa_Lvv"]("X,b,e") * eri["abab_oovv"]("m,n,b,f");

        // tmps_[103_aaaa_Lvvoo](X,f,e,n,m) = 1.00 eri[aaaa_oovv](m,n,b,e) * l2[aaaa_Loovv](X,i,j,a,f) * t2[aaaa_vvoo](b,a,i,j) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["103_aaaa_Lvvoo"]("X,f,e,n,m")  = eri["aaaa_oovv"]("m,n,b,e") * tmps_["32_aa_Lvv"]("X,b,f");
        tmps_["32_aa_Lvv"].~TArrayD();

        // sigmal2_aaaa += -0.50 P(e,f) <m,n||b,e>_aaaa l2_aaaa(i,j,a,f) t2_aaaa(b,a,i,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["103_aaaa_Lvvoo"]("X,f,e,n,m");
        sigmal2_aaaa("X,e,f,m,n") += tmps_["103_aaaa_Lvvoo"]("X,e,f,n,m");
        tmps_["103_aaaa_Lvvoo"].~TArrayD();

        // tmps_[33_bb_Lvv](X,f,a) = 1.00 eri[abab_oovv](j,i,b,a) * r2[abab_Lvvoo](X,b,f,j,i) // flops: o0v2L1 = o2v3L1 | mem: o0v2L1 = o0v2L1
        tmps_["33_bb_Lvv"]("X,f,a")  = eri["abab_oovv"]("j,i,b,a") * r2["abab_Lvvoo"]("X,b,f,j,i");

        // sigmar2_abab += -0.50 <j,i||b,a>_abab t2_abab(e,a,m,n) r2_abab(b,f,j,i) @ -0.50 <i,j||b,a>_abab t2_abab(e,a,m,n) r2_abab(b,f,i,j)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= tmps_["33_bb_Lvv"]("X,f,a") * t2["abab_vvoo"]("e,a,m,n");

        // tmps_[107_bbbb_Lvvoo](X,e,f,n,m) = 1.00 eri[abab_oovv](j,i,b,a) * r2[abab_Lvvoo](X,b,f,j,i) * t2[bbbb_vvoo](a,e,m,n) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["107_bbbb_Lvvoo"]("X,e,f,n,m")  = tmps_["33_bb_Lvv"]("X,f,a") * t2["bbbb_vvoo"]("a,e,m,n");
        tmps_["33_bb_Lvv"].~TArrayD();
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["107_bbbb_Lvvoo"]("X,f,e,n,m");

        // sigmar2_bbbb += +0.50 P(e,f) <j,i||b,a>_abab t2_bbbb(a,e,m,n) r2_abab(b,f,j,i) @ +0.50 P(e,f) <i,j||b,a>_abab t2_bbbb(a,e,m,n) r2_abab(b,f,i,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += tmps_["107_bbbb_Lvvoo"]("X,e,f,n,m");
        tmps_["107_bbbb_Lvvoo"].~TArrayD();

        // tmps_[34_bb_Lvv](X,f,a) = 0.50 eri[bbbb_oovv](j,i,a,b) * r2[bbbb_Lvvoo](X,b,f,j,i) // flops: o0v2L1 = o2v3L1 | mem: o0v2L1 = o0v2L1
        tmps_["34_bb_Lvv"]("X,f,a")  = 0.50 * eri["bbbb_oovv"]("j,i,a,b") * r2["bbbb_Lvvoo"]("X,b,f,j,i");

        // sigmar2_abab += +0.50 <j,i||a,b>_bbbb t2_abab(e,a,m,n) r2_bbbb(b,f,j,i)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += tmps_["34_bb_Lvv"]("X,f,a") * t2["abab_vvoo"]("e,a,m,n");

        // tmps_[106_bbbb_Lvvoo](X,e,f,n,m) = 1.00 t2[bbbb_vvoo](a,f,m,n) * eri[bbbb_oovv](j,i,a,b) * r2[bbbb_Lvvoo](X,b,e,j,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["106_bbbb_Lvvoo"]("X,e,f,n,m")  = t2["bbbb_vvoo"]("a,f,m,n") * tmps_["34_bb_Lvv"]("X,e,a");
        tmps_["34_bb_Lvv"].~TArrayD();
        sigmar2_bbbb("X,e,f,m,n") += tmps_["106_bbbb_Lvvoo"]("X,e,f,n,m");

        // sigmar2_bbbb += -0.50 P(e,f) <j,i||a,b>_bbbb t2_bbbb(a,e,m,n) r2_bbbb(b,f,j,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["106_bbbb_Lvvoo"]("X,f,e,n,m");
        tmps_["106_bbbb_Lvvoo"].~TArrayD();

        // tmps_[35_aa_Lvv](X,e,a) = 1.00 eri[abab_oovv](j,i,a,b) * r2[abab_Lvvoo](X,e,b,j,i) // flops: o0v2L1 = o2v3L1 | mem: o0v2L1 = o0v2L1
        tmps_["35_aa_Lvv"]("X,e,a")  = eri["abab_oovv"]("j,i,a,b") * r2["abab_Lvvoo"]("X,e,b,j,i");

        // sigmar2_abab += -0.50 <j,i||a,b>_abab t2_abab(a,f,m,n) r2_abab(e,b,j,i) @ -0.50 <i,j||a,b>_abab t2_abab(a,f,m,n) r2_abab(e,b,i,j)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= tmps_["35_aa_Lvv"]("X,e,a") * t2["abab_vvoo"]("a,f,m,n");

        // tmps_[109_aaaa_Lvvoo](X,e,f,n,m) = 1.00 t2[aaaa_vvoo](a,f,m,n) * eri[abab_oovv](j,i,a,b) * r2[abab_Lvvoo](X,e,b,j,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["109_aaaa_Lvvoo"]("X,e,f,n,m")  = t2["aaaa_vvoo"]("a,f,m,n") * tmps_["35_aa_Lvv"]("X,e,a");
        tmps_["35_aa_Lvv"].~TArrayD();

        // sigmar2_aaaa += +0.50 P(e,f) <j,i||a,b>_abab t2_aaaa(a,e,m,n) r2_abab(f,b,j,i) @ +0.50 P(e,f) <i,j||a,b>_abab t2_aaaa(a,e,m,n) r2_abab(f,b,i,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += tmps_["109_aaaa_Lvvoo"]("X,f,e,n,m");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["109_aaaa_Lvvoo"]("X,e,f,n,m");
        tmps_["109_aaaa_Lvvoo"].~TArrayD();

        // tmps_[36_aa_Lvv](X,f,a) = 0.50 eri[aaaa_oovv](j,i,a,b) * r2[aaaa_Lvvoo](X,b,f,j,i) // flops: o0v2L1 = o2v3L1 | mem: o0v2L1 = o0v2L1
        tmps_["36_aa_Lvv"]("X,f,a")  = 0.50 * eri["aaaa_oovv"]("j,i,a,b") * r2["aaaa_Lvvoo"]("X,b,f,j,i");

        // sigmar2_abab += +0.50 <j,i||a,b>_aaaa t2_abab(a,f,m,n) r2_aaaa(b,e,j,i)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += tmps_["36_aa_Lvv"]("X,e,a") * t2["abab_vvoo"]("a,f,m,n");

        // tmps_[110_aaaa_Lvvoo](X,f,e,n,m) = 1.00 t2[aaaa_vvoo](a,e,m,n) * eri[aaaa_oovv](j,i,a,b) * r2[aaaa_Lvvoo](X,b,f,j,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["110_aaaa_Lvvoo"]("X,f,e,n,m")  = t2["aaaa_vvoo"]("a,e,m,n") * tmps_["36_aa_Lvv"]("X,f,a");
        tmps_["36_aa_Lvv"].~TArrayD();

        // sigmar2_aaaa += -0.50 P(e,f) <j,i||a,b>_aaaa t2_aaaa(a,e,m,n) r2_aaaa(b,f,j,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["110_aaaa_Lvvoo"]("X,f,e,n,m");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["110_aaaa_Lvvoo"]("X,e,f,n,m");
        tmps_["110_aaaa_Lvvoo"].~TArrayD();

        // tmps_[37_bb_Lvv](X,c,b) = 0.50 l2[bbbb_Loovv](X,i,j,b,a) * t2[bbbb_vvoo](c,a,i,j) // flops: o0v2L1 = o2v3L1 | mem: o0v2L1 = o0v2L1
        tmps_["37_bb_Lvv"]("X,c,b")  = 0.50 * l2["bbbb_Loovv"]("X,i,j,b,a") * t2["bbbb_vvoo"]("c,a,i,j");

        // sigmal1_aa += +0.50 <m,b||e,c>_abab l2_bbbb(i,j,b,a) t2_bbbb(c,a,i,j)  // flops: o1v1L1 += o1v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= tmps_["37_bb_Lvv"]("X,c,b") * eri["baab_vovv"]("b,m,e,c");

        // sigmal1_bb += -0.50 <m,b||c,e>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(c,a,i,j)  // flops: o1v1L1 += o1v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += tmps_["37_bb_Lvv"]("X,c,b") * eri["bbbb_vovv"]("b,m,c,e");
        tmps_["37_bb_Lvv"].~TArrayD();

        // tmps_[38_aa_Lvv](X,c,b) = 0.50 l2[aaaa_Loovv](X,i,j,b,a) * t2[aaaa_vvoo](c,a,i,j) // flops: o0v2L1 = o2v3L1 | mem: o0v2L1 = o0v2L1
        tmps_["38_aa_Lvv"]("X,c,b")  = 0.50 * l2["aaaa_Loovv"]("X,i,j,b,a") * t2["aaaa_vvoo"]("c,a,i,j");

        // sigmal1_aa += -0.50 <m,b||c,e>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(c,a,i,j)  // flops: o1v1L1 += o1v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += tmps_["38_aa_Lvv"]("X,c,b") * eri["aaaa_vovv"]("b,m,c,e");

        // sigmal1_bb += +0.50 <b,m||c,e>_abab l2_aaaa(i,j,b,a) t2_aaaa(c,a,i,j)  // flops: o1v1L1 += o1v3L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += tmps_["38_aa_Lvv"]("X,c,b") * eri["abab_vovv"]("b,m,c,e");
        tmps_["38_aa_Lvv"].~TArrayD();

        // tmps_[39_bbbb_Lvvoo](X,f,e,m,n) = 1.00 eri[bbbb_vovv](a,n,e,f) * l1[bb_Lov](X,m,a) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["39_bbbb_Lvvoo"]("X,f,e,m,n")  = eri["bbbb_vovv"]("a,n,e,f") * l1["bb_Lov"]("X,m,a");

        // sigmal2_bbbb += -1.00 P(m,n) <n,a||e,f>_bbbb l1_bb(m,a)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") += tmps_["39_bbbb_Lvvoo"]("X,f,e,m,n");
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["39_bbbb_Lvvoo"]("X,f,e,n,m");
        tmps_["39_bbbb_Lvvoo"].~TArrayD();

        // tmps_[40_bbbb_Lvvoo](X,f,e,n,m) = 1.00 l2[bbbb_Loovv](X,m,n,a,e) * eri[bbbb_oovv](j,i,b,f) * t2[bbbb_vvoo](b,a,j,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["40_bbbb_Lvvoo"]("X,f,e,n,m")  = l2["bbbb_Loovv"]("X,m,n,a,e") * reuse_tmps_["16_bb_vv"]("a,f");
        sigmal2_bbbb("X,e,f,m,n") += tmps_["40_bbbb_Lvvoo"]("X,f,e,n,m");

        // sigmal2_bbbb += -0.50 P(e,f) <j,i||b,e>_bbbb l2_bbbb(m,n,a,f) t2_bbbb(b,a,j,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["40_bbbb_Lvvoo"]("X,e,f,n,m");
        tmps_["40_bbbb_Lvvoo"].~TArrayD();

        // tmps_[41_bbbb_Lvvoo](X,f,e,n,m) = 1.00 l2[bbbb_Loovv](X,m,n,a,e) * eri[abab_oovv](j,i,b,f) * t2[abab_vvoo](b,a,j,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["41_bbbb_Lvvoo"]("X,f,e,n,m")  = l2["bbbb_Loovv"]("X,m,n,a,e") * reuse_tmps_["15_bb_vv"]("a,f");
        sigmal2_bbbb("X,e,f,m,n") += tmps_["41_bbbb_Lvvoo"]("X,f,e,n,m");

        // sigmal2_bbbb += -0.50 P(e,f) <j,i||b,e>_abab l2_bbbb(m,n,a,f) t2_abab(b,a,j,i) @ -0.50 P(e,f) <i,j||b,e>_abab l2_bbbb(m,n,a,f) t2_abab(b,a,i,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["41_bbbb_Lvvoo"]("X,e,f,n,m");
        tmps_["41_bbbb_Lvvoo"].~TArrayD();

        // tmps_[42_aaaa_Lvvoo](X,e,f,n,m) = 1.00 l2[aaaa_Loovv](X,m,n,a,f) * eri[abab_oovv](j,i,e,b) * t2[abab_vvoo](a,b,j,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["42_aaaa_Lvvoo"]("X,e,f,n,m")  = l2["aaaa_Loovv"]("X,m,n,a,f") * reuse_tmps_["14_aa_vv"]("a,e");
        sigmal2_aaaa("X,e,f,m,n") += tmps_["42_aaaa_Lvvoo"]("X,f,e,n,m");

        // sigmal2_aaaa += -0.50 P(e,f) <j,i||e,b>_abab l2_aaaa(m,n,a,f) t2_abab(a,b,j,i) @ -0.50 P(e,f) <i,j||e,b>_abab l2_aaaa(m,n,a,f) t2_abab(a,b,i,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["42_aaaa_Lvvoo"]("X,e,f,n,m");
        tmps_["42_aaaa_Lvvoo"].~TArrayD();

        // tmps_[43_aaaa_Lvvoo](X,e,f,n,m) = 1.00 f[aa_vv](a,f) * l2[aaaa_Loovv](X,m,n,a,e) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["43_aaaa_Lvvoo"]("X,e,f,n,m")  = f["aa_vv"]("a,f") * l2["aaaa_Loovv"]("X,m,n,a,e");

        // sigmal2_aaaa += +1.00 P(e,f) f_aa(a,e) l2_aaaa(m,n,a,f)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") += tmps_["43_aaaa_Lvvoo"]("X,f,e,n,m");
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["43_aaaa_Lvvoo"]("X,e,f,n,m");
        tmps_["43_aaaa_Lvvoo"].~TArrayD();

        // tmps_[44_aaaa_Lvvoo](X,e,f,n,m) = 1.00 l2[aaaa_Loovv](X,m,n,a,f) * eri[aaaa_oovv](j,i,b,e) * t2[aaaa_vvoo](b,a,j,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["44_aaaa_Lvvoo"]("X,e,f,n,m")  = l2["aaaa_Loovv"]("X,m,n,a,f") * reuse_tmps_["17_aa_vv"]("a,e");
        sigmal2_aaaa("X,e,f,m,n") += tmps_["44_aaaa_Lvvoo"]("X,f,e,n,m");

        // sigmal2_aaaa += -0.50 P(e,f) <j,i||b,e>_aaaa l2_aaaa(m,n,a,f) t2_aaaa(b,a,j,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["44_aaaa_Lvvoo"]("X,e,f,n,m");
        tmps_["44_aaaa_Lvvoo"].~TArrayD();

        // tmps_[45_aaaa_Lvvoo](X,f,e,m,n) = 1.00 eri[aaaa_vovv](a,n,e,f) * l1[aa_Lov](X,m,a) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["45_aaaa_Lvvoo"]("X,f,e,m,n")  = eri["aaaa_vovv"]("a,n,e,f") * l1["aa_Lov"]("X,m,a");
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["45_aaaa_Lvvoo"]("X,f,e,n,m");

        // sigmal2_aaaa += -1.00 P(m,n) <n,a||e,f>_aaaa l1_aa(m,a)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") += tmps_["45_aaaa_Lvvoo"]("X,f,e,m,n");
        tmps_["45_aaaa_Lvvoo"].~TArrayD();

        // tmps_[46_bbbb_Lvvoo](X,f,e,m,n) = 1.00 eri[bbbb_vvvo](e,f,a,n) * r1[bb_Lvo](X,a,m) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["46_bbbb_Lvvoo"]("X,f,e,m,n")  = eri["bbbb_vvvo"]("e,f,a,n") * r1["bb_Lvo"]("X,a,m");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["46_bbbb_Lvvoo"]("X,f,e,n,m");

        // sigmar2_bbbb += +1.00 P(m,n) <e,f||a,n>_bbbb r1_bb(a,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += tmps_["46_bbbb_Lvvoo"]("X,f,e,m,n");
        tmps_["46_bbbb_Lvvoo"].~TArrayD();

        // tmps_[47_bbbb_Lvvoo](X,e,f,n,m) = 1.00 f[bb_vv](f,a) * r2[bbbb_Lvvoo](X,a,e,m,n) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["47_bbbb_Lvvoo"]("X,e,f,n,m")  = f["bb_vv"]("f,a") * r2["bbbb_Lvvoo"]("X,a,e,m,n");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["47_bbbb_Lvvoo"]("X,e,f,n,m");

        // sigmar2_bbbb += +1.00 P(e,f) f_bb(e,a) r2_bbbb(a,f,m,n)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += tmps_["47_bbbb_Lvvoo"]("X,f,e,n,m");
        tmps_["47_bbbb_Lvvoo"].~TArrayD();

        // tmps_[48_bbbb_Lvvoo](X,e,f,n,m) = 1.00 r2[bbbb_Lvvoo](X,b,f,m,n) * eri[abab_oovv](j,i,a,b) * t2[abab_vvoo](a,e,j,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["48_bbbb_Lvvoo"]("X,e,f,n,m")  = r2["bbbb_Lvvoo"]("X,b,f,m,n") * reuse_tmps_["15_bb_vv"]("e,b");

        // sigmar2_bbbb += -0.50 P(e,f) <j,i||a,b>_abab t2_abab(a,e,j,i) r2_bbbb(b,f,m,n) @ -0.50 P(e,f) <i,j||a,b>_abab t2_abab(a,e,i,j) r2_bbbb(b,f,m,n)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["48_bbbb_Lvvoo"]("X,e,f,n,m");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["48_bbbb_Lvvoo"]("X,f,e,n,m");
        tmps_["48_bbbb_Lvvoo"].~TArrayD();

        // tmps_[49_bbbb_Lvvoo](X,f,e,n,m) = 1.00 r2[bbbb_Lvvoo](X,b,e,m,n) * eri[bbbb_oovv](j,i,a,b) * t2[bbbb_vvoo](a,f,j,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["49_bbbb_Lvvoo"]("X,f,e,n,m")  = r2["bbbb_Lvvoo"]("X,b,e,m,n") * reuse_tmps_["16_bb_vv"]("f,b");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["49_bbbb_Lvvoo"]("X,f,e,n,m");

        // sigmar2_bbbb += -0.50 P(e,f) <j,i||a,b>_bbbb t2_bbbb(a,e,j,i) r2_bbbb(b,f,m,n)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["49_bbbb_Lvvoo"]("X,e,f,n,m");
        tmps_["49_bbbb_Lvvoo"].~TArrayD();

        // tmps_[50_bbbb_Lvvoo](X,e,f,n,m) = 1.00 r1[bb_Lvo](X,a,m) * eri[bbbb_oovo](j,i,a,n) * t2[bbbb_vvoo](e,f,j,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["50_bbbb_Lvvoo"]("X,e,f,n,m")  = r1["bb_Lvo"]("X,a,m") * reuse_tmps_["30_bbbb_vvvo"]("f,e,a,n");

        // sigmar2_bbbb += +0.50 P(m,n) <j,i||a,n>_bbbb t2_bbbb(e,f,j,i) r1_bb(a,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += tmps_["50_bbbb_Lvvoo"]("X,e,f,n,m");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["50_bbbb_Lvvoo"]("X,e,f,m,n");
        tmps_["50_bbbb_Lvvoo"].~TArrayD();

        // tmps_[51_bbbb_Lvvoo](X,e,f,n,m) = 1.00 f[bb_vv](a,f) * l2[bbbb_Loovv](X,m,n,a,e) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["51_bbbb_Lvvoo"]("X,e,f,n,m")  = f["bb_vv"]("a,f") * l2["bbbb_Loovv"]("X,m,n,a,e");

        // sigmal2_bbbb += +1.00 P(e,f) f_bb(a,e) l2_bbbb(m,n,a,f)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") += tmps_["51_bbbb_Lvvoo"]("X,f,e,n,m");
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["51_bbbb_Lvvoo"]("X,e,f,n,m");
        tmps_["51_bbbb_Lvvoo"].~TArrayD();

        // tmps_[52_aaaa_Lvvoo](X,e,f,n,m) = 1.00 f[aa_vv](f,a) * r2[aaaa_Lvvoo](X,a,e,m,n) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["52_aaaa_Lvvoo"]("X,e,f,n,m")  = f["aa_vv"]("f,a") * r2["aaaa_Lvvoo"]("X,a,e,m,n");

        // sigmar2_aaaa += +1.00 P(e,f) f_aa(e,a) r2_aaaa(a,f,m,n)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += tmps_["52_aaaa_Lvvoo"]("X,f,e,n,m");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["52_aaaa_Lvvoo"]("X,e,f,n,m");
        tmps_["52_aaaa_Lvvoo"].~TArrayD();

        // tmps_[53_aaaa_Lvvoo](X,e,f,n,m) = 1.00 r1[aa_Lvo](X,a,m) * eri[aaaa_oovo](j,i,a,n) * t2[aaaa_vvoo](e,f,j,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["53_aaaa_Lvvoo"]("X,e,f,n,m")  = r1["aa_Lvo"]("X,a,m") * reuse_tmps_["31_aaaa_vvvo"]("f,e,a,n");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["53_aaaa_Lvvoo"]("X,e,f,m,n");

        // sigmar2_aaaa += +0.50 P(m,n) <j,i||a,n>_aaaa t2_aaaa(e,f,j,i) r1_aa(a,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += tmps_["53_aaaa_Lvvoo"]("X,e,f,n,m");
        tmps_["53_aaaa_Lvvoo"].~TArrayD();

        // tmps_[54_aaaa_Lvvoo](X,e,f,n,m) = 1.00 r2[aaaa_Lvvoo](X,b,f,m,n) * eri[aaaa_oovv](j,i,a,b) * t2[aaaa_vvoo](a,e,j,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["54_aaaa_Lvvoo"]("X,e,f,n,m")  = r2["aaaa_Lvvoo"]("X,b,f,m,n") * reuse_tmps_["17_aa_vv"]("e,b");

        // sigmar2_aaaa += -0.50 P(e,f) <j,i||a,b>_aaaa t2_aaaa(a,e,j,i) r2_aaaa(b,f,m,n)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["54_aaaa_Lvvoo"]("X,e,f,n,m");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["54_aaaa_Lvvoo"]("X,f,e,n,m");
        tmps_["54_aaaa_Lvvoo"].~TArrayD();

        // tmps_[55_aaaa_Lvvoo](X,e,f,n,m) = 1.00 r2[aaaa_Lvvoo](X,b,f,m,n) * eri[abab_oovv](j,i,b,a) * t2[abab_vvoo](e,a,j,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["55_aaaa_Lvvoo"]("X,e,f,n,m")  = r2["aaaa_Lvvoo"]("X,b,f,m,n") * reuse_tmps_["14_aa_vv"]("e,b");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["55_aaaa_Lvvoo"]("X,f,e,n,m");

        // sigmar2_aaaa += -0.50 P(e,f) <j,i||b,a>_abab t2_abab(e,a,j,i) r2_aaaa(b,f,m,n) @ -0.50 P(e,f) <i,j||b,a>_abab t2_abab(e,a,i,j) r2_aaaa(b,f,m,n)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["55_aaaa_Lvvoo"]("X,e,f,n,m");
        tmps_["55_aaaa_Lvvoo"].~TArrayD();

        // tmps_[56_aaaa_Lvvoo](X,f,e,n,m) = 1.00 eri[aaaa_vvvo](e,f,a,m) * r1[aa_Lvo](X,a,n) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["56_aaaa_Lvvoo"]("X,f,e,n,m")  = eri["aaaa_vvvo"]("e,f,a,m") * r1["aa_Lvo"]("X,a,n");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["56_aaaa_Lvvoo"]("X,f,e,n,m");

        // sigmar2_aaaa += +1.00 P(m,n) <e,f||a,n>_aaaa r1_aa(a,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += tmps_["56_aaaa_Lvvoo"]("X,f,e,m,n");
        tmps_["56_aaaa_Lvvoo"].~TArrayD();

        // tmps_[57_aa_Loo](X,j,m) = 1.00 l2[abab_Loovv](X,m,i,b,a) * t2[abab_vvoo](b,a,j,i) // flops: o2v0L1 = o3v2L1 | mem: o2v0L1 = o2v0L1
        tmps_["57_aa_Loo"]("X,j,m")  = l2["abab_Loovv"]("X,m,i,b,a") * t2["abab_vvoo"]("b,a,j,i");

        // sigmal1_aa += -0.50 f_aa(j,e) l2_abab(m,i,b,a) t2_abab(b,a,j,i) @ -0.50 f_aa(j,e) l2_abab(m,i,a,b) t2_abab(a,b,j,i)  // flops: o1v1L1 += o2v1L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= tmps_["57_aa_Loo"]("X,j,m") * f["aa_ov"]("j,e");

        // sigmal1_aa += -0.50 <m,k||e,j>_aaaa l2_abab(j,i,b,a) t2_abab(b,a,k,i) @ -0.50 <m,k||e,j>_aaaa l2_abab(j,i,a,b) t2_abab(a,b,k,i)  // flops: o1v1L1 += o3v1L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= tmps_["57_aa_Loo"]("X,k,j") * eri["aaaa_oovo"]("m,k,e,j");

        // sigmal1_bb += -0.50 <k,m||j,e>_abab l2_abab(j,i,b,a) t2_abab(b,a,k,i) @ -0.50 <k,m||j,e>_abab l2_abab(j,i,a,b) t2_abab(a,b,k,i)  // flops: o1v1L1 += o3v1L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += tmps_["57_aa_Loo"]("X,k,j") * eri["abba_oovo"]("k,m,e,j");

        // sigmal2_abab += -0.50 <j,n||e,f>_abab l2_abab(m,i,b,a) t2_abab(b,a,j,i) @ -0.50 <j,n||e,f>_abab l2_abab(m,i,a,b) t2_abab(a,b,j,i)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= tmps_["57_aa_Loo"]("X,j,m") * eri["abab_oovv"]("j,n,e,f");

        // tmps_[115_aaaa_Lvvoo](X,f,e,n,m) = 1.00 eri[aaaa_oovv](m,j,e,f) * l2[abab_Loovv](X,n,i,b,a) * t2[abab_vvoo](b,a,j,i) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["115_aaaa_Lvvoo"]("X,f,e,n,m")  = eri["aaaa_oovv"]("m,j,e,f") * tmps_["57_aa_Loo"]("X,j,n");
        tmps_["57_aa_Loo"].~TArrayD();

        // sigmal2_aaaa += +0.50 P(m,n) <n,j||e,f>_aaaa l2_abab(m,i,b,a) t2_abab(b,a,j,i) @ +0.50 P(m,n) <n,j||e,f>_aaaa l2_abab(m,i,a,b) t2_abab(a,b,j,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") += tmps_["115_aaaa_Lvvoo"]("X,f,e,m,n");
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["115_aaaa_Lvvoo"]("X,f,e,n,m");
        tmps_["115_aaaa_Lvvoo"].~TArrayD();

        // tmps_[58_bb_Loo](X,k,j) = 1.00 l2[abab_Loovv](X,i,j,b,a) * t2[abab_vvoo](b,a,i,k) // flops: o2v0L1 = o3v2L1 | mem: o2v0L1 = o2v0L1
        tmps_["58_bb_Loo"]("X,k,j")  = l2["abab_Loovv"]("X,i,j,b,a") * t2["abab_vvoo"]("b,a,i,k");

        // sigmal1_bb += -0.50 f_bb(j,e) l2_abab(i,m,b,a) t2_abab(b,a,i,j) @ -0.50 f_bb(j,e) l2_abab(i,m,a,b) t2_abab(a,b,i,j)  // flops: o1v1L1 += o2v1L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= tmps_["58_bb_Loo"]("X,j,m") * f["bb_ov"]("j,e");

        // sigmal1_aa += -0.50 <m,k||e,j>_abab l2_abab(i,j,b,a) t2_abab(b,a,i,k) @ -0.50 <m,k||e,j>_abab l2_abab(i,j,a,b) t2_abab(a,b,i,k)  // flops: o1v1L1 += o3v1L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= tmps_["58_bb_Loo"]("X,k,j") * eri["abab_oovo"]("m,k,e,j");

        // sigmal1_bb += -0.50 <m,k||e,j>_bbbb l2_abab(i,j,b,a) t2_abab(b,a,i,k) @ -0.50 <m,k||e,j>_bbbb l2_abab(i,j,a,b) t2_abab(a,b,i,k)  // flops: o1v1L1 += o3v1L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= tmps_["58_bb_Loo"]("X,k,j") * eri["bbbb_oovo"]("m,k,e,j");

        // sigmal2_abab += -0.50 <m,j||e,f>_abab l2_abab(i,n,b,a) t2_abab(b,a,i,j) @ -0.50 <m,j||e,f>_abab l2_abab(i,n,a,b) t2_abab(a,b,i,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") -= tmps_["58_bb_Loo"]("X,j,n") * eri["abab_oovv"]("m,j,e,f");

        // tmps_[113_bbbb_Lvvoo](X,f,e,n,m) = 1.00 l2[abab_Loovv](X,i,m,b,a) * t2[abab_vvoo](b,a,i,j) * eri[bbbb_oovv](n,j,e,f) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["113_bbbb_Lvvoo"]("X,f,e,n,m")  = tmps_["58_bb_Loo"]("X,j,m") * eri["bbbb_oovv"]("n,j,e,f");
        tmps_["58_bb_Loo"].~TArrayD();

        // sigmal2_bbbb += +0.50 P(m,n) <n,j||e,f>_bbbb l2_abab(i,m,b,a) t2_abab(b,a,i,j) @ +0.50 P(m,n) <n,j||e,f>_bbbb l2_abab(i,m,a,b) t2_abab(a,b,i,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") += tmps_["113_bbbb_Lvvoo"]("X,f,e,n,m");
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["113_bbbb_Lvvoo"]("X,f,e,m,n");
        tmps_["113_bbbb_Lvvoo"].~TArrayD();

        // tmps_[59_bb_Loo](X,j,m) = 0.50 l2[bbbb_Loovv](X,m,i,b,a) * t2[bbbb_vvoo](b,a,i,j) // flops: o2v0L1 = o3v2L1 | mem: o2v0L1 = o2v0L1
        tmps_["59_bb_Loo"]("X,j,m")  = 0.50 * l2["bbbb_Loovv"]("X,m,i,b,a") * t2["bbbb_vvoo"]("b,a,i,j");

        // sigmal1_bb += +0.50 f_bb(j,e) l2_bbbb(m,i,b,a) t2_bbbb(b,a,i,j)  // flops: o1v1L1 += o2v1L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += tmps_["59_bb_Loo"]("X,j,m") * f["bb_ov"]("j,e");

        // sigmal2_abab += +0.50 <m,j||e,f>_abab l2_bbbb(n,i,b,a) t2_bbbb(b,a,i,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += tmps_["59_bb_Loo"]("X,j,n") * eri["abab_oovv"]("m,j,e,f");

        // tmps_[114_bbbb_Lvvoo](X,f,e,n,m) = 1.00 l2[bbbb_Loovv](X,m,i,b,a) * t2[bbbb_vvoo](b,a,i,j) * eri[bbbb_oovv](n,j,e,f) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["114_bbbb_Lvvoo"]("X,f,e,n,m")  = tmps_["59_bb_Loo"]("X,j,m") * eri["bbbb_oovv"]("n,j,e,f");
        tmps_["59_bb_Loo"].~TArrayD();
        sigmal2_bbbb("X,e,f,m,n") += tmps_["114_bbbb_Lvvoo"]("X,f,e,m,n");

        // sigmal2_bbbb += -0.50 P(m,n) <n,j||e,f>_bbbb l2_bbbb(m,i,b,a) t2_bbbb(b,a,i,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["114_bbbb_Lvvoo"]("X,f,e,n,m");
        tmps_["114_bbbb_Lvvoo"].~TArrayD();

        // tmps_[60_aa_Loo](X,j,m) = 0.50 l2[aaaa_Loovv](X,m,i,b,a) * t2[aaaa_vvoo](b,a,i,j) // flops: o2v0L1 = o3v2L1 | mem: o2v0L1 = o2v0L1
        tmps_["60_aa_Loo"]("X,j,m")  = 0.50 * l2["aaaa_Loovv"]("X,m,i,b,a") * t2["aaaa_vvoo"]("b,a,i,j");

        // sigmal1_aa += +0.50 f_aa(j,e) l2_aaaa(m,i,b,a) t2_aaaa(b,a,i,j)  // flops: o1v1L1 += o2v1L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") += tmps_["60_aa_Loo"]("X,j,m") * f["aa_ov"]("j,e");

        // sigmal2_abab += +0.50 <j,n||e,f>_abab l2_aaaa(m,i,b,a) t2_aaaa(b,a,i,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_abab("X,e,f,m,n") += tmps_["60_aa_Loo"]("X,j,m") * eri["abab_oovv"]("j,n,e,f");

        // tmps_[116_aaaa_Lvvoo](X,f,e,m,n) = 1.00 eri[aaaa_oovv](n,j,e,f) * l2[aaaa_Loovv](X,m,i,b,a) * t2[aaaa_vvoo](b,a,i,j) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["116_aaaa_Lvvoo"]("X,f,e,m,n")  = eri["aaaa_oovv"]("n,j,e,f") * tmps_["60_aa_Loo"]("X,j,m");
        tmps_["60_aa_Loo"].~TArrayD();
        sigmal2_aaaa("X,e,f,m,n") += tmps_["116_aaaa_Lvvoo"]("X,f,e,n,m");

        // sigmal2_aaaa += -0.50 P(m,n) <n,j||e,f>_aaaa l2_aaaa(m,i,b,a) t2_aaaa(b,a,i,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["116_aaaa_Lvvoo"]("X,f,e,m,n");
        tmps_["116_aaaa_Lvvoo"].~TArrayD();

        // tmps_[61_bbbb_Lvvoo](X,f,e,m,n) = 1.00 r1[bb_Lvo](X,e,j) * eri[abab_oovo](i,j,a,m) * t2[abab_vvoo](a,f,i,n) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["61_bbbb_Lvvoo"]("X,f,e,m,n")  = r1["bb_Lvo"]("X,e,j") * reuse_tmps_["49_bbbb_vooo"]("f,n,m,j");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["61_bbbb_Lvvoo"]("X,e,f,m,n");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["61_bbbb_Lvvoo"]("X,f,e,n,m");

        // sigmar2_bbbb += -1.00 P(m,n) P(e,f) <i,j||a,n>_abab t2_abab(a,e,i,m) r1_bb(f,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["61_bbbb_Lvvoo"]("X,e,f,n,m");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["61_bbbb_Lvvoo"]("X,f,e,m,n");
        tmps_["61_bbbb_Lvvoo"].~TArrayD();

        // tmps_[62_bbbb_Lvvoo](X,e,f,n,m) = 1.00 eri[bbbb_oovo](j,i,a,n) * t2[bbbb_vvoo](a,f,m,i) * r1[bb_Lvo](X,e,j) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["62_bbbb_Lvvoo"]("X,e,f,n,m")  = reuse_tmps_["46_bbbb_vooo"]("f,m,n,j") * r1["bb_Lvo"]("X,e,j");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["62_bbbb_Lvvoo"]("X,f,e,m,n");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["62_bbbb_Lvvoo"]("X,e,f,m,n");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["62_bbbb_Lvvoo"]("X,e,f,n,m");

        // sigmar2_bbbb += -1.00 P(m,n) P(e,f) <j,i||a,n>_bbbb t2_bbbb(a,e,m,i) r1_bb(f,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["62_bbbb_Lvvoo"]("X,f,e,n,m");
        tmps_["62_bbbb_Lvvoo"].~TArrayD();

        // tmps_[63_aaaa_Lvvoo](X,f,e,m,n) = 1.00 r1[aa_Lvo](X,e,j) * eri[aaaa_oovo](j,i,a,m) * t2[aaaa_vvoo](a,f,n,i) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["63_aaaa_Lvvoo"]("X,f,e,m,n")  = r1["aa_Lvo"]("X,e,j") * reuse_tmps_["47_aaaa_vooo"]("f,n,m,j");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["63_aaaa_Lvvoo"]("X,f,e,n,m");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["63_aaaa_Lvvoo"]("X,e,f,m,n");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["63_aaaa_Lvvoo"]("X,f,e,m,n");

        // sigmar2_aaaa += -1.00 P(m,n) P(e,f) <j,i||a,n>_aaaa t2_aaaa(a,e,m,i) r1_aa(f,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["63_aaaa_Lvvoo"]("X,e,f,n,m");
        tmps_["63_aaaa_Lvvoo"].~TArrayD();

        // tmps_[64_aaaa_Lvvoo](X,e,f,n,m) = 1.00 r1[aa_Lvo](X,f,j) * eri[abba_oovo](j,i,a,n) * t2[abab_vvoo](e,a,m,i) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["64_aaaa_Lvvoo"]("X,e,f,n,m")  = r1["aa_Lvo"]("X,f,j") * reuse_tmps_["48_aaaa_vooo"]("e,m,n,j");

        // sigmar2_aaaa += -1.00 P(m,n) P(e,f) <j,i||n,a>_abab t2_abab(e,a,m,i) r1_aa(f,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += tmps_["64_aaaa_Lvvoo"]("X,e,f,n,m");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["64_aaaa_Lvvoo"]("X,e,f,m,n");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["64_aaaa_Lvvoo"]("X,f,e,m,n");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["64_aaaa_Lvvoo"]("X,f,e,n,m");
        tmps_["64_aaaa_Lvvoo"].~TArrayD();

        // tmps_[65_bb_Loo](X,n,i) = 1.00 eri[abab_oovv](j,i,a,b) * r2[abab_Lvvoo](X,a,b,j,n) // flops: o2v0L1 = o3v2L1 | mem: o2v0L1 = o2v0L1
        tmps_["65_bb_Loo"]("X,n,i")  = eri["abab_oovv"]("j,i,a,b") * r2["abab_Lvvoo"]("X,a,b,j,n");

        // sigmar2_abab += -0.50 <j,i||a,b>_abab t2_abab(e,f,m,i) r2_abab(a,b,j,n) @ -0.50 <j,i||b,a>_abab t2_abab(e,f,m,i) r2_abab(b,a,j,n)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= tmps_["65_bb_Loo"]("X,n,i") * t2["abab_vvoo"]("e,f,m,i");

        // tmps_[121_bbbb_Lvvoo](X,f,e,n,m) = 1.00 eri[abab_oovv](j,i,a,b) * r2[abab_Lvvoo](X,a,b,j,m) * t2[bbbb_vvoo](e,f,n,i) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["121_bbbb_Lvvoo"]("X,f,e,n,m")  = tmps_["65_bb_Loo"]("X,m,i") * t2["bbbb_vvoo"]("e,f,n,i");
        tmps_["65_bb_Loo"].~TArrayD();

        // sigmar2_bbbb += +0.50 P(m,n) <j,i||a,b>_abab t2_bbbb(e,f,n,i) r2_abab(a,b,j,m) @ +0.50 P(m,n) <j,i||b,a>_abab t2_bbbb(e,f,n,i) r2_abab(b,a,j,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += tmps_["121_bbbb_Lvvoo"]("X,f,e,n,m");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["121_bbbb_Lvvoo"]("X,f,e,m,n");
        tmps_["121_bbbb_Lvvoo"].~TArrayD();

        // tmps_[66_bb_Loo](X,n,i) = 0.50 eri[bbbb_oovv](j,i,a,b) * r2[bbbb_Lvvoo](X,a,b,n,j) // flops: o2v0L1 = o3v2L1 | mem: o2v0L1 = o2v0L1
        tmps_["66_bb_Loo"]("X,n,i")  = 0.50 * eri["bbbb_oovv"]("j,i,a,b") * r2["bbbb_Lvvoo"]("X,a,b,n,j");

        // sigmar2_abab += +0.50 <j,i||a,b>_bbbb t2_abab(e,f,m,i) r2_bbbb(a,b,n,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += tmps_["66_bb_Loo"]("X,n,i") * t2["abab_vvoo"]("e,f,m,i");

        // tmps_[118_bbbb_Lvvoo](X,f,e,m,n) = 1.00 t2[bbbb_vvoo](e,f,n,i) * eri[bbbb_oovv](j,i,a,b) * r2[bbbb_Lvvoo](X,a,b,m,j) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["118_bbbb_Lvvoo"]("X,f,e,m,n")  = t2["bbbb_vvoo"]("e,f,n,i") * tmps_["66_bb_Loo"]("X,m,i");
        tmps_["66_bb_Loo"].~TArrayD();
        sigmar2_bbbb("X,e,f,m,n") += tmps_["118_bbbb_Lvvoo"]("X,f,e,n,m");

        // sigmar2_bbbb += -0.50 P(m,n) <j,i||a,b>_bbbb t2_bbbb(e,f,n,i) r2_bbbb(a,b,m,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["118_bbbb_Lvvoo"]("X,f,e,m,n");
        tmps_["118_bbbb_Lvvoo"].~TArrayD();

        // tmps_[67_aa_Loo](X,m,i) = 1.00 eri[abab_oovv](i,j,a,b) * r2[abab_Lvvoo](X,a,b,m,j) // flops: o2v0L1 = o3v2L1 | mem: o2v0L1 = o2v0L1
        tmps_["67_aa_Loo"]("X,m,i")  = eri["abab_oovv"]("i,j,a,b") * r2["abab_Lvvoo"]("X,a,b,m,j");

        // sigmar2_abab += -0.50 <i,j||a,b>_abab t2_abab(e,f,i,n) r2_abab(a,b,m,j) @ -0.50 <i,j||b,a>_abab t2_abab(e,f,i,n) r2_abab(b,a,m,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= tmps_["67_aa_Loo"]("X,m,i") * t2["abab_vvoo"]("e,f,i,n");

        // tmps_[125_aaaa_Lvvoo](X,f,e,n,m) = 1.00 t2[aaaa_vvoo](e,f,m,i) * eri[abab_oovv](i,j,a,b) * r2[abab_Lvvoo](X,a,b,n,j) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["125_aaaa_Lvvoo"]("X,f,e,n,m")  = t2["aaaa_vvoo"]("e,f,m,i") * tmps_["67_aa_Loo"]("X,n,i");
        tmps_["67_aa_Loo"].~TArrayD();

        // sigmar2_aaaa += +0.50 P(m,n) <i,j||a,b>_abab t2_aaaa(e,f,n,i) r2_abab(a,b,m,j) @ +0.50 P(m,n) <i,j||b,a>_abab t2_aaaa(e,f,n,i) r2_abab(b,a,m,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += tmps_["125_aaaa_Lvvoo"]("X,f,e,m,n");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["125_aaaa_Lvvoo"]("X,f,e,n,m");
        tmps_["125_aaaa_Lvvoo"].~TArrayD();

        // tmps_[68_aa_Loo](X,m,i) = 0.50 eri[aaaa_oovv](j,i,a,b) * r2[aaaa_Lvvoo](X,a,b,m,j) // flops: o2v0L1 = o3v2L1 | mem: o2v0L1 = o2v0L1
        tmps_["68_aa_Loo"]("X,m,i")  = 0.50 * eri["aaaa_oovv"]("j,i,a,b") * r2["aaaa_Lvvoo"]("X,a,b,m,j");

        // sigmar2_abab += +0.50 <j,i||a,b>_aaaa t2_abab(e,f,i,n) r2_aaaa(a,b,m,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += tmps_["68_aa_Loo"]("X,m,i") * t2["abab_vvoo"]("e,f,i,n");

        // tmps_[126_aaaa_Lvvoo](X,f,e,m,n) = 1.00 t2[aaaa_vvoo](e,f,n,i) * eri[aaaa_oovv](j,i,a,b) * r2[aaaa_Lvvoo](X,a,b,m,j) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["126_aaaa_Lvvoo"]("X,f,e,m,n")  = t2["aaaa_vvoo"]("e,f,n,i") * tmps_["68_aa_Loo"]("X,m,i");
        tmps_["68_aa_Loo"].~TArrayD();

        // sigmar2_aaaa += -0.50 P(m,n) <j,i||a,b>_aaaa t2_aaaa(e,f,n,i) r2_aaaa(a,b,m,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["126_aaaa_Lvvoo"]("X,f,e,m,n");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["126_aaaa_Lvvoo"]("X,f,e,n,m");
        tmps_["126_aaaa_Lvvoo"].~TArrayD();

        // tmps_[69_bb_Loo](X,k,j) = 0.50 l2[bbbb_Loovv](X,i,j,b,a) * t2[bbbb_vvoo](b,a,i,k) // flops: o2v0L1 = o3v2L1 | mem: o2v0L1 = o2v0L1
        tmps_["69_bb_Loo"]("X,k,j")  = 0.50 * l2["bbbb_Loovv"]("X,i,j,b,a") * t2["bbbb_vvoo"]("b,a,i,k");

        // sigmal1_aa += -0.50 <m,k||e,j>_abab l2_bbbb(i,j,b,a) t2_bbbb(b,a,i,k)  // flops: o1v1L1 += o3v1L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= tmps_["69_bb_Loo"]("X,k,j") * eri["abab_oovo"]("m,k,e,j");

        // sigmal1_bb += -0.50 <m,k||e,j>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(b,a,i,k)  // flops: o1v1L1 += o3v1L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") -= tmps_["69_bb_Loo"]("X,k,j") * eri["bbbb_oovo"]("m,k,e,j");
        tmps_["69_bb_Loo"].~TArrayD();

        // tmps_[70_aa_Loo](X,k,j) = 0.50 l2[aaaa_Loovv](X,i,j,b,a) * t2[aaaa_vvoo](b,a,i,k) // flops: o2v0L1 = o3v2L1 | mem: o2v0L1 = o2v0L1
        tmps_["70_aa_Loo"]("X,k,j")  = 0.50 * l2["aaaa_Loovv"]("X,i,j,b,a") * t2["aaaa_vvoo"]("b,a,i,k");

        // sigmal1_aa += -0.50 <m,k||e,j>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(b,a,i,k)  // flops: o1v1L1 += o3v1L1 | mem: o1v1L1 += o1v1L1
        sigmal1_aa("X,e,m") -= tmps_["70_aa_Loo"]("X,k,j") * eri["aaaa_oovo"]("m,k,e,j");

        // sigmal1_bb += -0.50 <k,m||j,e>_abab l2_aaaa(i,j,b,a) t2_aaaa(b,a,i,k)  // flops: o1v1L1 += o3v1L1 | mem: o1v1L1 += o1v1L1
        sigmal1_bb("X,e,m") += tmps_["70_aa_Loo"]("X,k,j") * eri["abba_oovo"]("k,m,e,j");
        tmps_["70_aa_Loo"].~TArrayD();

        // tmps_[71_bbbb_Lvvoo](X,f,e,m,n) = 1.00 eri[bbbb_oovv](n,j,a,b) * t2[bbbb_vvoo](a,b,i,j) * l2[bbbb_Loovv](X,m,i,e,f) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["71_bbbb_Lvvoo"]("X,f,e,m,n")  = reuse_tmps_["19_bb_oo"]("i,n") * l2["bbbb_Loovv"]("X,m,i,e,f");

        // sigmal2_bbbb += -0.50 P(m,n) <n,j||a,b>_bbbb l2_bbbb(m,i,e,f) t2_bbbb(a,b,i,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["71_bbbb_Lvvoo"]("X,f,e,m,n");
        sigmal2_bbbb("X,e,f,m,n") += tmps_["71_bbbb_Lvvoo"]("X,f,e,n,m");
        tmps_["71_bbbb_Lvvoo"].~TArrayD();

        // tmps_[72_bbbb_Lvvoo](X,f,e,m,n) = 1.00 l2[bbbb_Loovv](X,n,i,e,f) * eri[abab_oovv](j,m,a,b) * t2[abab_vvoo](a,b,j,i) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["72_bbbb_Lvvoo"]("X,f,e,m,n")  = l2["bbbb_Loovv"]("X,n,i,e,f") * reuse_tmps_["18_bb_oo"]("i,m");

        // sigmal2_bbbb += -0.50 P(m,n) <j,n||a,b>_abab l2_bbbb(m,i,e,f) t2_abab(a,b,j,i) @ -0.50 P(m,n) <j,n||b,a>_abab l2_bbbb(m,i,e,f) t2_abab(b,a,j,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["72_bbbb_Lvvoo"]("X,f,e,n,m");
        sigmal2_bbbb("X,e,f,m,n") += tmps_["72_bbbb_Lvvoo"]("X,f,e,m,n");
        tmps_["72_bbbb_Lvvoo"].~TArrayD();

        // tmps_[73_bbbb_Lvvoo](X,e,f,n,m) = 1.00 eri[bbbb_oovo](m,n,f,i) * l1[bb_Lov](X,i,e) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["73_bbbb_Lvvoo"]("X,e,f,n,m")  = eri["bbbb_oovo"]("m,n,f,i") * l1["bb_Lov"]("X,i,e");
        sigmal2_bbbb("X,e,f,m,n") += tmps_["73_bbbb_Lvvoo"]("X,e,f,n,m");

        // sigmal2_bbbb += -1.00 P(e,f) <m,n||e,i>_bbbb l1_bb(i,f)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["73_bbbb_Lvvoo"]("X,f,e,n,m");
        tmps_["73_bbbb_Lvvoo"].~TArrayD();

        // tmps_[74_aaaa_Lvvoo](X,f,e,m,n) = 1.00 l2[aaaa_Loovv](X,n,i,e,f) * eri[aaaa_oovv](m,j,a,b) * t2[aaaa_vvoo](a,b,i,j) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["74_aaaa_Lvvoo"]("X,f,e,m,n")  = l2["aaaa_Loovv"]("X,n,i,e,f") * reuse_tmps_["21_aa_oo"]("i,m");

        // sigmal2_aaaa += -0.50 P(m,n) <n,j||a,b>_aaaa l2_aaaa(m,i,e,f) t2_aaaa(a,b,i,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["74_aaaa_Lvvoo"]("X,f,e,n,m");
        sigmal2_aaaa("X,e,f,m,n") += tmps_["74_aaaa_Lvvoo"]("X,f,e,m,n");
        tmps_["74_aaaa_Lvvoo"].~TArrayD();

        // tmps_[75_aaaa_Lvvoo](X,f,e,n,m) = 1.00 f[aa_oo](m,i) * l2[aaaa_Loovv](X,n,i,e,f) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["75_aaaa_Lvvoo"]("X,f,e,n,m")  = f["aa_oo"]("m,i") * l2["aaaa_Loovv"]("X,n,i,e,f");

        // sigmal2_aaaa += -1.00 P(m,n) f_aa(n,i) l2_aaaa(m,i,e,f)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["75_aaaa_Lvvoo"]("X,f,e,m,n");
        sigmal2_aaaa("X,e,f,m,n") += tmps_["75_aaaa_Lvvoo"]("X,f,e,n,m");
        tmps_["75_aaaa_Lvvoo"].~TArrayD();

        // tmps_[76_aaaa_Lvvoo](X,f,e,m,n) = 1.00 l2[aaaa_Loovv](X,n,i,e,f) * eri[abab_oovv](m,j,a,b) * t2[abab_vvoo](a,b,i,j) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["76_aaaa_Lvvoo"]("X,f,e,m,n")  = l2["aaaa_Loovv"]("X,n,i,e,f") * reuse_tmps_["20_aa_oo"]("i,m");
        sigmal2_aaaa("X,e,f,m,n") += tmps_["76_aaaa_Lvvoo"]("X,f,e,m,n");

        // sigmal2_aaaa += -0.50 P(m,n) <n,j||a,b>_abab l2_aaaa(m,i,e,f) t2_abab(a,b,i,j) @ -0.50 P(m,n) <n,j||b,a>_abab l2_aaaa(m,i,e,f) t2_abab(b,a,i,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["76_aaaa_Lvvoo"]("X,f,e,n,m");
        tmps_["76_aaaa_Lvvoo"].~TArrayD();

        // tmps_[77_bbbb_Lvvoo](X,f,e,m,n) = 1.00 r2[bbbb_Lvvoo](X,e,f,n,j) * eri[bbbb_oovv](j,i,a,b) * t2[bbbb_vvoo](a,b,m,i) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["77_bbbb_Lvvoo"]("X,f,e,m,n")  = r2["bbbb_Lvvoo"]("X,e,f,n,j") * reuse_tmps_["19_bb_oo"]("m,j");

        // sigmar2_bbbb += -0.50 P(m,n) <j,i||a,b>_bbbb t2_bbbb(a,b,n,i) r2_bbbb(e,f,m,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["77_bbbb_Lvvoo"]("X,f,e,n,m");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["77_bbbb_Lvvoo"]("X,f,e,m,n");
        tmps_["77_bbbb_Lvvoo"].~TArrayD();

        // tmps_[78_bbbb_Lvvoo](X,f,e,n,m) = 1.00 f[bb_oo](i,m) * r2[bbbb_Lvvoo](X,e,f,n,i) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["78_bbbb_Lvvoo"]("X,f,e,n,m")  = f["bb_oo"]("i,m") * r2["bbbb_Lvvoo"]("X,e,f,n,i");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["78_bbbb_Lvvoo"]("X,f,e,n,m");

        // sigmar2_bbbb += -1.00 P(m,n) f_bb(i,n) r2_bbbb(e,f,m,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["78_bbbb_Lvvoo"]("X,f,e,m,n");
        tmps_["78_bbbb_Lvvoo"].~TArrayD();

        // tmps_[79_bbbb_Lvvoo](X,f,e,m,n) = 1.00 r1[bb_Lvo](X,e,i) * eri[bbbb_vovv](f,i,a,b) * t2[bbbb_vvoo](a,b,m,n) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["79_bbbb_Lvvoo"]("X,f,e,m,n")  = r1["bb_Lvo"]("X,e,i") * reuse_tmps_["32_bbbb_vooo"]("f,n,m,i");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["79_bbbb_Lvvoo"]("X,f,e,m,n");

        // sigmar2_bbbb += +0.50 P(e,f) <i,e||a,b>_bbbb t2_bbbb(a,b,m,n) r1_bb(f,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["79_bbbb_Lvvoo"]("X,e,f,m,n");
        tmps_["79_bbbb_Lvvoo"].~TArrayD();

        // tmps_[80_bbbb_Lvvoo](X,f,e,n,m) = 1.00 r2[bbbb_Lvvoo](X,e,f,m,j) * eri[abab_oovv](i,j,a,b) * t2[abab_vvoo](a,b,i,n) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["80_bbbb_Lvvoo"]("X,f,e,n,m")  = r2["bbbb_Lvvoo"]("X,e,f,m,j") * reuse_tmps_["18_bb_oo"]("n,j");

        // sigmar2_bbbb += -0.50 P(m,n) <i,j||a,b>_abab t2_abab(a,b,i,n) r2_bbbb(e,f,m,j) @ -0.50 P(m,n) <i,j||b,a>_abab t2_abab(b,a,i,n) r2_bbbb(e,f,m,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["80_bbbb_Lvvoo"]("X,f,e,n,m");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["80_bbbb_Lvvoo"]("X,f,e,m,n");
        tmps_["80_bbbb_Lvvoo"].~TArrayD();

        // tmps_[81_bbbb_Lvvoo](X,e,f,m,n) = 1.00 r1[bb_Lvo](X,f,i) * f[bb_ov](i,a) * t2[bbbb_vvoo](a,e,m,n) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["81_bbbb_Lvvoo"]("X,e,f,m,n")  = r1["bb_Lvo"]("X,f,i") * reuse_tmps_["50_bbbb_vooo"]("e,n,m,i");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["81_bbbb_Lvvoo"]("X,f,e,m,n");

        // sigmar2_bbbb += +1.00 P(e,f) f_bb(i,a) t2_bbbb(a,e,m,n) r1_bb(f,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += tmps_["81_bbbb_Lvvoo"]("X,e,f,m,n");
        tmps_["81_bbbb_Lvvoo"].~TArrayD();

        // tmps_[82_bbbb_Lvvoo](X,e,f,n,m) = 1.00 eri[bbbb_vooo](f,i,m,n) * r1[bb_Lvo](X,e,i) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["82_bbbb_Lvvoo"]("X,e,f,n,m")  = eri["bbbb_vooo"]("f,i,m,n") * r1["bb_Lvo"]("X,e,i");

        // sigmar2_bbbb += +1.00 P(e,f) <i,e||m,n>_bbbb r1_bb(f,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["82_bbbb_Lvvoo"]("X,f,e,n,m");
        sigmar2_bbbb("X,e,f,m,n") += tmps_["82_bbbb_Lvvoo"]("X,e,f,n,m");
        tmps_["82_bbbb_Lvvoo"].~TArrayD();

        // tmps_[83_aaaa_Lvvoo](X,f,e,n,m) = 1.00 r2[aaaa_Lvvoo](X,e,f,m,j) * eri[aaaa_oovv](j,i,a,b) * t2[aaaa_vvoo](a,b,n,i) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["83_aaaa_Lvvoo"]("X,f,e,n,m")  = r2["aaaa_Lvvoo"]("X,e,f,m,j") * reuse_tmps_["21_aa_oo"]("n,j");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["83_aaaa_Lvvoo"]("X,f,e,m,n");

        // sigmar2_aaaa += -0.50 P(m,n) <j,i||a,b>_aaaa t2_aaaa(a,b,n,i) r2_aaaa(e,f,m,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["83_aaaa_Lvvoo"]("X,f,e,n,m");
        tmps_["83_aaaa_Lvvoo"].~TArrayD();

        // tmps_[84_aaaa_Lvvoo](X,e,f,m,n) = 1.00 r1[aa_Lvo](X,f,i) * f[aa_ov](i,a) * t2[aaaa_vvoo](a,e,m,n) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["84_aaaa_Lvvoo"]("X,e,f,m,n")  = r1["aa_Lvo"]("X,f,i") * reuse_tmps_["51_aaaa_vooo"]("e,n,m,i");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["84_aaaa_Lvvoo"]("X,f,e,m,n");

        // sigmar2_aaaa += +1.00 P(e,f) f_aa(i,a) t2_aaaa(a,e,m,n) r1_aa(f,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += tmps_["84_aaaa_Lvvoo"]("X,e,f,m,n");
        tmps_["84_aaaa_Lvvoo"].~TArrayD();

        // tmps_[85_aaaa_Lvvoo](X,f,e,n,m) = 1.00 r2[aaaa_Lvvoo](X,e,f,m,j) * eri[abab_oovv](j,i,a,b) * t2[abab_vvoo](a,b,n,i) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["85_aaaa_Lvvoo"]("X,f,e,n,m")  = r2["aaaa_Lvvoo"]("X,e,f,m,j") * reuse_tmps_["20_aa_oo"]("n,j");

        // sigmar2_aaaa += -0.50 P(m,n) <j,i||a,b>_abab t2_abab(a,b,n,i) r2_aaaa(e,f,m,j) @ -0.50 P(m,n) <j,i||b,a>_abab t2_abab(b,a,n,i) r2_aaaa(e,f,m,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["85_aaaa_Lvvoo"]("X,f,e,n,m");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["85_aaaa_Lvvoo"]("X,f,e,m,n");
        tmps_["85_aaaa_Lvvoo"].~TArrayD();

        // tmps_[86_aaaa_Lvvoo](X,e,f,m,n) = 1.00 r1[aa_Lvo](X,f,i) * eri[aaaa_vovv](e,i,a,b) * t2[aaaa_vvoo](a,b,m,n) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["86_aaaa_Lvvoo"]("X,e,f,m,n")  = r1["aa_Lvo"]("X,f,i") * reuse_tmps_["33_aaaa_vooo"]("e,n,m,i");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["86_aaaa_Lvvoo"]("X,f,e,m,n");

        // sigmar2_aaaa += +0.50 P(e,f) <i,e||a,b>_aaaa t2_aaaa(a,b,m,n) r1_aa(f,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["86_aaaa_Lvvoo"]("X,e,f,m,n");
        tmps_["86_aaaa_Lvvoo"].~TArrayD();

        // tmps_[87_aaaa_Lvvoo](X,f,e,m,n) = 1.00 f[aa_oo](i,n) * r2[aaaa_Lvvoo](X,e,f,m,i) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["87_aaaa_Lvvoo"]("X,f,e,m,n")  = f["aa_oo"]("i,n") * r2["aaaa_Lvvoo"]("X,e,f,m,i");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["87_aaaa_Lvvoo"]("X,f,e,n,m");

        // sigmar2_aaaa += -1.00 P(m,n) f_aa(i,n) r2_aaaa(e,f,m,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["87_aaaa_Lvvoo"]("X,f,e,m,n");
        tmps_["87_aaaa_Lvvoo"].~TArrayD();

        // tmps_[88_aaaa_Lvvoo](X,e,f,n,m) = 1.00 eri[aaaa_vooo](f,i,m,n) * r1[aa_Lvo](X,e,i) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["88_aaaa_Lvvoo"]("X,e,f,n,m")  = eri["aaaa_vooo"]("f,i,m,n") * r1["aa_Lvo"]("X,e,i");
        sigmar2_aaaa("X,e,f,m,n") += tmps_["88_aaaa_Lvvoo"]("X,e,f,n,m");

        // sigmar2_aaaa += +1.00 P(e,f) <i,e||m,n>_aaaa r1_aa(f,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["88_aaaa_Lvvoo"]("X,f,e,n,m");
        tmps_["88_aaaa_Lvvoo"].~TArrayD();

        // tmps_[89_bbbb_Lvvoo](X,f,e,m,n) = 1.00 f[bb_oo](n,i) * l2[bbbb_Loovv](X,m,i,e,f) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["89_bbbb_Lvvoo"]("X,f,e,m,n")  = f["bb_oo"]("n,i") * l2["bbbb_Loovv"]("X,m,i,e,f");

        // sigmal2_bbbb += -1.00 P(m,n) f_bb(n,i) l2_bbbb(m,i,e,f)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_bbbb("X,e,f,m,n") -= tmps_["89_bbbb_Lvvoo"]("X,f,e,m,n");
        sigmal2_bbbb("X,e,f,m,n") += tmps_["89_bbbb_Lvvoo"]("X,f,e,n,m");
        tmps_["89_bbbb_Lvvoo"].~TArrayD();

        // tmps_[90_aaaa_Lvvoo](X,f,e,n,m) = 1.00 eri[aaaa_oovo](m,n,e,i) * l1[aa_Lov](X,i,f) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["90_aaaa_Lvvoo"]("X,f,e,n,m")  = eri["aaaa_oovo"]("m,n,e,i") * l1["aa_Lov"]("X,i,f");

        // sigmal2_aaaa += -1.00 P(e,f) <m,n||e,i>_aaaa l1_aa(i,f)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmal2_aaaa("X,e,f,m,n") -= tmps_["90_aaaa_Lvvoo"]("X,f,e,n,m");
        sigmal2_aaaa("X,e,f,m,n") += tmps_["90_aaaa_Lvvoo"]("X,e,f,n,m");
        tmps_["90_aaaa_Lvvoo"].~TArrayD();

        // tmps_[91_bb_Lvv](X,a,f) = 1.00 eri[baab_vovv](f,i,b,a) * r1[aa_Lvo](X,b,i) // flops: o0v2L1 = o1v3L1 | mem: o0v2L1 = o0v2L1
        tmps_["91_bb_Lvv"]("X,a,f")  = eri["baab_vovv"]("f,i,b,a") * r1["aa_Lvo"]("X,b,i");

        // sigmar2_abab += +1.00 <i,f||b,a>_abab t2_abab(e,a,m,n) r1_aa(b,i)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= tmps_["91_bb_Lvv"]("X,a,f") * t2["abab_vvoo"]("e,a,m,n");

        // tmps_[105_bbbb_Lvvoo](X,f,e,n,m) = 1.00 t2[bbbb_vvoo](a,e,m,n) * eri[baab_vovv](f,i,b,a) * r1[aa_Lvo](X,b,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["105_bbbb_Lvvoo"]("X,f,e,n,m")  = t2["bbbb_vvoo"]("a,e,m,n") * tmps_["91_bb_Lvv"]("X,a,f");
        tmps_["91_bb_Lvv"].~TArrayD();
        sigmar2_bbbb("X,e,f,m,n") += tmps_["105_bbbb_Lvvoo"]("X,f,e,n,m");

        // sigmar2_bbbb += +1.00 P(e,f) <i,e||b,a>_abab t2_bbbb(a,f,m,n) r1_aa(b,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["105_bbbb_Lvvoo"]("X,e,f,n,m");
        tmps_["105_bbbb_Lvvoo"].~TArrayD();

        // tmps_[92_bb_Lvv](X,a,f) = 1.00 eri[bbbb_vovv](f,i,a,b) * r1[bb_Lvo](X,b,i) // flops: o0v2L1 = o1v3L1 | mem: o0v2L1 = o0v2L1
        tmps_["92_bb_Lvv"]("X,a,f")  = eri["bbbb_vovv"]("f,i,a,b") * r1["bb_Lvo"]("X,b,i");

        // sigmar2_abab += -1.00 <i,f||a,b>_bbbb t2_abab(e,a,m,n) r1_bb(b,i)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += tmps_["92_bb_Lvv"]("X,a,f") * t2["abab_vvoo"]("e,a,m,n");

        // tmps_[108_bbbb_Lvvoo](X,f,e,n,m) = 1.00 t2[bbbb_vvoo](a,e,m,n) * eri[bbbb_vovv](f,i,a,b) * r1[bb_Lvo](X,b,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["108_bbbb_Lvvoo"]("X,f,e,n,m")  = t2["bbbb_vvoo"]("a,e,m,n") * tmps_["92_bb_Lvv"]("X,a,f");
        tmps_["92_bb_Lvv"].~TArrayD();

        // sigmar2_bbbb += -1.00 P(e,f) <i,e||a,b>_bbbb t2_bbbb(a,f,m,n) r1_bb(b,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += tmps_["108_bbbb_Lvvoo"]("X,e,f,n,m");
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["108_bbbb_Lvvoo"]("X,f,e,n,m");
        tmps_["108_bbbb_Lvvoo"].~TArrayD();

        // tmps_[93_aa_Lvv](X,a,f) = 1.00 eri[aaaa_vovv](f,i,a,b) * r1[aa_Lvo](X,b,i) // flops: o0v2L1 = o1v3L1 | mem: o0v2L1 = o0v2L1
        tmps_["93_aa_Lvv"]("X,a,f")  = eri["aaaa_vovv"]("f,i,a,b") * r1["aa_Lvo"]("X,b,i");

        // sigmar2_abab += -1.00 <i,e||a,b>_aaaa t2_abab(a,f,m,n) r1_aa(b,i)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += tmps_["93_aa_Lvv"]("X,a,e") * t2["abab_vvoo"]("a,f,m,n");

        // tmps_[112_aaaa_Lvvoo](X,f,e,n,m) = 1.00 t2[aaaa_vvoo](a,e,m,n) * eri[aaaa_vovv](f,i,a,b) * r1[aa_Lvo](X,b,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["112_aaaa_Lvvoo"]("X,f,e,n,m")  = t2["aaaa_vvoo"]("a,e,m,n") * tmps_["93_aa_Lvv"]("X,a,f");
        tmps_["93_aa_Lvv"].~TArrayD();
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["112_aaaa_Lvvoo"]("X,f,e,n,m");

        // sigmar2_aaaa += -1.00 P(e,f) <i,e||a,b>_aaaa t2_aaaa(a,f,m,n) r1_aa(b,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += tmps_["112_aaaa_Lvvoo"]("X,e,f,n,m");
        tmps_["112_aaaa_Lvvoo"].~TArrayD();

        // tmps_[94_aa_Lvv](X,a,f) = 1.00 eri[abab_vovv](f,i,a,b) * r1[bb_Lvo](X,b,i) // flops: o0v2L1 = o1v3L1 | mem: o0v2L1 = o0v2L1
        tmps_["94_aa_Lvv"]("X,a,f")  = eri["abab_vovv"]("f,i,a,b") * r1["bb_Lvo"]("X,b,i");

        // sigmar2_abab += +1.00 <e,i||a,b>_abab t2_abab(a,f,m,n) r1_bb(b,i)  // flops: o2v2L1 += o2v3L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += tmps_["94_aa_Lvv"]("X,a,e") * t2["abab_vvoo"]("a,f,m,n");

        // tmps_[111_aaaa_Lvvoo](X,e,f,n,m) = 1.00 t2[aaaa_vvoo](a,f,m,n) * eri[abab_vovv](e,i,a,b) * r1[bb_Lvo](X,b,i) // flops: o2v2L1 = o2v3L1 | mem: o2v2L1 = o2v2L1
        tmps_["111_aaaa_Lvvoo"]("X,e,f,n,m")  = t2["aaaa_vvoo"]("a,f,m,n") * tmps_["94_aa_Lvv"]("X,a,e");
        tmps_["94_aa_Lvv"].~TArrayD();
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["111_aaaa_Lvvoo"]("X,f,e,n,m");

        // sigmar2_aaaa += +1.00 P(e,f) <e,i||a,b>_abab t2_aaaa(a,f,m,n) r1_bb(b,i)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += tmps_["111_aaaa_Lvvoo"]("X,e,f,n,m");
        tmps_["111_aaaa_Lvvoo"].~TArrayD();

        // tmps_[95_bb_Loo](X,n,i) = 1.00 eri[abab_oovo](j,i,a,n) * r1[aa_Lvo](X,a,j) // flops: o2v0L1 = o3v1L1 | mem: o2v0L1 = o2v0L1
        tmps_["95_bb_Loo"]("X,n,i")  = eri["abab_oovo"]("j,i,a,n") * r1["aa_Lvo"]("X,a,j");

        // sigmar2_abab += -1.00 <j,i||a,n>_abab t2_abab(e,f,m,i) r1_aa(a,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= tmps_["95_bb_Loo"]("X,n,i") * t2["abab_vvoo"]("e,f,m,i");

        // tmps_[119_bbbb_Lvvoo](X,f,e,m,n) = 1.00 eri[abab_oovo](j,i,a,n) * r1[aa_Lvo](X,a,j) * t2[bbbb_vvoo](e,f,m,i) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["119_bbbb_Lvvoo"]("X,f,e,m,n")  = tmps_["95_bb_Loo"]("X,n,i") * t2["bbbb_vvoo"]("e,f,m,i");
        tmps_["95_bb_Loo"].~TArrayD();
        sigmar2_bbbb("X,e,f,m,n") += tmps_["119_bbbb_Lvvoo"]("X,f,e,n,m");

        // sigmar2_bbbb += -1.00 P(m,n) <j,i||a,n>_abab t2_bbbb(e,f,m,i) r1_aa(a,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["119_bbbb_Lvvoo"]("X,f,e,m,n");
        tmps_["119_bbbb_Lvvoo"].~TArrayD();

        // tmps_[96_bb_Loo](X,n,i) = 1.00 eri[bbbb_oovo](j,i,a,n) * r1[bb_Lvo](X,a,j) // flops: o2v0L1 = o3v1L1 | mem: o2v0L1 = o2v0L1
        tmps_["96_bb_Loo"]("X,n,i")  = eri["bbbb_oovo"]("j,i,a,n") * r1["bb_Lvo"]("X,a,j");

        // sigmar2_abab += -1.00 <j,i||a,n>_bbbb t2_abab(e,f,m,i) r1_bb(a,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= tmps_["96_bb_Loo"]("X,n,i") * t2["abab_vvoo"]("e,f,m,i");

        // tmps_[120_bbbb_Lvvoo](X,f,e,n,m) = 1.00 t2[bbbb_vvoo](e,f,m,i) * eri[bbbb_oovo](j,i,a,n) * r1[bb_Lvo](X,a,j) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["120_bbbb_Lvvoo"]("X,f,e,n,m")  = t2["bbbb_vvoo"]("e,f,m,i") * tmps_["96_bb_Loo"]("X,n,i");
        tmps_["96_bb_Loo"].~TArrayD();
        sigmar2_bbbb("X,e,f,m,n") += tmps_["120_bbbb_Lvvoo"]("X,f,e,m,n");

        // sigmar2_bbbb += -1.00 P(m,n) <j,i||a,n>_bbbb t2_bbbb(e,f,m,i) r1_bb(a,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["120_bbbb_Lvvoo"]("X,f,e,n,m");
        tmps_["120_bbbb_Lvvoo"].~TArrayD();

        // tmps_[97_aa_Loo](X,m,i) = 1.00 eri[aaaa_oovo](j,i,a,m) * r1[aa_Lvo](X,a,j) // flops: o2v0L1 = o3v1L1 | mem: o2v0L1 = o2v0L1
        tmps_["97_aa_Loo"]("X,m,i")  = eri["aaaa_oovo"]("j,i,a,m") * r1["aa_Lvo"]("X,a,j");

        // sigmar2_abab += -1.00 <j,i||a,m>_aaaa t2_abab(e,f,i,n) r1_aa(a,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= tmps_["97_aa_Loo"]("X,m,i") * t2["abab_vvoo"]("e,f,i,n");

        // tmps_[123_aaaa_Lvvoo](X,f,e,m,n) = 1.00 t2[aaaa_vvoo](e,f,n,i) * eri[aaaa_oovo](j,i,a,m) * r1[aa_Lvo](X,a,j) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["123_aaaa_Lvvoo"]("X,f,e,m,n")  = t2["aaaa_vvoo"]("e,f,n,i") * tmps_["97_aa_Loo"]("X,m,i");
        tmps_["97_aa_Loo"].~TArrayD();
        sigmar2_aaaa("X,e,f,m,n") += tmps_["123_aaaa_Lvvoo"]("X,f,e,m,n");

        // sigmar2_aaaa += -1.00 P(m,n) <j,i||a,n>_aaaa t2_aaaa(e,f,m,i) r1_aa(a,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["123_aaaa_Lvvoo"]("X,f,e,n,m");
        tmps_["123_aaaa_Lvvoo"].~TArrayD();

        // tmps_[98_aa_Loo](X,m,i) = 1.00 eri[abba_oovo](i,j,a,m) * r1[bb_Lvo](X,a,j) // flops: o2v0L1 = o3v1L1 | mem: o2v0L1 = o2v0L1
        tmps_["98_aa_Loo"]("X,m,i")  = eri["abba_oovo"]("i,j,a,m") * r1["bb_Lvo"]("X,a,j");

        // sigmar2_abab += -1.00 <i,j||m,a>_abab t2_abab(e,f,i,n) r1_bb(a,j)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") += tmps_["98_aa_Loo"]("X,m,i") * t2["abab_vvoo"]("e,f,i,n");

        // tmps_[124_aaaa_Lvvoo](X,f,e,m,n) = 1.00 t2[aaaa_vvoo](e,f,n,i) * eri[abba_oovo](i,j,a,m) * r1[bb_Lvo](X,a,j) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["124_aaaa_Lvvoo"]("X,f,e,m,n")  = t2["aaaa_vvoo"]("e,f,n,i") * tmps_["98_aa_Loo"]("X,m,i");
        tmps_["98_aa_Loo"].~TArrayD();

        // sigmar2_aaaa += -1.00 P(m,n) <i,j||n,a>_abab t2_aaaa(e,f,m,i) r1_bb(a,j)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += tmps_["124_aaaa_Lvvoo"]("X,f,e,n,m");
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["124_aaaa_Lvvoo"]("X,f,e,m,n");
        tmps_["124_aaaa_Lvvoo"].~TArrayD();

        // tmps_[99_bb_Loo](X,n,i) = 1.00 f[bb_ov](i,a) * r1[bb_Lvo](X,a,n) // flops: o2v0L1 = o2v1L1 | mem: o2v0L1 = o2v0L1
        tmps_["99_bb_Loo"]("X,n,i")  = f["bb_ov"]("i,a") * r1["bb_Lvo"]("X,a,n");

        // sigmar2_abab += -1.00 f_bb(i,a) t2_abab(e,f,m,i) r1_bb(a,n)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= tmps_["99_bb_Loo"]("X,n,i") * t2["abab_vvoo"]("e,f,m,i");

        // tmps_[117_bbbb_Lvvoo](X,f,e,n,m) = 1.00 t2[bbbb_vvoo](e,f,m,i) * f[bb_ov](i,a) * r1[bb_Lvo](X,a,n) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["117_bbbb_Lvvoo"]("X,f,e,n,m")  = t2["bbbb_vvoo"]("e,f,m,i") * tmps_["99_bb_Loo"]("X,n,i");
        tmps_["99_bb_Loo"].~TArrayD();
        sigmar2_bbbb("X,e,f,m,n") -= tmps_["117_bbbb_Lvvoo"]("X,f,e,n,m");

        // sigmar2_bbbb += +1.00 P(m,n) f_bb(i,a) t2_bbbb(e,f,n,i) r1_bb(a,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_bbbb("X,e,f,m,n") += tmps_["117_bbbb_Lvvoo"]("X,f,e,m,n");
        tmps_["117_bbbb_Lvvoo"].~TArrayD();

        // tmps_[100_aa_Loo](X,n,i) = 1.00 f[aa_ov](i,a) * r1[aa_Lvo](X,a,n) // flops: o2v0L1 = o2v1L1 | mem: o2v0L1 = o2v0L1
        tmps_["100_aa_Loo"]("X,n,i")  = f["aa_ov"]("i,a") * r1["aa_Lvo"]("X,a,n");

        // sigmar2_abab += -1.00 f_aa(i,a) t2_abab(e,f,i,n) r1_aa(a,m)  // flops: o2v2L1 += o3v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_abab("X,e,f,m,n") -= tmps_["100_aa_Loo"]("X,m,i") * t2["abab_vvoo"]("e,f,i,n");

        // tmps_[122_aaaa_Lvvoo](X,f,e,n,m) = 1.00 t2[aaaa_vvoo](e,f,m,i) * f[aa_ov](i,a) * r1[aa_Lvo](X,a,n) // flops: o2v2L1 = o3v2L1 | mem: o2v2L1 = o2v2L1
        tmps_["122_aaaa_Lvvoo"]("X,f,e,n,m")  = t2["aaaa_vvoo"]("e,f,m,i") * tmps_["100_aa_Loo"]("X,n,i");
        tmps_["100_aa_Loo"].~TArrayD();
        sigmar2_aaaa("X,e,f,m,n") -= tmps_["122_aaaa_Lvvoo"]("X,f,e,n,m");

        // sigmar2_aaaa += +1.00 P(m,n) f_aa(i,a) t2_aaaa(e,f,n,i) r1_aa(a,m)  // flops: o2v2L1 += o2v2L1 | mem: o2v2L1 += o2v2L1
        sigmar2_aaaa("X,e,f,m,n") += tmps_["122_aaaa_Lvvoo"]("X,f,e,m,n");
        tmps_["122_aaaa_Lvvoo"].~TArrayD();

    }

    double cc_energy = cc_wfn_->cc_energy_;
    sigmar0("I") += r0("I") * (cc_energy);
    sigmal0("I") += l0("I") * (cc_energy);

    sigmar1_aa("I,e,m") += r1["aa_Lvo"]("I,e,m") * (cc_energy);
    sigmar1_bb("I,e,m") += r1["bb_Lvo"]("I,e,m") * (cc_energy);

    sigmal1_aa("I,e,m") += l1["aa_Lov"]("I,m,e") * (cc_energy);
    sigmal1_bb("I,e,m") += l1["bb_Lov"]("I,m,e") * (cc_energy);

    sigmar2_aaaa("I,e,f,m,n") += r2["aaaa_Lvvoo"]("I,e,f,m,n") * (cc_energy);
    sigmar2_abab("I,e,f,m,n") += r2["abab_Lvvoo"]("I,e,f,m,n") * (cc_energy);
    sigmar2_bbbb("I,e,f,m,n") += r2["bbbb_Lvvoo"]("I,e,f,m,n") * (cc_energy);

    sigmal2_aaaa("I,e,f,m,n") += l2["aaaa_Loovv"]("I,m,n,e,f") * (cc_energy);
    sigmal2_abab("I,e,f,m,n") += l2["abab_Loovv"]("I,m,n,e,f") * (cc_energy);
    sigmal2_bbbb("I,e,f,m,n") += l2["bbbb_Loovv"]("I,m,n,e,f") * (cc_energy);

    world_.gop.fence();
}