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
        world_.gop.fence(); tempArray[0].~TArrayD();
        world_.gop.fence(); tempArray[1].~TArrayD();
        
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

    // one-body
    TA::TArrayD &t1_aa_vo =     cc_wfn_->amplitudes_["t1_aa"];
    TA::TArrayD &t1_bb_vo =     cc_wfn_->amplitudes_["t1_bb"];

    // two-body
    TA::TArrayD &t2_aaaa_vvoo = cc_wfn_->amplitudes_["t2_aaaa"];
    TA::TArrayD &t2_abab_vvoo = cc_wfn_->amplitudes_["t2_abab"];
    TA::TArrayD &t2_bbbb_vvoo = cc_wfn_->amplitudes_["t2_bbbb"];

    TArrayMap &V_blks_ = cc_wfn_->V_blks_;
    TArrayMap &F_blks_ = cc_wfn_->F_blks_;

    sigmaOps = vector<TA::TArrayD>(125);

    {

        // sigmaOps[0] += 1.000000 V_blks_["aaaa_vovv"]("e,i,a,b") t2_abab_vvoo("a,f,i,n")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[0]("e,b,f,n") = V_blks_["aaaa_vovv"]("e,i,a,b") * t2_abab_vvoo("a,f,i,n");

        // sigmaOps[1] += 1.000000 V_blks_["aaaa_vvvv"]("e,f,a,b") t2_aaaa_vvoo("a,b,m,n")
        // flops: o2v4: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[1]("e,f,m,n") = V_blks_["aaaa_vvvv"]("e,f,a,b") * t2_aaaa_vvoo("a,b,m,n");

        // sigmaOps[2] += 1.000000 t2_abab_vvoo("a,f,m,i") V_blks_["abab_vovv"]("e,i,a,b")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[2]("e,f,b,m") = t2_abab_vvoo("a,f,m,i") * V_blks_["abab_vovv"]("e,i,a,b");

        // sigmaOps[3] += 1.000000 t2_bbbb_vvoo("a,b,m,n") V_blks_["bbbb_vvvv"]("e,f,a,b")
        // flops: o2v4: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[3]("e,f,m,n") = t2_bbbb_vvoo("a,b,m,n") * V_blks_["bbbb_vvvv"]("e,f,a,b");

        // sigmaOps[4] += 1.000000 V_blks_["bbbb_vovv"]("e,i,a,b") t2_bbbb_vvoo("a,f,n,i")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[4]("e,b,f,n") = V_blks_["bbbb_vovv"]("e,i,a,b") * t2_bbbb_vvoo("a,f,n,i");

        // sigmaOps[5] += 1.000000 V_blks_["abab_vovv"]("e,i,b,a") t2_abab_vvoo("f,a,n,i")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[5]("e,b,f,n") = V_blks_["abab_vovv"]("e,i,b,a") * t2_abab_vvoo("f,a,n,i");

        // sigmaOps[6] += 1.000000 V_blks_["baab_vovv"]("f,i,a,b") t2_aaaa_vvoo("a,e,m,i")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[6]("e,f,b,m") = V_blks_["baab_vovv"]("f,i,a,b") * t2_aaaa_vvoo("a,e,m,i");

        // sigmaOps[7] += 1.000000 V_blks_["aaaa_vovv"]("e,i,a,b") t2_aaaa_vvoo("a,f,n,i")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[7]("e,b,f,n") = V_blks_["aaaa_vovv"]("e,i,a,b") * t2_aaaa_vvoo("a,f,n,i");

        // sigmaOps[8] += 1.000000 t2_abab_vvoo("a,f,i,n") V_blks_["baab_vovv"]("e,i,a,b")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[8]("f,e,b,n") = t2_abab_vvoo("a,f,i,n") * V_blks_["baab_vovv"]("e,i,a,b");

        // sigmaOps[9] += 1.000000 V_blks_["baab_vovv"]("f,i,b,a") t2_abab_vvoo("e,a,i,n")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[9]("b,e,f,n") = V_blks_["baab_vovv"]("f,i,b,a") * t2_abab_vvoo("e,a,i,n");

        // sigmaOps[10] += 1.000000 V_blks_["abab_vovv"]("e,i,b,a") t2_bbbb_vvoo("a,f,n,i")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[10]("e,b,f,n") = V_blks_["abab_vovv"]("e,i,b,a") * t2_bbbb_vvoo("a,f,n,i");

        // sigmaOps[11] += 1.000000 t2_abab_vvoo("b,a,m,n") V_blks_["abab_vvvv"]("e,f,b,a")
        // flops: o2v4: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[11]("e,f,m,n") = t2_abab_vvoo("b,a,m,n") * V_blks_["abab_vvvv"]("e,f,b,a");

        // sigmaOps[12] += 1.000000 t2_abab_vvoo("e,a,m,i") V_blks_["bbbb_vovv"]("f,i,a,b")
        // flops: o2v4: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[12]("e,f,b,m") = t2_abab_vvoo("e,a,m,i") * V_blks_["bbbb_vovv"]("f,i,a,b");

        // sigmaOps[13] += 1.000000 V_blks_["abab_oovv"]("j,i,b,a") t2_bbbb_vvoo("a,f,n,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[13]("b,f,j,n") = V_blks_["abab_oovv"]("j,i,b,a") * t2_bbbb_vvoo("a,f,n,i");

        // sigmaOps[14] += 1.000000 t2_abab_vvoo("e,a,n,i") V_blks_["abab_oovv"]("j,i,b,a")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[14]("e,b,n,j") = t2_abab_vvoo("e,a,n,i") * V_blks_["abab_oovv"]("j,i,b,a");

        // sigmaOps[15] += 1.000000 V_blks_["aaaa_oovv"]("j,i,a,b") t2_aaaa_vvoo("a,e,n,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[15]("b,e,j,n") = V_blks_["aaaa_oovv"]("j,i,a,b") * t2_aaaa_vvoo("a,e,n,i");

        // sigmaOps[16] += 1.000000 V_blks_["bbbb_oovv"]("j,i,a,b") t2_abab_vvoo("e,a,n,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[16]("e,b,n,j") = V_blks_["bbbb_oovv"]("j,i,a,b") * t2_abab_vvoo("e,a,n,i");

        // sigmaOps[17] += 1.000000 t2_bbbb_vvoo("a,f,n,i") V_blks_["bbbb_oovv"]("j,i,a,b")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[17]("f,b,n,j") = t2_bbbb_vvoo("a,f,n,i") * V_blks_["bbbb_oovv"]("j,i,a,b");

        // sigmaOps[18] += 1.000000 V_blks_["aaaa_oovv"]("j,i,a,b") t2_abab_vvoo("a,f,i,n")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[18]("b,f,j,n") = V_blks_["aaaa_oovv"]("j,i,a,b") * t2_abab_vvoo("a,f,i,n");

        // sigmaOps[19] += 1.000000 V_blks_["abab_oovv"]("i,j,b,a") t2_abab_vvoo("e,a,i,n") t2_abab_vvoo("b,f,m,j")
        // flops: o3v3: 2, o2v2: 1 | mem: o2v2: 3,
        sigmaOps[19]("e,f,m,n") = V_blks_["abab_oovv"]("i,j,b,a") * t2_abab_vvoo("e,a,i,n") * t2_abab_vvoo("b,f,m,j");

        // sigmaOps[20] += 1.000000 t2_aaaa_vvoo("a,e,n,i") V_blks_["abab_oovv"]("i,j,a,b")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[20]("e,b,n,j") = t2_aaaa_vvoo("a,e,n,i") * V_blks_["abab_oovv"]("i,j,a,b");

        // sigmaOps[21] += 1.000000 V_blks_["abab_oovv"]("i,j,a,b") t2_abab_vvoo("a,f,i,n")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[21]("b,f,j,n") = V_blks_["abab_oovv"]("i,j,a,b") * t2_abab_vvoo("a,f,i,n");

        // sigmaOps[22] += 1.000000 t2_abab_vvoo("a,b,m,n") V_blks_["baab_vovv"]("f,i,a,b")
        // flops: o3v3: 1, o3v1: 1 | mem: o3v1: 2,
        sigmaOps[22]("f,m,i,n") = t2_abab_vvoo("a,b,m,n") * V_blks_["baab_vovv"]("f,i,a,b");

        // sigmaOps[23] += 1.000000 V_blks_["aaaa_vovo"]("e,i,a,m") t2_abab_vvoo("a,f,i,n")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[23]("e,f,m,n") = V_blks_["aaaa_vovo"]("e,i,a,m") * t2_abab_vvoo("a,f,i,n");

        // sigmaOps[24] += 1.000000 V_blks_["aaaa_oovo"]("j,i,a,n") t2_aaaa_vvoo("e,f,j,i")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[24]("a,e,f,n") = V_blks_["aaaa_oovo"]("j,i,a,n") * t2_aaaa_vvoo("e,f,j,i");

        // sigmaOps[25] += 1.000000 t2_abab_vvoo("e,f,j,i") V_blks_["abab_oovo"]("j,i,a,n")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[25]("e,a,f,n") = t2_abab_vvoo("e,f,j,i") * V_blks_["abab_oovo"]("j,i,a,n");

        // sigmaOps[26] += 1.000000 V_blks_["aaaa_vovv"]("e,i,a,b") t2_aaaa_vvoo("a,b,m,n")
        // flops: o3v3: 1, o3v1: 1 | mem: o3v1: 2,
        sigmaOps[26]("e,i,m,n") = V_blks_["aaaa_vovv"]("e,i,a,b") * t2_aaaa_vvoo("a,b,m,n");

        // sigmaOps[27] += 1.000000 V_blks_["bbbb_vovv"]("e,i,a,b") t2_bbbb_vvoo("a,b,m,n")
        // flops: o3v3: 1, o3v1: 1 | mem: o3v1: 2,
        sigmaOps[27]("e,i,m,n") = V_blks_["bbbb_vovv"]("e,i,a,b") * t2_bbbb_vvoo("a,b,m,n");

        // sigmaOps[28] += 1.000000 t2_bbbb_vvoo("e,f,j,i") V_blks_["bbbb_oovo"]("j,i,a,n")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[28]("e,f,a,n") = t2_bbbb_vvoo("e,f,j,i") * V_blks_["bbbb_oovo"]("j,i,a,n");

        // sigmaOps[29] += 1.000000 t2_bbbb_vvoo("a,f,m,i") V_blks_["bbbb_vovo"]("e,i,a,n")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[29]("f,e,m,n") = t2_bbbb_vvoo("a,f,m,i") * V_blks_["bbbb_vovo"]("e,i,a,n");

        // sigmaOps[30] += 1.000000 t2_abab_vvoo("a,f,i,m") V_blks_["baab_vovo"]("e,i,a,n")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[30]("f,e,m,n") = t2_abab_vvoo("a,f,i,m") * V_blks_["baab_vovo"]("e,i,a,n");

        // sigmaOps[31] += 1.000000 V_blks_["abba_vovo"]("e,i,a,n") t2_abab_vvoo("f,a,m,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[31]("e,f,n,m") = V_blks_["abba_vovo"]("e,i,a,n") * t2_abab_vvoo("f,a,m,i");

        // sigmaOps[32] += 1.000000 V_blks_["baab_vovo"]("f,i,a,n") t2_aaaa_vvoo("a,e,m,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[32]("e,f,m,n") = V_blks_["baab_vovo"]("f,i,a,n") * t2_aaaa_vvoo("a,e,m,i");

        // sigmaOps[33] += 1.000000 t2_aaaa_vvoo("a,f,m,i") V_blks_["aaaa_vovo"]("e,i,a,n")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[33]("f,e,m,n") = t2_aaaa_vvoo("a,f,m,i") * V_blks_["aaaa_vovo"]("e,i,a,n");

        // sigmaOps[34] += 1.000000 V_blks_["abab_oovv"]("i,j,b,a") t2_abab_vvoo("e,a,i,n")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[34]("b,e,j,n") = V_blks_["abab_oovv"]("i,j,b,a") * t2_abab_vvoo("e,a,i,n");

        // sigmaOps[35] += 1.000000 V_blks_["abab_vovv"]("e,i,b,a") t2_abab_vvoo("b,a,m,n")
        // flops: o3v3: 1, o3v1: 1 | mem: o3v1: 2,
        sigmaOps[35]("e,m,i,n") = V_blks_["abab_vovv"]("e,i,b,a") * t2_abab_vvoo("b,a,m,n");

        // sigmaOps[36] += 1.000000 V_blks_["abab_oovv"]("j,i,a,b") t2_abab_vvoo("a,f,m,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[36]("b,f,j,m") = V_blks_["abab_oovv"]("j,i,a,b") * t2_abab_vvoo("a,f,m,i");

        // sigmaOps[37] += 1.000000 t2_bbbb_vvoo("a,f,n,i") V_blks_["abba_vovo"]("e,i,a,m")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[37]("e,f,m,n") = t2_bbbb_vvoo("a,f,n,i") * V_blks_["abba_vovo"]("e,i,a,m");

        // sigmaOps[38] += 1.000000 t2_abab_vvoo("e,a,m,i") V_blks_["bbbb_vovo"]("f,i,a,n")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[38]("e,f,m,n") = t2_abab_vvoo("e,a,m,i") * V_blks_["bbbb_vovo"]("f,i,a,n");

        // sigmaOps[39] += 1.000000 V_blks_["baba_vovo"]("f,i,a,m") t2_abab_vvoo("e,a,i,n")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[39]("e,f,m,n") = V_blks_["baba_vovo"]("f,i,a,m") * t2_abab_vvoo("e,a,i,n");

        // sigmaOps[40] += 1.000000 V_blks_["abba_oovo"]("j,i,a,m") t2_abab_vvoo("e,f,j,i")
        // flops: o3v3: 1, o1v3: 1 | mem: o1v3: 2,
        sigmaOps[40]("e,a,f,m") = V_blks_["abba_oovo"]("j,i,a,m") * t2_abab_vvoo("e,f,j,i");

        // sigmaOps[41] += 1.000000 V_blks_["abab_vovo"]("e,i,a,n") t2_abab_vvoo("a,f,m,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[41]("e,f,m,n") = V_blks_["abab_vovo"]("e,i,a,n") * t2_abab_vvoo("a,f,m,i");

        // sigmaOps[42] += 1.000000 t2_aaaa_vvoo("a,b,m,n") V_blks_["aaaa_oovv"]("j,i,a,b")
        // flops: o4v2: 1, o4v0: 1 | mem: o4v0: 2,
        sigmaOps[42]("m,n,j,i") = t2_aaaa_vvoo("a,b,m,n") * V_blks_["aaaa_oovv"]("j,i,a,b");

        // sigmaOps[43] += 1.000000 t2_bbbb_vvoo("a,b,m,n") V_blks_["bbbb_oovv"]("j,i,a,b")
        // flops: o4v2: 1, o4v0: 1 | mem: o4v0: 2,
        sigmaOps[43]("m,n,j,i") = t2_bbbb_vvoo("a,b,m,n") * V_blks_["bbbb_oovv"]("j,i,a,b");

        // sigmaOps[44] += 1.000000 V_blks_["abab_oovv"]("j,i,a,b") t2_abab_vvoo("a,b,m,n") t2_abab_vvoo("e,f,j,i")
        // flops: o4v2: 2, o2v2: 1 | mem: o2v2: 2, o4v0: 1,
        sigmaOps[44]("e,f,m,n") = V_blks_["abab_oovv"]("j,i,a,b") * t2_abab_vvoo("a,b,m,n") * t2_abab_vvoo("e,f,j,i");

        // sigmaOps[45] += 1.000000 t2_abab_vvoo("a,e,i,m") V_blks_["abab_oovo"]("i,j,a,n")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sigmaOps[45]("e,m,j,n") = t2_abab_vvoo("a,e,i,m") * V_blks_["abab_oovo"]("i,j,a,n");

        // sigmaOps[46] += 1.000000 V_blks_["abab_oooo"]("j,i,m,n") t2_abab_vvoo("e,f,j,i")
        // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[46]("e,f,m,n") = V_blks_["abab_oooo"]("j,i,m,n") * t2_abab_vvoo("e,f,j,i");

        // sigmaOps[47] += 1.000000 V_blks_["aaaa_oovo"]("j,i,a,n") t2_aaaa_vvoo("a,e,m,i")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sigmaOps[47]("e,j,n,m") = V_blks_["aaaa_oovo"]("j,i,a,n") * t2_aaaa_vvoo("a,e,m,i");

        // sigmaOps[48] += 1.000000 t2_bbbb_vvoo("e,f,j,i") V_blks_["bbbb_oooo"]("j,i,m,n")
        // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[48]("e,f,m,n") = t2_bbbb_vvoo("e,f,j,i") * V_blks_["bbbb_oooo"]("j,i,m,n");

        // sigmaOps[49] += 1.000000 t2_bbbb_vvoo("a,e,m,i") V_blks_["bbbb_oovo"]("j,i,a,n")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sigmaOps[49]("e,m,j,n") = t2_bbbb_vvoo("a,e,m,i") * V_blks_["bbbb_oovo"]("j,i,a,n");

        // sigmaOps[50] += 1.000000 V_blks_["abba_oovo"]("j,i,a,n") t2_abab_vvoo("e,a,m,i")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sigmaOps[50]("e,j,n,m") = V_blks_["abba_oovo"]("j,i,a,n") * t2_abab_vvoo("e,a,m,i");

        // sigmaOps[51] += 1.000000 V_blks_["aaaa_oooo"]("j,i,m,n") t2_aaaa_vvoo("e,f,j,i")
        // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[51]("e,f,m,n") = V_blks_["aaaa_oooo"]("j,i,m,n") * t2_aaaa_vvoo("e,f,j,i");

        // sigmaOps[52] += 1.000000 V_blks_["abba_oovo"]("j,i,a,m") t2_bbbb_vvoo("a,f,n,i")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sigmaOps[52]("f,j,m,n") = V_blks_["abba_oovo"]("j,i,a,m") * t2_bbbb_vvoo("a,f,n,i");

        // sigmaOps[53] += 1.000000 t2_abab_vvoo("a,f,m,i") V_blks_["abab_oovo"]("j,i,a,n")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sigmaOps[53]("f,m,j,n") = t2_abab_vvoo("a,f,m,i") * V_blks_["abab_oovo"]("j,i,a,n");

        // sigmaOps[54] += 1.000000 t2_abab_vvoo("a,b,m,n") V_blks_["abab_oovv"]("i,j,a,b")
        // flops: o4v2: 1, o4v0: 1 | mem: o4v0: 2,
        sigmaOps[54]("m,i,n,j") = t2_abab_vvoo("a,b,m,n") * V_blks_["abab_oovv"]("i,j,a,b");

        // sigmaOps[55] += 1.000000 V_blks_["bbbb_oovo"]("j,i,a,n") t2_abab_vvoo("e,a,m,i")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sigmaOps[55]("e,m,j,n") = V_blks_["bbbb_oovo"]("j,i,a,n") * t2_abab_vvoo("e,a,m,i");

        // sigmaOps[56] += 1.000000 t2_aaaa_vvoo("a,e,m,i") V_blks_["abab_oovo"]("i,j,a,n")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sigmaOps[56]("e,m,j,n") = t2_aaaa_vvoo("a,e,m,i") * V_blks_["abab_oovo"]("i,j,a,n");

        // sigmaOps[57] += 1.000000 V_blks_["aaaa_oovo"]("j,i,a,m") t2_abab_vvoo("a,f,i,n")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sigmaOps[57]("f,j,m,n") = V_blks_["aaaa_oovo"]("j,i,a,m") * t2_abab_vvoo("a,f,i,n");

        // sigmaOps[58] += 1.000000 V_blks_["abba_oovo"]("i,j,a,m") t2_abab_vvoo("e,a,i,n")
        // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2,
        sigmaOps[58]("e,m,j,n") = V_blks_["abba_oovo"]("i,j,a,m") * t2_abab_vvoo("e,a,i,n");

        // sigmaOps[59] += 1.000000 t2_abab_vvoo("a,e,j,i") V_blks_["abab_oovv"]("j,i,a,b")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sigmaOps[59]("e,b") = t2_abab_vvoo("a,e,j,i") * V_blks_["abab_oovv"]("j,i,a,b");

        // sigmaOps[60] += 1.000000 V_blks_["abab_oovv"]("i,j,b,a") t2_abab_vvoo("e,a,i,j")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sigmaOps[60]("b,e") = V_blks_["abab_oovv"]("i,j,b,a") * t2_abab_vvoo("e,a,i,j");

        // sigmaOps[61] += 1.000000 t2_aaaa_vvoo("a,e,j,i") V_blks_["aaaa_oovv"]("j,i,a,b")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sigmaOps[61]("e,b") = t2_aaaa_vvoo("a,e,j,i") * V_blks_["aaaa_oovv"]("j,i,a,b");

        // sigmaOps[62] += 1.000000 t2_bbbb_vvoo("a,e,j,i") V_blks_["bbbb_oovv"]("j,i,a,b")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sigmaOps[62]("e,b") = t2_bbbb_vvoo("a,e,j,i") * V_blks_["bbbb_oovv"]("j,i,a,b");

        // sigmaOps[63] += 1.000000 V_blks_["baab_vovv"]("e,i,a,b") t2_abab_vvoo("a,b,i,m")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[63]("e,m") = V_blks_["baab_vovv"]("e,i,a,b") * t2_abab_vvoo("a,b,i,m");

        // sigmaOps[64] += 1.000000 t2_abab_vvoo("a,b,m,i") V_blks_["abab_vovv"]("e,i,a,b")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[64]("e,m") = t2_abab_vvoo("a,b,m,i") * V_blks_["abab_vovv"]("e,i,a,b");

        // sigmaOps[65] += 1.000000 t2_bbbb_vvoo("a,b,m,i") V_blks_["bbbb_vovv"]("e,i,a,b")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[65]("e,m") = t2_bbbb_vvoo("a,b,m,i") * V_blks_["bbbb_vovv"]("e,i,a,b");

        // sigmaOps[66] += 1.000000 V_blks_["aaaa_vovv"]("e,i,a,b") t2_aaaa_vvoo("a,b,m,i")
        // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[66]("e,m") = V_blks_["aaaa_vovv"]("e,i,a,b") * t2_aaaa_vvoo("a,b,m,i");

        // sigmaOps[67] += 1.000000 V_blks_["aaaa_oovv"]("j,i,a,b") t2_aaaa_vvoo("b,f,j,i") t2_aaaa_vvoo("a,e,m,n")
        // flops: o2v3: 2, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
        sigmaOps[67]("f,e,m,n") = V_blks_["aaaa_oovv"]("j,i,a,b") * t2_aaaa_vvoo("b,f,j,i") * t2_aaaa_vvoo("a,e,m,n");

        // sigmaOps[68] += 1.000000 F_blks_["aa_vv"]("e,a") t2_abab_vvoo("a,f,m,n")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[68]("e,f,m,n") = F_blks_["aa_vv"]("e,a") * t2_abab_vvoo("a,f,m,n");

        // sigmaOps[69] += 1.000000 t2_bbbb_vvoo("a,f,m,n") F_blks_["bb_vv"]("e,a")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[69]("f,e,m,n") = t2_bbbb_vvoo("a,f,m,n") * F_blks_["bb_vv"]("e,a");

        // sigmaOps[70] += 1.000000 F_blks_["bb_vv"]("f,a") t2_abab_vvoo("e,a,m,n")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[70]("e,f,m,n") = F_blks_["bb_vv"]("f,a") * t2_abab_vvoo("e,a,m,n");

        // sigmaOps[71] += 1.000000 t2_bbbb_vvoo("b,f,j,i") V_blks_["bbbb_oovv"]("j,i,a,b")
        // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
        sigmaOps[71]("f,a") = t2_bbbb_vvoo("b,f,j,i") * V_blks_["bbbb_oovv"]("j,i,a,b");

        // sigmaOps[72] += 1.000000 F_blks_["aa_vv"]("e,a") t2_aaaa_vvoo("a,f,m,n")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[72]("e,f,m,n") = F_blks_["aa_vv"]("e,a") * t2_aaaa_vvoo("a,f,m,n");

        // sigmaOps[73] += 1.000000 V_blks_["abab_oovv"]("j,i,a,b") t2_abab_vvoo("a,b,m,i")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[73]("j,m") = V_blks_["abab_oovv"]("j,i,a,b") * t2_abab_vvoo("a,b,m,i");

        // sigmaOps[74] += 1.000000 V_blks_["abab_oovv"]("i,j,a,b") t2_abab_vvoo("a,b,i,m")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[74]("j,m") = V_blks_["abab_oovv"]("i,j,a,b") * t2_abab_vvoo("a,b,i,m");

        // sigmaOps[75] += 1.000000 V_blks_["aaaa_oovv"]("j,i,a,b") t2_aaaa_vvoo("a,b,m,i")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[75]("j,m") = V_blks_["aaaa_oovv"]("j,i,a,b") * t2_aaaa_vvoo("a,b,m,i");

        // sigmaOps[76] += 1.000000 V_blks_["bbbb_oovv"]("j,i,a,b") t2_bbbb_vvoo("a,b,m,i")
        // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
        sigmaOps[76]("j,m") = V_blks_["bbbb_oovv"]("j,i,a,b") * t2_bbbb_vvoo("a,b,m,i");

        // sigmaOps[77] += 1.000000 V_blks_["aaaa_oovo"]("j,i,a,m") t2_aaaa_vvoo("a,e,j,i")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[77]("e,m") = V_blks_["aaaa_oovo"]("j,i,a,m") * t2_aaaa_vvoo("a,e,j,i");

        // sigmaOps[78] += 1.000000 V_blks_["abab_oovo"]("j,i,a,m") t2_abab_vvoo("a,e,j,i")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[78]("e,m") = V_blks_["abab_oovo"]("j,i,a,m") * t2_abab_vvoo("a,e,j,i");

        // sigmaOps[79] += 1.000000 V_blks_["abba_oovo"]("j,i,a,m") t2_abab_vvoo("e,a,j,i")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[79]("e,m") = V_blks_["abba_oovo"]("j,i,a,m") * t2_abab_vvoo("e,a,j,i");

        // sigmaOps[80] += 1.000000 V_blks_["bbbb_oovo"]("j,i,a,m") t2_bbbb_vvoo("a,e,j,i")
        // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[80]("e,m") = V_blks_["bbbb_oovo"]("j,i,a,m") * t2_bbbb_vvoo("a,e,j,i");

        // sigmaOps[81] += 1.000000 F_blks_["aa_oo"]("i,m") t2_abab_vvoo("e,f,i,n")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[81]("e,f,m,n") = F_blks_["aa_oo"]("i,m") * t2_abab_vvoo("e,f,i,n");

        // sigmaOps[82] += 1.000000 F_blks_["aa_ov"]("i,a") t2_abab_vvoo("a,f,m,n")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sigmaOps[82]("f,i,m,n") = F_blks_["aa_ov"]("i,a") * t2_abab_vvoo("a,f,m,n");

        // sigmaOps[83] += 1.000000 F_blks_["bb_oo"]("i,n") t2_bbbb_vvoo("e,f,m,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[83]("e,f,n,m") = F_blks_["bb_oo"]("i,n") * t2_bbbb_vvoo("e,f,m,i");

        // sigmaOps[84] += 1.000000 t2_bbbb_vvoo("a,e,m,n") F_blks_["bb_ov"]("i,a")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sigmaOps[84]("e,m,n,i") = t2_bbbb_vvoo("a,e,m,n") * F_blks_["bb_ov"]("i,a");

        // sigmaOps[85] += 1.000000 t2_abab_vvoo("e,a,m,n") F_blks_["bb_ov"]("i,a")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sigmaOps[85]("e,m,n,i") = t2_abab_vvoo("e,a,m,n") * F_blks_["bb_ov"]("i,a");

        // sigmaOps[86] += 1.000000 t2_aaaa_vvoo("e,f,m,i") F_blks_["aa_oo"]("i,n")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[86]("e,f,m,n") = t2_aaaa_vvoo("e,f,m,i") * F_blks_["aa_oo"]("i,n");

        // sigmaOps[87] += 1.000000 F_blks_["bb_oo"]("i,n") t2_abab_vvoo("e,f,m,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[87]("e,f,m,n") = F_blks_["bb_oo"]("i,n") * t2_abab_vvoo("e,f,m,i");

        // sigmaOps[88] += 1.000000 F_blks_["aa_ov"]("i,a") t2_aaaa_vvoo("a,e,m,n")
        // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2,
        sigmaOps[88]("e,i,m,n") = F_blks_["aa_ov"]("i,a") * t2_aaaa_vvoo("a,e,m,n");

        // sigmaOps[89] += 1.000000 F_blks_["bb_ov"]("i,a") t2_bbbb_vvoo("a,e,m,i")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[89]("e,m") = F_blks_["bb_ov"]("i,a") * t2_bbbb_vvoo("a,e,m,i");

        // sigmaOps[90] += 1.000000 F_blks_["aa_ov"]("i,a") t2_abab_vvoo("a,e,i,m")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[90]("e,m") = F_blks_["aa_ov"]("i,a") * t2_abab_vvoo("a,e,i,m");

        // sigmaOps[91] += 1.000000 F_blks_["aa_ov"]("i,a") t2_aaaa_vvoo("a,e,m,i")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[91]("e,m") = F_blks_["aa_ov"]("i,a") * t2_aaaa_vvoo("a,e,m,i");

        // sigmaOps[92] += 1.000000 t2_abab_vvoo("e,a,m,i") F_blks_["bb_ov"]("i,a")
        // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
        sigmaOps[92]("e,m") = t2_abab_vvoo("e,a,m,i") * F_blks_["bb_ov"]("i,a");

        // sigmaOps[94] += 1.000000 tempOp14_abab_vvoo("a,f,i,m") t2_abab_vvoo("a,e,i,n") 1 V_blks_["abab_oovv"]("j,i,b,a") t2_bbbb_vvoo("a,f,n,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[94]("f,e,m,n") = sigmaOps[13]("a,f,i,m") * t2_abab_vvoo("a,e,i,n");

        // sigmaOps[95] += 1.000000 tempOp14_abab_vvoo("a,f,i,n") t2_aaaa_vvoo("a,e,m,i") 1 V_blks_["abab_oovv"]("j,i,b,a") t2_bbbb_vvoo("a,f,n,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[95]("e,f,m,n") = sigmaOps[13]("a,f,i,n") * t2_aaaa_vvoo("a,e,m,i");

        // sigmaOps[93] += 1.000000 t2_aaaa_vvoo("a,e,n,i") tempOp15_aaaa_vvoo("f,a,m,i") 1 t2_abab_vvoo("e,a,n,i") V_blks_["abab_oovv"]("j,i,b,a")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[93]("e,f,n,m") = t2_aaaa_vvoo("a,e,n,i") * sigmaOps[14]("f,a,m,i");

        // sigmaOps[96] += 1.000000 t2_abab_vvoo("b,f,j,n") tempOp15_aaaa_vvoo("e,b,m,j") 1 t2_abab_vvoo("e,a,n,i") V_blks_["abab_oovv"]("j,i,b,a")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[96]("e,f,m,n") = t2_abab_vvoo("b,f,j,n") * sigmaOps[14]("e,b,m,j");

        // sigmaOps[98] += 1.000000 t2_aaaa_vvoo("b,f,m,j") tempOp16_aaaa_vvoo("b,e,j,n") 1 V_blks_["aaaa_oovv"]("j,i,a,b") t2_aaaa_vvoo("a,e,n,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[98]("f,e,m,n") = t2_aaaa_vvoo("b,f,m,j") * sigmaOps[15]("b,e,j,n");

        // sigmaOps[100] += 1.000000 t2_abab_vvoo("b,f,j,n") tempOp16_aaaa_vvoo("b,e,j,m") 1 V_blks_["aaaa_oovv"]("j,i,a,b") t2_aaaa_vvoo("a,e,n,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[100]("e,f,m,n") = t2_abab_vvoo("b,f,j,n") * sigmaOps[15]("b,e,j,m");

        // sigmaOps[97] += 1.000000 t2_abab_vvoo("f,b,m,j") tempOp17_abab_vvoo("e,b,n,j") 1 V_blks_["bbbb_oovv"]("j,i,a,b") t2_abab_vvoo("e,a,n,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[97]("f,e,m,n") = t2_abab_vvoo("f,b,m,j") * sigmaOps[16]("e,b,n,j");

        // sigmaOps[99] += 1.000000 t2_bbbb_vvoo("b,f,n,j") tempOp17_abab_vvoo("e,b,m,j") 1 V_blks_["bbbb_oovv"]("j,i,a,b") t2_abab_vvoo("e,a,n,i")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[99]("e,f,m,n") = t2_bbbb_vvoo("b,f,n,j") * sigmaOps[16]("e,b,m,j");

        // sigmaOps[101] += 1.000000 t2_bbbb_vvoo("b,f,m,j") tempOp18_bbbb_vvoo("e,b,n,j") 1 t2_bbbb_vvoo("a,f,n,i") V_blks_["bbbb_oovv"]("j,i,a,b")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[101]("f,e,m,n") = t2_bbbb_vvoo("b,f,m,j") * sigmaOps[17]("e,b,n,j");

        // sigmaOps[103] += 1.000000 tempOp18_bbbb_vvoo("b,d,j,l") t2_abab_vvoo("a,d,i,l") 1 t2_bbbb_vvoo("a,f,n,i") V_blks_["bbbb_oovv"]("j,i,a,b")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[103]("a,b,i,j") = sigmaOps[17]("b,d,j,l") * t2_abab_vvoo("a,d,i,l");

        // sigmaOps[102] += 1.000000 t2_abab_vvoo("b,f,j,m") tempOp19_abab_vvoo("b,e,j,n") 1 V_blks_["aaaa_oovv"]("j,i,a,b") t2_abab_vvoo("a,f,i,n")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[102]("f,e,m,n") = t2_abab_vvoo("b,f,j,m") * sigmaOps[18]("b,e,j,n");

        // sigmaOps[104] += 1.000000 tempOp19_abab_vvoo("d,b,l,j") t2_aaaa_vvoo("d,a,i,l") 1 V_blks_["aaaa_oovv"]("j,i,a,b") t2_abab_vvoo("a,f,i,n")
        // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[104]("a,b,i,j") = sigmaOps[18]("d,b,l,j") * t2_aaaa_vvoo("d,a,i,l");

        // sigmaOps[105] += 1.000000 t2_aaaa_vvoo("e,f,j,i") tempOp43_aaaa_oooo("m,n,j,i") 1 t2_aaaa_vvoo("a,b,m,n") V_blks_["aaaa_oovv"]("j,i,a,b")
        // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[105]("e,f,m,n") = t2_aaaa_vvoo("e,f,j,i") * sigmaOps[42]("m,n,j,i");

        // sigmaOps[106] += 1.000000 t2_bbbb_vvoo("e,f,j,i") tempOp44_bbbb_oooo("m,n,j,i") 1 t2_bbbb_vvoo("a,b,m,n") V_blks_["bbbb_oovv"]("j,i,a,b")
        // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[106]("e,f,m,n") = t2_bbbb_vvoo("e,f,j,i") * sigmaOps[43]("m,n,j,i");

        // sigmaOps[108] += 1.000000 t2_bbbb_vvoo("a,e,m,n") tempOp60_bb_vv("f,a") 1 t2_abab_vvoo("a,e,j,i") V_blks_["abab_oovv"]("j,i,a,b")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[108]("e,f,m,n") = t2_bbbb_vvoo("a,e,m,n") * sigmaOps[59]("f,a");

        // sigmaOps[112] += 1.000000 tempOp60_bb_vv("f,a") t2_abab_vvoo("e,a,m,n") 1 t2_abab_vvoo("a,e,j,i") V_blks_["abab_oovv"]("j,i,a,b")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[112]("e,f,m,n") = sigmaOps[59]("f,a") * t2_abab_vvoo("e,a,m,n");

        // sigmaOps[107] += 1.000000 t2_aaaa_vvoo("a,e,m,n") tempOp61_aa_vv("a,f") 1 V_blks_["abab_oovv"]("i,j,b,a") t2_abab_vvoo("e,a,i,j")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[107]("e,f,m,n") = t2_aaaa_vvoo("a,e,m,n") * sigmaOps[60]("a,f");

        // sigmaOps[111] += 1.000000 tempOp61_aa_vv("b,e") t2_abab_vvoo("b,f,m,n") 1 V_blks_["abab_oovv"]("i,j,b,a") t2_abab_vvoo("e,a,i,j")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[111]("e,f,m,n") = sigmaOps[60]("b,e") * t2_abab_vvoo("b,f,m,n");

        // sigmaOps[109] += 1.000000 tempOp62_aa_vv("e,b") t2_aaaa_vvoo("b,f,m,n") 1 t2_aaaa_vvoo("a,e,j,i") V_blks_["aaaa_oovv"]("j,i,a,b")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[109]("e,f,m,n") = sigmaOps[61]("e,b") * t2_aaaa_vvoo("b,f,m,n");

        // sigmaOps[110] += 1.000000 tempOp62_aa_vv("e,b") t2_abab_vvoo("b,f,m,n") 1 t2_aaaa_vvoo("a,e,j,i") V_blks_["aaaa_oovv"]("j,i,a,b")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[110]("e,f,m,n") = sigmaOps[61]("e,b") * t2_abab_vvoo("b,f,m,n");

        // sigmaOps[113] += 1.000000 tempOp63_bb_vv("e,b") t2_bbbb_vvoo("b,f,m,n") 1 t2_bbbb_vvoo("a,e,j,i") V_blks_["bbbb_oovv"]("j,i,a,b")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[113]("e,f,m,n") = sigmaOps[62]("e,b") * t2_bbbb_vvoo("b,f,m,n");

        // sigmaOps[114] += 1.000000 tempOp63_bb_vv("b,d") t2_abab_vvoo("a,d,i,j") 1 t2_bbbb_vvoo("a,e,j,i") V_blks_["bbbb_oovv"]("j,i,a,b")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[114]("a,b,i,j") = sigmaOps[62]("b,d") * t2_abab_vvoo("a,d,i,j");

        // sigmaOps[115] += 1.000000 t2_bbbb_vvoo("a,e,m,n") tempOp72_bb_vv("f,a") 1 t2_bbbb_vvoo("b,f,j,i") V_blks_["bbbb_oovv"]("j,i,a,b")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[115]("e,f,m,n") = t2_bbbb_vvoo("a,e,m,n") * sigmaOps[71]("f,a");

        // sigmaOps[116] += 1.000000 t2_abab_vvoo("e,a,m,n") tempOp72_bb_vv("f,a") 1 t2_bbbb_vvoo("b,f,j,i") V_blks_["bbbb_oovv"]("j,i,a,b")
        // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[116]("e,f,m,n") = t2_abab_vvoo("e,a,m,n") * sigmaOps[71]("f,a");

        // sigmaOps[117] += 1.000000 t2_aaaa_vvoo("e,f,m,j") tempOp74_aa_oo("j,n") 1 V_blks_["abab_oovv"]("j,i,a,b") t2_abab_vvoo("a,b,m,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[117]("e,f,m,n") = t2_aaaa_vvoo("e,f,m,j") * sigmaOps[73]("j,n");

        // sigmaOps[121] += 1.000000 tempOp74_aa_oo("j,m") t2_abab_vvoo("e,f,j,n") 1 V_blks_["abab_oovv"]("j,i,a,b") t2_abab_vvoo("a,b,m,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[121]("e,f,m,n") = sigmaOps[73]("j,m") * t2_abab_vvoo("e,f,j,n");

        // sigmaOps[119] += 1.000000 tempOp75_bb_oo("j,n") t2_abab_vvoo("e,f,m,j") 1 V_blks_["abab_oovv"]("i,j,a,b") t2_abab_vvoo("a,b,i,m")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[119]("e,f,m,n") = sigmaOps[74]("j,n") * t2_abab_vvoo("e,f,m,j");

        // sigmaOps[124] += 1.000000 t2_bbbb_vvoo("e,f,m,j") tempOp75_bb_oo("j,n") 1 V_blks_["abab_oovv"]("i,j,a,b") t2_abab_vvoo("a,b,i,m")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[124]("e,f,m,n") = t2_bbbb_vvoo("e,f,m,j") * sigmaOps[74]("j,n");

        // sigmaOps[118] += 1.000000 tempOp76_aa_oo("j,n") t2_aaaa_vvoo("e,f,m,j") 1 V_blks_["aaaa_oovv"]("j,i,a,b") t2_aaaa_vvoo("a,b,m,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[118]("e,f,n,m") = sigmaOps[75]("j,n") * t2_aaaa_vvoo("e,f,m,j");

        // sigmaOps[120] += 1.000000 t2_abab_vvoo("e,f,j,n") tempOp76_aa_oo("j,m") 1 V_blks_["aaaa_oovv"]("j,i,a,b") t2_aaaa_vvoo("a,b,m,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[120]("e,f,m,n") = t2_abab_vvoo("e,f,j,n") * sigmaOps[75]("j,m");

        // sigmaOps[122] += 1.000000 tempOp77_bb_oo("j,n") t2_abab_vvoo("e,f,m,j") 1 V_blks_["bbbb_oovv"]("j,i,a,b") t2_bbbb_vvoo("a,b,m,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[122]("e,f,m,n") = sigmaOps[76]("j,n") * t2_abab_vvoo("e,f,m,j");

        // sigmaOps[123] += 1.000000 tempOp77_bb_oo("j,n") t2_bbbb_vvoo("e,f,m,j") 1 V_blks_["bbbb_oovv"]("j,i,a,b") t2_bbbb_vvoo("a,b,m,i")
        // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
        sigmaOps[123]("e,f,n,m") = sigmaOps[76]("j,n") * t2_bbbb_vvoo("e,f,m,j");
    }

    world_.gop.fence();
}

void hilbert::EOM_EE_CCSD::build_Hc_cH(size_t L) {

    // ensure that the integrals are transformed with t1 amplitudes
    if (!cc_wfn_->has_t1_integrals_) cc_wfn_->transform_integrals(true);

    /// build sigma vectors for this trial
    sigvec_blks_.clear();

    // reduce sigma into spins
    sigvec_blks_["sigmar_0"] = HelperD::makeTensor(world_, {L}, true);
    sigvec_blks_["sigmal_0"] = HelperD::makeTensor(world_, {L}, true);

    sigvec_blks_["sigmar_s_aa"] = HelperD::makeTensor(world_, {L, va_, oa_}, true);
    sigvec_blks_["sigmar_s_bb"] = HelperD::makeTensor(world_, {L, vb_, ob_}, true);
    sigvec_blks_["sigmal_s_aa"] = HelperD::makeTensor(world_, {L, va_, oa_}, true);
    sigvec_blks_["sigmal_s_bb"] = HelperD::makeTensor(world_, {L, vb_, ob_}, true);

    sigvec_blks_["sigmar_d_aaaa"] = HelperD::makeTensor(world_, {L, va_, va_, oa_, oa_}, true);
    sigvec_blks_["sigmar_d_abab"] = HelperD::makeTensor(world_, {L, va_, vb_, oa_, ob_}, true);
    sigvec_blks_["sigmar_d_bbbb"] = HelperD::makeTensor(world_, {L, vb_, vb_, ob_, ob_}, true);
    sigvec_blks_["sigmal_d_aaaa"] = HelperD::makeTensor(world_, {L, va_, va_, oa_, oa_}, true);
    sigvec_blks_["sigmal_d_abab"] = HelperD::makeTensor(world_, {L, va_, vb_, oa_, ob_}, true);
    sigvec_blks_["sigmal_d_bbbb"] = HelperD::makeTensor(world_, {L, vb_, vb_, ob_, ob_}, true);

    TA::TArrayD &t1_aa_vo =     cc_wfn_->amplitudes_["t1_aa"];
    TA::TArrayD &t1_bb_vo =     cc_wfn_->amplitudes_["t1_bb"];
    TA::TArrayD &t2_aaaa_vvoo = cc_wfn_->amplitudes_["t2_aaaa"];
    TA::TArrayD &t2_abab_vvoo = cc_wfn_->amplitudes_["t2_abab"];
    TA::TArrayD &t2_bbbb_vvoo = cc_wfn_->amplitudes_["t2_bbbb"];
    
    /// unpack the amplitude and sigma vectors
    // amplitude vector
    size_t idx = 0;
    auto &r0_L = evec_blks_["r0"];
    auto &r1_xaa_Lvo = evec_blks_["r1_aa"];
    auto &r1_xbb_Lvo = evec_blks_["r1_bb"];
    auto &r2_xaaaa_Lvvoo = evec_blks_["r2_aaaa"];
    auto &r2_xabab_Lvvoo = evec_blks_["r2_abab"];
    auto &r2_xbbbb_Lvvoo = evec_blks_["r2_bbbb"];

    // sigma vector
    idx = 0;
    auto &sigmar0_L = sigvec_blks_["sigmar_0"];
    auto &sigmal0_L = sigvec_blks_["sigmal_0"];

    auto &sigmar1_xaa_Lvo = sigvec_blks_["sigmar_s_aa"];
    auto &sigmar1_xbb_Lvo = sigvec_blks_["sigmar_s_bb"];
    auto &sigmal1_xaa_Lvo = sigvec_blks_["sigmal_s_aa"];
    auto &sigmal1_xbb_Lvo = sigvec_blks_["sigmal_s_bb"];

    auto &sigmar2_xaaaa_Lvvoo = sigvec_blks_["sigmar_d_aaaa"];
    auto &sigmar2_xabab_Lvvoo = sigvec_blks_["sigmar_d_abab"];
    auto &sigmar2_xbbbb_Lvvoo = sigvec_blks_["sigmar_d_bbbb"];
    auto &sigmal2_xaaaa_Lvvoo = sigvec_blks_["sigmal_d_aaaa"];
    auto &sigmal2_xabab_Lvvoo = sigvec_blks_["sigmal_d_abab"];
    auto &sigmal2_xbbbb_Lvvoo = sigvec_blks_["sigmal_d_bbbb"];

    /// initialize left operators from right
    auto &l0_L = evec_blks_["l0"];
    auto &l1_xaa_Lov = evec_blks_["l1_aa"];
    auto &l1_xbb_Lov = evec_blks_["l1_bb"];
    auto &l2_xaaaa_Loovv = evec_blks_["l2_aaaa"];
    auto &l2_xabab_Loovv = evec_blks_["l2_abab"];
    auto &l2_xbbbb_Loovv = evec_blks_["l2_bbbb"];

    TArrayMap &V_blks_ = cc_wfn_->V_blks_;
    TArrayMap &F_blks_ = cc_wfn_->F_blks_;
    TArrayMap &Id_blks_ = cc_wfn_->Id_blks_;

    {

        TA::TArrayD tempPerm_xaaaa_Lvvoo = HelperD::makeTensor(world_, {L, va_, va_, oa_, oa_}, true);
        TA::TArrayD tempPerm_xbbbb_Lvvoo = HelperD::makeTensor(world_, {L, vb_, vb_, ob_, ob_}, true);

        double scalar_0;
        double scalar_1;
        double scalar_2;
        double scalar_3;
        double scalar_4;
        double scalar_5;
        double scalar_6;
        double scalar_7;
        scalar_0 = dot(F_blks_["aa_oo"]("o0,o1"), Id_blks_["aa_oo"]("o0,o1"));
        scalar_1 = dot(F_blks_["bb_oo"]("o0,o1"), Id_blks_["bb_oo"]("o0,o1"));
        scalar_2 = dot(V_blks_["aaaa_oooo"]("o0,o1,o2,o3"), Id_blks_["aaaa_oooo"]("o0,o1,o2,o3"));
        scalar_3 = dot(V_blks_["abab_oooo"]("o0,o1,o2,o3"), Id_blks_["abab_oooo"]("o0,o1,o2,o3"));
        scalar_4 = dot(V_blks_["bbbb_oooo"]("o0,o1,o2,o3"), Id_blks_["bbbb_oooo"]("o0,o1,o2,o3"));
        scalar_5 = dot(V_blks_["aaaa_oovv"]("j,i,a,b"), t2_aaaa_vvoo("a,b,j,i"));
        scalar_6 = dot(V_blks_["abab_oovv"]("j,i,a,b"), t2_abab_vvoo("a,b,j,i"));
        scalar_7 = dot(V_blks_["bbbb_oovv"]("j,i,a,b"), t2_bbbb_vvoo("a,b,j,i"));
        auto tempArray = vector<TA::TArrayD>(34);

        /// ****** p†q ****** ///

        {

            // tempArray[0] += 1.000000 l2_xabab_Loovv("I,i,j,a,b") t2_abab_vvoo("a,c,i,j")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempArray[0]("I,b,c") = l2_xabab_Loovv("I,i,j,a,b") * t2_abab_vvoo("a,c,i,j");

            // sigmal1_xaa_Lvo += 0.500000 <m,b||e,c>_abab l2_abab(i,j,a,b) t2_abab(a,c,i,j)
            // +                  0.500000 <m,b||e,c>_abab l2_abab(j,i,a,b) t2_abab(a,c,j,i)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= tempArray[0]("I,b,c") * V_blks_["baab_vovv"]("b,m,e,c");

            // sigmal1_xbb_Lvo += -0.500000 <m,b||c,e>_bbbb l2_abab(i,j,a,b) t2_abab(a,c,i,j)
            // +                  -0.500000 <m,b||c,e>_bbbb l2_abab(j,i,a,b) t2_abab(a,c,j,i)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += tempArray[0]("I,b,c") * V_blks_["bbbb_vovv"]("b,m,c,e");

            // sigmal2_xabab_Lvvoo += -0.500000 <m,n||e,b>_abab l2_abab(i,j,a,f) t2_abab(a,b,i,j)
            // +                      -0.500000 <m,n||e,b>_abab l2_abab(j,i,a,f) t2_abab(a,b,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= tempArray[0]("I,f,b") * V_blks_["abab_oovv"]("m,n,e,b");

            // tempPerm_xbbbb_Lvvoo += 0.500000 P(e,f) <m,n||b,e>_bbbb l2_abab(i,j,a,f) t2_abab(a,b,i,j)
            // +                      0.500000 P(e,f) <m,n||b,e>_bbbb l2_abab(j,i,a,f) t2_abab(a,b,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") = tempArray[0]("I,f,b") * V_blks_["bbbb_oovv"]("m,n,b,e");

            // tempArray[1] += 1.000000 t2_abab_vvoo("c,a,i,j") l2_xabab_Loovv("I,i,j,b,a")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempArray[1]("I,c,b") = t2_abab_vvoo("c,a,i,j") * l2_xabab_Loovv("I,i,j,b,a");

            // sigmal1_xaa_Lvo += -0.500000 <m,b||c,e>_aaaa l2_abab(i,j,b,a) t2_abab(c,a,i,j)
            // +                  -0.500000 <m,b||c,e>_aaaa l2_abab(j,i,b,a) t2_abab(c,a,j,i)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += tempArray[1]("I,c,b") * V_blks_["aaaa_vovv"]("b,m,c,e");

            // sigmal1_xbb_Lvo += 0.500000 <b,m||c,e>_abab l2_abab(i,j,b,a) t2_abab(c,a,i,j)
            // +                  0.500000 <b,m||c,e>_abab l2_abab(j,i,b,a) t2_abab(c,a,j,i)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += tempArray[1]("I,c,b") * V_blks_["abab_vovv"]("b,m,c,e");

            // tempPerm_xaaaa_Lvvoo += 0.500000 P(e,f) <m,n||b,e>_aaaa l2_abab(i,j,f,a) t2_abab(b,a,i,j)
            // +                      0.500000 P(e,f) <m,n||b,e>_aaaa l2_abab(j,i,f,a) t2_abab(b,a,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") = tempArray[1]("I,b,f") * V_blks_["aaaa_oovv"]("m,n,b,e");

            // sigmal2_xabab_Lvvoo += -0.500000 <m,n||b,f>_abab l2_abab(i,j,e,a) t2_abab(b,a,i,j)
            // +                      -0.500000 <m,n||b,f>_abab l2_abab(j,i,e,a) t2_abab(b,a,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= tempArray[1]("I,b,e") * V_blks_["abab_oovv"]("m,n,b,f");

            // tempArray[2] += 0.500000 l2_xbbbb_Loovv("I,i,j,a,f") t2_bbbb_vvoo("b,a,i,j")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempArray[2]("I,f,b") = 0.500000 * l2_xbbbb_Loovv("I,i,j,a,f") * t2_bbbb_vvoo("b,a,i,j");

            // sigmal2_xabab_Lvvoo += 0.500000 <m,n||e,b>_abab l2_bbbb(i,j,a,f) t2_bbbb(b,a,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += tempArray[2]("I,f,b") * V_blks_["abab_oovv"]("m,n,e,b");

            // tempPerm_xbbbb_Lvvoo += -0.500000 P(e,f) <m,n||b,e>_bbbb l2_bbbb(i,j,a,f) t2_bbbb(b,a,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= tempArray[2]("I,f,b") * V_blks_["bbbb_oovv"]("m,n,b,e");

            // tempArray[3] += 0.500000 l2_xaaaa_Loovv("I,i,j,a,f") t2_aaaa_vvoo("b,a,i,j")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempArray[3]("I,f,b") = 0.500000 * l2_xaaaa_Loovv("I,i,j,a,f") * t2_aaaa_vvoo("b,a,i,j");

            // tempPerm_xaaaa_Lvvoo += -0.500000 P(e,f) <m,n||b,e>_aaaa l2_aaaa(i,j,a,f) t2_aaaa(b,a,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= tempArray[3]("I,f,b") * V_blks_["aaaa_oovv"]("m,n,b,e");

            // sigmal2_xabab_Lvvoo += 0.500000 <m,n||b,f>_abab l2_aaaa(i,j,a,e) t2_aaaa(b,a,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += tempArray[3]("I,e,b") * V_blks_["abab_oovv"]("m,n,b,f");

            // tempArray[4] += 0.500000 t2_aaaa_vvoo("c,a,i,j") l2_xaaaa_Loovv("I,i,j,b,a")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempArray[4]("I,c,b") = 0.500000 * t2_aaaa_vvoo("c,a,i,j") * l2_xaaaa_Loovv("I,i,j,b,a");

            // sigmal1_xaa_Lvo += -0.500000 <m,b||c,e>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(c,a,i,j)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += tempArray[4]("I,c,b") * V_blks_["aaaa_vovv"]("b,m,c,e");

            // sigmal1_xbb_Lvo += 0.500000 <b,m||c,e>_abab l2_aaaa(i,j,b,a) t2_aaaa(c,a,i,j)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += tempArray[4]("I,c,b") * V_blks_["abab_vovv"]("b,m,c,e");

            // tempArray[5] += 0.500000 t2_bbbb_vvoo("c,a,i,j") l2_xbbbb_Loovv("I,i,j,b,a")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempArray[5]("I,c,b") = 0.500000 * t2_bbbb_vvoo("c,a,i,j") * l2_xbbbb_Loovv("I,i,j,b,a");

            // sigmal1_xaa_Lvo += 0.500000 <m,b||e,c>_abab l2_bbbb(i,j,b,a) t2_bbbb(c,a,i,j)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= tempArray[5]("I,c,b") * V_blks_["baab_vovv"]("b,m,e,c");

            // sigmal1_xbb_Lvo += -0.500000 <m,b||c,e>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(c,a,i,j)
            // flops: o1v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += tempArray[5]("I,c,b") * V_blks_["bbbb_vovv"]("b,m,c,e");

            // tempArray[6] += 0.500000 r2_xbbbb_Lvvoo("I,b,f,j,i") V_blks_["bbbb_oovv"]("j,i,a,b")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempArray[6]("I,f,a") = 0.500000 * r2_xbbbb_Lvvoo("I,b,f,j,i") * V_blks_["bbbb_oovv"]("j,i,a,b");

            // sigmar2_xabab_Lvvoo += 0.500000 <j,i||a,b>_bbbb r2_bbbb(b,f,j,i) t2_abab(e,a,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += tempArray[6]("I,f,a") * t2_abab_vvoo("e,a,m,n");

            // sigmal2_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") += tempPerm_xbbbb_Lvvoo("I,e,f,m,n");
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") -= tempPerm_xbbbb_Lvvoo("I,f,e,m,n");

            // tempPerm_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xbbbb_Lvvoo = TArrayD(world_, tempPerm_xbbbb_Lvvoo.trange());
            tempPerm_xbbbb_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xbbbb_Lvvoo += -0.500000 P(e,f) <j,i||a,b>_bbbb r2_bbbb(b,f,j,i) t2_bbbb(a,e,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") = -1.0 * tempArray[6]("I,f,a") * t2_bbbb_vvoo("a,e,m,n");

            // tempArray[7] += 1.000000 V_blks_["abab_oovv"]("j,i,b,a") r2_xabab_Lvvoo("I,b,f,j,i")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempArray[7]("I,a,f") = V_blks_["abab_oovv"]("j,i,b,a") * r2_xabab_Lvvoo("I,b,f,j,i");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||b,a>_abab r2_abab(b,f,j,i) t2_abab(e,a,m,n)
            // +                      -0.500000 <i,j||b,a>_abab r2_abab(b,f,i,j) t2_abab(e,a,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= tempArray[7]("I,a,f") * t2_abab_vvoo("e,a,m,n");

            // tempPerm_xbbbb_Lvvoo += 0.500000 P(e,f) <j,i||b,a>_abab r2_abab(b,f,j,i) t2_bbbb(a,e,m,n)
            // +                      0.500000 P(e,f) <i,j||b,a>_abab r2_abab(b,f,i,j) t2_bbbb(a,e,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += tempArray[7]("I,a,f") * t2_bbbb_vvoo("a,e,m,n");

            // tempArray[8] += 1.000000 r2_xabab_Lvvoo("I,f,b,j,i") V_blks_["abab_oovv"]("j,i,a,b")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempArray[8]("I,f,a") = r2_xabab_Lvvoo("I,f,b,j,i") * V_blks_["abab_oovv"]("j,i,a,b");

            // sigmal2_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") += tempPerm_xaaaa_Lvvoo("I,e,f,m,n");
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") -= tempPerm_xaaaa_Lvvoo("I,f,e,m,n");

            // tempPerm_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xaaaa_Lvvoo = TArrayD(world_, tempPerm_xaaaa_Lvvoo.trange());
            tempPerm_xaaaa_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xaaaa_Lvvoo += 0.500000 P(e,f) <j,i||a,b>_abab r2_abab(f,b,j,i) t2_aaaa(a,e,m,n)
            // +                      0.500000 P(e,f) <i,j||a,b>_abab r2_abab(f,b,i,j) t2_aaaa(a,e,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") = tempArray[8]("I,f,a") * t2_aaaa_vvoo("a,e,m,n");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||a,b>_abab r2_abab(e,b,j,i) t2_abab(a,f,m,n)
            // +                      -0.500000 <i,j||a,b>_abab r2_abab(e,b,i,j) t2_abab(a,f,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= tempArray[8]("I,e,a") * t2_abab_vvoo("a,f,m,n");

            // tempArray[9] += 0.500000 r2_xaaaa_Lvvoo("I,b,f,j,i") V_blks_["aaaa_oovv"]("j,i,a,b")
            // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempArray[9]("I,f,a") = 0.500000 * r2_xaaaa_Lvvoo("I,b,f,j,i") * V_blks_["aaaa_oovv"]("j,i,a,b");

            // tempPerm_xaaaa_Lvvoo += -0.500000 P(e,f) <j,i||a,b>_aaaa r2_aaaa(b,f,j,i) t2_aaaa(a,e,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= tempArray[9]("I,f,a") * t2_aaaa_vvoo("a,e,m,n");

            // sigmar2_xabab_Lvvoo += 0.500000 <j,i||a,b>_aaaa r2_aaaa(b,e,j,i) t2_abab(a,f,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += tempArray[9]("I,e,a") * t2_abab_vvoo("a,f,m,n");

            // tempArray[10] += 1.000000 t2_abab_vvoo("b,a,j,i") l2_xabab_Loovv("I,m,i,b,a")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempArray[10]("I,j,m") = t2_abab_vvoo("b,a,j,i") * l2_xabab_Loovv("I,m,i,b,a");

            // sigmal1_xaa_Lvo += -0.500000 f_aa(j,e) l2_abab(m,i,b,a) t2_abab(b,a,j,i)
            // +                  -0.500000 f_aa(j,e) l2_abab(m,i,a,b) t2_abab(a,b,j,i)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= tempArray[10]("I,j,m") * F_blks_["aa_ov"]("j,e");

            // sigmal1_xaa_Lvo += -0.500000 <m,k||e,j>_aaaa l2_abab(j,i,b,a) t2_abab(b,a,k,i)
            // +                  -0.500000 <m,k||e,j>_aaaa l2_abab(j,i,a,b) t2_abab(a,b,k,i)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= tempArray[10]("I,k,j") * V_blks_["aaaa_oovo"]("m,k,e,j");

            // sigmal1_xbb_Lvo += -0.500000 <k,m||j,e>_abab l2_abab(j,i,b,a) t2_abab(b,a,k,i)
            // +                  -0.500000 <k,m||j,e>_abab l2_abab(j,i,a,b) t2_abab(a,b,k,i)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += tempArray[10]("I,k,j") * V_blks_["abba_oovo"]("k,m,e,j");

            // sigmar2_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += tempPerm_xaaaa_Lvvoo("I,e,f,m,n");
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") -= tempPerm_xaaaa_Lvvoo("I,f,e,m,n");

            // tempPerm_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xaaaa_Lvvoo = TArrayD(world_, tempPerm_xaaaa_Lvvoo.trange());
            tempPerm_xaaaa_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xaaaa_Lvvoo += 0.500000 P(m,n) <n,j||e,f>_aaaa l2_abab(m,i,b,a) t2_abab(b,a,j,i)
            // +                      0.500000 P(m,n) <n,j||e,f>_aaaa l2_abab(m,i,a,b) t2_abab(a,b,j,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") = tempArray[10]("I,j,m") * V_blks_["aaaa_oovv"]("n,j,e,f");

            // sigmal2_xabab_Lvvoo += -0.500000 <j,n||e,f>_abab l2_abab(m,i,b,a) t2_abab(b,a,j,i)
            // +                      -0.500000 <j,n||e,f>_abab l2_abab(m,i,a,b) t2_abab(a,b,j,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= tempArray[10]("I,j,m") * V_blks_["abab_oovv"]("j,n,e,f");

            // tempArray[11] += 1.000000 t2_abab_vvoo("b,a,i,k") l2_xabab_Loovv("I,i,j,b,a")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempArray[11]("I,k,j") = t2_abab_vvoo("b,a,i,k") * l2_xabab_Loovv("I,i,j,b,a");

            // sigmal1_xaa_Lvo += -0.500000 <m,k||e,j>_abab l2_abab(i,j,b,a) t2_abab(b,a,i,k)
            // +                  -0.500000 <m,k||e,j>_abab l2_abab(i,j,a,b) t2_abab(a,b,i,k)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= tempArray[11]("I,k,j") * V_blks_["abab_oovo"]("m,k,e,j");

            // sigmal1_xbb_Lvo += -0.500000 f_bb(j,e) l2_abab(i,m,b,a) t2_abab(b,a,i,j)
            // +                  -0.500000 f_bb(j,e) l2_abab(i,m,a,b) t2_abab(a,b,i,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= tempArray[11]("I,j,m") * F_blks_["bb_ov"]("j,e");

            // sigmal1_xbb_Lvo += -0.500000 <m,k||e,j>_bbbb l2_abab(i,j,b,a) t2_abab(b,a,i,k)
            // +                  -0.500000 <m,k||e,j>_bbbb l2_abab(i,j,a,b) t2_abab(a,b,i,k)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= tempArray[11]("I,k,j") * V_blks_["bbbb_oovo"]("m,k,e,j");

            // sigmal2_xabab_Lvvoo += -0.500000 <m,j||e,f>_abab l2_abab(i,n,b,a) t2_abab(b,a,i,j)
            // +                      -0.500000 <m,j||e,f>_abab l2_abab(i,n,a,b) t2_abab(a,b,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= tempArray[11]("I,j,n") * V_blks_["abab_oovv"]("m,j,e,f");

            // sigmar2_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += tempPerm_xbbbb_Lvvoo("I,e,f,m,n");
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") -= tempPerm_xbbbb_Lvvoo("I,f,e,m,n");

            // tempPerm_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xbbbb_Lvvoo = TArrayD(world_, tempPerm_xbbbb_Lvvoo.trange());
            tempPerm_xbbbb_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xbbbb_Lvvoo += 0.500000 P(m,n) <n,j||e,f>_bbbb l2_abab(i,m,b,a) t2_abab(b,a,i,j)
            // +                      0.500000 P(m,n) <n,j||e,f>_bbbb l2_abab(i,m,a,b) t2_abab(a,b,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") = tempArray[11]("I,j,m") * V_blks_["bbbb_oovv"]("n,j,e,f");

            // tempArray[12] += 0.500000 l2_xbbbb_Loovv("I,m,i,b,a") t2_bbbb_vvoo("b,a,i,j")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempArray[12]("I,m,j") = 0.500000 * l2_xbbbb_Loovv("I,m,i,b,a") * t2_bbbb_vvoo("b,a,i,j");

            // sigmal1_xbb_Lvo += 0.500000 f_bb(j,e) l2_bbbb(m,i,b,a) t2_bbbb(b,a,i,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += tempArray[12]("I,m,j") * F_blks_["bb_ov"]("j,e");

            // sigmal2_xabab_Lvvoo += 0.500000 <m,j||e,f>_abab l2_bbbb(n,i,b,a) t2_bbbb(b,a,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += tempArray[12]("I,n,j") * V_blks_["abab_oovv"]("m,j,e,f");

            // tempPerm_xbbbb_Lvvoo += -0.500000 P(m,n) <n,j||e,f>_bbbb l2_bbbb(m,i,b,a) t2_bbbb(b,a,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= tempArray[12]("I,m,j") * V_blks_["bbbb_oovv"]("n,j,e,f");

            // tempArray[13] += 0.500000 l2_xaaaa_Loovv("I,m,i,b,a") t2_aaaa_vvoo("b,a,i,j")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempArray[13]("I,m,j") = 0.500000 * l2_xaaaa_Loovv("I,m,i,b,a") * t2_aaaa_vvoo("b,a,i,j");

            // sigmal1_xaa_Lvo += 0.500000 f_aa(j,e) l2_aaaa(m,i,b,a) t2_aaaa(b,a,i,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += tempArray[13]("I,m,j") * F_blks_["aa_ov"]("j,e");

            // tempPerm_xaaaa_Lvvoo += -0.500000 P(m,n) <n,j||e,f>_aaaa l2_aaaa(m,i,b,a) t2_aaaa(b,a,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= tempArray[13]("I,m,j") * V_blks_["aaaa_oovv"]("n,j,e,f");

            // sigmal2_xabab_Lvvoo += 0.500000 <j,n||e,f>_abab l2_aaaa(m,i,b,a) t2_aaaa(b,a,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += tempArray[13]("I,m,j") * V_blks_["abab_oovv"]("j,n,e,f");

            // tempArray[14] += 0.500000 t2_aaaa_vvoo("b,a,i,k") l2_xaaaa_Loovv("I,i,j,b,a")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempArray[14]("I,k,j") = 0.500000 * t2_aaaa_vvoo("b,a,i,k") * l2_xaaaa_Loovv("I,i,j,b,a");

            // sigmal1_xaa_Lvo += -0.500000 <m,k||e,j>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(b,a,i,k)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= tempArray[14]("I,k,j") * V_blks_["aaaa_oovo"]("m,k,e,j");

            // sigmal1_xbb_Lvo += -0.500000 <k,m||j,e>_abab l2_aaaa(i,j,b,a) t2_aaaa(b,a,i,k)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += tempArray[14]("I,k,j") * V_blks_["abba_oovo"]("k,m,e,j");

            // tempArray[15] += 0.500000 t2_bbbb_vvoo("b,a,i,k") l2_xbbbb_Loovv("I,i,j,b,a")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempArray[15]("I,k,j") = 0.500000 * t2_bbbb_vvoo("b,a,i,k") * l2_xbbbb_Loovv("I,i,j,b,a");

            // sigmal1_xaa_Lvo += -0.500000 <m,k||e,j>_abab l2_bbbb(i,j,b,a) t2_bbbb(b,a,i,k)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= tempArray[15]("I,k,j") * V_blks_["abab_oovo"]("m,k,e,j");

            // sigmal1_xbb_Lvo += -0.500000 <m,k||e,j>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(b,a,i,k)
            // flops: o3v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= tempArray[15]("I,k,j") * V_blks_["bbbb_oovo"]("m,k,e,j");

            // tempArray[16] += 1.000000 V_blks_["abab_oovv"]("j,i,a,b") r2_xabab_Lvvoo("I,a,b,j,n")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempArray[16]("I,i,n") = V_blks_["abab_oovv"]("j,i,a,b") * r2_xabab_Lvvoo("I,a,b,j,n");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||a,b>_abab r2_abab(a,b,j,n) t2_abab(e,f,m,i)
            // +                      -0.500000 <j,i||b,a>_abab r2_abab(b,a,j,n) t2_abab(e,f,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= tempArray[16]("I,i,n") * t2_abab_vvoo("e,f,m,i");

            // sigmal2_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") += tempPerm_xbbbb_Lvvoo("I,e,f,m,n");
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") -= tempPerm_xbbbb_Lvvoo("I,e,f,n,m");

            // tempPerm_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xbbbb_Lvvoo = TArrayD(world_, tempPerm_xbbbb_Lvvoo.trange());
            tempPerm_xbbbb_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xbbbb_Lvvoo += 0.500000 P(m,n) <j,i||a,b>_abab r2_abab(a,b,j,m) t2_bbbb(e,f,n,i)
            // +                      0.500000 P(m,n) <j,i||b,a>_abab r2_abab(b,a,j,m) t2_bbbb(e,f,n,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") = tempArray[16]("I,i,m") * t2_bbbb_vvoo("e,f,n,i");

            // tempArray[17] += 0.500000 r2_xbbbb_Lvvoo("I,a,b,n,j") V_blks_["bbbb_oovv"]("j,i,a,b")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempArray[17]("I,n,i") = 0.500000 * r2_xbbbb_Lvvoo("I,a,b,n,j") * V_blks_["bbbb_oovv"]("j,i,a,b");

            // sigmar2_xabab_Lvvoo += 0.500000 <j,i||a,b>_bbbb r2_bbbb(a,b,n,j) t2_abab(e,f,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += tempArray[17]("I,n,i") * t2_abab_vvoo("e,f,m,i");

            // tempPerm_xbbbb_Lvvoo += -0.500000 P(m,n) <j,i||a,b>_bbbb r2_bbbb(a,b,m,j) t2_bbbb(e,f,n,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= tempArray[17]("I,m,i") * t2_bbbb_vvoo("e,f,n,i");

            // tempArray[18] += 1.000000 V_blks_["abab_oovv"]("i,j,a,b") r2_xabab_Lvvoo("I,a,b,m,j")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempArray[18]("I,i,m") = V_blks_["abab_oovv"]("i,j,a,b") * r2_xabab_Lvvoo("I,a,b,m,j");

            // sigmal2_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") += tempPerm_xaaaa_Lvvoo("I,e,f,m,n");
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") -= tempPerm_xaaaa_Lvvoo("I,e,f,n,m");

            // tempPerm_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xaaaa_Lvvoo = TArrayD(world_, tempPerm_xaaaa_Lvvoo.trange());
            tempPerm_xaaaa_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xaaaa_Lvvoo += 0.500000 P(m,n) <i,j||a,b>_abab r2_abab(a,b,m,j) t2_aaaa(e,f,n,i)
            // +                      0.500000 P(m,n) <i,j||b,a>_abab r2_abab(b,a,m,j) t2_aaaa(e,f,n,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") = tempArray[18]("I,i,m") * t2_aaaa_vvoo("e,f,n,i");

            // sigmar2_xabab_Lvvoo += -0.500000 <i,j||a,b>_abab r2_abab(a,b,m,j) t2_abab(e,f,i,n)
            // +                      -0.500000 <i,j||b,a>_abab r2_abab(b,a,m,j) t2_abab(e,f,i,n)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= tempArray[18]("I,i,m") * t2_abab_vvoo("e,f,i,n");

            // tempArray[19] += 0.500000 V_blks_["aaaa_oovv"]("j,i,a,b") r2_xaaaa_Lvvoo("I,a,b,m,j")
            // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
            tempArray[19]("I,i,m") = 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaaa_Lvvoo("I,a,b,m,j");

            // tempPerm_xaaaa_Lvvoo += -0.500000 P(m,n) <j,i||a,b>_aaaa r2_aaaa(a,b,m,j) t2_aaaa(e,f,n,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= tempArray[19]("I,i,m") * t2_aaaa_vvoo("e,f,n,i");

            // sigmar2_xabab_Lvvoo += 0.500000 <j,i||a,b>_aaaa r2_aaaa(a,b,m,j) t2_abab(e,f,i,n)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += tempArray[19]("I,i,m") * t2_abab_vvoo("e,f,i,n");

            // tempArray[20] += 1.000000 r1_xbb_Lvo("I,b,i") V_blks_["bbbb_vovv"]("f,i,a,b")
            // flops: o1v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempArray[20]("I,f,a") = r1_xbb_Lvo("I,b,i") * V_blks_["bbbb_vovv"]("f,i,a,b");

            // sigmar2_xabab_Lvvoo += -1.000000 <i,f||a,b>_bbbb r1_bb(b,i) t2_abab(e,a,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += tempArray[20]("I,f,a") * t2_abab_vvoo("e,a,m,n");

            // sigmar2_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += tempPerm_xbbbb_Lvvoo("I,e,f,m,n");
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") -= tempPerm_xbbbb_Lvvoo("I,e,f,n,m");

            // tempPerm_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xbbbb_Lvvoo = TArrayD(world_, tempPerm_xbbbb_Lvvoo.trange());
            tempPerm_xbbbb_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xbbbb_Lvvoo += -1.000000 P(e,f) <i,e||a,b>_bbbb r1_bb(b,i) t2_bbbb(a,f,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") = tempArray[20]("I,e,a") * t2_bbbb_vvoo("a,f,m,n");

            // tempArray[21] += 1.000000 V_blks_["baab_vovv"]("f,i,b,a") r1_xaa_Lvo("I,b,i")
            // flops: o1v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempArray[21]("I,f,a") = V_blks_["baab_vovv"]("f,i,b,a") * r1_xaa_Lvo("I,b,i");

            // sigmar2_xabab_Lvvoo += 1.000000 <i,f||b,a>_abab r1_aa(b,i) t2_abab(e,a,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= tempArray[21]("I,f,a") * t2_abab_vvoo("e,a,m,n");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(e,f) <i,e||b,a>_abab r1_aa(b,i) t2_bbbb(a,f,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= tempArray[21]("I,e,a") * t2_bbbb_vvoo("a,f,m,n");

            // tempArray[22] += 1.000000 V_blks_["aaaa_vovv"]("e,i,a,b") r1_xaa_Lvo("I,b,i")
            // flops: o1v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempArray[22]("I,e,a") = V_blks_["aaaa_vovv"]("e,i,a,b") * r1_xaa_Lvo("I,b,i");

            // sigmar2_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += tempPerm_xaaaa_Lvvoo("I,e,f,m,n");
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") -= tempPerm_xaaaa_Lvvoo("I,e,f,n,m");

            // tempPerm_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xaaaa_Lvvoo = TArrayD(world_, tempPerm_xaaaa_Lvvoo.trange());
            tempPerm_xaaaa_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xaaaa_Lvvoo += -1.000000 P(e,f) <i,e||a,b>_aaaa r1_aa(b,i) t2_aaaa(a,f,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") = tempArray[22]("I,e,a") * t2_aaaa_vvoo("a,f,m,n");

            // sigmar2_xabab_Lvvoo += -1.000000 <i,e||a,b>_aaaa r1_aa(b,i) t2_abab(a,f,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += tempArray[22]("I,e,a") * t2_abab_vvoo("a,f,m,n");

            // tempArray[23] += 1.000000 r1_xbb_Lvo("I,b,i") V_blks_["abab_vovv"]("e,i,a,b")
            // flops: o1v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
            tempArray[23]("I,e,a") = r1_xbb_Lvo("I,b,i") * V_blks_["abab_vovv"]("e,i,a,b");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(e,f) <e,i||a,b>_abab r1_bb(b,i) t2_aaaa(a,f,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += tempArray[23]("I,e,a") * t2_aaaa_vvoo("a,f,m,n");

            // sigmar2_xabab_Lvvoo += 1.000000 <e,i||a,b>_abab r1_bb(b,i) t2_abab(a,f,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += tempArray[23]("I,e,a") * t2_abab_vvoo("a,f,m,n");

            // tempArray[24] += 1.000000 r0_L("I") sigmaOps[94]("f,e,m,n")
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempArray[24]("I,f,e,m,n") = r0_L("I") * sigmaOps[94]("f,e,m,n");

            // sigmar2_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += tempPerm_xbbbb_Lvvoo("I,e,f,m,n");
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") -= tempPerm_xbbbb_Lvvoo("I,f,e,m,n");

            // tempPerm_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xbbbb_Lvvoo = TArrayD(world_, tempPerm_xbbbb_Lvvoo.trange());
            tempPerm_xbbbb_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) <j,i||b,a>_abab r0 t2_bbbb(a,e,n,i) t2_abab(b,f,j,m)
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") = tempArray[24]("I,e,f,n,m");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) <i,j||a,b>_abab r0 t2_abab(a,e,i,n) t2_bbbb(b,f,m,j)
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += tempArray[24]("I,f,e,m,n");

            // tempArray[25] += 1.000000 sigmaOps[93]("e,f,n,m") r0_L("I")
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempArray[25]("I,e,f,n,m") = sigmaOps[93]("e,f,n,m") * r0_L("I");

            // sigmar2_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += tempPerm_xaaaa_Lvvoo("I,e,f,m,n");
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") -= tempPerm_xaaaa_Lvvoo("I,f,e,m,n");

            // tempPerm_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xaaaa_Lvvoo = TArrayD(world_, tempPerm_xaaaa_Lvvoo.trange());
            tempPerm_xaaaa_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) <i,j||a,b>_abab r0 t2_aaaa(a,e,n,i) t2_abab(f,b,m,j)
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") = tempArray[25]("I,e,f,n,m");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) <j,i||b,a>_abab r0 t2_abab(e,a,n,i) t2_aaaa(b,f,m,j)
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += tempArray[25]("I,f,e,m,n");

            // tempArray[26] += 1.000000 l1_xaa_Lov("I,i,a") t2_abab_vvoo("a,b,i,j")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempArray[26]("I,b,j") = l1_xaa_Lov("I,i,a") * t2_abab_vvoo("a,b,i,j");

            // sigmal1_xaa_Lvo += 1.000000 <m,j||e,b>_abab l1_aa(i,a) t2_abab(a,b,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += tempArray[26]("I,b,j") * V_blks_["abab_oovv"]("m,j,e,b");

            // sigmal1_xbb_Lvo += -1.000000 <m,j||b,e>_bbbb l1_aa(i,a) t2_abab(a,b,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= tempArray[26]("I,b,j") * V_blks_["bbbb_oovv"]("m,j,b,e");

            // tempArray[27] += 1.000000 t2_aaaa_vvoo("b,a,i,j") l1_xaa_Lov("I,i,a")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempArray[27]("I,b,j") = t2_aaaa_vvoo("b,a,i,j") * l1_xaa_Lov("I,i,a");

            // sigmal1_xaa_Lvo += 1.000000 <m,j||b,e>_aaaa l1_aa(i,a) t2_aaaa(b,a,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += tempArray[27]("I,b,j") * V_blks_["aaaa_oovv"]("m,j,b,e");

            // sigmal1_xbb_Lvo += -1.000000 <j,m||b,e>_abab l1_aa(i,a) t2_aaaa(b,a,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= tempArray[27]("I,b,j") * V_blks_["abab_oovv"]("j,m,b,e");

            // tempArray[28] += 1.000000 t2_abab_vvoo("b,a,j,i") l1_xbb_Lov("I,i,a")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempArray[28]("I,b,j") = t2_abab_vvoo("b,a,j,i") * l1_xbb_Lov("I,i,a");

            // sigmal1_xaa_Lvo += -1.000000 <m,j||b,e>_aaaa l1_bb(i,a) t2_abab(b,a,j,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= tempArray[28]("I,b,j") * V_blks_["aaaa_oovv"]("m,j,b,e");

            // sigmal1_xbb_Lvo += 1.000000 <j,m||b,e>_abab l1_bb(i,a) t2_abab(b,a,j,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += tempArray[28]("I,b,j") * V_blks_["abab_oovv"]("j,m,b,e");

            // tempArray[29] += 1.000000 t2_bbbb_vvoo("b,a,i,j") l1_xbb_Lov("I,i,a")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempArray[29]("I,b,j") = t2_bbbb_vvoo("b,a,i,j") * l1_xbb_Lov("I,i,a");

            // sigmal1_xaa_Lvo += -1.000000 <m,j||e,b>_abab l1_bb(i,a) t2_bbbb(b,a,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= tempArray[29]("I,b,j") * V_blks_["abab_oovv"]("m,j,e,b");

            // sigmal1_xbb_Lvo += 1.000000 <m,j||b,e>_bbbb l1_bb(i,a) t2_bbbb(b,a,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += tempArray[29]("I,b,j") * V_blks_["bbbb_oovv"]("m,j,b,e");

            // tempArray[30] += 1.000000 V_blks_["abab_oovv"]("i,j,a,b") r1_xbb_Lvo("I,b,j")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempArray[30]("I,a,i") = V_blks_["abab_oovv"]("i,j,a,b") * r1_xbb_Lvo("I,b,j");

            // sigmar1_xaa_Lvo += -1.000000 <i,j||a,b>_abab r1_bb(b,j) t2_aaaa(a,e,m,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") -= tempArray[30]("I,a,i") * t2_aaaa_vvoo("a,e,m,i");

            // sigmar1_xbb_Lvo += 1.000000 <i,j||a,b>_abab r1_bb(b,j) t2_abab(a,e,i,m)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") += tempArray[30]("I,a,i") * t2_abab_vvoo("a,e,i,m");

            // tempArray[31] += 1.000000 V_blks_["abab_oovv"]("j,i,b,a") r1_xaa_Lvo("I,b,j")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempArray[31]("I,a,i") = V_blks_["abab_oovv"]("j,i,b,a") * r1_xaa_Lvo("I,b,j");

            // sigmar1_xaa_Lvo += 1.000000 <j,i||b,a>_abab r1_aa(b,j) t2_abab(e,a,m,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") += tempArray[31]("I,a,i") * t2_abab_vvoo("e,a,m,i");

            // sigmar1_xbb_Lvo += -1.000000 <j,i||b,a>_abab r1_aa(b,j) t2_bbbb(a,e,m,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") -= tempArray[31]("I,a,i") * t2_bbbb_vvoo("a,e,m,i");

            // tempArray[32] += 1.000000 r1_xbb_Lvo("I,b,j") V_blks_["bbbb_oovv"]("j,i,a,b")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempArray[32]("I,a,i") = r1_xbb_Lvo("I,b,j") * V_blks_["bbbb_oovv"]("j,i,a,b");

            // sigmar1_xaa_Lvo += -1.000000 <j,i||a,b>_bbbb r1_bb(b,j) t2_abab(e,a,m,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") -= tempArray[32]("I,a,i") * t2_abab_vvoo("e,a,m,i");

            // sigmar1_xbb_Lvo += 1.000000 <j,i||a,b>_bbbb r1_bb(b,j) t2_bbbb(a,e,m,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") += tempArray[32]("I,a,i") * t2_bbbb_vvoo("a,e,m,i");

            // tempArray[33] += 1.000000 r1_xaa_Lvo("I,b,j") V_blks_["aaaa_oovv"]("j,i,a,b")
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            tempArray[33]("I,a,i") = r1_xaa_Lvo("I,b,j") * V_blks_["aaaa_oovv"]("j,i,a,b");

            // sigmar1_xaa_Lvo += 1.000000 <j,i||a,b>_aaaa r1_aa(b,j) t2_aaaa(a,e,m,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") += tempArray[33]("I,a,i") * t2_aaaa_vvoo("a,e,m,i");

            // sigmar1_xbb_Lvo += -1.000000 <j,i||a,b>_aaaa r1_aa(b,j) t2_abab(a,e,i,m)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") -= tempArray[33]("I,a,i") * t2_abab_vvoo("a,e,i,m");

            // sigmar0_L += 1.000000 f_aa(i,i) r0
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            sigmar0_L("I") += scalar_0 * r0_L("I");

            // sigmar0_L += 1.000000 f_bb(i,i) r0
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            sigmar0_L("I") += scalar_1 * r0_L("I");

            // sigmar0_L += -0.500000 <j,i||j,i>_abab r0
            // +            -0.500000 <i,j||i,j>_abab r0
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            sigmar0_L("I") -= scalar_3 * r0_L("I");

            // sigmar0_L += 1.000000 f_aa(i,a) r1_aa(a,i)
            // flops: o1v1L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmar0_L("I") += F_blks_["aa_ov"]("i,a") * r1_xaa_Lvo("I,a,i");

            // sigmar0_L += 1.000000 f_bb(i,a) r1_bb(a,i)
            // flops: o1v1L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmar0_L("I") += F_blks_["bb_ov"]("i,a") * r1_xbb_Lvo("I,a,i");

            // sigmar0_L += 0.250000 <j,i||a,b>_aaaa r2_aaaa(a,b,j,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmar0_L("I") += 0.250000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaaa_Lvvoo("I,a,b,j,i");

            // sigmar0_L += -0.500000 <j,i||j,i>_aaaa r0
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            sigmar0_L("I") -= 0.500000 * scalar_2 * r0_L("I");

            // sigmar0_L += -0.500000 <j,i||j,i>_bbbb r0
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            sigmar0_L("I") -= 0.500000 * scalar_4 * r0_L("I");

            // sigmar0_L += 0.250000 <j,i||a,b>_abab r2_abab(a,b,j,i)
            // +            0.250000 <i,j||a,b>_abab r2_abab(a,b,i,j)
            // +            0.250000 <j,i||b,a>_abab r2_abab(b,a,j,i)
            // +            0.250000 <i,j||b,a>_abab r2_abab(b,a,i,j)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmar0_L("I") += V_blks_["abab_oovv"]("j,i,a,b") * r2_xabab_Lvvoo("I,a,b,j,i");

            // sigmar0_L += 0.250000 <j,i||a,b>_bbbb r2_bbbb(a,b,j,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmar0_L("I") += 0.250000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r2_xbbbb_Lvvoo("I,a,b,j,i");

            // sigmar0_L += 0.250000 <j,i||a,b>_aaaa r0 t2_aaaa(a,b,j,i)
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            sigmar0_L("I") += 0.250000 * scalar_5 * r0_L("I");

            // sigmar0_L += 0.250000 <j,i||a,b>_abab r0 t2_abab(a,b,j,i)
            // +            0.250000 <i,j||a,b>_abab r0 t2_abab(a,b,i,j)
            // +            0.250000 <j,i||b,a>_abab r0 t2_abab(b,a,j,i)
            // +            0.250000 <i,j||b,a>_abab r0 t2_abab(b,a,i,j)
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            sigmar0_L("I") += scalar_6 * r0_L("I");

            // sigmar0_L += 0.250000 <j,i||a,b>_bbbb r0 t2_bbbb(a,b,j,i)
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            sigmar0_L("I") += 0.250000 * scalar_7 * r0_L("I");

            // sigmar1_xaa_Lvo += -0.500000 <j,i||m,a>_abab r0 t2_abab(e,a,j,i)
            // +                  -0.500000 <i,j||m,a>_abab r0 t2_abab(e,a,i,j)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") += sigmaOps[79]("e,m") * r0_L("I");

            // sigmar1_xaa_Lvo += -0.500000 <j,i||a,b>_aaaa r1_aa(e,j) t2_aaaa(a,b,m,i)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") -= 0.500000 * sigmaOps[75]("j,m") * r1_xaa_Lvo("I,e,j");

            // sigmar1_xaa_Lvo += 1.000000 f_aa(i,i) r1_aa(e,m)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmar1_xaa_Lvo("I,e,m") += scalar_0 * r1_xaa_Lvo("I,e,m");

            // sigmar1_xaa_Lvo += 1.000000 f_aa(e,a) r1_aa(a,m)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") += F_blks_["aa_vv"]("e,a") * r1_xaa_Lvo("I,a,m");

            // sigmar1_xaa_Lvo += -0.500000 <j,i||j,i>_abab r1_aa(e,m)
            // +                  -0.500000 <i,j||i,j>_abab r1_aa(e,m)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmar1_xaa_Lvo("I,e,m") -= scalar_3 * r1_xaa_Lvo("I,e,m");

            // sigmar1_xaa_Lvo += -0.500000 <j,i||a,b>_abab r1_aa(e,j) t2_abab(a,b,m,i)
            // +                  -0.500000 <j,i||b,a>_abab r1_aa(e,j) t2_abab(b,a,m,i)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") -= sigmaOps[73]("j,m") * r1_xaa_Lvo("I,e,j");

            // sigmar1_xaa_Lvo += -0.500000 <j,i||m,a>_abab r2_abab(e,a,j,i)
            // +                  -0.500000 <i,j||m,a>_abab r2_abab(e,a,i,j)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") += V_blks_["abba_oovo"]("j,i,a,m") * r2_xabab_Lvvoo("I,e,a,j,i");

            // sigmar1_xaa_Lvo += 1.000000 f_aa(e,m) r0
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") += F_blks_["aa_vo"]("e,m") * r0_L("I");

            // sigmar1_xaa_Lvo += -0.500000 <j,i||a,m>_aaaa r0 t2_aaaa(a,e,j,i)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") -= 0.500000 * sigmaOps[77]("e,m") * r0_L("I");

            // sigmar1_xaa_Lvo += 0.250000 <j,i||a,b>_abab r1_aa(e,m) t2_abab(a,b,j,i)
            // +                  0.250000 <i,j||a,b>_abab r1_aa(e,m) t2_abab(a,b,i,j)
            // +                  0.250000 <j,i||b,a>_abab r1_aa(e,m) t2_abab(b,a,j,i)
            // +                  0.250000 <i,j||b,a>_abab r1_aa(e,m) t2_abab(b,a,i,j)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmar1_xaa_Lvo("I,e,m") += scalar_6 * r1_xaa_Lvo("I,e,m");

            // sigmar1_xaa_Lvo += -0.500000 <j,i||j,i>_aaaa r1_aa(e,m)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmar1_xaa_Lvo("I,e,m") -= 0.500000 * scalar_2 * r1_xaa_Lvo("I,e,m");

            // sigmar1_xaa_Lvo += 0.500000 <e,i||a,b>_abab r0 t2_abab(a,b,m,i)
            // +                  0.500000 <e,i||b,a>_abab r0 t2_abab(b,a,m,i)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") += sigmaOps[64]("e,m") * r0_L("I");

            // sigmar1_xaa_Lvo += 1.000000 <e,i||m,a>_abab r1_bb(a,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") -= V_blks_["abba_vovo"]("e,i,a,m") * r1_xbb_Lvo("I,a,i");

            // sigmar1_xaa_Lvo += 1.000000 <i,e||a,m>_aaaa r1_aa(a,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") -= V_blks_["aaaa_vovo"]("e,i,a,m") * r1_xaa_Lvo("I,a,i");

            // sigmar1_xaa_Lvo += 0.250000 <j,i||a,b>_bbbb r1_aa(e,m) t2_bbbb(a,b,j,i)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmar1_xaa_Lvo("I,e,m") += 0.250000 * scalar_7 * r1_xaa_Lvo("I,e,m");

            // sigmar1_xaa_Lvo += -0.500000 <j,i||a,b>_aaaa r1_aa(b,m) t2_aaaa(a,e,j,i)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") -= 0.500000 * sigmaOps[61]("e,b") * r1_xaa_Lvo("I,b,m");

            // sigmar1_xaa_Lvo += 1.000000 f_bb(i,i) r1_aa(e,m)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmar1_xaa_Lvo("I,e,m") += scalar_1 * r1_xaa_Lvo("I,e,m");

            // sigmar1_xaa_Lvo += -1.000000 f_aa(i,m) r1_aa(e,i)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") -= F_blks_["aa_oo"]("i,m") * r1_xaa_Lvo("I,e,i");

            // sigmar1_xaa_Lvo += -1.000000 f_aa(i,a) r2_aaaa(a,e,m,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") -= F_blks_["aa_ov"]("i,a") * r2_xaaaa_Lvvoo("I,a,e,m,i");

            // sigmar1_xaa_Lvo += 0.250000 <j,i||a,b>_aaaa r1_aa(e,m) t2_aaaa(a,b,j,i)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmar1_xaa_Lvo("I,e,m") += 0.250000 * scalar_5 * r1_xaa_Lvo("I,e,m");

            // sigmar1_xaa_Lvo += -0.500000 <j,i||a,m>_aaaa r2_aaaa(a,e,j,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") -= 0.500000 * V_blks_["aaaa_oovo"]("j,i,a,m") * r2_xaaaa_Lvvoo("I,a,e,j,i");

            // sigmar1_xaa_Lvo += -0.500000 <i,e||a,b>_aaaa r2_aaaa(a,b,m,i)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") += 0.500000 * V_blks_["aaaa_vovv"]("e,i,a,b") * r2_xaaaa_Lvvoo("I,a,b,m,i");

            // sigmar1_xaa_Lvo += -0.500000 <i,e||a,b>_aaaa r0 t2_aaaa(a,b,m,i)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") += 0.500000 * sigmaOps[66]("e,m") * r0_L("I");

            // sigmar1_xaa_Lvo += 0.500000 <e,i||a,b>_abab r2_abab(a,b,m,i)
            // +                  0.500000 <e,i||b,a>_abab r2_abab(b,a,m,i)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") += V_blks_["abab_vovv"]("e,i,a,b") * r2_xabab_Lvvoo("I,a,b,m,i");

            // sigmar1_xaa_Lvo += -1.000000 f_aa(i,a) r0 t2_aaaa(a,e,m,i)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") -= sigmaOps[91]("e,m") * r0_L("I");

            // sigmar1_xaa_Lvo += -0.500000 <j,i||b,a>_abab r1_aa(b,m) t2_abab(e,a,j,i)
            // +                  -0.500000 <i,j||b,a>_abab r1_aa(b,m) t2_abab(e,a,i,j)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") -= sigmaOps[60]("b,e") * r1_xaa_Lvo("I,b,m");

            // sigmar1_xaa_Lvo += -0.500000 <j,i||j,i>_bbbb r1_aa(e,m)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmar1_xaa_Lvo("I,e,m") -= 0.500000 * scalar_4 * r1_xaa_Lvo("I,e,m");

            // sigmar1_xaa_Lvo += 1.000000 f_bb(i,a) r0 t2_abab(e,a,m,i)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") += sigmaOps[92]("e,m") * r0_L("I");

            // sigmar1_xaa_Lvo += 1.000000 f_bb(i,a) r2_abab(e,a,m,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xaa_Lvo("I,e,m") += F_blks_["bb_ov"]("i,a") * r2_xabab_Lvvoo("I,e,a,m,i");

            // sigmar1_xbb_Lvo += -0.500000 <j,i||a,b>_bbbb r1_bb(e,j) t2_bbbb(a,b,m,i)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") -= 0.500000 * sigmaOps[76]("j,m") * r1_xbb_Lvo("I,e,j");

            // sigmar1_xbb_Lvo += 1.000000 <i,e||a,m>_bbbb r1_bb(a,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") -= V_blks_["bbbb_vovo"]("e,i,a,m") * r1_xbb_Lvo("I,a,i");

            // sigmar1_xbb_Lvo += 0.500000 <i,e||a,b>_abab r2_abab(a,b,i,m)
            // +                  0.500000 <i,e||b,a>_abab r2_abab(b,a,i,m)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") -= V_blks_["baab_vovv"]("e,i,a,b") * r2_xabab_Lvvoo("I,a,b,i,m");

            // sigmar1_xbb_Lvo += 1.000000 f_bb(e,a) r1_bb(a,m)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") += F_blks_["bb_vv"]("e,a") * r1_xbb_Lvo("I,a,m");

            // sigmar1_xbb_Lvo += -0.500000 <j,i||a,b>_bbbb r1_bb(b,m) t2_bbbb(a,e,j,i)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") -= 0.500000 * sigmaOps[62]("e,b") * r1_xbb_Lvo("I,b,m");

            // sigmar1_xbb_Lvo += -0.500000 <j,i||j,i>_bbbb r1_bb(e,m)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmar1_xbb_Lvo("I,e,m") -= 0.500000 * scalar_4 * r1_xbb_Lvo("I,e,m");

            // sigmar1_xbb_Lvo += 1.000000 f_aa(i,a) r0 t2_abab(a,e,i,m)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") += sigmaOps[90]("e,m") * r0_L("I");

            // sigmar1_xbb_Lvo += 1.000000 <i,e||a,m>_abab r1_aa(a,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") -= V_blks_["baab_vovo"]("e,i,a,m") * r1_xaa_Lvo("I,a,i");

            // sigmar1_xbb_Lvo += -0.500000 <i,e||a,b>_bbbb r0 t2_bbbb(a,b,m,i)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") += 0.500000 * sigmaOps[65]("e,m") * r0_L("I");

            // sigmar1_xbb_Lvo += 1.000000 f_bb(i,i) r1_bb(e,m)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmar1_xbb_Lvo("I,e,m") += scalar_1 * r1_xbb_Lvo("I,e,m");

            // sigmar1_xbb_Lvo += -0.500000 <j,i||a,b>_abab r1_bb(b,m) t2_abab(a,e,j,i)
            // +                  -0.500000 <i,j||a,b>_abab r1_bb(b,m) t2_abab(a,e,i,j)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") -= sigmaOps[59]("e,b") * r1_xbb_Lvo("I,b,m");

            // sigmar1_xbb_Lvo += -0.500000 <i,e||a,b>_bbbb r2_bbbb(a,b,m,i)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") += 0.500000 * V_blks_["bbbb_vovv"]("e,i,a,b") * r2_xbbbb_Lvvoo("I,a,b,m,i");

            // sigmar1_xbb_Lvo += -0.500000 <j,i||a,m>_bbbb r0 t2_bbbb(a,e,j,i)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") -= 0.500000 * sigmaOps[80]("e,m") * r0_L("I");

            // sigmar1_xbb_Lvo += 0.250000 <j,i||a,b>_abab r1_bb(e,m) t2_abab(a,b,j,i)
            // +                  0.250000 <i,j||a,b>_abab r1_bb(e,m) t2_abab(a,b,i,j)
            // +                  0.250000 <j,i||b,a>_abab r1_bb(e,m) t2_abab(b,a,j,i)
            // +                  0.250000 <i,j||b,a>_abab r1_bb(e,m) t2_abab(b,a,i,j)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmar1_xbb_Lvo("I,e,m") += scalar_6 * r1_xbb_Lvo("I,e,m");

            // sigmar1_xbb_Lvo += 1.000000 f_aa(i,a) r2_abab(a,e,i,m)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") += F_blks_["aa_ov"]("i,a") * r2_xabab_Lvvoo("I,a,e,i,m");

            // sigmar1_xbb_Lvo += 0.500000 <i,e||a,b>_abab r0 t2_abab(a,b,i,m)
            // +                  0.500000 <i,e||b,a>_abab r0 t2_abab(b,a,i,m)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") -= sigmaOps[63]("e,m") * r0_L("I");

            // sigmar1_xbb_Lvo += -0.500000 <j,i||a,m>_abab r2_abab(a,e,j,i)
            // +                  -0.500000 <i,j||a,m>_abab r2_abab(a,e,i,j)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") -= V_blks_["abab_oovo"]("j,i,a,m") * r2_xabab_Lvvoo("I,a,e,j,i");

            // sigmar1_xbb_Lvo += -0.500000 <i,j||a,b>_abab r1_bb(e,j) t2_abab(a,b,i,m)
            // +                  -0.500000 <i,j||b,a>_abab r1_bb(e,j) t2_abab(b,a,i,m)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") -= sigmaOps[74]("j,m") * r1_xbb_Lvo("I,e,j");

            // sigmar1_xbb_Lvo += 1.000000 f_aa(i,i) r1_bb(e,m)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmar1_xbb_Lvo("I,e,m") += scalar_0 * r1_xbb_Lvo("I,e,m");

            // sigmar1_xbb_Lvo += -1.000000 f_bb(i,a) r0 t2_bbbb(a,e,m,i)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") -= sigmaOps[89]("e,m") * r0_L("I");

            // sigmar1_xbb_Lvo += 1.000000 f_bb(e,m) r0
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") += F_blks_["bb_vo"]("e,m") * r0_L("I");

            // sigmar1_xbb_Lvo += -1.000000 f_bb(i,a) r2_bbbb(a,e,m,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") -= F_blks_["bb_ov"]("i,a") * r2_xbbbb_Lvvoo("I,a,e,m,i");

            // sigmar1_xbb_Lvo += -0.500000 <j,i||j,i>_abab r1_bb(e,m)
            // +                  -0.500000 <i,j||i,j>_abab r1_bb(e,m)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmar1_xbb_Lvo("I,e,m") -= scalar_3 * r1_xbb_Lvo("I,e,m");

            // sigmar1_xbb_Lvo += 0.250000 <j,i||a,b>_bbbb r1_bb(e,m) t2_bbbb(a,b,j,i)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmar1_xbb_Lvo("I,e,m") += 0.250000 * scalar_7 * r1_xbb_Lvo("I,e,m");

            // sigmar1_xbb_Lvo += -0.500000 <j,i||a,m>_abab r0 t2_abab(a,e,j,i)
            // +                  -0.500000 <i,j||a,m>_abab r0 t2_abab(a,e,i,j)
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") -= sigmaOps[78]("e,m") * r0_L("I");

            // sigmar1_xbb_Lvo += -0.500000 <j,i||a,m>_bbbb r2_bbbb(a,e,j,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") -= 0.500000 * V_blks_["bbbb_oovo"]("j,i,a,m") * r2_xbbbb_Lvvoo("I,a,e,j,i");

            // sigmar1_xbb_Lvo += 0.250000 <j,i||a,b>_aaaa r1_bb(e,m) t2_aaaa(a,b,j,i)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmar1_xbb_Lvo("I,e,m") += 0.250000 * scalar_5 * r1_xbb_Lvo("I,e,m");

            // sigmar1_xbb_Lvo += -1.000000 f_bb(i,m) r1_bb(e,i)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmar1_xbb_Lvo("I,e,m") -= F_blks_["bb_oo"]("i,m") * r1_xbb_Lvo("I,e,i");

            // sigmar1_xbb_Lvo += -0.500000 <j,i||j,i>_aaaa r1_bb(e,m)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmar1_xbb_Lvo("I,e,m") -= 0.500000 * scalar_2 * r1_xbb_Lvo("I,e,m");

            // sigmar2_xaaaa_Lvvoo += -0.500000 <j,i||j,i>_bbbb r2_aaaa(e,f,m,n)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") -= 0.500000 * scalar_4 * r2_xaaaa_Lvvoo("I,e,f,m,n");

            // sigmar2_xaaaa_Lvvoo += -0.500000 <j,i||b,a>_abab r0 t2_abab(e,a,j,i) t2_aaaa(b,f,m,n)
            // +                      -0.500000 <i,j||b,a>_abab r0 t2_abab(e,a,i,j) t2_aaaa(b,f,m,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") -= sigmaOps[107]("f,e,m,n") * r0_L("I");

            // sigmar2_xaaaa_Lvvoo += 0.500000 <j,i||a,b>_abab r0 t2_aaaa(a,e,m,n) t2_abab(f,b,j,i)
            // +                      0.500000 <i,j||a,b>_abab r0 t2_aaaa(a,e,m,n) t2_abab(f,b,i,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += sigmaOps[107]("e,f,m,n") * r0_L("I");

            // sigmar2_xaaaa_Lvvoo += 0.500000 <j,i||m,n>_aaaa r2_aaaa(e,f,j,i)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += 0.500000 * V_blks_["aaaa_oooo"]("j,i,m,n") * r2_xaaaa_Lvvoo("I,e,f,j,i");

            // sigmar2_xaaaa_Lvvoo += 0.250000 <j,i||a,b>_aaaa r0 t2_aaaa(a,b,m,n) t2_aaaa(e,f,j,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += 0.250000 * sigmaOps[105]("e,f,m,n") * r0_L("I");

            // sigmar2_xaaaa_Lvvoo += -0.500000 <j,i||j,i>_aaaa r2_aaaa(e,f,m,n)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") -= 0.500000 * scalar_2 * r2_xaaaa_Lvvoo("I,e,f,m,n");

            // sigmar2_xaaaa_Lvvoo += 0.500000 <e,f||a,b>_aaaa r0 t2_aaaa(a,b,m,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += 0.500000 * sigmaOps[1]("e,f,m,n") * r0_L("I");

            // sigmar2_xaaaa_Lvvoo += 0.500000 <j,i||m,n>_aaaa r0 t2_aaaa(e,f,j,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += 0.500000 * sigmaOps[51]("e,f,m,n") * r0_L("I");

            // sigmar2_xaaaa_Lvvoo += 0.250000 <j,i||a,b>_abab r2_aaaa(e,f,m,n) t2_abab(a,b,j,i)
            // +                      0.250000 <i,j||a,b>_abab r2_aaaa(e,f,m,n) t2_abab(a,b,i,j)
            // +                      0.250000 <j,i||b,a>_abab r2_aaaa(e,f,m,n) t2_abab(b,a,j,i)
            // +                      0.250000 <i,j||b,a>_abab r2_aaaa(e,f,m,n) t2_abab(b,a,i,j)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += scalar_6 * r2_xaaaa_Lvvoo("I,e,f,m,n");

            // sigmar2_xaaaa_Lvvoo += 0.250000 <j,i||a,b>_aaaa r2_aaaa(e,f,m,n) t2_aaaa(a,b,j,i)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += 0.250000 * scalar_5 * r2_xaaaa_Lvvoo("I,e,f,m,n");

            // sigmar2_xaaaa_Lvvoo += -0.500000 <j,i||j,i>_abab r2_aaaa(e,f,m,n)
            // +                      -0.500000 <i,j||i,j>_abab r2_aaaa(e,f,m,n)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") -= scalar_3 * r2_xaaaa_Lvvoo("I,e,f,m,n");

            // sigmar2_xaaaa_Lvvoo += 1.000000 <e,f||m,n>_aaaa r0
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += V_blks_["aaaa_vvoo"]("e,f,m,n") * r0_L("I");

            // sigmar2_xaaaa_Lvvoo += 0.500000 <e,f||a,b>_aaaa r2_aaaa(a,b,m,n)
            // flops: o2v4L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += 0.500000 * V_blks_["aaaa_vvvv"]("e,f,a,b") * r2_xaaaa_Lvvoo("I,a,b,m,n");

            // sigmar2_xaaaa_Lvvoo += 0.250000 <j,i||a,b>_aaaa r2_aaaa(a,b,m,n) t2_aaaa(e,f,j,i)
            // flops: o4v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o4v0L1: 1,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += 0.250000 * V_blks_["aaaa_oovv"]("j,i,a,b") * r2_xaaaa_Lvvoo("I,a,b,m,n") * t2_aaaa_vvoo("e,f,j,i");

            // sigmar2_xaaaa_Lvvoo += -0.500000 <j,i||a,b>_aaaa r0 t2_aaaa(a,e,j,i) t2_aaaa(b,f,m,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[109]("e,f,m,n") * r0_L("I");

            // sigmar2_xaaaa_Lvvoo += 0.250000 <j,i||a,b>_aaaa r2_aaaa(e,f,j,i) t2_aaaa(a,b,m,n)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += 0.250000 * sigmaOps[42]("m,n,j,i") * r2_xaaaa_Lvvoo("I,e,f,j,i");

            // sigmar2_xaaaa_Lvvoo += -0.500000 <j,i||a,b>_aaaa r0 t2_aaaa(a,e,m,n) t2_aaaa(b,f,j,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[67]("f,e,m,n") * r0_L("I");

            // sigmar2_xaaaa_Lvvoo += 1.000000 f_aa(i,i) r2_aaaa(e,f,m,n)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += scalar_0 * r2_xaaaa_Lvvoo("I,e,f,m,n");

            // sigmar2_xaaaa_Lvvoo += 1.000000 f_bb(i,i) r2_aaaa(e,f,m,n)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += scalar_1 * r2_xaaaa_Lvvoo("I,e,f,m,n");

            // sigmar2_xaaaa_Lvvoo += 0.250000 <j,i||a,b>_bbbb r2_aaaa(e,f,m,n) t2_bbbb(a,b,j,i)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += 0.250000 * scalar_7 * r2_xaaaa_Lvvoo("I,e,f,m,n");

            // sigmar2_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += tempPerm_xaaaa_Lvvoo("I,e,f,m,n");
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") -= tempPerm_xaaaa_Lvvoo("I,e,f,n,m");

            // tempPerm_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xaaaa_Lvvoo = TArrayD(world_, tempPerm_xaaaa_Lvvoo.trange());
            tempPerm_xaaaa_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xaaaa_Lvvoo += 1.000000 P(e,f) f_aa(e,a) r0 t2_aaaa(a,f,m,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") = sigmaOps[72]("e,f,m,n") * r0_L("I");

            // tempPerm_xaaaa_Lvvoo += -0.500000 P(e,f) <j,i||a,b>_aaaa r2_aaaa(b,f,m,n) t2_aaaa(a,e,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[61]("e,b") * r2_xaaaa_Lvvoo("I,b,f,m,n");

            // tempPerm_xaaaa_Lvvoo += 0.500000 P(e,f) <i,e||a,b>_aaaa r1_aa(f,i) t2_aaaa(a,b,m,n)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[26]("e,i,m,n") * r1_xaa_Lvo("I,f,i");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(e,f) f_aa(i,a) r1_aa(f,i) t2_aaaa(a,e,m,n)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += sigmaOps[88]("e,i,m,n") * r1_xaa_Lvo("I,f,i");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(e,f) <i,e||m,n>_aaaa r1_aa(f,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= V_blks_["aaaa_vooo"]("e,i,m,n") * r1_xaa_Lvo("I,f,i");

            // tempPerm_xaaaa_Lvvoo += -0.500000 P(e,f) <j,i||b,a>_abab r2_aaaa(b,f,m,n) t2_abab(e,a,j,i)
            // +                      -0.500000 P(e,f) <i,j||b,a>_abab r2_aaaa(b,f,m,n) t2_abab(e,a,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= sigmaOps[60]("b,e") * r2_xaaaa_Lvvoo("I,b,f,m,n");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(e,f) f_aa(e,a) r2_aaaa(a,f,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += F_blks_["aa_vv"]("e,a") * r2_xaaaa_Lvvoo("I,a,f,m,n");

            // sigmar2_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += tempPerm_xaaaa_Lvvoo("I,e,f,m,n");
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") -= tempPerm_xaaaa_Lvvoo("I,f,e,m,n");

            // tempPerm_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xaaaa_Lvvoo = TArrayD(world_, tempPerm_xaaaa_Lvvoo.trange());
            tempPerm_xaaaa_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) <j,i||a,b>_aaaa r0 t2_aaaa(a,e,n,i) t2_aaaa(b,f,m,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") = sigmaOps[98]("f,e,m,n") * r0_L("I");

            // tempPerm_xaaaa_Lvvoo += -1.000000 P(m,n) <j,i||a,n>_aaaa r1_aa(a,j) t2_aaaa(e,f,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1, o3v1L1: 1 | mem: o2v2L1: 2, o2v0L1: 1,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= r1_xaa_Lvo("I,a,j") * V_blks_["aaaa_oovo"]("j,i,a,n") * t2_aaaa_vvoo("e,f,m,i");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) <e,f||a,n>_aaaa r1_aa(a,m)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += V_blks_["aaaa_vvvo"]("e,f,a,n") * r1_xaa_Lvo("I,a,m");

            // tempPerm_xaaaa_Lvvoo += -0.500000 P(m,n) <j,i||a,b>_abab r2_aaaa(e,f,m,j) t2_abab(a,b,n,i)
            // +                      -0.500000 P(m,n) <j,i||b,a>_abab r2_aaaa(e,f,m,j) t2_abab(b,a,n,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= sigmaOps[73]("j,n") * r2_xaaaa_Lvvoo("I,e,f,m,j");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) <j,i||a,b>_bbbb r0 t2_abab(e,a,n,i) t2_abab(f,b,m,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += sigmaOps[97]("f,e,m,n") * r0_L("I");

            // tempPerm_xaaaa_Lvvoo += -1.000000 P(m,n) <i,j||n,a>_abab r1_bb(a,j) t2_aaaa(e,f,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1, o3v1L1: 1 | mem: o2v2L1: 2, o2v0L1: 1,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += r1_xbb_Lvo("I,a,j") * V_blks_["abba_oovo"]("i,j,a,n") * t2_aaaa_vvoo("e,f,m,i");

            // tempPerm_xaaaa_Lvvoo += -1.000000 P(m,n) f_aa(i,n) r0 t2_aaaa(e,f,m,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= sigmaOps[86]("e,f,m,n") * r0_L("I");

            // tempPerm_xaaaa_Lvvoo += -0.500000 P(m,n) <j,i||a,b>_aaaa r2_aaaa(e,f,m,j) t2_aaaa(a,b,n,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[75]("j,n") * r2_xaaaa_Lvvoo("I,e,f,m,j");

            // tempPerm_xaaaa_Lvvoo += -0.500000 P(m,n) <j,i||a,b>_abab r0 t2_abab(a,b,n,i) t2_aaaa(e,f,m,j)
            // +                      -0.500000 P(m,n) <j,i||b,a>_abab r0 t2_abab(b,a,n,i) t2_aaaa(e,f,m,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= sigmaOps[117]("e,f,m,n") * r0_L("I");

            // tempPerm_xaaaa_Lvvoo += -1.000000 P(m,n) f_aa(i,n) r2_aaaa(e,f,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= F_blks_["aa_oo"]("i,n") * r2_xaaaa_Lvvoo("I,e,f,m,i");

            // tempPerm_xaaaa_Lvvoo += 0.500000 P(m,n) <j,i||a,n>_aaaa r1_aa(a,m) t2_aaaa(e,f,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += 0.500000 * sigmaOps[24]("a,e,f,n") * r1_xaa_Lvo("I,a,m");

            // tempPerm_xaaaa_Lvvoo += -0.500000 P(m,n) <j,i||a,b>_aaaa r0 t2_aaaa(a,b,n,i) t2_aaaa(e,f,m,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[118]("e,f,n,m") * r0_L("I");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) f_aa(i,a) r1_aa(a,m) t2_aaaa(e,f,n,i)
            // flops: o3v2L1: 1, o2v2L1: 1, o2v1L1: 1 | mem: o2v2L1: 2, o2v0L1: 1,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += r1_xaa_Lvo("I,a,m") * F_blks_["aa_ov"]("i,a") * t2_aaaa_vvoo("e,f,n,i");

            // sigmar2_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += tempPerm_xaaaa_Lvvoo("I,e,f,m,n");
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") -= tempPerm_xaaaa_Lvvoo("I,e,f,n,m");

            // tempPerm_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xaaaa_Lvvoo = TArrayD(world_, tempPerm_xaaaa_Lvvoo.trange());
            tempPerm_xaaaa_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) P(e,f) <j,i||b,a>_abab r2_aaaa(b,f,m,j) t2_abab(e,a,n,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") = sigmaOps[14]("e,b,n,j") * r2_xaaaa_Lvvoo("I,b,f,m,j");

            // tempPerm_xaaaa_Lvvoo += -1.000000 P(m,n) P(e,f) <j,i||a,n>_aaaa r1_aa(f,j) t2_aaaa(a,e,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= sigmaOps[47]("e,j,n,m") * r1_xaa_Lvo("I,f,j");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) P(e,f) <i,j||a,b>_abab r2_abab(f,b,m,j) t2_aaaa(a,e,n,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += sigmaOps[20]("e,b,n,j") * r2_xabab_Lvvoo("I,f,b,m,j");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) P(e,f) <i,e||a,n>_aaaa r0 t2_aaaa(a,f,m,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= sigmaOps[33]("f,e,m,n") * r0_L("I");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) P(e,f) <i,e||a,n>_aaaa r2_aaaa(a,f,m,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= V_blks_["aaaa_vovo"]("e,i,a,n") * r2_xaaaa_Lvvoo("I,a,f,m,i");

            // tempPerm_xaaaa_Lvvoo += 0.500000 P(m,n) P(e,f) <j,i||n,a>_abab r1_aa(f,m) t2_abab(e,a,j,i)
            // +                      0.500000 P(m,n) P(e,f) <i,j||n,a>_abab r1_aa(f,m) t2_abab(e,a,i,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= sigmaOps[79]("e,n") * r1_xaa_Lvo("I,f,m");

            // tempPerm_xaaaa_Lvvoo += -1.000000 P(m,n) P(e,f) <e,i||n,a>_abab r0 t2_abab(f,a,m,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += sigmaOps[31]("e,f,n,m") * r0_L("I");

            // tempPerm_xaaaa_Lvvoo += -1.000000 P(m,n) P(e,f) <j,i||n,a>_abab r1_aa(f,j) t2_abab(e,a,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += sigmaOps[50]("e,j,n,m") * r1_xaa_Lvo("I,f,j");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) P(e,f) <j,i||a,b>_bbbb r2_abab(f,b,m,j) t2_abab(e,a,n,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += sigmaOps[16]("e,b,n,j") * r2_xabab_Lvvoo("I,f,b,m,j");

            // tempPerm_xaaaa_Lvvoo += -1.000000 P(m,n) P(e,f) <e,i||n,a>_abab r2_abab(f,a,m,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += V_blks_["abba_vovo"]("e,i,a,n") * r2_xabab_Lvvoo("I,f,a,m,i");

            // tempPerm_xaaaa_Lvvoo += -1.000000 P(m,n) P(e,f) f_aa(e,n) r1_aa(f,m)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= F_blks_["aa_vo"]("e,n") * r1_xaa_Lvo("I,f,m");

            // tempPerm_xaaaa_Lvvoo += -1.000000 P(m,n) P(e,f) <i,e||a,b>_aaaa r1_aa(b,m) t2_aaaa(a,f,n,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += sigmaOps[7]("e,b,f,n") * r1_xaa_Lvo("I,b,m");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) P(e,f) <j,i||a,b>_aaaa r2_aaaa(b,f,m,j) t2_aaaa(a,e,n,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += sigmaOps[15]("b,e,j,n") * r2_xaaaa_Lvvoo("I,b,f,m,j");

            // tempPerm_xaaaa_Lvvoo += -0.500000 P(m,n) P(e,f) <e,i||a,b>_abab r1_aa(f,m) t2_abab(a,b,n,i)
            // +                      -0.500000 P(m,n) P(e,f) <e,i||b,a>_abab r1_aa(f,m) t2_abab(b,a,n,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= sigmaOps[64]("e,n") * r1_xaa_Lvo("I,f,m");

            // tempPerm_xaaaa_Lvvoo += 0.500000 P(m,n) P(e,f) <i,e||a,b>_aaaa r1_aa(f,m) t2_aaaa(a,b,n,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[66]("e,n") * r1_xaa_Lvo("I,f,m");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) P(e,f) <e,i||b,a>_abab r1_aa(b,m) t2_abab(f,a,n,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += sigmaOps[5]("e,b,f,n") * r1_xaa_Lvo("I,b,m");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) P(e,f) f_aa(i,a) r1_aa(f,m) t2_aaaa(a,e,n,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += sigmaOps[91]("e,n") * r1_xaa_Lvo("I,f,m");

            // tempPerm_xaaaa_Lvvoo += 0.500000 P(m,n) P(e,f) <j,i||a,n>_aaaa r1_aa(f,m) t2_aaaa(a,e,j,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += 0.500000 * sigmaOps[77]("e,n") * r1_xaa_Lvo("I,f,m");

            // tempPerm_xaaaa_Lvvoo += -1.000000 P(m,n) P(e,f) f_bb(i,a) r1_aa(f,m) t2_abab(e,a,n,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= sigmaOps[92]("e,n") * r1_xaa_Lvo("I,f,m");

            // sigmar2_xabab_Lvvoo += 0.250000 <j,i||a,b>_abab r2_abab(e,f,j,i) t2_abab(a,b,m,n)
            // +                      0.250000 <j,i||b,a>_abab r2_abab(e,f,j,i) t2_abab(b,a,m,n)
            // +                      0.250000 <i,j||a,b>_abab r2_abab(e,f,i,j) t2_abab(a,b,m,n)
            // +                      0.250000 <i,j||b,a>_abab r2_abab(e,f,i,j) t2_abab(b,a,m,n)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[54]("m,j,n,i") * r2_xabab_Lvvoo("I,e,f,j,i");

            // sigmar2_xabab_Lvvoo += 1.000000 <j,i||a,b>_bbbb r0 t2_abab(e,a,m,i) t2_bbbb(b,f,n,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[99]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += -1.000000 <e,i||a,b>_abab r1_bb(b,n) t2_abab(a,f,m,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[2]("e,f,b,m") * r1_xbb_Lvo("I,b,n");

            // sigmar2_xabab_Lvvoo += 0.500000 <e,i||a,b>_abab r1_bb(f,n) t2_abab(a,b,m,i)
            // +                      0.500000 <e,i||b,a>_abab r1_bb(f,n) t2_abab(b,a,m,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[64]("e,m") * r1_xbb_Lvo("I,f,n");

            // sigmar2_xabab_Lvvoo += -1.000000 <i,f||a,n>_abab r0 t2_aaaa(a,e,m,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[32]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||a,b>_bbbb r0 t2_bbbb(a,b,n,i) t2_abab(e,f,m,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[122]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += 1.000000 <j,i||a,b>_aaaa r0 t2_aaaa(a,e,m,i) t2_abab(b,f,j,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[100]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += -1.000000 f_bb(i,n) r2_abab(e,f,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= F_blks_["bb_oo"]("i,n") * r2_xabab_Lvvoo("I,e,f,m,i");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||a,b>_bbbb r2_abab(e,b,m,n) t2_bbbb(a,f,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[62]("f,b") * r2_xabab_Lvvoo("I,e,b,m,n");

            // sigmar2_xabab_Lvvoo += 1.000000 <i,j||a,b>_abab r0 t2_aaaa(a,e,m,i) t2_bbbb(b,f,n,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[95]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += 1.000000 f_aa(e,a) r2_abab(a,f,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += F_blks_["aa_vv"]("e,a") * r2_xabab_Lvvoo("I,a,f,m,n");

            // sigmar2_xabab_Lvvoo += 1.000000 <i,j||b,a>_abab r2_abab(b,f,m,j) t2_abab(e,a,i,n)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[34]("b,e,j,n") * r2_xabab_Lvvoo("I,b,f,m,j");

            // sigmar2_xabab_Lvvoo += 1.000000 <i,f||a,n>_bbbb r2_abab(e,a,m,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= V_blks_["bbbb_vovo"]("f,i,a,n") * r2_xabab_Lvvoo("I,e,a,m,i");

            // sigmar2_xabab_Lvvoo += 1.000000 <i,j||m,a>_abab r1_bb(f,j) t2_abab(e,a,i,n)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[58]("e,m,j,n") * r1_xbb_Lvo("I,f,j");

            // sigmar2_xabab_Lvvoo += 0.500000 <e,f||a,b>_abab r2_abab(a,b,m,n)
            // +                      0.500000 <e,f||b,a>_abab r2_abab(b,a,m,n)
            // flops: o2v4L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += V_blks_["abab_vvvv"]("e,f,a,b") * r2_xabab_Lvvoo("I,a,b,m,n");

            // sigmar2_xabab_Lvvoo += -1.000000 <j,i||a,n>_bbbb r1_bb(a,j) t2_abab(e,f,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1, o3v1L1: 1 | mem: o2v2L1: 2, o2v0L1: 1,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= V_blks_["bbbb_oovo"]("j,i,a,n") * r1_xbb_Lvo("I,a,j") * t2_abab_vvoo("e,f,m,i");

            // sigmar2_xabab_Lvvoo += 1.000000 <j,i||a,b>_abab r2_abab(e,b,j,n) t2_abab(a,f,m,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[36]("b,f,j,m") * r2_xabab_Lvvoo("I,e,b,j,n");

            // sigmar2_xabab_Lvvoo += -1.000000 <e,i||a,n>_abab r0 t2_abab(a,f,m,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[41]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += 0.250000 <j,i||a,b>_abab r2_abab(a,b,m,n) t2_abab(e,f,j,i)
            // +                      0.250000 <i,j||a,b>_abab r2_abab(a,b,m,n) t2_abab(e,f,i,j)
            // +                      0.250000 <j,i||b,a>_abab r2_abab(b,a,m,n) t2_abab(e,f,j,i)
            // +                      0.250000 <i,j||b,a>_abab r2_abab(b,a,m,n) t2_abab(e,f,i,j)
            // flops: o4v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o4v0L1: 1,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += V_blks_["abab_oovv"]("j,i,a,b") * r2_xabab_Lvvoo("I,a,b,m,n") * t2_abab_vvoo("e,f,j,i");

            // sigmar2_xabab_Lvvoo += -1.000000 f_bb(i,a) r1_aa(e,m) t2_bbbb(a,f,n,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[89]("f,n") * r1_xaa_Lvo("I,e,m");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||a,b>_abab r2_abab(e,f,j,n) t2_abab(a,b,m,i)
            // +                      -0.500000 <j,i||b,a>_abab r2_abab(e,f,j,n) t2_abab(b,a,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[73]("j,m") * r2_xabab_Lvvoo("I,e,f,j,n");

            // sigmar2_xabab_Lvvoo += 0.500000 <j,i||a,b>_bbbb r0 t2_abab(e,a,m,n) t2_bbbb(b,f,j,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += 0.500000 * sigmaOps[116]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += 1.000000 <j,i||b,a>_abab r2_abab(b,f,j,n) t2_abab(e,a,m,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[14]("e,b,m,j") * r2_xabab_Lvvoo("I,b,f,j,n");

            // sigmar2_xabab_Lvvoo += -1.000000 <e,i||m,a>_abab r0 t2_bbbb(a,f,n,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[37]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += 1.000000 <i,j||a,b>_abab r2_bbbb(b,f,n,j) t2_aaaa(a,e,m,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[20]("e,b,m,j") * r2_xbbbb_Lvvoo("I,b,f,n,j");

            // sigmar2_xabab_Lvvoo += 1.000000 f_bb(f,a) r2_abab(e,a,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += F_blks_["bb_vv"]("f,a") * r2_xabab_Lvvoo("I,e,a,m,n");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||a,b>_aaaa r2_abab(b,f,m,n) t2_aaaa(a,e,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[61]("e,b") * r2_xabab_Lvvoo("I,b,f,m,n");

            // sigmar2_xabab_Lvvoo += 1.000000 <j,i||b,a>_abab r2_aaaa(b,e,m,j) t2_bbbb(a,f,n,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[13]("b,f,j,n") * r2_xaaaa_Lvvoo("I,b,e,m,j");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||b,a>_abab r2_abab(b,f,m,n) t2_abab(e,a,j,i)
            // +                      -0.500000 <i,j||b,a>_abab r2_abab(b,f,m,n) t2_abab(e,a,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[60]("b,e") * r2_xabab_Lvvoo("I,b,f,m,n");

            // sigmar2_xabab_Lvvoo += 0.250000 <j,i||a,b>_bbbb r2_abab(e,f,m,n) t2_bbbb(a,b,j,i)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += 0.250000 * scalar_7 * r2_xabab_Lvvoo("I,e,f,m,n");

            // sigmar2_xabab_Lvvoo += -1.000000 <j,i||a,m>_aaaa r1_aa(a,j) t2_abab(e,f,i,n)
            // flops: o3v2L1: 1, o2v2L1: 1, o3v1L1: 1 | mem: o2v2L1: 2, o2v0L1: 1,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= V_blks_["aaaa_oovo"]("j,i,a,m") * r1_xaa_Lvo("I,a,j") * t2_abab_vvoo("e,f,i,n");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||a,n>_abab r1_aa(e,m) t2_abab(a,f,j,i)
            // +                      -0.500000 <i,j||a,n>_abab r1_aa(e,m) t2_abab(a,f,i,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[78]("f,n") * r1_xaa_Lvo("I,e,m");

            // sigmar2_xabab_Lvvoo += 0.250000 <j,i||a,b>_abab r2_abab(e,f,m,n) t2_abab(a,b,j,i)
            // +                      0.250000 <i,j||a,b>_abab r2_abab(e,f,m,n) t2_abab(a,b,i,j)
            // +                      0.250000 <j,i||b,a>_abab r2_abab(e,f,m,n) t2_abab(b,a,j,i)
            // +                      0.250000 <i,j||b,a>_abab r2_abab(e,f,m,n) t2_abab(b,a,i,j)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += scalar_6 * r2_xabab_Lvvoo("I,e,f,m,n");

            // sigmar2_xabab_Lvvoo += 0.500000 <j,i||m,a>_abab r1_bb(a,n) t2_abab(e,f,j,i)
            // +                      0.500000 <i,j||m,a>_abab r1_bb(a,n) t2_abab(e,f,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[40]("e,a,f,m") * r1_xbb_Lvo("I,a,n");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||a,b>_aaaa r0 t2_aaaa(a,b,m,i) t2_abab(e,f,j,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[120]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += -1.000000 <e,i||b,a>_abab r1_aa(b,m) t2_bbbb(a,f,n,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[10]("e,b,f,n") * r1_xaa_Lvo("I,b,m");

            // sigmar2_xabab_Lvvoo += 1.000000 <j,i||a,n>_bbbb r1_bb(f,j) t2_abab(e,a,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[55]("e,m,j,n") * r1_xbb_Lvo("I,f,j");

            // sigmar2_xabab_Lvvoo += -1.000000 <i,f||m,a>_abab r2_abab(e,a,i,n)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= V_blks_["baba_vovo"]("f,i,a,m") * r2_xabab_Lvvoo("I,e,a,i,n");

            // sigmar2_xabab_Lvvoo += -1.000000 <i,f||b,a>_abab r1_aa(b,m) t2_abab(e,a,i,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[9]("b,e,f,n") * r1_xaa_Lvo("I,b,m");

            // sigmar2_xabab_Lvvoo += 1.000000 <j,i||b,a>_abab r0 t2_abab(e,a,m,i) t2_abab(b,f,j,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[96]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||a,b>_abab r0 t2_abab(a,b,m,i) t2_abab(e,f,j,n)
            // +                      -0.500000 <j,i||b,a>_abab r0 t2_abab(b,a,m,i) t2_abab(e,f,j,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[121]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += -1.000000 f_bb(i,n) r0 t2_abab(e,f,m,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[87]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += -1.000000 f_aa(i,m) r2_abab(e,f,i,n)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= F_blks_["aa_oo"]("i,m") * r2_xabab_Lvvoo("I,e,f,i,n");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||a,b>_abab r2_abab(e,b,m,n) t2_abab(a,f,j,i)
            // +                      -0.500000 <i,j||a,b>_abab r2_abab(e,b,m,n) t2_abab(a,f,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[59]("f,b") * r2_xabab_Lvvoo("I,e,b,m,n");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||j,i>_aaaa r2_abab(e,f,m,n)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= 0.500000 * scalar_2 * r2_xabab_Lvvoo("I,e,f,m,n");

            // sigmar2_xabab_Lvvoo += 1.000000 <i,j||a,b>_abab r2_abab(e,b,m,j) t2_abab(a,f,i,n)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[21]("b,f,j,n") * r2_xabab_Lvvoo("I,e,b,m,j");

            // sigmar2_xabab_Lvvoo += 1.000000 <j,i||a,b>_aaaa r2_abab(b,f,j,n) t2_aaaa(a,e,m,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[15]("b,e,j,m") * r2_xabab_Lvvoo("I,b,f,j,n");

            // sigmar2_xabab_Lvvoo += 1.000000 f_aa(e,a) r0 t2_abab(a,f,m,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[68]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += 1.000000 <i,j||b,a>_abab r0 t2_abab(e,a,i,n) t2_abab(b,f,m,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[19]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += 1.000000 f_bb(f,a) r0 t2_abab(e,a,m,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[70]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += -1.000000 <i,f||m,n>_abab r1_aa(e,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += V_blks_["baab_vooo"]("f,i,m,n") * r1_xaa_Lvo("I,e,i");

            // sigmar2_xabab_Lvvoo += -1.000000 <i,f||a,n>_abab r2_aaaa(a,e,m,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += V_blks_["baab_vovo"]("f,i,a,n") * r2_xaaaa_Lvvoo("I,a,e,m,i");

            // sigmar2_xabab_Lvvoo += -1.000000 f_aa(i,m) r0 t2_abab(e,f,i,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[81]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += 1.000000 <j,i||a,n>_abab r1_aa(e,j) t2_abab(a,f,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[53]("f,m,j,n") * r1_xaa_Lvo("I,e,j");

            // sigmar2_xabab_Lvvoo += 1.000000 <i,f||a,b>_bbbb r1_bb(b,n) t2_abab(e,a,m,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[12]("e,f,b,m") * r1_xbb_Lvo("I,b,n");

            // sigmar2_xabab_Lvvoo += 1.000000 <j,i||a,b>_aaaa r2_aaaa(b,e,m,j) t2_abab(a,f,i,n)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[18]("b,f,j,n") * r2_xaaaa_Lvvoo("I,b,e,m,j");

            // sigmar2_xabab_Lvvoo += -1.000000 f_aa(i,a) r1_aa(a,m) t2_abab(e,f,i,n)
            // flops: o3v2L1: 1, o2v2L1: 1, o2v1L1: 1 | mem: o2v2L1: 2, o2v0L1: 1,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= F_blks_["aa_ov"]("i,a") * r1_xaa_Lvo("I,a,m") * t2_abab_vvoo("e,f,i,n");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||a,n>_bbbb r1_aa(e,m) t2_bbbb(a,f,j,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[80]("f,n") * r1_xaa_Lvo("I,e,m");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||a,b>_aaaa r0 t2_aaaa(a,e,j,i) t2_abab(b,f,m,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[110]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += 1.000000 <i,f||a,n>_bbbb r0 t2_abab(e,a,m,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[38]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += 1.000000 <i,e||a,m>_aaaa r2_abab(a,f,i,n)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= V_blks_["aaaa_vovo"]("e,i,a,m") * r2_xabab_Lvvoo("I,a,f,i,n");

            // sigmar2_xabab_Lvvoo += 1.000000 <j,i||a,m>_aaaa r1_aa(e,j) t2_abab(a,f,i,n)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[57]("f,j,m,n") * r1_xaa_Lvo("I,e,j");

            // sigmar2_xabab_Lvvoo += -0.500000 <i,e||a,b>_aaaa r1_bb(f,n) t2_aaaa(a,b,m,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += 0.500000 * sigmaOps[66]("e,m") * r1_xbb_Lvo("I,f,n");

            // sigmar2_xabab_Lvvoo += 0.500000 <j,i||a,n>_abab r1_aa(a,m) t2_abab(e,f,j,i)
            // +                      0.500000 <i,j||a,n>_abab r1_aa(a,m) t2_abab(e,f,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[25]("e,a,f,n") * r1_xaa_Lvo("I,a,m");

            // sigmar2_xabab_Lvvoo += -1.000000 <i,f||a,b>_abab r1_bb(b,n) t2_aaaa(a,e,m,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[6]("e,f,b,m") * r1_xbb_Lvo("I,b,n");

            // sigmar2_xabab_Lvvoo += 0.500000 <j,i||m,n>_abab r0 t2_abab(e,f,j,i)
            // +                      0.500000 <i,j||m,n>_abab r0 t2_abab(e,f,i,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[46]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += 1.000000 <j,i||a,b>_bbbb r2_bbbb(b,f,n,j) t2_abab(e,a,m,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[16]("e,b,m,j") * r2_xbbbb_Lvvoo("I,b,f,n,j");

            // sigmar2_xabab_Lvvoo += -0.500000 <e,i||a,b>_abab r1_bb(f,i) t2_abab(a,b,m,n)
            // +                      -0.500000 <e,i||b,a>_abab r1_bb(f,i) t2_abab(b,a,m,n)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[35]("e,m,i,n") * r1_xbb_Lvo("I,f,i");

            // sigmar2_xabab_Lvvoo += 0.500000 <i,f||a,b>_abab r1_aa(e,m) t2_abab(a,b,i,n)
            // +                      0.500000 <i,f||b,a>_abab r1_aa(e,m) t2_abab(b,a,i,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[63]("f,n") * r1_xaa_Lvo("I,e,m");

            // sigmar2_xabab_Lvvoo += 1.000000 f_bb(i,i) r2_abab(e,f,m,n)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += scalar_1 * r2_xabab_Lvvoo("I,e,f,m,n");

            // sigmar2_xabab_Lvvoo += 1.000000 <e,f||m,a>_abab r1_bb(a,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= V_blks_["abba_vvvo"]("e,f,a,m") * r1_xbb_Lvo("I,a,n");

            // sigmar2_xabab_Lvvoo += 1.000000 f_bb(i,a) r1_bb(f,n) t2_abab(e,a,m,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[92]("e,m") * r1_xbb_Lvo("I,f,n");

            // sigmar2_xabab_Lvvoo += 1.000000 <i,j||a,n>_abab r1_bb(f,j) t2_aaaa(a,e,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[56]("e,m,j,n") * r1_xbb_Lvo("I,f,j");

            // sigmar2_xabab_Lvvoo += -1.000000 f_bb(i,a) r1_bb(f,i) t2_abab(e,a,m,n)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[85]("e,m,n,i") * r1_xbb_Lvo("I,f,i");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||b,a>_abab r0 t2_abab(e,a,m,n) t2_abab(b,f,j,i)
            // +                      -0.500000 <i,j||b,a>_abab r0 t2_abab(e,a,m,n) t2_abab(b,f,i,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[112]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += -0.500000 <i,j||a,b>_abab r2_abab(e,f,m,j) t2_abab(a,b,i,n)
            // +                      -0.500000 <i,j||b,a>_abab r2_abab(e,f,m,j) t2_abab(b,a,i,n)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[74]("j,n") * r2_xabab_Lvvoo("I,e,f,m,j");

            // sigmar2_xabab_Lvvoo += 1.000000 <j,i||a,b>_bbbb r2_abab(e,b,m,j) t2_bbbb(a,f,n,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[17]("f,b,n,j") * r2_xabab_Lvvoo("I,e,b,m,j");

            // sigmar2_xabab_Lvvoo += -1.000000 <i,f||m,a>_abab r0 t2_abab(e,a,i,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[39]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += -0.500000 <i,j||a,b>_abab r0 t2_abab(a,b,i,n) t2_abab(e,f,m,j)
            // +                      -0.500000 <i,j||b,a>_abab r0 t2_abab(b,a,i,n) t2_abab(e,f,m,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[119]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||j,i>_abab r2_abab(e,f,m,n)
            // +                      -0.500000 <i,j||i,j>_abab r2_abab(e,f,m,n)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= scalar_3 * r2_xabab_Lvvoo("I,e,f,m,n");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||j,i>_bbbb r2_abab(e,f,m,n)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= 0.500000 * scalar_4 * r2_xabab_Lvvoo("I,e,f,m,n");

            // sigmar2_xabab_Lvvoo += -1.000000 f_aa(i,a) r1_aa(e,i) t2_abab(a,f,m,n)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[82]("f,i,m,n") * r1_xaa_Lvo("I,e,i");

            // sigmar2_xabab_Lvvoo += 1.000000 <i,e||a,m>_aaaa r0 t2_abab(a,f,i,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[23]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += -1.000000 <j,i||a,n>_abab r1_aa(a,j) t2_abab(e,f,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1, o3v1L1: 1 | mem: o2v2L1: 2, o2v0L1: 1,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= V_blks_["abab_oovo"]("j,i,a,n") * r1_xaa_Lvo("I,a,j") * t2_abab_vvoo("e,f,m,i");

            // sigmar2_xabab_Lvvoo += -1.000000 <e,i||a,n>_abab r2_abab(a,f,m,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= V_blks_["abab_vovo"]("e,i,a,n") * r2_xabab_Lvvoo("I,a,f,m,i");

            // sigmar2_xabab_Lvvoo += 0.250000 <j,i||a,b>_abab r0 t2_abab(a,b,m,n) t2_abab(e,f,j,i)
            // +                      0.250000 <i,j||a,b>_abab r0 t2_abab(a,b,m,n) t2_abab(e,f,i,j)
            // +                      0.250000 <j,i||b,a>_abab r0 t2_abab(b,a,m,n) t2_abab(e,f,j,i)
            // +                      0.250000 <i,j||b,a>_abab r0 t2_abab(b,a,m,n) t2_abab(e,f,i,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[44]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += 1.000000 f_aa(i,i) r2_abab(e,f,m,n)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += scalar_0 * r2_xabab_Lvvoo("I,e,f,m,n");

            // sigmar2_xabab_Lvvoo += 1.000000 f_aa(i,a) r1_aa(e,m) t2_abab(a,f,i,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[90]("f,n") * r1_xaa_Lvo("I,e,m");

            // sigmar2_xabab_Lvvoo += 1.000000 <i,e||a,b>_aaaa r1_aa(b,m) t2_abab(a,f,i,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[0]("e,b,f,n") * r1_xaa_Lvo("I,b,m");

            // sigmar2_xabab_Lvvoo += -1.000000 f_aa(i,a) r1_bb(f,n) t2_aaaa(a,e,m,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[91]("e,m") * r1_xbb_Lvo("I,f,n");

            // sigmar2_xabab_Lvvoo += 1.000000 <e,f||m,n>_abab r0
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += V_blks_["abab_vvoo"]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += -0.500000 <i,f||a,b>_abab r1_aa(e,i) t2_abab(a,b,m,n)
            // +                      -0.500000 <i,f||b,a>_abab r1_aa(e,i) t2_abab(b,a,m,n)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[22]("f,m,i,n") * r1_xaa_Lvo("I,e,i");

            // sigmar2_xabab_Lvvoo += 1.000000 f_aa(e,m) r1_bb(f,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += F_blks_["aa_vo"]("e,m") * r1_xbb_Lvo("I,f,n");

            // sigmar2_xabab_Lvvoo += 1.000000 <e,f||a,n>_abab r1_aa(a,m)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += V_blks_["abab_vvvo"]("e,f,a,n") * r1_xaa_Lvo("I,a,m");

            // sigmar2_xabab_Lvvoo += -0.500000 <i,f||a,b>_bbbb r1_aa(e,m) t2_bbbb(a,b,n,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += 0.500000 * sigmaOps[65]("f,n") * r1_xaa_Lvo("I,e,m");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||b,a>_abab r0 t2_abab(e,a,j,i) t2_abab(b,f,m,n)
            // +                      -0.500000 <i,j||b,a>_abab r0 t2_abab(e,a,i,j) t2_abab(b,f,m,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[111]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += -1.000000 <e,i||m,a>_abab r2_bbbb(a,f,n,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += V_blks_["abba_vovo"]("e,i,a,m") * r2_xbbbb_Lvvoo("I,a,f,n,i");

            // sigmar2_xabab_Lvvoo += 1.000000 f_bb(f,n) r1_aa(e,m)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += F_blks_["bb_vo"]("f,n") * r1_xaa_Lvo("I,e,m");

            // sigmar2_xabab_Lvvoo += -1.000000 <i,j||m,a>_abab r1_bb(a,j) t2_abab(e,f,i,n)
            // flops: o3v2L1: 1, o2v2L1: 1, o3v1L1: 1 | mem: o2v2L1: 2, o2v0L1: 1,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += V_blks_["abba_oovo"]("i,j,a,m") * r1_xbb_Lvo("I,a,j") * t2_abab_vvoo("e,f,i,n");

            // sigmar2_xabab_Lvvoo += -1.000000 <e,i||m,n>_abab r1_bb(f,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= V_blks_["abab_vooo"]("e,i,m,n") * r1_xbb_Lvo("I,f,i");

            // sigmar2_xabab_Lvvoo += 1.000000 <j,i||m,a>_abab r1_aa(e,j) t2_bbbb(a,f,n,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[52]("f,j,m,n") * r1_xaa_Lvo("I,e,j");

            // sigmar2_xabab_Lvvoo += 0.250000 <j,i||a,b>_aaaa r2_abab(e,f,m,n) t2_aaaa(a,b,j,i)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += 0.250000 * scalar_5 * r2_xabab_Lvvoo("I,e,f,m,n");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||a,b>_bbbb r2_abab(e,f,m,j) t2_bbbb(a,b,n,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[76]("j,n") * r2_xabab_Lvvoo("I,e,f,m,j");

            // sigmar2_xabab_Lvvoo += 0.500000 <j,i||m,n>_abab r2_abab(e,f,j,i)
            // +                      0.500000 <i,j||m,n>_abab r2_abab(e,f,i,j)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += V_blks_["abab_oooo"]("j,i,m,n") * r2_xabab_Lvvoo("I,e,f,j,i");

            // sigmar2_xabab_Lvvoo += 0.500000 <e,f||a,b>_abab r0 t2_abab(a,b,m,n)
            // +                      0.500000 <e,f||b,a>_abab r0 t2_abab(b,a,m,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[11]("e,f,m,n") * r0_L("I");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||a,b>_aaaa r2_abab(e,f,j,n) t2_aaaa(a,b,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[75]("j,m") * r2_xabab_Lvvoo("I,e,f,j,n");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||m,a>_abab r1_bb(f,n) t2_abab(e,a,j,i)
            // +                      -0.500000 <i,j||m,a>_abab r1_bb(f,n) t2_abab(e,a,i,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[79]("e,m") * r1_xbb_Lvo("I,f,n");

            // sigmar2_xabab_Lvvoo += -0.500000 <j,i||a,m>_aaaa r1_bb(f,n) t2_aaaa(a,e,j,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[77]("e,m") * r1_xbb_Lvo("I,f,n");

            // sigmar2_xabab_Lvvoo += -1.000000 f_bb(i,a) r1_bb(a,n) t2_abab(e,f,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1, o2v1L1: 1 | mem: o2v2L1: 2, o2v0L1: 1,
            sigmar2_xabab_Lvvoo("I,e,f,m,n") -= F_blks_["bb_ov"]("i,a") * r1_xbb_Lvo("I,a,n") * t2_abab_vvoo("e,f,m,i");

            // sigmar2_xbbbb_Lvvoo += -0.500000 <j,i||j,i>_aaaa r2_bbbb(e,f,m,n)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") -= 0.500000 * scalar_2 * r2_xbbbb_Lvvoo("I,e,f,m,n");

            // sigmar2_xbbbb_Lvvoo += 0.250000 <j,i||a,b>_bbbb r0 t2_bbbb(a,b,m,n) t2_bbbb(e,f,j,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += 0.250000 * sigmaOps[106]("e,f,m,n") * r0_L("I");

            // sigmar2_xbbbb_Lvvoo += -0.500000 <j,i||j,i>_abab r2_bbbb(e,f,m,n)
            // +                      -0.500000 <i,j||i,j>_abab r2_bbbb(e,f,m,n)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") -= scalar_3 * r2_xbbbb_Lvvoo("I,e,f,m,n");

            // sigmar2_xbbbb_Lvvoo += 0.250000 <j,i||a,b>_aaaa r2_bbbb(e,f,m,n) t2_aaaa(a,b,j,i)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += 0.250000 * scalar_5 * r2_xbbbb_Lvvoo("I,e,f,m,n");

            // sigmar2_xbbbb_Lvvoo += 0.250000 <j,i||a,b>_bbbb r2_bbbb(e,f,m,n) t2_bbbb(a,b,j,i)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += 0.250000 * scalar_7 * r2_xbbbb_Lvvoo("I,e,f,m,n");

            // sigmar2_xbbbb_Lvvoo += 0.500000 <e,f||a,b>_bbbb r0 t2_bbbb(a,b,m,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += 0.500000 * sigmaOps[3]("e,f,m,n") * r0_L("I");

            // sigmar2_xbbbb_Lvvoo += 0.500000 <j,i||m,n>_bbbb r0 t2_bbbb(e,f,j,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += 0.500000 * sigmaOps[48]("e,f,m,n") * r0_L("I");

            // sigmar2_xbbbb_Lvvoo += -0.500000 <j,i||j,i>_bbbb r2_bbbb(e,f,m,n)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") -= 0.500000 * scalar_4 * r2_xbbbb_Lvvoo("I,e,f,m,n");

            // sigmar2_xbbbb_Lvvoo += 0.500000 <j,i||m,n>_bbbb r2_bbbb(e,f,j,i)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += 0.500000 * V_blks_["bbbb_oooo"]("j,i,m,n") * r2_xbbbb_Lvvoo("I,e,f,j,i");

            // sigmar2_xbbbb_Lvvoo += 0.250000 <j,i||a,b>_bbbb r2_bbbb(a,b,m,n) t2_bbbb(e,f,j,i)
            // flops: o4v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o4v0L1: 1,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += 0.250000 * V_blks_["bbbb_oovv"]("j,i,a,b") * r2_xbbbb_Lvvoo("I,a,b,m,n") * t2_bbbb_vvoo("e,f,j,i");

            // sigmar2_xbbbb_Lvvoo += 1.000000 f_bb(i,i) r2_bbbb(e,f,m,n)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += scalar_1 * r2_xbbbb_Lvvoo("I,e,f,m,n");

            // sigmar2_xbbbb_Lvvoo += 1.000000 <e,f||m,n>_bbbb r0
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += V_blks_["bbbb_vvoo"]("e,f,m,n") * r0_L("I");

            // sigmar2_xbbbb_Lvvoo += 0.250000 <j,i||a,b>_abab r2_bbbb(e,f,m,n) t2_abab(a,b,j,i)
            // +                      0.250000 <i,j||a,b>_abab r2_bbbb(e,f,m,n) t2_abab(a,b,i,j)
            // +                      0.250000 <j,i||b,a>_abab r2_bbbb(e,f,m,n) t2_abab(b,a,j,i)
            // +                      0.250000 <i,j||b,a>_abab r2_bbbb(e,f,m,n) t2_abab(b,a,i,j)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += scalar_6 * r2_xbbbb_Lvvoo("I,e,f,m,n");

            // sigmar2_xbbbb_Lvvoo += 0.500000 <j,i||b,a>_abab r0 t2_bbbb(a,e,m,n) t2_abab(b,f,j,i)
            // +                      0.500000 <i,j||b,a>_abab r0 t2_bbbb(a,e,m,n) t2_abab(b,f,i,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += sigmaOps[108]("e,f,m,n") * r0_L("I");

            // sigmar2_xbbbb_Lvvoo += -0.500000 <j,i||a,b>_abab r0 t2_abab(a,e,j,i) t2_bbbb(b,f,m,n)
            // +                      -0.500000 <i,j||a,b>_abab r0 t2_abab(a,e,i,j) t2_bbbb(b,f,m,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") -= sigmaOps[108]("f,e,m,n") * r0_L("I");

            // sigmar2_xbbbb_Lvvoo += 1.000000 f_aa(i,i) r2_bbbb(e,f,m,n)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += scalar_0 * r2_xbbbb_Lvvoo("I,e,f,m,n");

            // sigmar2_xbbbb_Lvvoo += -0.500000 <j,i||a,b>_bbbb r0 t2_bbbb(a,e,m,n) t2_bbbb(b,f,j,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[115]("e,f,m,n") * r0_L("I");

            // sigmar2_xbbbb_Lvvoo += 0.250000 <j,i||a,b>_bbbb r2_bbbb(e,f,j,i) t2_bbbb(a,b,m,n)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += 0.250000 * sigmaOps[43]("m,n,j,i") * r2_xbbbb_Lvvoo("I,e,f,j,i");

            // sigmar2_xbbbb_Lvvoo += -0.500000 <j,i||a,b>_bbbb r0 t2_bbbb(a,e,j,i) t2_bbbb(b,f,m,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[113]("e,f,m,n") * r0_L("I");

            // sigmar2_xbbbb_Lvvoo += 0.500000 <e,f||a,b>_bbbb r2_bbbb(a,b,m,n)
            // flops: o2v4L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += 0.500000 * V_blks_["bbbb_vvvv"]("e,f,a,b") * r2_xbbbb_Lvvoo("I,a,b,m,n");

            // sigmar2_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += tempPerm_xbbbb_Lvvoo("I,e,f,m,n");
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") -= tempPerm_xbbbb_Lvvoo("I,e,f,n,m");

            // tempPerm_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xbbbb_Lvvoo = TArrayD(world_, tempPerm_xbbbb_Lvvoo.trange());
            tempPerm_xbbbb_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xbbbb_Lvvoo += 1.000000 P(e,f) f_bb(e,a) r2_bbbb(a,f,m,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") = F_blks_["bb_vv"]("e,a") * r2_xbbbb_Lvvoo("I,a,f,m,n");

            // tempPerm_xbbbb_Lvvoo += -0.500000 P(e,f) <j,i||a,b>_abab r2_bbbb(b,f,m,n) t2_abab(a,e,j,i)
            // +                      -0.500000 P(e,f) <i,j||a,b>_abab r2_bbbb(b,f,m,n) t2_abab(a,e,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= sigmaOps[59]("e,b") * r2_xbbbb_Lvvoo("I,b,f,m,n");

            // tempPerm_xbbbb_Lvvoo += -0.500000 P(e,f) <j,i||a,b>_bbbb r2_bbbb(b,f,m,n) t2_bbbb(a,e,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[62]("e,b") * r2_xbbbb_Lvvoo("I,b,f,m,n");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(e,f) f_bb(i,a) r1_bb(f,i) t2_bbbb(a,e,m,n)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += sigmaOps[84]("e,m,n,i") * r1_xbb_Lvo("I,f,i");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(e,f) f_bb(e,a) r0 t2_bbbb(a,f,m,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += sigmaOps[69]("f,e,m,n") * r0_L("I");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(e,f) <i,e||m,n>_bbbb r1_bb(f,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= V_blks_["bbbb_vooo"]("e,i,m,n") * r1_xbb_Lvo("I,f,i");

            // tempPerm_xbbbb_Lvvoo += 0.500000 P(e,f) <i,e||a,b>_bbbb r1_bb(f,i) t2_bbbb(a,b,m,n)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[27]("e,i,m,n") * r1_xbb_Lvo("I,f,i");

            // sigmar2_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += tempPerm_xbbbb_Lvvoo("I,e,f,m,n");
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") -= tempPerm_xbbbb_Lvvoo("I,f,e,m,n");

            // tempPerm_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xbbbb_Lvvoo = TArrayD(world_, tempPerm_xbbbb_Lvvoo.trange());
            tempPerm_xbbbb_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xbbbb_Lvvoo += 0.500000 P(m,n) <j,i||a,n>_bbbb r1_bb(a,m) t2_bbbb(e,f,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") = 0.500000 * sigmaOps[28]("e,f,a,n") * r1_xbb_Lvo("I,a,m");

            // tempPerm_xbbbb_Lvvoo += -1.000000 P(m,n) f_bb(i,n) r2_bbbb(e,f,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= F_blks_["bb_oo"]("i,n") * r2_xbbbb_Lvvoo("I,e,f,m,i");

            // tempPerm_xbbbb_Lvvoo += -1.000000 P(m,n) <j,i||a,n>_bbbb r1_bb(a,j) t2_bbbb(e,f,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1, o3v1L1: 1 | mem: o2v2L1: 2, o2v0L1: 1,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= r1_xbb_Lvo("I,a,j") * V_blks_["bbbb_oovo"]("j,i,a,n") * t2_bbbb_vvoo("e,f,m,i");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) <j,i||a,b>_aaaa r0 t2_abab(a,e,i,n) t2_abab(b,f,j,m)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += sigmaOps[102]("f,e,m,n") * r0_L("I");

            // tempPerm_xbbbb_Lvvoo += -0.500000 P(m,n) <j,i||a,b>_bbbb r0 t2_bbbb(a,b,n,i) t2_bbbb(e,f,m,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[123]("e,f,n,m") * r0_L("I");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) <j,i||a,b>_bbbb r0 t2_bbbb(a,e,n,i) t2_bbbb(b,f,m,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += sigmaOps[101]("f,e,m,n") * r0_L("I");

            // tempPerm_xbbbb_Lvvoo += -0.500000 P(m,n) <j,i||a,b>_bbbb r2_bbbb(e,f,m,j) t2_bbbb(a,b,n,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[76]("j,n") * r2_xbbbb_Lvvoo("I,e,f,m,j");

            // tempPerm_xbbbb_Lvvoo += -0.500000 P(m,n) <i,j||a,b>_abab r2_bbbb(e,f,m,j) t2_abab(a,b,i,n)
            // +                      -0.500000 P(m,n) <i,j||b,a>_abab r2_bbbb(e,f,m,j) t2_abab(b,a,i,n)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= sigmaOps[74]("j,n") * r2_xbbbb_Lvvoo("I,e,f,m,j");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) <e,f||a,n>_bbbb r1_bb(a,m)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += V_blks_["bbbb_vvvo"]("e,f,a,n") * r1_xbb_Lvo("I,a,m");

            // tempPerm_xbbbb_Lvvoo += -1.000000 P(m,n) f_bb(i,n) r0 t2_bbbb(e,f,m,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= sigmaOps[83]("e,f,n,m") * r0_L("I");

            // tempPerm_xbbbb_Lvvoo += -1.000000 P(m,n) <j,i||a,n>_abab r1_aa(a,j) t2_bbbb(e,f,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1, o3v1L1: 1 | mem: o2v2L1: 2, o2v0L1: 1,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= r1_xaa_Lvo("I,a,j") * V_blks_["abab_oovo"]("j,i,a,n") * t2_bbbb_vvoo("e,f,m,i");

            // tempPerm_xbbbb_Lvvoo += -0.500000 P(m,n) <i,j||a,b>_abab r0 t2_abab(a,b,i,n) t2_bbbb(e,f,m,j)
            // +                      -0.500000 P(m,n) <i,j||b,a>_abab r0 t2_abab(b,a,i,n) t2_bbbb(e,f,m,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= sigmaOps[124]("e,f,m,n") * r0_L("I");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) f_bb(i,a) r1_bb(a,m) t2_bbbb(e,f,n,i)
            // flops: o3v2L1: 1, o2v2L1: 1, o2v1L1: 1 | mem: o2v2L1: 2, o2v0L1: 1,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += r1_xbb_Lvo("I,a,m") * F_blks_["bb_ov"]("i,a") * t2_bbbb_vvoo("e,f,n,i");

            // sigmar2_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += tempPerm_xbbbb_Lvvoo("I,e,f,m,n");
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") -= tempPerm_xbbbb_Lvvoo("I,e,f,n,m");

            // tempPerm_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xbbbb_Lvvoo = TArrayD(world_, tempPerm_xbbbb_Lvvoo.trange());
            tempPerm_xbbbb_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) P(e,f) <j,i||b,a>_abab r2_abab(b,f,j,m) t2_bbbb(a,e,n,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") = sigmaOps[13]("b,e,j,n") * r2_xabab_Lvvoo("I,b,f,j,m");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) P(e,f) <i,j||a,b>_abab r2_bbbb(b,f,m,j) t2_abab(a,e,i,n)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += sigmaOps[21]("b,e,j,n") * r2_xbbbb_Lvvoo("I,b,f,m,j");

            // tempPerm_xbbbb_Lvvoo += 0.500000 P(m,n) P(e,f) <j,i||a,n>_bbbb r1_bb(f,m) t2_bbbb(a,e,j,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += 0.500000 * sigmaOps[80]("e,n") * r1_xbb_Lvo("I,f,m");

            // tempPerm_xbbbb_Lvvoo += -0.500000 P(m,n) P(e,f) <i,e||a,b>_abab r1_bb(f,m) t2_abab(a,b,i,n)
            // +                      -0.500000 P(m,n) P(e,f) <i,e||b,a>_abab r1_bb(f,m) t2_abab(b,a,i,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += sigmaOps[63]("e,n") * r1_xbb_Lvo("I,f,m");

            // tempPerm_xbbbb_Lvvoo += -1.000000 P(m,n) P(e,f) f_aa(i,a) r1_bb(f,m) t2_abab(a,e,i,n)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= sigmaOps[90]("e,n") * r1_xbb_Lvo("I,f,m");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) P(e,f) <j,i||a,b>_bbbb r2_bbbb(b,f,m,j) t2_bbbb(a,e,n,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += sigmaOps[17]("e,b,n,j") * r2_xbbbb_Lvvoo("I,b,f,m,j");

            // tempPerm_xbbbb_Lvvoo += -1.000000 P(m,n) P(e,f) <i,e||a,b>_bbbb r1_bb(b,m) t2_bbbb(a,f,n,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += sigmaOps[4]("e,b,f,n") * r1_xbb_Lvo("I,b,m");

            // tempPerm_xbbbb_Lvvoo += -1.000000 P(m,n) P(e,f) f_bb(e,n) r1_bb(f,m)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= F_blks_["bb_vo"]("e,n") * r1_xbb_Lvo("I,f,m");

            // tempPerm_xbbbb_Lvvoo += 0.500000 P(m,n) P(e,f) <i,e||a,b>_bbbb r1_bb(f,m) t2_bbbb(a,b,n,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[65]("e,n") * r1_xbb_Lvo("I,f,m");

            // tempPerm_xbbbb_Lvvoo += 0.500000 P(m,n) P(e,f) <j,i||a,n>_abab r1_bb(f,m) t2_abab(a,e,j,i)
            // +                      0.500000 P(m,n) P(e,f) <i,j||a,n>_abab r1_bb(f,m) t2_abab(a,e,i,j)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += sigmaOps[78]("e,n") * r1_xbb_Lvo("I,f,m");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) P(e,f) <i,e||a,b>_abab r1_bb(b,m) t2_abab(a,f,i,n)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= sigmaOps[8]("f,e,b,n") * r1_xbb_Lvo("I,b,m");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) P(e,f) <j,i||a,b>_aaaa r2_abab(b,f,j,m) t2_abab(a,e,i,n)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += sigmaOps[18]("b,e,j,n") * r2_xabab_Lvvoo("I,b,f,j,m");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) P(e,f) <i,e||a,n>_bbbb r0 t2_bbbb(a,f,m,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= sigmaOps[29]("f,e,m,n") * r0_L("I");

            // tempPerm_xbbbb_Lvvoo += -1.000000 P(m,n) P(e,f) <i,j||a,n>_abab r1_bb(f,j) t2_abab(a,e,i,m)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= sigmaOps[45]("e,m,j,n") * r1_xbb_Lvo("I,f,j");

            // tempPerm_xbbbb_Lvvoo += -1.000000 P(m,n) P(e,f) <i,e||a,n>_abab r0 t2_abab(a,f,i,m)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += sigmaOps[30]("f,e,m,n") * r0_L("I");

            // tempPerm_xbbbb_Lvvoo += -1.000000 P(m,n) P(e,f) <j,i||a,n>_bbbb r1_bb(f,j) t2_bbbb(a,e,m,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= sigmaOps[49]("e,m,j,n") * r1_xbb_Lvo("I,f,j");

            // tempPerm_xbbbb_Lvvoo += -1.000000 P(m,n) P(e,f) <i,e||a,n>_abab r2_abab(a,f,i,m)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += V_blks_["baab_vovo"]("e,i,a,n") * r2_xabab_Lvvoo("I,a,f,i,m");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) P(e,f) f_bb(i,a) r1_bb(f,m) t2_bbbb(a,e,n,i)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += sigmaOps[89]("e,n") * r1_xbb_Lvo("I,f,m");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) P(e,f) <i,e||a,n>_bbbb r2_bbbb(a,f,m,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= V_blks_["bbbb_vovo"]("e,i,a,n") * r2_xbbbb_Lvvoo("I,a,f,m,i");

            // sigmal0_L += 0.125000 <b,a||c,d>_abab l2_abab(i,j,b,a) t2_abab(c,d,i,j)
            // +            0.125000 <b,a||d,c>_abab l2_abab(i,j,b,a) t2_abab(d,c,i,j)
            // +            0.125000 <a,b||c,d>_abab l2_abab(i,j,a,b) t2_abab(c,d,i,j)
            // +            0.125000 <a,b||d,c>_abab l2_abab(i,j,a,b) t2_abab(d,c,i,j)
            // +            0.125000 <b,a||c,d>_abab l2_abab(j,i,b,a) t2_abab(c,d,j,i)
            // +            0.125000 <b,a||d,c>_abab l2_abab(j,i,b,a) t2_abab(d,c,j,i)
            // +            0.125000 <a,b||c,d>_abab l2_abab(j,i,a,b) t2_abab(c,d,j,i)
            // +            0.125000 <a,b||d,c>_abab l2_abab(j,i,a,b) t2_abab(d,c,j,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += sigmaOps[11]("b,a,i,j") * l2_xabab_Loovv("I,i,j,b,a");

            // sigmal0_L += 0.500000 <l,k||c,d>_aaaa l2_abab(j,i,b,a) t2_aaaa(c,b,j,k) t2_abab(d,a,l,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.500000 * sigmaOps[100]("b,a,j,i") * l2_xabab_Loovv("I,j,i,b,a");

            // sigmal0_L += -0.500000 <j,a||b,c>_bbbb l1_bb(i,a) t2_bbbb(b,c,i,j)
            // flops: o1v1L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.500000 * sigmaOps[65]("a,i") * l1_xbb_Lov("I,i,a");

            // sigmal0_L += -0.250000 <l,k||c,d>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(c,b,l,k) t2_bbbb(d,a,i,j)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= 0.250000 * sigmaOps[113]("b,a,i,j") * l2_xbbbb_Loovv("I,i,j,b,a");

            // sigmal0_L += 0.500000 <l,k||c,d>_bbbb l2_aaaa(i,j,b,a) t2_abab(b,c,j,k) t2_abab(a,d,i,l)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.500000 * sigmaOps[97]("a,b,i,j") * l2_xaaaa_Loovv("I,i,j,b,a");

            // sigmal0_L += -0.500000 <j,i||j,i>_bbbb l0
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            sigmal0_L("I") -= 0.500000 * scalar_4 * l0_L("I");

            // sigmal0_L += -0.250000 <l,k||c,d>_aaaa l2_abab(j,i,b,a) t2_aaaa(c,d,j,k) t2_abab(b,a,l,i)
            // +            -0.250000 <l,k||c,d>_aaaa l2_abab(j,i,a,b) t2_aaaa(c,d,j,k) t2_abab(a,b,l,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= 0.500000 * sigmaOps[120]("b,a,j,i") * l2_xabab_Loovv("I,j,i,b,a");

            // sigmal0_L += 0.062500 <l,k||c,d>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(c,d,i,j) t2_bbbb(b,a,l,k)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.062500 * sigmaOps[106]("b,a,i,j") * l2_xbbbb_Loovv("I,i,j,b,a");

            // sigmal0_L += -1.000000 <k,b||c,j>_abab l2_abab(i,j,a,b) t2_aaaa(c,a,i,k)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += sigmaOps[32]("a,b,i,j") * l2_xabab_Loovv("I,i,j,a,b");

            // sigmal0_L += 0.125000 <b,a||c,d>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(c,d,i,j)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.125000 * sigmaOps[3]("b,a,i,j") * l2_xbbbb_Loovv("I,i,j,b,a");

            // sigmal0_L += -0.250000 <l,k||c,d>_bbbb l2_abab(i,j,a,b) t2_bbbb(c,b,l,k) t2_abab(a,d,i,j)
            // +            -0.250000 <l,k||c,d>_bbbb l2_abab(j,i,a,b) t2_bbbb(c,b,l,k) t2_abab(a,d,j,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= 0.500000 * sigmaOps[114]("a,b,i,j") * l2_xabab_Loovv("I,i,j,a,b");

            // sigmal0_L += -0.500000 f_bb(k,j) l2_abab(i,j,b,a) t2_abab(b,a,i,k)
            // +            -0.500000 f_bb(k,j) l2_abab(i,j,a,b) t2_abab(a,b,i,k)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= sigmaOps[87]("b,a,i,j") * l2_xabab_Loovv("I,i,j,b,a");

            // sigmal0_L += 1.000000 f_bb(j,b) l1_aa(i,a) t2_abab(a,b,i,j)
            // flops: o1v1L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += sigmaOps[92]("a,i") * l1_xaa_Lov("I,i,a");

            // sigmal0_L += -0.500000 <k,j||b,i>_bbbb l1_bb(i,a) t2_bbbb(b,a,k,j)
            // flops: o1v1L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= 0.500000 * sigmaOps[80]("a,i") * l1_xbb_Lov("I,i,a");

            // sigmal0_L += 0.500000 <l,k||c,d>_bbbb l2_abab(j,i,b,a) t2_abab(b,c,j,k) t2_bbbb(d,a,i,l)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.500000 * sigmaOps[99]("b,a,j,i") * l2_xabab_Loovv("I,j,i,b,a");

            // sigmal0_L += 0.500000 <l,k||c,d>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(c,b,j,k) t2_bbbb(d,a,i,l)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.500000 * sigmaOps[101]("a,b,i,j") * l2_xbbbb_Loovv("I,i,j,b,a");

            // sigmal0_L += -0.250000 <k,l||c,d>_abab l2_abab(i,j,b,a) t2_abab(c,d,k,j) t2_abab(b,a,i,l)
            // +            -0.250000 <k,l||d,c>_abab l2_abab(i,j,b,a) t2_abab(d,c,k,j) t2_abab(b,a,i,l)
            // +            -0.250000 <k,l||c,d>_abab l2_abab(i,j,a,b) t2_abab(c,d,k,j) t2_abab(a,b,i,l)
            // +            -0.250000 <k,l||d,c>_abab l2_abab(i,j,a,b) t2_abab(d,c,k,j) t2_abab(a,b,i,l)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= sigmaOps[119]("b,a,i,j") * l2_xabab_Loovv("I,i,j,b,a");

            // sigmal0_L += -0.250000 <l,k||d,c>_abab l2_aaaa(i,j,b,a) t2_abab(b,c,l,k) t2_aaaa(d,a,i,j)
            // +            -0.250000 <k,l||d,c>_abab l2_aaaa(i,j,b,a) t2_abab(b,c,k,l) t2_aaaa(d,a,i,j)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= 0.500000 * sigmaOps[107]("a,b,i,j") * l2_xaaaa_Loovv("I,i,j,b,a");

            // sigmal0_L += 0.062500 <l,k||c,d>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(c,d,i,j) t2_aaaa(b,a,l,k)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.062500 * sigmaOps[105]("b,a,i,j") * l2_xaaaa_Loovv("I,i,j,b,a");

            // sigmal0_L += -0.250000 <l,k||c,d>_abab l2_abab(j,i,b,a) t2_abab(c,d,j,k) t2_abab(b,a,l,i)
            // +            -0.250000 <l,k||d,c>_abab l2_abab(j,i,b,a) t2_abab(d,c,j,k) t2_abab(b,a,l,i)
            // +            -0.250000 <l,k||c,d>_abab l2_abab(j,i,a,b) t2_abab(c,d,j,k) t2_abab(a,b,l,i)
            // +            -0.250000 <l,k||d,c>_abab l2_abab(j,i,a,b) t2_abab(d,c,j,k) t2_abab(a,b,l,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= sigmaOps[121]("b,a,j,i") * l2_xabab_Loovv("I,j,i,b,a");

            // sigmal0_L += 0.250000 <j,i||a,b>_abab l0 t2_abab(a,b,j,i)
            // +            0.250000 <i,j||a,b>_abab l0 t2_abab(a,b,i,j)
            // +            0.250000 <j,i||b,a>_abab l0 t2_abab(b,a,j,i)
            // +            0.250000 <i,j||b,a>_abab l0 t2_abab(b,a,i,j)
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            sigmal0_L("I") += scalar_6 * l0_L("I");

            // sigmal0_L += 0.250000 <b,a||i,j>_aaaa l2_aaaa(i,j,b,a)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.250000 * V_blks_["aaaa_vvoo"]("b,a,i,j") * l2_xaaaa_Loovv("I,i,j,b,a");

            // sigmal0_L += 0.500000 <k,l||d,c>_abab l2_abab(i,j,b,a) t2_abab(b,c,k,j) t2_abab(d,a,i,l)
            // +            0.500000 <l,k||c,d>_abab l2_abab(j,i,a,b) t2_abab(c,b,j,k) t2_abab(a,d,l,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += sigmaOps[19]("b,a,i,j") * l2_xabab_Loovv("I,i,j,b,a");

            // sigmal0_L += 1.000000 f_bb(a,i) l1_bb(i,a)
            // flops: o1v1L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += F_blks_["bb_vo"]("a,i") * l1_xbb_Lov("I,i,a");

            // sigmal0_L += 0.500000 f_aa(b,c) l2_abab(i,j,b,a) t2_abab(c,a,i,j)
            // +            0.500000 f_aa(b,c) l2_abab(j,i,b,a) t2_abab(c,a,j,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += sigmaOps[68]("b,a,i,j") * l2_xabab_Loovv("I,i,j,b,a");

            // sigmal0_L += -1.000000 <b,k||c,j>_abab l2_abab(i,j,b,a) t2_abab(c,a,i,k)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= sigmaOps[41]("b,a,i,j") * l2_xabab_Loovv("I,i,j,b,a");

            // sigmal0_L += -0.500000 <j,i||j,i>_aaaa l0
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            sigmal0_L("I") -= 0.500000 * scalar_2 * l0_L("I");

            // sigmal0_L += -1.000000 <b,k||j,c>_abab l2_aaaa(i,j,b,a) t2_abab(a,c,i,k)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += sigmaOps[31]("b,a,j,i") * l2_xaaaa_Loovv("I,i,j,b,a");

            // sigmal0_L += 0.125000 <l,k||i,j>_abab l2_abab(i,j,b,a) t2_abab(b,a,l,k)
            // +            0.125000 <k,l||i,j>_abab l2_abab(i,j,b,a) t2_abab(b,a,k,l)
            // +            0.125000 <l,k||i,j>_abab l2_abab(i,j,a,b) t2_abab(a,b,l,k)
            // +            0.125000 <k,l||i,j>_abab l2_abab(i,j,a,b) t2_abab(a,b,k,l)
            // +            0.125000 <l,k||j,i>_abab l2_abab(j,i,b,a) t2_abab(b,a,l,k)
            // +            0.125000 <k,l||j,i>_abab l2_abab(j,i,b,a) t2_abab(b,a,k,l)
            // +            0.125000 <l,k||j,i>_abab l2_abab(j,i,a,b) t2_abab(a,b,l,k)
            // +            0.125000 <k,l||j,i>_abab l2_abab(j,i,a,b) t2_abab(a,b,k,l)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += sigmaOps[46]("b,a,i,j") * l2_xabab_Loovv("I,i,j,b,a");

            // sigmal0_L += 0.125000 <b,a||c,d>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(c,d,i,j)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.125000 * sigmaOps[1]("b,a,i,j") * l2_xaaaa_Loovv("I,i,j,b,a");

            // sigmal0_L += -1.000000 <b,k||j,c>_abab l2_abab(j,i,b,a) t2_bbbb(c,a,i,k)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += sigmaOps[37]("b,a,j,i") * l2_xabab_Loovv("I,j,i,b,a");

            // sigmal0_L += -0.250000 <k,l||c,d>_abab l2_bbbb(i,j,b,a) t2_abab(c,d,k,j) t2_bbbb(b,a,i,l)
            // +            -0.250000 <k,l||d,c>_abab l2_bbbb(i,j,b,a) t2_abab(d,c,k,j) t2_bbbb(b,a,i,l)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= 0.500000 * sigmaOps[124]("b,a,i,j") * l2_xbbbb_Loovv("I,i,j,b,a");

            // sigmal0_L += -0.500000 f_aa(k,j) l2_aaaa(i,j,b,a) t2_aaaa(b,a,i,k)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= 0.500000 * sigmaOps[86]("b,a,i,j") * l2_xaaaa_Loovv("I,i,j,b,a");

            // sigmal0_L += 1.000000 f_aa(i,i) l0
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            sigmal0_L("I") += scalar_0 * l0_L("I");

            // sigmal0_L += 0.125000 <l,k||i,j>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(b,a,l,k)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.125000 * sigmaOps[51]("b,a,i,j") * l2_xaaaa_Loovv("I,i,j,b,a");

            // sigmal0_L += -0.250000 <l,k||c,d>_abab l2_abab(i,j,a,b) t2_abab(c,b,l,k) t2_abab(a,d,i,j)
            // +            -0.250000 <k,l||c,d>_abab l2_abab(i,j,a,b) t2_abab(c,b,k,l) t2_abab(a,d,i,j)
            // +            -0.250000 <l,k||c,d>_abab l2_abab(j,i,a,b) t2_abab(c,b,l,k) t2_abab(a,d,j,i)
            // +            -0.250000 <k,l||c,d>_abab l2_abab(j,i,a,b) t2_abab(c,b,k,l) t2_abab(a,d,j,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= sigmaOps[112]("a,b,i,j") * l2_xabab_Loovv("I,i,j,a,b");

            // sigmal0_L += 0.250000 <b,a||i,j>_abab l2_abab(i,j,b,a)
            // +            0.250000 <a,b||i,j>_abab l2_abab(i,j,a,b)
            // +            0.250000 <b,a||j,i>_abab l2_abab(j,i,b,a)
            // +            0.250000 <a,b||j,i>_abab l2_abab(j,i,a,b)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += V_blks_["abab_vvoo"]("b,a,i,j") * l2_xabab_Loovv("I,i,j,b,a");

            // sigmal0_L += -0.500000 <k,j||b,i>_abab l1_bb(i,a) t2_abab(b,a,k,j)
            // +            -0.500000 <j,k||b,i>_abab l1_bb(i,a) t2_abab(b,a,j,k)
            // flops: o1v1L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= sigmaOps[78]("a,i") * l1_xbb_Lov("I,i,a");

            // sigmal0_L += -0.250000 <l,k||c,d>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(c,d,j,k) t2_bbbb(b,a,i,l)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= 0.250000 * sigmaOps[123]("b,a,j,i") * l2_xbbbb_Loovv("I,i,j,b,a");

            // sigmal0_L += -1.000000 <k,b||j,c>_abab l2_abab(j,i,a,b) t2_abab(a,c,k,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= sigmaOps[39]("a,b,j,i") * l2_xabab_Loovv("I,j,i,a,b");

            // sigmal0_L += -0.500000 <j,a||b,c>_aaaa l1_aa(i,a) t2_aaaa(b,c,i,j)
            // flops: o1v1L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.500000 * sigmaOps[66]("a,i") * l1_xaa_Lov("I,i,a");

            // sigmal0_L += 0.500000 <a,j||b,c>_abab l1_aa(i,a) t2_abab(b,c,i,j)
            // +            0.500000 <a,j||c,b>_abab l1_aa(i,a) t2_abab(c,b,i,j)
            // flops: o1v1L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += sigmaOps[64]("a,i") * l1_xaa_Lov("I,i,a");

            // sigmal0_L += 1.000000 f_aa(j,b) l1_bb(i,a) t2_abab(b,a,j,i)
            // flops: o1v1L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += sigmaOps[90]("a,i") * l1_xbb_Lov("I,i,a");

            // sigmal0_L += 0.500000 <l,k||c,d>_aaaa l2_bbbb(i,j,b,a) t2_abab(c,b,k,j) t2_abab(d,a,l,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.500000 * sigmaOps[102]("a,b,i,j") * l2_xbbbb_Loovv("I,i,j,b,a");

            // sigmal0_L += 1.000000 <k,b||c,j>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(c,a,i,k)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= sigmaOps[29]("a,b,i,j") * l2_xbbbb_Loovv("I,i,j,b,a");

            // sigmal0_L += 0.062500 <l,k||c,d>_abab l2_abab(i,j,b,a) t2_abab(c,d,i,j) t2_abab(b,a,l,k)
            // +            0.062500 <k,l||c,d>_abab l2_abab(i,j,b,a) t2_abab(c,d,i,j) t2_abab(b,a,k,l)
            // +            0.062500 <l,k||d,c>_abab l2_abab(i,j,b,a) t2_abab(d,c,i,j) t2_abab(b,a,l,k)
            // +            0.062500 <k,l||d,c>_abab l2_abab(i,j,b,a) t2_abab(d,c,i,j) t2_abab(b,a,k,l)
            // +            0.062500 <l,k||c,d>_abab l2_abab(i,j,a,b) t2_abab(c,d,i,j) t2_abab(a,b,l,k)
            // +            0.062500 <k,l||c,d>_abab l2_abab(i,j,a,b) t2_abab(c,d,i,j) t2_abab(a,b,k,l)
            // +            0.062500 <l,k||d,c>_abab l2_abab(i,j,a,b) t2_abab(d,c,i,j) t2_abab(a,b,l,k)
            // +            0.062500 <k,l||d,c>_abab l2_abab(i,j,a,b) t2_abab(d,c,i,j) t2_abab(a,b,k,l)
            // +            0.062500 <l,k||c,d>_abab l2_abab(j,i,b,a) t2_abab(c,d,j,i) t2_abab(b,a,l,k)
            // +            0.062500 <k,l||c,d>_abab l2_abab(j,i,b,a) t2_abab(c,d,j,i) t2_abab(b,a,k,l)
            // +            0.062500 <l,k||d,c>_abab l2_abab(j,i,b,a) t2_abab(d,c,j,i) t2_abab(b,a,l,k)
            // +            0.062500 <k,l||d,c>_abab l2_abab(j,i,b,a) t2_abab(d,c,j,i) t2_abab(b,a,k,l)
            // +            0.062500 <l,k||c,d>_abab l2_abab(j,i,a,b) t2_abab(c,d,j,i) t2_abab(a,b,l,k)
            // +            0.062500 <k,l||c,d>_abab l2_abab(j,i,a,b) t2_abab(c,d,j,i) t2_abab(a,b,k,l)
            // +            0.062500 <l,k||d,c>_abab l2_abab(j,i,a,b) t2_abab(d,c,j,i) t2_abab(a,b,l,k)
            // +            0.062500 <k,l||d,c>_abab l2_abab(j,i,a,b) t2_abab(d,c,j,i) t2_abab(a,b,k,l)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += sigmaOps[44]("b,a,i,j") * l2_xabab_Loovv("I,i,j,b,a");

            // sigmal0_L += -0.250000 <l,k||c,d>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(c,b,l,k) t2_aaaa(d,a,i,j)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= 0.250000 * sigmaOps[109]("b,a,i,j") * l2_xaaaa_Loovv("I,i,j,b,a");

            // sigmal0_L += 0.500000 <k,l||c,d>_abab l2_abab(i,j,a,b) t2_abab(c,b,k,j) t2_abab(a,d,i,l)
            // +            0.500000 <l,k||d,c>_abab l2_abab(j,i,b,a) t2_abab(b,c,j,k) t2_abab(d,a,l,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += sigmaOps[96]("a,b,i,j") * l2_xabab_Loovv("I,i,j,a,b");

            // sigmal0_L += -0.250000 <l,k||c,d>_bbbb l2_abab(i,j,b,a) t2_bbbb(c,d,j,k) t2_abab(b,a,i,l)
            // +            -0.250000 <l,k||c,d>_bbbb l2_abab(i,j,a,b) t2_bbbb(c,d,j,k) t2_abab(a,b,i,l)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= 0.500000 * sigmaOps[122]("b,a,i,j") * l2_xabab_Loovv("I,i,j,b,a");

            // sigmal0_L += 0.500000 f_bb(b,c) l2_abab(i,j,a,b) t2_abab(a,c,i,j)
            // +            0.500000 f_bb(b,c) l2_abab(j,i,a,b) t2_abab(a,c,j,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += sigmaOps[70]("a,b,i,j") * l2_xabab_Loovv("I,i,j,a,b");

            // sigmal0_L += 0.500000 <l,k||c,d>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(c,b,j,k) t2_aaaa(d,a,i,l)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.500000 * sigmaOps[98]("a,b,i,j") * l2_xaaaa_Loovv("I,i,j,b,a");

            // sigmal0_L += -1.000000 <k,b||c,j>_abab l2_bbbb(i,j,b,a) t2_abab(c,a,k,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += sigmaOps[30]("a,b,i,j") * l2_xbbbb_Loovv("I,i,j,b,a");

            // sigmal0_L += 0.250000 <j,i||a,b>_aaaa l0 t2_aaaa(a,b,j,i)
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            sigmal0_L("I") += 0.250000 * scalar_5 * l0_L("I");

            // sigmal0_L += 1.000000 f_bb(i,i) l0
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            sigmal0_L("I") += scalar_1 * l0_L("I");

            // sigmal0_L += -1.000000 f_bb(j,b) l1_bb(i,a) t2_bbbb(b,a,i,j)
            // flops: o1v1L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= sigmaOps[89]("a,i") * l1_xbb_Lov("I,i,a");

            // sigmal0_L += -0.500000 f_aa(k,j) l2_abab(j,i,b,a) t2_abab(b,a,k,i)
            // +            -0.500000 f_aa(k,j) l2_abab(j,i,a,b) t2_abab(a,b,k,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= sigmaOps[81]("b,a,j,i") * l2_xabab_Loovv("I,j,i,b,a");

            // sigmal0_L += 0.500000 f_bb(b,c) l2_bbbb(i,j,b,a) t2_bbbb(c,a,i,j)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.500000 * sigmaOps[69]("a,b,i,j") * l2_xbbbb_Loovv("I,i,j,b,a");

            // sigmal0_L += 1.000000 f_aa(a,i) l1_aa(i,a)
            // flops: o1v1L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += F_blks_["aa_vo"]("a,i") * l1_xaa_Lov("I,i,a");

            // sigmal0_L += -1.000000 f_aa(j,b) l1_aa(i,a) t2_aaaa(b,a,i,j)
            // flops: o1v1L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= sigmaOps[91]("a,i") * l1_xaa_Lov("I,i,a");

            // sigmal0_L += -0.250000 <l,k||c,d>_abab l2_aaaa(i,j,b,a) t2_abab(c,d,j,k) t2_aaaa(b,a,i,l)
            // +            -0.250000 <l,k||d,c>_abab l2_aaaa(i,j,b,a) t2_abab(d,c,j,k) t2_aaaa(b,a,i,l)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= 0.500000 * sigmaOps[117]("b,a,i,j") * l2_xaaaa_Loovv("I,i,j,b,a");

            // sigmal0_L += 0.250000 <j,i||a,b>_bbbb l0 t2_bbbb(a,b,j,i)
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            sigmal0_L("I") += 0.250000 * scalar_7 * l0_L("I");

            // sigmal0_L += 0.250000 <b,a||i,j>_bbbb l2_bbbb(i,j,b,a)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.250000 * V_blks_["bbbb_vvoo"]("b,a,i,j") * l2_xbbbb_Loovv("I,i,j,b,a");

            // sigmal0_L += 0.500000 <k,l||c,d>_abab l2_aaaa(i,j,b,a) t2_aaaa(c,b,j,k) t2_abab(a,d,i,l)
            // +            0.500000 <l,k||d,c>_abab l2_aaaa(i,j,b,a) t2_abab(b,c,j,k) t2_aaaa(d,a,i,l)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += sigmaOps[93]("b,a,j,i") * l2_xaaaa_Loovv("I,i,j,b,a");

            // sigmal0_L += -0.250000 <l,k||d,c>_abab l2_abab(i,j,b,a) t2_abab(b,c,l,k) t2_abab(d,a,i,j)
            // +            -0.250000 <k,l||d,c>_abab l2_abab(i,j,b,a) t2_abab(b,c,k,l) t2_abab(d,a,i,j)
            // +            -0.250000 <l,k||d,c>_abab l2_abab(j,i,b,a) t2_abab(b,c,l,k) t2_abab(d,a,j,i)
            // +            -0.250000 <k,l||d,c>_abab l2_abab(j,i,b,a) t2_abab(b,c,k,l) t2_abab(d,a,j,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= sigmaOps[111]("b,a,i,j") * l2_xabab_Loovv("I,i,j,b,a");

            // sigmal0_L += 0.500000 <l,k||c,d>_bbbb l2_abab(i,j,a,b) t2_bbbb(c,b,j,k) t2_abab(a,d,i,l)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.500000 * sigmaOps[103]("a,b,i,j") * l2_xabab_Loovv("I,i,j,a,b");

            // sigmal0_L += -0.250000 <l,k||c,d>_abab l2_bbbb(i,j,b,a) t2_abab(c,b,l,k) t2_bbbb(d,a,i,j)
            // +            -0.250000 <k,l||c,d>_abab l2_bbbb(i,j,b,a) t2_abab(c,b,k,l) t2_bbbb(d,a,i,j)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= 0.500000 * sigmaOps[108]("a,b,i,j") * l2_xbbbb_Loovv("I,i,j,b,a");

            // sigmal0_L += 0.125000 <l,k||i,j>_bbbb l2_bbbb(i,j,b,a) t2_bbbb(b,a,l,k)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.125000 * sigmaOps[48]("b,a,i,j") * l2_xbbbb_Loovv("I,i,j,b,a");

            // sigmal0_L += 0.500000 <l,k||d,c>_abab l2_bbbb(i,j,b,a) t2_bbbb(c,b,j,k) t2_abab(d,a,l,i)
            // +            0.500000 <k,l||c,d>_abab l2_bbbb(i,j,b,a) t2_abab(c,b,k,j) t2_bbbb(d,a,i,l)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += sigmaOps[94]("b,a,j,i") * l2_xbbbb_Loovv("I,i,j,b,a");

            // sigmal0_L += -0.500000 <k,j||i,b>_abab l1_aa(i,a) t2_abab(a,b,k,j)
            // +            -0.500000 <j,k||i,b>_abab l1_aa(i,a) t2_abab(a,b,j,k)
            // flops: o1v1L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += sigmaOps[79]("a,i") * l1_xaa_Lov("I,i,a");

            // sigmal0_L += 1.000000 <k,b||c,j>_bbbb l2_abab(i,j,a,b) t2_abab(a,c,i,k)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= sigmaOps[38]("a,b,i,j") * l2_xabab_Loovv("I,i,j,a,b");

            // sigmal0_L += 0.500000 f_aa(b,c) l2_aaaa(i,j,b,a) t2_aaaa(c,a,i,j)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.500000 * sigmaOps[72]("b,a,i,j") * l2_xaaaa_Loovv("I,i,j,b,a");

            // sigmal0_L += -0.250000 <l,k||c,d>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(c,d,j,k) t2_aaaa(b,a,i,l)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= 0.250000 * sigmaOps[118]("b,a,j,i") * l2_xaaaa_Loovv("I,i,j,b,a");

            // sigmal0_L += 1.000000 <k,b||c,j>_aaaa l2_abab(j,i,b,a) t2_abab(c,a,k,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= sigmaOps[23]("b,a,j,i") * l2_xabab_Loovv("I,j,i,b,a");

            // sigmal0_L += 0.500000 <l,k||c,d>_aaaa l2_abab(i,j,a,b) t2_abab(c,b,k,j) t2_aaaa(d,a,i,l)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += 0.500000 * sigmaOps[104]("a,b,i,j") * l2_xabab_Loovv("I,i,j,a,b");

            // sigmal0_L += -0.250000 <l,k||c,d>_aaaa l2_abab(i,j,b,a) t2_aaaa(c,b,l,k) t2_abab(d,a,i,j)
            // +            -0.250000 <l,k||c,d>_aaaa l2_abab(j,i,b,a) t2_aaaa(c,b,l,k) t2_abab(d,a,j,i)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= 0.500000 * sigmaOps[110]("b,a,i,j") * l2_xabab_Loovv("I,i,j,b,a");

            // sigmal0_L += 1.000000 <k,b||c,j>_aaaa l2_aaaa(i,j,b,a) t2_aaaa(c,a,i,k)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= sigmaOps[33]("a,b,i,j") * l2_xaaaa_Loovv("I,i,j,b,a");

            // sigmal0_L += -0.500000 <k,j||b,i>_aaaa l1_aa(i,a) t2_aaaa(b,a,k,j)
            // flops: o1v1L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= 0.500000 * sigmaOps[77]("a,i") * l1_xaa_Lov("I,i,a");

            // sigmal0_L += 0.500000 <k,l||c,d>_abab l2_abab(j,i,b,a) t2_aaaa(c,b,j,k) t2_bbbb(d,a,i,l)
            // +            0.500000 <l,k||d,c>_abab l2_abab(i,j,a,b) t2_bbbb(c,b,j,k) t2_aaaa(d,a,i,l)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") += sigmaOps[95]("b,a,j,i") * l2_xabab_Loovv("I,j,i,b,a");

            // sigmal0_L += 0.500000 <j,a||b,c>_abab l1_bb(i,a) t2_abab(b,c,j,i)
            // +            0.500000 <j,a||c,b>_abab l1_bb(i,a) t2_abab(c,b,j,i)
            // flops: o1v1L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= sigmaOps[63]("a,i") * l1_xbb_Lov("I,i,a");

            // sigmal0_L += -0.500000 <j,i||j,i>_abab l0
            // +            -0.500000 <i,j||i,j>_abab l0
            // flops: o0v0L1: 1, o0v0: 1 | mem: o0v0L1: 1, o0v0: 1,
            sigmal0_L("I") -= scalar_3 * l0_L("I");

            // sigmal0_L += -0.500000 f_bb(k,j) l2_bbbb(i,j,b,a) t2_bbbb(b,a,i,k)
            // flops: o2v2L1: 1, o0v0L1: 1 | mem: o0v0L1: 2,
            sigmal0_L("I") -= 0.500000 * sigmaOps[83]("b,a,j,i") * l2_xbbbb_Loovv("I,i,j,b,a");

            // sigmal1_xaa_Lvo += -1.000000 <j,b||c,e>_aaaa l2_aaaa(m,i,b,a) t2_aaaa(c,a,i,j)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += sigmaOps[7]("b,e,a,i") * l2_xaaaa_Loovv("I,m,i,b,a");

            // sigmal1_xaa_Lvo += -1.000000 f_aa(a,i) l2_aaaa(m,i,a,e)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= F_blks_["aa_vo"]("a,i") * l2_xaaaa_Loovv("I,m,i,a,e");

            // sigmal1_xaa_Lvo += 1.000000 f_aa(i,i) l1_aa(m,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmal1_xaa_Lvo("I,e,m") += scalar_0 * l1_xaa_Lov("I,m,e");

            // sigmal1_xaa_Lvo += -0.500000 <m,j||a,b>_abab l1_aa(i,e) t2_abab(a,b,i,j)
            // +                  -0.500000 <m,j||b,a>_abab l1_aa(i,e) t2_abab(b,a,i,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= sigmaOps[73]("m,i") * l1_xaa_Lov("I,i,e");

            // sigmal1_xaa_Lvo += 0.500000 <j,a||b,c>_abab l2_abab(m,i,e,a) t2_abab(b,c,j,i)
            // +                  0.500000 <j,a||c,b>_abab l2_abab(m,i,e,a) t2_abab(c,b,j,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= sigmaOps[63]("a,i") * l2_xabab_Loovv("I,m,i,e,a");

            // sigmal1_xaa_Lvo += -0.500000 <j,i||j,i>_bbbb l1_aa(m,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmal1_xaa_Lvo("I,e,m") -= 0.500000 * scalar_4 * l1_xaa_Lov("I,m,e");

            // sigmal1_xaa_Lvo += 0.250000 <m,a||b,c>_aaaa l2_aaaa(i,j,a,e) t2_aaaa(b,c,i,j)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= 0.250000 * sigmaOps[26]("a,m,i,j") * l2_xaaaa_Loovv("I,i,j,a,e");

            // sigmal1_xaa_Lvo += -1.000000 <m,k||j,b>_abab l2_aaaa(i,j,a,e) t2_abab(a,b,i,k)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += sigmaOps[50]("a,m,j,i") * l2_xaaaa_Loovv("I,i,j,a,e");

            // sigmal1_xaa_Lvo += 1.000000 <m,k||b,j>_aaaa l2_abab(j,i,e,a) t2_abab(b,a,k,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += sigmaOps[57]("a,m,j,i") * l2_xabab_Loovv("I,j,i,e,a");

            // sigmal1_xaa_Lvo += 0.500000 <j,a||b,c>_aaaa l2_aaaa(m,i,a,e) t2_aaaa(b,c,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= 0.500000 * sigmaOps[66]("a,i") * l2_xaaaa_Loovv("I,m,i,a,e");

            // sigmal1_xaa_Lvo += 0.250000 <j,i||a,b>_abab l1_aa(m,e) t2_abab(a,b,j,i)
            // +                  0.250000 <i,j||a,b>_abab l1_aa(m,e) t2_abab(a,b,i,j)
            // +                  0.250000 <j,i||b,a>_abab l1_aa(m,e) t2_abab(b,a,j,i)
            // +                  0.250000 <i,j||b,a>_abab l1_aa(m,e) t2_abab(b,a,i,j)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmal1_xaa_Lvo("I,e,m") += scalar_6 * l1_xaa_Lov("I,m,e");

            // sigmal1_xaa_Lvo += -1.000000 f_bb(j,b) l2_abab(m,i,e,a) t2_bbbb(b,a,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= sigmaOps[89]("a,i") * l2_xabab_Loovv("I,m,i,e,a");

            // sigmal1_xaa_Lvo += 0.250000 <k,j||e,i>_abab l2_abab(m,i,b,a) t2_abab(b,a,k,j)
            // +                  0.250000 <j,k||e,i>_abab l2_abab(m,i,b,a) t2_abab(b,a,j,k)
            // +                  0.250000 <k,j||e,i>_abab l2_abab(m,i,a,b) t2_abab(a,b,k,j)
            // +                  0.250000 <j,k||e,i>_abab l2_abab(m,i,a,b) t2_abab(a,b,j,k)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += sigmaOps[25]("b,e,a,i") * l2_xabab_Loovv("I,m,i,b,a");

            // sigmal1_xaa_Lvo += 0.500000 <k,j||i,b>_abab l2_aaaa(m,i,a,e) t2_abab(a,b,k,j)
            // +                  0.500000 <j,k||i,b>_abab l2_aaaa(m,i,a,e) t2_abab(a,b,j,k)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= sigmaOps[79]("a,i") * l2_xaaaa_Loovv("I,m,i,a,e");

            // sigmal1_xaa_Lvo += -1.000000 f_aa(m,i) l1_aa(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= F_blks_["aa_oo"]("m,i") * l1_xaa_Lov("I,i,e");

            // sigmal1_xaa_Lvo += -0.500000 <k,j||b,i>_bbbb l2_abab(m,i,e,a) t2_bbbb(b,a,k,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= 0.500000 * sigmaOps[80]("a,i") * l2_xabab_Loovv("I,m,i,e,a");

            // sigmal1_xaa_Lvo += 1.000000 <m,k||b,j>_abab l2_abab(i,j,e,a) t2_abab(b,a,i,k)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += sigmaOps[53]("a,i,m,j") * l2_xabab_Loovv("I,i,j,e,a");

            // sigmal1_xaa_Lvo += 0.250000 <j,i||a,b>_aaaa l1_aa(m,e) t2_aaaa(a,b,j,i)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmal1_xaa_Lvo("I,e,m") += 0.250000 * scalar_5 * l1_xaa_Lov("I,m,e");

            // sigmal1_xaa_Lvo += 1.000000 <m,k||j,b>_abab l2_abab(j,i,e,a) t2_bbbb(b,a,i,k)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= sigmaOps[52]("a,m,j,i") * l2_xabab_Loovv("I,j,i,e,a");

            // sigmal1_xaa_Lvo += -1.000000 <m,k||b,j>_aaaa l2_aaaa(i,j,a,e) t2_aaaa(b,a,i,k)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= sigmaOps[47]("a,m,j,i") * l2_xaaaa_Loovv("I,i,j,a,e");

            // sigmal1_xaa_Lvo += 1.000000 <j,b||c,e>_aaaa l2_abab(m,i,b,a) t2_abab(c,a,j,i)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= sigmaOps[0]("b,e,a,i") * l2_xabab_Loovv("I,m,i,b,a");

            // sigmal1_xaa_Lvo += 0.250000 <j,i||a,b>_bbbb l1_aa(m,e) t2_bbbb(a,b,j,i)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmal1_xaa_Lvo("I,e,m") += 0.250000 * scalar_7 * l1_xaa_Lov("I,m,e");

            // sigmal1_xaa_Lvo += -0.500000 <j,i||b,e>_aaaa l1_aa(m,a) t2_aaaa(b,a,j,i)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= 0.500000 * sigmaOps[61]("a,e") * l1_xaa_Lov("I,m,a");

            // sigmal1_xaa_Lvo += -0.500000 <j,i||e,b>_abab l1_aa(m,a) t2_abab(a,b,j,i)
            // +                  -0.500000 <i,j||e,b>_abab l1_aa(m,a) t2_abab(a,b,i,j)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= sigmaOps[60]("e,a") * l1_xaa_Lov("I,m,a");

            // sigmal1_xaa_Lvo += -0.500000 <k,j||b,i>_abab l2_abab(m,i,e,a) t2_abab(b,a,k,j)
            // +                  -0.500000 <j,k||b,i>_abab l2_abab(m,i,e,a) t2_abab(b,a,j,k)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= sigmaOps[78]("a,i") * l2_xabab_Loovv("I,m,i,e,a");

            // sigmal1_xaa_Lvo += 0.500000 <k,j||b,i>_aaaa l2_aaaa(m,i,a,e) t2_aaaa(b,a,k,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += 0.500000 * sigmaOps[77]("a,i") * l2_xaaaa_Loovv("I,m,i,a,e");

            // sigmal1_xaa_Lvo += 1.000000 <m,a||e,i>_abab l1_bb(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= V_blks_["baab_vovo"]("a,m,e,i") * l1_xbb_Lov("I,i,a");

            // sigmal1_xaa_Lvo += 1.000000 f_aa(a,e) l1_aa(m,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += F_blks_["aa_vv"]("a,e") * l1_xaa_Lov("I,m,a");

            // sigmal1_xaa_Lvo += -0.500000 <a,j||b,c>_abab l2_aaaa(m,i,a,e) t2_abab(b,c,i,j)
            // +                  -0.500000 <a,j||c,b>_abab l2_aaaa(m,i,a,e) t2_abab(c,b,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= sigmaOps[64]("a,i") * l2_xaaaa_Loovv("I,m,i,a,e");

            // sigmal1_xaa_Lvo += 0.500000 <b,a||e,i>_aaaa l2_aaaa(m,i,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += 0.500000 * V_blks_["aaaa_vvvo"]("b,a,e,i") * l2_xaaaa_Loovv("I,m,i,b,a");

            // sigmal1_xaa_Lvo += -1.000000 <j,b||e,c>_abab l2_abab(m,i,a,b) t2_abab(a,c,j,i)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += sigmaOps[9]("e,a,b,i") * l2_xabab_Loovv("I,m,i,a,b");

            // sigmal1_xaa_Lvo += 1.000000 f_aa(m,e) l0
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += F_blks_["aa_ov"]("m,e") * l0_L("I");

            // sigmal1_xaa_Lvo += -0.500000 f_aa(m,b) l2_abab(i,j,e,a) t2_abab(b,a,i,j)
            // +                  -0.500000 f_aa(m,b) l2_abab(j,i,e,a) t2_abab(b,a,j,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= sigmaOps[82]("a,m,i,j") * l2_xabab_Loovv("I,i,j,e,a");

            // sigmal1_xaa_Lvo += 1.000000 f_aa(j,b) l2_abab(m,i,e,a) t2_abab(b,a,j,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += sigmaOps[90]("a,i") * l2_xabab_Loovv("I,m,i,e,a");

            // sigmal1_xaa_Lvo += -0.500000 <m,a||i,j>_abab l2_abab(i,j,e,a)
            // +                  -0.500000 <m,a||j,i>_abab l2_abab(j,i,e,a)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += V_blks_["baab_vooo"]("a,m,i,j") * l2_xabab_Loovv("I,i,j,e,a");

            // sigmal1_xaa_Lvo += 1.000000 <m,a||e,i>_aaaa l1_aa(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= V_blks_["aaaa_vovo"]("a,m,e,i") * l1_xaa_Lov("I,i,a");

            // sigmal1_xaa_Lvo += -0.500000 <j,i||j,i>_aaaa l1_aa(m,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmal1_xaa_Lvo("I,e,m") -= 0.500000 * scalar_2 * l1_xaa_Lov("I,m,e");

            // sigmal1_xaa_Lvo += 0.500000 <b,a||e,i>_abab l2_abab(m,i,b,a)
            // +                  0.500000 <a,b||e,i>_abab l2_abab(m,i,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += V_blks_["abab_vvvo"]("b,a,e,i") * l2_xabab_Loovv("I,m,i,b,a");

            // sigmal1_xaa_Lvo += -0.500000 <m,j||a,b>_aaaa l1_aa(i,e) t2_aaaa(a,b,i,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= 0.500000 * sigmaOps[75]("m,i") * l1_xaa_Lov("I,i,e");

            // sigmal1_xaa_Lvo += -1.000000 <b,j||e,c>_abab l2_abab(m,i,b,a) t2_bbbb(c,a,i,j)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= sigmaOps[10]("b,e,a,i") * l2_xabab_Loovv("I,m,i,b,a");

            // sigmal1_xaa_Lvo += 0.500000 f_aa(m,b) l2_aaaa(i,j,a,e) t2_aaaa(b,a,i,j)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += 0.500000 * sigmaOps[88]("a,m,i,j") * l2_xaaaa_Loovv("I,i,j,a,e");

            // sigmal1_xaa_Lvo += -0.500000 <j,a||b,c>_bbbb l2_abab(m,i,e,a) t2_bbbb(b,c,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += 0.500000 * sigmaOps[65]("a,i") * l2_xabab_Loovv("I,m,i,e,a");

            // sigmal1_xaa_Lvo += 0.250000 <k,j||e,i>_aaaa l2_aaaa(m,i,b,a) t2_aaaa(b,a,k,j)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += 0.250000 * sigmaOps[24]("e,b,a,i") * l2_xaaaa_Loovv("I,m,i,b,a");

            // sigmal1_xaa_Lvo += 1.000000 f_aa(j,b) l2_aaaa(m,i,a,e) t2_aaaa(b,a,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += sigmaOps[91]("a,i") * l2_xaaaa_Loovv("I,m,i,a,e");

            // sigmal1_xaa_Lvo += -0.500000 <j,i||j,i>_abab l1_aa(m,e)
            // +                  -0.500000 <i,j||i,j>_abab l1_aa(m,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmal1_xaa_Lvo("I,e,m") -= scalar_3 * l1_xaa_Lov("I,m,e");

            // sigmal1_xaa_Lvo += 1.000000 <b,j||e,c>_abab l2_aaaa(m,i,b,a) t2_abab(a,c,i,j)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += sigmaOps[5]("b,e,a,i") * l2_xaaaa_Loovv("I,m,i,b,a");

            // sigmal1_xaa_Lvo += -0.250000 <m,a||b,c>_abab l2_abab(i,j,e,a) t2_abab(b,c,i,j)
            // +                  -0.250000 <m,a||c,b>_abab l2_abab(i,j,e,a) t2_abab(c,b,i,j)
            // +                  -0.250000 <m,a||b,c>_abab l2_abab(j,i,e,a) t2_abab(b,c,j,i)
            // +                  -0.250000 <m,a||c,b>_abab l2_abab(j,i,e,a) t2_abab(c,b,j,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += sigmaOps[22]("a,i,m,j") * l2_xabab_Loovv("I,i,j,e,a");

            // sigmal1_xaa_Lvo += -1.000000 f_bb(j,b) l2_aaaa(m,i,a,e) t2_abab(a,b,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= sigmaOps[92]("a,i") * l2_xaaaa_Loovv("I,m,i,a,e");

            // sigmal1_xaa_Lvo += 0.500000 <m,a||i,j>_aaaa l2_aaaa(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") -= 0.500000 * V_blks_["aaaa_vooo"]("a,m,i,j") * l2_xaaaa_Loovv("I,i,j,a,e");

            // sigmal1_xaa_Lvo += 1.000000 f_bb(a,i) l2_abab(m,i,e,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xaa_Lvo("I,e,m") += F_blks_["bb_vo"]("a,i") * l2_xabab_Loovv("I,m,i,e,a");

            // sigmal1_xaa_Lvo += 1.000000 f_bb(i,i) l1_aa(m,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmal1_xaa_Lvo("I,e,m") += scalar_1 * l1_xaa_Lov("I,m,e");

            // sigmal1_xbb_Lvo += -1.000000 f_aa(j,b) l2_bbbb(m,i,a,e) t2_abab(b,a,j,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= sigmaOps[90]("a,i") * l2_xbbbb_Loovv("I,m,i,a,e");

            // sigmal1_xbb_Lvo += -1.000000 <m,k||b,j>_bbbb l2_bbbb(i,j,a,e) t2_bbbb(b,a,i,k)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= sigmaOps[49]("a,i,m,j") * l2_xbbbb_Loovv("I,i,j,a,e");

            // sigmal1_xbb_Lvo += 0.500000 <k,j||b,i>_abab l2_bbbb(m,i,a,e) t2_abab(b,a,k,j)
            // +                  0.500000 <j,k||b,i>_abab l2_bbbb(m,i,a,e) t2_abab(b,a,j,k)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += sigmaOps[78]("a,i") * l2_xbbbb_Loovv("I,m,i,a,e");

            // sigmal1_xbb_Lvo += 1.000000 f_bb(i,i) l1_bb(m,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmal1_xbb_Lvo("I,e,m") += scalar_1 * l1_xbb_Lov("I,m,e");

            // sigmal1_xbb_Lvo += 1.000000 <a,m||i,e>_abab l1_aa(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= V_blks_["abba_vovo"]("a,m,e,i") * l1_xaa_Lov("I,i,a");

            // sigmal1_xbb_Lvo += -0.250000 <a,m||b,c>_abab l2_abab(i,j,a,e) t2_abab(b,c,i,j)
            // +                  -0.250000 <a,m||c,b>_abab l2_abab(i,j,a,e) t2_abab(c,b,i,j)
            // +                  -0.250000 <a,m||b,c>_abab l2_abab(j,i,a,e) t2_abab(b,c,j,i)
            // +                  -0.250000 <a,m||c,b>_abab l2_abab(j,i,a,e) t2_abab(c,b,j,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= sigmaOps[35]("a,i,m,j") * l2_xabab_Loovv("I,i,j,a,e");

            // sigmal1_xbb_Lvo += 0.500000 <m,a||i,j>_bbbb l2_bbbb(i,j,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= 0.500000 * V_blks_["bbbb_vooo"]("a,m,i,j") * l2_xbbbb_Loovv("I,i,j,a,e");

            // sigmal1_xbb_Lvo += 1.000000 <m,a||e,i>_bbbb l1_bb(i,a)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= V_blks_["bbbb_vovo"]("a,m,e,i") * l1_xbb_Lov("I,i,a");

            // sigmal1_xbb_Lvo += 0.500000 <j,a||b,c>_bbbb l2_bbbb(m,i,a,e) t2_bbbb(b,c,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= 0.500000 * sigmaOps[65]("a,i") * l2_xbbbb_Loovv("I,m,i,a,e");

            // sigmal1_xbb_Lvo += 0.250000 <m,a||b,c>_bbbb l2_bbbb(i,j,a,e) t2_bbbb(b,c,i,j)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= 0.250000 * sigmaOps[27]("a,m,i,j") * l2_xbbbb_Loovv("I,i,j,a,e");

            // sigmal1_xbb_Lvo += -1.000000 <k,m||b,j>_abab l2_bbbb(i,j,a,e) t2_abab(b,a,k,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= sigmaOps[45]("a,i,m,j") * l2_xbbbb_Loovv("I,i,j,a,e");

            // sigmal1_xbb_Lvo += 1.000000 <j,b||c,e>_bbbb l2_abab(i,m,a,b) t2_abab(a,c,i,j)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= sigmaOps[12]("a,b,e,i") * l2_xabab_Loovv("I,i,m,a,b");

            // sigmal1_xbb_Lvo += -0.500000 <a,m||i,j>_abab l2_abab(i,j,a,e)
            // +                  -0.500000 <a,m||j,i>_abab l2_abab(j,i,a,e)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= V_blks_["abab_vooo"]("a,m,i,j") * l2_xabab_Loovv("I,i,j,a,e");

            // sigmal1_xbb_Lvo += 0.250000 <j,i||a,b>_aaaa l1_bb(m,e) t2_aaaa(a,b,j,i)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmal1_xbb_Lvo("I,e,m") += 0.250000 * scalar_5 * l1_xbb_Lov("I,m,e");

            // sigmal1_xbb_Lvo += 1.000000 f_bb(m,e) l0
            // flops: o1v1L1: 2 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += F_blks_["bb_ov"]("m,e") * l0_L("I");

            // sigmal1_xbb_Lvo += 0.500000 f_bb(m,b) l2_bbbb(i,j,a,e) t2_bbbb(b,a,i,j)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += 0.500000 * sigmaOps[84]("a,i,j,m") * l2_xbbbb_Loovv("I,i,j,a,e");

            // sigmal1_xbb_Lvo += -1.000000 f_aa(j,b) l2_abab(i,m,a,e) t2_aaaa(b,a,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= sigmaOps[91]("a,i") * l2_xabab_Loovv("I,i,m,a,e");

            // sigmal1_xbb_Lvo += 0.500000 <a,j||b,c>_abab l2_abab(i,m,a,e) t2_abab(b,c,i,j)
            // +                  0.500000 <a,j||c,b>_abab l2_abab(i,m,a,e) t2_abab(c,b,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += sigmaOps[64]("a,i") * l2_xabab_Loovv("I,i,m,a,e");

            // sigmal1_xbb_Lvo += 1.000000 <j,b||c,e>_abab l2_bbbb(m,i,b,a) t2_abab(c,a,j,i)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= sigmaOps[8]("a,b,e,i") * l2_xbbbb_Loovv("I,m,i,b,a");

            // sigmal1_xbb_Lvo += 1.000000 f_bb(a,e) l1_bb(m,a)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += F_blks_["bb_vv"]("a,e") * l1_xbb_Lov("I,m,a");

            // sigmal1_xbb_Lvo += 0.500000 <b,a||e,i>_bbbb l2_bbbb(m,i,b,a)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += 0.500000 * V_blks_["bbbb_vvvo"]("b,a,e,i") * l2_xbbbb_Loovv("I,m,i,b,a");

            // sigmal1_xbb_Lvo += 1.000000 <k,m||b,j>_abab l2_abab(i,j,a,e) t2_aaaa(b,a,i,k)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += sigmaOps[56]("a,i,m,j") * l2_xabab_Loovv("I,i,j,a,e");

            // sigmal1_xbb_Lvo += -0.500000 <j,a||b,c>_aaaa l2_abab(i,m,a,e) t2_aaaa(b,c,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += 0.500000 * sigmaOps[66]("a,i") * l2_xabab_Loovv("I,i,m,a,e");

            // sigmal1_xbb_Lvo += -0.500000 f_bb(m,b) l2_abab(i,j,a,e) t2_abab(a,b,i,j)
            // +                  -0.500000 f_bb(m,b) l2_abab(j,i,a,e) t2_abab(a,b,j,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= sigmaOps[85]("a,i,j,m") * l2_xabab_Loovv("I,i,j,a,e");

            // sigmal1_xbb_Lvo += 1.000000 f_bb(j,b) l2_bbbb(m,i,a,e) t2_bbbb(b,a,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += sigmaOps[89]("a,i") * l2_xbbbb_Loovv("I,m,i,a,e");

            // sigmal1_xbb_Lvo += -0.500000 <j,m||a,b>_abab l1_bb(i,e) t2_abab(a,b,j,i)
            // +                  -0.500000 <j,m||b,a>_abab l1_bb(i,e) t2_abab(b,a,j,i)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= sigmaOps[74]("m,i") * l1_xbb_Lov("I,i,e");

            // sigmal1_xbb_Lvo += 1.000000 <k,m||j,b>_abab l2_abab(j,i,a,e) t2_abab(a,b,k,i)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= sigmaOps[58]("a,j,m,i") * l2_xabab_Loovv("I,j,i,a,e");

            // sigmal1_xbb_Lvo += 0.250000 <k,j||e,i>_bbbb l2_bbbb(m,i,b,a) t2_bbbb(b,a,k,j)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += 0.250000 * sigmaOps[28]("b,a,e,i") * l2_xbbbb_Loovv("I,m,i,b,a");

            // sigmal1_xbb_Lvo += 0.250000 <j,i||a,b>_abab l1_bb(m,e) t2_abab(a,b,j,i)
            // +                  0.250000 <i,j||a,b>_abab l1_bb(m,e) t2_abab(a,b,i,j)
            // +                  0.250000 <j,i||b,a>_abab l1_bb(m,e) t2_abab(b,a,j,i)
            // +                  0.250000 <i,j||b,a>_abab l1_bb(m,e) t2_abab(b,a,i,j)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmal1_xbb_Lvo("I,e,m") += scalar_6 * l1_xbb_Lov("I,m,e");

            // sigmal1_xbb_Lvo += 1.000000 f_aa(i,i) l1_bb(m,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmal1_xbb_Lvo("I,e,m") += scalar_0 * l1_xbb_Lov("I,m,e");

            // sigmal1_xbb_Lvo += 1.000000 <m,k||b,j>_bbbb l2_abab(i,j,a,e) t2_abab(a,b,i,k)
            // flops: o3v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += sigmaOps[55]("a,i,m,j") * l2_xabab_Loovv("I,i,j,a,e");

            // sigmal1_xbb_Lvo += 1.000000 f_bb(j,b) l2_abab(i,m,a,e) t2_abab(a,b,i,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += sigmaOps[92]("a,i") * l2_xabab_Loovv("I,i,m,a,e");

            // sigmal1_xbb_Lvo += -1.000000 <j,b||c,e>_abab l2_abab(i,m,a,b) t2_aaaa(c,a,i,j)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += sigmaOps[6]("a,b,e,i") * l2_xabab_Loovv("I,i,m,a,b");

            // sigmal1_xbb_Lvo += 0.500000 <b,a||i,e>_abab l2_abab(i,m,b,a)
            // +                  0.500000 <a,b||i,e>_abab l2_abab(i,m,a,b)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= V_blks_["abba_vvvo"]("b,a,e,i") * l2_xabab_Loovv("I,i,m,b,a");

            // sigmal1_xbb_Lvo += 0.250000 <k,j||i,e>_abab l2_abab(i,m,b,a) t2_abab(b,a,k,j)
            // +                  0.250000 <j,k||i,e>_abab l2_abab(i,m,b,a) t2_abab(b,a,j,k)
            // +                  0.250000 <k,j||i,e>_abab l2_abab(i,m,a,b) t2_abab(a,b,k,j)
            // +                  0.250000 <j,k||i,e>_abab l2_abab(i,m,a,b) t2_abab(a,b,j,k)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= sigmaOps[40]("b,e,a,i") * l2_xabab_Loovv("I,i,m,b,a");

            // sigmal1_xbb_Lvo += -0.500000 <k,j||i,b>_abab l2_abab(i,m,a,e) t2_abab(a,b,k,j)
            // +                  -0.500000 <j,k||i,b>_abab l2_abab(i,m,a,e) t2_abab(a,b,j,k)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += sigmaOps[79]("a,i") * l2_xabab_Loovv("I,i,m,a,e");

            // sigmal1_xbb_Lvo += -0.500000 <j,i||j,i>_aaaa l1_bb(m,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmal1_xbb_Lvo("I,e,m") -= 0.500000 * scalar_2 * l1_xbb_Lov("I,m,e");

            // sigmal1_xbb_Lvo += 0.500000 <k,j||b,i>_bbbb l2_bbbb(m,i,a,e) t2_bbbb(b,a,k,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += 0.500000 * sigmaOps[80]("a,i") * l2_xbbbb_Loovv("I,m,i,a,e");

            // sigmal1_xbb_Lvo += 1.000000 f_aa(a,i) l2_abab(i,m,a,e)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += F_blks_["aa_vo"]("a,i") * l2_xabab_Loovv("I,i,m,a,e");

            // sigmal1_xbb_Lvo += -1.000000 f_bb(m,i) l1_bb(i,e)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= F_blks_["bb_oo"]("m,i") * l1_xbb_Lov("I,i,e");

            // sigmal1_xbb_Lvo += -0.500000 <j,i||j,i>_bbbb l1_bb(m,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmal1_xbb_Lvo("I,e,m") -= 0.500000 * scalar_4 * l1_xbb_Lov("I,m,e");

            // sigmal1_xbb_Lvo += 0.250000 <j,i||a,b>_bbbb l1_bb(m,e) t2_bbbb(a,b,j,i)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmal1_xbb_Lvo("I,e,m") += 0.250000 * scalar_7 * l1_xbb_Lov("I,m,e");

            // sigmal1_xbb_Lvo += -0.500000 <j,a||b,c>_abab l2_bbbb(m,i,a,e) t2_abab(b,c,j,i)
            // +                  -0.500000 <j,a||c,b>_abab l2_bbbb(m,i,a,e) t2_abab(c,b,j,i)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += sigmaOps[63]("a,i") * l2_xbbbb_Loovv("I,m,i,a,e");

            // sigmal1_xbb_Lvo += -0.500000 <j,i||j,i>_abab l1_bb(m,e)
            // +                  -0.500000 <i,j||i,j>_abab l1_bb(m,e)
            // flops: o1v1L1: 1, o0v0: 1 | mem: o1v1L1: 1, o0v0: 1,
            sigmal1_xbb_Lvo("I,e,m") -= scalar_3 * l1_xbb_Lov("I,m,e");

            // sigmal1_xbb_Lvo += -1.000000 <b,j||c,e>_abab l2_abab(i,m,b,a) t2_abab(c,a,i,j)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= sigmaOps[2]("b,a,e,i") * l2_xabab_Loovv("I,i,m,b,a");

            // sigmal1_xbb_Lvo += -1.000000 f_bb(a,i) l2_bbbb(m,i,a,e)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= F_blks_["bb_vo"]("a,i") * l2_xbbbb_Loovv("I,m,i,a,e");

            // sigmal1_xbb_Lvo += -0.500000 <j,i||b,e>_abab l1_bb(m,a) t2_abab(b,a,j,i)
            // +                  -0.500000 <i,j||b,e>_abab l1_bb(m,a) t2_abab(b,a,i,j)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= sigmaOps[59]("a,e") * l1_xbb_Lov("I,m,a");

            // sigmal1_xbb_Lvo += -0.500000 <j,i||b,e>_bbbb l1_bb(m,a) t2_bbbb(b,a,j,i)
            // flops: o1v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= 0.500000 * sigmaOps[62]("a,e") * l1_xbb_Lov("I,m,a");

            // sigmal1_xbb_Lvo += -0.500000 <k,j||b,i>_aaaa l2_abab(i,m,a,e) t2_aaaa(b,a,k,j)
            // flops: o2v2L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= 0.500000 * sigmaOps[77]("a,i") * l2_xabab_Loovv("I,i,m,a,e");

            // sigmal1_xbb_Lvo += -0.500000 <m,j||a,b>_bbbb l1_bb(i,e) t2_bbbb(a,b,i,j)
            // flops: o2v1L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") -= 0.500000 * sigmaOps[76]("m,i") * l1_xbb_Lov("I,i,e");

            // sigmal1_xbb_Lvo += -1.000000 <j,b||c,e>_bbbb l2_bbbb(m,i,b,a) t2_bbbb(c,a,i,j)
            // flops: o2v3L1: 1, o1v1L1: 1 | mem: o1v1L1: 2,
            sigmal1_xbb_Lvo("I,e,m") += sigmaOps[4]("b,e,a,i") * l2_xbbbb_Loovv("I,m,i,b,a");

            // sigmal2_xaaaa_Lvvoo += 1.000000 f_aa(i,i) l2_aaaa(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") += scalar_0 * l2_xaaaa_Loovv("I,m,n,e,f");

            // sigmal2_xaaaa_Lvvoo += -0.500000 <j,i||j,i>_abab l2_aaaa(m,n,e,f)
            // +                      -0.500000 <i,j||i,j>_abab l2_aaaa(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") -= scalar_3 * l2_xaaaa_Loovv("I,m,n,e,f");

            // sigmal2_xaaaa_Lvvoo += 0.250000 <j,i||a,b>_bbbb l2_aaaa(m,n,e,f) t2_bbbb(a,b,j,i)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") += 0.250000 * scalar_7 * l2_xaaaa_Loovv("I,m,n,e,f");

            // sigmal2_xaaaa_Lvvoo += 0.500000 <b,a||e,f>_aaaa l2_aaaa(m,n,b,a)
            // flops: o2v4L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") += 0.500000 * V_blks_["aaaa_vvvv"]("b,a,e,f") * l2_xaaaa_Loovv("I,m,n,b,a");

            // sigmal2_xaaaa_Lvvoo += 0.500000 <m,n||i,j>_aaaa l2_aaaa(i,j,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") += 0.500000 * V_blks_["aaaa_oooo"]("m,n,i,j") * l2_xaaaa_Loovv("I,i,j,e,f");

            // sigmal2_xaaaa_Lvvoo += -0.500000 <j,i||j,i>_aaaa l2_aaaa(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") -= 0.500000 * scalar_2 * l2_xaaaa_Loovv("I,m,n,e,f");

            // sigmal2_xaaaa_Lvvoo += 0.250000 <j,i||e,f>_aaaa l2_aaaa(m,n,b,a) t2_aaaa(b,a,j,i)
            // flops: o4v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o4v0L1: 1,
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") += 0.250000 * l2_xaaaa_Loovv("I,m,n,b,a") * t2_aaaa_vvoo("b,a,j,i") * V_blks_["aaaa_oovv"]("j,i,e,f");

            // sigmal2_xaaaa_Lvvoo += -0.500000 <j,i||j,i>_bbbb l2_aaaa(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") -= 0.500000 * scalar_4 * l2_xaaaa_Loovv("I,m,n,e,f");

            // sigmal2_xaaaa_Lvvoo += 0.250000 <m,n||a,b>_aaaa l2_aaaa(i,j,e,f) t2_aaaa(a,b,i,j)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") += 0.250000 * sigmaOps[42]("i,j,m,n") * l2_xaaaa_Loovv("I,i,j,e,f");

            // sigmal2_xaaaa_Lvvoo += 1.000000 f_bb(i,i) l2_aaaa(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") += scalar_1 * l2_xaaaa_Loovv("I,m,n,e,f");

            // sigmal2_xaaaa_Lvvoo += 0.250000 <j,i||a,b>_abab l2_aaaa(m,n,e,f) t2_abab(a,b,j,i)
            // +                      0.250000 <i,j||a,b>_abab l2_aaaa(m,n,e,f) t2_abab(a,b,i,j)
            // +                      0.250000 <j,i||b,a>_abab l2_aaaa(m,n,e,f) t2_abab(b,a,j,i)
            // +                      0.250000 <i,j||b,a>_abab l2_aaaa(m,n,e,f) t2_abab(b,a,i,j)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") += scalar_6 * l2_xaaaa_Loovv("I,m,n,e,f");

            // sigmal2_xaaaa_Lvvoo += 0.250000 <j,i||a,b>_aaaa l2_aaaa(m,n,e,f) t2_aaaa(a,b,j,i)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") += 0.250000 * scalar_5 * l2_xaaaa_Loovv("I,m,n,e,f");

            // sigmal2_xaaaa_Lvvoo += 1.000000 <m,n||e,f>_aaaa l0
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") += V_blks_["aaaa_oovv"]("m,n,e,f") * l0_L("I");

            // sigmar2_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += tempPerm_xaaaa_Lvvoo("I,e,f,m,n");
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") -= tempPerm_xaaaa_Lvvoo("I,e,f,n,m");
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") -= tempPerm_xaaaa_Lvvoo("I,f,e,m,n");
            sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += tempPerm_xaaaa_Lvvoo("I,f,e,n,m");

            // tempPerm_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xaaaa_Lvvoo = TArrayD(world_, tempPerm_xaaaa_Lvvoo.trange());
            tempPerm_xaaaa_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xaaaa_Lvvoo += 1.000000 P(e,f) f_aa(a,e) l2_aaaa(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") = F_blks_["aa_vv"]("a,e") * l2_xaaaa_Loovv("I,m,n,a,f");

            // tempPerm_xaaaa_Lvvoo += -0.500000 P(e,f) <j,i||b,e>_aaaa l2_aaaa(m,n,a,f) t2_aaaa(b,a,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[61]("a,e") * l2_xaaaa_Loovv("I,m,n,a,f");

            // tempPerm_xaaaa_Lvvoo += -1.000000 P(e,f) <m,n||e,i>_aaaa l1_aa(i,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= V_blks_["aaaa_oovo"]("m,n,e,i") * l1_xaa_Lov("I,i,f");

            // tempPerm_xaaaa_Lvvoo += -0.500000 P(e,f) <j,i||e,b>_abab l2_aaaa(m,n,a,f) t2_abab(a,b,j,i)
            // +                      -0.500000 P(e,f) <i,j||e,b>_abab l2_aaaa(m,n,a,f) t2_abab(a,b,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= sigmaOps[60]("e,a") * l2_xaaaa_Loovv("I,m,n,a,f");

            // sigmal2_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") += tempPerm_xaaaa_Lvvoo("I,e,f,m,n");
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") -= tempPerm_xaaaa_Lvvoo("I,f,e,m,n");

            // tempPerm_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xaaaa_Lvvoo = TArrayD(world_, tempPerm_xaaaa_Lvvoo.trange());
            tempPerm_xaaaa_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xaaaa_Lvvoo += -0.500000 P(m,n) <n,j||a,b>_aaaa l2_aaaa(m,i,e,f) t2_aaaa(a,b,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") = -0.500000 * sigmaOps[75]("n,i") * l2_xaaaa_Loovv("I,m,i,e,f");

            // tempPerm_xaaaa_Lvvoo += -1.000000 P(m,n) f_aa(n,i) l2_aaaa(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= F_blks_["aa_oo"]("n,i") * l2_xaaaa_Loovv("I,m,i,e,f");

            // tempPerm_xaaaa_Lvvoo += -0.500000 P(m,n) <n,j||a,b>_abab l2_aaaa(m,i,e,f) t2_abab(a,b,i,j)
            // +                      -0.500000 P(m,n) <n,j||b,a>_abab l2_aaaa(m,i,e,f) t2_abab(b,a,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= sigmaOps[73]("n,i") * l2_xaaaa_Loovv("I,m,i,e,f");

            // tempPerm_xaaaa_Lvvoo += -1.000000 P(m,n) <n,a||e,f>_aaaa l1_aa(m,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += V_blks_["aaaa_vovv"]("a,n,e,f") * l1_xaa_Lov("I,m,a");

            // sigmal2_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") += tempPerm_xaaaa_Lvvoo("I,e,f,m,n");
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") -= tempPerm_xaaaa_Lvvoo("I,e,f,n,m");

            // tempPerm_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xaaaa_Lvvoo = TArrayD(world_, tempPerm_xaaaa_Lvvoo.trange());
            tempPerm_xaaaa_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) P(e,f) <n,j||e,b>_abab l2_abab(m,i,f,a) t2_bbbb(b,a,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") = sigmaOps[13]("e,a,n,i") * l2_xabab_Loovv("I,m,i,f,a");

            // tempPerm_xaaaa_Lvvoo += -1.000000 P(m,n) P(e,f) f_aa(n,e) l1_aa(m,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= F_blks_["aa_ov"]("n,e") * l1_xaa_Lov("I,m,f");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) P(e,f) <n,j||e,b>_abab l2_aaaa(m,i,a,f) t2_abab(a,b,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += sigmaOps[14]("a,e,i,n") * l2_xaaaa_Loovv("I,m,i,a,f");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) P(e,f) <n,j||b,e>_aaaa l2_aaaa(m,i,a,f) t2_aaaa(b,a,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += sigmaOps[15]("e,a,n,i") * l2_xaaaa_Loovv("I,m,i,a,f");

            // tempPerm_xaaaa_Lvvoo += -1.000000 P(m,n) P(e,f) <n,a||e,i>_abab l2_abab(m,i,f,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += V_blks_["baab_vovo"]("a,n,e,i") * l2_xabab_Loovv("I,m,i,f,a");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) P(e,f) <n,j||b,e>_aaaa l2_abab(m,i,f,a) t2_abab(b,a,j,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") += sigmaOps[18]("e,a,n,i") * l2_xabab_Loovv("I,m,i,f,a");

            // tempPerm_xaaaa_Lvvoo += 1.000000 P(m,n) P(e,f) <n,a||e,i>_aaaa l2_aaaa(m,i,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xaaaa_Lvvoo("I,e,f,m,n") -= V_blks_["aaaa_vovo"]("a,n,e,i") * l2_xaaaa_Loovv("I,m,i,a,f");

            // sigmal2_xabab_Lvvoo += -0.500000 <j,i||b,f>_abab l2_abab(m,n,e,a) t2_abab(b,a,j,i)
            // +                      -0.500000 <i,j||b,f>_abab l2_abab(m,n,e,a) t2_abab(b,a,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[59]("a,f") * l2_xabab_Loovv("I,m,n,e,a");

            // sigmal2_xabab_Lvvoo += -0.500000 <j,i||j,i>_bbbb l2_abab(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= 0.500000 * scalar_4 * l2_xabab_Loovv("I,m,n,e,f");

            // sigmal2_xabab_Lvvoo += 0.250000 <j,i||a,b>_bbbb l2_abab(m,n,e,f) t2_bbbb(a,b,j,i)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += 0.250000 * scalar_7 * l2_xabab_Loovv("I,m,n,e,f");

            // sigmal2_xabab_Lvvoo += 1.000000 <m,a||e,i>_aaaa l2_abab(i,n,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= V_blks_["aaaa_vovo"]("a,m,e,i") * l2_xabab_Loovv("I,i,n,a,f");

            // sigmal2_xabab_Lvvoo += 1.000000 f_bb(a,f) l2_abab(m,n,e,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += F_blks_["bb_vv"]("a,f") * l2_xabab_Loovv("I,m,n,e,a");

            // sigmal2_xabab_Lvvoo += 1.000000 <m,j||b,f>_abab l2_abab(i,n,e,a) t2_abab(b,a,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[36]("f,a,m,i") * l2_xabab_Loovv("I,i,n,e,a");

            // sigmal2_xabab_Lvvoo += 1.000000 <n,a||f,i>_bbbb l2_abab(m,i,e,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= V_blks_["bbbb_vovo"]("a,n,f,i") * l2_xabab_Loovv("I,m,i,e,a");

            // sigmal2_xabab_Lvvoo += -0.500000 <m,j||a,b>_aaaa l2_abab(i,n,e,f) t2_aaaa(a,b,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[75]("m,i") * l2_xabab_Loovv("I,i,n,e,f");

            // sigmal2_xabab_Lvvoo += 1.000000 <m,n||e,f>_abab l0
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += V_blks_["abab_oovv"]("m,n,e,f") * l0_L("I");

            // sigmal2_xabab_Lvvoo += -0.500000 <j,i||j,i>_abab l2_abab(m,n,e,f)
            // +                      -0.500000 <i,j||i,j>_abab l2_abab(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= scalar_3 * l2_xabab_Loovv("I,m,n,e,f");

            // sigmal2_xabab_Lvvoo += -1.000000 <m,n||i,f>_abab l1_aa(i,e)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += V_blks_["abba_oovo"]("m,n,f,i") * l1_xaa_Lov("I,i,e");

            // sigmal2_xabab_Lvvoo += -1.000000 f_aa(m,i) l2_abab(i,n,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= F_blks_["aa_oo"]("m,i") * l2_xabab_Loovv("I,i,n,e,f");

            // sigmal2_xabab_Lvvoo += 1.000000 <m,a||e,f>_abab l1_bb(n,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= V_blks_["baab_vovv"]("a,m,e,f") * l1_xbb_Lov("I,n,a");

            // sigmal2_xabab_Lvvoo += -1.000000 f_bb(n,i) l2_abab(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= F_blks_["bb_oo"]("n,i") * l2_xabab_Loovv("I,m,i,e,f");

            // sigmal2_xabab_Lvvoo += -1.000000 <m,a||i,f>_abab l2_abab(i,n,e,a)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= V_blks_["baba_vovo"]("a,m,f,i") * l2_xabab_Loovv("I,i,n,e,a");

            // sigmal2_xabab_Lvvoo += 1.000000 <m,j||b,e>_aaaa l2_bbbb(n,i,a,f) t2_abab(b,a,j,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[18]("e,a,m,i") * l2_xbbbb_Loovv("I,n,i,a,f");

            // sigmal2_xabab_Lvvoo += 1.000000 <j,n||b,f>_abab l2_abab(m,i,e,a) t2_abab(b,a,j,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[21]("f,a,n,i") * l2_xabab_Loovv("I,m,i,e,a");

            // sigmal2_xabab_Lvvoo += 0.250000 <j,i||e,f>_abab l2_abab(m,n,b,a) t2_abab(b,a,j,i)
            // +                      0.250000 <i,j||e,f>_abab l2_abab(m,n,b,a) t2_abab(b,a,i,j)
            // +                      0.250000 <j,i||e,f>_abab l2_abab(m,n,a,b) t2_abab(a,b,j,i)
            // +                      0.250000 <i,j||e,f>_abab l2_abab(m,n,a,b) t2_abab(a,b,i,j)
            // flops: o4v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o4v0L1: 1,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += l2_xabab_Loovv("I,m,n,b,a") * t2_abab_vvoo("b,a,j,i") * V_blks_["abab_oovv"]("j,i,e,f");

            // sigmal2_xabab_Lvvoo += 0.500000 <b,a||e,f>_abab l2_abab(m,n,b,a)
            // +                      0.500000 <a,b||e,f>_abab l2_abab(m,n,a,b)
            // flops: o2v4L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += V_blks_["abab_vvvv"]("b,a,e,f") * l2_xabab_Loovv("I,m,n,b,a");

            // sigmal2_xabab_Lvvoo += -1.000000 <a,n||i,f>_abab l2_aaaa(m,i,a,e)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += V_blks_["abba_vovo"]("a,n,f,i") * l2_xaaaa_Loovv("I,m,i,a,e");

            // sigmal2_xabab_Lvvoo += 1.000000 <j,n||e,b>_abab l2_abab(m,i,a,f) t2_abab(a,b,j,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[34]("e,a,n,i") * l2_xabab_Loovv("I,m,i,a,f");

            // sigmal2_xabab_Lvvoo += 1.000000 <n,j||b,f>_bbbb l2_aaaa(m,i,a,e) t2_abab(a,b,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[16]("a,f,i,n") * l2_xaaaa_Loovv("I,m,i,a,e");

            // sigmal2_xabab_Lvvoo += 1.000000 <m,j||e,b>_abab l2_bbbb(n,i,a,f) t2_bbbb(b,a,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[13]("e,a,m,i") * l2_xbbbb_Loovv("I,n,i,a,f");

            // sigmal2_xabab_Lvvoo += -0.500000 <n,j||a,b>_bbbb l2_abab(m,i,e,f) t2_bbbb(a,b,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[76]("n,i") * l2_xabab_Loovv("I,m,i,e,f");

            // sigmal2_xabab_Lvvoo += 1.000000 <a,n||e,f>_abab l1_aa(m,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += V_blks_["abab_vovv"]("a,n,e,f") * l1_xaa_Lov("I,m,a");

            // sigmal2_xabab_Lvvoo += 1.000000 <n,j||b,f>_bbbb l2_abab(m,i,e,a) t2_bbbb(b,a,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[17]("a,f,i,n") * l2_xabab_Loovv("I,m,i,e,a");

            // sigmal2_xabab_Lvvoo += -0.500000 <j,i||e,b>_abab l2_abab(m,n,a,f) t2_abab(a,b,j,i)
            // +                      -0.500000 <i,j||e,b>_abab l2_abab(m,n,a,f) t2_abab(a,b,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[60]("e,a") * l2_xabab_Loovv("I,m,n,a,f");

            // sigmal2_xabab_Lvvoo += -0.500000 <m,j||a,b>_abab l2_abab(i,n,e,f) t2_abab(a,b,i,j)
            // +                      -0.500000 <m,j||b,a>_abab l2_abab(i,n,e,f) t2_abab(b,a,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[73]("m,i") * l2_xabab_Loovv("I,i,n,e,f");

            // sigmal2_xabab_Lvvoo += 1.000000 <m,j||e,b>_abab l2_abab(i,n,a,f) t2_abab(a,b,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[14]("a,e,i,m") * l2_xabab_Loovv("I,i,n,a,f");

            // sigmal2_xabab_Lvvoo += 0.250000 <m,n||a,b>_abab l2_abab(i,j,e,f) t2_abab(a,b,i,j)
            // +                      0.250000 <m,n||b,a>_abab l2_abab(i,j,e,f) t2_abab(b,a,i,j)
            // +                      0.250000 <m,n||a,b>_abab l2_abab(j,i,e,f) t2_abab(a,b,j,i)
            // +                      0.250000 <m,n||b,a>_abab l2_abab(j,i,e,f) t2_abab(b,a,j,i)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[54]("i,m,j,n") * l2_xabab_Loovv("I,i,j,e,f");

            // sigmal2_xabab_Lvvoo += 1.000000 f_bb(n,f) l1_aa(m,e)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += F_blks_["bb_ov"]("n,f") * l1_xaa_Lov("I,m,e");

            // sigmal2_xabab_Lvvoo += -0.500000 <j,n||a,b>_abab l2_abab(m,i,e,f) t2_abab(a,b,j,i)
            // +                      -0.500000 <j,n||b,a>_abab l2_abab(m,i,e,f) t2_abab(b,a,j,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= sigmaOps[74]("n,i") * l2_xabab_Loovv("I,m,i,e,f");

            // sigmal2_xabab_Lvvoo += -0.500000 <j,i||j,i>_aaaa l2_abab(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= 0.500000 * scalar_2 * l2_xabab_Loovv("I,m,n,e,f");

            // sigmal2_xabab_Lvvoo += 1.000000 <j,n||b,f>_abab l2_aaaa(m,i,a,e) t2_aaaa(b,a,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[20]("a,f,i,n") * l2_xaaaa_Loovv("I,m,i,a,e");

            // sigmal2_xabab_Lvvoo += 1.000000 <m,j||b,e>_aaaa l2_abab(i,n,a,f) t2_aaaa(b,a,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += sigmaOps[15]("e,a,m,i") * l2_xabab_Loovv("I,i,n,a,f");

            // sigmal2_xabab_Lvvoo += 1.000000 f_aa(i,i) l2_abab(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += scalar_0 * l2_xabab_Loovv("I,m,n,e,f");

            // sigmal2_xabab_Lvvoo += 1.000000 f_bb(i,i) l2_abab(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += scalar_1 * l2_xabab_Loovv("I,m,n,e,f");

            // sigmal2_xabab_Lvvoo += 0.500000 <m,n||i,j>_abab l2_abab(i,j,e,f)
            // +                      0.500000 <m,n||j,i>_abab l2_abab(j,i,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += V_blks_["abab_oooo"]("m,n,i,j") * l2_xabab_Loovv("I,i,j,e,f");

            // sigmal2_xabab_Lvvoo += 0.250000 <j,i||a,b>_aaaa l2_abab(m,n,e,f) t2_aaaa(a,b,j,i)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += 0.250000 * scalar_5 * l2_xabab_Loovv("I,m,n,e,f");

            // sigmal2_xabab_Lvvoo += -0.500000 <j,i||b,e>_aaaa l2_abab(m,n,a,f) t2_aaaa(b,a,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[61]("a,e") * l2_xabab_Loovv("I,m,n,a,f");

            // sigmal2_xabab_Lvvoo += -0.500000 <j,i||b,f>_bbbb l2_abab(m,n,e,a) t2_bbbb(b,a,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[62]("a,f") * l2_xabab_Loovv("I,m,n,e,a");

            // sigmal2_xabab_Lvvoo += -1.000000 <m,a||e,i>_abab l2_bbbb(n,i,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += V_blks_["baab_vovo"]("a,m,e,i") * l2_xbbbb_Loovv("I,n,i,a,f");

            // sigmal2_xabab_Lvvoo += -1.000000 <m,n||e,i>_abab l1_bb(i,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= V_blks_["abab_oovo"]("m,n,e,i") * l1_xbb_Lov("I,i,f");

            // sigmal2_xabab_Lvvoo += -1.000000 <a,n||e,i>_abab l2_abab(m,i,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") -= V_blks_["abab_vovo"]("a,n,e,i") * l2_xabab_Loovv("I,m,i,a,f");

            // sigmal2_xabab_Lvvoo += 0.250000 <j,i||a,b>_abab l2_abab(m,n,e,f) t2_abab(a,b,j,i)
            // +                      0.250000 <i,j||a,b>_abab l2_abab(m,n,e,f) t2_abab(a,b,i,j)
            // +                      0.250000 <j,i||b,a>_abab l2_abab(m,n,e,f) t2_abab(b,a,j,i)
            // +                      0.250000 <i,j||b,a>_abab l2_abab(m,n,e,f) t2_abab(b,a,i,j)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += scalar_6 * l2_xabab_Loovv("I,m,n,e,f");

            // sigmal2_xabab_Lvvoo += 1.000000 f_aa(m,e) l1_bb(n,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += F_blks_["aa_ov"]("m,e") * l1_xbb_Lov("I,n,f");

            // sigmal2_xabab_Lvvoo += 1.000000 f_aa(a,e) l2_abab(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xabab_Lvvoo("I,e,f,m,n") += F_blks_["aa_vv"]("a,e") * l2_xabab_Loovv("I,m,n,a,f");

            // sigmal2_xbbbb_Lvvoo += 1.000000 f_aa(i,i) l2_bbbb(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") += scalar_0 * l2_xbbbb_Loovv("I,m,n,e,f");

            // sigmal2_xbbbb_Lvvoo += -0.500000 <j,i||j,i>_bbbb l2_bbbb(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") -= 0.500000 * scalar_4 * l2_xbbbb_Loovv("I,m,n,e,f");

            // sigmal2_xbbbb_Lvvoo += 0.500000 <m,n||i,j>_bbbb l2_bbbb(i,j,e,f)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") += 0.500000 * V_blks_["bbbb_oooo"]("m,n,i,j") * l2_xbbbb_Loovv("I,i,j,e,f");

            // sigmal2_xbbbb_Lvvoo += 0.250000 <m,n||a,b>_bbbb l2_bbbb(i,j,e,f) t2_bbbb(a,b,i,j)
            // flops: o4v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") += 0.250000 * sigmaOps[43]("i,j,m,n") * l2_xbbbb_Loovv("I,i,j,e,f");

            // sigmal2_xbbbb_Lvvoo += -0.500000 <j,i||j,i>_abab l2_bbbb(m,n,e,f)
            // +                      -0.500000 <i,j||i,j>_abab l2_bbbb(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") -= scalar_3 * l2_xbbbb_Loovv("I,m,n,e,f");

            // sigmal2_xbbbb_Lvvoo += 0.500000 <b,a||e,f>_bbbb l2_bbbb(m,n,b,a)
            // flops: o2v4L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") += 0.500000 * V_blks_["bbbb_vvvv"]("b,a,e,f") * l2_xbbbb_Loovv("I,m,n,b,a");

            // sigmal2_xbbbb_Lvvoo += 0.250000 <j,i||e,f>_bbbb l2_bbbb(m,n,b,a) t2_bbbb(b,a,j,i)
            // flops: o4v2L1: 2, o2v2L1: 1 | mem: o2v2L1: 2, o4v0L1: 1,
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") += 0.250000 * l2_xbbbb_Loovv("I,m,n,b,a") * t2_bbbb_vvoo("b,a,j,i") * V_blks_["bbbb_oovv"]("j,i,e,f");

            // sigmal2_xbbbb_Lvvoo += 1.000000 f_bb(i,i) l2_bbbb(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") += scalar_1 * l2_xbbbb_Loovv("I,m,n,e,f");

            // sigmal2_xbbbb_Lvvoo += 0.250000 <j,i||a,b>_bbbb l2_bbbb(m,n,e,f) t2_bbbb(a,b,j,i)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") += 0.250000 * scalar_7 * l2_xbbbb_Loovv("I,m,n,e,f");

            // sigmal2_xbbbb_Lvvoo += 1.000000 <m,n||e,f>_bbbb l0
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") += V_blks_["bbbb_oovv"]("m,n,e,f") * l0_L("I");

            // sigmal2_xbbbb_Lvvoo += -0.500000 <j,i||j,i>_aaaa l2_bbbb(m,n,e,f)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") -= 0.500000 * scalar_2 * l2_xbbbb_Loovv("I,m,n,e,f");

            // sigmal2_xbbbb_Lvvoo += 0.250000 <j,i||a,b>_aaaa l2_bbbb(m,n,e,f) t2_aaaa(a,b,j,i)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") += 0.250000 * scalar_5 * l2_xbbbb_Loovv("I,m,n,e,f");

            // sigmal2_xbbbb_Lvvoo += 0.250000 <j,i||a,b>_abab l2_bbbb(m,n,e,f) t2_abab(a,b,j,i)
            // +                      0.250000 <i,j||a,b>_abab l2_bbbb(m,n,e,f) t2_abab(a,b,i,j)
            // +                      0.250000 <j,i||b,a>_abab l2_bbbb(m,n,e,f) t2_abab(b,a,j,i)
            // +                      0.250000 <i,j||b,a>_abab l2_bbbb(m,n,e,f) t2_abab(b,a,i,j)
            // flops: o2v2L1: 1, o0v0: 1 | mem: o2v2L1: 1, o0v0: 1,
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") += scalar_6 * l2_xbbbb_Loovv("I,m,n,e,f");

            // sigmar2_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += tempPerm_xbbbb_Lvvoo("I,e,f,m,n");
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") -= tempPerm_xbbbb_Lvvoo("I,e,f,n,m");
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") -= tempPerm_xbbbb_Lvvoo("I,f,e,m,n");
            sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += tempPerm_xbbbb_Lvvoo("I,f,e,n,m");

            // tempPerm_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xbbbb_Lvvoo = TArrayD(world_, tempPerm_xbbbb_Lvvoo.trange());
            tempPerm_xbbbb_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xbbbb_Lvvoo += -0.500000 P(e,f) <j,i||b,e>_abab l2_bbbb(m,n,a,f) t2_abab(b,a,j,i)
            // +                      -0.500000 P(e,f) <i,j||b,e>_abab l2_bbbb(m,n,a,f) t2_abab(b,a,i,j)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") = -1.0 * sigmaOps[59]("a,e") * l2_xbbbb_Loovv("I,m,n,a,f");

            // tempPerm_xbbbb_Lvvoo += -0.500000 P(e,f) <j,i||b,e>_bbbb l2_bbbb(m,n,a,f) t2_bbbb(b,a,j,i)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= 0.500000 * sigmaOps[62]("a,e") * l2_xbbbb_Loovv("I,m,n,a,f");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(e,f) f_bb(a,e) l2_bbbb(m,n,a,f)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += F_blks_["bb_vv"]("a,e") * l2_xbbbb_Loovv("I,m,n,a,f");

            // tempPerm_xbbbb_Lvvoo += -1.000000 P(e,f) <m,n||e,i>_bbbb l1_bb(i,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= V_blks_["bbbb_oovo"]("m,n,e,i") * l1_xbb_Lov("I,i,f");

            // sigmal2_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") += tempPerm_xbbbb_Lvvoo("I,e,f,m,n");
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") -= tempPerm_xbbbb_Lvvoo("I,f,e,m,n");

            // tempPerm_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xbbbb_Lvvoo = TArrayD(world_, tempPerm_xbbbb_Lvvoo.trange());
            tempPerm_xbbbb_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xbbbb_Lvvoo += -0.500000 P(m,n) <n,j||a,b>_bbbb l2_bbbb(m,i,e,f) t2_bbbb(a,b,i,j)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") = -0.500000 * sigmaOps[76]("n,i") * l2_xbbbb_Loovv("I,m,i,e,f");

            // tempPerm_xbbbb_Lvvoo += -0.500000 P(m,n) <j,n||a,b>_abab l2_bbbb(m,i,e,f) t2_abab(a,b,j,i)
            // +                      -0.500000 P(m,n) <j,n||b,a>_abab l2_bbbb(m,i,e,f) t2_abab(b,a,j,i)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= sigmaOps[74]("n,i") * l2_xbbbb_Loovv("I,m,i,e,f");

            // tempPerm_xbbbb_Lvvoo += -1.000000 P(m,n) <n,a||e,f>_bbbb l1_bb(m,a)
            // flops: o2v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += V_blks_["bbbb_vovv"]("a,n,e,f") * l1_xbb_Lov("I,m,a");

            // tempPerm_xbbbb_Lvvoo += -1.000000 P(m,n) f_bb(n,i) l2_bbbb(m,i,e,f)
            // flops: o3v2L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= F_blks_["bb_oo"]("n,i") * l2_xbbbb_Loovv("I,m,i,e,f");

            // sigmal2_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") += tempPerm_xbbbb_Lvvoo("I,e,f,m,n");
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") -= tempPerm_xbbbb_Lvvoo("I,e,f,n,m");

            // tempPerm_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo --> Apply Perm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xbbbb_Lvvoo = TArrayD(world_, tempPerm_xbbbb_Lvvoo.trange());
            tempPerm_xbbbb_Lvvoo.fill(0.0); world_.gop.fence();


            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) P(e,f) <j,n||b,e>_abab l2_bbbb(m,i,a,f) t2_abab(b,a,j,i)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") = sigmaOps[21]("e,a,n,i") * l2_xbbbb_Loovv("I,m,i,a,f");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) P(e,f) <j,n||b,e>_abab l2_abab(i,m,a,f) t2_aaaa(b,a,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += sigmaOps[20]("a,e,i,n") * l2_xabab_Loovv("I,i,m,a,f");

            // tempPerm_xbbbb_Lvvoo += -1.000000 P(m,n) P(e,f) <a,n||i,e>_abab l2_abab(i,m,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += V_blks_["abba_vovo"]("a,n,e,i") * l2_xabab_Loovv("I,i,m,a,f");

            // tempPerm_xbbbb_Lvvoo += -1.000000 P(m,n) P(e,f) f_bb(n,e) l1_bb(m,f)
            // flops: o2v2L1: 2 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= F_blks_["bb_ov"]("n,e") * l1_xbb_Lov("I,m,f");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) P(e,f) <n,a||e,i>_bbbb l2_bbbb(m,i,a,f)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") -= V_blks_["bbbb_vovo"]("a,n,e,i") * l2_xbbbb_Loovv("I,m,i,a,f");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) P(e,f) <n,j||b,e>_bbbb l2_bbbb(m,i,a,f) t2_bbbb(b,a,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += sigmaOps[17]("a,e,i,n") * l2_xbbbb_Loovv("I,m,i,a,f");

            // tempPerm_xbbbb_Lvvoo += 1.000000 P(m,n) P(e,f) <n,j||b,e>_bbbb l2_abab(i,m,a,f) t2_abab(a,b,i,j)
            // flops: o3v3L1: 1, o2v2L1: 1 | mem: o2v2L1: 2,
            tempPerm_xbbbb_Lvvoo("I,e,f,m,n") += sigmaOps[16]("a,e,i,n") * l2_xabab_Loovv("I,i,m,a,f");

            // sigmal2_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo tempPerm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") += tempPerm_xaaaa_Lvvoo("I,e,f,m,n");
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") -= tempPerm_xaaaa_Lvvoo("I,e,f,n,m");
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") -= tempPerm_xaaaa_Lvvoo("I,f,e,m,n");
            sigmal2_xaaaa_Lvvoo("I,e,f,m,n") += tempPerm_xaaaa_Lvvoo("I,f,e,n,m");

            // tempPerm_xaaaa_Lvvoo += 1.000000 tempPerm_xaaaa_Lvvoo tempPerm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xaaaa_Lvvoo = TArrayD(world_, tempPerm_xaaaa_Lvvoo.trange());
            tempPerm_xaaaa_Lvvoo.fill(0.0); world_.gop.fence();


            // sigmal2_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo tempPerm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") += tempPerm_xbbbb_Lvvoo("I,e,f,m,n");
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") -= tempPerm_xbbbb_Lvvoo("I,e,f,n,m");
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") -= tempPerm_xbbbb_Lvvoo("I,f,e,m,n");
            sigmal2_xbbbb_Lvvoo("I,e,f,m,n") += tempPerm_xbbbb_Lvvoo("I,f,e,n,m");

            // tempPerm_xbbbb_Lvvoo += 1.000000 tempPerm_xbbbb_Lvvoo tempPerm
            // flops: o2v2L1: 1 | mem: o2v2L1: 1,
            tempPerm_xbbbb_Lvvoo = TArrayD(world_, tempPerm_xbbbb_Lvvoo.trange());
            tempPerm_xbbbb_Lvvoo.fill(0.0); world_.gop.fence();

        }
        world_.gop.fence(); tempArray[0].~TArrayD();
        world_.gop.fence(); tempArray[1].~TArrayD();
        world_.gop.fence(); tempArray[2].~TArrayD();
        world_.gop.fence(); tempArray[3].~TArrayD();
        world_.gop.fence(); tempArray[4].~TArrayD();
        world_.gop.fence(); tempArray[5].~TArrayD();
        world_.gop.fence(); tempArray[6].~TArrayD();
        world_.gop.fence(); tempArray[7].~TArrayD();
        world_.gop.fence(); tempArray[8].~TArrayD();
        world_.gop.fence(); tempArray[9].~TArrayD();
        world_.gop.fence(); tempArray[10].~TArrayD();
        world_.gop.fence(); tempArray[11].~TArrayD();
        world_.gop.fence(); tempArray[12].~TArrayD();
        world_.gop.fence(); tempArray[13].~TArrayD();
        world_.gop.fence(); tempArray[14].~TArrayD();
        world_.gop.fence(); tempArray[15].~TArrayD();
        world_.gop.fence(); tempArray[16].~TArrayD();
        world_.gop.fence(); tempArray[17].~TArrayD();
        world_.gop.fence(); tempArray[18].~TArrayD();
        world_.gop.fence(); tempArray[19].~TArrayD();
        world_.gop.fence(); tempArray[20].~TArrayD();
        world_.gop.fence(); tempArray[21].~TArrayD();
        world_.gop.fence(); tempArray[22].~TArrayD();
        world_.gop.fence(); tempArray[23].~TArrayD();
        world_.gop.fence(); tempArray[24].~TArrayD();
        world_.gop.fence(); tempArray[25].~TArrayD();
        world_.gop.fence(); tempArray[26].~TArrayD();
        world_.gop.fence(); tempArray[27].~TArrayD();
        world_.gop.fence(); tempArray[28].~TArrayD();
        world_.gop.fence(); tempArray[29].~TArrayD();
        world_.gop.fence(); tempArray[30].~TArrayD();
        world_.gop.fence(); tempArray[31].~TArrayD();
        world_.gop.fence(); tempArray[32].~TArrayD();
        world_.gop.fence(); tempArray[33].~TArrayD();

        /*

        Total Number of Terms: 746
        Number of Flops: (old) 2814 -> (new) 1711

        Total FLOP scaling:
        ------------------
           Scaling :      new |      old |     diff
          -------- : -------- | -------- | --------
            o2v4L1 :        6 |        8 |       -2
            o3v3L1 :       56 |       56 |        0
            o4v2L1 :       24 |       44 |      -20
            o2v3L1 :      102 |      156 |      -54
              o2v4 :       13 |       34 |      -21
            o3v2L1 :      114 |      176 |      -62
              o3v3 :       42 |      150 |     -108
              o4v2 :       20 |       94 |      -74
            o1v3L1 :       12 |       20 |       -8
            o2v2L1 :      646 |     1066 |     -420
              o2v3 :       25 |      154 |     -129
            o3v1L1 :       16 |       20 |       -4
              o3v2 :       24 |      152 |     -128
            o1v2L1 :       12 |       16 |       -4
              o1v3 :       14 |        0 |       14
            o2v1L1 :       20 |       26 |       -6
              o2v2 :       73 |       24 |       49
              o3v1 :       18 |        0 |       18
              o4v0 :        3 |        0 |        3
            o0v2L1 :       14 |        0 |       14
            o1v1L1 :      236 |      304 |      -68
            o2v0L1 :       10 |        0 |       10
              o0v2 :        5 |        0 |        5
              o1v1 :       12 |        0 |       12
              o2v0 :        4 |        0 |        4
            o0v0L1 :       94 |      170 |      -76
              o0v0 :       96 |      144 |      -48

        Total MEM scaling:
            o2v2L1 :      748 |      874 |     -126
            o4v0L1 :        6 |       12 |       -6
              o1v3 :       28 |       36 |       -8
              o2v2 :      139 |      256 |     -117
              o3v1 :       36 |       46 |      -10
              o4v0 :        7 |       36 |      -29
            o0v2L1 :       28 |       44 |      -16
            o1v1L1 :      376 |      496 |     -120
            o2v0L1 :       32 |       54 |      -22
              o0v2 :       11 |       72 |      -61
              o1v1 :       24 |       96 |      -72
              o2v0 :        8 |       66 |      -58
            o0v0L1 :      172 |      316 |     -144
              o0v0 :       96 |      144 |      -48
         */

    }

    double enuc_ = cc_wfn_->enuc_;
    sigmar0_L("I") += r0_L("I") * (enuc_);
    sigmal0_L("I") += l0_L("I") * (enuc_);

    sigmar1_xaa_Lvo("I,e,m") += r1_xaa_Lvo("I,e,m") * (enuc_);
    sigmar1_xbb_Lvo("I,e,m") += r1_xbb_Lvo("I,e,m") * (enuc_);

    sigmal1_xaa_Lvo("I,e,m") += l1_xaa_Lov("I,m,e") * (enuc_);
    sigmal1_xbb_Lvo("I,e,m") += l1_xbb_Lov("I,m,e") * (enuc_);

    sigmar2_xaaaa_Lvvoo("I,e,f,m,n") += r2_xaaaa_Lvvoo("I,e,f,m,n") * (enuc_);
    sigmar2_xabab_Lvvoo("I,e,f,m,n") += r2_xabab_Lvvoo("I,e,f,m,n") * (enuc_);
    sigmar2_xbbbb_Lvvoo("I,e,f,m,n") += r2_xbbbb_Lvvoo("I,e,f,m,n") * (enuc_);

    sigmal2_xaaaa_Lvvoo("I,e,f,m,n") += l2_xaaaa_Loovv("I,m,n,e,f") * (enuc_);
    sigmal2_xabab_Lvvoo("I,e,f,m,n") += l2_xabab_Loovv("I,m,n,e,f") * (enuc_);
    sigmal2_xbbbb_Lvvoo("I,e,f,m,n") += l2_xbbbb_Loovv("I,m,n,e,f") * (enuc_);

    world_.gop.fence();
}