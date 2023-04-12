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

#include "../../include/derived/eom_ee_rdm.h"

namespace hilbert {
    EOM_EE_RDM::EOM_EE_RDM(const shared_ptr<EOM_Driver>& eom_driver, Options & options) : EOM_RDM(eom_driver, options) {
        init_rdms();
    }

    void EOM_EE_RDM::init_rdms(){

        // initialize the dimensions
        vector<string> dims = {"oa", "va", "ob", "vb"};
        vector<size_t> dim_vec = {oa_, va_, ob_, vb_};

        // initialize the 1-RDMs
        for (int i = 0; i < dims.size(); i++) {
            for (int j = 0; j < dims.size(); j++) {
                string dim1 = dims[i], dim2 = dims[j];

                string ov = dim1.substr(0, 1) + dim2.substr(0, 1);
                string spin = dim1.substr(1, 1) + dim2.substr(1, 1);
                string rdm_str = "D1_" + spin + "_" + ov;
                RDM_blks_[rdm_str] = TArrayD(world_, HelperD::makeRange({M_, M_, dim_vec[i], dim_vec[j]}));
                RDM_blks_[rdm_str].fill(0.0);

                // initialize the 2-RDMs
                for (int k = 0; k < dims.size(); k++) {
                    for (int l = 0; l < dims.size(); l++) {
                        string dim3 = dims[k], dim4 = dims[l];

                        ov += dim3.substr(0, 1) + dim4.substr(0, 1);
                        spin += dim3.substr(1, 1) + dim4.substr(1, 1);
                        rdm_str = "D2_" + spin + "_" + ov;

                        RDM_blks_[rdm_str] = TArrayD(world_, HelperD::makeRange({M_, M_, dim_vec[i], dim_vec[j], dim_vec[k], dim_vec[l]}));
                        RDM_blks_[rdm_str].fill(0.0);
                    }
                }
            }
        }
        world_.gop.fence();
    }
    
    void EOM_EE_RDM::compute_eom_rdms(){

        Timer rdm_timer; rdm_timer.start();

        if (world_.rank() == 0) {
            outfile->Printf("\n  Building %s transition RDMs... ", eom_driver_->eom_type_.c_str());
        }
        // get cavity information
        double w0 = eom_driver_->cc_wfn_->cavity_frequency_[2];
        double coupling_factor_z = w0 * eom_driver_->cc_wfn_->cavity_coupling_strength_[2];
        
        TArrayMap & amplitudes = eom_driver_->cc_wfn_->amplitudes_;
        TArrayMap & evecs = eom_driver_->evec_blks_;
        TArrayMap & Id_blks = eom_driver_->cc_wfn_->Id_blks_;
        double u0 = eom_driver_->cc_wfn_->scalar_amps_["u0"];

        bool include_s0_ = include_u0_;
        bool include_m0_ = include_u0_;

        bool include_s1_ = include_u1_;
        bool include_m1_ = include_u1_;

        bool include_s2_ = include_u2_;
        bool include_m2_ = include_u2_;

        // compute 1RDM
        {
            // zero out 1RDM
            for (auto &rdm : RDM_blks_) {
                if (rdm.first.find("D1") == string::npos) continue;
                else HelperD::forall(rdm.second, [](auto &tile, auto &x) { tile[x] = 0.0; });
            }

            // initialize wavefunction amplitudes
            TA::TArrayD
                &t1_aa = amplitudes["t1_aa"], &t1_bb = amplitudes["t1_bb"], // t1 amplitudes
                &u1_aa = amplitudes["u1_aa"], &u1_bb = amplitudes["u1_bb"], // u1 amplitudes
                &t2_aaaa = amplitudes["t2_aaaa"], &t2_abab = amplitudes["t2_abab"], &t2_bbbb = amplitudes["t2_bbbb"],
                &u2_aaaa = amplitudes["u2_aaaa"], &u2_abab = amplitudes["u2_abab"], &u2_bbbb = amplitudes["u2_bbbb"];

            TA::TArrayD
                &r0 = evecs["r0"],
                &r1_aa = evecs["r1_aa"], &r1_bb = evecs["r1_bb"],
                &r2_aaaa = evecs["r2_aaaa"], &r2_abab = evecs["r2_abab"], &r2_bbbb = evecs["r2_bbbb"],
                &l0 = evecs["l0"],
                &l1_aa = evecs["l1_aa"], &l1_bb = evecs["l1_bb"],
                &l2_aaaa = evecs["l2_aaaa"], &l2_abab = evecs["l2_abab"], &l2_bbbb = evecs["l2_bbbb"];

            TA::TArrayD
                &s0 = evecs["s0"],
                &s1_aa = evecs["s1_aa"], &s1_bb = evecs["s1_bb"],
                &s2_aaaa = evecs["s2_aaaa"], &s2_abab = evecs["s2_abab"], &s2_bbbb = evecs["s2_bbbb"],
                &m0 = evecs["m0"],
                &m1_aa = evecs["m1_aa"], &m1_bb = evecs["m1_bb"],
                &m2_aaaa = evecs["m2_aaaa"], &m2_abab = evecs["m2_abab"], &m2_bbbb = evecs["m2_bbbb"];

            {
                TA::TArrayD tempOp_LF;
                TA::TArrayD tempOp_xaa_Loo;
                TA::TArrayD tempOp_xaa_Lvv;
                TA::TArrayD tempOp_xaaaa_Lvooo;
                TA::TArrayD tempOp_xaabb_Lvooo;
                TA::TArrayD tempOp_xbaab_Lvooo;
                TA::TArrayD tempOp_xbb_Loo;
                TA::TArrayD tempOp_xbb_Lvv;
                TA::TArrayD tempOp_xbbbb_Lvooo;
                TA::TArrayD tempOp_xxaa_LFoo;
                TA::TArrayD tempOp_xxaa_LFvo;
                TA::TArrayD tempOp_xxbb_LFoo;
                TA::TArrayD tempOp_xxbb_LFvo;

                /// ****** p†q ****** ///

                if (include_m2_ && include_s2_) {

                    // tempOp_LF += 1.000000 s2_abab("F,b,a,j,i") m2_abab("I,j,i,b,a") 
                    // flops: o2v2L1F1: 1, o0v0L1F1: 1 | mem: o0v0L1F1: 2, 
                    tempOp_LF("I,F") = s2_abab("F,b,a,j,i") * m2_abab("I,j,i,b,a");

                    // RDM_blks_["D1_aa_oo"] += 0.250000 d_aa(m,n) m2_abab(I,j,i,a,b) s2_abab(F,a,b,j,i) 
                    // +               0.250000 d_aa(m,n) m2_abab(I,j,i,b,a) s2_abab(F,b,a,j,i) 
                    // +               0.250000 d_aa(m,n) m2_abab(I,i,j,a,b) s2_abab(F,a,b,i,j) 
                    // +               0.250000 d_aa(m,n) m2_abab(I,i,j,b,a) s2_abab(F,b,a,i,j) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["aa_oo"]("m,n");

                    // RDM_blks_["D1_bb_oo"] += 0.250000 d_bb(m,n) m2_abab(I,j,i,a,b) s2_abab(F,a,b,j,i) 
                    // +               0.250000 d_bb(m,n) m2_abab(I,j,i,b,a) s2_abab(F,b,a,j,i) 
                    // +               0.250000 d_bb(m,n) m2_abab(I,i,j,a,b) s2_abab(F,a,b,i,j) 
                    // +               0.250000 d_bb(m,n) m2_abab(I,i,j,b,a) s2_abab(F,b,a,i,j) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["bb_oo"]("m,n");

                    // RDM_blks_["D1_aa_ov"] += 0.250000 t1_aa(e,m) m2_abab(I,j,i,a,b) s2_abab(F,a,b,j,i) 
                    // +               0.250000 t1_aa(e,m) m2_abab(I,j,i,b,a) s2_abab(F,b,a,j,i) 
                    // +               0.250000 t1_aa(e,m) m2_abab(I,i,j,a,b) s2_abab(F,a,b,i,j) 
                    // +               0.250000 t1_aa(e,m) m2_abab(I,i,j,b,a) s2_abab(F,b,a,i,j) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_aa("e,m");

                    // RDM_blks_["D1_bb_ov"] += 0.250000 t1_bb(e,m) m2_abab(I,j,i,a,b) s2_abab(F,a,b,j,i) 
                    // +               0.250000 t1_bb(e,m) m2_abab(I,j,i,b,a) s2_abab(F,b,a,j,i) 
                    // +               0.250000 t1_bb(e,m) m2_abab(I,i,j,a,b) s2_abab(F,a,b,i,j) 
                    // +               0.250000 t1_bb(e,m) m2_abab(I,i,j,b,a) s2_abab(F,b,a,i,j) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_bb("e,m");
                }

                if (include_m2_) {

                    // tempOp_xaa_Loo += 1.000000 t2_abab("b,a,m,i") m2_abab("I,n,i,b,a") 
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2, 
                    tempOp_xaa_Loo("I,m,n") = t2_abab("b,a,m,i") * m2_abab("I,n,i,b,a");
                }

                if (include_m2_ && include_s0_) {

                    // RDM_blks_["D1_aa_oo"] += -0.500000 t2_abab(b,a,m,i) m2_abab(I,n,i,b,a) s0(F) 
                    // +               -0.500000 t2_abab(a,b,m,i) m2_abab(I,n,i,a,b) s0(F) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempOp_xaa_Loo("I,m,n") * s0("F");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 t1_aa(e,j) t2_abab(b,a,m,i) m2_abab(I,j,i,b,a) s0(F) 
                    // +               -0.500000 t1_aa(e,j) t2_abab(a,b,m,i) m2_abab(I,j,i,a,b) s0(F) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t1_aa("e,j") * tempOp_xaa_Loo("I,m,j") * s0("F");
                }

                if (include_m2_ && include_s1_) {

                    // RDM_blks_["D1_aa_ov"] += -0.500000 t2_abab(b,a,m,i) m2_abab(I,j,i,b,a) s1_aa(F,e,j) 
                    // +               -0.500000 t2_abab(a,b,m,i) m2_abab(I,j,i,a,b) s1_aa(F,e,j) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xaa_Loo("I,m,j") * s1_aa("F,e,j");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_aa_ov"] += -0.500000 r0(F) t2_abab(b,a,m,i) u1_aa(e,j) m2_abab(I,j,i,b,a) 
                    // +               -0.500000 r0(F) t2_abab(a,b,m,i) u1_aa(e,j) m2_abab(I,j,i,a,b) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= u1_aa("e,j") * tempOp_xaa_Loo("I,m,j") * r0("F");
                }

                if (include_m2_ && include_u2_) {

                    // tempOp_xaa_Lvv += 1.000000 u2_abab("f,a,i,j") m2_abab("I,i,j,e,a") 
                    // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2, 
                    tempOp_xaa_Lvv("I,f,e") = u2_abab("f,a,i,j") * m2_abab("I,i,j,e,a");

                    // RDM_blks_["D1_aa_vv"] += 0.500000 r0(F) u2_abab(f,a,i,j) m2_abab(I,i,j,e,a) 
                    // +               0.500000 r0(F) u2_abab(f,a,j,i) m2_abab(I,j,i,e,a) 
                    // flops: o0v4L1F1: 1, o0v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += tempOp_xaa_Lvv("I,f,e") * r0("F");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 r1_aa(F,a,m) u2_abab(e,b,i,j) m2_abab(I,i,j,a,b) 
                    // +               -0.500000 r1_aa(F,a,m) u2_abab(e,b,j,i) m2_abab(I,j,i,a,b) 
                    // flops: o1v3L1F1: 1, o1v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xaa_Lvv("I,e,a") * r1_aa("F,a,m");
                }

                if (include_m2_) {

                    // tempOp_xaaaa_Lvooo += 1.000000 m2_aaaa("I,i,n,a,b") t1_aa("a,m") 
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2, 
                    tempOp_xaaaa_Lvooo("I,b,i,n,m") = m2_aaaa("I,i,n,a,b") * t1_aa("a,m");
                }

                if (include_m2_ && include_s0_) {

                    // RDM_blks_["D1_aa_ov"] += -0.500000 t1_aa(b,m) t2_aaaa(e,a,i,j) m2_aaaa(I,i,j,b,a) s0(F) 
                    // flops: o1v3L1F1: 1, o3v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= 0.500000 * t2_aaaa("e,a,i,j") * tempOp_xaaaa_Lvooo("I,a,i,j,m") * s0("F");
                }

                if (include_m2_ && include_s1_) {

                    // RDM_blks_["D1_aa_oo"] += 1.000000 t1_aa(a,m) m2_aaaa(I,i,n,a,b) s1_aa(F,b,i) 
                    // flops: o2v2L1F1: 1, o3v1L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempOp_xaaaa_Lvooo("I,b,i,n,m") * s1_aa("F,b,i");

                    // RDM_blks_["D1_aa_ov"] += 1.000000 t1_aa(e,i) t1_aa(a,m) m2_aaaa(I,j,i,a,b) s1_aa(F,b,j) 
                    // flops: o1v3L1F1: 1, o3v1L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += s1_aa("F,b,j") * tempOp_xaaaa_Lvooo("I,b,j,i,m") * t1_aa("e,i");
                }

                if (include_m2_ && include_s2_) {

                    // RDM_blks_["D1_aa_ov"] += -0.500000 t1_aa(a,m) m2_aaaa(I,j,i,a,b) s2_aaaa(F,e,b,j,i) 
                    // flops: o3v2L1F1: 1, o1v3L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= 0.500000 * tempOp_xaaaa_Lvooo("I,b,j,i,m") * s2_aaaa("F,e,b,j,i");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_aa_ov"] += 1.000000 r1_aa(F,a,i) t1_aa(b,m) u1_aa(e,j) m2_aaaa(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o3v1L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += r1_aa("F,a,i") * tempOp_xaaaa_Lvooo("I,a,i,j,m") * u1_aa("e,j");
                }

                if (include_m2_) {

                    // tempOp_xaabb_Lvooo += 1.000000 t1_bb("a,m") m2_abab("I,i,n,b,a") 
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2, 
                    tempOp_xaabb_Lvooo("I,b,i,m,n") = t1_bb("a,m") * m2_abab("I,i,n,b,a");
                }

                if (include_m2_ && include_s0_) {

                    // RDM_blks_["D1_bb_ov"] += -0.500000 t1_bb(b,m) t2_abab(a,e,i,j) m2_abab(I,i,j,a,b) s0(F) 
                    // +               -0.500000 t1_bb(b,m) t2_abab(a,e,j,i) m2_abab(I,j,i,a,b) s0(F) 
                    // flops: o1v3L1F1: 1, o3v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t2_abab("a,e,i,j") * tempOp_xaabb_Lvooo("I,a,i,m,j") * s0("F");
                }

                if (include_m2_ && include_s1_) {

                    // RDM_blks_["D1_bb_oo"] += -1.000000 t1_bb(a,m) m2_abab(I,i,n,b,a) s1_aa(F,b,i) 
                    // flops: o2v2L1F1: 1, o3v1L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempOp_xaabb_Lvooo("I,b,i,m,n") * s1_aa("F,b,i");

                    // RDM_blks_["D1_bb_ov"] += -1.000000 t1_bb(e,i) t1_bb(a,m) m2_abab(I,j,i,b,a) s1_aa(F,b,j) 
                    // flops: o1v3L1F1: 1, o3v1L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= s1_aa("F,b,j") * tempOp_xaabb_Lvooo("I,b,j,m,i") * t1_bb("e,i");
                }

                if (include_m2_ && include_s2_) {

                    // RDM_blks_["D1_bb_ov"] += -0.500000 t1_bb(a,m) m2_abab(I,j,i,b,a) s2_abab(F,b,e,j,i) 
                    // +               -0.500000 t1_bb(a,m) m2_abab(I,i,j,b,a) s2_abab(F,b,e,i,j) 
                    // flops: o3v2L1F1: 1, o1v3L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xaabb_Lvooo("I,b,j,m,i") * s2_abab("F,b,e,j,i");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_bb_ov"] += -1.000000 r1_aa(F,a,i) t1_bb(b,m) u1_bb(e,j) m2_abab(I,i,j,a,b) 
                    // flops: o1v3L1F1: 1, o3v1L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= r1_aa("F,a,i") * tempOp_xaabb_Lvooo("I,a,i,m,j") * u1_bb("e,j");
                }

                if (include_m2_ && include_u2_) {

                    // RDM_blks_["D1_bb_ov"] += -0.500000 r0(F) t1_bb(a,m) u2_abab(b,e,i,j) m2_abab(I,i,j,b,a) 
                    // +               -0.500000 r0(F) t1_bb(a,m) u2_abab(b,e,j,i) m2_abab(I,j,i,b,a) 
                    // flops: o1v3L1F1: 1, o3v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= u2_abab("b,e,i,j") * tempOp_xaabb_Lvooo("I,b,i,m,j") * r0("F");
                }

                if (include_m2_) {

                    // tempOp_xbaab_Lvooo += 1.000000 m2_abab("I,n,i,a,b") t1_aa("a,m") 
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2, 
                    tempOp_xbaab_Lvooo("I,b,n,m,i") = m2_abab("I,n,i,a,b") * t1_aa("a,m");
                }

                if (include_m2_ && include_s0_) {

                    // RDM_blks_["D1_aa_ov"] += -0.500000 t1_aa(b,m) t2_abab(e,a,i,j) m2_abab(I,i,j,b,a) s0(F) 
                    // +               -0.500000 t1_aa(b,m) t2_abab(e,a,j,i) m2_abab(I,j,i,b,a) s0(F) 
                    // flops: o1v3L1F1: 1, o3v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t2_abab("e,a,i,j") * tempOp_xbaab_Lvooo("I,a,i,m,j") * s0("F");
                }

                if (include_m2_ && include_s1_) {

                    // RDM_blks_["D1_aa_oo"] += -1.000000 t1_aa(a,m) m2_abab(I,n,i,a,b) s1_bb(F,b,i) 
                    // flops: o2v2L1F1: 1, o3v1L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempOp_xbaab_Lvooo("I,b,n,m,i") * s1_bb("F,b,i");

                    // RDM_blks_["D1_aa_ov"] += -1.000000 t1_aa(e,i) t1_aa(a,m) m2_abab(I,i,j,a,b) s1_bb(F,b,j) 
                    // flops: o1v3L1F1: 1, o3v1L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= s1_bb("F,b,j") * tempOp_xbaab_Lvooo("I,b,i,m,j") * t1_aa("e,i");
                }

                if (include_m2_ && include_s2_) {

                    // RDM_blks_["D1_aa_ov"] += -0.500000 t1_aa(a,m) m2_abab(I,j,i,a,b) s2_abab(F,e,b,j,i) 
                    // +               -0.500000 t1_aa(a,m) m2_abab(I,i,j,a,b) s2_abab(F,e,b,i,j) 
                    // flops: o3v2L1F1: 1, o1v3L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xbaab_Lvooo("I,b,j,m,i") * s2_abab("F,e,b,j,i");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_aa_ov"] += -1.000000 r1_bb(F,a,i) t1_aa(b,m) u1_aa(e,j) m2_abab(I,j,i,b,a) 
                    // flops: o1v3L1F1: 1, o3v1L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= r1_bb("F,a,i") * tempOp_xbaab_Lvooo("I,a,j,m,i") * u1_aa("e,j");
                }

                if (include_m2_ && include_u2_) {

                    // RDM_blks_["D1_aa_ov"] += -0.500000 r0(F) t1_aa(a,m) u2_abab(e,b,i,j) m2_abab(I,i,j,a,b) 
                    // +               -0.500000 r0(F) t1_aa(a,m) u2_abab(e,b,j,i) m2_abab(I,j,i,a,b) 
                    // flops: o1v3L1F1: 1, o3v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= u2_abab("e,b,i,j") * tempOp_xbaab_Lvooo("I,b,i,m,j") * r0("F");
                }

                if (include_m2_) {

                    // tempOp_xbb_Loo += 1.000000 m2_abab("I,i,n,b,a") t2_abab("b,a,i,m") 
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2, 
                    tempOp_xbb_Loo("I,n,m") = m2_abab("I,i,n,b,a") * t2_abab("b,a,i,m");
                }

                if (include_m2_ && include_s0_) {

                    // RDM_blks_["D1_bb_oo"] += -0.500000 t2_abab(b,a,i,m) m2_abab(I,i,n,b,a) s0(F) 
                    // +               -0.500000 t2_abab(a,b,i,m) m2_abab(I,i,n,a,b) s0(F) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempOp_xbb_Loo("I,n,m") * s0("F");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 t1_bb(e,j) t2_abab(b,a,i,m) m2_abab(I,i,j,b,a) s0(F) 
                    // +               -0.500000 t1_bb(e,j) t2_abab(a,b,i,m) m2_abab(I,i,j,a,b) s0(F) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t1_bb("e,j") * tempOp_xbb_Loo("I,j,m") * s0("F");
                }

                if (include_m2_ && include_s1_) {

                    // RDM_blks_["D1_bb_ov"] += -0.500000 t2_abab(b,a,i,m) m2_abab(I,i,j,b,a) s1_bb(F,e,j) 
                    // +               -0.500000 t2_abab(a,b,i,m) m2_abab(I,i,j,a,b) s1_bb(F,e,j) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xbb_Loo("I,j,m") * s1_bb("F,e,j");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_bb_ov"] += -0.500000 r0(F) t2_abab(b,a,i,m) u1_bb(e,j) m2_abab(I,i,j,b,a) 
                    // +               -0.500000 r0(F) t2_abab(a,b,i,m) u1_bb(e,j) m2_abab(I,i,j,a,b) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= u1_bb("e,j") * tempOp_xbb_Loo("I,j,m") * r0("F");
                }

                if (include_m2_) {

                    // tempOp_xbb_Lvv += 1.000000 m2_abab("I,i,j,a,e") t2_abab("a,f,i,j") 
                    // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2, 
                    tempOp_xbb_Lvv("I,e,f") = m2_abab("I,i,j,a,e") * t2_abab("a,f,i,j");
                }

                if (include_m2_ && include_s0_) {

                    // RDM_blks_["D1_bb_vv"] += 0.500000 t2_abab(a,f,i,j) m2_abab(I,i,j,a,e) s0(F) 
                    // +               0.500000 t2_abab(a,f,j,i) m2_abab(I,j,i,a,e) s0(F) 
                    // flops: o0v4L1F1: 1, o0v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += tempOp_xbb_Lvv("I,e,f") * s0("F");
                }

                if (include_m2_ && include_s1_) {

                    // RDM_blks_["D1_bb_ov"] += -0.500000 t2_abab(a,e,i,j) m2_abab(I,i,j,a,b) s1_bb(F,b,m) 
                    // +               -0.500000 t2_abab(a,e,j,i) m2_abab(I,j,i,a,b) s1_bb(F,b,m) 
                    // flops: o1v3L1F1: 1, o1v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xbb_Lvv("I,b,e") * s1_bb("F,b,m");
                }

                if (include_m2_) {

                    // tempOp_xbbbb_Lvooo += 1.000000 m2_bbbb("I,i,n,a,b") t1_bb("a,m") 
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2, 
                    tempOp_xbbbb_Lvooo("I,b,i,n,m") = m2_bbbb("I,i,n,a,b") * t1_bb("a,m");
                }

                if (include_m2_ && include_s0_) {

                    // RDM_blks_["D1_bb_ov"] += -0.500000 t1_bb(b,m) t2_bbbb(e,a,i,j) m2_bbbb(I,i,j,b,a) s0(F) 
                    // flops: o1v3L1F1: 1, o3v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= 0.500000 * t2_bbbb("e,a,i,j") * tempOp_xbbbb_Lvooo("I,a,i,j,m") * s0("F");
                }

                if (include_m2_ && include_s1_) {

                    // RDM_blks_["D1_bb_oo"] += 1.000000 t1_bb(a,m) m2_bbbb(I,i,n,a,b) s1_bb(F,b,i) 
                    // flops: o2v2L1F1: 1, o3v1L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempOp_xbbbb_Lvooo("I,b,i,n,m") * s1_bb("F,b,i");

                    // RDM_blks_["D1_bb_ov"] += 1.000000 t1_bb(e,i) t1_bb(a,m) m2_bbbb(I,j,i,a,b) s1_bb(F,b,j) 
                    // flops: o1v3L1F1: 1, o3v1L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += s1_bb("F,b,j") * tempOp_xbbbb_Lvooo("I,b,j,i,m") * t1_bb("e,i");
                }

                if (include_m2_ && include_s2_) {

                    // RDM_blks_["D1_bb_ov"] += -0.500000 t1_bb(a,m) m2_bbbb(I,j,i,a,b) s2_bbbb(F,e,b,j,i) 
                    // flops: o3v2L1F1: 1, o1v3L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= 0.500000 * tempOp_xbbbb_Lvooo("I,b,j,i,m") * s2_bbbb("F,e,b,j,i");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_bb_ov"] += 1.000000 r1_bb(F,a,i) t1_bb(b,m) u1_bb(e,j) m2_bbbb(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o3v1L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += r1_bb("F,a,i") * tempOp_xbbbb_Lvooo("I,a,i,j,m") * u1_bb("e,j");
                }

                if (include_m2_ && include_s2_) {

                    // tempOp_xxaa_LFoo += 0.500000 s2_aaaa("F,a,b,i,m") m2_aaaa("I,i,n,a,b") 
                    // flops: o3v2L1F1: 1, o2v0L1F1: 1 | mem: o2v0L1F1: 2, 
                    tempOp_xxaa_LFoo("I,F,m,n") = 0.500000 * s2_aaaa("F,a,b,i,m") * m2_aaaa("I,i,n,a,b");

                    // RDM_blks_["D1_aa_oo"] += -0.500000 m2_aaaa(I,i,n,a,b) s2_aaaa(F,a,b,i,m) 
                    // flops: o2v2L1F1: 1 | mem: o2v2L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempOp_xxaa_LFoo("I,F,m,n");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 t1_aa(e,i) m2_aaaa(I,j,i,a,b) s2_aaaa(F,a,b,j,m) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xxaa_LFoo("I,F,m,i") * t1_aa("e,i");
                }

                if (include_m2_) {

                    // tempOp_xxaa_LFvo += 1.000000 m2_aaaa("I,i,j,e,a") r1_aa("F,a,i") 
                    // flops: o2v2L1F1: 1, o1v1L1F1: 1 | mem: o1v1L1F1: 2, 
                    tempOp_xxaa_LFvo("I,F,e,j") = m2_aaaa("I,i,j,e,a") * r1_aa("F,a,i");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_aa_vv"] += -1.000000 r1_aa(F,a,i) u1_aa(f,j) m2_aaaa(I,i,j,e,a) 
                    // flops: o0v4L1F1: 1, o1v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") -= tempOp_xxaa_LFvo("I,F,e,j") * u1_aa("f,j");
                }

                if (include_m2_ && include_u2_) {

                    // RDM_blks_["D1_aa_ov"] += 1.000000 r1_aa(F,a,i) u2_aaaa(e,b,j,m) m2_aaaa(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_xxaa_LFvo("I,F,b,j") * u2_aaaa("e,b,j,m");

                    // RDM_blks_["D1_bb_ov"] += -1.000000 r1_aa(F,a,i) u2_abab(b,e,j,m) m2_aaaa(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xxaa_LFvo("I,F,b,j") * u2_abab("b,e,j,m");
                }

                {

                    // tempOp_xxbb_LFoo += 0.500000 r2_bbbb("F,a,b,i,m") l2_bbbb("I,i,n,a,b") 
                    // flops: o3v2L1F1: 1, o2v0L1F1: 1 | mem: o2v0L1F1: 2, 
                    tempOp_xxbb_LFoo("I,F,m,n") = 0.500000 * r2_bbbb("F,a,b,i,m") * l2_bbbb("I,i,n,a,b");

                    // RDM_blks_["D1_bb_oo"] += -0.500000 l2_bbbb(I,i,n,a,b) r2_bbbb(F,a,b,i,m) 
                    // flops: o2v2L1F1: 1 | mem: o2v2L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempOp_xxbb_LFoo("I,F,m,n");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 l2_bbbb(I,j,i,a,b) r2_bbbb(F,a,b,j,m) t1_bb(e,i) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xxbb_LFoo("I,F,m,i") * t1_bb("e,i");
                }

                if (include_m2_) {

                    // tempOp_xxbb_LFvo += 1.000000 r1_bb("F,a,i") m2_bbbb("I,i,j,e,a") 
                    // flops: o2v2L1F1: 1, o1v1L1F1: 1 | mem: o1v1L1F1: 2, 
                    tempOp_xxbb_LFvo("I,F,e,j") = r1_bb("F,a,i") * m2_bbbb("I,i,j,e,a");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_bb_vv"] += -1.000000 r1_bb(F,a,i) u1_bb(f,j) m2_bbbb(I,i,j,e,a) 
                    // flops: o0v4L1F1: 1, o1v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") -= tempOp_xxbb_LFvo("I,F,e,j") * u1_bb("f,j");
                }

                if (include_m2_ && include_u2_) {

                    // RDM_blks_["D1_aa_ov"] += -1.000000 r1_bb(F,a,i) u2_abab(e,b,m,j) m2_bbbb(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xxbb_LFvo("I,F,b,j") * u2_abab("e,b,m,j");

                    // RDM_blks_["D1_bb_ov"] += 1.000000 r1_bb(F,a,i) u2_bbbb(e,b,j,m) m2_bbbb(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_xxbb_LFvo("I,F,b,j") * u2_bbbb("e,b,j,m");
                }

                if (include_m2_ && include_s2_) {

                    // tempOp_LF += 0.250000 s2_aaaa("F,a,b,j,i") m2_aaaa("I,j,i,a,b") 
                    // flops: o2v2L1F1: 1, o0v0L1F1: 1 | mem: o0v0L1F1: 2, 
                    tempOp_LF("I,F") = 0.250000 * s2_aaaa("F,a,b,j,i") * m2_aaaa("I,j,i,a,b");

                    // RDM_blks_["D1_aa_oo"] += 0.250000 d_aa(m,n) m2_aaaa(I,j,i,a,b) s2_aaaa(F,a,b,j,i) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["aa_oo"]("m,n");

                    // RDM_blks_["D1_bb_oo"] += 0.250000 d_bb(m,n) m2_aaaa(I,j,i,a,b) s2_aaaa(F,a,b,j,i) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["bb_oo"]("m,n");

                    // RDM_blks_["D1_aa_ov"] += 0.250000 t1_aa(e,m) m2_aaaa(I,j,i,a,b) s2_aaaa(F,a,b,j,i) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_aa("e,m");

                    // RDM_blks_["D1_bb_ov"] += 0.250000 t1_bb(e,m) m2_aaaa(I,j,i,a,b) s2_aaaa(F,a,b,j,i) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_bb("e,m");
                }

                if (include_m2_) {

                    // tempOp_xaa_Loo += 0.500000 m2_aaaa("I,i,n,b,a") t2_aaaa("b,a,i,m") 
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2, 
                    tempOp_xaa_Loo("I,n,m") = 0.500000 * m2_aaaa("I,i,n,b,a") * t2_aaaa("b,a,i,m");
                }

                if (include_m2_ && include_s0_) {

                    // RDM_blks_["D1_aa_oo"] += -0.500000 t2_aaaa(b,a,i,m) m2_aaaa(I,i,n,b,a) s0(F) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempOp_xaa_Loo("I,n,m") * s0("F");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 t1_aa(e,j) t2_aaaa(b,a,i,m) m2_aaaa(I,i,j,b,a) s0(F) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t1_aa("e,j") * tempOp_xaa_Loo("I,j,m") * s0("F");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_aa_ov"] += -0.500000 r0(F) t2_aaaa(b,a,i,m) u1_aa(e,j) m2_aaaa(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= u1_aa("e,j") * tempOp_xaa_Loo("I,j,m") * r0("F");
                }

                {

                    // tempOp_xaa_Lvv += 1.000000 l2_abab("I,i,j,e,a") t2_abab("f,a,i,j") 
                    // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2, 
                    tempOp_xaa_Lvv("I,e,f") = l2_abab("I,i,j,e,a") * t2_abab("f,a,i,j");

                    // RDM_blks_["D1_aa_vv"] += 0.500000 l2_abab(I,i,j,e,a) r0(F) t2_abab(f,a,i,j) 
                    // +               0.500000 l2_abab(I,j,i,e,a) r0(F) t2_abab(f,a,j,i) 
                    // flops: o0v4L1F1: 1, o0v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += tempOp_xaa_Lvv("I,e,f") * r0("F");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 l2_abab(I,i,j,b,a) r1_aa(F,b,m) t2_abab(e,a,i,j) 
                    // +               -0.500000 l2_abab(I,j,i,b,a) r1_aa(F,b,m) t2_abab(e,a,j,i) 
                    // flops: o1v3L1F1: 1, o1v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xaa_Lvv("I,b,e") * r1_aa("F,b,m");

                    // tempOp_xaaaa_Lvooo += 1.000000 l2_aaaa("I,i,n,a,b") t1_aa("a,m") 
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2, 
                    tempOp_xaaaa_Lvooo("I,b,i,n,m") = l2_aaaa("I,i,n,a,b") * t1_aa("a,m");

                    // RDM_blks_["D1_aa_oo"] += 1.000000 l2_aaaa(I,i,n,a,b) r1_aa(F,b,i) t1_aa(a,m) 
                    // flops: o2v2L1F1: 1, o3v1L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempOp_xaaaa_Lvooo("I,b,i,n,m") * r1_aa("F,b,i");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 l2_aaaa(I,j,i,a,b) r2_aaaa(F,e,b,j,i) t1_aa(a,m) 
                    // flops: o3v2L1F1: 1, o1v3L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= 0.500000 * tempOp_xaaaa_Lvooo("I,b,j,i,m") * r2_aaaa("F,e,b,j,i");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 l2_aaaa(I,i,j,b,a) r0(F) t1_aa(b,m) t2_aaaa(e,a,i,j) 
                    // flops: o1v3L1F1: 1, o3v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= 0.500000 * t2_aaaa("e,a,i,j") * tempOp_xaaaa_Lvooo("I,a,i,j,m") * r0("F");

                    // RDM_blks_["D1_aa_ov"] += 1.000000 l2_aaaa(I,j,i,a,b) r1_aa(F,b,j) t1_aa(e,i) t1_aa(a,m) 
                    // flops: o1v3L1F1: 1, o3v1L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += r1_aa("F,b,j") * tempOp_xaaaa_Lvooo("I,b,j,i,m") * t1_aa("e,i");

                    // tempOp_xaabb_Lvooo += 1.000000 t1_bb("a,m") l2_abab("I,i,n,b,a") 
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2, 
                    tempOp_xaabb_Lvooo("I,b,i,m,n") = t1_bb("a,m") * l2_abab("I,i,n,b,a");

                    // RDM_blks_["D1_bb_oo"] += -1.000000 l2_abab(I,i,n,b,a) r1_aa(F,b,i) t1_bb(a,m) 
                    // flops: o2v2L1F1: 1, o3v1L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempOp_xaabb_Lvooo("I,b,i,m,n") * r1_aa("F,b,i");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 l2_abab(I,i,j,a,b) r0(F) t1_bb(b,m) t2_abab(a,e,i,j) 
                    // +               -0.500000 l2_abab(I,j,i,a,b) r0(F) t1_bb(b,m) t2_abab(a,e,j,i) 
                    // flops: o1v3L1F1: 1, o3v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t2_abab("a,e,i,j") * tempOp_xaabb_Lvooo("I,a,i,m,j") * r0("F");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 l2_abab(I,j,i,b,a) r2_abab(F,b,e,j,i) t1_bb(a,m) 
                    // +               -0.500000 l2_abab(I,i,j,b,a) r2_abab(F,b,e,i,j) t1_bb(a,m) 
                    // flops: o3v2L1F1: 1, o1v3L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xaabb_Lvooo("I,b,j,m,i") * r2_abab("F,b,e,j,i");

                    // RDM_blks_["D1_bb_ov"] += -1.000000 l2_abab(I,j,i,b,a) r1_aa(F,b,j) t1_bb(e,i) t1_bb(a,m) 
                    // flops: o1v3L1F1: 1, o3v1L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= r1_aa("F,b,j") * tempOp_xaabb_Lvooo("I,b,j,m,i") * t1_bb("e,i");
                }

                if (include_m2_ && include_u1_) {

                    // tempOp_xbaab_Lvooo += 1.000000 u1_aa("b,m") m2_abab("I,n,i,b,a") 
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2, 
                    tempOp_xbaab_Lvooo("I,a,m,n,i") = u1_aa("b,m") * m2_abab("I,n,i,b,a");

                    // RDM_blks_["D1_aa_oo"] += -1.000000 r1_bb(F,a,i) u1_aa(b,m) m2_abab(I,n,i,b,a) 
                    // flops: o2v2L1F1: 1, o3v1L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempOp_xbaab_Lvooo("I,a,m,n,i") * r1_bb("F,a,i");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 r2_abab(F,e,a,i,j) u1_aa(b,m) m2_abab(I,i,j,b,a) 
                    // +               -0.500000 r2_abab(F,e,a,j,i) u1_aa(b,m) m2_abab(I,j,i,b,a) 
                    // flops: o3v2L1F1: 1, o1v3L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xbaab_Lvooo("I,a,m,i,j") * r2_abab("F,e,a,i,j");

                    // RDM_blks_["D1_aa_ov"] += -1.000000 r1_bb(F,a,i) t1_aa(e,j) u1_aa(b,m) m2_abab(I,j,i,b,a) 
                    // flops: o1v3L1F1: 1, o3v1L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= r1_bb("F,a,i") * tempOp_xbaab_Lvooo("I,a,m,j,i") * t1_aa("e,j");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 r0(F) t2_abab(e,a,i,j) u1_aa(b,m) m2_abab(I,i,j,b,a) 
                    // +               -0.500000 r0(F) t2_abab(e,a,j,i) u1_aa(b,m) m2_abab(I,j,i,b,a) 
                    // flops: o1v3L1F1: 1, o3v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t2_abab("e,a,i,j") * tempOp_xbaab_Lvooo("I,a,m,i,j") * r0("F");
                }

                if (include_m2_) {

                    // tempOp_xbb_Loo += 0.500000 t2_bbbb("b,a,i,m") m2_bbbb("I,i,n,b,a") 
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2, 
                    tempOp_xbb_Loo("I,m,n") = 0.500000 * t2_bbbb("b,a,i,m") * m2_bbbb("I,i,n,b,a");
                }

                if (include_m2_ && include_s0_) {

                    // RDM_blks_["D1_bb_oo"] += -0.500000 t2_bbbb(b,a,i,m) m2_bbbb(I,i,n,b,a) s0(F) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempOp_xbb_Loo("I,m,n") * s0("F");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 t1_bb(e,j) t2_bbbb(b,a,i,m) m2_bbbb(I,i,j,b,a) s0(F) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t1_bb("e,j") * tempOp_xbb_Loo("I,m,j") * s0("F");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_bb_ov"] += -0.500000 r0(F) t2_bbbb(b,a,i,m) u1_bb(e,j) m2_bbbb(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= u1_bb("e,j") * tempOp_xbb_Loo("I,m,j") * r0("F");
                }

                {

                    // tempOp_xbb_Lvv += 1.000000 t2_abab("a,f,i,j") l2_abab("I,i,j,a,e") 
                    // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2, 
                    tempOp_xbb_Lvv("I,f,e") = t2_abab("a,f,i,j") * l2_abab("I,i,j,a,e");

                    // RDM_blks_["D1_bb_vv"] += 0.500000 l2_abab(I,i,j,a,e) r0(F) t2_abab(a,f,i,j) 
                    // +               0.500000 l2_abab(I,j,i,a,e) r0(F) t2_abab(a,f,j,i) 
                    // flops: o0v4L1F1: 1, o0v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += tempOp_xbb_Lvv("I,f,e") * r0("F");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 l2_abab(I,i,j,a,b) r1_bb(F,b,m) t2_abab(a,e,i,j) 
                    // +               -0.500000 l2_abab(I,j,i,a,b) r1_bb(F,b,m) t2_abab(a,e,j,i) 
                    // flops: o1v3L1F1: 1, o1v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xbb_Lvv("I,e,b") * r1_bb("F,b,m");
                }

                if (include_m2_ && include_u1_) {

                    // tempOp_xbbbb_Lvooo += 1.000000 u1_bb("b,m") m2_bbbb("I,i,n,b,a") 
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2, 
                    tempOp_xbbbb_Lvooo("I,a,m,i,n") = u1_bb("b,m") * m2_bbbb("I,i,n,b,a");

                    // RDM_blks_["D1_bb_oo"] += 1.000000 r1_bb(F,a,i) u1_bb(b,m) m2_bbbb(I,i,n,b,a) 
                    // flops: o2v2L1F1: 1, o3v1L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempOp_xbbbb_Lvooo("I,a,m,i,n") * r1_bb("F,a,i");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 r0(F) t2_bbbb(e,a,i,j) u1_bb(b,m) m2_bbbb(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o3v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= 0.500000 * t2_bbbb("e,a,i,j") * tempOp_xbbbb_Lvooo("I,a,m,i,j") * r0("F");

                    // RDM_blks_["D1_bb_ov"] += 1.000000 r1_bb(F,a,i) t1_bb(e,j) u1_bb(b,m) m2_bbbb(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o3v1L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += r1_bb("F,a,i") * tempOp_xbbbb_Lvooo("I,a,m,i,j") * t1_bb("e,j");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 r2_bbbb(F,e,a,i,j) u1_bb(b,m) m2_bbbb(I,i,j,b,a) 
                    // flops: o3v2L1F1: 1, o1v3L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= 0.500000 * tempOp_xbbbb_Lvooo("I,a,m,i,j") * r2_bbbb("F,e,a,i,j");
                }

                {

                    // tempOp_xxaa_LFoo += 1.000000 r2_abab("F,a,b,m,i") l2_abab("I,n,i,a,b") 
                    // flops: o3v2L1F1: 1, o2v0L1F1: 1 | mem: o2v0L1F1: 2, 
                    tempOp_xxaa_LFoo("I,F,m,n") = r2_abab("F,a,b,m,i") * l2_abab("I,n,i,a,b");

                    // RDM_blks_["D1_aa_oo"] += -0.500000 l2_abab(I,n,i,a,b) r2_abab(F,a,b,m,i) 
                    // +               -0.500000 l2_abab(I,n,i,b,a) r2_abab(F,b,a,m,i) 
                    // flops: o2v2L1F1: 1 | mem: o2v2L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempOp_xxaa_LFoo("I,F,m,n");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 l2_abab(I,i,j,a,b) r2_abab(F,a,b,m,j) t1_aa(e,i) 
                    // +               -0.500000 l2_abab(I,i,j,b,a) r2_abab(F,b,a,m,j) t1_aa(e,i) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xxaa_LFoo("I,F,m,i") * t1_aa("e,i");

                    // tempOp_xxaa_LFvo += 1.000000 l2_abab("I,i,j,e,a") r1_bb("F,a,j") 
                    // flops: o2v2L1F1: 1, o1v1L1F1: 1 | mem: o1v1L1F1: 2, 
                    tempOp_xxaa_LFvo("I,F,e,i") = l2_abab("I,i,j,e,a") * r1_bb("F,a,j");

                    // RDM_blks_["D1_aa_vv"] += 1.000000 l2_abab(I,i,j,e,a) r1_bb(F,a,j) t1_aa(f,i) 
                    // flops: o0v4L1F1: 1, o1v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += tempOp_xxaa_LFvo("I,F,e,i") * t1_aa("f,i");

                    // RDM_blks_["D1_aa_ov"] += -1.000000 l2_abab(I,i,j,a,b) r1_bb(F,b,j) t2_aaaa(e,a,i,m) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xxaa_LFvo("I,F,a,i") * t2_aaaa("e,a,i,m");

                    // RDM_blks_["D1_bb_ov"] += 1.000000 l2_abab(I,i,j,a,b) r1_bb(F,b,j) t2_abab(a,e,i,m) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_xxaa_LFvo("I,F,a,i") * t2_abab("a,e,i,m");

                    // RDM_blks_["D1_aa_vo"] += 1.000000 l2_abab(I,m,i,e,a) r1_bb(F,a,i) 
                    // flops: o1v3L1F1: 1 | mem: o1v3L1F1: 1, 
                    RDM_blks_["D1_aa_vo"]("I,F,e,m") += tempOp_xxaa_LFvo("I,F,e,m");
                }

                if (include_m2_ && include_s2_) {

                    // tempOp_xxbb_LFoo += 0.500000 m2_bbbb("I,i,n,a,b") s2_bbbb("F,a,b,i,m") 
                    // flops: o3v2L1F1: 1, o2v0L1F1: 1 | mem: o2v0L1F1: 2, 
                    tempOp_xxbb_LFoo("I,F,n,m") = 0.500000 * m2_bbbb("I,i,n,a,b") * s2_bbbb("F,a,b,i,m");

                    // RDM_blks_["D1_bb_oo"] += -0.500000 m2_bbbb(I,i,n,a,b) s2_bbbb(F,a,b,i,m) 
                    // flops: o2v2L1F1: 1 | mem: o2v2L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempOp_xxbb_LFoo("I,F,n,m");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 t1_bb(e,i) m2_bbbb(I,j,i,a,b) s2_bbbb(F,a,b,j,m) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xxbb_LFoo("I,F,i,m") * t1_bb("e,i");
                }

                if (include_m2_) {

                    // tempOp_xxbb_LFvo += 1.000000 m2_abab("I,i,j,a,e") r1_aa("F,a,i") 
                    // flops: o2v2L1F1: 1, o1v1L1F1: 1 | mem: o1v1L1F1: 2, 
                    tempOp_xxbb_LFvo("I,F,e,j") = m2_abab("I,i,j,a,e") * r1_aa("F,a,i");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_bb_vv"] += 1.000000 r1_aa(F,a,i) u1_bb(f,j) m2_abab(I,i,j,a,e) 
                    // flops: o0v4L1F1: 1, o1v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += tempOp_xxbb_LFvo("I,F,e,j") * u1_bb("f,j");
                }

                if (include_m2_ && include_u2_) {

                    // RDM_blks_["D1_aa_ov"] += 1.000000 r1_aa(F,a,i) u2_abab(e,b,m,j) m2_abab(I,i,j,a,b) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_xxbb_LFvo("I,F,b,j") * u2_abab("e,b,m,j");

                    // RDM_blks_["D1_bb_ov"] += -1.000000 r1_aa(F,a,i) u2_bbbb(e,b,j,m) m2_abab(I,i,j,a,b) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xxbb_LFvo("I,F,b,j") * u2_bbbb("e,b,j,m");
                }

                if (include_m2_ && include_s2_) {

                    // tempOp_LF += 0.250000 s2_bbbb("F,a,b,j,i") m2_bbbb("I,j,i,a,b") 
                    // flops: o2v2L1F1: 1, o0v0L1F1: 1 | mem: o0v0L1F1: 2, 
                    tempOp_LF("I,F") = 0.250000 * s2_bbbb("F,a,b,j,i") * m2_bbbb("I,j,i,a,b");

                    // RDM_blks_["D1_aa_oo"] += 0.250000 d_aa(m,n) m2_bbbb(I,j,i,a,b) s2_bbbb(F,a,b,j,i) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["aa_oo"]("m,n");

                    // RDM_blks_["D1_bb_oo"] += 0.250000 d_bb(m,n) m2_bbbb(I,j,i,a,b) s2_bbbb(F,a,b,j,i) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["bb_oo"]("m,n");

                    // RDM_blks_["D1_aa_ov"] += 0.250000 t1_aa(e,m) m2_bbbb(I,j,i,a,b) s2_bbbb(F,a,b,j,i) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_aa("e,m");

                    // RDM_blks_["D1_bb_ov"] += 0.250000 t1_bb(e,m) m2_bbbb(I,j,i,a,b) s2_bbbb(F,a,b,j,i) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_bb("e,m");
                }

                {

                    // tempOp_xaa_Loo += 1.000000 l2_abab("I,n,i,b,a") t2_abab("b,a,m,i") 
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2, 
                    tempOp_xaa_Loo("I,n,m") = l2_abab("I,n,i,b,a") * t2_abab("b,a,m,i");

                    // RDM_blks_["D1_aa_oo"] += -0.500000 l2_abab(I,n,i,b,a) r0(F) t2_abab(b,a,m,i) 
                    // +               -0.500000 l2_abab(I,n,i,a,b) r0(F) t2_abab(a,b,m,i) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempOp_xaa_Loo("I,n,m") * r0("F");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 l2_abab(I,j,i,b,a) r0(F) t1_aa(e,j) t2_abab(b,a,m,i) 
                    // +               -0.500000 l2_abab(I,j,i,a,b) r0(F) t1_aa(e,j) t2_abab(a,b,m,i) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t1_aa("e,j") * tempOp_xaa_Loo("I,j,m") * r0("F");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 l2_abab(I,j,i,b,a) r1_aa(F,e,j) t2_abab(b,a,m,i) 
                    // +               -0.500000 l2_abab(I,j,i,a,b) r1_aa(F,e,j) t2_abab(a,b,m,i) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xaa_Loo("I,j,m") * r1_aa("F,e,j");
                }

                if (include_m2_) {

                    // tempOp_xaa_Lvv += 1.000000 m2_abab("I,i,j,e,a") t2_abab("f,a,i,j") 
                    // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2, 
                    tempOp_xaa_Lvv("I,e,f") = m2_abab("I,i,j,e,a") * t2_abab("f,a,i,j");
                }

                if (include_m2_ && include_s0_) {

                    // RDM_blks_["D1_aa_vv"] += 0.500000 t2_abab(f,a,i,j) m2_abab(I,i,j,e,a) s0(F) 
                    // +               0.500000 t2_abab(f,a,j,i) m2_abab(I,j,i,e,a) s0(F) 
                    // flops: o0v4L1F1: 1, o0v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += tempOp_xaa_Lvv("I,e,f") * s0("F");
                }

                if (include_m2_ && include_s1_) {

                    // RDM_blks_["D1_aa_ov"] += -0.500000 t2_abab(e,a,i,j) m2_abab(I,i,j,b,a) s1_aa(F,b,m) 
                    // +               -0.500000 t2_abab(e,a,j,i) m2_abab(I,j,i,b,a) s1_aa(F,b,m) 
                    // flops: o1v3L1F1: 1, o1v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xaa_Lvv("I,b,e") * s1_aa("F,b,m");
                }

                if (include_m2_ && include_u1_) {

                    // tempOp_xaaaa_Lvooo += 1.000000 m2_aaaa("I,i,n,b,a") u1_aa("b,m") 
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2, 
                    tempOp_xaaaa_Lvooo("I,a,i,n,m") = m2_aaaa("I,i,n,b,a") * u1_aa("b,m");

                    // RDM_blks_["D1_aa_oo"] += 1.000000 r1_aa(F,a,i) u1_aa(b,m) m2_aaaa(I,i,n,b,a) 
                    // flops: o2v2L1F1: 1, o3v1L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempOp_xaaaa_Lvooo("I,a,i,n,m") * r1_aa("F,a,i");

                    // RDM_blks_["D1_aa_ov"] += 1.000000 r1_aa(F,a,i) t1_aa(e,j) u1_aa(b,m) m2_aaaa(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o3v1L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += r1_aa("F,a,i") * tempOp_xaaaa_Lvooo("I,a,i,j,m") * t1_aa("e,j");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 r0(F) t2_aaaa(e,a,i,j) u1_aa(b,m) m2_aaaa(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o3v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= 0.500000 * t2_aaaa("e,a,i,j") * tempOp_xaaaa_Lvooo("I,a,i,j,m") * r0("F");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 r2_aaaa(F,e,a,i,j) u1_aa(b,m) m2_aaaa(I,i,j,b,a) 
                    // flops: o3v2L1F1: 1, o1v3L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= 0.500000 * tempOp_xaaaa_Lvooo("I,a,i,j,m") * r2_aaaa("F,e,a,i,j");

                    // tempOp_xaabb_Lvooo += 1.000000 u1_bb("b,m") m2_abab("I,i,n,a,b") 
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2, 
                    tempOp_xaabb_Lvooo("I,a,i,m,n") = u1_bb("b,m") * m2_abab("I,i,n,a,b");

                    // RDM_blks_["D1_bb_oo"] += -1.000000 r1_aa(F,a,i) u1_bb(b,m) m2_abab(I,i,n,a,b) 
                    // flops: o2v2L1F1: 1, o3v1L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempOp_xaabb_Lvooo("I,a,i,m,n") * r1_aa("F,a,i");

                    // RDM_blks_["D1_bb_ov"] += -1.000000 r1_aa(F,a,i) t1_bb(e,j) u1_bb(b,m) m2_abab(I,i,j,a,b) 
                    // flops: o1v3L1F1: 1, o3v1L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= r1_aa("F,a,i") * tempOp_xaabb_Lvooo("I,a,i,m,j") * t1_bb("e,j");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 r2_abab(F,a,e,i,j) u1_bb(b,m) m2_abab(I,i,j,a,b) 
                    // +               -0.500000 r2_abab(F,a,e,j,i) u1_bb(b,m) m2_abab(I,j,i,a,b) 
                    // flops: o3v2L1F1: 1, o1v3L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xaabb_Lvooo("I,a,i,m,j") * r2_abab("F,a,e,i,j");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 r0(F) t2_abab(a,e,i,j) u1_bb(b,m) m2_abab(I,i,j,a,b) 
                    // +               -0.500000 r0(F) t2_abab(a,e,j,i) u1_bb(b,m) m2_abab(I,j,i,a,b) 
                    // flops: o1v3L1F1: 1, o3v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t2_abab("a,e,i,j") * tempOp_xaabb_Lvooo("I,a,i,m,j") * r0("F");
                }

                {

                    // tempOp_xbaab_Lvooo += 1.000000 l2_abab("I,n,i,a,b") t1_aa("a,m") 
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2, 
                    tempOp_xbaab_Lvooo("I,b,n,m,i") = l2_abab("I,n,i,a,b") * t1_aa("a,m");

                    // RDM_blks_["D1_aa_oo"] += -1.000000 l2_abab(I,n,i,a,b) r1_bb(F,b,i) t1_aa(a,m) 
                    // flops: o2v2L1F1: 1, o3v1L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempOp_xbaab_Lvooo("I,b,n,m,i") * r1_bb("F,b,i");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 l2_abab(I,i,j,b,a) r0(F) t1_aa(b,m) t2_abab(e,a,i,j) 
                    // +               -0.500000 l2_abab(I,j,i,b,a) r0(F) t1_aa(b,m) t2_abab(e,a,j,i) 
                    // flops: o1v3L1F1: 1, o3v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t2_abab("e,a,i,j") * tempOp_xbaab_Lvooo("I,a,i,m,j") * r0("F");

                    // RDM_blks_["D1_aa_ov"] += -1.000000 l2_abab(I,i,j,a,b) r1_bb(F,b,j) t1_aa(e,i) t1_aa(a,m) 
                    // flops: o1v3L1F1: 1, o3v1L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= r1_bb("F,b,j") * tempOp_xbaab_Lvooo("I,b,i,m,j") * t1_aa("e,i");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 l2_abab(I,j,i,a,b) r2_abab(F,e,b,j,i) t1_aa(a,m) 
                    // +               -0.500000 l2_abab(I,i,j,a,b) r2_abab(F,e,b,i,j) t1_aa(a,m) 
                    // flops: o3v2L1F1: 1, o1v3L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xbaab_Lvooo("I,b,j,m,i") * r2_abab("F,e,b,j,i");

                    // tempOp_xbb_Loo += 1.000000 t2_abab("b,a,i,m") l2_abab("I,i,n,b,a") 
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2, 
                    tempOp_xbb_Loo("I,m,n") = t2_abab("b,a,i,m") * l2_abab("I,i,n,b,a");

                    // RDM_blks_["D1_bb_oo"] += -0.500000 l2_abab(I,i,n,b,a) r0(F) t2_abab(b,a,i,m) 
                    // +               -0.500000 l2_abab(I,i,n,a,b) r0(F) t2_abab(a,b,i,m) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempOp_xbb_Loo("I,m,n") * r0("F");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 l2_abab(I,i,j,b,a) r1_bb(F,e,j) t2_abab(b,a,i,m) 
                    // +               -0.500000 l2_abab(I,i,j,a,b) r1_bb(F,e,j) t2_abab(a,b,i,m) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xbb_Loo("I,m,j") * r1_bb("F,e,j");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 l2_abab(I,i,j,b,a) r0(F) t1_bb(e,j) t2_abab(b,a,i,m) 
                    // +               -0.500000 l2_abab(I,i,j,a,b) r0(F) t1_bb(e,j) t2_abab(a,b,i,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t1_bb("e,j") * tempOp_xbb_Loo("I,m,j") * r0("F");
                }

                if (include_m2_ && include_u2_) {

                    // tempOp_xbb_Lvv += 1.000000 u2_abab("a,f,i,j") m2_abab("I,i,j,a,e") 
                    // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2, 
                    tempOp_xbb_Lvv("I,f,e") = u2_abab("a,f,i,j") * m2_abab("I,i,j,a,e");

                    // RDM_blks_["D1_bb_vv"] += 0.500000 r0(F) u2_abab(a,f,i,j) m2_abab(I,i,j,a,e) 
                    // +               0.500000 r0(F) u2_abab(a,f,j,i) m2_abab(I,j,i,a,e) 
                    // flops: o0v4L1F1: 1, o0v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += tempOp_xbb_Lvv("I,f,e") * r0("F");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 r1_bb(F,a,m) u2_abab(b,e,i,j) m2_abab(I,i,j,b,a) 
                    // +               -0.500000 r1_bb(F,a,m) u2_abab(b,e,j,i) m2_abab(I,j,i,b,a) 
                    // flops: o1v3L1F1: 1, o1v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xbb_Lvv("I,e,a") * r1_bb("F,a,m");
                }

                {

                    // tempOp_xbbbb_Lvooo += 1.000000 l2_bbbb("I,i,n,a,b") t1_bb("a,m") 
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2, 
                    tempOp_xbbbb_Lvooo("I,b,i,n,m") = l2_bbbb("I,i,n,a,b") * t1_bb("a,m");

                    // RDM_blks_["D1_bb_oo"] += 1.000000 l2_bbbb(I,i,n,a,b) r1_bb(F,b,i) t1_bb(a,m) 
                    // flops: o2v2L1F1: 1, o3v1L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempOp_xbbbb_Lvooo("I,b,i,n,m") * r1_bb("F,b,i");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 l2_bbbb(I,j,i,a,b) r2_bbbb(F,e,b,j,i) t1_bb(a,m) 
                    // flops: o3v2L1F1: 1, o1v3L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= 0.500000 * tempOp_xbbbb_Lvooo("I,b,j,i,m") * r2_bbbb("F,e,b,j,i");

                    // RDM_blks_["D1_bb_ov"] += 1.000000 l2_bbbb(I,j,i,a,b) r1_bb(F,b,j) t1_bb(e,i) t1_bb(a,m) 
                    // flops: o1v3L1F1: 1, o3v1L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += r1_bb("F,b,j") * tempOp_xbbbb_Lvooo("I,b,j,i,m") * t1_bb("e,i");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 l2_bbbb(I,i,j,b,a) r0(F) t1_bb(b,m) t2_bbbb(e,a,i,j) 
                    // flops: o1v3L1F1: 1, o3v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= 0.500000 * t2_bbbb("e,a,i,j") * tempOp_xbbbb_Lvooo("I,a,i,j,m") * r0("F");
                }

                if (include_m2_ && include_s2_) {

                    // tempOp_xxaa_LFoo += 1.000000 s2_abab("F,a,b,m,i") m2_abab("I,n,i,a,b") 
                    // flops: o3v2L1F1: 1, o2v0L1F1: 1 | mem: o2v0L1F1: 2, 
                    tempOp_xxaa_LFoo("I,F,m,n") = s2_abab("F,a,b,m,i") * m2_abab("I,n,i,a,b");

                    // RDM_blks_["D1_aa_oo"] += -0.500000 m2_abab(I,n,i,a,b) s2_abab(F,a,b,m,i) 
                    // +               -0.500000 m2_abab(I,n,i,b,a) s2_abab(F,b,a,m,i) 
                    // flops: o2v2L1F1: 1 | mem: o2v2L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempOp_xxaa_LFoo("I,F,m,n");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 t1_aa(e,i) m2_abab(I,i,j,a,b) s2_abab(F,a,b,m,j) 
                    // +               -0.500000 t1_aa(e,i) m2_abab(I,i,j,b,a) s2_abab(F,b,a,m,j) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xxaa_LFoo("I,F,m,i") * t1_aa("e,i");
                }

                if (include_m2_ && include_s1_) {

                    // tempOp_xxaa_LFvo += 1.000000 s1_aa("F,a,j") m2_aaaa("I,j,i,e,a") 
                    // flops: o2v2L1F1: 1, o1v1L1F1: 1 | mem: o1v1L1F1: 2, 
                    tempOp_xxaa_LFvo("I,F,e,i") = s1_aa("F,a,j") * m2_aaaa("I,j,i,e,a");

                    // RDM_blks_["D1_aa_vv"] += -1.000000 t1_aa(f,i) m2_aaaa(I,j,i,e,a) s1_aa(F,a,j) 
                    // flops: o0v4L1F1: 1, o1v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") -= tempOp_xxaa_LFvo("I,F,e,i") * t1_aa("f,i");

                    // RDM_blks_["D1_aa_ov"] += 1.000000 t2_aaaa(e,a,i,m) m2_aaaa(I,j,i,a,b) s1_aa(F,b,j) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_xxaa_LFvo("I,F,a,i") * t2_aaaa("e,a,i,m");

                    // RDM_blks_["D1_bb_ov"] += -1.000000 t2_abab(a,e,i,m) m2_aaaa(I,j,i,a,b) s1_aa(F,b,j) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xxaa_LFvo("I,F,a,i") * t2_abab("a,e,i,m");

                    // RDM_blks_["D1_aa_vo"] += -1.000000 m2_aaaa(I,i,m,e,a) s1_aa(F,a,i) 
                    // flops: o1v3L1F1: 1 | mem: o1v3L1F1: 1, 
                    RDM_blks_["D1_aa_vo"]("I,F,e,m") -= tempOp_xxaa_LFvo("I,F,e,m");
                }

                {

                    // tempOp_xxbb_LFoo += 1.000000 l2_abab("I,i,n,a,b") r2_abab("F,a,b,i,m") 
                    // flops: o3v2L1F1: 1, o2v0L1F1: 1 | mem: o2v0L1F1: 2, 
                    tempOp_xxbb_LFoo("I,F,n,m") = l2_abab("I,i,n,a,b") * r2_abab("F,a,b,i,m");

                    // RDM_blks_["D1_bb_oo"] += -0.500000 l2_abab(I,i,n,a,b) r2_abab(F,a,b,i,m) 
                    // +               -0.500000 l2_abab(I,i,n,b,a) r2_abab(F,b,a,i,m) 
                    // flops: o2v2L1F1: 1 | mem: o2v2L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempOp_xxbb_LFoo("I,F,n,m");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 l2_abab(I,j,i,a,b) r2_abab(F,a,b,j,m) t1_bb(e,i) 
                    // +               -0.500000 l2_abab(I,j,i,b,a) r2_abab(F,b,a,j,m) t1_bb(e,i) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xxbb_LFoo("I,F,i,m") * t1_bb("e,i");
                }

                if (include_m2_ && include_s1_) {

                    // tempOp_xxbb_LFvo += 1.000000 s1_bb("F,a,j") m2_bbbb("I,j,i,e,a") 
                    // flops: o2v2L1F1: 1, o1v1L1F1: 1 | mem: o1v1L1F1: 2, 
                    tempOp_xxbb_LFvo("I,F,e,i") = s1_bb("F,a,j") * m2_bbbb("I,j,i,e,a");

                    // RDM_blks_["D1_bb_vv"] += -1.000000 t1_bb(f,i) m2_bbbb(I,j,i,e,a) s1_bb(F,a,j) 
                    // flops: o0v4L1F1: 1, o1v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") -= tempOp_xxbb_LFvo("I,F,e,i") * t1_bb("f,i");

                    // RDM_blks_["D1_aa_ov"] += -1.000000 t2_abab(e,a,m,i) m2_bbbb(I,j,i,a,b) s1_bb(F,b,j) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xxbb_LFvo("I,F,a,i") * t2_abab("e,a,m,i");

                    // RDM_blks_["D1_bb_ov"] += 1.000000 t2_bbbb(e,a,i,m) m2_bbbb(I,j,i,a,b) s1_bb(F,b,j) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_xxbb_LFvo("I,F,a,i") * t2_bbbb("e,a,i,m");

                    // RDM_blks_["D1_bb_vo"] += -1.000000 m2_bbbb(I,i,m,e,a) s1_bb(F,a,i) 
                    // flops: o1v3L1F1: 1 | mem: o1v3L1F1: 1, 
                    RDM_blks_["D1_bb_vo"]("I,F,e,m") -= tempOp_xxbb_LFvo("I,F,e,m");
                }

                {

                    // tempOp_LF += 0.250000 l2_aaaa("I,j,i,a,b") r2_aaaa("F,a,b,j,i") 
                    // flops: o2v2L1F1: 1, o0v0L1F1: 1 | mem: o0v0L1F1: 2, 
                    tempOp_LF("I,F") = 0.250000 * l2_aaaa("I,j,i,a,b") * r2_aaaa("F,a,b,j,i");

                    // RDM_blks_["D1_aa_oo"] += 0.250000 d_aa(m,n) l2_aaaa(I,j,i,a,b) r2_aaaa(F,a,b,j,i) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["aa_oo"]("m,n");

                    // RDM_blks_["D1_bb_oo"] += 0.250000 d_bb(m,n) l2_aaaa(I,j,i,a,b) r2_aaaa(F,a,b,j,i) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["bb_oo"]("m,n");

                    // RDM_blks_["D1_aa_ov"] += 0.250000 l2_aaaa(I,j,i,a,b) r2_aaaa(F,a,b,j,i) t1_aa(e,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_aa("e,m");

                    // RDM_blks_["D1_bb_ov"] += 0.250000 l2_aaaa(I,j,i,a,b) r2_aaaa(F,a,b,j,i) t1_bb(e,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_bb("e,m");
                }

                if (include_m2_ && include_u2_) {

                    // tempOp_xaa_Loo += 1.000000 m2_abab("I,n,i,b,a") u2_abab("b,a,m,i") 
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2, 
                    tempOp_xaa_Loo("I,n,m") = m2_abab("I,n,i,b,a") * u2_abab("b,a,m,i");

                    // RDM_blks_["D1_aa_oo"] += -0.500000 r0(F) u2_abab(b,a,m,i) m2_abab(I,n,i,b,a) 
                    // +               -0.500000 r0(F) u2_abab(a,b,m,i) m2_abab(I,n,i,a,b) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempOp_xaa_Loo("I,n,m") * r0("F");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 r1_aa(F,e,i) u2_abab(b,a,m,j) m2_abab(I,i,j,b,a) 
                    // +               -0.500000 r1_aa(F,e,i) u2_abab(a,b,m,j) m2_abab(I,i,j,a,b) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xaa_Loo("I,i,m") * r1_aa("F,e,i");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 r0(F) t1_aa(e,i) u2_abab(b,a,m,j) m2_abab(I,i,j,b,a) 
                    // +               -0.500000 r0(F) t1_aa(e,i) u2_abab(a,b,m,j) m2_abab(I,i,j,a,b) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t1_aa("e,i") * tempOp_xaa_Loo("I,i,m") * r0("F");

                    // tempOp_xbb_Loo += 1.000000 m2_abab("I,i,n,b,a") u2_abab("b,a,i,m") 
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2, 
                    tempOp_xbb_Loo("I,n,m") = m2_abab("I,i,n,b,a") * u2_abab("b,a,i,m");

                    // RDM_blks_["D1_bb_oo"] += -0.500000 r0(F) u2_abab(b,a,i,m) m2_abab(I,i,n,b,a) 
                    // +               -0.500000 r0(F) u2_abab(a,b,i,m) m2_abab(I,i,n,a,b) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempOp_xbb_Loo("I,n,m") * r0("F");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 r1_bb(F,e,i) u2_abab(b,a,j,m) m2_abab(I,j,i,b,a) 
                    // +               -0.500000 r1_bb(F,e,i) u2_abab(a,b,j,m) m2_abab(I,j,i,a,b) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xbb_Loo("I,i,m") * r1_bb("F,e,i");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 r0(F) t1_bb(e,i) u2_abab(b,a,j,m) m2_abab(I,j,i,b,a) 
                    // +               -0.500000 r0(F) t1_bb(e,i) u2_abab(a,b,j,m) m2_abab(I,j,i,a,b) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t1_bb("e,i") * tempOp_xbb_Loo("I,i,m") * r0("F");
                }

                {

                    // tempOp_xxaa_LFoo += 0.500000 r2_aaaa("F,a,b,i,m") l2_aaaa("I,i,n,a,b") 
                    // flops: o3v2L1F1: 1, o2v0L1F1: 1 | mem: o2v0L1F1: 2, 
                    tempOp_xxaa_LFoo("I,F,m,n") = 0.500000 * r2_aaaa("F,a,b,i,m") * l2_aaaa("I,i,n,a,b");

                    // RDM_blks_["D1_aa_oo"] += -0.500000 l2_aaaa(I,i,n,a,b) r2_aaaa(F,a,b,i,m) 
                    // flops: o2v2L1F1: 1 | mem: o2v2L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempOp_xxaa_LFoo("I,F,m,n");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 l2_aaaa(I,j,i,a,b) r2_aaaa(F,a,b,j,m) t1_aa(e,i) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xxaa_LFoo("I,F,m,i") * t1_aa("e,i");
                }

                if (include_m2_ && include_s1_) {

                    // tempOp_xxaa_LFvo += 1.000000 s1_bb("F,a,j") m2_abab("I,i,j,e,a") 
                    // flops: o2v2L1F1: 1, o1v1L1F1: 1 | mem: o1v1L1F1: 2, 
                    tempOp_xxaa_LFvo("I,F,e,i") = s1_bb("F,a,j") * m2_abab("I,i,j,e,a");

                    // RDM_blks_["D1_aa_vv"] += 1.000000 t1_aa(f,i) m2_abab(I,i,j,e,a) s1_bb(F,a,j) 
                    // flops: o0v4L1F1: 1, o1v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += tempOp_xxaa_LFvo("I,F,e,i") * t1_aa("f,i");

                    // RDM_blks_["D1_aa_ov"] += -1.000000 t2_aaaa(e,a,i,m) m2_abab(I,i,j,a,b) s1_bb(F,b,j) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xxaa_LFvo("I,F,a,i") * t2_aaaa("e,a,i,m");

                    // RDM_blks_["D1_bb_ov"] += 1.000000 t2_abab(a,e,i,m) m2_abab(I,i,j,a,b) s1_bb(F,b,j) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_xxaa_LFvo("I,F,a,i") * t2_abab("a,e,i,m");

                    // RDM_blks_["D1_aa_vo"] += 1.000000 m2_abab(I,m,i,e,a) s1_bb(F,a,i) 
                    // flops: o1v3L1F1: 1 | mem: o1v3L1F1: 1, 
                    RDM_blks_["D1_aa_vo"]("I,F,e,m") += tempOp_xxaa_LFvo("I,F,e,m");
                }

                if (include_m2_ && include_s2_) {

                    // tempOp_xxbb_LFoo += 1.000000 m2_abab("I,i,n,a,b") s2_abab("F,a,b,i,m") 
                    // flops: o3v2L1F1: 1, o2v0L1F1: 1 | mem: o2v0L1F1: 2, 
                    tempOp_xxbb_LFoo("I,F,n,m") = m2_abab("I,i,n,a,b") * s2_abab("F,a,b,i,m");

                    // RDM_blks_["D1_bb_oo"] += -0.500000 m2_abab(I,i,n,a,b) s2_abab(F,a,b,i,m) 
                    // +               -0.500000 m2_abab(I,i,n,b,a) s2_abab(F,b,a,i,m) 
                    // flops: o2v2L1F1: 1 | mem: o2v2L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempOp_xxbb_LFoo("I,F,n,m");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 t1_bb(e,i) m2_abab(I,j,i,a,b) s2_abab(F,a,b,j,m) 
                    // +               -0.500000 t1_bb(e,i) m2_abab(I,j,i,b,a) s2_abab(F,b,a,j,m) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xxbb_LFoo("I,F,i,m") * t1_bb("e,i");
                }

                {

                    // tempOp_xxbb_LFvo += 1.000000 r1_bb("F,a,j") l2_bbbb("I,j,i,e,a") 
                    // flops: o2v2L1F1: 1, o1v1L1F1: 1 | mem: o1v1L1F1: 2, 
                    tempOp_xxbb_LFvo("I,F,e,i") = r1_bb("F,a,j") * l2_bbbb("I,j,i,e,a");

                    // RDM_blks_["D1_bb_vv"] += -1.000000 l2_bbbb(I,j,i,e,a) r1_bb(F,a,j) t1_bb(f,i) 
                    // flops: o0v4L1F1: 1, o1v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") -= tempOp_xxbb_LFvo("I,F,e,i") * t1_bb("f,i");

                    // RDM_blks_["D1_aa_ov"] += -1.000000 l2_bbbb(I,j,i,a,b) r1_bb(F,b,j) t2_abab(e,a,m,i) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xxbb_LFvo("I,F,a,i") * t2_abab("e,a,m,i");

                    // RDM_blks_["D1_bb_ov"] += 1.000000 l2_bbbb(I,j,i,a,b) r1_bb(F,b,j) t2_bbbb(e,a,i,m) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_xxbb_LFvo("I,F,a,i") * t2_bbbb("e,a,i,m");

                    // RDM_blks_["D1_bb_vo"] += -1.000000 l2_bbbb(I,i,m,e,a) r1_bb(F,a,i) 
                    // flops: o1v3L1F1: 1 | mem: o1v3L1F1: 1, 
                    RDM_blks_["D1_bb_vo"]("I,F,e,m") -= tempOp_xxbb_LFvo("I,F,e,m");

                    // tempOp_LF += 1.000000 l2_abab("I,i,j,b,a") r2_abab("F,b,a,i,j") 
                    // flops: o2v2L1F1: 1, o0v0L1F1: 1 | mem: o0v0L1F1: 2, 
                    tempOp_LF("I,F") = l2_abab("I,i,j,b,a") * r2_abab("F,b,a,i,j");

                    // RDM_blks_["D1_aa_oo"] += 0.250000 d_aa(m,n) l2_abab(I,j,i,a,b) r2_abab(F,a,b,j,i) 
                    // +               0.250000 d_aa(m,n) l2_abab(I,j,i,b,a) r2_abab(F,b,a,j,i) 
                    // +               0.250000 d_aa(m,n) l2_abab(I,i,j,a,b) r2_abab(F,a,b,i,j) 
                    // +               0.250000 d_aa(m,n) l2_abab(I,i,j,b,a) r2_abab(F,b,a,i,j) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["aa_oo"]("m,n");

                    // RDM_blks_["D1_bb_oo"] += 0.250000 d_bb(m,n) l2_abab(I,j,i,a,b) r2_abab(F,a,b,j,i) 
                    // +               0.250000 d_bb(m,n) l2_abab(I,j,i,b,a) r2_abab(F,b,a,j,i) 
                    // +               0.250000 d_bb(m,n) l2_abab(I,i,j,a,b) r2_abab(F,a,b,i,j) 
                    // +               0.250000 d_bb(m,n) l2_abab(I,i,j,b,a) r2_abab(F,b,a,i,j) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["bb_oo"]("m,n");

                    // RDM_blks_["D1_aa_ov"] += 0.250000 l2_abab(I,j,i,a,b) r2_abab(F,a,b,j,i) t1_aa(e,m) 
                    // +               0.250000 l2_abab(I,j,i,b,a) r2_abab(F,b,a,j,i) t1_aa(e,m) 
                    // +               0.250000 l2_abab(I,i,j,a,b) r2_abab(F,a,b,i,j) t1_aa(e,m) 
                    // +               0.250000 l2_abab(I,i,j,b,a) r2_abab(F,b,a,i,j) t1_aa(e,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_aa("e,m");

                    // RDM_blks_["D1_bb_ov"] += 0.250000 l2_abab(I,j,i,a,b) r2_abab(F,a,b,j,i) t1_bb(e,m) 
                    // +               0.250000 l2_abab(I,j,i,b,a) r2_abab(F,b,a,j,i) t1_bb(e,m) 
                    // +               0.250000 l2_abab(I,i,j,a,b) r2_abab(F,a,b,i,j) t1_bb(e,m) 
                    // +               0.250000 l2_abab(I,i,j,b,a) r2_abab(F,b,a,i,j) t1_bb(e,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_bb("e,m");

                    // tempOp_xaa_Loo += 0.500000 t2_aaaa("b,a,i,m") l2_aaaa("I,i,n,b,a") 
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2, 
                    tempOp_xaa_Loo("I,m,n") = 0.500000 * t2_aaaa("b,a,i,m") * l2_aaaa("I,i,n,b,a");

                    // RDM_blks_["D1_aa_oo"] += -0.500000 l2_aaaa(I,i,n,b,a) r0(F) t2_aaaa(b,a,i,m) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempOp_xaa_Loo("I,m,n") * r0("F");

                    // RDM_blks_["D1_aa_ov"] += -0.500000 l2_aaaa(I,i,j,b,a) r0(F) t1_aa(e,j) t2_aaaa(b,a,i,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t1_aa("e,j") * tempOp_xaa_Loo("I,m,j") * r0("F");

                    // tempOp_xbb_Loo += 0.500000 t2_bbbb("b,a,i,m") l2_bbbb("I,i,n,b,a") 
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2, 
                    tempOp_xbb_Loo("I,m,n") = 0.500000 * t2_bbbb("b,a,i,m") * l2_bbbb("I,i,n,b,a");

                    // RDM_blks_["D1_bb_oo"] += -0.500000 l2_bbbb(I,i,n,b,a) r0(F) t2_bbbb(b,a,i,m) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempOp_xbb_Loo("I,m,n") * r0("F");

                    // RDM_blks_["D1_bb_ov"] += -0.500000 l2_bbbb(I,i,j,b,a) r0(F) t1_bb(e,j) t2_bbbb(b,a,i,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t1_bb("e,j") * tempOp_xbb_Loo("I,m,j") * r0("F");

                    // tempOp_xxaa_LFvo += 1.000000 r1_aa("F,a,j") l2_aaaa("I,j,i,e,a") 
                    // flops: o2v2L1F1: 1, o1v1L1F1: 1 | mem: o1v1L1F1: 2, 
                    tempOp_xxaa_LFvo("I,F,e,i") = r1_aa("F,a,j") * l2_aaaa("I,j,i,e,a");

                    // RDM_blks_["D1_aa_vv"] += -1.000000 l2_aaaa(I,j,i,e,a) r1_aa(F,a,j) t1_aa(f,i) 
                    // flops: o0v4L1F1: 1, o1v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") -= tempOp_xxaa_LFvo("I,F,e,i") * t1_aa("f,i");

                    // RDM_blks_["D1_aa_ov"] += 1.000000 l2_aaaa(I,j,i,a,b) r1_aa(F,b,j) t2_aaaa(e,a,i,m) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_xxaa_LFvo("I,F,a,i") * t2_aaaa("e,a,i,m");

                    // RDM_blks_["D1_bb_ov"] += -1.000000 l2_aaaa(I,j,i,a,b) r1_aa(F,b,j) t2_abab(a,e,i,m) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xxaa_LFvo("I,F,a,i") * t2_abab("a,e,i,m");

                    // RDM_blks_["D1_aa_vo"] += -1.000000 l2_aaaa(I,i,m,e,a) r1_aa(F,a,i) 
                    // flops: o1v3L1F1: 1 | mem: o1v3L1F1: 1, 
                    RDM_blks_["D1_aa_vo"]("I,F,e,m") -= tempOp_xxaa_LFvo("I,F,e,m");

                    // tempOp_xxbb_LFvo += 1.000000 r1_aa("F,a,j") l2_abab("I,j,i,a,e") 
                    // flops: o2v2L1F1: 1, o1v1L1F1: 1 | mem: o1v1L1F1: 2, 
                    tempOp_xxbb_LFvo("I,F,e,i") = r1_aa("F,a,j") * l2_abab("I,j,i,a,e");

                    // RDM_blks_["D1_bb_vv"] += 1.000000 l2_abab(I,j,i,a,e) r1_aa(F,a,j) t1_bb(f,i) 
                    // flops: o0v4L1F1: 1, o1v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += tempOp_xxbb_LFvo("I,F,e,i") * t1_bb("f,i");

                    // RDM_blks_["D1_aa_ov"] += 1.000000 l2_abab(I,j,i,b,a) r1_aa(F,b,j) t2_abab(e,a,m,i) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_xxbb_LFvo("I,F,a,i") * t2_abab("e,a,m,i");

                    // RDM_blks_["D1_bb_ov"] += -1.000000 l2_abab(I,j,i,b,a) r1_aa(F,b,j) t2_bbbb(e,a,i,m) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xxbb_LFvo("I,F,a,i") * t2_bbbb("e,a,i,m");

                    // RDM_blks_["D1_bb_vo"] += 1.000000 l2_abab(I,i,m,a,e) r1_aa(F,a,i) 
                    // flops: o1v3L1F1: 1 | mem: o1v3L1F1: 1, 
                    RDM_blks_["D1_bb_vo"]("I,F,e,m") += tempOp_xxbb_LFvo("I,F,e,m");

                    // tempOp_LF += 0.250000 r2_bbbb("F,a,b,j,i") l2_bbbb("I,j,i,a,b") 
                    // flops: o2v2L1F1: 1, o0v0L1F1: 1 | mem: o0v0L1F1: 2, 
                    tempOp_LF("I,F") = 0.250000 * r2_bbbb("F,a,b,j,i") * l2_bbbb("I,j,i,a,b");

                    // RDM_blks_["D1_aa_oo"] += 0.250000 d_aa(m,n) l2_bbbb(I,j,i,a,b) r2_bbbb(F,a,b,j,i) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["aa_oo"]("m,n");

                    // RDM_blks_["D1_bb_oo"] += 0.250000 d_bb(m,n) l2_bbbb(I,j,i,a,b) r2_bbbb(F,a,b,j,i) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["bb_oo"]("m,n");

                    // RDM_blks_["D1_aa_ov"] += 0.250000 l2_bbbb(I,j,i,a,b) r2_bbbb(F,a,b,j,i) t1_aa(e,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_aa("e,m");

                    // RDM_blks_["D1_bb_ov"] += 0.250000 l2_bbbb(I,j,i,a,b) r2_bbbb(F,a,b,j,i) t1_bb(e,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_bb("e,m");
                }

                if (include_m2_ && include_u2_) {

                    // tempOp_xaa_Loo += 0.500000 u2_aaaa("b,a,j,m") m2_aaaa("I,i,j,b,a") 
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2, 
                    tempOp_xaa_Loo("I,m,i") = 0.500000 * u2_aaaa("b,a,j,m") * m2_aaaa("I,i,j,b,a");

                    // RDM_blks_["D1_aa_ov"] += 0.500000 r1_aa(F,e,i) u2_aaaa(b,a,j,m) m2_aaaa(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_xaa_Loo("I,m,i") * r1_aa("F,e,i");

                    // RDM_blks_["D1_aa_ov"] += 0.500000 r0(F) t1_aa(e,i) u2_aaaa(b,a,j,m) m2_aaaa(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += t1_aa("e,i") * tempOp_xaa_Loo("I,m,i") * r0("F");

                    // tempOp_xbb_Loo += 0.500000 u2_bbbb("b,a,j,m") m2_bbbb("I,i,j,b,a") 
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2, 
                    tempOp_xbb_Loo("I,m,i") = 0.500000 * u2_bbbb("b,a,j,m") * m2_bbbb("I,i,j,b,a");

                    // RDM_blks_["D1_bb_ov"] += 0.500000 r0(F) t1_bb(e,i) u2_bbbb(b,a,j,m) m2_bbbb(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += t1_bb("e,i") * tempOp_xbb_Loo("I,m,i") * r0("F");

                    // RDM_blks_["D1_bb_ov"] += 0.500000 r1_bb(F,e,i) u2_bbbb(b,a,j,m) m2_bbbb(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_xbb_Loo("I,m,i") * r1_bb("F,e,i");
                }

                if (include_m2_) {

                    // tempOp_xxaa_LFvo += 1.000000 r1_bb("F,a,i") m2_abab("I,j,i,e,a") 
                    // flops: o2v2L1F1: 1, o1v1L1F1: 1 | mem: o1v1L1F1: 2, 
                    tempOp_xxaa_LFvo("I,F,e,j") = r1_bb("F,a,i") * m2_abab("I,j,i,e,a");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_aa_vv"] += 1.000000 r1_bb(F,a,i) u1_aa(f,j) m2_abab(I,j,i,e,a) 
                    // flops: o0v4L1F1: 1, o1v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += tempOp_xxaa_LFvo("I,F,e,j") * u1_aa("f,j");
                }

                if (include_m2_ && include_u2_) {

                    // RDM_blks_["D1_aa_ov"] += -1.000000 r1_bb(F,a,i) u2_aaaa(e,b,j,m) m2_abab(I,j,i,b,a) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xxaa_LFvo("I,F,b,j") * u2_aaaa("e,b,j,m");

                    // RDM_blks_["D1_bb_ov"] += 1.000000 r1_bb(F,a,i) u2_abab(b,e,j,m) m2_abab(I,j,i,b,a) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_xxaa_LFvo("I,F,b,j") * u2_abab("b,e,j,m");
                }

                if (include_m2_ && include_s1_) {

                    // tempOp_xxbb_LFvo += 1.000000 s1_aa("F,a,j") m2_abab("I,j,i,a,e") 
                    // flops: o2v2L1F1: 1, o1v1L1F1: 1 | mem: o1v1L1F1: 2, 
                    tempOp_xxbb_LFvo("I,F,e,i") = s1_aa("F,a,j") * m2_abab("I,j,i,a,e");

                    // RDM_blks_["D1_bb_vv"] += 1.000000 t1_bb(f,i) m2_abab(I,j,i,a,e) s1_aa(F,a,j) 
                    // flops: o0v4L1F1: 1, o1v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += tempOp_xxbb_LFvo("I,F,e,i") * t1_bb("f,i");

                    // RDM_blks_["D1_aa_ov"] += 1.000000 t2_abab(e,a,m,i) m2_abab(I,j,i,b,a) s1_aa(F,b,j) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_xxbb_LFvo("I,F,a,i") * t2_abab("e,a,m,i");

                    // RDM_blks_["D1_bb_ov"] += -1.000000 t2_bbbb(e,a,i,m) m2_abab(I,j,i,b,a) s1_aa(F,b,j) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xxbb_LFvo("I,F,a,i") * t2_bbbb("e,a,i,m");

                    // RDM_blks_["D1_bb_vo"] += 1.000000 m2_abab(I,i,m,a,e) s1_aa(F,a,i) 
                    // flops: o1v3L1F1: 1 | mem: o1v3L1F1: 1, 
                    RDM_blks_["D1_bb_vo"]("I,F,e,m") += tempOp_xxbb_LFvo("I,F,e,m");
                }

                if (include_m2_) {

                    // tempOp_LF += 0.250000 m2_bbbb("I,i,j,b,a") r2_bbbb("F,b,a,i,j") 
                    // flops: o2v2L1F1: 1, o0v0L1F1: 1 | mem: o0v0L1F1: 2, 
                    tempOp_LF("I,F") = 0.250000 * m2_bbbb("I,i,j,b,a") * r2_bbbb("F,b,a,i,j");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_aa_ov"] += 0.250000 r2_bbbb(F,b,a,i,j) u1_aa(e,m) m2_bbbb(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_LF("I,F") * u1_aa("e,m");

                    // RDM_blks_["D1_bb_ov"] += 0.250000 r2_bbbb(F,b,a,i,j) u1_bb(e,m) m2_bbbb(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_LF("I,F") * u1_bb("e,m");
                }

                if (include_m1_) {

                    // tempOp_xaa_Loo += 1.000000 m1_aa("I,n,a") t1_aa("a,m") 
                    // flops: o2v1L1: 1, o2v0L1: 1 | mem: o2v0L1: 2, 
                    tempOp_xaa_Loo("I,n,m") = m1_aa("I,n,a") * t1_aa("a,m");
                }

                if (include_m1_ && include_s0_) {

                    // RDM_blks_["D1_aa_oo"] += -1.000000 t1_aa(a,m) m1_aa(I,n,a) s0(F) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempOp_xaa_Loo("I,n,m") * s0("F");

                    // RDM_blks_["D1_aa_ov"] += -1.000000 t1_aa(e,i) t1_aa(a,m) m1_aa(I,i,a) s0(F) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t1_aa("e,i") * tempOp_xaa_Loo("I,i,m") * s0("F");
                }

                if (include_m1_ && include_s1_) {

                    // RDM_blks_["D1_aa_ov"] += -1.000000 t1_aa(a,m) m1_aa(I,i,a) s1_aa(F,e,i) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempOp_xaa_Loo("I,i,m") * s1_aa("F,e,i");
                }

                if (include_m1_ && include_u1_) {

                    // RDM_blks_["D1_aa_ov"] += -1.000000 r0(F) t1_aa(a,m) u1_aa(e,i) m1_aa(I,i,a) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= u1_aa("e,i") * tempOp_xaa_Loo("I,i,m") * r0("F");
                }

                if (include_m1_) {

                    // tempOp_xbb_Loo += 1.000000 m1_bb("I,n,a") t1_bb("a,m") 
                    // flops: o2v1L1: 1, o2v0L1: 1 | mem: o2v0L1: 2, 
                    tempOp_xbb_Loo("I,n,m") = m1_bb("I,n,a") * t1_bb("a,m");
                }

                if (include_m1_ && include_s0_) {

                    // RDM_blks_["D1_bb_oo"] += -1.000000 t1_bb(a,m) m1_bb(I,n,a) s0(F) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempOp_xbb_Loo("I,n,m") * s0("F");

                    // RDM_blks_["D1_bb_ov"] += -1.000000 t1_bb(e,i) t1_bb(a,m) m1_bb(I,i,a) s0(F) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t1_bb("e,i") * tempOp_xbb_Loo("I,i,m") * s0("F");
                }

                if (include_m1_ && include_s1_) {

                    // RDM_blks_["D1_bb_ov"] += -1.000000 t1_bb(a,m) m1_bb(I,i,a) s1_bb(F,e,i) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempOp_xbb_Loo("I,i,m") * s1_bb("F,e,i");
                }

                if (include_m1_ && include_u1_) {

                    // RDM_blks_["D1_bb_ov"] += -1.000000 r0(F) t1_bb(a,m) u1_bb(e,i) m1_bb(I,i,a) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= u1_bb("e,i") * tempOp_xbb_Loo("I,i,m") * r0("F");
                }

                if (include_m2_) {

                    // tempOp_LF += 0.250000 m2_aaaa("I,i,j,b,a") r2_aaaa("F,b,a,i,j") 
                    // flops: o2v2L1F1: 1, o0v0L1F1: 1 | mem: o0v0L1F1: 2, 
                    tempOp_LF("I,F") = 0.250000 * m2_aaaa("I,i,j,b,a") * r2_aaaa("F,b,a,i,j");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_aa_ov"] += 0.250000 r2_aaaa(F,b,a,i,j) u1_aa(e,m) m2_aaaa(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_LF("I,F") * u1_aa("e,m");

                    // RDM_blks_["D1_bb_ov"] += 0.250000 r2_aaaa(F,b,a,i,j) u1_bb(e,m) m2_aaaa(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_LF("I,F") * u1_bb("e,m");
                }

                if (include_m2_) {

                    // tempOp_LF += 1.000000 r2_abab("F,a,b,j,i") m2_abab("I,j,i,a,b") 
                    // flops: o2v2L1F1: 1, o0v0L1F1: 1 | mem: o0v0L1F1: 2, 
                    tempOp_LF("I,F") = r2_abab("F,a,b,j,i") * m2_abab("I,j,i,a,b");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_aa_ov"] += 0.250000 r2_abab(F,b,a,i,j) u1_aa(e,m) m2_abab(I,i,j,b,a) 
                    // +               0.250000 r2_abab(F,b,a,j,i) u1_aa(e,m) m2_abab(I,j,i,b,a) 
                    // +               0.250000 r2_abab(F,a,b,i,j) u1_aa(e,m) m2_abab(I,i,j,a,b) 
                    // +               0.250000 r2_abab(F,a,b,j,i) u1_aa(e,m) m2_abab(I,j,i,a,b) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_LF("I,F") * u1_aa("e,m");

                    // RDM_blks_["D1_bb_ov"] += 0.250000 r2_abab(F,b,a,i,j) u1_bb(e,m) m2_abab(I,i,j,b,a) 
                    // +               0.250000 r2_abab(F,b,a,j,i) u1_bb(e,m) m2_abab(I,j,i,b,a) 
                    // +               0.250000 r2_abab(F,a,b,i,j) u1_bb(e,m) m2_abab(I,i,j,a,b) 
                    // +               0.250000 r2_abab(F,a,b,j,i) u1_bb(e,m) m2_abab(I,j,i,a,b) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_LF("I,F") * u1_bb("e,m");
                }

                {

                    // tempOp_LF += 1.000000 l1_aa("I,i,a") r1_aa("F,a,i") 
                    // flops: o1v1L1F1: 1, o0v0L1F1: 1 | mem: o0v0L1F1: 2, 
                    tempOp_LF("I,F") = l1_aa("I,i,a") * r1_aa("F,a,i");

                    // RDM_blks_["D1_aa_oo"] += 1.000000 d_aa(m,n) l1_aa(I,i,a) r1_aa(F,a,i) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["aa_oo"]("m,n");

                    // RDM_blks_["D1_bb_oo"] += 1.000000 d_bb(m,n) l1_aa(I,i,a) r1_aa(F,a,i) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["bb_oo"]("m,n");

                    // RDM_blks_["D1_aa_ov"] += 1.000000 l1_aa(I,i,a) r1_aa(F,a,i) t1_aa(e,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_aa("e,m");

                    // RDM_blks_["D1_bb_ov"] += 1.000000 l1_aa(I,i,a) r1_aa(F,a,i) t1_bb(e,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_bb("e,m");
                }

                if (include_m1_ && include_s1_) {

                    // tempOp_LF += 1.000000 m1_bb("I,i,a") s1_bb("F,a,i") 
                    // flops: o1v1L1F1: 1, o0v0L1F1: 1 | mem: o0v0L1F1: 2, 
                    tempOp_LF("I,F") = m1_bb("I,i,a") * s1_bb("F,a,i");

                    // RDM_blks_["D1_aa_oo"] += 1.000000 d_aa(m,n) m1_bb(I,i,a) s1_bb(F,a,i) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["aa_oo"]("m,n");

                    // RDM_blks_["D1_bb_oo"] += 1.000000 d_bb(m,n) m1_bb(I,i,a) s1_bb(F,a,i) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["bb_oo"]("m,n");

                    // RDM_blks_["D1_aa_ov"] += 1.000000 t1_aa(e,m) m1_bb(I,i,a) s1_bb(F,a,i) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_aa("e,m");

                    // RDM_blks_["D1_bb_ov"] += 1.000000 t1_bb(e,m) m1_bb(I,i,a) s1_bb(F,a,i) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_bb("e,m");

                    // tempOp_LF += 1.000000 m1_aa("I,i,a") s1_aa("F,a,i") 
                    // flops: o1v1L1F1: 1, o0v0L1F1: 1 | mem: o0v0L1F1: 2, 
                    tempOp_LF("I,F") = m1_aa("I,i,a") * s1_aa("F,a,i");

                    // RDM_blks_["D1_aa_oo"] += 1.000000 d_aa(m,n) m1_aa(I,i,a) s1_aa(F,a,i) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["aa_oo"]("m,n");

                    // RDM_blks_["D1_bb_oo"] += 1.000000 d_bb(m,n) m1_aa(I,i,a) s1_aa(F,a,i) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["bb_oo"]("m,n");

                    // RDM_blks_["D1_aa_ov"] += 1.000000 t1_aa(e,m) m1_aa(I,i,a) s1_aa(F,a,i) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_aa("e,m");

                    // RDM_blks_["D1_bb_ov"] += 1.000000 t1_bb(e,m) m1_aa(I,i,a) s1_aa(F,a,i) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_bb("e,m");
                }

                {

                    // tempOp_LF += 1.000000 r1_bb("F,a,i") l1_bb("I,i,a") 
                    // flops: o1v1L1F1: 1, o0v0L1F1: 1 | mem: o0v0L1F1: 2, 
                    tempOp_LF("I,F") = r1_bb("F,a,i") * l1_bb("I,i,a");

                    // RDM_blks_["D1_aa_oo"] += 1.000000 d_aa(m,n) l1_bb(I,i,a) r1_bb(F,a,i) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["aa_oo"]("m,n");

                    // RDM_blks_["D1_bb_oo"] += 1.000000 d_bb(m,n) l1_bb(I,i,a) r1_bb(F,a,i) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempOp_LF("I,F") * Id_blks["bb_oo"]("m,n");

                    // RDM_blks_["D1_aa_ov"] += 1.000000 l1_bb(I,i,a) r1_bb(F,a,i) t1_aa(e,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_aa("e,m");

                    // RDM_blks_["D1_bb_ov"] += 1.000000 l1_bb(I,i,a) r1_bb(F,a,i) t1_bb(e,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_LF("I,F") * t1_bb("e,m");
                }

                if (include_m1_) {

                    // tempOp_LF += 1.000000 m1_aa("I,i,a") r1_aa("F,a,i") 
                    // flops: o1v1L1F1: 1, o0v0L1F1: 1 | mem: o0v0L1F1: 2, 
                    tempOp_LF("I,F") = m1_aa("I,i,a") * r1_aa("F,a,i");
                }

                if (include_m1_ && include_u1_) {

                    // RDM_blks_["D1_aa_ov"] += 1.000000 r1_aa(F,a,i) u1_aa(e,m) m1_aa(I,i,a) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_LF("I,F") * u1_aa("e,m");

                    // RDM_blks_["D1_bb_ov"] += 1.000000 r1_aa(F,a,i) u1_bb(e,m) m1_aa(I,i,a) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_LF("I,F") * u1_bb("e,m");
                }

                if (include_m1_) {

                    // tempOp_LF += 1.000000 r1_bb("F,a,i") m1_bb("I,i,a") 
                    // flops: o1v1L1F1: 1, o0v0L1F1: 1 | mem: o0v0L1F1: 2, 
                    tempOp_LF("I,F") = r1_bb("F,a,i") * m1_bb("I,i,a");
                }

                if (include_m1_ && include_u1_) {

                    // RDM_blks_["D1_aa_ov"] += 1.000000 r1_bb(F,a,i) u1_aa(e,m) m1_bb(I,i,a) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempOp_LF("I,F") * u1_aa("e,m");

                    // RDM_blks_["D1_bb_ov"] += 1.000000 r1_bb(F,a,i) u1_bb(e,m) m1_bb(I,i,a) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempOp_LF("I,F") * u1_bb("e,m");
                }

                {

                    // RDM_blks_["D1_aa_oo"] += -1.000000 l1_aa(I,n,a) r1_aa(F,a,m) 
                    // flops: o2v2L1F1: 1, o2v1L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= l1_aa("I,n,a") * r1_aa("F,a,m");
                }

                if (include_m2_ && include_u2_) {

                    // RDM_blks_["D1_aa_oo"] += -0.500000 r0(F) u2_aaaa(b,a,i,m) m2_aaaa(I,i,n,b,a) 
                    // flops: o2v2L1F1: 1, o3v2L1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, o2v0L1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= 0.500000 * u2_aaaa("b,a,i,m") * m2_aaaa("I,i,n,b,a") * r0("F");
                }

                if (include_m1_ && include_s1_) {

                    // RDM_blks_["D1_aa_oo"] += -1.000000 m1_aa(I,n,a) s1_aa(F,a,m) 
                    // flops: o2v2L1F1: 1, o2v1L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= m1_aa("I,n,a") * s1_aa("F,a,m");
                }

                {

                    // RDM_blks_["D1_aa_oo"] += 1.000000 d_aa(m,n) l0(I) r0(F) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1, o0v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, o0v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += l0("I") * r0("F") * Id_blks["aa_oo"]("m,n");
                }

                if (include_m1_ && include_u1_) {

                    // RDM_blks_["D1_aa_oo"] += -1.000000 r0(F) u1_aa(a,m) m1_aa(I,n,a) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1, o2v1L1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, o2v0L1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= u1_aa("a,m") * m1_aa("I,n,a") * r0("F");
                }

                {

                    // RDM_blks_["D1_aa_oo"] += -1.000000 l1_aa(I,n,a) r0(F) t1_aa(a,m) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1, o2v1L1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, o2v0L1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= l1_aa("I,n,a") * t1_aa("a,m") * r0("F");
                }

                if (include_m0_ && include_s0_) {

                    // RDM_blks_["D1_aa_oo"] += 1.000000 d_aa(m,n) m0(I) s0(F) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1, o0v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, o0v0L1F1: 1, 
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += m0("I") * s0("F") * Id_blks["aa_oo"]("m,n");
                }

                {

                    // RDM_blks_["D1_bb_oo"] += -1.000000 l1_bb(I,n,a) r0(F) t1_bb(a,m) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1, o2v1L1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, o2v0L1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= l1_bb("I,n,a") * t1_bb("a,m") * r0("F");

                    // RDM_blks_["D1_bb_oo"] += 1.000000 d_bb(m,n) l0(I) r0(F) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1, o0v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, o0v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += l0("I") * r0("F") * Id_blks["bb_oo"]("m,n");
                }

                if (include_m1_ && include_u1_) {

                    // RDM_blks_["D1_bb_oo"] += -1.000000 r0(F) u1_bb(a,m) m1_bb(I,n,a) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1, o2v1L1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, o2v0L1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= u1_bb("a,m") * m1_bb("I,n,a") * r0("F");
                }

                if (include_m1_ && include_s1_) {

                    // RDM_blks_["D1_bb_oo"] += -1.000000 m1_bb(I,n,a) s1_bb(F,a,m) 
                    // flops: o2v2L1F1: 1, o2v1L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= m1_bb("I,n,a") * s1_bb("F,a,m");
                }

                if (include_m2_ && include_u2_) {

                    // RDM_blks_["D1_bb_oo"] += -0.500000 r0(F) u2_bbbb(b,a,i,m) m2_bbbb(I,i,n,b,a) 
                    // flops: o2v2L1F1: 1, o3v2L1: 1, o2v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, o2v0L1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= 0.500000 * u2_bbbb("b,a,i,m") * m2_bbbb("I,i,n,b,a") * r0("F");
                }

                {

                    // RDM_blks_["D1_bb_oo"] += -1.000000 l1_bb(I,n,a) r1_bb(F,a,m) 
                    // flops: o2v2L1F1: 1, o2v1L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= l1_bb("I,n,a") * r1_bb("F,a,m");
                }

                if (include_m0_ && include_s0_) {

                    // RDM_blks_["D1_bb_oo"] += 1.000000 d_bb(m,n) m0(I) s0(F) 
                    // flops: o2v2L1F1: 1, o2v0L1F1: 1, o0v0L1F1: 1 | mem: o2v2L1F1: 1, o2v0L1F1: 1, o0v0L1F1: 1, 
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += m0("I") * s0("F") * Id_blks["bb_oo"]("m,n");
                }

                {

                    // RDM_blks_["D1_aa_vv"] += 0.500000 l2_abab(I,j,i,e,a) r2_abab(F,f,a,j,i) 
                    // +               0.500000 l2_abab(I,i,j,e,a) r2_abab(F,f,a,i,j) 
                    // flops: o2v3L1F1: 1, o0v4L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += l2_abab("I,j,i,e,a") * r2_abab("F,f,a,j,i");
                }

                if (include_m2_ && include_s2_) {

                    // RDM_blks_["D1_aa_vv"] += 0.500000 m2_abab(I,j,i,e,a) s2_abab(F,f,a,j,i) 
                    // +               0.500000 m2_abab(I,i,j,e,a) s2_abab(F,f,a,i,j) 
                    // flops: o2v3L1F1: 1, o0v4L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += m2_abab("I,j,i,e,a") * s2_abab("F,f,a,j,i");
                }

                if (include_m2_ && include_u2_) {

                    // RDM_blks_["D1_aa_vv"] += 0.500000 r0(F) u2_aaaa(f,a,i,j) m2_aaaa(I,i,j,e,a) 
                    // flops: o0v4L1F1: 1, o2v3L1: 1, o0v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, o0v2L1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += 0.500000 * u2_aaaa("f,a,i,j") * m2_aaaa("I,i,j,e,a") * r0("F");
                }

                {

                    // RDM_blks_["D1_aa_vv"] += 1.000000 l1_aa(I,i,e) r0(F) t1_aa(f,i) 
                    // flops: o0v4L1F1: 1, o0v2L1F1: 1, o1v2L1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, o0v2L1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += l1_aa("I,i,e") * t1_aa("f,i") * r0("F");
                }

                if (include_m2_ && include_s0_) {

                    // RDM_blks_["D1_aa_vv"] += 0.500000 t2_aaaa(f,a,i,j) m2_aaaa(I,i,j,e,a) s0(F) 
                    // flops: o0v4L1F1: 1, o2v3L1: 1, o0v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, o0v2L1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += 0.500000 * t2_aaaa("f,a,i,j") * m2_aaaa("I,i,j,e,a") * s0("F");
                }

                if (include_m2_ && include_s2_) {

                    // RDM_blks_["D1_aa_vv"] += 0.500000 m2_aaaa(I,j,i,e,a) s2_aaaa(F,f,a,j,i) 
                    // flops: o2v3L1F1: 1, o0v4L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += 0.500000 * m2_aaaa("I,j,i,e,a") * s2_aaaa("F,f,a,j,i");
                }

                {

                    // RDM_blks_["D1_aa_vv"] += 0.500000 l2_aaaa(I,j,i,e,a) r2_aaaa(F,f,a,j,i) 
                    // flops: o2v3L1F1: 1, o0v4L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += 0.500000 * l2_aaaa("I,j,i,e,a") * r2_aaaa("F,f,a,j,i");

                    // RDM_blks_["D1_aa_vv"] += 0.500000 l2_aaaa(I,i,j,e,a) r0(F) t2_aaaa(f,a,i,j) 
                    // flops: o0v4L1F1: 1, o2v3L1: 1, o0v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, o0v2L1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += 0.500000 * l2_aaaa("I,i,j,e,a") * t2_aaaa("f,a,i,j") * r0("F");
                }

                if (include_m1_ && include_s1_) {

                    // RDM_blks_["D1_aa_vv"] += 1.000000 m1_aa(I,i,e) s1_aa(F,f,i) 
                    // flops: o0v4L1F1: 1, o1v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += m1_aa("I,i,e") * s1_aa("F,f,i");
                }

                {

                    // RDM_blks_["D1_aa_vv"] += 1.000000 l1_aa(I,i,e) r1_aa(F,f,i) 
                    // flops: o0v4L1F1: 1, o1v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += l1_aa("I,i,e") * r1_aa("F,f,i");
                }

                if (include_m1_ && include_u1_) {

                    // RDM_blks_["D1_aa_vv"] += 1.000000 r0(F) u1_aa(f,i) m1_aa(I,i,e) 
                    // flops: o0v4L1F1: 1, o0v2L1F1: 1, o1v2L1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, o0v2L1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += u1_aa("f,i") * m1_aa("I,i,e") * r0("F");
                }

                if (include_m1_ && include_s0_) {

                    // RDM_blks_["D1_aa_vv"] += 1.000000 t1_aa(f,i) m1_aa(I,i,e) s0(F) 
                    // flops: o0v4L1F1: 1, o0v2L1F1: 1, o1v2L1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, o0v2L1: 1, 
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += t1_aa("f,i") * m1_aa("I,i,e") * s0("F");
                }

                if (include_m2_ && include_u2_) {

                    // RDM_blks_["D1_bb_vv"] += 0.500000 r0(F) u2_bbbb(f,a,i,j) m2_bbbb(I,i,j,e,a) 
                    // flops: o0v4L1F1: 1, o2v3L1: 1, o0v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, o0v2L1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += 0.500000 * u2_bbbb("f,a,i,j") * m2_bbbb("I,i,j,e,a") * r0("F");
                }

                if (include_m1_ && include_s0_) {

                    // RDM_blks_["D1_bb_vv"] += 1.000000 t1_bb(f,i) m1_bb(I,i,e) s0(F) 
                    // flops: o0v4L1F1: 1, o0v2L1F1: 1, o1v2L1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, o0v2L1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += t1_bb("f,i") * m1_bb("I,i,e") * s0("F");
                }

                {

                    // RDM_blks_["D1_bb_vv"] += 0.500000 l2_bbbb(I,i,j,e,a) r0(F) t2_bbbb(f,a,i,j) 
                    // flops: o0v4L1F1: 1, o2v3L1: 1, o0v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, o0v2L1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += 0.500000 * l2_bbbb("I,i,j,e,a") * t2_bbbb("f,a,i,j") * r0("F");
                }

                if (include_m2_ && include_s0_) {

                    // RDM_blks_["D1_bb_vv"] += 0.500000 t2_bbbb(f,a,i,j) m2_bbbb(I,i,j,e,a) s0(F) 
                    // flops: o0v4L1F1: 1, o2v3L1: 1, o0v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, o0v2L1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += 0.500000 * t2_bbbb("f,a,i,j") * m2_bbbb("I,i,j,e,a") * s0("F");
                }

                if (include_m2_ && include_s2_) {

                    // RDM_blks_["D1_bb_vv"] += 0.500000 m2_abab(I,j,i,a,e) s2_abab(F,a,f,j,i) 
                    // +               0.500000 m2_abab(I,i,j,a,e) s2_abab(F,a,f,i,j) 
                    // flops: o2v3L1F1: 1, o0v4L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += m2_abab("I,j,i,a,e") * s2_abab("F,a,f,j,i");
                }

                {

                    // RDM_blks_["D1_bb_vv"] += 0.500000 l2_abab(I,j,i,a,e) r2_abab(F,a,f,j,i) 
                    // +               0.500000 l2_abab(I,i,j,a,e) r2_abab(F,a,f,i,j) 
                    // flops: o2v3L1F1: 1, o0v4L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += l2_abab("I,j,i,a,e") * r2_abab("F,a,f,j,i");
                }

                if (include_m2_ && include_s2_) {

                    // RDM_blks_["D1_bb_vv"] += 0.500000 m2_bbbb(I,j,i,e,a) s2_bbbb(F,f,a,j,i) 
                    // flops: o2v3L1F1: 1, o0v4L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += 0.500000 * m2_bbbb("I,j,i,e,a") * s2_bbbb("F,f,a,j,i");
                }

                {

                    // RDM_blks_["D1_bb_vv"] += 1.000000 l1_bb(I,i,e) r0(F) t1_bb(f,i) 
                    // flops: o0v4L1F1: 1, o0v2L1F1: 1, o1v2L1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, o0v2L1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += l1_bb("I,i,e") * t1_bb("f,i") * r0("F");

                    // RDM_blks_["D1_bb_vv"] += 1.000000 l1_bb(I,i,e) r1_bb(F,f,i) 
                    // flops: o0v4L1F1: 1, o1v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += l1_bb("I,i,e") * r1_bb("F,f,i");
                }

                if (include_m1_ && include_s1_) {

                    // RDM_blks_["D1_bb_vv"] += 1.000000 m1_bb(I,i,e) s1_bb(F,f,i) 
                    // flops: o0v4L1F1: 1, o1v2L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += m1_bb("I,i,e") * s1_bb("F,f,i");
                }

                if (include_m1_ && include_u1_) {

                    // RDM_blks_["D1_bb_vv"] += 1.000000 r0(F) u1_bb(f,i) m1_bb(I,i,e) 
                    // flops: o0v4L1F1: 1, o0v2L1F1: 1, o1v2L1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, o0v2L1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += u1_bb("f,i") * m1_bb("I,i,e") * r0("F");
                }

                {

                    // RDM_blks_["D1_bb_vv"] += 0.500000 l2_bbbb(I,j,i,e,a) r2_bbbb(F,f,a,j,i) 
                    // flops: o2v3L1F1: 1, o0v4L1F1: 1 | mem: o0v4L1F1: 1, o0v2L1F1: 1, 
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += 0.500000 * l2_bbbb("I,j,i,e,a") * r2_bbbb("F,f,a,j,i");

                    // RDM_blks_["D1_aa_ov"] += 1.000000 l1_bb(I,i,a) r0(F) t2_abab(e,a,m,i) 
                    // flops: o1v3L1F1: 1, o2v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += l1_bb("I,i,a") * t2_abab("e,a,m,i") * r0("F");
                }

                if (include_m1_ && include_u1_) {

                    // RDM_blks_["D1_aa_ov"] += -1.000000 r1_aa(F,e,i) u1_aa(a,m) m1_aa(I,i,a) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= u1_aa("a,m") * m1_aa("I,i,a") * r1_aa("F,e,i");
                }

                {

                    // RDM_blks_["D1_aa_ov"] += 1.000000 l0(I) r0(F) t1_aa(e,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o0v0L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o0v0L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += l0("I") * r0("F") * t1_aa("e,m");

                    // RDM_blks_["D1_aa_ov"] += -1.000000 l1_aa(I,i,a) r0(F) t1_aa(e,i) t1_aa(a,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 2 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, o2v0L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= l1_aa("I,i,a") * t1_aa("a,m") * t1_aa("e,i") * r0("F");
                }

                if (include_m1_ && include_u2_) {

                    // RDM_blks_["D1_aa_ov"] += -1.000000 r0(F) u2_aaaa(e,a,i,m) m1_aa(I,i,a) 
                    // flops: o1v3L1F1: 1, o2v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= u2_aaaa("e,a,i,m") * m1_aa("I,i,a") * r0("F");
                }

                if (include_m2_ && include_u2_) {

                    // RDM_blks_["D1_aa_ov"] += 0.500000 r0(F) t1_aa(a,m) u2_aaaa(e,b,i,j) m2_aaaa(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o3v2L1: 2, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o3v1L1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += 0.500000 * t1_aa("a,m") * m2_aaaa("I,i,j,b,a") * u2_aaaa("e,b,i,j") * r0("F");
                }

                if (include_m1_ && include_u1_) {

                    // RDM_blks_["D1_aa_ov"] += -1.000000 r0(F) t1_aa(e,i) u1_aa(a,m) m1_aa(I,i,a) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 2 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, o2v0L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= u1_aa("a,m") * m1_aa("I,i,a") * t1_aa("e,i") * r0("F");
                }

                if (include_m0_ && include_u1_) {

                    // RDM_blks_["D1_aa_ov"] += 1.000000 r0(F) u1_aa(e,m) m0(I) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o0v0L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o0v0L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += r0("F") * m0("I") * u1_aa("e,m");
                }

                if (include_m1_ && include_s1_) {

                    // RDM_blks_["D1_aa_ov"] += -1.000000 t1_aa(e,i) m1_aa(I,i,a) s1_aa(F,a,m) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 2 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= m1_aa("I,i,a") * s1_aa("F,a,m") * t1_aa("e,i");
                }

                if (include_m2_ && include_s1_) {

                    // RDM_blks_["D1_aa_ov"] += 0.500000 t2_aaaa(b,a,i,m) m2_aaaa(I,j,i,b,a) s1_aa(F,e,j) 
                    // flops: o1v3L1F1: 1, o3v2L1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += 0.500000 * t2_aaaa("b,a,i,m") * m2_aaaa("I,j,i,b,a") * s1_aa("F,e,j");
                }

                {

                    // RDM_blks_["D1_aa_ov"] += 0.500000 l2_aaaa(I,j,i,b,a) r1_aa(F,e,j) t2_aaaa(b,a,i,m) 
                    // flops: o1v3L1F1: 1, o3v2L1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += 0.500000 * l2_aaaa("I,j,i,b,a") * t2_aaaa("b,a,i,m") * r1_aa("F,e,j");
                }

                if (include_m1_ && include_s0_) {

                    // RDM_blks_["D1_aa_ov"] += -1.000000 t2_aaaa(e,a,i,m) m1_aa(I,i,a) s0(F) 
                    // flops: o1v3L1F1: 1, o2v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t2_aaaa("e,a,i,m") * m1_aa("I,i,a") * s0("F");
                }

                if (include_m1_ && include_u1_) {

                    // RDM_blks_["D1_aa_ov"] += -1.000000 r1_aa(F,a,m) u1_aa(e,i) m1_aa(I,i,a) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 2 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= r1_aa("F,a,m") * m1_aa("I,i,a") * u1_aa("e,i");
                }

                {

                    // RDM_blks_["D1_aa_ov"] += 1.000000 l0(I) r1_aa(F,e,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += l0("I") * r1_aa("F,e,m");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_aa_ov"] += -0.500000 r2_aaaa(F,b,a,i,m) u1_aa(e,j) m2_aaaa(I,i,j,b,a) 
                    // flops: o3v2L1F1: 1, o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= 0.500000 * r2_aaaa("F,b,a,i,m") * m2_aaaa("I,i,j,b,a") * u1_aa("e,j");
                }

                {

                    // RDM_blks_["D1_aa_ov"] += 0.500000 l2_aaaa(I,i,j,a,b) r1_aa(F,b,m) t2_aaaa(e,a,i,j) 
                    // flops: o1v3L1F1: 1, o2v3L1: 1, o1v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o0v2L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += 0.500000 * l2_aaaa("I,i,j,a,b") * t2_aaaa("e,a,i,j") * r1_aa("F,b,m");

                    // RDM_blks_["D1_aa_ov"] += -1.000000 l1_aa(I,i,a) r2_aaaa(F,e,a,i,m) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= l1_aa("I,i,a") * r2_aaaa("F,e,a,i,m");
                }

                if (include_m2_ && include_u2_) {

                    // RDM_blks_["D1_aa_ov"] += 0.500000 r1_aa(F,a,m) u2_aaaa(e,b,i,j) m2_aaaa(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o2v3L1: 1, o1v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o0v2L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += 0.500000 * u2_aaaa("e,b,i,j") * m2_aaaa("I,i,j,b,a") * r1_aa("F,a,m");
                }

                if (include_m1_ && include_s0_) {

                    // RDM_blks_["D1_aa_ov"] += 1.000000 t2_abab(e,a,m,i) m1_bb(I,i,a) s0(F) 
                    // flops: o1v3L1F1: 1, o2v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += t2_abab("e,a,m,i") * m1_bb("I,i,a") * s0("F");
                }

                {

                    // RDM_blks_["D1_aa_ov"] += -1.000000 l1_aa(I,i,a) r1_aa(F,a,m) t1_aa(e,i) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 2 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= l1_aa("I,i,a") * r1_aa("F,a,m") * t1_aa("e,i");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_aa_ov"] += -0.500000 r2_abab(F,b,a,m,i) u1_aa(e,j) m2_abab(I,j,i,b,a) 
                    // +               -0.500000 r2_abab(F,a,b,m,i) u1_aa(e,j) m2_abab(I,j,i,a,b) 
                    // flops: o3v2L1F1: 1, o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= r2_abab("F,b,a,m,i") * m2_abab("I,j,i,b,a") * u1_aa("e,j");
                }

                if (include_m0_ && include_s0_) {

                    // RDM_blks_["D1_aa_ov"] += 1.000000 t1_aa(e,m) m0(I) s0(F) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o0v0L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o0v0L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += m0("I") * s0("F") * t1_aa("e,m");
                }

                {

                    // RDM_blks_["D1_aa_ov"] += -1.000000 l1_aa(I,i,a) r1_aa(F,e,i) t1_aa(a,m) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= l1_aa("I,i,a") * t1_aa("a,m") * r1_aa("F,e,i");
                }

                if (include_m0_ && include_s1_) {

                    // RDM_blks_["D1_aa_ov"] += 1.000000 m0(I) s1_aa(F,e,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += m0("I") * s1_aa("F,e,m");
                }

                {

                    // RDM_blks_["D1_aa_ov"] += -1.000000 l1_aa(I,i,a) r0(F) t2_aaaa(e,a,i,m) 
                    // flops: o1v3L1F1: 1, o2v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= l1_aa("I,i,a") * t2_aaaa("e,a,i,m") * r0("F");
                }

                if (include_m1_ && include_s2_) {

                    // RDM_blks_["D1_aa_ov"] += 1.000000 m1_bb(I,i,a) s2_abab(F,e,a,m,i) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += m1_bb("I,i,a") * s2_abab("F,e,a,m,i");
                }

                if (include_m1_ && include_u2_) {

                    // RDM_blks_["D1_aa_ov"] += 1.000000 r0(F) u2_abab(e,a,m,i) m1_bb(I,i,a) 
                    // flops: o1v3L1F1: 1, o2v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += u2_abab("e,a,m,i") * m1_bb("I,i,a") * r0("F");
                }

                {

                    // RDM_blks_["D1_aa_ov"] += 1.000000 l1_bb(I,i,a) r2_abab(F,e,a,m,i) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += l1_bb("I,i,a") * r2_abab("F,e,a,m,i");
                }

                if (include_m1_ && include_s2_) {

                    // RDM_blks_["D1_aa_ov"] += -1.000000 m1_aa(I,i,a) s2_aaaa(F,e,a,i,m) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= m1_aa("I,i,a") * s2_aaaa("F,e,a,i,m");
                }

                if (include_m2_ && include_s1_) {

                    // RDM_blks_["D1_aa_ov"] += 0.500000 t2_aaaa(e,a,i,j) m2_aaaa(I,i,j,a,b) s1_aa(F,b,m) 
                    // flops: o1v3L1F1: 1, o2v3L1: 1, o1v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o0v2L1: 1, 
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += 0.500000 * t2_aaaa("e,a,i,j") * m2_aaaa("I,i,j,a,b") * s1_aa("F,b,m");
                }

                {

                    // RDM_blks_["D1_bb_ov"] += -1.000000 l1_bb(I,i,a) r1_bb(F,e,i) t1_bb(a,m) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= l1_bb("I,i,a") * t1_bb("a,m") * r1_bb("F,e,i");
                }

                if (include_m1_ && include_s0_) {

                    // RDM_blks_["D1_bb_ov"] += 1.000000 t2_abab(a,e,i,m) m1_aa(I,i,a) s0(F) 
                    // flops: o1v3L1F1: 1, o2v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += t2_abab("a,e,i,m") * m1_aa("I,i,a") * s0("F");

                    // RDM_blks_["D1_bb_ov"] += -1.000000 t2_bbbb(e,a,i,m) m1_bb(I,i,a) s0(F) 
                    // flops: o1v3L1F1: 1, o2v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t2_bbbb("e,a,i,m") * m1_bb("I,i,a") * s0("F");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_bb_ov"] += -0.500000 r2_abab(F,b,a,i,m) u1_bb(e,j) m2_abab(I,i,j,b,a) 
                    // +               -0.500000 r2_abab(F,a,b,i,m) u1_bb(e,j) m2_abab(I,i,j,a,b) 
                    // flops: o3v2L1F1: 1, o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= r2_abab("F,b,a,i,m") * m2_abab("I,i,j,b,a") * u1_bb("e,j");
                }

                if (include_m1_ && include_u2_) {

                    // RDM_blks_["D1_bb_ov"] += 1.000000 r0(F) u2_abab(a,e,i,m) m1_aa(I,i,a) 
                    // flops: o1v3L1F1: 1, o2v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += u2_abab("a,e,i,m") * m1_aa("I,i,a") * r0("F");
                }

                {

                    // RDM_blks_["D1_bb_ov"] += 1.000000 l1_aa(I,i,a) r0(F) t2_abab(a,e,i,m) 
                    // flops: o1v3L1F1: 1, o2v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += l1_aa("I,i,a") * t2_abab("a,e,i,m") * r0("F");
                }

                if (include_m2_ && include_u1_) {

                    // RDM_blks_["D1_bb_ov"] += -0.500000 r2_bbbb(F,b,a,i,m) u1_bb(e,j) m2_bbbb(I,i,j,b,a) 
                    // flops: o3v2L1F1: 1, o1v3L1F1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= 0.500000 * r2_bbbb("F,b,a,i,m") * m2_bbbb("I,i,j,b,a") * u1_bb("e,j");
                }

                if (include_m1_ && include_s1_) {

                    // RDM_blks_["D1_bb_ov"] += -1.000000 t1_bb(e,i) m1_bb(I,i,a) s1_bb(F,a,m) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 2 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= m1_bb("I,i,a") * s1_bb("F,a,m") * t1_bb("e,i");
                }

                if (include_m2_ && include_u2_) {

                    // RDM_blks_["D1_bb_ov"] += 0.500000 r0(F) t1_bb(a,m) u2_bbbb(e,b,i,j) m2_bbbb(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o3v2L1: 2, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o3v1L1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += 0.500000 * t1_bb("a,m") * m2_bbbb("I,i,j,b,a") * u2_bbbb("e,b,i,j") * r0("F");
                }

                {

                    // RDM_blks_["D1_bb_ov"] += -1.000000 l1_bb(I,i,a) r0(F) t2_bbbb(e,a,i,m) 
                    // flops: o1v3L1F1: 1, o2v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= l1_bb("I,i,a") * t2_bbbb("e,a,i,m") * r0("F");

                    // RDM_blks_["D1_bb_ov"] += 0.500000 l2_bbbb(I,j,i,b,a) r1_bb(F,e,j) t2_bbbb(b,a,i,m) 
                    // flops: o1v3L1F1: 1, o3v2L1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += 0.500000 * l2_bbbb("I,j,i,b,a") * t2_bbbb("b,a,i,m") * r1_bb("F,e,j");
                }

                if (include_m1_ && include_s2_) {

                    // RDM_blks_["D1_bb_ov"] += 1.000000 m1_aa(I,i,a) s2_abab(F,a,e,i,m) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += m1_aa("I,i,a") * s2_abab("F,a,e,i,m");
                }

                if (include_m1_ && include_u2_) {

                    // RDM_blks_["D1_bb_ov"] += -1.000000 r0(F) u2_bbbb(e,a,i,m) m1_bb(I,i,a) 
                    // flops: o1v3L1F1: 1, o2v2L1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= u2_bbbb("e,a,i,m") * m1_bb("I,i,a") * r0("F");
                }

                {

                    // RDM_blks_["D1_bb_ov"] += -1.000000 l1_bb(I,i,a) r1_bb(F,a,m) t1_bb(e,i) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 2 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= l1_bb("I,i,a") * r1_bb("F,a,m") * t1_bb("e,i");
                }

                if (include_m1_ && include_s2_) {

                    // RDM_blks_["D1_bb_ov"] += -1.000000 m1_bb(I,i,a) s2_bbbb(F,e,a,i,m) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= m1_bb("I,i,a") * s2_bbbb("F,e,a,i,m");
                }

                {

                    // RDM_blks_["D1_bb_ov"] += 0.500000 l2_bbbb(I,i,j,a,b) r1_bb(F,b,m) t2_bbbb(e,a,i,j) 
                    // flops: o1v3L1F1: 1, o2v3L1: 1, o1v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o0v2L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += 0.500000 * l2_bbbb("I,i,j,a,b") * t2_bbbb("e,a,i,j") * r1_bb("F,b,m");
                }

                if (include_m2_ && include_s1_) {

                    // RDM_blks_["D1_bb_ov"] += 0.500000 t2_bbbb(b,a,i,m) m2_bbbb(I,j,i,b,a) s1_bb(F,e,j) 
                    // flops: o1v3L1F1: 1, o3v2L1: 1, o2v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += 0.500000 * t2_bbbb("b,a,i,m") * m2_bbbb("I,j,i,b,a") * s1_bb("F,e,j");
                }

                {

                    // RDM_blks_["D1_bb_ov"] += 1.000000 l0(I) r0(F) t1_bb(e,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o0v0L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o0v0L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += l0("I") * r0("F") * t1_bb("e,m");
                }

                if (include_m1_ && include_u1_) {

                    // RDM_blks_["D1_bb_ov"] += -1.000000 r1_bb(F,a,m) u1_bb(e,i) m1_bb(I,i,a) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 2 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= r1_bb("F,a,m") * m1_bb("I,i,a") * u1_bb("e,i");

                    // RDM_blks_["D1_bb_ov"] += -1.000000 r1_bb(F,e,i) u1_bb(a,m) m1_bb(I,i,a) 
                    // flops: o1v3L1F1: 1, o2v1L1F1: 1, o2v1L1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o2v0L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= u1_bb("a,m") * m1_bb("I,i,a") * r1_bb("F,e,i");
                }

                {

                    // RDM_blks_["D1_bb_ov"] += 1.000000 l0(I) r1_bb(F,e,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += l0("I") * r1_bb("F,e,m");
                }

                if (include_m0_ && include_u1_) {

                    // RDM_blks_["D1_bb_ov"] += 1.000000 r0(F) u1_bb(e,m) m0(I) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o0v0L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o0v0L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += r0("F") * m0("I") * u1_bb("e,m");
                }

                if (include_m0_ && include_s0_) {

                    // RDM_blks_["D1_bb_ov"] += 1.000000 t1_bb(e,m) m0(I) s0(F) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o0v0L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o0v0L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += m0("I") * s0("F") * t1_bb("e,m");
                }

                {

                    // RDM_blks_["D1_bb_ov"] += -1.000000 l1_bb(I,i,a) r2_bbbb(F,e,a,i,m) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= l1_bb("I,i,a") * r2_bbbb("F,e,a,i,m");
                }

                if (include_m1_ && include_u1_) {

                    // RDM_blks_["D1_bb_ov"] += -1.000000 r0(F) t1_bb(e,i) u1_bb(a,m) m1_bb(I,i,a) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 2 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, o2v0L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= u1_bb("a,m") * m1_bb("I,i,a") * t1_bb("e,i") * r0("F");
                }

                if (include_m2_ && include_u2_) {

                    // RDM_blks_["D1_bb_ov"] += 0.500000 r1_bb(F,a,m) u2_bbbb(e,b,i,j) m2_bbbb(I,i,j,b,a) 
                    // flops: o1v3L1F1: 1, o2v3L1: 1, o1v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o0v2L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += 0.500000 * u2_bbbb("e,b,i,j") * m2_bbbb("I,i,j,b,a") * r1_bb("F,a,m");
                }

                if (include_m2_ && include_s1_) {

                    // RDM_blks_["D1_bb_ov"] += 0.500000 t2_bbbb(e,a,i,j) m2_bbbb(I,i,j,a,b) s1_bb(F,b,m) 
                    // flops: o1v3L1F1: 1, o2v3L1: 1, o1v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o0v2L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += 0.500000 * t2_bbbb("e,a,i,j") * m2_bbbb("I,i,j,a,b") * s1_bb("F,b,m");
                }

                {

                    // RDM_blks_["D1_bb_ov"] += 1.000000 l1_aa(I,i,a) r2_abab(F,a,e,i,m) 
                    // flops: o1v3L1F1: 1, o2v2L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += l1_aa("I,i,a") * r2_abab("F,a,e,i,m");
                }

                if (include_m0_ && include_s1_) {

                    // RDM_blks_["D1_bb_ov"] += 1.000000 m0(I) s1_bb(F,e,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += m0("I") * s1_bb("F,e,m");
                }

                {

                    // RDM_blks_["D1_bb_ov"] += -1.000000 l1_bb(I,i,a) r0(F) t1_bb(e,i) t1_bb(a,m) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1, o2v1L1: 2 | mem: o1v3L1F1: 1, o1v1L1F1: 1, o1v1L1: 1, o2v0L1: 1, 
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= l1_bb("I,i,a") * t1_bb("a,m") * t1_bb("e,i") * r0("F");
                }

                if (include_m1_ && include_s0_) {

                    // RDM_blks_["D1_aa_vo"] += 1.000000 m1_aa(I,m,e) s0(F) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_vo"]("I,F,e,m") += m1_aa("I,m,e") * s0("F");
                }

                {

                    // RDM_blks_["D1_aa_vo"] += 1.000000 l1_aa(I,m,e) r0(F) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_aa_vo"]("I,F,e,m") += l1_aa("I,m,e") * r0("F");
                }

                if (include_m1_ && include_s0_) {

                    // RDM_blks_["D1_bb_vo"] += 1.000000 m1_bb(I,m,e) s0(F) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_vo"]("I,F,e,m") += m1_bb("I,m,e") * s0("F");
                }

                {

                    // RDM_blks_["D1_bb_vo"] += 1.000000 l1_bb(I,m,e) r0(F) 
                    // flops: o1v3L1F1: 1, o1v1L1F1: 1 | mem: o1v3L1F1: 1, o1v1L1F1: 1, 
                    RDM_blks_["D1_bb_vo"]("I,F,e,m") += l1_bb("I,m,e") * r0("F");
                }

/*

                Total Number of Terms: 387
                Number of Flops: (old) 1246 -> (new) 884

                Total FLOP scaling:
                ------------------
                Scaling :      new |      old |     diff
                -------- : -------- | -------- | --------
                o2v3L1F1 :        8 |       12 |       -4
                o3v2L1F1 :       24 |       48 |      -24
                o0v4L1F1 :       42 |       52 |      -10
                o1v3L1F1 :      212 |      270 |      -58
                o2v2L1F1 :      119 |      224 |     -105
                o2v3L1   :       18 |       36 |      -18
                o3v1L1F1 :       28 |       28 |        0
                o3v2L1   :       48 |      154 |     -106
                o1v2L1F1 :       28 |       34 |       -6
                o2v1L1F1 :       62 |       74 |      -12
                o2v2L1   :       12 |       12 |        0
                o3v1L1   :       12 |        0 |       12
                o0v2L1F1 :       18 |       24 |       -6
                o1v1L1F1 :      114 |      150 |      -36
                o1v2L1   :        6 |        6 |        0
                o2v0L1F1 :       50 |       60 |      -10
                o2v1L1   :       38 |       52 |      -14
                o0v2L1   :        6 |        0 |        6
                o2v0L1   :       14 |        0 |       14
                o0v0L1F1 :       25 |       10 |       15

                Total MEM scaling:
                o0v4L1F1 :       42 |       52 |      -10
                o1v3L1F1 :      212 |      270 |      -58
                o2v2L1F1 :       66 |       88 |      -22
                o3v1L1   :       26 |       70 |      -44
                o0v2L1F1 :       42 |       52 |      -10
                o1v1L1F1 :      228 |      306 |      -78
                o2v0L1F1 :      100 |      128 |      -28
                o0v2L1   :       30 |       42 |      -12
                o1v1L1   :       52 |       68 |      -16
                o2v0L1   :       46 |       80 |      -34
                o0v0L1F1 :       40 |       90 |      -50
*/

            }
        }

        rdm_timer.stop();
        if (world_.rank() == 0) {
//            outfile->Printf("\n  Building EOM-EE-CC RDMs... ");
            outfile->Printf("Done --> Time: %s\n", rdm_timer.elapsed().c_str());
        }
    }

} // cc_cavity