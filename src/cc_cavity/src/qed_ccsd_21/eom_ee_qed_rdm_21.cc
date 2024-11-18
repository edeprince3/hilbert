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

#include "cc_cavity/include/qed_ccsd_21/eom_ee_qed_rdm_21.h"
#include <psi4/libmints/writer.h>
#include "misc/nonsym_davidson_solver.h"

namespace hilbert {
    EOM_EE_QED_RDM_21::EOM_EE_QED_RDM_21(const shared_ptr<EOM_Driver>& eom_driver, Options & options) : EOM_EE_RDM(eom_driver, options) {
    }

    void EOM_EE_QED_RDM_21::compute_eom_1rdm() {

        // remove t1 transforms if they exist (shouldn't)
        if (eom_driver_->cc_wfn_->has_t1_integrals_)
            eom_driver_->cc_wfn_->transform_integrals(false);

        // reinitialize the 1-RDMs if not already done

        RDM_blks_["D1_aa_vv"] = makeTensor(world_, {M_, M_, va_, va_}, true);
        RDM_blks_["D1_bb_vv"] = makeTensor(world_, {M_, M_, vb_, vb_}, true);
        RDM_blks_["D1_aa_vo"] = makeTensor(world_, {M_, M_, va_, oa_}, true);
        RDM_blks_["D1_bb_vo"] = makeTensor(world_, {M_, M_, vb_, ob_}, true);
        RDM_blks_["D1_aa_ov"] = makeTensor(world_, {M_, M_, oa_, va_}, true);
        RDM_blks_["D1_bb_ov"] = makeTensor(world_, {M_, M_, ob_, vb_}, true);
        RDM_blks_["D1_aa_oo"] = makeTensor(world_, {M_, M_, oa_, oa_}, true);
        RDM_blks_["D1_bb_oo"] = makeTensor(world_, {M_, M_, ob_, ob_}, true);

        TArrayD &D_vv_aa = RDM_blks_["D1_aa_vv"];
        TArrayD &D_vv_bb = RDM_blks_["D1_bb_vv"];
        TArrayD &D_vo_aa = RDM_blks_["D1_aa_vo"];
        TArrayD &D_vo_bb = RDM_blks_["D1_bb_vo"];
        TArrayD &D_ov_aa = RDM_blks_["D1_aa_ov"];
        TArrayD &D_ov_bb = RDM_blks_["D1_bb_ov"];
        TArrayD &D_oo_aa = RDM_blks_["D1_aa_oo"];
        TArrayD &D_oo_bb = RDM_blks_["D1_bb_oo"];

        Timer rdm_timer; rdm_timer.start();

        // get cavity information
        double w0 = eom_driver_->cc_wfn_->cavity_frequency_[2];
        double coupling_factor_z = w0 * eom_driver_->cc_wfn_->cavity_coupling_strength_[2];

        // extract 0-body amplitudes
        TArrayD r0 = eom_driver_->evec_blks_["r0"];
        TArrayD l0 = eom_driver_->evec_blks_["l0"];

        double t0_1 = eom_driver_->cc_wfn_->scalars_["t0_1"];
        TArrayD r0_1 = eom_driver_->evec_blks_["l0_1"];
        TArrayD l0_1 = eom_driver_->evec_blks_["r0_1"];

        // extract 1-body amplitudes
        std::map<std::string, TA::TArrayD> t1 {
                {"aa", eom_driver_->cc_wfn_->amplitudes_["t1_aa"]},
                {"bb", eom_driver_->cc_wfn_->amplitudes_["t1_bb"]}
        };
        std::map<std::string, TA::TArrayD> r1 {
                {"aa", eom_driver_->evec_blks_["r1_aa"]},
                {"bb", eom_driver_->evec_blks_["r1_bb"]}
        };
        std::map<std::string, TA::TArrayD> l1 {
                {"aa", eom_driver_->evec_blks_["l1_aa"]},
                {"bb", eom_driver_->evec_blks_["l1_bb"]}
        };


        std::map<std::string, TA::TArrayD> t1_1 {
                {"aa", eom_driver_->cc_wfn_->amplitudes_["t1_1_aa"]},
                {"bb", eom_driver_->cc_wfn_->amplitudes_["t1_1_bb"]}
        };
        std::map<std::string, TA::TArrayD> r1_1 {
                {"aa", eom_driver_->evec_blks_["r1_1_aa"]},
                {"bb", eom_driver_->evec_blks_["r1_1_bb"]}
        };
        std::map<std::string, TA::TArrayD> l1_1 {
                {"aa", eom_driver_->evec_blks_["l1_1_aa"]},
                {"bb", eom_driver_->evec_blks_["l1_1_bb"]}
        };


        // extract 2-body amplitudes
        std::map<std::string, TA::TArrayD> t2 {
                {"aaaa", eom_driver_->cc_wfn_->amplitudes_["t2_aaaa"]},
                {"abab", eom_driver_->cc_wfn_->amplitudes_["t2_abab"]},
                {"bbbb", eom_driver_->cc_wfn_->amplitudes_["t2_bbbb"]}
        };
        std::map<std::string, TA::TArrayD> r2 {
                {"aaaa", eom_driver_->evec_blks_["r2_aaaa"]},
                {"abab", eom_driver_->evec_blks_["r2_abab"]},
                {"bbbb", eom_driver_->evec_blks_["r2_bbbb"]}
        };
        std::map<std::string, TA::TArrayD> l2 {
                {"aaaa", eom_driver_->evec_blks_["l2_aaaa"]},
                {"abab", eom_driver_->evec_blks_["l2_abab"]},
                {"bbbb", eom_driver_->evec_blks_["l2_bbbb"]}
        };

        std::map<std::string, TA::TArrayD> t2_1 {
                {"aaaa", eom_driver_->cc_wfn_->amplitudes_["t2_1_aaaa"]},
                {"abab", eom_driver_->cc_wfn_->amplitudes_["t2_1_abab"]},
                {"bbbb", eom_driver_->cc_wfn_->amplitudes_["t2_1_bbbb"]}
        };
        std::map<std::string, TA::TArrayD> r2_1 {
                {"aaaa", eom_driver_->evec_blks_["r2_1_aaaa"]},
                {"abab", eom_driver_->evec_blks_["r2_1_abab"]},
                {"bbbb", eom_driver_->evec_blks_["r2_1_bbbb"]}
        };
        std::map<std::string, TA::TArrayD> l2_1 {
                {"aaaa", eom_driver_->evec_blks_["l2_1_aaaa"]},
                {"abab", eom_driver_->evec_blks_["l2_1_abab"]},
                {"bbbb", eom_driver_->evec_blks_["l2_1_bbbb"]}
        };

        // get integrals
        TArrayMap &Id  = eom_driver_->cc_wfn_->Id_blks_;
        TArrayMap &eri = eom_driver_->cc_wfn_->V_blks_;
        TArrayMap &f   = eom_driver_->cc_wfn_->F_blks_;

        rdm_timer.start();
        {

            // D_ov_bb  = +1.00 r1_1_bb(a,i) l0_1
            // flops: o1v1L2  = o1v1L2
            //  mems: o1v1L2  = o1v1L2
            D_ov_bb("L,R,i,a")  = l0_1("L") * r1_1["bb"]("R,a,i");

            // D_vo_aa  = +1.00 r0 l1_aa(i,a)
            // flops: o1v1L2  = o1v1L2
            //  mems: o1v1L2  = o1v1L2
            D_vo_aa("L,R,a,i")  = l1["aa"]("L,i,a") * r0("R");

            // D_ov_aa  = +1.00 r1_1_aa(a,i) l0_1
            // flops: o1v1L2  = o1v1L2
            //  mems: o1v1L2  = o1v1L2
            D_ov_aa("L,R,i,a")  = l0_1("L") * r1_1["aa"]("R,a,i");

            // D_vv_bb  = +1.00 r1_1_bb(b,i) l1_1_bb(i,a)
            // flops: o0v2L2  = o1v2L2
            //  mems: o0v2L2  = o0v2L2
            D_vv_bb("L,R,a,b")  = l1_1["bb"]("L,i,a") * r1_1["bb"]("R,b,i");

            // D_vo_bb  = +1.00 r0 l1_bb(i,a)
            // flops: o1v1L2  = o1v1L2
            //  mems: o1v1L2  = o1v1L2
            D_vo_bb("L,R,a,i")  = l1["bb"]("L,i,a") * r0("R");

            // D_vv_aa  = +1.00 r1_1_aa(b,i) l1_1_aa(i,a)
            // flops: o0v2L2  = o1v2L2
            //  mems: o0v2L2  = o0v2L2
            D_vv_aa("L,R,a,b")  = l1_1["aa"]("L,i,a") * r1_1["aa"]("R,b,i");

            // D_ov_aa += +1.00 r1_aa(a,i) l0
            // flops: o1v1L2 += o1v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") += l0("L") * r1["aa"]("R,a,i");

            // D_ov_aa += +1.00 r2_abab(a,b,i,j) l1_bb(j,b)
            // flops: o1v1L2 += o2v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") += l1["bb"]("L,j,b") * r2["abab"]("R,a,b,i,j");

            // D_ov_aa += +1.00 r2_1_abab(a,b,i,j) l1_1_bb(j,b)
            // flops: o1v1L2 += o2v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") += l1_1["bb"]("L,j,b") * r2_1["abab"]("R,a,b,i,j");

            // D_ov_aa += -1.00 r2_aaaa(a,b,j,i) l1_aa(j,b)
            // flops: o1v1L2 += o2v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") -= l1["aa"]("L,j,b") * r2["aaaa"]("R,a,b,j,i");

            // D_ov_aa += -1.00 r2_1_aaaa(a,b,j,i) l1_1_aa(j,b)
            // flops: o1v1L2 += o2v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") -= l1_1["aa"]("L,j,b") * r2_1["aaaa"]("R,a,b,j,i");

            // D_ov_aa += +1.00 r0 t2_abab(a,b,i,j) l1_bb(j,b)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") += l1["bb"]("L,j,b") * t2["abab"]("a,b,i,j") * r0("R");

            // D_ov_aa += +1.00 r0 t2_1_abab(a,b,i,j) l1_1_bb(j,b)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") += l1_1["bb"]("L,j,b") * t2_1["abab"]("a,b,i,j") * r0("R");

            // D_ov_aa += +1.00 r0_1 t2_abab(a,b,i,j) l1_1_bb(j,b)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") += l1_1["bb"]("L,j,b") * t2["abab"]("a,b,i,j") * r0_1("R");

            // D_ov_aa += -1.00 r0 t2_1_aaaa(a,b,j,i) l1_1_aa(j,b)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= l1_1["aa"]("L,j,b") * t2_1["aaaa"]("a,b,j,i") * r0("R");

            // D_ov_aa += -1.00 r0_1 t2_aaaa(a,b,j,i) l1_1_aa(j,b)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= l1_1["aa"]("L,j,b") * t2["aaaa"]("a,b,j,i") * r0_1("R");

            // D_ov_aa += -1.00 r0 t2_aaaa(a,b,j,i) l1_aa(j,b)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= l1["aa"]("L,j,b") * t2["aaaa"]("a,b,j,i") * r0("R");

            // D_ov_aa += -1.00 r1_aa(b,i) t1_1_aa(a,j) l1_1_aa(j,b)
            // flops: o1v1L2 += o2v1L2 o2v1L2
            //  mems: o1v1L2 += o2v0L2 o1v1L2
            D_ov_aa("L,R,i,a") -= l1_1["aa"]("L,j,b") * r1["aa"]("R,b,i") * t1_1["aa"]("a,j");

            // D_ov_aa += +0.50 r1_aa(a,k) t2_1_aaaa(c,b,j,i) l2_1_aaaa(k,j,c,b)
            // flops: o1v1L2 += o3v2L1 o2v1L2
            //  mems: o1v1L2 += o2v0L1 o1v1L2
            D_ov_aa("L,R,i,a") += 0.50 * l2_1["aaaa"]("L,k,j,c,b") * t2_1["aaaa"]("c,b,j,i") * r1["aa"]("R,a,k");

            // D_ov_aa += +0.50 r1_aa(a,k) t2_aaaa(c,b,j,i) l2_aaaa(k,j,c,b)
            // flops: o1v1L2 += o3v2L1 o2v1L2
            //  mems: o1v1L2 += o2v0L1 o1v1L2
            D_ov_aa("L,R,i,a") += 0.50 * l2["aaaa"]("L,k,j,c,b") * t2["aaaa"]("c,b,j,i") * r1["aa"]("R,a,k");

            // D_ov_aa += +0.50 r1_aa(c,i) t2_aaaa(a,b,j,k) l2_aaaa(j,k,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o0v2L1 o1v1L2
            D_ov_aa("L,R,i,a") += 0.50 * l2["aaaa"]("L,j,k,b,c") * t2["aaaa"]("a,b,j,k") * r1["aa"]("R,c,i");

            // D_ov_aa += +0.50 r1_aa(c,i) t2_1_aaaa(a,b,j,k) l2_1_aaaa(j,k,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o0v2L1 o1v1L2
            D_ov_aa("L,R,i,a") += 0.50 * l2_1["aaaa"]("L,j,k,b,c") * t2_1["aaaa"]("a,b,j,k") * r1["aa"]("R,c,i");

            // D_ov_aa += -0.50 r2_abab(b,c,i,k) t1_1_aa(a,j) l2_1_abab(j,k,b,c)
            //         += -0.50 r2_abab(c,b,i,k) t1_1_aa(a,j) l2_1_abab(j,k,c,b)
            // flops: o1v1L2 += o3v2L2 o2v1L2
            //  mems: o1v1L2 += o2v0L2 o1v1L2
            D_ov_aa("L,R,i,a") -= l2_1["abab"]("L,j,k,b,c") * r2["abab"]("R,b,c,i,k") * t1_1["aa"]("a,j");

            // D_ov_aa += -0.50 r2_abab(a,c,k,j) t1_1_aa(b,i) l2_1_abab(k,j,b,c)
            //         += -0.50 r2_abab(a,c,j,k) t1_1_aa(b,i) l2_1_abab(j,k,b,c)
            // flops: o1v1L2 += o3v2L1 o3v2L2
            //  mems: o1v1L2 += o3v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= l2_1["abab"]("L,k,j,b,c") * t1_1["aa"]("b,i") * r2["abab"]("R,a,c,k,j");

            // D_ov_aa += -0.50 r2_aaaa(b,c,k,i) t1_1_aa(a,j) l2_1_aaaa(k,j,b,c)
            // flops: o1v1L2 += o3v2L2 o2v1L2
            //  mems: o1v1L2 += o2v0L2 o1v1L2
            D_ov_aa("L,R,i,a") -= 0.50 * l2_1["aaaa"]("L,k,j,b,c") * r2["aaaa"]("R,b,c,k,i") * t1_1["aa"]("a,j");

            // D_ov_aa += -0.50 r2_aaaa(a,c,k,j) t1_1_aa(b,i) l2_1_aaaa(k,j,b,c)
            // flops: o1v1L2 += o3v2L1 o3v2L2
            //  mems: o1v1L2 += o3v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= 0.50 * l2_1["aaaa"]("L,k,j,b,c") * t1_1["aa"]("b,i") * r2["aaaa"]("R,a,c,k,j");

            // D_ov_bb += +1.00 r1_bb(a,i) l0
            // flops: o1v1L2 += o1v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") += l0("L") * r1["bb"]("R,a,i");

            // D_ov_bb += -1.00 r2_1_bbbb(a,b,j,i) l1_1_bb(j,b)
            // flops: o1v1L2 += o2v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") -= l1_1["bb"]("L,j,b") * r2_1["bbbb"]("R,a,b,j,i");

            // D_ov_bb += -1.00 r2_bbbb(a,b,j,i) l1_bb(j,b)
            // flops: o1v1L2 += o2v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") -= l1["bb"]("L,j,b") * r2["bbbb"]("R,a,b,j,i");

            // D_ov_bb += +1.00 r2_abab(b,a,j,i) l1_aa(j,b)
            // flops: o1v1L2 += o2v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") += l1["aa"]("L,j,b") * r2["abab"]("R,b,a,j,i");

            // D_ov_bb += +1.00 r2_1_abab(b,a,j,i) l1_1_aa(j,b)
            // flops: o1v1L2 += o2v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") += l1_1["aa"]("L,j,b") * r2_1["abab"]("R,b,a,j,i");

            // D_ov_bb += -1.00 r0 t2_1_bbbb(a,b,j,i) l1_1_bb(j,b)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l1_1["bb"]("L,j,b") * t2_1["bbbb"]("a,b,j,i") * r0("R");

            // D_ov_bb += -1.00 r0_1 t2_bbbb(a,b,j,i) l1_1_bb(j,b)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l1_1["bb"]("L,j,b") * t2["bbbb"]("a,b,j,i") * r0_1("R");

            // D_ov_bb += -1.00 r0 t2_bbbb(a,b,j,i) l1_bb(j,b)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l1["bb"]("L,j,b") * t2["bbbb"]("a,b,j,i") * r0("R");

            // D_ov_bb += +1.00 r0 t2_1_abab(b,a,j,i) l1_1_aa(j,b)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") += l1_1["aa"]("L,j,b") * t2_1["abab"]("b,a,j,i") * r0("R");

            // D_ov_bb += +1.00 r0 t2_abab(b,a,j,i) l1_aa(j,b)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") += l1["aa"]("L,j,b") * t2["abab"]("b,a,j,i") * r0("R");

            // D_ov_bb += +1.00 r0_1 t2_abab(b,a,j,i) l1_1_aa(j,b)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") += l1_1["aa"]("L,j,b") * t2["abab"]("b,a,j,i") * r0_1("R");

            // D_ov_bb += -1.00 r1_bb(b,i) t1_1_bb(a,j) l1_1_bb(j,b)
            // flops: o1v1L2 += o2v1L2 o2v1L2
            //  mems: o1v1L2 += o2v0L2 o1v1L2
            D_ov_bb("L,R,i,a") -= l1_1["bb"]("L,j,b") * r1["bb"]("R,b,i") * t1_1["bb"]("a,j");

            // D_ov_bb += +0.50 r1_bb(a,k) t2_bbbb(c,b,j,i) l2_bbbb(k,j,c,b)
            // flops: o1v1L2 += o3v2L1 o2v1L2
            //  mems: o1v1L2 += o2v0L1 o1v1L2
            D_ov_bb("L,R,i,a") += 0.50 * l2["bbbb"]("L,k,j,c,b") * t2["bbbb"]("c,b,j,i") * r1["bb"]("R,a,k");

            // D_ov_bb += +0.50 r1_bb(a,k) t2_1_bbbb(c,b,j,i) l2_1_bbbb(k,j,c,b)
            // flops: o1v1L2 += o3v2L1 o2v1L2
            //  mems: o1v1L2 += o2v0L1 o1v1L2
            D_ov_bb("L,R,i,a") += 0.50 * l2_1["bbbb"]("L,k,j,c,b") * t2_1["bbbb"]("c,b,j,i") * r1["bb"]("R,a,k");

            // D_ov_bb += +0.50 r1_bb(c,i) t2_bbbb(a,b,j,k) l2_bbbb(j,k,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o0v2L1 o1v1L2
            D_ov_bb("L,R,i,a") += 0.50 * l2["bbbb"]("L,j,k,b,c") * t2["bbbb"]("a,b,j,k") * r1["bb"]("R,c,i");

            // D_ov_bb += +0.50 r1_bb(c,i) t2_1_bbbb(a,b,j,k) l2_1_bbbb(j,k,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o0v2L1 o1v1L2
            D_ov_bb("L,R,i,a") += 0.50 * l2_1["bbbb"]("L,j,k,b,c") * t2_1["bbbb"]("a,b,j,k") * r1["bb"]("R,c,i");

            // D_ov_bb += -0.50 r2_bbbb(b,c,k,i) t1_1_bb(a,j) l2_1_bbbb(k,j,b,c)
            // flops: o1v1L2 += o3v2L2 o2v1L2
            //  mems: o1v1L2 += o2v0L2 o1v1L2
            D_ov_bb("L,R,i,a") -= 0.50 * l2_1["bbbb"]("L,k,j,b,c") * r2["bbbb"]("R,b,c,k,i") * t1_1["bb"]("a,j");

            // D_ov_bb += -0.50 r2_bbbb(a,c,k,j) t1_1_bb(b,i) l2_1_bbbb(k,j,b,c)
            // flops: o1v1L2 += o3v2L1 o3v2L2
            //  mems: o1v1L2 += o3v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= 0.50 * l2_1["bbbb"]("L,k,j,b,c") * t1_1["bb"]("b,i") * r2["bbbb"]("R,a,c,k,j");

            // D_ov_bb += -0.50 r2_abab(b,c,k,i) t1_1_bb(a,j) l2_1_abab(k,j,b,c)
            //         += -0.50 r2_abab(c,b,k,i) t1_1_bb(a,j) l2_1_abab(k,j,c,b)
            // flops: o1v1L2 += o3v2L2 o2v1L2
            //  mems: o1v1L2 += o2v0L2 o1v1L2
            D_ov_bb("L,R,i,a") -= l2_1["abab"]("L,k,j,b,c") * r2["abab"]("R,b,c,k,i") * t1_1["bb"]("a,j");

            // D_ov_bb += -0.50 r2_abab(c,a,k,j) t1_1_bb(b,i) l2_1_abab(k,j,c,b)
            //         += -0.50 r2_abab(c,a,j,k) t1_1_bb(b,i) l2_1_abab(j,k,c,b)
            // flops: o1v1L2 += o3v2L1 o3v2L2
            //  mems: o1v1L2 += o3v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l2_1["abab"]("L,k,j,c,b") * t1_1["bb"]("b,i") * r2["abab"]("R,c,a,k,j");

            // D_vo_aa += +1.00 r0_1 l1_1_aa(i,a)
            // flops: o1v1L2 += o1v1L2
            //  mems: o1v1L2 += o1v1L2
            D_vo_aa("L,R,a,i") += l1_1["aa"]("L,i,a") * r0_1("R");

            // D_vo_bb += +1.00 r0_1 l1_1_bb(i,a)
            // flops: o1v1L2 += o1v1L2
            //  mems: o1v1L2 += o1v1L2
            D_vo_bb("L,R,a,i") += l1_1["bb"]("L,i,a") * r0_1("R");

            // D_vv_aa += +1.00 r1_aa(b,i) l1_aa(i,a)
            // flops: o0v2L2 += o1v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_aa("L,R,a,b") += l1["aa"]("L,i,a") * r1["aa"]("R,b,i");

            // D_vv_aa += +1.00 r0_1 t1_aa(b,i) l1_1_aa(i,a)
            // flops: o0v2L2 += o1v2L1 o0v2L2
            //  mems: o0v2L2 += o0v2L1 o0v2L2
            D_vv_aa("L,R,a,b") += l1_1["aa"]("L,i,a") * t1["aa"]("b,i") * r0_1("R");

            // D_vv_aa += +1.00 r0 t1_aa(b,i) l1_aa(i,a)
            // flops: o0v2L2 += o1v2L1 o0v2L2
            //  mems: o0v2L2 += o0v2L1 o0v2L2
            D_vv_aa("L,R,a,b") += l1["aa"]("L,i,a") * t1["aa"]("b,i") * r0("R");

            // D_vv_aa += +1.00 r0 t1_1_aa(b,i) l1_1_aa(i,a)
            // flops: o0v2L2 += o1v2L1 o0v2L2
            //  mems: o0v2L2 += o0v2L1 o0v2L2
            D_vv_aa("L,R,a,b") += l1_1["aa"]("L,i,a") * t1_1["aa"]("b,i") * r0("R");

            // D_vv_bb += +1.00 r1_bb(b,i) l1_bb(i,a)
            // flops: o0v2L2 += o1v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_bb("L,R,a,b") += l1["bb"]("L,i,a") * r1["bb"]("R,b,i");

            // D_vv_bb += +1.00 r0 t1_bb(b,i) l1_bb(i,a)
            // flops: o0v2L2 += o1v2L1 o0v2L2
            //  mems: o0v2L2 += o0v2L1 o0v2L2
            D_vv_bb("L,R,a,b") += l1["bb"]("L,i,a") * t1["bb"]("b,i") * r0("R");

            // D_vv_bb += +1.00 r0_1 t1_bb(b,i) l1_1_bb(i,a)
            // flops: o0v2L2 += o1v2L1 o0v2L2
            //  mems: o0v2L2 += o0v2L1 o0v2L2
            D_vv_bb("L,R,a,b") += l1_1["bb"]("L,i,a") * t1["bb"]("b,i") * r0_1("R");

            // D_vv_bb += +1.00 r0 t1_1_bb(b,i) l1_1_bb(i,a)
            // flops: o0v2L2 += o1v2L1 o0v2L2
            //  mems: o0v2L2 += o0v2L1 o0v2L2
            D_vv_bb("L,R,a,b") += l1_1["bb"]("L,i,a") * t1_1["bb"]("b,i") * r0("R");

            // flops: o0v2L1  = o2v3L1
            //  mems: o0v2L1  = o0v2L1
            tmps_["1_aa_Lvv"]("L,a,b")  = l2_1["abab"]("L,i,j,a,c") * t2["abab"]("b,c,i,j");

            // D_ov_aa += -0.50 r1_1_aa(c,i) t2_abab(a,b,j,k) l2_1_abab(j,k,c,b)
            //         += -0.50 r1_1_aa(c,i) t2_abab(a,b,k,j) l2_1_abab(k,j,c,b)
            // flops: o1v1L2 += o1v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") -= r1_1["aa"]("R,c,i") * tmps_["1_aa_Lvv"]("L,c,a");

            // D_ov_aa += -0.50 r0 t2_abab(a,c,j,k) t1_1_aa(b,i) l2_1_abab(j,k,b,c)
            //         += -0.50 r0 t2_abab(a,c,k,j) t1_1_aa(b,i) l2_1_abab(k,j,b,c)
            // flops: o1v1L2 += o1v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= t1_1["aa"]("b,i") * tmps_["1_aa_Lvv"]("L,b,a") * r0("R");

            // D_ov_aa += -0.50 r0_1 t2_abab(a,b,j,k) t1_aa(c,i) l2_1_abab(j,k,c,b)
            //         += -0.50 r0_1 t2_abab(a,b,k,j) t1_aa(c,i) l2_1_abab(k,j,c,b)
            // flops: o1v1L2 += o1v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= t1["aa"]("c,i") * tmps_["1_aa_Lvv"]("L,c,a") * r0_1("R");

            // D_vv_aa += +0.50 r0_1 t2_abab(b,c,i,j) l2_1_abab(i,j,a,c)
            //         += +0.50 r0_1 t2_abab(b,c,j,i) l2_1_abab(j,i,a,c)
            // flops: o0v2L2 += o0v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_aa("L,R,a,b") += r0_1("R") * tmps_["1_aa_Lvv"]("L,a,b");
            tmps_["1_aa_Lvv"].~TArrayD();

            // flops: o0v2L1  = o2v3L1
            //  mems: o0v2L1  = o0v2L1
            tmps_["2_bb_Lvv"]("L,a,b")  = l2_1["abab"]("L,i,j,c,b") * t2["abab"]("c,a,i,j");

            // D_ov_bb += -0.50 r1_1_bb(c,i) t2_abab(b,a,j,k) l2_1_abab(j,k,b,c)
            //         += -0.50 r1_1_bb(c,i) t2_abab(b,a,k,j) l2_1_abab(k,j,b,c)
            // flops: o1v1L2 += o1v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") -= r1_1["bb"]("R,c,i") * tmps_["2_bb_Lvv"]("L,a,c");

            // D_ov_bb += -0.50 r0 t2_abab(c,a,j,k) t1_1_bb(b,i) l2_1_abab(j,k,c,b)
            //         += -0.50 r0 t2_abab(c,a,k,j) t1_1_bb(b,i) l2_1_abab(k,j,c,b)
            // flops: o1v1L2 += o1v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= t1_1["bb"]("b,i") * tmps_["2_bb_Lvv"]("L,a,b") * r0("R");

            // D_ov_bb += -0.50 r0_1 t2_abab(b,a,j,k) t1_bb(c,i) l2_1_abab(j,k,b,c)
            //         += -0.50 r0_1 t2_abab(b,a,k,j) t1_bb(c,i) l2_1_abab(k,j,b,c)
            // flops: o1v1L2 += o1v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= t1["bb"]("c,i") * tmps_["2_bb_Lvv"]("L,a,c") * r0_1("R");

            // D_vv_bb += +0.50 r0_1 t2_abab(c,b,i,j) l2_1_abab(i,j,c,a)
            //         += +0.50 r0_1 t2_abab(c,b,j,i) l2_1_abab(j,i,c,a)
            // flops: o0v2L2 += o0v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_bb("L,R,a,b") += r0_1("R") * tmps_["2_bb_Lvv"]("L,b,a");
            tmps_["2_bb_Lvv"].~TArrayD();

            // flops: o0v2L1  = o2v3L1
            //  mems: o0v2L1  = o0v2L1
            tmps_["3_aa_Lvv"]("L,a,b")  = l2_1["aaaa"]("L,i,j,c,b") * t2["aaaa"]("a,c,i,j");

            // D_ov_aa += +0.50 r1_1_aa(c,i) t2_aaaa(a,b,j,k) l2_1_aaaa(j,k,b,c)
            // flops: o1v1L2 += o1v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") += 0.50 * r1_1["aa"]("R,c,i") * tmps_["3_aa_Lvv"]("L,a,c");

            // D_ov_aa += +0.50 r0 t2_aaaa(a,c,j,k) t1_1_aa(b,i) l2_1_aaaa(j,k,c,b)
            // flops: o1v1L2 += o1v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") += 0.50 * t1_1["aa"]("b,i") * tmps_["3_aa_Lvv"]("L,a,b") * r0("R");
            tmps_["3_aa_Lvv"].~TArrayD();

            // flops: o0v2L1  = o2v3L1
            //  mems: o0v2L1  = o0v2L1
            tmps_["4_aa_Lvv"]("L,a,b")  = l2_1["aaaa"]("L,i,j,a,c") * t2["aaaa"]("b,c,i,j");

            // D_ov_aa += -0.50 r0_1 t2_aaaa(a,b,j,k) t1_aa(c,i) l2_1_aaaa(j,k,c,b)
            // flops: o1v1L2 += o1v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= 0.50 * t1["aa"]("c,i") * tmps_["4_aa_Lvv"]("L,c,a") * r0_1("R");

            // D_vv_aa += +0.50 r0_1 t2_aaaa(b,c,i,j) l2_1_aaaa(i,j,a,c)
            // flops: o0v2L2 += o0v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_aa("L,R,a,b") += 0.50 * r0_1("R") * tmps_["4_aa_Lvv"]("L,a,b");
            tmps_["4_aa_Lvv"].~TArrayD();

            // flops: o2v0L1  = o3v2L1
            //  mems: o2v0L1  = o2v0L1
            tmps_["5_aa_Loo"]("L,i,j")  = l2_1["aaaa"]("L,j,k,a,b") * t2["aaaa"]("a,b,k,i");

            // D_ov_aa += +0.50 r1_1_aa(a,k) t2_aaaa(c,b,j,i) l2_1_aaaa(k,j,c,b)
            // flops: o1v1L2 += o2v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") += 0.50 * r1_1["aa"]("R,a,k") * tmps_["5_aa_Loo"]("L,i,k");

            // D_ov_aa += +0.50 r0 t2_aaaa(c,b,k,i) t1_1_aa(a,j) l2_1_aaaa(j,k,c,b)
            // flops: o1v1L2 += o2v1L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") += 0.50 * t1_1["aa"]("a,j") * tmps_["5_aa_Loo"]("L,i,j") * r0("R");
            tmps_["5_aa_Loo"].~TArrayD();

            // flops: o2v0L1  = o3v2L1
            //  mems: o2v0L1  = o2v0L1
            tmps_["6_aa_Loo"]("L,i,j")  = l2_1["aaaa"]("L,k,j,a,b") * t2["aaaa"]("a,b,k,i");

            // D_oo_aa  = -0.50 r0_1 t2_aaaa(b,a,k,i) l2_1_aaaa(k,j,b,a)
            // flops: o2v0L2  = o2v0L2
            //  mems: o2v0L2  = o2v0L2
            D_oo_aa("L,R,i,j")  = -0.50 * r0_1("R") * tmps_["6_aa_Loo"]("L,i,j");

            // D_ov_aa += -0.50 r0_1 t1_aa(a,k) t2_aaaa(c,b,j,i) l2_1_aaaa(j,k,c,b)
            // flops: o1v1L2 += o2v1L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= 0.50 * t1["aa"]("a,k") * tmps_["6_aa_Loo"]("L,i,k") * r0_1("R");
            tmps_["6_aa_Loo"].~TArrayD();

            // flops: o2v0L1  = o3v2L1
            //  mems: o2v0L1  = o2v0L1
            tmps_["7_bb_Loo"]("L,i,j")  = l2_1["abab"]("L,k,j,a,b") * t2["abab"]("a,b,k,i");

            // D_oo_bb  = -0.50 r0_1 t2_abab(b,a,k,i) l2_1_abab(k,j,b,a)
            //         += -0.50 r0_1 t2_abab(a,b,k,i) l2_1_abab(k,j,a,b)
            // flops: o2v0L2  = o2v0L2
            //  mems: o2v0L2  = o2v0L2
            D_oo_bb("L,R,i,j")  = -1.00 * r0_1("R") * tmps_["7_bb_Loo"]("L,i,j");

            // D_ov_bb += -0.50 r1_1_bb(a,k) t2_abab(c,b,j,i) l2_1_abab(j,k,c,b)
            //         += -0.50 r1_1_bb(a,k) t2_abab(b,c,j,i) l2_1_abab(j,k,b,c)
            // flops: o1v1L2 += o2v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") -= r1_1["bb"]("R,a,k") * tmps_["7_bb_Loo"]("L,i,k");

            // D_ov_bb += -0.50 r0_1 t1_bb(a,k) t2_abab(c,b,j,i) l2_1_abab(j,k,c,b)
            //         += -0.50 r0_1 t1_bb(a,k) t2_abab(b,c,j,i) l2_1_abab(j,k,b,c)
            // flops: o1v1L2 += o2v1L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= t1["bb"]("a,k") * tmps_["7_bb_Loo"]("L,i,k") * r0_1("R");

            // D_ov_bb += -0.50 r0 t2_abab(c,b,k,i) t1_1_bb(a,j) l2_1_abab(k,j,c,b)
            //         += -0.50 r0 t2_abab(b,c,k,i) t1_1_bb(a,j) l2_1_abab(k,j,b,c)
            // flops: o1v1L2 += o2v1L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= t1_1["bb"]("a,j") * tmps_["7_bb_Loo"]("L,i,j") * r0("R");
            tmps_["7_bb_Loo"].~TArrayD();

            // flops: o2v0L1  = o3v2L1
            //  mems: o2v0L1  = o2v0L1
            tmps_["8_aa_Loo"]("L,i,j")  = l2_1["abab"]("L,j,k,a,b") * t2["abab"]("a,b,i,k");

            // D_oo_aa += -0.50 r0_1 t2_abab(b,a,i,k) l2_1_abab(j,k,b,a)
            //         += -0.50 r0_1 t2_abab(a,b,i,k) l2_1_abab(j,k,a,b)
            // flops: o2v0L2 += o2v0L2
            //  mems: o2v0L2 += o2v0L2
            D_oo_aa("L,R,i,j") -= r0_1("R") * tmps_["8_aa_Loo"]("L,i,j");

            // D_ov_aa += -0.50 r1_1_aa(a,k) t2_abab(c,b,i,j) l2_1_abab(k,j,c,b)
            //         += -0.50 r1_1_aa(a,k) t2_abab(b,c,i,j) l2_1_abab(k,j,b,c)
            // flops: o1v1L2 += o2v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") -= r1_1["aa"]("R,a,k") * tmps_["8_aa_Loo"]("L,i,k");

            // D_ov_aa += -0.50 r0 t2_abab(c,b,i,k) t1_1_aa(a,j) l2_1_abab(j,k,c,b)
            //         += -0.50 r0 t2_abab(b,c,i,k) t1_1_aa(a,j) l2_1_abab(j,k,b,c)
            // flops: o1v1L2 += o2v1L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= t1_1["aa"]("a,j") * tmps_["8_aa_Loo"]("L,i,j") * r0("R");

            // D_ov_aa += -0.50 r0_1 t1_aa(a,k) t2_abab(c,b,i,j) l2_1_abab(k,j,c,b)
            //         += -0.50 r0_1 t1_aa(a,k) t2_abab(b,c,i,j) l2_1_abab(k,j,b,c)
            // flops: o1v1L2 += o2v1L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= t1["aa"]("a,k") * tmps_["8_aa_Loo"]("L,i,k") * r0_1("R");
            tmps_["8_aa_Loo"].~TArrayD();

            // flops: o0v2L1  = o2v3L1
            //  mems: o0v2L1  = o0v2L1
            tmps_["9_bb_Lvv"]("L,a,b")  = l2_1["bbbb"]("L,i,j,c,b") * t2["bbbb"]("a,c,i,j");

            // D_ov_bb += +0.50 r1_1_bb(c,i) t2_bbbb(a,b,j,k) l2_1_bbbb(j,k,b,c)
            // flops: o1v1L2 += o1v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") += 0.50 * r1_1["bb"]("R,c,i") * tmps_["9_bb_Lvv"]("L,a,c");

            // D_ov_bb += +0.50 r0 t2_bbbb(a,c,j,k) t1_1_bb(b,i) l2_1_bbbb(j,k,c,b)
            // flops: o1v1L2 += o1v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") += 0.50 * t1_1["bb"]("b,i") * tmps_["9_bb_Lvv"]("L,a,b") * r0("R");
            tmps_["9_bb_Lvv"].~TArrayD();

            // flops: o0v2L1  = o2v3L1
            //  mems: o0v2L1  = o0v2L1
            tmps_["10_bb_Lvv"]("L,a,b")  = l2_1["bbbb"]("L,i,j,a,c") * t2["bbbb"]("b,c,i,j");

            // D_ov_bb += -0.50 r0_1 t2_bbbb(a,b,j,k) t1_bb(c,i) l2_1_bbbb(j,k,c,b)
            // flops: o1v1L2 += o1v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= 0.50 * t1["bb"]("c,i") * tmps_["10_bb_Lvv"]("L,c,a") * r0_1("R");

            // D_vv_bb += +0.50 r0_1 t2_bbbb(b,c,i,j) l2_1_bbbb(i,j,a,c)
            // flops: o0v2L2 += o0v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_bb("L,R,a,b") += 0.50 * r0_1("R") * tmps_["10_bb_Lvv"]("L,a,b");
            tmps_["10_bb_Lvv"].~TArrayD();

            // flops: o2v0L1  = o3v2L1
            //  mems: o2v0L1  = o2v0L1
            tmps_["11_bb_Loo"]("L,i,j")  = l2_1["bbbb"]("L,j,k,a,b") * t2["bbbb"]("a,b,k,i");

            // D_ov_bb += +0.50 r1_1_bb(a,k) t2_bbbb(c,b,j,i) l2_1_bbbb(k,j,c,b)
            // flops: o1v1L2 += o2v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") += 0.50 * r1_1["bb"]("R,a,k") * tmps_["11_bb_Loo"]("L,i,k");

            // D_ov_bb += +0.50 r0 t2_bbbb(c,b,k,i) t1_1_bb(a,j) l2_1_bbbb(j,k,c,b)
            // flops: o1v1L2 += o2v1L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") += 0.50 * t1_1["bb"]("a,j") * tmps_["11_bb_Loo"]("L,i,j") * r0("R");
            tmps_["11_bb_Loo"].~TArrayD();

            // flops: o2v0L1  = o3v2L1
            //  mems: o2v0L1  = o2v0L1
            tmps_["12_bb_Loo"]("L,i,j")  = l2_1["bbbb"]("L,k,j,a,b") * t2["bbbb"]("a,b,k,i");

            // D_oo_bb += -0.50 r0_1 t2_bbbb(b,a,k,i) l2_1_bbbb(k,j,b,a)
            // flops: o2v0L2 += o2v0L2
            //  mems: o2v0L2 += o2v0L2
            D_oo_bb("L,R,i,j") -= 0.50 * r0_1("R") * tmps_["12_bb_Loo"]("L,i,j");

            // D_ov_bb += -0.50 r0_1 t1_bb(a,k) t2_bbbb(c,b,j,i) l2_1_bbbb(j,k,c,b)
            // flops: o1v1L2 += o2v1L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= 0.50 * t1["bb"]("a,k") * tmps_["12_bb_Loo"]("L,i,k") * r0_1("R");
            tmps_["12_bb_Loo"].~TArrayD();

            // flops: o2v0L1  = o2v1L1
            //  mems: o2v0L1  = o2v0L1
            tmps_["13_aa_Loo"]("L,i,j")  = l1_1["aa"]("L,j,a") * t1["aa"]("a,i");

            // D_oo_aa += -1.00 r0_1 t1_aa(a,i) l1_1_aa(j,a)
            // flops: o2v0L2 += o2v0L2
            //  mems: o2v0L2 += o2v0L2
            D_oo_aa("L,R,i,j") -= r0_1("R") * tmps_["13_aa_Loo"]("L,i,j");

            // D_ov_aa += -1.00 r1_1_aa(a,j) t1_aa(b,i) l1_1_aa(j,b)
            // flops: o1v1L2 += o2v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") -= r1_1["aa"]("R,a,j") * tmps_["13_aa_Loo"]("L,i,j");

            // D_ov_aa += -1.00 r0_1 t1_aa(a,j) t1_aa(b,i) l1_1_aa(j,b)
            // flops: o1v1L2 += o2v1L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= t1["aa"]("a,j") * tmps_["13_aa_Loo"]("L,i,j") * r0_1("R");

            // D_ov_aa += -1.00 r0 t1_aa(b,i) t1_1_aa(a,j) l1_1_aa(j,b)
            // flops: o1v1L2 += o2v1L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= t1_1["aa"]("a,j") * tmps_["13_aa_Loo"]("L,i,j") * r0("R");
            tmps_["13_aa_Loo"].~TArrayD();

            // flops: o2v0L1  = o2v1L1
            //  mems: o2v0L1  = o2v0L1
            tmps_["14_bb_Loo"]("L,i,j")  = l1_1["bb"]("L,j,a") * t1["bb"]("a,i");

            // D_oo_bb += -1.00 r0_1 t1_bb(a,i) l1_1_bb(j,a)
            // flops: o2v0L2 += o2v0L2
            //  mems: o2v0L2 += o2v0L2
            D_oo_bb("L,R,i,j") -= r0_1("R") * tmps_["14_bb_Loo"]("L,i,j");

            // D_ov_bb += -1.00 r1_1_bb(a,j) t1_bb(b,i) l1_1_bb(j,b)
            // flops: o1v1L2 += o2v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") -= r1_1["bb"]("R,a,j") * tmps_["14_bb_Loo"]("L,i,j");

            // D_ov_bb += -1.00 r0 t1_bb(b,i) t1_1_bb(a,j) l1_1_bb(j,b)
            // flops: o1v1L2 += o2v1L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= t1_1["bb"]("a,j") * tmps_["14_bb_Loo"]("L,i,j") * r0("R");

            // D_ov_bb += -1.00 r0_1 t1_bb(a,j) t1_bb(b,i) l1_1_bb(j,b)
            // flops: o1v1L2 += o2v1L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= t1["bb"]("a,j") * tmps_["14_bb_Loo"]("L,i,j") * r0_1("R");
            tmps_["14_bb_Loo"].~TArrayD();

            // flops: o0v2L2  = o2v3L2 o2v3L2 o2v3L2 o2v3L2 o0v2L2 o0v2L2 o0v2L2
            //  mems: o0v2L2  = o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2
            tmps_["15_bb_LLvv"]("L,R,a,b")  = 0.50 * l2["bbbb"]("L,i,j,a,c") * r2["bbbb"]("R,b,c,i,j");
            tmps_["15_bb_LLvv"]("L,R,a,b") += 0.50 * l2_1["bbbb"]("L,i,j,a,c") * r2_1["bbbb"]("R,b,c,i,j");
            tmps_["15_bb_LLvv"]("L,R,a,b") += l2["abab"]("L,k,j,d,a") * r2["abab"]("R,d,b,k,j");
            tmps_["15_bb_LLvv"]("L,R,a,b") += l2_1["abab"]("L,k,j,d,a") * r2_1["abab"]("R,d,b,k,j");

            // D_ov_bb += -0.50 r2_1_abab(c,a,k,j) t1_bb(b,i) l2_1_abab(k,j,c,b)
            //         += -0.50 r2_1_abab(c,a,j,k) t1_bb(b,i) l2_1_abab(j,k,c,b)
            //         += -0.50 r2_abab(c,a,k,j) t1_bb(b,i) l2_abab(k,j,c,b)
            //         += -0.50 r2_abab(c,a,j,k) t1_bb(b,i) l2_abab(j,k,c,b)
            //         += -0.50 r2_1_bbbb(a,c,k,j) t1_bb(b,i) l2_1_bbbb(k,j,b,c)
            //         += -0.50 r2_bbbb(a,c,k,j) t1_bb(b,i) l2_bbbb(k,j,b,c)
            // flops: o1v1L2 += o1v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") -= t1["bb"]("b,i") * tmps_["15_bb_LLvv"]("L,R,b,a");

            // D_vv_bb += +0.50 r2_1_abab(c,b,j,i) l2_1_abab(j,i,c,a)
            //         += +0.50 r2_1_abab(c,b,i,j) l2_1_abab(i,j,c,a)
            //         += +0.50 r2_abab(c,b,j,i) l2_abab(j,i,c,a)
            //         += +0.50 r2_abab(c,b,i,j) l2_abab(i,j,c,a)
            //         += +0.50 r2_1_bbbb(b,c,j,i) l2_1_bbbb(j,i,a,c)
            //         += +0.50 r2_bbbb(b,c,j,i) l2_bbbb(j,i,a,c)
            D_vv_bb("L,R,a,b") += tmps_["15_bb_LLvv"]("L,R,a,b");
            tmps_["15_bb_LLvv"].~TArrayD();

            // flops: o0v2L2  = o2v3L2 o2v3L2 o0v2L2 o2v3L2 o0v2L2 o2v3L2 o0v2L2
            //  mems: o0v2L2  = o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2
            tmps_["16_aa_LLvv"]("L,R,a,b")  = 0.50 * l2["aaaa"]("L,i,k,a,d") * r2["aaaa"]("R,b,d,i,k");
            tmps_["16_aa_LLvv"]("L,R,a,b") += l2_1["abab"]("L,i,j,a,c") * r2_1["abab"]("R,b,c,i,j");
            tmps_["16_aa_LLvv"]("L,R,a,b") += 0.50 * l2_1["aaaa"]("L,i,k,a,d") * r2_1["aaaa"]("R,b,d,i,k");
            tmps_["16_aa_LLvv"]("L,R,a,b") += l2["abab"]("L,i,j,a,c") * r2["abab"]("R,b,c,i,j");

            // D_ov_aa += -0.50 r2_1_abab(a,c,k,j) t1_aa(b,i) l2_1_abab(k,j,b,c)
            //         += -0.50 r2_1_abab(a,c,j,k) t1_aa(b,i) l2_1_abab(j,k,b,c)
            //         += -0.50 r2_aaaa(a,c,k,j) t1_aa(b,i) l2_aaaa(k,j,b,c)
            //         += -0.50 r2_1_aaaa(a,c,k,j) t1_aa(b,i) l2_1_aaaa(k,j,b,c)
            //         += -0.50 r2_abab(a,c,k,j) t1_aa(b,i) l2_abab(k,j,b,c)
            //         += -0.50 r2_abab(a,c,j,k) t1_aa(b,i) l2_abab(j,k,b,c)
            // flops: o1v1L2 += o1v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") -= t1["aa"]("b,i") * tmps_["16_aa_LLvv"]("L,R,b,a");

            // D_vv_aa += +0.50 r2_1_abab(b,c,j,i) l2_1_abab(j,i,a,c)
            //         += +0.50 r2_1_abab(b,c,i,j) l2_1_abab(i,j,a,c)
            //         += +0.50 r2_aaaa(b,c,j,i) l2_aaaa(j,i,a,c)
            //         += +0.50 r2_1_aaaa(b,c,j,i) l2_1_aaaa(j,i,a,c)
            //         += +0.50 r2_abab(b,c,j,i) l2_abab(j,i,a,c)
            //         += +0.50 r2_abab(b,c,i,j) l2_abab(i,j,a,c)
            D_vv_aa("L,R,a,b") += tmps_["16_aa_LLvv"]("L,R,a,b");
            tmps_["16_aa_LLvv"].~TArrayD();

            // flops: o0v0L2  = o0v0L2 o1v1L2 o1v1L2 o2v2L2 o2v2L2 o0v0L2 o0v0L2 o2v2L2 o0v0L2 o0v0L2 o0v0L2
            //  mems: o0v0L2  = o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2
            tmps_["17_LL"]("L,R")  = 4.00 * l0_1("L") * r0("R");
            tmps_["17_LL"]("L,R") += 4.00 * l1_1["aa"]("L,k,c") * r1["aa"]("R,c,k");
            tmps_["17_LL"]("L,R") += 4.00 * l1_1["bb"]("L,j,a") * r1["bb"]("R,a,j");
            tmps_["17_LL"]("L,R") += 4.00 * l2_1["abab"]("L,l,j,c,b") * r2["abab"]("R,c,b,l,j");
            tmps_["17_LL"]("L,R") += l2_1["aaaa"]("L,l,k,c,d") * r2["aaaa"]("R,c,d,l,k");
            tmps_["17_LL"]("L,R") += l2_1["bbbb"]("L,i,j,a,b") * r2["bbbb"]("R,a,b,i,j");

            // D_ov_aa += +0.25 r2_aaaa(b,c,k,j) t1_1_aa(a,i) l2_1_aaaa(k,j,b,c)
            //         += +0.25 r2_abab(b,c,k,j) t1_1_aa(a,i) l2_1_abab(k,j,b,c)
            //         += +0.25 r2_abab(b,c,j,k) t1_1_aa(a,i) l2_1_abab(j,k,b,c)
            //         += +0.25 r2_abab(c,b,k,j) t1_1_aa(a,i) l2_1_abab(k,j,c,b)
            //         += +0.25 r2_abab(c,b,j,k) t1_1_aa(a,i) l2_1_abab(j,k,c,b)
            //         += +1.00 r1_bb(b,j) t1_1_aa(a,i) l1_1_bb(j,b)
            //         += +0.25 r2_bbbb(b,c,k,j) t1_1_aa(a,i) l2_1_bbbb(k,j,b,c)
            //         += +1.00 r1_aa(b,j) t1_1_aa(a,i) l1_1_aa(j,b)
            //         += +1.00 r0 t1_1_aa(a,i) l0_1
            // flops: o1v1L2 += o1v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") += 0.25 * t1_1["aa"]("a,i") * tmps_["17_LL"]("L,R");

            // D_ov_bb += +0.25 r2_aaaa(b,c,k,j) t1_1_bb(a,i) l2_1_aaaa(k,j,b,c)
            //         += +0.25 r2_abab(b,c,k,j) t1_1_bb(a,i) l2_1_abab(k,j,b,c)
            //         += +0.25 r2_abab(b,c,j,k) t1_1_bb(a,i) l2_1_abab(j,k,b,c)
            //         += +0.25 r2_abab(c,b,k,j) t1_1_bb(a,i) l2_1_abab(k,j,c,b)
            //         += +0.25 r2_abab(c,b,j,k) t1_1_bb(a,i) l2_1_abab(j,k,c,b)
            //         += +1.00 r1_bb(b,j) t1_1_bb(a,i) l1_1_bb(j,b)
            //         += +0.25 r2_bbbb(b,c,k,j) t1_1_bb(a,i) l2_1_bbbb(k,j,b,c)
            //         += +1.00 r1_aa(b,j) t1_1_bb(a,i) l1_1_aa(j,b)
            //         += +1.00 r0 t1_1_bb(a,i) l0_1
            // flops: o1v1L2 += o1v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") += 0.25 * t1_1["bb"]("a,i") * tmps_["17_LL"]("L,R");
            tmps_["17_LL"].~TArrayD();

            // flops: o0v2L1  = o2v3L1 o2v3L1 o0v2L1
            //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
            tmps_["18_aa_Lvv"]("L,a,b")  = l2["abab"]("L,i,j,a,c") * t2["abab"]("b,c,i,j");
            tmps_["18_aa_Lvv"]("L,a,b") += l2_1["abab"]("L,i,j,a,c") * t2_1["abab"]("b,c,i,j");

            // D_ov_aa += -0.50 r1_aa(c,i) t2_1_abab(a,b,j,k) l2_1_abab(j,k,c,b)
            //         += -0.50 r1_aa(c,i) t2_1_abab(a,b,k,j) l2_1_abab(k,j,c,b)
            //         += -0.50 r1_aa(c,i) t2_abab(a,b,j,k) l2_abab(j,k,c,b)
            //         += -0.50 r1_aa(c,i) t2_abab(a,b,k,j) l2_abab(k,j,c,b)
            // flops: o1v1L2 += o1v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") -= r1["aa"]("R,c,i") * tmps_["18_aa_Lvv"]("L,c,a");

            // D_ov_aa += -0.50 r0 t1_aa(c,i) t2_1_abab(a,b,j,k) l2_1_abab(j,k,c,b)
            //         += -0.50 r0 t1_aa(c,i) t2_1_abab(a,b,k,j) l2_1_abab(k,j,c,b)
            //         += -0.50 r0 t2_abab(a,b,j,k) t1_aa(c,i) l2_abab(j,k,c,b)
            //         += -0.50 r0 t2_abab(a,b,k,j) t1_aa(c,i) l2_abab(k,j,c,b)
            // flops: o1v1L2 += o1v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= t1["aa"]("c,i") * tmps_["18_aa_Lvv"]("L,c,a") * r0("R");

            // D_vv_aa += +0.50 r0 t2_1_abab(b,c,i,j) l2_1_abab(i,j,a,c)
            //         += +0.50 r0 t2_1_abab(b,c,j,i) l2_1_abab(j,i,a,c)
            //         += +0.50 r0 t2_abab(b,c,i,j) l2_abab(i,j,a,c)
            //         += +0.50 r0 t2_abab(b,c,j,i) l2_abab(j,i,a,c)
            // flops: o0v2L2 += o0v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_aa("L,R,a,b") += r0("R") * tmps_["18_aa_Lvv"]("L,a,b");
            tmps_["18_aa_Lvv"].~TArrayD();

            // flops: o0v2L1  = o2v3L1 o2v3L1 o0v2L1
            //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
            tmps_["19_bb_Lvv"]("L,a,b")  = l2["abab"]("L,i,j,c,a") * t2["abab"]("c,b,i,j");
            tmps_["19_bb_Lvv"]("L,a,b") += l2_1["abab"]("L,i,j,c,a") * t2_1["abab"]("c,b,i,j");

            // D_ov_bb += -0.50 r1_bb(c,i) t2_abab(b,a,j,k) l2_abab(j,k,b,c)
            //         += -0.50 r1_bb(c,i) t2_abab(b,a,k,j) l2_abab(k,j,b,c)
            //         += -0.50 r1_bb(c,i) t2_1_abab(b,a,j,k) l2_1_abab(j,k,b,c)
            //         += -0.50 r1_bb(c,i) t2_1_abab(b,a,k,j) l2_1_abab(k,j,b,c)
            // flops: o1v1L2 += o1v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") -= r1["bb"]("R,c,i") * tmps_["19_bb_Lvv"]("L,c,a");

            // D_ov_bb += -0.50 r0 t2_abab(b,a,j,k) t1_bb(c,i) l2_abab(j,k,b,c)
            //         += -0.50 r0 t2_abab(b,a,k,j) t1_bb(c,i) l2_abab(k,j,b,c)
            //         += -0.50 r0 t1_bb(c,i) t2_1_abab(b,a,j,k) l2_1_abab(j,k,b,c)
            //         += -0.50 r0 t1_bb(c,i) t2_1_abab(b,a,k,j) l2_1_abab(k,j,b,c)
            // flops: o1v1L2 += o1v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= t1["bb"]("c,i") * tmps_["19_bb_Lvv"]("L,c,a") * r0("R");

            // D_vv_bb += +0.50 r0 t2_abab(c,b,i,j) l2_abab(i,j,c,a)
            //         += +0.50 r0 t2_abab(c,b,j,i) l2_abab(j,i,c,a)
            //         += +0.50 r0 t2_1_abab(c,b,i,j) l2_1_abab(i,j,c,a)
            //         += +0.50 r0 t2_1_abab(c,b,j,i) l2_1_abab(j,i,c,a)
            // flops: o0v2L2 += o0v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_bb("L,R,a,b") += r0("R") * tmps_["19_bb_Lvv"]("L,a,b");
            tmps_["19_bb_Lvv"].~TArrayD();

            // flops: o0v2L1  = o2v3L1 o2v3L1 o0v2L1
            //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
            tmps_["20_aa_Lvv"]("L,a,b")  = l2["aaaa"]("L,i,j,a,c") * t2["aaaa"]("b,c,i,j");
            tmps_["20_aa_Lvv"]("L,a,b") += l2_1["aaaa"]("L,i,j,a,c") * t2_1["aaaa"]("b,c,i,j");

            // D_ov_aa += -0.50 r0 t2_aaaa(a,b,j,k) t1_aa(c,i) l2_aaaa(j,k,c,b)
            //         += -0.50 r0 t1_aa(c,i) t2_1_aaaa(a,b,j,k) l2_1_aaaa(j,k,c,b)
            // flops: o1v1L2 += o1v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= 0.50 * t1["aa"]("c,i") * tmps_["20_aa_Lvv"]("L,c,a") * r0("R");

            // D_vv_aa += +0.50 r0 t2_aaaa(b,c,i,j) l2_aaaa(i,j,a,c)
            //         += +0.50 r0 t2_1_aaaa(b,c,i,j) l2_1_aaaa(i,j,a,c)
            // flops: o0v2L2 += o0v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_aa("L,R,a,b") += 0.50 * r0("R") * tmps_["20_aa_Lvv"]("L,a,b");
            tmps_["20_aa_Lvv"].~TArrayD();

            // flops: o2v0L1  = o3v2L1 o3v2L1 o2v0L1
            //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
            tmps_["21_aa_Loo"]("L,i,j")  = l2["aaaa"]("L,k,j,a,b") * t2["aaaa"]("a,b,k,i");
            tmps_["21_aa_Loo"]("L,i,j") += l2_1["aaaa"]("L,k,j,a,b") * t2_1["aaaa"]("a,b,k,i");

            // D_oo_aa += -0.50 r0 t2_aaaa(b,a,k,i) l2_aaaa(k,j,b,a)
            //         += -0.50 r0 t2_1_aaaa(b,a,k,i) l2_1_aaaa(k,j,b,a)
            // flops: o2v0L2 += o2v0L2
            //  mems: o2v0L2 += o2v0L2
            D_oo_aa("L,R,i,j") -= 0.50 * r0("R") * tmps_["21_aa_Loo"]("L,i,j");

            // D_ov_aa += -0.50 r0 t1_aa(a,k) t2_aaaa(c,b,j,i) l2_aaaa(j,k,c,b)
            //         += -0.50 r0 t1_aa(a,k) t2_1_aaaa(c,b,j,i) l2_1_aaaa(j,k,c,b)
            // flops: o1v1L2 += o2v1L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= 0.50 * t1["aa"]("a,k") * tmps_["21_aa_Loo"]("L,i,k") * r0("R");
            tmps_["21_aa_Loo"].~TArrayD();

            // flops: o2v0L1  = o3v2L1 o3v2L1 o2v0L1
            //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
            tmps_["22_bb_Loo"]("L,i,j")  = l2["abab"]("L,k,j,a,b") * t2["abab"]("a,b,k,i");
            tmps_["22_bb_Loo"]("L,i,j") += l2_1["abab"]("L,k,j,a,b") * t2_1["abab"]("a,b,k,i");

            // D_oo_bb += -0.50 r0 t2_1_abab(b,a,k,i) l2_1_abab(k,j,b,a)
            //         += -0.50 r0 t2_1_abab(a,b,k,i) l2_1_abab(k,j,a,b)
            //         += -0.50 r0 t2_abab(b,a,k,i) l2_abab(k,j,b,a)
            //         += -0.50 r0 t2_abab(a,b,k,i) l2_abab(k,j,a,b)
            // flops: o2v0L2 += o2v0L2
            //  mems: o2v0L2 += o2v0L2
            D_oo_bb("L,R,i,j") -= r0("R") * tmps_["22_bb_Loo"]("L,i,j");

            // D_ov_bb += -0.50 r1_bb(a,k) t2_1_abab(c,b,j,i) l2_1_abab(j,k,c,b)
            //         += -0.50 r1_bb(a,k) t2_1_abab(b,c,j,i) l2_1_abab(j,k,b,c)
            //         += -0.50 r1_bb(a,k) t2_abab(c,b,j,i) l2_abab(j,k,c,b)
            //         += -0.50 r1_bb(a,k) t2_abab(b,c,j,i) l2_abab(j,k,b,c)
            // flops: o1v1L2 += o2v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") -= r1["bb"]("R,a,k") * tmps_["22_bb_Loo"]("L,i,k");

            // D_ov_bb += -0.50 r0 t1_bb(a,k) t2_1_abab(c,b,j,i) l2_1_abab(j,k,c,b)
            //         += -0.50 r0 t1_bb(a,k) t2_1_abab(b,c,j,i) l2_1_abab(j,k,b,c)
            //         += -0.50 r0 t1_bb(a,k) t2_abab(c,b,j,i) l2_abab(j,k,c,b)
            //         += -0.50 r0 t1_bb(a,k) t2_abab(b,c,j,i) l2_abab(j,k,b,c)
            // flops: o1v1L2 += o2v1L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= t1["bb"]("a,k") * tmps_["22_bb_Loo"]("L,i,k") * r0("R");
            tmps_["22_bb_Loo"].~TArrayD();

            // flops: o2v0L1  = o3v2L1 o3v2L1 o2v0L1
            //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
            tmps_["23_aa_Loo"]("L,i,j")  = l2["abab"]("L,j,k,a,b") * t2["abab"]("a,b,i,k");
            tmps_["23_aa_Loo"]("L,i,j") += l2_1["abab"]("L,j,k,a,b") * t2_1["abab"]("a,b,i,k");

            // D_oo_aa += -0.50 r0 t2_1_abab(b,a,i,k) l2_1_abab(j,k,b,a)
            //         += -0.50 r0 t2_1_abab(a,b,i,k) l2_1_abab(j,k,a,b)
            //         += -0.50 r0 t2_abab(b,a,i,k) l2_abab(j,k,b,a)
            //         += -0.50 r0 t2_abab(a,b,i,k) l2_abab(j,k,a,b)
            // flops: o2v0L2 += o2v0L2
            //  mems: o2v0L2 += o2v0L2
            D_oo_aa("L,R,i,j") -= r0("R") * tmps_["23_aa_Loo"]("L,i,j");

            // D_ov_aa += -0.50 r1_aa(a,k) t2_1_abab(c,b,i,j) l2_1_abab(k,j,c,b)
            //         += -0.50 r1_aa(a,k) t2_1_abab(b,c,i,j) l2_1_abab(k,j,b,c)
            //         += -0.50 r1_aa(a,k) t2_abab(c,b,i,j) l2_abab(k,j,c,b)
            //         += -0.50 r1_aa(a,k) t2_abab(b,c,i,j) l2_abab(k,j,b,c)
            // flops: o1v1L2 += o2v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") -= r1["aa"]("R,a,k") * tmps_["23_aa_Loo"]("L,i,k");

            // D_ov_aa += -0.50 r0 t1_aa(a,k) t2_1_abab(c,b,i,j) l2_1_abab(k,j,c,b)
            //         += -0.50 r0 t1_aa(a,k) t2_1_abab(b,c,i,j) l2_1_abab(k,j,b,c)
            //         += -0.50 r0 t1_aa(a,k) t2_abab(c,b,i,j) l2_abab(k,j,c,b)
            //         += -0.50 r0 t1_aa(a,k) t2_abab(b,c,i,j) l2_abab(k,j,b,c)
            // flops: o1v1L2 += o2v1L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= t1["aa"]("a,k") * tmps_["23_aa_Loo"]("L,i,k") * r0("R");
            tmps_["23_aa_Loo"].~TArrayD();

            // flops: o0v2L1  = o2v3L1 o2v3L1 o0v2L1
            //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
            tmps_["24_bb_Lvv"]("L,a,b")  = l2["bbbb"]("L,i,j,a,c") * t2["bbbb"]("b,c,i,j");
            tmps_["24_bb_Lvv"]("L,a,b") += l2_1["bbbb"]("L,i,j,a,c") * t2_1["bbbb"]("b,c,i,j");

            // D_ov_bb += -0.50 r0 t1_bb(c,i) t2_1_bbbb(a,b,j,k) l2_1_bbbb(j,k,c,b)
            //         += -0.50 r0 t2_bbbb(a,b,j,k) t1_bb(c,i) l2_bbbb(j,k,c,b)
            // flops: o1v1L2 += o1v2L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= 0.50 * t1["bb"]("c,i") * tmps_["24_bb_Lvv"]("L,c,a") * r0("R");

            // D_vv_bb += +0.50 r0 t2_1_bbbb(b,c,i,j) l2_1_bbbb(i,j,a,c)
            //         += +0.50 r0 t2_bbbb(b,c,i,j) l2_bbbb(i,j,a,c)
            // flops: o0v2L2 += o0v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_bb("L,R,a,b") += 0.50 * r0("R") * tmps_["24_bb_Lvv"]("L,a,b");
            tmps_["24_bb_Lvv"].~TArrayD();

            // flops: o2v0L1  = o3v2L1 o3v2L1 o2v0L1
            //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
            tmps_["25_bb_Loo"]("L,i,j")  = l2["bbbb"]("L,k,j,a,b") * t2["bbbb"]("a,b,k,i");
            tmps_["25_bb_Loo"]("L,i,j") += l2_1["bbbb"]("L,k,j,a,b") * t2_1["bbbb"]("a,b,k,i");

            // D_oo_bb += -0.50 r0 t2_1_bbbb(b,a,k,i) l2_1_bbbb(k,j,b,a)
            //         += -0.50 r0 t2_bbbb(b,a,k,i) l2_bbbb(k,j,b,a)
            // flops: o2v0L2 += o2v0L2
            //  mems: o2v0L2 += o2v0L2
            D_oo_bb("L,R,i,j") -= 0.50 * r0("R") * tmps_["25_bb_Loo"]("L,i,j");

            // D_ov_bb += -0.50 r0 t1_bb(a,k) t2_1_bbbb(c,b,j,i) l2_1_bbbb(j,k,c,b)
            //         += -0.50 r0 t1_bb(a,k) t2_bbbb(c,b,j,i) l2_bbbb(j,k,c,b)
            // flops: o1v1L2 += o2v1L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= 0.50 * t1["bb"]("a,k") * tmps_["25_bb_Loo"]("L,i,k") * r0("R");
            tmps_["25_bb_Loo"].~TArrayD();

            // flops: o2v0L1  = o2v1L1 o2v1L1 o2v0L1
            //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
            tmps_["26_aa_Loo"]("L,i,j")  = l1["aa"]("L,j,a") * t1["aa"]("a,i");
            tmps_["26_aa_Loo"]("L,i,j") += l1_1["aa"]("L,j,a") * t1_1["aa"]("a,i");

            // D_oo_aa += -1.00 r0 t1_aa(a,i) l1_aa(j,a)
            //         += -1.00 r0 t1_1_aa(a,i) l1_1_aa(j,a)
            // flops: o2v0L2 += o2v0L2
            //  mems: o2v0L2 += o2v0L2
            D_oo_aa("L,R,i,j") -= r0("R") * tmps_["26_aa_Loo"]("L,i,j");

            // D_ov_aa += -1.00 r1_aa(a,j) t1_aa(b,i) l1_aa(j,b)
            //         += -1.00 r1_aa(a,j) t1_1_aa(b,i) l1_1_aa(j,b)
            // flops: o1v1L2 += o2v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") -= r1["aa"]("R,a,j") * tmps_["26_aa_Loo"]("L,i,j");

            // D_ov_aa += -1.00 r0 t1_aa(a,j) t1_aa(b,i) l1_aa(j,b)
            //         += -1.00 r0 t1_aa(a,j) t1_1_aa(b,i) l1_1_aa(j,b)
            // flops: o1v1L2 += o2v1L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= t1["aa"]("a,j") * tmps_["26_aa_Loo"]("L,i,j") * r0("R");
            tmps_["26_aa_Loo"].~TArrayD();

            // flops: o2v0L1  = o2v1L1 o2v1L1 o2v0L1
            //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
            tmps_["27_bb_Loo"]("L,i,j")  = l1["bb"]("L,j,a") * t1["bb"]("a,i");
            tmps_["27_bb_Loo"]("L,i,j") += l1_1["bb"]("L,j,a") * t1_1["bb"]("a,i");

            // D_oo_bb += -1.00 r0 t1_bb(a,i) l1_bb(j,a)
            //         += -1.00 r0 t1_1_bb(a,i) l1_1_bb(j,a)
            // flops: o2v0L2 += o2v0L2
            //  mems: o2v0L2 += o2v0L2
            D_oo_bb("L,R,i,j") -= r0("R") * tmps_["27_bb_Loo"]("L,i,j");

            // D_ov_bb += -1.00 r1_bb(a,j) t1_bb(b,i) l1_bb(j,b)
            //         += -1.00 r1_bb(a,j) t1_1_bb(b,i) l1_1_bb(j,b)
            // flops: o1v1L2 += o2v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") -= r1["bb"]("R,a,j") * tmps_["27_bb_Loo"]("L,i,j");

            // D_ov_bb += -1.00 r0 t1_bb(a,j) t1_bb(b,i) l1_bb(j,b)
            //         += -1.00 r0 t1_bb(a,j) t1_1_bb(b,i) l1_1_bb(j,b)
            // flops: o1v1L2 += o2v1L1 o1v1L2
            //  mems: o1v1L2 += o1v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= t1["bb"]("a,j") * tmps_["27_bb_Loo"]("L,i,j") * r0("R");
            tmps_["27_bb_Loo"].~TArrayD();

            // flops: o0v0L2  = o2v2L2 o0v0L2 o2v2L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o2v2L2 o0v0L2 o1v1L2 o0v0L2 o1v1L2 o0v0L2 o1v1L2 o0v0L2 o2v2L2 o0v0L2 o1v1L2 o0v0L2 o2v2L2 o0v0L2 o2v2L2 o0v0L2
            //  mems: o0v0L2  = o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2
            tmps_["28_LL"]("L,R")  = 0.25 * l2["aaaa"]("L,i,k,a,d") * r2["aaaa"]("R,a,d,i,k");
            tmps_["28_LL"]("L,R") += l0("L") * r0("R");
            tmps_["28_LL"]("L,R") += l2_1["abab"]("L,i,j,a,b") * r2_1["abab"]("R,a,b,i,j");
            tmps_["28_LL"]("L,R") += l0_1("L") * r0_1("R");
            tmps_["28_LL"]("L,R") += 0.25 * l2["bbbb"]("L,l,j,c,b") * r2["bbbb"]("R,c,b,l,j");
            tmps_["28_LL"]("L,R") += l1_1["aa"]("L,k,a") * r1_1["aa"]("R,a,k");
            tmps_["28_LL"]("L,R") += l1["aa"]("L,k,a") * r1["aa"]("R,a,k");
            tmps_["28_LL"]("L,R") += l1_1["bb"]("L,j,c") * r1_1["bb"]("R,c,j");
            tmps_["28_LL"]("L,R") += 0.25 * l2_1["bbbb"]("L,l,j,c,b") * r2_1["bbbb"]("R,c,b,l,j");
            tmps_["28_LL"]("L,R") += l1["bb"]("L,j,c") * r1["bb"]("R,c,j");
            tmps_["28_LL"]("L,R") += 0.25 * l2_1["aaaa"]("L,i,k,a,d") * r2_1["aaaa"]("R,a,d,i,k");
            tmps_["28_LL"]("L,R") += l2["abab"]("L,i,j,a,b") * r2["abab"]("R,a,b,i,j");

            // D_oo_aa += +1.00 d_aa(i,j) r0 l0
            //         += +0.25 d_aa(i,j) r2_1_abab(a,b,l,k) l2_1_abab(l,k,a,b)
            //         += +0.25 d_aa(i,j) r2_1_abab(a,b,k,l) l2_1_abab(k,l,a,b)
            //         += +0.25 d_aa(i,j) r2_1_abab(b,a,l,k) l2_1_abab(l,k,b,a)
            //         += +0.25 d_aa(i,j) r2_1_abab(b,a,k,l) l2_1_abab(k,l,b,a)
            //         += +0.25 d_aa(i,j) r2_aaaa(a,b,l,k) l2_aaaa(l,k,a,b)
            //         += +1.00 d_aa(i,j) r0_1 l0_1
            //         += +0.25 d_aa(i,j) r2_bbbb(a,b,l,k) l2_bbbb(l,k,a,b)
            //         += +1.00 d_aa(i,j) r1_1_aa(a,k) l1_1_aa(k,a)
            //         += +1.00 d_aa(i,j) r1_aa(a,k) l1_aa(k,a)
            //         += +1.00 d_aa(i,j) r1_1_bb(a,k) l1_1_bb(k,a)
            //         += +0.25 d_aa(i,j) r2_1_bbbb(a,b,l,k) l2_1_bbbb(l,k,a,b)
            //         += +1.00 d_aa(i,j) r1_bb(a,k) l1_bb(k,a)
            //         += +0.25 d_aa(i,j) r2_1_aaaa(a,b,l,k) l2_1_aaaa(l,k,a,b)
            //         += +0.25 d_aa(i,j) r2_abab(a,b,l,k) l2_abab(l,k,a,b)
            //         += +0.25 d_aa(i,j) r2_abab(a,b,k,l) l2_abab(k,l,a,b)
            //         += +0.25 d_aa(i,j) r2_abab(b,a,l,k) l2_abab(l,k,b,a)
            //         += +0.25 d_aa(i,j) r2_abab(b,a,k,l) l2_abab(k,l,b,a)
            // flops: o2v0L2 += o2v0L2
            //  mems: o2v0L2 += o2v0L2
            D_oo_aa("L,R,i,j") += Id["aa_oo"]("i,j") * tmps_["28_LL"]("L,R");

            // D_oo_bb += +1.00 d_bb(i,j) r0 l0
            //         += +0.25 d_bb(i,j) r2_1_abab(a,b,l,k) l2_1_abab(l,k,a,b)
            //         += +0.25 d_bb(i,j) r2_1_abab(a,b,k,l) l2_1_abab(k,l,a,b)
            //         += +0.25 d_bb(i,j) r2_1_abab(b,a,l,k) l2_1_abab(l,k,b,a)
            //         += +0.25 d_bb(i,j) r2_1_abab(b,a,k,l) l2_1_abab(k,l,b,a)
            //         += +0.25 d_bb(i,j) r2_aaaa(a,b,l,k) l2_aaaa(l,k,a,b)
            //         += +1.00 d_bb(i,j) r0_1 l0_1
            //         += +0.25 d_bb(i,j) r2_bbbb(a,b,l,k) l2_bbbb(l,k,a,b)
            //         += +1.00 d_bb(i,j) r1_1_aa(a,k) l1_1_aa(k,a)
            //         += +1.00 d_bb(i,j) r1_aa(a,k) l1_aa(k,a)
            //         += +1.00 d_bb(i,j) r1_1_bb(a,k) l1_1_bb(k,a)
            //         += +0.25 d_bb(i,j) r2_1_bbbb(a,b,l,k) l2_1_bbbb(l,k,a,b)
            //         += +1.00 d_bb(i,j) r1_bb(a,k) l1_bb(k,a)
            //         += +0.25 d_bb(i,j) r2_1_aaaa(a,b,l,k) l2_1_aaaa(l,k,a,b)
            //         += +0.25 d_bb(i,j) r2_abab(a,b,l,k) l2_abab(l,k,a,b)
            //         += +0.25 d_bb(i,j) r2_abab(a,b,k,l) l2_abab(k,l,a,b)
            //         += +0.25 d_bb(i,j) r2_abab(b,a,l,k) l2_abab(l,k,b,a)
            //         += +0.25 d_bb(i,j) r2_abab(b,a,k,l) l2_abab(k,l,b,a)
            // flops: o2v0L2 += o2v0L2
            //  mems: o2v0L2 += o2v0L2
            D_oo_bb("L,R,i,j") += Id["bb_oo"]("i,j") * tmps_["28_LL"]("L,R");

            // D_ov_aa += +1.00 r0 t1_aa(a,i) l0
            //         += +0.25 r2_1_abab(b,c,k,j) t1_aa(a,i) l2_1_abab(k,j,b,c)
            //         += +0.25 r2_1_abab(b,c,j,k) t1_aa(a,i) l2_1_abab(j,k,b,c)
            //         += +0.25 r2_1_abab(c,b,k,j) t1_aa(a,i) l2_1_abab(k,j,c,b)
            //         += +0.25 r2_1_abab(c,b,j,k) t1_aa(a,i) l2_1_abab(j,k,c,b)
            //         += +0.25 r2_aaaa(b,c,k,j) t1_aa(a,i) l2_aaaa(k,j,b,c)
            //         += +1.00 r0_1 t1_aa(a,i) l0_1
            //         += +0.25 r2_bbbb(b,c,k,j) t1_aa(a,i) l2_bbbb(k,j,b,c)
            //         += +1.00 r1_1_aa(b,j) t1_aa(a,i) l1_1_aa(j,b)
            //         += +1.00 r1_aa(b,j) t1_aa(a,i) l1_aa(j,b)
            //         += +1.00 r1_1_bb(b,j) t1_aa(a,i) l1_1_bb(j,b)
            //         += +0.25 r2_1_bbbb(b,c,k,j) t1_aa(a,i) l2_1_bbbb(k,j,b,c)
            //         += +1.00 r1_bb(b,j) t1_aa(a,i) l1_bb(j,b)
            //         += +0.25 r2_1_aaaa(b,c,k,j) t1_aa(a,i) l2_1_aaaa(k,j,b,c)
            //         += +0.25 r2_abab(b,c,k,j) t1_aa(a,i) l2_abab(k,j,b,c)
            //         += +0.25 r2_abab(b,c,j,k) t1_aa(a,i) l2_abab(j,k,b,c)
            //         += +0.25 r2_abab(c,b,k,j) t1_aa(a,i) l2_abab(k,j,c,b)
            //         += +0.25 r2_abab(c,b,j,k) t1_aa(a,i) l2_abab(j,k,c,b)
            // flops: o1v1L2 += o1v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") += t1["aa"]("a,i") * tmps_["28_LL"]("L,R");

            // D_ov_bb += +1.00 r0 t1_bb(a,i) l0
            //         += +0.25 r2_1_abab(b,c,k,j) t1_bb(a,i) l2_1_abab(k,j,b,c)
            //         += +0.25 r2_1_abab(b,c,j,k) t1_bb(a,i) l2_1_abab(j,k,b,c)
            //         += +0.25 r2_1_abab(c,b,k,j) t1_bb(a,i) l2_1_abab(k,j,c,b)
            //         += +0.25 r2_1_abab(c,b,j,k) t1_bb(a,i) l2_1_abab(j,k,c,b)
            //         += +0.25 r2_aaaa(b,c,k,j) t1_bb(a,i) l2_aaaa(k,j,b,c)
            //         += +1.00 r0_1 t1_bb(a,i) l0_1
            //         += +0.25 r2_bbbb(b,c,k,j) t1_bb(a,i) l2_bbbb(k,j,b,c)
            //         += +1.00 r1_1_aa(b,j) t1_bb(a,i) l1_1_aa(j,b)
            //         += +1.00 r1_aa(b,j) t1_bb(a,i) l1_aa(j,b)
            //         += +1.00 r1_1_bb(b,j) t1_bb(a,i) l1_1_bb(j,b)
            //         += +0.25 r2_1_bbbb(b,c,k,j) t1_bb(a,i) l2_1_bbbb(k,j,b,c)
            //         += +1.00 r1_bb(b,j) t1_bb(a,i) l1_bb(j,b)
            //         += +0.25 r2_1_aaaa(b,c,k,j) t1_bb(a,i) l2_1_aaaa(k,j,b,c)
            //         += +0.25 r2_abab(b,c,k,j) t1_bb(a,i) l2_abab(k,j,b,c)
            //         += +0.25 r2_abab(b,c,j,k) t1_bb(a,i) l2_abab(j,k,b,c)
            //         += +0.25 r2_abab(c,b,k,j) t1_bb(a,i) l2_abab(k,j,c,b)
            //         += +0.25 r2_abab(c,b,j,k) t1_bb(a,i) l2_abab(j,k,c,b)
            // flops: o1v1L2 += o1v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") += t1["bb"]("a,i") * tmps_["28_LL"]("L,R");
            tmps_["28_LL"].~TArrayD();

            // flops: o1v1L2  = o2v2L2 o2v2L2 o2v2L2 o2v2L2 o1v1L2 o1v1L2 o1v1L2
            //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2
            tmps_["29_bb_LLov"]("L,R,i,a")  = -1.00 * l2["bbbb"]("L,j,i,a,b") * r1["bb"]("R,b,j");
            tmps_["29_bb_LLov"]("L,R,i,a") -= l2_1["bbbb"]("L,j,i,a,b") * r1_1["bb"]("R,b,j");
            tmps_["29_bb_LLov"]("L,R,i,a") += l2["abab"]("L,k,i,c,a") * r1["aa"]("R,c,k");
            tmps_["29_bb_LLov"]("L,R,i,a") += l2_1["abab"]("L,k,i,c,a") * r1_1["aa"]("R,c,k");

            // D_ov_aa += +1.00 r1_aa(c,k) t2_abab(a,b,i,j) l2_abab(k,j,c,b)
            //         += +1.00 r1_1_aa(c,k) t2_abab(a,b,i,j) l2_1_abab(k,j,c,b)
            //         += -1.00 r1_1_bb(c,k) t2_abab(a,b,i,j) l2_1_bbbb(k,j,b,c)
            //         += -1.00 r1_bb(c,k) t2_abab(a,b,i,j) l2_bbbb(k,j,b,c)
            // flops: o1v1L2 += o2v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") += t2["abab"]("a,b,i,j") * tmps_["29_bb_LLov"]("L,R,j,b");

            // D_ov_bb += -1.00 r1_aa(c,k) t2_bbbb(a,b,j,i) l2_abab(k,j,c,b)
            //         += -1.00 r1_1_aa(c,k) t2_bbbb(a,b,j,i) l2_1_abab(k,j,c,b)
            //         += +1.00 r1_1_bb(c,k) t2_bbbb(a,b,j,i) l2_1_bbbb(k,j,b,c)
            //         += +1.00 r1_bb(c,k) t2_bbbb(a,b,j,i) l2_bbbb(k,j,b,c)
            // flops: o1v1L2 += o2v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") -= t2["bbbb"]("a,b,j,i") * tmps_["29_bb_LLov"]("L,R,j,b");

            // D_vo_bb += +1.00 r1_aa(b,j) l2_abab(j,i,b,a)
            //         += +1.00 r1_1_aa(b,j) l2_1_abab(j,i,b,a)
            //         += -1.00 r1_1_bb(b,j) l2_1_bbbb(j,i,a,b)
            //         += -1.00 r1_bb(b,j) l2_bbbb(j,i,a,b)
            D_vo_bb("L,R,a,i") += tmps_["29_bb_LLov"]("L,R,i,a");

            // D_vv_bb += +1.00 r1_aa(c,j) t1_bb(b,i) l2_abab(j,i,c,a)
            //         += +1.00 r1_1_aa(c,j) t1_bb(b,i) l2_1_abab(j,i,c,a)
            //         += -1.00 r1_1_bb(c,j) t1_bb(b,i) l2_1_bbbb(j,i,a,c)
            //         += -1.00 r1_bb(c,j) t1_bb(b,i) l2_bbbb(j,i,a,c)
            // flops: o0v2L2 += o1v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_bb("L,R,a,b") += t1["bb"]("b,i") * tmps_["29_bb_LLov"]("L,R,i,a");

            // flops: o1v1L2  = o2v2L2 o2v2L2 o1v1L2
            //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2
            tmps_["30_bb_LLov"]("L,R,i,a")  = -1.00 * l2_1["bbbb"]("L,k,i,a,c") * r1["bb"]("R,c,k");
            tmps_["30_bb_LLov"]("L,R,i,a") += l2_1["abab"]("L,j,i,b,a") * r1["aa"]("R,b,j");

            // D_ov_aa += +1.00 r1_aa(c,k) t2_1_abab(a,b,i,j) l2_1_abab(k,j,c,b)
            //         += -1.00 r1_bb(c,k) t2_1_abab(a,b,i,j) l2_1_bbbb(k,j,b,c)
            // flops: o1v1L2 += o2v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") += t2_1["abab"]("a,b,i,j") * tmps_["30_bb_LLov"]("L,R,j,b");

            // D_ov_bb += -1.00 r1_aa(c,k) t2_1_bbbb(a,b,j,i) l2_1_abab(k,j,c,b)
            //         += +1.00 r1_bb(c,k) t2_1_bbbb(a,b,j,i) l2_1_bbbb(k,j,b,c)
            // flops: o1v1L2 += o2v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") -= t2_1["bbbb"]("a,b,j,i") * tmps_["30_bb_LLov"]("L,R,j,b");

            // D_ov_bb += -1.00 r1_aa(c,k) t1_bb(b,i) t1_1_bb(a,j) l2_1_abab(k,j,c,b)
            //         += +1.00 r1_bb(c,k) t1_bb(b,i) t1_1_bb(a,j) l2_1_bbbb(k,j,b,c)
            // flops: o1v1L2 += o2v1L2 o2v1L2
            //  mems: o1v1L2 += o2v0L2 o1v1L2
            D_ov_bb("L,R,i,a") -= t1["bb"]("b,i") * tmps_["30_bb_LLov"]("L,R,j,b") * t1_1["bb"]("a,j");

            // D_vv_bb += +1.00 r1_aa(c,j) t1_1_bb(b,i) l2_1_abab(j,i,c,a)
            //         += -1.00 r1_bb(c,j) t1_1_bb(b,i) l2_1_bbbb(j,i,a,c)
            // flops: o0v2L2 += o1v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_bb("L,R,a,b") += t1_1["bb"]("b,i") * tmps_["30_bb_LLov"]("L,R,i,a");

            // flops: o2v0L2  = o2v1L2 o3v2L2 o2v0L2 o3v2L2 o2v0L2 o3v2L2 o2v0L2 o2v1L2 o2v0L2 o3v2L2 o2v0L2 o2v1L2 o2v0L2 o2v1L2 o2v0L2
            //  mems: o2v0L2  = o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2
            tmps_["31_bb_LLoo"]("L,R,i,j")  = 2.00 * l1["bb"]("L,i,a") * r1["bb"]("R,a,j");
            tmps_["31_bb_LLoo"]("L,R,i,j") += l2_1["bbbb"]("L,l,i,a,c") * r2_1["bbbb"]("R,a,c,l,j");
            tmps_["31_bb_LLoo"]("L,R,i,j") += l2["bbbb"]("L,l,i,a,c") * r2["bbbb"]("R,a,c,l,j");
            tmps_["31_bb_LLoo"]("L,R,i,j") += 2.00 * l2_1["abab"]("L,k,i,b,c") * r2_1["abab"]("R,b,c,k,j");
            tmps_["31_bb_LLoo"]("L,R,i,j") += 2.00 * l1_1["bb"]("L,i,a") * r1_1["bb"]("R,a,j");
            tmps_["31_bb_LLoo"]("L,R,i,j") += 2.00 * l2["abab"]("L,k,i,b,c") * r2["abab"]("R,b,c,k,j");
            tmps_["31_bb_LLoo"]("L,R,i,j") += 2.00 * t1_1["bb"]("a,j") * tmps_["30_bb_LLov"]("L,R,i,a");
            tmps_["31_bb_LLoo"]("L,R,i,j") += 2.00 * t1["bb"]("a,j") * tmps_["29_bb_LLov"]("L,R,i,a");
            tmps_["30_bb_LLov"].~TArrayD();
            tmps_["29_bb_LLov"].~TArrayD();

            // D_oo_bb += -0.50 r2_1_bbbb(a,b,k,i) l2_1_bbbb(k,j,a,b)
            //         += -1.00 r1_bb(a,i) l1_bb(j,a)
            //         += -0.50 r2_bbbb(a,b,k,i) l2_bbbb(k,j,a,b)
            //         += -0.50 r2_1_abab(a,b,k,i) l2_1_abab(k,j,a,b)
            //         += -0.50 r2_1_abab(b,a,k,i) l2_1_abab(k,j,b,a)
            //         += -1.00 r1_1_bb(a,i) l1_1_bb(j,a)
            //         += -0.50 r2_abab(a,b,k,i) l2_abab(k,j,a,b)
            //         += -0.50 r2_abab(b,a,k,i) l2_abab(k,j,b,a)
            //         += -1.00 r1_aa(b,k) t1_1_bb(a,i) l2_1_abab(k,j,b,a)
            //         += +1.00 r1_bb(b,k) t1_1_bb(a,i) l2_1_bbbb(k,j,a,b)
            //         += -1.00 r1_aa(b,k) t1_bb(a,i) l2_abab(k,j,b,a)
            //         += -1.00 r1_1_aa(b,k) t1_bb(a,i) l2_1_abab(k,j,b,a)
            //         += +1.00 r1_1_bb(b,k) t1_bb(a,i) l2_1_bbbb(k,j,a,b)
            //         += +1.00 r1_bb(b,k) t1_bb(a,i) l2_bbbb(k,j,a,b)
            D_oo_bb("L,R,i,j") -= 0.50 * tmps_["31_bb_LLoo"]("L,R,j,i");

            // D_ov_bb += -0.50 r2_1_bbbb(b,c,k,i) t1_bb(a,j) l2_1_bbbb(k,j,b,c)
            //         += -1.00 r1_bb(b,i) t1_bb(a,j) l1_bb(j,b)
            //         += -0.50 r2_bbbb(b,c,k,i) t1_bb(a,j) l2_bbbb(k,j,b,c)
            //         += -0.50 r2_1_abab(b,c,k,i) t1_bb(a,j) l2_1_abab(k,j,b,c)
            //         += -0.50 r2_1_abab(c,b,k,i) t1_bb(a,j) l2_1_abab(k,j,c,b)
            //         += -1.00 r1_1_bb(b,i) t1_bb(a,j) l1_1_bb(j,b)
            //         += -0.50 r2_abab(b,c,k,i) t1_bb(a,j) l2_abab(k,j,b,c)
            //         += -0.50 r2_abab(c,b,k,i) t1_bb(a,j) l2_abab(k,j,c,b)
            //         += -1.00 r1_aa(c,k) t1_bb(a,j) t1_1_bb(b,i) l2_1_abab(k,j,c,b)
            //         += +1.00 r1_bb(c,k) t1_bb(a,j) t1_1_bb(b,i) l2_1_bbbb(k,j,b,c)
            //         += -1.00 r1_aa(c,k) t1_bb(a,j) t1_bb(b,i) l2_abab(k,j,c,b)
            //         += -1.00 r1_1_aa(c,k) t1_bb(a,j) t1_bb(b,i) l2_1_abab(k,j,c,b)
            //         += +1.00 r1_1_bb(c,k) t1_bb(a,j) t1_bb(b,i) l2_1_bbbb(k,j,b,c)
            //         += +1.00 r1_bb(c,k) t1_bb(a,j) t1_bb(b,i) l2_bbbb(k,j,b,c)
            // flops: o1v1L2 += o2v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") -= 0.50 * t1["bb"]("a,j") * tmps_["31_bb_LLoo"]("L,R,j,i");
            tmps_["31_bb_LLoo"].~TArrayD();

            // flops: o1v1L2  = o2v2L2 o2v2L2 o1v1L2
            //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2
            tmps_["32_aa_LLov"]("L,R,i,a")  = -1.00 * l2_1["aaaa"]("L,k,i,a,c") * r1["aa"]("R,c,k");
            tmps_["32_aa_LLov"]("L,R,i,a") += l2_1["abab"]("L,i,j,a,b") * r1["bb"]("R,b,j");

            // D_ov_aa += -1.00 r1_bb(c,k) t2_1_aaaa(a,b,j,i) l2_1_abab(j,k,b,c)
            //         += +1.00 r1_aa(c,k) t2_1_aaaa(a,b,j,i) l2_1_aaaa(k,j,b,c)
            // flops: o1v1L2 += o2v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") -= t2_1["aaaa"]("a,b,j,i") * tmps_["32_aa_LLov"]("L,R,j,b");

            // D_ov_aa += -1.00 r1_bb(c,k) t1_aa(b,i) t1_1_aa(a,j) l2_1_abab(j,k,b,c)
            //         += +1.00 r1_aa(c,k) t1_aa(b,i) t1_1_aa(a,j) l2_1_aaaa(k,j,b,c)
            // flops: o1v1L2 += o2v1L2 o2v1L2
            //  mems: o1v1L2 += o2v0L2 o1v1L2
            D_ov_aa("L,R,i,a") -= t1["aa"]("b,i") * tmps_["32_aa_LLov"]("L,R,j,b") * t1_1["aa"]("a,j");

            // D_ov_bb += +1.00 r1_bb(c,k) t2_1_abab(b,a,j,i) l2_1_abab(j,k,b,c)
            //         += -1.00 r1_aa(c,k) t2_1_abab(b,a,j,i) l2_1_aaaa(k,j,b,c)
            // flops: o1v1L2 += o2v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") += t2_1["abab"]("b,a,j,i") * tmps_["32_aa_LLov"]("L,R,j,b");

            // D_vv_aa += +1.00 r1_bb(c,j) t1_1_aa(b,i) l2_1_abab(i,j,a,c)
            //         += -1.00 r1_aa(c,j) t1_1_aa(b,i) l2_1_aaaa(j,i,a,c)
            // flops: o0v2L2 += o1v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_aa("L,R,a,b") += t1_1["aa"]("b,i") * tmps_["32_aa_LLov"]("L,R,i,a");

            // flops: o1v1L2  = o2v2L2 o2v2L2 o2v2L2 o1v1L2 o1v1L2 o2v2L2 o1v1L2
            //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2
            tmps_["33_aa_LLov"]("L,R,i,a")  = -1.00 * l2["aaaa"]("L,j,i,a,b") * r1["aa"]("R,b,j");
            tmps_["33_aa_LLov"]("L,R,i,a") += l2["abab"]("L,i,k,a,c") * r1["bb"]("R,c,k");
            tmps_["33_aa_LLov"]("L,R,i,a") += l2_1["abab"]("L,i,k,a,c") * r1_1["bb"]("R,c,k");
            tmps_["33_aa_LLov"]("L,R,i,a") -= l2_1["aaaa"]("L,j,i,a,b") * r1_1["aa"]("R,b,j");

            // D_ov_aa += -1.00 r1_1_bb(c,k) t2_aaaa(a,b,j,i) l2_1_abab(j,k,b,c)
            //         += -1.00 r1_bb(c,k) t2_aaaa(a,b,j,i) l2_abab(j,k,b,c)
            //         += +1.00 r1_aa(c,k) t2_aaaa(a,b,j,i) l2_aaaa(k,j,b,c)
            //         += +1.00 r1_1_aa(c,k) t2_aaaa(a,b,j,i) l2_1_aaaa(k,j,b,c)
            // flops: o1v1L2 += o2v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") -= t2["aaaa"]("a,b,j,i") * tmps_["33_aa_LLov"]("L,R,j,b");

            // D_ov_bb += +1.00 r1_1_bb(c,k) t2_abab(b,a,j,i) l2_1_abab(j,k,b,c)
            //         += +1.00 r1_bb(c,k) t2_abab(b,a,j,i) l2_abab(j,k,b,c)
            //         += -1.00 r1_aa(c,k) t2_abab(b,a,j,i) l2_aaaa(k,j,b,c)
            //         += -1.00 r1_1_aa(c,k) t2_abab(b,a,j,i) l2_1_aaaa(k,j,b,c)
            // flops: o1v1L2 += o2v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") += t2["abab"]("b,a,j,i") * tmps_["33_aa_LLov"]("L,R,j,b");

            // D_vo_aa += +1.00 r1_1_bb(b,j) l2_1_abab(i,j,a,b)
            //         += +1.00 r1_bb(b,j) l2_abab(i,j,a,b)
            //         += -1.00 r1_aa(b,j) l2_aaaa(j,i,a,b)
            //         += -1.00 r1_1_aa(b,j) l2_1_aaaa(j,i,a,b)
            D_vo_aa("L,R,a,i") += tmps_["33_aa_LLov"]("L,R,i,a");

            // D_vv_aa += +1.00 r1_1_bb(c,j) t1_aa(b,i) l2_1_abab(i,j,a,c)
            //         += +1.00 r1_bb(c,j) t1_aa(b,i) l2_abab(i,j,a,c)
            //         += -1.00 r1_aa(c,j) t1_aa(b,i) l2_aaaa(j,i,a,c)
            //         += -1.00 r1_1_aa(c,j) t1_aa(b,i) l2_1_aaaa(j,i,a,c)
            // flops: o0v2L2 += o1v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_aa("L,R,a,b") += t1["aa"]("b,i") * tmps_["33_aa_LLov"]("L,R,i,a");

            // flops: o2v0L2  = o3v2L2 o2v1L2 o3v2L2 o3v2L2 o2v0L2 o2v0L2 o2v0L2 o2v1L2 o2v0L2 o3v2L2 o2v0L2 o2v1L2 o2v1L2 o2v0L2 o2v0L2
            //  mems: o2v0L2  = o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2
            tmps_["34_aa_LLoo"]("L,R,i,j")  = 0.50 * l2["aaaa"]("L,k,j,a,b") * r2["aaaa"]("R,a,b,k,i");
            tmps_["34_aa_LLoo"]("L,R,i,j") += l1_1["aa"]("L,j,a") * r1_1["aa"]("R,a,i");
            tmps_["34_aa_LLoo"]("L,R,i,j") += l2["abab"]("L,j,l,a,c") * r2["abab"]("R,a,c,i,l");
            tmps_["34_aa_LLoo"]("L,R,i,j") += l2_1["abab"]("L,j,l,a,c") * r2_1["abab"]("R,a,c,i,l");
            tmps_["34_aa_LLoo"]("L,R,i,j") += l1["aa"]("L,j,a") * r1["aa"]("R,a,i");
            tmps_["34_aa_LLoo"]("L,R,i,j") += 0.50 * l2_1["aaaa"]("L,k,j,a,b") * r2_1["aaaa"]("R,a,b,k,i");
            tmps_["34_aa_LLoo"]("L,R,i,j") += t1["aa"]("a,i") * tmps_["33_aa_LLov"]("L,R,j,a");
            tmps_["34_aa_LLoo"]("L,R,i,j") += t1_1["aa"]("a,i") * tmps_["32_aa_LLov"]("L,R,j,a");
            tmps_["33_aa_LLov"].~TArrayD();
            tmps_["32_aa_LLov"].~TArrayD();

            // D_oo_aa += -1.00 r1_bb(b,k) t1_1_aa(a,i) l2_1_abab(j,k,a,b)
            //         += +1.00 r1_aa(b,k) t1_1_aa(a,i) l2_1_aaaa(k,j,a,b)
            //         += -1.00 r1_1_bb(b,k) t1_aa(a,i) l2_1_abab(j,k,a,b)
            //         += -1.00 r1_bb(b,k) t1_aa(a,i) l2_abab(j,k,a,b)
            //         += +1.00 r1_aa(b,k) t1_aa(a,i) l2_aaaa(k,j,a,b)
            //         += +1.00 r1_1_aa(b,k) t1_aa(a,i) l2_1_aaaa(k,j,a,b)
            //         += -0.50 r2_1_abab(a,b,i,k) l2_1_abab(j,k,a,b)
            //         += -0.50 r2_1_abab(b,a,i,k) l2_1_abab(j,k,b,a)
            //         += -0.50 r2_abab(a,b,i,k) l2_abab(j,k,a,b)
            //         += -0.50 r2_abab(b,a,i,k) l2_abab(j,k,b,a)
            //         += -1.00 r1_1_aa(a,i) l1_1_aa(j,a)
            //         += -0.50 r2_aaaa(a,b,k,i) l2_aaaa(k,j,a,b)
            //         += -1.00 r1_aa(a,i) l1_aa(j,a)
            //         += -0.50 r2_1_aaaa(a,b,k,i) l2_1_aaaa(k,j,a,b)
            D_oo_aa("L,R,i,j") -= tmps_["34_aa_LLoo"]("L,R,i,j");

            // D_ov_aa += -1.00 r1_bb(c,k) t1_aa(a,j) t1_1_aa(b,i) l2_1_abab(j,k,b,c)
            //         += +1.00 r1_aa(c,k) t1_aa(a,j) t1_1_aa(b,i) l2_1_aaaa(k,j,b,c)
            //         += -1.00 r1_1_bb(c,k) t1_aa(a,j) t1_aa(b,i) l2_1_abab(j,k,b,c)
            //         += -1.00 r1_bb(c,k) t1_aa(a,j) t1_aa(b,i) l2_abab(j,k,b,c)
            //         += +1.00 r1_aa(c,k) t1_aa(a,j) t1_aa(b,i) l2_aaaa(k,j,b,c)
            //         += +1.00 r1_1_aa(c,k) t1_aa(a,j) t1_aa(b,i) l2_1_aaaa(k,j,b,c)
            //         += -0.50 r2_1_abab(b,c,i,k) t1_aa(a,j) l2_1_abab(j,k,b,c)
            //         += -0.50 r2_1_abab(c,b,i,k) t1_aa(a,j) l2_1_abab(j,k,c,b)
            //         += -0.50 r2_abab(b,c,i,k) t1_aa(a,j) l2_abab(j,k,b,c)
            //         += -0.50 r2_abab(c,b,i,k) t1_aa(a,j) l2_abab(j,k,c,b)
            //         += -1.00 r1_1_aa(b,i) t1_aa(a,j) l1_1_aa(j,b)
            //         += -0.50 r2_aaaa(b,c,k,i) t1_aa(a,j) l2_aaaa(k,j,b,c)
            //         += -1.00 r1_aa(b,i) t1_aa(a,j) l1_aa(j,b)
            //         += -0.50 r2_1_aaaa(b,c,k,i) t1_aa(a,j) l2_1_aaaa(k,j,b,c)
            // flops: o1v1L2 += o2v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") -= t1["aa"]("a,j") * tmps_["34_aa_LLoo"]("L,R,i,j");
            tmps_["34_aa_LLoo"].~TArrayD();
        }

        rdm_timer.stop();
        Printf("Done --> Time: %s\n", rdm_timer.elapsed().c_str());

    }

#ifndef BUILD_QED_2RDM
    void EOM_EE_QED_RDM_21::compute_eom_2rdm() {
        throw PsiException("EOM_EE_QED_RDM_21::compute_eom_2rdm() should not be called when BUILD_QED_2RDM is not defined.",
                           __FILE__, __LINE__);
    }
#else
    void EOM_EE_QED_RDM_21::compute_eom_2rdm() {
        // reinitialize the 2-RDMs if not already done (only allowed for a single state or transition state)
        vector<string> dims = {"oa", "va", "ob", "vb"};
        vector<size_t> dim_vec = {oa_, va_, ob_, vb_};
        for (int i = 0; i < dims.size(); i++) {
            for (int j = 0; j < dims.size(); j++) {
                string dim1 = dims[i], dim2 = dims[j];

                string ov = dim1.substr(0, 1) + dim2.substr(0, 1);
                string spin = dim1.substr(1, 1) + dim2.substr(1, 1);
                for (int k = 0; k < dims.size(); k++) {
                    for (int l = 0; l < dims.size(); l++) {
                        string dim3 = dims[k], dim4 = dims[l];

                        ov += dim3.substr(0, 1) + dim4.substr(0, 1);
                        spin += dim3.substr(1, 1) + dim4.substr(1, 1);
                        string rdm_str = "D2_" + spin + "_" + ov;

                        // check that spin is either aaaa, abab, or bbbb.
                        if (spin[0] != spin[2] || spin[1] != spin[3])
                            continue; // if not, continue

                        RDM_blks_[rdm_str] = TArrayD(world_, makeRange({dim_vec[i], dim_vec[j], dim_vec[k], dim_vec[l]}));
                        RDM_blks_[rdm_str].fill(0.0);
                    }
                }
            }
        }
        world_.gop.fence();

        Timer rdm_timer; rdm_timer.start();

        rdm2_21_1();
        rdm2_21_2();
        rdm2_21_3();
        rdm2_21_4();
        rdm2_21_5();
        rdm2_21_6();
        world.gop.fence();

        rdm_timer.stop();
        Printf("Done --> Time: %s\n", rdm_timer.elapsed().c_str());
    }
#endif

} // cc_cavity
