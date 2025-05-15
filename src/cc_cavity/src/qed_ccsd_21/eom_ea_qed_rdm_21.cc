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

#include "cc_cavity/include/qed_ccsd_21/eom_ea_qed_rdm_21.h"
#include <psi4/libmints/writer.h>
#include "misc/nonsym_davidson_solver.h"

namespace hilbert {
    EOM_EA_QED_RDM_21::EOM_EA_QED_RDM_21(const shared_ptr<EOM_Driver>& eom_driver, Options & options) : EOM_EA_RDM(eom_driver, options) {
    }

    void EOM_EA_QED_RDM_21::compute_eom_1rdm() {

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



        // extract amplitudes
        std::map<std::string, TA::TArrayD> t1 {
                {"aa", eom_driver_->cc_wfn_->amplitudes_["t1_aa"]},
                {"bb", eom_driver_->cc_wfn_->amplitudes_["t1_bb"]}
        };
        std::map<std::string, TA::TArrayD> t2 {
                {"aaaa", eom_driver_->cc_wfn_->amplitudes_["t2_aaaa"]},
                {"abab", eom_driver_->cc_wfn_->amplitudes_["t2_abab"]},
                {"bbbb", eom_driver_->cc_wfn_->amplitudes_["t2_bbbb"]}
        };

        double t0_1;
        foreach_inplace(eom_driver_->cc_wfn_->amplitudes_["t0_1"], [&t0_1](auto &tile){
            for(auto &x : tile.range())
                t0_1 = tile[x];
        });

        std::map<std::string, TA::TArrayD> t1_1 {
                {"aa", eom_driver_->cc_wfn_->amplitudes_["t1_1_aa"]},
                {"bb", eom_driver_->cc_wfn_->amplitudes_["t1_1_bb"]}
        };
        std::map<std::string, TA::TArrayD> t2_1 {
                {"aaaa", eom_driver_->cc_wfn_->amplitudes_["t2_1_aaaa"]},
                {"abab", eom_driver_->cc_wfn_->amplitudes_["t2_1_abab"]},
                {"bbbb", eom_driver_->cc_wfn_->amplitudes_["t2_1_bbbb"]}
        };


        std::map<std::string, TA::TArrayD> r1 {
                {"a", eom_driver_->evec_blks_["r1_a"]},
        };
        std::map<std::string, TA::TArrayD> l1 {
                {"a", eom_driver_->evec_blks_["l1_a"]},
        };
        std::map<std::string, TA::TArrayD> r1_1 {
                {"a", eom_driver_->evec_blks_["r1_1_a"]},
        };
        std::map<std::string, TA::TArrayD> l1_1 {
                {"a", eom_driver_->evec_blks_["l1_1_a"]},
        };


        // extract 2-body amplitudes
        std::map<std::string, TA::TArrayD> r2 {
                {"aaa", eom_driver_->evec_blks_["r2_aaa"]},
                {"abb", eom_driver_->evec_blks_["r2_abb"]},
        };
        std::map<std::string, TA::TArrayD> l2 {
                {"aaa", eom_driver_->evec_blks_["l2_aaa"]},
                {"bab", eom_driver_->evec_blks_["l2_bab"]},
        };
        std::map<std::string, TA::TArrayD> r2_1 {
                {"aaa", eom_driver_->evec_blks_["r2_1_aaa"]},
                {"abb", eom_driver_->evec_blks_["r2_1_abb"]},
        };
        std::map<std::string, TA::TArrayD> l2_1 {
                {"aaa", eom_driver_->evec_blks_["l2_1_aaa"]},
                {"bab", eom_driver_->evec_blks_["l2_1_bab"]},
        };

        // assuming singlet reference, build beta blocks from alpha
        r1["b"] = r1["a"].clone();
        r2["bbb"] = r2["aaa"].clone();
        r2["aba"] = r2["abb"].clone();

        l1["b"] = l1["a"].clone();
        l2["bbb"] = l2["aaa"].clone();
        l2["aab"] = l2["bab"].clone();

        r1_1["b"] = r1_1["a"].clone();
        r2_1["bbb"] = r2_1["aaa"].clone();
        r2_1["aba"] = r2_1["abb"].clone();

        l1_1["b"] = l1_1["a"].clone();
        l2_1["bbb"] = l2_1["aaa"].clone();
        l2_1["aab"] = l2_1["bab"].clone();

        // get integrals
        TArrayMap &Id  = eom_driver_->cc_wfn_->Id_blks_;

        rdm_timer.start();

        {
            // D_ov_bb  = -1.00 r2_1_bbb(a,b,i) l1_1_b(b)
            // flops: o1v1L2  = o1v2L2
            //  mems: o1v1L2  = o1v1L2
            D_ov_bb("L,R,i,a")  = -1.00 * l1_1["b"]("L,b") * r2_1["bbb"]("R,a,b,i");

            // D_vv_bb  = +1.00 r1_b(b) l1_b(a)
            // flops: o0v2L2  = o0v2L2
            //  mems: o0v2L2  = o0v2L2
            D_vv_bb("L,R,a,b")  = l1["b"]("L,a") * r1["b"]("R,b");

            // D_vv_aa  = +1.00 r1_a(b) l1_a(a)
            // flops: o0v2L2  = o0v2L2
            //  mems: o0v2L2  = o0v2L2
            D_vv_aa("L,R,a,b")  = l1["a"]("L,a") * r1["a"]("R,b");

            // D_ov_aa  = -1.00 r2_1_aba(a,b,i) l1_1_b(b)
            // flops: o1v1L2  = o1v2L2
            //  mems: o1v1L2  = o1v1L2
            D_ov_aa("L,R,i,a")  = -1.00 * l1_1["b"]("L,b") * r2_1["aba"]("R,a,b,i");

            // D_ov_aa += -1.00 r2_aba(a,b,i) l1_b(b)
            // flops: o1v1L2 += o1v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") -= l1["b"]("L,b") * r2["aba"]("R,a,b,i");

            // D_ov_aa += -1.00 r2_aaa(a,b,i) l1_a(b)
            // flops: o1v1L2 += o1v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") -= l1["a"]("L,b") * r2["aaa"]("R,a,b,i");

            // D_ov_aa += -1.00 r2_1_aaa(a,b,i) l1_1_a(b)
            // flops: o1v1L2 += o1v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") -= l1_1["a"]("L,b") * r2_1["aaa"]("R,a,b,i");

            // D_ov_aa += -1.00 r1_a(a) t1_aa(b,i) l1_a(b)
            // flops: o1v1L2 += o1v1L1 o1v1L2
            //  mems: o1v1L2 += o1v0L1 o1v1L2
            D_ov_aa("L,R,i,a") -= l1["a"]("L,b") * t1["aa"]("b,i") * r1["a"]("R,a");

            // D_ov_aa += -0.50 r1_a(a) t2_abab(c,b,i,j) l2_bab(j,c,b)
            //         += -0.50 r1_a(a) t2_abab(b,c,i,j) l2_bab(j,b,c)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v0L1 o1v1L2
            D_ov_aa("L,R,i,a") -= l2["bab"]("L,j,c,b") * t2["abab"]("c,b,i,j") * r1["a"]("R,a");

            // D_ov_aa += -0.50 r1_a(a) t2_1_abab(c,b,i,j) l2_1_bab(j,c,b)
            //         += -0.50 r1_a(a) t2_1_abab(b,c,i,j) l2_1_bab(j,b,c)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v0L1 o1v1L2
            D_ov_aa("L,R,i,a") -= l2_1["bab"]("L,j,c,b") * t2_1["abab"]("c,b,i,j") * r1["a"]("R,a");

            // D_ov_aa += -0.50 r1_1_a(a) t2_abab(c,b,i,j) l2_1_bab(j,c,b)
            //         += -0.50 r1_1_a(a) t2_abab(b,c,i,j) l2_1_bab(j,b,c)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v0L1 o1v1L2
            D_ov_aa("L,R,i,a") -= l2_1["bab"]("L,j,c,b") * t2["abab"]("c,b,i,j") * r1_1["a"]("R,a");

            // D_ov_aa += +0.50 r1_a(a) t2_aaaa(c,b,j,i) l2_aaa(j,c,b)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v0L1 o1v1L2
            D_ov_aa("L,R,i,a") += 0.50 * l2["aaa"]("L,j,c,b") * t2["aaaa"]("c,b,j,i") * r1["a"]("R,a");

            // D_ov_aa += +0.50 r1_1_a(a) t2_aaaa(c,b,j,i) l2_1_aaa(j,c,b)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v0L1 o1v1L2
            D_ov_aa("L,R,i,a") += 0.50 * l2_1["aaa"]("L,j,c,b") * t2["aaaa"]("c,b,j,i") * r1_1["a"]("R,a");

            // D_ov_aa += +0.50 r1_a(a) t2_1_aaaa(c,b,j,i) l2_1_aaa(j,c,b)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v0L1 o1v1L2
            D_ov_aa("L,R,i,a") += 0.50 * l2_1["aaa"]("L,j,c,b") * t2_1["aaaa"]("c,b,j,i") * r1["a"]("R,a");

            // D_ov_aa += -1.00 r1_b(c) t2_1_abab(a,b,i,j) l2_1_bbb(j,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_aa("L,R,i,a") -= l2_1["bbb"]("L,j,b,c") * t2_1["abab"]("a,b,i,j") * r1["b"]("R,c");

            // D_ov_aa += -1.00 r1_b(c) t2_abab(a,b,i,j) l2_bbb(j,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_aa("L,R,i,a") -= l2["bbb"]("L,j,b,c") * t2["abab"]("a,b,i,j") * r1["b"]("R,c");

            // D_ov_aa += -1.00 r1_1_b(c) t2_abab(a,b,i,j) l2_1_bbb(j,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_aa("L,R,i,a") -= l2_1["bbb"]("L,j,b,c") * t2["abab"]("a,b,i,j") * r1_1["b"]("R,c");

            // D_ov_aa += +1.00 r1_b(c) t2_aaaa(a,b,j,i) l2_aab(j,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_aa("L,R,i,a") += l2["aab"]("L,j,b,c") * t2["aaaa"]("a,b,j,i") * r1["b"]("R,c");

            // D_ov_aa += +1.00 r1_1_b(c) t2_aaaa(a,b,j,i) l2_1_aab(j,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_aa("L,R,i,a") += l2_1["aab"]("L,j,b,c") * t2["aaaa"]("a,b,j,i") * r1_1["b"]("R,c");

            // D_ov_aa += +1.00 r1_b(c) t2_1_aaaa(a,b,j,i) l2_1_aab(j,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_aa("L,R,i,a") += l2_1["aab"]("L,j,b,c") * t2_1["aaaa"]("a,b,j,i") * r1["b"]("R,c");

            // D_ov_aa += +1.00 r1_a(c) t2_1_abab(a,b,i,j) l2_1_bab(j,c,b)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_aa("L,R,i,a") += l2_1["bab"]("L,j,c,b") * t2_1["abab"]("a,b,i,j") * r1["a"]("R,c");

            // D_ov_aa += +1.00 r1_a(c) t2_abab(a,b,i,j) l2_bab(j,c,b)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_aa("L,R,i,a") += l2["bab"]("L,j,c,b") * t2["abab"]("a,b,i,j") * r1["a"]("R,c");

            // D_ov_aa += +1.00 r1_1_a(c) t2_abab(a,b,i,j) l2_1_bab(j,c,b)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_aa("L,R,i,a") += l2_1["bab"]("L,j,c,b") * t2["abab"]("a,b,i,j") * r1_1["a"]("R,c");

            // D_ov_aa += +1.00 r1_a(c) t2_1_aaaa(a,b,j,i) l2_1_aaa(j,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_aa("L,R,i,a") += l2_1["aaa"]("L,j,b,c") * t2_1["aaaa"]("a,b,j,i") * r1["a"]("R,c");

            // D_ov_aa += +1.00 r1_1_a(c) t2_aaaa(a,b,j,i) l2_1_aaa(j,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_aa("L,R,i,a") += l2_1["aaa"]("L,j,b,c") * t2["aaaa"]("a,b,j,i") * r1_1["a"]("R,c");

            // D_ov_aa += +1.00 r1_a(c) t2_aaaa(a,b,j,i) l2_aaa(j,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_aa("L,R,i,a") += l2["aaa"]("L,j,b,c") * t2["aaaa"]("a,b,j,i") * r1["a"]("R,c");

            // D_ov_aa += -1.00 r2_aba(a,c,j) t1_1_aa(b,i) l2_1_aab(j,b,c)
            // flops: o1v1L2 += o2v2L1 o2v2L2
            //  mems: o1v1L2 += o2v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= l2_1["aab"]("L,j,b,c") * t1_1["aa"]("b,i") * r2["aba"]("R,a,c,j");

            // D_ov_aa += -0.50 r2_aba(b,c,i) t1_1_aa(a,j) l2_1_aab(j,b,c)
            //         += -0.50 r2_aba(c,b,i) t1_1_aa(a,j) l2_1_aab(j,c,b)
            // flops: o1v1L2 += o2v2L2 o2v1L2
            //  mems: o1v1L2 += o2v0L2 o1v1L2
            D_ov_aa("L,R,i,a") -= l2_1["aab"]("L,j,b,c") * r2["aba"]("R,b,c,i") * t1_1["aa"]("a,j");

            // D_ov_aa += -1.00 r2_abb(a,c,j) t1_1_aa(b,i) l2_1_bab(j,b,c)
            // flops: o1v1L2 += o2v2L1 o2v2L2
            //  mems: o1v1L2 += o2v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= l2_1["bab"]("L,j,b,c") * t1_1["aa"]("b,i") * r2["abb"]("R,a,c,j");

            // D_ov_aa += -1.00 r2_aaa(a,c,j) t1_1_aa(b,i) l2_1_aaa(j,b,c)
            // flops: o1v1L2 += o2v2L1 o2v2L2
            //  mems: o1v1L2 += o2v1L1 o1v1L2
            D_ov_aa("L,R,i,a") -= l2_1["aaa"]("L,j,b,c") * t1_1["aa"]("b,i") * r2["aaa"]("R,a,c,j");

            // D_ov_aa += -0.50 r2_aaa(b,c,i) t1_1_aa(a,j) l2_1_aaa(j,b,c)
            // flops: o1v1L2 += o2v2L2 o2v1L2
            //  mems: o1v1L2 += o2v0L2 o1v1L2
            D_ov_aa("L,R,i,a") -= 0.50 * l2_1["aaa"]("L,j,b,c") * r2["aaa"]("R,b,c,i") * t1_1["aa"]("a,j");

            // D_ov_aa += +1.00 r1_b(c) t1_aa(b,i) t1_1_aa(a,j) l2_1_aab(j,b,c)
            // flops: o1v1L2 += o2v2L1 o2v1L2 o2v1L2
            //  mems: o1v1L2 += o2v1L1 o2v0L2 o1v1L2
            D_ov_aa("L,R,i,a") += l2_1["aab"]("L,j,b,c") * t1["aa"]("b,i") * r1["b"]("R,c") * t1_1["aa"]("a,j");

            // D_ov_bb += -1.00 r2_bbb(a,b,i) l1_b(b)
            // flops: o1v1L2 += o1v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") -= l1["b"]("L,b") * r2["bbb"]("R,a,b,i");

            // D_ov_bb += +1.00 r2_1_abb(b,a,i) l1_1_a(b)
            // flops: o1v1L2 += o1v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") += l1_1["a"]("L,b") * r2_1["abb"]("R,b,a,i");

            // D_ov_bb += +1.00 r2_abb(b,a,i) l1_a(b)
            // flops: o1v1L2 += o1v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") += l1["a"]("L,b") * r2["abb"]("R,b,a,i");

            // D_ov_bb += -1.00 r1_b(a) t1_1_bb(b,i) l1_1_b(b)
            // flops: o1v1L2 += o1v1L1 o1v1L2
            //  mems: o1v1L2 += o1v0L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l1_1["b"]("L,b") * t1_1["bb"]("b,i") * r1["b"]("R,a");

            // D_ov_bb += -1.00 r1_1_b(a) t1_bb(b,i) l1_1_b(b)
            // flops: o1v1L2 += o1v1L1 o1v1L2
            //  mems: o1v1L2 += o1v0L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l1_1["b"]("L,b") * t1["bb"]("b,i") * r1_1["b"]("R,a");

            // D_ov_bb += -1.00 r1_b(a) t1_bb(b,i) l1_b(b)
            // flops: o1v1L2 += o1v1L1 o1v1L2
            //  mems: o1v1L2 += o1v0L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l1["b"]("L,b") * t1["bb"]("b,i") * r1["b"]("R,a");

            // D_ov_bb += +0.50 r1_b(a) t2_1_bbbb(c,b,j,i) l2_1_bbb(j,c,b)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v0L1 o1v1L2
            D_ov_bb("L,R,i,a") += 0.50 * l2_1["bbb"]("L,j,c,b") * t2_1["bbbb"]("c,b,j,i") * r1["b"]("R,a");

            // D_ov_bb += +0.50 r1_1_b(a) t2_bbbb(c,b,j,i) l2_1_bbb(j,c,b)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v0L1 o1v1L2
            D_ov_bb("L,R,i,a") += 0.50 * l2_1["bbb"]("L,j,c,b") * t2["bbbb"]("c,b,j,i") * r1_1["b"]("R,a");

            // D_ov_bb += +0.50 r1_b(a) t2_bbbb(c,b,j,i) l2_bbb(j,c,b)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v0L1 o1v1L2
            D_ov_bb("L,R,i,a") += 0.50 * l2["bbb"]("L,j,c,b") * t2["bbbb"]("c,b,j,i") * r1["b"]("R,a");

            // D_ov_bb += +0.50 r1_b(a) t2_1_abab(c,b,j,i) l2_1_aab(j,c,b)
            //         += +0.50 r1_b(a) t2_1_abab(b,c,j,i) l2_1_aab(j,b,c)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v0L1 o1v1L2
            D_ov_bb("L,R,i,a") += l2_1["aab"]("L,j,c,b") * t2_1["abab"]("c,b,j,i") * r1["b"]("R,a");

            // D_ov_bb += +0.50 r1_b(a) t2_abab(c,b,j,i) l2_aab(j,c,b)
            //         += +0.50 r1_b(a) t2_abab(b,c,j,i) l2_aab(j,b,c)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v0L1 o1v1L2
            D_ov_bb("L,R,i,a") += l2["aab"]("L,j,c,b") * t2["abab"]("c,b,j,i") * r1["b"]("R,a");

            // D_ov_bb += +0.50 r1_1_b(a) t2_abab(c,b,j,i) l2_1_aab(j,c,b)
            //         += +0.50 r1_1_b(a) t2_abab(b,c,j,i) l2_1_aab(j,b,c)
            // flops: o1v1L2 += o2v2L1 o1v1L2
            //  mems: o1v1L2 += o1v0L1 o1v1L2
            D_ov_bb("L,R,i,a") += l2_1["aab"]("L,j,c,b") * t2["abab"]("c,b,j,i") * r1_1["b"]("R,a");

            // D_ov_bb += +1.00 r1_b(c) t2_bbbb(a,b,j,i) l2_bbb(j,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_bb("L,R,i,a") += l2["bbb"]("L,j,b,c") * t2["bbbb"]("a,b,j,i") * r1["b"]("R,c");

            // D_ov_bb += +1.00 r1_1_b(c) t2_bbbb(a,b,j,i) l2_1_bbb(j,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_bb("L,R,i,a") += l2_1["bbb"]("L,j,b,c") * t2["bbbb"]("a,b,j,i") * r1_1["b"]("R,c");

            // D_ov_bb += +1.00 r1_b(c) t2_1_bbbb(a,b,j,i) l2_1_bbb(j,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_bb("L,R,i,a") += l2_1["bbb"]("L,j,b,c") * t2_1["bbbb"]("a,b,j,i") * r1["b"]("R,c");

            // D_ov_bb += -1.00 r1_b(c) t2_abab(b,a,j,i) l2_aab(j,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l2["aab"]("L,j,b,c") * t2["abab"]("b,a,j,i") * r1["b"]("R,c");

            // D_ov_bb += -1.00 r1_1_b(c) t2_abab(b,a,j,i) l2_1_aab(j,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l2_1["aab"]("L,j,b,c") * t2["abab"]("b,a,j,i") * r1_1["b"]("R,c");

            // D_ov_bb += -1.00 r1_b(c) t2_1_abab(b,a,j,i) l2_1_aab(j,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l2_1["aab"]("L,j,b,c") * t2_1["abab"]("b,a,j,i") * r1["b"]("R,c");

            // D_ov_bb += -1.00 r1_a(c) t2_bbbb(a,b,j,i) l2_bab(j,c,b)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l2["bab"]("L,j,c,b") * t2["bbbb"]("a,b,j,i") * r1["a"]("R,c");

            // D_ov_bb += -1.00 r1_1_a(c) t2_bbbb(a,b,j,i) l2_1_bab(j,c,b)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l2_1["bab"]("L,j,c,b") * t2["bbbb"]("a,b,j,i") * r1_1["a"]("R,c");

            // D_ov_bb += -1.00 r1_a(c) t2_1_bbbb(a,b,j,i) l2_1_bab(j,c,b)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l2_1["bab"]("L,j,c,b") * t2_1["bbbb"]("a,b,j,i") * r1["a"]("R,c");

            // D_ov_bb += -1.00 r1_a(c) t2_abab(b,a,j,i) l2_aaa(j,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l2["aaa"]("L,j,b,c") * t2["abab"]("b,a,j,i") * r1["a"]("R,c");

            // D_ov_bb += -1.00 r1_a(c) t2_1_abab(b,a,j,i) l2_1_aaa(j,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l2_1["aaa"]("L,j,b,c") * t2_1["abab"]("b,a,j,i") * r1["a"]("R,c");

            // D_ov_bb += -1.00 r1_1_a(c) t2_abab(b,a,j,i) l2_1_aaa(j,b,c)
            // flops: o1v1L2 += o2v3L1 o1v2L2
            //  mems: o1v1L2 += o1v2L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l2_1["aaa"]("L,j,b,c") * t2["abab"]("b,a,j,i") * r1_1["a"]("R,c");

            // D_ov_bb += -1.00 r2_bbb(a,c,j) t1_1_bb(b,i) l2_1_bbb(j,b,c)
            // flops: o1v1L2 += o2v2L1 o2v2L2
            //  mems: o1v1L2 += o2v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l2_1["bbb"]("L,j,b,c") * t1_1["bb"]("b,i") * r2["bbb"]("R,a,c,j");

            // D_ov_bb += -0.50 r2_bbb(b,c,i) t1_1_bb(a,j) l2_1_bbb(j,b,c)
            // flops: o1v1L2 += o2v2L2 o2v1L2
            //  mems: o1v1L2 += o2v0L2 o1v1L2
            D_ov_bb("L,R,i,a") -= 0.50 * l2_1["bbb"]("L,j,b,c") * r2["bbb"]("R,b,c,i") * t1_1["bb"]("a,j");

            // D_ov_bb += -1.00 r2_aba(c,a,j) t1_1_bb(b,i) l2_1_aab(j,c,b)
            // flops: o1v1L2 += o2v2L1 o2v2L2
            //  mems: o1v1L2 += o2v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l2_1["aab"]("L,j,c,b") * t1_1["bb"]("b,i") * r2["aba"]("R,c,a,j");

            // D_ov_bb += -1.00 r2_abb(c,a,j) t1_1_bb(b,i) l2_1_bab(j,c,b)
            // flops: o1v1L2 += o2v2L1 o2v2L2
            //  mems: o1v1L2 += o2v1L1 o1v1L2
            D_ov_bb("L,R,i,a") -= l2_1["bab"]("L,j,c,b") * t1_1["bb"]("b,i") * r2["abb"]("R,c,a,j");

            // D_ov_bb += -0.50 r2_abb(b,c,i) t1_1_bb(a,j) l2_1_bab(j,b,c)
            //         += -0.50 r2_abb(c,b,i) t1_1_bb(a,j) l2_1_bab(j,c,b)
            // flops: o1v1L2 += o2v2L2 o2v1L2
            //  mems: o1v1L2 += o2v0L2 o1v1L2
            D_ov_bb("L,R,i,a") -= l2_1["bab"]("L,j,b,c") * r2["abb"]("R,b,c,i") * t1_1["bb"]("a,j");

            // D_vv_aa += +1.00 r1_1_a(b) l1_1_a(a)
            // flops: o0v2L2 += o0v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_aa("L,R,a,b") += l1_1["a"]("L,a") * r1_1["a"]("R,b");

            // D_vv_aa += -1.00 r1_b(c) t1_1_aa(b,i) l2_1_aab(i,a,c)
            // flops: o0v2L2 += o1v2L2 o1v2L2
            //  mems: o0v2L2 += o1v1L2 o0v2L2
            D_vv_aa("L,R,a,b") -= l2_1["aab"]("L,i,a,c") * r1["b"]("R,c") * t1_1["aa"]("b,i");

            // D_vv_bb += +1.00 r1_1_b(b) l1_1_b(a)
            // flops: o0v2L2 += o0v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_bb("L,R,a,b") += l1_1["b"]("L,a") * r1_1["b"]("R,b");

            // flops: o1v1L2  = o1v2L2
            //  mems: o1v1L2  = o1v1L2
            tmps_["1_aa_LLov"]("L,R,i,a")  = l2_1["aaa"]("L,i,a,b") * r1["a"]("R,b");

            // D_ov_aa += +1.00 r1_a(c) t1_aa(b,i) t1_1_aa(a,j) l2_1_aaa(j,b,c)
            // flops: o1v1L2 += o2v1L2 o2v1L2
            //  mems: o1v1L2 += o2v0L2 o1v1L2
            D_ov_aa("L,R,i,a") += t1["aa"]("b,i") * tmps_["1_aa_LLov"]("L,R,j,b") * t1_1["aa"]("a,j");

            // D_vv_aa += -1.00 r1_a(c) t1_1_aa(b,i) l2_1_aaa(i,a,c)
            // flops: o0v2L2 += o1v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_aa("L,R,a,b") -= t1_1["aa"]("b,i") * tmps_["1_aa_LLov"]("L,R,i,a");
            tmps_["1_aa_LLov"].~TArrayD();

            // flops: o1v1L2  = o1v2L2
            //  mems: o1v1L2  = o1v1L2
            tmps_["2_bb_LLov"]("L,R,i,a")  = l2_1["bbb"]("L,i,a,b") * r1["b"]("R,b");

            // D_ov_bb += +1.00 r1_b(c) t1_bb(b,i) t1_1_bb(a,j) l2_1_bbb(j,b,c)
            // flops: o1v1L2 += o2v1L2 o2v1L2
            //  mems: o1v1L2 += o2v0L2 o1v1L2
            D_ov_bb("L,R,i,a") += t1["bb"]("b,i") * tmps_["2_bb_LLov"]("L,R,j,b") * t1_1["bb"]("a,j");

            // D_vv_bb += -1.00 r1_b(c) t1_1_bb(b,i) l2_1_bbb(i,a,c)
            // flops: o0v2L2 += o1v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_bb("L,R,a,b") -= t1_1["bb"]("b,i") * tmps_["2_bb_LLov"]("L,R,i,a");
            tmps_["2_bb_LLov"].~TArrayD();

            // flops: o1v0L1  = o1v1L1
            //  mems: o1v0L1  = o1v0L1
            tmps_["3_a_Lo"]("L,i")  = l1_1["a"]("L,a") * t1_1["aa"]("a,i");

            // D_ov_aa += -1.00 r1_a(a) t1_1_aa(b,i) l1_1_a(b)
            // flops: o1v1L2 += o1v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") -= r1["a"]("R,a") * tmps_["3_a_Lo"]("L,i");
            tmps_["3_a_Lo"].~TArrayD();

            // flops: o1v0L1  = o1v1L1
            //  mems: o1v0L1  = o1v0L1
            tmps_["4_a_Lo"]("L,i")  = l1_1["a"]("L,a") * t1["aa"]("a,i");

            // D_ov_aa += -1.00 r1_1_a(a) t1_aa(b,i) l1_1_a(b)
            // flops: o1v1L2 += o1v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") -= r1_1["a"]("R,a") * tmps_["4_a_Lo"]("L,i");
            tmps_["4_a_Lo"].~TArrayD();

            // flops: o0v0L2  = o0v1L2 o1v2L2 o0v0L2 o1v2L2 o0v0L2 o1v2L2 o0v0L2 o1v2L2 o0v0L2 o0v1L2 o0v0L2 o0v1L2 o0v0L2 o1v2L2 o0v0L2 o1v2L2 o0v0L2 o1v2L2 o0v0L2 o0v1L2 o0v0L2 o1v2L2 o0v0L2
            //  mems: o0v0L2  = o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2
            tmps_["5_LL"]("L,R")  = 2.00 * l1["a"]("L,a") * r1["a"]("R,a");
            tmps_["5_LL"]("L,R") += l2_1["aaa"]("L,i,a,c") * r2_1["aaa"]("R,a,c,i");
            tmps_["5_LL"]("L,R") += 2.00 * l2_1["aab"]("L,i,a,b") * r2_1["aba"]("R,a,b,i");
            tmps_["5_LL"]("L,R") += 2.00 * l2["bab"]("L,j,a,b") * r2["abb"]("R,a,b,j");
            tmps_["5_LL"]("L,R") += l2["bbb"]("L,j,d,b") * r2["bbb"]("R,d,b,j");
            tmps_["5_LL"]("L,R") += 2.00 * l1["b"]("L,d") * r1["b"]("R,d");
            tmps_["5_LL"]("L,R") += 2.00 * l1_1["a"]("L,a") * r1_1["a"]("R,a");
            tmps_["5_LL"]("L,R") += 2.00 * l2_1["bab"]("L,j,a,b") * r2_1["abb"]("R,a,b,j");
            tmps_["5_LL"]("L,R") += l2_1["bbb"]("L,j,d,b") * r2_1["bbb"]("R,d,b,j");
            tmps_["5_LL"]("L,R") += l2["aaa"]("L,i,a,c") * r2["aaa"]("R,a,c,i");
            tmps_["5_LL"]("L,R") += 2.00 * l1_1["b"]("L,d") * r1_1["b"]("R,d");
            tmps_["5_LL"]("L,R") += 2.00 * l2["aab"]("L,i,a,b") * r2["aba"]("R,a,b,i");

            // D_oo_aa  = +0.50 d_aa(i,j) r2_1_aaa(a,b,k) l2_1_aaa(k,a,b)
            //         += +1.00 d_aa(i,j) r1_a(a) l1_a(a)
            //         += +0.50 d_aa(i,j) r2_1_aba(a,b,k) l2_1_aab(k,a,b)
            //         += +0.50 d_aa(i,j) r2_1_aba(b,a,k) l2_1_aab(k,b,a)
            //         += +0.50 d_aa(i,j) r2_abb(a,b,k) l2_bab(k,a,b)
            //         += +0.50 d_aa(i,j) r2_abb(b,a,k) l2_bab(k,b,a)
            //         += +0.50 d_aa(i,j) r2_bbb(a,b,k) l2_bbb(k,a,b)
            //         += +1.00 d_aa(i,j) r1_b(a) l1_b(a)
            //         += +1.00 d_aa(i,j) r1_1_a(a) l1_1_a(a)
            //         += +0.50 d_aa(i,j) r2_1_abb(a,b,k) l2_1_bab(k,a,b)
            //         += +0.50 d_aa(i,j) r2_1_abb(b,a,k) l2_1_bab(k,b,a)
            //         += +0.50 d_aa(i,j) r2_1_bbb(a,b,k) l2_1_bbb(k,a,b)
            //         += +0.50 d_aa(i,j) r2_aaa(a,b,k) l2_aaa(k,a,b)
            //         += +1.00 d_aa(i,j) r1_1_b(a) l1_1_b(a)
            //         += +0.50 d_aa(i,j) r2_aba(a,b,k) l2_aab(k,a,b)
            //         += +0.50 d_aa(i,j) r2_aba(b,a,k) l2_aab(k,b,a)
            // flops: o2v0L2  = o2v0L2
            //  mems: o2v0L2  = o2v0L2
            D_oo_aa("L,R,i,j")  = 0.50 * Id["aa_oo"]("i,j") * tmps_["5_LL"]("L,R");

            // D_oo_bb  = +0.50 d_bb(i,j) r2_1_aaa(a,b,k) l2_1_aaa(k,a,b)
            //         += +1.00 d_bb(i,j) r1_a(a) l1_a(a)
            //         += +0.50 d_bb(i,j) r2_1_aba(a,b,k) l2_1_aab(k,a,b)
            //         += +0.50 d_bb(i,j) r2_1_aba(b,a,k) l2_1_aab(k,b,a)
            //         += +0.50 d_bb(i,j) r2_abb(a,b,k) l2_bab(k,a,b)
            //         += +0.50 d_bb(i,j) r2_abb(b,a,k) l2_bab(k,b,a)
            //         += +0.50 d_bb(i,j) r2_bbb(a,b,k) l2_bbb(k,a,b)
            //         += +1.00 d_bb(i,j) r1_b(a) l1_b(a)
            //         += +1.00 d_bb(i,j) r1_1_a(a) l1_1_a(a)
            //         += +0.50 d_bb(i,j) r2_1_abb(a,b,k) l2_1_bab(k,a,b)
            //         += +0.50 d_bb(i,j) r2_1_abb(b,a,k) l2_1_bab(k,b,a)
            //         += +0.50 d_bb(i,j) r2_1_bbb(a,b,k) l2_1_bbb(k,a,b)
            //         += +0.50 d_bb(i,j) r2_aaa(a,b,k) l2_aaa(k,a,b)
            //         += +1.00 d_bb(i,j) r1_1_b(a) l1_1_b(a)
            //         += +0.50 d_bb(i,j) r2_aba(a,b,k) l2_aab(k,a,b)
            //         += +0.50 d_bb(i,j) r2_aba(b,a,k) l2_aab(k,b,a)
            // flops: o2v0L2  = o2v0L2
            //  mems: o2v0L2  = o2v0L2
            D_oo_bb("L,R,i,j")  = 0.50 * Id["bb_oo"]("i,j") * tmps_["5_LL"]("L,R");

            // D_ov_aa += +0.50 r2_1_aaa(b,c,j) t1_aa(a,i) l2_1_aaa(j,b,c)
            //         += +1.00 r1_a(b) t1_aa(a,i) l1_a(b)
            //         += +0.50 r2_1_aba(b,c,j) t1_aa(a,i) l2_1_aab(j,b,c)
            //         += +0.50 r2_1_aba(c,b,j) t1_aa(a,i) l2_1_aab(j,c,b)
            //         += +0.50 r2_abb(b,c,j) t1_aa(a,i) l2_bab(j,b,c)
            //         += +0.50 r2_abb(c,b,j) t1_aa(a,i) l2_bab(j,c,b)
            //         += +0.50 r2_bbb(b,c,j) t1_aa(a,i) l2_bbb(j,b,c)
            //         += +1.00 r1_b(b) t1_aa(a,i) l1_b(b)
            //         += +1.00 r1_1_a(b) t1_aa(a,i) l1_1_a(b)
            //         += +0.50 r2_1_abb(b,c,j) t1_aa(a,i) l2_1_bab(j,b,c)
            //         += +0.50 r2_1_abb(c,b,j) t1_aa(a,i) l2_1_bab(j,c,b)
            //         += +0.50 r2_1_bbb(b,c,j) t1_aa(a,i) l2_1_bbb(j,b,c)
            //         += +0.50 r2_aaa(b,c,j) t1_aa(a,i) l2_aaa(j,b,c)
            //         += +1.00 r1_1_b(b) t1_aa(a,i) l1_1_b(b)
            //         += +0.50 r2_aba(b,c,j) t1_aa(a,i) l2_aab(j,b,c)
            //         += +0.50 r2_aba(c,b,j) t1_aa(a,i) l2_aab(j,c,b)
            // flops: o1v1L2 += o1v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") += 0.50 * t1["aa"]("a,i") * tmps_["5_LL"]("L,R");

            // D_ov_bb += +0.50 r2_1_aaa(b,c,j) t1_bb(a,i) l2_1_aaa(j,b,c)
            //         += +1.00 r1_a(b) t1_bb(a,i) l1_a(b)
            //         += +0.50 r2_1_aba(b,c,j) t1_bb(a,i) l2_1_aab(j,b,c)
            //         += +0.50 r2_1_aba(c,b,j) t1_bb(a,i) l2_1_aab(j,c,b)
            //         += +0.50 r2_abb(b,c,j) t1_bb(a,i) l2_bab(j,b,c)
            //         += +0.50 r2_abb(c,b,j) t1_bb(a,i) l2_bab(j,c,b)
            //         += +0.50 r2_bbb(b,c,j) t1_bb(a,i) l2_bbb(j,b,c)
            //         += +1.00 r1_b(b) t1_bb(a,i) l1_b(b)
            //         += +1.00 r1_1_a(b) t1_bb(a,i) l1_1_a(b)
            //         += +0.50 r2_1_abb(b,c,j) t1_bb(a,i) l2_1_bab(j,b,c)
            //         += +0.50 r2_1_abb(c,b,j) t1_bb(a,i) l2_1_bab(j,c,b)
            //         += +0.50 r2_1_bbb(b,c,j) t1_bb(a,i) l2_1_bbb(j,b,c)
            //         += +0.50 r2_aaa(b,c,j) t1_bb(a,i) l2_aaa(j,b,c)
            //         += +1.00 r1_1_b(b) t1_bb(a,i) l1_1_b(b)
            //         += +0.50 r2_aba(b,c,j) t1_bb(a,i) l2_aab(j,b,c)
            //         += +0.50 r2_aba(c,b,j) t1_bb(a,i) l2_aab(j,c,b)
            // flops: o1v1L2 += o1v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") += 0.50 * t1["bb"]("a,i") * tmps_["5_LL"]("L,R");
            tmps_["5_LL"].~TArrayD();

            // flops: o0v0L2  = o1v2L2 o1v2L2 o0v1L2 o1v2L2 o0v0L2 o0v0L2 o0v0L2 o0v1L2 o0v0L2 o1v2L2 o0v0L2
            //  mems: o0v0L2  = o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2
            tmps_["6_LL"]("L,R")  = 0.50 * l2_1["aaa"]("L,j,a,d") * r2["aaa"]("R,a,d,j");
            tmps_["6_LL"]("L,R") += 0.50 * l2_1["bbb"]("L,i,c,b") * r2["bbb"]("R,c,b,i");
            tmps_["6_LL"]("L,R") += l1_1["a"]("L,a") * r1["a"]("R,a");
            tmps_["6_LL"]("L,R") += l2_1["aab"]("L,j,a,b") * r2["aba"]("R,a,b,j");
            tmps_["6_LL"]("L,R") += l1_1["b"]("L,c") * r1["b"]("R,c");
            tmps_["6_LL"]("L,R") += l2_1["bab"]("L,i,a,b") * r2["abb"]("R,a,b,i");

            // D_ov_aa += +0.50 r2_aba(b,c,j) t1_1_aa(a,i) l2_1_aab(j,b,c)
            //         += +0.50 r2_aba(c,b,j) t1_1_aa(a,i) l2_1_aab(j,c,b)
            //         += +1.00 r1_a(b) t1_1_aa(a,i) l1_1_a(b)
            //         += +0.50 r2_bbb(b,c,j) t1_1_aa(a,i) l2_1_bbb(j,b,c)
            //         += +0.50 r2_aaa(b,c,j) t1_1_aa(a,i) l2_1_aaa(j,b,c)
            //         += +1.00 r1_b(b) t1_1_aa(a,i) l1_1_b(b)
            //         += +0.50 r2_abb(b,c,j) t1_1_aa(a,i) l2_1_bab(j,b,c)
            //         += +0.50 r2_abb(c,b,j) t1_1_aa(a,i) l2_1_bab(j,c,b)
            // flops: o1v1L2 += o1v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") += t1_1["aa"]("a,i") * tmps_["6_LL"]("L,R");

            // D_ov_bb += +0.50 r2_aba(b,c,j) t1_1_bb(a,i) l2_1_aab(j,b,c)
            //         += +0.50 r2_aba(c,b,j) t1_1_bb(a,i) l2_1_aab(j,c,b)
            //         += +1.00 r1_a(b) t1_1_bb(a,i) l1_1_a(b)
            //         += +0.50 r2_bbb(b,c,j) t1_1_bb(a,i) l2_1_bbb(j,b,c)
            //         += +0.50 r2_aaa(b,c,j) t1_1_bb(a,i) l2_1_aaa(j,b,c)
            //         += +1.00 r1_b(b) t1_1_bb(a,i) l1_1_b(b)
            //         += +0.50 r2_abb(b,c,j) t1_1_bb(a,i) l2_1_bab(j,b,c)
            //         += +0.50 r2_abb(c,b,j) t1_1_bb(a,i) l2_1_bab(j,c,b)
            // flops: o1v1L2 += o1v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") += t1_1["bb"]("a,i") * tmps_["6_LL"]("L,R");
            tmps_["6_LL"].~TArrayD();

            // flops: o0v2L2  = o1v3L2 o1v3L2 o1v3L2 o0v2L2 o1v3L2 o0v2L2 o1v3L2 o0v2L2 o0v2L2 o1v3L2 o0v2L2
            //  mems: o0v2L2  = o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2
            tmps_["7_aa_LLvv"]("L,R,a,b")  = l2["aaa"]("L,i,a,d") * r2["aaa"]("R,b,d,i");
            tmps_["7_aa_LLvv"]("L,R,a,b") += l2["bab"]("L,j,a,c") * r2["abb"]("R,b,c,j");
            tmps_["7_aa_LLvv"]("L,R,a,b") += l2_1["aaa"]("L,i,a,d") * r2_1["aaa"]("R,b,d,i");
            tmps_["7_aa_LLvv"]("L,R,a,b") += l2_1["bab"]("L,j,a,c") * r2_1["abb"]("R,b,c,j");
            tmps_["7_aa_LLvv"]("L,R,a,b") += l2_1["aab"]("L,i,a,c") * r2_1["aba"]("R,b,c,i");
            tmps_["7_aa_LLvv"]("L,R,a,b") += l2["aab"]("L,i,a,c") * r2["aba"]("R,b,c,i");

            // D_ov_aa += -1.00 r2_1_aaa(a,c,j) t1_aa(b,i) l2_1_aaa(j,b,c)
            //         += -1.00 r2_abb(a,c,j) t1_aa(b,i) l2_bab(j,b,c)
            //         += -1.00 r2_1_abb(a,c,j) t1_aa(b,i) l2_1_bab(j,b,c)
            //         += -1.00 r2_1_aba(a,c,j) t1_aa(b,i) l2_1_aab(j,b,c)
            //         += -1.00 r2_aaa(a,c,j) t1_aa(b,i) l2_aaa(j,b,c)
            //         += -1.00 r2_aba(a,c,j) t1_aa(b,i) l2_aab(j,b,c)
            // flops: o1v1L2 += o1v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") -= t1["aa"]("b,i") * tmps_["7_aa_LLvv"]("L,R,b,a");

            // D_vv_aa += +1.00 r2_1_aaa(b,c,i) l2_1_aaa(i,a,c)
            //         += +1.00 r2_abb(b,c,i) l2_bab(i,a,c)
            //         += +1.00 r2_1_abb(b,c,i) l2_1_bab(i,a,c)
            //         += +1.00 r2_1_aba(b,c,i) l2_1_aab(i,a,c)
            //         += +1.00 r2_aaa(b,c,i) l2_aaa(i,a,c)
            //         += +1.00 r2_aba(b,c,i) l2_aab(i,a,c)
            D_vv_aa("L,R,a,b") += tmps_["7_aa_LLvv"]("L,R,a,b");
            tmps_["7_aa_LLvv"].~TArrayD();

            // flops: o0v2L2  = o1v3L2 o1v3L2 o1v3L2 o0v2L2 o1v3L2 o0v2L2 o1v3L2 o0v2L2 o1v3L2 o0v2L2 o0v2L2
            //  mems: o0v2L2  = o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2
            tmps_["8_bb_LLvv"]("L,R,a,b")  = l2["aab"]("L,i,c,a") * r2["aba"]("R,c,b,i");
            tmps_["8_bb_LLvv"]("L,R,a,b") += l2["bab"]("L,j,c,a") * r2["abb"]("R,c,b,j");
            tmps_["8_bb_LLvv"]("L,R,a,b") += l2_1["bab"]("L,j,c,a") * r2_1["abb"]("R,c,b,j");
            tmps_["8_bb_LLvv"]("L,R,a,b") += l2_1["aab"]("L,i,c,a") * r2_1["aba"]("R,c,b,i");
            tmps_["8_bb_LLvv"]("L,R,a,b") += l2["bbb"]("L,j,a,d") * r2["bbb"]("R,b,d,j");
            tmps_["8_bb_LLvv"]("L,R,a,b") += l2_1["bbb"]("L,j,a,d") * r2_1["bbb"]("R,b,d,j");

            // D_ov_bb += -1.00 r2_abb(c,a,j) t1_bb(b,i) l2_bab(j,c,b)
            //         += -1.00 r2_1_abb(c,a,j) t1_bb(b,i) l2_1_bab(j,c,b)
            //         += -1.00 r2_1_aba(c,a,j) t1_bb(b,i) l2_1_aab(j,c,b)
            //         += -1.00 r2_bbb(a,c,j) t1_bb(b,i) l2_bbb(j,b,c)
            //         += -1.00 r2_1_bbb(a,c,j) t1_bb(b,i) l2_1_bbb(j,b,c)
            //         += -1.00 r2_aba(c,a,j) t1_bb(b,i) l2_aab(j,c,b)
            // flops: o1v1L2 += o1v2L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") -= t1["bb"]("b,i") * tmps_["8_bb_LLvv"]("L,R,b,a");

            // D_vv_bb += +1.00 r2_abb(c,b,i) l2_bab(i,c,a)
            //         += +1.00 r2_1_abb(c,b,i) l2_1_bab(i,c,a)
            //         += +1.00 r2_1_aba(c,b,i) l2_1_aab(i,c,a)
            //         += +1.00 r2_bbb(b,c,i) l2_bbb(i,a,c)
            //         += +1.00 r2_1_bbb(b,c,i) l2_1_bbb(i,a,c)
            //         += +1.00 r2_aba(c,b,i) l2_aab(i,c,a)
            D_vv_bb("L,R,a,b") += tmps_["8_bb_LLvv"]("L,R,a,b");
            tmps_["8_bb_LLvv"].~TArrayD();

            // flops: o1v1L2  = o1v2L2 o1v2L2 o1v1L2
            //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2
            tmps_["9_aa_LLov"]("L,R,i,a")  = l2["aab"]("L,i,a,b") * r1["b"]("R,b");
            tmps_["9_aa_LLov"]("L,R,i,a") += l2_1["aab"]("L,i,a,b") * r1_1["b"]("R,b");

            // D_vo_aa  = -1.00 r1_1_b(b) l2_1_aab(i,a,b)
            //         += -1.00 r1_b(b) l2_aab(i,a,b)
            D_vo_aa("L,R,a,i")  = -1.00 * tmps_["9_aa_LLov"]("L,R,i,a");

            // D_vv_aa += -1.00 r1_1_b(c) t1_aa(b,i) l2_1_aab(i,a,c)
            //         += -1.00 r1_b(c) t1_aa(b,i) l2_aab(i,a,c)
            // flops: o0v2L2 += o1v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_aa("L,R,a,b") -= t1["aa"]("b,i") * tmps_["9_aa_LLov"]("L,R,i,a");
            tmps_["9_aa_LLov"].~TArrayD();

            // flops: o1v1L2  = o1v2L2 o1v2L2 o1v1L2
            //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2
            tmps_["10_aa_LLov"]("L,R,i,a")  = l2["aaa"]("L,i,a,b") * r1["a"]("R,b");
            tmps_["10_aa_LLov"]("L,R,i,a") += l2_1["aaa"]("L,i,a,b") * r1_1["a"]("R,b");

            // D_vo_aa += -1.00 r1_1_a(b) l2_1_aaa(i,a,b)
            //         += -1.00 r1_a(b) l2_aaa(i,a,b)
            D_vo_aa("L,R,a,i") -= tmps_["10_aa_LLov"]("L,R,i,a");

            // D_vv_aa += -1.00 r1_1_a(c) t1_aa(b,i) l2_1_aaa(i,a,c)
            //         += -1.00 r1_a(c) t1_aa(b,i) l2_aaa(i,a,c)
            // flops: o0v2L2 += o1v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_aa("L,R,a,b") -= t1["aa"]("b,i") * tmps_["10_aa_LLov"]("L,R,i,a");

            // flops: o2v0L2  = o2v2L2 o2v2L2 o2v0L2 o2v2L2 o2v0L2 o2v2L1 o2v1L2 o2v0L2 o2v2L1 o2v1L2 o2v0L2 o2v2L1 o2v1L2 o2v0L2 o2v2L2 o2v0L2 o2v2L1 o2v1L2 o2v0L2 o2v1L2 o2v0L2
            //  mems: o2v0L2  = o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v1L1 o2v0L2 o2v0L2 o2v1L1 o2v0L2 o2v0L2 o2v1L1 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v1L1 o2v0L2 o2v0L2 o2v0L2 o2v0L2
            tmps_["11_aa_LLoo"]("L,R,i,j")  = -0.50 * l2_1["aaa"]("L,j,a,b") * r2_1["aaa"]("R,a,b,i");
            tmps_["11_aa_LLoo"]("L,R,i,j") -= l2["aab"]("L,j,a,c") * r2["aba"]("R,a,c,i");
            tmps_["11_aa_LLoo"]("L,R,i,j") -= 0.50 * l2["aaa"]("L,j,a,b") * r2["aaa"]("R,a,b,i");
            tmps_["11_aa_LLoo"]("L,R,i,j") += t1_1["aa"]("a,i") * l2_1["aaa"]("L,j,a,b") * r1["a"]("R,b");
            tmps_["11_aa_LLoo"]("L,R,i,j") += t1_1["aa"]("a,i") * l2_1["aab"]("L,j,a,c") * r1["b"]("R,c");
            tmps_["11_aa_LLoo"]("L,R,i,j") += t1["aa"]("a,i") * l2_1["aab"]("L,j,a,c") * r1_1["b"]("R,c");
            tmps_["11_aa_LLoo"]("L,R,i,j") -= l2_1["aab"]("L,j,a,c") * r2_1["aba"]("R,a,c,i");
            tmps_["11_aa_LLoo"]("L,R,i,j") += t1["aa"]("a,i") * l2["aab"]("L,j,a,c") * r1["b"]("R,c");
            tmps_["11_aa_LLoo"]("L,R,i,j") += t1["aa"]("a,i") * tmps_["10_aa_LLov"]("L,R,j,a");
            tmps_["10_aa_LLov"].~TArrayD();

            // D_oo_aa += +1.00 r1_1_a(b) t1_aa(a,i) l2_1_aaa(j,a,b)
            //         += +1.00 r1_a(b) t1_aa(a,i) l2_aaa(j,a,b)
            //         += -0.50 r2_1_aaa(a,b,i) l2_1_aaa(j,a,b)
            //         += -0.50 r2_aba(a,b,i) l2_aab(j,a,b)
            //         += -0.50 r2_aba(b,a,i) l2_aab(j,b,a)
            //         += -0.50 r2_aaa(a,b,i) l2_aaa(j,a,b)
            //         += +1.00 r1_a(b) t1_1_aa(a,i) l2_1_aaa(j,a,b)
            //         += +1.00 r1_b(b) t1_1_aa(a,i) l2_1_aab(j,a,b)
            //         += +1.00 r1_1_b(b) t1_aa(a,i) l2_1_aab(j,a,b)
            //         += -0.50 r2_1_aba(a,b,i) l2_1_aab(j,a,b)
            //         += -0.50 r2_1_aba(b,a,i) l2_1_aab(j,b,a)
            //         += +1.00 r1_b(b) t1_aa(a,i) l2_aab(j,a,b)
            D_oo_aa("L,R,i,j") += tmps_["11_aa_LLoo"]("L,R,i,j");

            // D_ov_aa += +1.00 r1_1_a(c) t1_aa(a,j) t1_aa(b,i) l2_1_aaa(j,b,c)
            //         += +1.00 r1_a(c) t1_aa(a,j) t1_aa(b,i) l2_aaa(j,b,c)
            //         += -0.50 r2_1_aaa(b,c,i) t1_aa(a,j) l2_1_aaa(j,b,c)
            //         += -0.50 r2_aba(b,c,i) t1_aa(a,j) l2_aab(j,b,c)
            //         += -0.50 r2_aba(c,b,i) t1_aa(a,j) l2_aab(j,c,b)
            //         += -0.50 r2_aaa(b,c,i) t1_aa(a,j) l2_aaa(j,b,c)
            //         += +1.00 r1_a(c) t1_aa(a,j) t1_1_aa(b,i) l2_1_aaa(j,b,c)
            //         += +1.00 r1_b(c) t1_aa(a,j) t1_1_aa(b,i) l2_1_aab(j,b,c)
            //         += +1.00 r1_1_b(c) t1_aa(a,j) t1_aa(b,i) l2_1_aab(j,b,c)
            //         += -0.50 r2_1_aba(b,c,i) t1_aa(a,j) l2_1_aab(j,b,c)
            //         += -0.50 r2_1_aba(c,b,i) t1_aa(a,j) l2_1_aab(j,c,b)
            //         += +1.00 r1_b(c) t1_aa(a,j) t1_aa(b,i) l2_aab(j,b,c)
            // flops: o1v1L2 += o2v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_aa("L,R,i,a") += t1["aa"]("a,j") * tmps_["11_aa_LLoo"]("L,R,i,j");
            tmps_["11_aa_LLoo"].~TArrayD();

            // flops: o1v1L2  = o1v2L2
            //  mems: o1v1L2  = o1v1L2
            tmps_["12_bb_LLov"]("L,R,i,a")  = l2_1["bab"]("L,i,b,a") * r1["a"]("R,b");

            // D_ov_bb += -1.00 r1_a(c) t1_bb(b,i) t1_1_bb(a,j) l2_1_bab(j,c,b)
            // flops: o1v1L2 += o2v1L2 o2v1L2
            //  mems: o1v1L2 += o2v0L2 o1v1L2
            D_ov_bb("L,R,i,a") -= t1["bb"]("b,i") * tmps_["12_bb_LLov"]("L,R,j,b") * t1_1["bb"]("a,j");

            // D_vv_bb += +1.00 r1_a(c) t1_1_bb(b,i) l2_1_bab(i,c,a)
            // flops: o0v2L2 += o1v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_bb("L,R,a,b") += t1_1["bb"]("b,i") * tmps_["12_bb_LLov"]("L,R,i,a");

            // flops: o1v1L2  = o1v2L2 o1v2L2 o1v2L2 o1v1L2 o1v2L2 o1v1L2 o1v1L2
            //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2
            tmps_["13_bb_LLov"]("L,R,i,a")  = -1.00 * l2["bab"]("L,i,c,a") * r1["a"]("R,c");
            tmps_["13_bb_LLov"]("L,R,i,a") -= l2_1["bab"]("L,i,c,a") * r1_1["a"]("R,c");
            tmps_["13_bb_LLov"]("L,R,i,a") += l2_1["bbb"]("L,i,a,b") * r1_1["b"]("R,b");
            tmps_["13_bb_LLov"]("L,R,i,a") += l2["bbb"]("L,i,a,b") * r1["b"]("R,b");

            // D_vo_bb  = -1.00 r1_1_b(b) l2_1_bbb(i,a,b)
            //         += +1.00 r1_1_a(b) l2_1_bab(i,b,a)
            //         += -1.00 r1_b(b) l2_bbb(i,a,b)
            //         += +1.00 r1_a(b) l2_bab(i,b,a)
            D_vo_bb("L,R,a,i")  = -1.00 * tmps_["13_bb_LLov"]("L,R,i,a");

            // D_vv_bb += -1.00 r1_1_b(c) t1_bb(b,i) l2_1_bbb(i,a,c)
            //         += +1.00 r1_1_a(c) t1_bb(b,i) l2_1_bab(i,c,a)
            //         += -1.00 r1_b(c) t1_bb(b,i) l2_bbb(i,a,c)
            //         += +1.00 r1_a(c) t1_bb(b,i) l2_bab(i,c,a)
            // flops: o0v2L2 += o1v2L2
            //  mems: o0v2L2 += o0v2L2
            D_vv_bb("L,R,a,b") -= t1["bb"]("b,i") * tmps_["13_bb_LLov"]("L,R,i,a");

            // flops: o2v0L2  = o2v1L2 o2v2L2 o2v2L1 o2v1L2 o2v0L2 o2v2L2 o2v0L2 o2v2L2 o2v0L2 o2v2L2 o2v0L2 o2v1L2 o2v0L2 o2v0L2
            //  mems: o2v0L2  = o2v0L2 o2v0L2 o2v1L1 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2
            tmps_["14_bb_LLoo"]("L,R,i,j")  = -1.00 * t1_1["bb"]("a,i") * tmps_["12_bb_LLov"]("L,R,j,a");
            tmps_["14_bb_LLoo"]("L,R,i,j") -= 0.50 * l2["bbb"]("L,j,a,d") * r2["bbb"]("R,a,d,i");
            tmps_["14_bb_LLoo"]("L,R,i,j") += l2_1["bbb"]("L,j,a,d") * t1_1["bb"]("a,i") * r1["b"]("R,d");
            tmps_["14_bb_LLoo"]("L,R,i,j") -= 0.50 * l2_1["bbb"]("L,j,a,d") * r2_1["bbb"]("R,a,d,i");
            tmps_["14_bb_LLoo"]("L,R,i,j") -= l2_1["bab"]("L,j,c,d") * r2_1["abb"]("R,c,d,i");
            tmps_["14_bb_LLoo"]("L,R,i,j") -= l2["bab"]("L,j,c,d") * r2["abb"]("R,c,d,i");
            tmps_["14_bb_LLoo"]("L,R,i,j") += t1["bb"]("a,i") * tmps_["13_bb_LLov"]("L,R,j,a");
            tmps_["13_bb_LLov"].~TArrayD();
            tmps_["12_bb_LLov"].~TArrayD();

            // D_oo_bb += +1.00 r1_1_b(b) t1_bb(a,i) l2_1_bbb(j,a,b)
            //         += -1.00 r1_1_a(b) t1_bb(a,i) l2_1_bab(j,b,a)
            //         += +1.00 r1_b(b) t1_bb(a,i) l2_bbb(j,a,b)
            //         += -1.00 r1_a(b) t1_bb(a,i) l2_bab(j,b,a)
            //         += +1.00 r1_b(b) t1_1_bb(a,i) l2_1_bbb(j,a,b)
            //         += -0.50 r2_bbb(a,b,i) l2_bbb(j,a,b)
            //         += -0.50 r2_1_bbb(a,b,i) l2_1_bbb(j,a,b)
            //         += -0.50 r2_1_abb(a,b,i) l2_1_bab(j,a,b)
            //         += -0.50 r2_1_abb(b,a,i) l2_1_bab(j,b,a)
            //         += -0.50 r2_abb(a,b,i) l2_bab(j,a,b)
            //         += -0.50 r2_abb(b,a,i) l2_bab(j,b,a)
            //         += -1.00 r1_a(b) t1_1_bb(a,i) l2_1_bab(j,b,a)
            D_oo_bb("L,R,i,j") += tmps_["14_bb_LLoo"]("L,R,i,j");

            // D_ov_bb += +1.00 r1_1_b(c) t1_bb(a,j) t1_bb(b,i) l2_1_bbb(j,b,c)
            //         += -1.00 r1_1_a(c) t1_bb(a,j) t1_bb(b,i) l2_1_bab(j,c,b)
            //         += +1.00 r1_b(c) t1_bb(a,j) t1_bb(b,i) l2_bbb(j,b,c)
            //         += -1.00 r1_a(c) t1_bb(a,j) t1_bb(b,i) l2_bab(j,c,b)
            //         += +1.00 r1_b(c) t1_bb(a,j) t1_1_bb(b,i) l2_1_bbb(j,b,c)
            //         += -0.50 r2_bbb(b,c,i) t1_bb(a,j) l2_bbb(j,b,c)
            //         += -0.50 r2_1_bbb(b,c,i) t1_bb(a,j) l2_1_bbb(j,b,c)
            //         += -0.50 r2_1_abb(b,c,i) t1_bb(a,j) l2_1_bab(j,b,c)
            //         += -0.50 r2_1_abb(c,b,i) t1_bb(a,j) l2_1_bab(j,c,b)
            //         += -0.50 r2_abb(b,c,i) t1_bb(a,j) l2_bab(j,b,c)
            //         += -0.50 r2_abb(c,b,i) t1_bb(a,j) l2_bab(j,c,b)
            //         += -1.00 r1_a(c) t1_bb(a,j) t1_1_bb(b,i) l2_1_bab(j,c,b)
            // flops: o1v1L2 += o2v1L2
            //  mems: o1v1L2 += o1v1L2
            D_ov_bb("L,R,i,a") += t1["bb"]("a,j") * tmps_["14_bb_LLoo"]("L,R,i,j");
            tmps_["14_bb_LLoo"].~TArrayD();
        }

        rdm_timer.stop();
        Printf("Done --> Time: %s\n", rdm_timer.elapsed().c_str());

    }

} // cc_cavity
