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

#ifdef BUILD_QED_2RDM
#include <psi4/libmints/writer.h>
#include "cc_cavity/include/qed_ccsd_21/eom_ee_qed_rdm_21.h"
#include "misc/nonsym_davidson_solver.h"

namespace hilbert {


    void EOM_EE_QED_RDM_21::rdm2_21_3() {

        auto &cc_wfn_ = eom_driver_->cc_wfn_;

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

        TArrayMap &Id = cc_wfn_->Id_blks_;
        TArrayMap &evecs = eom_driver_->evec_blks_;


        /// unpack left/right operators

        TArrayMap   r1,   r2,   l1,   l2;
        TArrayMap r1_1, r2_1, l1_1, l2_1;


        // slice out the relevant states
        TArrayD &l0 = evecs["l0"];
        TArrayD &r0 = evecs["r0"];

        l1["aa"] = evecs["l1_aa"];
        l1["bb"] = evecs["l1_bb"];
        r1["aa"] = evecs["r1_aa"];
        r1["bb"] = evecs["r1_bb"];

        l2["aaaa"] = evecs["l2_aaaa"];
        l2["abab"] = evecs["l2_abab"];
        l2["bbbb"] = evecs["l2_bbbb"];
        r2["aaaa"] = evecs["r2_aaaa"];
        r2["abab"] = evecs["r2_abab"];
        r2["bbbb"] = evecs["r2_bbbb"];

        TArrayD &l0_1 = evecs["l0_1"];
        TArrayD &r0_1 = evecs["r0_1"];

        l1_1["aa"] = evecs["l1_1_aa"];
        l1_1["bb"] = evecs["l1_1_bb"];
        r1_1["aa"] = evecs["r1_1_aa"];
        r1_1["bb"] = evecs["r1_1_bb"];

        l2_1["aaaa"] = evecs["l2_1_aaaa"];
        l2_1["abab"] = evecs["l2_1_abab"];
        l2_1["bbbb"] = evecs["l2_1_bbbb"];
        r2_1["aaaa"] = evecs["r2_1_aaaa"];
        r2_1["abab"] = evecs["r2_1_abab"];
        r2_1["bbbb"] = evecs["r2_1_bbbb"];

        TArrayD &D_oooo_aaaa = RDM_blks_["D2_oooo_aaaa"];
        TArrayD &D_oooo_abab = RDM_blks_["D2_oooo_abab"];
        TArrayD &D_oooo_bbbb = RDM_blks_["D2_oooo_bbbb"];
        TArrayD &D_ooov_aaaa = RDM_blks_["D2_ooov_aaaa"];
        TArrayD &D_ooov_abab = RDM_blks_["D2_ooov_abab"];
        TArrayD &D_ooov_bbbb = RDM_blks_["D2_ooov_bbbb"];
        TArrayD &D_oovo_aaaa = RDM_blks_["D2_oovo_aaaa"];
        TArrayD &D_oovo_abab = RDM_blks_["D2_oovo_abab"];
        TArrayD &D_oovo_bbbb = RDM_blks_["D2_oovo_bbbb"];
        TArrayD &D_oovv_aaaa = RDM_blks_["D2_oovv_aaaa"];
        TArrayD &D_oovv_abab = RDM_blks_["D2_oovv_abab"];
        TArrayD &D_oovv_bbbb = RDM_blks_["D2_oovv_bbbb"];
        TArrayD &D_ovov_aaaa = RDM_blks_["D2_ovov_aaaa"];
        TArrayD &D_ovov_abab = RDM_blks_["D2_ovov_abab"];
        TArrayD &D_ovov_bbbb = RDM_blks_["D2_ovov_bbbb"];
        TArrayD &D_ovvo_aaaa = RDM_blks_["D2_ovvo_aaaa"];
        TArrayD &D_ovvo_abab = RDM_blks_["D2_ovvo_abab"];
        TArrayD &D_ovvo_bbbb = RDM_blks_["D2_ovvo_bbbb"];
        TArrayD &D_ovvv_aaaa = RDM_blks_["D2_ovvv_aaaa"];
        TArrayD &D_ovvv_abab = RDM_blks_["D2_ovvv_abab"];
        TArrayD &D_ovvv_bbbb = RDM_blks_["D2_ovvv_bbbb"];
        TArrayD &D_vovv_aaaa = RDM_blks_["D2_vovv_aaaa"];
        TArrayD &D_vovv_abab = RDM_blks_["D2_vovv_abab"];
        TArrayD &D_vovv_bbbb = RDM_blks_["D2_vovv_bbbb"];
        TArrayD &D_vvov_aaaa = RDM_blks_["D2_vvov_aaaa"];
        TArrayD &D_vvov_abab = RDM_blks_["D2_vvov_abab"];
        TArrayD &D_vvov_bbbb = RDM_blks_["D2_vvov_bbbb"];
        TArrayD &D_vvvo_aaaa = RDM_blks_["D2_vvvo_aaaa"];
        TArrayD &D_vvvo_abab = RDM_blks_["D2_vvvo_abab"];
        TArrayD &D_vvvo_bbbb = RDM_blks_["D2_vvvo_bbbb"];
        TArrayD &D_vvvv_aaaa = RDM_blks_["D2_vvvv_aaaa"];
        TArrayD &D_vvvv_abab = RDM_blks_["D2_vvvv_abab"];
        TArrayD &D_vvvv_bbbb = RDM_blks_["D2_vvvv_bbbb"];


        // D_oovv_bbbb += -0.50 P(a,b) r2_abab(d,a,l,k) t2_1_bbbb(b,c,i,j) l2_1_abab(l,k,d,c)
        //             += -0.50 P(a,b) r2_abab(d,a,k,l) t2_1_bbbb(b,c,i,j) l2_1_abab(k,l,d,c)
        //             += -0.50 P(a,b) r2_bbbb(a,d,l,k) t2_1_bbbb(b,c,i,j) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["98_bb_LLvv"]("L,R,c,a") * t2_1["bbbb"]("b,c,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v1L2  = o1v2L2 o2v1L2 o2v1L1 o2v1L2 o1v1L2 o2v1L2 o1v1L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o1v1L2 o2v0L1 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2
        tmps_["99_bb_LLov"]("L,R,i,a")  = tmps_["98_bb_LLvv"]("L,R,b,a") * t1["bb"]("b,i");
        tmps_["99_bb_LLov"]("L,R,i,a") -= 0.50 * tmps_["97_bb_Loo"]("L,i,j") * r1["bb"]("R,a,j");
        tmps_["99_bb_LLov"]("L,R,i,a") += l1_1["bb"]("L,k,b") * t1["bb"]("b,i") * r1["bb"]("R,a,k");
        tmps_["99_bb_LLov"]("L,R,i,a") += tmps_["96_bb_Loo"]("L,i,j") * r1["bb"]("R,a,j");

        // D_oovv_abab += +1.00 r1_bb(b,k) t1_bb(c,j) t1_1_aa(a,i) l1_1_bb(k,c)
        //             += -0.50 r1_bb(b,l) t2_bbbb(d,c,k,j) t1_1_aa(a,i) l2_1_bbbb(l,k,d,c)
        //             += +0.50 r1_bb(b,l) t2_abab(d,c,k,j) t1_1_aa(a,i) l2_1_abab(k,l,d,c)
        //             += +0.50 r1_bb(b,l) t2_abab(c,d,k,j) t1_1_aa(a,i) l2_1_abab(k,l,c,d)
        //             += +0.50 r2_abab(d,b,l,k) t1_bb(c,j) t1_1_aa(a,i) l2_1_abab(l,k,d,c)
        //             += +0.50 r2_abab(d,b,k,l) t1_bb(c,j) t1_1_aa(a,i) l2_1_abab(k,l,d,c)
        //             += +0.50 r2_bbbb(b,d,l,k) t1_bb(c,j) t1_1_aa(a,i) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1_1["aa"]("a,i") * tmps_["99_bb_LLov"]("L,R,j,b");

        // D_oovv_bbbb += -1.00 P(i,j) P(a,b) r1_bb(a,k) t1_bb(c,j) t1_1_bb(b,i) l1_1_bb(k,c)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,l) t2_bbbb(d,c,k,j) t1_1_bb(b,i) l2_1_bbbb(l,k,d,c)
        //             += -0.50 P(i,j) P(a,b) r1_bb(a,l) t2_abab(d,c,k,j) t1_1_bb(b,i) l2_1_abab(k,l,d,c)
        //             += -0.50 P(i,j) P(a,b) r1_bb(a,l) t2_abab(c,d,k,j) t1_1_bb(b,i) l2_1_abab(k,l,c,d)
        //             += -0.50 P(i,j) P(a,b) r2_abab(d,a,l,k) t1_bb(c,j) t1_1_bb(b,i) l2_1_abab(l,k,d,c)
        //             += -0.50 P(i,j) P(a,b) r2_abab(d,a,k,l) t1_bb(c,j) t1_1_bb(b,i) l2_1_abab(k,l,d,c)
        //             += -0.50 P(i,j) P(a,b) r2_bbbb(a,d,l,k) t1_bb(c,j) t1_1_bb(b,i) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1_1["bb"]("b,i") * tmps_["99_bb_LLov"]("L,R,j,a");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
        tmps_["99_bb_LLov"].~TArrayD();

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["100_bb_Lvv"]("L,a,b")  = t2["bbbb"]("a,c,i,j") * l2_1["bbbb"]("L,i,j,c,b");

        // D_oovv_abab += -0.50 r2_1_abab(a,d,i,j) t2_bbbb(b,c,k,l) l2_1_bbbb(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * tmps_["100_bb_Lvv"]("L,b,d") * r2_1["abab"]("R,a,d,i,j");

        // D_oovv_abab += -0.50 r0_1 t2_abab(a,c,i,j) t2_bbbb(b,d,k,l) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * t2["abab"]("a,c,i,j") * tmps_["100_bb_Lvv"]("L,b,c") * r0_1("R");

        // D_oovv_abab += -0.50 r0 t2_bbbb(b,d,k,l) t2_1_abab(a,c,i,j) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * t2_1["abab"]("a,c,i,j") * tmps_["100_bb_Lvv"]("L,b,c") * r0("R");

        // D_oovv_bbbb += -0.50 P(a,b) r2_1_bbbb(a,d,i,j) t2_bbbb(b,c,k,l) l2_1_bbbb(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["100_bb_Lvv"]("L,b,d") * r2_1["bbbb"]("R,a,d,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += -0.50 r0_1 t2_bbbb(a,c,i,j) t2_bbbb(b,d,k,l) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") -= 0.50 * t2["bbbb"]("a,c,i,j") * tmps_["100_bb_Lvv"]("L,b,c") * r0_1("R");

        // D_oovv_bbbb += -0.50 P(a,b) r0 t2_bbbb(b,d,k,l) t2_1_bbbb(a,c,i,j) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t2_1["bbbb"]("a,c,i,j") * tmps_["100_bb_Lvv"]("L,b,c") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v0L2  = o3v2L2 o3v2L2 o2v0L2 o2v1L2 o2v0L2
        //  mems: o2v0L2  = o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2
        tmps_["101_bb_LLoo"]("L,R,i,j")  = 0.50 * l2_1["bbbb"]("L,l,i,a,c") * r2["bbbb"]("R,a,c,l,j");
        tmps_["101_bb_LLoo"]("L,R,i,j") += l2_1["abab"]("L,k,i,b,c") * r2["abab"]("R,b,c,k,j");
        tmps_["101_bb_LLoo"]("L,R,i,j") += l1_1["bb"]("L,i,a") * r1["bb"]("R,a,j");

        // D_oovv_abab += +0.50 r2_abab(c,d,l,j) t2_1_abab(a,b,i,k) l2_1_abab(l,k,c,d)
        //             += +0.50 r2_abab(d,c,l,j) t2_1_abab(a,b,i,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r2_bbbb(c,d,l,j) t2_1_abab(a,b,i,k) l2_1_bbbb(l,k,c,d)
        //             += +1.00 r1_bb(c,j) t2_1_abab(a,b,i,k) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2_1["abab"]("a,b,i,k") * tmps_["101_bb_LLoo"]("L,R,k,j");

        // D_oovv_bbbb += -0.50 P(i,j) r2_abab(c,d,l,i) t2_1_bbbb(b,a,k,j) l2_1_abab(l,k,c,d)
        //             += -0.50 P(i,j) r2_abab(d,c,l,i) t2_1_bbbb(b,a,k,j) l2_1_abab(l,k,d,c)
        //             += -0.50 P(i,j) r2_bbbb(c,d,l,i) t2_1_bbbb(b,a,k,j) l2_1_bbbb(l,k,c,d)
        //             += -1.00 P(i,j) r1_bb(c,i) t2_1_bbbb(b,a,k,j) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2_1["bbbb"]("b,a,k,j") * tmps_["101_bb_LLoo"]("L,R,k,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v1L2  = o2v1L2 o1v2L2 o1v1L2 o1v2L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2
        tmps_["102_bb_LLvo"]("L,R,a,i")  = -2.00 * t1["bb"]("a,j") * tmps_["101_bb_LLoo"]("L,R,j,i");
        tmps_["102_bb_LLvo"]("L,R,a,i") += r1["bb"]("R,b,i") * tmps_["100_bb_Lvv"]("L,a,b");
        tmps_["102_bb_LLvo"]("L,R,a,i") -= 2.00 * tmps_["84_bb_Lvv"]("L,a,b") * r1["bb"]("R,b,i");

        // D_oovv_abab += -0.50 r1_bb(d,j) t2_bbbb(b,c,k,l) t1_1_aa(a,i) l2_1_bbbb(k,l,c,d)
        //             += +0.50 r2_abab(c,d,l,j) t1_bb(b,k) t1_1_aa(a,i) l2_1_abab(l,k,c,d)
        //             += +0.50 r2_abab(d,c,l,j) t1_bb(b,k) t1_1_aa(a,i) l2_1_abab(l,k,d,c)
        //             += +0.50 r2_bbbb(c,d,l,j) t1_bb(b,k) t1_1_aa(a,i) l2_1_bbbb(l,k,c,d)
        //             += +1.00 r1_bb(c,j) t1_bb(b,k) t1_1_aa(a,i) l1_1_bb(k,c)
        //             += +0.50 r1_bb(d,j) t2_abab(c,b,k,l) t1_1_aa(a,i) l2_1_abab(k,l,c,d)
        //             += +0.50 r1_bb(d,j) t2_abab(c,b,l,k) t1_1_aa(a,i) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * t1_1["aa"]("a,i") * tmps_["102_bb_LLvo"]("L,R,b,j");

        // D_oovv_bbbb += +0.50 P(i,j) P(a,b) r1_bb(d,i) t2_bbbb(b,c,k,l) t1_1_bb(a,j) l2_1_bbbb(k,l,c,d)
        //             += -0.50 P(i,j) P(a,b) r2_abab(c,d,l,i) t1_bb(b,k) t1_1_bb(a,j) l2_1_abab(l,k,c,d)
        //             += -0.50 P(i,j) P(a,b) r2_abab(d,c,l,i) t1_bb(b,k) t1_1_bb(a,j) l2_1_abab(l,k,d,c)
        //             += -0.50 P(i,j) P(a,b) r2_bbbb(c,d,l,i) t1_bb(b,k) t1_1_bb(a,j) l2_1_bbbb(l,k,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_bb(c,i) t1_bb(b,k) t1_1_bb(a,j) l1_1_bb(k,c)
        //             += -0.50 P(i,j) P(a,b) r1_bb(d,i) t2_abab(c,b,k,l) t1_1_bb(a,j) l2_1_abab(k,l,c,d)
        //             += -0.50 P(i,j) P(a,b) r1_bb(d,i) t2_abab(c,b,l,k) t1_1_bb(a,j) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t1_1["bb"]("a,j") * tmps_["102_bb_LLvo"]("L,R,b,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
        tmps_["102_bb_LLvo"].~TArrayD();

        // flops: o2v0L1  = o2v0L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["103_aa_Loo"]("R,i,j")  = Id["aa_oo"]("i,j") * r0("R");

        // D_ooov_aaaa += -0.50 P(i,j) d_aa(j,k) r0 t1_aa(a,m) t2_1_abab(c,b,i,l) l2_1_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t1_aa(a,m) t2_1_abab(b,c,i,l) l2_1_abab(m,l,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t1_aa(a,m) t2_abab(c,b,i,l) l2_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t1_aa(a,m) t2_abab(b,c,i,l) l2_abab(m,l,b,c)
        // flops: o3v1L2 += o3v2L1 o3v2L1 o2v0L1 o2v1L1 o3v1L2
        //  mems: o3v1L2 += o2v0L1 o2v0L1 o2v0L1 o1v1L1 o3v1L2
        tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k")  = t1["aa"]("a,m") * (t2_1["abab"]("c,b,i,l") * l2_1["abab"]("L,m,l,c,b") + -1.00 * t2["abab"]("c,b,i,l") * l2["abab"]("L,m,l,c,b")) * tmps_["103_aa_Loo"]("R,j,k");
        D_ooov_aaaa("L,R,i,j,k,a") -= tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k");
        D_ooov_aaaa("L,R,i,j,k,a") += tmps_["perm_aaaa_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_aaaa_LLvooo"].~TArrayD();

        // D_ooov_aaaa += -0.50 P(i,j) d_aa(j,k) r0 t1_aa(a,m) t2_1_aaaa(c,b,l,i) l2_1_aaaa(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t1_aa(a,m) t2_aaaa(c,b,l,i) l2_aaaa(l,m,c,b)
        // flops: o3v1L2 += o3v2L1 o3v2L1 o2v0L1 o2v1L1 o3v1L2
        //  mems: o3v1L2 += o2v0L1 o2v0L1 o2v0L1 o1v1L1 o3v1L2
        tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k")  = 0.50 * t1["aa"]("a,m") * (t2_1["aaaa"]("c,b,l,i") * l2_1["aaaa"]("L,l,m,c,b") + -1.00 * t2["aaaa"]("c,b,l,i") * l2["aaaa"]("L,l,m,c,b")) * tmps_["103_aa_Loo"]("R,j,k");
        D_ooov_aaaa("L,R,i,j,k,a") -= tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k");
        D_ooov_aaaa("L,R,i,j,k,a") += tmps_["perm_aaaa_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_aaaa_LLvooo"].~TArrayD();

        // D_oooo_abab  = +0.50 d_aa(i,k) r0 t2_abab(b,a,m,j) l2_abab(m,l,b,a)
        //             += +0.50 d_aa(i,k) r0 t2_abab(a,b,m,j) l2_abab(m,l,a,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t1_bb(a,m) t2_1_abab(c,b,l,i) l2_1_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t1_bb(a,m) t2_1_abab(b,c,l,i) l2_1_abab(l,m,b,c)
        // flops: o4v0L2  = o6v0L2
        //  mems: o4v0L2  = o6v0L2
        D_oooo_abab("L,R,i,j,k,l")  = tmps_["60_bbbb_Loooo"]("L,j,l,i,m") * tmps_["103_aa_Loo"]("R,i,k");

        // flops: o5v0L2  = o5v0L2
        //  mems: o5v0L2  = o5v0L2
        tmps_["104_aaaaa_LLooooo"]("L,R,i,j,k,l,m")  = tmps_["49_aaa_Looo"]("L,i,j,k") * tmps_["103_aa_Loo"]("R,l,m");

        // D_oooo_aaaa += +0.50 P(i,j) d_aa(j,l) r0 t2_1_abab(b,a,i,m) l2_1_abab(k,m,b,a)
        //             += +0.50 P(i,j) d_aa(j,l) r0 t2_1_abab(a,b,i,m) l2_1_abab(k,m,a,b)
        //             += -0.50 P(i,j) d_aa(j,k) r1_aa(a,m) t2_abab(c,b,i,l) l2_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r1_aa(a,m) t2_abab(b,c,i,l) l2_abab(m,l,b,c)
        D_oooo_aaaa("L,R,i,j,k,l") += tmps_["104_aaaaa_LLooooo"]("L,R,i,k,m,j,l");
        D_oooo_aaaa("L,R,i,j,k,l") -= tmps_["104_aaaaa_LLooooo"]("L,R,j,k,m,i,l");

        // D_oooo_aaaa += -0.50 P(i,j) d_aa(j,k) r0 t2_1_abab(b,a,i,m) l2_1_abab(l,m,b,a)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_1_abab(a,b,i,m) l2_1_abab(l,m,a,b)
        //             += +0.50 P(i,j) d_aa(j,k) r1_aa(a,m) t2_abab(c,b,i,l) l2_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r1_aa(a,m) t2_abab(b,c,i,l) l2_abab(m,l,b,c)
        D_oooo_aaaa("L,R,i,j,k,l") -= tmps_["104_aaaaa_LLooooo"]("L,R,i,l,m,j,k");
        D_oooo_aaaa("L,R,i,j,k,l") += tmps_["104_aaaaa_LLooooo"]("L,R,j,l,m,i,k");
        tmps_["104_aaaaa_LLooooo"].~TArrayD();

        // flops: o3v0L1  = o3v2L1 o3v2L1 o3v0L1
        //  mems: o3v0L1  = o2v0L1 o2v0L1 o3v0L1
        tmps_["105_aaa_Looo"]("L,i,j,k")  = t2_1["aaaa"]("a,b,k,i") * l2_1["aaaa"]("L,k,j,a,b");
        tmps_["105_aaa_Looo"]("L,i,j,k") -= t2["aaaa"]("c,a,l,i") * l2["aaaa"]("L,l,k,c,a");

        // flops: o2v1L1  = o3v1L1
        //  mems: o2v1L1  = o2v1L1
        tmps_["106_aaa_Lvoo"]("L,a,i,j")  = t1["aa"]("a,k") * tmps_["105_aaa_Looo"]("L,i,k,j");

        // D_oovv_aaaa += +0.50 P(i,j) P(a,b) r1_aa(a,i) t1_aa(b,l) t2_1_aaaa(d,c,k,j) l2_1_aaaa(k,l,d,c)
        //             += -0.50 d_bb(i,k) r0 t1_aa(a,m) t2_aaaa(c,b,l,j) l2_aaaa(l,m,c,b)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o3v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["106_aaa_Lvoo"]("L,b,j,m") * r1["aa"]("R,a,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r1_bb(b,j) t1_aa(a,l) t2_1_aaaa(d,c,k,i) l2_1_aaaa(k,l,d,c)
        //             += -0.50 P(i,j) r0 t1_aa(a,j) t2_aaaa(c,b,l,i) l2_aaaa(l,k,c,b)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o3v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * tmps_["106_aaa_Lvoo"]("L,a,i,k") * r1["bb"]("R,b,j");

        // flops: o2v0L1  = o3v2L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["107_bb_Loo"]("L,i,j")  = t2["bbbb"]("a,b,k,i") * l2_1["bbbb"]("L,k,j,a,b");

        // D_oovv_abab += +0.50 r0_1 t2_abab(a,b,i,l) t2_bbbb(d,c,k,j) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * t2["abab"]("a,b,i,l") * tmps_["107_bb_Loo"]("L,j,l") * r0_1("R");

        // D_oovv_bbbb += -0.50 P(i,j) r0_1 t2_bbbb(b,a,l,j) t2_bbbb(d,c,k,i) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t2["bbbb"]("b,a,l,j") * tmps_["107_bb_Loo"]("L,i,l") * r0_1("R");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v0L1  = o2v1L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["108_bb_Loo"]("L,i,j")  = t1["bb"]("a,i") * l1_1["bb"]("L,j,a");

        // D_oovv_abab += +1.00 r2_1_abab(a,b,i,k) t1_bb(c,j) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["108_bb_Loo"]("L,j,k") * r2_1["abab"]("R,a,b,i,k");

        // D_oovv_abab += +1.00 r0_1 t2_abab(a,b,i,k) t1_bb(c,j) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,b,i,k") * tmps_["108_bb_Loo"]("L,j,k") * r0_1("R");

        // D_oovv_abab += +1.00 r0 t1_bb(c,j) t2_1_abab(a,b,i,k) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2_1["abab"]("a,b,i,k") * tmps_["108_bb_Loo"]("L,j,k") * r0("R");

        // D_oovv_bbbb += +1.00 P(i,j) r2_1_bbbb(b,a,k,i) t1_bb(c,j) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["108_bb_Loo"]("L,j,k") * r2_1["bbbb"]("R,b,a,k,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) r0 t1_bb(c,j) t2_1_bbbb(b,a,k,i) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2_1["bbbb"]("b,a,k,i") * tmps_["108_bb_Loo"]("L,j,k") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) r0_1 t2_bbbb(b,a,k,i) t1_bb(c,j) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["bbbb"]("b,a,k,i") * tmps_["108_bb_Loo"]("L,j,k") * r0_1("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v0L1  = o2v1L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["109_bb_Loo"]("L,i,j")  = t1["bb"]("a,i") * l1["bb"]("L,j,a");

        // D_oovv_abab += +1.00 r2_abab(a,b,i,k) t1_bb(c,j) l1_bb(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["109_bb_Loo"]("L,j,k") * r2["abab"]("R,a,b,i,k");

        // D_oovv_abab += +1.00 r0 t2_abab(a,b,i,k) t1_bb(c,j) l1_bb(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,b,i,k") * tmps_["109_bb_Loo"]("L,j,k") * r0("R");

        // D_oovv_bbbb += +1.00 P(i,j) r2_bbbb(b,a,k,i) t1_bb(c,j) l1_bb(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["109_bb_Loo"]("L,j,k") * r2["bbbb"]("R,b,a,k,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oooo_abab += +1.00 d_aa(i,k) r0 t1_bb(a,j) l1_bb(l,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += tmps_["109_bb_Loo"]("L,j,l") * tmps_["103_aa_Loo"]("R,i,k");

        // D_oovv_bbbb += +1.00 P(i,j) r0 t2_bbbb(b,a,k,i) t1_bb(c,j) l1_bb(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["bbbb"]("b,a,k,i") * tmps_["109_bb_Loo"]("L,j,k") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v0L1  = o2v1L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["110_bb_Loo"]("L,i,j")  = t1_1["bb"]("a,i") * l1_1["bb"]("L,j,a");

        // D_oovv_abab += +1.00 r0 t2_abab(a,b,i,k) t1_1_bb(c,j) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,b,i,k") * tmps_["110_bb_Loo"]("L,j,k") * r0("R");

        // D_oovv_bbbb += -1.00 P(i,j) r0 t2_bbbb(b,a,k,j) t1_1_bb(c,i) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["bbbb"]("b,a,k,j") * tmps_["110_bb_Loo"]("L,i,k") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oooo_abab += +1.00 d_aa(i,k) r0 t1_1_bb(a,j) l1_1_bb(l,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += tmps_["110_bb_Loo"]("L,j,l") * tmps_["103_aa_Loo"]("R,i,k");

        // D_oovv_abab += +1.00 r2_abab(a,b,i,k) t1_1_bb(c,j) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["110_bb_Loo"]("L,j,k") * r2["abab"]("R,a,b,i,k");

        // D_oovv_bbbb += +1.00 P(i,j) r2_bbbb(b,a,k,i) t1_1_bb(c,j) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["110_bb_Loo"]("L,j,k") * r2["bbbb"]("R,b,a,k,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v0L2  = o2v1L2 o3v2L2 o3v2L2 o2v0L2 o3v2L2 o2v0L2 o2v0L2 o3v2L2 o2v0L2 o2v1L2 o2v0L2
        //  mems: o2v0L2  = o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2
        tmps_["111_bb_LLoo"]("L,R,i,j")  = 2.00 * l1["bb"]("L,i,a") * r1["bb"]("R,a,j");
        tmps_["111_bb_LLoo"]("L,R,i,j") += 2.00 * l2["abab"]("L,k,i,b,c") * r2["abab"]("R,b,c,k,j");
        tmps_["111_bb_LLoo"]("L,R,i,j") += l2_1["bbbb"]("L,l,i,a,c") * r2_1["bbbb"]("R,a,c,l,j");
        tmps_["111_bb_LLoo"]("L,R,i,j") += l2["bbbb"]("L,l,i,a,c") * r2["bbbb"]("R,a,c,l,j");
        tmps_["111_bb_LLoo"]("L,R,i,j") += 2.00 * l2_1["abab"]("L,k,i,b,c") * r2_1["abab"]("R,b,c,k,j");
        tmps_["111_bb_LLoo"]("L,R,i,j") += 2.00 * l1_1["bb"]("L,i,a") * r1_1["bb"]("R,a,j");

        // D_oooo_abab += +0.50 d_aa(i,k) r2_1_bbbb(a,b,m,j) l2_1_bbbb(m,l,a,b)
        //             += +0.50 d_aa(i,k) r2_abab(a,b,m,j) l2_abab(m,l,a,b)
        //             += +0.50 d_aa(i,k) r2_abab(b,a,m,j) l2_abab(m,l,b,a)
        //             += +0.50 d_aa(i,k) r2_bbbb(a,b,m,j) l2_bbbb(m,l,a,b)
        //             += +1.00 d_aa(i,k) r1_bb(a,j) l1_bb(l,a)
        //             += +0.50 d_aa(i,k) r2_1_abab(a,b,m,j) l2_1_abab(m,l,a,b)
        //             += +0.50 d_aa(i,k) r2_1_abab(b,a,m,j) l2_1_abab(m,l,b,a)
        //             += +1.00 d_aa(i,k) r1_1_bb(a,j) l1_1_bb(l,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += 0.50 * Id["aa_oo"]("i,k") * tmps_["111_bb_LLoo"]("L,R,l,j");

        // D_oovv_abab += +0.50 r2_1_bbbb(c,d,l,j) t2_abab(a,b,i,k) l2_1_bbbb(l,k,c,d)
        //             += +0.50 r2_abab(c,d,l,j) t2_abab(a,b,i,k) l2_abab(l,k,c,d)
        //             += +0.50 r2_abab(d,c,l,j) t2_abab(a,b,i,k) l2_abab(l,k,d,c)
        //             += +0.50 r2_bbbb(c,d,l,j) t2_abab(a,b,i,k) l2_bbbb(l,k,c,d)
        //             += +1.00 r1_bb(c,j) t2_abab(a,b,i,k) l1_bb(k,c)
        //             += +0.50 r2_1_abab(c,d,l,j) t2_abab(a,b,i,k) l2_1_abab(l,k,c,d)
        //             += +0.50 r2_1_abab(d,c,l,j) t2_abab(a,b,i,k) l2_1_abab(l,k,d,c)
        //             += +1.00 r1_1_bb(c,j) t2_abab(a,b,i,k) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * t2["abab"]("a,b,i,k") * tmps_["111_bb_LLoo"]("L,R,k,j");

        // D_oovv_bbbb += -0.50 P(i,j) r2_1_bbbb(c,d,l,i) t2_bbbb(b,a,k,j) l2_1_bbbb(l,k,c,d)
        //             += -0.50 P(i,j) r2_abab(c,d,l,i) t2_bbbb(b,a,k,j) l2_abab(l,k,c,d)
        //             += -0.50 P(i,j) r2_abab(d,c,l,i) t2_bbbb(b,a,k,j) l2_abab(l,k,d,c)
        //             += -0.50 P(i,j) r2_bbbb(c,d,l,i) t2_bbbb(b,a,k,j) l2_bbbb(l,k,c,d)
        //             += -1.00 P(i,j) r1_bb(c,i) t2_bbbb(b,a,k,j) l1_bb(k,c)
        //             += -0.50 P(i,j) r2_1_abab(c,d,l,i) t2_bbbb(b,a,k,j) l2_1_abab(l,k,c,d)
        //             += -0.50 P(i,j) r2_1_abab(d,c,l,i) t2_bbbb(b,a,k,j) l2_1_abab(l,k,d,c)
        //             += -1.00 P(i,j) r1_1_bb(c,i) t2_bbbb(b,a,k,j) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t2["bbbb"]("b,a,k,j") * tmps_["111_bb_LLoo"]("L,R,k,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v0L1  = o2v0L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["112_bb_Loo"]("R,i,j")  = Id["bb_oo"]("i,j") * r0("R");

        // D_oooo_abab += +0.50 d_bb(j,l) r0 t2_1_abab(b,a,i,m) l2_1_abab(k,m,b,a)
        //             += +0.50 d_bb(j,l) r0 t2_1_abab(a,b,i,m) l2_1_abab(k,m,a,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t1_aa(a,m) t2_abab(c,b,i,l) l2_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t1_aa(a,m) t2_abab(b,c,i,l) l2_abab(m,l,b,c)
        // flops: o4v0L2 += o5v0L2
        //  mems: o4v0L2 += o5v0L2
        D_oooo_abab("L,R,i,j,k,l") += tmps_["49_aaa_Looo"]("L,i,k,m") * tmps_["112_bb_Loo"]("R,j,l");

        // D_oooo_abab += +0.50 d_bb(j,l) r0 t2_1_aaaa(b,a,m,i) l2_1_aaaa(m,k,b,a)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t1_aa(a,m) t2_aaaa(c,b,l,i) l2_aaaa(l,m,c,b)
        // flops: o4v0L2 += o5v0L2
        //  mems: o4v0L2 += o5v0L2
        D_oooo_abab("L,R,i,j,k,l") += 0.50 * tmps_["105_aaa_Looo"]("L,i,k,m") * tmps_["112_bb_Loo"]("R,j,l");

        // D_ooov_bbbb += -0.50 P(i,j) d_bb(j,k) r0 t1_bb(a,m) t2_abab(c,b,l,i) l2_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t1_bb(a,m) t2_abab(b,c,l,i) l2_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t1_bb(a,m) t2_1_abab(c,b,l,i) l2_1_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t1_bb(a,m) t2_1_abab(b,c,l,i) l2_1_abab(l,m,b,c)
        // flops: o3v1L2 += o3v2L1 o3v2L1 o2v0L1 o2v1L1 o3v1L2
        //  mems: o3v1L2 += o2v0L1 o2v0L1 o2v0L1 o1v1L1 o3v1L2
        tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k")  = t1["bb"]("a,m") * (t2["abab"]("c,b,l,i") * l2["abab"]("L,l,m,c,b") + -1.00 * t2_1["abab"]("c,b,l,i") * l2_1["abab"]("L,l,m,c,b")) * tmps_["112_bb_Loo"]("R,j,k");
        D_ooov_bbbb("L,R,i,j,k,a") -= tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k");
        D_ooov_bbbb("L,R,i,j,k,a") += tmps_["perm_bbbb_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_bbbb_LLvooo"].~TArrayD();

        // D_oovo_abab += -0.50 d_bb(i,k) r0 t1_aa(a,m) t2_1_aaaa(c,b,l,j) l2_1_aaaa(l,m,c,b)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,i) t1_aa(b,l) t2_aaaa(d,c,k,j) l2_aaaa(k,l,d,c)
        // flops: o3v1L2 += o4v1L2
        //  mems: o3v1L2 += o4v1L2
        D_oovo_abab("L,R,i,j,a,k") -= 0.50 * tmps_["106_aaa_Lvoo"]("L,a,j,l") * tmps_["112_bb_Loo"]("R,i,k");
        tmps_["106_aaa_Lvoo"].~TArrayD();

        // flops: o2v0L1  = o2v0L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["113_bb_Loo"]("R,i,j")  = Id["bb_oo"]("i,j") * r0_1("R");

        // D_oooo_abab += +0.50 d_bb(j,l) r0_1 t2_abab(b,a,i,m) l2_1_abab(k,m,b,a)
        //             += +0.50 d_bb(j,l) r0_1 t2_abab(a,b,i,m) l2_1_abab(k,m,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += tmps_["63_aa_Loo"]("L,i,k") * tmps_["113_bb_Loo"]("R,j,l");

        // flops: o2v0L2  = o1v1L2 o2v1L2
        //  mems: o2v0L2  = o1v1L2 o2v0L2
        tmps_["114_bb_LLoo"]("L,R,i,j")  = (tmps_["72_bb_LLov"]("L,R,i,a") + -1.00 * tmps_["69_bb_LLov"]("L,R,i,a")) * t1_1["bb"]("a,j");

        // D_oooo_abab += -1.00 d_aa(i,k) r1_bb(b,m) t1_1_bb(a,j) l2_1_bbbb(m,l,a,b)
        //             += +1.00 d_aa(i,k) r1_aa(b,m) t1_1_bb(a,j) l2_1_abab(m,l,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") -= Id["aa_oo"]("i,k") * tmps_["114_bb_LLoo"]("L,R,l,j");

        // D_oovv_abab += -1.00 r1_bb(d,l) t2_abab(a,b,i,k) t1_1_bb(c,j) l2_1_bbbb(l,k,c,d)
        //             += +1.00 r1_aa(d,l) t2_abab(a,b,i,k) t1_1_bb(c,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,b,i,k") * tmps_["114_bb_LLoo"]("L,R,k,j");

        // D_oovv_bbbb += +1.00 P(i,j) r1_bb(d,l) t2_bbbb(b,a,k,j) t1_1_bb(c,i) l2_1_bbbb(l,k,c,d)
        //             += -1.00 P(i,j) r1_aa(d,l) t2_bbbb(b,a,k,j) t1_1_bb(c,i) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["bbbb"]("b,a,k,j") * tmps_["114_bb_LLoo"]("L,R,k,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v0L2  = o1v1L2 o1v1L2 o1v1L2 o2v1L2
        //  mems: o2v0L2  = o1v1L2 o1v1L2 o1v1L2 o2v0L2
        tmps_["115_bb_LLoo"]("L,R,i,j")  = (-1.00 * tmps_["70_bb_LLov"]("L,R,i,a") + tmps_["74_bb_LLov"]("L,R,i,a") + -1.00 * tmps_["71_bb_LLov"]("L,R,i,a") + tmps_["73_bb_LLov"]("L,R,i,a")) * t1["bb"]("a,j");

        // D_oooo_abab += -1.00 d_aa(i,k) r1_bb(b,m) t1_bb(a,j) l2_bbbb(m,l,a,b)
        //             += +1.00 d_aa(i,k) r1_aa(b,m) t1_bb(a,j) l2_abab(m,l,b,a)
        //             += +1.00 d_aa(i,k) r1_1_aa(b,m) t1_bb(a,j) l2_1_abab(m,l,b,a)
        //             += -1.00 d_aa(i,k) r1_1_bb(b,m) t1_bb(a,j) l2_1_bbbb(m,l,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") -= Id["aa_oo"]("i,k") * tmps_["115_bb_LLoo"]("L,R,l,j");

        // D_oovv_abab += -1.00 r1_bb(d,l) t2_abab(a,b,i,k) t1_bb(c,j) l2_bbbb(l,k,c,d)
        //             += +1.00 r1_aa(d,l) t2_abab(a,b,i,k) t1_bb(c,j) l2_abab(l,k,d,c)
        //             += +1.00 r1_1_aa(d,l) t2_abab(a,b,i,k) t1_bb(c,j) l2_1_abab(l,k,d,c)
        //             += -1.00 r1_1_bb(d,l) t2_abab(a,b,i,k) t1_bb(c,j) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,b,i,k") * tmps_["115_bb_LLoo"]("L,R,k,j");

        // D_oovv_bbbb += -1.00 P(i,j) r1_bb(d,l) t2_bbbb(b,a,k,i) t1_bb(c,j) l2_bbbb(l,k,c,d)
        //             += +1.00 P(i,j) r1_aa(d,l) t2_bbbb(b,a,k,i) t1_bb(c,j) l2_abab(l,k,d,c)
        //             += +1.00 P(i,j) r1_1_aa(d,l) t2_bbbb(b,a,k,i) t1_bb(c,j) l2_1_abab(l,k,d,c)
        //             += -1.00 P(i,j) r1_1_bb(d,l) t2_bbbb(b,a,k,i) t1_bb(c,j) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["bbbb"]("b,a,k,i") * tmps_["115_bb_LLoo"]("L,R,k,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o4v0L2  = o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2
        //  mems: o4v0L2  = o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2
        tmps_["116_bbbb_LLoooo"]("L,R,i,j,k,l")  = 2.00 * Id["bb_oo"]("i,k") * tmps_["114_bb_LLoo"]("L,R,j,l");
        tmps_["116_bbbb_LLoooo"]("L,R,i,j,k,l") += 2.00 * Id["bb_oo"]("i,k") * tmps_["115_bb_LLoo"]("L,R,j,l");
        tmps_["116_bbbb_LLoooo"]("L,R,i,j,k,l") += Id["bb_oo"]("i,j") * tmps_["111_bb_LLoo"]("L,R,k,l");
        tmps_["116_bbbb_LLoooo"]("L,R,i,j,k,l") += 2.00 * tmps_["96_bb_Loo"]("L,l,k") * tmps_["113_bb_Loo"]("R,i,j");
        tmps_["116_bbbb_LLoooo"]("L,R,i,j,k,l") += 2.00 * tmps_["110_bb_Loo"]("L,l,k") * tmps_["112_bb_Loo"]("R,i,j");
        tmps_["116_bbbb_LLoooo"]("L,R,i,j,k,l") += tmps_["107_bb_Loo"]("L,l,k") * tmps_["113_bb_Loo"]("R,i,j");
        tmps_["116_bbbb_LLoooo"]("L,R,i,j,k,l") += 2.00 * tmps_["108_bb_Loo"]("L,l,k") * tmps_["113_bb_Loo"]("R,i,j");
        tmps_["116_bbbb_LLoooo"]("L,R,i,j,k,l") += 2.00 * tmps_["109_bb_Loo"]("L,l,k") * tmps_["112_bb_Loo"]("R,i,j");

        // D_oooo_bbbb += -0.50 P(i,j) d_bb(j,k) r2_1_bbbb(a,b,m,i) l2_1_bbbb(m,l,a,b)
        //             += -0.50 P(i,j) d_bb(j,k) r2_abab(a,b,m,i) l2_abab(m,l,a,b)
        //             += -0.50 P(i,j) d_bb(j,k) r2_abab(b,a,m,i) l2_abab(m,l,b,a)
        //             += -0.50 P(i,j) d_bb(j,k) r2_bbbb(a,b,m,i) l2_bbbb(m,l,a,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(a,i) l1_bb(l,a)
        //             += -0.50 P(i,j) d_bb(j,k) r2_1_abab(a,b,m,i) l2_1_abab(m,l,a,b)
        //             += -0.50 P(i,j) d_bb(j,k) r2_1_abab(b,a,m,i) l2_1_abab(m,l,b,a)
        //             += -1.00 P(i,j) d_bb(j,k) r1_1_bb(a,i) l1_1_bb(l,a)
        //             += -1.00 P(i,j) d_bb(j,l) r1_bb(b,m) t1_bb(a,i) l2_bbbb(m,k,a,b)
        //             += +1.00 P(i,j) d_bb(j,l) r1_aa(b,m) t1_bb(a,i) l2_abab(m,k,b,a)
        //             += +1.00 P(i,j) d_bb(j,l) r1_1_aa(b,m) t1_bb(a,i) l2_1_abab(m,k,b,a)
        //             += -1.00 P(i,j) d_bb(j,l) r1_1_bb(b,m) t1_bb(a,i) l2_1_bbbb(m,k,a,b)
        //             += -1.00 P(i,j) d_bb(j,l) r1_bb(b,m) t1_1_bb(a,i) l2_1_bbbb(m,k,a,b)
        //             += +1.00 P(i,j) d_bb(j,l) r1_aa(b,m) t1_1_bb(a,i) l2_1_abab(m,k,b,a)
        //             += -0.50 P(i,j) d_bb(j,k) r0_1 t2_abab(b,a,m,i) l2_1_abab(m,l,b,a)
        //             += -0.50 P(i,j) d_bb(j,k) r0_1 t2_abab(a,b,m,i) l2_1_abab(m,l,a,b)
        //             += -1.00 P(i,j) d_bb(j,k) r0 t1_1_bb(a,i) l1_1_bb(l,a)
        //             += -0.50 P(i,j) d_bb(j,k) r0_1 t2_bbbb(b,a,m,i) l2_1_bbbb(m,l,b,a)
        //             += -1.00 P(i,j) d_bb(j,k) r0_1 t1_bb(a,i) l1_1_bb(l,a)
        //             += -1.00 P(i,j) d_bb(j,k) r0 t1_bb(a,i) l1_bb(l,a)
        D_oooo_bbbb("L,R,i,j,k,l") -= 0.50 * tmps_["116_bbbb_LLoooo"]("L,R,j,k,l,i");
        D_oooo_bbbb("L,R,i,j,k,l") += 0.50 * tmps_["116_bbbb_LLoooo"]("L,R,i,k,l,j");

        // D_oooo_bbbb += +0.50 P(i,j) d_bb(j,l) r2_1_bbbb(a,b,m,i) l2_1_bbbb(m,k,a,b)
        //             += +0.50 P(i,j) d_bb(j,l) r2_abab(a,b,m,i) l2_abab(m,k,a,b)
        //             += +0.50 P(i,j) d_bb(j,l) r2_abab(b,a,m,i) l2_abab(m,k,b,a)
        //             += +0.50 P(i,j) d_bb(j,l) r2_bbbb(a,b,m,i) l2_bbbb(m,k,a,b)
        //             += +1.00 P(i,j) d_bb(j,l) r1_bb(a,i) l1_bb(k,a)
        //             += +0.50 P(i,j) d_bb(j,l) r2_1_abab(a,b,m,i) l2_1_abab(m,k,a,b)
        //             += +0.50 P(i,j) d_bb(j,l) r2_1_abab(b,a,m,i) l2_1_abab(m,k,b,a)
        //             += +1.00 P(i,j) d_bb(j,l) r1_1_bb(a,i) l1_1_bb(k,a)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(b,m) t1_bb(a,i) l2_bbbb(m,l,a,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_aa(b,m) t1_bb(a,i) l2_abab(m,l,b,a)
        //             += -1.00 P(i,j) d_bb(j,k) r1_1_aa(b,m) t1_bb(a,i) l2_1_abab(m,l,b,a)
        //             += +1.00 P(i,j) d_bb(j,k) r1_1_bb(b,m) t1_bb(a,i) l2_1_bbbb(m,l,a,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(b,m) t1_1_bb(a,i) l2_1_bbbb(m,l,a,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_aa(b,m) t1_1_bb(a,i) l2_1_abab(m,l,b,a)
        //             += +0.50 P(i,j) d_bb(j,l) r0_1 t2_abab(b,a,m,i) l2_1_abab(m,k,b,a)
        //             += +0.50 P(i,j) d_bb(j,l) r0_1 t2_abab(a,b,m,i) l2_1_abab(m,k,a,b)
        //             += +1.00 P(i,j) d_bb(j,l) r0 t1_1_bb(a,i) l1_1_bb(k,a)
        //             += +0.50 P(i,j) d_bb(j,l) r0_1 t2_bbbb(b,a,m,i) l2_1_bbbb(m,k,b,a)
        //             += +1.00 P(i,j) d_bb(j,l) r0_1 t1_bb(a,i) l1_1_bb(k,a)
        //             += +1.00 P(i,j) d_bb(j,l) r0 t1_bb(a,i) l1_bb(k,a)
        D_oooo_bbbb("L,R,i,j,k,l") += 0.50 * tmps_["116_bbbb_LLoooo"]("L,R,j,l,k,i");
        D_oooo_bbbb("L,R,i,j,k,l") -= 0.50 * tmps_["116_bbbb_LLoooo"]("L,R,i,l,k,j");
        tmps_["116_bbbb_LLoooo"].~TArrayD();

        // flops: o1v1L2  = o2v2L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["117_aa_LLov"]("L,R,i,a")  = l2_1["aaaa"]("L,j,i,a,b") * r1["aa"]("R,b,j");

        // D_oovv_aaaa += +1.00 P(a,b) r1_aa(d,l) t2_aaaa(b,c,i,j) t1_1_aa(a,k) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["aaaa"]("b,c,i,j") * tmps_["117_aa_LLov"]("L,R,k,c") * t1_1["aa"]("a,k");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r1_aa(d,l) t2_abab(c,b,i,j) t1_1_aa(a,k) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("c,b,i,j") * tmps_["117_aa_LLov"]("L,R,k,c") * t1_1["aa"]("a,k");

        // D_oovv_abab += -1.00 r1_aa(d,l) t1_aa(a,k) t2_1_abab(c,b,i,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2_1["abab"]("c,b,i,j") * tmps_["117_aa_LLov"]("L,R,k,c") * t1["aa"]("a,k");

        // flops: o1v1L2  = o2v2L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["118_aa_LLov"]("L,R,i,a")  = l2_1["abab"]("L,i,j,a,b") * r1["bb"]("R,b,j");

        // D_oovv_aaaa += -1.00 P(a,b) r1_bb(d,l) t2_aaaa(b,c,i,j) t1_1_aa(a,k) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["aaaa"]("b,c,i,j") * tmps_["118_aa_LLov"]("L,R,k,c") * t1_1["aa"]("a,k");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +1.00 r1_bb(d,l) t1_aa(a,k) t2_1_abab(c,b,i,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2_1["abab"]("c,b,i,j") * tmps_["118_aa_LLov"]("L,R,k,c") * t1["aa"]("a,k");

        // D_oovv_abab += +1.00 r1_bb(d,l) t2_abab(c,b,i,j) t1_1_aa(a,k) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("c,b,i,j") * tmps_["118_aa_LLov"]("L,R,k,c") * t1_1["aa"]("a,k");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["119_aaaa_Lvoov"]("L,a,i,j,b")  = t2["aaaa"]("a,c,k,i") * l2_1["aaaa"]("L,j,k,c,b");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r1_aa(a,l) t2_aaaa(b,d,k,j) t1_1_aa(c,i) l2_1_aaaa(l,k,d,c)
        // flops: o2v2L2 += o3v2L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["119_aaaa_Lvoov"]("L,b,j,l,c") * t1_1["aa"]("c,i") * r1["aa"]("R,a,l");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r1_aa(d,i) t2_aaaa(b,c,l,j) t1_1_aa(a,k) l2_1_aaaa(k,l,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["119_aaaa_Lvoov"]("L,b,j,k,d") * r1["aa"]("R,d,i") * t1_1["aa"]("a,k");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r2_1_abab(d,b,l,j) t2_aaaa(a,c,k,i) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["119_aaaa_Lvoov"]("L,a,i,l,d") * r2_1["abab"]("R,d,b,l,j");

        // D_oovv_abab += -1.00 r0 t2_aaaa(a,d,l,i) t2_1_abab(c,b,k,j) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["119_aaaa_Lvoov"]("L,a,i,k,c") * t2_1["abab"]("c,b,k,j") * r0("R");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["120_aabb_Lvoov"]("L,a,i,j,b")  = t2["abab"]("a,c,i,k") * l2_1["bbbb"]("L,j,k,c,b");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_1_abab(a,d,i,l) t2_abab(b,c,j,k) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["120_aabb_Lvoov"]("L,b,j,l,d") * r2_1["abab"]("R,a,d,i,l");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) r0_1 t2_abab(a,c,i,k) t2_abab(b,d,j,l) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["abab"]("a,c,i,k") * tmps_["120_aabb_Lvoov"]("L,b,j,k,c") * r0_1("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r0 t2_abab(b,d,j,l) t2_1_abab(a,c,i,k) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2_1["abab"]("a,c,i,k") * tmps_["120_aabb_Lvoov"]("L,b,j,k,c") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r0 t2_abab(a,d,i,l) t2_1_bbbb(b,c,k,j) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["120_aabb_Lvoov"]("L,a,i,k,c") * t2_1["bbbb"]("b,c,k,j") * r0("R");

        // D_oovv_abab += -1.00 r1_bb(b,l) t2_abab(a,d,i,k) t1_1_bb(c,j) l2_1_bbbb(l,k,d,c)
        // flops: o2v2L2 += o3v2L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["120_aabb_Lvoov"]("L,a,i,l,c") * t1_1["bb"]("c,j") * r1["bb"]("R,b,l");

        // D_oovv_abab += -1.00 r1_bb(d,j) t2_abab(a,c,i,l) t1_1_bb(b,k) l2_1_bbbb(k,l,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["120_aabb_Lvoov"]("L,a,i,k,d") * r1["bb"]("R,d,j") * t1_1["bb"]("b,k");

        // flops: o2v0L2  = o1v1L2 o2v1L2
        //  mems: o2v0L2  = o1v1L2 o2v0L2
        tmps_["121_aa_LLoo"]("L,R,i,j")  = (tmps_["118_aa_LLov"]("L,R,i,a") + -1.00 * tmps_["117_aa_LLov"]("L,R,i,a")) * t1["aa"]("a,j");

        // D_oovv_aaaa += +1.00 P(i,j) r1_bb(d,l) t1_aa(c,j) t2_1_aaaa(b,a,k,i) l2_1_abab(k,l,c,d)
        //             += -1.00 P(i,j) r1_aa(d,l) t1_aa(c,j) t2_1_aaaa(b,a,k,i) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2_1["aaaa"]("b,a,k,i") * tmps_["121_aa_LLoo"]("L,R,k,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +1.00 r1_bb(d,l) t1_aa(c,i) t2_1_abab(a,b,k,j) l2_1_abab(k,l,c,d)
        //             += -1.00 r1_aa(d,l) t1_aa(c,i) t2_1_abab(a,b,k,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2_1["abab"]("a,b,k,j") * tmps_["121_aa_LLoo"]("L,R,k,i");

        // flops: o1v1L2  = o2v2L2 o2v2L2 o1v1L2 o2v2L2 o1v1L2 o2v2L2 o1v1L2 o2v1L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2
        tmps_["122_aa_LLvo"]("L,R,a,i")  = -1.00 * t2["abab"]("a,e,i,m") * tmps_["69_bb_LLov"]("L,R,m,e");
        tmps_["122_aa_LLvo"]("L,R,a,i") += tmps_["16_aabb_Lvoov"]("L,a,i,l,d") * r1["bb"]("R,d,l");
        tmps_["122_aa_LLvo"]("L,R,a,i") += tmps_["120_aabb_Lvoov"]("L,a,i,l,d") * r1["bb"]("R,d,l");
        tmps_["122_aa_LLvo"]("L,R,a,i") -= tmps_["119_aaaa_Lvoov"]("L,a,i,k,c") * r1["aa"]("R,c,k");
        tmps_["122_aa_LLvo"]("L,R,a,i") += t1["aa"]("a,j") * tmps_["121_aa_LLoo"]("L,R,j,i");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r1_bb(d,l) t1_aa(b,k) t1_aa(c,j) t1_1_aa(a,i) l2_1_abab(k,l,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_aa(d,l) t1_aa(b,k) t1_aa(c,j) t1_1_aa(a,i) l2_1_aaaa(l,k,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_aa(d,l) t2_abab(b,c,j,k) t1_1_aa(a,i) l2_1_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r1_bb(d,l) t2_aaaa(b,c,k,j) t1_1_aa(a,i) l2_1_abab(k,l,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_bb(d,l) t2_abab(b,c,j,k) t1_1_aa(a,i) l2_1_bbbb(l,k,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_aa(d,l) t2_aaaa(b,c,k,j) t1_1_aa(a,i) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1_1["aa"]("a,i") * tmps_["122_aa_LLvo"]("L,R,b,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +1.00 r1_bb(d,l) t1_aa(a,k) t1_aa(c,i) t1_1_bb(b,j) l2_1_abab(k,l,c,d)
        //             += -1.00 r1_aa(d,l) t1_aa(a,k) t1_aa(c,i) t1_1_bb(b,j) l2_1_aaaa(l,k,c,d)
        //             += -1.00 r1_aa(d,l) t2_abab(a,c,i,k) t1_1_bb(b,j) l2_1_abab(l,k,d,c)
        //             += +1.00 r1_bb(d,l) t2_aaaa(a,c,k,i) t1_1_bb(b,j) l2_1_abab(k,l,c,d)
        //             += +1.00 r1_bb(d,l) t2_abab(a,c,i,k) t1_1_bb(b,j) l2_1_bbbb(l,k,c,d)
        //             += -1.00 r1_aa(d,l) t2_aaaa(a,c,k,i) t1_1_bb(b,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["122_aa_LLvo"]("L,R,a,i") * t1_1["bb"]("b,j");
        tmps_["122_aa_LLvo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["123_bbaa_Lvoov"]("L,a,i,j,b")  = t2["abab"]("c,a,k,i") * l2_1["aaaa"]("L,j,k,c,b");

        // D_oovv_abab += -1.00 r1_aa(a,l) t2_abab(d,b,k,j) t1_1_aa(c,i) l2_1_aaaa(l,k,d,c)
        // flops: o2v2L2 += o3v2L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["123_bbaa_Lvoov"]("L,b,j,l,c") * t1_1["aa"]("c,i") * r1["aa"]("R,a,l");

        // D_oovv_abab += -1.00 r1_aa(d,i) t2_abab(c,b,l,j) t1_1_aa(a,k) l2_1_aaaa(k,l,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["123_bbaa_Lvoov"]("L,b,j,k,d") * r1["aa"]("R,d,i") * t1_1["aa"]("a,k");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r2_1_abab(d,a,l,i) t2_abab(c,b,k,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["123_bbaa_Lvoov"]("L,b,j,l,d") * r2_1["abab"]("R,d,a,l,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) r0_1 t2_abab(c,a,k,i) t2_abab(d,b,l,j) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["abab"]("c,a,k,i") * tmps_["123_bbaa_Lvoov"]("L,b,j,k,c") * r0_1("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r0 t2_abab(d,b,l,j) t2_1_abab(c,a,k,i) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2_1["abab"]("c,a,k,i") * tmps_["123_bbaa_Lvoov"]("L,b,j,k,c") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["124_bbbb_Lvoov"]("L,a,i,j,b")  = t2["bbbb"]("a,c,k,i") * l2_1["bbbb"]("L,j,k,c,b");

        // D_oovv_abab += -1.00 r2_1_abab(a,d,i,l) t2_bbbb(b,c,k,j) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["124_bbbb_Lvoov"]("L,b,j,l,d") * r2_1["abab"]("R,a,d,i,l");

        // D_oovv_abab += -1.00 r0_1 t2_abab(a,c,i,k) t2_bbbb(b,d,l,j) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,c,i,k") * tmps_["124_bbbb_Lvoov"]("L,b,j,k,c") * r0_1("R");

        // D_oovv_abab += -1.00 r0 t2_bbbb(b,d,l,j) t2_1_abab(a,c,i,k) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2_1["abab"]("a,c,i,k") * tmps_["124_bbbb_Lvoov"]("L,b,j,k,c") * r0("R");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r0 t2_bbbb(b,d,l,j) t2_1_bbbb(a,c,k,i) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2_1["bbbb"]("a,c,k,i") * tmps_["124_bbbb_Lvoov"]("L,b,j,k,c") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) r0_1 t2_bbbb(a,c,k,i) t2_bbbb(b,d,l,j) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["bbbb"]("a,c,k,i") * tmps_["124_bbbb_Lvoov"]("L,b,j,k,c") * r0_1("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r1_bb(a,l) t2_bbbb(b,d,k,j) t1_1_bb(c,i) l2_1_bbbb(l,k,d,c)
        // flops: o2v2L2 += o3v2L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["124_bbbb_Lvoov"]("L,b,j,l,c") * t1_1["bb"]("c,i") * r1["bb"]("R,a,l");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r1_bb(d,i) t2_bbbb(b,c,l,j) t1_1_bb(a,k) l2_1_bbbb(k,l,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["124_bbbb_Lvoov"]("L,b,j,k,d") * r1["bb"]("R,d,i") * t1_1["bb"]("a,k");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v0L2  = o1v1L2 o2v1L2
        //  mems: o2v0L2  = o1v1L2 o2v0L2
        tmps_["125_bb_LLoo"]("L,R,i,j")  = (tmps_["69_bb_LLov"]("L,R,i,a") + -1.00 * tmps_["72_bb_LLov"]("L,R,i,a")) * t1["bb"]("a,j");

        // D_oovv_abab += +1.00 r1_aa(d,l) t1_bb(c,j) t2_1_abab(a,b,i,k) l2_1_abab(l,k,d,c)
        //             += -1.00 r1_bb(d,l) t1_bb(c,j) t2_1_abab(a,b,i,k) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2_1["abab"]("a,b,i,k") * tmps_["125_bb_LLoo"]("L,R,k,j");

        // D_oovv_bbbb += +1.00 P(i,j) r1_aa(d,l) t1_bb(c,j) t2_1_bbbb(b,a,k,i) l2_1_abab(l,k,d,c)
        //             += -1.00 P(i,j) r1_bb(d,l) t1_bb(c,j) t2_1_bbbb(b,a,k,i) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2_1["bbbb"]("b,a,k,i") * tmps_["125_bb_LLoo"]("L,R,k,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v1L2  = o2v2L2 o2v2L2 o1v1L2 o2v2L2 o1v1L2 o2v2L2 o1v1L2 o2v1L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2
        tmps_["126_bb_LLvo"]("L,R,a,i")  = -1.00 * tmps_["81_bbbb_Lvoov"]("L,a,i,k,c") * r1["bb"]("R,c,k");
        tmps_["126_bb_LLvo"]("L,R,a,i") += tmps_["93_bbaa_Lvoov"]("L,a,i,l,d") * r1["aa"]("R,d,l");
        tmps_["126_bb_LLvo"]("L,R,a,i") += tmps_["123_bbaa_Lvoov"]("L,a,i,l,d") * r1["aa"]("R,d,l");
        tmps_["126_bb_LLvo"]("L,R,a,i") -= tmps_["124_bbbb_Lvoov"]("L,a,i,k,c") * r1["bb"]("R,c,k");
        tmps_["126_bb_LLvo"]("L,R,a,i") += t1["bb"]("a,j") * tmps_["125_bb_LLoo"]("L,R,j,i");

        // D_oovv_abab += +1.00 r1_aa(d,l) t1_bb(b,k) t1_bb(c,j) t1_1_aa(a,i) l2_1_abab(l,k,d,c)
        //             += -1.00 r1_bb(d,l) t1_bb(b,k) t1_bb(c,j) t1_1_aa(a,i) l2_1_bbbb(l,k,c,d)
        //             += -1.00 r1_bb(d,l) t2_abab(c,b,k,j) t1_1_aa(a,i) l2_1_abab(k,l,c,d)
        //             += +1.00 r1_aa(d,l) t2_bbbb(b,c,k,j) t1_1_aa(a,i) l2_1_abab(l,k,d,c)
        //             += +1.00 r1_aa(d,l) t2_abab(c,b,k,j) t1_1_aa(a,i) l2_1_aaaa(l,k,c,d)
        //             += -1.00 r1_bb(d,l) t2_bbbb(b,c,k,j) t1_1_aa(a,i) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1_1["aa"]("a,i") * tmps_["126_bb_LLvo"]("L,R,b,j");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r1_aa(d,l) t1_bb(b,k) t1_bb(c,j) t1_1_bb(a,i) l2_1_abab(l,k,d,c)
        //             += -1.00 P(i,j) P(a,b) r1_bb(d,l) t1_bb(b,k) t1_bb(c,j) t1_1_bb(a,i) l2_1_bbbb(l,k,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_bb(d,l) t2_abab(c,b,k,j) t1_1_bb(a,i) l2_1_abab(k,l,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_aa(d,l) t2_bbbb(b,c,k,j) t1_1_bb(a,i) l2_1_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r1_aa(d,l) t2_abab(c,b,k,j) t1_1_bb(a,i) l2_1_aaaa(l,k,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_bb(d,l) t2_bbbb(b,c,k,j) t1_1_bb(a,i) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1_1["bb"]("a,i") * tmps_["126_bb_LLvo"]("L,R,b,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
        tmps_["126_bb_LLvo"].~TArrayD();

        // flops: o1v1L2  = o2v2L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["127_aa_LLov"]("L,R,i,a")  = l2["aaaa"]("L,j,i,a,b") * r1["aa"]("R,b,j");

        // D_oovv_abab += -1.00 r1_aa(d,l) t1_aa(a,k) t2_abab(c,b,i,j) l2_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("c,b,i,j") * tmps_["127_aa_LLov"]("L,R,k,c") * t1["aa"]("a,k");

        // flops: o1v1L2  = o2v2L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["128_aa_LLov"]("L,R,i,a")  = l2_1["aaaa"]("L,j,i,a,b") * r1_1["aa"]("R,b,j");

        // D_oovv_abab += -1.00 r1_1_aa(d,l) t1_aa(a,k) t2_abab(c,b,i,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("c,b,i,j") * tmps_["128_aa_LLov"]("L,R,k,c") * t1["aa"]("a,k");

        // flops: o1v1L2  = o2v2L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["129_aa_LLov"]("L,R,i,a")  = l2["abab"]("L,i,j,a,b") * r1["bb"]("R,b,j");

        // D_oovv_abab += +1.00 r1_bb(d,l) t1_aa(a,k) t2_abab(c,b,i,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("c,b,i,j") * tmps_["129_aa_LLov"]("L,R,k,c") * t1["aa"]("a,k");

        // flops: o1v1L2  = o2v2L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["130_aa_LLov"]("L,R,i,a")  = l2_1["abab"]("L,i,j,a,b") * r1_1["bb"]("R,b,j");

        // D_oovv_abab += +1.00 r1_1_bb(d,l) t1_aa(a,k) t2_abab(c,b,i,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("c,b,i,j") * tmps_["130_aa_LLov"]("L,R,k,c") * t1["aa"]("a,k");

        // flops: o2v0L1  = o3v2L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["131_aa_Loo"]("L,i,j")  = t2["aaaa"]("a,b,k,i") * l2_1["aaaa"]("L,k,j,a,b");

        // D_oovv_aaaa += -0.50 P(i,j) r0_1 t2_aaaa(b,a,l,j) t2_aaaa(d,c,k,i) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t2["aaaa"]("b,a,l,j") * tmps_["131_aa_Loo"]("L,i,l") * r0_1("R");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r0_1 t2_abab(a,b,l,j) t2_aaaa(d,c,k,i) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * t2["abab"]("a,b,l,j") * tmps_["131_aa_Loo"]("L,i,l") * r0_1("R");

        // D_oooo_abab += +0.50 d_bb(j,l) r0_1 t2_aaaa(b,a,m,i) l2_1_aaaa(m,k,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += 0.50 * tmps_["131_aa_Loo"]("L,i,k") * tmps_["113_bb_Loo"]("R,j,l");

        // flops: o2v0L1  = o2v1L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["132_aa_Loo"]("L,i,j")  = t1["aa"]("a,i") * l1_1["aa"]("L,j,a");

        // D_oovv_aaaa += +1.00 P(i,j) r2_1_aaaa(b,a,k,i) t1_aa(c,j) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["132_aa_Loo"]("L,j,k") * r2_1["aaaa"]("R,b,a,k,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) r0_1 t2_aaaa(b,a,k,i) t1_aa(c,j) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["aaaa"]("b,a,k,i") * tmps_["132_aa_Loo"]("L,j,k") * r0_1("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) r0 t1_aa(c,j) t2_1_aaaa(b,a,k,i) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2_1["aaaa"]("b,a,k,i") * tmps_["132_aa_Loo"]("L,j,k") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oooo_abab += +1.00 d_bb(j,l) r0_1 t1_aa(a,i) l1_1_aa(k,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += tmps_["132_aa_Loo"]("L,i,k") * tmps_["113_bb_Loo"]("R,j,l");

        // D_oovv_abab += +1.00 r2_1_abab(a,b,k,j) t1_aa(c,i) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["132_aa_Loo"]("L,i,k") * r2_1["abab"]("R,a,b,k,j");

        // D_oovv_abab += +1.00 r0 t1_aa(c,i) t2_1_abab(a,b,k,j) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2_1["abab"]("a,b,k,j") * tmps_["132_aa_Loo"]("L,i,k") * r0("R");

        // D_oovv_abab += +1.00 r0_1 t2_abab(a,b,k,j) t1_aa(c,i) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,b,k,j") * tmps_["132_aa_Loo"]("L,i,k") * r0_1("R");

        // flops: o2v0L1  = o2v1L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["133_aa_Loo"]("L,i,j")  = t1["aa"]("a,i") * l1["aa"]("L,j,a");

        // D_oovv_aaaa += +1.00 P(i,j) r2_aaaa(b,a,k,i) t1_aa(c,j) l1_aa(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["133_aa_Loo"]("L,j,k") * r2["aaaa"]("R,b,a,k,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) r0 t2_aaaa(b,a,k,i) t1_aa(c,j) l1_aa(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["aaaa"]("b,a,k,i") * tmps_["133_aa_Loo"]("L,j,k") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oooo_abab += +1.00 d_bb(j,l) r0 t1_aa(a,i) l1_aa(k,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += tmps_["133_aa_Loo"]("L,i,k") * tmps_["112_bb_Loo"]("R,j,l");

        // D_oovv_abab += +1.00 r2_abab(a,b,k,j) t1_aa(c,i) l1_aa(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["133_aa_Loo"]("L,i,k") * r2["abab"]("R,a,b,k,j");

        // D_oovv_abab += +1.00 r0 t2_abab(a,b,k,j) t1_aa(c,i) l1_aa(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,b,k,j") * tmps_["133_aa_Loo"]("L,i,k") * r0("R");

        // flops: o2v0L1  = o2v1L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["134_aa_Loo"]("L,i,j")  = t1_1["aa"]("a,i") * l1_1["aa"]("L,j,a");

        // D_oovv_aaaa += +1.00 P(i,j) r2_aaaa(b,a,k,i) t1_1_aa(c,j) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["134_aa_Loo"]("L,j,k") * r2["aaaa"]("R,b,a,k,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += -1.00 P(i,j) r0 t2_aaaa(b,a,k,j) t1_1_aa(c,i) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["aaaa"]("b,a,k,j") * tmps_["134_aa_Loo"]("L,i,k") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +1.00 r2_abab(a,b,k,j) t1_1_aa(c,i) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["134_aa_Loo"]("L,i,k") * r2["abab"]("R,a,b,k,j");

        // D_oovv_abab += +1.00 r0 t2_abab(a,b,k,j) t1_1_aa(c,i) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,b,k,j") * tmps_["134_aa_Loo"]("L,i,k") * r0("R");

        // D_oooo_abab += +1.00 d_bb(j,l) r0 t1_1_aa(a,i) l1_1_aa(k,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += tmps_["134_aa_Loo"]("L,i,k") * tmps_["112_bb_Loo"]("R,j,l");

        // flops: o2v0L2  = o2v1L2 o3v2L2 o2v0L2 o3v2L2 o2v0L2 o2v1L2 o2v0L2 o3v2L2 o2v0L2 o3v2L2 o2v0L2
        //  mems: o2v0L2  = o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2
        tmps_["135_aa_LLoo"]("L,R,i,j")  = 2.00 * l1["aa"]("L,i,a") * r1["aa"]("R,a,j");
        tmps_["135_aa_LLoo"]("L,R,i,j") += l2["aaaa"]("L,l,i,a,c") * r2["aaaa"]("R,a,c,l,j");
        tmps_["135_aa_LLoo"]("L,R,i,j") += l2_1["aaaa"]("L,l,i,a,c") * r2_1["aaaa"]("R,a,c,l,j");
        tmps_["135_aa_LLoo"]("L,R,i,j") += 2.00 * l1_1["aa"]("L,i,a") * r1_1["aa"]("R,a,j");
        tmps_["135_aa_LLoo"]("L,R,i,j") += 2.00 * l2_1["abab"]("L,i,k,a,b") * r2_1["abab"]("R,a,b,j,k");
        tmps_["135_aa_LLoo"]("L,R,i,j") += 2.00 * l2["abab"]("L,i,k,a,b") * r2["abab"]("R,a,b,j,k");

        // D_oooo_abab += +0.50 d_bb(j,l) r2_aaaa(a,b,m,i) l2_aaaa(m,k,a,b)
        //             += +1.00 d_bb(j,l) r1_aa(a,i) l1_aa(k,a)
        //             += +0.50 d_bb(j,l) r2_1_aaaa(a,b,m,i) l2_1_aaaa(m,k,a,b)
        //             += +1.00 d_bb(j,l) r1_1_aa(a,i) l1_1_aa(k,a)
        //             += +0.50 d_bb(j,l) r2_1_abab(a,b,i,m) l2_1_abab(k,m,a,b)
        //             += +0.50 d_bb(j,l) r2_1_abab(b,a,i,m) l2_1_abab(k,m,b,a)
        //             += +0.50 d_bb(j,l) r2_abab(a,b,i,m) l2_abab(k,m,a,b)
        //             += +0.50 d_bb(j,l) r2_abab(b,a,i,m) l2_abab(k,m,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += 0.50 * tmps_["135_aa_LLoo"]("L,R,k,i") * Id["bb_oo"]("j,l");

        // D_oovv_aaaa += -0.50 P(i,j) r2_aaaa(c,d,l,i) t2_aaaa(b,a,k,j) l2_aaaa(l,k,c,d)
        //             += -1.00 P(i,j) r1_aa(c,i) t2_aaaa(b,a,k,j) l1_aa(k,c)
        //             += -0.50 P(i,j) r2_1_aaaa(c,d,l,i) t2_aaaa(b,a,k,j) l2_1_aaaa(l,k,c,d)
        //             += -1.00 P(i,j) r1_1_aa(c,i) t2_aaaa(b,a,k,j) l1_1_aa(k,c)
        //             += -0.50 P(i,j) r2_1_abab(c,d,i,l) t2_aaaa(b,a,k,j) l2_1_abab(k,l,c,d)
        //             += -0.50 P(i,j) r2_1_abab(d,c,i,l) t2_aaaa(b,a,k,j) l2_1_abab(k,l,d,c)
        //             += -0.50 P(i,j) r2_abab(c,d,i,l) t2_aaaa(b,a,k,j) l2_abab(k,l,c,d)
        //             += -0.50 P(i,j) r2_abab(d,c,i,l) t2_aaaa(b,a,k,j) l2_abab(k,l,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t2["aaaa"]("b,a,k,j") * tmps_["135_aa_LLoo"]("L,R,k,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r2_aaaa(c,d,l,i) t2_abab(a,b,k,j) l2_aaaa(l,k,c,d)
        //             += +1.00 r1_aa(c,i) t2_abab(a,b,k,j) l1_aa(k,c)
        //             += +0.50 r2_1_aaaa(c,d,l,i) t2_abab(a,b,k,j) l2_1_aaaa(l,k,c,d)
        //             += +1.00 r1_1_aa(c,i) t2_abab(a,b,k,j) l1_1_aa(k,c)
        //             += +0.50 r2_1_abab(c,d,i,l) t2_abab(a,b,k,j) l2_1_abab(k,l,c,d)
        //             += +0.50 r2_1_abab(d,c,i,l) t2_abab(a,b,k,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r2_abab(c,d,i,l) t2_abab(a,b,k,j) l2_abab(k,l,c,d)
        //             += +0.50 r2_abab(d,c,i,l) t2_abab(a,b,k,j) l2_abab(k,l,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * t2["abab"]("a,b,k,j") * tmps_["135_aa_LLoo"]("L,R,k,i");

        // flops: o2v0L1  = o2v0L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["136_aa_Loo"]("R,i,j")  = Id["aa_oo"]("i,j") * r0_1("R");

        // D_oooo_abab += +0.50 d_aa(i,k) r0_1 t2_abab(b,a,m,j) l2_1_abab(m,l,b,a)
        //             += +0.50 d_aa(i,k) r0_1 t2_abab(a,b,m,j) l2_1_abab(m,l,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += tmps_["96_bb_Loo"]("L,j,l") * tmps_["136_aa_Loo"]("R,i,k");

        // D_oooo_abab += +1.00 d_aa(i,k) r0_1 t1_bb(a,j) l1_1_bb(l,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += tmps_["108_bb_Loo"]("L,j,l") * tmps_["136_aa_Loo"]("R,i,k");

        // D_oooo_abab += +0.50 d_aa(i,k) r0_1 t2_bbbb(b,a,m,j) l2_1_bbbb(m,l,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += 0.50 * tmps_["107_bb_Loo"]("L,j,l") * tmps_["136_aa_Loo"]("R,i,k");

        // flops: o2v0L2  = o1v1L2 o1v1L2 o1v1L2 o2v1L2
        //  mems: o2v0L2  = o1v1L2 o1v1L2 o1v1L2 o2v0L2
        tmps_["137_aa_LLoo"]("L,R,i,j")  = (-1.00 * tmps_["127_aa_LLov"]("L,R,i,a") + tmps_["129_aa_LLov"]("L,R,i,a") + tmps_["130_aa_LLov"]("L,R,i,a") + -1.00 * tmps_["128_aa_LLov"]("L,R,i,a")) * t1["aa"]("a,j");

        // D_oooo_abab += +1.00 d_bb(j,l) r1_bb(b,m) t1_aa(a,i) l2_abab(k,m,a,b)
        //             += +1.00 d_bb(j,l) r1_1_bb(b,m) t1_aa(a,i) l2_1_abab(k,m,a,b)
        //             += -1.00 d_bb(j,l) r1_aa(b,m) t1_aa(a,i) l2_aaaa(m,k,a,b)
        //             += -1.00 d_bb(j,l) r1_1_aa(b,m) t1_aa(a,i) l2_1_aaaa(m,k,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += tmps_["137_aa_LLoo"]("L,R,k,i") * Id["bb_oo"]("j,l");

        // D_oovv_aaaa += +1.00 P(i,j) r1_bb(d,l) t2_aaaa(b,a,k,i) t1_aa(c,j) l2_abab(k,l,c,d)
        //             += +1.00 P(i,j) r1_1_bb(d,l) t2_aaaa(b,a,k,i) t1_aa(c,j) l2_1_abab(k,l,c,d)
        //             += -1.00 P(i,j) r1_aa(d,l) t2_aaaa(b,a,k,i) t1_aa(c,j) l2_aaaa(l,k,c,d)
        //             += -1.00 P(i,j) r1_1_aa(d,l) t2_aaaa(b,a,k,i) t1_aa(c,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["aaaa"]("b,a,k,i") * tmps_["137_aa_LLoo"]("L,R,k,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +1.00 r1_bb(d,l) t2_abab(a,b,k,j) t1_aa(c,i) l2_abab(k,l,c,d)
        //             += +1.00 r1_1_bb(d,l) t2_abab(a,b,k,j) t1_aa(c,i) l2_1_abab(k,l,c,d)
        //             += -1.00 r1_aa(d,l) t2_abab(a,b,k,j) t1_aa(c,i) l2_aaaa(l,k,c,d)
        //             += -1.00 r1_1_aa(d,l) t2_abab(a,b,k,j) t1_aa(c,i) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,b,k,j") * tmps_["137_aa_LLoo"]("L,R,k,i");

        // flops: o2v0L2  = o1v1L2 o2v1L2
        //  mems: o2v0L2  = o1v1L2 o2v0L2
        tmps_["138_aa_LLoo"]("L,R,i,j")  = (tmps_["117_aa_LLov"]("L,R,i,a") + -1.00 * tmps_["118_aa_LLov"]("L,R,i,a")) * t1_1["aa"]("a,j");

        // D_oooo_abab += -1.00 d_bb(j,l) r1_aa(b,m) t1_1_aa(a,i) l2_1_aaaa(m,k,a,b)
        //             += +1.00 d_bb(j,l) r1_bb(b,m) t1_1_aa(a,i) l2_1_abab(k,m,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") -= tmps_["138_aa_LLoo"]("L,R,k,i") * Id["bb_oo"]("j,l");

        // D_oovv_aaaa += +1.00 P(i,j) r1_aa(d,l) t2_aaaa(b,a,k,j) t1_1_aa(c,i) l2_1_aaaa(l,k,c,d)
        //             += -1.00 P(i,j) r1_bb(d,l) t2_aaaa(b,a,k,j) t1_1_aa(c,i) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["aaaa"]("b,a,k,j") * tmps_["138_aa_LLoo"]("L,R,k,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r1_aa(d,l) t2_abab(a,b,k,j) t1_1_aa(c,i) l2_1_aaaa(l,k,c,d)
        //             += +1.00 r1_bb(d,l) t2_abab(a,b,k,j) t1_1_aa(c,i) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,b,k,j") * tmps_["138_aa_LLoo"]("L,R,k,i");

        // flops: o4v0L2  = o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2
        //  mems: o4v0L2  = o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2
        tmps_["139_aaaa_LLoooo"]("L,R,i,j,k,l")  = 0.50 * tmps_["131_aa_Loo"]("L,i,j") * tmps_["136_aa_Loo"]("R,k,l");
        tmps_["139_aaaa_LLoooo"]("L,R,i,j,k,l") += 0.50 * tmps_["135_aa_LLoo"]("L,R,j,i") * Id["aa_oo"]("k,l");
        tmps_["139_aaaa_LLoooo"]("L,R,i,j,k,l") += tmps_["137_aa_LLoo"]("L,R,j,i") * Id["aa_oo"]("k,l");
        tmps_["139_aaaa_LLoooo"]("L,R,i,j,k,l") += tmps_["133_aa_Loo"]("L,i,j") * tmps_["103_aa_Loo"]("R,k,l");
        tmps_["139_aaaa_LLoooo"]("L,R,i,j,k,l") += tmps_["134_aa_Loo"]("L,i,j") * tmps_["103_aa_Loo"]("R,k,l");
        tmps_["139_aaaa_LLoooo"]("L,R,i,j,k,l") += tmps_["63_aa_Loo"]("L,i,j") * tmps_["136_aa_Loo"]("R,k,l");
        tmps_["139_aaaa_LLoooo"]("L,R,i,j,k,l") += Id["aa_oo"]("k,j") * tmps_["138_aa_LLoo"]("L,R,l,i");
        tmps_["139_aaaa_LLoooo"]("L,R,i,j,k,l") += tmps_["132_aa_Loo"]("L,i,j") * tmps_["136_aa_Loo"]("R,k,l");

        // D_oooo_aaaa += +1.00 P(i,j) d_aa(j,l) r0 t1_aa(a,i) l1_aa(k,a)
        //             += +1.00 P(i,j) d_aa(j,l) r1_bb(b,m) t1_aa(a,i) l2_abab(k,m,a,b)
        //             += +1.00 P(i,j) d_aa(j,l) r1_1_bb(b,m) t1_aa(a,i) l2_1_abab(k,m,a,b)
        //             += -1.00 P(i,j) d_aa(j,l) r1_aa(b,m) t1_aa(a,i) l2_aaaa(m,k,a,b)
        //             += -1.00 P(i,j) d_aa(j,l) r1_1_aa(b,m) t1_aa(a,i) l2_1_aaaa(m,k,a,b)
        //             += +0.50 P(i,j) d_aa(j,l) r2_aaaa(a,b,m,i) l2_aaaa(m,k,a,b)
        //             += +1.00 P(i,j) d_aa(j,l) r1_aa(a,i) l1_aa(k,a)
        //             += +0.50 P(i,j) d_aa(j,l) r2_1_aaaa(a,b,m,i) l2_1_aaaa(m,k,a,b)
        //             += +1.00 P(i,j) d_aa(j,l) r1_1_aa(a,i) l1_1_aa(k,a)
        //             += +0.50 P(i,j) d_aa(j,l) r2_1_abab(a,b,i,m) l2_1_abab(k,m,a,b)
        //             += +0.50 P(i,j) d_aa(j,l) r2_1_abab(b,a,i,m) l2_1_abab(k,m,b,a)
        //             += +0.50 P(i,j) d_aa(j,l) r2_abab(a,b,i,m) l2_abab(k,m,a,b)
        //             += +0.50 P(i,j) d_aa(j,l) r2_abab(b,a,i,m) l2_abab(k,m,b,a)
        //             += +0.50 P(i,j) d_aa(j,l) r0_1 t2_aaaa(b,a,m,i) l2_1_aaaa(m,k,b,a)
        //             += +1.00 P(i,j) d_aa(j,l) r0 t1_1_aa(a,i) l1_1_aa(k,a)
        //             += +0.50 P(i,j) d_aa(j,l) r0_1 t2_abab(b,a,i,m) l2_1_abab(k,m,b,a)
        //             += +0.50 P(i,j) d_aa(j,l) r0_1 t2_abab(a,b,i,m) l2_1_abab(k,m,a,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(b,m) t1_1_aa(a,i) l2_1_aaaa(m,l,a,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_bb(b,m) t1_1_aa(a,i) l2_1_abab(l,m,a,b)
        //             += +1.00 P(i,j) d_aa(j,l) r0_1 t1_aa(a,i) l1_1_aa(k,a)
        D_oooo_aaaa("L,R,i,j,k,l") += tmps_["139_aaaa_LLoooo"]("L,R,i,k,j,l");
        D_oooo_aaaa("L,R,i,j,k,l") -= tmps_["139_aaaa_LLoooo"]("L,R,j,k,i,l");

        // D_oooo_aaaa += -1.00 P(i,j) d_aa(j,k) r0 t1_aa(a,i) l1_aa(l,a)
        //             += -1.00 P(i,j) d_aa(j,k) r1_bb(b,m) t1_aa(a,i) l2_abab(l,m,a,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_1_bb(b,m) t1_aa(a,i) l2_1_abab(l,m,a,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(b,m) t1_aa(a,i) l2_aaaa(m,l,a,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_1_aa(b,m) t1_aa(a,i) l2_1_aaaa(m,l,a,b)
        //             += -0.50 P(i,j) d_aa(j,k) r2_aaaa(a,b,m,i) l2_aaaa(m,l,a,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(a,i) l1_aa(l,a)
        //             += -0.50 P(i,j) d_aa(j,k) r2_1_aaaa(a,b,m,i) l2_1_aaaa(m,l,a,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_1_aa(a,i) l1_1_aa(l,a)
        //             += -0.50 P(i,j) d_aa(j,k) r2_1_abab(a,b,i,m) l2_1_abab(l,m,a,b)
        //             += -0.50 P(i,j) d_aa(j,k) r2_1_abab(b,a,i,m) l2_1_abab(l,m,b,a)
        //             += -0.50 P(i,j) d_aa(j,k) r2_abab(a,b,i,m) l2_abab(l,m,a,b)
        //             += -0.50 P(i,j) d_aa(j,k) r2_abab(b,a,i,m) l2_abab(l,m,b,a)
        //             += -0.50 P(i,j) d_aa(j,k) r0_1 t2_aaaa(b,a,m,i) l2_1_aaaa(m,l,b,a)
        //             += -1.00 P(i,j) d_aa(j,k) r0 t1_1_aa(a,i) l1_1_aa(l,a)
        //             += -0.50 P(i,j) d_aa(j,k) r0_1 t2_abab(b,a,i,m) l2_1_abab(l,m,b,a)
        //             += -0.50 P(i,j) d_aa(j,k) r0_1 t2_abab(a,b,i,m) l2_1_abab(l,m,a,b)
        //             += -1.00 P(i,j) d_aa(j,l) r1_aa(b,m) t1_1_aa(a,i) l2_1_aaaa(m,k,a,b)
        //             += +1.00 P(i,j) d_aa(j,l) r1_bb(b,m) t1_1_aa(a,i) l2_1_abab(k,m,a,b)
        //             += -1.00 P(i,j) d_aa(j,k) r0_1 t1_aa(a,i) l1_1_aa(l,a)
        D_oooo_aaaa("L,R,i,j,k,l") -= tmps_["139_aaaa_LLoooo"]("L,R,i,l,j,k");
        D_oooo_aaaa("L,R,i,j,k,l") += tmps_["139_aaaa_LLoooo"]("L,R,j,l,i,k");
        tmps_["139_aaaa_LLoooo"].~TArrayD();

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["140_aa_Lvv"]("L,a,b")  = t2["aaaa"]("a,c,i,j") * l2_1["aaaa"]("L,i,j,c,b");

        // D_oovv_aaaa += -0.50 P(a,b) r2_1_aaaa(a,d,i,j) t2_aaaa(b,c,k,l) l2_1_aaaa(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["140_aa_Lvv"]("L,b,d") * r2_1["aaaa"]("R,a,d,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += -0.50 r0_1 t2_aaaa(a,c,i,j) t2_aaaa(b,d,k,l) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") -= 0.50 * t2["aaaa"]("a,c,i,j") * tmps_["140_aa_Lvv"]("L,b,c") * r0_1("R");

        // D_oovv_aaaa += -0.50 P(a,b) r0 t2_aaaa(b,d,k,l) t2_1_aaaa(a,c,i,j) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t2_1["aaaa"]("a,c,i,j") * tmps_["140_aa_Lvv"]("L,b,c") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -0.50 r2_1_abab(d,b,i,j) t2_aaaa(a,c,k,l) l2_1_aaaa(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * tmps_["140_aa_Lvv"]("L,a,d") * r2_1["abab"]("R,d,b,i,j");

        // D_oovv_abab += -0.50 r0 t2_aaaa(a,d,k,l) t2_1_abab(c,b,i,j) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * tmps_["140_aa_Lvv"]("L,a,c") * t2_1["abab"]("c,b,i,j") * r0("R");

        // flops: o2v0L2  = o3v2L2 o2v1L2 o3v2L2 o2v0L2 o2v0L2
        //  mems: o2v0L2  = o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2
        tmps_["141_aa_LLoo"]("L,R,i,j")  = 0.50 * l2_1["aaaa"]("L,l,i,a,c") * r2["aaaa"]("R,a,c,l,j");
        tmps_["141_aa_LLoo"]("L,R,i,j") += l1_1["aa"]("L,i,a") * r1["aa"]("R,a,j");
        tmps_["141_aa_LLoo"]("L,R,i,j") += l2_1["abab"]("L,i,k,a,b") * r2["abab"]("R,a,b,j,k");

        // D_oovv_aaaa += -0.50 P(i,j) r2_abab(c,d,i,l) t2_1_aaaa(b,a,k,j) l2_1_abab(k,l,c,d)
        //             += -0.50 P(i,j) r2_abab(d,c,i,l) t2_1_aaaa(b,a,k,j) l2_1_abab(k,l,d,c)
        //             += -1.00 P(i,j) r1_aa(c,i) t2_1_aaaa(b,a,k,j) l1_1_aa(k,c)
        //             += -0.50 P(i,j) r2_aaaa(c,d,l,i) t2_1_aaaa(b,a,k,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2_1["aaaa"]("b,a,k,j") * tmps_["141_aa_LLoo"]("L,R,k,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r2_abab(c,d,i,l) t2_1_abab(a,b,k,j) l2_1_abab(k,l,c,d)
        //             += +0.50 r2_abab(d,c,i,l) t2_1_abab(a,b,k,j) l2_1_abab(k,l,d,c)
        //             += +1.00 r1_aa(c,i) t2_1_abab(a,b,k,j) l1_1_aa(k,c)
        //             += +0.50 r2_aaaa(c,d,l,i) t2_1_abab(a,b,k,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2_1["abab"]("a,b,k,j") * tmps_["141_aa_LLoo"]("L,R,k,i");

        // flops: o1v1L2  = o1v2L2 o2v1L2 o1v1L2 o1v2L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2
        tmps_["142_aa_LLov"]("R,L,i,a")  = -0.50 * tmps_["140_aa_Lvv"]("L,a,b") * r1["aa"]("R,b,i");
        tmps_["142_aa_LLov"]("R,L,i,a") += t1["aa"]("a,j") * tmps_["141_aa_LLoo"]("L,R,j,i");
        tmps_["142_aa_LLov"]("R,L,i,a") += r1["aa"]("R,b,i") * tmps_["29_aa_Lvv"]("L,b,a");

        // D_oovv_aaaa += -0.50 P(i,j) P(a,b) r1_aa(d,i) t2_abab(b,c,k,l) t1_1_aa(a,j) l2_1_abab(k,l,d,c)
        //             += -0.50 P(i,j) P(a,b) r1_aa(d,i) t2_abab(b,c,l,k) t1_1_aa(a,j) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(d,i) t2_aaaa(b,c,k,l) t1_1_aa(a,j) l2_1_aaaa(k,l,c,d)
        //             += -0.50 P(i,j) P(a,b) r2_abab(c,d,i,l) t1_aa(b,k) t1_1_aa(a,j) l2_1_abab(k,l,c,d)
        //             += -0.50 P(i,j) P(a,b) r2_abab(d,c,i,l) t1_aa(b,k) t1_1_aa(a,j) l2_1_abab(k,l,d,c)
        //             += -1.00 P(i,j) P(a,b) r1_aa(c,i) t1_aa(b,k) t1_1_aa(a,j) l1_1_aa(k,c)
        //             += -0.50 P(i,j) P(a,b) r2_aaaa(c,d,l,i) t1_aa(b,k) t1_1_aa(a,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1_1["aa"]("a,j") * tmps_["142_aa_LLov"]("R,L,i,b");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r1_aa(d,i) t2_abab(a,c,k,l) t1_1_bb(b,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r1_aa(d,i) t2_abab(a,c,l,k) t1_1_bb(b,j) l2_1_abab(l,k,d,c)
        //             += -0.50 r1_aa(d,i) t2_aaaa(a,c,k,l) t1_1_bb(b,j) l2_1_aaaa(k,l,c,d)
        //             += +0.50 r2_abab(c,d,i,l) t1_aa(a,k) t1_1_bb(b,j) l2_1_abab(k,l,c,d)
        //             += +0.50 r2_abab(d,c,i,l) t1_aa(a,k) t1_1_bb(b,j) l2_1_abab(k,l,d,c)
        //             += +1.00 r1_aa(c,i) t1_aa(a,k) t1_1_bb(b,j) l1_1_aa(k,c)
        //             += +0.50 r2_aaaa(c,d,l,i) t1_aa(a,k) t1_1_bb(b,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1_1["bb"]("b,j") * tmps_["142_aa_LLov"]("R,L,i,a");
        tmps_["142_aa_LLov"].~TArrayD();
    }
} // hilbert
#endif