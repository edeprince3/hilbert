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


    void EOM_EE_QED_RDM_21::rdm2_21_5() {

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


        // flops: o1v1L2  = o2v1L2 o1v2L2 o2v1L2 o2v2L2 o2v2L2 o2v1L2 o2v2L2 o1v1L2 o2v1L2 o1v1L2 o2v1L2 o1v1L2 o1v1L2 o1v1L2 o2v2L2 o1v1L2 o1v1L2 o1v1L2 o2v2L2 o1v1L2 o1v2L2 o1v1L2 o1v2L2 o1v1L2 o2v2L2 o1v1L2 o1v2L2 o1v1L2 o2v2L2 o1v1L2 o1v2L2 o1v1L2 o2v2L2 o1v1L2 o2v2L2 o1v1L2 o1v1L2 o1v2L2 o1v1L2 o2v2L2 o1v1L2 o2v2L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2
        tmps_["211_bb_LLvo"]("L,R,a,i")  = -0.50 * t1["bb"]("a,m") * tmps_["111_bb_LLoo"]("L,R,m,i");
        tmps_["211_bb_LLvo"]("L,R,a,i") -= 0.50 * t1["bb"]("e,i") * tmps_["210_bb_LLvv"]("L,R,e,a");
        tmps_["211_bb_LLvo"]("L,R,a,i") -= t1_1["bb"]("a,m") * tmps_["101_bb_LLoo"]("L,R,m,i");
        tmps_["211_bb_LLvo"]("L,R,a,i") -= t2["abab"]("c,a,k,i") * tmps_["127_aa_LLov"]("L,R,k,c");
        tmps_["211_bb_LLvo"]("L,R,a,i") -= r1["aa"]("R,d,l") * tmps_["209_bbaa_Lvoov"]("L,a,i,l,d");
        tmps_["211_bb_LLvo"]("L,R,a,i") -= r1_1["bb"]("R,a,j") * tmps_["96_bb_Loo"]("L,i,j");
        tmps_["211_bb_LLvo"]("L,R,a,i") += t2_1["bbbb"]("a,e,m,i") * tmps_["72_bb_LLov"]("L,R,m,e");
        tmps_["211_bb_LLvo"]("L,R,a,i") += 0.50 * r1_1["bb"]("R,a,j") * tmps_["97_bb_Loo"]("L,i,j");
        tmps_["211_bb_LLvo"]("L,R,a,i") += 0.50 * r1["bb"]("R,a,j") * tmps_["202_bb_Loo"]("L,i,j");
        tmps_["211_bb_LLvo"]("L,R,a,i") += r1["bb"]("R,b,j") * tmps_["199_bbbb_Lvoov"]("L,a,i,j,b");
        tmps_["211_bb_LLvo"]("L,R,a,i") += r1["bb"]("R,b,j") * tmps_["208_bbbb_Lvoov"]("L,a,i,j,b");
        tmps_["211_bb_LLvo"]("L,R,a,i") -= r1_1["bb"]("R,b,i") * tmps_["84_bb_Lvv"]("L,a,b");
        tmps_["211_bb_LLvo"]("L,R,a,i") -= t1_1["bb"]("e,i") * tmps_["98_bb_LLvv"]("L,R,e,a");
        tmps_["211_bb_LLvo"]("L,R,a,i") += r1_1["bb"]("R,b,j") * tmps_["124_bbbb_Lvoov"]("L,a,i,j,b");
        tmps_["211_bb_LLvo"]("L,R,a,i") += 0.50 * r1["bb"]("R,b,i") * tmps_["201_bb_Lvv"]("L,a,b");
        tmps_["211_bb_LLvo"]("L,R,a,i") -= r1_1["aa"]("R,d,l") * tmps_["93_bbaa_Lvoov"]("L,a,i,l,d");
        tmps_["211_bb_LLvo"]("L,R,a,i") += 0.50 * r1["bb"]("R,b,i") * tmps_["200_bb_Lvv"]("L,a,b");
        tmps_["211_bb_LLvo"]("L,R,a,i") += r1_1["bb"]("R,b,j") * tmps_["81_bbbb_Lvoov"]("L,a,i,j,b");
        tmps_["211_bb_LLvo"]("L,R,a,i") -= r1_1["aa"]("R,d,l") * tmps_["123_bbaa_Lvoov"]("L,a,i,l,d");
        tmps_["211_bb_LLvo"]("L,R,a,i") += 0.50 * r1_1["bb"]("R,b,i") * tmps_["100_bb_Lvv"]("L,a,b");
        tmps_["211_bb_LLvo"]("L,R,a,i") -= t2_1["abab"]("c,a,k,i") * tmps_["117_aa_LLov"]("L,R,k,c");
        tmps_["211_bb_LLvo"]("L,R,a,i") += r1["bb"]("R,b,j") * tmps_["82_bbbb_Lvoov"]("L,a,i,j,b");
        tmps_["202_bb_Loo"].~TArrayD();
        tmps_["201_bb_Lvv"].~TArrayD();
        tmps_["200_bb_Lvv"].~TArrayD();
        tmps_["199_bbbb_Lvoov"].~TArrayD();
        tmps_["124_bbbb_Lvoov"].~TArrayD();
        tmps_["123_bbaa_Lvoov"].~TArrayD();

        // D_oovv_abab += -1.00 r1_bb(d,l) t1_aa(a,i) t2_1_bbbb(b,c,k,j) l2_1_bbbb(l,k,c,d)
        //             += +0.50 r1_1_bb(b,l) t1_aa(a,i) t2_abab(d,c,k,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r1_1_bb(b,l) t1_aa(a,i) t2_abab(c,d,k,j) l2_1_abab(k,l,c,d)
        //             += -0.50 r1_1_bb(b,l) t1_aa(a,i) t2_bbbb(d,c,k,j) l2_1_bbbb(l,k,d,c)
        //             += -0.50 r1_bb(b,l) t1_aa(a,i) t2_1_bbbb(d,c,k,j) l2_1_bbbb(l,k,d,c)
        //             += -0.50 r1_bb(b,l) t1_aa(a,i) t2_bbbb(d,c,k,j) l2_bbbb(l,k,d,c)
        //             += +1.00 r1_aa(d,l) t1_aa(a,i) t2_bbbb(b,c,k,j) l2_abab(l,k,d,c)
        //             += +1.00 r1_aa(d,l) t1_aa(a,i) t2_1_bbbb(b,c,k,j) l2_1_abab(l,k,d,c)
        //             += +1.00 r1_aa(d,l) t1_aa(a,i) t2_abab(c,b,k,j) l2_aaaa(l,k,c,d)
        //             += -1.00 r1_bb(d,l) t1_aa(a,i) t2_bbbb(b,c,k,j) l2_bbbb(l,k,c,d)
        //             += +0.50 r2_abab(c,d,l,j) t1_aa(a,i) t1_1_bb(b,k) l2_1_abab(l,k,c,d)
        //             += +0.50 r2_abab(d,c,l,j) t1_aa(a,i) t1_1_bb(b,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r2_bbbb(c,d,l,j) t1_aa(a,i) t1_1_bb(b,k) l2_1_bbbb(l,k,c,d)
        //             += +1.00 r1_bb(c,j) t1_aa(a,i) t1_1_bb(b,k) l1_1_bb(k,c)
        //             += +0.50 r2_bbbb(b,d,l,k) t1_aa(a,i) t1_bb(c,j) l2_bbbb(l,k,c,d)
        //             += +0.50 r2_1_abab(d,b,l,k) t1_aa(a,i) t1_bb(c,j) l2_1_abab(l,k,d,c)
        //             += +0.50 r2_1_abab(d,b,k,l) t1_aa(a,i) t1_bb(c,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r2_1_bbbb(b,d,l,k) t1_aa(a,i) t1_bb(c,j) l2_1_bbbb(l,k,c,d)
        //             += +0.50 r2_abab(d,b,l,k) t1_aa(a,i) t1_bb(c,j) l2_abab(l,k,d,c)
        //             += +0.50 r2_abab(d,b,k,l) t1_aa(a,i) t1_bb(c,j) l2_abab(k,l,d,c)
        //             += -1.00 r1_bb(d,l) t1_aa(a,i) t2_1_abab(c,b,k,j) l2_1_abab(k,l,c,d)
        //             += +0.50 r1_1_bb(d,j) t1_aa(a,i) t2_abab(c,b,k,l) l2_1_abab(k,l,c,d)
        //             += +0.50 r1_1_bb(d,j) t1_aa(a,i) t2_abab(c,b,l,k) l2_1_abab(l,k,c,d)
        //             += +0.50 r2_abab(d,b,l,k) t1_aa(a,i) t1_1_bb(c,j) l2_1_abab(l,k,d,c)
        //             += +0.50 r2_abab(d,b,k,l) t1_aa(a,i) t1_1_bb(c,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r2_bbbb(b,d,l,k) t1_aa(a,i) t1_1_bb(c,j) l2_1_bbbb(l,k,c,d)
        //             += -1.00 r1_1_bb(d,l) t1_aa(a,i) t2_bbbb(b,c,k,j) l2_1_bbbb(l,k,c,d)
        //             += -0.50 r1_bb(d,j) t1_aa(a,i) t2_1_bbbb(b,c,k,l) l2_1_bbbb(k,l,c,d)
        //             += +1.00 r1_1_aa(d,l) t1_aa(a,i) t2_bbbb(b,c,k,j) l2_1_abab(l,k,d,c)
        //             += -0.50 r1_bb(d,j) t1_aa(a,i) t2_bbbb(b,c,k,l) l2_bbbb(k,l,c,d)
        //             += -1.00 r1_1_bb(d,l) t1_aa(a,i) t2_abab(c,b,k,j) l2_1_abab(k,l,c,d)
        //             += +1.00 r1_1_aa(d,l) t1_aa(a,i) t2_abab(c,b,k,j) l2_1_aaaa(l,k,c,d)
        //             += +0.50 r2_1_bbbb(c,d,l,j) t1_aa(a,i) t1_bb(b,k) l2_1_bbbb(l,k,c,d)
        //             += +0.50 r2_abab(c,d,l,j) t1_aa(a,i) t1_bb(b,k) l2_abab(l,k,c,d)
        //             += +0.50 r2_abab(d,c,l,j) t1_aa(a,i) t1_bb(b,k) l2_abab(l,k,d,c)
        //             += +0.50 r2_bbbb(c,d,l,j) t1_aa(a,i) t1_bb(b,k) l2_bbbb(l,k,c,d)
        //             += +1.00 r1_bb(c,j) t1_aa(a,i) t1_bb(b,k) l1_bb(k,c)
        //             += +0.50 r2_1_abab(c,d,l,j) t1_aa(a,i) t1_bb(b,k) l2_1_abab(l,k,c,d)
        //             += +0.50 r2_1_abab(d,c,l,j) t1_aa(a,i) t1_bb(b,k) l2_1_abab(l,k,d,c)
        //             += +1.00 r1_1_bb(c,j) t1_aa(a,i) t1_bb(b,k) l1_1_bb(k,c)
        //             += -0.50 r1_1_bb(d,j) t1_aa(a,i) t2_bbbb(b,c,k,l) l2_1_bbbb(k,l,c,d)
        //             += +1.00 r1_aa(d,l) t1_aa(a,i) t2_1_abab(c,b,k,j) l2_1_aaaa(l,k,c,d)
        //             += -1.00 r1_bb(d,l) t1_aa(a,i) t2_abab(c,b,k,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["aa"]("a,i") * tmps_["211_bb_LLvo"]("L,R,b,j");

        // D_oovv_bbbb += -1.00 P(i,j) P(a,b) r1_bb(d,l) t1_bb(b,j) t2_1_bbbb(a,c,k,i) l2_1_bbbb(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_1_bb(a,l) t1_bb(b,j) t2_abab(d,c,k,i) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_1_bb(a,l) t1_bb(b,j) t2_abab(c,d,k,i) l2_1_abab(k,l,c,d)
        //             += -0.50 P(i,j) P(a,b) r1_1_bb(a,l) t1_bb(b,j) t2_bbbb(d,c,k,i) l2_1_bbbb(l,k,d,c)
        //             += -0.50 P(i,j) P(a,b) r1_bb(a,l) t1_bb(b,j) t2_1_bbbb(d,c,k,i) l2_1_bbbb(l,k,d,c)
        //             += -0.50 P(i,j) P(a,b) r1_bb(a,l) t1_bb(b,j) t2_bbbb(d,c,k,i) l2_bbbb(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r1_aa(d,l) t2_bbbb(a,c,k,i) t1_bb(b,j) l2_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r1_aa(d,l) t1_bb(b,j) t2_1_bbbb(a,c,k,i) l2_1_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r1_aa(d,l) t2_abab(c,a,k,i) t1_bb(b,j) l2_aaaa(l,k,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_bb(d,l) t2_bbbb(a,c,k,i) t1_bb(b,j) l2_bbbb(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_abab(c,d,l,i) t1_bb(b,j) t1_1_bb(a,k) l2_1_abab(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_abab(d,c,l,i) t1_bb(b,j) t1_1_bb(a,k) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r2_bbbb(c,d,l,i) t1_bb(b,j) t1_1_bb(a,k) l2_1_bbbb(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_bb(c,i) t1_bb(b,j) t1_1_bb(a,k) l1_1_bb(k,c)
        //             += +0.50 P(i,j) P(a,b) r2_bbbb(a,d,l,k) t1_bb(b,j) t1_bb(c,i) l2_bbbb(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_1_abab(d,a,l,k) t1_bb(b,j) t1_bb(c,i) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r2_1_abab(d,a,k,l) t1_bb(b,j) t1_bb(c,i) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r2_1_bbbb(a,d,l,k) t1_bb(b,j) t1_bb(c,i) l2_1_bbbb(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_abab(d,a,l,k) t1_bb(b,j) t1_bb(c,i) l2_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r2_abab(d,a,k,l) t1_bb(b,j) t1_bb(c,i) l2_abab(k,l,d,c)
        //             += -1.00 P(i,j) P(a,b) r1_bb(d,l) t1_bb(b,j) t2_1_abab(c,a,k,i) l2_1_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_1_bb(d,i) t2_abab(c,a,k,l) t1_bb(b,j) l2_1_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_1_bb(d,i) t2_abab(c,a,l,k) t1_bb(b,j) l2_1_abab(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_abab(d,a,l,k) t1_bb(b,j) t1_1_bb(c,i) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r2_abab(d,a,k,l) t1_bb(b,j) t1_1_bb(c,i) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r2_bbbb(a,d,l,k) t1_bb(b,j) t1_1_bb(c,i) l2_1_bbbb(l,k,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_1_bb(d,l) t2_bbbb(a,c,k,i) t1_bb(b,j) l2_1_bbbb(l,k,c,d)
        //             += -0.50 P(i,j) P(a,b) r1_bb(d,i) t1_bb(b,j) t2_1_bbbb(a,c,k,l) l2_1_bbbb(k,l,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_1_aa(d,l) t2_bbbb(a,c,k,i) t1_bb(b,j) l2_1_abab(l,k,d,c)
        //             += -0.50 P(i,j) P(a,b) r1_bb(d,i) t2_bbbb(a,c,k,l) t1_bb(b,j) l2_bbbb(k,l,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_1_bb(d,l) t2_abab(c,a,k,i) t1_bb(b,j) l2_1_abab(k,l,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_1_aa(d,l) t2_abab(c,a,k,i) t1_bb(b,j) l2_1_aaaa(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_1_bbbb(c,d,l,i) t1_bb(a,k) t1_bb(b,j) l2_1_bbbb(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_abab(c,d,l,i) t1_bb(a,k) t1_bb(b,j) l2_abab(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_abab(d,c,l,i) t1_bb(a,k) t1_bb(b,j) l2_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r2_bbbb(c,d,l,i) t1_bb(a,k) t1_bb(b,j) l2_bbbb(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_bb(c,i) t1_bb(a,k) t1_bb(b,j) l1_bb(k,c)
        //             += +0.50 P(i,j) P(a,b) r2_1_abab(c,d,l,i) t1_bb(a,k) t1_bb(b,j) l2_1_abab(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_1_abab(d,c,l,i) t1_bb(a,k) t1_bb(b,j) l2_1_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r1_1_bb(c,i) t1_bb(a,k) t1_bb(b,j) l1_1_bb(k,c)
        //             += -0.50 P(i,j) P(a,b) r1_1_bb(d,i) t2_bbbb(a,c,k,l) t1_bb(b,j) l2_1_bbbb(k,l,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_aa(d,l) t1_bb(b,j) t2_1_abab(c,a,k,i) l2_1_aaaa(l,k,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_bb(d,l) t2_abab(c,a,k,i) t1_bb(b,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["211_bb_LLvo"]("L,R,a,i") * t1["bb"]("b,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o3v2L1 o3v2L1 o3v2L1 o3v1L1 o3v1L1
        //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1
        tmps_["212_bbbb_Lvooo"]("L,a,i,j,k")  = tmps_["82_bbbb_Lvoov"]("L,a,i,j,b") * t1["bb"]("b,k");
        tmps_["212_bbbb_Lvooo"]("L,a,i,j,k") += tmps_["208_bbbb_Lvoov"]("L,a,i,j,b") * t1["bb"]("b,k");
        tmps_["212_bbbb_Lvooo"]("L,a,i,j,k") -= tmps_["81_bbbb_Lvoov"]("L,a,k,j,d") * t1_1["bb"]("d,i");

        // D_oovv_bbbb += -1.00 P(i,j) P(a,b) r1_bb(a,l) t1_bb(d,j) t2_1_abab(c,b,k,i) l2_1_abab(k,l,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_bb(a,l) t2_abab(d,b,k,j) t1_1_bb(c,i) l2_1_abab(k,l,d,c)
        //             += -1.00 P(i,j) P(a,b) r1_bb(a,l) t2_abab(c,b,k,i) t1_bb(d,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["212_bbbb_Lvooo"]("L,b,i,l,j") * r1["bb"]("R,a,l");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r0 t1_bb(b,l) t1_bb(d,j) t2_1_abab(c,a,k,i) l2_1_abab(k,l,c,d)
        //             += -1.00 P(i,j) P(a,b) r0 t2_abab(d,a,k,j) t1_bb(b,l) t1_1_bb(c,i) l2_1_abab(k,l,d,c)
        //             += +1.00 P(i,j) P(a,b) r0 t2_abab(c,a,k,i) t1_bb(b,l) t1_bb(d,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["212_bbbb_Lvooo"]("L,a,i,l,j") * t1["bb"]("b,l") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["213_bbbb_Lvoov"]("L,a,i,j,b")  = t2["bbbb"]("a,c,k,i") * l2["bbbb"]("L,k,j,b,c");

        // D_oovv_bbbb += +1.00 P(i,j) r0 t2_bbbb(a,c,k,i) t2_bbbb(b,d,l,j) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["213_bbbb_Lvoov"]("L,a,i,l,d") * t2["bbbb"]("b,d,l,j") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["214_bbbb_Lvoov"]("L,a,i,j,b")  = t2_1["bbbb"]("a,c,k,i") * l2_1["bbbb"]("L,k,j,b,c");

        // flops: o3v1L1  = o3v2L1 o3v2L1 o3v2L1 o4v2L1 o3v1L1 o3v1L1
        //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1
        tmps_["215_bbbb_Lvooo"]("L,a,i,j,k")  = tmps_["213_bbbb_Lvoov"]("L,a,i,j,b") * t1["bb"]("b,k");
        tmps_["215_bbbb_Lvooo"]("L,a,i,j,k") += tmps_["214_bbbb_Lvoov"]("L,a,i,j,b") * t1["bb"]("b,k");
        tmps_["215_bbbb_Lvooo"]("L,a,i,j,k") += t1_1["bb"]("c,i") * l2_1["bbbb"]("L,l,j,b,c") * t2["bbbb"]("a,b,l,k");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r0 t1_bb(b,l) t1_bb(d,j) t2_1_bbbb(a,c,k,i) l2_1_bbbb(k,l,d,c)
        //             += +1.00 P(i,j) P(a,b) r0 t2_bbbb(a,d,k,j) t1_bb(b,l) t1_1_bb(c,i) l2_1_bbbb(k,l,d,c)
        //             += +1.00 P(i,j) P(a,b) r0 t2_bbbb(a,c,k,i) t1_bb(b,l) t1_bb(d,j) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["215_bbbb_Lvooo"]("L,a,i,l,j") * t1["bb"]("b,l") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v1L2  = o2v1L2 o2v1L2 o1v1L2 o2v1L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2
        tmps_["216_bb_LLvo"]("L,R,a,i")  = -1.00 * t1_1["bb"]("a,j") * tmps_["125_bb_LLoo"]("L,R,j,i");
        tmps_["216_bb_LLvo"]("L,R,a,i") += t1["bb"]("a,j") * tmps_["115_bb_LLoo"]("L,R,j,i");
        tmps_["216_bb_LLvo"]("L,R,a,i") += t1["bb"]("a,j") * tmps_["114_bb_LLoo"]("L,R,j,i");

        // D_oovv_abab += -1.00 r1_bb(d,l) t1_aa(a,i) t1_bb(b,k) t1_bb(c,j) l2_bbbb(l,k,c,d)
        //             += +1.00 r1_aa(d,l) t1_aa(a,i) t1_bb(b,k) t1_bb(c,j) l2_abab(l,k,d,c)
        //             += +1.00 r1_1_aa(d,l) t1_aa(a,i) t1_bb(b,k) t1_bb(c,j) l2_1_abab(l,k,d,c)
        //             += -1.00 r1_1_bb(d,l) t1_aa(a,i) t1_bb(b,k) t1_bb(c,j) l2_1_bbbb(l,k,c,d)
        //             += +1.00 r1_aa(d,l) t1_aa(a,i) t1_bb(c,j) t1_1_bb(b,k) l2_1_abab(l,k,d,c)
        //             += -1.00 r1_bb(d,l) t1_aa(a,i) t1_bb(c,j) t1_1_bb(b,k) l2_1_bbbb(l,k,c,d)
        //             += -1.00 r1_bb(d,l) t1_aa(a,i) t1_bb(b,k) t1_1_bb(c,j) l2_1_bbbb(l,k,c,d)
        //             += +1.00 r1_aa(d,l) t1_aa(a,i) t1_bb(b,k) t1_1_bb(c,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["aa"]("a,i") * tmps_["216_bb_LLvo"]("L,R,b,j");

        // D_oovv_bbbb += -1.00 P(i,j) P(a,b) r1_bb(d,l) t1_bb(a,k) t1_bb(b,j) t1_bb(c,i) l2_bbbb(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_aa(d,l) t1_bb(a,k) t1_bb(b,j) t1_bb(c,i) l2_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r1_1_aa(d,l) t1_bb(a,k) t1_bb(b,j) t1_bb(c,i) l2_1_abab(l,k,d,c)
        //             += -1.00 P(i,j) P(a,b) r1_1_bb(d,l) t1_bb(a,k) t1_bb(b,j) t1_bb(c,i) l2_1_bbbb(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_aa(d,l) t1_bb(b,j) t1_bb(c,i) t1_1_bb(a,k) l2_1_abab(l,k,d,c)
        //             += -1.00 P(i,j) P(a,b) r1_bb(d,l) t1_bb(b,j) t1_bb(c,i) t1_1_bb(a,k) l2_1_bbbb(l,k,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_bb(d,l) t1_bb(a,k) t1_bb(b,j) t1_1_bb(c,i) l2_1_bbbb(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_aa(d,l) t1_bb(a,k) t1_bb(b,j) t1_1_bb(c,i) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["216_bb_LLvo"]("L,R,a,i") * t1["bb"]("b,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["217_bb_Lvv"]("L,a,b")  = l2_1["abab"]("L,i,j,c,a") * t2_1["abab"]("c,b,i,j");

        // D_oovv_abab += +0.50 r2_abab(a,d,i,j) t2_1_abab(c,b,k,l) l2_1_abab(k,l,c,d)
        //             += +0.50 r2_abab(a,d,i,j) t2_1_abab(c,b,l,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["217_bb_Lvv"]("L,d,b") * r2["abab"]("R,a,d,i,j");

        // D_oovv_abab += +0.50 r0 t2_abab(a,d,i,j) t2_1_abab(c,b,k,l) l2_1_abab(k,l,c,d)
        //             += +0.50 r0 t2_abab(a,d,i,j) t2_1_abab(c,b,l,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,d,i,j") * tmps_["217_bb_Lvv"]("L,d,b") * r0("R");

        // D_oovv_bbbb += +0.50 P(a,b) r2_bbbb(a,d,i,j) t2_1_abab(c,b,k,l) l2_1_abab(k,l,c,d)
        //             += +0.50 P(a,b) r2_bbbb(a,d,i,j) t2_1_abab(c,b,l,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["217_bb_Lvv"]("L,d,b") * r2["bbbb"]("R,a,d,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += -0.50 P(a,b) r0 t2_bbbb(b,d,i,j) t2_1_abab(c,a,k,l) l2_1_abab(k,l,c,d)
        //             += -0.50 P(a,b) r0 t2_bbbb(b,d,i,j) t2_1_abab(c,a,l,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["217_bb_Lvv"]("L,d,a") * t2["bbbb"]("b,d,i,j") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v1L2  = o0v2L1 o1v2L2
        //  mems: o1v1L2  = o0v2L1 o1v1L2
        tmps_["218_bb_LLvo"]("L,R,a,i")  = (tmps_["217_bb_Lvv"]("L,b,a") + tmps_["85_bb_Lvv"]("L,b,a")) * r1["bb"]("R,b,i");

        // D_oovv_abab += +0.50 r1_bb(d,j) t1_aa(a,i) t2_1_abab(c,b,k,l) l2_1_abab(k,l,c,d)
        //             += +0.50 r1_bb(d,j) t1_aa(a,i) t2_1_abab(c,b,l,k) l2_1_abab(l,k,c,d)
        //             += +0.50 r1_bb(d,j) t1_aa(a,i) t2_abab(c,b,k,l) l2_abab(k,l,c,d)
        //             += +0.50 r1_bb(d,j) t1_aa(a,i) t2_abab(c,b,l,k) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1["aa"]("a,i") * tmps_["218_bb_LLvo"]("L,R,b,j");

        // D_oovv_bbbb += +0.50 P(i,j) P(a,b) r1_bb(d,i) t1_bb(b,j) t2_1_abab(c,a,k,l) l2_1_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_bb(d,i) t1_bb(b,j) t2_1_abab(c,a,l,k) l2_1_abab(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_bb(d,i) t2_abab(c,a,k,l) t1_bb(b,j) l2_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_bb(d,i) t2_abab(c,a,l,k) t1_bb(b,j) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("b,j") * tmps_["218_bb_LLvo"]("L,R,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["219_bb_Lvv"]("L,a,b")  = t2["bbbb"]("a,c,i,j") * l2_1["bbbb"]("L,i,j,b,c");

        // D_oovv_bbbb += -0.50 r0_1 t2_bbbb(a,c,k,l) t2_bbbb(b,d,i,j) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") -= 0.50 * tmps_["219_bb_Lvv"]("L,a,d") * t2["bbbb"]("b,d,i,j") * r0_1("R");

        // flops: o1v1L1  = o2v2 o2v2 o2v2L1 o2v2L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o1v2L1 o1v1L1 o1v2L1 o1v1L1
        //  mems: o1v1L1  = o2v2 o2v2 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
        tmps_["220_bb_Lvo"]("L,a,i")  = (t2["bbbb"]("a,d,l,i") + t1["bb"]("a,l") * t1["bb"]("d,i")) * l1_1["bb"]("L,l,d");
        tmps_["220_bb_Lvo"]("L,a,i") -= l1_1["aa"]("L,j,c") * t2["abab"]("c,a,j,i");
        tmps_["220_bb_Lvo"]("L,a,i") += t1["bb"]("a,k") * tmps_["96_bb_Loo"]("L,i,k");
        tmps_["220_bb_Lvo"]("L,a,i") += 0.50 * t1["bb"]("a,k") * tmps_["107_bb_Loo"]("L,i,k");
        tmps_["220_bb_Lvo"]("L,a,i") += 0.50 * t1["bb"]("b,i") * tmps_["219_bb_Lvv"]("L,a,b");
        tmps_["220_bb_Lvo"]("L,a,i") += t1["bb"]("b,i") * tmps_["84_bb_Lvv"]("L,a,b");

        // D_oovv_abab += +1.00 r1_1_aa(a,i) t2_bbbb(b,c,k,j) l1_1_bb(k,c)
        //             += +1.00 r1_1_aa(a,i) t1_bb(b,k) t1_bb(c,j) l1_1_bb(k,c)
        //             += -1.00 r1_1_aa(a,i) t2_abab(c,b,k,j) l1_1_aa(k,c)
        //             += +0.50 r1_1_aa(a,i) t1_bb(b,l) t2_abab(d,c,k,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r1_1_aa(a,i) t1_bb(b,l) t2_abab(c,d,k,j) l2_1_abab(k,l,c,d)
        //             += +0.50 r1_1_aa(a,i) t1_bb(b,l) t2_bbbb(d,c,k,j) l2_1_bbbb(k,l,d,c)
        //             += +0.50 r1_1_aa(a,i) t2_bbbb(b,c,k,l) t1_bb(d,j) l2_1_bbbb(k,l,d,c)
        //             += +0.50 r1_1_aa(a,i) t2_abab(c,b,k,l) t1_bb(d,j) l2_1_abab(k,l,c,d)
        //             += +0.50 r1_1_aa(a,i) t2_abab(c,b,l,k) t1_bb(d,j) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["220_bb_Lvo"]("L,b,j") * r1_1["aa"]("R,a,i");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r1_1_bb(a,i) t2_bbbb(b,c,k,j) l1_1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_1_bb(a,i) t1_bb(b,k) t1_bb(c,j) l1_1_bb(k,c)
        //             += -1.00 P(i,j) P(a,b) r1_1_bb(a,i) t2_abab(c,b,k,j) l1_1_aa(k,c)
        //             += +0.50 P(i,j) P(a,b) r1_1_bb(a,i) t1_bb(b,l) t2_abab(d,c,k,j) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_1_bb(a,i) t1_bb(b,l) t2_abab(c,d,k,j) l2_1_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_1_bb(a,i) t1_bb(b,l) t2_bbbb(d,c,k,j) l2_1_bbbb(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_1_bb(a,i) t2_bbbb(b,c,k,l) t1_bb(d,j) l2_1_bbbb(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_1_bb(a,i) t2_abab(c,b,k,l) t1_bb(d,j) l2_1_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_1_bb(a,i) t2_abab(c,b,l,k) t1_bb(d,j) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["220_bb_Lvo"]("L,b,j") * r1_1["bb"]("R,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["221_bb_Lvv"]("L,a,b")  = l2_1["bbbb"]("L,i,j,a,c") * t2_1["bbbb"]("b,c,i,j");

        // D_oovv_abab += +0.50 r0 t2_abab(a,d,i,j) t2_1_bbbb(b,c,k,l) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * t2["abab"]("a,d,i,j") * tmps_["221_bb_Lvv"]("L,d,b") * r0("R");

        // D_oovv_bbbb += -0.50 P(a,b) r0 t2_bbbb(b,d,i,j) t2_1_bbbb(a,c,k,l) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["221_bb_Lvv"]("L,d,a") * t2["bbbb"]("b,d,i,j") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["222_bb_Lvv"]("L,a,b")  = t2["bbbb"]("a,c,i,j") * l2["bbbb"]("L,i,j,b,c");

        // D_oovv_bbbb += -0.50 r0 t2_bbbb(a,c,k,l) t2_bbbb(b,d,i,j) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") -= 0.50 * tmps_["222_bb_Lvv"]("L,a,d") * t2["bbbb"]("b,d,i,j") * r0("R");

        // flops: o1v1L1  = o2v1L1 o1v2L1 o2v2L1 o2v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o1v2L1 o1v2L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1 o1v2L1 o1v1L1 o1v1L1 o1v2L1 o1v1L1 o1v1L1 o1v2L1 o1v1L1
        //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o2v0L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o2v0L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o2v0L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
        tmps_["223_bb_Lov"]("L,i,a")  = -0.50 * t1_1["bb"]("a,j") * tmps_["97_bb_Loo"]("L,i,j");
        tmps_["223_bb_Lov"]("L,i,a") -= 0.50 * tmps_["100_bb_Lvv"]("L,a,c") * t1_1["bb"]("c,i");
        tmps_["223_bb_Lov"]("L,i,a") -= l1["aa"]("L,l,d") * t2["abab"]("d,a,l,i");
        tmps_["223_bb_Lov"]("L,i,a") += l1_1["bb"]("L,j,c") * t1["bb"]("c,i") * t1_1["bb"]("a,j");
        tmps_["223_bb_Lov"]("L,i,a") -= l1_1["aa"]("L,l,d") * t2_1["abab"]("d,a,l,i");
        tmps_["223_bb_Lov"]("L,i,a") += l1["bb"]("L,j,c") * t1["bb"]("c,i") * t1["bb"]("a,j");
        tmps_["223_bb_Lov"]("L,i,a") += l1["bb"]("L,j,c") * t2["bbbb"]("a,c,j,i");
        tmps_["223_bb_Lov"]("L,i,a") += l1_1["bb"]("L,j,c") * t1_1["bb"]("c,i") * t1["bb"]("a,j");
        tmps_["223_bb_Lov"]("L,i,a") += l1_1["bb"]("L,j,c") * t2_1["bbbb"]("a,c,j,i");
        tmps_["223_bb_Lov"]("L,i,a") += t1["bb"]("b,i") * tmps_["217_bb_Lvv"]("L,b,a");
        tmps_["223_bb_Lov"]("L,i,a") += t1_1["bb"]("c,i") * tmps_["84_bb_Lvv"]("L,a,c");
        tmps_["223_bb_Lov"]("L,i,a") += t1_1["bb"]("a,j") * tmps_["96_bb_Loo"]("L,i,j");
        tmps_["223_bb_Lov"]("L,i,a") += t1["bb"]("b,i") * tmps_["85_bb_Lvv"]("L,b,a");
        tmps_["223_bb_Lov"]("L,i,a") += 0.50 * t1["bb"]("b,i") * tmps_["222_bb_Lvv"]("L,a,b");
        tmps_["223_bb_Lov"]("L,i,a") += 0.50 * t1["bb"]("b,i") * tmps_["221_bb_Lvv"]("L,b,a");
        tmps_["100_bb_Lvv"].~TArrayD();
        tmps_["97_bb_Loo"].~TArrayD();

        // D_oovv_abab += +0.50 r1_aa(a,i) t1_bb(d,j) t2_1_abab(c,b,k,l) l2_1_abab(k,l,c,d)
        //             += +0.50 r1_aa(a,i) t1_bb(d,j) t2_1_abab(c,b,l,k) l2_1_abab(l,k,c,d)
        //             += +0.50 r1_aa(a,i) t2_abab(d,b,k,l) t1_1_bb(c,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r1_aa(a,i) t2_abab(d,b,l,k) t1_1_bb(c,j) l2_1_abab(l,k,d,c)
        //             += +1.00 r1_aa(a,i) t1_bb(c,j) t1_1_bb(b,k) l1_1_bb(k,c)
        //             += -1.00 r1_aa(a,i) t2_abab(c,b,k,j) l1_aa(k,c)
        //             += -1.00 r1_aa(a,i) t2_1_abab(c,b,k,j) l1_1_aa(k,c)
        //             += +1.00 r1_aa(a,i) t1_bb(b,k) t1_bb(c,j) l1_bb(k,c)
        //             += +1.00 r1_aa(a,i) t2_bbbb(b,c,k,j) l1_bb(k,c)
        //             += +1.00 r1_aa(a,i) t1_bb(b,k) t1_1_bb(c,j) l1_1_bb(k,c)
        //             += +1.00 r1_aa(a,i) t2_1_bbbb(b,c,k,j) l1_1_bb(k,c)
        //             += +0.50 r1_aa(a,i) t2_abab(d,c,l,j) t1_1_bb(b,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r1_aa(a,i) t2_abab(c,d,l,j) t1_1_bb(b,k) l2_1_abab(l,k,c,d)
        //             += +0.50 r1_aa(a,i) t2_abab(c,b,k,l) t1_bb(d,j) l2_abab(k,l,c,d)
        //             += +0.50 r1_aa(a,i) t2_abab(c,b,l,k) t1_bb(d,j) l2_abab(l,k,c,d)
        //             += -0.50 r1_aa(a,i) t2_bbbb(b,d,k,l) t1_1_bb(c,j) l2_1_bbbb(k,l,d,c)
        //             += +0.50 r1_aa(a,i) t2_bbbb(b,c,k,l) t1_bb(d,j) l2_bbbb(k,l,d,c)
        //             += -0.50 r1_aa(a,i) t2_bbbb(d,c,l,j) t1_1_bb(b,k) l2_1_bbbb(k,l,d,c)
        //             += +0.50 r1_aa(a,i) t1_bb(d,j) t2_1_bbbb(b,c,k,l) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["223_bb_Lov"]("L,j,b") * r1["aa"]("R,a,i");

        // D_oovv_bbbb += +0.50 P(i,j) P(a,b) r1_bb(a,i) t1_bb(d,j) t2_1_abab(c,b,k,l) l2_1_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,i) t1_bb(d,j) t2_1_abab(c,b,l,k) l2_1_abab(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,i) t2_abab(d,b,k,l) t1_1_bb(c,j) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,i) t2_abab(d,b,l,k) t1_1_bb(c,j) l2_1_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r1_bb(a,i) t1_bb(c,j) t1_1_bb(b,k) l1_1_bb(k,c)
        //             += -1.00 P(i,j) P(a,b) r1_bb(a,i) t2_abab(c,b,k,j) l1_aa(k,c)
        //             += -1.00 P(i,j) P(a,b) r1_bb(a,i) t2_1_abab(c,b,k,j) l1_1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_bb(a,i) t1_bb(b,k) t1_bb(c,j) l1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_bb(a,i) t2_bbbb(b,c,k,j) l1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_bb(a,i) t1_bb(b,k) t1_1_bb(c,j) l1_1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_bb(a,i) t2_1_bbbb(b,c,k,j) l1_1_bb(k,c)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,i) t2_abab(d,c,l,j) t1_1_bb(b,k) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,i) t2_abab(c,d,l,j) t1_1_bb(b,k) l2_1_abab(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,i) t2_abab(c,b,k,l) t1_bb(d,j) l2_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,i) t2_abab(c,b,l,k) t1_bb(d,j) l2_abab(l,k,c,d)
        //             += -0.50 P(i,j) P(a,b) r1_bb(a,i) t2_bbbb(b,d,k,l) t1_1_bb(c,j) l2_1_bbbb(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,i) t2_bbbb(b,c,k,l) t1_bb(d,j) l2_bbbb(k,l,d,c)
        //             += -0.50 P(i,j) P(a,b) r1_bb(a,i) t2_bbbb(d,c,l,j) t1_1_bb(b,k) l2_1_bbbb(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,i) t1_bb(d,j) t2_1_bbbb(b,c,k,l) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["223_bb_Lov"]("L,j,b") * r1["bb"]("R,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v1L1  = o1v1L1
        //  mems: o1v1L1  = o1v1L1
        tmps_["224_bb_Lvo"]("R,a,i")  = r0("R") * t1["bb"]("a,i");

        // D_oovv_bbbb += +0.50 P(i,j) P(a,b) r0 t1_bb(b,j) t1_bb(d,i) t2_1_abab(c,a,k,l) l2_1_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r0 t1_bb(b,j) t1_bb(d,i) t2_1_abab(c,a,l,k) l2_1_abab(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r0 t2_abab(d,a,k,l) t1_bb(b,j) t1_1_bb(c,i) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t2_abab(d,a,l,k) t1_bb(b,j) t1_1_bb(c,i) l2_1_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r0 t1_bb(b,j) t1_bb(c,i) t1_1_bb(a,k) l1_1_bb(k,c)
        //             += -1.00 P(i,j) P(a,b) r0 t2_abab(c,a,k,i) t1_bb(b,j) l1_aa(k,c)
        //             += -1.00 P(i,j) P(a,b) r0 t1_bb(b,j) t2_1_abab(c,a,k,i) l1_1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r0 t1_bb(a,k) t1_bb(b,j) t1_bb(c,i) l1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r0 t2_bbbb(a,c,k,i) t1_bb(b,j) l1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r0 t1_bb(a,k) t1_bb(b,j) t1_1_bb(c,i) l1_1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r0 t1_bb(b,j) t2_1_bbbb(a,c,k,i) l1_1_bb(k,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_bb(b,j) t2_abab(d,c,l,i) t1_1_bb(a,k) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_bb(b,j) t2_abab(c,d,l,i) t1_1_bb(a,k) l2_1_abab(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r0 t2_abab(c,a,k,l) t1_bb(b,j) t1_bb(d,i) l2_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r0 t2_abab(c,a,l,k) t1_bb(b,j) t1_bb(d,i) l2_abab(l,k,c,d)
        //             += -0.50 P(i,j) P(a,b) r0 t2_bbbb(a,d,k,l) t1_bb(b,j) t1_1_bb(c,i) l2_1_bbbb(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t2_bbbb(a,c,k,l) t1_bb(b,j) t1_bb(d,i) l2_bbbb(k,l,d,c)
        //             += -0.50 P(i,j) P(a,b) r0 t1_bb(b,j) t2_bbbb(d,c,l,i) t1_1_bb(a,k) l2_1_bbbb(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_bb(b,j) t1_bb(d,i) t2_1_bbbb(a,c,k,l) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["223_bb_Lov"]("L,i,a") * tmps_["224_bb_Lvo"]("R,b,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += -1.00 P(i,j) r0 t1_bb(a,i) t1_bb(b,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["192_bb_Lvo"]("L,b,j") * tmps_["224_bb_Lvo"]("R,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
        tmps_["192_bb_Lvo"].~TArrayD();

        // D_oovo_bbbb += -0.50 P(i,j) r0 t1_bb(a,j) t2_abab(c,b,l,i) l2_abab(l,k,c,b)
        //             += -0.50 P(i,j) r0 t1_bb(a,j) t2_abab(b,c,l,i) l2_abab(l,k,b,c)
        //             += +0.50 r1_aa(a,i) t1_bb(b,l) t2_1_abab(d,c,k,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r1_aa(a,i) t1_bb(b,l) t2_1_abab(c,d,k,j) l2_1_abab(k,l,c,d)
        // flops: o3v1L2 += o4v1L2
        //  mems: o3v1L2 += o3v1L2
        tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k")  = tmps_["60_bbbb_Loooo"]("L,i,k,j,l") * tmps_["224_bb_Lvo"]("R,a,j");
        D_oovo_bbbb("L,R,i,j,a,k") -= tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k");
        D_oovo_bbbb("L,R,i,j,a,k") += tmps_["perm_bbbb_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_bbbb_LLvooo"].~TArrayD();

        // D_oovv_bbbb += +0.50 P(i,j) P(a,b) r0 t1_bb(a,l) t1_bb(b,j) t2_bbbb(d,c,k,i) l2_bbbb(k,l,d,c)
        //             += +0.50 r1_aa(a,j) t2_1_bbbb(c,b,l,i) l2_1_bbbb(l,k,c,b)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o3v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["195_bbb_Loov"]("L,i,k,a") * tmps_["224_bb_Lvo"]("R,b,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
        tmps_["195_bbb_Loov"].~TArrayD();

        // D_ooov_bbbb += +0.50 P(i,j) r0 t1_bb(a,j) t2_abab(c,b,l,i) l2_abab(l,k,c,b)
        //             += +0.50 P(i,j) r0 t1_bb(a,j) t2_abab(b,c,l,i) l2_abab(l,k,b,c)
        //             += +0.50 r1_bb(b,l) t1_aa(a,i) t2_1_abab(d,c,k,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r1_bb(b,l) t1_aa(a,i) t2_1_abab(c,d,k,j) l2_1_abab(k,l,c,d)
        // flops: o3v1L2 += o3v2L1 o3v2L1 o4v0L1 o4v1L2
        //  mems: o3v1L2 += o2v0L1 o2v0L1 o4v0L1 o3v1L2
        tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k")  = (t2["abab"]("c,b,l,i") * l2["abab"]("L,l,k,c,b") + t2_1["abab"]("d,c,k,j") * l2_1["abab"]("L,k,l,d,c")) * tmps_["224_bb_Lvo"]("R,a,j");
        D_ooov_bbbb("L,R,i,j,k,a") += tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k");
        D_ooov_bbbb("L,R,i,j,k,a") -= tmps_["perm_bbbb_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_bbbb_LLvooo"].~TArrayD();

        // D_oovv_abab += +0.50 r0 t1_aa(a,l) t1_bb(b,j) t2_1_aaaa(d,c,k,i) l2_1_aaaa(k,l,d,c)
        //             += +0.50 P(i,j) r0 t1_aa(a,j) t2_aaaa(c,b,l,i) l2_aaaa(l,k,c,b)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o3v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * tmps_["193_aaa_Loov"]("L,i,k,a") * tmps_["224_bb_Lvo"]("R,b,j");
        tmps_["193_aaa_Loov"].~TArrayD();

        // D_oovo_bbbb += -0.50 P(i,j) r0 t1_bb(a,j) t2_bbbb(c,b,l,i) l2_bbbb(l,k,c,b)
        //             += -0.50 r0 t1_aa(a,j) t2_1_bbbb(c,b,l,i) l2_1_bbbb(l,k,c,b)
        // flops: o3v1L2 += o3v2L1 o3v2L1 o2v0L1 o3v1L2
        //  mems: o3v1L2 += o2v0L1 o2v0L1 o2v0L1 o3v1L2
        tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k")  = 0.50 * (t2["bbbb"]("c,b,l,i") * l2["bbbb"]("L,l,k,c,b") + t2_1["bbbb"]("c,b,l,i") * l2_1["bbbb"]("L,l,k,c,b")) * tmps_["224_bb_Lvo"]("R,a,j");
        D_oovo_bbbb("L,R,i,j,a,k") -= tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k");
        D_oovo_bbbb("L,R,i,j,a,k") += tmps_["perm_bbbb_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_bbbb_LLvooo"].~TArrayD();

        // D_oovv_abab += +0.50 r0 t1_aa(a,l) t1_bb(b,j) t2_1_abab(d,c,i,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r0 t1_aa(a,l) t1_bb(b,j) t2_1_abab(c,d,i,k) l2_1_abab(l,k,c,d)
        //             += -0.50 P(i,j) r1_aa(a,i) t2_abab(c,b,j,l) l2_abab(k,l,c,b)
        //             += -0.50 P(i,j) r1_aa(a,i) t2_abab(b,c,j,l) l2_abab(k,l,b,c)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o4v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["194_aaaa_Looov"]("L,i,j,k,a") * tmps_["224_bb_Lvo"]("R,b,j");
        tmps_["194_aaaa_Looov"].~TArrayD();

        // D_oovv_bbbb += +0.50 P(i,j) P(a,b) r0 t1_bb(a,l) t1_bb(b,j) t2_abab(d,c,k,i) l2_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_bb(a,l) t1_bb(b,j) t2_abab(c,d,k,i) l2_abab(k,l,c,d)
        //             += +0.50 r1_aa(a,j) t2_1_abab(c,b,l,i) l2_1_abab(l,k,c,b)
        //             += +0.50 r1_aa(a,j) t2_1_abab(b,c,l,i) l2_1_abab(l,k,b,c)
        // flops: o2v2L2 += o3v2L1 o3v2L1 o3v0L1 o3v1L1 o3v2L2
        //  mems: o2v2L2 += o2v0L1 o2v0L1 o3v0L1 o2v1L1 o3v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = (t2["abab"]("d,c,k,i") * l2["abab"]("L,k,l,d,c") + t2_1["abab"]("c,b,l,i") * l2_1["abab"]("L,l,k,c,b")) * t1["bb"]("a,l") * tmps_["224_bb_Lvo"]("R,b,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v1L1  = o1v1L1
        //  mems: o1v1L1  = o1v1L1
        tmps_["225_bb_Lvo"]("R,a,i")  = r0_1("R") * t1["bb"]("a,i");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r0_1 t2_bbbb(a,c,k,i) t1_bb(b,j) l1_1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r0_1 t1_bb(a,k) t1_bb(b,j) t1_bb(c,i) l1_1_bb(k,c)
        //             += -1.00 P(i,j) P(a,b) r0_1 t2_abab(c,a,k,i) t1_bb(b,j) l1_1_aa(k,c)
        //             += +0.50 P(i,j) P(a,b) r0_1 t1_bb(a,l) t1_bb(b,j) t2_abab(d,c,k,i) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0_1 t1_bb(a,l) t1_bb(b,j) t2_abab(c,d,k,i) l2_1_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r0_1 t1_bb(a,l) t1_bb(b,j) t2_bbbb(d,c,k,i) l2_1_bbbb(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0_1 t2_bbbb(a,c,k,l) t1_bb(b,j) t1_bb(d,i) l2_1_bbbb(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0_1 t2_abab(c,a,k,l) t1_bb(b,j) t1_bb(d,i) l2_1_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r0_1 t2_abab(c,a,l,k) t1_bb(b,j) t1_bb(d,i) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["220_bb_Lvo"]("L,a,i") * tmps_["225_bb_Lvo"]("R,b,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v1L1  = o1v1L1
        //  mems: o1v1L1  = o1v1L1
        tmps_["226_bb_Lvo"]("R,a,i")  = r0("R") * t1_1["bb"]("a,i");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r0 t2_bbbb(b,c,k,j) t1_1_bb(a,i) l1_1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r0 t1_bb(b,k) t1_bb(c,j) t1_1_bb(a,i) l1_1_bb(k,c)
        //             += -1.00 P(i,j) P(a,b) r0 t2_abab(c,b,k,j) t1_1_bb(a,i) l1_1_aa(k,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_bb(b,l) t2_abab(d,c,k,j) t1_1_bb(a,i) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_bb(b,l) t2_abab(c,d,k,j) t1_1_bb(a,i) l2_1_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r0 t1_bb(b,l) t2_bbbb(d,c,k,j) t1_1_bb(a,i) l2_1_bbbb(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t2_bbbb(b,c,k,l) t1_bb(d,j) t1_1_bb(a,i) l2_1_bbbb(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t2_abab(c,b,k,l) t1_bb(d,j) t1_1_bb(a,i) l2_1_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r0 t2_abab(c,b,l,k) t1_bb(d,j) t1_1_bb(a,i) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["220_bb_Lvo"]("L,b,j") * tmps_["226_bb_Lvo"]("R,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o0v0L2  = o2v2L2 o2v2L2 o0v0L2 o2v2L2 o0v0L2 o1v1L2 o0v0L2 o1v1L2 o0v0L2 o1v1L2 o0v0L2 o2v2L2 o0v0L2 o1v1L2 o0v0L2 o2v2L2 o0v0L2 o2v2L2 o0v0L2
        //  mems: o0v0L2  = o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2
        tmps_["229_LL"]("L,R")  = 0.25 * l2["aaaa"]("L,i,k,a,d") * r2["aaaa"]("R,a,d,i,k");
        tmps_["229_LL"]("L,R") += l2_1["abab"]("L,i,j,a,b") * r2_1["abab"]("R,a,b,i,j");
        tmps_["229_LL"]("L,R") += 0.25 * l2["bbbb"]("L,l,j,c,b") * r2["bbbb"]("R,c,b,l,j");
        tmps_["229_LL"]("L,R") += l1_1["aa"]("L,k,a") * r1_1["aa"]("R,a,k");
        tmps_["229_LL"]("L,R") += l1_1["bb"]("L,j,c") * r1_1["bb"]("R,c,j");
        tmps_["229_LL"]("L,R") += l1["aa"]("L,k,a") * r1["aa"]("R,a,k");
        tmps_["229_LL"]("L,R") += 0.25 * l2_1["bbbb"]("L,l,j,c,b") * r2_1["bbbb"]("R,c,b,l,j");
        tmps_["229_LL"]("L,R") += l1["bb"]("L,j,c") * r1["bb"]("R,c,j");
        tmps_["229_LL"]("L,R") += 0.25 * l2_1["aaaa"]("L,i,k,a,d") * r2_1["aaaa"]("R,a,d,i,k");
        tmps_["229_LL"]("L,R") += l2["abab"]("L,i,j,a,b") * r2["abab"]("R,a,b,i,j");

        // D_oovv_aaaa += +0.25 r2_1_abab(c,d,l,k) t2_aaaa(b,a,i,j) l2_1_abab(l,k,c,d)
        //             += +0.25 r2_1_abab(c,d,k,l) t2_aaaa(b,a,i,j) l2_1_abab(k,l,c,d)
        //             += +0.25 r2_1_abab(d,c,l,k) t2_aaaa(b,a,i,j) l2_1_abab(l,k,d,c)
        //             += +0.25 r2_1_abab(d,c,k,l) t2_aaaa(b,a,i,j) l2_1_abab(k,l,d,c)
        //             += +0.25 r2_aaaa(c,d,l,k) t2_aaaa(b,a,i,j) l2_aaaa(l,k,c,d)
        //             += +0.25 r2_bbbb(c,d,l,k) t2_aaaa(b,a,i,j) l2_bbbb(l,k,c,d)
        //             += +1.00 r1_1_aa(c,k) t2_aaaa(b,a,i,j) l1_1_aa(k,c)
        //             += +1.00 r1_1_bb(c,k) t2_aaaa(b,a,i,j) l1_1_bb(k,c)
        //             += +1.00 r1_aa(c,k) t2_aaaa(b,a,i,j) l1_aa(k,c)
        //             += +0.25 r2_1_bbbb(c,d,l,k) t2_aaaa(b,a,i,j) l2_1_bbbb(l,k,c,d)
        //             += +1.00 r1_bb(c,k) t2_aaaa(b,a,i,j) l1_bb(k,c)
        //             += +0.25 r2_1_aaaa(c,d,l,k) t2_aaaa(b,a,i,j) l2_1_aaaa(l,k,c,d)
        //             += +0.25 r2_abab(c,d,l,k) t2_aaaa(b,a,i,j) l2_abab(l,k,c,d)
        //             += +0.25 r2_abab(c,d,k,l) t2_aaaa(b,a,i,j) l2_abab(k,l,c,d)
        //             += +0.25 r2_abab(d,c,l,k) t2_aaaa(b,a,i,j) l2_abab(l,k,d,c)
        //             += +0.25 r2_abab(d,c,k,l) t2_aaaa(b,a,i,j) l2_abab(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") += t2["aaaa"]("b,a,i,j") * tmps_["229_LL"]("L,R");

        // D_oovv_abab += -0.25 r2_1_abab(c,d,l,k) t2_abab(a,b,i,j) l2_1_abab(l,k,c,d)
        //             += -0.25 r2_1_abab(c,d,k,l) t2_abab(a,b,i,j) l2_1_abab(k,l,c,d)
        //             += -0.25 r2_1_abab(d,c,l,k) t2_abab(a,b,i,j) l2_1_abab(l,k,d,c)
        //             += -0.25 r2_1_abab(d,c,k,l) t2_abab(a,b,i,j) l2_1_abab(k,l,d,c)
        //             += -0.25 r2_aaaa(c,d,l,k) t2_abab(a,b,i,j) l2_aaaa(l,k,c,d)
        //             += -0.25 r2_bbbb(c,d,l,k) t2_abab(a,b,i,j) l2_bbbb(l,k,c,d)
        //             += -1.00 r1_1_aa(c,k) t2_abab(a,b,i,j) l1_1_aa(k,c)
        //             += -1.00 r1_1_bb(c,k) t2_abab(a,b,i,j) l1_1_bb(k,c)
        //             += -1.00 r1_aa(c,k) t2_abab(a,b,i,j) l1_aa(k,c)
        //             += -0.25 r2_1_bbbb(c,d,l,k) t2_abab(a,b,i,j) l2_1_bbbb(l,k,c,d)
        //             += -1.00 r1_bb(c,k) t2_abab(a,b,i,j) l1_bb(k,c)
        //             += -0.25 r2_1_aaaa(c,d,l,k) t2_abab(a,b,i,j) l2_1_aaaa(l,k,c,d)
        //             += -0.25 r2_abab(c,d,l,k) t2_abab(a,b,i,j) l2_abab(l,k,c,d)
        //             += -0.25 r2_abab(c,d,k,l) t2_abab(a,b,i,j) l2_abab(k,l,c,d)
        //             += -0.25 r2_abab(d,c,l,k) t2_abab(a,b,i,j) l2_abab(l,k,d,c)
        //             += -0.25 r2_abab(d,c,k,l) t2_abab(a,b,i,j) l2_abab(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,b,i,j") * tmps_["229_LL"]("L,R");

        // D_oovv_abab += -0.25 r2_1_abab(c,d,l,k) t1_aa(a,i) t1_bb(b,j) l2_1_abab(l,k,c,d)
        //             += -0.25 r2_1_abab(c,d,k,l) t1_aa(a,i) t1_bb(b,j) l2_1_abab(k,l,c,d)
        //             += -0.25 r2_1_abab(d,c,l,k) t1_aa(a,i) t1_bb(b,j) l2_1_abab(l,k,d,c)
        //             += -0.25 r2_1_abab(d,c,k,l) t1_aa(a,i) t1_bb(b,j) l2_1_abab(k,l,d,c)
        //             += -0.25 r2_aaaa(c,d,l,k) t1_aa(a,i) t1_bb(b,j) l2_aaaa(l,k,c,d)
        //             += -0.25 r2_bbbb(c,d,l,k) t1_aa(a,i) t1_bb(b,j) l2_bbbb(l,k,c,d)
        //             += -1.00 r1_1_aa(c,k) t1_aa(a,i) t1_bb(b,j) l1_1_aa(k,c)
        //             += -1.00 r1_1_bb(c,k) t1_aa(a,i) t1_bb(b,j) l1_1_bb(k,c)
        //             += -1.00 r1_aa(c,k) t1_aa(a,i) t1_bb(b,j) l1_aa(k,c)
        //             += -0.25 r2_1_bbbb(c,d,l,k) t1_aa(a,i) t1_bb(b,j) l2_1_bbbb(l,k,c,d)
        //             += -1.00 r1_bb(c,k) t1_aa(a,i) t1_bb(b,j) l1_bb(k,c)
        //             += -0.25 r2_1_aaaa(c,d,l,k) t1_aa(a,i) t1_bb(b,j) l2_1_aaaa(l,k,c,d)
        //             += -0.25 r2_abab(c,d,l,k) t1_aa(a,i) t1_bb(b,j) l2_abab(l,k,c,d)
        //             += -0.25 r2_abab(c,d,k,l) t1_aa(a,i) t1_bb(b,j) l2_abab(k,l,c,d)
        //             += -0.25 r2_abab(d,c,l,k) t1_aa(a,i) t1_bb(b,j) l2_abab(l,k,d,c)
        //             += -0.25 r2_abab(d,c,k,l) t1_aa(a,i) t1_bb(b,j) l2_abab(k,l,d,c)
        // flops: o2v2L2 += o2v2 o2v2L2
        //  mems: o2v2L2 += o2v2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["aa"]("a,i") * t1["bb"]("b,j") * tmps_["229_LL"]("L,R");

        // D_oovv_bbbb += +0.25 r2_1_abab(c,d,l,k) t2_bbbb(b,a,i,j) l2_1_abab(l,k,c,d)
        //             += +0.25 r2_1_abab(c,d,k,l) t2_bbbb(b,a,i,j) l2_1_abab(k,l,c,d)
        //             += +0.25 r2_1_abab(d,c,l,k) t2_bbbb(b,a,i,j) l2_1_abab(l,k,d,c)
        //             += +0.25 r2_1_abab(d,c,k,l) t2_bbbb(b,a,i,j) l2_1_abab(k,l,d,c)
        //             += +0.25 r2_aaaa(c,d,l,k) t2_bbbb(b,a,i,j) l2_aaaa(l,k,c,d)
        //             += +0.25 r2_bbbb(c,d,l,k) t2_bbbb(b,a,i,j) l2_bbbb(l,k,c,d)
        //             += +1.00 r1_1_aa(c,k) t2_bbbb(b,a,i,j) l1_1_aa(k,c)
        //             += +1.00 r1_1_bb(c,k) t2_bbbb(b,a,i,j) l1_1_bb(k,c)
        //             += +1.00 r1_aa(c,k) t2_bbbb(b,a,i,j) l1_aa(k,c)
        //             += +0.25 r2_1_bbbb(c,d,l,k) t2_bbbb(b,a,i,j) l2_1_bbbb(l,k,c,d)
        //             += +1.00 r1_bb(c,k) t2_bbbb(b,a,i,j) l1_bb(k,c)
        //             += +0.25 r2_1_aaaa(c,d,l,k) t2_bbbb(b,a,i,j) l2_1_aaaa(l,k,c,d)
        //             += +0.25 r2_abab(c,d,l,k) t2_bbbb(b,a,i,j) l2_abab(l,k,c,d)
        //             += +0.25 r2_abab(c,d,k,l) t2_bbbb(b,a,i,j) l2_abab(k,l,c,d)
        //             += +0.25 r2_abab(d,c,l,k) t2_bbbb(b,a,i,j) l2_abab(l,k,d,c)
        //             += +0.25 r2_abab(d,c,k,l) t2_bbbb(b,a,i,j) l2_abab(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += t2["bbbb"]("b,a,i,j") * tmps_["229_LL"]("L,R");

        // D_oovv_bbbb += -0.25 P(i,j) r2_1_abab(c,d,l,k) t1_bb(a,i) t1_bb(b,j) l2_1_abab(l,k,c,d)
        //             += -0.25 P(i,j) r2_1_abab(c,d,k,l) t1_bb(a,i) t1_bb(b,j) l2_1_abab(k,l,c,d)
        //             += -0.25 P(i,j) r2_1_abab(d,c,l,k) t1_bb(a,i) t1_bb(b,j) l2_1_abab(l,k,d,c)
        //             += -0.25 P(i,j) r2_1_abab(d,c,k,l) t1_bb(a,i) t1_bb(b,j) l2_1_abab(k,l,d,c)
        //             += -0.25 P(i,j) r2_aaaa(c,d,l,k) t1_bb(a,i) t1_bb(b,j) l2_aaaa(l,k,c,d)
        //             += -0.25 P(i,j) r2_bbbb(c,d,l,k) t1_bb(a,i) t1_bb(b,j) l2_bbbb(l,k,c,d)
        //             += -1.00 P(i,j) r1_1_aa(c,k) t1_bb(a,i) t1_bb(b,j) l1_1_aa(k,c)
        //             += -1.00 P(i,j) r1_1_bb(c,k) t1_bb(a,i) t1_bb(b,j) l1_1_bb(k,c)
        //             += -1.00 P(i,j) r1_aa(c,k) t1_bb(a,i) t1_bb(b,j) l1_aa(k,c)
        //             += -0.25 P(i,j) r2_1_bbbb(c,d,l,k) t1_bb(a,i) t1_bb(b,j) l2_1_bbbb(l,k,c,d)
        //             += -1.00 P(i,j) r1_bb(c,k) t1_bb(a,i) t1_bb(b,j) l1_bb(k,c)
        //             += -0.25 P(i,j) r2_1_aaaa(c,d,l,k) t1_bb(a,i) t1_bb(b,j) l2_1_aaaa(l,k,c,d)
        //             += -0.25 P(i,j) r2_abab(c,d,l,k) t1_bb(a,i) t1_bb(b,j) l2_abab(l,k,c,d)
        //             += -0.25 P(i,j) r2_abab(c,d,k,l) t1_bb(a,i) t1_bb(b,j) l2_abab(k,l,c,d)
        //             += -0.25 P(i,j) r2_abab(d,c,l,k) t1_bb(a,i) t1_bb(b,j) l2_abab(l,k,d,c)
        //             += -0.25 P(i,j) r2_abab(d,c,k,l) t1_bb(a,i) t1_bb(b,j) l2_abab(k,l,d,c)
        // flops: o2v2L2 += o2v2 o2v2L2
        //  mems: o2v2L2 += o2v2 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("a,i") * t1["bb"]("b,j") * tmps_["229_LL"]("L,R");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oooo_abab += -0.25 d_aa(i,k) d_bb(j,l) r2_1_abab(a,b,n,m) l2_1_abab(n,m,a,b)
        //             += -0.25 d_aa(i,k) d_bb(j,l) r2_1_abab(a,b,m,n) l2_1_abab(m,n,a,b)
        //             += -0.25 d_aa(i,k) d_bb(j,l) r2_1_abab(b,a,n,m) l2_1_abab(n,m,b,a)
        //             += -0.25 d_aa(i,k) d_bb(j,l) r2_1_abab(b,a,m,n) l2_1_abab(m,n,b,a)
        //             += -0.25 d_aa(i,k) d_bb(j,l) r2_aaaa(a,b,n,m) l2_aaaa(n,m,a,b)
        //             += -0.25 d_aa(i,k) d_bb(j,l) r2_bbbb(a,b,n,m) l2_bbbb(n,m,a,b)
        //             += -1.00 d_aa(i,k) d_bb(j,l) r1_1_aa(a,m) l1_1_aa(m,a)
        //             += -1.00 d_aa(i,k) d_bb(j,l) r1_1_bb(a,m) l1_1_bb(m,a)
        //             += -1.00 d_aa(i,k) d_bb(j,l) r1_aa(a,m) l1_aa(m,a)
        //             += -0.25 d_aa(i,k) d_bb(j,l) r2_1_bbbb(a,b,n,m) l2_1_bbbb(n,m,a,b)
        //             += -1.00 d_aa(i,k) d_bb(j,l) r1_bb(a,m) l1_bb(m,a)
        //             += -0.25 d_aa(i,k) d_bb(j,l) r2_1_aaaa(a,b,n,m) l2_1_aaaa(n,m,a,b)
        //             += -0.25 d_aa(i,k) d_bb(j,l) r2_abab(a,b,n,m) l2_abab(n,m,a,b)
        //             += -0.25 d_aa(i,k) d_bb(j,l) r2_abab(a,b,m,n) l2_abab(m,n,a,b)
        //             += -0.25 d_aa(i,k) d_bb(j,l) r2_abab(b,a,n,m) l2_abab(n,m,b,a)
        //             += -0.25 d_aa(i,k) d_bb(j,l) r2_abab(b,a,m,n) l2_abab(m,n,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") -= tmps_["229_LL"]("L,R") * tmps_["197_aabb_oooo"]("i,k,j,l");
        tmps_["197_aabb_oooo"].~TArrayD();

        // D_oooo_bbbb += +0.25 P(i,j) d_bb(i,l) d_bb(j,k) r2_1_abab(a,b,n,m) l2_1_abab(n,m,a,b)
        //             += +0.25 P(i,j) d_bb(i,l) d_bb(j,k) r2_1_abab(a,b,m,n) l2_1_abab(m,n,a,b)
        //             += +0.25 P(i,j) d_bb(i,l) d_bb(j,k) r2_1_abab(b,a,n,m) l2_1_abab(n,m,b,a)
        //             += +0.25 P(i,j) d_bb(i,l) d_bb(j,k) r2_1_abab(b,a,m,n) l2_1_abab(m,n,b,a)
        //             += +0.25 P(i,j) d_bb(i,l) d_bb(j,k) r2_aaaa(a,b,n,m) l2_aaaa(n,m,a,b)
        //             += +0.25 P(i,j) d_bb(i,l) d_bb(j,k) r2_bbbb(a,b,n,m) l2_bbbb(n,m,a,b)
        //             += +1.00 P(i,j) d_bb(i,l) d_bb(j,k) r1_1_aa(a,m) l1_1_aa(m,a)
        //             += +1.00 P(i,j) d_bb(i,l) d_bb(j,k) r1_1_bb(a,m) l1_1_bb(m,a)
        //             += +1.00 P(i,j) d_bb(i,l) d_bb(j,k) r1_aa(a,m) l1_aa(m,a)
        //             += +0.25 P(i,j) d_bb(i,l) d_bb(j,k) r2_1_bbbb(a,b,n,m) l2_1_bbbb(n,m,a,b)
        //             += +1.00 P(i,j) d_bb(i,l) d_bb(j,k) r1_bb(a,m) l1_bb(m,a)
        //             += +0.25 P(i,j) d_bb(i,l) d_bb(j,k) r2_1_aaaa(a,b,n,m) l2_1_aaaa(n,m,a,b)
        //             += +0.25 P(i,j) d_bb(i,l) d_bb(j,k) r2_abab(a,b,n,m) l2_abab(n,m,a,b)
        //             += +0.25 P(i,j) d_bb(i,l) d_bb(j,k) r2_abab(a,b,m,n) l2_abab(m,n,a,b)
        //             += +0.25 P(i,j) d_bb(i,l) d_bb(j,k) r2_abab(b,a,n,m) l2_abab(n,m,b,a)
        //             += +0.25 P(i,j) d_bb(i,l) d_bb(j,k) r2_abab(b,a,m,n) l2_abab(m,n,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        tmps_["perm_bbbb_LLoooo"]("L,R,i,j,k,l")  = tmps_["229_LL"]("L,R") * tmps_["198_bbbb_oooo"]("i,l,j,k");
        D_oooo_bbbb("L,R,i,j,k,l") += tmps_["perm_bbbb_LLoooo"]("L,R,i,j,k,l");
        D_oooo_bbbb("L,R,i,j,k,l") -= tmps_["perm_bbbb_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_bbbb_LLoooo"].~TArrayD();
        tmps_["198_bbbb_oooo"].~TArrayD();

        // D_oovv_aaaa += -0.25 P(i,j) r2_1_abab(c,d,l,k) t1_aa(a,i) t1_aa(b,j) l2_1_abab(l,k,c,d)
        //             += -0.25 P(i,j) r2_1_abab(c,d,k,l) t1_aa(a,i) t1_aa(b,j) l2_1_abab(k,l,c,d)
        //             += -0.25 P(i,j) r2_1_abab(d,c,l,k) t1_aa(a,i) t1_aa(b,j) l2_1_abab(l,k,d,c)
        //             += -0.25 P(i,j) r2_1_abab(d,c,k,l) t1_aa(a,i) t1_aa(b,j) l2_1_abab(k,l,d,c)
        //             += -0.25 P(i,j) r2_aaaa(c,d,l,k) t1_aa(a,i) t1_aa(b,j) l2_aaaa(l,k,c,d)
        //             += -0.25 P(i,j) r2_bbbb(c,d,l,k) t1_aa(a,i) t1_aa(b,j) l2_bbbb(l,k,c,d)
        //             += -1.00 P(i,j) r1_1_aa(c,k) t1_aa(a,i) t1_aa(b,j) l1_1_aa(k,c)
        //             += -1.00 P(i,j) r1_1_bb(c,k) t1_aa(a,i) t1_aa(b,j) l1_1_bb(k,c)
        //             += -1.00 P(i,j) r1_aa(c,k) t1_aa(a,i) t1_aa(b,j) l1_aa(k,c)
        //             += -0.25 P(i,j) r2_1_bbbb(c,d,l,k) t1_aa(a,i) t1_aa(b,j) l2_1_bbbb(l,k,c,d)
        //             += -1.00 P(i,j) r1_bb(c,k) t1_aa(a,i) t1_aa(b,j) l1_bb(k,c)
        //             += -0.25 P(i,j) r2_1_aaaa(c,d,l,k) t1_aa(a,i) t1_aa(b,j) l2_1_aaaa(l,k,c,d)
        //             += -0.25 P(i,j) r2_abab(c,d,l,k) t1_aa(a,i) t1_aa(b,j) l2_abab(l,k,c,d)
        //             += -0.25 P(i,j) r2_abab(c,d,k,l) t1_aa(a,i) t1_aa(b,j) l2_abab(k,l,c,d)
        //             += -0.25 P(i,j) r2_abab(d,c,l,k) t1_aa(a,i) t1_aa(b,j) l2_abab(l,k,d,c)
        //             += -0.25 P(i,j) r2_abab(d,c,k,l) t1_aa(a,i) t1_aa(b,j) l2_abab(k,l,d,c)
        // flops: o2v2L2 += o2v2 o2v2L2
        //  mems: o2v2L2 += o2v2 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("a,i") * t1["aa"]("b,j") * tmps_["229_LL"]("L,R");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oooo_aaaa += +0.25 P(i,j) d_aa(i,l) d_aa(j,k) r2_1_abab(a,b,n,m) l2_1_abab(n,m,a,b)
        //             += +0.25 P(i,j) d_aa(i,l) d_aa(j,k) r2_1_abab(a,b,m,n) l2_1_abab(m,n,a,b)
        //             += +0.25 P(i,j) d_aa(i,l) d_aa(j,k) r2_1_abab(b,a,n,m) l2_1_abab(n,m,b,a)
        //             += +0.25 P(i,j) d_aa(i,l) d_aa(j,k) r2_1_abab(b,a,m,n) l2_1_abab(m,n,b,a)
        //             += +0.25 P(i,j) d_aa(i,l) d_aa(j,k) r2_aaaa(a,b,n,m) l2_aaaa(n,m,a,b)
        //             += +0.25 P(i,j) d_aa(i,l) d_aa(j,k) r2_bbbb(a,b,n,m) l2_bbbb(n,m,a,b)
        //             += +1.00 P(i,j) d_aa(i,l) d_aa(j,k) r1_1_aa(a,m) l1_1_aa(m,a)
        //             += +1.00 P(i,j) d_aa(i,l) d_aa(j,k) r1_1_bb(a,m) l1_1_bb(m,a)
        //             += +1.00 P(i,j) d_aa(i,l) d_aa(j,k) r1_aa(a,m) l1_aa(m,a)
        //             += +0.25 P(i,j) d_aa(i,l) d_aa(j,k) r2_1_bbbb(a,b,n,m) l2_1_bbbb(n,m,a,b)
        //             += +1.00 P(i,j) d_aa(i,l) d_aa(j,k) r1_bb(a,m) l1_bb(m,a)
        //             += +0.25 P(i,j) d_aa(i,l) d_aa(j,k) r2_1_aaaa(a,b,n,m) l2_1_aaaa(n,m,a,b)
        //             += +0.25 P(i,j) d_aa(i,l) d_aa(j,k) r2_abab(a,b,n,m) l2_abab(n,m,a,b)
        //             += +0.25 P(i,j) d_aa(i,l) d_aa(j,k) r2_abab(a,b,m,n) l2_abab(m,n,a,b)
        //             += +0.25 P(i,j) d_aa(i,l) d_aa(j,k) r2_abab(b,a,n,m) l2_abab(n,m,b,a)
        //             += +0.25 P(i,j) d_aa(i,l) d_aa(j,k) r2_abab(b,a,m,n) l2_abab(m,n,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l")  = tmps_["229_LL"]("L,R") * tmps_["196_aaaa_oooo"]("i,l,j,k");
        D_oooo_aaaa("L,R,i,j,k,l") += tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l");
        D_oooo_aaaa("L,R,i,j,k,l") -= tmps_["perm_aaaa_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_aaaa_LLoooo"].~TArrayD();
        tmps_["196_aaaa_oooo"].~TArrayD();

        // flops: o0v0L2  = o2v2L2 o1v1L2 o1v1L2 o0v0L2 o0v0L2 o2v2L2 o0v0L2 o2v2L2 o0v0L2
        //  mems: o0v0L2  = o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2
        tmps_["230_LL"]("L,R")  = 0.25 * l2_1["aaaa"]("L,i,l,a,d") * r2["aaaa"]("R,a,d,i,l");
        tmps_["230_LL"]("L,R") += l1_1["aa"]("L,l,a") * r1["aa"]("R,a,l");
        tmps_["230_LL"]("L,R") += l1_1["bb"]("L,j,c") * r1["bb"]("R,c,j");
        tmps_["230_LL"]("L,R") += 0.25 * l2_1["bbbb"]("L,k,j,c,b") * r2["bbbb"]("R,c,b,k,j");
        tmps_["230_LL"]("L,R") += l2_1["abab"]("L,i,j,a,b") * r2["abab"]("R,a,b,i,j");

        // D_oovv_aaaa += +1.00 r1_bb(c,k) t2_1_aaaa(b,a,i,j) l1_1_bb(k,c)
        //             += +1.00 r1_aa(c,k) t2_1_aaaa(b,a,i,j) l1_1_aa(k,c)
        //             += +0.25 r2_aaaa(c,d,l,k) t2_1_aaaa(b,a,i,j) l2_1_aaaa(l,k,c,d)
        //             += +0.25 r2_bbbb(c,d,l,k) t2_1_aaaa(b,a,i,j) l2_1_bbbb(l,k,c,d)
        //             += +0.25 r2_abab(c,d,l,k) t2_1_aaaa(b,a,i,j) l2_1_abab(l,k,c,d)
        //             += +0.25 r2_abab(c,d,k,l) t2_1_aaaa(b,a,i,j) l2_1_abab(k,l,c,d)
        //             += +0.25 r2_abab(d,c,l,k) t2_1_aaaa(b,a,i,j) l2_1_abab(l,k,d,c)
        //             += +0.25 r2_abab(d,c,k,l) t2_1_aaaa(b,a,i,j) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") += t2_1["aaaa"]("b,a,i,j") * tmps_["230_LL"]("L,R");

        // D_oovv_aaaa += -1.00 P(i,j) P(a,b) r1_bb(c,k) t1_aa(b,j) t1_1_aa(a,i) l1_1_bb(k,c)
        //             += -1.00 P(i,j) P(a,b) r1_aa(c,k) t1_aa(b,j) t1_1_aa(a,i) l1_1_aa(k,c)
        //             += -0.25 P(i,j) P(a,b) r2_aaaa(c,d,l,k) t1_aa(b,j) t1_1_aa(a,i) l2_1_aaaa(l,k,c,d)
        //             += -0.25 P(i,j) P(a,b) r2_bbbb(c,d,l,k) t1_aa(b,j) t1_1_aa(a,i) l2_1_bbbb(l,k,c,d)
        //             += -0.25 P(i,j) P(a,b) r2_abab(c,d,l,k) t1_aa(b,j) t1_1_aa(a,i) l2_1_abab(l,k,c,d)
        //             += -0.25 P(i,j) P(a,b) r2_abab(c,d,k,l) t1_aa(b,j) t1_1_aa(a,i) l2_1_abab(k,l,c,d)
        //             += -0.25 P(i,j) P(a,b) r2_abab(d,c,l,k) t1_aa(b,j) t1_1_aa(a,i) l2_1_abab(l,k,d,c)
        //             += -0.25 P(i,j) P(a,b) r2_abab(d,c,k,l) t1_aa(b,j) t1_1_aa(a,i) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o2v2 o2v2L2
        //  mems: o2v2L2 += o2v2 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("b,j") * t1_1["aa"]("a,i") * tmps_["230_LL"]("L,R");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r1_bb(c,k) t2_1_abab(a,b,i,j) l1_1_bb(k,c)
        //             += -1.00 r1_aa(c,k) t2_1_abab(a,b,i,j) l1_1_aa(k,c)
        //             += -0.25 r2_aaaa(c,d,l,k) t2_1_abab(a,b,i,j) l2_1_aaaa(l,k,c,d)
        //             += -0.25 r2_bbbb(c,d,l,k) t2_1_abab(a,b,i,j) l2_1_bbbb(l,k,c,d)
        //             += -0.25 r2_abab(c,d,l,k) t2_1_abab(a,b,i,j) l2_1_abab(l,k,c,d)
        //             += -0.25 r2_abab(c,d,k,l) t2_1_abab(a,b,i,j) l2_1_abab(k,l,c,d)
        //             += -0.25 r2_abab(d,c,l,k) t2_1_abab(a,b,i,j) l2_1_abab(l,k,d,c)
        //             += -0.25 r2_abab(d,c,k,l) t2_1_abab(a,b,i,j) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2_1["abab"]("a,b,i,j") * tmps_["230_LL"]("L,R");

        // D_oovv_abab += -1.00 r1_bb(c,k) t1_aa(a,i) t1_1_bb(b,j) l1_1_bb(k,c)
        //             += -1.00 r1_aa(c,k) t1_aa(a,i) t1_1_bb(b,j) l1_1_aa(k,c)
        //             += -0.25 r2_aaaa(c,d,l,k) t1_aa(a,i) t1_1_bb(b,j) l2_1_aaaa(l,k,c,d)
        //             += -0.25 r2_bbbb(c,d,l,k) t1_aa(a,i) t1_1_bb(b,j) l2_1_bbbb(l,k,c,d)
        //             += -0.25 r2_abab(c,d,l,k) t1_aa(a,i) t1_1_bb(b,j) l2_1_abab(l,k,c,d)
        //             += -0.25 r2_abab(c,d,k,l) t1_aa(a,i) t1_1_bb(b,j) l2_1_abab(k,l,c,d)
        //             += -0.25 r2_abab(d,c,l,k) t1_aa(a,i) t1_1_bb(b,j) l2_1_abab(l,k,d,c)
        //             += -0.25 r2_abab(d,c,k,l) t1_aa(a,i) t1_1_bb(b,j) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o2v2 o2v2L2
        //  mems: o2v2L2 += o2v2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["aa"]("a,i") * t1_1["bb"]("b,j") * tmps_["230_LL"]("L,R");

        // D_oovv_abab += -1.00 r1_bb(c,k) t1_bb(b,j) t1_1_aa(a,i) l1_1_bb(k,c)
        //             += -1.00 r1_aa(c,k) t1_bb(b,j) t1_1_aa(a,i) l1_1_aa(k,c)
        //             += -0.25 r2_aaaa(c,d,l,k) t1_bb(b,j) t1_1_aa(a,i) l2_1_aaaa(l,k,c,d)
        //             += -0.25 r2_bbbb(c,d,l,k) t1_bb(b,j) t1_1_aa(a,i) l2_1_bbbb(l,k,c,d)
        //             += -0.25 r2_abab(c,d,l,k) t1_bb(b,j) t1_1_aa(a,i) l2_1_abab(l,k,c,d)
        //             += -0.25 r2_abab(c,d,k,l) t1_bb(b,j) t1_1_aa(a,i) l2_1_abab(k,l,c,d)
        //             += -0.25 r2_abab(d,c,l,k) t1_bb(b,j) t1_1_aa(a,i) l2_1_abab(l,k,d,c)
        //             += -0.25 r2_abab(d,c,k,l) t1_bb(b,j) t1_1_aa(a,i) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o2v2 o2v2L2
        //  mems: o2v2L2 += o2v2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["bb"]("b,j") * t1_1["aa"]("a,i") * tmps_["230_LL"]("L,R");

        // D_oovv_bbbb += +1.00 r1_bb(c,k) t2_1_bbbb(b,a,i,j) l1_1_bb(k,c)
        //             += +1.00 r1_aa(c,k) t2_1_bbbb(b,a,i,j) l1_1_aa(k,c)
        //             += +0.25 r2_aaaa(c,d,l,k) t2_1_bbbb(b,a,i,j) l2_1_aaaa(l,k,c,d)
        //             += +0.25 r2_bbbb(c,d,l,k) t2_1_bbbb(b,a,i,j) l2_1_bbbb(l,k,c,d)
        //             += +0.25 r2_abab(c,d,l,k) t2_1_bbbb(b,a,i,j) l2_1_abab(l,k,c,d)
        //             += +0.25 r2_abab(c,d,k,l) t2_1_bbbb(b,a,i,j) l2_1_abab(k,l,c,d)
        //             += +0.25 r2_abab(d,c,l,k) t2_1_bbbb(b,a,i,j) l2_1_abab(l,k,d,c)
        //             += +0.25 r2_abab(d,c,k,l) t2_1_bbbb(b,a,i,j) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += t2_1["bbbb"]("b,a,i,j") * tmps_["230_LL"]("L,R");

        // D_oovv_bbbb += -1.00 P(i,j) P(a,b) r1_bb(c,k) t1_bb(b,j) t1_1_bb(a,i) l1_1_bb(k,c)
        //             += -1.00 P(i,j) P(a,b) r1_aa(c,k) t1_bb(b,j) t1_1_bb(a,i) l1_1_aa(k,c)
        //             += -0.25 P(i,j) P(a,b) r2_aaaa(c,d,l,k) t1_bb(b,j) t1_1_bb(a,i) l2_1_aaaa(l,k,c,d)
        //             += -0.25 P(i,j) P(a,b) r2_bbbb(c,d,l,k) t1_bb(b,j) t1_1_bb(a,i) l2_1_bbbb(l,k,c,d)
        //             += -0.25 P(i,j) P(a,b) r2_abab(c,d,l,k) t1_bb(b,j) t1_1_bb(a,i) l2_1_abab(l,k,c,d)
        //             += -0.25 P(i,j) P(a,b) r2_abab(c,d,k,l) t1_bb(b,j) t1_1_bb(a,i) l2_1_abab(k,l,c,d)
        //             += -0.25 P(i,j) P(a,b) r2_abab(d,c,l,k) t1_bb(b,j) t1_1_bb(a,i) l2_1_abab(l,k,d,c)
        //             += -0.25 P(i,j) P(a,b) r2_abab(d,c,k,l) t1_bb(b,j) t1_1_bb(a,i) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o2v2 o2v2L2
        //  mems: o2v2L2 += o2v2 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("b,j") * t1_1["bb"]("a,i") * tmps_["230_LL"]("L,R");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v0L1  = o2v0L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["228_bb_Loo"]("L,i,j")  = Id["bb_oo"]("i,j") * l0("L");

        // flops: o2v0L1  = o2v0L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["227_bb_Loo"]("L,i,j")  = Id["bb_oo"]("i,j") * l0_1("L");

        // flops: o3v1L2  = o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v3L1 o3v3L1 o3v2L1 o3v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L2 o3v2L2 o3v1L2 o3v3L1 o3v2L2 o3v1L2 o3v1L2 o3v1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        //  mems: o3v1L2  = o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o2v2L1 o2v2L1 o3v1L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v1L2 o3v1L2 o3v1L2 o2v2L1 o3v1L2 o3v1L2 o3v1L2 o3v1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k")  = -1.00 * Id["bb_oo"]("j,k") * tmps_["211_bb_LLvo"]("L,R,a,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += t1["bb"]("a,j") * tmps_["114_bb_LLoo"]("L,R,k,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += t1_1["bb"]("a,i") * tmps_["125_bb_LLoo"]("L,R,k,j");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += tmps_["220_bb_Lvo"]("L,a,i") * tmps_["113_bb_Loo"]("R,j,k");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += tmps_["223_bb_Lov"]("L,i,a") * tmps_["112_bb_Loo"]("R,j,k");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") -= tmps_["227_bb_Loo"]("L,j,k") * tmps_["226_bb_Lvo"]("R,a,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") -= 0.50 * t1["bb"]("a,j") * tmps_["111_bb_LLoo"]("L,R,k,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") -= t1_1["bb"]("a,j") * tmps_["101_bb_LLoo"]("L,R,k,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") -= tmps_["96_bb_Loo"]("L,i,k") * tmps_["225_bb_Lvo"]("R,a,j");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += tmps_["203_bb_LLov"]("L,R,i,a") * Id["bb_oo"]("j,k");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") -= Id["bb_oo"]("j,k") * t1["bb"]("a,i") * tmps_["229_LL"]("L,R");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") -= tmps_["110_bb_Loo"]("L,i,k") * tmps_["224_bb_Lvo"]("R,a,j");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") -= tmps_["108_bb_Loo"]("L,i,k") * tmps_["225_bb_Lvo"]("R,a,j");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") -= tmps_["227_bb_Loo"]("L,j,k") * tmps_["225_bb_Lvo"]("R,a,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") -= (-1.00 * t2["bbbb"]("a,c,n,j") * l2["bbbb"]("L,n,k,c,b") + -1.00 * l2_1["bbbb"]("L,n,k,c,b") * t2_1["bbbb"]("a,c,n,j") + tmps_["208_bbbb_Lvoov"]("L,a,j,k,b") + -1.00 * l2_1["bbbb"]("L,n,k,c,b") * t1["bb"]("c,j") * t1_1["bb"]("a,n") + tmps_["82_bbbb_Lvoov"]("L,a,j,k,b")) * r1["bb"]("R,b,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") -= tmps_["81_bbbb_Lvoov"]("L,a,j,k,b") * r1_1["bb"]("R,b,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += t2["bbbb"]("a,c,n,j") * l2_1["bbbb"]("L,n,k,c,b") * r1_1["bb"]("R,b,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") -= Id["bb_oo"]("j,k") * t1_1["bb"]("a,i") * tmps_["230_LL"]("L,R");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += tmps_["108_bb_Loo"]("L,j,k") * tmps_["226_bb_Lvo"]("R,a,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") -= tmps_["109_bb_Loo"]("L,i,k") * tmps_["224_bb_Lvo"]("R,a,j");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += 0.50 * tmps_["107_bb_Loo"]("L,j,k") * tmps_["226_bb_Lvo"]("R,a,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += tmps_["96_bb_Loo"]("L,j,k") * tmps_["226_bb_Lvo"]("R,a,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") -= tmps_["228_bb_Loo"]("L,j,k") * tmps_["224_bb_Lvo"]("R,a,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") -= 0.50 * tmps_["107_bb_Loo"]("L,i,k") * tmps_["225_bb_Lvo"]("R,a,j");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += r0("R") * tmps_["212_bbbb_Lvooo"]("L,a,i,k,j");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += t1["bb"]("a,j") * tmps_["115_bb_LLoo"]("L,R,k,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += r0("R") * tmps_["215_bbbb_Lvooo"]("L,a,i,k,j");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") -= Id["bb_oo"]("j,k") * tmps_["216_bb_LLvo"]("L,R,a,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") -= tmps_["228_bb_Loo"]("L,j,k") * r1["bb"]("R,a,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += tmps_["204_bbbb_Lvooo"]("L,a,i,k,j") * r0_1("R");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += 0.50 * tmps_["107_bb_Loo"]("L,j,k") * r1_1["bb"]("R,a,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += tmps_["206_bbbb_Lvooo"]("L,a,i,k,j") * r0_1("R");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += tmps_["96_bb_Loo"]("L,j,k") * r1_1["bb"]("R,a,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += tmps_["207_bbbb_Looov"]("L,i,k,j,a") * r0("R");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += tmps_["110_bb_Loo"]("L,j,k") * r1["bb"]("R,a,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += tmps_["108_bb_Loo"]("L,j,k") * r1_1["bb"]("R,a,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += tmps_["109_bb_Loo"]("L,j,k") * r1["bb"]("R,a,i");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") -= r1_1["bb"]("R,a,i") * tmps_["227_bb_Loo"]("L,j,k");
        tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k") += Id["bb_oo"]("j,k") * tmps_["218_bb_LLvo"]("L,R,a,i");
        tmps_["218_bb_LLvo"].~TArrayD();
        tmps_["216_bb_LLvo"].~TArrayD();
        tmps_["215_bbbb_Lvooo"].~TArrayD();
        tmps_["212_bbbb_Lvooo"].~TArrayD();
        tmps_["211_bb_LLvo"].~TArrayD();
        tmps_["207_bbbb_Looov"].~TArrayD();
        tmps_["206_bbbb_Lvooo"].~TArrayD();
        tmps_["204_bbbb_Lvooo"].~TArrayD();
        tmps_["203_bb_LLov"].~TArrayD();

        // D_ooov_bbbb += -0.50 P(i,j) r1_1_bb(a,i) t2_abab(c,b,l,j) l2_1_abab(l,k,c,b)
        //             += -0.50 P(i,j) r1_1_bb(a,i) t2_abab(b,c,l,j) l2_1_abab(l,k,b,c)
        //             += -1.00 P(i,j) r0_1 t2_bbbb(a,b,l,i) t1_bb(c,j) l2_1_bbbb(l,k,c,b)
        //             += -0.50 P(i,j) r1_1_bb(a,i) t2_bbbb(c,b,l,j) l2_1_bbbb(l,k,c,b)
        //             += -1.00 P(i,j) r0_1 t2_abab(b,a,l,i) t1_bb(c,j) l2_1_abab(l,k,b,c)
        //             += -1.00 P(i,j) r0 t1_bb(a,l) t1_bb(c,j) t1_1_bb(b,i) l2_1_bbbb(l,k,c,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(a,i) l0
        //             += -1.00 P(i,j) r1_bb(a,i) t1_1_bb(b,j) l1_1_bb(k,b)
        //             += -1.00 P(i,j) r1_1_bb(a,i) t1_bb(b,j) l1_1_bb(k,b)
        //             += -1.00 P(i,j) r1_bb(c,l) t1_bb(a,j) t1_bb(b,i) l2_bbbb(l,k,b,c)
        //             += +1.00 P(i,j) r1_aa(c,l) t1_bb(a,j) t1_bb(b,i) l2_abab(l,k,c,b)
        //             += +1.00 P(i,j) r1_1_aa(c,l) t1_bb(a,j) t1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += -1.00 P(i,j) r1_1_bb(c,l) t1_bb(a,j) t1_bb(b,i) l2_1_bbbb(l,k,b,c)
        //             += -1.00 P(i,j) r0 t1_bb(c,j) t2_1_abab(b,a,l,i) l2_1_abab(l,k,b,c)
        //             += +1.00 P(i,j) r0 t2_abab(c,a,l,j) t1_1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += -1.00 P(i,j) r0 t2_abab(b,a,l,i) t1_bb(c,j) l2_abab(l,k,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r0 t1_1_bb(a,i) l0_1
        //             += +0.50 P(i,j) r2_1_bbbb(b,c,l,i) t1_bb(a,j) l2_1_bbbb(l,k,b,c)
        //             += +0.50 P(i,j) r2_abab(b,c,l,i) t1_bb(a,j) l2_abab(l,k,b,c)
        //             += +0.50 P(i,j) r2_abab(c,b,l,i) t1_bb(a,j) l2_abab(l,k,c,b)
        //             += +0.50 P(i,j) r2_bbbb(b,c,l,i) t1_bb(a,j) l2_bbbb(l,k,b,c)
        //             += +1.00 P(i,j) r1_bb(b,i) t1_bb(a,j) l1_bb(k,b)
        //             += +0.50 P(i,j) r2_1_abab(b,c,l,i) t1_bb(a,j) l2_1_abab(l,k,b,c)
        //             += +0.50 P(i,j) r2_1_abab(c,b,l,i) t1_bb(a,j) l2_1_abab(l,k,c,b)
        //             += +1.00 P(i,j) r1_1_bb(b,i) t1_bb(a,j) l1_1_bb(k,b)
        //             += +0.50 P(i,j) r2_abab(b,c,l,i) t1_1_bb(a,j) l2_1_abab(l,k,b,c)
        //             += +0.50 P(i,j) r2_abab(c,b,l,i) t1_1_bb(a,j) l2_1_abab(l,k,c,b)
        //             += +0.50 P(i,j) r2_bbbb(b,c,l,i) t1_1_bb(a,j) l2_1_bbbb(l,k,b,c)
        //             += +1.00 P(i,j) r1_bb(b,i) t1_1_bb(a,j) l1_1_bb(k,b)
        //             += +0.50 P(i,j) r0_1 t1_bb(a,j) t2_abab(c,b,l,i) l2_1_abab(l,k,c,b)
        //             += +0.50 P(i,j) r0_1 t1_bb(a,j) t2_abab(b,c,l,i) l2_1_abab(l,k,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_1_bb(a,l) t1_bb(b,i) l1_1_bb(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r2_abab(b,a,l,i) l1_aa(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r2_1_abab(b,a,l,i) l1_1_aa(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r2_1_bbbb(a,b,l,i) l1_1_bb(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r2_bbbb(a,b,l,i) l1_bb(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(a,l) t1_1_bb(b,i) l1_1_bb(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(a,l) t1_bb(b,i) l1_bb(l,b)
        //             += +0.25 P(i,j) d_bb(j,k) r2_1_abab(b,c,m,l) t1_bb(a,i) l2_1_abab(m,l,b,c)
        //             += +0.25 P(i,j) d_bb(j,k) r2_1_abab(b,c,l,m) t1_bb(a,i) l2_1_abab(l,m,b,c)
        //             += +0.25 P(i,j) d_bb(j,k) r2_1_abab(c,b,m,l) t1_bb(a,i) l2_1_abab(m,l,c,b)
        //             += +0.25 P(i,j) d_bb(j,k) r2_1_abab(c,b,l,m) t1_bb(a,i) l2_1_abab(l,m,c,b)
        //             += +0.25 P(i,j) d_bb(j,k) r2_aaaa(b,c,m,l) t1_bb(a,i) l2_aaaa(m,l,b,c)
        //             += +0.25 P(i,j) d_bb(j,k) r2_bbbb(b,c,m,l) t1_bb(a,i) l2_bbbb(m,l,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_1_aa(b,l) t1_bb(a,i) l1_1_aa(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_1_bb(b,l) t1_bb(a,i) l1_1_bb(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_aa(b,l) t1_bb(a,i) l1_aa(l,b)
        //             += +0.25 P(i,j) d_bb(j,k) r2_1_bbbb(b,c,m,l) t1_bb(a,i) l2_1_bbbb(m,l,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(b,l) t1_bb(a,i) l1_bb(l,b)
        //             += +0.25 P(i,j) d_bb(j,k) r2_1_aaaa(b,c,m,l) t1_bb(a,i) l2_1_aaaa(m,l,b,c)
        //             += +0.25 P(i,j) d_bb(j,k) r2_abab(b,c,m,l) t1_bb(a,i) l2_abab(m,l,b,c)
        //             += +0.25 P(i,j) d_bb(j,k) r2_abab(b,c,l,m) t1_bb(a,i) l2_abab(l,m,b,c)
        //             += +0.25 P(i,j) d_bb(j,k) r2_abab(c,b,m,l) t1_bb(a,i) l2_abab(m,l,c,b)
        //             += +0.25 P(i,j) d_bb(j,k) r2_abab(c,b,l,m) t1_bb(a,i) l2_abab(l,m,c,b)
        //             += +1.00 P(i,j) r0 t1_bb(a,j) t1_1_bb(b,i) l1_1_bb(k,b)
        //             += +1.00 P(i,j) r0_1 t1_bb(a,j) t1_bb(b,i) l1_1_bb(k,b)
        //             += +1.00 P(i,j) d_bb(j,k) r0_1 t1_bb(a,i) l0_1
        //             += +1.00 P(i,j) r1_bb(c,i) t2_abab(b,a,l,j) l2_abab(l,k,b,c)
        //             += -1.00 P(i,j) r1_bb(c,i) t1_bb(b,j) t1_1_bb(a,l) l2_1_bbbb(l,k,b,c)
        //             += +1.00 P(i,j) r1_bb(c,i) t2_1_abab(b,a,l,j) l2_1_abab(l,k,b,c)
        //             += -1.00 P(i,j) r1_bb(c,i) t2_1_bbbb(a,b,l,j) l2_1_bbbb(l,k,b,c)
        //             += -1.00 P(i,j) r1_bb(c,i) t2_bbbb(a,b,l,j) l2_bbbb(l,k,b,c)
        //             += +1.00 P(i,j) r1_1_bb(c,i) t2_abab(b,a,l,j) l2_1_abab(l,k,b,c)
        //             += -1.00 P(i,j) r1_1_bb(c,i) t2_bbbb(a,b,l,j) l2_1_bbbb(l,k,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(b,l) t1_1_bb(a,i) l1_1_bb(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_aa(b,l) t1_1_bb(a,i) l1_1_aa(l,b)
        //             += +0.25 P(i,j) d_bb(j,k) r2_aaaa(b,c,m,l) t1_1_bb(a,i) l2_1_aaaa(m,l,b,c)
        //             += +0.25 P(i,j) d_bb(j,k) r2_bbbb(b,c,m,l) t1_1_bb(a,i) l2_1_bbbb(m,l,b,c)
        //             += +0.25 P(i,j) d_bb(j,k) r2_abab(b,c,m,l) t1_1_bb(a,i) l2_1_abab(m,l,b,c)
        //             += +0.25 P(i,j) d_bb(j,k) r2_abab(b,c,l,m) t1_1_bb(a,i) l2_1_abab(l,m,b,c)
        //             += +0.25 P(i,j) d_bb(j,k) r2_abab(c,b,m,l) t1_1_bb(a,i) l2_1_abab(m,l,c,b)
        //             += +0.25 P(i,j) d_bb(j,k) r2_abab(c,b,l,m) t1_1_bb(a,i) l2_1_abab(l,m,c,b)
        //             += -1.00 P(i,j) r0 t1_bb(b,j) t1_1_bb(a,i) l1_1_bb(k,b)
        //             += +1.00 P(i,j) r0 t1_bb(a,j) t1_bb(b,i) l1_bb(k,b)
        //             += -0.50 P(i,j) r0 t2_bbbb(c,b,l,j) t1_1_bb(a,i) l2_1_bbbb(l,k,c,b)
        //             += -0.50 P(i,j) r0 t2_abab(c,b,l,j) t1_1_bb(a,i) l2_1_abab(l,k,c,b)
        //             += -0.50 P(i,j) r0 t2_abab(b,c,l,j) t1_1_bb(a,i) l2_1_abab(l,k,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r0 t1_bb(a,i) l0
        //             += +0.50 P(i,j) r0_1 t1_bb(a,j) t2_bbbb(c,b,l,i) l2_1_bbbb(l,k,c,b)
        //             += -1.00 P(i,j) r0 t1_bb(c,j) t2_1_bbbb(a,b,l,i) l2_1_bbbb(l,k,c,b)
        //             += -1.00 P(i,j) r0 t2_bbbb(a,c,l,j) t1_1_bb(b,i) l2_1_bbbb(l,k,c,b)
        //             += -1.00 P(i,j) r0 t2_bbbb(a,b,l,i) t1_bb(c,j) l2_bbbb(l,k,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t1_bb(c,i) t2_1_abab(b,a,l,m) l2_1_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t1_bb(c,i) t2_1_abab(b,a,m,l) l2_1_abab(m,l,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_abab(c,a,l,m) t1_1_bb(b,i) l2_1_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_abab(c,a,m,l) t1_1_bb(b,i) l2_1_abab(m,l,c,b)
        //             += -1.00 P(i,j) d_bb(j,k) r0 t1_bb(b,i) t1_1_bb(a,l) l1_1_bb(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r0 t2_abab(b,a,l,i) l1_aa(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r0 t2_1_abab(b,a,l,i) l1_1_aa(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r0 t1_bb(a,l) t1_bb(b,i) l1_bb(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r0 t2_bbbb(a,b,l,i) l1_bb(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r0 t1_bb(a,l) t1_1_bb(b,i) l1_1_bb(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r0 t2_1_bbbb(a,b,l,i) l1_1_bb(l,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_abab(c,b,m,i) t1_1_bb(a,l) l2_1_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_abab(b,c,m,i) t1_1_bb(a,l) l2_1_abab(m,l,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_abab(b,a,l,m) t1_bb(c,i) l2_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_abab(b,a,m,l) t1_bb(c,i) l2_abab(m,l,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t2_bbbb(a,c,l,m) t1_1_bb(b,i) l2_1_bbbb(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_bbbb(a,b,l,m) t1_bb(c,i) l2_bbbb(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t2_bbbb(c,b,m,i) t1_1_bb(a,l) l2_1_bbbb(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t1_bb(c,i) t2_1_bbbb(a,b,l,m) l2_1_bbbb(l,m,c,b)
        //             += -1.00 P(i,j) d_bb(j,k) r0_1 t2_bbbb(a,b,l,i) l1_1_bb(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r0_1 t1_bb(a,l) t1_bb(b,i) l1_1_bb(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r0_1 t2_abab(b,a,l,i) l1_1_aa(l,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0_1 t1_bb(a,m) t2_abab(c,b,l,i) l2_1_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0_1 t1_bb(a,m) t2_abab(b,c,l,i) l2_1_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r0_1 t1_bb(a,m) t2_bbbb(c,b,l,i) l2_1_bbbb(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0_1 t2_bbbb(a,b,l,m) t1_bb(c,i) l2_1_bbbb(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0_1 t2_abab(b,a,l,m) t1_bb(c,i) l2_1_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r0_1 t2_abab(b,a,m,l) t1_bb(c,i) l2_1_abab(m,l,b,c)
        //             += -1.00 P(i,j) r1_aa(c,l) t1_bb(b,j) t1_1_bb(a,i) l2_1_abab(l,k,c,b)
        //             += +1.00 P(i,j) r1_bb(c,l) t1_bb(b,j) t1_1_bb(a,i) l2_1_bbbb(l,k,b,c)
        //             += -1.00 P(i,j) r1_bb(c,l) t1_bb(a,j) t1_1_bb(b,i) l2_1_bbbb(l,k,b,c)
        //             += +1.00 P(i,j) r1_aa(c,l) t1_bb(a,j) t1_1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t2_1_bbbb(a,b,l,i) l2_1_bbbb(m,l,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r1_1_bb(a,m) t2_abab(c,b,l,i) l2_1_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r1_1_bb(a,m) t2_abab(b,c,l,i) l2_1_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r1_1_bb(a,m) t2_bbbb(c,b,l,i) l2_1_bbbb(m,l,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r1_bb(a,m) t2_1_bbbb(c,b,l,i) l2_1_bbbb(m,l,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r1_bb(a,m) t2_bbbb(c,b,l,i) l2_bbbb(m,l,c,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t2_bbbb(a,b,l,i) l2_abab(m,l,c,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t2_1_bbbb(a,b,l,i) l2_1_abab(m,l,c,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t2_abab(b,a,l,i) l2_aaaa(m,l,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t2_bbbb(a,b,l,i) l2_bbbb(m,l,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r2_abab(b,c,m,i) t1_1_bb(a,l) l2_1_abab(m,l,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r2_abab(c,b,m,i) t1_1_bb(a,l) l2_1_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r2_bbbb(b,c,m,i) t1_1_bb(a,l) l2_1_bbbb(m,l,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(b,i) t1_1_bb(a,l) l1_1_bb(l,b)
        //             += -0.50 P(i,j) d_bb(j,k) r2_bbbb(a,c,m,l) t1_bb(b,i) l2_bbbb(m,l,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r2_1_abab(c,a,m,l) t1_bb(b,i) l2_1_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r2_1_abab(c,a,l,m) t1_bb(b,i) l2_1_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r2_1_bbbb(a,c,m,l) t1_bb(b,i) l2_1_bbbb(m,l,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r2_abab(c,a,m,l) t1_bb(b,i) l2_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r2_abab(c,a,l,m) t1_bb(b,i) l2_abab(l,m,c,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t2_1_abab(b,a,l,i) l2_1_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r1_1_bb(c,i) t2_abab(b,a,l,m) l2_1_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r1_1_bb(c,i) t2_abab(b,a,m,l) l2_1_abab(m,l,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r2_abab(c,a,m,l) t1_1_bb(b,i) l2_1_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r2_abab(c,a,l,m) t1_1_bb(b,i) l2_1_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r2_bbbb(a,c,m,l) t1_1_bb(b,i) l2_1_bbbb(m,l,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_1_bb(c,m) t2_bbbb(a,b,l,i) l2_1_bbbb(m,l,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r1_bb(c,i) t2_1_bbbb(a,b,l,m) l2_1_bbbb(l,m,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_1_aa(c,m) t2_bbbb(a,b,l,i) l2_1_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r1_bb(c,i) t2_bbbb(a,b,l,m) l2_bbbb(l,m,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_1_bb(c,m) t2_abab(b,a,l,i) l2_1_abab(l,m,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_1_aa(c,m) t2_abab(b,a,l,i) l2_1_aaaa(m,l,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r2_1_bbbb(b,c,m,i) t1_bb(a,l) l2_1_bbbb(m,l,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r2_abab(b,c,m,i) t1_bb(a,l) l2_abab(m,l,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r2_abab(c,b,m,i) t1_bb(a,l) l2_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r2_bbbb(b,c,m,i) t1_bb(a,l) l2_bbbb(m,l,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(b,i) t1_bb(a,l) l1_bb(l,b)
        //             += -0.50 P(i,j) d_bb(j,k) r2_1_abab(b,c,m,i) t1_bb(a,l) l2_1_abab(m,l,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r2_1_abab(c,b,m,i) t1_bb(a,l) l2_1_abab(m,l,c,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_1_bb(b,i) t1_bb(a,l) l1_1_bb(l,b)
        //             += +0.50 P(i,j) d_bb(j,k) r1_1_bb(c,i) t2_bbbb(a,b,l,m) l2_1_bbbb(l,m,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t2_1_abab(b,a,l,i) l2_1_aaaa(m,l,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t2_abab(b,a,l,i) l2_abab(l,m,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t1_bb(a,l) t1_bb(b,i) l2_bbbb(m,l,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t1_bb(a,l) t1_bb(b,i) l2_abab(m,l,c,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_1_aa(c,m) t1_bb(a,l) t1_bb(b,i) l2_1_abab(m,l,c,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_1_bb(c,m) t1_bb(a,l) t1_bb(b,i) l2_1_bbbb(m,l,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t1_bb(b,i) t1_1_bb(a,l) l2_1_abab(m,l,c,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t1_bb(b,i) t1_1_bb(a,l) l2_1_bbbb(m,l,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t1_bb(a,l) t1_1_bb(b,i) l2_1_bbbb(m,l,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t1_bb(a,l) t1_1_bb(b,i) l2_1_abab(m,l,c,b)
        //             += -1.00 P(i,j) r1_bb(a,i) t1_bb(b,j) l1_bb(k,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_1_bb(a,i) l0_1
        //             += -0.50 P(i,j) d_bb(j,k) r1_bb(c,i) t2_1_abab(b,a,l,m) l2_1_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r1_bb(c,i) t2_1_abab(b,a,m,l) l2_1_abab(m,l,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r1_bb(c,i) t2_abab(b,a,l,m) l2_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r1_bb(c,i) t2_abab(b,a,m,l) l2_abab(m,l,b,c)
        D_ooov_bbbb("L,R,i,j,k,a") -= tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k");
        D_ooov_bbbb("L,R,i,j,k,a") += tmps_["231_bbbb_LLvooo"]("R,L,a,j,i,k");

        // D_oovo_bbbb += +0.50 P(i,j) r1_1_bb(a,i) t2_abab(c,b,l,j) l2_1_abab(l,k,c,b)
        //             += +0.50 P(i,j) r1_1_bb(a,i) t2_abab(b,c,l,j) l2_1_abab(l,k,b,c)
        //             += +1.00 P(i,j) r0_1 t2_bbbb(a,b,l,i) t1_bb(c,j) l2_1_bbbb(l,k,c,b)
        //             += +0.50 P(i,j) r1_1_bb(a,i) t2_bbbb(c,b,l,j) l2_1_bbbb(l,k,c,b)
        //             += +1.00 P(i,j) r0_1 t2_abab(b,a,l,i) t1_bb(c,j) l2_1_abab(l,k,b,c)
        //             += +1.00 P(i,j) r0 t1_bb(a,l) t1_bb(c,j) t1_1_bb(b,i) l2_1_bbbb(l,k,c,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(a,i) l0
        //             += +1.00 P(i,j) r1_bb(a,i) t1_1_bb(b,j) l1_1_bb(k,b)
        //             += +1.00 P(i,j) r1_1_bb(a,i) t1_bb(b,j) l1_1_bb(k,b)
        //             += +1.00 P(i,j) r1_bb(c,l) t1_bb(a,j) t1_bb(b,i) l2_bbbb(l,k,b,c)
        //             += -1.00 P(i,j) r1_aa(c,l) t1_bb(a,j) t1_bb(b,i) l2_abab(l,k,c,b)
        //             += -1.00 P(i,j) r1_1_aa(c,l) t1_bb(a,j) t1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += +1.00 P(i,j) r1_1_bb(c,l) t1_bb(a,j) t1_bb(b,i) l2_1_bbbb(l,k,b,c)
        //             += +1.00 P(i,j) r0 t1_bb(c,j) t2_1_abab(b,a,l,i) l2_1_abab(l,k,b,c)
        //             += -1.00 P(i,j) r0 t2_abab(c,a,l,j) t1_1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += +1.00 P(i,j) r0 t2_abab(b,a,l,i) t1_bb(c,j) l2_abab(l,k,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r0 t1_1_bb(a,i) l0_1
        //             += -0.50 P(i,j) r2_1_bbbb(b,c,l,i) t1_bb(a,j) l2_1_bbbb(l,k,b,c)
        //             += -0.50 P(i,j) r2_abab(b,c,l,i) t1_bb(a,j) l2_abab(l,k,b,c)
        //             += -0.50 P(i,j) r2_abab(c,b,l,i) t1_bb(a,j) l2_abab(l,k,c,b)
        //             += -0.50 P(i,j) r2_bbbb(b,c,l,i) t1_bb(a,j) l2_bbbb(l,k,b,c)
        //             += -1.00 P(i,j) r1_bb(b,i) t1_bb(a,j) l1_bb(k,b)
        //             += -0.50 P(i,j) r2_1_abab(b,c,l,i) t1_bb(a,j) l2_1_abab(l,k,b,c)
        //             += -0.50 P(i,j) r2_1_abab(c,b,l,i) t1_bb(a,j) l2_1_abab(l,k,c,b)
        //             += -1.00 P(i,j) r1_1_bb(b,i) t1_bb(a,j) l1_1_bb(k,b)
        //             += -0.50 P(i,j) r2_abab(b,c,l,i) t1_1_bb(a,j) l2_1_abab(l,k,b,c)
        //             += -0.50 P(i,j) r2_abab(c,b,l,i) t1_1_bb(a,j) l2_1_abab(l,k,c,b)
        //             += -0.50 P(i,j) r2_bbbb(b,c,l,i) t1_1_bb(a,j) l2_1_bbbb(l,k,b,c)
        //             += -1.00 P(i,j) r1_bb(b,i) t1_1_bb(a,j) l1_1_bb(k,b)
        //             += -0.50 P(i,j) r0_1 t1_bb(a,j) t2_abab(c,b,l,i) l2_1_abab(l,k,c,b)
        //             += -0.50 P(i,j) r0_1 t1_bb(a,j) t2_abab(b,c,l,i) l2_1_abab(l,k,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_1_bb(a,l) t1_bb(b,i) l1_1_bb(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r2_abab(b,a,l,i) l1_aa(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r2_1_abab(b,a,l,i) l1_1_aa(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r2_1_bbbb(a,b,l,i) l1_1_bb(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r2_bbbb(a,b,l,i) l1_bb(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(a,l) t1_1_bb(b,i) l1_1_bb(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(a,l) t1_bb(b,i) l1_bb(l,b)
        //             += -0.25 P(i,j) d_bb(j,k) r2_1_abab(b,c,m,l) t1_bb(a,i) l2_1_abab(m,l,b,c)
        //             += -0.25 P(i,j) d_bb(j,k) r2_1_abab(b,c,l,m) t1_bb(a,i) l2_1_abab(l,m,b,c)
        //             += -0.25 P(i,j) d_bb(j,k) r2_1_abab(c,b,m,l) t1_bb(a,i) l2_1_abab(m,l,c,b)
        //             += -0.25 P(i,j) d_bb(j,k) r2_1_abab(c,b,l,m) t1_bb(a,i) l2_1_abab(l,m,c,b)
        //             += -0.25 P(i,j) d_bb(j,k) r2_aaaa(b,c,m,l) t1_bb(a,i) l2_aaaa(m,l,b,c)
        //             += -0.25 P(i,j) d_bb(j,k) r2_bbbb(b,c,m,l) t1_bb(a,i) l2_bbbb(m,l,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_1_aa(b,l) t1_bb(a,i) l1_1_aa(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_1_bb(b,l) t1_bb(a,i) l1_1_bb(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_aa(b,l) t1_bb(a,i) l1_aa(l,b)
        //             += -0.25 P(i,j) d_bb(j,k) r2_1_bbbb(b,c,m,l) t1_bb(a,i) l2_1_bbbb(m,l,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(b,l) t1_bb(a,i) l1_bb(l,b)
        //             += -0.25 P(i,j) d_bb(j,k) r2_1_aaaa(b,c,m,l) t1_bb(a,i) l2_1_aaaa(m,l,b,c)
        //             += -0.25 P(i,j) d_bb(j,k) r2_abab(b,c,m,l) t1_bb(a,i) l2_abab(m,l,b,c)
        //             += -0.25 P(i,j) d_bb(j,k) r2_abab(b,c,l,m) t1_bb(a,i) l2_abab(l,m,b,c)
        //             += -0.25 P(i,j) d_bb(j,k) r2_abab(c,b,m,l) t1_bb(a,i) l2_abab(m,l,c,b)
        //             += -0.25 P(i,j) d_bb(j,k) r2_abab(c,b,l,m) t1_bb(a,i) l2_abab(l,m,c,b)
        //             += -1.00 P(i,j) r0 t1_bb(a,j) t1_1_bb(b,i) l1_1_bb(k,b)
        //             += -1.00 P(i,j) r0_1 t1_bb(a,j) t1_bb(b,i) l1_1_bb(k,b)
        //             += -1.00 P(i,j) d_bb(j,k) r0_1 t1_bb(a,i) l0_1
        //             += -1.00 P(i,j) r1_bb(c,i) t2_abab(b,a,l,j) l2_abab(l,k,b,c)
        //             += +1.00 P(i,j) r1_bb(c,i) t1_bb(b,j) t1_1_bb(a,l) l2_1_bbbb(l,k,b,c)
        //             += -1.00 P(i,j) r1_bb(c,i) t2_1_abab(b,a,l,j) l2_1_abab(l,k,b,c)
        //             += +1.00 P(i,j) r1_bb(c,i) t2_1_bbbb(a,b,l,j) l2_1_bbbb(l,k,b,c)
        //             += +1.00 P(i,j) r1_bb(c,i) t2_bbbb(a,b,l,j) l2_bbbb(l,k,b,c)
        //             += -1.00 P(i,j) r1_1_bb(c,i) t2_abab(b,a,l,j) l2_1_abab(l,k,b,c)
        //             += +1.00 P(i,j) r1_1_bb(c,i) t2_bbbb(a,b,l,j) l2_1_bbbb(l,k,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(b,l) t1_1_bb(a,i) l1_1_bb(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_aa(b,l) t1_1_bb(a,i) l1_1_aa(l,b)
        //             += -0.25 P(i,j) d_bb(j,k) r2_aaaa(b,c,m,l) t1_1_bb(a,i) l2_1_aaaa(m,l,b,c)
        //             += -0.25 P(i,j) d_bb(j,k) r2_bbbb(b,c,m,l) t1_1_bb(a,i) l2_1_bbbb(m,l,b,c)
        //             += -0.25 P(i,j) d_bb(j,k) r2_abab(b,c,m,l) t1_1_bb(a,i) l2_1_abab(m,l,b,c)
        //             += -0.25 P(i,j) d_bb(j,k) r2_abab(b,c,l,m) t1_1_bb(a,i) l2_1_abab(l,m,b,c)
        //             += -0.25 P(i,j) d_bb(j,k) r2_abab(c,b,m,l) t1_1_bb(a,i) l2_1_abab(m,l,c,b)
        //             += -0.25 P(i,j) d_bb(j,k) r2_abab(c,b,l,m) t1_1_bb(a,i) l2_1_abab(l,m,c,b)
        //             += +1.00 P(i,j) r0 t1_bb(b,j) t1_1_bb(a,i) l1_1_bb(k,b)
        //             += -1.00 P(i,j) r0 t1_bb(a,j) t1_bb(b,i) l1_bb(k,b)
        //             += +0.50 P(i,j) r0 t2_bbbb(c,b,l,j) t1_1_bb(a,i) l2_1_bbbb(l,k,c,b)
        //             += +0.50 P(i,j) r0 t2_abab(c,b,l,j) t1_1_bb(a,i) l2_1_abab(l,k,c,b)
        //             += +0.50 P(i,j) r0 t2_abab(b,c,l,j) t1_1_bb(a,i) l2_1_abab(l,k,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r0 t1_bb(a,i) l0
        //             += -0.50 P(i,j) r0_1 t1_bb(a,j) t2_bbbb(c,b,l,i) l2_1_bbbb(l,k,c,b)
        //             += +1.00 P(i,j) r0 t1_bb(c,j) t2_1_bbbb(a,b,l,i) l2_1_bbbb(l,k,c,b)
        //             += +1.00 P(i,j) r0 t2_bbbb(a,c,l,j) t1_1_bb(b,i) l2_1_bbbb(l,k,c,b)
        //             += +1.00 P(i,j) r0 t2_bbbb(a,b,l,i) t1_bb(c,j) l2_bbbb(l,k,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t1_bb(c,i) t2_1_abab(b,a,l,m) l2_1_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t1_bb(c,i) t2_1_abab(b,a,m,l) l2_1_abab(m,l,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t2_abab(c,a,l,m) t1_1_bb(b,i) l2_1_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t2_abab(c,a,m,l) t1_1_bb(b,i) l2_1_abab(m,l,c,b)
        //             += +1.00 P(i,j) d_bb(j,k) r0 t1_bb(b,i) t1_1_bb(a,l) l1_1_bb(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r0 t2_abab(b,a,l,i) l1_aa(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r0 t2_1_abab(b,a,l,i) l1_1_aa(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r0 t1_bb(a,l) t1_bb(b,i) l1_bb(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r0 t2_bbbb(a,b,l,i) l1_bb(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r0 t1_bb(a,l) t1_1_bb(b,i) l1_1_bb(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r0 t2_1_bbbb(a,b,l,i) l1_1_bb(l,b)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t2_abab(c,b,m,i) t1_1_bb(a,l) l2_1_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t2_abab(b,c,m,i) t1_1_bb(a,l) l2_1_abab(m,l,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t2_abab(b,a,l,m) t1_bb(c,i) l2_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t2_abab(b,a,m,l) t1_bb(c,i) l2_abab(m,l,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_bbbb(a,c,l,m) t1_1_bb(b,i) l2_1_bbbb(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t2_bbbb(a,b,l,m) t1_bb(c,i) l2_bbbb(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_bbbb(c,b,m,i) t1_1_bb(a,l) l2_1_bbbb(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t1_bb(c,i) t2_1_bbbb(a,b,l,m) l2_1_bbbb(l,m,c,b)
        //             += +1.00 P(i,j) d_bb(j,k) r0_1 t2_bbbb(a,b,l,i) l1_1_bb(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r0_1 t1_bb(a,l) t1_bb(b,i) l1_1_bb(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r0_1 t2_abab(b,a,l,i) l1_1_aa(l,b)
        //             += +0.50 P(i,j) d_bb(j,k) r0_1 t1_bb(a,m) t2_abab(c,b,l,i) l2_1_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r0_1 t1_bb(a,m) t2_abab(b,c,l,i) l2_1_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r0_1 t1_bb(a,m) t2_bbbb(c,b,l,i) l2_1_bbbb(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r0_1 t2_bbbb(a,b,l,m) t1_bb(c,i) l2_1_bbbb(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r0_1 t2_abab(b,a,l,m) t1_bb(c,i) l2_1_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r0_1 t2_abab(b,a,m,l) t1_bb(c,i) l2_1_abab(m,l,b,c)
        //             += +1.00 P(i,j) r1_aa(c,l) t1_bb(b,j) t1_1_bb(a,i) l2_1_abab(l,k,c,b)
        //             += -1.00 P(i,j) r1_bb(c,l) t1_bb(b,j) t1_1_bb(a,i) l2_1_bbbb(l,k,b,c)
        //             += +1.00 P(i,j) r1_bb(c,l) t1_bb(a,j) t1_1_bb(b,i) l2_1_bbbb(l,k,b,c)
        //             += -1.00 P(i,j) r1_aa(c,l) t1_bb(a,j) t1_1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t2_1_bbbb(a,b,l,i) l2_1_bbbb(m,l,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r1_1_bb(a,m) t2_abab(c,b,l,i) l2_1_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r1_1_bb(a,m) t2_abab(b,c,l,i) l2_1_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r1_1_bb(a,m) t2_bbbb(c,b,l,i) l2_1_bbbb(m,l,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r1_bb(a,m) t2_1_bbbb(c,b,l,i) l2_1_bbbb(m,l,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r1_bb(a,m) t2_bbbb(c,b,l,i) l2_bbbb(m,l,c,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t2_bbbb(a,b,l,i) l2_abab(m,l,c,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t2_1_bbbb(a,b,l,i) l2_1_abab(m,l,c,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t2_abab(b,a,l,i) l2_aaaa(m,l,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t2_bbbb(a,b,l,i) l2_bbbb(m,l,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r2_abab(b,c,m,i) t1_1_bb(a,l) l2_1_abab(m,l,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r2_abab(c,b,m,i) t1_1_bb(a,l) l2_1_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r2_bbbb(b,c,m,i) t1_1_bb(a,l) l2_1_bbbb(m,l,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(b,i) t1_1_bb(a,l) l1_1_bb(l,b)
        //             += +0.50 P(i,j) d_bb(j,k) r2_bbbb(a,c,m,l) t1_bb(b,i) l2_bbbb(m,l,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r2_1_abab(c,a,m,l) t1_bb(b,i) l2_1_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r2_1_abab(c,a,l,m) t1_bb(b,i) l2_1_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r2_1_bbbb(a,c,m,l) t1_bb(b,i) l2_1_bbbb(m,l,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r2_abab(c,a,m,l) t1_bb(b,i) l2_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r2_abab(c,a,l,m) t1_bb(b,i) l2_abab(l,m,c,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t2_1_abab(b,a,l,i) l2_1_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r1_1_bb(c,i) t2_abab(b,a,l,m) l2_1_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r1_1_bb(c,i) t2_abab(b,a,m,l) l2_1_abab(m,l,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r2_abab(c,a,m,l) t1_1_bb(b,i) l2_1_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r2_abab(c,a,l,m) t1_1_bb(b,i) l2_1_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r2_bbbb(a,c,m,l) t1_1_bb(b,i) l2_1_bbbb(m,l,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_1_bb(c,m) t2_bbbb(a,b,l,i) l2_1_bbbb(m,l,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r1_bb(c,i) t2_1_bbbb(a,b,l,m) l2_1_bbbb(l,m,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_1_aa(c,m) t2_bbbb(a,b,l,i) l2_1_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r1_bb(c,i) t2_bbbb(a,b,l,m) l2_bbbb(l,m,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_1_bb(c,m) t2_abab(b,a,l,i) l2_1_abab(l,m,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_1_aa(c,m) t2_abab(b,a,l,i) l2_1_aaaa(m,l,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r2_1_bbbb(b,c,m,i) t1_bb(a,l) l2_1_bbbb(m,l,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r2_abab(b,c,m,i) t1_bb(a,l) l2_abab(m,l,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r2_abab(c,b,m,i) t1_bb(a,l) l2_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r2_bbbb(b,c,m,i) t1_bb(a,l) l2_bbbb(m,l,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(b,i) t1_bb(a,l) l1_bb(l,b)
        //             += +0.50 P(i,j) d_bb(j,k) r2_1_abab(b,c,m,i) t1_bb(a,l) l2_1_abab(m,l,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r2_1_abab(c,b,m,i) t1_bb(a,l) l2_1_abab(m,l,c,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_1_bb(b,i) t1_bb(a,l) l1_1_bb(l,b)
        //             += -0.50 P(i,j) d_bb(j,k) r1_1_bb(c,i) t2_bbbb(a,b,l,m) l2_1_bbbb(l,m,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t2_1_abab(b,a,l,i) l2_1_aaaa(m,l,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t2_abab(b,a,l,i) l2_abab(l,m,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t1_bb(a,l) t1_bb(b,i) l2_bbbb(m,l,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t1_bb(a,l) t1_bb(b,i) l2_abab(m,l,c,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_1_aa(c,m) t1_bb(a,l) t1_bb(b,i) l2_1_abab(m,l,c,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_1_bb(c,m) t1_bb(a,l) t1_bb(b,i) l2_1_bbbb(m,l,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t1_bb(b,i) t1_1_bb(a,l) l2_1_abab(m,l,c,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t1_bb(b,i) t1_1_bb(a,l) l2_1_bbbb(m,l,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t1_bb(a,l) t1_1_bb(b,i) l2_1_bbbb(m,l,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t1_bb(a,l) t1_1_bb(b,i) l2_1_abab(m,l,c,b)
        //             += +1.00 P(i,j) r1_bb(a,i) t1_bb(b,j) l1_bb(k,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_1_bb(a,i) l0_1
        //             += +0.50 P(i,j) d_bb(j,k) r1_bb(c,i) t2_1_abab(b,a,l,m) l2_1_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r1_bb(c,i) t2_1_abab(b,a,m,l) l2_1_abab(m,l,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r1_bb(c,i) t2_abab(b,a,l,m) l2_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r1_bb(c,i) t2_abab(b,a,m,l) l2_abab(m,l,b,c)
        D_oovo_bbbb("L,R,i,j,a,k") += tmps_["231_bbbb_LLvooo"]("R,L,a,i,j,k");
        D_oovo_bbbb("L,R,i,j,a,k") -= tmps_["231_bbbb_LLvooo"]("R,L,a,j,i,k");
        tmps_["231_bbbb_LLvooo"].~TArrayD();

        // flops: o2v1L1  = o3v2L1 o3v2L1 o3v0L1 o3v1L1
        //  mems: o2v1L1  = o2v0L1 o2v0L1 o3v0L1 o2v1L1
        tmps_["232_bbb_Loov"]("L,i,j,a")  = (t2["bbbb"]("b,c,j,i") * l2["bbbb"]("L,j,k,b,c") + -1.00 * t2_1["bbbb"]("c,a,k,i") * l2_1["bbbb"]("L,k,j,c,a")) * t1["bb"]("a,k");

        // D_oovo_bbbb += +0.50 P(i,j) d_bb(j,k) r0 t1_bb(a,m) t2_bbbb(c,b,l,i) l2_bbbb(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_1_bbbb(b,a,m,i) l2_1_bbbb(m,l,b,a)
        // flops: o3v1L2 += o4v1L2
        //  mems: o3v1L2 += o4v1L2
        tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k")  = 0.50 * tmps_["232_bbb_Loov"]("L,i,l,a") * tmps_["112_bb_Loo"]("R,j,k");
        D_oovo_bbbb("L,R,i,j,a,k") += tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k");
        D_oovo_bbbb("L,R,i,j,a,k") -= tmps_["perm_bbbb_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_bbbb_LLvooo"].~TArrayD();

        // flops: o4v0L1  = o3v2L1 o3v2L1 o4v0L1
        //  mems: o4v0L1  = o2v0L1 o2v0L1 o4v0L1
        tmps_["233_bbbb_Loooo"]("L,i,j,k,l")  = t2["bbbb"]("a,b,m,i") * l2["bbbb"]("L,m,j,a,b");
        tmps_["233_bbbb_Loooo"]("L,i,j,k,l") += t2_1["bbbb"]("a,b,m,k") * l2_1["bbbb"]("L,m,l,a,b");

        // D_ooov_bbbb += +0.50 P(i,j) r0 t1_bb(a,j) t2_bbbb(c,b,l,i) l2_bbbb(l,k,c,b)
        //             += +0.50 r0 t2_abab(a,b,i,l) t2_1_bbbb(d,c,k,j) l2_1_bbbb(k,l,d,c)
        // flops: o3v1L2 += o4v1L2
        //  mems: o3v1L2 += o3v1L2
        tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k")  = 0.50 * tmps_["233_bbbb_Loooo"]("L,i,k,j,l") * tmps_["224_bb_Lvo"]("R,a,j");
        D_ooov_bbbb("L,R,i,j,k,a") += tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k");
        D_ooov_bbbb("L,R,i,j,k,a") -= tmps_["perm_bbbb_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_bbbb_LLvooo"].~TArrayD();

    }
} // hilbert
#endif