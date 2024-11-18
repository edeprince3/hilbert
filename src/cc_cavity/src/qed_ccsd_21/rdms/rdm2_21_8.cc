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


    void EOM_EE_QED_RDM_21::rdm2_21_8() {

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


        // flops: o1v1L2  = o1v2L2 o1v2L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2
        tmps_["308_aa_LLov"]("R,L,i,a")  = tmps_["29_aa_Lvv"]("L,b,a") * r1_1["aa"]("R,b,i");
        tmps_["308_aa_LLov"]("R,L,i,a") += tmps_["249_aa_Lvv"]("L,b,a") * r1["aa"]("R,b,i");
        tmps_["249_aa_Lvv"].~TArrayD();
        tmps_["29_aa_Lvv"].~TArrayD();

        // D_oovv_aaaa += +0.50 P(i,j) P(a,b) r1_aa(d,i) t1_aa(b,j) t2_1_abab(a,c,k,l) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(d,i) t1_aa(b,j) t2_1_abab(a,c,l,k) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_1_aa(d,i) t2_abab(a,c,k,l) t1_aa(b,j) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_1_aa(d,i) t2_abab(a,c,l,k) t1_aa(b,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("b,j") * tmps_["308_aa_LLov"]("R,L,i,a");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r1_aa(d,i) t1_bb(b,j) t2_1_abab(a,c,k,l) l2_1_abab(k,l,d,c)
        //             += +0.50 r1_aa(d,i) t1_bb(b,j) t2_1_abab(a,c,l,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r1_1_aa(d,i) t2_abab(a,c,k,l) t1_bb(b,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r1_1_aa(d,i) t2_abab(a,c,l,k) t1_bb(b,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1["bb"]("b,j") * tmps_["308_aa_LLov"]("R,L,i,a");

        // flops: o2v0L1  = o2v0L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["293_aa_Loo"]("L,i,j")  = Id["aa_oo"]("i,j") * l0("L");

        // flops: o2v0L1  = o2v0L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["292_aa_Loo"]("L,i,j")  = Id["aa_oo"]("i,j") * l0_1("L");

        // flops: o3v1L2  = o2v0L1 o2v0L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o1v1L2 o1v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v3L1 o2v2L1 o3v2L2 o3v3L1 o3v2L2 o3v1L2 o3v2L2 o3v1L2 o4v1L2 o4v1L2 o3v1L2 o3v2L2 o3v1L2 o3v3L1 o3v2L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        //  mems: o3v1L2  = o2v0L1 o2v0L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o1v1L2 o1v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o2v2L1 o2v2L1 o3v1L2 o2v2L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o4v0L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o2v2L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k")  = (-1.00 * tmps_["293_aa_Loo"]("L,i,j") + tmps_["133_aa_Loo"]("L,i,j") + tmps_["134_aa_Loo"]("L,i,j")) * r1["aa"]("R,a,k");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") -= r1_1["aa"]("R,a,k") * tmps_["292_aa_Loo"]("L,i,j");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += r1_1["aa"]("R,a,k") * tmps_["63_aa_Loo"]("L,i,j");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += r0_1("R") * tmps_["294_aaaa_Lvooo"]("L,a,k,j,i");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += tmps_["307_aa_LLvo"]("L,R,a,k") * Id["aa_oo"]("i,j");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += tmps_["297_aaaa_Lvooo"]("L,a,k,j,i") * r0("R");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") -= t1["aa"]("a,i") * tmps_["137_aa_LLoo"]("L,R,j,k");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += tmps_["305_aa_Lov"]("L,k,a") * tmps_["136_aa_Loo"]("R,i,j");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += t1_1["aa"]("a,k") * tmps_["121_aa_LLoo"]("L,R,j,i");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") -= (t1["aa"]("a,k") * tmps_["229_LL"]("L,R") + -1.00 * tmps_["302_aa_LLvo"]("R,L,a,k")) * Id["aa_oo"]("i,j");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") -= 0.50 * tmps_["131_aa_Loo"]("L,k,j") * tmps_["254_aa_Lvo"]("R,a,i");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") -= t1_1["aa"]("a,i") * tmps_["141_aa_LLoo"]("L,R,j,k");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") -= tmps_["292_aa_Loo"]("L,i,j") * tmps_["255_aa_Lvo"]("R,a,k");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += tmps_["63_aa_Loo"]("L,i,j") * tmps_["255_aa_Lvo"]("R,a,k");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") -= 0.50 * t1["aa"]("a,i") * tmps_["135_aa_LLoo"]("L,R,j,k");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") -= tmps_["293_aa_Loo"]("L,i,j") * tmps_["253_aa_Lvo"]("R,a,k");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") -= tmps_["63_aa_Loo"]("L,k,j") * tmps_["254_aa_Lvo"]("R,a,i");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += 0.50 * tmps_["131_aa_Loo"]("L,i,j") * tmps_["255_aa_Lvo"]("R,a,k");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += tmps_["132_aa_Loo"]("L,i,j") * tmps_["255_aa_Lvo"]("R,a,k");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += (t2["aaaa"]("a,b,l,i") * l2_1["aaaa"]("L,l,j,b,c") + -1.00 * tmps_["243_aaaa_Lvoov"]("L,a,i,j,c")) * r1_1["aa"]("R,c,k");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += t2_1["aaaa"]("a,b,l,i") * l2_1["aaaa"]("L,l,j,b,c") * r1["aa"]("R,c,k");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") -= tmps_["245_aaaa_Lvoov"]("L,a,i,j,c") * r1["aa"]("R,c,k");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += tmps_["22_aaaa_Looov"]("L,i,l,j,c") * r1["aa"]("R,c,k") * t1_1["aa"]("a,l");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") -= tmps_["244_aaaa_Lvoov"]("L,a,i,j,c") * r1["aa"]("R,c,k");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += t2["aaaa"]("a,b,l,i") * l2["aaaa"]("L,l,j,b,c") * r1["aa"]("R,c,k");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") -= tmps_["134_aa_Loo"]("L,k,j") * tmps_["253_aa_Lvo"]("R,a,i");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") -= tmps_["132_aa_Loo"]("L,k,j") * tmps_["254_aa_Lvo"]("R,a,i");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") -= tmps_["133_aa_Loo"]("L,k,j") * tmps_["253_aa_Lvo"]("R,a,i");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") -= tmps_["292_aa_Loo"]("L,i,j") * tmps_["254_aa_Lvo"]("R,a,k");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") -= Id["aa_oo"]("i,j") * t1_1["aa"]("a,k") * tmps_["230_LL"]("L,R");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += tmps_["304_aa_LLvo"]("L,R,a,k") * Id["aa_oo"]("i,j");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += 0.50 * tmps_["306_aa_Lov"]("L,k,a") * tmps_["103_aa_Loo"]("R,i,j");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += t1["aa"]("a,i") * tmps_["138_aa_LLoo"]("L,R,j,k");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += tmps_["298_aaaa_Lvooo"]("L,a,k,j,i") * r0("R");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += r0_1("R") * tmps_["295_aaaa_Lvooo"]("L,a,k,j,i");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += 0.50 * r1_1["aa"]("R,a,k") * tmps_["131_aa_Loo"]("L,i,j");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += r0("R") * tmps_["296_aaaa_Looov"]("L,k,i,j,a");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += r1_1["aa"]("R,a,k") * tmps_["132_aa_Loo"]("L,i,j");
        tmps_["309_aaaa_LLoovo"]("L,R,i,j,a,k") += Id["aa_oo"]("i,j") * tmps_["308_aa_LLov"]("R,L,k,a");
        tmps_["298_aaaa_Lvooo"].~TArrayD();
        tmps_["297_aaaa_Lvooo"].~TArrayD();
        tmps_["296_aaaa_Looov"].~TArrayD();
        tmps_["295_aaaa_Lvooo"].~TArrayD();
        tmps_["294_aaaa_Lvooo"].~TArrayD();
        tmps_["293_aa_Loo"].~TArrayD();
        tmps_["292_aa_Loo"].~TArrayD();
        tmps_["245_aaaa_Lvoov"].~TArrayD();
        tmps_["244_aaaa_Lvoov"].~TArrayD();
        tmps_["243_aaaa_Lvoov"].~TArrayD();
        tmps_["141_aa_LLoo"].~TArrayD();
        tmps_["138_aa_LLoo"].~TArrayD();
        tmps_["137_aa_LLoo"].~TArrayD();
        tmps_["136_aa_Loo"].~TArrayD();
        tmps_["135_aa_LLoo"].~TArrayD();
        tmps_["134_aa_Loo"].~TArrayD();
        tmps_["133_aa_Loo"].~TArrayD();
        tmps_["132_aa_Loo"].~TArrayD();
        tmps_["131_aa_Loo"].~TArrayD();
        tmps_["121_aa_LLoo"].~TArrayD();
        tmps_["103_aa_Loo"].~TArrayD();
        tmps_["63_aa_Loo"].~TArrayD();
        tmps_["22_aaaa_Looov"].~TArrayD();

        // D_ooov_aaaa += -1.00 P(i,j) r1_aa(a,i) t1_aa(b,j) l1_aa(k,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(a,i) l0
        //             += -1.00 P(i,j) r1_aa(a,i) t1_1_aa(b,j) l1_1_aa(k,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_1_aa(a,i) l0_1
        //             += -0.50 P(i,j) r1_1_aa(a,i) t2_abab(c,b,j,l) l2_1_abab(k,l,c,b)
        //             += -0.50 P(i,j) r1_1_aa(a,i) t2_abab(b,c,j,l) l2_1_abab(k,l,b,c)
        //             += -1.00 P(i,j) r0_1 t2_abab(a,b,i,l) t1_aa(c,j) l2_1_abab(k,l,c,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t1_aa(b,i) t1_1_aa(a,l) l2_1_abab(l,m,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t1_aa(b,i) t1_1_aa(a,l) l2_1_aaaa(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t1_aa(a,l) t1_aa(b,i) l2_abab(l,m,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_1_bb(c,m) t1_aa(a,l) t1_aa(b,i) l2_1_abab(l,m,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t1_aa(a,l) t1_aa(b,i) l2_aaaa(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_1_aa(c,m) t1_aa(a,l) t1_aa(b,i) l2_1_aaaa(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t1_aa(a,l) t1_1_aa(b,i) l2_1_aaaa(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t1_aa(a,l) t1_1_aa(b,i) l2_1_abab(l,m,b,c)
        //             += -1.00 P(i,j) r0 t1_aa(c,j) t2_1_abab(a,b,i,l) l2_1_abab(k,l,c,b)
        //             += -1.00 P(i,j) r0 t2_abab(a,b,i,l) t1_aa(c,j) l2_abab(k,l,c,b)
        //             += +1.00 P(i,j) r0 t2_abab(a,c,j,l) t1_1_aa(b,i) l2_1_abab(k,l,b,c)
        //             += +1.00 P(i,j) r1_bb(c,l) t1_aa(a,j) t1_aa(b,i) l2_abab(k,l,b,c)
        //             += +1.00 P(i,j) r1_1_bb(c,l) t1_aa(a,j) t1_aa(b,i) l2_1_abab(k,l,b,c)
        //             += -1.00 P(i,j) r1_aa(c,l) t1_aa(a,j) t1_aa(b,i) l2_aaaa(l,k,b,c)
        //             += -1.00 P(i,j) r1_1_aa(c,l) t1_aa(a,j) t1_aa(b,i) l2_1_aaaa(l,k,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r0_1 t2_abab(a,b,l,m) t1_aa(c,i) l2_1_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0_1 t2_abab(a,b,m,l) t1_aa(c,i) l2_1_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0_1 t1_aa(a,m) t2_aaaa(c,b,l,i) l2_1_aaaa(l,m,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0_1 t1_aa(a,m) t2_abab(c,b,i,l) l2_1_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0_1 t1_aa(a,m) t2_abab(b,c,i,l) l2_1_abab(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r0_1 t1_aa(a,l) t1_aa(b,i) l1_1_aa(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r0_1 t2_aaaa(a,b,l,i) l1_1_aa(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r0_1 t2_abab(a,b,i,l) l1_1_bb(l,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0_1 t2_aaaa(a,b,l,m) t1_aa(c,i) l2_1_aaaa(l,m,c,b)
        //             += -1.00 P(i,j) r1_bb(c,l) t1_aa(b,j) t1_1_aa(a,i) l2_1_abab(k,l,b,c)
        //             += +1.00 P(i,j) r1_aa(c,l) t1_aa(b,j) t1_1_aa(a,i) l2_1_aaaa(l,k,b,c)
        //             += +0.25 P(i,j) d_aa(j,k) r2_1_abab(b,c,m,l) t1_aa(a,i) l2_1_abab(m,l,b,c)
        //             += +0.25 P(i,j) d_aa(j,k) r2_1_abab(b,c,l,m) t1_aa(a,i) l2_1_abab(l,m,b,c)
        //             += +0.25 P(i,j) d_aa(j,k) r2_1_abab(c,b,m,l) t1_aa(a,i) l2_1_abab(m,l,c,b)
        //             += +0.25 P(i,j) d_aa(j,k) r2_1_abab(c,b,l,m) t1_aa(a,i) l2_1_abab(l,m,c,b)
        //             += +0.25 P(i,j) d_aa(j,k) r2_aaaa(b,c,m,l) t1_aa(a,i) l2_aaaa(m,l,b,c)
        //             += +0.25 P(i,j) d_aa(j,k) r2_bbbb(b,c,m,l) t1_aa(a,i) l2_bbbb(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_1_aa(b,l) t1_aa(a,i) l1_1_aa(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_1_bb(b,l) t1_aa(a,i) l1_1_bb(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(b,l) t1_aa(a,i) l1_aa(l,b)
        //             += +0.25 P(i,j) d_aa(j,k) r2_1_bbbb(b,c,m,l) t1_aa(a,i) l2_1_bbbb(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_bb(b,l) t1_aa(a,i) l1_bb(l,b)
        //             += +0.25 P(i,j) d_aa(j,k) r2_1_aaaa(b,c,m,l) t1_aa(a,i) l2_1_aaaa(m,l,b,c)
        //             += +0.25 P(i,j) d_aa(j,k) r2_abab(b,c,m,l) t1_aa(a,i) l2_abab(m,l,b,c)
        //             += +0.25 P(i,j) d_aa(j,k) r2_abab(b,c,l,m) t1_aa(a,i) l2_abab(l,m,b,c)
        //             += +0.25 P(i,j) d_aa(j,k) r2_abab(c,b,m,l) t1_aa(a,i) l2_abab(m,l,c,b)
        //             += +0.25 P(i,j) d_aa(j,k) r2_abab(c,b,l,m) t1_aa(a,i) l2_abab(l,m,c,b)
        //             += -1.00 P(i,j) d_aa(j,k) r2_1_aaaa(a,b,l,i) l1_1_aa(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(a,l) t1_1_aa(b,i) l1_1_aa(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r2_aaaa(a,b,l,i) l1_aa(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t2_1_abab(a,b,i,l) l2_1_bbbb(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(a,l) t1_aa(b,i) l1_aa(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r2_1_abab(a,b,i,l) l1_1_bb(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r2_abab(a,b,i,l) l1_bb(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_1_aa(a,l) t1_aa(b,i) l1_1_aa(l,b)
        //             += +0.50 P(i,j) r0_1 t1_aa(a,j) t2_aaaa(c,b,l,i) l2_1_aaaa(l,k,c,b)
        //             += +0.50 P(i,j) r2_abab(b,c,i,l) t1_1_aa(a,j) l2_1_abab(k,l,b,c)
        //             += +0.50 P(i,j) r2_abab(c,b,i,l) t1_1_aa(a,j) l2_1_abab(k,l,c,b)
        //             += +1.00 P(i,j) r1_aa(b,i) t1_1_aa(a,j) l1_1_aa(k,b)
        //             += +0.50 P(i,j) r2_aaaa(b,c,l,i) t1_1_aa(a,j) l2_1_aaaa(l,k,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r0 t1_1_aa(a,i) l0_1
        //             += -0.50 P(i,j) r0 t2_abab(c,b,j,l) t1_1_aa(a,i) l2_1_abab(k,l,c,b)
        //             += -0.50 P(i,j) r0 t2_abab(b,c,j,l) t1_1_aa(a,i) l2_1_abab(k,l,b,c)
        //             += +0.50 P(i,j) r2_aaaa(b,c,l,i) t1_aa(a,j) l2_aaaa(l,k,b,c)
        //             += +1.00 P(i,j) r1_aa(b,i) t1_aa(a,j) l1_aa(k,b)
        //             += +0.50 P(i,j) r2_1_aaaa(b,c,l,i) t1_aa(a,j) l2_1_aaaa(l,k,b,c)
        //             += +1.00 P(i,j) r1_1_aa(b,i) t1_aa(a,j) l1_1_aa(k,b)
        //             += +0.50 P(i,j) r2_1_abab(b,c,i,l) t1_aa(a,j) l2_1_abab(k,l,b,c)
        //             += +0.50 P(i,j) r2_1_abab(c,b,i,l) t1_aa(a,j) l2_1_abab(k,l,c,b)
        //             += +0.50 P(i,j) r2_abab(b,c,i,l) t1_aa(a,j) l2_abab(k,l,b,c)
        //             += +0.50 P(i,j) r2_abab(c,b,i,l) t1_aa(a,j) l2_abab(k,l,c,b)
        //             += +1.00 P(i,j) d_aa(j,k) r0 t1_aa(a,i) l0
        //             += +0.50 P(i,j) r0_1 t1_aa(a,j) t2_abab(c,b,i,l) l2_1_abab(k,l,c,b)
        //             += +0.50 P(i,j) r0_1 t1_aa(a,j) t2_abab(b,c,i,l) l2_1_abab(k,l,b,c)
        //             += -0.50 P(i,j) r0 t2_aaaa(c,b,l,j) t1_1_aa(a,i) l2_1_aaaa(l,k,c,b)
        //             += -1.00 P(i,j) r0 t1_aa(b,j) t1_1_aa(a,i) l1_1_aa(k,b)
        //             += -1.00 P(i,j) r1_1_aa(c,i) t2_aaaa(a,b,l,j) l2_1_aaaa(l,k,b,c)
        //             += +1.00 P(i,j) r1_1_aa(c,i) t2_abab(a,b,j,l) l2_1_abab(k,l,c,b)
        //             += -1.00 P(i,j) r1_aa(c,i) t2_1_aaaa(a,b,l,j) l2_1_aaaa(l,k,b,c)
        //             += +1.00 P(i,j) r1_aa(c,i) t2_abab(a,b,j,l) l2_abab(k,l,c,b)
        //             += -1.00 P(i,j) r1_aa(c,i) t1_aa(b,j) t1_1_aa(a,l) l2_1_aaaa(l,k,b,c)
        //             += +1.00 P(i,j) r1_aa(c,i) t2_1_abab(a,b,j,l) l2_1_abab(k,l,c,b)
        //             += -1.00 P(i,j) r1_aa(c,i) t2_aaaa(a,b,l,j) l2_aaaa(l,k,b,c)
        //             += +1.00 P(i,j) r0 t1_aa(a,j) t1_1_aa(b,i) l1_1_aa(k,b)
        //             += +1.00 P(i,j) r0_1 t1_aa(a,j) t1_aa(b,i) l1_1_aa(k,b)
        //             += +1.00 P(i,j) r0 t1_aa(a,j) t1_aa(b,i) l1_aa(k,b)
        //             += +1.00 P(i,j) d_aa(j,k) r0_1 t1_aa(a,i) l0_1
        //             += +1.00 P(i,j) d_aa(j,k) r1_bb(b,l) t1_1_aa(a,i) l1_1_bb(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(b,l) t1_1_aa(a,i) l1_1_aa(l,b)
        //             += +0.25 P(i,j) d_aa(j,k) r2_aaaa(b,c,m,l) t1_1_aa(a,i) l2_1_aaaa(m,l,b,c)
        //             += +0.25 P(i,j) d_aa(j,k) r2_bbbb(b,c,m,l) t1_1_aa(a,i) l2_1_bbbb(m,l,b,c)
        //             += +0.25 P(i,j) d_aa(j,k) r2_abab(b,c,m,l) t1_1_aa(a,i) l2_1_abab(m,l,b,c)
        //             += +0.25 P(i,j) d_aa(j,k) r2_abab(b,c,l,m) t1_1_aa(a,i) l2_1_abab(l,m,b,c)
        //             += +0.25 P(i,j) d_aa(j,k) r2_abab(c,b,m,l) t1_1_aa(a,i) l2_1_abab(m,l,c,b)
        //             += +0.25 P(i,j) d_aa(j,k) r2_abab(c,b,l,m) t1_1_aa(a,i) l2_1_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r2_abab(a,c,m,l) t1_1_aa(b,i) l2_1_abab(m,l,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r2_abab(a,c,l,m) t1_1_aa(b,i) l2_1_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r2_aaaa(a,c,m,l) t1_1_aa(b,i) l2_1_aaaa(m,l,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r1_aa(c,i) t2_aaaa(a,b,l,m) l2_aaaa(l,m,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t2_abab(a,b,i,l) l2_bbbb(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t2_aaaa(a,b,l,i) l2_abab(l,m,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t2_1_abab(a,b,i,l) l2_1_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r2_aaaa(b,c,m,i) t1_aa(a,l) l2_aaaa(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(b,i) t1_aa(a,l) l1_aa(l,b)
        //             += -0.50 P(i,j) d_aa(j,k) r2_1_aaaa(b,c,m,i) t1_aa(a,l) l2_1_aaaa(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_1_aa(b,i) t1_aa(a,l) l1_1_aa(l,b)
        //             += -0.50 P(i,j) d_aa(j,k) r2_1_abab(b,c,i,m) t1_aa(a,l) l2_1_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r2_1_abab(c,b,i,m) t1_aa(a,l) l2_1_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r2_abab(b,c,i,m) t1_aa(a,l) l2_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r2_abab(c,b,i,m) t1_aa(a,l) l2_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r1_aa(a,m) t2_aaaa(c,b,l,i) l2_aaaa(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r1_aa(a,m) t2_1_aaaa(c,b,l,i) l2_1_aaaa(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r1_aa(c,i) t2_abab(a,b,l,m) l2_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r1_aa(c,i) t2_abab(a,b,m,l) l2_abab(m,l,c,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_1_bb(c,m) t2_abab(a,b,i,l) l2_1_bbbb(m,l,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r1_1_aa(a,m) t2_aaaa(c,b,l,i) l2_1_aaaa(m,l,c,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_1_bb(c,m) t2_aaaa(a,b,l,i) l2_1_abab(l,m,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_1_aa(c,m) t2_aaaa(a,b,l,i) l2_1_aaaa(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t2_1_aaaa(a,b,l,i) l2_1_abab(l,m,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t2_abab(a,b,i,l) l2_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r2_abab(b,c,i,m) t1_1_aa(a,l) l2_1_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r2_abab(c,b,i,m) t1_1_aa(a,l) l2_1_abab(l,m,c,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(b,i) t1_1_aa(a,l) l1_1_aa(l,b)
        //             += -0.50 P(i,j) d_aa(j,k) r2_aaaa(b,c,m,i) t1_1_aa(a,l) l2_1_aaaa(m,l,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r1_1_aa(c,i) t2_aaaa(a,b,l,m) l2_1_aaaa(l,m,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t2_aaaa(a,b,l,i) l2_aaaa(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t2_1_aaaa(a,b,l,i) l2_1_aaaa(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_1_aa(c,m) t2_abab(a,b,i,l) l2_1_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r1_1_aa(a,m) t2_abab(c,b,i,l) l2_1_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r1_1_aa(a,m) t2_abab(b,c,i,l) l2_1_abab(m,l,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r1_aa(c,i) t2_1_aaaa(a,b,l,m) l2_1_aaaa(l,m,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r2_1_abab(a,c,m,l) t1_aa(b,i) l2_1_abab(m,l,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r2_1_abab(a,c,l,m) t1_aa(b,i) l2_1_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r2_aaaa(a,c,m,l) t1_aa(b,i) l2_aaaa(m,l,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r2_abab(a,c,m,l) t1_aa(b,i) l2_abab(m,l,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r2_abab(a,c,l,m) t1_aa(b,i) l2_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r2_1_aaaa(a,c,m,l) t1_aa(b,i) l2_1_aaaa(m,l,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t1_aa(c,i) t2_1_aaaa(a,b,l,m) l2_1_aaaa(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t2_aaaa(a,c,l,m) t1_1_aa(b,i) l2_1_aaaa(l,m,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_abab(a,c,l,m) t1_1_aa(b,i) l2_1_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_abab(a,c,m,l) t1_1_aa(b,i) l2_1_abab(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r0 t1_aa(a,l) t1_aa(b,i) l1_aa(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r0 t2_abab(a,b,i,l) l1_bb(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r0 t1_aa(a,l) t1_1_aa(b,i) l1_1_aa(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r0 t2_aaaa(a,b,l,i) l1_aa(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r0 t1_aa(b,i) t1_1_aa(a,l) l1_1_aa(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r0 t2_1_abab(a,b,i,l) l1_1_bb(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r0 t2_1_aaaa(a,b,l,i) l1_1_aa(l,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t2_aaaa(c,b,m,i) t1_1_aa(a,l) l2_1_aaaa(l,m,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_abab(a,b,l,m) t1_aa(c,i) l2_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_abab(a,b,m,l) t1_aa(c,i) l2_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_aaaa(a,b,l,m) t1_aa(c,i) l2_aaaa(l,m,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t1_aa(c,i) t2_1_abab(a,b,l,m) l2_1_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t1_aa(c,i) t2_1_abab(a,b,m,l) l2_1_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_abab(c,b,i,m) t1_1_aa(a,l) l2_1_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_abab(b,c,i,m) t1_1_aa(a,l) l2_1_abab(l,m,b,c)
        //             += -1.00 P(i,j) r1_aa(c,l) t1_aa(a,j) t1_1_aa(b,i) l2_1_aaaa(l,k,b,c)
        //             += +1.00 P(i,j) r1_bb(c,l) t1_aa(a,j) t1_1_aa(b,i) l2_1_abab(k,l,b,c)
        //             += -1.00 P(i,j) r0 t1_aa(c,j) t2_1_aaaa(a,b,l,i) l2_1_aaaa(l,k,c,b)
        //             += -1.00 P(i,j) r0 t2_aaaa(a,b,l,i) t1_aa(c,j) l2_aaaa(l,k,c,b)
        //             += -1.00 P(i,j) r0 t2_aaaa(a,c,l,j) t1_1_aa(b,i) l2_1_aaaa(l,k,c,b)
        //             += -1.00 P(i,j) r0_1 t2_aaaa(a,b,l,i) t1_aa(c,j) l2_1_aaaa(l,k,c,b)
        //             += -0.50 P(i,j) r1_1_aa(a,i) t2_aaaa(c,b,l,j) l2_1_aaaa(l,k,c,b)
        //             += -1.00 P(i,j) r0 t1_aa(a,l) t1_aa(c,j) t1_1_aa(b,i) l2_1_aaaa(l,k,c,b)
        //             += -1.00 P(i,j) r1_1_aa(a,i) t1_aa(b,j) l1_1_aa(k,b)
        //             += -0.50 P(i,j) d_aa(j,k) r1_aa(c,i) t2_1_abab(a,b,l,m) l2_1_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r1_aa(c,i) t2_1_abab(a,b,m,l) l2_1_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r1_1_aa(c,i) t2_abab(a,b,l,m) l2_1_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r1_1_aa(c,i) t2_abab(a,b,m,l) l2_1_abab(m,l,c,b)
        D_ooov_aaaa("L,R,i,j,k,a") -= tmps_["309_aaaa_LLoovo"]("L,R,j,k,a,i");
        D_ooov_aaaa("L,R,i,j,k,a") += tmps_["309_aaaa_LLoovo"]("L,R,i,k,a,j");

        // D_oovo_aaaa += +1.00 P(i,j) r1_aa(a,i) t1_aa(b,j) l1_aa(k,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(a,i) l0
        //             += +1.00 P(i,j) r1_aa(a,i) t1_1_aa(b,j) l1_1_aa(k,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_1_aa(a,i) l0_1
        //             += +0.50 P(i,j) r1_1_aa(a,i) t2_abab(c,b,j,l) l2_1_abab(k,l,c,b)
        //             += +0.50 P(i,j) r1_1_aa(a,i) t2_abab(b,c,j,l) l2_1_abab(k,l,b,c)
        //             += +1.00 P(i,j) r0_1 t2_abab(a,b,i,l) t1_aa(c,j) l2_1_abab(k,l,c,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t1_aa(b,i) t1_1_aa(a,l) l2_1_abab(l,m,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t1_aa(b,i) t1_1_aa(a,l) l2_1_aaaa(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t1_aa(a,l) t1_aa(b,i) l2_abab(l,m,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_1_bb(c,m) t1_aa(a,l) t1_aa(b,i) l2_1_abab(l,m,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t1_aa(a,l) t1_aa(b,i) l2_aaaa(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_1_aa(c,m) t1_aa(a,l) t1_aa(b,i) l2_1_aaaa(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t1_aa(a,l) t1_1_aa(b,i) l2_1_aaaa(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t1_aa(a,l) t1_1_aa(b,i) l2_1_abab(l,m,b,c)
        //             += +1.00 P(i,j) r0 t1_aa(c,j) t2_1_abab(a,b,i,l) l2_1_abab(k,l,c,b)
        //             += +1.00 P(i,j) r0 t2_abab(a,b,i,l) t1_aa(c,j) l2_abab(k,l,c,b)
        //             += -1.00 P(i,j) r0 t2_abab(a,c,j,l) t1_1_aa(b,i) l2_1_abab(k,l,b,c)
        //             += -1.00 P(i,j) r1_bb(c,l) t1_aa(a,j) t1_aa(b,i) l2_abab(k,l,b,c)
        //             += -1.00 P(i,j) r1_1_bb(c,l) t1_aa(a,j) t1_aa(b,i) l2_1_abab(k,l,b,c)
        //             += +1.00 P(i,j) r1_aa(c,l) t1_aa(a,j) t1_aa(b,i) l2_aaaa(l,k,b,c)
        //             += +1.00 P(i,j) r1_1_aa(c,l) t1_aa(a,j) t1_aa(b,i) l2_1_aaaa(l,k,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r0_1 t2_abab(a,b,l,m) t1_aa(c,i) l2_1_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0_1 t2_abab(a,b,m,l) t1_aa(c,i) l2_1_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0_1 t1_aa(a,m) t2_aaaa(c,b,l,i) l2_1_aaaa(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0_1 t1_aa(a,m) t2_abab(c,b,i,l) l2_1_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0_1 t1_aa(a,m) t2_abab(b,c,i,l) l2_1_abab(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r0_1 t1_aa(a,l) t1_aa(b,i) l1_1_aa(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r0_1 t2_aaaa(a,b,l,i) l1_1_aa(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r0_1 t2_abab(a,b,i,l) l1_1_bb(l,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0_1 t2_aaaa(a,b,l,m) t1_aa(c,i) l2_1_aaaa(l,m,c,b)
        //             += +1.00 P(i,j) r1_bb(c,l) t1_aa(b,j) t1_1_aa(a,i) l2_1_abab(k,l,b,c)
        //             += -1.00 P(i,j) r1_aa(c,l) t1_aa(b,j) t1_1_aa(a,i) l2_1_aaaa(l,k,b,c)
        //             += -0.25 P(i,j) d_aa(j,k) r2_1_abab(b,c,m,l) t1_aa(a,i) l2_1_abab(m,l,b,c)
        //             += -0.25 P(i,j) d_aa(j,k) r2_1_abab(b,c,l,m) t1_aa(a,i) l2_1_abab(l,m,b,c)
        //             += -0.25 P(i,j) d_aa(j,k) r2_1_abab(c,b,m,l) t1_aa(a,i) l2_1_abab(m,l,c,b)
        //             += -0.25 P(i,j) d_aa(j,k) r2_1_abab(c,b,l,m) t1_aa(a,i) l2_1_abab(l,m,c,b)
        //             += -0.25 P(i,j) d_aa(j,k) r2_aaaa(b,c,m,l) t1_aa(a,i) l2_aaaa(m,l,b,c)
        //             += -0.25 P(i,j) d_aa(j,k) r2_bbbb(b,c,m,l) t1_aa(a,i) l2_bbbb(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_1_aa(b,l) t1_aa(a,i) l1_1_aa(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_1_bb(b,l) t1_aa(a,i) l1_1_bb(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(b,l) t1_aa(a,i) l1_aa(l,b)
        //             += -0.25 P(i,j) d_aa(j,k) r2_1_bbbb(b,c,m,l) t1_aa(a,i) l2_1_bbbb(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_bb(b,l) t1_aa(a,i) l1_bb(l,b)
        //             += -0.25 P(i,j) d_aa(j,k) r2_1_aaaa(b,c,m,l) t1_aa(a,i) l2_1_aaaa(m,l,b,c)
        //             += -0.25 P(i,j) d_aa(j,k) r2_abab(b,c,m,l) t1_aa(a,i) l2_abab(m,l,b,c)
        //             += -0.25 P(i,j) d_aa(j,k) r2_abab(b,c,l,m) t1_aa(a,i) l2_abab(l,m,b,c)
        //             += -0.25 P(i,j) d_aa(j,k) r2_abab(c,b,m,l) t1_aa(a,i) l2_abab(m,l,c,b)
        //             += -0.25 P(i,j) d_aa(j,k) r2_abab(c,b,l,m) t1_aa(a,i) l2_abab(l,m,c,b)
        //             += +1.00 P(i,j) d_aa(j,k) r2_1_aaaa(a,b,l,i) l1_1_aa(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(a,l) t1_1_aa(b,i) l1_1_aa(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r2_aaaa(a,b,l,i) l1_aa(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t2_1_abab(a,b,i,l) l2_1_bbbb(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(a,l) t1_aa(b,i) l1_aa(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r2_1_abab(a,b,i,l) l1_1_bb(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r2_abab(a,b,i,l) l1_bb(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_1_aa(a,l) t1_aa(b,i) l1_1_aa(l,b)
        //             += -0.50 P(i,j) r0_1 t1_aa(a,j) t2_aaaa(c,b,l,i) l2_1_aaaa(l,k,c,b)
        //             += -0.50 P(i,j) r2_abab(b,c,i,l) t1_1_aa(a,j) l2_1_abab(k,l,b,c)
        //             += -0.50 P(i,j) r2_abab(c,b,i,l) t1_1_aa(a,j) l2_1_abab(k,l,c,b)
        //             += -1.00 P(i,j) r1_aa(b,i) t1_1_aa(a,j) l1_1_aa(k,b)
        //             += -0.50 P(i,j) r2_aaaa(b,c,l,i) t1_1_aa(a,j) l2_1_aaaa(l,k,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r0 t1_1_aa(a,i) l0_1
        //             += +0.50 P(i,j) r0 t2_abab(c,b,j,l) t1_1_aa(a,i) l2_1_abab(k,l,c,b)
        //             += +0.50 P(i,j) r0 t2_abab(b,c,j,l) t1_1_aa(a,i) l2_1_abab(k,l,b,c)
        //             += -0.50 P(i,j) r2_aaaa(b,c,l,i) t1_aa(a,j) l2_aaaa(l,k,b,c)
        //             += -1.00 P(i,j) r1_aa(b,i) t1_aa(a,j) l1_aa(k,b)
        //             += -0.50 P(i,j) r2_1_aaaa(b,c,l,i) t1_aa(a,j) l2_1_aaaa(l,k,b,c)
        //             += -1.00 P(i,j) r1_1_aa(b,i) t1_aa(a,j) l1_1_aa(k,b)
        //             += -0.50 P(i,j) r2_1_abab(b,c,i,l) t1_aa(a,j) l2_1_abab(k,l,b,c)
        //             += -0.50 P(i,j) r2_1_abab(c,b,i,l) t1_aa(a,j) l2_1_abab(k,l,c,b)
        //             += -0.50 P(i,j) r2_abab(b,c,i,l) t1_aa(a,j) l2_abab(k,l,b,c)
        //             += -0.50 P(i,j) r2_abab(c,b,i,l) t1_aa(a,j) l2_abab(k,l,c,b)
        //             += -1.00 P(i,j) d_aa(j,k) r0 t1_aa(a,i) l0
        //             += -0.50 P(i,j) r0_1 t1_aa(a,j) t2_abab(c,b,i,l) l2_1_abab(k,l,c,b)
        //             += -0.50 P(i,j) r0_1 t1_aa(a,j) t2_abab(b,c,i,l) l2_1_abab(k,l,b,c)
        //             += +0.50 P(i,j) r0 t2_aaaa(c,b,l,j) t1_1_aa(a,i) l2_1_aaaa(l,k,c,b)
        //             += +1.00 P(i,j) r0 t1_aa(b,j) t1_1_aa(a,i) l1_1_aa(k,b)
        //             += +1.00 P(i,j) r1_1_aa(c,i) t2_aaaa(a,b,l,j) l2_1_aaaa(l,k,b,c)
        //             += -1.00 P(i,j) r1_1_aa(c,i) t2_abab(a,b,j,l) l2_1_abab(k,l,c,b)
        //             += +1.00 P(i,j) r1_aa(c,i) t2_1_aaaa(a,b,l,j) l2_1_aaaa(l,k,b,c)
        //             += -1.00 P(i,j) r1_aa(c,i) t2_abab(a,b,j,l) l2_abab(k,l,c,b)
        //             += +1.00 P(i,j) r1_aa(c,i) t1_aa(b,j) t1_1_aa(a,l) l2_1_aaaa(l,k,b,c)
        //             += -1.00 P(i,j) r1_aa(c,i) t2_1_abab(a,b,j,l) l2_1_abab(k,l,c,b)
        //             += +1.00 P(i,j) r1_aa(c,i) t2_aaaa(a,b,l,j) l2_aaaa(l,k,b,c)
        //             += -1.00 P(i,j) r0 t1_aa(a,j) t1_1_aa(b,i) l1_1_aa(k,b)
        //             += -1.00 P(i,j) r0_1 t1_aa(a,j) t1_aa(b,i) l1_1_aa(k,b)
        //             += -1.00 P(i,j) r0 t1_aa(a,j) t1_aa(b,i) l1_aa(k,b)
        //             += -1.00 P(i,j) d_aa(j,k) r0_1 t1_aa(a,i) l0_1
        //             += -1.00 P(i,j) d_aa(j,k) r1_bb(b,l) t1_1_aa(a,i) l1_1_bb(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(b,l) t1_1_aa(a,i) l1_1_aa(l,b)
        //             += -0.25 P(i,j) d_aa(j,k) r2_aaaa(b,c,m,l) t1_1_aa(a,i) l2_1_aaaa(m,l,b,c)
        //             += -0.25 P(i,j) d_aa(j,k) r2_bbbb(b,c,m,l) t1_1_aa(a,i) l2_1_bbbb(m,l,b,c)
        //             += -0.25 P(i,j) d_aa(j,k) r2_abab(b,c,m,l) t1_1_aa(a,i) l2_1_abab(m,l,b,c)
        //             += -0.25 P(i,j) d_aa(j,k) r2_abab(b,c,l,m) t1_1_aa(a,i) l2_1_abab(l,m,b,c)
        //             += -0.25 P(i,j) d_aa(j,k) r2_abab(c,b,m,l) t1_1_aa(a,i) l2_1_abab(m,l,c,b)
        //             += -0.25 P(i,j) d_aa(j,k) r2_abab(c,b,l,m) t1_1_aa(a,i) l2_1_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r2_abab(a,c,m,l) t1_1_aa(b,i) l2_1_abab(m,l,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r2_abab(a,c,l,m) t1_1_aa(b,i) l2_1_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r2_aaaa(a,c,m,l) t1_1_aa(b,i) l2_1_aaaa(m,l,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r1_aa(c,i) t2_aaaa(a,b,l,m) l2_aaaa(l,m,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t2_abab(a,b,i,l) l2_bbbb(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t2_aaaa(a,b,l,i) l2_abab(l,m,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t2_1_abab(a,b,i,l) l2_1_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r2_aaaa(b,c,m,i) t1_aa(a,l) l2_aaaa(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(b,i) t1_aa(a,l) l1_aa(l,b)
        //             += +0.50 P(i,j) d_aa(j,k) r2_1_aaaa(b,c,m,i) t1_aa(a,l) l2_1_aaaa(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_1_aa(b,i) t1_aa(a,l) l1_1_aa(l,b)
        //             += +0.50 P(i,j) d_aa(j,k) r2_1_abab(b,c,i,m) t1_aa(a,l) l2_1_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r2_1_abab(c,b,i,m) t1_aa(a,l) l2_1_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r2_abab(b,c,i,m) t1_aa(a,l) l2_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r2_abab(c,b,i,m) t1_aa(a,l) l2_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r1_aa(a,m) t2_aaaa(c,b,l,i) l2_aaaa(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r1_aa(a,m) t2_1_aaaa(c,b,l,i) l2_1_aaaa(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r1_aa(c,i) t2_abab(a,b,l,m) l2_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r1_aa(c,i) t2_abab(a,b,m,l) l2_abab(m,l,c,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_1_bb(c,m) t2_abab(a,b,i,l) l2_1_bbbb(m,l,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r1_1_aa(a,m) t2_aaaa(c,b,l,i) l2_1_aaaa(m,l,c,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_1_bb(c,m) t2_aaaa(a,b,l,i) l2_1_abab(l,m,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_1_aa(c,m) t2_aaaa(a,b,l,i) l2_1_aaaa(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t2_1_aaaa(a,b,l,i) l2_1_abab(l,m,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t2_abab(a,b,i,l) l2_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r2_abab(b,c,i,m) t1_1_aa(a,l) l2_1_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r2_abab(c,b,i,m) t1_1_aa(a,l) l2_1_abab(l,m,c,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(b,i) t1_1_aa(a,l) l1_1_aa(l,b)
        //             += +0.50 P(i,j) d_aa(j,k) r2_aaaa(b,c,m,i) t1_1_aa(a,l) l2_1_aaaa(m,l,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r1_1_aa(c,i) t2_aaaa(a,b,l,m) l2_1_aaaa(l,m,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t2_aaaa(a,b,l,i) l2_aaaa(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t2_1_aaaa(a,b,l,i) l2_1_aaaa(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_1_aa(c,m) t2_abab(a,b,i,l) l2_1_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r1_1_aa(a,m) t2_abab(c,b,i,l) l2_1_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r1_1_aa(a,m) t2_abab(b,c,i,l) l2_1_abab(m,l,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r1_aa(c,i) t2_1_aaaa(a,b,l,m) l2_1_aaaa(l,m,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r2_1_abab(a,c,m,l) t1_aa(b,i) l2_1_abab(m,l,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r2_1_abab(a,c,l,m) t1_aa(b,i) l2_1_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r2_aaaa(a,c,m,l) t1_aa(b,i) l2_aaaa(m,l,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r2_abab(a,c,m,l) t1_aa(b,i) l2_abab(m,l,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r2_abab(a,c,l,m) t1_aa(b,i) l2_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r2_1_aaaa(a,c,m,l) t1_aa(b,i) l2_1_aaaa(m,l,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t1_aa(c,i) t2_1_aaaa(a,b,l,m) l2_1_aaaa(l,m,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_aaaa(a,c,l,m) t1_1_aa(b,i) l2_1_aaaa(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t2_abab(a,c,l,m) t1_1_aa(b,i) l2_1_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t2_abab(a,c,m,l) t1_1_aa(b,i) l2_1_abab(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r0 t1_aa(a,l) t1_aa(b,i) l1_aa(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r0 t2_abab(a,b,i,l) l1_bb(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r0 t1_aa(a,l) t1_1_aa(b,i) l1_1_aa(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r0 t2_aaaa(a,b,l,i) l1_aa(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r0 t1_aa(b,i) t1_1_aa(a,l) l1_1_aa(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r0 t2_1_abab(a,b,i,l) l1_1_bb(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r0 t2_1_aaaa(a,b,l,i) l1_1_aa(l,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_aaaa(c,b,m,i) t1_1_aa(a,l) l2_1_aaaa(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t2_abab(a,b,l,m) t1_aa(c,i) l2_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t2_abab(a,b,m,l) t1_aa(c,i) l2_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t2_aaaa(a,b,l,m) t1_aa(c,i) l2_aaaa(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t1_aa(c,i) t2_1_abab(a,b,l,m) l2_1_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t1_aa(c,i) t2_1_abab(a,b,m,l) l2_1_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t2_abab(c,b,i,m) t1_1_aa(a,l) l2_1_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t2_abab(b,c,i,m) t1_1_aa(a,l) l2_1_abab(l,m,b,c)
        //             += +1.00 P(i,j) r1_aa(c,l) t1_aa(a,j) t1_1_aa(b,i) l2_1_aaaa(l,k,b,c)
        //             += -1.00 P(i,j) r1_bb(c,l) t1_aa(a,j) t1_1_aa(b,i) l2_1_abab(k,l,b,c)
        //             += +1.00 P(i,j) r0 t1_aa(c,j) t2_1_aaaa(a,b,l,i) l2_1_aaaa(l,k,c,b)
        //             += +1.00 P(i,j) r0 t2_aaaa(a,b,l,i) t1_aa(c,j) l2_aaaa(l,k,c,b)
        //             += +1.00 P(i,j) r0 t2_aaaa(a,c,l,j) t1_1_aa(b,i) l2_1_aaaa(l,k,c,b)
        //             += +1.00 P(i,j) r0_1 t2_aaaa(a,b,l,i) t1_aa(c,j) l2_1_aaaa(l,k,c,b)
        //             += +0.50 P(i,j) r1_1_aa(a,i) t2_aaaa(c,b,l,j) l2_1_aaaa(l,k,c,b)
        //             += +1.00 P(i,j) r0 t1_aa(a,l) t1_aa(c,j) t1_1_aa(b,i) l2_1_aaaa(l,k,c,b)
        //             += +1.00 P(i,j) r1_1_aa(a,i) t1_aa(b,j) l1_1_aa(k,b)
        //             += +0.50 P(i,j) d_aa(j,k) r1_aa(c,i) t2_1_abab(a,b,l,m) l2_1_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r1_aa(c,i) t2_1_abab(a,b,m,l) l2_1_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r1_1_aa(c,i) t2_abab(a,b,l,m) l2_1_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r1_1_aa(c,i) t2_abab(a,b,m,l) l2_1_abab(m,l,c,b)
        D_oovo_aaaa("L,R,i,j,a,k") += tmps_["309_aaaa_LLoovo"]("L,R,j,k,a,i");
        D_oovo_aaaa("L,R,i,j,a,k") -= tmps_["309_aaaa_LLoovo"]("L,R,i,k,a,j");
        tmps_["309_aaaa_LLoovo"].~TArrayD();

        // flops: o2v2L2  = o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o3v2L1 o1v1L1 o1v1L1 o3v2L2 o3v2L1 o3v2L2 o2v2L2 o3v2L2 o2v2L2 o2v2L2 o2v2L2 o3v2L2 o3v2L2 o2v2L2 o3v2L1 o3v2L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o3v2L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o3v2L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2
        //  mems: o2v2L2  = o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o3v1L1 o1v1L1 o1v1L1 o2v2L2 o3v1L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o3v1L2 o2v2L2 o2v2L2 o3v1L1 o2v2L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j")  = -1.00 * tmps_["85_bb_Lvv"]("L,a,b") * tmps_["112_bb_Loo"]("R,i,j");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") -= 0.50 * tmps_["219_bb_Lvv"]("L,b,a") * tmps_["113_bb_Loo"]("R,i,j");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") -= tmps_["276_bb_LLvv"]("L,R,a,b") * Id["bb_oo"]("i,j");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += t1_1["bb"]("b,i") * tmps_["69_bb_LLov"]("L,R,j,a");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += t1["bb"]("b,i") * tmps_["70_bb_LLov"]("L,R,j,a");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += l2["bbbb"]("L,k,j,a,c") * t1["bb"]("c,i") * (r0("R") * t1["bb"]("b,k") + r1["bb"]("R,b,k"));
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += t1_1["bb"]("c,i") * l2_1["bbbb"]("L,k,j,a,c") * r1["bb"]("R,b,k");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += tmps_["171_bbbb_Looov"]("L,i,k,j,a") * r1_1["bb"]("R,b,k");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += l1["bb"]("L,j,a") * r1["bb"]("R,b,i");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += l2_1["bbbb"]("L,k,j,a,c") * r1["bb"]("R,c,i") * t1_1["bb"]("b,k");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += t1_1["bb"]("c,i") * l2_1["bbbb"]("L,k,j,a,c") * t1["bb"]("b,k") * r0("R");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += l1_1["bb"]("L,j,a") * r1_1["bb"]("R,b,i");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += t1["bb"]("b,i") * tmps_["71_bb_LLov"]("L,R,j,a");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += l1_1["bb"]("L,j,a") * tmps_["225_bb_Lvo"]("R,b,i");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") -= t1["bb"]("b,i") * tmps_["73_bb_LLov"]("L,R,j,a");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") -= t1["bb"]("b,i") * tmps_["74_bb_LLov"]("L,R,j,a");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += t1_1["bb"]("b,k") * tmps_["171_bbbb_Looov"]("L,i,k,j,a") * r0("R");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += l1_1["bb"]("L,j,a") * tmps_["226_bb_Lvo"]("R,b,i");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") -= t1_1["bb"]("b,i") * tmps_["72_bb_LLov"]("L,R,j,a");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += l1["bb"]("L,j,a") * tmps_["224_bb_Lvo"]("R,b,i");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") -= 0.50 * tmps_["210_bb_LLvv"]("L,R,a,b") * Id["bb_oo"]("i,j");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += t1["bb"]("b,k") * tmps_["171_bbbb_Looov"]("L,i,k,j,a") * r0_1("R");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") -= 0.50 * tmps_["222_bb_Lvv"]("L,b,a") * tmps_["112_bb_Loo"]("R,i,j");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") -= tmps_["217_bb_Lvv"]("L,a,b") * tmps_["112_bb_Loo"]("R,i,j");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") -= tmps_["84_bb_Lvv"]("L,b,a") * tmps_["113_bb_Loo"]("R,i,j");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") -= tmps_["277_bb_Lvv"]("L,a,b") * tmps_["112_bb_Loo"]("R,i,j");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") -= 0.50 * tmps_["221_bb_Lvv"]("L,a,b") * tmps_["112_bb_Loo"]("R,i,j");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") -= tmps_["275_bb_Lvv"]("L,a,b") * tmps_["113_bb_Loo"]("R,i,j");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += r0("R") * tmps_["208_bbbb_Lvoov"]("L,b,i,j,a");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += r0("R") * tmps_["213_bbbb_Lvoov"]("L,b,i,j,a");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += r0("R") * tmps_["82_bbbb_Lvoov"]("L,b,i,j,a");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += r0_1("R") * tmps_["205_bbbb_Lvoov"]("L,b,i,j,a");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += r0_1("R") * tmps_["81_bbbb_Lvoov"]("L,b,i,j,a");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += r0("R") * tmps_["214_bbbb_Lvoov"]("L,b,i,j,a");
        tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j") += Id["bb_oo"]("i,j") * tmps_["278_bb_LLvv"]("L,R,a,b");
        tmps_["278_bb_LLvv"].~TArrayD();
        tmps_["277_bb_Lvv"].~TArrayD();
        tmps_["276_bb_LLvv"].~TArrayD();
        tmps_["275_bb_Lvv"].~TArrayD();
        tmps_["226_bb_Lvo"].~TArrayD();
        tmps_["225_bb_Lvo"].~TArrayD();
        tmps_["224_bb_Lvo"].~TArrayD();
        tmps_["222_bb_Lvv"].~TArrayD();
        tmps_["221_bb_Lvv"].~TArrayD();
        tmps_["219_bb_Lvv"].~TArrayD();
        tmps_["217_bb_Lvv"].~TArrayD();
        tmps_["214_bbbb_Lvoov"].~TArrayD();
        tmps_["213_bbbb_Lvoov"].~TArrayD();
        tmps_["210_bb_LLvv"].~TArrayD();
        tmps_["208_bbbb_Lvoov"].~TArrayD();
        tmps_["205_bbbb_Lvoov"].~TArrayD();
        tmps_["171_bbbb_Looov"].~TArrayD();
        tmps_["85_bb_Lvv"].~TArrayD();
        tmps_["84_bb_Lvv"].~TArrayD();
        tmps_["82_bbbb_Lvoov"].~TArrayD();
        tmps_["81_bbbb_Lvoov"].~TArrayD();
        tmps_["74_bb_LLov"].~TArrayD();
        tmps_["73_bb_LLov"].~TArrayD();
        tmps_["72_bb_LLov"].~TArrayD();
        tmps_["71_bb_LLov"].~TArrayD();
        tmps_["70_bb_LLov"].~TArrayD();
        tmps_["69_bb_LLov"].~TArrayD();

        // D_ovov_bbbb += +1.00 d_bb(i,j) r1_bb(c,l) t1_1_bb(b,k) l2_1_bbbb(l,k,a,c)
        //             += -1.00 d_bb(i,j) r1_aa(c,l) t1_1_bb(b,k) l2_1_abab(l,k,c,a)
        //             += +1.00 d_bb(i,j) r1_1_bb(c,l) t1_bb(b,k) l2_1_bbbb(l,k,a,c)
        //             += -1.00 d_bb(i,j) r1_1_aa(c,l) t1_bb(b,k) l2_1_abab(l,k,c,a)
        //             += +1.00 d_bb(i,j) r1_bb(c,l) t1_bb(b,k) l2_bbbb(l,k,a,c)
        //             += -1.00 d_bb(i,j) r1_aa(c,l) t1_bb(b,k) l2_abab(l,k,c,a)
        //             += +1.00 r0 t2_1_abab(c,b,k,i) l2_1_abab(k,j,c,a)
        //             += -0.50 d_bb(i,j) r0 t2_abab(c,b,k,l) l2_abab(k,l,c,a)
        //             += -0.50 d_bb(i,j) r0 t2_abab(c,b,l,k) l2_abab(l,k,c,a)
        //             += -0.50 d_bb(i,j) r0_1 t2_bbbb(b,c,k,l) l2_1_bbbb(k,l,a,c)
        //             += -1.00 d_bb(i,j) r1_1_bb(b,k) l1_1_bb(k,a)
        //             += -1.00 d_bb(i,j) r1_bb(b,k) l1_bb(k,a)
        //             += +1.00 r1_aa(c,k) t1_1_bb(b,i) l2_1_abab(k,j,c,a)
        //             += +1.00 r1_aa(c,k) t1_bb(b,i) l2_abab(k,j,c,a)
        //             += +1.00 r0 t1_bb(b,k) t1_bb(c,i) l2_bbbb(k,j,a,c)
        //             += +1.00 r1_bb(b,k) t1_bb(c,i) l2_bbbb(k,j,a,c)
        //             += +1.00 r1_bb(b,k) t1_1_bb(c,i) l2_1_bbbb(k,j,a,c)
        //             += +1.00 r1_1_bb(b,k) t1_bb(c,i) l2_1_bbbb(k,j,a,c)
        //             += +1.00 r1_bb(b,i) l1_bb(j,a)
        //             += +1.00 r1_bb(c,i) t1_1_bb(b,k) l2_1_bbbb(k,j,a,c)
        //             += +1.00 r0 t1_bb(b,k) t1_1_bb(c,i) l2_1_bbbb(k,j,a,c)
        //             += +1.00 r1_1_bb(b,i) l1_1_bb(j,a)
        //             += +1.00 r1_1_aa(c,k) t1_bb(b,i) l2_1_abab(k,j,c,a)
        //             += +1.00 r0_1 t1_bb(b,i) l1_1_bb(j,a)
        //             += -1.00 r1_1_bb(c,k) t1_bb(b,i) l2_1_bbbb(k,j,a,c)
        //             += -1.00 r1_bb(c,k) t1_bb(b,i) l2_bbbb(k,j,a,c)
        //             += +1.00 r0 t1_bb(c,i) t1_1_bb(b,k) l2_1_bbbb(k,j,a,c)
        //             += +1.00 r0 t1_1_bb(b,i) l1_1_bb(j,a)
        //             += -1.00 r1_bb(c,k) t1_1_bb(b,i) l2_1_bbbb(k,j,a,c)
        //             += +1.00 r0 t1_bb(b,i) l1_bb(j,a)
        //             += -0.50 d_bb(i,j) r2_bbbb(b,c,l,k) l2_bbbb(l,k,a,c)
        //             += -0.50 d_bb(i,j) r2_1_abab(c,b,l,k) l2_1_abab(l,k,c,a)
        //             += -0.50 d_bb(i,j) r2_1_abab(c,b,k,l) l2_1_abab(k,l,c,a)
        //             += -0.50 d_bb(i,j) r2_1_bbbb(b,c,l,k) l2_1_bbbb(l,k,a,c)
        //             += -0.50 d_bb(i,j) r2_abab(c,b,l,k) l2_abab(l,k,c,a)
        //             += -0.50 d_bb(i,j) r2_abab(c,b,k,l) l2_abab(k,l,c,a)
        //             += +1.00 r0_1 t1_bb(b,k) t1_bb(c,i) l2_1_bbbb(k,j,a,c)
        //             += -0.50 d_bb(i,j) r0 t2_bbbb(b,c,k,l) l2_bbbb(k,l,a,c)
        //             += -0.50 d_bb(i,j) r0 t2_1_abab(c,b,k,l) l2_1_abab(k,l,c,a)
        //             += -0.50 d_bb(i,j) r0 t2_1_abab(c,b,l,k) l2_1_abab(l,k,c,a)
        //             += -0.50 d_bb(i,j) r0_1 t2_abab(c,b,k,l) l2_1_abab(k,l,c,a)
        //             += -0.50 d_bb(i,j) r0_1 t2_abab(c,b,l,k) l2_1_abab(l,k,c,a)
        //             += -1.00 d_bb(i,j) r0 t1_1_bb(b,k) l1_1_bb(k,a)
        //             += -1.00 d_bb(i,j) r0 t1_bb(b,k) l1_bb(k,a)
        //             += -0.50 d_bb(i,j) r0 t2_1_bbbb(b,c,k,l) l2_1_bbbb(k,l,a,c)
        //             += -1.00 d_bb(i,j) r0_1 t1_bb(b,k) l1_1_bb(k,a)
        //             += +1.00 r0 t2_bbbb(b,c,k,i) l2_bbbb(k,j,a,c)
        //             += +1.00 r0 t2_abab(c,b,k,i) l2_abab(k,j,c,a)
        //             += +1.00 r0_1 t2_bbbb(b,c,k,i) l2_1_bbbb(k,j,a,c)
        //             += +1.00 r0_1 t2_abab(c,b,k,i) l2_1_abab(k,j,c,a)
        //             += +1.00 r0 t2_1_bbbb(b,c,k,i) l2_1_bbbb(k,j,a,c)
        D_ovov_bbbb("L,R,i,a,j,b") += tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j");

        // D_ovvo_bbbb += -1.00 d_bb(i,j) r1_bb(c,l) t1_1_bb(b,k) l2_1_bbbb(l,k,a,c)
        //             += +1.00 d_bb(i,j) r1_aa(c,l) t1_1_bb(b,k) l2_1_abab(l,k,c,a)
        //             += -1.00 d_bb(i,j) r1_1_bb(c,l) t1_bb(b,k) l2_1_bbbb(l,k,a,c)
        //             += +1.00 d_bb(i,j) r1_1_aa(c,l) t1_bb(b,k) l2_1_abab(l,k,c,a)
        //             += -1.00 d_bb(i,j) r1_bb(c,l) t1_bb(b,k) l2_bbbb(l,k,a,c)
        //             += +1.00 d_bb(i,j) r1_aa(c,l) t1_bb(b,k) l2_abab(l,k,c,a)
        //             += -1.00 r0 t2_1_abab(c,b,k,i) l2_1_abab(k,j,c,a)
        //             += +0.50 d_bb(i,j) r0 t2_abab(c,b,k,l) l2_abab(k,l,c,a)
        //             += +0.50 d_bb(i,j) r0 t2_abab(c,b,l,k) l2_abab(l,k,c,a)
        //             += +0.50 d_bb(i,j) r0_1 t2_bbbb(b,c,k,l) l2_1_bbbb(k,l,a,c)
        //             += +1.00 d_bb(i,j) r1_1_bb(b,k) l1_1_bb(k,a)
        //             += +1.00 d_bb(i,j) r1_bb(b,k) l1_bb(k,a)
        //             += -1.00 r1_aa(c,k) t1_1_bb(b,i) l2_1_abab(k,j,c,a)
        //             += -1.00 r1_aa(c,k) t1_bb(b,i) l2_abab(k,j,c,a)
        //             += -1.00 r0 t1_bb(b,k) t1_bb(c,i) l2_bbbb(k,j,a,c)
        //             += -1.00 r1_bb(b,k) t1_bb(c,i) l2_bbbb(k,j,a,c)
        //             += -1.00 r1_bb(b,k) t1_1_bb(c,i) l2_1_bbbb(k,j,a,c)
        //             += -1.00 r1_1_bb(b,k) t1_bb(c,i) l2_1_bbbb(k,j,a,c)
        //             += -1.00 r1_bb(b,i) l1_bb(j,a)
        //             += -1.00 r1_bb(c,i) t1_1_bb(b,k) l2_1_bbbb(k,j,a,c)
        //             += -1.00 r0 t1_bb(b,k) t1_1_bb(c,i) l2_1_bbbb(k,j,a,c)
        //             += -1.00 r1_1_bb(b,i) l1_1_bb(j,a)
        //             += -1.00 r1_1_aa(c,k) t1_bb(b,i) l2_1_abab(k,j,c,a)
        //             += -1.00 r0_1 t1_bb(b,i) l1_1_bb(j,a)
        //             += +1.00 r1_1_bb(c,k) t1_bb(b,i) l2_1_bbbb(k,j,a,c)
        //             += +1.00 r1_bb(c,k) t1_bb(b,i) l2_bbbb(k,j,a,c)
        //             += -1.00 r0 t1_bb(c,i) t1_1_bb(b,k) l2_1_bbbb(k,j,a,c)
        //             += -1.00 r0 t1_1_bb(b,i) l1_1_bb(j,a)
        //             += +1.00 r1_bb(c,k) t1_1_bb(b,i) l2_1_bbbb(k,j,a,c)
        //             += -1.00 r0 t1_bb(b,i) l1_bb(j,a)
        //             += +0.50 d_bb(i,j) r2_bbbb(b,c,l,k) l2_bbbb(l,k,a,c)
        //             += +0.50 d_bb(i,j) r2_1_abab(c,b,l,k) l2_1_abab(l,k,c,a)
        //             += +0.50 d_bb(i,j) r2_1_abab(c,b,k,l) l2_1_abab(k,l,c,a)
        //             += +0.50 d_bb(i,j) r2_1_bbbb(b,c,l,k) l2_1_bbbb(l,k,a,c)
        //             += +0.50 d_bb(i,j) r2_abab(c,b,l,k) l2_abab(l,k,c,a)
        //             += +0.50 d_bb(i,j) r2_abab(c,b,k,l) l2_abab(k,l,c,a)
        //             += -1.00 r0_1 t1_bb(b,k) t1_bb(c,i) l2_1_bbbb(k,j,a,c)
        //             += +0.50 d_bb(i,j) r0 t2_bbbb(b,c,k,l) l2_bbbb(k,l,a,c)
        //             += +0.50 d_bb(i,j) r0 t2_1_abab(c,b,k,l) l2_1_abab(k,l,c,a)
        //             += +0.50 d_bb(i,j) r0 t2_1_abab(c,b,l,k) l2_1_abab(l,k,c,a)
        //             += +0.50 d_bb(i,j) r0_1 t2_abab(c,b,k,l) l2_1_abab(k,l,c,a)
        //             += +0.50 d_bb(i,j) r0_1 t2_abab(c,b,l,k) l2_1_abab(l,k,c,a)
        //             += +1.00 d_bb(i,j) r0 t1_1_bb(b,k) l1_1_bb(k,a)
        //             += +1.00 d_bb(i,j) r0 t1_bb(b,k) l1_bb(k,a)
        //             += +0.50 d_bb(i,j) r0 t2_1_bbbb(b,c,k,l) l2_1_bbbb(k,l,a,c)
        //             += +1.00 d_bb(i,j) r0_1 t1_bb(b,k) l1_1_bb(k,a)
        //             += -1.00 r0 t2_bbbb(b,c,k,i) l2_bbbb(k,j,a,c)
        //             += -1.00 r0 t2_abab(c,b,k,i) l2_abab(k,j,c,a)
        //             += -1.00 r0_1 t2_bbbb(b,c,k,i) l2_1_bbbb(k,j,a,c)
        //             += -1.00 r0_1 t2_abab(c,b,k,i) l2_1_abab(k,j,c,a)
        //             += -1.00 r0 t2_1_bbbb(b,c,k,i) l2_1_bbbb(k,j,a,c)
        D_ovvo_bbbb("L,R,i,a,b,j") -= tmps_["310_bbbb_LLvvoo"]("L,R,a,b,i,j");
        tmps_["310_bbbb_LLvvoo"].~TArrayD();

        // flops: o4v0L2  = o4v2L2
        //  mems: o4v0L2  = o4v0L2
        tmps_["311_abab_LLoooo"]("L,R,i,j,k,l")  = l2_1["abab"]("L,i,j,a,b") * r2["abab"]("R,a,b,k,l");

        // D_oovv_abab += -0.25 r2_abab(c,d,i,j) t2_1_abab(a,b,k,l) l2_1_abab(k,l,c,d)
        //             += -0.25 r2_abab(c,d,i,j) t2_1_abab(a,b,l,k) l2_1_abab(l,k,c,d)
        //             += -0.25 r2_abab(d,c,i,j) t2_1_abab(a,b,k,l) l2_1_abab(k,l,d,c)
        //             += -0.25 r2_abab(d,c,i,j) t2_1_abab(a,b,l,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2_1["abab"]("a,b,k,l") * tmps_["311_abab_LLoooo"]("L,R,k,l,i,j");

        // D_oovv_abab += -0.50 r2_abab(c,d,i,j) t1_aa(a,l) t1_1_bb(b,k) l2_1_abab(l,k,c,d)
        //             += -0.50 r2_abab(d,c,i,j) t1_aa(a,l) t1_1_bb(b,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o4v1L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1_1["bb"]("b,k") * tmps_["311_abab_LLoooo"]("L,R,l,k,i,j") * t1["aa"]("a,l");

        // D_oovv_abab += -0.50 r2_abab(c,d,i,j) t1_bb(b,l) t1_1_aa(a,k) l2_1_abab(k,l,c,d)
        //             += -0.50 r2_abab(d,c,i,j) t1_bb(b,l) t1_1_aa(a,k) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o4v1L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["bb"]("b,l") * tmps_["311_abab_LLoooo"]("L,R,k,l,i,j") * t1_1["aa"]("a,k");

        // flops: o3v1L1  = o3v2L1 o4v2L1
        //  mems: o3v1L1  = o3v1L1 o3v1L1
        tmps_["312_bbaa_Loovo"]("L,i,j,a,k")  = t1["bb"]("c,i") * l2_1["bbbb"]("L,l,j,c,b") * t2["abab"]("a,b,k,l");

        // D_oovv_abab += -1.00 r0_1 t2_abab(a,c,i,k) t1_bb(b,l) t1_bb(d,j) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["bb"]("b,l") * tmps_["312_bbaa_Loovo"]("L,j,l,a,i") * r0_1("R");

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["313_aabb_Lvooo"]("L,a,i,j,k")  = t2["abab"]("a,b,i,j") * l1_1["bb"]("L,k,b");

        // D_oovv_abab += +1.00 r1_1_bb(b,k) t2_abab(a,c,i,j) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["313_aabb_Lvooo"]("L,a,i,j,k") * r1_1["bb"]("R,b,k");

        // D_oovv_abab += +1.00 r0 t2_abab(a,c,i,j) t1_1_bb(b,k) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["313_aabb_Lvooo"]("L,a,i,j,k") * t1_1["bb"]("b,k") * r0("R");

        // D_oovv_abab += +1.00 r0_1 t2_abab(a,c,i,j) t1_bb(b,k) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["313_aabb_Lvooo"]("L,a,i,j,k") * t1["bb"]("b,k") * r0_1("R");

        // flops: o4v0L2  = o4v2L2 o4v2L2 o4v0L2
        //  mems: o4v0L2  = o4v0L2 o4v0L2 o4v0L2
        tmps_["314_abab_LLoooo"]("L,R,i,j,k,l")  = l2["abab"]("L,i,j,a,b") * r2["abab"]("R,a,b,k,l");
        tmps_["314_abab_LLoooo"]("L,R,i,j,k,l") += l2_1["abab"]("L,i,j,a,b") * r2_1["abab"]("R,a,b,k,l");

        // D_oooo_abab += -0.50 r2_abab(a,b,i,j) l2_abab(k,l,a,b)
        //             += -0.50 r2_abab(b,a,i,j) l2_abab(k,l,b,a)
        //             += -0.50 r2_1_abab(a,b,i,j) l2_1_abab(k,l,a,b)
        //             += -0.50 r2_1_abab(b,a,i,j) l2_1_abab(k,l,b,a)
        D_oooo_abab("L,R,i,j,k,l") -= tmps_["314_abab_LLoooo"]("L,R,k,l,i,j");

        // D_oovv_abab += -0.25 r2_abab(c,d,i,j) t2_abab(a,b,k,l) l2_abab(k,l,c,d)
        //             += -0.25 r2_abab(c,d,i,j) t2_abab(a,b,l,k) l2_abab(l,k,c,d)
        //             += -0.25 r2_abab(d,c,i,j) t2_abab(a,b,k,l) l2_abab(k,l,d,c)
        //             += -0.25 r2_abab(d,c,i,j) t2_abab(a,b,l,k) l2_abab(l,k,d,c)
        //             += -0.25 r2_1_abab(c,d,i,j) t2_abab(a,b,k,l) l2_1_abab(k,l,c,d)
        //             += -0.25 r2_1_abab(c,d,i,j) t2_abab(a,b,l,k) l2_1_abab(l,k,c,d)
        //             += -0.25 r2_1_abab(d,c,i,j) t2_abab(a,b,k,l) l2_1_abab(k,l,d,c)
        //             += -0.25 r2_1_abab(d,c,i,j) t2_abab(a,b,l,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,b,k,l") * tmps_["314_abab_LLoooo"]("L,R,k,l,i,j");

        // D_oovv_abab += -0.50 r2_abab(c,d,i,j) t1_aa(a,k) t1_bb(b,l) l2_abab(k,l,c,d)
        //             += -0.50 r2_abab(d,c,i,j) t1_aa(a,k) t1_bb(b,l) l2_abab(k,l,d,c)
        //             += -0.50 r2_1_abab(c,d,i,j) t1_aa(a,k) t1_bb(b,l) l2_1_abab(k,l,c,d)
        //             += -0.50 r2_1_abab(d,c,i,j) t1_aa(a,k) t1_bb(b,l) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o4v1L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["bb"]("b,l") * tmps_["314_abab_LLoooo"]("L,R,k,l,i,j") * t1["aa"]("a,k");

        // flops: o4v0L2  = o3v2L1 o4v1L2 o3v2L1 o4v1L2 o4v0L2
        //  mems: o4v0L2  = o3v1L1 o4v0L2 o3v1L1 o4v0L2 o4v0L2
        tmps_["315_abba_LLoooo"]("L,R,i,j,k,l")  = l2_1["abab"]("L,i,j,a,c") * t1["bb"]("c,k") * r1["aa"]("R,a,l");
        tmps_["315_abba_LLoooo"]("L,R,i,j,k,l") += l2_1["abab"]("L,i,j,d,b") * t1["aa"]("d,l") * r1["bb"]("R,b,k");

        // D_oovv_abab += -1.00 r1_aa(d,i) t1_aa(a,l) t1_bb(c,j) t1_1_bb(b,k) l2_1_abab(l,k,d,c)
        //             += -1.00 r1_bb(d,j) t1_aa(a,l) t1_aa(c,i) t1_1_bb(b,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o4v1L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1_1["bb"]("b,k") * tmps_["315_abba_LLoooo"]("L,R,l,k,j,i") * t1["aa"]("a,l");

        // D_oovv_abab += -1.00 r1_aa(d,i) t1_bb(b,l) t1_bb(c,j) t1_1_aa(a,k) l2_1_abab(k,l,d,c)
        //             += -1.00 r1_bb(d,j) t1_bb(b,l) t1_aa(c,i) t1_1_aa(a,k) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o4v1L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["bb"]("b,l") * tmps_["315_abba_LLoooo"]("L,R,k,l,j,i") * t1_1["aa"]("a,k");

        // flops: o4v0L2  = o3v2L1 o4v1L2 o4v1L2 o3v2L1 o4v1L2 o4v0L2 o3v2L1 o4v1L2 o4v0L2 o4v0L2 o3v2L1 o4v1L2 o4v0L2 o3v2L1 o4v1L2 o4v0L2
        //  mems: o4v0L2  = o3v1L1 o4v0L2 o4v0L2 o3v1L1 o4v0L2 o4v0L2 o3v1L1 o4v0L2 o4v0L2 o4v0L2 o3v1L1 o4v0L2 o4v0L2 o3v1L1 o4v0L2 o4v0L2
        tmps_["316_aabb_LLoooo"]("L,R,i,j,k,l")  = l2_1["abab"]("L,j,k,c,d") * t1_1["bb"]("d,l") * r1["aa"]("R,c,i");
        tmps_["316_aabb_LLoooo"]("L,R,i,j,k,l") += tmps_["264_aabb_Looov"]("L,i,j,k,a") * r1_1["bb"]("R,a,l");
        tmps_["316_aabb_LLoooo"]("L,R,i,j,k,l") += l2_1["abab"]("L,j,k,c,d") * t1["bb"]("d,l") * r1_1["aa"]("R,c,i");
        tmps_["316_aabb_LLoooo"]("L,R,i,j,k,l") += l2["abab"]("L,j,k,c,d") * t1["bb"]("d,l") * r1["aa"]("R,c,i");
        tmps_["316_aabb_LLoooo"]("L,R,i,j,k,l") += t1_1["aa"]("b,i") * l2_1["abab"]("L,j,k,b,a") * r1["bb"]("R,a,l");
        tmps_["316_aabb_LLoooo"]("L,R,i,j,k,l") += t1["aa"]("b,i") * l2["abab"]("L,j,k,b,a") * r1["bb"]("R,a,l");
        tmps_["264_aabb_Looov"].~TArrayD();

        // D_oooo_abab += -1.00 r1_1_bb(b,j) t1_aa(a,i) l2_1_abab(k,l,a,b)
        //             += -1.00 r1_1_aa(b,i) t1_bb(a,j) l2_1_abab(k,l,b,a)
        //             += -1.00 r1_aa(b,i) t1_bb(a,j) l2_abab(k,l,b,a)
        //             += -1.00 r1_aa(b,i) t1_1_bb(a,j) l2_1_abab(k,l,b,a)
        //             += -1.00 r1_bb(b,j) t1_1_aa(a,i) l2_1_abab(k,l,a,b)
        //             += -1.00 r1_bb(b,j) t1_aa(a,i) l2_abab(k,l,a,b)
        D_oooo_abab("L,R,i,j,k,l") -= tmps_["316_aabb_LLoooo"]("L,R,i,k,l,j");

        // D_oovv_abab += -1.00 r1_1_bb(d,j) t1_aa(a,k) t1_bb(b,l) t1_aa(c,i) l2_1_abab(k,l,c,d)
        //             += -1.00 r1_1_aa(d,i) t1_aa(a,k) t1_bb(b,l) t1_bb(c,j) l2_1_abab(k,l,d,c)
        //             += -1.00 r1_aa(d,i) t1_aa(a,k) t1_bb(b,l) t1_bb(c,j) l2_abab(k,l,d,c)
        //             += -1.00 r1_aa(d,i) t1_aa(a,k) t1_bb(b,l) t1_1_bb(c,j) l2_1_abab(k,l,d,c)
        //             += -1.00 r1_bb(d,j) t1_aa(a,k) t1_bb(b,l) t1_1_aa(c,i) l2_1_abab(k,l,c,d)
        //             += -1.00 r1_bb(d,j) t1_aa(a,k) t1_bb(b,l) t1_aa(c,i) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o4v1L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["bb"]("b,l") * tmps_["316_aabb_LLoooo"]("L,R,i,k,l,j") * t1["aa"]("a,k");

        // flops: o3v1L1  = o3v2L1 o3v2L1 o3v1L1
        //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1
        tmps_["317_aabb_Lvooo"]("L,a,i,j,k")  = t2["abab"]("a,b,i,j") * l1["bb"]("L,k,b");
        tmps_["317_aabb_Lvooo"]("L,a,i,j,k") += t2_1["abab"]("a,b,i,j") * l1_1["bb"]("L,k,b");

        // D_oovv_abab += +1.00 r1_bb(b,k) t2_abab(a,c,i,j) l1_bb(k,c)
        //             += +1.00 r1_bb(b,k) t2_1_abab(a,c,i,j) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["317_aabb_Lvooo"]("L,a,i,j,k") * r1["bb"]("R,b,k");

        // D_oovv_abab += +1.00 r0 t2_abab(a,c,i,j) t1_bb(b,k) l1_bb(k,c)
        //             += +1.00 r0 t1_bb(b,k) t2_1_abab(a,c,i,j) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["317_aabb_Lvooo"]("L,a,i,j,k") * t1["bb"]("b,k") * r0("R");

        // flops: o3v1L1  = o3v2L1 o4v2L1 o3v2L1 o4v2L1 o3v2L1 o4v2L1 o3v1L1 o3v1L1
        //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1
        tmps_["318_bbaa_Loovo"]("L,i,j,a,k")  = -1.00 * t1["bb"]("b,i") * l2["bbbb"]("L,l,j,b,c") * t2["abab"]("a,c,k,l");
        tmps_["318_bbaa_Loovo"]("L,i,j,a,k") -= t1["bb"]("b,i") * l2_1["bbbb"]("L,l,j,b,c") * t2_1["abab"]("a,c,k,l");
        tmps_["318_bbaa_Loovo"]("L,i,j,a,k") += t1_1["bb"]("c,i") * l2_1["bbbb"]("L,l,j,b,c") * t2["abab"]("a,b,k,l");

        // D_oovv_abab += +1.00 r0 t2_abab(a,d,i,k) t1_bb(b,l) t1_1_bb(c,j) l2_1_bbbb(k,l,d,c)
        //             += -1.00 r0 t1_bb(b,l) t1_bb(d,j) t2_1_abab(a,c,i,k) l2_1_bbbb(k,l,d,c)
        //             += -1.00 r0 t2_abab(a,c,i,k) t1_bb(b,l) t1_bb(d,j) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1["bb"]("b,l") * tmps_["318_bbaa_Loovo"]("L,j,l,a,i") * r0("R");

        // flops: o3v1L1  = o3v2L1 o3v2L1 o3v1L1
        //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1
        tmps_["319_aabb_Lvooo"]("L,a,i,j,k")  = tmps_["16_aabb_Lvoov"]("L,a,i,j,b") * t1["bb"]("b,k");
        tmps_["319_aabb_Lvooo"]("L,a,i,j,k") += tmps_["13_abba_Lvoov"]("L,a,k,j,c") * t1["aa"]("c,i");

        // D_oovv_abab += -1.00 r1_1_bb(b,l) t2_aaaa(a,c,k,i) t1_bb(d,j) l2_1_abab(k,l,c,d)
        //             += -1.00 r1_1_bb(b,l) t2_abab(a,c,k,j) t1_aa(d,i) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["319_aabb_Lvooo"]("L,a,i,l,j") * r1_1["bb"]("R,b,l");

        // D_oovv_abab += -1.00 r0_1 t2_aaaa(a,c,k,i) t1_bb(b,l) t1_bb(d,j) l2_1_abab(k,l,c,d)
        //             += -1.00 r0_1 t2_abab(a,c,k,j) t1_bb(b,l) t1_aa(d,i) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["319_aabb_Lvooo"]("L,a,i,l,j") * t1["bb"]("b,l") * r0_1("R");

        // D_oovv_abab += -1.00 r0 t2_aaaa(a,c,l,i) t1_bb(d,j) t1_1_bb(b,k) l2_1_abab(l,k,c,d)
        //             += -1.00 r0 t2_abab(a,c,l,j) t1_aa(d,i) t1_1_bb(b,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["319_aabb_Lvooo"]("L,a,i,k,j") * t1_1["bb"]("b,k") * r0("R");

        // flops: o3v1L1  = o3v2L1 o3v2L1 o3v2L1 o3v2L1 o3v2L1 o3v1L1 o3v2L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1
        //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1
        tmps_["320_aabb_Lvooo"]("L,a,i,j,k")  = tmps_["17_aabb_Lvoov"]("L,a,i,j,b") * t1["bb"]("b,k");
        tmps_["320_aabb_Lvooo"]("L,a,i,j,k") += tmps_["13_abba_Lvoov"]("L,a,k,j,c") * t1_1["aa"]("c,i");
        tmps_["320_aabb_Lvooo"]("L,a,i,j,k") += tmps_["15_abba_Lvoov"]("L,a,k,j,d") * t1["aa"]("d,i");
        tmps_["320_aabb_Lvooo"]("L,a,i,j,k") += tmps_["16_aabb_Lvoov"]("L,a,i,j,e") * t1_1["bb"]("e,k");
        tmps_["320_aabb_Lvooo"]("L,a,i,j,k") += tmps_["18_aabb_Lvoov"]("L,a,i,j,b") * t1["bb"]("b,k");
        tmps_["320_aabb_Lvooo"]("L,a,i,j,k") += tmps_["14_abba_Lvoov"]("L,a,k,j,d") * t1["aa"]("d,i");
        tmps_["18_aabb_Lvoov"].~TArrayD();
        tmps_["17_aabb_Lvoov"].~TArrayD();
        tmps_["16_aabb_Lvoov"].~TArrayD();
        tmps_["15_abba_Lvoov"].~TArrayD();
        tmps_["14_abba_Lvoov"].~TArrayD();
        tmps_["13_abba_Lvoov"].~TArrayD();

        // D_oovv_abab += -1.00 r1_bb(b,l) t1_bb(d,j) t2_1_aaaa(a,c,k,i) l2_1_abab(k,l,c,d)
        //             += -1.00 r1_bb(b,l) t2_aaaa(a,d,k,i) t1_1_bb(c,j) l2_1_abab(k,l,d,c)
        //             += -1.00 r1_bb(b,l) t2_abab(a,c,k,j) t1_aa(d,i) l2_abab(k,l,d,c)
        //             += -1.00 r1_bb(b,l) t1_aa(d,i) t2_1_abab(a,c,k,j) l2_1_abab(k,l,d,c)
        //             += -1.00 r1_bb(b,l) t2_abab(a,d,k,j) t1_1_aa(c,i) l2_1_abab(k,l,c,d)
        //             += -1.00 r1_bb(b,l) t2_aaaa(a,c,k,i) t1_bb(d,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["320_aabb_Lvooo"]("L,a,i,l,j") * r1["bb"]("R,b,l");

        // D_oovv_abab += -1.00 r0 t1_bb(b,l) t1_bb(d,j) t2_1_aaaa(a,c,k,i) l2_1_abab(k,l,c,d)
        //             += -1.00 r0 t2_aaaa(a,d,k,i) t1_bb(b,l) t1_1_bb(c,j) l2_1_abab(k,l,d,c)
        //             += -1.00 r0 t2_abab(a,c,k,j) t1_bb(b,l) t1_aa(d,i) l2_abab(k,l,d,c)
        //             += -1.00 r0 t1_bb(b,l) t1_aa(d,i) t2_1_abab(a,c,k,j) l2_1_abab(k,l,d,c)
        //             += -1.00 r0 t2_abab(a,d,k,j) t1_bb(b,l) t1_1_aa(c,i) l2_1_abab(k,l,c,d)
        //             += -1.00 r0 t2_aaaa(a,c,k,i) t1_bb(b,l) t1_bb(d,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["320_aabb_Lvooo"]("L,a,i,l,j") * t1["bb"]("b,l") * r0("R");

        // flops: o3v1L2  = o4v1L1 o3v1L1 o3v1L2 o4v1L1 o3v1L2 o3v1L2 o4v1L2 o3v1L2 o3v1L2 o3v1L2 o4v1L2 o3v1L2 o4v1L2 o3v1L2 o3v1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o4v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1 o3v1L2 o3v1L2 o4v1L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o4v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o4v1L1 o3v1L2 o3v1L2 o4v1L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o4v1L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o4v1L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o4v1L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        //  mems: o3v1L2  = o3v1L1 o3v1L1 o3v1L2 o3v1L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1 o3v1L2 o3v1L2 o3v1L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L1 o3v1L2 o3v1L2 o3v1L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k")  = -1.00 * (t1["aa"]("a,l") * tmps_["267_abab_Loooo"]("L,i,j,l,k") + -1.00 * tmps_["318_bbaa_Loovo"]("L,j,k,a,i")) * r0("R");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= t1["aa"]("a,l") * tmps_["265_aabb_Loooo"]("L,i,l,k,j") * r0("R");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= t1_1["aa"]("a,l") * tmps_["311_abab_LLoooo"]("L,R,l,k,i,j");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += 0.50 * t1["aa"]("a,i") * tmps_["111_bb_LLoo"]("L,R,k,j");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= t1["aa"]("a,l") * tmps_["314_abab_LLoooo"]("L,R,l,k,i,j");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= r1["aa"]("R,a,l") * tmps_["267_abab_Loooo"]("L,i,j,l,k");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= Id["bb_oo"]("j,k") * t1["aa"]("a,i") * tmps_["229_LL"]("L,R");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += tmps_["108_bb_Loo"]("L,j,k") * tmps_["254_aa_Lvo"]("R,a,i");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= t1_1["aa"]("a,l") * tmps_["315_abba_LLoooo"]("L,R,l,k,j,i");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += tmps_["302_aa_LLvo"]("R,L,a,i") * Id["bb_oo"]("j,k");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= Id["bb_oo"]("j,k") * t1_1["aa"]("a,i") * tmps_["230_LL"]("L,R");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= t1_1["aa"]("a,l") * tmps_["272_abab_Loooo"]("L,i,j,l,k") * r0("R");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += tmps_["317_aabb_Lvooo"]("L,a,i,j,k") * r0("R");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += t1_1["aa"]("a,i") * tmps_["101_bb_LLoo"]("L,R,k,j");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= t1["aa"]("a,l") * tmps_["316_aabb_LLoooo"]("L,R,i,l,k,j");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= tmps_["227_bb_Loo"]("L,j,k") * tmps_["255_aa_Lvo"]("R,a,i");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= tmps_["227_bb_Loo"]("L,j,k") * tmps_["254_aa_Lvo"]("R,a,i");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += 0.50 * tmps_["107_bb_Loo"]("L,j,k") * tmps_["255_aa_Lvo"]("R,a,i");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += tmps_["96_bb_Loo"]("L,j,k") * tmps_["255_aa_Lvo"]("R,a,i");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += 0.50 * tmps_["107_bb_Loo"]("L,j,k") * tmps_["254_aa_Lvo"]("R,a,i");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += tmps_["108_bb_Loo"]("L,j,k") * tmps_["255_aa_Lvo"]("R,a,i");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= tmps_["228_bb_Loo"]("L,j,k") * tmps_["253_aa_Lvo"]("R,a,i");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += tmps_["96_bb_Loo"]("L,j,k") * tmps_["254_aa_Lvo"]("R,a,i");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= t1["aa"]("a,l") * tmps_["266_aabb_Loooo"]("L,i,l,k,j") * r0("R");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= t1["aa"]("a,l") * tmps_["272_abab_Loooo"]("L,i,j,l,k") * r0_1("R");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += tmps_["109_bb_Loo"]("L,j,k") * tmps_["253_aa_Lvo"]("R,a,i");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += tmps_["110_bb_Loo"]("L,j,k") * tmps_["253_aa_Lvo"]("R,a,i");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= r0("R") * tmps_["320_aabb_Lvooo"]("L,a,i,k,j");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += tmps_["305_aa_Lov"]("L,i,a") * tmps_["113_bb_Loo"]("R,j,k");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= t1["aa"]("a,i") * tmps_["114_bb_LLoo"]("L,R,k,j");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += Id["bb_oo"]("j,k") * tmps_["304_aa_LLvo"]("L,R,a,i");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= r0_1("R") * tmps_["319_aabb_Lvooo"]("L,a,i,k,j");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += t1_1["aa"]("a,i") * tmps_["125_bb_LLoo"]("L,R,k,j");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += 0.50 * tmps_["306_aa_Lov"]("L,i,a") * tmps_["112_bb_Loo"]("R,j,k");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= t1["aa"]("a,i") * tmps_["115_bb_LLoo"]("L,R,k,j");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += Id["bb_oo"]("j,k") * tmps_["307_aa_LLvo"]("L,R,a,i");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= r1_1["aa"]("R,a,i") * tmps_["227_bb_Loo"]("L,j,k");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += r0_1("R") * tmps_["313_aabb_Lvooo"]("L,a,i,j,k");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= t1_1["aa"]("a,l") * tmps_["273_baab_Loooo"]("L,j,i,l,k") * r0("R");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += r1_1["aa"]("R,a,i") * tmps_["96_bb_Loo"]("L,j,k");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= t1["aa"]("a,l") * tmps_["268_baab_Loooo"]("L,j,i,l,k") * r0("R");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += 0.50 * r1_1["aa"]("R,a,i") * tmps_["107_bb_Loo"]("L,j,k");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= r1["aa"]("R,a,i") * tmps_["228_bb_Loo"]("L,j,k");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += r1_1["aa"]("R,a,i") * tmps_["108_bb_Loo"]("L,j,k");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= t1["aa"]("a,l") * tmps_["273_baab_Loooo"]("L,j,i,l,k") * r0_1("R");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += r1["aa"]("R,a,i") * tmps_["109_bb_Loo"]("L,j,k");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") -= r0_1("R") * tmps_["312_bbaa_Loovo"]("L,j,k,a,i");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += r1["aa"]("R,a,i") * tmps_["110_bb_Loo"]("L,j,k");
        tmps_["321_aabb_LLovoo"]("R,L,i,a,j,k") += Id["bb_oo"]("j,k") * tmps_["308_aa_LLov"]("R,L,i,a");
        tmps_["320_aabb_Lvooo"].~TArrayD();
        tmps_["319_aabb_Lvooo"].~TArrayD();
        tmps_["318_bbaa_Loovo"].~TArrayD();
        tmps_["317_aabb_Lvooo"].~TArrayD();
        tmps_["316_aabb_LLoooo"].~TArrayD();
        tmps_["315_abba_LLoooo"].~TArrayD();
        tmps_["314_abab_LLoooo"].~TArrayD();
        tmps_["313_aabb_Lvooo"].~TArrayD();
        tmps_["312_bbaa_Loovo"].~TArrayD();
        tmps_["311_abab_LLoooo"].~TArrayD();
        tmps_["308_aa_LLov"].~TArrayD();
        tmps_["307_aa_LLvo"].~TArrayD();
        tmps_["306_aa_Lov"].~TArrayD();
        tmps_["305_aa_Lov"].~TArrayD();
        tmps_["304_aa_LLvo"].~TArrayD();
        tmps_["302_aa_LLvo"].~TArrayD();
        tmps_["273_baab_Loooo"].~TArrayD();
        tmps_["272_abab_Loooo"].~TArrayD();
        tmps_["268_baab_Loooo"].~TArrayD();
        tmps_["267_abab_Loooo"].~TArrayD();
        tmps_["266_aabb_Loooo"].~TArrayD();
        tmps_["265_aabb_Loooo"].~TArrayD();
        tmps_["255_aa_Lvo"].~TArrayD();
        tmps_["254_aa_Lvo"].~TArrayD();
        tmps_["253_aa_Lvo"].~TArrayD();
        tmps_["230_LL"].~TArrayD();
        tmps_["229_LL"].~TArrayD();
        tmps_["228_bb_Loo"].~TArrayD();
        tmps_["227_bb_Loo"].~TArrayD();
        tmps_["125_bb_LLoo"].~TArrayD();
        tmps_["115_bb_LLoo"].~TArrayD();
        tmps_["114_bb_LLoo"].~TArrayD();
        tmps_["113_bb_Loo"].~TArrayD();
        tmps_["112_bb_Loo"].~TArrayD();
        tmps_["111_bb_LLoo"].~TArrayD();
        tmps_["110_bb_Loo"].~TArrayD();
        tmps_["109_bb_Loo"].~TArrayD();
        tmps_["108_bb_Loo"].~TArrayD();
        tmps_["107_bb_Loo"].~TArrayD();
        tmps_["101_bb_LLoo"].~TArrayD();
        tmps_["96_bb_Loo"].~TArrayD();

        // D_ooov_abab += +0.50 d_bb(i,k) r1_aa(c,j) t2_1_abab(a,b,l,m) l2_1_abab(l,m,c,b)
        //             += +0.50 d_bb(i,k) r1_aa(c,j) t2_1_abab(a,b,m,l) l2_1_abab(m,l,c,b)
        //             += +0.50 d_bb(i,k) r1_1_aa(c,j) t2_abab(a,b,l,m) l2_1_abab(l,m,c,b)
        //             += +0.50 d_bb(i,k) r1_1_aa(c,j) t2_abab(a,b,m,l) l2_1_abab(m,l,c,b)
        //             += +1.00 r0_1 t2_abab(a,b,j,i) l1_1_bb(k,b)
        //             += -1.00 d_bb(i,k) r1_1_aa(a,j) l0_1
        //             += -1.00 r0 t1_bb(b,i) t1_aa(c,j) t1_1_aa(a,l) l2_1_abab(l,k,c,b)
        //             += +0.50 r1_1_aa(a,j) t2_abab(c,b,l,i) l2_1_abab(l,k,c,b)
        //             += +0.50 r1_1_aa(a,j) t2_abab(b,c,l,i) l2_1_abab(l,k,b,c)
        //             += +0.50 d_bb(i,k) r0_1 t2_abab(a,b,l,m) t1_aa(c,j) l2_1_abab(l,m,c,b)
        //             += +0.50 d_bb(i,k) r0_1 t2_abab(a,b,m,l) t1_aa(c,j) l2_1_abab(m,l,c,b)
        //             += +0.50 d_bb(i,k) r0_1 t1_aa(a,m) t2_aaaa(c,b,l,j) l2_1_aaaa(l,m,c,b)
        //             += +0.50 d_bb(i,k) r0_1 t1_aa(a,m) t2_abab(c,b,j,l) l2_1_abab(m,l,c,b)
        //             += +0.50 d_bb(i,k) r0_1 t1_aa(a,m) t2_abab(b,c,j,l) l2_1_abab(m,l,b,c)
        //             += +1.00 d_bb(i,k) r0_1 t1_aa(a,l) t1_aa(b,j) l1_1_aa(l,b)
        //             += +1.00 d_bb(i,k) r0_1 t2_aaaa(a,b,l,j) l1_1_aa(l,b)
        //             += -1.00 d_bb(i,k) r0_1 t2_abab(a,b,j,l) l1_1_bb(l,b)
        //             += +0.50 d_bb(i,k) r0_1 t2_aaaa(a,b,l,m) t1_aa(c,j) l2_1_aaaa(l,m,c,b)
        //             += -1.00 r0 t1_bb(c,i) t2_1_aaaa(a,b,l,j) l2_1_abab(l,k,b,c)
        //             += -1.00 r0 t2_aaaa(a,c,l,j) t1_1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += -1.00 r0 t2_abab(a,b,l,i) t1_aa(c,j) l2_abab(l,k,c,b)
        //             += -1.00 r0 t1_aa(c,j) t2_1_abab(a,b,l,i) l2_1_abab(l,k,c,b)
        //             += -1.00 r0 t2_abab(a,c,l,i) t1_1_aa(b,j) l2_1_abab(l,k,b,c)
        //             += -1.00 r0 t2_aaaa(a,b,l,j) t1_bb(c,i) l2_abab(l,k,b,c)
        //             += -0.50 r0 t1_aa(a,l) t2_abab(c,b,j,i) l2_abab(l,k,c,b)
        //             += -0.50 r0 t1_aa(a,l) t2_abab(b,c,j,i) l2_abab(l,k,b,c)
        //             += -0.50 r0 t1_aa(a,l) t2_1_abab(c,b,j,i) l2_1_abab(l,k,c,b)
        //             += -0.50 r0 t1_aa(a,l) t2_1_abab(b,c,j,i) l2_1_abab(l,k,b,c)
        //             += +1.00 r0 t2_abab(a,c,j,l) t1_1_bb(b,i) l2_1_bbbb(l,k,c,b)
        //             += -1.00 r0 t1_bb(c,i) t2_1_abab(a,b,j,l) l2_1_bbbb(l,k,c,b)
        //             += -1.00 r0 t2_abab(a,b,j,l) t1_bb(c,i) l2_bbbb(l,k,c,b)
        //             += -1.00 r0 t1_aa(a,l) t1_bb(c,i) t1_1_aa(b,j) l2_1_abab(l,k,b,c)
        //             += -0.50 r2_abab(b,c,j,i) t1_1_aa(a,l) l2_1_abab(l,k,b,c)
        //             += -0.50 r2_abab(c,b,j,i) t1_1_aa(a,l) l2_1_abab(l,k,c,b)
        //             += +0.50 r2_1_bbbb(b,c,l,i) t1_aa(a,j) l2_1_bbbb(l,k,b,c)
        //             += +0.50 r2_abab(b,c,l,i) t1_aa(a,j) l2_abab(l,k,b,c)
        //             += +0.50 r2_abab(c,b,l,i) t1_aa(a,j) l2_abab(l,k,c,b)
        //             += +0.50 r2_bbbb(b,c,l,i) t1_aa(a,j) l2_bbbb(l,k,b,c)
        //             += +1.00 r1_bb(b,i) t1_aa(a,j) l1_bb(k,b)
        //             += +0.50 r2_1_abab(b,c,l,i) t1_aa(a,j) l2_1_abab(l,k,b,c)
        //             += +0.50 r2_1_abab(c,b,l,i) t1_aa(a,j) l2_1_abab(l,k,c,b)
        //             += +1.00 r1_1_bb(b,i) t1_aa(a,j) l1_1_bb(k,b)
        //             += -0.50 r2_abab(b,c,j,i) t1_aa(a,l) l2_abab(l,k,b,c)
        //             += -0.50 r2_abab(c,b,j,i) t1_aa(a,l) l2_abab(l,k,c,b)
        //             += -0.50 r2_1_abab(b,c,j,i) t1_aa(a,l) l2_1_abab(l,k,b,c)
        //             += -0.50 r2_1_abab(c,b,j,i) t1_aa(a,l) l2_1_abab(l,k,c,b)
        //             += -0.50 r1_aa(a,l) t2_abab(c,b,j,i) l2_abab(l,k,c,b)
        //             += -0.50 r1_aa(a,l) t2_abab(b,c,j,i) l2_abab(l,k,b,c)
        //             += -0.50 r1_aa(a,l) t2_1_abab(c,b,j,i) l2_1_abab(l,k,c,b)
        //             += -0.50 r1_aa(a,l) t2_1_abab(b,c,j,i) l2_1_abab(l,k,b,c)
        //             += -0.25 d_bb(i,k) r2_1_abab(b,c,m,l) t1_aa(a,j) l2_1_abab(m,l,b,c)
        //             += -0.25 d_bb(i,k) r2_1_abab(b,c,l,m) t1_aa(a,j) l2_1_abab(l,m,b,c)
        //             += -0.25 d_bb(i,k) r2_1_abab(c,b,m,l) t1_aa(a,j) l2_1_abab(m,l,c,b)
        //             += -0.25 d_bb(i,k) r2_1_abab(c,b,l,m) t1_aa(a,j) l2_1_abab(l,m,c,b)
        //             += -0.25 d_bb(i,k) r2_aaaa(b,c,m,l) t1_aa(a,j) l2_aaaa(m,l,b,c)
        //             += -0.25 d_bb(i,k) r2_bbbb(b,c,m,l) t1_aa(a,j) l2_bbbb(m,l,b,c)
        //             += -1.00 d_bb(i,k) r1_1_aa(b,l) t1_aa(a,j) l1_1_aa(l,b)
        //             += -1.00 d_bb(i,k) r1_1_bb(b,l) t1_aa(a,j) l1_1_bb(l,b)
        //             += -1.00 d_bb(i,k) r1_aa(b,l) t1_aa(a,j) l1_aa(l,b)
        //             += -0.25 d_bb(i,k) r2_1_bbbb(b,c,m,l) t1_aa(a,j) l2_1_bbbb(m,l,b,c)
        //             += -1.00 d_bb(i,k) r1_bb(b,l) t1_aa(a,j) l1_bb(l,b)
        //             += -0.25 d_bb(i,k) r2_1_aaaa(b,c,m,l) t1_aa(a,j) l2_1_aaaa(m,l,b,c)
        //             += -0.25 d_bb(i,k) r2_abab(b,c,m,l) t1_aa(a,j) l2_abab(m,l,b,c)
        //             += -0.25 d_bb(i,k) r2_abab(b,c,l,m) t1_aa(a,j) l2_abab(l,m,b,c)
        //             += -0.25 d_bb(i,k) r2_abab(c,b,m,l) t1_aa(a,j) l2_abab(m,l,c,b)
        //             += -0.25 d_bb(i,k) r2_abab(c,b,l,m) t1_aa(a,j) l2_abab(l,m,c,b)
        //             += +1.00 r0_1 t1_aa(a,j) t1_bb(b,i) l1_1_bb(k,b)
        //             += -1.00 r1_aa(c,j) t1_bb(b,i) t1_1_aa(a,l) l2_1_abab(l,k,c,b)
        //             += -1.00 r1_bb(c,i) t1_aa(b,j) t1_1_aa(a,l) l2_1_abab(l,k,b,c)
        //             += +1.00 d_bb(i,k) r2_1_aaaa(a,b,l,j) l1_1_aa(l,b)
        //             += +1.00 d_bb(i,k) r1_aa(a,l) t1_1_aa(b,j) l1_1_aa(l,b)
        //             += +1.00 d_bb(i,k) r2_aaaa(a,b,l,j) l1_aa(l,b)
        //             += +1.00 d_bb(i,k) r1_bb(c,m) t2_1_abab(a,b,j,l) l2_1_bbbb(m,l,b,c)
        //             += +1.00 d_bb(i,k) r1_aa(a,l) t1_aa(b,j) l1_aa(l,b)
        //             += -1.00 d_bb(i,k) r2_1_abab(a,b,j,l) l1_1_bb(l,b)
        //             += -1.00 d_bb(i,k) r2_abab(a,b,j,l) l1_bb(l,b)
        //             += +1.00 d_bb(i,k) r1_1_aa(a,l) t1_aa(b,j) l1_1_aa(l,b)
        //             += -1.00 d_bb(i,k) r1_bb(b,l) t1_1_aa(a,j) l1_1_bb(l,b)
        //             += -1.00 d_bb(i,k) r1_aa(b,l) t1_1_aa(a,j) l1_1_aa(l,b)
        //             += -0.25 d_bb(i,k) r2_aaaa(b,c,m,l) t1_1_aa(a,j) l2_1_aaaa(m,l,b,c)
        //             += -0.25 d_bb(i,k) r2_bbbb(b,c,m,l) t1_1_aa(a,j) l2_1_bbbb(m,l,b,c)
        //             += -0.25 d_bb(i,k) r2_abab(b,c,m,l) t1_1_aa(a,j) l2_1_abab(m,l,b,c)
        //             += -0.25 d_bb(i,k) r2_abab(b,c,l,m) t1_1_aa(a,j) l2_1_abab(l,m,b,c)
        //             += -0.25 d_bb(i,k) r2_abab(c,b,m,l) t1_1_aa(a,j) l2_1_abab(m,l,c,b)
        //             += -0.25 d_bb(i,k) r2_abab(c,b,l,m) t1_1_aa(a,j) l2_1_abab(l,m,c,b)
        //             += -0.50 r0 t2_abab(c,b,j,i) t1_1_aa(a,l) l2_1_abab(l,k,c,b)
        //             += -0.50 r0 t2_abab(b,c,j,i) t1_1_aa(a,l) l2_1_abab(l,k,b,c)
        //             += +1.00 r0 t2_abab(a,b,j,i) l1_bb(k,b)
        //             += +1.00 r0 t2_1_abab(a,b,j,i) l1_1_bb(k,b)
        //             += +0.50 r2_abab(b,c,l,i) t1_1_aa(a,j) l2_1_abab(l,k,b,c)
        //             += +0.50 r2_abab(c,b,l,i) t1_1_aa(a,j) l2_1_abab(l,k,c,b)
        //             += +0.50 r2_bbbb(b,c,l,i) t1_1_aa(a,j) l2_1_bbbb(l,k,b,c)
        //             += +1.00 r1_bb(b,i) t1_1_aa(a,j) l1_1_bb(k,b)
        //             += -1.00 r1_1_bb(c,i) t1_aa(a,l) t1_aa(b,j) l2_1_abab(l,k,b,c)
        //             += -1.00 r1_1_aa(c,j) t1_aa(a,l) t1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += -1.00 r1_aa(c,j) t1_aa(a,l) t1_bb(b,i) l2_abab(l,k,c,b)
        //             += -1.00 r1_aa(c,j) t1_aa(a,l) t1_1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += -1.00 r1_bb(c,i) t1_aa(a,l) t1_1_aa(b,j) l2_1_abab(l,k,b,c)
        //             += -1.00 r1_bb(c,i) t1_aa(a,l) t1_aa(b,j) l2_abab(l,k,b,c)
        //             += -0.50 r1_1_aa(a,l) t2_abab(c,b,j,i) l2_1_abab(l,k,c,b)
        //             += -0.50 r1_1_aa(a,l) t2_abab(b,c,j,i) l2_1_abab(l,k,b,c)
        //             += -1.00 r1_aa(a,l) t1_bb(b,i) t1_aa(c,j) l2_abab(l,k,c,b)
        //             += -1.00 r1_1_aa(a,l) t1_bb(b,i) t1_aa(c,j) l2_1_abab(l,k,c,b)
        //             += -1.00 r1_aa(a,l) t1_bb(c,i) t1_1_aa(b,j) l2_1_abab(l,k,b,c)
        //             += -1.00 r1_aa(a,l) t1_aa(c,j) t1_1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += -1.00 d_bb(i,k) r0 t1_1_aa(a,j) l0_1
        //             += -1.00 d_bb(i,k) r0_1 t1_aa(a,j) l0_1
        //             += +0.50 r0 t2_bbbb(c,b,l,i) t1_1_aa(a,j) l2_1_bbbb(l,k,c,b)
        //             += +0.50 r0 t2_abab(c,b,l,i) t1_1_aa(a,j) l2_1_abab(l,k,c,b)
        //             += +0.50 r0 t2_abab(b,c,l,i) t1_1_aa(a,j) l2_1_abab(l,k,b,c)
        //             += +0.50 r0_1 t1_aa(a,j) t2_bbbb(c,b,l,i) l2_1_bbbb(l,k,c,b)
        //             += +1.00 r0 t1_bb(b,i) t1_1_aa(a,j) l1_1_bb(k,b)
        //             += -1.00 d_bb(i,k) r0 t1_aa(a,j) l0
        //             += +0.50 r0_1 t1_aa(a,j) t2_abab(c,b,l,i) l2_1_abab(l,k,c,b)
        //             += +0.50 r0_1 t1_aa(a,j) t2_abab(b,c,l,i) l2_1_abab(l,k,b,c)
        //             += -1.00 r0 t1_aa(a,l) t1_bb(b,i) t1_aa(c,j) l2_abab(l,k,c,b)
        //             += -0.50 r0_1 t1_aa(a,l) t2_abab(c,b,j,i) l2_1_abab(l,k,c,b)
        //             += -0.50 r0_1 t1_aa(a,l) t2_abab(b,c,j,i) l2_1_abab(l,k,b,c)
        //             += +1.00 r0 t1_aa(a,j) t1_bb(b,i) l1_bb(k,b)
        //             += +1.00 r0 t1_aa(a,j) t1_1_bb(b,i) l1_1_bb(k,b)
        //             += -1.00 r1_bb(c,l) t1_aa(a,j) t1_1_bb(b,i) l2_1_bbbb(l,k,b,c)
        //             += +1.00 r1_aa(c,l) t1_aa(a,j) t1_1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += +0.50 d_bb(i,k) r2_abab(a,c,m,l) t1_1_aa(b,j) l2_1_abab(m,l,b,c)
        //             += +0.50 d_bb(i,k) r2_abab(a,c,l,m) t1_1_aa(b,j) l2_1_abab(l,m,b,c)
        //             += +0.50 d_bb(i,k) r2_aaaa(a,c,m,l) t1_1_aa(b,j) l2_1_aaaa(m,l,b,c)
        //             += -0.50 d_bb(i,k) r1_aa(c,j) t2_aaaa(a,b,l,m) l2_aaaa(l,m,b,c)
        //             += +1.00 d_bb(i,k) r1_bb(c,m) t2_abab(a,b,j,l) l2_bbbb(m,l,b,c)
        //             += +1.00 d_bb(i,k) r1_bb(c,m) t2_aaaa(a,b,l,j) l2_abab(l,m,b,c)
        //             += -1.00 d_bb(i,k) r1_aa(c,m) t2_1_abab(a,b,j,l) l2_1_abab(m,l,c,b)
        //             += +0.50 d_bb(i,k) r2_aaaa(b,c,m,j) t1_aa(a,l) l2_aaaa(m,l,b,c)
        //             += +1.00 d_bb(i,k) r1_aa(b,j) t1_aa(a,l) l1_aa(l,b)
        //             += +0.50 d_bb(i,k) r2_1_aaaa(b,c,m,j) t1_aa(a,l) l2_1_aaaa(m,l,b,c)
        //             += +1.00 d_bb(i,k) r1_1_aa(b,j) t1_aa(a,l) l1_1_aa(l,b)
        //             += +0.50 d_bb(i,k) r2_1_abab(b,c,j,m) t1_aa(a,l) l2_1_abab(l,m,b,c)
        //             += +0.50 d_bb(i,k) r2_1_abab(c,b,j,m) t1_aa(a,l) l2_1_abab(l,m,c,b)
        //             += +0.50 d_bb(i,k) r2_abab(b,c,j,m) t1_aa(a,l) l2_abab(l,m,b,c)
        //             += +0.50 d_bb(i,k) r2_abab(c,b,j,m) t1_aa(a,l) l2_abab(l,m,c,b)
        //             += -0.50 d_bb(i,k) r1_aa(a,m) t2_aaaa(c,b,l,j) l2_aaaa(m,l,c,b)
        //             += -0.50 d_bb(i,k) r1_aa(a,m) t2_1_aaaa(c,b,l,j) l2_1_aaaa(m,l,c,b)
        //             += +0.50 d_bb(i,k) r1_aa(c,j) t2_abab(a,b,l,m) l2_abab(l,m,c,b)
        //             += +0.50 d_bb(i,k) r1_aa(c,j) t2_abab(a,b,m,l) l2_abab(m,l,c,b)
        //             += +1.00 d_bb(i,k) r1_1_bb(c,m) t2_abab(a,b,j,l) l2_1_bbbb(m,l,b,c)
        //             += -0.50 d_bb(i,k) r1_1_aa(a,m) t2_aaaa(c,b,l,j) l2_1_aaaa(m,l,c,b)
        //             += +1.00 d_bb(i,k) r1_1_bb(c,m) t2_aaaa(a,b,l,j) l2_1_abab(l,m,b,c)
        //             += -1.00 d_bb(i,k) r1_1_aa(c,m) t2_aaaa(a,b,l,j) l2_1_aaaa(m,l,b,c)
        //             += +1.00 d_bb(i,k) r1_bb(c,m) t2_1_aaaa(a,b,l,j) l2_1_abab(l,m,b,c)
        //             += -1.00 d_bb(i,k) r1_aa(c,m) t2_abab(a,b,j,l) l2_abab(m,l,c,b)
        //             += +0.50 d_bb(i,k) r2_abab(b,c,j,m) t1_1_aa(a,l) l2_1_abab(l,m,b,c)
        //             += +0.50 d_bb(i,k) r2_abab(c,b,j,m) t1_1_aa(a,l) l2_1_abab(l,m,c,b)
        //             += +1.00 d_bb(i,k) r1_aa(b,j) t1_1_aa(a,l) l1_1_aa(l,b)
        //             += +0.50 d_bb(i,k) r2_aaaa(b,c,m,j) t1_1_aa(a,l) l2_1_aaaa(m,l,b,c)
        //             += -0.50 d_bb(i,k) r1_1_aa(c,j) t2_aaaa(a,b,l,m) l2_1_aaaa(l,m,b,c)
        //             += -1.00 d_bb(i,k) r1_aa(c,m) t2_aaaa(a,b,l,j) l2_aaaa(m,l,b,c)
        //             += -1.00 d_bb(i,k) r1_aa(c,m) t2_1_aaaa(a,b,l,j) l2_1_aaaa(m,l,b,c)
        //             += -1.00 d_bb(i,k) r1_1_aa(c,m) t2_abab(a,b,j,l) l2_1_abab(m,l,c,b)
        //             += +0.50 d_bb(i,k) r1_1_aa(a,m) t2_abab(c,b,j,l) l2_1_abab(m,l,c,b)
        //             += +0.50 d_bb(i,k) r1_1_aa(a,m) t2_abab(b,c,j,l) l2_1_abab(m,l,b,c)
        //             += -0.50 d_bb(i,k) r1_aa(c,j) t2_1_aaaa(a,b,l,m) l2_1_aaaa(l,m,b,c)
        //             += +0.50 d_bb(i,k) r2_1_abab(a,c,m,l) t1_aa(b,j) l2_1_abab(m,l,b,c)
        //             += +0.50 d_bb(i,k) r2_1_abab(a,c,l,m) t1_aa(b,j) l2_1_abab(l,m,b,c)
        //             += +0.50 d_bb(i,k) r2_aaaa(a,c,m,l) t1_aa(b,j) l2_aaaa(m,l,b,c)
        //             += +0.50 d_bb(i,k) r2_abab(a,c,m,l) t1_aa(b,j) l2_abab(m,l,b,c)
        //             += +0.50 d_bb(i,k) r2_abab(a,c,l,m) t1_aa(b,j) l2_abab(l,m,b,c)
        //             += +0.50 d_bb(i,k) r2_1_aaaa(a,c,m,l) t1_aa(b,j) l2_1_aaaa(m,l,b,c)
        //             += -1.00 r0_1 t2_aaaa(a,b,l,j) t1_bb(c,i) l2_1_abab(l,k,b,c)
        //             += -1.00 r0_1 t2_abab(a,b,l,i) t1_aa(c,j) l2_1_abab(l,k,c,b)
        //             += +1.00 r1_aa(c,l) t1_bb(b,i) t1_1_aa(a,j) l2_1_abab(l,k,c,b)
        //             += -1.00 r1_bb(c,l) t1_bb(b,i) t1_1_aa(a,j) l2_1_bbbb(l,k,b,c)
        //             += +0.50 d_bb(i,k) r0 t1_aa(c,j) t2_1_aaaa(a,b,l,m) l2_1_aaaa(l,m,c,b)
        //             += -0.50 d_bb(i,k) r0 t2_aaaa(a,c,l,m) t1_1_aa(b,j) l2_1_aaaa(l,m,c,b)
        //             += +0.50 d_bb(i,k) r0 t2_abab(a,c,l,m) t1_1_aa(b,j) l2_1_abab(l,m,b,c)
        //             += +0.50 d_bb(i,k) r0 t2_abab(a,c,m,l) t1_1_aa(b,j) l2_1_abab(m,l,b,c)
        //             += +1.00 d_bb(i,k) r0 t1_aa(a,l) t1_aa(b,j) l1_aa(l,b)
        //             += -1.00 d_bb(i,k) r0 t2_abab(a,b,j,l) l1_bb(l,b)
        //             += +1.00 d_bb(i,k) r0 t1_aa(a,l) t1_1_aa(b,j) l1_1_aa(l,b)
        //             += +1.00 d_bb(i,k) r0 t2_aaaa(a,b,l,j) l1_aa(l,b)
        //             += +1.00 d_bb(i,k) r0 t1_aa(b,j) t1_1_aa(a,l) l1_1_aa(l,b)
        //             += -1.00 d_bb(i,k) r0 t2_1_abab(a,b,j,l) l1_1_bb(l,b)
        //             += +1.00 d_bb(i,k) r0 t2_1_aaaa(a,b,l,j) l1_1_aa(l,b)
        //             += -0.50 d_bb(i,k) r0 t2_aaaa(c,b,m,j) t1_1_aa(a,l) l2_1_aaaa(l,m,c,b)
        //             += +0.50 d_bb(i,k) r0 t2_abab(a,b,l,m) t1_aa(c,j) l2_abab(l,m,c,b)
        //             += +0.50 d_bb(i,k) r0 t2_abab(a,b,m,l) t1_aa(c,j) l2_abab(m,l,c,b)
        //             += +0.50 d_bb(i,k) r0 t2_aaaa(a,b,l,m) t1_aa(c,j) l2_aaaa(l,m,c,b)
        //             += +0.50 d_bb(i,k) r0 t1_aa(c,j) t2_1_abab(a,b,l,m) l2_1_abab(l,m,c,b)
        //             += +0.50 d_bb(i,k) r0 t1_aa(c,j) t2_1_abab(a,b,m,l) l2_1_abab(m,l,c,b)
        //             += +0.50 d_bb(i,k) r0 t2_abab(c,b,j,m) t1_1_aa(a,l) l2_1_abab(l,m,c,b)
        //             += +0.50 d_bb(i,k) r0 t2_abab(b,c,j,m) t1_1_aa(a,l) l2_1_abab(l,m,b,c)
        //             += -1.00 r1_bb(c,l) t1_aa(a,j) t1_bb(b,i) l2_bbbb(l,k,b,c)
        //             += +1.00 r1_aa(c,l) t1_aa(a,j) t1_bb(b,i) l2_abab(l,k,c,b)
        //             += +1.00 r1_1_aa(c,l) t1_aa(a,j) t1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += -1.00 r1_1_bb(c,l) t1_aa(a,j) t1_bb(b,i) l2_1_bbbb(l,k,b,c)
        //             += +1.00 d_bb(i,k) r1_bb(c,m) t1_aa(b,j) t1_1_aa(a,l) l2_1_abab(l,m,b,c)
        //             += -1.00 d_bb(i,k) r1_aa(c,m) t1_aa(b,j) t1_1_aa(a,l) l2_1_aaaa(m,l,b,c)
        //             += +1.00 d_bb(i,k) r1_bb(c,m) t1_aa(a,l) t1_aa(b,j) l2_abab(l,m,b,c)
        //             += +1.00 d_bb(i,k) r1_1_bb(c,m) t1_aa(a,l) t1_aa(b,j) l2_1_abab(l,m,b,c)
        //             += -1.00 d_bb(i,k) r1_aa(c,m) t1_aa(a,l) t1_aa(b,j) l2_aaaa(m,l,b,c)
        //             += -1.00 d_bb(i,k) r1_1_aa(c,m) t1_aa(a,l) t1_aa(b,j) l2_1_aaaa(m,l,b,c)
        //             += -1.00 d_bb(i,k) r1_aa(c,m) t1_aa(a,l) t1_1_aa(b,j) l2_1_aaaa(m,l,b,c)
        //             += +1.00 d_bb(i,k) r1_bb(c,m) t1_aa(a,l) t1_1_aa(b,j) l2_1_abab(l,m,b,c)
        //             += -1.00 r0 t1_aa(a,l) t1_aa(c,j) t1_1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += +0.50 r1_1_aa(a,j) t2_bbbb(c,b,l,i) l2_1_bbbb(l,k,c,b)
        //             += -1.00 d_bb(i,k) r1_aa(a,j) l0
        //             += +1.00 r1_1_aa(a,j) t1_bb(b,i) l1_1_bb(k,b)
        //             += -1.00 r0_1 t1_aa(a,l) t1_bb(b,i) t1_aa(c,j) l2_1_abab(l,k,c,b)
        //             += +1.00 r1_aa(a,j) t1_bb(b,i) l1_bb(k,b)
        //             += -1.00 r0_1 t2_abab(a,b,j,l) t1_bb(c,i) l2_1_bbbb(l,k,c,b)
        //             += +1.00 r1_aa(a,j) t1_1_bb(b,i) l1_1_bb(k,b)
        D_ooov_abab("L,R,i,j,k,a") += tmps_["321_aabb_LLovoo"]("R,L,j,a,i,k");

        // D_oovo_abab += -0.50 d_bb(i,k) r1_aa(c,j) t2_1_abab(a,b,l,m) l2_1_abab(l,m,c,b)
        //             += -0.50 d_bb(i,k) r1_aa(c,j) t2_1_abab(a,b,m,l) l2_1_abab(m,l,c,b)
        //             += -0.50 d_bb(i,k) r1_1_aa(c,j) t2_abab(a,b,l,m) l2_1_abab(l,m,c,b)
        //             += -0.50 d_bb(i,k) r1_1_aa(c,j) t2_abab(a,b,m,l) l2_1_abab(m,l,c,b)
        //             += -1.00 r0_1 t2_abab(a,b,j,i) l1_1_bb(k,b)
        //             += +1.00 d_bb(i,k) r1_1_aa(a,j) l0_1
        //             += +1.00 r0 t1_bb(b,i) t1_aa(c,j) t1_1_aa(a,l) l2_1_abab(l,k,c,b)
        //             += -0.50 r1_1_aa(a,j) t2_abab(c,b,l,i) l2_1_abab(l,k,c,b)
        //             += -0.50 r1_1_aa(a,j) t2_abab(b,c,l,i) l2_1_abab(l,k,b,c)
        //             += -0.50 d_bb(i,k) r0_1 t2_abab(a,b,l,m) t1_aa(c,j) l2_1_abab(l,m,c,b)
        //             += -0.50 d_bb(i,k) r0_1 t2_abab(a,b,m,l) t1_aa(c,j) l2_1_abab(m,l,c,b)
        //             += -0.50 d_bb(i,k) r0_1 t1_aa(a,m) t2_aaaa(c,b,l,j) l2_1_aaaa(l,m,c,b)
        //             += -0.50 d_bb(i,k) r0_1 t1_aa(a,m) t2_abab(c,b,j,l) l2_1_abab(m,l,c,b)
        //             += -0.50 d_bb(i,k) r0_1 t1_aa(a,m) t2_abab(b,c,j,l) l2_1_abab(m,l,b,c)
        //             += -1.00 d_bb(i,k) r0_1 t1_aa(a,l) t1_aa(b,j) l1_1_aa(l,b)
        //             += -1.00 d_bb(i,k) r0_1 t2_aaaa(a,b,l,j) l1_1_aa(l,b)
        //             += +1.00 d_bb(i,k) r0_1 t2_abab(a,b,j,l) l1_1_bb(l,b)
        //             += -0.50 d_bb(i,k) r0_1 t2_aaaa(a,b,l,m) t1_aa(c,j) l2_1_aaaa(l,m,c,b)
        //             += +1.00 r0 t1_bb(c,i) t2_1_aaaa(a,b,l,j) l2_1_abab(l,k,b,c)
        //             += +1.00 r0 t2_aaaa(a,c,l,j) t1_1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += +1.00 r0 t2_abab(a,b,l,i) t1_aa(c,j) l2_abab(l,k,c,b)
        //             += +1.00 r0 t1_aa(c,j) t2_1_abab(a,b,l,i) l2_1_abab(l,k,c,b)
        //             += +1.00 r0 t2_abab(a,c,l,i) t1_1_aa(b,j) l2_1_abab(l,k,b,c)
        //             += +1.00 r0 t2_aaaa(a,b,l,j) t1_bb(c,i) l2_abab(l,k,b,c)
        //             += +0.50 r0 t1_aa(a,l) t2_abab(c,b,j,i) l2_abab(l,k,c,b)
        //             += +0.50 r0 t1_aa(a,l) t2_abab(b,c,j,i) l2_abab(l,k,b,c)
        //             += +0.50 r0 t1_aa(a,l) t2_1_abab(c,b,j,i) l2_1_abab(l,k,c,b)
        //             += +0.50 r0 t1_aa(a,l) t2_1_abab(b,c,j,i) l2_1_abab(l,k,b,c)
        //             += -1.00 r0 t2_abab(a,c,j,l) t1_1_bb(b,i) l2_1_bbbb(l,k,c,b)
        //             += +1.00 r0 t1_bb(c,i) t2_1_abab(a,b,j,l) l2_1_bbbb(l,k,c,b)
        //             += +1.00 r0 t2_abab(a,b,j,l) t1_bb(c,i) l2_bbbb(l,k,c,b)
        //             += +1.00 r0 t1_aa(a,l) t1_bb(c,i) t1_1_aa(b,j) l2_1_abab(l,k,b,c)
        //             += +0.50 r2_abab(b,c,j,i) t1_1_aa(a,l) l2_1_abab(l,k,b,c)
        //             += +0.50 r2_abab(c,b,j,i) t1_1_aa(a,l) l2_1_abab(l,k,c,b)
        //             += -0.50 r2_1_bbbb(b,c,l,i) t1_aa(a,j) l2_1_bbbb(l,k,b,c)
        //             += -0.50 r2_abab(b,c,l,i) t1_aa(a,j) l2_abab(l,k,b,c)
        //             += -0.50 r2_abab(c,b,l,i) t1_aa(a,j) l2_abab(l,k,c,b)
        //             += -0.50 r2_bbbb(b,c,l,i) t1_aa(a,j) l2_bbbb(l,k,b,c)
        //             += -1.00 r1_bb(b,i) t1_aa(a,j) l1_bb(k,b)
        //             += -0.50 r2_1_abab(b,c,l,i) t1_aa(a,j) l2_1_abab(l,k,b,c)
        //             += -0.50 r2_1_abab(c,b,l,i) t1_aa(a,j) l2_1_abab(l,k,c,b)
        //             += -1.00 r1_1_bb(b,i) t1_aa(a,j) l1_1_bb(k,b)
        //             += +0.50 r2_abab(b,c,j,i) t1_aa(a,l) l2_abab(l,k,b,c)
        //             += +0.50 r2_abab(c,b,j,i) t1_aa(a,l) l2_abab(l,k,c,b)
        //             += +0.50 r2_1_abab(b,c,j,i) t1_aa(a,l) l2_1_abab(l,k,b,c)
        //             += +0.50 r2_1_abab(c,b,j,i) t1_aa(a,l) l2_1_abab(l,k,c,b)
        //             += +0.50 r1_aa(a,l) t2_abab(c,b,j,i) l2_abab(l,k,c,b)
        //             += +0.50 r1_aa(a,l) t2_abab(b,c,j,i) l2_abab(l,k,b,c)
        //             += +0.50 r1_aa(a,l) t2_1_abab(c,b,j,i) l2_1_abab(l,k,c,b)
        //             += +0.50 r1_aa(a,l) t2_1_abab(b,c,j,i) l2_1_abab(l,k,b,c)
        //             += +0.25 d_bb(i,k) r2_1_abab(b,c,m,l) t1_aa(a,j) l2_1_abab(m,l,b,c)
        //             += +0.25 d_bb(i,k) r2_1_abab(b,c,l,m) t1_aa(a,j) l2_1_abab(l,m,b,c)
        //             += +0.25 d_bb(i,k) r2_1_abab(c,b,m,l) t1_aa(a,j) l2_1_abab(m,l,c,b)
        //             += +0.25 d_bb(i,k) r2_1_abab(c,b,l,m) t1_aa(a,j) l2_1_abab(l,m,c,b)
        //             += +0.25 d_bb(i,k) r2_aaaa(b,c,m,l) t1_aa(a,j) l2_aaaa(m,l,b,c)
        //             += +0.25 d_bb(i,k) r2_bbbb(b,c,m,l) t1_aa(a,j) l2_bbbb(m,l,b,c)
        //             += +1.00 d_bb(i,k) r1_1_aa(b,l) t1_aa(a,j) l1_1_aa(l,b)
        //             += +1.00 d_bb(i,k) r1_1_bb(b,l) t1_aa(a,j) l1_1_bb(l,b)
        //             += +1.00 d_bb(i,k) r1_aa(b,l) t1_aa(a,j) l1_aa(l,b)
        //             += +0.25 d_bb(i,k) r2_1_bbbb(b,c,m,l) t1_aa(a,j) l2_1_bbbb(m,l,b,c)
        //             += +1.00 d_bb(i,k) r1_bb(b,l) t1_aa(a,j) l1_bb(l,b)
        //             += +0.25 d_bb(i,k) r2_1_aaaa(b,c,m,l) t1_aa(a,j) l2_1_aaaa(m,l,b,c)
        //             += +0.25 d_bb(i,k) r2_abab(b,c,m,l) t1_aa(a,j) l2_abab(m,l,b,c)
        //             += +0.25 d_bb(i,k) r2_abab(b,c,l,m) t1_aa(a,j) l2_abab(l,m,b,c)
        //             += +0.25 d_bb(i,k) r2_abab(c,b,m,l) t1_aa(a,j) l2_abab(m,l,c,b)
        //             += +0.25 d_bb(i,k) r2_abab(c,b,l,m) t1_aa(a,j) l2_abab(l,m,c,b)
        //             += -1.00 r0_1 t1_aa(a,j) t1_bb(b,i) l1_1_bb(k,b)
        //             += +1.00 r1_aa(c,j) t1_bb(b,i) t1_1_aa(a,l) l2_1_abab(l,k,c,b)
        //             += +1.00 r1_bb(c,i) t1_aa(b,j) t1_1_aa(a,l) l2_1_abab(l,k,b,c)
        //             += -1.00 d_bb(i,k) r2_1_aaaa(a,b,l,j) l1_1_aa(l,b)
        //             += -1.00 d_bb(i,k) r1_aa(a,l) t1_1_aa(b,j) l1_1_aa(l,b)
        //             += -1.00 d_bb(i,k) r2_aaaa(a,b,l,j) l1_aa(l,b)
        //             += -1.00 d_bb(i,k) r1_bb(c,m) t2_1_abab(a,b,j,l) l2_1_bbbb(m,l,b,c)
        //             += -1.00 d_bb(i,k) r1_aa(a,l) t1_aa(b,j) l1_aa(l,b)
        //             += +1.00 d_bb(i,k) r2_1_abab(a,b,j,l) l1_1_bb(l,b)
        //             += +1.00 d_bb(i,k) r2_abab(a,b,j,l) l1_bb(l,b)
        //             += -1.00 d_bb(i,k) r1_1_aa(a,l) t1_aa(b,j) l1_1_aa(l,b)
        //             += +1.00 d_bb(i,k) r1_bb(b,l) t1_1_aa(a,j) l1_1_bb(l,b)
        //             += +1.00 d_bb(i,k) r1_aa(b,l) t1_1_aa(a,j) l1_1_aa(l,b)
        //             += +0.25 d_bb(i,k) r2_aaaa(b,c,m,l) t1_1_aa(a,j) l2_1_aaaa(m,l,b,c)
        //             += +0.25 d_bb(i,k) r2_bbbb(b,c,m,l) t1_1_aa(a,j) l2_1_bbbb(m,l,b,c)
        //             += +0.25 d_bb(i,k) r2_abab(b,c,m,l) t1_1_aa(a,j) l2_1_abab(m,l,b,c)
        //             += +0.25 d_bb(i,k) r2_abab(b,c,l,m) t1_1_aa(a,j) l2_1_abab(l,m,b,c)
        //             += +0.25 d_bb(i,k) r2_abab(c,b,m,l) t1_1_aa(a,j) l2_1_abab(m,l,c,b)
        //             += +0.25 d_bb(i,k) r2_abab(c,b,l,m) t1_1_aa(a,j) l2_1_abab(l,m,c,b)
        //             += +0.50 r0 t2_abab(c,b,j,i) t1_1_aa(a,l) l2_1_abab(l,k,c,b)
        //             += +0.50 r0 t2_abab(b,c,j,i) t1_1_aa(a,l) l2_1_abab(l,k,b,c)
        //             += -1.00 r0 t2_abab(a,b,j,i) l1_bb(k,b)
        //             += -1.00 r0 t2_1_abab(a,b,j,i) l1_1_bb(k,b)
        //             += -0.50 r2_abab(b,c,l,i) t1_1_aa(a,j) l2_1_abab(l,k,b,c)
        //             += -0.50 r2_abab(c,b,l,i) t1_1_aa(a,j) l2_1_abab(l,k,c,b)
        //             += -0.50 r2_bbbb(b,c,l,i) t1_1_aa(a,j) l2_1_bbbb(l,k,b,c)
        //             += -1.00 r1_bb(b,i) t1_1_aa(a,j) l1_1_bb(k,b)
        //             += +1.00 r1_1_bb(c,i) t1_aa(a,l) t1_aa(b,j) l2_1_abab(l,k,b,c)
        //             += +1.00 r1_1_aa(c,j) t1_aa(a,l) t1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += +1.00 r1_aa(c,j) t1_aa(a,l) t1_bb(b,i) l2_abab(l,k,c,b)
        //             += +1.00 r1_aa(c,j) t1_aa(a,l) t1_1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += +1.00 r1_bb(c,i) t1_aa(a,l) t1_1_aa(b,j) l2_1_abab(l,k,b,c)
        //             += +1.00 r1_bb(c,i) t1_aa(a,l) t1_aa(b,j) l2_abab(l,k,b,c)
        //             += +0.50 r1_1_aa(a,l) t2_abab(c,b,j,i) l2_1_abab(l,k,c,b)
        //             += +0.50 r1_1_aa(a,l) t2_abab(b,c,j,i) l2_1_abab(l,k,b,c)
        //             += +1.00 r1_aa(a,l) t1_bb(b,i) t1_aa(c,j) l2_abab(l,k,c,b)
        //             += +1.00 r1_1_aa(a,l) t1_bb(b,i) t1_aa(c,j) l2_1_abab(l,k,c,b)
        //             += +1.00 r1_aa(a,l) t1_bb(c,i) t1_1_aa(b,j) l2_1_abab(l,k,b,c)
        //             += +1.00 r1_aa(a,l) t1_aa(c,j) t1_1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += +1.00 d_bb(i,k) r0 t1_1_aa(a,j) l0_1
        //             += +1.00 d_bb(i,k) r0_1 t1_aa(a,j) l0_1
        //             += -0.50 r0 t2_bbbb(c,b,l,i) t1_1_aa(a,j) l2_1_bbbb(l,k,c,b)
        //             += -0.50 r0 t2_abab(c,b,l,i) t1_1_aa(a,j) l2_1_abab(l,k,c,b)
        //             += -0.50 r0 t2_abab(b,c,l,i) t1_1_aa(a,j) l2_1_abab(l,k,b,c)
        //             += -0.50 r0_1 t1_aa(a,j) t2_bbbb(c,b,l,i) l2_1_bbbb(l,k,c,b)
        //             += -1.00 r0 t1_bb(b,i) t1_1_aa(a,j) l1_1_bb(k,b)
        //             += +1.00 d_bb(i,k) r0 t1_aa(a,j) l0
        //             += -0.50 r0_1 t1_aa(a,j) t2_abab(c,b,l,i) l2_1_abab(l,k,c,b)
        //             += -0.50 r0_1 t1_aa(a,j) t2_abab(b,c,l,i) l2_1_abab(l,k,b,c)
        //             += +1.00 r0 t1_aa(a,l) t1_bb(b,i) t1_aa(c,j) l2_abab(l,k,c,b)
        //             += +0.50 r0_1 t1_aa(a,l) t2_abab(c,b,j,i) l2_1_abab(l,k,c,b)
        //             += +0.50 r0_1 t1_aa(a,l) t2_abab(b,c,j,i) l2_1_abab(l,k,b,c)
        //             += -1.00 r0 t1_aa(a,j) t1_bb(b,i) l1_bb(k,b)
        //             += -1.00 r0 t1_aa(a,j) t1_1_bb(b,i) l1_1_bb(k,b)
        //             += +1.00 r1_bb(c,l) t1_aa(a,j) t1_1_bb(b,i) l2_1_bbbb(l,k,b,c)
        //             += -1.00 r1_aa(c,l) t1_aa(a,j) t1_1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += -0.50 d_bb(i,k) r2_abab(a,c,m,l) t1_1_aa(b,j) l2_1_abab(m,l,b,c)
        //             += -0.50 d_bb(i,k) r2_abab(a,c,l,m) t1_1_aa(b,j) l2_1_abab(l,m,b,c)
        //             += -0.50 d_bb(i,k) r2_aaaa(a,c,m,l) t1_1_aa(b,j) l2_1_aaaa(m,l,b,c)
        //             += +0.50 d_bb(i,k) r1_aa(c,j) t2_aaaa(a,b,l,m) l2_aaaa(l,m,b,c)
        //             += -1.00 d_bb(i,k) r1_bb(c,m) t2_abab(a,b,j,l) l2_bbbb(m,l,b,c)
        //             += -1.00 d_bb(i,k) r1_bb(c,m) t2_aaaa(a,b,l,j) l2_abab(l,m,b,c)
        //             += +1.00 d_bb(i,k) r1_aa(c,m) t2_1_abab(a,b,j,l) l2_1_abab(m,l,c,b)
        //             += -0.50 d_bb(i,k) r2_aaaa(b,c,m,j) t1_aa(a,l) l2_aaaa(m,l,b,c)
        //             += -1.00 d_bb(i,k) r1_aa(b,j) t1_aa(a,l) l1_aa(l,b)
        //             += -0.50 d_bb(i,k) r2_1_aaaa(b,c,m,j) t1_aa(a,l) l2_1_aaaa(m,l,b,c)
        //             += -1.00 d_bb(i,k) r1_1_aa(b,j) t1_aa(a,l) l1_1_aa(l,b)
        //             += -0.50 d_bb(i,k) r2_1_abab(b,c,j,m) t1_aa(a,l) l2_1_abab(l,m,b,c)
        //             += -0.50 d_bb(i,k) r2_1_abab(c,b,j,m) t1_aa(a,l) l2_1_abab(l,m,c,b)
        //             += -0.50 d_bb(i,k) r2_abab(b,c,j,m) t1_aa(a,l) l2_abab(l,m,b,c)
        //             += -0.50 d_bb(i,k) r2_abab(c,b,j,m) t1_aa(a,l) l2_abab(l,m,c,b)
        //             += +0.50 d_bb(i,k) r1_aa(a,m) t2_aaaa(c,b,l,j) l2_aaaa(m,l,c,b)
        //             += +0.50 d_bb(i,k) r1_aa(a,m) t2_1_aaaa(c,b,l,j) l2_1_aaaa(m,l,c,b)
        //             += -0.50 d_bb(i,k) r1_aa(c,j) t2_abab(a,b,l,m) l2_abab(l,m,c,b)
        //             += -0.50 d_bb(i,k) r1_aa(c,j) t2_abab(a,b,m,l) l2_abab(m,l,c,b)
        //             += -1.00 d_bb(i,k) r1_1_bb(c,m) t2_abab(a,b,j,l) l2_1_bbbb(m,l,b,c)
        //             += +0.50 d_bb(i,k) r1_1_aa(a,m) t2_aaaa(c,b,l,j) l2_1_aaaa(m,l,c,b)
        //             += -1.00 d_bb(i,k) r1_1_bb(c,m) t2_aaaa(a,b,l,j) l2_1_abab(l,m,b,c)
        //             += +1.00 d_bb(i,k) r1_1_aa(c,m) t2_aaaa(a,b,l,j) l2_1_aaaa(m,l,b,c)
        //             += -1.00 d_bb(i,k) r1_bb(c,m) t2_1_aaaa(a,b,l,j) l2_1_abab(l,m,b,c)
        //             += +1.00 d_bb(i,k) r1_aa(c,m) t2_abab(a,b,j,l) l2_abab(m,l,c,b)
        //             += -0.50 d_bb(i,k) r2_abab(b,c,j,m) t1_1_aa(a,l) l2_1_abab(l,m,b,c)
        //             += -0.50 d_bb(i,k) r2_abab(c,b,j,m) t1_1_aa(a,l) l2_1_abab(l,m,c,b)
        //             += -1.00 d_bb(i,k) r1_aa(b,j) t1_1_aa(a,l) l1_1_aa(l,b)
        //             += -0.50 d_bb(i,k) r2_aaaa(b,c,m,j) t1_1_aa(a,l) l2_1_aaaa(m,l,b,c)
        //             += +0.50 d_bb(i,k) r1_1_aa(c,j) t2_aaaa(a,b,l,m) l2_1_aaaa(l,m,b,c)
        //             += +1.00 d_bb(i,k) r1_aa(c,m) t2_aaaa(a,b,l,j) l2_aaaa(m,l,b,c)
        //             += +1.00 d_bb(i,k) r1_aa(c,m) t2_1_aaaa(a,b,l,j) l2_1_aaaa(m,l,b,c)
        //             += +1.00 d_bb(i,k) r1_1_aa(c,m) t2_abab(a,b,j,l) l2_1_abab(m,l,c,b)
        //             += -0.50 d_bb(i,k) r1_1_aa(a,m) t2_abab(c,b,j,l) l2_1_abab(m,l,c,b)
        //             += -0.50 d_bb(i,k) r1_1_aa(a,m) t2_abab(b,c,j,l) l2_1_abab(m,l,b,c)
        //             += +0.50 d_bb(i,k) r1_aa(c,j) t2_1_aaaa(a,b,l,m) l2_1_aaaa(l,m,b,c)
        //             += -0.50 d_bb(i,k) r2_1_abab(a,c,m,l) t1_aa(b,j) l2_1_abab(m,l,b,c)
        //             += -0.50 d_bb(i,k) r2_1_abab(a,c,l,m) t1_aa(b,j) l2_1_abab(l,m,b,c)
        //             += -0.50 d_bb(i,k) r2_aaaa(a,c,m,l) t1_aa(b,j) l2_aaaa(m,l,b,c)
        //             += -0.50 d_bb(i,k) r2_abab(a,c,m,l) t1_aa(b,j) l2_abab(m,l,b,c)
        //             += -0.50 d_bb(i,k) r2_abab(a,c,l,m) t1_aa(b,j) l2_abab(l,m,b,c)
        //             += -0.50 d_bb(i,k) r2_1_aaaa(a,c,m,l) t1_aa(b,j) l2_1_aaaa(m,l,b,c)
        //             += +1.00 r0_1 t2_aaaa(a,b,l,j) t1_bb(c,i) l2_1_abab(l,k,b,c)
        //             += +1.00 r0_1 t2_abab(a,b,l,i) t1_aa(c,j) l2_1_abab(l,k,c,b)
        //             += -1.00 r1_aa(c,l) t1_bb(b,i) t1_1_aa(a,j) l2_1_abab(l,k,c,b)
        //             += +1.00 r1_bb(c,l) t1_bb(b,i) t1_1_aa(a,j) l2_1_bbbb(l,k,b,c)
        //             += -0.50 d_bb(i,k) r0 t1_aa(c,j) t2_1_aaaa(a,b,l,m) l2_1_aaaa(l,m,c,b)
        //             += +0.50 d_bb(i,k) r0 t2_aaaa(a,c,l,m) t1_1_aa(b,j) l2_1_aaaa(l,m,c,b)
        //             += -0.50 d_bb(i,k) r0 t2_abab(a,c,l,m) t1_1_aa(b,j) l2_1_abab(l,m,b,c)
        //             += -0.50 d_bb(i,k) r0 t2_abab(a,c,m,l) t1_1_aa(b,j) l2_1_abab(m,l,b,c)
        //             += -1.00 d_bb(i,k) r0 t1_aa(a,l) t1_aa(b,j) l1_aa(l,b)
        //             += +1.00 d_bb(i,k) r0 t2_abab(a,b,j,l) l1_bb(l,b)
        //             += -1.00 d_bb(i,k) r0 t1_aa(a,l) t1_1_aa(b,j) l1_1_aa(l,b)
        //             += -1.00 d_bb(i,k) r0 t2_aaaa(a,b,l,j) l1_aa(l,b)
        //             += -1.00 d_bb(i,k) r0 t1_aa(b,j) t1_1_aa(a,l) l1_1_aa(l,b)
        //             += +1.00 d_bb(i,k) r0 t2_1_abab(a,b,j,l) l1_1_bb(l,b)
        //             += -1.00 d_bb(i,k) r0 t2_1_aaaa(a,b,l,j) l1_1_aa(l,b)
        //             += +0.50 d_bb(i,k) r0 t2_aaaa(c,b,m,j) t1_1_aa(a,l) l2_1_aaaa(l,m,c,b)
        //             += -0.50 d_bb(i,k) r0 t2_abab(a,b,l,m) t1_aa(c,j) l2_abab(l,m,c,b)
        //             += -0.50 d_bb(i,k) r0 t2_abab(a,b,m,l) t1_aa(c,j) l2_abab(m,l,c,b)
        //             += -0.50 d_bb(i,k) r0 t2_aaaa(a,b,l,m) t1_aa(c,j) l2_aaaa(l,m,c,b)
        //             += -0.50 d_bb(i,k) r0 t1_aa(c,j) t2_1_abab(a,b,l,m) l2_1_abab(l,m,c,b)
        //             += -0.50 d_bb(i,k) r0 t1_aa(c,j) t2_1_abab(a,b,m,l) l2_1_abab(m,l,c,b)
        //             += -0.50 d_bb(i,k) r0 t2_abab(c,b,j,m) t1_1_aa(a,l) l2_1_abab(l,m,c,b)
        //             += -0.50 d_bb(i,k) r0 t2_abab(b,c,j,m) t1_1_aa(a,l) l2_1_abab(l,m,b,c)
        //             += +1.00 r1_bb(c,l) t1_aa(a,j) t1_bb(b,i) l2_bbbb(l,k,b,c)
        //             += -1.00 r1_aa(c,l) t1_aa(a,j) t1_bb(b,i) l2_abab(l,k,c,b)
        //             += -1.00 r1_1_aa(c,l) t1_aa(a,j) t1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += +1.00 r1_1_bb(c,l) t1_aa(a,j) t1_bb(b,i) l2_1_bbbb(l,k,b,c)
        //             += -1.00 d_bb(i,k) r1_bb(c,m) t1_aa(b,j) t1_1_aa(a,l) l2_1_abab(l,m,b,c)
        //             += +1.00 d_bb(i,k) r1_aa(c,m) t1_aa(b,j) t1_1_aa(a,l) l2_1_aaaa(m,l,b,c)
        //             += -1.00 d_bb(i,k) r1_bb(c,m) t1_aa(a,l) t1_aa(b,j) l2_abab(l,m,b,c)
        //             += -1.00 d_bb(i,k) r1_1_bb(c,m) t1_aa(a,l) t1_aa(b,j) l2_1_abab(l,m,b,c)
        //             += +1.00 d_bb(i,k) r1_aa(c,m) t1_aa(a,l) t1_aa(b,j) l2_aaaa(m,l,b,c)
        //             += +1.00 d_bb(i,k) r1_1_aa(c,m) t1_aa(a,l) t1_aa(b,j) l2_1_aaaa(m,l,b,c)
        //             += +1.00 d_bb(i,k) r1_aa(c,m) t1_aa(a,l) t1_1_aa(b,j) l2_1_aaaa(m,l,b,c)
        //             += -1.00 d_bb(i,k) r1_bb(c,m) t1_aa(a,l) t1_1_aa(b,j) l2_1_abab(l,m,b,c)
        //             += +1.00 r0 t1_aa(a,l) t1_aa(c,j) t1_1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += -0.50 r1_1_aa(a,j) t2_bbbb(c,b,l,i) l2_1_bbbb(l,k,c,b)
        //             += +1.00 d_bb(i,k) r1_aa(a,j) l0
        //             += -1.00 r1_1_aa(a,j) t1_bb(b,i) l1_1_bb(k,b)
        //             += +1.00 r0_1 t1_aa(a,l) t1_bb(b,i) t1_aa(c,j) l2_1_abab(l,k,c,b)
        //             += -1.00 r1_aa(a,j) t1_bb(b,i) l1_bb(k,b)
        //             += +1.00 r0_1 t2_abab(a,b,j,l) t1_bb(c,i) l2_1_bbbb(l,k,c,b)
        //             += -1.00 r1_aa(a,j) t1_1_bb(b,i) l1_1_bb(k,b)
        D_oovo_abab("L,R,i,j,a,k") -= tmps_["321_aabb_LLovoo"]("R,L,j,a,i,k");
        tmps_["321_aabb_LLovoo"].~TArrayD();
    }
} // hilbert
#endif