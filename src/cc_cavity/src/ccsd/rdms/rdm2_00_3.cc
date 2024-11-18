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
#include "cc_cavity/include/ccsd/eom_ee_rdm.h"
#include "misc/nonsym_davidson_solver.h"

namespace hilbert {


    void EOM_EE_RDM::rdm2_00_3() {

        auto &cc_wfn_ = eom_driver_->cc_wfn_;

        // extract 1-body amplitudes
        std::map<std::string, TA::TArrayD> t1 {
                {"aa", cc_wfn_->amplitudes_["t1_aa"]},
                {"bb", cc_wfn_->amplitudes_["t1_bb"]}
        };

        // extract 2-body amplitudes
        std::map<std::string, TA::TArrayD> t2 {
                {"aaaa", cc_wfn_->amplitudes_["t2_aaaa"]},
                {"abab", cc_wfn_->amplitudes_["t2_abab"]},
                {"bbbb", cc_wfn_->amplitudes_["t2_bbbb"]}
        };

        world_.gop.fence();

        TArrayMap &Id = cc_wfn_->Id_blks_;
        TArrayMap &evecs = eom_driver_->evec_blks_;


        /// unpack left/right operators

        TArrayMap r1, r2, l1, l2;

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


        // D_oovv_abab += -0.50 r1_aa(d,i) t2_abab(a,b,k,l) t1_bb(c,j) l2_abab(k,l,d,c)
        //             += -0.50 r1_aa(d,i) t2_abab(a,b,l,k) t1_bb(c,j) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r1["aa"]("R,d,i") * tmps_["96_aabb_Lvvvo"]("L,d,a,b,j");

        // D_oovv_abab += -0.50 r0 t2_abab(a,b,k,l) t1_aa(c,i) t1_bb(d,j) l2_abab(k,l,c,d)
        //             += -0.50 r0 t2_abab(a,b,l,k) t1_aa(c,i) t1_bb(d,j) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["aa"]("c,i") * tmps_["96_aabb_Lvvvo"]("L,c,a,b,j") * r0("R");

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["97_aa_Lvv"]("L,a,b")  = l2["aaaa"]("L,i,j,a,c") * t2["aaaa"]("b,c,i,j");

        // D_oovv_aaaa += -0.50 r0 t2_aaaa(a,c,k,l) t2_aaaa(b,d,i,j) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") -= 0.50 * t2["aaaa"]("b,d,i,j") * tmps_["97_aa_Lvv"]("L,d,a") * r0("R");

        // D_oovv_abab += +0.50 r0 t2_aaaa(a,c,k,l) t2_abab(d,b,i,j) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * t2["abab"]("d,b,i,j") * tmps_["97_aa_Lvv"]("L,d,a") * r0("R");

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["98_aa_Lvv"]("L,a,b")  = l2["abab"]("L,i,j,a,c") * t2["abab"]("b,c,i,j");

        // D_oovv_aaaa += +0.50 P(a,b) r2_aaaa(a,d,i,j) t2_abab(b,c,k,l) l2_abab(k,l,d,c)
        //             += +0.50 P(a,b) r2_aaaa(a,d,i,j) t2_abab(b,c,l,k) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = r2["aaaa"]("R,a,d,i,j") * tmps_["98_aa_Lvv"]("L,d,b");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r2_abab(d,b,i,j) t2_abab(a,c,k,l) l2_abab(k,l,d,c)
        //             += +0.50 r2_abab(d,b,i,j) t2_abab(a,c,l,k) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += r2["abab"]("R,d,b,i,j") * tmps_["98_aa_Lvv"]("L,d,a");

        // D_oovv_abab += +0.50 r0 t2_abab(a,c,k,l) t2_abab(d,b,i,j) l2_abab(k,l,d,c)
        //             += +0.50 r0 t2_abab(a,c,l,k) t2_abab(d,b,i,j) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("d,b,i,j") * tmps_["98_aa_Lvv"]("L,d,a") * r0("R");

        // flops: o0v2L2  = o2v3L2 o2v3L2 o0v2L2
        //  mems: o0v2L2  = o0v2L2 o0v2L2 o0v2L2
        tmps_["101_aa_LLvv"]("L,R,a,b")  = 0.50 * l2["aaaa"]("L,i,k,a,d") * r2["aaaa"]("R,b,d,i,k");
        tmps_["101_aa_LLvv"]("L,R,a,b") += l2["abab"]("L,i,j,a,c") * r2["abab"]("R,b,c,i,j");

        // D_oovv_aaaa += -0.50 P(a,b) r2_abab(a,d,l,k) t2_aaaa(b,c,i,j) l2_abab(l,k,c,d)
        //             += -0.50 P(a,b) r2_abab(a,d,k,l) t2_aaaa(b,c,i,j) l2_abab(k,l,c,d)
        //             += -0.50 P(a,b) r2_aaaa(a,d,l,k) t2_aaaa(b,c,i,j) l2_aaaa(l,k,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["aaaa"]("b,c,i,j") * tmps_["101_aa_LLvv"]("L,R,c,a");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r2_abab(a,d,l,k) t2_abab(c,b,i,j) l2_abab(l,k,c,d)
        //             += +0.50 r2_abab(a,d,k,l) t2_abab(c,b,i,j) l2_abab(k,l,c,d)
        //             += +0.50 r2_aaaa(a,d,l,k) t2_abab(c,b,i,j) l2_aaaa(l,k,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("c,b,i,j") * tmps_["101_aa_LLvv"]("L,R,c,a");

        // flops: o0v2L2  = o1v1L2 o1v2L2
        //  mems: o0v2L2  = o1v1L2 o0v2L2
        tmps_["102_aa_LLvv"]("L,R,a,b")  = (tmps_["56_aa_LLov"]("L,R,i,a") + -1.00 * tmps_["55_aa_LLov"]("L,R,i,a")) * t1["aa"]("b,i");

        // flops: o0v2L1  = o1v2L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["100_aa_Lvv"]("L,a,b")  = l1["aa"]("L,i,a") * t1["aa"]("b,i");

        // flops: o0v2L2  = o1v2L2
        //  mems: o0v2L2  = o0v2L2
        tmps_["99_aa_LLvv"]("L,R,a,b")  = l1["aa"]("L,i,a") * r1["aa"]("R,b,i");

        // flops: o1v3L2  = o1v3L2 o1v3L2 o2v3L2 o1v4L2 o1v3L2 o2v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o2v3L2 o1v3L2 o3v2L1 o3v2L2 o2v3L2 o1v3L2 o3v2L2 o3v2L2 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L1 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        //  mems: o1v3L2  = o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o3v1L1 o2v2L2 o1v3L2 o1v3L2 o3v1L2 o2v2L2 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L1 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i")  = -0.50 * tmps_["97_aa_Lvv"]("L,a,b") * r1["bb"]("R,c,i");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") -= 0.50 * tmps_["97_aa_Lvv"]("L,a,b") * tmps_["90_bb_Lvo"]("R,c,i");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") -= tmps_["56_aa_LLov"]("L,R,j,a") * t2["abab"]("b,c,j,i");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") += tmps_["94_abab_Lvvvv"]("L,a,e,b,c") * r1["bb"]("R,e,i");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") += tmps_["55_aa_LLov"]("L,R,j,a") * t2["abab"]("b,c,j,i");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") -= l2["aaaa"]("L,l,j,a,d") * t2["abab"]("d,c,j,i") * r1["aa"]("R,b,l");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") -= l1["aa"]("L,j,a") * t2["abab"]("b,c,j,i") * r0("R");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") += l2["abab"]("L,j,m,a,e") * t2["abab"]("b,e,j,i") * r1["bb"]("R,c,m");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") += l2["abab"]("L,l,k,a,e") * t1["bb"]("e,i") * t1["bb"]("c,k") * r1["aa"]("R,b,l");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") += l2["abab"]("L,l,k,a,e") * t2["bbbb"]("c,e,k,i") * r1["aa"]("R,b,l");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") -= l1["aa"]("L,j,a") * r2["abab"]("R,b,c,j,i");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") += l2["abab"]("L,j,m,a,e") * t1["bb"]("e,i") * r1["bb"]("R,c,m") * t1["aa"]("b,j");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") += l2["abab"]("L,l,k,a,e") * r1["bb"]("R,e,i") * t1["bb"]("c,k") * t1["aa"]("b,l");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") += l2["abab"]("L,l,k,a,e") * t1["bb"]("e,i") * t1["bb"]("c,k") * t1["aa"]("b,l") * r0("R");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") += t1["bb"]("c,m") * tmps_["12_abba_Lvoov"]("L,b,i,m,a") * r0("R");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") += t1["aa"]("b,l") * tmps_["43_bbaa_Lvoov"]("L,c,i,l,a") * r0("R");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") -= t1["bb"]("c,i") * tmps_["101_aa_LLvv"]("L,R,a,b");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") -= tmps_["100_aa_Lvv"]("L,a,b") * tmps_["90_bb_Lvo"]("R,c,i");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") += t1["aa"]("b,l") * tmps_["95_bbaa_Lvoov"]("L,c,i,l,a") * r0("R");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") -= tmps_["98_aa_Lvv"]("L,a,b") * tmps_["90_bb_Lvo"]("R,c,i");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") += tmps_["96_aabb_Lvvvo"]("L,a,b,c,i") * r0("R");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") -= tmps_["99_aa_LLvv"]("L,R,a,b") * t1["bb"]("c,i");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") -= tmps_["98_aa_Lvv"]("L,a,b") * r1["bb"]("R,c,i");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") += tmps_["38_baab_LLovvo"]("L,R,k,a,b,i") * t1["bb"]("c,k");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") -= tmps_["100_aa_Lvv"]("L,a,b") * r1["bb"]("R,c,i");
        tmps_["103_aabb_LLvvvo"]("R,L,a,b,c,i") -= t1["bb"]("c,i") * tmps_["102_aa_LLvv"]("L,R,a,b");
        tmps_["96_aabb_Lvvvo"].~TArrayD();
        tmps_["95_bbaa_Lvoov"].~TArrayD();
        tmps_["94_abab_Lvvvv"].~TArrayD();
        tmps_["38_baab_LLovvo"].~TArrayD();

        // D_ovvv_abab += +0.50 r0 t2_abab(c,b,j,k) t1_bb(d,i) l2_abab(j,k,a,d)
        //             += +0.50 r0 t2_abab(c,b,k,j) t1_bb(d,i) l2_abab(k,j,a,d)
        //             += +1.00 r0 t1_bb(b,k) t2_abab(c,d,j,i) l2_abab(j,k,a,d)
        //             += -1.00 r1_bb(d,k) t2_abab(c,b,j,i) l2_abab(j,k,a,d)
        //             += +0.50 r1_bb(d,i) t2_abab(c,b,j,k) l2_abab(j,k,a,d)
        //             += +0.50 r1_bb(d,i) t2_abab(c,b,k,j) l2_abab(k,j,a,d)
        //             += +1.00 r1_aa(d,k) t2_abab(c,b,j,i) l2_aaaa(k,j,a,d)
        //             += -1.00 r1_aa(c,k) t2_abab(d,b,j,i) l2_aaaa(k,j,a,d)
        //             += -1.00 r0 t2_abab(c,b,j,i) l1_aa(j,a)
        //             += +1.00 r1_bb(b,k) t2_abab(c,d,j,i) l2_abab(j,k,a,d)
        //             += +1.00 r1_aa(c,k) t1_bb(b,j) t1_bb(d,i) l2_abab(k,j,a,d)
        //             += +1.00 r1_aa(c,k) t2_bbbb(b,d,j,i) l2_abab(k,j,a,d)
        //             += -1.00 r2_abab(c,b,j,i) l1_aa(j,a)
        //             += +1.00 r1_bb(b,k) t1_aa(c,j) t1_bb(d,i) l2_abab(j,k,a,d)
        //             += +1.00 r1_bb(d,i) t1_bb(b,j) t1_aa(c,k) l2_abab(k,j,a,d)
        //             += +1.00 r0 t1_bb(b,j) t1_aa(c,k) t1_bb(d,i) l2_abab(k,j,a,d)
        //             += +1.00 r0 t2_bbbb(b,d,j,i) t1_aa(c,k) l2_abab(k,j,a,d)
        //             += -0.50 r0 t1_bb(b,i) t2_aaaa(c,d,j,k) l2_aaaa(j,k,a,d)
        //             += -0.50 r2_abab(c,d,k,j) t1_bb(b,i) l2_abab(k,j,a,d)
        //             += -0.50 r2_abab(c,d,j,k) t1_bb(b,i) l2_abab(j,k,a,d)
        //             += -0.50 r2_aaaa(c,d,k,j) t1_bb(b,i) l2_aaaa(k,j,a,d)
        //             += -1.00 r0 t1_bb(b,i) t1_aa(c,j) l1_aa(j,a)
        //             += +1.00 r0 t2_abab(d,b,j,i) t1_aa(c,k) l2_aaaa(j,k,a,d)
        //             += -0.50 r0 t1_bb(b,i) t2_abab(c,d,j,k) l2_abab(j,k,a,d)
        //             += -0.50 r0 t1_bb(b,i) t2_abab(c,d,k,j) l2_abab(k,j,a,d)
        //             += -1.00 r1_aa(c,j) t1_bb(b,i) l1_aa(j,a)
        //             += -0.50 r1_bb(b,i) t2_abab(c,d,j,k) l2_abab(j,k,a,d)
        //             += -0.50 r1_bb(b,i) t2_abab(c,d,k,j) l2_abab(k,j,a,d)
        //             += +1.00 r2_abab(c,d,k,i) t1_bb(b,j) l2_abab(k,j,a,d)
        //             += -1.00 r1_bb(b,i) t1_aa(c,j) l1_aa(j,a)
        //             += -0.50 r1_bb(b,i) t2_aaaa(c,d,j,k) l2_aaaa(j,k,a,d)
        //             += -1.00 r1_bb(d,k) t1_bb(b,i) t1_aa(c,j) l2_abab(j,k,a,d)
        //             += +1.00 r1_aa(d,k) t1_bb(b,i) t1_aa(c,j) l2_aaaa(k,j,a,d)
        D_ovvv_abab("L,R,i,a,b,c") += tmps_["103_aabb_LLvvvo"]("R,L,a,c,b,i");

        // D_vovv_abab += -0.50 r0 t2_abab(c,b,j,k) t1_bb(d,i) l2_abab(j,k,a,d)
        //             += -0.50 r0 t2_abab(c,b,k,j) t1_bb(d,i) l2_abab(k,j,a,d)
        //             += -1.00 r0 t1_bb(b,k) t2_abab(c,d,j,i) l2_abab(j,k,a,d)
        //             += +1.00 r1_bb(d,k) t2_abab(c,b,j,i) l2_abab(j,k,a,d)
        //             += -0.50 r1_bb(d,i) t2_abab(c,b,j,k) l2_abab(j,k,a,d)
        //             += -0.50 r1_bb(d,i) t2_abab(c,b,k,j) l2_abab(k,j,a,d)
        //             += -1.00 r1_aa(d,k) t2_abab(c,b,j,i) l2_aaaa(k,j,a,d)
        //             += +1.00 r1_aa(c,k) t2_abab(d,b,j,i) l2_aaaa(k,j,a,d)
        //             += +1.00 r0 t2_abab(c,b,j,i) l1_aa(j,a)
        //             += -1.00 r1_bb(b,k) t2_abab(c,d,j,i) l2_abab(j,k,a,d)
        //             += -1.00 r1_aa(c,k) t1_bb(b,j) t1_bb(d,i) l2_abab(k,j,a,d)
        //             += -1.00 r1_aa(c,k) t2_bbbb(b,d,j,i) l2_abab(k,j,a,d)
        //             += +1.00 r2_abab(c,b,j,i) l1_aa(j,a)
        //             += -1.00 r1_bb(b,k) t1_aa(c,j) t1_bb(d,i) l2_abab(j,k,a,d)
        //             += -1.00 r1_bb(d,i) t1_bb(b,j) t1_aa(c,k) l2_abab(k,j,a,d)
        //             += -1.00 r0 t1_bb(b,j) t1_aa(c,k) t1_bb(d,i) l2_abab(k,j,a,d)
        //             += -1.00 r0 t2_bbbb(b,d,j,i) t1_aa(c,k) l2_abab(k,j,a,d)
        //             += +0.50 r0 t1_bb(b,i) t2_aaaa(c,d,j,k) l2_aaaa(j,k,a,d)
        //             += +0.50 r2_abab(c,d,k,j) t1_bb(b,i) l2_abab(k,j,a,d)
        //             += +0.50 r2_abab(c,d,j,k) t1_bb(b,i) l2_abab(j,k,a,d)
        //             += +0.50 r2_aaaa(c,d,k,j) t1_bb(b,i) l2_aaaa(k,j,a,d)
        //             += +1.00 r0 t1_bb(b,i) t1_aa(c,j) l1_aa(j,a)
        //             += -1.00 r0 t2_abab(d,b,j,i) t1_aa(c,k) l2_aaaa(j,k,a,d)
        //             += +0.50 r0 t1_bb(b,i) t2_abab(c,d,j,k) l2_abab(j,k,a,d)
        //             += +0.50 r0 t1_bb(b,i) t2_abab(c,d,k,j) l2_abab(k,j,a,d)
        //             += +1.00 r1_aa(c,j) t1_bb(b,i) l1_aa(j,a)
        //             += +0.50 r1_bb(b,i) t2_abab(c,d,j,k) l2_abab(j,k,a,d)
        //             += +0.50 r1_bb(b,i) t2_abab(c,d,k,j) l2_abab(k,j,a,d)
        //             += -1.00 r2_abab(c,d,k,i) t1_bb(b,j) l2_abab(k,j,a,d)
        //             += +1.00 r1_bb(b,i) t1_aa(c,j) l1_aa(j,a)
        //             += +0.50 r1_bb(b,i) t2_aaaa(c,d,j,k) l2_aaaa(j,k,a,d)
        //             += +1.00 r1_bb(d,k) t1_bb(b,i) t1_aa(c,j) l2_abab(j,k,a,d)
        //             += -1.00 r1_aa(d,k) t1_bb(b,i) t1_aa(c,j) l2_aaaa(k,j,a,d)
        D_vovv_abab("L,R,a,i,b,c") -= tmps_["103_aabb_LLvvvo"]("R,L,a,c,b,i");
        tmps_["103_aabb_LLvvvo"].~TArrayD();

        // flops: o2v2L2  = o2v2L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o3v2L1 o3v2L2 o2v2L2 o2v2L2 o3v2L1 o3v2L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2
        //  mems: o2v2L2  = o2v2L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o3v1L1 o2v2L2 o2v2L2 o2v2L2 o3v1L1 o2v2L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2
        tmps_["104_bbbb_LLvoov"]("L,R,a,i,j,b")  = (tmps_["24_bbbb_Lvoov"]("L,a,i,j,b") + tmps_["86_bbbb_Lvoov"]("L,a,i,j,b")) * r0("R");
        tmps_["104_bbbb_LLvoov"]("L,R,a,i,j,b") -= Id["bb_oo"]("i,j") * tmps_["88_bb_LLvv"]("L,R,b,a");
        tmps_["104_bbbb_LLvoov"]("L,R,a,i,j,b") -= tmps_["30_bb_Lvv"]("L,b,a") * tmps_["52_bb_Loo"]("R,i,j");
        tmps_["104_bbbb_LLvoov"]("L,R,a,i,j,b") += t1["bb"]("d,i") * l2["bbbb"]("L,k,j,b,d") * r1["bb"]("R,a,k");
        tmps_["104_bbbb_LLvoov"]("L,R,a,i,j,b") += l1["bb"]("L,j,b") * r1["bb"]("R,a,i");
        tmps_["104_bbbb_LLvoov"]("L,R,a,i,j,b") += t1["bb"]("d,i") * l2["bbbb"]("L,k,j,b,d") * t1["bb"]("a,k") * r0("R");
        tmps_["104_bbbb_LLvoov"]("L,R,a,i,j,b") += t1["bb"]("a,i") * tmps_["39_bb_LLov"]("L,R,j,b");
        tmps_["104_bbbb_LLvoov"]("L,R,a,i,j,b") -= t1["bb"]("a,i") * tmps_["40_bb_LLov"]("L,R,j,b");
        tmps_["104_bbbb_LLvoov"]("L,R,a,i,j,b") -= tmps_["91_bb_LLvv"]("L,R,b,a") * Id["bb_oo"]("i,j");
        tmps_["104_bbbb_LLvoov"]("L,R,a,i,j,b") += l1["bb"]("L,j,b") * tmps_["90_bb_Lvo"]("R,a,i");
        tmps_["104_bbbb_LLvoov"]("L,R,a,i,j,b") -= tmps_["89_bb_Lvv"]("L,b,a") * tmps_["52_bb_Loo"]("R,i,j");
        tmps_["104_bbbb_LLvoov"]("L,R,a,i,j,b") -= 0.50 * tmps_["87_bb_Lvv"]("L,b,a") * tmps_["52_bb_Loo"]("R,i,j");
        tmps_["104_bbbb_LLvoov"]("L,R,a,i,j,b") -= Id["bb_oo"]("i,j") * tmps_["92_bb_LLvv"]("L,R,b,a");
        tmps_["92_bb_LLvv"].~TArrayD();
        tmps_["89_bb_Lvv"].~TArrayD();
        tmps_["88_bb_LLvv"].~TArrayD();
        tmps_["40_bb_LLov"].~TArrayD();

        // D_ovov_bbbb += +1.00 r0 t2_abab(c,b,k,i) l2_abab(k,j,c,a)
        //             += +1.00 r0 t2_bbbb(b,c,k,i) l2_bbbb(k,j,a,c)
        //             += -1.00 d_bb(i,j) r1_bb(b,k) l1_bb(k,a)
        //             += -0.50 d_bb(i,j) r0 t2_abab(c,b,k,l) l2_abab(k,l,c,a)
        //             += -0.50 d_bb(i,j) r0 t2_abab(c,b,l,k) l2_abab(l,k,c,a)
        //             += +1.00 r1_bb(b,k) t1_bb(c,i) l2_bbbb(k,j,a,c)
        //             += +1.00 r1_bb(b,i) l1_bb(j,a)
        //             += +1.00 r0 t1_bb(b,k) t1_bb(c,i) l2_bbbb(k,j,a,c)
        //             += +1.00 r1_aa(c,k) t1_bb(b,i) l2_abab(k,j,c,a)
        //             += -1.00 r1_bb(c,k) t1_bb(b,i) l2_bbbb(k,j,a,c)
        //             += -0.50 d_bb(i,j) r2_abab(c,b,l,k) l2_abab(l,k,c,a)
        //             += -0.50 d_bb(i,j) r2_abab(c,b,k,l) l2_abab(k,l,c,a)
        //             += -0.50 d_bb(i,j) r2_bbbb(b,c,l,k) l2_bbbb(l,k,a,c)
        //             += +1.00 r0 t1_bb(b,i) l1_bb(j,a)
        //             += -1.00 d_bb(i,j) r0 t1_bb(b,k) l1_bb(k,a)
        //             += -0.50 d_bb(i,j) r0 t2_bbbb(b,c,k,l) l2_bbbb(k,l,a,c)
        //             += -1.00 d_bb(i,j) r1_aa(c,l) t1_bb(b,k) l2_abab(l,k,c,a)
        //             += +1.00 d_bb(i,j) r1_bb(c,l) t1_bb(b,k) l2_bbbb(l,k,a,c)
        D_ovov_bbbb("L,R,i,a,j,b") += tmps_["104_bbbb_LLvoov"]("L,R,b,i,j,a");

        // D_ovvo_bbbb += -1.00 r0 t2_abab(c,b,k,i) l2_abab(k,j,c,a)
        //             += -1.00 r0 t2_bbbb(b,c,k,i) l2_bbbb(k,j,a,c)
        //             += +1.00 d_bb(i,j) r1_bb(b,k) l1_bb(k,a)
        //             += +0.50 d_bb(i,j) r0 t2_abab(c,b,k,l) l2_abab(k,l,c,a)
        //             += +0.50 d_bb(i,j) r0 t2_abab(c,b,l,k) l2_abab(l,k,c,a)
        //             += -1.00 r1_bb(b,k) t1_bb(c,i) l2_bbbb(k,j,a,c)
        //             += -1.00 r1_bb(b,i) l1_bb(j,a)
        //             += -1.00 r0 t1_bb(b,k) t1_bb(c,i) l2_bbbb(k,j,a,c)
        //             += -1.00 r1_aa(c,k) t1_bb(b,i) l2_abab(k,j,c,a)
        //             += +1.00 r1_bb(c,k) t1_bb(b,i) l2_bbbb(k,j,a,c)
        //             += +0.50 d_bb(i,j) r2_abab(c,b,l,k) l2_abab(l,k,c,a)
        //             += +0.50 d_bb(i,j) r2_abab(c,b,k,l) l2_abab(k,l,c,a)
        //             += +0.50 d_bb(i,j) r2_bbbb(b,c,l,k) l2_bbbb(l,k,a,c)
        //             += -1.00 r0 t1_bb(b,i) l1_bb(j,a)
        //             += +1.00 d_bb(i,j) r0 t1_bb(b,k) l1_bb(k,a)
        //             += +0.50 d_bb(i,j) r0 t2_bbbb(b,c,k,l) l2_bbbb(k,l,a,c)
        //             += +1.00 d_bb(i,j) r1_aa(c,l) t1_bb(b,k) l2_abab(l,k,c,a)
        //             += -1.00 d_bb(i,j) r1_bb(c,l) t1_bb(b,k) l2_bbbb(l,k,a,c)
        D_ovvo_bbbb("L,R,i,a,b,j") -= tmps_["104_bbbb_LLvoov"]("L,R,b,i,j,a");
        tmps_["104_bbbb_LLvoov"].~TArrayD();

        // flops: o1v1L1  = o1v1L1
        //  mems: o1v1L1  = o1v1L1
        tmps_["105_bb_Lvo"]("L,a,i")  = l0("L") * t1["bb"]("a,i");

        // D_oovv_abab += -1.00 r1_aa(a,i) t1_bb(b,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r1["aa"]("R,a,i") * tmps_["105_bb_Lvo"]("L,b,j");

        // D_oovv_bbbb += -1.00 P(i,j) P(a,b) r1_bb(a,i) t1_bb(b,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = r1["bb"]("R,a,i") * tmps_["105_bb_Lvo"]("L,b,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += -1.00 P(i,j) r0 t1_bb(a,i) t1_bb(b,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["90_bb_Lvo"]("R,a,i") * tmps_["105_bb_Lvo"]("L,b,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v1L1  = o1v1L1
        //  mems: o1v1L1  = o1v1L1
        tmps_["106_aa_Lvo"]("L,a,i")  = l0("L") * t1["aa"]("a,i");

        // D_oovv_aaaa += -1.00 P(i,j) P(a,b) r1_aa(a,i) t1_aa(b,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = r1["aa"]("R,a,i") * tmps_["106_aa_Lvo"]("L,b,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r1_bb(b,j) t1_aa(a,i) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r1["bb"]("R,b,j") * tmps_["106_aa_Lvo"]("L,a,i");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["107_aaaa_Lvoov"]("L,a,i,j,b")  = l2["aaaa"]("L,k,j,b,c") * t2["aaaa"]("a,c,k,i");

        // D_oovv_aaaa += +1.00 P(i,j) r0 t2_aaaa(a,c,k,i) t2_aaaa(b,d,l,j) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["aaaa"]("b,d,l,j") * tmps_["107_aaaa_Lvoov"]("L,a,i,l,d") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r0 t2_aaaa(a,c,k,i) t2_abab(d,b,l,j) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("d,b,l,j") * tmps_["107_aaaa_Lvoov"]("L,a,i,l,d") * r0("R");

        // flops: o1v1L1  = o1v1L1
        //  mems: o1v1L1  = o1v1L1
        tmps_["109_aa_Lvo"]("R,a,i")  = r0("R") * t1["aa"]("a,i");

        // D_oovv_abab += -1.00 r0 t1_aa(a,i) t1_bb(b,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["109_aa_Lvo"]("R,a,i") * tmps_["105_bb_Lvo"]("L,b,j");
        tmps_["105_bb_Lvo"].~TArrayD();

        // D_oovv_aaaa += -1.00 P(i,j) r0 t1_aa(a,i) t1_aa(b,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["109_aa_Lvo"]("R,a,i") * tmps_["106_aa_Lvo"]("L,b,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();
        tmps_["106_aa_Lvo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["108_aaaa_Lvoov"]("L,a,i,j,b")  = l2["abab"]("L,j,k,b,c") * t2["abab"]("a,c,i,k");

        // flops: o1v3L2  = o1v3L2 o2v3L1 o1v3L2 o1v3L2 o3v3L1 o2v3L2 o3v3L1 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        //  mems: o1v3L2  = o1v3L2 o1v3L1 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        tmps_["110_aaaa_LLvvvo"]("L,R,a,b,c,i")  = 0.50 * tmps_["97_aa_Lvv"]("L,a,b") * tmps_["109_aa_Lvo"]("R,c,i");
        tmps_["110_aaaa_LLvvvo"]("L,R,a,b,c,i") -= tmps_["107_aaaa_Lvoov"]("L,b,i,l,a") * t1["aa"]("c,l") * r0("R");
        tmps_["110_aaaa_LLvvvo"]("L,R,a,b,c,i") += l2["abab"]("L,l,m,a,d") * t2["abab"]("c,d,i,m") * r1["aa"]("R,b,l");
        tmps_["110_aaaa_LLvvvo"]("L,R,a,b,c,i") -= l2["aaaa"]("L,l,j,a,e") * t2["aaaa"]("c,e,j,i") * r1["aa"]("R,b,l");
        tmps_["110_aaaa_LLvvvo"]("L,R,a,b,c,i") -= l2["aaaa"]("L,l,j,a,e") * t1["aa"]("e,i") * t1["aa"]("c,j") * r1["aa"]("R,b,l");
        tmps_["110_aaaa_LLvvvo"]("L,R,a,b,c,i") -= tmps_["108_aaaa_Lvoov"]("L,b,i,l,a") * t1["aa"]("c,l") * r0("R");
        tmps_["110_aaaa_LLvvvo"]("L,R,a,b,c,i") += tmps_["100_aa_Lvv"]("L,a,b") * tmps_["109_aa_Lvo"]("R,c,i");
        tmps_["110_aaaa_LLvvvo"]("L,R,a,b,c,i") += tmps_["98_aa_Lvv"]("L,a,b") * tmps_["109_aa_Lvo"]("R,c,i");
        tmps_["110_aaaa_LLvvvo"]("L,R,a,b,c,i") += tmps_["101_aa_LLvv"]("L,R,a,b") * t1["aa"]("c,i");
        tmps_["110_aaaa_LLvvvo"]("L,R,a,b,c,i") -= t1["aa"]("c,j") * tmps_["35_aaaa_LLovvo"]("L,R,j,a,b,i");
        tmps_["110_aaaa_LLvvvo"]("L,R,a,b,c,i") -= t1["aa"]("c,j") * tmps_["34_aaaa_LLovvo"]("L,R,j,a,b,i");
        tmps_["110_aaaa_LLvvvo"]("L,R,a,b,c,i") -= 0.50 * r1["aa"]("R,b,i") * tmps_["97_aa_Lvv"]("L,a,c");
        tmps_["110_aaaa_LLvvvo"]("L,R,a,b,c,i") -= r1["aa"]("R,b,i") * tmps_["100_aa_Lvv"]("L,a,c");
        tmps_["110_aaaa_LLvvvo"]("L,R,a,b,c,i") += t1["aa"]("c,i") * tmps_["99_aa_LLvv"]("L,R,a,b");
        tmps_["110_aaaa_LLvvvo"]("L,R,a,b,c,i") -= r1["aa"]("R,b,i") * tmps_["98_aa_Lvv"]("L,a,c");
        tmps_["110_aaaa_LLvvvo"]("L,R,a,b,c,i") += t1["aa"]("c,i") * tmps_["102_aa_LLvv"]("L,R,a,b");
        tmps_["35_aaaa_LLovvo"].~TArrayD();
        tmps_["34_aaaa_LLovvo"].~TArrayD();

        // D_ovvv_aaaa += +1.00 P(b,c) r1_bb(d,k) t1_aa(b,j) t1_aa(c,i) l2_abab(j,k,a,d)
        //             += -1.00 P(b,c) r1_aa(d,k) t1_aa(b,j) t1_aa(c,i) l2_aaaa(k,j,a,d)
        //             += +0.50 P(b,c) r0 t2_aaaa(b,d,j,k) t1_aa(c,i) l2_aaaa(j,k,a,d)
        //             += -1.00 P(b,c) r0 t2_aaaa(b,d,j,i) t1_aa(c,k) l2_aaaa(j,k,a,d)
        //             += +1.00 P(b,c) r1_aa(b,k) t2_abab(c,d,i,j) l2_abab(k,j,a,d)
        //             += -1.00 P(b,c) r1_aa(b,k) t2_aaaa(c,d,j,i) l2_aaaa(k,j,a,d)
        //             += -1.00 P(b,c) r1_aa(b,k) t1_aa(c,j) t1_aa(d,i) l2_aaaa(k,j,a,d)
        //             += -1.00 P(b,c) r0 t2_abab(b,d,i,j) t1_aa(c,k) l2_abab(k,j,a,d)
        //             += +1.00 P(b,c) r0 t1_aa(b,j) t1_aa(c,i) l1_aa(j,a)
        //             += +0.50 P(b,c) r0 t2_abab(b,d,j,k) t1_aa(c,i) l2_abab(j,k,a,d)
        //             += +0.50 P(b,c) r0 t2_abab(b,d,k,j) t1_aa(c,i) l2_abab(k,j,a,d)
        //             += +0.50 P(b,c) r2_abab(b,d,k,j) t1_aa(c,i) l2_abab(k,j,a,d)
        //             += +0.50 P(b,c) r2_abab(b,d,j,k) t1_aa(c,i) l2_abab(j,k,a,d)
        //             += +0.50 P(b,c) r2_aaaa(b,d,k,j) t1_aa(c,i) l2_aaaa(k,j,a,d)
        //             += -1.00 P(b,c) r2_abab(b,d,i,k) t1_aa(c,j) l2_abab(j,k,a,d)
        //             += -1.00 P(b,c) r2_aaaa(b,d,k,i) t1_aa(c,j) l2_aaaa(k,j,a,d)
        //             += -0.50 P(b,c) r1_aa(b,i) t2_aaaa(c,d,j,k) l2_aaaa(j,k,a,d)
        //             += -1.00 P(b,c) r1_aa(b,i) t1_aa(c,j) l1_aa(j,a)
        //             += +1.00 P(b,c) r1_aa(b,j) t1_aa(c,i) l1_aa(j,a)
        //             += -0.50 P(b,c) r1_aa(b,i) t2_abab(c,d,j,k) l2_abab(j,k,a,d)
        //             += -0.50 P(b,c) r1_aa(b,i) t2_abab(c,d,k,j) l2_abab(k,j,a,d)
        D_ovvv_aaaa("L,R,i,a,b,c") += tmps_["110_aaaa_LLvvvo"]("L,R,a,b,c,i");
        D_ovvv_aaaa("L,R,i,a,b,c") -= tmps_["110_aaaa_LLvvvo"]("L,R,a,c,b,i");

        // D_vovv_aaaa += -1.00 P(b,c) r1_bb(d,k) t1_aa(b,j) t1_aa(c,i) l2_abab(j,k,a,d)
        //             += +1.00 P(b,c) r1_aa(d,k) t1_aa(b,j) t1_aa(c,i) l2_aaaa(k,j,a,d)
        //             += -0.50 P(b,c) r0 t2_aaaa(b,d,j,k) t1_aa(c,i) l2_aaaa(j,k,a,d)
        //             += +1.00 P(b,c) r0 t2_aaaa(b,d,j,i) t1_aa(c,k) l2_aaaa(j,k,a,d)
        //             += -1.00 P(b,c) r1_aa(b,k) t2_abab(c,d,i,j) l2_abab(k,j,a,d)
        //             += +1.00 P(b,c) r1_aa(b,k) t2_aaaa(c,d,j,i) l2_aaaa(k,j,a,d)
        //             += +1.00 P(b,c) r1_aa(b,k) t1_aa(c,j) t1_aa(d,i) l2_aaaa(k,j,a,d)
        //             += +1.00 P(b,c) r0 t2_abab(b,d,i,j) t1_aa(c,k) l2_abab(k,j,a,d)
        //             += -1.00 P(b,c) r0 t1_aa(b,j) t1_aa(c,i) l1_aa(j,a)
        //             += -0.50 P(b,c) r0 t2_abab(b,d,j,k) t1_aa(c,i) l2_abab(j,k,a,d)
        //             += -0.50 P(b,c) r0 t2_abab(b,d,k,j) t1_aa(c,i) l2_abab(k,j,a,d)
        //             += -0.50 P(b,c) r2_abab(b,d,k,j) t1_aa(c,i) l2_abab(k,j,a,d)
        //             += -0.50 P(b,c) r2_abab(b,d,j,k) t1_aa(c,i) l2_abab(j,k,a,d)
        //             += -0.50 P(b,c) r2_aaaa(b,d,k,j) t1_aa(c,i) l2_aaaa(k,j,a,d)
        //             += +1.00 P(b,c) r2_abab(b,d,i,k) t1_aa(c,j) l2_abab(j,k,a,d)
        //             += +1.00 P(b,c) r2_aaaa(b,d,k,i) t1_aa(c,j) l2_aaaa(k,j,a,d)
        //             += +0.50 P(b,c) r1_aa(b,i) t2_aaaa(c,d,j,k) l2_aaaa(j,k,a,d)
        //             += +1.00 P(b,c) r1_aa(b,i) t1_aa(c,j) l1_aa(j,a)
        //             += -1.00 P(b,c) r1_aa(b,j) t1_aa(c,i) l1_aa(j,a)
        //             += +0.50 P(b,c) r1_aa(b,i) t2_abab(c,d,j,k) l2_abab(j,k,a,d)
        //             += +0.50 P(b,c) r1_aa(b,i) t2_abab(c,d,k,j) l2_abab(k,j,a,d)
        D_vovv_aaaa("L,R,a,i,b,c") -= tmps_["110_aaaa_LLvvvo"]("L,R,a,b,c,i");
        D_vovv_aaaa("L,R,a,i,b,c") += tmps_["110_aaaa_LLvvvo"]("L,R,a,c,b,i");
        tmps_["110_aaaa_LLvvvo"].~TArrayD();

        // flops: o2v2L2  = o0v2L1 o2v2L2 o3v2L1 o1v1L1 o1v1L1 o3v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2
        //  mems: o2v2L2  = o0v2L1 o2v2L2 o3v1L1 o1v1L1 o1v1L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2
        tmps_["111_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * (tmps_["97_aa_Lvv"]("L,a,b") + 2.00 * tmps_["98_aa_Lvv"]("L,a,b")) * tmps_["61_aa_Loo"]("R,i,j");
        tmps_["111_aaaa_LLvvoo"]("L,R,a,b,i,j") -= l2["aaaa"]("L,k,j,a,d") * t1["aa"]("d,i") * (r0("R") * t1["aa"]("b,k") + r1["aa"]("R,b,k"));
        tmps_["111_aaaa_LLvvoo"]("L,R,a,b,i,j") -= l1["aa"]("L,j,a") * r1["aa"]("R,b,i");
        tmps_["111_aaaa_LLvvoo"]("L,R,a,b,i,j") += tmps_["101_aa_LLvv"]("L,R,a,b") * Id["aa_oo"]("i,j");
        tmps_["111_aaaa_LLvvoo"]("L,R,a,b,i,j") -= l1["aa"]("L,j,a") * tmps_["109_aa_Lvo"]("R,b,i");
        tmps_["111_aaaa_LLvvoo"]("L,R,a,b,i,j") += t1["aa"]("b,i") * tmps_["55_aa_LLov"]("L,R,j,a");
        tmps_["111_aaaa_LLvvoo"]("L,R,a,b,i,j") -= t1["aa"]("b,i") * tmps_["56_aa_LLov"]("L,R,j,a");
        tmps_["111_aaaa_LLvvoo"]("L,R,a,b,i,j") += tmps_["100_aa_Lvv"]("L,a,b") * tmps_["61_aa_Loo"]("R,i,j");
        tmps_["111_aaaa_LLvvoo"]("L,R,a,b,i,j") += Id["aa_oo"]("i,j") * tmps_["99_aa_LLvv"]("L,R,a,b");
        tmps_["111_aaaa_LLvvoo"]("L,R,a,b,i,j") -= r0("R") * tmps_["107_aaaa_Lvoov"]("L,b,i,j,a");
        tmps_["111_aaaa_LLvvoo"]("L,R,a,b,i,j") -= r0("R") * tmps_["108_aaaa_Lvoov"]("L,b,i,j,a");
        tmps_["111_aaaa_LLvvoo"]("L,R,a,b,i,j") += Id["aa_oo"]("i,j") * tmps_["102_aa_LLvv"]("L,R,a,b");
        tmps_["102_aa_LLvv"].~TArrayD();
        tmps_["100_aa_Lvv"].~TArrayD();
        tmps_["99_aa_LLvv"].~TArrayD();
        tmps_["56_aa_LLov"].~TArrayD();

        // D_ovov_aaaa += -1.00 d_aa(i,j) r1_bb(c,l) t1_aa(b,k) l2_abab(k,l,a,c)
        //             += +1.00 d_aa(i,j) r1_aa(c,l) t1_aa(b,k) l2_aaaa(l,k,a,c)
        //             += -0.50 d_aa(i,j) r0 t2_aaaa(b,c,k,l) l2_aaaa(k,l,a,c)
        //             += -0.50 d_aa(i,j) r0 t2_abab(b,c,k,l) l2_abab(k,l,a,c)
        //             += -0.50 d_aa(i,j) r0 t2_abab(b,c,l,k) l2_abab(l,k,a,c)
        //             += +1.00 r0 t1_aa(b,k) t1_aa(c,i) l2_aaaa(k,j,a,c)
        //             += +1.00 r1_aa(b,k) t1_aa(c,i) l2_aaaa(k,j,a,c)
        //             += +1.00 r1_aa(b,i) l1_aa(j,a)
        //             += -0.50 d_aa(i,j) r2_abab(b,c,l,k) l2_abab(l,k,a,c)
        //             += -0.50 d_aa(i,j) r2_abab(b,c,k,l) l2_abab(k,l,a,c)
        //             += -0.50 d_aa(i,j) r2_aaaa(b,c,l,k) l2_aaaa(l,k,a,c)
        //             += +1.00 r0 t1_aa(b,i) l1_aa(j,a)
        //             += -1.00 r1_aa(c,k) t1_aa(b,i) l2_aaaa(k,j,a,c)
        //             += +1.00 r1_bb(c,k) t1_aa(b,i) l2_abab(j,k,a,c)
        //             += -1.00 d_aa(i,j) r0 t1_aa(b,k) l1_aa(k,a)
        //             += -1.00 d_aa(i,j) r1_aa(b,k) l1_aa(k,a)
        //             += +1.00 r0 t2_aaaa(b,c,k,i) l2_aaaa(k,j,a,c)
        //             += +1.00 r0 t2_abab(b,c,i,k) l2_abab(j,k,a,c)
        D_ovov_aaaa("L,R,i,a,j,b") -= tmps_["111_aaaa_LLvvoo"]("L,R,a,b,i,j");

        // D_ovvo_aaaa += +1.00 d_aa(i,j) r1_bb(c,l) t1_aa(b,k) l2_abab(k,l,a,c)
        //             += -1.00 d_aa(i,j) r1_aa(c,l) t1_aa(b,k) l2_aaaa(l,k,a,c)
        //             += +0.50 d_aa(i,j) r0 t2_aaaa(b,c,k,l) l2_aaaa(k,l,a,c)
        //             += +0.50 d_aa(i,j) r0 t2_abab(b,c,k,l) l2_abab(k,l,a,c)
        //             += +0.50 d_aa(i,j) r0 t2_abab(b,c,l,k) l2_abab(l,k,a,c)
        //             += -1.00 r0 t1_aa(b,k) t1_aa(c,i) l2_aaaa(k,j,a,c)
        //             += -1.00 r1_aa(b,k) t1_aa(c,i) l2_aaaa(k,j,a,c)
        //             += -1.00 r1_aa(b,i) l1_aa(j,a)
        //             += +0.50 d_aa(i,j) r2_abab(b,c,l,k) l2_abab(l,k,a,c)
        //             += +0.50 d_aa(i,j) r2_abab(b,c,k,l) l2_abab(k,l,a,c)
        //             += +0.50 d_aa(i,j) r2_aaaa(b,c,l,k) l2_aaaa(l,k,a,c)
        //             += -1.00 r0 t1_aa(b,i) l1_aa(j,a)
        //             += +1.00 r1_aa(c,k) t1_aa(b,i) l2_aaaa(k,j,a,c)
        //             += -1.00 r1_bb(c,k) t1_aa(b,i) l2_abab(j,k,a,c)
        //             += +1.00 d_aa(i,j) r0 t1_aa(b,k) l1_aa(k,a)
        //             += +1.00 d_aa(i,j) r1_aa(b,k) l1_aa(k,a)
        //             += -1.00 r0 t2_aaaa(b,c,k,i) l2_aaaa(k,j,a,c)
        //             += -1.00 r0 t2_abab(b,c,i,k) l2_abab(j,k,a,c)
        D_ovvo_aaaa("L,R,i,a,b,j") += tmps_["111_aaaa_LLvvoo"]("L,R,a,b,i,j");
        tmps_["111_aaaa_LLvvoo"].~TArrayD();

        // flops: o0v0L2  = o0v0L2
        //  mems: o0v0L2  = o0v0L2
        tmps_["112_LL"]("L,R")  = l0("L") * r0("R");

        // D_oovv_aaaa += +1.00 r0 t2_aaaa(b,a,i,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") += t2["aaaa"]("b,a,i,j") * tmps_["112_LL"]("L,R");

        // D_oovv_abab += -1.00 r0 t2_abab(a,b,i,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,b,i,j") * tmps_["112_LL"]("L,R");

        // D_oovv_bbbb += +1.00 r0 t2_bbbb(b,a,i,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += t2["bbbb"]("b,a,i,j") * tmps_["112_LL"]("L,R");

        // flops: o4v0  = o4v0
        //  mems: o4v0  = o4v0
        tmps_["113_aaaa_oooo"]("i,j,k,l")  = Id["aa_oo"]("i,j") * Id["aa_oo"]("k,l");

        // D_oooo_aaaa += +1.00 P(i,j) d_aa(i,l) d_aa(j,k) r0 l0
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l")  = tmps_["112_LL"]("L,R") * tmps_["113_aaaa_oooo"]("i,l,j,k");
        D_oooo_aaaa("L,R,i,j,k,l") += tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l");
        D_oooo_aaaa("L,R,i,j,k,l") -= tmps_["perm_aaaa_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_aaaa_LLoooo"].~TArrayD();

        // flops: o4v0  = o4v0
        //  mems: o4v0  = o4v0
        tmps_["114_bbbb_oooo"]("i,j,k,l")  = Id["bb_oo"]("i,j") * Id["bb_oo"]("k,l");

        // D_oooo_bbbb += +1.00 P(i,j) d_bb(i,l) d_bb(j,k) r0 l0
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        tmps_["perm_bbbb_LLoooo"]("L,R,i,j,k,l")  = tmps_["112_LL"]("L,R") * tmps_["114_bbbb_oooo"]("i,l,j,k");
        D_oooo_bbbb("L,R,i,j,k,l") += tmps_["perm_bbbb_LLoooo"]("L,R,i,j,k,l");
        D_oooo_bbbb("L,R,i,j,k,l") -= tmps_["perm_bbbb_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_bbbb_LLoooo"].~TArrayD();

        // flops: o4v0  = o4v0
        //  mems: o4v0  = o4v0
        tmps_["115_aabb_oooo"]("i,j,k,l")  = Id["aa_oo"]("i,j") * Id["bb_oo"]("k,l");

        // D_oooo_abab += -1.00 d_aa(i,k) d_bb(j,l) r0 l0
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") -= tmps_["112_LL"]("L,R") * tmps_["115_aabb_oooo"]("i,k,j,l");
        tmps_["112_LL"].~TArrayD();

        // flops: o4v0L2  = o4v2L2
        //  mems: o4v0L2  = o4v0L2
        tmps_["116_abab_LLoooo"]("L,R,i,j,k,l")  = l2["abab"]("L,i,j,a,b") * r2["abab"]("R,a,b,k,l");

        // D_oooo_abab += -0.50 r2_abab(a,b,i,j) l2_abab(k,l,a,b)
        //             += -0.50 r2_abab(b,a,i,j) l2_abab(k,l,b,a)
        D_oooo_abab("L,R,i,j,k,l") -= tmps_["116_abab_LLoooo"]("L,R,k,l,i,j");

        // D_oovv_abab += -0.25 r2_abab(c,d,i,j) t2_abab(a,b,k,l) l2_abab(k,l,c,d)
        //             += -0.25 r2_abab(c,d,i,j) t2_abab(a,b,l,k) l2_abab(l,k,c,d)
        //             += -0.25 r2_abab(d,c,i,j) t2_abab(a,b,k,l) l2_abab(k,l,d,c)
        //             += -0.25 r2_abab(d,c,i,j) t2_abab(a,b,l,k) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,b,k,l") * tmps_["116_abab_LLoooo"]("L,R,k,l,i,j");

        // D_oovv_abab += -0.50 r2_abab(c,d,i,j) t1_aa(a,k) t1_bb(b,l) l2_abab(k,l,c,d)
        //             += -0.50 r2_abab(d,c,i,j) t1_aa(a,k) t1_bb(b,l) l2_abab(k,l,d,c)
        // flops: o2v2L2 += o4v1L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["bb"]("b,l") * tmps_["116_abab_LLoooo"]("L,R,k,l,i,j") * t1["aa"]("a,k");

        // flops: o3v1L1  = o3v2L1 o4v2L1
        //  mems: o3v1L1  = o3v1L1 o3v1L1
        tmps_["117_bbaa_Loovo"]("L,i,j,a,k")  = l2["bbbb"]("L,l,j,c,b") * t1["bb"]("c,i") * t2["abab"]("a,b,k,l");

        // D_oovv_abab += -1.00 r0 t2_abab(a,c,i,k) t1_bb(b,l) t1_bb(d,j) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["bb"]("b,l") * tmps_["117_bbaa_Loovo"]("L,j,l,a,i") * r0("R");

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["118_aabb_Lvooo"]("L,a,i,j,k")  = l1["bb"]("L,k,b") * t2["abab"]("a,b,i,j");

        // D_oovv_abab += +1.00 r1_bb(b,k) t2_abab(a,c,i,j) l1_bb(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += r1["bb"]("R,b,k") * tmps_["118_aabb_Lvooo"]("L,a,i,j,k");

        // D_oovv_abab += +1.00 r0 t2_abab(a,c,i,j) t1_bb(b,k) l1_bb(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1["bb"]("b,k") * tmps_["118_aabb_Lvooo"]("L,a,i,j,k") * r0("R");

        // flops: o4v0L2  = o3v2L1 o4v1L2 o3v2L1 o4v1L2 o4v0L2
        //  mems: o4v0L2  = o3v1L1 o4v0L2 o3v1L1 o4v0L2 o4v0L2
        tmps_["119_aabb_LLoooo"]("L,R,i,j,k,l")  = l2["abab"]("L,j,k,b,d") * t1["bb"]("d,l") * r1["aa"]("R,b,i");
        tmps_["119_aabb_LLoooo"]("L,R,i,j,k,l") += l2["abab"]("L,j,k,c,a") * t1["aa"]("c,i") * r1["bb"]("R,a,l");

        // D_oooo_abab += -1.00 r1_bb(b,j) t1_aa(a,i) l2_abab(k,l,a,b)
        //             += -1.00 r1_aa(b,i) t1_bb(a,j) l2_abab(k,l,b,a)
        D_oooo_abab("L,R,i,j,k,l") -= tmps_["119_aabb_LLoooo"]("L,R,i,k,l,j");

        // D_oovv_abab += -1.00 r1_bb(d,j) t1_aa(a,k) t1_bb(b,l) t1_aa(c,i) l2_abab(k,l,c,d)
        //             += -1.00 r1_aa(d,i) t1_aa(a,k) t1_bb(b,l) t1_bb(c,j) l2_abab(k,l,d,c)
        // flops: o2v2L2 += o4v1L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["bb"]("b,l") * tmps_["119_aabb_LLoooo"]("L,R,i,k,l,j") * t1["aa"]("a,k");

        // flops: o3v1L1  = o3v2L1 o3v2L1 o3v1L1
        //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1
        tmps_["120_aabb_Lvooo"]("L,a,i,j,k")  = t1["aa"]("c,i") * tmps_["12_abba_Lvoov"]("L,a,k,j,c");
        tmps_["120_aabb_Lvooo"]("L,a,i,j,k") += t1["bb"]("b,k") * tmps_["13_aabb_Lvoov"]("L,a,i,j,b");
        tmps_["12_abba_Lvoov"].~TArrayD();

        // D_oovv_abab += -1.00 r1_bb(b,l) t2_aaaa(a,c,k,i) t1_bb(d,j) l2_abab(k,l,c,d)
        //             += -1.00 r1_bb(b,l) t2_abab(a,c,k,j) t1_aa(d,i) l2_abab(k,l,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r1["bb"]("R,b,l") * tmps_["120_aabb_Lvooo"]("L,a,i,l,j");

        // D_oovv_abab += -1.00 r0 t2_aaaa(a,c,k,i) t1_bb(b,l) t1_bb(d,j) l2_abab(k,l,c,d)
        //             += -1.00 r0 t2_abab(a,c,k,j) t1_bb(b,l) t1_aa(d,i) l2_abab(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["bb"]("b,l") * tmps_["120_aabb_Lvooo"]("L,a,i,l,j") * r0("R");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["121_aabb_Lvoov"]("L,a,i,j,b")  = l2["bbbb"]("L,j,k,c,b") * t2["abab"]("a,c,i,k");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_abab(a,d,i,l) t2_abab(b,c,j,k) l2_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = r2["abab"]("R,a,d,i,l") * tmps_["121_aabb_Lvoov"]("L,b,j,l,d");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) r0 t2_abab(a,c,i,k) t2_abab(b,d,j,l) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["abab"]("a,c,i,k") * tmps_["121_aabb_Lvoov"]("L,b,j,k,c") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["122_aa_Lvv"]("L,a,b")  = l2["aaaa"]("L,i,j,c,b") * t2["aaaa"]("a,c,i,j");

        // D_oovv_aaaa += -0.50 P(a,b) r2_aaaa(a,d,i,j) t2_aaaa(b,c,k,l) l2_aaaa(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * r2["aaaa"]("R,a,d,i,j") * tmps_["122_aa_Lvv"]("L,b,d");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += -0.50 r0 t2_aaaa(a,c,i,j) t2_aaaa(b,d,k,l) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") -= 0.50 * t2["aaaa"]("a,c,i,j") * tmps_["122_aa_Lvv"]("L,b,c") * r0("R");

        // D_oovv_abab += -0.50 r2_abab(d,b,i,j) t2_aaaa(a,c,k,l) l2_aaaa(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * r2["abab"]("R,d,b,i,j") * tmps_["122_aa_Lvv"]("L,a,d");

        // flops: o2v0L1  = o3v2L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["123_aa_Loo"]("L,i,j")  = l2["aaaa"]("L,j,k,a,b") * t2["aaaa"]("a,b,k,i");

        // D_oovv_aaaa += -0.50 P(i,j) r2_aaaa(b,a,l,i) t2_aaaa(d,c,k,j) l2_aaaa(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * r2["aaaa"]("R,b,a,l,i") * tmps_["123_aa_Loo"]("L,j,l");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -0.50 r2_abab(a,b,l,j) t2_aaaa(d,c,k,i) l2_aaaa(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * r2["abab"]("R,a,b,l,j") * tmps_["123_aa_Loo"]("L,i,l");

        // flops: o1v1L2  = o2v2L2 o2v1L1 o2v1L2 o1v1L2 o2v2L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o2v0L1 o1v1L2 o1v1L2 o1v1L2 o1v1L2
        tmps_["125_aa_LLov"]("L,R,i,a")  = -1.00 * l1["bb"]("L,k,c") * r2["abab"]("R,a,c,i,k");
        tmps_["125_aa_LLov"]("L,R,i,a") += l1["aa"]("L,j,b") * t1["aa"]("b,i") * r1["aa"]("R,a,j");
        tmps_["125_aa_LLov"]("L,R,i,a") += l1["aa"]("L,j,b") * r2["aaaa"]("R,a,b,j,i");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r1_aa(a,k) t1_aa(b,j) t1_aa(c,i) l1_aa(k,c)
        //             += -1.00 P(i,j) P(a,b) r2_abab(a,c,i,k) t1_aa(b,j) l1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r2_aaaa(a,c,k,i) t1_aa(b,j) l1_aa(k,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("b,j") * tmps_["125_aa_LLov"]("L,R,i,a");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +1.00 r1_aa(a,k) t1_bb(b,j) t1_aa(c,i) l1_aa(k,c)
        //             += -1.00 r2_abab(a,c,i,k) t1_bb(b,j) l1_bb(k,c)
        //             += +1.00 r2_aaaa(a,c,k,i) t1_bb(b,j) l1_aa(k,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1["bb"]("b,j") * tmps_["125_aa_LLov"]("L,R,i,a");

        // flops: o0v0L2  = o2v2L2 o1v1L2 o1v1L2 o0v0L2 o0v0L2 o2v2L2 o0v0L2 o2v2L2 o0v0L2
        //  mems: o0v0L2  = o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2 o0v0L2
        tmps_["126_LL"]("L,R")  = 0.25 * l2["aaaa"]("L,i,l,a,d") * r2["aaaa"]("R,a,d,i,l");
        tmps_["126_LL"]("L,R") += l1["aa"]("L,l,a") * r1["aa"]("R,a,l");
        tmps_["126_LL"]("L,R") += l1["bb"]("L,j,c") * r1["bb"]("R,c,j");
        tmps_["126_LL"]("L,R") += l2["abab"]("L,i,j,a,b") * r2["abab"]("R,a,b,i,j");
        tmps_["126_LL"]("L,R") += 0.25 * l2["bbbb"]("L,k,j,c,b") * r2["bbbb"]("R,c,b,k,j");

        // D_oooo_bbbb += +1.00 P(i,j) d_bb(i,l) d_bb(j,k) r1_bb(a,m) l1_bb(m,a)
        //             += +1.00 P(i,j) d_bb(i,l) d_bb(j,k) r1_aa(a,m) l1_aa(m,a)
        //             += +0.25 P(i,j) d_bb(i,l) d_bb(j,k) r2_aaaa(a,b,n,m) l2_aaaa(n,m,a,b)
        //             += +0.25 P(i,j) d_bb(i,l) d_bb(j,k) r2_abab(a,b,n,m) l2_abab(n,m,a,b)
        //             += +0.25 P(i,j) d_bb(i,l) d_bb(j,k) r2_abab(a,b,m,n) l2_abab(m,n,a,b)
        //             += +0.25 P(i,j) d_bb(i,l) d_bb(j,k) r2_abab(b,a,n,m) l2_abab(n,m,b,a)
        //             += +0.25 P(i,j) d_bb(i,l) d_bb(j,k) r2_abab(b,a,m,n) l2_abab(m,n,b,a)
        //             += +0.25 P(i,j) d_bb(i,l) d_bb(j,k) r2_bbbb(a,b,n,m) l2_bbbb(n,m,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        tmps_["perm_bbbb_LLoooo"]("L,R,i,j,k,l")  = tmps_["126_LL"]("L,R") * tmps_["114_bbbb_oooo"]("i,l,j,k");
        D_oooo_bbbb("L,R,i,j,k,l") += tmps_["perm_bbbb_LLoooo"]("L,R,i,j,k,l");
        D_oooo_bbbb("L,R,i,j,k,l") -= tmps_["perm_bbbb_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_bbbb_LLoooo"].~TArrayD();
        tmps_["114_bbbb_oooo"].~TArrayD();

        // D_oooo_abab += -1.00 d_aa(i,k) d_bb(j,l) r1_bb(a,m) l1_bb(m,a)
        //             += -1.00 d_aa(i,k) d_bb(j,l) r1_aa(a,m) l1_aa(m,a)
        //             += -0.25 d_aa(i,k) d_bb(j,l) r2_aaaa(a,b,n,m) l2_aaaa(n,m,a,b)
        //             += -0.25 d_aa(i,k) d_bb(j,l) r2_abab(a,b,n,m) l2_abab(n,m,a,b)
        //             += -0.25 d_aa(i,k) d_bb(j,l) r2_abab(a,b,m,n) l2_abab(m,n,a,b)
        //             += -0.25 d_aa(i,k) d_bb(j,l) r2_abab(b,a,n,m) l2_abab(n,m,b,a)
        //             += -0.25 d_aa(i,k) d_bb(j,l) r2_abab(b,a,m,n) l2_abab(m,n,b,a)
        //             += -0.25 d_aa(i,k) d_bb(j,l) r2_bbbb(a,b,n,m) l2_bbbb(n,m,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") -= tmps_["126_LL"]("L,R") * tmps_["115_aabb_oooo"]("i,k,j,l");
        tmps_["115_aabb_oooo"].~TArrayD();

        // D_oovv_aaaa += +1.00 r1_bb(c,k) t2_aaaa(b,a,i,j) l1_bb(k,c)
        //             += +1.00 r1_aa(c,k) t2_aaaa(b,a,i,j) l1_aa(k,c)
        //             += +0.25 r2_aaaa(c,d,l,k) t2_aaaa(b,a,i,j) l2_aaaa(l,k,c,d)
        //             += +0.25 r2_abab(c,d,l,k) t2_aaaa(b,a,i,j) l2_abab(l,k,c,d)
        //             += +0.25 r2_abab(c,d,k,l) t2_aaaa(b,a,i,j) l2_abab(k,l,c,d)
        //             += +0.25 r2_abab(d,c,l,k) t2_aaaa(b,a,i,j) l2_abab(l,k,d,c)
        //             += +0.25 r2_abab(d,c,k,l) t2_aaaa(b,a,i,j) l2_abab(k,l,d,c)
        //             += +0.25 r2_bbbb(c,d,l,k) t2_aaaa(b,a,i,j) l2_bbbb(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") += t2["aaaa"]("b,a,i,j") * tmps_["126_LL"]("L,R");

        // D_oovv_abab += -1.00 r1_bb(c,k) t2_abab(a,b,i,j) l1_bb(k,c)
        //             += -1.00 r1_aa(c,k) t2_abab(a,b,i,j) l1_aa(k,c)
        //             += -0.25 r2_aaaa(c,d,l,k) t2_abab(a,b,i,j) l2_aaaa(l,k,c,d)
        //             += -0.25 r2_abab(c,d,l,k) t2_abab(a,b,i,j) l2_abab(l,k,c,d)
        //             += -0.25 r2_abab(c,d,k,l) t2_abab(a,b,i,j) l2_abab(k,l,c,d)
        //             += -0.25 r2_abab(d,c,l,k) t2_abab(a,b,i,j) l2_abab(l,k,d,c)
        //             += -0.25 r2_abab(d,c,k,l) t2_abab(a,b,i,j) l2_abab(k,l,d,c)
        //             += -0.25 r2_bbbb(c,d,l,k) t2_abab(a,b,i,j) l2_bbbb(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,b,i,j") * tmps_["126_LL"]("L,R");

        // D_oovv_abab += -1.00 r1_bb(c,k) t1_aa(a,i) t1_bb(b,j) l1_bb(k,c)
        //             += -1.00 r1_aa(c,k) t1_aa(a,i) t1_bb(b,j) l1_aa(k,c)
        //             += -0.25 r2_aaaa(c,d,l,k) t1_aa(a,i) t1_bb(b,j) l2_aaaa(l,k,c,d)
        //             += -0.25 r2_abab(c,d,l,k) t1_aa(a,i) t1_bb(b,j) l2_abab(l,k,c,d)
        //             += -0.25 r2_abab(c,d,k,l) t1_aa(a,i) t1_bb(b,j) l2_abab(k,l,c,d)
        //             += -0.25 r2_abab(d,c,l,k) t1_aa(a,i) t1_bb(b,j) l2_abab(l,k,d,c)
        //             += -0.25 r2_abab(d,c,k,l) t1_aa(a,i) t1_bb(b,j) l2_abab(k,l,d,c)
        //             += -0.25 r2_bbbb(c,d,l,k) t1_aa(a,i) t1_bb(b,j) l2_bbbb(l,k,c,d)
        // flops: o2v2L2 += o2v2 o2v2L2
        //  mems: o2v2L2 += o2v2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["aa"]("a,i") * t1["bb"]("b,j") * tmps_["126_LL"]("L,R");

        // D_oovv_bbbb += +1.00 r1_bb(c,k) t2_bbbb(b,a,i,j) l1_bb(k,c)
        //             += +1.00 r1_aa(c,k) t2_bbbb(b,a,i,j) l1_aa(k,c)
        //             += +0.25 r2_aaaa(c,d,l,k) t2_bbbb(b,a,i,j) l2_aaaa(l,k,c,d)
        //             += +0.25 r2_abab(c,d,l,k) t2_bbbb(b,a,i,j) l2_abab(l,k,c,d)
        //             += +0.25 r2_abab(c,d,k,l) t2_bbbb(b,a,i,j) l2_abab(k,l,c,d)
        //             += +0.25 r2_abab(d,c,l,k) t2_bbbb(b,a,i,j) l2_abab(l,k,d,c)
        //             += +0.25 r2_abab(d,c,k,l) t2_bbbb(b,a,i,j) l2_abab(k,l,d,c)
        //             += +0.25 r2_bbbb(c,d,l,k) t2_bbbb(b,a,i,j) l2_bbbb(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += t2["bbbb"]("b,a,i,j") * tmps_["126_LL"]("L,R");

        // D_oovv_bbbb += -1.00 P(i,j) r1_bb(c,k) t1_bb(a,i) t1_bb(b,j) l1_bb(k,c)
        //             += -1.00 P(i,j) r1_aa(c,k) t1_bb(a,i) t1_bb(b,j) l1_aa(k,c)
        //             += -0.25 P(i,j) r2_aaaa(c,d,l,k) t1_bb(a,i) t1_bb(b,j) l2_aaaa(l,k,c,d)
        //             += -0.25 P(i,j) r2_abab(c,d,l,k) t1_bb(a,i) t1_bb(b,j) l2_abab(l,k,c,d)
        //             += -0.25 P(i,j) r2_abab(c,d,k,l) t1_bb(a,i) t1_bb(b,j) l2_abab(k,l,c,d)
        //             += -0.25 P(i,j) r2_abab(d,c,l,k) t1_bb(a,i) t1_bb(b,j) l2_abab(l,k,d,c)
        //             += -0.25 P(i,j) r2_abab(d,c,k,l) t1_bb(a,i) t1_bb(b,j) l2_abab(k,l,d,c)
        //             += -0.25 P(i,j) r2_bbbb(c,d,l,k) t1_bb(a,i) t1_bb(b,j) l2_bbbb(l,k,c,d)
        // flops: o2v2L2 += o2v2 o2v2L2
        //  mems: o2v2L2 += o2v2 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("a,i") * t1["bb"]("b,j") * tmps_["126_LL"]("L,R");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oooo_aaaa += +1.00 P(i,j) d_aa(i,l) d_aa(j,k) r1_bb(a,m) l1_bb(m,a)
        //             += +1.00 P(i,j) d_aa(i,l) d_aa(j,k) r1_aa(a,m) l1_aa(m,a)
        //             += +0.25 P(i,j) d_aa(i,l) d_aa(j,k) r2_aaaa(a,b,n,m) l2_aaaa(n,m,a,b)
        //             += +0.25 P(i,j) d_aa(i,l) d_aa(j,k) r2_abab(a,b,n,m) l2_abab(n,m,a,b)
        //             += +0.25 P(i,j) d_aa(i,l) d_aa(j,k) r2_abab(a,b,m,n) l2_abab(m,n,a,b)
        //             += +0.25 P(i,j) d_aa(i,l) d_aa(j,k) r2_abab(b,a,n,m) l2_abab(n,m,b,a)
        //             += +0.25 P(i,j) d_aa(i,l) d_aa(j,k) r2_abab(b,a,m,n) l2_abab(m,n,b,a)
        //             += +0.25 P(i,j) d_aa(i,l) d_aa(j,k) r2_bbbb(a,b,n,m) l2_bbbb(n,m,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l")  = tmps_["126_LL"]("L,R") * tmps_["113_aaaa_oooo"]("i,l,j,k");
        D_oooo_aaaa("L,R,i,j,k,l") += tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l");
        D_oooo_aaaa("L,R,i,j,k,l") -= tmps_["perm_aaaa_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_aaaa_LLoooo"].~TArrayD();
        tmps_["113_aaaa_oooo"].~TArrayD();

        // D_oovv_aaaa += -1.00 P(i,j) r1_bb(c,k) t1_aa(a,i) t1_aa(b,j) l1_bb(k,c)
        //             += -1.00 P(i,j) r1_aa(c,k) t1_aa(a,i) t1_aa(b,j) l1_aa(k,c)
        //             += -0.25 P(i,j) r2_aaaa(c,d,l,k) t1_aa(a,i) t1_aa(b,j) l2_aaaa(l,k,c,d)
        //             += -0.25 P(i,j) r2_abab(c,d,l,k) t1_aa(a,i) t1_aa(b,j) l2_abab(l,k,c,d)
        //             += -0.25 P(i,j) r2_abab(c,d,k,l) t1_aa(a,i) t1_aa(b,j) l2_abab(k,l,c,d)
        //             += -0.25 P(i,j) r2_abab(d,c,l,k) t1_aa(a,i) t1_aa(b,j) l2_abab(l,k,d,c)
        //             += -0.25 P(i,j) r2_abab(d,c,k,l) t1_aa(a,i) t1_aa(b,j) l2_abab(k,l,d,c)
        //             += -0.25 P(i,j) r2_bbbb(c,d,l,k) t1_aa(a,i) t1_aa(b,j) l2_bbbb(l,k,c,d)
        // flops: o2v2L2 += o2v2 o2v2L2
        //  mems: o2v2L2 += o2v2 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("a,i") * t1["aa"]("b,j") * tmps_["126_LL"]("L,R");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o1v1L2  = o2v2L2 o2v1L2 o2v2L2 o2v1L2 o2v2L2 o2v1L2 o1v1L2 o1v1L2 o2v2L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v2L2 o1v1L2 o1v2L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2
        tmps_["127_aa_LLvo"]("L,R,a,i")  = -1.00 * t2["aaaa"]("a,b,j,i") * tmps_["55_aa_LLov"]("L,R,j,b");
        tmps_["127_aa_LLvo"]("L,R,a,i") -= 0.50 * r1["aa"]("R,a,k") * tmps_["123_aa_Loo"]("L,i,k");
        tmps_["127_aa_LLvo"]("L,R,a,i") -= t2["abab"]("a,e,i,l") * tmps_["39_bb_LLov"]("L,R,l,e");
        tmps_["127_aa_LLvo"]("L,R,a,i") += r1["aa"]("R,a,k") * tmps_["58_aa_Loo"]("L,i,k");
        tmps_["127_aa_LLvo"]("L,R,a,i") += r1["bb"]("R,d,m") * tmps_["13_aabb_Lvoov"]("L,a,i,m,d");
        tmps_["127_aa_LLvo"]("L,R,a,i") += t1["aa"]("a,j") * tmps_["60_aa_LLoo"]("L,R,j,i");
        tmps_["127_aa_LLvo"]("L,R,a,i") += r1["bb"]("R,d,m") * tmps_["121_aabb_Lvoov"]("L,a,i,m,d");
        tmps_["127_aa_LLvo"]("L,R,a,i") += t1["aa"]("b,i") * tmps_["101_aa_LLvv"]("L,R,b,a");
        tmps_["127_aa_LLvo"]("L,R,a,i") -= 0.50 * r1["aa"]("R,c,i") * tmps_["122_aa_Lvv"]("L,a,c");
        tmps_["123_aa_Loo"].~TArrayD();
        tmps_["122_aa_Lvv"].~TArrayD();
        tmps_["121_aabb_Lvoov"].~TArrayD();
        tmps_["101_aa_LLvv"].~TArrayD();
        tmps_["39_bb_LLov"].~TArrayD();
        tmps_["13_aabb_Lvoov"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r1_aa(c,i) t1_aa(a,k) t1_aa(b,j) l1_aa(k,c)
        //             += +0.50 P(i,j) P(a,b) r2_abab(c,d,i,l) t1_aa(a,k) t1_aa(b,j) l2_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_abab(d,c,i,l) t1_aa(a,k) t1_aa(b,j) l2_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r2_aaaa(c,d,l,i) t1_aa(a,k) t1_aa(b,j) l2_aaaa(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_bb(d,l) t2_aaaa(a,c,k,i) t1_aa(b,j) l2_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,l) t1_aa(b,j) t2_abab(d,c,i,k) l2_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,l) t1_aa(b,j) t2_abab(c,d,i,k) l2_abab(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_bb(d,l) t2_abab(a,c,i,k) t1_aa(b,j) l2_bbbb(l,k,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_aa(d,l) t2_abab(a,c,i,k) t1_aa(b,j) l2_abab(l,k,d,c)
        //             += -0.50 P(i,j) P(a,b) r1_aa(a,l) t1_aa(b,j) t2_aaaa(d,c,k,i) l2_aaaa(l,k,d,c)
        //             += -1.00 P(i,j) P(a,b) r1_aa(d,l) t2_aaaa(a,c,k,i) t1_aa(b,j) l2_aaaa(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_abab(a,d,l,k) t1_aa(b,j) t1_aa(c,i) l2_abab(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_abab(a,d,k,l) t1_aa(b,j) t1_aa(c,i) l2_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_aaaa(a,d,l,k) t1_aa(b,j) t1_aa(c,i) l2_aaaa(l,k,c,d)
        //             += -0.50 P(i,j) P(a,b) r1_aa(d,i) t2_aaaa(a,c,k,l) t1_aa(b,j) l2_aaaa(k,l,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("b,j") * tmps_["127_aa_LLvo"]("L,R,a,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +1.00 r1_aa(c,i) t1_aa(a,k) t1_bb(b,j) l1_aa(k,c)
        //             += +0.50 r2_abab(c,d,i,l) t1_aa(a,k) t1_bb(b,j) l2_abab(k,l,c,d)
        //             += +0.50 r2_abab(d,c,i,l) t1_aa(a,k) t1_bb(b,j) l2_abab(k,l,d,c)
        //             += +0.50 r2_aaaa(c,d,l,i) t1_aa(a,k) t1_bb(b,j) l2_aaaa(l,k,c,d)
        //             += +1.00 r1_bb(d,l) t2_aaaa(a,c,k,i) t1_bb(b,j) l2_abab(k,l,c,d)
        //             += +0.50 r1_aa(a,l) t1_bb(b,j) t2_abab(d,c,i,k) l2_abab(l,k,d,c)
        //             += +0.50 r1_aa(a,l) t1_bb(b,j) t2_abab(c,d,i,k) l2_abab(l,k,c,d)
        //             += +1.00 r1_bb(d,l) t2_abab(a,c,i,k) t1_bb(b,j) l2_bbbb(l,k,c,d)
        //             += -1.00 r1_aa(d,l) t2_abab(a,c,i,k) t1_bb(b,j) l2_abab(l,k,d,c)
        //             += -0.50 r1_aa(a,l) t1_bb(b,j) t2_aaaa(d,c,k,i) l2_aaaa(l,k,d,c)
        //             += -1.00 r1_aa(d,l) t2_aaaa(a,c,k,i) t1_bb(b,j) l2_aaaa(l,k,c,d)
        //             += +0.50 r2_abab(a,d,l,k) t1_bb(b,j) t1_aa(c,i) l2_abab(l,k,c,d)
        //             += +0.50 r2_abab(a,d,k,l) t1_bb(b,j) t1_aa(c,i) l2_abab(k,l,c,d)
        //             += +0.50 r2_aaaa(a,d,l,k) t1_bb(b,j) t1_aa(c,i) l2_aaaa(l,k,c,d)
        //             += -0.50 r1_aa(d,i) t2_aaaa(a,c,k,l) t1_bb(b,j) l2_aaaa(k,l,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1["bb"]("b,j") * tmps_["127_aa_LLvo"]("L,R,a,i");

        // flops: o1v1L1  = o2v2L1 o2v2L1 o1v1L1 o2v1L1 o2v1L1 o1v1L1 o1v2L1 o2v1L1 o1v2L1 o1v1L1 o1v1L1 o1v1L1 o2v1L1 o1v1L1
        //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o2v0L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
        tmps_["128_aa_Lov"]("L,i,a")  = -1.00 * l1["bb"]("L,m,d") * t2["abab"]("a,d,i,m");
        tmps_["128_aa_Lov"]("L,i,a") += l1["aa"]("L,k,c") * t2["aaaa"]("a,c,k,i");
        tmps_["128_aa_Lov"]("L,i,a") += l1["aa"]("L,k,c") * t1["aa"]("c,i") * t1["aa"]("a,k");
        tmps_["128_aa_Lov"]("L,i,a") += 0.50 * tmps_["97_aa_Lvv"]("L,b,a") * t1["aa"]("b,i");
        tmps_["128_aa_Lov"]("L,i,a") += t1["aa"]("a,j") * tmps_["58_aa_Loo"]("L,i,j");
        tmps_["128_aa_Lov"]("L,i,a") += tmps_["98_aa_Lvv"]("L,b,a") * t1["aa"]("b,i");
        tmps_["128_aa_Lov"]("L,i,a") += 0.50 * t1["aa"]("a,j") * tmps_["57_aa_Loo"]("L,i,j");
        tmps_["97_aa_Lvv"].~TArrayD();

        // D_oovv_aaaa += +0.50 P(i,j) P(a,b) r0 t2_abab(a,c,k,l) t1_aa(b,j) t1_aa(d,i) l2_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t2_abab(a,c,l,k) t1_aa(b,j) t1_aa(d,i) l2_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_aa(a,l) t1_aa(b,j) t2_abab(d,c,i,k) l2_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_aa(a,l) t1_aa(b,j) t2_abab(c,d,i,k) l2_abab(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r0 t2_aaaa(a,c,k,l) t1_aa(b,j) t1_aa(d,i) l2_aaaa(k,l,d,c)
        //             += -1.00 P(i,j) P(a,b) r0 t2_abab(a,c,i,k) t1_aa(b,j) l1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r0 t2_aaaa(a,c,k,i) t1_aa(b,j) l1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r0 t1_aa(a,k) t1_aa(b,j) t1_aa(c,i) l1_aa(k,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_aa(a,l) t1_aa(b,j) t2_aaaa(d,c,k,i) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["109_aa_Lvo"]("R,b,j") * tmps_["128_aa_Lov"]("L,i,a");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +0.50 P(i,j) P(a,b) r1_aa(a,i) t2_abab(b,c,k,l) t1_aa(d,j) l2_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,i) t2_abab(b,c,l,k) t1_aa(d,j) l2_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,i) t1_aa(b,l) t2_abab(d,c,j,k) l2_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,i) t1_aa(b,l) t2_abab(c,d,j,k) l2_abab(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,i) t2_aaaa(b,c,k,l) t1_aa(d,j) l2_aaaa(k,l,d,c)
        //             += -1.00 P(i,j) P(a,b) r1_aa(a,i) t2_abab(b,c,j,k) l1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_aa(a,i) t2_aaaa(b,c,k,j) l1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_aa(a,i) t1_aa(b,k) t1_aa(c,j) l1_aa(k,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,i) t1_aa(b,l) t2_aaaa(d,c,k,j) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = r1["aa"]("R,a,i") * tmps_["128_aa_Lov"]("L,j,b");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();
    }
} // hilbert
#endif