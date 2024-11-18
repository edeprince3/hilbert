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


    void EOM_EE_QED_RDM_21::rdm2_21_4() {

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


        // flops: o1v3L1  = o2v3L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["143_aaaa_Lvovv"]("L,a,i,b,c")  = t1["aa"]("a,j") * l2_1["aaaa"]("L,j,i,b,c");

        // D_vvvv_aaaa += -1.00 r0_1 t1_aa(c,i) t1_aa(d,j) l2_1_aaaa(i,j,a,b)
        // flops: o0v4L2 += o1v4L1 o0v4L2
        //  mems: o0v4L2 += o0v4L1 o0v4L2
        D_vvvv_aaaa("L,R,a,b,c,d") -= tmps_["143_aaaa_Lvovv"]("L,c,j,a,b") * t1["aa"]("d,j") * r0_1("R");

        // flops: o1v3L1  = o2v3L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["144_aaaa_Lvovv"]("L,a,i,b,c")  = t1["aa"]("a,j") * l2["aaaa"]("L,j,i,b,c");

        // D_vvvv_aaaa += -1.00 r0 t1_aa(c,i) t1_aa(d,j) l2_aaaa(i,j,a,b)
        // flops: o0v4L2 += o1v4L1 o0v4L2
        //  mems: o0v4L2 += o0v4L1 o0v4L2
        D_vvvv_aaaa("L,R,a,b,c,d") -= tmps_["144_aaaa_Lvovv"]("L,c,j,a,b") * t1["aa"]("d,j") * r0("R");

        // flops: o1v3L1  = o2v3L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["145_aaaa_Lvovv"]("L,a,i,b,c")  = t1_1["aa"]("a,j") * l2_1["aaaa"]("L,j,i,b,c");

        // D_vvvv_aaaa += -1.00 P(c,d) r0 t1_aa(d,j) t1_1_aa(c,i) l2_1_aaaa(i,j,a,b)
        // flops: o0v4L2 += o1v4L1 o0v4L2
        //  mems: o0v4L2 += o0v4L1 o0v4L2
        tmps_["perm_aaaa_LLvvvv"]("L,R,a,b,c,d")  = tmps_["145_aaaa_Lvovv"]("L,c,j,a,b") * t1["aa"]("d,j") * r0("R");
        D_vvvv_aaaa("L,R,a,b,c,d") -= tmps_["perm_aaaa_LLvvvv"]("L,R,a,b,c,d");
        D_vvvv_aaaa("L,R,a,b,c,d") += tmps_["perm_aaaa_LLvvvv"]("L,R,a,b,d,c");
        tmps_["perm_aaaa_LLvvvv"].~TArrayD();

        // flops: o1v3L2  = o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        //  mems: o1v3L2  = o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        tmps_["146_aaaa_LLvovv"]("R,L,a,i,b,c")  = tmps_["143_aaaa_Lvovv"]("L,a,i,b,c") * r0_1("R");
        tmps_["146_aaaa_LLvovv"]("R,L,a,i,b,c") += tmps_["145_aaaa_Lvovv"]("L,a,i,b,c") * r0("R");
        tmps_["146_aaaa_LLvovv"]("R,L,a,i,b,c") += tmps_["144_aaaa_Lvovv"]("L,a,i,b,c") * r0("R");
        tmps_["145_aaaa_Lvovv"].~TArrayD();
        tmps_["144_aaaa_Lvovv"].~TArrayD();
        tmps_["143_aaaa_Lvovv"].~TArrayD();

        // D_vvov_aaaa += +1.00 r0_1 t1_aa(c,j) l2_1_aaaa(j,i,a,b)
        //             += +1.00 r0 t1_1_aa(c,j) l2_1_aaaa(j,i,a,b)
        //             += +1.00 r0 t1_aa(c,j) l2_aaaa(j,i,a,b)
        D_vvov_aaaa("L,R,a,b,i,c") += tmps_["146_aaaa_LLvovv"]("R,L,c,i,a,b");

        // D_vvvo_aaaa += -1.00 r0_1 t1_aa(c,j) l2_1_aaaa(j,i,a,b)
        //             += -1.00 r0 t1_1_aa(c,j) l2_1_aaaa(j,i,a,b)
        //             += -1.00 r0 t1_aa(c,j) l2_aaaa(j,i,a,b)
        D_vvvo_aaaa("L,R,a,b,c,i") -= tmps_["146_aaaa_LLvovv"]("R,L,c,i,a,b");
        tmps_["146_aaaa_LLvovv"].~TArrayD();

        // flops: o1v3L1  = o2v3L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["147_bbbb_Lvovv"]("L,a,i,b,c")  = t1["bb"]("a,j") * l2["bbbb"]("L,j,i,b,c");

        // D_vvvv_bbbb += -1.00 r0 t1_bb(c,i) t1_bb(d,j) l2_bbbb(i,j,a,b)
        // flops: o0v4L2 += o1v4L1 o0v4L2
        //  mems: o0v4L2 += o0v4L1 o0v4L2
        D_vvvv_bbbb("L,R,a,b,c,d") -= tmps_["147_bbbb_Lvovv"]("L,c,j,a,b") * t1["bb"]("d,j") * r0("R");

        // flops: o1v3L1  = o2v3L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["148_bbbb_Lvovv"]("L,a,i,b,c")  = t1_1["bb"]("a,j") * l2_1["bbbb"]("L,j,i,b,c");

        // D_vvvv_bbbb += -1.00 P(c,d) r0 t1_bb(d,j) t1_1_bb(c,i) l2_1_bbbb(i,j,a,b)
        // flops: o0v4L2 += o1v4L1 o0v4L2
        //  mems: o0v4L2 += o0v4L1 o0v4L2
        tmps_["perm_bbbb_LLvvvv"]("L,R,a,b,c,d")  = tmps_["148_bbbb_Lvovv"]("L,c,j,a,b") * t1["bb"]("d,j") * r0("R");
        D_vvvv_bbbb("L,R,a,b,c,d") -= tmps_["perm_bbbb_LLvvvv"]("L,R,a,b,c,d");
        D_vvvv_bbbb("L,R,a,b,c,d") += tmps_["perm_bbbb_LLvvvv"]("L,R,a,b,d,c");
        tmps_["perm_bbbb_LLvvvv"].~TArrayD();

        // flops: o1v3L1  = o2v3L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["149_bbbb_Lvovv"]("L,a,i,b,c")  = t1["bb"]("a,j") * l2_1["bbbb"]("L,j,i,b,c");

        // D_vvvv_bbbb += -1.00 r0_1 t1_bb(c,i) t1_bb(d,j) l2_1_bbbb(i,j,a,b)
        // flops: o0v4L2 += o1v4L1 o0v4L2
        //  mems: o0v4L2 += o0v4L1 o0v4L2
        D_vvvv_bbbb("L,R,a,b,c,d") -= tmps_["149_bbbb_Lvovv"]("L,c,j,a,b") * t1["bb"]("d,j") * r0_1("R");

        // flops: o1v3L2  = o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        //  mems: o1v3L2  = o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        tmps_["150_bbbb_LLvovv"]("R,L,a,i,b,c")  = tmps_["147_bbbb_Lvovv"]("L,a,i,b,c") * r0("R");
        tmps_["150_bbbb_LLvovv"]("R,L,a,i,b,c") += tmps_["148_bbbb_Lvovv"]("L,a,i,b,c") * r0("R");
        tmps_["150_bbbb_LLvovv"]("R,L,a,i,b,c") += tmps_["149_bbbb_Lvovv"]("L,a,i,b,c") * r0_1("R");
        tmps_["149_bbbb_Lvovv"].~TArrayD();
        tmps_["148_bbbb_Lvovv"].~TArrayD();
        tmps_["147_bbbb_Lvovv"].~TArrayD();

        // D_vvov_bbbb += +1.00 r0_1 t1_bb(c,j) l2_1_bbbb(j,i,a,b)
        //             += +1.00 r0 t1_1_bb(c,j) l2_1_bbbb(j,i,a,b)
        //             += +1.00 r0 t1_bb(c,j) l2_bbbb(j,i,a,b)
        D_vvov_bbbb("L,R,a,b,i,c") += tmps_["150_bbbb_LLvovv"]("R,L,c,i,a,b");

        // D_vvvo_bbbb += -1.00 r0_1 t1_bb(c,j) l2_1_bbbb(j,i,a,b)
        //             += -1.00 r0 t1_1_bb(c,j) l2_1_bbbb(j,i,a,b)
        //             += -1.00 r0 t1_bb(c,j) l2_bbbb(j,i,a,b)
        D_vvvo_bbbb("L,R,a,b,c,i") -= tmps_["150_bbbb_LLvovv"]("R,L,c,i,a,b");
        tmps_["150_bbbb_LLvovv"].~TArrayD();

        // flops: o0v0L2  = o0v0L2
        //  mems: o0v0L2  = o0v0L2
        tmps_["151_LL"]("L,R")  = l0_1("L") * r0("R");

        // D_oovv_aaaa += +1.00 r0 t2_1_aaaa(b,a,i,j) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") += t2_1["aaaa"]("b,a,i,j") * tmps_["151_LL"]("L,R");

        // D_oovv_abab += -1.00 r0 t2_1_abab(a,b,i,j) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2_1["abab"]("a,b,i,j") * tmps_["151_LL"]("L,R");

        // D_oovv_bbbb += +1.00 r0 t2_1_bbbb(b,a,i,j) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += t2_1["bbbb"]("b,a,i,j") * tmps_["151_LL"]("L,R");
        tmps_["151_LL"].~TArrayD();

        // flops: o2v2L2  = o3v2L2 o3v2L2
        //  mems: o2v2L2  = o3v1L2 o2v2L2
        tmps_["152_bbbb_LLovov"]("L,R,i,a,j,b")  = l2["bbbb"]("L,k,i,a,c") * r1["bb"]("R,c,j") * t1["bb"]("b,k");

        // D_ovov_bbbb += +1.00 r1_bb(c,i) t1_bb(b,k) l2_bbbb(k,j,a,c)
        D_ovov_bbbb("L,R,i,a,j,b") += tmps_["152_bbbb_LLovov"]("L,R,j,a,i,b");

        // D_ovvo_bbbb += -1.00 r1_bb(c,i) t1_bb(b,k) l2_bbbb(k,j,a,c)
        D_ovvo_bbbb("L,R,i,a,b,j") -= tmps_["152_bbbb_LLovov"]("L,R,j,a,i,b");

        // flops: o2v2L2  = o3v2L2 o3v2L2
        //  mems: o2v2L2  = o3v1L2 o2v2L2
        tmps_["153_bbbb_LLovov"]("L,R,i,a,j,b")  = l2_1["bbbb"]("L,k,i,a,c") * r1_1["bb"]("R,c,j") * t1["bb"]("b,k");

        // D_ovov_bbbb += +1.00 r1_1_bb(c,i) t1_bb(b,k) l2_1_bbbb(k,j,a,c)
        D_ovov_bbbb("L,R,i,a,j,b") += tmps_["153_bbbb_LLovov"]("L,R,j,a,i,b");

        // D_ovvo_bbbb += -1.00 r1_1_bb(c,i) t1_bb(b,k) l2_1_bbbb(k,j,a,c)
        D_ovvo_bbbb("L,R,i,a,b,j") -= tmps_["153_bbbb_LLovov"]("L,R,j,a,i,b");

        // flops: o0v4L1  = o2v4L1
        //  mems: o0v4L1  = o0v4L1
        tmps_["154_bbbb_Lvvvv"]("L,a,b,c,d")  = l2_1["bbbb"]("L,i,j,a,b") * t2["bbbb"]("c,d,i,j");

        // D_oovv_bbbb += -0.50 P(i,j) r1_bb(d,i) t2_bbbb(b,a,k,l) t1_1_bb(c,j) l2_1_bbbb(k,l,c,d)
        // flops: o2v2L2 += o1v4L1 o2v3L2
        //  mems: o2v2L2 += o1v3L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["154_bbbb_Lvvvv"]("L,c,d,b,a") * t1_1["bb"]("c,j") * r1["bb"]("R,d,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_vvvv_bbbb += +0.50 r0_1 t2_bbbb(d,c,i,j) l2_1_bbbb(i,j,a,b)
        // flops: o0v4L2 += o0v4L2
        //  mems: o0v4L2 += o0v4L2
        D_vvvv_bbbb("L,R,a,b,c,d") += 0.50 * tmps_["154_bbbb_Lvvvv"]("L,a,b,d,c") * r0_1("R");

        // flops: o0v4L1  = o2v4L1
        //  mems: o0v4L1  = o0v4L1
        tmps_["155_bbbb_Lvvvv"]("L,a,b,c,d")  = l2_1["bbbb"]("L,i,j,a,b") * t2_1["bbbb"]("c,d,i,j");

        // D_vvvv_bbbb += +0.50 r0 t2_1_bbbb(d,c,i,j) l2_1_bbbb(i,j,a,b)
        // flops: o0v4L2 += o0v4L2
        //  mems: o0v4L2 += o0v4L2
        D_vvvv_bbbb("L,R,a,b,c,d") += 0.50 * tmps_["155_bbbb_Lvvvv"]("L,a,b,d,c") * r0("R");

        // flops: o0v4L1  = o2v4L1
        //  mems: o0v4L1  = o0v4L1
        tmps_["156_bbbb_Lvvvv"]("L,a,b,c,d")  = l2["bbbb"]("L,i,j,a,b") * t2["bbbb"]("c,d,i,j");

        // D_vvvv_bbbb += +0.50 r0 t2_bbbb(d,c,i,j) l2_bbbb(i,j,a,b)
        // flops: o0v4L2 += o0v4L2
        //  mems: o0v4L2 += o0v4L2
        D_vvvv_bbbb("L,R,a,b,c,d") += 0.50 * tmps_["156_bbbb_Lvvvv"]("L,a,b,d,c") * r0("R");

        // flops: o1v3L1  = o1v4L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["157_bbbb_Lvvvo"]("L,a,b,c,i")  = tmps_["154_bbbb_Lvvvv"]("L,a,d,b,c") * t1_1["bb"]("d,i");

        // D_oovv_bbbb += -0.50 P(i,j) r0 t2_bbbb(b,a,k,l) t1_bb(d,j) t1_1_bb(c,i) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["157_bbbb_Lvvvo"]("L,d,b,a,i") * t1["bb"]("d,j") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v3L1  = o1v4L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["158_bbbb_Lvvvo"]("L,a,b,c,i")  = tmps_["154_bbbb_Lvvvv"]("L,a,d,b,c") * t1["bb"]("d,i");

        // D_oovv_bbbb += -0.50 r0_1 t2_bbbb(b,a,k,l) t1_bb(c,i) t1_bb(d,j) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") -= 0.50 * tmps_["158_bbbb_Lvvvo"]("L,d,b,a,i") * t1["bb"]("d,j") * r0_1("R");

        // flops: o1v3L1  = o0v4L1 o1v4L1
        //  mems: o1v3L1  = o0v4L1 o1v3L1
        tmps_["159_bbbb_Lvvvo"]("L,a,b,c,i")  = (tmps_["155_bbbb_Lvvvv"]("L,a,d,b,c") + tmps_["156_bbbb_Lvvvv"]("L,a,d,b,c")) * t1["bb"]("d,i");

        // D_oovv_bbbb += -0.50 r0 t1_bb(c,i) t1_bb(d,j) t2_1_bbbb(b,a,k,l) l2_1_bbbb(k,l,d,c)
        //             += -0.50 r0 t2_bbbb(b,a,k,l) t1_bb(c,i) t1_bb(d,j) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") -= 0.50 * tmps_["159_bbbb_Lvvvo"]("L,d,b,a,i") * t1["bb"]("d,j") * r0("R");

        // flops: o1v3L2  = o2v3L2 o1v4L2 o1v4L2 o1v3L2 o1v4L2 o1v3L2 o1v3L2 o2v3L2 o2v3L1 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L1 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L1 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        //  mems: o1v3L2  = o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L1 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i")  = -2.00 * tmps_["71_bb_LLov"]("L,R,j,a") * t2["bbbb"]("b,c,j,i");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") += tmps_["154_bbbb_Lvvvv"]("L,a,d,b,c") * r1_1["bb"]("R,d,i");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") += tmps_["155_bbbb_Lvvvv"]("L,a,d,b,c") * r1["bb"]("R,d,i");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") += tmps_["156_bbbb_Lvvvv"]("L,a,d,b,c") * r1["bb"]("R,d,i");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") += 2.00 * tmps_["72_bb_LLov"]("L,R,j,a") * t2_1["bbbb"]("b,c,j,i");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l1_1["bb"]("L,j,a") * t2["bbbb"]("b,c,j,i") * r0_1("R");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * tmps_["152_bbbb_LLovov"]("L,R,k,a,i,c") * t1["bb"]("b,k");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l2_1["bbbb"]("L,j,k,a,d") * t1_1["bb"]("d,i") * t1["bb"]("c,j") * t1["bb"]("b,k") * r0("R");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l1["bb"]("L,j,a") * t2["bbbb"]("b,c,j,i") * r0("R");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * tmps_["70_bb_LLov"]("L,R,j,a") * t2["bbbb"]("b,c,j,i");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l2["bbbb"]("L,j,k,a,d") * t1["bb"]("d,i") * t1["bb"]("c,j") * t1["bb"]("b,k") * r0("R");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * tmps_["69_bb_LLov"]("L,R,j,a") * t2_1["bbbb"]("b,c,j,i");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") += 2.00 * tmps_["74_bb_LLov"]("L,R,j,a") * t2["bbbb"]("b,c,j,i");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l1_1["bb"]("L,j,a") * t2_1["bbbb"]("b,c,j,i") * r0("R");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * tmps_["153_bbbb_LLovov"]("L,R,k,a,i,c") * t1["bb"]("b,k");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l1_1["bb"]("L,j,a") * r2_1["bbbb"]("R,b,c,j,i");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l1["bb"]("L,j,a") * r2["bbbb"]("R,b,c,j,i");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l2_1["bbbb"]("L,j,k,a,d") * t1["bb"]("d,i") * t1["bb"]("c,j") * t1["bb"]("b,k") * r0_1("R");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") += 2.00 * tmps_["73_bb_LLov"]("L,R,j,a") * t2["bbbb"]("b,c,j,i");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") += r0("R") * tmps_["159_bbbb_Lvvvo"]("L,a,b,c,i");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") += r0_1("R") * tmps_["158_bbbb_Lvvvo"]("L,a,b,c,i");
        tmps_["160_bbbb_LLvvvo"]("R,L,a,b,c,i") += r0("R") * tmps_["157_bbbb_Lvvvo"]("L,a,b,c,i");
        tmps_["159_bbbb_Lvvvo"].~TArrayD();
        tmps_["158_bbbb_Lvvvo"].~TArrayD();
        tmps_["157_bbbb_Lvvvo"].~TArrayD();
        tmps_["156_bbbb_Lvvvv"].~TArrayD();
        tmps_["155_bbbb_Lvvvv"].~TArrayD();
        tmps_["154_bbbb_Lvvvv"].~TArrayD();
        tmps_["153_bbbb_LLovov"].~TArrayD();
        tmps_["152_bbbb_LLovov"].~TArrayD();

        // D_ovvv_bbbb += +0.50 r0_1 t2_bbbb(c,b,j,k) t1_bb(d,i) l2_1_bbbb(j,k,a,d)
        //             += +0.50 r1_1_bb(d,i) t2_bbbb(c,b,j,k) l2_1_bbbb(j,k,a,d)
        //             += +0.50 r1_bb(d,i) t2_1_bbbb(c,b,j,k) l2_1_bbbb(j,k,a,d)
        //             += +0.50 r1_bb(d,i) t2_bbbb(c,b,j,k) l2_bbbb(j,k,a,d)
        //             += -1.00 r1_1_aa(d,k) t2_bbbb(c,b,j,i) l2_1_abab(k,j,d,a)
        //             += +1.00 r1_bb(d,k) t2_1_bbbb(c,b,j,i) l2_1_bbbb(k,j,a,d)
        //             += -1.00 r0_1 t2_bbbb(c,b,j,i) l1_1_bb(j,a)
        //             += -1.00 r1_bb(d,i) t1_bb(b,j) t1_bb(c,k) l2_bbbb(j,k,a,d)
        //             += -1.00 r0 t1_bb(b,j) t1_bb(c,k) t1_1_bb(d,i) l2_1_bbbb(j,k,a,d)
        //             += -1.00 r0 t2_bbbb(c,b,j,i) l1_bb(j,a)
        //             += -1.00 r1_aa(d,k) t2_bbbb(c,b,j,i) l2_abab(k,j,d,a)
        //             += -1.00 r0 t1_bb(b,j) t1_bb(c,k) t1_bb(d,i) l2_bbbb(j,k,a,d)
        //             += -1.00 r1_aa(d,k) t2_1_bbbb(c,b,j,i) l2_1_abab(k,j,d,a)
        //             += +1.00 r1_bb(d,k) t2_bbbb(c,b,j,i) l2_bbbb(k,j,a,d)
        //             += -1.00 r0 t2_1_bbbb(c,b,j,i) l1_1_bb(j,a)
        //             += -1.00 r1_1_bb(d,i) t1_bb(b,j) t1_bb(c,k) l2_1_bbbb(j,k,a,d)
        //             += -1.00 r2_1_bbbb(c,b,j,i) l1_1_bb(j,a)
        //             += -1.00 r2_bbbb(c,b,j,i) l1_bb(j,a)
        //             += -1.00 r0_1 t1_bb(b,j) t1_bb(c,k) t1_bb(d,i) l2_1_bbbb(j,k,a,d)
        //             += +1.00 r1_1_bb(d,k) t2_bbbb(c,b,j,i) l2_1_bbbb(k,j,a,d)
        //             += +0.50 r0 t1_bb(d,i) t2_1_bbbb(c,b,j,k) l2_1_bbbb(j,k,a,d)
        //             += +0.50 r0 t2_bbbb(c,b,j,k) t1_bb(d,i) l2_bbbb(j,k,a,d)
        //             += +0.50 r0 t2_bbbb(c,b,j,k) t1_1_bb(d,i) l2_1_bbbb(j,k,a,d)
        D_ovvv_bbbb("L,R,i,a,b,c") += 0.50 * tmps_["160_bbbb_LLvvvo"]("R,L,a,c,b,i");

        // D_vovv_bbbb += -0.50 r0_1 t2_bbbb(c,b,j,k) t1_bb(d,i) l2_1_bbbb(j,k,a,d)
        //             += -0.50 r1_1_bb(d,i) t2_bbbb(c,b,j,k) l2_1_bbbb(j,k,a,d)
        //             += -0.50 r1_bb(d,i) t2_1_bbbb(c,b,j,k) l2_1_bbbb(j,k,a,d)
        //             += -0.50 r1_bb(d,i) t2_bbbb(c,b,j,k) l2_bbbb(j,k,a,d)
        //             += +1.00 r1_1_aa(d,k) t2_bbbb(c,b,j,i) l2_1_abab(k,j,d,a)
        //             += -1.00 r1_bb(d,k) t2_1_bbbb(c,b,j,i) l2_1_bbbb(k,j,a,d)
        //             += +1.00 r0_1 t2_bbbb(c,b,j,i) l1_1_bb(j,a)
        //             += +1.00 r1_bb(d,i) t1_bb(b,j) t1_bb(c,k) l2_bbbb(j,k,a,d)
        //             += +1.00 r0 t1_bb(b,j) t1_bb(c,k) t1_1_bb(d,i) l2_1_bbbb(j,k,a,d)
        //             += +1.00 r0 t2_bbbb(c,b,j,i) l1_bb(j,a)
        //             += +1.00 r1_aa(d,k) t2_bbbb(c,b,j,i) l2_abab(k,j,d,a)
        //             += +1.00 r0 t1_bb(b,j) t1_bb(c,k) t1_bb(d,i) l2_bbbb(j,k,a,d)
        //             += +1.00 r1_aa(d,k) t2_1_bbbb(c,b,j,i) l2_1_abab(k,j,d,a)
        //             += -1.00 r1_bb(d,k) t2_bbbb(c,b,j,i) l2_bbbb(k,j,a,d)
        //             += +1.00 r0 t2_1_bbbb(c,b,j,i) l1_1_bb(j,a)
        //             += +1.00 r1_1_bb(d,i) t1_bb(b,j) t1_bb(c,k) l2_1_bbbb(j,k,a,d)
        //             += +1.00 r2_1_bbbb(c,b,j,i) l1_1_bb(j,a)
        //             += +1.00 r2_bbbb(c,b,j,i) l1_bb(j,a)
        //             += +1.00 r0_1 t1_bb(b,j) t1_bb(c,k) t1_bb(d,i) l2_1_bbbb(j,k,a,d)
        //             += -1.00 r1_1_bb(d,k) t2_bbbb(c,b,j,i) l2_1_bbbb(k,j,a,d)
        //             += -0.50 r0 t1_bb(d,i) t2_1_bbbb(c,b,j,k) l2_1_bbbb(j,k,a,d)
        //             += -0.50 r0 t2_bbbb(c,b,j,k) t1_bb(d,i) l2_bbbb(j,k,a,d)
        //             += -0.50 r0 t2_bbbb(c,b,j,k) t1_1_bb(d,i) l2_1_bbbb(j,k,a,d)
        D_vovv_bbbb("L,R,a,i,b,c") -= 0.50 * tmps_["160_bbbb_LLvvvo"]("R,L,a,c,b,i");
        tmps_["160_bbbb_LLvvvo"].~TArrayD();

        // flops: o4v0L1  = o4v2L1
        //  mems: o4v0L1  = o4v0L1
        tmps_["161_bbbb_Loooo"]("L,i,j,k,l")  = t2["bbbb"]("a,b,i,j") * l2_1["bbbb"]("L,k,l,a,b");

        // D_oooo_bbbb += +0.50 r0_1 t2_bbbb(b,a,i,j) l2_1_bbbb(l,k,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_bbbb("L,R,i,j,k,l") += 0.50 * tmps_["161_bbbb_Loooo"]("L,i,j,l,k") * r0_1("R");

        // D_oovv_bbbb += +0.25 r2_1_bbbb(b,a,l,k) t2_bbbb(d,c,i,j) l2_1_bbbb(l,k,d,c)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += 0.25 * tmps_["161_bbbb_Loooo"]("L,i,j,l,k") * r2_1["bbbb"]("R,b,a,l,k");

        // D_oovv_bbbb += +0.25 r0_1 t2_bbbb(b,a,k,l) t2_bbbb(d,c,i,j) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o4v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += 0.25 * t2["bbbb"]("b,a,k,l") * tmps_["161_bbbb_Loooo"]("L,i,j,k,l") * r0_1("R");

        // D_oovv_bbbb += +0.25 r0 t2_bbbb(d,c,i,j) t2_1_bbbb(b,a,k,l) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o4v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += 0.25 * t2_1["bbbb"]("b,a,k,l") * tmps_["161_bbbb_Loooo"]("L,i,j,k,l") * r0("R");

        // D_oovv_bbbb += -0.50 P(a,b) r1_bb(a,l) t2_bbbb(d,c,i,j) t1_1_bb(b,k) l2_1_bbbb(l,k,d,c)
        // flops: o2v2L2 += o4v1L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t1_1["bb"]("b,k") * tmps_["161_bbbb_Loooo"]("L,i,j,l,k") * r1["bb"]("R,a,l");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["162_bbbb_Loovo"]("L,i,j,a,k")  = l2["bbbb"]("L,i,j,b,a") * t1["bb"]("b,k");

        // D_oooo_bbbb += -1.00 P(i,j) r1_bb(b,i) t1_bb(a,j) l2_bbbb(l,k,a,b)
        // flops: o4v0L2 += o4v1L2
        //  mems: o4v0L2 += o4v0L2
        tmps_["perm_bbbb_LLoooo"]("L,R,i,j,k,l")  = tmps_["162_bbbb_Loovo"]("L,l,k,b,j") * r1["bb"]("R,b,i");
        D_oooo_bbbb("L,R,i,j,k,l") -= tmps_["perm_bbbb_LLoooo"]("L,R,i,j,k,l");
        D_oooo_bbbb("L,R,i,j,k,l") += tmps_["perm_bbbb_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_bbbb_LLoooo"].~TArrayD();

        // D_oooo_bbbb += -1.00 r0 t1_bb(a,i) t1_bb(b,j) l2_bbbb(l,k,b,a)
        // flops: o4v0L2 += o4v1L1 o4v0L2
        //  mems: o4v0L2 += o4v0L1 o4v0L2
        D_oooo_bbbb("L,R,i,j,k,l") -= t1["bb"]("a,i") * tmps_["162_bbbb_Loovo"]("L,l,k,a,j") * r0("R");

        // D_oovv_abab += +1.00 r1_bb(b,l) t2_abab(a,c,i,k) t1_bb(d,j) l2_bbbb(l,k,d,c)
        // flops: o2v2L2 += o4v2L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,c,i,k") * tmps_["162_bbbb_Loovo"]("L,l,k,c,j") * r1["bb"]("R,b,l");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r1_bb(a,l) t2_bbbb(b,c,k,i) t1_bb(d,j) l2_bbbb(l,k,d,c)
        // flops: o2v2L2 += o4v2L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["bbbb"]("b,c,k,i") * tmps_["162_bbbb_Loovo"]("L,l,k,c,j") * r1["bb"]("R,a,l");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["163_bbbb_Lvooo"]("L,a,i,j,k")  = t2["bbbb"]("a,b,i,j") * l1_1["bb"]("L,k,b");

        // D_oovv_bbbb += -1.00 P(a,b) r1_1_bb(a,k) t2_bbbb(b,c,i,j) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["163_bbbb_Lvooo"]("L,b,i,j,k") * r1_1["bb"]("R,a,k");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(a,b) r0_1 t2_bbbb(a,c,i,j) t1_bb(b,k) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["163_bbbb_Lvooo"]("L,a,i,j,k") * t1["bb"]("b,k") * r0_1("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += -1.00 P(a,b) r0 t2_bbbb(b,c,i,j) t1_1_bb(a,k) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1_1["bb"]("a,k") * tmps_["163_bbbb_Lvooo"]("L,b,i,j,k") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o4v0L1  = o4v2L1 o4v2L1 o4v0L1
        //  mems: o4v0L1  = o4v0L1 o4v0L1 o4v0L1
        tmps_["164_bbbb_Loooo"]("L,i,j,k,l")  = t2["bbbb"]("a,b,i,j") * l2["bbbb"]("L,k,l,a,b");
        tmps_["164_bbbb_Loooo"]("L,i,j,k,l") += t2_1["bbbb"]("a,b,i,j") * l2_1["bbbb"]("L,k,l,a,b");

        // D_oooo_bbbb += +0.50 r0 t2_bbbb(b,a,i,j) l2_bbbb(l,k,b,a)
        //             += +0.50 r0 t2_1_bbbb(b,a,i,j) l2_1_bbbb(l,k,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_bbbb("L,R,i,j,k,l") += 0.50 * tmps_["164_bbbb_Loooo"]("L,i,j,l,k") * r0("R");

        // D_oovv_bbbb += +0.25 r2_bbbb(b,a,l,k) t2_bbbb(d,c,i,j) l2_bbbb(l,k,d,c)
        //             += +0.25 r2_bbbb(b,a,l,k) t2_1_bbbb(d,c,i,j) l2_1_bbbb(l,k,d,c)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += 0.25 * tmps_["164_bbbb_Loooo"]("L,i,j,l,k") * r2["bbbb"]("R,b,a,l,k");

        // D_oovv_bbbb += +0.25 r0 t2_bbbb(b,a,k,l) t2_bbbb(d,c,i,j) l2_bbbb(k,l,d,c)
        //             += +0.25 r0 t2_bbbb(b,a,k,l) t2_1_bbbb(d,c,i,j) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o4v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += 0.25 * t2["bbbb"]("b,a,k,l") * tmps_["164_bbbb_Loooo"]("L,i,j,k,l") * r0("R");

        // flops: o3v1L1  = o3v2L1 o3v2L1 o3v1L1
        //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1
        tmps_["165_bbbb_Lvooo"]("L,a,i,j,k")  = t2["bbbb"]("a,b,i,j") * l1["bb"]("L,k,b");
        tmps_["165_bbbb_Lvooo"]("L,a,i,j,k") += t2_1["bbbb"]("a,b,i,j") * l1_1["bb"]("L,k,b");

        // D_oovv_bbbb += -1.00 P(a,b) r1_bb(a,k) t2_bbbb(b,c,i,j) l1_bb(k,c)
        //             += -1.00 P(a,b) r1_bb(a,k) t2_1_bbbb(b,c,i,j) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["165_bbbb_Lvooo"]("L,b,i,j,k") * r1["bb"]("R,a,k");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(a,b) r0 t2_bbbb(a,c,i,j) t1_bb(b,k) l1_bb(k,c)
        //             += +1.00 P(a,b) r0 t1_bb(b,k) t2_1_bbbb(a,c,i,j) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["165_bbbb_Lvooo"]("L,a,i,j,k") * t1["bb"]("b,k") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["166_bbbb_Loovo"]("L,i,j,a,k")  = l2_1["bbbb"]("L,i,j,b,a") * t1["bb"]("b,k");

        // D_oooo_bbbb += -1.00 P(i,j) r1_1_bb(b,i) t1_bb(a,j) l2_1_bbbb(l,k,a,b)
        // flops: o4v0L2 += o4v1L2
        //  mems: o4v0L2 += o4v0L2
        tmps_["perm_bbbb_LLoooo"]("L,R,i,j,k,l")  = tmps_["166_bbbb_Loovo"]("L,l,k,b,j") * r1_1["bb"]("R,b,i");
        D_oooo_bbbb("L,R,i,j,k,l") -= tmps_["perm_bbbb_LLoooo"]("L,R,i,j,k,l");
        D_oooo_bbbb("L,R,i,j,k,l") += tmps_["perm_bbbb_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_bbbb_LLoooo"].~TArrayD();

        // D_oooo_bbbb += -1.00 P(i,j) r0 t1_bb(b,j) t1_1_bb(a,i) l2_1_bbbb(l,k,b,a)
        // flops: o4v0L2 += o4v1L1 o4v0L2
        //  mems: o4v0L2 += o4v0L1 o4v0L2
        tmps_["perm_bbbb_LLoooo"]("L,R,i,j,k,l")  = t1_1["bb"]("a,i") * tmps_["166_bbbb_Loovo"]("L,l,k,a,j") * r0("R");
        D_oooo_bbbb("L,R,i,j,k,l") -= tmps_["perm_bbbb_LLoooo"]("L,R,i,j,k,l");
        D_oooo_bbbb("L,R,i,j,k,l") += tmps_["perm_bbbb_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_bbbb_LLoooo"].~TArrayD();

        // D_oovv_abab += +1.00 r1_bb(b,l) t1_bb(d,j) t2_1_abab(a,c,i,k) l2_1_bbbb(l,k,d,c)
        // flops: o2v2L2 += o4v2L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2_1["abab"]("a,c,i,k") * tmps_["166_bbbb_Loovo"]("L,l,k,c,j") * r1["bb"]("R,b,l");

        // D_oovv_abab += -1.00 r2_abab(a,d,i,l) t1_bb(c,j) t1_1_bb(b,k) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o4v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["166_bbbb_Loovo"]("L,l,k,d,j") * r2["abab"]("R,a,d,i,l") * t1_1["bb"]("b,k");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r1_bb(a,l) t1_bb(d,j) t2_1_bbbb(b,c,k,i) l2_1_bbbb(l,k,d,c)
        // flops: o2v2L2 += o4v2L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2_1["bbbb"]("b,c,k,i") * tmps_["166_bbbb_Loovo"]("L,l,k,c,j") * r1["bb"]("R,a,l");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o4v0L1  = o4v1L1
        //  mems: o4v0L1  = o4v0L1
        tmps_["167_bbbb_Loooo"]("L,i,j,k,l")  = t1["bb"]("a,i") * tmps_["166_bbbb_Loovo"]("L,j,k,a,l");

        // D_oooo_bbbb += -1.00 r0_1 t1_bb(a,i) t1_bb(b,j) l2_1_bbbb(l,k,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_bbbb("L,R,i,j,k,l") -= tmps_["167_bbbb_Loooo"]("L,i,l,k,j") * r0_1("R");

        // D_oovv_bbbb += +1.00 P(a,b) r1_bb(a,l) t1_bb(c,i) t1_bb(d,j) t1_1_bb(b,k) l2_1_bbbb(l,k,d,c)
        // flops: o2v2L2 += o4v1L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1_1["bb"]("b,k") * tmps_["167_bbbb_Loooo"]("L,i,l,k,j") * r1["bb"]("R,a,l");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o4v1L1 o4v0L1 o4v1L1
        //  mems: o3v1L1  = o4v0L1 o4v0L1 o3v1L1
        tmps_["168_bbbb_Lvooo"]("L,a,i,j,k")  = t1["bb"]("a,l") * (t1["bb"]("b,i") * tmps_["162_bbbb_Loovo"]("L,l,j,b,k") + -0.50 * tmps_["164_bbbb_Loooo"]("L,i,k,l,j"));
        tmps_["162_bbbb_Loovo"].~TArrayD();

        // D_oovv_bbbb += +1.00 r0 t1_bb(a,k) t1_bb(b,l) t1_bb(c,i) t1_bb(d,j) l2_bbbb(k,l,d,c)
        //             += -0.50 r0 t1_bb(a,k) t1_bb(b,l) t2_bbbb(d,c,i,j) l2_bbbb(k,l,d,c)
        //             += -0.50 r0 t1_bb(a,k) t1_bb(b,l) t2_1_bbbb(d,c,i,j) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["168_bbbb_Lvooo"]("L,a,i,l,j") * t1["bb"]("b,l") * r0("R");

        // flops: o3v1L1  = o4v0L1 o4v1L1
        //  mems: o3v1L1  = o4v0L1 o3v1L1
        tmps_["169_bbbb_Looov"]("L,i,j,k,a")  = (tmps_["161_bbbb_Loooo"]("L,i,j,l,k") + -2.00 * tmps_["167_bbbb_Loooo"]("L,i,l,k,j")) * t1_1["bb"]("a,l");

        // D_oovv_bbbb += -0.50 P(a,b) r0 t1_bb(b,l) t2_bbbb(d,c,i,j) t1_1_bb(a,k) l2_1_bbbb(k,l,d,c)
        //             += +1.00 P(a,b) r0 t1_bb(b,l) t1_bb(c,i) t1_bb(d,j) t1_1_bb(a,k) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t1["bb"]("b,l") * tmps_["169_bbbb_Looov"]("L,i,j,l,a") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o4v0L1 o4v1L1
        //  mems: o3v1L1  = o4v0L1 o3v1L1
        tmps_["170_bbbb_Looov"]("L,i,j,k,a")  = (tmps_["167_bbbb_Loooo"]("L,i,l,j,k") + -0.50 * tmps_["161_bbbb_Loooo"]("L,i,k,l,j")) * t1["bb"]("a,l");
        tmps_["167_bbbb_Loooo"].~TArrayD();

        // D_oovv_bbbb += +1.00 r0_1 t1_bb(a,k) t1_bb(b,l) t1_bb(c,i) t1_bb(d,j) l2_1_bbbb(k,l,d,c)
        //             += -0.50 r0_1 t1_bb(a,k) t1_bb(b,l) t2_bbbb(d,c,i,j) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += t1["bb"]("b,l") * tmps_["170_bbbb_Looov"]("L,i,l,j,a") * r0_1("R");

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["171_bbbb_Looov"]("L,i,j,k,a")  = t1["bb"]("b,i") * l2_1["bbbb"]("L,j,k,a,b");

        // flops: o3v1L2  = o3v1L2 o4v1L2 o1v1L2 o3v2L2 o3v2L2 o4v1L1 o4v1L2 o3v1L2 o3v2L1 o4v1L1 o4v1L2 o3v1L2 o4v1L2 o3v1L2 o3v2L2 o3v1L2 o3v2L2 o3v1L2 o3v2L2 o3v1L2 o3v2L2 o3v1L2 o3v2L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        //  mems: o3v1L2  = o3v1L2 o3v1L2 o1v1L2 o3v1L2 o3v1L2 o4v0L1 o3v1L2 o3v1L2 o3v1L1 o4v0L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        tmps_["172_bbbb_LLvooo"]("R,L,a,i,j,k")  = tmps_["163_bbbb_Lvooo"]("L,a,i,j,k") * r0_1("R");
        tmps_["172_bbbb_LLvooo"]("R,L,a,i,j,k") -= 0.50 * r1["bb"]("R,a,l") * tmps_["164_bbbb_Loooo"]("L,i,j,l,k");
        tmps_["172_bbbb_LLvooo"]("R,L,a,i,j,k") += t2["bbbb"]("a,b,i,j") * (tmps_["71_bb_LLov"]("L,R,k,b") + -1.00 * tmps_["73_bb_LLov"]("L,R,k,b"));
        tmps_["172_bbbb_LLvooo"]("R,L,a,i,j,k") -= tmps_["74_bb_LLov"]("L,R,k,b") * t2["bbbb"]("a,b,i,j");
        tmps_["172_bbbb_LLvooo"]("R,L,a,i,j,k") += tmps_["171_bbbb_Looov"]("L,i,l,k,c") * t1["bb"]("c,j") * r1_1["bb"]("R,a,l");
        tmps_["172_bbbb_LLvooo"]("R,L,a,i,j,k") += t1["bb"]("b,i") * l2["bbbb"]("L,l,k,c,b") * t1["bb"]("c,j") * r1["bb"]("R,a,l");
        tmps_["172_bbbb_LLvooo"]("R,L,a,i,j,k") -= 0.50 * tmps_["161_bbbb_Loooo"]("L,i,j,l,k") * r1_1["bb"]("R,a,l");
        tmps_["172_bbbb_LLvooo"]("R,L,a,i,j,k") += tmps_["69_bb_LLov"]("L,R,k,b") * t2_1["bbbb"]("a,b,i,j");
        tmps_["172_bbbb_LLvooo"]("R,L,a,i,j,k") += tmps_["70_bb_LLov"]("L,R,k,b") * t2["bbbb"]("a,b,i,j");
        tmps_["172_bbbb_LLvooo"]("R,L,a,i,j,k") += l1["bb"]("L,k,b") * r2["bbbb"]("R,a,b,i,j");
        tmps_["172_bbbb_LLvooo"]("R,L,a,i,j,k") += l1_1["bb"]("L,k,b") * r2_1["bbbb"]("R,a,b,i,j");
        tmps_["172_bbbb_LLvooo"]("R,L,a,i,j,k") -= tmps_["72_bb_LLov"]("L,R,k,b") * t2_1["bbbb"]("a,b,i,j");
        tmps_["172_bbbb_LLvooo"]("R,L,a,i,j,k") += tmps_["165_bbbb_Lvooo"]("L,a,i,j,k") * r0("R");
        tmps_["172_bbbb_LLvooo"]("R,L,a,i,j,k") += r0("R") * tmps_["168_bbbb_Lvooo"]("L,a,i,k,j");
        tmps_["172_bbbb_LLvooo"]("R,L,a,i,j,k") -= 0.50 * r0("R") * tmps_["169_bbbb_Looov"]("L,i,j,k,a");
        tmps_["172_bbbb_LLvooo"]("R,L,a,i,j,k") += r0_1("R") * tmps_["170_bbbb_Looov"]("L,i,k,j,a");
        tmps_["170_bbbb_Looov"].~TArrayD();
        tmps_["169_bbbb_Looov"].~TArrayD();
        tmps_["168_bbbb_Lvooo"].~TArrayD();
        tmps_["165_bbbb_Lvooo"].~TArrayD();
        tmps_["164_bbbb_Loooo"].~TArrayD();
        tmps_["163_bbbb_Lvooo"].~TArrayD();
        tmps_["161_bbbb_Loooo"].~TArrayD();

        // D_ooov_bbbb += -1.00 r0_1 t2_bbbb(a,b,i,j) l1_1_bb(k,b)
        //             += -1.00 r0_1 t1_bb(a,l) t1_bb(b,i) t1_bb(c,j) l2_1_bbbb(l,k,c,b)
        //             += +0.50 r0_1 t1_bb(a,l) t2_bbbb(c,b,i,j) l2_1_bbbb(l,k,c,b)
        //             += +0.50 r0 t2_bbbb(c,b,i,j) t1_1_bb(a,l) l2_1_bbbb(l,k,c,b)
        //             += -1.00 r0 t1_bb(b,i) t1_bb(c,j) t1_1_bb(a,l) l2_1_bbbb(l,k,c,b)
        //             += -1.00 r0 t1_bb(a,l) t1_bb(b,i) t1_bb(c,j) l2_bbbb(l,k,c,b)
        //             += +0.50 r0 t1_bb(a,l) t2_bbbb(c,b,i,j) l2_bbbb(l,k,c,b)
        //             += +0.50 r0 t1_bb(a,l) t2_1_bbbb(c,b,i,j) l2_1_bbbb(l,k,c,b)
        //             += +0.50 r1_bb(a,l) t2_bbbb(c,b,i,j) l2_bbbb(l,k,c,b)
        //             += +0.50 r1_bb(a,l) t2_1_bbbb(c,b,i,j) l2_1_bbbb(l,k,c,b)
        //             += -1.00 r1_1_aa(c,l) t2_bbbb(a,b,i,j) l2_1_abab(l,k,c,b)
        //             += +1.00 r1_1_bb(c,l) t2_bbbb(a,b,i,j) l2_1_bbbb(l,k,b,c)
        //             += +1.00 r1_bb(c,l) t2_bbbb(a,b,i,j) l2_bbbb(l,k,b,c)
        //             += -1.00 r1_1_bb(a,l) t1_bb(b,i) t1_bb(c,j) l2_1_bbbb(l,k,c,b)
        //             += -1.00 r1_bb(a,l) t1_bb(b,i) t1_bb(c,j) l2_bbbb(l,k,c,b)
        //             += +0.50 r1_1_bb(a,l) t2_bbbb(c,b,i,j) l2_1_bbbb(l,k,c,b)
        //             += -1.00 r1_aa(c,l) t2_1_bbbb(a,b,i,j) l2_1_abab(l,k,c,b)
        //             += -1.00 r1_aa(c,l) t2_bbbb(a,b,i,j) l2_abab(l,k,c,b)
        //             += -1.00 r2_bbbb(a,b,i,j) l1_bb(k,b)
        //             += -1.00 r2_1_bbbb(a,b,i,j) l1_1_bb(k,b)
        //             += +1.00 r1_bb(c,l) t2_1_bbbb(a,b,i,j) l2_1_bbbb(l,k,b,c)
        //             += -1.00 r0 t2_bbbb(a,b,i,j) l1_bb(k,b)
        //             += -1.00 r0 t2_1_bbbb(a,b,i,j) l1_1_bb(k,b)
        D_ooov_bbbb("L,R,i,j,k,a") -= tmps_["172_bbbb_LLvooo"]("R,L,a,i,j,k");

        // D_oovo_bbbb += +1.00 r0_1 t2_bbbb(a,b,i,j) l1_1_bb(k,b)
        //             += +1.00 r0_1 t1_bb(a,l) t1_bb(b,i) t1_bb(c,j) l2_1_bbbb(l,k,c,b)
        //             += -0.50 r0_1 t1_bb(a,l) t2_bbbb(c,b,i,j) l2_1_bbbb(l,k,c,b)
        //             += -0.50 r0 t2_bbbb(c,b,i,j) t1_1_bb(a,l) l2_1_bbbb(l,k,c,b)
        //             += +1.00 r0 t1_bb(b,i) t1_bb(c,j) t1_1_bb(a,l) l2_1_bbbb(l,k,c,b)
        //             += +1.00 r0 t1_bb(a,l) t1_bb(b,i) t1_bb(c,j) l2_bbbb(l,k,c,b)
        //             += -0.50 r0 t1_bb(a,l) t2_bbbb(c,b,i,j) l2_bbbb(l,k,c,b)
        //             += -0.50 r0 t1_bb(a,l) t2_1_bbbb(c,b,i,j) l2_1_bbbb(l,k,c,b)
        //             += -0.50 r1_bb(a,l) t2_bbbb(c,b,i,j) l2_bbbb(l,k,c,b)
        //             += -0.50 r1_bb(a,l) t2_1_bbbb(c,b,i,j) l2_1_bbbb(l,k,c,b)
        //             += +1.00 r1_1_aa(c,l) t2_bbbb(a,b,i,j) l2_1_abab(l,k,c,b)
        //             += -1.00 r1_1_bb(c,l) t2_bbbb(a,b,i,j) l2_1_bbbb(l,k,b,c)
        //             += -1.00 r1_bb(c,l) t2_bbbb(a,b,i,j) l2_bbbb(l,k,b,c)
        //             += +1.00 r1_1_bb(a,l) t1_bb(b,i) t1_bb(c,j) l2_1_bbbb(l,k,c,b)
        //             += +1.00 r1_bb(a,l) t1_bb(b,i) t1_bb(c,j) l2_bbbb(l,k,c,b)
        //             += -0.50 r1_1_bb(a,l) t2_bbbb(c,b,i,j) l2_1_bbbb(l,k,c,b)
        //             += +1.00 r1_aa(c,l) t2_1_bbbb(a,b,i,j) l2_1_abab(l,k,c,b)
        //             += +1.00 r1_aa(c,l) t2_bbbb(a,b,i,j) l2_abab(l,k,c,b)
        //             += +1.00 r2_bbbb(a,b,i,j) l1_bb(k,b)
        //             += +1.00 r2_1_bbbb(a,b,i,j) l1_1_bb(k,b)
        //             += -1.00 r1_bb(c,l) t2_1_bbbb(a,b,i,j) l2_1_bbbb(l,k,b,c)
        //             += +1.00 r0 t2_bbbb(a,b,i,j) l1_bb(k,b)
        //             += +1.00 r0 t2_1_bbbb(a,b,i,j) l1_1_bb(k,b)
        D_oovo_bbbb("L,R,i,j,a,k") += tmps_["172_bbbb_LLvooo"]("R,L,a,i,j,k");
        tmps_["172_bbbb_LLvooo"].~TArrayD();

        // flops: o4v0L1  = o4v2L1
        //  mems: o4v0L1  = o4v0L1
        tmps_["173_aaaa_Loooo"]("L,i,j,k,l")  = t2["aaaa"]("a,b,i,j") * l2_1["aaaa"]("L,k,l,a,b");

        // D_oooo_aaaa += +0.50 r0_1 t2_aaaa(b,a,i,j) l2_1_aaaa(l,k,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_aaaa("L,R,i,j,k,l") += 0.50 * tmps_["173_aaaa_Loooo"]("L,i,j,l,k") * r0_1("R");

        // D_oovv_aaaa += +0.25 r2_1_aaaa(b,a,l,k) t2_aaaa(d,c,i,j) l2_1_aaaa(l,k,d,c)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") += 0.25 * tmps_["173_aaaa_Loooo"]("L,i,j,l,k") * r2_1["aaaa"]("R,b,a,l,k");

        // D_oovv_aaaa += +0.25 r0_1 t2_aaaa(b,a,k,l) t2_aaaa(d,c,i,j) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o4v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") += 0.25 * t2["aaaa"]("b,a,k,l") * tmps_["173_aaaa_Loooo"]("L,i,j,k,l") * r0_1("R");

        // D_oovv_aaaa += +0.25 r0 t2_aaaa(d,c,i,j) t2_1_aaaa(b,a,k,l) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o4v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") += 0.25 * t2_1["aaaa"]("b,a,k,l") * tmps_["173_aaaa_Loooo"]("L,i,j,k,l") * r0("R");

        // D_oovv_aaaa += -0.50 P(a,b) r1_aa(a,l) t2_aaaa(d,c,i,j) t1_1_aa(b,k) l2_1_aaaa(l,k,d,c)
        // flops: o2v2L2 += o4v1L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t1_1["aa"]("b,k") * tmps_["173_aaaa_Loooo"]("L,i,j,l,k") * r1["aa"]("R,a,l");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["174_aaaa_Lvooo"]("L,a,i,j,k")  = t2["aaaa"]("a,b,i,j") * l1_1["aa"]("L,k,b");

        // D_oovv_aaaa += -1.00 P(a,b) r1_1_aa(a,k) t2_aaaa(b,c,i,j) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["174_aaaa_Lvooo"]("L,b,i,j,k") * r1_1["aa"]("R,a,k");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(a,b) r0_1 t2_aaaa(a,c,i,j) t1_aa(b,k) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["174_aaaa_Lvooo"]("L,a,i,j,k") * t1["aa"]("b,k") * r0_1("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += -1.00 P(a,b) r0 t2_aaaa(b,c,i,j) t1_1_aa(a,k) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1_1["aa"]("a,k") * tmps_["174_aaaa_Lvooo"]("L,b,i,j,k") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o4v0L1  = o4v2L1 o4v2L1 o4v0L1
        //  mems: o4v0L1  = o4v0L1 o4v0L1 o4v0L1
        tmps_["175_aaaa_Loooo"]("L,i,j,k,l")  = t2["aaaa"]("a,b,i,j") * l2["aaaa"]("L,k,l,a,b");
        tmps_["175_aaaa_Loooo"]("L,i,j,k,l") += t2_1["aaaa"]("a,b,i,j") * l2_1["aaaa"]("L,k,l,a,b");

        // D_oooo_aaaa += +0.50 r0 t2_aaaa(b,a,i,j) l2_aaaa(l,k,b,a)
        //             += +0.50 r0 t2_1_aaaa(b,a,i,j) l2_1_aaaa(l,k,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_aaaa("L,R,i,j,k,l") += 0.50 * tmps_["175_aaaa_Loooo"]("L,i,j,l,k") * r0("R");

        // D_oovv_aaaa += +0.25 r2_aaaa(b,a,l,k) t2_aaaa(d,c,i,j) l2_aaaa(l,k,d,c)
        //             += +0.25 r2_aaaa(b,a,l,k) t2_1_aaaa(d,c,i,j) l2_1_aaaa(l,k,d,c)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") += 0.25 * tmps_["175_aaaa_Loooo"]("L,i,j,l,k") * r2["aaaa"]("R,b,a,l,k");

        // D_oovv_aaaa += +0.25 r0 t2_aaaa(b,a,k,l) t2_aaaa(d,c,i,j) l2_aaaa(k,l,d,c)
        //             += +0.25 r0 t2_aaaa(b,a,k,l) t2_1_aaaa(d,c,i,j) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o4v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") += 0.25 * t2["aaaa"]("b,a,k,l") * tmps_["175_aaaa_Loooo"]("L,i,j,k,l") * r0("R");

        // flops: o3v1L1  = o3v2L1 o3v2L1 o3v1L1
        //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1
        tmps_["176_aaaa_Lvooo"]("L,a,i,j,k")  = t2["aaaa"]("a,b,i,j") * l1["aa"]("L,k,b");
        tmps_["176_aaaa_Lvooo"]("L,a,i,j,k") += t2_1["aaaa"]("a,b,i,j") * l1_1["aa"]("L,k,b");

        // D_oovv_aaaa += -1.00 P(a,b) r1_aa(a,k) t2_aaaa(b,c,i,j) l1_aa(k,c)
        //             += -1.00 P(a,b) r1_aa(a,k) t2_1_aaaa(b,c,i,j) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["176_aaaa_Lvooo"]("L,b,i,j,k") * r1["aa"]("R,a,k");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(a,b) r0 t2_aaaa(a,c,i,j) t1_aa(b,k) l1_aa(k,c)
        //             += +1.00 P(a,b) r0 t1_aa(b,k) t2_1_aaaa(a,c,i,j) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["176_aaaa_Lvooo"]("L,a,i,j,k") * t1["aa"]("b,k") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o4v0L1  = o4v1L1
        //  mems: o4v0L1  = o4v0L1
        tmps_["177_aaaa_Loooo"]("L,i,j,k,l")  = t1["aa"]("a,i") * tmps_["22_aaaa_Looov"]("L,j,k,l,a");

        // D_oooo_aaaa += -1.00 r0_1 t1_aa(a,i) t1_aa(b,j) l2_1_aaaa(l,k,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_aaaa("L,R,i,j,k,l") -= tmps_["177_aaaa_Loooo"]("L,i,j,l,k") * r0_1("R");

        // D_oovv_aaaa += +1.00 P(a,b) r1_aa(a,l) t1_aa(c,i) t1_aa(d,j) t1_1_aa(b,k) l2_1_aaaa(l,k,d,c)
        // flops: o2v2L2 += o4v1L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1_1["aa"]("b,k") * tmps_["177_aaaa_Loooo"]("L,i,j,l,k") * r1["aa"]("R,a,l");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o4v1L1 o4v0L1 o4v1L1
        //  mems: o3v1L1  = o4v0L1 o4v0L1 o3v1L1
        tmps_["178_aaaa_Looov"]("L,i,j,k,a")  = (tmps_["175_aaaa_Loooo"]("L,i,j,l,k") + -2.00 * t1["aa"]("b,i") * tmps_["23_aaaa_Looov"]("L,j,l,k,b")) * t1["aa"]("a,l");
        tmps_["23_aaaa_Looov"].~TArrayD();

        // D_oovv_aaaa += -0.50 r0 t1_aa(a,k) t1_aa(b,l) t2_aaaa(d,c,i,j) l2_aaaa(k,l,d,c)
        //             += -0.50 r0 t1_aa(a,k) t1_aa(b,l) t2_1_aaaa(d,c,i,j) l2_1_aaaa(k,l,d,c)
        //             += +1.00 r0 t1_aa(a,k) t1_aa(b,l) t1_aa(c,i) t1_aa(d,j) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") -= 0.50 * t1["aa"]("b,l") * tmps_["178_aaaa_Looov"]("L,i,j,l,a") * r0("R");

        // flops: o3v1L1  = o4v0L1 o4v1L1
        //  mems: o3v1L1  = o4v0L1 o3v1L1
        tmps_["179_aaaa_Looov"]("L,i,j,k,a")  = (tmps_["177_aaaa_Loooo"]("L,i,j,l,k") + -0.50 * tmps_["173_aaaa_Loooo"]("L,i,j,l,k")) * t1["aa"]("a,l");

        // D_oovv_aaaa += +1.00 r0_1 t1_aa(a,k) t1_aa(b,l) t1_aa(c,i) t1_aa(d,j) l2_1_aaaa(k,l,d,c)
        //             += -0.50 r0_1 t1_aa(a,k) t1_aa(b,l) t2_aaaa(d,c,i,j) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") += t1["aa"]("b,l") * tmps_["179_aaaa_Looov"]("L,i,j,l,a") * r0_1("R");

        // flops: o3v1L1  = o4v0L1 o4v1L1
        //  mems: o3v1L1  = o4v0L1 o3v1L1
        tmps_["180_aaaa_Looov"]("L,i,j,k,a")  = (tmps_["177_aaaa_Loooo"]("L,i,j,l,k") + -0.50 * tmps_["173_aaaa_Loooo"]("L,i,j,l,k")) * t1_1["aa"]("a,l");
        tmps_["177_aaaa_Loooo"].~TArrayD();
        tmps_["173_aaaa_Loooo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(a,b) r0 t1_aa(b,l) t1_aa(c,i) t1_aa(d,j) t1_1_aa(a,k) l2_1_aaaa(k,l,d,c)
        //             += -0.50 P(a,b) r0 t1_aa(b,l) t2_aaaa(d,c,i,j) t1_1_aa(a,k) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("b,l") * tmps_["180_aaaa_Looov"]("L,i,j,l,a") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o3v1L2  = o3v1L2 o3v1L2 o4v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        //  mems: o3v1L2  = o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        tmps_["181_aaaa_LLvooo"]("R,L,a,i,j,k")  = tmps_["174_aaaa_Lvooo"]("L,a,i,j,k") * r0_1("R");
        tmps_["181_aaaa_LLvooo"]("R,L,a,i,j,k") -= 0.50 * tmps_["178_aaaa_Looov"]("L,i,j,k,a") * r0("R");
        tmps_["181_aaaa_LLvooo"]("R,L,a,i,j,k") -= 0.50 * r1["aa"]("R,a,l") * tmps_["175_aaaa_Loooo"]("L,i,j,l,k");
        tmps_["181_aaaa_LLvooo"]("R,L,a,i,j,k") += tmps_["176_aaaa_Lvooo"]("L,a,i,j,k") * r0("R");
        tmps_["181_aaaa_LLvooo"]("R,L,a,i,j,k") += r0("R") * tmps_["180_aaaa_Looov"]("L,i,j,k,a");
        tmps_["181_aaaa_LLvooo"]("R,L,a,i,j,k") += r0_1("R") * tmps_["179_aaaa_Looov"]("L,i,j,k,a");
        tmps_["180_aaaa_Looov"].~TArrayD();
        tmps_["179_aaaa_Looov"].~TArrayD();
        tmps_["178_aaaa_Looov"].~TArrayD();
        tmps_["176_aaaa_Lvooo"].~TArrayD();
        tmps_["175_aaaa_Loooo"].~TArrayD();
        tmps_["174_aaaa_Lvooo"].~TArrayD();

        // D_ooov_aaaa += -1.00 r0_1 t2_aaaa(a,b,i,j) l1_1_aa(k,b)
        //             += -1.00 r0_1 t1_aa(a,l) t1_aa(b,i) t1_aa(c,j) l2_1_aaaa(l,k,c,b)
        //             += +0.50 r0_1 t1_aa(a,l) t2_aaaa(c,b,i,j) l2_1_aaaa(l,k,c,b)
        //             += -1.00 r0 t1_aa(b,i) t1_aa(c,j) t1_1_aa(a,l) l2_1_aaaa(l,k,c,b)
        //             += +0.50 r0 t2_aaaa(c,b,i,j) t1_1_aa(a,l) l2_1_aaaa(l,k,c,b)
        //             += +0.50 r0 t1_aa(a,l) t2_aaaa(c,b,i,j) l2_aaaa(l,k,c,b)
        //             += +0.50 r0 t1_aa(a,l) t2_1_aaaa(c,b,i,j) l2_1_aaaa(l,k,c,b)
        //             += -1.00 r0 t1_aa(a,l) t1_aa(b,i) t1_aa(c,j) l2_aaaa(l,k,c,b)
        //             += +0.50 r1_aa(a,l) t2_aaaa(c,b,i,j) l2_aaaa(l,k,c,b)
        //             += +0.50 r1_aa(a,l) t2_1_aaaa(c,b,i,j) l2_1_aaaa(l,k,c,b)
        //             += -1.00 r0 t2_aaaa(a,b,i,j) l1_aa(k,b)
        //             += -1.00 r0 t2_1_aaaa(a,b,i,j) l1_1_aa(k,b)
        D_ooov_aaaa("L,R,i,j,k,a") -= tmps_["181_aaaa_LLvooo"]("R,L,a,i,j,k");

        // D_oovo_aaaa += +1.00 r0_1 t2_aaaa(a,b,i,j) l1_1_aa(k,b)
        //             += +1.00 r0_1 t1_aa(a,l) t1_aa(b,i) t1_aa(c,j) l2_1_aaaa(l,k,c,b)
        //             += -0.50 r0_1 t1_aa(a,l) t2_aaaa(c,b,i,j) l2_1_aaaa(l,k,c,b)
        //             += +1.00 r0 t1_aa(b,i) t1_aa(c,j) t1_1_aa(a,l) l2_1_aaaa(l,k,c,b)
        //             += -0.50 r0 t2_aaaa(c,b,i,j) t1_1_aa(a,l) l2_1_aaaa(l,k,c,b)
        //             += -0.50 r0 t1_aa(a,l) t2_aaaa(c,b,i,j) l2_aaaa(l,k,c,b)
        //             += -0.50 r0 t1_aa(a,l) t2_1_aaaa(c,b,i,j) l2_1_aaaa(l,k,c,b)
        //             += +1.00 r0 t1_aa(a,l) t1_aa(b,i) t1_aa(c,j) l2_aaaa(l,k,c,b)
        //             += -0.50 r1_aa(a,l) t2_aaaa(c,b,i,j) l2_aaaa(l,k,c,b)
        //             += -0.50 r1_aa(a,l) t2_1_aaaa(c,b,i,j) l2_1_aaaa(l,k,c,b)
        //             += +1.00 r0 t2_aaaa(a,b,i,j) l1_aa(k,b)
        //             += +1.00 r0 t2_1_aaaa(a,b,i,j) l1_1_aa(k,b)
        D_oovo_aaaa("L,R,i,j,a,k") += tmps_["181_aaaa_LLvooo"]("R,L,a,i,j,k");
        tmps_["181_aaaa_LLvooo"].~TArrayD();

        // flops: o2v2L2  = o3v2L2 o3v2L2
        //  mems: o2v2L2  = o3v1L2 o2v2L2
        tmps_["182_aaaa_LLovov"]("L,R,i,a,j,b")  = l2_1["aaaa"]("L,k,i,a,c") * r1_1["aa"]("R,c,j") * t1["aa"]("b,k");

        // D_ovov_aaaa += +1.00 r1_1_aa(c,i) t1_aa(b,k) l2_1_aaaa(k,j,a,c)
        D_ovov_aaaa("L,R,i,a,j,b") += tmps_["182_aaaa_LLovov"]("L,R,j,a,i,b");

        // D_ovvo_aaaa += -1.00 r1_1_aa(c,i) t1_aa(b,k) l2_1_aaaa(k,j,a,c)
        D_ovvo_aaaa("L,R,i,a,b,j") -= tmps_["182_aaaa_LLovov"]("L,R,j,a,i,b");

        // flops: o2v2L2  = o3v2L2 o3v2L2
        //  mems: o2v2L2  = o3v1L2 o2v2L2
        tmps_["183_aaaa_LLovov"]("L,R,i,a,j,b")  = l2["aaaa"]("L,k,i,a,c") * r1["aa"]("R,c,j") * t1["aa"]("b,k");

        // D_ovov_aaaa += +1.00 r1_aa(c,i) t1_aa(b,k) l2_aaaa(k,j,a,c)
        D_ovov_aaaa("L,R,i,a,j,b") += tmps_["183_aaaa_LLovov"]("L,R,j,a,i,b");

        // D_ovvo_aaaa += -1.00 r1_aa(c,i) t1_aa(b,k) l2_aaaa(k,j,a,c)
        D_ovvo_aaaa("L,R,i,a,b,j") -= tmps_["183_aaaa_LLovov"]("L,R,j,a,i,b");

        // flops: o0v4L1  = o2v4L1
        //  mems: o0v4L1  = o0v4L1
        tmps_["184_aaaa_Lvvvv"]("L,a,b,c,d")  = l2_1["aaaa"]("L,i,j,a,b") * t2["aaaa"]("c,d,i,j");

        // D_oovv_aaaa += -0.50 P(i,j) r1_aa(d,i) t2_aaaa(b,a,k,l) t1_1_aa(c,j) l2_1_aaaa(k,l,c,d)
        // flops: o2v2L2 += o1v4L1 o2v3L2
        //  mems: o2v2L2 += o1v3L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["184_aaaa_Lvvvv"]("L,c,d,b,a") * t1_1["aa"]("c,j") * r1["aa"]("R,d,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_vvvv_aaaa += +0.50 r0_1 t2_aaaa(d,c,i,j) l2_1_aaaa(i,j,a,b)
        // flops: o0v4L2 += o0v4L2
        //  mems: o0v4L2 += o0v4L2
        D_vvvv_aaaa("L,R,a,b,c,d") += 0.50 * tmps_["184_aaaa_Lvvvv"]("L,a,b,d,c") * r0_1("R");

        // flops: o0v4L1  = o2v4L1
        //  mems: o0v4L1  = o0v4L1
        tmps_["185_aaaa_Lvvvv"]("L,a,b,c,d")  = l2["aaaa"]("L,i,j,a,b") * t2["aaaa"]("c,d,i,j");

        // D_vvvv_aaaa += +0.50 r0 t2_aaaa(d,c,i,j) l2_aaaa(i,j,a,b)
        // flops: o0v4L2 += o0v4L2
        //  mems: o0v4L2 += o0v4L2
        D_vvvv_aaaa("L,R,a,b,c,d") += 0.50 * tmps_["185_aaaa_Lvvvv"]("L,a,b,d,c") * r0("R");

        // flops: o0v4L1  = o2v4L1
        //  mems: o0v4L1  = o0v4L1
        tmps_["186_aaaa_Lvvvv"]("L,a,b,c,d")  = l2_1["aaaa"]("L,i,j,a,b") * t2_1["aaaa"]("c,d,i,j");

        // D_vvvv_aaaa += +0.50 r0 t2_1_aaaa(d,c,i,j) l2_1_aaaa(i,j,a,b)
        // flops: o0v4L2 += o0v4L2
        //  mems: o0v4L2 += o0v4L2
        D_vvvv_aaaa("L,R,a,b,c,d") += 0.50 * tmps_["186_aaaa_Lvvvv"]("L,a,b,d,c") * r0("R");

        // flops: o1v3L1  = o1v4L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["187_aaaa_Lvvvo"]("L,a,b,c,i")  = tmps_["184_aaaa_Lvvvv"]("L,a,d,b,c") * t1["aa"]("d,i");

        // D_oovv_aaaa += -0.50 r0_1 t2_aaaa(b,a,k,l) t1_aa(c,i) t1_aa(d,j) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") -= 0.50 * tmps_["187_aaaa_Lvvvo"]("L,d,b,a,i") * t1["aa"]("d,j") * r0_1("R");

        // flops: o1v3L1  = o1v4L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["188_aaaa_Lvvvo"]("L,a,b,c,i")  = tmps_["184_aaaa_Lvvvv"]("L,a,d,b,c") * t1_1["aa"]("d,i");

        // D_oovv_aaaa += -0.50 P(i,j) r0 t2_aaaa(b,a,k,l) t1_aa(d,j) t1_1_aa(c,i) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["188_aaaa_Lvvvo"]("L,d,b,a,i") * t1["aa"]("d,j") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o1v3L1  = o0v4L1 o1v4L1
        //  mems: o1v3L1  = o0v4L1 o1v3L1
        tmps_["189_aaaa_Lvvvo"]("L,a,b,c,i")  = (tmps_["186_aaaa_Lvvvv"]("L,a,d,b,c") + tmps_["185_aaaa_Lvvvv"]("L,a,d,b,c")) * t1["aa"]("d,i");

        // D_oovv_aaaa += -0.50 r0 t1_aa(c,i) t1_aa(d,j) t2_1_aaaa(b,a,k,l) l2_1_aaaa(k,l,d,c)
        //             += -0.50 r0 t2_aaaa(b,a,k,l) t1_aa(c,i) t1_aa(d,j) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") -= 0.50 * tmps_["189_aaaa_Lvvvo"]("L,d,b,a,i") * t1["aa"]("d,j") * r0("R");

        // flops: o1v3L2  = o2v3L2 o1v4L2 o1v4L2 o1v3L2 o1v4L2 o1v3L2 o1v3L2 o3v2L1 o3v2L1 o2v2L2 o2v2L2 o2v3L2 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L1 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        //  mems: o1v3L2  = o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o2v2L2 o2v2L2 o1v3L2 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L1 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i")  = -2.00 * tmps_["129_aa_LLov"]("L,R,j,a") * t2["aaaa"]("b,c,j,i");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") += tmps_["184_aaaa_Lvvvv"]("L,a,d,b,c") * r1_1["aa"]("R,d,i");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") += tmps_["186_aaaa_Lvvvv"]("L,a,d,b,c") * r1["aa"]("R,d,i");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") += tmps_["185_aaaa_Lvvvv"]("L,a,d,b,c") * r1["aa"]("R,d,i");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * (tmps_["183_aaaa_LLovov"]("L,R,k,a,i,c") + t1["aa"]("d,i") * l2_1["aaaa"]("L,j,k,a,d") * t1["aa"]("c,j") * r0_1("R")) * t1["aa"]("b,k");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * tmps_["118_aa_LLov"]("L,R,j,a") * t2_1["aaaa"]("b,c,j,i");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l2["aaaa"]("L,j,k,a,d") * t1["aa"]("d,i") * t1["aa"]("c,j") * t1["aa"]("b,k") * r0("R");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l1_1["aa"]("L,j,a") * t2["aaaa"]("b,c,j,i") * r0_1("R");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") += 2.00 * tmps_["117_aa_LLov"]("L,R,j,a") * t2_1["aaaa"]("b,c,j,i");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") += 2.00 * tmps_["128_aa_LLov"]("L,R,j,a") * t2["aaaa"]("b,c,j,i");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l1["aa"]("L,j,a") * t2["aaaa"]("b,c,j,i") * r0("R");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l1["aa"]("L,j,a") * r2["aaaa"]("R,b,c,j,i");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * tmps_["182_aaaa_LLovov"]("L,R,k,a,i,c") * t1["aa"]("b,k");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l1_1["aa"]("L,j,a") * r2_1["aaaa"]("R,b,c,j,i");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * tmps_["130_aa_LLov"]("L,R,j,a") * t2["aaaa"]("b,c,j,i");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") += 2.00 * tmps_["127_aa_LLov"]("L,R,j,a") * t2["aaaa"]("b,c,j,i");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l1_1["aa"]("L,j,a") * t2_1["aaaa"]("b,c,j,i") * r0("R");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l2_1["aaaa"]("L,j,k,a,d") * t1_1["aa"]("d,i") * t1["aa"]("c,j") * t1["aa"]("b,k") * r0("R");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") += tmps_["189_aaaa_Lvvvo"]("L,a,b,c,i") * r0("R");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") += r0("R") * tmps_["188_aaaa_Lvvvo"]("L,a,b,c,i");
        tmps_["190_aaaa_LLvvvo"]("R,L,a,b,c,i") += r0_1("R") * tmps_["187_aaaa_Lvvvo"]("L,a,b,c,i");
        tmps_["189_aaaa_Lvvvo"].~TArrayD();
        tmps_["188_aaaa_Lvvvo"].~TArrayD();
        tmps_["187_aaaa_Lvvvo"].~TArrayD();
        tmps_["186_aaaa_Lvvvv"].~TArrayD();
        tmps_["185_aaaa_Lvvvv"].~TArrayD();
        tmps_["184_aaaa_Lvvvv"].~TArrayD();
        tmps_["183_aaaa_LLovov"].~TArrayD();
        tmps_["182_aaaa_LLovov"].~TArrayD();

        // D_ovvv_aaaa += +0.50 r0_1 t2_aaaa(c,b,j,k) t1_aa(d,i) l2_1_aaaa(j,k,a,d)
        //             += +0.50 r0 t2_aaaa(c,b,j,k) t1_1_aa(d,i) l2_1_aaaa(j,k,a,d)
        //             += -1.00 r1_bb(d,k) t2_aaaa(c,b,j,i) l2_abab(j,k,a,d)
        //             += +0.50 r1_1_aa(d,i) t2_aaaa(c,b,j,k) l2_1_aaaa(j,k,a,d)
        //             += +0.50 r1_aa(d,i) t2_1_aaaa(c,b,j,k) l2_1_aaaa(j,k,a,d)
        //             += +0.50 r1_aa(d,i) t2_aaaa(c,b,j,k) l2_aaaa(j,k,a,d)
        //             += -1.00 r1_aa(d,i) t1_aa(b,j) t1_aa(c,k) l2_aaaa(j,k,a,d)
        //             += -1.00 r0_1 t1_aa(b,j) t1_aa(c,k) t1_aa(d,i) l2_1_aaaa(j,k,a,d)
        //             += -1.00 r1_bb(d,k) t2_1_aaaa(c,b,j,i) l2_1_abab(j,k,a,d)
        //             += -1.00 r0 t1_aa(b,j) t1_aa(c,k) t1_aa(d,i) l2_aaaa(j,k,a,d)
        //             += -1.00 r0_1 t2_aaaa(c,b,j,i) l1_1_aa(j,a)
        //             += +1.00 r1_aa(d,k) t2_1_aaaa(c,b,j,i) l2_1_aaaa(k,j,a,d)
        //             += +1.00 r1_1_aa(d,k) t2_aaaa(c,b,j,i) l2_1_aaaa(k,j,a,d)
        //             += -1.00 r0 t2_aaaa(c,b,j,i) l1_aa(j,a)
        //             += -1.00 r2_aaaa(c,b,j,i) l1_aa(j,a)
        //             += -1.00 r1_1_aa(d,i) t1_aa(b,j) t1_aa(c,k) l2_1_aaaa(j,k,a,d)
        //             += -1.00 r2_1_aaaa(c,b,j,i) l1_1_aa(j,a)
        //             += -1.00 r1_1_bb(d,k) t2_aaaa(c,b,j,i) l2_1_abab(j,k,a,d)
        //             += +1.00 r1_aa(d,k) t2_aaaa(c,b,j,i) l2_aaaa(k,j,a,d)
        //             += -1.00 r0 t2_1_aaaa(c,b,j,i) l1_1_aa(j,a)
        //             += -1.00 r0 t1_aa(b,j) t1_aa(c,k) t1_1_aa(d,i) l2_1_aaaa(j,k,a,d)
        //             += +0.50 r0 t1_aa(d,i) t2_1_aaaa(c,b,j,k) l2_1_aaaa(j,k,a,d)
        //             += +0.50 r0 t2_aaaa(c,b,j,k) t1_aa(d,i) l2_aaaa(j,k,a,d)
        D_ovvv_aaaa("L,R,i,a,b,c") += 0.50 * tmps_["190_aaaa_LLvvvo"]("R,L,a,c,b,i");

        // D_vovv_aaaa += -0.50 r0_1 t2_aaaa(c,b,j,k) t1_aa(d,i) l2_1_aaaa(j,k,a,d)
        //             += -0.50 r0 t2_aaaa(c,b,j,k) t1_1_aa(d,i) l2_1_aaaa(j,k,a,d)
        //             += +1.00 r1_bb(d,k) t2_aaaa(c,b,j,i) l2_abab(j,k,a,d)
        //             += -0.50 r1_1_aa(d,i) t2_aaaa(c,b,j,k) l2_1_aaaa(j,k,a,d)
        //             += -0.50 r1_aa(d,i) t2_1_aaaa(c,b,j,k) l2_1_aaaa(j,k,a,d)
        //             += -0.50 r1_aa(d,i) t2_aaaa(c,b,j,k) l2_aaaa(j,k,a,d)
        //             += +1.00 r1_aa(d,i) t1_aa(b,j) t1_aa(c,k) l2_aaaa(j,k,a,d)
        //             += +1.00 r0_1 t1_aa(b,j) t1_aa(c,k) t1_aa(d,i) l2_1_aaaa(j,k,a,d)
        //             += +1.00 r1_bb(d,k) t2_1_aaaa(c,b,j,i) l2_1_abab(j,k,a,d)
        //             += +1.00 r0 t1_aa(b,j) t1_aa(c,k) t1_aa(d,i) l2_aaaa(j,k,a,d)
        //             += +1.00 r0_1 t2_aaaa(c,b,j,i) l1_1_aa(j,a)
        //             += -1.00 r1_aa(d,k) t2_1_aaaa(c,b,j,i) l2_1_aaaa(k,j,a,d)
        //             += -1.00 r1_1_aa(d,k) t2_aaaa(c,b,j,i) l2_1_aaaa(k,j,a,d)
        //             += +1.00 r0 t2_aaaa(c,b,j,i) l1_aa(j,a)
        //             += +1.00 r2_aaaa(c,b,j,i) l1_aa(j,a)
        //             += +1.00 r1_1_aa(d,i) t1_aa(b,j) t1_aa(c,k) l2_1_aaaa(j,k,a,d)
        //             += +1.00 r2_1_aaaa(c,b,j,i) l1_1_aa(j,a)
        //             += +1.00 r1_1_bb(d,k) t2_aaaa(c,b,j,i) l2_1_abab(j,k,a,d)
        //             += -1.00 r1_aa(d,k) t2_aaaa(c,b,j,i) l2_aaaa(k,j,a,d)
        //             += +1.00 r0 t2_1_aaaa(c,b,j,i) l1_1_aa(j,a)
        //             += +1.00 r0 t1_aa(b,j) t1_aa(c,k) t1_1_aa(d,i) l2_1_aaaa(j,k,a,d)
        //             += -0.50 r0 t1_aa(d,i) t2_1_aaaa(c,b,j,k) l2_1_aaaa(j,k,a,d)
        //             += -0.50 r0 t2_aaaa(c,b,j,k) t1_aa(d,i) l2_aaaa(j,k,a,d)
        D_vovv_aaaa("L,R,a,i,b,c") -= 0.50 * tmps_["190_aaaa_LLvvvo"]("R,L,a,c,b,i");
        tmps_["190_aaaa_LLvvvo"].~TArrayD();

        // flops: o0v0L2  = o0v0L2 o0v0L2 o0v0L2
        //  mems: o0v0L2  = o0v0L2 o0v0L2 o0v0L2
        tmps_["191_LL"]("L,R")  = l0("L") * r0("R");
        tmps_["191_LL"]("L,R") += l0_1("L") * r0_1("R");

        // D_oovv_aaaa += +1.00 r0_1 t2_aaaa(b,a,i,j) l0_1
        //             += +1.00 r0 t2_aaaa(b,a,i,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") += t2["aaaa"]("b,a,i,j") * tmps_["191_LL"]("L,R");

        // D_oovv_abab += -1.00 r0_1 t2_abab(a,b,i,j) l0_1
        //             += -1.00 r0 t2_abab(a,b,i,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,b,i,j") * tmps_["191_LL"]("L,R");

        // D_oovv_bbbb += +1.00 r0_1 t2_bbbb(b,a,i,j) l0_1
        //             += +1.00 r0 t2_bbbb(b,a,i,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += t2["bbbb"]("b,a,i,j") * tmps_["191_LL"]("L,R");

        // flops: o1v1L1  = o1v1L1
        //  mems: o1v1L1  = o1v1L1
        tmps_["192_bb_Lvo"]("L,a,i")  = l0("L") * t1["bb"]("a,i");

        // D_oovv_abab += -1.00 r1_aa(a,i) t1_bb(b,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["192_bb_Lvo"]("L,b,j") * r1["aa"]("R,a,i");

        // D_oovv_bbbb += -1.00 P(i,j) P(a,b) r1_bb(a,i) t1_bb(b,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["192_bb_Lvo"]("L,b,j") * r1["bb"]("R,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v1L1  = o3v2L1 o3v2L1 o3v0L1 o3v1L1
        //  mems: o2v1L1  = o2v0L1 o2v0L1 o3v0L1 o2v1L1
        tmps_["193_aaa_Loov"]("L,i,j,a")  = (t2_1["aaaa"]("b,c,j,i") * l2_1["aaaa"]("L,j,k,b,c") + t2["aaaa"]("c,d,k,i") * l2["aaaa"]("L,k,j,c,d")) * t1["aa"]("a,k");

        // D_oovo_aaaa += +0.50 P(i,j) d_aa(j,k) r0 t1_aa(a,m) t2_1_aaaa(c,b,l,i) l2_1_aaaa(l,m,c,b)
        //             += +0.50 d_bb(j,l) r0 t2_aaaa(b,a,m,i) l2_aaaa(m,k,b,a)
        // flops: o3v1L2 += o3v1L1 o2v1L2
        //  mems: o3v1L2 += o2v1L1 o2v1L2
        tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k")  = 0.50 * tmps_["193_aaa_Loov"]("L,i,k,a") * Id["aa_oo"]("j,k") * r0("R");
        D_oovo_aaaa("L,R,i,j,a,k") += tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k");
        D_oovo_aaaa("L,R,i,j,a,k") -= tmps_["perm_aaaa_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_aaaa_LLvooo"].~TArrayD();

        // flops: o3v1L1  = o3v2L1 o3v2L1 o4v0L1 o4v1L1
        //  mems: o3v1L1  = o2v0L1 o2v0L1 o4v0L1 o3v1L1
        tmps_["194_aaaa_Looov"]("L,i,j,k,a")  = (t2_1["abab"]("b,c,i,m") * l2_1["abab"]("L,l,m,b,c") + -1.00 * t2["abab"]("d,e,j,n") * l2["abab"]("L,k,n,d,e")) * t1["aa"]("a,l");

        // D_oovo_abab += -0.50 d_bb(i,k) r0 t1_aa(a,m) t2_1_abab(c,b,j,l) l2_1_abab(m,l,c,b)
        //             += -0.50 d_bb(i,k) r0 t1_aa(a,m) t2_1_abab(b,c,j,l) l2_1_abab(m,l,b,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,l) t1_aa(b,j) t2_abab(d,c,i,k) l2_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,l) t1_aa(b,j) t2_abab(c,d,i,k) l2_abab(l,k,c,d)
        // flops: o3v1L2 += o5v1L2
        //  mems: o3v1L2 += o5v1L2
        D_oovo_abab("L,R,i,j,a,k") -= tmps_["194_aaaa_Looov"]("L,j,i,l,a") * tmps_["112_bb_Loo"]("R,i,k");

        // flops: o2v1L1  = o3v2L1 o3v2L1 o3v0L1 o3v1L1
        //  mems: o2v1L1  = o2v0L1 o2v0L1 o3v0L1 o2v1L1
        tmps_["195_bbb_Loov"]("L,i,j,a")  = (t2["bbbb"]("b,c,j,i") * l2["bbbb"]("L,j,k,b,c") + t2_1["bbbb"]("c,a,k,i") * l2_1["bbbb"]("L,k,j,c,a")) * t1["bb"]("a,k");

        // D_oovv_abab += +0.50 r1_aa(a,i) t1_bb(b,l) t2_bbbb(d,c,k,j) l2_bbbb(k,l,d,c)
        //             += +0.50 P(i,j) r1_bb(a,i) t2_1_bbbb(c,b,l,j) l2_1_bbbb(l,k,c,b)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o3v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * tmps_["195_bbb_Loov"]("L,j,k,b") * r1["aa"]("R,a,i");

        // flops: o4v0  = o4v0
        //  mems: o4v0  = o4v0
        tmps_["196_aaaa_oooo"]("i,j,k,l")  = Id["aa_oo"]("i,j") * Id["aa_oo"]("k,l");

        // D_oooo_aaaa += +1.00 P(i,j) d_aa(i,l) d_aa(j,k) r0_1 l0_1
        //             += +1.00 P(i,j) d_aa(i,l) d_aa(j,k) r0 l0
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l")  = tmps_["191_LL"]("L,R") * tmps_["196_aaaa_oooo"]("i,l,j,k");
        D_oooo_aaaa("L,R,i,j,k,l") += tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l");
        D_oooo_aaaa("L,R,i,j,k,l") -= tmps_["perm_aaaa_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_aaaa_LLoooo"].~TArrayD();

        // flops: o4v0  = o4v0
        //  mems: o4v0  = o4v0
        tmps_["197_aabb_oooo"]("i,j,k,l")  = Id["aa_oo"]("i,j") * Id["bb_oo"]("k,l");

        // D_oooo_abab += -1.00 d_aa(i,k) d_bb(j,l) r0_1 l0_1
        //             += -1.00 d_aa(i,k) d_bb(j,l) r0 l0
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") -= tmps_["191_LL"]("L,R") * tmps_["197_aabb_oooo"]("i,k,j,l");

        // flops: o4v0  = o4v0
        //  mems: o4v0  = o4v0
        tmps_["198_bbbb_oooo"]("i,j,k,l")  = Id["bb_oo"]("i,j") * Id["bb_oo"]("k,l");

        // D_oooo_bbbb += +1.00 P(i,j) d_bb(i,l) d_bb(j,k) r0_1 l0_1
        //             += +1.00 P(i,j) d_bb(i,l) d_bb(j,k) r0 l0
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        tmps_["perm_bbbb_LLoooo"]("L,R,i,j,k,l")  = tmps_["191_LL"]("L,R") * tmps_["198_bbbb_oooo"]("i,l,j,k");
        D_oooo_bbbb("L,R,i,j,k,l") += tmps_["perm_bbbb_LLoooo"]("L,R,i,j,k,l");
        D_oooo_bbbb("L,R,i,j,k,l") -= tmps_["perm_bbbb_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_bbbb_LLoooo"].~TArrayD();
        tmps_["191_LL"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["199_bbbb_Lvoov"]("L,a,i,j,b")  = t2["bbbb"]("a,c,k,i") * l2["bbbb"]("L,j,k,c,b");

        // D_oovv_abab += -1.00 r2_abab(a,d,i,l) t2_bbbb(b,c,k,j) l2_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["199_bbbb_Lvoov"]("L,b,j,l,d") * r2["abab"]("R,a,d,i,l");

        // D_oovv_abab += -1.00 r0 t2_abab(a,c,i,k) t2_bbbb(b,d,l,j) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,c,i,k") * tmps_["199_bbbb_Lvoov"]("L,b,j,k,c") * r0("R");

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["200_bb_Lvv"]("L,a,b")  = t2["bbbb"]("a,c,i,j") * l2["bbbb"]("L,i,j,c,b");

        // D_oovv_abab += -0.50 r2_abab(a,d,i,j) t2_bbbb(b,c,k,l) l2_bbbb(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * tmps_["200_bb_Lvv"]("L,b,d") * r2["abab"]("R,a,d,i,j");

        // D_oovv_abab += -0.50 r0 t2_abab(a,c,i,j) t2_bbbb(b,d,k,l) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * t2["abab"]("a,c,i,j") * tmps_["200_bb_Lvv"]("L,b,c") * r0("R");

        // D_oovv_bbbb += -0.50 P(a,b) r2_bbbb(a,d,i,j) t2_bbbb(b,c,k,l) l2_bbbb(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["200_bb_Lvv"]("L,b,d") * r2["bbbb"]("R,a,d,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += -0.50 r0 t2_bbbb(a,c,i,j) t2_bbbb(b,d,k,l) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") -= 0.50 * t2["bbbb"]("a,c,i,j") * tmps_["200_bb_Lvv"]("L,b,c") * r0("R");

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["201_bb_Lvv"]("L,a,b")  = t2_1["bbbb"]("a,c,i,j") * l2_1["bbbb"]("L,i,j,c,b");

        // D_oovv_abab += -0.50 r2_abab(a,d,i,j) t2_1_bbbb(b,c,k,l) l2_1_bbbb(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * tmps_["201_bb_Lvv"]("L,b,d") * r2["abab"]("R,a,d,i,j");

        // D_oovv_bbbb += -0.50 P(a,b) r2_bbbb(a,d,i,j) t2_1_bbbb(b,c,k,l) l2_1_bbbb(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["201_bb_Lvv"]("L,b,d") * r2["bbbb"]("R,a,d,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v0L1  = o3v2L1 o3v2L1 o2v0L1
        //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
        tmps_["202_bb_Loo"]("L,i,j")  = t2["bbbb"]("a,b,k,i") * l2["bbbb"]("L,j,k,a,b");
        tmps_["202_bb_Loo"]("L,i,j") += t2_1["bbbb"]("a,b,k,i") * l2_1["bbbb"]("L,j,k,a,b");

        // D_oovv_abab += -0.50 r2_abab(a,b,i,l) t2_1_bbbb(d,c,k,j) l2_1_bbbb(l,k,d,c)
        //             += -0.50 r2_abab(a,b,i,l) t2_bbbb(d,c,k,j) l2_bbbb(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * tmps_["202_bb_Loo"]("L,j,l") * r2["abab"]("R,a,b,i,l");

        // D_oovv_bbbb += -0.50 P(i,j) r2_bbbb(b,a,l,i) t2_1_bbbb(d,c,k,j) l2_1_bbbb(l,k,d,c)
        //             += -0.50 P(i,j) r2_bbbb(b,a,l,i) t2_bbbb(d,c,k,j) l2_bbbb(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["202_bb_Loo"]("L,j,l") * r2["bbbb"]("R,b,a,l,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v1L2  = o2v2L2 o2v1L1 o2v1L2 o1v1L2 o2v2L2 o1v1L2 o2v2L2 o1v1L2 o2v2L2 o1v1L2 o2v1L1 o2v1L2 o1v1L2 o2v1L1 o2v1L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o2v0L1 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o2v0L1 o1v1L2 o1v1L2 o2v0L1 o1v1L2 o1v1L2
        tmps_["203_bb_LLov"]("L,R,i,a")  = -1.00 * l1["aa"]("L,k,c") * r2["abab"]("R,c,a,k,i");
        tmps_["203_bb_LLov"]("L,R,i,a") += l1_1["bb"]("L,j,b") * t1["bb"]("b,i") * r1_1["bb"]("R,a,j");
        tmps_["203_bb_LLov"]("L,R,i,a") -= l1_1["aa"]("L,k,c") * r2_1["abab"]("R,c,a,k,i");
        tmps_["203_bb_LLov"]("L,R,i,a") += l1_1["bb"]("L,j,b") * r2_1["bbbb"]("R,a,b,j,i");
        tmps_["203_bb_LLov"]("L,R,i,a") += l1["bb"]("L,j,b") * r2["bbbb"]("R,a,b,j,i");
        tmps_["203_bb_LLov"]("L,R,i,a") += l1_1["bb"]("L,j,b") * t1_1["bb"]("b,i") * r1["bb"]("R,a,j");
        tmps_["203_bb_LLov"]("L,R,i,a") += l1["bb"]("L,j,b") * t1["bb"]("b,i") * r1["bb"]("R,a,j");

        // D_oovv_abab += +1.00 r1_1_bb(b,k) t1_aa(a,i) t1_bb(c,j) l1_1_bb(k,c)
        //             += -1.00 r2_abab(c,b,k,j) t1_aa(a,i) l1_aa(k,c)
        //             += -1.00 r2_1_abab(c,b,k,j) t1_aa(a,i) l1_1_aa(k,c)
        //             += +1.00 r2_1_bbbb(b,c,k,j) t1_aa(a,i) l1_1_bb(k,c)
        //             += +1.00 r2_bbbb(b,c,k,j) t1_aa(a,i) l1_bb(k,c)
        //             += +1.00 r1_bb(b,k) t1_aa(a,i) t1_1_bb(c,j) l1_1_bb(k,c)
        //             += +1.00 r1_bb(b,k) t1_aa(a,i) t1_bb(c,j) l1_bb(k,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1["aa"]("a,i") * tmps_["203_bb_LLov"]("L,R,j,b");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r1_1_bb(a,k) t1_bb(b,j) t1_bb(c,i) l1_1_bb(k,c)
        //             += -1.00 P(i,j) P(a,b) r2_abab(c,a,k,i) t1_bb(b,j) l1_aa(k,c)
        //             += -1.00 P(i,j) P(a,b) r2_1_abab(c,a,k,i) t1_bb(b,j) l1_1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r2_1_bbbb(a,c,k,i) t1_bb(b,j) l1_1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r2_bbbb(a,c,k,i) t1_bb(b,j) l1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_bb(a,k) t1_bb(b,j) t1_1_bb(c,i) l1_1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_bb(a,k) t1_bb(b,j) t1_bb(c,i) l1_bb(k,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("b,j") * tmps_["203_bb_LLov"]("L,R,i,a");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["204_bbbb_Lvooo"]("L,a,i,j,k")  = tmps_["81_bbbb_Lvoov"]("L,a,i,j,b") * t1["bb"]("b,k");

        // D_oovv_bbbb += -1.00 P(i,j) P(a,b) r1_1_bb(a,l) t2_abab(c,b,k,i) t1_bb(d,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["204_bbbb_Lvooo"]("L,b,i,l,j") * r1_1["bb"]("R,a,l");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += -1.00 P(i,j) P(a,b) r0 t2_abab(c,b,l,i) t1_bb(d,j) t1_1_bb(a,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1_1["bb"]("a,k") * tmps_["204_bbbb_Lvooo"]("L,b,i,k,j") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r0_1 t2_abab(c,a,k,i) t1_bb(b,l) t1_bb(d,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["204_bbbb_Lvooo"]("L,a,i,l,j") * t1["bb"]("b,l") * r0_1("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["205_bbbb_Lvoov"]("L,a,i,j,b")  = t2["bbbb"]("a,c,k,i") * l2_1["bbbb"]("L,k,j,b,c");

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["206_bbbb_Lvooo"]("L,a,i,j,k")  = tmps_["205_bbbb_Lvoov"]("L,a,i,j,b") * t1["bb"]("b,k");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r0_1 t2_bbbb(a,c,k,i) t1_bb(b,l) t1_bb(d,j) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["206_bbbb_Lvooo"]("L,a,i,l,j") * t1["bb"]("b,l") * r0_1("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o4v1L1 o4v1L1
        //  mems: o3v1L1  = o4v0L1 o3v1L1
        tmps_["207_bbbb_Looov"]("L,i,j,k,a")  = t1_1["bb"]("b,i") * tmps_["166_bbbb_Loovo"]("L,l,j,b,k") * t1["bb"]("a,l");
        tmps_["166_bbbb_Loovo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) r0 t1_bb(a,k) t1_bb(b,l) t1_bb(d,j) t1_1_bb(c,i) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("b,l") * tmps_["207_bbbb_Looov"]("L,i,l,j,a") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["208_bbbb_Lvoov"]("L,a,i,j,b")  = t2_1["abab"]("c,a,k,i") * l2_1["abab"]("L,k,j,c,b");

        // D_oovv_abab += -1.00 r2_abab(a,d,i,l) t2_1_abab(c,b,k,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["208_bbbb_Lvoov"]("L,b,j,l,d") * r2["abab"]("R,a,d,i,l");

        // D_oovv_abab += -1.00 r0 t2_abab(a,d,i,l) t2_1_abab(c,b,k,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,d,i,l") * tmps_["208_bbbb_Lvoov"]("L,b,j,l,d") * r0("R");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r2_bbbb(a,d,l,i) t2_1_abab(c,b,k,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["208_bbbb_Lvoov"]("L,b,j,l,d") * r2["bbbb"]("R,a,d,l,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r0 t2_bbbb(b,d,l,j) t2_1_abab(c,a,k,i) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["208_bbbb_Lvoov"]("L,a,i,l,d") * t2["bbbb"]("b,d,l,j") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1 o3v3L1 o2v2L1
        //  mems: o2v2L1  = o2v2L1 o2v2L1 o2v2L1
        tmps_["209_bbaa_Lvoov"]("L,a,i,j,b")  = t2["bbbb"]("a,c,k,i") * l2["abab"]("L,j,k,b,c");
        tmps_["209_bbaa_Lvoov"]("L,a,i,j,b") += t2_1["bbbb"]("a,c,k,i") * l2_1["abab"]("L,j,k,b,c");

        // D_oovv_abab += -1.00 r1_aa(d,i) t1_aa(a,l) t2_bbbb(b,c,k,j) l2_abab(l,k,d,c)
        //             += -1.00 r1_aa(d,i) t1_aa(a,l) t2_1_bbbb(b,c,k,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["209_bbaa_Lvoov"]("L,b,j,l,d") * r1["aa"]("R,d,i") * t1["aa"]("a,l");

        // flops: o0v2L2  = o2v3L2 o2v3L2 o2v3L2 o0v2L2 o2v3L2 o0v2L2 o0v2L2
        //  mems: o0v2L2  = o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2
        tmps_["210_bb_LLvv"]("L,R,a,b")  = 2.00 * l2["abab"]("L,k,j,d,a") * r2["abab"]("R,d,b,k,j");
        tmps_["210_bb_LLvv"]("L,R,a,b") += 2.00 * l2_1["abab"]("L,k,j,d,a") * r2_1["abab"]("R,d,b,k,j");
        tmps_["210_bb_LLvv"]("L,R,a,b") += l2["bbbb"]("L,i,j,a,c") * r2["bbbb"]("R,b,c,i,j");
        tmps_["210_bb_LLvv"]("L,R,a,b") += l2_1["bbbb"]("L,i,j,a,c") * r2_1["bbbb"]("R,b,c,i,j");

        // D_oovv_abab += +0.50 r2_bbbb(b,d,l,k) t2_abab(a,c,i,j) l2_bbbb(l,k,c,d)
        //             += +0.50 r2_1_abab(d,b,l,k) t2_abab(a,c,i,j) l2_1_abab(l,k,d,c)
        //             += +0.50 r2_1_abab(d,b,k,l) t2_abab(a,c,i,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r2_1_bbbb(b,d,l,k) t2_abab(a,c,i,j) l2_1_bbbb(l,k,c,d)
        //             += +0.50 r2_abab(d,b,l,k) t2_abab(a,c,i,j) l2_abab(l,k,d,c)
        //             += +0.50 r2_abab(d,b,k,l) t2_abab(a,c,i,j) l2_abab(k,l,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * t2["abab"]("a,c,i,j") * tmps_["210_bb_LLvv"]("L,R,c,b");

        // D_oovv_bbbb += -0.50 P(a,b) r2_bbbb(a,d,l,k) t2_bbbb(b,c,i,j) l2_bbbb(l,k,c,d)
        //             += -0.50 P(a,b) r2_1_abab(d,a,l,k) t2_bbbb(b,c,i,j) l2_1_abab(l,k,d,c)
        //             += -0.50 P(a,b) r2_1_abab(d,a,k,l) t2_bbbb(b,c,i,j) l2_1_abab(k,l,d,c)
        //             += -0.50 P(a,b) r2_1_bbbb(a,d,l,k) t2_bbbb(b,c,i,j) l2_1_bbbb(l,k,c,d)
        //             += -0.50 P(a,b) r2_abab(d,a,l,k) t2_bbbb(b,c,i,j) l2_abab(l,k,d,c)
        //             += -0.50 P(a,b) r2_abab(d,a,k,l) t2_bbbb(b,c,i,j) l2_abab(k,l,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["210_bb_LLvv"]("L,R,c,a") * t2["bbbb"]("b,c,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
    }
} // hilbert
#endif