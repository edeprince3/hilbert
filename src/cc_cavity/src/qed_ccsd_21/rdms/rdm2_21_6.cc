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


    void EOM_EE_QED_RDM_21::rdm2_21_6() {

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


        // D_oovv_abab += +0.50 r0 t2_abab(a,b,i,l) t2_bbbb(d,c,k,j) l2_bbbb(k,l,d,c)
        //             += +0.50 P(i,j) r0 t1_bb(a,j) t2_1_bbbb(c,b,l,i) l2_1_bbbb(l,k,c,b)
        // flops: o2v2L2 += o5v2L1 o4v2L2
        //  mems: o2v2L2 += o4v2L1 o4v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * t2["abab"]("a,b,i,l") * tmps_["233_bbbb_Loooo"]("L,j,l,i,k") * r0("R");

        // D_oovv_bbbb += +0.50 P(i,j) P(a,b) r1_bb(a,i) t1_bb(b,l) t2_bbbb(d,c,k,j) l2_bbbb(k,l,d,c)
        //             += +0.50 r0 t1_aa(a,j) t2_1_bbbb(c,b,l,i) l2_1_bbbb(l,k,c,b)
        // flops: o2v2L2 += o4v1L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t1["bb"]("b,l") * tmps_["233_bbbb_Loooo"]("L,j,l,i,k") * r1["bb"]("R,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oooo_abab += +0.50 d_aa(i,k) r0 t2_bbbb(b,a,m,j) l2_bbbb(m,l,b,a)
        //             += +0.50 P(i,j) d_bb(j,l) r0 t2_1_bbbb(b,a,m,i) l2_1_bbbb(m,k,b,a)
        // flops: o4v0L2 += o6v0L2
        //  mems: o4v0L2 += o6v0L2
        D_oooo_abab("L,R,i,j,k,l") += 0.50 * tmps_["233_bbbb_Loooo"]("L,j,l,i,k") * tmps_["103_aa_Loo"]("R,i,k");

        // flops: o1v1L1  = o1v1L1
        //  mems: o1v1L1  = o1v1L1
        tmps_["234_aa_Lvo"]("L,a,i")  = l0_1("L") * t1_1["aa"]("a,i");

        // D_oovv_aaaa += -1.00 P(i,j) P(a,b) r1_aa(a,i) t1_1_aa(b,j) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["234_aa_Lvo"]("L,b,j") * r1["aa"]("R,a,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r1_bb(b,j) t1_1_aa(a,i) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["234_aa_Lvo"]("L,a,i") * r1["bb"]("R,b,j");

        // D_oovv_abab += -1.00 r0 t1_bb(b,j) t1_1_aa(a,i) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["234_aa_Lvo"]("L,a,i") * tmps_["224_bb_Lvo"]("R,b,j");

        // flops: o1v1L1  = o1v1L1
        //  mems: o1v1L1  = o1v1L1
        tmps_["235_aa_Lvo"]("L,a,i")  = l0("L") * t1["aa"]("a,i");

        // D_oovv_aaaa += -1.00 P(i,j) P(a,b) r1_aa(a,i) t1_aa(b,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["235_aa_Lvo"]("L,b,j") * r1["aa"]("R,a,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r1_bb(b,j) t1_aa(a,i) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["235_aa_Lvo"]("L,a,i") * r1["bb"]("R,b,j");

        // D_oovv_abab += -1.00 r0 t1_aa(a,i) t1_bb(b,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["235_aa_Lvo"]("L,a,i") * tmps_["224_bb_Lvo"]("R,b,j");

        // flops: o3v1L1  = o3v2L1 o3v2L1 o4v0L1 o3v2L1 o3v2L1 o4v0L1 o4v0L1 o4v1L1
        //  mems: o3v1L1  = o2v0L1 o2v0L1 o4v0L1 o2v0L1 o2v0L1 o4v0L1 o4v0L1 o3v1L1
        tmps_["236_aaaa_Lvooo"]("L,a,i,j,k")  = t1["aa"]("a,l") * (0.50 * (t2_1["aaaa"]("d,f,k,i") * l2_1["aaaa"]("L,k,l,d,f") + t2["aaaa"]("b,d,o,j") * l2["aaaa"]("L,o,k,b,d")) + l2["abab"]("L,k,m,b,c") * t2["abab"]("b,c,j,m") + l2_1["abab"]("L,l,n,d,e") * t2_1["abab"]("d,e,i,n"));

        // D_ooov_abab += +0.50 d_bb(i,k) r0 t1_aa(a,m) t2_1_abab(c,b,j,l) l2_1_abab(m,l,c,b)
        //             += +0.50 d_bb(i,k) r0 t1_aa(a,m) t2_1_abab(b,c,j,l) l2_1_abab(m,l,b,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_aa(a,l) t1_aa(b,j) t2_abab(d,c,i,k) l2_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_aa(a,l) t1_aa(b,j) t2_abab(c,d,i,k) l2_abab(l,k,c,d)
        //             += +0.50 d_bb(i,k) r0 t1_aa(a,m) t2_1_aaaa(c,b,l,j) l2_1_aaaa(l,m,c,b)
        //             += +0.50 P(i,j) P(a,b) r0 t1_aa(a,l) t1_aa(b,j) t2_aaaa(d,c,k,i) l2_aaaa(k,l,d,c)
        // flops: o3v1L2 += o5v1L2
        //  mems: o3v1L2 += o5v1L2
        D_ooov_abab("L,R,i,j,k,a") += tmps_["236_aaaa_Lvooo"]("L,a,j,i,l") * tmps_["112_bb_Loo"]("R,i,k");

        // flops: o1v1L1  = o1v1L1
        //  mems: o1v1L1  = o1v1L1
        tmps_["237_bb_Lvo"]("L,a,i")  = l0_1("L") * t1_1["bb"]("a,i");

        // D_oovv_bbbb += -1.00 P(i,j) P(a,b) r0 t1_bb(b,j) t1_1_bb(a,i) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["237_bb_Lvo"]("L,a,i") * tmps_["224_bb_Lvo"]("R,b,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r1_aa(a,i) t1_1_bb(b,j) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["237_bb_Lvo"]("L,b,j") * r1["aa"]("R,a,i");

        // D_oovv_bbbb += -1.00 P(i,j) P(a,b) r1_bb(a,i) t1_1_bb(b,j) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["237_bb_Lvo"]("L,b,j") * r1["bb"]("R,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v1L1  = o1v1L1
        //  mems: o1v1L1  = o1v1L1
        tmps_["238_aa_Lvo"]("L,a,i")  = l0_1("L") * t1["aa"]("a,i");

        // D_oovv_aaaa += -1.00 P(i,j) P(a,b) r1_1_aa(a,i) t1_aa(b,j) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["238_aa_Lvo"]("L,b,j") * r1_1["aa"]("R,a,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r1_1_bb(b,j) t1_aa(a,i) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["238_aa_Lvo"]("L,a,i") * r1_1["bb"]("R,b,j");

        // flops: o1v1L1  = o1v1L1
        //  mems: o1v1L1  = o1v1L1
        tmps_["239_bb_Lvo"]("L,a,i")  = l0_1("L") * t1["bb"]("a,i");

        // D_oovv_bbbb += -1.00 P(i,j) r0_1 t1_bb(a,i) t1_bb(b,j) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["239_bb_Lvo"]("L,b,j") * tmps_["225_bb_Lvo"]("R,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r1_1_aa(a,i) t1_bb(b,j) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["239_bb_Lvo"]("L,b,j") * r1_1["aa"]("R,a,i");

        // D_oovv_bbbb += -1.00 P(i,j) P(a,b) r1_1_bb(a,i) t1_bb(b,j) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["239_bb_Lvo"]("L,b,j") * r1_1["bb"]("R,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["240_aaaa_Lvoov"]("L,a,i,j,b")  = t2["aaaa"]("a,c,k,i") * l2_1["aaaa"]("L,k,j,b,c");

        // D_oovv_aaaa += +1.00 P(i,j) r0_1 t2_aaaa(a,c,k,i) t2_aaaa(b,d,l,j) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["240_aaaa_Lvoov"]("L,a,i,l,d") * t2["aaaa"]("b,d,l,j") * r0_1("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r0_1 t2_aaaa(a,c,k,i) t2_abab(d,b,l,j) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["240_aaaa_Lvoov"]("L,a,i,l,d") * t2["abab"]("d,b,l,j") * r0_1("R");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["241_aaaa_Lvoov"]("L,a,i,j,b")  = t2_1["aaaa"]("a,c,k,i") * l2_1["aaaa"]("L,k,j,b,c");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r0 t2_aaaa(b,d,l,j) t2_1_aaaa(a,c,k,i) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["241_aaaa_Lvoov"]("L,a,i,l,d") * t2["aaaa"]("b,d,l,j") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r0 t2_abab(d,b,l,j) t2_1_aaaa(a,c,k,i) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["241_aaaa_Lvoov"]("L,a,i,l,d") * t2["abab"]("d,b,l,j") * r0("R");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["242_aaaa_Lvoov"]("L,a,i,j,b")  = t2["aaaa"]("a,c,k,i") * l2["aaaa"]("L,k,j,b,c");

        // D_oovv_aaaa += +1.00 P(i,j) r0 t2_aaaa(a,c,k,i) t2_aaaa(b,d,l,j) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["242_aaaa_Lvoov"]("L,a,i,l,d") * t2["aaaa"]("b,d,l,j") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r0 t2_aaaa(a,c,k,i) t2_abab(d,b,l,j) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["242_aaaa_Lvoov"]("L,a,i,l,d") * t2["abab"]("d,b,l,j") * r0("R");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["243_aaaa_Lvoov"]("L,a,i,j,b")  = t2["abab"]("a,c,i,k") * l2_1["abab"]("L,j,k,b,c");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r1_aa(d,i) t2_abab(b,c,j,l) t1_1_aa(a,k) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["243_aaaa_Lvoov"]("L,b,j,k,d") * r1["aa"]("R,d,i") * t1_1["aa"]("a,k");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["246_aa_Lvv"]("L,a,b")  = l2_1["aaaa"]("L,i,j,a,c") * t2["aaaa"]("b,c,i,j");

        // D_oovv_aaaa += -0.50 r0_1 t2_aaaa(a,c,k,l) t2_aaaa(b,d,i,j) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") -= 0.50 * tmps_["246_aa_Lvv"]("L,d,a") * t2["aaaa"]("b,d,i,j") * r0_1("R");

        // D_oovv_abab += +0.50 r0_1 t2_aaaa(a,c,k,l) t2_abab(d,b,i,j) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * tmps_["246_aa_Lvv"]("L,d,a") * t2["abab"]("d,b,i,j") * r0_1("R");

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["247_aa_Lvv"]("L,a,b")  = l2_1["aaaa"]("L,i,j,a,c") * t2_1["aaaa"]("b,c,i,j");

        // D_oovv_aaaa += -0.50 P(a,b) r0 t2_aaaa(b,d,i,j) t2_1_aaaa(a,c,k,l) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["247_aa_Lvv"]("L,d,a") * t2["aaaa"]("b,d,i,j") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r0 t2_abab(d,b,i,j) t2_1_aaaa(a,c,k,l) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * tmps_["247_aa_Lvv"]("L,d,a") * t2["abab"]("d,b,i,j") * r0("R");

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["248_aa_Lvv"]("L,a,b")  = t2["aaaa"]("a,c,i,j") * l2["aaaa"]("L,i,j,b,c");

        // D_oovv_aaaa += -0.50 r0 t2_aaaa(a,c,k,l) t2_aaaa(b,d,i,j) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") -= 0.50 * tmps_["248_aa_Lvv"]("L,a,d") * t2["aaaa"]("b,d,i,j") * r0("R");

        // D_oovv_abab += +0.50 r0 t2_aaaa(a,c,k,l) t2_abab(d,b,i,j) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * tmps_["248_aa_Lvv"]("L,a,d") * t2["abab"]("d,b,i,j") * r0("R");

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["249_aa_Lvv"]("L,a,b")  = l2_1["abab"]("L,i,j,a,c") * t2_1["abab"]("b,c,i,j");

        // D_oovv_aaaa += +0.50 P(a,b) r2_aaaa(a,d,i,j) t2_1_abab(b,c,k,l) l2_1_abab(k,l,d,c)
        //             += +0.50 P(a,b) r2_aaaa(a,d,i,j) t2_1_abab(b,c,l,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["249_aa_Lvv"]("L,d,b") * r2["aaaa"]("R,a,d,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += -0.50 P(a,b) r0 t2_aaaa(b,d,i,j) t2_1_abab(a,c,k,l) l2_1_abab(k,l,d,c)
        //             += -0.50 P(a,b) r0 t2_aaaa(b,d,i,j) t2_1_abab(a,c,l,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["249_aa_Lvv"]("L,d,a") * t2["aaaa"]("b,d,i,j") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r2_abab(d,b,i,j) t2_1_abab(a,c,k,l) l2_1_abab(k,l,d,c)
        //             += +0.50 r2_abab(d,b,i,j) t2_1_abab(a,c,l,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["249_aa_Lvv"]("L,d,a") * r2["abab"]("R,d,b,i,j");

        // D_oovv_abab += +0.50 r0 t2_abab(d,b,i,j) t2_1_abab(a,c,k,l) l2_1_abab(k,l,d,c)
        //             += +0.50 r0 t2_abab(d,b,i,j) t2_1_abab(a,c,l,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["249_aa_Lvv"]("L,d,a") * t2["abab"]("d,b,i,j") * r0("R");

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["250_aa_Lvv"]("L,a,b")  = t2["abab"]("a,c,i,j") * l2["abab"]("L,i,j,b,c");

        // D_oovv_aaaa += +0.50 P(a,b) r2_aaaa(a,d,i,j) t2_abab(b,c,k,l) l2_abab(k,l,d,c)
        //             += +0.50 P(a,b) r2_aaaa(a,d,i,j) t2_abab(b,c,l,k) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["250_aa_Lvv"]("L,b,d") * r2["aaaa"]("R,a,d,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r2_abab(d,b,i,j) t2_abab(a,c,k,l) l2_abab(k,l,d,c)
        //             += +0.50 r2_abab(d,b,i,j) t2_abab(a,c,l,k) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["250_aa_Lvv"]("L,a,d") * r2["abab"]("R,d,b,i,j");

        // D_oovv_abab += +0.50 r0 t2_abab(a,c,k,l) t2_abab(d,b,i,j) l2_abab(k,l,d,c)
        //             += +0.50 r0 t2_abab(a,c,l,k) t2_abab(d,b,i,j) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["250_aa_Lvv"]("L,a,d") * t2["abab"]("d,b,i,j") * r0("R");

        // flops: o1v1L1  = o1v1L1
        //  mems: o1v1L1  = o1v1L1
        tmps_["253_aa_Lvo"]("R,a,i")  = r0("R") * t1["aa"]("a,i");

        // D_ooov_aaaa += +0.50 P(i,j) r0 t1_aa(a,j) t2_1_abab(c,b,i,l) l2_1_abab(k,l,c,b)
        //             += +0.50 P(i,j) r0 t1_aa(a,j) t2_1_abab(b,c,i,l) l2_1_abab(k,l,b,c)
        //             += +0.50 r1_bb(b,j) t1_aa(a,l) t2_abab(d,c,i,k) l2_abab(l,k,d,c)
        //             += +0.50 r1_bb(b,j) t1_aa(a,l) t2_abab(c,d,i,k) l2_abab(l,k,c,d)
        // flops: o3v1L2 += o3v2L1 o3v2L1 o3v0L1 o4v1L2
        //  mems: o3v1L2 += o2v0L1 o2v0L1 o3v0L1 o4v1L2
        tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k")  = (t2_1["abab"]("c,b,i,l") * l2_1["abab"]("L,k,l,c,b") + t2["abab"]("d,c,i,k") * l2["abab"]("L,l,k,d,c")) * tmps_["253_aa_Lvo"]("R,a,j");
        D_ooov_aaaa("L,R,i,j,k,a") += tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k");
        D_ooov_aaaa("L,R,i,j,k,a") -= tmps_["perm_aaaa_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_aaaa_LLvooo"].~TArrayD();

        // D_ooov_aaaa += +0.50 P(i,j) r0 t1_aa(a,j) t2_1_aaaa(c,b,l,i) l2_1_aaaa(l,k,c,b)
        //             += +0.50 r0 t1_aa(a,l) t1_bb(b,j) t2_aaaa(d,c,k,i) l2_aaaa(k,l,d,c)
        // flops: o3v1L2 += o3v2L1 o3v2L1 o3v0L1 o4v1L2
        //  mems: o3v1L2 += o2v0L1 o2v0L1 o3v0L1 o4v1L2
        tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k")  = 0.50 * (t2_1["aaaa"]("c,b,l,i") * l2_1["aaaa"]("L,l,k,c,b") + t2["aaaa"]("d,c,k,i") * l2["aaaa"]("L,k,l,d,c")) * tmps_["253_aa_Lvo"]("R,a,j");
        D_ooov_aaaa("L,R,i,j,k,a") += tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k");
        D_ooov_aaaa("L,R,i,j,k,a") -= tmps_["perm_aaaa_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_aaaa_LLvooo"].~TArrayD();

        // D_oovv_abab += +0.50 r0 t1_aa(a,i) t1_bb(b,l) t2_bbbb(d,c,k,j) l2_bbbb(k,l,d,c)
        //             += -0.50 P(i,j) r1_bb(a,i) t2_1_bbbb(c,b,l,j) l2_1_bbbb(l,k,c,b)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o3v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * tmps_["232_bbb_Loov"]("L,j,k,b") * tmps_["253_aa_Lvo"]("R,a,i");
        tmps_["232_bbb_Loov"].~TArrayD();

        // D_oovv_abab += -1.00 r0 t1_aa(a,i) t1_1_bb(b,j) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["237_bb_Lvo"]("L,b,j") * tmps_["253_aa_Lvo"]("R,a,i");
        tmps_["237_bb_Lvo"].~TArrayD();

        // D_oovv_abab += +0.50 r0 t1_aa(a,i) t1_bb(d,j) t2_1_abab(c,b,k,l) l2_1_abab(k,l,c,d)
        //             += +0.50 r0 t1_aa(a,i) t1_bb(d,j) t2_1_abab(c,b,l,k) l2_1_abab(l,k,c,d)
        //             += +0.50 r0 t1_aa(a,i) t2_abab(d,b,k,l) t1_1_bb(c,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r0 t1_aa(a,i) t2_abab(d,b,l,k) t1_1_bb(c,j) l2_1_abab(l,k,d,c)
        //             += +1.00 r0 t1_aa(a,i) t1_bb(c,j) t1_1_bb(b,k) l1_1_bb(k,c)
        //             += -1.00 r0 t1_aa(a,i) t2_abab(c,b,k,j) l1_aa(k,c)
        //             += -1.00 r0 t1_aa(a,i) t2_1_abab(c,b,k,j) l1_1_aa(k,c)
        //             += +1.00 r0 t1_aa(a,i) t1_bb(b,k) t1_bb(c,j) l1_bb(k,c)
        //             += +1.00 r0 t1_aa(a,i) t2_bbbb(b,c,k,j) l1_bb(k,c)
        //             += +1.00 r0 t1_aa(a,i) t1_bb(b,k) t1_1_bb(c,j) l1_1_bb(k,c)
        //             += +1.00 r0 t1_aa(a,i) t2_1_bbbb(b,c,k,j) l1_1_bb(k,c)
        //             += +0.50 r0 t1_aa(a,i) t2_abab(d,c,l,j) t1_1_bb(b,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r0 t1_aa(a,i) t2_abab(c,d,l,j) t1_1_bb(b,k) l2_1_abab(l,k,c,d)
        //             += +0.50 r0 t1_aa(a,i) t2_abab(c,b,k,l) t1_bb(d,j) l2_abab(k,l,c,d)
        //             += +0.50 r0 t1_aa(a,i) t2_abab(c,b,l,k) t1_bb(d,j) l2_abab(l,k,c,d)
        //             += -0.50 r0 t1_aa(a,i) t2_bbbb(b,d,k,l) t1_1_bb(c,j) l2_1_bbbb(k,l,d,c)
        //             += +0.50 r0 t1_aa(a,i) t2_bbbb(b,c,k,l) t1_bb(d,j) l2_bbbb(k,l,d,c)
        //             += -0.50 r0 t1_aa(a,i) t2_bbbb(d,c,l,j) t1_1_bb(b,k) l2_1_bbbb(k,l,d,c)
        //             += +0.50 r0 t1_aa(a,i) t1_bb(d,j) t2_1_bbbb(b,c,k,l) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["223_bb_Lov"]("L,j,b") * tmps_["253_aa_Lvo"]("R,a,i");
        tmps_["223_bb_Lov"].~TArrayD();

        // D_ooov_abab += +0.50 r0 t1_aa(a,j) t2_bbbb(c,b,l,i) l2_bbbb(l,k,c,b)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,i) t1_bb(b,l) t2_1_bbbb(d,c,k,j) l2_1_bbbb(k,l,d,c)
        // flops: o3v1L2 += o5v1L2
        //  mems: o3v1L2 += o5v1L2
        D_ooov_abab("L,R,i,j,k,a") += 0.50 * tmps_["233_bbbb_Loooo"]("L,i,k,j,l") * tmps_["253_aa_Lvo"]("R,a,j");
        tmps_["233_bbbb_Loooo"].~TArrayD();

        // D_oovv_aaaa += -1.00 P(i,j) P(a,b) r0 t1_aa(b,j) t1_1_aa(a,i) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["234_aa_Lvo"]("L,a,i") * tmps_["253_aa_Lvo"]("R,b,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();
        tmps_["234_aa_Lvo"].~TArrayD();

        // D_oovv_aaaa += +0.50 P(i,j) P(a,b) r0 t1_aa(a,l) t1_aa(b,j) t2_1_abab(d,c,i,k) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_aa(a,l) t1_aa(b,j) t2_1_abab(c,d,i,k) l2_1_abab(l,k,c,d)
        //             += +0.50 d_bb(i,k) r0 t1_aa(a,m) t2_abab(c,b,j,l) l2_abab(m,l,c,b)
        //             += +0.50 d_bb(i,k) r0 t1_aa(a,m) t2_abab(b,c,j,l) l2_abab(m,l,b,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_aa(a,l) t1_aa(b,j) t2_1_aaaa(d,c,k,i) l2_1_aaaa(k,l,d,c)
        //             += +0.50 d_bb(i,k) r0 t1_aa(a,m) t2_aaaa(c,b,l,j) l2_aaaa(l,m,c,b)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["236_aaaa_Lvooo"]("L,a,i,j,m") * tmps_["253_aa_Lvo"]("R,b,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();
        tmps_["236_aaaa_Lvooo"].~TArrayD();

        // D_ooov_abab += +0.50 r0 t1_aa(a,j) t2_abab(c,b,l,i) l2_abab(l,k,c,b)
        //             += +0.50 r0 t1_aa(a,j) t2_abab(b,c,l,i) l2_abab(l,k,b,c)
        //             += +0.50 P(i,j) r2_bbbb(b,a,l,i) t2_1_abab(d,c,k,j) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) r2_bbbb(b,a,l,i) t2_1_abab(c,d,k,j) l2_1_abab(k,l,c,d)
        // flops: o3v1L2 += o3v2L1 o3v2L1 o4v0L1 o5v1L2
        //  mems: o3v1L2 += o2v0L1 o2v0L1 o4v0L1 o5v1L2
        D_ooov_abab("L,R,i,j,k,a") += (t2["abab"]("c,b,l,i") * l2["abab"]("L,l,k,c,b") + t2_1["abab"]("d,c,k,j") * l2_1["abab"]("L,k,l,d,c")) * tmps_["253_aa_Lvo"]("R,a,j");

        // D_oovv_aaaa += -1.00 P(i,j) r0 t1_aa(a,i) t1_aa(b,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["235_aa_Lvo"]("L,b,j") * tmps_["253_aa_Lvo"]("R,a,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();
        tmps_["235_aa_Lvo"].~TArrayD();

        // D_oovo_abab += -0.50 r0 t1_aa(a,j) t2_bbbb(c,b,l,i) l2_bbbb(l,k,c,b)
        //             += -0.50 P(i,j) r0 t1_bb(a,j) t2_1_bbbb(c,b,l,i) l2_1_bbbb(l,k,c,b)
        // flops: o3v1L2 += o3v2L1 o3v2L1 o2v0L1 o3v1L2
        //  mems: o3v1L2 += o2v0L1 o2v0L1 o2v0L1 o3v1L2
        D_oovo_abab("L,R,i,j,a,k") -= 0.50 * (t2["bbbb"]("c,b,l,i") * l2["bbbb"]("L,l,k,c,b") + t2_1["bbbb"]("c,b,l,i") * l2_1["bbbb"]("L,l,k,c,b")) * tmps_["253_aa_Lvo"]("R,a,j");

        // D_oovv_abab += +0.50 r0 t1_aa(a,i) t1_bb(b,l) t2_abab(d,c,k,j) l2_abab(k,l,d,c)
        //             += +0.50 r0 t1_aa(a,i) t1_bb(b,l) t2_abab(c,d,k,j) l2_abab(k,l,c,d)
        //             += -0.50 P(i,j) r1_bb(a,i) t2_1_abab(c,b,l,j) l2_1_abab(l,k,c,b)
        //             += -0.50 P(i,j) r1_bb(a,i) t2_1_abab(b,c,l,j) l2_1_abab(l,k,b,c)
        // flops: o2v2L2 += o3v2L1 o3v2L1 o3v0L1 o3v1L1 o3v2L2
        //  mems: o2v2L2 += o2v0L1 o2v0L1 o3v0L1 o2v1L1 o3v2L2
        D_oovv_abab("L,R,i,j,a,b") += (t2["abab"]("d,c,k,j") * l2["abab"]("L,k,l,d,c") + -1.00 * t2_1["abab"]("c,b,l,j") * l2_1["abab"]("L,l,k,c,b")) * t1["bb"]("b,l") * tmps_["253_aa_Lvo"]("R,a,i");

        // D_oovo_abab += -0.50 r0 t1_aa(a,j) t2_abab(c,b,l,i) l2_abab(l,k,c,b)
        //             += -0.50 r0 t1_aa(a,j) t2_abab(b,c,l,i) l2_abab(l,k,b,c)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,i) t1_bb(b,l) t2_1_abab(d,c,k,j) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,i) t1_bb(b,l) t2_1_abab(c,d,k,j) l2_1_abab(k,l,c,d)
        // flops: o3v1L2 += o5v1L2
        //  mems: o3v1L2 += o5v1L2
        D_oovo_abab("L,R,i,j,a,k") -= tmps_["60_bbbb_Loooo"]("L,i,k,j,l") * tmps_["253_aa_Lvo"]("R,a,j");
        tmps_["60_bbbb_Loooo"].~TArrayD();

        // D_oovo_aaaa += -0.50 P(i,j) r0 t1_aa(a,j) t2_1_abab(c,b,i,l) l2_1_abab(k,l,c,b)
        //             += -0.50 P(i,j) r0 t1_aa(a,j) t2_1_abab(b,c,i,l) l2_1_abab(k,l,b,c)
        //             += +0.50 r1_aa(a,l) t1_bb(b,j) t2_abab(d,c,i,k) l2_abab(l,k,d,c)
        //             += +0.50 r1_aa(a,l) t1_bb(b,j) t2_abab(c,d,i,k) l2_abab(l,k,c,d)
        // flops: o3v1L2 += o4v1L2
        //  mems: o3v1L2 += o4v1L2
        tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k")  = tmps_["49_aaa_Looo"]("L,i,k,l") * tmps_["253_aa_Lvo"]("R,a,j");
        D_oovo_aaaa("L,R,i,j,a,k") -= tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k");
        D_oovo_aaaa("L,R,i,j,a,k") += tmps_["perm_aaaa_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_aaaa_LLvooo"].~TArrayD();
        tmps_["49_aaa_Looo"].~TArrayD();

        // D_oovo_aaaa += -0.50 P(i,j) r0 t1_aa(a,j) t2_1_aaaa(c,b,l,i) l2_1_aaaa(l,k,c,b)
        //             += +0.50 r1_bb(b,j) t1_aa(a,l) t2_aaaa(d,c,k,i) l2_aaaa(k,l,d,c)
        // flops: o3v1L2 += o4v1L2
        //  mems: o3v1L2 += o4v1L2
        tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k")  = 0.50 * tmps_["105_aaa_Looo"]("L,i,k,l") * tmps_["253_aa_Lvo"]("R,a,j");
        D_oovo_aaaa("L,R,i,j,a,k") -= tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k");
        D_oovo_aaaa("L,R,i,j,a,k") += tmps_["perm_aaaa_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_aaaa_LLvooo"].~TArrayD();
        tmps_["105_aaa_Looo"].~TArrayD();

        // flops: o1v1L1  = o1v1L1
        //  mems: o1v1L1  = o1v1L1
        tmps_["254_aa_Lvo"]("R,a,i")  = r0_1("R") * t1["aa"]("a,i");

        // D_oovv_abab += +1.00 r0_1 t1_aa(a,i) t2_bbbb(b,c,k,j) l1_1_bb(k,c)
        //             += +1.00 r0_1 t1_aa(a,i) t1_bb(b,k) t1_bb(c,j) l1_1_bb(k,c)
        //             += -1.00 r0_1 t1_aa(a,i) t2_abab(c,b,k,j) l1_1_aa(k,c)
        //             += +0.50 r0_1 t1_aa(a,i) t1_bb(b,l) t2_abab(d,c,k,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r0_1 t1_aa(a,i) t1_bb(b,l) t2_abab(c,d,k,j) l2_1_abab(k,l,c,d)
        //             += +0.50 r0_1 t1_aa(a,i) t1_bb(b,l) t2_bbbb(d,c,k,j) l2_1_bbbb(k,l,d,c)
        //             += +0.50 r0_1 t1_aa(a,i) t2_bbbb(b,c,k,l) t1_bb(d,j) l2_1_bbbb(k,l,d,c)
        //             += +0.50 r0_1 t1_aa(a,i) t2_abab(c,b,k,l) t1_bb(d,j) l2_1_abab(k,l,c,d)
        //             += +0.50 r0_1 t1_aa(a,i) t2_abab(c,b,l,k) t1_bb(d,j) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["220_bb_Lvo"]("L,b,j") * tmps_["254_aa_Lvo"]("R,a,i");

        // D_oovv_aaaa += -1.00 P(i,j) r0_1 t1_aa(a,i) t1_aa(b,j) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["238_aa_Lvo"]("L,b,j") * tmps_["254_aa_Lvo"]("R,a,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();
        tmps_["238_aa_Lvo"].~TArrayD();

        // D_oovv_abab += -1.00 r0_1 t1_aa(a,i) t1_bb(b,j) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["239_bb_Lvo"]("L,b,j") * tmps_["254_aa_Lvo"]("R,a,i");
        tmps_["239_bb_Lvo"].~TArrayD();

        // flops: o1v1L1  = o1v1L1
        //  mems: o1v1L1  = o1v1L1
        tmps_["255_aa_Lvo"]("R,a,i")  = r0("R") * t1_1["aa"]("a,i");

        // D_oovv_abab += +1.00 r0 t2_bbbb(b,c,k,j) t1_1_aa(a,i) l1_1_bb(k,c)
        //             += +1.00 r0 t1_bb(b,k) t1_bb(c,j) t1_1_aa(a,i) l1_1_bb(k,c)
        //             += -1.00 r0 t2_abab(c,b,k,j) t1_1_aa(a,i) l1_1_aa(k,c)
        //             += +0.50 r0 t1_bb(b,l) t2_abab(d,c,k,j) t1_1_aa(a,i) l2_1_abab(k,l,d,c)
        //             += +0.50 r0 t1_bb(b,l) t2_abab(c,d,k,j) t1_1_aa(a,i) l2_1_abab(k,l,c,d)
        //             += +0.50 r0 t1_bb(b,l) t2_bbbb(d,c,k,j) t1_1_aa(a,i) l2_1_bbbb(k,l,d,c)
        //             += +0.50 r0 t2_bbbb(b,c,k,l) t1_bb(d,j) t1_1_aa(a,i) l2_1_bbbb(k,l,d,c)
        //             += +0.50 r0 t2_abab(c,b,k,l) t1_bb(d,j) t1_1_aa(a,i) l2_1_abab(k,l,c,d)
        //             += +0.50 r0 t2_abab(c,b,l,k) t1_bb(d,j) t1_1_aa(a,i) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["220_bb_Lvo"]("L,b,j") * tmps_["255_aa_Lvo"]("R,a,i");
        tmps_["220_bb_Lvo"].~TArrayD();

        // flops: o0v2L2  = o2v3L2 o2v3L2 o0v2L2 o2v3L2 o0v2L2 o2v3L2 o0v2L2
        //  mems: o0v2L2  = o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2
        tmps_["256_aa_LLvv"]("L,R,a,b")  = 0.50 * l2["aaaa"]("L,i,k,a,d") * r2["aaaa"]("R,b,d,i,k");
        tmps_["256_aa_LLvv"]("L,R,a,b") += l2_1["abab"]("L,i,j,a,c") * r2_1["abab"]("R,b,c,i,j");
        tmps_["256_aa_LLvv"]("L,R,a,b") += l2["abab"]("L,i,j,a,c") * r2["abab"]("R,b,c,i,j");
        tmps_["256_aa_LLvv"]("L,R,a,b") += 0.50 * l2_1["aaaa"]("L,i,k,a,d") * r2_1["aaaa"]("R,b,d,i,k");

        // D_oovv_aaaa += -0.50 P(a,b) r2_1_abab(a,d,l,k) t2_aaaa(b,c,i,j) l2_1_abab(l,k,c,d)
        //             += -0.50 P(a,b) r2_1_abab(a,d,k,l) t2_aaaa(b,c,i,j) l2_1_abab(k,l,c,d)
        //             += -0.50 P(a,b) r2_aaaa(a,d,l,k) t2_aaaa(b,c,i,j) l2_aaaa(l,k,c,d)
        //             += -0.50 P(a,b) r2_abab(a,d,l,k) t2_aaaa(b,c,i,j) l2_abab(l,k,c,d)
        //             += -0.50 P(a,b) r2_abab(a,d,k,l) t2_aaaa(b,c,i,j) l2_abab(k,l,c,d)
        //             += -0.50 P(a,b) r2_1_aaaa(a,d,l,k) t2_aaaa(b,c,i,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["256_aa_LLvv"]("L,R,c,a") * t2["aaaa"]("b,c,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r2_1_abab(a,d,l,k) t2_abab(c,b,i,j) l2_1_abab(l,k,c,d)
        //             += +0.50 r2_1_abab(a,d,k,l) t2_abab(c,b,i,j) l2_1_abab(k,l,c,d)
        //             += +0.50 r2_aaaa(a,d,l,k) t2_abab(c,b,i,j) l2_aaaa(l,k,c,d)
        //             += +0.50 r2_abab(a,d,l,k) t2_abab(c,b,i,j) l2_abab(l,k,c,d)
        //             += +0.50 r2_abab(a,d,k,l) t2_abab(c,b,i,j) l2_abab(k,l,c,d)
        //             += +0.50 r2_1_aaaa(a,d,l,k) t2_abab(c,b,i,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["256_aa_LLvv"]("L,R,c,a") * t2["abab"]("c,b,i,j");

        // flops: o0v2L2  = o1v1L2 o1v2L2
        //  mems: o0v2L2  = o1v1L2 o0v2L2
        tmps_["260_aa_LLvv"]("L,R,a,b")  = (tmps_["118_aa_LLov"]("L,R,i,a") + -1.00 * tmps_["117_aa_LLov"]("L,R,i,a")) * t1["aa"]("b,i");

        // flops: o0v2L2  = o1v2L2 o1v2L2 o1v2L2 o1v2L2 o0v2L2 o0v2L2 o0v2L2 o1v2L2 o0v2L2 o1v2L2 o0v2L2
        //  mems: o0v2L2  = o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2
        tmps_["259_aa_LLvv"]("L,R,a,b")  = -1.00 * t1["aa"]("a,i") * tmps_["129_aa_LLov"]("L,R,i,b");
        tmps_["259_aa_LLvv"]("L,R,a,b") -= t1["aa"]("a,i") * tmps_["130_aa_LLov"]("L,R,i,b");
        tmps_["259_aa_LLvv"]("L,R,a,b") += t1["aa"]("a,i") * tmps_["127_aa_LLov"]("L,R,i,b");
        tmps_["259_aa_LLvv"]("L,R,a,b") += t1_1["aa"]("a,i") * tmps_["117_aa_LLov"]("L,R,i,b");
        tmps_["259_aa_LLvv"]("L,R,a,b") += t1["aa"]("a,i") * tmps_["128_aa_LLov"]("L,R,i,b");
        tmps_["259_aa_LLvv"]("L,R,a,b") -= t1_1["aa"]("a,i") * tmps_["118_aa_LLov"]("L,R,i,b");

        // flops: o0v2L1  = o1v2L1 o1v2L1 o0v2L1
        //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
        tmps_["258_aa_Lvv"]("L,a,b")  = l1["aa"]("L,i,a") * t1["aa"]("b,i");
        tmps_["258_aa_Lvv"]("L,a,b") += l1_1["aa"]("L,i,a") * t1_1["aa"]("b,i");

        // flops: o0v2L2  = o1v2L2 o1v2L2 o0v2L2
        //  mems: o0v2L2  = o0v2L2 o0v2L2 o0v2L2
        tmps_["257_aa_LLvv"]("L,R,a,b")  = l1["aa"]("L,i,a") * r1["aa"]("R,b,i");
        tmps_["257_aa_LLvv"]("L,R,a,b") += l1_1["aa"]("L,i,a") * r1_1["aa"]("R,b,i");

        // flops: o0v2L1  = o1v2L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["252_aa_Lvv"]("L,a,b")  = l1_1["aa"]("L,i,a") * t1["aa"]("b,i");

        // flops: o0v2L2  = o1v2L2
        //  mems: o0v2L2  = o0v2L2
        tmps_["251_aa_LLvv"]("L,R,a,b")  = l1_1["aa"]("L,i,a") * r1["aa"]("R,b,i");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["245_aaaa_Lvoov"]("L,a,i,j,b")  = t2["abab"]("a,c,i,k") * l2["abab"]("L,j,k,b,c");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["244_aaaa_Lvoov"]("L,a,i,j,b")  = t2_1["abab"]("a,c,i,k") * l2_1["abab"]("L,j,k,b,c");

        // flops: o1v3L2  = o1v3L2 o1v3L2 o2v3L2 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L2 o2v3L1 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o3v3L1 o2v3L2 o3v3L1 o2v3L2 o3v2L1 o3v2L1 o2v3L2 o3v3L1 o2v3L2 o3v2L1 o3v2L1 o2v3L2 o3v3L1 o2v3L2 o3v2L1 o3v2L1 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o3v2L2 o3v2L2 o2v3L2 o1v3L2 o1v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L2 o1v3L2 o2v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        //  mems: o1v3L2  = o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o2v2L1 o1v3L2 o2v2L1 o1v3L2 o3v1L1 o2v2L1 o1v3L2 o2v2L1 o1v3L2 o3v1L1 o2v2L1 o1v3L2 o2v2L1 o1v3L2 o3v1L1 o2v2L1 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o3v1L2 o2v2L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i")  = t1["aa"]("a,i") * tmps_["259_aa_LLvv"]("L,R,c,b");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += 0.50 * tmps_["246_aa_Lvv"]("L,b,a") * r1_1["aa"]("R,c,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += t1_1["aa"]("a,j") * tmps_["87_aaaa_LLovvo"]("L,R,j,b,c,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") -= t1_1["aa"]("a,i") * tmps_["251_aa_LLvv"]("L,R,b,c");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += 0.50 * tmps_["247_aa_Lvv"]("L,b,a") * r1["aa"]("R,c,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += t1_1["aa"]("a,j") * tmps_["88_aaaa_LLovvo"]("L,R,j,b,c,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += r1["aa"]("R,c,i") * tmps_["249_aa_Lvv"]("L,b,a");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += t1["aa"]("a,l") * tmps_["245_aaaa_Lvoov"]("L,c,i,l,b") * r0("R");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") -= t1["aa"]("a,i") * tmps_["257_aa_LLvv"]("L,R,b,c");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") -= 0.50 * tmps_["246_aa_Lvv"]("L,b,c") * tmps_["254_aa_Lvo"]("R,a,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") -= 0.50 * tmps_["247_aa_Lvv"]("L,b,c") * tmps_["253_aa_Lvo"]("R,a,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += 0.50 * tmps_["248_aa_Lvv"]("L,a,b") * r1["aa"]("R,c,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += t1["aa"]("a,j") * tmps_["89_aaaa_LLovvo"]("L,R,j,b,c,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") -= tmps_["243_aaaa_Lvoov"]("L,a,i,j,b") * t1_1["aa"]("c,j") * r0("R");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") -= tmps_["253_aa_Lvo"]("R,a,i") * tmps_["258_aa_Lvv"]("L,b,c");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += t1["aa"]("a,l") * tmps_["240_aaaa_Lvoov"]("L,c,i,l,b") * r0_1("R");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") -= t2["abab"]("a,d,i,m") * l2["abab"]("L,l,m,b,d") * r1["aa"]("R,c,l");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") -= t2_1["abab"]("a,d,i,m") * l2_1["abab"]("L,l,m,b,d") * r1["aa"]("R,c,l");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += l2["aaaa"]("L,l,j,b,e") * t1["aa"]("e,i") * t1["aa"]("a,j") * r1["aa"]("R,c,l");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += t2_1["aaaa"]("a,e,j,i") * l2_1["aaaa"]("L,l,j,b,e") * r1["aa"]("R,c,l");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += l2_1["aaaa"]("L,l,j,b,e") * t1["aa"]("e,i") * t1_1["aa"]("a,j") * r1["aa"]("R,c,l");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += t2["aaaa"]("a,e,j,i") * l2_1["aaaa"]("L,l,j,b,e") * r1_1["aa"]("R,c,l");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += l2_1["aaaa"]("L,j,l,b,e") * t1["aa"]("e,i") * t1["aa"]("a,l") * t1_1["aa"]("c,j") * r0("R");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += l2_1["aaaa"]("L,j,l,b,e") * r1["aa"]("R,e,i") * t1["aa"]("a,l") * t1_1["aa"]("c,j");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += l2_1["aaaa"]("L,l,j,b,e") * t1_1["aa"]("e,i") * t1["aa"]("a,j") * r1["aa"]("R,c,l");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") -= t2["abab"]("a,d,i,m") * l2_1["abab"]("L,l,m,b,d") * r1_1["aa"]("R,c,l");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += t2["aaaa"]("a,e,j,i") * l2["aaaa"]("L,l,j,b,e") * r1["aa"]("R,c,l");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += l2_1["aaaa"]("L,l,j,b,e") * t1["aa"]("e,i") * t1["aa"]("a,j") * r1_1["aa"]("R,c,l");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += t1["aa"]("a,j") * tmps_["90_aaaa_LLovvo"]("L,R,j,b,c,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") -= tmps_["250_aa_Lvv"]("L,c,b") * tmps_["253_aa_Lvo"]("R,a,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") -= tmps_["252_aa_Lvv"]("L,b,c") * tmps_["254_aa_Lvo"]("R,a,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += 0.50 * tmps_["246_aa_Lvv"]("L,b,a") * tmps_["255_aa_Lvo"]("R,c,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") -= tmps_["29_aa_Lvv"]("L,b,c") * tmps_["254_aa_Lvo"]("R,a,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += r1["aa"]("R,c,i") * tmps_["250_aa_Lvv"]("L,a,b");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += t1["aa"]("a,l") * tmps_["243_aaaa_Lvoov"]("L,c,i,l,b") * r0_1("R");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += r1["aa"]("R,c,i") * tmps_["258_aa_Lvv"]("L,b,a");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += tmps_["52_aaaa_Lvoov"]("L,a,i,j,b") * t1_1["aa"]("c,j") * r0("R");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") -= 0.50 * tmps_["248_aa_Lvv"]("L,c,b") * tmps_["253_aa_Lvo"]("R,a,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") -= tmps_["249_aa_Lvv"]("L,b,c") * tmps_["253_aa_Lvo"]("R,a,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") -= t1_1["aa"]("a,i") * tmps_["64_aa_LLvv"]("L,R,b,c");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += t1["aa"]("a,l") * tmps_["244_aaaa_Lvoov"]("L,c,i,l,b") * r0("R");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += tmps_["29_aa_Lvv"]("L,b,a") * tmps_["255_aa_Lvo"]("R,c,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") -= t1["aa"]("a,i") * tmps_["256_aa_LLvv"]("L,R,b,c");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += t1["aa"]("a,l") * tmps_["242_aaaa_Lvoov"]("L,c,i,l,b") * r0("R");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += t1["aa"]("a,l") * tmps_["241_aaaa_Lvoov"]("L,c,i,l,b") * r0("R");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += tmps_["252_aa_Lvv"]("L,b,a") * tmps_["255_aa_Lvo"]("R,c,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += tmps_["29_aa_Lvv"]("L,b,a") * r1_1["aa"]("R,c,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += tmps_["252_aa_Lvv"]("L,b,a") * r1_1["aa"]("R,c,i");
        tmps_["261_aaaa_LLvvvo"]("L,R,a,b,c,i") += t1_1["aa"]("c,i") * tmps_["260_aa_LLvv"]("L,R,b,a");
        tmps_["90_aaaa_LLovvo"].~TArrayD();
        tmps_["89_aaaa_LLovvo"].~TArrayD();
        tmps_["88_aaaa_LLovvo"].~TArrayD();
        tmps_["87_aaaa_LLovvo"].~TArrayD();
        tmps_["52_aaaa_Lvoov"].~TArrayD();

        // D_ovvv_aaaa += -1.00 P(b,c) r2_abab(b,d,i,k) t1_1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += -0.50 P(b,c) r1_aa(b,i) t2_1_abab(c,d,j,k) l2_1_abab(j,k,a,d)
        //             += -0.50 P(b,c) r1_aa(b,i) t2_1_abab(c,d,k,j) l2_1_abab(k,j,a,d)
        //             += -0.50 P(b,c) r1_aa(b,i) t2_1_aaaa(c,d,j,k) l2_1_aaaa(j,k,a,d)
        //             += +1.00 P(b,c) r1_aa(b,j) t1_1_aa(c,i) l1_1_aa(j,a)
        //             += -1.00 P(b,c) r0_1 t2_abab(b,d,i,j) t1_aa(c,k) l2_1_abab(k,j,a,d)
        //             += -0.50 P(b,c) r1_aa(b,i) t2_abab(c,d,j,k) l2_abab(j,k,a,d)
        //             += -0.50 P(b,c) r1_aa(b,i) t2_abab(c,d,k,j) l2_abab(k,j,a,d)
        //             += +0.50 P(b,c) r0_1 t2_abab(b,d,j,k) t1_aa(c,i) l2_1_abab(j,k,a,d)
        //             += +0.50 P(b,c) r0_1 t2_abab(b,d,k,j) t1_aa(c,i) l2_1_abab(k,j,a,d)
        //             += -0.50 P(b,c) r0 t2_aaaa(c,d,j,k) t1_1_aa(b,i) l2_1_aaaa(j,k,a,d)
        //             += -1.00 P(b,c) r1_aa(b,i) t1_1_aa(c,j) l1_1_aa(j,a)
        //             += -1.00 P(b,c) r1_aa(b,i) t1_aa(c,j) l1_aa(j,a)
        //             += +1.00 P(b,c) r0_1 t1_aa(b,j) t1_aa(c,i) l1_1_aa(j,a)
        //             += +0.50 P(b,c) r0 t2_abab(b,d,j,k) t1_aa(c,i) l2_abab(j,k,a,d)
        //             += +0.50 P(b,c) r0 t2_abab(b,d,k,j) t1_aa(c,i) l2_abab(k,j,a,d)
        //             += -1.00 P(b,c) r2_1_abab(b,d,i,k) t1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += -1.00 P(b,c) r2_abab(b,d,i,k) t1_aa(c,j) l2_abab(j,k,a,d)
        //             += -1.00 P(b,c) r0 t1_aa(c,k) t1_aa(d,i) t1_1_aa(b,j) l2_1_aaaa(j,k,a,d)
        //             += -1.00 P(b,c) r1_1_aa(b,k) t2_aaaa(c,d,j,i) l2_1_aaaa(k,j,a,d)
        //             += -1.00 P(b,c) r1_aa(b,k) t1_aa(d,i) t1_1_aa(c,j) l2_1_aaaa(k,j,a,d)
        //             += -1.00 P(b,c) r1_aa(b,k) t2_1_aaaa(c,d,j,i) l2_1_aaaa(k,j,a,d)
        //             += -1.00 P(b,c) r1_aa(b,k) t1_aa(c,j) t1_aa(d,i) l2_aaaa(k,j,a,d)
        //             += -1.00 P(b,c) r1_aa(d,i) t1_aa(c,k) t1_1_aa(b,j) l2_1_aaaa(j,k,a,d)
        //             += +1.00 P(b,c) r1_aa(b,k) t2_1_abab(c,d,i,j) l2_1_abab(k,j,a,d)
        //             += +1.00 P(b,c) r1_aa(b,k) t2_abab(c,d,i,j) l2_abab(k,j,a,d)
        //             += -1.00 P(b,c) r1_aa(b,k) t1_aa(c,j) t1_1_aa(d,i) l2_1_aaaa(k,j,a,d)
        //             += +1.00 P(b,c) r1_1_aa(b,k) t2_abab(c,d,i,j) l2_1_abab(k,j,a,d)
        //             += -1.00 P(b,c) r1_aa(b,k) t2_aaaa(c,d,j,i) l2_aaaa(k,j,a,d)
        //             += -1.00 P(b,c) r1_1_aa(b,k) t1_aa(c,j) t1_aa(d,i) l2_1_aaaa(k,j,a,d)
        //             += -1.00 P(b,c) r0 t2_aaaa(c,d,k,i) t1_1_aa(b,j) l2_1_aaaa(j,k,a,d)
        //             += -1.00 P(b,c) r0_1 t2_aaaa(b,d,j,i) t1_aa(c,k) l2_1_aaaa(j,k,a,d)
        //             += +1.00 P(b,c) r0 t1_aa(c,i) t1_1_aa(b,j) l1_1_aa(j,a)
        //             += +1.00 P(b,c) r0 t1_aa(b,j) t1_aa(c,i) l1_aa(j,a)
        //             += +1.00 P(b,c) r0 t2_abab(c,d,i,k) t1_1_aa(b,j) l2_1_abab(j,k,a,d)
        //             += -1.00 P(b,c) r2_1_aaaa(b,d,k,i) t1_aa(c,j) l2_1_aaaa(k,j,a,d)
        //             += -1.00 P(b,c) r2_aaaa(b,d,k,i) t1_aa(c,j) l2_aaaa(k,j,a,d)
        //             += +0.50 P(b,c) r0 t2_aaaa(b,d,j,k) t1_aa(c,i) l2_aaaa(j,k,a,d)
        //             += +0.50 P(b,c) r0 t1_aa(c,i) t2_1_abab(b,d,j,k) l2_1_abab(j,k,a,d)
        //             += +0.50 P(b,c) r0 t1_aa(c,i) t2_1_abab(b,d,k,j) l2_1_abab(k,j,a,d)
        //             += +0.50 P(b,c) r2_abab(b,d,k,j) t1_1_aa(c,i) l2_1_abab(k,j,a,d)
        //             += +0.50 P(b,c) r2_abab(b,d,j,k) t1_1_aa(c,i) l2_1_abab(j,k,a,d)
        //             += +0.50 P(b,c) r2_aaaa(b,d,k,j) t1_1_aa(c,i) l2_1_aaaa(k,j,a,d)
        //             += -1.00 P(b,c) r0 t1_aa(c,k) t2_1_abab(b,d,i,j) l2_1_abab(k,j,a,d)
        //             += -0.50 P(b,c) r0 t2_abab(c,d,j,k) t1_1_aa(b,i) l2_1_abab(j,k,a,d)
        //             += -0.50 P(b,c) r0 t2_abab(c,d,k,j) t1_1_aa(b,i) l2_1_abab(k,j,a,d)
        //             += +0.50 P(b,c) r2_1_abab(b,d,k,j) t1_aa(c,i) l2_1_abab(k,j,a,d)
        //             += +0.50 P(b,c) r2_1_abab(b,d,j,k) t1_aa(c,i) l2_1_abab(j,k,a,d)
        //             += +0.50 P(b,c) r2_aaaa(b,d,k,j) t1_aa(c,i) l2_aaaa(k,j,a,d)
        //             += +0.50 P(b,c) r2_abab(b,d,k,j) t1_aa(c,i) l2_abab(k,j,a,d)
        //             += +0.50 P(b,c) r2_abab(b,d,j,k) t1_aa(c,i) l2_abab(j,k,a,d)
        //             += +0.50 P(b,c) r2_1_aaaa(b,d,k,j) t1_aa(c,i) l2_1_aaaa(k,j,a,d)
        //             += -1.00 P(b,c) r0 t2_aaaa(b,d,j,i) t1_aa(c,k) l2_aaaa(j,k,a,d)
        //             += -0.50 P(b,c) r1_aa(b,i) t2_aaaa(c,d,j,k) l2_aaaa(j,k,a,d)
        //             += -1.00 P(b,c) r0 t1_aa(c,k) t2_1_aaaa(b,d,j,i) l2_1_aaaa(j,k,a,d)
        //             += -1.00 P(b,c) r0 t1_aa(c,j) t1_1_aa(b,i) l1_1_aa(j,a)
        //             += +0.50 P(b,c) r0 t1_aa(c,i) t2_1_aaaa(b,d,j,k) l2_1_aaaa(j,k,a,d)
        //             += +0.50 P(b,c) r0_1 t2_aaaa(b,d,j,k) t1_aa(c,i) l2_1_aaaa(j,k,a,d)
        //             += +1.00 P(b,c) r1_1_aa(b,j) t1_aa(c,i) l1_1_aa(j,a)
        //             += +1.00 P(b,c) r1_aa(b,j) t1_aa(c,i) l1_aa(j,a)
        //             += -1.00 P(b,c) r0 t2_abab(b,d,i,j) t1_aa(c,k) l2_abab(k,j,a,d)
        //             += -0.50 P(b,c) r1_1_aa(b,i) t2_abab(c,d,j,k) l2_1_abab(j,k,a,d)
        //             += -0.50 P(b,c) r1_1_aa(b,i) t2_abab(c,d,k,j) l2_1_abab(k,j,a,d)
        //             += -1.00 P(b,c) r2_aaaa(b,d,k,i) t1_1_aa(c,j) l2_1_aaaa(k,j,a,d)
        //             += -0.50 P(b,c) r1_1_aa(b,i) t2_aaaa(c,d,j,k) l2_1_aaaa(j,k,a,d)
        //             += -1.00 P(b,c) r1_1_aa(b,i) t1_aa(c,j) l1_1_aa(j,a)
        //             += -1.00 P(b,c) r1_bb(d,k) t1_aa(c,j) t1_1_aa(b,i) l2_1_abab(j,k,a,d)
        //             += +1.00 P(b,c) r1_aa(d,k) t1_aa(c,j) t1_1_aa(b,i) l2_1_aaaa(k,j,a,d)
        //             += -1.00 P(b,c) r1_aa(d,k) t1_aa(b,j) t1_aa(c,i) l2_aaaa(k,j,a,d)
        //             += -1.00 P(b,c) r1_aa(d,k) t1_aa(c,i) t1_1_aa(b,j) l2_1_aaaa(k,j,a,d)
        //             += +1.00 P(b,c) r1_1_bb(d,k) t1_aa(b,j) t1_aa(c,i) l2_1_abab(j,k,a,d)
        //             += +1.00 P(b,c) r1_bb(d,k) t1_aa(b,j) t1_aa(c,i) l2_abab(j,k,a,d)
        //             += -1.00 P(b,c) r1_1_aa(d,k) t1_aa(b,j) t1_aa(c,i) l2_1_aaaa(k,j,a,d)
        //             += +1.00 P(b,c) r1_bb(d,k) t1_aa(c,i) t1_1_aa(b,j) l2_1_abab(j,k,a,d)
        D_ovvv_aaaa("L,R,i,a,b,c") -= tmps_["261_aaaa_LLvvvo"]("L,R,c,a,b,i");
        D_ovvv_aaaa("L,R,i,a,b,c") += tmps_["261_aaaa_LLvvvo"]("L,R,b,a,c,i");

        // D_vovv_aaaa += +1.00 P(b,c) r2_abab(b,d,i,k) t1_1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += +0.50 P(b,c) r1_aa(b,i) t2_1_abab(c,d,j,k) l2_1_abab(j,k,a,d)
        //             += +0.50 P(b,c) r1_aa(b,i) t2_1_abab(c,d,k,j) l2_1_abab(k,j,a,d)
        //             += +0.50 P(b,c) r1_aa(b,i) t2_1_aaaa(c,d,j,k) l2_1_aaaa(j,k,a,d)
        //             += -1.00 P(b,c) r1_aa(b,j) t1_1_aa(c,i) l1_1_aa(j,a)
        //             += +1.00 P(b,c) r0_1 t2_abab(b,d,i,j) t1_aa(c,k) l2_1_abab(k,j,a,d)
        //             += +0.50 P(b,c) r1_aa(b,i) t2_abab(c,d,j,k) l2_abab(j,k,a,d)
        //             += +0.50 P(b,c) r1_aa(b,i) t2_abab(c,d,k,j) l2_abab(k,j,a,d)
        //             += -0.50 P(b,c) r0_1 t2_abab(b,d,j,k) t1_aa(c,i) l2_1_abab(j,k,a,d)
        //             += -0.50 P(b,c) r0_1 t2_abab(b,d,k,j) t1_aa(c,i) l2_1_abab(k,j,a,d)
        //             += +0.50 P(b,c) r0 t2_aaaa(c,d,j,k) t1_1_aa(b,i) l2_1_aaaa(j,k,a,d)
        //             += +1.00 P(b,c) r1_aa(b,i) t1_1_aa(c,j) l1_1_aa(j,a)
        //             += +1.00 P(b,c) r1_aa(b,i) t1_aa(c,j) l1_aa(j,a)
        //             += -1.00 P(b,c) r0_1 t1_aa(b,j) t1_aa(c,i) l1_1_aa(j,a)
        //             += -0.50 P(b,c) r0 t2_abab(b,d,j,k) t1_aa(c,i) l2_abab(j,k,a,d)
        //             += -0.50 P(b,c) r0 t2_abab(b,d,k,j) t1_aa(c,i) l2_abab(k,j,a,d)
        //             += +1.00 P(b,c) r2_1_abab(b,d,i,k) t1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += +1.00 P(b,c) r2_abab(b,d,i,k) t1_aa(c,j) l2_abab(j,k,a,d)
        //             += +1.00 P(b,c) r0 t1_aa(c,k) t1_aa(d,i) t1_1_aa(b,j) l2_1_aaaa(j,k,a,d)
        //             += +1.00 P(b,c) r1_1_aa(b,k) t2_aaaa(c,d,j,i) l2_1_aaaa(k,j,a,d)
        //             += +1.00 P(b,c) r1_aa(b,k) t1_aa(d,i) t1_1_aa(c,j) l2_1_aaaa(k,j,a,d)
        //             += +1.00 P(b,c) r1_aa(b,k) t2_1_aaaa(c,d,j,i) l2_1_aaaa(k,j,a,d)
        //             += +1.00 P(b,c) r1_aa(b,k) t1_aa(c,j) t1_aa(d,i) l2_aaaa(k,j,a,d)
        //             += +1.00 P(b,c) r1_aa(d,i) t1_aa(c,k) t1_1_aa(b,j) l2_1_aaaa(j,k,a,d)
        //             += -1.00 P(b,c) r1_aa(b,k) t2_1_abab(c,d,i,j) l2_1_abab(k,j,a,d)
        //             += -1.00 P(b,c) r1_aa(b,k) t2_abab(c,d,i,j) l2_abab(k,j,a,d)
        //             += +1.00 P(b,c) r1_aa(b,k) t1_aa(c,j) t1_1_aa(d,i) l2_1_aaaa(k,j,a,d)
        //             += -1.00 P(b,c) r1_1_aa(b,k) t2_abab(c,d,i,j) l2_1_abab(k,j,a,d)
        //             += +1.00 P(b,c) r1_aa(b,k) t2_aaaa(c,d,j,i) l2_aaaa(k,j,a,d)
        //             += +1.00 P(b,c) r1_1_aa(b,k) t1_aa(c,j) t1_aa(d,i) l2_1_aaaa(k,j,a,d)
        //             += +1.00 P(b,c) r0 t2_aaaa(c,d,k,i) t1_1_aa(b,j) l2_1_aaaa(j,k,a,d)
        //             += +1.00 P(b,c) r0_1 t2_aaaa(b,d,j,i) t1_aa(c,k) l2_1_aaaa(j,k,a,d)
        //             += -1.00 P(b,c) r0 t1_aa(c,i) t1_1_aa(b,j) l1_1_aa(j,a)
        //             += -1.00 P(b,c) r0 t1_aa(b,j) t1_aa(c,i) l1_aa(j,a)
        //             += -1.00 P(b,c) r0 t2_abab(c,d,i,k) t1_1_aa(b,j) l2_1_abab(j,k,a,d)
        //             += +1.00 P(b,c) r2_1_aaaa(b,d,k,i) t1_aa(c,j) l2_1_aaaa(k,j,a,d)
        //             += +1.00 P(b,c) r2_aaaa(b,d,k,i) t1_aa(c,j) l2_aaaa(k,j,a,d)
        //             += -0.50 P(b,c) r0 t2_aaaa(b,d,j,k) t1_aa(c,i) l2_aaaa(j,k,a,d)
        //             += -0.50 P(b,c) r0 t1_aa(c,i) t2_1_abab(b,d,j,k) l2_1_abab(j,k,a,d)
        //             += -0.50 P(b,c) r0 t1_aa(c,i) t2_1_abab(b,d,k,j) l2_1_abab(k,j,a,d)
        //             += -0.50 P(b,c) r2_abab(b,d,k,j) t1_1_aa(c,i) l2_1_abab(k,j,a,d)
        //             += -0.50 P(b,c) r2_abab(b,d,j,k) t1_1_aa(c,i) l2_1_abab(j,k,a,d)
        //             += -0.50 P(b,c) r2_aaaa(b,d,k,j) t1_1_aa(c,i) l2_1_aaaa(k,j,a,d)
        //             += +1.00 P(b,c) r0 t1_aa(c,k) t2_1_abab(b,d,i,j) l2_1_abab(k,j,a,d)
        //             += +0.50 P(b,c) r0 t2_abab(c,d,j,k) t1_1_aa(b,i) l2_1_abab(j,k,a,d)
        //             += +0.50 P(b,c) r0 t2_abab(c,d,k,j) t1_1_aa(b,i) l2_1_abab(k,j,a,d)
        //             += -0.50 P(b,c) r2_1_abab(b,d,k,j) t1_aa(c,i) l2_1_abab(k,j,a,d)
        //             += -0.50 P(b,c) r2_1_abab(b,d,j,k) t1_aa(c,i) l2_1_abab(j,k,a,d)
        //             += -0.50 P(b,c) r2_aaaa(b,d,k,j) t1_aa(c,i) l2_aaaa(k,j,a,d)
        //             += -0.50 P(b,c) r2_abab(b,d,k,j) t1_aa(c,i) l2_abab(k,j,a,d)
        //             += -0.50 P(b,c) r2_abab(b,d,j,k) t1_aa(c,i) l2_abab(j,k,a,d)
        //             += -0.50 P(b,c) r2_1_aaaa(b,d,k,j) t1_aa(c,i) l2_1_aaaa(k,j,a,d)
        //             += +1.00 P(b,c) r0 t2_aaaa(b,d,j,i) t1_aa(c,k) l2_aaaa(j,k,a,d)
        //             += +0.50 P(b,c) r1_aa(b,i) t2_aaaa(c,d,j,k) l2_aaaa(j,k,a,d)
        //             += +1.00 P(b,c) r0 t1_aa(c,k) t2_1_aaaa(b,d,j,i) l2_1_aaaa(j,k,a,d)
        //             += +1.00 P(b,c) r0 t1_aa(c,j) t1_1_aa(b,i) l1_1_aa(j,a)
        //             += -0.50 P(b,c) r0 t1_aa(c,i) t2_1_aaaa(b,d,j,k) l2_1_aaaa(j,k,a,d)
        //             += -0.50 P(b,c) r0_1 t2_aaaa(b,d,j,k) t1_aa(c,i) l2_1_aaaa(j,k,a,d)
        //             += -1.00 P(b,c) r1_1_aa(b,j) t1_aa(c,i) l1_1_aa(j,a)
        //             += -1.00 P(b,c) r1_aa(b,j) t1_aa(c,i) l1_aa(j,a)
        //             += +1.00 P(b,c) r0 t2_abab(b,d,i,j) t1_aa(c,k) l2_abab(k,j,a,d)
        //             += +0.50 P(b,c) r1_1_aa(b,i) t2_abab(c,d,j,k) l2_1_abab(j,k,a,d)
        //             += +0.50 P(b,c) r1_1_aa(b,i) t2_abab(c,d,k,j) l2_1_abab(k,j,a,d)
        //             += +1.00 P(b,c) r2_aaaa(b,d,k,i) t1_1_aa(c,j) l2_1_aaaa(k,j,a,d)
        //             += +0.50 P(b,c) r1_1_aa(b,i) t2_aaaa(c,d,j,k) l2_1_aaaa(j,k,a,d)
        //             += +1.00 P(b,c) r1_1_aa(b,i) t1_aa(c,j) l1_1_aa(j,a)
        //             += +1.00 P(b,c) r1_bb(d,k) t1_aa(c,j) t1_1_aa(b,i) l2_1_abab(j,k,a,d)
        //             += -1.00 P(b,c) r1_aa(d,k) t1_aa(c,j) t1_1_aa(b,i) l2_1_aaaa(k,j,a,d)
        //             += +1.00 P(b,c) r1_aa(d,k) t1_aa(b,j) t1_aa(c,i) l2_aaaa(k,j,a,d)
        //             += +1.00 P(b,c) r1_aa(d,k) t1_aa(c,i) t1_1_aa(b,j) l2_1_aaaa(k,j,a,d)
        //             += -1.00 P(b,c) r1_1_bb(d,k) t1_aa(b,j) t1_aa(c,i) l2_1_abab(j,k,a,d)
        //             += -1.00 P(b,c) r1_bb(d,k) t1_aa(b,j) t1_aa(c,i) l2_abab(j,k,a,d)
        //             += +1.00 P(b,c) r1_1_aa(d,k) t1_aa(b,j) t1_aa(c,i) l2_1_aaaa(k,j,a,d)
        //             += -1.00 P(b,c) r1_bb(d,k) t1_aa(c,i) t1_1_aa(b,j) l2_1_abab(j,k,a,d)
        D_vovv_aaaa("L,R,a,i,b,c") += tmps_["261_aaaa_LLvvvo"]("L,R,c,a,b,i");
        D_vovv_aaaa("L,R,a,i,b,c") -= tmps_["261_aaaa_LLvvvo"]("L,R,b,a,c,i");
        tmps_["261_aaaa_LLvvvo"].~TArrayD();

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["262_aaaa_Looov"]("L,i,j,k,a")  = t1["aa"]("b,i") * l2_1["aaaa"]("L,j,k,a,b");

        // flops: o2v2L2  = o2v2L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o3v2L2 o2v2L2 o3v2L2 o3v2L2 o2v2L2 o3v2L1 o3v2L2 o2v2L2 o3v2L1 o3v2L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o3v2L1 o3v2L1 o2v2L2 o2v2L2 o3v2L1 o3v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o3v2L1 o2v2L2 o2v2L2 o3v2L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2
        //  mems: o2v2L2  = o2v2L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o3v1L2 o2v2L2 o2v2L2 o3v1L1 o2v2L2 o2v2L2 o3v1L1 o2v2L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o3v1L1 o2v2L1 o2v2L2 o2v2L2 o3v1L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L1 o2v2L2 o2v2L2 o2v2L1 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b")  = (tmps_["245_aaaa_Lvoov"]("L,a,i,j,b") + tmps_["242_aaaa_Lvoov"]("L,a,i,j,b")) * r0("R");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += r0_1("R") * tmps_["240_aaaa_Lvoov"]("L,a,i,j,b");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += r0("R") * tmps_["241_aaaa_Lvoov"]("L,a,i,j,b");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += r0("R") * tmps_["244_aaaa_Lvoov"]("L,a,i,j,b");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += r0_1("R") * tmps_["243_aaaa_Lvoov"]("L,a,i,j,b");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") -= tmps_["257_aa_LLvv"]("L,R,b,a") * Id["aa_oo"]("i,j");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += l1["aa"]("L,j,b") * r1["aa"]("R,a,i");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += tmps_["262_aaaa_Looov"]("L,i,l,j,b") * r1_1["aa"]("R,a,l");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += l2_1["aaaa"]("L,l,j,b,d") * r1["aa"]("R,d,i") * t1_1["aa"]("a,l");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += t1_1["aa"]("d,i") * l2_1["aaaa"]("L,l,j,b,d") * r1["aa"]("R,a,l");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += t1_1["aa"]("d,i") * l2_1["aaaa"]("L,l,j,b,d") * t1["aa"]("a,l") * r0("R");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += l1_1["aa"]("L,j,b") * r1_1["aa"]("R,a,i");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += t1["aa"]("d,i") * l2["aaaa"]("L,l,j,b,d") * t1["aa"]("a,l") * r0("R");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += t1["aa"]("d,i") * l2["aaaa"]("L,l,j,b,d") * r1["aa"]("R,a,l");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += t1["aa"]("a,i") * tmps_["129_aa_LLov"]("L,R,j,b");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") -= t1["aa"]("a,i") * tmps_["127_aa_LLov"]("L,R,j,b");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") -= t1_1["aa"]("a,i") * tmps_["117_aa_LLov"]("L,R,j,b");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += t1["aa"]("a,i") * tmps_["130_aa_LLov"]("L,R,j,b");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") -= t1["aa"]("a,i") * tmps_["128_aa_LLov"]("L,R,j,b");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += l1["aa"]("L,j,b") * tmps_["253_aa_Lvo"]("R,a,i");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += l1_1["aa"]("L,j,b") * tmps_["255_aa_Lvo"]("R,a,i");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") -= tmps_["256_aa_LLvv"]("L,R,b,a") * Id["aa_oo"]("i,j");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += t1_1["aa"]("a,i") * tmps_["118_aa_LLov"]("L,R,j,b");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += l1_1["aa"]("L,j,b") * tmps_["254_aa_Lvo"]("R,a,i");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += t1_1["aa"]("a,l") * tmps_["262_aaaa_Looov"]("L,i,l,j,b") * r0("R");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += t1["aa"]("a,l") * tmps_["262_aaaa_Looov"]("L,i,l,j,b") * r0_1("R");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") -= tmps_["250_aa_Lvv"]("L,a,b") * tmps_["103_aa_Loo"]("R,i,j");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") -= tmps_["249_aa_Lvv"]("L,b,a") * tmps_["103_aa_Loo"]("R,i,j");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") -= tmps_["252_aa_Lvv"]("L,b,a") * tmps_["136_aa_Loo"]("R,i,j");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") -= 0.50 * tmps_["247_aa_Lvv"]("L,b,a") * tmps_["103_aa_Loo"]("R,i,j");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") -= 0.50 * tmps_["246_aa_Lvv"]("L,b,a") * tmps_["136_aa_Loo"]("R,i,j");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") -= tmps_["258_aa_Lvv"]("L,b,a") * tmps_["103_aa_Loo"]("R,i,j");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") -= tmps_["29_aa_Lvv"]("L,b,a") * tmps_["136_aa_Loo"]("R,i,j");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") -= 0.50 * tmps_["248_aa_Lvv"]("L,a,b") * tmps_["103_aa_Loo"]("R,i,j");
        tmps_["263_aaaa_LLvoov"]("L,R,a,i,j,b") += Id["aa_oo"]("i,j") * tmps_["259_aa_LLvv"]("L,R,a,b");
        tmps_["262_aaaa_Looov"].~TArrayD();

        // D_ovov_aaaa += +1.00 r0 t2_abab(b,c,i,k) l2_abab(j,k,a,c)
        //             += +1.00 r0 t2_aaaa(b,c,k,i) l2_aaaa(k,j,a,c)
        //             += +1.00 r0_1 t2_aaaa(b,c,k,i) l2_1_aaaa(k,j,a,c)
        //             += +1.00 r0 t2_1_aaaa(b,c,k,i) l2_1_aaaa(k,j,a,c)
        //             += +1.00 r0 t2_1_abab(b,c,i,k) l2_1_abab(j,k,a,c)
        //             += +1.00 r0_1 t2_abab(b,c,i,k) l2_1_abab(j,k,a,c)
        //             += -1.00 d_aa(i,j) r1_1_aa(b,k) l1_1_aa(k,a)
        //             += -1.00 d_aa(i,j) r1_aa(b,k) l1_aa(k,a)
        //             += +1.00 r1_aa(b,i) l1_aa(j,a)
        //             += +1.00 r1_1_aa(b,k) t1_aa(c,i) l2_1_aaaa(k,j,a,c)
        //             += +1.00 r1_aa(c,i) t1_1_aa(b,k) l2_1_aaaa(k,j,a,c)
        //             += +1.00 r1_aa(b,k) t1_1_aa(c,i) l2_1_aaaa(k,j,a,c)
        //             += +1.00 r0 t1_aa(b,k) t1_1_aa(c,i) l2_1_aaaa(k,j,a,c)
        //             += +1.00 r1_1_aa(b,i) l1_1_aa(j,a)
        //             += +1.00 r0 t1_aa(b,k) t1_aa(c,i) l2_aaaa(k,j,a,c)
        //             += +1.00 r1_aa(b,k) t1_aa(c,i) l2_aaaa(k,j,a,c)
        //             += +1.00 r1_bb(c,k) t1_aa(b,i) l2_abab(j,k,a,c)
        //             += -1.00 r1_aa(c,k) t1_aa(b,i) l2_aaaa(k,j,a,c)
        //             += -1.00 r1_aa(c,k) t1_1_aa(b,i) l2_1_aaaa(k,j,a,c)
        //             += +1.00 r1_1_bb(c,k) t1_aa(b,i) l2_1_abab(j,k,a,c)
        //             += -1.00 r1_1_aa(c,k) t1_aa(b,i) l2_1_aaaa(k,j,a,c)
        //             += +1.00 r0 t1_aa(b,i) l1_aa(j,a)
        //             += +1.00 r0 t1_1_aa(b,i) l1_1_aa(j,a)
        //             += -0.50 d_aa(i,j) r2_1_abab(b,c,l,k) l2_1_abab(l,k,a,c)
        //             += -0.50 d_aa(i,j) r2_1_abab(b,c,k,l) l2_1_abab(k,l,a,c)
        //             += -0.50 d_aa(i,j) r2_aaaa(b,c,l,k) l2_aaaa(l,k,a,c)
        //             += -0.50 d_aa(i,j) r2_abab(b,c,l,k) l2_abab(l,k,a,c)
        //             += -0.50 d_aa(i,j) r2_abab(b,c,k,l) l2_abab(k,l,a,c)
        //             += -0.50 d_aa(i,j) r2_1_aaaa(b,c,l,k) l2_1_aaaa(l,k,a,c)
        //             += +1.00 r1_bb(c,k) t1_1_aa(b,i) l2_1_abab(j,k,a,c)
        //             += +1.00 r0_1 t1_aa(b,i) l1_1_aa(j,a)
        //             += +1.00 r0 t1_aa(c,i) t1_1_aa(b,k) l2_1_aaaa(k,j,a,c)
        //             += +1.00 r0_1 t1_aa(b,k) t1_aa(c,i) l2_1_aaaa(k,j,a,c)
        //             += -0.50 d_aa(i,j) r0 t2_abab(b,c,k,l) l2_abab(k,l,a,c)
        //             += -0.50 d_aa(i,j) r0 t2_abab(b,c,l,k) l2_abab(l,k,a,c)
        //             += -0.50 d_aa(i,j) r0 t2_1_abab(b,c,k,l) l2_1_abab(k,l,a,c)
        //             += -0.50 d_aa(i,j) r0 t2_1_abab(b,c,l,k) l2_1_abab(l,k,a,c)
        //             += -1.00 d_aa(i,j) r0_1 t1_aa(b,k) l1_1_aa(k,a)
        //             += -0.50 d_aa(i,j) r0 t2_1_aaaa(b,c,k,l) l2_1_aaaa(k,l,a,c)
        //             += -0.50 d_aa(i,j) r0_1 t2_aaaa(b,c,k,l) l2_1_aaaa(k,l,a,c)
        //             += -1.00 d_aa(i,j) r0 t1_1_aa(b,k) l1_1_aa(k,a)
        //             += -1.00 d_aa(i,j) r0 t1_aa(b,k) l1_aa(k,a)
        //             += -0.50 d_aa(i,j) r0_1 t2_abab(b,c,k,l) l2_1_abab(k,l,a,c)
        //             += -0.50 d_aa(i,j) r0_1 t2_abab(b,c,l,k) l2_1_abab(l,k,a,c)
        //             += -0.50 d_aa(i,j) r0 t2_aaaa(b,c,k,l) l2_aaaa(k,l,a,c)
        //             += +1.00 d_aa(i,j) r1_aa(c,l) t1_aa(b,k) l2_aaaa(l,k,a,c)
        //             += +1.00 d_aa(i,j) r1_aa(c,l) t1_1_aa(b,k) l2_1_aaaa(l,k,a,c)
        //             += -1.00 d_aa(i,j) r1_1_bb(c,l) t1_aa(b,k) l2_1_abab(k,l,a,c)
        //             += -1.00 d_aa(i,j) r1_bb(c,l) t1_aa(b,k) l2_abab(k,l,a,c)
        //             += +1.00 d_aa(i,j) r1_1_aa(c,l) t1_aa(b,k) l2_1_aaaa(l,k,a,c)
        //             += -1.00 d_aa(i,j) r1_bb(c,l) t1_1_aa(b,k) l2_1_abab(k,l,a,c)
        D_ovov_aaaa("L,R,i,a,j,b") += tmps_["263_aaaa_LLvoov"]("L,R,b,i,j,a");

        // D_ovvo_aaaa += -1.00 r0 t2_abab(b,c,i,k) l2_abab(j,k,a,c)
        //             += -1.00 r0 t2_aaaa(b,c,k,i) l2_aaaa(k,j,a,c)
        //             += -1.00 r0_1 t2_aaaa(b,c,k,i) l2_1_aaaa(k,j,a,c)
        //             += -1.00 r0 t2_1_aaaa(b,c,k,i) l2_1_aaaa(k,j,a,c)
        //             += -1.00 r0 t2_1_abab(b,c,i,k) l2_1_abab(j,k,a,c)
        //             += -1.00 r0_1 t2_abab(b,c,i,k) l2_1_abab(j,k,a,c)
        //             += +1.00 d_aa(i,j) r1_1_aa(b,k) l1_1_aa(k,a)
        //             += +1.00 d_aa(i,j) r1_aa(b,k) l1_aa(k,a)
        //             += -1.00 r1_aa(b,i) l1_aa(j,a)
        //             += -1.00 r1_1_aa(b,k) t1_aa(c,i) l2_1_aaaa(k,j,a,c)
        //             += -1.00 r1_aa(c,i) t1_1_aa(b,k) l2_1_aaaa(k,j,a,c)
        //             += -1.00 r1_aa(b,k) t1_1_aa(c,i) l2_1_aaaa(k,j,a,c)
        //             += -1.00 r0 t1_aa(b,k) t1_1_aa(c,i) l2_1_aaaa(k,j,a,c)
        //             += -1.00 r1_1_aa(b,i) l1_1_aa(j,a)
        //             += -1.00 r0 t1_aa(b,k) t1_aa(c,i) l2_aaaa(k,j,a,c)
        //             += -1.00 r1_aa(b,k) t1_aa(c,i) l2_aaaa(k,j,a,c)
        //             += -1.00 r1_bb(c,k) t1_aa(b,i) l2_abab(j,k,a,c)
        //             += +1.00 r1_aa(c,k) t1_aa(b,i) l2_aaaa(k,j,a,c)
        //             += +1.00 r1_aa(c,k) t1_1_aa(b,i) l2_1_aaaa(k,j,a,c)
        //             += -1.00 r1_1_bb(c,k) t1_aa(b,i) l2_1_abab(j,k,a,c)
        //             += +1.00 r1_1_aa(c,k) t1_aa(b,i) l2_1_aaaa(k,j,a,c)
        //             += -1.00 r0 t1_aa(b,i) l1_aa(j,a)
        //             += -1.00 r0 t1_1_aa(b,i) l1_1_aa(j,a)
        //             += +0.50 d_aa(i,j) r2_1_abab(b,c,l,k) l2_1_abab(l,k,a,c)
        //             += +0.50 d_aa(i,j) r2_1_abab(b,c,k,l) l2_1_abab(k,l,a,c)
        //             += +0.50 d_aa(i,j) r2_aaaa(b,c,l,k) l2_aaaa(l,k,a,c)
        //             += +0.50 d_aa(i,j) r2_abab(b,c,l,k) l2_abab(l,k,a,c)
        //             += +0.50 d_aa(i,j) r2_abab(b,c,k,l) l2_abab(k,l,a,c)
        //             += +0.50 d_aa(i,j) r2_1_aaaa(b,c,l,k) l2_1_aaaa(l,k,a,c)
        //             += -1.00 r1_bb(c,k) t1_1_aa(b,i) l2_1_abab(j,k,a,c)
        //             += -1.00 r0_1 t1_aa(b,i) l1_1_aa(j,a)
        //             += -1.00 r0 t1_aa(c,i) t1_1_aa(b,k) l2_1_aaaa(k,j,a,c)
        //             += -1.00 r0_1 t1_aa(b,k) t1_aa(c,i) l2_1_aaaa(k,j,a,c)
        //             += +0.50 d_aa(i,j) r0 t2_abab(b,c,k,l) l2_abab(k,l,a,c)
        //             += +0.50 d_aa(i,j) r0 t2_abab(b,c,l,k) l2_abab(l,k,a,c)
        //             += +0.50 d_aa(i,j) r0 t2_1_abab(b,c,k,l) l2_1_abab(k,l,a,c)
        //             += +0.50 d_aa(i,j) r0 t2_1_abab(b,c,l,k) l2_1_abab(l,k,a,c)
        //             += +1.00 d_aa(i,j) r0_1 t1_aa(b,k) l1_1_aa(k,a)
        //             += +0.50 d_aa(i,j) r0 t2_1_aaaa(b,c,k,l) l2_1_aaaa(k,l,a,c)
        //             += +0.50 d_aa(i,j) r0_1 t2_aaaa(b,c,k,l) l2_1_aaaa(k,l,a,c)
        //             += +1.00 d_aa(i,j) r0 t1_1_aa(b,k) l1_1_aa(k,a)
        //             += +1.00 d_aa(i,j) r0 t1_aa(b,k) l1_aa(k,a)
        //             += +0.50 d_aa(i,j) r0_1 t2_abab(b,c,k,l) l2_1_abab(k,l,a,c)
        //             += +0.50 d_aa(i,j) r0_1 t2_abab(b,c,l,k) l2_1_abab(l,k,a,c)
        //             += +0.50 d_aa(i,j) r0 t2_aaaa(b,c,k,l) l2_aaaa(k,l,a,c)
        //             += -1.00 d_aa(i,j) r1_aa(c,l) t1_aa(b,k) l2_aaaa(l,k,a,c)
        //             += -1.00 d_aa(i,j) r1_aa(c,l) t1_1_aa(b,k) l2_1_aaaa(l,k,a,c)
        //             += +1.00 d_aa(i,j) r1_1_bb(c,l) t1_aa(b,k) l2_1_abab(k,l,a,c)
        //             += +1.00 d_aa(i,j) r1_bb(c,l) t1_aa(b,k) l2_abab(k,l,a,c)
        //             += -1.00 d_aa(i,j) r1_1_aa(c,l) t1_aa(b,k) l2_1_aaaa(l,k,a,c)
        //             += +1.00 d_aa(i,j) r1_bb(c,l) t1_1_aa(b,k) l2_1_abab(k,l,a,c)
        D_ovvo_aaaa("L,R,i,a,b,j") -= tmps_["263_aaaa_LLvoov"]("L,R,b,i,j,a");
        tmps_["263_aaaa_LLvoov"].~TArrayD();

        // flops: o4v0L1  = o3v2L1 o4v1L1
        //  mems: o4v0L1  = o3v1L1 o4v0L1
        tmps_["265_aabb_Loooo"]("L,i,j,k,l")  = t1_1["aa"]("b,i") * l2_1["abab"]("L,j,k,b,a") * t1["bb"]("a,l");

        // D_oooo_abab += -1.00 r0 t1_bb(b,j) t1_1_aa(a,i) l2_1_abab(k,l,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") -= tmps_["265_aabb_Loooo"]("L,i,k,l,j") * r0("R");

        // D_oovv_abab += -0.50 r2_abab(a,b,l,k) t1_bb(d,j) t1_1_aa(c,i) l2_1_abab(l,k,c,d)
        //             += -0.50 r2_abab(a,b,k,l) t1_bb(d,j) t1_1_aa(c,i) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["265_aabb_Loooo"]("L,i,l,k,j") * r2["abab"]("R,a,b,l,k");

        // D_oovv_abab += -1.00 r1_bb(b,l) t1_aa(a,k) t1_bb(d,j) t1_1_aa(c,i) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o4v1L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["265_aabb_Loooo"]("L,i,k,l,j") * r1["bb"]("R,b,l") * t1["aa"]("a,k");

        // flops: o4v0L1  = o3v2L1 o4v1L1
        //  mems: o4v0L1  = o3v1L1 o4v0L1
        tmps_["266_aabb_Loooo"]("L,i,j,k,l")  = t1["aa"]("b,i") * l2["abab"]("L,j,k,b,a") * t1["bb"]("a,l");

        // D_oooo_abab += -1.00 r0 t1_aa(a,i) t1_bb(b,j) l2_abab(k,l,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") -= tmps_["266_aabb_Loooo"]("L,i,k,l,j") * r0("R");

        // D_oovv_abab += -1.00 r1_bb(b,l) t1_aa(a,k) t1_aa(c,i) t1_bb(d,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o4v1L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["266_aabb_Loooo"]("L,i,k,l,j") * r1["bb"]("R,b,l") * t1["aa"]("a,k");

        // flops: o4v0L1  = o4v2L1 o4v2L1 o4v0L1
        //  mems: o4v0L1  = o4v0L1 o4v0L1 o4v0L1
        tmps_["267_abab_Loooo"]("L,i,j,k,l")  = t2["abab"]("a,b,i,j") * l2["abab"]("L,k,l,a,b");
        tmps_["267_abab_Loooo"]("L,i,j,k,l") += t2_1["abab"]("a,b,i,j") * l2_1["abab"]("L,k,l,a,b");

        // D_oooo_abab += -0.50 r0 t2_abab(b,a,i,j) l2_abab(k,l,b,a)
        //             += -0.50 r0 t2_abab(a,b,i,j) l2_abab(k,l,a,b)
        //             += -0.50 r0 t2_1_abab(b,a,i,j) l2_1_abab(k,l,b,a)
        //             += -0.50 r0 t2_1_abab(a,b,i,j) l2_1_abab(k,l,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") -= tmps_["267_abab_Loooo"]("L,i,j,k,l") * r0("R");

        // D_oovv_abab += -0.25 r2_abab(a,b,l,k) t2_abab(d,c,i,j) l2_abab(l,k,d,c)
        //             += -0.25 r2_abab(a,b,l,k) t2_abab(c,d,i,j) l2_abab(l,k,c,d)
        //             += -0.25 r2_abab(a,b,k,l) t2_abab(d,c,i,j) l2_abab(k,l,d,c)
        //             += -0.25 r2_abab(a,b,k,l) t2_abab(c,d,i,j) l2_abab(k,l,c,d)
        //             += -0.25 r2_abab(a,b,l,k) t2_1_abab(d,c,i,j) l2_1_abab(l,k,d,c)
        //             += -0.25 r2_abab(a,b,l,k) t2_1_abab(c,d,i,j) l2_1_abab(l,k,c,d)
        //             += -0.25 r2_abab(a,b,k,l) t2_1_abab(d,c,i,j) l2_1_abab(k,l,d,c)
        //             += -0.25 r2_abab(a,b,k,l) t2_1_abab(c,d,i,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["267_abab_Loooo"]("L,i,j,l,k") * r2["abab"]("R,a,b,l,k");

        // D_oovv_abab += -0.25 r0 t2_abab(a,b,k,l) t2_abab(d,c,i,j) l2_abab(k,l,d,c)
        //             += -0.25 r0 t2_abab(a,b,k,l) t2_abab(c,d,i,j) l2_abab(k,l,c,d)
        //             += -0.25 r0 t2_abab(a,b,l,k) t2_abab(d,c,i,j) l2_abab(l,k,d,c)
        //             += -0.25 r0 t2_abab(a,b,l,k) t2_abab(c,d,i,j) l2_abab(l,k,c,d)
        //             += -0.25 r0 t2_abab(a,b,k,l) t2_1_abab(d,c,i,j) l2_1_abab(k,l,d,c)
        //             += -0.25 r0 t2_abab(a,b,k,l) t2_1_abab(c,d,i,j) l2_1_abab(k,l,c,d)
        //             += -0.25 r0 t2_abab(a,b,l,k) t2_1_abab(d,c,i,j) l2_1_abab(l,k,d,c)
        //             += -0.25 r0 t2_abab(a,b,l,k) t2_1_abab(c,d,i,j) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o4v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,b,k,l") * tmps_["267_abab_Loooo"]("L,i,j,k,l") * r0("R");

        // D_oovv_abab += -0.50 r1_bb(b,l) t1_aa(a,k) t2_abab(d,c,i,j) l2_abab(k,l,d,c)
        //             += -0.50 r1_bb(b,l) t1_aa(a,k) t2_abab(c,d,i,j) l2_abab(k,l,c,d)
        //             += -0.50 r1_bb(b,l) t1_aa(a,k) t2_1_abab(d,c,i,j) l2_1_abab(k,l,d,c)
        //             += -0.50 r1_bb(b,l) t1_aa(a,k) t2_1_abab(c,d,i,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o4v1L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["267_abab_Loooo"]("L,i,j,k,l") * r1["bb"]("R,b,l") * t1["aa"]("a,k");

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["264_aabb_Looov"]("L,i,j,k,a")  = t1["aa"]("b,i") * l2_1["abab"]("L,j,k,b,a");

        // flops: o4v0L1  = o4v1L1
        //  mems: o4v0L1  = o4v0L1
        tmps_["268_baab_Loooo"]("L,i,j,k,l")  = tmps_["264_aabb_Looov"]("L,j,k,l,a") * t1_1["bb"]("a,i");

        // D_oooo_abab += -1.00 r0 t1_aa(b,i) t1_1_bb(a,j) l2_1_abab(k,l,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") -= r0("R") * tmps_["268_baab_Loooo"]("L,j,i,k,l");

        // D_oovv_abab += -1.00 r1_bb(b,l) t1_aa(a,k) t1_aa(d,i) t1_1_bb(c,j) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o4v1L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r1["bb"]("R,b,l") * tmps_["268_baab_Loooo"]("L,j,i,k,l") * t1["aa"]("a,k");

        // flops: o3v1L1  = o4v0L1 o4v0L1 o4v0L1 o4v1L1
        //  mems: o3v1L1  = o4v0L1 o4v0L1 o4v0L1 o3v1L1
        tmps_["269_baab_Looov"]("L,i,j,k,a")  = t1["bb"]("a,l") * (tmps_["268_baab_Loooo"]("L,i,j,k,l") + tmps_["265_aabb_Loooo"]("L,j,k,l,i") + tmps_["266_aabb_Loooo"]("L,j,k,l,i") + tmps_["267_abab_Loooo"]("L,j,i,k,l"));

        // D_oovv_abab += -1.00 r1_aa(a,l) t1_bb(b,k) t1_aa(d,i) t1_1_bb(c,j) l2_1_abab(l,k,d,c)
        //             += -1.00 r1_aa(a,l) t1_bb(b,k) t1_aa(c,i) t1_bb(d,j) l2_abab(l,k,c,d)
        //             += -0.50 r1_aa(a,l) t1_bb(b,k) t2_abab(d,c,i,j) l2_abab(l,k,d,c)
        //             += -0.50 r1_aa(a,l) t1_bb(b,k) t2_abab(c,d,i,j) l2_abab(l,k,c,d)
        //             += -0.50 r1_aa(a,l) t1_bb(b,k) t2_1_abab(d,c,i,j) l2_1_abab(l,k,d,c)
        //             += -0.50 r1_aa(a,l) t1_bb(b,k) t2_1_abab(c,d,i,j) l2_1_abab(l,k,c,d)
        //             += -1.00 r1_aa(a,l) t1_bb(b,k) t1_bb(d,j) t1_1_aa(c,i) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r1["aa"]("R,a,l") * tmps_["269_baab_Looov"]("L,j,i,l,b");

        // D_oovv_abab += -1.00 r0 t1_aa(a,k) t1_bb(b,l) t1_aa(d,i) t1_1_bb(c,j) l2_1_abab(k,l,d,c)
        //             += -1.00 r0 t1_aa(a,k) t1_bb(b,l) t1_aa(c,i) t1_bb(d,j) l2_abab(k,l,c,d)
        //             += -0.50 r0 t1_aa(a,k) t1_bb(b,l) t2_abab(d,c,i,j) l2_abab(k,l,d,c)
        //             += -0.50 r0 t1_aa(a,k) t1_bb(b,l) t2_abab(c,d,i,j) l2_abab(k,l,c,d)
        //             += -0.50 r0 t1_aa(a,k) t1_bb(b,l) t2_1_abab(d,c,i,j) l2_1_abab(k,l,d,c)
        //             += -0.50 r0 t1_aa(a,k) t1_bb(b,l) t2_1_abab(c,d,i,j) l2_1_abab(k,l,c,d)
        //             += -1.00 r0 t1_aa(a,k) t1_bb(b,l) t1_bb(d,j) t1_1_aa(c,i) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["aa"]("a,k") * tmps_["269_baab_Looov"]("L,j,i,k,b") * r0("R");
        tmps_["269_baab_Looov"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["270_baab_Lvoov"]("L,a,i,j,b")  = t2["abab"]("c,a,i,k") * l2["abab"]("L,j,k,c,b");

        // D_oovv_abab += -1.00 r2_abab(a,d,l,j) t2_abab(c,b,i,k) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["270_baab_Lvoov"]("L,b,i,l,d") * r2["abab"]("R,a,d,l,j");

        // D_oovv_abab += -1.00 r0 t2_abab(a,c,k,j) t2_abab(d,b,i,l) l2_abab(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,c,k,j") * tmps_["270_baab_Lvoov"]("L,b,i,k,c") * r0("R");

        // D_oovv_abab += -1.00 r1_bb(d,j) t1_aa(a,l) t2_abab(c,b,i,k) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["270_baab_Lvoov"]("L,b,i,l,d") * r1["bb"]("R,d,j") * t1["aa"]("a,l");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["271_baab_Lvoov"]("L,a,i,j,b")  = t2_1["abab"]("c,a,i,k") * l2_1["abab"]("L,j,k,c,b");

        // D_oovv_abab += -1.00 r2_abab(a,d,l,j) t2_1_abab(c,b,i,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["271_baab_Lvoov"]("L,b,i,l,d") * r2["abab"]("R,a,d,l,j");

        // D_oovv_abab += -1.00 r0 t2_abab(a,d,l,j) t2_1_abab(c,b,i,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,d,l,j") * tmps_["271_baab_Lvoov"]("L,b,i,l,d") * r0("R");

        // D_oovv_abab += -1.00 r1_bb(d,j) t1_aa(a,l) t2_1_abab(c,b,i,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["271_baab_Lvoov"]("L,b,i,l,d") * r1["bb"]("R,d,j") * t1["aa"]("a,l");

        // flops: o4v0L1  = o4v2L1
        //  mems: o4v0L1  = o4v0L1
        tmps_["272_abab_Loooo"]("L,i,j,k,l")  = t2["abab"]("a,b,i,j") * l2_1["abab"]("L,k,l,a,b");

        // D_oooo_abab += -0.50 r0_1 t2_abab(b,a,i,j) l2_1_abab(k,l,b,a)
        //             += -0.50 r0_1 t2_abab(a,b,i,j) l2_1_abab(k,l,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") -= tmps_["272_abab_Loooo"]("L,i,j,k,l") * r0_1("R");

        // D_oovv_abab += -0.25 r2_1_abab(a,b,l,k) t2_abab(d,c,i,j) l2_1_abab(l,k,d,c)
        //             += -0.25 r2_1_abab(a,b,l,k) t2_abab(c,d,i,j) l2_1_abab(l,k,c,d)
        //             += -0.25 r2_1_abab(a,b,k,l) t2_abab(d,c,i,j) l2_1_abab(k,l,d,c)
        //             += -0.25 r2_1_abab(a,b,k,l) t2_abab(c,d,i,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["272_abab_Loooo"]("L,i,j,l,k") * r2_1["abab"]("R,a,b,l,k");

        // D_oovv_abab += -0.25 r0 t2_abab(d,c,i,j) t2_1_abab(a,b,k,l) l2_1_abab(k,l,d,c)
        //             += -0.25 r0 t2_abab(d,c,i,j) t2_1_abab(a,b,l,k) l2_1_abab(l,k,d,c)
        //             += -0.25 r0 t2_abab(c,d,i,j) t2_1_abab(a,b,k,l) l2_1_abab(k,l,c,d)
        //             += -0.25 r0 t2_abab(c,d,i,j) t2_1_abab(a,b,l,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o4v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2_1["abab"]("a,b,k,l") * tmps_["272_abab_Loooo"]("L,i,j,k,l") * r0("R");

        // D_oovv_abab += -0.25 r0_1 t2_abab(a,b,k,l) t2_abab(d,c,i,j) l2_1_abab(k,l,d,c)
        //             += -0.25 r0_1 t2_abab(a,b,k,l) t2_abab(c,d,i,j) l2_1_abab(k,l,c,d)
        //             += -0.25 r0_1 t2_abab(a,b,l,k) t2_abab(d,c,i,j) l2_1_abab(l,k,d,c)
        //             += -0.25 r0_1 t2_abab(a,b,l,k) t2_abab(c,d,i,j) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o4v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,b,k,l") * tmps_["272_abab_Loooo"]("L,i,j,k,l") * r0_1("R");

        // D_oovv_abab += -0.50 r1_bb(b,l) t2_abab(d,c,i,j) t1_1_aa(a,k) l2_1_abab(k,l,d,c)
        //             += -0.50 r1_bb(b,l) t2_abab(c,d,i,j) t1_1_aa(a,k) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o4v1L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["272_abab_Loooo"]("L,i,j,k,l") * r1["bb"]("R,b,l") * t1_1["aa"]("a,k");

        // D_oovv_abab += -0.50 r1_1_bb(b,l) t1_aa(a,k) t2_abab(d,c,i,j) l2_1_abab(k,l,d,c)
        //             += -0.50 r1_1_bb(b,l) t1_aa(a,k) t2_abab(c,d,i,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o4v1L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["272_abab_Loooo"]("L,i,j,k,l") * r1_1["bb"]("R,b,l") * t1["aa"]("a,k");

        // flops: o4v0L1  = o4v1L1
        //  mems: o4v0L1  = o4v0L1
        tmps_["273_baab_Loooo"]("L,i,j,k,l")  = tmps_["264_aabb_Looov"]("L,j,k,l,a") * t1["bb"]("a,i");

        // D_oooo_abab += -1.00 r0_1 t1_aa(a,i) t1_bb(b,j) l2_1_abab(k,l,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") -= r0_1("R") * tmps_["273_baab_Loooo"]("L,j,i,k,l");

        // D_oovv_abab += -1.00 r1_bb(b,l) t1_aa(c,i) t1_bb(d,j) t1_1_aa(a,k) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o4v1L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r1["bb"]("R,b,l") * tmps_["273_baab_Loooo"]("L,j,i,k,l") * t1_1["aa"]("a,k");

        // D_oovv_abab += -1.00 r1_1_bb(b,l) t1_aa(a,k) t1_aa(c,i) t1_bb(d,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o4v1L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r1_1["bb"]("R,b,l") * tmps_["273_baab_Loooo"]("L,j,i,k,l") * t1["aa"]("a,k");

        // flops: o3v1L1  = o3v2L1 o4v1L1 o3v2L1 o3v1L1 o3v1L1 o3v2L1 o3v1L1 o3v2L1 o3v1L1 o3v2L1 o3v1L1 o4v1L1 o3v1L1
        //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1
        tmps_["274_bbaa_Lvooo"]("L,a,i,j,k")  = t1["aa"]("e,j") * tmps_["209_bbaa_Lvoov"]("L,a,i,k,e");
        tmps_["274_bbaa_Lvooo"]("L,a,i,j,k") += t1_1["bb"]("a,l") * tmps_["272_abab_Loooo"]("L,j,i,k,l");
        tmps_["274_bbaa_Lvooo"]("L,a,i,j,k") += t1["bb"]("b,i") * tmps_["270_baab_Lvoov"]("L,a,j,k,b");
        tmps_["274_bbaa_Lvooo"]("L,a,i,j,k") += t1_1["aa"]("d,j") * tmps_["93_bbaa_Lvoov"]("L,a,i,k,d");
        tmps_["274_bbaa_Lvooo"]("L,a,i,j,k") += t1["bb"]("b,i") * tmps_["271_baab_Lvoov"]("L,a,j,k,b");
        tmps_["274_bbaa_Lvooo"]("L,a,i,j,k") += t1_1["bb"]("c,i") * tmps_["94_baab_Lvoov"]("L,a,j,k,c");
        tmps_["274_bbaa_Lvooo"]("L,a,i,j,k") += t1_1["bb"]("a,l") * tmps_["273_baab_Loooo"]("L,i,j,k,l");
        tmps_["271_baab_Lvoov"].~TArrayD();
        tmps_["270_baab_Lvoov"].~TArrayD();
        tmps_["94_baab_Lvoov"].~TArrayD();

        // D_oovv_abab += -1.00 r1_aa(a,l) t1_aa(c,i) t1_bb(d,j) t1_1_bb(b,k) l2_1_abab(l,k,c,d)
        //             += -1.00 r1_aa(a,l) t2_abab(c,b,i,k) t1_bb(d,j) l2_abab(l,k,c,d)
        //             += -0.50 r1_aa(a,l) t2_abab(d,c,i,j) t1_1_bb(b,k) l2_1_abab(l,k,d,c)
        //             += -0.50 r1_aa(a,l) t2_abab(c,d,i,j) t1_1_bb(b,k) l2_1_abab(l,k,c,d)
        //             += -1.00 r1_aa(a,l) t2_bbbb(b,c,k,j) t1_aa(d,i) l2_abab(l,k,d,c)
        //             += -1.00 r1_aa(a,l) t1_aa(d,i) t2_1_bbbb(b,c,k,j) l2_1_abab(l,k,d,c)
        //             += -1.00 r1_aa(a,l) t2_bbbb(b,d,k,j) t1_1_aa(c,i) l2_1_abab(l,k,c,d)
        //             += -1.00 r1_aa(a,l) t1_bb(d,j) t2_1_abab(c,b,i,k) l2_1_abab(l,k,c,d)
        //             += -1.00 r1_aa(a,l) t2_abab(d,b,i,k) t1_1_bb(c,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r1["aa"]("R,a,l") * tmps_["274_bbaa_Lvooo"]("L,b,j,i,l");

        // D_oovv_abab += -1.00 r0 t1_aa(a,l) t1_aa(c,i) t1_bb(d,j) t1_1_bb(b,k) l2_1_abab(l,k,c,d)
        //             += -1.00 r0 t1_aa(a,l) t2_abab(c,b,i,k) t1_bb(d,j) l2_abab(l,k,c,d)
        //             += -0.50 r0 t1_aa(a,l) t2_abab(d,c,i,j) t1_1_bb(b,k) l2_1_abab(l,k,d,c)
        //             += -0.50 r0 t1_aa(a,l) t2_abab(c,d,i,j) t1_1_bb(b,k) l2_1_abab(l,k,c,d)
        //             += -1.00 r0 t1_aa(a,l) t2_bbbb(b,c,k,j) t1_aa(d,i) l2_abab(l,k,d,c)
        //             += -1.00 r0 t1_aa(a,l) t1_aa(d,i) t2_1_bbbb(b,c,k,j) l2_1_abab(l,k,d,c)
        //             += -1.00 r0 t1_aa(a,l) t2_bbbb(b,d,k,j) t1_1_aa(c,i) l2_1_abab(l,k,c,d)
        //             += -1.00 r0 t1_aa(a,l) t1_bb(d,j) t2_1_abab(c,b,i,k) l2_1_abab(l,k,c,d)
        //             += -1.00 r0 t1_aa(a,l) t2_abab(d,b,i,k) t1_1_bb(c,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["aa"]("a,l") * tmps_["274_bbaa_Lvooo"]("L,b,j,i,l") * r0("R");
        tmps_["274_bbaa_Lvooo"].~TArrayD();

        // flops: o0v2L2  = o1v1L2 o1v2L2 o1v2L2 o0v2L2 o1v2L2 o0v2L2 o1v2L2 o0v2L2 o1v2L2 o0v2L2
        //  mems: o0v2L2  = o1v1L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2 o0v2L2
        tmps_["278_bb_LLvv"]("L,R,a,b")  = (tmps_["72_bb_LLov"]("L,R,i,a") + -1.00 * tmps_["69_bb_LLov"]("L,R,i,a")) * t1_1["bb"]("b,i");
        tmps_["278_bb_LLvv"]("L,R,a,b") += t1["bb"]("b,i") * tmps_["73_bb_LLov"]("L,R,i,a");
        tmps_["278_bb_LLvv"]("L,R,a,b") -= t1["bb"]("b,i") * tmps_["71_bb_LLov"]("L,R,i,a");
        tmps_["278_bb_LLvv"]("L,R,a,b") += t1["bb"]("b,i") * tmps_["74_bb_LLov"]("L,R,i,a");
        tmps_["278_bb_LLvv"]("L,R,a,b") -= t1["bb"]("b,i") * tmps_["70_bb_LLov"]("L,R,i,a");

        // flops: o0v2L1  = o1v2L1 o1v2L1 o0v2L1
        //  mems: o0v2L1  = o0v2L1 o0v2L1 o0v2L1
        tmps_["277_bb_Lvv"]("L,a,b")  = l1["bb"]("L,i,a") * t1["bb"]("b,i");
        tmps_["277_bb_Lvv"]("L,a,b") += l1_1["bb"]("L,i,a") * t1_1["bb"]("b,i");

        // flops: o0v2L2  = o1v2L2 o1v2L2 o0v2L2
        //  mems: o0v2L2  = o0v2L2 o0v2L2 o0v2L2
        tmps_["276_bb_LLvv"]("L,R,a,b")  = l1["bb"]("L,i,a") * r1["bb"]("R,b,i");
        tmps_["276_bb_LLvv"]("L,R,a,b") += l1_1["bb"]("L,i,a") * r1_1["bb"]("R,b,i");

        // flops: o0v2L1  = o1v2L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["275_bb_Lvv"]("L,a,b")  = l1_1["bb"]("L,i,a") * t1["bb"]("b,i");

        // flops: o1v3L2  = o1v2L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L2 o3v2L1 o3v2L1 o2v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L1 o1v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o3v2L2 o3v2L2 o2v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o1v2L2 o1v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o1v3L2 o1v2L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        //  mems: o1v3L2  = o0v2L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L2 o3v1L1 o2v2L1 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L1 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o3v1L2 o2v2L2 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o0v2L2 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o0v2L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i")  = -1.00 * tmps_["72_bb_LLov"]("L,R,l,a") * t1["bb"]("c,l") * t1_1["bb"]("b,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += 0.50 * tmps_["219_bb_Lvv"]("L,c,a") * r1_1["bb"]("R,b,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["275_bb_Lvv"]("L,a,c") * tmps_["226_bb_Lvo"]("R,b,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["81_bbbb_Lvoov"]("L,b,i,k,a") * t1["bb"]("c,k") * r0_1("R");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["84_bb_Lvv"]("L,b,a") * tmps_["225_bb_Lvo"]("R,c,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["277_bb_Lvv"]("L,a,b") * tmps_["224_bb_Lvo"]("R,c,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["82_bbbb_Lvoov"]("L,b,i,k,a") * t1["bb"]("c,k") * r0("R");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["208_bbbb_Lvoov"]("L,b,i,k,a") * t1["bb"]("c,k") * r0("R");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += l2_1["bbbb"]("L,k,l,a,e") * t1["bb"]("e,i") * t1["bb"]("c,l") * r1_1["bb"]("R,b,k");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += l2_1["bbbb"]("L,k,l,a,e") * t1_1["bb"]("e,i") * t1["bb"]("c,l") * r1["bb"]("R,b,k");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += l2_1["bbbb"]("L,k,l,a,e") * t2["bbbb"]("c,e,l,i") * r1_1["bb"]("R,b,k");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += l2_1["bbbb"]("L,k,l,a,e") * t1["bb"]("e,i") * t1_1["bb"]("c,l") * r1["bb"]("R,b,k");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") -= l2["abab"]("L,j,k,d,a") * t2["abab"]("d,c,j,i") * r1["bb"]("R,b,k");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += l2_1["bbbb"]("L,k,l,a,e") * t2_1["bbbb"]("c,e,l,i") * r1["bb"]("R,b,k");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += l2["bbbb"]("L,k,l,a,e") * t1["bb"]("e,i") * t1["bb"]("c,l") * r1["bb"]("R,b,k");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += l2_1["bbbb"]("L,l,k,a,e") * t1["bb"]("e,i") * t1["bb"]("c,k") * t1_1["bb"]("b,l") * r0("R");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,j,k,d,a") * t2["abab"]("d,c,j,i") * r1_1["bb"]("R,b,k");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += l2_1["bbbb"]("L,l,k,a,e") * r1["bb"]("R,e,i") * t1["bb"]("c,k") * t1_1["bb"]("b,l");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += l2["bbbb"]("L,k,l,a,e") * t2["bbbb"]("c,e,l,i") * r1["bb"]("R,b,k");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") -= l1_1["bb"]("L,l,a") * r1["bb"]("R,b,l") * t1_1["bb"]("c,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,j,k,d,a") * t2_1["abab"]("d,c,j,i") * r1["bb"]("R,b,k");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["277_bb_Lvv"]("L,a,c") * r1["bb"]("R,b,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["205_bbbb_Lvoov"]("L,b,i,k,a") * t1["bb"]("c,k") * r0_1("R");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["217_bb_Lvv"]("L,a,b") * tmps_["224_bb_Lvo"]("R,c,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["98_bb_LLvv"]("L,R,a,b") * t1_1["bb"]("c,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["214_bbbb_Lvoov"]("L,b,i,k,a") * t1["bb"]("c,k") * r0("R");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["275_bb_Lvv"]("L,a,b") * tmps_["225_bb_Lvo"]("R,c,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") -= 0.50 * tmps_["219_bb_Lvv"]("L,b,a") * tmps_["225_bb_Lvo"]("R,c,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") -= 0.50 * tmps_["221_bb_Lvv"]("L,a,b") * tmps_["224_bb_Lvo"]("R,c,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += t1_1["bb"]("b,l") * tmps_["58_bbbb_Lvoov"]("L,c,i,l,a") * r0("R");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") -= 0.50 * tmps_["210_bb_LLvv"]("L,R,a,b") * t1["bb"]("c,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["84_bb_Lvv"]("L,c,a") * r1_1["bb"]("R,b,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += 0.50 * tmps_["222_bb_Lvv"]("L,c,a") * r1["bb"]("R,b,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") -= 0.50 * tmps_["222_bb_Lvv"]("L,b,a") * tmps_["224_bb_Lvo"]("R,c,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += 0.50 * tmps_["219_bb_Lvv"]("L,c,a") * tmps_["226_bb_Lvo"]("R,b,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["276_bb_LLvv"]("L,R,a,b") * t1["bb"]("c,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["85_bb_Lvv"]("L,a,b") * tmps_["224_bb_Lvo"]("R,c,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["77_bbbb_LLovvo"]("L,R,l,a,b,i") * t1_1["bb"]("c,l");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["84_bb_Lvv"]("L,c,a") * tmps_["226_bb_Lvo"]("R,b,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["213_bbbb_Lvoov"]("L,b,i,k,a") * t1["bb"]("c,k") * r0("R");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") -= t1_1["bb"]("b,l") * tmps_["81_bbbb_Lvoov"]("L,c,i,l,a") * r0("R");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["78_bbbb_LLovvo"]("L,R,l,a,b,i") * t1["bb"]("c,l");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += 0.50 * tmps_["221_bb_Lvv"]("L,a,c") * r1["bb"]("R,b,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["69_bb_LLov"]("L,R,l,a") * t1["bb"]("c,l") * t1_1["bb"]("b,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["275_bb_Lvv"]("L,a,c") * r1_1["bb"]("R,b,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["217_bb_Lvv"]("L,a,c") * r1["bb"]("R,b,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["85_bb_Lvv"]("L,a,c") * r1["bb"]("R,b,i");
        tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i") += t1["bb"]("c,i") * tmps_["278_bb_LLvv"]("L,R,a,b");
        tmps_["98_bb_LLvv"].~TArrayD();
        tmps_["78_bbbb_LLovvo"].~TArrayD();
        tmps_["77_bbbb_LLovvo"].~TArrayD();
        tmps_["58_bbbb_Lvoov"].~TArrayD();
    }
} // hilbert
#endif