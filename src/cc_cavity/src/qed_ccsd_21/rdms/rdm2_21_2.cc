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


    void EOM_EE_QED_RDM_21::rdm2_21_2() {

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


        // D_vvvv_bbbb += +0.50 r2_bbbb(d,c,j,i) l2_bbbb(j,i,a,b)
        //             += +0.50 r2_1_bbbb(d,c,j,i) l2_1_bbbb(j,i,a,b)
        D_vvvv_bbbb("L,R,a,b,c,d") += 0.50 * tmps_["39_bbbb_LLvvvv"]("L,R,a,b,d,c");

        // flops: o1v3L2  = o1v4L2
        //  mems: o1v3L2  = o1v3L2
        tmps_["40_bbbb_LLvvvo"]("L,R,a,b,c,i")  = tmps_["39_bbbb_LLvvvv"]("L,R,a,d,b,c") * t1["bb"]("d,i");
        tmps_["39_bbbb_LLvvvv"].~TArrayD();

        // D_oovv_bbbb += -0.50 r2_bbbb(b,a,l,k) t1_bb(c,i) t1_bb(d,j) l2_bbbb(l,k,d,c)
        //             += -0.50 r2_1_bbbb(b,a,l,k) t1_bb(c,i) t1_bb(d,j) l2_1_bbbb(l,k,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") -= 0.50 * tmps_["40_bbbb_LLvvvo"]("L,R,d,b,a,i") * t1["bb"]("d,j");

        // D_ovvv_bbbb += +0.50 r2_bbbb(c,b,k,j) t1_bb(d,i) l2_bbbb(k,j,a,d)
        //             += +0.50 r2_1_bbbb(c,b,k,j) t1_bb(d,i) l2_1_bbbb(k,j,a,d)
        D_ovvv_bbbb("L,R,i,a,b,c") += 0.50 * tmps_["40_bbbb_LLvvvo"]("L,R,a,c,b,i");

        // D_vovv_bbbb += -0.50 r2_bbbb(c,b,k,j) t1_bb(d,i) l2_bbbb(k,j,a,d)
        //             += -0.50 r2_1_bbbb(c,b,k,j) t1_bb(d,i) l2_1_bbbb(k,j,a,d)
        D_vovv_bbbb("L,R,a,i,b,c") -= 0.50 * tmps_["40_bbbb_LLvvvo"]("L,R,a,c,b,i");
        tmps_["40_bbbb_LLvvvo"].~TArrayD();

        // flops: o4v0L2  = o4v2L2
        //  mems: o4v0L2  = o4v0L2
        tmps_["41_aaaa_LLoooo"]("L,R,i,j,k,l")  = l2_1["aaaa"]("L,i,j,a,b") * r2["aaaa"]("R,a,b,k,l");

        // D_oovv_aaaa += +0.25 r2_aaaa(c,d,i,j) t2_1_aaaa(b,a,k,l) l2_1_aaaa(k,l,c,d)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") += 0.25 * t2_1["aaaa"]("b,a,k,l") * tmps_["41_aaaa_LLoooo"]("L,R,k,l,i,j");

        // flops: o3v1L2  = o4v1L2
        //  mems: o3v1L2  = o3v1L2
        tmps_["42_aaaa_LLvooo"]("L,R,a,i,j,k")  = t1_1["aa"]("a,l") * tmps_["41_aaaa_LLoooo"]("L,R,l,i,j,k");
        tmps_["41_aaaa_LLoooo"].~TArrayD();

        // D_ooov_aaaa += +0.50 r2_aaaa(b,c,i,j) t1_1_aa(a,l) l2_1_aaaa(l,k,b,c)
        D_ooov_aaaa("L,R,i,j,k,a") += 0.50 * tmps_["42_aaaa_LLvooo"]("L,R,a,k,i,j");

        // D_oovo_aaaa += -0.50 r2_aaaa(b,c,i,j) t1_1_aa(a,l) l2_1_aaaa(l,k,b,c)
        D_oovo_aaaa("L,R,i,j,a,k") -= 0.50 * tmps_["42_aaaa_LLvooo"]("L,R,a,k,i,j");

        // D_oovv_aaaa += -0.50 P(a,b) r2_aaaa(c,d,i,j) t1_aa(b,l) t1_1_aa(a,k) l2_1_aaaa(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["42_aaaa_LLvooo"]("L,R,a,l,i,j") * t1["aa"]("b,l");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();
        tmps_["42_aaaa_LLvooo"].~TArrayD();

        // flops: o4v0L2  = o4v2L2 o4v2L2 o4v0L2
        //  mems: o4v0L2  = o4v0L2 o4v0L2 o4v0L2
        tmps_["43_aaaa_LLoooo"]("L,R,i,j,k,l")  = l2["aaaa"]("L,i,j,a,b") * r2["aaaa"]("R,a,b,k,l");
        tmps_["43_aaaa_LLoooo"]("L,R,i,j,k,l") += l2_1["aaaa"]("L,i,j,a,b") * r2_1["aaaa"]("R,a,b,k,l");

        // D_oooo_aaaa += +0.50 r2_aaaa(a,b,i,j) l2_aaaa(l,k,a,b)
        //             += +0.50 r2_1_aaaa(a,b,i,j) l2_1_aaaa(l,k,a,b)
        D_oooo_aaaa("L,R,i,j,k,l") += 0.50 * tmps_["43_aaaa_LLoooo"]("L,R,l,k,i,j");

        // D_oovv_aaaa += +0.25 r2_aaaa(c,d,i,j) t2_aaaa(b,a,k,l) l2_aaaa(k,l,c,d)
        //             += +0.25 r2_1_aaaa(c,d,i,j) t2_aaaa(b,a,k,l) l2_1_aaaa(k,l,c,d)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") += 0.25 * t2["aaaa"]("b,a,k,l") * tmps_["43_aaaa_LLoooo"]("L,R,k,l,i,j");

        // flops: o3v1L2  = o4v1L2
        //  mems: o3v1L2  = o3v1L2
        tmps_["44_aaaa_LLvooo"]("L,R,a,i,j,k")  = t1["aa"]("a,l") * tmps_["43_aaaa_LLoooo"]("L,R,l,i,j,k");
        tmps_["43_aaaa_LLoooo"].~TArrayD();

        // D_ooov_aaaa += +0.50 r2_aaaa(b,c,i,j) t1_aa(a,l) l2_aaaa(l,k,b,c)
        //             += +0.50 r2_1_aaaa(b,c,i,j) t1_aa(a,l) l2_1_aaaa(l,k,b,c)
        D_ooov_aaaa("L,R,i,j,k,a") += 0.50 * tmps_["44_aaaa_LLvooo"]("L,R,a,k,i,j");

        // D_oovo_aaaa += -0.50 r2_aaaa(b,c,i,j) t1_aa(a,l) l2_aaaa(l,k,b,c)
        //             += -0.50 r2_1_aaaa(b,c,i,j) t1_aa(a,l) l2_1_aaaa(l,k,b,c)
        D_oovo_aaaa("L,R,i,j,a,k") -= 0.50 * tmps_["44_aaaa_LLvooo"]("L,R,a,k,i,j");

        // D_oovv_aaaa += -0.50 r2_aaaa(c,d,i,j) t1_aa(a,k) t1_aa(b,l) l2_aaaa(k,l,c,d)
        //             += -0.50 r2_1_aaaa(c,d,i,j) t1_aa(a,k) t1_aa(b,l) l2_1_aaaa(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") -= 0.50 * tmps_["44_aaaa_LLvooo"]("L,R,a,l,i,j") * t1["aa"]("b,l");
        tmps_["44_aaaa_LLvooo"].~TArrayD();

        // flops: o4v0L2  = o4v2L2 o4v2L2 o4v0L2
        //  mems: o4v0L2  = o4v0L2 o4v0L2 o4v0L2
        tmps_["45_bbbb_LLoooo"]("L,R,i,j,k,l")  = l2["bbbb"]("L,i,j,a,b") * r2["bbbb"]("R,a,b,k,l");
        tmps_["45_bbbb_LLoooo"]("L,R,i,j,k,l") += l2_1["bbbb"]("L,i,j,a,b") * r2_1["bbbb"]("R,a,b,k,l");

        // D_oooo_bbbb += +0.50 r2_1_bbbb(a,b,i,j) l2_1_bbbb(l,k,a,b)
        //             += +0.50 r2_bbbb(a,b,i,j) l2_bbbb(l,k,a,b)
        D_oooo_bbbb("L,R,i,j,k,l") += 0.50 * tmps_["45_bbbb_LLoooo"]("L,R,l,k,i,j");

        // D_oovv_bbbb += +0.25 r2_1_bbbb(c,d,i,j) t2_bbbb(b,a,k,l) l2_1_bbbb(k,l,c,d)
        //             += +0.25 r2_bbbb(c,d,i,j) t2_bbbb(b,a,k,l) l2_bbbb(k,l,c,d)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += 0.25 * t2["bbbb"]("b,a,k,l") * tmps_["45_bbbb_LLoooo"]("L,R,k,l,i,j");

        // flops: o3v1L2  = o4v1L2
        //  mems: o3v1L2  = o3v1L2
        tmps_["46_bbbb_LLvooo"]("L,R,a,i,j,k")  = t1["bb"]("a,l") * tmps_["45_bbbb_LLoooo"]("L,R,l,i,j,k");
        tmps_["45_bbbb_LLoooo"].~TArrayD();

        // D_ooov_bbbb += +0.50 r2_1_bbbb(b,c,i,j) t1_bb(a,l) l2_1_bbbb(l,k,b,c)
        //             += +0.50 r2_bbbb(b,c,i,j) t1_bb(a,l) l2_bbbb(l,k,b,c)
        D_ooov_bbbb("L,R,i,j,k,a") += 0.50 * tmps_["46_bbbb_LLvooo"]("L,R,a,k,i,j");

        // D_oovo_bbbb += -0.50 r2_1_bbbb(b,c,i,j) t1_bb(a,l) l2_1_bbbb(l,k,b,c)
        //             += -0.50 r2_bbbb(b,c,i,j) t1_bb(a,l) l2_bbbb(l,k,b,c)
        D_oovo_bbbb("L,R,i,j,a,k") -= 0.50 * tmps_["46_bbbb_LLvooo"]("L,R,a,k,i,j");

        // D_oovv_bbbb += -0.50 r2_1_bbbb(c,d,i,j) t1_bb(a,k) t1_bb(b,l) l2_1_bbbb(k,l,c,d)
        //             += -0.50 r2_bbbb(c,d,i,j) t1_bb(a,k) t1_bb(b,l) l2_bbbb(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") -= 0.50 * tmps_["46_bbbb_LLvooo"]("L,R,a,l,i,j") * t1["bb"]("b,l");
        tmps_["46_bbbb_LLvooo"].~TArrayD();

        // flops: o4v0L2  = o4v2L2
        //  mems: o4v0L2  = o4v0L2
        tmps_["47_bbbb_LLoooo"]("L,R,i,j,k,l")  = l2_1["bbbb"]("L,i,j,a,b") * r2["bbbb"]("R,a,b,k,l");

        // D_oovv_bbbb += +0.25 r2_bbbb(c,d,i,j) t2_1_bbbb(b,a,k,l) l2_1_bbbb(k,l,c,d)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += 0.25 * t2_1["bbbb"]("b,a,k,l") * tmps_["47_bbbb_LLoooo"]("L,R,k,l,i,j");

        // flops: o3v1L2  = o4v1L2
        //  mems: o3v1L2  = o3v1L2
        tmps_["48_bbbb_LLvooo"]("L,R,a,i,j,k")  = t1_1["bb"]("a,l") * tmps_["47_bbbb_LLoooo"]("L,R,l,i,j,k");
        tmps_["47_bbbb_LLoooo"].~TArrayD();

        // D_ooov_bbbb += +0.50 r2_bbbb(b,c,i,j) t1_1_bb(a,l) l2_1_bbbb(l,k,b,c)
        D_ooov_bbbb("L,R,i,j,k,a") += 0.50 * tmps_["48_bbbb_LLvooo"]("L,R,a,k,i,j");

        // D_oovo_bbbb += -0.50 r2_bbbb(b,c,i,j) t1_1_bb(a,l) l2_1_bbbb(l,k,b,c)
        D_oovo_bbbb("L,R,i,j,a,k") -= 0.50 * tmps_["48_bbbb_LLvooo"]("L,R,a,k,i,j");

        // D_oovv_bbbb += -0.50 P(a,b) r2_bbbb(c,d,i,j) t1_bb(b,l) t1_1_bb(a,k) l2_1_bbbb(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["48_bbbb_LLvooo"]("L,R,a,l,i,j") * t1["bb"]("b,l");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
        tmps_["48_bbbb_LLvooo"].~TArrayD();

        // flops: o3v0L1  = o3v2L1 o3v2L1 o3v0L1
        //  mems: o3v0L1  = o2v0L1 o2v0L1 o3v0L1
        tmps_["49_aaa_Looo"]("L,i,j,k")  = t2_1["abab"]("a,b,i,l") * l2_1["abab"]("L,j,l,a,b");
        tmps_["49_aaa_Looo"]("L,i,j,k") -= t2["abab"]("c,d,i,m") * l2["abab"]("L,k,m,c,d");

        // D_oovv_aaaa += +0.50 P(i,j) P(a,b) r1_aa(a,i) t1_aa(b,l) t2_1_abab(d,c,j,k) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,i) t1_aa(b,l) t2_1_abab(c,d,j,k) l2_1_abab(l,k,c,d)
        //             += -0.50 d_bb(i,k) r1_aa(a,m) t2_abab(c,b,j,l) l2_abab(m,l,c,b)
        //             += -0.50 d_bb(i,k) r1_aa(a,m) t2_abab(b,c,j,l) l2_abab(m,l,b,c)
        // flops: o2v2L2 += o3v1L1 o3v2L2
        //  mems: o2v2L2 += o2v1L1 o3v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("b,l") * tmps_["49_aaa_Looo"]("L,j,l,m") * r1["aa"]("R,a,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o2v1L2  = o3v1L2
        //  mems: o2v1L2  = o2v1L2
        tmps_["50_aaa_LLoov"]("L,R,i,j,a")  = tmps_["49_aaa_Looo"]("L,i,k,j") * r1["aa"]("R,a,k");

        // D_ooov_aaaa += -0.50 P(i,j) d_aa(j,k) r1_aa(a,m) t2_1_abab(c,b,i,l) l2_1_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r1_aa(a,m) t2_1_abab(b,c,i,l) l2_1_abab(m,l,b,c)
        //             += +0.50 P(i,j) d_aa(j,l) r0 t2_abab(b,a,i,m) l2_abab(k,m,b,a)
        //             += +0.50 P(i,j) d_aa(j,l) r0 t2_abab(a,b,i,m) l2_abab(k,m,a,b)
        // flops: o3v1L2 += o3v1L2
        //  mems: o3v1L2 += o2v1L2
        tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k")  = tmps_["50_aaa_LLoov"]("L,R,i,k,a") * Id["aa_oo"]("j,k");
        D_ooov_aaaa("L,R,i,j,k,a") -= tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k");
        D_ooov_aaaa("L,R,i,j,k,a") += tmps_["perm_aaaa_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_aaaa_LLvooo"].~TArrayD();

        // D_oovo_aaaa += +0.50 P(i,j) d_aa(j,k) r1_aa(a,m) t2_1_abab(c,b,i,l) l2_1_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r1_aa(a,m) t2_1_abab(b,c,i,l) l2_1_abab(m,l,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_abab(b,a,i,m) l2_abab(l,m,b,a)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_abab(a,b,i,m) l2_abab(l,m,a,b)
        // flops: o3v1L2 += o4v1L2
        //  mems: o3v1L2 += o4v1L2
        tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k")  = tmps_["50_aaa_LLoov"]("L,R,i,l,a") * Id["aa_oo"]("j,k");
        D_oovo_aaaa("L,R,i,j,a,k") += tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k");
        D_oovo_aaaa("L,R,i,j,a,k") -= tmps_["perm_aaaa_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_aaaa_LLvooo"].~TArrayD();

        // D_oovo_abab += -0.50 d_bb(i,k) r1_aa(a,m) t2_1_abab(c,b,j,l) l2_1_abab(m,l,c,b)
        //             += -0.50 d_bb(i,k) r1_aa(a,m) t2_1_abab(b,c,j,l) l2_1_abab(m,l,b,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,i) t1_aa(b,l) t2_abab(d,c,j,k) l2_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,i) t1_aa(b,l) t2_abab(c,d,j,k) l2_abab(l,k,c,d)
        // flops: o3v1L2 += o4v1L2
        //  mems: o3v1L2 += o4v1L2
        D_oovo_abab("L,R,i,j,a,k") -= tmps_["50_aaa_LLoov"]("L,R,j,l,a") * Id["bb_oo"]("i,k");

        // D_oovv_abab += +0.50 r1_aa(a,l) t1_bb(b,j) t2_1_abab(d,c,i,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r1_aa(a,l) t1_bb(b,j) t2_1_abab(c,d,i,k) l2_1_abab(l,k,c,d)
        //             += -0.50 P(i,j) r0 t1_aa(a,j) t2_abab(c,b,i,l) l2_abab(k,l,c,b)
        //             += -0.50 P(i,j) r0 t1_aa(a,j) t2_abab(b,c,i,l) l2_abab(k,l,b,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o3v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1["bb"]("b,j") * tmps_["50_aaa_LLoov"]("L,R,i,k,a");
        tmps_["50_aaa_LLoov"].~TArrayD();

        // flops: o2v1L2  = o3v2L1 o3v2L1 o3v0L1 o3v1L2
        //  mems: o2v1L2  = o2v0L1 o2v0L1 o3v0L1 o2v1L2
        tmps_["51_bbb_LLoov"]("L,R,i,j,a")  = (t2["abab"]("b,c,l,i") * l2["abab"]("L,l,k,b,c") + -1.00 * t2_1["abab"]("d,a,m,i") * l2_1["abab"]("L,m,j,d,a")) * r1["bb"]("R,a,k");

        // D_ooov_bbbb += -0.50 P(i,j) d_bb(j,k) r1_bb(a,m) t2_abab(c,b,l,i) l2_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r1_bb(a,m) t2_abab(b,c,l,i) l2_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_bb(j,l) r0 t2_1_abab(b,a,m,i) l2_1_abab(m,k,b,a)
        //             += +0.50 P(i,j) d_bb(j,l) r0 t2_1_abab(a,b,m,i) l2_1_abab(m,k,a,b)
        // flops: o3v1L2 += o3v1L2
        //  mems: o3v1L2 += o2v1L2
        tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k")  = tmps_["51_bbb_LLoov"]("L,R,i,k,a") * Id["bb_oo"]("j,k");
        D_ooov_bbbb("L,R,i,j,k,a") -= tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k");
        D_ooov_bbbb("L,R,i,j,k,a") += tmps_["perm_bbbb_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_bbbb_LLvooo"].~TArrayD();

        // D_oovo_bbbb += +0.50 P(i,j) d_bb(j,k) r1_bb(a,m) t2_abab(c,b,l,i) l2_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r1_bb(a,m) t2_abab(b,c,l,i) l2_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_1_abab(b,a,m,i) l2_1_abab(m,l,b,a)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_1_abab(a,b,m,i) l2_1_abab(m,l,a,b)
        // flops: o3v1L2 += o4v1L2
        //  mems: o3v1L2 += o4v1L2
        tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k")  = tmps_["51_bbb_LLoov"]("L,R,i,l,a") * Id["bb_oo"]("j,k");
        D_oovo_bbbb("L,R,i,j,a,k") += tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k");
        D_oovo_bbbb("L,R,i,j,a,k") -= tmps_["perm_bbbb_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_bbbb_LLvooo"].~TArrayD();

        // D_oovv_bbbb += +0.50 P(i,j) P(a,b) r1_bb(a,l) t1_bb(b,j) t2_abab(d,c,k,i) l2_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,l) t1_bb(b,j) t2_abab(c,d,k,i) l2_abab(k,l,c,d)
        //             += -0.50 r1_aa(a,j) t2_1_abab(c,b,l,i) l2_1_abab(l,k,c,b)
        //             += -0.50 r1_aa(a,j) t2_1_abab(b,c,l,i) l2_1_abab(l,k,b,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o3v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("b,j") * tmps_["51_bbb_LLoov"]("L,R,i,k,a");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
        tmps_["51_bbb_LLoov"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["52_aaaa_Lvoov"]("L,a,i,j,b")  = t2["aaaa"]("a,c,k,i") * l2_1["aaaa"]("L,j,k,b,c");

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["53_aaaa_Lvooo"]("L,a,i,j,k")  = tmps_["52_aaaa_Lvoov"]("L,a,i,j,b") * t1["aa"]("b,k");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r1_1_aa(a,l) t2_aaaa(b,c,k,i) t1_aa(d,j) l2_1_aaaa(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["53_aaaa_Lvooo"]("L,b,i,l,j") * r1_1["aa"]("R,a,l");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r0 t2_aaaa(b,c,l,i) t1_aa(d,j) t1_1_aa(a,k) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1_1["aa"]("a,k") * tmps_["53_aaaa_Lvooo"]("L,b,i,k,j") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();
        tmps_["53_aaaa_Lvooo"].~TArrayD();

        // flops: o2v1L1  = o3v2L1 o3v2L1 o3v0L1 o3v1L1
        //  mems: o2v1L1  = o2v0L1 o2v0L1 o3v0L1 o2v1L1
        tmps_["54_aaa_Loov"]("L,i,j,a")  = (t2_1["abab"]("b,c,i,l") * l2_1["abab"]("L,k,l,b,c") + t2["abab"]("d,e,i,m") * l2["abab"]("L,j,m,d,e")) * t1["aa"]("a,k");

        // D_oovo_aaaa += +0.50 P(i,j) d_aa(j,k) r0 t1_aa(a,m) t2_1_abab(c,b,i,l) l2_1_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t1_aa(a,m) t2_1_abab(b,c,i,l) l2_1_abab(m,l,b,c)
        //             += +0.50 d_bb(j,l) r0 t2_abab(b,a,i,m) l2_abab(k,m,b,a)
        //             += +0.50 d_bb(j,l) r0 t2_abab(a,b,i,m) l2_abab(k,m,a,b)
        // flops: o3v1L2 += o3v1L1 o2v1L2
        //  mems: o3v1L2 += o2v1L1 o2v1L2
        tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k")  = tmps_["54_aaa_Loov"]("L,i,k,a") * Id["aa_oo"]("j,k") * r0("R");
        D_oovo_aaaa("L,R,i,j,a,k") += tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k");
        D_oovo_aaaa("L,R,i,j,a,k") -= tmps_["perm_aaaa_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_aaaa_LLvooo"].~TArrayD();

        // D_oovv_abab += +0.50 r1_bb(b,j) t1_aa(a,l) t2_1_abab(d,c,i,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r1_bb(b,j) t1_aa(a,l) t2_1_abab(c,d,i,k) l2_1_abab(l,k,c,d)
        //             += +0.50 P(i,j) r0 t1_aa(a,j) t2_abab(c,b,i,l) l2_abab(k,l,c,b)
        //             += +0.50 P(i,j) r0 t1_aa(a,j) t2_abab(b,c,i,l) l2_abab(k,l,b,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o3v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["54_aaa_Loov"]("L,i,k,a") * r1["bb"]("R,b,j");
        tmps_["54_aaa_Loov"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["55_bbaa_Lvoov"]("L,a,i,j,b")  = t2["abab"]("c,a,k,i") * l2_1["aaaa"]("L,j,k,b,c");

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["56_bbaa_Lvooo"]("L,a,i,j,k")  = tmps_["55_bbaa_Lvoov"]("L,a,i,j,b") * t1["aa"]("b,k");

        // D_oovv_abab += +1.00 r1_1_aa(a,l) t2_abab(c,b,k,j) t1_aa(d,i) l2_1_aaaa(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["56_bbaa_Lvooo"]("L,b,j,l,i") * r1_1["aa"]("R,a,l");

        // D_oovv_abab += +1.00 r0 t2_abab(c,b,l,j) t1_aa(d,i) t1_1_aa(a,k) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1_1["aa"]("a,k") * tmps_["56_bbaa_Lvooo"]("L,b,j,k,i") * r0("R");
        tmps_["56_bbaa_Lvooo"].~TArrayD();

        // flops: o4v1L2  = o3v2L1 o3v2L1 o3v0L1 o4v1L2
        //  mems: o4v1L2  = o2v0L1 o2v0L1 o3v0L1 o4v1L2
        tmps_["57_bbbaa_LLooovo"]("L,R,i,j,k,a,l")  = (t2["bbbb"]("b,c,k,i") * l2["bbbb"]("L,k,j,b,c") + t2_1["bbbb"]("d,b,j,i") * l2_1["bbbb"]("L,j,k,d,b")) * r1["aa"]("R,a,l");

        // D_ooov_abab += +0.50 r1_aa(a,j) t2_bbbb(c,b,l,i) l2_bbbb(l,k,c,b)
        //             += +0.50 P(i,j) P(a,b) r0 t1_bb(a,l) t1_bb(b,j) t2_1_bbbb(d,c,k,i) l2_1_bbbb(k,l,d,c)
        D_ooov_abab("L,R,i,j,k,a") += 0.50 * tmps_["57_bbbaa_LLooovo"]("L,R,i,k,l,a,j");

        // D_oovo_abab += -0.50 r1_aa(a,j) t2_bbbb(c,b,l,i) l2_bbbb(l,k,c,b)
        //             += -0.50 P(i,j) r0 t2_bbbb(b,a,l,j) t2_1_bbbb(d,c,k,i) l2_1_bbbb(k,l,d,c)
        D_oovo_abab("L,R,i,j,a,k") -= 0.50 * tmps_["57_bbbaa_LLooovo"]("L,R,i,k,l,a,j");
        tmps_["57_bbbaa_LLooovo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["58_bbbb_Lvoov"]("L,a,i,j,b")  = t2["bbbb"]("a,c,k,i") * l2_1["bbbb"]("L,j,k,b,c");

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["59_bbbb_Lvooo"]("L,a,i,j,k")  = tmps_["58_bbbb_Lvoov"]("L,a,i,j,b") * t1["bb"]("b,k");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r1_1_bb(a,l) t2_bbbb(b,c,k,i) t1_bb(d,j) l2_1_bbbb(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["59_bbbb_Lvooo"]("L,b,i,l,j") * r1_1["bb"]("R,a,l");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r0 t2_bbbb(b,c,l,i) t1_bb(d,j) t1_1_bb(a,k) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1_1["bb"]("a,k") * tmps_["59_bbbb_Lvooo"]("L,b,i,k,j") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
        tmps_["59_bbbb_Lvooo"].~TArrayD();

        // flops: o4v0L1  = o3v2L1 o3v2L1 o4v0L1
        //  mems: o4v0L1  = o2v0L1 o2v0L1 o4v0L1
        tmps_["60_bbbb_Loooo"]("L,i,j,k,l")  = t2["abab"]("a,b,m,i") * l2["abab"]("L,m,j,a,b");
        tmps_["60_bbbb_Loooo"]("L,i,j,k,l") -= t2_1["abab"]("c,d,n,k") * l2_1["abab"]("L,n,l,c,d");

        // flops: o3v1L1  = o4v1L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["61_bbbb_Lvooo"]("L,a,i,j,k")  = t1["bb"]("a,l") * tmps_["60_bbbb_Loooo"]("L,i,l,j,k");

        // D_oovv_abab += +0.50 r1_aa(a,i) t1_bb(b,l) t2_abab(d,c,k,j) l2_abab(k,l,d,c)
        //             += +0.50 r1_aa(a,i) t1_bb(b,l) t2_abab(c,d,k,j) l2_abab(k,l,c,d)
        //             += -0.50 P(i,j) r0 t1_bb(a,j) t2_1_abab(c,b,l,i) l2_1_abab(l,k,c,b)
        //             += -0.50 P(i,j) r0 t1_bb(a,j) t2_1_abab(b,c,l,i) l2_1_abab(l,k,b,c)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o4v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["61_bbbb_Lvooo"]("L,b,j,i,k") * r1["aa"]("R,a,i");

        // D_oovv_bbbb += +0.50 P(i,j) P(a,b) r1_bb(a,i) t1_bb(b,l) t2_abab(d,c,k,j) l2_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,i) t1_bb(b,l) t2_abab(c,d,k,j) l2_abab(k,l,c,d)
        //             += -0.50 r0 t1_aa(a,j) t2_1_abab(c,b,l,i) l2_1_abab(l,k,c,b)
        //             += -0.50 r0 t1_aa(a,j) t2_1_abab(b,c,l,i) l2_1_abab(l,k,b,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["61_bbbb_Lvooo"]("L,b,j,i,k") * r1["bb"]("R,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
        tmps_["61_bbbb_Lvooo"].~TArrayD();

        // flops: o2v0L1  = o3v2L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["62_aa_Loo"]("L,i,j")  = t2["aaaa"]("a,b,k,i") * l2_1["aaaa"]("L,j,k,a,b");

        // D_oovv_aaaa += -0.50 P(i,j) r2_1_aaaa(b,a,l,i) t2_aaaa(d,c,k,j) l2_1_aaaa(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["62_aa_Loo"]("L,j,l") * r2_1["aaaa"]("R,b,a,l,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += -0.50 P(i,j) r0 t2_aaaa(d,c,l,j) t2_1_aaaa(b,a,k,i) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t2_1["aaaa"]("b,a,k,i") * tmps_["62_aa_Loo"]("L,j,k") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -0.50 r2_1_abab(a,b,l,j) t2_aaaa(d,c,k,i) l2_1_aaaa(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * tmps_["62_aa_Loo"]("L,i,l") * r2_1["abab"]("R,a,b,l,j");

        // D_oovv_abab += -0.50 r0 t2_aaaa(d,c,l,i) t2_1_abab(a,b,k,j) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * t2_1["abab"]("a,b,k,j") * tmps_["62_aa_Loo"]("L,i,k") * r0("R");

        // flops: o2v0L1  = o3v2L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["63_aa_Loo"]("L,i,j")  = t2["abab"]("a,b,i,k") * l2_1["abab"]("L,j,k,a,b");

        // D_oovv_aaaa += +0.50 P(i,j) r2_1_aaaa(b,a,l,i) t2_abab(d,c,j,k) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) r2_1_aaaa(b,a,l,i) t2_abab(c,d,j,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["63_aa_Loo"]("L,j,l") * r2_1["aaaa"]("R,b,a,l,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += -0.50 P(i,j) r0_1 t2_aaaa(b,a,l,j) t2_abab(d,c,i,k) l2_1_abab(l,k,d,c)
        //             += -0.50 P(i,j) r0_1 t2_aaaa(b,a,l,j) t2_abab(c,d,i,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["aaaa"]("b,a,l,j") * tmps_["63_aa_Loo"]("L,i,l") * r0_1("R");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +0.50 P(i,j) r0 t2_abab(d,c,j,l) t2_1_aaaa(b,a,k,i) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) r0 t2_abab(c,d,j,l) t2_1_aaaa(b,a,k,i) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2_1["aaaa"]("b,a,k,i") * tmps_["63_aa_Loo"]("L,j,k") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r2_1_abab(a,b,l,j) t2_abab(d,c,i,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r2_1_abab(a,b,l,j) t2_abab(c,d,i,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["63_aa_Loo"]("L,i,l") * r2_1["abab"]("R,a,b,l,j");

        // D_oovv_abab += +0.50 r0 t2_abab(d,c,i,l) t2_1_abab(a,b,k,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r0 t2_abab(c,d,i,l) t2_1_abab(a,b,k,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2_1["abab"]("a,b,k,j") * tmps_["63_aa_Loo"]("L,i,k") * r0("R");

        // D_oovv_abab += +0.50 r0_1 t2_abab(a,b,l,j) t2_abab(d,c,i,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r0_1 t2_abab(a,b,l,j) t2_abab(c,d,i,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,b,l,j") * tmps_["63_aa_Loo"]("L,i,l") * r0_1("R");

        // flops: o0v2L2  = o2v3L2 o2v3L2 o0v2L2
        //  mems: o0v2L2  = o0v2L2 o0v2L2 o0v2L2
        tmps_["64_aa_LLvv"]("L,R,a,b")  = 0.50 * l2_1["aaaa"]("L,i,k,a,d") * r2["aaaa"]("R,b,d,i,k");
        tmps_["64_aa_LLvv"]("L,R,a,b") += l2_1["abab"]("L,i,j,a,c") * r2["abab"]("R,b,c,i,j");

        // D_oovv_aaaa += -0.50 P(a,b) r2_abab(a,d,l,k) t2_1_aaaa(b,c,i,j) l2_1_abab(l,k,c,d)
        //             += -0.50 P(a,b) r2_abab(a,d,k,l) t2_1_aaaa(b,c,i,j) l2_1_abab(k,l,c,d)
        //             += -0.50 P(a,b) r2_aaaa(a,d,l,k) t2_1_aaaa(b,c,i,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["64_aa_LLvv"]("L,R,c,a") * t2_1["aaaa"]("b,c,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r2_abab(a,d,l,k) t2_1_abab(c,b,i,j) l2_1_abab(l,k,c,d)
        //             += +0.50 r2_abab(a,d,k,l) t2_1_abab(c,b,i,j) l2_1_abab(k,l,c,d)
        //             += +0.50 r2_aaaa(a,d,l,k) t2_1_abab(c,b,i,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["64_aa_LLvv"]("L,R,c,a") * t2_1["abab"]("c,b,i,j");

        // flops: o1v1L2  = o1v2L2 o2v1L1 o2v1L2 o2v1L2 o1v1L2 o1v1L2 o2v1L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o2v0L1 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2
        tmps_["65_aa_LLov"]("L,R,i,a")  = tmps_["64_aa_LLvv"]("L,R,b,a") * t1["aa"]("b,i");
        tmps_["65_aa_LLov"]("L,R,i,a") += l1_1["aa"]("L,k,b") * t1["aa"]("b,i") * r1["aa"]("R,a,k");
        tmps_["65_aa_LLov"]("L,R,i,a") += r1["aa"]("R,a,j") * tmps_["63_aa_Loo"]("L,i,j");
        tmps_["65_aa_LLov"]("L,R,i,a") -= 0.50 * tmps_["62_aa_Loo"]("L,i,j") * r1["aa"]("R,a,j");

        // D_oovv_aaaa += -1.00 P(i,j) P(a,b) r1_aa(a,k) t1_aa(c,j) t1_1_aa(b,i) l1_1_aa(k,c)
        //             += -0.50 P(i,j) P(a,b) r1_aa(a,l) t2_abab(d,c,j,k) t1_1_aa(b,i) l2_1_abab(l,k,d,c)
        //             += -0.50 P(i,j) P(a,b) r1_aa(a,l) t2_abab(c,d,j,k) t1_1_aa(b,i) l2_1_abab(l,k,c,d)
        //             += -0.50 P(i,j) P(a,b) r2_abab(a,d,l,k) t1_aa(c,j) t1_1_aa(b,i) l2_1_abab(l,k,c,d)
        //             += -0.50 P(i,j) P(a,b) r2_abab(a,d,k,l) t1_aa(c,j) t1_1_aa(b,i) l2_1_abab(k,l,c,d)
        //             += -0.50 P(i,j) P(a,b) r2_aaaa(a,d,l,k) t1_aa(c,j) t1_1_aa(b,i) l2_1_aaaa(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,l) t2_aaaa(d,c,k,j) t1_1_aa(b,i) l2_1_aaaa(l,k,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1_1["aa"]("b,i") * tmps_["65_aa_LLov"]("L,R,j,a");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +1.00 r1_aa(a,k) t1_aa(c,i) t1_1_bb(b,j) l1_1_aa(k,c)
        //             += +0.50 r1_aa(a,l) t2_abab(d,c,i,k) t1_1_bb(b,j) l2_1_abab(l,k,d,c)
        //             += +0.50 r1_aa(a,l) t2_abab(c,d,i,k) t1_1_bb(b,j) l2_1_abab(l,k,c,d)
        //             += +0.50 r2_abab(a,d,l,k) t1_aa(c,i) t1_1_bb(b,j) l2_1_abab(l,k,c,d)
        //             += +0.50 r2_abab(a,d,k,l) t1_aa(c,i) t1_1_bb(b,j) l2_1_abab(k,l,c,d)
        //             += +0.50 r2_aaaa(a,d,l,k) t1_aa(c,i) t1_1_bb(b,j) l2_1_aaaa(l,k,c,d)
        //             += -0.50 r1_aa(a,l) t2_aaaa(d,c,k,i) t1_1_bb(b,j) l2_1_aaaa(l,k,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1_1["bb"]("b,j") * tmps_["65_aa_LLov"]("L,R,i,a");
        tmps_["65_aa_LLov"].~TArrayD();

        // flops: o2v2L2  = o3v3L2
        //  mems: o2v2L2  = o2v2L2
        tmps_["66_bbaa_LLovvo"]("L,R,i,a,b,j")  = l2_1["abab"]("L,k,i,c,a") * r2["aaaa"]("R,b,c,k,j");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_aaaa(a,d,l,i) t2_1_abab(b,c,j,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["66_bbaa_LLovvo"]("L,R,k,c,a,i") * t2_1["abab"]("b,c,j,k");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r2_aaaa(a,d,l,i) t2_1_bbbb(b,c,k,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["66_bbaa_LLovvo"]("L,R,k,c,a,i") * t2_1["bbbb"]("b,c,k,j");

        // D_oovv_abab += -1.00 r2_aaaa(a,d,l,i) t1_bb(c,j) t1_1_bb(b,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["66_bbaa_LLovvo"]("L,R,k,c,a,i") * t1["bb"]("c,j") * t1_1["bb"]("b,k");

        // flops: o2v2L2  = o3v3L2 o3v3L2 o2v2L2
        //  mems: o2v2L2  = o2v2L2 o2v2L2 o2v2L2
        tmps_["67_bbaa_LLovvo"]("L,R,i,a,b,j")  = l2["abab"]("L,k,i,c,a") * r2["aaaa"]("R,b,c,k,j");
        tmps_["67_bbaa_LLovvo"]("L,R,i,a,b,j") += l2_1["abab"]("L,k,i,c,a") * r2_1["aaaa"]("R,b,c,k,j");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_aaaa(a,d,l,i) t2_abab(b,c,j,k) l2_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r2_1_aaaa(a,d,l,i) t2_abab(b,c,j,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["67_bbaa_LLovvo"]("L,R,k,c,a,i") * t2["abab"]("b,c,j,k");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r2_aaaa(a,d,l,i) t2_bbbb(b,c,k,j) l2_abab(l,k,d,c)
        //             += -1.00 r2_1_aaaa(a,d,l,i) t2_bbbb(b,c,k,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["67_bbaa_LLovvo"]("L,R,k,c,a,i") * t2["bbbb"]("b,c,k,j");

        // flops: o2v2L2  = o3v3L2
        //  mems: o2v2L2  = o2v2L2
        tmps_["68_baab_LLovvo"]("L,R,i,a,b,j")  = l2_1["abab"]("L,k,i,a,c") * r2["abab"]("R,b,c,k,j");

        // D_oovv_abab += -1.00 r2_abab(a,d,l,j) t1_aa(c,i) t1_1_bb(b,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["68_baab_LLovvo"]("L,R,k,c,a,j") * t1["aa"]("c,i") * t1_1["bb"]("b,k");

        // flops: o1v1L2  = o2v2L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["69_bb_LLov"]("L,R,i,a")  = l2_1["abab"]("L,j,i,b,a") * r1["aa"]("R,b,j");

        // D_oovv_abab += +1.00 r1_aa(d,l) t2_abab(a,c,i,j) t1_1_bb(b,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,c,i,j") * tmps_["69_bb_LLov"]("L,R,k,c") * t1_1["bb"]("b,k");

        // D_oovv_bbbb += -1.00 P(a,b) r1_aa(d,l) t2_bbbb(b,c,i,j) t1_1_bb(a,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["bbbb"]("b,c,i,j") * tmps_["69_bb_LLov"]("L,R,k,c") * t1_1["bb"]("a,k");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v1L2  = o2v2L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["72_bb_LLov"]("L,R,i,a")  = l2_1["bbbb"]("L,j,i,a,b") * r1["bb"]("R,b,j");

        // D_oovv_abab += -1.00 r1_bb(d,l) t2_abab(a,c,i,j) t1_1_bb(b,k) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,c,i,j") * tmps_["72_bb_LLov"]("L,R,k,c") * t1_1["bb"]("b,k");

        // D_oovv_bbbb += +1.00 P(a,b) r1_bb(d,l) t2_bbbb(b,c,i,j) t1_1_bb(a,k) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["bbbb"]("b,c,i,j") * tmps_["72_bb_LLov"]("L,R,k,c") * t1_1["bb"]("a,k");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v2L2  = o3v3L2 o3v3L2 o2v2L2
        //  mems: o2v2L2  = o2v2L2 o2v2L2 o2v2L2
        tmps_["75_baab_LLovvo"]("L,R,i,a,b,j")  = l2["abab"]("L,k,i,a,c") * r2["abab"]("R,b,c,k,j");
        tmps_["75_baab_LLovvo"]("L,R,i,a,b,j") += l2_1["abab"]("L,k,i,a,c") * r2_1["abab"]("R,b,c,k,j");

        // flops: o1v1L2  = o2v2L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["74_bb_LLov"]("L,R,i,a")  = l2["bbbb"]("L,j,i,a,b") * r1["bb"]("R,b,j");

        // flops: o1v1L2  = o2v2L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["73_bb_LLov"]("L,R,i,a")  = l2_1["bbbb"]("L,j,i,a,b") * r1_1["bb"]("R,b,j");

        // flops: o1v1L2  = o2v2L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["71_bb_LLov"]("L,R,i,a")  = l2_1["abab"]("L,j,i,b,a") * r1_1["aa"]("R,b,j");

        // flops: o1v1L2  = o2v2L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["70_bb_LLov"]("L,R,i,a")  = l2["abab"]("L,j,i,b,a") * r1["aa"]("R,b,j");

        // flops: o3v1L2  = o3v2L2 o3v2L2 o3v2L2 o3v2L2 o3v2L2 o3v2L1 o4v2L2 o3v1L2 o3v1L2 o3v2L2 o3v1L2 o3v1L2 o3v1L2 o3v2L1 o4v2L2 o3v1L2 o3v2L1 o4v2L2 o3v1L2 o3v1L2 o3v2L2 o3v1L2 o3v2L2 o3v1L2 o3v2L2 o3v1L2 o3v2L2 o3v1L2 o3v2L2 o3v1L2 o3v2L2 o3v1L2
        //  mems: o3v1L2  = o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L1 o3v1L2 o3v1L2 o3v1L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        tmps_["76_bbaa_LLoovo"]("L,R,i,j,a,k")  = -1.00 * t2["abab"]("a,c,k,i") * tmps_["70_bb_LLov"]("L,R,j,c");
        tmps_["76_bbaa_LLoovo"]("L,R,i,j,a,k") -= t2["abab"]("a,c,k,i") * tmps_["71_bb_LLov"]("L,R,j,c");
        tmps_["76_bbaa_LLoovo"]("L,R,i,j,a,k") -= t2_1["abab"]("a,c,k,i") * tmps_["69_bb_LLov"]("L,R,j,c");
        tmps_["76_bbaa_LLoovo"]("L,R,i,j,a,k") += t2["abab"]("a,c,k,i") * tmps_["74_bb_LLov"]("L,R,j,c");
        tmps_["76_bbaa_LLoovo"]("L,R,i,j,a,k") -= l1_1["bb"]("L,j,c") * r2_1["abab"]("R,a,c,k,i");
        tmps_["76_bbaa_LLoovo"]("L,R,i,j,a,k") += l2_1["bbbb"]("L,m,j,c,d") * t1_1["bb"]("c,i") * r2["abab"]("R,a,d,k,m");
        tmps_["76_bbaa_LLoovo"]("L,R,i,j,a,k") += t2_1["abab"]("a,c,k,i") * tmps_["72_bb_LLov"]("L,R,j,c");
        tmps_["76_bbaa_LLoovo"]("L,R,i,j,a,k") += l2_1["bbbb"]("L,m,j,c,d") * t1["bb"]("c,i") * r2_1["abab"]("R,a,d,k,m");
        tmps_["76_bbaa_LLoovo"]("L,R,i,j,a,k") += l2["bbbb"]("L,m,j,c,d") * t1["bb"]("c,i") * r2["abab"]("R,a,d,k,m");
        tmps_["76_bbaa_LLoovo"]("L,R,i,j,a,k") -= l1["bb"]("L,j,c") * r2["abab"]("R,a,c,k,i");
        tmps_["76_bbaa_LLoovo"]("L,R,i,j,a,k") += t2["abab"]("a,c,k,i") * tmps_["73_bb_LLov"]("L,R,j,c");
        tmps_["76_bbaa_LLoovo"]("L,R,i,j,a,k") += t1_1["bb"]("c,i") * tmps_["66_bbaa_LLovvo"]("L,R,j,c,a,k");
        tmps_["76_bbaa_LLoovo"]("L,R,i,j,a,k") += t1_1["aa"]("b,k") * tmps_["68_baab_LLovvo"]("L,R,j,b,a,i");
        tmps_["76_bbaa_LLoovo"]("L,R,i,j,a,k") += t1["bb"]("c,i") * tmps_["67_bbaa_LLovvo"]("L,R,j,c,a,k");
        tmps_["76_bbaa_LLoovo"]("L,R,i,j,a,k") += t1["aa"]("b,k") * tmps_["75_baab_LLovvo"]("L,R,j,b,a,i");
        tmps_["67_bbaa_LLovvo"].~TArrayD();
        tmps_["66_bbaa_LLovvo"].~TArrayD();

        // D_ooov_abab += -1.00 r2_abab(a,c,j,l) t1_1_bb(b,i) l2_1_bbbb(l,k,b,c)
        //             += +1.00 r2_1_abab(a,b,j,i) l1_1_bb(k,b)
        //             += -1.00 r1_bb(c,l) t2_abab(a,b,j,i) l2_bbbb(l,k,b,c)
        //             += -1.00 r1_bb(c,l) t2_1_abab(a,b,j,i) l2_1_bbbb(l,k,b,c)
        //             += +1.00 r1_aa(c,l) t2_1_abab(a,b,j,i) l2_1_abab(l,k,c,b)
        //             += +1.00 r1_1_aa(c,l) t2_abab(a,b,j,i) l2_1_abab(l,k,c,b)
        //             += -1.00 r2_1_abab(a,c,j,l) t1_bb(b,i) l2_1_bbbb(l,k,b,c)
        //             += -1.00 r2_abab(a,c,j,l) t1_bb(b,i) l2_bbbb(l,k,b,c)
        //             += +1.00 r1_aa(c,l) t2_abab(a,b,j,i) l2_abab(l,k,c,b)
        //             += +1.00 r2_abab(a,b,j,i) l1_bb(k,b)
        //             += -1.00 r1_1_bb(c,l) t2_abab(a,b,j,i) l2_1_bbbb(l,k,b,c)
        //             += -1.00 r2_aaaa(a,c,l,j) t1_1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += -1.00 r2_abab(a,c,l,i) t1_1_aa(b,j) l2_1_abab(l,k,b,c)
        //             += -1.00 r2_aaaa(a,c,l,j) t1_bb(b,i) l2_abab(l,k,c,b)
        //             += -1.00 r2_1_aaaa(a,c,l,j) t1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += -1.00 r2_abab(a,c,l,i) t1_aa(b,j) l2_abab(l,k,b,c)
        //             += -1.00 r2_1_abab(a,c,l,i) t1_aa(b,j) l2_1_abab(l,k,b,c)
        D_ooov_abab("L,R,i,j,k,a") -= tmps_["76_bbaa_LLoovo"]("L,R,i,k,a,j");

        // D_oovo_abab += +1.00 r2_abab(a,c,j,l) t1_1_bb(b,i) l2_1_bbbb(l,k,b,c)
        //             += -1.00 r2_1_abab(a,b,j,i) l1_1_bb(k,b)
        //             += +1.00 r1_bb(c,l) t2_abab(a,b,j,i) l2_bbbb(l,k,b,c)
        //             += +1.00 r1_bb(c,l) t2_1_abab(a,b,j,i) l2_1_bbbb(l,k,b,c)
        //             += -1.00 r1_aa(c,l) t2_1_abab(a,b,j,i) l2_1_abab(l,k,c,b)
        //             += -1.00 r1_1_aa(c,l) t2_abab(a,b,j,i) l2_1_abab(l,k,c,b)
        //             += +1.00 r2_1_abab(a,c,j,l) t1_bb(b,i) l2_1_bbbb(l,k,b,c)
        //             += +1.00 r2_abab(a,c,j,l) t1_bb(b,i) l2_bbbb(l,k,b,c)
        //             += -1.00 r1_aa(c,l) t2_abab(a,b,j,i) l2_abab(l,k,c,b)
        //             += -1.00 r2_abab(a,b,j,i) l1_bb(k,b)
        //             += +1.00 r1_1_bb(c,l) t2_abab(a,b,j,i) l2_1_bbbb(l,k,b,c)
        //             += +1.00 r2_aaaa(a,c,l,j) t1_1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += +1.00 r2_abab(a,c,l,i) t1_1_aa(b,j) l2_1_abab(l,k,b,c)
        //             += +1.00 r2_aaaa(a,c,l,j) t1_bb(b,i) l2_abab(l,k,c,b)
        //             += +1.00 r2_1_aaaa(a,c,l,j) t1_bb(b,i) l2_1_abab(l,k,c,b)
        //             += +1.00 r2_abab(a,c,l,i) t1_aa(b,j) l2_abab(l,k,b,c)
        //             += +1.00 r2_1_abab(a,c,l,i) t1_aa(b,j) l2_1_abab(l,k,b,c)
        D_oovo_abab("L,R,i,j,a,k") += tmps_["76_bbaa_LLoovo"]("L,R,i,k,a,j");

        // D_oovv_abab += -1.00 r2_abab(a,d,i,l) t1_bb(b,k) t1_1_bb(c,j) l2_1_bbbb(l,k,c,d)
        //             += +1.00 r2_1_abab(a,c,i,j) t1_bb(b,k) l1_1_bb(k,c)
        //             += -1.00 r1_bb(d,l) t2_abab(a,c,i,j) t1_bb(b,k) l2_bbbb(l,k,c,d)
        //             += -1.00 r1_bb(d,l) t1_bb(b,k) t2_1_abab(a,c,i,j) l2_1_bbbb(l,k,c,d)
        //             += +1.00 r1_aa(d,l) t1_bb(b,k) t2_1_abab(a,c,i,j) l2_1_abab(l,k,d,c)
        //             += +1.00 r1_1_aa(d,l) t2_abab(a,c,i,j) t1_bb(b,k) l2_1_abab(l,k,d,c)
        //             += -1.00 r2_1_abab(a,d,i,l) t1_bb(b,k) t1_bb(c,j) l2_1_bbbb(l,k,c,d)
        //             += -1.00 r2_abab(a,d,i,l) t1_bb(b,k) t1_bb(c,j) l2_bbbb(l,k,c,d)
        //             += +1.00 r1_aa(d,l) t2_abab(a,c,i,j) t1_bb(b,k) l2_abab(l,k,d,c)
        //             += +1.00 r2_abab(a,c,i,j) t1_bb(b,k) l1_bb(k,c)
        //             += -1.00 r1_1_bb(d,l) t2_abab(a,c,i,j) t1_bb(b,k) l2_1_bbbb(l,k,c,d)
        //             += -1.00 r2_aaaa(a,d,l,i) t1_bb(b,k) t1_1_bb(c,j) l2_1_abab(l,k,d,c)
        //             += -1.00 r2_abab(a,d,l,j) t1_bb(b,k) t1_1_aa(c,i) l2_1_abab(l,k,c,d)
        //             += -1.00 r2_aaaa(a,d,l,i) t1_bb(b,k) t1_bb(c,j) l2_abab(l,k,d,c)
        //             += -1.00 r2_1_aaaa(a,d,l,i) t1_bb(b,k) t1_bb(c,j) l2_1_abab(l,k,d,c)
        //             += -1.00 r2_abab(a,d,l,j) t1_bb(b,k) t1_aa(c,i) l2_abab(l,k,c,d)
        //             += -1.00 r2_1_abab(a,d,l,j) t1_bb(b,k) t1_aa(c,i) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["bb"]("b,k") * tmps_["76_bbaa_LLoovo"]("L,R,j,k,a,i");
        tmps_["76_bbaa_LLoovo"].~TArrayD();

        // flops: o2v2L2  = o3v3L2 o3v3L2 o2v2L2
        //  mems: o2v2L2  = o2v2L2 o2v2L2 o2v2L2
        tmps_["77_bbbb_LLovvo"]("L,R,i,a,b,j")  = l2_1["abab"]("L,k,i,c,a") * r2["abab"]("R,c,b,k,j");
        tmps_["77_bbbb_LLovvo"]("L,R,i,a,b,j") += l2_1["bbbb"]("L,l,i,a,d") * r2["bbbb"]("R,b,d,l,j");

        // D_oovv_abab += -1.00 r2_abab(d,b,l,j) t2_1_abab(a,c,i,k) l2_1_abab(l,k,d,c)
        //             += -1.00 r2_bbbb(b,d,l,j) t2_1_abab(a,c,i,k) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2_1["abab"]("a,c,i,k") * tmps_["77_bbbb_LLovvo"]("L,R,k,c,b,j");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r2_abab(d,a,l,i) t2_1_bbbb(b,c,k,j) l2_1_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r2_bbbb(a,d,l,i) t2_1_bbbb(b,c,k,j) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["77_bbbb_LLovvo"]("L,R,k,c,a,i") * t2_1["bbbb"]("b,c,k,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r2_abab(d,a,l,i) t1_bb(c,j) t1_1_bb(b,k) l2_1_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r2_bbbb(a,d,l,i) t1_bb(c,j) t1_1_bb(b,k) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["77_bbbb_LLovvo"]("L,R,k,c,a,i") * t1["bb"]("c,j") * t1_1["bb"]("b,k");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v2L2  = o3v3L2 o3v3L2 o2v2L2 o3v3L2 o2v2L2 o3v3L2 o2v2L2
        //  mems: o2v2L2  = o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2
        tmps_["78_bbbb_LLovvo"]("L,R,i,a,b,j")  = l2["abab"]("L,l,i,d,a") * r2["abab"]("R,d,b,l,j");
        tmps_["78_bbbb_LLovvo"]("L,R,i,a,b,j") += l2_1["bbbb"]("L,k,i,a,c") * r2_1["bbbb"]("R,b,c,k,j");
        tmps_["78_bbbb_LLovvo"]("L,R,i,a,b,j") += l2_1["abab"]("L,l,i,d,a") * r2_1["abab"]("R,d,b,l,j");
        tmps_["78_bbbb_LLovvo"]("L,R,i,a,b,j") += l2["bbbb"]("L,k,i,a,c") * r2["bbbb"]("R,b,c,k,j");

        // D_ovvo_bbbb  = -1.00 r2_1_bbbb(b,c,k,i) l2_1_bbbb(k,j,a,c)
        //             += -1.00 r2_abab(c,b,k,i) l2_abab(k,j,c,a)
        //             += -1.00 r2_1_abab(c,b,k,i) l2_1_abab(k,j,c,a)
        //             += -1.00 r2_bbbb(b,c,k,i) l2_bbbb(k,j,a,c)
        D_ovvo_bbbb("L,R,i,a,b,j")  = -1.00 * tmps_["78_bbbb_LLovvo"]("L,R,j,a,b,i");

        // D_ovov_bbbb  = +1.00 r2_1_bbbb(b,c,k,i) l2_1_bbbb(k,j,a,c)
        //             += +1.00 r2_abab(c,b,k,i) l2_abab(k,j,c,a)
        //             += +1.00 r2_1_abab(c,b,k,i) l2_1_abab(k,j,c,a)
        //             += +1.00 r2_bbbb(b,c,k,i) l2_bbbb(k,j,a,c)
        D_ovov_bbbb("L,R,i,a,j,b")  = tmps_["78_bbbb_LLovvo"]("L,R,j,a,b,i");

        // D_oovv_abab += -1.00 r2_1_bbbb(b,d,l,j) t2_abab(a,c,i,k) l2_1_bbbb(l,k,c,d)
        //             += -1.00 r2_abab(d,b,l,j) t2_abab(a,c,i,k) l2_abab(l,k,d,c)
        //             += -1.00 r2_1_abab(d,b,l,j) t2_abab(a,c,i,k) l2_1_abab(l,k,d,c)
        //             += -1.00 r2_bbbb(b,d,l,j) t2_abab(a,c,i,k) l2_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,c,i,k") * tmps_["78_bbbb_LLovvo"]("L,R,k,c,b,j");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r2_1_bbbb(a,d,l,i) t2_bbbb(b,c,k,j) l2_1_bbbb(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r2_abab(d,a,l,i) t2_bbbb(b,c,k,j) l2_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r2_1_abab(d,a,l,i) t2_bbbb(b,c,k,j) l2_1_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r2_bbbb(a,d,l,i) t2_bbbb(b,c,k,j) l2_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["78_bbbb_LLovvo"]("L,R,k,c,a,i") * t2["bbbb"]("b,c,k,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o3v1L2  = o3v2L2 o3v2L2 o3v1L2 o3v2L1 o4v1L1 o4v1L2 o3v1L2
        //  mems: o3v1L2  = o3v1L2 o3v1L2 o3v1L2 o3v1L1 o4v0L1 o3v1L2 o3v1L2
        tmps_["79_bbbb_LLoovo"]("L,R,i,j,a,k")  = t1["bb"]("b,i") * tmps_["78_bbbb_LLovvo"]("L,R,j,b,a,k");
        tmps_["79_bbbb_LLoovo"]("L,R,i,j,a,k") += t1_1["bb"]("b,i") * tmps_["77_bbbb_LLovvo"]("L,R,j,b,a,k");
        tmps_["79_bbbb_LLoovo"]("L,R,i,j,a,k") += t1["bb"]("c,i") * l2_1["bbbb"]("L,l,j,c,b") * t1_1["bb"]("b,k") * r1["bb"]("R,a,l");

        // D_ooov_bbbb += -1.00 P(i,j) r2_abab(c,a,l,i) t1_1_bb(b,j) l2_1_abab(l,k,c,b)
        //             += -1.00 P(i,j) r2_bbbb(a,c,l,i) t1_1_bb(b,j) l2_1_bbbb(l,k,b,c)
        //             += -1.00 P(i,j) r2_1_bbbb(a,c,l,i) t1_bb(b,j) l2_1_bbbb(l,k,b,c)
        //             += -1.00 P(i,j) r2_abab(c,a,l,i) t1_bb(b,j) l2_abab(l,k,c,b)
        //             += -1.00 P(i,j) r2_1_abab(c,a,l,i) t1_bb(b,j) l2_1_abab(l,k,c,b)
        //             += -1.00 P(i,j) r2_bbbb(a,c,l,i) t1_bb(b,j) l2_bbbb(l,k,b,c)
        //             += -1.00 P(i,j) r1_bb(a,l) t1_bb(c,j) t1_1_bb(b,i) l2_1_bbbb(l,k,c,b)
        D_ooov_bbbb("L,R,i,j,k,a") -= tmps_["79_bbbb_LLoovo"]("L,R,j,k,a,i");
        D_ooov_bbbb("L,R,i,j,k,a") += tmps_["79_bbbb_LLoovo"]("L,R,i,k,a,j");

        // D_oovo_bbbb += +1.00 P(i,j) r2_abab(c,a,l,i) t1_1_bb(b,j) l2_1_abab(l,k,c,b)
        //             += +1.00 P(i,j) r2_bbbb(a,c,l,i) t1_1_bb(b,j) l2_1_bbbb(l,k,b,c)
        //             += +1.00 P(i,j) r2_1_bbbb(a,c,l,i) t1_bb(b,j) l2_1_bbbb(l,k,b,c)
        //             += +1.00 P(i,j) r2_abab(c,a,l,i) t1_bb(b,j) l2_abab(l,k,c,b)
        //             += +1.00 P(i,j) r2_1_abab(c,a,l,i) t1_bb(b,j) l2_1_abab(l,k,c,b)
        //             += +1.00 P(i,j) r2_bbbb(a,c,l,i) t1_bb(b,j) l2_bbbb(l,k,b,c)
        //             += +1.00 P(i,j) r1_bb(a,l) t1_bb(c,j) t1_1_bb(b,i) l2_1_bbbb(l,k,c,b)
        D_oovo_bbbb("L,R,i,j,a,k") += tmps_["79_bbbb_LLoovo"]("L,R,j,k,a,i");
        D_oovo_bbbb("L,R,i,j,a,k") -= tmps_["79_bbbb_LLoovo"]("L,R,i,k,a,j");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r2_abab(d,a,l,i) t1_bb(b,k) t1_1_bb(c,j) l2_1_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r2_bbbb(a,d,l,i) t1_bb(b,k) t1_1_bb(c,j) l2_1_bbbb(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r2_1_bbbb(a,d,l,i) t1_bb(b,k) t1_bb(c,j) l2_1_bbbb(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r2_abab(d,a,l,i) t1_bb(b,k) t1_bb(c,j) l2_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r2_1_abab(d,a,l,i) t1_bb(b,k) t1_bb(c,j) l2_1_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r2_bbbb(a,d,l,i) t1_bb(b,k) t1_bb(c,j) l2_bbbb(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_bb(a,l) t1_bb(b,k) t1_bb(d,j) t1_1_bb(c,i) l2_1_bbbb(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("b,k") * tmps_["79_bbbb_LLoovo"]("L,R,j,k,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
        tmps_["79_bbbb_LLoovo"].~TArrayD();

        // flops: o2v2L2  = o3v3L1 o2v2L2 o3v3L1 o2v2L2 o2v2L2
        //  mems: o2v2L2  = o2v2L1 o2v2L2 o2v2L1 o2v2L2 o2v2L2
        tmps_["80_aaaa_LLvovo"]("L,R,a,i,b,j")  = t2["abab"]("a,c,i,k") * tmps_["17_aabb_Lvoov"]("L,b,j,k,c") * r0("R");
        tmps_["80_aaaa_LLvovo"]("L,R,a,i,b,j") += t2["abab"]("b,d,j,l") * tmps_["16_aabb_Lvoov"]("L,a,i,l,d") * r0_1("R");

        // D_oovv_aaaa += +1.00 P(i,j) r0 t2_abab(a,c,i,k) t2_aaaa(b,d,l,j) l2_abab(l,k,d,c)
        //             += +1.00 P(i,j) r0_1 t2_aaaa(a,c,k,i) t2_abab(b,d,j,l) l2_1_abab(k,l,c,d)
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["80_aaaa_LLvovo"]("L,R,a,i,b,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["80_aaaa_LLvovo"]("L,R,a,j,b,i");

        // D_oovv_aaaa += +1.00 P(i,j) r0 t2_aaaa(a,c,k,i) t2_abab(b,d,j,l) l2_abab(k,l,c,d)
        //             += +1.00 P(i,j) r0_1 t2_abab(a,c,i,k) t2_aaaa(b,d,l,j) l2_1_abab(l,k,d,c)
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["80_aaaa_LLvovo"]("L,R,b,j,a,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["80_aaaa_LLvovo"]("L,R,b,i,a,j");
        tmps_["80_aaaa_LLvovo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["81_bbbb_Lvoov"]("L,a,i,j,b")  = t2["abab"]("c,a,k,i") * l2_1["abab"]("L,k,j,c,b");

        // D_oovv_abab += -1.00 r2_1_abab(a,d,i,l) t2_abab(c,b,k,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["81_bbbb_Lvoov"]("L,b,j,l,d") * r2_1["abab"]("R,a,d,i,l");

        // D_oovv_abab += -1.00 r0 t2_abab(d,b,l,j) t2_1_abab(a,c,i,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2_1["abab"]("a,c,i,k") * tmps_["81_bbbb_Lvoov"]("L,b,j,k,c") * r0("R");

        // D_oovv_abab += -1.00 r0_1 t2_abab(a,c,i,k) t2_abab(d,b,l,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,c,i,k") * tmps_["81_bbbb_Lvoov"]("L,b,j,k,c") * r0_1("R");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r2_1_bbbb(a,d,l,i) t2_abab(c,b,k,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["81_bbbb_Lvoov"]("L,b,j,l,d") * r2_1["bbbb"]("R,a,d,l,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r0 t2_abab(d,b,l,j) t2_1_bbbb(a,c,k,i) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2_1["bbbb"]("a,c,k,i") * tmps_["81_bbbb_Lvoov"]("L,b,j,k,c") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r1_bb(d,i) t2_abab(c,b,l,j) t1_1_bb(a,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["81_bbbb_Lvoov"]("L,b,j,k,d") * r1["bb"]("R,d,i") * t1_1["bb"]("a,k");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["82_bbbb_Lvoov"]("L,a,i,j,b")  = t2["abab"]("c,a,k,i") * l2["abab"]("L,k,j,c,b");

        // D_oovv_abab += -1.00 r2_abab(a,d,i,l) t2_abab(c,b,k,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["82_bbbb_Lvoov"]("L,b,j,l,d") * r2["abab"]("R,a,d,i,l");

        // D_oovv_abab += -1.00 r0 t2_abab(a,c,i,k) t2_abab(d,b,l,j) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,c,i,k") * tmps_["82_bbbb_Lvoov"]("L,b,j,k,c") * r0("R");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r2_bbbb(a,d,l,i) t2_abab(c,b,k,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["82_bbbb_Lvoov"]("L,b,j,l,d") * r2["bbbb"]("R,a,d,l,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v2L2  = o3v3L1 o2v2L2 o3v3L1 o2v2L2 o2v2L2
        //  mems: o2v2L2  = o2v2L1 o2v2L2 o2v2L1 o2v2L2 o2v2L2
        tmps_["83_bbbb_LLvovo"]("L,R,a,i,b,j")  = t2["bbbb"]("a,c,k,i") * tmps_["81_bbbb_Lvoov"]("L,b,j,k,c") * r0_1("R");
        tmps_["83_bbbb_LLvovo"]("L,R,a,i,b,j") += t2["bbbb"]("b,d,l,j") * tmps_["82_bbbb_Lvoov"]("L,a,i,l,d") * r0("R");

        // D_oovv_bbbb += +1.00 P(i,j) r0_1 t2_abab(c,a,k,i) t2_bbbb(b,d,l,j) l2_1_abab(k,l,c,d)
        //             += +1.00 P(i,j) r0 t2_bbbb(a,c,k,i) t2_abab(d,b,l,j) l2_abab(l,k,d,c)
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["83_bbbb_LLvovo"]("L,R,b,j,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["83_bbbb_LLvovo"]("L,R,b,i,a,j");

        // D_oovv_bbbb += +1.00 P(i,j) r0_1 t2_bbbb(a,c,k,i) t2_abab(d,b,l,j) l2_1_abab(l,k,d,c)
        //             += +1.00 P(i,j) r0 t2_abab(c,a,k,i) t2_bbbb(b,d,l,j) l2_abab(k,l,c,d)
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["83_bbbb_LLvovo"]("L,R,a,i,b,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["83_bbbb_LLvovo"]("L,R,a,j,b,i");
        tmps_["83_bbbb_LLvovo"].~TArrayD();

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["84_bb_Lvv"]("L,a,b")  = t2["abab"]("c,a,i,j") * l2_1["abab"]("L,i,j,c,b");

        // D_oovv_abab += +0.50 r2_1_abab(a,d,i,j) t2_abab(c,b,k,l) l2_1_abab(k,l,c,d)
        //             += +0.50 r2_1_abab(a,d,i,j) t2_abab(c,b,l,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["84_bb_Lvv"]("L,b,d") * r2_1["abab"]("R,a,d,i,j");

        // D_oovv_abab += +0.50 r0_1 t2_abab(a,c,i,j) t2_abab(d,b,k,l) l2_1_abab(k,l,d,c)
        //             += +0.50 r0_1 t2_abab(a,c,i,j) t2_abab(d,b,l,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,c,i,j") * tmps_["84_bb_Lvv"]("L,b,c") * r0_1("R");

        // D_oovv_abab += +0.50 r0 t2_abab(d,b,k,l) t2_1_abab(a,c,i,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r0 t2_abab(d,b,l,k) t2_1_abab(a,c,i,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2_1["abab"]("a,c,i,j") * tmps_["84_bb_Lvv"]("L,b,c") * r0("R");

        // D_oovv_bbbb += +0.50 P(a,b) r2_1_bbbb(a,d,i,j) t2_abab(c,b,k,l) l2_1_abab(k,l,c,d)
        //             += +0.50 P(a,b) r2_1_bbbb(a,d,i,j) t2_abab(c,b,l,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["84_bb_Lvv"]("L,b,d") * r2_1["bbbb"]("R,a,d,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +0.50 P(a,b) r0 t2_abab(d,b,k,l) t2_1_bbbb(a,c,i,j) l2_1_abab(k,l,d,c)
        //             += +0.50 P(a,b) r0 t2_abab(d,b,l,k) t2_1_bbbb(a,c,i,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2_1["bbbb"]("a,c,i,j") * tmps_["84_bb_Lvv"]("L,b,c") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["85_bb_Lvv"]("L,a,b")  = l2["abab"]("L,i,j,c,a") * t2["abab"]("c,b,i,j");

        // D_oovv_abab += +0.50 r2_abab(a,d,i,j) t2_abab(c,b,k,l) l2_abab(k,l,c,d)
        //             += +0.50 r2_abab(a,d,i,j) t2_abab(c,b,l,k) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["85_bb_Lvv"]("L,d,b") * r2["abab"]("R,a,d,i,j");

        // D_oovv_abab += +0.50 r0 t2_abab(a,c,i,j) t2_abab(d,b,k,l) l2_abab(k,l,d,c)
        //             += +0.50 r0 t2_abab(a,c,i,j) t2_abab(d,b,l,k) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,c,i,j") * tmps_["85_bb_Lvv"]("L,c,b") * r0("R");

        // D_oovv_bbbb += +0.50 P(a,b) r2_bbbb(a,d,i,j) t2_abab(c,b,k,l) l2_abab(k,l,c,d)
        //             += +0.50 P(a,b) r2_bbbb(a,d,i,j) t2_abab(c,b,l,k) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["85_bb_Lvv"]("L,d,b") * r2["bbbb"]("R,a,d,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v2L2  = o0v2L2 o0v2L2 o0v2L2 o2v3L2
        //  mems: o2v2L2  = o0v2L2 o0v2L2 o0v2L2 o2v2L2
        tmps_["86_bbbb_LLvvoo"]("L,R,a,b,i,j")  = (tmps_["85_bb_Lvv"]("L,c,a") * r0("R") + tmps_["84_bb_Lvv"]("L,a,c") * r0_1("R")) * t2["bbbb"]("b,c,i,j");

        // D_oovv_bbbb += -0.50 r0 t2_abab(c,a,k,l) t2_bbbb(b,d,i,j) l2_abab(k,l,c,d)
        //             += -0.50 r0 t2_abab(c,a,l,k) t2_bbbb(b,d,i,j) l2_abab(l,k,c,d)
        //             += -0.50 r0_1 t2_abab(c,a,k,l) t2_bbbb(b,d,i,j) l2_1_abab(k,l,c,d)
        //             += -0.50 r0_1 t2_abab(c,a,l,k) t2_bbbb(b,d,i,j) l2_1_abab(l,k,c,d)
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["86_bbbb_LLvvoo"]("L,R,a,b,i,j");

        // D_oovv_bbbb += +0.50 r0 t2_bbbb(a,c,i,j) t2_abab(d,b,k,l) l2_abab(k,l,d,c)
        //             += +0.50 r0 t2_bbbb(a,c,i,j) t2_abab(d,b,l,k) l2_abab(l,k,d,c)
        //             += +0.50 r0_1 t2_bbbb(a,c,i,j) t2_abab(d,b,k,l) l2_1_abab(k,l,d,c)
        //             += +0.50 r0_1 t2_bbbb(a,c,i,j) t2_abab(d,b,l,k) l2_1_abab(l,k,d,c)
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["86_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["86_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v2L2  = o3v3L2
        //  mems: o2v2L2  = o2v2L2
        tmps_["87_aaaa_LLovvo"]("L,R,i,a,b,j")  = l2_1["aaaa"]("L,k,i,a,c") * r2["aaaa"]("R,b,c,k,j");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_aaaa(a,d,l,i) t2_1_aaaa(b,c,k,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["87_aaaa_LLovvo"]("L,R,k,c,a,i") * t2_1["aaaa"]("b,c,k,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_aaaa(a,d,l,i) t1_aa(c,j) t1_1_aa(b,k) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["87_aaaa_LLovvo"]("L,R,k,c,a,i") * t1["aa"]("c,j") * t1_1["aa"]("b,k");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r2_aaaa(a,d,l,i) t2_1_abab(c,b,k,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["87_aaaa_LLovvo"]("L,R,k,c,a,i") * t2_1["abab"]("c,b,k,j");

        // flops: o2v2L2  = o3v3L2
        //  mems: o2v2L2  = o2v2L2
        tmps_["88_aaaa_LLovvo"]("L,R,i,a,b,j")  = l2_1["abab"]("L,i,k,a,c") * r2["abab"]("R,b,c,j,k");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_abab(a,d,i,l) t1_aa(c,j) t1_1_aa(b,k) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["88_aaaa_LLovvo"]("L,R,k,c,a,i") * t1["aa"]("c,j") * t1_1["aa"]("b,k");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o2v2L2  = o3v3L2 o3v3L2 o2v2L2
        //  mems: o2v2L2  = o2v2L2 o2v2L2 o2v2L2
        tmps_["89_aaaa_LLovvo"]("L,R,i,a,b,j")  = l2["aaaa"]("L,k,i,a,c") * r2["aaaa"]("R,b,c,k,j");
        tmps_["89_aaaa_LLovvo"]("L,R,i,a,b,j") += l2_1["aaaa"]("L,k,i,a,c") * r2_1["aaaa"]("R,b,c,k,j");

        // D_ovvo_aaaa  = -1.00 r2_1_aaaa(b,c,k,i) l2_1_aaaa(k,j,a,c)
        //             += -1.00 r2_aaaa(b,c,k,i) l2_aaaa(k,j,a,c)
        D_ovvo_aaaa("L,R,i,a,b,j")  = -1.00 * tmps_["89_aaaa_LLovvo"]("L,R,j,a,b,i");

        // D_ovov_aaaa  = +1.00 r2_1_aaaa(b,c,k,i) l2_1_aaaa(k,j,a,c)
        //             += +1.00 r2_aaaa(b,c,k,i) l2_aaaa(k,j,a,c)
        D_ovov_aaaa("L,R,i,a,j,b")  = tmps_["89_aaaa_LLovvo"]("L,R,j,a,b,i");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_1_aaaa(a,d,l,i) t2_aaaa(b,c,k,j) l2_1_aaaa(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r2_aaaa(a,d,l,i) t2_aaaa(b,c,k,j) l2_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["89_aaaa_LLovvo"]("L,R,k,c,a,i") * t2["aaaa"]("b,c,k,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r2_1_aaaa(a,d,l,i) t2_abab(c,b,k,j) l2_1_aaaa(l,k,c,d)
        //             += -1.00 r2_aaaa(a,d,l,i) t2_abab(c,b,k,j) l2_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["89_aaaa_LLovvo"]("L,R,k,c,a,i") * t2["abab"]("c,b,k,j");

        // flops: o2v2L2  = o3v3L2 o3v3L2 o2v2L2
        //  mems: o2v2L2  = o2v2L2 o2v2L2 o2v2L2
        tmps_["90_aaaa_LLovvo"]("L,R,i,a,b,j")  = l2["abab"]("L,i,k,a,c") * r2["abab"]("R,b,c,j,k");
        tmps_["90_aaaa_LLovvo"]("L,R,i,a,b,j") += l2_1["abab"]("L,i,k,a,c") * r2_1["abab"]("R,b,c,j,k");

        // D_ovov_aaaa += +1.00 r2_1_abab(b,c,i,k) l2_1_abab(j,k,a,c)
        //             += +1.00 r2_abab(b,c,i,k) l2_abab(j,k,a,c)
        D_ovov_aaaa("L,R,i,a,j,b") += tmps_["90_aaaa_LLovvo"]("L,R,j,a,b,i");

        // D_ovvo_aaaa += -1.00 r2_1_abab(b,c,i,k) l2_1_abab(j,k,a,c)
        //             += -1.00 r2_abab(b,c,i,k) l2_abab(j,k,a,c)
        D_ovvo_aaaa("L,R,i,a,b,j") -= tmps_["90_aaaa_LLovvo"]("L,R,j,a,b,i");

        // flops: o3v1L2  = o3v2L2 o3v2L2 o3v1L2 o4v1L1 o4v1L2 o3v1L2 o3v2L2 o3v1L2 o3v2L2 o3v1L2
        //  mems: o3v1L2  = o3v1L2 o3v1L2 o3v1L2 o4v0L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        tmps_["91_aaaa_LLoovo"]("L,R,i,j,a,k")  = t1["aa"]("b,i") * tmps_["89_aaaa_LLovvo"]("L,R,j,b,a,k");
        tmps_["91_aaaa_LLoovo"]("L,R,i,j,a,k") += t1_1["aa"]("b,i") * tmps_["88_aaaa_LLovvo"]("L,R,j,b,a,k");
        tmps_["91_aaaa_LLoovo"]("L,R,i,j,a,k") += t1_1["aa"]("b,k") * tmps_["22_aaaa_Looov"]("L,i,l,j,b") * r1["aa"]("R,a,l");
        tmps_["91_aaaa_LLoovo"]("L,R,i,j,a,k") += t1["aa"]("b,i") * tmps_["90_aaaa_LLovvo"]("L,R,j,b,a,k");
        tmps_["91_aaaa_LLoovo"]("L,R,i,j,a,k") += t1_1["aa"]("b,i") * tmps_["87_aaaa_LLovvo"]("L,R,j,b,a,k");

        // D_ooov_aaaa += -1.00 P(i,j) r2_1_aaaa(a,c,l,i) t1_aa(b,j) l2_1_aaaa(l,k,b,c)
        //             += -1.00 P(i,j) r2_aaaa(a,c,l,i) t1_aa(b,j) l2_aaaa(l,k,b,c)
        //             += -1.00 P(i,j) r2_abab(a,c,i,l) t1_1_aa(b,j) l2_1_abab(k,l,b,c)
        //             += -1.00 P(i,j) r1_aa(a,l) t1_aa(c,j) t1_1_aa(b,i) l2_1_aaaa(l,k,c,b)
        //             += -1.00 P(i,j) r2_1_abab(a,c,i,l) t1_aa(b,j) l2_1_abab(k,l,b,c)
        //             += -1.00 P(i,j) r2_abab(a,c,i,l) t1_aa(b,j) l2_abab(k,l,b,c)
        //             += -1.00 P(i,j) r2_aaaa(a,c,l,i) t1_1_aa(b,j) l2_1_aaaa(l,k,b,c)
        D_ooov_aaaa("L,R,i,j,k,a") -= tmps_["91_aaaa_LLoovo"]("L,R,j,k,a,i");
        D_ooov_aaaa("L,R,i,j,k,a") += tmps_["91_aaaa_LLoovo"]("L,R,i,k,a,j");

        // D_oovo_aaaa += +1.00 P(i,j) r2_1_aaaa(a,c,l,i) t1_aa(b,j) l2_1_aaaa(l,k,b,c)
        //             += +1.00 P(i,j) r2_aaaa(a,c,l,i) t1_aa(b,j) l2_aaaa(l,k,b,c)
        //             += +1.00 P(i,j) r2_abab(a,c,i,l) t1_1_aa(b,j) l2_1_abab(k,l,b,c)
        //             += +1.00 P(i,j) r1_aa(a,l) t1_aa(c,j) t1_1_aa(b,i) l2_1_aaaa(l,k,c,b)
        //             += +1.00 P(i,j) r2_1_abab(a,c,i,l) t1_aa(b,j) l2_1_abab(k,l,b,c)
        //             += +1.00 P(i,j) r2_abab(a,c,i,l) t1_aa(b,j) l2_abab(k,l,b,c)
        //             += +1.00 P(i,j) r2_aaaa(a,c,l,i) t1_1_aa(b,j) l2_1_aaaa(l,k,b,c)
        D_oovo_aaaa("L,R,i,j,a,k") += tmps_["91_aaaa_LLoovo"]("L,R,j,k,a,i");
        D_oovo_aaaa("L,R,i,j,a,k") -= tmps_["91_aaaa_LLoovo"]("L,R,i,k,a,j");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_1_aaaa(a,d,l,i) t1_aa(b,k) t1_aa(c,j) l2_1_aaaa(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r2_aaaa(a,d,l,i) t1_aa(b,k) t1_aa(c,j) l2_aaaa(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r2_abab(a,d,i,l) t1_aa(b,k) t1_1_aa(c,j) l2_1_abab(k,l,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_aa(a,l) t1_aa(b,k) t1_aa(d,j) t1_1_aa(c,i) l2_1_aaaa(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r2_1_abab(a,d,i,l) t1_aa(b,k) t1_aa(c,j) l2_1_abab(k,l,c,d)
        //             += +1.00 P(i,j) P(a,b) r2_abab(a,d,i,l) t1_aa(b,k) t1_aa(c,j) l2_abab(k,l,c,d)
        //             += +1.00 P(i,j) P(a,b) r2_aaaa(a,d,l,i) t1_aa(b,k) t1_1_aa(c,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("b,k") * tmps_["91_aaaa_LLoovo"]("L,R,j,k,a,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();
        tmps_["91_aaaa_LLoovo"].~TArrayD();

        // flops: o5v0L2  = o2v0L1 o3v2L1 o3v2L1 o3v0L1 o3v2L1 o3v2L1 o3v0L1 o3v0L1 o5v0L2
        //  mems: o5v0L2  = o2v0L1 o2v0L1 o2v0L1 o3v0L1 o2v0L1 o2v0L1 o3v0L1 o3v0L1 o5v0L2
        tmps_["92_bbbbb_LLooooo"]("R,L,i,j,k,l,m")  = Id["bb_oo"]("i,j") * r0("R") * (2.00 * (t2["abab"]("d,c,n,k") * l2["abab"]("L,n,l,d,c") + -1.00 * t2_1["abab"]("e,b,o,k") * l2_1["abab"]("L,o,m,e,b")) + -1.00 * l2_1["bbbb"]("L,j,m,a,b") * t2_1["bbbb"]("a,b,j,k") + l2["bbbb"]("L,m,l,b,c") * t2["bbbb"]("b,c,m,k"));

        // D_oooo_bbbb += +0.50 P(i,j) d_bb(j,l) r0 t2_bbbb(b,a,m,i) l2_bbbb(m,k,b,a)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t1_bb(a,m) t2_1_bbbb(c,b,l,i) l2_1_bbbb(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,l) r0 t2_abab(b,a,m,i) l2_abab(m,k,b,a)
        //             += +0.50 P(i,j) d_bb(j,l) r0 t2_abab(a,b,m,i) l2_abab(m,k,a,b)
        //             += -0.50 P(i,j) d_bb(j,k) r1_bb(a,m) t2_1_abab(c,b,l,i) l2_1_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r1_bb(a,m) t2_1_abab(b,c,l,i) l2_1_abab(l,m,b,c)
        D_oooo_bbbb("L,R,i,j,k,l") += 0.50 * tmps_["92_bbbbb_LLooooo"]("R,L,j,l,i,k,m");
        D_oooo_bbbb("L,R,i,j,k,l") -= 0.50 * tmps_["92_bbbbb_LLooooo"]("R,L,i,l,j,k,m");

        // D_oooo_bbbb += -0.50 P(i,j) d_bb(j,k) r0 t2_bbbb(b,a,m,i) l2_bbbb(m,l,b,a)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t1_bb(a,m) t2_1_bbbb(c,b,l,i) l2_1_bbbb(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_abab(b,a,m,i) l2_abab(m,l,b,a)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_abab(a,b,m,i) l2_abab(m,l,a,b)
        //             += +0.50 P(i,j) d_bb(j,k) r1_bb(a,m) t2_1_abab(c,b,l,i) l2_1_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r1_bb(a,m) t2_1_abab(b,c,l,i) l2_1_abab(l,m,b,c)
        D_oooo_bbbb("L,R,i,j,k,l") -= 0.50 * tmps_["92_bbbbb_LLooooo"]("R,L,j,k,i,l,m");
        D_oooo_bbbb("L,R,i,j,k,l") += 0.50 * tmps_["92_bbbbb_LLooooo"]("R,L,i,k,j,l,m");
        tmps_["92_bbbbb_LLooooo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["93_bbaa_Lvoov"]("L,a,i,j,b")  = t2["bbbb"]("a,c,k,i") * l2_1["abab"]("L,j,k,b,c");

        // D_oovv_abab += -1.00 r1_aa(d,i) t2_bbbb(b,c,l,j) t1_1_aa(a,k) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["93_bbaa_Lvoov"]("L,b,j,k,d") * r1["aa"]("R,d,i") * t1_1["aa"]("a,k");

        // D_oovv_abab += -1.00 r1_1_aa(d,i) t1_aa(a,l) t2_bbbb(b,c,k,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["93_bbaa_Lvoov"]("L,b,j,l,d") * r1_1["aa"]("R,d,i") * t1["aa"]("a,l");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["94_baab_Lvoov"]("L,a,i,j,b")  = t2["abab"]("c,a,i,k") * l2_1["abab"]("L,j,k,c,b");

        // D_oovv_abab += -1.00 r2_1_abab(a,d,l,j) t2_abab(c,b,i,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["94_baab_Lvoov"]("L,b,i,l,d") * r2_1["abab"]("R,a,d,l,j");

        // D_oovv_abab += -1.00 r0 t2_abab(d,b,i,l) t2_1_abab(a,c,k,j) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2_1["abab"]("a,c,k,j") * tmps_["94_baab_Lvoov"]("L,b,i,k,c") * r0("R");

        // D_oovv_abab += -1.00 r0_1 t2_abab(a,c,k,j) t2_abab(d,b,i,l) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,c,k,j") * tmps_["94_baab_Lvoov"]("L,b,i,k,c") * r0_1("R");

        // D_oovv_abab += -1.00 r1_bb(d,j) t2_abab(c,b,i,l) t1_1_aa(a,k) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["94_baab_Lvoov"]("L,b,i,k,d") * r1["bb"]("R,d,j") * t1_1["aa"]("a,k");

        // D_oovv_abab += -1.00 r1_1_bb(d,j) t1_aa(a,l) t2_abab(c,b,i,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["94_baab_Lvoov"]("L,b,i,l,d") * r1_1["bb"]("R,d,j") * t1["aa"]("a,l");

        // flops: o3v1L1  = o3v2L1 o3v2L1 o3v1L1
        //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1
        tmps_["95_baab_Lvooo"]("L,a,i,j,k")  = tmps_["94_baab_Lvoov"]("L,a,i,j,b") * t1["bb"]("b,k");
        tmps_["95_baab_Lvooo"]("L,a,i,j,k") += tmps_["93_bbaa_Lvoov"]("L,a,k,j,c") * t1["aa"]("c,i");

        // D_oovv_abab += -1.00 r1_1_aa(a,l) t2_abab(c,b,i,k) t1_bb(d,j) l2_1_abab(l,k,c,d)
        //             += -1.00 r1_1_aa(a,l) t2_bbbb(b,c,k,j) t1_aa(d,i) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["95_baab_Lvooo"]("L,b,i,l,j") * r1_1["aa"]("R,a,l");

        // D_oovv_abab += -1.00 r0_1 t1_aa(a,l) t2_abab(c,b,i,k) t1_bb(d,j) l2_1_abab(l,k,c,d)
        //             += -1.00 r0_1 t1_aa(a,l) t2_bbbb(b,c,k,j) t1_aa(d,i) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["aa"]("a,l") * tmps_["95_baab_Lvooo"]("L,b,i,l,j") * r0_1("R");

        // D_oovv_abab += -1.00 r0 t2_abab(c,b,i,l) t1_bb(d,j) t1_1_aa(a,k) l2_1_abab(k,l,c,d)
        //             += -1.00 r0 t2_bbbb(b,c,l,j) t1_aa(d,i) t1_1_aa(a,k) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1_1["aa"]("a,k") * tmps_["95_baab_Lvooo"]("L,b,i,k,j") * r0("R");
        tmps_["95_baab_Lvooo"].~TArrayD();

        // flops: o2v0L1  = o3v2L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["96_bb_Loo"]("L,i,j")  = t2["abab"]("a,b,k,i") * l2_1["abab"]("L,k,j,a,b");

        // D_oovv_abab += +0.50 r2_1_abab(a,b,i,l) t2_abab(d,c,k,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r2_1_abab(a,b,i,l) t2_abab(c,d,k,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["96_bb_Loo"]("L,j,l") * r2_1["abab"]("R,a,b,i,l");

        // D_oovv_abab += +0.50 r0 t2_abab(d,c,l,j) t2_1_abab(a,b,i,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r0 t2_abab(c,d,l,j) t2_1_abab(a,b,i,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2_1["abab"]("a,b,i,k") * tmps_["96_bb_Loo"]("L,j,k") * r0("R");

        // D_oovv_abab += +0.50 r0_1 t2_abab(a,b,i,l) t2_abab(d,c,k,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r0_1 t2_abab(a,b,i,l) t2_abab(c,d,k,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,b,i,l") * tmps_["96_bb_Loo"]("L,j,l") * r0_1("R");

        // D_oovv_bbbb += +0.50 P(i,j) r2_1_bbbb(b,a,l,i) t2_abab(d,c,k,j) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) r2_1_bbbb(b,a,l,i) t2_abab(c,d,k,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["96_bb_Loo"]("L,j,l") * r2_1["bbbb"]("R,b,a,l,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += -0.50 P(i,j) r0_1 t2_bbbb(b,a,l,j) t2_abab(d,c,k,i) l2_1_abab(k,l,d,c)
        //             += -0.50 P(i,j) r0_1 t2_bbbb(b,a,l,j) t2_abab(c,d,k,i) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["bbbb"]("b,a,l,j") * tmps_["96_bb_Loo"]("L,i,l") * r0_1("R");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +0.50 P(i,j) r0 t2_abab(d,c,l,j) t2_1_bbbb(b,a,k,i) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) r0 t2_abab(c,d,l,j) t2_1_bbbb(b,a,k,i) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2_1["bbbb"]("b,a,k,i") * tmps_["96_bb_Loo"]("L,j,k") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v0L1  = o3v2L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["97_bb_Loo"]("L,i,j")  = t2["bbbb"]("a,b,k,i") * l2_1["bbbb"]("L,j,k,a,b");

        // D_oovv_abab += -0.50 r2_1_abab(a,b,i,l) t2_bbbb(d,c,k,j) l2_1_bbbb(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * tmps_["97_bb_Loo"]("L,j,l") * r2_1["abab"]("R,a,b,i,l");

        // D_oovv_abab += -0.50 r0 t2_bbbb(d,c,l,j) t2_1_abab(a,b,i,k) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * t2_1["abab"]("a,b,i,k") * tmps_["97_bb_Loo"]("L,j,k") * r0("R");

        // D_oovv_bbbb += -0.50 P(i,j) r2_1_bbbb(b,a,l,i) t2_bbbb(d,c,k,j) l2_1_bbbb(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["97_bb_Loo"]("L,j,l") * r2_1["bbbb"]("R,b,a,l,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += -0.50 P(i,j) r0 t2_bbbb(d,c,l,j) t2_1_bbbb(b,a,k,i) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t2_1["bbbb"]("b,a,k,i") * tmps_["97_bb_Loo"]("L,j,k") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o0v2L2  = o2v3L2 o2v3L2 o0v2L2
        //  mems: o0v2L2  = o0v2L2 o0v2L2 o0v2L2
        tmps_["98_bb_LLvv"]("L,R,a,b")  = 0.50 * l2_1["bbbb"]("L,k,j,a,d") * r2["bbbb"]("R,b,d,k,j");
        tmps_["98_bb_LLvv"]("L,R,a,b") += l2_1["abab"]("L,i,j,c,a") * r2["abab"]("R,c,b,i,j");

        // D_oovv_abab += +0.50 r2_abab(d,b,l,k) t2_1_abab(a,c,i,j) l2_1_abab(l,k,d,c)
        //             += +0.50 r2_abab(d,b,k,l) t2_1_abab(a,c,i,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r2_bbbb(b,d,l,k) t2_1_abab(a,c,i,j) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2_1["abab"]("a,c,i,j") * tmps_["98_bb_LLvv"]("L,R,c,b");
    }
} // hilbert
#endif