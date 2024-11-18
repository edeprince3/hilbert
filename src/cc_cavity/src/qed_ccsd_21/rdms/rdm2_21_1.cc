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


    void EOM_EE_QED_RDM_21::rdm2_21_1() {

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


        // D_ooov_aaaa  = -0.50 P(i,j) r1_aa(a,i) t2_1_aaaa(c,b,l,j) l2_1_aaaa(l,k,c,b)
        //             += +0.50 P(i,j) d_aa(j,l) r0 t2_aaaa(b,a,m,i) l2_aaaa(m,k,b,a)
        // flops: o3v1L2  = o3v2L1 o3v2L1 o3v0L1 o3v1L2
        //  mems: o3v1L2  = o2v0L1 o2v0L1 o3v0L1 o2v1L2
        tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k")  = 0.50 * (t2_1["aaaa"]("c,b,l,j") * l2_1["aaaa"]("L,l,k,c,b") + -1.00 * t2["aaaa"]("b,a,m,i") * l2["aaaa"]("L,m,k,b,a")) * r1["aa"]("R,a,i");
        D_ooov_aaaa("L,R,i,j,k,a")  = -1.00 * tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k");
        D_ooov_aaaa("L,R,i,j,k,a") += tmps_["perm_aaaa_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_aaaa_LLvooo"].~TArrayD();

        // D_oovo_aaaa  = +0.50 P(i,j) r1_aa(a,i) t2_1_abab(c,b,j,l) l2_1_abab(k,l,c,b)
        //             += +0.50 P(i,j) r1_aa(a,i) t2_1_abab(b,c,j,l) l2_1_abab(k,l,b,c)
        //             += +0.50 r2_abab(a,b,l,j) t2_abab(d,c,i,k) l2_abab(l,k,d,c)
        //             += +0.50 r2_abab(a,b,l,j) t2_abab(c,d,i,k) l2_abab(l,k,c,d)
        // flops: o3v1L2  = o3v2L1 o3v2L1 o4v0L1 o4v1L2
        //  mems: o3v1L2  = o2v0L1 o2v0L1 o4v0L1 o3v1L2
        tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k")  = (t2_1["abab"]("c,b,j,l") * l2_1["abab"]("L,k,l,c,b") + t2["abab"]("d,c,i,k") * l2["abab"]("L,l,k,d,c")) * r1["aa"]("R,a,i");
        D_oovo_aaaa("L,R,i,j,a,k")  = tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k");
        D_oovo_aaaa("L,R,i,j,a,k") -= tmps_["perm_aaaa_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_aaaa_LLvooo"].~TArrayD();

        // D_oooo_aaaa  = -1.00 P(i,j) r1_aa(b,i) t1_1_aa(a,j) l2_1_aaaa(l,k,a,b)
        // flops: o4v0L2  = o3v2L1 o4v1L2
        //  mems: o4v0L2  = o3v1L1 o4v0L2
        tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l")  = t1_1["aa"]("a,j") * l2_1["aaaa"]("L,l,k,a,b") * r1["aa"]("R,b,i");
        D_oooo_aaaa("L,R,i,j,k,l")  = -1.00 * tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l");
        D_oooo_aaaa("L,R,i,j,k,l") += tmps_["perm_aaaa_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_aaaa_LLoooo"].~TArrayD();

        // D_ooov_bbbb  = -0.50 P(i,j) r1_bb(a,i) t2_bbbb(c,b,l,j) l2_bbbb(l,k,c,b)
        //             += +0.50 r0 t1_aa(a,i) t1_bb(b,l) t2_1_bbbb(d,c,k,j) l2_1_bbbb(k,l,d,c)
        // flops: o3v1L2  = o3v2L1 o3v2L1 o3v0L1 o4v1L2
        //  mems: o3v1L2  = o2v0L1 o2v0L1 o3v0L1 o4v1L2
        tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k")  = 0.50 * (t2["bbbb"]("c,b,l,j") * l2["bbbb"]("L,l,k,c,b") + -1.00 * t2_1["bbbb"]("d,c,k,j") * l2_1["bbbb"]("L,k,l,d,c")) * r1["bb"]("R,a,i");
        D_ooov_bbbb("L,R,i,j,k,a")  = -1.00 * tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k");
        D_ooov_bbbb("L,R,i,j,k,a") += tmps_["perm_bbbb_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_bbbb_LLvooo"].~TArrayD();

        // D_vvvv_aaaa  = -1.00 P(c,d) r1_aa(c,j) t1_1_aa(d,i) l2_1_aaaa(j,i,a,b)
        // flops: o0v4L2  = o2v3L1 o1v4L2
        //  mems: o0v4L2  = o1v3L1 o0v4L2
        tmps_["perm_aaaa_LLvvvv"]("L,R,a,b,c,d")  = l2_1["aaaa"]("L,j,i,a,b") * t1_1["aa"]("d,i") * r1["aa"]("R,c,j");
        D_vvvv_aaaa("L,R,a,b,c,d")  = -1.00 * tmps_["perm_aaaa_LLvvvv"]("L,R,a,b,c,d");
        D_vvvv_aaaa("L,R,a,b,c,d") += tmps_["perm_aaaa_LLvvvv"]("L,R,a,b,d,c");
        tmps_["perm_aaaa_LLvvvv"].~TArrayD();

        // D_oovv_abab  = -1.00 r2_abab(a,b,i,j) l0
        // flops: o2v2L2  = o2v2L2
        //  mems: o2v2L2  = o2v2L2
        D_oovv_abab("L,R,i,j,a,b")  = -1.00 * l0("L") * r2["abab"]("R,a,b,i,j");

        // D_oovo_abab  = -0.50 r1_aa(a,j) t2_abab(c,b,l,i) l2_abab(l,k,c,b)
        //             += -0.50 r1_aa(a,j) t2_abab(b,c,l,i) l2_abab(l,k,b,c)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,l) t1_bb(b,j) t2_1_abab(d,c,k,i) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,l) t1_bb(b,j) t2_1_abab(c,d,k,i) l2_1_abab(k,l,c,d)
        // flops: o3v1L2  = o3v2L1 o3v2L1 o3v0L1 o4v1L2
        //  mems: o3v1L2  = o2v0L1 o2v0L1 o3v0L1 o4v1L2
        D_oovo_abab("L,R,i,j,a,k")  = -1.00 * (t2["abab"]("c,b,l,i") * l2["abab"]("L,l,k,c,b") + -1.00 * t2_1["abab"]("d,c,k,i") * l2_1["abab"]("L,k,l,d,c")) * r1["aa"]("R,a,j");

        // D_vvvv_abab  = -1.00 r1_bb(d,j) t1_1_aa(c,i) l2_1_abab(i,j,a,b)
        // flops: o0v4L2  = o2v3L2 o1v4L2
        //  mems: o0v4L2  = o1v3L2 o0v4L2
        D_vvvv_abab("L,R,a,b,c,d")  = -1.00 * l2_1["abab"]("L,i,j,a,b") * r1["bb"]("R,d,j") * t1_1["aa"]("c,i");

        // D_oovo_bbbb  = +0.50 P(i,j) r1_bb(a,i) t2_bbbb(c,b,l,j) l2_bbbb(l,k,c,b)
        //             += +0.50 r1_aa(a,i) t1_bb(b,l) t2_1_bbbb(d,c,k,j) l2_1_bbbb(k,l,d,c)
        // flops: o3v1L2  = o3v2L1 o3v2L1 o3v0L1 o4v1L2
        //  mems: o3v1L2  = o2v0L1 o2v0L1 o3v0L1 o4v1L2
        tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k")  = 0.50 * (t2["bbbb"]("c,b,l,j") * l2["bbbb"]("L,l,k,c,b") + t2_1["bbbb"]("d,c,k,j") * l2_1["bbbb"]("L,k,l,d,c")) * r1["bb"]("R,a,i");
        D_oovo_bbbb("L,R,i,j,a,k")  = tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k");
        D_oovo_bbbb("L,R,i,j,a,k") -= tmps_["perm_bbbb_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_bbbb_LLvooo"].~TArrayD();

        // D_oovv_aaaa  = +1.00 r2_aaaa(b,a,i,j) l0
        // flops: o2v2L2  = o2v2L2
        //  mems: o2v2L2  = o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b")  = l0("L") * r2["aaaa"]("R,b,a,i,j");

        // D_oovv_bbbb  = +1.00 r2_1_bbbb(b,a,i,j) l0_1
        // flops: o2v2L2  = o2v2L2
        //  mems: o2v2L2  = o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b")  = l0_1("L") * r2_1["bbbb"]("R,b,a,i,j");

        // D_ooov_abab  = +0.50 r1_aa(a,j) t2_abab(c,b,l,i) l2_abab(l,k,c,b)
        //             += +0.50 r1_aa(a,j) t2_abab(b,c,l,i) l2_abab(l,k,b,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_bb(a,l) t1_bb(b,j) t2_1_abab(d,c,k,i) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_bb(a,l) t1_bb(b,j) t2_1_abab(c,d,k,i) l2_1_abab(k,l,c,d)
        // flops: o3v1L2  = o3v2L1 o3v2L1 o3v0L1 o4v1L2
        //  mems: o3v1L2  = o2v0L1 o2v0L1 o3v0L1 o4v1L2
        D_ooov_abab("L,R,i,j,k,a")  = (t2["abab"]("c,b,l,i") * l2["abab"]("L,l,k,c,b") + t2_1["abab"]("d,c,k,i") * l2_1["abab"]("L,k,l,d,c")) * r1["aa"]("R,a,j");

        // D_vvvv_bbbb  = -1.00 P(c,d) r1_bb(c,j) t1_1_bb(d,i) l2_1_bbbb(j,i,a,b)
        // flops: o0v4L2  = o2v3L1 o1v4L2
        //  mems: o0v4L2  = o1v3L1 o0v4L2
        tmps_["perm_bbbb_LLvvvv"]("L,R,a,b,c,d")  = l2_1["bbbb"]("L,j,i,a,b") * t1_1["bb"]("d,i") * r1["bb"]("R,c,j");
        D_vvvv_bbbb("L,R,a,b,c,d")  = -1.00 * tmps_["perm_bbbb_LLvvvv"]("L,R,a,b,c,d");
        D_vvvv_bbbb("L,R,a,b,c,d") += tmps_["perm_bbbb_LLvvvv"]("L,R,a,b,d,c");
        tmps_["perm_bbbb_LLvvvv"].~TArrayD();

        // D_oooo_bbbb  = -1.00 P(i,j) r1_bb(b,i) t1_1_bb(a,j) l2_1_bbbb(l,k,a,b)
        // flops: o4v0L2  = o3v2L1 o4v1L2
        //  mems: o4v0L2  = o3v1L1 o4v0L2
        tmps_["perm_bbbb_LLoooo"]("L,R,i,j,k,l")  = t1_1["bb"]("a,j") * l2_1["bbbb"]("L,l,k,a,b") * r1["bb"]("R,b,i");
        D_oooo_bbbb("L,R,i,j,k,l")  = -1.00 * tmps_["perm_bbbb_LLoooo"]("L,R,i,j,k,l");
        D_oooo_bbbb("L,R,i,j,k,l") += tmps_["perm_bbbb_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_bbbb_LLoooo"].~TArrayD();

        // D_oooo_aaaa += -0.50 P(i,j) d_aa(j,k) r0 t2_1_aaaa(b,a,m,i) l2_1_aaaa(m,l,b,a)
        //             += +0.50 P(i,j) r1_aa(a,i) t2_aaaa(c,b,l,j) l2_aaaa(l,k,c,b)
        // flops: o4v0L2 += o3v2L1 o3v2L1 o4v0L1 o4v0L1 o2v0L2
        //  mems: o4v0L2 += o2v0L1 o2v0L1 o4v0L1 o2v0L1 o2v0L2
        tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l")  = 0.50 * (t2_1["aaaa"]("b,a,m,i") * l2_1["aaaa"]("L,m,l,b,a") + -1.00 * t2["aaaa"]("c,b,l,j") * l2["aaaa"]("L,l,k,c,b")) * Id["aa_oo"]("j,k") * r0("R");
        D_oooo_aaaa("L,R,i,j,k,l") -= tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l");
        D_oooo_aaaa("L,R,i,j,k,l") += tmps_["perm_aaaa_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_aaaa_LLoooo"].~TArrayD();

        // D_oooo_aaaa += +0.50 P(i,j) d_aa(j,l) r0 t2_1_aaaa(b,a,m,i) l2_1_aaaa(m,k,b,a)
        //             += -0.50 P(i,j) r1_aa(a,i) t2_aaaa(c,b,l,j) l2_aaaa(l,k,c,b)
        // flops: o4v0L2 += o3v2L1 o3v2L1 o3v0L1 o4v0L1 o3v0L2
        //  mems: o4v0L2 += o2v0L1 o2v0L1 o3v0L1 o3v0L1 o3v0L2
        tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l")  = 0.50 * (t2_1["aaaa"]("b,a,m,i") * l2_1["aaaa"]("L,m,k,b,a") + -1.00 * t2["aaaa"]("c,b,l,j") * l2["aaaa"]("L,l,k,c,b")) * Id["aa_oo"]("j,l") * r0("R");
        D_oooo_aaaa("L,R,i,j,k,l") += tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l");
        D_oooo_aaaa("L,R,i,j,k,l") -= tmps_["perm_aaaa_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_aaaa_LLoooo"].~TArrayD();

        // D_ooov_aaaa += -0.50 P(i,j) r1_aa(a,i) t2_1_abab(c,b,j,l) l2_1_abab(k,l,c,b)
        //             += -0.50 P(i,j) r1_aa(a,i) t2_1_abab(b,c,j,l) l2_1_abab(k,l,b,c)
        //             += +0.50 r0 t1_aa(a,l) t1_bb(b,j) t2_abab(d,c,i,k) l2_abab(l,k,d,c)
        //             += +0.50 r0 t1_aa(a,l) t1_bb(b,j) t2_abab(c,d,i,k) l2_abab(l,k,c,d)
        // flops: o3v1L2 += o3v2L1 o3v2L1 o4v0L1 o4v1L2
        //  mems: o3v1L2 += o2v0L1 o2v0L1 o4v0L1 o3v1L2
        tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k")  = (t2_1["abab"]("c,b,j,l") * l2_1["abab"]("L,k,l,c,b") + -1.00 * t2["abab"]("d,c,i,k") * l2["abab"]("L,l,k,d,c")) * r1["aa"]("R,a,i");
        D_ooov_aaaa("L,R,i,j,k,a") -= tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k");
        D_ooov_aaaa("L,R,i,j,k,a") += tmps_["perm_aaaa_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_aaaa_LLvooo"].~TArrayD();

        // D_ooov_abab += +0.50 d_bb(i,k) r1_aa(a,m) t2_1_abab(c,b,j,l) l2_1_abab(m,l,c,b)
        //             += +0.50 d_bb(i,k) r1_aa(a,m) t2_1_abab(b,c,j,l) l2_1_abab(m,l,b,c)
        //             += +0.50 P(i,j) r2_aaaa(b,a,l,i) t2_abab(d,c,j,k) l2_abab(l,k,d,c)
        //             += +0.50 P(i,j) r2_aaaa(b,a,l,i) t2_abab(c,d,j,k) l2_abab(l,k,c,d)
        // flops: o3v1L2 += o3v2L1 o3v2L1 o3v0L1 o3v1L2 o4v1L2
        //  mems: o3v1L2 += o2v0L1 o2v0L1 o3v0L1 o2v1L2 o4v1L2
        D_ooov_abab("L,R,i,j,k,a") += (t2_1["abab"]("c,b,j,l") * l2_1["abab"]("L,m,l,c,b") + t2["abab"]("d,c,j,k") * l2["abab"]("L,l,k,d,c")) * r1["aa"]("R,a,m") * Id["bb_oo"]("i,k");

        // D_ooov_bbbb += -0.50 P(i,j) r1_bb(a,i) t2_abab(c,b,l,j) l2_abab(l,k,c,b)
        //             += -0.50 P(i,j) r1_bb(a,i) t2_abab(b,c,l,j) l2_abab(l,k,b,c)
        //             += +0.50 r0 t1_aa(a,i) t1_bb(b,l) t2_1_abab(d,c,k,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r0 t1_aa(a,i) t1_bb(b,l) t2_1_abab(c,d,k,j) l2_1_abab(k,l,c,d)
        // flops: o3v1L2 += o3v2L1 o3v2L1 o3v0L1 o4v1L2
        //  mems: o3v1L2 += o2v0L1 o2v0L1 o3v0L1 o4v1L2
        tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k")  = (t2["abab"]("c,b,l,j") * l2["abab"]("L,l,k,c,b") + -1.00 * t2_1["abab"]("d,c,k,j") * l2_1["abab"]("L,k,l,d,c")) * r1["bb"]("R,a,i");
        D_ooov_bbbb("L,R,i,j,k,a") -= tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k");
        D_ooov_bbbb("L,R,i,j,k,a") += tmps_["perm_bbbb_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_bbbb_LLvooo"].~TArrayD();

        // D_ooov_bbbb += -0.50 P(i,j) d_bb(j,k) r0 t1_bb(a,m) t2_bbbb(c,b,l,i) l2_bbbb(l,m,c,b)
        //             += +0.50 d_aa(i,k) r0 t2_1_bbbb(b,a,m,j) l2_1_bbbb(m,l,b,a)
        // flops: o3v1L2 += o3v2L1 o3v2L1 o4v0L1 o5v0L1 o4v1L1 o3v1L2
        //  mems: o3v1L2 += o2v0L1 o2v0L1 o4v0L1 o4v0L1 o3v1L1 o3v1L2
        tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k")  = 0.50 * (t2["bbbb"]("c,b,l,i") * l2["bbbb"]("L,l,m,c,b") + -1.00 * t2_1["bbbb"]("b,a,m,j") * l2_1["bbbb"]("L,m,l,b,a")) * Id["bb_oo"]("j,k") * t1["bb"]("a,m") * r0("R");
        D_ooov_bbbb("L,R,i,j,k,a") -= tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k");
        D_ooov_bbbb("L,R,i,j,k,a") += tmps_["perm_bbbb_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_bbbb_LLvooo"].~TArrayD();

        // D_oovo_aaaa += +0.50 P(i,j) r1_aa(a,i) t2_1_aaaa(c,b,l,j) l2_1_aaaa(l,k,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_aaaa(b,a,m,i) l2_aaaa(m,l,b,a)
        // flops: o3v1L2 += o3v2L1 o3v2L1 o4v0L1 o4v1L2
        //  mems: o3v1L2 += o2v0L1 o2v0L1 o4v0L1 o3v1L2
        tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k")  = 0.50 * (t2_1["aaaa"]("c,b,l,j") * l2_1["aaaa"]("L,l,k,c,b") + -1.00 * t2["aaaa"]("b,a,m,i") * l2["aaaa"]("L,m,l,b,a")) * r1["aa"]("R,a,i");
        D_oovo_aaaa("L,R,i,j,a,k") += tmps_["perm_aaaa_LLvooo"]("L,R,a,i,j,k");
        D_oovo_aaaa("L,R,i,j,a,k") -= tmps_["perm_aaaa_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_aaaa_LLvooo"].~TArrayD();

        // D_oovo_bbbb += +0.50 P(i,j) r1_bb(a,i) t2_abab(c,b,l,j) l2_abab(l,k,c,b)
        //             += +0.50 P(i,j) r1_bb(a,i) t2_abab(b,c,l,j) l2_abab(l,k,b,c)
        //             += +0.50 r2_abab(a,b,i,l) t2_1_abab(d,c,k,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r2_abab(a,b,i,l) t2_1_abab(c,d,k,j) l2_1_abab(k,l,c,d)
        // flops: o3v1L2 += o3v2L1 o3v2L1 o3v0L1 o4v1L2
        //  mems: o3v1L2 += o2v0L1 o2v0L1 o3v0L1 o4v1L2
        tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k")  = (t2["abab"]("c,b,l,j") * l2["abab"]("L,l,k,c,b") + t2_1["abab"]("d,c,k,j") * l2_1["abab"]("L,k,l,d,c")) * r1["bb"]("R,a,i");
        D_oovo_bbbb("L,R,i,j,a,k") += tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k");
        D_oovo_bbbb("L,R,i,j,a,k") -= tmps_["perm_bbbb_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_bbbb_LLvooo"].~TArrayD();

        // D_oovo_bbbb += +0.50 P(i,j) d_bb(j,k) r0 t1_bb(a,m) t2_abab(c,b,l,i) l2_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t1_bb(a,m) t2_abab(b,c,l,i) l2_abab(l,m,b,c)
        //             += +0.50 d_aa(i,k) r0 t2_1_abab(b,a,m,j) l2_1_abab(m,l,b,a)
        //             += +0.50 d_aa(i,k) r0 t2_1_abab(a,b,m,j) l2_1_abab(m,l,a,b)
        // flops: o3v1L2 += o3v2L1 o3v2L1 o4v0L1 o5v0L1 o4v1L1 o3v1L2
        //  mems: o3v1L2 += o2v0L1 o2v0L1 o4v0L1 o4v0L1 o3v1L1 o3v1L2
        tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k")  = (t2["abab"]("c,b,l,i") * l2["abab"]("L,l,m,c,b") + t2_1["abab"]("b,a,m,j") * l2_1["abab"]("L,m,l,b,a")) * Id["bb_oo"]("j,k") * t1["bb"]("a,m") * r0("R");
        D_oovo_bbbb("L,R,i,j,a,k") += tmps_["perm_bbbb_LLvooo"]("L,R,a,i,j,k");
        D_oovo_bbbb("L,R,i,j,a,k") -= tmps_["perm_bbbb_LLvooo"]("L,R,a,j,i,k");
        tmps_["perm_bbbb_LLvooo"].~TArrayD();

        // D_oovv_aaaa += +1.00 r2_1_aaaa(b,a,i,j) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") += l0_1("L") * r2_1["aaaa"]("R,b,a,i,j");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r1_1_aa(d,i) t2_aaaa(a,c,k,j) t1_aa(b,l) l2_1_aaaa(k,l,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_1_aa(d,i) t2_abab(a,c,j,k) t1_aa(b,l) l2_1_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r1_aa(d,i) t1_aa(b,l) t2_1_aaaa(a,c,k,j) l2_1_aaaa(k,l,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_aa(d,i) t2_abab(a,c,j,k) t1_aa(b,l) l2_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r1_aa(d,i) t1_aa(b,l) t1_aa(c,j) t1_1_aa(a,k) l2_1_aaaa(k,l,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_aa(d,i) t1_aa(b,l) t2_1_abab(a,c,j,k) l2_1_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r1_aa(d,i) t2_aaaa(a,c,k,j) t1_aa(b,l) l2_aaaa(k,l,c,d)
        // flops: o2v2L2 += o3v3L1 o3v3L1 o2v2L1 o3v2L2 o3v3L1 o3v2L2 o3v1L2 o3v3L1 o3v2L2 o3v1L2 o3v2L1 o4v1L2 o4v1L2 o3v1L2 o3v3L1 o3v2L2 o3v1L2 o3v3L1 o3v2L2 o3v1L2 o3v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L1 o2v2L1 o3v1L2 o2v2L1 o3v1L2 o3v1L2 o2v2L1 o3v1L2 o3v1L2 o3v1L1 o4v0L2 o3v1L2 o3v1L2 o2v2L1 o3v1L2 o3v1L2 o2v2L1 o3v1L2 o3v1L2 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = ((t2["aaaa"]("a,c,k,j") * l2_1["aaaa"]("L,k,l,c,d") + -1.00 * t2["abab"]("a,c,j,k") * l2_1["abab"]("L,l,k,d,c")) * r1_1["aa"]("R,d,i") + t2_1["aaaa"]("a,c,k,j") * l2_1["aaaa"]("L,k,l,c,d") * r1["aa"]("R,d,i") + -1.00 * t2["abab"]("a,c,j,k") * l2["abab"]("L,l,k,d,c") * r1["aa"]("R,d,i") + t1["aa"]("c,j") * l2_1["aaaa"]("L,k,l,c,d") * r1["aa"]("R,d,i") * t1_1["aa"]("a,k") + -1.00 * t2_1["abab"]("a,c,j,k") * l2_1["abab"]("L,l,k,d,c") * r1["aa"]("R,d,i") + t2["aaaa"]("a,c,k,j") * l2["aaaa"]("L,k,l,c,d") * r1["aa"]("R,d,i")) * t1["aa"]("b,l");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +0.50 P(i,j) r2_aaaa(b,a,l,i) t2_1_abab(d,c,j,k) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) r2_aaaa(b,a,l,i) t2_1_abab(c,d,j,k) l2_1_abab(l,k,c,d)
        //             += +0.50 d_bb(i,k) r1_aa(a,m) t2_abab(c,b,j,l) l2_abab(m,l,c,b)
        //             += +0.50 d_bb(i,k) r1_aa(a,m) t2_abab(b,c,j,l) l2_abab(m,l,b,c)
        // flops: o2v2L2 += o3v2L1 o3v2L1 o3v0L1 o4v2L2
        //  mems: o2v2L2 += o2v0L1 o2v0L1 o3v0L1 o3v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = (t2_1["abab"]("d,c,j,k") * l2_1["abab"]("L,l,k,d,c") + t2["abab"]("c,b,j,l") * l2["abab"]("L,m,l,c,b")) * r2["aaaa"]("R,b,a,l,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += -0.50 P(i,j) r1_1_aa(d,i) t2_aaaa(b,a,k,l) t1_aa(c,j) l2_1_aaaa(k,l,c,d)
        //             += -0.50 P(i,j) r1_aa(d,i) t1_aa(c,j) t2_1_aaaa(b,a,k,l) l2_1_aaaa(k,l,c,d)
        //             += -0.50 P(i,j) r1_aa(d,i) t2_aaaa(b,a,k,l) t1_aa(c,j) l2_aaaa(k,l,c,d)
        // flops: o2v2L2 += o3v3L1 o3v3L1 o3v3L1 o3v4L2 o2v4L1 o1v4L2 o1v3L2 o2v3L2
        //  mems: o2v2L2 += o3v3L1 o3v3L1 o3v3L1 o1v3L2 o0v4L1 o1v3L2 o1v3L2 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * ((t2["aaaa"]("b,a,k,l") * r1_1["aa"]("R,d,i") + t2_1["aaaa"]("b,a,k,l") * r1["aa"]("R,d,i")) * l2_1["aaaa"]("L,k,l,c,d") + l2["aaaa"]("L,k,l,c,d") * t2["aaaa"]("b,a,k,l") * r1["aa"]("R,d,i")) * t1["aa"]("c,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += -0.50 P(i,j) r0 t2_aaaa(b,a,l,j) t2_1_abab(d,c,i,k) l2_1_abab(l,k,d,c)
        //             += -0.50 P(i,j) r0 t2_aaaa(b,a,l,j) t2_1_abab(c,d,i,k) l2_1_abab(l,k,c,d)
        //             += +0.50 r0 t2_abab(a,b,l,j) t2_abab(d,c,i,k) l2_abab(l,k,d,c)
        //             += +0.50 r0 t2_abab(a,b,l,j) t2_abab(c,d,i,k) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L1 o3v2L1 o2v0L1 o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v0L1 o2v0L1 o2v0L1 o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = (t2_1["abab"]("d,c,i,k") * l2_1["abab"]("L,l,k,d,c") + -1.00 * t2["abab"]("d,c,i,k") * l2["abab"]("L,l,k,d,c")) * t2["aaaa"]("b,a,l,j") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += -0.50 P(i,j) r0 t2_aaaa(b,a,l,j) t2_1_aaaa(d,c,k,i) l2_1_aaaa(k,l,d,c)
        //             += +0.50 r0 t2_abab(a,b,l,j) t2_aaaa(d,c,k,i) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o3v2L1 o2v0L1 o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v0L1 o2v0L1 o2v0L1 o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * (t2_1["aaaa"]("d,c,k,i") * l2_1["aaaa"]("L,k,l,d,c") + -1.00 * t2["aaaa"]("d,c,k,i") * l2["aaaa"]("L,k,l,d,c")) * t2["aaaa"]("b,a,l,j") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +0.50 P(i,j) P(a,b) r1_aa(a,l) t1_aa(b,j) t2_1_abab(d,c,i,k) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,l) t1_aa(b,j) t2_1_abab(c,d,i,k) l2_1_abab(l,k,c,d)
        //             += -0.50 d_bb(i,k) r0 t1_aa(a,m) t2_abab(c,b,j,l) l2_abab(m,l,c,b)
        //             += -0.50 d_bb(i,k) r0 t1_aa(a,m) t2_abab(b,c,j,l) l2_abab(m,l,b,c)
        // flops: o2v2L2 += o3v2L1 o3v2L1 o4v0L1 o4v1L1 o3v2L2
        //  mems: o2v2L2 += o2v0L1 o2v0L1 o4v0L1 o3v1L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = (t2_1["abab"]("d,c,i,k") * l2_1["abab"]("L,l,k,d,c") + -1.00 * t2["abab"]("c,b,j,l") * l2["abab"]("L,m,l,c,b")) * t1["aa"]("b,j") * r1["aa"]("R,a,l");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(a,b) r2_aaaa(a,c,i,j) t1_1_aa(b,k) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = l1_1["aa"]("L,k,c") * r2["aaaa"]("R,a,c,i,j") * t1_1["aa"]("b,k");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_abab(a,d,i,l) t2_1_abab(b,c,j,k) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v3L1 o3v3L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = l2_1["bbbb"]("L,l,k,c,d") * t2_1["abab"]("b,c,j,k") * r2["abab"]("R,a,d,i,l");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r2_1_abab(a,b,i,j) l0_1
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= l0_1("L") * r2_1["abab"]("R,a,b,i,j");

        // D_oovv_abab += +0.50 r2_abab(a,b,i,l) t2_abab(d,c,k,j) l2_abab(k,l,d,c)
        //             += +0.50 r2_abab(a,b,i,l) t2_abab(c,d,k,j) l2_abab(k,l,c,d)
        //             += +0.50 P(i,j) r1_bb(a,i) t2_1_abab(c,b,l,j) l2_1_abab(l,k,c,b)
        //             += +0.50 P(i,j) r1_bb(a,i) t2_1_abab(b,c,l,j) l2_1_abab(l,k,b,c)
        // flops: o2v2L2 += o3v2L1 o3v2L1 o3v0L1 o4v2L2
        //  mems: o2v2L2 += o2v0L1 o2v0L1 o3v0L1 o3v2L2
        D_oovv_abab("L,R,i,j,a,b") += (t2["abab"]("d,c,k,j") * l2["abab"]("L,k,l,d,c") + t2_1["abab"]("c,b,l,j") * l2_1["abab"]("L,l,k,c,b")) * r2["abab"]("R,a,b,i,l");

        // D_oovv_abab += +0.50 r2_abab(a,b,l,j) t2_1_abab(d,c,i,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r2_abab(a,b,l,j) t2_1_abab(c,d,i,k) l2_1_abab(l,k,c,d)
        //             += +0.50 P(i,j) r1_aa(a,i) t2_abab(c,b,j,l) l2_abab(k,l,c,b)
        //             += +0.50 P(i,j) r1_aa(a,i) t2_abab(b,c,j,l) l2_abab(k,l,b,c)
        // flops: o2v2L2 += o3v2L1 o3v2L1 o4v0L1 o5v2L2
        //  mems: o2v2L2 += o2v0L1 o2v0L1 o4v0L1 o4v2L2
        D_oovv_abab("L,R,i,j,a,b") += (t2_1["abab"]("d,c,i,k") * l2_1["abab"]("L,l,k,d,c") + t2["abab"]("c,b,j,l") * l2["abab"]("L,k,l,c,b")) * r2["abab"]("R,a,b,l,j");

        // D_oovv_abab += +0.50 r0 t2_abab(a,b,l,j) t2_1_abab(d,c,i,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r0 t2_abab(a,b,l,j) t2_1_abab(c,d,i,k) l2_1_abab(l,k,c,d)
        //             += -0.50 P(i,j) r0 t2_aaaa(b,a,l,j) t2_abab(d,c,i,k) l2_abab(l,k,d,c)
        //             += -0.50 P(i,j) r0 t2_aaaa(b,a,l,j) t2_abab(c,d,i,k) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L1 o3v2L1 o2v0L1 o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v0L1 o2v0L1 o2v0L1 o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += (t2_1["abab"]("d,c,i,k") * l2_1["abab"]("L,l,k,d,c") + -1.00 * t2["abab"]("d,c,i,k") * l2["abab"]("L,l,k,d,c")) * t2["abab"]("a,b,l,j") * r0("R");

        // D_oovv_abab += +0.50 r0 t2_abab(a,b,l,j) t2_1_aaaa(d,c,k,i) l2_1_aaaa(k,l,d,c)
        //             += -0.50 P(i,j) r0 t2_aaaa(b,a,l,j) t2_aaaa(d,c,k,i) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o3v2L1 o2v0L1 o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v0L1 o2v0L1 o2v0L1 o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * (t2_1["aaaa"]("d,c,k,i") * l2_1["aaaa"]("L,k,l,d,c") + -1.00 * t2["aaaa"]("d,c,k,i") * l2["aaaa"]("L,k,l,d,c")) * t2["abab"]("a,b,l,j") * r0("R");

        // D_oovv_abab += +1.00 r2_abab(c,b,i,j) t1_1_aa(a,k) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += l1_1["aa"]("L,k,c") * r2["abab"]("R,c,b,i,j") * t1_1["aa"]("a,k");

        // D_oovv_abab += +1.00 r2_abab(c,b,i,j) t1_aa(a,k) l1_aa(k,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += l1["aa"]("L,k,c") * r2["abab"]("R,c,b,i,j") * t1["aa"]("a,k");

        // D_oovv_abab += +1.00 r2_1_abab(c,b,i,j) t1_aa(a,k) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += l1_1["aa"]("L,k,c") * r2_1["abab"]("R,c,b,i,j") * t1["aa"]("a,k");

        // D_oovv_abab += +0.50 r0 t2_abab(a,b,i,l) t2_abab(d,c,k,j) l2_abab(k,l,d,c)
        //             += +0.50 r0 t2_abab(a,b,i,l) t2_abab(c,d,k,j) l2_abab(k,l,c,d)
        //             += -0.50 P(i,j) r0 t2_bbbb(b,a,l,j) t2_1_abab(d,c,k,i) l2_1_abab(k,l,d,c)
        //             += -0.50 P(i,j) r0 t2_bbbb(b,a,l,j) t2_1_abab(c,d,k,i) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L1 o3v2L1 o3v0L1 o4v2L1 o3v2L2
        //  mems: o2v2L2 += o2v0L1 o2v0L1 o3v0L1 o3v2L1 o3v2L2
        D_oovv_abab("L,R,i,j,a,b") += (t2["abab"]("d,c,k,j") * l2["abab"]("L,k,l,d,c") + -1.00 * t2_1["abab"]("d,c,k,i") * l2_1["abab"]("L,k,l,d,c")) * t2["abab"]("a,b,i,l") * r0("R");

        // D_oovv_abab += +1.00 r2_abab(a,c,i,j) t1_1_bb(b,k) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += l1_1["bb"]("L,k,c") * r2["abab"]("R,a,c,i,j") * t1_1["bb"]("b,k");

        // D_oovv_abab += +0.50 r1_bb(b,l) t1_aa(a,i) t2_abab(d,c,k,j) l2_abab(k,l,d,c)
        //             += +0.50 r1_bb(b,l) t1_aa(a,i) t2_abab(c,d,k,j) l2_abab(k,l,c,d)
        //             += +0.50 P(i,j) r0 t1_bb(a,j) t2_1_abab(c,b,l,i) l2_1_abab(l,k,c,b)
        //             += +0.50 P(i,j) r0 t1_bb(a,j) t2_1_abab(b,c,l,i) l2_1_abab(l,k,b,c)
        // flops: o2v2L2 += o3v2L1 o3v2L1 o4v0L1 o4v1L2 o4v2L2
        //  mems: o2v2L2 += o2v0L1 o2v0L1 o4v0L1 o3v1L2 o4v2L2
        D_oovv_abab("L,R,i,j,a,b") += (t2["abab"]("d,c,k,j") * l2["abab"]("L,k,l,d,c") + t2_1["abab"]("c,b,l,i") * l2_1["abab"]("L,l,k,c,b")) * r1["bb"]("R,b,l") * t1["aa"]("a,i");

        // D_oovv_abab += -1.00 r2_abab(a,d,i,l) t2_1_bbbb(b,c,k,j) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v3L1 o3v3L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= l2_1["bbbb"]("L,l,k,c,d") * t2_1["bbbb"]("b,c,k,j") * r2["abab"]("R,a,d,i,l");

        // D_oovv_abab += -1.00 r2_abab(d,b,l,j) t2_aaaa(a,c,k,i) l2_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v3L1 o3v3L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= l2["aaaa"]("L,l,k,c,d") * t2["aaaa"]("a,c,k,i") * r2["abab"]("R,d,b,l,j");

        // D_oovv_abab += +1.00 r1_aa(d,i) t1_aa(a,l) t2_abab(c,b,k,j) l2_aaaa(k,l,c,d)
        // flops: o2v2L2 += o3v3L1 o3v2L2 o3v2L2
        //  mems: o2v2L2 += o2v2L1 o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += l2["aaaa"]("L,k,l,c,d") * t2["abab"]("c,b,k,j") * r1["aa"]("R,d,i") * t1["aa"]("a,l");

        // D_oovv_abab += +1.00 r1_aa(d,i) t1_aa(a,l) t2_1_abab(c,b,k,j) l2_1_aaaa(k,l,c,d)
        // flops: o2v2L2 += o3v3L1 o3v2L2 o3v2L2
        //  mems: o2v2L2 += o2v2L1 o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += l2_1["aaaa"]("L,k,l,c,d") * t2_1["abab"]("c,b,k,j") * r1["aa"]("R,d,i") * t1["aa"]("a,l");

        // D_oovv_bbbb += +1.00 r2_bbbb(b,a,i,j) l0
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += l0("L") * r2["bbbb"]("R,b,a,i,j");

        // D_oovv_bbbb += -1.00 P(i,j) P(a,b) r1_bb(d,i) t2_abab(c,a,k,j) t1_bb(b,l) l2_abab(k,l,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_bb(d,i) t1_bb(b,l) t1_bb(c,j) t1_1_bb(a,k) l2_1_bbbb(k,l,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_bb(d,i) t1_bb(b,l) t2_1_abab(c,a,k,j) l2_1_abab(k,l,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_bb(d,i) t1_bb(b,l) t2_1_bbbb(a,c,k,j) l2_1_bbbb(k,l,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_bb(d,i) t2_bbbb(a,c,k,j) t1_bb(b,l) l2_bbbb(k,l,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_1_bb(d,i) t2_abab(c,a,k,j) t1_bb(b,l) l2_1_abab(k,l,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_1_bb(d,i) t2_bbbb(a,c,k,j) t1_bb(b,l) l2_1_bbbb(k,l,c,d)
        // flops: o2v2L2 += o3v3L1 o3v3L1 o3v3L1 o3v2L1 o3v2L1 o3v3L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v2L2 o3v3L1 o3v2L2 o3v1L2 o3v3L1 o3v2L2 o3v1L2 o3v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L1 o2v2L1 o3v1L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o2v2L1 o3v1L2 o2v2L1 o3v1L2 o3v1L2 o2v2L1 o3v1L2 o3v1L2 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = ((-1.00 * t2["bbbb"]("a,c,k,j") * l2["bbbb"]("L,k,l,c,d") + -1.00 * l2_1["bbbb"]("L,k,l,c,d") * t2_1["bbbb"]("a,c,k,j") + t2_1["abab"]("c,a,k,j") * l2_1["abab"]("L,k,l,c,d") + -1.00 * l2_1["bbbb"]("L,k,l,c,d") * t1["bb"]("c,j") * t1_1["bb"]("a,k") + l2["abab"]("L,k,l,c,d") * t2["abab"]("c,a,k,j")) * r1["bb"]("R,d,i") + t2["abab"]("c,a,k,j") * l2_1["abab"]("L,k,l,c,d") * r1_1["bb"]("R,d,i") + -1.00 * t2["bbbb"]("a,c,k,j") * l2_1["bbbb"]("L,k,l,c,d") * r1_1["bb"]("R,d,i")) * t1["bb"]("b,l");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(a,b) r1_1_aa(d,l) t2_bbbb(a,c,i,j) t1_bb(b,k) l2_1_abab(l,k,d,c)
        //             += -1.00 P(a,b) r1_1_bb(d,l) t2_bbbb(a,c,i,j) t1_bb(b,k) l2_1_bbbb(l,k,c,d)
        //             += -1.00 P(a,b) r1_bb(d,l) t2_bbbb(a,c,i,j) t1_bb(b,k) l2_bbbb(l,k,c,d)
        //             += +1.00 P(a,b) r1_1_bb(a,l) t1_bb(b,k) t1_bb(c,i) t1_bb(d,j) l2_1_bbbb(l,k,d,c)
        //             += +1.00 P(a,b) r1_bb(a,l) t1_bb(b,k) t1_bb(c,i) t1_bb(d,j) l2_bbbb(l,k,d,c)
        //             += -0.50 P(a,b) r1_1_bb(a,l) t1_bb(b,k) t2_bbbb(d,c,i,j) l2_1_bbbb(l,k,d,c)
        //             += +1.00 P(a,b) r1_aa(d,l) t1_bb(b,k) t2_1_bbbb(a,c,i,j) l2_1_abab(l,k,d,c)
        //             += +1.00 P(a,b) r1_aa(d,l) t2_bbbb(a,c,i,j) t1_bb(b,k) l2_abab(l,k,d,c)
        //             += +1.00 P(a,b) r2_bbbb(a,c,i,j) t1_bb(b,k) l1_bb(k,c)
        //             += +1.00 P(a,b) r2_1_bbbb(a,c,i,j) t1_bb(b,k) l1_1_bb(k,c)
        //             += -1.00 P(a,b) r1_bb(d,l) t1_bb(b,k) t2_1_bbbb(a,c,i,j) l2_1_bbbb(l,k,c,d)
        // flops: o2v2L2 += o2v2L2 o2v2L2 o1v1L2 o3v2L2 o2v2L2 o3v2L2 o3v2L1 o4v1L1 o4v1L2 o3v1L2 o3v2L1 o4v1L1 o4v1L2 o3v1L2 o4v2L1 o4v1L2 o3v1L2 o2v2L2 o3v2L2 o3v1L2 o2v2L2 o3v2L2 o3v1L2 o3v2L2 o3v1L2 o3v2L2 o3v1L2 o2v2L2 o3v2L2 o3v1L2 o3v1L2 o3v2L2
        //  mems: o2v2L2 += o1v1L2 o1v1L2 o1v1L2 o3v1L2 o1v1L2 o3v1L2 o3v1L1 o4v0L1 o3v1L2 o3v1L2 o3v1L1 o4v0L1 o3v1L2 o3v1L2 o4v0L1 o3v1L2 o3v1L2 o1v1L2 o3v1L2 o3v1L2 o1v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o1v1L2 o3v1L2 o3v1L2 o3v1L2 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = (t2["bbbb"]("a,c,i,j") * (l2_1["abab"]("L,l,k,d,c") * r1_1["aa"]("R,d,l") + -1.00 * l2_1["bbbb"]("L,l,k,c,d") * r1_1["bb"]("R,d,l")) + -1.00 * (l2["bbbb"]("L,l,k,c,d") * r1["bb"]("R,d,l") * t2["bbbb"]("a,c,i,j") + -1.00 * t1["bb"]("c,i") * l2_1["bbbb"]("L,l,k,d,c") * t1["bb"]("d,j") * r1_1["bb"]("R,a,l") + -1.00 * t1["bb"]("c,i") * l2["bbbb"]("L,l,k,d,c") * t1["bb"]("d,j") * r1["bb"]("R,a,l") + 0.50 * t2["bbbb"]("d,c,i,j") * l2_1["bbbb"]("L,l,k,d,c") * r1_1["bb"]("R,a,l") + -1.00 * l2_1["abab"]("L,l,k,d,c") * r1["aa"]("R,d,l") * t2_1["bbbb"]("a,c,i,j") + -1.00 * l2["abab"]("L,l,k,d,c") * r1["aa"]("R,d,l") * t2["bbbb"]("a,c,i,j") + -1.00 * l1["bb"]("L,k,c") * r2["bbbb"]("R,a,c,i,j") + -1.00 * l1_1["bb"]("L,k,c") * r2_1["bbbb"]("R,a,c,i,j") + l2_1["bbbb"]("L,l,k,c,d") * r1["bb"]("R,d,l") * t2_1["bbbb"]("a,c,i,j"))) * t1["bb"]("b,k");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +0.50 P(i,j) r2_bbbb(b,a,l,i) t2_abab(d,c,k,j) l2_abab(k,l,d,c)
        //             += +0.50 P(i,j) r2_bbbb(b,a,l,i) t2_abab(c,d,k,j) l2_abab(k,l,c,d)
        //             += +0.50 r0 t1_aa(a,j) t2_1_abab(c,b,l,i) l2_1_abab(l,k,c,b)
        //             += +0.50 r0 t1_aa(a,j) t2_1_abab(b,c,l,i) l2_1_abab(l,k,b,c)
        // flops: o2v2L2 += o3v2L1 o3v2L1 o4v0L1 o4v2L2
        //  mems: o2v2L2 += o2v0L1 o2v0L1 o4v0L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = (t2["abab"]("d,c,k,j") * l2["abab"]("L,k,l,d,c") + t2_1["abab"]("c,b,l,i") * l2_1["abab"]("L,l,k,c,b")) * r2["bbbb"]("R,b,a,l,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += -0.50 P(i,j) r1_1_bb(d,i) t2_bbbb(b,a,k,l) t1_bb(c,j) l2_1_bbbb(k,l,c,d)
        //             += -0.50 P(i,j) r1_bb(d,i) t1_bb(c,j) t2_1_bbbb(b,a,k,l) l2_1_bbbb(k,l,c,d)
        //             += -0.50 P(i,j) r1_bb(d,i) t2_bbbb(b,a,k,l) t1_bb(c,j) l2_bbbb(k,l,c,d)
        // flops: o2v2L2 += o3v3L1 o3v3L1 o3v3L1 o3v4L2 o2v4L1 o1v4L2 o1v3L2 o2v3L2
        //  mems: o2v2L2 += o3v3L1 o3v3L1 o3v3L1 o1v3L2 o0v4L1 o1v3L2 o1v3L2 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * ((t2["bbbb"]("b,a,k,l") * r1_1["bb"]("R,d,i") + t2_1["bbbb"]("b,a,k,l") * r1["bb"]("R,d,i")) * l2_1["bbbb"]("L,k,l,c,d") + l2["bbbb"]("L,k,l,c,d") * t2["bbbb"]("b,a,k,l") * r1["bb"]("R,d,i")) * t1["bb"]("c,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += -0.50 P(i,j) r0 t2_bbbb(b,a,l,j) t2_abab(d,c,k,i) l2_abab(k,l,d,c)
        //             += -0.50 P(i,j) r0 t2_bbbb(b,a,l,j) t2_abab(c,d,k,i) l2_abab(k,l,c,d)
        //             += +0.50 r0 t2_abab(a,b,i,l) t2_1_abab(d,c,k,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r0 t2_abab(a,b,i,l) t2_1_abab(c,d,k,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L1 o3v2L1 o3v0L1 o3v2L1 o1v2L2
        //  mems: o2v2L2 += o2v0L1 o2v0L1 o3v0L1 o1v2L1 o1v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = (t2["abab"]("d,c,k,i") * l2["abab"]("L,k,l,d,c") + -1.00 * t2_1["abab"]("d,c,k,j") * l2_1["abab"]("L,k,l,d,c")) * t2["bbbb"]("b,a,l,j") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += -0.50 P(i,j) r0 t2_bbbb(b,a,l,j) t2_bbbb(d,c,k,i) l2_bbbb(k,l,d,c)
        //             += -0.50 r1_aa(a,j) t2_1_bbbb(c,b,l,i) l2_1_bbbb(l,k,c,b)
        // flops: o2v2L2 += o3v2L1 o3v2L1 o3v0L1 o4v2L1 o3v2L2
        //  mems: o2v2L2 += o2v0L1 o2v0L1 o3v0L1 o3v2L1 o3v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * (t2["bbbb"]("d,c,k,i") * l2["bbbb"]("L,k,l,d,c") + t2_1["bbbb"]("c,b,l,i") * l2_1["bbbb"]("L,l,k,c,b")) * t2["bbbb"]("b,a,l,j") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(a,b) r2_bbbb(a,c,i,j) t1_1_bb(b,k) l1_1_bb(k,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = l1_1["bb"]("L,k,c") * r2["bbbb"]("R,a,c,i,j") * t1_1["bb"]("b,k");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r2_abab(d,a,l,i) t2_abab(c,b,k,j) l2_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v3L1 o3v3L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = l2["aaaa"]("L,l,k,c,d") * t2["abab"]("c,b,k,j") * r2["abab"]("R,d,a,l,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_vvvv_abab += -1.00 r1_1_bb(d,j) t1_aa(c,i) l2_1_abab(i,j,a,b)
        // flops: o0v4L2 += o2v3L2 o1v4L2
        //  mems: o0v4L2 += o1v3L2 o0v4L2
        D_vvvv_abab("L,R,a,b,c,d") -= l2_1["abab"]("L,i,j,a,b") * r1_1["bb"]("R,d,j") * t1["aa"]("c,i");

        // D_vvvv_abab += -1.00 r1_bb(d,j) t1_aa(c,i) l2_abab(i,j,a,b)
        // flops: o0v4L2 += o2v3L2 o1v4L2
        //  mems: o0v4L2 += o1v3L2 o0v4L2
        D_vvvv_abab("L,R,a,b,c,d") -= l2["abab"]("L,i,j,a,b") * r1["bb"]("R,d,j") * t1["aa"]("c,i");

        // flops: o1v3L2  = o3v2L1 o3v3L2
        //  mems: o1v3L2  = o3v1L1 o1v3L2
        tmps_["1_aaaa_LLvovv"]("L,R,a,i,b,c")  = l2_1["aaaa"]("L,j,k,a,d") * t1_1["aa"]("d,i") * r2["aaaa"]("R,b,c,j,k");

        // D_ovvv_aaaa  = +0.50 r2_aaaa(c,b,k,j) t1_1_aa(d,i) l2_1_aaaa(k,j,a,d)
        D_ovvv_aaaa("L,R,i,a,b,c")  = 0.50 * tmps_["1_aaaa_LLvovv"]("L,R,a,i,c,b");

        // D_vovv_aaaa  = -0.50 r2_aaaa(c,b,k,j) t1_1_aa(d,i) l2_1_aaaa(k,j,a,d)
        D_vovv_aaaa("L,R,a,i,b,c")  = -0.50 * tmps_["1_aaaa_LLvovv"]("L,R,a,i,c,b");

        // D_oovv_aaaa += -0.50 P(i,j) r2_aaaa(b,a,l,k) t1_aa(d,j) t1_1_aa(c,i) l2_1_aaaa(l,k,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["1_aaaa_LLvovv"]("L,R,d,i,b,a") * t1["aa"]("d,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();
        tmps_["1_aaaa_LLvovv"].~TArrayD();

        // flops: o1v3L2  = o3v2L1 o3v3L2
        //  mems: o1v3L2  = o3v1L1 o1v3L2
        tmps_["2_abab_LLvovv"]("L,R,a,i,b,c")  = l2_1["abab"]("L,j,k,a,d") * t1_1["bb"]("d,i") * r2["abab"]("R,b,c,j,k");

        // D_ovvv_abab  = +0.50 r2_abab(c,b,k,j) t1_1_bb(d,i) l2_1_abab(k,j,a,d)
        //             += +0.50 r2_abab(c,b,j,k) t1_1_bb(d,i) l2_1_abab(j,k,a,d)
        D_ovvv_abab("L,R,i,a,b,c")  = tmps_["2_abab_LLvovv"]("L,R,a,i,c,b");

        // D_vovv_abab  = -0.50 r2_abab(c,b,k,j) t1_1_bb(d,i) l2_1_abab(k,j,a,d)
        //             += -0.50 r2_abab(c,b,j,k) t1_1_bb(d,i) l2_1_abab(j,k,a,d)
        D_vovv_abab("L,R,a,i,b,c")  = -1.00 * tmps_["2_abab_LLvovv"]("L,R,a,i,c,b");

        // D_oovv_abab += -0.50 r2_abab(a,b,l,k) t1_aa(d,i) t1_1_bb(c,j) l2_1_abab(l,k,d,c)
        //             += -0.50 r2_abab(a,b,k,l) t1_aa(d,i) t1_1_bb(c,j) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["aa"]("d,i") * tmps_["2_abab_LLvovv"]("L,R,d,j,a,b");
        tmps_["2_abab_LLvovv"].~TArrayD();

        // flops: o1v3L2  = o3v2L1 o3v3L2
        //  mems: o1v3L2  = o3v1L1 o1v3L2
        tmps_["3_bbbb_LLvovv"]("L,R,a,i,b,c")  = l2_1["bbbb"]("L,j,k,a,d") * t1_1["bb"]("d,i") * r2["bbbb"]("R,b,c,j,k");

        // D_ovvv_bbbb  = +0.50 r2_bbbb(c,b,k,j) t1_1_bb(d,i) l2_1_bbbb(k,j,a,d)
        D_ovvv_bbbb("L,R,i,a,b,c")  = 0.50 * tmps_["3_bbbb_LLvovv"]("L,R,a,i,c,b");

        // D_vovv_bbbb  = -0.50 r2_bbbb(c,b,k,j) t1_1_bb(d,i) l2_1_bbbb(k,j,a,d)
        D_vovv_bbbb("L,R,a,i,b,c")  = -0.50 * tmps_["3_bbbb_LLvovv"]("L,R,a,i,c,b");

        // D_oovv_bbbb += -0.50 P(i,j) r2_bbbb(b,a,l,k) t1_bb(d,j) t1_1_bb(c,i) l2_1_bbbb(l,k,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["3_bbbb_LLvovv"]("L,R,d,i,b,a") * t1["bb"]("d,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
        tmps_["3_bbbb_LLvovv"].~TArrayD();

        // flops: o2v2L2  = o3v3L2
        //  mems: o2v2L2  = o2v2L2
        tmps_["4_abba_LLovvo"]("L,R,i,a,b,j")  = l2_1["abab"]("L,i,k,c,a") * r2["abab"]("R,c,b,j,k");

        // D_oovv_abab += -1.00 r2_abab(d,b,i,l) t2_1_abab(a,c,k,j) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2_1["abab"]("a,c,k,j") * tmps_["4_abba_LLovvo"]("L,R,k,c,b,i");

        // D_oovv_abab += -1.00 r2_abab(d,b,i,l) t1_aa(a,k) t1_1_bb(c,j) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["4_abba_LLovvo"]("L,R,k,c,b,i") * t1_1["bb"]("c,j") * t1["aa"]("a,k");

        // D_oovv_abab += -1.00 r2_abab(d,b,i,l) t1_bb(c,j) t1_1_aa(a,k) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["4_abba_LLovvo"]("L,R,k,c,b,i") * t1["bb"]("c,j") * t1_1["aa"]("a,k");
        tmps_["4_abba_LLovvo"].~TArrayD();

        // flops: o3v1L1  = o3v2L1 o4v2L1
        //  mems: o3v1L1  = o3v1L1 o3v1L1
        tmps_["5_bbaa_Loovo"]("L,i,j,a,k")  = t1["bb"]("c,i") * l2_1["bbbb"]("L,j,l,c,b") * t2["abab"]("a,b,k,l");

        // D_oovv_abab += +1.00 r1_1_bb(b,l) t2_abab(a,c,i,k) t1_bb(d,j) l2_1_bbbb(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["5_bbaa_Loovo"]("L,j,l,a,i") * r1_1["bb"]("R,b,l");

        // D_oovv_abab += +1.00 r0 t2_abab(a,c,i,l) t1_bb(d,j) t1_1_bb(b,k) l2_1_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1_1["bb"]("b,k") * tmps_["5_bbaa_Loovo"]("L,j,k,a,i") * r0("R");
        tmps_["5_bbaa_Loovo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["6_bbaa_Lvoov"]("L,a,i,j,b")  = t2["abab"]("c,a,k,i") * l2_1["aaaa"]("L,k,j,c,b");

        // D_oovv_abab += +1.00 r1_1_aa(d,i) t1_aa(a,l) t2_abab(c,b,k,j) l2_1_aaaa(k,l,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["6_bbaa_Lvoov"]("L,b,j,l,d") * r1_1["aa"]("R,d,i") * t1["aa"]("a,l");

        // D_oovv_abab += +1.00 r0 t1_aa(a,l) t2_abab(d,b,k,j) t1_1_aa(c,i) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o3v2L1 o2v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["6_bbaa_Lvoov"]("L,b,j,l,c") * t1_1["aa"]("c,i") * t1["aa"]("a,l") * r0("R");
        tmps_["6_bbaa_Lvoov"].~TArrayD();

        // flops: o1v3L1  = o2v3L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["7_baab_Lvovv"]("L,a,i,b,c")  = t1["bb"]("a,j") * l2_1["abab"]("L,i,j,b,c");

        // D_vvvv_abab += -1.00 r1_1_aa(c,j) t1_bb(d,i) l2_1_abab(j,i,a,b)
        // flops: o0v4L2 += o1v4L2
        //  mems: o0v4L2 += o0v4L2
        D_vvvv_abab("L,R,a,b,c,d") -= tmps_["7_baab_Lvovv"]("L,d,j,a,b") * r1_1["aa"]("R,c,j");

        // D_vvvv_abab += -1.00 r0_1 t1_aa(c,i) t1_bb(d,j) l2_1_abab(i,j,a,b)
        // flops: o0v4L2 += o1v4L1 o0v4L2
        //  mems: o0v4L2 += o0v4L1 o0v4L2
        D_vvvv_abab("L,R,a,b,c,d") -= t1["aa"]("c,i") * tmps_["7_baab_Lvovv"]("L,d,i,a,b") * r0_1("R");

        // D_vvvv_abab += -1.00 r0 t1_bb(d,j) t1_1_aa(c,i) l2_1_abab(i,j,a,b)
        // flops: o0v4L2 += o1v4L1 o0v4L2
        //  mems: o0v4L2 += o0v4L1 o0v4L2
        D_vvvv_abab("L,R,a,b,c,d") -= t1_1["aa"]("c,i") * tmps_["7_baab_Lvovv"]("L,d,i,a,b") * r0("R");
        tmps_["7_baab_Lvovv"].~TArrayD();

        // flops: o1v3L1  = o2v3L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["8_baab_Lvovv"]("L,a,i,b,c")  = t1_1["bb"]("a,j") * l2_1["abab"]("L,i,j,b,c");

        // D_vvvv_abab += -1.00 r1_aa(c,j) t1_1_bb(d,i) l2_1_abab(j,i,a,b)
        // flops: o0v4L2 += o1v4L2
        //  mems: o0v4L2 += o0v4L2
        D_vvvv_abab("L,R,a,b,c,d") -= tmps_["8_baab_Lvovv"]("L,d,j,a,b") * r1["aa"]("R,c,j");

        // D_vvvv_abab += -1.00 r0 t1_aa(c,j) t1_1_bb(d,i) l2_1_abab(j,i,a,b)
        // flops: o0v4L2 += o1v4L1 o0v4L2
        //  mems: o0v4L2 += o0v4L1 o0v4L2
        D_vvvv_abab("L,R,a,b,c,d") -= t1["aa"]("c,j") * tmps_["8_baab_Lvovv"]("L,d,j,a,b") * r0("R");
        tmps_["8_baab_Lvovv"].~TArrayD();

        // flops: o1v3L1  = o2v3L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["9_baab_Lvovv"]("L,a,i,b,c")  = t1["bb"]("a,j") * l2["abab"]("L,i,j,b,c");

        // D_vvvv_abab += -1.00 r1_aa(c,j) t1_bb(d,i) l2_abab(j,i,a,b)
        // flops: o0v4L2 += o1v4L2
        //  mems: o0v4L2 += o0v4L2
        D_vvvv_abab("L,R,a,b,c,d") -= tmps_["9_baab_Lvovv"]("L,d,j,a,b") * r1["aa"]("R,c,j");

        // D_vvvv_abab += -1.00 r0 t1_aa(c,i) t1_bb(d,j) l2_abab(i,j,a,b)
        // flops: o0v4L2 += o1v4L1 o0v4L2
        //  mems: o0v4L2 += o0v4L1 o0v4L2
        D_vvvv_abab("L,R,a,b,c,d") -= t1["aa"]("c,i") * tmps_["9_baab_Lvovv"]("L,d,i,a,b") * r0("R");
        tmps_["9_baab_Lvovv"].~TArrayD();

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["10_baba_Lvooo"]("L,a,i,j,k")  = t2["abab"]("b,a,i,j") * l1_1["aa"]("L,k,b");

        // D_oovv_abab += +1.00 r1_1_aa(a,k) t2_abab(c,b,i,j) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["10_baba_Lvooo"]("L,b,i,j,k") * r1_1["aa"]("R,a,k");

        // D_oovv_abab += +1.00 r0 t2_abab(c,b,i,j) t1_1_aa(a,k) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1_1["aa"]("a,k") * tmps_["10_baba_Lvooo"]("L,b,i,j,k") * r0("R");

        // D_oovv_abab += +1.00 r0_1 t1_aa(a,k) t2_abab(c,b,i,j) l1_1_aa(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1["aa"]("a,k") * tmps_["10_baba_Lvooo"]("L,b,i,j,k") * r0_1("R");
        tmps_["10_baba_Lvooo"].~TArrayD();

        // flops: o2v2L2  = o3v3L2 o3v3L2 o2v2L2
        //  mems: o2v2L2  = o2v2L2 o2v2L2 o2v2L2
        tmps_["11_abba_LLovvo"]("L,R,i,a,b,j")  = l2["abab"]("L,i,k,c,a") * r2["abab"]("R,c,b,j,k");
        tmps_["11_abba_LLovvo"]("L,R,i,a,b,j") += l2_1["abab"]("L,i,k,c,a") * r2_1["abab"]("R,c,b,j,k");

        // D_oovv_abab += -1.00 r2_abab(d,b,i,l) t2_abab(a,c,k,j) l2_abab(k,l,d,c)
        //             += -1.00 r2_1_abab(d,b,i,l) t2_abab(a,c,k,j) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,c,k,j") * tmps_["11_abba_LLovvo"]("L,R,k,c,b,i");

        // D_oovv_abab += -1.00 r2_abab(d,b,i,l) t1_aa(a,k) t1_bb(c,j) l2_abab(k,l,d,c)
        //             += -1.00 r2_1_abab(d,b,i,l) t1_aa(a,k) t1_bb(c,j) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["11_abba_LLovvo"]("L,R,k,c,b,i") * t1["bb"]("c,j") * t1["aa"]("a,k");
        tmps_["11_abba_LLovvo"].~TArrayD();

        // flops: o1v3L2  = o2v3L2 o2v3L2 o1v3L2
        //  mems: o1v3L2  = o1v3L2 o1v3L2 o1v3L2
        tmps_["12_bbbb_LLovvv"]("L,R,i,a,b,c")  = l2["bbbb"]("L,j,i,a,b") * r1["bb"]("R,c,j");
        tmps_["12_bbbb_LLovvv"]("L,R,i,a,b,c") += l2_1["bbbb"]("L,j,i,a,b") * r1_1["bb"]("R,c,j");

        // D_vvov_bbbb  = +1.00 r1_bb(c,j) l2_bbbb(j,i,a,b)
        //             += +1.00 r1_1_bb(c,j) l2_1_bbbb(j,i,a,b)
        D_vvov_bbbb("L,R,a,b,i,c")  = tmps_["12_bbbb_LLovvv"]("L,R,i,a,b,c");

        // D_vvvo_bbbb  = -1.00 r1_bb(c,j) l2_bbbb(j,i,a,b)
        //             += -1.00 r1_1_bb(c,j) l2_1_bbbb(j,i,a,b)
        D_vvvo_bbbb("L,R,a,b,c,i")  = -1.00 * tmps_["12_bbbb_LLovvv"]("L,R,i,a,b,c");

        // D_vvvv_bbbb += -1.00 P(c,d) r1_bb(c,j) t1_bb(d,i) l2_bbbb(j,i,a,b)
        //             += -1.00 P(c,d) r1_1_bb(c,j) t1_bb(d,i) l2_1_bbbb(j,i,a,b)
        // flops: o0v4L2 += o1v4L2
        //  mems: o0v4L2 += o0v4L2
        tmps_["perm_bbbb_LLvvvv"]("L,R,a,b,c,d")  = tmps_["12_bbbb_LLovvv"]("L,R,i,a,b,c") * t1["bb"]("d,i");
        D_vvvv_bbbb("L,R,a,b,c,d") -= tmps_["perm_bbbb_LLvvvv"]("L,R,a,b,c,d");
        D_vvvv_bbbb("L,R,a,b,c,d") += tmps_["perm_bbbb_LLvvvv"]("L,R,a,b,d,c");
        tmps_["perm_bbbb_LLvvvv"].~TArrayD();
        tmps_["12_bbbb_LLovvv"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["13_abba_Lvoov"]("L,a,i,j,b")  = t2["abab"]("a,c,k,i") * l2_1["abab"]("L,k,j,b,c");

        // D_oovv_abab += -1.00 r1_aa(d,i) t2_abab(a,c,l,j) t1_1_bb(b,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["13_abba_Lvoov"]("L,a,j,k,d") * r1["aa"]("R,d,i") * t1_1["bb"]("b,k");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["16_aabb_Lvoov"]("L,a,i,j,b")  = t2["aaaa"]("a,c,k,i") * l2_1["abab"]("L,k,j,c,b");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_1_abab(a,d,i,l) t2_aaaa(b,c,k,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["16_aabb_Lvoov"]("L,b,j,l,d") * r2_1["abab"]("R,a,d,i,l");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r0 t2_aaaa(b,d,l,j) t2_1_abab(a,c,i,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2_1["abab"]("a,c,i,k") * tmps_["16_aabb_Lvoov"]("L,b,j,k,c") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r2_1_bbbb(b,d,l,j) t2_aaaa(a,c,k,i) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["16_aabb_Lvoov"]("L,a,i,l,d") * r2_1["bbbb"]("R,b,d,l,j");

        // D_oovv_abab += -1.00 r0 t2_aaaa(a,d,l,i) t2_1_bbbb(b,c,k,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["16_aabb_Lvoov"]("L,a,i,k,c") * t2_1["bbbb"]("b,c,k,j") * r0("R");

        // D_oovv_abab += -1.00 r0_1 t2_aaaa(a,c,k,i) t2_bbbb(b,d,l,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["16_aabb_Lvoov"]("L,a,i,l,d") * t2["bbbb"]("b,d,l,j") * r0_1("R");

        // D_oovv_abab += -1.00 r1_bb(d,j) t2_aaaa(a,c,l,i) t1_1_bb(b,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["16_aabb_Lvoov"]("L,a,i,k,d") * r1["bb"]("R,d,j") * t1_1["bb"]("b,k");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["17_aabb_Lvoov"]("L,a,i,j,b")  = t2["aaaa"]("a,c,k,i") * l2["abab"]("L,k,j,c,b");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_abab(a,d,i,l) t2_aaaa(b,c,k,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["17_aabb_Lvoov"]("L,b,j,l,d") * r2["abab"]("R,a,d,i,l");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r2_bbbb(b,d,l,j) t2_aaaa(a,c,k,i) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["17_aabb_Lvoov"]("L,a,i,l,d") * r2["bbbb"]("R,b,d,l,j");

        // D_oovv_abab += -1.00 r0 t2_aaaa(a,c,k,i) t2_bbbb(b,d,l,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["17_aabb_Lvoov"]("L,a,i,l,d") * t2["bbbb"]("b,d,l,j") * r0("R");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["18_aabb_Lvoov"]("L,a,i,j,b")  = t2_1["aaaa"]("a,c,k,i") * l2_1["abab"]("L,k,j,c,b");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_abab(a,d,i,l) t2_1_aaaa(b,c,k,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["18_aabb_Lvoov"]("L,b,j,l,d") * r2["abab"]("R,a,d,i,l");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r0 t2_abab(b,d,j,l) t2_1_aaaa(a,c,k,i) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["18_aabb_Lvoov"]("L,a,i,l,d") * t2["abab"]("b,d,j,l") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r2_bbbb(b,d,l,j) t2_1_aaaa(a,c,k,i) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["18_aabb_Lvoov"]("L,a,i,l,d") * r2["bbbb"]("R,b,d,l,j");

        // D_oovv_abab += -1.00 r0 t2_bbbb(b,d,l,j) t2_1_aaaa(a,c,k,i) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["18_aabb_Lvoov"]("L,a,i,l,d") * t2["bbbb"]("b,d,l,j") * r0("R");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["15_abba_Lvoov"]("L,a,i,j,b")  = t2_1["abab"]("a,c,k,i") * l2_1["abab"]("L,k,j,b,c");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["14_abba_Lvoov"]("L,a,i,j,b")  = t2["abab"]("a,c,k,i") * l2["abab"]("L,k,j,b,c");

        // flops: o3v1L2  = o3v3L1 o3v2L2 o3v2L2 o3v1L2 o3v2L2 o3v1L2 o3v2L2 o3v1L2 o3v2L2 o3v1L2 o3v3L1 o3v2L2 o3v1L2 o3v2L2 o3v1L2 o3v2L2 o3v1L2 o3v3L1 o3v2L2 o3v1L2
        //  mems: o3v1L2  = o2v2L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o2v2L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o2v2L1 o3v1L2 o3v1L2
        tmps_["19_abba_LLvooo"]("L,R,a,i,j,k")  = -1.00 * t2["abab"]("a,d,k,m") * l2["bbbb"]("L,m,j,d,c") * r1["bb"]("R,c,i");
        tmps_["19_abba_LLvooo"]("L,R,a,i,j,k") += tmps_["13_abba_Lvoov"]("L,a,i,j,b") * r1_1["aa"]("R,b,k");
        tmps_["19_abba_LLvooo"]("L,R,a,i,j,k") += tmps_["18_aabb_Lvoov"]("L,a,k,j,c") * r1["bb"]("R,c,i");
        tmps_["19_abba_LLvooo"]("L,R,a,i,j,k") += tmps_["16_aabb_Lvoov"]("L,a,k,j,c") * r1_1["bb"]("R,c,i");
        tmps_["19_abba_LLvooo"]("L,R,a,i,j,k") += tmps_["14_abba_Lvoov"]("L,a,i,j,b") * r1["aa"]("R,b,k");
        tmps_["19_abba_LLvooo"]("L,R,a,i,j,k") -= t2["abab"]("a,d,k,m") * l2_1["bbbb"]("L,m,j,d,c") * r1_1["bb"]("R,c,i");
        tmps_["19_abba_LLvooo"]("L,R,a,i,j,k") += tmps_["17_aabb_Lvoov"]("L,a,k,j,c") * r1["bb"]("R,c,i");
        tmps_["19_abba_LLvooo"]("L,R,a,i,j,k") += tmps_["15_abba_Lvoov"]("L,a,i,j,b") * r1["aa"]("R,b,k");
        tmps_["19_abba_LLvooo"]("L,R,a,i,j,k") -= t2_1["abab"]("a,d,k,m") * l2_1["bbbb"]("L,m,j,d,c") * r1["bb"]("R,c,i");

        // D_ooov_abab += -1.00 r1_1_aa(c,j) t2_abab(a,b,l,i) l2_1_abab(l,k,c,b)
        //             += +1.00 r1_bb(c,i) t2_abab(a,b,j,l) l2_bbbb(l,k,b,c)
        //             += -1.00 r1_bb(c,i) t2_1_aaaa(a,b,l,j) l2_1_abab(l,k,b,c)
        //             += -1.00 r1_1_bb(c,i) t2_aaaa(a,b,l,j) l2_1_abab(l,k,b,c)
        //             += -1.00 r1_aa(c,j) t2_abab(a,b,l,i) l2_abab(l,k,c,b)
        //             += +1.00 r1_1_bb(c,i) t2_abab(a,b,j,l) l2_1_bbbb(l,k,b,c)
        //             += -1.00 r1_bb(c,i) t2_aaaa(a,b,l,j) l2_abab(l,k,b,c)
        //             += -1.00 r1_aa(c,j) t2_1_abab(a,b,l,i) l2_1_abab(l,k,c,b)
        //             += +1.00 r1_bb(c,i) t2_1_abab(a,b,j,l) l2_1_bbbb(l,k,b,c)
        D_ooov_abab("L,R,i,j,k,a") -= tmps_["19_abba_LLvooo"]("L,R,a,i,k,j");

        // D_oovo_abab += +1.00 r1_1_aa(c,j) t2_abab(a,b,l,i) l2_1_abab(l,k,c,b)
        //             += -1.00 r1_bb(c,i) t2_abab(a,b,j,l) l2_bbbb(l,k,b,c)
        //             += +1.00 r1_bb(c,i) t2_1_aaaa(a,b,l,j) l2_1_abab(l,k,b,c)
        //             += +1.00 r1_1_bb(c,i) t2_aaaa(a,b,l,j) l2_1_abab(l,k,b,c)
        //             += +1.00 r1_aa(c,j) t2_abab(a,b,l,i) l2_abab(l,k,c,b)
        //             += -1.00 r1_1_bb(c,i) t2_abab(a,b,j,l) l2_1_bbbb(l,k,b,c)
        //             += +1.00 r1_bb(c,i) t2_aaaa(a,b,l,j) l2_abab(l,k,b,c)
        //             += +1.00 r1_aa(c,j) t2_1_abab(a,b,l,i) l2_1_abab(l,k,c,b)
        //             += -1.00 r1_bb(c,i) t2_1_abab(a,b,j,l) l2_1_bbbb(l,k,b,c)
        D_oovo_abab("L,R,i,j,a,k") += tmps_["19_abba_LLvooo"]("L,R,a,i,k,j");

        // D_oovv_abab += -1.00 r1_1_aa(d,i) t2_abab(a,c,k,j) t1_bb(b,l) l2_1_abab(k,l,d,c)
        //             += +1.00 r1_bb(d,j) t2_abab(a,c,i,k) t1_bb(b,l) l2_bbbb(k,l,c,d)
        //             += -1.00 r1_bb(d,j) t1_bb(b,l) t2_1_aaaa(a,c,k,i) l2_1_abab(k,l,c,d)
        //             += -1.00 r1_1_bb(d,j) t2_aaaa(a,c,k,i) t1_bb(b,l) l2_1_abab(k,l,c,d)
        //             += -1.00 r1_aa(d,i) t2_abab(a,c,k,j) t1_bb(b,l) l2_abab(k,l,d,c)
        //             += +1.00 r1_1_bb(d,j) t2_abab(a,c,i,k) t1_bb(b,l) l2_1_bbbb(k,l,c,d)
        //             += -1.00 r1_bb(d,j) t2_aaaa(a,c,k,i) t1_bb(b,l) l2_abab(k,l,c,d)
        //             += -1.00 r1_aa(d,i) t1_bb(b,l) t2_1_abab(a,c,k,j) l2_1_abab(k,l,d,c)
        //             += +1.00 r1_bb(d,j) t1_bb(b,l) t2_1_abab(a,c,i,k) l2_1_bbbb(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["19_abba_LLvooo"]("L,R,a,j,l,i") * t1["bb"]("b,l");
        tmps_["19_abba_LLvooo"].~TArrayD();

        // flops: o3v1L2  = o3v2L2 o3v2L2 o3v2L2 o3v1L2 o3v2L2 o3v1L2 o3v2L1 o4v1L1 o4v1L2 o3v1L2 o3v2L2 o3v1L2 o3v1L2 o3v2L2 o3v1L2 o4v1L1 o4v1L2 o3v1L2 o3v2L2 o3v1L2 o4v1L2 o3v1L2 o3v2L2 o3v1L2
        //  mems: o3v1L2  = o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L1 o4v0L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o4v0L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        tmps_["20_aaaa_LLovoo"]("L,R,i,a,j,k")  = -1.00 * t2["aaaa"]("a,b,j,k") * tmps_["127_aa_LLov"]("L,R,i,b");
        tmps_["20_aaaa_LLovoo"]("L,R,i,a,j,k") += l1["aa"]("L,i,b") * r2["aaaa"]("R,a,b,j,k");
        tmps_["20_aaaa_LLovoo"]("L,R,i,a,j,k") += l1_1["aa"]("L,i,b") * r2_1["aaaa"]("R,a,b,j,k");
        tmps_["20_aaaa_LLovoo"]("L,R,i,a,j,k") += t2_1["aaaa"]("a,b,j,k") * tmps_["118_aa_LLov"]("L,R,i,b");
        tmps_["20_aaaa_LLovoo"]("L,R,i,a,j,k") += l2["aaaa"]("L,l,i,c,b") * t1["aa"]("b,j") * t1["aa"]("c,k") * r1["aa"]("R,a,l");
        tmps_["20_aaaa_LLovoo"]("L,R,i,a,j,k") += t2["aaaa"]("a,b,j,k") * tmps_["129_aa_LLov"]("L,R,i,b");
        tmps_["20_aaaa_LLovoo"]("L,R,i,a,j,k") += t2["aaaa"]("a,b,j,k") * tmps_["130_aa_LLov"]("L,R,i,b");
        tmps_["20_aaaa_LLovoo"]("L,R,i,a,j,k") += t1["aa"]("c,k") * tmps_["262_aaaa_Looov"]("L,j,l,i,c") * r1_1["aa"]("R,a,l");
        tmps_["20_aaaa_LLovoo"]("L,R,i,a,j,k") -= t2["aaaa"]("a,b,j,k") * tmps_["128_aa_LLov"]("L,R,i,b");
        tmps_["20_aaaa_LLovoo"]("L,R,i,a,j,k") -= 0.50 * tmps_["173_aaaa_Loooo"]("L,j,k,l,i") * r1_1["aa"]("R,a,l");
        tmps_["20_aaaa_LLovoo"]("L,R,i,a,j,k") -= t2_1["aaaa"]("a,b,j,k") * tmps_["117_aa_LLov"]("L,R,i,b");

        // D_ooov_aaaa += -1.00 r2_1_aaaa(a,b,i,j) l1_1_aa(k,b)
        //             += -1.00 r2_aaaa(a,b,i,j) l1_aa(k,b)
        //             += -1.00 r1_bb(c,l) t2_1_aaaa(a,b,i,j) l2_1_abab(k,l,b,c)
        //             += -1.00 r1_aa(a,l) t1_aa(b,i) t1_aa(c,j) l2_aaaa(l,k,c,b)
        //             += -1.00 r1_bb(c,l) t2_aaaa(a,b,i,j) l2_abab(k,l,b,c)
        //             += +1.00 r1_aa(c,l) t2_aaaa(a,b,i,j) l2_aaaa(l,k,b,c)
        //             += -1.00 r1_1_bb(c,l) t2_aaaa(a,b,i,j) l2_1_abab(k,l,b,c)
        //             += -1.00 r1_1_aa(a,l) t1_aa(b,i) t1_aa(c,j) l2_1_aaaa(l,k,c,b)
        //             += +1.00 r1_1_aa(c,l) t2_aaaa(a,b,i,j) l2_1_aaaa(l,k,b,c)
        //             += +0.50 r1_1_aa(a,l) t2_aaaa(c,b,i,j) l2_1_aaaa(l,k,c,b)
        //             += +1.00 r1_aa(c,l) t2_1_aaaa(a,b,i,j) l2_1_aaaa(l,k,b,c)
        D_ooov_aaaa("L,R,i,j,k,a") -= tmps_["20_aaaa_LLovoo"]("L,R,k,a,i,j");

        // D_oovo_aaaa += +1.00 r2_1_aaaa(a,b,i,j) l1_1_aa(k,b)
        //             += +1.00 r2_aaaa(a,b,i,j) l1_aa(k,b)
        //             += +1.00 r1_bb(c,l) t2_1_aaaa(a,b,i,j) l2_1_abab(k,l,b,c)
        //             += +1.00 r1_aa(a,l) t1_aa(b,i) t1_aa(c,j) l2_aaaa(l,k,c,b)
        //             += +1.00 r1_bb(c,l) t2_aaaa(a,b,i,j) l2_abab(k,l,b,c)
        //             += -1.00 r1_aa(c,l) t2_aaaa(a,b,i,j) l2_aaaa(l,k,b,c)
        //             += +1.00 r1_1_bb(c,l) t2_aaaa(a,b,i,j) l2_1_abab(k,l,b,c)
        //             += +1.00 r1_1_aa(a,l) t1_aa(b,i) t1_aa(c,j) l2_1_aaaa(l,k,c,b)
        //             += -1.00 r1_1_aa(c,l) t2_aaaa(a,b,i,j) l2_1_aaaa(l,k,b,c)
        //             += -0.50 r1_1_aa(a,l) t2_aaaa(c,b,i,j) l2_1_aaaa(l,k,c,b)
        //             += -1.00 r1_aa(c,l) t2_1_aaaa(a,b,i,j) l2_1_aaaa(l,k,b,c)
        D_oovo_aaaa("L,R,i,j,a,k") += tmps_["20_aaaa_LLovoo"]("L,R,k,a,i,j");

        // D_oovv_aaaa += +1.00 P(a,b) r2_1_aaaa(a,c,i,j) t1_aa(b,k) l1_1_aa(k,c)
        //             += +1.00 P(a,b) r2_aaaa(a,c,i,j) t1_aa(b,k) l1_aa(k,c)
        //             += +1.00 P(a,b) r1_bb(d,l) t1_aa(b,k) t2_1_aaaa(a,c,i,j) l2_1_abab(k,l,c,d)
        //             += +1.00 P(a,b) r1_aa(a,l) t1_aa(b,k) t1_aa(c,i) t1_aa(d,j) l2_aaaa(l,k,d,c)
        //             += +1.00 P(a,b) r1_bb(d,l) t2_aaaa(a,c,i,j) t1_aa(b,k) l2_abab(k,l,c,d)
        //             += -1.00 P(a,b) r1_aa(d,l) t2_aaaa(a,c,i,j) t1_aa(b,k) l2_aaaa(l,k,c,d)
        //             += +1.00 P(a,b) r1_1_bb(d,l) t2_aaaa(a,c,i,j) t1_aa(b,k) l2_1_abab(k,l,c,d)
        //             += +1.00 P(a,b) r1_1_aa(a,l) t1_aa(b,k) t1_aa(c,i) t1_aa(d,j) l2_1_aaaa(l,k,d,c)
        //             += -1.00 P(a,b) r1_1_aa(d,l) t2_aaaa(a,c,i,j) t1_aa(b,k) l2_1_aaaa(l,k,c,d)
        //             += -0.50 P(a,b) r1_1_aa(a,l) t1_aa(b,k) t2_aaaa(d,c,i,j) l2_1_aaaa(l,k,d,c)
        //             += -1.00 P(a,b) r1_aa(d,l) t1_aa(b,k) t2_1_aaaa(a,c,i,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["20_aaaa_LLovoo"]("L,R,k,a,i,j") * t1["aa"]("b,k");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();
        tmps_["20_aaaa_LLovoo"].~TArrayD();

        // flops: o1v3L2  = o2v3L2 o2v3L2 o1v3L2
        //  mems: o1v3L2  = o1v3L2 o1v3L2 o1v3L2
        tmps_["21_aaaa_LLovvv"]("L,R,i,a,b,c")  = l2["aaaa"]("L,j,i,a,b") * r1["aa"]("R,c,j");
        tmps_["21_aaaa_LLovvv"]("L,R,i,a,b,c") += l2_1["aaaa"]("L,j,i,a,b") * r1_1["aa"]("R,c,j");

        // D_vvvo_aaaa  = -1.00 r1_1_aa(c,j) l2_1_aaaa(j,i,a,b)
        //             += -1.00 r1_aa(c,j) l2_aaaa(j,i,a,b)
        D_vvvo_aaaa("L,R,a,b,c,i")  = -1.00 * tmps_["21_aaaa_LLovvv"]("L,R,i,a,b,c");

        // D_vvov_aaaa  = +1.00 r1_1_aa(c,j) l2_1_aaaa(j,i,a,b)
        //             += +1.00 r1_aa(c,j) l2_aaaa(j,i,a,b)
        D_vvov_aaaa("L,R,a,b,i,c")  = tmps_["21_aaaa_LLovvv"]("L,R,i,a,b,c");

        // D_vvvv_aaaa += -1.00 P(c,d) r1_1_aa(c,j) t1_aa(d,i) l2_1_aaaa(j,i,a,b)
        //             += -1.00 P(c,d) r1_aa(c,j) t1_aa(d,i) l2_aaaa(j,i,a,b)
        // flops: o0v4L2 += o1v4L2
        //  mems: o0v4L2 += o0v4L2
        tmps_["perm_aaaa_LLvvvv"]("L,R,a,b,c,d")  = tmps_["21_aaaa_LLovvv"]("L,R,i,a,b,c") * t1["aa"]("d,i");
        D_vvvv_aaaa("L,R,a,b,c,d") -= tmps_["perm_aaaa_LLvvvv"]("L,R,a,b,c,d");
        D_vvvv_aaaa("L,R,a,b,c,d") += tmps_["perm_aaaa_LLvvvv"]("L,R,a,b,d,c");
        tmps_["perm_aaaa_LLvvvv"].~TArrayD();
        tmps_["21_aaaa_LLovvv"].~TArrayD();

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["22_aaaa_Looov"]("L,i,j,k,a")  = t1["aa"]("b,i") * l2_1["aaaa"]("L,j,k,b,a");

        // D_oooo_aaaa += -1.00 P(i,j) r1_1_aa(b,i) t1_aa(a,j) l2_1_aaaa(l,k,a,b)
        // flops: o4v0L2 += o4v1L2
        //  mems: o4v0L2 += o4v0L2
        tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l")  = tmps_["22_aaaa_Looov"]("L,j,l,k,b") * r1_1["aa"]("R,b,i");
        D_oooo_aaaa("L,R,i,j,k,l") -= tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l");
        D_oooo_aaaa("L,R,i,j,k,l") += tmps_["perm_aaaa_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_aaaa_LLoooo"].~TArrayD();

        // D_oooo_aaaa += -1.00 P(i,j) r0 t1_aa(b,j) t1_1_aa(a,i) l2_1_aaaa(l,k,b,a)
        // flops: o4v0L2 += o4v1L1 o4v0L2
        //  mems: o4v0L2 += o4v0L1 o4v0L2
        tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l")  = t1_1["aa"]("a,i") * tmps_["22_aaaa_Looov"]("L,j,l,k,a") * r0("R");
        D_oooo_aaaa("L,R,i,j,k,l") -= tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l");
        D_oooo_aaaa("L,R,i,j,k,l") += tmps_["perm_aaaa_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_aaaa_LLoooo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r1_aa(a,l) t1_aa(d,j) t2_1_aaaa(b,c,k,i) l2_1_aaaa(l,k,d,c)
        // flops: o2v2L2 += o4v2L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2_1["aaaa"]("b,c,k,i") * tmps_["22_aaaa_Looov"]("L,j,l,k,c") * r1["aa"]("R,a,l");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +1.00 r1_aa(a,l) t1_aa(d,i) t2_1_abab(c,b,k,j) l2_1_aaaa(l,k,d,c)
        // flops: o2v2L2 += o4v2L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2_1["abab"]("c,b,k,j") * tmps_["22_aaaa_Looov"]("L,i,l,k,c") * r1["aa"]("R,a,l");

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["23_aaaa_Looov"]("L,i,j,k,a")  = t1["aa"]("b,i") * l2["aaaa"]("L,j,k,b,a");

        // D_oooo_aaaa += -1.00 P(i,j) r1_aa(b,i) t1_aa(a,j) l2_aaaa(l,k,a,b)
        // flops: o4v0L2 += o4v1L2
        //  mems: o4v0L2 += o4v0L2
        tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l")  = tmps_["23_aaaa_Looov"]("L,j,l,k,b") * r1["aa"]("R,b,i");
        D_oooo_aaaa("L,R,i,j,k,l") -= tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l");
        D_oooo_aaaa("L,R,i,j,k,l") += tmps_["perm_aaaa_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_aaaa_LLoooo"].~TArrayD();

        // D_oooo_aaaa += -1.00 r0 t1_aa(a,i) t1_aa(b,j) l2_aaaa(l,k,b,a)
        // flops: o4v0L2 += o4v1L1 o4v0L2
        //  mems: o4v0L2 += o4v0L1 o4v0L2
        D_oooo_aaaa("L,R,i,j,k,l") -= t1["aa"]("a,i") * tmps_["23_aaaa_Looov"]("L,j,l,k,a") * r0("R");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r1_aa(a,l) t2_aaaa(b,c,k,i) t1_aa(d,j) l2_aaaa(l,k,d,c)
        // flops: o2v2L2 += o4v2L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["aaaa"]("b,c,k,i") * tmps_["23_aaaa_Looov"]("L,j,l,k,c") * r1["aa"]("R,a,l");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +1.00 r1_aa(a,l) t2_abab(c,b,k,j) t1_aa(d,i) l2_aaaa(l,k,d,c)
        // flops: o2v2L2 += o4v2L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("c,b,k,j") * tmps_["23_aaaa_Looov"]("L,i,l,k,c") * r1["aa"]("R,a,l");

        // flops: o3v1L2  = o3v2L1 o4v1L2 o4v1L2 o4v1L2 o4v1L2 o3v1L2 o4v1L2 o4v1L2 o3v1L2
        //  mems: o3v1L2  = o3v1L1 o4v0L2 o3v1L2 o4v0L2 o3v1L2 o3v1L2 o4v0L2 o3v1L2 o3v1L2
        tmps_["24_aaaa_LLooov"]("L,R,i,j,k,a")  = t1_1["aa"]("c,i") * l2_1["aaaa"]("L,l,j,c,b") * r1["aa"]("R,b,k") * t1["aa"]("a,l");
        tmps_["24_aaaa_LLooov"]("L,R,i,j,k,a") += tmps_["23_aaaa_Looov"]("L,i,l,j,b") * r1["aa"]("R,b,k") * t1["aa"]("a,l");
        tmps_["24_aaaa_LLooov"]("L,R,i,j,k,a") += tmps_["22_aaaa_Looov"]("L,i,l,j,b") * r1_1["aa"]("R,b,k") * t1["aa"]("a,l");

        // D_ooov_aaaa += -1.00 P(i,j) r1_aa(c,i) t1_aa(a,l) t1_1_aa(b,j) l2_1_aaaa(l,k,b,c)
        //             += -1.00 P(i,j) r1_aa(c,i) t1_aa(a,l) t1_aa(b,j) l2_aaaa(l,k,b,c)
        //             += -1.00 P(i,j) r1_1_aa(c,i) t1_aa(a,l) t1_aa(b,j) l2_1_aaaa(l,k,b,c)
        D_ooov_aaaa("L,R,i,j,k,a") -= tmps_["24_aaaa_LLooov"]("L,R,j,k,i,a");
        D_ooov_aaaa("L,R,i,j,k,a") += tmps_["24_aaaa_LLooov"]("L,R,i,k,j,a");

        // D_oovo_aaaa += +1.00 P(i,j) r1_aa(c,i) t1_aa(a,l) t1_1_aa(b,j) l2_1_aaaa(l,k,b,c)
        //             += +1.00 P(i,j) r1_aa(c,i) t1_aa(a,l) t1_aa(b,j) l2_aaaa(l,k,b,c)
        //             += +1.00 P(i,j) r1_1_aa(c,i) t1_aa(a,l) t1_aa(b,j) l2_1_aaaa(l,k,b,c)
        D_oovo_aaaa("L,R,i,j,a,k") += tmps_["24_aaaa_LLooov"]("L,R,j,k,i,a");
        D_oovo_aaaa("L,R,i,j,a,k") -= tmps_["24_aaaa_LLooov"]("L,R,i,k,j,a");

        // D_oovv_aaaa += +1.00 P(i,j) r1_aa(d,i) t1_aa(a,k) t1_aa(b,l) t1_1_aa(c,j) l2_1_aaaa(k,l,c,d)
        //             += +1.00 P(i,j) r1_aa(d,i) t1_aa(a,k) t1_aa(b,l) t1_aa(c,j) l2_aaaa(k,l,c,d)
        //             += +1.00 P(i,j) r1_1_aa(d,i) t1_aa(a,k) t1_aa(b,l) t1_aa(c,j) l2_1_aaaa(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("b,l") * tmps_["24_aaaa_LLooov"]("L,R,j,l,i,a");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();
        tmps_["24_aaaa_LLooov"].~TArrayD();

        // flops: o3v1L2  = o3v2L1 o4v1L2 o4v1L2 o3v2L1 o4v1L2 o4v1L2 o3v1L2 o3v2L1 o4v1L2 o4v1L2 o3v1L2
        //  mems: o3v1L2  = o3v1L1 o4v0L2 o3v1L2 o3v1L1 o4v0L2 o3v1L2 o3v1L2 o3v1L1 o4v0L2 o3v1L2 o3v1L2
        tmps_["25_bbbb_LLooov"]("L,R,i,j,k,a")  = t1["bb"]("c,i") * l2["bbbb"]("L,l,j,c,b") * r1["bb"]("R,b,k") * t1["bb"]("a,l");
        tmps_["25_bbbb_LLooov"]("L,R,i,j,k,a") += t1["bb"]("c,i") * l2_1["bbbb"]("L,l,j,c,b") * r1_1["bb"]("R,b,k") * t1["bb"]("a,l");
        tmps_["25_bbbb_LLooov"]("L,R,i,j,k,a") += t1_1["bb"]("c,i") * l2_1["bbbb"]("L,l,j,c,b") * r1["bb"]("R,b,k") * t1["bb"]("a,l");

        // D_ooov_bbbb += -1.00 P(i,j) r1_1_bb(c,i) t1_bb(a,l) t1_bb(b,j) l2_1_bbbb(l,k,b,c)
        //             += -1.00 P(i,j) r1_bb(c,i) t1_bb(a,l) t1_bb(b,j) l2_bbbb(l,k,b,c)
        //             += -1.00 P(i,j) r1_bb(c,i) t1_bb(a,l) t1_1_bb(b,j) l2_1_bbbb(l,k,b,c)
        D_ooov_bbbb("L,R,i,j,k,a") -= tmps_["25_bbbb_LLooov"]("L,R,j,k,i,a");
        D_ooov_bbbb("L,R,i,j,k,a") += tmps_["25_bbbb_LLooov"]("L,R,i,k,j,a");

        // D_oovo_bbbb += +1.00 P(i,j) r1_1_bb(c,i) t1_bb(a,l) t1_bb(b,j) l2_1_bbbb(l,k,b,c)
        //             += +1.00 P(i,j) r1_bb(c,i) t1_bb(a,l) t1_bb(b,j) l2_bbbb(l,k,b,c)
        //             += +1.00 P(i,j) r1_bb(c,i) t1_bb(a,l) t1_1_bb(b,j) l2_1_bbbb(l,k,b,c)
        D_oovo_bbbb("L,R,i,j,a,k") += tmps_["25_bbbb_LLooov"]("L,R,j,k,i,a");
        D_oovo_bbbb("L,R,i,j,a,k") -= tmps_["25_bbbb_LLooov"]("L,R,i,k,j,a");

        // D_oovv_bbbb += +1.00 P(i,j) r1_1_bb(d,i) t1_bb(a,k) t1_bb(b,l) t1_bb(c,j) l2_1_bbbb(k,l,c,d)
        //             += +1.00 P(i,j) r1_bb(d,i) t1_bb(a,k) t1_bb(b,l) t1_bb(c,j) l2_bbbb(k,l,c,d)
        //             += +1.00 P(i,j) r1_bb(d,i) t1_bb(a,k) t1_bb(b,l) t1_1_bb(c,j) l2_1_bbbb(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("b,l") * tmps_["25_bbbb_LLooov"]("L,R,j,l,i,a");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
        tmps_["25_bbbb_LLooov"].~TArrayD();

        // flops: o1v3L2  = o3v3L2 o3v3L2 o2v2L2 o3v3L2 o2v2L2 o3v3L2 o2v2L2 o2v3L2
        //  mems: o1v3L2  = o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o1v3L2
        tmps_["26_aabb_LLvvvo"]("L,R,a,b,c,i")  = t1["aa"]("a,j") * (l2["aaaa"]("L,l,j,b,e") * r2["abab"]("R,e,c,l,i") + l2["abab"]("L,j,k,b,d") * r2["bbbb"]("R,c,d,k,i") + l2_1["aaaa"]("L,l,j,b,e") * r2_1["abab"]("R,e,c,l,i") + l2_1["abab"]("L,j,k,b,d") * r2_1["bbbb"]("R,c,d,k,i"));

        // D_oovv_abab += -1.00 r2_abab(d,b,l,j) t1_aa(a,k) t1_aa(c,i) l2_aaaa(l,k,c,d)
        //             += -1.00 r2_bbbb(b,d,l,j) t1_aa(a,k) t1_aa(c,i) l2_abab(k,l,c,d)
        //             += -1.00 r2_1_abab(d,b,l,j) t1_aa(a,k) t1_aa(c,i) l2_1_aaaa(l,k,c,d)
        //             += -1.00 r2_1_bbbb(b,d,l,j) t1_aa(a,k) t1_aa(c,i) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["26_aabb_LLvvvo"]("L,R,a,c,b,j") * t1["aa"]("c,i");
        tmps_["26_aabb_LLvvvo"].~TArrayD();

        // flops: o1v1L2  = o2v2L2 o2v2L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2
        tmps_["27_aa_LLvo"]("L,R,a,i")  = -1.00 * l1_1["bb"]("L,k,c") * r2["abab"]("R,a,c,i,k");
        tmps_["27_aa_LLvo"]("L,R,a,i") += l1_1["aa"]("L,j,b") * r2["aaaa"]("R,a,b,j,i");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_aaaa(a,c,k,i) t1_1_aa(b,j) l1_1_aa(k,c)
        //             += -1.00 P(i,j) P(a,b) r2_abab(a,c,i,k) t1_1_aa(b,j) l1_1_bb(k,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["27_aa_LLvo"]("L,R,a,i") * t1_1["aa"]("b,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +1.00 r2_aaaa(a,c,k,i) t1_1_bb(b,j) l1_1_aa(k,c)
        //             += -1.00 r2_abab(a,c,i,k) t1_1_bb(b,j) l1_1_bb(k,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["27_aa_LLvo"]("L,R,a,i") * t1_1["bb"]("b,j");
        tmps_["27_aa_LLvo"].~TArrayD();

        // flops: o1v1L2  = o2v2L2 o2v2L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2
        tmps_["28_bb_LLvo"]("L,R,a,i")  = -1.00 * l1_1["bb"]("L,k,c") * r2["bbbb"]("R,a,c,k,i");
        tmps_["28_bb_LLvo"]("L,R,a,i") += l1_1["aa"]("L,j,b") * r2["abab"]("R,b,a,j,i");

        // D_oovv_abab += -1.00 r2_abab(c,b,k,j) t1_1_aa(a,i) l1_1_aa(k,c)
        //             += +1.00 r2_bbbb(b,c,k,j) t1_1_aa(a,i) l1_1_bb(k,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1_1["aa"]("a,i") * tmps_["28_bb_LLvo"]("L,R,b,j");

        // D_oovv_bbbb += -1.00 P(i,j) P(a,b) r2_abab(c,a,k,i) t1_1_bb(b,j) l1_1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r2_bbbb(a,c,k,i) t1_1_bb(b,j) l1_1_bb(k,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["28_bb_LLvo"]("L,R,a,i") * t1_1["bb"]("b,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
        tmps_["28_bb_LLvo"].~TArrayD();

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["29_aa_Lvv"]("L,a,b")  = l2_1["abab"]("L,i,j,a,c") * t2["abab"]("b,c,i,j");

        // D_oovv_aaaa += +0.50 P(a,b) r2_1_aaaa(a,d,i,j) t2_abab(b,c,k,l) l2_1_abab(k,l,d,c)
        //             += +0.50 P(a,b) r2_1_aaaa(a,d,i,j) t2_abab(b,c,l,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["29_aa_Lvv"]("L,d,b") * r2_1["aaaa"]("R,a,d,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += -0.50 r0 t2_abab(a,c,k,l) t2_aaaa(b,d,i,j) l2_abab(k,l,d,c)
        //             += -0.50 r0 t2_abab(a,c,l,k) t2_aaaa(b,d,i,j) l2_abab(l,k,d,c)
        //             += -0.50 r0_1 t2_abab(a,c,k,l) t2_aaaa(b,d,i,j) l2_1_abab(k,l,d,c)
        //             += -0.50 r0_1 t2_abab(a,c,l,k) t2_aaaa(b,d,i,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L1 o0v2L2 o0v2L2 o0v2L2 o2v3L2
        //  mems: o2v2L2 += o0v2L1 o0v2L2 o0v2L2 o0v2L2 o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") -= (t2["abab"]("a,c,k,l") * l2["abab"]("L,k,l,d,c") * r0("R") + tmps_["29_aa_Lvv"]("L,d,a") * r0_1("R")) * t2["aaaa"]("b,d,i,j");

        // D_oovv_aaaa += +0.50 P(a,b) r0 t2_abab(b,d,k,l) t2_1_aaaa(a,c,i,j) l2_1_abab(k,l,c,d)
        //             += +0.50 P(a,b) r0 t2_abab(b,d,l,k) t2_1_aaaa(a,c,i,j) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2_1["aaaa"]("a,c,i,j") * tmps_["29_aa_Lvv"]("L,c,b") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r2_1_abab(d,b,i,j) t2_abab(a,c,k,l) l2_1_abab(k,l,d,c)
        //             += +0.50 r2_1_abab(d,b,i,j) t2_abab(a,c,l,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["29_aa_Lvv"]("L,d,a") * r2_1["abab"]("R,d,b,i,j");

        // D_oovv_abab += +0.50 r0 t2_abab(a,d,k,l) t2_1_abab(c,b,i,j) l2_1_abab(k,l,c,d)
        //             += +0.50 r0 t2_abab(a,d,l,k) t2_1_abab(c,b,i,j) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["29_aa_Lvv"]("L,c,a") * t2_1["abab"]("c,b,i,j") * r0("R");

        // D_oovv_abab += +0.50 r0_1 t2_abab(a,c,k,l) t2_abab(d,b,i,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r0_1 t2_abab(a,c,l,k) t2_abab(d,b,i,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["29_aa_Lvv"]("L,d,a") * t2["abab"]("d,b,i,j") * r0_1("R");

        // flops: o2v2L2  = o2v3L1 o0v2L2 o0v2L2 o0v2L2 o2v3L2
        //  mems: o2v2L2  = o0v2L1 o0v2L2 o0v2L2 o0v2L2 o2v2L2
        tmps_["30_aaaa_LLvoov"]("L,R,a,i,j,b")  = t2["aaaa"]("a,c,i,j") * (t2["abab"]("b,d,k,l") * l2["abab"]("L,k,l,c,d") * r0("R") + tmps_["29_aa_Lvv"]("L,c,b") * r0_1("R"));

        // D_oovv_aaaa += +0.50 r0 t2_aaaa(a,c,i,j) t2_abab(b,d,k,l) l2_abab(k,l,c,d)
        //             += +0.50 r0 t2_aaaa(a,c,i,j) t2_abab(b,d,l,k) l2_abab(l,k,c,d)
        //             += +0.50 r0_1 t2_aaaa(a,c,i,j) t2_abab(b,d,k,l) l2_1_abab(k,l,c,d)
        //             += +0.50 r0_1 t2_aaaa(a,c,i,j) t2_abab(b,d,l,k) l2_1_abab(l,k,c,d)
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["30_aaaa_LLvoov"]("L,R,a,i,j,b");
        tmps_["30_aaaa_LLvoov"].~TArrayD();

        // flops: o1v3L2  = o2v3L1 o1v3L2 o2v3L2 o2v3L1 o1v3L2 o2v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        //  mems: o1v3L2  = o1v3L1 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        tmps_["31_abab_LLvovv"]("L,R,a,i,b,c")  = t1_1["aa"]("a,j") * l2_1["abab"]("L,j,i,b,c") * r0("R");
        tmps_["31_abab_LLvovv"]("L,R,a,i,b,c") += l2_1["abab"]("L,j,i,b,c") * r1_1["aa"]("R,a,j");
        tmps_["31_abab_LLvovv"]("L,R,a,i,b,c") += t1["aa"]("a,j") * l2["abab"]("L,j,i,b,c") * r0("R");
        tmps_["31_abab_LLvovv"]("L,R,a,i,b,c") += l2["abab"]("L,j,i,b,c") * r1["aa"]("R,a,j");
        tmps_["31_abab_LLvovv"]("L,R,a,i,b,c") += t1["aa"]("a,j") * l2_1["abab"]("L,j,i,b,c") * r0_1("R");

        // D_vvov_abab  = +1.00 r0 t1_aa(c,j) l2_abab(j,i,a,b)
        //             += +1.00 r1_aa(c,j) l2_abab(j,i,a,b)
        //             += +1.00 r0_1 t1_aa(c,j) l2_1_abab(j,i,a,b)
        //             += +1.00 r1_1_aa(c,j) l2_1_abab(j,i,a,b)
        //             += +1.00 r0 t1_1_aa(c,j) l2_1_abab(j,i,a,b)
        D_vvov_abab("L,R,a,b,i,c")  = tmps_["31_abab_LLvovv"]("L,R,c,i,a,b");

        // D_vvvo_abab  = -1.00 r0 t1_aa(c,j) l2_abab(j,i,a,b)
        //             += -1.00 r1_aa(c,j) l2_abab(j,i,a,b)
        //             += -1.00 r0_1 t1_aa(c,j) l2_1_abab(j,i,a,b)
        //             += -1.00 r1_1_aa(c,j) l2_1_abab(j,i,a,b)
        //             += -1.00 r0 t1_1_aa(c,j) l2_1_abab(j,i,a,b)
        D_vvvo_abab("L,R,a,b,c,i")  = -1.00 * tmps_["31_abab_LLvovv"]("L,R,c,i,a,b");
        tmps_["31_abab_LLvovv"].~TArrayD();

        // flops: o3v1L1  = o3v2L1 o3v2L1 o3v1L1
        //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1
        tmps_["32_baba_Lvooo"]("L,a,i,j,k")  = t2["abab"]("b,a,i,j") * l1["aa"]("L,k,b");
        tmps_["32_baba_Lvooo"]("L,a,i,j,k") += t2_1["abab"]("b,a,i,j") * l1_1["aa"]("L,k,b");

        // D_oovv_abab += +1.00 r1_aa(a,k) t2_1_abab(c,b,i,j) l1_1_aa(k,c)
        //             += +1.00 r1_aa(a,k) t2_abab(c,b,i,j) l1_aa(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["32_baba_Lvooo"]("L,b,i,j,k") * r1["aa"]("R,a,k");

        // D_oovv_abab += +1.00 r0 t1_aa(a,k) t2_1_abab(c,b,i,j) l1_1_aa(k,c)
        //             += +1.00 r0 t1_aa(a,k) t2_abab(c,b,i,j) l1_aa(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1["aa"]("a,k") * tmps_["32_baba_Lvooo"]("L,b,i,j,k") * r0("R");
        tmps_["32_baba_Lvooo"].~TArrayD();

        // flops: o3v1L2  = o4v2L1 o4v2L1 o4v0L1 o4v1L2
        //  mems: o3v1L2  = o4v0L1 o4v0L1 o4v0L1 o3v1L2
        tmps_["33_aaaa_LLooov"]("L,R,i,j,k,a")  = (t2_1["aaaa"]("b,c,i,j") * l2_1["aaaa"]("L,l,k,b,c") + t2["aaaa"]("b,c,i,j") * l2["aaaa"]("L,l,k,b,c")) * r1["aa"]("R,a,l");

        // D_oovv_aaaa += -0.50 P(a,b) r1_aa(a,l) t1_aa(b,k) t2_1_aaaa(d,c,i,j) l2_1_aaaa(l,k,d,c)
        //             += -0.50 P(a,b) r1_aa(a,l) t1_aa(b,k) t2_aaaa(d,c,i,j) l2_aaaa(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t1["aa"]("b,k") * tmps_["33_aaaa_LLooov"]("L,R,i,j,k,a");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();
        tmps_["33_aaaa_LLooov"].~TArrayD();

        // flops: o3v1L2  = o4v2L1 o4v2L1 o4v0L1 o4v1L2
        //  mems: o3v1L2  = o4v0L1 o4v0L1 o4v0L1 o3v1L2
        tmps_["34_bbbb_LLooov"]("L,R,i,j,k,a")  = (t2_1["bbbb"]("b,c,i,j") * l2_1["bbbb"]("L,l,k,b,c") + t2["bbbb"]("b,c,i,j") * l2["bbbb"]("L,l,k,b,c")) * r1["bb"]("R,a,l");

        // D_oovv_bbbb += -0.50 P(a,b) r1_bb(a,l) t1_bb(b,k) t2_1_bbbb(d,c,i,j) l2_1_bbbb(l,k,d,c)
        //             += -0.50 P(a,b) r1_bb(a,l) t1_bb(b,k) t2_bbbb(d,c,i,j) l2_bbbb(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t1["bb"]("b,k") * tmps_["34_bbbb_LLooov"]("L,R,i,j,k,a");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
        tmps_["34_bbbb_LLooov"].~TArrayD();

        // flops: o0v4L2  = o2v4L2 o2v4L2 o0v4L2
        //  mems: o0v4L2  = o0v4L2 o0v4L2 o0v4L2
        tmps_["35_aaaa_LLvvvv"]("L,R,a,b,c,d")  = l2["aaaa"]("L,i,j,a,b") * r2["aaaa"]("R,c,d,i,j");
        tmps_["35_aaaa_LLvvvv"]("L,R,a,b,c,d") += l2_1["aaaa"]("L,i,j,a,b") * r2_1["aaaa"]("R,c,d,i,j");

        // D_vvvv_aaaa += +0.50 r2_1_aaaa(d,c,j,i) l2_1_aaaa(j,i,a,b)
        //             += +0.50 r2_aaaa(d,c,j,i) l2_aaaa(j,i,a,b)
        D_vvvv_aaaa("L,R,a,b,c,d") += 0.50 * tmps_["35_aaaa_LLvvvv"]("L,R,a,b,d,c");

        // flops: o1v3L2  = o1v4L2
        //  mems: o1v3L2  = o1v3L2
        tmps_["36_aaaa_LLvvvo"]("L,R,a,b,c,i")  = tmps_["35_aaaa_LLvvvv"]("L,R,a,d,b,c") * t1["aa"]("d,i");
        tmps_["35_aaaa_LLvvvv"].~TArrayD();

        // D_oovv_aaaa += -0.50 r2_1_aaaa(b,a,l,k) t1_aa(c,i) t1_aa(d,j) l2_1_aaaa(l,k,d,c)
        //             += -0.50 r2_aaaa(b,a,l,k) t1_aa(c,i) t1_aa(d,j) l2_aaaa(l,k,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") -= 0.50 * tmps_["36_aaaa_LLvvvo"]("L,R,d,b,a,i") * t1["aa"]("d,j");

        // D_ovvv_aaaa += +0.50 r2_1_aaaa(c,b,k,j) t1_aa(d,i) l2_1_aaaa(k,j,a,d)
        //             += +0.50 r2_aaaa(c,b,k,j) t1_aa(d,i) l2_aaaa(k,j,a,d)
        D_ovvv_aaaa("L,R,i,a,b,c") += 0.50 * tmps_["36_aaaa_LLvvvo"]("L,R,a,c,b,i");

        // D_vovv_aaaa += -0.50 r2_1_aaaa(c,b,k,j) t1_aa(d,i) l2_1_aaaa(k,j,a,d)
        //             += -0.50 r2_aaaa(c,b,k,j) t1_aa(d,i) l2_aaaa(k,j,a,d)
        D_vovv_aaaa("L,R,a,i,b,c") -= 0.50 * tmps_["36_aaaa_LLvvvo"]("L,R,a,c,b,i");
        tmps_["36_aaaa_LLvvvo"].~TArrayD();

        // flops: o0v4L2  = o2v4L2 o2v4L2 o0v4L2
        //  mems: o0v4L2  = o0v4L2 o0v4L2 o0v4L2
        tmps_["37_abab_LLvvvv"]("L,R,a,b,c,d")  = l2["abab"]("L,i,j,a,b") * r2["abab"]("R,c,d,i,j");
        tmps_["37_abab_LLvvvv"]("L,R,a,b,c,d") += l2_1["abab"]("L,i,j,a,b") * r2_1["abab"]("R,c,d,i,j");

        // D_vvvv_abab += -0.50 r2_1_abab(c,d,j,i) l2_1_abab(j,i,a,b)
        //             += -0.50 r2_1_abab(c,d,i,j) l2_1_abab(i,j,a,b)
        //             += -0.50 r2_abab(c,d,j,i) l2_abab(j,i,a,b)
        //             += -0.50 r2_abab(c,d,i,j) l2_abab(i,j,a,b)
        D_vvvv_abab("L,R,a,b,c,d") -= tmps_["37_abab_LLvvvv"]("L,R,a,b,c,d");

        // flops: o1v3L2  = o1v4L2
        //  mems: o1v3L2  = o1v3L2
        tmps_["38_aabb_LLvvvo"]("L,R,a,b,c,i")  = tmps_["37_abab_LLvvvv"]("L,R,a,d,b,c") * t1["bb"]("d,i");
        tmps_["37_abab_LLvvvv"].~TArrayD();

        // D_oovv_abab += -0.50 r2_1_abab(a,b,l,k) t1_aa(c,i) t1_bb(d,j) l2_1_abab(l,k,c,d)
        //             += -0.50 r2_1_abab(a,b,k,l) t1_aa(c,i) t1_bb(d,j) l2_1_abab(k,l,c,d)
        //             += -0.50 r2_abab(a,b,l,k) t1_aa(c,i) t1_bb(d,j) l2_abab(l,k,c,d)
        //             += -0.50 r2_abab(a,b,k,l) t1_aa(c,i) t1_bb(d,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["38_aabb_LLvvvo"]("L,R,c,a,b,j") * t1["aa"]("c,i");

        // D_ovvv_abab += +0.50 r2_1_abab(c,b,k,j) t1_bb(d,i) l2_1_abab(k,j,a,d)
        //             += +0.50 r2_1_abab(c,b,j,k) t1_bb(d,i) l2_1_abab(j,k,a,d)
        //             += +0.50 r2_abab(c,b,k,j) t1_bb(d,i) l2_abab(k,j,a,d)
        //             += +0.50 r2_abab(c,b,j,k) t1_bb(d,i) l2_abab(j,k,a,d)
        D_ovvv_abab("L,R,i,a,b,c") += tmps_["38_aabb_LLvvvo"]("L,R,a,c,b,i");

        // D_vovv_abab += -0.50 r2_1_abab(c,b,k,j) t1_bb(d,i) l2_1_abab(k,j,a,d)
        //             += -0.50 r2_1_abab(c,b,j,k) t1_bb(d,i) l2_1_abab(j,k,a,d)
        //             += -0.50 r2_abab(c,b,k,j) t1_bb(d,i) l2_abab(k,j,a,d)
        //             += -0.50 r2_abab(c,b,j,k) t1_bb(d,i) l2_abab(j,k,a,d)
        D_vovv_abab("L,R,a,i,b,c") -= tmps_["38_aabb_LLvvvo"]("L,R,a,c,b,i");
        tmps_["38_aabb_LLvvvo"].~TArrayD();

        // flops: o0v4L2  = o2v4L2 o2v4L2 o0v4L2
        //  mems: o0v4L2  = o0v4L2 o0v4L2 o0v4L2
        tmps_["39_bbbb_LLvvvv"]("L,R,a,b,c,d")  = l2["bbbb"]("L,i,j,a,b") * r2["bbbb"]("R,c,d,i,j");
        tmps_["39_bbbb_LLvvvv"]("L,R,a,b,c,d") += l2_1["bbbb"]("L,i,j,a,b") * r2_1["bbbb"]("R,c,d,i,j");
    }
} // hilbert
#endif