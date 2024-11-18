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


    void EOM_EE_RDM::rdm2_00_2() {

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


        // flops: o2v0L1  = o2v1L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["50_bb_Loo"]("L,i,j")  = l1["bb"]("L,j,a") * t1["bb"]("a,i");

        // D_oovv_abab += +1.00 r2_abab(a,b,i,k) t1_bb(c,j) l1_bb(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += r2["abab"]("R,a,b,i,k") * tmps_["50_bb_Loo"]("L,j,k");

        // D_oovv_abab += +1.00 r0 t2_abab(a,b,i,k) t1_bb(c,j) l1_bb(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,b,i,k") * tmps_["50_bb_Loo"]("L,j,k") * r0("R");

        // D_oovv_bbbb += +1.00 P(i,j) r2_bbbb(b,a,k,i) t1_bb(c,j) l1_bb(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = r2["bbbb"]("R,b,a,k,i") * tmps_["50_bb_Loo"]("L,j,k");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) r0 t2_bbbb(b,a,k,i) t1_bb(c,j) l1_bb(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["bbbb"]("b,a,k,i") * tmps_["50_bb_Loo"]("L,j,k") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v0L2  = o3v2L2 o2v1L2 o3v2L2 o2v0L2 o2v0L2
        //  mems: o2v0L2  = o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2
        tmps_["51_bb_LLoo"]("L,R,i,j")  = 0.50 * l2["bbbb"]("L,l,i,a,c") * r2["bbbb"]("R,a,c,l,j");
        tmps_["51_bb_LLoo"]("L,R,i,j") += l1["bb"]("L,i,a") * r1["bb"]("R,a,j");
        tmps_["51_bb_LLoo"]("L,R,i,j") += l2["abab"]("L,k,i,b,c") * r2["abab"]("R,b,c,k,j");

        // D_oooo_abab += +1.00 d_aa(i,k) r1_bb(a,j) l1_bb(l,a)
        //             += +0.50 d_aa(i,k) r2_abab(a,b,m,j) l2_abab(m,l,a,b)
        //             += +0.50 d_aa(i,k) r2_abab(b,a,m,j) l2_abab(m,l,b,a)
        //             += +0.50 d_aa(i,k) r2_bbbb(a,b,m,j) l2_bbbb(m,l,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += Id["aa_oo"]("i,k") * tmps_["51_bb_LLoo"]("L,R,l,j");

        // D_oovv_abab += +1.00 r1_bb(c,j) t2_abab(a,b,i,k) l1_bb(k,c)
        //             += +0.50 r2_abab(c,d,l,j) t2_abab(a,b,i,k) l2_abab(l,k,c,d)
        //             += +0.50 r2_abab(d,c,l,j) t2_abab(a,b,i,k) l2_abab(l,k,d,c)
        //             += +0.50 r2_bbbb(c,d,l,j) t2_abab(a,b,i,k) l2_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,b,i,k") * tmps_["51_bb_LLoo"]("L,R,k,j");

        // D_oovv_bbbb += -1.00 P(i,j) r1_bb(c,i) t2_bbbb(b,a,k,j) l1_bb(k,c)
        //             += -0.50 P(i,j) r2_abab(c,d,l,i) t2_bbbb(b,a,k,j) l2_abab(l,k,c,d)
        //             += -0.50 P(i,j) r2_abab(d,c,l,i) t2_bbbb(b,a,k,j) l2_abab(l,k,d,c)
        //             += -0.50 P(i,j) r2_bbbb(c,d,l,i) t2_bbbb(b,a,k,j) l2_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["bbbb"]("b,a,k,j") * tmps_["51_bb_LLoo"]("L,R,k,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v0L2  = o1v1L2 o2v1L2
        //  mems: o2v0L2  = o1v1L2 o2v0L2
        tmps_["53_bb_LLoo"]("L,R,i,j")  = (tmps_["40_bb_LLov"]("L,R,i,a") + -1.00 * tmps_["39_bb_LLov"]("L,R,i,a")) * t1["bb"]("a,j");

        // D_oooo_abab += -1.00 d_aa(i,k) r1_bb(b,m) t1_bb(a,j) l2_bbbb(m,l,a,b)
        //             += +1.00 d_aa(i,k) r1_aa(b,m) t1_bb(a,j) l2_abab(m,l,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") -= Id["aa_oo"]("i,k") * tmps_["53_bb_LLoo"]("L,R,l,j");

        // D_oovv_abab += -1.00 r1_bb(d,l) t2_abab(a,b,i,k) t1_bb(c,j) l2_bbbb(l,k,c,d)
        //             += +1.00 r1_aa(d,l) t2_abab(a,b,i,k) t1_bb(c,j) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,b,i,k") * tmps_["53_bb_LLoo"]("L,R,k,j");

        // D_oovv_bbbb += -1.00 P(i,j) r1_bb(d,l) t2_bbbb(b,a,k,i) t1_bb(c,j) l2_bbbb(l,k,c,d)
        //             += +1.00 P(i,j) r1_aa(d,l) t2_bbbb(b,a,k,i) t1_bb(c,j) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["bbbb"]("b,a,k,i") * tmps_["53_bb_LLoo"]("L,R,k,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v0L1  = o2v0L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["52_bb_Loo"]("R,i,j")  = Id["bb_oo"]("i,j") * r0("R");

        // flops: o4v0L2  = o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2
        //  mems: o4v0L2  = o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2
        tmps_["54_bbbb_LLoooo"]("L,R,i,j,k,l")  = 0.50 * tmps_["49_bb_Loo"]("L,l,j") * tmps_["52_bb_Loo"]("R,i,k");
        tmps_["54_bbbb_LLoooo"]("L,R,i,j,k,l") += Id["bb_oo"]("i,j") * tmps_["53_bb_LLoo"]("L,R,k,l");
        tmps_["54_bbbb_LLoooo"]("L,R,i,j,k,l") += tmps_["48_bb_Loo"]("L,l,j") * tmps_["52_bb_Loo"]("R,i,k");
        tmps_["54_bbbb_LLoooo"]("L,R,i,j,k,l") += Id["bb_oo"]("i,k") * tmps_["51_bb_LLoo"]("L,R,j,l");
        tmps_["54_bbbb_LLoooo"]("L,R,i,j,k,l") += tmps_["50_bb_Loo"]("L,l,j") * tmps_["52_bb_Loo"]("R,i,k");

        // D_oooo_bbbb += -1.00 P(i,j) d_bb(j,l) r1_bb(b,m) t1_bb(a,i) l2_bbbb(m,k,a,b)
        //             += +1.00 P(i,j) d_bb(j,l) r1_aa(b,m) t1_bb(a,i) l2_abab(m,k,b,a)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_abab(b,a,m,i) l2_abab(m,l,b,a)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_abab(a,b,m,i) l2_abab(m,l,a,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(a,i) l1_bb(l,a)
        //             += -0.50 P(i,j) d_bb(j,k) r2_abab(a,b,m,i) l2_abab(m,l,a,b)
        //             += -0.50 P(i,j) d_bb(j,k) r2_abab(b,a,m,i) l2_abab(m,l,b,a)
        //             += -0.50 P(i,j) d_bb(j,k) r2_bbbb(a,b,m,i) l2_bbbb(m,l,a,b)
        //             += -1.00 P(i,j) d_bb(j,k) r0 t1_bb(a,i) l1_bb(l,a)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_bbbb(b,a,m,i) l2_bbbb(m,l,b,a)
        D_oooo_bbbb("L,R,i,j,k,l") -= tmps_["54_bbbb_LLoooo"]("L,R,j,l,k,i");
        D_oooo_bbbb("L,R,i,j,k,l") += tmps_["54_bbbb_LLoooo"]("L,R,i,l,k,j");

        // D_oooo_bbbb += +1.00 P(i,j) d_bb(j,k) r1_bb(b,m) t1_bb(a,i) l2_bbbb(m,l,a,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_aa(b,m) t1_bb(a,i) l2_abab(m,l,b,a)
        //             += +0.50 P(i,j) d_bb(j,l) r0 t2_abab(b,a,m,i) l2_abab(m,k,b,a)
        //             += +0.50 P(i,j) d_bb(j,l) r0 t2_abab(a,b,m,i) l2_abab(m,k,a,b)
        //             += +1.00 P(i,j) d_bb(j,l) r1_bb(a,i) l1_bb(k,a)
        //             += +0.50 P(i,j) d_bb(j,l) r2_abab(a,b,m,i) l2_abab(m,k,a,b)
        //             += +0.50 P(i,j) d_bb(j,l) r2_abab(b,a,m,i) l2_abab(m,k,b,a)
        //             += +0.50 P(i,j) d_bb(j,l) r2_bbbb(a,b,m,i) l2_bbbb(m,k,a,b)
        //             += +1.00 P(i,j) d_bb(j,l) r0 t1_bb(a,i) l1_bb(k,a)
        //             += +0.50 P(i,j) d_bb(j,l) r0 t2_bbbb(b,a,m,i) l2_bbbb(m,k,b,a)
        D_oooo_bbbb("L,R,i,j,k,l") += tmps_["54_bbbb_LLoooo"]("L,R,j,k,l,i");
        D_oooo_bbbb("L,R,i,j,k,l") -= tmps_["54_bbbb_LLoooo"]("L,R,i,k,l,j");
        tmps_["54_bbbb_LLoooo"].~TArrayD();

        // flops: o1v1L2  = o2v2L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["55_aa_LLov"]("L,R,i,a")  = l2["aaaa"]("L,j,i,a,b") * r1["aa"]("R,b,j");

        // D_oovv_abab += -1.00 r1_aa(d,l) t1_aa(a,k) t2_abab(c,b,i,j) l2_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("c,b,i,j") * tmps_["55_aa_LLov"]("L,R,k,c") * t1["aa"]("a,k");

        // flops: o1v1L2  = o2v2L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["56_aa_LLov"]("L,R,i,a")  = l2["abab"]("L,i,j,a,b") * r1["bb"]("R,b,j");

        // D_oovv_abab += +1.00 r1_bb(d,l) t1_aa(a,k) t2_abab(c,b,i,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("c,b,i,j") * tmps_["56_aa_LLov"]("L,R,k,c") * t1["aa"]("a,k");

        // flops: o2v0L1  = o3v2L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["57_aa_Loo"]("L,i,j")  = l2["aaaa"]("L,k,j,a,b") * t2["aaaa"]("a,b,k,i");

        // D_oovv_aaaa += -0.50 P(i,j) r0 t2_aaaa(b,a,l,j) t2_aaaa(d,c,k,i) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t2["aaaa"]("b,a,l,j") * tmps_["57_aa_Loo"]("L,i,l") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oooo_abab += +0.50 d_bb(j,l) r0 t2_aaaa(b,a,m,i) l2_aaaa(m,k,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += 0.50 * tmps_["57_aa_Loo"]("L,i,k") * tmps_["52_bb_Loo"]("R,j,l");

        // D_oovv_abab += +0.50 r0 t2_abab(a,b,l,j) t2_aaaa(d,c,k,i) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * t2["abab"]("a,b,l,j") * tmps_["57_aa_Loo"]("L,i,l") * r0("R");

        // flops: o2v0L1  = o3v2L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["58_aa_Loo"]("L,i,j")  = l2["abab"]("L,j,k,a,b") * t2["abab"]("a,b,i,k");

        // D_oooo_abab += +0.50 d_bb(j,l) r0 t2_abab(b,a,i,m) l2_abab(k,m,b,a)
        //             += +0.50 d_bb(j,l) r0 t2_abab(a,b,i,m) l2_abab(k,m,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += tmps_["58_aa_Loo"]("L,i,k") * tmps_["52_bb_Loo"]("R,j,l");

        // D_oovv_aaaa += +0.50 P(i,j) r2_aaaa(b,a,l,i) t2_abab(d,c,j,k) l2_abab(l,k,d,c)
        //             += +0.50 P(i,j) r2_aaaa(b,a,l,i) t2_abab(c,d,j,k) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = r2["aaaa"]("R,b,a,l,i") * tmps_["58_aa_Loo"]("L,j,l");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += -0.50 P(i,j) r0 t2_aaaa(b,a,l,j) t2_abab(d,c,i,k) l2_abab(l,k,d,c)
        //             += -0.50 P(i,j) r0 t2_aaaa(b,a,l,j) t2_abab(c,d,i,k) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["aaaa"]("b,a,l,j") * tmps_["58_aa_Loo"]("L,i,l") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r2_abab(a,b,l,j) t2_abab(d,c,i,k) l2_abab(l,k,d,c)
        //             += +0.50 r2_abab(a,b,l,j) t2_abab(c,d,i,k) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += r2["abab"]("R,a,b,l,j") * tmps_["58_aa_Loo"]("L,i,l");

        // D_oovv_abab += +0.50 r0 t2_abab(a,b,l,j) t2_abab(d,c,i,k) l2_abab(l,k,d,c)
        //             += +0.50 r0 t2_abab(a,b,l,j) t2_abab(c,d,i,k) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,b,l,j") * tmps_["58_aa_Loo"]("L,i,l") * r0("R");

        // flops: o2v0L1  = o2v1L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["59_aa_Loo"]("L,i,j")  = l1["aa"]("L,j,a") * t1["aa"]("a,i");

        // D_oovv_aaaa += +1.00 P(i,j) r2_aaaa(b,a,k,i) t1_aa(c,j) l1_aa(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = r2["aaaa"]("R,b,a,k,i") * tmps_["59_aa_Loo"]("L,j,k");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) r0 t2_aaaa(b,a,k,i) t1_aa(c,j) l1_aa(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["aaaa"]("b,a,k,i") * tmps_["59_aa_Loo"]("L,j,k") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +1.00 r2_abab(a,b,k,j) t1_aa(c,i) l1_aa(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += r2["abab"]("R,a,b,k,j") * tmps_["59_aa_Loo"]("L,i,k");

        // D_oovv_abab += +1.00 r0 t2_abab(a,b,k,j) t1_aa(c,i) l1_aa(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,b,k,j") * tmps_["59_aa_Loo"]("L,i,k") * r0("R");

        // D_oooo_abab += +1.00 d_bb(j,l) r0 t1_aa(a,i) l1_aa(k,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += tmps_["59_aa_Loo"]("L,i,k") * tmps_["52_bb_Loo"]("R,j,l");

        // flops: o2v0L2  = o3v2L2 o2v1L2 o3v2L2 o2v0L2 o2v0L2
        //  mems: o2v0L2  = o2v0L2 o2v0L2 o2v0L2 o2v0L2 o2v0L2
        tmps_["60_aa_LLoo"]("L,R,i,j")  = 0.50 * l2["aaaa"]("L,l,i,a,c") * r2["aaaa"]("R,a,c,l,j");
        tmps_["60_aa_LLoo"]("L,R,i,j") += l1["aa"]("L,i,a") * r1["aa"]("R,a,j");
        tmps_["60_aa_LLoo"]("L,R,i,j") += l2["abab"]("L,i,k,a,b") * r2["abab"]("R,a,b,j,k");

        // D_oooo_abab += +1.00 d_bb(j,l) r1_aa(a,i) l1_aa(k,a)
        //             += +0.50 d_bb(j,l) r2_abab(a,b,i,m) l2_abab(k,m,a,b)
        //             += +0.50 d_bb(j,l) r2_abab(b,a,i,m) l2_abab(k,m,b,a)
        //             += +0.50 d_bb(j,l) r2_aaaa(a,b,m,i) l2_aaaa(m,k,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += Id["bb_oo"]("j,l") * tmps_["60_aa_LLoo"]("L,R,k,i");

        // D_oovv_aaaa += -1.00 P(i,j) r1_aa(c,i) t2_aaaa(b,a,k,j) l1_aa(k,c)
        //             += -0.50 P(i,j) r2_abab(c,d,i,l) t2_aaaa(b,a,k,j) l2_abab(k,l,c,d)
        //             += -0.50 P(i,j) r2_abab(d,c,i,l) t2_aaaa(b,a,k,j) l2_abab(k,l,d,c)
        //             += -0.50 P(i,j) r2_aaaa(c,d,l,i) t2_aaaa(b,a,k,j) l2_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["aaaa"]("b,a,k,j") * tmps_["60_aa_LLoo"]("L,R,k,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +1.00 r1_aa(c,i) t2_abab(a,b,k,j) l1_aa(k,c)
        //             += +0.50 r2_abab(c,d,i,l) t2_abab(a,b,k,j) l2_abab(k,l,c,d)
        //             += +0.50 r2_abab(d,c,i,l) t2_abab(a,b,k,j) l2_abab(k,l,d,c)
        //             += +0.50 r2_aaaa(c,d,l,i) t2_abab(a,b,k,j) l2_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,b,k,j") * tmps_["60_aa_LLoo"]("L,R,k,i");

        // flops: o2v0L1  = o2v0L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["61_aa_Loo"]("R,i,j")  = Id["aa_oo"]("i,j") * r0("R");

        // D_oooo_abab += +1.00 d_aa(i,k) r0 t1_bb(a,j) l1_bb(l,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += tmps_["50_bb_Loo"]("L,j,l") * tmps_["61_aa_Loo"]("R,i,k");

        // D_oooo_abab += +0.50 d_aa(i,k) r0 t2_abab(b,a,m,j) l2_abab(m,l,b,a)
        //             += +0.50 d_aa(i,k) r0 t2_abab(a,b,m,j) l2_abab(m,l,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += tmps_["48_bb_Loo"]("L,j,l") * tmps_["61_aa_Loo"]("R,i,k");

        // D_oooo_abab += +0.50 d_aa(i,k) r0 t2_bbbb(b,a,m,j) l2_bbbb(m,l,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") += 0.50 * tmps_["49_bb_Loo"]("L,j,l") * tmps_["61_aa_Loo"]("R,i,k");

        // flops: o2v0L2  = o1v1L2 o2v1L2
        //  mems: o2v0L2  = o1v1L2 o2v0L2
        tmps_["62_aa_LLoo"]("L,R,i,j")  = (tmps_["55_aa_LLov"]("L,R,i,a") + -1.00 * tmps_["56_aa_LLov"]("L,R,i,a")) * t1["aa"]("a,j");

        // D_oooo_abab += -1.00 d_bb(j,l) r1_aa(b,m) t1_aa(a,i) l2_aaaa(m,k,a,b)
        //             += +1.00 d_bb(j,l) r1_bb(b,m) t1_aa(a,i) l2_abab(k,m,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") -= Id["bb_oo"]("j,l") * tmps_["62_aa_LLoo"]("L,R,k,i");

        // D_oovv_aaaa += -1.00 P(i,j) r1_aa(d,l) t2_aaaa(b,a,k,i) t1_aa(c,j) l2_aaaa(l,k,c,d)
        //             += +1.00 P(i,j) r1_bb(d,l) t2_aaaa(b,a,k,i) t1_aa(c,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["aaaa"]("b,a,k,i") * tmps_["62_aa_LLoo"]("L,R,k,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r1_aa(d,l) t2_abab(a,b,k,j) t1_aa(c,i) l2_aaaa(l,k,c,d)
        //             += +1.00 r1_bb(d,l) t2_abab(a,b,k,j) t1_aa(c,i) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,b,k,j") * tmps_["62_aa_LLoo"]("L,R,k,i");

        // flops: o4v0L2  = o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2
        //  mems: o4v0L2  = o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2 o4v0L2
        tmps_["63_aaaa_LLoooo"]("L,R,i,j,k,l")  = 0.50 * tmps_["57_aa_Loo"]("L,l,j") * tmps_["61_aa_Loo"]("R,i,k");
        tmps_["63_aaaa_LLoooo"]("L,R,i,j,k,l") += Id["aa_oo"]("i,j") * tmps_["62_aa_LLoo"]("L,R,k,l");
        tmps_["63_aaaa_LLoooo"]("L,R,i,j,k,l") += Id["aa_oo"]("i,k") * tmps_["60_aa_LLoo"]("L,R,j,l");
        tmps_["63_aaaa_LLoooo"]("L,R,i,j,k,l") += tmps_["59_aa_Loo"]("L,l,j") * tmps_["61_aa_Loo"]("R,i,k");
        tmps_["63_aaaa_LLoooo"]("L,R,i,j,k,l") += tmps_["58_aa_Loo"]("L,l,j") * tmps_["61_aa_Loo"]("R,i,k");

        // D_oooo_aaaa += +1.00 P(i,j) d_aa(j,k) r1_aa(b,m) t1_aa(a,i) l2_aaaa(m,l,a,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_bb(b,m) t1_aa(a,i) l2_abab(l,m,a,b)
        //             += +1.00 P(i,j) d_aa(j,l) r1_aa(a,i) l1_aa(k,a)
        //             += +0.50 P(i,j) d_aa(j,l) r2_abab(a,b,i,m) l2_abab(k,m,a,b)
        //             += +0.50 P(i,j) d_aa(j,l) r2_abab(b,a,i,m) l2_abab(k,m,b,a)
        //             += +0.50 P(i,j) d_aa(j,l) r2_aaaa(a,b,m,i) l2_aaaa(m,k,a,b)
        //             += +0.50 P(i,j) d_aa(j,l) r0 t2_aaaa(b,a,m,i) l2_aaaa(m,k,b,a)
        //             += +1.00 P(i,j) d_aa(j,l) r0 t1_aa(a,i) l1_aa(k,a)
        //             += +0.50 P(i,j) d_aa(j,l) r0 t2_abab(b,a,i,m) l2_abab(k,m,b,a)
        //             += +0.50 P(i,j) d_aa(j,l) r0 t2_abab(a,b,i,m) l2_abab(k,m,a,b)
        D_oooo_aaaa("L,R,i,j,k,l") += tmps_["63_aaaa_LLoooo"]("L,R,j,k,l,i");
        D_oooo_aaaa("L,R,i,j,k,l") -= tmps_["63_aaaa_LLoooo"]("L,R,i,k,l,j");

        // D_oooo_aaaa += -1.00 P(i,j) d_aa(j,l) r1_aa(b,m) t1_aa(a,i) l2_aaaa(m,k,a,b)
        //             += +1.00 P(i,j) d_aa(j,l) r1_bb(b,m) t1_aa(a,i) l2_abab(k,m,a,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(a,i) l1_aa(l,a)
        //             += -0.50 P(i,j) d_aa(j,k) r2_abab(a,b,i,m) l2_abab(l,m,a,b)
        //             += -0.50 P(i,j) d_aa(j,k) r2_abab(b,a,i,m) l2_abab(l,m,b,a)
        //             += -0.50 P(i,j) d_aa(j,k) r2_aaaa(a,b,m,i) l2_aaaa(m,l,a,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_aaaa(b,a,m,i) l2_aaaa(m,l,b,a)
        //             += -1.00 P(i,j) d_aa(j,k) r0 t1_aa(a,i) l1_aa(l,a)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_abab(b,a,i,m) l2_abab(l,m,b,a)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_abab(a,b,i,m) l2_abab(l,m,a,b)
        D_oooo_aaaa("L,R,i,j,k,l") -= tmps_["63_aaaa_LLoooo"]("L,R,j,l,k,i");
        D_oooo_aaaa("L,R,i,j,k,l") += tmps_["63_aaaa_LLoooo"]("L,R,i,l,k,j");
        tmps_["63_aaaa_LLoooo"].~TArrayD();

        // flops: o1v3L1  = o2v3L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["64_bbbb_Lvovv"]("L,a,i,b,c")  = l2["bbbb"]("L,j,i,b,c") * t1["bb"]("a,j");

        // D_vvvv_bbbb += -1.00 r0 t1_bb(c,i) t1_bb(d,j) l2_bbbb(i,j,a,b)
        // flops: o0v4L2 += o1v4L1 o0v4L2
        //  mems: o0v4L2 += o0v4L1 o0v4L2
        D_vvvv_bbbb("L,R,a,b,c,d") -= t1["bb"]("d,j") * tmps_["64_bbbb_Lvovv"]("L,c,j,a,b") * r0("R");

        // flops: o1v3L2  = o1v3L2
        //  mems: o1v3L2  = o1v3L2
        tmps_["65_bbbb_LLvovv"]("R,L,a,i,b,c")  = r0("R") * tmps_["64_bbbb_Lvovv"]("L,a,i,b,c");
        tmps_["64_bbbb_Lvovv"].~TArrayD();

        // D_vvov_bbbb += +1.00 r0 t1_bb(c,j) l2_bbbb(j,i,a,b)
        D_vvov_bbbb("L,R,a,b,i,c") += tmps_["65_bbbb_LLvovv"]("R,L,c,i,a,b");

        // D_vvvo_bbbb += -1.00 r0 t1_bb(c,j) l2_bbbb(j,i,a,b)
        D_vvvo_bbbb("L,R,a,b,c,i") -= tmps_["65_bbbb_LLvovv"]("R,L,c,i,a,b");
        tmps_["65_bbbb_LLvovv"].~TArrayD();

        // flops: o1v3L1  = o2v3L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["66_aaaa_Lvovv"]("L,a,i,b,c")  = l2["aaaa"]("L,j,i,b,c") * t1["aa"]("a,j");

        // D_vvvv_aaaa += -1.00 r0 t1_aa(c,i) t1_aa(d,j) l2_aaaa(i,j,a,b)
        // flops: o0v4L2 += o1v4L1 o0v4L2
        //  mems: o0v4L2 += o0v4L1 o0v4L2
        D_vvvv_aaaa("L,R,a,b,c,d") -= t1["aa"]("d,j") * tmps_["66_aaaa_Lvovv"]("L,c,j,a,b") * r0("R");

        // flops: o1v3L2  = o1v3L2
        //  mems: o1v3L2  = o1v3L2
        tmps_["67_aaaa_LLvovv"]("R,L,a,i,b,c")  = r0("R") * tmps_["66_aaaa_Lvovv"]("L,a,i,b,c");
        tmps_["66_aaaa_Lvovv"].~TArrayD();

        // D_vvov_aaaa += +1.00 r0 t1_aa(c,j) l2_aaaa(j,i,a,b)
        D_vvov_aaaa("L,R,a,b,i,c") += tmps_["67_aaaa_LLvovv"]("R,L,c,i,a,b");

        // D_vvvo_aaaa += -1.00 r0 t1_aa(c,j) l2_aaaa(j,i,a,b)
        D_vvvo_aaaa("L,R,a,b,c,i") -= tmps_["67_aaaa_LLvovv"]("R,L,c,i,a,b");
        tmps_["67_aaaa_LLvovv"].~TArrayD();

        // flops: o3v1L2  = o3v2L2 o3v2L2 o3v2L1 o4v1L1 o4v1L2 o3v2L2 o3v1L2 o3v1L2 o3v1L2 o4v1L2 o3v1L2
        //  mems: o3v1L2  = o3v1L2 o3v1L2 o3v1L1 o4v0L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        tmps_["68_aaaa_LLovoo"]("L,R,i,a,j,k")  = -1.00 * t2["aaaa"]("a,b,j,k") * tmps_["56_aa_LLov"]("L,R,i,b");
        tmps_["68_aaaa_LLovoo"]("L,R,i,a,j,k") -= l1["aa"]("L,i,b") * r2["aaaa"]("R,a,b,j,k");
        tmps_["68_aaaa_LLovoo"]("L,R,i,a,j,k") -= l2["aaaa"]("L,l,i,c,b") * t1["aa"]("b,j") * t1["aa"]("c,k") * r1["aa"]("R,a,l");
        tmps_["68_aaaa_LLovoo"]("L,R,i,a,j,k") += t2["aaaa"]("a,b,j,k") * tmps_["55_aa_LLov"]("L,R,i,b");
        tmps_["68_aaaa_LLovoo"]("L,R,i,a,j,k") += 0.50 * r1["aa"]("R,a,l") * tmps_["69_aaaa_Loooo"]("L,j,k,l,i");

        // D_ooov_aaaa += +1.00 r1_aa(c,l) t2_aaaa(a,b,i,j) l2_aaaa(l,k,b,c)
        //             += -1.00 r1_aa(a,l) t1_aa(b,i) t1_aa(c,j) l2_aaaa(l,k,c,b)
        //             += -1.00 r2_aaaa(a,b,i,j) l1_aa(k,b)
        //             += -1.00 r1_bb(c,l) t2_aaaa(a,b,i,j) l2_abab(k,l,b,c)
        //             += +0.50 r1_aa(a,l) t2_aaaa(c,b,i,j) l2_aaaa(l,k,c,b)
        D_ooov_aaaa("L,R,i,j,k,a") += tmps_["68_aaaa_LLovoo"]("L,R,k,a,i,j");

        // D_oovo_aaaa += -1.00 r1_aa(c,l) t2_aaaa(a,b,i,j) l2_aaaa(l,k,b,c)
        //             += +1.00 r1_aa(a,l) t1_aa(b,i) t1_aa(c,j) l2_aaaa(l,k,c,b)
        //             += +1.00 r2_aaaa(a,b,i,j) l1_aa(k,b)
        //             += +1.00 r1_bb(c,l) t2_aaaa(a,b,i,j) l2_abab(k,l,b,c)
        //             += -0.50 r1_aa(a,l) t2_aaaa(c,b,i,j) l2_aaaa(l,k,c,b)
        D_oovo_aaaa("L,R,i,j,a,k") -= tmps_["68_aaaa_LLovoo"]("L,R,k,a,i,j");

        // D_oovv_aaaa += -1.00 P(a,b) r1_aa(d,l) t2_aaaa(a,c,i,j) t1_aa(b,k) l2_aaaa(l,k,c,d)
        //             += +1.00 P(a,b) r1_aa(a,l) t1_aa(b,k) t1_aa(c,i) t1_aa(d,j) l2_aaaa(l,k,d,c)
        //             += +1.00 P(a,b) r2_aaaa(a,c,i,j) t1_aa(b,k) l1_aa(k,c)
        //             += +1.00 P(a,b) r1_bb(d,l) t2_aaaa(a,c,i,j) t1_aa(b,k) l2_abab(k,l,c,d)
        //             += -0.50 P(a,b) r1_aa(a,l) t1_aa(b,k) t2_aaaa(d,c,i,j) l2_aaaa(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("b,k") * tmps_["68_aaaa_LLovoo"]("L,R,k,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();
        tmps_["68_aaaa_LLovoo"].~TArrayD();

        // flops: o4v0L1  = o4v2L1
        //  mems: o4v0L1  = o4v0L1
        tmps_["69_aaaa_Loooo"]("L,i,j,k,l")  = l2["aaaa"]("L,k,l,a,b") * t2["aaaa"]("a,b,i,j");

        // D_oooo_aaaa += +0.50 r0 t2_aaaa(b,a,i,j) l2_aaaa(l,k,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_aaaa("L,R,i,j,k,l") += 0.50 * r0("R") * tmps_["69_aaaa_Loooo"]("L,i,j,l,k");

        // D_oovv_aaaa += +0.25 r2_aaaa(b,a,l,k) t2_aaaa(d,c,i,j) l2_aaaa(l,k,d,c)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") += 0.25 * r2["aaaa"]("R,b,a,l,k") * tmps_["69_aaaa_Loooo"]("L,i,j,l,k");

        // D_oovv_aaaa += +0.25 r0 t2_aaaa(b,a,k,l) t2_aaaa(d,c,i,j) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o4v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") += 0.25 * t2["aaaa"]("b,a,k,l") * tmps_["69_aaaa_Loooo"]("L,i,j,k,l") * r0("R");

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["70_aaaa_Looov"]("L,i,j,k,a")  = l2["aaaa"]("L,j,k,b,a") * t1["aa"]("b,i");

        // D_oooo_aaaa += -1.00 P(i,j) r1_aa(b,i) t1_aa(a,j) l2_aaaa(l,k,a,b)
        // flops: o4v0L2 += o4v1L2
        //  mems: o4v0L2 += o4v0L2
        tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l")  = r1["aa"]("R,b,i") * tmps_["70_aaaa_Looov"]("L,j,l,k,b");
        D_oooo_aaaa("L,R,i,j,k,l") -= tmps_["perm_aaaa_LLoooo"]("L,R,i,j,k,l");
        D_oooo_aaaa("L,R,i,j,k,l") += tmps_["perm_aaaa_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_aaaa_LLoooo"].~TArrayD();

        // D_oooo_aaaa += -1.00 r0 t1_aa(a,i) t1_aa(b,j) l2_aaaa(l,k,b,a)
        // flops: o4v0L2 += o4v1L1 o4v0L2
        //  mems: o4v0L2 += o4v0L1 o4v0L2
        D_oooo_aaaa("L,R,i,j,k,l") -= t1["aa"]("a,i") * tmps_["70_aaaa_Looov"]("L,j,l,k,a") * r0("R");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r1_aa(a,l) t2_aaaa(b,c,k,i) t1_aa(d,j) l2_aaaa(l,k,d,c)
        // flops: o2v2L2 += o4v2L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["aaaa"]("b,c,k,i") * tmps_["70_aaaa_Looov"]("L,j,l,k,c") * r1["aa"]("R,a,l");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +1.00 r1_aa(a,l) t2_abab(c,b,k,j) t1_aa(d,i) l2_aaaa(l,k,d,c)
        // flops: o2v2L2 += o4v2L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("c,b,k,j") * tmps_["70_aaaa_Looov"]("L,i,l,k,c") * r1["aa"]("R,a,l");

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["71_aaaa_Lvooo"]("L,a,i,j,k")  = l1["aa"]("L,k,b") * t2["aaaa"]("a,b,i,j");

        // D_oovv_aaaa += -1.00 P(a,b) r1_aa(a,k) t2_aaaa(b,c,i,j) l1_aa(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = r1["aa"]("R,a,k") * tmps_["71_aaaa_Lvooo"]("L,b,i,j,k");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(a,b) r0 t2_aaaa(a,c,i,j) t1_aa(b,k) l1_aa(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("b,k") * tmps_["71_aaaa_Lvooo"]("L,a,i,j,k") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o4v1L1 o4v0L1 o4v1L1
        //  mems: o3v1L1  = o4v0L1 o4v0L1 o3v1L1
        tmps_["72_aaaa_Looov"]("L,i,j,k,a")  = t1["aa"]("a,l") * (tmps_["69_aaaa_Loooo"]("L,i,j,l,k") + -2.00 * t1["aa"]("c,i") * tmps_["70_aaaa_Looov"]("L,j,l,k,c"));
        tmps_["70_aaaa_Looov"].~TArrayD();
        tmps_["69_aaaa_Loooo"].~TArrayD();

        // D_oovv_aaaa += -0.50 r0 t1_aa(a,k) t1_aa(b,l) t2_aaaa(d,c,i,j) l2_aaaa(k,l,d,c)
        //             += +1.00 r0 t1_aa(a,k) t1_aa(b,l) t1_aa(c,i) t1_aa(d,j) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") -= 0.50 * t1["aa"]("b,l") * tmps_["72_aaaa_Looov"]("L,i,j,l,a") * r0("R");

        // flops: o3v1L2  = o3v1L1 o3v1L2
        //  mems: o3v1L2  = o3v1L1 o3v1L2
        tmps_["73_aaaa_LLooov"]("L,R,i,j,k,a")  = (tmps_["72_aaaa_Looov"]("L,i,j,k,a") + -2.00 * tmps_["71_aaaa_Lvooo"]("L,a,i,j,k")) * r0("R");
        tmps_["72_aaaa_Looov"].~TArrayD();
        tmps_["71_aaaa_Lvooo"].~TArrayD();

        // D_ooov_aaaa += +0.50 r0 t1_aa(a,l) t2_aaaa(c,b,i,j) l2_aaaa(l,k,c,b)
        //             += -1.00 r0 t1_aa(a,l) t1_aa(b,i) t1_aa(c,j) l2_aaaa(l,k,c,b)
        //             += -1.00 r0 t2_aaaa(a,b,i,j) l1_aa(k,b)
        D_ooov_aaaa("L,R,i,j,k,a") += 0.50 * tmps_["73_aaaa_LLooov"]("L,R,i,j,k,a");

        // D_oovo_aaaa += -0.50 r0 t1_aa(a,l) t2_aaaa(c,b,i,j) l2_aaaa(l,k,c,b)
        //             += +1.00 r0 t1_aa(a,l) t1_aa(b,i) t1_aa(c,j) l2_aaaa(l,k,c,b)
        //             += +1.00 r0 t2_aaaa(a,b,i,j) l1_aa(k,b)
        D_oovo_aaaa("L,R,i,j,a,k") -= 0.50 * tmps_["73_aaaa_LLooov"]("L,R,i,j,k,a");
        tmps_["73_aaaa_LLooov"].~TArrayD();

        // flops: o0v4L1  = o2v4L1
        //  mems: o0v4L1  = o0v4L1
        tmps_["74_aaaa_Lvvvv"]("L,a,b,c,d")  = l2["aaaa"]("L,i,j,a,b") * t2["aaaa"]("c,d,i,j");

        // D_vvvv_aaaa += +0.50 r0 t2_aaaa(d,c,i,j) l2_aaaa(i,j,a,b)
        // flops: o0v4L2 += o0v4L2
        //  mems: o0v4L2 += o0v4L2
        D_vvvv_aaaa("L,R,a,b,c,d") += 0.50 * r0("R") * tmps_["74_aaaa_Lvvvv"]("L,a,b,d,c");

        // flops: o2v2L2  = o3v2L2 o3v2L2
        //  mems: o2v2L2  = o3v1L2 o2v2L2
        tmps_["75_aaaa_LLovov"]("L,R,i,a,j,b")  = l2["aaaa"]("L,k,i,a,c") * r1["aa"]("R,c,j") * t1["aa"]("b,k");

        // D_ovov_aaaa += +1.00 r1_aa(c,i) t1_aa(b,k) l2_aaaa(k,j,a,c)
        D_ovov_aaaa("L,R,i,a,j,b") += tmps_["75_aaaa_LLovov"]("L,R,j,a,i,b");

        // D_ovvo_aaaa += -1.00 r1_aa(c,i) t1_aa(b,k) l2_aaaa(k,j,a,c)
        D_ovvo_aaaa("L,R,i,a,b,j") -= tmps_["75_aaaa_LLovov"]("L,R,j,a,i,b");

        // flops: o1v3L1  = o1v4L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["76_aaaa_Lvvvo"]("L,a,b,c,i")  = t1["aa"]("d,i") * tmps_["74_aaaa_Lvvvv"]("L,a,d,b,c");
        tmps_["74_aaaa_Lvvvv"].~TArrayD();

        // D_oovv_aaaa += -0.50 r0 t2_aaaa(b,a,k,l) t1_aa(c,i) t1_aa(d,j) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") -= 0.50 * t1["aa"]("d,j") * tmps_["76_aaaa_Lvvvo"]("L,d,b,a,i") * r0("R");

        // flops: o1v3L2  = o3v2L1 o3v2L1 o2v3L1 o1v3L2 o2v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o1v3L2 o1v3L2
        //  mems: o1v3L2  = o3v1L1 o2v2L1 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        tmps_["77_aaaa_LLvvvo"]("R,L,a,b,c,i")  = -2.00 * l2["aaaa"]("L,j,k,a,d") * t1["aa"]("d,i") * t1["aa"]("c,j") * t1["aa"]("b,k") * r0("R");
        tmps_["77_aaaa_LLvvvo"]("R,L,a,b,c,i") += 2.00 * tmps_["55_aa_LLov"]("L,R,j,a") * t2["aaaa"]("b,c,j,i");
        tmps_["77_aaaa_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * tmps_["75_aaaa_LLovov"]("L,R,k,a,i,c") * t1["aa"]("b,k");
        tmps_["77_aaaa_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * tmps_["56_aa_LLov"]("L,R,j,a") * t2["aaaa"]("b,c,j,i");
        tmps_["77_aaaa_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l1["aa"]("L,j,a") * t2["aaaa"]("b,c,j,i") * r0("R");
        tmps_["77_aaaa_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l1["aa"]("L,j,a") * r2["aaaa"]("R,b,c,j,i");
        tmps_["77_aaaa_LLvvvo"]("R,L,a,b,c,i") += r0("R") * tmps_["76_aaaa_Lvvvo"]("L,a,b,c,i");
        tmps_["76_aaaa_Lvvvo"].~TArrayD();
        tmps_["75_aaaa_LLovov"].~TArrayD();

        // D_ovvv_aaaa += +0.50 r0 t2_aaaa(c,b,j,k) t1_aa(d,i) l2_aaaa(j,k,a,d)
        //             += -1.00 r0 t1_aa(b,j) t1_aa(c,k) t1_aa(d,i) l2_aaaa(j,k,a,d)
        //             += +1.00 r1_aa(d,k) t2_aaaa(c,b,j,i) l2_aaaa(k,j,a,d)
        //             += -1.00 r1_aa(d,i) t1_aa(b,j) t1_aa(c,k) l2_aaaa(j,k,a,d)
        //             += -1.00 r1_bb(d,k) t2_aaaa(c,b,j,i) l2_abab(j,k,a,d)
        //             += -1.00 r0 t2_aaaa(c,b,j,i) l1_aa(j,a)
        //             += -1.00 r2_aaaa(c,b,j,i) l1_aa(j,a)
        D_ovvv_aaaa("L,R,i,a,b,c") += 0.50 * tmps_["77_aaaa_LLvvvo"]("R,L,a,c,b,i");

        // D_vovv_aaaa += -0.50 r0 t2_aaaa(c,b,j,k) t1_aa(d,i) l2_aaaa(j,k,a,d)
        //             += +1.00 r0 t1_aa(b,j) t1_aa(c,k) t1_aa(d,i) l2_aaaa(j,k,a,d)
        //             += -1.00 r1_aa(d,k) t2_aaaa(c,b,j,i) l2_aaaa(k,j,a,d)
        //             += +1.00 r1_aa(d,i) t1_aa(b,j) t1_aa(c,k) l2_aaaa(j,k,a,d)
        //             += +1.00 r1_bb(d,k) t2_aaaa(c,b,j,i) l2_abab(j,k,a,d)
        //             += +1.00 r0 t2_aaaa(c,b,j,i) l1_aa(j,a)
        //             += +1.00 r2_aaaa(c,b,j,i) l1_aa(j,a)
        D_vovv_aaaa("L,R,a,i,b,c") -= 0.50 * tmps_["77_aaaa_LLvvvo"]("R,L,a,c,b,i");
        tmps_["77_aaaa_LLvvvo"].~TArrayD();

        // flops: o4v0L1  = o4v2L1
        //  mems: o4v0L1  = o4v0L1
        tmps_["78_bbbb_Loooo"]("L,i,j,k,l")  = l2["bbbb"]("L,k,l,a,b") * t2["bbbb"]("a,b,i,j");

        // D_oooo_bbbb += +0.50 r0 t2_bbbb(b,a,i,j) l2_bbbb(l,k,b,a)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_bbbb("L,R,i,j,k,l") += 0.50 * r0("R") * tmps_["78_bbbb_Loooo"]("L,i,j,l,k");

        // D_oovv_bbbb += +0.25 r2_bbbb(b,a,l,k) t2_bbbb(d,c,i,j) l2_bbbb(l,k,d,c)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += 0.25 * r2["bbbb"]("R,b,a,l,k") * tmps_["78_bbbb_Loooo"]("L,i,j,l,k");

        // D_oovv_bbbb += +0.25 r0 t2_bbbb(b,a,k,l) t2_bbbb(d,c,i,j) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o4v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += 0.25 * t2["bbbb"]("b,a,k,l") * tmps_["78_bbbb_Loooo"]("L,i,j,k,l") * r0("R");

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["79_bbbb_Loovo"]("L,i,j,a,k")  = l2["bbbb"]("L,i,j,b,a") * t1["bb"]("b,k");

        // D_oooo_bbbb += -1.00 P(i,j) r1_bb(b,i) t1_bb(a,j) l2_bbbb(l,k,a,b)
        // flops: o4v0L2 += o4v1L2
        //  mems: o4v0L2 += o4v0L2
        tmps_["perm_bbbb_LLoooo"]("L,R,i,j,k,l")  = r1["bb"]("R,b,i") * tmps_["79_bbbb_Loovo"]("L,l,k,b,j");
        D_oooo_bbbb("L,R,i,j,k,l") -= tmps_["perm_bbbb_LLoooo"]("L,R,i,j,k,l");
        D_oooo_bbbb("L,R,i,j,k,l") += tmps_["perm_bbbb_LLoooo"]("L,R,j,i,k,l");
        tmps_["perm_bbbb_LLoooo"].~TArrayD();

        // D_oooo_bbbb += -1.00 r0 t1_bb(a,i) t1_bb(b,j) l2_bbbb(l,k,b,a)
        // flops: o4v0L2 += o4v1L1 o4v0L2
        //  mems: o4v0L2 += o4v0L1 o4v0L2
        D_oooo_bbbb("L,R,i,j,k,l") -= t1["bb"]("a,i") * tmps_["79_bbbb_Loovo"]("L,l,k,a,j") * r0("R");

        // D_oovv_abab += +1.00 r1_bb(b,l) t2_abab(a,c,i,k) t1_bb(d,j) l2_bbbb(l,k,d,c)
        // flops: o2v2L2 += o4v2L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,c,i,k") * tmps_["79_bbbb_Loovo"]("L,l,k,c,j") * r1["bb"]("R,b,l");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r1_bb(a,l) t2_bbbb(b,c,k,i) t1_bb(d,j) l2_bbbb(l,k,d,c)
        // flops: o2v2L2 += o4v2L1 o3v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["bbbb"]("b,c,k,i") * tmps_["79_bbbb_Loovo"]("L,l,k,c,j") * r1["bb"]("R,a,l");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["80_bbbb_Lvooo"]("L,a,i,j,k")  = l1["bb"]("L,k,b") * t2["bbbb"]("a,b,i,j");

        // D_oovv_bbbb += -1.00 P(a,b) r1_bb(a,k) t2_bbbb(b,c,i,j) l1_bb(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = r1["bb"]("R,a,k") * tmps_["80_bbbb_Lvooo"]("L,b,i,j,k");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(a,b) r0 t2_bbbb(a,c,i,j) t1_bb(b,k) l1_bb(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("b,k") * tmps_["80_bbbb_Lvooo"]("L,a,i,j,k") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o4v1L1 o4v0L1 o4v1L1
        //  mems: o3v1L1  = o4v0L1 o4v0L1 o3v1L1
        tmps_["81_bbbb_Looov"]("L,i,j,k,a")  = t1["bb"]("a,l") * (tmps_["78_bbbb_Loooo"]("L,i,j,l,k") + -2.00 * t1["bb"]("c,i") * tmps_["79_bbbb_Loovo"]("L,l,k,c,j"));
        tmps_["79_bbbb_Loovo"].~TArrayD();

        // D_oovv_bbbb += -0.50 r0 t1_bb(a,k) t1_bb(b,l) t2_bbbb(d,c,i,j) l2_bbbb(k,l,d,c)
        //             += +1.00 r0 t1_bb(a,k) t1_bb(b,l) t1_bb(c,i) t1_bb(d,j) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") -= 0.50 * t1["bb"]("b,l") * tmps_["81_bbbb_Looov"]("L,i,j,l,a") * r0("R");

        // flops: o3v1L2  = o1v1L2 o3v2L2 o3v2L2 o3v1L2 o3v2L1 o4v1L1 o4v1L2 o3v1L2 o4v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        //  mems: o3v1L2  = o1v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L1 o4v0L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        tmps_["82_bbbb_LLvooo"]("R,L,a,i,j,k")  = -1.00 * t2["bbbb"]("a,b,i,j") * (tmps_["40_bb_LLov"]("L,R,k,b") + -1.00 * tmps_["39_bb_LLov"]("L,R,k,b"));
        tmps_["82_bbbb_LLvooo"]("R,L,a,i,j,k") += l1["bb"]("L,k,b") * r2["bbbb"]("R,a,b,i,j");
        tmps_["82_bbbb_LLvooo"]("R,L,a,i,j,k") += t1["bb"]("b,i") * l2["bbbb"]("L,l,k,c,b") * t1["bb"]("c,j") * r1["bb"]("R,a,l");
        tmps_["82_bbbb_LLvooo"]("R,L,a,i,j,k") -= 0.50 * tmps_["78_bbbb_Loooo"]("L,i,j,l,k") * r1["bb"]("R,a,l");
        tmps_["82_bbbb_LLvooo"]("R,L,a,i,j,k") -= 0.50 * tmps_["81_bbbb_Looov"]("L,i,j,k,a") * r0("R");
        tmps_["82_bbbb_LLvooo"]("R,L,a,i,j,k") += r0("R") * tmps_["80_bbbb_Lvooo"]("L,a,i,j,k");
        tmps_["81_bbbb_Looov"].~TArrayD();
        tmps_["80_bbbb_Lvooo"].~TArrayD();
        tmps_["78_bbbb_Loooo"].~TArrayD();

        // D_ooov_bbbb += -1.00 r0 t2_bbbb(a,b,i,j) l1_bb(k,b)
        //             += +1.00 r1_bb(c,l) t2_bbbb(a,b,i,j) l2_bbbb(l,k,b,c)
        //             += -1.00 r1_aa(c,l) t2_bbbb(a,b,i,j) l2_abab(l,k,c,b)
        //             += -1.00 r2_bbbb(a,b,i,j) l1_bb(k,b)
        //             += -1.00 r1_bb(a,l) t1_bb(b,i) t1_bb(c,j) l2_bbbb(l,k,c,b)
        //             += +0.50 r1_bb(a,l) t2_bbbb(c,b,i,j) l2_bbbb(l,k,c,b)
        //             += +0.50 r0 t1_bb(a,l) t2_bbbb(c,b,i,j) l2_bbbb(l,k,c,b)
        //             += -1.00 r0 t1_bb(a,l) t1_bb(b,i) t1_bb(c,j) l2_bbbb(l,k,c,b)
        D_ooov_bbbb("L,R,i,j,k,a") -= tmps_["82_bbbb_LLvooo"]("R,L,a,i,j,k");

        // D_oovo_bbbb += +1.00 r0 t2_bbbb(a,b,i,j) l1_bb(k,b)
        //             += -1.00 r1_bb(c,l) t2_bbbb(a,b,i,j) l2_bbbb(l,k,b,c)
        //             += +1.00 r1_aa(c,l) t2_bbbb(a,b,i,j) l2_abab(l,k,c,b)
        //             += +1.00 r2_bbbb(a,b,i,j) l1_bb(k,b)
        //             += +1.00 r1_bb(a,l) t1_bb(b,i) t1_bb(c,j) l2_bbbb(l,k,c,b)
        //             += -0.50 r1_bb(a,l) t2_bbbb(c,b,i,j) l2_bbbb(l,k,c,b)
        //             += -0.50 r0 t1_bb(a,l) t2_bbbb(c,b,i,j) l2_bbbb(l,k,c,b)
        //             += +1.00 r0 t1_bb(a,l) t1_bb(b,i) t1_bb(c,j) l2_bbbb(l,k,c,b)
        D_oovo_bbbb("L,R,i,j,a,k") += tmps_["82_bbbb_LLvooo"]("R,L,a,i,j,k");
        tmps_["82_bbbb_LLvooo"].~TArrayD();

        // flops: o2v2L2  = o3v2L2 o3v2L2
        //  mems: o2v2L2  = o3v1L2 o2v2L2
        tmps_["83_bbbb_LLovov"]("L,R,i,a,j,b")  = l2["bbbb"]("L,k,i,a,c") * r1["bb"]("R,c,j") * t1["bb"]("b,k");

        // D_ovov_bbbb += +1.00 r1_bb(c,i) t1_bb(b,k) l2_bbbb(k,j,a,c)
        D_ovov_bbbb("L,R,i,a,j,b") += tmps_["83_bbbb_LLovov"]("L,R,j,a,i,b");

        // D_ovvo_bbbb += -1.00 r1_bb(c,i) t1_bb(b,k) l2_bbbb(k,j,a,c)
        D_ovvo_bbbb("L,R,i,a,b,j") -= tmps_["83_bbbb_LLovov"]("L,R,j,a,i,b");

        // flops: o1v3L1  = o1v4L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["84_bbbb_Lvvvo"]("L,a,b,c,i")  = t1["bb"]("d,i") * tmps_["4_bbbb_Lvvvv"]("L,a,d,b,c");
        tmps_["4_bbbb_Lvvvv"].~TArrayD();

        // D_oovv_bbbb += -0.50 r0 t2_bbbb(b,a,k,l) t1_bb(c,i) t1_bb(d,j) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") -= 0.50 * t1["bb"]("d,j") * tmps_["84_bbbb_Lvvvo"]("L,d,b,a,i") * r0("R");

        // flops: o1v3L2  = o2v3L2 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L1 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        //  mems: o1v3L2  = o1v3L2 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        tmps_["85_bbbb_LLvvvo"]("R,L,a,b,c,i")  = -2.00 * l1["bb"]("L,j,a") * r2["bbbb"]("R,b,c,j,i");
        tmps_["85_bbbb_LLvvvo"]("R,L,a,b,c,i") += 2.00 * tmps_["40_bb_LLov"]("L,R,j,a") * t2["bbbb"]("b,c,j,i");
        tmps_["85_bbbb_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l2["bbbb"]("L,j,k,a,d") * t1["bb"]("d,i") * t1["bb"]("c,j") * t1["bb"]("b,k") * r0("R");
        tmps_["85_bbbb_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * tmps_["39_bb_LLov"]("L,R,j,a") * t2["bbbb"]("b,c,j,i");
        tmps_["85_bbbb_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * tmps_["83_bbbb_LLovov"]("L,R,k,a,i,c") * t1["bb"]("b,k");
        tmps_["85_bbbb_LLvvvo"]("R,L,a,b,c,i") -= 2.00 * l1["bb"]("L,j,a") * t2["bbbb"]("b,c,j,i") * r0("R");
        tmps_["85_bbbb_LLvvvo"]("R,L,a,b,c,i") += r0("R") * tmps_["84_bbbb_Lvvvo"]("L,a,b,c,i");
        tmps_["84_bbbb_Lvvvo"].~TArrayD();
        tmps_["83_bbbb_LLovov"].~TArrayD();

        // D_ovvv_bbbb += +0.50 r0 t2_bbbb(c,b,j,k) t1_bb(d,i) l2_bbbb(j,k,a,d)
        //             += -1.00 r2_bbbb(c,b,j,i) l1_bb(j,a)
        //             += +1.00 r1_bb(d,k) t2_bbbb(c,b,j,i) l2_bbbb(k,j,a,d)
        //             += -1.00 r0 t1_bb(b,j) t1_bb(c,k) t1_bb(d,i) l2_bbbb(j,k,a,d)
        //             += -1.00 r1_aa(d,k) t2_bbbb(c,b,j,i) l2_abab(k,j,d,a)
        //             += -1.00 r1_bb(d,i) t1_bb(b,j) t1_bb(c,k) l2_bbbb(j,k,a,d)
        //             += -1.00 r0 t2_bbbb(c,b,j,i) l1_bb(j,a)
        D_ovvv_bbbb("L,R,i,a,b,c") += 0.50 * tmps_["85_bbbb_LLvvvo"]("R,L,a,c,b,i");

        // D_vovv_bbbb += -0.50 r0 t2_bbbb(c,b,j,k) t1_bb(d,i) l2_bbbb(j,k,a,d)
        //             += +1.00 r2_bbbb(c,b,j,i) l1_bb(j,a)
        //             += -1.00 r1_bb(d,k) t2_bbbb(c,b,j,i) l2_bbbb(k,j,a,d)
        //             += +1.00 r0 t1_bb(b,j) t1_bb(c,k) t1_bb(d,i) l2_bbbb(j,k,a,d)
        //             += +1.00 r1_aa(d,k) t2_bbbb(c,b,j,i) l2_abab(k,j,d,a)
        //             += +1.00 r1_bb(d,i) t1_bb(b,j) t1_bb(c,k) l2_bbbb(j,k,a,d)
        //             += +1.00 r0 t2_bbbb(c,b,j,i) l1_bb(j,a)
        D_vovv_bbbb("L,R,a,i,b,c") -= 0.50 * tmps_["85_bbbb_LLvvvo"]("R,L,a,c,b,i");
        tmps_["85_bbbb_LLvvvo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["86_bbbb_Lvoov"]("L,a,i,j,b")  = l2["bbbb"]("L,k,j,b,c") * t2["bbbb"]("a,c,k,i");

        // D_oovv_bbbb += +1.00 P(i,j) r0 t2_bbbb(a,c,k,i) t2_bbbb(b,d,l,j) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["bbbb"]("b,d,l,j") * tmps_["86_bbbb_Lvoov"]("L,a,i,l,d") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["87_bb_Lvv"]("L,a,b")  = l2["bbbb"]("L,i,j,a,c") * t2["bbbb"]("b,c,i,j");

        // D_oovv_bbbb += -0.50 r0 t2_bbbb(a,c,k,l) t2_bbbb(b,d,i,j) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") -= 0.50 * t2["bbbb"]("b,d,i,j") * tmps_["87_bb_Lvv"]("L,d,a") * r0("R");

        // flops: o0v2L2  = o2v3L2 o2v3L2 o0v2L2
        //  mems: o0v2L2  = o0v2L2 o0v2L2 o0v2L2
        tmps_["91_bb_LLvv"]("L,R,a,b")  = 0.50 * l2["bbbb"]("L,k,j,a,d") * r2["bbbb"]("R,b,d,k,j");
        tmps_["91_bb_LLvv"]("L,R,a,b") += l2["abab"]("L,i,j,c,a") * r2["abab"]("R,c,b,i,j");

        // D_oovv_abab += +0.50 r2_abab(d,b,l,k) t2_abab(a,c,i,j) l2_abab(l,k,d,c)
        //             += +0.50 r2_abab(d,b,k,l) t2_abab(a,c,i,j) l2_abab(k,l,d,c)
        //             += +0.50 r2_bbbb(b,d,l,k) t2_abab(a,c,i,j) l2_bbbb(l,k,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,c,i,j") * tmps_["91_bb_LLvv"]("L,R,c,b");

        // D_oovv_bbbb += -0.50 P(a,b) r2_abab(d,a,l,k) t2_bbbb(b,c,i,j) l2_abab(l,k,d,c)
        //             += -0.50 P(a,b) r2_abab(d,a,k,l) t2_bbbb(b,c,i,j) l2_abab(k,l,d,c)
        //             += -0.50 P(a,b) r2_bbbb(a,d,l,k) t2_bbbb(b,c,i,j) l2_bbbb(l,k,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["bbbb"]("b,c,i,j") * tmps_["91_bb_LLvv"]("L,R,c,a");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o0v2L2  = o1v1L2 o1v2L2
        //  mems: o0v2L2  = o1v1L2 o0v2L2
        tmps_["92_bb_LLvv"]("L,R,a,b")  = (tmps_["39_bb_LLov"]("L,R,i,a") + -1.00 * tmps_["40_bb_LLov"]("L,R,i,a")) * t1["bb"]("b,i");

        // flops: o1v1L1  = o1v1L1
        //  mems: o1v1L1  = o1v1L1
        tmps_["90_bb_Lvo"]("R,a,i")  = r0("R") * t1["bb"]("a,i");

        // flops: o0v2L1  = o1v2L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["89_bb_Lvv"]("L,a,b")  = l1["bb"]("L,i,a") * t1["bb"]("b,i");

        // flops: o0v2L2  = o1v2L2
        //  mems: o0v2L2  = o0v2L2
        tmps_["88_bb_LLvv"]("L,R,a,b")  = l1["bb"]("L,i,a") * r1["bb"]("R,b,i");

        // flops: o1v3L2  = o0v2L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o3v3L1 o2v3L2 o3v2L1 o3v2L1 o2v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        //  mems: o1v3L2  = o0v2L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o3v1L1 o2v2L1 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        tmps_["93_bbbb_LLvvvo"]("L,R,a,b,c,i")  = (tmps_["89_bb_Lvv"]("L,a,b") + tmps_["30_bb_Lvv"]("L,a,b")) * r1["bb"]("R,c,i");
        tmps_["93_bbbb_LLvvvo"]("L,R,a,b,c,i") -= t1["bb"]("b,i") * tmps_["88_bb_LLvv"]("L,R,a,c");
        tmps_["93_bbbb_LLvvvo"]("L,R,a,b,c,i") += 0.50 * r1["bb"]("R,c,i") * tmps_["87_bb_Lvv"]("L,a,b");
        tmps_["93_bbbb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["91_bb_LLvv"]("L,R,a,c") * t1["bb"]("b,i");
        tmps_["93_bbbb_LLvvvo"]("L,R,a,b,c,i") -= 0.50 * tmps_["87_bb_Lvv"]("L,a,c") * tmps_["90_bb_Lvo"]("R,b,i");
        tmps_["93_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["24_bbbb_Lvoov"]("L,c,i,k,a") * t1["bb"]("b,k") * r0("R");
        tmps_["93_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["20_bbbb_LLovvo"]("L,R,j,a,c,i") * t1["bb"]("b,j");
        tmps_["93_bbbb_LLvvvo"]("L,R,a,b,c,i") += l2["bbbb"]("L,k,j,a,d") * t2["bbbb"]("b,d,j,i") * r1["bb"]("R,c,k");
        tmps_["93_bbbb_LLvvvo"]("L,R,a,b,c,i") += l2["bbbb"]("L,k,j,a,d") * t1["bb"]("d,i") * t1["bb"]("b,j") * r1["bb"]("R,c,k");
        tmps_["93_bbbb_LLvvvo"]("L,R,a,b,c,i") -= l2["abab"]("L,l,k,e,a") * t2["abab"]("e,b,l,i") * r1["bb"]("R,c,k");
        tmps_["93_bbbb_LLvvvo"]("L,R,a,b,c,i") += tmps_["86_bbbb_Lvoov"]("L,c,i,k,a") * t1["bb"]("b,k") * r0("R");
        tmps_["93_bbbb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["89_bb_Lvv"]("L,a,c") * tmps_["90_bb_Lvo"]("R,b,i");
        tmps_["93_bbbb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["30_bb_Lvv"]("L,a,c") * tmps_["90_bb_Lvo"]("R,b,i");
        tmps_["93_bbbb_LLvvvo"]("L,R,a,b,c,i") -= t1["bb"]("b,i") * tmps_["92_bb_LLvv"]("L,R,a,c");
        tmps_["20_bbbb_LLovvo"].~TArrayD();

        // D_ovvv_bbbb += -1.00 P(b,c) r1_bb(b,i) t1_bb(c,j) l1_bb(j,a)
        //             += -0.50 P(b,c) r1_bb(b,i) t2_abab(d,c,j,k) l2_abab(j,k,d,a)
        //             += -0.50 P(b,c) r1_bb(b,i) t2_abab(d,c,k,j) l2_abab(k,j,d,a)
        //             += +1.00 P(b,c) r1_bb(b,j) t1_bb(c,i) l1_bb(j,a)
        //             += -0.50 P(b,c) r1_bb(b,i) t2_bbbb(c,d,j,k) l2_bbbb(j,k,a,d)
        //             += +0.50 P(b,c) r2_abab(d,b,k,j) t1_bb(c,i) l2_abab(k,j,d,a)
        //             += +0.50 P(b,c) r2_abab(d,b,j,k) t1_bb(c,i) l2_abab(j,k,d,a)
        //             += +0.50 P(b,c) r2_bbbb(b,d,k,j) t1_bb(c,i) l2_bbbb(k,j,a,d)
        //             += +0.50 P(b,c) r0 t2_bbbb(b,d,j,k) t1_bb(c,i) l2_bbbb(j,k,a,d)
        //             += -1.00 P(b,c) r0 t2_abab(d,b,j,i) t1_bb(c,k) l2_abab(j,k,d,a)
        //             += -1.00 P(b,c) r2_bbbb(b,d,k,i) t1_bb(c,j) l2_bbbb(k,j,a,d)
        //             += -1.00 P(b,c) r2_abab(d,b,k,i) t1_bb(c,j) l2_abab(k,j,d,a)
        //             += -1.00 P(b,c) r1_bb(b,k) t2_bbbb(c,d,j,i) l2_bbbb(k,j,a,d)
        //             += -1.00 P(b,c) r1_bb(b,k) t1_bb(c,j) t1_bb(d,i) l2_bbbb(k,j,a,d)
        //             += +1.00 P(b,c) r1_bb(b,k) t2_abab(d,c,j,i) l2_abab(j,k,d,a)
        //             += -1.00 P(b,c) r0 t2_bbbb(b,d,j,i) t1_bb(c,k) l2_bbbb(j,k,a,d)
        //             += +1.00 P(b,c) r0 t1_bb(b,j) t1_bb(c,i) l1_bb(j,a)
        //             += +0.50 P(b,c) r0 t2_abab(d,b,j,k) t1_bb(c,i) l2_abab(j,k,d,a)
        //             += +0.50 P(b,c) r0 t2_abab(d,b,k,j) t1_bb(c,i) l2_abab(k,j,d,a)
        //             += +1.00 P(b,c) r1_aa(d,k) t1_bb(b,j) t1_bb(c,i) l2_abab(k,j,d,a)
        //             += -1.00 P(b,c) r1_bb(d,k) t1_bb(b,j) t1_bb(c,i) l2_bbbb(k,j,a,d)
        D_ovvv_bbbb("L,R,i,a,b,c") -= tmps_["93_bbbb_LLvvvo"]("L,R,a,c,b,i");
        D_ovvv_bbbb("L,R,i,a,b,c") += tmps_["93_bbbb_LLvvvo"]("L,R,a,b,c,i");

        // D_vovv_bbbb += +1.00 P(b,c) r1_bb(b,i) t1_bb(c,j) l1_bb(j,a)
        //             += +0.50 P(b,c) r1_bb(b,i) t2_abab(d,c,j,k) l2_abab(j,k,d,a)
        //             += +0.50 P(b,c) r1_bb(b,i) t2_abab(d,c,k,j) l2_abab(k,j,d,a)
        //             += -1.00 P(b,c) r1_bb(b,j) t1_bb(c,i) l1_bb(j,a)
        //             += +0.50 P(b,c) r1_bb(b,i) t2_bbbb(c,d,j,k) l2_bbbb(j,k,a,d)
        //             += -0.50 P(b,c) r2_abab(d,b,k,j) t1_bb(c,i) l2_abab(k,j,d,a)
        //             += -0.50 P(b,c) r2_abab(d,b,j,k) t1_bb(c,i) l2_abab(j,k,d,a)
        //             += -0.50 P(b,c) r2_bbbb(b,d,k,j) t1_bb(c,i) l2_bbbb(k,j,a,d)
        //             += -0.50 P(b,c) r0 t2_bbbb(b,d,j,k) t1_bb(c,i) l2_bbbb(j,k,a,d)
        //             += +1.00 P(b,c) r0 t2_abab(d,b,j,i) t1_bb(c,k) l2_abab(j,k,d,a)
        //             += +1.00 P(b,c) r2_bbbb(b,d,k,i) t1_bb(c,j) l2_bbbb(k,j,a,d)
        //             += +1.00 P(b,c) r2_abab(d,b,k,i) t1_bb(c,j) l2_abab(k,j,d,a)
        //             += +1.00 P(b,c) r1_bb(b,k) t2_bbbb(c,d,j,i) l2_bbbb(k,j,a,d)
        //             += +1.00 P(b,c) r1_bb(b,k) t1_bb(c,j) t1_bb(d,i) l2_bbbb(k,j,a,d)
        //             += -1.00 P(b,c) r1_bb(b,k) t2_abab(d,c,j,i) l2_abab(j,k,d,a)
        //             += +1.00 P(b,c) r0 t2_bbbb(b,d,j,i) t1_bb(c,k) l2_bbbb(j,k,a,d)
        //             += -1.00 P(b,c) r0 t1_bb(b,j) t1_bb(c,i) l1_bb(j,a)
        //             += -0.50 P(b,c) r0 t2_abab(d,b,j,k) t1_bb(c,i) l2_abab(j,k,d,a)
        //             += -0.50 P(b,c) r0 t2_abab(d,b,k,j) t1_bb(c,i) l2_abab(k,j,d,a)
        //             += -1.00 P(b,c) r1_aa(d,k) t1_bb(b,j) t1_bb(c,i) l2_abab(k,j,d,a)
        //             += +1.00 P(b,c) r1_bb(d,k) t1_bb(b,j) t1_bb(c,i) l2_bbbb(k,j,a,d)
        D_vovv_bbbb("L,R,a,i,b,c") += tmps_["93_bbbb_LLvvvo"]("L,R,a,c,b,i");
        D_vovv_bbbb("L,R,a,i,b,c") -= tmps_["93_bbbb_LLvvvo"]("L,R,a,b,c,i");
        tmps_["93_bbbb_LLvvvo"].~TArrayD();

        // flops: o0v4L1  = o2v4L1
        //  mems: o0v4L1  = o0v4L1
        tmps_["94_abab_Lvvvv"]("L,a,b,c,d")  = l2["abab"]("L,i,j,a,b") * t2["abab"]("c,d,i,j");

        // D_oovv_abab += -0.50 r1_bb(d,j) t2_abab(a,b,k,l) t1_aa(c,i) l2_abab(k,l,c,d)
        //             += -0.50 r1_bb(d,j) t2_abab(a,b,l,k) t1_aa(c,i) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o1v4L1 o2v3L2
        //  mems: o2v2L2 += o1v3L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["aa"]("c,i") * tmps_["94_abab_Lvvvv"]("L,c,d,a,b") * r1["bb"]("R,d,j");

        // D_vvvv_abab += -0.50 r0 t2_abab(c,d,i,j) l2_abab(i,j,a,b)
        //             += -0.50 r0 t2_abab(c,d,j,i) l2_abab(j,i,a,b)
        // flops: o0v4L2 += o0v4L2
        //  mems: o0v4L2 += o0v4L2
        D_vvvv_abab("L,R,a,b,c,d") -= r0("R") * tmps_["94_abab_Lvvvv"]("L,a,b,c,d");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["95_bbaa_Lvoov"]("L,a,i,j,b")  = l2["aaaa"]("L,k,j,b,c") * t2["abab"]("c,a,k,i");

        // D_oovv_abab += -1.00 r0 t1_aa(a,l) t2_abab(c,b,k,j) t1_aa(d,i) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o3v2L1 o2v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["aa"]("d,i") * tmps_["95_bbaa_Lvoov"]("L,b,j,l,d") * t1["aa"]("a,l") * r0("R");

        // D_oovv_bbbb += +1.00 P(i,j) r0 t2_abab(c,a,k,i) t2_abab(d,b,l,j) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["abab"]("d,b,l,j") * tmps_["95_bbaa_Lvoov"]("L,a,i,l,d") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v3L1  = o1v4L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["96_aabb_Lvvvo"]("L,a,b,c,i")  = t1["bb"]("d,i") * tmps_["94_abab_Lvvvv"]("L,a,d,b,c");
    }
} // hilbert
#endif