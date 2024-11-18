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


    void EOM_EE_RDM::rdm2_00_4() {

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


        // D_oovv_abab += +0.50 r0 t2_abab(a,c,k,l) t1_bb(b,j) t1_aa(d,i) l2_abab(k,l,d,c)
        //             += +0.50 r0 t2_abab(a,c,l,k) t1_bb(b,j) t1_aa(d,i) l2_abab(l,k,d,c)
        //             += +0.50 r0 t1_aa(a,l) t1_bb(b,j) t2_abab(d,c,i,k) l2_abab(l,k,d,c)
        //             += +0.50 r0 t1_aa(a,l) t1_bb(b,j) t2_abab(c,d,i,k) l2_abab(l,k,c,d)
        //             += +0.50 r0 t2_aaaa(a,c,k,l) t1_bb(b,j) t1_aa(d,i) l2_aaaa(k,l,d,c)
        //             += -1.00 r0 t2_abab(a,c,i,k) t1_bb(b,j) l1_bb(k,c)
        //             += +1.00 r0 t2_aaaa(a,c,k,i) t1_bb(b,j) l1_aa(k,c)
        //             += +1.00 r0 t1_aa(a,k) t1_bb(b,j) t1_aa(c,i) l1_aa(k,c)
        //             += +0.50 r0 t1_aa(a,l) t1_bb(b,j) t2_aaaa(d,c,k,i) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["90_bb_Lvo"]("R,b,j") * tmps_["128_aa_Lov"]("L,i,a");

        // D_oovv_abab += +0.50 r1_bb(b,j) t2_abab(a,c,k,l) t1_aa(d,i) l2_abab(k,l,d,c)
        //             += +0.50 r1_bb(b,j) t2_abab(a,c,l,k) t1_aa(d,i) l2_abab(l,k,d,c)
        //             += +0.50 r1_bb(b,j) t1_aa(a,l) t2_abab(d,c,i,k) l2_abab(l,k,d,c)
        //             += +0.50 r1_bb(b,j) t1_aa(a,l) t2_abab(c,d,i,k) l2_abab(l,k,c,d)
        //             += +0.50 r1_bb(b,j) t2_aaaa(a,c,k,l) t1_aa(d,i) l2_aaaa(k,l,d,c)
        //             += -1.00 r1_bb(b,j) t2_abab(a,c,i,k) l1_bb(k,c)
        //             += +1.00 r1_bb(b,j) t2_aaaa(a,c,k,i) l1_aa(k,c)
        //             += +1.00 r1_bb(b,j) t1_aa(a,k) t1_aa(c,i) l1_aa(k,c)
        //             += +0.50 r1_bb(b,j) t1_aa(a,l) t2_aaaa(d,c,k,i) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += r1["bb"]("R,b,j") * tmps_["128_aa_Lov"]("L,i,a");

        // flops: o1v1L2  = o2v1L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["129_aa_LLvo"]("L,R,a,i")  = t1["aa"]("a,j") * tmps_["62_aa_LLoo"]("L,R,j,i");

        // D_oovv_aaaa += -1.00 P(i,j) P(a,b) r1_aa(d,l) t1_aa(a,k) t1_aa(b,j) t1_aa(c,i) l2_aaaa(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_bb(d,l) t1_aa(a,k) t1_aa(b,j) t1_aa(c,i) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("b,j") * tmps_["129_aa_LLvo"]("L,R,a,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r1_aa(d,l) t1_aa(a,k) t1_bb(b,j) t1_aa(c,i) l2_aaaa(l,k,c,d)
        //             += +1.00 r1_bb(d,l) t1_aa(a,k) t1_bb(b,j) t1_aa(c,i) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["bb"]("b,j") * tmps_["129_aa_LLvo"]("L,R,a,i");

        // flops: o1v1L2  = o1v2L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["130_aa_LLov"]("R,L,i,a")  = r1["aa"]("R,b,i") * tmps_["98_aa_Lvv"]("L,b,a");
        tmps_["98_aa_Lvv"].~TArrayD();

        // D_oovv_aaaa += +0.50 P(i,j) P(a,b) r1_aa(d,i) t2_abab(a,c,k,l) t1_aa(b,j) l2_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(d,i) t2_abab(a,c,l,k) t1_aa(b,j) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("b,j") * tmps_["130_aa_LLov"]("R,L,i,a");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r1_aa(d,i) t2_abab(a,c,k,l) t1_bb(b,j) l2_abab(k,l,d,c)
        //             += +0.50 r1_aa(d,i) t2_abab(a,c,l,k) t1_bb(b,j) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1["bb"]("b,j") * tmps_["130_aa_LLov"]("R,L,i,a");

        // flops: o2v0L1  = o2v0L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["124_bb_Loo"]("L,i,j")  = Id["bb_oo"]("i,j") * l0("L");

        // flops: o3v1L2  = o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o2v0L1 o3v1L2 o3v1L2 o3v1L2 o3v1 o3v1L2 o3v1L2 o4v1L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o4v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o4v1L2 o3v1L2 o4v1L1 o3v1L2 o3v1L2 o4v0L1 o4v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        //  mems: o3v1L2  = o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o2v0L1 o3v1L2 o3v1L2 o3v1L2 o3v1 o3v1L2 o3v1L2 o3v1L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L1 o3v1L2 o3v1L2 o4v0L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k")  = -1.00 * tmps_["127_aa_LLvo"]("L,R,a,k") * Id["bb_oo"]("i,j");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") -= tmps_["128_aa_Lov"]("L,k,a") * tmps_["52_bb_Loo"]("R,i,j");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") += tmps_["120_aabb_Lvooo"]("L,a,k,j,i") * r0("R");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") += t1["aa"]("a,k") * tmps_["53_bb_LLoo"]("L,R,j,i");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") -= 0.50 * (tmps_["49_bb_Loo"]("L,i,j") + 2.00 * tmps_["48_bb_Loo"]("L,i,j")) * tmps_["109_aa_Lvo"]("R,a,k");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") -= tmps_["125_aa_LLov"]("L,R,k,a") * Id["bb_oo"]("i,j");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") += Id["bb_oo"]("i,j") * t1["aa"]("a,k") * tmps_["126_LL"]("L,R");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") += t1["aa"]("a,l") * tmps_["46_aabb_Loooo"]("L,k,l,j,i") * r0("R");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") += tmps_["124_bb_Loo"]("L,i,j") * tmps_["109_aa_Lvo"]("R,a,k");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") += t1["aa"]("a,l") * tmps_["116_abab_LLoooo"]("L,R,l,j,k,i");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") -= t1["aa"]("a,k") * tmps_["51_bb_LLoo"]("L,R,j,i");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") -= tmps_["50_bb_Loo"]("L,i,j") * tmps_["109_aa_Lvo"]("R,a,k");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") += t1["aa"]("a,l") * tmps_["119_aabb_LLoooo"]("L,R,k,l,j,i");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") += t1["aa"]("a,l") * tmps_["45_abab_Loooo"]("L,k,i,l,j") * r0("R");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") += (tmps_["45_abab_Loooo"]("L,k,i,l,j") + tmps_["46_aabb_Loooo"]("L,k,l,j,i")) * r1["aa"]("R,a,l");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") -= tmps_["118_aabb_Lvooo"]("L,a,k,i,j") * r0("R");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") -= 0.50 * tmps_["49_bb_Loo"]("L,i,j") * r1["aa"]("R,a,k");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") -= tmps_["50_bb_Loo"]("L,i,j") * r1["aa"]("R,a,k");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") += tmps_["129_aa_LLvo"]("L,R,a,k") * Id["bb_oo"]("i,j");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") += tmps_["117_bbaa_Loovo"]("L,i,j,a,k") * r0("R");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") += tmps_["124_bb_Loo"]("L,i,j") * r1["aa"]("R,a,k");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") -= tmps_["48_bb_Loo"]("L,i,j") * r1["aa"]("R,a,k");
        tmps_["131_bbaa_LLoovo"]("L,R,i,j,a,k") -= tmps_["130_aa_LLov"]("R,L,k,a") * Id["bb_oo"]("i,j");
        tmps_["120_aabb_Lvooo"].~TArrayD();
        tmps_["119_aabb_LLoooo"].~TArrayD();
        tmps_["118_aabb_Lvooo"].~TArrayD();
        tmps_["117_bbaa_Loovo"].~TArrayD();
        tmps_["116_abab_LLoooo"].~TArrayD();
        tmps_["46_aabb_Loooo"].~TArrayD();
        tmps_["45_abab_Loooo"].~TArrayD();

        // D_ooov_abab += -1.00 d_bb(i,k) r1_aa(c,m) t1_aa(a,l) t1_aa(b,j) l2_aaaa(m,l,b,c)
        //             += +1.00 d_bb(i,k) r1_bb(c,m) t1_aa(a,l) t1_aa(b,j) l2_abab(l,m,b,c)
        //             += +1.00 r1_aa(a,j) t1_bb(b,i) l1_bb(k,b)
        //             += +0.50 r1_aa(a,j) t2_bbbb(c,b,l,i) l2_bbbb(l,k,c,b)
        //             += -1.00 r0 t2_abab(a,b,j,l) t1_bb(c,i) l2_bbbb(l,k,c,b)
        //             += +1.00 r0 t2_abab(a,b,j,i) l1_bb(k,b)
        //             += -1.00 d_bb(i,k) r1_aa(a,j) l0
        //             += +1.00 d_bb(i,k) r1_aa(b,j) t1_aa(a,l) l1_aa(l,b)
        //             += +0.50 d_bb(i,k) r2_abab(b,c,j,m) t1_aa(a,l) l2_abab(l,m,b,c)
        //             += +0.50 d_bb(i,k) r2_abab(c,b,j,m) t1_aa(a,l) l2_abab(l,m,c,b)
        //             += +0.50 d_bb(i,k) r2_aaaa(b,c,m,j) t1_aa(a,l) l2_aaaa(m,l,b,c)
        //             += +1.00 d_bb(i,k) r1_bb(c,m) t2_aaaa(a,b,l,j) l2_abab(l,m,b,c)
        //             += +0.50 d_bb(i,k) r1_aa(a,m) t2_abab(c,b,j,l) l2_abab(m,l,c,b)
        //             += +0.50 d_bb(i,k) r1_aa(a,m) t2_abab(b,c,j,l) l2_abab(m,l,b,c)
        //             += +1.00 d_bb(i,k) r1_bb(c,m) t2_abab(a,b,j,l) l2_bbbb(m,l,b,c)
        //             += -1.00 d_bb(i,k) r1_aa(c,m) t2_abab(a,b,j,l) l2_abab(m,l,c,b)
        //             += -0.50 d_bb(i,k) r1_aa(a,m) t2_aaaa(c,b,l,j) l2_aaaa(m,l,c,b)
        //             += -1.00 d_bb(i,k) r1_aa(c,m) t2_aaaa(a,b,l,j) l2_aaaa(m,l,b,c)
        //             += +0.50 d_bb(i,k) r2_abab(a,c,m,l) t1_aa(b,j) l2_abab(m,l,b,c)
        //             += +0.50 d_bb(i,k) r2_abab(a,c,l,m) t1_aa(b,j) l2_abab(l,m,b,c)
        //             += +0.50 d_bb(i,k) r2_aaaa(a,c,m,l) t1_aa(b,j) l2_aaaa(m,l,b,c)
        //             += -0.50 d_bb(i,k) r1_aa(c,j) t2_aaaa(a,b,l,m) l2_aaaa(l,m,b,c)
        //             += +0.50 d_bb(i,k) r0 t2_abab(a,b,l,m) t1_aa(c,j) l2_abab(l,m,c,b)
        //             += +0.50 d_bb(i,k) r0 t2_abab(a,b,m,l) t1_aa(c,j) l2_abab(m,l,c,b)
        //             += +0.50 d_bb(i,k) r0 t1_aa(a,m) t2_abab(c,b,j,l) l2_abab(m,l,c,b)
        //             += +0.50 d_bb(i,k) r0 t1_aa(a,m) t2_abab(b,c,j,l) l2_abab(m,l,b,c)
        //             += +0.50 d_bb(i,k) r0 t2_aaaa(a,b,l,m) t1_aa(c,j) l2_aaaa(l,m,c,b)
        //             += -1.00 d_bb(i,k) r0 t2_abab(a,b,j,l) l1_bb(l,b)
        //             += +1.00 d_bb(i,k) r0 t2_aaaa(a,b,l,j) l1_aa(l,b)
        //             += +1.00 d_bb(i,k) r0 t1_aa(a,l) t1_aa(b,j) l1_aa(l,b)
        //             += +0.50 d_bb(i,k) r0 t1_aa(a,m) t2_aaaa(c,b,l,j) l2_aaaa(l,m,c,b)
        //             += -1.00 r0 t2_aaaa(a,b,l,j) t1_bb(c,i) l2_abab(l,k,b,c)
        //             += -1.00 r0 t2_abab(a,b,l,i) t1_aa(c,j) l2_abab(l,k,c,b)
        //             += -1.00 r1_bb(c,l) t1_aa(a,j) t1_bb(b,i) l2_bbbb(l,k,b,c)
        //             += +1.00 r1_aa(c,l) t1_aa(a,j) t1_bb(b,i) l2_abab(l,k,c,b)
        //             += +0.50 r0 t1_aa(a,j) t2_bbbb(c,b,l,i) l2_bbbb(l,k,c,b)
        //             += +0.50 r0 t1_aa(a,j) t2_abab(c,b,l,i) l2_abab(l,k,c,b)
        //             += +0.50 r0 t1_aa(a,j) t2_abab(b,c,l,i) l2_abab(l,k,b,c)
        //             += +1.00 d_bb(i,k) r1_aa(a,l) t1_aa(b,j) l1_aa(l,b)
        //             += -1.00 d_bb(i,k) r2_abab(a,b,j,l) l1_bb(l,b)
        //             += +1.00 d_bb(i,k) r2_aaaa(a,b,l,j) l1_aa(l,b)
        //             += -1.00 d_bb(i,k) r1_bb(b,l) t1_aa(a,j) l1_bb(l,b)
        //             += -1.00 d_bb(i,k) r1_aa(b,l) t1_aa(a,j) l1_aa(l,b)
        //             += -0.25 d_bb(i,k) r2_aaaa(b,c,m,l) t1_aa(a,j) l2_aaaa(m,l,b,c)
        //             += -0.25 d_bb(i,k) r2_abab(b,c,m,l) t1_aa(a,j) l2_abab(m,l,b,c)
        //             += -0.25 d_bb(i,k) r2_abab(b,c,l,m) t1_aa(a,j) l2_abab(l,m,b,c)
        //             += -0.25 d_bb(i,k) r2_abab(c,b,m,l) t1_aa(a,j) l2_abab(m,l,c,b)
        //             += -0.25 d_bb(i,k) r2_abab(c,b,l,m) t1_aa(a,j) l2_abab(l,m,c,b)
        //             += -0.25 d_bb(i,k) r2_bbbb(b,c,m,l) t1_aa(a,j) l2_bbbb(m,l,b,c)
        //             += -1.00 r0 t1_aa(a,l) t1_bb(b,i) t1_aa(c,j) l2_abab(l,k,c,b)
        //             += -1.00 d_bb(i,k) r0 t1_aa(a,j) l0
        //             += -0.50 r2_abab(b,c,j,i) t1_aa(a,l) l2_abab(l,k,b,c)
        //             += -0.50 r2_abab(c,b,j,i) t1_aa(a,l) l2_abab(l,k,c,b)
        //             += +1.00 r1_bb(b,i) t1_aa(a,j) l1_bb(k,b)
        //             += +0.50 r2_abab(b,c,l,i) t1_aa(a,j) l2_abab(l,k,b,c)
        //             += +0.50 r2_abab(c,b,l,i) t1_aa(a,j) l2_abab(l,k,c,b)
        //             += +0.50 r2_bbbb(b,c,l,i) t1_aa(a,j) l2_bbbb(l,k,b,c)
        //             += +1.00 r0 t1_aa(a,j) t1_bb(b,i) l1_bb(k,b)
        //             += -1.00 r1_bb(c,i) t1_aa(a,l) t1_aa(b,j) l2_abab(l,k,b,c)
        //             += -1.00 r1_aa(c,j) t1_aa(a,l) t1_bb(b,i) l2_abab(l,k,c,b)
        //             += -0.50 r0 t1_aa(a,l) t2_abab(c,b,j,i) l2_abab(l,k,c,b)
        //             += -0.50 r0 t1_aa(a,l) t2_abab(b,c,j,i) l2_abab(l,k,b,c)
        //             += -0.50 r1_aa(a,l) t2_abab(c,b,j,i) l2_abab(l,k,c,b)
        //             += -0.50 r1_aa(a,l) t2_abab(b,c,j,i) l2_abab(l,k,b,c)
        //             += -1.00 r1_aa(a,l) t1_bb(b,i) t1_aa(c,j) l2_abab(l,k,c,b)
        //             += +0.50 r1_aa(a,j) t2_abab(c,b,l,i) l2_abab(l,k,c,b)
        //             += +0.50 r1_aa(a,j) t2_abab(b,c,l,i) l2_abab(l,k,b,c)
        //             += +0.50 d_bb(i,k) r1_aa(c,j) t2_abab(a,b,l,m) l2_abab(l,m,c,b)
        //             += +0.50 d_bb(i,k) r1_aa(c,j) t2_abab(a,b,m,l) l2_abab(m,l,c,b)
        D_ooov_abab("L,R,i,j,k,a") -= tmps_["131_bbaa_LLoovo"]("L,R,i,k,a,j");

        // D_oovo_abab += +1.00 d_bb(i,k) r1_aa(c,m) t1_aa(a,l) t1_aa(b,j) l2_aaaa(m,l,b,c)
        //             += -1.00 d_bb(i,k) r1_bb(c,m) t1_aa(a,l) t1_aa(b,j) l2_abab(l,m,b,c)
        //             += -1.00 r1_aa(a,j) t1_bb(b,i) l1_bb(k,b)
        //             += -0.50 r1_aa(a,j) t2_bbbb(c,b,l,i) l2_bbbb(l,k,c,b)
        //             += +1.00 r0 t2_abab(a,b,j,l) t1_bb(c,i) l2_bbbb(l,k,c,b)
        //             += -1.00 r0 t2_abab(a,b,j,i) l1_bb(k,b)
        //             += +1.00 d_bb(i,k) r1_aa(a,j) l0
        //             += -1.00 d_bb(i,k) r1_aa(b,j) t1_aa(a,l) l1_aa(l,b)
        //             += -0.50 d_bb(i,k) r2_abab(b,c,j,m) t1_aa(a,l) l2_abab(l,m,b,c)
        //             += -0.50 d_bb(i,k) r2_abab(c,b,j,m) t1_aa(a,l) l2_abab(l,m,c,b)
        //             += -0.50 d_bb(i,k) r2_aaaa(b,c,m,j) t1_aa(a,l) l2_aaaa(m,l,b,c)
        //             += -1.00 d_bb(i,k) r1_bb(c,m) t2_aaaa(a,b,l,j) l2_abab(l,m,b,c)
        //             += -0.50 d_bb(i,k) r1_aa(a,m) t2_abab(c,b,j,l) l2_abab(m,l,c,b)
        //             += -0.50 d_bb(i,k) r1_aa(a,m) t2_abab(b,c,j,l) l2_abab(m,l,b,c)
        //             += -1.00 d_bb(i,k) r1_bb(c,m) t2_abab(a,b,j,l) l2_bbbb(m,l,b,c)
        //             += +1.00 d_bb(i,k) r1_aa(c,m) t2_abab(a,b,j,l) l2_abab(m,l,c,b)
        //             += +0.50 d_bb(i,k) r1_aa(a,m) t2_aaaa(c,b,l,j) l2_aaaa(m,l,c,b)
        //             += +1.00 d_bb(i,k) r1_aa(c,m) t2_aaaa(a,b,l,j) l2_aaaa(m,l,b,c)
        //             += -0.50 d_bb(i,k) r2_abab(a,c,m,l) t1_aa(b,j) l2_abab(m,l,b,c)
        //             += -0.50 d_bb(i,k) r2_abab(a,c,l,m) t1_aa(b,j) l2_abab(l,m,b,c)
        //             += -0.50 d_bb(i,k) r2_aaaa(a,c,m,l) t1_aa(b,j) l2_aaaa(m,l,b,c)
        //             += +0.50 d_bb(i,k) r1_aa(c,j) t2_aaaa(a,b,l,m) l2_aaaa(l,m,b,c)
        //             += -0.50 d_bb(i,k) r0 t2_abab(a,b,l,m) t1_aa(c,j) l2_abab(l,m,c,b)
        //             += -0.50 d_bb(i,k) r0 t2_abab(a,b,m,l) t1_aa(c,j) l2_abab(m,l,c,b)
        //             += -0.50 d_bb(i,k) r0 t1_aa(a,m) t2_abab(c,b,j,l) l2_abab(m,l,c,b)
        //             += -0.50 d_bb(i,k) r0 t1_aa(a,m) t2_abab(b,c,j,l) l2_abab(m,l,b,c)
        //             += -0.50 d_bb(i,k) r0 t2_aaaa(a,b,l,m) t1_aa(c,j) l2_aaaa(l,m,c,b)
        //             += +1.00 d_bb(i,k) r0 t2_abab(a,b,j,l) l1_bb(l,b)
        //             += -1.00 d_bb(i,k) r0 t2_aaaa(a,b,l,j) l1_aa(l,b)
        //             += -1.00 d_bb(i,k) r0 t1_aa(a,l) t1_aa(b,j) l1_aa(l,b)
        //             += -0.50 d_bb(i,k) r0 t1_aa(a,m) t2_aaaa(c,b,l,j) l2_aaaa(l,m,c,b)
        //             += +1.00 r0 t2_aaaa(a,b,l,j) t1_bb(c,i) l2_abab(l,k,b,c)
        //             += +1.00 r0 t2_abab(a,b,l,i) t1_aa(c,j) l2_abab(l,k,c,b)
        //             += +1.00 r1_bb(c,l) t1_aa(a,j) t1_bb(b,i) l2_bbbb(l,k,b,c)
        //             += -1.00 r1_aa(c,l) t1_aa(a,j) t1_bb(b,i) l2_abab(l,k,c,b)
        //             += -0.50 r0 t1_aa(a,j) t2_bbbb(c,b,l,i) l2_bbbb(l,k,c,b)
        //             += -0.50 r0 t1_aa(a,j) t2_abab(c,b,l,i) l2_abab(l,k,c,b)
        //             += -0.50 r0 t1_aa(a,j) t2_abab(b,c,l,i) l2_abab(l,k,b,c)
        //             += -1.00 d_bb(i,k) r1_aa(a,l) t1_aa(b,j) l1_aa(l,b)
        //             += +1.00 d_bb(i,k) r2_abab(a,b,j,l) l1_bb(l,b)
        //             += -1.00 d_bb(i,k) r2_aaaa(a,b,l,j) l1_aa(l,b)
        //             += +1.00 d_bb(i,k) r1_bb(b,l) t1_aa(a,j) l1_bb(l,b)
        //             += +1.00 d_bb(i,k) r1_aa(b,l) t1_aa(a,j) l1_aa(l,b)
        //             += +0.25 d_bb(i,k) r2_aaaa(b,c,m,l) t1_aa(a,j) l2_aaaa(m,l,b,c)
        //             += +0.25 d_bb(i,k) r2_abab(b,c,m,l) t1_aa(a,j) l2_abab(m,l,b,c)
        //             += +0.25 d_bb(i,k) r2_abab(b,c,l,m) t1_aa(a,j) l2_abab(l,m,b,c)
        //             += +0.25 d_bb(i,k) r2_abab(c,b,m,l) t1_aa(a,j) l2_abab(m,l,c,b)
        //             += +0.25 d_bb(i,k) r2_abab(c,b,l,m) t1_aa(a,j) l2_abab(l,m,c,b)
        //             += +0.25 d_bb(i,k) r2_bbbb(b,c,m,l) t1_aa(a,j) l2_bbbb(m,l,b,c)
        //             += +1.00 r0 t1_aa(a,l) t1_bb(b,i) t1_aa(c,j) l2_abab(l,k,c,b)
        //             += +1.00 d_bb(i,k) r0 t1_aa(a,j) l0
        //             += +0.50 r2_abab(b,c,j,i) t1_aa(a,l) l2_abab(l,k,b,c)
        //             += +0.50 r2_abab(c,b,j,i) t1_aa(a,l) l2_abab(l,k,c,b)
        //             += -1.00 r1_bb(b,i) t1_aa(a,j) l1_bb(k,b)
        //             += -0.50 r2_abab(b,c,l,i) t1_aa(a,j) l2_abab(l,k,b,c)
        //             += -0.50 r2_abab(c,b,l,i) t1_aa(a,j) l2_abab(l,k,c,b)
        //             += -0.50 r2_bbbb(b,c,l,i) t1_aa(a,j) l2_bbbb(l,k,b,c)
        //             += -1.00 r0 t1_aa(a,j) t1_bb(b,i) l1_bb(k,b)
        //             += +1.00 r1_bb(c,i) t1_aa(a,l) t1_aa(b,j) l2_abab(l,k,b,c)
        //             += +1.00 r1_aa(c,j) t1_aa(a,l) t1_bb(b,i) l2_abab(l,k,c,b)
        //             += +0.50 r0 t1_aa(a,l) t2_abab(c,b,j,i) l2_abab(l,k,c,b)
        //             += +0.50 r0 t1_aa(a,l) t2_abab(b,c,j,i) l2_abab(l,k,b,c)
        //             += +0.50 r1_aa(a,l) t2_abab(c,b,j,i) l2_abab(l,k,c,b)
        //             += +0.50 r1_aa(a,l) t2_abab(b,c,j,i) l2_abab(l,k,b,c)
        //             += +1.00 r1_aa(a,l) t1_bb(b,i) t1_aa(c,j) l2_abab(l,k,c,b)
        //             += -0.50 r1_aa(a,j) t2_abab(c,b,l,i) l2_abab(l,k,c,b)
        //             += -0.50 r1_aa(a,j) t2_abab(b,c,l,i) l2_abab(l,k,b,c)
        //             += -0.50 d_bb(i,k) r1_aa(c,j) t2_abab(a,b,l,m) l2_abab(l,m,c,b)
        //             += -0.50 d_bb(i,k) r1_aa(c,j) t2_abab(a,b,m,l) l2_abab(m,l,c,b)
        D_oovo_abab("L,R,i,j,a,k") += tmps_["131_bbaa_LLoovo"]("L,R,i,k,a,j");
        tmps_["131_bbaa_LLoovo"].~TArrayD();

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["133_aaaa_Lvooo"]("L,a,i,j,k")  = t1["aa"]("b,k") * tmps_["108_aaaa_Lvoov"]("L,a,i,j,b");

        // D_oovv_aaaa += -1.00 P(i,j) P(a,b) r1_aa(a,l) t2_abab(b,c,i,k) t1_aa(d,j) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = r1["aa"]("R,a,l") * tmps_["133_aaaa_Lvooo"]("L,b,i,l,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r0 t2_abab(a,c,i,k) t1_aa(b,l) t1_aa(d,j) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("b,l") * tmps_["133_aaaa_Lvooo"]("L,a,i,l,j") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["134_aaaa_Lvooo"]("L,a,i,j,k")  = t1["aa"]("b,k") * tmps_["107_aaaa_Lvoov"]("L,a,i,j,b");
        tmps_["107_aaaa_Lvoov"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r0 t2_aaaa(a,c,k,i) t1_aa(b,l) t1_aa(d,j) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("b,l") * tmps_["134_aaaa_Lvooo"]("L,a,i,l,j") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o2v0L1  = o2v0L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["132_aa_Loo"]("L,i,j")  = Id["aa_oo"]("i,j") * l0("L");

        // flops: o3v1L2  = o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v3L1 o2v2L1 o3v2L2 o3v1L2 o3v1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        //  mems: o3v1L2  = o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o2v2L1 o2v2L1 o3v1L2 o3v1L2 o3v1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k")  = -1.00 * Id["aa_oo"]("k,j") * tmps_["129_aa_LLvo"]("L,R,a,i");
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k") += r0("R") * tmps_["133_aaaa_Lvooo"]("L,a,i,j,k");
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k") += 0.50 * r1["aa"]("R,a,i") * tmps_["57_aa_Loo"]("L,k,j");
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k") += tmps_["127_aa_LLvo"]("L,R,a,i") * Id["aa_oo"]("k,j");
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k") += t1["aa"]("a,k") * tmps_["62_aa_LLoo"]("L,R,j,i");
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k") -= tmps_["132_aa_Loo"]("L,k,j") * tmps_["109_aa_Lvo"]("R,a,i");
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k") -= tmps_["59_aa_Loo"]("L,i,j") * tmps_["109_aa_Lvo"]("R,a,k");
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k") += (t2["aaaa"]("a,c,l,k") * l2["aaaa"]("L,l,j,c,b") + -1.00 * tmps_["108_aaaa_Lvoov"]("L,a,k,j,b")) * r1["aa"]("R,b,i");
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k") -= Id["aa_oo"]("k,j") * t1["aa"]("a,i") * tmps_["126_LL"]("L,R");
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k") += tmps_["125_aa_LLov"]("L,R,i,a") * Id["aa_oo"]("k,j");
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k") -= 0.50 * tmps_["57_aa_Loo"]("L,i,j") * tmps_["109_aa_Lvo"]("R,a,k");
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k") -= t1["aa"]("a,k") * tmps_["60_aa_LLoo"]("L,R,j,i");
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k") -= tmps_["58_aa_Loo"]("L,i,j") * tmps_["109_aa_Lvo"]("R,a,k");
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k") += tmps_["128_aa_Lov"]("L,i,a") * tmps_["61_aa_Loo"]("R,k,j");
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k") += r0("R") * tmps_["134_aaaa_Lvooo"]("L,a,i,j,k");
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k") -= r1["aa"]("R,a,i") * tmps_["132_aa_Loo"]("L,k,j");
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k") += r1["aa"]("R,a,i") * tmps_["58_aa_Loo"]("L,k,j");
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k") += r1["aa"]("R,a,i") * tmps_["59_aa_Loo"]("L,k,j");
        tmps_["135_aaaa_LLvooo"]("R,L,a,i,j,k") += Id["aa_oo"]("k,j") * tmps_["130_aa_LLov"]("R,L,i,a");
        tmps_["134_aaaa_Lvooo"].~TArrayD();
        tmps_["133_aaaa_Lvooo"].~TArrayD();
        tmps_["132_aa_Loo"].~TArrayD();
        tmps_["130_aa_LLov"].~TArrayD();
        tmps_["129_aa_LLvo"].~TArrayD();
        tmps_["128_aa_Lov"].~TArrayD();
        tmps_["127_aa_LLvo"].~TArrayD();
        tmps_["125_aa_LLov"].~TArrayD();
        tmps_["108_aaaa_Lvoov"].~TArrayD();
        tmps_["62_aa_LLoo"].~TArrayD();
        tmps_["61_aa_Loo"].~TArrayD();
        tmps_["60_aa_LLoo"].~TArrayD();
        tmps_["59_aa_Loo"].~TArrayD();
        tmps_["58_aa_Loo"].~TArrayD();
        tmps_["57_aa_Loo"].~TArrayD();

        // D_ooov_aaaa += -1.00 P(i,j) r0 t2_abab(a,b,i,l) t1_aa(c,j) l2_abab(k,l,c,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t1_aa(a,l) t1_aa(b,i) l2_aaaa(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t1_aa(a,l) t1_aa(b,i) l2_abab(l,m,b,c)
        //             += -0.50 P(i,j) r1_aa(a,i) t2_aaaa(c,b,l,j) l2_aaaa(l,k,c,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(b,i) t1_aa(a,l) l1_aa(l,b)
        //             += -0.50 P(i,j) d_aa(j,k) r2_abab(b,c,i,m) t1_aa(a,l) l2_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r2_abab(c,b,i,m) t1_aa(a,l) l2_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r2_aaaa(b,c,m,i) t1_aa(a,l) l2_aaaa(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t2_aaaa(a,b,l,i) l2_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r1_aa(a,m) t2_abab(c,b,i,l) l2_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r1_aa(a,m) t2_abab(b,c,i,l) l2_abab(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t2_abab(a,b,i,l) l2_bbbb(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t2_abab(a,b,i,l) l2_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r1_aa(a,m) t2_aaaa(c,b,l,i) l2_aaaa(m,l,c,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t2_aaaa(a,b,l,i) l2_aaaa(m,l,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r2_abab(a,c,m,l) t1_aa(b,i) l2_abab(m,l,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r2_abab(a,c,l,m) t1_aa(b,i) l2_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r2_aaaa(a,c,m,l) t1_aa(b,i) l2_aaaa(m,l,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r1_aa(c,i) t2_aaaa(a,b,l,m) l2_aaaa(l,m,b,c)
        //             += -1.00 P(i,j) r1_aa(c,l) t1_aa(a,j) t1_aa(b,i) l2_aaaa(l,k,b,c)
        //             += +1.00 P(i,j) r1_bb(c,l) t1_aa(a,j) t1_aa(b,i) l2_abab(k,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r0 t1_aa(a,i) l0
        //             += +1.00 P(i,j) r0 t1_aa(a,j) t1_aa(b,i) l1_aa(k,b)
        //             += -1.00 P(i,j) r1_aa(c,i) t2_aaaa(a,b,l,j) l2_aaaa(l,k,b,c)
        //             += +1.00 P(i,j) r1_aa(c,i) t2_abab(a,b,j,l) l2_abab(k,l,c,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_bb(b,l) t1_aa(a,i) l1_bb(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(b,l) t1_aa(a,i) l1_aa(l,b)
        //             += +0.25 P(i,j) d_aa(j,k) r2_aaaa(b,c,m,l) t1_aa(a,i) l2_aaaa(m,l,b,c)
        //             += +0.25 P(i,j) d_aa(j,k) r2_abab(b,c,m,l) t1_aa(a,i) l2_abab(m,l,b,c)
        //             += +0.25 P(i,j) d_aa(j,k) r2_abab(b,c,l,m) t1_aa(a,i) l2_abab(l,m,b,c)
        //             += +0.25 P(i,j) d_aa(j,k) r2_abab(c,b,m,l) t1_aa(a,i) l2_abab(m,l,c,b)
        //             += +0.25 P(i,j) d_aa(j,k) r2_abab(c,b,l,m) t1_aa(a,i) l2_abab(l,m,c,b)
        //             += +0.25 P(i,j) d_aa(j,k) r2_bbbb(b,c,m,l) t1_aa(a,i) l2_bbbb(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(a,l) t1_aa(b,i) l1_aa(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r2_abab(a,b,i,l) l1_bb(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r2_aaaa(a,b,l,i) l1_aa(l,b)
        //             += +0.50 P(i,j) r0 t1_aa(a,j) t2_aaaa(c,b,l,i) l2_aaaa(l,k,c,b)
        //             += +1.00 P(i,j) r1_aa(b,i) t1_aa(a,j) l1_aa(k,b)
        //             += +0.50 P(i,j) r2_abab(b,c,i,l) t1_aa(a,j) l2_abab(k,l,b,c)
        //             += +0.50 P(i,j) r2_abab(c,b,i,l) t1_aa(a,j) l2_abab(k,l,c,b)
        //             += +0.50 P(i,j) r2_aaaa(b,c,l,i) t1_aa(a,j) l2_aaaa(l,k,b,c)
        //             += +0.50 P(i,j) r0 t1_aa(a,j) t2_abab(c,b,i,l) l2_abab(k,l,c,b)
        //             += +0.50 P(i,j) r0 t1_aa(a,j) t2_abab(b,c,i,l) l2_abab(k,l,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_abab(a,b,l,m) t1_aa(c,i) l2_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_abab(a,b,m,l) t1_aa(c,i) l2_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t1_aa(a,m) t2_abab(c,b,i,l) l2_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t1_aa(a,m) t2_abab(b,c,i,l) l2_abab(m,l,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t2_aaaa(a,b,l,m) t1_aa(c,i) l2_aaaa(l,m,c,b)
        //             += +1.00 P(i,j) d_aa(j,k) r0 t2_abab(a,b,i,l) l1_bb(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r0 t2_aaaa(a,b,l,i) l1_aa(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r0 t1_aa(a,l) t1_aa(b,i) l1_aa(l,b)
        //             += -0.50 P(i,j) d_aa(j,k) r0 t1_aa(a,m) t2_aaaa(c,b,l,i) l2_aaaa(l,m,c,b)
        //             += -1.00 P(i,j) r0 t2_aaaa(a,b,l,i) t1_aa(c,j) l2_aaaa(l,k,c,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(a,i) l0
        //             += -0.50 P(i,j) r1_aa(a,i) t2_abab(c,b,j,l) l2_abab(k,l,c,b)
        //             += -0.50 P(i,j) r1_aa(a,i) t2_abab(b,c,j,l) l2_abab(k,l,b,c)
        //             += -1.00 P(i,j) r1_aa(a,i) t1_aa(b,j) l1_aa(k,b)
        //             += -0.50 P(i,j) d_aa(j,k) r1_aa(c,i) t2_abab(a,b,l,m) l2_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r1_aa(c,i) t2_abab(a,b,m,l) l2_abab(m,l,c,b)
        D_ooov_aaaa("L,R,i,j,k,a") -= tmps_["135_aaaa_LLvooo"]("R,L,a,i,k,j");
        D_ooov_aaaa("L,R,i,j,k,a") += tmps_["135_aaaa_LLvooo"]("R,L,a,j,k,i");

        // D_oovo_aaaa += +1.00 P(i,j) r0 t2_abab(a,b,i,l) t1_aa(c,j) l2_abab(k,l,c,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t1_aa(a,l) t1_aa(b,i) l2_aaaa(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t1_aa(a,l) t1_aa(b,i) l2_abab(l,m,b,c)
        //             += +0.50 P(i,j) r1_aa(a,i) t2_aaaa(c,b,l,j) l2_aaaa(l,k,c,b)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(b,i) t1_aa(a,l) l1_aa(l,b)
        //             += +0.50 P(i,j) d_aa(j,k) r2_abab(b,c,i,m) t1_aa(a,l) l2_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r2_abab(c,b,i,m) t1_aa(a,l) l2_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r2_aaaa(b,c,m,i) t1_aa(a,l) l2_aaaa(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t2_aaaa(a,b,l,i) l2_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r1_aa(a,m) t2_abab(c,b,i,l) l2_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r1_aa(a,m) t2_abab(b,c,i,l) l2_abab(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_bb(c,m) t2_abab(a,b,i,l) l2_bbbb(m,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t2_abab(a,b,i,l) l2_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_aa(j,k) r1_aa(a,m) t2_aaaa(c,b,l,i) l2_aaaa(m,l,c,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(c,m) t2_aaaa(a,b,l,i) l2_aaaa(m,l,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r2_abab(a,c,m,l) t1_aa(b,i) l2_abab(m,l,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r2_abab(a,c,l,m) t1_aa(b,i) l2_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r2_aaaa(a,c,m,l) t1_aa(b,i) l2_aaaa(m,l,b,c)
        //             += -0.50 P(i,j) d_aa(j,k) r1_aa(c,i) t2_aaaa(a,b,l,m) l2_aaaa(l,m,b,c)
        //             += +1.00 P(i,j) r1_aa(c,l) t1_aa(a,j) t1_aa(b,i) l2_aaaa(l,k,b,c)
        //             += -1.00 P(i,j) r1_bb(c,l) t1_aa(a,j) t1_aa(b,i) l2_abab(k,l,b,c)
        //             += -1.00 P(i,j) d_aa(j,k) r0 t1_aa(a,i) l0
        //             += -1.00 P(i,j) r0 t1_aa(a,j) t1_aa(b,i) l1_aa(k,b)
        //             += +1.00 P(i,j) r1_aa(c,i) t2_aaaa(a,b,l,j) l2_aaaa(l,k,b,c)
        //             += -1.00 P(i,j) r1_aa(c,i) t2_abab(a,b,j,l) l2_abab(k,l,c,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_bb(b,l) t1_aa(a,i) l1_bb(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(b,l) t1_aa(a,i) l1_aa(l,b)
        //             += -0.25 P(i,j) d_aa(j,k) r2_aaaa(b,c,m,l) t1_aa(a,i) l2_aaaa(m,l,b,c)
        //             += -0.25 P(i,j) d_aa(j,k) r2_abab(b,c,m,l) t1_aa(a,i) l2_abab(m,l,b,c)
        //             += -0.25 P(i,j) d_aa(j,k) r2_abab(b,c,l,m) t1_aa(a,i) l2_abab(l,m,b,c)
        //             += -0.25 P(i,j) d_aa(j,k) r2_abab(c,b,m,l) t1_aa(a,i) l2_abab(m,l,c,b)
        //             += -0.25 P(i,j) d_aa(j,k) r2_abab(c,b,l,m) t1_aa(a,i) l2_abab(l,m,c,b)
        //             += -0.25 P(i,j) d_aa(j,k) r2_bbbb(b,c,m,l) t1_aa(a,i) l2_bbbb(m,l,b,c)
        //             += +1.00 P(i,j) d_aa(j,k) r1_aa(a,l) t1_aa(b,i) l1_aa(l,b)
        //             += -1.00 P(i,j) d_aa(j,k) r2_abab(a,b,i,l) l1_bb(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r2_aaaa(a,b,l,i) l1_aa(l,b)
        //             += -0.50 P(i,j) r0 t1_aa(a,j) t2_aaaa(c,b,l,i) l2_aaaa(l,k,c,b)
        //             += -1.00 P(i,j) r1_aa(b,i) t1_aa(a,j) l1_aa(k,b)
        //             += -0.50 P(i,j) r2_abab(b,c,i,l) t1_aa(a,j) l2_abab(k,l,b,c)
        //             += -0.50 P(i,j) r2_abab(c,b,i,l) t1_aa(a,j) l2_abab(k,l,c,b)
        //             += -0.50 P(i,j) r2_aaaa(b,c,l,i) t1_aa(a,j) l2_aaaa(l,k,b,c)
        //             += -0.50 P(i,j) r0 t1_aa(a,j) t2_abab(c,b,i,l) l2_abab(k,l,c,b)
        //             += -0.50 P(i,j) r0 t1_aa(a,j) t2_abab(b,c,i,l) l2_abab(k,l,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t2_abab(a,b,l,m) t1_aa(c,i) l2_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t2_abab(a,b,m,l) t1_aa(c,i) l2_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t1_aa(a,m) t2_abab(c,b,i,l) l2_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t1_aa(a,m) t2_abab(b,c,i,l) l2_abab(m,l,b,c)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t2_aaaa(a,b,l,m) t1_aa(c,i) l2_aaaa(l,m,c,b)
        //             += -1.00 P(i,j) d_aa(j,k) r0 t2_abab(a,b,i,l) l1_bb(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r0 t2_aaaa(a,b,l,i) l1_aa(l,b)
        //             += +1.00 P(i,j) d_aa(j,k) r0 t1_aa(a,l) t1_aa(b,i) l1_aa(l,b)
        //             += +0.50 P(i,j) d_aa(j,k) r0 t1_aa(a,m) t2_aaaa(c,b,l,i) l2_aaaa(l,m,c,b)
        //             += +1.00 P(i,j) r0 t2_aaaa(a,b,l,i) t1_aa(c,j) l2_aaaa(l,k,c,b)
        //             += -1.00 P(i,j) d_aa(j,k) r1_aa(a,i) l0
        //             += +0.50 P(i,j) r1_aa(a,i) t2_abab(c,b,j,l) l2_abab(k,l,c,b)
        //             += +0.50 P(i,j) r1_aa(a,i) t2_abab(b,c,j,l) l2_abab(k,l,b,c)
        //             += +1.00 P(i,j) r1_aa(a,i) t1_aa(b,j) l1_aa(k,b)
        //             += +0.50 P(i,j) d_aa(j,k) r1_aa(c,i) t2_abab(a,b,l,m) l2_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_aa(j,k) r1_aa(c,i) t2_abab(a,b,m,l) l2_abab(m,l,c,b)
        D_oovo_aaaa("L,R,i,j,a,k") += tmps_["135_aaaa_LLvooo"]("R,L,a,i,k,j");
        D_oovo_aaaa("L,R,i,j,a,k") -= tmps_["135_aaaa_LLvooo"]("R,L,a,j,k,i");
        tmps_["135_aaaa_LLvooo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["136_bbbb_Lvoov"]("L,a,i,j,b")  = l2["bbbb"]("L,j,k,c,b") * t2["bbbb"]("a,c,k,i");

        // D_oovv_abab += -1.00 r2_abab(a,d,i,l) t2_bbbb(b,c,k,j) l2_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r2["abab"]("R,a,d,i,l") * tmps_["136_bbbb_Lvoov"]("L,b,j,l,d");

        // D_oovv_abab += -1.00 r0 t2_abab(a,c,i,k) t2_bbbb(b,d,l,j) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,c,i,k") * tmps_["136_bbbb_Lvoov"]("L,b,j,k,c") * r0("R");

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["137_bb_Lvv"]("L,a,b")  = l2["bbbb"]("L,i,j,c,b") * t2["bbbb"]("a,c,i,j");

        // D_oovv_abab += -0.50 r2_abab(a,d,i,j) t2_bbbb(b,c,k,l) l2_bbbb(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * r2["abab"]("R,a,d,i,j") * tmps_["137_bb_Lvv"]("L,b,d");

        // D_oovv_abab += -0.50 r0 t2_abab(a,c,i,j) t2_bbbb(b,d,k,l) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * t2["abab"]("a,c,i,j") * tmps_["137_bb_Lvv"]("L,b,c") * r0("R");

        // D_oovv_bbbb += -0.50 P(a,b) r2_bbbb(a,d,i,j) t2_bbbb(b,c,k,l) l2_bbbb(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * r2["bbbb"]("R,a,d,i,j") * tmps_["137_bb_Lvv"]("L,b,d");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += -0.50 r0 t2_bbbb(a,c,i,j) t2_bbbb(b,d,k,l) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") -= 0.50 * t2["bbbb"]("a,c,i,j") * tmps_["137_bb_Lvv"]("L,b,c") * r0("R");

        // flops: o2v0L1  = o3v2L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["138_bb_Loo"]("L,i,j")  = l2["bbbb"]("L,j,k,a,b") * t2["bbbb"]("a,b,k,i");

        // D_oovv_abab += -0.50 r2_abab(a,b,i,l) t2_bbbb(d,c,k,j) l2_bbbb(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * r2["abab"]("R,a,b,i,l") * tmps_["138_bb_Loo"]("L,j,l");

        // D_oovv_bbbb += -0.50 P(i,j) r2_bbbb(b,a,l,i) t2_bbbb(d,c,k,j) l2_bbbb(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * r2["bbbb"]("R,b,a,l,i") * tmps_["138_bb_Loo"]("L,j,l");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v1L2  = o2v2L2 o2v2L2 o1v1L2 o2v1L1 o2v1L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2 o2v0L1 o1v1L2 o1v1L2
        tmps_["139_bb_LLvo"]("L,R,a,i")  = -1.00 * l1["aa"]("L,k,c") * r2["abab"]("R,c,a,k,i");
        tmps_["139_bb_LLvo"]("L,R,a,i") += l1["bb"]("L,j,b") * r2["bbbb"]("R,a,b,j,i");
        tmps_["139_bb_LLvo"]("L,R,a,i") += l1["bb"]("L,j,b") * t1["bb"]("b,i") * r1["bb"]("R,a,j");

        // D_oovv_abab += +1.00 r2_bbbb(b,c,k,j) t1_aa(a,i) l1_bb(k,c)
        //             += -1.00 r2_abab(c,b,k,j) t1_aa(a,i) l1_aa(k,c)
        //             += +1.00 r1_bb(b,k) t1_aa(a,i) t1_bb(c,j) l1_bb(k,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1["aa"]("a,i") * tmps_["139_bb_LLvo"]("L,R,b,j");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r2_bbbb(a,c,k,i) t1_bb(b,j) l1_bb(k,c)
        //             += -1.00 P(i,j) P(a,b) r2_abab(c,a,k,i) t1_bb(b,j) l1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_bb(a,k) t1_bb(b,j) t1_bb(c,i) l1_bb(k,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("b,j") * tmps_["139_bb_LLvo"]("L,R,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["140_bbbb_Lvooo"]("L,a,i,j,k")  = t1["bb"]("b,k") * tmps_["24_bbbb_Lvoov"]("L,a,i,j,b");

        // D_oovv_bbbb += -1.00 P(i,j) P(a,b) r1_bb(a,l) t2_abab(c,b,k,i) t1_bb(d,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = r1["bb"]("R,a,l") * tmps_["140_bbbb_Lvooo"]("L,b,i,l,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r0 t2_abab(c,a,k,i) t1_bb(b,l) t1_bb(d,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("b,l") * tmps_["140_bbbb_Lvooo"]("L,a,i,l,j") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["141_bbbb_Lvooo"]("L,a,i,j,k")  = t1["bb"]("b,k") * tmps_["86_bbbb_Lvoov"]("L,a,i,j,b");
        tmps_["86_bbbb_Lvoov"].~TArrayD();

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r0 t2_bbbb(a,c,k,i) t1_bb(b,l) t1_bb(d,j) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("b,l") * tmps_["141_bbbb_Lvooo"]("L,a,i,l,j") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v1L2  = o1v2L2 o2v2L2 o2v2L2 o1v1L2 o1v1L2 o2v1L2 o1v1L2 o1v2L2 o1v1L2 o2v1L2 o1v1L2 o2v1L2 o1v1L2 o2v2L2 o1v1L2 o2v2L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2
        tmps_["142_bb_LLvo"]("L,R,a,i")  = -0.50 * tmps_["137_bb_Lvv"]("L,a,c") * r1["bb"]("R,c,i");
        tmps_["142_bb_LLvo"]("L,R,a,i") -= tmps_["136_bbbb_Lvoov"]("L,a,i,k,c") * r1["bb"]("R,c,k");
        tmps_["142_bb_LLvo"]("L,R,a,i") += t2["abab"]("e,a,m,i") * tmps_["55_aa_LLov"]("L,R,m,e");
        tmps_["142_bb_LLvo"]("L,R,a,i") -= 0.50 * tmps_["138_bb_Loo"]("L,i,k") * r1["bb"]("R,a,k");
        tmps_["142_bb_LLvo"]("L,R,a,i") += tmps_["91_bb_LLvv"]("L,R,d,a") * t1["bb"]("d,i");
        tmps_["142_bb_LLvo"]("L,R,a,i") += t1["bb"]("a,l") * tmps_["51_bb_LLoo"]("L,R,l,i");
        tmps_["142_bb_LLvo"]("L,R,a,i") += r1["bb"]("R,a,k") * tmps_["48_bb_Loo"]("L,i,k");
        tmps_["142_bb_LLvo"]("L,R,a,i") += r1["aa"]("R,b,j") * tmps_["43_bbaa_Lvoov"]("L,a,i,j,b");
        tmps_["142_bb_LLvo"]("L,R,a,i") -= r1["bb"]("R,c,k") * tmps_["24_bbbb_Lvoov"]("L,a,i,k,c");
        tmps_["138_bb_Loo"].~TArrayD();
        tmps_["137_bb_Lvv"].~TArrayD();
        tmps_["136_bbbb_Lvoov"].~TArrayD();
        tmps_["91_bb_LLvv"].~TArrayD();
        tmps_["55_aa_LLov"].~TArrayD();
        tmps_["43_bbaa_Lvoov"].~TArrayD();
        tmps_["24_bbbb_Lvoov"].~TArrayD();

        // D_oovv_abab += +1.00 r1_aa(d,l) t1_aa(a,i) t2_abab(c,b,k,j) l2_aaaa(l,k,c,d)
        //             += -1.00 r1_bb(d,l) t1_aa(a,i) t2_bbbb(b,c,k,j) l2_bbbb(l,k,c,d)
        //             += -0.50 r1_bb(d,j) t1_aa(a,i) t2_bbbb(b,c,k,l) l2_bbbb(k,l,c,d)
        //             += -0.50 r1_bb(b,l) t1_aa(a,i) t2_bbbb(d,c,k,j) l2_bbbb(l,k,d,c)
        //             += +0.50 r2_abab(d,b,l,k) t1_aa(a,i) t1_bb(c,j) l2_abab(l,k,d,c)
        //             += +0.50 r2_abab(d,b,k,l) t1_aa(a,i) t1_bb(c,j) l2_abab(k,l,d,c)
        //             += +0.50 r2_bbbb(b,d,l,k) t1_aa(a,i) t1_bb(c,j) l2_bbbb(l,k,c,d)
        //             += +1.00 r1_bb(c,j) t1_aa(a,i) t1_bb(b,k) l1_bb(k,c)
        //             += +0.50 r2_abab(c,d,l,j) t1_aa(a,i) t1_bb(b,k) l2_abab(l,k,c,d)
        //             += +0.50 r2_abab(d,c,l,j) t1_aa(a,i) t1_bb(b,k) l2_abab(l,k,d,c)
        //             += +0.50 r2_bbbb(c,d,l,j) t1_aa(a,i) t1_bb(b,k) l2_bbbb(l,k,c,d)
        //             += +0.50 r1_bb(b,l) t1_aa(a,i) t2_abab(d,c,k,j) l2_abab(k,l,d,c)
        //             += +0.50 r1_bb(b,l) t1_aa(a,i) t2_abab(c,d,k,j) l2_abab(k,l,c,d)
        //             += +1.00 r1_aa(d,l) t1_aa(a,i) t2_bbbb(b,c,k,j) l2_abab(l,k,d,c)
        //             += -1.00 r1_bb(d,l) t1_aa(a,i) t2_abab(c,b,k,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1["aa"]("a,i") * tmps_["142_bb_LLvo"]("L,R,b,j");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r1_aa(d,l) t2_abab(c,a,k,i) t1_bb(b,j) l2_aaaa(l,k,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_bb(d,l) t2_bbbb(a,c,k,i) t1_bb(b,j) l2_bbbb(l,k,c,d)
        //             += -0.50 P(i,j) P(a,b) r1_bb(d,i) t2_bbbb(a,c,k,l) t1_bb(b,j) l2_bbbb(k,l,c,d)
        //             += -0.50 P(i,j) P(a,b) r1_bb(a,l) t1_bb(b,j) t2_bbbb(d,c,k,i) l2_bbbb(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r2_abab(d,a,l,k) t1_bb(b,j) t1_bb(c,i) l2_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r2_abab(d,a,k,l) t1_bb(b,j) t1_bb(c,i) l2_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r2_bbbb(a,d,l,k) t1_bb(b,j) t1_bb(c,i) l2_bbbb(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_bb(c,i) t1_bb(a,k) t1_bb(b,j) l1_bb(k,c)
        //             += +0.50 P(i,j) P(a,b) r2_abab(c,d,l,i) t1_bb(a,k) t1_bb(b,j) l2_abab(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_abab(d,c,l,i) t1_bb(a,k) t1_bb(b,j) l2_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r2_bbbb(c,d,l,i) t1_bb(a,k) t1_bb(b,j) l2_bbbb(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,l) t1_bb(b,j) t2_abab(d,c,k,i) l2_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,l) t1_bb(b,j) t2_abab(c,d,k,i) l2_abab(k,l,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_aa(d,l) t2_bbbb(a,c,k,i) t1_bb(b,j) l2_abab(l,k,d,c)
        //             += -1.00 P(i,j) P(a,b) r1_bb(d,l) t2_abab(c,a,k,i) t1_bb(b,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("b,j") * tmps_["142_bb_LLvo"]("L,R,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v1L1  = o2v2L1 o2v2L1 o1v1L1 o2v1L1 o2v1L1 o1v1L1 o1v2L1 o1v1L1 o2v1L1 o1v1L1 o2v1L1 o1v1L1 o1v2L1 o1v1L1
        //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o2v0L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
        tmps_["143_bb_Lov"]("L,i,a")  = -2.00 * l1["aa"]("L,k,c") * t2["abab"]("c,a,k,i");
        tmps_["143_bb_Lov"]("L,i,a") += 2.00 * l1["bb"]("L,l,d") * t2["bbbb"]("a,d,l,i");
        tmps_["143_bb_Lov"]("L,i,a") += 2.00 * l1["bb"]("L,l,d") * t1["bb"]("d,i") * t1["bb"]("a,l");
        tmps_["143_bb_Lov"]("L,i,a") += t1["bb"]("b,i") * tmps_["87_bb_Lvv"]("L,b,a");
        tmps_["143_bb_Lov"]("L,i,a") += t1["bb"]("a,j") * tmps_["49_bb_Loo"]("L,i,j");
        tmps_["143_bb_Lov"]("L,i,a") += 2.00 * t1["bb"]("a,j") * tmps_["48_bb_Loo"]("L,i,j");
        tmps_["143_bb_Lov"]("L,i,a") += 2.00 * t1["bb"]("b,i") * tmps_["30_bb_Lvv"]("L,b,a");
        tmps_["87_bb_Lvv"].~TArrayD();

        // D_oovv_bbbb += +0.50 P(i,j) P(a,b) r0 t2_bbbb(a,c,k,l) t1_bb(b,j) t1_bb(d,i) l2_bbbb(k,l,d,c)
        //             += -1.00 P(i,j) P(a,b) r0 t2_abab(c,a,k,i) t1_bb(b,j) l1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r0 t2_bbbb(a,c,k,i) t1_bb(b,j) l1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r0 t1_bb(a,k) t1_bb(b,j) t1_bb(c,i) l1_bb(k,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_bb(a,l) t1_bb(b,j) t2_bbbb(d,c,k,i) l2_bbbb(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_bb(a,l) t1_bb(b,j) t2_abab(d,c,k,i) l2_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_bb(a,l) t1_bb(b,j) t2_abab(c,d,k,i) l2_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r0 t2_abab(c,a,k,l) t1_bb(b,j) t1_bb(d,i) l2_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r0 t2_abab(c,a,l,k) t1_bb(b,j) t1_bb(d,i) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["90_bb_Lvo"]("R,b,j") * tmps_["143_bb_Lov"]("L,i,a");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += +0.50 P(i,j) P(a,b) r1_bb(a,i) t2_bbbb(b,c,k,l) t1_bb(d,j) l2_bbbb(k,l,d,c)
        //             += -1.00 P(i,j) P(a,b) r1_bb(a,i) t2_abab(c,b,k,j) l1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_bb(a,i) t2_bbbb(b,c,k,j) l1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_bb(a,i) t1_bb(b,k) t1_bb(c,j) l1_bb(k,c)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,i) t1_bb(b,l) t2_bbbb(d,c,k,j) l2_bbbb(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,i) t1_bb(b,l) t2_abab(d,c,k,j) l2_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,i) t1_bb(b,l) t2_abab(c,d,k,j) l2_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,i) t2_abab(c,b,k,l) t1_bb(d,j) l2_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_bb(a,i) t2_abab(c,b,l,k) t1_bb(d,j) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * r1["bb"]("R,a,i") * tmps_["143_bb_Lov"]("L,j,b");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r0 t1_aa(a,i) t2_bbbb(b,c,k,l) t1_bb(d,j) l2_bbbb(k,l,d,c)
        //             += -1.00 r0 t1_aa(a,i) t2_abab(c,b,k,j) l1_aa(k,c)
        //             += +1.00 r0 t1_aa(a,i) t2_bbbb(b,c,k,j) l1_bb(k,c)
        //             += +1.00 r0 t1_aa(a,i) t1_bb(b,k) t1_bb(c,j) l1_bb(k,c)
        //             += +0.50 r0 t1_aa(a,i) t1_bb(b,l) t2_bbbb(d,c,k,j) l2_bbbb(k,l,d,c)
        //             += +0.50 r0 t1_aa(a,i) t1_bb(b,l) t2_abab(d,c,k,j) l2_abab(k,l,d,c)
        //             += +0.50 r0 t1_aa(a,i) t1_bb(b,l) t2_abab(c,d,k,j) l2_abab(k,l,c,d)
        //             += +0.50 r0 t1_aa(a,i) t2_abab(c,b,k,l) t1_bb(d,j) l2_abab(k,l,c,d)
        //             += +0.50 r0 t1_aa(a,i) t2_abab(c,b,l,k) t1_bb(d,j) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * tmps_["109_aa_Lvo"]("R,a,i") * tmps_["143_bb_Lov"]("L,j,b");
        tmps_["109_aa_Lvo"].~TArrayD();

        // D_oovv_abab += +0.50 r1_aa(a,i) t2_bbbb(b,c,k,l) t1_bb(d,j) l2_bbbb(k,l,d,c)
        //             += -1.00 r1_aa(a,i) t2_abab(c,b,k,j) l1_aa(k,c)
        //             += +1.00 r1_aa(a,i) t2_bbbb(b,c,k,j) l1_bb(k,c)
        //             += +1.00 r1_aa(a,i) t1_bb(b,k) t1_bb(c,j) l1_bb(k,c)
        //             += +0.50 r1_aa(a,i) t1_bb(b,l) t2_bbbb(d,c,k,j) l2_bbbb(k,l,d,c)
        //             += +0.50 r1_aa(a,i) t1_bb(b,l) t2_abab(d,c,k,j) l2_abab(k,l,d,c)
        //             += +0.50 r1_aa(a,i) t1_bb(b,l) t2_abab(c,d,k,j) l2_abab(k,l,c,d)
        //             += +0.50 r1_aa(a,i) t2_abab(c,b,k,l) t1_bb(d,j) l2_abab(k,l,c,d)
        //             += +0.50 r1_aa(a,i) t2_abab(c,b,l,k) t1_bb(d,j) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * r1["aa"]("R,a,i") * tmps_["143_bb_Lov"]("L,j,b");

        // flops: o1v1L2  = o2v1L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["144_bb_LLvo"]("L,R,a,i")  = t1["bb"]("a,j") * tmps_["53_bb_LLoo"]("L,R,j,i");

        // D_oovv_abab += -1.00 r1_bb(d,l) t1_aa(a,i) t1_bb(b,k) t1_bb(c,j) l2_bbbb(l,k,c,d)
        //             += +1.00 r1_aa(d,l) t1_aa(a,i) t1_bb(b,k) t1_bb(c,j) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["aa"]("a,i") * tmps_["144_bb_LLvo"]("L,R,b,j");

        // D_oovv_bbbb += -1.00 P(i,j) P(a,b) r1_bb(d,l) t1_bb(a,k) t1_bb(b,j) t1_bb(c,i) l2_bbbb(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_aa(d,l) t1_bb(a,k) t1_bb(b,j) t1_bb(c,i) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("b,j") * tmps_["144_bb_LLvo"]("L,R,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v1L2  = o1v2L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["145_bb_LLov"]("R,L,i,a")  = r1["bb"]("R,b,i") * tmps_["30_bb_Lvv"]("L,b,a");
        tmps_["30_bb_Lvv"].~TArrayD();

        // D_oovv_abab += +0.50 r1_bb(d,j) t1_aa(a,i) t2_abab(c,b,k,l) l2_abab(k,l,c,d)
        //             += +0.50 r1_bb(d,j) t1_aa(a,i) t2_abab(c,b,l,k) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1["aa"]("a,i") * tmps_["145_bb_LLov"]("R,L,j,b");

        // D_oovv_bbbb += +0.50 P(i,j) P(a,b) r1_bb(d,i) t2_abab(c,a,k,l) t1_bb(b,j) l2_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_bb(d,i) t2_abab(c,a,l,k) t1_bb(b,j) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("b,j") * tmps_["145_bb_LLov"]("R,L,i,a");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o3v1L2  = o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        //  mems: o3v1L2  = o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        tmps_["146_bbbb_LLooov"]("R,L,i,j,k,a")  = -1.00 * Id["bb_oo"]("i,j") * tmps_["144_bb_LLvo"]("L,R,a,k");
        tmps_["146_bbbb_LLooov"]("R,L,i,j,k,a") += r1["bb"]("R,a,k") * tmps_["48_bb_Loo"]("L,i,j");
        tmps_["146_bbbb_LLooov"]("R,L,i,j,k,a") += t1["bb"]("a,i") * tmps_["53_bb_LLoo"]("L,R,j,k");
        tmps_["146_bbbb_LLooov"]("R,L,i,j,k,a") += 0.50 * tmps_["143_bb_Lov"]("L,k,a") * tmps_["52_bb_Loo"]("R,i,j");
        tmps_["146_bbbb_LLooov"]("R,L,i,j,k,a") += tmps_["142_bb_LLvo"]("L,R,a,k") * Id["bb_oo"]("i,j");
        tmps_["146_bbbb_LLooov"]("R,L,i,j,k,a") -= t1["bb"]("a,i") * tmps_["51_bb_LLoo"]("L,R,j,k");
        tmps_["146_bbbb_LLooov"]("R,L,i,j,k,a") += tmps_["139_bb_LLvo"]("L,R,a,k") * Id["bb_oo"]("i,j");
        tmps_["146_bbbb_LLooov"]("R,L,i,j,k,a") -= tmps_["124_bb_Loo"]("L,i,j") * tmps_["90_bb_Lvo"]("R,a,k");
        tmps_["146_bbbb_LLooov"]("R,L,i,j,k,a") -= 0.50 * tmps_["49_bb_Loo"]("L,k,j") * tmps_["90_bb_Lvo"]("R,a,i");
        tmps_["146_bbbb_LLooov"]("R,L,i,j,k,a") -= tmps_["50_bb_Loo"]("L,k,j") * tmps_["90_bb_Lvo"]("R,a,i");
        tmps_["146_bbbb_LLooov"]("R,L,i,j,k,a") -= Id["bb_oo"]("i,j") * t1["bb"]("a,k") * tmps_["126_LL"]("L,R");
        tmps_["146_bbbb_LLooov"]("R,L,i,j,k,a") -= tmps_["48_bb_Loo"]("L,k,j") * tmps_["90_bb_Lvo"]("R,a,i");
        tmps_["146_bbbb_LLooov"]("R,L,i,j,k,a") -= r1["bb"]("R,a,k") * tmps_["124_bb_Loo"]("L,i,j");
        tmps_["146_bbbb_LLooov"]("R,L,i,j,k,a") += r0("R") * tmps_["141_bbbb_Lvooo"]("L,a,k,j,i");
        tmps_["146_bbbb_LLooov"]("R,L,i,j,k,a") += 0.50 * r1["bb"]("R,a,k") * tmps_["49_bb_Loo"]("L,i,j");
        tmps_["146_bbbb_LLooov"]("R,L,i,j,k,a") += r1["bb"]("R,a,k") * tmps_["50_bb_Loo"]("L,i,j");
        tmps_["146_bbbb_LLooov"]("R,L,i,j,k,a") += r0("R") * tmps_["140_bbbb_Lvooo"]("L,a,k,j,i");
        tmps_["146_bbbb_LLooov"]("R,L,i,j,k,a") += Id["bb_oo"]("i,j") * tmps_["145_bb_LLov"]("R,L,k,a");
        tmps_["145_bb_LLov"].~TArrayD();
        tmps_["144_bb_LLvo"].~TArrayD();
        tmps_["143_bb_Lov"].~TArrayD();
        tmps_["142_bb_LLvo"].~TArrayD();
        tmps_["141_bbbb_Lvooo"].~TArrayD();
        tmps_["140_bbbb_Lvooo"].~TArrayD();
        tmps_["139_bb_LLvo"].~TArrayD();
        tmps_["126_LL"].~TArrayD();
        tmps_["124_bb_Loo"].~TArrayD();
        tmps_["90_bb_Lvo"].~TArrayD();
        tmps_["53_bb_LLoo"].~TArrayD();
        tmps_["52_bb_Loo"].~TArrayD();
        tmps_["51_bb_LLoo"].~TArrayD();
        tmps_["50_bb_Loo"].~TArrayD();
        tmps_["49_bb_Loo"].~TArrayD();
        tmps_["48_bb_Loo"].~TArrayD();

        // D_ooov_bbbb += -0.50 P(i,j) d_bb(j,k) r1_bb(c,i) t2_abab(b,a,l,m) l2_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r1_bb(c,i) t2_abab(b,a,m,l) l2_abab(m,l,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t1_bb(a,l) t1_bb(b,i) l2_bbbb(m,l,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t1_bb(a,l) t1_bb(b,i) l2_abab(m,l,c,b)
        //             += -0.50 P(i,j) r1_bb(a,i) t2_abab(c,b,l,j) l2_abab(l,k,c,b)
        //             += -0.50 P(i,j) r1_bb(a,i) t2_abab(b,c,l,j) l2_abab(l,k,b,c)
        //             += -1.00 P(i,j) r1_bb(c,l) t1_bb(a,j) t1_bb(b,i) l2_bbbb(l,k,b,c)
        //             += +1.00 P(i,j) r1_aa(c,l) t1_bb(a,j) t1_bb(b,i) l2_abab(l,k,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_bbbb(a,b,l,m) t1_bb(c,i) l2_bbbb(l,m,c,b)
        //             += +1.00 P(i,j) d_bb(j,k) r0 t2_abab(b,a,l,i) l1_aa(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r0 t2_bbbb(a,b,l,i) l1_bb(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r0 t1_bb(a,l) t1_bb(b,i) l1_bb(l,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t1_bb(a,m) t2_bbbb(c,b,l,i) l2_bbbb(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t1_bb(a,m) t2_abab(c,b,l,i) l2_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t1_bb(a,m) t2_abab(b,c,l,i) l2_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_abab(b,a,l,m) t1_bb(c,i) l2_abab(l,m,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r0 t2_abab(b,a,m,l) t1_bb(c,i) l2_abab(m,l,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t2_abab(b,a,l,i) l2_aaaa(m,l,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t2_bbbb(a,b,l,i) l2_bbbb(m,l,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r1_bb(c,i) t2_bbbb(a,b,l,m) l2_bbbb(l,m,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r1_bb(a,m) t2_bbbb(c,b,l,i) l2_bbbb(m,l,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r2_abab(c,a,m,l) t1_bb(b,i) l2_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r2_abab(c,a,l,m) t1_bb(b,i) l2_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r2_bbbb(a,c,m,l) t1_bb(b,i) l2_bbbb(m,l,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(b,i) t1_bb(a,l) l1_bb(l,b)
        //             += -0.50 P(i,j) d_bb(j,k) r2_abab(b,c,m,i) t1_bb(a,l) l2_abab(m,l,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r2_abab(c,b,m,i) t1_bb(a,l) l2_abab(m,l,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r2_bbbb(b,c,m,i) t1_bb(a,l) l2_bbbb(m,l,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r1_bb(a,m) t2_abab(c,b,l,i) l2_abab(l,m,c,b)
        //             += -0.50 P(i,j) d_bb(j,k) r1_bb(a,m) t2_abab(b,c,l,i) l2_abab(l,m,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t2_bbbb(a,b,l,i) l2_abab(m,l,c,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t2_abab(b,a,l,i) l2_abab(l,m,b,c)
        //             += +1.00 P(i,j) r1_bb(b,i) t1_bb(a,j) l1_bb(k,b)
        //             += +0.50 P(i,j) r2_abab(b,c,l,i) t1_bb(a,j) l2_abab(l,k,b,c)
        //             += +0.50 P(i,j) r2_abab(c,b,l,i) t1_bb(a,j) l2_abab(l,k,c,b)
        //             += +0.50 P(i,j) r2_bbbb(b,c,l,i) t1_bb(a,j) l2_bbbb(l,k,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r2_bbbb(a,b,l,i) l1_bb(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r2_abab(b,a,l,i) l1_aa(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(a,l) t1_bb(b,i) l1_bb(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r0 t1_bb(a,i) l0
        //             += +0.50 P(i,j) r0 t1_bb(a,j) t2_bbbb(c,b,l,i) l2_bbbb(l,k,c,b)
        //             += +1.00 P(i,j) r0 t1_bb(a,j) t1_bb(b,i) l1_bb(k,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(b,l) t1_bb(a,i) l1_bb(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_aa(b,l) t1_bb(a,i) l1_aa(l,b)
        //             += +0.25 P(i,j) d_bb(j,k) r2_aaaa(b,c,m,l) t1_bb(a,i) l2_aaaa(m,l,b,c)
        //             += +0.25 P(i,j) d_bb(j,k) r2_abab(b,c,m,l) t1_bb(a,i) l2_abab(m,l,b,c)
        //             += +0.25 P(i,j) d_bb(j,k) r2_abab(b,c,l,m) t1_bb(a,i) l2_abab(l,m,b,c)
        //             += +0.25 P(i,j) d_bb(j,k) r2_abab(c,b,m,l) t1_bb(a,i) l2_abab(m,l,c,b)
        //             += +0.25 P(i,j) d_bb(j,k) r2_abab(c,b,l,m) t1_bb(a,i) l2_abab(l,m,c,b)
        //             += +0.25 P(i,j) d_bb(j,k) r2_bbbb(b,c,m,l) t1_bb(a,i) l2_bbbb(m,l,b,c)
        //             += +0.50 P(i,j) r0 t1_bb(a,j) t2_abab(c,b,l,i) l2_abab(l,k,c,b)
        //             += +0.50 P(i,j) r0 t1_bb(a,j) t2_abab(b,c,l,i) l2_abab(l,k,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(a,i) l0
        //             += -1.00 P(i,j) r0 t2_bbbb(a,b,l,i) t1_bb(c,j) l2_bbbb(l,k,c,b)
        //             += -0.50 P(i,j) r1_bb(a,i) t2_bbbb(c,b,l,j) l2_bbbb(l,k,c,b)
        //             += -1.00 P(i,j) r1_bb(a,i) t1_bb(b,j) l1_bb(k,b)
        //             += -1.00 P(i,j) r0 t2_abab(b,a,l,i) t1_bb(c,j) l2_abab(l,k,b,c)
        D_ooov_bbbb("L,R,i,j,k,a") -= tmps_["146_bbbb_LLooov"]("R,L,j,k,i,a");
        D_ooov_bbbb("L,R,i,j,k,a") += tmps_["146_bbbb_LLooov"]("R,L,i,k,j,a");

        // D_oovo_bbbb += +0.50 P(i,j) d_bb(j,k) r1_bb(c,i) t2_abab(b,a,l,m) l2_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r1_bb(c,i) t2_abab(b,a,m,l) l2_abab(m,l,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t1_bb(a,l) t1_bb(b,i) l2_bbbb(m,l,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t1_bb(a,l) t1_bb(b,i) l2_abab(m,l,c,b)
        //             += +0.50 P(i,j) r1_bb(a,i) t2_abab(c,b,l,j) l2_abab(l,k,c,b)
        //             += +0.50 P(i,j) r1_bb(a,i) t2_abab(b,c,l,j) l2_abab(l,k,b,c)
        //             += +1.00 P(i,j) r1_bb(c,l) t1_bb(a,j) t1_bb(b,i) l2_bbbb(l,k,b,c)
        //             += -1.00 P(i,j) r1_aa(c,l) t1_bb(a,j) t1_bb(b,i) l2_abab(l,k,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t2_bbbb(a,b,l,m) t1_bb(c,i) l2_bbbb(l,m,c,b)
        //             += -1.00 P(i,j) d_bb(j,k) r0 t2_abab(b,a,l,i) l1_aa(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r0 t2_bbbb(a,b,l,i) l1_bb(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r0 t1_bb(a,l) t1_bb(b,i) l1_bb(l,b)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t1_bb(a,m) t2_bbbb(c,b,l,i) l2_bbbb(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t1_bb(a,m) t2_abab(c,b,l,i) l2_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t1_bb(a,m) t2_abab(b,c,l,i) l2_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t2_abab(b,a,l,m) t1_bb(c,i) l2_abab(l,m,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r0 t2_abab(b,a,m,l) t1_bb(c,i) l2_abab(m,l,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t2_abab(b,a,l,i) l2_aaaa(m,l,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t2_bbbb(a,b,l,i) l2_bbbb(m,l,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r1_bb(c,i) t2_bbbb(a,b,l,m) l2_bbbb(l,m,b,c)
        //             += -0.50 P(i,j) d_bb(j,k) r1_bb(a,m) t2_bbbb(c,b,l,i) l2_bbbb(m,l,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r2_abab(c,a,m,l) t1_bb(b,i) l2_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r2_abab(c,a,l,m) t1_bb(b,i) l2_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r2_bbbb(a,c,m,l) t1_bb(b,i) l2_bbbb(m,l,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(b,i) t1_bb(a,l) l1_bb(l,b)
        //             += +0.50 P(i,j) d_bb(j,k) r2_abab(b,c,m,i) t1_bb(a,l) l2_abab(m,l,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r2_abab(c,b,m,i) t1_bb(a,l) l2_abab(m,l,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r2_bbbb(b,c,m,i) t1_bb(a,l) l2_bbbb(m,l,b,c)
        //             += +0.50 P(i,j) d_bb(j,k) r1_bb(a,m) t2_abab(c,b,l,i) l2_abab(l,m,c,b)
        //             += +0.50 P(i,j) d_bb(j,k) r1_bb(a,m) t2_abab(b,c,l,i) l2_abab(l,m,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r1_aa(c,m) t2_bbbb(a,b,l,i) l2_abab(m,l,c,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(c,m) t2_abab(b,a,l,i) l2_abab(l,m,b,c)
        //             += -1.00 P(i,j) r1_bb(b,i) t1_bb(a,j) l1_bb(k,b)
        //             += -0.50 P(i,j) r2_abab(b,c,l,i) t1_bb(a,j) l2_abab(l,k,b,c)
        //             += -0.50 P(i,j) r2_abab(c,b,l,i) t1_bb(a,j) l2_abab(l,k,c,b)
        //             += -0.50 P(i,j) r2_bbbb(b,c,l,i) t1_bb(a,j) l2_bbbb(l,k,b,c)
        //             += +1.00 P(i,j) d_bb(j,k) r2_bbbb(a,b,l,i) l1_bb(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r2_abab(b,a,l,i) l1_aa(l,b)
        //             += +1.00 P(i,j) d_bb(j,k) r1_bb(a,l) t1_bb(b,i) l1_bb(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r0 t1_bb(a,i) l0
        //             += -0.50 P(i,j) r0 t1_bb(a,j) t2_bbbb(c,b,l,i) l2_bbbb(l,k,c,b)
        //             += -1.00 P(i,j) r0 t1_bb(a,j) t1_bb(b,i) l1_bb(k,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(b,l) t1_bb(a,i) l1_bb(l,b)
        //             += -1.00 P(i,j) d_bb(j,k) r1_aa(b,l) t1_bb(a,i) l1_aa(l,b)
        //             += -0.25 P(i,j) d_bb(j,k) r2_aaaa(b,c,m,l) t1_bb(a,i) l2_aaaa(m,l,b,c)
        //             += -0.25 P(i,j) d_bb(j,k) r2_abab(b,c,m,l) t1_bb(a,i) l2_abab(m,l,b,c)
        //             += -0.25 P(i,j) d_bb(j,k) r2_abab(b,c,l,m) t1_bb(a,i) l2_abab(l,m,b,c)
        //             += -0.25 P(i,j) d_bb(j,k) r2_abab(c,b,m,l) t1_bb(a,i) l2_abab(m,l,c,b)
        //             += -0.25 P(i,j) d_bb(j,k) r2_abab(c,b,l,m) t1_bb(a,i) l2_abab(l,m,c,b)
        //             += -0.25 P(i,j) d_bb(j,k) r2_bbbb(b,c,m,l) t1_bb(a,i) l2_bbbb(m,l,b,c)
        //             += -0.50 P(i,j) r0 t1_bb(a,j) t2_abab(c,b,l,i) l2_abab(l,k,c,b)
        //             += -0.50 P(i,j) r0 t1_bb(a,j) t2_abab(b,c,l,i) l2_abab(l,k,b,c)
        //             += -1.00 P(i,j) d_bb(j,k) r1_bb(a,i) l0
        //             += +1.00 P(i,j) r0 t2_bbbb(a,b,l,i) t1_bb(c,j) l2_bbbb(l,k,c,b)
        //             += +0.50 P(i,j) r1_bb(a,i) t2_bbbb(c,b,l,j) l2_bbbb(l,k,c,b)
        //             += +1.00 P(i,j) r1_bb(a,i) t1_bb(b,j) l1_bb(k,b)
        //             += +1.00 P(i,j) r0 t2_abab(b,a,l,i) t1_bb(c,j) l2_abab(l,k,b,c)
        D_oovo_bbbb("L,R,i,j,a,k") += tmps_["146_bbbb_LLooov"]("R,L,j,k,i,a");
        D_oovo_bbbb("L,R,i,j,a,k") -= tmps_["146_bbbb_LLooov"]("R,L,i,k,j,a");
        tmps_["146_bbbb_LLooov"].~TArrayD();
    }
} // hilbert
#endif