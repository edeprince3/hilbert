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


    void EOM_EE_RDM::rdm2_00_1() {

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

                // D_oovv_aaaa  = +1.00 r2_aaaa(b,a,i,j) l0
        // flops: o2v2L2  = o2v2L2
        //  mems: o2v2L2  = o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b")  = l0("L") * r2["aaaa"]("R,b,a,i,j");

        // D_oovv_bbbb  = +1.00 r2_bbbb(b,a,i,j) l0
        // flops: o2v2L2  = o2v2L2
        //  mems: o2v2L2  = o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b")  = l0("L") * r2["bbbb"]("R,b,a,i,j");

        // D_vvvv_abab  = -1.00 r1_bb(d,j) t1_aa(c,i) l2_abab(i,j,a,b)
        // flops: o0v4L2  = o2v3L2 o1v4L2
        //  mems: o0v4L2  = o1v3L2 o0v4L2
        D_vvvv_abab("L,R,a,b,c,d")  = -1.00 * l2["abab"]("L,i,j,a,b") * r1["bb"]("R,d,j") * t1["aa"]("c,i");

        // D_oovv_abab  = -1.00 r2_abab(a,b,i,j) l0
        // flops: o2v2L2  = o2v2L2
        //  mems: o2v2L2  = o2v2L2
        D_oovv_abab("L,R,i,j,a,b")  = -1.00 * l0("L") * r2["abab"]("R,a,b,i,j");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r1_aa(d,i) t2_aaaa(a,c,k,j) t1_aa(b,l) l2_aaaa(k,l,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_aa(d,i) t2_abab(a,c,j,k) t1_aa(b,l) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o3v3L1 o3v3L1 o2v2L1 o3v2L2 o3v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L1 o2v2L1 o3v1L2 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = (t2["aaaa"]("a,c,k,j") * l2["aaaa"]("L,k,l,c,d") + -1.00 * t2["abab"]("a,c,j,k") * l2["abab"]("L,l,k,d,c")) * r1["aa"]("R,d,i") * t1["aa"]("b,l");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +1.00 r2_abab(c,b,i,j) t1_aa(a,k) l1_aa(k,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += l1["aa"]("L,k,c") * r2["abab"]("R,c,b,i,j") * t1["aa"]("a,k");

        // D_oovv_abab += -1.00 r2_abab(d,b,l,j) t2_aaaa(a,c,k,i) l2_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v3L1 o3v3L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= l2["aaaa"]("L,l,k,c,d") * t2["aaaa"]("a,c,k,i") * r2["abab"]("R,d,b,l,j");

        // D_oovv_abab += +1.00 r1_aa(d,i) t1_aa(a,l) t2_abab(c,b,k,j) l2_aaaa(k,l,c,d)
        // flops: o2v2L2 += o3v3L1 o3v2L2 o3v2L2
        //  mems: o2v2L2 += o2v2L1 o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += l2["aaaa"]("L,k,l,c,d") * t2["abab"]("c,b,k,j") * r1["aa"]("R,d,i") * t1["aa"]("a,l");

        // D_oovv_bbbb += -1.00 P(a,b) r1_bb(d,l) t2_bbbb(a,c,i,j) t1_bb(b,k) l2_bbbb(l,k,c,d)
        //             += +1.00 P(a,b) r1_aa(d,l) t2_bbbb(a,c,i,j) t1_bb(b,k) l2_abab(l,k,d,c)
        //             += +1.00 P(a,b) r2_bbbb(a,c,i,j) t1_bb(b,k) l1_bb(k,c)
        //             += +1.00 P(a,b) r1_bb(a,l) t1_bb(b,k) t1_bb(c,i) t1_bb(d,j) l2_bbbb(l,k,d,c)
        //             += -0.50 P(a,b) r1_bb(a,l) t1_bb(b,k) t2_bbbb(d,c,i,j) l2_bbbb(l,k,d,c)
        // flops: o2v2L2 += o2v2L2 o2v2L2 o1v1L2 o3v2L2 o3v2L2 o3v1L2 o3v2L1 o4v1L1 o4v1L2 o3v1L2 o4v2L1 o4v1L2 o3v1L2 o3v2L2
        //  mems: o2v2L2 += o1v1L2 o1v1L2 o1v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L1 o4v0L1 o3v1L2 o3v1L2 o4v0L1 o3v1L2 o3v1L2 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = (t2["bbbb"]("a,c,i,j") * (l2["bbbb"]("L,l,k,c,d") * r1["bb"]("R,d,l") + -1.00 * l2["abab"]("L,l,k,d,c") * r1["aa"]("R,d,l")) + -1.00 * l1["bb"]("L,k,c") * r2["bbbb"]("R,a,c,i,j") + -1.00 * t1["bb"]("c,i") * l2["bbbb"]("L,l,k,d,c") * t1["bb"]("d,j") * r1["bb"]("R,a,l") + 0.50 * t2["bbbb"]("d,c,i,j") * l2["bbbb"]("L,l,k,d,c") * r1["bb"]("R,a,l")) * t1["bb"]("b,k");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
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

        // flops: o2v2L2  = o3v3L2
        //  mems: o2v2L2  = o2v2L2
        tmps_["1_abba_LLovvo"]("L,R,i,a,b,j")  = l2["abab"]("L,i,k,c,a") * r2["abab"]("R,c,b,j,k");

        // D_oovv_abab += -1.00 r2_abab(d,b,i,l) t2_abab(a,c,k,j) l2_abab(k,l,d,c)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,c,k,j") * tmps_["1_abba_LLovvo"]("L,R,k,c,b,i");

        // D_oovv_abab += -1.00 r2_abab(d,b,i,l) t1_aa(a,k) t1_bb(c,j) l2_abab(k,l,d,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["bb"]("c,j") * tmps_["1_abba_LLovvo"]("L,R,k,c,b,i") * t1["aa"]("a,k");
        tmps_["1_abba_LLovvo"].~TArrayD();

        // flops: o1v3L2  = o2v4L1 o1v4L2
        //  mems: o1v3L2  = o0v4L1 o1v3L2
        tmps_["2_aaaa_LLvvvo"]("L,R,a,b,c,i")  = l2["aaaa"]("L,j,k,a,d") * t2["aaaa"]("b,c,j,k") * r1["aa"]("R,d,i");

        // D_vovv_aaaa  = -0.50 r1_aa(d,i) t2_aaaa(c,b,j,k) l2_aaaa(j,k,a,d)
        D_vovv_aaaa("L,R,a,i,b,c")  = -0.50 * tmps_["2_aaaa_LLvvvo"]("L,R,a,c,b,i");

        // D_ovvv_aaaa  = +0.50 r1_aa(d,i) t2_aaaa(c,b,j,k) l2_aaaa(j,k,a,d)
        D_ovvv_aaaa("L,R,i,a,b,c")  = 0.50 * tmps_["2_aaaa_LLvvvo"]("L,R,a,c,b,i");

        // D_oovv_aaaa += -0.50 P(i,j) r1_aa(d,i) t2_aaaa(b,a,k,l) t1_aa(c,j) l2_aaaa(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t1["aa"]("c,j") * tmps_["2_aaaa_LLvvvo"]("L,R,c,b,a,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();
        tmps_["2_aaaa_LLvvvo"].~TArrayD();

        // flops: o1v3L2  = o2v4L1 o1v4L2
        //  mems: o1v3L2  = o0v4L1 o1v3L2
        tmps_["3_bbbb_LLvvvo"]("L,R,a,b,c,i")  = l2["bbbb"]("L,j,k,a,d") * t2["bbbb"]("b,c,j,k") * r1["bb"]("R,d,i");

        // D_vovv_bbbb  = -0.50 r1_bb(d,i) t2_bbbb(c,b,j,k) l2_bbbb(j,k,a,d)
        D_vovv_bbbb("L,R,a,i,b,c")  = -0.50 * tmps_["3_bbbb_LLvvvo"]("L,R,a,c,b,i");

        // D_ovvv_bbbb  = +0.50 r1_bb(d,i) t2_bbbb(c,b,j,k) l2_bbbb(j,k,a,d)
        D_ovvv_bbbb("L,R,i,a,b,c")  = 0.50 * tmps_["3_bbbb_LLvvvo"]("L,R,a,c,b,i");

        // D_oovv_bbbb += -0.50 P(i,j) r1_bb(d,i) t2_bbbb(b,a,k,l) t1_bb(c,j) l2_bbbb(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t1["bb"]("c,j") * tmps_["3_bbbb_LLvvvo"]("L,R,c,b,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
        tmps_["3_bbbb_LLvvvo"].~TArrayD();

        // flops: o0v4L1  = o2v4L1
        //  mems: o0v4L1  = o0v4L1
        tmps_["4_bbbb_Lvvvv"]("L,a,b,c,d")  = l2["bbbb"]("L,i,j,a,b") * t2["bbbb"]("c,d,i,j");

        // D_vvvv_bbbb  = +0.50 r0 t2_bbbb(d,c,i,j) l2_bbbb(i,j,a,b)
        // flops: o0v4L2  = o0v4L2
        //  mems: o0v4L2  = o0v4L2
        D_vvvv_bbbb("L,R,a,b,c,d")  = 0.50 * r0("R") * tmps_["4_bbbb_Lvvvv"]("L,a,b,d,c");

        // flops: o1v3L2  = o2v3L2
        //  mems: o1v3L2  = o1v3L2
        tmps_["5_aaaa_LLovvv"]("L,R,i,a,b,c")  = l2["aaaa"]("L,j,i,a,b") * r1["aa"]("R,c,j");

        // D_vvov_aaaa  = +1.00 r1_aa(c,j) l2_aaaa(j,i,a,b)
        D_vvov_aaaa("L,R,a,b,i,c")  = tmps_["5_aaaa_LLovvv"]("L,R,i,a,b,c");

        // D_vvvo_aaaa  = -1.00 r1_aa(c,j) l2_aaaa(j,i,a,b)
        D_vvvo_aaaa("L,R,a,b,c,i")  = -1.00 * tmps_["5_aaaa_LLovvv"]("L,R,i,a,b,c");

        // D_vvvv_aaaa  = -1.00 P(c,d) r1_aa(c,j) t1_aa(d,i) l2_aaaa(j,i,a,b)
        // flops: o0v4L2  = o1v4L2
        //  mems: o0v4L2  = o0v4L2
        tmps_["perm_aaaa_LLvvvv"]("L,R,a,b,c,d")  = t1["aa"]("d,i") * tmps_["5_aaaa_LLovvv"]("L,R,i,a,b,c");
        D_vvvv_aaaa("L,R,a,b,c,d")  = -1.00 * tmps_["perm_aaaa_LLvvvv"]("L,R,a,b,c,d");
        D_vvvv_aaaa("L,R,a,b,c,d") += tmps_["perm_aaaa_LLvvvv"]("L,R,a,b,d,c");
        tmps_["perm_aaaa_LLvvvv"].~TArrayD();
        tmps_["5_aaaa_LLovvv"].~TArrayD();

        // flops: o1v3L2  = o2v3L2
        //  mems: o1v3L2  = o1v3L2
        tmps_["6_bbbb_LLovvv"]("L,R,i,a,b,c")  = l2["bbbb"]("L,j,i,a,b") * r1["bb"]("R,c,j");

        // D_vvov_bbbb  = +1.00 r1_bb(c,j) l2_bbbb(j,i,a,b)
        D_vvov_bbbb("L,R,a,b,i,c")  = tmps_["6_bbbb_LLovvv"]("L,R,i,a,b,c");

        // D_vvvo_bbbb  = -1.00 r1_bb(c,j) l2_bbbb(j,i,a,b)
        D_vvvo_bbbb("L,R,a,b,c,i")  = -1.00 * tmps_["6_bbbb_LLovvv"]("L,R,i,a,b,c");

        // D_vvvv_bbbb += -1.00 P(c,d) r1_bb(c,j) t1_bb(d,i) l2_bbbb(j,i,a,b)
        // flops: o0v4L2 += o1v4L2
        //  mems: o0v4L2 += o0v4L2
        tmps_["perm_bbbb_LLvvvv"]("L,R,a,b,c,d")  = t1["bb"]("d,i") * tmps_["6_bbbb_LLovvv"]("L,R,i,a,b,c");
        D_vvvv_bbbb("L,R,a,b,c,d") -= tmps_["perm_bbbb_LLvvvv"]("L,R,a,b,c,d");
        D_vvvv_bbbb("L,R,a,b,c,d") += tmps_["perm_bbbb_LLvvvv"]("L,R,a,b,d,c");
        tmps_["perm_bbbb_LLvvvv"].~TArrayD();
        tmps_["6_bbbb_LLovvv"].~TArrayD();

        // flops: o3v1L2  = o3v2L1 o4v1L2 o4v1L2
        //  mems: o3v1L2  = o3v1L1 o4v0L2 o3v1L2
        tmps_["7_aaaa_LLooov"]("L,R,i,j,k,a")  = l2["aaaa"]("L,l,j,c,b") * t1["aa"]("c,i") * r1["aa"]("R,b,k") * t1["aa"]("a,l");

        // D_oovo_aaaa  = +1.00 P(i,j) r1_aa(c,i) t1_aa(a,l) t1_aa(b,j) l2_aaaa(l,k,b,c)
        D_oovo_aaaa("L,R,i,j,a,k")  = tmps_["7_aaaa_LLooov"]("L,R,j,k,i,a");
        D_oovo_aaaa("L,R,i,j,a,k") -= tmps_["7_aaaa_LLooov"]("L,R,i,k,j,a");

        // D_ooov_aaaa  = -1.00 P(i,j) r1_aa(c,i) t1_aa(a,l) t1_aa(b,j) l2_aaaa(l,k,b,c)
        D_ooov_aaaa("L,R,i,j,k,a")  = -1.00 * tmps_["7_aaaa_LLooov"]("L,R,j,k,i,a");
        D_ooov_aaaa("L,R,i,j,k,a") += tmps_["7_aaaa_LLooov"]("L,R,i,k,j,a");

        // D_oovv_aaaa += +1.00 P(i,j) r1_aa(d,i) t1_aa(a,k) t1_aa(b,l) t1_aa(c,j) l2_aaaa(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("b,l") * tmps_["7_aaaa_LLooov"]("L,R,j,l,i,a");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();
        tmps_["7_aaaa_LLooov"].~TArrayD();

        // flops: o3v1L2  = o3v2L1 o4v1L2 o4v1L2
        //  mems: o3v1L2  = o3v1L1 o4v0L2 o3v1L2
        tmps_["8_bbbb_LLooov"]("L,R,i,j,k,a")  = l2["bbbb"]("L,l,j,c,b") * t1["bb"]("c,i") * r1["bb"]("R,b,k") * t1["bb"]("a,l");

        // D_oovo_bbbb  = +1.00 P(i,j) r1_bb(c,i) t1_bb(a,l) t1_bb(b,j) l2_bbbb(l,k,b,c)
        D_oovo_bbbb("L,R,i,j,a,k")  = tmps_["8_bbbb_LLooov"]("L,R,j,k,i,a");
        D_oovo_bbbb("L,R,i,j,a,k") -= tmps_["8_bbbb_LLooov"]("L,R,i,k,j,a");

        // D_ooov_bbbb  = -1.00 P(i,j) r1_bb(c,i) t1_bb(a,l) t1_bb(b,j) l2_bbbb(l,k,b,c)
        D_ooov_bbbb("L,R,i,j,k,a")  = -1.00 * tmps_["8_bbbb_LLooov"]("L,R,j,k,i,a");
        D_ooov_bbbb("L,R,i,j,k,a") += tmps_["8_bbbb_LLooov"]("L,R,i,k,j,a");

        // D_oovv_bbbb += +1.00 P(i,j) r1_bb(d,i) t1_bb(a,k) t1_bb(b,l) t1_bb(c,j) l2_bbbb(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("b,l") * tmps_["8_bbbb_LLooov"]("L,R,j,l,i,a");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
        tmps_["8_bbbb_LLooov"].~TArrayD();

        // flops: o2v2L2  = o2v3L1 o2v3L1 o2v2L2
        //  mems: o2v2L2  = o0v2L1 o2v2L1 o2v2L2
        tmps_["9_aaaa_LLvvoo"]("L,R,a,b,i,j")  = l2["abab"]("L,k,l,c,d") * t2["abab"]("a,d,k,l") * t2["aaaa"]("b,c,i,j") * r0("R");

        // D_oovv_aaaa += -0.50 r0 t2_abab(a,c,k,l) t2_aaaa(b,d,i,j) l2_abab(k,l,d,c)
        //             += -0.50 r0 t2_abab(a,c,l,k) t2_aaaa(b,d,i,j) l2_abab(l,k,d,c)
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["9_aaaa_LLvvoo"]("L,R,a,b,i,j");

        // D_oovv_aaaa += +0.50 r0 t2_aaaa(a,c,i,j) t2_abab(b,d,k,l) l2_abab(k,l,c,d)
        //             += +0.50 r0 t2_aaaa(a,c,i,j) t2_abab(b,d,l,k) l2_abab(l,k,c,d)
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["9_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["9_aaaa_LLvvoo"].~TArrayD();

        // flops: o1v3L1  = o2v3L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["10_baab_Lvovv"]("L,a,i,b,c")  = l2["abab"]("L,i,j,b,c") * t1["bb"]("a,j");

        // D_vvvv_abab += -1.00 r1_aa(c,j) t1_bb(d,i) l2_abab(j,i,a,b)
        // flops: o0v4L2 += o1v4L2
        //  mems: o0v4L2 += o0v4L2
        D_vvvv_abab("L,R,a,b,c,d") -= r1["aa"]("R,c,j") * tmps_["10_baab_Lvovv"]("L,d,j,a,b");

        // D_vvvv_abab += -1.00 r0 t1_aa(c,i) t1_bb(d,j) l2_abab(i,j,a,b)
        // flops: o0v4L2 += o1v4L1 o0v4L2
        //  mems: o0v4L2 += o0v4L1 o0v4L2
        D_vvvv_abab("L,R,a,b,c,d") -= t1["aa"]("c,i") * tmps_["10_baab_Lvovv"]("L,d,i,a,b") * r0("R");
        tmps_["10_baab_Lvovv"].~TArrayD();

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["11_baba_Lvooo"]("L,a,i,j,k")  = l1["aa"]("L,k,b") * t2["abab"]("b,a,i,j");

        // D_oovv_abab += +1.00 r1_aa(a,k) t2_abab(c,b,i,j) l1_aa(k,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += r1["aa"]("R,a,k") * tmps_["11_baba_Lvooo"]("L,b,i,j,k");

        // D_oovv_abab += +1.00 r0 t1_aa(a,k) t2_abab(c,b,i,j) l1_aa(k,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t1["aa"]("a,k") * tmps_["11_baba_Lvooo"]("L,b,i,j,k") * r0("R");
        tmps_["11_baba_Lvooo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["13_aabb_Lvoov"]("L,a,i,j,b")  = l2["abab"]("L,k,j,c,b") * t2["aaaa"]("a,c,k,i");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_abab(a,d,i,l) t2_aaaa(b,c,k,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = r2["abab"]("R,a,d,i,l") * tmps_["13_aabb_Lvoov"]("L,b,j,l,d");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r2_bbbb(b,d,l,j) t2_aaaa(a,c,k,i) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r2["bbbb"]("R,b,d,l,j") * tmps_["13_aabb_Lvoov"]("L,a,i,l,d");

        // D_oovv_abab += -1.00 r0 t2_aaaa(a,c,k,i) t2_bbbb(b,d,l,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["bbbb"]("b,d,l,j") * tmps_["13_aabb_Lvoov"]("L,a,i,l,d") * r0("R");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["12_abba_Lvoov"]("L,a,i,j,b")  = l2["abab"]("L,k,j,b,c") * t2["abab"]("a,c,k,i");

        // flops: o3v1L2  = o3v3L1 o3v2L2 o3v2L2 o3v2L2 o3v1L2 o3v1L2
        //  mems: o3v1L2  = o2v2L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        tmps_["14_abba_LLvooo"]("L,R,a,i,j,k")  = -1.00 * l2["bbbb"]("L,m,j,d,c") * t2["abab"]("a,d,k,m") * r1["bb"]("R,c,i");
        tmps_["14_abba_LLvooo"]("L,R,a,i,j,k") += r1["aa"]("R,b,k") * tmps_["12_abba_Lvoov"]("L,a,i,j,b");
        tmps_["14_abba_LLvooo"]("L,R,a,i,j,k") += r1["bb"]("R,c,i") * tmps_["13_aabb_Lvoov"]("L,a,k,j,c");

        // D_ooov_abab  = -1.00 r1_aa(c,j) t2_abab(a,b,l,i) l2_abab(l,k,c,b)
        //             += -1.00 r1_bb(c,i) t2_aaaa(a,b,l,j) l2_abab(l,k,b,c)
        //             += +1.00 r1_bb(c,i) t2_abab(a,b,j,l) l2_bbbb(l,k,b,c)
        D_ooov_abab("L,R,i,j,k,a")  = -1.00 * tmps_["14_abba_LLvooo"]("L,R,a,i,k,j");

        // D_oovo_abab  = +1.00 r1_aa(c,j) t2_abab(a,b,l,i) l2_abab(l,k,c,b)
        //             += +1.00 r1_bb(c,i) t2_aaaa(a,b,l,j) l2_abab(l,k,b,c)
        //             += -1.00 r1_bb(c,i) t2_abab(a,b,j,l) l2_bbbb(l,k,b,c)
        D_oovo_abab("L,R,i,j,a,k")  = tmps_["14_abba_LLvooo"]("L,R,a,i,k,j");

        // D_oovv_abab += -1.00 r1_aa(d,i) t2_abab(a,c,k,j) t1_bb(b,l) l2_abab(k,l,d,c)
        //             += -1.00 r1_bb(d,j) t2_aaaa(a,c,k,i) t1_bb(b,l) l2_abab(k,l,c,d)
        //             += +1.00 r1_bb(d,j) t2_abab(a,c,i,k) t1_bb(b,l) l2_bbbb(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["bb"]("b,l") * tmps_["14_abba_LLvooo"]("L,R,a,j,l,i");
        tmps_["14_abba_LLvooo"].~TArrayD();

        // flops: o1v3L2  = o2v3L2 o2v3L1 o1v3L2 o1v3L2
        //  mems: o1v3L2  = o1v3L2 o1v3L1 o1v3L2 o1v3L2
        tmps_["15_abab_LLvovv"]("L,R,a,i,b,c")  = l2["abab"]("L,j,i,b,c") * r1["aa"]("R,a,j");
        tmps_["15_abab_LLvovv"]("L,R,a,i,b,c") += l2["abab"]("L,j,i,b,c") * t1["aa"]("a,j") * r0("R");

        // D_vvvo_abab  = -1.00 r0 t1_aa(c,j) l2_abab(j,i,a,b)
        //             += -1.00 r1_aa(c,j) l2_abab(j,i,a,b)
        D_vvvo_abab("L,R,a,b,c,i")  = -1.00 * tmps_["15_abab_LLvovv"]("L,R,c,i,a,b");

        // D_vvov_abab  = +1.00 r0 t1_aa(c,j) l2_abab(j,i,a,b)
        //             += +1.00 r1_aa(c,j) l2_abab(j,i,a,b)
        D_vvov_abab("L,R,a,b,i,c")  = tmps_["15_abab_LLvovv"]("L,R,c,i,a,b");
        tmps_["15_abab_LLvovv"].~TArrayD();

        // flops: o0v4L2  = o2v4L2
        //  mems: o0v4L2  = o0v4L2
        tmps_["16_aaaa_LLvvvv"]("L,R,a,b,c,d")  = l2["aaaa"]("L,i,j,a,b") * r2["aaaa"]("R,c,d,i,j");

        // D_vvvv_aaaa += +0.50 r2_aaaa(d,c,j,i) l2_aaaa(j,i,a,b)
        D_vvvv_aaaa("L,R,a,b,c,d") += 0.50 * tmps_["16_aaaa_LLvvvv"]("L,R,a,b,d,c");

        // flops: o1v3L2  = o1v4L2
        //  mems: o1v3L2  = o1v3L2
        tmps_["17_aaaa_LLvvvo"]("L,R,a,b,c,i")  = t1["aa"]("d,i") * tmps_["16_aaaa_LLvvvv"]("L,R,a,d,b,c");
        tmps_["16_aaaa_LLvvvv"].~TArrayD();

        // D_oovv_aaaa += -0.50 r2_aaaa(b,a,l,k) t1_aa(c,i) t1_aa(d,j) l2_aaaa(l,k,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") -= 0.50 * t1["aa"]("d,j") * tmps_["17_aaaa_LLvvvo"]("L,R,d,b,a,i");

        // D_ovvv_aaaa += +0.50 r2_aaaa(c,b,k,j) t1_aa(d,i) l2_aaaa(k,j,a,d)
        D_ovvv_aaaa("L,R,i,a,b,c") += 0.50 * tmps_["17_aaaa_LLvvvo"]("L,R,a,c,b,i");

        // D_vovv_aaaa += -0.50 r2_aaaa(c,b,k,j) t1_aa(d,i) l2_aaaa(k,j,a,d)
        D_vovv_aaaa("L,R,a,i,b,c") -= 0.50 * tmps_["17_aaaa_LLvvvo"]("L,R,a,c,b,i");
        tmps_["17_aaaa_LLvvvo"].~TArrayD();

        // flops: o0v4L2  = o2v4L2
        //  mems: o0v4L2  = o0v4L2
        tmps_["18_bbbb_LLvvvv"]("L,R,a,b,c,d")  = l2["bbbb"]("L,i,j,a,b") * r2["bbbb"]("R,c,d,i,j");

        // D_vvvv_bbbb += +0.50 r2_bbbb(d,c,j,i) l2_bbbb(j,i,a,b)
        D_vvvv_bbbb("L,R,a,b,c,d") += 0.50 * tmps_["18_bbbb_LLvvvv"]("L,R,a,b,d,c");

        // flops: o1v3L2  = o1v4L2
        //  mems: o1v3L2  = o1v3L2
        tmps_["19_bbbb_LLvvvo"]("L,R,a,b,c,i")  = t1["bb"]("d,i") * tmps_["18_bbbb_LLvvvv"]("L,R,a,d,b,c");
        tmps_["18_bbbb_LLvvvv"].~TArrayD();

        // D_oovv_bbbb += -0.50 r2_bbbb(b,a,l,k) t1_bb(c,i) t1_bb(d,j) l2_bbbb(l,k,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") -= 0.50 * t1["bb"]("d,j") * tmps_["19_bbbb_LLvvvo"]("L,R,d,b,a,i");

        // D_ovvv_bbbb += +0.50 r2_bbbb(c,b,k,j) t1_bb(d,i) l2_bbbb(k,j,a,d)
        D_ovvv_bbbb("L,R,i,a,b,c") += 0.50 * tmps_["19_bbbb_LLvvvo"]("L,R,a,c,b,i");

        // D_vovv_bbbb += -0.50 r2_bbbb(c,b,k,j) t1_bb(d,i) l2_bbbb(k,j,a,d)
        D_vovv_bbbb("L,R,a,i,b,c") -= 0.50 * tmps_["19_bbbb_LLvvvo"]("L,R,a,c,b,i");
        tmps_["19_bbbb_LLvvvo"].~TArrayD();

        // flops: o2v2L2  = o3v3L2 o3v3L2 o2v2L2
        //  mems: o2v2L2  = o2v2L2 o2v2L2 o2v2L2
        tmps_["20_bbbb_LLovvo"]("L,R,i,a,b,j")  = l2["abab"]("L,l,i,d,a") * r2["abab"]("R,d,b,l,j");
        tmps_["20_bbbb_LLovvo"]("L,R,i,a,b,j") += l2["bbbb"]("L,k,i,a,c") * r2["bbbb"]("R,b,c,k,j");

        // D_ovov_bbbb  = +1.00 r2_bbbb(b,c,k,i) l2_bbbb(k,j,a,c)
        //             += +1.00 r2_abab(c,b,k,i) l2_abab(k,j,c,a)
        D_ovov_bbbb("L,R,i,a,j,b")  = tmps_["20_bbbb_LLovvo"]("L,R,j,a,b,i");

        // D_ovvo_bbbb  = -1.00 r2_bbbb(b,c,k,i) l2_bbbb(k,j,a,c)
        //             += -1.00 r2_abab(c,b,k,i) l2_abab(k,j,c,a)
        D_ovvo_bbbb("L,R,i,a,b,j")  = -1.00 * tmps_["20_bbbb_LLovvo"]("L,R,j,a,b,i");

        // D_oovv_abab += -1.00 r2_bbbb(b,d,l,j) t2_abab(a,c,i,k) l2_bbbb(l,k,c,d)
        //             += -1.00 r2_abab(d,b,l,j) t2_abab(a,c,i,k) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,c,i,k") * tmps_["20_bbbb_LLovvo"]("L,R,k,c,b,j");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r2_bbbb(a,d,l,i) t2_bbbb(b,c,k,j) l2_bbbb(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r2_abab(d,a,l,i) t2_bbbb(b,c,k,j) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["bbbb"]("b,c,k,j") * tmps_["20_bbbb_LLovvo"]("L,R,k,c,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o3v1L2  = o3v2L2
        //  mems: o3v1L2  = o3v1L2
        tmps_["21_bbbb_LLoovo"]("L,R,i,j,a,k")  = t1["bb"]("b,i") * tmps_["20_bbbb_LLovvo"]("L,R,j,b,a,k");

        // D_ooov_bbbb += -1.00 P(i,j) r2_bbbb(a,c,l,i) t1_bb(b,j) l2_bbbb(l,k,b,c)
        //             += -1.00 P(i,j) r2_abab(c,a,l,i) t1_bb(b,j) l2_abab(l,k,c,b)
        D_ooov_bbbb("L,R,i,j,k,a") -= tmps_["21_bbbb_LLoovo"]("L,R,j,k,a,i");
        D_ooov_bbbb("L,R,i,j,k,a") += tmps_["21_bbbb_LLoovo"]("L,R,i,k,a,j");

        // D_oovo_bbbb += +1.00 P(i,j) r2_bbbb(a,c,l,i) t1_bb(b,j) l2_bbbb(l,k,b,c)
        //             += +1.00 P(i,j) r2_abab(c,a,l,i) t1_bb(b,j) l2_abab(l,k,c,b)
        D_oovo_bbbb("L,R,i,j,a,k") += tmps_["21_bbbb_LLoovo"]("L,R,j,k,a,i");
        D_oovo_bbbb("L,R,i,j,a,k") -= tmps_["21_bbbb_LLoovo"]("L,R,i,k,a,j");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r2_bbbb(a,d,l,i) t1_bb(b,k) t1_bb(c,j) l2_bbbb(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r2_abab(d,a,l,i) t1_bb(b,k) t1_bb(c,j) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("b,k") * tmps_["21_bbbb_LLoovo"]("L,R,j,k,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
        tmps_["21_bbbb_LLoovo"].~TArrayD();

        // flops: o2v2L2  = o3v3L1 o2v2L2
        //  mems: o2v2L2  = o2v2L1 o2v2L2
        tmps_["22_aaaa_LLvovo"]("L,R,a,i,b,j")  = t2["abab"]("a,c,i,k") * tmps_["13_aabb_Lvoov"]("L,b,j,k,c") * r0("R");

        // D_oovv_aaaa += +1.00 P(i,j) r0 t2_abab(a,c,i,k) t2_aaaa(b,d,l,j) l2_abab(l,k,d,c)
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["22_aaaa_LLvovo"]("L,R,a,i,b,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["22_aaaa_LLvovo"]("L,R,a,j,b,i");

        // D_oovv_aaaa += +1.00 P(i,j) r0 t2_aaaa(a,c,k,i) t2_abab(b,d,j,l) l2_abab(k,l,c,d)
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["22_aaaa_LLvovo"]("L,R,b,j,a,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["22_aaaa_LLvovo"]("L,R,b,i,a,j");
        tmps_["22_aaaa_LLvovo"].~TArrayD();

        // flops: o3v1L2  = o3v3L1 o2v2L1 o3v2L2
        //  mems: o3v1L2  = o2v2L1 o2v2L1 o3v1L2
        tmps_["23_bbbb_LLvooo"]("L,R,a,i,j,k")  = (tmps_["24_bbbb_Lvoov"]("L,a,i,j,b") + -1.00 * t2["bbbb"]("a,d,m,i") * l2["bbbb"]("L,m,j,d,b")) * r1["bb"]("R,b,k");

        // D_ooov_bbbb += +1.00 P(i,j) r1_bb(c,i) t2_abab(b,a,l,j) l2_abab(l,k,b,c)
        //             += -1.00 P(i,j) r1_bb(c,i) t2_bbbb(a,b,l,j) l2_bbbb(l,k,b,c)
        D_ooov_bbbb("L,R,i,j,k,a") += tmps_["23_bbbb_LLvooo"]("L,R,a,j,k,i");
        D_ooov_bbbb("L,R,i,j,k,a") -= tmps_["23_bbbb_LLvooo"]("L,R,a,i,k,j");

        // D_oovo_bbbb += -1.00 P(i,j) r1_bb(c,i) t2_abab(b,a,l,j) l2_abab(l,k,b,c)
        //             += +1.00 P(i,j) r1_bb(c,i) t2_bbbb(a,b,l,j) l2_bbbb(l,k,b,c)
        D_oovo_bbbb("L,R,i,j,a,k") -= tmps_["23_bbbb_LLvooo"]("L,R,a,j,k,i");
        D_oovo_bbbb("L,R,i,j,a,k") += tmps_["23_bbbb_LLvooo"]("L,R,a,i,k,j");

        // D_oovv_bbbb += -1.00 P(i,j) P(a,b) r1_bb(d,i) t2_abab(c,a,k,j) t1_bb(b,l) l2_abab(k,l,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_bb(d,i) t2_bbbb(a,c,k,j) t1_bb(b,l) l2_bbbb(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t1["bb"]("b,l") * tmps_["23_bbbb_LLvooo"]("L,R,a,j,l,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
        tmps_["23_bbbb_LLvooo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["24_bbbb_Lvoov"]("L,a,i,j,b")  = l2["abab"]("L,k,j,c,b") * t2["abab"]("c,a,k,i");

        // D_oovv_abab += -1.00 r2_abab(a,d,i,l) t2_abab(c,b,k,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r2["abab"]("R,a,d,i,l") * tmps_["24_bbbb_Lvoov"]("L,b,j,l,d");

        // D_oovv_abab += -1.00 r0 t2_abab(a,c,i,k) t2_abab(d,b,l,j) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,c,i,k") * tmps_["24_bbbb_Lvoov"]("L,b,j,k,c") * r0("R");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r2_bbbb(a,d,l,i) t2_abab(c,b,k,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = r2["bbbb"]("R,a,d,l,i") * tmps_["24_bbbb_Lvoov"]("L,b,j,l,d");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v2L2  = o3v3L1 o2v2L2
        //  mems: o2v2L2  = o2v2L1 o2v2L2
        tmps_["25_bbbb_LLvovo"]("L,R,a,i,b,j")  = t2["bbbb"]("a,c,k,i") * tmps_["24_bbbb_Lvoov"]("L,b,j,k,c") * r0("R");

        // D_oovv_bbbb += +1.00 P(i,j) r0 t2_bbbb(a,c,k,i) t2_abab(d,b,l,j) l2_abab(l,k,d,c)
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["25_bbbb_LLvovo"]("L,R,a,i,b,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["25_bbbb_LLvovo"]("L,R,a,j,b,i");

        // D_oovv_bbbb += +1.00 P(i,j) r0 t2_abab(c,a,k,i) t2_bbbb(b,d,l,j) l2_abab(k,l,c,d)
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["25_bbbb_LLvovo"]("L,R,b,j,a,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["25_bbbb_LLvovo"]("L,R,b,i,a,j");
        tmps_["25_bbbb_LLvovo"].~TArrayD();

        // flops: o4v0L2  = o4v2L2
        //  mems: o4v0L2  = o4v0L2
        tmps_["26_aaaa_LLoooo"]("L,R,i,j,k,l")  = l2["aaaa"]("L,i,j,a,b") * r2["aaaa"]("R,a,b,k,l");

        // D_oooo_aaaa  = +0.50 r2_aaaa(a,b,i,j) l2_aaaa(l,k,a,b)
        D_oooo_aaaa("L,R,i,j,k,l")  = 0.50 * tmps_["26_aaaa_LLoooo"]("L,R,l,k,i,j");

        // D_oovv_aaaa += +0.25 r2_aaaa(c,d,i,j) t2_aaaa(b,a,k,l) l2_aaaa(k,l,c,d)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") += 0.25 * t2["aaaa"]("b,a,k,l") * tmps_["26_aaaa_LLoooo"]("L,R,k,l,i,j");

        // flops: o3v1L2  = o4v1L2
        //  mems: o3v1L2  = o3v1L2
        tmps_["27_aaaa_LLvooo"]("L,R,a,i,j,k")  = t1["aa"]("a,l") * tmps_["26_aaaa_LLoooo"]("L,R,l,i,j,k");
        tmps_["26_aaaa_LLoooo"].~TArrayD();

        // D_ooov_aaaa += +0.50 r2_aaaa(b,c,i,j) t1_aa(a,l) l2_aaaa(l,k,b,c)
        D_ooov_aaaa("L,R,i,j,k,a") += 0.50 * tmps_["27_aaaa_LLvooo"]("L,R,a,k,i,j");

        // D_oovo_aaaa += -0.50 r2_aaaa(b,c,i,j) t1_aa(a,l) l2_aaaa(l,k,b,c)
        D_oovo_aaaa("L,R,i,j,a,k") -= 0.50 * tmps_["27_aaaa_LLvooo"]("L,R,a,k,i,j");

        // D_oovv_aaaa += -0.50 r2_aaaa(c,d,i,j) t1_aa(a,k) t1_aa(b,l) l2_aaaa(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") -= 0.50 * t1["aa"]("b,l") * tmps_["27_aaaa_LLvooo"]("L,R,a,l,i,j");
        tmps_["27_aaaa_LLvooo"].~TArrayD();

        // flops: o4v0L2  = o4v2L2
        //  mems: o4v0L2  = o4v0L2
        tmps_["28_bbbb_LLoooo"]("L,R,i,j,k,l")  = l2["bbbb"]("L,i,j,a,b") * r2["bbbb"]("R,a,b,k,l");

        // D_oooo_bbbb  = +0.50 r2_bbbb(a,b,i,j) l2_bbbb(l,k,a,b)
        D_oooo_bbbb("L,R,i,j,k,l")  = 0.50 * tmps_["28_bbbb_LLoooo"]("L,R,l,k,i,j");

        // D_oovv_bbbb += +0.25 r2_bbbb(c,d,i,j) t2_bbbb(b,a,k,l) l2_bbbb(k,l,c,d)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") += 0.25 * t2["bbbb"]("b,a,k,l") * tmps_["28_bbbb_LLoooo"]("L,R,k,l,i,j");

        // flops: o3v1L2  = o4v1L2
        //  mems: o3v1L2  = o3v1L2
        tmps_["29_bbbb_LLvooo"]("L,R,a,i,j,k")  = t1["bb"]("a,l") * tmps_["28_bbbb_LLoooo"]("L,R,l,i,j,k");
        tmps_["28_bbbb_LLoooo"].~TArrayD();

        // D_ooov_bbbb += +0.50 r2_bbbb(b,c,i,j) t1_bb(a,l) l2_bbbb(l,k,b,c)
        D_ooov_bbbb("L,R,i,j,k,a") += 0.50 * tmps_["29_bbbb_LLvooo"]("L,R,a,k,i,j");

        // D_oovo_bbbb += -0.50 r2_bbbb(b,c,i,j) t1_bb(a,l) l2_bbbb(l,k,b,c)
        D_oovo_bbbb("L,R,i,j,a,k") -= 0.50 * tmps_["29_bbbb_LLvooo"]("L,R,a,k,i,j");

        // D_oovv_bbbb += -0.50 r2_bbbb(c,d,i,j) t1_bb(a,k) t1_bb(b,l) l2_bbbb(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_bbbb("L,R,i,j,a,b") -= 0.50 * t1["bb"]("b,l") * tmps_["29_bbbb_LLvooo"]("L,R,a,l,i,j");
        tmps_["29_bbbb_LLvooo"].~TArrayD();

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["30_bb_Lvv"]("L,a,b")  = l2["abab"]("L,i,j,c,a") * t2["abab"]("c,b,i,j");

        // D_oovv_abab += +0.50 r2_abab(a,d,i,j) t2_abab(c,b,k,l) l2_abab(k,l,c,d)
        //             += +0.50 r2_abab(a,d,i,j) t2_abab(c,b,l,k) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += r2["abab"]("R,a,d,i,j") * tmps_["30_bb_Lvv"]("L,d,b");

        // D_oovv_abab += +0.50 r0 t2_abab(a,c,i,j) t2_abab(d,b,k,l) l2_abab(k,l,d,c)
        //             += +0.50 r0 t2_abab(a,c,i,j) t2_abab(d,b,l,k) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,c,i,j") * tmps_["30_bb_Lvv"]("L,c,b") * r0("R");

        // D_oovv_bbbb += +0.50 P(a,b) r2_bbbb(a,d,i,j) t2_abab(c,b,k,l) l2_abab(k,l,c,d)
        //             += +0.50 P(a,b) r2_bbbb(a,d,i,j) t2_abab(c,b,l,k) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = r2["bbbb"]("R,a,d,i,j") * tmps_["30_bb_Lvv"]("L,d,b");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v2L2  = o2v3L1 o2v2L2
        //  mems: o2v2L2  = o2v2L1 o2v2L2
        tmps_["31_bbbb_LLvoov"]("L,R,a,i,j,b")  = t2["bbbb"]("a,c,i,j") * tmps_["30_bb_Lvv"]("L,c,b") * r0("R");

        // D_oovv_bbbb += +0.50 r0 t2_bbbb(a,c,i,j) t2_abab(d,b,k,l) l2_abab(k,l,d,c)
        //             += +0.50 r0 t2_bbbb(a,c,i,j) t2_abab(d,b,l,k) l2_abab(l,k,d,c)
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["31_bbbb_LLvoov"]("L,R,a,i,j,b");

        // D_oovv_bbbb += -0.50 r0 t2_abab(c,a,k,l) t2_bbbb(b,d,i,j) l2_abab(k,l,c,d)
        //             += -0.50 r0 t2_abab(c,a,l,k) t2_bbbb(b,d,i,j) l2_abab(l,k,c,d)
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["31_bbbb_LLvoov"]("L,R,b,i,j,a");
        tmps_["31_bbbb_LLvoov"].~TArrayD();

        // flops: o0v4L2  = o2v4L2
        //  mems: o0v4L2  = o0v4L2
        tmps_["32_abab_LLvvvv"]("L,R,a,b,c,d")  = l2["abab"]("L,i,j,a,b") * r2["abab"]("R,c,d,i,j");

        // D_vvvv_abab += -0.50 r2_abab(c,d,j,i) l2_abab(j,i,a,b)
        //             += -0.50 r2_abab(c,d,i,j) l2_abab(i,j,a,b)
        D_vvvv_abab("L,R,a,b,c,d") -= tmps_["32_abab_LLvvvv"]("L,R,a,b,c,d");

        // flops: o1v3L2  = o3v3L2 o3v3L2 o2v2L2 o2v3L2 o1v4L2 o1v3L2
        //  mems: o1v3L2  = o2v2L2 o2v2L2 o2v2L2 o1v3L2 o1v3L2 o1v3L2
        tmps_["33_aabb_LLvvvo"]("L,R,a,b,c,i")  = (l2["aaaa"]("L,k,j,b,e") * r2["abab"]("R,e,c,k,i") + l2["abab"]("L,j,m,b,d") * r2["bbbb"]("R,c,d,m,i")) * t1["aa"]("a,j");
        tmps_["33_aabb_LLvvvo"]("L,R,a,b,c,i") += t1["bb"]("d,i") * tmps_["32_abab_LLvvvv"]("L,R,b,d,a,c");
        tmps_["32_abab_LLvvvv"].~TArrayD();

        // D_vovv_abab  = -1.00 r2_abab(d,b,k,i) t1_aa(c,j) l2_aaaa(k,j,a,d)
        //             += -1.00 r2_bbbb(b,d,k,i) t1_aa(c,j) l2_abab(j,k,a,d)
        //             += -0.50 r2_abab(c,b,k,j) t1_bb(d,i) l2_abab(k,j,a,d)
        //             += -0.50 r2_abab(c,b,j,k) t1_bb(d,i) l2_abab(j,k,a,d)
        D_vovv_abab("L,R,a,i,b,c")  = -1.00 * tmps_["33_aabb_LLvvvo"]("L,R,c,a,b,i");

        // D_ovvv_abab  = +1.00 r2_abab(d,b,k,i) t1_aa(c,j) l2_aaaa(k,j,a,d)
        //             += +1.00 r2_bbbb(b,d,k,i) t1_aa(c,j) l2_abab(j,k,a,d)
        //             += +0.50 r2_abab(c,b,k,j) t1_bb(d,i) l2_abab(k,j,a,d)
        //             += +0.50 r2_abab(c,b,j,k) t1_bb(d,i) l2_abab(j,k,a,d)
        D_ovvv_abab("L,R,i,a,b,c")  = tmps_["33_aabb_LLvvvo"]("L,R,c,a,b,i");

        // D_oovv_abab += -1.00 r2_abab(d,b,l,j) t1_aa(a,k) t1_aa(c,i) l2_aaaa(l,k,c,d)
        //             += -1.00 r2_bbbb(b,d,l,j) t1_aa(a,k) t1_aa(c,i) l2_abab(k,l,c,d)
        //             += -0.50 r2_abab(a,b,l,k) t1_aa(c,i) t1_bb(d,j) l2_abab(l,k,c,d)
        //             += -0.50 r2_abab(a,b,k,l) t1_aa(c,i) t1_bb(d,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["aa"]("c,i") * tmps_["33_aabb_LLvvvo"]("L,R,a,c,b,j");
        tmps_["33_aabb_LLvvvo"].~TArrayD();

        // flops: o2v2L2  = o3v3L2
        //  mems: o2v2L2  = o2v2L2
        tmps_["34_aaaa_LLovvo"]("L,R,i,a,b,j")  = l2["aaaa"]("L,k,i,a,c") * r2["aaaa"]("R,b,c,k,j");

        // D_ovvo_aaaa  = -1.00 r2_aaaa(b,c,k,i) l2_aaaa(k,j,a,c)
        D_ovvo_aaaa("L,R,i,a,b,j")  = -1.00 * tmps_["34_aaaa_LLovvo"]("L,R,j,a,b,i");

        // D_ovov_aaaa  = +1.00 r2_aaaa(b,c,k,i) l2_aaaa(k,j,a,c)
        D_ovov_aaaa("L,R,i,a,j,b")  = tmps_["34_aaaa_LLovvo"]("L,R,j,a,b,i");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_aaaa(a,d,l,i) t2_aaaa(b,c,k,j) l2_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["aaaa"]("b,c,k,j") * tmps_["34_aaaa_LLovvo"]("L,R,k,c,a,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r2_aaaa(a,d,l,i) t2_abab(c,b,k,j) l2_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("c,b,k,j") * tmps_["34_aaaa_LLovvo"]("L,R,k,c,a,i");

        // flops: o2v2L2  = o3v3L2
        //  mems: o2v2L2  = o2v2L2
        tmps_["35_aaaa_LLovvo"]("L,R,i,a,b,j")  = l2["abab"]("L,i,k,a,c") * r2["abab"]("R,b,c,j,k");

        // D_ovov_aaaa += +1.00 r2_abab(b,c,i,k) l2_abab(j,k,a,c)
        D_ovov_aaaa("L,R,i,a,j,b") += tmps_["35_aaaa_LLovvo"]("L,R,j,a,b,i");

        // D_ovvo_aaaa += -1.00 r2_abab(b,c,i,k) l2_abab(j,k,a,c)
        D_ovvo_aaaa("L,R,i,a,b,j") -= tmps_["35_aaaa_LLovvo"]("L,R,j,a,b,i");

        // flops: o3v1L2  = o2v2L2 o3v2L2
        //  mems: o3v1L2  = o2v2L2 o3v1L2
        tmps_["36_aaaa_LLovoo"]("L,R,i,a,j,k")  = (tmps_["34_aaaa_LLovvo"]("L,R,i,b,a,j") + tmps_["35_aaaa_LLovvo"]("L,R,i,b,a,j")) * t1["aa"]("b,k");

        // D_ooov_aaaa += -1.00 P(i,j) r2_aaaa(a,c,l,i) t1_aa(b,j) l2_aaaa(l,k,b,c)
        //             += -1.00 P(i,j) r2_abab(a,c,i,l) t1_aa(b,j) l2_abab(k,l,b,c)
        D_ooov_aaaa("L,R,i,j,k,a") -= tmps_["36_aaaa_LLovoo"]("L,R,k,a,i,j");
        D_ooov_aaaa("L,R,i,j,k,a") += tmps_["36_aaaa_LLovoo"]("L,R,k,a,j,i");

        // D_oovo_aaaa += +1.00 P(i,j) r2_aaaa(a,c,l,i) t1_aa(b,j) l2_aaaa(l,k,b,c)
        //             += +1.00 P(i,j) r2_abab(a,c,i,l) t1_aa(b,j) l2_abab(k,l,b,c)
        D_oovo_aaaa("L,R,i,j,a,k") += tmps_["36_aaaa_LLovoo"]("L,R,k,a,i,j");
        D_oovo_aaaa("L,R,i,j,a,k") -= tmps_["36_aaaa_LLovoo"]("L,R,k,a,j,i");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_aaaa(a,d,l,i) t1_aa(b,k) t1_aa(c,j) l2_aaaa(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r2_abab(a,d,i,l) t1_aa(b,k) t1_aa(c,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("b,k") * tmps_["36_aaaa_LLovoo"]("L,R,k,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();
        tmps_["36_aaaa_LLovoo"].~TArrayD();

        // flops: o2v2L2  = o3v3L2
        //  mems: o2v2L2  = o2v2L2
        tmps_["37_bbaa_LLovvo"]("L,R,i,a,b,j")  = l2["abab"]("L,k,i,c,a") * r2["aaaa"]("R,b,c,k,j");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_aaaa(a,d,l,i) t2_abab(b,c,j,k) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["abab"]("b,c,j,k") * tmps_["37_bbaa_LLovvo"]("L,R,k,c,a,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -1.00 r2_aaaa(a,d,l,i) t2_bbbb(b,c,k,j) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["bbbb"]("b,c,k,j") * tmps_["37_bbaa_LLovvo"]("L,R,k,c,a,i");

        // flops: o1v1L2  = o2v2L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["40_bb_LLov"]("L,R,i,a")  = l2["bbbb"]("L,j,i,a,b") * r1["bb"]("R,b,j");

        // flops: o1v1L2  = o2v2L2
        //  mems: o1v1L2  = o1v1L2
        tmps_["39_bb_LLov"]("L,R,i,a")  = l2["abab"]("L,j,i,b,a") * r1["aa"]("R,b,j");

        // flops: o2v2L2  = o3v3L2
        //  mems: o2v2L2  = o2v2L2
        tmps_["38_baab_LLovvo"]("L,R,i,a,b,j")  = l2["abab"]("L,k,i,a,c") * r2["abab"]("R,b,c,k,j");

        // flops: o3v1L2  = o3v2L2 o3v2L2 o3v2L2 o3v1L2 o3v2L1 o4v2L2 o3v1L2 o3v2L2 o3v1L2 o3v2L2 o3v1L2 o3v1L2
        //  mems: o3v1L2  = o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L1 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2 o3v1L2
        tmps_["41_bbaa_LLoovo"]("L,R,i,j,a,k")  = t1["aa"]("b,k") * tmps_["38_baab_LLovvo"]("L,R,j,b,a,i");
        tmps_["41_bbaa_LLoovo"]("L,R,i,j,a,k") -= l1["bb"]("L,j,c") * r2["abab"]("R,a,c,k,i");
        tmps_["41_bbaa_LLoovo"]("L,R,i,j,a,k") += tmps_["40_bb_LLov"]("L,R,j,c") * t2["abab"]("a,c,k,i");
        tmps_["41_bbaa_LLoovo"]("L,R,i,j,a,k") += t1["bb"]("c,i") * l2["bbbb"]("L,m,j,c,d") * r2["abab"]("R,a,d,k,m");
        tmps_["41_bbaa_LLoovo"]("L,R,i,j,a,k") -= tmps_["39_bb_LLov"]("L,R,j,c") * t2["abab"]("a,c,k,i");
        tmps_["41_bbaa_LLoovo"]("L,R,i,j,a,k") += t1["bb"]("c,i") * tmps_["37_bbaa_LLovvo"]("L,R,j,c,a,k");
        tmps_["37_bbaa_LLovvo"].~TArrayD();

        // D_ooov_abab += -1.00 r2_aaaa(a,c,l,j) t1_bb(b,i) l2_abab(l,k,c,b)
        //             += +1.00 r2_abab(a,b,j,i) l1_bb(k,b)
        //             += -1.00 r1_bb(c,l) t2_abab(a,b,j,i) l2_bbbb(l,k,b,c)
        //             += -1.00 r2_abab(a,c,j,l) t1_bb(b,i) l2_bbbb(l,k,b,c)
        //             += +1.00 r1_aa(c,l) t2_abab(a,b,j,i) l2_abab(l,k,c,b)
        //             += -1.00 r2_abab(a,c,l,i) t1_aa(b,j) l2_abab(l,k,b,c)
        D_ooov_abab("L,R,i,j,k,a") -= tmps_["41_bbaa_LLoovo"]("L,R,i,k,a,j");

        // D_oovo_abab += +1.00 r2_aaaa(a,c,l,j) t1_bb(b,i) l2_abab(l,k,c,b)
        //             += -1.00 r2_abab(a,b,j,i) l1_bb(k,b)
        //             += +1.00 r1_bb(c,l) t2_abab(a,b,j,i) l2_bbbb(l,k,b,c)
        //             += +1.00 r2_abab(a,c,j,l) t1_bb(b,i) l2_bbbb(l,k,b,c)
        //             += -1.00 r1_aa(c,l) t2_abab(a,b,j,i) l2_abab(l,k,c,b)
        //             += +1.00 r2_abab(a,c,l,i) t1_aa(b,j) l2_abab(l,k,b,c)
        D_oovo_abab("L,R,i,j,a,k") += tmps_["41_bbaa_LLoovo"]("L,R,i,k,a,j");

        // D_oovv_abab += -1.00 r2_aaaa(a,d,l,i) t1_bb(b,k) t1_bb(c,j) l2_abab(l,k,d,c)
        //             += +1.00 r2_abab(a,c,i,j) t1_bb(b,k) l1_bb(k,c)
        //             += -1.00 r1_bb(d,l) t2_abab(a,c,i,j) t1_bb(b,k) l2_bbbb(l,k,c,d)
        //             += -1.00 r2_abab(a,d,i,l) t1_bb(b,k) t1_bb(c,j) l2_bbbb(l,k,c,d)
        //             += +1.00 r1_aa(d,l) t2_abab(a,c,i,j) t1_bb(b,k) l2_abab(l,k,d,c)
        //             += -1.00 r2_abab(a,d,l,j) t1_bb(b,k) t1_aa(c,i) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["bb"]("b,k") * tmps_["41_bbaa_LLoovo"]("L,R,j,k,a,i");
        tmps_["41_bbaa_LLoovo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["42_baab_Lvoov"]("L,a,i,j,b")  = l2["abab"]("L,j,k,c,b") * t2["abab"]("c,a,i,k");

        // D_oovv_abab += -1.00 r2_abab(a,d,l,j) t2_abab(c,b,i,k) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r2["abab"]("R,a,d,l,j") * tmps_["42_baab_Lvoov"]("L,b,i,l,d");

        // D_oovv_abab += -1.00 r0 t2_abab(a,c,k,j) t2_abab(d,b,i,l) l2_abab(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,c,k,j") * tmps_["42_baab_Lvoov"]("L,b,i,k,c") * r0("R");

        // D_oovv_abab += -1.00 r1_bb(d,j) t1_aa(a,l) t2_abab(c,b,i,k) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r1["bb"]("R,d,j") * tmps_["42_baab_Lvoov"]("L,b,i,l,d") * t1["aa"]("a,l");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["43_bbaa_Lvoov"]("L,a,i,j,b")  = l2["abab"]("L,j,k,b,c") * t2["bbbb"]("a,c,k,i");

        // D_oovv_abab += -1.00 r1_aa(d,i) t1_aa(a,l) t2_bbbb(b,c,k,j) l2_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r1["aa"]("R,d,i") * tmps_["43_bbaa_Lvoov"]("L,b,j,l,d") * t1["aa"]("a,l");

        // flops: o3v1L1  = o3v2L1 o3v2L1 o3v1L1
        //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1
        tmps_["44_bbaa_Lvooo"]("L,a,i,j,k")  = t1["aa"]("b,k") * tmps_["43_bbaa_Lvoov"]("L,a,i,j,b");
        tmps_["44_bbaa_Lvooo"]("L,a,i,j,k") += t1["bb"]("c,i") * tmps_["42_baab_Lvoov"]("L,a,k,j,c");
        tmps_["42_baab_Lvoov"].~TArrayD();

        // D_oovv_abab += -1.00 r1_aa(a,l) t2_bbbb(b,c,k,j) t1_aa(d,i) l2_abab(l,k,d,c)
        //             += -1.00 r1_aa(a,l) t2_abab(c,b,i,k) t1_bb(d,j) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r1["aa"]("R,a,l") * tmps_["44_bbaa_Lvooo"]("L,b,j,l,i");

        // D_oovv_abab += -1.00 r0 t1_aa(a,l) t2_bbbb(b,c,k,j) t1_aa(d,i) l2_abab(l,k,d,c)
        //             += -1.00 r0 t1_aa(a,l) t2_abab(c,b,i,k) t1_bb(d,j) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["aa"]("a,l") * tmps_["44_bbaa_Lvooo"]("L,b,j,l,i") * r0("R");
        tmps_["44_bbaa_Lvooo"].~TArrayD();

        // flops: o4v0L1  = o4v2L1
        //  mems: o4v0L1  = o4v0L1
        tmps_["45_abab_Loooo"]("L,i,j,k,l")  = l2["abab"]("L,k,l,a,b") * t2["abab"]("a,b,i,j");

        // D_oooo_abab  = -0.50 r0 t2_abab(b,a,i,j) l2_abab(k,l,b,a)
        //             += -0.50 r0 t2_abab(a,b,i,j) l2_abab(k,l,a,b)
        // flops: o4v0L2  = o4v0L2
        //  mems: o4v0L2  = o4v0L2
        D_oooo_abab("L,R,i,j,k,l")  = -1.00 * r0("R") * tmps_["45_abab_Loooo"]("L,i,j,k,l");

        // D_oovv_abab += -0.25 r2_abab(a,b,l,k) t2_abab(d,c,i,j) l2_abab(l,k,d,c)
        //             += -0.25 r2_abab(a,b,l,k) t2_abab(c,d,i,j) l2_abab(l,k,c,d)
        //             += -0.25 r2_abab(a,b,k,l) t2_abab(d,c,i,j) l2_abab(k,l,d,c)
        //             += -0.25 r2_abab(a,b,k,l) t2_abab(c,d,i,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o4v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r2["abab"]("R,a,b,l,k") * tmps_["45_abab_Loooo"]("L,i,j,l,k");

        // D_oovv_abab += -0.25 r0 t2_abab(a,b,k,l) t2_abab(d,c,i,j) l2_abab(k,l,d,c)
        //             += -0.25 r0 t2_abab(a,b,k,l) t2_abab(c,d,i,j) l2_abab(k,l,c,d)
        //             += -0.25 r0 t2_abab(a,b,l,k) t2_abab(d,c,i,j) l2_abab(l,k,d,c)
        //             += -0.25 r0 t2_abab(a,b,l,k) t2_abab(c,d,i,j) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o4v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2["abab"]("a,b,k,l") * tmps_["45_abab_Loooo"]("L,i,j,k,l") * r0("R");

        // D_oovv_abab += -0.50 r1_bb(b,l) t1_aa(a,k) t2_abab(d,c,i,j) l2_abab(k,l,d,c)
        //             += -0.50 r1_bb(b,l) t1_aa(a,k) t2_abab(c,d,i,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o4v1L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r1["bb"]("R,b,l") * tmps_["45_abab_Loooo"]("L,i,j,k,l") * t1["aa"]("a,k");

        // flops: o4v0L1  = o3v2L1 o4v1L1
        //  mems: o4v0L1  = o3v1L1 o4v0L1
        tmps_["46_aabb_Loooo"]("L,i,j,k,l")  = l2["abab"]("L,j,k,b,a") * t1["aa"]("b,i") * t1["bb"]("a,l");

        // D_oooo_abab += -1.00 r0 t1_aa(a,i) t1_bb(b,j) l2_abab(k,l,a,b)
        // flops: o4v0L2 += o4v0L2
        //  mems: o4v0L2 += o4v0L2
        D_oooo_abab("L,R,i,j,k,l") -= r0("R") * tmps_["46_aabb_Loooo"]("L,i,k,l,j");

        // D_oovv_abab += -1.00 r1_bb(b,l) t1_aa(a,k) t1_aa(c,i) t1_bb(d,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o4v1L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r1["bb"]("R,b,l") * tmps_["46_aabb_Loooo"]("L,i,k,l,j") * t1["aa"]("a,k");

        // flops: o3v1L1  = o4v0L1 o4v1L1
        //  mems: o3v1L1  = o4v0L1 o3v1L1
        tmps_["47_aabb_Looov"]("L,i,j,k,a")  = t1["bb"]("a,l") * (tmps_["46_aabb_Loooo"]("L,i,j,l,k") + tmps_["45_abab_Loooo"]("L,i,k,j,l"));

        // D_oovv_abab += -1.00 r1_aa(a,l) t1_bb(b,k) t1_aa(c,i) t1_bb(d,j) l2_abab(l,k,c,d)
        //             += -0.50 r1_aa(a,l) t1_bb(b,k) t2_abab(d,c,i,j) l2_abab(l,k,d,c)
        //             += -0.50 r1_aa(a,l) t1_bb(b,k) t2_abab(c,d,i,j) l2_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r1["aa"]("R,a,l") * tmps_["47_aabb_Looov"]("L,i,l,j,b");

        // D_oovv_abab += -1.00 r0 t1_aa(a,k) t1_bb(b,l) t1_aa(c,i) t1_bb(d,j) l2_abab(k,l,c,d)
        //             += -0.50 r0 t1_aa(a,k) t1_bb(b,l) t2_abab(d,c,i,j) l2_abab(k,l,d,c)
        //             += -0.50 r0 t1_aa(a,k) t1_bb(b,l) t2_abab(c,d,i,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["aa"]("a,k") * tmps_["47_aabb_Looov"]("L,i,k,j,b") * r0("R");
        tmps_["47_aabb_Looov"].~TArrayD();

        // flops: o2v0L1  = o3v2L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["48_bb_Loo"]("L,i,j")  = l2["abab"]("L,k,j,a,b") * t2["abab"]("a,b,k,i");

        // D_oovv_abab += +0.50 r2_abab(a,b,i,l) t2_abab(d,c,k,j) l2_abab(k,l,d,c)
        //             += +0.50 r2_abab(a,b,i,l) t2_abab(c,d,k,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += r2["abab"]("R,a,b,i,l") * tmps_["48_bb_Loo"]("L,j,l");

        // D_oovv_abab += +0.50 r0 t2_abab(a,b,i,l) t2_abab(d,c,k,j) l2_abab(k,l,d,c)
        //             += +0.50 r0 t2_abab(a,b,i,l) t2_abab(c,d,k,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += t2["abab"]("a,b,i,l") * tmps_["48_bb_Loo"]("L,j,l") * r0("R");

        // D_oovv_bbbb += +0.50 P(i,j) r2_bbbb(b,a,l,i) t2_abab(d,c,k,j) l2_abab(k,l,d,c)
        //             += +0.50 P(i,j) r2_bbbb(b,a,l,i) t2_abab(c,d,k,j) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = r2["bbbb"]("R,b,a,l,i") * tmps_["48_bb_Loo"]("L,j,l");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // D_oovv_bbbb += -0.50 P(i,j) r0 t2_bbbb(b,a,l,j) t2_abab(d,c,k,i) l2_abab(k,l,d,c)
        //             += -0.50 P(i,j) r0 t2_bbbb(b,a,l,j) t2_abab(c,d,k,i) l2_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = t2["bbbb"]("b,a,l,j") * tmps_["48_bb_Loo"]("L,i,l") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v0L1  = o3v2L1
        //  mems: o2v0L1  = o2v0L1
        tmps_["49_bb_Loo"]("L,i,j")  = l2["bbbb"]("L,k,j,a,b") * t2["bbbb"]("a,b,k,i");

        // D_oovv_abab += +0.50 r0 t2_abab(a,b,i,l) t2_bbbb(d,c,k,j) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * t2["abab"]("a,b,i,l") * tmps_["49_bb_Loo"]("L,j,l") * r0("R");

        // D_oovv_bbbb += -0.50 P(i,j) r0 t2_bbbb(b,a,l,j) t2_bbbb(d,c,k,i) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * t2["bbbb"]("b,a,l,j") * tmps_["49_bb_Loo"]("L,i,l") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();
    }
} // hilbert
#endif