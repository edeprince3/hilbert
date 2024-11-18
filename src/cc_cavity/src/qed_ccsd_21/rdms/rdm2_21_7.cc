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


    void EOM_EE_QED_RDM_21::rdm2_21_7() {

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


        // D_ovvv_bbbb += -1.00 P(b,c) r1_bb(d,k) t1_bb(c,i) t1_1_bb(b,j) l2_1_bbbb(k,j,a,d)
        //             += +1.00 P(b,c) r1_aa(d,k) t1_bb(c,i) t1_1_bb(b,j) l2_1_abab(k,j,d,a)
        //             += -1.00 P(b,c) r1_1_bb(d,k) t1_bb(b,j) t1_bb(c,i) l2_1_bbbb(k,j,a,d)
        //             += +1.00 P(b,c) r1_1_aa(d,k) t1_bb(b,j) t1_bb(c,i) l2_1_abab(k,j,d,a)
        //             += -1.00 P(b,c) r1_bb(d,k) t1_bb(b,j) t1_bb(c,i) l2_bbbb(k,j,a,d)
        //             += +1.00 P(b,c) r1_aa(d,k) t1_bb(b,j) t1_bb(c,i) l2_abab(k,j,d,a)
        //             += -1.00 P(b,c) r1_aa(d,k) t1_bb(c,j) t1_1_bb(b,i) l2_1_abab(k,j,d,a)
        //             += -0.50 P(b,c) r1_bb(b,i) t2_1_bbbb(c,d,j,k) l2_1_bbbb(j,k,a,d)
        //             += -0.50 P(b,c) r1_1_bb(b,i) t2_bbbb(c,d,j,k) l2_1_bbbb(j,k,a,d)
        //             += -1.00 P(b,c) r0 t1_bb(c,j) t1_1_bb(b,i) l1_1_bb(j,a)
        //             += -1.00 P(b,c) r0_1 t2_abab(d,b,j,i) t1_bb(c,k) l2_1_abab(j,k,d,a)
        //             += +0.50 P(b,c) r0_1 t2_abab(d,b,j,k) t1_bb(c,i) l2_1_abab(j,k,d,a)
        //             += +0.50 P(b,c) r0_1 t2_abab(d,b,k,j) t1_bb(c,i) l2_1_abab(k,j,d,a)
        //             += +1.00 P(b,c) r0 t1_bb(c,i) t1_1_bb(b,j) l1_1_bb(j,a)
        //             += +1.00 P(b,c) r0 t1_bb(b,j) t1_bb(c,i) l1_bb(j,a)
        //             += -1.00 P(b,c) r0 t2_abab(d,b,j,i) t1_bb(c,k) l2_abab(j,k,d,a)
        //             += -1.00 P(b,c) r0 t1_bb(c,k) t2_1_abab(d,b,j,i) l2_1_abab(j,k,d,a)
        //             += -1.00 P(b,c) r1_1_bb(b,k) t1_bb(c,j) t1_bb(d,i) l2_1_bbbb(k,j,a,d)
        //             += -1.00 P(b,c) r1_bb(b,k) t1_bb(c,j) t1_1_bb(d,i) l2_1_bbbb(k,j,a,d)
        //             += -1.00 P(b,c) r1_1_bb(b,k) t2_bbbb(c,d,j,i) l2_1_bbbb(k,j,a,d)
        //             += -1.00 P(b,c) r1_bb(b,k) t1_bb(d,i) t1_1_bb(c,j) l2_1_bbbb(k,j,a,d)
        //             += +1.00 P(b,c) r1_bb(b,k) t2_abab(d,c,j,i) l2_abab(j,k,d,a)
        //             += -1.00 P(b,c) r1_bb(b,k) t2_1_bbbb(c,d,j,i) l2_1_bbbb(k,j,a,d)
        //             += -1.00 P(b,c) r1_bb(b,k) t1_bb(c,j) t1_bb(d,i) l2_bbbb(k,j,a,d)
        //             += -1.00 P(b,c) r0 t1_bb(c,k) t1_bb(d,i) t1_1_bb(b,j) l2_1_bbbb(j,k,a,d)
        //             += +1.00 P(b,c) r1_1_bb(b,k) t2_abab(d,c,j,i) l2_1_abab(j,k,d,a)
        //             += -1.00 P(b,c) r1_bb(d,i) t1_bb(c,k) t1_1_bb(b,j) l2_1_bbbb(j,k,a,d)
        //             += -1.00 P(b,c) r1_bb(b,k) t2_bbbb(c,d,j,i) l2_bbbb(k,j,a,d)
        //             += +1.00 P(b,c) r1_bb(b,j) t1_1_bb(c,i) l1_1_bb(j,a)
        //             += +1.00 P(b,c) r1_bb(b,k) t2_1_abab(d,c,j,i) l2_1_abab(j,k,d,a)
        //             += -1.00 P(b,c) r1_bb(b,i) t1_1_bb(c,j) l1_1_bb(j,a)
        //             += -1.00 P(b,c) r1_bb(b,i) t1_bb(c,j) l1_bb(j,a)
        //             += -1.00 P(b,c) r0_1 t2_bbbb(b,d,j,i) t1_bb(c,k) l2_1_bbbb(j,k,a,d)
        //             += +0.50 P(b,c) r0 t1_bb(c,i) t2_1_abab(d,b,j,k) l2_1_abab(j,k,d,a)
        //             += +0.50 P(b,c) r0 t1_bb(c,i) t2_1_abab(d,b,k,j) l2_1_abab(k,j,d,a)
        //             += +0.50 P(b,c) r2_abab(d,b,k,j) t1_1_bb(c,i) l2_1_abab(k,j,d,a)
        //             += +0.50 P(b,c) r2_abab(d,b,j,k) t1_1_bb(c,i) l2_1_abab(j,k,d,a)
        //             += +0.50 P(b,c) r2_bbbb(b,d,k,j) t1_1_bb(c,i) l2_1_bbbb(k,j,a,d)
        //             += -1.00 P(b,c) r0 t1_bb(c,k) t2_1_bbbb(b,d,j,i) l2_1_bbbb(j,k,a,d)
        //             += +1.00 P(b,c) r0_1 t1_bb(b,j) t1_bb(c,i) l1_1_bb(j,a)
        //             += +0.50 P(b,c) r0_1 t2_bbbb(b,d,j,k) t1_bb(c,i) l2_1_bbbb(j,k,a,d)
        //             += +0.50 P(b,c) r0 t1_bb(c,i) t2_1_bbbb(b,d,j,k) l2_1_bbbb(j,k,a,d)
        //             += -1.00 P(b,c) r0 t2_bbbb(c,d,k,i) t1_1_bb(b,j) l2_1_bbbb(j,k,a,d)
        //             += +0.50 P(b,c) r2_bbbb(b,d,k,j) t1_bb(c,i) l2_bbbb(k,j,a,d)
        //             += +0.50 P(b,c) r2_1_abab(d,b,k,j) t1_bb(c,i) l2_1_abab(k,j,d,a)
        //             += +0.50 P(b,c) r2_1_abab(d,b,j,k) t1_bb(c,i) l2_1_abab(j,k,d,a)
        //             += +0.50 P(b,c) r2_1_bbbb(b,d,k,j) t1_bb(c,i) l2_1_bbbb(k,j,a,d)
        //             += +0.50 P(b,c) r2_abab(d,b,k,j) t1_bb(c,i) l2_abab(k,j,d,a)
        //             += +0.50 P(b,c) r2_abab(d,b,j,k) t1_bb(c,i) l2_abab(j,k,d,a)
        //             += -0.50 P(b,c) r1_1_bb(b,i) t2_abab(d,c,j,k) l2_1_abab(j,k,d,a)
        //             += -0.50 P(b,c) r1_1_bb(b,i) t2_abab(d,c,k,j) l2_1_abab(k,j,d,a)
        //             += -0.50 P(b,c) r1_bb(b,i) t2_bbbb(c,d,j,k) l2_bbbb(j,k,a,d)
        //             += +0.50 P(b,c) r0 t2_bbbb(b,d,j,k) t1_bb(c,i) l2_bbbb(j,k,a,d)
        //             += -0.50 P(b,c) r0 t2_bbbb(c,d,j,k) t1_1_bb(b,i) l2_1_bbbb(j,k,a,d)
        //             += +1.00 P(b,c) r1_1_bb(b,j) t1_bb(c,i) l1_1_bb(j,a)
        //             += +1.00 P(b,c) r1_bb(b,j) t1_bb(c,i) l1_bb(j,a)
        //             += +0.50 P(b,c) r0 t2_abab(d,b,j,k) t1_bb(c,i) l2_abab(j,k,d,a)
        //             += +0.50 P(b,c) r0 t2_abab(d,b,k,j) t1_bb(c,i) l2_abab(k,j,d,a)
        //             += -1.00 P(b,c) r2_abab(d,b,k,i) t1_1_bb(c,j) l2_1_abab(k,j,d,a)
        //             += -1.00 P(b,c) r2_bbbb(b,d,k,i) t1_1_bb(c,j) l2_1_bbbb(k,j,a,d)
        //             += -0.50 P(b,c) r0 t2_abab(d,c,j,k) t1_1_bb(b,i) l2_1_abab(j,k,d,a)
        //             += -0.50 P(b,c) r0 t2_abab(d,c,k,j) t1_1_bb(b,i) l2_1_abab(k,j,d,a)
        //             += -1.00 P(b,c) r0 t2_bbbb(b,d,j,i) t1_bb(c,k) l2_bbbb(j,k,a,d)
        //             += +1.00 P(b,c) r0 t2_abab(d,c,k,i) t1_1_bb(b,j) l2_1_abab(k,j,d,a)
        //             += -1.00 P(b,c) r2_1_bbbb(b,d,k,i) t1_bb(c,j) l2_1_bbbb(k,j,a,d)
        //             += -1.00 P(b,c) r2_abab(d,b,k,i) t1_bb(c,j) l2_abab(k,j,d,a)
        //             += -1.00 P(b,c) r2_1_abab(d,b,k,i) t1_bb(c,j) l2_1_abab(k,j,d,a)
        //             += -1.00 P(b,c) r2_bbbb(b,d,k,i) t1_bb(c,j) l2_bbbb(k,j,a,d)
        //             += +1.00 P(b,c) r1_bb(d,k) t1_bb(c,j) t1_1_bb(b,i) l2_1_bbbb(k,j,a,d)
        //             += -1.00 P(b,c) r1_1_bb(b,i) t1_bb(c,j) l1_1_bb(j,a)
        //             += -0.50 P(b,c) r1_bb(b,i) t2_1_abab(d,c,j,k) l2_1_abab(j,k,d,a)
        //             += -0.50 P(b,c) r1_bb(b,i) t2_1_abab(d,c,k,j) l2_1_abab(k,j,d,a)
        //             += -0.50 P(b,c) r1_bb(b,i) t2_abab(d,c,j,k) l2_abab(j,k,d,a)
        //             += -0.50 P(b,c) r1_bb(b,i) t2_abab(d,c,k,j) l2_abab(k,j,d,a)
        D_ovvv_bbbb("L,R,i,a,b,c") -= tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i");
        D_ovvv_bbbb("L,R,i,a,b,c") += tmps_["279_bbbb_LLvvvo"]("L,R,a,c,b,i");

        // D_vovv_bbbb += +1.00 P(b,c) r1_bb(d,k) t1_bb(c,i) t1_1_bb(b,j) l2_1_bbbb(k,j,a,d)
        //             += -1.00 P(b,c) r1_aa(d,k) t1_bb(c,i) t1_1_bb(b,j) l2_1_abab(k,j,d,a)
        //             += +1.00 P(b,c) r1_1_bb(d,k) t1_bb(b,j) t1_bb(c,i) l2_1_bbbb(k,j,a,d)
        //             += -1.00 P(b,c) r1_1_aa(d,k) t1_bb(b,j) t1_bb(c,i) l2_1_abab(k,j,d,a)
        //             += +1.00 P(b,c) r1_bb(d,k) t1_bb(b,j) t1_bb(c,i) l2_bbbb(k,j,a,d)
        //             += -1.00 P(b,c) r1_aa(d,k) t1_bb(b,j) t1_bb(c,i) l2_abab(k,j,d,a)
        //             += +1.00 P(b,c) r1_aa(d,k) t1_bb(c,j) t1_1_bb(b,i) l2_1_abab(k,j,d,a)
        //             += +0.50 P(b,c) r1_bb(b,i) t2_1_bbbb(c,d,j,k) l2_1_bbbb(j,k,a,d)
        //             += +0.50 P(b,c) r1_1_bb(b,i) t2_bbbb(c,d,j,k) l2_1_bbbb(j,k,a,d)
        //             += +1.00 P(b,c) r0 t1_bb(c,j) t1_1_bb(b,i) l1_1_bb(j,a)
        //             += +1.00 P(b,c) r0_1 t2_abab(d,b,j,i) t1_bb(c,k) l2_1_abab(j,k,d,a)
        //             += -0.50 P(b,c) r0_1 t2_abab(d,b,j,k) t1_bb(c,i) l2_1_abab(j,k,d,a)
        //             += -0.50 P(b,c) r0_1 t2_abab(d,b,k,j) t1_bb(c,i) l2_1_abab(k,j,d,a)
        //             += -1.00 P(b,c) r0 t1_bb(c,i) t1_1_bb(b,j) l1_1_bb(j,a)
        //             += -1.00 P(b,c) r0 t1_bb(b,j) t1_bb(c,i) l1_bb(j,a)
        //             += +1.00 P(b,c) r0 t2_abab(d,b,j,i) t1_bb(c,k) l2_abab(j,k,d,a)
        //             += +1.00 P(b,c) r0 t1_bb(c,k) t2_1_abab(d,b,j,i) l2_1_abab(j,k,d,a)
        //             += +1.00 P(b,c) r1_1_bb(b,k) t1_bb(c,j) t1_bb(d,i) l2_1_bbbb(k,j,a,d)
        //             += +1.00 P(b,c) r1_bb(b,k) t1_bb(c,j) t1_1_bb(d,i) l2_1_bbbb(k,j,a,d)
        //             += +1.00 P(b,c) r1_1_bb(b,k) t2_bbbb(c,d,j,i) l2_1_bbbb(k,j,a,d)
        //             += +1.00 P(b,c) r1_bb(b,k) t1_bb(d,i) t1_1_bb(c,j) l2_1_bbbb(k,j,a,d)
        //             += -1.00 P(b,c) r1_bb(b,k) t2_abab(d,c,j,i) l2_abab(j,k,d,a)
        //             += +1.00 P(b,c) r1_bb(b,k) t2_1_bbbb(c,d,j,i) l2_1_bbbb(k,j,a,d)
        //             += +1.00 P(b,c) r1_bb(b,k) t1_bb(c,j) t1_bb(d,i) l2_bbbb(k,j,a,d)
        //             += +1.00 P(b,c) r0 t1_bb(c,k) t1_bb(d,i) t1_1_bb(b,j) l2_1_bbbb(j,k,a,d)
        //             += -1.00 P(b,c) r1_1_bb(b,k) t2_abab(d,c,j,i) l2_1_abab(j,k,d,a)
        //             += +1.00 P(b,c) r1_bb(d,i) t1_bb(c,k) t1_1_bb(b,j) l2_1_bbbb(j,k,a,d)
        //             += +1.00 P(b,c) r1_bb(b,k) t2_bbbb(c,d,j,i) l2_bbbb(k,j,a,d)
        //             += -1.00 P(b,c) r1_bb(b,j) t1_1_bb(c,i) l1_1_bb(j,a)
        //             += -1.00 P(b,c) r1_bb(b,k) t2_1_abab(d,c,j,i) l2_1_abab(j,k,d,a)
        //             += +1.00 P(b,c) r1_bb(b,i) t1_1_bb(c,j) l1_1_bb(j,a)
        //             += +1.00 P(b,c) r1_bb(b,i) t1_bb(c,j) l1_bb(j,a)
        //             += +1.00 P(b,c) r0_1 t2_bbbb(b,d,j,i) t1_bb(c,k) l2_1_bbbb(j,k,a,d)
        //             += -0.50 P(b,c) r0 t1_bb(c,i) t2_1_abab(d,b,j,k) l2_1_abab(j,k,d,a)
        //             += -0.50 P(b,c) r0 t1_bb(c,i) t2_1_abab(d,b,k,j) l2_1_abab(k,j,d,a)
        //             += -0.50 P(b,c) r2_abab(d,b,k,j) t1_1_bb(c,i) l2_1_abab(k,j,d,a)
        //             += -0.50 P(b,c) r2_abab(d,b,j,k) t1_1_bb(c,i) l2_1_abab(j,k,d,a)
        //             += -0.50 P(b,c) r2_bbbb(b,d,k,j) t1_1_bb(c,i) l2_1_bbbb(k,j,a,d)
        //             += +1.00 P(b,c) r0 t1_bb(c,k) t2_1_bbbb(b,d,j,i) l2_1_bbbb(j,k,a,d)
        //             += -1.00 P(b,c) r0_1 t1_bb(b,j) t1_bb(c,i) l1_1_bb(j,a)
        //             += -0.50 P(b,c) r0_1 t2_bbbb(b,d,j,k) t1_bb(c,i) l2_1_bbbb(j,k,a,d)
        //             += -0.50 P(b,c) r0 t1_bb(c,i) t2_1_bbbb(b,d,j,k) l2_1_bbbb(j,k,a,d)
        //             += +1.00 P(b,c) r0 t2_bbbb(c,d,k,i) t1_1_bb(b,j) l2_1_bbbb(j,k,a,d)
        //             += -0.50 P(b,c) r2_bbbb(b,d,k,j) t1_bb(c,i) l2_bbbb(k,j,a,d)
        //             += -0.50 P(b,c) r2_1_abab(d,b,k,j) t1_bb(c,i) l2_1_abab(k,j,d,a)
        //             += -0.50 P(b,c) r2_1_abab(d,b,j,k) t1_bb(c,i) l2_1_abab(j,k,d,a)
        //             += -0.50 P(b,c) r2_1_bbbb(b,d,k,j) t1_bb(c,i) l2_1_bbbb(k,j,a,d)
        //             += -0.50 P(b,c) r2_abab(d,b,k,j) t1_bb(c,i) l2_abab(k,j,d,a)
        //             += -0.50 P(b,c) r2_abab(d,b,j,k) t1_bb(c,i) l2_abab(j,k,d,a)
        //             += +0.50 P(b,c) r1_1_bb(b,i) t2_abab(d,c,j,k) l2_1_abab(j,k,d,a)
        //             += +0.50 P(b,c) r1_1_bb(b,i) t2_abab(d,c,k,j) l2_1_abab(k,j,d,a)
        //             += +0.50 P(b,c) r1_bb(b,i) t2_bbbb(c,d,j,k) l2_bbbb(j,k,a,d)
        //             += -0.50 P(b,c) r0 t2_bbbb(b,d,j,k) t1_bb(c,i) l2_bbbb(j,k,a,d)
        //             += +0.50 P(b,c) r0 t2_bbbb(c,d,j,k) t1_1_bb(b,i) l2_1_bbbb(j,k,a,d)
        //             += -1.00 P(b,c) r1_1_bb(b,j) t1_bb(c,i) l1_1_bb(j,a)
        //             += -1.00 P(b,c) r1_bb(b,j) t1_bb(c,i) l1_bb(j,a)
        //             += -0.50 P(b,c) r0 t2_abab(d,b,j,k) t1_bb(c,i) l2_abab(j,k,d,a)
        //             += -0.50 P(b,c) r0 t2_abab(d,b,k,j) t1_bb(c,i) l2_abab(k,j,d,a)
        //             += +1.00 P(b,c) r2_abab(d,b,k,i) t1_1_bb(c,j) l2_1_abab(k,j,d,a)
        //             += +1.00 P(b,c) r2_bbbb(b,d,k,i) t1_1_bb(c,j) l2_1_bbbb(k,j,a,d)
        //             += +0.50 P(b,c) r0 t2_abab(d,c,j,k) t1_1_bb(b,i) l2_1_abab(j,k,d,a)
        //             += +0.50 P(b,c) r0 t2_abab(d,c,k,j) t1_1_bb(b,i) l2_1_abab(k,j,d,a)
        //             += +1.00 P(b,c) r0 t2_bbbb(b,d,j,i) t1_bb(c,k) l2_bbbb(j,k,a,d)
        //             += -1.00 P(b,c) r0 t2_abab(d,c,k,i) t1_1_bb(b,j) l2_1_abab(k,j,d,a)
        //             += +1.00 P(b,c) r2_1_bbbb(b,d,k,i) t1_bb(c,j) l2_1_bbbb(k,j,a,d)
        //             += +1.00 P(b,c) r2_abab(d,b,k,i) t1_bb(c,j) l2_abab(k,j,d,a)
        //             += +1.00 P(b,c) r2_1_abab(d,b,k,i) t1_bb(c,j) l2_1_abab(k,j,d,a)
        //             += +1.00 P(b,c) r2_bbbb(b,d,k,i) t1_bb(c,j) l2_bbbb(k,j,a,d)
        //             += -1.00 P(b,c) r1_bb(d,k) t1_bb(c,j) t1_1_bb(b,i) l2_1_bbbb(k,j,a,d)
        //             += +1.00 P(b,c) r1_1_bb(b,i) t1_bb(c,j) l1_1_bb(j,a)
        //             += +0.50 P(b,c) r1_bb(b,i) t2_1_abab(d,c,j,k) l2_1_abab(j,k,d,a)
        //             += +0.50 P(b,c) r1_bb(b,i) t2_1_abab(d,c,k,j) l2_1_abab(k,j,d,a)
        //             += +0.50 P(b,c) r1_bb(b,i) t2_abab(d,c,j,k) l2_abab(j,k,d,a)
        //             += +0.50 P(b,c) r1_bb(b,i) t2_abab(d,c,k,j) l2_abab(k,j,d,a)
        D_vovv_bbbb("L,R,a,i,b,c") += tmps_["279_bbbb_LLvvvo"]("L,R,a,b,c,i");
        D_vovv_bbbb("L,R,a,i,b,c") -= tmps_["279_bbbb_LLvvvo"]("L,R,a,c,b,i");
        tmps_["279_bbbb_LLvvvo"].~TArrayD();

        // flops: o3v1L1  = o4v0L1 o4v1L1
        //  mems: o3v1L1  = o4v0L1 o3v1L1
        tmps_["280_abab_Looov"]("L,i,j,k,a")  = t1["bb"]("a,l") * (tmps_["272_abab_Loooo"]("L,i,j,k,l") + tmps_["273_baab_Loooo"]("L,j,i,k,l"));

        // D_oovv_abab += -0.50 r1_1_aa(a,l) t1_bb(b,k) t2_abab(d,c,i,j) l2_1_abab(l,k,d,c)
        //             += -0.50 r1_1_aa(a,l) t1_bb(b,k) t2_abab(c,d,i,j) l2_1_abab(l,k,c,d)
        //             += -1.00 r1_1_aa(a,l) t1_bb(b,k) t1_aa(c,i) t1_bb(d,j) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= r1_1["aa"]("R,a,l") * tmps_["280_abab_Looov"]("L,i,j,l,b");

        // D_oovv_abab += -0.50 r0 t1_bb(b,l) t2_abab(d,c,i,j) t1_1_aa(a,k) l2_1_abab(k,l,d,c)
        //             += -0.50 r0 t1_bb(b,l) t2_abab(c,d,i,j) t1_1_aa(a,k) l2_1_abab(k,l,c,d)
        //             += -1.00 r0 t1_bb(b,l) t1_aa(c,i) t1_bb(d,j) t1_1_aa(a,k) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1_1["aa"]("a,k") * tmps_["280_abab_Looov"]("L,i,j,k,b") * r0("R");

        // D_oovv_abab += -0.50 r0_1 t1_aa(a,k) t1_bb(b,l) t2_abab(d,c,i,j) l2_1_abab(k,l,d,c)
        //             += -0.50 r0_1 t1_aa(a,k) t1_bb(b,l) t2_abab(c,d,i,j) l2_1_abab(k,l,c,d)
        //             += -1.00 r0_1 t1_aa(a,k) t1_bb(b,l) t1_aa(c,i) t1_bb(d,j) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t1["aa"]("a,k") * tmps_["280_abab_Looov"]("L,i,j,k,b") * r0_1("R");
        tmps_["280_abab_Looov"].~TArrayD();

        // flops: o2v2L2  = o3v3L2
        //  mems: o2v2L2  = o2v2L2
        tmps_["281_aabb_LLovvo"]("L,R,i,a,b,j")  = l2_1["aaaa"]("L,k,i,a,c") * r2["abab"]("R,c,b,k,j");

        // D_oovv_abab += -1.00 r2_abab(d,b,l,j) t2_1_aaaa(a,c,k,i) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= t2_1["aaaa"]("a,c,k,i") * tmps_["281_aabb_LLovvo"]("L,R,k,c,b,j");

        // D_oovv_abab += -1.00 r2_abab(d,b,l,j) t1_aa(a,k) t1_1_aa(c,i) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["281_aabb_LLovvo"]("L,R,k,c,b,j") * t1_1["aa"]("c,i") * t1["aa"]("a,k");

        // D_oovv_abab += -1.00 r2_abab(d,b,l,j) t1_aa(c,i) t1_1_aa(a,k) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["281_aabb_LLovvo"]("L,R,k,c,b,j") * t1["aa"]("c,i") * t1_1["aa"]("a,k");

        // D_oovv_bbbb += +1.00 P(i,j) P(a,b) r2_abab(d,a,l,i) t2_1_abab(c,b,k,j) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["281_aabb_LLovvo"]("L,R,k,c,a,i") * t2_1["abab"]("c,b,k,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o2v2L2  = o3v3L2
        //  mems: o2v2L2  = o2v2L2
        tmps_["282_aabb_LLovvo"]("L,R,i,a,b,j")  = l2_1["abab"]("L,i,k,a,c") * r2["bbbb"]("R,b,c,k,j");

        // D_oovv_abab += -1.00 r2_bbbb(b,d,l,j) t1_aa(a,k) t1_1_aa(c,i) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["282_aabb_LLovvo"]("L,R,k,c,b,j") * t1_1["aa"]("c,i") * t1["aa"]("a,k");

        // D_oovv_abab += -1.00 r2_bbbb(b,d,l,j) t1_aa(c,i) t1_1_aa(a,k) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o3v2L2 o3v2L2
        //  mems: o2v2L2 += o3v1L2 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["282_aabb_LLovvo"]("L,R,k,c,b,j") * t1["aa"]("c,i") * t1_1["aa"]("a,k");

        // flops: o0v4L1  = o2v4L1
        //  mems: o0v4L1  = o0v4L1
        tmps_["283_abab_Lvvvv"]("L,a,b,c,d")  = l2_1["abab"]("L,i,j,a,b") * t2["abab"]("c,d,i,j");

        // D_oovv_abab += -0.50 r1_1_bb(d,j) t2_abab(a,b,k,l) t1_aa(c,i) l2_1_abab(k,l,c,d)
        //             += -0.50 r1_1_bb(d,j) t2_abab(a,b,l,k) t1_aa(c,i) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o1v4L1 o2v3L2
        //  mems: o2v2L2 += o1v3L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["283_abab_Lvvvv"]("L,c,d,a,b") * t1["aa"]("c,i") * r1_1["bb"]("R,d,j");

        // D_oovv_abab += -0.50 r1_bb(d,j) t2_abab(a,b,k,l) t1_1_aa(c,i) l2_1_abab(k,l,c,d)
        //             += -0.50 r1_bb(d,j) t2_abab(a,b,l,k) t1_1_aa(c,i) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o1v4L1 o2v3L2
        //  mems: o2v2L2 += o1v3L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["283_abab_Lvvvv"]("L,c,d,a,b") * t1_1["aa"]("c,i") * r1["bb"]("R,d,j");

        // D_vvvv_abab += -0.50 r0_1 t2_abab(c,d,i,j) l2_1_abab(i,j,a,b)
        //             += -0.50 r0_1 t2_abab(c,d,j,i) l2_1_abab(j,i,a,b)
        // flops: o0v4L2 += o0v4L2
        //  mems: o0v4L2 += o0v4L2
        D_vvvv_abab("L,R,a,b,c,d") -= tmps_["283_abab_Lvvvv"]("L,a,b,c,d") * r0_1("R");

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["284_bbaa_Lvoov"]("L,a,i,j,b")  = t2["abab"]("c,a,k,i") * l2["aaaa"]("L,k,j,b,c");

        // D_oovv_abab += -1.00 r0 t1_aa(a,l) t2_abab(c,b,k,j) t1_aa(d,i) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o3v2L1 o2v2L2
        //  mems: o2v2L2 += o3v1L1 o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["284_bbaa_Lvoov"]("L,b,j,l,d") * t1["aa"]("d,i") * t1["aa"]("a,l") * r0("R");

        // D_oovv_bbbb += +1.00 P(i,j) r0 t2_abab(c,a,k,i) t2_abab(d,b,l,j) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j")  = tmps_["284_bbaa_Lvoov"]("L,a,i,l,d") * t2["abab"]("d,b,l,j") * r0("R");
        D_oovv_bbbb("L,R,i,j,a,b") += tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_bbbb("L,R,i,j,a,b") -= tmps_["perm_bbbb_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_bbbb_LLvvoo"].~TArrayD();

        // flops: o1v3L1  = o3v3L1 o2v3L1
        //  mems: o1v3L1  = o2v2L1 o1v3L1
        tmps_["285_abba_Lvvov"]("L,a,b,i,c")  = l2_1["aaaa"]("L,k,j,a,d") * t2_1["abab"]("d,b,k,i") * t1["aa"]("c,j");

        // D_oovv_abab += -1.00 r0 t1_aa(a,l) t1_aa(d,i) t2_1_abab(c,b,k,j) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["285_abba_Lvvov"]("L,d,b,j,a") * t1["aa"]("d,i") * r0("R");

        // flops: o1v3L1  = o3v3L1 o2v3L1
        //  mems: o1v3L1  = o2v2L1 o1v3L1
        tmps_["286_abba_Lvvov"]("L,a,b,i,c")  = l2_1["aaaa"]("L,k,j,a,d") * t2["abab"]("d,b,k,i") * t1["aa"]("c,j");

        // D_oovv_abab += -1.00 r0_1 t1_aa(a,l) t2_abab(c,b,k,j) t1_aa(d,i) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["286_abba_Lvvov"]("L,d,b,j,a") * t1["aa"]("d,i") * r0_1("R");

        // flops: o0v4L1  = o2v4L1 o2v4L1 o0v4L1
        //  mems: o0v4L1  = o0v4L1 o0v4L1 o0v4L1
        tmps_["287_abab_Lvvvv"]("L,a,b,c,d")  = l2["abab"]("L,i,j,a,b") * t2["abab"]("c,d,i,j");
        tmps_["287_abab_Lvvvv"]("L,a,b,c,d") += l2_1["abab"]("L,i,j,a,b") * t2_1["abab"]("c,d,i,j");

        // D_oovv_abab += -0.50 r1_bb(d,j) t2_abab(a,b,k,l) t1_aa(c,i) l2_abab(k,l,c,d)
        //             += -0.50 r1_bb(d,j) t2_abab(a,b,l,k) t1_aa(c,i) l2_abab(l,k,c,d)
        //             += -0.50 r1_bb(d,j) t1_aa(c,i) t2_1_abab(a,b,k,l) l2_1_abab(k,l,c,d)
        //             += -0.50 r1_bb(d,j) t1_aa(c,i) t2_1_abab(a,b,l,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o1v4L1 o2v3L2
        //  mems: o2v2L2 += o1v3L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["287_abab_Lvvvv"]("L,c,d,a,b") * t1["aa"]("c,i") * r1["bb"]("R,d,j");

        // D_vvvv_abab += -0.50 r0 t2_abab(c,d,i,j) l2_abab(i,j,a,b)
        //             += -0.50 r0 t2_abab(c,d,j,i) l2_abab(j,i,a,b)
        //             += -0.50 r0 t2_1_abab(c,d,i,j) l2_1_abab(i,j,a,b)
        //             += -0.50 r0 t2_1_abab(c,d,j,i) l2_1_abab(j,i,a,b)
        // flops: o0v4L2 += o0v4L2
        //  mems: o0v4L2 += o0v4L2
        D_vvvv_abab("L,R,a,b,c,d") -= tmps_["287_abab_Lvvvv"]("L,a,b,c,d") * r0("R");

        // flops: o1v3L1  = o1v4L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["288_aabb_Lvvvo"]("L,a,b,c,i")  = tmps_["283_abab_Lvvvv"]("L,a,d,b,c") * t1["bb"]("d,i");

        // D_oovv_abab += -0.50 r1_1_aa(d,i) t2_abab(a,b,k,l) t1_bb(c,j) l2_1_abab(k,l,d,c)
        //             += -0.50 r1_1_aa(d,i) t2_abab(a,b,l,k) t1_bb(c,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["288_aabb_Lvvvo"]("L,d,a,b,j") * r1_1["aa"]("R,d,i");

        // D_oovv_abab += -0.50 r0_1 t2_abab(a,b,k,l) t1_aa(c,i) t1_bb(d,j) l2_1_abab(k,l,c,d)
        //             += -0.50 r0_1 t2_abab(a,b,l,k) t1_aa(c,i) t1_bb(d,j) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["288_aabb_Lvvvo"]("L,c,a,b,j") * t1["aa"]("c,i") * r0_1("R");

        // D_oovv_abab += -0.50 r0 t2_abab(a,b,k,l) t1_bb(d,j) t1_1_aa(c,i) l2_1_abab(k,l,c,d)
        //             += -0.50 r0 t2_abab(a,b,l,k) t1_bb(d,j) t1_1_aa(c,i) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["288_aabb_Lvvvo"]("L,c,a,b,j") * t1_1["aa"]("c,i") * r0("R");

        // flops: o1v3L1  = o1v4L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["289_aabb_Lvvvo"]("L,a,b,c,i")  = tmps_["287_abab_Lvvvv"]("L,a,d,b,c") * t1["bb"]("d,i");

        // D_oovv_abab += -0.50 r1_aa(d,i) t2_abab(a,b,k,l) t1_bb(c,j) l2_abab(k,l,d,c)
        //             += -0.50 r1_aa(d,i) t2_abab(a,b,l,k) t1_bb(c,j) l2_abab(l,k,d,c)
        //             += -0.50 r1_aa(d,i) t1_bb(c,j) t2_1_abab(a,b,k,l) l2_1_abab(k,l,d,c)
        //             += -0.50 r1_aa(d,i) t1_bb(c,j) t2_1_abab(a,b,l,k) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["289_aabb_Lvvvo"]("L,d,a,b,j") * r1["aa"]("R,d,i");

        // D_oovv_abab += -0.50 r0 t2_abab(a,b,k,l) t1_aa(c,i) t1_bb(d,j) l2_abab(k,l,c,d)
        //             += -0.50 r0 t2_abab(a,b,l,k) t1_aa(c,i) t1_bb(d,j) l2_abab(l,k,c,d)
        //             += -0.50 r0 t1_aa(c,i) t1_bb(d,j) t2_1_abab(a,b,k,l) l2_1_abab(k,l,c,d)
        //             += -0.50 r0 t1_aa(c,i) t1_bb(d,j) t2_1_abab(a,b,l,k) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["289_aabb_Lvvvo"]("L,c,a,b,j") * t1["aa"]("c,i") * r0("R");

        // flops: o1v3L1  = o1v4L1
        //  mems: o1v3L1  = o1v3L1
        tmps_["290_aabb_Lvvvo"]("L,a,b,c,i")  = tmps_["283_abab_Lvvvv"]("L,a,d,b,c") * t1_1["bb"]("d,i");

        // D_oovv_abab += -0.50 r1_aa(d,i) t2_abab(a,b,k,l) t1_1_bb(c,j) l2_1_abab(k,l,d,c)
        //             += -0.50 r1_aa(d,i) t2_abab(a,b,l,k) t1_1_bb(c,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["290_aabb_Lvvvo"]("L,d,a,b,j") * r1["aa"]("R,d,i");

        // D_oovv_abab += -0.50 r0 t2_abab(a,b,k,l) t1_aa(d,i) t1_1_bb(c,j) l2_1_abab(k,l,d,c)
        //             += -0.50 r0 t2_abab(a,b,l,k) t1_aa(d,i) t1_1_bb(c,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= tmps_["290_aabb_Lvvvo"]("L,d,a,b,j") * t1["aa"]("d,i") * r0("R");

        // flops: o1v3L2  = o0v2L1 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L2 o3v2L2 o3v2L2 o2v3L2 o1v3L2 o3v2L1 o3v2L2 o2v3L2 o1v3L2 o1v4L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L1 o1v3L2 o1v3L2 o3v2L2 o3v2L2 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L2 o1v3L2 o3v2L1 o3v2L2 o2v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L1 o1v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L2 o1v3L2 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L1 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o3v2L1 o3v2L2 o2v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L1 o1v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o3v2L1 o3v2L2 o2v3L2 o1v3L2 o3v2L2 o3v2L2 o2v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o3v3L1 o2v3L2 o1v3L2 o3v2L2 o3v2L2 o2v3L2 o1v3L2 o3v2L1 o3v2L1 o2v3L1 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o1v3L2 o1v3L2 o3v3L2 o3v3L2 o2v2L2 o3v3L2 o2v2L2 o3v3L2 o2v2L2 o2v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v4L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        //  mems: o1v3L2  = o0v2L1 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o3v1L2 o2v2L2 o1v3L2 o1v3L2 o3v1L1 o2v2L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L1 o1v3L2 o1v3L2 o3v1L2 o2v2L2 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L2 o1v3L2 o3v1L1 o2v2L2 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L1 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o3v1L1 o2v2L2 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L1 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o3v1L1 o2v2L2 o1v3L2 o1v3L2 o3v1L2 o2v2L2 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o2v2L1 o1v3L2 o1v3L2 o3v1L2 o2v2L2 o1v3L2 o1v3L2 o3v1L1 o2v2L1 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o2v2L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L1 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2 o1v3L2
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i")  = 0.50 * (tmps_["246_aa_Lvv"]("L,a,b") + 2.00 * tmps_["252_aa_Lvv"]("L,a,b")) * tmps_["225_bb_Lvo"]("R,c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += tmps_["250_aa_Lvv"]("L,b,a") * r1["bb"]("R,c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["13_abba_Lvoov"]("L,b,i,l,a") * t1["bb"]("c,l") * r0_1("R");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += tmps_["250_aa_Lvv"]("L,b,a") * tmps_["224_bb_Lvo"]("R,c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += tmps_["249_aa_Lvv"]("L,a,b") * tmps_["224_bb_Lvo"]("R,c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += tmps_["257_aa_LLvv"]("L,R,a,b") * t1["bb"]("c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += l1["aa"]("L,j,a") * r2["abab"]("R,b,c,j,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,j,l,a,e") * r1["bb"]("R,e,i") * t1["bb"]("c,l") * t1_1["aa"]("b,j");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,j,l,a,e") * t1["bb"]("e,i") * r1["bb"]("R,c,l") * t1_1["aa"]("b,j");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["283_abab_Lvvvv"]("L,a,e,b,c") * r1_1["bb"]("R,e,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += l1_1["aa"]("L,j,a") * t2["abab"]("b,c,j,i") * r0_1("R");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,k,m,a,e") * t1["bb"]("e,i") * t1_1["bb"]("c,m") * r1["aa"]("R,b,k");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2["abab"]("L,k,m,a,e") * t2["bbbb"]("c,e,m,i") * r1["aa"]("R,b,k");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2["abab"]("L,k,m,a,e") * t1["bb"]("e,i") * t1["bb"]("c,m") * t1["aa"]("b,k") * r0("R");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,k,m,a,e") * r1_1["bb"]("R,e,i") * t1["bb"]("c,m") * t1["aa"]("b,k");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,k,m,a,e") * t1["bb"]("e,i") * t1["bb"]("c,m") * r1_1["aa"]("R,b,k");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,j,l,a,e") * t1["bb"]("e,i") * r1_1["bb"]("R,c,l") * t1["aa"]("b,j");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,k,m,a,e") * t2["bbbb"]("c,e,m,i") * r1_1["aa"]("R,b,k");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,k,m,a,e") * t1["bb"]("e,i") * t1["bb"]("c,m") * t1["aa"]("b,k") * r0_1("R");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2["abab"]("L,j,l,a,e") * t2["abab"]("b,e,j,i") * r1["bb"]("R,c,l");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += l2["aaaa"]("L,k,j,a,d") * t2["abab"]("d,c,j,i") * r1["aa"]("R,b,k");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2["abab"]("L,k,m,a,e") * t1["bb"]("e,i") * t1["bb"]("c,m") * r1["aa"]("R,b,k");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += l1["aa"]("L,j,a") * t2["abab"]("b,c,j,i") * r0("R");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["127_aa_LLov"]("L,R,j,a") * t2["abab"]("b,c,j,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += l2_1["aaaa"]("L,k,j,a,d") * t2["abab"]("d,c,j,i") * r1_1["aa"]("R,b,k");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += l1_1["aa"]("L,j,a") * r2_1["abab"]("R,b,c,j,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,k,m,a,e") * t1_1["bb"]("e,i") * t1["bb"]("c,m") * r1["aa"]("R,b,k");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["128_aa_LLov"]("L,R,j,a") * t2["abab"]("b,c,j,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,k,m,a,e") * t1_1["bb"]("e,i") * t1["bb"]("c,m") * t1["aa"]("b,k") * r0("R");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += tmps_["130_aa_LLov"]("L,R,j,a") * t2["abab"]("b,c,j,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += tmps_["118_aa_LLov"]("L,R,j,a") * t2_1["abab"]("b,c,j,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["117_aa_LLov"]("L,R,j,a") * t2_1["abab"]("b,c,j,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += l1_1["aa"]("L,j,a") * t2_1["abab"]("b,c,j,i") * r0("R");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,j,l,a,e") * t1_1["bb"]("e,i") * r1["bb"]("R,c,l") * t1["aa"]("b,j");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,j,l,a,e") * t2_1["abab"]("b,e,j,i") * r1["bb"]("R,c,l");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,j,l,a,e") * t1["bb"]("e,i") * t1["bb"]("c,l") * t1_1["aa"]("b,j") * r0("R");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,j,l,a,e") * t2["abab"]("b,e,j,i") * r1_1["bb"]("R,c,l");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2["abab"]("L,j,l,a,e") * t1["bb"]("e,i") * r1["bb"]("R,c,l") * t1["aa"]("b,j");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,k,m,a,e") * r1["bb"]("R,e,i") * t1_1["bb"]("c,m") * t1["aa"]("b,k");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += l2_1["aaaa"]("L,k,j,a,d") * t2_1["abab"]("d,c,j,i") * r1["aa"]("R,b,k");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,k,m,a,e") * t2_1["bbbb"]("c,e,m,i") * r1["aa"]("R,b,k");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2["abab"]("L,k,m,a,e") * r1["bb"]("R,e,i") * t1["bb"]("c,m") * t1["aa"]("b,k");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= l2_1["abab"]("L,k,m,a,e") * t1["bb"]("e,i") * t1_1["bb"]("c,m") * t1["aa"]("b,k") * r0("R");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += tmps_["129_aa_LLov"]("L,R,j,a") * t2["abab"]("b,c,j,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += 0.50 * tmps_["247_aa_Lvv"]("L,a,b") * tmps_["224_bb_Lvo"]("R,c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += t1_1["aa"]("b,j") * tmps_["55_bbaa_Lvoov"]("L,c,i,j,a") * r0("R");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["75_baab_LLovvo"]("L,R,m,a,b,i") * t1["bb"]("c,m");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += 0.50 * tmps_["248_aa_Lvv"]("L,b,a") * tmps_["224_bb_Lvo"]("R,c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= (l2["aaaa"]("L,k,j,a,d") * r2["abab"]("R,d,c,k,i") + l2["abab"]("L,j,l,a,e") * r2["bbbb"]("R,c,e,l,i") + l2_1["aaaa"]("L,k,j,a,d") * r2_1["abab"]("R,d,c,k,i") + l2_1["abab"]("L,j,l,a,e") * r2_1["bbbb"]("R,c,e,l,i")) * t1["aa"]("b,j");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= t1["aa"]("b,k") * tmps_["209_bbaa_Lvoov"]("L,c,i,k,a") * r0("R");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += tmps_["258_aa_Lvv"]("L,a,b") * r1["bb"]("R,c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += tmps_["64_aa_LLvv"]("L,R,a,b") * t1_1["bb"]("c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["15_abba_Lvoov"]("L,b,i,l,a") * t1["bb"]("c,l") * r0("R");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += 0.50 * tmps_["246_aa_Lvv"]("L,a,b") * tmps_["226_bb_Lvo"]("R,c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= t1["aa"]("b,k") * tmps_["284_bbaa_Lvoov"]("L,c,i,k,a") * r0("R");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += 0.50 * tmps_["248_aa_Lvv"]("L,b,a") * r1["bb"]("R,c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += tmps_["258_aa_Lvv"]("L,a,b") * tmps_["224_bb_Lvo"]("R,c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= r1["bb"]("R,e,i") * tmps_["287_abab_Lvvvv"]("L,a,e,b,c");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["13_abba_Lvoov"]("L,b,i,m,a") * t1_1["bb"]("c,m") * r0("R");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += tmps_["29_aa_Lvv"]("L,a,b") * tmps_["226_bb_Lvo"]("R,c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= tmps_["14_abba_Lvoov"]("L,b,i,l,a") * t1["bb"]("c,l") * r0("R");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= t1_1["aa"]("b,j") * tmps_["93_bbaa_Lvoov"]("L,c,i,j,a") * r0("R");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += tmps_["252_aa_Lvv"]("L,a,b") * tmps_["226_bb_Lvo"]("R,c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= t1["aa"]("b,k") * tmps_["93_bbaa_Lvoov"]("L,c,i,k,a") * r0_1("R");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += tmps_["29_aa_Lvv"]("L,a,b") * tmps_["225_bb_Lvo"]("R,c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += tmps_["256_aa_LLvv"]("L,R,a,b") * t1["bb"]("c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= r0("R") * tmps_["289_aabb_Lvvvo"]("L,a,b,c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= r0("R") * tmps_["290_aabb_Lvvvo"]("L,a,b,c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += r1_1["bb"]("R,c,i") * tmps_["252_aa_Lvv"]("L,a,b");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= t1_1["aa"]("b,j") * tmps_["281_aabb_LLovvo"]("L,R,j,a,c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= r0_1("R") * tmps_["286_abba_Lvvov"]("L,a,c,i,b");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += 0.50 * r1_1["bb"]("R,c,i") * tmps_["246_aa_Lvv"]("L,a,b");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= t1_1["aa"]("b,j") * tmps_["282_aabb_LLovvo"]("L,R,j,a,c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= r0("R") * tmps_["285_abba_Lvvov"]("L,a,c,i,b");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += r1["bb"]("R,c,i") * tmps_["249_aa_Lvv"]("L,a,b");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += r1_1["bb"]("R,c,i") * tmps_["29_aa_Lvv"]("L,a,b");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += t1_1["bb"]("c,i") * tmps_["251_aa_LLvv"]("L,R,a,b");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= t1_1["bb"]("c,m") * tmps_["68_baab_LLovvo"]("L,R,m,a,b,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= r0_1("R") * tmps_["288_aabb_Lvvvo"]("L,a,b,c,i");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += 0.50 * r1["bb"]("R,c,i") * tmps_["247_aa_Lvv"]("L,a,b");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") -= t1["bb"]("c,i") * tmps_["259_aa_LLvv"]("L,R,b,a");
        tmps_["291_aabb_LLvvvo"]("L,R,a,b,c,i") += t1_1["bb"]("c,i") * tmps_["260_aa_LLvv"]("L,R,a,b");
        tmps_["290_aabb_Lvvvo"].~TArrayD();
        tmps_["289_aabb_Lvvvo"].~TArrayD();
        tmps_["288_aabb_Lvvvo"].~TArrayD();
        tmps_["287_abab_Lvvvv"].~TArrayD();
        tmps_["286_abba_Lvvov"].~TArrayD();
        tmps_["285_abba_Lvvov"].~TArrayD();
        tmps_["284_bbaa_Lvoov"].~TArrayD();
        tmps_["283_abab_Lvvvv"].~TArrayD();
        tmps_["282_aabb_LLovvo"].~TArrayD();
        tmps_["281_aabb_LLovvo"].~TArrayD();
        tmps_["260_aa_LLvv"].~TArrayD();
        tmps_["259_aa_LLvv"].~TArrayD();
        tmps_["258_aa_Lvv"].~TArrayD();
        tmps_["257_aa_LLvv"].~TArrayD();
        tmps_["252_aa_Lvv"].~TArrayD();
        tmps_["251_aa_LLvv"].~TArrayD();
        tmps_["209_bbaa_Lvoov"].~TArrayD();
        tmps_["130_aa_LLov"].~TArrayD();
        tmps_["129_aa_LLov"].~TArrayD();
        tmps_["128_aa_LLov"].~TArrayD();
        tmps_["118_aa_LLov"].~TArrayD();
        tmps_["93_bbaa_Lvoov"].~TArrayD();
        tmps_["75_baab_LLovvo"].~TArrayD();
        tmps_["68_baab_LLovvo"].~TArrayD();
        tmps_["55_bbaa_Lvoov"].~TArrayD();

        // D_ovvv_abab += -1.00 r1_bb(d,k) t1_aa(c,j) t1_1_bb(b,i) l2_1_abab(j,k,a,d)
        //             += +1.00 r1_aa(d,k) t1_aa(c,j) t1_1_bb(b,i) l2_1_aaaa(k,j,a,d)
        //             += +1.00 r1_aa(d,k) t1_bb(b,i) t1_aa(c,j) l2_aaaa(k,j,a,d)
        //             += +1.00 r1_aa(d,k) t1_bb(b,i) t1_1_aa(c,j) l2_1_aaaa(k,j,a,d)
        //             += -1.00 r1_1_bb(d,k) t1_bb(b,i) t1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += -1.00 r1_bb(d,k) t1_bb(b,i) t1_aa(c,j) l2_abab(j,k,a,d)
        //             += +1.00 r1_1_aa(d,k) t1_bb(b,i) t1_aa(c,j) l2_1_aaaa(k,j,a,d)
        //             += -1.00 r1_bb(d,k) t1_bb(b,i) t1_1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += -0.50 r0_1 t1_bb(b,i) t2_aaaa(c,d,j,k) l2_1_aaaa(j,k,a,d)
        //             += -1.00 r0_1 t1_bb(b,i) t1_aa(c,j) l1_1_aa(j,a)
        //             += -0.50 r1_bb(b,i) t2_abab(c,d,j,k) l2_abab(j,k,a,d)
        //             += -0.50 r1_bb(b,i) t2_abab(c,d,k,j) l2_abab(k,j,a,d)
        //             += +1.00 r0_1 t1_bb(b,k) t2_abab(c,d,j,i) l2_1_abab(j,k,a,d)
        //             += -0.50 r0 t1_bb(b,i) t2_abab(c,d,j,k) l2_abab(j,k,a,d)
        //             += -0.50 r0 t1_bb(b,i) t2_abab(c,d,k,j) l2_abab(k,j,a,d)
        //             += -0.50 r0 t1_bb(b,i) t2_1_abab(c,d,j,k) l2_1_abab(j,k,a,d)
        //             += -0.50 r0 t1_bb(b,i) t2_1_abab(c,d,k,j) l2_1_abab(k,j,a,d)
        //             += -1.00 r1_1_aa(c,j) t1_bb(b,i) l1_1_aa(j,a)
        //             += -1.00 r1_aa(c,j) t1_bb(b,i) l1_aa(j,a)
        //             += -1.00 r2_abab(c,b,j,i) l1_aa(j,a)
        //             += +1.00 r1_bb(d,i) t1_bb(b,k) t1_1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += +1.00 r1_bb(b,k) t1_bb(d,i) t1_1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += +0.50 r1_1_bb(d,i) t2_abab(c,b,j,k) l2_1_abab(j,k,a,d)
        //             += +0.50 r1_1_bb(d,i) t2_abab(c,b,k,j) l2_1_abab(k,j,a,d)
        //             += -1.00 r0_1 t2_abab(c,b,j,i) l1_1_aa(j,a)
        //             += +1.00 r1_aa(c,k) t1_bb(d,i) t1_1_bb(b,j) l2_1_abab(k,j,a,d)
        //             += +1.00 r1_aa(c,k) t2_bbbb(b,d,j,i) l2_abab(k,j,a,d)
        //             += +1.00 r0 t1_bb(b,j) t1_aa(c,k) t1_bb(d,i) l2_abab(k,j,a,d)
        //             += +1.00 r1_1_bb(d,i) t1_bb(b,j) t1_aa(c,k) l2_1_abab(k,j,a,d)
        //             += +1.00 r1_1_aa(c,k) t1_bb(b,j) t1_bb(d,i) l2_1_abab(k,j,a,d)
        //             += +1.00 r1_1_bb(b,k) t1_aa(c,j) t1_bb(d,i) l2_1_abab(j,k,a,d)
        //             += +1.00 r1_1_aa(c,k) t2_bbbb(b,d,j,i) l2_1_abab(k,j,a,d)
        //             += +1.00 r0_1 t1_bb(b,j) t1_aa(c,k) t1_bb(d,i) l2_1_abab(k,j,a,d)
        //             += +1.00 r1_bb(b,k) t2_abab(c,d,j,i) l2_abab(j,k,a,d)
        //             += -1.00 r1_aa(c,k) t2_abab(d,b,j,i) l2_aaaa(k,j,a,d)
        //             += +1.00 r1_aa(c,k) t1_bb(b,j) t1_bb(d,i) l2_abab(k,j,a,d)
        //             += -1.00 r0 t2_abab(c,b,j,i) l1_aa(j,a)
        //             += +1.00 r1_aa(d,k) t2_abab(c,b,j,i) l2_aaaa(k,j,a,d)
        //             += -1.00 r1_1_aa(c,k) t2_abab(d,b,j,i) l2_1_aaaa(k,j,a,d)
        //             += -1.00 r2_1_abab(c,b,j,i) l1_1_aa(j,a)
        //             += +1.00 r1_aa(c,k) t1_bb(b,j) t1_1_bb(d,i) l2_1_abab(k,j,a,d)
        //             += +1.00 r1_1_aa(d,k) t2_abab(c,b,j,i) l2_1_aaaa(k,j,a,d)
        //             += +1.00 r0 t1_bb(b,j) t1_aa(c,k) t1_1_bb(d,i) l2_1_abab(k,j,a,d)
        //             += -1.00 r1_1_bb(d,k) t2_abab(c,b,j,i) l2_1_abab(j,k,a,d)
        //             += -1.00 r1_bb(d,k) t2_1_abab(c,b,j,i) l2_1_abab(j,k,a,d)
        //             += +1.00 r1_aa(d,k) t2_1_abab(c,b,j,i) l2_1_aaaa(k,j,a,d)
        //             += -1.00 r0 t2_1_abab(c,b,j,i) l1_1_aa(j,a)
        //             += +1.00 r1_bb(b,k) t1_aa(c,j) t1_1_bb(d,i) l2_1_abab(j,k,a,d)
        //             += +1.00 r1_bb(b,k) t2_1_abab(c,d,j,i) l2_1_abab(j,k,a,d)
        //             += +1.00 r0 t1_bb(b,k) t1_bb(d,i) t1_1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += +1.00 r1_1_bb(b,k) t2_abab(c,d,j,i) l2_1_abab(j,k,a,d)
        //             += +1.00 r1_bb(b,k) t1_aa(c,j) t1_bb(d,i) l2_abab(j,k,a,d)
        //             += +1.00 r1_bb(d,i) t1_aa(c,k) t1_1_bb(b,j) l2_1_abab(k,j,a,d)
        //             += -1.00 r1_aa(c,k) t2_1_abab(d,b,j,i) l2_1_aaaa(k,j,a,d)
        //             += +1.00 r1_aa(c,k) t2_1_bbbb(b,d,j,i) l2_1_abab(k,j,a,d)
        //             += +1.00 r1_bb(d,i) t1_bb(b,j) t1_aa(c,k) l2_abab(k,j,a,d)
        //             += +1.00 r0 t1_aa(c,k) t1_bb(d,i) t1_1_bb(b,j) l2_1_abab(k,j,a,d)
        //             += -1.00 r1_bb(d,k) t2_abab(c,b,j,i) l2_abab(j,k,a,d)
        //             += -0.50 r0 t1_bb(b,i) t2_1_aaaa(c,d,j,k) l2_1_aaaa(j,k,a,d)
        //             += -1.00 r0 t2_abab(d,b,k,i) t1_1_aa(c,j) l2_1_aaaa(j,k,a,d)
        //             += +1.00 r2_abab(c,d,k,i) t1_bb(b,j) l2_abab(k,j,a,d)
        //             += +1.00 r2_1_abab(c,d,k,i) t1_bb(b,j) l2_1_abab(k,j,a,d)
        //             += -0.50 r0 t1_bb(b,i) t2_aaaa(c,d,j,k) l2_aaaa(j,k,a,d)
        //             += +1.00 r2_abab(d,b,k,i) t1_aa(c,j) l2_aaaa(k,j,a,d)
        //             += +1.00 r2_bbbb(b,d,k,i) t1_aa(c,j) l2_abab(j,k,a,d)
        //             += +1.00 r2_1_abab(d,b,k,i) t1_aa(c,j) l2_1_aaaa(k,j,a,d)
        //             += +1.00 r2_1_bbbb(b,d,k,i) t1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += +1.00 r0 t2_bbbb(b,d,j,i) t1_aa(c,k) l2_abab(k,j,a,d)
        //             += +1.00 r0 t1_aa(c,k) t2_1_bbbb(b,d,j,i) l2_1_abab(k,j,a,d)
        //             += -1.00 r1_bb(b,i) t1_1_aa(c,j) l1_1_aa(j,a)
        //             += -1.00 r1_bb(b,i) t1_aa(c,j) l1_aa(j,a)
        //             += -0.50 r2_abab(c,d,k,j) t1_1_bb(b,i) l2_1_abab(k,j,a,d)
        //             += -0.50 r2_abab(c,d,j,k) t1_1_bb(b,i) l2_1_abab(j,k,a,d)
        //             += -0.50 r2_aaaa(c,d,k,j) t1_1_bb(b,i) l2_1_aaaa(k,j,a,d)
        //             += +1.00 r0 t1_bb(b,k) t2_1_abab(c,d,j,i) l2_1_abab(j,k,a,d)
        //             += -0.50 r0 t2_aaaa(c,d,j,k) t1_1_bb(b,i) l2_1_aaaa(j,k,a,d)
        //             += +1.00 r0 t2_abab(d,b,j,i) t1_aa(c,k) l2_aaaa(j,k,a,d)
        //             += -0.50 r1_bb(b,i) t2_aaaa(c,d,j,k) l2_aaaa(j,k,a,d)
        //             += -1.00 r0 t1_bb(b,i) t1_1_aa(c,j) l1_1_aa(j,a)
        //             += -1.00 r0 t1_bb(b,i) t1_aa(c,j) l1_aa(j,a)
        //             += +0.50 r1_bb(d,i) t2_abab(c,b,j,k) l2_abab(j,k,a,d)
        //             += +0.50 r1_bb(d,i) t2_abab(c,b,k,j) l2_abab(k,j,a,d)
        //             += +0.50 r1_bb(d,i) t2_1_abab(c,b,j,k) l2_1_abab(j,k,a,d)
        //             += +0.50 r1_bb(d,i) t2_1_abab(c,b,k,j) l2_1_abab(k,j,a,d)
        //             += +1.00 r0 t2_abab(c,d,k,i) t1_1_bb(b,j) l2_1_abab(k,j,a,d)
        //             += -0.50 r0 t2_abab(c,d,j,k) t1_1_bb(b,i) l2_1_abab(j,k,a,d)
        //             += -0.50 r0 t2_abab(c,d,k,j) t1_1_bb(b,i) l2_1_abab(k,j,a,d)
        //             += +1.00 r0 t1_bb(b,k) t2_abab(c,d,j,i) l2_abab(j,k,a,d)
        //             += +1.00 r0 t2_bbbb(b,d,k,i) t1_1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += -1.00 r0 t1_aa(c,j) t1_1_bb(b,i) l1_1_aa(j,a)
        //             += +1.00 r0_1 t2_bbbb(b,d,j,i) t1_aa(c,k) l2_1_abab(k,j,a,d)
        //             += -0.50 r0_1 t1_bb(b,i) t2_abab(c,d,j,k) l2_1_abab(j,k,a,d)
        //             += -0.50 r0_1 t1_bb(b,i) t2_abab(c,d,k,j) l2_1_abab(k,j,a,d)
        //             += -0.50 r2_1_abab(c,d,k,j) t1_bb(b,i) l2_1_abab(k,j,a,d)
        //             += -0.50 r2_1_abab(c,d,j,k) t1_bb(b,i) l2_1_abab(j,k,a,d)
        //             += -0.50 r2_aaaa(c,d,k,j) t1_bb(b,i) l2_aaaa(k,j,a,d)
        //             += -0.50 r2_abab(c,d,k,j) t1_bb(b,i) l2_abab(k,j,a,d)
        //             += -0.50 r2_abab(c,d,j,k) t1_bb(b,i) l2_abab(j,k,a,d)
        //             += -0.50 r2_1_aaaa(c,d,k,j) t1_bb(b,i) l2_1_aaaa(k,j,a,d)
        //             += +0.50 r0 t2_abab(c,b,j,k) t1_bb(d,i) l2_abab(j,k,a,d)
        //             += +0.50 r0 t2_abab(c,b,k,j) t1_bb(d,i) l2_abab(k,j,a,d)
        //             += +0.50 r0 t1_bb(d,i) t2_1_abab(c,b,j,k) l2_1_abab(j,k,a,d)
        //             += +0.50 r0 t1_bb(d,i) t2_1_abab(c,b,k,j) l2_1_abab(k,j,a,d)
        //             += +0.50 r0 t2_abab(c,b,j,k) t1_1_bb(d,i) l2_1_abab(j,k,a,d)
        //             += +0.50 r0 t2_abab(c,b,k,j) t1_1_bb(d,i) l2_1_abab(k,j,a,d)
        //             += -1.00 r1_1_bb(b,i) t1_aa(c,j) l1_1_aa(j,a)
        //             += +1.00 r2_abab(d,b,k,i) t1_1_aa(c,j) l2_1_aaaa(k,j,a,d)
        //             += +1.00 r0_1 t2_abab(d,b,j,i) t1_aa(c,k) l2_1_aaaa(j,k,a,d)
        //             += -0.50 r1_1_bb(b,i) t2_aaaa(c,d,j,k) l2_1_aaaa(j,k,a,d)
        //             += +1.00 r2_bbbb(b,d,k,i) t1_1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += +1.00 r0 t1_aa(c,k) t2_1_abab(d,b,j,i) l2_1_aaaa(j,k,a,d)
        //             += -0.50 r1_bb(b,i) t2_1_abab(c,d,j,k) l2_1_abab(j,k,a,d)
        //             += -0.50 r1_bb(b,i) t2_1_abab(c,d,k,j) l2_1_abab(k,j,a,d)
        //             += -0.50 r1_1_bb(b,i) t2_abab(c,d,j,k) l2_1_abab(j,k,a,d)
        //             += -0.50 r1_1_bb(b,i) t2_abab(c,d,k,j) l2_1_abab(k,j,a,d)
        //             += -1.00 r1_aa(c,j) t1_1_bb(b,i) l1_1_aa(j,a)
        //             += +1.00 r2_abab(c,d,k,i) t1_1_bb(b,j) l2_1_abab(k,j,a,d)
        //             += +0.50 r0_1 t2_abab(c,b,j,k) t1_bb(d,i) l2_1_abab(j,k,a,d)
        //             += +0.50 r0_1 t2_abab(c,b,k,j) t1_bb(d,i) l2_1_abab(k,j,a,d)
        //             += -0.50 r1_bb(b,i) t2_1_aaaa(c,d,j,k) l2_1_aaaa(j,k,a,d)
        D_ovvv_abab("L,R,i,a,b,c") -= tmps_["291_aabb_LLvvvo"]("L,R,a,c,b,i");

        // D_vovv_abab += +1.00 r1_bb(d,k) t1_aa(c,j) t1_1_bb(b,i) l2_1_abab(j,k,a,d)
        //             += -1.00 r1_aa(d,k) t1_aa(c,j) t1_1_bb(b,i) l2_1_aaaa(k,j,a,d)
        //             += -1.00 r1_aa(d,k) t1_bb(b,i) t1_aa(c,j) l2_aaaa(k,j,a,d)
        //             += -1.00 r1_aa(d,k) t1_bb(b,i) t1_1_aa(c,j) l2_1_aaaa(k,j,a,d)
        //             += +1.00 r1_1_bb(d,k) t1_bb(b,i) t1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += +1.00 r1_bb(d,k) t1_bb(b,i) t1_aa(c,j) l2_abab(j,k,a,d)
        //             += -1.00 r1_1_aa(d,k) t1_bb(b,i) t1_aa(c,j) l2_1_aaaa(k,j,a,d)
        //             += +1.00 r1_bb(d,k) t1_bb(b,i) t1_1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += +0.50 r0_1 t1_bb(b,i) t2_aaaa(c,d,j,k) l2_1_aaaa(j,k,a,d)
        //             += +1.00 r0_1 t1_bb(b,i) t1_aa(c,j) l1_1_aa(j,a)
        //             += +0.50 r1_bb(b,i) t2_abab(c,d,j,k) l2_abab(j,k,a,d)
        //             += +0.50 r1_bb(b,i) t2_abab(c,d,k,j) l2_abab(k,j,a,d)
        //             += -1.00 r0_1 t1_bb(b,k) t2_abab(c,d,j,i) l2_1_abab(j,k,a,d)
        //             += +0.50 r0 t1_bb(b,i) t2_abab(c,d,j,k) l2_abab(j,k,a,d)
        //             += +0.50 r0 t1_bb(b,i) t2_abab(c,d,k,j) l2_abab(k,j,a,d)
        //             += +0.50 r0 t1_bb(b,i) t2_1_abab(c,d,j,k) l2_1_abab(j,k,a,d)
        //             += +0.50 r0 t1_bb(b,i) t2_1_abab(c,d,k,j) l2_1_abab(k,j,a,d)
        //             += +1.00 r1_1_aa(c,j) t1_bb(b,i) l1_1_aa(j,a)
        //             += +1.00 r1_aa(c,j) t1_bb(b,i) l1_aa(j,a)
        //             += +1.00 r2_abab(c,b,j,i) l1_aa(j,a)
        //             += -1.00 r1_bb(d,i) t1_bb(b,k) t1_1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += -1.00 r1_bb(b,k) t1_bb(d,i) t1_1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += -0.50 r1_1_bb(d,i) t2_abab(c,b,j,k) l2_1_abab(j,k,a,d)
        //             += -0.50 r1_1_bb(d,i) t2_abab(c,b,k,j) l2_1_abab(k,j,a,d)
        //             += +1.00 r0_1 t2_abab(c,b,j,i) l1_1_aa(j,a)
        //             += -1.00 r1_aa(c,k) t1_bb(d,i) t1_1_bb(b,j) l2_1_abab(k,j,a,d)
        //             += -1.00 r1_aa(c,k) t2_bbbb(b,d,j,i) l2_abab(k,j,a,d)
        //             += -1.00 r0 t1_bb(b,j) t1_aa(c,k) t1_bb(d,i) l2_abab(k,j,a,d)
        //             += -1.00 r1_1_bb(d,i) t1_bb(b,j) t1_aa(c,k) l2_1_abab(k,j,a,d)
        //             += -1.00 r1_1_aa(c,k) t1_bb(b,j) t1_bb(d,i) l2_1_abab(k,j,a,d)
        //             += -1.00 r1_1_bb(b,k) t1_aa(c,j) t1_bb(d,i) l2_1_abab(j,k,a,d)
        //             += -1.00 r1_1_aa(c,k) t2_bbbb(b,d,j,i) l2_1_abab(k,j,a,d)
        //             += -1.00 r0_1 t1_bb(b,j) t1_aa(c,k) t1_bb(d,i) l2_1_abab(k,j,a,d)
        //             += -1.00 r1_bb(b,k) t2_abab(c,d,j,i) l2_abab(j,k,a,d)
        //             += +1.00 r1_aa(c,k) t2_abab(d,b,j,i) l2_aaaa(k,j,a,d)
        //             += -1.00 r1_aa(c,k) t1_bb(b,j) t1_bb(d,i) l2_abab(k,j,a,d)
        //             += +1.00 r0 t2_abab(c,b,j,i) l1_aa(j,a)
        //             += -1.00 r1_aa(d,k) t2_abab(c,b,j,i) l2_aaaa(k,j,a,d)
        //             += +1.00 r1_1_aa(c,k) t2_abab(d,b,j,i) l2_1_aaaa(k,j,a,d)
        //             += +1.00 r2_1_abab(c,b,j,i) l1_1_aa(j,a)
        //             += -1.00 r1_aa(c,k) t1_bb(b,j) t1_1_bb(d,i) l2_1_abab(k,j,a,d)
        //             += -1.00 r1_1_aa(d,k) t2_abab(c,b,j,i) l2_1_aaaa(k,j,a,d)
        //             += -1.00 r0 t1_bb(b,j) t1_aa(c,k) t1_1_bb(d,i) l2_1_abab(k,j,a,d)
        //             += +1.00 r1_1_bb(d,k) t2_abab(c,b,j,i) l2_1_abab(j,k,a,d)
        //             += +1.00 r1_bb(d,k) t2_1_abab(c,b,j,i) l2_1_abab(j,k,a,d)
        //             += -1.00 r1_aa(d,k) t2_1_abab(c,b,j,i) l2_1_aaaa(k,j,a,d)
        //             += +1.00 r0 t2_1_abab(c,b,j,i) l1_1_aa(j,a)
        //             += -1.00 r1_bb(b,k) t1_aa(c,j) t1_1_bb(d,i) l2_1_abab(j,k,a,d)
        //             += -1.00 r1_bb(b,k) t2_1_abab(c,d,j,i) l2_1_abab(j,k,a,d)
        //             += -1.00 r0 t1_bb(b,k) t1_bb(d,i) t1_1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += -1.00 r1_1_bb(b,k) t2_abab(c,d,j,i) l2_1_abab(j,k,a,d)
        //             += -1.00 r1_bb(b,k) t1_aa(c,j) t1_bb(d,i) l2_abab(j,k,a,d)
        //             += -1.00 r1_bb(d,i) t1_aa(c,k) t1_1_bb(b,j) l2_1_abab(k,j,a,d)
        //             += +1.00 r1_aa(c,k) t2_1_abab(d,b,j,i) l2_1_aaaa(k,j,a,d)
        //             += -1.00 r1_aa(c,k) t2_1_bbbb(b,d,j,i) l2_1_abab(k,j,a,d)
        //             += -1.00 r1_bb(d,i) t1_bb(b,j) t1_aa(c,k) l2_abab(k,j,a,d)
        //             += -1.00 r0 t1_aa(c,k) t1_bb(d,i) t1_1_bb(b,j) l2_1_abab(k,j,a,d)
        //             += +1.00 r1_bb(d,k) t2_abab(c,b,j,i) l2_abab(j,k,a,d)
        //             += +0.50 r0 t1_bb(b,i) t2_1_aaaa(c,d,j,k) l2_1_aaaa(j,k,a,d)
        //             += +1.00 r0 t2_abab(d,b,k,i) t1_1_aa(c,j) l2_1_aaaa(j,k,a,d)
        //             += -1.00 r2_abab(c,d,k,i) t1_bb(b,j) l2_abab(k,j,a,d)
        //             += -1.00 r2_1_abab(c,d,k,i) t1_bb(b,j) l2_1_abab(k,j,a,d)
        //             += +0.50 r0 t1_bb(b,i) t2_aaaa(c,d,j,k) l2_aaaa(j,k,a,d)
        //             += -1.00 r2_abab(d,b,k,i) t1_aa(c,j) l2_aaaa(k,j,a,d)
        //             += -1.00 r2_bbbb(b,d,k,i) t1_aa(c,j) l2_abab(j,k,a,d)
        //             += -1.00 r2_1_abab(d,b,k,i) t1_aa(c,j) l2_1_aaaa(k,j,a,d)
        //             += -1.00 r2_1_bbbb(b,d,k,i) t1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += -1.00 r0 t2_bbbb(b,d,j,i) t1_aa(c,k) l2_abab(k,j,a,d)
        //             += -1.00 r0 t1_aa(c,k) t2_1_bbbb(b,d,j,i) l2_1_abab(k,j,a,d)
        //             += +1.00 r1_bb(b,i) t1_1_aa(c,j) l1_1_aa(j,a)
        //             += +1.00 r1_bb(b,i) t1_aa(c,j) l1_aa(j,a)
        //             += +0.50 r2_abab(c,d,k,j) t1_1_bb(b,i) l2_1_abab(k,j,a,d)
        //             += +0.50 r2_abab(c,d,j,k) t1_1_bb(b,i) l2_1_abab(j,k,a,d)
        //             += +0.50 r2_aaaa(c,d,k,j) t1_1_bb(b,i) l2_1_aaaa(k,j,a,d)
        //             += -1.00 r0 t1_bb(b,k) t2_1_abab(c,d,j,i) l2_1_abab(j,k,a,d)
        //             += +0.50 r0 t2_aaaa(c,d,j,k) t1_1_bb(b,i) l2_1_aaaa(j,k,a,d)
        //             += -1.00 r0 t2_abab(d,b,j,i) t1_aa(c,k) l2_aaaa(j,k,a,d)
        //             += +0.50 r1_bb(b,i) t2_aaaa(c,d,j,k) l2_aaaa(j,k,a,d)
        //             += +1.00 r0 t1_bb(b,i) t1_1_aa(c,j) l1_1_aa(j,a)
        //             += +1.00 r0 t1_bb(b,i) t1_aa(c,j) l1_aa(j,a)
        //             += -0.50 r1_bb(d,i) t2_abab(c,b,j,k) l2_abab(j,k,a,d)
        //             += -0.50 r1_bb(d,i) t2_abab(c,b,k,j) l2_abab(k,j,a,d)
        //             += -0.50 r1_bb(d,i) t2_1_abab(c,b,j,k) l2_1_abab(j,k,a,d)
        //             += -0.50 r1_bb(d,i) t2_1_abab(c,b,k,j) l2_1_abab(k,j,a,d)
        //             += -1.00 r0 t2_abab(c,d,k,i) t1_1_bb(b,j) l2_1_abab(k,j,a,d)
        //             += +0.50 r0 t2_abab(c,d,j,k) t1_1_bb(b,i) l2_1_abab(j,k,a,d)
        //             += +0.50 r0 t2_abab(c,d,k,j) t1_1_bb(b,i) l2_1_abab(k,j,a,d)
        //             += -1.00 r0 t1_bb(b,k) t2_abab(c,d,j,i) l2_abab(j,k,a,d)
        //             += -1.00 r0 t2_bbbb(b,d,k,i) t1_1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += +1.00 r0 t1_aa(c,j) t1_1_bb(b,i) l1_1_aa(j,a)
        //             += -1.00 r0_1 t2_bbbb(b,d,j,i) t1_aa(c,k) l2_1_abab(k,j,a,d)
        //             += +0.50 r0_1 t1_bb(b,i) t2_abab(c,d,j,k) l2_1_abab(j,k,a,d)
        //             += +0.50 r0_1 t1_bb(b,i) t2_abab(c,d,k,j) l2_1_abab(k,j,a,d)
        //             += +0.50 r2_1_abab(c,d,k,j) t1_bb(b,i) l2_1_abab(k,j,a,d)
        //             += +0.50 r2_1_abab(c,d,j,k) t1_bb(b,i) l2_1_abab(j,k,a,d)
        //             += +0.50 r2_aaaa(c,d,k,j) t1_bb(b,i) l2_aaaa(k,j,a,d)
        //             += +0.50 r2_abab(c,d,k,j) t1_bb(b,i) l2_abab(k,j,a,d)
        //             += +0.50 r2_abab(c,d,j,k) t1_bb(b,i) l2_abab(j,k,a,d)
        //             += +0.50 r2_1_aaaa(c,d,k,j) t1_bb(b,i) l2_1_aaaa(k,j,a,d)
        //             += -0.50 r0 t2_abab(c,b,j,k) t1_bb(d,i) l2_abab(j,k,a,d)
        //             += -0.50 r0 t2_abab(c,b,k,j) t1_bb(d,i) l2_abab(k,j,a,d)
        //             += -0.50 r0 t1_bb(d,i) t2_1_abab(c,b,j,k) l2_1_abab(j,k,a,d)
        //             += -0.50 r0 t1_bb(d,i) t2_1_abab(c,b,k,j) l2_1_abab(k,j,a,d)
        //             += -0.50 r0 t2_abab(c,b,j,k) t1_1_bb(d,i) l2_1_abab(j,k,a,d)
        //             += -0.50 r0 t2_abab(c,b,k,j) t1_1_bb(d,i) l2_1_abab(k,j,a,d)
        //             += +1.00 r1_1_bb(b,i) t1_aa(c,j) l1_1_aa(j,a)
        //             += -1.00 r2_abab(d,b,k,i) t1_1_aa(c,j) l2_1_aaaa(k,j,a,d)
        //             += -1.00 r0_1 t2_abab(d,b,j,i) t1_aa(c,k) l2_1_aaaa(j,k,a,d)
        //             += +0.50 r1_1_bb(b,i) t2_aaaa(c,d,j,k) l2_1_aaaa(j,k,a,d)
        //             += -1.00 r2_bbbb(b,d,k,i) t1_1_aa(c,j) l2_1_abab(j,k,a,d)
        //             += -1.00 r0 t1_aa(c,k) t2_1_abab(d,b,j,i) l2_1_aaaa(j,k,a,d)
        //             += +0.50 r1_bb(b,i) t2_1_abab(c,d,j,k) l2_1_abab(j,k,a,d)
        //             += +0.50 r1_bb(b,i) t2_1_abab(c,d,k,j) l2_1_abab(k,j,a,d)
        //             += +0.50 r1_1_bb(b,i) t2_abab(c,d,j,k) l2_1_abab(j,k,a,d)
        //             += +0.50 r1_1_bb(b,i) t2_abab(c,d,k,j) l2_1_abab(k,j,a,d)
        //             += +1.00 r1_aa(c,j) t1_1_bb(b,i) l1_1_aa(j,a)
        //             += -1.00 r2_abab(c,d,k,i) t1_1_bb(b,j) l2_1_abab(k,j,a,d)
        //             += -0.50 r0_1 t2_abab(c,b,j,k) t1_bb(d,i) l2_1_abab(j,k,a,d)
        //             += -0.50 r0_1 t2_abab(c,b,k,j) t1_bb(d,i) l2_1_abab(k,j,a,d)
        //             += +0.50 r1_bb(b,i) t2_1_aaaa(c,d,j,k) l2_1_aaaa(j,k,a,d)
        D_vovv_abab("L,R,a,i,b,c") += tmps_["291_aabb_LLvvvo"]("L,R,a,c,b,i");
        tmps_["291_aabb_LLvvvo"].~TArrayD();

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["294_aaaa_Lvooo"]("L,a,i,j,k")  = tmps_["243_aaaa_Lvoov"]("L,a,i,j,b") * t1["aa"]("b,k");

        // D_oovv_aaaa += -1.00 P(i,j) P(a,b) r1_1_aa(a,l) t2_abab(b,c,i,k) t1_aa(d,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["294_aaaa_Lvooo"]("L,b,i,l,j") * r1_1["aa"]("R,a,l");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += -1.00 P(i,j) P(a,b) r0 t2_abab(b,c,i,l) t1_aa(d,j) t1_1_aa(a,k) l2_1_abab(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1_1["aa"]("a,k") * tmps_["294_aaaa_Lvooo"]("L,b,i,k,j") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r0_1 t2_abab(a,c,i,k) t1_aa(b,l) t1_aa(d,j) l2_1_abab(l,k,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["294_aaaa_Lvooo"]("L,a,i,l,j") * t1["aa"]("b,l") * r0_1("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o3v2L1
        //  mems: o3v1L1  = o3v1L1
        tmps_["295_aaaa_Lvooo"]("L,a,i,j,k")  = tmps_["240_aaaa_Lvoov"]("L,a,i,j,b") * t1["aa"]("b,k");
        tmps_["240_aaaa_Lvoov"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r0_1 t2_aaaa(a,c,k,i) t1_aa(b,l) t1_aa(d,j) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["295_aaaa_Lvooo"]("L,a,i,l,j") * t1["aa"]("b,l") * r0_1("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o4v1L1 o4v1L1
        //  mems: o3v1L1  = o4v0L1 o3v1L1
        tmps_["296_aaaa_Looov"]("L,i,j,k,a")  = t1_1["aa"]("b,i") * tmps_["22_aaaa_Looov"]("L,j,l,k,b") * t1["aa"]("a,l");

        // D_oovv_aaaa += +1.00 P(i,j) r0 t1_aa(a,k) t1_aa(b,l) t1_aa(d,j) t1_1_aa(c,i) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t1["aa"]("b,l") * tmps_["296_aaaa_Looov"]("L,i,j,l,a") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o2v2L1 o3v2L1 o3v2L1 o3v1L1
        //  mems: o3v1L1  = o2v2L1 o3v1L1 o3v1L1 o3v1L1
        tmps_["297_aaaa_Lvooo"]("L,a,i,j,k")  = (tmps_["244_aaaa_Lvoov"]("L,a,i,j,b") + tmps_["245_aaaa_Lvoov"]("L,a,i,j,b")) * t1["aa"]("b,k");
        tmps_["297_aaaa_Lvooo"]("L,a,i,j,k") -= t1_1["aa"]("c,i") * tmps_["243_aaaa_Lvoov"]("L,a,k,j,c");

        // D_oovv_aaaa += -1.00 P(i,j) P(a,b) r1_aa(a,l) t1_aa(d,j) t2_1_abab(b,c,i,k) l2_1_abab(l,k,d,c)
        //             += -1.00 P(i,j) P(a,b) r1_aa(a,l) t2_abab(b,c,i,k) t1_aa(d,j) l2_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r1_aa(a,l) t2_abab(b,d,j,k) t1_1_aa(c,i) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["297_aaaa_Lvooo"]("L,b,i,l,j") * r1["aa"]("R,a,l");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r0 t1_aa(b,l) t1_aa(d,j) t2_1_abab(a,c,i,k) l2_1_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r0 t2_abab(a,c,i,k) t1_aa(b,l) t1_aa(d,j) l2_abab(l,k,d,c)
        //             += -1.00 P(i,j) P(a,b) r0 t2_abab(a,d,j,k) t1_aa(b,l) t1_1_aa(c,i) l2_1_abab(l,k,c,d)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["297_aaaa_Lvooo"]("L,a,i,l,j") * t1["aa"]("b,l") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o3v1L1  = o2v2L1 o3v2L1 o3v2L1 o4v2L1 o3v1L1
        //  mems: o3v1L1  = o2v2L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1
        tmps_["298_aaaa_Lvooo"]("L,a,i,j,k")  = (tmps_["241_aaaa_Lvoov"]("L,a,i,j,b") + tmps_["242_aaaa_Lvoov"]("L,a,i,j,b")) * t1["aa"]("b,k");
        tmps_["298_aaaa_Lvooo"]("L,a,i,j,k") += l2_1["aaaa"]("L,l,j,b,c") * t1_1["aa"]("c,i") * t2["aaaa"]("a,b,l,k");
        tmps_["242_aaaa_Lvoov"].~TArrayD();
        tmps_["241_aaaa_Lvoov"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r0 t1_aa(b,l) t1_aa(d,j) t2_1_aaaa(a,c,k,i) l2_1_aaaa(k,l,d,c)
        //             += +1.00 P(i,j) P(a,b) r0 t2_aaaa(a,c,k,i) t1_aa(b,l) t1_aa(d,j) l2_aaaa(k,l,d,c)
        //             += +1.00 P(i,j) P(a,b) r0 t2_aaaa(a,d,k,j) t1_aa(b,l) t1_1_aa(c,i) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o3v2L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["298_aaaa_Lvooo"]("L,a,i,l,j") * t1["aa"]("b,l") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o2v2L1  = o3v3L1
        //  mems: o2v2L1  = o2v2L1
        tmps_["299_aabb_Lvoov"]("L,a,i,j,b")  = t2["abab"]("a,c,i,k") * l2["bbbb"]("L,j,k,c,b");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_abab(a,d,i,l) t2_abab(b,c,j,k) l2_bbbb(l,k,c,d)
        // flops: o2v2L2 += o3v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["299_aabb_Lvoov"]("L,b,j,l,d") * r2["abab"]("R,a,d,i,l");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +1.00 P(i,j) r0 t2_abab(a,c,i,k) t2_abab(b,d,j,l) l2_bbbb(k,l,d,c)
        // flops: o2v2L2 += o3v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = t2["abab"]("a,c,i,k") * tmps_["299_aabb_Lvoov"]("L,b,j,k,c") * r0("R");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["300_aa_Lvv"]("L,a,b")  = t2["aaaa"]("a,c,i,j") * l2["aaaa"]("L,i,j,c,b");

        // D_oovv_aaaa += -0.50 P(a,b) r2_aaaa(a,d,i,j) t2_aaaa(b,c,k,l) l2_aaaa(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["300_aa_Lvv"]("L,b,d") * r2["aaaa"]("R,a,d,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += -0.50 r0 t2_aaaa(a,c,i,j) t2_aaaa(b,d,k,l) l2_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v3L1 o2v2L2
        //  mems: o2v2L2 += o2v2L1 o2v2L2
        D_oovv_aaaa("L,R,i,j,a,b") -= 0.50 * t2["aaaa"]("a,c,i,j") * tmps_["300_aa_Lvv"]("L,b,c") * r0("R");

        // D_oovv_abab += -0.50 r2_abab(d,b,i,j) t2_aaaa(a,c,k,l) l2_aaaa(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * tmps_["300_aa_Lvv"]("L,a,d") * r2["abab"]("R,d,b,i,j");

        // flops: o0v2L1  = o2v3L1
        //  mems: o0v2L1  = o0v2L1
        tmps_["301_aa_Lvv"]("L,a,b")  = t2_1["aaaa"]("a,c,i,j") * l2_1["aaaa"]("L,i,j,c,b");

        // D_oovv_aaaa += -0.50 P(a,b) r2_aaaa(a,d,i,j) t2_1_aaaa(b,c,k,l) l2_1_aaaa(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["301_aa_Lvv"]("L,b,d") * r2["aaaa"]("R,a,d,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -0.50 r2_abab(d,b,i,j) t2_1_aaaa(a,c,k,l) l2_1_aaaa(k,l,c,d)
        // flops: o2v2L2 += o2v3L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * tmps_["301_aa_Lvv"]("L,a,d") * r2["abab"]("R,d,b,i,j");

        // flops: o1v1L2  = o2v2L2 o2v2L2 o2v2L2 o2v2L1 o2v2L1 o2v2L2 o1v1L2 o2v2L2 o2v2L2 o1v1L2 o2v1L1 o2v1L2 o1v1L2 o1v1L2 o1v1L2 o2v1L1 o2v1L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2 o2v2L1 o2v2L1 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o2v0L1 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o2v0L1 o1v1L2 o1v1L2
        tmps_["302_aa_LLvo"]("R,L,a,i")  = -1.00 * l1["bb"]("L,k,c") * r2["abab"]("R,a,c,i,k");
        tmps_["302_aa_LLvo"]("R,L,a,i") -= l1_1["bb"]("L,k,c") * r2_1["abab"]("R,a,c,i,k");
        tmps_["302_aa_LLvo"]("R,L,a,i") += l1["aa"]("L,j,b") * r2["aaaa"]("R,a,b,j,i");
        tmps_["302_aa_LLvo"]("R,L,a,i") += l1_1["aa"]("L,j,b") * (r2_1["aaaa"]("R,a,b,j,i") + r1["aa"]("R,a,j") * t1_1["aa"]("b,i"));
        tmps_["302_aa_LLvo"]("R,L,a,i") += l2_1["bbbb"]("L,l,k,c,d") * r1["bb"]("R,d,l") * t2_1["abab"]("a,c,i,k");
        tmps_["302_aa_LLvo"]("R,L,a,i") += l1["aa"]("L,j,b") * t1["aa"]("b,i") * r1["aa"]("R,a,j");
        tmps_["302_aa_LLvo"]("R,L,a,i") += l1_1["aa"]("L,j,b") * t1["aa"]("b,i") * r1_1["aa"]("R,a,j");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r2_1_aaaa(a,c,k,i) t1_aa(b,j) l1_1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_aa(a,k) t1_aa(b,j) t1_1_aa(c,i) l1_1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r2_aaaa(a,c,k,i) t1_aa(b,j) l1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_bb(d,l) t1_aa(b,j) t2_1_abab(a,c,i,k) l2_1_bbbb(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_aa(a,k) t1_aa(b,j) t1_aa(c,i) l1_aa(k,c)
        //             += -1.00 P(i,j) P(a,b) r2_1_abab(a,c,i,k) t1_aa(b,j) l1_1_bb(k,c)
        //             += -1.00 P(i,j) P(a,b) r2_abab(a,c,i,k) t1_aa(b,j) l1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_1_aa(a,k) t1_aa(b,j) t1_aa(c,i) l1_1_aa(k,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["302_aa_LLvo"]("R,L,a,i") * t1["aa"]("b,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +1.00 r2_1_aaaa(a,c,k,i) t1_bb(b,j) l1_1_aa(k,c)
        //             += +1.00 r1_aa(a,k) t1_bb(b,j) t1_1_aa(c,i) l1_1_aa(k,c)
        //             += +1.00 r2_aaaa(a,c,k,i) t1_bb(b,j) l1_aa(k,c)
        //             += +1.00 r1_bb(d,l) t1_bb(b,j) t2_1_abab(a,c,i,k) l2_1_bbbb(l,k,c,d)
        //             += +1.00 r1_aa(a,k) t1_bb(b,j) t1_aa(c,i) l1_aa(k,c)
        //             += -1.00 r2_1_abab(a,c,i,k) t1_bb(b,j) l1_1_bb(k,c)
        //             += -1.00 r2_abab(a,c,i,k) t1_bb(b,j) l1_bb(k,c)
        //             += +1.00 r1_1_aa(a,k) t1_bb(b,j) t1_aa(c,i) l1_1_aa(k,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["302_aa_LLvo"]("R,L,a,i") * t1["bb"]("b,j");

        // flops: o2v0L1  = o3v2L1 o3v2L1 o2v0L1
        //  mems: o2v0L1  = o2v0L1 o2v0L1 o2v0L1
        tmps_["303_aa_Loo"]("L,i,j")  = t2["aaaa"]("a,b,k,i") * l2["aaaa"]("L,j,k,a,b");
        tmps_["303_aa_Loo"]("L,i,j") += t2_1["aaaa"]("a,b,k,i") * l2_1["aaaa"]("L,j,k,a,b");

        // D_oovv_aaaa += -0.50 P(i,j) r2_aaaa(b,a,l,i) t2_aaaa(d,c,k,j) l2_aaaa(l,k,d,c)
        //             += -0.50 P(i,j) r2_aaaa(b,a,l,i) t2_1_aaaa(d,c,k,j) l2_1_aaaa(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["303_aa_Loo"]("L,j,l") * r2["aaaa"]("R,b,a,l,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += -0.50 r2_abab(a,b,l,j) t2_aaaa(d,c,k,i) l2_aaaa(l,k,d,c)
        //             += -0.50 r2_abab(a,b,l,j) t2_1_aaaa(d,c,k,i) l2_1_aaaa(l,k,d,c)
        // flops: o2v2L2 += o3v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") -= 0.50 * tmps_["303_aa_Loo"]("L,i,l") * r2["abab"]("R,a,b,l,j");

        // flops: o1v1L2  = o1v2L2 o1v2L2 o1v2L2 o1v1L2 o2v2L2 o1v1L2 o2v2L2 o1v1L2 o2v2L2 o1v1L2 o2v1L2 o1v1L2 o2v1L2 o1v1L2 o1v2L2 o1v1L2 o2v2L2 o1v1L2 o2v1L2 o1v1L2 o2v2L2 o1v1L2 o2v2L2 o1v1L2 o2v2L2 o1v1L2 o2v2L2 o1v1L2 o2v1L2 o1v1L2 o1v1L2 o2v2L2 o1v1L2 o2v2L2 o1v1L2 o2v2L2 o1v1L2 o2v1L2 o1v1L2 o1v2L2 o1v1L2 o1v2L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2
        tmps_["304_aa_LLvo"]("L,R,a,i")  = -0.50 * tmps_["140_aa_Lvv"]("L,a,c") * r1_1["aa"]("R,c,i");
        tmps_["304_aa_LLvo"]("L,R,a,i") -= 0.50 * tmps_["300_aa_Lvv"]("L,a,c") * r1["aa"]("R,c,i");
        tmps_["304_aa_LLvo"]("L,R,a,i") += tmps_["64_aa_LLvv"]("L,R,b,a") * t1_1["aa"]("b,i");
        tmps_["304_aa_LLvo"]("L,R,a,i") += tmps_["299_aabb_Lvoov"]("L,a,i,m,e") * r1["bb"]("R,e,m");
        tmps_["304_aa_LLvo"]("L,R,a,i") += tmps_["17_aabb_Lvoov"]("L,a,i,m,e") * r1["bb"]("R,e,m");
        tmps_["304_aa_LLvo"]("L,R,a,i") -= t2_1["abab"]("a,d,i,k") * tmps_["69_bb_LLov"]("L,R,k,d");
        tmps_["304_aa_LLvo"]("L,R,a,i") += 0.50 * t1["aa"]("a,l") * tmps_["135_aa_LLoo"]("L,R,l,i");
        tmps_["304_aa_LLvo"]("L,R,a,i") -= 0.50 * tmps_["303_aa_Loo"]("L,i,j") * r1["aa"]("R,a,j");
        tmps_["304_aa_LLvo"]("L,R,a,i") += tmps_["250_aa_Lvv"]("L,a,c") * r1["aa"]("R,c,i");
        tmps_["304_aa_LLvo"]("L,R,a,i") += tmps_["120_aabb_Lvoov"]("L,a,i,m,e") * r1_1["bb"]("R,e,m");
        tmps_["304_aa_LLvo"]("L,R,a,i") -= 0.50 * tmps_["62_aa_Loo"]("L,i,j") * r1_1["aa"]("R,a,j");
        tmps_["304_aa_LLvo"]("L,R,a,i") += tmps_["16_aabb_Lvoov"]("L,a,i,m,e") * r1_1["bb"]("R,e,m");
        tmps_["304_aa_LLvo"]("L,R,a,i") -= tmps_["119_aaaa_Lvoov"]("L,a,i,j,c") * r1_1["aa"]("R,c,j");
        tmps_["304_aa_LLvo"]("L,R,a,i") += tmps_["18_aabb_Lvoov"]("L,a,i,m,e") * r1["bb"]("R,e,m");
        tmps_["304_aa_LLvo"]("L,R,a,i") -= t2["abab"]("a,d,i,k") * tmps_["70_bb_LLov"]("L,R,k,d");
        tmps_["304_aa_LLvo"]("L,R,a,i") += t1_1["aa"]("a,l") * tmps_["141_aa_LLoo"]("L,R,l,i");
        tmps_["304_aa_LLvo"]("L,R,a,i") -= t2["aaaa"]("a,b,l,i") * tmps_["127_aa_LLov"]("L,R,l,b");
        tmps_["304_aa_LLvo"]("L,R,a,i") -= t2_1["aaaa"]("a,b,l,i") * tmps_["117_aa_LLov"]("L,R,l,b");
        tmps_["304_aa_LLvo"]("L,R,a,i") -= t2["abab"]("a,d,i,k") * tmps_["71_bb_LLov"]("L,R,k,d");
        tmps_["304_aa_LLvo"]("L,R,a,i") += r1_1["aa"]("R,a,j") * tmps_["63_aa_Loo"]("L,i,j");
        tmps_["304_aa_LLvo"]("L,R,a,i") -= 0.50 * tmps_["301_aa_Lvv"]("L,a,c") * r1["aa"]("R,c,i");
        tmps_["304_aa_LLvo"]("L,R,a,i") += tmps_["256_aa_LLvv"]("L,R,b,a") * t1["aa"]("b,i");
        tmps_["303_aa_Loo"].~TArrayD();
        tmps_["301_aa_Lvv"].~TArrayD();
        tmps_["300_aa_Lvv"].~TArrayD();
        tmps_["299_aabb_Lvoov"].~TArrayD();
        tmps_["256_aa_LLvv"].~TArrayD();
        tmps_["127_aa_LLov"].~TArrayD();
        tmps_["120_aabb_Lvoov"].~TArrayD();
        tmps_["119_aaaa_Lvoov"].~TArrayD();
        tmps_["117_aa_LLov"].~TArrayD();
        tmps_["64_aa_LLvv"].~TArrayD();

        // D_oovv_aaaa += +0.50 P(i,j) P(a,b) r2_abab(a,d,l,k) t1_aa(b,j) t1_1_aa(c,i) l2_1_abab(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_abab(a,d,k,l) t1_aa(b,j) t1_1_aa(c,i) l2_1_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_aaaa(a,d,l,k) t1_aa(b,j) t1_1_aa(c,i) l2_1_aaaa(l,k,c,d)
        //             += -0.50 P(i,j) P(a,b) r1_aa(d,i) t2_aaaa(a,c,k,l) t1_aa(b,j) l2_aaaa(k,l,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_bb(d,l) t2_abab(a,c,i,k) t1_aa(b,j) l2_bbbb(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_bb(d,l) t2_aaaa(a,c,k,i) t1_aa(b,j) l2_abab(k,l,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_aa(d,l) t1_aa(b,j) t2_1_abab(a,c,i,k) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r2_aaaa(c,d,l,i) t1_aa(a,k) t1_aa(b,j) l2_aaaa(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_aa(c,i) t1_aa(a,k) t1_aa(b,j) l1_aa(k,c)
        //             += +0.50 P(i,j) P(a,b) r2_1_aaaa(c,d,l,i) t1_aa(a,k) t1_aa(b,j) l2_1_aaaa(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_1_aa(c,i) t1_aa(a,k) t1_aa(b,j) l1_1_aa(k,c)
        //             += +0.50 P(i,j) P(a,b) r2_1_abab(c,d,i,l) t1_aa(a,k) t1_aa(b,j) l2_1_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_1_abab(d,c,i,l) t1_aa(a,k) t1_aa(b,j) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r2_abab(c,d,i,l) t1_aa(a,k) t1_aa(b,j) l2_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_abab(d,c,i,l) t1_aa(a,k) t1_aa(b,j) l2_abab(k,l,d,c)
        //             += -0.50 P(i,j) P(a,b) r1_aa(a,l) t1_aa(b,j) t2_aaaa(d,c,k,i) l2_aaaa(l,k,d,c)
        //             += -0.50 P(i,j) P(a,b) r1_aa(a,l) t1_aa(b,j) t2_1_aaaa(d,c,k,i) l2_1_aaaa(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(d,i) t2_abab(a,c,k,l) t1_aa(b,j) l2_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(d,i) t2_abab(a,c,l,k) t1_aa(b,j) l2_abab(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r1_1_bb(d,l) t2_abab(a,c,i,k) t1_aa(b,j) l2_1_bbbb(l,k,c,d)
        //             += -0.50 P(i,j) P(a,b) r1_1_aa(a,l) t1_aa(b,j) t2_aaaa(d,c,k,i) l2_1_aaaa(l,k,d,c)
        //             += +1.00 P(i,j) P(a,b) r1_1_bb(d,l) t2_aaaa(a,c,k,i) t1_aa(b,j) l2_1_abab(k,l,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_1_aa(d,l) t2_aaaa(a,c,k,i) t1_aa(b,j) l2_1_aaaa(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_bb(d,l) t1_aa(b,j) t2_1_aaaa(a,c,k,i) l2_1_abab(k,l,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_aa(d,l) t2_abab(a,c,i,k) t1_aa(b,j) l2_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r2_abab(c,d,i,l) t1_aa(b,j) t1_1_aa(a,k) l2_1_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_abab(d,c,i,l) t1_aa(b,j) t1_1_aa(a,k) l2_1_abab(k,l,d,c)
        //             += +1.00 P(i,j) P(a,b) r1_aa(c,i) t1_aa(b,j) t1_1_aa(a,k) l1_1_aa(k,c)
        //             += +0.50 P(i,j) P(a,b) r2_aaaa(c,d,l,i) t1_aa(b,j) t1_1_aa(a,k) l2_1_aaaa(l,k,c,d)
        //             += -0.50 P(i,j) P(a,b) r1_1_aa(d,i) t2_aaaa(a,c,k,l) t1_aa(b,j) l2_1_aaaa(k,l,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_aa(d,l) t2_aaaa(a,c,k,i) t1_aa(b,j) l2_aaaa(l,k,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_aa(d,l) t1_aa(b,j) t2_1_aaaa(a,c,k,i) l2_1_aaaa(l,k,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_1_aa(d,l) t2_abab(a,c,i,k) t1_aa(b,j) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_1_aa(a,l) t1_aa(b,j) t2_abab(d,c,i,k) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_1_aa(a,l) t1_aa(b,j) t2_abab(c,d,i,k) l2_1_abab(l,k,c,d)
        //             += -0.50 P(i,j) P(a,b) r1_aa(d,i) t1_aa(b,j) t2_1_aaaa(a,c,k,l) l2_1_aaaa(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_1_abab(a,d,l,k) t1_aa(b,j) t1_aa(c,i) l2_1_abab(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_1_abab(a,d,k,l) t1_aa(b,j) t1_aa(c,i) l2_1_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_aaaa(a,d,l,k) t1_aa(b,j) t1_aa(c,i) l2_aaaa(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_abab(a,d,l,k) t1_aa(b,j) t1_aa(c,i) l2_abab(l,k,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_abab(a,d,k,l) t1_aa(b,j) t1_aa(c,i) l2_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r2_1_aaaa(a,d,l,k) t1_aa(b,j) t1_aa(c,i) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["304_aa_LLvo"]("L,R,a,i") * t1["aa"]("b,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r2_abab(a,d,l,k) t1_bb(b,j) t1_1_aa(c,i) l2_1_abab(l,k,c,d)
        //             += +0.50 r2_abab(a,d,k,l) t1_bb(b,j) t1_1_aa(c,i) l2_1_abab(k,l,c,d)
        //             += +0.50 r2_aaaa(a,d,l,k) t1_bb(b,j) t1_1_aa(c,i) l2_1_aaaa(l,k,c,d)
        //             += -0.50 r1_aa(d,i) t2_aaaa(a,c,k,l) t1_bb(b,j) l2_aaaa(k,l,c,d)
        //             += +1.00 r1_bb(d,l) t2_abab(a,c,i,k) t1_bb(b,j) l2_bbbb(l,k,c,d)
        //             += +1.00 r1_bb(d,l) t2_aaaa(a,c,k,i) t1_bb(b,j) l2_abab(k,l,c,d)
        //             += -1.00 r1_aa(d,l) t1_bb(b,j) t2_1_abab(a,c,i,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r2_aaaa(c,d,l,i) t1_aa(a,k) t1_bb(b,j) l2_aaaa(l,k,c,d)
        //             += +1.00 r1_aa(c,i) t1_aa(a,k) t1_bb(b,j) l1_aa(k,c)
        //             += +0.50 r2_1_aaaa(c,d,l,i) t1_aa(a,k) t1_bb(b,j) l2_1_aaaa(l,k,c,d)
        //             += +1.00 r1_1_aa(c,i) t1_aa(a,k) t1_bb(b,j) l1_1_aa(k,c)
        //             += +0.50 r2_1_abab(c,d,i,l) t1_aa(a,k) t1_bb(b,j) l2_1_abab(k,l,c,d)
        //             += +0.50 r2_1_abab(d,c,i,l) t1_aa(a,k) t1_bb(b,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r2_abab(c,d,i,l) t1_aa(a,k) t1_bb(b,j) l2_abab(k,l,c,d)
        //             += +0.50 r2_abab(d,c,i,l) t1_aa(a,k) t1_bb(b,j) l2_abab(k,l,d,c)
        //             += -0.50 r1_aa(a,l) t1_bb(b,j) t2_aaaa(d,c,k,i) l2_aaaa(l,k,d,c)
        //             += -0.50 r1_aa(a,l) t1_bb(b,j) t2_1_aaaa(d,c,k,i) l2_1_aaaa(l,k,d,c)
        //             += +0.50 r1_aa(d,i) t2_abab(a,c,k,l) t1_bb(b,j) l2_abab(k,l,d,c)
        //             += +0.50 r1_aa(d,i) t2_abab(a,c,l,k) t1_bb(b,j) l2_abab(l,k,d,c)
        //             += +1.00 r1_1_bb(d,l) t2_abab(a,c,i,k) t1_bb(b,j) l2_1_bbbb(l,k,c,d)
        //             += -0.50 r1_1_aa(a,l) t1_bb(b,j) t2_aaaa(d,c,k,i) l2_1_aaaa(l,k,d,c)
        //             += +1.00 r1_1_bb(d,l) t2_aaaa(a,c,k,i) t1_bb(b,j) l2_1_abab(k,l,c,d)
        //             += -1.00 r1_1_aa(d,l) t2_aaaa(a,c,k,i) t1_bb(b,j) l2_1_aaaa(l,k,c,d)
        //             += +1.00 r1_bb(d,l) t1_bb(b,j) t2_1_aaaa(a,c,k,i) l2_1_abab(k,l,c,d)
        //             += -1.00 r1_aa(d,l) t2_abab(a,c,i,k) t1_bb(b,j) l2_abab(l,k,d,c)
        //             += +0.50 r2_abab(c,d,i,l) t1_bb(b,j) t1_1_aa(a,k) l2_1_abab(k,l,c,d)
        //             += +0.50 r2_abab(d,c,i,l) t1_bb(b,j) t1_1_aa(a,k) l2_1_abab(k,l,d,c)
        //             += +1.00 r1_aa(c,i) t1_bb(b,j) t1_1_aa(a,k) l1_1_aa(k,c)
        //             += +0.50 r2_aaaa(c,d,l,i) t1_bb(b,j) t1_1_aa(a,k) l2_1_aaaa(l,k,c,d)
        //             += -0.50 r1_1_aa(d,i) t2_aaaa(a,c,k,l) t1_bb(b,j) l2_1_aaaa(k,l,c,d)
        //             += -1.00 r1_aa(d,l) t2_aaaa(a,c,k,i) t1_bb(b,j) l2_aaaa(l,k,c,d)
        //             += -1.00 r1_aa(d,l) t1_bb(b,j) t2_1_aaaa(a,c,k,i) l2_1_aaaa(l,k,c,d)
        //             += -1.00 r1_1_aa(d,l) t2_abab(a,c,i,k) t1_bb(b,j) l2_1_abab(l,k,d,c)
        //             += +0.50 r1_1_aa(a,l) t1_bb(b,j) t2_abab(d,c,i,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r1_1_aa(a,l) t1_bb(b,j) t2_abab(c,d,i,k) l2_1_abab(l,k,c,d)
        //             += -0.50 r1_aa(d,i) t1_bb(b,j) t2_1_aaaa(a,c,k,l) l2_1_aaaa(k,l,c,d)
        //             += +0.50 r2_1_abab(a,d,l,k) t1_bb(b,j) t1_aa(c,i) l2_1_abab(l,k,c,d)
        //             += +0.50 r2_1_abab(a,d,k,l) t1_bb(b,j) t1_aa(c,i) l2_1_abab(k,l,c,d)
        //             += +0.50 r2_aaaa(a,d,l,k) t1_bb(b,j) t1_aa(c,i) l2_aaaa(l,k,c,d)
        //             += +0.50 r2_abab(a,d,l,k) t1_bb(b,j) t1_aa(c,i) l2_abab(l,k,c,d)
        //             += +0.50 r2_abab(a,d,k,l) t1_bb(b,j) t1_aa(c,i) l2_abab(k,l,c,d)
        //             += +0.50 r2_1_aaaa(a,d,l,k) t1_bb(b,j) t1_aa(c,i) l2_1_aaaa(l,k,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["304_aa_LLvo"]("L,R,a,i") * t1["bb"]("b,j");

        // flops: o1v1L1  = o2v2L1 o2v1L1 o2v1L1 o2v2L1 o1v1L1 o1v1L1 o2v1L1 o1v2L1 o1v1L1 o2v1L1 o1v1L1 o1v1L1 o1v2L1 o1v1L1
        //  mems: o1v1L1  = o1v1L1 o2v0L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
        tmps_["305_aa_Lov"]("L,i,a")  = -1.00 * l1_1["bb"]("L,l,d") * t2["abab"]("a,d,i,l");
        tmps_["305_aa_Lov"]("L,i,a") += l1_1["aa"]("L,k,c") * t1["aa"]("c,i") * t1["aa"]("a,k");
        tmps_["305_aa_Lov"]("L,i,a") += l1_1["aa"]("L,k,c") * t2["aaaa"]("a,c,k,i");
        tmps_["305_aa_Lov"]("L,i,a") += 0.50 * t1["aa"]("a,j") * tmps_["131_aa_Loo"]("L,i,j");
        tmps_["305_aa_Lov"]("L,i,a") += t1["aa"]("b,i") * tmps_["29_aa_Lvv"]("L,b,a");
        tmps_["305_aa_Lov"]("L,i,a") += t1["aa"]("a,j") * tmps_["63_aa_Loo"]("L,i,j");
        tmps_["305_aa_Lov"]("L,i,a") += 0.50 * t1["aa"]("b,i") * tmps_["246_aa_Lvv"]("L,b,a");
        tmps_["246_aa_Lvv"].~TArrayD();

        // D_oovv_abab += +0.50 r0 t2_abab(a,c,k,l) t1_aa(d,i) t1_1_bb(b,j) l2_1_abab(k,l,d,c)
        //             += +0.50 r0 t2_abab(a,c,l,k) t1_aa(d,i) t1_1_bb(b,j) l2_1_abab(l,k,d,c)
        //             += +0.50 r0 t1_aa(a,l) t2_aaaa(d,c,k,i) t1_1_bb(b,j) l2_1_aaaa(k,l,d,c)
        //             += +0.50 r0 t1_aa(a,l) t2_abab(d,c,i,k) t1_1_bb(b,j) l2_1_abab(l,k,d,c)
        //             += +0.50 r0 t1_aa(a,l) t2_abab(c,d,i,k) t1_1_bb(b,j) l2_1_abab(l,k,c,d)
        //             += +1.00 r0 t1_aa(a,k) t1_aa(c,i) t1_1_bb(b,j) l1_1_aa(k,c)
        //             += +1.00 r0 t2_aaaa(a,c,k,i) t1_1_bb(b,j) l1_1_aa(k,c)
        //             += -1.00 r0 t2_abab(a,c,i,k) t1_1_bb(b,j) l1_1_bb(k,c)
        //             += +0.50 r0 t2_aaaa(a,c,k,l) t1_aa(d,i) t1_1_bb(b,j) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["305_aa_Lov"]("L,i,a") * tmps_["226_bb_Lvo"]("R,b,j");

        // D_oovv_aaaa += +0.50 P(i,j) P(a,b) r0_1 t2_abab(a,c,k,l) t1_aa(b,j) t1_aa(d,i) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0_1 t2_abab(a,c,l,k) t1_aa(b,j) t1_aa(d,i) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r0_1 t1_aa(a,l) t1_aa(b,j) t2_aaaa(d,c,k,i) l2_1_aaaa(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0_1 t1_aa(a,l) t1_aa(b,j) t2_abab(d,c,i,k) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r0_1 t1_aa(a,l) t1_aa(b,j) t2_abab(c,d,i,k) l2_1_abab(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r0_1 t1_aa(a,k) t1_aa(b,j) t1_aa(c,i) l1_1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r0_1 t2_aaaa(a,c,k,i) t1_aa(b,j) l1_1_aa(k,c)
        //             += -1.00 P(i,j) P(a,b) r0_1 t2_abab(a,c,i,k) t1_aa(b,j) l1_1_bb(k,c)
        //             += +0.50 P(i,j) P(a,b) r0_1 t2_aaaa(a,c,k,l) t1_aa(b,j) t1_aa(d,i) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["305_aa_Lov"]("L,i,a") * tmps_["254_aa_Lvo"]("R,b,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_aaaa += +0.50 P(i,j) P(a,b) r0 t2_abab(b,c,k,l) t1_aa(d,j) t1_1_aa(a,i) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t2_abab(b,c,l,k) t1_aa(d,j) t1_1_aa(a,i) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_aa(b,l) t2_aaaa(d,c,k,j) t1_1_aa(a,i) l2_1_aaaa(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_aa(b,l) t2_abab(d,c,j,k) t1_1_aa(a,i) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_aa(b,l) t2_abab(c,d,j,k) t1_1_aa(a,i) l2_1_abab(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r0 t1_aa(b,k) t1_aa(c,j) t1_1_aa(a,i) l1_1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r0 t2_aaaa(b,c,k,j) t1_1_aa(a,i) l1_1_aa(k,c)
        //             += -1.00 P(i,j) P(a,b) r0 t2_abab(b,c,j,k) t1_1_aa(a,i) l1_1_bb(k,c)
        //             += +0.50 P(i,j) P(a,b) r0 t2_aaaa(b,c,k,l) t1_aa(d,j) t1_1_aa(a,i) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["305_aa_Lov"]("L,j,b") * tmps_["255_aa_Lvo"]("R,a,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r0_1 t2_abab(a,c,k,l) t1_bb(b,j) t1_aa(d,i) l2_1_abab(k,l,d,c)
        //             += +0.50 r0_1 t2_abab(a,c,l,k) t1_bb(b,j) t1_aa(d,i) l2_1_abab(l,k,d,c)
        //             += +0.50 r0_1 t1_aa(a,l) t1_bb(b,j) t2_aaaa(d,c,k,i) l2_1_aaaa(k,l,d,c)
        //             += +0.50 r0_1 t1_aa(a,l) t1_bb(b,j) t2_abab(d,c,i,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r0_1 t1_aa(a,l) t1_bb(b,j) t2_abab(c,d,i,k) l2_1_abab(l,k,c,d)
        //             += +1.00 r0_1 t1_aa(a,k) t1_bb(b,j) t1_aa(c,i) l1_1_aa(k,c)
        //             += +1.00 r0_1 t2_aaaa(a,c,k,i) t1_bb(b,j) l1_1_aa(k,c)
        //             += -1.00 r0_1 t2_abab(a,c,i,k) t1_bb(b,j) l1_1_bb(k,c)
        //             += +0.50 r0_1 t2_aaaa(a,c,k,l) t1_bb(b,j) t1_aa(d,i) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["305_aa_Lov"]("L,i,a") * tmps_["225_bb_Lvo"]("R,b,j");

        // D_oovv_aaaa += +0.50 P(i,j) P(a,b) r1_1_aa(a,i) t2_abab(b,c,k,l) t1_aa(d,j) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_1_aa(a,i) t2_abab(b,c,l,k) t1_aa(d,j) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_1_aa(a,i) t1_aa(b,l) t2_aaaa(d,c,k,j) l2_1_aaaa(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_1_aa(a,i) t1_aa(b,l) t2_abab(d,c,j,k) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_1_aa(a,i) t1_aa(b,l) t2_abab(c,d,j,k) l2_1_abab(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_1_aa(a,i) t1_aa(b,k) t1_aa(c,j) l1_1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_1_aa(a,i) t2_aaaa(b,c,k,j) l1_1_aa(k,c)
        //             += -1.00 P(i,j) P(a,b) r1_1_aa(a,i) t2_abab(b,c,j,k) l1_1_bb(k,c)
        //             += +0.50 P(i,j) P(a,b) r1_1_aa(a,i) t2_aaaa(b,c,k,l) t1_aa(d,j) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["305_aa_Lov"]("L,j,b") * r1_1["aa"]("R,a,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r1_1_bb(b,j) t2_abab(a,c,k,l) t1_aa(d,i) l2_1_abab(k,l,d,c)
        //             += +0.50 r1_1_bb(b,j) t2_abab(a,c,l,k) t1_aa(d,i) l2_1_abab(l,k,d,c)
        //             += +0.50 r1_1_bb(b,j) t1_aa(a,l) t2_aaaa(d,c,k,i) l2_1_aaaa(k,l,d,c)
        //             += +0.50 r1_1_bb(b,j) t1_aa(a,l) t2_abab(d,c,i,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r1_1_bb(b,j) t1_aa(a,l) t2_abab(c,d,i,k) l2_1_abab(l,k,c,d)
        //             += +1.00 r1_1_bb(b,j) t1_aa(a,k) t1_aa(c,i) l1_1_aa(k,c)
        //             += +1.00 r1_1_bb(b,j) t2_aaaa(a,c,k,i) l1_1_aa(k,c)
        //             += -1.00 r1_1_bb(b,j) t2_abab(a,c,i,k) l1_1_bb(k,c)
        //             += +0.50 r1_1_bb(b,j) t2_aaaa(a,c,k,l) t1_aa(d,i) l2_1_aaaa(k,l,d,c)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["305_aa_Lov"]("L,i,a") * r1_1["bb"]("R,b,j");

        // flops: o1v1L1  = o2v1L1 o1v2L1 o1v2L1 o1v1L1 o1v2L1 o1v1L1 o2v1L1 o2v1L1 o2v2L1 o1v1L1 o2v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o2v1L1 o2v1L1 o1v1L1 o2v2L1 o1v1L1 o2v2L1 o1v1L1 o1v1L1 o1v1L1 o1v2L1 o1v1L1 o1v2L1 o1v1L1 o1v2L1 o1v1L1 o2v1L1 o1v1L1
        //  mems: o1v1L1  = o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o2v0L1 o1v1L1 o1v1L1 o1v1L1 o2v0L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o2v0L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1 o1v1L1
        tmps_["306_aa_Lov"]("L,i,a")  = -1.00 * t1_1["aa"]("a,j") * tmps_["62_aa_Loo"]("L,i,j");
        tmps_["306_aa_Lov"]("L,i,a") -= t1_1["aa"]("d,i") * tmps_["140_aa_Lvv"]("L,a,d");
        tmps_["306_aa_Lov"]("L,i,a") += t1["aa"]("b,i") * tmps_["247_aa_Lvv"]("L,b,a");
        tmps_["306_aa_Lov"]("L,i,a") += 2.00 * t1_1["aa"]("d,i") * tmps_["29_aa_Lvv"]("L,d,a");
        tmps_["306_aa_Lov"]("L,i,a") += 2.00 * l1["aa"]("L,j,d") * t1["aa"]("d,i") * t1["aa"]("a,j");
        tmps_["306_aa_Lov"]("L,i,a") -= 2.00 * l1["bb"]("L,m,c") * t2["abab"]("a,c,i,m");
        tmps_["306_aa_Lov"]("L,i,a") += 2.00 * l1_1["aa"]("L,j,d") * t1_1["aa"]("d,i") * t1["aa"]("a,j");
        tmps_["306_aa_Lov"]("L,i,a") += 2.00 * l1["aa"]("L,j,d") * t2["aaaa"]("a,d,j,i");
        tmps_["306_aa_Lov"]("L,i,a") += 2.00 * l1_1["aa"]("L,j,d") * t1["aa"]("d,i") * t1_1["aa"]("a,j");
        tmps_["306_aa_Lov"]("L,i,a") -= 2.00 * l1_1["bb"]("L,m,c") * t2_1["abab"]("a,c,i,m");
        tmps_["306_aa_Lov"]("L,i,a") += 2.00 * l1_1["aa"]("L,j,d") * t2_1["aaaa"]("a,d,j,i");
        tmps_["306_aa_Lov"]("L,i,a") += 2.00 * t1["aa"]("b,i") * tmps_["250_aa_Lvv"]("L,a,b");
        tmps_["306_aa_Lov"]("L,i,a") += t1["aa"]("b,i") * tmps_["248_aa_Lvv"]("L,a,b");
        tmps_["306_aa_Lov"]("L,i,a") += 2.00 * t1["aa"]("b,i") * tmps_["249_aa_Lvv"]("L,b,a");
        tmps_["306_aa_Lov"]("L,i,a") += 2.00 * t1_1["aa"]("a,j") * tmps_["63_aa_Loo"]("L,i,j");
        tmps_["250_aa_Lvv"].~TArrayD();
        tmps_["248_aa_Lvv"].~TArrayD();
        tmps_["247_aa_Lvv"].~TArrayD();
        tmps_["140_aa_Lvv"].~TArrayD();
        tmps_["62_aa_Loo"].~TArrayD();

        // D_oovv_aaaa += +0.50 P(i,j) P(a,b) r1_aa(a,i) t1_aa(d,j) t2_1_aaaa(b,c,k,l) l2_1_aaaa(k,l,d,c)
        //             += -0.50 P(i,j) P(a,b) r1_aa(a,i) t2_aaaa(b,d,k,l) t1_1_aa(c,j) l2_1_aaaa(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,i) t2_abab(b,d,k,l) t1_1_aa(c,j) l2_1_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,i) t2_abab(b,d,l,k) t1_1_aa(c,j) l2_1_abab(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_aa(a,i) t1_aa(b,k) t1_aa(c,j) l1_aa(k,c)
        //             += -1.00 P(i,j) P(a,b) r1_aa(a,i) t2_abab(b,c,j,k) l1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_aa(a,i) t1_aa(b,k) t1_1_aa(c,j) l1_1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_aa(a,i) t2_aaaa(b,c,k,j) l1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_aa(a,i) t1_aa(c,j) t1_1_aa(b,k) l1_1_aa(k,c)
        //             += -1.00 P(i,j) P(a,b) r1_aa(a,i) t2_1_abab(b,c,j,k) l1_1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r1_aa(a,i) t2_1_aaaa(b,c,k,j) l1_1_aa(k,c)
        //             += -0.50 P(i,j) P(a,b) r1_aa(a,i) t2_aaaa(d,c,l,j) t1_1_aa(b,k) l2_1_aaaa(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,i) t2_abab(b,c,k,l) t1_aa(d,j) l2_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,i) t2_abab(b,c,l,k) t1_aa(d,j) l2_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,i) t2_aaaa(b,c,k,l) t1_aa(d,j) l2_aaaa(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,i) t1_aa(d,j) t2_1_abab(b,c,k,l) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,i) t1_aa(d,j) t2_1_abab(b,c,l,k) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,i) t2_abab(d,c,j,l) t1_1_aa(b,k) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r1_aa(a,i) t2_abab(c,d,j,l) t1_1_aa(b,k) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["306_aa_Lov"]("L,j,b") * r1["aa"]("R,a,i");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +0.50 r0 t1_bb(b,j) t1_aa(d,i) t2_1_aaaa(a,c,k,l) l2_1_aaaa(k,l,d,c)
        //             += -0.50 r0 t2_aaaa(a,d,k,l) t1_bb(b,j) t1_1_aa(c,i) l2_1_aaaa(k,l,d,c)
        //             += +0.50 r0 t2_abab(a,d,k,l) t1_bb(b,j) t1_1_aa(c,i) l2_1_abab(k,l,c,d)
        //             += +0.50 r0 t2_abab(a,d,l,k) t1_bb(b,j) t1_1_aa(c,i) l2_1_abab(l,k,c,d)
        //             += +1.00 r0 t1_aa(a,k) t1_bb(b,j) t1_aa(c,i) l1_aa(k,c)
        //             += -1.00 r0 t2_abab(a,c,i,k) t1_bb(b,j) l1_bb(k,c)
        //             += +1.00 r0 t1_aa(a,k) t1_bb(b,j) t1_1_aa(c,i) l1_1_aa(k,c)
        //             += +1.00 r0 t2_aaaa(a,c,k,i) t1_bb(b,j) l1_aa(k,c)
        //             += +1.00 r0 t1_bb(b,j) t1_aa(c,i) t1_1_aa(a,k) l1_1_aa(k,c)
        //             += -1.00 r0 t1_bb(b,j) t2_1_abab(a,c,i,k) l1_1_bb(k,c)
        //             += +1.00 r0 t1_bb(b,j) t2_1_aaaa(a,c,k,i) l1_1_aa(k,c)
        //             += -0.50 r0 t1_bb(b,j) t2_aaaa(d,c,l,i) t1_1_aa(a,k) l2_1_aaaa(k,l,d,c)
        //             += +0.50 r0 t2_abab(a,c,k,l) t1_bb(b,j) t1_aa(d,i) l2_abab(k,l,d,c)
        //             += +0.50 r0 t2_abab(a,c,l,k) t1_bb(b,j) t1_aa(d,i) l2_abab(l,k,d,c)
        //             += +0.50 r0 t2_aaaa(a,c,k,l) t1_bb(b,j) t1_aa(d,i) l2_aaaa(k,l,d,c)
        //             += +0.50 r0 t1_bb(b,j) t1_aa(d,i) t2_1_abab(a,c,k,l) l2_1_abab(k,l,d,c)
        //             += +0.50 r0 t1_bb(b,j) t1_aa(d,i) t2_1_abab(a,c,l,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r0 t1_bb(b,j) t2_abab(d,c,i,l) t1_1_aa(a,k) l2_1_abab(k,l,d,c)
        //             += +0.50 r0 t1_bb(b,j) t2_abab(c,d,i,l) t1_1_aa(a,k) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * tmps_["306_aa_Lov"]("L,i,a") * tmps_["224_bb_Lvo"]("R,b,j");

        // D_oovv_abab += +0.50 r1_bb(b,j) t1_aa(d,i) t2_1_aaaa(a,c,k,l) l2_1_aaaa(k,l,d,c)
        //             += -0.50 r1_bb(b,j) t2_aaaa(a,d,k,l) t1_1_aa(c,i) l2_1_aaaa(k,l,d,c)
        //             += +0.50 r1_bb(b,j) t2_abab(a,d,k,l) t1_1_aa(c,i) l2_1_abab(k,l,c,d)
        //             += +0.50 r1_bb(b,j) t2_abab(a,d,l,k) t1_1_aa(c,i) l2_1_abab(l,k,c,d)
        //             += +1.00 r1_bb(b,j) t1_aa(a,k) t1_aa(c,i) l1_aa(k,c)
        //             += -1.00 r1_bb(b,j) t2_abab(a,c,i,k) l1_bb(k,c)
        //             += +1.00 r1_bb(b,j) t1_aa(a,k) t1_1_aa(c,i) l1_1_aa(k,c)
        //             += +1.00 r1_bb(b,j) t2_aaaa(a,c,k,i) l1_aa(k,c)
        //             += +1.00 r1_bb(b,j) t1_aa(c,i) t1_1_aa(a,k) l1_1_aa(k,c)
        //             += -1.00 r1_bb(b,j) t2_1_abab(a,c,i,k) l1_1_bb(k,c)
        //             += +1.00 r1_bb(b,j) t2_1_aaaa(a,c,k,i) l1_1_aa(k,c)
        //             += -0.50 r1_bb(b,j) t2_aaaa(d,c,l,i) t1_1_aa(a,k) l2_1_aaaa(k,l,d,c)
        //             += +0.50 r1_bb(b,j) t2_abab(a,c,k,l) t1_aa(d,i) l2_abab(k,l,d,c)
        //             += +0.50 r1_bb(b,j) t2_abab(a,c,l,k) t1_aa(d,i) l2_abab(l,k,d,c)
        //             += +0.50 r1_bb(b,j) t2_aaaa(a,c,k,l) t1_aa(d,i) l2_aaaa(k,l,d,c)
        //             += +0.50 r1_bb(b,j) t1_aa(d,i) t2_1_abab(a,c,k,l) l2_1_abab(k,l,d,c)
        //             += +0.50 r1_bb(b,j) t1_aa(d,i) t2_1_abab(a,c,l,k) l2_1_abab(l,k,d,c)
        //             += +0.50 r1_bb(b,j) t2_abab(d,c,i,l) t1_1_aa(a,k) l2_1_abab(k,l,d,c)
        //             += +0.50 r1_bb(b,j) t2_abab(c,d,i,l) t1_1_aa(a,k) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += 0.50 * tmps_["306_aa_Lov"]("L,i,a") * r1["bb"]("R,b,j");

        // D_oovv_aaaa += +0.50 P(i,j) P(a,b) r0 t1_aa(b,j) t1_aa(d,i) t2_1_aaaa(a,c,k,l) l2_1_aaaa(k,l,d,c)
        //             += -0.50 P(i,j) P(a,b) r0 t2_aaaa(a,d,k,l) t1_aa(b,j) t1_1_aa(c,i) l2_1_aaaa(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t2_abab(a,d,k,l) t1_aa(b,j) t1_1_aa(c,i) l2_1_abab(k,l,c,d)
        //             += +0.50 P(i,j) P(a,b) r0 t2_abab(a,d,l,k) t1_aa(b,j) t1_1_aa(c,i) l2_1_abab(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r0 t1_aa(a,k) t1_aa(b,j) t1_aa(c,i) l1_aa(k,c)
        //             += -1.00 P(i,j) P(a,b) r0 t2_abab(a,c,i,k) t1_aa(b,j) l1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r0 t1_aa(a,k) t1_aa(b,j) t1_1_aa(c,i) l1_1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r0 t2_aaaa(a,c,k,i) t1_aa(b,j) l1_aa(k,c)
        //             += +1.00 P(i,j) P(a,b) r0 t1_aa(b,j) t1_aa(c,i) t1_1_aa(a,k) l1_1_aa(k,c)
        //             += -1.00 P(i,j) P(a,b) r0 t1_aa(b,j) t2_1_abab(a,c,i,k) l1_1_bb(k,c)
        //             += +1.00 P(i,j) P(a,b) r0 t1_aa(b,j) t2_1_aaaa(a,c,k,i) l1_1_aa(k,c)
        //             += -0.50 P(i,j) P(a,b) r0 t1_aa(b,j) t2_aaaa(d,c,l,i) t1_1_aa(a,k) l2_1_aaaa(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t2_abab(a,c,k,l) t1_aa(b,j) t1_aa(d,i) l2_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t2_abab(a,c,l,k) t1_aa(b,j) t1_aa(d,i) l2_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t2_aaaa(a,c,k,l) t1_aa(b,j) t1_aa(d,i) l2_aaaa(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_aa(b,j) t1_aa(d,i) t2_1_abab(a,c,k,l) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_aa(b,j) t1_aa(d,i) t2_1_abab(a,c,l,k) l2_1_abab(l,k,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_aa(b,j) t2_abab(d,c,i,l) t1_1_aa(a,k) l2_1_abab(k,l,d,c)
        //             += +0.50 P(i,j) P(a,b) r0 t1_aa(b,j) t2_abab(c,d,i,l) t1_1_aa(a,k) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = 0.50 * tmps_["306_aa_Lov"]("L,i,a") * tmps_["253_aa_Lvo"]("R,b,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // flops: o1v1L2  = o2v1L2 o2v1L2 o2v1L2 o1v1L2 o1v1L2
        //  mems: o1v1L2  = o1v1L2 o1v1L2 o1v1L2 o1v1L2 o1v1L2
        tmps_["307_aa_LLvo"]("L,R,a,i")  = -1.00 * t1["aa"]("a,j") * tmps_["138_aa_LLoo"]("L,R,j,i");
        tmps_["307_aa_LLvo"]("L,R,a,i") += t1["aa"]("a,j") * tmps_["137_aa_LLoo"]("L,R,j,i");
        tmps_["307_aa_LLvo"]("L,R,a,i") += t1_1["aa"]("a,j") * tmps_["121_aa_LLoo"]("L,R,j,i");

        // D_oovv_aaaa += +1.00 P(i,j) P(a,b) r1_bb(d,l) t1_aa(b,j) t1_aa(c,i) t1_1_aa(a,k) l2_1_abab(k,l,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_aa(d,l) t1_aa(b,j) t1_aa(c,i) t1_1_aa(a,k) l2_1_aaaa(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_bb(d,l) t1_aa(a,k) t1_aa(b,j) t1_aa(c,i) l2_abab(k,l,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_1_bb(d,l) t1_aa(a,k) t1_aa(b,j) t1_aa(c,i) l2_1_abab(k,l,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_aa(d,l) t1_aa(a,k) t1_aa(b,j) t1_aa(c,i) l2_aaaa(l,k,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_1_aa(d,l) t1_aa(a,k) t1_aa(b,j) t1_aa(c,i) l2_1_aaaa(l,k,c,d)
        //             += -1.00 P(i,j) P(a,b) r1_aa(d,l) t1_aa(a,k) t1_aa(b,j) t1_1_aa(c,i) l2_1_aaaa(l,k,c,d)
        //             += +1.00 P(i,j) P(a,b) r1_bb(d,l) t1_aa(a,k) t1_aa(b,j) t1_1_aa(c,i) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j")  = tmps_["307_aa_LLvo"]("L,R,a,i") * t1["aa"]("b,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,a,b,j,i");
        D_oovv_aaaa("L,R,i,j,a,b") -= tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,i,j");
        D_oovv_aaaa("L,R,i,j,a,b") += tmps_["perm_aaaa_LLvvoo"]("L,R,b,a,j,i");
        tmps_["perm_aaaa_LLvvoo"].~TArrayD();

        // D_oovv_abab += +1.00 r1_bb(d,l) t1_bb(b,j) t1_aa(c,i) t1_1_aa(a,k) l2_1_abab(k,l,c,d)
        //             += -1.00 r1_aa(d,l) t1_bb(b,j) t1_aa(c,i) t1_1_aa(a,k) l2_1_aaaa(l,k,c,d)
        //             += +1.00 r1_bb(d,l) t1_aa(a,k) t1_bb(b,j) t1_aa(c,i) l2_abab(k,l,c,d)
        //             += +1.00 r1_1_bb(d,l) t1_aa(a,k) t1_bb(b,j) t1_aa(c,i) l2_1_abab(k,l,c,d)
        //             += -1.00 r1_aa(d,l) t1_aa(a,k) t1_bb(b,j) t1_aa(c,i) l2_aaaa(l,k,c,d)
        //             += -1.00 r1_1_aa(d,l) t1_aa(a,k) t1_bb(b,j) t1_aa(c,i) l2_1_aaaa(l,k,c,d)
        //             += -1.00 r1_aa(d,l) t1_aa(a,k) t1_bb(b,j) t1_1_aa(c,i) l2_1_aaaa(l,k,c,d)
        //             += +1.00 r1_bb(d,l) t1_aa(a,k) t1_bb(b,j) t1_1_aa(c,i) l2_1_abab(k,l,c,d)
        // flops: o2v2L2 += o2v2L2
        //  mems: o2v2L2 += o2v2L2
        D_oovv_abab("L,R,i,j,a,b") += tmps_["307_aa_LLvo"]("L,R,a,i") * t1["bb"]("b,j");
    }
} // hilbert
#endif