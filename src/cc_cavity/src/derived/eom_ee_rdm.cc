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

#include <psi4/libmints/writer.h>
#include "cc_cavity/include/derived/eom_ee_rdm.h"
#include "misc/nonsym_davidson_solver.h"

namespace hilbert {
    EOM_EE_RDM::EOM_EE_RDM(const shared_ptr<EOM_Driver>& eom_driver, Options & options) : EOM_RDM(eom_driver, options) {
    }

    void EOM_EE_RDM::compute_eom_1rdm() {

        // reinitialize the 1-RDMs if not already done
        vector<string> dims = {"oa", "va", "ob", "vb"};
        vector<size_t> dim_vec = {oa_, va_, ob_, vb_};
        for (int i = 0; i < dims.size(); i++) {
            for (int j = 0; j < dims.size(); j++) {
                string dim1 = dims[i], dim2 = dims[j];

                string ov = dim1.substr(0, 1) + dim2.substr(0, 1);
                string spin = dim1.substr(1, 1) + dim2.substr(1, 1);
                string rdm_str = "D1_" + spin + "_" + ov;
                RDM_blks_[rdm_str] = TArrayD(world_, makeRange({M_, M_, dim_vec[i], dim_vec[j]}));
                RDM_blks_[rdm_str].fill(0.0);
            }
        }
        world_.gop.fence();

        TArrayD &rdm_vv_aa = RDM_blks_["D1_aa_vv"];
        TArrayD &rdm_vv_bb = RDM_blks_["D1_bb_vv"];
        TArrayD &rdm_vo_aa = RDM_blks_["D1_aa_vo"];
        TArrayD &rdm_vo_bb = RDM_blks_["D1_bb_vo"];
        TArrayD &rdm_ov_aa = RDM_blks_["D1_aa_ov"];
        TArrayD &rdm_ov_bb = RDM_blks_["D1_bb_ov"];
        TArrayD &rdm_oo_aa = RDM_blks_["D1_aa_oo"];
        TArrayD &rdm_oo_bb = RDM_blks_["D1_bb_oo"];

        Timer rdm_timer; rdm_timer.start();

        // get cavity information
        double w0 = eom_driver_->cc_wfn_->cavity_frequency_[2];
        double coupling_factor_z = w0 * eom_driver_->cc_wfn_->cavity_coupling_strength_[2];

        // extract 0-body amplitudes
        TArrayD r0 = eom_driver_->evec_blks_["r0"];
        TArrayD l0 = eom_driver_->evec_blks_["l0"];

        double u0 = eom_driver_->cc_wfn_->scalar_amps_["u0"];
        TArrayD s0 = eom_driver_->evec_blks_["m0"];
        TArrayD m0 = eom_driver_->evec_blks_["s0"];

        // extract 1-body amplitudes
        std::map<std::string, TA::TArrayD> t1 {
                {"aa_vo", eom_driver_->cc_wfn_->amplitudes_["t1_aa"]},
                {"bb_vo", eom_driver_->cc_wfn_->amplitudes_["t1_bb"]}
        };
        std::map<std::string, TA::TArrayD> r1 {
                {"aa_Lvo", eom_driver_->evec_blks_["r1_aa"]},
                {"bb_Lvo", eom_driver_->evec_blks_["r1_bb"]}
        };
        std::map<std::string, TA::TArrayD> l1 {
                {"aa_Lov", eom_driver_->evec_blks_["l1_aa"]},
                {"bb_Lov", eom_driver_->evec_blks_["l1_bb"]}
        };


        std::map<std::string, TA::TArrayD> u1 {
                {"aa_vo", eom_driver_->cc_wfn_->amplitudes_["u1_aa"]},
                {"bb_vo", eom_driver_->cc_wfn_->amplitudes_["u1_bb"]}
        };
        std::map<std::string, TA::TArrayD> s1 {
                {"aa_Lvo", eom_driver_->evec_blks_["s1_aa"]},
                {"bb_Lvo", eom_driver_->evec_blks_["s1_bb"]}
        };
        std::map<std::string, TA::TArrayD> m1 {
                {"aa_Lov", eom_driver_->evec_blks_["m1_aa"]},
                {"bb_Lov", eom_driver_->evec_blks_["m1_bb"]}
        };


        // extract 2-body amplitudes
        std::map<std::string, TA::TArrayD> t2 {
                {"aaaa_vvoo", eom_driver_->cc_wfn_->amplitudes_["t2_aaaa"]},
                {"abab_vvoo", eom_driver_->cc_wfn_->amplitudes_["t2_abab"]},
                {"bbbb_vvoo", eom_driver_->cc_wfn_->amplitudes_["t2_bbbb"]}
        };
        std::map<std::string, TA::TArrayD> r2 {
                {"aaaa_Lvvoo", eom_driver_->evec_blks_["r2_aaaa"]},
                {"abab_Lvvoo", eom_driver_->evec_blks_["r2_abab"]},
                {"bbbb_Lvvoo", eom_driver_->evec_blks_["r2_bbbb"]}
        };
        std::map<std::string, TA::TArrayD> l2 {
                {"aaaa_Loovv", eom_driver_->evec_blks_["l2_aaaa"]},
                {"abab_Loovv", eom_driver_->evec_blks_["l2_abab"]},
                {"bbbb_Loovv", eom_driver_->evec_blks_["l2_bbbb"]}
        };

        std::map<std::string, TA::TArrayD> u2 {
                {"aaaa_vvoo", eom_driver_->cc_wfn_->amplitudes_["u2_aaaa"]},
                {"abab_vvoo", eom_driver_->cc_wfn_->amplitudes_["u2_abab"]},
                {"bbbb_vvoo", eom_driver_->cc_wfn_->amplitudes_["u2_bbbb"]}
        };
        std::map<std::string, TA::TArrayD> s2 {
                {"aaaa_Lvvoo", eom_driver_->evec_blks_["s2_aaaa"]},
                {"abab_Lvvoo", eom_driver_->evec_blks_["s2_abab"]},
                {"bbbb_Lvvoo", eom_driver_->evec_blks_["s2_bbbb"]}
        };
        std::map<std::string, TA::TArrayD> m2 {
                {"aaaa_Loovv", eom_driver_->evec_blks_["m2_aaaa"]},
                {"abab_Loovv", eom_driver_->evec_blks_["m2_abab"]},
                {"bbbb_Loovv", eom_driver_->evec_blks_["m2_bbbb"]}
        };

        // get integrals
        TArrayMap &Id  = eom_driver_->cc_wfn_->Id_blks_;
        TArrayMap &eri = eom_driver_->cc_wfn_->V_blks_;
        TArrayMap &f   = eom_driver_->cc_wfn_->F_blks_;

        rdm_timer.start();
        {
            map<string, double> scalars_;
            TArrayMap tmps_;

            scalars_["0"]  = 1.00 * dot(l0("X"), r0("X"));
            scalars_["1"]  = 1.00 * dot(l1["aa_Lov"]("X,k,a"), r1["aa_Lvo"]("X,a,k"));
            scalars_["2"]  = 1.00 * dot(l2["aaaa_Loovv"]("X,l,k,a,b"), r2["aaaa_Lvvoo"]("X,a,b,l,k"));

            if (includes_["u0"]) {
                scalars_["3"]  = 1.00 * dot(m0("X"), s0("X"));
            }

            if (includes_["u1"]) {
                scalars_["4"]  = 1.00 * dot(m1["aa_Lov"]("X,k,a"), s1["aa_Lvo"]("X,a,k"));
            }

            if (includes_["u2"]) {
                scalars_["5"]  = 1.00 * dot(m2["aaaa_Loovv"]("X,l,k,a,b"), s2["aaaa_Lvvoo"]("X,a,b,l,k"));
            }
            scalars_["6"]  = 1.00 * dot(l1["bb_Lov"]("X,k,a"), r1["bb_Lvo"]("X,a,k"));
            scalars_["7"]  = 1.00 * dot(l2["bbbb_Loovv"]("X,l,k,a,b"), r2["bbbb_Lvvoo"]("X,a,b,l,k"));

            if (includes_["u1"]) {
                scalars_["8"]  = 1.00 * dot(m1["bb_Lov"]("X,k,a"), s1["bb_Lvo"]("X,a,k"));
            }

            if (includes_["u2"]) {
                scalars_["9"]  = 1.00 * dot(m2["bbbb_Loovv"]("X,l,k,a,b"), s2["bbbb_Lvvoo"]("X,a,b,l,k"));
            }

            if (includes_["u0"]) {
                scalars_["10"]  = 1.00 * dot(m0("X"), r0("X"));
            }

            if (includes_["u1"]) {
                scalars_["11"]  = 1.00 * dot(m1["aa_Lov"]("X,j,b"), r1["aa_Lvo"]("X,b,j"));
            }

            if (includes_["u2"]) {
                scalars_["12"]  = 1.00 * dot(m2["aaaa_Loovv"]("X,j,k,c,b"), r2["aaaa_Lvvoo"]("X,c,b,j,k"));
                scalars_["13"]  = 1.00 * dot(m2["abab_Loovv"]("X,j,k,b,c"), r2["abab_Lvvoo"]("X,b,c,j,k"));
            }

            if (includes_["u1"]) {
                scalars_["14"]  = 1.00 * dot(m1["bb_Lov"]("X,j,b"), r1["bb_Lvo"]("X,b,j"));
            }

            if (includes_["u2"]) {
                scalars_["15"]  = 1.00 * dot(m2["abab_Loovv"]("X,k,j,c,b"), r2["abab_Lvvoo"]("X,c,b,k,j"));
                scalars_["16"]  = 1.00 * dot(m2["bbbb_Loovv"]("X,j,k,c,b"), r2["bbbb_Lvvoo"]("X,c,b,j,k"));
            }

            // rdm_oo_aa = +1.00 d_aa(i,j) r0 l0  // flops: o2v0 = o2v0 | mem: o2v0 = o2v0
            rdm_oo_aa("i,j")  = scalars_["0"] * Id["aa_oo"]("i,j");

            // rdm_oo_aa += +0.25 d_aa(i,j) r2_aaaa(a,b,l,k) l2_aaaa(l,k,a,b)  // flops: o2v0 += o2v0 | mem: o2v0 += o2v0
            rdm_oo_aa("i,j") += 0.25 * scalars_["2"] * Id["aa_oo"]("i,j");

            // rdm_oo_aa += +1.00 d_aa(i,j) r1_aa(a,k) l1_aa(k,a)  // flops: o2v0 += o2v0 | mem: o2v0 += o2v0
            rdm_oo_aa("i,j") += scalars_["1"] * Id["aa_oo"]("i,j");

            // rdm_oo_bb = +1.00 d_bb(i,j) r0 l0  // flops: o2v0 = o2v0 | mem: o2v0 = o2v0
            rdm_oo_bb("i,j")  = scalars_["0"] * Id["bb_oo"]("i,j");

            // rdm_oo_bb += +1.00 d_bb(i,j) r1_bb(a,k) l1_bb(k,a)  // flops: o2v0 += o2v0 | mem: o2v0 += o2v0
            rdm_oo_bb("i,j") += scalars_["6"] * Id["bb_oo"]("i,j");

            // rdm_oo_bb += +0.25 d_bb(i,j) r2_bbbb(a,b,l,k) l2_bbbb(l,k,a,b)  // flops: o2v0 += o2v0 | mem: o2v0 += o2v0
            rdm_oo_bb("i,j") += 0.25 * scalars_["7"] * Id["bb_oo"]("i,j");

            // rdm_ov_aa = +1.00 r1_aa(a,i) l0  // flops: o1v1 = o1v1L1 | mem: o1v1 = o1v1
            rdm_ov_aa("a,i")  = r1["aa_Lvo"]("X,a,i") * l0("X");

            // rdm_ov_bb = +1.00 r1_bb(a,i) l0  // flops: o1v1 = o1v1L1 | mem: o1v1 = o1v1
            rdm_ov_bb("a,i")  = r1["bb_Lvo"]("X,a,i") * l0("X");

            // rdm_oo_aa += -1.00 r1_aa(a,i) l1_aa(j,a)  // flops: o2v0 += o2v1L1 | mem: o2v0 += o2v0
            rdm_oo_aa("i,j") -= r1["aa_Lvo"]("X,a,i") * l1["aa_Lov"]("X,j,a");

            // rdm_oo_bb += -1.00 r1_bb(a,i) l1_bb(j,a)  // flops: o2v0 += o2v1L1 | mem: o2v0 += o2v0
            rdm_oo_bb("i,j") -= r1["bb_Lvo"]("X,a,i") * l1["bb_Lov"]("X,j,a");

            // rdm_vv_aa = +1.00 r1_aa(b,i) l1_aa(i,a)  // flops: o0v2 = o1v2L1 | mem: o0v2 = o0v2
            rdm_vv_aa("b,a")  = r1["aa_Lvo"]("X,b,i") * l1["aa_Lov"]("X,i,a");

            // rdm_vv_bb = +1.00 r1_bb(b,i) l1_bb(i,a)  // flops: o0v2 = o1v2L1 | mem: o0v2 = o0v2
            rdm_vv_bb("b,a")  = r1["bb_Lvo"]("X,b,i") * l1["bb_Lov"]("X,i,a");

            // rdm_ov_aa += -1.00 r2_aaaa(a,b,j,i) l1_aa(j,b)  // flops: o1v1 += o2v2L1 | mem: o1v1 += o1v1
            rdm_ov_aa("a,i") -= r2["aaaa_Lvvoo"]("X,a,b,j,i") * l1["aa_Lov"]("X,j,b");

            // rdm_ov_bb += -1.00 r2_bbbb(a,b,j,i) l1_bb(j,b)  // flops: o1v1 += o2v2L1 | mem: o1v1 += o1v1
            rdm_ov_bb("a,i") -= r2["bbbb_Lvvoo"]("X,a,b,j,i") * l1["bb_Lov"]("X,j,b");

            // rdm_oo_aa += -0.50 r2_aaaa(a,b,k,i) l2_aaaa(k,j,a,b)  // flops: o2v0 += o3v2L1 | mem: o2v0 += o2v0
            rdm_oo_aa("i,j") -= 0.50 * r2["aaaa_Lvvoo"]("X,a,b,k,i") * l2["aaaa_Loovv"]("X,k,j,a,b");

            // rdm_oo_bb += -0.50 r2_bbbb(a,b,k,i) l2_bbbb(k,j,a,b)  // flops: o2v0 += o3v2L1 | mem: o2v0 += o2v0
            rdm_oo_bb("i,j") -= 0.50 * r2["bbbb_Lvvoo"]("X,a,b,k,i") * l2["bbbb_Loovv"]("X,k,j,a,b");

            // rdm_vv_aa += +0.50 r2_aaaa(b,c,j,i) l2_aaaa(j,i,a,c)  // flops: o0v2 += o2v3L1 | mem: o0v2 += o0v2
            rdm_vv_aa("b,a") += 0.50 * r2["aaaa_Lvvoo"]("X,b,c,j,i") * l2["aaaa_Loovv"]("X,j,i,a,c");

            // rdm_vv_bb += +0.50 r2_bbbb(b,c,j,i) l2_bbbb(j,i,a,c)  // flops: o0v2 += o2v3L1 | mem: o0v2 += o0v2
            rdm_vv_bb("b,a") += 0.50 * r2["bbbb_Lvvoo"]("X,b,c,j,i") * l2["bbbb_Loovv"]("X,j,i,a,c");

            // rdm_ov_aa += -1.00 r1_bb(c,k) t2_aaaa(a,b,j,i) l2_abab(j,k,b,c)  // flops: o1v1 += o2v2L1 o2v2 | mem: o1v1 += o1v1 o1v1
            rdm_ov_aa("a,i") -= r1["bb_Lvo"]("X,c,k") * l2["abab_Loovv"]("X,j,k,b,c") * t2["aaaa_vvoo"]("a,b,j,i");

            // rdm_ov_bb += -1.00 r1_aa(c,k) t2_bbbb(a,b,j,i) l2_abab(k,j,c,b)  // flops: o1v1 += o2v2L1 o2v2 | mem: o1v1 += o1v1 o1v1
            rdm_ov_bb("a,i") -= r1["aa_Lvo"]("X,c,k") * l2["abab_Loovv"]("X,k,j,c,b") * t2["bbbb_vvoo"]("a,b,j,i");

            // rdm_ov_aa += +0.50 r1_aa(a,k) t2_aaaa(c,b,j,i) l2_aaaa(k,j,c,b)  // flops: o1v1 += o3v2L1 o2v1L1 | mem: o1v1 += o2v0L1 o1v1
            rdm_ov_aa("a,i") += 0.50 * t2["aaaa_vvoo"]("c,b,j,i") * l2["aaaa_Loovv"]("X,k,j,c,b") * r1["aa_Lvo"]("X,a,k");

            // rdm_ov_bb += +0.50 r1_bb(a,k) t2_bbbb(c,b,j,i) l2_bbbb(k,j,c,b)  // flops: o1v1 += o3v2L1 o2v1L1 | mem: o1v1 += o2v0L1 o1v1
            rdm_ov_bb("a,i") += 0.50 * t2["bbbb_vvoo"]("c,b,j,i") * l2["bbbb_Loovv"]("X,k,j,c,b") * r1["bb_Lvo"]("X,a,k");

            // rdm_ov_aa += +0.50 r1_aa(c,i) t2_aaaa(a,b,j,k) l2_aaaa(j,k,b,c)  // flops: o1v1 += o3v2L1 o3v2 | mem: o1v1 += o3v1 o1v1
            rdm_ov_aa("a,i") += 0.50 * r1["aa_Lvo"]("X,c,i") * l2["aaaa_Loovv"]("X,j,k,b,c") * t2["aaaa_vvoo"]("a,b,j,k");

            // rdm_ov_bb += +0.50 r1_bb(c,i) t2_bbbb(a,b,j,k) l2_bbbb(j,k,b,c)  // flops: o1v1 += o3v2L1 o3v2 | mem: o1v1 += o3v1 o1v1
            rdm_ov_bb("a,i") += 0.50 * r1["bb_Lvo"]("X,c,i") * l2["bbbb_Loovv"]("X,j,k,b,c") * t2["bbbb_vvoo"]("a,b,j,k");

            if (includes_["u0"]) {

                // rdm_oo_aa += +1.00 d_aa(i,j) s0 m0  // flops: o2v0 += o2v0 | mem: o2v0 += o2v0
                rdm_oo_aa("i,j") += scalars_["3"] * Id["aa_oo"]("i,j");

                // rdm_oo_bb += +1.00 d_bb(i,j) s0 m0  // flops: o2v0 += o2v0 | mem: o2v0 += o2v0
                rdm_oo_bb("i,j") += scalars_["3"] * Id["bb_oo"]("i,j");
            }

            if (includes_["u0"] && includes_["u1"]) {

                // rdm_ov_aa += +1.00 r0 u1_aa(a,i) m0  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
                rdm_ov_aa("a,i") += scalars_["10"] * u1["aa_vo"]("a,i");

                // rdm_ov_bb += +1.00 r0 u1_bb(a,i) m0  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
                rdm_ov_bb("a,i") += scalars_["10"] * u1["bb_vo"]("a,i");

                // rdm_ov_aa += +1.00 s1_aa(a,i) m0  // flops: o1v1 += o1v1L1 | mem: o1v1 += o1v1
                rdm_ov_aa("a,i") += s1["aa_Lvo"]("X,a,i") * m0("X");

                // rdm_ov_bb += +1.00 s1_bb(a,i) m0  // flops: o1v1 += o1v1L1 | mem: o1v1 += o1v1
                rdm_ov_bb("a,i") += s1["bb_Lvo"]("X,a,i") * m0("X");
            }

            if (includes_["u1"]) {

                // rdm_oo_aa += +1.00 d_aa(i,j) s1_aa(a,k) m1_aa(k,a)  // flops: o2v0 += o2v0 | mem: o2v0 += o2v0
                rdm_oo_aa("i,j") += scalars_["4"] * Id["aa_oo"]("i,j");

                // rdm_oo_bb += +1.00 d_bb(i,j) s1_bb(a,k) m1_bb(k,a)  // flops: o2v0 += o2v0 | mem: o2v0 += o2v0
                rdm_oo_bb("i,j") += scalars_["8"] * Id["bb_oo"]("i,j");

                // rdm_ov_aa += +1.00 r1_aa(b,j) u1_aa(a,i) m1_aa(j,b)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
                rdm_ov_aa("a,i") += scalars_["11"] * u1["aa_vo"]("a,i");

                // rdm_ov_bb += +1.00 r1_bb(b,j) u1_bb(a,i) m1_bb(j,b)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
                rdm_ov_bb("a,i") += scalars_["14"] * u1["bb_vo"]("a,i");

                // rdm_oo_aa += -1.00 s1_aa(a,i) m1_aa(j,a)  // flops: o2v0 += o2v1L1 | mem: o2v0 += o2v0
                rdm_oo_aa("i,j") -= s1["aa_Lvo"]("X,a,i") * m1["aa_Lov"]("X,j,a");

                // rdm_oo_bb += -1.00 s1_bb(a,i) m1_bb(j,a)  // flops: o2v0 += o2v1L1 | mem: o2v0 += o2v0
                rdm_oo_bb("i,j") -= s1["bb_Lvo"]("X,a,i") * m1["bb_Lov"]("X,j,a");

                // rdm_vv_aa += +1.00 s1_aa(b,i) m1_aa(i,a)  // flops: o0v2 += o1v2L1 | mem: o0v2 += o0v2
                rdm_vv_aa("b,a") += s1["aa_Lvo"]("X,b,i") * m1["aa_Lov"]("X,i,a");

                // rdm_vv_bb += +1.00 s1_bb(b,i) m1_bb(i,a)  // flops: o0v2 += o1v2L1 | mem: o0v2 += o0v2
                rdm_vv_bb("b,a") += s1["bb_Lvo"]("X,b,i") * m1["bb_Lov"]("X,i,a");

                // rdm_ov_aa += -1.00 r1_aa(b,i) u1_aa(a,j) m1_aa(j,b)  // flops: o1v1 += o2v1L1 o2v1 | mem: o1v1 += o2v0 o1v1
                rdm_ov_aa("a,i") -= r1["aa_Lvo"]("X,b,i") * m1["aa_Lov"]("X,j,b") * u1["aa_vo"]("a,j");

                // rdm_ov_bb += -1.00 r1_bb(b,i) u1_bb(a,j) m1_bb(j,b)  // flops: o1v1 += o2v1L1 o2v1 | mem: o1v1 += o2v0 o1v1
                rdm_ov_bb("a,i") -= r1["bb_Lvo"]("X,b,i") * m1["bb_Lov"]("X,j,b") * u1["bb_vo"]("a,j");

                // rdm_ov_aa += -1.00 r1_aa(a,j) u1_aa(b,i) m1_aa(j,b)  // flops: o1v1 += o2v1L1 o2v1L1 | mem: o1v1 += o2v0L1 o1v1
                rdm_ov_aa("a,i") -= u1["aa_vo"]("b,i") * m1["aa_Lov"]("X,j,b") * r1["aa_Lvo"]("X,a,j");

                // rdm_ov_bb += -1.00 r1_bb(a,j) u1_bb(b,i) m1_bb(j,b)  // flops: o1v1 += o2v1L1 o2v1L1 | mem: o1v1 += o2v0L1 o1v1
                rdm_ov_bb("a,i") -= u1["bb_vo"]("b,i") * m1["bb_Lov"]("X,j,b") * r1["bb_Lvo"]("X,a,j");
            }

            if (includes_["u1"] && includes_["u2"]) {

                // rdm_ov_aa += +0.25 r2_abab(b,c,j,k) u1_aa(a,i) m2_abab(j,k,b,c)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
                rdm_ov_aa("a,i") += 0.25 * scalars_["13"] * u1["aa_vo"]("a,i");

                // rdm_ov_aa += +0.25 r2_aaaa(c,b,j,k) u1_aa(a,i) m2_aaaa(j,k,c,b)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
                rdm_ov_aa("a,i") += 0.25 * scalars_["12"] * u1["aa_vo"]("a,i");

                // rdm_ov_bb += +0.25 r2_bbbb(c,b,j,k) u1_bb(a,i) m2_bbbb(j,k,c,b)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
                rdm_ov_bb("a,i") += 0.25 * scalars_["16"] * u1["bb_vo"]("a,i");

                // rdm_ov_bb += +0.25 r2_abab(c,b,k,j) u1_bb(a,i) m2_abab(k,j,c,b)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
                rdm_ov_bb("a,i") += 0.25 * scalars_["15"] * u1["bb_vo"]("a,i");

                // rdm_ov_aa += -1.00 s2_aaaa(a,b,j,i) m1_aa(j,b)  // flops: o1v1 += o2v2L1 | mem: o1v1 += o1v1
                rdm_ov_aa("a,i") -= s2["aaaa_Lvvoo"]("X,a,b,j,i") * m1["aa_Lov"]("X,j,b");

                // rdm_ov_bb += -1.00 s2_bbbb(a,b,j,i) m1_bb(j,b)  // flops: o1v1 += o2v2L1 | mem: o1v1 += o1v1
                rdm_ov_bb("a,i") -= s2["bbbb_Lvvoo"]("X,a,b,j,i") * m1["bb_Lov"]("X,j,b");

                // rdm_ov_aa += -0.50 r2_aaaa(c,b,j,i) u1_aa(a,k) m2_aaaa(j,k,c,b)  // flops: o1v1 += o3v2L1 o2v1 | mem: o1v1 += o2v0 o1v1
                rdm_ov_aa("a,i") -= 0.50 * r2["aaaa_Lvvoo"]("X,c,b,j,i") * m2["aaaa_Loovv"]("X,j,k,c,b") * u1["aa_vo"]("a,k");

                // rdm_ov_bb += -0.50 r2_bbbb(c,b,j,i) u1_bb(a,k) m2_bbbb(j,k,c,b)  // flops: o1v1 += o3v2L1 o2v1 | mem: o1v1 += o2v0 o1v1
                rdm_ov_bb("a,i") -= 0.50 * r2["bbbb_Lvvoo"]("X,c,b,j,i") * m2["bbbb_Loovv"]("X,j,k,c,b") * u1["bb_vo"]("a,k");

                // rdm_ov_aa += -1.00 s1_bb(c,k) t2_aaaa(a,b,j,i) m2_abab(j,k,b,c)  // flops: o1v1 += o2v2L1 o2v2 | mem: o1v1 += o1v1 o1v1
                rdm_ov_aa("a,i") -= s1["bb_Lvo"]("X,c,k") * m2["abab_Loovv"]("X,j,k,b,c") * t2["aaaa_vvoo"]("a,b,j,i");

                // rdm_ov_bb += -1.00 s1_aa(c,k) t2_bbbb(a,b,j,i) m2_abab(k,j,c,b)  // flops: o1v1 += o2v2L1 o2v2 | mem: o1v1 += o1v1 o1v1
                rdm_ov_bb("a,i") -= s1["aa_Lvo"]("X,c,k") * m2["abab_Loovv"]("X,k,j,c,b") * t2["bbbb_vvoo"]("a,b,j,i");

                // rdm_ov_aa += +0.50 s1_aa(a,k) t2_aaaa(c,b,j,i) m2_aaaa(k,j,c,b)  // flops: o1v1 += o3v2L1 o2v1L1 | mem: o1v1 += o2v0L1 o1v1
                rdm_ov_aa("a,i") += 0.50 * t2["aaaa_vvoo"]("c,b,j,i") * m2["aaaa_Loovv"]("X,k,j,c,b") * s1["aa_Lvo"]("X,a,k");

                // rdm_ov_bb += +0.50 s1_bb(a,k) t2_bbbb(c,b,j,i) m2_bbbb(k,j,c,b)  // flops: o1v1 += o3v2L1 o2v1L1 | mem: o1v1 += o2v0L1 o1v1
                rdm_ov_bb("a,i") += 0.50 * t2["bbbb_vvoo"]("c,b,j,i") * m2["bbbb_Loovv"]("X,k,j,c,b") * s1["bb_Lvo"]("X,a,k");

                // rdm_ov_aa += +0.50 s1_aa(c,i) t2_aaaa(a,b,j,k) m2_aaaa(j,k,b,c)  // flops: o1v1 += o3v2L1 o3v2 | mem: o1v1 += o3v1 o1v1
                rdm_ov_aa("a,i") += 0.50 * s1["aa_Lvo"]("X,c,i") * m2["aaaa_Loovv"]("X,j,k,b,c") * t2["aaaa_vvoo"]("a,b,j,k");

                // rdm_ov_bb += +0.50 s1_bb(c,i) t2_bbbb(a,b,j,k) m2_bbbb(j,k,b,c)  // flops: o1v1 += o3v2L1 o3v2 | mem: o1v1 += o3v1 o1v1
                rdm_ov_bb("a,i") += 0.50 * s1["bb_Lvo"]("X,c,i") * m2["bbbb_Loovv"]("X,j,k,b,c") * t2["bbbb_vvoo"]("a,b,j,k");

                // rdm_ov_aa += -0.50 r2_aaaa(a,b,j,k) u1_aa(c,i) m2_aaaa(j,k,c,b)  // flops: o1v1 += o3v2L1 o3v2L1 | mem: o1v1 += o3v1L1 o1v1
                rdm_ov_aa("a,i") -= 0.50 * u1["aa_vo"]("c,i") * m2["aaaa_Loovv"]("X,j,k,c,b") * r2["aaaa_Lvvoo"]("X,a,b,j,k");

                // rdm_ov_bb += -0.50 r2_bbbb(a,b,j,k) u1_bb(c,i) m2_bbbb(j,k,c,b)  // flops: o1v1 += o3v2L1 o3v2L1 | mem: o1v1 += o3v1L1 o1v1
                rdm_ov_bb("a,i") -= 0.50 * u1["bb_vo"]("c,i") * m2["bbbb_Loovv"]("X,j,k,c,b") * r2["bbbb_Lvvoo"]("X,a,b,j,k");
            }

            if (includes_["u2"]) {

                // rdm_oo_aa += +0.25 d_aa(i,j) s2_aaaa(a,b,l,k) m2_aaaa(l,k,a,b)  // flops: o2v0 += o2v0 | mem: o2v0 += o2v0
                rdm_oo_aa("i,j") += 0.25 * scalars_["5"] * Id["aa_oo"]("i,j");

                // rdm_oo_bb += +0.25 d_bb(i,j) s2_bbbb(a,b,l,k) m2_bbbb(l,k,a,b)  // flops: o2v0 += o2v0 | mem: o2v0 += o2v0
                rdm_oo_bb("i,j") += 0.25 * scalars_["9"] * Id["bb_oo"]("i,j");

                // rdm_oo_aa += -0.50 s2_aaaa(a,b,k,i) m2_aaaa(k,j,a,b)  // flops: o2v0 += o3v2L1 | mem: o2v0 += o2v0
                rdm_oo_aa("i,j") -= 0.50 * s2["aaaa_Lvvoo"]("X,a,b,k,i") * m2["aaaa_Loovv"]("X,k,j,a,b");

                // rdm_oo_bb += -0.50 s2_bbbb(a,b,k,i) m2_bbbb(k,j,a,b)  // flops: o2v0 += o3v2L1 | mem: o2v0 += o2v0
                rdm_oo_bb("i,j") -= 0.50 * s2["bbbb_Lvvoo"]("X,a,b,k,i") * m2["bbbb_Loovv"]("X,k,j,a,b");

                // rdm_vv_aa += +0.50 s2_aaaa(b,c,j,i) m2_aaaa(j,i,a,c)  // flops: o0v2 += o2v3L1 | mem: o0v2 += o0v2
                rdm_vv_aa("b,a") += 0.50 * s2["aaaa_Lvvoo"]("X,b,c,j,i") * m2["aaaa_Loovv"]("X,j,i,a,c");

                // rdm_vv_bb += +0.50 s2_bbbb(b,c,j,i) m2_bbbb(j,i,a,c)  // flops: o0v2 += o2v3L1 | mem: o0v2 += o0v2
                rdm_vv_bb("b,a") += 0.50 * s2["bbbb_Lvvoo"]("X,b,c,j,i") * m2["bbbb_Loovv"]("X,j,i,a,c");

                // rdm_ov_aa += +1.00 r1_aa(b,j) u2_abab(a,c,i,k) m2_abab(j,k,b,c)  // flops: o1v1 += o2v2L1 o2v2 | mem: o1v1 += o1v1 o1v1
                rdm_ov_aa("a,i") += r1["aa_Lvo"]("X,b,j") * m2["abab_Loovv"]("X,j,k,b,c") * u2["abab_vvoo"]("a,c,i,k");

                // rdm_ov_bb += +1.00 r1_bb(b,j) u2_abab(c,a,k,i) m2_abab(k,j,c,b)  // flops: o1v1 += o2v2L1 o2v2 | mem: o1v1 += o1v1 o1v1
                rdm_ov_bb("a,i") += r1["bb_Lvo"]("X,b,j") * m2["abab_Loovv"]("X,k,j,c,b") * u2["abab_vvoo"]("c,a,k,i");

                // rdm_ov_aa += -0.50 r1_aa(a,j) u2_abab(b,c,i,k) m2_abab(j,k,b,c)  // flops: o1v1 += o3v2L1 o2v1L1 | mem: o1v1 += o2v0L1 o1v1
                rdm_ov_aa("a,i") -= 0.50 * u2["abab_vvoo"]("b,c,i,k") * m2["abab_Loovv"]("X,j,k,b,c") * r1["aa_Lvo"]("X,a,j");

                // rdm_ov_aa += +0.50 r1_aa(a,j) u2_aaaa(c,b,k,i) m2_aaaa(j,k,c,b)  // flops: o1v1 += o3v2L1 o2v1L1 | mem: o1v1 += o2v0L1 o1v1
                rdm_ov_aa("a,i") += 0.50 * u2["aaaa_vvoo"]("c,b,k,i") * m2["aaaa_Loovv"]("X,j,k,c,b") * r1["aa_Lvo"]("X,a,j");

                // rdm_ov_bb += -0.50 r1_bb(a,j) u2_abab(c,b,k,i) m2_abab(k,j,c,b)  // flops: o1v1 += o3v2L1 o2v1L1 | mem: o1v1 += o2v0L1 o1v1
                rdm_ov_bb("a,i") -= 0.50 * u2["abab_vvoo"]("c,b,k,i") * m2["abab_Loovv"]("X,k,j,c,b") * r1["bb_Lvo"]("X,a,j");

                // rdm_ov_bb += +0.50 r1_bb(a,j) u2_bbbb(c,b,k,i) m2_bbbb(j,k,c,b)  // flops: o1v1 += o3v2L1 o2v1L1 | mem: o1v1 += o2v0L1 o1v1
                rdm_ov_bb("a,i") += 0.50 * u2["bbbb_vvoo"]("c,b,k,i") * m2["bbbb_Loovv"]("X,j,k,c,b") * r1["bb_Lvo"]("X,a,j");

                // rdm_ov_aa += -0.50 r1_aa(b,i) u2_abab(a,c,j,k) m2_abab(j,k,b,c)  // flops: o1v1 += o3v2L1 o3v2 | mem: o1v1 += o3v1 o1v1
                rdm_ov_aa("a,i") -= 0.50 * r1["aa_Lvo"]("X,b,i") * m2["abab_Loovv"]("X,j,k,b,c") * u2["abab_vvoo"]("a,c,j,k");

                // rdm_ov_aa += +0.50 r1_aa(b,i) u2_aaaa(a,c,j,k) m2_aaaa(j,k,c,b)  // flops: o1v1 += o3v2L1 o3v2 | mem: o1v1 += o3v1 o1v1
                rdm_ov_aa("a,i") += 0.50 * r1["aa_Lvo"]("X,b,i") * m2["aaaa_Loovv"]("X,j,k,c,b") * u2["aaaa_vvoo"]("a,c,j,k");

                // rdm_ov_bb += -0.50 r1_bb(b,i) u2_abab(c,a,k,j) m2_abab(k,j,c,b)  // flops: o1v1 += o3v2L1 o3v2 | mem: o1v1 += o3v1 o1v1
                rdm_ov_bb("a,i") -= 0.50 * r1["bb_Lvo"]("X,b,i") * m2["abab_Loovv"]("X,k,j,c,b") * u2["abab_vvoo"]("c,a,k,j");

                // rdm_ov_bb += +0.50 r1_bb(b,i) u2_bbbb(a,c,j,k) m2_bbbb(j,k,c,b)  // flops: o1v1 += o3v2L1 o3v2 | mem: o1v1 += o3v1 o1v1
                rdm_ov_bb("a,i") += 0.50 * r1["bb_Lvo"]("X,b,i") * m2["bbbb_Loovv"]("X,j,k,c,b") * u2["bbbb_vvoo"]("a,c,j,k");

                // tmps_[1_bbbb_vvoo](c,a,j,i) = 0.50 m2[bbbb_Loovv](X,i,j,a,c) * r0(X) // flops: o2v2 = o2v2L1 | mem: o2v2 = o2v2
                tmps_["1_bbbb_vvoo"]("c,a,j,i")  = 0.50 * m2["bbbb_Loovv"]("X,i,j,a,c") * r0("X");
            }

            if (includes_["u1"] && includes_["u2"]) {

                // rdm_ov_bb += +0.50 r0 u1_bb(a,j) t2_bbbb(c,b,k,i) m2_bbbb(j,k,c,b)  // flops: o1v1 += o3v2 o2v1 | mem: o1v1 += o2v0 o1v1
                rdm_ov_bb("a,i") += tmps_["1_bbbb_vvoo"]("b,c,k,j") * t2["bbbb_vvoo"]("c,b,k,i") * u1["bb_vo"]("a,j");

                // rdm_ov_bb += +0.50 r0 u1_bb(b,i) t2_bbbb(a,c,j,k) m2_bbbb(j,k,c,b)  // flops: o1v1 += o3v2 o3v2 | mem: o1v1 += o3v1 o1v1
                rdm_ov_bb("a,i") += tmps_["1_bbbb_vvoo"]("b,c,k,j") * u1["bb_vo"]("b,i") * t2["bbbb_vvoo"]("a,c,j,k");
            }

            if (includes_["u2"]) {

                // rdm_oo_bb += -0.50 r0 u2_bbbb(b,a,k,i) m2_bbbb(k,j,b,a)  // flops: o2v0 += o3v2 | mem: o2v0 += o2v0
                rdm_oo_bb("j,i") -= tmps_["1_bbbb_vvoo"]("a,b,j,k") * u2["bbbb_vvoo"]("b,a,k,i");

                // rdm_vv_bb += +0.50 r0 u2_bbbb(b,c,i,j) m2_bbbb(i,j,a,c)  // flops: o0v2 += o2v3 | mem: o0v2 += o0v2
                rdm_vv_bb("a,b") += tmps_["1_bbbb_vvoo"]("c,a,j,i") * u2["bbbb_vvoo"]("b,c,i,j");
                tmps_["1_bbbb_vvoo"].~TArrayD();

                // tmps_[2_aaaa_vvoo](c,a,j,i) = 0.50 m2[aaaa_Loovv](X,i,j,a,c) * r0(X) // flops: o2v2 = o2v2L1 | mem: o2v2 = o2v2
                tmps_["2_aaaa_vvoo"]("c,a,j,i")  = 0.50 * m2["aaaa_Loovv"]("X,i,j,a,c") * r0("X");
            }

            if (includes_["u1"] && includes_["u2"]) {

                // rdm_ov_aa += +0.50 r0 u1_aa(a,j) t2_aaaa(c,b,k,i) m2_aaaa(j,k,c,b)  // flops: o1v1 += o3v2 o2v1 | mem: o1v1 += o2v0 o1v1
                rdm_ov_aa("a,i") += tmps_["2_aaaa_vvoo"]("b,c,k,j") * t2["aaaa_vvoo"]("c,b,k,i") * u1["aa_vo"]("a,j");

                // rdm_ov_aa += +0.50 r0 u1_aa(b,i) t2_aaaa(a,c,j,k) m2_aaaa(j,k,c,b)  // flops: o1v1 += o3v2 o3v2 | mem: o1v1 += o3v1 o1v1
                rdm_ov_aa("a,i") += tmps_["2_aaaa_vvoo"]("b,c,k,j") * u1["aa_vo"]("b,i") * t2["aaaa_vvoo"]("a,c,j,k");
            }

            if (includes_["u2"]) {

                // rdm_oo_aa += -0.50 r0 u2_aaaa(b,a,k,i) m2_aaaa(k,j,b,a)  // flops: o2v0 += o3v2 | mem: o2v0 += o2v0
                rdm_oo_aa("j,i") -= tmps_["2_aaaa_vvoo"]("a,b,j,k") * u2["aaaa_vvoo"]("b,a,k,i");

                // rdm_vv_aa += +0.50 r0 u2_aaaa(b,c,i,j) m2_aaaa(i,j,a,c)  // flops: o0v2 += o2v3 | mem: o0v2 += o0v2
                rdm_vv_aa("a,b") += tmps_["2_aaaa_vvoo"]("c,a,j,i") * u2["aaaa_vvoo"]("b,c,i,j");
                tmps_["2_aaaa_vvoo"].~TArrayD();

                // tmps_[3_abab_vvoo](c,b,k,j) = 0.50 m2[abab_Loovv](X,k,j,c,b) * r0(X) // flops: o2v2 = o2v2L1 | mem: o2v2 = o2v2
                tmps_["3_abab_vvoo"]("c,b,k,j")  = 0.50 * m2["abab_Loovv"]("X,k,j,c,b") * r0("X");
            }

            if (includes_["u1"] && includes_["u2"]) {

                // rdm_ov_aa += -0.50 r0 u1_aa(a,j) t2_abab(b,c,i,k) m2_abab(j,k,b,c)  // flops: o1v1 += o3v2 o2v1 | mem: o1v1 += o2v0 o1v1
                rdm_ov_aa("a,i") -= tmps_["3_abab_vvoo"]("b,c,j,k") * t2["abab_vvoo"]("b,c,i,k") * u1["aa_vo"]("a,j");

                // rdm_ov_bb += -0.50 r0 u1_bb(a,j) t2_abab(c,b,k,i) m2_abab(k,j,c,b)  // flops: o1v1 += o3v2 o2v1 | mem: o1v1 += o2v0 o1v1
                rdm_ov_bb("a,i") -= tmps_["3_abab_vvoo"]("c,b,k,j") * t2["abab_vvoo"]("c,b,k,i") * u1["bb_vo"]("a,j");

                // rdm_ov_aa += -0.50 r0 u1_aa(b,i) t2_abab(a,c,j,k) m2_abab(j,k,b,c)  // flops: o1v1 += o3v2 o3v2 | mem: o1v1 += o3v1 o1v1
                rdm_ov_aa("a,i") -= tmps_["3_abab_vvoo"]("b,c,j,k") * u1["aa_vo"]("b,i") * t2["abab_vvoo"]("a,c,j,k");

                // rdm_ov_bb += -0.50 r0 u1_bb(b,i) t2_abab(c,a,k,j) m2_abab(k,j,c,b)  // flops: o1v1 += o3v2 o3v2 | mem: o1v1 += o3v1 o1v1
                rdm_ov_bb("a,i") -= tmps_["3_abab_vvoo"]("c,b,k,j") * u1["bb_vo"]("b,i") * t2["abab_vvoo"]("c,a,k,j");
            }

            if (includes_["u2"]) {
                tmps_["3_abab_vvoo"].~TArrayD();

                // tmps_[4_bb_vo](a,j) = 1.00 m2[bbbb_Loovv](X,i,j,a,c) * r1[bb_Lvo](X,c,i) // flops: o1v1 = o2v2L1 | mem: o1v1 = o1v1
                tmps_["4_bb_vo"]("a,j")  = m2["bbbb_Loovv"]("X,i,j,a,c") * r1["bb_Lvo"]("X,c,i");
            }

            if (includes_["u1"] && includes_["u2"]) {

                // rdm_oo_bb += +1.00 r1_bb(a,k) u1_bb(b,i) m2_bbbb(k,j,b,a)  // flops: o2v0 += o2v1 | mem: o2v0 += o2v0
                rdm_oo_bb("j,i") += tmps_["4_bb_vo"]("b,j") * u1["bb_vo"]("b,i");

                // rdm_vv_bb += -1.00 r1_bb(c,i) u1_bb(b,j) m2_bbbb(i,j,a,c)  // flops: o0v2 += o1v2 | mem: o0v2 += o0v2
                rdm_vv_bb("a,b") -= tmps_["4_bb_vo"]("a,j") * u1["bb_vo"]("b,j");
            }

            if (includes_["u2"]) {

                // rdm_ov_bb += +1.00 r1_bb(b,j) u2_bbbb(a,c,k,i) m2_bbbb(j,k,c,b)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
                rdm_ov_bb("a,i") += tmps_["4_bb_vo"]("c,k") * u2["bbbb_vvoo"]("a,c,k,i");
                tmps_["4_bb_vo"].~TArrayD();

                // tmps_[5_aa_vo](a,j) = 1.00 m2[aaaa_Loovv](X,i,j,a,c) * r1[aa_Lvo](X,c,i) // flops: o1v1 = o2v2L1 | mem: o1v1 = o1v1
                tmps_["5_aa_vo"]("a,j")  = m2["aaaa_Loovv"]("X,i,j,a,c") * r1["aa_Lvo"]("X,c,i");
            }

            if (includes_["u1"] && includes_["u2"]) {

                // rdm_oo_aa += +1.00 r1_aa(a,k) u1_aa(b,i) m2_aaaa(k,j,b,a)  // flops: o2v0 += o2v1 | mem: o2v0 += o2v0
                rdm_oo_aa("j,i") += tmps_["5_aa_vo"]("b,j") * u1["aa_vo"]("b,i");

                // rdm_vv_aa += -1.00 r1_aa(c,i) u1_aa(b,j) m2_aaaa(i,j,a,c)  // flops: o0v2 += o1v2 | mem: o0v2 += o0v2
                rdm_vv_aa("a,b") -= tmps_["5_aa_vo"]("a,j") * u1["aa_vo"]("b,j");
            }

            if (includes_["u2"]) {

                // rdm_ov_aa += +1.00 r1_aa(b,j) u2_aaaa(a,c,k,i) m2_aaaa(j,k,c,b)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
                rdm_ov_aa("a,i") += tmps_["5_aa_vo"]("c,k") * u2["aaaa_vvoo"]("a,c,k,i");
                tmps_["5_aa_vo"].~TArrayD();
            }

            // tmps_[6_bb_vo](b,j) = 1.00 l2[bbbb_Loovv](X,k,j,b,c) * r1[bb_Lvo](X,c,k) // flops: o1v1 = o2v2L1 | mem: o1v1 = o1v1
            tmps_["6_bb_vo"]("b,j")  = l2["bbbb_Loovv"]("X,k,j,b,c") * r1["bb_Lvo"]("X,c,k");

            // rdm_vo_bb = -1.00 r1_bb(b,j) l2_bbbb(j,i,a,b)  // flops: o1v1 = o1v1 | mem: o1v1 = o1v1
            rdm_vo_bb("a,i")  = -1.00 * tmps_["6_bb_vo"]("a,i");

            // rdm_ov_bb += +1.00 r1_bb(c,k) t2_bbbb(a,b,j,i) l2_bbbb(k,j,b,c)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
            rdm_ov_bb("a,i") += tmps_["6_bb_vo"]("b,j") * t2["bbbb_vvoo"]("a,b,j,i");
            tmps_["6_bb_vo"].~TArrayD();

            if (includes_["u1"] && includes_["u2"]) {

                // tmps_[7_bb_vo](b,j) = 1.00 m2[bbbb_Loovv](X,k,j,b,c) * s1[bb_Lvo](X,c,k) // flops: o1v1 = o2v2L1 | mem: o1v1 = o1v1
                tmps_["7_bb_vo"]("b,j")  = m2["bbbb_Loovv"]("X,k,j,b,c") * s1["bb_Lvo"]("X,c,k");

                // rdm_vo_bb += -1.00 s1_bb(b,j) m2_bbbb(j,i,a,b)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
                rdm_vo_bb("a,i") -= tmps_["7_bb_vo"]("a,i");

                // rdm_ov_bb += +1.00 s1_bb(c,k) t2_bbbb(a,b,j,i) m2_bbbb(k,j,b,c)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
                rdm_ov_bb("a,i") += tmps_["7_bb_vo"]("b,j") * t2["bbbb_vvoo"]("a,b,j,i");
                tmps_["7_bb_vo"].~TArrayD();
            }

            // tmps_[8_aa_vo](b,j) = 1.00 l2[aaaa_Loovv](X,k,j,b,c) * r1[aa_Lvo](X,c,k) // flops: o1v1 = o2v2L1 | mem: o1v1 = o1v1
            tmps_["8_aa_vo"]("b,j")  = l2["aaaa_Loovv"]("X,k,j,b,c") * r1["aa_Lvo"]("X,c,k");

            // rdm_vo_aa = -1.00 r1_aa(b,j) l2_aaaa(j,i,a,b)  // flops: o1v1 = o1v1 | mem: o1v1 = o1v1
            rdm_vo_aa("a,i")  = -1.00 * tmps_["8_aa_vo"]("a,i");

            // rdm_ov_aa += +1.00 r1_aa(c,k) t2_aaaa(a,b,j,i) l2_aaaa(k,j,b,c)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
            rdm_ov_aa("a,i") += tmps_["8_aa_vo"]("b,j") * t2["aaaa_vvoo"]("a,b,j,i");
            tmps_["8_aa_vo"].~TArrayD();

            if (includes_["u1"] && includes_["u2"]) {

                // tmps_[9_aa_vo](b,j) = 1.00 m2[aaaa_Loovv](X,k,j,b,c) * s1[aa_Lvo](X,c,k) // flops: o1v1 = o2v2L1 | mem: o1v1 = o1v1
                tmps_["9_aa_vo"]("b,j")  = m2["aaaa_Loovv"]("X,k,j,b,c") * s1["aa_Lvo"]("X,c,k");

                // rdm_vo_aa += -1.00 s1_aa(b,j) m2_aaaa(j,i,a,b)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
                rdm_vo_aa("a,i") -= tmps_["9_aa_vo"]("a,i");

                // rdm_ov_aa += +1.00 s1_aa(c,k) t2_aaaa(a,b,j,i) m2_aaaa(k,j,b,c)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
                rdm_ov_aa("a,i") += tmps_["9_aa_vo"]("b,j") * t2["aaaa_vvoo"]("a,b,j,i");
                tmps_["9_aa_vo"].~TArrayD();
            }

            // tmps_[10_bbbb_vvoo](c,a,j,i) = 0.50 l2[bbbb_Loovv](X,i,j,a,c) * r0(X) // flops: o2v2 = o2v2L1 | mem: o2v2 = o2v2
            tmps_["10_bbbb_vvoo"]("c,a,j,i")  = 0.50 * l2["bbbb_Loovv"]("X,i,j,a,c") * r0("X");

            // rdm_oo_bb += -0.50 r0 t2_bbbb(b,a,k,i) l2_bbbb(k,j,b,a)  // flops: o2v0 += o3v2 | mem: o2v0 += o2v0
            rdm_oo_bb("j,i") -= tmps_["10_bbbb_vvoo"]("a,b,j,k") * t2["bbbb_vvoo"]("b,a,k,i");

            // rdm_vv_bb += +0.50 r0 t2_bbbb(b,c,i,j) l2_bbbb(i,j,a,c)  // flops: o0v2 += o2v3 | mem: o0v2 += o0v2
            rdm_vv_bb("a,b") += tmps_["10_bbbb_vvoo"]("c,a,j,i") * t2["bbbb_vvoo"]("b,c,i,j");
            tmps_["10_bbbb_vvoo"].~TArrayD();

            if (includes_["u0"] && includes_["u2"]) {

                // tmps_[11_bbbb_vvoo](c,a,j,i) = 0.50 m2[bbbb_Loovv](X,i,j,a,c) * s0(X) // flops: o2v2 = o2v2L1 | mem: o2v2 = o2v2
                tmps_["11_bbbb_vvoo"]("c,a,j,i")  = 0.50 * m2["bbbb_Loovv"]("X,i,j,a,c") * s0("X");

                // rdm_oo_bb += -0.50 s0 t2_bbbb(b,a,k,i) m2_bbbb(k,j,b,a)  // flops: o2v0 += o3v2 | mem: o2v0 += o2v0
                rdm_oo_bb("j,i") -= tmps_["11_bbbb_vvoo"]("a,b,j,k") * t2["bbbb_vvoo"]("b,a,k,i");

                // rdm_vv_bb += +0.50 s0 t2_bbbb(b,c,i,j) m2_bbbb(i,j,a,c)  // flops: o0v2 += o2v3 | mem: o0v2 += o0v2
                rdm_vv_bb("a,b") += tmps_["11_bbbb_vvoo"]("c,a,j,i") * t2["bbbb_vvoo"]("b,c,i,j");
                tmps_["11_bbbb_vvoo"].~TArrayD();
            }

            // tmps_[12_aaaa_vvoo](c,a,j,i) = 0.50 l2[aaaa_Loovv](X,i,j,a,c) * r0(X) // flops: o2v2 = o2v2L1 | mem: o2v2 = o2v2
            tmps_["12_aaaa_vvoo"]("c,a,j,i")  = 0.50 * l2["aaaa_Loovv"]("X,i,j,a,c") * r0("X");

            // rdm_oo_aa += -0.50 r0 t2_aaaa(b,a,k,i) l2_aaaa(k,j,b,a)  // flops: o2v0 += o3v2 | mem: o2v0 += o2v0
            rdm_oo_aa("j,i") -= tmps_["12_aaaa_vvoo"]("a,b,j,k") * t2["aaaa_vvoo"]("b,a,k,i");

            // rdm_vv_aa += +0.50 r0 t2_aaaa(b,c,i,j) l2_aaaa(i,j,a,c)  // flops: o0v2 += o2v3 | mem: o0v2 += o0v2
            rdm_vv_aa("a,b") += tmps_["12_aaaa_vvoo"]("c,a,j,i") * t2["aaaa_vvoo"]("b,c,i,j");
            tmps_["12_aaaa_vvoo"].~TArrayD();

            if (includes_["u0"] && includes_["u2"]) {

                // tmps_[13_aaaa_vvoo](c,a,j,i) = 0.50 m2[aaaa_Loovv](X,i,j,a,c) * s0(X) // flops: o2v2 = o2v2L1 | mem: o2v2 = o2v2
                tmps_["13_aaaa_vvoo"]("c,a,j,i")  = 0.50 * m2["aaaa_Loovv"]("X,i,j,a,c") * s0("X");

                // rdm_oo_aa += -0.50 s0 t2_aaaa(b,a,k,i) m2_aaaa(k,j,b,a)  // flops: o2v0 += o3v2 | mem: o2v0 += o2v0
                rdm_oo_aa("j,i") -= tmps_["13_aaaa_vvoo"]("a,b,j,k") * t2["aaaa_vvoo"]("b,a,k,i");

                // rdm_vv_aa += +0.50 s0 t2_aaaa(b,c,i,j) m2_aaaa(i,j,a,c)  // flops: o0v2 += o2v3 | mem: o0v2 += o0v2
                rdm_vv_aa("a,b") += tmps_["13_aaaa_vvoo"]("c,a,j,i") * t2["aaaa_vvoo"]("b,c,i,j");
                tmps_["13_aaaa_vvoo"].~TArrayD();
            }

            if (includes_["u1"]) {

                // tmps_[14_bb_vo](a,j) = 1.00 m1[bb_Lov](X,j,a) * r0(X) // flops: o1v1 = o1v1L1 | mem: o1v1 = o1v1
                tmps_["14_bb_vo"]("a,j")  = m1["bb_Lov"]("X,j,a") * r0("X");

                // rdm_oo_bb += -1.00 r0 u1_bb(a,i) m1_bb(j,a)  // flops: o2v0 += o2v1 | mem: o2v0 += o2v0
                rdm_oo_bb("j,i") -= tmps_["14_bb_vo"]("a,j") * u1["bb_vo"]("a,i");

                // rdm_vv_bb += +1.00 r0 u1_bb(b,i) m1_bb(i,a)  // flops: o0v2 += o1v2 | mem: o0v2 += o0v2
                rdm_vv_bb("a,b") += tmps_["14_bb_vo"]("a,i") * u1["bb_vo"]("b,i");
            }

            if (includes_["u1"] && includes_["u2"]) {

                // rdm_ov_bb += -1.00 r0 u2_bbbb(a,b,j,i) m1_bb(j,b)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
                rdm_ov_bb("a,i") -= tmps_["14_bb_vo"]("b,j") * u2["bbbb_vvoo"]("a,b,j,i");
            }

            if (includes_["u1"]) {
                tmps_["14_bb_vo"].~TArrayD();

                // tmps_[15_aa_vo](a,j) = 1.00 m1[aa_Lov](X,j,a) * r0(X) // flops: o1v1 = o1v1L1 | mem: o1v1 = o1v1
                tmps_["15_aa_vo"]("a,j")  = m1["aa_Lov"]("X,j,a") * r0("X");

                // rdm_oo_aa += -1.00 r0 u1_aa(a,i) m1_aa(j,a)  // flops: o2v0 += o2v1 | mem: o2v0 += o2v0
                rdm_oo_aa("j,i") -= tmps_["15_aa_vo"]("a,j") * u1["aa_vo"]("a,i");

                // rdm_vv_aa += +1.00 r0 u1_aa(b,i) m1_aa(i,a)  // flops: o0v2 += o1v2 | mem: o0v2 += o0v2
                rdm_vv_aa("a,b") += tmps_["15_aa_vo"]("a,i") * u1["aa_vo"]("b,i");
            }

            if (includes_["u1"] && includes_["u2"]) {

                // rdm_ov_aa += -1.00 r0 u2_aaaa(a,b,j,i) m1_aa(j,b)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
                rdm_ov_aa("a,i") -= tmps_["15_aa_vo"]("b,j") * u2["aaaa_vvoo"]("a,b,j,i");
            }

            if (includes_["u1"]) {
                tmps_["15_aa_vo"].~TArrayD();
            }

            // tmps_[16_bb_vo](b,j) = 1.00 l1[bb_Lov](X,j,b) * r0(X) // flops: o1v1 = o1v1L1 | mem: o1v1 = o1v1
            tmps_["16_bb_vo"]("b,j")  = l1["bb_Lov"]("X,j,b") * r0("X");

            // rdm_vo_bb += +1.00 r0 l1_bb(i,a)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
            rdm_vo_bb("a,i") += tmps_["16_bb_vo"]("a,i");

            // rdm_ov_bb += -1.00 r0 t2_bbbb(a,b,j,i) l1_bb(j,b)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
            rdm_ov_bb("a,i") -= tmps_["16_bb_vo"]("b,j") * t2["bbbb_vvoo"]("a,b,j,i");
            tmps_["16_bb_vo"].~TArrayD();

            // tmps_[17_aa_vo](b,j) = 1.00 l1[aa_Lov](X,j,b) * r0(X) // flops: o1v1 = o1v1L1 | mem: o1v1 = o1v1
            tmps_["17_aa_vo"]("b,j")  = l1["aa_Lov"]("X,j,b") * r0("X");

            // rdm_vo_aa += +1.00 r0 l1_aa(i,a)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
            rdm_vo_aa("a,i") += tmps_["17_aa_vo"]("a,i");

            // rdm_ov_aa += -1.00 r0 t2_aaaa(a,b,j,i) l1_aa(j,b)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
            rdm_ov_aa("a,i") -= tmps_["17_aa_vo"]("b,j") * t2["aaaa_vvoo"]("a,b,j,i");
            tmps_["17_aa_vo"].~TArrayD();

            if (includes_["u0"] && includes_["u1"]) {

                // tmps_[18_aa_vo](b,j) = 1.00 m1[aa_Lov](X,j,b) * s0(X) // flops: o1v1 = o1v1L1 | mem: o1v1 = o1v1
                tmps_["18_aa_vo"]("b,j")  = m1["aa_Lov"]("X,j,b") * s0("X");

                // rdm_vo_aa += +1.00 s0 m1_aa(i,a)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
                rdm_vo_aa("a,i") += tmps_["18_aa_vo"]("a,i");

                // rdm_ov_aa += -1.00 s0 t2_aaaa(a,b,j,i) m1_aa(j,b)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
                rdm_ov_aa("a,i") -= tmps_["18_aa_vo"]("b,j") * t2["aaaa_vvoo"]("a,b,j,i");
                tmps_["18_aa_vo"].~TArrayD();

                // tmps_[19_bb_vo](b,j) = 1.00 m1[bb_Lov](X,j,b) * s0(X) // flops: o1v1 = o1v1L1 | mem: o1v1 = o1v1
                tmps_["19_bb_vo"]("b,j")  = m1["bb_Lov"]("X,j,b") * s0("X");

                // rdm_vo_bb += +1.00 s0 m1_bb(i,a)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
                rdm_vo_bb("a,i") += tmps_["19_bb_vo"]("a,i");

                // rdm_ov_bb += -1.00 s0 t2_bbbb(a,b,j,i) m1_bb(j,b)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
                rdm_ov_bb("a,i") -= tmps_["19_bb_vo"]("b,j") * t2["bbbb_vvoo"]("a,b,j,i");
                tmps_["19_bb_vo"].~TArrayD();
            }
        }

        rdm_timer.stop();
        if (world_.rank() == 0) {
//            outfile->Printf("\n  Building EOM-EE-CC RDMs... ");
            outfile->Printf("Done --> Time: %s\n", rdm_timer.elapsed().c_str());
        }
    }

    void EOM_EE_RDM::save_density(vector<int> rdm_states) {

        // if only single state, we do not use a transition density matrix
        if (rdm_states.size() == 1) {
            rdm_states.push_back(rdm_states[0]);
        }

        // remove t1 transforms if they exist (shouldn't)
        if (eom_driver_->cc_wfn_->has_t1_integrals_)
            eom_driver_->cc_wfn_->transform_integrals(false);

        // get reference wavefunction
        const auto & cc_wfn = eom_driver_->cc_wfn_;
        size_t nso = cc_wfn->nso(); // could have done this sooner...

        // extract the RDMs
        SharedMatrix D1a(new Matrix(nso, nso));
        SharedMatrix D1b(new Matrix(nso, nso));
        
        double** D1a_p = D1a->pointer();
        double** D1b_p = D1b->pointer();
        
        // helper function to extract density from blocks
        auto extract_density = [rdm_states](TArrayD& D1_blk, double** D1_p, pair<size_t, size_t> offs){
            forall(D1_blk, [D1_p, offs, rdm_states](auto &tile, auto &x){
                size_t mu = x[2] + offs.first, nu = x[3] + offs.second;
                if (x[0] == rdm_states[0] && x[1] == rdm_states[1])
                    D1_p[mu][nu] = tile[x];
            });
        };

        extract_density(RDM_blks_["D1_aa_oo"], D1a_p, {0, 0});
        extract_density(RDM_blks_["D1_aa_ov"], D1a_p, {0, oa_});
        extract_density(RDM_blks_["D1_aa_vo"], D1a_p, {oa_, 0});
        extract_density(RDM_blks_["D1_aa_vv"], D1a_p, {oa_, oa_});

        extract_density(RDM_blks_["D1_bb_oo"], D1b_p, {0, 0});
        extract_density(RDM_blks_["D1_bb_ov"], D1b_p, {0, ob_});
        extract_density(RDM_blks_["D1_bb_vo"], D1b_p, {ob_, 0});
        extract_density(RDM_blks_["D1_bb_vv"], D1b_p, {ob_, ob_});
        world_.gop.fence();

        // transform MO basis back to AO basis
        D1a->back_transform(cc_wfn->Ca());
        D1b->back_transform(cc_wfn->Cb());

        // set objects as members of the cc_wfn
        cc_wfn->save_density(D1a, D1b);
    }

} // cc_cavity
