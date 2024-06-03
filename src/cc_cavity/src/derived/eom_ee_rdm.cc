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

        Timer rdm_timer; rdm_timer.start();

        // get cavity information
        double w0 = eom_driver_->cc_wfn_->cavity_frequency_[2];
        double coupling_factor_z = w0 * eom_driver_->cc_wfn_->cavity_coupling_strength_[2];

        TArrayMap & amplitudes = eom_driver_->cc_wfn_->amplitudes_;
        TArrayMap & evecs = eom_driver_->evec_blks_;
        TArrayMap & Id_blks = eom_driver_->cc_wfn_->Id_blks_;
        double u0 = eom_driver_->cc_wfn_->scalar_amps_["u0"];

        bool include_s0_ = include_u0_;
        bool include_m0_ = include_u0_;

        bool include_s1_ = include_u1_;
        bool include_m1_ = include_u1_;

        bool include_s2_ = include_u2_;
        bool include_m2_ = include_u2_;

        // compute 1RDM
        {

            // unpack wavefunction amplitudes
            TA::TArrayD
                    &t1_aa = amplitudes["t1_aa"], &t1_bb = amplitudes["t1_bb"], // t1 amplitudes
            &u1_aa = amplitudes["u1_aa"], &u1_bb = amplitudes["u1_bb"], // u1 amplitudes
            &t2_aaaa = amplitudes["t2_aaaa"], &t2_abab = amplitudes["t2_abab"], &t2_bbbb = amplitudes["t2_bbbb"],
                    &u2_aaaa = amplitudes["u2_aaaa"], &u2_abab = amplitudes["u2_abab"], &u2_bbbb = amplitudes["u2_bbbb"];

            // unpack eigenvectors
            TA::TArrayD
                    &r0 = evecs["r0"],
                    &r1_aa = evecs["r1_aa"], &r1_bb = evecs["r1_bb"],
                    &r2_aaaa = evecs["r2_aaaa"], &r2_abab = evecs["r2_abab"], &r2_bbbb = evecs["r2_bbbb"],
                    &l0 = evecs["l0"],
                    &l1_aa = evecs["l1_aa"], &l1_bb = evecs["l1_bb"],
                    &l2_aaaa = evecs["l2_aaaa"], &l2_abab = evecs["l2_abab"], &l2_bbbb = evecs["l2_bbbb"];

            TA::TArrayD
                    &s0 = evecs["s0"],
                    &s1_aa = evecs["s1_aa"], &s1_bb = evecs["s1_bb"],
                    &s2_aaaa = evecs["s2_aaaa"], &s2_abab = evecs["s2_abab"], &s2_bbbb = evecs["s2_bbbb"],
                    &m0 = evecs["m0"],
                    &m1_aa = evecs["m1_aa"], &m1_bb = evecs["m1_bb"],
                    &m2_aaaa = evecs["m2_aaaa"], &m2_abab = evecs["m2_abab"], &m2_bbbb = evecs["m2_bbbb"];

            {

                auto tempArray = vector<TA::TArrayD>(67);

                /// ****** p†q ****** ///

                {

                    // tempArray[0] += 0.500000 r2_bbbb("F,a,b,i,m") l2_bbbb("I,i,n,a,b")
                    // flops: o3v2L2: 1, o2v0L2: 1 | mem: o2v0L2: 2,
                    tempArray[0]("I,F,m,n") = 0.500000 * r2_bbbb("F,a,b,i,m") * l2_bbbb("I,i,n,a,b");

                    // dm_xxbb_LLoo += -0.500000 l2_bbbb(I,i,n,a,b) r2_bbbb(F,a,b,i,m)
                    // flops: o2v2L2: 1 | mem: o2v2L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempArray[0]("I,F,m,n");

                    // dm_xxbb_LLov += -0.500000 l2_bbbb(I,j,i,a,b) r2_bbbb(F,a,b,j,m) t1_bb(e,i)
                    // flops: o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[0]("I,F,m,i") * t1_bb("e,i");
                }
                tempArray[0].~TArrayD();

                if (include_u2_) {

                    // tempArray[1] += 0.500000 m2_bbbb("I,i,n,a,b") s2_bbbb("F,a,b,i,m")
                    // flops: o3v2L2: 1, o2v0L2: 1 | mem: o2v0L2: 2,
                    tempArray[1]("I,F,n,m") = 0.500000 * m2_bbbb("I,i,n,a,b") * s2_bbbb("F,a,b,i,m");

                    // dm_xxbb_LLoo += -0.500000 m2_bbbb(I,i,n,a,b) s2_bbbb(F,a,b,i,m)
                    // flops: o2v2L2: 1 | mem: o2v2L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempArray[1]("I,F,n,m");

                    // dm_xxbb_LLov += -0.500000 t1_bb(e,i) m2_bbbb(I,j,i,a,b) s2_bbbb(F,a,b,j,m)
                    // flops: o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[1]("I,F,i,m") * t1_bb("e,i");
                }
                tempArray[1].~TArrayD();

                {

                    // tempArray[2] += 1.000000 r2_abab("F,a,b,i,m") l2_abab("I,i,n,a,b")
                    // flops: o3v2L2: 1, o2v0L2: 1 | mem: o2v0L2: 2,
                    tempArray[2]("I,F,m,n") = r2_abab("F,a,b,i,m") * l2_abab("I,i,n,a,b");

                    // dm_xxbb_LLoo += -0.500000 l2_abab(I,i,n,a,b) r2_abab(F,a,b,i,m)
                    // +               -0.500000 l2_abab(I,i,n,b,a) r2_abab(F,b,a,i,m)
                    // flops: o2v2L2: 1 | mem: o2v2L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempArray[2]("I,F,m,n");

                    // dm_xxbb_LLov += -0.500000 l2_abab(I,j,i,a,b) r2_abab(F,a,b,j,m) t1_bb(e,i)
                    // +               -0.500000 l2_abab(I,j,i,b,a) r2_abab(F,b,a,j,m) t1_bb(e,i)
                    // flops: o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[2]("I,F,m,i") * t1_bb("e,i");
                }
                tempArray[2].~TArrayD();

                if (include_u2_) {

                    // tempArray[3] += 0.500000 s2_aaaa("F,a,b,i,m") m2_aaaa("I,i,n,a,b")
                    // flops: o3v2L2: 1, o2v0L2: 1 | mem: o2v0L2: 2,
                    tempArray[3]("I,F,m,n") = 0.500000 * s2_aaaa("F,a,b,i,m") * m2_aaaa("I,i,n,a,b");

                    // dm_xxaa_LLoo += -0.500000 m2_aaaa(I,i,n,a,b) s2_aaaa(F,a,b,i,m)
                    // flops: o2v2L2: 1 | mem: o2v2L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempArray[3]("I,F,m,n");

                    // dm_xxaa_LLov += -0.500000 t1_aa(e,i) m2_aaaa(I,j,i,a,b) s2_aaaa(F,a,b,j,m)
                    // flops: o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[3]("I,F,m,i") * t1_aa("e,i");
                }
                tempArray[3].~TArrayD();

                {

                    // tempArray[4] += 1.000000 r2_abab("F,a,b,m,i") l2_abab("I,n,i,a,b")
                    // flops: o3v2L2: 1, o2v0L2: 1 | mem: o2v0L2: 2,
                    tempArray[4]("I,F,m,n") = r2_abab("F,a,b,m,i") * l2_abab("I,n,i,a,b");

                    // dm_xxaa_LLoo += -0.500000 l2_abab(I,n,i,a,b) r2_abab(F,a,b,m,i)
                    // +               -0.500000 l2_abab(I,n,i,b,a) r2_abab(F,b,a,m,i)
                    // flops: o2v2L2: 1 | mem: o2v2L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempArray[4]("I,F,m,n");

                    // dm_xxaa_LLov += -0.500000 l2_abab(I,i,j,a,b) r2_abab(F,a,b,m,j) t1_aa(e,i)
                    // +               -0.500000 l2_abab(I,i,j,b,a) r2_abab(F,b,a,m,j) t1_aa(e,i)
                    // flops: o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[4]("I,F,m,i") * t1_aa("e,i");
                }
                tempArray[4].~TArrayD();

                if (include_u2_) {

                    // tempArray[5] += 1.000000 s2_abab("F,a,b,i,m") m2_abab("I,i,n,a,b")
                    // flops: o3v2L2: 1, o2v0L2: 1 | mem: o2v0L2: 2,
                    tempArray[5]("I,F,m,n") = s2_abab("F,a,b,i,m") * m2_abab("I,i,n,a,b");

                    // dm_xxbb_LLoo += -0.500000 m2_abab(I,i,n,a,b) s2_abab(F,a,b,i,m)
                    // +               -0.500000 m2_abab(I,i,n,b,a) s2_abab(F,b,a,i,m)
                    // flops: o2v2L2: 1 | mem: o2v2L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempArray[5]("I,F,m,n");

                    // dm_xxbb_LLov += -0.500000 t1_bb(e,i) m2_abab(I,j,i,a,b) s2_abab(F,a,b,j,m)
                    // +               -0.500000 t1_bb(e,i) m2_abab(I,j,i,b,a) s2_abab(F,b,a,j,m)
                    // flops: o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[5]("I,F,m,i") * t1_bb("e,i");

                    // tempArray[6] += 1.000000 s2_abab("F,b,a,m,i") m2_abab("I,n,i,b,a")
                    // flops: o3v2L2: 1, o2v0L2: 1 | mem: o2v0L2: 2,
                    tempArray[6]("I,F,m,n") = s2_abab("F,b,a,m,i") * m2_abab("I,n,i,b,a");

                    // dm_xxaa_LLoo += -0.500000 m2_abab(I,n,i,a,b) s2_abab(F,a,b,m,i)
                    // +               -0.500000 m2_abab(I,n,i,b,a) s2_abab(F,b,a,m,i)
                    // flops: o2v2L2: 1 | mem: o2v2L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempArray[6]("I,F,m,n");

                    // dm_xxaa_LLov += -0.500000 t1_aa(e,i) m2_abab(I,i,j,a,b) s2_abab(F,a,b,m,j)
                    // +               -0.500000 t1_aa(e,i) m2_abab(I,i,j,b,a) s2_abab(F,b,a,m,j)
                    // flops: o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[6]("I,F,m,i") * t1_aa("e,i");
                }
                tempArray[5].~TArrayD();
                tempArray[6].~TArrayD();

                {

                    // tempArray[7] += 0.500000 r2_aaaa("F,a,b,i,m") l2_aaaa("I,i,n,a,b")
                    // flops: o3v2L2: 1, o2v0L2: 1 | mem: o2v0L2: 2,
                    tempArray[7]("I,F,m,n") = 0.500000 * r2_aaaa("F,a,b,i,m") * l2_aaaa("I,i,n,a,b");

                    // dm_xxaa_LLoo += -0.500000 l2_aaaa(I,i,n,a,b) r2_aaaa(F,a,b,i,m)
                    // flops: o2v2L2: 1 | mem: o2v2L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempArray[7]("I,F,m,n");

                    // dm_xxaa_LLov += -0.500000 l2_aaaa(I,j,i,a,b) r2_aaaa(F,a,b,j,m) t1_aa(e,i)
                    // flops: o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[7]("I,F,m,i") * t1_aa("e,i");
                }
                tempArray[7].~TArrayD();

                if (include_u2_) {

                    // tempArray[8] += 1.000000 s2_abab("F,a,b,j,i") m2_abab("I,j,i,a,b")
                    // flops: o2v2L2: 1, o0v0L2: 1 | mem: o0v0L2: 2,
                    tempArray[8]("I,F") = s2_abab("F,a,b,j,i") * m2_abab("I,j,i,a,b");

                    // dm_xxaa_LLoo += 0.250000 d_aa(m,n) m2_abab(I,j,i,a,b) s2_abab(F,a,b,j,i)
                    // +               0.250000 d_aa(m,n) m2_abab(I,j,i,b,a) s2_abab(F,b,a,j,i)
                    // +               0.250000 d_aa(m,n) m2_abab(I,i,j,a,b) s2_abab(F,a,b,i,j)
                    // +               0.250000 d_aa(m,n) m2_abab(I,i,j,b,a) s2_abab(F,b,a,i,j)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempArray[8]("I,F") * Id_blks["aa_oo"]("m,n");

                    // dm_xxbb_LLoo += 0.250000 d_bb(m,n) m2_abab(I,j,i,a,b) s2_abab(F,a,b,j,i)
                    // +               0.250000 d_bb(m,n) m2_abab(I,j,i,b,a) s2_abab(F,b,a,j,i)
                    // +               0.250000 d_bb(m,n) m2_abab(I,i,j,a,b) s2_abab(F,a,b,i,j)
                    // +               0.250000 d_bb(m,n) m2_abab(I,i,j,b,a) s2_abab(F,b,a,i,j)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempArray[8]("I,F") * Id_blks["bb_oo"]("m,n");

                    // dm_xxaa_LLov += 0.250000 t1_aa(e,m) m2_abab(I,j,i,a,b) s2_abab(F,a,b,j,i)
                    // +               0.250000 t1_aa(e,m) m2_abab(I,j,i,b,a) s2_abab(F,b,a,j,i)
                    // +               0.250000 t1_aa(e,m) m2_abab(I,i,j,a,b) s2_abab(F,a,b,i,j)
                    // +               0.250000 t1_aa(e,m) m2_abab(I,i,j,b,a) s2_abab(F,b,a,i,j)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[8]("I,F") * t1_aa("e,m");

                    // dm_xxbb_LLov += 0.250000 t1_bb(e,m) m2_abab(I,j,i,a,b) s2_abab(F,a,b,j,i)
                    // +               0.250000 t1_bb(e,m) m2_abab(I,j,i,b,a) s2_abab(F,b,a,j,i)
                    // +               0.250000 t1_bb(e,m) m2_abab(I,i,j,a,b) s2_abab(F,a,b,i,j)
                    // +               0.250000 t1_bb(e,m) m2_abab(I,i,j,b,a) s2_abab(F,b,a,i,j)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[8]("I,F") * t1_bb("e,m");
                }
                tempArray[8].~TArrayD();

                {

                    // tempArray[9] += 0.250000 r2_aaaa("F,a,b,j,i") l2_aaaa("I,j,i,a,b")
                    // flops: o2v2L2: 1, o0v0L2: 1 | mem: o0v0L2: 2,
                    tempArray[9]("I,F") = 0.250000 * r2_aaaa("F,a,b,j,i") * l2_aaaa("I,j,i,a,b");

                    // dm_xxaa_LLoo += 0.250000 d_aa(m,n) l2_aaaa(I,j,i,a,b) r2_aaaa(F,a,b,j,i)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempArray[9]("I,F") * Id_blks["aa_oo"]("m,n");

                    // dm_xxbb_LLoo += 0.250000 d_bb(m,n) l2_aaaa(I,j,i,a,b) r2_aaaa(F,a,b,j,i)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempArray[9]("I,F") * Id_blks["bb_oo"]("m,n");

                    // dm_xxaa_LLov += 0.250000 l2_aaaa(I,j,i,a,b) r2_aaaa(F,a,b,j,i) t1_aa(e,m)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[9]("I,F") * t1_aa("e,m");

                    // dm_xxbb_LLov += 0.250000 l2_aaaa(I,j,i,a,b) r2_aaaa(F,a,b,j,i) t1_bb(e,m)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[9]("I,F") * t1_bb("e,m");
                }
                tempArray[9].~TArrayD();

                if (include_u2_) {

                    // tempArray[10] += 0.250000 s2_aaaa("F,a,b,j,i") m2_aaaa("I,j,i,a,b")
                    // flops: o2v2L2: 1, o0v0L2: 1 | mem: o0v0L2: 2,
                    tempArray[10]("I,F") = 0.250000 * s2_aaaa("F,a,b,j,i") * m2_aaaa("I,j,i,a,b");

                    // dm_xxaa_LLoo += 0.250000 d_aa(m,n) m2_aaaa(I,j,i,a,b) s2_aaaa(F,a,b,j,i)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempArray[10]("I,F") * Id_blks["aa_oo"]("m,n");

                    // dm_xxbb_LLoo += 0.250000 d_bb(m,n) m2_aaaa(I,j,i,a,b) s2_aaaa(F,a,b,j,i)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempArray[10]("I,F") * Id_blks["bb_oo"]("m,n");

                    // dm_xxaa_LLov += 0.250000 t1_aa(e,m) m2_aaaa(I,j,i,a,b) s2_aaaa(F,a,b,j,i)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[10]("I,F") * t1_aa("e,m");

                    // dm_xxbb_LLov += 0.250000 t1_bb(e,m) m2_aaaa(I,j,i,a,b) s2_aaaa(F,a,b,j,i)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[10]("I,F") * t1_bb("e,m");

                    // tempArray[11] += 0.250000 s2_bbbb("F,a,b,j,i") m2_bbbb("I,j,i,a,b")
                    // flops: o2v2L2: 1, o0v0L2: 1 | mem: o0v0L2: 2,
                    tempArray[11]("I,F") = 0.250000 * s2_bbbb("F,a,b,j,i") * m2_bbbb("I,j,i,a,b");

                    // dm_xxaa_LLoo += 0.250000 d_aa(m,n) m2_bbbb(I,j,i,a,b) s2_bbbb(F,a,b,j,i)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempArray[11]("I,F") * Id_blks["aa_oo"]("m,n");

                    // dm_xxbb_LLoo += 0.250000 d_bb(m,n) m2_bbbb(I,j,i,a,b) s2_bbbb(F,a,b,j,i)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempArray[11]("I,F") * Id_blks["bb_oo"]("m,n");

                    // dm_xxaa_LLov += 0.250000 t1_aa(e,m) m2_bbbb(I,j,i,a,b) s2_bbbb(F,a,b,j,i)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[11]("I,F") * t1_aa("e,m");

                    // dm_xxbb_LLov += 0.250000 t1_bb(e,m) m2_bbbb(I,j,i,a,b) s2_bbbb(F,a,b,j,i)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[11]("I,F") * t1_bb("e,m");
                }
                tempArray[10].~TArrayD();
                tempArray[11].~TArrayD();

                {

                    // tempArray[12] += 1.000000 l2_abab("I,j,i,b,a") r2_abab("F,b,a,j,i")
                    // flops: o2v2L2: 1, o0v0L2: 1 | mem: o0v0L2: 2,
                    tempArray[12]("I,F") = l2_abab("I,j,i,b,a") * r2_abab("F,b,a,j,i");

                    // dm_xxaa_LLoo += 0.250000 d_aa(m,n) l2_abab(I,j,i,a,b) r2_abab(F,a,b,j,i)
                    // +               0.250000 d_aa(m,n) l2_abab(I,j,i,b,a) r2_abab(F,b,a,j,i)
                    // +               0.250000 d_aa(m,n) l2_abab(I,i,j,a,b) r2_abab(F,a,b,i,j)
                    // +               0.250000 d_aa(m,n) l2_abab(I,i,j,b,a) r2_abab(F,b,a,i,j)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempArray[12]("I,F") * Id_blks["aa_oo"]("m,n");

                    // dm_xxbb_LLoo += 0.250000 d_bb(m,n) l2_abab(I,j,i,a,b) r2_abab(F,a,b,j,i)
                    // +               0.250000 d_bb(m,n) l2_abab(I,j,i,b,a) r2_abab(F,b,a,j,i)
                    // +               0.250000 d_bb(m,n) l2_abab(I,i,j,a,b) r2_abab(F,a,b,i,j)
                    // +               0.250000 d_bb(m,n) l2_abab(I,i,j,b,a) r2_abab(F,b,a,i,j)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempArray[12]("I,F") * Id_blks["bb_oo"]("m,n");

                    // dm_xxaa_LLov += 0.250000 l2_abab(I,j,i,a,b) r2_abab(F,a,b,j,i) t1_aa(e,m)
                    // +               0.250000 l2_abab(I,j,i,b,a) r2_abab(F,b,a,j,i) t1_aa(e,m)
                    // +               0.250000 l2_abab(I,i,j,a,b) r2_abab(F,a,b,i,j) t1_aa(e,m)
                    // +               0.250000 l2_abab(I,i,j,b,a) r2_abab(F,b,a,i,j) t1_aa(e,m)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[12]("I,F") * t1_aa("e,m");

                    // dm_xxbb_LLov += 0.250000 l2_abab(I,j,i,a,b) r2_abab(F,a,b,j,i) t1_bb(e,m)
                    // +               0.250000 l2_abab(I,j,i,b,a) r2_abab(F,b,a,j,i) t1_bb(e,m)
                    // +               0.250000 l2_abab(I,i,j,a,b) r2_abab(F,a,b,i,j) t1_bb(e,m)
                    // +               0.250000 l2_abab(I,i,j,b,a) r2_abab(F,b,a,i,j) t1_bb(e,m)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[12]("I,F") * t1_bb("e,m");

                    // tempArray[13] += 0.250000 r2_bbbb("F,a,b,j,i") l2_bbbb("I,j,i,a,b")
                    // flops: o2v2L2: 1, o0v0L2: 1 | mem: o0v0L2: 2,
                    tempArray[13]("I,F") = 0.250000 * r2_bbbb("F,a,b,j,i") * l2_bbbb("I,j,i,a,b");

                    // dm_xxaa_LLoo += 0.250000 d_aa(m,n) l2_bbbb(I,j,i,a,b) r2_bbbb(F,a,b,j,i)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempArray[13]("I,F") * Id_blks["aa_oo"]("m,n");

                    // dm_xxbb_LLoo += 0.250000 d_bb(m,n) l2_bbbb(I,j,i,a,b) r2_bbbb(F,a,b,j,i)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempArray[13]("I,F") * Id_blks["bb_oo"]("m,n");

                    // dm_xxaa_LLov += 0.250000 l2_bbbb(I,j,i,a,b) r2_bbbb(F,a,b,j,i) t1_aa(e,m)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[13]("I,F") * t1_aa("e,m");

                    // dm_xxbb_LLov += 0.250000 l2_bbbb(I,j,i,a,b) r2_bbbb(F,a,b,j,i) t1_bb(e,m)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[13]("I,F") * t1_bb("e,m");
                }
                tempArray[12].~TArrayD();
                tempArray[13].~TArrayD();

                if (include_u2_) {

                    // tempArray[14] += 1.000000 r1_bb("F,a,i") m2_bbbb("I,i,j,e,a")
                    // flops: o2v2L2: 1, o1v1L2: 1 | mem: o1v1L2: 2,
                    tempArray[14]("I,F,e,j") = r1_bb("F,a,i") * m2_bbbb("I,i,j,e,a");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxbb_LLvv += -1.000000 r1_bb(F,a,i) u1_bb(f,j) m2_bbbb(I,i,j,e,a)
                    // flops: o0v4L2: 1, o1v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") -= tempArray[14]("I,F,e,j") * u1_bb("f,j");
                }

                if (include_u2_) {

                    // dm_xxaa_LLov += -1.000000 r1_bb(F,a,i) u2_abab(e,b,m,j) m2_bbbb(I,i,j,b,a)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[14]("I,F,b,j") * u2_abab("e,b,m,j");

                    // dm_xxbb_LLov += 1.000000 r1_bb(F,a,i) u2_bbbb(e,b,j,m) m2_bbbb(I,i,j,b,a)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[14]("I,F,b,j") * u2_bbbb("e,b,j,m");

                    // tempArray[15] += 1.000000 m2_aaaa("I,i,j,e,a") r1_aa("F,a,i")
                    // flops: o2v2L2: 1, o1v1L2: 1 | mem: o1v1L2: 2,
                    tempArray[15]("I,F,e,j") = m2_aaaa("I,i,j,e,a") * r1_aa("F,a,i");
                }
                tempArray[14].~TArrayD();

                if (include_u1_ && include_u2_) {

                    // dm_xxaa_LLvv += -1.000000 r1_aa(F,a,i) u1_aa(f,j) m2_aaaa(I,i,j,e,a)
                    // flops: o0v4L2: 1, o1v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") -= tempArray[15]("I,F,e,j") * u1_aa("f,j");
                }

                if (include_u2_) {

                    // dm_xxaa_LLov += 1.000000 r1_aa(F,a,i) u2_aaaa(e,b,j,m) m2_aaaa(I,i,j,b,a)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[15]("I,F,b,j") * u2_aaaa("e,b,j,m");

                    // dm_xxbb_LLov += -1.000000 r1_aa(F,a,i) u2_abab(b,e,j,m) m2_aaaa(I,i,j,b,a)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[15]("I,F,b,j") * u2_abab("b,e,j,m");
                }
                tempArray[15].~TArrayD();

                if (include_u1_ && include_u2_) {

                    // tempArray[16] += 1.000000 s1_bb("F,a,j") m2_abab("I,i,j,e,a")
                    // flops: o2v2L2: 1, o1v1L2: 1 | mem: o1v1L2: 2,
                    tempArray[16]("I,F,e,i") = s1_bb("F,a,j") * m2_abab("I,i,j,e,a");

                    // dm_xxaa_LLvv += 1.000000 t1_aa(f,i) m2_abab(I,i,j,e,a) s1_bb(F,a,j)
                    // flops: o0v4L2: 1, o1v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += tempArray[16]("I,F,e,i") * t1_aa("f,i");

                    // dm_xxaa_LLov += -1.000000 t2_aaaa(e,a,i,m) m2_abab(I,i,j,a,b) s1_bb(F,b,j)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[16]("I,F,a,i") * t2_aaaa("e,a,i,m");

                    // dm_xxbb_LLov += 1.000000 t2_abab(a,e,i,m) m2_abab(I,i,j,a,b) s1_bb(F,b,j)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[16]("I,F,a,i") * t2_abab("a,e,i,m");

                    // dm_xxaa_LLvo += 1.000000 m2_abab(I,m,i,e,a) s1_bb(F,a,i)
                    // flops: o1v3L2: 1 | mem: o1v3L2: 1,
                    RDM_blks_["D1_aa_vo"]("I,F,e,m") += tempArray[16]("I,F,e,m");
                }
                tempArray[16].~TArrayD();

                {

                    // tempArray[17] += 1.000000 l2_abab("I,i,j,e,a") r1_bb("F,a,j")
                    // flops: o2v2L2: 1, o1v1L2: 1 | mem: o1v1L2: 2,
                    tempArray[17]("I,F,e,i") = l2_abab("I,i,j,e,a") * r1_bb("F,a,j");

                    // dm_xxaa_LLvv += 1.000000 l2_abab(I,i,j,e,a) r1_bb(F,a,j) t1_aa(f,i)
                    // flops: o0v4L2: 1, o1v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += tempArray[17]("I,F,e,i") * t1_aa("f,i");

                    // dm_xxaa_LLov += -1.000000 l2_abab(I,i,j,a,b) r1_bb(F,b,j) t2_aaaa(e,a,i,m)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[17]("I,F,a,i") * t2_aaaa("e,a,i,m");

                    // dm_xxbb_LLov += 1.000000 l2_abab(I,i,j,a,b) r1_bb(F,b,j) t2_abab(a,e,i,m)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[17]("I,F,a,i") * t2_abab("a,e,i,m");

                    // dm_xxaa_LLvo += 1.000000 l2_abab(I,m,i,e,a) r1_bb(F,a,i)
                    // flops: o1v3L2: 1 | mem: o1v3L2: 1,
                    RDM_blks_["D1_aa_vo"]("I,F,e,m") += tempArray[17]("I,F,e,m");
                }
                tempArray[17].~TArrayD();

                if (include_u1_ && include_u2_) {

                    // tempArray[18] += 1.000000 s1_aa("F,a,j") m2_aaaa("I,j,i,e,a")
                    // flops: o2v2L2: 1, o1v1L2: 1 | mem: o1v1L2: 2,
                    tempArray[18]("I,F,e,i") = s1_aa("F,a,j") * m2_aaaa("I,j,i,e,a");

                    // dm_xxaa_LLvv += -1.000000 t1_aa(f,i) m2_aaaa(I,j,i,e,a) s1_aa(F,a,j)
                    // flops: o0v4L2: 1, o1v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") -= tempArray[18]("I,F,e,i") * t1_aa("f,i");

                    // dm_xxaa_LLov += 1.000000 t2_aaaa(e,a,i,m) m2_aaaa(I,j,i,a,b) s1_aa(F,b,j)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[18]("I,F,a,i") * t2_aaaa("e,a,i,m");

                    // dm_xxbb_LLov += -1.000000 t2_abab(a,e,i,m) m2_aaaa(I,j,i,a,b) s1_aa(F,b,j)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[18]("I,F,a,i") * t2_abab("a,e,i,m");

                    // dm_xxaa_LLvo += -1.000000 m2_aaaa(I,i,m,e,a) s1_aa(F,a,i)
                    // flops: o1v3L2: 1 | mem: o1v3L2: 1,
                    RDM_blks_["D1_aa_vo"]("I,F,e,m") -= tempArray[18]("I,F,e,m");
                }
                tempArray[18].~TArrayD();

                if (include_u2_) {

                    // tempArray[19] += 1.000000 m2_abab("I,i,j,a,e") r1_aa("F,a,i")
                    // flops: o2v2L2: 1, o1v1L2: 1 | mem: o1v1L2: 2,
                    tempArray[19]("I,F,e,j") = m2_abab("I,i,j,a,e") * r1_aa("F,a,i");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxbb_LLvv += 1.000000 r1_aa(F,a,i) u1_bb(f,j) m2_abab(I,i,j,a,e)
                    // flops: o0v4L2: 1, o1v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += tempArray[19]("I,F,e,j") * u1_bb("f,j");
                }

                if (include_u2_) {

                    // dm_xxaa_LLov += 1.000000 r1_aa(F,a,i) u2_abab(e,b,m,j) m2_abab(I,i,j,a,b)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[19]("I,F,b,j") * u2_abab("e,b,m,j");

                    // dm_xxbb_LLov += -1.000000 r1_aa(F,a,i) u2_bbbb(e,b,j,m) m2_abab(I,i,j,a,b)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[19]("I,F,b,j") * u2_bbbb("e,b,j,m");
                }
                tempArray[19].~TArrayD();

                {

                    // tempArray[20] += 1.000000 r1_bb("F,a,j") l2_bbbb("I,j,i,e,a")
                    // flops: o2v2L2: 1, o1v1L2: 1 | mem: o1v1L2: 2,
                    tempArray[20]("I,F,e,i") = r1_bb("F,a,j") * l2_bbbb("I,j,i,e,a");

                    // dm_xxbb_LLvv += -1.000000 l2_bbbb(I,j,i,e,a) r1_bb(F,a,j) t1_bb(f,i)
                    // flops: o0v4L2: 1, o1v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") -= tempArray[20]("I,F,e,i") * t1_bb("f,i");

                    // dm_xxaa_LLov += -1.000000 l2_bbbb(I,j,i,a,b) r1_bb(F,b,j) t2_abab(e,a,m,i)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[20]("I,F,a,i") * t2_abab("e,a,m,i");

                    // dm_xxbb_LLov += 1.000000 l2_bbbb(I,j,i,a,b) r1_bb(F,b,j) t2_bbbb(e,a,i,m)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[20]("I,F,a,i") * t2_bbbb("e,a,i,m");

                    // dm_xxbb_LLvo += -1.000000 l2_bbbb(I,i,m,e,a) r1_bb(F,a,i)
                    // flops: o1v3L2: 1 | mem: o1v3L2: 1,
                    RDM_blks_["D1_bb_vo"]("I,F,e,m") -= tempArray[20]("I,F,e,m");

                    // tempArray[21] += 1.000000 r1_aa("F,a,j") l2_aaaa("I,j,i,e,a")
                    // flops: o2v2L2: 1, o1v1L2: 1 | mem: o1v1L2: 2,
                    tempArray[21]("I,F,e,i") = r1_aa("F,a,j") * l2_aaaa("I,j,i,e,a");

                    // dm_xxaa_LLvv += -1.000000 l2_aaaa(I,j,i,e,a) r1_aa(F,a,j) t1_aa(f,i)
                    // flops: o0v4L2: 1, o1v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") -= tempArray[21]("I,F,e,i") * t1_aa("f,i");

                    // dm_xxaa_LLov += 1.000000 l2_aaaa(I,j,i,a,b) r1_aa(F,b,j) t2_aaaa(e,a,i,m)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[21]("I,F,a,i") * t2_aaaa("e,a,i,m");

                    // dm_xxbb_LLov += -1.000000 l2_aaaa(I,j,i,a,b) r1_aa(F,b,j) t2_abab(a,e,i,m)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[21]("I,F,a,i") * t2_abab("a,e,i,m");

                    // dm_xxaa_LLvo += -1.000000 l2_aaaa(I,i,m,e,a) r1_aa(F,a,i)
                    // flops: o1v3L2: 1 | mem: o1v3L2: 1,
                    RDM_blks_["D1_aa_vo"]("I,F,e,m") -= tempArray[21]("I,F,e,m");
                }
                tempArray[20].~TArrayD();
                tempArray[21].~TArrayD();

                if (include_u1_ && include_u2_) {

                    // tempArray[22] += 1.000000 s1_bb("F,a,j") m2_bbbb("I,j,i,e,a")
                    // flops: o2v2L2: 1, o1v1L2: 1 | mem: o1v1L2: 2,
                    tempArray[22]("I,F,e,i") = s1_bb("F,a,j") * m2_bbbb("I,j,i,e,a");

                    // dm_xxbb_LLvv += -1.000000 t1_bb(f,i) m2_bbbb(I,j,i,e,a) s1_bb(F,a,j)
                    // flops: o0v4L2: 1, o1v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") -= tempArray[22]("I,F,e,i") * t1_bb("f,i");

                    // dm_xxaa_LLov += -1.000000 t2_abab(e,a,m,i) m2_bbbb(I,j,i,a,b) s1_bb(F,b,j)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[22]("I,F,a,i") * t2_abab("e,a,m,i");

                    // dm_xxbb_LLov += 1.000000 t2_bbbb(e,a,i,m) m2_bbbb(I,j,i,a,b) s1_bb(F,b,j)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[22]("I,F,a,i") * t2_bbbb("e,a,i,m");

                    // dm_xxbb_LLvo += -1.000000 m2_bbbb(I,i,m,e,a) s1_bb(F,a,i)
                    // flops: o1v3L2: 1 | mem: o1v3L2: 1,
                    RDM_blks_["D1_bb_vo"]("I,F,e,m") -= tempArray[22]("I,F,e,m");
                }
                tempArray[22].~TArrayD();

                if (include_u2_) {

                    // tempArray[23] += 1.000000 m2_abab("I,j,i,e,a") r1_bb("F,a,i")
                    // flops: o2v2L2: 1, o1v1L2: 1 | mem: o1v1L2: 2,
                    tempArray[23]("I,F,e,j") = m2_abab("I,j,i,e,a") * r1_bb("F,a,i");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxaa_LLvv += 1.000000 r1_bb(F,a,i) u1_aa(f,j) m2_abab(I,j,i,e,a)
                    // flops: o0v4L2: 1, o1v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += tempArray[23]("I,F,e,j") * u1_aa("f,j");
                }

                if (include_u2_) {

                    // dm_xxaa_LLov += -1.000000 r1_bb(F,a,i) u2_aaaa(e,b,j,m) m2_abab(I,j,i,b,a)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[23]("I,F,b,j") * u2_aaaa("e,b,j,m");

                    // dm_xxbb_LLov += 1.000000 r1_bb(F,a,i) u2_abab(b,e,j,m) m2_abab(I,j,i,b,a)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[23]("I,F,b,j") * u2_abab("b,e,j,m");
                }
                tempArray[23].~TArrayD();

                {

                    // tempArray[24] += 1.000000 r1_aa("F,a,j") l2_abab("I,j,i,a,e")
                    // flops: o2v2L2: 1, o1v1L2: 1 | mem: o1v1L2: 2,
                    tempArray[24]("I,F,e,i") = r1_aa("F,a,j") * l2_abab("I,j,i,a,e");

                    // dm_xxbb_LLvv += 1.000000 l2_abab(I,j,i,a,e) r1_aa(F,a,j) t1_bb(f,i)
                    // flops: o0v4L2: 1, o1v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += tempArray[24]("I,F,e,i") * t1_bb("f,i");

                    // dm_xxaa_LLov += 1.000000 l2_abab(I,j,i,b,a) r1_aa(F,b,j) t2_abab(e,a,m,i)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[24]("I,F,a,i") * t2_abab("e,a,m,i");

                    // dm_xxbb_LLov += -1.000000 l2_abab(I,j,i,b,a) r1_aa(F,b,j) t2_bbbb(e,a,i,m)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[24]("I,F,a,i") * t2_bbbb("e,a,i,m");

                    // dm_xxbb_LLvo += 1.000000 l2_abab(I,i,m,a,e) r1_aa(F,a,i)
                    // flops: o1v3L2: 1 | mem: o1v3L2: 1,
                    RDM_blks_["D1_bb_vo"]("I,F,e,m") += tempArray[24]("I,F,e,m");
                }
                tempArray[24].~TArrayD();

                if (include_u1_ && include_u2_) {

                    // tempArray[25] += 1.000000 m2_abab("I,j,i,a,e") s1_aa("F,a,j")
                    // flops: o2v2L2: 1, o1v1L2: 1 | mem: o1v1L2: 2,
                    tempArray[25]("I,F,e,i") = m2_abab("I,j,i,a,e") * s1_aa("F,a,j");

                    // dm_xxbb_LLvv += 1.000000 t1_bb(f,i) m2_abab(I,j,i,a,e) s1_aa(F,a,j)
                    // flops: o0v4L2: 1, o1v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += tempArray[25]("I,F,e,i") * t1_bb("f,i");

                    // dm_xxaa_LLov += 1.000000 t2_abab(e,a,m,i) m2_abab(I,j,i,b,a) s1_aa(F,b,j)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[25]("I,F,a,i") * t2_abab("e,a,m,i");

                    // dm_xxbb_LLov += -1.000000 t2_bbbb(e,a,i,m) m2_abab(I,j,i,b,a) s1_aa(F,b,j)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[25]("I,F,a,i") * t2_bbbb("e,a,i,m");

                    // dm_xxbb_LLvo += 1.000000 m2_abab(I,i,m,a,e) s1_aa(F,a,i)
                    // flops: o1v3L2: 1 | mem: o1v3L2: 1,
                    RDM_blks_["D1_bb_vo"]("I,F,e,m") += tempArray[25]("I,F,e,m");
                }
                tempArray[25].~TArrayD();

                if (include_u2_) {

                    // tempArray[26] += 0.250000 m2_bbbb("I,i,j,b,a") r2_bbbb("F,b,a,i,j")
                    // flops: o2v2L2: 1, o0v0L2: 1 | mem: o0v0L2: 2,
                    tempArray[26]("I,F") = 0.250000 * m2_bbbb("I,i,j,b,a") * r2_bbbb("F,b,a,i,j");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxaa_LLov += 0.250000 r2_bbbb(F,b,a,i,j) u1_aa(e,m) m2_bbbb(I,i,j,b,a)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[26]("I,F") * u1_aa("e,m");

                    // dm_xxbb_LLov += 0.250000 r2_bbbb(F,b,a,i,j) u1_bb(e,m) m2_bbbb(I,i,j,b,a)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[26]("I,F") * u1_bb("e,m");
                }
                tempArray[26].~TArrayD();

                if (include_u2_) {

                    // tempArray[27] += 0.250000 m2_aaaa("I,i,j,b,a") r2_aaaa("F,b,a,i,j")
                    // flops: o2v2L2: 1, o0v0L2: 1 | mem: o0v0L2: 2,
                    tempArray[27]("I,F") = 0.250000 * m2_aaaa("I,i,j,b,a") * r2_aaaa("F,b,a,i,j");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxaa_LLov += 0.250000 r2_aaaa(F,b,a,i,j) u1_aa(e,m) m2_aaaa(I,i,j,b,a)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[27]("I,F") * u1_aa("e,m");

                    // dm_xxbb_LLov += 0.250000 r2_aaaa(F,b,a,i,j) u1_bb(e,m) m2_aaaa(I,i,j,b,a)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[27]("I,F") * u1_bb("e,m");
                }
                tempArray[27].~TArrayD();

                if (include_u2_) {

                    // tempArray[28] += 1.000000 r2_abab("F,b,a,i,j") m2_abab("I,i,j,b,a")
                    // flops: o2v2L2: 1, o0v0L2: 1 | mem: o0v0L2: 2,
                    tempArray[28]("I,F") = r2_abab("F,b,a,i,j") * m2_abab("I,i,j,b,a");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxaa_LLov += 0.250000 r2_abab(F,b,a,i,j) u1_aa(e,m) m2_abab(I,i,j,b,a)
                    // +               0.250000 r2_abab(F,b,a,j,i) u1_aa(e,m) m2_abab(I,j,i,b,a)
                    // +               0.250000 r2_abab(F,a,b,i,j) u1_aa(e,m) m2_abab(I,i,j,a,b)
                    // +               0.250000 r2_abab(F,a,b,j,i) u1_aa(e,m) m2_abab(I,j,i,a,b)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[28]("I,F") * u1_aa("e,m");

                    // dm_xxbb_LLov += 0.250000 r2_abab(F,b,a,i,j) u1_bb(e,m) m2_abab(I,i,j,b,a)
                    // +               0.250000 r2_abab(F,b,a,j,i) u1_bb(e,m) m2_abab(I,j,i,b,a)
                    // +               0.250000 r2_abab(F,a,b,i,j) u1_bb(e,m) m2_abab(I,i,j,a,b)
                    // +               0.250000 r2_abab(F,a,b,j,i) u1_bb(e,m) m2_abab(I,j,i,a,b)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[28]("I,F") * u1_bb("e,m");
                }
                tempArray[28].~TArrayD();

                if (include_u2_) {

                    // tempArray[29] += 1.000000 m2_abab("I,j,i,a,e") t2_abab("a,f,j,i")
                    // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
                    tempArray[29]("I,e,f") = m2_abab("I,j,i,a,e") * t2_abab("a,f,j,i");
                }

                if (include_u0_ && include_u2_) {

                    // dm_xxbb_LLvv += 0.500000 t2_abab(a,f,i,j) m2_abab(I,i,j,a,e) s0(F)
                    // +               0.500000 t2_abab(a,f,j,i) m2_abab(I,j,i,a,e) s0(F)
                    // flops: o0v4L2: 1, o0v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += tempArray[29]("I,e,f") * s0("F");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxbb_LLov += -0.500000 t2_abab(a,e,i,j) m2_abab(I,i,j,a,b) s1_bb(F,b,m)
                    // +               -0.500000 t2_abab(a,e,j,i) m2_abab(I,j,i,a,b) s1_bb(F,b,m)
                    // flops: o1v3L2: 1, o1v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[29]("I,b,e") * s1_bb("F,b,m");
                }
                tempArray[29].~TArrayD();

                {

                    // tempArray[30] += 1.000000 t2_abab("a,f,i,j") l2_abab("I,i,j,a,e")
                    // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
                    tempArray[30]("I,f,e") = t2_abab("a,f,i,j") * l2_abab("I,i,j,a,e");

                    // dm_xxbb_LLvv += 0.500000 l2_abab(I,i,j,a,e) r0(F) t2_abab(a,f,i,j)
                    // +               0.500000 l2_abab(I,j,i,a,e) r0(F) t2_abab(a,f,j,i)
                    // flops: o0v4L2: 1, o0v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += tempArray[30]("I,f,e") * r0("F");

                    // dm_xxbb_LLov += -0.500000 l2_abab(I,i,j,a,b) r1_bb(F,b,m) t2_abab(a,e,i,j)
                    // +               -0.500000 l2_abab(I,j,i,a,b) r1_bb(F,b,m) t2_abab(a,e,j,i)
                    // flops: o1v3L2: 1, o1v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[30]("I,e,b") * r1_bb("F,b,m");
                }
                tempArray[30].~TArrayD();

                if (include_u2_) {

                    // tempArray[31] += 1.000000 u2_abab("f,a,i,j") m2_abab("I,i,j,e,a")
                    // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
                    tempArray[31]("I,f,e") = u2_abab("f,a,i,j") * m2_abab("I,i,j,e,a");

                    // RDM_blks_["D1_aa_vv"](F) u2_abab(f,a,i,j) m2_abab(I,i,j,e,a)
                    // +               0.500000 r0(F) u2_abab(f,a,j,i) m2_abab(I,j,i,e,a)
                    // flops: o0v4L2: 1, o0v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += tempArray[31]("I,f,e") * r0("F");

                    // dm_xxaa_LLov += -0.500000 r1_aa(F,a,m) u2_abab(e,b,i,j) m2_abab(I,i,j,a,b)
                    // +               -0.500000 r1_aa(F,a,m) u2_abab(e,b,j,i) m2_abab(I,j,i,a,b)
                    // flops: o1v3L2: 1, o1v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[31]("I,e,a") * r1_aa("F,a,m");
                }
                tempArray[31].~TArrayD();

                {

                    // tempArray[32] += 1.000000 l2_abab("I,i,j,e,a") t2_abab("f,a,i,j")
                    // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
                    tempArray[32]("I,e,f") = l2_abab("I,i,j,e,a") * t2_abab("f,a,i,j");

                    // dm_xxaa_LLvv += 0.500000 l2_abab(I,i,j,e,a) r0(F) t2_abab(f,a,i,j)
                    // +               0.500000 l2_abab(I,j,i,e,a) r0(F) t2_abab(f,a,j,i)
                    // flops: o0v4L2: 1, o0v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += tempArray[32]("I,e,f") * r0("F");

                    // dm_xxaa_LLov += -0.500000 l2_abab(I,i,j,b,a) r1_aa(F,b,m) t2_abab(e,a,i,j)
                    // +               -0.500000 l2_abab(I,j,i,b,a) r1_aa(F,b,m) t2_abab(e,a,j,i)
                    // flops: o1v3L2: 1, o1v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[32]("I,b,e") * r1_aa("F,b,m");
                }
                tempArray[32].~TArrayD();

                if (include_u2_) {

                    // tempArray[33] += 1.000000 u2_abab("a,f,j,i") m2_abab("I,j,i,a,e")
                    // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
                    tempArray[33]("I,f,e") = u2_abab("a,f,j,i") * m2_abab("I,j,i,a,e");

                    // RDM_blks_["D1_bb_vv"](F) u2_abab(a,f,i,j) m2_abab(I,i,j,a,e)
                    // +               0.500000 r0(F) u2_abab(a,f,j,i) m2_abab(I,j,i,a,e)
                    // flops: o0v4L2: 1, o0v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += tempArray[33]("I,f,e") * r0("F");

                    // dm_xxbb_LLov += -0.500000 r1_bb(F,a,m) u2_abab(b,e,i,j) m2_abab(I,i,j,b,a)
                    // +               -0.500000 r1_bb(F,a,m) u2_abab(b,e,j,i) m2_abab(I,j,i,b,a)
                    // flops: o1v3L2: 1, o1v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[33]("I,e,a") * r1_bb("F,a,m");

                    // tempArray[34] += 1.000000 m2_abab("I,i,j,e,a") t2_abab("f,a,i,j")
                    // flops: o2v3L1: 1, o0v2L1: 1 | mem: o0v2L1: 2,
                    tempArray[34]("I,e,f") = m2_abab("I,i,j,e,a") * t2_abab("f,a,i,j");
                }
                tempArray[33].~TArrayD();

                if (include_u0_ && include_u2_) {

                    // dm_xxaa_LLvv += 0.500000 t2_abab(f,a,i,j) m2_abab(I,i,j,e,a) s0(F)
                    // +               0.500000 t2_abab(f,a,j,i) m2_abab(I,j,i,e,a) s0(F)
                    // flops: o0v4L2: 1, o0v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += tempArray[34]("I,e,f") * s0("F");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxaa_LLov += -0.500000 t2_abab(e,a,i,j) m2_abab(I,i,j,b,a) s1_aa(F,b,m)
                    // +               -0.500000 t2_abab(e,a,j,i) m2_abab(I,j,i,b,a) s1_aa(F,b,m)
                    // flops: o1v3L2: 1, o1v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[34]("I,b,e") * s1_aa("F,b,m");
                }
                tempArray[34].~TArrayD();

                if (include_u2_) {

                    // tempArray[35] += 1.000000 t2_abab("b,a,m,i") m2_abab("I,n,i,b,a")
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
                    tempArray[35]("I,m,n") = t2_abab("b,a,m,i") * m2_abab("I,n,i,b,a");
                }

                if (include_u0_ && include_u2_) {

                    // dm_xxaa_LLoo += -0.500000 t2_abab(b,a,m,i) m2_abab(I,n,i,b,a) s0(F)
                    // +               -0.500000 t2_abab(a,b,m,i) m2_abab(I,n,i,a,b) s0(F)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempArray[35]("I,m,n") * s0("F");
                }

                if (include_u1_ && include_u2_) {

                    // RDM_blks_["D1_aa_ov"](F) t2_abab(b,a,m,i) u1_aa(e,j) m2_abab(I,j,i,b,a)
                    // +               -0.500000 r0(F) t2_abab(a,b,m,i) u1_aa(e,j) m2_abab(I,j,i,a,b)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= u1_aa("e,j") * tempArray[35]("I,m,j") * r0("F");

                    // dm_xxaa_LLov += -0.500000 t2_abab(b,a,m,i) m2_abab(I,j,i,b,a) s1_aa(F,e,j)
                    // +               -0.500000 t2_abab(a,b,m,i) m2_abab(I,j,i,a,b) s1_aa(F,e,j)
                    // flops: o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[35]("I,m,j") * s1_aa("F,e,j");
                }

                if (include_u0_ && include_u2_) {

                    // dm_xxaa_LLov += -0.500000 t1_aa(e,j) t2_abab(b,a,m,i) m2_abab(I,j,i,b,a) s0(F)
                    // +               -0.500000 t1_aa(e,j) t2_abab(a,b,m,i) m2_abab(I,j,i,a,b) s0(F)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t1_aa("e,j") * tempArray[35]("I,m,j") * s0("F");
                }
                tempArray[35].~TArrayD();

                if (include_u2_) {

                    // tempArray[36] += 1.000000 m2_abab("I,i,n,b,a") t2_abab("b,a,i,m")
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
                    tempArray[36]("I,n,m") = m2_abab("I,i,n,b,a") * t2_abab("b,a,i,m");
                }

                if (include_u0_ && include_u2_) {

                    // dm_xxbb_LLoo += -0.500000 t2_abab(b,a,i,m) m2_abab(I,i,n,b,a) s0(F)
                    // +               -0.500000 t2_abab(a,b,i,m) m2_abab(I,i,n,a,b) s0(F)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempArray[36]("I,n,m") * s0("F");
                }

                if (include_u1_ && include_u2_) {

                    // RDM_blks_["D1_bb_ov"](F) t2_abab(b,a,i,m) u1_bb(e,j) m2_abab(I,i,j,b,a)
                    // +               -0.500000 r0(F) t2_abab(a,b,i,m) u1_bb(e,j) m2_abab(I,i,j,a,b)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= u1_bb("e,j") * tempArray[36]("I,j,m") * r0("F");

                    // dm_xxbb_LLov += -0.500000 t2_abab(b,a,i,m) m2_abab(I,i,j,b,a) s1_bb(F,e,j)
                    // +               -0.500000 t2_abab(a,b,i,m) m2_abab(I,i,j,a,b) s1_bb(F,e,j)
                    // flops: o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[36]("I,j,m") * s1_bb("F,e,j");
                }

                if (include_u0_ && include_u2_) {

                    // dm_xxbb_LLov += -0.500000 t1_bb(e,j) t2_abab(b,a,i,m) m2_abab(I,i,j,b,a) s0(F)
                    // +               -0.500000 t1_bb(e,j) t2_abab(a,b,i,m) m2_abab(I,i,j,a,b) s0(F)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t1_bb("e,j") * tempArray[36]("I,j,m") * s0("F");
                }
                tempArray[36].~TArrayD();

                if (include_u2_) {

                    // tempArray[37] += 1.000000 m2_abab("I,n,i,a,b") t1_aa("a,m")
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
                    tempArray[37]("I,b,n,m,i") = m2_abab("I,n,i,a,b") * t1_aa("a,m");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxaa_LLoo += -1.000000 t1_aa(a,m) m2_abab(I,n,i,a,b) s1_bb(F,b,i)
                    // flops: o2v2L2: 1, o3v1L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempArray[37]("I,b,n,m,i") * s1_bb("F,b,i");
                }

                if (include_u0_ && include_u2_) {

                    // dm_xxaa_LLov += -0.500000 t1_aa(b,m) t2_abab(e,a,i,j) m2_abab(I,i,j,b,a) s0(F)
                    // +               -0.500000 t1_aa(b,m) t2_abab(e,a,j,i) m2_abab(I,j,i,b,a) s0(F)
                    // flops: o1v3L2: 1, o3v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t2_abab("e,a,i,j") * tempArray[37]("I,a,i,m,j") * s0("F");
                }

                if (include_u2_) {

                    // dm_xxaa_LLov += -0.500000 t1_aa(a,m) m2_abab(I,j,i,a,b) s2_abab(F,e,b,j,i)
                    // +               -0.500000 t1_aa(a,m) m2_abab(I,i,j,a,b) s2_abab(F,e,b,i,j)
                    // flops: o3v2L2: 1, o1v3L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[37]("I,b,j,m,i") * s2_abab("F,e,b,j,i");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxaa_LLov += -1.000000 t1_aa(e,i) t1_aa(a,m) m2_abab(I,i,j,a,b) s1_bb(F,b,j)
                    // flops: o1v3L2: 1, o3v1L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= s1_bb("F,b,j") * tempArray[37]("I,b,i,m,j") * t1_aa("e,i");
                }

                if (include_u2_) {

                    // RDM_blks_["D1_aa_ov"](F) t1_aa(a,m) u2_abab(e,b,i,j) m2_abab(I,i,j,a,b)
                    // +               -0.500000 r0(F) t1_aa(a,m) u2_abab(e,b,j,i) m2_abab(I,j,i,a,b)
                    // flops: o1v3L2: 1, o3v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= u2_abab("e,b,i,j") * tempArray[37]("I,b,i,m,j") * r0("F");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxaa_LLov += -1.000000 r1_bb(F,a,i) t1_aa(b,m) u1_aa(e,j) m2_abab(I,j,i,b,a)
                    // flops: o1v3L2: 1, o3v1L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= r1_bb("F,a,i") * tempArray[37]("I,a,j,m,i") * u1_aa("e,j");
                }
                tempArray[37].~TArrayD();

                if (include_u2_) {

                    // tempArray[38] += 1.000000 t1_bb("a,m") m2_abab("I,i,n,b,a")
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
                    tempArray[38]("I,b,i,m,n") = t1_bb("a,m") * m2_abab("I,i,n,b,a");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxbb_LLoo += -1.000000 t1_bb(a,m) m2_abab(I,i,n,b,a) s1_aa(F,b,i)
                    // flops: o2v2L2: 1, o3v1L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempArray[38]("I,b,i,m,n") * s1_aa("F,b,i");
                }

                if (include_u2_) {

                    // dm_xxbb_LLov += -0.500000 t1_bb(a,m) m2_abab(I,j,i,b,a) s2_abab(F,b,e,j,i)
                    // +               -0.500000 t1_bb(a,m) m2_abab(I,i,j,b,a) s2_abab(F,b,e,i,j)
                    // flops: o3v2L2: 1, o1v3L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[38]("I,b,j,m,i") * s2_abab("F,b,e,j,i");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxbb_LLov += -1.000000 t1_bb(e,i) t1_bb(a,m) m2_abab(I,j,i,b,a) s1_aa(F,b,j)
                    // flops: o1v3L2: 1, o3v1L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= s1_aa("F,b,j") * tempArray[38]("I,b,j,m,i") * t1_bb("e,i");
                }

                if (include_u0_ && include_u2_) {

                    // dm_xxbb_LLov += -0.500000 t1_bb(b,m) t2_abab(a,e,i,j) m2_abab(I,i,j,a,b) s0(F)
                    // +               -0.500000 t1_bb(b,m) t2_abab(a,e,j,i) m2_abab(I,j,i,a,b) s0(F)
                    // flops: o1v3L2: 1, o3v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t2_abab("a,e,i,j") * tempArray[38]("I,a,i,m,j") * s0("F");
                }

                if (include_u2_) {

                    // RDM_blks_["D1_bb_ov"](F) t1_bb(a,m) u2_abab(b,e,i,j) m2_abab(I,i,j,b,a)
                    // +               -0.500000 r0(F) t1_bb(a,m) u2_abab(b,e,j,i) m2_abab(I,j,i,b,a)
                    // flops: o1v3L2: 1, o3v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= u2_abab("b,e,i,j") * tempArray[38]("I,b,i,m,j") * r0("F");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxbb_LLov += -1.000000 r1_aa(F,a,i) t1_bb(b,m) u1_bb(e,j) m2_abab(I,i,j,a,b)
                    // flops: o1v3L2: 1, o3v1L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= r1_aa("F,a,i") * tempArray[38]("I,a,i,m,j") * u1_bb("e,j");
                }
                tempArray[38].~TArrayD();

                if (include_u2_) {

                    // tempArray[39] += 0.500000 m2_aaaa("I,i,n,b,a") t2_aaaa("b,a,i,m")
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
                    tempArray[39]("I,n,m") = 0.500000 * m2_aaaa("I,i,n,b,a") * t2_aaaa("b,a,i,m");
                }

                if (include_u0_ && include_u2_) {

                    // dm_xxaa_LLoo += -0.500000 t2_aaaa(b,a,i,m) m2_aaaa(I,i,n,b,a) s0(F)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempArray[39]("I,n,m") * s0("F");
                }

                if (include_u1_ && include_u2_) {

                    // RDM_blks_["D1_aa_ov"](F) t2_aaaa(b,a,i,m) u1_aa(e,j) m2_aaaa(I,i,j,b,a)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= u1_aa("e,j") * tempArray[39]("I,j,m") * r0("F");
                }

                if (include_u0_ && include_u2_) {

                    // dm_xxaa_LLov += -0.500000 t1_aa(e,j) t2_aaaa(b,a,i,m) m2_aaaa(I,i,j,b,a) s0(F)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t1_aa("e,j") * tempArray[39]("I,j,m") * s0("F");
                }
                tempArray[39].~TArrayD();

                if (include_u2_) {

                    // tempArray[40] += 0.500000 t2_bbbb("b,a,i,m") m2_bbbb("I,i,n,b,a")
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
                    tempArray[40]("I,m,n") = 0.500000 * t2_bbbb("b,a,i,m") * m2_bbbb("I,i,n,b,a");
                }

                if (include_u0_ && include_u2_) {

                    // dm_xxbb_LLoo += -0.500000 t2_bbbb(b,a,i,m) m2_bbbb(I,i,n,b,a) s0(F)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempArray[40]("I,m,n") * s0("F");
                }

                if (include_u1_ && include_u2_) {

                    // RDM_blks_["D1_bb_ov"](F) t2_bbbb(b,a,i,m) u1_bb(e,j) m2_bbbb(I,i,j,b,a)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= u1_bb("e,j") * tempArray[40]("I,m,j") * r0("F");
                }

                if (include_u0_ && include_u2_) {

                    // dm_xxbb_LLov += -0.500000 t1_bb(e,j) t2_bbbb(b,a,i,m) m2_bbbb(I,i,j,b,a) s0(F)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t1_bb("e,j") * tempArray[40]("I,m,j") * s0("F");
                }
                tempArray[40].~TArrayD();

                {

                    // tempArray[41] += 1.000000 l2_abab("I,n,i,b,a") t2_abab("b,a,m,i")
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
                    tempArray[41]("I,n,m") = l2_abab("I,n,i,b,a") * t2_abab("b,a,m,i");

                    // dm_xxaa_LLoo += -0.500000 l2_abab(I,n,i,b,a) r0(F) t2_abab(b,a,m,i)
                    // +               -0.500000 l2_abab(I,n,i,a,b) r0(F) t2_abab(a,b,m,i)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempArray[41]("I,n,m") * r0("F");

                    // dm_xxaa_LLov += -0.500000 l2_abab(I,j,i,b,a) r1_aa(F,e,j) t2_abab(b,a,m,i)
                    // +               -0.500000 l2_abab(I,j,i,a,b) r1_aa(F,e,j) t2_abab(a,b,m,i)
                    // flops: o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[41]("I,j,m") * r1_aa("F,e,j");

                    // dm_xxaa_LLov += -0.500000 l2_abab(I,j,i,b,a) r0(F) t1_aa(e,j) t2_abab(b,a,m,i)
                    // +               -0.500000 l2_abab(I,j,i,a,b) r0(F) t1_aa(e,j) t2_abab(a,b,m,i)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t1_aa("e,j") * tempArray[41]("I,j,m") * r0("F");
                }
                tempArray[41].~TArrayD();

                if (include_u2_) {

                    // tempArray[42] += 1.000000 m2_abab("I,n,i,b,a") u2_abab("b,a,m,i")
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
                    tempArray[42]("I,n,m") = m2_abab("I,n,i,b,a") * u2_abab("b,a,m,i");

                    // RDM_blks_["D1_aa_oo"](F) u2_abab(b,a,m,i) m2_abab(I,n,i,b,a)
                    // +               -0.500000 r0(F) u2_abab(a,b,m,i) m2_abab(I,n,i,a,b)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempArray[42]("I,n,m") * r0("F");

                    // RDM_blks_["D1_aa_ov"](F) t1_aa(e,i) u2_abab(b,a,m,j) m2_abab(I,i,j,b,a)
                    // +               -0.500000 r0(F) t1_aa(e,i) u2_abab(a,b,m,j) m2_abab(I,i,j,a,b)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t1_aa("e,i") * tempArray[42]("I,i,m") * r0("F");

                    // dm_xxaa_LLov += -0.500000 r1_aa(F,e,i) u2_abab(b,a,m,j) m2_abab(I,i,j,b,a)
                    // +               -0.500000 r1_aa(F,e,i) u2_abab(a,b,m,j) m2_abab(I,i,j,a,b)
                    // flops: o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[42]("I,i,m") * r1_aa("F,e,i");
                }
                tempArray[42].~TArrayD();

                {

                    // tempArray[43] += 1.000000 l2_abab("I,i,n,b,a") t2_abab("b,a,i,m")
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
                    tempArray[43]("I,n,m") = l2_abab("I,i,n,b,a") * t2_abab("b,a,i,m");

                    // dm_xxbb_LLoo += -0.500000 l2_abab(I,i,n,b,a) r0(F) t2_abab(b,a,i,m)
                    // +               -0.500000 l2_abab(I,i,n,a,b) r0(F) t2_abab(a,b,i,m)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempArray[43]("I,n,m") * r0("F");

                    // dm_xxbb_LLov += -0.500000 l2_abab(I,i,j,b,a) r1_bb(F,e,j) t2_abab(b,a,i,m)
                    // +               -0.500000 l2_abab(I,i,j,a,b) r1_bb(F,e,j) t2_abab(a,b,i,m)
                    // flops: o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[43]("I,j,m") * r1_bb("F,e,j");

                    // dm_xxbb_LLov += -0.500000 l2_abab(I,i,j,b,a) r0(F) t1_bb(e,j) t2_abab(b,a,i,m)
                    // +               -0.500000 l2_abab(I,i,j,a,b) r0(F) t1_bb(e,j) t2_abab(a,b,i,m)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t1_bb("e,j") * tempArray[43]("I,j,m") * r0("F");
                }
                tempArray[43].~TArrayD();

                if (include_u2_) {

                    // tempArray[44] += 1.000000 u2_abab("a,b,i,m") m2_abab("I,i,n,a,b")
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
                    tempArray[44]("I,m,n") = u2_abab("a,b,i,m") * m2_abab("I,i,n,a,b");

                    // RDM_blks_["D1_bb_oo"](F) u2_abab(b,a,i,m) m2_abab(I,i,n,b,a)
                    // +               -0.500000 r0(F) u2_abab(a,b,i,m) m2_abab(I,i,n,a,b)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempArray[44]("I,m,n") * r0("F");

                    // RDM_blks_["D1_bb_ov"](F) t1_bb(e,i) u2_abab(b,a,j,m) m2_abab(I,j,i,b,a)
                    // +               -0.500000 r0(F) t1_bb(e,i) u2_abab(a,b,j,m) m2_abab(I,j,i,a,b)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t1_bb("e,i") * tempArray[44]("I,m,i") * r0("F");

                    // dm_xxbb_LLov += -0.500000 r1_bb(F,e,i) u2_abab(b,a,j,m) m2_abab(I,j,i,b,a)
                    // +               -0.500000 r1_bb(F,e,i) u2_abab(a,b,j,m) m2_abab(I,j,i,a,b)
                    // flops: o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[44]("I,m,i") * r1_bb("F,e,i");

                    // tempArray[45] += 1.000000 m2_aaaa("I,i,n,a,b") t1_aa("a,m")
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
                    tempArray[45]("I,b,i,n,m") = m2_aaaa("I,i,n,a,b") * t1_aa("a,m");
                }
                tempArray[44].~TArrayD();

                if (include_u1_ && include_u2_) {

                    // dm_xxaa_LLoo += 1.000000 t1_aa(a,m) m2_aaaa(I,i,n,a,b) s1_aa(F,b,i)
                    // flops: o2v2L2: 1, o3v1L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempArray[45]("I,b,i,n,m") * s1_aa("F,b,i");

                    // dm_xxaa_LLov += 1.000000 t1_aa(e,i) t1_aa(a,m) m2_aaaa(I,j,i,a,b) s1_aa(F,b,j)
                    // flops: o1v3L2: 1, o3v1L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += s1_aa("F,b,j") * tempArray[45]("I,b,j,i,m") * t1_aa("e,i");
                }

                if (include_u0_ && include_u2_) {

                    // dm_xxaa_LLov += -0.500000 t1_aa(b,m) t2_aaaa(e,a,i,j) m2_aaaa(I,i,j,b,a) s0(F)
                    // flops: o1v3L2: 1, o3v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= 0.500000 * t2_aaaa("e,a,i,j") * tempArray[45]("I,a,i,j,m") * s0("F");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxaa_LLov += 1.000000 r1_aa(F,a,i) t1_aa(b,m) u1_aa(e,j) m2_aaaa(I,i,j,b,a)
                    // flops: o1v3L2: 1, o3v1L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += r1_aa("F,a,i") * tempArray[45]("I,a,i,j,m") * u1_aa("e,j");
                }

                if (include_u2_) {

                    // dm_xxaa_LLov += -0.500000 t1_aa(a,m) m2_aaaa(I,j,i,a,b) s2_aaaa(F,e,b,j,i)
                    // flops: o3v2L2: 1, o1v3L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= 0.500000 * tempArray[45]("I,b,j,i,m") * s2_aaaa("F,e,b,j,i");

                    // tempArray[46] += 1.000000 m2_bbbb("I,i,n,a,b") t1_bb("a,m")
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
                    tempArray[46]("I,b,i,n,m") = m2_bbbb("I,i,n,a,b") * t1_bb("a,m");
                }
                tempArray[45].~TArrayD();

                if (include_u1_ && include_u2_) {

                    // dm_xxbb_LLoo += 1.000000 t1_bb(a,m) m2_bbbb(I,i,n,a,b) s1_bb(F,b,i)
                    // flops: o2v2L2: 1, o3v1L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempArray[46]("I,b,i,n,m") * s1_bb("F,b,i");

                    // dm_xxbb_LLov += 1.000000 t1_bb(e,i) t1_bb(a,m) m2_bbbb(I,j,i,a,b) s1_bb(F,b,j)
                    // flops: o1v3L2: 1, o3v1L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += s1_bb("F,b,j") * tempArray[46]("I,b,j,i,m") * t1_bb("e,i");
                }

                if (include_u0_ && include_u2_) {

                    // dm_xxbb_LLov += -0.500000 t1_bb(b,m) t2_bbbb(e,a,i,j) m2_bbbb(I,i,j,b,a) s0(F)
                    // flops: o1v3L2: 1, o3v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= 0.500000 * t2_bbbb("e,a,i,j") * tempArray[46]("I,a,i,j,m") * s0("F");
                }

                if (include_u2_) {

                    // dm_xxbb_LLov += -0.500000 t1_bb(a,m) m2_bbbb(I,j,i,a,b) s2_bbbb(F,e,b,j,i)
                    // flops: o3v2L2: 1, o1v3L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= 0.500000 * tempArray[46]("I,b,j,i,m") * s2_bbbb("F,e,b,j,i");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxbb_LLov += 1.000000 r1_bb(F,a,i) t1_bb(b,m) u1_bb(e,j) m2_bbbb(I,i,j,b,a)
                    // flops: o1v3L2: 1, o3v1L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += r1_bb("F,a,i") * tempArray[46]("I,a,i,j,m") * u1_bb("e,j");
                }
                tempArray[46].~TArrayD();

                {

                    // tempArray[47] += 0.500000 t2_bbbb("b,a,i,m") l2_bbbb("I,i,n,b,a")
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
                    tempArray[47]("I,m,n") = 0.500000 * t2_bbbb("b,a,i,m") * l2_bbbb("I,i,n,b,a");

                    // dm_xxbb_LLoo += -0.500000 l2_bbbb(I,i,n,b,a) r0(F) t2_bbbb(b,a,i,m)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempArray[47]("I,m,n") * r0("F");

                    // dm_xxbb_LLov += -0.500000 l2_bbbb(I,i,j,b,a) r0(F) t1_bb(e,j) t2_bbbb(b,a,i,m)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t1_bb("e,j") * tempArray[47]("I,m,j") * r0("F");

                    // tempArray[48] += 0.500000 t2_aaaa("b,a,i,m") l2_aaaa("I,i,n,b,a")
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
                    tempArray[48]("I,m,n") = 0.500000 * t2_aaaa("b,a,i,m") * l2_aaaa("I,i,n,b,a");

                    // dm_xxaa_LLoo += -0.500000 l2_aaaa(I,i,n,b,a) r0(F) t2_aaaa(b,a,i,m)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempArray[48]("I,m,n") * r0("F");

                    // dm_xxaa_LLov += -0.500000 l2_aaaa(I,i,j,b,a) r0(F) t1_aa(e,j) t2_aaaa(b,a,i,m)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t1_aa("e,j") * tempArray[48]("I,m,j") * r0("F");
                }
                tempArray[47].~TArrayD();
                tempArray[48].~TArrayD();

                if (include_u2_) {

                    // tempArray[49] += 0.500000 u2_aaaa("b,a,j,m") m2_aaaa("I,i,j,b,a")
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
                    tempArray[49]("I,m,i") = 0.500000 * u2_aaaa("b,a,j,m") * m2_aaaa("I,i,j,b,a");

                    // dm_xxaa_LLov += 0.500000 r1_aa(F,e,i) u2_aaaa(b,a,j,m) m2_aaaa(I,i,j,b,a)
                    // flops: o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[49]("I,m,i") * r1_aa("F,e,i");

                    // RDM_blks_["D1_aa_ov"](F) t1_aa(e,i) u2_aaaa(b,a,j,m) m2_aaaa(I,i,j,b,a)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += t1_aa("e,i") * tempArray[49]("I,m,i") * r0("F");

                    // tempArray[50] += 0.500000 u2_bbbb("b,a,j,m") m2_bbbb("I,i,j,b,a")
                    // flops: o3v2L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
                    tempArray[50]("I,m,i") = 0.500000 * u2_bbbb("b,a,j,m") * m2_bbbb("I,i,j,b,a");

                    // dm_xxbb_LLov += 0.500000 r1_bb(F,e,i) u2_bbbb(b,a,j,m) m2_bbbb(I,i,j,b,a)
                    // flops: o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[50]("I,m,i") * r1_bb("F,e,i");

                    // RDM_blks_["D1_bb_ov"](F) t1_bb(e,i) u2_bbbb(b,a,j,m) m2_bbbb(I,i,j,b,a)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += t1_bb("e,i") * tempArray[50]("I,m,i") * r0("F");
                }
                tempArray[49].~TArrayD();
                tempArray[50].~TArrayD();

                {

                    // tempArray[51] += 1.000000 t1_bb("a,m") l2_abab("I,i,n,b,a")
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
                    tempArray[51]("I,b,i,m,n") = t1_bb("a,m") * l2_abab("I,i,n,b,a");

                    // dm_xxbb_LLoo += -1.000000 l2_abab(I,i,n,b,a) r1_aa(F,b,i) t1_bb(a,m)
                    // flops: o2v2L2: 1, o3v1L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempArray[51]("I,b,i,m,n") * r1_aa("F,b,i");

                    // dm_xxbb_LLov += -0.500000 l2_abab(I,i,j,a,b) r0(F) t1_bb(b,m) t2_abab(a,e,i,j)
                    // +               -0.500000 l2_abab(I,j,i,a,b) r0(F) t1_bb(b,m) t2_abab(a,e,j,i)
                    // flops: o1v3L2: 1, o3v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t2_abab("a,e,i,j") * tempArray[51]("I,a,i,m,j") * r0("F");

                    // dm_xxbb_LLov += -1.000000 l2_abab(I,j,i,b,a) r1_aa(F,b,j) t1_bb(e,i) t1_bb(a,m)
                    // flops: o1v3L2: 1, o3v1L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= r1_aa("F,b,j") * tempArray[51]("I,b,j,m,i") * t1_bb("e,i");

                    // dm_xxbb_LLov += -0.500000 l2_abab(I,j,i,b,a) r2_abab(F,b,e,j,i) t1_bb(a,m)
                    // +               -0.500000 l2_abab(I,i,j,b,a) r2_abab(F,b,e,i,j) t1_bb(a,m)
                    // flops: o3v2L2: 1, o1v3L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[51]("I,b,j,m,i") * r2_abab("F,b,e,j,i");

                    // tempArray[52] += 1.000000 l2_aaaa("I,i,n,a,b") t1_aa("a,m")
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
                    tempArray[52]("I,b,i,n,m") = l2_aaaa("I,i,n,a,b") * t1_aa("a,m");

                    // dm_xxaa_LLoo += 1.000000 l2_aaaa(I,i,n,a,b) r1_aa(F,b,i) t1_aa(a,m)
                    // flops: o2v2L2: 1, o3v1L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempArray[52]("I,b,i,n,m") * r1_aa("F,b,i");

                    // dm_xxaa_LLov += -0.500000 l2_aaaa(I,j,i,a,b) r2_aaaa(F,e,b,j,i) t1_aa(a,m)
                    // flops: o3v2L2: 1, o1v3L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= 0.500000 * tempArray[52]("I,b,j,i,m") * r2_aaaa("F,e,b,j,i");

                    // dm_xxaa_LLov += -0.500000 l2_aaaa(I,i,j,b,a) r0(F) t1_aa(b,m) t2_aaaa(e,a,i,j)
                    // flops: o1v3L2: 1, o3v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= 0.500000 * t2_aaaa("e,a,i,j") * tempArray[52]("I,a,i,j,m") * r0("F");

                    // dm_xxaa_LLov += 1.000000 l2_aaaa(I,j,i,a,b) r1_aa(F,b,j) t1_aa(e,i) t1_aa(a,m)
                    // flops: o1v3L2: 1, o3v1L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += r1_aa("F,b,j") * tempArray[52]("I,b,j,i,m") * t1_aa("e,i");
                }
                tempArray[51].~TArrayD();
                tempArray[52].~TArrayD();

                if (include_u1_ && include_u2_) {

                    // tempArray[53] += 1.000000 u1_aa("b,m") m2_aaaa("I,i,n,b,a")
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
                    tempArray[53]("I,a,m,i,n") = u1_aa("b,m") * m2_aaaa("I,i,n,b,a");

                    // dm_xxaa_LLoo += 1.000000 r1_aa(F,a,i) u1_aa(b,m) m2_aaaa(I,i,n,b,a)
                    // flops: o2v2L2: 1, o3v1L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempArray[53]("I,a,m,i,n") * r1_aa("F,a,i");

                    // RDM_blks_["D1_aa_ov"](F) t2_aaaa(e,a,i,j) u1_aa(b,m) m2_aaaa(I,i,j,b,a)
                    // flops: o1v3L2: 1, o3v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= 0.500000 * t2_aaaa("e,a,i,j") * tempArray[53]("I,a,m,i,j") * r0("F");

                    // dm_xxaa_LLov += -0.500000 r2_aaaa(F,e,a,i,j) u1_aa(b,m) m2_aaaa(I,i,j,b,a)
                    // flops: o3v2L2: 1, o1v3L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= 0.500000 * tempArray[53]("I,a,m,i,j") * r2_aaaa("F,e,a,i,j");

                    // dm_xxaa_LLov += 1.000000 r1_aa(F,a,i) t1_aa(e,j) u1_aa(b,m) m2_aaaa(I,i,j,b,a)
                    // flops: o1v3L2: 1, o3v1L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += r1_aa("F,a,i") * tempArray[53]("I,a,m,i,j") * t1_aa("e,j");

                    // tempArray[54] += 1.000000 m2_abab("I,n,i,b,a") u1_aa("b,m")
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
                    tempArray[54]("I,a,n,m,i") = m2_abab("I,n,i,b,a") * u1_aa("b,m");

                    // dm_xxaa_LLoo += -1.000000 r1_bb(F,a,i) u1_aa(b,m) m2_abab(I,n,i,b,a)
                    // flops: o2v2L2: 1, o3v1L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempArray[54]("I,a,n,m,i") * r1_bb("F,a,i");

                    // dm_xxaa_LLov += -1.000000 r1_bb(F,a,i) t1_aa(e,j) u1_aa(b,m) m2_abab(I,j,i,b,a)
                    // flops: o1v3L2: 1, o3v1L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= r1_bb("F,a,i") * tempArray[54]("I,a,j,m,i") * t1_aa("e,j");

                    // RDM_blks_["D1_aa_ov"](F) t2_abab(e,a,i,j) u1_aa(b,m) m2_abab(I,i,j,b,a)
                    // +               -0.500000 r0(F) t2_abab(e,a,j,i) u1_aa(b,m) m2_abab(I,j,i,b,a)
                    // flops: o1v3L2: 1, o3v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t2_abab("e,a,i,j") * tempArray[54]("I,a,i,m,j") * r0("F");

                    // dm_xxaa_LLov += -0.500000 r2_abab(F,e,a,i,j) u1_aa(b,m) m2_abab(I,i,j,b,a)
                    // +               -0.500000 r2_abab(F,e,a,j,i) u1_aa(b,m) m2_abab(I,j,i,b,a)
                    // flops: o3v2L2: 1, o1v3L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[54]("I,a,i,m,j") * r2_abab("F,e,a,i,j");

                    // tempArray[55] += 1.000000 u1_bb("b,m") m2_bbbb("I,i,n,b,a")
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
                    tempArray[55]("I,a,m,i,n") = u1_bb("b,m") * m2_bbbb("I,i,n,b,a");

                    // dm_xxbb_LLoo += 1.000000 r1_bb(F,a,i) u1_bb(b,m) m2_bbbb(I,i,n,b,a)
                    // flops: o2v2L2: 1, o3v1L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempArray[55]("I,a,m,i,n") * r1_bb("F,a,i");

                    // RDM_blks_["D1_bb_ov"](F) t2_bbbb(e,a,i,j) u1_bb(b,m) m2_bbbb(I,i,j,b,a)
                    // flops: o1v3L2: 1, o3v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= 0.500000 * t2_bbbb("e,a,i,j") * tempArray[55]("I,a,m,i,j") * r0("F");

                    // dm_xxbb_LLov += 1.000000 r1_bb(F,a,i) t1_bb(e,j) u1_bb(b,m) m2_bbbb(I,i,j,b,a)
                    // flops: o1v3L2: 1, o3v1L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += r1_bb("F,a,i") * tempArray[55]("I,a,m,i,j") * t1_bb("e,j");

                    // dm_xxbb_LLov += -0.500000 r2_bbbb(F,e,a,i,j) u1_bb(b,m) m2_bbbb(I,i,j,b,a)
                    // flops: o3v2L2: 1, o1v3L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= 0.500000 * tempArray[55]("I,a,m,i,j") * r2_bbbb("F,e,a,i,j");
                }
                tempArray[53].~TArrayD();
                tempArray[54].~TArrayD();
                tempArray[55].~TArrayD();

                {

                    // tempArray[56] += 1.000000 l2_bbbb("I,i,n,a,b") t1_bb("a,m")
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
                    tempArray[56]("I,b,i,n,m") = l2_bbbb("I,i,n,a,b") * t1_bb("a,m");

                    // dm_xxbb_LLoo += 1.000000 l2_bbbb(I,i,n,a,b) r1_bb(F,b,i) t1_bb(a,m)
                    // flops: o2v2L2: 1, o3v1L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempArray[56]("I,b,i,n,m") * r1_bb("F,b,i");

                    // dm_xxbb_LLov += -0.500000 l2_bbbb(I,j,i,a,b) r2_bbbb(F,e,b,j,i) t1_bb(a,m)
                    // flops: o3v2L2: 1, o1v3L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= 0.500000 * tempArray[56]("I,b,j,i,m") * r2_bbbb("F,e,b,j,i");

                    // dm_xxbb_LLov += -0.500000 l2_bbbb(I,i,j,b,a) r0(F) t1_bb(b,m) t2_bbbb(e,a,i,j)
                    // flops: o1v3L2: 1, o3v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= 0.500000 * t2_bbbb("e,a,i,j") * tempArray[56]("I,a,i,j,m") * r0("F");

                    // dm_xxbb_LLov += 1.000000 l2_bbbb(I,j,i,a,b) r1_bb(F,b,j) t1_bb(e,i) t1_bb(a,m)
                    // flops: o1v3L2: 1, o3v1L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += r1_bb("F,b,j") * tempArray[56]("I,b,j,i,m") * t1_bb("e,i");
                }
                tempArray[56].~TArrayD();

                if (include_u1_ && include_u2_) {

                    // tempArray[57] += 1.000000 u1_bb("b,m") m2_abab("I,i,n,a,b")
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
                    tempArray[57]("I,a,i,m,n") = u1_bb("b,m") * m2_abab("I,i,n,a,b");

                    // dm_xxbb_LLoo += -1.000000 r1_aa(F,a,i) u1_bb(b,m) m2_abab(I,i,n,a,b)
                    // flops: o2v2L2: 1, o3v1L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempArray[57]("I,a,i,m,n") * r1_aa("F,a,i");

                    // RDM_blks_["D1_bb_ov"](F) t2_abab(a,e,i,j) u1_bb(b,m) m2_abab(I,i,j,a,b)
                    // +               -0.500000 r0(F) t2_abab(a,e,j,i) u1_bb(b,m) m2_abab(I,j,i,a,b)
                    // flops: o1v3L2: 1, o3v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t2_abab("a,e,i,j") * tempArray[57]("I,a,i,m,j") * r0("F");

                    // dm_xxbb_LLov += -1.000000 r1_aa(F,a,i) t1_bb(e,j) u1_bb(b,m) m2_abab(I,i,j,a,b)
                    // flops: o1v3L2: 1, o3v1L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= r1_aa("F,a,i") * tempArray[57]("I,a,i,m,j") * t1_bb("e,j");

                    // dm_xxbb_LLov += -0.500000 r2_abab(F,a,e,i,j) u1_bb(b,m) m2_abab(I,i,j,a,b)
                    // +               -0.500000 r2_abab(F,a,e,j,i) u1_bb(b,m) m2_abab(I,j,i,a,b)
                    // flops: o3v2L2: 1, o1v3L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[57]("I,a,i,m,j") * r2_abab("F,a,e,i,j");
                }
                tempArray[57].~TArrayD();

                {

                    // tempArray[58] += 1.000000 l2_abab("I,n,i,a,b") t1_aa("a,m")
                    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
                    tempArray[58]("I,b,n,m,i") = l2_abab("I,n,i,a,b") * t1_aa("a,m");

                    // dm_xxaa_LLoo += -1.000000 l2_abab(I,n,i,a,b) r1_bb(F,b,i) t1_aa(a,m)
                    // flops: o2v2L2: 1, o3v1L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempArray[58]("I,b,n,m,i") * r1_bb("F,b,i");

                    // dm_xxaa_LLov += -0.500000 l2_abab(I,j,i,a,b) r2_abab(F,e,b,j,i) t1_aa(a,m)
                    // +               -0.500000 l2_abab(I,i,j,a,b) r2_abab(F,e,b,i,j) t1_aa(a,m)
                    // flops: o3v2L2: 1, o1v3L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[58]("I,b,j,m,i") * r2_abab("F,e,b,j,i");

                    // dm_xxaa_LLov += -1.000000 l2_abab(I,i,j,a,b) r1_bb(F,b,j) t1_aa(e,i) t1_aa(a,m)
                    // flops: o1v3L2: 1, o3v1L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= r1_bb("F,b,j") * tempArray[58]("I,b,i,m,j") * t1_aa("e,i");

                    // dm_xxaa_LLov += -0.500000 l2_abab(I,i,j,b,a) r0(F) t1_aa(b,m) t2_abab(e,a,i,j)
                    // +               -0.500000 l2_abab(I,j,i,b,a) r0(F) t1_aa(b,m) t2_abab(e,a,j,i)
                    // flops: o1v3L2: 1, o3v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t2_abab("e,a,i,j") * tempArray[58]("I,a,i,m,j") * r0("F");
                }
                tempArray[58].~TArrayD();

                if (include_u1_) {

                    // tempArray[59] += 1.000000 m1_aa("I,i,a") s1_aa("F,a,i")
                    // flops: o1v1L2: 1, o0v0L2: 1 | mem: o0v0L2: 2,
                    tempArray[59]("I,F") = m1_aa("I,i,a") * s1_aa("F,a,i");

                    // dm_xxaa_LLoo += 1.000000 d_aa(m,n) m1_aa(I,i,a) s1_aa(F,a,i)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempArray[59]("I,F") * Id_blks["aa_oo"]("m,n");

                    // dm_xxbb_LLoo += 1.000000 d_bb(m,n) m1_aa(I,i,a) s1_aa(F,a,i)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempArray[59]("I,F") * Id_blks["bb_oo"]("m,n");

                    // dm_xxaa_LLov += 1.000000 t1_aa(e,m) m1_aa(I,i,a) s1_aa(F,a,i)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[59]("I,F") * t1_aa("e,m");

                    // dm_xxbb_LLov += 1.000000 t1_bb(e,m) m1_aa(I,i,a) s1_aa(F,a,i)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[59]("I,F") * t1_bb("e,m");
                }
                tempArray[59].~TArrayD();

                {

                    // tempArray[60] += 1.000000 l1_bb("I,i,a") r1_bb("F,a,i")
                    // flops: o1v1L2: 1, o0v0L2: 1 | mem: o0v0L2: 2,
                    tempArray[60]("I,F") = l1_bb("I,i,a") * r1_bb("F,a,i");

                    // dm_xxaa_LLoo += 1.000000 d_aa(m,n) l1_bb(I,i,a) r1_bb(F,a,i)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempArray[60]("I,F") * Id_blks["aa_oo"]("m,n");

                    // dm_xxbb_LLoo += 1.000000 d_bb(m,n) l1_bb(I,i,a) r1_bb(F,a,i)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempArray[60]("I,F") * Id_blks["bb_oo"]("m,n");

                    // dm_xxaa_LLov += 1.000000 l1_bb(I,i,a) r1_bb(F,a,i) t1_aa(e,m)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[60]("I,F") * t1_aa("e,m");

                    // dm_xxbb_LLov += 1.000000 l1_bb(I,i,a) r1_bb(F,a,i) t1_bb(e,m)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[60]("I,F") * t1_bb("e,m");
                }
                tempArray[60].~TArrayD();

                if (include_u1_) {

                    // tempArray[61] += 1.000000 m1_bb("I,i,a") s1_bb("F,a,i")
                    // flops: o1v1L2: 1, o0v0L2: 1 | mem: o0v0L2: 2,
                    tempArray[61]("I,F") = m1_bb("I,i,a") * s1_bb("F,a,i");

                    // dm_xxaa_LLoo += 1.000000 d_aa(m,n) m1_bb(I,i,a) s1_bb(F,a,i)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempArray[61]("I,F") * Id_blks["aa_oo"]("m,n");

                    // dm_xxbb_LLoo += 1.000000 d_bb(m,n) m1_bb(I,i,a) s1_bb(F,a,i)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempArray[61]("I,F") * Id_blks["bb_oo"]("m,n");

                    // dm_xxaa_LLov += 1.000000 t1_aa(e,m) m1_bb(I,i,a) s1_bb(F,a,i)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[61]("I,F") * t1_aa("e,m");

                    // dm_xxbb_LLov += 1.000000 t1_bb(e,m) m1_bb(I,i,a) s1_bb(F,a,i)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[61]("I,F") * t1_bb("e,m");
                }
                tempArray[61].~TArrayD();

                {

                    // tempArray[62] += 1.000000 l1_aa("I,i,a") r1_aa("F,a,i")
                    // flops: o1v1L2: 1, o0v0L2: 1 | mem: o0v0L2: 2,
                    tempArray[62]("I,F") = l1_aa("I,i,a") * r1_aa("F,a,i");

                    // dm_xxaa_LLoo += 1.000000 d_aa(m,n) l1_aa(I,i,a) r1_aa(F,a,i)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += tempArray[62]("I,F") * Id_blks["aa_oo"]("m,n");

                    // dm_xxbb_LLoo += 1.000000 d_bb(m,n) l1_aa(I,i,a) r1_aa(F,a,i)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += tempArray[62]("I,F") * Id_blks["bb_oo"]("m,n");

                    // dm_xxaa_LLov += 1.000000 l1_aa(I,i,a) r1_aa(F,a,i) t1_aa(e,m)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[62]("I,F") * t1_aa("e,m");

                    // dm_xxbb_LLov += 1.000000 l1_aa(I,i,a) r1_aa(F,a,i) t1_bb(e,m)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[62]("I,F") * t1_bb("e,m");
                }
                tempArray[62].~TArrayD();

                if (include_u1_) {

                    // tempArray[63] += 1.000000 m1_aa("I,i,a") r1_aa("F,a,i")
                    // flops: o1v1L2: 1, o0v0L2: 1 | mem: o0v0L2: 2,
                    tempArray[63]("I,F") = m1_aa("I,i,a") * r1_aa("F,a,i");

                    // dm_xxaa_LLov += 1.000000 r1_aa(F,a,i) u1_aa(e,m) m1_aa(I,i,a)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[63]("I,F") * u1_aa("e,m");

                    // dm_xxbb_LLov += 1.000000 r1_aa(F,a,i) u1_bb(e,m) m1_aa(I,i,a)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[63]("I,F") * u1_bb("e,m");

                    // tempArray[64] += 1.000000 r1_bb("F,a,i") m1_bb("I,i,a")
                    // flops: o1v1L2: 1, o0v0L2: 1 | mem: o0v0L2: 2,
                    tempArray[64]("I,F") = r1_bb("F,a,i") * m1_bb("I,i,a");

                    // dm_xxaa_LLov += 1.000000 r1_bb(F,a,i) u1_aa(e,m) m1_bb(I,i,a)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += tempArray[64]("I,F") * u1_aa("e,m");

                    // dm_xxbb_LLov += 1.000000 r1_bb(F,a,i) u1_bb(e,m) m1_bb(I,i,a)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += tempArray[64]("I,F") * u1_bb("e,m");

                    // tempArray[65] += 1.000000 m1_aa("I,n,a") t1_aa("a,m")
                    // flops: o2v1L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
                    tempArray[65]("I,n,m") = m1_aa("I,n,a") * t1_aa("a,m");
                }
                tempArray[63].~TArrayD();
                tempArray[64].~TArrayD();

                if (include_u0_ && include_u1_) {

                    // dm_xxaa_LLoo += -1.000000 t1_aa(a,m) m1_aa(I,n,a) s0(F)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= tempArray[65]("I,n,m") * s0("F");
                }

                if (include_u1_) {

                    // dm_xxaa_LLov += -1.000000 t1_aa(a,m) m1_aa(I,i,a) s1_aa(F,e,i)
                    // flops: o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= tempArray[65]("I,i,m") * s1_aa("F,e,i");

                    // RDM_blks_["D1_aa_ov"](F) t1_aa(a,m) u1_aa(e,i) m1_aa(I,i,a)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= u1_aa("e,i") * tempArray[65]("I,i,m") * r0("F");
                }

                if (include_u0_ && include_u1_) {

                    // dm_xxaa_LLov += -1.000000 t1_aa(e,i) t1_aa(a,m) m1_aa(I,i,a) s0(F)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t1_aa("e,i") * tempArray[65]("I,i,m") * s0("F");
                }
                tempArray[65].~TArrayD();

                if (include_u1_) {

                    // tempArray[66] += 1.000000 m1_bb("I,n,a") t1_bb("a,m")
                    // flops: o2v1L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
                    tempArray[66]("I,n,m") = m1_bb("I,n,a") * t1_bb("a,m");
                }

                if (include_u0_ && include_u1_) {

                    // dm_xxbb_LLoo += -1.000000 t1_bb(a,m) m1_bb(I,n,a) s0(F)
                    // flops: o2v2L2: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= tempArray[66]("I,n,m") * s0("F");
                }

                if (include_u1_) {

                    // dm_xxbb_LLov += -1.000000 t1_bb(a,m) m1_bb(I,i,a) s1_bb(F,e,i)
                    // flops: o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= tempArray[66]("I,i,m") * s1_bb("F,e,i");
                }

                if (include_u0_ && include_u1_) {

                    // dm_xxbb_LLov += -1.000000 t1_bb(e,i) t1_bb(a,m) m1_bb(I,i,a) s0(F)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t1_bb("e,i") * tempArray[66]("I,i,m") * s0("F");
                }

                if (include_u1_) {

                    // RDM_blks_["D1_bb_ov"](F) t1_bb(a,m) u1_bb(e,i) m1_bb(I,i,a)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= u1_bb("e,i") * tempArray[66]("I,i,m") * r0("F");
                }
                tempArray[66].~TArrayD();

                if (include_u0_) {

                    // dm_xxaa_LLoo += 1.000000 d_aa(m,n) m0(I) s0(F)
                    // flops: o2v2L2: 1, o2v0L2: 1, o0v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1, o0v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += m0("I") * s0("F") * Id_blks["aa_oo"]("m,n");
                }

                if (include_u2_) {

                    // RDM_blks_["D1_aa_oo"](F) u2_aaaa(b,a,i,m) m2_aaaa(I,i,n,b,a)
                    // flops: o2v2L2: 1, o3v2L1: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1, o2v0L1: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= 0.500000 * u2_aaaa("b,a,i,m") * m2_aaaa("I,i,n,b,a") * r0("F");
                }

                {

                    // dm_xxaa_LLoo += 1.000000 d_aa(m,n) l0(I) r0(F)
                    // flops: o2v2L2: 1, o2v0L2: 1, o0v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1, o0v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") += l0("I") * r0("F") * Id_blks["aa_oo"]("m,n");
                }

                if (include_u1_) {

                    // RDM_blks_["D1_aa_oo"](F) u1_aa(a,m) m1_aa(I,n,a)
                    // flops: o2v2L2: 1, o2v0L2: 1, o2v1L1: 1 | mem: o2v2L2: 1, o2v0L2: 1, o2v0L1: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= u1_aa("a,m") * m1_aa("I,n,a") * r0("F");
                }

                {

                    // dm_xxaa_LLoo += -1.000000 l1_aa(I,n,a) r0(F) t1_aa(a,m)
                    // flops: o2v2L2: 1, o2v0L2: 1, o2v1L1: 1 | mem: o2v2L2: 1, o2v0L2: 1, o2v0L1: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= l1_aa("I,n,a") * t1_aa("a,m") * r0("F");
                }

                if (include_u1_) {

                    // dm_xxaa_LLoo += -1.000000 m1_aa(I,n,a) s1_aa(F,a,m)
                    // flops: o2v2L2: 1, o2v1L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= m1_aa("I,n,a") * s1_aa("F,a,m");
                }

                {

                    // dm_xxaa_LLoo += -1.000000 l1_aa(I,n,a) r1_aa(F,a,m)
                    // flops: o2v2L2: 1, o2v1L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_oo"]("I,F,m,n") -= l1_aa("I,n,a") * r1_aa("F,a,m");

                    // dm_xxbb_LLoo += -1.000000 l1_bb(I,n,a) r0(F) t1_bb(a,m)
                    // flops: o2v2L2: 1, o2v0L2: 1, o2v1L1: 1 | mem: o2v2L2: 1, o2v0L2: 1, o2v0L1: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= l1_bb("I,n,a") * t1_bb("a,m") * r0("F");
                }

                if (include_u1_) {

                    // dm_xxbb_LLoo += -1.000000 m1_bb(I,n,a) s1_bb(F,a,m)
                    // flops: o2v2L2: 1, o2v1L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= m1_bb("I,n,a") * s1_bb("F,a,m");
                }

                if (include_u0_) {

                    // dm_xxbb_LLoo += 1.000000 d_bb(m,n) m0(I) s0(F)
                    // flops: o2v2L2: 1, o2v0L2: 1, o0v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1, o0v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += m0("I") * s0("F") * Id_blks["bb_oo"]("m,n");
                }

                {

                    // dm_xxbb_LLoo += 1.000000 d_bb(m,n) l0(I) r0(F)
                    // flops: o2v2L2: 1, o2v0L2: 1, o0v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1, o0v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") += l0("I") * r0("F") * Id_blks["bb_oo"]("m,n");
                }

                if (include_u1_) {

                    // RDM_blks_["D1_bb_oo"](F) u1_bb(a,m) m1_bb(I,n,a)
                    // flops: o2v2L2: 1, o2v0L2: 1, o2v1L1: 1 | mem: o2v2L2: 1, o2v0L2: 1, o2v0L1: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= u1_bb("a,m") * m1_bb("I,n,a") * r0("F");
                }

                if (include_u2_) {

                    // RDM_blks_["D1_bb_oo"](F) u2_bbbb(b,a,i,m) m2_bbbb(I,i,n,b,a)
                    // flops: o2v2L2: 1, o3v2L1: 1, o2v0L2: 1 | mem: o2v2L2: 1, o2v0L2: 1, o2v0L1: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= 0.500000 * u2_bbbb("b,a,i,m") * m2_bbbb("I,i,n,b,a") * r0("F");
                }

                {

                    // dm_xxbb_LLoo += -1.000000 l1_bb(I,n,a) r1_bb(F,a,m)
                    // flops: o2v2L2: 1, o2v1L2: 1 | mem: o2v2L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_oo"]("I,F,m,n") -= l1_bb("I,n,a") * r1_bb("F,a,m");

                    // dm_xxaa_LLvv += 0.500000 l2_abab(I,j,i,e,a) r2_abab(F,f,a,j,i)
                    // +               0.500000 l2_abab(I,i,j,e,a) r2_abab(F,f,a,i,j)
                    // flops: o2v3L2: 1, o0v4L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += l2_abab("I,j,i,e,a") * r2_abab("F,f,a,j,i");
                }

                if (include_u2_) {

                    // RDM_blks_["D1_aa_vv"](F) u2_aaaa(f,a,i,j) m2_aaaa(I,i,j,e,a)
                    // flops: o0v4L2: 1, o2v3L1: 1, o0v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1, o0v2L1: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += 0.500000 * u2_aaaa("f,a,i,j") * m2_aaaa("I,i,j,e,a") * r0("F");

                    // dm_xxaa_LLvv += 0.500000 m2_abab(I,j,i,e,a) s2_abab(F,f,a,j,i)
                    // +               0.500000 m2_abab(I,i,j,e,a) s2_abab(F,f,a,i,j)
                    // flops: o2v3L2: 1, o0v4L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += m2_abab("I,j,i,e,a") * s2_abab("F,f,a,j,i");
                }

                {

                    // dm_xxaa_LLvv += 0.500000 l2_aaaa(I,j,i,e,a) r2_aaaa(F,f,a,j,i)
                    // flops: o2v3L2: 1, o0v4L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += 0.500000 * l2_aaaa("I,j,i,e,a") * r2_aaaa("F,f,a,j,i");

                    // dm_xxaa_LLvv += 0.500000 l2_aaaa(I,i,j,e,a) r0(F) t2_aaaa(f,a,i,j)
                    // flops: o0v4L2: 1, o2v3L1: 1, o0v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1, o0v2L1: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += 0.500000 * l2_aaaa("I,i,j,e,a") * t2_aaaa("f,a,i,j") * r0("F");
                }

                if (include_u0_ && include_u1_) {

                    // dm_xxaa_LLvv += 1.000000 t1_aa(f,i) m1_aa(I,i,e) s0(F)
                    // flops: o0v4L2: 1, o0v2L2: 1, o1v2L1: 1 | mem: o0v4L2: 1, o0v2L2: 1, o0v2L1: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += t1_aa("f,i") * m1_aa("I,i,e") * s0("F");
                }

                {

                    // dm_xxaa_LLvv += 1.000000 l1_aa(I,i,e) r0(F) t1_aa(f,i)
                    // flops: o0v4L2: 1, o0v2L2: 1, o1v2L1: 1 | mem: o0v4L2: 1, o0v2L2: 1, o0v2L1: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += l1_aa("I,i,e") * t1_aa("f,i") * r0("F");
                }

                if (include_u1_) {

                    // dm_xxaa_LLvv += 1.000000 m1_aa(I,i,e) s1_aa(F,f,i)
                    // flops: o0v4L2: 1, o1v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += m1_aa("I,i,e") * s1_aa("F,f,i");
                }

                {

                    // dm_xxaa_LLvv += 1.000000 l1_aa(I,i,e) r1_aa(F,f,i)
                    // flops: o0v4L2: 1, o1v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += l1_aa("I,i,e") * r1_aa("F,f,i");
                }

                if (include_u1_) {

                    // RDM_blks_["D1_aa_vv"](F) u1_aa(f,i) m1_aa(I,i,e)
                    // flops: o0v4L2: 1, o0v2L2: 1, o1v2L1: 1 | mem: o0v4L2: 1, o0v2L2: 1, o0v2L1: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += u1_aa("f,i") * m1_aa("I,i,e") * r0("F");
                }

                if (include_u0_ && include_u2_) {

                    // dm_xxaa_LLvv += 0.500000 t2_aaaa(f,a,i,j) m2_aaaa(I,i,j,e,a) s0(F)
                    // flops: o0v4L2: 1, o2v3L1: 1, o0v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1, o0v2L1: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += 0.500000 * t2_aaaa("f,a,i,j") * m2_aaaa("I,i,j,e,a") * s0("F");
                }

                if (include_u2_) {

                    // dm_xxaa_LLvv += 0.500000 m2_aaaa(I,j,i,e,a) s2_aaaa(F,f,a,j,i)
                    // flops: o2v3L2: 1, o0v4L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_aa_vv"]("I,F,e,f") += 0.500000 * m2_aaaa("I,j,i,e,a") * s2_aaaa("F,f,a,j,i");

                    // dm_xxbb_LLvv += 0.500000 m2_bbbb(I,j,i,e,a) s2_bbbb(F,f,a,j,i)
                    // flops: o2v3L2: 1, o0v4L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += 0.500000 * m2_bbbb("I,j,i,e,a") * s2_bbbb("F,f,a,j,i");
                }

                {

                    // dm_xxbb_LLvv += 1.000000 l1_bb(I,i,e) r0(F) t1_bb(f,i)
                    // flops: o0v4L2: 1, o0v2L2: 1, o1v2L1: 1 | mem: o0v4L2: 1, o0v2L2: 1, o0v2L1: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += l1_bb("I,i,e") * t1_bb("f,i") * r0("F");
                }

                if (include_u0_ && include_u2_) {

                    // dm_xxbb_LLvv += 0.500000 t2_bbbb(f,a,i,j) m2_bbbb(I,i,j,e,a) s0(F)
                    // flops: o0v4L2: 1, o2v3L1: 1, o0v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1, o0v2L1: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += 0.500000 * t2_bbbb("f,a,i,j") * m2_bbbb("I,i,j,e,a") * s0("F");
                }

                {

                    // dm_xxbb_LLvv += 0.500000 l2_bbbb(I,i,j,e,a) r0(F) t2_bbbb(f,a,i,j)
                    // flops: o0v4L2: 1, o2v3L1: 1, o0v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1, o0v2L1: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += 0.500000 * l2_bbbb("I,i,j,e,a") * t2_bbbb("f,a,i,j") * r0("F");
                }

                if (include_u2_) {

                    // RDM_blks_["D1_bb_vv"](F) u2_bbbb(f,a,i,j) m2_bbbb(I,i,j,e,a)
                    // flops: o0v4L2: 1, o2v3L1: 1, o0v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1, o0v2L1: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += 0.500000 * u2_bbbb("f,a,i,j") * m2_bbbb("I,i,j,e,a") * r0("F");
                }

                if (include_u0_ && include_u1_) {

                    // dm_xxbb_LLvv += 1.000000 t1_bb(f,i) m1_bb(I,i,e) s0(F)
                    // flops: o0v4L2: 1, o0v2L2: 1, o1v2L1: 1 | mem: o0v4L2: 1, o0v2L2: 1, o0v2L1: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += t1_bb("f,i") * m1_bb("I,i,e") * s0("F");
                }

                {

                    // dm_xxbb_LLvv += 0.500000 l2_bbbb(I,j,i,e,a) r2_bbbb(F,f,a,j,i)
                    // flops: o2v3L2: 1, o0v4L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += 0.500000 * l2_bbbb("I,j,i,e,a") * r2_bbbb("F,f,a,j,i");
                }

                if (include_u1_) {

                    // dm_xxbb_LLvv += 1.000000 m1_bb(I,i,e) s1_bb(F,f,i)
                    // flops: o0v4L2: 1, o1v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += m1_bb("I,i,e") * s1_bb("F,f,i");

                    // RDM_blks_["D1_bb_vv"](F) u1_bb(f,i) m1_bb(I,i,e)
                    // flops: o0v4L2: 1, o0v2L2: 1, o1v2L1: 1 | mem: o0v4L2: 1, o0v2L2: 1, o0v2L1: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += u1_bb("f,i") * m1_bb("I,i,e") * r0("F");
                }

                if (include_u2_) {

                    // dm_xxbb_LLvv += 0.500000 m2_abab(I,j,i,a,e) s2_abab(F,a,f,j,i)
                    // +               0.500000 m2_abab(I,i,j,a,e) s2_abab(F,a,f,i,j)
                    // flops: o2v3L2: 1, o0v4L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += m2_abab("I,j,i,a,e") * s2_abab("F,a,f,j,i");
                }

                {

                    // dm_xxbb_LLvv += 1.000000 l1_bb(I,i,e) r1_bb(F,f,i)
                    // flops: o0v4L2: 1, o1v2L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += l1_bb("I,i,e") * r1_bb("F,f,i");

                    // dm_xxbb_LLvv += 0.500000 l2_abab(I,j,i,a,e) r2_abab(F,a,f,j,i)
                    // +               0.500000 l2_abab(I,i,j,a,e) r2_abab(F,a,f,i,j)
                    // flops: o2v3L2: 1, o0v4L2: 1 | mem: o0v4L2: 1, o0v2L2: 1,
                    RDM_blks_["D1_bb_vv"]("I,F,e,f") += l2_abab("I,j,i,a,e") * r2_abab("F,a,f,j,i");

                    // dm_xxaa_LLov += -1.000000 l1_aa(I,i,a) r2_aaaa(F,e,a,i,m)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= l1_aa("I,i,a") * r2_aaaa("F,e,a,i,m");
                }

                if (include_u1_) {

                    // dm_xxaa_LLov += -1.000000 r1_aa(F,e,i) u1_aa(a,m) m1_aa(I,i,a)
                    // flops: o1v3L2: 1, o2v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= u1_aa("a,m") * m1_aa("I,i,a") * r1_aa("F,e,i");
                }

                if (include_u2_) {

                    // dm_xxaa_LLov += 0.500000 r1_aa(F,a,m) u2_aaaa(e,b,i,j) m2_aaaa(I,i,j,b,a)
                    // flops: o1v3L2: 1, o2v3L1: 1, o1v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o0v2L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += 0.500000 * u2_aaaa("e,b,i,j") * m2_aaaa("I,i,j,b,a") * r1_aa("F,a,m");

                    // RDM_blks_["D1_aa_ov"](F) t1_aa(a,m) u2_aaaa(e,b,i,j) m2_aaaa(I,i,j,b,a)
                    // flops: o1v3L2: 1, o3v2L1: 2, o1v1L2: 1 | mem: o1v3L2: 1, o3v1L1: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += 0.500000 * t1_aa("a,m") * m2_aaaa("I,i,j,b,a") * u2_aaaa("e,b,i,j") * r0("F");
                }

                if (include_u1_) {

                    // RDM_blks_["D1_aa_ov"](F) t1_aa(e,i) u1_aa(a,m) m1_aa(I,i,a)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 2 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1, o2v0L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= u1_aa("a,m") * m1_aa("I,i,a") * t1_aa("e,i") * r0("F");
                }

                if (include_u0_ && include_u1_) {

                    // dm_xxaa_LLov += 1.000000 t2_abab(e,a,m,i) m1_bb(I,i,a) s0(F)
                    // flops: o1v3L2: 1, o2v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += t2_abab("e,a,m,i") * m1_bb("I,i,a") * s0("F");
                }

                if (include_u0_) {

                    // dm_xxaa_LLov += 1.000000 t1_aa(e,m) m0(I) s0(F)
                    // flops: o1v3L2: 1, o1v1L2: 1, o0v0L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o0v0L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += m0("I") * s0("F") * t1_aa("e,m");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxaa_LLov += 0.500000 t2_aaaa(b,a,i,m) m2_aaaa(I,j,i,b,a) s1_aa(F,e,j)
                    // flops: o1v3L2: 1, o3v2L1: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += 0.500000 * t2_aaaa("b,a,i,m") * m2_aaaa("I,j,i,b,a") * s1_aa("F,e,j");
                }

                {

                    // dm_xxaa_LLov += 0.500000 l2_aaaa(I,j,i,b,a) r1_aa(F,e,j) t2_aaaa(b,a,i,m)
                    // flops: o1v3L2: 1, o3v2L1: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += 0.500000 * l2_aaaa("I,j,i,b,a") * t2_aaaa("b,a,i,m") * r1_aa("F,e,j");
                }

                if (include_u0_ && include_u1_) {

                    // dm_xxaa_LLov += -1.000000 t2_aaaa(e,a,i,m) m1_aa(I,i,a) s0(F)
                    // flops: o1v3L2: 1, o2v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= t2_aaaa("e,a,i,m") * m1_aa("I,i,a") * s0("F");
                }

                {

                    // dm_xxaa_LLov += -1.000000 l1_aa(I,i,a) r1_aa(F,e,i) t1_aa(a,m)
                    // flops: o1v3L2: 1, o2v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= l1_aa("I,i,a") * t1_aa("a,m") * r1_aa("F,e,i");
                }

                if (include_u0_ && include_u1_) {

                    // RDM_blks_["D1_aa_ov"](I) s1_aa(F,e,m)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += m0("I") * s1_aa("F,e,m");

                    // RDM_blks_["D1_aa_ov"](F) u1_aa(e,m) m0(I)
                    // flops: o1v3L2: 1, o1v1L2: 1, o0v0L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o0v0L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += r0("F") * m0("I") * u1_aa("e,m");
                }

                if (include_u1_ && include_u2_) {

                    // RDM_blks_["D1_aa_ov"](F) u2_abab(e,a,m,i) m1_bb(I,i,a)
                    // flops: o1v3L2: 1, o2v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += u2_abab("e,a,m,i") * m1_bb("I,i,a") * r0("F");
                }

                {

                    // RDM_blks_["D1_aa_ov"](I) r0(F) t1_aa(e,m)
                    // flops: o1v3L2: 1, o1v1L2: 1, o0v0L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o0v0L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += l0("I") * r0("F") * t1_aa("e,m");
                }

                if (include_u1_) {

                    // dm_xxaa_LLov += -1.000000 t1_aa(e,i) m1_aa(I,i,a) s1_aa(F,a,m)
                    // flops: o1v3L2: 1, o2v1L2: 2 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= m1_aa("I,i,a") * s1_aa("F,a,m") * t1_aa("e,i");

                    // dm_xxaa_LLov += -1.000000 r1_aa(F,a,m) u1_aa(e,i) m1_aa(I,i,a)
                    // flops: o1v3L2: 1, o2v1L2: 2 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= r1_aa("F,a,m") * m1_aa("I,i,a") * u1_aa("e,i");
                }

                {

                    // RDM_blks_["D1_aa_ov"](I) r1_aa(F,e,m)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += l0("I") * r1_aa("F,e,m");

                    // dm_xxaa_LLov += 1.000000 l1_bb(I,i,a) r0(F) t2_abab(e,a,m,i)
                    // flops: o1v3L2: 1, o2v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += l1_bb("I,i,a") * t2_abab("e,a,m,i") * r0("F");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxaa_LLov += -0.500000 r2_abab(F,b,a,m,i) u1_aa(e,j) m2_abab(I,j,i,b,a)
                    // +               -0.500000 r2_abab(F,a,b,m,i) u1_aa(e,j) m2_abab(I,j,i,a,b)
                    // flops: o3v2L2: 1, o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= r2_abab("F,b,a,m,i") * m2_abab("I,j,i,b,a") * u1_aa("e,j");

                    // RDM_blks_["D1_aa_ov"](F) u2_aaaa(e,a,i,m) m1_aa(I,i,a)
                    // flops: o1v3L2: 1, o2v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= u2_aaaa("e,a,i,m") * m1_aa("I,i,a") * r0("F");
                }

                {

                    // dm_xxaa_LLov += 1.000000 l1_bb(I,i,a) r2_abab(F,e,a,m,i)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += l1_bb("I,i,a") * r2_abab("F,e,a,m,i");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxaa_LLov += 0.500000 t2_aaaa(e,a,i,j) m2_aaaa(I,i,j,a,b) s1_aa(F,b,m)
                    // flops: o1v3L2: 1, o2v3L1: 1, o1v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o0v2L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += 0.500000 * t2_aaaa("e,a,i,j") * m2_aaaa("I,i,j,a,b") * s1_aa("F,b,m");

                    // dm_xxaa_LLov += -1.000000 m1_aa(I,i,a) s2_aaaa(F,e,a,i,m)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= m1_aa("I,i,a") * s2_aaaa("F,e,a,i,m");
                }

                {

                    // dm_xxaa_LLov += 0.500000 l2_aaaa(I,i,j,a,b) r1_aa(F,b,m) t2_aaaa(e,a,i,j)
                    // flops: o1v3L2: 1, o2v3L1: 1, o1v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o0v2L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += 0.500000 * l2_aaaa("I,i,j,a,b") * t2_aaaa("e,a,i,j") * r1_aa("F,b,m");

                    // dm_xxaa_LLov += -1.000000 l1_aa(I,i,a) r0(F) t2_aaaa(e,a,i,m)
                    // flops: o1v3L2: 1, o2v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= l1_aa("I,i,a") * t2_aaaa("e,a,i,m") * r0("F");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxaa_LLov += 1.000000 m1_bb(I,i,a) s2_abab(F,e,a,m,i)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") += m1_bb("I,i,a") * s2_abab("F,e,a,m,i");
                }

                {

                    // dm_xxaa_LLov += -1.000000 l1_aa(I,i,a) r1_aa(F,a,m) t1_aa(e,i)
                    // flops: o1v3L2: 1, o2v1L2: 2 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= l1_aa("I,i,a") * r1_aa("F,a,m") * t1_aa("e,i");

                    // dm_xxaa_LLov += -1.000000 l1_aa(I,i,a) r0(F) t1_aa(e,i) t1_aa(a,m)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 2 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1, o2v0L1: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= l1_aa("I,i,a") * t1_aa("a,m") * t1_aa("e,i") * r0("F");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxaa_LLov += -0.500000 r2_aaaa(F,b,a,i,m) u1_aa(e,j) m2_aaaa(I,i,j,b,a)
                    // flops: o3v2L2: 1, o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_aa_ov"]("I,F,m,e") -= 0.500000 * r2_aaaa("F,b,a,i,m") * m2_aaaa("I,i,j,b,a") * u1_aa("e,j");

                    // dm_xxbb_LLov += -1.000000 m1_bb(I,i,a) s2_bbbb(F,e,a,i,m)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= m1_bb("I,i,a") * s2_bbbb("F,e,a,i,m");
                }

                if (include_u0_ && include_u1_) {

                    // RDM_blks_["D1_bb_ov"](F) u1_bb(e,m) m0(I)
                    // flops: o1v3L2: 1, o1v1L2: 1, o0v0L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o0v0L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += r0("F") * m0("I") * u1_bb("e,m");

                    // dm_xxbb_LLov += 1.000000 t2_abab(a,e,i,m) m1_aa(I,i,a) s0(F)
                    // flops: o1v3L2: 1, o2v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += t2_abab("a,e,i,m") * m1_aa("I,i,a") * s0("F");
                }

                if (include_u1_) {

                    // dm_xxbb_LLov += -1.000000 r1_bb(F,e,i) u1_bb(a,m) m1_bb(I,i,a)
                    // flops: o1v3L2: 1, o2v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= u1_bb("a,m") * m1_bb("I,i,a") * r1_bb("F,e,i");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxbb_LLov += 1.000000 m1_aa(I,i,a) s2_abab(F,a,e,i,m)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += m1_aa("I,i,a") * s2_abab("F,a,e,i,m");
                }

                {

                    // RDM_blks_["D1_bb_ov"](I) r1_bb(F,e,m)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += l0("I") * r1_bb("F,e,m");
                }

                if (include_u1_) {

                    // dm_xxbb_LLov += -1.000000 r1_bb(F,a,m) u1_bb(e,i) m1_bb(I,i,a)
                    // flops: o1v3L2: 1, o2v1L2: 2 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= r1_bb("F,a,m") * m1_bb("I,i,a") * u1_bb("e,i");
                }

                {

                    // dm_xxbb_LLov += -1.000000 l1_bb(I,i,a) r0(F) t2_bbbb(e,a,i,m)
                    // flops: o1v3L2: 1, o2v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= l1_bb("I,i,a") * t2_bbbb("e,a,i,m") * r0("F");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxbb_LLov += 0.500000 t2_bbbb(e,a,i,j) m2_bbbb(I,i,j,a,b) s1_bb(F,b,m)
                    // flops: o1v3L2: 1, o2v3L1: 1, o1v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o0v2L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += 0.500000 * t2_bbbb("e,a,i,j") * m2_bbbb("I,i,j,a,b") * s1_bb("F,b,m");
                }

                {

                    // dm_xxbb_LLov += 1.000000 l1_aa(I,i,a) r2_abab(F,a,e,i,m)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += l1_aa("I,i,a") * r2_abab("F,a,e,i,m");

                    // dm_xxbb_LLov += -1.000000 l1_bb(I,i,a) r1_bb(F,e,i) t1_bb(a,m)
                    // flops: o1v3L2: 1, o2v1L2: 1, o2v1L1: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= l1_bb("I,i,a") * t1_bb("a,m") * r1_bb("F,e,i");
                }

                if (include_u1_ && include_u2_) {

                    // RDM_blks_["D1_bb_ov"](F) u2_abab(a,e,i,m) m1_aa(I,i,a)
                    // flops: o1v3L2: 1, o2v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += u2_abab("a,e,i,m") * m1_aa("I,i,a") * r0("F");

                    // dm_xxbb_LLov += -0.500000 r2_abab(F,b,a,i,m) u1_bb(e,j) m2_abab(I,i,j,b,a)
                    // +               -0.500000 r2_abab(F,a,b,i,m) u1_bb(e,j) m2_abab(I,i,j,a,b)
                    // flops: o3v2L2: 1, o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= r2_abab("F,b,a,i,m") * m2_abab("I,i,j,b,a") * u1_bb("e,j");
                }

                if (include_u0_ && include_u1_) {

                    // dm_xxbb_LLov += -1.000000 t2_bbbb(e,a,i,m) m1_bb(I,i,a) s0(F)
                    // flops: o1v3L2: 1, o2v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= t2_bbbb("e,a,i,m") * m1_bb("I,i,a") * s0("F");
                }

                if (include_u0_) {

                    // dm_xxbb_LLov += 1.000000 t1_bb(e,m) m0(I) s0(F)
                    // flops: o1v3L2: 1, o1v1L2: 1, o0v0L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o0v0L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += m0("I") * s0("F") * t1_bb("e,m");
                }

                if (include_u1_ && include_u2_) {

                    // RDM_blks_["D1_bb_ov"](F) u2_bbbb(e,a,i,m) m1_bb(I,i,a)
                    // flops: o1v3L2: 1, o2v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= u2_bbbb("e,a,i,m") * m1_bb("I,i,a") * r0("F");
                }

                {

                    // dm_xxbb_LLov += -1.000000 l1_bb(I,i,a) r1_bb(F,a,m) t1_bb(e,i)
                    // flops: o1v3L2: 1, o2v1L2: 2 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= l1_bb("I,i,a") * r1_bb("F,a,m") * t1_bb("e,i");
                }

                if (include_u1_ && include_u2_) {

                    // dm_xxbb_LLov += 0.500000 t2_bbbb(b,a,i,m) m2_bbbb(I,j,i,b,a) s1_bb(F,e,j)
                    // flops: o1v3L2: 1, o3v2L1: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += 0.500000 * t2_bbbb("b,a,i,m") * m2_bbbb("I,j,i,b,a") * s1_bb("F,e,j");

                    // dm_xxbb_LLov += -0.500000 r2_bbbb(F,b,a,i,m) u1_bb(e,j) m2_bbbb(I,i,j,b,a)
                    // flops: o3v2L2: 1, o1v3L2: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= 0.500000 * r2_bbbb("F,b,a,i,m") * m2_bbbb("I,i,j,b,a") * u1_bb("e,j");
                }

                if (include_u1_) {

                    // dm_xxbb_LLov += -1.000000 t1_bb(e,i) m1_bb(I,i,a) s1_bb(F,a,m)
                    // flops: o1v3L2: 1, o2v1L2: 2 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= m1_bb("I,i,a") * s1_bb("F,a,m") * t1_bb("e,i");
                }

                if (include_u2_) {

                    // RDM_blks_["D1_bb_ov"](F) t1_bb(a,m) u2_bbbb(e,b,i,j) m2_bbbb(I,i,j,b,a)
                    // flops: o1v3L2: 1, o3v2L1: 2, o1v1L2: 1 | mem: o1v3L2: 1, o3v1L1: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += 0.500000 * t1_bb("a,m") * m2_bbbb("I,i,j,b,a") * u2_bbbb("e,b,i,j") * r0("F");
                }

                {

                    // dm_xxbb_LLov += 0.500000 l2_bbbb(I,i,j,a,b) r1_bb(F,b,m) t2_bbbb(e,a,i,j)
                    // flops: o1v3L2: 1, o2v3L1: 1, o1v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o0v2L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += 0.500000 * l2_bbbb("I,i,j,a,b") * t2_bbbb("e,a,i,j") * r1_bb("F,b,m");

                    // dm_xxbb_LLov += 1.000000 l1_aa(I,i,a) r0(F) t2_abab(a,e,i,m)
                    // flops: o1v3L2: 1, o2v2L1: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += l1_aa("I,i,a") * t2_abab("a,e,i,m") * r0("F");

                    // dm_xxbb_LLov += 0.500000 l2_bbbb(I,j,i,b,a) r1_bb(F,e,j) t2_bbbb(b,a,i,m)
                    // flops: o1v3L2: 1, o3v2L1: 1, o2v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o2v0L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += 0.500000 * l2_bbbb("I,j,i,b,a") * t2_bbbb("b,a,i,m") * r1_bb("F,e,j");
                }

                if (include_u0_ && include_u1_) {

                    // RDM_blks_["D1_bb_ov"](I) s1_bb(F,e,m)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += m0("I") * s1_bb("F,e,m");
                }

                {

                    // dm_xxbb_LLov += -1.000000 l1_bb(I,i,a) r0(F) t1_bb(e,i) t1_bb(a,m)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 2 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1, o2v0L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= l1_bb("I,i,a") * t1_bb("a,m") * t1_bb("e,i") * r0("F");

                    // RDM_blks_["D1_bb_ov"](I) r0(F) t1_bb(e,m)
                    // flops: o1v3L2: 1, o1v1L2: 1, o0v0L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o0v0L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += l0("I") * r0("F") * t1_bb("e,m");

                    // dm_xxbb_LLov += -1.000000 l1_bb(I,i,a) r2_bbbb(F,e,a,i,m)
                    // flops: o1v3L2: 1, o2v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= l1_bb("I,i,a") * r2_bbbb("F,e,a,i,m");
                }

                if (include_u1_) {

                    // RDM_blks_["D1_bb_ov"](F) t1_bb(e,i) u1_bb(a,m) m1_bb(I,i,a)
                    // flops: o1v3L2: 1, o1v1L2: 1, o2v1L1: 2 | mem: o1v3L2: 1, o1v1L2: 1, o1v1L1: 1, o2v0L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") -= u1_bb("a,m") * m1_bb("I,i,a") * t1_bb("e,i") * r0("F");
                }

                if (include_u2_) {

                    // dm_xxbb_LLov += 0.500000 r1_bb(F,a,m) u2_bbbb(e,b,i,j) m2_bbbb(I,i,j,b,a)
                    // flops: o1v3L2: 1, o2v3L1: 1, o1v2L2: 1 | mem: o1v3L2: 1, o1v1L2: 1, o0v2L1: 1,
                    RDM_blks_["D1_bb_ov"]("I,F,m,e") += 0.500000 * u2_bbbb("e,b,i,j") * m2_bbbb("I,i,j,b,a") * r1_bb("F,a,m");
                }

                if (include_u0_ && include_u1_) {

                    // dm_xxaa_LLvo += 1.000000 m1_aa(I,m,e) s0(F)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_vo"]("I,F,e,m") += m1_aa("I,m,e") * s0("F");
                }

                {

                    // dm_xxaa_LLvo += 1.000000 l1_aa(I,m,e) r0(F)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_aa_vo"]("I,F,e,m") += l1_aa("I,m,e") * r0("F");
                }

                if (include_u0_ && include_u1_) {

                    // dm_xxbb_LLvo += 1.000000 m1_bb(I,m,e) s0(F)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_vo"]("I,F,e,m") += m1_bb("I,m,e") * s0("F");
                }

                {

                    // dm_xxbb_LLvo += 1.000000 l1_bb(I,m,e) r0(F)
                    // flops: o1v3L2: 1, o1v1L2: 1 | mem: o1v3L2: 1, o1v1L2: 1,
                    RDM_blks_["D1_bb_vo"]("I,F,e,m") += l1_bb("I,m,e") * r0("F");
                }

                /*
                              Total Number of Terms: 387
                Number of Flops: (old) 1246 -> (new) 884

                Total FLOP scaling:
                ------------------
                   Scaling :      new |      old |     diff
                  -------- : -------- | -------- | --------
                  o2v3L2 :        8 |       12 |       -4
                  o3v2L2 :       24 |       48 |      -24
                  o0v4L2 :       42 |       52 |      -10
                  o1v3L2 :      212 |      270 |      -58
                  o2v2L2 :      119 |      224 |     -105
                  o2v3L1 :       18 |       36 |      -18
                  o3v1L2 :       28 |       28 |        0
                  o3v2L1 :       48 |      154 |     -106
                  o1v2L2 :       28 |       34 |       -6
                  o2v1L2 :       62 |       74 |      -12
                  o2v2L1 :       12 |       12 |        0
                  o3v1L1 :       12 |        0 |       12
                  o0v2L2 :       18 |       24 |       -6
                  o1v1L2 :      114 |      150 |      -36
                  o1v2L1 :        6 |        6 |        0
                  o2v0L2 :       50 |       60 |      -10
                  o2v1L1 :       38 |       52 |      -14
                  o0v2L1 :        6 |        0 |        6
                  o2v0L1 :       14 |        0 |       14
                  o0v0L2 :       25 |       10 |       15

                Total MEM scaling:
                  o0v4L2 :       42 |       52 |      -10
                  o1v3L2 :      212 |      270 |      -58
                  o2v2L2 :       66 |       88 |      -22
                  o3v1L1 :       26 |       70 |      -44
                  o0v2L2 :       42 |       52 |      -10
                  o1v1L2 :      228 |      306 |      -78
                  o2v0L2 :      100 |      128 |      -28
                  o0v2L1 :       30 |       42 |      -12
                  o1v1L1 :       52 |       68 |      -16
                  o2v0L1 :       46 |       80 |      -34
                  o0v0L2 :       40 |       90 |      -50

                */

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
