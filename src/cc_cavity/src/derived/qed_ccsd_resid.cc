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
 *  You should have renergyeived a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 *  @END LICENSE
 */

#include "../../include/derived/qed_ccsd.h"


using namespace std;
using namespace TA;
using namespace hilbert;

double QED_CCSD::build_residuals() {
    if ( !has_t1_integrals_ ) transform_integrals(true);


    if ( options_.get_bool("ZERO_U0") ) {
        include_u0_ = false;
        scalar_amps_["u0"] = 0.0;
    }

    // Get cavity information
    double w0 = cavity_frequency_[2];
    double coupling_factor_z = w0 * cavity_coupling_strength_[2];

    TA::TArrayD dp_aa_oo, dp_bb_oo, dp_aa_ov, dp_bb_ov, dp_aa_vo, dp_bb_vo, dp_aa_vv, dp_bb_vv;

    // get dipole integrals
    dp_aa_oo("i, j") = coupling_factor_z * Dip_blks_["dz_aa_oo"]("i, j");
    dp_aa_ov("i, a") = coupling_factor_z * Dip_blks_["dz_aa_ov"]("i, a");
    dp_aa_vo("a, i") = coupling_factor_z * Dip_blks_["dz_aa_vo"]("a, i");
    dp_aa_vv("a, b") = coupling_factor_z * Dip_blks_["dz_aa_vv"]("a, b");

    dp_bb_oo("i, j") = coupling_factor_z * Dip_blks_["dz_bb_oo"]("i, j");
    dp_bb_ov("i, a") = coupling_factor_z * Dip_blks_["dz_bb_ov"]("i, a");
    dp_bb_vo("a, i") = coupling_factor_z * Dip_blks_["dz_bb_vo"]("a, i");
    dp_bb_vv("a, b") = coupling_factor_z * Dip_blks_["dz_bb_vv"]("a, b");

    // get residual tensors
    TArrayD &rt1_aa_vo = residuals_["t1_aa"];
    TArrayD &rt1_bb_vo = residuals_["t1_bb"];
    TArrayD &ru1_aa_vo = residuals_["u1_aa"];
    TArrayD &ru1_bb_vo = residuals_["u1_bb"];
    TArrayD &rt2_aaaa_vvoo = residuals_["t2_aaaa"];
    TArrayD &rt2_abab_vvoo = residuals_["t2_abab"];
    TArrayD &rt2_bbbb_vvoo = residuals_["t2_bbbb"];
    TArrayD &ru2_aaaa_vvoo = residuals_["u2_aaaa"];
    TArrayD &ru2_abab_vvoo = residuals_["u2_abab"];
    TArrayD &ru2_bbbb_vvoo = residuals_["u2_bbbb"];
    double ru0 = 0.0;

    // initialize to zero
    double energy = 0;
    zero_tiles(rt1_aa_vo);
    zero_tiles(rt1_bb_vo);
    zero_tiles(rt2_aaaa_vvoo);
    zero_tiles(rt2_abab_vvoo);
    zero_tiles(rt2_bbbb_vvoo);

    if ( include_u1_){
        zero_tiles(ru1_aa_vo);
        zero_tiles(ru1_bb_vo);
    }
    if ( include_u2_){
        zero_tiles(ru2_aaaa_vvoo);
        zero_tiles(ru2_abab_vvoo);
        zero_tiles(ru2_bbbb_vvoo);
    }

    // initialize coherent state amplitudes to zero
    TA::TArrayD crt1_aa_vo = HelperD::makeTensor(world_, {va_, oa_}, true);
    TA::TArrayD crt1_bb_vo = HelperD::makeTensor(world_, {vb_, ob_}, true);
    TA::TArrayD crt2_aaaa_vvoo = HelperD::makeTensor(world_, {va_, va_, oa_, oa_}, true);
    TA::TArrayD crt2_abab_vvoo = HelperD::makeTensor(world_, {va_, vb_, oa_, ob_}, true);
    TA::TArrayD crt2_bbbb_vvoo = HelperD::makeTensor(world_, {vb_, vb_, ob_, ob_}, true);

    double cenergy = 0.0;
    double cru0 = 0.0;
    TA::TArrayD cru1_aa_vo;
    TA::TArrayD cru1_bb_vo;
    TA::TArrayD cru2_aaaa_vvoo;
    TA::TArrayD cru2_abab_vvoo;
    TA::TArrayD cru2_bbbb_vvoo;

    if (include_u1_){
        cru1_aa_vo = HelperD::makeTensor(world_, {va_, oa_}, true);
        cru1_bb_vo = HelperD::makeTensor(world_, {vb_, ob_}, true);
    }
    if (include_u2_){
        cru2_aaaa_vvoo = HelperD::makeTensor(world_, {va_, va_, oa_, oa_}, true);
        cru2_abab_vvoo = HelperD::makeTensor(world_, {va_, vb_, oa_, ob_}, true);
        cru2_bbbb_vvoo = HelperD::makeTensor(world_, {vb_, vb_, ob_, ob_}, true);
    }

    /// get reference amplitudes

    // extract scalar amplitude
    double& u0 = scalar_amps_["u0"];
    if(include_u0_) {
        world_.gop.fence();
        HelperD::forall(amplitudes_["u0"], [&u0](auto &tile, auto &x){
            u0 = tile[x];
        });
        world_.gop.fence();
    }

    // extract 1-body amplitudes
    TA::TArrayD &t1_aa_vo = amplitudes_["t1_aa"];
    TA::TArrayD &t1_bb_vo = amplitudes_["t1_bb"];
    TA::TArrayD &u1_aa_vo = amplitudes_["u1_aa"];
    TA::TArrayD &u1_bb_vo = amplitudes_["u1_bb"];

    // extract 2-body amplitudes
    TA::TArrayD &t2_aaaa_vvoo = amplitudes_["t2_aaaa"];
    TA::TArrayD &t2_abab_vvoo = amplitudes_["t2_abab"];
    TA::TArrayD &t2_bbbb_vvoo = amplitudes_["t2_bbbb"];
    TA::TArrayD &u2_aaaa_vvoo = amplitudes_["u2_aaaa"];
    TA::TArrayD &u2_abab_vvoo = amplitudes_["u2_abab"];
    TA::TArrayD &u2_bbbb_vvoo = amplitudes_["u2_bbbb"];

    world_.gop.fence();

    // compute energy and residuals
    {

        TA::TArrayD tempPerm_aaaa_vvoo;
        TA::TArrayD tempPerm_bbbb_vvoo;
        double scalar0;
        double scalar1;
        double scalar2;
        double scalar3;
        double scalar4;
        double scalar5;
        double scalar6;
        double scalar7;
        double scalar8;
        double scalar9;
        double scalar10;
        double scalar11;
        double scalar12;
        double scalar13;
        double scalar14;
        double scalar15;
        double scalar16;
        double scalar17;
        double scalar18;
        scalar0 = dot(F_blks_["aa_oo"]("o0,o1"), Id_blks_["aa_oo"]("o0,o1"));
        scalar1 = dot(V_blks_["aaaa_oovv"]("j,i,a,b"), t2_aaaa_vvoo("a,b,j,i"));
        scalar2 = dot(V_blks_["aaaa_oooo"]("o0,o1,o2,o3"), Id_blks_["aaaa_oooo"]("o0,o1,o2,o3"));
        scalar3 = dot(dp_aa_oo("o0,o1"), Id_blks_["aa_oo"]("o0,o1"));
        scalar4 = dot(V_blks_["bbbb_oooo"]("o0,o1,o2,o3"), Id_blks_["bbbb_oooo"]("o0,o1,o2,o3"));
        scalar5 = dot(dp_aa_ov("i,a"), u1_aa_vo("a,i"));
        scalar6 = dot(V_blks_["abab_oovv"]("i,j,a,b"), t2_abab_vvoo("a,b,i,j"));
        scalar7 = dot(F_blks_["bb_oo"]("o0,o1"), Id_blks_["bb_oo"]("o0,o1"));
        scalar8 = dot(V_blks_["abab_oooo"]("o0,o1,o2,o3"), Id_blks_["abab_oooo"]("o0,o1,o2,o3"));
        scalar9 = dot(dp_bb_ov("i,a"), u1_bb_vo("a,i"));
        scalar10 = dot(dp_bb_oo("o0,o1"), Id_blks_["bb_oo"]("o0,o1"));
        scalar11 = dot(V_blks_["bbbb_oovv"]("j,i,a,b"), t2_bbbb_vvoo("a,b,j,i"));
        scalar12 = dot(F_blks_["bb_ov"]("i,a"), u1_bb_vo("a,i"));
        scalar13 = dot(V_blks_["aaaa_oovv"]("j,i,a,b"), u2_aaaa_vvoo("a,b,j,i"));
        scalar14 = dot(F_blks_["aa_ov"]("i,a"), u1_aa_vo("a,i"));
        scalar15 = dot(dp_bb_oo("o0,o1"), Id_blks_["bb_oo"]("o0,o1"));
        scalar16 = dot(V_blks_["abab_oovv"]("j,i,b,a"), u2_abab_vvoo("b,a,j,i"));
        scalar17 = dot(V_blks_["bbbb_oovv"]("j,i,a,b"), u2_bbbb_vvoo("a,b,j,i"));
        scalar18 = dot(dp_aa_oo("o0,o1"), Id_blks_["aa_oo"]("o0,o1"));
        auto tempArray = vector<TA::TArrayD>(71);



/// ****** p†q ****** ///



        {

            // tempArray[0] += 1.000000 V_blks_["abab_oovv"]("i,j,a,b") t2_abab_vvoo("f,b,m,j")
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[0]("a,f,i,m") = V_blks_["abab_oovv"]("i,j,a,b") * t2_abab_vvoo("f,b,m,j");

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) <j,i||b,a>_abab t2_abab(e,a,n,i) t2_aaaa(b,f,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempArray[0]("b,e,j,n") * t2_aaaa_vvoo("b,f,m,j");

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) <i,j||a,b>_abab t2_aaaa(a,e,n,i) t2_abab(f,b,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += tempArray[0]("a,f,i,m") * t2_aaaa_vvoo("a,e,n,i");

            // rt2_abab_vvoo += 1.000000 <j,i||b,a>_abab t2_abab(e,a,m,i) t2_abab(b,f,j,n)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += tempArray[0]("b,e,j,m") * t2_abab_vvoo("b,f,j,n");

            // rt2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <j,i||b,a>_abab t2_abab(e,a,n,i) u2_aaaa(b,f,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempArray[0]("b,e,j,n") * u2_aaaa_vvoo("b,f,m,j");

            // ru2_abab_vvoo += 1.000000 <j,i||b,a>_abab t2_abab(e,a,m,i) u2_abab(b,f,j,n)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[0]("b,e,j,m") * u2_abab_vvoo("b,f,j,n");
        }
        world_.gop.fence(); tempArray[0].~TArrayD();

        {

            // tempArray[1] += 1.000000 V_blks_["abab_oovv"]("i,j,a,b") t2_bbbb_vvoo("b,f,n,j")
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[1]("a,f,i,n") = V_blks_["abab_oovv"]("i,j,a,b") * t2_bbbb_vvoo("b,f,n,j");

            // rt2_abab_vvoo += 1.000000 <i,j||a,b>_abab t2_aaaa(a,e,m,i) t2_bbbb(b,f,n,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += tempArray[1]("a,f,i,n") * t2_aaaa_vvoo("a,e,m,i");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) <i,j||a,b>_abab t2_abab(a,e,i,n) t2_bbbb(b,f,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempArray[1]("a,f,i,m") * t2_abab_vvoo("a,e,i,n");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) <j,i||b,a>_abab t2_bbbb(a,e,n,i) t2_abab(b,f,j,m)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += tempArray[1]("b,e,j,n") * t2_abab_vvoo("b,f,j,m");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += 1.000000 <j,i||b,a>_abab t2_bbbb(a,f,n,i) u2_aaaa(b,e,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[1]("b,f,j,n") * u2_aaaa_vvoo("b,e,m,j");
        }

        {

            // rt2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <j,i||b,a>_abab t2_bbbb(a,e,n,i) u2_abab(b,f,j,m)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempArray[1]("b,e,j,n") * u2_abab_vvoo("b,f,j,m");
        }
        world_.gop.fence(); tempArray[1].~TArrayD();

        {

            // tempArray[2] += 1.000000 t2_aaaa_vvoo("a,e,n,i") V_blks_["aaaa_oovv"]("j,i,a,b")
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[2]("e,b,n,j") = t2_aaaa_vvoo("a,e,n,i") * V_blks_["aaaa_oovv"]("j,i,a,b");
        }

        if (include_u2_) {

            // ru2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("f,e,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) <j,i||a,b>_aaaa t2_aaaa(a,e,n,i) t2_aaaa(b,f,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempArray[2]("e,b,n,j") * t2_aaaa_vvoo("b,f,m,j");

            // rt2_abab_vvoo += 1.000000 <j,i||a,b>_aaaa t2_aaaa(a,e,m,i) t2_abab(b,f,j,n)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += tempArray[2]("e,b,m,j") * t2_abab_vvoo("b,f,j,n");

            // rt2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <j,i||a,b>_aaaa t2_aaaa(a,e,n,i) u2_aaaa(b,f,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempArray[2]("e,b,n,j") * u2_aaaa_vvoo("b,f,m,j");

            // ru2_abab_vvoo += 1.000000 <j,i||a,b>_aaaa t2_aaaa(a,e,m,i) u2_abab(b,f,j,n)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[2]("e,b,m,j") * u2_abab_vvoo("b,f,j,n");
        }
        world_.gop.fence(); tempArray[2].~TArrayD();

        {

            // tempArray[3] += 1.000000 t2_abab_vvoo("e,a,n,i") V_blks_["bbbb_oovv"]("j,i,a,b")
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[3]("e,b,n,j") = t2_abab_vvoo("e,a,n,i") * V_blks_["bbbb_oovv"]("j,i,a,b");
        }

        if (include_u2_) {

            // ru2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("f,e,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) <j,i||a,b>_bbbb t2_abab(e,a,n,i) t2_abab(f,b,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempArray[3]("e,b,n,j") * t2_abab_vvoo("f,b,m,j");

            // rt2_abab_vvoo += 1.000000 <j,i||a,b>_bbbb t2_abab(e,a,m,i) t2_bbbb(b,f,n,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += tempArray[3]("e,b,m,j") * t2_bbbb_vvoo("b,f,n,j");

            // rt2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <j,i||a,b>_bbbb t2_abab(e,a,n,i) u2_abab(f,b,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempArray[3]("e,b,n,j") * u2_abab_vvoo("f,b,m,j");

            // ru2_abab_vvoo += 1.000000 <j,i||a,b>_bbbb t2_abab(e,a,m,i) u2_bbbb(b,f,n,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[3]("e,b,m,j") * u2_bbbb_vvoo("b,f,n,j");
        }
        world_.gop.fence(); tempArray[3].~TArrayD();

        {

            // tempArray[4] += 1.000000 t2_bbbb_vvoo("a,e,n,i") V_blks_["bbbb_oovv"]("j,i,a,b")
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[4]("e,b,n,j") = t2_bbbb_vvoo("a,e,n,i") * V_blks_["bbbb_oovv"]("j,i,a,b");
        }

        if (include_u2_) {

            // ru2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("f,e,n,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) <j,i||a,b>_bbbb t2_bbbb(a,e,n,i) t2_bbbb(b,f,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempArray[4]("e,b,n,j") * t2_bbbb_vvoo("b,f,m,j");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += 1.000000 <j,i||a,b>_bbbb t2_bbbb(a,f,n,i) u2_abab(e,b,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[4]("f,b,n,j") * u2_abab_vvoo("e,b,m,j");
        }

        {

            // rt2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <j,i||a,b>_bbbb t2_bbbb(a,e,n,i) u2_bbbb(b,f,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempArray[4]("e,b,n,j") * u2_bbbb_vvoo("b,f,m,j");
        }
        world_.gop.fence(); tempArray[4].~TArrayD();

        {

            // tempArray[5] += 1.000000 V_blks_["aaaa_oovv"]("j,i,a,b") t2_abab_vvoo("a,e,i,n")
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[5]("b,e,j,n") = V_blks_["aaaa_oovv"]("j,i,a,b") * t2_abab_vvoo("a,e,i,n");
        }

        if (include_u2_) {

            // ru2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("f,e,n,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) <j,i||a,b>_aaaa t2_abab(a,e,i,n) t2_abab(b,f,j,m)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempArray[5]("b,e,j,n") * t2_abab_vvoo("b,f,j,m");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += 1.000000 <j,i||a,b>_aaaa t2_abab(a,f,i,n) u2_aaaa(b,e,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[5]("b,f,j,n") * u2_aaaa_vvoo("b,e,m,j");
        }

        {

            // rt2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <j,i||a,b>_aaaa t2_abab(a,e,i,n) u2_abab(b,f,j,m)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempArray[5]("b,e,j,n") * u2_abab_vvoo("b,f,j,m");
        }
        world_.gop.fence(); tempArray[5].~TArrayD();

        {

            // tempArray[6] += 1.000000 V_blks_["abab_oovv"]("i,j,a,b") t2_aaaa_vvoo("a,e,n,i")
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[6]("e,b,n,j") = V_blks_["abab_oovv"]("i,j,a,b") * t2_aaaa_vvoo("a,e,n,i");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <i,j||a,b>_abab t2_aaaa(a,e,n,i) u2_abab(f,b,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += tempArray[6]("e,b,n,j") * u2_abab_vvoo("f,b,m,j");

            // ru2_abab_vvoo += 1.000000 <i,j||a,b>_abab t2_aaaa(a,e,m,i) u2_bbbb(b,f,n,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[6]("e,b,m,j") * u2_bbbb_vvoo("b,f,n,j");
        }
        world_.gop.fence(); tempArray[6].~TArrayD();

        {

            // tempArray[7] += 1.000000 V_blks_["abab_oovv"]("i,j,b,a") t2_abab_vvoo("e,a,i,n")
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[7]("b,e,j,n") = V_blks_["abab_oovv"]("i,j,b,a") * t2_abab_vvoo("e,a,i,n");

            // rt2_abab_vvoo += 1.000000 <i,j||b,a>_abab t2_abab(e,a,i,n) t2_abab(b,f,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += tempArray[7]("b,e,j,n") * t2_abab_vvoo("b,f,m,j");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += 1.000000 <i,j||b,a>_abab t2_abab(e,a,i,n) u2_abab(b,f,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[7]("b,e,j,n") * u2_abab_vvoo("b,f,m,j");
        }
        world_.gop.fence(); tempArray[7].~TArrayD();

        {

            // tempArray[8] += 1.000000 V_blks_["abab_oovv"]("j,i,b,a") t2_abab_vvoo("b,f,j,n")
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[8]("a,f,i,n") = V_blks_["abab_oovv"]("j,i,b,a") * t2_abab_vvoo("b,f,j,n");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += 1.000000 <i,j||a,b>_abab t2_abab(a,f,i,n) u2_abab(e,b,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[8]("b,f,j,n") * u2_abab_vvoo("e,b,m,j");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <i,j||a,b>_abab t2_abab(a,e,i,n) u2_bbbb(b,f,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += tempArray[8]("b,e,j,n") * u2_bbbb_vvoo("b,f,m,j");
        }
        world_.gop.fence(); tempArray[8].~TArrayD();

        {

            // tempArray[9] += 0.250000 t2_aaaa_vvoo("a,b,m,n") V_blks_["aaaa_oovv"]("j,i,a,b")
            // flops: o4v2: 1, o4v0: 1 | mem: o4v0: 2,
            tempArray[9]("m,n,j,i") = 0.250000 * t2_aaaa_vvoo("a,b,m,n") * V_blks_["aaaa_oovv"]("j,i,a,b");

            // rt2_aaaa_vvoo += 0.250000 <j,i||a,b>_aaaa t2_aaaa(a,b,m,n) t2_aaaa(e,f,j,i)
            // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_aaaa_vvoo("e,f,m,n") += tempArray[9]("m,n,j,i") * t2_aaaa_vvoo("e,f,j,i");
        }

        if (include_u2_) {

            // ru2_aaaa_vvoo += 0.250000 <j,i||a,b>_aaaa t2_aaaa(a,b,m,n) u2_aaaa(e,f,j,i)
            // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_aaaa_vvoo("e,f,m,n") += tempArray[9]("m,n,j,i") * u2_aaaa_vvoo("e,f,j,i");
        }
        world_.gop.fence(); tempArray[9].~TArrayD();

        {

            // tempArray[10] += 0.250000 V_blks_["bbbb_oovv"]("j,i,a,b") t2_bbbb_vvoo("a,b,m,n")
            // flops: o4v2: 1, o4v0: 1 | mem: o4v0: 2,
            tempArray[10]("j,i,m,n") = 0.250000 * V_blks_["bbbb_oovv"]("j,i,a,b") * t2_bbbb_vvoo("a,b,m,n");

            // rt2_bbbb_vvoo += 0.250000 <j,i||a,b>_bbbb t2_bbbb(a,b,m,n) t2_bbbb(e,f,j,i)
            // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_bbbb_vvoo("e,f,m,n") += tempArray[10]("j,i,m,n") * t2_bbbb_vvoo("e,f,j,i");
        }

        if (include_u2_) {

            // ru2_bbbb_vvoo += 0.250000 <j,i||a,b>_bbbb t2_bbbb(a,b,m,n) u2_bbbb(e,f,j,i)
            // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_bbbb_vvoo("e,f,m,n") += tempArray[10]("j,i,m,n") * u2_bbbb_vvoo("e,f,j,i");
        }
        world_.gop.fence(); tempArray[10].~TArrayD();

        {

            // tempArray[11] += 1.000000 V_blks_["abab_oovv"]("i,j,a,b") t2_abab_vvoo("a,b,m,n")
            // flops: o4v2: 1, o4v0: 1 | mem: o4v0: 2,
            tempArray[11]("i,m,j,n") = V_blks_["abab_oovv"]("i,j,a,b") * t2_abab_vvoo("a,b,m,n");

            // rt2_abab_vvoo += 0.250000 <j,i||a,b>_abab t2_abab(a,b,m,n) t2_abab(e,f,j,i)
            // +                0.250000 <i,j||a,b>_abab t2_abab(a,b,m,n) t2_abab(e,f,i,j)
            // +                0.250000 <j,i||b,a>_abab t2_abab(b,a,m,n) t2_abab(e,f,j,i)
            // +                0.250000 <i,j||b,a>_abab t2_abab(b,a,m,n) t2_abab(e,f,i,j)
            // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += tempArray[11]("j,m,i,n") * t2_abab_vvoo("e,f,j,i");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += 0.250000 <j,i||a,b>_abab t2_abab(a,b,m,n) u2_abab(e,f,j,i)
            // +                0.250000 <i,j||a,b>_abab t2_abab(a,b,m,n) u2_abab(e,f,i,j)
            // +                0.250000 <j,i||b,a>_abab t2_abab(b,a,m,n) u2_abab(e,f,j,i)
            // +                0.250000 <i,j||b,a>_abab t2_abab(b,a,m,n) u2_abab(e,f,i,j)
            // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[11]("j,m,i,n") * u2_abab_vvoo("e,f,j,i");
        }
        world_.gop.fence(); tempArray[11].~TArrayD();

        {

            // tempArray[12] += 1.000000 V_blks_["abab_oovv"]("i,j,a,b") t2_abab_vvoo("a,e,i,j") t2_bbbb_vvoo("b,f,m,n")
            // flops: o2v3: 2, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
            tempArray[12]("e,f,m,n") = V_blks_["abab_oovv"]("i,j,a,b") * t2_abab_vvoo("a,e,i,j") * t2_bbbb_vvoo("b,f,m,n");

            // rt2_bbbb_vvoo += 0.500000 <j,i||b,a>_abab t2_bbbb(a,e,m,n) t2_abab(b,f,j,i)
            // +                0.500000 <i,j||b,a>_abab t2_bbbb(a,e,m,n) t2_abab(b,f,i,j)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += tempArray[12]("f,e,m,n");

            // rt2_bbbb_vvoo += -0.500000 <j,i||a,b>_abab t2_abab(a,e,j,i) t2_bbbb(b,f,m,n)
            // +                -0.500000 <i,j||a,b>_abab t2_abab(a,e,i,j) t2_bbbb(b,f,m,n)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") -= tempArray[12]("e,f,m,n");

            // tempArray[13] += 1.000000 V_blks_["abab_oovv"]("j,i,b,a") t2_abab_vvoo("e,a,j,i") t2_aaaa_vvoo("b,f,m,n")
            // flops: o2v3: 2, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
            tempArray[13]("e,f,m,n") = V_blks_["abab_oovv"]("j,i,b,a") * t2_abab_vvoo("e,a,j,i") * t2_aaaa_vvoo("b,f,m,n");

            // rt2_aaaa_vvoo += -0.500000 <j,i||b,a>_abab t2_abab(e,a,j,i) t2_aaaa(b,f,m,n)
            // +                -0.500000 <i,j||b,a>_abab t2_abab(e,a,i,j) t2_aaaa(b,f,m,n)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") -= tempArray[13]("e,f,m,n");

            // rt2_aaaa_vvoo += 0.500000 <j,i||a,b>_abab t2_aaaa(a,e,m,n) t2_abab(f,b,j,i)
            // +                0.500000 <i,j||a,b>_abab t2_aaaa(a,e,m,n) t2_abab(f,b,i,j)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += tempArray[13]("f,e,m,n");

            // tempArray[14] += 0.500000 t2_aaaa_vvoo("a,e,j,i") V_blks_["aaaa_oovv"]("j,i,a,b")
            // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
            tempArray[14]("e,b") = 0.500000 * t2_aaaa_vvoo("a,e,j,i") * V_blks_["aaaa_oovv"]("j,i,a,b");

            // rt2_aaaa_vvoo += -0.500000 <j,i||a,b>_aaaa t2_aaaa(a,e,j,i) t2_aaaa(b,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_aaaa_vvoo("e,f,m,n") -= tempArray[14]("e,b") * t2_aaaa_vvoo("b,f,m,n");

            // rt2_abab_vvoo += -0.500000 <j,i||a,b>_aaaa t2_aaaa(a,e,j,i) t2_abab(b,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= tempArray[14]("e,b") * t2_abab_vvoo("b,f,m,n");
        }
        world_.gop.fence(); tempArray[12].~TArrayD();
        world_.gop.fence(); tempArray[13].~TArrayD();

        if (include_u2_) {

            // ru2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("f,e,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += -0.500000 P(e,f) <j,i||a,b>_aaaa t2_aaaa(a,e,j,i) u2_aaaa(b,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempArray[14]("e,b") * u2_aaaa_vvoo("b,f,m,n");

            // ru2_abab_vvoo += -0.500000 <j,i||a,b>_aaaa t2_aaaa(a,e,j,i) u2_abab(b,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[14]("e,b") * u2_abab_vvoo("b,f,m,n");
        }
        world_.gop.fence(); tempArray[14].~TArrayD();

        {

            // tempArray[15] += 0.500000 t2_bbbb_vvoo("a,e,j,i") V_blks_["bbbb_oovv"]("j,i,a,b")
            // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
            tempArray[15]("e,b") = 0.500000 * t2_bbbb_vvoo("a,e,j,i") * V_blks_["bbbb_oovv"]("j,i,a,b");

            // rt2_bbbb_vvoo += -0.500000 <j,i||a,b>_bbbb t2_bbbb(a,e,j,i) t2_bbbb(b,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_bbbb_vvoo("e,f,m,n") -= tempArray[15]("e,b") * t2_bbbb_vvoo("b,f,m,n");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += -0.500000 <j,i||a,b>_bbbb t2_bbbb(a,f,j,i) u2_abab(e,b,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[15]("f,b") * u2_abab_vvoo("e,b,m,n");

            // ru2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("f,e,n,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += -0.500000 P(e,f) <j,i||a,b>_bbbb t2_bbbb(a,e,j,i) u2_bbbb(b,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempArray[15]("e,b") * u2_bbbb_vvoo("b,f,m,n");
        }
        world_.gop.fence(); tempArray[15].~TArrayD();

        {

            // tempArray[16] += 1.000000 V_blks_["abab_oovv"]("j,i,b,a") t2_abab_vvoo("e,a,j,i")
            // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
            tempArray[16]("b,e") = V_blks_["abab_oovv"]("j,i,b,a") * t2_abab_vvoo("e,a,j,i");

            // rt2_abab_vvoo += -0.500000 <j,i||b,a>_abab t2_abab(e,a,j,i) t2_abab(b,f,m,n)
            // +                -0.500000 <i,j||b,a>_abab t2_abab(e,a,i,j) t2_abab(b,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= tempArray[16]("b,e") * t2_abab_vvoo("b,f,m,n");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += -0.500000 P(e,f) <j,i||b,a>_abab t2_abab(e,a,j,i) u2_aaaa(b,f,m,n)
            // +                -0.500000 P(e,f) <i,j||b,a>_abab t2_abab(e,a,i,j) u2_aaaa(b,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= tempArray[16]("b,e") * u2_aaaa_vvoo("b,f,m,n");

            // ru2_abab_vvoo += -0.500000 <j,i||b,a>_abab t2_abab(e,a,j,i) u2_abab(b,f,m,n)
            // +                -0.500000 <i,j||b,a>_abab t2_abab(e,a,i,j) u2_abab(b,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[16]("b,e") * u2_abab_vvoo("b,f,m,n");
        }
        world_.gop.fence(); tempArray[16].~TArrayD();

        {

            // tempArray[17] += 1.000000 t2_abab_vvoo("b,f,j,i") V_blks_["abab_oovv"]("j,i,b,a")
            // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
            tempArray[17]("f,a") = t2_abab_vvoo("b,f,j,i") * V_blks_["abab_oovv"]("j,i,b,a");

            // rt2_abab_vvoo += -0.500000 <j,i||b,a>_abab t2_abab(e,a,m,n) t2_abab(b,f,j,i)
            // +                -0.500000 <i,j||b,a>_abab t2_abab(e,a,m,n) t2_abab(b,f,i,j)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= tempArray[17]("f,a") * t2_abab_vvoo("e,a,m,n");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += -0.500000 <j,i||a,b>_abab t2_abab(a,f,j,i) u2_abab(e,b,m,n)
            // +                -0.500000 <i,j||a,b>_abab t2_abab(a,f,i,j) u2_abab(e,b,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[17]("f,b") * u2_abab_vvoo("e,b,m,n");

            // tempPerm_bbbb_vvoo += -0.500000 P(e,f) <j,i||a,b>_abab t2_abab(a,e,j,i) u2_bbbb(b,f,m,n)
            // +                -0.500000 P(e,f) <i,j||a,b>_abab t2_abab(a,e,i,j) u2_bbbb(b,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= tempArray[17]("e,b") * u2_bbbb_vvoo("b,f,m,n");

            // tempArray[18] += 1.000000 dp_bb_vv("e,a") u2_bbbb_vvoo("a,f,m,n")
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[18]("e,f,m,n") = dp_bb_vv("e,a") * u2_bbbb_vvoo("a,f,m,n");

            // ru2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
        }
        world_.gop.fence(); tempArray[17].~TArrayD();

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(e,f) d-_bb(e,a) u2_bbbb(a,f,m,n)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempArray[18]("e,f,m,n");
        }

        {

            // rt2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u0_ && include_u2_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(e,f) d-_bb(e,a) u0 u2_bbbb(a,f,m,n)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempArray[18]("e,f,m,n") * u0;
        }
        world_.gop.fence(); tempArray[18].~TArrayD();

        if (include_u2_) {

            // tempArray[19] += 1.000000 dp_aa_vv("e,a") u2_aaaa_vvoo("a,f,m,n")
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[19]("e,f,m,n") = dp_aa_vv("e,a") * u2_aaaa_vvoo("a,f,m,n");

            // ru2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(e,f) d-_aa(e,a) u2_aaaa(a,f,m,n)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempArray[19]("e,f,m,n");
        }

        {

            // rt2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u0_ && include_u2_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(e,f) d-_aa(e,a) u0 u2_aaaa(a,f,m,n)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempArray[19]("e,f,m,n") * u0;
        }
        world_.gop.fence(); tempArray[19].~TArrayD();

        if (include_u2_) {

            // tempArray[20] += 0.500000 V_blks_["aaaa_oovv"]("j,i,a,b") u2_aaaa_vvoo("b,f,j,i")
            // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
            tempArray[20]("a,f") = 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * u2_aaaa_vvoo("b,f,j,i");

            // tempPerm_aaaa_vvoo += -0.500000 P(e,f) <j,i||a,b>_aaaa t2_aaaa(a,e,m,n) u2_aaaa(b,f,j,i)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= tempArray[20]("a,f") * t2_aaaa_vvoo("a,e,m,n");

            // ru2_abab_vvoo += 0.500000 <j,i||a,b>_aaaa t2_abab(a,f,m,n) u2_aaaa(b,e,j,i)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[20]("a,e") * t2_abab_vvoo("a,f,m,n");

            // tempArray[21] += 1.000000 V_blks_["abab_oovv"]("j,i,a,b") u2_abab_vvoo("f,b,j,i")
            // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
            tempArray[21]("a,f") = V_blks_["abab_oovv"]("j,i,a,b") * u2_abab_vvoo("f,b,j,i");

            // tempPerm_aaaa_vvoo += 0.500000 P(e,f) <j,i||a,b>_abab t2_aaaa(a,e,m,n) u2_abab(f,b,j,i)
            // +                0.500000 P(e,f) <i,j||a,b>_abab t2_aaaa(a,e,m,n) u2_abab(f,b,i,j)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += tempArray[21]("a,f") * t2_aaaa_vvoo("a,e,m,n");

            // ru2_abab_vvoo += -0.500000 <j,i||a,b>_abab t2_abab(a,f,m,n) u2_abab(e,b,j,i)
            // +                -0.500000 <i,j||a,b>_abab t2_abab(a,f,m,n) u2_abab(e,b,i,j)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[21]("a,e") * t2_abab_vvoo("a,f,m,n");

            // tempArray[22] += 1.000000 V_blks_["abab_oovv"]("i,j,b,a") u2_abab_vvoo("b,f,i,j")
            // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
            tempArray[22]("a,f") = V_blks_["abab_oovv"]("i,j,b,a") * u2_abab_vvoo("b,f,i,j");

            // ru2_abab_vvoo += -0.500000 <j,i||b,a>_abab t2_abab(e,a,m,n) u2_abab(b,f,j,i)
            // +                -0.500000 <i,j||b,a>_abab t2_abab(e,a,m,n) u2_abab(b,f,i,j)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[22]("a,f") * t2_abab_vvoo("e,a,m,n");

            // tempPerm_bbbb_vvoo += 0.500000 P(e,f) <j,i||b,a>_abab t2_bbbb(a,e,m,n) u2_abab(b,f,j,i)
            // +                0.500000 P(e,f) <i,j||b,a>_abab t2_bbbb(a,e,m,n) u2_abab(b,f,i,j)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += tempArray[22]("a,f") * t2_bbbb_vvoo("a,e,m,n");
        }
        world_.gop.fence(); tempArray[20].~TArrayD();
        world_.gop.fence(); tempArray[21].~TArrayD();
        world_.gop.fence(); tempArray[22].~TArrayD();

        {

            // tempArray[23] += 0.500000 t2_bbbb_vvoo("b,f,j,i") V_blks_["bbbb_oovv"]("j,i,a,b")
            // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
            tempArray[23]("f,a") = 0.500000 * t2_bbbb_vvoo("b,f,j,i") * V_blks_["bbbb_oovv"]("j,i,a,b");

            // rt2_abab_vvoo += 0.500000 <j,i||a,b>_bbbb t2_abab(e,a,m,n) t2_bbbb(b,f,j,i)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += tempArray[23]("f,a") * t2_abab_vvoo("e,a,m,n");

            // rt2_bbbb_vvoo += -0.500000 <j,i||a,b>_bbbb t2_bbbb(a,e,m,n) t2_bbbb(b,f,j,i)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_bbbb_vvoo("e,f,m,n") -= tempArray[23]("f,a") * t2_bbbb_vvoo("a,e,m,n");
        }
        world_.gop.fence(); tempArray[23].~TArrayD();

        if (include_u2_) {

            // tempArray[24] += 0.500000 V_blks_["bbbb_oovv"]("j,i,a,b") u2_bbbb_vvoo("b,f,j,i")
            // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
            tempArray[24]("a,f") = 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * u2_bbbb_vvoo("b,f,j,i");

            // ru2_abab_vvoo += 0.500000 <j,i||a,b>_bbbb t2_abab(e,a,m,n) u2_bbbb(b,f,j,i)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[24]("a,f") * t2_abab_vvoo("e,a,m,n");

            // tempPerm_bbbb_vvoo += -0.500000 P(e,f) <j,i||a,b>_bbbb t2_bbbb(a,e,m,n) u2_bbbb(b,f,j,i)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= tempArray[24]("a,f") * t2_bbbb_vvoo("a,e,m,n");
        }
        world_.gop.fence(); tempArray[24].~TArrayD();

        if (include_u1_) {

            // tempArray[25] += 1.000000 u1_bb_vo("b,n") V_blks_["bbbb_vovv"]("f,i,a,b")
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[25]("f,a,n,i") = u1_bb_vo("b,n") * V_blks_["bbbb_vovv"]("f,i,a,b");
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += 1.000000 <i,f||a,b>_bbbb t2_abab(e,a,m,i) u1_bb(b,n)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[25]("f,a,n,i") * t2_abab_vvoo("e,a,m,i");
        }

        if (include_u2_) {

            // ru2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) <i,e||a,b>_bbbb t2_bbbb(a,f,n,i) u1_bb(b,m)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempArray[25]("e,a,m,i") * t2_bbbb_vvoo("a,f,n,i");

            // tempArray[26] += 1.000000 V_blks_["baab_vovv"]("f,i,a,b") u1_bb_vo("b,n")
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[26]("a,f,i,n") = V_blks_["baab_vovv"]("f,i,a,b") * u1_bb_vo("b,n");
        }
        world_.gop.fence(); tempArray[25].~TArrayD();

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += -1.000000 <i,f||a,b>_abab t2_aaaa(a,e,m,i) u1_bb(b,n)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[26]("a,f,i,n") * t2_aaaa_vvoo("a,e,m,i");
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <i,e||a,b>_abab t2_abab(a,f,i,n) u1_bb(b,m)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= tempArray[26]("a,e,i,m") * t2_abab_vvoo("a,f,i,n");

            // tempArray[27] += 1.000000 V_blks_["aaaa_vovv"]("e,i,a,b") u1_aa_vo("b,m")
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[27]("e,a,i,m") = V_blks_["aaaa_vovv"]("e,i,a,b") * u1_aa_vo("b,m");
        }
        world_.gop.fence(); tempArray[26].~TArrayD();

        if (include_u2_) {

            // ru2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) <i,e||a,b>_aaaa t2_aaaa(a,f,n,i) u1_aa(b,m)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempArray[27]("e,a,i,m") * t2_aaaa_vvoo("a,f,n,i");
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += 1.000000 <i,e||a,b>_aaaa t2_abab(a,f,i,n) u1_aa(b,m)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[27]("e,a,i,m") * t2_abab_vvoo("a,f,i,n");
        }
        world_.gop.fence(); tempArray[27].~TArrayD();

        if (include_u1_) {

            // tempArray[28] += 1.000000 u1_aa_vo("b,m") V_blks_["abab_vovv"]("e,i,b,a")
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[28]("e,a,m,i") = u1_aa_vo("b,m") * V_blks_["abab_vovv"]("e,i,b,a");

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <e,i||b,a>_abab t2_abab(f,a,n,i) u1_aa(b,m)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += tempArray[28]("e,a,m,i") * t2_abab_vvoo("f,a,n,i");
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += -1.000000 <e,i||b,a>_abab t2_bbbb(a,f,n,i) u1_aa(b,m)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[28]("e,a,m,i") * t2_bbbb_vvoo("a,f,n,i");
        }
        world_.gop.fence(); tempArray[28].~TArrayD();

        if (include_u2_) {

            // tempArray[29] += 1.000000 dp_aa_vv("e,a") u2_abab_vvoo("a,f,m,n")
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[29]("e,f,m,n") = dp_aa_vv("e,a") * u2_abab_vvoo("a,f,m,n");

            // rt2_abab_vvoo += -1.000000 d-_aa(e,a) u2_abab(a,f,m,n)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_abab_vvoo("e,f,m,n") -= tempArray[29]("e,f,m,n");
        }

        if (include_u0_ && include_u2_) {

            // ru2_abab_vvoo += -1.000000 d-_aa(e,a) u0 u2_abab(a,f,m,n)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[29]("e,f,m,n") * u0;
        }
        world_.gop.fence(); tempArray[29].~TArrayD();

        if (include_u2_) {

            // tempArray[30] += 1.000000 dp_bb_vv("f,a") u2_abab_vvoo("e,a,m,n")
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[30]("e,f,m,n") = dp_bb_vv("f,a") * u2_abab_vvoo("e,a,m,n");

            // rt2_abab_vvoo += -1.000000 d-_bb(f,a) u2_abab(e,a,m,n)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_abab_vvoo("e,f,m,n") -= tempArray[30]("e,f,m,n");
        }

        if (include_u0_ && include_u2_) {

            // ru2_abab_vvoo += -1.000000 d-_bb(f,a) u0 u2_abab(e,a,m,n)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[30]("e,f,m,n") * u0;
        }
        world_.gop.fence(); tempArray[30].~TArrayD();

        if (include_u1_) {

            // tempArray[31] += 1.000000 dp_aa_ov("i,a") t2_aaaa_vvoo("a,e,m,n") u1_aa_vo("f,i")
            // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            tempArray[31]("e,f,m,n") = dp_aa_ov("i,a") * t2_aaaa_vvoo("a,e,m,n") * u1_aa_vo("f,i");
        }

        if (include_u2_) {

            // ru2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("f,e,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(e,f) d-_aa(i,a) t2_aaaa(a,e,m,n) u1_aa(f,i)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempArray[31]("e,f,m,n");
        }

        {

            // rt2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u0_ && include_u1_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(e,f) d-_aa(i,a) t2_aaaa(a,e,m,n) u0 u1_aa(f,i)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempArray[31]("e,f,m,n") * u0;
        }
        world_.gop.fence(); tempArray[31].~TArrayD();

        if (include_u1_) {

            // tempArray[32] += 1.000000 dp_bb_ov("i,a") t2_bbbb_vvoo("a,e,m,n") u1_bb_vo("f,i")
            // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            tempArray[32]("e,f,m,n") = dp_bb_ov("i,a") * t2_bbbb_vvoo("a,e,m,n") * u1_bb_vo("f,i");
        }

        if (include_u2_) {

            // ru2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("f,e,n,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(e,f) d-_bb(i,a) t2_bbbb(a,e,m,n) u1_bb(f,i)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempArray[32]("e,f,m,n");
        }

        {

            // rt2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u0_ && include_u1_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(e,f) d-_bb(i,a) t2_bbbb(a,e,m,n) u0 u1_bb(f,i)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempArray[32]("e,f,m,n") * u0;
        }
        world_.gop.fence(); tempArray[32].~TArrayD();

        if (include_u1_) {

            // tempArray[33] += 1.000000 dp_aa_ov("i,a") t2_abab_vvoo("a,f,m,n") u1_aa_vo("e,i")
            // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            tempArray[33]("e,f,m,n") = dp_aa_ov("i,a") * t2_abab_vvoo("a,f,m,n") * u1_aa_vo("e,i");

            // rt2_abab_vvoo += 1.000000 d-_aa(i,a) t2_abab(a,f,m,n) u1_aa(e,i)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_abab_vvoo("e,f,m,n") += tempArray[33]("e,f,m,n");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // ru2_abab_vvoo += 1.000000 d-_aa(i,a) t2_abab(a,f,m,n) u0 u1_aa(e,i)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            ru2_abab_vvoo("e,f,m,n") += tempArray[33]("e,f,m,n") * u0;
        }
        world_.gop.fence(); tempArray[33].~TArrayD();

        if (include_u1_) {

            // tempArray[34] += 1.000000 dp_bb_ov("i,a") t2_abab_vvoo("e,a,m,n") u1_bb_vo("f,i")
            // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            tempArray[34]("e,f,m,n") = dp_bb_ov("i,a") * t2_abab_vvoo("e,a,m,n") * u1_bb_vo("f,i");

            // rt2_abab_vvoo += 1.000000 d-_bb(i,a) t2_abab(e,a,m,n) u1_bb(f,i)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_abab_vvoo("e,f,m,n") += tempArray[34]("e,f,m,n");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // ru2_abab_vvoo += 1.000000 d-_bb(i,a) t2_abab(e,a,m,n) u0 u1_bb(f,i)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            ru2_abab_vvoo("e,f,m,n") += tempArray[34]("e,f,m,n") * u0;
        }
        world_.gop.fence(); tempArray[34].~TArrayD();

        if (include_u1_) {

            // tempArray[35] += 1.000000 dp_aa_ov("i,a") u1_aa_vo("a,m") t2_aaaa_vvoo("e,f,n,i")
            // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
            tempArray[35]("e,f,m,n") = dp_aa_ov("i,a") * u1_aa_vo("a,m") * t2_aaaa_vvoo("e,f,n,i");
        }

        if (include_u2_) {

            // ru2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) d-_aa(i,a) t2_aaaa(e,f,n,i) u1_aa(a,m)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempArray[35]("e,f,m,n");
        }

        {

            // rt2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u0_ && include_u1_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) d-_aa(i,a) t2_aaaa(e,f,n,i) u0 u1_aa(a,m)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempArray[35]("e,f,m,n") * u0;
        }
        world_.gop.fence(); tempArray[35].~TArrayD();

        if (include_u1_) {

            // tempArray[36] += 1.000000 dp_bb_ov("i,a") u1_bb_vo("a,m") t2_bbbb_vvoo("e,f,n,i")
            // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
            tempArray[36]("e,f,m,n") = dp_bb_ov("i,a") * u1_bb_vo("a,m") * t2_bbbb_vvoo("e,f,n,i");
        }

        if (include_u2_) {

            // ru2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) d-_bb(i,a) t2_bbbb(e,f,n,i) u1_bb(a,m)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempArray[36]("e,f,m,n");
        }

        {

            // rt2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u0_ && include_u1_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) d-_bb(i,a) t2_bbbb(e,f,n,i) u0 u1_bb(a,m)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempArray[36]("e,f,m,n") * u0;
        }
        world_.gop.fence(); tempArray[36].~TArrayD();

        if (include_u1_) {

            // tempArray[37] += 1.000000 dp_aa_ov("i,a") u1_aa_vo("a,m") t2_abab_vvoo("e,f,i,n")
            // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
            tempArray[37]("e,f,m,n") = dp_aa_ov("i,a") * u1_aa_vo("a,m") * t2_abab_vvoo("e,f,i,n");

            // rt2_abab_vvoo += 1.000000 d-_aa(i,a) t2_abab(e,f,i,n) u1_aa(a,m)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_abab_vvoo("e,f,m,n") += tempArray[37]("e,f,m,n");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // ru2_abab_vvoo += 1.000000 d-_aa(i,a) t2_abab(e,f,i,n) u0 u1_aa(a,m)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            ru2_abab_vvoo("e,f,m,n") += tempArray[37]("e,f,m,n") * u0;
        }
        world_.gop.fence(); tempArray[37].~TArrayD();

        if (include_u1_) {

            // tempArray[38] += 1.000000 dp_bb_ov("i,a") u1_bb_vo("a,n") t2_abab_vvoo("e,f,m,i")
            // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
            tempArray[38]("e,f,m,n") = dp_bb_ov("i,a") * u1_bb_vo("a,n") * t2_abab_vvoo("e,f,m,i");

            // rt2_abab_vvoo += 1.000000 d-_bb(i,a) t2_abab(e,f,m,i) u1_bb(a,n)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_abab_vvoo("e,f,m,n") += tempArray[38]("e,f,m,n");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // ru2_abab_vvoo += 1.000000 d-_bb(i,a) t2_abab(e,f,m,i) u0 u1_bb(a,n)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            ru2_abab_vvoo("e,f,m,n") += tempArray[38]("e,f,m,n") * u0;
        }
        world_.gop.fence(); tempArray[38].~TArrayD();

        {

            // tempArray[39] += 0.500000 V_blks_["bbbb_oovv"]("j,i,a,b") t2_bbbb_vvoo("a,b,n,i")
            // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
            tempArray[39]("j,n") = 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * t2_bbbb_vvoo("a,b,n,i");

            // rt2_abab_vvoo += -0.500000 <j,i||a,b>_bbbb t2_bbbb(a,b,n,i) t2_abab(e,f,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= tempArray[39]("j,n") * t2_abab_vvoo("e,f,m,j");
        }

        if (include_u2_) {

            // ru2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();

            // tempPerm_bbbb_vvoo += -0.500000 P(m,n) <j,i||a,b>_bbbb t2_bbbb(a,b,n,i) t2_bbbb(e,f,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempArray[39]("j,n") * t2_bbbb_vvoo("e,f,m,j");
        }

        if (include_u1_) {

            // ru1_bb_vo += -0.500000 <j,i||a,b>_bbbb t2_bbbb(a,b,m,i) u1_bb(e,j)
            // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_bb_vo("e,m") -= tempArray[39]("j,m") * u1_bb_vo("e,j");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += -0.500000 <j,i||a,b>_bbbb t2_bbbb(a,b,n,i) u2_abab(e,f,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[39]("j,n") * u2_abab_vvoo("e,f,m,j");
        }

        {

            // rt2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += -0.500000 P(m,n) <j,i||a,b>_bbbb t2_bbbb(a,b,n,i) u2_bbbb(e,f,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempArray[39]("j,n") * u2_bbbb_vvoo("e,f,m,j");
        }
        world_.gop.fence(); tempArray[39].~TArrayD();

        {

            // tempArray[40] += 1.000000 t2_abab_vvoo("a,b,n,i") V_blks_["abab_oovv"]("j,i,a,b")
            // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
            tempArray[40]("n,j") = t2_abab_vvoo("a,b,n,i") * V_blks_["abab_oovv"]("j,i,a,b");
        }

        if (include_u2_) {

            // ru2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();

            // tempPerm_aaaa_vvoo += -0.500000 P(m,n) <j,i||a,b>_abab t2_abab(a,b,n,i) t2_aaaa(e,f,m,j)
            // +                -0.500000 P(m,n) <j,i||b,a>_abab t2_abab(b,a,n,i) t2_aaaa(e,f,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempArray[40]("n,j") * t2_aaaa_vvoo("e,f,m,j");

            // rt2_abab_vvoo += -0.500000 <j,i||a,b>_abab t2_abab(a,b,m,i) t2_abab(e,f,j,n)
            // +                -0.500000 <j,i||b,a>_abab t2_abab(b,a,m,i) t2_abab(e,f,j,n)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= tempArray[40]("m,j") * t2_abab_vvoo("e,f,j,n");
        }

        if (include_u1_) {

            // ru1_aa_vo += -0.500000 <j,i||a,b>_abab t2_abab(a,b,m,i) u1_aa(e,j)
            // +            -0.500000 <j,i||b,a>_abab t2_abab(b,a,m,i) u1_aa(e,j)
            // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_aa_vo("e,m") -= tempArray[40]("m,j") * u1_aa_vo("e,j");
        }

        {

            // rt2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += -0.500000 P(m,n) <j,i||a,b>_abab t2_abab(a,b,n,i) u2_aaaa(e,f,m,j)
            // +                -0.500000 P(m,n) <j,i||b,a>_abab t2_abab(b,a,n,i) u2_aaaa(e,f,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempArray[40]("n,j") * u2_aaaa_vvoo("e,f,m,j");

            // ru2_abab_vvoo += -0.500000 <j,i||a,b>_abab t2_abab(a,b,m,i) u2_abab(e,f,j,n)
            // +                -0.500000 <j,i||b,a>_abab t2_abab(b,a,m,i) u2_abab(e,f,j,n)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[40]("m,j") * u2_abab_vvoo("e,f,j,n");
        }
        world_.gop.fence(); tempArray[40].~TArrayD();

        {

            // tempArray[41] += 0.500000 t2_aaaa_vvoo("a,b,n,i") V_blks_["aaaa_oovv"]("j,i,a,b")
            // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
            tempArray[41]("n,j") = 0.500000 * t2_aaaa_vvoo("a,b,n,i") * V_blks_["aaaa_oovv"]("j,i,a,b");
        }

        if (include_u2_) {

            // ru2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();

            // tempPerm_aaaa_vvoo += -0.500000 P(m,n) <j,i||a,b>_aaaa t2_aaaa(a,b,n,i) t2_aaaa(e,f,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempArray[41]("n,j") * t2_aaaa_vvoo("e,f,m,j");

            // rt2_abab_vvoo += -0.500000 <j,i||a,b>_aaaa t2_aaaa(a,b,m,i) t2_abab(e,f,j,n)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= tempArray[41]("m,j") * t2_abab_vvoo("e,f,j,n");
        }

        if (include_u1_) {

            // ru1_aa_vo += -0.500000 <j,i||a,b>_aaaa t2_aaaa(a,b,m,i) u1_aa(e,j)
            // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_aa_vo("e,m") -= tempArray[41]("m,j") * u1_aa_vo("e,j");
        }

        {

            // rt2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += -0.500000 P(m,n) <j,i||a,b>_aaaa t2_aaaa(a,b,n,i) u2_aaaa(e,f,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempArray[41]("n,j") * u2_aaaa_vvoo("e,f,m,j");

            // ru2_abab_vvoo += -0.500000 <j,i||a,b>_aaaa t2_aaaa(a,b,m,i) u2_abab(e,f,j,n)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[41]("m,j") * u2_abab_vvoo("e,f,j,n");
        }
        world_.gop.fence(); tempArray[41].~TArrayD();

        {

            // tempArray[42] += 1.000000 t2_abab_vvoo("a,b,i,n") V_blks_["abab_oovv"]("i,j,a,b")
            // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
            tempArray[42]("n,j") = t2_abab_vvoo("a,b,i,n") * V_blks_["abab_oovv"]("i,j,a,b");

            // rt2_abab_vvoo += -0.500000 <i,j||a,b>_abab t2_abab(a,b,i,n) t2_abab(e,f,m,j)
            // +                -0.500000 <i,j||b,a>_abab t2_abab(b,a,i,n) t2_abab(e,f,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= tempArray[42]("n,j") * t2_abab_vvoo("e,f,m,j");
        }

        if (include_u2_) {

            // ru2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();

            // tempPerm_bbbb_vvoo += -0.500000 P(m,n) <i,j||a,b>_abab t2_abab(a,b,i,n) t2_bbbb(e,f,m,j)
            // +                -0.500000 P(m,n) <i,j||b,a>_abab t2_abab(b,a,i,n) t2_bbbb(e,f,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempArray[42]("n,j") * t2_bbbb_vvoo("e,f,m,j");
        }

        if (include_u1_) {

            // ru1_bb_vo += -0.500000 <i,j||a,b>_abab t2_abab(a,b,i,m) u1_bb(e,j)
            // +            -0.500000 <i,j||b,a>_abab t2_abab(b,a,i,m) u1_bb(e,j)
            // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_bb_vo("e,m") -= tempArray[42]("m,j") * u1_bb_vo("e,j");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += -0.500000 <i,j||a,b>_abab t2_abab(a,b,i,n) u2_abab(e,f,m,j)
            // +                -0.500000 <i,j||b,a>_abab t2_abab(b,a,i,n) u2_abab(e,f,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[42]("n,j") * u2_abab_vvoo("e,f,m,j");
        }

        {

            // rt2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += -0.500000 P(m,n) <i,j||a,b>_abab t2_abab(a,b,i,n) u2_bbbb(e,f,m,j)
            // +                -0.500000 P(m,n) <i,j||b,a>_abab t2_abab(b,a,i,n) u2_bbbb(e,f,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempArray[42]("n,j") * u2_bbbb_vvoo("e,f,m,j");

            // tempArray[43] += 1.000000 u2_bbbb_vvoo("e,f,m,i") dp_bb_oo("i,n")
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[43]("e,f,m,n") = u2_bbbb_vvoo("e,f,m,i") * dp_bb_oo("i,n");

            // ru2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
        }
        world_.gop.fence(); tempArray[42].~TArrayD();

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) d-_bb(i,n) u2_bbbb(e,f,m,i)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempArray[43]("e,f,m,n");
        }

        {

            // rt2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u0_ && include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) d-_bb(i,n) u0 u2_bbbb(e,f,m,i)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempArray[43]("e,f,m,n") * u0;
        }
        world_.gop.fence(); tempArray[43].~TArrayD();

        if (include_u2_) {

            // tempArray[44] += 1.000000 dp_aa_oo("i,n") u2_aaaa_vvoo("e,f,m,i")
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[44]("e,f,n,m") = dp_aa_oo("i,n") * u2_aaaa_vvoo("e,f,m,i");

            // ru2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) d-_aa(i,n) u2_aaaa(e,f,m,i)
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempArray[44]("e,f,n,m");
        }

        {

            // rt2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u0_ && include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) d-_aa(i,n) u0 u2_aaaa(e,f,m,i)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempArray[44]("e,f,n,m") * u0;
        }
        world_.gop.fence(); tempArray[44].~TArrayD();

        if (include_u2_) {

            // tempArray[45] += 1.000000 u2_abab_vvoo("a,b,j,n") V_blks_["abab_oovv"]("j,i,a,b")
            // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
            tempArray[45]("n,i") = u2_abab_vvoo("a,b,j,n") * V_blks_["abab_oovv"]("j,i,a,b");

            // ru2_abab_vvoo += -0.500000 <j,i||a,b>_abab t2_abab(e,f,m,i) u2_abab(a,b,j,n)
            // +                -0.500000 <j,i||b,a>_abab t2_abab(e,f,m,i) u2_abab(b,a,j,n)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[45]("n,i") * t2_abab_vvoo("e,f,m,i");

            // tempPerm_bbbb_vvoo += 0.500000 P(m,n) <j,i||a,b>_abab t2_bbbb(e,f,n,i) u2_abab(a,b,j,m)
            // +                0.500000 P(m,n) <j,i||b,a>_abab t2_bbbb(e,f,n,i) u2_abab(b,a,j,m)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += tempArray[45]("m,i") * t2_bbbb_vvoo("e,f,n,i");

            // tempArray[46] += 1.000000 u2_abab_vvoo("b,a,m,j") V_blks_["abab_oovv"]("i,j,b,a")
            // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
            tempArray[46]("m,i") = u2_abab_vvoo("b,a,m,j") * V_blks_["abab_oovv"]("i,j,b,a");

            // tempPerm_aaaa_vvoo += 0.500000 P(m,n) <i,j||a,b>_abab t2_aaaa(e,f,n,i) u2_abab(a,b,m,j)
            // +                0.500000 P(m,n) <i,j||b,a>_abab t2_aaaa(e,f,n,i) u2_abab(b,a,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += tempArray[46]("m,i") * t2_aaaa_vvoo("e,f,n,i");

            // ru2_abab_vvoo += -0.500000 <i,j||a,b>_abab t2_abab(e,f,i,n) u2_abab(a,b,m,j)
            // +                -0.500000 <i,j||b,a>_abab t2_abab(e,f,i,n) u2_abab(b,a,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[46]("m,i") * t2_abab_vvoo("e,f,i,n");

            // tempArray[47] += 0.500000 V_blks_["aaaa_oovv"]("j,i,a,b") u2_aaaa_vvoo("a,b,m,j")
            // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
            tempArray[47]("i,m") = 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * u2_aaaa_vvoo("a,b,m,j");

            // tempPerm_aaaa_vvoo += -0.500000 P(m,n) <j,i||a,b>_aaaa t2_aaaa(e,f,n,i) u2_aaaa(a,b,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= tempArray[47]("i,m") * t2_aaaa_vvoo("e,f,n,i");

            // ru2_abab_vvoo += 0.500000 <j,i||a,b>_aaaa t2_abab(e,f,i,n) u2_aaaa(a,b,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[47]("i,m") * t2_abab_vvoo("e,f,i,n");

            // tempArray[48] += 0.500000 u2_bbbb_vvoo("a,b,n,j") V_blks_["bbbb_oovv"]("j,i,a,b")
            // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
            tempArray[48]("n,i") = 0.500000 * u2_bbbb_vvoo("a,b,n,j") * V_blks_["bbbb_oovv"]("j,i,a,b");

            // ru2_abab_vvoo += 0.500000 <j,i||a,b>_bbbb t2_abab(e,f,m,i) u2_bbbb(a,b,n,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[48]("n,i") * t2_abab_vvoo("e,f,m,i");

            // tempPerm_bbbb_vvoo += -0.500000 P(m,n) <j,i||a,b>_bbbb t2_bbbb(e,f,n,i) u2_bbbb(a,b,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= tempArray[48]("m,i") * t2_bbbb_vvoo("e,f,n,i");

            // tempArray[49] += 1.000000 dp_bb_oo("i,n") u2_abab_vvoo("e,f,m,i")
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[49]("e,f,m,n") = dp_bb_oo("i,n") * u2_abab_vvoo("e,f,m,i");

            // rt2_abab_vvoo += 1.000000 d-_bb(i,n) u2_abab(e,f,m,i)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_abab_vvoo("e,f,m,n") += tempArray[49]("e,f,m,n");
        }
        world_.gop.fence(); tempArray[45].~TArrayD();
        world_.gop.fence(); tempArray[46].~TArrayD();
        world_.gop.fence(); tempArray[47].~TArrayD();
        world_.gop.fence(); tempArray[48].~TArrayD();

        if (include_u0_ && include_u2_) {

            // ru2_abab_vvoo += 1.000000 d-_bb(i,n) u0 u2_abab(e,f,m,i)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            ru2_abab_vvoo("e,f,m,n") += tempArray[49]("e,f,m,n") * u0;
        }
        world_.gop.fence(); tempArray[49].~TArrayD();

        if (include_u2_) {

            // tempArray[50] += 1.000000 dp_aa_oo("i,m") u2_abab_vvoo("e,f,i,n")
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[50]("e,f,m,n") = dp_aa_oo("i,m") * u2_abab_vvoo("e,f,i,n");

            // rt2_abab_vvoo += 1.000000 d-_aa(i,m) u2_abab(e,f,i,n)
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_abab_vvoo("e,f,m,n") += tempArray[50]("e,f,m,n");
        }

        if (include_u0_ && include_u2_) {

            // ru2_abab_vvoo += 1.000000 d-_aa(i,m) u0 u2_abab(e,f,i,n)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            ru2_abab_vvoo("e,f,m,n") += tempArray[50]("e,f,m,n") * u0;
        }
        world_.gop.fence(); tempArray[50].~TArrayD();

        if (include_u1_) {

            // tempArray[51] += 1.000000 V_blks_["baab_vovv"]("f,i,b,a") u1_aa_vo("b,i")
            // flops: o1v3: 1, o0v2: 1 | mem: o0v2: 2,
            tempArray[51]("f,a") = V_blks_["baab_vovv"]("f,i,b,a") * u1_aa_vo("b,i");
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += 1.000000 <i,f||b,a>_abab t2_abab(e,a,m,n) u1_aa(b,i)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[51]("f,a") * t2_abab_vvoo("e,a,m,n");
        }

        if (include_u2_) {

            // ru2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(e,f) <i,e||b,a>_abab t2_bbbb(a,f,m,n) u1_aa(b,i)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * tempArray[51]("e,a") * t2_bbbb_vvoo("a,f,m,n");

            // tempArray[52] += 1.000000 u1_bb_vo("b,i") V_blks_["bbbb_vovv"]("f,i,a,b")
            // flops: o1v3: 1, o0v2: 1 | mem: o0v2: 2,
            tempArray[52]("f,a") = u1_bb_vo("b,i") * V_blks_["bbbb_vovv"]("f,i,a,b");
        }
        world_.gop.fence(); tempArray[51].~TArrayD();

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += -1.000000 <i,f||a,b>_bbbb t2_abab(e,a,m,n) u1_bb(b,i)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[52]("f,a") * t2_abab_vvoo("e,a,m,n");
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(e,f) <i,e||a,b>_bbbb t2_bbbb(a,f,m,n) u1_bb(b,i)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += tempArray[52]("e,a") * t2_bbbb_vvoo("a,f,m,n");

            // tempArray[53] += 1.000000 u1_aa_vo("b,i") V_blks_["aaaa_vovv"]("e,i,a,b")
            // flops: o1v3: 1, o0v2: 1 | mem: o0v2: 2,
            tempArray[53]("e,a") = u1_aa_vo("b,i") * V_blks_["aaaa_vovv"]("e,i,a,b");
        }
        world_.gop.fence(); tempArray[52].~TArrayD();

        if (include_u2_) {

            // ru2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(e,f) <i,e||a,b>_aaaa t2_aaaa(a,f,m,n) u1_aa(b,i)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempArray[53]("e,a") * t2_aaaa_vvoo("a,f,m,n");
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += -1.000000 <i,e||a,b>_aaaa t2_abab(a,f,m,n) u1_aa(b,i)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[53]("e,a") * t2_abab_vvoo("a,f,m,n");
        }
        world_.gop.fence(); tempArray[53].~TArrayD();

        if (include_u1_) {

            // tempArray[54] += 1.000000 V_blks_["abab_vovv"]("e,i,a,b") u1_bb_vo("b,i")
            // flops: o1v3: 1, o0v2: 1 | mem: o0v2: 2,
            tempArray[54]("e,a") = V_blks_["abab_vovv"]("e,i,a,b") * u1_bb_vo("b,i");

            // tempPerm_aaaa_vvoo += 1.000000 P(e,f) <e,i||a,b>_abab t2_aaaa(a,f,m,n) u1_bb(b,i)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += tempArray[54]("e,a") * t2_aaaa_vvoo("a,f,m,n");
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += 1.000000 <e,i||a,b>_abab t2_abab(a,f,m,n) u1_bb(b,i)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[54]("e,a") * t2_abab_vvoo("a,f,m,n");
        }
        world_.gop.fence(); tempArray[54].~TArrayD();

        if (include_u2_) {

            // tempArray[55] += 1.000000 u2_aaaa_vvoo("a,e,m,i") dp_aa_ov("i,a")
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            tempArray[55]("e,m") = u2_aaaa_vvoo("a,e,m,i") * dp_aa_ov("i,a");

            // rt1_aa_vo += 1.000000 d-_aa(i,a) u2_aaaa(a,e,m,i)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rt1_aa_vo("e,m") += tempArray[55]("e,m");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // ru1_aa_vo += 1.000000 d-_aa(i,a) u0 u2_aaaa(a,e,m,i)
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            ru1_aa_vo("e,m") += tempArray[55]("e,m") * u0;
        }

        if (include_u2_) {

            // ru2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) d-_aa(i,a) u1_aa(e,n) u2_aaaa(a,f,m,i)
            // flops: o2v2: 2 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempArray[55]("f,m") * u1_aa_vo("e,n");

            // ru2_abab_vvoo += 1.000000 d-_aa(i,a) u1_bb(f,n) u2_aaaa(a,e,m,i)
            // flops: o2v2: 2 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[55]("e,m") * u1_bb_vo("f,n");
        }
        world_.gop.fence(); tempArray[55].~TArrayD();

        if (include_u2_) {

            // tempArray[56] += 1.000000 dp_aa_ov("i,a") u2_abab_vvoo("a,e,i,m")
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            tempArray[56]("e,m") = dp_aa_ov("i,a") * u2_abab_vvoo("a,e,i,m");

            // rt1_bb_vo += -1.000000 d-_aa(i,a) u2_abab(a,e,i,m)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rt1_bb_vo("e,m") -= tempArray[56]("e,m");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // ru1_bb_vo += -1.000000 d-_aa(i,a) u0 u2_abab(a,e,i,m)
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            ru1_bb_vo("e,m") -= tempArray[56]("e,m") * u0;
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += -1.000000 d-_aa(i,a) u1_aa(e,m) u2_abab(a,f,i,n)
            // flops: o2v2: 2 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[56]("f,n") * u1_aa_vo("e,m");
        }

        if (include_u2_) {

            // ru2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) d-_aa(i,a) u1_bb(e,n) u2_abab(a,f,i,m)
            // flops: o2v2: 2 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempArray[56]("f,m") * u1_bb_vo("e,n");
        }
        world_.gop.fence(); tempArray[56].~TArrayD();

        if (include_u2_) {

            // tempArray[57] += 1.000000 dp_bb_ov("i,a") u2_bbbb_vvoo("a,e,m,i")
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            tempArray[57]("e,m") = dp_bb_ov("i,a") * u2_bbbb_vvoo("a,e,m,i");

            // rt1_bb_vo += 1.000000 d-_bb(i,a) u2_bbbb(a,e,m,i)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rt1_bb_vo("e,m") += tempArray[57]("e,m");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // ru1_bb_vo += 1.000000 d-_bb(i,a) u0 u2_bbbb(a,e,m,i)
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            ru1_bb_vo("e,m") += tempArray[57]("e,m") * u0;
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += 1.000000 d-_bb(i,a) u1_aa(e,m) u2_bbbb(a,f,n,i)
            // flops: o2v2: 2 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[57]("f,n") * u1_aa_vo("e,m");

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) d-_bb(i,a) u1_bb(e,n) u2_bbbb(a,f,m,i)
            // flops: o2v2: 2 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= tempArray[57]("f,m") * u1_bb_vo("e,n");
        }
        world_.gop.fence(); tempArray[57].~TArrayD();

        if (include_u2_) {

            // tempArray[58] += 1.000000 u2_abab_vvoo("e,a,m,i") dp_bb_ov("i,a")
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            tempArray[58]("e,m") = u2_abab_vvoo("e,a,m,i") * dp_bb_ov("i,a");

            // rt1_aa_vo += -1.000000 d-_bb(i,a) u2_abab(e,a,m,i)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rt1_aa_vo("e,m") -= tempArray[58]("e,m");
        }

        if (include_u0_ && include_u1_ && include_u2_) {

            // ru1_aa_vo += -1.000000 d-_bb(i,a) u0 u2_abab(e,a,m,i)
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            ru1_aa_vo("e,m") -= tempArray[58]("e,m") * u0;
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) d-_bb(i,a) u1_aa(e,n) u2_abab(f,a,m,i)
            // flops: o2v2: 2 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += tempArray[58]("f,m") * u1_aa_vo("e,n");

            // ru2_abab_vvoo += -1.000000 d-_bb(i,a) u1_bb(f,n) u2_abab(e,a,m,i)
            // flops: o2v2: 2 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[58]("e,m") * u1_bb_vo("f,n");
        }
        world_.gop.fence(); tempArray[58].~TArrayD();

        {

            // tempArray[59] += 1.000000 dp_aa_ov("i,a") t2_aaaa_vvoo("a,e,m,i")
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            tempArray[59]("e,m") = dp_aa_ov("i,a") * t2_aaaa_vvoo("a,e,m,i");
        }

        if (include_u0_) {

            // rt1_aa_vo += 1.000000 d-_aa(i,a) t2_aaaa(a,e,m,i) u0
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            rt1_aa_vo("e,m") += tempArray[59]("e,m") * u0;
        }

        if (include_u2_) {

            // ru2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("f,e,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) d-_aa(i,a) t2_aaaa(a,e,n,i) u1_aa(f,m)
            // flops: o2v2: 2 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * tempArray[59]("e,n") * u1_aa_vo("f,m");

            // rt2_abab_vvoo += 1.000000 d-_aa(i,a) t2_aaaa(a,e,m,i) u1_bb(f,n)
            // flops: o2v2: 2 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += tempArray[59]("e,m") * u1_bb_vo("f,n");
        }
        world_.gop.fence(); tempArray[59].~TArrayD();

        {

            // tempArray[60] += 1.000000 t2_abab_vvoo("a,e,i,m") dp_aa_ov("i,a")
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            tempArray[60]("e,m") = t2_abab_vvoo("a,e,i,m") * dp_aa_ov("i,a");
        }

        if (include_u0_) {

            // rt1_bb_vo += -1.000000 d-_aa(i,a) t2_abab(a,e,i,m) u0
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            rt1_bb_vo("e,m") -= tempArray[60]("e,m") * u0;
        }

        if (include_u1_) {

            // rt2_abab_vvoo += -1.000000 d-_aa(i,a) t2_abab(a,f,i,n) u1_aa(e,m)
            // flops: o2v2: 2 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= tempArray[60]("f,n") * u1_aa_vo("e,m");
        }

        if (include_u2_) {

            // ru2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("f,e,n,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) d-_aa(i,a) t2_abab(a,e,i,n) u1_bb(f,m)
            // flops: o2v2: 2 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempArray[60]("e,n") * u1_bb_vo("f,m");
        }
        world_.gop.fence(); tempArray[60].~TArrayD();

        {

            // tempArray[61] += 1.000000 t2_bbbb_vvoo("a,e,m,i") dp_bb_ov("i,a")
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            tempArray[61]("e,m") = t2_bbbb_vvoo("a,e,m,i") * dp_bb_ov("i,a");
        }

        if (include_u0_) {

            // rt1_bb_vo += 1.000000 d-_bb(i,a) t2_bbbb(a,e,m,i) u0
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            rt1_bb_vo("e,m") += tempArray[61]("e,m") * u0;
        }

        if (include_u1_) {

            // rt2_abab_vvoo += 1.000000 d-_bb(i,a) t2_bbbb(a,f,n,i) u1_aa(e,m)
            // flops: o2v2: 2 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += tempArray[61]("f,n") * u1_aa_vo("e,m");

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) d-_bb(i,a) t2_bbbb(a,e,n,i) u1_bb(f,m)
            // flops: o2v2: 2 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= tempArray[61]("e,n") * u1_bb_vo("f,m");
        }
        world_.gop.fence(); tempArray[61].~TArrayD();

        {

            // tempArray[62] += 1.000000 t2_abab_vvoo("e,a,m,i") dp_bb_ov("i,a")
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            tempArray[62]("e,m") = t2_abab_vvoo("e,a,m,i") * dp_bb_ov("i,a");
        }

        if (include_u0_) {

            // rt1_aa_vo += -1.000000 d-_bb(i,a) t2_abab(e,a,m,i) u0
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            rt1_aa_vo("e,m") -= tempArray[62]("e,m") * u0;
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) d-_bb(i,a) t2_abab(e,a,n,i) u1_aa(f,m)
            // flops: o2v2: 2 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += tempArray[62]("e,n") * u1_aa_vo("f,m");

            // rt2_abab_vvoo += -1.000000 d-_bb(i,a) t2_abab(e,a,m,i) u1_bb(f,n)
            // flops: o2v2: 2 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= tempArray[62]("e,m") * u1_bb_vo("f,n");

            // tempArray[63] += 1.000000 V_blks_["aaaa_oovv"]("j,i,a,b") u1_aa_vo("b,j")
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            tempArray[63]("a,i") = V_blks_["aaaa_oovv"]("j,i,a,b") * u1_aa_vo("b,j");

            // ru1_aa_vo += 1.000000 <j,i||a,b>_aaaa t2_aaaa(a,e,m,i) u1_aa(b,j)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_aa_vo("e,m") += tempArray[63]("a,i") * t2_aaaa_vvoo("a,e,m,i");

            // ru1_bb_vo += -1.000000 <j,i||a,b>_aaaa t2_abab(a,e,i,m) u1_aa(b,j)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_bb_vo("e,m") -= tempArray[63]("a,i") * t2_abab_vvoo("a,e,i,m");

            // tempArray[64] += 1.000000 u1_bb_vo("b,j") V_blks_["abab_oovv"]("i,j,a,b")
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            tempArray[64]("a,i") = u1_bb_vo("b,j") * V_blks_["abab_oovv"]("i,j,a,b");

            // ru1_aa_vo += -1.000000 <i,j||a,b>_abab t2_aaaa(a,e,m,i) u1_bb(b,j)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_aa_vo("e,m") -= tempArray[64]("a,i") * t2_aaaa_vvoo("a,e,m,i");

            // ru1_bb_vo += 1.000000 <i,j||a,b>_abab t2_abab(a,e,i,m) u1_bb(b,j)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_bb_vo("e,m") += tempArray[64]("a,i") * t2_abab_vvoo("a,e,i,m");

            // tempArray[65] += 1.000000 V_blks_["abab_oovv"]("j,i,b,a") u1_aa_vo("b,j")
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            tempArray[65]("a,i") = V_blks_["abab_oovv"]("j,i,b,a") * u1_aa_vo("b,j");

            // ru1_aa_vo += 1.000000 <j,i||b,a>_abab t2_abab(e,a,m,i) u1_aa(b,j)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_aa_vo("e,m") += tempArray[65]("a,i") * t2_abab_vvoo("e,a,m,i");

            // ru1_bb_vo += -1.000000 <j,i||b,a>_abab t2_bbbb(a,e,m,i) u1_aa(b,j)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_bb_vo("e,m") -= tempArray[65]("a,i") * t2_bbbb_vvoo("a,e,m,i");

            // tempArray[66] += 1.000000 V_blks_["bbbb_oovv"]("j,i,a,b") u1_bb_vo("b,j")
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            tempArray[66]("a,i") = V_blks_["bbbb_oovv"]("j,i,a,b") * u1_bb_vo("b,j");

            // ru1_aa_vo += -1.000000 <j,i||a,b>_bbbb t2_abab(e,a,m,i) u1_bb(b,j)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_aa_vo("e,m") -= tempArray[66]("a,i") * t2_abab_vvoo("e,a,m,i");

            // ru1_bb_vo += 1.000000 <j,i||a,b>_bbbb t2_bbbb(a,e,m,i) u1_bb(b,j)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_bb_vo("e,m") += tempArray[66]("a,i") * t2_bbbb_vvoo("a,e,m,i");

            // tempArray[67] += 1.000000 dp_bb_vv("e,a") u1_bb_vo("a,m")
            // flops: o1v2: 1, o1v1: 1 | mem: o1v1: 2,
            tempArray[67]("e,m") = dp_bb_vv("e,a") * u1_bb_vo("a,m");

            // rt1_bb_vo += -1.000000 d-_bb(e,a) u1_bb(a,m)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rt1_bb_vo("e,m") -= tempArray[67]("e,m");
        }
        world_.gop.fence(); tempArray[62].~TArrayD();
        world_.gop.fence(); tempArray[63].~TArrayD();
        world_.gop.fence(); tempArray[64].~TArrayD();
        world_.gop.fence(); tempArray[65].~TArrayD();
        world_.gop.fence(); tempArray[66].~TArrayD();

        if (include_u0_ && include_u1_) {

            // ru1_bb_vo += -1.000000 d-_bb(e,a) u0 u1_bb(a,m)
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            ru1_bb_vo("e,m") -= tempArray[67]("e,m") * u0;
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += -1.000000 d-_bb(f,a) u1_bb(a,n) u1_aa(e,m)
            // flops: o2v2: 2 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[67]("f,n") * u1_aa_vo("e,m");
        }

        {

            // rt2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("f,e,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) d-_bb(e,a) u1_bb(a,n) u1_bb(f,m)
            // flops: o2v2: 2 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempArray[67]("e,n") * u1_bb_vo("f,m");

            // tempArray[68] += 1.000000 dp_aa_vv("e,a") u1_aa_vo("a,m")
            // flops: o1v2: 1, o1v1: 1 | mem: o1v1: 2,
            tempArray[68]("e,m") = dp_aa_vv("e,a") * u1_aa_vo("a,m");

            // rt1_aa_vo += -1.000000 d-_aa(e,a) u1_aa(a,m)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rt1_aa_vo("e,m") -= tempArray[68]("e,m");
        }
        world_.gop.fence(); tempArray[67].~TArrayD();

        if (include_u0_ && include_u1_) {

            // ru1_aa_vo += -1.000000 d-_aa(e,a) u0 u1_aa(a,m)
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            ru1_aa_vo("e,m") -= tempArray[68]("e,m") * u0;
        }

        {

            // rt2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("f,e,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) d-_aa(e,a) u1_aa(a,n) u1_aa(f,m)
            // flops: o2v2: 2 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempArray[68]("e,n") * u1_aa_vo("f,m");
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += -1.000000 d-_aa(e,a) u1_aa(a,m) u1_bb(f,n)
            // flops: o2v2: 2 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= tempArray[68]("e,m") * u1_bb_vo("f,n");
        }
        world_.gop.fence(); tempArray[68].~TArrayD();

        if (include_u1_) {

            // tempArray[69] += 1.000000 dp_aa_oo("i,m") u1_aa_vo("e,i")
            // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
            tempArray[69]("e,m") = dp_aa_oo("i,m") * u1_aa_vo("e,i");

            // rt1_aa_vo += 1.000000 d-_aa(i,m) u1_aa(e,i)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rt1_aa_vo("e,m") += tempArray[69]("e,m");
        }

        if (include_u0_ && include_u1_) {

            // ru1_aa_vo += 1.000000 d-_aa(i,m) u0 u1_aa(e,i)
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            ru1_aa_vo("e,m") += tempArray[69]("e,m") * u0;
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) d-_aa(i,n) u1_aa(e,i) u1_aa(f,m)
            // flops: o2v2: 2 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= tempArray[69]("e,n") * u1_aa_vo("f,m");
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += 1.000000 d-_aa(i,m) u1_aa(e,i) u1_bb(f,n)
            // flops: o2v2: 2 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[69]("e,m") * u1_bb_vo("f,n");
        }
        world_.gop.fence(); tempArray[69].~TArrayD();

        if (include_u1_) {

            // tempArray[70] += 1.000000 u1_bb_vo("e,i") dp_bb_oo("i,m")
            // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
            tempArray[70]("e,m") = u1_bb_vo("e,i") * dp_bb_oo("i,m");

            // rt1_bb_vo += 1.000000 d-_bb(i,m) u1_bb(e,i)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rt1_bb_vo("e,m") += tempArray[70]("e,m");
        }

        if (include_u0_ && include_u1_) {

            // ru1_bb_vo += 1.000000 d-_bb(i,m) u0 u1_bb(e,i)
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            ru1_bb_vo("e,m") += tempArray[70]("e,m") * u0;
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += 1.000000 d-_bb(i,n) u1_bb(f,i) u1_aa(e,m)
            // flops: o2v2: 2 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += tempArray[70]("f,n") * u1_aa_vo("e,m");
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) d-_bb(i,n) u1_bb(e,i) u1_bb(f,m)
            // flops: o2v2: 2 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= tempArray[70]("e,n") * u1_bb_vo("f,m");
        }
        world_.gop.fence(); tempArray[70].~TArrayD();

        {

            // energy+= 1.000000 f_aa(i,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy+= scalar0;

            // energy+= 1.000000 f_bb(i,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy+= scalar7;
        }

        if (include_u0_) {

            // energy+= -1.000000 d-_aa(i,i) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy-= scalar3 * u0;

            // energy+= -1.000000 d-_bb(i,i) u0
            // flops: o0v0: 2 | mem: o0v0: 2,
            energy-= scalar10 * u0;
        }

        if (include_u1_) {

            // energy+= -1.000000 d-_aa(i,a) u1_aa(a,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy-= scalar5;

            // energy+= -1.000000 d-_bb(i,a) u1_bb(a,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy-= scalar9;
        }

        {

            // energy+= -0.500000 <j,i||j,i>_aaaa
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy-= 0.500000 * scalar2;

            // energy+= -0.500000 <j,i||j,i>_abab
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy-= 0.500000 * scalar8;

            // energy+= -0.500000 <i,j||i,j>_abab
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy-= 0.500000 * scalar8;

            // energy+= -0.500000 <j,i||j,i>_bbbb
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy-= 0.500000 * scalar4;

            // energy+= 0.250000 <j,i||a,b>_aaaa t2_aaaa(a,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy+= 0.250000 * scalar1;

            // energy+= 0.250000 <j,i||a,b>_abab t2_abab(a,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy+= 0.250000 * scalar6;

            // energy+= 0.250000 <i,j||a,b>_abab t2_abab(a,b,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy+= 0.250000 * scalar6;

            // energy+= 0.250000 <j,i||b,a>_abab t2_abab(b,a,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy+= 0.250000 * scalar6;

            // energy+= 0.250000 <i,j||b,a>_abab t2_abab(b,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy+= 0.250000 * scalar6;

            // energy+= 0.250000 <j,i||a,b>_bbbb t2_bbbb(a,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy+= 0.250000 * scalar11;

            // rt1_aa_vo += 1.000000 f_bb(i,a) t2_abab(e,a,m,i)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_aa_vo("e,m") += F_blks_["bb_ov"]("i,a") * t2_abab_vvoo("e,a,m,i");

            // rt1_aa_vo += -0.500000 <i,e||a,b>_aaaa t2_aaaa(a,b,m,i)
            // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_aa_vo("e,m") += 0.500000 * V_blks_["aaaa_vovv"]("e,i,a,b") * t2_aaaa_vvoo("a,b,m,i");

            // rt1_aa_vo += 0.500000 <e,i||a,b>_abab t2_abab(a,b,m,i)
            // +            0.500000 <e,i||b,a>_abab t2_abab(b,a,m,i)
            // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_aa_vo("e,m") += V_blks_["abab_vovv"]("e,i,a,b") * t2_abab_vvoo("a,b,m,i");
        }

        if (include_u0_) {

            // rt1_aa_vo += -1.000000 d-_aa(e,m) u0
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            rt1_aa_vo("e,m") -= dp_aa_vo("e,m") * u0;
        }

        if (include_u1_) {

            // rt1_aa_vo += -1.000000 d-_aa(i,i) u1_aa(e,m)
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            rt1_aa_vo("e,m") -= scalar3 * u1_aa_vo("e,m");
        }

        {

            // rt1_aa_vo += -0.500000 <j,i||a,m>_aaaa t2_aaaa(a,e,j,i)
            // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_aa_vo("e,m") -= 0.500000 * V_blks_["aaaa_oovo"]("j,i,a,m") * t2_aaaa_vvoo("a,e,j,i");

            // rt1_aa_vo += -1.000000 f_aa(i,a) t2_aaaa(a,e,m,i)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_aa_vo("e,m") -= F_blks_["aa_ov"]("i,a") * t2_aaaa_vvoo("a,e,m,i");

            // rt1_aa_vo += 1.000000 f_aa(e,m)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rt1_aa_vo("e,m") += F_blks_["aa_vo"]("e,m");
        }

        if (include_u1_) {

            // rt1_aa_vo += -1.000000 d-_bb(i,i) u1_aa(e,m)
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            rt1_aa_vo("e,m") -= scalar10 * u1_aa_vo("e,m");
        }

        {

            // rt1_aa_vo += -0.500000 <j,i||m,a>_abab t2_abab(e,a,j,i)
            // +            -0.500000 <i,j||m,a>_abab t2_abab(e,a,i,j)
            // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_aa_vo("e,m") += V_blks_["abba_oovo"]("j,i,a,m") * t2_abab_vvoo("e,a,j,i");
        }

        if (include_u1_) {

            // rt1_bb_vo += -1.000000 d-_bb(i,i) u1_bb(e,m)
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            rt1_bb_vo("e,m") -= scalar10 * u1_bb_vo("e,m");
        }

        {

            // rt1_bb_vo += -0.500000 <i,e||a,b>_bbbb t2_bbbb(a,b,m,i)
            // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_bb_vo("e,m") += 0.500000 * V_blks_["bbbb_vovv"]("e,i,a,b") * t2_bbbb_vvoo("a,b,m,i");
        }

        if (include_u1_) {

            // rt1_bb_vo += -1.000000 d-_aa(i,i) u1_bb(e,m)
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            rt1_bb_vo("e,m") -= scalar3 * u1_bb_vo("e,m");
        }

        {

            // rt1_bb_vo += 1.000000 f_aa(i,a) t2_abab(a,e,i,m)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_bb_vo("e,m") += F_blks_["aa_ov"]("i,a") * t2_abab_vvoo("a,e,i,m");

            // rt1_bb_vo += 1.000000 f_bb(e,m)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rt1_bb_vo("e,m") += F_blks_["bb_vo"]("e,m");
        }

        if (include_u0_) {

            // rt1_bb_vo += -1.000000 d-_bb(e,m) u0
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            rt1_bb_vo("e,m") -= dp_bb_vo("e,m") * u0;
        }

        {

            // rt1_bb_vo += -0.500000 <j,i||a,m>_abab t2_abab(a,e,j,i)
            // +            -0.500000 <i,j||a,m>_abab t2_abab(a,e,i,j)
            // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_bb_vo("e,m") -= V_blks_["abab_oovo"]("j,i,a,m") * t2_abab_vvoo("a,e,j,i");

            // rt1_bb_vo += -0.500000 <j,i||a,m>_bbbb t2_bbbb(a,e,j,i)
            // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_bb_vo("e,m") -= 0.500000 * V_blks_["bbbb_oovo"]("j,i,a,m") * t2_bbbb_vvoo("a,e,j,i");

            // rt1_bb_vo += -1.000000 f_bb(i,a) t2_bbbb(a,e,m,i)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_bb_vo("e,m") -= F_blks_["bb_ov"]("i,a") * t2_bbbb_vvoo("a,e,m,i");

            // rt1_bb_vo += 0.500000 <i,e||a,b>_abab t2_abab(a,b,i,m)
            // +            0.500000 <i,e||b,a>_abab t2_abab(b,a,i,m)
            // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_bb_vo("e,m") -= V_blks_["baab_vovv"]("e,i,a,b") * t2_abab_vvoo("a,b,i,m");

            // rt2_aaaa_vvoo += -0.500000 <j,i||a,b>_aaaa t2_aaaa(a,e,m,n) t2_aaaa(b,f,j,i)
            // flops: o2v3: 2, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") -= 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * t2_aaaa_vvoo("b,f,j,i") * t2_aaaa_vvoo("a,e,m,n");

            // rt2_aaaa_vvoo += 1.000000 <e,f||m,n>_aaaa
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += V_blks_["aaaa_vvoo"]("e,f,m,n");

            // rt2_aaaa_vvoo += 0.500000 <j,i||m,n>_aaaa t2_aaaa(e,f,j,i)
            // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_aaaa_vvoo("e,f,m,n") += 0.500000 * V_blks_["aaaa_oooo"]("j,i,m,n") * t2_aaaa_vvoo("e,f,j,i");
        }

        if (include_u2_) {

            // rt2_aaaa_vvoo += -1.000000 d-_aa(i,i) u2_aaaa(e,f,m,n)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            rt2_aaaa_vvoo("e,f,m,n") -= scalar3 * u2_aaaa_vvoo("e,f,m,n");
        }

        {

            // rt2_aaaa_vvoo += 0.500000 <e,f||a,b>_aaaa t2_aaaa(a,b,m,n)
            // flops: o2v4: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_aaaa_vvoo("e,f,m,n") += 0.500000 * V_blks_["aaaa_vvvv"]("e,f,a,b") * t2_aaaa_vvoo("a,b,m,n");
        }

        if (include_u2_) {

            // rt2_aaaa_vvoo += -1.000000 d-_bb(i,i) u2_aaaa(e,f,m,n)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            rt2_aaaa_vvoo("e,f,m,n") -= scalar10 * u2_aaaa_vvoo("e,f,m,n");

            // ru2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("f,e,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();

            // tempPerm_aaaa_vvoo += 1.000000 P(e,f) f_aa(e,a) t2_aaaa(a,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = F_blks_["aa_vv"]("e,a") * t2_aaaa_vvoo("a,f,m,n");
        }

        if (include_u0_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(e,f) d-_aa(e,a) t2_aaaa(a,f,m,n) u0
            // flops: o2v3: 1, o2v2: 1, o0v0: 1 | mem: o2v2: 2, o0v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") -= dp_aa_vv("e,a") * u0 * t2_aaaa_vvoo("a,f,m,n");
        }

        {

            // rt2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u0_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) d-_aa(i,n) t2_aaaa(e,f,m,i) u0
            // flops: o3v2: 1, o2v2: 1, o0v0: 1 | mem: o2v2: 2, o0v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = dp_aa_oo("i,n") * u0 * t2_aaaa_vvoo("e,f,m,i");
        }

        {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) f_aa(i,n) t2_aaaa(e,f,m,i)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= F_blks_["aa_oo"]("i,n") * t2_aaaa_vvoo("e,f,m,i");

            // rt2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <i,e||a,n>_aaaa t2_aaaa(a,f,m,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * V_blks_["aaaa_vovo"]("e,i,a,n") * t2_aaaa_vvoo("a,f,m,i");

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) <e,i||n,a>_abab t2_abab(f,a,m,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += V_blks_["abba_vovo"]("e,i,a,n") * t2_abab_vvoo("f,a,m,i");
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) d-_aa(e,n) u1_aa(f,m)
            // flops: o2v2: 2 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += dp_aa_vo("e,n") * u1_aa_vo("f,m");

            // rt2_abab_vvoo += -1.000000 d-_bb(f,n) u1_aa(e,m)
            // flops: o2v2: 2 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= dp_bb_vo("f,n") * u1_aa_vo("e,m");
        }

        {

            // rt2_abab_vvoo += -1.000000 <e,i||m,a>_abab t2_bbbb(a,f,n,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += V_blks_["abba_vovo"]("e,i,a,m") * t2_bbbb_vvoo("a,f,n,i");
        }

        if (include_u0_) {

            // rt2_abab_vvoo += -1.000000 d-_bb(f,a) t2_abab(e,a,m,n) u0
            // flops: o2v3: 1, o2v2: 1, o0v0: 1 | mem: o2v2: 2, o0v0: 1,
            rt2_abab_vvoo("e,f,m,n") -= dp_bb_vv("f,a") * t2_abab_vvoo("e,a,m,n") * u0;

            // rt2_abab_vvoo += -1.000000 d-_aa(e,a) t2_abab(a,f,m,n) u0
            // flops: o2v3: 1, o2v2: 1, o0v0: 1 | mem: o2v2: 2, o0v0: 1,
            rt2_abab_vvoo("e,f,m,n") -= dp_aa_vv("e,a") * t2_abab_vvoo("a,f,m,n") * u0;
        }

        {

            // rt2_abab_vvoo += -1.000000 <e,i||a,n>_abab t2_abab(a,f,m,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= V_blks_["abab_vovo"]("e,i,a,n") * t2_abab_vvoo("a,f,m,i");
        }

        if (include_u2_) {

            // rt2_abab_vvoo += -1.000000 d-_aa(i,i) u2_abab(e,f,m,n)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            rt2_abab_vvoo("e,f,m,n") -= scalar3 * u2_abab_vvoo("e,f,m,n");
        }

        {

            // rt2_abab_vvoo += -1.000000 f_bb(i,n) t2_abab(e,f,m,i)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= F_blks_["bb_oo"]("i,n") * t2_abab_vvoo("e,f,m,i");

            // rt2_abab_vvoo += 1.000000 <i,f||a,n>_bbbb t2_abab(e,a,m,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= V_blks_["bbbb_vovo"]("f,i,a,n") * t2_abab_vvoo("e,a,m,i");
        }

        if (include_u0_) {

            // rt2_abab_vvoo += 1.000000 d-_bb(i,n) t2_abab(e,f,m,i) u0
            // flops: o3v2: 1, o2v2: 1, o0v0: 1 | mem: o2v2: 2, o0v0: 1,
            rt2_abab_vvoo("e,f,m,n") += dp_bb_oo("i,n") * t2_abab_vvoo("e,f,m,i") * u0;
        }

        {

            // rt2_abab_vvoo += 1.000000 <e,f||m,n>_abab
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_abab_vvoo("e,f,m,n") += V_blks_["abab_vvoo"]("e,f,m,n");

            // rt2_abab_vvoo += 1.000000 f_aa(e,a) t2_abab(a,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += F_blks_["aa_vv"]("e,a") * t2_abab_vvoo("a,f,m,n");
        }

        if (include_u1_) {

            // rt2_abab_vvoo += -1.000000 d-_aa(e,m) u1_bb(f,n)
            // flops: o2v2: 2 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= dp_aa_vo("e,m") * u1_bb_vo("f,n");
        }

        if (include_u2_) {

            // rt2_abab_vvoo += -1.000000 d-_bb(i,i) u2_abab(e,f,m,n)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            rt2_abab_vvoo("e,f,m,n") -= scalar10 * u2_abab_vvoo("e,f,m,n");
        }

        {

            // rt2_abab_vvoo += 1.000000 f_bb(f,a) t2_abab(e,a,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += F_blks_["bb_vv"]("f,a") * t2_abab_vvoo("e,a,m,n");

            // rt2_abab_vvoo += 0.500000 <e,f||a,b>_abab t2_abab(a,b,m,n)
            // +                0.500000 <e,f||b,a>_abab t2_abab(b,a,m,n)
            // flops: o2v4: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += V_blks_["abab_vvvv"]("e,f,a,b") * t2_abab_vvoo("a,b,m,n");
        }

        if (include_u0_) {

            // rt2_abab_vvoo += 1.000000 d-_aa(i,m) t2_abab(e,f,i,n) u0
            // flops: o3v2: 1, o2v2: 1, o0v0: 1 | mem: o2v2: 2, o0v0: 1,
            rt2_abab_vvoo("e,f,m,n") += dp_aa_oo("i,m") * t2_abab_vvoo("e,f,i,n") * u0;
        }

        {

            // rt2_abab_vvoo += -1.000000 f_aa(i,m) t2_abab(e,f,i,n)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= F_blks_["aa_oo"]("i,m") * t2_abab_vvoo("e,f,i,n");

            // rt2_abab_vvoo += -1.000000 <i,f||a,n>_abab t2_aaaa(a,e,m,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += V_blks_["baab_vovo"]("f,i,a,n") * t2_aaaa_vvoo("a,e,m,i");

            // rt2_abab_vvoo += 0.500000 <j,i||m,n>_abab t2_abab(e,f,j,i)
            // +                0.500000 <i,j||m,n>_abab t2_abab(e,f,i,j)
            // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += V_blks_["abab_oooo"]("j,i,m,n") * t2_abab_vvoo("e,f,j,i");

            // rt2_abab_vvoo += 1.000000 <i,e||a,m>_aaaa t2_abab(a,f,i,n)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= V_blks_["aaaa_vovo"]("e,i,a,m") * t2_abab_vvoo("a,f,i,n");

            // rt2_abab_vvoo += -1.000000 <i,f||m,a>_abab t2_abab(e,a,i,n)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= V_blks_["baba_vovo"]("f,i,a,m") * t2_abab_vvoo("e,a,i,n");

            // rt2_bbbb_vvoo += 0.500000 <j,i||m,n>_bbbb t2_bbbb(e,f,j,i)
            // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_bbbb_vvoo("e,f,m,n") += 0.500000 * V_blks_["bbbb_oooo"]("j,i,m,n") * t2_bbbb_vvoo("e,f,j,i");
        }

        if (include_u2_) {

            // rt2_bbbb_vvoo += -1.000000 d-_bb(i,i) u2_bbbb(e,f,m,n)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            rt2_bbbb_vvoo("e,f,m,n") -= scalar10 * u2_bbbb_vvoo("e,f,m,n");

            // rt2_bbbb_vvoo += -1.000000 d-_aa(i,i) u2_bbbb(e,f,m,n)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            rt2_bbbb_vvoo("e,f,m,n") -= scalar3 * u2_bbbb_vvoo("e,f,m,n");
        }

        {

            // rt2_bbbb_vvoo += 1.000000 <e,f||m,n>_bbbb
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += V_blks_["bbbb_vvoo"]("e,f,m,n");

            // rt2_bbbb_vvoo += 0.500000 <e,f||a,b>_bbbb t2_bbbb(a,b,m,n)
            // flops: o2v4: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_bbbb_vvoo("e,f,m,n") += 0.500000 * V_blks_["bbbb_vvvv"]("e,f,a,b") * t2_bbbb_vvoo("a,b,m,n");
        }

        if (include_u2_) {

            // ru2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("f,e,n,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u0_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(e,f) d-_bb(e,a) t2_bbbb(a,f,m,n) u0
            // flops: o2v3: 1, o2v2: 1, o0v0: 1 | mem: o2v2: 2, o0v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * dp_bb_vv("e,a") * u0 * t2_bbbb_vvoo("a,f,m,n");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 P(e,f) f_bb(e,a) t2_bbbb(a,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += F_blks_["bb_vv"]("e,a") * t2_bbbb_vvoo("a,f,m,n");

            // rt2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) f_bb(i,n) t2_bbbb(e,f,m,i)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * F_blks_["bb_oo"]("i,n") * t2_bbbb_vvoo("e,f,m,i");
        }

        if (include_u0_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) d-_bb(i,n) t2_bbbb(e,f,m,i) u0
            // flops: o3v2: 1, o2v2: 1, o0v0: 1 | mem: o2v2: 2, o0v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") += dp_bb_oo("i,n") * u0 * t2_bbbb_vvoo("e,f,m,i");
        }

        {

            // rt2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) d-_bb(e,n) u1_bb(f,m)
            // flops: o2v2: 2 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = dp_bb_vo("e,n") * u1_bb_vo("f,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <i,e||a,n>_bbbb t2_bbbb(a,f,m,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= V_blks_["bbbb_vovo"]("e,i,a,n") * t2_bbbb_vvoo("a,f,m,i");

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) <i,e||a,n>_abab t2_abab(a,f,i,m)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += V_blks_["baab_vovo"]("e,i,a,n") * t2_abab_vvoo("a,f,i,m");
        }

        if (include_u0_) {

            // ru0 += 1.000000 u0 w0
            // flops: o0v0: 2 | mem: o0v0: 2,
            ru0 += u0 * w0;
        }

        if (include_u0_ && include_u1_) {

            // ru0 += 1.000000 f_aa(i,a) u1_aa(a,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            ru0 += scalar14;

            // ru0 += 1.000000 f_bb(i,a) u1_bb(a,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            ru0 += scalar12;
        }

        if (include_u0_) {

            // ru0 += -1.000000 d+_aa(i,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            ru0 -= scalar18;

            // ru0 += -1.000000 d+_bb(i,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            ru0 -= scalar15;
        }

        if (include_u0_ && include_u1_) {

            // ru0 += -1.000000 d-_aa(i,a) u0 u1_aa(a,i)
            // flops: o0v0: 2 | mem: o0v0: 2,
            ru0 -= scalar5 * u0;

            // ru0 += -1.000000 d-_bb(i,a) u0 u1_bb(a,i)
            // flops: o0v0: 2 | mem: o0v0: 2,
            ru0 -= scalar9 * u0;
        }

        if (include_u0_ && include_u2_) {

            // ru0 += 0.250000 <j,i||a,b>_aaaa u2_aaaa(a,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            ru0 += 0.250000 * scalar13;

            // ru0 += 0.250000 <j,i||a,b>_abab u2_abab(a,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            ru0 += 0.250000 * scalar16;

            // ru0 += 0.250000 <i,j||a,b>_abab u2_abab(a,b,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            ru0 += 0.250000 * scalar16;

            // ru0 += 0.250000 <j,i||b,a>_abab u2_abab(b,a,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            ru0 += 0.250000 * scalar16;

            // ru0 += 0.250000 <i,j||b,a>_abab u2_abab(b,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            ru0 += 0.250000 * scalar16;

            // ru0 += 0.250000 <j,i||a,b>_bbbb u2_bbbb(a,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            ru0 += 0.250000 * scalar17;
        }

        if (include_u1_ && include_u2_) {

            // ru1_aa_vo += -0.500000 <j,i||m,a>_abab u2_abab(e,a,j,i)
            // +            -0.500000 <i,j||m,a>_abab u2_abab(e,a,i,j)
            // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_aa_vo("e,m") += V_blks_["abba_oovo"]("j,i,a,m") * u2_abab_vvoo("e,a,j,i");
        }

        if (include_u1_) {

            // ru1_aa_vo += -1.000000 d-_aa(i,a) u1_aa(a,i) u1_aa(e,m)
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            ru1_aa_vo("e,m") -= scalar5 * u1_aa_vo("e,m");
        }

        if (include_u1_ && include_u2_) {

            // ru1_aa_vo += -0.500000 <i,e||a,b>_aaaa u2_aaaa(a,b,m,i)
            // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_aa_vo("e,m") += 0.500000 * V_blks_["aaaa_vovv"]("e,i,a,b") * u2_aaaa_vvoo("a,b,m,i");
        }

        if (include_u1_) {

            // ru1_aa_vo += -1.000000 d-_bb(i,a) u1_bb(a,i) u1_aa(e,m)
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            ru1_aa_vo("e,m") -= scalar9 * u1_aa_vo("e,m");

            // ru1_aa_vo += 1.000000 f_aa(e,a) u1_aa(a,m)
            // flops: o1v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_aa_vo("e,m") += F_blks_["aa_vv"]("e,a") * u1_aa_vo("a,m");

            // ru1_aa_vo += 1.000000 u1_aa(e,m) w0
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            ru1_aa_vo("e,m") += u1_aa_vo("e,m") * w0;
        }

        if (include_u1_ && include_u2_) {

            // ru1_aa_vo += -1.000000 f_aa(i,a) u2_aaaa(a,e,m,i)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_aa_vo("e,m") -= F_blks_["aa_ov"]("i,a") * u2_aaaa_vvoo("a,e,m,i");
        }

        if (include_u1_) {

            // ru1_aa_vo += 2.000000 d-_aa(i,a) u1_aa(a,m) u1_aa(e,i)
            // flops: o2v1: 2, o1v1: 1 | mem: o1v1: 2, o2v0: 1,
            ru1_aa_vo("e,m") += 2.000000 * dp_aa_ov("i,a") * u1_aa_vo("a,m") * u1_aa_vo("e,i");

            // ru1_aa_vo += -1.000000 d+_aa(e,m)
            // flops: o1v1: 1 | mem: o1v1: 1,
            ru1_aa_vo("e,m") -= dp_aa_vo("e,m");

            // ru1_aa_vo += -1.000000 f_aa(i,m) u1_aa(e,i)
            // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_aa_vo("e,m") -= F_blks_["aa_oo"]("i,m") * u1_aa_vo("e,i");

            // ru1_aa_vo += -0.500000 <j,i||b,a>_abab t2_abab(e,a,j,i) u1_aa(b,m)
            // +            -0.500000 <i,j||b,a>_abab t2_abab(e,a,i,j) u1_aa(b,m)
            // flops: o3v2: 2, o1v1: 1 | mem: o3v1: 1, o1v1: 2,
            ru1_aa_vo("e,m") -= V_blks_["abab_oovv"]("j,i,b,a") * u1_aa_vo("b,m") * t2_abab_vvoo("e,a,j,i");

            // ru1_aa_vo += -0.500000 <j,i||a,b>_aaaa t2_aaaa(a,e,j,i) u1_aa(b,m)
            // flops: o3v2: 2, o1v1: 1 | mem: o3v1: 1, o1v1: 2,
            ru1_aa_vo("e,m") -= 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * u1_aa_vo("b,m") * t2_aaaa_vvoo("a,e,j,i");

            // ru1_aa_vo += 1.000000 <i,e||a,m>_aaaa u1_aa(a,i)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_aa_vo("e,m") -= V_blks_["aaaa_vovo"]("e,i,a,m") * u1_aa_vo("a,i");
        }

        if (include_u1_ && include_u2_) {

            // ru1_aa_vo += 1.000000 f_bb(i,a) u2_abab(e,a,m,i)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_aa_vo("e,m") += F_blks_["bb_ov"]("i,a") * u2_abab_vvoo("e,a,m,i");
        }

        if (include_u1_) {

            // ru1_aa_vo += 1.000000 d+_aa(i,a) t2_aaaa(a,e,m,i)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_aa_vo("e,m") += dp_aa_ov("i,a") * t2_aaaa_vvoo("a,e,m,i");
        }

        if (include_u1_ && include_u2_) {

            // ru1_aa_vo += 0.500000 <e,i||a,b>_abab u2_abab(a,b,m,i)
            // +            0.500000 <e,i||b,a>_abab u2_abab(b,a,m,i)
            // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_aa_vo("e,m") += V_blks_["abab_vovv"]("e,i,a,b") * u2_abab_vvoo("a,b,m,i");
        }

        if (include_u1_) {

            // ru1_aa_vo += 1.000000 <e,i||m,a>_abab u1_bb(a,i)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_aa_vo("e,m") -= V_blks_["abba_vovo"]("e,i,a,m") * u1_bb_vo("a,i");

            // ru1_aa_vo += -1.000000 d+_bb(i,a) t2_abab(e,a,m,i)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_aa_vo("e,m") -= dp_bb_ov("i,a") * t2_abab_vvoo("e,a,m,i");
        }

        if (include_u1_ && include_u2_) {

            // ru1_aa_vo += -0.500000 <j,i||a,m>_aaaa u2_aaaa(a,e,j,i)
            // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_aa_vo("e,m") -= 0.500000 * V_blks_["aaaa_oovo"]("j,i,a,m") * u2_aaaa_vvoo("a,e,j,i");

            // ru1_bb_vo += -0.500000 <i,e||a,b>_bbbb u2_bbbb(a,b,m,i)
            // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_bb_vo("e,m") += 0.500000 * V_blks_["bbbb_vovv"]("e,i,a,b") * u2_bbbb_vvoo("a,b,m,i");
        }

        if (include_u1_) {

            // ru1_bb_vo += -1.000000 d-_aa(i,a) u1_aa(a,i) u1_bb(e,m)
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            ru1_bb_vo("e,m") -= scalar5 * u1_bb_vo("e,m");
        }

        if (include_u1_ && include_u2_) {

            // ru1_bb_vo += -1.000000 f_bb(i,a) u2_bbbb(a,e,m,i)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_bb_vo("e,m") -= F_blks_["bb_ov"]("i,a") * u2_bbbb_vvoo("a,e,m,i");
        }

        if (include_u1_) {

            // ru1_bb_vo += -1.000000 d+_aa(i,a) t2_abab(a,e,i,m)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_bb_vo("e,m") -= dp_aa_ov("i,a") * t2_abab_vvoo("a,e,i,m");
        }

        if (include_u1_ && include_u2_) {

            // ru1_bb_vo += 0.500000 <i,e||a,b>_abab u2_abab(a,b,i,m)
            // +            0.500000 <i,e||b,a>_abab u2_abab(b,a,i,m)
            // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_bb_vo("e,m") -= V_blks_["baab_vovv"]("e,i,a,b") * u2_abab_vvoo("a,b,i,m");

            // ru1_bb_vo += 1.000000 f_aa(i,a) u2_abab(a,e,i,m)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_bb_vo("e,m") += F_blks_["aa_ov"]("i,a") * u2_abab_vvoo("a,e,i,m");
        }

        if (include_u1_) {

            // ru1_bb_vo += 1.000000 d+_bb(i,a) t2_bbbb(a,e,m,i)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_bb_vo("e,m") += dp_bb_ov("i,a") * t2_bbbb_vvoo("a,e,m,i");

            // ru1_bb_vo += 1.000000 u1_bb(e,m) w0
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            ru1_bb_vo("e,m") += u1_bb_vo("e,m") * w0;

            // ru1_bb_vo += -0.500000 <j,i||a,b>_bbbb t2_bbbb(a,e,j,i) u1_bb(b,m)
            // flops: o3v2: 2, o1v1: 1 | mem: o3v1: 1, o1v1: 2,
            ru1_bb_vo("e,m") -= 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * u1_bb_vo("b,m") * t2_bbbb_vvoo("a,e,j,i");
        }

        if (include_u1_ && include_u2_) {

            // ru1_bb_vo += -0.500000 <j,i||a,m>_abab u2_abab(a,e,j,i)
            // +            -0.500000 <i,j||a,m>_abab u2_abab(a,e,i,j)
            // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_bb_vo("e,m") -= V_blks_["abab_oovo"]("j,i,a,m") * u2_abab_vvoo("a,e,j,i");
        }

        if (include_u1_) {

            // ru1_bb_vo += -0.500000 <j,i||a,b>_abab t2_abab(a,e,j,i) u1_bb(b,m)
            // +            -0.500000 <i,j||a,b>_abab t2_abab(a,e,i,j) u1_bb(b,m)
            // flops: o3v2: 2, o1v1: 1 | mem: o3v1: 1, o1v1: 2,
            ru1_bb_vo("e,m") -= V_blks_["abab_oovv"]("j,i,a,b") * u1_bb_vo("b,m") * t2_abab_vvoo("a,e,j,i");

            // ru1_bb_vo += -1.000000 f_bb(i,m) u1_bb(e,i)
            // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_bb_vo("e,m") -= F_blks_["bb_oo"]("i,m") * u1_bb_vo("e,i");

            // ru1_bb_vo += 2.000000 d-_bb(i,a) u1_bb(a,m) u1_bb(e,i)
            // flops: o2v1: 2, o1v1: 1 | mem: o1v1: 2, o2v0: 1,
            ru1_bb_vo("e,m") += 2.000000 * dp_bb_ov("i,a") * u1_bb_vo("a,m") * u1_bb_vo("e,i");

            // ru1_bb_vo += -1.000000 d+_bb(e,m)
            // flops: o1v1: 1 | mem: o1v1: 1,
            ru1_bb_vo("e,m") -= dp_bb_vo("e,m");
        }

        if (include_u1_ && include_u2_) {

            // ru1_bb_vo += -0.500000 <j,i||a,m>_bbbb u2_bbbb(a,e,j,i)
            // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_bb_vo("e,m") -= 0.500000 * V_blks_["bbbb_oovo"]("j,i,a,m") * u2_bbbb_vvoo("a,e,j,i");
        }

        if (include_u1_) {

            // ru1_bb_vo += 1.000000 <i,e||a,m>_abab u1_aa(a,i)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_bb_vo("e,m") -= V_blks_["baab_vovo"]("e,i,a,m") * u1_aa_vo("a,i");

            // ru1_bb_vo += -1.000000 d-_bb(i,a) u1_bb(a,i) u1_bb(e,m)
            // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1,
            ru1_bb_vo("e,m") -= scalar9 * u1_bb_vo("e,m");

            // ru1_bb_vo += 1.000000 <i,e||a,m>_bbbb u1_bb(a,i)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_bb_vo("e,m") -= V_blks_["bbbb_vovo"]("e,i,a,m") * u1_bb_vo("a,i");

            // ru1_bb_vo += 1.000000 f_bb(e,a) u1_bb(a,m)
            // flops: o1v2: 1, o1v1: 1 | mem: o1v1: 2,
            ru1_bb_vo("e,m") += F_blks_["bb_vv"]("e,a") * u1_bb_vo("a,m");
        }

        if (include_u1_ && include_u2_) {

            // ru2_aaaa_vvoo += -1.000000 d-_aa(i,a) u1_aa(a,i) u2_aaaa(e,f,m,n)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            ru2_aaaa_vvoo("e,f,m,n") -= scalar5 * u2_aaaa_vvoo("e,f,m,n");

            // ru2_aaaa_vvoo += -1.000000 d-_bb(i,a) u1_bb(a,i) u2_aaaa(e,f,m,n)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            ru2_aaaa_vvoo("e,f,m,n") -= scalar9 * u2_aaaa_vvoo("e,f,m,n");
        }

        if (include_u2_) {

            // ru2_aaaa_vvoo += 0.500000 <j,i||m,n>_aaaa u2_aaaa(e,f,j,i)
            // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_aaaa_vvoo("e,f,m,n") += 0.500000 * V_blks_["aaaa_oooo"]("j,i,m,n") * u2_aaaa_vvoo("e,f,j,i");

            // ru2_aaaa_vvoo += 0.250000 <j,i||a,b>_aaaa t2_aaaa(e,f,j,i) u2_aaaa(a,b,m,n)
            // flops: o4v2: 2, o2v2: 1 | mem: o2v2: 2, o4v0: 1,
            ru2_aaaa_vvoo("e,f,m,n") += 0.250000 * V_blks_["aaaa_oovv"]("j,i,a,b") * u2_aaaa_vvoo("a,b,m,n") * t2_aaaa_vvoo("e,f,j,i");

            // ru2_aaaa_vvoo += 1.000000 u2_aaaa(e,f,m,n) w0
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            ru2_aaaa_vvoo("e,f,m,n") += u2_aaaa_vvoo("e,f,m,n") * w0;

            // ru2_aaaa_vvoo += 0.500000 <e,f||a,b>_aaaa u2_aaaa(a,b,m,n)
            // flops: o2v4: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_aaaa_vvoo("e,f,m,n") += 0.500000 * V_blks_["aaaa_vvvv"]("e,f,a,b") * u2_aaaa_vvoo("a,b,m,n");
        }

        {

            // rt2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("f,e,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += 0.500000 P(e,f) <i,e||a,b>_aaaa t2_aaaa(a,b,m,n) u1_aa(f,i)
            // flops: o3v3: 1, o3v2: 1, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = -0.500000 * t2_aaaa_vvoo("a,b,m,n") * V_blks_["aaaa_vovv"]("e,i,a,b") * u1_aa_vo("f,i");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += 2.000000 P(e,f) d-_aa(i,a) u1_aa(e,i) u2_aaaa(a,f,m,n)
            // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") += 2.000000 * u2_aaaa_vvoo("a,f,m,n") * dp_aa_ov("i,a") * u1_aa_vo("e,i");
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(e,f) <i,e||m,n>_aaaa u1_aa(f,i)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= V_blks_["aaaa_vooo"]("e,i,m,n") * u1_aa_vo("f,i");

            // tempPerm_aaaa_vvoo += 1.000000 P(e,f) f_aa(i,a) t2_aaaa(a,e,m,n) u1_aa(f,i)
            // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") += t2_aaaa_vvoo("a,e,m,n") * F_blks_["aa_ov"]("i,a") * u1_aa_vo("f,i");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(e,f) f_aa(e,a) u2_aaaa(a,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += F_blks_["aa_vv"]("e,a") * u2_aaaa_vvoo("a,f,m,n");
        }

        {

            // tempPerm_aaaa_vvoo += -1.000000 P(e,f) d+_aa(e,a) t2_aaaa(a,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= dp_aa_vv("e,a") * t2_aaaa_vvoo("a,f,m,n");
        }

        if (include_u2_) {

            // ru2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) <i,j||n,a>_abab t2_aaaa(e,f,m,i) u1_bb(a,j)
            // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o2v2: 2, o2v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") = u1_bb_vo("a,j") * V_blks_["abba_oovo"]("i,j,a,n") * t2_aaaa_vvoo("e,f,m,i");
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) f_aa(i,n) u2_aaaa(e,f,m,i)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= F_blks_["aa_oo"]("i,n") * u2_aaaa_vvoo("e,f,m,i");
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) f_aa(i,a) t2_aaaa(e,f,n,i) u1_aa(a,m)
            // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") += u1_aa_vo("a,m") * F_blks_["aa_ov"]("i,a") * t2_aaaa_vvoo("e,f,n,i");

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) <e,f||a,n>_aaaa u1_aa(a,m)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += V_blks_["aaaa_vvvo"]("e,f,a,n") * u1_aa_vo("a,m");

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) <j,i||a,n>_aaaa t2_aaaa(e,f,m,i) u1_aa(a,j)
            // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o2v2: 2, o2v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") -= u1_aa_vo("a,j") * V_blks_["aaaa_oovo"]("j,i,a,n") * t2_aaaa_vvoo("e,f,m,i");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_aaaa_vvoo += 2.000000 P(m,n) d-_aa(i,a) u1_aa(a,n) u2_aaaa(e,f,m,i)
            // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") += 2.000000 * u1_aa_vo("a,n") * dp_aa_ov("i,a") * u2_aaaa_vvoo("e,f,m,i");
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += 0.500000 P(m,n) <j,i||a,n>_aaaa t2_aaaa(e,f,j,i) u1_aa(a,m)
            // flops: o4v2: 1, o4v1: 1, o2v2: 1 | mem: o2v2: 2, o4v0: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") += 0.500000 * u1_aa_vo("a,m") * V_blks_["aaaa_oovo"]("j,i,a,n") * t2_aaaa_vvoo("e,f,j,i");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) d+_aa(i,n) t2_aaaa(e,f,m,i)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += dp_aa_oo("i,n") * t2_aaaa_vvoo("e,f,m,i");
        }

        if (include_u2_) {

            // ru2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u2_) {

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <i,e||a,n>_aaaa u2_aaaa(a,f,m,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * V_blks_["aaaa_vovo"]("e,i,a,n") * u2_aaaa_vvoo("a,f,m,i");

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) <e,i||n,a>_abab u2_abab(f,a,m,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += V_blks_["abba_vovo"]("e,i,a,n") * u2_abab_vvoo("f,a,m,i");
        }

        if (include_u1_) {

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) <j,i||a,n>_aaaa t2_aaaa(a,e,m,i) u1_aa(f,j)
            // flops: o4v2: 1, o3v2: 1, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") -= t2_aaaa_vvoo("a,e,m,i") * V_blks_["aaaa_oovo"]("j,i,a,n") * u1_aa_vo("f,j");

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) <j,i||n,a>_abab t2_abab(e,a,m,i) u1_aa(f,j)
            // flops: o4v2: 1, o3v2: 1, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            tempPerm_aaaa_vvoo("e,f,m,n") += t2_abab_vvoo("e,a,m,i") * V_blks_["abba_oovo"]("j,i,a,n") * u1_aa_vo("f,j");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += 1.000000 f_aa(e,a) u2_abab(a,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += F_blks_["aa_vv"]("e,a") * u2_abab_vvoo("a,f,m,n");
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += -1.000000 <j,i||a,n>_abab t2_abab(e,f,m,i) u1_aa(a,j)
            // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o2v2: 2, o2v0: 1,
            ru2_abab_vvoo("e,f,m,n") -= V_blks_["abab_oovo"]("j,i,a,n") * u1_aa_vo("a,j") * t2_abab_vvoo("e,f,m,i");

            // ru2_abab_vvoo += 1.000000 <e,f||m,a>_abab u1_bb(a,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= V_blks_["abba_vvvo"]("e,f,a,m") * u1_bb_vo("a,n");

            // ru2_abab_vvoo += 1.000000 <j,i||a,m>_aaaa t2_abab(a,f,i,n) u1_aa(e,j)
            // flops: o4v2: 1, o3v2: 1, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            ru2_abab_vvoo("e,f,m,n") += V_blks_["aaaa_oovo"]("j,i,a,m") * t2_abab_vvoo("a,f,i,n") * u1_aa_vo("e,j");

            // ru2_abab_vvoo += 1.000000 <j,i||a,n>_abab t2_abab(a,f,m,i) u1_aa(e,j)
            // flops: o4v2: 1, o3v2: 1, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            ru2_abab_vvoo("e,f,m,n") += V_blks_["abab_oovo"]("j,i,a,n") * t2_abab_vvoo("a,f,m,i") * u1_aa_vo("e,j");

            // ru2_abab_vvoo += -1.000000 <j,i||a,n>_bbbb t2_abab(e,f,m,i) u1_bb(a,j)
            // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o2v2: 2, o2v0: 1,
            ru2_abab_vvoo("e,f,m,n") -= V_blks_["bbbb_oovo"]("j,i,a,n") * u1_bb_vo("a,j") * t2_abab_vvoo("e,f,m,i");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += 0.500000 <e,f||a,b>_abab u2_abab(a,b,m,n)
            // +                0.500000 <e,f||b,a>_abab u2_abab(b,a,m,n)
            // flops: o2v4: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += V_blks_["abab_vvvv"]("e,f,a,b") * u2_abab_vvoo("a,b,m,n");

            // ru2_abab_vvoo += -1.000000 <i,f||m,a>_abab u2_abab(e,a,i,n)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= V_blks_["baba_vovo"]("f,i,a,m") * u2_abab_vvoo("e,a,i,n");

            // ru2_abab_vvoo += -1.000000 d+_aa(e,a) t2_abab(a,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= dp_aa_vv("e,a") * t2_abab_vvoo("a,f,m,n");

            // ru2_abab_vvoo += -1.000000 <i,f||a,n>_abab u2_aaaa(a,e,m,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += V_blks_["baab_vovo"]("f,i,a,n") * u2_aaaa_vvoo("a,e,m,i");

            // ru2_abab_vvoo += -1.000000 <e,i||a,n>_abab u2_abab(a,f,m,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= V_blks_["abab_vovo"]("e,i,a,n") * u2_abab_vvoo("a,f,m,i");
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += 1.000000 <j,i||a,n>_bbbb t2_abab(e,a,m,i) u1_bb(f,j)
            // flops: o4v2: 1, o3v2: 1, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            ru2_abab_vvoo("e,f,m,n") += V_blks_["bbbb_oovo"]("j,i,a,n") * t2_abab_vvoo("e,a,m,i") * u1_bb_vo("f,j");

            // ru2_abab_vvoo += -1.000000 <i,j||m,a>_abab t2_abab(e,f,i,n) u1_bb(a,j)
            // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o2v2: 2, o2v0: 1,
            ru2_abab_vvoo("e,f,m,n") += V_blks_["abba_oovo"]("i,j,a,m") * u1_bb_vo("a,j") * t2_abab_vvoo("e,f,i,n");

            // ru2_abab_vvoo += -1.000000 d-_bb(i,a) u1_bb(a,i) u2_abab(e,f,m,n)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            ru2_abab_vvoo("e,f,m,n") -= scalar9 * u2_abab_vvoo("e,f,m,n");

            // ru2_abab_vvoo += 0.500000 <j,i||a,n>_abab t2_abab(e,f,j,i) u1_aa(a,m)
            // +                0.500000 <i,j||a,n>_abab t2_abab(e,f,i,j) u1_aa(a,m)
            // flops: o4v2: 1, o4v1: 1, o2v2: 1 | mem: o2v2: 2, o4v0: 1,
            ru2_abab_vvoo("e,f,m,n") += V_blks_["abab_oovo"]("j,i,a,n") * u1_aa_vo("a,m") * t2_abab_vvoo("e,f,j,i");

            // ru2_abab_vvoo += -1.000000 f_aa(i,a) t2_abab(e,f,i,n) u1_aa(a,m)
            // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
            ru2_abab_vvoo("e,f,m,n") -= F_blks_["aa_ov"]("i,a") * u1_aa_vo("a,m") * t2_abab_vvoo("e,f,i,n");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += 0.250000 <j,i||a,b>_abab t2_abab(e,f,j,i) u2_abab(a,b,m,n)
            // +                0.250000 <j,i||b,a>_abab t2_abab(e,f,j,i) u2_abab(b,a,m,n)
            // +                0.250000 <i,j||a,b>_abab t2_abab(e,f,i,j) u2_abab(a,b,m,n)
            // +                0.250000 <i,j||b,a>_abab t2_abab(e,f,i,j) u2_abab(b,a,m,n)
            // flops: o4v2: 2, o2v2: 1 | mem: o2v2: 2, o4v0: 1,
            ru2_abab_vvoo("e,f,m,n") += V_blks_["abab_oovv"]("j,i,a,b") * u2_abab_vvoo("a,b,m,n") * t2_abab_vvoo("e,f,j,i");
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += 2.000000 d-_aa(i,a) u1_aa(e,i) u2_abab(a,f,m,n)
            // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            ru2_abab_vvoo("e,f,m,n") += 2.000000 * dp_aa_ov("i,a") * u2_abab_vvoo("a,f,m,n") * u1_aa_vo("e,i");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += 1.000000 f_bb(f,a) u2_abab(e,a,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += F_blks_["bb_vv"]("f,a") * u2_abab_vvoo("e,a,m,n");

            // ru2_abab_vvoo += 1.000000 u2_abab(e,f,m,n) w0
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            ru2_abab_vvoo("e,f,m,n") += u2_abab_vvoo("e,f,m,n") * w0;
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += 2.000000 d-_bb(i,a) u1_bb(f,i) u2_abab(e,a,m,n)
            // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            ru2_abab_vvoo("e,f,m,n") += 2.000000 * dp_bb_ov("i,a") * u2_abab_vvoo("e,a,m,n") * u1_bb_vo("f,i");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += 1.000000 <i,e||a,m>_aaaa u2_abab(a,f,i,n)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= V_blks_["aaaa_vovo"]("e,i,a,m") * u2_abab_vvoo("a,f,i,n");
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += 1.000000 <i,j||a,n>_abab t2_aaaa(a,e,m,i) u1_bb(f,j)
            // flops: o4v2: 1, o3v2: 1, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            ru2_abab_vvoo("e,f,m,n") += V_blks_["abab_oovo"]("i,j,a,n") * t2_aaaa_vvoo("a,e,m,i") * u1_bb_vo("f,j");

            // ru2_abab_vvoo += -1.000000 f_aa(i,a) t2_abab(a,f,m,n) u1_aa(e,i)
            // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            ru2_abab_vvoo("e,f,m,n") -= F_blks_["aa_ov"]("i,a") * t2_abab_vvoo("a,f,m,n") * u1_aa_vo("e,i");

            // ru2_abab_vvoo += 2.000000 d-_aa(i,a) u1_aa(a,m) u2_abab(e,f,i,n)
            // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
            ru2_abab_vvoo("e,f,m,n") += 2.000000 * dp_aa_ov("i,a") * u1_aa_vo("a,m") * u2_abab_vvoo("e,f,i,n");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += -1.000000 <e,i||m,a>_abab u2_bbbb(a,f,n,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += V_blks_["abba_vovo"]("e,i,a,m") * u2_bbbb_vvoo("a,f,n,i");

            // ru2_abab_vvoo += 1.000000 <j,i||a,b>_abab t2_abab(a,f,m,i) u2_abab(e,b,j,n)
            // flops: o3v3: 2, o2v2: 1 | mem: o2v2: 3,
            ru2_abab_vvoo("e,f,m,n") += V_blks_["abab_oovv"]("j,i,a,b") * t2_abab_vvoo("a,f,m,i") * u2_abab_vvoo("e,b,j,n");
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += -0.500000 <i,f||a,b>_abab t2_abab(a,b,m,n) u1_aa(e,i)
            // +                -0.500000 <i,f||b,a>_abab t2_abab(b,a,m,n) u1_aa(e,i)
            // flops: o3v3: 1, o3v2: 1, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            ru2_abab_vvoo("e,f,m,n") += V_blks_["baab_vovv"]("f,i,a,b") * t2_abab_vvoo("a,b,m,n") * u1_aa_vo("e,i");

            // ru2_abab_vvoo += 0.500000 <j,i||m,a>_abab t2_abab(e,f,j,i) u1_bb(a,n)
            // +                0.500000 <i,j||m,a>_abab t2_abab(e,f,i,j) u1_bb(a,n)
            // flops: o4v2: 1, o4v1: 1, o2v2: 1 | mem: o2v2: 2, o4v0: 1,
            ru2_abab_vvoo("e,f,m,n") -= V_blks_["abba_oovo"]("j,i,a,m") * u1_bb_vo("a,n") * t2_abab_vvoo("e,f,j,i");

            // ru2_abab_vvoo += -1.000000 <j,i||a,m>_aaaa t2_abab(e,f,i,n) u1_aa(a,j)
            // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o2v2: 2, o2v0: 1,
            ru2_abab_vvoo("e,f,m,n") -= V_blks_["aaaa_oovo"]("j,i,a,m") * u1_aa_vo("a,j") * t2_abab_vvoo("e,f,i,n");

            // ru2_abab_vvoo += -1.000000 f_bb(i,a) t2_abab(e,a,m,n) u1_bb(f,i)
            // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            ru2_abab_vvoo("e,f,m,n") -= F_blks_["bb_ov"]("i,a") * t2_abab_vvoo("e,a,m,n") * u1_bb_vo("f,i");

            // ru2_abab_vvoo += -0.500000 <e,i||a,b>_abab t2_abab(a,b,m,n) u1_bb(f,i)
            // +                -0.500000 <e,i||b,a>_abab t2_abab(b,a,m,n) u1_bb(f,i)
            // flops: o3v3: 1, o3v2: 1, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            ru2_abab_vvoo("e,f,m,n") -= V_blks_["abab_vovv"]("e,i,a,b") * t2_abab_vvoo("a,b,m,n") * u1_bb_vo("f,i");

            // ru2_abab_vvoo += -1.000000 <e,i||m,n>_abab u1_bb(f,i)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= V_blks_["abab_vooo"]("e,i,m,n") * u1_bb_vo("f,i");

            // ru2_abab_vvoo += -1.000000 <e,i||a,b>_abab t2_abab(a,f,m,i) u1_bb(b,n)
            // flops: o3v3: 1, o2v3: 1, o2v2: 1 | mem: o2v2: 3,
            ru2_abab_vvoo("e,f,m,n") -= V_blks_["abab_vovv"]("e,i,a,b") * u1_bb_vo("b,n") * t2_abab_vvoo("a,f,m,i");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += 1.000000 <i,f||a,n>_bbbb u2_abab(e,a,m,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= V_blks_["bbbb_vovo"]("f,i,a,n") * u2_abab_vvoo("e,a,m,i");

            // ru2_abab_vvoo += 0.500000 <j,i||m,n>_abab u2_abab(e,f,j,i)
            // +                0.500000 <i,j||m,n>_abab u2_abab(e,f,i,j)
            // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += V_blks_["abab_oooo"]("j,i,m,n") * u2_abab_vvoo("e,f,j,i");
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += -1.000000 f_bb(i,a) t2_abab(e,f,m,i) u1_bb(a,n)
            // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
            ru2_abab_vvoo("e,f,m,n") -= F_blks_["bb_ov"]("i,a") * u1_bb_vo("a,n") * t2_abab_vvoo("e,f,m,i");

            // ru2_abab_vvoo += 2.000000 d-_bb(i,a) u1_bb(a,n) u2_abab(e,f,m,i)
            // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
            ru2_abab_vvoo("e,f,m,n") += 2.000000 * dp_bb_ov("i,a") * u1_bb_vo("a,n") * u2_abab_vvoo("e,f,m,i");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += -1.000000 d+_bb(f,a) t2_abab(e,a,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= dp_bb_vv("f,a") * t2_abab_vvoo("e,a,m,n");

            // ru2_abab_vvoo += 1.000000 d+_bb(i,n) t2_abab(e,f,m,i)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += dp_bb_oo("i,n") * t2_abab_vvoo("e,f,m,i");
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += 1.000000 <e,f||a,n>_abab u1_aa(a,m)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += V_blks_["abab_vvvo"]("e,f,a,n") * u1_aa_vo("a,m");

            // ru2_abab_vvoo += -1.000000 <i,f||m,n>_abab u1_aa(e,i)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += V_blks_["baab_vooo"]("f,i,m,n") * u1_aa_vo("e,i");

            // ru2_abab_vvoo += 1.000000 <i,j||m,a>_abab t2_abab(e,a,i,n) u1_bb(f,j)
            // flops: o4v2: 1, o3v2: 1, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            ru2_abab_vvoo("e,f,m,n") -= V_blks_["abba_oovo"]("i,j,a,m") * t2_abab_vvoo("e,a,i,n") * u1_bb_vo("f,j");

            // ru2_abab_vvoo += -1.000000 <i,f||b,a>_abab t2_abab(e,a,i,n) u1_aa(b,m)
            // flops: o3v3: 1, o2v3: 1, o2v2: 1 | mem: o2v2: 3,
            ru2_abab_vvoo("e,f,m,n") += V_blks_["baab_vovv"]("f,i,b,a") * u1_aa_vo("b,m") * t2_abab_vvoo("e,a,i,n");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += 1.000000 d+_aa(i,m) t2_abab(e,f,i,n)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") += dp_aa_oo("i,m") * t2_abab_vvoo("e,f,i,n");
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += -1.000000 d-_aa(i,a) u1_aa(a,i) u2_abab(e,f,m,n)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            ru2_abab_vvoo("e,f,m,n") -= scalar5 * u2_abab_vvoo("e,f,m,n");
        }

        if (include_u2_) {

            // ru2_abab_vvoo += -1.000000 f_aa(i,m) u2_abab(e,f,i,n)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= F_blks_["aa_oo"]("i,m") * u2_abab_vvoo("e,f,i,n");

            // ru2_abab_vvoo += -1.000000 f_bb(i,n) u2_abab(e,f,m,i)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_abab_vvoo("e,f,m,n") -= F_blks_["bb_oo"]("i,n") * u2_abab_vvoo("e,f,m,i");
        }

        if (include_u1_ && include_u2_) {

            // ru2_abab_vvoo += 1.000000 <j,i||m,a>_abab t2_bbbb(a,f,n,i) u1_aa(e,j)
            // flops: o4v2: 1, o3v2: 1, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            ru2_abab_vvoo("e,f,m,n") -= V_blks_["abba_oovo"]("j,i,a,m") * t2_bbbb_vvoo("a,f,n,i") * u1_aa_vo("e,j");
        }

        if (include_u2_) {

            // ru2_bbbb_vvoo += 1.000000 u2_bbbb(e,f,m,n) w0
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            ru2_bbbb_vvoo("e,f,m,n") += u2_bbbb_vvoo("e,f,m,n") * w0;

            // ru2_bbbb_vvoo += 0.250000 <j,i||a,b>_bbbb t2_bbbb(e,f,j,i) u2_bbbb(a,b,m,n)
            // flops: o4v2: 2, o2v2: 1 | mem: o2v2: 2, o4v0: 1,
            ru2_bbbb_vvoo("e,f,m,n") += 0.250000 * V_blks_["bbbb_oovv"]("j,i,a,b") * u2_bbbb_vvoo("a,b,m,n") * t2_bbbb_vvoo("e,f,j,i");
        }

        if (include_u1_ && include_u2_) {

            // ru2_bbbb_vvoo += -1.000000 d-_bb(i,a) u1_bb(a,i) u2_bbbb(e,f,m,n)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            ru2_bbbb_vvoo("e,f,m,n") -= scalar9 * u2_bbbb_vvoo("e,f,m,n");

            // ru2_bbbb_vvoo += -1.000000 d-_aa(i,a) u1_aa(a,i) u2_bbbb(e,f,m,n)
            // flops: o2v2: 1, o0v0: 1 | mem: o2v2: 1, o0v0: 1,
            ru2_bbbb_vvoo("e,f,m,n") -= scalar5 * u2_bbbb_vvoo("e,f,m,n");
        }

        if (include_u2_) {

            // ru2_bbbb_vvoo += 0.500000 <j,i||m,n>_bbbb u2_bbbb(e,f,j,i)
            // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_bbbb_vvoo("e,f,m,n") += 0.500000 * V_blks_["bbbb_oooo"]("j,i,m,n") * u2_bbbb_vvoo("e,f,j,i");

            // ru2_bbbb_vvoo += 0.500000 <e,f||a,b>_bbbb u2_bbbb(a,b,m,n)
            // flops: o2v4: 1, o2v2: 1 | mem: o2v2: 2,
            ru2_bbbb_vvoo("e,f,m,n") += 0.500000 * V_blks_["bbbb_vvvv"]("e,f,a,b") * u2_bbbb_vvoo("a,b,m,n");
        }

        {

            // rt2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("f,e,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_bbbb_vvoo += 2.000000 P(e,f) d-_bb(i,a) u1_bb(e,i) u2_bbbb(a,f,m,n)
            // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = 2.000000 * u2_bbbb_vvoo("a,f,m,n") * dp_bb_ov("i,a") * u1_bb_vo("e,i");
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(e,f) <i,e||m,n>_bbbb u1_bb(f,i)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= V_blks_["bbbb_vooo"]("e,i,m,n") * u1_bb_vo("f,i");
        }

        {

            // tempPerm_bbbb_vvoo += -1.000000 P(e,f) d+_bb(e,a) t2_bbbb(a,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= dp_bb_vv("e,a") * t2_bbbb_vvoo("a,f,m,n");
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(e,f) f_bb(e,a) u2_bbbb(a,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += F_blks_["bb_vv"]("e,a") * u2_bbbb_vvoo("a,f,m,n");
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(e,f) f_bb(i,a) t2_bbbb(a,e,m,n) u1_bb(f,i)
            // flops: o3v2: 2, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") += t2_bbbb_vvoo("a,e,m,n") * F_blks_["bb_ov"]("i,a") * u1_bb_vo("f,i");

            // tempPerm_bbbb_vvoo += 0.500000 P(e,f) <i,e||a,b>_bbbb t2_bbbb(a,b,m,n) u1_bb(f,i)
            // flops: o3v3: 1, o3v2: 1, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") -= 0.500000 * t2_bbbb_vvoo("a,b,m,n") * V_blks_["bbbb_vovv"]("e,i,a,b") * u1_bb_vo("f,i");
        }

        if (include_u2_) {

            // ru2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) d+_bb(i,n) t2_bbbb(e,f,m,i)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = dp_bb_oo("i,n") * t2_bbbb_vvoo("e,f,m,i");
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) <j,i||a,n>_abab t2_bbbb(e,f,m,i) u1_aa(a,j)
            // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o2v2: 2, o2v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") -= u1_aa_vo("a,j") * V_blks_["abab_oovo"]("j,i,a,n") * t2_bbbb_vvoo("e,f,m,i");

            // tempPerm_bbbb_vvoo += 0.500000 P(m,n) <j,i||a,n>_bbbb t2_bbbb(e,f,j,i) u1_bb(a,m)
            // flops: o4v2: 1, o4v1: 1, o2v2: 1 | mem: o2v2: 2, o4v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") += 0.500000 * u1_bb_vo("a,m") * V_blks_["bbbb_oovo"]("j,i,a,n") * t2_bbbb_vvoo("e,f,j,i");
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) f_bb(i,n) u2_bbbb(e,f,m,i)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= F_blks_["bb_oo"]("i,n") * u2_bbbb_vvoo("e,f,m,i");
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) <j,i||a,n>_bbbb t2_bbbb(e,f,m,i) u1_bb(a,j)
            // flops: o3v2: 1, o2v2: 1, o3v1: 1 | mem: o2v2: 2, o2v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") -= u1_bb_vo("a,j") * V_blks_["bbbb_oovo"]("j,i,a,n") * t2_bbbb_vvoo("e,f,m,i");
        }

        if (include_u1_ && include_u2_) {

            // tempPerm_bbbb_vvoo += 2.000000 P(m,n) d-_bb(i,a) u1_bb(a,n) u2_bbbb(e,f,m,i)
            // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") += 2.000000 * u1_bb_vo("a,n") * dp_bb_ov("i,a") * u2_bbbb_vvoo("e,f,m,i");
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) f_bb(i,a) t2_bbbb(e,f,n,i) u1_bb(a,m)
            // flops: o3v2: 1, o2v2: 1, o2v1: 1 | mem: o2v2: 2, o2v0: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") += u1_bb_vo("a,m") * F_blks_["bb_ov"]("i,a") * t2_bbbb_vvoo("e,f,n,i");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) <e,f||a,n>_bbbb u1_bb(a,m)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += V_blks_["bbbb_vvvo"]("e,f,a,n") * u1_bb_vo("a,m");
        }

        if (include_u2_) {

            // ru2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

        if (include_u1_) {

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) <i,j||a,n>_abab t2_abab(a,e,i,m) u1_bb(f,j)
            // flops: o4v2: 1, o3v2: 1, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * t2_abab_vvoo("a,e,i,m") * V_blks_["abab_oovo"]("i,j,a,n") * u1_bb_vo("f,j");

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) <j,i||a,n>_bbbb t2_bbbb(a,e,m,i) u1_bb(f,j)
            // flops: o4v2: 1, o3v2: 1, o2v2: 1 | mem: o2v2: 2, o3v1: 1,
            tempPerm_bbbb_vvoo("e,f,m,n") -= t2_bbbb_vvoo("a,e,m,i") * V_blks_["bbbb_oovo"]("j,i,a,n") * u1_bb_vo("f,j");
        }

        if (include_u2_) {

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <i,e||a,n>_bbbb u2_bbbb(a,f,m,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= V_blks_["bbbb_vovo"]("e,i,a,n") * u2_bbbb_vvoo("a,f,m,i");

            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) <i,e||a,n>_abab u2_abab(a,f,i,m)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += V_blks_["baab_vovo"]("e,i,a,n") * u2_abab_vvoo("a,f,i,m");
        }

        if (include_u0_) {

            // cenergy+= 1.000000 u0
            // flops: o0v0: 1 | mem: o0v0: 1,
            cenergy+= u0;
        }

        if (include_u1_) {

            // crt1_aa_vo += 1.000000 u1_aa(e,m)
            // flops: o1v1: 1 | mem: o1v1: 1,
            crt1_aa_vo("e,m") += u1_aa_vo("e,m");

            // crt1_bb_vo += 1.000000 u1_bb(e,m)
            // flops: o1v1: 1 | mem: o1v1: 1,
            crt1_bb_vo("e,m") += u1_bb_vo("e,m");
        }

        if (include_u2_) {

            // crt2_aaaa_vvoo += 1.000000 u2_aaaa(e,f,m,n)
            // flops: o2v2: 1 | mem: o2v2: 1,
            crt2_aaaa_vvoo("e,f,m,n") += u2_aaaa_vvoo("e,f,m,n");

            // crt2_abab_vvoo += 1.000000 u2_abab(e,f,m,n)
            // flops: o2v2: 1 | mem: o2v2: 1,
            crt2_abab_vvoo("e,f,m,n") += u2_abab_vvoo("e,f,m,n");

            // crt2_bbbb_vvoo += 1.000000 u2_bbbb(e,f,m,n)
            // flops: o2v2: 1 | mem: o2v2: 1,
            crt2_bbbb_vvoo("e,f,m,n") += u2_bbbb_vvoo("e,f,m,n");
        }

        if (include_u0_) {

            // cru0 += 1.000000
            // flops: o0v0: 1 | mem: o0v0: 1,
            cru0 += 1.000000;
        }

        if (include_u2_) {

            // ru2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo tempPerm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
            ru2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
            ru2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("f,e,n,m");
        }

        {

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo tempPerm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_aaaa_vvoo) TArrayD();
        }

        if (include_u2_) {

            // ru2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo tempPerm
            // flops: o2v2: 1 | mem: o2v2: 1,
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
            ru2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
            ru2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("f,e,n,m");
        }

        {

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo tempPerm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo.~TArrayD(); world_.gop.fence(); new (&tempPerm_bbbb_vvoo) TArrayD();
        }

/*

        Total Number of Terms: 617
        Number of Flops: (old) 1425 -> (new) 1123

            Total FLOP scaling:
            ------------------
         Scaling :      new |      old |     diff
        -------- : -------- | -------- | --------
            o2v4 :        6 |        8 |       -2
            o3v3 :       75 |       98 |      -23
            o4v2 :       35 |       60 |      -25
            o2v3 :       84 |      142 |      -58
            o3v2 :      130 |      202 |      -72
            o4v1 :        4 |        6 |       -2
            o1v3 :        4 |        8 |       -4
            o2v2 :      518 |      644 |     -126
            o3v1 :        8 |        8 |        0
            o4v0 :        3 |        0 |        3
            o1v2 :        4 |       10 |       -6
            o2v1 :       24 |       36 |      -12
            o0v2 :       13 |        0 |       13
            o1v1 :      108 |      104 |        4
            o2v0 :        8 |        0 |        8
            o0v0 :       99 |       99 |        0
    
            Total MEM scaling:
            o2v2 :      738 |      716 |       22
            o3v1 :       30 |       38 |       -8
            o4v0 :       13 |       24 |      -11
            o0v2 :       27 |       50 |      -23
            o1v1 :      178 |      222 |      -44
            o2v0 :       38 |       68 |      -30
            o0v0 :       99 |       99 |        0

*/


        double coherent_scalar;
        if ( options_.get_bool("QED_USE_RELAXED_ORBITALS")) {
            coherent_scalar = coupling_factor_z * e_dip_z_;
        } else {
            coherent_scalar = -coupling_factor_z * nuc_dip_z_;
        }

        energy += coherent_scalar * cenergy;
        rt1_aa_vo(idx_map_[2]) += coherent_scalar * crt1_aa_vo(idx_map_[2]);
        rt1_bb_vo(idx_map_[2]) += coherent_scalar * crt1_bb_vo(idx_map_[2]);
        rt2_aaaa_vvoo(idx_map_[4]) += coherent_scalar * crt2_aaaa_vvoo(idx_map_[4]);
        rt2_abab_vvoo(idx_map_[4]) += coherent_scalar * crt2_abab_vvoo(idx_map_[4]);
        rt2_bbbb_vvoo(idx_map_[4]) += coherent_scalar * crt2_bbbb_vvoo(idx_map_[4]);
        if (include_u0_) ru0 = ru0 + coherent_scalar * cru0;
        if (include_u1_) {
            ru1_aa_vo(idx_map_[2]) += coherent_scalar * cru1_aa_vo(idx_map_[2]);
            ru1_bb_vo(idx_map_[2]) += coherent_scalar * cru1_bb_vo(idx_map_[2]);
        }
        if (include_u2_) {
            ru2_aaaa_vvoo(idx_map_[4]) += coherent_scalar * cru2_aaaa_vvoo(idx_map_[4]);
            ru2_abab_vvoo(idx_map_[4]) += coherent_scalar * cru2_abab_vvoo(idx_map_[4]);
            ru2_bbbb_vvoo(idx_map_[4]) += coherent_scalar * cru2_bbbb_vvoo(idx_map_[4]);
        }

        if ( options_.get_bool("ZERO_U0") ){
            ru0 = 0.0;
            u0 = 0.0;
            include_u0_ = true;
        }

        scalar_resids_["u0"] = ru0;

    }
    world_.gop.fence();
    return energy + enuc_ + average_electric_dipole_self_energy_;
}
