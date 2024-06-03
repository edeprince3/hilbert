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

#include "cc_cavity/include/derived/qed_ccsd.h"


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
    std::map<std::string, TA::TArrayD> dp;
    
    dp["aa_oo"]("i, j") = coupling_factor_z * Dip_blks_["dz_aa_oo"]("i, j");
    dp["aa_ov"]("i, a") = coupling_factor_z * Dip_blks_["dz_aa_ov"]("i, a");
    dp["aa_vo"]("a, i") = coupling_factor_z * Dip_blks_["dz_aa_vo"]("a, i");
    dp["aa_vv"]("a, b") = coupling_factor_z * Dip_blks_["dz_aa_vv"]("a, b");

    dp["bb_oo"]("i, j") = coupling_factor_z * Dip_blks_["dz_bb_oo"]("i, j");
    dp["bb_ov"]("i, a") = coupling_factor_z * Dip_blks_["dz_bb_ov"]("i, a");
    dp["bb_vo"]("a, i") = coupling_factor_z * Dip_blks_["dz_bb_vo"]("a, i");
    dp["bb_vv"]("a, b") = coupling_factor_z * Dip_blks_["dz_bb_vv"]("a, b");

    // get residual tensors
    TArrayD &rt1_aa   = residuals_["t1_aa"];
    TArrayD &rt1_bb   = residuals_["t1_bb"];
    TArrayD &ru1_aa   = residuals_["u1_aa"];
    TArrayD &ru1_bb   = residuals_["u1_bb"];
    TArrayD &rt2_aaaa = residuals_["t2_aaaa"];
    TArrayD &rt2_abab = residuals_["t2_abab"];
    TArrayD &rt2_bbbb = residuals_["t2_bbbb"];
    TArrayD &ru2_aaaa = residuals_["u2_aaaa"];
    TArrayD &ru2_abab = residuals_["u2_abab"];
    TArrayD &ru2_bbbb = residuals_["u2_bbbb"];
    double ru0 = 0.0;

    // initialize to zero
    double energy = 0;
    zero_tiles(rt1_aa);
    zero_tiles(rt1_bb);
    zero_tiles(rt2_aaaa);
    zero_tiles(rt2_abab);
    zero_tiles(rt2_bbbb);

    if ( include_u1_){
        zero_tiles(ru1_aa);
        zero_tiles(ru1_bb);
    }
    if ( include_u2_){
        zero_tiles(ru2_aaaa);
        zero_tiles(ru2_abab);
        zero_tiles(ru2_bbbb);
    }

    // initialize coherent state amplitudes to zero
    TA::TArrayD crt1_aa = makeTensor(world_, {va_, oa_}, true);
    TA::TArrayD crt1_bb = makeTensor(world_, {vb_, ob_}, true);
    TA::TArrayD crt2_aaaa = makeTensor(world_, {va_, va_, oa_, oa_}, true);
    TA::TArrayD crt2_abab = makeTensor(world_, {va_, vb_, oa_, ob_}, true);
    TA::TArrayD crt2_bbbb = makeTensor(world_, {vb_, vb_, ob_, ob_}, true);

    double cenergy = 0.0;
    double cru0 = 0.0;
    TA::TArrayD cru1_aa;
    TA::TArrayD cru1_bb;
    TA::TArrayD cru2_aaaa;
    TA::TArrayD cru2_abab;
    TA::TArrayD cru2_bbbb;

    if (include_u1_){
        cru1_aa = makeTensor(world_, {va_, oa_}, true);
        cru1_bb = makeTensor(world_, {vb_, ob_}, true);
    }
    if (include_u2_){
        cru2_aaaa = makeTensor(world_, {va_, va_, oa_, oa_}, true);
        cru2_abab = makeTensor(world_, {va_, vb_, oa_, ob_}, true);
        cru2_bbbb = makeTensor(world_, {vb_, vb_, ob_, ob_}, true);
    }

    /// get reference amplitudes

    // extract scalar amplitude
    double& u0 = scalar_amps_["u0"];
    if(include_u0_) {
        world_.gop.fence();
        forall(amplitudes_["u0"], [&u0](auto &tile, auto &x){
            u0 = tile[x];
        });
        world_.gop.fence();
    }

    // extract 1-body amplitudes
    std::map<std::string, TA::TArrayD> t1 {
            {"aa_vo", amplitudes_["t1_aa"]},
            {"bb_vo", amplitudes_["t1_bb"]}
    };
    std::map<std::string, TA::TArrayD> u1 {
            {"aa_vo", amplitudes_["u1_aa"]},
            {"bb_vo", amplitudes_["u1_bb"]}
    };

    // extract 2-body amplitudes
    std::map<std::string, TA::TArrayD> t2 {
            {"aaaa_vvoo", amplitudes_["t2_aaaa"]},
            {"abab_vvoo", amplitudes_["t2_abab"]},
            {"bbbb_vvoo", amplitudes_["t2_bbbb"]}
    };
    std::map<std::string, TA::TArrayD> u2 {
            {"aaaa_vvoo", amplitudes_["u2_aaaa"]},
            {"abab_vvoo", amplitudes_["u2_abab"]},
            {"bbbb_vvoo", amplitudes_["u2_bbbb"]}
    };

    world_.gop.fence();

    // compute energy and residuals
    {

        std::map<std::string, double> scalars_;
        std::map<std::string, TA::TArrayD> tmps_;
        std::map<std::string, TA::TArrayD> perm_tmps;

    {
        scalars_["1"]  = 0.25 * dot(V_blks_["abab_oovv"]("i,j,a,b"), t2["abab_vvoo"]("a,b,i,j"));
        scalars_["2"]  = 0.25 * dot(V_blks_["abab_oovv"]("j,i,b,a"), u2["abab_vvoo"]("b,a,j,i"));
        scalars_["3"]  = 1.00 * dot(V_blks_["bbbb_oovv"]("j,i,a,b"), t2["bbbb_vvoo"]("a,b,j,i"));
        scalars_["4"]  = 1.00 * dot(V_blks_["aaaa_oovv"]("j,i,a,b"), t2["aaaa_vvoo"]("a,b,j,i"));
        scalars_["5"]  = 1.00 * dot(V_blks_["aaaa_oovv"]("j,i,a,b"), u2["aaaa_vvoo"]("a,b,j,i"));
        scalars_["6"]  = 1.00 * dot(V_blks_["bbbb_oovv"]("j,i,a,b"), u2["bbbb_vvoo"]("a,b,j,i"));
        scalars_["7"]  = 1.00 * dot(dp["aa_ov"]("i,a"), u1["aa_vo"]("a,i"));
        scalars_["8"]  = 1.00 * dot(dp["bb_ov"]("i,a"), u1["bb_vo"]("a,i"));
        scalars_["9"]  = 1.00 * dot(F_blks_["bb_ov"]("i,a"), u1["bb_vo"]("a,i"));
        scalars_["10"]  = 1.00 * dot(F_blks_["aa_ov"]("i,a"), u1["aa_vo"]("a,i"));
        scalars_["11"]  = 1.00 * dot(Id_blks_["bb_oo"]("i0,i1"), dp["bb_oo"]("i0,i1"));
        scalars_["12"]  = 1.00 * dot(Id_blks_["aa_oo"]("i0,i1"), dp["aa_oo"]("i0,i1"));
        scalars_["13"]  = 1.00 * dot(Id_blks_["aa_oo"]("i0,i1"), F_blks_["aa_oo"]("i0,i1"));
        scalars_["14"]  = 1.00 * dot(Id_blks_["bb_oo"]("i0,i1"), F_blks_["bb_oo"]("i0,i1"));
    }
 
    {
        // cenergy = +1.00 u0  // flops: o0v0 = o0v0 | mem: o0v0 = o0v0
        cenergy  = u0;
    
        // crt1_aa = +1.00 u1_aa(a,i)  // flops: o1v1 = o1v1 | mem: o1v1 = o1v1
        crt1_aa("a,i")  = u1["aa_vo"]("a,i");
    
        // crt1_bb = +1.00 u1_bb(a,i)  // flops: o1v1 = o1v1 | mem: o1v1 = o1v1
        crt1_bb("a,i")  = u1["bb_vo"]("a,i");
    
        // crt2_aaaa = +1.00 u2_aaaa(a,b,i,j)  // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
        crt2_aaaa("a,b,i,j")  = u2["aaaa_vvoo"]("a,b,i,j");
    
        // crt2_abab = +1.00 u2_abab(a,b,i,j)  // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
        crt2_abab("a,b,i,j")  = u2["abab_vvoo"]("a,b,i,j");
    
        // crt2_bbbb = +1.00 u2_bbbb(a,b,i,j)  // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
        crt2_bbbb("a,b,i,j")  = u2["bbbb_vvoo"]("a,b,i,j");
    
        // cru0 = +1.00  // flops: o0v0 | mem: o0v0
        cru0  = 1.00;
    
        // energy = +1.00 f_aa(i,i)  // flops: o0v0 = o0v0 | mem: o0v0 = o0v0
        energy  = scalars_["13"];
    
        // energy += +1.00 f_bb(i,i)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        energy += scalars_["14"];
    
        // energy += -1.00 d-_aa(i,i) u0  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        energy -= scalars_["12"] * u0;
    
        // energy += -1.00 d-_bb(i,i) u0  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        energy -= scalars_["11"] * u0;
    
        // energy += -1.00 d-_aa(i,a) u1_aa(a,i)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        energy -= scalars_["7"];
    
        // energy += -1.00 d-_bb(i,a) u1_bb(a,i)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        energy -= scalars_["8"];
    
        // energy += -0.50 <j,i||j,i>_aaaa  // flops: o0v0 += o4v0 o2v0 | mem: o0v0 += o2v0 o0v0
        energy -= 0.50 * dot(Id_blks_["aa_oo"]("i0,i1") * V_blks_["aaaa_oooo"]("j0,i0,j1,i1"), Id_blks_["aa_oo"]("j0,j1"));
    
        // energy += -0.50 <j,i||j,i>_bbbb  // flops: o0v0 += o4v0 o2v0 | mem: o0v0 += o2v0 o0v0
        energy -= 0.50 * dot(Id_blks_["bb_oo"]("i0,i1") * V_blks_["bbbb_oooo"]("j0,i0,j1,i1"), Id_blks_["bb_oo"]("j0,j1"));
    
        // energy += +0.25 <j,i||a,b>_aaaa t2_aaaa(a,b,j,i)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        energy += 0.25 * scalars_["4"];
    
        // energy += +0.25 <j,i||a,b>_abab t2_abab(a,b,j,i)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        energy += scalars_["1"];
    
        // energy += +0.25 <i,j||a,b>_abab t2_abab(a,b,i,j)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        energy += scalars_["1"];
    
        // energy += +0.25 <j,i||b,a>_abab t2_abab(b,a,j,i)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        energy += scalars_["1"];
    
        // energy += +0.25 <i,j||b,a>_abab t2_abab(b,a,i,j)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        energy += scalars_["1"];
    
        // energy += +0.25 <j,i||a,b>_bbbb t2_bbbb(a,b,j,i)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        energy += 0.25 * scalars_["3"];
    
        // rt1_aa = +1.00 f_aa(a,i)  // flops: o1v1 = o1v1 | mem: o1v1 = o1v1
        rt1_aa("a,i")  = F_blks_["aa_vo"]("a,i");
    
        // rt1_aa += -1.00 f_aa(j,b) t2_aaaa(b,a,i,j)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        rt1_aa("a,i") -= F_blks_["aa_ov"]("j,b") * t2["aaaa_vvoo"]("b,a,i,j");
    
        // rt1_aa += +1.00 f_bb(j,b) t2_abab(a,b,i,j)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        rt1_aa("a,i") += F_blks_["bb_ov"]("j,b") * t2["abab_vvoo"]("a,b,i,j");
    
        // rt1_aa += -1.00 d-_aa(a,i) u0  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_aa("a,i") -= dp["aa_vo"]("a,i") * u0;
    
        // rt1_aa += -1.00 d-_aa(j,j) u1_aa(a,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_aa("a,i") -= scalars_["12"] * u1["aa_vo"]("a,i");
    
        // rt1_aa += -1.00 d-_bb(j,j) u1_aa(a,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_aa("a,i") -= scalars_["11"] * u1["aa_vo"]("a,i");
    
        // rt1_aa += -0.50 <k,j||b,i>_aaaa t2_aaaa(b,a,k,j)  // flops: o1v1 += o3v2 | mem: o1v1 += o1v1
        rt1_aa("a,i") -= 0.50 * V_blks_["aaaa_oovo"]("k,j,b,i") * t2["aaaa_vvoo"]("b,a,k,j");
    
        // rt1_aa += -0.50 <j,a||b,c>_aaaa t2_aaaa(b,c,i,j)  // flops: o1v1 += o2v3 | mem: o1v1 += o1v1
        rt1_aa("a,i") += 0.50 * V_blks_["aaaa_vovv"]("a,j,b,c") * t2["aaaa_vvoo"]("b,c,i,j");
    
        // rt1_bb = +1.00 f_bb(a,i)  // flops: o1v1 = o1v1 | mem: o1v1 = o1v1
        rt1_bb("a,i")  = F_blks_["bb_vo"]("a,i");
    
        // rt1_bb += +1.00 f_aa(j,b) t2_abab(b,a,j,i)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        rt1_bb("a,i") += F_blks_["aa_ov"]("j,b") * t2["abab_vvoo"]("b,a,j,i");
    
        // rt1_bb += -1.00 f_bb(j,b) t2_bbbb(b,a,i,j)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        rt1_bb("a,i") -= F_blks_["bb_ov"]("j,b") * t2["bbbb_vvoo"]("b,a,i,j");
    
        // rt1_bb += -1.00 d-_bb(a,i) u0  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_bb("a,i") -= dp["bb_vo"]("a,i") * u0;
    
        // rt1_bb += -1.00 d-_aa(j,j) u1_bb(a,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_bb("a,i") -= scalars_["12"] * u1["bb_vo"]("a,i");
    
        // rt1_bb += -1.00 d-_bb(j,j) u1_bb(a,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_bb("a,i") -= scalars_["11"] * u1["bb_vo"]("a,i");
    
        // rt1_bb += -0.50 <k,j||b,i>_bbbb t2_bbbb(b,a,k,j)  // flops: o1v1 += o3v2 | mem: o1v1 += o1v1
        rt1_bb("a,i") -= 0.50 * V_blks_["bbbb_oovo"]("k,j,b,i") * t2["bbbb_vvoo"]("b,a,k,j");
    
        // rt1_bb += -0.50 <j,a||b,c>_bbbb t2_bbbb(b,c,i,j)  // flops: o1v1 += o2v3 | mem: o1v1 += o1v1
        rt1_bb("a,i") += 0.50 * V_blks_["bbbb_vovv"]("a,j,b,c") * t2["bbbb_vvoo"]("b,c,i,j");
    
        // rt2_aaaa = +1.00 P(i,j) P(a,b) d-_aa(a,j) u1_aa(b,i)  // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
        rt2_aaaa("a,b,i,j")  = dp["aa_vo"]("a,j") * u1["aa_vo"]("b,i");
        rt2_aaaa("a,b,i,j") -= dp["aa_vo"]("a,i") * u1["aa_vo"]("b,j");
        rt2_aaaa("a,b,i,j") -= dp["aa_vo"]("b,j") * u1["aa_vo"]("a,i");
        rt2_aaaa("a,b,i,j") += dp["aa_vo"]("b,i") * u1["aa_vo"]("a,j");
    
        // rt2_aaaa += -1.00 d-_aa(k,k) u2_aaaa(a,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") -= scalars_["12"] * u2["aaaa_vvoo"]("a,b,i,j");
    
        // rt2_aaaa += -1.00 d-_bb(k,k) u2_aaaa(a,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") -= scalars_["11"] * u2["aaaa_vvoo"]("a,b,i,j");
    
        // rt2_aaaa += +1.00 <a,b||i,j>_aaaa  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") += V_blks_["aaaa_vvoo"]("a,b,i,j");
    
        // rt2_aaaa += +0.50 <l,k||i,j>_aaaa t2_aaaa(a,b,l,k)  // flops: o2v2 += o4v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") += 0.50 * V_blks_["aaaa_oooo"]("l,k,i,j") * t2["aaaa_vvoo"]("a,b,l,k");
    
        // rt2_aaaa += +0.50 <a,b||c,d>_aaaa t2_aaaa(c,d,i,j)  // flops: o2v2 += o2v4 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") += 0.50 * V_blks_["aaaa_vvvv"]("a,b,c,d") * t2["aaaa_vvoo"]("c,d,i,j");
    
        // rt2_aaaa += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,a,i,j) t2_aaaa(d,b,l,k)  // flops: o2v2 += o2v3 o2v3 | mem: o2v2 += o0v2 o2v2
        rt2_aaaa("a,b,i,j") -= 0.50 * V_blks_["aaaa_oovv"]("l,k,c,d") * t2["aaaa_vvoo"]("d,b,l,k") * t2["aaaa_vvoo"]("c,a,i,j");
    
        // rt2_abab = -1.00 f_bb(k,j) t2_abab(a,b,i,k)  // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        rt2_abab("a,b,i,j")  = -1.00 * F_blks_["bb_oo"]("k,j") * t2["abab_vvoo"]("a,b,i,k");
    
        // rt2_abab += -1.00 f_aa(k,i) t2_abab(a,b,k,j)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= F_blks_["aa_oo"]("k,i") * t2["abab_vvoo"]("a,b,k,j");
    
        // rt2_abab += +1.00 f_aa(a,c) t2_abab(c,b,i,j)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += F_blks_["aa_vv"]("a,c") * t2["abab_vvoo"]("c,b,i,j");
    
        // rt2_abab += +1.00 f_bb(b,c) t2_abab(a,c,i,j)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += F_blks_["bb_vv"]("b,c") * t2["abab_vvoo"]("a,c,i,j");
    
        // rt2_abab += -1.00 d-_bb(b,j) u1_aa(a,i)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= dp["bb_vo"]("b,j") * u1["aa_vo"]("a,i");
    
        // rt2_abab += -1.00 d-_aa(a,i) u1_bb(b,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= dp["aa_vo"]("a,i") * u1["bb_vo"]("b,j");
    
        // rt2_abab += -1.00 d-_aa(k,k) u2_abab(a,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= scalars_["12"] * u2["abab_vvoo"]("a,b,i,j");
    
        // rt2_abab += -1.00 d-_bb(k,k) u2_abab(a,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= scalars_["11"] * u2["abab_vvoo"]("a,b,i,j");
    
        // rt2_abab += +1.00 <a,b||i,j>_abab  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += V_blks_["abab_vvoo"]("a,b,i,j");
    
        // rt2_abab += -1.00 <a,k||c,j>_abab t2_abab(c,b,i,k)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= V_blks_["abab_vovo"]("a,k,c,j") * t2["abab_vvoo"]("c,b,i,k");
    
        // rt2_abab += -1.00 <k,b||c,j>_abab t2_aaaa(c,a,i,k)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += V_blks_["baab_vovo"]("b,k,c,j") * t2["aaaa_vvoo"]("c,a,i,k");
    
        // rt2_abab += +1.00 <k,b||c,j>_bbbb t2_abab(a,c,i,k)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= V_blks_["bbbb_vovo"]("b,k,c,j") * t2["abab_vvoo"]("a,c,i,k");
    
        // rt2_abab += +1.00 <k,a||c,i>_aaaa t2_abab(c,b,k,j)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= V_blks_["aaaa_vovo"]("a,k,c,i") * t2["abab_vvoo"]("c,b,k,j");
    
        // rt2_abab += -1.00 <a,k||i,c>_abab t2_bbbb(c,b,j,k)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += V_blks_["abba_vovo"]("a,k,c,i") * t2["bbbb_vvoo"]("c,b,j,k");
    
        // rt2_abab += -1.00 <k,b||i,c>_abab t2_abab(a,c,k,j)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= V_blks_["baba_vovo"]("b,k,c,i") * t2["abab_vvoo"]("a,c,k,j");
    
        // rt2_bbbb = +1.00 P(i,j) P(a,b) d-_bb(a,j) u1_bb(b,i)  // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
        rt2_bbbb("a,b,i,j")  = dp["bb_vo"]("a,j") * u1["bb_vo"]("b,i");
        rt2_bbbb("a,b,i,j") -= dp["bb_vo"]("a,i") * u1["bb_vo"]("b,j");
        rt2_bbbb("a,b,i,j") -= dp["bb_vo"]("b,j") * u1["bb_vo"]("a,i");
        rt2_bbbb("a,b,i,j") += dp["bb_vo"]("b,i") * u1["bb_vo"]("a,j");
    
        // rt2_bbbb += -1.00 d-_aa(k,k) u2_bbbb(a,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") -= scalars_["12"] * u2["bbbb_vvoo"]("a,b,i,j");
    
        // rt2_bbbb += -1.00 d-_bb(k,k) u2_bbbb(a,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") -= scalars_["11"] * u2["bbbb_vvoo"]("a,b,i,j");
    
        // rt2_bbbb += +1.00 <a,b||i,j>_bbbb  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") += V_blks_["bbbb_vvoo"]("a,b,i,j");
    
        // rt2_bbbb += +0.50 <l,k||i,j>_bbbb t2_bbbb(a,b,l,k)  // flops: o2v2 += o4v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") += 0.50 * V_blks_["bbbb_oooo"]("l,k,i,j") * t2["bbbb_vvoo"]("a,b,l,k");
    
        // rt2_bbbb += +0.50 <a,b||c,d>_bbbb t2_bbbb(c,d,i,j)  // flops: o2v2 += o2v4 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") += 0.50 * V_blks_["bbbb_vvvv"]("a,b,c,d") * t2["bbbb_vvoo"]("c,d,i,j");
    
        // ru0 = +1.00 u0 w0  // flops: o0v0 = o0v0 | mem: o0v0 = o0v0
        ru0  = u0 * w0;
    
        // ru0 += +1.00 f_aa(i,a) u1_aa(a,i)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        ru0 += scalars_["10"];
    
        // ru0 += +1.00 f_bb(i,a) u1_bb(a,i)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        ru0 += scalars_["9"];
    
        // ru0 += -1.00 d+_aa(i,i)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        ru0 -= scalars_["12"];
    
        // ru0 += -1.00 d+_bb(i,i)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        ru0 -= scalars_["11"];
    
        // ru0 += -1.00 d-_aa(i,a) u0 u1_aa(a,i)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        ru0 -= scalars_["7"] * u0;
    
        // ru0 += -1.00 d-_bb(i,a) u0 u1_bb(a,i)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        ru0 -= scalars_["8"] * u0;
    
        // ru0 += +0.25 <j,i||a,b>_aaaa u2_aaaa(a,b,j,i)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        ru0 += 0.25 * scalars_["5"];
    
        // ru0 += +0.25 <j,i||a,b>_abab u2_abab(a,b,j,i)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        ru0 += scalars_["2"];
    
        // ru0 += +0.25 <i,j||a,b>_abab u2_abab(a,b,i,j)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        ru0 += scalars_["2"];
    
        // ru0 += +0.25 <j,i||b,a>_abab u2_abab(b,a,j,i)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        ru0 += scalars_["2"];
    
        // ru0 += +0.25 <i,j||b,a>_abab u2_abab(b,a,i,j)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        ru0 += scalars_["2"];
    
        // ru0 += +0.25 <j,i||a,b>_bbbb u2_bbbb(a,b,j,i)  // flops: o0v0 += o0v0 | mem: o0v0 += o0v0
        ru0 += 0.25 * scalars_["6"];
    
        // ru1_aa = +1.00 u1_aa(a,i) w0  // flops: o1v1 = o1v1 | mem: o1v1 = o1v1
        ru1_aa("a,i")  = u1["aa_vo"]("a,i") * w0;
    
        // ru1_aa += -1.00 f_aa(j,i) u1_aa(a,j)  // flops: o1v1 += o2v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= F_blks_["aa_oo"]("j,i") * u1["aa_vo"]("a,j");
    
        // ru1_aa += +1.00 f_aa(a,b) u1_aa(b,i)  // flops: o1v1 += o1v2 | mem: o1v1 += o1v1
        ru1_aa("a,i") += F_blks_["aa_vv"]("a,b") * u1["aa_vo"]("b,i");
    
        // ru1_aa += -1.00 f_aa(j,b) u2_aaaa(b,a,i,j)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= F_blks_["aa_ov"]("j,b") * u2["aaaa_vvoo"]("b,a,i,j");
    
        // ru1_aa += +1.00 f_bb(j,b) u2_abab(a,b,i,j)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        ru1_aa("a,i") += F_blks_["bb_ov"]("j,b") * u2["abab_vvoo"]("a,b,i,j");
    
        // ru1_aa += -1.00 d+_aa(a,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= dp["aa_vo"]("a,i");
    
        // ru1_aa += -1.00 d-_aa(j,b) u1_aa(b,j) u1_aa(a,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= scalars_["7"] * u1["aa_vo"]("a,i");
    
        // ru1_aa += -1.00 d-_bb(j,b) u1_bb(b,j) u1_aa(a,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= scalars_["8"] * u1["aa_vo"]("a,i");
    
        // ru1_aa += +1.00 <j,a||b,i>_aaaa u1_aa(b,j)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= V_blks_["aaaa_vovo"]("a,j,b,i") * u1["aa_vo"]("b,j");
    
        // ru1_aa += +1.00 <a,j||i,b>_abab u1_bb(b,j)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= V_blks_["abba_vovo"]("a,j,b,i") * u1["bb_vo"]("b,j");
    
        // ru1_aa += -0.50 <k,j||b,i>_aaaa u2_aaaa(b,a,k,j)  // flops: o1v1 += o3v2 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= 0.50 * V_blks_["aaaa_oovo"]("k,j,b,i") * u2["aaaa_vvoo"]("b,a,k,j");
    
        // ru1_aa += -0.50 <j,a||b,c>_aaaa u2_aaaa(b,c,i,j)  // flops: o1v1 += o2v3 | mem: o1v1 += o1v1
        ru1_aa("a,i") += 0.50 * V_blks_["aaaa_vovv"]("a,j,b,c") * u2["aaaa_vvoo"]("b,c,i,j");
    
        // ru1_bb = +1.00 u1_bb(a,i) w0  // flops: o1v1 = o1v1 | mem: o1v1 = o1v1
        ru1_bb("a,i")  = u1["bb_vo"]("a,i") * w0;
    
        // ru1_bb += -1.00 f_bb(j,i) u1_bb(a,j)  // flops: o1v1 += o2v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= F_blks_["bb_oo"]("j,i") * u1["bb_vo"]("a,j");
    
        // ru1_bb += +1.00 f_bb(a,b) u1_bb(b,i)  // flops: o1v1 += o1v2 | mem: o1v1 += o1v1
        ru1_bb("a,i") += F_blks_["bb_vv"]("a,b") * u1["bb_vo"]("b,i");
    
        // ru1_bb += +1.00 f_aa(j,b) u2_abab(b,a,j,i)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        ru1_bb("a,i") += F_blks_["aa_ov"]("j,b") * u2["abab_vvoo"]("b,a,j,i");
    
        // ru1_bb += -1.00 f_bb(j,b) u2_bbbb(b,a,i,j)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= F_blks_["bb_ov"]("j,b") * u2["bbbb_vvoo"]("b,a,i,j");
    
        // ru1_bb += -1.00 d+_bb(a,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= dp["bb_vo"]("a,i");
    
        // ru1_bb += -1.00 d-_aa(j,b) u1_aa(b,j) u1_bb(a,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= scalars_["7"] * u1["bb_vo"]("a,i");
    
        // ru1_bb += -1.00 d-_bb(j,b) u1_bb(b,j) u1_bb(a,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= scalars_["8"] * u1["bb_vo"]("a,i");
    
        // ru1_bb += +1.00 <j,a||b,i>_abab u1_aa(b,j)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= V_blks_["baab_vovo"]("a,j,b,i") * u1["aa_vo"]("b,j");
    
        // ru1_bb += +1.00 <j,a||b,i>_bbbb u1_bb(b,j)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= V_blks_["bbbb_vovo"]("a,j,b,i") * u1["bb_vo"]("b,j");
    
        // ru1_bb += -0.50 <k,j||b,i>_bbbb u2_bbbb(b,a,k,j)  // flops: o1v1 += o3v2 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= 0.50 * V_blks_["bbbb_oovo"]("k,j,b,i") * u2["bbbb_vvoo"]("b,a,k,j");
    
        // ru1_bb += -0.50 <j,a||b,c>_bbbb u2_bbbb(b,c,i,j)  // flops: o1v1 += o2v3 | mem: o1v1 += o1v1
        ru1_bb("a,i") += 0.50 * V_blks_["bbbb_vovv"]("a,j,b,c") * u2["bbbb_vvoo"]("b,c,i,j");
    
        // ru2_aaaa = +1.00 u2_aaaa(a,b,i,j) w0  // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
        ru2_aaaa("a,b,i,j")  = u2["aaaa_vvoo"]("a,b,i,j") * w0;
    
        // ru2_aaaa += -1.00 d-_aa(k,c) u1_aa(c,k) u2_aaaa(a,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= scalars_["7"] * u2["aaaa_vvoo"]("a,b,i,j");
    
        // ru2_aaaa += -1.00 d-_bb(k,c) u1_bb(c,k) u2_aaaa(a,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= scalars_["8"] * u2["aaaa_vvoo"]("a,b,i,j");
    
        // ru2_aaaa += +0.50 <l,k||i,j>_aaaa u2_aaaa(a,b,l,k)  // flops: o2v2 += o4v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += 0.50 * V_blks_["aaaa_oooo"]("l,k,i,j") * u2["aaaa_vvoo"]("a,b,l,k");
    
        // ru2_aaaa += +0.50 <a,b||c,d>_aaaa u2_aaaa(c,d,i,j)  // flops: o2v2 += o2v4 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += 0.50 * V_blks_["aaaa_vvvv"]("a,b,c,d") * u2["aaaa_vvoo"]("c,d,i,j");
    
        // ru2_aaaa += +0.25 <l,k||c,d>_aaaa u2_aaaa(c,d,i,j) t2_aaaa(a,b,l,k)  // flops: o2v2 += o4v2 o4v2 | mem: o2v2 += o4v0 o2v2
        ru2_aaaa("a,b,i,j") += 0.25 * V_blks_["aaaa_oovv"]("l,k,c,d") * u2["aaaa_vvoo"]("c,d,i,j") * t2["aaaa_vvoo"]("a,b,l,k");
    
        // ru2_abab = +1.00 u2_abab(a,b,i,j) w0  // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
        ru2_abab("a,b,i,j")  = u2["abab_vvoo"]("a,b,i,j") * w0;
    
        // ru2_abab += -1.00 f_bb(k,j) u2_abab(a,b,i,k)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= F_blks_["bb_oo"]("k,j") * u2["abab_vvoo"]("a,b,i,k");
    
        // ru2_abab += -1.00 f_aa(k,i) u2_abab(a,b,k,j)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= F_blks_["aa_oo"]("k,i") * u2["abab_vvoo"]("a,b,k,j");
    
        // ru2_abab += +1.00 f_aa(a,c) u2_abab(c,b,i,j)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += F_blks_["aa_vv"]("a,c") * u2["abab_vvoo"]("c,b,i,j");
    
        // ru2_abab += +1.00 f_bb(b,c) u2_abab(a,c,i,j)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += F_blks_["bb_vv"]("b,c") * u2["abab_vvoo"]("a,c,i,j");
    
        // ru2_abab += -1.00 f_bb(k,c) u1_bb(b,k) t2_abab(a,c,i,j)  // flops: o2v2 += o3v2 o3v2 | mem: o2v2 += o3v1 o2v2
        ru2_abab("a,b,i,j") -= F_blks_["bb_ov"]("k,c") * t2["abab_vvoo"]("a,c,i,j") * u1["bb_vo"]("b,k");
    
        // ru2_abab += -1.00 f_aa(k,c) u1_aa(a,k) t2_abab(c,b,i,j)  // flops: o2v2 += o3v2 o3v2 | mem: o2v2 += o3v1 o2v2
        ru2_abab("a,b,i,j") -= F_blks_["aa_ov"]("k,c") * t2["abab_vvoo"]("c,b,i,j") * u1["aa_vo"]("a,k");
    
        // ru2_abab += -1.00 d-_aa(k,c) u1_aa(c,k) u2_abab(a,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= scalars_["7"] * u2["abab_vvoo"]("a,b,i,j");
    
        // ru2_abab += -1.00 d-_bb(k,c) u1_bb(c,k) u2_abab(a,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= scalars_["8"] * u2["abab_vvoo"]("a,b,i,j");
    
        // ru2_abab += +2.00 d-_aa(k,c) u1_aa(a,k) u2_abab(c,b,i,j)  // flops: o2v2 += o3v2 o3v2 | mem: o2v2 += o3v1 o2v2
        ru2_abab("a,b,i,j") += 2.00 * dp["aa_ov"]("k,c") * u2["abab_vvoo"]("c,b,i,j") * u1["aa_vo"]("a,k");
    
        // ru2_abab += +2.00 d-_bb(k,c) u1_bb(b,k) u2_abab(a,c,i,j)  // flops: o2v2 += o3v2 o3v2 | mem: o2v2 += o3v1 o2v2
        ru2_abab("a,b,i,j") += 2.00 * dp["bb_ov"]("k,c") * u2["abab_vvoo"]("a,c,i,j") * u1["bb_vo"]("b,k");
    
        // ru2_abab += -1.00 <a,k||i,j>_abab u1_bb(b,k)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= V_blks_["abab_vooo"]("a,k,i,j") * u1["bb_vo"]("b,k");
    
        // ru2_abab += -1.00 <k,b||i,j>_abab u1_aa(a,k)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += V_blks_["baab_vooo"]("b,k,i,j") * u1["aa_vo"]("a,k");
    
        // ru2_abab += +1.00 <a,b||c,j>_abab u1_aa(c,i)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += V_blks_["abab_vvvo"]("a,b,c,j") * u1["aa_vo"]("c,i");
    
        // ru2_abab += +1.00 <a,b||i,c>_abab u1_bb(c,j)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= V_blks_["abba_vvvo"]("a,b,c,i") * u1["bb_vo"]("c,j");
    
        // ru2_abab += -1.00 <a,k||c,j>_abab u2_abab(c,b,i,k)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= V_blks_["abab_vovo"]("a,k,c,j") * u2["abab_vvoo"]("c,b,i,k");
    
        // ru2_abab += -1.00 <k,b||c,j>_abab u2_aaaa(c,a,i,k)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += V_blks_["baab_vovo"]("b,k,c,j") * u2["aaaa_vvoo"]("c,a,i,k");
    
        // ru2_abab += +1.00 <k,b||c,j>_bbbb u2_abab(a,c,i,k)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= V_blks_["bbbb_vovo"]("b,k,c,j") * u2["abab_vvoo"]("a,c,i,k");
    
        // ru2_abab += +1.00 <k,a||c,i>_aaaa u2_abab(c,b,k,j)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= V_blks_["aaaa_vovo"]("a,k,c,i") * u2["abab_vvoo"]("c,b,k,j");
    
        // ru2_abab += -1.00 <a,k||i,c>_abab u2_bbbb(c,b,j,k)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += V_blks_["abba_vovo"]("a,k,c,i") * u2["bbbb_vvoo"]("c,b,j,k");
    
        // ru2_abab += -1.00 <k,b||i,c>_abab u2_abab(a,c,k,j)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= V_blks_["baba_vovo"]("b,k,c,i") * u2["abab_vvoo"]("a,c,k,j");
    
        // ru2_abab += +1.00 <k,l||c,j>_abab u1_bb(b,l) t2_aaaa(c,a,i,k)  // flops: o2v2 += o4v2 o3v2 | mem: o2v2 += o3v1 o2v2
        ru2_abab("a,b,i,j") += V_blks_["abab_oovo"]("k,l,c,j") * t2["aaaa_vvoo"]("c,a,i,k") * u1["bb_vo"]("b,l");
    
        // ru2_abab += +1.00 <l,k||c,j>_bbbb u1_bb(b,l) t2_abab(a,c,i,k)  // flops: o2v2 += o4v2 o3v2 | mem: o2v2 += o3v1 o2v2
        ru2_abab("a,b,i,j") += V_blks_["bbbb_oovo"]("l,k,c,j") * t2["abab_vvoo"]("a,c,i,k") * u1["bb_vo"]("b,l");
    
        // ru2_abab += +1.00 <l,k||c,j>_abab u1_aa(a,l) t2_abab(c,b,i,k)  // flops: o2v2 += o4v2 o3v2 | mem: o2v2 += o3v1 o2v2
        ru2_abab("a,b,i,j") += V_blks_["abab_oovo"]("l,k,c,j") * t2["abab_vvoo"]("c,b,i,k") * u1["aa_vo"]("a,l");
    
        // ru2_abab += +1.00 <k,l||i,c>_abab u1_bb(b,l) t2_abab(a,c,k,j)  // flops: o2v2 += o4v2 o3v2 | mem: o2v2 += o3v1 o2v2
        ru2_abab("a,b,i,j") -= V_blks_["abba_oovo"]("k,l,c,i") * t2["abab_vvoo"]("a,c,k,j") * u1["bb_vo"]("b,l");
    
        // ru2_abab += +1.00 <l,k||c,i>_aaaa u1_aa(a,l) t2_abab(c,b,k,j)  // flops: o2v2 += o4v2 o3v2 | mem: o2v2 += o3v1 o2v2
        ru2_abab("a,b,i,j") += V_blks_["aaaa_oovo"]("l,k,c,i") * t2["abab_vvoo"]("c,b,k,j") * u1["aa_vo"]("a,l");
    
        // ru2_abab += +1.00 <l,k||i,c>_abab u1_aa(a,l) t2_bbbb(c,b,j,k)  // flops: o2v2 += o4v2 o3v2 | mem: o2v2 += o3v1 o2v2
        ru2_abab("a,b,i,j") -= V_blks_["abba_oovo"]("l,k,c,i") * t2["bbbb_vvoo"]("c,b,j,k") * u1["aa_vo"]("a,l");
    
        // ru2_abab += -1.00 <k,b||d,c>_abab u1_aa(d,i) t2_abab(a,c,k,j)  // flops: o2v2 += o2v3 o3v3 | mem: o2v2 += o2v2 o2v2
        ru2_abab("a,b,i,j") += V_blks_["baab_vovv"]("b,k,d,c") * u1["aa_vo"]("d,i") * t2["abab_vvoo"]("a,c,k,j");
    
        // ru2_abab += -1.00 <a,k||c,d>_abab u1_bb(d,j) t2_abab(c,b,i,k)  // flops: o2v2 += o2v3 o3v3 | mem: o2v2 += o2v2 o2v2
        ru2_abab("a,b,i,j") -= V_blks_["abab_vovv"]("a,k,c,d") * u1["bb_vo"]("d,j") * t2["abab_vvoo"]("c,b,i,k");
    
        // ru2_abab += +1.00 <k,l||d,c>_abab u2_abab(d,b,i,l) t2_abab(a,c,k,j)  // flops: o2v2 += o3v3 o3v3 | mem: o2v2 += o2v2 o2v2
        ru2_abab("a,b,i,j") += V_blks_["abab_oovv"]("k,l,d,c") * t2["abab_vvoo"]("a,c,k,j") * u2["abab_vvoo"]("d,b,i,l");
    
        // ru2_bbbb = +1.00 u2_bbbb(a,b,i,j) w0  // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
        ru2_bbbb("a,b,i,j")  = u2["bbbb_vvoo"]("a,b,i,j") * w0;
    
        // ru2_bbbb += -1.00 d-_aa(k,c) u1_aa(c,k) u2_bbbb(a,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= scalars_["7"] * u2["bbbb_vvoo"]("a,b,i,j");
    
        // ru2_bbbb += -1.00 d-_bb(k,c) u1_bb(c,k) u2_bbbb(a,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= scalars_["8"] * u2["bbbb_vvoo"]("a,b,i,j");
    
        // ru2_bbbb += +0.50 <l,k||i,j>_bbbb u2_bbbb(a,b,l,k)  // flops: o2v2 += o4v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += 0.50 * V_blks_["bbbb_oooo"]("l,k,i,j") * u2["bbbb_vvoo"]("a,b,l,k");
    
        // ru2_bbbb += +0.50 <a,b||c,d>_bbbb u2_bbbb(c,d,i,j)  // flops: o2v2 += o2v4 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += 0.50 * V_blks_["bbbb_vvvv"]("a,b,c,d") * u2["bbbb_vvoo"]("c,d,i,j");
    
        // ru2_bbbb += +0.25 <l,k||c,d>_bbbb u2_bbbb(c,d,i,j) t2_bbbb(a,b,l,k)  // flops: o2v2 += o4v2 o4v2 | mem: o2v2 += o4v0 o2v2
        ru2_bbbb("a,b,i,j") += 0.25 * V_blks_["bbbb_oovv"]("l,k,c,d") * u2["bbbb_vvoo"]("c,d,i,j") * t2["bbbb_vvoo"]("a,b,l,k");
    
        // tmps_[abab_vvoo_1](a,b,i,j) = 0.50 V_blks_[abab_vvvv](a,b,d,c) * u2[abab_vvoo](d,c,i,j) // flops: o2v2 = o2v4 | mem: o2v2 = o2v2
        tmps_["abab_vvoo_1"]("a,b,i,j")  = 0.50 * V_blks_["abab_vvvv"]("a,b,d,c") * u2["abab_vvoo"]("d,c,i,j");
    
        // ru2_abab += +0.50 <a,b||c,d>_abab u2_abab(c,d,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["abab_vvoo_1"]("a,b,i,j");
    
        // ru2_abab += +0.50 <a,b||d,c>_abab u2_abab(d,c,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["abab_vvoo_1"]("a,b,i,j");
        tmps_["abab_vvoo_1"].~TArrayD();
    
        // tmps_[abab_vvoo_2](a,b,i,j) = 0.50 V_blks_[abab_vvvv](a,b,c,d) * t2[abab_vvoo](c,d,i,j) // flops: o2v2 = o2v4 | mem: o2v2 = o2v2
        tmps_["abab_vvoo_2"]("a,b,i,j")  = 0.50 * V_blks_["abab_vvvv"]("a,b,c,d") * t2["abab_vvoo"]("c,d,i,j");
    
        // rt2_abab += +0.50 <a,b||c,d>_abab t2_abab(c,d,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["abab_vvoo_2"]("a,b,i,j");
    
        // rt2_abab += +0.50 <a,b||d,c>_abab t2_abab(d,c,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["abab_vvoo_2"]("a,b,i,j");
        tmps_["abab_vvoo_2"].~TArrayD();
    
        // tmps_[aabb_ovvo_3](l,d,a,j) = 1.00 V_blks_[abab_oovv](l,k,d,c) * t2[bbbb_vvoo](c,a,j,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["aabb_ovvo_3"]("l,d,a,j")  = V_blks_["abab_oovv"]("l,k,d,c") * t2["bbbb_vvoo"]("c,a,j,k");
    
        // rt2_abab += +1.00 <k,l||c,d>_abab t2_aaaa(c,a,i,k) t2_bbbb(d,b,j,l)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["aabb_ovvo_3"]("k,c,b,j") * t2["aaaa_vvoo"]("c,a,i,k");
    
        // ru2_abab += +1.00 <l,k||d,c>_abab u2_aaaa(d,a,i,l) t2_bbbb(c,b,j,k)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["aabb_ovvo_3"]("l,d,b,j") * u2["aaaa_vvoo"]("d,a,i,l");
    
        // tmps_[bbbb_vovo_91](b,j,a,i) = 1.00 u2[abab_vvoo](d,b,l,j) * V_blks_[abab_oovv](l,k,d,c) * t2[bbbb_vvoo](c,a,i,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vovo_91"]("b,j,a,i")  = u2["abab_vvoo"]("d,b,l,j") * tmps_["aabb_ovvo_3"]("l,d,a,i");
    
        // ru2_bbbb += +1.00 P(i,j) P(a,b) <l,k||d,c>_abab u2_abab(d,b,l,i) t2_bbbb(c,a,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_91"]("b,i,a,j");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_91"]("b,j,a,i");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_91"]("a,i,b,j");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_91"]("a,j,b,i");
        tmps_["bbbb_vovo_91"].~TArrayD();
    
        // tmps_[bbbb_vovo_103](a,j,b,i) = 1.00 t2[abab_vvoo](c,a,k,j) * V_blks_[abab_oovv](k,l,c,d) * t2[bbbb_vvoo](d,b,i,l) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vovo_103"]("a,j,b,i")  = t2["abab_vvoo"]("c,a,k,j") * tmps_["aabb_ovvo_3"]("k,c,b,i");
        tmps_["aabb_ovvo_3"].~TArrayD();
    
        // rt2_bbbb += +1.00 P(i,j) <k,l||c,d>_abab t2_abab(c,a,k,j) t2_bbbb(d,b,i,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_103"]("a,j,b,i");
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_103"]("a,i,b,j");
    
        // rt2_bbbb += +1.00 P(i,j) <l,k||d,c>_abab t2_bbbb(c,a,j,k) t2_abab(d,b,l,i)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_103"]("b,i,a,j");
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_103"]("b,j,a,i");
        tmps_["bbbb_vovo_103"].~TArrayD();
    
        // tmps_[bbbb_ovvo_4](k,c,b,i) = 1.00 V_blks_[abab_oovv](l,k,d,c) * t2[abab_vvoo](d,b,l,i) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["bbbb_ovvo_4"]("k,c,b,i")  = V_blks_["abab_oovv"]("l,k,d,c") * t2["abab_vvoo"]("d,b,l,i");
    
        // rt2_abab += +1.00 <l,k||d,c>_abab t2_abab(a,c,i,k) t2_abab(d,b,l,j)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["bbbb_ovvo_4"]("k,c,b,j") * t2["abab_vvoo"]("a,c,i,k");
    
        // ru2_abab += +1.00 <k,l||c,d>_abab u2_abab(a,d,i,l) t2_abab(c,b,k,j)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["bbbb_ovvo_4"]("l,d,b,j") * u2["abab_vvoo"]("a,d,i,l");
    
        // tmps_[bbbb_vovo_92](b,i,a,j) = 1.00 u2[bbbb_vvoo](d,b,i,l) * V_blks_[abab_oovv](k,l,c,d) * t2[abab_vvoo](c,a,k,j) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vovo_92"]("b,i,a,j")  = u2["bbbb_vvoo"]("d,b,i,l") * tmps_["bbbb_ovvo_4"]("l,d,a,j");
        tmps_["bbbb_ovvo_4"].~TArrayD();
    
        // ru2_bbbb += +1.00 P(i,j) P(a,b) <k,l||c,d>_abab u2_bbbb(d,b,i,l) t2_abab(c,a,k,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_92"]("b,i,a,j");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_92"]("b,j,a,i");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_92"]("a,i,b,j");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_92"]("a,j,b,i");
        tmps_["bbbb_vovo_92"].~TArrayD();
    
        // tmps_[aaaa_ovvo_5](k,c,b,i) = 1.00 V_blks_[abab_oovv](k,l,c,d) * t2[abab_vvoo](b,d,i,l) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["aaaa_ovvo_5"]("k,c,b,i")  = V_blks_["abab_oovv"]("k,l,c,d") * t2["abab_vvoo"]("b,d,i,l");
    
        // ru2_abab += +1.00 <l,k||d,c>_abab u2_abab(d,b,l,j) t2_abab(a,c,i,k)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["aaaa_ovvo_5"]("l,d,a,i") * u2["abab_vvoo"]("d,b,l,j");
    
        // tmps_[aaaa_vovo_97](b,j,a,i) = 1.00 u2[aaaa_vvoo](d,b,j,l) * V_blks_[abab_oovv](l,k,d,c) * t2[abab_vvoo](a,c,i,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vovo_97"]("b,j,a,i")  = u2["aaaa_vvoo"]("d,b,j,l") * tmps_["aaaa_ovvo_5"]("l,d,a,i");
    
        // ru2_aaaa += +1.00 P(i,j) P(a,b) <l,k||d,c>_abab u2_aaaa(d,b,i,l) t2_abab(a,c,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_97"]("b,i,a,j");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_97"]("b,j,a,i");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_97"]("a,i,b,j");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_97"]("a,j,b,i");
        tmps_["aaaa_vovo_97"].~TArrayD();
    
        // tmps_[aaaa_vovo_104](a,j,b,i) = 1.00 t2[aaaa_vvoo](c,a,j,k) * V_blks_[abab_oovv](k,l,c,d) * t2[abab_vvoo](b,d,i,l) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vovo_104"]("a,j,b,i")  = t2["aaaa_vvoo"]("c,a,j,k") * tmps_["aaaa_ovvo_5"]("k,c,b,i");
        tmps_["aaaa_ovvo_5"].~TArrayD();
    
        // rt2_aaaa += +1.00 P(i,j) <k,l||c,d>_abab t2_aaaa(c,a,j,k) t2_abab(b,d,i,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_104"]("a,j,b,i");
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_104"]("a,i,b,j");
    
        // rt2_aaaa += +1.00 P(i,j) <l,k||d,c>_abab t2_abab(a,c,j,k) t2_aaaa(d,b,i,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_104"]("b,i,a,j");
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_104"]("b,j,a,i");
        tmps_["aaaa_vovo_104"].~TArrayD();
    
        // tmps_[bbaa_ovvo_6](l,d,a,i) = 1.00 V_blks_[abab_oovv](k,l,c,d) * t2[aaaa_vvoo](c,a,i,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["bbaa_ovvo_6"]("l,d,a,i")  = V_blks_["abab_oovv"]("k,l,c,d") * t2["aaaa_vvoo"]("c,a,i,k");
    
        // ru2_abab += +1.00 <k,l||c,d>_abab u2_bbbb(d,b,j,l) t2_aaaa(c,a,i,k)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["bbaa_ovvo_6"]("l,d,a,i") * u2["bbbb_vvoo"]("d,b,j,l");
    
        // tmps_[aaaa_vovo_98](a,j,b,i) = 1.00 u2[abab_vvoo](a,d,j,l) * V_blks_[abab_oovv](k,l,c,d) * t2[aaaa_vvoo](c,b,i,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vovo_98"]("a,j,b,i")  = u2["abab_vvoo"]("a,d,j,l") * tmps_["bbaa_ovvo_6"]("l,d,b,i");
        tmps_["bbaa_ovvo_6"].~TArrayD();
    
        // ru2_aaaa += +1.00 P(i,j) P(a,b) <k,l||c,d>_abab u2_abab(b,d,i,l) t2_aaaa(c,a,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_98"]("b,i,a,j");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_98"]("b,j,a,i");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_98"]("a,i,b,j");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_98"]("a,j,b,i");
        tmps_["aaaa_vovo_98"].~TArrayD();
    
        // tmps_[aaaa_ovvo_7](l,d,a,i) = 1.00 V_blks_[aaaa_oovv](l,k,c,d) * t2[aaaa_vvoo](c,a,i,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["aaaa_ovvo_7"]("l,d,a,i")  = V_blks_["aaaa_oovv"]("l,k,c,d") * t2["aaaa_vvoo"]("c,a,i,k");
    
        // rt2_abab += +1.00 <l,k||c,d>_aaaa t2_aaaa(c,a,i,k) t2_abab(d,b,l,j)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["aaaa_ovvo_7"]("l,d,a,i") * t2["abab_vvoo"]("d,b,l,j");
    
        // ru2_abab += +1.00 <l,k||c,d>_aaaa u2_abab(d,b,l,j) t2_aaaa(c,a,i,k)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["aaaa_ovvo_7"]("l,d,a,i") * u2["abab_vvoo"]("d,b,l,j");
    
        // tmps_[aaaa_vovo_100](b,i,a,j) = 1.00 V_blks_[aaaa_oovv](l,k,c,d) * t2[aaaa_vvoo](c,b,i,k) * u2[aaaa_vvoo](d,a,j,l) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vovo_100"]("b,i,a,j")  = tmps_["aaaa_ovvo_7"]("l,d,b,i") * u2["aaaa_vvoo"]("d,a,j,l");
    
        // ru2_aaaa += +1.00 P(i,j) P(a,b) <l,k||c,d>_aaaa u2_aaaa(d,b,i,l) t2_aaaa(c,a,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_100"]("a,j,b,i");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_100"]("a,i,b,j");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_100"]("b,j,a,i");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_100"]("b,i,a,j");
        tmps_["aaaa_vovo_100"].~TArrayD();
    
        // tmps_[aaaa_vovo_107](b,i,a,j) = 1.00 t2[aaaa_vvoo](d,b,i,l) * V_blks_[aaaa_oovv](l,k,c,d) * t2[aaaa_vvoo](c,a,j,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vovo_107"]("b,i,a,j")  = t2["aaaa_vvoo"]("d,b,i,l") * tmps_["aaaa_ovvo_7"]("l,d,a,j");
        tmps_["aaaa_ovvo_7"].~TArrayD();
    
        // rt2_aaaa += +1.00 P(i,j) <l,k||c,d>_aaaa t2_aaaa(c,a,j,k) t2_aaaa(d,b,i,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_107"]("b,i,a,j");
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_107"]("b,j,a,i");
        tmps_["aaaa_vovo_107"].~TArrayD();
    
        // tmps_[bbaa_ovvo_8](l,d,a,i) = 1.00 V_blks_[bbbb_oovv](l,k,c,d) * t2[abab_vvoo](a,c,i,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["bbaa_ovvo_8"]("l,d,a,i")  = V_blks_["bbbb_oovv"]("l,k,c,d") * t2["abab_vvoo"]("a,c,i,k");
    
        // rt2_abab += +1.00 <l,k||c,d>_bbbb t2_abab(a,c,i,k) t2_bbbb(d,b,j,l)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["bbaa_ovvo_8"]("l,d,a,i") * t2["bbbb_vvoo"]("d,b,j,l");
    
        // ru2_abab += +1.00 <l,k||c,d>_bbbb u2_bbbb(d,b,j,l) t2_abab(a,c,i,k)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["bbaa_ovvo_8"]("l,d,a,i") * u2["bbbb_vvoo"]("d,b,j,l");
    
        // tmps_[aaaa_vovo_99](b,j,a,i) = 1.00 u2[abab_vvoo](b,d,j,l) * V_blks_[bbbb_oovv](l,k,c,d) * t2[abab_vvoo](a,c,i,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vovo_99"]("b,j,a,i")  = u2["abab_vvoo"]("b,d,j,l") * tmps_["bbaa_ovvo_8"]("l,d,a,i");
    
        // ru2_aaaa += +1.00 P(i,j) P(a,b) <l,k||c,d>_bbbb u2_abab(b,d,i,l) t2_abab(a,c,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_99"]("b,i,a,j");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_99"]("b,j,a,i");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_99"]("a,i,b,j");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_99"]("a,j,b,i");
        tmps_["aaaa_vovo_99"].~TArrayD();
    
        // tmps_[aaaa_vovo_108](b,i,a,j) = 1.00 t2[abab_vvoo](b,d,i,l) * V_blks_[bbbb_oovv](l,k,c,d) * t2[abab_vvoo](a,c,j,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vovo_108"]("b,i,a,j")  = t2["abab_vvoo"]("b,d,i,l") * tmps_["bbaa_ovvo_8"]("l,d,a,j");
        tmps_["bbaa_ovvo_8"].~TArrayD();
    
        // rt2_aaaa += +1.00 P(i,j) <l,k||c,d>_bbbb t2_abab(a,c,j,k) t2_abab(b,d,i,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_108"]("b,i,a,j");
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_108"]("b,j,a,i");
        tmps_["aaaa_vovo_108"].~TArrayD();
    
        // tmps_[aabb_ovvo_9](l,d,b,i) = 1.00 V_blks_[aaaa_oovv](l,k,c,d) * t2[abab_vvoo](c,b,k,i) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["aabb_ovvo_9"]("l,d,b,i")  = V_blks_["aaaa_oovv"]("l,k,c,d") * t2["abab_vvoo"]("c,b,k,i");
    
        // ru2_abab += +1.00 <l,k||c,d>_aaaa u2_aaaa(d,a,i,l) t2_abab(c,b,k,j)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["aabb_ovvo_9"]("l,d,b,j") * u2["aaaa_vvoo"]("d,a,i,l");
    
        // tmps_[bbbb_vovo_94](b,j,a,i) = 1.00 V_blks_[aaaa_oovv](l,k,c,d) * t2[abab_vvoo](c,b,k,j) * u2[abab_vvoo](d,a,l,i) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vovo_94"]("b,j,a,i")  = tmps_["aabb_ovvo_9"]("l,d,b,j") * u2["abab_vvoo"]("d,a,l,i");
    
        // ru2_bbbb += +1.00 P(i,j) P(a,b) <l,k||c,d>_aaaa u2_abab(d,b,l,i) t2_abab(c,a,k,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_94"]("a,j,b,i");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_94"]("a,i,b,j");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_94"]("b,j,a,i");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_94"]("b,i,a,j");
        tmps_["bbbb_vovo_94"].~TArrayD();
    
        // tmps_[bbbb_vovo_105](b,i,a,j) = 1.00 t2[abab_vvoo](d,b,l,i) * V_blks_[aaaa_oovv](l,k,c,d) * t2[abab_vvoo](c,a,k,j) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vovo_105"]("b,i,a,j")  = t2["abab_vvoo"]("d,b,l,i") * tmps_["aabb_ovvo_9"]("l,d,a,j");
        tmps_["aabb_ovvo_9"].~TArrayD();
    
        // rt2_bbbb += +1.00 P(i,j) <l,k||c,d>_aaaa t2_abab(c,a,k,j) t2_abab(d,b,l,i)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_105"]("b,i,a,j");
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_105"]("b,j,a,i");
        tmps_["bbbb_vovo_105"].~TArrayD();
    
        // tmps_[bbbb_ovvo_10](l,d,a,i) = 1.00 V_blks_[bbbb_oovv](l,k,c,d) * t2[bbbb_vvoo](c,a,i,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["bbbb_ovvo_10"]("l,d,a,i")  = V_blks_["bbbb_oovv"]("l,k,c,d") * t2["bbbb_vvoo"]("c,a,i,k");
    
        // ru2_abab += +1.00 <l,k||c,d>_bbbb u2_abab(a,d,i,l) t2_bbbb(c,b,j,k)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["bbbb_ovvo_10"]("l,d,b,j") * u2["abab_vvoo"]("a,d,i,l");
    
        // tmps_[bbbb_vovo_93](a,i,b,j) = 1.00 V_blks_[bbbb_oovv](l,k,c,d) * t2[bbbb_vvoo](c,a,i,k) * u2[bbbb_vvoo](d,b,j,l) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vovo_93"]("a,i,b,j")  = tmps_["bbbb_ovvo_10"]("l,d,a,i") * u2["bbbb_vvoo"]("d,b,j,l");
    
        // ru2_bbbb += +1.00 P(i,j) P(a,b) <l,k||c,d>_bbbb u2_bbbb(d,b,i,l) t2_bbbb(c,a,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_93"]("a,j,b,i");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_93"]("a,i,b,j");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_93"]("b,j,a,i");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_93"]("b,i,a,j");
        tmps_["bbbb_vovo_93"].~TArrayD();
    
        // tmps_[bbbb_vovo_106](a,j,b,i) = 1.00 V_blks_[bbbb_oovv](l,k,c,d) * t2[bbbb_vvoo](c,a,j,k) * t2[bbbb_vvoo](d,b,i,l) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vovo_106"]("a,j,b,i")  = tmps_["bbbb_ovvo_10"]("l,d,a,j") * t2["bbbb_vvoo"]("d,b,i,l");
        tmps_["bbbb_ovvo_10"].~TArrayD();
    
        // rt2_bbbb += +1.00 P(i,j) <l,k||c,d>_bbbb t2_bbbb(c,a,j,k) t2_bbbb(d,b,i,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_106"]("a,j,b,i");
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_106"]("a,i,b,j");
        tmps_["bbbb_vovo_106"].~TArrayD();
    
        // tmps_[bbbb_vovo_11](b,j,a,i) = 1.00 V_blks_[baab_vovo](b,k,c,j) * u2[abab_vvoo](c,a,k,i) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vovo_11"]("b,j,a,i")  = V_blks_["baab_vovo"]("b,k,c,j") * u2["abab_vvoo"]("c,a,k,i");
    
        // ru2_bbbb += -1.00 P(i,j) P(a,b) <k,a||c,j>_abab u2_abab(c,b,k,i)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_11"]("a,j,b,i");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_11"]("a,i,b,j");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_11"]("b,j,a,i");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_11"]("b,i,a,j");
        tmps_["bbbb_vovo_11"].~TArrayD();
    
        // tmps_[bbbb_vovo_12](a,j,b,i) = 1.00 V_blks_[bbbb_vovo](a,k,c,j) * u2[bbbb_vvoo](c,b,i,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vovo_12"]("a,j,b,i")  = V_blks_["bbbb_vovo"]("a,k,c,j") * u2["bbbb_vvoo"]("c,b,i,k");
    
        // ru2_bbbb += +1.00 P(i,j) P(a,b) <k,a||c,j>_bbbb u2_bbbb(c,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_12"]("a,j,b,i");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_12"]("a,i,b,j");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_12"]("b,j,a,i");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_12"]("b,i,a,j");
        tmps_["bbbb_vovo_12"].~TArrayD();
    
        // tmps_[bbbb_vovo_13](b,j,a,i) = 1.00 V_blks_[bbbb_vovo](b,k,c,j) * t2[bbbb_vvoo](c,a,i,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vovo_13"]("b,j,a,i")  = V_blks_["bbbb_vovo"]("b,k,c,j") * t2["bbbb_vvoo"]("c,a,i,k");
    
        // rt2_bbbb += +1.00 P(i,j) P(a,b) <k,a||c,j>_bbbb t2_bbbb(c,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_13"]("a,j,b,i");
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_13"]("a,i,b,j");
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_13"]("b,j,a,i");
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_13"]("b,i,a,j");
        tmps_["bbbb_vovo_13"].~TArrayD();
    
        // tmps_[aaaa_vovo_14](b,i,a,j) = 1.00 V_blks_[aaaa_vovo](b,k,c,i) * t2[aaaa_vvoo](c,a,j,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vovo_14"]("b,i,a,j")  = V_blks_["aaaa_vovo"]("b,k,c,i") * t2["aaaa_vvoo"]("c,a,j,k");
    
        // rt2_aaaa += +1.00 P(i,j) P(a,b) <k,a||c,j>_aaaa t2_aaaa(c,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_14"]("a,j,b,i");
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_14"]("a,i,b,j");
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_14"]("b,j,a,i");
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_14"]("b,i,a,j");
        tmps_["aaaa_vovo_14"].~TArrayD();
    
        // tmps_[aaaa_vovo_15](a,j,b,i) = 1.00 V_blks_[abba_vovo](a,k,c,j) * t2[abab_vvoo](b,c,i,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vovo_15"]("a,j,b,i")  = V_blks_["abba_vovo"]("a,k,c,j") * t2["abab_vvoo"]("b,c,i,k");
    
        // rt2_aaaa += -1.00 P(i,j) P(a,b) <a,k||j,c>_abab t2_abab(b,c,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_15"]("a,j,b,i");
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_15"]("a,i,b,j");
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_15"]("b,j,a,i");
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_15"]("b,i,a,j");
        tmps_["aaaa_vovo_15"].~TArrayD();
    
        // tmps_[aaaa_vovo_16](a,j,b,i) = 1.00 V_blks_[abba_vovo](a,k,c,j) * u2[abab_vvoo](b,c,i,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vovo_16"]("a,j,b,i")  = V_blks_["abba_vovo"]("a,k,c,j") * u2["abab_vvoo"]("b,c,i,k");
    
        // ru2_aaaa += -1.00 P(i,j) P(a,b) <a,k||j,c>_abab u2_abab(b,c,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_16"]("a,j,b,i");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_16"]("a,i,b,j");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_16"]("b,j,a,i");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_16"]("b,i,a,j");
        tmps_["aaaa_vovo_16"].~TArrayD();
    
        // tmps_[bbbb_vovo_17](b,j,a,i) = 1.00 V_blks_[baab_vovo](b,k,c,j) * t2[abab_vvoo](c,a,k,i) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vovo_17"]("b,j,a,i")  = V_blks_["baab_vovo"]("b,k,c,j") * t2["abab_vvoo"]("c,a,k,i");
    
        // rt2_bbbb += -1.00 P(i,j) P(a,b) <k,a||c,j>_abab t2_abab(c,b,k,i)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_17"]("a,j,b,i");
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_17"]("a,i,b,j");
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_17"]("b,j,a,i");
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_17"]("b,i,a,j");
        tmps_["bbbb_vovo_17"].~TArrayD();
    
        // tmps_[aaaa_vovo_18](a,i,b,j) = 1.00 V_blks_[aaaa_vovo](a,k,c,i) * u2[aaaa_vvoo](c,b,j,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vovo_18"]("a,i,b,j")  = V_blks_["aaaa_vovo"]("a,k,c,i") * u2["aaaa_vvoo"]("c,b,j,k");
    
        // ru2_aaaa += +1.00 P(i,j) P(a,b) <k,a||c,j>_aaaa u2_aaaa(c,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_18"]("a,j,b,i");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_18"]("a,i,b,j");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_18"]("b,j,a,i");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_18"]("b,i,a,j");
        tmps_["aaaa_vovo_18"].~TArrayD();
    
        // tmps_[abba_ovvo_23](l,d,b,i) = 1.00 V_blks_[abab_oovv](l,k,c,d) * t2[abab_vvoo](c,b,i,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["abba_ovvo_23"]("l,d,b,i")  = V_blks_["abab_oovv"]("l,k,c,d") * t2["abab_vvoo"]("c,b,i,k");
    
        // rt2_abab += +1.00 <k,l||d,c>_abab t2_abab(a,c,k,j) t2_abab(d,b,i,l)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["abba_ovvo_23"]("k,c,b,i") * t2["abab_vvoo"]("a,c,k,j");
    
        // ru2_abab += +1.00 <l,k||c,d>_abab u2_abab(a,d,l,j) t2_abab(c,b,i,k)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["abba_ovvo_23"]("l,d,b,i") * u2["abab_vvoo"]("a,d,l,j");
        tmps_["abba_ovvo_23"].~TArrayD();
    
        // tmps_[abab_oooo_25](k,l,i,j) = 0.25 V_blks_[abab_oovv](k,l,c,d) * t2[abab_vvoo](c,d,i,j) // flops: o4v0 = o4v2 | mem: o4v0 = o4v0
        tmps_["abab_oooo_25"]("k,l,i,j")  = 0.25 * V_blks_["abab_oovv"]("k,l,c,d") * t2["abab_vvoo"]("c,d,i,j");
    
        // tmps_[abab_vvoo_110](a,b,i,j) = 1.00 u2[abab_vvoo](a,b,l,k) * V_blks_[abab_oovv](l,k,c,d) * t2[abab_vvoo](c,d,i,j) // flops: o2v2 = o4v2 | mem: o2v2 = o2v2
        tmps_["abab_vvoo_110"]("a,b,i,j")  = u2["abab_vvoo"]("a,b,l,k") * tmps_["abab_oooo_25"]("l,k,i,j");
    
        // ru2_abab += +0.25 <l,k||c,d>_abab u2_abab(a,b,l,k) t2_abab(c,d,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["abab_vvoo_110"]("a,b,i,j");
    
        // ru2_abab += +0.25 <k,l||c,d>_abab u2_abab(a,b,k,l) t2_abab(c,d,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["abab_vvoo_110"]("a,b,i,j");
    
        // ru2_abab += +0.25 <l,k||d,c>_abab u2_abab(a,b,l,k) t2_abab(d,c,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["abab_vvoo_110"]("a,b,i,j");
    
        // ru2_abab += +0.25 <k,l||d,c>_abab u2_abab(a,b,k,l) t2_abab(d,c,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["abab_vvoo_110"]("a,b,i,j");
        tmps_["abab_vvoo_110"].~TArrayD();
    
        // tmps_[abab_vvoo_111](a,b,i,j) = 1.00 t2[abab_vvoo](a,b,l,k) * V_blks_[abab_oovv](l,k,c,d) * t2[abab_vvoo](c,d,i,j) // flops: o2v2 = o4v2 | mem: o2v2 = o2v2
        tmps_["abab_vvoo_111"]("a,b,i,j")  = t2["abab_vvoo"]("a,b,l,k") * tmps_["abab_oooo_25"]("l,k,i,j");
        tmps_["abab_oooo_25"].~TArrayD();
    
        // rt2_abab += +0.25 <l,k||c,d>_abab t2_abab(c,d,i,j) t2_abab(a,b,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["abab_vvoo_111"]("a,b,i,j");
    
        // rt2_abab += +0.25 <k,l||c,d>_abab t2_abab(c,d,i,j) t2_abab(a,b,k,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["abab_vvoo_111"]("a,b,i,j");
    
        // rt2_abab += +0.25 <l,k||d,c>_abab t2_abab(d,c,i,j) t2_abab(a,b,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["abab_vvoo_111"]("a,b,i,j");
    
        // rt2_abab += +0.25 <k,l||d,c>_abab t2_abab(d,c,i,j) t2_abab(a,b,k,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["abab_vvoo_111"]("a,b,i,j");
        tmps_["abab_vvoo_111"].~TArrayD();
    
        // tmps_[aaaa_oooo_31](l,k,i,j) = 0.25 V_blks_[aaaa_oovv](l,k,c,d) * t2[aaaa_vvoo](c,d,i,j) // flops: o4v0 = o4v2 | mem: o4v0 = o4v0
        tmps_["aaaa_oooo_31"]("l,k,i,j")  = 0.25 * V_blks_["aaaa_oovv"]("l,k,c,d") * t2["aaaa_vvoo"]("c,d,i,j");
    
        // rt2_aaaa += +0.25 <l,k||c,d>_aaaa t2_aaaa(c,d,i,j) t2_aaaa(a,b,l,k)  // flops: o2v2 += o4v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_oooo_31"]("l,k,i,j") * t2["aaaa_vvoo"]("a,b,l,k");
    
        // ru2_aaaa += +0.25 <l,k||c,d>_aaaa u2_aaaa(a,b,l,k) t2_aaaa(c,d,i,j)  // flops: o2v2 += o4v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_oooo_31"]("l,k,i,j") * u2["aaaa_vvoo"]("a,b,l,k");
        tmps_["aaaa_oooo_31"].~TArrayD();
    
        // tmps_[bbbb_oooo_32](l,k,i,j) = 0.25 V_blks_[bbbb_oovv](l,k,c,d) * t2[bbbb_vvoo](c,d,i,j) // flops: o4v0 = o4v2 | mem: o4v0 = o4v0
        tmps_["bbbb_oooo_32"]("l,k,i,j")  = 0.25 * V_blks_["bbbb_oovv"]("l,k,c,d") * t2["bbbb_vvoo"]("c,d,i,j");
    
        // rt2_bbbb += +0.25 <l,k||c,d>_bbbb t2_bbbb(c,d,i,j) t2_bbbb(a,b,l,k)  // flops: o2v2 += o4v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_oooo_32"]("l,k,i,j") * t2["bbbb_vvoo"]("a,b,l,k");
    
        // ru2_bbbb += +0.25 <l,k||c,d>_bbbb u2_bbbb(a,b,l,k) t2_bbbb(c,d,i,j)  // flops: o2v2 += o4v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_oooo_32"]("l,k,i,j") * u2["bbbb_vvoo"]("a,b,l,k");
        tmps_["bbbb_oooo_32"].~TArrayD();
    
        // tmps_[abab_oovv_33](i,j,a,b) = 0.50 V_blks_[abab_oooo](k,l,i,j) * u2[abab_vvoo](a,b,k,l) // flops: o2v2 = o4v2 | mem: o2v2 = o2v2
        tmps_["abab_oovv_33"]("i,j,a,b")  = 0.50 * V_blks_["abab_oooo"]("k,l,i,j") * u2["abab_vvoo"]("a,b,k,l");
    
        // ru2_abab += +0.50 <l,k||i,j>_abab u2_abab(a,b,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["abab_oovv_33"]("i,j,a,b");
    
        // ru2_abab += +0.50 <k,l||i,j>_abab u2_abab(a,b,k,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["abab_oovv_33"]("i,j,a,b");
        tmps_["abab_oovv_33"].~TArrayD();
    
        // tmps_[abab_oovv_34](i,j,a,b) = 0.50 V_blks_[abab_oooo](k,l,i,j) * t2[abab_vvoo](a,b,k,l) // flops: o2v2 = o4v2 | mem: o2v2 = o2v2
        tmps_["abab_oovv_34"]("i,j,a,b")  = 0.50 * V_blks_["abab_oooo"]("k,l,i,j") * t2["abab_vvoo"]("a,b,k,l");
    
        // rt2_abab += +0.50 <l,k||i,j>_abab t2_abab(a,b,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["abab_oovv_34"]("i,j,a,b");
    
        // rt2_abab += +0.50 <k,l||i,j>_abab t2_abab(a,b,k,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["abab_oovv_34"]("i,j,a,b");
        tmps_["abab_oovv_34"].~TArrayD();
    
        // tmps_[bb_vv_35](c,b) = 0.50 V_blks_[abab_oovv](l,k,d,c) * t2[abab_vvoo](d,b,l,k) // flops: o0v2 = o2v3 | mem: o0v2 = o0v2
        tmps_["bb_vv_35"]("c,b")  = 0.50 * V_blks_["abab_oovv"]("l,k,d,c") * t2["abab_vvoo"]("d,b,l,k");
    
        // tmps_[bbbb_voov_112](a,i,j,b) = 1.00 u2[bbbb_vvoo](d,a,i,j) * V_blks_[abab_oovv](l,k,c,d) * t2[abab_vvoo](c,b,l,k) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["bbbb_voov_112"]("a,i,j,b")  = u2["bbbb_vvoo"]("d,a,i,j") * tmps_["bb_vv_35"]("d,b");
    
        // ru2_bbbb += -0.50 P(a,b) <l,k||c,d>_abab u2_bbbb(d,b,i,j) t2_abab(c,a,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_voov_112"]("b,i,j,a");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_voov_112"]("a,i,j,b");
    
        // ru2_bbbb += -0.50 P(a,b) <k,l||c,d>_abab u2_bbbb(d,b,i,j) t2_abab(c,a,k,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_voov_112"]("b,i,j,a");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_voov_112"]("a,i,j,b");
        tmps_["bbbb_voov_112"].~TArrayD();
    
        // tmps_[bbbb_voov_116](b,i,j,a) = 1.00 t2[bbbb_vvoo](d,b,i,j) * V_blks_[abab_oovv](l,k,c,d) * t2[abab_vvoo](c,a,l,k) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["bbbb_voov_116"]("b,i,j,a")  = t2["bbbb_vvoo"]("d,b,i,j") * tmps_["bb_vv_35"]("d,a");
    
        // rt2_bbbb += -0.50 <l,k||c,d>_abab t2_abab(c,a,l,k) t2_bbbb(d,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_voov_116"]("b,i,j,a");
    
        // rt2_bbbb += -0.50 <k,l||c,d>_abab t2_abab(c,a,k,l) t2_bbbb(d,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_voov_116"]("b,i,j,a");
    
        // rt2_bbbb += +0.50 <l,k||d,c>_abab t2_bbbb(c,a,i,j) t2_abab(d,b,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_voov_116"]("a,i,j,b");
    
        // rt2_bbbb += +0.50 <k,l||d,c>_abab t2_bbbb(c,a,i,j) t2_abab(d,b,k,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_voov_116"]("a,i,j,b");
        tmps_["bbbb_voov_116"].~TArrayD();
    
        // tmps_[aabb_voov_123](a,i,j,b) = 1.00 u2[abab_vvoo](a,d,i,j) * V_blks_[abab_oovv](k,l,c,d) * t2[abab_vvoo](c,b,k,l) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["aabb_voov_123"]("a,i,j,b")  = u2["abab_vvoo"]("a,d,i,j") * tmps_["bb_vv_35"]("d,b");
    
        // ru2_abab += -0.50 <l,k||c,d>_abab u2_abab(a,d,i,j) t2_abab(c,b,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["aabb_voov_123"]("a,i,j,b");
    
        // ru2_abab += -0.50 <k,l||c,d>_abab u2_abab(a,d,i,j) t2_abab(c,b,k,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["aabb_voov_123"]("a,i,j,b");
        tmps_["aabb_voov_123"].~TArrayD();
    
        // tmps_[aabb_voov_127](a,i,j,b) = 1.00 t2[abab_vvoo](a,c,i,j) * V_blks_[abab_oovv](l,k,d,c) * t2[abab_vvoo](d,b,l,k) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["aabb_voov_127"]("a,i,j,b")  = t2["abab_vvoo"]("a,c,i,j") * tmps_["bb_vv_35"]("c,b");
    
        // rt2_abab += -0.50 <l,k||d,c>_abab t2_abab(a,c,i,j) t2_abab(d,b,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= tmps_["aabb_voov_127"]("a,i,j,b");
    
        // rt2_abab += -0.50 <k,l||d,c>_abab t2_abab(a,c,i,j) t2_abab(d,b,k,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= tmps_["aabb_voov_127"]("a,i,j,b");
        tmps_["aabb_voov_127"].~TArrayD();
    
        // tmps_[bb_ov_206](i,a) = 1.00 u1[bb_vo](c,i) * V_blks_[abab_oovv](k,j,b,c) * t2[abab_vvoo](b,a,k,j) // flops: o1v1 = o1v2 | mem: o1v1 = o1v1
        tmps_["bb_ov_206"]("i,a")  = u1["bb_vo"]("c,i") * tmps_["bb_vv_35"]("c,a");
        tmps_["bb_vv_35"].~TArrayD();
    
        // ru1_bb += -0.50 <k,j||b,c>_abab u1_bb(c,i) t2_abab(b,a,k,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= tmps_["bb_ov_206"]("i,a");
    
        // ru1_bb += -0.50 <j,k||b,c>_abab u1_bb(c,i) t2_abab(b,a,j,k)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= tmps_["bb_ov_206"]("i,a");
        tmps_["bb_ov_206"].~TArrayD();
    
        // tmps_[aa_vv_36](d,a) = 0.50 V_blks_[abab_oovv](k,l,d,c) * t2[abab_vvoo](a,c,k,l) // flops: o0v2 = o2v3 | mem: o0v2 = o0v2
        tmps_["aa_vv_36"]("d,a")  = 0.50 * V_blks_["abab_oovv"]("k,l,d,c") * t2["abab_vvoo"]("a,c,k,l");
    
        // tmps_[aaaa_vvoo_114](b,a,i,j) = 1.00 V_blks_[abab_oovv](l,k,d,c) * t2[abab_vvoo](b,c,l,k) * u2[aaaa_vvoo](d,a,i,j) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vvoo_114"]("b,a,i,j")  = tmps_["aa_vv_36"]("d,b") * u2["aaaa_vvoo"]("d,a,i,j");
    
        // ru2_aaaa += -0.50 P(a,b) <l,k||d,c>_abab u2_aaaa(d,b,i,j) t2_abab(a,c,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_114"]("a,b,i,j");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_114"]("b,a,i,j");
    
        // ru2_aaaa += -0.50 P(a,b) <k,l||d,c>_abab u2_aaaa(d,b,i,j) t2_abab(a,c,k,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_114"]("a,b,i,j");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_114"]("b,a,i,j");
        tmps_["aaaa_vvoo_114"].~TArrayD();
    
        // tmps_[aaaa_voov_117](b,i,j,a) = 1.00 t2[aaaa_vvoo](d,b,i,j) * V_blks_[abab_oovv](l,k,d,c) * t2[abab_vvoo](a,c,l,k) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["aaaa_voov_117"]("b,i,j,a")  = t2["aaaa_vvoo"]("d,b,i,j") * tmps_["aa_vv_36"]("d,a");
    
        // rt2_aaaa += -0.50 <l,k||d,c>_abab t2_abab(a,c,l,k) t2_aaaa(d,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_voov_117"]("b,i,j,a");
    
        // rt2_aaaa += -0.50 <k,l||d,c>_abab t2_abab(a,c,k,l) t2_aaaa(d,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_voov_117"]("b,i,j,a");
    
        // rt2_aaaa += +0.50 <l,k||c,d>_abab t2_aaaa(c,a,i,j) t2_abab(b,d,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_voov_117"]("a,i,j,b");
    
        // rt2_aaaa += +0.50 <k,l||c,d>_abab t2_aaaa(c,a,i,j) t2_abab(b,d,k,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_voov_117"]("a,i,j,b");
        tmps_["aaaa_voov_117"].~TArrayD();
    
        // tmps_[baba_voov_120](b,i,j,a) = 1.00 u2[abab_vvoo](d,b,i,j) * V_blks_[abab_oovv](k,l,d,c) * t2[abab_vvoo](a,c,k,l) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["baba_voov_120"]("b,i,j,a")  = u2["abab_vvoo"]("d,b,i,j") * tmps_["aa_vv_36"]("d,a");
    
        // ru2_abab += -0.50 <l,k||d,c>_abab u2_abab(d,b,i,j) t2_abab(a,c,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["baba_voov_120"]("b,i,j,a");
    
        // ru2_abab += -0.50 <k,l||d,c>_abab u2_abab(d,b,i,j) t2_abab(a,c,k,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["baba_voov_120"]("b,i,j,a");
        tmps_["baba_voov_120"].~TArrayD();
    
        // tmps_[baba_voov_126](b,i,j,a) = 1.00 t2[abab_vvoo](d,b,i,j) * V_blks_[abab_oovv](k,l,d,c) * t2[abab_vvoo](a,c,k,l) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["baba_voov_126"]("b,i,j,a")  = t2["abab_vvoo"]("d,b,i,j") * tmps_["aa_vv_36"]("d,a");
    
        // rt2_abab += -0.50 <l,k||d,c>_abab t2_abab(a,c,l,k) t2_abab(d,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= tmps_["baba_voov_126"]("b,i,j,a");
    
        // rt2_abab += -0.50 <k,l||d,c>_abab t2_abab(a,c,k,l) t2_abab(d,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= tmps_["baba_voov_126"]("b,i,j,a");
        tmps_["baba_voov_126"].~TArrayD();
    
        // tmps_[aa_ov_207](i,a) = 1.00 u1[aa_vo](c,i) * V_blks_[abab_oovv](k,j,c,b) * t2[abab_vvoo](a,b,k,j) // flops: o1v1 = o1v2 | mem: o1v1 = o1v1
        tmps_["aa_ov_207"]("i,a")  = u1["aa_vo"]("c,i") * tmps_["aa_vv_36"]("c,a");
        tmps_["aa_vv_36"].~TArrayD();
    
        // ru1_aa += -0.50 <k,j||c,b>_abab u1_aa(c,i) t2_abab(a,b,k,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= tmps_["aa_ov_207"]("i,a");
    
        // ru1_aa += -0.50 <j,k||c,b>_abab u1_aa(c,i) t2_abab(a,b,j,k)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= tmps_["aa_ov_207"]("i,a");
        tmps_["aa_ov_207"].~TArrayD();
    
        // tmps_[bb_vv_37](c,a) = 0.50 V_blks_[abab_oovv](k,l,d,c) * u2[abab_vvoo](d,a,k,l) // flops: o0v2 = o2v3 | mem: o0v2 = o0v2
        tmps_["bb_vv_37"]("c,a")  = 0.50 * V_blks_["abab_oovv"]("k,l,d,c") * u2["abab_vvoo"]("d,a,k,l");
    
        // tmps_[bbbb_voov_113](b,i,j,a) = 1.00 t2[bbbb_vvoo](c,b,i,j) * V_blks_[abab_oovv](l,k,d,c) * u2[abab_vvoo](d,a,l,k) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["bbbb_voov_113"]("b,i,j,a")  = t2["bbbb_vvoo"]("c,b,i,j") * tmps_["bb_vv_37"]("c,a");
    
        // ru2_bbbb += +0.50 P(a,b) <l,k||d,c>_abab u2_abab(d,b,l,k) t2_bbbb(c,a,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_voov_113"]("a,i,j,b");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_voov_113"]("b,i,j,a");
    
        // ru2_bbbb += +0.50 P(a,b) <k,l||d,c>_abab u2_abab(d,b,k,l) t2_bbbb(c,a,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_voov_113"]("a,i,j,b");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_voov_113"]("b,i,j,a");
        tmps_["bbbb_voov_113"].~TArrayD();
    
        // tmps_[aabb_voov_121](a,i,j,b) = 1.00 t2[abab_vvoo](a,c,i,j) * V_blks_[abab_oovv](l,k,d,c) * u2[abab_vvoo](d,b,l,k) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["aabb_voov_121"]("a,i,j,b")  = t2["abab_vvoo"]("a,c,i,j") * tmps_["bb_vv_37"]("c,b");
        tmps_["bb_vv_37"].~TArrayD();
    
        // ru2_abab += -0.50 <l,k||d,c>_abab u2_abab(d,b,l,k) t2_abab(a,c,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["aabb_voov_121"]("a,i,j,b");
    
        // ru2_abab += -0.50 <k,l||d,c>_abab u2_abab(d,b,k,l) t2_abab(a,c,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["aabb_voov_121"]("a,i,j,b");
        tmps_["aabb_voov_121"].~TArrayD();
    
        // tmps_[aa_vv_38](c,a) = 0.50 V_blks_[abab_oovv](l,k,c,d) * u2[abab_vvoo](a,d,l,k) // flops: o0v2 = o2v3 | mem: o0v2 = o0v2
        tmps_["aa_vv_38"]("c,a")  = 0.50 * V_blks_["abab_oovv"]("l,k,c,d") * u2["abab_vvoo"]("a,d,l,k");
    
        // tmps_[aaaa_voov_115](a,i,j,b) = 1.00 t2[aaaa_vvoo](c,a,i,j) * V_blks_[abab_oovv](l,k,c,d) * u2[abab_vvoo](b,d,l,k) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["aaaa_voov_115"]("a,i,j,b")  = t2["aaaa_vvoo"]("c,a,i,j") * tmps_["aa_vv_38"]("c,b");
    
        // ru2_aaaa += +0.50 P(a,b) <l,k||c,d>_abab u2_abab(b,d,l,k) t2_aaaa(c,a,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_voov_115"]("a,i,j,b");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_voov_115"]("b,i,j,a");
    
        // ru2_aaaa += +0.50 P(a,b) <k,l||c,d>_abab u2_abab(b,d,k,l) t2_aaaa(c,a,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_voov_115"]("a,i,j,b");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_voov_115"]("b,i,j,a");
        tmps_["aaaa_voov_115"].~TArrayD();
    
        // tmps_[abab_vvoo_122](a,b,i,j) = 1.00 V_blks_[abab_oovv](l,k,c,d) * u2[abab_vvoo](a,d,l,k) * t2[abab_vvoo](c,b,i,j) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["abab_vvoo_122"]("a,b,i,j")  = tmps_["aa_vv_38"]("c,a") * t2["abab_vvoo"]("c,b,i,j");
        tmps_["aa_vv_38"].~TArrayD();
    
        // ru2_abab += -0.50 <l,k||c,d>_abab u2_abab(a,d,l,k) t2_abab(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["abab_vvoo_122"]("a,b,i,j");
    
        // ru2_abab += -0.50 <k,l||c,d>_abab u2_abab(a,d,k,l) t2_abab(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["abab_vvoo_122"]("a,b,i,j");
        tmps_["abab_vvoo_122"].~TArrayD();
    
        // tmps_[aa_vv_39](d,a) = 0.50 V_blks_[aaaa_oovv](l,k,c,d) * t2[aaaa_vvoo](c,a,l,k) // flops: o0v2 = o2v3 | mem: o0v2 = o0v2
        tmps_["aa_vv_39"]("d,a")  = 0.50 * V_blks_["aaaa_oovv"]("l,k,c,d") * t2["aaaa_vvoo"]("c,a,l,k");
    
        // rt2_aaaa += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,a,l,k) t2_aaaa(d,b,i,j)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") -= tmps_["aa_vv_39"]("d,a") * t2["aaaa_vvoo"]("d,b,i,j");
    
        // rt2_abab += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,a,l,k) t2_abab(d,b,i,j)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= tmps_["aa_vv_39"]("d,a") * t2["abab_vvoo"]("d,b,i,j");
    
        // ru1_aa += -0.50 <k,j||b,c>_aaaa u1_aa(c,i) t2_aaaa(b,a,k,j)  // flops: o1v1 += o1v2 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= tmps_["aa_vv_39"]("c,a") * u1["aa_vo"]("c,i");
    
        // ru2_abab += -0.50 <l,k||c,d>_aaaa u2_abab(d,b,i,j) t2_aaaa(c,a,l,k)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["aa_vv_39"]("d,a") * u2["abab_vvoo"]("d,b,i,j");
    
        // tmps_[aaaa_voov_124](a,i,j,b) = 1.00 u2[aaaa_vvoo](d,a,i,j) * V_blks_[aaaa_oovv](l,k,c,d) * t2[aaaa_vvoo](c,b,l,k) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["aaaa_voov_124"]("a,i,j,b")  = u2["aaaa_vvoo"]("d,a,i,j") * tmps_["aa_vv_39"]("d,b");
        tmps_["aa_vv_39"].~TArrayD();
    
        // ru2_aaaa += -0.50 P(a,b) <l,k||c,d>_aaaa u2_aaaa(d,b,i,j) t2_aaaa(c,a,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_voov_124"]("b,i,j,a");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_voov_124"]("a,i,j,b");
        tmps_["aaaa_voov_124"].~TArrayD();
    
        // tmps_[bbbb_vovo_40](b,k,c,i) = 1.00 V_blks_[bbbb_vovv](b,k,c,d) * u1[bb_vo](d,i) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vovo_40"]("b,k,c,i")  = V_blks_["bbbb_vovv"]("b,k,c,d") * u1["bb_vo"]("d,i");
    
        // ru2_abab += +1.00 <k,b||c,d>_bbbb u1_bb(d,j) t2_abab(a,c,i,k)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["bbbb_vovo_40"]("b,k,c,j") * t2["abab_vvoo"]("a,c,i,k");
    
        // tmps_[bbbb_vovo_96](a,j,b,i) = 1.00 t2[bbbb_vvoo](c,a,j,k) * V_blks_[bbbb_vovv](b,k,c,d) * u1[bb_vo](d,i) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vovo_96"]("a,j,b,i")  = t2["bbbb_vvoo"]("c,a,j,k") * tmps_["bbbb_vovo_40"]("b,k,c,i");
        tmps_["bbbb_vovo_40"].~TArrayD();
    
        // ru2_bbbb += -1.00 P(i,j) P(a,b) <k,a||c,d>_bbbb u1_bb(d,i) t2_bbbb(c,b,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_96"]("b,j,a,i");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_96"]("b,i,a,j");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_96"]("a,j,b,i");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_96"]("a,i,b,j");
        tmps_["bbbb_vovo_96"].~TArrayD();
    
        // tmps_[baab_vovo_41](b,k,c,j) = 1.00 V_blks_[baab_vovv](b,k,c,d) * u1[bb_vo](d,j) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["baab_vovo_41"]("b,k,c,j")  = V_blks_["baab_vovv"]("b,k,c,d") * u1["bb_vo"]("d,j");
    
        // ru2_abab += -1.00 <k,b||c,d>_abab u1_bb(d,j) t2_aaaa(c,a,i,k)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["baab_vovo_41"]("b,k,c,j") * t2["aaaa_vvoo"]("c,a,i,k");
    
        // tmps_[bbbb_vovo_95](b,j,a,i) = 1.00 t2[abab_vvoo](c,b,k,j) * V_blks_[baab_vovv](a,k,c,d) * u1[bb_vo](d,i) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vovo_95"]("b,j,a,i")  = t2["abab_vvoo"]("c,b,k,j") * tmps_["baab_vovo_41"]("a,k,c,i");
        tmps_["baab_vovo_41"].~TArrayD();
    
        // ru2_bbbb += +1.00 P(i,j) P(a,b) <k,a||c,d>_abab u1_bb(d,i) t2_abab(c,b,k,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_95"]("b,j,a,i");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_95"]("b,i,a,j");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vovo_95"]("a,j,b,i");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vovo_95"]("a,i,b,j");
        tmps_["bbbb_vovo_95"].~TArrayD();
    
        // tmps_[aaaa_vovo_42](a,k,c,i) = 1.00 V_blks_[aaaa_vovv](a,k,c,d) * u1[aa_vo](d,i) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vovo_42"]("a,k,c,i")  = V_blks_["aaaa_vovv"]("a,k,c,d") * u1["aa_vo"]("d,i");
    
        // ru2_abab += +1.00 <k,a||c,d>_aaaa u1_aa(d,i) t2_abab(c,b,k,j)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["aaaa_vovo_42"]("a,k,c,i") * t2["abab_vvoo"]("c,b,k,j");
    
        // tmps_[aaaa_vovo_102](a,i,b,j) = 1.00 t2[aaaa_vvoo](c,a,i,k) * V_blks_[aaaa_vovv](b,k,c,d) * u1[aa_vo](d,j) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vovo_102"]("a,i,b,j")  = t2["aaaa_vvoo"]("c,a,i,k") * tmps_["aaaa_vovo_42"]("b,k,c,j");
        tmps_["aaaa_vovo_42"].~TArrayD();
    
        // ru2_aaaa += -1.00 P(i,j) P(a,b) <k,a||c,d>_aaaa u1_aa(d,i) t2_aaaa(c,b,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_102"]("b,j,a,i");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_102"]("b,i,a,j");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_102"]("a,j,b,i");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_102"]("a,i,b,j");
        tmps_["aaaa_vovo_102"].~TArrayD();
    
        // tmps_[abba_vovo_43](a,k,c,j) = 1.00 V_blks_[abab_vovv](a,k,d,c) * u1[aa_vo](d,j) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["abba_vovo_43"]("a,k,c,j")  = V_blks_["abab_vovv"]("a,k,d,c") * u1["aa_vo"]("d,j");
    
        // ru2_abab += -1.00 <a,k||d,c>_abab u1_aa(d,i) t2_bbbb(c,b,j,k)  // flops: o2v2 += o3v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["abba_vovo_43"]("a,k,c,i") * t2["bbbb_vvoo"]("c,b,j,k");
    
        // tmps_[aaaa_vovo_101](b,j,a,i) = 1.00 t2[abab_vvoo](b,c,j,k) * V_blks_[abab_vovv](a,k,d,c) * u1[aa_vo](d,i) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vovo_101"]("b,j,a,i")  = t2["abab_vvoo"]("b,c,j,k") * tmps_["abba_vovo_43"]("a,k,c,i");
        tmps_["abba_vovo_43"].~TArrayD();
    
        // ru2_aaaa += +1.00 P(i,j) P(a,b) <a,k||d,c>_abab u1_aa(d,i) t2_abab(b,c,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_101"]("b,j,a,i");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_101"]("b,i,a,j");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vovo_101"]("a,j,b,i");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vovo_101"]("a,i,b,j");
        tmps_["aaaa_vovo_101"].~TArrayD();
    
        // tmps_[bb_vv_44](d,a) = 0.50 V_blks_[bbbb_oovv](l,k,c,d) * t2[bbbb_vvoo](c,a,l,k) // flops: o0v2 = o2v3 | mem: o0v2 = o0v2
        tmps_["bb_vv_44"]("d,a")  = 0.50 * V_blks_["bbbb_oovv"]("l,k,c,d") * t2["bbbb_vvoo"]("c,a,l,k");
    
        // rt2_bbbb += -0.50 <l,k||c,d>_bbbb t2_bbbb(c,a,l,k) t2_bbbb(d,b,i,j)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") -= tmps_["bb_vv_44"]("d,a") * t2["bbbb_vvoo"]("d,b,i,j");
    
        // ru1_bb += -0.50 <k,j||b,c>_bbbb u1_bb(c,i) t2_bbbb(b,a,k,j)  // flops: o1v1 += o1v2 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= tmps_["bb_vv_44"]("c,a") * u1["bb_vo"]("c,i");
    
        // ru2_abab += -0.50 <l,k||c,d>_bbbb u2_abab(a,d,i,j) t2_bbbb(c,b,l,k)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["bb_vv_44"]("d,b") * u2["abab_vvoo"]("a,d,i,j");
    
        // tmps_[bbbb_vvoo_118](b,a,i,j) = 1.00 V_blks_[bbbb_oovv](l,k,c,d) * t2[bbbb_vvoo](c,b,l,k) * u2[bbbb_vvoo](d,a,i,j) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vvoo_118"]("b,a,i,j")  = tmps_["bb_vv_44"]("d,b") * u2["bbbb_vvoo"]("d,a,i,j");
        tmps_["bb_vv_44"].~TArrayD();
    
        // ru2_bbbb += -0.50 P(a,b) <l,k||c,d>_bbbb u2_bbbb(d,b,i,j) t2_bbbb(c,a,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_118"]("a,b,i,j");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_118"]("b,a,i,j");
        tmps_["bbbb_vvoo_118"].~TArrayD();
    
        // tmps_[aaaa_vvoo_45](b,a,i,j) = 1.00 dp[aa_vv](b,c) * t2[aaaa_vvoo](c,a,i,j) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vvoo_45"]("b,a,i,j")  = dp["aa_vv"]("b,c") * t2["aaaa_vvoo"]("c,a,i,j");
    
        // rt2_aaaa += -1.00 P(a,b) d-_aa(a,c) u0 t2_aaaa(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_45"]("a,b,i,j") * u0;
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_45"]("b,a,i,j") * u0;
    
        // ru2_aaaa += -1.00 P(a,b) d+_aa(a,c) t2_aaaa(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_45"]("a,b,i,j");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_45"]("b,a,i,j");
        tmps_["aaaa_vvoo_45"].~TArrayD();
    
        // tmps_[aaaa_vvoo_46](a,b,i,j) = 1.00 dp[aa_vv](a,c) * u2[aaaa_vvoo](c,b,i,j) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vvoo_46"]("a,b,i,j")  = dp["aa_vv"]("a,c") * u2["aaaa_vvoo"]("c,b,i,j");
    
        // rt2_aaaa += -1.00 P(a,b) d-_aa(a,c) u2_aaaa(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_46"]("a,b,i,j");
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_46"]("b,a,i,j");
    
        // ru2_aaaa += -1.00 P(a,b) d-_aa(a,c) u0 u2_aaaa(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_46"]("a,b,i,j") * u0;
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_46"]("b,a,i,j") * u0;
        tmps_["aaaa_vvoo_46"].~TArrayD();
    
        // tmps_[bbbb_vvoo_47](a,b,i,j) = 1.00 dp[bb_vv](a,c) * t2[bbbb_vvoo](c,b,i,j) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vvoo_47"]("a,b,i,j")  = dp["bb_vv"]("a,c") * t2["bbbb_vvoo"]("c,b,i,j");
    
        // rt2_bbbb += -1.00 P(a,b) d-_bb(a,c) u0 t2_bbbb(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_47"]("a,b,i,j") * u0;
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_47"]("b,a,i,j") * u0;
    
        // ru2_bbbb += -1.00 P(a,b) d+_bb(a,c) t2_bbbb(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_47"]("a,b,i,j");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_47"]("b,a,i,j");
        tmps_["bbbb_vvoo_47"].~TArrayD();
    
        // tmps_[bbbb_vvoo_48](a,b,i,j) = 1.00 dp[bb_vv](a,c) * u2[bbbb_vvoo](c,b,i,j) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vvoo_48"]("a,b,i,j")  = dp["bb_vv"]("a,c") * u2["bbbb_vvoo"]("c,b,i,j");
    
        // rt2_bbbb += -1.00 P(a,b) d-_bb(a,c) u2_bbbb(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_48"]("a,b,i,j");
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_48"]("b,a,i,j");
    
        // ru2_bbbb += -1.00 P(a,b) d-_bb(a,c) u0 u2_bbbb(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_48"]("a,b,i,j") * u0;
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_48"]("b,a,i,j") * u0;
        tmps_["bbbb_vvoo_48"].~TArrayD();
    
        // tmps_[aa_vv_49](c,b) = 0.50 V_blks_[aaaa_oovv](l,k,c,d) * u2[aaaa_vvoo](d,b,l,k) // flops: o0v2 = o2v3 | mem: o0v2 = o0v2
        tmps_["aa_vv_49"]("c,b")  = 0.50 * V_blks_["aaaa_oovv"]("l,k,c,d") * u2["aaaa_vvoo"]("d,b,l,k");
    
        // ru2_abab += +0.50 <l,k||c,d>_aaaa u2_aaaa(d,a,l,k) t2_abab(c,b,i,j)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["aa_vv_49"]("c,a") * t2["abab_vvoo"]("c,b,i,j");
    
        // tmps_[aaaa_voov_125](a,i,j,b) = 1.00 t2[aaaa_vvoo](c,a,i,j) * V_blks_[aaaa_oovv](l,k,c,d) * u2[aaaa_vvoo](d,b,l,k) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["aaaa_voov_125"]("a,i,j,b")  = t2["aaaa_vvoo"]("c,a,i,j") * tmps_["aa_vv_49"]("c,b");
        tmps_["aa_vv_49"].~TArrayD();
    
        // ru2_aaaa += -0.50 P(a,b) <l,k||c,d>_aaaa u2_aaaa(d,b,l,k) t2_aaaa(c,a,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_voov_125"]("a,i,j,b");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_voov_125"]("b,i,j,a");
        tmps_["aaaa_voov_125"].~TArrayD();
    
        // tmps_[bb_vv_50](c,b) = 0.50 V_blks_[bbbb_oovv](l,k,c,d) * u2[bbbb_vvoo](d,b,l,k) // flops: o0v2 = o2v3 | mem: o0v2 = o0v2
        tmps_["bb_vv_50"]("c,b")  = 0.50 * V_blks_["bbbb_oovv"]("l,k,c,d") * u2["bbbb_vvoo"]("d,b,l,k");
    
        // ru2_abab += +0.50 <l,k||c,d>_bbbb u2_bbbb(d,b,l,k) t2_abab(a,c,i,j)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["bb_vv_50"]("c,b") * t2["abab_vvoo"]("a,c,i,j");
    
        // tmps_[bbbb_voov_119](a,i,j,b) = 1.00 t2[bbbb_vvoo](c,a,i,j) * V_blks_[bbbb_oovv](l,k,c,d) * u2[bbbb_vvoo](d,b,l,k) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["bbbb_voov_119"]("a,i,j,b")  = t2["bbbb_vvoo"]("c,a,i,j") * tmps_["bb_vv_50"]("c,b");
        tmps_["bb_vv_50"].~TArrayD();
    
        // ru2_bbbb += -0.50 P(a,b) <l,k||c,d>_bbbb u2_bbbb(d,b,l,k) t2_bbbb(c,a,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_voov_119"]("a,i,j,b");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_voov_119"]("b,i,j,a");
        tmps_["bbbb_voov_119"].~TArrayD();
    
        // tmps_[aa_vo_51](a,i) = 0.50 V_blks_[abab_vovv](a,j,c,b) * t2[abab_vvoo](c,b,i,j) // flops: o1v1 = o2v3 | mem: o1v1 = o1v1
        tmps_["aa_vo_51"]("a,i")  = 0.50 * V_blks_["abab_vovv"]("a,j,c,b") * t2["abab_vvoo"]("c,b,i,j");
    
        // rt1_aa += +0.50 <a,j||b,c>_abab t2_abab(b,c,i,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_aa("a,i") += tmps_["aa_vo_51"]("a,i");
    
        // rt1_aa += +0.50 <a,j||c,b>_abab t2_abab(c,b,i,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_aa("a,i") += tmps_["aa_vo_51"]("a,i");
        tmps_["aa_vo_51"].~TArrayD();
    
        // tmps_[bb_vo_52](a,i) = 0.50 V_blks_[baab_vovv](a,j,b,c) * t2[abab_vvoo](b,c,j,i) // flops: o1v1 = o2v3 | mem: o1v1 = o1v1
        tmps_["bb_vo_52"]("a,i")  = 0.50 * V_blks_["baab_vovv"]("a,j,b,c") * t2["abab_vvoo"]("b,c,j,i");
    
        // rt1_bb += +0.50 <j,a||b,c>_abab t2_abab(b,c,j,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_bb("a,i") -= tmps_["bb_vo_52"]("a,i");
    
        // rt1_bb += +0.50 <j,a||c,b>_abab t2_abab(c,b,j,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_bb("a,i") -= tmps_["bb_vo_52"]("a,i");
        tmps_["bb_vo_52"].~TArrayD();
    
        // tmps_[aa_vo_53](a,i) = 0.50 V_blks_[abab_vovv](a,j,c,b) * u2[abab_vvoo](c,b,i,j) // flops: o1v1 = o2v3 | mem: o1v1 = o1v1
        tmps_["aa_vo_53"]("a,i")  = 0.50 * V_blks_["abab_vovv"]("a,j,c,b") * u2["abab_vvoo"]("c,b,i,j");
    
        // ru1_aa += +0.50 <a,j||b,c>_abab u2_abab(b,c,i,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") += tmps_["aa_vo_53"]("a,i");
    
        // ru1_aa += +0.50 <a,j||c,b>_abab u2_abab(c,b,i,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") += tmps_["aa_vo_53"]("a,i");
        tmps_["aa_vo_53"].~TArrayD();
    
        // tmps_[bb_vo_54](a,i) = 0.50 V_blks_[baab_vovv](a,j,b,c) * u2[abab_vvoo](b,c,j,i) // flops: o1v1 = o2v3 | mem: o1v1 = o1v1
        tmps_["bb_vo_54"]("a,i")  = 0.50 * V_blks_["baab_vovv"]("a,j,b,c") * u2["abab_vvoo"]("b,c,j,i");
    
        // ru1_bb += +0.50 <j,a||b,c>_abab u2_abab(b,c,j,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= tmps_["bb_vo_54"]("a,i");
    
        // ru1_bb += +0.50 <j,a||c,b>_abab u2_abab(c,b,j,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= tmps_["bb_vo_54"]("a,i");
        tmps_["bb_vo_54"].~TArrayD();
    
        // tmps_[bb_vv_55](c,b) = 0.50 V_blks_[bbbb_oovv](l,k,c,d) * t2[bbbb_vvoo](d,b,l,k) // flops: o0v2 = o2v3 | mem: o0v2 = o0v2
        tmps_["bb_vv_55"]("c,b")  = 0.50 * V_blks_["bbbb_oovv"]("l,k,c,d") * t2["bbbb_vvoo"]("d,b,l,k");
    
        // rt2_abab += +0.50 <l,k||c,d>_bbbb t2_abab(a,c,i,j) t2_bbbb(d,b,l,k)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["bb_vv_55"]("c,b") * t2["abab_vvoo"]("a,c,i,j");
    
        // rt2_bbbb += -0.50 <l,k||c,d>_bbbb t2_bbbb(c,a,i,j) t2_bbbb(d,b,l,k)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") -= tmps_["bb_vv_55"]("c,b") * t2["bbbb_vvoo"]("c,a,i,j");
        tmps_["bb_vv_55"].~TArrayD();
    
        // tmps_[abab_vvoo_56](a,b,i,j) = 1.00 dp[aa_vv](a,c) * u2[abab_vvoo](c,b,i,j) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["abab_vvoo_56"]("a,b,i,j")  = dp["aa_vv"]("a,c") * u2["abab_vvoo"]("c,b,i,j");
    
        // rt2_abab += -1.00 d-_aa(a,c) u2_abab(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= tmps_["abab_vvoo_56"]("a,b,i,j");
    
        // ru2_abab += -1.00 d-_aa(a,c) u0 u2_abab(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["abab_vvoo_56"]("a,b,i,j") * u0;
        tmps_["abab_vvoo_56"].~TArrayD();
    
        // tmps_[baab_vvoo_57](b,a,i,j) = 1.00 dp[bb_vv](b,c) * u2[abab_vvoo](a,c,i,j) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["baab_vvoo_57"]("b,a,i,j")  = dp["bb_vv"]("b,c") * u2["abab_vvoo"]("a,c,i,j");
    
        // rt2_abab += -1.00 d-_bb(b,c) u2_abab(a,c,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= tmps_["baab_vvoo_57"]("b,a,i,j");
    
        // ru2_abab += -1.00 d-_bb(b,c) u0 u2_abab(a,c,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["baab_vvoo_57"]("b,a,i,j") * u0;
        tmps_["baab_vvoo_57"].~TArrayD();
    
        // tmps_[baab_vvoo_58](b,a,i,j) = 1.00 dp[bb_vv](b,c) * t2[abab_vvoo](a,c,i,j) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["baab_vvoo_58"]("b,a,i,j")  = dp["bb_vv"]("b,c") * t2["abab_vvoo"]("a,c,i,j");
    
        // rt2_abab += -1.00 d-_bb(b,c) u0 t2_abab(a,c,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= tmps_["baab_vvoo_58"]("b,a,i,j") * u0;
    
        // ru2_abab += -1.00 d+_bb(b,c) t2_abab(a,c,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["baab_vvoo_58"]("b,a,i,j");
        tmps_["baab_vvoo_58"].~TArrayD();
    
        // tmps_[abab_vvoo_59](a,b,i,j) = 1.00 dp[aa_vv](a,c) * t2[abab_vvoo](c,b,i,j) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["abab_vvoo_59"]("a,b,i,j")  = dp["aa_vv"]("a,c") * t2["abab_vvoo"]("c,b,i,j");
    
        // rt2_abab += -1.00 d-_aa(a,c) u0 t2_abab(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= tmps_["abab_vvoo_59"]("a,b,i,j") * u0;
    
        // ru2_abab += -1.00 d+_aa(a,c) t2_abab(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["abab_vvoo_59"]("a,b,i,j");
        tmps_["abab_vvoo_59"].~TArrayD();
    
        // tmps_[bbbb_vvoo_60](a,b,i,j) = 1.00 F_blks_[bb_vv](a,c) * u2[bbbb_vvoo](c,b,i,j) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vvoo_60"]("a,b,i,j")  = F_blks_["bb_vv"]("a,c") * u2["bbbb_vvoo"]("c,b,i,j");
    
        // ru2_bbbb += +1.00 P(a,b) f_bb(a,c) u2_bbbb(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_60"]("a,b,i,j");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_60"]("b,a,i,j");
        tmps_["bbbb_vvoo_60"].~TArrayD();
    
        // tmps_[bbbb_vvoo_61](a,b,j,i) = 1.00 V_blks_[bbbb_vvvo](a,b,c,j) * u1[bb_vo](c,i) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vvoo_61"]("a,b,j,i")  = V_blks_["bbbb_vvvo"]("a,b,c,j") * u1["bb_vo"]("c,i");
    
        // ru2_bbbb += +1.00 P(i,j) <a,b||c,j>_bbbb u1_bb(c,i)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_61"]("a,b,j,i");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_61"]("a,b,i,j");
        tmps_["bbbb_vvoo_61"].~TArrayD();
    
        // tmps_[aaaa_vvoo_62](a,b,j,i) = 1.00 V_blks_[aaaa_vvvo](a,b,c,j) * u1[aa_vo](c,i) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vvoo_62"]("a,b,j,i")  = V_blks_["aaaa_vvvo"]("a,b,c,j") * u1["aa_vo"]("c,i");
    
        // ru2_aaaa += +1.00 P(i,j) <a,b||c,j>_aaaa u1_aa(c,i)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_62"]("a,b,j,i");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_62"]("a,b,i,j");
        tmps_["aaaa_vvoo_62"].~TArrayD();
    
        // tmps_[bbbb_vvoo_63](a,b,i,j) = 1.00 F_blks_[bb_vv](a,c) * t2[bbbb_vvoo](c,b,i,j) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["bbbb_vvoo_63"]("a,b,i,j")  = F_blks_["bb_vv"]("a,c") * t2["bbbb_vvoo"]("c,b,i,j");
    
        // rt2_bbbb += +1.00 P(a,b) f_bb(a,c) t2_bbbb(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_63"]("a,b,i,j");
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_63"]("b,a,i,j");
        tmps_["bbbb_vvoo_63"].~TArrayD();
    
        // tmps_[aaaa_vvoo_64](a,b,i,j) = 1.00 F_blks_[aa_vv](a,c) * u2[aaaa_vvoo](c,b,i,j) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vvoo_64"]("a,b,i,j")  = F_blks_["aa_vv"]("a,c") * u2["aaaa_vvoo"]("c,b,i,j");
    
        // ru2_aaaa += +1.00 P(a,b) f_aa(a,c) u2_aaaa(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_64"]("a,b,i,j");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_64"]("b,a,i,j");
        tmps_["aaaa_vvoo_64"].~TArrayD();
    
        // tmps_[aaaa_vvoo_65](a,b,i,j) = 1.00 F_blks_[aa_vv](a,c) * t2[aaaa_vvoo](c,b,i,j) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vvoo_65"]("a,b,i,j")  = F_blks_["aa_vv"]("a,c") * t2["aaaa_vvoo"]("c,b,i,j");
    
        // rt2_aaaa += +1.00 P(a,b) f_aa(a,c) t2_aaaa(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_65"]("a,b,i,j");
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_65"]("b,a,i,j");
        tmps_["aaaa_vvoo_65"].~TArrayD();
    
        // tmps_[aa_oo_66](l,j) = 0.50 V_blks_[abab_oovv](l,k,c,d) * t2[abab_vvoo](c,d,j,k) // flops: o2v0 = o3v2 | mem: o2v0 = o2v0
        tmps_["aa_oo_66"]("l,j")  = 0.50 * V_blks_["abab_oovv"]("l,k,c,d") * t2["abab_vvoo"]("c,d,j,k");
    
        // tmps_[aaaa_vvoo_134](a,b,i,j) = 1.00 u2[aaaa_vvoo](a,b,i,l) * V_blks_[abab_oovv](l,k,c,d) * t2[abab_vvoo](c,d,j,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["aaaa_vvoo_134"]("a,b,i,j")  = u2["aaaa_vvoo"]("a,b,i,l") * tmps_["aa_oo_66"]("l,j");
    
        // ru2_aaaa += -0.50 P(i,j) <l,k||c,d>_abab u2_aaaa(a,b,i,l) t2_abab(c,d,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_134"]("a,b,i,j");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_134"]("a,b,j,i");
    
        // ru2_aaaa += -0.50 P(i,j) <l,k||d,c>_abab u2_aaaa(a,b,i,l) t2_abab(d,c,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_134"]("a,b,i,j");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_134"]("a,b,j,i");
        tmps_["aaaa_vvoo_134"].~TArrayD();
    
        // tmps_[aaaa_vvoo_139](a,b,i,j) = 1.00 t2[aaaa_vvoo](a,b,i,l) * V_blks_[abab_oovv](l,k,c,d) * t2[abab_vvoo](c,d,j,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["aaaa_vvoo_139"]("a,b,i,j")  = t2["aaaa_vvoo"]("a,b,i,l") * tmps_["aa_oo_66"]("l,j");
    
        // rt2_aaaa += -0.50 P(i,j) <l,k||c,d>_abab t2_abab(c,d,j,k) t2_aaaa(a,b,i,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_139"]("a,b,i,j");
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_139"]("a,b,j,i");
    
        // rt2_aaaa += -0.50 P(i,j) <l,k||d,c>_abab t2_abab(d,c,j,k) t2_aaaa(a,b,i,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_139"]("a,b,i,j");
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_139"]("a,b,j,i");
        tmps_["aaaa_vvoo_139"].~TArrayD();
    
        // tmps_[abba_vvoo_156](a,b,j,i) = 1.00 u2[abab_vvoo](a,b,l,j) * V_blks_[abab_oovv](l,k,c,d) * t2[abab_vvoo](c,d,i,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["abba_vvoo_156"]("a,b,j,i")  = u2["abab_vvoo"]("a,b,l,j") * tmps_["aa_oo_66"]("l,i");
    
        // ru2_abab += -0.50 <l,k||c,d>_abab u2_abab(a,b,l,j) t2_abab(c,d,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["abba_vvoo_156"]("a,b,j,i");
    
        // ru2_abab += -0.50 <l,k||d,c>_abab u2_abab(a,b,l,j) t2_abab(d,c,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["abba_vvoo_156"]("a,b,j,i");
        tmps_["abba_vvoo_156"].~TArrayD();
    
        // tmps_[abba_vvoo_166](a,b,j,i) = 1.00 t2[abab_vvoo](a,b,l,j) * V_blks_[abab_oovv](l,k,c,d) * t2[abab_vvoo](c,d,i,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["abba_vvoo_166"]("a,b,j,i")  = t2["abab_vvoo"]("a,b,l,j") * tmps_["aa_oo_66"]("l,i");
    
        // rt2_abab += -0.50 <l,k||c,d>_abab t2_abab(c,d,i,k) t2_abab(a,b,l,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= tmps_["abba_vvoo_166"]("a,b,j,i");
    
        // rt2_abab += -0.50 <l,k||d,c>_abab t2_abab(d,c,i,k) t2_abab(a,b,l,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= tmps_["abba_vvoo_166"]("a,b,j,i");
        tmps_["abba_vvoo_166"].~TArrayD();
    
        // tmps_[aa_vo_215](a,i) = 1.00 u1[aa_vo](a,k) * V_blks_[abab_oovv](k,j,b,c) * t2[abab_vvoo](b,c,i,j) // flops: o1v1 = o2v1 | mem: o1v1 = o1v1
        tmps_["aa_vo_215"]("a,i")  = u1["aa_vo"]("a,k") * tmps_["aa_oo_66"]("k,i");
        tmps_["aa_oo_66"].~TArrayD();
    
        // ru1_aa += -0.50 <k,j||b,c>_abab u1_aa(a,k) t2_abab(b,c,i,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= tmps_["aa_vo_215"]("a,i");
    
        // ru1_aa += -0.50 <k,j||c,b>_abab u1_aa(a,k) t2_abab(c,b,i,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= tmps_["aa_vo_215"]("a,i");
        tmps_["aa_vo_215"].~TArrayD();
    
        // tmps_[bb_oo_67](l,j) = 0.50 V_blks_[abab_oovv](k,l,d,c) * t2[abab_vvoo](d,c,k,j) // flops: o2v0 = o3v2 | mem: o2v0 = o2v0
        tmps_["bb_oo_67"]("l,j")  = 0.50 * V_blks_["abab_oovv"]("k,l,d,c") * t2["abab_vvoo"]("d,c,k,j");
    
        // tmps_[bbbb_ovvo_131](j,a,b,i) = 1.00 V_blks_[abab_oovv](k,l,c,d) * t2[abab_vvoo](c,d,k,j) * u2[bbbb_vvoo](a,b,i,l) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["bbbb_ovvo_131"]("j,a,b,i")  = tmps_["bb_oo_67"]("l,j") * u2["bbbb_vvoo"]("a,b,i,l");
    
        // ru2_bbbb += -0.50 P(i,j) <k,l||c,d>_abab u2_bbbb(a,b,i,l) t2_abab(c,d,k,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_ovvo_131"]("j,a,b,i");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_ovvo_131"]("i,a,b,j");
    
        // ru2_bbbb += -0.50 P(i,j) <k,l||d,c>_abab u2_bbbb(a,b,i,l) t2_abab(d,c,k,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_ovvo_131"]("j,a,b,i");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_ovvo_131"]("i,a,b,j");
        tmps_["bbbb_ovvo_131"].~TArrayD();
    
        // tmps_[bbbb_vvoo_137](a,b,i,j) = 1.00 t2[bbbb_vvoo](a,b,i,l) * V_blks_[abab_oovv](k,l,c,d) * t2[abab_vvoo](c,d,k,j) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["bbbb_vvoo_137"]("a,b,i,j")  = t2["bbbb_vvoo"]("a,b,i,l") * tmps_["bb_oo_67"]("l,j");
    
        // rt2_bbbb += -0.50 P(i,j) <k,l||c,d>_abab t2_abab(c,d,k,j) t2_bbbb(a,b,i,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_137"]("a,b,i,j");
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_137"]("a,b,j,i");
    
        // rt2_bbbb += -0.50 P(i,j) <k,l||d,c>_abab t2_abab(d,c,k,j) t2_bbbb(a,b,i,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_137"]("a,b,i,j");
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_137"]("a,b,j,i");
        tmps_["bbbb_vvoo_137"].~TArrayD();
    
        // tmps_[abab_vvoo_155](a,b,i,j) = 1.00 u2[abab_vvoo](a,b,i,l) * V_blks_[abab_oovv](k,l,c,d) * t2[abab_vvoo](c,d,k,j) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["abab_vvoo_155"]("a,b,i,j")  = u2["abab_vvoo"]("a,b,i,l") * tmps_["bb_oo_67"]("l,j");
    
        // ru2_abab += -0.50 <k,l||c,d>_abab u2_abab(a,b,i,l) t2_abab(c,d,k,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["abab_vvoo_155"]("a,b,i,j");
    
        // ru2_abab += -0.50 <k,l||d,c>_abab u2_abab(a,b,i,l) t2_abab(d,c,k,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["abab_vvoo_155"]("a,b,i,j");
        tmps_["abab_vvoo_155"].~TArrayD();
    
        // tmps_[abab_vvoo_167](a,b,i,j) = 1.00 t2[abab_vvoo](a,b,i,l) * V_blks_[abab_oovv](k,l,c,d) * t2[abab_vvoo](c,d,k,j) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["abab_vvoo_167"]("a,b,i,j")  = t2["abab_vvoo"]("a,b,i,l") * tmps_["bb_oo_67"]("l,j");
    
        // rt2_abab += -0.50 <k,l||c,d>_abab t2_abab(c,d,k,j) t2_abab(a,b,i,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= tmps_["abab_vvoo_167"]("a,b,i,j");
    
        // rt2_abab += -0.50 <k,l||d,c>_abab t2_abab(d,c,k,j) t2_abab(a,b,i,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= tmps_["abab_vvoo_167"]("a,b,i,j");
        tmps_["abab_vvoo_167"].~TArrayD();
    
        // tmps_[bb_ov_214](i,a) = 1.00 V_blks_[abab_oovv](j,k,c,b) * t2[abab_vvoo](c,b,j,i) * u1[bb_vo](a,k) // flops: o1v1 = o2v1 | mem: o1v1 = o1v1
        tmps_["bb_ov_214"]("i,a")  = tmps_["bb_oo_67"]("k,i") * u1["bb_vo"]("a,k");
        tmps_["bb_oo_67"].~TArrayD();
    
        // ru1_bb += -0.50 <j,k||b,c>_abab u1_bb(a,k) t2_abab(b,c,j,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= tmps_["bb_ov_214"]("i,a");
    
        // ru1_bb += -0.50 <j,k||c,b>_abab u1_bb(a,k) t2_abab(c,b,j,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= tmps_["bb_ov_214"]("i,a");
        tmps_["bb_ov_214"].~TArrayD();
    
        // tmps_[bb_oo_68](l,j) = 0.50 V_blks_[bbbb_oovv](l,k,c,d) * t2[bbbb_vvoo](c,d,j,k) // flops: o2v0 = o3v2 | mem: o2v0 = o2v0
        tmps_["bb_oo_68"]("l,j")  = 0.50 * V_blks_["bbbb_oovv"]("l,k,c,d") * t2["bbbb_vvoo"]("c,d,j,k");
    
        // rt2_abab += -0.50 <l,k||c,d>_bbbb t2_bbbb(c,d,j,k) t2_abab(a,b,i,l)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= tmps_["bb_oo_68"]("l,j") * t2["abab_vvoo"]("a,b,i,l");
    
        // ru1_bb += -0.50 <k,j||b,c>_bbbb u1_bb(a,k) t2_bbbb(b,c,i,j)  // flops: o1v1 += o2v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= tmps_["bb_oo_68"]("k,i") * u1["bb_vo"]("a,k");
    
        // ru2_abab += -0.50 <l,k||c,d>_bbbb u2_abab(a,b,i,l) t2_bbbb(c,d,j,k)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["bb_oo_68"]("l,j") * u2["abab_vvoo"]("a,b,i,l");
    
        // tmps_[bbbb_vvoo_148](a,b,j,i) = 1.00 u2[bbbb_vvoo](a,b,j,l) * V_blks_[bbbb_oovv](l,k,c,d) * t2[bbbb_vvoo](c,d,i,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["bbbb_vvoo_148"]("a,b,j,i")  = u2["bbbb_vvoo"]("a,b,j,l") * tmps_["bb_oo_68"]("l,i");
    
        // ru2_bbbb += -0.50 P(i,j) <l,k||c,d>_bbbb u2_bbbb(a,b,i,l) t2_bbbb(c,d,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_148"]("a,b,i,j");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_148"]("a,b,j,i");
        tmps_["bbbb_vvoo_148"].~TArrayD();
    
        // tmps_[bbbb_vvoo_165](a,b,i,j) = 1.00 t2[bbbb_vvoo](a,b,i,l) * V_blks_[bbbb_oovv](l,k,c,d) * t2[bbbb_vvoo](c,d,j,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["bbbb_vvoo_165"]("a,b,i,j")  = t2["bbbb_vvoo"]("a,b,i,l") * tmps_["bb_oo_68"]("l,j");
        tmps_["bb_oo_68"].~TArrayD();
    
        // rt2_bbbb += -0.50 P(i,j) <l,k||c,d>_bbbb t2_bbbb(c,d,j,k) t2_bbbb(a,b,i,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_165"]("a,b,i,j");
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_165"]("a,b,j,i");
        tmps_["bbbb_vvoo_165"].~TArrayD();
    
        // tmps_[aa_oo_69](l,j) = 0.50 V_blks_[aaaa_oovv](l,k,c,d) * t2[aaaa_vvoo](c,d,j,k) // flops: o2v0 = o3v2 | mem: o2v0 = o2v0
        tmps_["aa_oo_69"]("l,j")  = 0.50 * V_blks_["aaaa_oovv"]("l,k,c,d") * t2["aaaa_vvoo"]("c,d,j,k");
    
        // rt2_abab += -0.50 <l,k||c,d>_aaaa t2_aaaa(c,d,i,k) t2_abab(a,b,l,j)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= tmps_["aa_oo_69"]("l,i") * t2["abab_vvoo"]("a,b,l,j");
    
        // ru1_aa += -0.50 <k,j||b,c>_aaaa u1_aa(a,k) t2_aaaa(b,c,i,j)  // flops: o1v1 += o2v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= tmps_["aa_oo_69"]("k,i") * u1["aa_vo"]("a,k");
    
        // ru2_abab += -0.50 <l,k||c,d>_aaaa u2_abab(a,b,l,j) t2_aaaa(c,d,i,k)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["aa_oo_69"]("l,i") * u2["abab_vvoo"]("a,b,l,j");
    
        // tmps_[aaaa_vvoo_161](a,b,i,j) = 1.00 u2[aaaa_vvoo](a,b,i,l) * V_blks_[aaaa_oovv](l,k,c,d) * t2[aaaa_vvoo](c,d,j,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["aaaa_vvoo_161"]("a,b,i,j")  = u2["aaaa_vvoo"]("a,b,i,l") * tmps_["aa_oo_69"]("l,j");
    
        // ru2_aaaa += -0.50 P(i,j) <l,k||c,d>_aaaa u2_aaaa(a,b,i,l) t2_aaaa(c,d,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_161"]("a,b,i,j");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_161"]("a,b,j,i");
        tmps_["aaaa_vvoo_161"].~TArrayD();
    
        // tmps_[aaaa_vvoo_170](a,b,i,j) = 1.00 t2[aaaa_vvoo](a,b,i,l) * V_blks_[aaaa_oovv](l,k,c,d) * t2[aaaa_vvoo](c,d,j,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["aaaa_vvoo_170"]("a,b,i,j")  = t2["aaaa_vvoo"]("a,b,i,l") * tmps_["aa_oo_69"]("l,j");
        tmps_["aa_oo_69"].~TArrayD();
    
        // rt2_aaaa += -0.50 P(i,j) <l,k||c,d>_aaaa t2_aaaa(c,d,j,k) t2_aaaa(a,b,i,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_170"]("a,b,i,j");
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_170"]("a,b,j,i");
        tmps_["aaaa_vvoo_170"].~TArrayD();
    
        // tmps_[bb_oo_70](k,i) = 0.50 V_blks_[abab_oovv](l,k,c,d) * u2[abab_vvoo](c,d,l,i) // flops: o2v0 = o3v2 | mem: o2v0 = o2v0
        tmps_["bb_oo_70"]("k,i")  = 0.50 * V_blks_["abab_oovv"]("l,k,c,d") * u2["abab_vvoo"]("c,d,l,i");
    
        // tmps_[bbbb_vvoo_132](a,b,i,j) = 1.00 t2[bbbb_vvoo](a,b,i,k) * V_blks_[abab_oovv](l,k,d,c) * u2[abab_vvoo](d,c,l,j) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["bbbb_vvoo_132"]("a,b,i,j")  = t2["bbbb_vvoo"]("a,b,i,k") * tmps_["bb_oo_70"]("k,j");
    
        // ru2_bbbb += +0.50 P(i,j) <l,k||c,d>_abab u2_abab(c,d,l,i) t2_bbbb(a,b,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_132"]("a,b,j,i");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_132"]("a,b,i,j");
    
        // ru2_bbbb += +0.50 P(i,j) <l,k||d,c>_abab u2_abab(d,c,l,i) t2_bbbb(a,b,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_132"]("a,b,j,i");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_132"]("a,b,i,j");
        tmps_["bbbb_vvoo_132"].~TArrayD();
    
        // tmps_[baba_ovvo_154](j,a,b,i) = 1.00 V_blks_[abab_oovv](l,k,d,c) * u2[abab_vvoo](d,c,l,j) * t2[abab_vvoo](a,b,i,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["baba_ovvo_154"]("j,a,b,i")  = tmps_["bb_oo_70"]("k,j") * t2["abab_vvoo"]("a,b,i,k");
        tmps_["bb_oo_70"].~TArrayD();
    
        // ru2_abab += -0.50 <l,k||c,d>_abab u2_abab(c,d,l,j) t2_abab(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["baba_ovvo_154"]("j,a,b,i");
    
        // ru2_abab += -0.50 <l,k||d,c>_abab u2_abab(d,c,l,j) t2_abab(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["baba_ovvo_154"]("j,a,b,i");
        tmps_["baba_ovvo_154"].~TArrayD();
    
        // tmps_[aa_oo_71](k,j) = 0.50 V_blks_[abab_oovv](k,l,c,d) * u2[abab_vvoo](c,d,j,l) // flops: o2v0 = o3v2 | mem: o2v0 = o2v0
        tmps_["aa_oo_71"]("k,j")  = 0.50 * V_blks_["abab_oovv"]("k,l,c,d") * u2["abab_vvoo"]("c,d,j,l");
    
        // tmps_[aaaa_vvoo_135](a,b,i,j) = 1.00 t2[aaaa_vvoo](a,b,i,k) * V_blks_[abab_oovv](k,l,c,d) * u2[abab_vvoo](c,d,j,l) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["aaaa_vvoo_135"]("a,b,i,j")  = t2["aaaa_vvoo"]("a,b,i,k") * tmps_["aa_oo_71"]("k,j");
    
        // ru2_aaaa += +0.50 P(i,j) <k,l||c,d>_abab u2_abab(c,d,i,l) t2_aaaa(a,b,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_135"]("a,b,j,i");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_135"]("a,b,i,j");
    
        // ru2_aaaa += +0.50 P(i,j) <k,l||d,c>_abab u2_abab(d,c,i,l) t2_aaaa(a,b,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_135"]("a,b,j,i");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_135"]("a,b,i,j");
        tmps_["aaaa_vvoo_135"].~TArrayD();
    
        // tmps_[abba_vvoo_157](a,b,j,i) = 1.00 t2[abab_vvoo](a,b,k,j) * V_blks_[abab_oovv](k,l,d,c) * u2[abab_vvoo](d,c,i,l) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["abba_vvoo_157"]("a,b,j,i")  = t2["abab_vvoo"]("a,b,k,j") * tmps_["aa_oo_71"]("k,i");
        tmps_["aa_oo_71"].~TArrayD();
    
        // ru2_abab += -0.50 <k,l||c,d>_abab u2_abab(c,d,i,l) t2_abab(a,b,k,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["abba_vvoo_157"]("a,b,j,i");
    
        // ru2_abab += -0.50 <k,l||d,c>_abab u2_abab(d,c,i,l) t2_abab(a,b,k,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["abba_vvoo_157"]("a,b,j,i");
        tmps_["abba_vvoo_157"].~TArrayD();
    
        // tmps_[bbbb_ovvo_74](j,a,b,i) = 1.00 dp[bb_oo](k,j) * t2[bbbb_vvoo](a,b,i,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["bbbb_ovvo_74"]("j,a,b,i")  = dp["bb_oo"]("k,j") * t2["bbbb_vvoo"]("a,b,i,k");
    
        // rt2_bbbb += +1.00 P(i,j) d-_bb(k,j) u0 t2_bbbb(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_ovvo_74"]("j,a,b,i") * u0;
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_ovvo_74"]("i,a,b,j") * u0;
    
        // ru2_bbbb += +1.00 P(i,j) d+_bb(k,j) t2_bbbb(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_ovvo_74"]("j,a,b,i");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_ovvo_74"]("i,a,b,j");
        tmps_["bbbb_ovvo_74"].~TArrayD();
    
        // tmps_[aaaa_ovvo_75](j,a,b,i) = 1.00 dp[aa_oo](k,j) * u2[aaaa_vvoo](a,b,i,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["aaaa_ovvo_75"]("j,a,b,i")  = dp["aa_oo"]("k,j") * u2["aaaa_vvoo"]("a,b,i,k");
    
        // rt2_aaaa += +1.00 P(i,j) d-_aa(k,j) u2_aaaa(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_ovvo_75"]("j,a,b,i");
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_ovvo_75"]("i,a,b,j");
    
        // ru2_aaaa += +1.00 P(i,j) d-_aa(k,j) u0 u2_aaaa(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_ovvo_75"]("j,a,b,i") * u0;
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_ovvo_75"]("i,a,b,j") * u0;
        tmps_["aaaa_ovvo_75"].~TArrayD();
    
        // tmps_[bbbb_ovvo_76](j,a,b,i) = 1.00 dp[bb_oo](k,j) * u2[bbbb_vvoo](a,b,i,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["bbbb_ovvo_76"]("j,a,b,i")  = dp["bb_oo"]("k,j") * u2["bbbb_vvoo"]("a,b,i,k");
    
        // rt2_bbbb += +1.00 P(i,j) d-_bb(k,j) u2_bbbb(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_ovvo_76"]("j,a,b,i");
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_ovvo_76"]("i,a,b,j");
    
        // ru2_bbbb += +1.00 P(i,j) d-_bb(k,j) u0 u2_bbbb(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_ovvo_76"]("j,a,b,i") * u0;
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_ovvo_76"]("i,a,b,j") * u0;
        tmps_["bbbb_ovvo_76"].~TArrayD();
    
        // tmps_[aaaa_ovvo_77](j,a,b,i) = 1.00 dp[aa_oo](k,j) * t2[aaaa_vvoo](a,b,i,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["aaaa_ovvo_77"]("j,a,b,i")  = dp["aa_oo"]("k,j") * t2["aaaa_vvoo"]("a,b,i,k");
    
        // rt2_aaaa += +1.00 P(i,j) d-_aa(k,j) u0 t2_aaaa(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_ovvo_77"]("j,a,b,i") * u0;
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_ovvo_77"]("i,a,b,j") * u0;
    
        // ru2_aaaa += +1.00 P(i,j) d+_aa(k,j) t2_aaaa(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_ovvo_77"]("j,a,b,i");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_ovvo_77"]("i,a,b,j");
        tmps_["aaaa_ovvo_77"].~TArrayD();
    
        // tmps_[bb_oo_78](k,i) = 0.50 V_blks_[bbbb_oovv](l,k,c,d) * u2[bbbb_vvoo](c,d,i,l) // flops: o2v0 = o3v2 | mem: o2v0 = o2v0
        tmps_["bb_oo_78"]("k,i")  = 0.50 * V_blks_["bbbb_oovv"]("l,k,c,d") * u2["bbbb_vvoo"]("c,d,i,l");
    
        // ru2_abab += +0.50 <l,k||c,d>_bbbb u2_bbbb(c,d,j,l) t2_abab(a,b,i,k)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["bb_oo_78"]("k,j") * t2["abab_vvoo"]("a,b,i,k");
    
        // tmps_[bbbb_vvoo_149](a,b,j,i) = 1.00 t2[bbbb_vvoo](a,b,j,k) * V_blks_[bbbb_oovv](l,k,c,d) * u2[bbbb_vvoo](c,d,i,l) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["bbbb_vvoo_149"]("a,b,j,i")  = t2["bbbb_vvoo"]("a,b,j,k") * tmps_["bb_oo_78"]("k,i");
        tmps_["bb_oo_78"].~TArrayD();
    
        // ru2_bbbb += -0.50 P(i,j) <l,k||c,d>_bbbb u2_bbbb(c,d,i,l) t2_bbbb(a,b,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_149"]("a,b,j,i");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_149"]("a,b,i,j");
        tmps_["bbbb_vvoo_149"].~TArrayD();
    
        // tmps_[aa_oo_79](k,i) = 0.50 V_blks_[aaaa_oovv](l,k,c,d) * u2[aaaa_vvoo](c,d,i,l) // flops: o2v0 = o3v2 | mem: o2v0 = o2v0
        tmps_["aa_oo_79"]("k,i")  = 0.50 * V_blks_["aaaa_oovv"]("l,k,c,d") * u2["aaaa_vvoo"]("c,d,i,l");
    
        // ru2_abab += +0.50 <l,k||c,d>_aaaa u2_aaaa(c,d,i,l) t2_abab(a,b,k,j)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["aa_oo_79"]("k,i") * t2["abab_vvoo"]("a,b,k,j");
    
        // tmps_[aaaa_vvoo_162](a,b,i,j) = 1.00 t2[aaaa_vvoo](a,b,i,k) * V_blks_[aaaa_oovv](l,k,c,d) * u2[aaaa_vvoo](c,d,j,l) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["aaaa_vvoo_162"]("a,b,i,j")  = t2["aaaa_vvoo"]("a,b,i,k") * tmps_["aa_oo_79"]("k,j");
        tmps_["aa_oo_79"].~TArrayD();
    
        // ru2_aaaa += -0.50 P(i,j) <l,k||c,d>_aaaa u2_aaaa(c,d,i,l) t2_aaaa(a,b,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_162"]("a,b,j,i");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_162"]("a,b,i,j");
        tmps_["aaaa_vvoo_162"].~TArrayD();
    
        // tmps_[bb_ov_80](i,a) = 0.50 V_blks_[abab_oovo](k,j,b,i) * t2[abab_vvoo](b,a,k,j) // flops: o1v1 = o3v2 | mem: o1v1 = o1v1
        tmps_["bb_ov_80"]("i,a")  = 0.50 * V_blks_["abab_oovo"]("k,j,b,i") * t2["abab_vvoo"]("b,a,k,j");
    
        // rt1_bb += -0.50 <k,j||b,i>_abab t2_abab(b,a,k,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_bb("a,i") -= tmps_["bb_ov_80"]("i,a");
    
        // rt1_bb += -0.50 <j,k||b,i>_abab t2_abab(b,a,j,k)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_bb("a,i") -= tmps_["bb_ov_80"]("i,a");
        tmps_["bb_ov_80"].~TArrayD();
    
        // tmps_[aa_ov_81](i,a) = 0.50 V_blks_[abba_oovo](j,k,b,i) * t2[abab_vvoo](a,b,j,k) // flops: o1v1 = o3v2 | mem: o1v1 = o1v1
        tmps_["aa_ov_81"]("i,a")  = 0.50 * V_blks_["abba_oovo"]("j,k,b,i") * t2["abab_vvoo"]("a,b,j,k");
    
        // rt1_aa += -0.50 <k,j||i,b>_abab t2_abab(a,b,k,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_aa("a,i") += tmps_["aa_ov_81"]("i,a");
    
        // rt1_aa += -0.50 <j,k||i,b>_abab t2_abab(a,b,j,k)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_aa("a,i") += tmps_["aa_ov_81"]("i,a");
        tmps_["aa_ov_81"].~TArrayD();
    
        // tmps_[aa_ov_82](i,a) = 0.50 V_blks_[abba_oovo](k,j,b,i) * u2[abab_vvoo](a,b,k,j) // flops: o1v1 = o3v2 | mem: o1v1 = o1v1
        tmps_["aa_ov_82"]("i,a")  = 0.50 * V_blks_["abba_oovo"]("k,j,b,i") * u2["abab_vvoo"]("a,b,k,j");
    
        // ru1_aa += -0.50 <k,j||i,b>_abab u2_abab(a,b,k,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") += tmps_["aa_ov_82"]("i,a");
    
        // ru1_aa += -0.50 <j,k||i,b>_abab u2_abab(a,b,j,k)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") += tmps_["aa_ov_82"]("i,a");
        tmps_["aa_ov_82"].~TArrayD();
    
        // tmps_[bb_ov_83](i,a) = 0.50 V_blks_[abab_oovo](j,k,b,i) * u2[abab_vvoo](b,a,j,k) // flops: o1v1 = o3v2 | mem: o1v1 = o1v1
        tmps_["bb_ov_83"]("i,a")  = 0.50 * V_blks_["abab_oovo"]("j,k,b,i") * u2["abab_vvoo"]("b,a,j,k");
    
        // ru1_bb += -0.50 <k,j||b,i>_abab u2_abab(b,a,k,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= tmps_["bb_ov_83"]("i,a");
    
        // ru1_bb += -0.50 <j,k||b,i>_abab u2_abab(b,a,j,k)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= tmps_["bb_ov_83"]("i,a");
        tmps_["bb_ov_83"].~TArrayD();
    
        // tmps_[baba_ovvo_90](j,a,b,i) = 1.00 dp[bb_oo](k,j) * u2[abab_vvoo](a,b,i,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["baba_ovvo_90"]("j,a,b,i")  = dp["bb_oo"]("k,j") * u2["abab_vvoo"]("a,b,i,k");
    
        // rt2_abab += +1.00 d-_bb(k,j) u2_abab(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["baba_ovvo_90"]("j,a,b,i");
    
        // ru2_abab += +1.00 d-_bb(k,j) u0 u2_abab(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["baba_ovvo_90"]("j,a,b,i") * u0;
        tmps_["baba_ovvo_90"].~TArrayD();
    
        // tmps_[abab_oovv_109](i,j,a,b) = 0.25 V_blks_[abab_oovv](k,l,c,d) * u2[abab_vvoo](c,d,i,j) * t2[abab_vvoo](a,b,k,l) // flops: o2v2 = o4v2 o4v2 | mem: o2v2 = o4v0 o2v2
        tmps_["abab_oovv_109"]("i,j,a,b")  = 0.25 * V_blks_["abab_oovv"]("k,l,c,d") * u2["abab_vvoo"]("c,d,i,j") * t2["abab_vvoo"]("a,b,k,l");
    
        // ru2_abab += +0.25 <l,k||c,d>_abab u2_abab(c,d,i,j) t2_abab(a,b,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["abab_oovv_109"]("i,j,a,b");
    
        // ru2_abab += +0.25 <l,k||d,c>_abab u2_abab(d,c,i,j) t2_abab(a,b,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["abab_oovv_109"]("i,j,a,b");
    
        // ru2_abab += +0.25 <k,l||c,d>_abab u2_abab(c,d,i,j) t2_abab(a,b,k,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["abab_oovv_109"]("i,j,a,b");
    
        // ru2_abab += +0.25 <k,l||d,c>_abab u2_abab(d,c,i,j) t2_abab(a,b,k,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["abab_oovv_109"]("i,j,a,b");
        tmps_["abab_oovv_109"].~TArrayD();
    
        // tmps_[bbbb_vvoo_128](a,b,i,j) = 1.00 dp[bb_ov](k,c) * t2[bbbb_vvoo](c,b,i,j) * u1[bb_vo](a,k) // flops: o2v2 = o3v2 o3v2 | mem: o2v2 = o3v1 o2v2
        tmps_["bbbb_vvoo_128"]("a,b,i,j")  = dp["bb_ov"]("k,c") * t2["bbbb_vvoo"]("c,b,i,j") * u1["bb_vo"]("a,k");
    
        // rt2_bbbb += -1.00 P(a,b) d-_bb(k,c) u1_bb(b,k) t2_bbbb(c,a,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_128"]("b,a,i,j");
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_128"]("a,b,i,j");
    
        // ru2_bbbb += -1.00 P(a,b) d-_bb(k,c) u0 u1_bb(b,k) t2_bbbb(c,a,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= u0 * tmps_["bbbb_vvoo_128"]("b,a,i,j");
        ru2_bbbb("a,b,i,j") += u0 * tmps_["bbbb_vvoo_128"]("a,b,i,j");
        tmps_["bbbb_vvoo_128"].~TArrayD();
    
        // tmps_[aaaa_vvoo_129](a,b,i,j) = 1.00 dp[aa_ov](k,c) * t2[aaaa_vvoo](c,b,i,j) * u1[aa_vo](a,k) // flops: o2v2 = o3v2 o3v2 | mem: o2v2 = o3v1 o2v2
        tmps_["aaaa_vvoo_129"]("a,b,i,j")  = dp["aa_ov"]("k,c") * t2["aaaa_vvoo"]("c,b,i,j") * u1["aa_vo"]("a,k");
    
        // rt2_aaaa += -1.00 P(a,b) d-_aa(k,c) u1_aa(b,k) t2_aaaa(c,a,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_129"]("b,a,i,j");
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_129"]("a,b,i,j");
    
        // ru2_aaaa += -1.00 P(a,b) d-_aa(k,c) u0 u1_aa(b,k) t2_aaaa(c,a,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= u0 * tmps_["aaaa_vvoo_129"]("b,a,i,j");
        ru2_aaaa("a,b,i,j") += u0 * tmps_["aaaa_vvoo_129"]("a,b,i,j");
        tmps_["aaaa_vvoo_129"].~TArrayD();
    
        // tmps_[bbbb_ovvo_130](j,b,a,i) = 1.00 V_blks_[abab_oovo](k,l,c,j) * t2[abab_vvoo](c,a,k,i) * u1[bb_vo](b,l) // flops: o2v2 = o4v2 o3v2 | mem: o2v2 = o3v1 o2v2
        tmps_["bbbb_ovvo_130"]("j,b,a,i")  = V_blks_["abab_oovo"]("k,l,c,j") * t2["abab_vvoo"]("c,a,k,i") * u1["bb_vo"]("b,l");
    
        // ru2_bbbb += -1.00 P(i,j) P(a,b) <k,l||c,j>_abab u1_bb(b,l) t2_abab(c,a,k,i)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_ovvo_130"]("j,b,a,i");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_ovvo_130"]("i,b,a,j");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_ovvo_130"]("j,a,b,i");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_ovvo_130"]("i,a,b,j");
        tmps_["bbbb_ovvo_130"].~TArrayD();
    
        // tmps_[aaaa_ovov_133](j,b,i,a) = 1.00 V_blks_[aaaa_oovo](l,k,c,j) * t2[aaaa_vvoo](c,b,i,k) * u1[aa_vo](a,l) // flops: o2v2 = o4v2 o3v2 | mem: o2v2 = o3v1 o2v2
        tmps_["aaaa_ovov_133"]("j,b,i,a")  = V_blks_["aaaa_oovo"]("l,k,c,j") * t2["aaaa_vvoo"]("c,b,i,k") * u1["aa_vo"]("a,l");
    
        // ru2_aaaa += -1.00 P(i,j) P(a,b) <l,k||c,j>_aaaa u1_aa(b,l) t2_aaaa(c,a,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_ovov_133"]("j,a,i,b");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_ovov_133"]("i,a,j,b");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_ovov_133"]("j,b,i,a");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_ovov_133"]("i,b,j,a");
        tmps_["aaaa_ovov_133"].~TArrayD();
    
        // tmps_[bbbb_ovov_136](j,a,i,b) = 1.00 V_blks_[bbbb_oovo](l,k,c,j) * t2[bbbb_vvoo](c,a,i,k) * u1[bb_vo](b,l) // flops: o2v2 = o4v2 o3v2 | mem: o2v2 = o3v1 o2v2
        tmps_["bbbb_ovov_136"]("j,a,i,b")  = V_blks_["bbbb_oovo"]("l,k,c,j") * t2["bbbb_vvoo"]("c,a,i,k") * u1["bb_vo"]("b,l");
    
        // ru2_bbbb += -1.00 P(i,j) P(a,b) <l,k||c,j>_bbbb u1_bb(b,l) t2_bbbb(c,a,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_ovov_136"]("j,a,i,b");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_ovov_136"]("i,a,j,b");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_ovov_136"]("j,b,i,a");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_ovov_136"]("i,b,j,a");
        tmps_["bbbb_ovov_136"].~TArrayD();
    
        // tmps_[aaaa_ovov_138](j,b,i,a) = 1.00 V_blks_[abba_oovo](l,k,c,j) * t2[abab_vvoo](b,c,i,k) * u1[aa_vo](a,l) // flops: o2v2 = o4v2 o3v2 | mem: o2v2 = o3v1 o2v2
        tmps_["aaaa_ovov_138"]("j,b,i,a")  = V_blks_["abba_oovo"]("l,k,c,j") * t2["abab_vvoo"]("b,c,i,k") * u1["aa_vo"]("a,l");
    
        // ru2_aaaa += -1.00 P(i,j) P(a,b) <l,k||j,c>_abab u1_aa(b,l) t2_abab(a,c,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_ovov_138"]("j,a,i,b");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_ovov_138"]("i,a,j,b");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_ovov_138"]("j,b,i,a");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_ovov_138"]("i,b,j,a");
        tmps_["aaaa_ovov_138"].~TArrayD();
    
        // tmps_[abab_vvoo_140](a,b,i,j) = 1.00 dp[aa_ov](k,c) * t2[abab_vvoo](c,b,i,j) * u1[aa_vo](a,k) // flops: o2v2 = o3v2 o3v2 | mem: o2v2 = o3v1 o2v2
        tmps_["abab_vvoo_140"]("a,b,i,j")  = dp["aa_ov"]("k,c") * t2["abab_vvoo"]("c,b,i,j") * u1["aa_vo"]("a,k");
    
        // rt2_abab += +1.00 d-_aa(k,c) u1_aa(a,k) t2_abab(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["abab_vvoo_140"]("a,b,i,j");
    
        // ru2_abab += +1.00 d-_aa(k,c) u0 u1_aa(a,k) t2_abab(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += u0 * tmps_["abab_vvoo_140"]("a,b,i,j");
        tmps_["abab_vvoo_140"].~TArrayD();
    
        // tmps_[baab_vvoo_141](b,a,i,j) = 1.00 dp[bb_ov](k,c) * t2[abab_vvoo](a,c,i,j) * u1[bb_vo](b,k) // flops: o2v2 = o3v2 o3v2 | mem: o2v2 = o3v1 o2v2
        tmps_["baab_vvoo_141"]("b,a,i,j")  = dp["bb_ov"]("k,c") * t2["abab_vvoo"]("a,c,i,j") * u1["bb_vo"]("b,k");
    
        // rt2_abab += +1.00 d-_bb(k,c) u1_bb(b,k) t2_abab(a,c,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["baab_vvoo_141"]("b,a,i,j");
    
        // ru2_abab += +1.00 d-_bb(k,c) u0 u1_bb(b,k) t2_abab(a,c,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += u0 * tmps_["baab_vvoo_141"]("b,a,i,j");
        tmps_["baab_vvoo_141"].~TArrayD();
    
        // tmps_[aabb_ovvo_142](i,a,b,j) = 1.00 dp[aa_oo](k,i) * u2[abab_vvoo](a,b,k,j) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["aabb_ovvo_142"]("i,a,b,j")  = dp["aa_oo"]("k,i") * u2["abab_vvoo"]("a,b,k,j");
    
        // rt2_abab += +1.00 d-_aa(k,i) u2_abab(a,b,k,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["aabb_ovvo_142"]("i,a,b,j");
    
        // ru2_abab += +1.00 d-_aa(k,i) u0 u2_abab(a,b,k,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["aabb_ovvo_142"]("i,a,b,j") * u0;
        tmps_["aabb_ovvo_142"].~TArrayD();
    
        // tmps_[aabb_ovvo_143](i,a,b,j) = 1.00 dp[aa_oo](k,i) * t2[abab_vvoo](a,b,k,j) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["aabb_ovvo_143"]("i,a,b,j")  = dp["aa_oo"]("k,i") * t2["abab_vvoo"]("a,b,k,j");
    
        // rt2_abab += +1.00 d-_aa(k,i) u0 t2_abab(a,b,k,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["aabb_ovvo_143"]("i,a,b,j") * u0;
    
        // ru2_abab += +1.00 d+_aa(k,i) t2_abab(a,b,k,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["aabb_ovvo_143"]("i,a,b,j");
        tmps_["aabb_ovvo_143"].~TArrayD();
    
        // tmps_[baba_ovvo_144](j,a,b,i) = 1.00 dp[bb_oo](k,j) * t2[abab_vvoo](a,b,i,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["baba_ovvo_144"]("j,a,b,i")  = dp["bb_oo"]("k,j") * t2["abab_vvoo"]("a,b,i,k");
    
        // rt2_abab += +1.00 d-_bb(k,j) u0 t2_abab(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["baba_ovvo_144"]("j,a,b,i") * u0;
    
        // ru2_abab += +1.00 d+_bb(k,j) t2_abab(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["baba_ovvo_144"]("j,a,b,i");
        tmps_["baba_ovvo_144"].~TArrayD();
    
        // tmps_[bbbb_ovvo_145](j,a,b,i) = 1.00 F_blks_[bb_oo](k,j) * u2[bbbb_vvoo](a,b,i,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["bbbb_ovvo_145"]("j,a,b,i")  = F_blks_["bb_oo"]("k,j") * u2["bbbb_vvoo"]("a,b,i,k");
    
        // ru2_bbbb += -1.00 P(i,j) f_bb(k,j) u2_bbbb(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_ovvo_145"]("j,a,b,i");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_ovvo_145"]("i,a,b,j");
        tmps_["bbbb_ovvo_145"].~TArrayD();
    
        // tmps_[bbbb_vvoo_146](a,b,i,j) = 0.50 V_blks_[bbbb_vovv](a,k,c,d) * t2[bbbb_vvoo](c,d,i,j) * u1[bb_vo](b,k) // flops: o2v2 = o3v3 o3v2 | mem: o2v2 = o3v1 o2v2
        tmps_["bbbb_vvoo_146"]("a,b,i,j")  = 0.50 * V_blks_["bbbb_vovv"]("a,k,c,d") * t2["bbbb_vvoo"]("c,d,i,j") * u1["bb_vo"]("b,k");
    
        // ru2_bbbb += +0.50 P(a,b) <k,a||c,d>_bbbb u1_bb(b,k) t2_bbbb(c,d,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_146"]("a,b,i,j");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_146"]("b,a,i,j");
        tmps_["bbbb_vvoo_146"].~TArrayD();
    
        // tmps_[bbbb_voov_147](a,i,j,b) = 1.00 V_blks_[bbbb_vooo](a,k,i,j) * u1[bb_vo](b,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["bbbb_voov_147"]("a,i,j,b")  = V_blks_["bbbb_vooo"]("a,k,i,j") * u1["bb_vo"]("b,k");
    
        // ru2_bbbb += +1.00 P(a,b) <k,a||i,j>_bbbb u1_bb(b,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_voov_147"]("a,i,j,b");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_voov_147"]("b,i,j,a");
        tmps_["bbbb_voov_147"].~TArrayD();
    
        // tmps_[bbbb_vvoo_150](b,a,i,j) = 1.00 F_blks_[bb_ov](k,c) * t2[bbbb_vvoo](c,a,i,j) * u1[bb_vo](b,k) // flops: o2v2 = o3v2 o3v2 | mem: o2v2 = o3v1 o2v2
        tmps_["bbbb_vvoo_150"]("b,a,i,j")  = F_blks_["bb_ov"]("k,c") * t2["bbbb_vvoo"]("c,a,i,j") * u1["bb_vo"]("b,k");
    
        // ru2_bbbb += +1.00 P(a,b) f_bb(k,c) u1_bb(b,k) t2_bbbb(c,a,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_150"]("b,a,i,j");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_150"]("a,b,i,j");
        tmps_["bbbb_vvoo_150"].~TArrayD();
    
        // tmps_[bbbb_vvoo_151](a,b,i,j) = 1.00 dp[bb_ov](k,c) * u2[bbbb_vvoo](c,b,i,j) * u1[bb_vo](a,k) // flops: o2v2 = o3v2 o3v2 | mem: o2v2 = o3v1 o2v2
        tmps_["bbbb_vvoo_151"]("a,b,i,j")  = dp["bb_ov"]("k,c") * u2["bbbb_vvoo"]("c,b,i,j") * u1["bb_vo"]("a,k");
    
        // ru2_bbbb += +2.00 P(a,b) d-_bb(k,c) u1_bb(a,k) u2_bbbb(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += 2.00 * tmps_["bbbb_vvoo_151"]("a,b,i,j");
        ru2_bbbb("a,b,i,j") -= 2.00 * tmps_["bbbb_vvoo_151"]("b,a,i,j");
        tmps_["bbbb_vvoo_151"].~TArrayD();
    
        // tmps_[baba_voov_152](b,i,j,a) = 0.50 V_blks_[baab_vovv](b,k,c,d) * t2[abab_vvoo](c,d,i,j) * u1[aa_vo](a,k) // flops: o2v2 = o3v3 o3v2 | mem: o2v2 = o3v1 o2v2
        tmps_["baba_voov_152"]("b,i,j,a")  = 0.50 * V_blks_["baab_vovv"]("b,k,c,d") * t2["abab_vvoo"]("c,d,i,j") * u1["aa_vo"]("a,k");
    
        // ru2_abab += -0.50 <k,b||c,d>_abab u1_aa(a,k) t2_abab(c,d,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["baba_voov_152"]("b,i,j,a");
    
        // ru2_abab += -0.50 <k,b||d,c>_abab u1_aa(a,k) t2_abab(d,c,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["baba_voov_152"]("b,i,j,a");
        tmps_["baba_voov_152"].~TArrayD();
    
        // tmps_[abab_vvoo_153](a,b,i,j) = 0.50 V_blks_[abab_vovv](a,k,c,d) * t2[abab_vvoo](c,d,i,j) * u1[bb_vo](b,k) // flops: o2v2 = o3v3 o3v2 | mem: o2v2 = o3v1 o2v2
        tmps_["abab_vvoo_153"]("a,b,i,j")  = 0.50 * V_blks_["abab_vovv"]("a,k,c,d") * t2["abab_vvoo"]("c,d,i,j") * u1["bb_vo"]("b,k");
    
        // ru2_abab += -0.50 <a,k||c,d>_abab u1_bb(b,k) t2_abab(c,d,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["abab_vvoo_153"]("a,b,i,j");
    
        // ru2_abab += -0.50 <a,k||d,c>_abab u1_bb(b,k) t2_abab(d,c,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["abab_vvoo_153"]("a,b,i,j");
        tmps_["abab_vvoo_153"].~TArrayD();
    
        // tmps_[aaaa_ovvo_158](j,a,b,i) = 1.00 F_blks_[aa_oo](k,j) * u2[aaaa_vvoo](a,b,i,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["aaaa_ovvo_158"]("j,a,b,i")  = F_blks_["aa_oo"]("k,j") * u2["aaaa_vvoo"]("a,b,i,k");
    
        // ru2_aaaa += -1.00 P(i,j) f_aa(k,j) u2_aaaa(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_ovvo_158"]("j,a,b,i");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_ovvo_158"]("i,a,b,j");
        tmps_["aaaa_ovvo_158"].~TArrayD();
    
        // tmps_[aaaa_vvoo_159](a,b,i,j) = 0.50 V_blks_[aaaa_vovv](a,k,c,d) * t2[aaaa_vvoo](c,d,i,j) * u1[aa_vo](b,k) // flops: o2v2 = o3v3 o3v2 | mem: o2v2 = o3v1 o2v2
        tmps_["aaaa_vvoo_159"]("a,b,i,j")  = 0.50 * V_blks_["aaaa_vovv"]("a,k,c,d") * t2["aaaa_vvoo"]("c,d,i,j") * u1["aa_vo"]("b,k");
    
        // ru2_aaaa += +0.50 P(a,b) <k,a||c,d>_aaaa u1_aa(b,k) t2_aaaa(c,d,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_159"]("a,b,i,j");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_159"]("b,a,i,j");
        tmps_["aaaa_vvoo_159"].~TArrayD();
    
        // tmps_[aaaa_voov_160](b,i,j,a) = 1.00 dp[aa_ov](k,c) * u2[aaaa_vvoo](c,b,i,j) * u1[aa_vo](a,k) // flops: o2v2 = o3v2 o3v2 | mem: o2v2 = o3v1 o2v2
        tmps_["aaaa_voov_160"]("b,i,j,a")  = dp["aa_ov"]("k,c") * u2["aaaa_vvoo"]("c,b,i,j") * u1["aa_vo"]("a,k");
    
        // ru2_aaaa += +2.00 P(a,b) d-_aa(k,c) u1_aa(a,k) u2_aaaa(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += 2.00 * tmps_["aaaa_voov_160"]("b,i,j,a");
        ru2_aaaa("a,b,i,j") -= 2.00 * tmps_["aaaa_voov_160"]("a,i,j,b");
        tmps_["aaaa_voov_160"].~TArrayD();
    
        // tmps_[aaaa_voov_163](b,i,j,a) = 1.00 F_blks_[aa_ov](k,c) * t2[aaaa_vvoo](c,b,i,j) * u1[aa_vo](a,k) // flops: o2v2 = o3v2 o3v2 | mem: o2v2 = o3v1 o2v2
        tmps_["aaaa_voov_163"]("b,i,j,a")  = F_blks_["aa_ov"]("k,c") * t2["aaaa_vvoo"]("c,b,i,j") * u1["aa_vo"]("a,k");
    
        // ru2_aaaa += +1.00 P(a,b) f_aa(k,c) u1_aa(b,k) t2_aaaa(c,a,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_voov_163"]("a,i,j,b");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_voov_163"]("b,i,j,a");
        tmps_["aaaa_voov_163"].~TArrayD();
    
        // tmps_[bbbb_ovvo_164](j,a,b,i) = 1.00 F_blks_[bb_oo](k,j) * t2[bbbb_vvoo](a,b,i,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["bbbb_ovvo_164"]("j,a,b,i")  = F_blks_["bb_oo"]("k,j") * t2["bbbb_vvoo"]("a,b,i,k");
    
        // rt2_bbbb += -1.00 P(i,j) f_bb(k,j) t2_bbbb(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_ovvo_164"]("j,a,b,i");
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_ovvo_164"]("i,a,b,j");
        tmps_["bbbb_ovvo_164"].~TArrayD();
    
        // tmps_[aaaa_ovvo_168](j,a,b,i) = 1.00 F_blks_[aa_oo](k,j) * t2[aaaa_vvoo](a,b,i,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["aaaa_ovvo_168"]("j,a,b,i")  = F_blks_["aa_oo"]("k,j") * t2["aaaa_vvoo"]("a,b,i,k");
    
        // rt2_aaaa += -1.00 P(i,j) f_aa(k,j) t2_aaaa(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_ovvo_168"]("j,a,b,i");
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_ovvo_168"]("i,a,b,j");
        tmps_["aaaa_ovvo_168"].~TArrayD();
    
        // tmps_[aaaa_voov_169](a,i,j,b) = 1.00 V_blks_[aaaa_vooo](a,k,i,j) * u1[aa_vo](b,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["aaaa_voov_169"]("a,i,j,b")  = V_blks_["aaaa_vooo"]("a,k,i,j") * u1["aa_vo"]("b,k");
    
        // ru2_aaaa += +1.00 P(a,b) <k,a||i,j>_aaaa u1_aa(b,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_voov_169"]("a,i,j,b");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_voov_169"]("b,i,j,a");
        tmps_["aaaa_voov_169"].~TArrayD();
    
        // tmps_[bb_vv_175](a,c) = 1.00 V_blks_[baab_vovv](a,k,d,c) * u1[aa_vo](d,k) // flops: o0v2 = o1v3 | mem: o0v2 = o0v2
        tmps_["bb_vv_175"]("a,c")  = V_blks_["baab_vovv"]("a,k,d,c") * u1["aa_vo"]("d,k");
    
        // ru2_abab += +1.00 <k,b||d,c>_abab u1_aa(d,k) t2_abab(a,c,i,j)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["bb_vv_175"]("b,c") * t2["abab_vvoo"]("a,c,i,j");
    
        // tmps_[bbbb_voov_196](a,i,j,b) = 1.00 t2[bbbb_vvoo](c,a,i,j) * V_blks_[baab_vovv](b,k,d,c) * u1[aa_vo](d,k) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["bbbb_voov_196"]("a,i,j,b")  = t2["bbbb_vvoo"]("c,a,i,j") * tmps_["bb_vv_175"]("b,c");
        tmps_["bb_vv_175"].~TArrayD();
    
        // ru2_bbbb += +1.00 P(a,b) <k,a||d,c>_abab u1_aa(d,k) t2_bbbb(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_voov_196"]("b,i,j,a");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_voov_196"]("a,i,j,b");
        tmps_["bbbb_voov_196"].~TArrayD();
    
        // tmps_[bb_vv_176](b,c) = 1.00 V_blks_[bbbb_vovv](b,k,c,d) * u1[bb_vo](d,k) // flops: o0v2 = o1v3 | mem: o0v2 = o0v2
        tmps_["bb_vv_176"]("b,c")  = V_blks_["bbbb_vovv"]("b,k,c,d") * u1["bb_vo"]("d,k");
    
        // ru2_abab += -1.00 <k,b||c,d>_bbbb u1_bb(d,k) t2_abab(a,c,i,j)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["bb_vv_176"]("b,c") * t2["abab_vvoo"]("a,c,i,j");
    
        // tmps_[bbbb_voov_195](b,i,j,a) = 1.00 t2[bbbb_vvoo](c,b,i,j) * V_blks_[bbbb_vovv](a,k,c,d) * u1[bb_vo](d,k) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["bbbb_voov_195"]("b,i,j,a")  = t2["bbbb_vvoo"]("c,b,i,j") * tmps_["bb_vv_176"]("a,c");
        tmps_["bb_vv_176"].~TArrayD();
    
        // ru2_bbbb += -1.00 P(a,b) <k,a||c,d>_bbbb u1_bb(d,k) t2_bbbb(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_voov_195"]("b,i,j,a");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_voov_195"]("a,i,j,b");
        tmps_["bbbb_voov_195"].~TArrayD();
    
        // tmps_[aa_vv_177](a,c) = 1.00 V_blks_[aaaa_vovv](a,k,c,d) * u1[aa_vo](d,k) // flops: o0v2 = o1v3 | mem: o0v2 = o0v2
        tmps_["aa_vv_177"]("a,c")  = V_blks_["aaaa_vovv"]("a,k,c,d") * u1["aa_vo"]("d,k");
    
        // ru2_abab += -1.00 <k,a||c,d>_aaaa u1_aa(d,k) t2_abab(c,b,i,j)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["aa_vv_177"]("a,c") * t2["abab_vvoo"]("c,b,i,j");
    
        // tmps_[aaaa_vvoo_194](a,b,i,j) = 1.00 V_blks_[aaaa_vovv](a,k,c,d) * u1[aa_vo](d,k) * t2[aaaa_vvoo](c,b,i,j) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["aaaa_vvoo_194"]("a,b,i,j")  = tmps_["aa_vv_177"]("a,c") * t2["aaaa_vvoo"]("c,b,i,j");
        tmps_["aa_vv_177"].~TArrayD();
    
        // ru2_aaaa += -1.00 P(a,b) <k,a||c,d>_aaaa u1_aa(d,k) t2_aaaa(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_194"]("a,b,i,j");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_194"]("b,a,i,j");
        tmps_["aaaa_vvoo_194"].~TArrayD();
    
        // tmps_[aa_vv_178](b,c) = 1.00 V_blks_[abab_vovv](b,k,c,d) * u1[bb_vo](d,k) // flops: o0v2 = o1v3 | mem: o0v2 = o0v2
        tmps_["aa_vv_178"]("b,c")  = V_blks_["abab_vovv"]("b,k,c,d") * u1["bb_vo"]("d,k");
    
        // ru2_abab += +1.00 <a,k||c,d>_abab u1_bb(d,k) t2_abab(c,b,i,j)  // flops: o2v2 += o2v3 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["aa_vv_178"]("a,c") * t2["abab_vvoo"]("c,b,i,j");
    
        // tmps_[aaaa_voov_193](a,i,j,b) = 1.00 t2[aaaa_vvoo](c,a,i,j) * V_blks_[abab_vovv](b,k,c,d) * u1[bb_vo](d,k) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
        tmps_["aaaa_voov_193"]("a,i,j,b")  = t2["aaaa_vvoo"]("c,a,i,j") * tmps_["aa_vv_178"]("b,c");
        tmps_["aa_vv_178"].~TArrayD();
    
        // ru2_aaaa += +1.00 P(a,b) <a,k||c,d>_abab u1_bb(d,k) t2_aaaa(c,b,i,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_voov_193"]("b,i,j,a");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_voov_193"]("a,i,j,b");
        tmps_["aaaa_voov_193"].~TArrayD();
    
        // tmps_[bb_vo_179](a,i) = 1.00 dp[bb_ov](k,c) * u2[bbbb_vvoo](c,a,i,k) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
        tmps_["bb_vo_179"]("a,i")  = dp["bb_ov"]("k,c") * u2["bbbb_vvoo"]("c,a,i,k");
    
        // rt1_bb += +1.00 d-_bb(j,b) u2_bbbb(b,a,i,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_bb("a,i") += tmps_["bb_vo_179"]("a,i");
    
        // ru1_bb += +1.00 d-_bb(j,b) u0 u2_bbbb(b,a,i,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") += tmps_["bb_vo_179"]("a,i") * u0;
    
        // ru2_abab += +1.00 d-_bb(k,c) u1_aa(a,i) u2_bbbb(c,b,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["bb_vo_179"]("b,j") * u1["aa_vo"]("a,i");
    
        // ru2_bbbb += -1.00 P(i,j) P(a,b) d-_bb(k,c) u1_bb(a,j) u2_bbbb(c,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bb_vo_179"]("b,i") * u1["bb_vo"]("a,j");
        ru2_bbbb("a,b,i,j") += tmps_["bb_vo_179"]("b,j") * u1["bb_vo"]("a,i");
        ru2_bbbb("a,b,i,j") += tmps_["bb_vo_179"]("a,i") * u1["bb_vo"]("b,j");
        ru2_bbbb("a,b,i,j") -= tmps_["bb_vo_179"]("a,j") * u1["bb_vo"]("b,i");
        tmps_["bb_vo_179"].~TArrayD();
    
        // tmps_[aa_vo_180](b,i) = 1.00 dp[aa_ov](k,c) * u2[aaaa_vvoo](c,b,i,k) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
        tmps_["aa_vo_180"]("b,i")  = dp["aa_ov"]("k,c") * u2["aaaa_vvoo"]("c,b,i,k");
    
        // rt1_aa += +1.00 d-_aa(j,b) u2_aaaa(b,a,i,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_aa("a,i") += tmps_["aa_vo_180"]("a,i");
    
        // ru1_aa += +1.00 d-_aa(j,b) u0 u2_aaaa(b,a,i,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") += tmps_["aa_vo_180"]("a,i") * u0;
    
        // ru2_aaaa += -1.00 P(i,j) P(a,b) d-_aa(k,c) u1_aa(a,j) u2_aaaa(c,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= tmps_["aa_vo_180"]("b,i") * u1["aa_vo"]("a,j");
        ru2_aaaa("a,b,i,j") += tmps_["aa_vo_180"]("b,j") * u1["aa_vo"]("a,i");
        ru2_aaaa("a,b,i,j") += tmps_["aa_vo_180"]("a,i") * u1["aa_vo"]("b,j");
        ru2_aaaa("a,b,i,j") -= tmps_["aa_vo_180"]("a,j") * u1["aa_vo"]("b,i");
    
        // ru2_abab += +1.00 d-_aa(k,c) u1_bb(b,j) u2_aaaa(c,a,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["aa_vo_180"]("a,i") * u1["bb_vo"]("b,j");
        tmps_["aa_vo_180"].~TArrayD();
    
        // tmps_[bb_vo_181](a,j) = 1.00 dp[aa_ov](k,c) * t2[abab_vvoo](c,a,k,j) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
        tmps_["bb_vo_181"]("a,j")  = dp["aa_ov"]("k,c") * t2["abab_vvoo"]("c,a,k,j");
    
        // rt1_bb += -1.00 d-_aa(j,b) u0 t2_abab(b,a,j,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_bb("a,i") -= tmps_["bb_vo_181"]("a,i") * u0;
    
        // rt2_abab += -1.00 d-_aa(k,c) u1_aa(a,i) t2_abab(c,b,k,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= tmps_["bb_vo_181"]("b,j") * u1["aa_vo"]("a,i");
    
        // rt2_bbbb += +1.00 P(i,j) P(a,b) d-_aa(k,c) u1_bb(b,i) t2_abab(c,a,k,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") += tmps_["bb_vo_181"]("a,j") * u1["bb_vo"]("b,i");
        rt2_bbbb("a,b,i,j") -= tmps_["bb_vo_181"]("a,i") * u1["bb_vo"]("b,j");
        rt2_bbbb("a,b,i,j") -= tmps_["bb_vo_181"]("b,j") * u1["bb_vo"]("a,i");
        rt2_bbbb("a,b,i,j") += tmps_["bb_vo_181"]("b,i") * u1["bb_vo"]("a,j");
    
        // ru1_bb += -1.00 d+_aa(j,b) t2_abab(b,a,j,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= tmps_["bb_vo_181"]("a,i");
        tmps_["bb_vo_181"].~TArrayD();
    
        // tmps_[bb_vo_182](a,i) = 1.00 dp[bb_ov](k,c) * t2[bbbb_vvoo](c,a,i,k) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
        tmps_["bb_vo_182"]("a,i")  = dp["bb_ov"]("k,c") * t2["bbbb_vvoo"]("c,a,i,k");
    
        // rt1_bb += +1.00 d-_bb(j,b) u0 t2_bbbb(b,a,i,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_bb("a,i") += tmps_["bb_vo_182"]("a,i") * u0;
    
        // rt2_abab += +1.00 d-_bb(k,c) u1_aa(a,i) t2_bbbb(c,b,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["bb_vo_182"]("b,j") * u1["aa_vo"]("a,i");
    
        // rt2_bbbb += -1.00 P(i,j) P(a,b) d-_bb(k,c) u1_bb(b,i) t2_bbbb(c,a,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") -= tmps_["bb_vo_182"]("a,j") * u1["bb_vo"]("b,i");
        rt2_bbbb("a,b,i,j") += tmps_["bb_vo_182"]("a,i") * u1["bb_vo"]("b,j");
        rt2_bbbb("a,b,i,j") += tmps_["bb_vo_182"]("b,j") * u1["bb_vo"]("a,i");
        rt2_bbbb("a,b,i,j") -= tmps_["bb_vo_182"]("b,i") * u1["bb_vo"]("a,j");
    
        // ru1_bb += +1.00 d+_bb(j,b) t2_bbbb(b,a,i,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") += tmps_["bb_vo_182"]("a,i");
        tmps_["bb_vo_182"].~TArrayD();
    
        // tmps_[bb_vo_183](a,i) = 1.00 dp[aa_ov](k,c) * u2[abab_vvoo](c,a,k,i) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
        tmps_["bb_vo_183"]("a,i")  = dp["aa_ov"]("k,c") * u2["abab_vvoo"]("c,a,k,i");
    
        // rt1_bb += -1.00 d-_aa(j,b) u2_abab(b,a,j,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_bb("a,i") -= tmps_["bb_vo_183"]("a,i");
    
        // ru1_bb += -1.00 d-_aa(j,b) u0 u2_abab(b,a,j,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= tmps_["bb_vo_183"]("a,i") * u0;
    
        // ru2_abab += -1.00 d-_aa(k,c) u1_aa(a,i) u2_abab(c,b,k,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["bb_vo_183"]("b,j") * u1["aa_vo"]("a,i");
    
        // ru2_bbbb += +1.00 P(i,j) P(a,b) d-_aa(k,c) u1_bb(a,j) u2_abab(c,b,k,i)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bb_vo_183"]("b,i") * u1["bb_vo"]("a,j");
        ru2_bbbb("a,b,i,j") -= tmps_["bb_vo_183"]("b,j") * u1["bb_vo"]("a,i");
        ru2_bbbb("a,b,i,j") -= tmps_["bb_vo_183"]("a,i") * u1["bb_vo"]("b,j");
        ru2_bbbb("a,b,i,j") += tmps_["bb_vo_183"]("a,j") * u1["bb_vo"]("b,i");
        tmps_["bb_vo_183"].~TArrayD();
    
        // tmps_[aa_vo_184](b,i) = 1.00 dp[aa_ov](k,c) * t2[aaaa_vvoo](c,b,i,k) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
        tmps_["aa_vo_184"]("b,i")  = dp["aa_ov"]("k,c") * t2["aaaa_vvoo"]("c,b,i,k");
    
        // rt1_aa += +1.00 d-_aa(j,b) u0 t2_aaaa(b,a,i,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_aa("a,i") += tmps_["aa_vo_184"]("a,i") * u0;
    
        // rt2_aaaa += -1.00 P(i,j) P(a,b) d-_aa(k,c) u1_aa(b,i) t2_aaaa(c,a,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") -= tmps_["aa_vo_184"]("a,j") * u1["aa_vo"]("b,i");
        rt2_aaaa("a,b,i,j") += tmps_["aa_vo_184"]("a,i") * u1["aa_vo"]("b,j");
        rt2_aaaa("a,b,i,j") += tmps_["aa_vo_184"]("b,j") * u1["aa_vo"]("a,i");
        rt2_aaaa("a,b,i,j") -= tmps_["aa_vo_184"]("b,i") * u1["aa_vo"]("a,j");
    
        // rt2_abab += +1.00 d-_aa(k,c) u1_bb(b,j) t2_aaaa(c,a,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["aa_vo_184"]("a,i") * u1["bb_vo"]("b,j");
    
        // ru1_aa += +1.00 d+_aa(j,b) t2_aaaa(b,a,i,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") += tmps_["aa_vo_184"]("a,i");
        tmps_["aa_vo_184"].~TArrayD();
    
        // tmps_[aa_vo_185](a,i) = 1.00 dp[bb_ov](j,b) * t2[abab_vvoo](a,b,i,j) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
        tmps_["aa_vo_185"]("a,i")  = dp["bb_ov"]("j,b") * t2["abab_vvoo"]("a,b,i,j");
    
        // rt1_aa += -1.00 d-_bb(j,b) u0 t2_abab(a,b,i,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_aa("a,i") -= tmps_["aa_vo_185"]("a,i") * u0;
    
        // rt2_aaaa += +1.00 P(i,j) P(a,b) d-_bb(k,c) u1_aa(b,i) t2_abab(a,c,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") += tmps_["aa_vo_185"]("a,j") * u1["aa_vo"]("b,i");
        rt2_aaaa("a,b,i,j") -= tmps_["aa_vo_185"]("a,i") * u1["aa_vo"]("b,j");
        rt2_aaaa("a,b,i,j") -= tmps_["aa_vo_185"]("b,j") * u1["aa_vo"]("a,i");
        rt2_aaaa("a,b,i,j") += tmps_["aa_vo_185"]("b,i") * u1["aa_vo"]("a,j");
    
        // rt2_abab += -1.00 d-_bb(k,c) u1_bb(b,j) t2_abab(a,c,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") -= tmps_["aa_vo_185"]("a,i") * u1["bb_vo"]("b,j");
    
        // ru1_aa += -1.00 d+_bb(j,b) t2_abab(a,b,i,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= tmps_["aa_vo_185"]("a,i");
        tmps_["aa_vo_185"].~TArrayD();
    
        // tmps_[aa_vo_186](a,i) = 1.00 dp[bb_ov](j,b) * u2[abab_vvoo](a,b,i,j) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
        tmps_["aa_vo_186"]("a,i")  = dp["bb_ov"]("j,b") * u2["abab_vvoo"]("a,b,i,j");
    
        // rt1_aa += -1.00 d-_bb(j,b) u2_abab(a,b,i,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_aa("a,i") -= tmps_["aa_vo_186"]("a,i");
    
        // ru1_aa += -1.00 d-_bb(j,b) u0 u2_abab(a,b,i,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= tmps_["aa_vo_186"]("a,i") * u0;
    
        // ru2_aaaa += +1.00 P(i,j) P(a,b) d-_bb(k,c) u1_aa(a,j) u2_abab(b,c,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aa_vo_186"]("b,i") * u1["aa_vo"]("a,j");
        ru2_aaaa("a,b,i,j") -= tmps_["aa_vo_186"]("b,j") * u1["aa_vo"]("a,i");
        ru2_aaaa("a,b,i,j") -= tmps_["aa_vo_186"]("a,i") * u1["aa_vo"]("b,j");
        ru2_aaaa("a,b,i,j") += tmps_["aa_vo_186"]("a,j") * u1["aa_vo"]("b,i");
    
        // ru2_abab += -1.00 d-_bb(k,c) u1_bb(b,j) u2_abab(a,c,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["aa_vo_186"]("a,i") * u1["bb_vo"]("b,j");
        tmps_["aa_vo_186"].~TArrayD();
    
        // tmps_[aa_ov_187](j,b) = 1.00 V_blks_[abab_oovv](j,k,b,c) * u1[bb_vo](c,k) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
        tmps_["aa_ov_187"]("j,b")  = V_blks_["abab_oovv"]("j,k,b,c") * u1["bb_vo"]("c,k");
    
        // ru1_aa += -1.00 <j,k||b,c>_abab u1_bb(c,k) t2_aaaa(b,a,i,j)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= tmps_["aa_ov_187"]("j,b") * t2["aaaa_vvoo"]("b,a,i,j");
    
        // ru1_bb += +1.00 <j,k||b,c>_abab u1_bb(c,k) t2_abab(b,a,j,i)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        ru1_bb("a,i") += tmps_["aa_ov_187"]("j,b") * t2["abab_vvoo"]("b,a,j,i");
        tmps_["aa_ov_187"].~TArrayD();
    
        // tmps_[bb_ov_188](j,b) = 1.00 V_blks_[bbbb_oovv](k,j,b,c) * u1[bb_vo](c,k) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
        tmps_["bb_ov_188"]("j,b")  = V_blks_["bbbb_oovv"]("k,j,b,c") * u1["bb_vo"]("c,k");
    
        // ru1_aa += -1.00 <k,j||b,c>_bbbb u1_bb(c,k) t2_abab(a,b,i,j)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= tmps_["bb_ov_188"]("j,b") * t2["abab_vvoo"]("a,b,i,j");
    
        // ru1_bb += +1.00 <k,j||b,c>_bbbb u1_bb(c,k) t2_bbbb(b,a,i,j)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        ru1_bb("a,i") += tmps_["bb_ov_188"]("j,b") * t2["bbbb_vvoo"]("b,a,i,j");
        tmps_["bb_ov_188"].~TArrayD();
    
        // tmps_[aaaa_ovvo_189](i,a,b,j) = 0.50 V_blks_[aaaa_oovo](l,k,c,i) * u1[aa_vo](c,j) * t2[aaaa_vvoo](a,b,l,k) // flops: o2v2 = o4v1 o4v2 | mem: o2v2 = o4v0 o2v2
        tmps_["aaaa_ovvo_189"]("i,a,b,j")  = 0.50 * V_blks_["aaaa_oovo"]("l,k,c,i") * u1["aa_vo"]("c,j") * t2["aaaa_vvoo"]("a,b,l,k");
    
        // ru2_aaaa += +0.50 P(i,j) <l,k||c,j>_aaaa u1_aa(c,i) t2_aaaa(a,b,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_ovvo_189"]("j,a,b,i");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_ovvo_189"]("i,a,b,j");
        tmps_["aaaa_ovvo_189"].~TArrayD();
    
        // tmps_[abab_oovv_190](i,j,a,b) = 0.50 V_blks_[abba_oovo](l,k,c,i) * u1[bb_vo](c,j) * t2[abab_vvoo](a,b,l,k) // flops: o2v2 = o4v1 o4v2 | mem: o2v2 = o4v0 o2v2
        tmps_["abab_oovv_190"]("i,j,a,b")  = 0.50 * V_blks_["abba_oovo"]("l,k,c,i") * u1["bb_vo"]("c,j") * t2["abab_vvoo"]("a,b,l,k");
    
        // ru2_abab += +0.50 <l,k||i,c>_abab u1_bb(c,j) t2_abab(a,b,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["abab_oovv_190"]("i,j,a,b");
    
        // ru2_abab += +0.50 <k,l||i,c>_abab u1_bb(c,j) t2_abab(a,b,k,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["abab_oovv_190"]("i,j,a,b");
        tmps_["abab_oovv_190"].~TArrayD();
    
        // tmps_[baba_ovvo_191](j,a,b,i) = 0.50 V_blks_[abab_oovo](l,k,c,j) * u1[aa_vo](c,i) * t2[abab_vvoo](a,b,l,k) // flops: o2v2 = o4v1 o4v2 | mem: o2v2 = o4v0 o2v2
        tmps_["baba_ovvo_191"]("j,a,b,i")  = 0.50 * V_blks_["abab_oovo"]("l,k,c,j") * u1["aa_vo"]("c,i") * t2["abab_vvoo"]("a,b,l,k");
    
        // ru2_abab += +0.50 <l,k||c,j>_abab u1_aa(c,i) t2_abab(a,b,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["baba_ovvo_191"]("j,a,b,i");
    
        // ru2_abab += +0.50 <k,l||c,j>_abab u1_aa(c,i) t2_abab(a,b,k,l)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["baba_ovvo_191"]("j,a,b,i");
        tmps_["baba_ovvo_191"].~TArrayD();
    
        // tmps_[bbbb_ovvo_192](i,a,b,j) = 0.50 V_blks_[bbbb_oovo](l,k,c,i) * u1[bb_vo](c,j) * t2[bbbb_vvoo](a,b,l,k) // flops: o2v2 = o4v1 o4v2 | mem: o2v2 = o4v0 o2v2
        tmps_["bbbb_ovvo_192"]("i,a,b,j")  = 0.50 * V_blks_["bbbb_oovo"]("l,k,c,i") * u1["bb_vo"]("c,j") * t2["bbbb_vvoo"]("a,b,l,k");
    
        // ru2_bbbb += +0.50 P(i,j) <l,k||c,j>_bbbb u1_bb(c,i) t2_bbbb(a,b,l,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_ovvo_192"]("j,a,b,i");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_ovvo_192"]("i,a,b,j");
        tmps_["bbbb_ovvo_192"].~TArrayD();
    
        // tmps_[aa_ov_197](j,b) = 1.00 V_blks_[aaaa_oovv](k,j,b,c) * u1[aa_vo](c,k) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
        tmps_["aa_ov_197"]("j,b")  = V_blks_["aaaa_oovv"]("k,j,b,c") * u1["aa_vo"]("c,k");
    
        // ru1_aa += +1.00 <k,j||b,c>_aaaa u1_aa(c,k) t2_aaaa(b,a,i,j)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        ru1_aa("a,i") += tmps_["aa_ov_197"]("j,b") * t2["aaaa_vvoo"]("b,a,i,j");
    
        // ru1_bb += -1.00 <k,j||b,c>_aaaa u1_aa(c,k) t2_abab(b,a,j,i)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= tmps_["aa_ov_197"]("j,b") * t2["abab_vvoo"]("b,a,j,i");
        tmps_["aa_ov_197"].~TArrayD();
    
        // tmps_[bb_ov_198](j,b) = 1.00 V_blks_[abab_oovv](k,j,c,b) * u1[aa_vo](c,k) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
        tmps_["bb_ov_198"]("j,b")  = V_blks_["abab_oovv"]("k,j,c,b") * u1["aa_vo"]("c,k");
    
        // ru1_aa += +1.00 <k,j||c,b>_abab u1_aa(c,k) t2_abab(a,b,i,j)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        ru1_aa("a,i") += tmps_["bb_ov_198"]("j,b") * t2["abab_vvoo"]("a,b,i,j");
    
        // ru1_bb += -1.00 <k,j||c,b>_abab u1_aa(c,k) t2_bbbb(b,a,i,j)  // flops: o1v1 += o2v2 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= tmps_["bb_ov_198"]("j,b") * t2["bbbb_vvoo"]("b,a,i,j");
        tmps_["bb_ov_198"].~TArrayD();
    
        // tmps_[bb_oo_199](k,i) = 1.00 V_blks_[abab_oovo](l,k,c,i) * u1[aa_vo](c,l) // flops: o2v0 = o3v1 | mem: o2v0 = o2v0
        tmps_["bb_oo_199"]("k,i")  = V_blks_["abab_oovo"]("l,k,c,i") * u1["aa_vo"]("c,l");
    
        // ru2_abab += -1.00 <l,k||c,j>_abab u1_aa(c,l) t2_abab(a,b,i,k)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["bb_oo_199"]("k,j") * t2["abab_vvoo"]("a,b,i,k");
    
        // tmps_[bbbb_vvoo_227](a,b,i,j) = 1.00 t2[bbbb_vvoo](a,b,i,k) * V_blks_[abab_oovo](l,k,c,j) * u1[aa_vo](c,l) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["bbbb_vvoo_227"]("a,b,i,j")  = t2["bbbb_vvoo"]("a,b,i,k") * tmps_["bb_oo_199"]("k,j");
        tmps_["bb_oo_199"].~TArrayD();
    
        // ru2_bbbb += -1.00 P(i,j) <l,k||c,j>_abab u1_aa(c,l) t2_bbbb(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_227"]("a,b,i,j");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_227"]("a,b,j,i");
        tmps_["bbbb_vvoo_227"].~TArrayD();
    
        // tmps_[aa_oo_200](k,i) = 1.00 V_blks_[aaaa_oovo](l,k,c,i) * u1[aa_vo](c,l) // flops: o2v0 = o3v1 | mem: o2v0 = o2v0
        tmps_["aa_oo_200"]("k,i")  = V_blks_["aaaa_oovo"]("l,k,c,i") * u1["aa_vo"]("c,l");
    
        // ru2_abab += -1.00 <l,k||c,i>_aaaa u1_aa(c,l) t2_abab(a,b,k,j)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["aa_oo_200"]("k,i") * t2["abab_vvoo"]("a,b,k,j");
    
        // tmps_[aaaa_vvoo_226](a,b,i,j) = 1.00 t2[aaaa_vvoo](a,b,i,k) * V_blks_[aaaa_oovo](l,k,c,j) * u1[aa_vo](c,l) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["aaaa_vvoo_226"]("a,b,i,j")  = t2["aaaa_vvoo"]("a,b,i,k") * tmps_["aa_oo_200"]("k,j");
        tmps_["aa_oo_200"].~TArrayD();
    
        // ru2_aaaa += -1.00 P(i,j) <l,k||c,j>_aaaa u1_aa(c,l) t2_aaaa(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_226"]("a,b,i,j");
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_226"]("a,b,j,i");
        tmps_["aaaa_vvoo_226"].~TArrayD();
    
        // tmps_[aa_oo_201](k,j) = 1.00 V_blks_[abba_oovo](k,l,c,j) * u1[bb_vo](c,l) // flops: o2v0 = o3v1 | mem: o2v0 = o2v0
        tmps_["aa_oo_201"]("k,j")  = V_blks_["abba_oovo"]("k,l,c,j") * u1["bb_vo"]("c,l");
    
        // ru2_abab += -1.00 <k,l||i,c>_abab u1_bb(c,l) t2_abab(a,b,k,j)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["aa_oo_201"]("k,i") * t2["abab_vvoo"]("a,b,k,j");
    
        // tmps_[aaaa_vvoo_225](a,b,i,j) = 1.00 t2[aaaa_vvoo](a,b,i,k) * V_blks_[abba_oovo](k,l,c,j) * u1[bb_vo](c,l) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["aaaa_vvoo_225"]("a,b,i,j")  = t2["aaaa_vvoo"]("a,b,i,k") * tmps_["aa_oo_201"]("k,j");
        tmps_["aa_oo_201"].~TArrayD();
    
        // ru2_aaaa += -1.00 P(i,j) <k,l||j,c>_abab u1_bb(c,l) t2_aaaa(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_225"]("a,b,i,j");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_225"]("a,b,j,i");
        tmps_["aaaa_vvoo_225"].~TArrayD();
    
        // tmps_[bb_oo_202](k,j) = 1.00 V_blks_[bbbb_oovo](l,k,c,j) * u1[bb_vo](c,l) // flops: o2v0 = o3v1 | mem: o2v0 = o2v0
        tmps_["bb_oo_202"]("k,j")  = V_blks_["bbbb_oovo"]("l,k,c,j") * u1["bb_vo"]("c,l");
    
        // ru2_abab += -1.00 <l,k||c,j>_bbbb u1_bb(c,l) t2_abab(a,b,i,k)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["bb_oo_202"]("k,j") * t2["abab_vvoo"]("a,b,i,k");
    
        // tmps_[bbbb_vvoo_224](a,b,j,i) = 1.00 t2[bbbb_vvoo](a,b,j,k) * V_blks_[bbbb_oovo](l,k,c,i) * u1[bb_vo](c,l) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["bbbb_vvoo_224"]("a,b,j,i")  = t2["bbbb_vvoo"]("a,b,j,k") * tmps_["bb_oo_202"]("k,i");
        tmps_["bb_oo_202"].~TArrayD();
    
        // ru2_bbbb += -1.00 P(i,j) <l,k||c,j>_bbbb u1_bb(c,l) t2_bbbb(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_224"]("a,b,i,j");
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_224"]("a,b,j,i");
        tmps_["bbbb_vvoo_224"].~TArrayD();
    
        // tmps_[aa_oo_203](j0,j1) = 0.50 Id_blks_[bb_oo](i0,i1) * V_blks_[abab_oooo](j0,i0,j1,i1) // flops: o2v0 = o4v0 | mem: o2v0 = o2v0
        tmps_["aa_oo_203"]("j0,j1")  = 0.50 * Id_blks_["bb_oo"]("i0,i1") * V_blks_["abab_oooo"]("j0,i0,j1,i1");
    
        // energy += -0.50 <j,i||j,i>_abab  // flops: o0v0 += o2v0 | mem: o0v0 += o0v0
        energy -= 1.00 * dot(Id_blks_["aa_oo"]("j0,j1"), tmps_["aa_oo_203"]("j0,j1"));
    
        // energy += -0.50 <i,j||i,j>_abab  // flops: o0v0 += o2v0 | mem: o0v0 += o0v0
        energy -= 1.00 * dot(Id_blks_["aa_oo"]("i0,i1"), tmps_["aa_oo_203"]("i0,i1"));
        tmps_["aa_oo_203"].~TArrayD();
    
        // tmps_[aa_vo_204](a,i) = 1.00 dp[aa_vv](a,b) * u1[aa_vo](b,i) // flops: o1v1 = o1v2 | mem: o1v1 = o1v1
        tmps_["aa_vo_204"]("a,i")  = dp["aa_vv"]("a,b") * u1["aa_vo"]("b,i");
    
        // rt1_aa += -1.00 d-_aa(a,b) u1_aa(b,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_aa("a,i") -= tmps_["aa_vo_204"]("a,i");
    
        // ru1_aa += -1.00 d-_aa(a,b) u0 u1_aa(b,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") -= tmps_["aa_vo_204"]("a,i") * u0;
    
        // ru2_aaaa += +1.00 P(i,j) P(a,b) d-_aa(a,c) u1_aa(c,j) u1_aa(b,i)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aa_vo_204"]("a,j") * u1["aa_vo"]("b,i");
        ru2_aaaa("a,b,i,j") -= tmps_["aa_vo_204"]("a,i") * u1["aa_vo"]("b,j");
        ru2_aaaa("a,b,i,j") -= tmps_["aa_vo_204"]("b,j") * u1["aa_vo"]("a,i");
        ru2_aaaa("a,b,i,j") += tmps_["aa_vo_204"]("b,i") * u1["aa_vo"]("a,j");
    
        // ru2_abab += -1.00 d-_aa(a,c) u1_aa(c,i) u1_bb(b,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["aa_vo_204"]("a,i") * u1["bb_vo"]("b,j");
        tmps_["aa_vo_204"].~TArrayD();
    
        // tmps_[bb_vo_205](a,i) = 1.00 dp[bb_vv](a,b) * u1[bb_vo](b,i) // flops: o1v1 = o1v2 | mem: o1v1 = o1v1
        tmps_["bb_vo_205"]("a,i")  = dp["bb_vv"]("a,b") * u1["bb_vo"]("b,i");
    
        // rt1_bb += -1.00 d-_bb(a,b) u1_bb(b,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_bb("a,i") -= tmps_["bb_vo_205"]("a,i");
    
        // ru1_bb += -1.00 d-_bb(a,b) u0 u1_bb(b,i)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") -= tmps_["bb_vo_205"]("a,i") * u0;
    
        // ru2_abab += -1.00 d-_bb(b,c) u1_bb(c,j) u1_aa(a,i)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["bb_vo_205"]("b,j") * u1["aa_vo"]("a,i");
    
        // ru2_bbbb += +1.00 P(i,j) P(a,b) d-_bb(a,c) u1_bb(c,j) u1_bb(b,i)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bb_vo_205"]("a,j") * u1["bb_vo"]("b,i");
        ru2_bbbb("a,b,i,j") -= tmps_["bb_vo_205"]("a,i") * u1["bb_vo"]("b,j");
        ru2_bbbb("a,b,i,j") -= tmps_["bb_vo_205"]("b,j") * u1["bb_vo"]("a,i");
        ru2_bbbb("a,b,i,j") += tmps_["bb_vo_205"]("b,i") * u1["bb_vo"]("a,j");
        tmps_["bb_vo_205"].~TArrayD();
    
        // tmps_[bb_oo_208](k,j) = 1.00 dp[bb_ov](k,c) * u1[bb_vo](c,j) // flops: o2v0 = o2v1 | mem: o2v0 = o2v0
        tmps_["bb_oo_208"]("k,j")  = dp["bb_ov"]("k,c") * u1["bb_vo"]("c,j");
    
        // ru1_bb += +2.00 d-_bb(j,b) u1_bb(b,i) u1_bb(a,j)  // flops: o1v1 += o2v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") += 2.00 * tmps_["bb_oo_208"]("j,i") * u1["bb_vo"]("a,j");
    
        // ru2_abab += +2.00 d-_bb(k,c) u1_bb(c,j) u2_abab(a,b,i,k)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += 2.00 * tmps_["bb_oo_208"]("k,j") * u2["abab_vvoo"]("a,b,i,k");
    
        // tmps_[bbbb_vvoo_217](a,b,i,j) = 1.00 t2[bbbb_vvoo](a,b,i,k) * dp[bb_ov](k,c) * u1[bb_vo](c,j) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["bbbb_vvoo_217"]("a,b,i,j")  = t2["bbbb_vvoo"]("a,b,i,k") * tmps_["bb_oo_208"]("k,j");
    
        // rt2_bbbb += -1.00 P(i,j) d-_bb(k,c) u1_bb(c,i) t2_bbbb(a,b,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_217"]("a,b,j,i");
        rt2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_217"]("a,b,i,j");
    
        // ru2_bbbb += -1.00 P(i,j) d-_bb(k,c) u0 u1_bb(c,i) t2_bbbb(a,b,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= u0 * tmps_["bbbb_vvoo_217"]("a,b,j,i");
        ru2_bbbb("a,b,i,j") += u0 * tmps_["bbbb_vvoo_217"]("a,b,i,j");
        tmps_["bbbb_vvoo_217"].~TArrayD();
    
        // tmps_[abab_vvoo_219](a,b,i,j) = 1.00 t2[abab_vvoo](a,b,i,k) * dp[bb_ov](k,c) * u1[bb_vo](c,j) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["abab_vvoo_219"]("a,b,i,j")  = t2["abab_vvoo"]("a,b,i,k") * tmps_["bb_oo_208"]("k,j");
    
        // rt2_abab += +1.00 d-_bb(k,c) u1_bb(c,j) t2_abab(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["abab_vvoo_219"]("a,b,i,j");
    
        // ru2_abab += +1.00 d-_bb(k,c) u0 u1_bb(c,j) t2_abab(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += u0 * tmps_["abab_vvoo_219"]("a,b,i,j");
        tmps_["abab_vvoo_219"].~TArrayD();
    
        // tmps_[bbbb_vvoo_223](a,b,i,j) = 1.00 u2[bbbb_vvoo](a,b,i,k) * dp[bb_ov](k,c) * u1[bb_vo](c,j) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["bbbb_vvoo_223"]("a,b,i,j")  = u2["bbbb_vvoo"]("a,b,i,k") * tmps_["bb_oo_208"]("k,j");
        tmps_["bb_oo_208"].~TArrayD();
    
        // ru2_bbbb += +2.00 P(i,j) d-_bb(k,c) u1_bb(c,j) u2_bbbb(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += 2.00 * tmps_["bbbb_vvoo_223"]("a,b,i,j");
        ru2_bbbb("a,b,i,j") -= 2.00 * tmps_["bbbb_vvoo_223"]("a,b,j,i");
        tmps_["bbbb_vvoo_223"].~TArrayD();
    
        // tmps_[aa_oo_209](k,j) = 1.00 dp[aa_ov](k,c) * u1[aa_vo](c,j) // flops: o2v0 = o2v1 | mem: o2v0 = o2v0
        tmps_["aa_oo_209"]("k,j")  = dp["aa_ov"]("k,c") * u1["aa_vo"]("c,j");
    
        // ru1_aa += +2.00 d-_aa(j,b) u1_aa(b,i) u1_aa(a,j)  // flops: o1v1 += o2v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") += 2.00 * tmps_["aa_oo_209"]("j,i") * u1["aa_vo"]("a,j");
    
        // ru2_abab += +2.00 d-_aa(k,c) u1_aa(c,i) u2_abab(a,b,k,j)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += 2.00 * tmps_["aa_oo_209"]("k,i") * u2["abab_vvoo"]("a,b,k,j");
    
        // tmps_[aaaa_vvoo_216](a,b,i,j) = 1.00 t2[aaaa_vvoo](a,b,i,k) * dp[aa_ov](k,c) * u1[aa_vo](c,j) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["aaaa_vvoo_216"]("a,b,i,j")  = t2["aaaa_vvoo"]("a,b,i,k") * tmps_["aa_oo_209"]("k,j");
    
        // rt2_aaaa += -1.00 P(i,j) d-_aa(k,c) u1_aa(c,i) t2_aaaa(a,b,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_216"]("a,b,j,i");
        rt2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_216"]("a,b,i,j");
    
        // ru2_aaaa += -1.00 P(i,j) d-_aa(k,c) u0 u1_aa(c,i) t2_aaaa(a,b,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= u0 * tmps_["aaaa_vvoo_216"]("a,b,j,i");
        ru2_aaaa("a,b,i,j") += u0 * tmps_["aaaa_vvoo_216"]("a,b,i,j");
        tmps_["aaaa_vvoo_216"].~TArrayD();
    
        // tmps_[abba_vvoo_218](a,b,j,i) = 1.00 t2[abab_vvoo](a,b,k,j) * dp[aa_ov](k,c) * u1[aa_vo](c,i) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["abba_vvoo_218"]("a,b,j,i")  = t2["abab_vvoo"]("a,b,k,j") * tmps_["aa_oo_209"]("k,i");
    
        // rt2_abab += +1.00 d-_aa(k,c) u1_aa(c,i) t2_abab(a,b,k,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        rt2_abab("a,b,i,j") += tmps_["abba_vvoo_218"]("a,b,j,i");
    
        // ru2_abab += +1.00 d-_aa(k,c) u0 u1_aa(c,i) t2_abab(a,b,k,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += u0 * tmps_["abba_vvoo_218"]("a,b,j,i");
        tmps_["abba_vvoo_218"].~TArrayD();
    
        // tmps_[aaaa_vvoo_222](a,b,i,j) = 1.00 u2[aaaa_vvoo](a,b,i,k) * dp[aa_ov](k,c) * u1[aa_vo](c,j) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["aaaa_vvoo_222"]("a,b,i,j")  = u2["aaaa_vvoo"]("a,b,i,k") * tmps_["aa_oo_209"]("k,j");
        tmps_["aa_oo_209"].~TArrayD();
    
        // ru2_aaaa += +2.00 P(i,j) d-_aa(k,c) u1_aa(c,j) u2_aaaa(a,b,i,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += 2.00 * tmps_["aaaa_vvoo_222"]("a,b,i,j");
        ru2_aaaa("a,b,i,j") -= 2.00 * tmps_["aaaa_vvoo_222"]("a,b,j,i");
        tmps_["aaaa_vvoo_222"].~TArrayD();
    
        // tmps_[bb_ov_210](i,a) = 1.00 dp[bb_oo](j,i) * u1[bb_vo](a,j) // flops: o1v1 = o2v1 | mem: o1v1 = o1v1
        tmps_["bb_ov_210"]("i,a")  = dp["bb_oo"]("j,i") * u1["bb_vo"]("a,j");
    
        // rt1_bb += +1.00 d-_bb(j,i) u1_bb(a,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_bb("a,i") += tmps_["bb_ov_210"]("i,a");
    
        // ru1_bb += +1.00 d-_bb(j,i) u0 u1_bb(a,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_bb("a,i") += tmps_["bb_ov_210"]("i,a") * u0;
    
        // ru2_abab += +1.00 d-_bb(k,j) u1_bb(b,k) u1_aa(a,i)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["bb_ov_210"]("j,b") * u1["aa_vo"]("a,i");
    
        // ru2_bbbb += -1.00 P(i,j) P(a,b) d-_bb(k,j) u1_bb(a,k) u1_bb(b,i)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") -= tmps_["bb_ov_210"]("j,a") * u1["bb_vo"]("b,i");
        ru2_bbbb("a,b,i,j") += tmps_["bb_ov_210"]("i,a") * u1["bb_vo"]("b,j");
        ru2_bbbb("a,b,i,j") += tmps_["bb_ov_210"]("j,b") * u1["bb_vo"]("a,i");
        ru2_bbbb("a,b,i,j") -= tmps_["bb_ov_210"]("i,b") * u1["bb_vo"]("a,j");
        tmps_["bb_ov_210"].~TArrayD();
    
        // tmps_[aa_ov_211](i,a) = 1.00 dp[aa_oo](j,i) * u1[aa_vo](a,j) // flops: o1v1 = o2v1 | mem: o1v1 = o1v1
        tmps_["aa_ov_211"]("i,a")  = dp["aa_oo"]("j,i") * u1["aa_vo"]("a,j");
    
        // rt1_aa += +1.00 d-_aa(j,i) u1_aa(a,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        rt1_aa("a,i") += tmps_["aa_ov_211"]("i,a");
    
        // ru1_aa += +1.00 d-_aa(j,i) u0 u1_aa(a,j)  // flops: o1v1 += o1v1 | mem: o1v1 += o1v1
        ru1_aa("a,i") += tmps_["aa_ov_211"]("i,a") * u0;
    
        // ru2_aaaa += -1.00 P(i,j) P(a,b) d-_aa(k,j) u1_aa(a,k) u1_aa(b,i)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") -= tmps_["aa_ov_211"]("j,a") * u1["aa_vo"]("b,i");
        ru2_aaaa("a,b,i,j") += tmps_["aa_ov_211"]("i,a") * u1["aa_vo"]("b,j");
        ru2_aaaa("a,b,i,j") += tmps_["aa_ov_211"]("j,b") * u1["aa_vo"]("a,i");
        ru2_aaaa("a,b,i,j") -= tmps_["aa_ov_211"]("i,b") * u1["aa_vo"]("a,j");
    
        // ru2_abab += +1.00 d-_aa(k,i) u1_aa(a,k) u1_bb(b,j)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") += tmps_["aa_ov_211"]("i,a") * u1["bb_vo"]("b,j");
        tmps_["aa_ov_211"].~TArrayD();
    
        // tmps_[bb_oo_212](k,j) = 1.00 F_blks_[bb_ov](k,c) * u1[bb_vo](c,j) // flops: o2v0 = o2v1 | mem: o2v0 = o2v0
        tmps_["bb_oo_212"]("k,j")  = F_blks_["bb_ov"]("k,c") * u1["bb_vo"]("c,j");
    
        // ru2_abab += -1.00 f_bb(k,c) u1_bb(c,j) t2_abab(a,b,i,k)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["bb_oo_212"]("k,j") * t2["abab_vvoo"]("a,b,i,k");
    
        // tmps_[bbbb_vvoo_221](a,b,j,i) = 1.00 t2[bbbb_vvoo](a,b,j,k) * F_blks_[bb_ov](k,c) * u1[bb_vo](c,i) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["bbbb_vvoo_221"]("a,b,j,i")  = t2["bbbb_vvoo"]("a,b,j,k") * tmps_["bb_oo_212"]("k,i");
        tmps_["bb_oo_212"].~TArrayD();
    
        // ru2_bbbb += +1.00 P(i,j) f_bb(k,c) u1_bb(c,i) t2_bbbb(a,b,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_bbbb("a,b,i,j") += tmps_["bbbb_vvoo_221"]("a,b,j,i");
        ru2_bbbb("a,b,i,j") -= tmps_["bbbb_vvoo_221"]("a,b,i,j");
        tmps_["bbbb_vvoo_221"].~TArrayD();
    
        // tmps_[aa_oo_213](k,i) = 1.00 F_blks_[aa_ov](k,c) * u1[aa_vo](c,i) // flops: o2v0 = o2v1 | mem: o2v0 = o2v0
        tmps_["aa_oo_213"]("k,i")  = F_blks_["aa_ov"]("k,c") * u1["aa_vo"]("c,i");
    
        // ru2_abab += -1.00 f_aa(k,c) u1_aa(c,i) t2_abab(a,b,k,j)  // flops: o2v2 += o3v2 | mem: o2v2 += o2v2
        ru2_abab("a,b,i,j") -= tmps_["aa_oo_213"]("k,i") * t2["abab_vvoo"]("a,b,k,j");
    
        // tmps_[aaaa_vvoo_220](a,b,j,i) = 1.00 t2[aaaa_vvoo](a,b,j,k) * F_blks_[aa_ov](k,c) * u1[aa_vo](c,i) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
        tmps_["aaaa_vvoo_220"]("a,b,j,i")  = t2["aaaa_vvoo"]("a,b,j,k") * tmps_["aa_oo_213"]("k,i");
        tmps_["aa_oo_213"].~TArrayD();
    
        // ru2_aaaa += +1.00 P(i,j) f_aa(k,c) u1_aa(c,i) t2_aaaa(a,b,j,k)  // flops: o2v2 += o2v2 | mem: o2v2 += o2v2
        ru2_aaaa("a,b,i,j") += tmps_["aaaa_vvoo_220"]("a,b,j,i");
        ru2_aaaa("a,b,i,j") -= tmps_["aaaa_vvoo_220"]("a,b,i,j");
        tmps_["aaaa_vvoo_220"].~TArrayD();
    }




/*

        Total Number of Terms: 2074
        Total Contractions: (last) 4521 -> (new) 4211

        Total FLOP scaling:
        ------------------
           Scaling :  initial |  reorder | optimize ||  init diff |   opt diff
          -------- : -------- | -------- | -------- || ---------- | ----------
              o3v4 :       60 |        0 |        0 ||        -60 |          0
              o4v3 :       78 |        0 |        0 ||        -78 |          0
          -------- : -------- | -------- | -------- || ---------- | ----------
              o2v4 :       86 |        8 |        6 ||        -80 |         -2
              o3v3 :      446 |      260 |      111 ||       -335 |       -149
              o4v2 :       98 |      254 |       75 ||        -23 |       -179
          -------- : -------- | -------- | -------- || ---------- | ----------
              o1v4 :       53 |        9 |        4 ||        -49 |         -5
              o2v3 :      537 |      303 |      113 ||       -424 |       -190
              o3v2 :      765 |     1063 |      463 ||       -302 |       -600
              o4v1 :       83 |      143 |       44 ||        -39 |        -99
          -------- : -------- | -------- | -------- || ---------- | ----------
              o1v3 :       68 |       48 |       12 ||        -56 |        -36
              o2v2 :     1462 |     1596 |     2549 ||       1087 |        953
              o3v1 :       68 |       48 |      133 ||         65 |         85
              o4v0 :        8 |        4 |       24 ||         16 |         20
          -------- : -------- | -------- | -------- || ---------- | ----------
              o1v2 :      128 |       44 |       19 ||       -109 |        -25
              o2v1 :      160 |      302 |       72 ||        -88 |       -230
          -------- : -------- | -------- | -------- || ---------- | ----------
              o0v2 :       16 |       12 |       17 ||          1 |          5
              o1v1 :      324 |      314 |      394 ||         70 |         80
              o2v0 :       32 |       60 |       43 ||         11 |        -17
          -------- : -------- | -------- | -------- || ---------- | ----------
              o0v0 :       49 |       53 |      132 ||         83 |         79
          -------- : -------- | -------- | -------- || ---------- | ----------
             Total :     4521 |     4521 |     4211 ||       -310 |       -310

*/


        double coherent_scalar = coupling_factor_z;
        if ( options_.get_bool("QED_USE_RELAXED_ORBITALS")) {
            coherent_scalar *= e_dip_z_;
        } else {
            coherent_scalar *= -nuc_dip_z_;
        }

        energy += coherent_scalar * cenergy;
        rt1_aa(idx_map_[2]) += coherent_scalar * crt1_aa(idx_map_[2]);
        rt1_bb(idx_map_[2]) += coherent_scalar * crt1_bb(idx_map_[2]);
        rt2_aaaa(idx_map_[4]) += coherent_scalar * crt2_aaaa(idx_map_[4]);
        rt2_abab(idx_map_[4]) += coherent_scalar * crt2_abab(idx_map_[4]);
        rt2_bbbb(idx_map_[4]) += coherent_scalar * crt2_bbbb(idx_map_[4]);
        if (include_u0_) ru0 = ru0 + coherent_scalar * cru0;
        if (include_u1_) {
            ru1_aa(idx_map_[2]) += coherent_scalar * cru1_aa(idx_map_[2]);
            ru1_bb(idx_map_[2]) += coherent_scalar * cru1_bb(idx_map_[2]);
        }
        if (include_u2_) {
            ru2_aaaa(idx_map_[4]) += coherent_scalar * cru2_aaaa(idx_map_[4]);
            ru2_abab(idx_map_[4]) += coherent_scalar * cru2_abab(idx_map_[4]);
            ru2_bbbb(idx_map_[4]) += coherent_scalar * cru2_bbbb(idx_map_[4]);
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
