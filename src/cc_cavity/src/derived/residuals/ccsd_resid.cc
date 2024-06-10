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

#include "cc_cavity/include/derived//ccsd.h"


using namespace std;
using namespace TA;
using namespace hilbert;

double CCSD::build_residuals() {
    if ( !has_t1_integrals_ ) transform_integrals(true);

        // get residual tensors
    TArrayD &rt1_aa_vo = residuals_["t1_aa"];
    TArrayD &rt1_bb_vo = residuals_["t1_bb"];
    TArrayD &rt2_aaaa_vvoo = residuals_["t2_aaaa"];
    TArrayD &rt2_abab_vvoo = residuals_["t2_abab"];
    TArrayD &rt2_bbbb_vvoo = residuals_["t2_bbbb"];

    // initialize to zero
    double energy = 0;
    zero_tiles(rt1_aa_vo);
    zero_tiles(rt1_bb_vo);
    zero_tiles(rt2_aaaa_vvoo);
    zero_tiles(rt2_abab_vvoo);
    zero_tiles(rt2_bbbb_vvoo);

    // extract 1-body amplitudes
    TA::TArrayD &t1_aa_vo = amplitudes_["t1_aa"];
    TA::TArrayD &t1_bb_vo = amplitudes_["t1_bb"];

    // extract 2-body amplitudes
    TA::TArrayD &t2_aaaa_vvoo = amplitudes_["t2_aaaa"];
    TA::TArrayD &t2_abab_vvoo = amplitudes_["t2_abab"];
    TA::TArrayD &t2_bbbb_vvoo = amplitudes_["t2_bbbb"];

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
        scalar0 = dot(F_blks_["aa_oo"]("o0,o1"), Id_blks_["aa_oo"]("o0,o1"));
        scalar1 = dot(V_blks_["abab_oovv"]("j,i,b,a"), t2_abab_vvoo("b,a,j,i"));
        scalar2 = dot(V_blks_["abab_oooo"]("o0,o1,o2,o3"), Id_blks_["abab_oooo"]("o0,o1,o2,o3"));
        scalar3 = dot(V_blks_["bbbb_oovv"]("j,i,a,b"), t2_bbbb_vvoo("a,b,j,i"));
        scalar4 = dot(F_blks_["bb_oo"]("o0,o1"), Id_blks_["bb_oo"]("o0,o1"));
        scalar5 = dot(V_blks_["aaaa_oooo"]("o0,o1,o2,o3"), Id_blks_["aaaa_oooo"]("o0,o1,o2,o3"));
        scalar6 = dot(V_blks_["aaaa_oovv"]("j,i,a,b"), t2_aaaa_vvoo("a,b,j,i"));
        scalar7 = dot(V_blks_["bbbb_oooo"]("o0,o1,o2,o3"), Id_blks_["bbbb_oooo"]("o0,o1,o2,o3"));
        auto tempArray = vector<TA::TArrayD>(12);



/// ****** p†q ****** ///



        {

            // tempArray[0] += 1.000000 V_blks_["abab_oovv"]("i,j,a,b") t2_bbbb_vvoo("b,f,n,j")
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[0]("a,f,i,n") = V_blks_["abab_oovv"]("i,j,a,b") * t2_bbbb_vvoo("b,f,n,j");

            // rt2_abab_vvoo += 1.000000 <i,j||a,b>_abab t2_aaaa(a,e,m,i) t2_bbbb(b,f,n,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += tempArray[0]("a,f,i,n") * t2_aaaa_vvoo("a,e,m,i");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) <i,j||a,b>_abab t2_abab(a,e,i,n) t2_bbbb(b,f,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = tempArray[0]("a,f,i,m") * t2_abab_vvoo("a,e,i,n");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) <j,i||b,a>_abab t2_bbbb(a,e,n,i) t2_abab(b,f,j,m)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += tempArray[0]("b,e,j,n") * t2_abab_vvoo("b,f,j,m");

            // tempArray[1] += 1.000000 t2_abab_vvoo("f,b,m,j") V_blks_["abab_oovv"]("i,j,a,b")
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[1]("f,a,m,i") = t2_abab_vvoo("f,b,m,j") * V_blks_["abab_oovv"]("i,j,a,b");

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) <j,i||b,a>_abab t2_abab(e,a,n,i) t2_aaaa(b,f,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = tempArray[1]("e,b,n,j") * t2_aaaa_vvoo("b,f,m,j");

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) <i,j||a,b>_abab t2_aaaa(a,e,n,i) t2_abab(f,b,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += tempArray[1]("f,a,m,i") * t2_aaaa_vvoo("a,e,n,i");

            // rt2_abab_vvoo += 1.000000 <j,i||b,a>_abab t2_abab(e,a,m,i) t2_abab(b,f,j,n)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += tempArray[1]("e,b,m,j") * t2_abab_vvoo("b,f,j,n");

            // tempArray[2] += 1.000000 t2_abab_vvoo("b,f,j,n") V_blks_["aaaa_oovv"]("j,i,a,b")
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[2]("a,f,i,n") = t2_abab_vvoo("b,f,j,n") * V_blks_["aaaa_oovv"]("j,i,a,b");

            // rt2_abab_vvoo += 1.000000 <j,i||a,b>_aaaa t2_aaaa(a,e,m,i) t2_abab(b,f,j,n)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += tempArray[2]("a,f,i,n") * t2_aaaa_vvoo("a,e,m,i");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) <j,i||a,b>_aaaa t2_abab(a,e,i,n) t2_abab(b,f,j,m)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += tempArray[2]("a,f,i,m") * t2_abab_vvoo("a,e,i,n");

            // tempArray[3] += 1.000000 V_blks_["bbbb_oovv"]("j,i,a,b") t2_bbbb_vvoo("b,f,n,j")
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempArray[3]("a,f,i,n") = V_blks_["bbbb_oovv"]("j,i,a,b") * t2_bbbb_vvoo("b,f,n,j");

            // rt2_abab_vvoo += 1.000000 <j,i||a,b>_bbbb t2_abab(e,a,m,i) t2_bbbb(b,f,n,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += tempArray[3]("a,f,i,n") * t2_abab_vvoo("e,a,m,i");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) <j,i||a,b>_bbbb t2_bbbb(a,e,n,i) t2_bbbb(b,f,m,j)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") += tempArray[3]("a,f,i,m") * t2_bbbb_vvoo("a,e,n,i");

            // tempArray[4] += 1.000000 V_blks_["abab_oovv"]("j,i,b,a") t2_abab_vvoo("b,f,j,i")
            // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
            tempArray[4]("a,f") = V_blks_["abab_oovv"]("j,i,b,a") * t2_abab_vvoo("b,f,j,i");

            // rt2_abab_vvoo += -0.500000 <j,i||b,a>_abab t2_abab(e,a,m,n) t2_abab(b,f,j,i)
            // +                -0.500000 <i,j||b,a>_abab t2_abab(e,a,m,n) t2_abab(b,f,i,j)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= tempArray[4]("a,f") * t2_abab_vvoo("e,a,m,n");

            // rt2_bbbb_vvoo += 0.500000 <j,i||b,a>_abab t2_bbbb(a,e,m,n) t2_abab(b,f,j,i)
            // +                0.500000 <i,j||b,a>_abab t2_bbbb(a,e,m,n) t2_abab(b,f,i,j)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_bbbb_vvoo("e,f,m,n") += tempArray[4]("a,f") * t2_bbbb_vvoo("a,e,m,n");

            // rt2_bbbb_vvoo += -0.500000 <j,i||a,b>_abab t2_abab(a,e,j,i) t2_bbbb(b,f,m,n)
            // +                -0.500000 <i,j||a,b>_abab t2_abab(a,e,i,j) t2_bbbb(b,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_bbbb_vvoo("e,f,m,n") -= tempArray[4]("b,e") * t2_bbbb_vvoo("b,f,m,n");

            // tempArray[5] += 1.000000 t2_abab_vvoo("e,a,j,i") V_blks_["abab_oovv"]("j,i,b,a")
            // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
            tempArray[5]("e,b") = t2_abab_vvoo("e,a,j,i") * V_blks_["abab_oovv"]("j,i,b,a");

            // rt2_aaaa_vvoo += 0.500000 <j,i||a,b>_abab t2_aaaa(a,e,m,n) t2_abab(f,b,j,i)
            // +                0.500000 <i,j||a,b>_abab t2_aaaa(a,e,m,n) t2_abab(f,b,i,j)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_aaaa_vvoo("e,f,m,n") += tempArray[5]("f,a") * t2_aaaa_vvoo("a,e,m,n");

            // rt2_aaaa_vvoo += -0.500000 <j,i||b,a>_abab t2_abab(e,a,j,i) t2_aaaa(b,f,m,n)
            // +                -0.500000 <i,j||b,a>_abab t2_abab(e,a,i,j) t2_aaaa(b,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_aaaa_vvoo("e,f,m,n") -= tempArray[5]("e,b") * t2_aaaa_vvoo("b,f,m,n");

            // rt2_abab_vvoo += -0.500000 <j,i||b,a>_abab t2_abab(e,a,j,i) t2_abab(b,f,m,n)
            // +                -0.500000 <i,j||b,a>_abab t2_abab(e,a,i,j) t2_abab(b,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= tempArray[5]("e,b") * t2_abab_vvoo("b,f,m,n");

            // tempArray[6] += 0.500000 V_blks_["bbbb_oovv"]("j,i,a,b") t2_bbbb_vvoo("b,f,j,i")
            // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
            tempArray[6]("a,f") = 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * t2_bbbb_vvoo("b,f,j,i");

            // rt2_abab_vvoo += 0.500000 <j,i||a,b>_bbbb t2_abab(e,a,m,n) t2_bbbb(b,f,j,i)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += tempArray[6]("a,f") * t2_abab_vvoo("e,a,m,n");

            // rt2_bbbb_vvoo += -0.500000 <j,i||a,b>_bbbb t2_bbbb(a,e,m,n) t2_bbbb(b,f,j,i)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_bbbb_vvoo("e,f,m,n") -= tempArray[6]("a,f") * t2_bbbb_vvoo("a,e,m,n");

            // tempArray[7] += 0.500000 t2_aaaa_vvoo("a,e,j,i") V_blks_["aaaa_oovv"]("j,i,a,b")
            // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2,
            tempArray[7]("e,b") = 0.500000 * t2_aaaa_vvoo("a,e,j,i") * V_blks_["aaaa_oovv"]("j,i,a,b");

            // rt2_aaaa_vvoo += -0.500000 <j,i||a,b>_aaaa t2_aaaa(a,e,j,i) t2_aaaa(b,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_aaaa_vvoo("e,f,m,n") -= tempArray[7]("e,b") * t2_aaaa_vvoo("b,f,m,n");

            // rt2_abab_vvoo += -0.500000 <j,i||a,b>_aaaa t2_aaaa(a,e,j,i) t2_abab(b,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= tempArray[7]("e,b") * t2_abab_vvoo("b,f,m,n");

            // tempArray[8] += 1.000000 t2_abab_vvoo("a,b,i,n") V_blks_["abab_oovv"]("i,j,a,b")
            // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
            tempArray[8]("n,j") = t2_abab_vvoo("a,b,i,n") * V_blks_["abab_oovv"]("i,j,a,b");

            // rt2_abab_vvoo += -0.500000 <i,j||a,b>_abab t2_abab(a,b,i,n) t2_abab(e,f,m,j)
            // +                -0.500000 <i,j||b,a>_abab t2_abab(b,a,i,n) t2_abab(e,f,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= tempArray[8]("n,j") * t2_abab_vvoo("e,f,m,j");

            // tempPerm_bbbb_vvoo += -0.500000 P(m,n) <i,j||a,b>_abab t2_abab(a,b,i,n) t2_bbbb(e,f,m,j)
            // +                -0.500000 P(m,n) <i,j||b,a>_abab t2_abab(b,a,i,n) t2_bbbb(e,f,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= tempArray[8]("n,j") * t2_bbbb_vvoo("e,f,m,j");

            // tempArray[9] += 0.500000 V_blks_["bbbb_oovv"]("j,i,a,b") t2_bbbb_vvoo("a,b,n,i")
            // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
            tempArray[9]("j,n") = 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * t2_bbbb_vvoo("a,b,n,i");

            // rt2_abab_vvoo += -0.500000 <j,i||a,b>_bbbb t2_bbbb(a,b,n,i) t2_abab(e,f,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= tempArray[9]("j,n") * t2_abab_vvoo("e,f,m,j");

            // tempPerm_bbbb_vvoo += -0.500000 P(m,n) <j,i||a,b>_bbbb t2_bbbb(a,b,n,i) t2_bbbb(e,f,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= tempArray[9]("j,n") * t2_bbbb_vvoo("e,f,m,j");

            // tempArray[10] += 1.000000 t2_abab_vvoo("b,a,n,i") V_blks_["abab_oovv"]("j,i,b,a")
            // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
            tempArray[10]("n,j") = t2_abab_vvoo("b,a,n,i") * V_blks_["abab_oovv"]("j,i,b,a");

            // tempPerm_aaaa_vvoo += -0.500000 P(m,n) <j,i||a,b>_abab t2_abab(a,b,n,i) t2_aaaa(e,f,m,j)
            // +                -0.500000 P(m,n) <j,i||b,a>_abab t2_abab(b,a,n,i) t2_aaaa(e,f,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= tempArray[10]("n,j") * t2_aaaa_vvoo("e,f,m,j");

            // rt2_abab_vvoo += -0.500000 <j,i||a,b>_abab t2_abab(a,b,m,i) t2_abab(e,f,j,n)
            // +                -0.500000 <j,i||b,a>_abab t2_abab(b,a,m,i) t2_abab(e,f,j,n)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= tempArray[10]("m,j") * t2_abab_vvoo("e,f,j,n");

            // tempArray[11] += 0.500000 t2_aaaa_vvoo("a,b,n,i") V_blks_["aaaa_oovv"]("j,i,a,b")
            // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2,
            tempArray[11]("n,j") = 0.500000 * t2_aaaa_vvoo("a,b,n,i") * V_blks_["aaaa_oovv"]("j,i,a,b");

            // tempPerm_aaaa_vvoo += -0.500000 P(m,n) <j,i||a,b>_aaaa t2_aaaa(a,b,n,i) t2_aaaa(e,f,m,j)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") -= tempArray[11]("n,j") * t2_aaaa_vvoo("e,f,m,j");

            // rt2_abab_vvoo += -0.500000 <j,i||a,b>_aaaa t2_aaaa(a,b,m,i) t2_abab(e,f,j,n)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= tempArray[11]("m,j") * t2_abab_vvoo("e,f,j,n");

            // energy += 1.000000 f_aa(i,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar0;

            // energy += 1.000000 f_bb(i,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += scalar4;

            // energy += -0.500000 <j,i||j,i>_aaaa
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar5;

            // energy += -0.500000 <j,i||j,i>_abab
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar2;

            // energy += -0.500000 <i,j||i,j>_abab
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar2;

            // energy += -0.500000 <j,i||j,i>_bbbb
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy -= 0.500000 * scalar7;

            // energy += 0.250000 <j,i||a,b>_aaaa t2_aaaa(a,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.250000 * scalar6;

            // energy += 0.250000 <j,i||a,b>_abab t2_abab(a,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.250000 * scalar1;

            // energy += 0.250000 <i,j||a,b>_abab t2_abab(a,b,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.250000 * scalar1;

            // energy += 0.250000 <j,i||b,a>_abab t2_abab(b,a,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.250000 * scalar1;

            // energy += 0.250000 <i,j||b,a>_abab t2_abab(b,a,i,j)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.250000 * scalar1;

            // energy += 0.250000 <j,i||a,b>_bbbb t2_bbbb(a,b,j,i)
            // flops: o0v0: 1 | mem: o0v0: 1,
            energy += 0.250000 * scalar3;

            // rt1_aa_vo += 0.500000 <e,i||a,b>_abab t2_abab(a,b,m,i)
            // +            0.500000 <e,i||b,a>_abab t2_abab(b,a,m,i)
            // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_aa_vo("e,m") += V_blks_["abab_vovv"]("e,i,a,b") * t2_abab_vvoo("a,b,m,i");

            // rt1_aa_vo += -0.500000 <i,e||a,b>_aaaa t2_aaaa(a,b,m,i)
            // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_aa_vo("e,m") += 0.500000 * V_blks_["aaaa_vovv"]("e,i,a,b") * t2_aaaa_vvoo("a,b,m,i");

            // rt1_aa_vo += -0.500000 <j,i||m,a>_abab t2_abab(e,a,j,i)
            // +            -0.500000 <i,j||m,a>_abab t2_abab(e,a,i,j)
            // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_aa_vo("e,m") += V_blks_["abba_oovo"]("j,i,a,m") * t2_abab_vvoo("e,a,j,i");

            // rt1_aa_vo += -0.500000 <j,i||a,m>_aaaa t2_aaaa(a,e,j,i)
            // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_aa_vo("e,m") -= 0.500000 * V_blks_["aaaa_oovo"]("j,i,a,m") * t2_aaaa_vvoo("a,e,j,i");

            // rt1_aa_vo += 1.000000 f_bb(i,a) t2_abab(e,a,m,i)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_aa_vo("e,m") += F_blks_["bb_ov"]("i,a") * t2_abab_vvoo("e,a,m,i");

            // rt1_aa_vo += -1.000000 f_aa(i,a) t2_aaaa(a,e,m,i)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_aa_vo("e,m") -= F_blks_["aa_ov"]("i,a") * t2_aaaa_vvoo("a,e,m,i");

            // rt1_aa_vo += 1.000000 f_aa(e,m)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rt1_aa_vo("e,m") += F_blks_["aa_vo"]("e,m");

            // rt1_bb_vo += -0.500000 <i,e||a,b>_bbbb t2_bbbb(a,b,m,i)
            // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_bb_vo("e,m") += 0.500000 * V_blks_["bbbb_vovv"]("e,i,a,b") * t2_bbbb_vvoo("a,b,m,i");

            // rt1_bb_vo += -0.500000 <j,i||a,m>_bbbb t2_bbbb(a,e,j,i)
            // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_bb_vo("e,m") -= 0.500000 * V_blks_["bbbb_oovo"]("j,i,a,m") * t2_bbbb_vvoo("a,e,j,i");

            // rt1_bb_vo += -0.500000 <j,i||a,m>_abab t2_abab(a,e,j,i)
            // +            -0.500000 <i,j||a,m>_abab t2_abab(a,e,i,j)
            // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_bb_vo("e,m") -= V_blks_["abab_oovo"]("j,i,a,m") * t2_abab_vvoo("a,e,j,i");

            // rt1_bb_vo += 0.500000 <i,e||a,b>_abab t2_abab(a,b,i,m)
            // +            0.500000 <i,e||b,a>_abab t2_abab(b,a,i,m)
            // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_bb_vo("e,m") -= V_blks_["baab_vovv"]("e,i,a,b") * t2_abab_vvoo("a,b,i,m");

            // rt1_bb_vo += -1.000000 f_bb(i,a) t2_bbbb(a,e,m,i)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_bb_vo("e,m") -= F_blks_["bb_ov"]("i,a") * t2_bbbb_vvoo("a,e,m,i");

            // rt1_bb_vo += 1.000000 f_aa(i,a) t2_abab(a,e,i,m)
            // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2,
            rt1_bb_vo("e,m") += F_blks_["aa_ov"]("i,a") * t2_abab_vvoo("a,e,i,m");

            // rt1_bb_vo += 1.000000 f_bb(e,m)
            // flops: o1v1: 1 | mem: o1v1: 1,
            rt1_bb_vo("e,m") += F_blks_["bb_vo"]("e,m");

            // rt2_aaaa_vvoo += 0.500000 <e,f||a,b>_aaaa t2_aaaa(a,b,m,n)
            // flops: o2v4: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_aaaa_vvoo("e,f,m,n") += 0.500000 * V_blks_["aaaa_vvvv"]("e,f,a,b") * t2_aaaa_vvoo("a,b,m,n");

            // rt2_aaaa_vvoo += 0.250000 <j,i||a,b>_aaaa t2_aaaa(a,b,m,n) t2_aaaa(e,f,j,i)
            // flops: o4v2: 2, o2v2: 1 | mem: o2v2: 2, o4v0: 1,
            rt2_aaaa_vvoo("e,f,m,n") += 0.250000 * V_blks_["aaaa_oovv"]("j,i,a,b") * t2_aaaa_vvoo("a,b,m,n") * t2_aaaa_vvoo("e,f,j,i");

            // rt2_aaaa_vvoo += 1.000000 <e,f||m,n>_aaaa
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += V_blks_["aaaa_vvoo"]("e,f,m,n");

            // rt2_aaaa_vvoo += -0.500000 <j,i||a,b>_aaaa t2_aaaa(a,e,m,n) t2_aaaa(b,f,j,i)
            // flops: o2v3: 2, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") -= 0.500000 * V_blks_["aaaa_oovv"]("j,i,a,b") * t2_aaaa_vvoo("b,f,j,i") * t2_aaaa_vvoo("a,e,m,n");

            // rt2_aaaa_vvoo += 0.500000 <j,i||m,n>_aaaa t2_aaaa(e,f,j,i)
            // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_aaaa_vvoo("e,f,m,n") += 0.500000 * V_blks_["aaaa_oooo"]("j,i,m,n") * t2_aaaa_vvoo("e,f,j,i");

            // rt2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo = TArrayD(world_, tempPerm_aaaa_vvoo.trange());
            tempPerm_aaaa_vvoo.fill(0.0); world_.gop.fence();


            // tempPerm_aaaa_vvoo += 1.000000 P(e,f) f_aa(e,a) t2_aaaa(a,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = F_blks_["aa_vv"]("e,a") * t2_aaaa_vvoo("a,f,m,n");

            // rt2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo = TArrayD(world_, tempPerm_aaaa_vvoo.trange());
            tempPerm_aaaa_vvoo.fill(0.0); world_.gop.fence();


            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) f_aa(i,n) t2_aaaa(e,f,m,i)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * F_blks_["aa_oo"]("i,n") * t2_aaaa_vvoo("e,f,m,i");

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) <j,i||a,b>_aaaa t2_aaaa(a,e,n,i) t2_aaaa(b,f,m,j)
            // flops: o3v3: 2, o2v2: 1 | mem: o2v2: 3,
            tempPerm_aaaa_vvoo("e,f,m,n") += V_blks_["aaaa_oovv"]("j,i,a,b") * t2_aaaa_vvoo("b,f,m,j") * t2_aaaa_vvoo("a,e,n,i");

            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) <j,i||a,b>_bbbb t2_abab(e,a,n,i) t2_abab(f,b,m,j)
            // flops: o3v3: 2, o2v2: 1 | mem: o2v2: 3,
            tempPerm_aaaa_vvoo("e,f,m,n") += V_blks_["bbbb_oovv"]("j,i,a,b") * t2_abab_vvoo("f,b,m,j") * t2_abab_vvoo("e,a,n,i");

            // rt2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo = TArrayD(world_, tempPerm_aaaa_vvoo.trange());
            tempPerm_aaaa_vvoo.fill(0.0); world_.gop.fence();


            // tempPerm_aaaa_vvoo += 1.000000 P(m,n) P(e,f) <i,e||a,n>_aaaa t2_aaaa(a,f,m,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") = -1.0 * V_blks_["aaaa_vovo"]("e,i,a,n") * t2_aaaa_vvoo("a,f,m,i");

            // tempPerm_aaaa_vvoo += -1.000000 P(m,n) P(e,f) <e,i||n,a>_abab t2_abab(f,a,m,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_aaaa_vvoo("e,f,m,n") += V_blks_["abba_vovo"]("e,i,a,n") * t2_abab_vvoo("f,a,m,i");

            // rt2_abab_vvoo += 0.250000 <j,i||a,b>_abab t2_abab(a,b,m,n) t2_abab(e,f,j,i)
            // +                0.250000 <i,j||a,b>_abab t2_abab(a,b,m,n) t2_abab(e,f,i,j)
            // +                0.250000 <j,i||b,a>_abab t2_abab(b,a,m,n) t2_abab(e,f,j,i)
            // +                0.250000 <i,j||b,a>_abab t2_abab(b,a,m,n) t2_abab(e,f,i,j)
            // flops: o4v2: 2, o2v2: 1 | mem: o2v2: 2, o4v0: 1,
            rt2_abab_vvoo("e,f,m,n") += V_blks_["abab_oovv"]("j,i,a,b") * t2_abab_vvoo("a,b,m,n") * t2_abab_vvoo("e,f,j,i");

            // rt2_abab_vvoo += -1.000000 f_aa(i,m) t2_abab(e,f,i,n)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= F_blks_["aa_oo"]("i,m") * t2_abab_vvoo("e,f,i,n");

            // rt2_abab_vvoo += 1.000000 <e,f||m,n>_abab
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_abab_vvoo("e,f,m,n") += V_blks_["abab_vvoo"]("e,f,m,n");

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

            // rt2_abab_vvoo += 0.500000 <e,f||a,b>_abab t2_abab(a,b,m,n)
            // +                0.500000 <e,f||b,a>_abab t2_abab(b,a,m,n)
            // flops: o2v4: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += V_blks_["abab_vvvv"]("e,f,a,b") * t2_abab_vvoo("a,b,m,n");

            // rt2_abab_vvoo += -1.000000 <e,i||m,a>_abab t2_bbbb(a,f,n,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += V_blks_["abba_vovo"]("e,i,a,m") * t2_bbbb_vvoo("a,f,n,i");

            // rt2_abab_vvoo += 1.000000 f_bb(f,a) t2_abab(e,a,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += F_blks_["bb_vv"]("f,a") * t2_abab_vvoo("e,a,m,n");

            // rt2_abab_vvoo += -1.000000 <i,f||m,a>_abab t2_abab(e,a,i,n)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= V_blks_["baba_vovo"]("f,i,a,m") * t2_abab_vvoo("e,a,i,n");

            // rt2_abab_vvoo += 1.000000 <i,j||b,a>_abab t2_abab(e,a,i,n) t2_abab(b,f,m,j)
            // flops: o3v3: 2, o2v2: 1 | mem: o2v2: 3,
            rt2_abab_vvoo("e,f,m,n") += V_blks_["abab_oovv"]("i,j,b,a") * t2_abab_vvoo("e,a,i,n") * t2_abab_vvoo("b,f,m,j");

            // rt2_abab_vvoo += -1.000000 <e,i||a,n>_abab t2_abab(a,f,m,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= V_blks_["abab_vovo"]("e,i,a,n") * t2_abab_vvoo("a,f,m,i");

            // rt2_abab_vvoo += 1.000000 <i,f||a,n>_bbbb t2_abab(e,a,m,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= V_blks_["bbbb_vovo"]("f,i,a,n") * t2_abab_vvoo("e,a,m,i");

            // rt2_abab_vvoo += 1.000000 f_aa(e,a) t2_abab(a,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") += F_blks_["aa_vv"]("e,a") * t2_abab_vvoo("a,f,m,n");

            // rt2_abab_vvoo += -1.000000 f_bb(i,n) t2_abab(e,f,m,i)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_abab_vvoo("e,f,m,n") -= F_blks_["bb_oo"]("i,n") * t2_abab_vvoo("e,f,m,i");

            // rt2_bbbb_vvoo += 1.000000 <e,f||m,n>_bbbb
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += V_blks_["bbbb_vvoo"]("e,f,m,n");

            // rt2_bbbb_vvoo += -0.500000 <j,i||a,b>_bbbb t2_bbbb(a,e,j,i) t2_bbbb(b,f,m,n)
            // flops: o2v3: 2, o2v2: 1 | mem: o2v2: 2, o0v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") -= 0.500000 * V_blks_["bbbb_oovv"]("j,i,a,b") * t2_bbbb_vvoo("a,e,j,i") * t2_bbbb_vvoo("b,f,m,n");

            // rt2_bbbb_vvoo += 0.500000 <e,f||a,b>_bbbb t2_bbbb(a,b,m,n)
            // flops: o2v4: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_bbbb_vvoo("e,f,m,n") += 0.500000 * V_blks_["bbbb_vvvv"]("e,f,a,b") * t2_bbbb_vvoo("a,b,m,n");

            // rt2_bbbb_vvoo += 0.250000 <j,i||a,b>_bbbb t2_bbbb(a,b,m,n) t2_bbbb(e,f,j,i)
            // flops: o4v2: 2, o2v2: 1 | mem: o2v2: 2, o4v0: 1,
            rt2_bbbb_vvoo("e,f,m,n") += 0.250000 * V_blks_["bbbb_oovv"]("j,i,a,b") * t2_bbbb_vvoo("a,b,m,n") * t2_bbbb_vvoo("e,f,j,i");

            // rt2_bbbb_vvoo += 0.500000 <j,i||m,n>_bbbb t2_bbbb(e,f,j,i)
            // flops: o4v2: 1, o2v2: 1 | mem: o2v2: 2,
            rt2_bbbb_vvoo("e,f,m,n") += 0.500000 * V_blks_["bbbb_oooo"]("j,i,m,n") * t2_bbbb_vvoo("e,f,j,i");

            // rt2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo = TArrayD(world_, tempPerm_bbbb_vvoo.trange());
            tempPerm_bbbb_vvoo.fill(0.0); world_.gop.fence();


            // tempPerm_bbbb_vvoo += 1.000000 P(e,f) f_bb(e,a) t2_bbbb(a,f,m,n)
            // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = F_blks_["bb_vv"]("e,a") * t2_bbbb_vvoo("a,f,m,n");

            // rt2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo = TArrayD(world_, tempPerm_bbbb_vvoo.trange());
            tempPerm_bbbb_vvoo.fill(0.0); world_.gop.fence();


            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) f_bb(i,n) t2_bbbb(e,f,m,i)
            // flops: o3v2: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = -1.0 * F_blks_["bb_oo"]("i,n") * t2_bbbb_vvoo("e,f,m,i");

            // rt2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo --> Apply Perm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo = TArrayD(world_, tempPerm_bbbb_vvoo.trange());
            tempPerm_bbbb_vvoo.fill(0.0); world_.gop.fence();


            // tempPerm_bbbb_vvoo += -1.000000 P(m,n) P(e,f) <i,e||a,n>_abab t2_abab(a,f,i,m)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") = V_blks_["baab_vovo"]("e,i,a,n") * t2_abab_vvoo("a,f,i,m");

            // tempPerm_bbbb_vvoo += 1.000000 P(m,n) P(e,f) <i,e||a,n>_bbbb t2_bbbb(a,f,m,i)
            // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2,
            tempPerm_bbbb_vvoo("e,f,m,n") -= V_blks_["bbbb_vovo"]("e,i,a,n") * t2_bbbb_vvoo("a,f,m,i");

            // rt2_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo tempPerm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("e,f,m,n");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("e,f,n,m");
            rt2_aaaa_vvoo("e,f,m,n") -= tempPerm_aaaa_vvoo("f,e,m,n");
            rt2_aaaa_vvoo("e,f,m,n") += tempPerm_aaaa_vvoo("f,e,n,m");

            // tempPerm_aaaa_vvoo += 1.000000 tempPerm_aaaa_vvoo tempPerm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_aaaa_vvoo = TArrayD(world_, tempPerm_aaaa_vvoo.trange());
            tempPerm_aaaa_vvoo.fill(0.0); world_.gop.fence();


            // rt2_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo tempPerm
            // flops: o2v2: 1 | mem: o2v2: 1,
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("e,f,m,n");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("e,f,n,m");
            rt2_bbbb_vvoo("e,f,m,n") -= tempPerm_bbbb_vvoo("f,e,m,n");
            rt2_bbbb_vvoo("e,f,m,n") += tempPerm_bbbb_vvoo("f,e,n,m");

            // tempPerm_bbbb_vvoo += 1.000000 tempPerm_bbbb_vvoo tempPerm
            // flops: o2v2: 1 | mem: o2v2: 1,
            tempPerm_bbbb_vvoo = TArrayD(world_, tempPerm_bbbb_vvoo.trange());
            tempPerm_bbbb_vvoo.fill(0.0); world_.gop.fence();

        }
        world_.gop.fence(); tempArray[0].~TArrayD();
        world_.gop.fence(); tempArray[1].~TArrayD();
        world_.gop.fence(); tempArray[2].~TArrayD();
        world_.gop.fence(); tempArray[3].~TArrayD();
        world_.gop.fence(); tempArray[4].~TArrayD();
        world_.gop.fence(); tempArray[5].~TArrayD();
        world_.gop.fence(); tempArray[6].~TArrayD();
        world_.gop.fence(); tempArray[7].~TArrayD();
        world_.gop.fence(); tempArray[8].~TArrayD();
        world_.gop.fence(); tempArray[9].~TArrayD();
        world_.gop.fence(); tempArray[10].~TArrayD();
        world_.gop.fence(); tempArray[11].~TArrayD();

        /*
            Total Number of Terms: 117
            Number of Flops: (old) 278 -> (new) 209

            Total FLOP scaling:
            ------------------
               Scaling :      new |      old |     diff
              -------- : -------- | -------- | --------
                  o2v4 :        3 |        4 |       -1
                  o3v3 :       30 |       36 |       -6
                  o4v2 :        9 |       16 |       -7
                  o2v3 :       26 |       46 |      -20
                  o3v2 :       20 |       34 |      -14
                  o2v2 :       87 |      112 |      -25
                  o0v2 :        4 |        0 |        4
                  o1v1 :       14 |       18 |       -4
                  o2v0 :        4 |        0 |        4
                  o0v0 :       12 |       12 |        0

            Total MEM scaling:
                  o2v2 :      150 |      166 |      -16
                  o4v0 :        3 |        6 |       -3
                  o0v2 :       10 |       18 |       -8
                  o1v1 :       26 |       34 |       -8
                  o2v0 :        8 |       12 |       -4
                  o0v0 :       12 |       12 |        0
         */
    }
    world_.gop.fence();
    return energy + enuc_;
}
