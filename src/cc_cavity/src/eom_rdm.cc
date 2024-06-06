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

#include "cc_cavity/include/eom_rdm.h"

namespace hilbert {
    EOM_RDM::EOM_RDM(const shared_ptr<EOM_Driver> &eom_driver, Options &options) :
            eom_driver_(eom_driver), options_(options) {}

    void EOM_RDM::compute_oscillators() {

        // Remove t1 amplitudes from integrals if they are present
        if (eom_driver_->cc_wfn_->has_t1_integrals_) eom_driver_->cc_wfn_->transform_integrals(false);

        // Compute transition dipole tensors
        properties_["OSCILLATOR STRENGTHS"] = makeTensor(world_, {M_, M_}, true),
        properties_["X TRANSITION DIPOLES"] = makeTensor(world_, {M_, M_}, true),
        properties_["Y TRANSITION DIPOLES"] = makeTensor(world_, {M_, M_}, true),
        properties_["Z TRANSITION DIPOLES"] = makeTensor(world_, {M_, M_}, true);

        TArrayMap &Dip_blks = eom_driver_->cc_wfn_->Dip_blks_;

        // Compute transition dipole tensors
        for (auto &rdm: RDM_blks_) {
            string rdm_blk = rdm.first;

            // check if this is a 1RDM
            if (rdm_blk.find("D1") == string::npos) continue;

            // strip off the D1_ prefix and keep the rest
            rdm_blk = rdm_blk.substr(3);

            properties_["X TRANSITION DIPOLES"]("I,F") += rdm.second("I,F,p,q") * Dip_blks["dx_" + rdm_blk]("p,q");
            properties_["Y TRANSITION DIPOLES"]("I,F") += rdm.second("I,F,p,q") * Dip_blks["dy_" + rdm_blk]("p,q");
            properties_["Z TRANSITION DIPOLES"]("I,F") += rdm.second("I,F,p,q") * Dip_blks["dz_" + rdm_blk]("p,q");
        }

        // Compute oscillator strengths
        properties_["OSCILLATOR STRENGTHS"]("I,F") =
                properties_["X TRANSITION DIPOLES"]("I,F") * properties_["X TRANSITION DIPOLES"]("F,I") +
                properties_["Y TRANSITION DIPOLES"]("I,F") * properties_["Y TRANSITION DIPOLES"]("F,I") +
                properties_["Z TRANSITION DIPOLES"]("I,F") * properties_["Z TRANSITION DIPOLES"]("F,I");

        // Scale by energy difference
        double *eigval = eom_driver_->eigvals_->pointer();
        forall(properties_["OSCILLATOR STRENGTHS"], [eigval](auto &tile, auto &x) {
            double w = eigval[x[1]] - eigval[x[0]];
            tile[x] *= 2.0 / 3.0 * w;
        });
        world_.gop.fence();
    }

    void EOM_RDM::print_oscillators() {
        TArrayD &oscillators = properties_["OSCILLATOR STRENGTHS"];
        TArrayD &x_transitions = properties_["X TRANSITION DIPOLES"];
        TArrayD &y_transitions = properties_["Y TRANSITION DIPOLES"];
        TArrayD &z_transitions = properties_["Z TRANSITION DIPOLES"];


        //TODO: this is a hack to get the data out of the TArrayD and needs to be broadcast to all nodes.
        // However, these arrays are small enough that it is not a big deal.

        // convert to SharedMatrix
        SharedMatrix osc(new Matrix(M_, M_)), x(new Matrix(M_, M_)), y(new Matrix(M_, M_)), z(new Matrix(M_, M_));
        double **oscp = osc->pointer(), **xp, **yp, **zp;

        if (options_.get_str("ROTATE_POLARIZATION_AXIS") == "XYZ"){
            xp = x->pointer();
            yp = y->pointer();
            zp = z->pointer();
        } else if (options_.get_str("ROTATE_POLARIZATION_AXIS") == "ZXY"){
            xp = z->pointer();
            yp = x->pointer();
            zp = y->pointer();
        } else if (options_.get_str("ROTATE_POLARIZATION_AXIS") == "YZX"){
            xp = y->pointer();
            yp = z->pointer();
            zp = x->pointer();
        } else {
            throw PsiException("ROTATE_POLARIZATION_AXIS must be XYZ, ZXY, or YZX", __FILE__, __LINE__);
        }

        forall(properties_["OSCILLATOR STRENGTHS"],
                        [oscp](auto &tile, auto &x) { oscp[x[0]][x[1]] = tile[x]; });
        forall(properties_["X TRANSITION DIPOLES"], [xp](auto &tile, auto &x) { xp[x[0]][x[1]] = tile[x]; });
        forall(properties_["Y TRANSITION DIPOLES"], [yp](auto &tile, auto &x) { yp[x[0]][x[1]] = tile[x]; });
        forall(properties_["Z TRANSITION DIPOLES"], [zp](auto &tile, auto &x) { zp[x[0]][x[1]] = tile[x]; });
        world_.gop.fence();

        double cc_energy = eom_driver_->cc_wfn_->cc_energy_;
        double *eigval = eom_driver_->eigvals_->pointer();

        // Print oscillator strengths
        world_.gop.serial_invoke([=] () {
            const char *eom_type = eom_driver_->eom_type_.c_str();
            for (int ref = 0; ref < M_; ++ref) {
                if (ref == 0) outfile->Printf("\n  ==> %s Initial State Transition Dipoles <==\n\n", eom_type);
                else if (ref == 1) outfile->Printf("\n\n  ==> %s Excited State Transition Dipoles <==\n\n", eom_type);
                else outfile->Printf("  -- Excited State %d --\n", ref);

                outfile->Printf("%5s %20s", "state", "omega (Eh)");
                outfile->Printf(" %20s | %20s | %12s  %12s | %12s  %12s | %12s  %12s\n", "Oscillator Strength", "<0|mu|n><n|mu|0>",
                                "<0|mu_x|n>","<n|mu_x|0>", "<0|mu_y|n>","<n|mu_y|0>", "<0|mu_z|n>","<n|mu_z|0>");
                for (int state = ref; state < M_; ++state) {

                    double oscilator_strength = osc->get(ref, state);
                    double xket = x->get(state, ref),
                           yket = y->get(state, ref),
                           zket = z->get(state, ref);
                    double xbra = x->get(ref, state),
                           ybra = y->get(ref, state),
                           zbra = z->get(ref, state);

                    if (state == ref) oscilator_strength = 0.0; // same state oscillator strength must be zero

                    // get excitation energy
                    double w = eigval[state] - eigval[ref];

                    // print oscillator strength
                    outfile->Printf("%5d %20.12lf ", state, w);
                    if (fabs(oscilator_strength) > 1e-12) outfile->Printf("%20.12lf |", oscilator_strength);
                    else outfile->Printf("-------------------- |");

                    // print dipole strengths
                    double dip_strength = xbra * xket + ybra * yket + zbra * zket;
                    if (fabs(dip_strength) > 1e-12) outfile->Printf(" %20.12lf", dip_strength);
                    else outfile->Printf(" --------------------");

                    // print transition dipole moments
                    if (fabs(xbra) > 1e-8) outfile->Printf(" | %12.8lf", xbra);
                    else outfile->Printf(" | ------------");
                    if (fabs(xket) > 1e-8) outfile->Printf("  %12.8lf", xket);
                    else outfile->Printf("  ------------");

                    if (fabs(ybra) > 1e-8) outfile->Printf(" | %12.8lf", ybra);
                    else outfile->Printf(" | ------------");
                    if (fabs(yket) > 1e-8) outfile->Printf("  %12.8lf", yket);
                    else outfile->Printf("  ------------");

                    if (fabs(zbra) > 1e-8) outfile->Printf(" | %12.8lf", zbra);
                    else outfile->Printf(" | ------------");
                    if (fabs(zket) > 1e-8) outfile->Printf("  %12.8lf", zket);
                    else outfile->Printf("  ------------");
                    outfile->Printf("\n");
                }
                outfile->Printf("\n");
            }
        });
    }
}