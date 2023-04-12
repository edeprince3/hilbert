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

#include "../../include/derived/eom_ea_driver.h"


namespace hilbert {
    hilbert::EOM_EA_Driver::EOM_EA_Driver(const shared_ptr <CC_Cavity> &cc_wfn, Options &options) : EOM_Driver(cc_wfn,
                                                                                                                 options) {

        //TODO: replace this will SCF call through psi4's python interface

        old_e_dip_z_ = cc_wfn_->e_dip_z_;
        old_dipole_self_energy_ = cc_wfn_->average_electric_dipole_self_energy_;

        // reset cavity parameters to get dipole of the anion
        if ( options_["ANION_EDIPZ"].has_changed() && options_["ANION_DSE"].has_changed()) {
            cc_wfn_->e_dip_z_ = options_.get_double("ANION_EDIPZ");
            cc_wfn_->average_electric_dipole_self_energy_ = options_.get_double("ANION_DSE");
        } else {
            outfile->Printf("Warning, settings for dipole moment and dipole self energy not set for QED-EOM-EA-CC");
        }

        // remove T1 transformation
        if ( cc_wfn_->has_t1_integrals_ ) cc_wfn_->transform_integrals(false);

        // some one-electron quantities must be rebuilt
        // to account for updated electric dipole
        //
        // normally, the following happens in update_cavity_terms().
        //

        cc_wfn_->tot_dip_x_ = cc_wfn_->e_dip_x_ + cc_wfn_->nuc_dip_x_;
        cc_wfn_->tot_dip_y_ = cc_wfn_->e_dip_y_ + cc_wfn_->nuc_dip_y_;
        cc_wfn_->tot_dip_z_ = cc_wfn_->e_dip_z_ + cc_wfn_->nuc_dip_z_;

        // e-(n-<d>) contribution 0.5 * 2 (lambda . de) ( lambda . (dn - <d>) )

        double lambda_x = cc_wfn_->cavity_coupling_strength_[0] * sqrt(2.0 * cc_wfn_->cavity_frequency_[0]);
        double lambda_y = cc_wfn_->cavity_coupling_strength_[1] * sqrt(2.0 * cc_wfn_->cavity_frequency_[1]);
        double lambda_z = cc_wfn_->cavity_coupling_strength_[2] * sqrt(2.0 * cc_wfn_->cavity_frequency_[2]);

        size_t nso_ = cc_wfn_->ns_/2;

        // e contribution: lambda . de
        std::shared_ptr<Matrix> el_dipdot (new Matrix(nso_,nso_));
        el_dipdot->zero();
        el_dipdot->axpy(lambda_x,cc_wfn_->dipole_[0]);
        el_dipdot->axpy(lambda_y,cc_wfn_->dipole_[1]);
        el_dipdot->axpy(lambda_z,cc_wfn_->dipole_[2]);

        // n-<d> contribution: lambda . (dn - <d>)
        double nuc_dipdot = 0.0;
        nuc_dipdot += lambda_x * ( cc_wfn_->nuc_dip_x_ - cc_wfn_->tot_dip_x_ );
        nuc_dipdot += lambda_y * ( cc_wfn_->nuc_dip_y_ - cc_wfn_->tot_dip_y_ );
        nuc_dipdot += lambda_z * ( cc_wfn_->nuc_dip_z_ - cc_wfn_->tot_dip_z_ );

        // e-(n-<d>) contribution 0.5 * 2 (lambda . de) ( lambda . (dn - <d>) )
        cc_wfn_->scaled_e_n_dipole_squared_ = std::make_shared<Matrix>(new Matrix(nso_,nso_));
        cc_wfn_->scaled_e_n_dipole_squared_->copy(el_dipdot);
        cc_wfn_->scaled_e_n_dipole_squared_->scale(nuc_dipdot);

        // now, rebuild oe_cavity_terms_

        // e-(n-<d>) contribution to dipole self energy
        double ** scaled_e_n_dipole_squared_p = cc_wfn_->scaled_e_n_dipole_squared_->pointer();

        // one-electron part of electron-electron contribution to dipole self energy
        double ** quadrupole_scaled_sum_p = cc_wfn_->quadrupole_scaled_sum_->pointer();

        cc_wfn_->oe_cavity_terms_ = HelperD::makeTensor(world_, {2L*nso_, 2L*nso_}, false);
        cc_wfn_->oe_cavity_terms_.init_elements([scaled_e_n_dipole_squared_p, quadrupole_scaled_sum_p, nso_](auto &I) {
            if (I[0] < nso_ && I[1] < nso_)
                return scaled_e_n_dipole_squared_p[I[0]][I[1]] - quadrupole_scaled_sum_p[I[0]][I[1]];
            else if (I[0] >= nso_ && I[1] >= nso_)
                return scaled_e_n_dipole_squared_p[I[0]-nso_][I[1]-nso_] - quadrupole_scaled_sum_p[I[0]-nso_][I[1]-nso_];
            else return 0.0;
        });
        world_.gop.fence();

        // retransform T1 integrals
        cc_wfn_->transform_integrals(true);

    }

    void hilbert::EOM_EA_Driver::set_problem_size() {
        singleDim_ = va_ + vb_;
        doubleDim_ = va_*(va_-1)/2*oa_ + vb_*(vb_-1)/2*ob_ + ob_*vb_*va_ + oa_*va_*vb_;
        dim_e_ = singleDim_ + doubleDim_;
        dim_p_ = 0;
        if (include_u1_) dim_p_ += singleDim_;
        if (include_u2_) dim_p_ += doubleDim_;

        N_ = dim_e_ + dim_p_;
    }

    void hilbert::EOM_EA_Driver::print_eom_summary() const {
        // print out energies and amplitude norms
        Printf("\n%5s", "state");
        Printf(" %20s", "total energy (Eh)");
        Printf(" %17s", "ex. energy (Eh)");
        Printf(" %13s", "|l1*r1|");
        Printf(" %13s", "|l2*r2|");
        if (include_u1_) Printf(" %13s", "|m1*s1|");
        if (include_u2_) Printf(" %13s", "|m2*s2|");
        Printf("\n");

        // count number of operators
        int nid = 2; // r1, r2
        if (include_u1_) nid++; // s1
        if (include_u2_) nid++; // s2

        // loop over roots
        for (int i = 0; i < M_; i++) {
            // calculate norms of amplitudes
            double *norms = get_state_norms(i);

            // print energies and excitation energies
            double ee_energy = eigvals_->get(i);
            Printf("%5d %20.12lf %17.12lf ", i, ee_energy,
                   ee_energy - cc_wfn_->cc_energy_);

            // print out the norm of the amplitudes. Ignore if less than 1e-10
            for (int j = 0; j < nid; j++) {
                if (fabs(norms[j]) > 1e-10) Printf("%13.10lf ", norms[j]);
                else Printf("------------- ");
            }
            Printf("\n");
            free(norms);
        }
    }

    void hilbert::EOM_EA_Driver::build_hamiltonian() {

    }

    double *hilbert::EOM_EA_Driver::build_guess() {

        // get properties from the CC_Cavity object
        double const *cavity_frequency_ = cc_wfn_->cavity_frequency_;
        double &enuc_ = cc_wfn_->enuc_;
        double &average_electric_dipole_self_energy_ = cc_wfn_->average_electric_dipole_self_energy_;
        double const *epsilon_ = cc_wfn_->epsilon_;

        // adjust the energy of the CCSD wavefunction to include the anion dipole self energy
        double cc_energy_ = cc_wfn_->cc_energy_;
        cc_energy_ += average_electric_dipole_self_energy_ - old_dipole_self_energy_;

        int oa = (int)oa_, ob = (int)ob_, // occupied alpha/beta
        va = (int)va_, vb = (int)vb_; // virtual alpha/beta

        auto * Hdiag = (double*) calloc(N_, sizeof(double));

        bool singleGuess = false;
        double* ss_diag = nullptr;

        // singles
        //alpha
        int id = 0;
        int id_s = 0;
        for (int a = 0; a < va; a++) {
            if (singleGuess) {
                Hdiag[id] = ss_diag[id_s];
                Hdiag[id] += average_electric_dipole_self_energy_ + enuc_;

                if (include_u1_) {
                    Hdiag[id + dim_e_] = ss_diag[id_s + oa * va + ob * vb];
                    Hdiag[id + dim_e_] += average_electric_dipole_self_energy_ + enuc_;
                }
            } else {
                Hdiag[id] = cc_energy_ + epsilon_[a + o_];
                if (include_u1_) {
                    Hdiag[id + dim_e_] = cc_energy_ + epsilon_[a + o_] + cavity_frequency_[2];
                }

            }
            id++; id_s++;
        }

        //beta
        for (int a = va; a < va + vb; a++) {
            if (singleGuess) {
                Hdiag[id] = ss_diag[id_s];
                Hdiag[id] += average_electric_dipole_self_energy_ + enuc_;
                if (include_u1_) {
                    Hdiag[id + dim_e_] = ss_diag[id_s + oa * va + ob * vb];
                    Hdiag[id + dim_e_] += average_electric_dipole_self_energy_ + enuc_;
                }
            } else {
                Hdiag[id] = cc_energy_ + epsilon_[a + o_];
                if (include_u1_) {
                    Hdiag[id + dim_e_] = cc_energy_ + epsilon_[a + o_] + cavity_frequency_[2];
                }
            }
            id++; id_s++;
        }

        // doubles
        //aaa
        for (int a = 0; a < va; a++) {
            for (int b = a + 1; b < va; b++) {
                for (int j = 0; j < oa; j++) {
                    Hdiag[id] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[j];
                    if ( include_u2_ ){
                        Hdiag[id + dim_e_] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[j] + cavity_frequency_[2];
                    }
                    id++;
                }
            }
        }

        // abb
        for (int a = 0; a < va; a++) {
            for (int b = vb; b < vb + vb; b++) {
                for (int j = oa; j < oa+ob; j++) {
                    Hdiag[id] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[j];
                    if ( include_u2_ ){
                        Hdiag[id + dim_e_] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[j] + cavity_frequency_[2];
                    }
                    id++;
                }
            }
        }
        // baa
        for (int a = vb; a < vb + vb; a++) {
            for (int b = 0; b < va; b++) {
                for (int j = 0; j < oa; j++) {
                    Hdiag[id] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[j];
                    if ( include_u2_ ){
                        Hdiag[id + dim_e_] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[j] + cavity_frequency_[2];
                    }
                    id++;
                }
            }
        }

        //bbb
        for (int a = vb; a < va+vb; a++) {
            for (int b = a + 1; b < va+vb; b++) {
                for (int j = ob; j < oa+ob; j++) {
                    Hdiag[id] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[j];
                    if ( include_u2_ ){
                        Hdiag[id + dim_e_] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[j] + cavity_frequency_[2];
                    }
                    id++;
                }
            }
        }

        if (singleGuess) { free(ss_diag); }
        return Hdiag;
    }

    void hilbert::EOM_EA_Driver::unpack_trial_vectors(size_t L, double **Q) {
        /// unpack the trial vectors from Q into the appropriate blocks of the left/right eigenvectors

        evec_blks_.clear(); // clear the map of trial vectors

        // initialize operators
        evec_blks_["r1_a"] = HelperD::makeTensor(world_, {L, va_}, false); // alpha singles excitations
        evec_blks_["r1_b"] = HelperD::makeTensor(world_, {L, vb_}, false); // beta singles excitations
        evec_blks_["r2_aaa"] = HelperD::makeTensor(world_, {L, va_, va_, oa_}, false); // alpha doubles excitations
        evec_blks_["r2_abb"] = HelperD::makeTensor(world_, {L, va_, vb_, ob_}, false); // alpha-beta doubles excitations
        evec_blks_["r2_baa"] = HelperD::makeTensor(world_, {L, vb_, va_, oa_}, false); // beta-alpha doubles excitations
        evec_blks_["r2_bbb"] = HelperD::makeTensor(world_, {L, vb_, vb_, ob_}, false); // beta doubles excitations

        if (include_u1_){
            evec_blks_["s1_a"] = HelperD::makeTensor(world_, {L, va_}, false); // alpha singles excitations + hw
            evec_blks_["s1_b"] = HelperD::makeTensor(world_, {L, vb_}, false); // beta singles excitations + hw
        }
        if (include_u2_){
            evec_blks_["s2_aaa"] = HelperD::makeTensor(world_, {L, va_, va_, oa_}, false); // alpha doubles excitations + hw
            evec_blks_["s2_abb"] = HelperD::makeTensor(world_, {L, va_, vb_, ob_}, false); // alpha-beta doubles excitations + hw
            evec_blks_["s2_baa"] = HelperD::makeTensor(world_, {L, vb_, va_, oa_}, false); // beta-alpha doubles excitations + hw
            evec_blks_["s2_bbb"] = HelperD::makeTensor(world_, {L, vb_, vb_, ob_}, false); // beta doubles excitations + hw
        }

        size_t oa = oa_, ob = ob_, // number of occupied orbitals
        va = va_, vb = vb_; // number of virtual beta orbitals

        size_t offset = 0;
        /// pack electronic trial vectors into TA::TArrayD

        // pack singles into tensors
        {
            evec_blks_["r1_a"].init_elements([Q, oa, offset](auto &I) {
                return Q[I[0]][I[1] + offset];
            }); offset += va;

            evec_blks_["r1_b"].init_elements([Q, ob, offset](auto &I) {
                return Q[I[0]][I[1] + offset];
            }); offset += vb;
        }

        // pack redundant doubles into tensors
        {
            evec_blks_["r2_aaa"].init_elements([Q, oa, ob, vb, va, offset, this](auto &I) {
                size_t a = I[1], b = I[2], j = I[3];
                if (a == b) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b);
                if (a > b) std::swap(a, b);

                size_t ab_off = this->sqr_2_tri_idx(a, b, va);
                size_t upper_tri_index = ab_off*oa + j;

                double value = Q[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;
            }); offset += va*(va-1)/2*oa;

            evec_blks_["r2_abb"].init_elements([Q, oa, ob, vb, va, offset](auto &I) {
                return Q[I[0]][I[1]*vb*ob + I[2]*ob + I[3] + offset];
            }); offset += va*vb*ob;

            evec_blks_["r2_baa"].init_elements([Q, oa, ob, vb, va, offset](auto &I) {
                return Q[I[0]][I[1]*va*oa + I[2]*oa + I[3] + offset];
            }); offset += vb*va*oa;

            evec_blks_["r2_bbb"].init_elements([Q, oa, ob, vb, va, offset, this](auto &I) {
                size_t a = I[1], b = I[2], j = I[3];
                if (a == b) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b);
                if (a > b) std::swap(a, b);

                size_t ab_off = this->sqr_2_tri_idx(a, b, vb);
                size_t upper_tri_index = ab_off*ob + j;

                double value = Q[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;
            }); offset += vb*(vb-1)/2*ob;
        }

        /// pack photonic trial vectors into TA::TArrayD
        if (include_u1_){
            // pack singles + hw operator
            evec_blks_["s1_a"].init_elements([Q, oa, ob, vb, va, offset](auto &I) {
                return Q[I[0]][I[1] + offset];
            }); offset += va;

            evec_blks_["s1_b"].init_elements([Q, oa, ob, vb, va, offset](auto &I) {
                return Q[I[0]][I[1] + offset];
            }); offset += vb;
        }
        if (include_u2_){
            // pack redundant doubles + hw operator
            evec_blks_["s2_aaa"].init_elements([Q, oa, ob, vb, va, offset, this](auto &I) {
                size_t a = I[1], b = I[2], j = I[3];
                if (a == b) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b);
                if (a > b) std::swap(a, b);

                size_t ab_off = this->sqr_2_tri_idx(a, b, va);
                size_t upper_tri_index = ab_off*oa + j;

                double value = Q[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;
            }); offset += va*(va-1)/2*oa;

            evec_blks_["s2_abb"].init_elements([Q, oa, ob, vb, va, offset](auto &I) {
                return Q[I[0]][I[1]*vb*ob + I[2]*ob + I[3] + offset];
            }); offset += va*vb*ob;

            evec_blks_["s2_baa"].init_elements([Q, oa, ob, vb, va, offset](auto &I) {
                return Q[I[0]][I[1]*va*oa + I[2]*oa + I[3] + offset];
            }); offset += vb*va*oa;

            evec_blks_["s2_bbb"].init_elements([Q, oa, ob, vb, va, offset, this](auto &I) {
                size_t a = I[1], b = I[2], j = I[3];
                if (a == b) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b);
                if (a > b) std::swap(a, b);

                size_t ab_off = this->sqr_2_tri_idx(a, b, vb);
                size_t upper_tri_index = ab_off*ob + j;

                double value = Q[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;
            }); offset += vb*(vb-1)/2*ob;
        }
        world_.gop.fence();

        /// initialize left operators from right
        evec_blks_["l1_a"]("I, e") = evec_blks_["r1_a"]("I, e");
        evec_blks_["l1_b"]("I, e") = evec_blks_["r1_b"]("I, e");
        evec_blks_["l2_aaa"]("I, n, f, e") = evec_blks_["r2_aaa"]("I, e, f, n");
        evec_blks_["l2_bba"]("I, n, f, e") = evec_blks_["r2_abb"]("I, e, f, n");
        evec_blks_["l2_aab"]("I, n, f, e") = evec_blks_["r2_baa"]("I, e, f, n");
        evec_blks_["l2_bbb"]("I, n, f, e") = evec_blks_["r2_bbb"]("I, e, f, n");

        if ( include_u1_ ){
            evec_blks_["m1_a"]("I, e") = evec_blks_["s1_a"]("I, e");
            evec_blks_["m1_b"]("I, e") = evec_blks_["s1_b"]("I, e");
        }
        if ( include_u2_ ){
            evec_blks_["m2_aaa"]("I, n, f, e") = evec_blks_["s2_aaa"]("I, e, f, n");
            evec_blks_["m2_bba"]("I, n, f, e") = evec_blks_["s2_abb"]("I, e, f, n");
            evec_blks_["m2_aab"]("I, n, f, e") = evec_blks_["s2_baa"]("I, e, f, n");
            evec_blks_["m2_bbb"]("I, n, f, e") = evec_blks_["s2_bbb"]("I, e, f, n");
        }
        world_.gop.fence();
    }

    void hilbert::EOM_EA_Driver::pack_sigma_vectors(size_t L, double **sigmar, double **sigmal) {
        size_t oa = oa_, ob = ob_, // number of occupied orbitals
        va = va_, vb = vb_; // number of virtual beta orbitals

        /// unpack electronic part of right sigma vectors
        size_t offset = 0;
        // singles
        {
            HelperD::forall(sigvec_blks_["sigmar1_a"], [sigmar, oa, offset](auto &tile, auto &x) {
                sigmar[x[1] + offset][x[0]] = tile[x];
            }); offset += va;

            HelperD::forall(sigvec_blks_["sigmar1_b"], [sigmar, ob, offset](auto &tile, auto &x) {
                sigmar[x[1] + offset][x[0]] = tile[x];
            }); offset += vb;
        }

        // doubles
        {
            HelperD::forall(sigvec_blks_["sigmar2_aaa"], [sigmar, oa, va, ob, vb, offset, this](auto &tile, auto &x) {
                if (x[1] < x[2]) {
                    size_t ab_off = this->sqr_2_tri_idx(x[1], x[2], va);
                    size_t sigmar_idx = ab_off*oa + x[3] + offset;

                    sigmar[sigmar_idx][x[0]] = tile[x];
                }
            }); offset += va*(va-1)/2*oa;

            HelperD::forall(sigvec_blks_["sigmar2_abb"], [sigmar, oa, ob, vb, offset](auto &tile, auto &x) {
                sigmar[x[1] * vb * ob + x[2] * ob + x[3] + offset][x[0]] = tile[x];
            }); offset += va*vb*ob;

            HelperD::forall(sigvec_blks_["sigmar2_baa"], [sigmar, oa, ob, va, vb, offset](auto &tile, auto &x) {
                sigmar[x[1] * va * oa + x[2] * oa + x[3] + offset][x[0]] = tile[x];
            }); offset += vb*va*oa;

            HelperD::forall(sigvec_blks_["sigmar2_bbb"], [sigmar, oa, ob, va, vb, offset, this](auto &tile, auto &x) {
                if (x[1] < x[2]) {
                    size_t ab_off = this->sqr_2_tri_idx(x[1], x[2], vb);
                    size_t sigmar_idx = ab_off*ob + x[3] + offset;

                    sigmar[sigmar_idx][x[0]] = tile[x];
                }
            }); offset += vb*(vb-1)/2*ob;
        }

        /// unpack photonic part of right sigma vectors
        // singles + hw
        if (include_u1_) {
            HelperD::forall(sigvec_blks_["sigmas1_a"], [sigmar, oa, offset](auto &tile, auto &x) {
                sigmar[x[1] + offset][x[0]] = tile[x];
            }); offset += va;

            HelperD::forall(sigvec_blks_["sigmas1_b"], [sigmar, ob, offset](auto &tile, auto &x) {
                sigmar[x[1] + offset][x[0]] = tile[x];
            }); offset += vb;
        }

        // doubles + hw
        if (include_u2_){
            HelperD::forall(sigvec_blks_["sigmas2_aaa"], [sigmar, oa, ob, va, vb, offset, this](auto &tile, auto &x) {
                if (x[1] < x[2]) {
                    size_t ab_off = this->sqr_2_tri_idx(x[1], x[2], va);
                    size_t sigmar_idx = ab_off*oa + x[3] + offset;

                    sigmar[sigmar_idx][x[0]] = tile[x];
                }
            }); offset += va*(va-1)/2*oa;

            HelperD::forall(sigvec_blks_["sigmas2_abb"], [sigmar, oa, ob, vb, offset](auto &tile, auto &x) {
                sigmar[x[1] * vb * ob + x[2] * ob + x[3] + offset][x[0]] = tile[x];
            }); offset += va*vb*ob;

            HelperD::forall(sigvec_blks_["sigmas2_baa"], [sigmar, oa, ob, va, vb, offset](auto &tile, auto &x) {
                sigmar[x[1] * va * oa + x[2] * oa + x[3] + offset][x[0]] = tile[x];
            }); offset += vb*va*oa;

            HelperD::forall(sigvec_blks_["sigmas2_bbb"], [sigmar, oa, ob, va, vb, offset, this](auto &tile, auto &x) {
                if (x[1] < x[2]) {
                    size_t ab_off = this->sqr_2_tri_idx(x[1], x[2], vb);
                    size_t sigmar_idx = ab_off*ob + x[3] + offset;

                    sigmar[sigmar_idx][x[0]] = tile[x];
                }
            }); offset += vb*(vb-1)/2*ob;
        }

        /// unpack electronic part of left sigma vectors
        offset = 0;
        // singles
        {
            HelperD::forall(sigvec_blks_["sigmal1_a"], [sigmal, oa, offset](auto &tile, auto &x) {
                sigmal[x[1] + offset][x[0]] = tile[x];
            }); offset += va;

            HelperD::forall(sigvec_blks_["sigmal1_b"], [sigmal, ob, offset](auto &tile, auto &x) {
                sigmal[x[1] + offset][x[0]] = tile[x];
            }); offset += vb;
        }

        // doubles
        {
            HelperD::forall(sigvec_blks_["sigmal2_aaa"], [sigmal, oa, ob, va, vb, offset, this](auto &tile, auto &x) {
                if (x[1] < x[2]) {
                    size_t ab_off = this->sqr_2_tri_idx(x[1], x[2], va);
                    size_t sigmal_idx = ab_off*oa + x[3] + offset;

                    sigmal[sigmal_idx][x[0]] = tile[x];
                }
            }); offset += va*(va-1)/2*oa;

            HelperD::forall(sigvec_blks_["sigmal2_abb"], [sigmal, oa, ob, vb, offset](auto &tile, auto &x) {
                sigmal[x[1] * vb * ob + x[2] * ob + x[3] + offset][x[0]] = tile[x];
            }); offset += va*vb*ob;

            HelperD::forall(sigvec_blks_["sigmal2_baa"], [sigmal, oa, ob, va, vb, offset](auto &tile, auto &x) {
                sigmal[x[1] * va * oa + x[2] * oa + x[3] + offset][x[0]] = tile[x];
            }); offset += vb*va*oa;

            HelperD::forall(sigvec_blks_["sigmal2_bbb"], [sigmal, oa, ob, va, vb, offset, this](auto &tile, auto &x) {
                if (x[1] < x[2]) {
                    size_t ab_off = this->sqr_2_tri_idx(x[1], x[2], vb);
                    size_t sigmal_idx = ab_off*ob + x[3] + offset;

                    sigmal[sigmal_idx][x[0]] = tile[x];
                }
            }); offset += vb*(vb-1)/2*ob;
        }

        /// unpack photonic part of left sigma vectors
        // singles + hw
        if (include_u1_) {
            HelperD::forall(sigvec_blks_["sigmam1_a"], [sigmal, oa, offset](auto &tile, auto &x) {
                sigmal[x[1] + offset][x[0]] = tile[x];
            }); offset += va;

            HelperD::forall(sigvec_blks_["sigmam1_b"], [sigmal, ob, offset](auto &tile, auto &x) {
                sigmal[x[1] + offset][x[0]] = tile[x];
            }); offset += vb;
        }

        // doubles + hw
        if (include_u2_) {
            HelperD::forall(sigvec_blks_["sigmam2_aaa"], [sigmal, oa, ob, va, vb, offset, this](auto &tile, auto &x) {
                if (x[1] < x[2]) {
                    size_t ab_off = this->sqr_2_tri_idx(x[1], x[2], va);
                    size_t sigmal_idx = ab_off*oa + x[3] + offset;

                    sigmal[sigmal_idx][x[0]] = tile[x];
                }
            }); offset += va*(va-1)/2*oa;

            HelperD::forall(sigvec_blks_["sigmam2_abb"], [sigmal, oa, ob, vb, offset](auto &tile, auto &x) {
                sigmal[x[1] * vb * ob + x[2] * ob + x[3] + offset][x[0]] = tile[x];
            }); offset += va*vb*ob;

            HelperD::forall(sigvec_blks_["sigmam2_baa"], [sigmal, oa, ob, va, vb, offset](auto &tile, auto &x) {
                sigmal[x[1] * va * oa + x[2] * oa + x[3] + offset][x[0]] = tile[x];
            }); offset += vb*va*oa;

            HelperD::forall(sigvec_blks_["sigmam2_bbb"], [sigmal, oa, ob, va, vb, offset, this](auto &tile, auto &x) {
                if (x[1] < x[2]) {
                    size_t ab_off = this->sqr_2_tri_idx(x[1], x[2], vb);
                    size_t sigmal_idx = ab_off*ob + x[3] + offset;

                    sigmal[sigmal_idx][x[0]] = tile[x];
                }
            }); offset += vb*(vb-1)/2*ob;
        }
    }

    double *hilbert::EOM_EA_Driver::get_state_norms(size_t i) const {
        // get number of amplitudes
        int numNorms = 2;
        if (include_u1_) numNorms++;
        if (include_u2_) numNorms++;

        double** rerp_ = revec_->pointer();
        double** relp_ = levec_->pointer();

        auto *norms = (double *) calloc(numNorms, sizeof(double));
        int nid = 0;

        double totalNorm = fabs(C_DDOT(N_, rerp_[i], 1, relp_[i], 1));
        norms[nid] = C_DDOT(singleDim_, rerp_[i] + 1, 1, relp_[i] + 1, 1) / totalNorm;
        norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0; nid++;

        norms[nid] = C_DDOT(doubleDim_, rerp_[i] + singleDim_ + 1, 1, relp_[i] + singleDim_ + 1, 1) / totalNorm;
        norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0; nid++;

        if (include_u1_) {
            norms[nid] = C_DDOT(singleDim_, rerp_[i] + dim_e_ + 1, 1, relp_[i] + dim_e_ + 1, 1) / totalNorm;
            norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0; nid++;
        }

        if (include_u2_) {
            norms[nid] = C_DDOT(doubleDim_, rerp_[i] + dim_e_ + singleDim_ + 1, 1, relp_[i] + dim_e_ + singleDim_ + 1, 1) / totalNorm;
            norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0;
        }

        return norms;
    }

    void hilbert::EOM_EA_Driver::unpack_eigenvectors() {
        exit(0);
    }

    void EOM_EA_Driver::print_dominant_transitions() {

    }
} // cc_cavity