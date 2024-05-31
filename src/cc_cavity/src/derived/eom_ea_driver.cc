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

#include "cc_cavity/include/derived/eom_ea_driver.h"


namespace hilbert {
    hilbert::EOM_EA_CCSD::EOM_EA_CCSD(const shared_ptr <CC_Cavity> &cc_wfn, Options &options) : EOM_Driver(cc_wfn,
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

    void hilbert::EOM_EA_CCSD::set_problem_size() {
        singleDim_ = va_ + vb_;
        doubleDim_ = (va_*va_-va_)/2*oa_  // aaa
                   +  va_*vb_*ob_         // abb
                   +  va_*vb_*oa_         // aba
                   + (vb_*vb_-vb_)/2*ob_; // bbb
        dim_e_ = singleDim_ + doubleDim_;
        nops_ = 2;

        dim_p_ = 0;
//        if (include_u1_) {
//            dim_p_ += singleDim_;
//            nops_++;
//        }
//        if (include_u2_) {
//            dim_p_ += doubleDim_;
//            nops_++;
//        }

        N_ = dim_e_ + dim_p_;
    }

    void hilbert::EOM_EA_CCSD::print_eom_header() {
        // print out energies and amplitude norms
        Printf("\n%5s", "state");
        Printf(" %20s", "total energy (Eh)");
        Printf(" %17s", "ex. energy (Eh)");
        Printf(" %13s", "|l1*r1|");
        Printf(" %13s", "|l2*r2|");
//        if (include_u1_) Printf(" %13s", "|m1*s1|");
//        if (include_u2_) Printf(" %13s", "|m2*s2|");
        Printf("\n");

        // set reference energy
        ground_energy_ref_ = cc_wfn_->cc_energy_; //eigvals_->get(0);
    }

    double *hilbert::EOM_EA_CCSD::build_preconditioner() {

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
        size_t aaa = (va*va-va)/2*oa, abb = va*vb*ob, baa = vb*va*oa, bbb = (vb*vb-vb)/2*ob;

        auto* Hdiag = (double*) calloc(N_, sizeof(double));

        bool singleGuess = false;
        double* ss_diag = nullptr;

        // singles
        //alpha
        int id = 0;
        int id_s = 0;
        for (int a = 0; a < va; a++) {
            Hdiag[id] = cc_energy_ + epsilon_[a + o_];
            id++; id_s++;
        }

        //beta
        for (int a = va; a < va + vb; a++) {
            Hdiag[id] = cc_energy_ + epsilon_[a + o_];
            id++; id_s++;
        }

        // doubles
        //aaa
        for (int a = 0; a < va; a++) {
            for (int b = a+1; b < va; b++) {
                for (int j = 0; j < oa; j++) {
                    if (a == b) Hdiag[id] = 0;
                    else
                        Hdiag[id] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[j];
                    id++;
                }
            }
        }

        // abb
        for (int a = 0; a < va; a++) {
            for (int b = va; b < va + vb; b++) {
                for (int j = oa; j < oa+ob; j++) {
                    Hdiag[id] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[j];
                    id++;
                }
            }
        }

        // aba
        for (int a = 0; a < va; a++) {
            for (int b = va; b < va+vb; b++) {
                for (int j = 0; j < oa; j++) {
                    Hdiag[id] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[j];
                    id++;
                }
            }
        }

        //bbb
        for (int a = va; a < va+vb; a++) {
            for (int b = a+1; b < va+vb; b++) {
                for (int j = oa; j < oa+ob; j++) {
                    if (a == b) Hdiag[id] = 0;
                    else
                        Hdiag[id] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[j];
                    id++;
                }
            }
        }

        if (singleGuess) { free(ss_diag); }
        return Hdiag;
    }

    void hilbert::EOM_EA_CCSD::unpack_trial_vectors(size_t L, double **Q) {
        /// unpack the trial vectors from Q into the appropriate blocks of the left/right eigenvectors

        evec_blks_.clear(); // clear the map of trial vectors

        // initialize operators
        evec_blks_["r1_a"] = HelperD::makeTensor(world_, {L, va_}, false); // alpha singles excitations
        evec_blks_["r1_b"] = HelperD::makeTensor(world_, {L, vb_}, false); // beta singles excitations
        evec_blks_["r2_aaa"] = HelperD::makeTensor(world_, {L, va_, va_, oa_}, false); // alpha doubles excitations
        evec_blks_["r2_abb"] = HelperD::makeTensor(world_, {L, va_, vb_, ob_}, false); // alpha-beta doubles excitations
        evec_blks_["r2_aba"] = HelperD::makeTensor(world_, {L, va_, vb_, oa_}, false); // beta-alpha doubles excitations
        evec_blks_["r2_bbb"] = HelperD::makeTensor(world_, {L, vb_, vb_, ob_}, false); // beta doubles excitations

        size_t oa = oa_, ob = ob_, // number of occupied orbitals
        va = va_, vb = vb_; // number of virtual beta orbitals
        size_t aaa = (va*va-va)/2*oa, abb = va*vb*ob, aba = va*vb*oa, bbb = (vb*vb-vb)/2*ob; // offsets for indexing into Q

        size_t offset = 0;
        /// pack electronic trial vectors into TA::TArrayD

        // pack singles into tensors
        {
            evec_blks_["r1_a"].init_elements([Q, offset](auto &I) {
                size_t trial = I[0], e = I[1];
                return Q[trial][e + offset];
            }); offset += va;

            evec_blks_["r1_b"].init_elements([Q, offset](auto &I) {
                size_t trial = I[0], e = I[1];
                return Q[trial][e + offset];
            }); offset += vb;
        }

        // pack redundant doubles into tensors
        {
            evec_blks_["r2_aaa"].init_elements([Q, oa, ob, vb, va, offset](auto &I) {
                size_t trial = I[0], e = I[1], f = I[2], m = I[3];

                if (e == f) return 0.0; // return 0 for diagonal elements

                bool change_sign = (e > f);
                if (change_sign) std::swap(e, f);

                size_t ef_off = EOM_Driver::sqr_2_tri_idx(e, f, va);
                size_t upper_tri_index = ef_off*oa + m + offset;

                double value = Q[trial][upper_tri_index];
                return change_sign ? -value : value;
            }); offset += aaa;

            evec_blks_["r2_abb"].init_elements([Q, oa, ob, vb, va, offset](auto &I) {
                size_t trial = I[0], e = I[1], f = I[2], m = I[3];
                return Q[trial][e*vb*ob + f*ob + m + offset];
            }); offset += abb;

            evec_blks_["r2_aba"].init_elements([Q, oa, ob, vb, va, offset](auto &I) {
                size_t trial = I[0], e = I[1], f = I[2], m = I[3];
                return Q[trial][e*vb*oa + f*oa + m + offset];
            }); offset += aba;

            evec_blks_["r2_bbb"].init_elements([Q, oa, ob, vb, va, offset](auto &I) {
                size_t trial = I[0], e = I[1], f = I[2], m = I[3];

                if (e == f) return 0.0; // return 0 for diagonal elements

                bool change_sign = (e > f);
                if (change_sign) std::swap(e, f);

                size_t ef_off = EOM_Driver::sqr_2_tri_idx(e, f, vb);
                size_t upper_tri_index = ef_off*ob + m + offset;

                double value = Q[trial][upper_tri_index];
                return change_sign ? -value : value;
            }); offset += bbb;
        }
        world_.gop.fence();

        /// initialize left operators from right
        evec_blks_["l1_a"]("I, e") = evec_blks_["r1_a"]("I, e");
        evec_blks_["l1_b"]("I, e") = evec_blks_["r1_b"]("I, e");
        evec_blks_["l2_aaa"]("I, m, e, f") = evec_blks_["r2_aaa"]("I, e, f, m");
        evec_blks_["l2_bab"]("I, m, e, f") = evec_blks_["r2_abb"]("I, e, f, m");
        evec_blks_["l2_aab"]("I, m, e, f") = evec_blks_["r2_aba"]("I, e, f, m");
        evec_blks_["l2_bbb"]("I, m, e, f") = evec_blks_["r2_bbb"]("I, e, f, m");
        world_.gop.fence();
    }

    void hilbert::EOM_EA_CCSD::pack_sigma_vectors(size_t L, double **sigmar, double **sigmal) {
        size_t oa = oa_, ob = ob_, // number of occupied orbitals
        va = va_, vb = vb_; // number of virtual beta orbitals

        size_t aaa = (va*va-va)/2*oa, abb = va*vb*ob, aba = va*vb*oa, bbb = (vb*vb-vb)/2*ob;

        /// unpack electronic part of right sigma vectors
        size_t offset = 0;

        // singles
        {
            HelperD::forall(sigvec_blks_["sigmar1_a"], [sigmar, oa, offset](auto &tile, auto &x) {
                size_t trial = x[0], e = x[1];
                sigmar[e + offset][trial] = tile[x];
            }); offset += va;

            HelperD::forall(sigvec_blks_["sigmar1_b"], [sigmar, ob, offset](auto &tile, auto &x) {
                size_t trial = x[0], e = x[1];
                sigmar[e + offset][trial] = tile[x];
            }); offset += vb;
        }

        // doubles
        {
            HelperD::forall(sigvec_blks_["sigmar2_aaa"], [sigmar, oa, va, ob, vb, offset, this](auto &tile, auto &x) {
                size_t trial = x[0], e = x[1], f = x[2], m = x[3];
                if (e < f) {
                    size_t ef_off = EOM_Driver::sqr_2_tri_idx(e, f, va);
                    size_t sigmar_idx = ef_off*oa + m + offset;
                    sigmar[sigmar_idx][trial] = tile[x];
                }
            }); offset += aaa;

            HelperD::forall(sigvec_blks_["sigmar2_abb"], [sigmar, oa, ob, vb, offset](auto &tile, auto &x) {
                size_t trial = x[0], e = x[1], f = x[2], m = x[3];
                sigmar[e*vb*ob + f*ob + m + offset][trial] = tile[x];
            }); offset += abb;

            HelperD::forall(sigvec_blks_["sigmar2_aba"], [sigmar, oa, ob, va, vb, offset](auto &tile, auto &x) {
                size_t trial = x[0], e = x[1], f = x[2], m = x[3];
                sigmar[e*vb*oa + f*oa + m + offset][trial] = tile[x];
            }); offset += aba;

            HelperD::forall(sigvec_blks_["sigmar2_bbb"], [sigmar, oa, ob, va, vb, offset, this](auto &tile, auto &x) {
                size_t trial = x[0], e = x[1], f = x[2], m = x[3];
                if (e < f) {
                    size_t ef_off = EOM_Driver::sqr_2_tri_idx(e, f, vb);
                    size_t sigmar_idx = ef_off*ob + m + offset;
                    sigmar[sigmar_idx][trial] = tile[x];
                }
            }); offset += bbb;
        }

        /// unpack electronic part of left sigma vectors
        offset = 0;
        // singles
        {
            HelperD::forall(sigvec_blks_["sigmal1_a"], [sigmal, oa, offset](auto &tile, auto &x) {
                size_t trial = x[0], e = x[1];
                sigmal[e + offset][trial] = tile[x];
            }); offset += va;

            HelperD::forall(sigvec_blks_["sigmal1_b"], [sigmal, ob, offset](auto &tile, auto &x) {
                size_t trial = x[0], e = x[1];
                sigmal[e + offset][trial] = tile[x];
            }); offset += vb;
        }

        // doubles
        {
            HelperD::forall(sigvec_blks_["sigmal2_aaa"], [sigmal, oa, ob, va, vb, offset, this](auto &tile, auto &x) {
                size_t trial = x[0], e = x[1], f = x[2], m = x[3];
                if (e < f) {
                    size_t ef_off = EOM_Driver::sqr_2_tri_idx(e, f, va);
                    size_t sigmal_idx = ef_off*oa + m + offset;
                    sigmal[sigmal_idx][trial] = tile[x];
                }
            }); offset += aaa;

            HelperD::forall(sigvec_blks_["sigmal2_abb"], [sigmal, oa, ob, vb, offset](auto &tile, auto &x) {
                size_t trial = x[0], e = x[1], f = x[2], m = x[3];
                sigmal[e*vb*ob + f*ob + m + offset][trial] = tile[x];
            }); offset += abb;

            HelperD::forall(sigvec_blks_["sigmal2_aba"], [sigmal, oa, ob, va, vb, offset](auto &tile, auto &x) {
                size_t trial = x[0], e = x[1], f = x[2], m = x[3];
                sigmal[e*vb*oa + f*oa + m + offset][trial] = tile[x];
            }); offset += aba;

            HelperD::forall(sigvec_blks_["sigmal2_bbb"], [sigmal, oa, ob, va, vb, offset, this](auto &tile, auto &x) {
                size_t trial = x[0], e = x[1], f = x[2], m = x[3];
                if (e < f) {
                    size_t ef_off = EOM_Driver::sqr_2_tri_idx(e, f, vb);
                    size_t sigmal_idx = ef_off*ob + m + offset;
                    sigmal[sigmal_idx][trial] = tile[x];
                }
            }); offset += bbb;
        }

    }

    double *hilbert::EOM_EA_CCSD::get_state_norms(size_t i) const {
        // get number of amplitudes
        int numNorms = 2;
        if (include_u1_) numNorms++;
        if (include_u2_) numNorms++;

        double** rerp_ = revec_->pointer();
        double** relp_ = levec_->pointer();

        auto *norms = (double *) calloc(numNorms, sizeof(double));
        int nid = 0;

        double totalNorm = fabs(C_DDOT(N_, rerp_[i], 1, relp_[i], 1));
        norms[nid] = C_DDOT(singleDim_, rerp_[i], 1, relp_[i], 1) / totalNorm;
        norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0; nid++;

        norms[nid] = C_DDOT(doubleDim_, rerp_[i] + singleDim_, 1, relp_[i] + singleDim_, 1) / totalNorm;
        norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0; nid++;

//        if (include_u1_) {
//            norms[nid] = C_DDOT(singleDim_, rerp_[i] + dim_e_, 1, relp_[i] + dim_e_, 1) / totalNorm;
//            norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0; nid++;
//        }
//
//        if (include_u2_) {
//            norms[nid] = C_DDOT(doubleDim_, rerp_[i] + dim_e_ + singleDim_, 1, relp_[i] + dim_e_ + singleDim_, 1) / totalNorm;
//            norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0; nid++;
//        }

        return norms;
    }

    void hilbert::EOM_EA_CCSD::unpack_eigenvectors() {
        return;
    }

    EOM_Driver::DominantTransitionsType EOM_EA_CCSD::find_dominant_transitions(size_t I) {
        // this is not implemented for EA-EOM-CC. Return empty vector
        return {};
    }
} // cc_cavity