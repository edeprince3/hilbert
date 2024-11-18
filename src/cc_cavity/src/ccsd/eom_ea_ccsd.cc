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

#include "cc_cavity/include/ccsd/eom_ea_ccsd.h"

namespace hilbert {
    EOM_EA_CCSD::EOM_EA_CCSD(shared_ptr <CC_Cavity> &cc_wfn, Options &options) : EOM_Driver(cc_wfn, options) {
        if (cc_wfn_->nalpha() != cc_wfn_->nbeta()) {
            throw PsiException("EOM-EA-CCSD requires a singlet reference wavefunction",__FILE__,__LINE__);
        }
    }

    void EOM_EA_CCSD::set_problem_size() {
        singleDim_ = va_;
        doubleDim_ = (va_*va_-va_)*oa_ / 2  // aaa
                   +  va_*vb_*ob_;         // abb
        dim_e_ = singleDim_ + doubleDim_;
        nops_ = 2;

        dim_p_ = 0;
        N_ = dim_e_ + dim_p_;
    }

    void EOM_EA_CCSD::print_eom_header() {
        // print out energies and amplitude norms
        Printf("\n%5s", "state");
        Printf(" %20s", "total energy (Eh)");
        Printf(" %17s", "ex. energy (Eh)");
        Printf(" %13s", "|l1*r1|");
        Printf(" %13s", "|l2*r2|");
        Printf("\n");
    }

    double *EOM_EA_CCSD::build_preconditioner() {

        // get properties from the CC_Cavity object
        double const *cavity_frequency_ = cc_wfn_->cavity_frequency_;
        double &enuc_ = cc_wfn_->enuc_;
        double const *epsilon_ = cc_wfn_->epsilon_;

        // adjust the energy of the CCSD wavefunction to include the anion dipole self energy
        double cc_energy_ = cc_wfn_->cc_energy_;

        int oa = (int)oa_, ob = (int)ob_, // occupied alpha/beta
        va = (int)va_, vb = (int)vb_; // virtual alpha/beta
        size_t aaa = (va*va-va)*oa/2, abb = va*vb*ob, baa = vb*va*oa, bbb = (vb*vb-vb)*ob/2;

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

        if (singleGuess) { free(ss_diag); }
        return Hdiag;
    }

    void EOM_EA_CCSD::unpack_trial_vectors(size_t L, double **Q) {
        /// unpack the trial vectors from Q into the appropriate blocks of the left/right eigenvectors

        evec_blks_.clear(); // clear the map of trial vectors

        // initialize operators
        evec_blks_["r1_a"] = makeTensor(world_, {L, va_}, false); // alpha singles excitations
        evec_blks_["r2_aaa"] = makeTensor(world_, {L, va_, va_, oa_}, false); // alpha doubles excitations
        evec_blks_["r2_abb"] = makeTensor(world_, {L, va_, vb_, ob_}, false); // alpha-beta doubles excitations

        size_t oa = oa_, ob = ob_, // number of occupied orbitals
        va = va_, vb = vb_; // number of virtual beta orbitals
        size_t aaa = (va*va-va)*oa/2, abb = va*vb*ob; // offsets for indexing into Q

        size_t offset = 0;
        /// pack electronic trial vectors into TA::TArrayD

        // pack singles into tensors
        {
            evec_blks_["r1_a"].init_elements([Q, offset](auto &I) {
                size_t trial = I[0], e = I[1];
                return Q[trial][e + offset];
            }); offset += va;
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
        }
        world_.gop.fence();

        /// initialize left operators from right
        evec_blks_["l1_a"]("I, e") = evec_blks_["r1_a"]("I, e");
        evec_blks_["l2_aaa"]("I, m, e, f") = evec_blks_["r2_aaa"]("I, e, f, m");
        evec_blks_["l2_bab"]("I, m, e, f") = evec_blks_["r2_abb"]("I, e, f, m");

        world_.gop.fence();
    }

    void EOM_EA_CCSD::build_Hc_cH(size_t L) {
        // transform integrals if needed with t1 amplitudes
        if (!cc_wfn_->has_t1_integrals_)
            cc_wfn_->transform_integrals(true);

        /// initialize sigma vectors for this trial
        sigvec_blks_.clear();

        sigvec_blks_["sigmar1_a"]   = makeTensor(world_, {L, va_}, true);
        sigvec_blks_["sigmal1_a"]   = makeTensor(world_, {L, va_}, true);

        sigvec_blks_["sigmar2_aaa"]   = makeTensor(world_, {L, va_, va_, oa_}, true);
        sigvec_blks_["sigmal2_aaa"]   = makeTensor(world_, {L, va_, va_, oa_}, true);

        sigvec_blks_["sigmar2_abb"]   = makeTensor(world_, {L, va_, vb_, ob_}, true);
        sigvec_blks_["sigmal2_abb"]   = makeTensor(world_, {L, va_, vb_, ob_}, true);

        // run sigma vector builds
        sigma_ea_00_1();
        sigma_ea_00_2();
        sigma_ea_00_3();
        sigma_ea_00_4();
        world_.gop.fence();

        /// add constant terms to sigma vectors
        double nuclear = cc_wfn_->enuc_;

        // set singles
        sigvec_blks_["sigmar1_a"]("I, a") += evec_blks_["r1_a"]("I, a") * nuclear;
        sigvec_blks_["sigmal1_a"]("I, a") += evec_blks_["l1_a"]("I, a") * nuclear;

        // set doubles
        sigvec_blks_["sigmar2_aaa"]("I, a, b, i") += evec_blks_["r2_aaa"]("I, a, b, i") * nuclear;
        sigvec_blks_["sigmal2_aaa"]("I, a, b, i") += evec_blks_["l2_aaa"]("I, i, a, b") * nuclear;

        sigvec_blks_["sigmar2_abb"]("I, a, b, i") += evec_blks_["r2_abb"]("I, a, b, i") * nuclear;
        sigvec_blks_["sigmal2_abb"]("I, a, b, i") += evec_blks_["l2_bab"]("I, i, a, b") * nuclear;

        // clear temporary arrays
        tmps_.clear();

        world_.gop.fence();
    }

    void EOM_EA_CCSD::pack_sigma_vectors(size_t L, double **sigmar, double **sigmal) {
        size_t oa = oa_, ob = ob_, // number of occupied orbitals
        va = va_, vb = vb_; // number of virtual beta orbitals

        size_t aaa = (va*va-va)*oa/2, abb = va*vb*ob;

        /// unpack electronic part of right sigma vectors
        size_t offset = 0;

        // singles
        {
            foreach_inplace(sigvec_blks_["sigmar1_a"], [sigmar, oa, offset](auto &tile) {
                for (auto &x : tile.range()) {
                    size_t trial = x[0], e = x[1];
                    sigmar[e + offset][trial] = tile[x];
                }
            }); offset += va;
        }

        // doubles
        {
            foreach_inplace(sigvec_blks_["sigmar2_aaa"], [sigmar, oa, va, ob, vb, offset, this](auto &tile) {
                for (auto &x : tile.range()) {
                    size_t trial = x[0], e = x[1], f = x[2], m = x[3];
                    if (e < f) {
                        size_t ef_off = EOM_Driver::sqr_2_tri_idx(e, f, va);
                        size_t sigmar_idx = ef_off * oa + m + offset;
                        sigmar[sigmar_idx][trial] = tile[x];
                    }
                }
            }); offset += aaa;

            foreach_inplace(sigvec_blks_["sigmar2_abb"], [sigmar, oa, ob, vb, offset](auto &tile) {
                for (auto &x : tile.range()) {
                    size_t trial = x[0], e = x[1], f = x[2], m = x[3];
                    sigmar[e * vb * ob + f * ob + m + offset][trial] = tile[x];
                }
            }); offset += abb;
        }

        /// unpack electronic part of left sigma vectors
        offset = 0;
        // singles
        {
            foreach_inplace(sigvec_blks_["sigmal1_a"], [sigmal, oa, offset](auto &tile) {
                for (auto &x : tile.range()) {
                    size_t trial = x[0], e = x[1];
                    sigmal[e + offset][trial] = tile[x];
                }
            }); offset += va;
        }

        // doubles
        {
            foreach_inplace(sigvec_blks_["sigmal2_aaa"], [sigmal, oa, ob, va, vb, offset, this](auto &tile) {
                for (auto &x : tile.range()) {
                    size_t trial = x[0], e = x[1], f = x[2], m = x[3];
                    if (e < f) {
                        size_t ef_off = EOM_Driver::sqr_2_tri_idx(e, f, va);
                        size_t sigmal_idx = ef_off * oa + m + offset;
                        sigmal[sigmal_idx][trial] = tile[x];
                    }
                }
            }); offset += aaa;

            foreach_inplace(sigvec_blks_["sigmal2_abb"], [sigmal, oa, ob, vb, offset](auto &tile) {
                for (auto &x : tile.range()) {
                    size_t trial = x[0], e = x[1], f = x[2], m = x[3];
                    sigmal[e * vb * ob + f * ob + m + offset][trial] = tile[x];
                }
            }); offset += abb;
        }

    }

    double *EOM_EA_CCSD::get_state_norms(size_t i) const {
        // get pointers to the eigenvectors
        double** rerp_ = revec_->pointer();
        double** relp_ = levec_->pointer();

        // allocate memory for the norms
        auto *norms = (double *) calloc(nops_, sizeof(double));

        size_t nid = 0; // norm id

        { // singles
            norms[nid] = C_DDOT(singleDim_, rerp_[i], 1, relp_[i], 1);
            norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0;
            nid++;
        }
        { // doubles
            norms[nid] = C_DDOT(doubleDim_, rerp_[i] + singleDim_, 1, relp_[i] + singleDim_, 1);
            norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0;
            nid++;
        }

        return norms;
    }

    void EOM_EA_CCSD::unpack_eigenvectors() {
        return;
    }

    EOM_Driver::DominantTransitionsType EOM_EA_CCSD::find_dominant_transitions(size_t I) {
        /// get dominant transition for each root in each block

        // get pointers to eigenvectors
        double **rerp = revec_->pointer();
        double **relp = levec_->pointer();

        /// get dominant transitions for each state in each block
        static constexpr double threshold = 1e-10;

        // a map of the operator name to a priority queue of its dominant transitions with:
        //     the magnitude of the transition
        //     the spin of the transition
        //     the indicies of the transition
        EOM_Driver::DominantTransitionsType dominant_transitions;

        size_t off = 0;
        size_t id = 0;

        // singles a
        double l, r, lr;
        for (size_t a = 0; a < va_; a++) {
            l = relp[I][id + off]; // get l
            r = rerp[I][id + off]; // get r
            lr = l*r; // get lr
            if (fabs(lr) > threshold) {
                dominant_transitions["l1*r1"].push({lr, l, r, "a", {a+1 + oa_}});
            }
            id++;
        }
        off += id;

        // doubles aaa
        for (size_t a = 0; a < va_; a++) {
            for (size_t b = a + 1; b < va_; b++) {
                for (size_t i = 0; i < oa_; i++) {
                    l = relp[I][id + off]; // get l
                    r = rerp[I][id + off]; // get r
                    lr = l*r; // get lr
                    if (fabs(lr) > threshold) {
                        dominant_transitions["l2*r2"].push({lr, l, r, "aaa", {a+1 + oa_, b+1 + oa_, oa_ - i}});
                    }
                    id++;
                }
            }
        }
        off += id;
        id = 0;

        // doubles abb
        for (size_t a = 0; a < va_; a++) {
            for (size_t b = 0; b < vb_; b++) {
                for (size_t i = 0; i < oa_; i++) {
                    l = relp[I][id + off]; // get l
                    r = rerp[I][id + off]; // get r
                    lr = l*r; // get lr
                    if (fabs(lr) > threshold) {
                        dominant_transitions["l2*r2"].push({lr, l, r, "abb", {a+1 + oa_, b+1 + ob_, ob_ - i}});
                    }
                    id++;
                }
            }
        }

        return dominant_transitions;
    }
} // cc_cavity