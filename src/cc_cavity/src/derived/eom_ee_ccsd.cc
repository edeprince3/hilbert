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

#include "cc_cavity/include/derived/eom_ee_ccsd.h"
#include <memory>

namespace hilbert {

    EOM_EE_CCSD::EOM_EE_CCSD(const std::shared_ptr<CC_Cavity> &cc_wfn, Options &options) :
            EOM_Driver(cc_wfn, options) {}

    void EOM_EE_CCSD::set_problem_size() {

        // count number of operators
        nops_ = 3; // r0, r1, r2

        // compute the number of amplitudes
        singleDim_ = oa_*va_ + ob_*vb_;
        doubleDim_ = oa_*(oa_-1)/2*va_*(va_-1)/2 + ob_*(ob_-1)/2*vb_*(vb_-1)/2 + oa_*ob_*va_*vb_;

        dim_e_ = 1 + singleDim_ + doubleDim_;
        dim_p_ = 0;

        N_ = dim_e_;

    }

    void EOM_EE_CCSD::build_hamiltonian() {
        throw PsiException("EOM_EE_CCSD::build_hamiltonian not implemented", __FILE__, __LINE__);
    }

    double* EOM_EE_CCSD::build_preconditioner() {

        // get properties from the CC_Cavity object
        double &cc_energy_ = cc_wfn_->cc_energy_;
        double const *cavity_frequency_ = cc_wfn_->cavity_frequency_;
        double &enuc_ = cc_wfn_->enuc_;
        double &average_electric_dipole_self_energy_ = cc_wfn_->average_electric_dipole_self_energy_;
        double const *epsilon_ = cc_wfn_->epsilon_;

        int oa = (int)oa_, ob = (int)ob_, // occupied alpha/beta
            va = (int)va_, vb = (int)vb_; // virtual alpha/beta

        // allocate memory for the preconditioner
        auto * precond = (double*) calloc(N_, sizeof(double));

        double* ss_diag;
        int id_s = 0;
        if (eom_ss_guess_) // if we want to use the singles block in the preconditioner
            ss_diag = build_ss_diagonal();

        /// reference
        precond[0] = cc_energy_;

        /// singles
        //alpha
        int id = 1;
        for (int a = 0; a < va; a++)
            for (int i = 0; i < oa; ++i) {
                if (eom_ss_guess_) precond[id++] = ss_diag[id_s++] + enuc_;
                else precond[id++] = cc_energy_ + epsilon_[a + o_] - epsilon_[i];
            }

        //beta
        for (int a = va; a < va + vb; a++)
            for (int i = oa; i < oa + ob; ++i) {
                if (eom_ss_guess_) precond[id++] = ss_diag[id_s++] + enuc_;
                else precond[id++] = cc_energy_ + epsilon_[a + o_] - epsilon_[i];
            }

        /// doubles
        //alpha/alpha
        for (int a = 0; a < va; a++)
            for (int b = a + 1; b < va; b++)
                for (int i = 0; i < oa; i++)
                    for (int j = i + 1; j < oa; ++j)
                        precond[id++] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[i] - epsilon_[j];

        //alpha/beta
        for (int a = 0; a < va; a++)
            for (int b = va; b < va + vb; b++)
                for (int i = 0; i < oa; i++)
                    for (int j = oa; j < oa + ob; ++j)
                        precond[id++] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[i] - epsilon_[j];

        //beta/beta
        for (int a = va; a < va + vb; a++)
            for (int b = a + 1; b < va + vb; b++)
                for (int i = oa; i < oa + ob; i++)
                    for (int j = i + 1; j < oa + ob; ++j)
                        precond[id++] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[i] - epsilon_[j];

        // free memory if needed
        if (eom_ss_guess_) free(ss_diag);

        // return the preconditioner
        return precond;
    }

    void EOM_EE_CCSD::print_eom_header() {
        // print out energies and amplitude norms
        Printf("\n%5s", "state");
        Printf(" %20s", "total energy (Eh)");
        Printf(" %17s", "ex. energy (Eh)");
        Printf(" %13s", "l0*r0");
        Printf(" %13s", "l1*r1");
        Printf(" %13s", "l2*r2");
        Printf("\n");

        // set reference energy
        ground_energy_ref_ = eigvals_->get(0);
    }

    double* EOM_EE_CCSD::get_state_norms(size_t i) const {

        // get pointers to the eigenvectors
        double** rerp_ = revec_->pointer();
        double** relp_ = levec_->pointer();

        // allocate memory for the norms
        auto *norms = (double *) calloc(nops_, sizeof(double));

        size_t nid = 0; // norm id

        { // reference
            norms[nid] = rerp_[i][0] * relp_[i][0];
            norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0;
            nid++;
        }
        { // singles
            norms[nid] = C_DDOT(singleDim_, rerp_[i] + 1, 1, relp_[i] + 1, 1);
            norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0;
            nid++;
        }
        { // doubles
            norms[nid] = C_DDOT(doubleDim_, rerp_[i] + singleDim_ + 1, 1, relp_[i] + singleDim_ + 1, 1);
            norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0;
            nid++;
        }

        return norms;
    };

    void EOM_EE_CCSD::unpack_trial_vectors(size_t L, double **Q) {
        /// unpack the trial vectors from Q into the appropriate blocks of the left/right eigenvectors
        
        evec_blks_.clear(); // clear the map of trial vectors
        
        /// initialize electronic operators
        evec_blks_["r0"] = makeTensor(world_, {L}, false); // reference
        evec_blks_["r1_aa"] = makeTensor(world_, {L, va_, oa_}, false); // alpha singles excitations
        evec_blks_["r1_bb"] = makeTensor(world_, {L, vb_, ob_}, false); // beta singles excitations
        evec_blks_["r2_aaaa"] = makeTensor(world_, {L, va_, va_, oa_, oa_}, false); // alpha doubles excitations
        evec_blks_["r2_abab"] = makeTensor(world_, {L, va_, vb_, oa_, ob_}, false); // alpha-beta doubles excitations
        evec_blks_["r2_bbbb"] = makeTensor(world_, {L, vb_, vb_, ob_, ob_}, false); // beta doubles excitations

        size_t oa = oa_, ob = ob_, // number of occupied orbitals
               va = va_, vb = vb_; // number of virtual beta orbitals
        size_t aa = va*oa, bb = vb*ob; // singles offsets
        size_t aaaa = (va*va-va)*(oa*oa-oa)/4, abab = va*vb*oa*ob, bbbb = (vb*vb-vb)*(ob*ob-ob)/4; // doubles offsets

        size_t offset;
        /// pack electronic trial vectors into TA::TArrayD
        // pack ground state operator
        evec_blks_["r0"].init_elements([Q](auto &I) {
            return Q[I[0]][0];
        }); offset = 1;

        // pack singles into tensors
        {
            evec_blks_["r1_aa"].init_elements([Q, oa, offset](auto &I) {
                return Q[I[0]][I[1]*oa + I[2] + offset];
            }); offset += aa;

            evec_blks_["r1_bb"].init_elements([Q, ob, offset](auto &I) {
                return Q[I[0]][I[1]*ob + I[2] + offset];
            }); offset += bb;
        }

        // pack redundant doubles into tensors
        {
            evec_blks_["r2_aaaa"].init_elements([Q, va, oa, offset](auto &I) {
                size_t a = I[1], b = I[2], i = I[3], j = I[4];
                if (a == b || i == j) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b) ^ (i > j);
                if (a > b) std::swap(a, b);
                if (i > j) std::swap(i, j);

                size_t ab_off = EOM_Driver::sqr_2_tri_idx(a, b, va);
                size_t ij_off = EOM_Driver::sqr_2_tri_idx(i, j, oa);
                size_t upper_tri_index = ab_off*oa*(oa-1)/2 + ij_off;

                double value = Q[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;


            }); offset += aaaa;

            evec_blks_["r2_abab"].init_elements([Q, oa, ob, vb, offset](auto &I) {
                return Q[I[0]][I[1]*vb*oa*ob + I[2]*oa*ob + I[3]*ob + I[4] + offset];
            }); offset += abab;

            evec_blks_["r2_bbbb"].init_elements([Q, vb, ob, offset](auto &I) {
                size_t a = I[1], b = I[2], i = I[3], j = I[4];
                if (a == b || i == j) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b) ^ (i > j);
                if (a > b) std::swap(a, b);
                if (i > j) std::swap(i, j);

                size_t ab_off = EOM_Driver::sqr_2_tri_idx(a, b, vb);
                size_t ij_off = EOM_Driver::sqr_2_tri_idx(i, j, ob);
                size_t upper_tri_index = ab_off*ob*(ob-1)/2 + ij_off;
                double value = Q[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;

            }); offset += bbbb;
        }
        world_.gop.fence();

        /// initialize left operators from right
        evec_blks_["l0"]("I") = evec_blks_["r0"]("I");
        evec_blks_["l1_aa"]("I, m, e") = evec_blks_["r1_aa"]("I, e, m");
        evec_blks_["l1_bb"]("I, m, e") = evec_blks_["r1_bb"]("I, e, m");
        evec_blks_["l2_aaaa"]("I, m, n, e, f") = evec_blks_["r2_aaaa"]("I, e, f, m, n");
        evec_blks_["l2_abab"]("I, m, n, e, f") = evec_blks_["r2_abab"]("I, e, f, m, n");
        evec_blks_["l2_bbbb"]("I, m, n, e, f") = evec_blks_["r2_bbbb"]("I, e, f, m, n");
        world_.gop.fence();
    }

    void EOM_EE_CCSD::pack_sigma_vectors(size_t L, double **sigmar, double **sigmal) {

        /// fill elements of all sigma vectors from their respective blocks

        size_t oa = oa_, ob = ob_, // number of occupied orbitals
        va = va_, vb = vb_; // number of virtual beta orbitals
        size_t aa = va*oa, bb = vb*ob; // singles offsets
        size_t aaaa = (va*va-va)*(oa*oa-oa)/4, abab = va*vb*oa*ob, bbbb = (vb*vb-vb)*(ob*ob-ob)/4; // doubles offsets

        /// unpack electronic part of right sigma vectors
        size_t offset = 0;
        // ground state
        {
            forall(sigvec_blks_["sigmar0"], [sigmar](auto &tile, auto &x) {
                sigmar[0][x[0]] = tile[x];
            }); offset += 1;
        }
        // singles
        {
            forall(sigvec_blks_["sigmar1_aa"], [sigmar, oa, offset](auto &tile, auto &x) {
                sigmar[x[1]*oa + x[2] + offset][x[0]] = tile[x];
            }); offset += aa;
            
            forall(sigvec_blks_["sigmar1_bb"], [sigmar, ob, offset](auto &tile, auto &x) {
                sigmar[x[1]*ob + x[2] + offset][x[0]] = tile[x];
            }); offset += bb;
        }

        // doubles
        {

            forall(sigvec_blks_["sigmar2_aaaa"], [sigmar, oa, va, offset](auto &tile, auto &x) {
                if (x[1] < x[2] && x[3] < x[4]) {
                    size_t ab_off = EOM_Driver::sqr_2_tri_idx(x[1], x[2], va);
                    size_t ij_off = EOM_Driver::sqr_2_tri_idx(x[3], x[4], oa);
                    size_t sigmar_idx = ab_off*oa*(oa-1)/2 + ij_off + offset;
                    sigmar[sigmar_idx][x[0]] = tile[x];
                }
            }); offset += aaaa;
            
            forall(sigvec_blks_["sigmar2_abab"], [sigmar, oa, ob, vb, offset](auto &tile, auto &x) {
                sigmar[x[1]*vb*oa*ob + x[2]*oa*ob + x[3]*ob + x[4] + offset][x[0]] = tile[x];
            }); offset += abab;
            
            forall(sigvec_blks_["sigmar2_bbbb"], [sigmar, ob, vb, offset](auto &tile, auto &x) {
                if (x[1] < x[2] && x[3] < x[4]) {
                    size_t ab_off = EOM_Driver::sqr_2_tri_idx(x[1], x[2], vb);
                    size_t ij_off = EOM_Driver::sqr_2_tri_idx(x[3], x[4], ob);
                    size_t sigmar_idx = ab_off*ob*(ob-1)/2 + ij_off + offset;
                    sigmar[sigmar_idx][x[0]] = tile[x];
                }
            }); offset += bbbb;
        }

        /// unpack electronic part of left sigma vectors
        offset = 0;
        // ground state
        {
            forall(sigvec_blks_["sigmal0"], [sigmal](auto &tile, auto &x) {
                sigmal[0][x[0]] = tile[x];
            }); offset += 1;
        }
        // singles
        {
            forall(sigvec_blks_["sigmal1_aa"], [sigmal, oa, offset](auto &tile, auto &x) {
                sigmal[x[1]*oa + x[2] + offset][x[0]] = tile[x];
            }); offset += aa;

            forall(sigvec_blks_["sigmal1_bb"], [sigmal, ob, offset](auto &tile, auto &x) {
                sigmal[x[1]*ob + x[2] + offset][x[0]] = tile[x];
            }); offset += bb;
        }

        // doubles
        {
            forall(sigvec_blks_["sigmal2_aaaa"], [sigmal, oa, va, offset](auto &tile, auto &x) {
                if (x[1] < x[2] && x[3] < x[4]) {
                    size_t ab_off = EOM_Driver::sqr_2_tri_idx(x[1], x[2], va);
                    size_t ij_off = EOM_Driver::sqr_2_tri_idx(x[3], x[4], oa);
                    size_t sigmal_idx = ab_off*oa*(oa-1)/2 + ij_off + offset;
                    sigmal[sigmal_idx][x[0]] = tile[x];
                }
            }); offset += aaaa;

            forall(sigvec_blks_["sigmal2_abab"], [sigmal, oa, ob, vb, offset](auto &tile, auto &x) {
                sigmal[x[1] * vb * oa * ob + x[2] * oa * ob + x[3] * ob + x[4] +
                       offset][x[0]] = tile[x];
            }); offset += abab;

            forall(sigvec_blks_["sigmal2_bbbb"], [sigmal, ob, vb, offset](auto &tile, auto &x) {
                if (x[1] < x[2] && x[3] < x[4]) {
                    size_t ab_off = EOM_Driver::sqr_2_tri_idx(x[1], x[2], vb);
                    size_t ij_off = EOM_Driver::sqr_2_tri_idx(x[3], x[4], ob);
                    size_t sigmal_idx = ab_off*ob*(ob-1)/2 + ij_off + offset;
                    sigmal[sigmal_idx][x[0]] = tile[x];
                }
            }); offset += bbbb;
        }
        world_.gop.fence();
    }

    void EOM_EE_CCSD::unpack_eigenvectors() {

        evec_blks_.clear(); // clear the map of trial vectors

        size_t L = M_;
        double **rerp = revec_->pointer(), **relp = levec_->pointer();
        
        // initialize operators
        evec_blks_["r0"] = makeTensor(world_, {L}, false); // reference
        evec_blks_["r1_aa"] = makeTensor(world_, {L, va_, oa_}, false); // alpha singles excitations
        evec_blks_["r1_bb"] = makeTensor(world_, {L, vb_, ob_}, false); // beta singles excitations
        evec_blks_["r2_aaaa"] = makeTensor(world_, {L, va_, va_, oa_, oa_}, false); // alpha doubles excitations
        evec_blks_["r2_abab"] = makeTensor(world_, {L, va_, vb_, oa_, ob_}, false); // alpha-beta doubles excitations
        evec_blks_["r2_bbbb"] = makeTensor(world_, {L, vb_, vb_, ob_, ob_}, false); // beta doubles excitations

        size_t oa = oa_, ob = ob_, // number of occupied orbitals
        va = va_, vb = vb_; // number of virtual beta orbitals
        size_t aa = va*oa, bb = vb*ob; // singles offsets
        size_t aaaa = va*oa*(oa-1)*(va-1)/4, abab = va*vb*oa*ob, bbbb = vb*ob*(ob-1)*(vb-1)/4; // doubles offsets

        size_t offset;
        
        /// pack electronic trial vectors into TA::TArrayD
        // pack ground state operator
        evec_blks_["r0"].init_elements([rerp](auto &I) {
            return rerp[I[0]][0];
        }); offset = 1;

        // pack singles into tensors
        {
            evec_blks_["r1_aa"].init_elements([rerp, oa, offset](auto &I) {
                return rerp[I[0]][I[1]*oa + I[2] + offset];
            }); offset += aa;

            evec_blks_["r1_bb"].init_elements([rerp, ob, offset](auto &I) {
                return rerp[I[0]][I[1]*ob + I[2] + offset];
            }); offset += bb;
        }

        // pack redundant doubles into tensors
        {
            evec_blks_["r2_aaaa"].init_elements([rerp, va, oa, offset](auto &I) {
                size_t a = I[1], b = I[2], i = I[3], j = I[4];
                if (a == b || i == j) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b) ^ (i > j);
                if (a > b) std::swap(a, b);
                if (i > j) std::swap(i, j);

                size_t ab_off = EOM_Driver::sqr_2_tri_idx(a, b, va);
                size_t ij_off = EOM_Driver::sqr_2_tri_idx(i, j, oa);
                size_t upper_tri_index = ab_off*oa*(oa-1)/2 + ij_off;
                double value = rerp[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;

            }); offset += aaaa;

            evec_blks_["r2_abab"].init_elements([rerp, oa, ob, vb, offset](auto &I) {
                return rerp[I[0]][I[1]*vb*oa*ob + I[2]*oa*ob + I[3]*ob + I[4] + offset];
            }); offset += abab;
            
            evec_blks_["r2_bbbb"].init_elements([rerp, vb, ob, offset](auto &I) {
                size_t a = I[1], b = I[2], i = I[3], j = I[4];
                if (a == b || i == j) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b) ^ (i > j);
                if (a > b) std::swap(a, b);
                if (i > j) std::swap(i, j);

                size_t ab_off = EOM_Driver::sqr_2_tri_idx(a, b, vb);
                size_t ij_off = EOM_Driver::sqr_2_tri_idx(i, j, ob);
                size_t upper_tri_index = ab_off*ob*(ob-1)/2 + ij_off;
                double value = rerp[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;

            }); offset += bbbb;
        }

        /// initialize left operators
        evec_blks_["l0"] = makeTensor(world_, {L}, false); // reference
        evec_blks_["l1_aa"] = makeTensor(world_, {L, oa_, va_}, false); // alpha singles excitations
        evec_blks_["l1_bb"] = makeTensor(world_, {L, ob_, vb_}, false); // beta singles excitations
        evec_blks_["l2_aaaa"] = makeTensor(world_, {L, oa_, oa_, va_, va_}, false); // alpha doubles excitations
        evec_blks_["l2_abab"] = makeTensor(world_, {L, oa_, ob_, va_, vb_}, false); // alpha-beta doubles excitations
        evec_blks_["l2_bbbb"] = makeTensor(world_, {L, ob_, ob_, vb_, vb_}, false); // beta doubles excitations
        
        /// pack electronic trial vectors into TA::TArrayD
        // pack ground state operator
        evec_blks_["l0"].init_elements([relp](auto &I) {
            return relp[I[0]][0];
        }); offset = 1;

        // pack singles into tensors
        {
            evec_blks_["l1_aa"].init_elements([relp, oa, offset](auto &I) {
                return relp[I[0]][I[2]*oa + I[1] + offset];
            }); offset += aa;

            evec_blks_["l1_bb"].init_elements([relp, ob, offset](auto &I) {
                return relp[I[0]][I[2]*ob + I[1] + offset];
            }); offset += bb;
        }

        // pack redundant doubles into tensors
        {
            evec_blks_["l2_aaaa"].init_elements([relp, va, oa, offset](auto &I) {
                size_t a = I[3], b = I[4], i = I[1], j = I[2];
                if (a == b || i == j) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b) ^ (i > j);
                if (a > b) std::swap(a, b);
                if (i > j) std::swap(i, j);

                size_t ab_off = EOM_Driver::sqr_2_tri_idx(a, b, va);
                size_t ij_off = EOM_Driver::sqr_2_tri_idx(i, j, oa);
                size_t upper_tri_index = ab_off*oa*(oa-1)/2 + ij_off;
                double value = relp[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;

            }); offset += aaaa;

            evec_blks_["l2_abab"].init_elements([relp, oa, ob, vb, offset](auto &I) {
                return relp[I[0]][I[3]*vb*oa*ob + I[4]*oa*ob + I[1]*ob + I[2] + offset];
            }); offset += abab;

            evec_blks_["l2_bbbb"].init_elements([relp, vb, ob, offset](auto &I) {
                size_t a = I[3], b = I[4], i = I[1], j = I[2];
                if (a == b || i == j) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b) ^ (i > j);
                if (a > b) std::swap(a, b);
                if (i > j) std::swap(i, j);

                size_t ab_off = EOM_Driver::sqr_2_tri_idx(a, b, vb);
                size_t ij_off = EOM_Driver::sqr_2_tri_idx(i, j, ob);
                size_t upper_tri_index = ab_off*ob*(ob-1)/2 + ij_off;
                double value = relp[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;

            }); offset += bbbb;
        }

        world_.gop.fence();
    }

    EOM_Driver::DominantTransitionsType EOM_EE_CCSD::find_dominant_transitions(size_t I) {
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
        double l = relp[I][0]; // get l
        double r = rerp[I][0]; // get r
        double lr = l*r; // get lr

        // ground state
        if (fabs(lr) > threshold) {
            dominant_transitions["l0*r0"].push({lr, l, r, "", {}});
        }
        off++;


        size_t id = 0;
        // singles aa
        for (size_t a = 0; a < va_; a++) {
            for (size_t i = 0; i < oa_; i++) {
                l = relp[I][id + off]; // get l
                r = rerp[I][id + off]; // get r
                lr = l*r; // get lr
                if (fabs(lr) > threshold) {
                    dominant_transitions["l1*r1"].push({lr, l, r, "a", {a + oa_, i}});
                }
                id++;
            }
        }
        off += id;
        id = 0;

        // singles bb
        for (size_t a = 0; a < vb_; a++) {
            for (size_t i = 0; i < ob_; i++) {
                l = relp[I][id + off]; // get l
                r = rerp[I][id + off]; // get r
                lr = l*r; // get lr
                if (fabs(lr) > threshold) {
                    dominant_transitions["l1*r1"].push({lr, l, r, "b", {a + ob_, i}});
                }
                id++;
            }
        }
        off += id;
        id = 0;

        // doubles aaaa
        for (size_t a = 0; a < va_; a++) {
            for (size_t b = a + 1; b < va_; b++) {
                for (size_t i = 0; i < oa_; i++) {
                    for (size_t j = i + 1; j < oa_; j++) {
                        l = relp[I][id + off]; // get l
                        r = rerp[I][id + off]; // get r
                        lr = l*r; // get lr
                        if (fabs(lr) > threshold) {
                            dominant_transitions["l2*r2"].push({lr, l, r, "aa", {a + oa_, b + oa_, i, j}});
                        }
                        id++;
                    }
                }
            }
        }
        off += id;
        id = 0;

        // doubles abab
        for (size_t a = 0; a < va_; a++) {
            for (size_t b = 0; b < vb_; b++) {
                for (size_t i = 0; i < oa_; i++) {
                    for (size_t j = 0; j < ob_; j++) {
                        l = relp[I][id + off]; // get l
                        r = rerp[I][id + off]; // get r
                        lr = l*r; // get lr
                        if (fabs(lr) > threshold) {
                            dominant_transitions["l2*r2"].push({lr, l, r, "ab", {a + oa_, b + ob_, i, j}});
                        }
                        id++;
                    }
                }
            }
        }
        off += id;
        id = 0;

        // doubles bbbb
        for (size_t a = 0; a < vb_; a++) {
            for (size_t b = a + 1; b < vb_; b++) {
                for (size_t i = 0; i < ob_; i++) {
                    for (size_t j = i + 1; j < ob_; j++) {
                        l = relp[I][id + off]; // get l
                        r = rerp[I][id + off]; // get r
                        lr = l*r; // get lr
                        if (fabs(lr) > threshold) {
                            dominant_transitions["l2*r2"].push({lr, l, r, "bb", {a + ob_, b + ob_, i, j}});
                        }
                        id++;
                    }
                }
            }
        }

        return dominant_transitions;
    }

} // cc_cavity
