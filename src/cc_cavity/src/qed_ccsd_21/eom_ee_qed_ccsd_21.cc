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

#include "cc_cavity/include/qed_ccsd_21/eom_ee_qed_ccsd_21.h"
#include "cc_cavity/include/qed_ccsd_21/qed_ccsd_21.h"

namespace hilbert {

    EOM_EE_QED_CCSD::EOM_EE_QED_CCSD(std::shared_ptr<CC_Cavity> &cc_wfn, Options &options) :
            EOM_EE_CCSD(cc_wfn, options) {}


    void EOM_EE_QED_CCSD::set_problem_size() {

        // get problem size from the EOM_EE_CCSD object
        EOM_EE_CCSD::set_problem_size();

        // count number of operators
        nops_ += 3; // r0_1, r1_1, r2_1

        // compute the number of photon amplitudes
        dim_p_  = 1; // r0_1
        dim_p_ += singleDim_; // r1_1
        dim_p_ += doubleDim_; // r2_1

        N_ = dim_e_ + dim_p_;

    }

    void EOM_EE_QED_CCSD::build_hamiltonian() {
        throw PsiException("EOM_EE_QEDCCSD::build_hamiltonian not implemented", __FILE__, __LINE__);
    }

    double* EOM_EE_QED_CCSD::build_preconditioner() {

        // get properties from the CC_Cavity object
        double cc_energy_ = cc_wfn_->cc_energy_;
        double enuc = cc_wfn_->enuc_;
        double dse = cc_wfn_->average_electric_dipole_self_energy_;
        double w0 = cc_wfn_->cavity_frequency_[2];
        double const *epsilon_ = cc_wfn_->epsilon_;

        int oa = (int)oa_, ob = (int)ob_, // occupied alpha/beta
        va = (int)va_, vb = (int)vb_; // virtual alpha/beta

        // allocate memory for the preconditioner
        auto * precond = (double*) calloc(N_, sizeof(double));

        double* ss_diag;
        if (eom_ss_guess_)
            ss_diag = build_ss_diagonal();

        // reference
        precond[0] = cc_energy_;

        // reference + hw
        int id_s = 0;
        precond[dim_e_] = cc_energy_ + w0;
        if (eom_ss_guess_)
            precond[dim_e_] = ss_diag[id_s++] + dse + enuc;

        // singles
        //alpha
        int id = 1;
        for (int a = 0; a < va; a++)
            for (int i = 0; i < oa; ++i) {
                if (eom_ss_guess_) {
                    precond[id] = ss_diag[id_s] + dse + enuc;
                    precond[id + dim_e_] = ss_diag[id_s + singleDim_] + dse + enuc;
                } else {
                    precond[id] = cc_energy_ + epsilon_[a + o_] - epsilon_[i];
                    precond[id + dim_e_] = cc_energy_ + epsilon_[a + o_] - epsilon_[i] + w0;
                } id++; id_s++;
            }

        //beta
        for (int a = va; a < va + vb; a++)
            for (int i = oa; i < oa + ob; ++i) {
                if (eom_ss_guess_) {
                    precond[id] = ss_diag[id_s] + dse + enuc;
                    precond[id + dim_e_] = ss_diag[id_s + singleDim_] + dse + enuc;
                } else {
                    precond[id] = cc_energy_ + epsilon_[a + o_] - epsilon_[i];
                    precond[id + dim_e_] = cc_energy_ + epsilon_[a + o_] - epsilon_[i] + w0;
                } id++; id_s++;
            }

        // doubles
        //alpha/alpha
        for (int a = 0; a < va; a++)
            for (int b = a + 1; b < va; b++)
                for (int i = 0; i < oa; i++)
                    for (int j = i + 1; j < oa; ++j) {
                        precond[id] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[i] - epsilon_[j];
                        precond[id + dim_e_] = cc_energy_
                                               + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[i] - epsilon_[j]
                                               + w0;
                        id++;
                    }

        //alpha/beta
        for (int a = 0; a < va; a++)
            for (int b = va; b < va + vb; b++)
                for (int i = 0; i < oa; i++)
                    for (int j = oa; j < oa + ob; ++j) {
                        precond[id] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[i] - epsilon_[j];
                        precond[id + dim_e_] = cc_energy_
                                               + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[i] - epsilon_[j]
                                               + w0;
                        id++;
                    }

        //beta/beta
        for (int a = va; a < va + vb; a++)
            for (int b = a + 1; b < va + vb; b++)
                for (int i = oa; i < oa + ob; i++)
                    for (int j = i + 1; j < oa + ob; ++j) {
                        precond[id] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[i] - epsilon_[j];
                        precond[id + dim_e_] = cc_energy_
                                               + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[i] - epsilon_[j]
                                               + w0;
                        id++;
                    }

        // free memory if needed
        if (eom_ss_guess_) free(ss_diag);

        // return the preconditioner
        return precond;
    }

    void EOM_EE_QED_CCSD::print_eom_header() {
        // print out energies and amplitude norms
        Printf("\n%5s", "state");
        Printf(" %20s", "total energy (Eh)");
        Printf(" %17s", "ex. energy (Eh)");
        Printf(" %13s", "l0*r0");
        Printf(" %13s", "l1*r1");
        Printf(" %13s", "l2*r2");
        Printf(" %13s", "l0,1*r0,1");
        Printf(" %13s", "l1,1*r1,1");
        Printf(" %13s", "l2,1*r2,1");
        Printf("\n");

        // set reference energy
//        ground_energy_ref_ = eigvals_->get(0);
    }

    double* EOM_EE_QED_CCSD::get_state_norms(size_t i) const {

        // get pointers to the eigenvectors
        double** rerp_ = revec_->pointer();
        double** relp_ = levec_->pointer();

        // grab the norms from super class (should have correct size since nops_ is set in this constructor)
        auto *norms = EOM_EE_CCSD::get_state_norms(i);

        size_t nid = 3; // norm id

        size_t off = dim_e_;
        norms[nid] = (rerp_[i][off]*relp_[i][off]);
        norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0;
        off++; nid++;
        
        norms[nid] = C_DDOT(singleDim_, rerp_[i] + off, 1, relp_[i] + off, 1);
        norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0;
        off += singleDim_; nid++;
        
        norms[nid] = C_DDOT(doubleDim_, rerp_[i] + off, 1, relp_[i] + off, 1);
        norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0;
        off += doubleDim_; nid++;

        return norms;
    };

    void EOM_EE_QED_CCSD::unpack_trial_vectors(size_t L, double **Q) {

        // initialize electronic operators
        EOM_EE_CCSD::unpack_trial_vectors(L, Q);

        /// initialize photon operators
        evec_blks_["r0_1"] = makeTensor(world_, {L}, false); // reference + hw
        evec_blks_["r1_1_aa"] = makeTensor(world_, {L, va_, oa_}, false); // alpha singles excitations + hw
        evec_blks_["r1_1_bb"] = makeTensor(world_, {L, vb_, ob_}, false); // beta singles excitations + hw
        evec_blks_["r2_1_aaaa"] = makeTensor(world_, {L, va_, va_, oa_, oa_}, false); // alpha doubles excitations + hw
        evec_blks_["r2_1_abab"] = makeTensor(world_, {L, va_, vb_, oa_, ob_}, false); // alpha-beta doubles excitations + hw
        evec_blks_["r2_1_bbbb"] = makeTensor(world_, {L, vb_, vb_, ob_, ob_}, false); // beta doubles excitations + hw

        size_t oa = oa_, ob = ob_, // number of occupied orbitals
        va = va_, vb = vb_; // number of virtual beta orbitals
        size_t aa = va*oa, bb = vb*ob; // singles offsets
        size_t aaaa = (va*va-va)*(oa*oa-oa)/4, abab = va*vb*oa*ob, bbbb = (vb*vb-vb)*(ob*ob-ob)/4; // doubles offsets

        size_t offset = dim_e_;

        /// pack photonic trial vectors into TA::TArrayD
        // pack ground state + hw operator
        evec_blks_["r0_1"].init_elements([Q, offset](auto &I) {
            return Q[I[0]][offset];
        }); offset++;
        
        // pack singles + hw operator
        evec_blks_["r1_1_aa"].init_elements([Q, oa, offset](auto &I) {
            return Q[I[0]][I[1]*oa + I[2] + offset];
        }); offset += aa;

        evec_blks_["r1_1_bb"].init_elements([Q, ob, offset](auto &I) {
            return Q[I[0]][I[1]*ob + I[2] + offset];
        }); offset += bb;
        
        // pack redundant doubles + hw operator
        evec_blks_["r2_1_aaaa"].init_elements([Q, va, oa, offset](auto &I) {
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

        evec_blks_["r2_1_abab"].init_elements([Q, oa, ob, vb, offset](auto &I) {
            return Q[I[0]][I[1]*vb*oa*ob + I[2]*oa*ob + I[3]*ob + I[4] + offset];
        }); offset += abab;

        evec_blks_["r2_1_bbbb"].init_elements([Q, vb, ob, offset](auto &I) {
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
        world_.gop.fence();

        /// initialize left operators from right
        evec_blks_["l0_1"]("I") = evec_blks_["r0_1"]("I");
        evec_blks_["l1_1_aa"]("I, m, e") = evec_blks_["r1_1_aa"]("I, e, m");
        evec_blks_["l1_1_bb"]("I, m, e") = evec_blks_["r1_1_bb"]("I, e, m");
        evec_blks_["l2_1_aaaa"]("I, m, n, e, f") = evec_blks_["r2_1_aaaa"]("I, e, f, m, n");
        evec_blks_["l2_1_abab"]("I, m, n, e, f") = evec_blks_["r2_1_abab"]("I, e, f, m, n");
        evec_blks_["l2_1_bbbb"]("I, m, n, e, f") = evec_blks_["r2_1_bbbb"]("I, e, f, m, n");
        world_.gop.fence();
    }

    void EOM_EE_QED_CCSD::build_Hc_cH(size_t L) {

        // transform integrals if needed with t1 amplitudes
        if (!cc_wfn_->has_t1_integrals_)
            cc_wfn_->transform_integrals(true);

        /// initialize sigma vectors for this trial
        sigvec_blks_.clear();

        // reduce sigma into spins
        sigvec_blks_["sigmar0"]   = makeTensor(world_, {L}, true);
        sigvec_blks_["sigmar0_1"] = makeTensor(world_, {L}, true);

        sigvec_blks_["sigmal0"]   = makeTensor(world_, {L}, true);
        sigvec_blks_["sigmal0_1"] = makeTensor(world_, {L}, true);

        sigvec_blks_["sigmar1_aa"]   = makeTensor(world_, {L, va_, oa_}, true);
        sigvec_blks_["sigmar1_bb"]   = makeTensor(world_, {L, vb_, ob_}, true);
        sigvec_blks_["sigmar1_1_aa"] = makeTensor(world_, {L, va_, oa_}, true);
        sigvec_blks_["sigmar1_1_bb"] = makeTensor(world_, {L, vb_, ob_}, true);

        sigvec_blks_["sigmal1_aa"]   = makeTensor(world_, {L, va_, oa_}, true);
        sigvec_blks_["sigmal1_bb"]   = makeTensor(world_, {L, vb_, ob_}, true);
        sigvec_blks_["sigmal1_1_aa"] = makeTensor(world_, {L, va_, oa_}, true);
        sigvec_blks_["sigmal1_1_bb"] = makeTensor(world_, {L, vb_, ob_}, true);

        sigvec_blks_["sigmar2_aaaa"]   = makeTensor(world_, {L, va_, va_, oa_, oa_}, true);
        sigvec_blks_["sigmar2_abab"]   = makeTensor(world_, {L, va_, vb_, oa_, ob_}, true);
        sigvec_blks_["sigmar2_bbbb"]   = makeTensor(world_, {L, vb_, vb_, ob_, ob_}, true);
        sigvec_blks_["sigmar2_1_aaaa"] = makeTensor(world_, {L, va_, va_, oa_, oa_}, true);
        sigvec_blks_["sigmar2_1_abab"] = makeTensor(world_, {L, va_, vb_, oa_, ob_}, true);
        sigvec_blks_["sigmar2_1_bbbb"] = makeTensor(world_, {L, vb_, vb_, ob_, ob_}, true);

        sigvec_blks_["sigmal2_aaaa"]   = makeTensor(world_, {L, va_, va_, oa_, oa_}, true);
        sigvec_blks_["sigmal2_abab"]   = makeTensor(world_, {L, va_, vb_, oa_, ob_}, true);
        sigvec_blks_["sigmal2_bbbb"]   = makeTensor(world_, {L, vb_, vb_, ob_, ob_}, true);
        sigvec_blks_["sigmal2_1_aaaa"] = makeTensor(world_, {L, va_, va_, oa_, oa_}, true);
        sigvec_blks_["sigmal2_1_abab"] = makeTensor(world_, {L, va_, vb_, oa_, ob_}, true);
        sigvec_blks_["sigmal2_1_bbbb"] = makeTensor(world_, {L, vb_, vb_, ob_, ob_}, true);


        // initialize containers for coherent state basis terms
        for (auto &[name, resid] : sigvec_blks_) {
            tmps_["c" + name] = resid.clone();
            cc_wfn_->zero_tiles(tmps_["c" + name]);
        }

        // run sigma vector builds
        sigma_ee_21_1();
        sigma_ee_21_2();
        sigma_ee_21_3();
        sigma_ee_21_4();
        sigma_ee_21_5();
        sigma_ee_21_6();
        sigma_ee_21_7();
        sigma_ee_21_8();
        world_.gop.fence();

        /// add coherent state basis terms

        // Get cavity information
        double w0 = cc_wfn_->cavity_frequency_[2];
        double coupling_factor_z = w0 * cc_wfn_->cavity_coupling_strength_[2];

        double coherent_scalar;
        double e_dip_z = cc_wfn_->e_dip_z_;
        double nuc_dip_z = cc_wfn_->nuc_dip_z_;
        if ( options_.get_bool("QED_USE_RELAXED_ORBITALS")) {
            coherent_scalar = coupling_factor_z * e_dip_z;
        } else {
            coherent_scalar = -coupling_factor_z * nuc_dip_z;
        }

        // update sigma vectors with coherent state basis terms
        for (auto &[name, resid] : sigvec_blks_) {
            string idx = CC_Cavity::get_index(resid);
            resid(idx) += coherent_scalar * tmps_["c" + name](idx);
        }

        /// add constant terms to sigma vectors
        double nuclear = cc_wfn_->enuc_ + cc_wfn_->average_electric_dipole_self_energy_;

        // set ground state
        sigvec_blks_["sigmar0"]("I") = evec_blks_["r0"]("I") * cc_wfn_->cc_energy_;
        sigvec_blks_["sigmal0"]("I") = evec_blks_["l0"]("I") * cc_wfn_->cc_energy_;

        // set singles
        sigvec_blks_["sigmar1_aa"]("I, e, m") += evec_blks_["r1_aa"]("I, e, m") * nuclear;
        sigvec_blks_["sigmar1_bb"]("I, e, m") += evec_blks_["r1_bb"]("I, e, m") * nuclear;
        sigvec_blks_["sigmal1_aa"]("I, e, m") += evec_blks_["l1_aa"]("I, m, e") * nuclear;
        sigvec_blks_["sigmal1_bb"]("I, e, m") += evec_blks_["l1_bb"]("I, m, e") * nuclear;

        // set doubles
        sigvec_blks_["sigmar2_aaaa"]("I, e, f, m, n") += evec_blks_["r2_aaaa"]("I, e, f, m, n") * nuclear;
        sigvec_blks_["sigmar2_abab"]("I, e, f, m, n") += evec_blks_["r2_abab"]("I, e, f, m, n") * nuclear;
        sigvec_blks_["sigmar2_bbbb"]("I, e, f, m, n") += evec_blks_["r2_bbbb"]("I, e, f, m, n") * nuclear;
        sigvec_blks_["sigmal2_aaaa"]("I, e, f, m, n") += evec_blks_["l2_aaaa"]("I, m, n, e, f") * nuclear;
        sigvec_blks_["sigmal2_abab"]("I, e, f, m, n") += evec_blks_["l2_abab"]("I, m, n, e, f") * nuclear;
        sigvec_blks_["sigmal2_bbbb"]("I, e, f, m, n") += evec_blks_["l2_bbbb"]("I, m, n, e, f") * nuclear;

        // set ground state + hw
        sigvec_blks_["sigmar0_1"]("I") += evec_blks_["r0_1"]("I") * nuclear;
        sigvec_blks_["sigmal0_1"]("I") += evec_blks_["l0_1"]("I") * nuclear;

        // set singles + hw
        sigvec_blks_["sigmar1_1_aa"]("I, e, m") += evec_blks_["r1_1_aa"]("I, e, m") * nuclear;
        sigvec_blks_["sigmar1_1_bb"]("I, e, m") += evec_blks_["r1_1_bb"]("I, e, m") * nuclear;
        sigvec_blks_["sigmal1_1_aa"]("I, e, m") += evec_blks_["l1_1_aa"]("I, m, e") * nuclear;
        sigvec_blks_["sigmal1_1_bb"]("I, e, m") += evec_blks_["l1_1_bb"]("I, m, e") * nuclear;

        // set doubles + hw
        sigvec_blks_["sigmar2_1_aaaa"]("I, e, f, m, n") += evec_blks_["r2_1_aaaa"]("I, e, f, m, n") * nuclear;
        sigvec_blks_["sigmar2_1_abab"]("I, e, f, m, n") += evec_blks_["r2_1_abab"]("I, e, f, m, n") * nuclear;
        sigvec_blks_["sigmar2_1_bbbb"]("I, e, f, m, n") += evec_blks_["r2_1_bbbb"]("I, e, f, m, n") * nuclear;
        sigvec_blks_["sigmal2_1_aaaa"]("I, e, f, m, n") += evec_blks_["l2_1_aaaa"]("I, m, n, e, f") * nuclear;
        sigvec_blks_["sigmal2_1_abab"]("I, e, f, m, n") += evec_blks_["l2_1_abab"]("I, m, n, e, f") * nuclear;
        sigvec_blks_["sigmal2_1_bbbb"]("I, e, f, m, n") += evec_blks_["l2_1_bbbb"]("I, m, n, e, f") * nuclear;

        // clear temporary arrays
        tmps_.clear();

        world_.gop.fence();
    }

    void EOM_EE_QED_CCSD::pack_sigma_vectors(size_t L, double **sigmar, double **sigmal) {

        /// initialize electronic sigma vectors
        EOM_EE_CCSD::pack_sigma_vectors(L, sigmar, sigmal);

        /// fill elements of all sigma vectors from their respective blocks
        size_t oa = oa_, ob = ob_, // number of occupied orbitals
        va = va_, vb = vb_; // number of virtual beta orbitals
        size_t aa = va*oa, bb = vb*ob; // singles offsets
        size_t aaaa = (va*va-va)*(oa*oa-oa)/4, abab = va*vb*oa*ob, bbbb = (vb*vb-vb)*(ob*ob-ob)/4; // doubles offsets


        /// unpack photonic part of right sigma vectors
        size_t offset = dim_e_;

        // ground state + hw
        foreach_inplace(sigvec_blks_["sigmar0_1"], [sigmar, offset](auto &tile) {
            for (auto &x : tile.range())
                sigmar[offset][x[0]] = tile[x];
        }); offset += 1;

        // singles + hw
        foreach_inplace(sigvec_blks_["sigmar1_1_aa"], [sigmar, oa, offset](auto &tile) {
            for (auto &x : tile.range())
                sigmar[x[1] * oa + x[2] + offset][x[0]] = tile[x];
        }); offset += aa;

        foreach_inplace(sigvec_blks_["sigmar1_1_bb"], [sigmar, ob, offset](auto &tile) {
            for (auto &x : tile.range())
                sigmar[x[1] * ob + x[2] + offset][x[0]] = tile[x];
        }); offset += bb;

        foreach_inplace(sigvec_blks_["sigmar2_1_aaaa"], [sigmar, oa, va, offset](auto &tile) {
            for (auto &x : tile.range()) {
                if (x[1] < x[2] && x[3] < x[4]) {
                    size_t ab_off = EOM_Driver::sqr_2_tri_idx(x[1], x[2], va);
                    size_t ij_off = EOM_Driver::sqr_2_tri_idx(x[3], x[4], oa);
                    size_t sigmar_idx = ab_off * (oa*oa - oa)/2 + ij_off + offset;
                    sigmar[sigmar_idx][x[0]] = tile[x];
                }
            }
        }); offset += aaaa;

        foreach_inplace(sigvec_blks_["sigmar2_1_abab"], [sigmar, oa, ob, vb, offset](auto &tile) {
            for (auto &x : tile.range())
                sigmar[x[1]*vb*oa*ob + x[2]*oa*ob + x[3]*ob + x[4] + offset][x[0]] = tile[x];
        }); offset += abab;

        foreach_inplace(sigvec_blks_["sigmar2_1_bbbb"], [sigmar, ob, vb, offset](auto &tile) {
            for (auto &x : tile.range()) {
                if (x[1] < x[2] && x[3] < x[4]) {
                    size_t ab_off = EOM_Driver::sqr_2_tri_idx(x[1], x[2], vb);
                    size_t ij_off = EOM_Driver::sqr_2_tri_idx(x[3], x[4], ob);
                    size_t sigmar_idx = ab_off * (ob*ob - ob)/2 + ij_off + offset;
                    sigmar[sigmar_idx][x[0]] = tile[x];
                }
            }
        }); offset += bbbb;

        /// unpack photonic part of left sigma vectors
        offset = dim_e_;
        // ground state + hw
        foreach_inplace(sigvec_blks_["sigmal0_1"], [sigmal, offset](auto &tile) {
            for (auto &x : tile.range())
                sigmal[offset][x[0]] = tile[x];
        }); offset += 1;

        // singles + hw
        foreach_inplace(sigvec_blks_["sigmal1_1_aa"], [sigmal, oa, offset](auto &tile) {
            for (auto &x : tile.range())
                sigmal[x[1] * oa + x[2] + offset][x[0]] = tile[x];
        }); offset += aa;

        foreach_inplace(sigvec_blks_["sigmal1_1_bb"], [sigmal, ob, offset](auto &tile) {
            for (auto &x : tile.range())
                sigmal[x[1] * ob + x[2] + offset][x[0]] = tile[x];
        }); offset += bb;

        foreach_inplace(sigvec_blks_["sigmal2_1_aaaa"], [sigmal, oa, va, offset](auto &tile) {
            for (auto &x : tile.range()) {
                if (x[1] < x[2] && x[3] < x[4]) {
                    size_t ab_off = EOM_Driver::sqr_2_tri_idx(x[1], x[2], va);
                    size_t ij_off = EOM_Driver::sqr_2_tri_idx(x[3], x[4], oa);
                    size_t sigmal_idx = ab_off * (oa*oa - oa)/2 + ij_off + offset;
                    sigmal[sigmal_idx][x[0]] = tile[x];
                }
            }
        }); offset += aaaa;

        foreach_inplace(sigvec_blks_["sigmal2_1_abab"], [sigmal, oa, ob, vb, offset](auto &tile) {
            for (auto &x : tile.range())
                sigmal[x[1]*vb*oa*ob + x[2]*oa*ob + x[3]*ob + x[4] + offset][x[0]] = tile[x];
        }); offset += abab;

        foreach_inplace(sigvec_blks_["sigmal2_1_bbbb"], [sigmal, ob, vb, offset](auto &tile) {
            for (auto &x : tile.range()) {
                if (x[1] < x[2] && x[3] < x[4]) {
                    size_t ab_off = EOM_Driver::sqr_2_tri_idx(x[1], x[2], vb);
                    size_t ij_off = EOM_Driver::sqr_2_tri_idx(x[3], x[4], ob);
                    size_t sigmal_idx = ab_off * (ob*ob - ob)/2 + ij_off + offset;
                    sigmal[sigmal_idx][x[0]] = tile[x];
                }
            }
        }); offset += bbbb;
        world_.gop.fence();
    }

    void EOM_EE_QED_CCSD::unpack_eigenvectors() {

        // initialize electronic operators
        EOM_EE_CCSD::unpack_eigenvectors();

        size_t L = M_;
        double **rerp = revec_->pointer(), **relp = levec_->pointer();

        // initialize photonic operators
        evec_blks_["r0_1"] = makeTensor(world_, {L}, false); // reference + hw
        evec_blks_["r1_1_aa"] = makeTensor(world_, {L, va_, oa_}, false); // alpha singles excitations + hw
        evec_blks_["r1_1_bb"] = makeTensor(world_, {L, vb_, ob_}, false); // beta singles excitations + hw
        evec_blks_["r2_1_aaaa"] = makeTensor(world_, {L, va_, va_, oa_, oa_}, false); // alpha doubles excitations + hw
        evec_blks_["r2_1_abab"] = makeTensor(world_, {L, va_, vb_, oa_, ob_}, false); // alpha-beta doubles excitations + hw
        evec_blks_["r2_1_bbbb"] = makeTensor(world_, {L, vb_, vb_, ob_, ob_}, false); // beta doubles excitations + hw

        size_t oa = oa_, ob = ob_, // number of occupied orbitals
        va = va_, vb = vb_; // number of virtual beta orbitals
        size_t aa = va*oa, bb = vb*ob; // singles offsets
        size_t aaaa = (va*va-va)*(oa*oa-oa)/4, abab = va*vb*oa*ob, bbbb = (vb*vb-vb)*(ob*ob-ob)/4; // doubles offsets

        size_t offset = dim_e_;

        /// pack photonic eigenvectors into TA::TArrayD
        // pack ground state + hw operator
        evec_blks_["r0_1"].init_elements([rerp, offset](auto &I) {
            return rerp[I[0]][offset];
        }); offset++;

        // pack singles + hw operator
        evec_blks_["r1_1_aa"].init_elements([rerp, oa, offset](auto &I) {
            return rerp[I[0]][I[1]*oa + I[2] + offset];
        }); offset += aa;

        evec_blks_["r1_1_bb"].init_elements([rerp, ob, offset](auto &I) {
            return rerp[I[0]][I[1]*ob + I[2] + offset];
        }); offset += bb;

        // pack redundant doubles + hw operator
        evec_blks_["r2_1_aaaa"].init_elements([rerp, va, oa, offset](auto &I) {
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

        evec_blks_["r2_1_abab"].init_elements([rerp, oa, ob, vb, offset](auto &I) {
            return rerp[I[0]][I[1]*vb*oa*ob + I[2]*oa*ob + I[3]*ob + I[4] + offset];
        }); offset += abab;

        evec_blks_["r2_1_bbbb"].init_elements([rerp, vb, ob, offset](auto &I) {
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

        /// initialize left photonic operators
        evec_blks_["l0_1"] = makeTensor(world_, {L}, false); // reference + hw
        evec_blks_["l1_1_aa"] = makeTensor(world_, {L, oa_, va_}, false); // alpha singles excitations + hw
        evec_blks_["l1_1_bb"] = makeTensor(world_, {L, ob_, vb_}, false); // beta singles excitations + hw
        evec_blks_["l2_1_aaaa"] = makeTensor(world_, {L, oa_, oa_, va_, va_}, false); // alpha doubles excitations + hw
        evec_blks_["l2_1_abab"] = makeTensor(world_, {L, oa_, ob_, va_, vb_}, false); // alpha-beta doubles excitations + hw
        evec_blks_["l2_1_bbbb"] = makeTensor(world_, {L, ob_, ob_, vb_, vb_}, false); // beta doubles excitations + hw

        // reset offset
        offset = dim_e_;

        /// pack photonic trial vectors into TA::TArrayD
        // pack ground state + hw operator
        evec_blks_["l0_1"].init_elements([relp, offset](auto &I) {
            return relp[I[0]][offset];
        }); offset++;
        // pack singles + hw operator
        evec_blks_["l1_1_aa"].init_elements([relp, oa, offset](auto &I) {
            return relp[I[0]][I[2]*oa + I[1] + offset];
        }); offset += aa;

        evec_blks_["l1_1_bb"].init_elements([relp, ob, offset](auto &I) {
            return relp[I[0]][I[2]*ob + I[1] + offset];
        }); offset += bb;

        // pack redundant doubles + hw operator
        evec_blks_["l2_1_aaaa"].init_elements([relp, va, oa, offset](auto &I) {
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

        evec_blks_["l2_1_abab"].init_elements([relp, oa, ob, vb, offset](auto &I) {
            return relp[I[0]][I[3]*vb*oa*ob + I[4]*oa*ob + I[1]*ob + I[2] + offset];
        }); offset += abab;

        evec_blks_["l2_1_bbbb"].init_elements([relp, vb, ob, offset](auto &I) {
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

        world_.gop.fence();
    }

    EOM_Driver::DominantTransitionsType EOM_EE_QED_CCSD::find_dominant_transitions(size_t I) {
        /// get dominant transition for each root in each block

        // get pointers to eigenvectors
        double **rerp = revec_->pointer();
        double **relp = levec_->pointer();

        /// get dominant transitions for each state in each block
        static constexpr double threshold = 1e-10;

        // get electronic dominant transitions from EOM_EE_CCSD
        EOM_Driver::DominantTransitionsType dominant_transitions = EOM_EE_CCSD::find_dominant_transitions(I);

        size_t off = dim_e_;
        size_t id = 0;

        // ground state + hw
        double  l, r, lr;
        l = relp[I][id + off]; // get l
        r = rerp[I][id + off]; // get r
        lr = l*r; // get lr
        if (fabs(lr) > threshold) {
            dominant_transitions["l0_1*r0_1"].push({lr, l, r, "", {}});
        }
        off++;

        // singles aa + hw
        for (size_t a = 0; a < va_; a++) {
            for (size_t i = 0; i < oa_; i++) {
                l = relp[I][id + off]; // get l
                r = rerp[I][id + off]; // get r
                lr = l*r; // get lr
                if (fabs(lr) > threshold) {
                    dominant_transitions["l1_1*r1_1"].push({lr, l, r, "aa", {a+1 + oa_, oa_ - i}});
                }
                id++;
            }
        }
        off += id;
        id = 0;

        // singles bb + hw
        for (size_t a = 0; a < vb_; a++) {
            for (size_t i = 0; i < ob_; i++) {
                l = relp[I][id + off]; // get l
                r = rerp[I][id + off]; // get r
                lr = l*r; // get lr
                if (fabs(lr) > threshold) {
                    dominant_transitions["l1_1*r1_1"].push({lr, l, r, "bb", {a+1 + ob_, ob_ - i}});
                }
                id++;
            }
        }
        off += id;
        id = 0;

        // doubles aaaa + hw
        for (size_t a = 0; a < va_; a++) {
            for (size_t b = a + 1; b < va_; b++) {
                for (size_t i = 0; i < oa_; i++) {
                    for (size_t j = i + 1; j < oa_; j++) {
                        l = relp[I][id + off]; // get l
                        r = rerp[I][id + off]; // get r
                        lr = l*r; // get lr
                        if (fabs(lr) > threshold) {
                            dominant_transitions["l2_1*r2_1"].push({lr, l, r, "aaaa", {a+1 + oa_, b+1 + oa_, oa_ - i, oa_ - j}});
                        }
                        id++;
                    }
                }
            }
        }
        off += id;
        id = 0;

        // doubles abab + hw
        for (size_t a = 0; a < va_; a++) {
            for (size_t b = 0; b < vb_; b++) {
                for (size_t i = 0; i < oa_; i++) {
                    for (size_t j = 0; j < ob_; j++) {
                        l = relp[I][id + off]; // get l
                        r = rerp[I][id + off]; // get r
                        lr = l*r; // get lr
                        if (fabs(lr) > threshold) {
                            dominant_transitions["l2_1*r2_1"].push({lr, l, r, "abab", {a+1 + oa_, b+1 + ob_, oa_ - i, ob_ - j}});
                        }
                        id++;
                    }
                }
            }
        }
        off += id;
        id = 0;

        // doubles bbbb + hw
        for (size_t a = 0; a < vb_; a++) {
            for (size_t b = a + 1; b < vb_; b++) {
                for (size_t i = 0; i < ob_; i++) {
                    for (size_t j = i + 1; j < ob_; j++) {
                        l = relp[I][id + off]; // get l
                        r = rerp[I][id + off]; // get r
                        lr = l*r; // get lr
                        if (fabs(lr) > threshold) {
                            dominant_transitions["l2_1*r2_1"].push({lr, l, r, "bbbb", {a+1 + ob_, b+1 + ob_, ob_ - i, ob_ - j}});
                        }
                        id++;
                    }
                }
            }
        }

        return dominant_transitions;
    }

} // hilbert