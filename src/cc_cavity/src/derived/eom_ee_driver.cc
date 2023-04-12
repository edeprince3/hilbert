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

#include "../../include/derived/eom_ee_driver.h"
#include <memory>

namespace hilbert {

    EOM_EE_Driver::EOM_EE_Driver(const std::shared_ptr<CC_Cavity> &cc_wfn, Options &options) :
            EOM_Driver(cc_wfn, options) {}

    void EOM_EE_Driver::set_problem_size() {

        // compute the number of amplitudes
        singleDim_ = oa_*va_ + ob_*vb_;
        doubleDim_ = oa_*(oa_-1)/2*va_*(va_-1)/2 + ob_*(ob_-1)/2*vb_*(vb_-1)/2 + oa_*ob_*va_*vb_;

        dim_e_ = 1 + singleDim_ + doubleDim_;
        dim_p_ = 0;
        if (include_u0_) dim_p_ += 1;
        if (include_u1_) dim_p_ += singleDim_;
        if (include_u2_) dim_p_ += doubleDim_;

        N_ = dim_e_ + dim_p_;

    }

    void EOM_EE_Driver::build_hamiltonian() {

    }

    double* EOM_EE_Driver::build_guess() {

        // get properties from the CC_Cavity object
        double &cc_energy_ = cc_wfn_->cc_energy_;
        double const *cavity_frequency_ = cc_wfn_->cavity_frequency_;
        double &enuc_ = cc_wfn_->enuc_;
        double &average_electric_dipole_self_energy_ = cc_wfn_->average_electric_dipole_self_energy_;
        double const *epsilon_ = cc_wfn_->epsilon_;

        int oa = (int)oa_, ob = (int)ob_, // occupied alpha/beta
            va = (int)va_, vb = (int)vb_; // virtual alpha/beta

        auto * Hdiag = (double*) calloc(N_, sizeof(double));

        double* ss_diag = nullptr;
        if (eom_ss_guess_) ss_diag = build_diagonal();

        // reference
        Hdiag[0] = cc_energy_;

        // reference + hw
        int id_s = 0;
        if ( include_u0_ ){
            Hdiag[dim_e_] = cc_energy_ + cavity_frequency_[2];
            if (eom_ss_guess_) {
                Hdiag[dim_e_] = ss_diag[id_s++] + enuc_;
                Hdiag[dim_e_] += average_electric_dipole_self_energy_;
            }
        }

        // singles
        //alpha
        int id = 1;
        for (int a = 0; a < va; a++) {
            for (int i = 0; i < oa; ++i) {
                if (eom_ss_guess_) {
                    Hdiag[id] = ss_diag[id_s];
                    Hdiag[id] += average_electric_dipole_self_energy_ + enuc_;

                    if (include_u1_) {
                        Hdiag[id + dim_e_] = ss_diag[id_s + singleDim_];
                        Hdiag[id + dim_e_] += average_electric_dipole_self_energy_ + enuc_;
                    }
                } else {
                    Hdiag[id] = cc_energy_ + epsilon_[a + o_] - epsilon_[i];
                    if (include_u1_) {
                        Hdiag[id + dim_e_] = cc_energy_ + epsilon_[a + o_] - epsilon_[i] + cavity_frequency_[2];
                    }

                }
                id++; id_s++;
            }
        }
        //beta
        for (int a = va; a < va + vb; a++) {
            for (int i = oa; i < oa + ob; ++i) {
                if (eom_ss_guess_) {
                    Hdiag[id] = ss_diag[id_s];
                    Hdiag[id] += average_electric_dipole_self_energy_ + enuc_;
                    if (include_u1_) {
                        Hdiag[id + dim_e_] = ss_diag[id_s + singleDim_];
                        Hdiag[id + dim_e_] += average_electric_dipole_self_energy_ + enuc_;
                    }
                } else {
                    Hdiag[id] = cc_energy_ + epsilon_[a + o_] - epsilon_[i];
                    if (include_u1_) {
                        Hdiag[id + dim_e_] = cc_energy_ + epsilon_[a + o_] - epsilon_[i] + cavity_frequency_[2];
                    }
                }
                id++; id_s++;
            }
        }

        // doubles
        //alpha/alpha
        for (int a = 0; a < va; a++) {
            for (int b = a + 1; b < va; b++) {
                for (int i = 0; i < oa; i++) {
                    for (int j = i + 1; j < oa; ++j) {
                        Hdiag[id] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[i] - epsilon_[j];
                        if ( include_u2_ ){
                            Hdiag[id + dim_e_] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[i] - epsilon_[j] + cavity_frequency_[2];
                        }
                        id++;
                    }
                }
            }
        }
        //alpha/beta
        for (int a = 0; a < va; a++) {
            for (int b = va; b < va + vb; b++) {
                for (int i = 0; i < oa; i++) {
                    for (int j = oa; j < oa + ob; ++j) {
                        Hdiag[id] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[i] - epsilon_[j];
                        if ( include_u2_ ){
                            Hdiag[id + dim_e_] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[i] - epsilon_[j] + cavity_frequency_[2];
                        }
                        id++;
                    }
                }
            }
        }
        //beta/beta
        for (int a = va; a < va + vb; a++) {
            for (int b = a + 1; b < va + vb; b++) {
                for (int i = oa; i < oa + ob; i++) {
                    for (int j = i + 1; j < oa + ob; ++j) {
                        Hdiag[id] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[i] - epsilon_[j];
                        if ( include_u2_ ){
                            Hdiag[id + dim_e_] = cc_energy_ + epsilon_[a + o_] + epsilon_[b + o_] - epsilon_[i] - epsilon_[j] + cavity_frequency_[2];
                        }
                        id++;
                    }
                }
            }
        }

        if (eom_ss_guess_) { free(ss_diag); }
        return Hdiag;
    }

    void EOM_EE_Driver::print_eom_summary() const {
        // print out energies and amplitude norms
        Printf("\n%5s", "state");
        Printf(" %20s", "total energy (Eh)");
        Printf(" %17s", "ex. energy (Eh)");
        Printf(" %13s", "l0*r0");
        Printf(" %13s", "l1*r1");
        Printf(" %13s", "l2*r2");
        if (include_u0_) Printf(" %13s", "m0*s0");
        if (include_u1_) Printf(" %13s", "m1*s1");
        if (include_u2_) Printf(" %13s", "m2*s2");
        Printf("\n");


        // count number of operators
        int nid = 3; // r0, r1, r2
        if (include_u0_) nid++; // s0
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
                if (fabs(norms[j]) > 1e-10) {
                    Printf("%13.10lf ", norms[j]);
                } else {
                    Printf("------------- ");
                }
            }
            Printf("\n");
            free(norms);
        }
    }

    double* EOM_EE_Driver::get_state_norms(size_t i) const {
        // get number of amplitudes
        int numNorms = 3;
        if (include_u0_) numNorms++;
        if (include_u1_) numNorms++;
        if (include_u2_) numNorms++;

        double** rerp_ = revec_->pointer();
        double** relp_ = levec_->pointer();

        auto *norms = (double *) calloc(numNorms, sizeof(double));
        int nid = 0;

        norms[nid] = rerp_[i][0]*relp_[i][0];
        norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0;
        nid++;
        norms[nid] = C_DDOT(singleDim_, rerp_[i] + 1, 1, relp_[i] + 1, 1);
        norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0;
        nid++;
        norms[nid] = C_DDOT(doubleDim_, rerp_[i] + singleDim_ + 1, 1, relp_[i] + singleDim_ + 1, 1);
        norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0;
        nid++;
        if (include_u0_) {
            norms[nid] = (rerp_[i][dim_e_]*relp_[i][dim_e_]);
            norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0;
            nid++;
        }
        if (include_u1_) {
            norms[nid] = C_DDOT(singleDim_, rerp_[i] + dim_e_ + 1, 1, relp_[i] + dim_e_ + 1, 1);
            norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0;
            nid++;
        }
        if (include_u2_) {
            norms[nid] =
                    C_DDOT(doubleDim_, rerp_[i] + dim_e_ + singleDim_ + 1, 1, relp_[i] + dim_e_ + singleDim_ + 1, 1);
            norms[nid] = fabs(norms[nid]) > 1e-12 ? norms[nid] : 0.0;
            nid++;
        }

        return norms;
    };

    void EOM_EE_Driver::unpack_trial_vectors(size_t L, double **Q) {
        /// unpack the trial vectors from Q into the appropriate blocks of the left/right eigenvectors
        
        evec_blks_.clear(); // clear the map of trial vectors
        
        // initialize operators
        evec_blks_["r0"] = HelperD::makeTensor(world_, {L}, false); // reference
        evec_blks_["r1_aa"] = HelperD::makeTensor(world_, {L, va_, oa_}, false); // alpha singles excitations
        evec_blks_["r1_bb"] = HelperD::makeTensor(world_, {L, vb_, ob_}, false); // beta singles excitations
        evec_blks_["r2_aaaa"] = HelperD::makeTensor(world_, {L, va_, va_, oa_, oa_}, false); // alpha doubles excitations
        evec_blks_["r2_abab"] = HelperD::makeTensor(world_, {L, va_, vb_, oa_, ob_}, false); // alpha-beta doubles excitations
        evec_blks_["r2_bbbb"] = HelperD::makeTensor(world_, {L, vb_, vb_, ob_, ob_}, false); // beta doubles excitations

        if (include_u0_){
            evec_blks_["s0"] = HelperD::makeTensor(world_, {L}, false); // reference + hw
        }
        if (include_u1_){
            evec_blks_["s1_aa"] = HelperD::makeTensor(world_, {L, va_, oa_}, false); // alpha singles excitations + hw
            evec_blks_["s1_bb"] = HelperD::makeTensor(world_, {L, vb_, ob_}, false); // beta singles excitations + hw
        }
        if (include_u2_){
            evec_blks_["s2_aaaa"] = HelperD::makeTensor(world_, {L, va_, va_, oa_, oa_}, false); // alpha doubles excitations + hw
            evec_blks_["s2_abab"] = HelperD::makeTensor(world_, {L, va_, vb_, oa_, ob_}, false); // alpha-beta doubles excitations + hw
            evec_blks_["s2_bbbb"] = HelperD::makeTensor(world_, {L, vb_, vb_, ob_, ob_}, false); // beta doubles excitations + hw
        }

        size_t oa = oa_, ob = ob_, // number of occupied orbitals
               va = va_, vb = vb_; // number of virtual beta orbitals

        size_t offset;
        /// pack electronic trial vectors into TA::TArrayD
        // pack ground state operator
        evec_blks_["r0"].init_elements([Q](auto &I) {
//            return I[0] == 0 ? Q[I[0]][0] : 0.0;
            return Q[I[0]][0];
        }); offset = 1;

        // pack singles into tensors
        {
            evec_blks_["r1_aa"].init_elements([Q, oa, offset](auto &I) {
                return Q[I[0]][I[1]*oa + I[2] + offset];
            }); offset += oa*va;

            evec_blks_["r1_bb"].init_elements([Q, ob, offset](auto &I) {
                return Q[I[0]][I[1]*ob + I[2] + offset];
            }); offset += vb*ob;
        }

        // pack redundant doubles into tensors
        {
            evec_blks_["r2_aaaa"].init_elements([Q, va, oa, offset, this](auto &I) {
                size_t a = I[1], b = I[2], i = I[3], j = I[4];
                if (a == b || i == j) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b) ^ (i > j);
                if (a > b) std::swap(a, b);
                if (i > j) std::swap(i, j);

                size_t ab_off = this->sqr_2_tri_idx(a, b, va);
                size_t ij_off = this->sqr_2_tri_idx(i, j, oa);
                size_t upper_tri_index = ab_off*oa*(oa-1)/2 + ij_off;

                double value = Q[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;


            }); offset += va*oa*(oa-1)*(va-1)/4;

            evec_blks_["r2_abab"].init_elements([Q, oa, ob, vb, offset](auto &I) {
                return Q[I[0]][I[1]*vb*oa*ob + I[2]*oa*ob + I[3]*ob + I[4] + offset];
            }); offset += va*vb*oa*ob;

            evec_blks_["r2_bbbb"].init_elements([Q, vb, ob, offset, this](auto &I) {
                size_t a = I[1], b = I[2], i = I[3], j = I[4];
                if (a == b || i == j) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b) ^ (i > j);
                if (a > b) std::swap(a, b);
                if (i > j) std::swap(i, j);

                size_t ab_off = this->sqr_2_tri_idx(a, b, vb);
                size_t ij_off = this->sqr_2_tri_idx(i, j, ob);
                size_t upper_tri_index = ab_off*ob*(ob-1)/2 + ij_off;
                double value = Q[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;

            }); offset += vb*ob*(vb-1)*(ob-1)/4;
        }

        /// pack photonic trial vectors into TA::TArrayD
        if (include_u0_){
            // pack ground state + hw operator
            evec_blks_["s0"].init_elements([Q, offset](auto &I) {
                return Q[I[0]][offset];
            }); offset++;
        }
        if (include_u1_){
            // pack singles + hw operator
            evec_blks_["s1_aa"].init_elements([Q, oa, offset](auto &I) {
                return Q[I[0]][I[1]*oa + I[2] + offset];
            }); offset += oa*va;

            evec_blks_["s1_bb"].init_elements([Q, ob, offset](auto &I) {
                return Q[I[0]][I[1]*ob + I[2] + offset];
            }); offset += ob*vb;
        }
        if (include_u2_){
            // pack redundant doubles + hw operator
            evec_blks_["s2_aaaa"].init_elements([Q, va, oa, offset, this](auto &I) {
                size_t a = I[1], b = I[2], i = I[3], j = I[4];
                if (a == b || i == j) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b) ^ (i > j);
                if (a > b) std::swap(a, b);
                if (i > j) std::swap(i, j);

                size_t ab_off = this->sqr_2_tri_idx(a, b, va);
                size_t ij_off = this->sqr_2_tri_idx(i, j, oa);
                size_t upper_tri_index = ab_off*oa*(oa-1)/2 + ij_off;
                double value = Q[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;

            }); offset += va*oa*(oa-1)*(va-1)/4;

            evec_blks_["s2_abab"].init_elements([Q, oa, ob, vb, offset](auto &I) {
                return Q[I[0]][I[1]*vb*oa*ob + I[2]*oa*ob + I[3]*ob + I[4] + offset];
            }); offset += va*vb*oa*ob;

            evec_blks_["s2_bbbb"].init_elements([Q, vb, ob, offset, this](auto &I) {
                size_t a = I[1], b = I[2], i = I[3], j = I[4];
                if (a == b || i == j) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b) ^ (i > j);
                if (a > b) std::swap(a, b);
                if (i > j) std::swap(i, j);

                size_t ab_off = this->sqr_2_tri_idx(a, b, vb);
                size_t ij_off = this->sqr_2_tri_idx(i, j, ob);
                size_t upper_tri_index = ab_off*ob*(ob-1)/2 + ij_off;
                double value = Q[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;

            }); offset += vb*ob*(ob-1)*(vb-1)/4;
        }
        world_.gop.fence();

        /// initialize left operators from right
        evec_blks_["l0"]("I") = evec_blks_["r0"]("I");
        evec_blks_["l1_aa"]("I, m, e") = evec_blks_["r1_aa"]("I, e, m");
        evec_blks_["l1_bb"]("I, m, e") = evec_blks_["r1_bb"]("I, e, m");
        evec_blks_["l2_aaaa"]("I, m, n, e, f") = evec_blks_["r2_aaaa"]("I, e, f, m, n");
        evec_blks_["l2_abab"]("I, m, n, e, f") = evec_blks_["r2_abab"]("I, e, f, m, n");
        evec_blks_["l2_bbbb"]("I, m, n, e, f") = evec_blks_["r2_bbbb"]("I, e, f, m, n");
        
        if ( include_u0_ ) evec_blks_["m0"]("I") = evec_blks_["s0"]("I");
        if ( include_u1_ ){
            evec_blks_["m1_aa"]("I, m, e") = evec_blks_["s1_aa"]("I, e, m");
            evec_blks_["m1_bb"]("I, m, e") = evec_blks_["s1_bb"]("I, e, m");
        }
        if ( include_u2_ ){
            evec_blks_["m2_aaaa"]("I, m, n, e, f") = evec_blks_["s2_aaaa"]("I, e, f, m, n");
            evec_blks_["m2_abab"]("I, m, n, e, f") = evec_blks_["s2_abab"]("I, e, f, m, n");
            evec_blks_["m2_bbbb"]("I, m, n, e, f") = evec_blks_["s2_bbbb"]("I, e, f, m, n");
        }
        world_.gop.fence();
    }

    void EOM_EE_Driver::pack_sigma_vectors(size_t L, double **sigmar, double **sigmal) {

        /// fill elements of all sigma vectors from their respective blocks

        size_t oa = oa_, ob = ob_, // number of occupied orbitals
        va = va_, vb = vb_; // number of virtual beta orbitals

        double cc_energy = cc_wfn_->cc_energy_;
        /// unpack electronic part of right sigma vectors
        size_t offset = 0;
        // ground state
        {
            HelperD::forall(sigvec_blks_["sigmar_0"], [sigmar](auto &tile, auto &x) {
                sigmar[0][x[0]] = tile[x];
            }); offset += 1;
        }
        // singles
        {
            HelperD::forall(sigvec_blks_["sigmar_s_aa"], [sigmar, oa, offset](auto &tile, auto &x) {
                sigmar[x[1]*oa + x[2] + offset][x[0]] = tile[x];
            }); offset += oa*va;
            
            HelperD::forall(sigvec_blks_["sigmar_s_bb"], [sigmar, ob, offset](auto &tile, auto &x) {
                sigmar[x[1]*ob + x[2] + offset][x[0]] = tile[x];
            }); offset += ob*vb;
        }

        // doubles
        {

            HelperD::forall(sigvec_blks_["sigmar_d_aaaa"], [sigmar, oa, va, offset, this](auto &tile, auto &x) {
                if (x[1] < x[2] && x[3] < x[4]) {
                    size_t ab_off = this->sqr_2_tri_idx(x[1], x[2], va);
                    size_t ij_off = this->sqr_2_tri_idx(x[3], x[4], oa);
                    size_t sigmar_idx = ab_off*oa*(oa-1)/2 + ij_off + offset;
                    sigmar[sigmar_idx][x[0]] = tile[x];
                }
            }); offset += va*oa*(oa-1)*(va-1)/4;
            
            HelperD::forall(sigvec_blks_["sigmar_d_abab"], [sigmar, oa, ob, vb, offset](auto &tile, auto &x) {
                sigmar[x[1]*vb*oa*ob + x[2]*oa*ob + x[3]*ob + x[4] + offset][x[0]] = tile[x];
            }); offset += va*vb*oa*ob;
            
            HelperD::forall(sigvec_blks_["sigmar_d_bbbb"], [sigmar, ob, vb, offset, this](auto &tile, auto &x) {
                if (x[1] < x[2] && x[3] < x[4]) {
                    size_t ab_off = this->sqr_2_tri_idx(x[1], x[2], vb);
                    size_t ij_off = this->sqr_2_tri_idx(x[3], x[4], ob);
                    size_t sigmar_idx = ab_off*ob*(ob-1)/2 + ij_off + offset;
                    sigmar[sigmar_idx][x[0]] = tile[x];
                }
            }); offset += vb*ob*(ob-1)*(vb-1)/4;
        }

        /// unpack photonic part of right sigma vectors
        // ground state + hw
        if (include_u0_){
            HelperD::forall(sigvec_blks_["sigmar_u0"], [sigmar, offset](auto &tile, auto &x) {
                sigmar[offset][x[0]] = tile[x];
            }); offset += 1;
        }

        // singles + hw
        if (include_u1_){
            HelperD::forall(sigvec_blks_["sigmar_u1_aa"], [sigmar, oa, offset](auto &tile, auto &x) {
                sigmar[x[1] * oa + x[2] + offset][x[0]] = tile[x];
            }); offset += oa*va;

            HelperD::forall(sigvec_blks_["sigmar_u1_bb"], [sigmar, ob, offset](auto &tile, auto &x) {
                sigmar[x[1] * ob + x[2] + offset][x[0]] = tile[x];
            }); offset += ob*vb;
            
        }

        if (include_u2_){
            HelperD::forall(sigvec_blks_["sigmar_u2_aaaa"], [sigmar, oa, va, offset, this](auto &tile, auto &x) {
                if (x[1] < x[2] && x[3] < x[4]) {
                    size_t ab_off = this->sqr_2_tri_idx(x[1], x[2], va);
                    size_t ij_off = this->sqr_2_tri_idx(x[3], x[4], oa);
                    size_t sigmar_idx = ab_off*oa*(oa-1)/2 + ij_off + offset;
                    sigmar[sigmar_idx][x[0]] = tile[x];
                }
            }); offset += va*oa*(oa-1)*(va-1)/4;
            
            HelperD::forall(sigvec_blks_["sigmar_u2_abab"], [sigmar, oa, ob, vb, offset](auto &tile, auto &x) {
                sigmar[x[1]*vb*oa*ob + x[2]*oa*ob + x[3]*ob + x[4] + offset][x[0]] = tile[x];
            }); offset += va*vb*oa*ob;
            
            HelperD::forall(sigvec_blks_["sigmar_u2_bbbb"], [sigmar, ob, vb, offset, this](auto &tile, auto &x) {
                if (x[1] < x[2] && x[3] < x[4]) {
                    size_t ab_off = this->sqr_2_tri_idx(x[1], x[2], vb);
                    size_t ij_off = this->sqr_2_tri_idx(x[3], x[4], ob);
                    size_t sigmar_idx = ab_off*ob*(ob-1)/2 + ij_off + offset;
                    sigmar[sigmar_idx][x[0]] = tile[x];
                }
            }); offset += vb*ob*(ob-1)*(vb-1)/4;
        }

        /// unpack electronic part of left sigma vectors
        offset = 0;
        // ground state
        {
            HelperD::forall(sigvec_blks_["sigmal_0"], [sigmal](auto &tile, auto &x) {
                sigmal[0][x[0]] = tile[x];
            }); offset += 1;
        }
        // singles
        {
            HelperD::forall(sigvec_blks_["sigmal_s_aa"], [sigmal, oa, offset](auto &tile, auto &x) {
                sigmal[x[1]*oa + x[2] + offset][x[0]] = tile[x];
            }); offset += oa*va;

            HelperD::forall(sigvec_blks_["sigmal_s_bb"], [sigmal, ob, offset](auto &tile, auto &x) {
                sigmal[x[1]*ob + x[2] + offset][x[0]] = tile[x];
            }); offset += ob*vb;
        }

        // doubles
        {
            HelperD::forall(sigvec_blks_["sigmal_d_aaaa"], [sigmal, oa, va, offset, this](auto &tile, auto &x) {
                if (x[1] < x[2] && x[3] < x[4]) {
                    size_t ab_off = this->sqr_2_tri_idx(x[1], x[2], va);
                    size_t ij_off = this->sqr_2_tri_idx(x[3], x[4], oa);
                    size_t sigmal_idx = ab_off*oa*(oa-1)/2 + ij_off + offset;
                    sigmal[sigmal_idx][x[0]] = tile[x];
                }
            }); offset += va*oa*(oa-1)*(va-1)/4;

            HelperD::forall(sigvec_blks_["sigmal_d_abab"], [sigmal, oa, ob, vb, offset](auto &tile, auto &x) {
                sigmal[x[1] * vb * oa * ob + x[2] * oa * ob + x[3] * ob + x[4] +
                       offset][x[0]] = tile[x];
            }); offset += va*vb*oa*ob;

            HelperD::forall(sigvec_blks_["sigmal_d_bbbb"], [sigmal, ob, vb, offset, this](auto &tile, auto &x) {
                if (x[1] < x[2] && x[3] < x[4]) {
                    size_t ab_off = this->sqr_2_tri_idx(x[1], x[2], vb);
                    size_t ij_off = this->sqr_2_tri_idx(x[3], x[4], ob);
                    size_t sigmal_idx = ab_off*ob*(ob-1)/2 + ij_off + offset;
                    sigmal[sigmal_idx][x[0]] = tile[x];
                }
            }); offset += vb*ob*(ob-1)*(vb-1)/4;
        }

        /// unpack photonic part of left sigma vectors
        // ground state + hw
        if (include_u0_){
            HelperD::forall(sigvec_blks_["sigmal_u0"], [sigmal, offset](auto &tile, auto &x) {
                sigmal[offset][x[0]] = tile[x];
            }); offset += 1;
        }

        // singles + hw
        if (include_u1_){
            HelperD::forall(sigvec_blks_["sigmal_u1_aa"], [sigmal, oa, offset](auto &tile, auto &x) {
                sigmal[x[1] * oa + x[2] + offset][x[0]] = tile[x];
            }); offset += oa*va;

            HelperD::forall(sigvec_blks_["sigmal_u1_bb"], [sigmal, ob, offset](auto &tile, auto &x) {
                sigmal[x[1] * ob + x[2] + offset][x[0]] = tile[x];
            }); offset += ob*vb;

        }

        if (include_u2_){
            HelperD::forall(sigvec_blks_["sigmal_u2_aaaa"], [sigmal, oa, va, offset, this](auto &tile, auto &x) {
                if (x[1] < x[2] && x[3] < x[4]) {
                    size_t ab_off = this->sqr_2_tri_idx(x[1], x[2], va);
                    size_t ij_off = this->sqr_2_tri_idx(x[3], x[4], oa);
                    size_t sigmal_idx = ab_off*oa*(oa-1)/2 + ij_off + offset;
                    sigmal[sigmal_idx][x[0]] = tile[x];
                }
            }); offset += va*oa*(oa-1)*(va-1)/4;

            HelperD::forall(sigvec_blks_["sigmal_u2_abab"], [sigmal, oa, ob, vb, offset](auto &tile, auto &x) {
                sigmal[x[1]*vb*oa*ob + x[2]*oa*ob + x[3]*ob + x[4] + offset][x[0]] = tile[x];
            }); offset += va*vb*oa*ob;

            HelperD::forall(sigvec_blks_["sigmal_u2_bbbb"], [sigmal, ob, vb, offset, this](auto &tile, auto &x) {
                if (x[1] < x[2] && x[3] < x[4]) {
                    size_t ab_off = this->sqr_2_tri_idx(x[1], x[2], vb);
                    size_t ij_off = this->sqr_2_tri_idx(x[3], x[4], ob);
                    size_t sigmal_idx = ab_off*ob*(ob-1)/2 + ij_off + offset;
                    sigmal[sigmal_idx][x[0]] = tile[x];
                }
            }); offset += vb*ob*(ob-1)*(vb-1)/4;
        }

    }

    void EOM_EE_Driver::unpack_eigenvectors() {

        evec_blks_.clear(); // clear the map of trial vectors

        size_t L = M_;
        double **rerp = revec_->pointer(), **relp = levec_->pointer();
        
        // initialize operators
        evec_blks_["r0"] = HelperD::makeTensor(world_, {L}, false); // reference
        evec_blks_["r1_aa"] = HelperD::makeTensor(world_, {L, va_, oa_}, false); // alpha singles excitations
        evec_blks_["r1_bb"] = HelperD::makeTensor(world_, {L, vb_, ob_}, false); // beta singles excitations
        evec_blks_["r2_aaaa"] = HelperD::makeTensor(world_, {L, va_, va_, oa_, oa_}, false); // alpha doubles excitations
        evec_blks_["r2_abab"] = HelperD::makeTensor(world_, {L, va_, vb_, oa_, ob_}, false); // alpha-beta doubles excitations
        evec_blks_["r2_bbbb"] = HelperD::makeTensor(world_, {L, vb_, vb_, ob_, ob_}, false); // beta doubles excitations

        if (include_u0_){
            evec_blks_["s0"] = HelperD::makeTensor(world_, {L}, false); // reference + hw
        }
        if (include_u1_){
            evec_blks_["s1_aa"] = HelperD::makeTensor(world_, {L, va_, oa_}, false); // alpha singles excitations + hw
            evec_blks_["s1_bb"] = HelperD::makeTensor(world_, {L, vb_, ob_}, false); // beta singles excitations + hw
        }
        if (include_u2_){
            evec_blks_["s2_aaaa"] = HelperD::makeTensor(world_, {L, va_, va_, oa_, oa_}, false); // alpha doubles excitations + hw
            evec_blks_["s2_abab"] = HelperD::makeTensor(world_, {L, va_, vb_, oa_, ob_}, false); // alpha-beta doubles excitations + hw
            evec_blks_["s2_bbbb"] = HelperD::makeTensor(world_, {L, vb_, vb_, ob_, ob_}, false); // beta doubles excitations + hw
        }

        size_t oa = oa_, ob = ob_, // number of occupied orbitals
        va = va_, vb = vb_; // number of virtual beta orbitals

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
            }); offset += oa*va;

            evec_blks_["r1_bb"].init_elements([rerp, ob, offset](auto &I) {
                return rerp[I[0]][I[1]*ob + I[2] + offset];
            }); offset += vb*ob;
        }

        // pack redundant doubles into tensors
        {
            evec_blks_["r2_aaaa"].init_elements([rerp, va, oa, offset, this](auto &I) {
                size_t a = I[1], b = I[2], i = I[3], j = I[4];
                if (a == b || i == j) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b) ^ (i > j);
                if (a > b) std::swap(a, b);
                if (i > j) std::swap(i, j);

                size_t ab_off = this->sqr_2_tri_idx(a, b, va);
                size_t ij_off = this->sqr_2_tri_idx(i, j, oa);
                size_t upper_tri_index = ab_off*oa*(oa-1)/2 + ij_off;
                double value = rerp[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;

            }); offset += va*oa*(va-1)/2*(oa-1)/2;

            evec_blks_["r2_abab"].init_elements([rerp, oa, ob, vb, offset](auto &I) {
                return rerp[I[0]][I[1]*vb*oa*ob + I[2]*oa*ob + I[3]*ob + I[4] + offset];
            }); offset += va*vb*oa*ob;
            
            evec_blks_["r2_bbbb"].init_elements([rerp, vb, ob, offset, this](auto &I) {
                size_t a = I[1], b = I[2], i = I[3], j = I[4];
                if (a == b || i == j) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b) ^ (i > j);
                if (a > b) std::swap(a, b);
                if (i > j) std::swap(i, j);

                size_t ab_off = this->sqr_2_tri_idx(a, b, vb);
                size_t ij_off = this->sqr_2_tri_idx(i, j, ob);
                size_t upper_tri_index = ab_off*ob*(ob-1)/2 + ij_off;
                double value = rerp[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;

            }); offset += vb*ob*(vb-1)/2*(ob-1)/2;
        }

        /// pack photonic eigenvectors into TA::TArrayD
        if (include_u0_){
            // pack ground state + hw operator
            evec_blks_["s0"].init_elements([rerp, offset](auto &I) {
                return rerp[I[0]][offset];
            }); offset++;
        }
        if (include_u1_){
            // pack singles + hw operator
            evec_blks_["s1_aa"].init_elements([rerp, oa, offset](auto &I) {
                return rerp[I[0]][I[1]*oa + I[2] + offset];
            }); offset += oa*va;

            evec_blks_["s1_bb"].init_elements([rerp, ob, offset](auto &I) {
                return rerp[I[0]][I[1]*ob + I[2] + offset];
            }); offset += ob*vb;
        }
        if (include_u2_){
            // pack redundant doubles + hw operator
            evec_blks_["s2_aaaa"].init_elements([rerp, va, oa, offset, this](auto &I) {
                size_t a = I[1], b = I[2], i = I[3], j = I[4];
                if (a == b || i == j) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b) ^ (i > j);
                if (a > b) std::swap(a, b);
                if (i > j) std::swap(i, j);

                size_t ab_off = this->sqr_2_tri_idx(a, b, va);
                size_t ij_off = this->sqr_2_tri_idx(i, j, oa);
                size_t upper_tri_index = ab_off*oa*(oa-1)/2 + ij_off;
                double value = rerp[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;

            }); offset += va*oa*(va-1)/2*(oa-1)/2;

            evec_blks_["s2_abab"].init_elements([rerp, oa, ob, vb, offset](auto &I) {
                return rerp[I[0]][I[1]*vb*oa*ob + I[2]*oa*ob + I[3]*ob + I[4] + offset];
            }); offset += va*vb*oa*ob;

            evec_blks_["s2_bbbb"].init_elements([rerp, vb, ob, offset, this](auto &I) {
                size_t a = I[1], b = I[2], i = I[3], j = I[4];
                if (a == b || i == j) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b) ^ (i > j);
                if (a > b) std::swap(a, b);
                if (i > j) std::swap(i, j);

                size_t ab_off = this->sqr_2_tri_idx(a, b, vb);
                size_t ij_off = this->sqr_2_tri_idx(i, j, ob);
                size_t upper_tri_index = ab_off*ob*(ob-1)/2 + ij_off;
                double value = rerp[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;

            }); offset += vb*ob*(vb-1)/2*(ob-1)/2;
        }

        /// initialize left operators
        evec_blks_["l0"] = HelperD::makeTensor(world_, {L}, false); // reference
        evec_blks_["l1_aa"] = HelperD::makeTensor(world_, {L, oa_, va_}, false); // alpha singles excitations
        evec_blks_["l1_bb"] = HelperD::makeTensor(world_, {L, ob_, vb_}, false); // beta singles excitations
        evec_blks_["l2_aaaa"] = HelperD::makeTensor(world_, {L, oa_, oa_, va_, va_}, false); // alpha doubles excitations
        evec_blks_["l2_abab"] = HelperD::makeTensor(world_, {L, oa_, ob_, va_, vb_}, false); // alpha-beta doubles excitations
        evec_blks_["l2_bbbb"] = HelperD::makeTensor(world_, {L, ob_, ob_, vb_, vb_}, false); // beta doubles excitations

        if (include_u0_){
            evec_blks_["m0"] = HelperD::makeTensor(world_, {L}, false); // reference + hw
        }
        if (include_u1_){
            evec_blks_["m1_aa"] = HelperD::makeTensor(world_, {L, oa_, va_}, false); // alpha singles excitations + hw
            evec_blks_["m1_bb"] = HelperD::makeTensor(world_, {L, ob_, vb_}, false); // beta singles excitations + hw
        }
        if (include_u2_){
            evec_blks_["m2_aaaa"] = HelperD::makeTensor(world_, {L, oa_, oa_, va_, va_}, false); // alpha doubles excitations + hw
            evec_blks_["m2_abab"] = HelperD::makeTensor(world_, {L, oa_, ob_, va_, vb_}, false); // alpha-beta doubles excitations + hw
            evec_blks_["m2_bbbb"] = HelperD::makeTensor(world_, {L, ob_, ob_, vb_, vb_}, false); // beta doubles excitations + hw
        }
        
        /// pack electronic trial vectors into TA::TArrayD
        // pack ground state operator
        evec_blks_["l0"].init_elements([relp](auto &I) {
            return relp[I[0]][0];
        }); offset = 1;

        // pack singles into tensors
        {
            evec_blks_["l1_aa"].init_elements([relp, oa, offset](auto &I) {
                return relp[I[0]][I[2]*oa + I[1] + offset];
            }); offset += oa*va;

            evec_blks_["l1_bb"].init_elements([relp, ob, offset](auto &I) {
                return relp[I[0]][I[2]*ob + I[1] + offset];
            }); offset += vb*ob;
        }

        // pack redundant doubles into tensors
        {
            evec_blks_["l2_aaaa"].init_elements([relp, va, oa, offset, this](auto &I) {
                size_t a = I[3], b = I[4], i = I[1], j = I[2];
                if (a == b || i == j) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b) ^ (i > j);
                if (a > b) std::swap(a, b);
                if (i > j) std::swap(i, j);

                size_t ab_off = this->sqr_2_tri_idx(a, b, va);
                size_t ij_off = this->sqr_2_tri_idx(i, j, oa);
                size_t upper_tri_index = ab_off*oa*(oa-1)/2 + ij_off;
                double value = relp[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;

            }); offset += va*oa*(va-1)/2*(oa-1)/2;

            evec_blks_["l2_abab"].init_elements([relp, oa, ob, vb, offset](auto &I) {
                return relp[I[0]][I[3]*vb*oa*ob + I[4]*oa*ob + I[1]*ob + I[2] + offset];
            }); offset += va*vb*oa*ob;

            evec_blks_["l2_bbbb"].init_elements([relp, vb, ob, offset, this](auto &I) {
                size_t a = I[3], b = I[4], i = I[1], j = I[2];
                if (a == b || i == j) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b) ^ (i > j);
                if (a > b) std::swap(a, b);
                if (i > j) std::swap(i, j);

                size_t ab_off = this->sqr_2_tri_idx(a, b, vb);
                size_t ij_off = this->sqr_2_tri_idx(i, j, ob);
                size_t upper_tri_index = ab_off*ob*(ob-1)/2 + ij_off;
                double value = relp[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;

            }); offset += vb*ob*(vb-1)/2*(ob-1)/2;
        }

        /// pack photonic trial vectors into TA::TArrayD
        if (include_u0_){
            // pack ground state + hw operator
            evec_blks_["m0"].init_elements([relp, offset](auto &I) {
                return relp[I[0]][offset];
            }); offset++;
        }
        if (include_u1_){
            // pack singles + hw operator
            evec_blks_["m1_aa"].init_elements([relp, oa, offset](auto &I) {
                return relp[I[0]][I[2]*oa + I[1] + offset];
            }); offset += oa*va;

            evec_blks_["m1_bb"].init_elements([relp, ob, offset](auto &I) {
                return relp[I[0]][I[2]*ob + I[1] + offset];
            }); offset += ob*vb;
        }
        if (include_u2_){
            // pack redundant doubles + hw operator
            evec_blks_["m2_aaaa"].init_elements([relp, va, oa, offset, this](auto &I) {
                size_t a = I[3], b = I[4], i = I[1], j = I[2];
                if (a == b || i == j) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b) ^ (i > j);
                if (a > b) std::swap(a, b);
                if (i > j) std::swap(i, j);

                size_t ab_off = this->sqr_2_tri_idx(a, b, va);
                size_t ij_off = this->sqr_2_tri_idx(i, j, oa);
                size_t upper_tri_index = ab_off*oa*(oa-1)/2 + ij_off;
                double value = relp[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;

            }); offset += va*oa*(va-1)/2*(oa-1)/2;

            evec_blks_["m2_abab"].init_elements([relp, oa, ob, vb, offset](auto &I) {
                return relp[I[0]][I[3]*vb*oa*ob + I[4]*oa*ob + I[1]*ob + I[2] + offset];
            }); offset += va*vb*oa*ob;

            evec_blks_["m2_bbbb"].init_elements([relp, vb, ob, offset, this](auto &I) {
                size_t a = I[3], b = I[4], i = I[1], j = I[2];
                if (a == b || i == j) return 0.0; // return 0 for diagonal elements

                bool change_sign = (a > b) ^ (i > j);
                if (a > b) std::swap(a, b);
                if (i > j) std::swap(i, j);

                size_t ab_off = this->sqr_2_tri_idx(a, b, vb);
                size_t ij_off = this->sqr_2_tri_idx(i, j, ob);
                size_t upper_tri_index = ab_off*ob*(ob-1)/2 + ij_off;
                double value = relp[I[0]][upper_tri_index + offset];
                return change_sign ? -value : value;

            }); offset += vb*ob*(vb-1)/2*(ob-1)/2;
        }

        world_.gop.fence();
    }

    void EOM_EE_Driver::print_dominant_transitions() {
        /// get dominant transition for each root in each block

        // get pointers to eigenvectors
        double **rerp = revec_->pointer();
        double **relp = levec_->pointer();

        /// get dominant transitions for each state in each block
        double threshold = 1e-10;

        for (size_t I = 0; I < M_; I++) { // loop over roots

            // a map of the operator name to a priority queue of its dominant transitions with:
            //     the magnitude of the transition
            //     the spin of the transition
            //     the indicies of the transition
            map<string, priority_queue<pair<double, pair<string, vector<size_t>>>>> dominant_transitions;

            size_t off = 0;

            // ground state
            double value = rerp[I][0]*relp[I][0]; // get value
            if (fabs(value) > threshold) {
                dominant_transitions["l0*r0"].push({value, {"", {}}});
            }
            off++;


            size_t id = 0;
            // singles aa
            for (size_t a = 0; a < va_; a++) {
                for (size_t i = 0; i < oa_; i++) {
                    value = rerp[I][id + off]*relp[I][id + off]; // get value
                    if (fabs(value) > threshold) {
                        dominant_transitions["l1*r1"].push({value, {"a", {a+oa_, i}}});
                    }
                    id++;
                }
            }
            off += id;
            id = 0;

            // singles bb
            for (size_t a = 0; a < vb_; a++) {
                for (size_t i = 0; i < ob_; i++) {
                    value = rerp[I][id + off]*relp[I][id + off]; // get value
                    if (fabs(value) > threshold) {
                        dominant_transitions["l1*r1"].push({value, {"b", {a+ob_, i}}});
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
                            value = rerp[I][id + off] * relp[I][id + off]; // get value
                            if (fabs(value) > threshold) {
                                dominant_transitions["l2*r2"].push({value, {"aa", {a+oa_, b+oa_, i, j}}});
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
                            value = rerp[I][id + off] * relp[I][id + off]; // get value
                            if (fabs(value) > threshold) {
                                dominant_transitions["l2*r2"].push({value, {"ab", {a+oa_, b+ob_, i, j}}});
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
                            value = rerp[I][id + off] * relp[I][id + off]; // get value
                            if (fabs(value) > threshold) {
                                dominant_transitions["l2*r2"].push({value, {"bb", {a+ob_, b+ob_, i, j}}});
                            }
                            id++;
                        }
                    }
                }
            }
            off += id;
            id = 0;

            // ground state + hw
            if (include_u0_) {
                value = rerp[I][id + off] * relp[I][id + off]; // get value
                if (fabs(value) > threshold) {
                    dominant_transitions["m0*s0"].push({value, {"",{}}});
                }
                off++;
            }

            if (include_u1_) {
                // singles aa + hw
                for (size_t a = 0; a < va_; a++) {
                    for (size_t i = 0; i < oa_; i++) {
                        value = rerp[I][id + off]*relp[I][id + off]; // get value
                        if (fabs(value) > threshold) {
                            dominant_transitions["m1*s1"].push({value, {"a", {a+oa_, i}}});
                        }
                        id++;
                    }
                }
                off += id;
                id = 0;

                // singles bb + hw
                for (size_t a = 0; a < vb_; a++) {
                    for (size_t i = 0; i < ob_; i++) {
                        value = rerp[I][id + off]*relp[I][id + off]; // get value
                        if (fabs(value) > threshold) {
                            dominant_transitions["m1*s1"].push({value, {"b", {a+ob_, i}}});
                        }
                        id++;
                    }
                }
                off += id;
                id = 0;
            }


            if (include_u2_) {

                // doubles aaaa + hw
                for (size_t a = 0; a < va_; a++) {
                    for (size_t b = a + 1; b < va_; b++) {
                        for (size_t i = 0; i < oa_; i++) {
                            for (size_t j = i + 1; j < oa_; j++) {
                                value = rerp[I][id + off] * relp[I][id + off]; // get value
                                if (fabs(value) > threshold) {
                                    dominant_transitions["m2*s2"].push({value, {"aa", {a+oa_, b+oa_, i, j}}});
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
                                value = rerp[I][id + off] * relp[I][id + off]; // get value
                                if (fabs(value) > threshold) {
                                    dominant_transitions["m2*s2"].push({value, {"ab", {a+oa_, b+ob_, i, j}}});
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
                                value = rerp[I][id + off] * relp[I][id + off]; // get value
                                if (fabs(value) > threshold) {
                                    dominant_transitions["m2*s2"].push({value, {"bb", {a+ob_, b+ob_, i, j}}});
                                }
                                id++;
                            }
                        }
                    }
                }
            }
            off += id;
            id = 0;

            /// Print the dominant transitions
            Printf("\n\n  --> Root %zu <--\n", I);

            for(auto blk : dominant_transitions) { // print top 5 transitions for each block
                // grab first 5 values from priority queue
                vector<pair<double, pair<string,vector<size_t>>>> top5;
                for (size_t i = 0; i < 5; i++) {
                    if (blk.second.size() > 0) {
                        top5.push_back(blk.second.top());
                        blk.second.pop();
                    }
                }

                // print top 5 in following format:
                //     --> l1*r1_bbbb <--
                //          1: (  i  )->(  a  ), (+-)0.1234
                //          2: (  i  )->(  a  ), (+-)0.1234
                //          3: (  i  )->(  a  ), (+-)0.1234
                //          4: (  i  )->(  a  ), (+-)0.1234
                //          5: (  i  )->(  a  ), (+-)0.1234
                //     --> l2*r2_aaaa <--
                //          1: (  i,  j  )->(  a,  b  ), (+-)0.1234
                //          2: (  i,  j  )->(  a,  b  ), (+-)0.1234
                //          3: (  i,  j  )->(  a,  b  ), (+-)0.1234
                //          4: (  i,  j  )->(  a,  b  ), (+-)0.1234
                //          5: (  i,  j  )->(  a,  b  ), (+-)0.1234

                size_t top5_size = top5.size();
                if (top5_size == 0) continue;

                Printf("\n    %s", blk.first.c_str());
                for (size_t i = 0; i < top5_size; i++) {
                    pair<double, pair<string,vector<size_t>>> &tpair = top5[i];
                    string &spin = tpair.second.first; // spin
                    vector<size_t> &labels = tpair.second.second; // transition pairs

                    double val = tpair.first; // value
                    size_t size = labels.size(); // size of transition pairs
                    if (size == 0) { // ground state
                        Printf(": %15.12lf", i+1, val);
                    } else if (size == 2) { // singles
                        Printf("\n        %d: (%2d%c ) -> (%2d%c ), %15.12lf", i+1,
                               labels[1], spin[0], labels[0], spin[0], val);
                    } else if (size == 4) { // doubles
                        Printf("\n        %d: (%2d%c,%2d%c ) -> (%2d%c,%2d%c ), %15.12lf", i+1,
                               labels[2], spin[0], labels[3], spin[1],
                               labels[0], spin[0], labels[1], spin[1], val);
                    }
                }
            }
            // free memory
            dominant_transitions.clear();
        }
        Printf("\n\n");
    }


} // cc_cavity
