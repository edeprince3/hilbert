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

#include <psi4/psi4-dec.h>
#include <psi4/physconst.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsio/psio.hpp>

#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/mintshelper.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/basisset.h>
#include <psi4/lib3index/dftensor.h>
#include <psi4/libqt/qt.h>

// jk object
#include <psi4/libfock/jk.h>

// for dft
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libscf_solver/hf.h"

#include <psi4/psifiles.h>
#include <psi4/libtrans/integraltransform.h>

#include "utddft.h"

#include <misc/blas.h>
#include <misc/hilbert_psifiles.h>
#include <misc/omp.h>
#include <misc/threeindexintegrals.h>
#include <misc/nonsym_davidson_solver.h>

// diis solver
#include <misc/diis.h>

#include <algorithm>

using namespace psi;
using namespace fnocc;

namespace hilbert{ 

PolaritonicUTDDFT::PolaritonicUTDDFT(std::shared_ptr<Wavefunction> reference_wavefunction, Options& options_, std::shared_ptr<Wavefunction> dummy_wfn):
    PolaritonicHF(reference_wavefunction,options_) {
    common_init(dummy_wfn);
}

PolaritonicUTDDFT::~PolaritonicUTDDFT() {
}

void PolaritonicUTDDFT::common_init(std::shared_ptr<Wavefunction> dummy_wfn) {

    outfile->Printf("\n");
    outfile->Printf("\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    Polaritonic Unrestricted TDDFT                   *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf("\n");

    // UTDDFT only works with TDA for now
    if ( !options_.get_bool("TDSCF_TDA") ) {
        //throw PsiException("polaritonic rtddft only works with TDA df for now",__FILE__,__LINE__);
    }

    // check SCF type
    if ( options_.get_str("SCF_TYPE") != "DF" && options_.get_str("SCF_TYPE") != "CD" && options_.get_str("SCF_TYPE") != "PK") {
        throw PsiException("invalid SCF_TYPE for qed-utddft",__FILE__,__LINE__);
    }

    // ensure running in c1 symmetry
    if ( reference_wavefunction_->nirrep() > 1 ) {
        throw PsiException("qed-utddft only works with c1 symmetry for now.",__FILE__,__LINE__);
    }

    // get primary basis:
    std::shared_ptr<BasisSet> primary = reference_wavefunction_->get_basisset("ORBITAL");

    // determine the DFT functional and initialize the potential object
    psi::scf::HF* scfwfn = (psi::scf::HF*)dummy_wfn.get();
    std::shared_ptr<SuperFunctional> functional = scfwfn->functional();
    potential_ = (std::shared_ptr<VBase>)VBase::build_V(primary,functional,options_,(options_.get_str("REFERENCE") == "RKS" ? "RV" : "UV"));

    // initialize potential object
    potential_->initialize();

    // apparently compute_Vx wants me to set the density
    potential_->set_D({Da_,Db_});

    if ( options_.get_str("SCF_TYPE") == "DF" ) {

        // get auxiliary basis:
        std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");

        jk_ = (std::shared_ptr<DiskDFJK>)(new DiskDFJK(primary,auxiliary,options_));

    }else if ( options_.get_str("SCF_TYPE") == "CD" ) {

        jk_ = (std::shared_ptr<CDJK>)(new CDJK(primary,options_,options_.get_double("CHOLESKY_TOLERANCE")));

    }else if ( options_.get_str("SCF_TYPE") == "PK" ) {

        jk_ = (std::shared_ptr<PKJK>)(new PKJK(primary,options_));
    }

    // memory for jk (say, 80% of what is available)
    jk_->set_memory(0.8 * memory_);

    // integral cutoff
    jk_->set_cutoff(options_.get_double("INTS_TOLERANCE"));

    is_hf_ = functional->name() == "HF" ? true : false;

    is_x_lrc_    = functional->is_x_lrc();
    is_x_hybrid_ = functional->is_x_hybrid();
    x_omega_     = functional->x_omega();
    if ( options_["DFT_OMEGA"].has_changed() ) {
        x_omega_ = options_.get_double("DFT_OMEGA");
    }
    x_alpha_     = functional->x_alpha();
    needs_xc_    = functional->needs_xc();

    jk_->set_do_J(true);
    jk_->set_do_K(is_x_hybrid_);
    jk_->set_do_wK(is_x_lrc_);
    jk_->set_omega(x_omega_);

    jk_->initialize();

    // alpha + beta MO transformation matrix
    C_ = (std::shared_ptr<Matrix>)(new Matrix(2*nso_,2*nmo_));
    C_->zero();
    double ** cp = C_->pointer();
    double ** ca = Ca_->pointer();
    double ** cb = Cb_->pointer();

    for (size_t mu = 0; mu < nso_; mu++) {
        size_t count = 0;
        for (size_t i = 0; i < nalpha_; i++) {
            cp[mu][count++] = ca[mu][i];
        }
        for (size_t i = 0; i < nbeta_; i++) {
            cp[mu+nso_][count++] = cb[mu][i];
        }
        for (size_t i = nalpha_; i < nmo_; i++) {
            cp[mu][count++] = ca[mu][i];
        }
        for (size_t i = nbeta_; i < nmo_; i++) {
            cp[mu+nso_][count++] = cb[mu][i];
        }
    }

    // dipole integrals in spin-orbital basis
    Dipole_x_ = (std::shared_ptr<Matrix>)(new Matrix(2*nso_,2*nso_));
    Dipole_y_ = (std::shared_ptr<Matrix>)(new Matrix(2*nso_,2*nso_));
    Dipole_z_ = (std::shared_ptr<Matrix>)(new Matrix(2*nso_,2*nso_));
    Dipole_x_->zero();
    Dipole_y_->zero();
    Dipole_z_->zero();
    double ** dx = Dipole_x_->pointer();
    double ** dy = Dipole_y_->pointer();
    double ** dz = Dipole_z_->pointer();
    double ** dipole_x_p = dipole_[0]->pointer();
    double ** dipole_y_p = dipole_[1]->pointer();
    double ** dipole_z_p = dipole_[2]->pointer();
    for (size_t mu = 0; mu < nso_; mu++) {
        for (size_t nu = 0; nu < nso_; nu++) {

            dx[mu][nu]           = dipole_x_p[mu][nu];
            dx[mu+nso_][nu+nso_] = dipole_x_p[mu][nu];

            dy[mu][nu]           = dipole_y_p[mu][nu];
            dy[mu+nso_][nu+nso_] = dipole_y_p[mu][nu];

            dz[mu][nu]           = dipole_z_p[mu][nu];
            dz[mu+nso_][nu+nso_] = dipole_z_p[mu][nu];

        }
    }
    Dipole_x_->transform(C_);
    Dipole_y_->transform(C_);
    Dipole_z_->transform(C_);
}

void PolaritonicUTDDFT::compute_static_responses() {

    outfile->Printf("\n");
    outfile->Printf("    No. basis functions:            %5i\n",nso_);
    outfile->Printf("    No. alpha electrons:            %5i\n",nalpha_);
    outfile->Printf("    No. beta electrons:             %5i\n",nbeta_);

    if ( n_photon_states_ > 1 ) {
        update_cavity_terms();
    }

    double coupling_factor_x = cavity_frequency_[0] * cavity_coupling_strength_[0];
    double coupling_factor_y = cavity_frequency_[1] * cavity_coupling_strength_[1];
    double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];

    double lambda_x = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
    double lambda_y = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

    double ** dx = Dipole_x_->pointer();
    double ** dy = Dipole_y_->pointer();
    double ** dz = Dipole_z_->pointer();

    std::shared_ptr<Matrix> HCavity_z (new Matrix(n_photon_states_,n_photon_states_));
    HCavity_z->zero();
    if ( n_photon_states_ > 1 ) {
        HCavity_z->pointer()[1][1] = cavity_frequency_[2];
    }
    if ( n_photon_states_ > 2 ) {
        throw PsiException("polaritonic response properties do not work with N_PHOTON_STATES > 2",__FILE__,__LINE__);
    }

    // dimension of the problem
    int oa = nalpha_;
    int ob = nbeta_;
    int va = nso_ - oa;
    int vb = nso_ - ob;

    int N = 2*(oa*va+ob*vb) + 2;

    double * amps1 = (double*)malloc(3*N*sizeof(double));
    memset((void*)amps1, '\0', 3*N * sizeof(double));

    double * old_amps1 = (double*)malloc(3*N*sizeof(double));
    memset((void*)old_amps1, '\0', 3*N * sizeof(double));

    double * amps1_error = (double*)malloc(3*N*sizeof(double));
    memset((void*)amps1_error, '\0', 3*N * sizeof(double));

    double * ABu = (double*)malloc(3*2*nmo_*nmo_*sizeof(double));
    memset((void*)ABu, '\0', 3*2*nmo_*nmo_ * sizeof(double));

    double * ea = epsilon_a_->pointer();
    double * eb = epsilon_b_->pointer();

    std::shared_ptr<MintsHelper> mints (new MintsHelper(reference_wavefunction_));

    // dipole integrals
    std::vector<std::shared_ptr<Matrix>> mua = mints->so_dipole();
    mua[0]->transform(Ca_);
    mua[1]->transform(Ca_);
    mua[2]->transform(Ca_);

    std::vector<std::shared_ptr<Matrix>> mub = mints->so_dipole();
    mub[0]->transform(Cb_);
    mub[1]->transform(Cb_);
    mub[2]->transform(Cb_);

    // photon parts
    int L = 1;
    double *gm = (double*)malloc(3*N*L*sizeof(double));
    double *sigma_m_r = (double*)malloc(3*L*sizeof(double));
    double *sigma_m_l = (double*)malloc(3*L*sizeof(double));

    memset((void*)gm,'\0',3*N*L*sizeof(double));
    memset((void*)sigma_m_r,'\0',3*L*sizeof(double));
    memset((void*)sigma_m_l,'\0',3*L*sizeof(double));

    double * alpha = (double*)malloc(9*sizeof(double));
    memset((void*)alpha, '\0', 9*sizeof(double));

    std::shared_ptr<DIIS> diis (new DIIS(3*N));

    outfile->Printf("\n");
    outfile->Printf("    ==> First-Order Response Equations <==\n");
    outfile->Printf("\n");

    for (int iter = 0; iter < 500; iter++) {
        double damp = 0.5;
        if (iter > 10) {
            damp = 0.0;
        }
        double err = 0.0;
        C_DCOPY(3*N, amps1, 1, old_amps1, 1);
        for (int p = 0; p < 3; p++) {

            build_Au_Bu_response(N, &amps1[p*N], &ABu[p*2*nmo_*nmo_]);
            build_gm(N, 1, &amps1[p*N + (oa*va+ob*vb)], &gm[p*N]);

            for (int i = 0; i < oa; i++) {
                for (int a = 0; a < va; a++) {
                    double precon = 1.0 / (ea[a+oa] - ea[i]);
                    int ia = i * va + a;
                    amps1[p*N+ia] = damp * amps1[p*N+ia] + (1.0 - damp) * precon * (mua[p]->pointer()[i][a+oa] - ABu[p*2*nmo_*nmo_ + i*nmo_+(a+oa)] - ABu[p*2*nmo_*nmo_ + (a+oa)*nmo_+i] - 2.0 * gm[p*N+ia]);
                }
            }
            for (int i = 0; i < ob; i++) {
                for (int a = 0; a < vb; a++) {
                    int ia = oa*va + i * vb + a;
                    double precon = 1.0 / (eb[a+ob] - eb[i]);
                    amps1[p*N+ia] = damp * amps1[p*N+ia] + (1.0 - damp) * precon * (mub[p]->pointer()[i][a+ob] - ABu[p*2*nmo_*nmo_ + nmo_*nmo_ + i*nmo_+(a+ob)] - ABu[p*2*nmo_*nmo_ + nmo_*nmo_ + (a+ob)*nmo_+i] - 2.0 * gm[p*N+ia]);
                }
            }
            C_DCOPY(N, &amps1[p*N], 1, &amps1_error[p*N], 1);
            C_DAXPY(N, -1.0, &old_amps1[p*N], 1, &amps1_error[p*N], 1);

            // photon part
            build_sigma_m(N, L, &amps1[p*N], &amps1[p*N], &amps1[p*N + (oa*va+ob*vb)], &sigma_m_r[p], &sigma_m_l[p]);
            amps1[p*N + (oa*va+ob*vb)] = -sigma_m_r[p] / cavity_frequency_[2];
            
            alpha[p*3 + p] = 0.0;
            for (int i = 0; i < oa; i++) {
                for (int a = 0; a < va; a++) {
                    int ia = i * va + a;
                    alpha[p*3 + p] += 2.0 * mua[p]->pointer()[i][a+oa] * amps1[p*N+ia];
                }
            }
            for (int i = 0; i < ob; i++) {
                for (int a = 0; a < vb; a++) {
                    int ia = oa*va + i * vb + a;
                    alpha[p*3 + p] += 2.0 * mub[p]->pointer()[i][a+ob] * amps1[p*N+ia];
                }
            }
        }

        err += C_DDOT(N, amps1_error, 1, amps1_error, 1);
        err = sqrt(err);

        // DIIS extrapolation
        diis->WriteVector(amps1);
        diis->WriteErrorVector(amps1_error);
        diis->Extrapolate(amps1);

        outfile->Printf("%5i %20.12lf %20.12lf %20.12lf %20.12lf\n", iter, alpha[0*3+0], alpha[1*3+1], alpha[2*3+2], err);
        if ( err < options_.get_double("D_CONVERGENCE") ) {
            break;
        }
    }

    for (int p = 0; p < 3; p++) {
        build_Au_Bu_response(N, &amps1[p*N], &ABu[p*2*nmo_*nmo_]);

        // photon part
        build_sigma_m(N, L, &amps1[p*N], &amps1[p*N], &amps1[p*N + (oa*va+ob*vb)], &sigma_m_r[p], &sigma_m_l[p]);
        amps1[p*N + (oa*va+ob*vb)] = - sigma_m_r[p] / cavity_frequency_[2];
    }

    outfile->Printf("\n");
    outfile->Printf("    ==> Static Polarizability <==\n");
    outfile->Printf("\n");
    std::vector<std::string> dir {"X", "Y", "Z"};
    for (int p = 0; p < 3; p++) {
        for (int q = p; q < 3; q++) {
            alpha[p*3 + q] = 0.0;
            for (int i = 0; i < oa; i++) {
                for (int a = 0; a < va; a++) {
                    int ia = i * va + a;
                    alpha[p*3 + q] += 2.0 * mua[p]->pointer()[i][a+oa] * amps1[q*N+ia];
                }
            }
            for (int i = 0; i < ob; i++) {
                for (int a = 0; a < vb; a++) {
                    int ia = oa*va + i * vb + a;
                    alpha[p*3 + q] += 2.0 * mub[p]->pointer()[i][a+ob] * amps1[q*N+ia];
                }
            }
            outfile->Printf("    ALPHA(%s%s) %20.12lf\n", dir[p].c_str(), dir[q].c_str(), alpha[p*3 + q]);

            // add polarizabilities to psi variables
            std::string label = "QED-DFT ALPHA(";
            label += dir[p] + dir[q] + ")";
            std::transform(label.begin(), label.end(), label.begin(),
                [](unsigned char c) { return std::toupper(c); });
            Process::environment.globals[label] = alpha[p*3 + q];
        }
    }

    free(alpha);
    free(gm);
    free(sigma_m_r);
    free(sigma_m_l);
    free(old_amps1);
    free(amps1_error);

    // static hyperpolarizability is only implemented for QED-HF
    if ( !is_hf_ ) {
        free(amps1);
        free(ABu);
        return;
    }

    outfile->Printf("\n");
    outfile->Printf("    ==> Static Hyperpolarizability <==\n");
    outfile->Printf("\n");

    for (int p = 0; p < 3; p++) {
        for (int q = p; q < 3; q++) {
            for (int r = q; r < 3; r++) {
                double beta = 0.0;
                // alpha-spin
                for (int i = 0; i < oa; i++) {
                    for (int j = 0; j < oa; j++) {
                        double dum_pq = 0.0;
                        double dum_pr = 0.0;
                        double dum_qr = 0.0;
                        for (int a = 0; a < va; a++) {
                            int ia = i * va + a;
                            int ja = j * va + a;
                            dum_pq += amps1[p*N + ia] * amps1[q*N + ja];
                            dum_pr += amps1[p*N + ia] * amps1[r*N + ja];
                            dum_qr += amps1[q*N + ia] * amps1[r*N + ja];
                        }
                        int ij = i * nmo_ + j;
                        int ji = j * nmo_ + i;
                        beta -= 2 * (ABu[p*2*nmo_*nmo_ + ij] + ABu[p*2*nmo_*nmo_ + ji]) * dum_qr; // symmetrize to pick up all six combinations
                        beta -= 2 * (ABu[q*2*nmo_*nmo_ + ij] + ABu[q*2*nmo_*nmo_ + ji]) * dum_pr;
                        beta -= 2 * (ABu[r*2*nmo_*nmo_ + ij] + ABu[r*2*nmo_*nmo_ + ji]) * dum_pq;
                    }
                }
                for (int a = 0; a < va; a++) {
                    for (int b = 0; b < va; b++) {
                        double dum_pq = 0.0;
                        double dum_pr = 0.0;
                        double dum_qr = 0.0;
                        for (int i = 0; i < oa; i++) {
                            int ia = i * va + a;
                            int ib = i * va + b;
                            dum_pq += amps1[p*N + ia] * amps1[q*N + ib];
                            dum_pr += amps1[p*N + ia] * amps1[r*N + ib];
                            dum_qr += amps1[q*N + ia] * amps1[r*N + ib];
                        }
                        int ab = (a+oa) * nmo_ + (b+oa);
                        int ba = (b+oa) * nmo_ + (a+oa);
                        beta += 2 * (ABu[p*2*nmo_*nmo_ + ab] + ABu[p*2*nmo_*nmo_ + ba]) * dum_qr; // symmetrize to pick up all six combinations
                        beta += 2 * (ABu[q*2*nmo_*nmo_ + ab] + ABu[q*2*nmo_*nmo_ + ba]) * dum_pr;
                        beta += 2 * (ABu[r*2*nmo_*nmo_ + ab] + ABu[r*2*nmo_*nmo_ + ba]) * dum_pq;
                    }
                }
                // beta-spin
                for (int i = 0; i < ob; i++) {
                    for (int j = 0; j < ob; j++) {
                        double dum_pq = 0.0;
                        double dum_pr = 0.0;
                        double dum_qr = 0.0;
                        for (int a = 0; a < vb; a++) {
                            int ia = i * vb + a + oa*va;
                            int ja = j * vb + a + oa*va;
                            dum_pq += amps1[p*N + ia] * amps1[q*N + ja];
                            dum_pr += amps1[p*N + ia] * amps1[r*N + ja];
                            dum_qr += amps1[q*N + ia] * amps1[r*N + ja];
                        }
                        int ij = i * nmo_ + j + nmo_*nmo_;
                        int ji = j * nmo_ + i + nmo_*nmo_;
                        beta -= 2 * (ABu[p*2*nmo_*nmo_ + ij] +ABu[p*2*nmo_*nmo_ + ji]) * dum_qr; // symmetrize to pick up all six combinations
                        beta -= 2 * (ABu[q*2*nmo_*nmo_ + ij] +ABu[q*2*nmo_*nmo_ + ji]) * dum_pr;
                        beta -= 2 * (ABu[r*2*nmo_*nmo_ + ij] +ABu[r*2*nmo_*nmo_ + ji]) * dum_pq;
                    }
                }
                for (int a = 0; a < vb; a++) {
                    for (int b = 0; b < vb; b++) {
                        double dum_pq = 0.0;
                        double dum_pr = 0.0;
                        double dum_qr = 0.0;
                        for (int i = 0; i < ob; i++) {
                            int ia = i * vb + a + oa*va;
                            int ib = i * vb + b + oa*va;
                            dum_pq += amps1[p*N + ia] * amps1[q*N + ib];
                            dum_pr += amps1[p*N + ia] * amps1[r*N + ib];
                            dum_qr += amps1[q*N + ia] * amps1[r*N + ib];
                        }
                        int ab = (a+ob) * nmo_ + (b+ob) + nmo_*nmo_;
                        int ba = (b+ob) * nmo_ + (a+ob) + nmo_*nmo_;
                        beta += 2 * (ABu[p*2*nmo_*nmo_ + ab] +ABu[p*2*nmo_*nmo_ + ba]) * dum_qr; // symmetrize to pick up all six combinations
                        beta += 2 * (ABu[q*2*nmo_*nmo_ + ab] +ABu[q*2*nmo_*nmo_ + ba]) * dum_pr;
                        beta += 2 * (ABu[r*2*nmo_*nmo_ + ab] +ABu[r*2*nmo_*nmo_ + ba]) * dum_pq;
                    }
                }

                // derivative wrt field plus second derivative wrt wfn parameters
                //beta += -2.00 * einsum('ji,ai,aj', f[o, o], t1, t1, optimize=['einsum_path', (0, 1), (0, 1)])
                //beta +=  2.00 * einsum('ab,ai,bi', f[v, v], t1, t1, optimize=['einsum_path', (0, 1), (0, 1)])

                // also, third derivative wrt wfn parameters has similar structure
                //beta += -2.00 * einsum('ji,ai,aj,', dipole[o, o], t1, t1, t0_1p, optimize=['einsum_path', (0, 1), (0, 2), (0, 1)])
                //beta +=  2.00 * einsum('ab,ai,bi,', dipole[v, v], t1, t1, t0_1p, optimize=['einsum_path', (0, 1), (0, 2), (0, 1)])

                // alpha-spin
                for (int i = 0; i < oa; i++) {
                    for (int j = 0; j < oa; j++) {
                        double dum_pq = 0.0;
                        double dum_pr = 0.0;
                        double dum_qr = 0.0;
                        for (int a = 0; a < va; a++) {
                            int ia = i * va + a;
                            int ja = j * va + a;
                            dum_pq += amps1[p*N + ja] * amps1[q*N + ia];
                            dum_pr += amps1[p*N + ja] * amps1[r*N + ia];
                            dum_qr += amps1[q*N + ja] * amps1[r*N + ia];
                        }
                        beta += 2 * dum_qr * mua[p]->pointer()[j][i];
                        beta += 2 * dum_pr * mua[q]->pointer()[j][i];
                        beta += 2 * dum_pq * mua[r]->pointer()[j][i];

                        // extra factor of 2 because we get terms from b and b^ in H
                        beta += +4 * dum_qr * amps1[p*N + (oa*va+ob*vb)] * mua[2]->pointer()[j][i] * coupling_factor_z;
                        beta += +4 * dum_pr * amps1[q*N + (oa*va+ob*vb)] * mua[2]->pointer()[j][i] * coupling_factor_z;
                        beta += +4 * dum_pq * amps1[r*N + (oa*va+ob*vb)] * mua[2]->pointer()[j][i] * coupling_factor_z;
                    }
                }

                for (int a = 0; a < va; a++) {
                    for (int b = 0; b < va; b++) {
                        double dum_pq = 0.0;
                        double dum_pr = 0.0;
                        double dum_qr = 0.0;
                        for (int i = 0; i < oa; i++) {
                            int ia = i * va + a;
                            int ib = i * va + b;
                            dum_pq += amps1[p*N + ia] * amps1[q*N + ib];
                            dum_pr += amps1[p*N + ia] * amps1[r*N + ib];
                            dum_qr += amps1[q*N + ia] * amps1[r*N + ib];
                        }
                        beta -= 2 * dum_qr * mua[p]->pointer()[a+oa][b+oa];
                        beta -= 2 * dum_pr * mua[q]->pointer()[a+oa][b+oa];
                        beta -= 2 * dum_pq * mua[r]->pointer()[a+oa][b+oa];

                        // extra factor of 2 because we get terms from b and b^ in H
                        beta -= +4 * dum_qr * amps1[p*N + (oa*va+ob*vb)] * mua[2]->pointer()[a+oa][b+oa] * coupling_factor_z;
                        beta -= +4 * dum_pr * amps1[q*N + (oa*va+ob*vb)] * mua[2]->pointer()[a+oa][b+oa] * coupling_factor_z;
                        beta -= +4 * dum_pq * amps1[r*N + (oa*va+ob*vb)] * mua[2]->pointer()[a+oa][b+oa] * coupling_factor_z;
                    }
                }

                // beta-spin
                for (int i = 0; i < ob; i++) {
                    for (int j = 0; j < ob; j++) {
                        double dum_pq = 0.0;
                        double dum_pr = 0.0;
                        double dum_qr = 0.0;
                        for (int a = 0; a < vb; a++) {
                            int ia = i * vb + a + oa*va;
                            int ja = j * vb + a + oa*va;
                            dum_pq += amps1[p*N + ja] * amps1[q*N + ia];
                            dum_pr += amps1[p*N + ja] * amps1[r*N + ia];
                            dum_qr += amps1[q*N + ja] * amps1[r*N + ia];
                        }
                        beta += 2 * dum_qr * mub[p]->pointer()[j][i];
                        beta += 2 * dum_pr * mub[q]->pointer()[j][i];
                        beta += 2 * dum_pq * mub[r]->pointer()[j][i];

                        // extra factor of 2 because we get terms from b and b^ in H
                        beta += +4 * dum_qr * amps1[p*N + (oa*va+ob*vb)] * mub[2]->pointer()[j][i] * coupling_factor_z;
                        beta += +4 * dum_pr * amps1[q*N + (oa*va+ob*vb)] * mub[2]->pointer()[j][i] * coupling_factor_z;
                        beta += +4 * dum_pq * amps1[r*N + (oa*va+ob*vb)] * mub[2]->pointer()[j][i] * coupling_factor_z;
                    }
                }
                for (int a = 0; a < vb; a++) {
                    for (int b = 0; b < vb; b++) {
                        double dum_pq = 0.0;
                        double dum_pr = 0.0;
                        double dum_qr = 0.0;
                        for (int i = 0; i < ob; i++) {
                            int ia = i * vb + a + oa*va;
                            int ib = i * vb + b + oa*va;
                            dum_pq += amps1[p*N + ia] * amps1[q*N + ib];
                            dum_pr += amps1[p*N + ia] * amps1[r*N + ib];
                            dum_qr += amps1[q*N + ia] * amps1[r*N + ib];
                        }
                        beta -= 2 * dum_qr * mub[p]->pointer()[a+ob][b+ob];
                        beta -= 2 * dum_pr * mub[q]->pointer()[a+ob][b+ob];
                        beta -= 2 * dum_pq * mub[r]->pointer()[a+ob][b+ob];

                        // extra factor of 2 because we get terms from b and b^ in H
                        beta -= +4 * dum_qr * amps1[p*N + (oa*va+ob*vb)] * mub[2]->pointer()[a+ob][b+ob] * coupling_factor_z;
                        beta -= +4 * dum_pr * amps1[q*N + (oa*va+ob*vb)] * mub[2]->pointer()[a+ob][b+ob] * coupling_factor_z;
                        beta -= +4 * dum_pq * amps1[r*N + (oa*va+ob*vb)] * mub[2]->pointer()[a+ob][b+ob] * coupling_factor_z;
                    }
                }

                // third derivative of dipole-self energy 
              
                if ( options_.get_bool("QED_USE_RELAXED_ORBITALS") ) {

                    double dipole_Ja_p = 0.0;
                    double dipole_Ja_q = 0.0;
                    double dipole_Ja_r = 0.0;
                    for (int j = 0; j < oa; j++) {
                        for (int a = 0; a < va; a++) {
                            double dja = mua[2]->pointer()[j][a+oa];
                            int aj = j * va + a;
                            dipole_Ja_p += 2.0 * dja * lambda_z * lambda_z * amps1[p*N + aj];
                            dipole_Ja_q += 2.0 * dja * lambda_z * lambda_z * amps1[q*N + aj];
                            dipole_Ja_r += 2.0 * dja * lambda_z * lambda_z * amps1[r*N + aj];
                        }
                    }
                    double dipole_Jb_p = 0.0;
                    double dipole_Jb_q = 0.0;
                    double dipole_Jb_r = 0.0;
                    for (int j = 0; j < ob; j++) {
                        for (int a = 0; a < vb; a++) {
                            double dja = mub[2]->pointer()[j][a+ob];
                            int aj = j * vb + a + oa*va;
                            dipole_Jb_p += 2.0 * dja * lambda_z * lambda_z * amps1[p*N + aj];
                            dipole_Jb_q += 2.0 * dja * lambda_z * lambda_z * amps1[q*N + aj];
                            dipole_Jb_r += 2.0 * dja * lambda_z * lambda_z * amps1[r*N + aj];
                        }
                    }
                    // alpha-alpha (exchange)
                    for (int i = 0; i < oa; i++) {
                        for (int k = 0; k < oa; k++) {
                            double tmp_ik_pq = 0.0;
                            double tmp_ik_qp = 0.0;
                            double tmp_ik_qr = 0.0;
                            double tmp_ik_rq = 0.0;
                            double tmp_ik_pr = 0.0;
                            double tmp_ik_rp = 0.0;
                            for (int b = 0; b < va; b++) {
                                int bi = i * va + b;
                                int bk = k * va + b;
                                tmp_ik_qr += 2 * amps1[q*N + bi] * amps1[r*N + bk] * lambda_z * lambda_z;
                                tmp_ik_rq += 2 * amps1[r*N + bi] * amps1[q*N + bk] * lambda_z * lambda_z;
                                tmp_ik_pq += 2 * amps1[p*N + bi] * amps1[q*N + bk] * lambda_z * lambda_z;
                                tmp_ik_qp += 2 * amps1[q*N + bi] * amps1[p*N + bk] * lambda_z * lambda_z;
                                tmp_ik_pr += 2 * amps1[p*N + bi] * amps1[r*N + bk] * lambda_z * lambda_z;
                                tmp_ik_rp += 2 * amps1[r*N + bi] * amps1[p*N + bk] * lambda_z * lambda_z;
                            }
                            for (int j = 0; j < oa; j++) {
                                double dji = mua[2]->pointer()[j][i];
                                for (int a = 0; a < va; a++) {
                                    double dip = dji * mua[2]->pointer()[k][a+oa];
                                    int aj = j * va + a;
                                    beta += dip * amps1[p*N + aj] * (tmp_ik_qr + tmp_ik_rq); 
                                    beta += dip * amps1[q*N + aj] * (tmp_ik_pr + tmp_ik_rp); 
                                    beta += dip * amps1[r*N + aj] * (tmp_ik_pq + tmp_ik_qp);
                                }
                            }
                        }
                    }
                    // beta-beta (exchange)
                    for (int i = 0; i < ob; i++) {
                        for (int k = 0; k < ob; k++) {
                            double tmp_ik_pq = 0.0;
                            double tmp_ik_qp = 0.0;
                            double tmp_ik_qr = 0.0;
                            double tmp_ik_rq = 0.0;
                            double tmp_ik_pr = 0.0;
                            double tmp_ik_rp = 0.0;
                            for (int b = 0; b < vb; b++) {
                                int bi = i * vb + b + oa*va;
                                int bk = k * vb + b + oa*va;
                                tmp_ik_qr += 2 * amps1[q*N + bi] * amps1[r*N + bk] * lambda_z * lambda_z;
                                tmp_ik_rq += 2 * amps1[r*N + bi] * amps1[q*N + bk] * lambda_z * lambda_z;
                                tmp_ik_pq += 2 * amps1[p*N + bi] * amps1[q*N + bk] * lambda_z * lambda_z;
                                tmp_ik_qp += 2 * amps1[q*N + bi] * amps1[p*N + bk] * lambda_z * lambda_z;
                                tmp_ik_pr += 2 * amps1[p*N + bi] * amps1[r*N + bk] * lambda_z * lambda_z;
                                tmp_ik_rp += 2 * amps1[r*N + bi] * amps1[p*N + bk] * lambda_z * lambda_z;
                            }
                            for (int j = 0; j < ob; j++) {
                                double dji = mub[2]->pointer()[j][i];
                                for (int a = 0; a < vb; a++) {
                                    double dip = dji * mub[2]->pointer()[k][a+ob];
                                    int aj = j * vb + a + oa*va;
                                    beta += dip * amps1[p*N + aj] * (tmp_ik_qr + tmp_ik_rq); 
                                    beta += dip * amps1[q*N + aj] * (tmp_ik_pr + tmp_ik_rp); 
                                    beta += dip * amps1[r*N + aj] * (tmp_ik_pq + tmp_ik_qp); 
                                }
                            }
                        }
                    }
                    // alpha-alpha (coulomb) and alpha-beta (coulomb)
                    for (int k = 0; k < oa; k++) {
                        for (int i = 0; i < oa; i++) {
                            double dki = mua[2]->pointer()[k][i];
                            double tmp_pq = 0.0;
                            double tmp_qp = 0.0;
                            double tmp_qr = 0.0;
                            double tmp_rq = 0.0;
                            double tmp_pr = 0.0;
                            double tmp_rp = 0.0;
                            for (int b = 0; b < va; b++) {
                                int bi = i * va + b;
                                int bk = k * va + b;
                                tmp_pq += amps1[p*N + bi] * amps1[q*N + bk];
                                tmp_qp += amps1[q*N + bi] * amps1[p*N + bk];
                                tmp_pr += amps1[p*N + bi] * amps1[r*N + bk];
                                tmp_rp += amps1[r*N + bi] * amps1[p*N + bk];
                                tmp_qr += amps1[q*N + bi] * amps1[r*N + bk];
                                tmp_rq += amps1[r*N + bi] * amps1[q*N + bk];
                            }
                            beta += -(dipole_Ja_p + dipole_Jb_p) * dki * (tmp_qr + tmp_rq);
                            beta += -(dipole_Ja_q + dipole_Jb_q) * dki * (tmp_pr + tmp_rp);
                            beta += -(dipole_Ja_r + dipole_Jb_r) * dki * (tmp_pq + tmp_qp);
                        }
                    }
                    // beta-beta (coulomb) and beta-alpha (coulomb)
                    for (int k = 0; k < ob; k++) {
                        for (int i = 0; i < ob; i++) {
                            double dki = mub[2]->pointer()[k][i];
                            double tmp_pq = 0.0;
                            double tmp_qp = 0.0;
                            double tmp_qr = 0.0;
                            double tmp_rq = 0.0;
                            double tmp_pr = 0.0;
                            double tmp_rp = 0.0;
                            for (int b = 0; b < vb; b++) {
                                int bi = i * vb + b + oa*va;
                                int bk = k * vb + b + oa*va;
                                tmp_pq += amps1[p*N + bi] * amps1[q*N + bk];
                                tmp_qp += amps1[q*N + bi] * amps1[p*N + bk];
                                tmp_pr += amps1[p*N + bi] * amps1[r*N + bk];
                                tmp_rp += amps1[r*N + bi] * amps1[p*N + bk];
                                tmp_qr += amps1[q*N + bi] * amps1[r*N + bk];
                                tmp_rq += amps1[r*N + bi] * amps1[q*N + bk];
                            }
                            beta += -(dipole_Ja_p + dipole_Jb_p) * dki * (tmp_qr + tmp_rq);
                            beta += -(dipole_Ja_q + dipole_Jb_q) * dki * (tmp_pr + tmp_rp);
                            beta += -(dipole_Ja_r + dipole_Jb_r) * dki * (tmp_pq + tmp_qp);
                        }
                    }
                    dipole_Ja_p = 0.0;
                    dipole_Ja_q = 0.0;
                    dipole_Ja_r = 0.0;
                    for (int i = 0; i < oa; i++) {
                        for (int b = 0; b < va; b++) {
                            int bi = i * va + b;
                            double dib = mua[2]->pointer()[i][b+oa];
                            dipole_Ja_p += 2.0 * dib * amps1[p*N + bi] * lambda_z * lambda_z;
                            dipole_Ja_q += 2.0 * dib * amps1[q*N + bi] * lambda_z * lambda_z;
                            dipole_Ja_r += 2.0 * dib * amps1[r*N + bi] * lambda_z * lambda_z;
                        }
                    }
                    dipole_Jb_p = 0.0;
                    dipole_Jb_q = 0.0;
                    dipole_Jb_r = 0.0;
                    for (int i = 0; i < ob; i++) {
                        for (int b = 0; b < vb; b++) {
                            int bi = i * vb + b + oa*va;
                            double dib = mub[2]->pointer()[i][b+ob];
                            dipole_Jb_p += 2.0 * dib * amps1[p*N + bi] * lambda_z * lambda_z;
                            dipole_Jb_q += 2.0 * dib * amps1[q*N + bi] * lambda_z * lambda_z;
                            dipole_Jb_r += 2.0 * dib * amps1[r*N + bi] * lambda_z * lambda_z;
                        }
                    }
                    // alpha-alpha (exchange)
                    for (int a = 0; a < va; a++) {
                        for (int c = 0; c < va; c++) {
                            double tmp_ac_pq = 0.0;
                            double tmp_ac_qp = 0.0;
                            double tmp_ac_pr = 0.0;
                            double tmp_ac_rp = 0.0;
                            double tmp_ac_qr = 0.0;
                            double tmp_ac_rq = 0.0;
                            for (int j = 0; j < oa; j++) {
                                int aj = j * va + a;
                                int cj = j * va + c;
                                tmp_ac_pq += 2.0 * lambda_z * lambda_z * amps1[p*N + aj] * amps1[q*N + cj];
                                tmp_ac_qp += 2.0 * lambda_z * lambda_z * amps1[q*N + aj] * amps1[p*N + cj];
                                tmp_ac_pr += 2.0 * lambda_z * lambda_z * amps1[p*N + aj] * amps1[r*N + cj];
                                tmp_ac_rp += 2.0 * lambda_z * lambda_z * amps1[r*N + aj] * amps1[p*N + cj];
                                tmp_ac_qr += 2.0 * lambda_z * lambda_z * amps1[q*N + aj] * amps1[r*N + cj];
                                tmp_ac_rq += 2.0 * lambda_z * lambda_z * amps1[r*N + aj] * amps1[q*N + cj];
                            }
                            for (int i = 0; i < oa; i++) {
                                double dic = mua[2]->pointer()[i][c+oa];
                                for (int b = 0; b < va; b++) {
                                    double dip =  dic * mua[2]->pointer()[a+oa][b+oa];
                                    int bi = i * va + b;
                                    beta -= dip * amps1[p*N + bi] * tmp_ac_qr;
                                    beta -= dip * amps1[p*N + bi] * tmp_ac_rq;
                                    beta -= dip * amps1[q*N + bi] * tmp_ac_pr;
                                    beta -= dip * amps1[q*N + bi] * tmp_ac_rp;
                                    beta -= dip * amps1[r*N + bi] * tmp_ac_pq;
                                    beta -= dip * amps1[r*N + bi] * tmp_ac_qp;
                                }
                            }
                        }
                    }
                    // beta-beta (exchange)
                    for (int a = 0; a < vb; a++) {
                        for (int c = 0; c < vb; c++) {
                            double tmp_ac_pq = 0.0;
                            double tmp_ac_qp = 0.0;
                            double tmp_ac_pr = 0.0;
                            double tmp_ac_rp = 0.0;
                            double tmp_ac_qr = 0.0;
                            double tmp_ac_rq = 0.0;
                            for (int j = 0; j < ob; j++) {
                                int aj = j * vb + a + oa*va;
                                int cj = j * vb + c + oa*va;
                                tmp_ac_pq += 2.0 * lambda_z * lambda_z * amps1[p*N + aj] * amps1[q*N + cj];
                                tmp_ac_qp += 2.0 * lambda_z * lambda_z * amps1[q*N + aj] * amps1[p*N + cj];
                                tmp_ac_pr += 2.0 * lambda_z * lambda_z * amps1[p*N + aj] * amps1[r*N + cj];
                                tmp_ac_rp += 2.0 * lambda_z * lambda_z * amps1[r*N + aj] * amps1[p*N + cj];
                                tmp_ac_qr += 2.0 * lambda_z * lambda_z * amps1[q*N + aj] * amps1[r*N + cj];
                                tmp_ac_rq += 2.0 * lambda_z * lambda_z * amps1[r*N + aj] * amps1[q*N + cj];
                            }
                            for (int i = 0; i < ob; i++) {
                                double dic = mub[2]->pointer()[i][c+ob];
                                for (int b = 0; b < vb; b++) {
                                    double dip =  dic * mub[2]->pointer()[a+ob][b+ob];
                                    int bi = i * vb + b + oa*va;
                                    beta -= dip * amps1[p*N + bi] * tmp_ac_qr;
                                    beta -= dip * amps1[p*N + bi] * tmp_ac_rq;
                                    beta -= dip * amps1[q*N + bi] * tmp_ac_pr;
                                    beta -= dip * amps1[q*N + bi] * tmp_ac_rp;
                                    beta -= dip * amps1[r*N + bi] * tmp_ac_pq;
                                    beta -= dip * amps1[r*N + bi] * tmp_ac_qp;
                                }
                            }
                        }
                    }
                    // alpha-alpha (coulomb) + alpha-beta (coulomb)
                    for (int a = 0; a < va; a++) {
                        for (int c = 0; c < va; c++) {
                            double dac = mua[2]->pointer()[a+oa][c+oa];
                            double tmp_pq = 0.0;
                            double tmp_qp = 0.0;
                            double tmp_qr = 0.0;
                            double tmp_rq = 0.0;
                            double tmp_pr = 0.0;
                            double tmp_rp = 0.0;
                            for (int j = 0; j < oa; j++) {
                                int aj = j * va + a;
                                int cj = j * va + c;
                                tmp_pq += amps1[p*N + aj] * amps1[q*N + cj];
                                tmp_qp += amps1[q*N + aj] * amps1[p*N + cj];
                                tmp_pr += amps1[p*N + aj] * amps1[r*N + cj];
                                tmp_rp += amps1[r*N + aj] * amps1[p*N + cj];
                                tmp_qr += amps1[q*N + aj] * amps1[r*N + cj];
                                tmp_rq += amps1[r*N + aj] * amps1[q*N + cj];
                            }
                            beta -= -(dipole_Ja_p + dipole_Jb_p) * dac * (tmp_qr + tmp_rq);
                            beta -= -(dipole_Ja_q + dipole_Jb_q) * dac * (tmp_pr + tmp_rp);
                            beta -= -(dipole_Ja_r + dipole_Jb_r) * dac * (tmp_pq + tmp_qp);
                        }
                    }
                    // beta-alpha (coulomb) + beta-beta (coulomb)
                    for (int a = 0; a < vb; a++) {
                        for (int c = 0; c < vb; c++) {
                            double dac = mub[2]->pointer()[a+ob][c+ob];
                            double tmp_pq = 0.0;
                            double tmp_qp = 0.0;
                            double tmp_qr = 0.0;
                            double tmp_rq = 0.0;
                            double tmp_pr = 0.0;
                            double tmp_rp = 0.0;
                            for (int j = 0; j < ob; j++) {
                                int aj = j * vb + a + oa*va;
                                int cj = j * vb + c + oa*va;
                                tmp_pq += amps1[p*N + aj] * amps1[q*N + cj];
                                tmp_qp += amps1[q*N + aj] * amps1[p*N + cj];
                                tmp_pr += amps1[p*N + aj] * amps1[r*N + cj];
                                tmp_rp += amps1[r*N + aj] * amps1[p*N + cj];
                                tmp_qr += amps1[q*N + aj] * amps1[r*N + cj];
                                tmp_rq += amps1[r*N + aj] * amps1[q*N + cj];
                            }
                            beta -= -(dipole_Ja_p + dipole_Jb_p) * dac * (tmp_qr + tmp_rq);
                            beta -= -(dipole_Ja_q + dipole_Jb_q) * dac * (tmp_pr + tmp_rp);
                            beta -= -(dipole_Ja_r + dipole_Jb_r) * dac * (tmp_pq + tmp_qp);
                        }
                    }
                }

                outfile->Printf("    BETA(%s%s%s) %20.12lf\n", dir[p].c_str(), dir[q].c_str(), dir[r].c_str(), beta);

                // add hyperpolarizabilities to psi variables
                std::string label = "QED-DFT BETA(";
                label += dir[p] + dir[q] + dir[r] + ")";
                std::transform(label.begin(), label.end(), label.begin(),
                    [](unsigned char c) { return std::toupper(c); });
                Process::environment.globals[label.c_str()] = beta;
            }
        }
    }

    free(amps1);
    free(ABu);

    return;
}

double PolaritonicUTDDFT::compute_energy() {

    outfile->Printf("\n");
    outfile->Printf("    No. basis functions:            %5i\n",nso_);
    outfile->Printf("    No. alpha electrons:            %5i\n",nalpha_);
    outfile->Printf("    No. beta electrons:             %5i\n",nbeta_);

    if ( n_photon_states_ > 1 ) {
        update_cavity_terms();
    }

    double coupling_factor_x = cavity_frequency_[0] * cavity_coupling_strength_[0];
    double coupling_factor_y = cavity_frequency_[1] * cavity_coupling_strength_[1];
    double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];

    double lambda_x = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
    double lambda_y = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

    double ** dx = Dipole_x_->pointer();
    double ** dy = Dipole_y_->pointer();
    double ** dz = Dipole_z_->pointer();

    std::shared_ptr<Matrix> HCavity_z (new Matrix(n_photon_states_,n_photon_states_));
    HCavity_z->zero();
    if ( n_photon_states_ > 1 ) {
        HCavity_z->pointer()[1][1] = cavity_frequency_[2];
    }
    if ( n_photon_states_ > 2 ) {
        throw PsiException("polaritonic TDDFT does not work with N_PHOTON_STATES > 2",__FILE__,__LINE__);
    }

    // ok, try solving iteratively with davison
    //std::shared_ptr<Nonsym_DavidsonSolver> david (new Nonsym_DavidsonSolver());
    std::shared_ptr<Nonsym_DavidsonSolver> david (new Nonsym_DavidsonSolver());

    // dimension of the problem
    int oa = nalpha_;
    int ob = nbeta_;
    int va = nso_ - oa;
    int vb = nso_ - ob;

    int N = 2*(oa*va+ob*vb) + 2;

    // number of desired roots
    int M = options_.get_int("NUM_ROOTS");

    // maximum subspace size
    size_t MD = options_.get_int("MAXDIM");
    size_t ID = options_.get_int("INDIM");
    size_t maxdim = MD*M;
    size_t init_dim = ID*M;
    if ( M > N ) {
        outfile->Printf("    WARNING: NUM_ROOTS > N, NUM_ROOTS = %i, N = %i\n",M,N);
        outfile->Printf("        setting NUM_ROOTS = N\n");
        M = N;
    }
    if ( init_dim < M ) {
        outfile->Printf("    WARNING: INDIM < NUM_ROOTS, INDIM = %i, NUM_ROOTS = %i\n",init_dim,M);
        outfile->Printf("        setting INDIM = NUM_ROOTS\n");
        init_dim = M;
    }
    if ( init_dim > N ) {
        outfile->Printf("    WARNING: INDIM > N, INDIM = %i, N = %i\n",init_dim,N);
        outfile->Printf("        setting INDIM = N\n");
        init_dim = N;
    }
    if ( maxdim < init_dim ) { 
        outfile->Printf("    WARNING: MAXDIM < INDIM, MAXDIM = %i, INDIM = %i\n",maxdim,init_dim);
        outfile->Printf("        setting MAXDIM = INDIM\n");
        maxdim = init_dim;
    }
    if ( maxdim > N ) {
        outfile->Printf("    WARNING: MAXDIM > N, MAXDIM = %i, N = %i\n",maxdim,N);
        outfile->Printf("        setting MAXDIM = N\n");
        maxdim = N;
    }
    outfile->Printf("    No. requested roots:            %5i\n",M);
    outfile->Printf("    Problem dimension:              %5i\n",N);
    outfile->Printf("    Max subspace dimension:         %5i\n",maxdim);
    outfile->Printf("    Initial subspace dimension:     %5i\n",init_dim);

    // (approximate) diagonal of Hamiltonian
    double * Hdiag = build_hamiltonian_diagonals();
    double * Sdiag = build_overlap_diagonals();

    std::shared_ptr<Matrix> rel ( new Matrix(M,N) );
    std::shared_ptr<Matrix> rer ( new Matrix(M,N) );
    std::shared_ptr<Vector> revals ( new Vector(M) );

    double * revalp = revals->pointer();
    double ** rerp = rer->pointer();
    double ** relp = rel->pointer();

    bool use_residual_norm = true;
    double residual_norm = options_.get_double("RESIDUAL_NORM");

    // call eigensolver
    david->real_generalized_eigenvalue_problem(Hdiag,
                                               Sdiag,
                                               N,
                                               M,
                                               revalp,
                                               rerp,
                                               [&](int N, int maxdim, int L, double **Q, double **sigmar, double **sigmal) {
                                                      build_sigma_generalized(N,maxdim,L, Q, sigmar,sigmal);
                                               }, maxdim, init_dim, residual_norm, use_residual_norm);

    free(Hdiag);

    outfile->Printf("\n");
    outfile->Printf("    ==> QED-TDDFT energies <==\n");
    outfile->Printf("\n");
    outfile->Printf("    ");
    outfile->Printf("%5s","state");
    outfile->Printf(" %5s","type");
    outfile->Printf(" %20s","ex energy (Eh)");
    outfile->Printf(" %20s","energy (Eh)");
    outfile->Printf(" %20s","photon weight");
    outfile->Printf(" %10s","mu_x");
    outfile->Printf(" %10s","mu_y");
    outfile->Printf(" %10s","mu_z");
    outfile->Printf(" %10s","f");
    outfile->Printf("\n");

    for (int state = 0; state < M; state++) {

        // normalization
        double m_weight = rerp[state][2 * (oa*va + ob*vb)    ] * rerp[state][2 * (oa*va + ob*vb)    ];
        double n_weight = rerp[state][2 * (oa*va + ob*vb) + 1] * rerp[state][2 * (oa*va + ob*vb) + 1];
        double photon_weight = m_weight - n_weight;
        double x_weight = 0.0;
        double y_weight = 0.0;
        for (int j = 0; j < oa*va; j++) {
            x_weight += rerp[state][j        ] * rerp[state][j        ];
            y_weight += rerp[state][j + oa*va] * rerp[state][j + oa*va];
        }
        for (int j = 0; j < ob*vb; j++) {
            x_weight += rerp[state][j + 2*oa*va        ] * rerp[state][j + 2*oa*va        ];
            y_weight += rerp[state][j + 2*oa*va + ob*vb] * rerp[state][j + 2*oa*va + ob*vb];
        }
        double electron_weight = x_weight - y_weight;
        double nrm = photon_weight + electron_weight;
        std::string type;
        if ( nrm > 0.0 ) {
            type = "ex";
            electron_weight /= nrm;
            photon_weight /= nrm;
        }else {
            type = "de-ex";
            electron_weight /= nrm;
            photon_weight /= nrm;
            nrm = -nrm;
        }

        // oscillator strengths
        double mu_x_r = 0.0;
        double mu_y_r = 0.0;
        double mu_z_r = 0.0;
        double mu_x_l = 0.0;
        double mu_y_l = 0.0;
        double mu_z_l = 0.0;
        nrm = sqrt(nrm);
        for (int j = 0; j < oa; j++) {
            for (int b = 0; b < va; b++) {
                double cr = (rerp[state][j*va + b] + rerp[state][j*va + b + oa*va])/nrm;
                mu_x_r += dx[j][b + oa + ob] * cr;
                mu_y_r += dy[j][b + oa + ob] * cr;
                mu_z_r += dz[j][b + oa + ob] * cr;
            }
        }
        for (int j = 0; j < ob; j++) {
            for (int b = 0; b < vb; b++) {
                double cr = (rerp[state][j*vb + b + 2 * oa*va] + rerp[state][j*vb + b + ob*vb + 2 * oa*va])/nrm;
                mu_x_r += dx[j + oa][b + oa + ob + va] * cr;
                mu_y_r += dy[j + oa][b + oa + ob + va] * cr;
                mu_z_r += dz[j + oa][b + oa + ob + va] * cr;
            }
        }
        double w = revalp[state];
        double f = 2.0/3.0*w*(mu_x_r * mu_x_r + mu_y_r * mu_y_r + mu_z_r * mu_z_r);
        outfile->Printf("    %5i %5s %20.12lf %20.12lf %20.12lf %10.6lf %10.6lf %10.6lf %10.6lf\n", state, type.c_str(), w, energy_ + w, photon_weight, mu_x_r, mu_y_r, mu_z_r, f);
    }

    outfile->Printf("\n");
    outfile->Printf("    ==> QED-TDDFT significant amplitudes <==\n");
    for (int state = 0; state < M; state++) {

        double w = revalp[state];
        outfile->Printf("\n");
        outfile->Printf("    %5s","state");
        outfile->Printf(" %20s","ex energy (eV)");
        outfile->Printf(" %4s","id");
        outfile->Printf(" %12s","value");
        outfile->Printf(" %12s","transition");
        outfile->Printf(" %5s","type");
        outfile->Printf("\n");

        bool print = true;
        for (size_t  p = 0; p < N; p++) {
             double dum = rerp[state][p];
             //print only transitions' contribution with amplitude larger than 0.1
             if (fabs(dum)  > 0.1) {
                 if (print) {
                     outfile->Printf("    %5i", state);
                     outfile->Printf(" %20.12lf", w * pc_hartree2ev);
                     print = false;
                 }else {
                     outfile->Printf("         ");
                     outfile->Printf("                     ");
                 }
                 outfile->Printf(" %4i", p);
                 outfile->Printf(" %12.8lf", dum);
                 if (p < oa*va) {
                    size_t a = p % va;
                    size_t i = (p-a) / va;
                    outfile->Printf(" %4d -> %4d",i + 1, a + oa + 1);
                    outfile->Printf(" %5s\n", "Xa");
                 }else if ( p < 2 * oa*va) {
                    size_t a = p % va;
                    size_t i = (p-a) / va;
                    outfile->Printf(" %4d -> %4d",i + 1, a + oa + 1);
                    outfile->Printf(" %5s\n", "Ya");
                 }else if (p < 2*oa*va + ob*vb) {
                    size_t a = (p-2*oa*va) % vb;
                    size_t i = (p-2*oa*va-a) / vb;
                    outfile->Printf(" %4d -> %4d",i + 1, a + ob + 1);
                    outfile->Printf(" %5s\n", "Xb");
                 }else if ( p < 2 * (oa*va + ob*vb) ) {
                    size_t a = (p-2*oa*va) % vb;
                    size_t i = (p-2*oa*va-a) / vb;
                    outfile->Printf(" %4d -> %4d",i + 1, a + ob + 1);
                    outfile->Printf(" %5s\n", "Yb");
                 }else if (p == 2*(oa*va+ob*vb)) {
                     outfile->Printf("             ");
                     outfile->Printf(" %5s\n", "M");
                 }else {
                     outfile->Printf("             ");
                     outfile->Printf(" %5s\n", "N");
                 }
             }
        }
    }

    // add excitation energies to psi variables
    for (int state = 0; state < M; state++) {
        Process::environment.globals["QED-TDDFT ROOT 0 -> ROOT " + std::to_string(state+1) + " EXCITATION ENERGY"] = revalp[state];
    }

    return 0.0;
}

// QED-TDDFT:
// 
// |  A  B  g*  g* | ( X )     |  1  0  0  0 | ( X )
// |  B  A  g*  g* | ( Y ) = W |  0 -1  0  0 | ( Y )
// |  g  g  w   0  | ( M )     |  0  0  1  0 | ( M )
// |  g  g  0   w  | ( N )     |  0  0  0 -1 | ( N )
// 
// 
// 
// this function solves the generalized eigenvalue problem S c = 1/w H c
void PolaritonicUTDDFT::build_sigma_generalized(int N, int maxdim, int L, double **Q, double **sigmah, double **sigmas){

    for (int i = 0; i < N; i++) {
        for (int I = 0; I < maxdim; I++) {
            sigmah[i][I] = 0.0;
            sigmas[i][I] = 0.0;
        }
    }

    // TDDFT-PF:

    // sigmah_x =  Ax + By + g*m + g*n
    // sigmah_y =  Bx + Ay + g*m + g*n
    // sigmah_m =  gx + gy +  wm
    // sigmah_n =  gx + gy       +  wn

    double * Ax = (double*)malloc(N*L*sizeof(double));
    double * Bx = (double*)malloc(N*L*sizeof(double));
    double * Ay = (double*)malloc(N*L*sizeof(double));
    double * By = (double*)malloc(N*L*sizeof(double));
    double * x  = (double*)malloc(N*L*sizeof(double));
    double * y  = (double*)malloc(N*L*sizeof(double));

    memset((void*)Ax,'\0',N*L*sizeof(double));
    memset((void*)Bx,'\0',N*L*sizeof(double));
    memset((void*)Ay,'\0',N*L*sizeof(double));
    memset((void*)By,'\0',N*L*sizeof(double));
    memset((void*)x, '\0',N*L*sizeof(double));
    memset((void*)y, '\0',N*L*sizeof(double));

    double *m = (double*)malloc(L*sizeof(double));
    double *gm = (double*)malloc(N*L*sizeof(double));
    double *sigma_m_h = (double*)malloc(L*sizeof(double));
    double *sigma_m_s = (double*)malloc(L*sizeof(double));
    double *n = (double*)malloc(L*sizeof(double));
    double *gn = (double*)malloc(N*L*sizeof(double));

    memset((void*)m,'\0',L*sizeof(double));
    memset((void*)gm,'\0',N*L*sizeof(double));
    memset((void*)sigma_m_h,'\0',L*sizeof(double));
    memset((void*)sigma_m_s,'\0',L*sizeof(double));
    memset((void*)n,'\0',L*sizeof(double));
    memset((void*)gn,'\0',N*L*sizeof(double));

    int oa = nalpha_;
    int va = nso_ - oa;
    int ob = nbeta_;
    int vb = nso_ - ob;

    // unpack input vectors
    for (int I = 0; I < L; I++) {
        for (int j = 0; j < oa*va; j++) {
            x[I*N+j] = Q[I][j];
            y[I*N+j] = Q[I][j + oa*va];
        }
        for (int j = 0; j < ob*vb; j++) {
            x[I*N + oa*va + j] = Q[I][2*oa*va + j];
            y[I*N + oa*va + j] = Q[I][2*oa*va + j + ob*vb];
        }
        m[I] = Q[I][2*(oa*va+ob*vb)];
        n[I] = Q[I][2*(oa*va+ob*vb)+1];
    }

    build_Au_Bu(N, L, x, Ax, Bx);
    build_Au_Bu(N, L, y, Ay, By);
    build_gm(N, L, m, gm);
    build_gm(N, L, n, gn);
    build_sigma_m(N, L, x, y, m, sigma_m_h, sigma_m_s); // sigma_m_s need to be overwritten

    // x, y
    for (int i = 0; i < oa*va; i++) {
        for (int I = 0; I < L; I++) {
            sigmah[i        ][I] =  Ax[I*N + i] + By[I*N + i] + gm[I*N + i] + gn[I*N + i];
            sigmah[i + oa*va][I] =  Bx[I*N + i] + Ay[I*N + i] + gm[I*N + i] + gn[I*N + i];

            sigmas[i        ][I] =  x[I*N + i];
            sigmas[i + oa*va][I] = -y[I*N + i];
        }
    }
    for (int i = 0; i < ob*vb; i++) {
        for (int I = 0; I < L; I++) {
            sigmah[2*oa*va + i        ][I] =  Ax[I*N + oa*va + i] + By[I*N + oa*va + i] + gm[I*N + oa*va + i] + gn[I*N + oa*va + i];
            sigmah[2*oa*va + i + ob*vb][I] =  Bx[I*N + oa*va + i] + Ay[I*N + oa*va + i] + gm[I*N + oa*va + i] + gn[I*N + oa*va + i];

            sigmas[2*oa*va + i        ][I] =  x[I*N + oa*va + i];
            sigmas[2*oa*va + i + ob*vb][I] = -y[I*N + oa*va + i];
        }
    }
    // m, n
    for (int I = 0; I < L; I++) {
        sigmah[2*(oa*va + ob*vb)    ][I] =  sigma_m_h[I] + cavity_frequency_[2] * m[I];
        sigmas[2*(oa*va + ob*vb)    ][I] =  m[I]; 

        sigmah[2*(oa*va + ob*vb) + 1][I] =  sigma_m_h[I] + cavity_frequency_[2] * n[I]; 
        sigmas[2*(oa*va + ob*vb) + 1][I] = -n[I];
    }

    free(x);
    free(y);
    free(Ax);
    free(Bx);
    free(Ay);
    free(By);
    free(m);
    free(sigma_m_h);
    free(sigma_m_s);
    free(gm);
    free(n);
    free(gn);
}

// build Au and Bu according to Yihan's formalism that includes only (1) single 
// excitations and de-excitations and (2) photon excitations and de-excitations, 
// but no coupled photon / electron excitations. Also, this formalism only works
// for n_photon_states_ = 2

// Yihan calls this TDDFT-JC
// note that Au will contain contribution from m

void PolaritonicUTDDFT::build_Au_Bu(int N, int L, double *u, double *Au, double *Bu){

    if ( n_photon_states_ > 2 ) {
        throw PsiException("qed-utddft only works for n_photon_states <= 2",__FILE__,__LINE__);
    }

    int oa = nalpha_;
    int ob = nbeta_;
    int va = nso_ - oa;
    int vb = nso_ - ob;

    double * c = (double*)malloc(N*sizeof(double));

    double ** cap = Ca_->pointer();
    double ** cbp = Cb_->pointer();

    std::vector<SharedMatrix>& C_left  = jk_->C_left();
    std::vector<SharedMatrix>& C_right = jk_->C_right();
    C_left.clear();
    C_right.clear();

    double coupling_factor_x = cavity_frequency_[0] * cavity_coupling_strength_[0];
    double coupling_factor_y = cavity_frequency_[1] * cavity_coupling_strength_[1];
    double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];

    double lambda_x = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
    double lambda_y = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

    double ** dx = Dipole_x_->pointer();
    double ** dy = Dipole_y_->pointer();
    double ** dz = Dipole_z_->pointer();

    std::vector<std::shared_ptr<Matrix> > Vx;
    std::vector<std::shared_ptr<Matrix> > Dx;

    // J/K-like contributions
    for (int I = 0; I < L; I++) {

        // unpack a vector
        for (int j = 0; j < N; j++) {
            c[j] = u[I*N+j];
        }

        // push density matrices on JK object
        // we need density matrices for alpha as well as
        // alpha plus a photon

        // Da(mu,nu) = ca(j,b) Ca(mu,j) Ca(nu,b) = Ca(mu,j) Ca'(nu,j)
        // Ca'(nu,j) = ca(j,b) Ca(nu,b)

        // etc.

        // singles

        // alpha
        std::shared_ptr<Matrix> cra (new Matrix(Ca_) );
        std::shared_ptr<Matrix> cla (new Matrix(Ca_) );

        double ** clap = cla->pointer();
        double ** crap = cra->pointer();

        cra->zero();
        cla->zero();

        for (int mu = 0; mu < nso_; mu++) {
            for (int i = 0; i < oa; i++) {

                // left is plain orbitals
                clap[mu][i] = cap[mu][i];

                // right is modified orbitals
                double dum = 0.0;
                for (int a = 0; a < va; a++) {
                    int ia = i * va + a;
                    dum += cap[mu][a+oa] * c[ia];
                }
                crap[mu][i] = dum;

            }
        }

        // push alpha orbitals onto JK object
        C_left.push_back(cla);
        C_right.push_back(cra);

        // build pseudo densities for xc contribution
        if (needs_xc_) {

            auto Dx_a = linalg::doublet(cla, cra, false, true);

            Vx.push_back(std::make_shared<Matrix>("Vax temp", Dx_a->rowspi(), Dx_a->colspi(), Dx_a->symmetry()));

            Dx.push_back(Dx_a);
        }

        // beta
        std::shared_ptr<Matrix> crb (new Matrix(Cb_) );
        std::shared_ptr<Matrix> clb (new Matrix(Cb_) );

        double ** clbp = clb->pointer();
        double ** crbp = crb->pointer();

        crb->zero();
        clb->zero();

        for (int mu = 0; mu < nso_; mu++) {
            for (int i = 0; i < ob; i++) {

                // left is plain orbitals
                clbp[mu][i] = cbp[mu][i];

                // right is modified orbitals
                double dum = 0.0;
                for (int a = 0; a < vb; a++) {
                    int ia = i * vb + a;
                    dum += cbp[mu][a+ob] * c[oa*va + ia];
                }
                crbp[mu][i] = dum;

            }
        }

        // push beta orbitals onto JK object
        C_left.push_back(clb);
        C_right.push_back(crb);

        // build pseudo densities for xc contribution
        if (needs_xc_) {

            auto Dx_b = linalg::doublet(clb, crb, false, true);

            Vx.push_back(std::make_shared<Matrix>("Vbx temp", Dx_b->rowspi(), Dx_b->colspi(), Dx_b->symmetry()));

            Dx.push_back(Dx_b);
        }

    }

    // form J/K
    jk_->compute();

    // form xc contributions
    if ( needs_xc_ ) {
        potential_->compute_Vx(Dx, Vx);
    }

    // now, accumulate sigma vectors

    // offset for j/k matrices
    int count = 0;

    double * ea = epsilon_a_->pointer();
    double * eb = epsilon_b_->pointer();

    double * tmpa_a = (double*)malloc(oa*va*sizeof(double));
    double * tmpa_b = (double*)malloc(oa*oa*sizeof(double));
    double * tmpb_a = (double*)malloc(ob*vb*sizeof(double));
    double * tmpb_b = (double*)malloc(ob*ob*sizeof(double));

    for (int I = 0; I < L; I++) {

        // unpack a vector
        for (int j = 0; j < N; j++) {
            c[j] = u[I*N+j];
        }

        // a <- a+b coulomb
        std::shared_ptr<Matrix> sa (new Matrix(jk_->J()[count]));
        sa->add(jk_->J()[count+1]);

        // b <- a+b coulomb
        std::shared_ptr<Matrix> sb (new Matrix(jk_->J()[count]));
        sb->add(jk_->J()[count+1]);

        // xc?
        // a <- a
        // b <- b
        if ( needs_xc_ ) {
            sa->axpy(1.0, Vx[count]);
            sb->axpy(1.0, Vx[count+1]);
        }

        // exact exchange?
        // a <- a
        // b <- b
        if (is_x_hybrid_) {
            sa->axpy(-x_alpha_,jk_->K()[count]);
            sb->axpy(-x_alpha_,jk_->K()[count+1]);
        }

        // LRC functional?
        // a <- a
        // b <- b
        if (is_x_lrc_) {
            double beta = 1.0 - x_alpha_;
            sa->axpy(-beta,jk_->wK()[count]);
            sb->axpy(-beta,jk_->wK()[count+1]);
        }

        // update counter
        count += 2;

        // transform jk, e.g., j(a,i) = j(mu,nu) c(mu,a) c(nu,i)

        sa->transform(Ca_);
        sb->transform(Cb_);

        double ** sap = sa->pointer();
        double ** sbp = sb->pointer();

        // alpha
        for (int i = 0; i < oa; i++) {
            for (int a = 0; a < va; a++) {
                int ia = i * va + a;

                Au[I*N+ia] = sap[i][a+oa];
                Bu[I*N+ia] = sap[a+oa][i];
            }
        }
        // beta
        for (int i = 0; i < ob; i++) {
            for (int a = 0; a < vb; a++) {
                int ia = i * vb + a;

                Au[I*N+ia + oa*va] = sbp[i][a+ob];
                Bu[I*N+ia + oa*va] = sbp[a+ob][i];
            }
        }

        // singles diagonals (alpha)
        for (int i = 0; i < oa; i++) {
            for (int a = 0; a < va; a++) {
                int ia = i * va + a;
                Au[I*N+ia] += c[ia] * (ea[a+oa] - ea[i]);
            }
        }
        // singles diagonals (beta)
        for (int i = 0; i < ob; i++) {
            for (int a = 0; a < vb; a++) {
                int ia = i * vb + a;
                Au[I*N+ia + oa*va] += c[ia + oa*va] * (eb[a+ob] - eb[i]);
            }
        }

        // dipole self energy (for singles)
        // J-like contribution from dipole self energy (alpha)
        double dipole_Ja = 0.0;
        for (int i = 0; i < oa; i++) {
            for (int a = 0; a < va; a++) {
                int ia = i * va + a;
                dipole_Ja += c[ia] * dz[i][a + oa + ob];
            }
        }
        double dipole_Jb = 0.0;
        for (int i = 0; i < ob; i++) {
            for (int a = 0; a < vb; a++) {
                int ia = i * vb + a;
                dipole_Jb += c[ia + oa*va] * dz[i + oa][a + oa + ob + va];
            }
        }

        // intermediate for K-like contribution from dipole self energy

        // ignore if following QED-TDDFT outlined in J. Chem. Phys. 155, 064107 (2021)
        memset((void*)tmpa_a,'\0',oa*va*sizeof(double));
        memset((void*)tmpa_b,'\0',oa*oa*sizeof(double));
        memset((void*)tmpb_a,'\0',ob*vb*sizeof(double));
        memset((void*)tmpb_b,'\0',ob*ob*sizeof(double));

        if ( options_.get_bool("QED_USE_RELAXED_ORBITALS") ) {

            // A term (alpha)
            for (int i = 0; i < oa; i++) {
                for (int a = 0; a < va; a++) {
                    double dum_a = 0.0;
                    for (int j = 0; j < oa; j++) {
                        int ja = j * va + a;
                        dum_a += c[ja] * dz[i][j];
                    }
                    tmpa_a[a*oa+i] = dum_a;
                }
            }

            // A term (beta)
            for (int i = 0; i < ob; i++) {
                for (int a = 0; a < vb; a++) {
                    double dum_a = 0.0;
                    for (int j = 0; j < ob; j++) {
                        int ja = j * vb + a;
                        dum_a += c[ja + oa*va] * dz[i + oa][j + oa];
                    }
                    tmpb_a[a*ob+i] = dum_a;
                }
            }

            // B term (alpha)
            for (int i = 0; i < oa; i++) {
                for (int j = 0; j < oa; j++) {
                    double dum_b = 0.0;
                    for (int b = 0; b < va; b++) {
                        int jb = j * va + b;
                        dum_b += c[jb] * dz[b + oa + ob][i];
                    }
                    tmpa_b[i*oa+j] = dum_b;
                }
            }

            // B term (beta)
            for (int i = 0; i < ob; i++) {
                for (int j = 0; j < ob; j++) {
                    double dum_b = 0.0;
                    for (int b = 0; b < vb; b++) {
                        int jb = j * vb + b;
                        dum_b += c[jb + oa*va] * dz[b + oa + ob + va][i + oa];
                    }
                    tmpb_b[i*ob+j] = dum_b;
                }
            }
        }

        // alpha
        for (int i = 0; i < oa; i++) {
            for (int a = 0; a < va; a++) {

                double dipole_J_ia = (dipole_Ja + dipole_Jb) * dz[i][a + oa + ob];

                // A term
                double dipole_Ka_A = 0.0;
                for (int b = 0; b < va; b++) {
                    dipole_Ka_A += tmpa_a[b*oa + i] * dz[a + oa + ob][b + oa + ob];
                }
                // B term
                double dipole_Ka_B = 0.0;
                for (int j = 0; j < oa; j++) {
                    dipole_Ka_B += tmpa_b[i*oa + j] * dz[a + oa + ob][j];
                }

                int ia = i * va + a;

                Au[I*N + ia] += lambda_z * lambda_z * (dipole_J_ia - dipole_Ka_A);
                Bu[I*N + ia] += lambda_z * lambda_z * (dipole_J_ia - dipole_Ka_B);

            }
        }
        // beta
        for (int i = 0; i < ob; i++) {
            for (int a = 0; a < vb; a++) {

                double dipole_J_ia = (dipole_Ja + dipole_Jb) * dz[i + oa][a + oa + ob + va];

                // A term
                double dipole_Kb_A = 0.0;
                for (int b = 0; b < vb; b++) {
                    dipole_Kb_A += tmpb_a[b*ob + i] * dz[a + oa + ob + va][b + oa + ob + va];
                }
                // B term
                double dipole_Kb_B = 0.0;
                for (int j = 0; j < ob; j++) {
                    dipole_Kb_B += tmpb_b[i*ob + j] * dz[a + oa + ob + va][j + oa];
                }

                int ia = i * vb + a;

                Au[I*N + ia + oa*va] += lambda_z * lambda_z * (dipole_J_ia - dipole_Kb_A);
                Bu[I*N + ia + oa*va] += lambda_z * lambda_z * (dipole_J_ia - dipole_Kb_B);
            }
        }
    }
    free(tmpa_a);
    free(tmpa_b);
    free(tmpb_a);
    free(tmpb_b);
    free(c);
}

void PolaritonicUTDDFT::build_Au_Bu_response(int N, double *u, double *ABu){

    if ( n_photon_states_ > 2 ) {
        throw PsiException("qed-utddft only works for n_photon_states <= 2",__FILE__,__LINE__);
    }

    int oa = nalpha_;
    int ob = nbeta_;
    int va = nso_ - oa;
    int vb = nso_ - ob;

    double ** cap = Ca_->pointer();
    double ** cbp = Cb_->pointer();

    std::vector<SharedMatrix>& C_left  = jk_->C_left();
    std::vector<SharedMatrix>& C_right = jk_->C_right();
    C_left.clear();
    C_right.clear();

    double coupling_factor_x = cavity_frequency_[0] * cavity_coupling_strength_[0];
    double coupling_factor_y = cavity_frequency_[1] * cavity_coupling_strength_[1];
    double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];

    double lambda_x = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
    double lambda_y = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

    double ** dx = Dipole_x_->pointer();
    double ** dy = Dipole_y_->pointer();
    double ** dz = Dipole_z_->pointer();

    std::vector<std::shared_ptr<Matrix> > Vx;
    std::vector<std::shared_ptr<Matrix> > Dx;

    // J/K-like contributions

    // push density matrices on JK object
    // we need density matrices for alpha as well as
    // alpha plus a photon

    // Da(mu,nu) = ca(j,b) Ca(mu,j) Ca(nu,b) = Ca(mu,j) Ca'(nu,j)
    // Ca'(nu,j) = ca(j,b) Ca(nu,b)

    // etc.

    // singles

    // alpha
    std::shared_ptr<Matrix> cra (new Matrix(Ca_) );
    std::shared_ptr<Matrix> cla (new Matrix(Ca_) );

    double ** clap = cla->pointer();
    double ** crap = cra->pointer();

    cra->zero();
    cla->zero();

    for (int mu = 0; mu < nso_; mu++) {
        for (int i = 0; i < oa; i++) {

            // left is plain orbitals
            clap[mu][i] = cap[mu][i];

            // right is modified orbitals
            double dum = 0.0;
            for (int a = 0; a < va; a++) {
                int ia = i * va + a;
                dum += cap[mu][a+oa] * u[ia];
            }
            crap[mu][i] = dum;

        }
    }

    // push alpha orbitals onto JK object
    C_left.push_back(cla);
    C_right.push_back(cra);

    // build pseudo densities for xc contribution
    if (needs_xc_) {

        auto Dx_a = linalg::doublet(cla, cra, false, true);

        Vx.push_back(std::make_shared<Matrix>("Vax temp", Dx_a->rowspi(), Dx_a->colspi(), Dx_a->symmetry()));

        Dx.push_back(Dx_a);
    }

    // beta
    std::shared_ptr<Matrix> crb (new Matrix(Cb_) );
    std::shared_ptr<Matrix> clb (new Matrix(Cb_) );

    double ** clbp = clb->pointer();
    double ** crbp = crb->pointer();

    crb->zero();
    clb->zero();

    for (int mu = 0; mu < nso_; mu++) {
        for (int i = 0; i < ob; i++) {

            // left is plain orbitals
            clbp[mu][i] = cbp[mu][i];

            // right is modified orbitals
            double dum = 0.0;
            for (int a = 0; a < vb; a++) {
                int ia = i * vb + a;
                dum += cbp[mu][a+ob] * u[oa*va + ia];
            }
            crbp[mu][i] = dum;

        }
    }

    // push beta orbitals onto JK object
    C_left.push_back(clb);
    C_right.push_back(crb);

    // build pseudo densities for xc contribution
    if (needs_xc_) {

        auto Dx_b = linalg::doublet(clb, crb, false, true);

        Vx.push_back(std::make_shared<Matrix>("Vbx temp", Dx_b->rowspi(), Dx_b->colspi(), Dx_b->symmetry()));

        Dx.push_back(Dx_b);
    }

    // form J/K
    jk_->compute();

    // form xc contributions
    if ( needs_xc_ ) {
        potential_->compute_Vx(Dx, Vx);
    }

    // now, accumulate sigma vectors

    // offset for j/k matrices
    int count = 0;

    double * ea = epsilon_a_->pointer();
    double * eb = epsilon_b_->pointer();

    double * tmpa_a = (double*)malloc(oa*va*sizeof(double));
    double * tmpa_b = (double*)malloc(oa*oa*sizeof(double));
    double * tmpb_a = (double*)malloc(ob*vb*sizeof(double));
    double * tmpb_b = (double*)malloc(ob*ob*sizeof(double));

    // a <- a+b coulomb
    std::shared_ptr<Matrix> sa (new Matrix(jk_->J()[count]));
    sa->add(jk_->J()[count+1]);

    // b <- a+b coulomb
    std::shared_ptr<Matrix> sb (new Matrix(jk_->J()[count]));
    sb->add(jk_->J()[count+1]);

    // xc?
    // a <- a
    // b <- b
    if ( needs_xc_ ) {
        sa->axpy(1.0, Vx[count]);
        sb->axpy(1.0, Vx[count+1]);
    }

    // exact exchange?
    // a <- a
    // b <- b
    if (is_x_hybrid_) {
        sa->axpy(-x_alpha_,jk_->K()[count]);
        sb->axpy(-x_alpha_,jk_->K()[count+1]);
    }

    // LRC functional?
    // a <- a
    // b <- b
    if (is_x_lrc_) {
        double beta = 1.0 - x_alpha_;
        sa->axpy(-beta,jk_->wK()[count]);
        sb->axpy(-beta,jk_->wK()[count+1]);
    }

    // update counter
    count += 2;

    // transform jk, e.g., j(a,i) = j(mu,nu) c(mu,a) c(nu,i)

    sa->transform(Ca_);
    sb->transform(Cb_);

    double ** sap = sa->pointer();
    double ** sbp = sb->pointer();

    // alpha
    C_DCOPY(nmo_*nmo_, &sap[0][0], 1, &ABu[0], 1);

    // beta
    C_DCOPY(nmo_*nmo_, &sbp[0][0], 1, &ABu[nmo_*nmo_], 1);

    // dipole self energy (for singles)

    // J-like contribution from dipole self energy (alpha)
    double dipole_Ja = 0.0;
    for (int i = 0; i < oa; i++) {
        for (int a = 0; a < va; a++) {
            int ia = i * va + a;
            dipole_Ja += u[ia] * dz[i][a + oa + ob];
        }
    }
    double dipole_Jb = 0.0;
    for (int i = 0; i < ob; i++) {
        for (int a = 0; a < vb; a++) {
            int ia = i * vb + a;
            dipole_Jb += u[ia + oa*va] * dz[i + oa][a + oa + ob + va];
        }
    }

    // intermediate for K-like contribution from dipole self energy

    // ignore if following QED-TDDFT outlined in J. Chem. Phys. 155, 064107 (2021)
    memset((void*)tmpa_a,'\0',oa*va*sizeof(double));
    memset((void*)tmpa_b,'\0',oa*oa*sizeof(double));
    memset((void*)tmpb_a,'\0',ob*vb*sizeof(double));
    memset((void*)tmpb_b,'\0',ob*ob*sizeof(double));

    if ( options_.get_bool("QED_USE_RELAXED_ORBITALS") ) {

        // A term (alpha)
        for (int i = 0; i < oa; i++) {
            for (int a = 0; a < va; a++) {
                double dum_a = 0.0;
                for (int j = 0; j < oa; j++) {
                    int ja = j * va + a;
                    dum_a += u[ja] * dz[i][j];
                }
                tmpa_a[a*oa+i] = dum_a;
            }
        }

        // A term (beta)
        for (int i = 0; i < ob; i++) {
            for (int a = 0; a < vb; a++) {
                double dum_a = 0.0;
                for (int j = 0; j < ob; j++) {
                    int ja = j * vb + a;
                    dum_a += u[ja + oa*va] * dz[i + oa][j + oa];
                }
                tmpb_a[a*ob+i] = dum_a;
            }
        }

        // B term (alpha)
        for (int i = 0; i < oa; i++) {
            for (int j = 0; j < oa; j++) {
                double dum_b = 0.0;
                for (int b = 0; b < va; b++) {
                    int jb = j * va + b;
                    dum_b += u[jb] * dz[b + oa + ob][i];
                }
                tmpa_b[i*oa+j] = dum_b;
            }
        }

        // B term (beta)
        for (int i = 0; i < ob; i++) {
            for (int j = 0; j < ob; j++) {
                double dum_b = 0.0;
                for (int b = 0; b < vb; b++) {
                    int jb = j * vb + b;
                    dum_b += u[jb + oa*va] * dz[b + oa + ob + va][i + oa];
                }
                tmpb_b[i*ob+j] = dum_b;
            }
        }
    }

    // ia and ai parts for second derivatives

    // alpha
    for (int i = 0; i < oa; i++) {
        for (int a = 0; a < va; a++) {

            double dipole_J_ia = (dipole_Ja + dipole_Jb) * dz[i][a + oa + ob];

            // A term
            double dipole_Ka_A = 0.0;
            for (int b = 0; b < va; b++) {
                dipole_Ka_A += tmpa_a[b*oa + i] * dz[a + oa + ob][b + oa + ob];
            }
            // B term
            double dipole_Ka_B = 0.0;
            for (int j = 0; j < oa; j++) {
                dipole_Ka_B += tmpa_b[i*oa + j] * dz[a + oa + ob][j];
            }

            int ia = i * nmo_ + (a + oa);
            int ai = (a + oa) * nmo_ + i;

            ABu[ia] += lambda_z * lambda_z * (dipole_J_ia - dipole_Ka_A);
            ABu[ai] += lambda_z * lambda_z * (dipole_J_ia - dipole_Ka_B);
        }
    }
    // beta
    for (int i = 0; i < ob; i++) {
        for (int a = 0; a < vb; a++) {

            double dipole_J_ia = (dipole_Ja + dipole_Jb) * dz[i + oa][a + oa + ob + va];

            // A term
            double dipole_Kb_A = 0.0;
            for (int b = 0; b < vb; b++) {
                dipole_Kb_A += tmpb_a[b*ob + i] * dz[a + oa + ob + va][b + oa + ob + va];
            }
            // B term
            double dipole_Kb_B = 0.0;
            for (int j = 0; j < ob; j++) {
                dipole_Kb_B += tmpb_b[i*ob + j] * dz[a + oa + ob + va][j + oa];
            }

            int ia = i * nmo_ + (a + ob) + nmo_*nmo_;
            int ai = (a + ob) * nmo_ + i + nmo_*nmo_;

            ABu[ia] += lambda_z * lambda_z * (dipole_J_ia - dipole_Kb_A);
            ABu[ai] += lambda_z * lambda_z * (dipole_J_ia - dipole_Kb_B);
        }
    }

    free(tmpa_a);
    free(tmpa_b);
    free(tmpb_a);
    free(tmpb_b);
}

void PolaritonicUTDDFT::build_gm(int N, int L, double *m, double *gm) {

    if ( n_photon_states_ > 2 ) {
        throw PsiException("qed-utddft only works for n_photon_states <= 2",__FILE__,__LINE__);
    }

    int oa = nalpha_;
    int ob = nbeta_;
    int va = nso_ - oa;
    int vb = nso_ - ob;

    double coupling_factor_x = cavity_frequency_[0] * cavity_coupling_strength_[0];
    double coupling_factor_y = cavity_frequency_[1] * cavity_coupling_strength_[1];
    double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];

    double lambda_x = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
    double lambda_y = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

    double ** dx = Dipole_x_->pointer();
    double ** dy = Dipole_y_->pointer();
    double ** dz = Dipole_z_->pointer();

    for (int I = 0; I < L; I++) {

        // couple |0,1> to |ia,0> (alpha)
        for (int i = 0; i < oa; i++) {
            for (int a = 0; a < va; a++) {
                int ia = i * va + a;

                // <ia| H |0,1>
                gm[I*N + ia] = -coupling_factor_z * dz[i][a + oa + ob] * m[I];
            }
        }
        // couple |0,1> to |ia,0> (beta)
        for (int i = 0; i < ob; i++) {
            for (int a = 0; a < vb; a++) {
                int ia = i * vb + a;

                // <ia| H |0,1>
                gm[I*N + ia + oa*va] = -coupling_factor_z * dz[i + oa][a + oa + ob + va] * m[I];
            }
        }
    }
}

void PolaritonicUTDDFT::build_sigma_m(int N, int L, double *x, double *y, double *m, double *sigma_m_r, double *sigma_m_l) {

    if ( n_photon_states_ > 2 ) {
        throw PsiException("qed-utddft only works for n_photon_states <= 2",__FILE__,__LINE__);
    }

    int oa = nalpha_;
    int ob = nbeta_;
    int va = nso_ - oa;
    int vb = nso_ - ob;

    double coupling_factor_x = cavity_frequency_[0] * cavity_coupling_strength_[0];
    double coupling_factor_y = cavity_frequency_[1] * cavity_coupling_strength_[1];
    double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];

    double lambda_x = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
    double lambda_y = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

    double ** dx = Dipole_x_->pointer();
    double ** dy = Dipole_y_->pointer();
    double ** dz = Dipole_z_->pointer();

    for (int I = 0; I < L; I++) {

        // |0,1> diagonal
        sigma_m_r[I] = 0.0;//cavity_frequency_[2] * m[I];
        sigma_m_l[I] = 0.0;//cavity_frequency_[2] * m[I];

        // couple |0,1> to |ia,0> (alpha)
        for (int i = 0; i < oa; i++) {
            for (int a = 0; a < va; a++) {

                int ia = i * va + a;

                // <ia| H |0,1>
                double factor = -coupling_factor_z * dz[i][a + oa + ob];
                    
                sigma_m_r[I] += factor * ( x[I*N + ia] + y[I*N + ia] );
                sigma_m_l[I] += factor * ( x[I*N + ia] - y[I*N + ia] );
            }
        }
        // couple |0,1> to |ia,0> (beta)
        for (int i = 0; i < ob; i++) {
            for (int a = 0; a < vb; a++) {

                int ia = i * vb + a;

                // <ia| H |0,1>
                double factor = -coupling_factor_z * dz[i + oa][a + oa + ob + va];
                    
                sigma_m_r[I] += factor * ( x[I*N + ia + oa*va] + y[I*N + ia + oa*va] );
                sigma_m_l[I] += factor * ( x[I*N + ia + oa*va] - y[I*N + ia + oa*va] );
            }
        }
    }
}

// build approximation to diagonal of A matrix in basis that excludes |0,0> and de-excitation |0,-n>
double * PolaritonicUTDDFT::build_hamiltonian_diagonals(){

    int oa = nalpha_;
    int ob = nbeta_;
    int va = nso_ - oa;
    int vb = nso_ - ob;

    int dim = 2*(oa*va+ob*vb)+2;

    double * H = (double*)malloc(dim*sizeof(double));
    memset((void*)H,'\0',dim*sizeof(double));

    double * ea = epsilon_a_->pointer();
    double * eb = epsilon_b_->pointer();

    for (int i = 0; i < oa; i++) {
        for (int a = 0; a < va; a++) {
            int ia = i * va + a;
            H[ia]         = ea[a+oa] - ea[i];
            H[oa*va + ia] = ea[a+oa] - ea[i];
        }
    }

    for (int i = 0; i < ob; i++) {
        for (int a = 0; a < vb; a++) {
            int ia = i * vb + a;
            H[2*oa*va + ia]         = eb[a+ob] - eb[i];
            H[2*oa*va + ob*vb + ia] = eb[a+ob] - eb[i];
        }
    }

    // now, |0,1> diagonals
    if ( n_photon_states_ > 2 ) {
        throw PsiException("qed-utddft only works for n_photon_states <= 2",__FILE__,__LINE__);
    }

    int off = 2 * (oa * va + ob * vb);
    H[off  ] = cavity_frequency_[2];
    H[off+1] = cavity_frequency_[2];

    return H;
}

double * PolaritonicUTDDFT::build_overlap_diagonals(){

    int oa = nalpha_;
    int ob = nbeta_;
    int va = nso_ - oa;
    int vb = nso_ - ob;

    int dim = 2*(oa*va+ob*vb)+2;

    double * S = (double*)malloc(dim*sizeof(double));
    memset((void*)S,'\0',dim*sizeof(double));

    double * ea = epsilon_a_->pointer();

    for (int i = 0; i < oa; i++) {
        for (int a = 0; a < va; a++) {
            int ia = i * va + a;
            S[ia]         =  1.0;
            S[oa*va + ia] = -1.0;
        }
    }
    for (int i = 0; i < ob; i++) {
        for (int a = 0; a < vb; a++) {
            int ia = i * vb + a;
            S[2*oa*va + ia]         =  1.0;
            S[2*oa*va + ob*vb + ia] = -1.0;
        }
    }

    // now, |0,1> diagonals
    if ( n_photon_states_ > 2 ) {
        throw PsiException("qed-utddft only works for n_photon_states <= 2",__FILE__,__LINE__);
    }

    int off = 2 * (oa * va + ob * vb);
    S[off  ] =  1.0;
    S[off+1] = -1.0;

    return S;
}

} // End namespaces

