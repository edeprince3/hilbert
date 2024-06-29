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

    // ensure scf_type df
    if ( options_.get_str("SCF_TYPE") != "DF" && options_.get_str("SCF_TYPE") != "CD" ) {
        throw PsiException("polaritonic utddft only works with scf_type df for now",__FILE__,__LINE__);
    }

    // ensure running in c1 symmetry
    if ( reference_wavefunction_->nirrep() > 1 ) {
        throw PsiException("polaritonic utddft only works with c1 symmetry for now.",__FILE__,__LINE__);
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

    // print the ks information
    //potential_->print_header();

    is_x_lrc_    = functional->is_x_lrc();
    is_x_hybrid_ = functional->is_x_hybrid();
    x_omega_     = functional->x_omega();
    x_alpha_     = functional->x_alpha();
    needs_xc_    = functional->needs_xc();

    if ( options_.get_str("SCF_TYPE") == "DF" ) {

        // get auxiliary basis:
        std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");

        // total number of auxiliary basis functions
        //nQ = auxiliary->nbf();

        std::shared_ptr<DiskDFJK> myjk = (std::shared_ptr<DiskDFJK>)(new DiskDFJK(primary,auxiliary,options_));

        // memory for jk (say, 80% of what is available)
        myjk->set_memory(0.8 * memory_);

        // integral cutoff
        myjk->set_cutoff(options_.get_double("INTS_TOLERANCE"));

        // Do J/K/wK?
        is_x_lrc_  = functional->is_x_lrc();
        is_x_hybrid_  = functional->is_x_hybrid();
        //if ( options_["IP_FITTING"].has_changed() ) {
        //    if ( options_.get_bool("IP_FITTING") ) {
        //        is_x_lrc = true;
        //    }
        //}
        x_omega_ = functional->x_omega();
        if ( options_["DFT_OMEGA"].has_changed() ) {
            x_omega_ = options_.get_double("DFT_OMEGA");
        }

        myjk->set_do_J(true);
        myjk->set_do_K(is_x_hybrid_);
        myjk->set_do_wK(is_x_lrc_);
        myjk->set_omega(x_omega_);

        myjk->initialize();

        jk_ = myjk;

    }else if ( options_.get_str("SCF_TYPE") == "CD" ) {

        std::shared_ptr<CDJK> myjk = (std::shared_ptr<CDJK>)(new CDJK(primary,options_,options_.get_double("CHOLESKY_TOLERANCE")));

        // memory for jk (say, 80% of what is available)
        myjk->set_memory(0.8 * memory_);

        // integral cutoff
        myjk->set_cutoff(options_.get_double("INTS_TOLERANCE"));

        // Do J/K/wK?
        is_x_lrc_  = functional->is_x_lrc();
        is_x_hybrid_  = functional->is_x_hybrid();
        //if ( options_["IP_FITTING"].has_changed() ) {
        //    if ( options_.get_bool("IP_FITTING") ) {
        //        is_x_lrc = true;
        //    }
        //}
        x_omega_ = functional->x_omega();
        if ( options_["DFT_OMEGA"].has_changed() ) {
            x_omega_ = options_.get_double("DFT_OMEGA");
        }

        myjk->set_do_J(true);
        myjk->set_do_K(is_x_hybrid_);
        myjk->set_do_wK(is_x_lrc_);
        myjk->set_omega(x_omega_);

        myjk->initialize();

        jk_ = myjk;
    }

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
    outfile->Printf("    No. states:                     %5i\n",N);
    outfile->Printf("    No. roots:                      %5i\n",M);
    outfile->Printf("    Max subspace dim:               %5i\n",maxdim);
    outfile->Printf("    Initial subspace dim:           %5i\n",init_dim);

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
    outfile->Printf("    QED-TDDFT energies:\n");
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
/*
// TODO: fix
    for (int state = 0; state < M; state++) {

        double w = revalp[state];//sqrt(revalp[state]);
        outfile->Printf("state =%5i eig = %20.12lf eV\n", state,w* 27.21138);

        for (size_t  p = 0; p < N; p++) {
             double dum = rerp[state][p];
             //print only transitions' contribution with amplitude larger than 0.1
             if (fabs(dum)  > 0.1) {
             outfile->Printf("%4d %20.12lf\t", p, dum);
             if (p<o*v) {
                outfile->Printf("electron    excitation");
                size_t a =p%v;
                size_t i = (p-a)/v;
                outfile->Printf("%4d   ->%4d\n",i+1,a+o+1);
             }
             else if ((o*v-1< p) && (p < 2*o*v)) {
                outfile->Printf("electron de-excitation");
                size_t a =p%v;
                size_t i = (p-a)/v;
                outfile->Printf("%4d   ->%4d\n",i+1,a+o+1);
             }
             else if (p==2*o*v) {outfile->Printf("photon excitation\n");
             }
             else {outfile->Printf("photon de-excitation\n");
             }
             }
        }
*/
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
    // TODO
    //build_gm(N, L, m, gm);
    // TODO
    //build_gm(N, L, n, gn);
    // TODO
    //build_sigma_m(N, L, x, y, m, sigma_m_h, sigma_m_s); // sigma_m_s need to be overwritten

    // x, y
    for (int i = 0; i < oa*va; i++) {
        for (int I = 0; I < L; I++) {
            sigmah[i        ][I] =  Ax[I*N + i] + By[I*N + i]; // TODO + gm[I*N+i] + gn[I*N+i];
            sigmah[i + oa*va][I] =  Bx[I*N + i] + Ay[I*N + i]; // TODO + gm[I*N+i] + gn[I*N+i];

            sigmas[i        ][I] =  x[I*N + i];
            sigmas[i + oa*va][I] = -y[I*N + i];
        }
    }
    for (int i = 0; i < ob*vb; i++) {
        for (int I = 0; I < L; I++) {
            sigmah[2*oa*va + i        ][I] =  Ax[I*N + oa*va + i] + By[I*N + oa*va + i]; // TODO + gm[I*N+i] + gn[I*N+i];
            sigmah[2*oa*va + i + ob*vb][I] =  Bx[I*N + oa*va + i] + Ay[I*N + oa*va + i]; // TODO + gm[I*N+i] + gn[I*N+i];

            sigmas[2*oa*va + i        ][I] =  x[I*N + oa*va + i];
            sigmas[2*oa*va + i + ob*vb][I] = -y[I*N + oa*va + i];
        }
    }
    // m, n
    // TODO:
/*
    for (int I = 0; I < L; I++) {
        sigmah[2*o*v  ][I] =  sigma_m_h[I] + cavity_frequency_[2] * m[I];
        sigmas[2*o*v  ][I] =  m[I]; 

        sigmah[2*o*v+1][I] =  sigma_m_h[I] + cavity_frequency_[2] * n[I]; 
        sigmas[2*o*v+1][I] = -n[I];
    }
*/

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

    // TODO: need separate dipole_a and dipole_b
    double ** dx = dipole_[0]->pointer();
    double ** dy = dipole_[1]->pointer();
    double ** dz = dipole_[2]->pointer();

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

// TODO:
/*
    double * tmp_a = (double*)malloc(o*v*sizeof(double));
    double * tmp_b = (double*)malloc(o*o*sizeof(double));
*/

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
// TODO
/*
        // J-like contribution from dipole self energy
        double dipole_Ja = 0.0;
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                int ia = i * v + a;

                dipole_Ja += c[ia] * dz[i][a+o];
            }
        }

        // intermediate for K-like contribution from dipole self energy

        // ignore if following QED-TDDFT outlined in J. Chem. Phys. 155, 064107 (2021)
        memset((void*)tmp_a,'\0',o*v*sizeof(double));
        memset((void*)tmp_b,'\0',o*o*sizeof(double));

        if ( options_.get_bool("QED_USE_RELAXED_ORBITALS") ) {

            // A term
            for (int i = 0; i < o; i++) {
                for (int a = 0; a < v; a++) {
                    double dum_a = 0.0;
                    for (int j = 0; j < o; j++) {
                        int ja = j * v + a;

                        dum_a += c[ja] * dz[i][j];
                    }

                    tmp_a[a*o+i] = dum_a;
                }
            }

            // B term
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    double dum_b = 0.0;
                    for (int b = 0; b < v; b++) {
                        int jb = j * v + b;

                        dum_b += c[jb] * dz[b+o][i];
                    }

                    tmp_b[i*o+j] = dum_b;
                }
            }
        }

        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {

                double dipole_Ja_ia = dipole_Ja * dz[i][a+o];

                // A term
                double dipole_Ka_A = 0.0;
                for (int b = 0; b < v; b++) {
                    dipole_Ka_A += tmp_a[b*o+i] * dz[a+o][b+o];
                }
                // B term
                double dipole_Ka_B = 0.0;
                for (int j = 0; j < o; j++) {
                    dipole_Ka_B += tmp_b[i*o+j] * dz[a+o][j];
                }

                int ia = i * v + a;

                Au[I*N+ia] += lambda_z * lambda_z * (dipole_Ja_ia + dipole_Ja_ia - dipole_Ka_A);
                Bu[I*N+ia] += lambda_z * lambda_z * (dipole_Ja_ia + dipole_Ja_ia - dipole_Ka_B);
            }
        }
*/

    }
// TODO:
/*
    free(tmp_a);
    free(tmp_b);
*/
    free(c);
}

void PolaritonicUTDDFT::build_gm(int N, int L, double *m, double *gm) {

    if ( n_photon_states_ > 2 ) {
        throw PsiException("qed-utddft only works for n_photon_states <= 2",__FILE__,__LINE__);
    }

    int o = nalpha_;
    int v = nso_ - nalpha_;

    double coupling_factor_x = cavity_frequency_[0] * cavity_coupling_strength_[0];
    double coupling_factor_y = cavity_frequency_[1] * cavity_coupling_strength_[1];
    double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];

    double lambda_x = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
    double lambda_y = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

    double ** dx = dipole_[0]->pointer();
    double ** dy = dipole_[1]->pointer();
    double ** dz = dipole_[2]->pointer();

    for (int I = 0; I < L; I++) {

        // couple |0,1> to |ia,0>
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                int ia = i * v + a;

                // <ia| H |0,1>
                gm[I*N+ia] = -sqrt(2.0) * coupling_factor_z * dz[i][a+o] * m[I];
            }
        }
    }
}

void PolaritonicUTDDFT::build_sigma_m(int N, int L, double *x, double *y, double *m, double *sigma_m_r, double *sigma_m_l) {

    if ( n_photon_states_ > 2 ) {
        throw PsiException("qed-utddft only works for n_photon_states <= 2",__FILE__,__LINE__);
    }

    int o = nalpha_;
    int v = nso_ - nalpha_;

    double coupling_factor_x = cavity_frequency_[0] * cavity_coupling_strength_[0];
    double coupling_factor_y = cavity_frequency_[1] * cavity_coupling_strength_[1];
    double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];

    double lambda_x = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
    double lambda_y = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

    double ** dx = dipole_[0]->pointer();
    double ** dy = dipole_[1]->pointer();
    double ** dz = dipole_[2]->pointer();

    for (int I = 0; I < L; I++) {

        // |0,1> diagonal
        sigma_m_r[I] = 0.0;//cavity_frequency_[2] * m[I];
        sigma_m_l[I] = 0.0;//cavity_frequency_[2] * m[I];

        // couple |0,1> to |ia,0>

        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {

                int ia = i * v + a;

                // <ia| H |0,1>
                double factor = -sqrt(2.0)*coupling_factor_z * dz[i][a+o];
                    
                sigma_m_r[I] += factor * ( x[I*N+ia] + y[I*N+ia] );
                sigma_m_l[I] += factor * ( x[I*N+ia] - y[I*N+ia] );
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

