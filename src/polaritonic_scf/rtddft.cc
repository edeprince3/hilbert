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

#include "rtddft.h"

#include <misc/blas.h>
#include <misc/hilbert_psifiles.h>
#include <misc/omp.h>
#include <misc/threeindexintegrals.h>
#include <misc/nonsym_davidson_solver.h>

using namespace psi;
using namespace fnocc;

namespace hilbert{ 

PolaritonicRTDDFT::PolaritonicRTDDFT(std::shared_ptr<Wavefunction> reference_wavefunction, Options& options_, std::shared_ptr<Wavefunction> dummy_wfn):
    PolaritonicHF(reference_wavefunction,options_) {
    common_init(dummy_wfn);
}

PolaritonicRTDDFT::~PolaritonicRTDDFT() {

/*
    free(int1_);
    free(int2_);
*/

}

void PolaritonicRTDDFT::common_init(std::shared_ptr<Wavefunction> dummy_wfn) {

    outfile->Printf("\n");
    outfile->Printf("\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    Polaritonic Restricted TDDFT                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf("\n");

    // RTTDFT only works with TDA for now
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

    // ensure closed shell
    if ( nalpha_ != nbeta_ ) {
        throw PsiException("polaritonic TDDFT only works with nalpha = nbeta (for now)",__FILE__,__LINE__);
    }

    // get primary basis:
    std::shared_ptr<BasisSet> primary = reference_wavefunction_->get_basisset("ORBITAL");

    int o = nalpha_;
    int v = nso_ - nalpha_;

/*
    nQ_ = 0;
    if ( options_.get_str("SCF_TYPE") == "DF" ) {

        // get auxiliary basis:
        std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");

        // total number of auxiliary basis functions
        nQ_ = auxiliary->nbf();

        std::shared_ptr<DFTensor> DF (new DFTensor(primary,auxiliary,Ca_,o,v,o,v,options_));

        std::shared_ptr<Matrix> tmpoo = DF->Qoo();
        std::shared_ptr<Matrix> tmpov = DF->Qov();
        std::shared_ptr<Matrix> tmpvv = DF->Qvv();

        double ** Qoo = tmpoo->pointer();
        double ** Qov = tmpov->pointer();
        double ** Qvv = tmpvv->pointer();

        int1_ = (double*)malloc(o*o*v*v*sizeof(double));
        int2_ = (double*)malloc(o*o*v*v*sizeof(double));

        memset((void*)int1_,'\0',o*o*v*v*sizeof(double));
        memset((void*)int2_,'\0',o*o*v*v*sizeof(double));

        F_DGEMM('n','t',o*v,o*v,nQ_,1.0,&(Qov[0][0]),o*v,&(Qov[0][0]),o*v,0.0,int1_,o*v);
        F_DGEMM('n','t',v*v,o*o,nQ_,1.0,&(Qvv[0][0]),v*v,&(Qoo[0][0]),o*o,0.0,int2_,v*v);

    }else if ( options_.get_str("SCF_TYPE") == "CD" ) {

        //outfile->Printf("    ==> Transform three-index integrals <==\n");
        //outfile->Printf("\n");

        double start = omp_get_wtime();
        ThreeIndexIntegrals(reference_wavefunction_,nQ_,memory_);

        double * Qmo = (double*)malloc(nmo_*(nmo_+1)/2*nQ_*sizeof(double));
        memset((void*)Qmo,'\0',nmo_*(nmo_+1)/2*nQ_*sizeof(double));

        std::shared_ptr<PSIO> psio(new PSIO());
        psio->open(PSIF_DCC_QMO,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_QMO,"(Q|mn) Integrals",(char*)Qmo,sizeof(double)*nQ_ * nmo_*(nmo_+1)/2);
        psio->close(PSIF_DCC_QMO,1);

        double end = omp_get_wtime();

        double * Qoo = (double*)malloc(o*o*nQ_*sizeof(double));
        double * Qov = (double*)malloc(o*v*nQ_*sizeof(double));
        double * Qvv = (double*)malloc(v*v*nQ_*sizeof(double));

        for (size_t Q = 0; Q < nQ_; Q++) {
            for (size_t i = 0; i < o; i++) {
                for (size_t j = 0; j < o; j++) {
                    Qoo[Q*o*o+i*o+j] = Qmo[Q*nmo_*(nmo_+1)/2+INDEX(i,j)];
                }
            }
            for (size_t i = 0; i < o; i++) {
                for (size_t a = 0; a < v; a++) {
                    Qov[Q*o*v+i*v+a] = Qmo[Q*nmo_*(nmo_+1)/2+INDEX(i,a+o)];
                }
            }
            for (size_t a = 0; a < v; a++) {
                for (size_t b = 0; b < v; b++) {
                    Qvv[Q*v*v+a*v+b] = Qmo[Q*nmo_*(nmo_+1)/2+INDEX(a+o,b+o)];
                }
            }
        }

        free(Qmo);

        //outfile->Printf("\n");
        //outfile->Printf("        Time for integral transformation:  %7.2lf s\n",end-start);
        //outfile->Printf("\n");

        int1_ = (double*)malloc(o*o*v*v*sizeof(double));
        int2_ = (double*)malloc(o*o*v*v*sizeof(double));

        memset((void*)int1_,'\0',o*o*v*v*sizeof(double));
        memset((void*)int2_,'\0',o*o*v*v*sizeof(double));

        F_DGEMM('n','t',o*v,o*v,nQ_,1.0,Qov,o*v,Qov,o*v,0.0,int1_,o*v);
        F_DGEMM('n','t',v*v,o*o,nQ_,1.0,Qvv,v*v,Qoo,o*o,0.0,int2_,v*v);

        free(Qoo);
        free(Qov);
        free(Qvv);

    }
*/

    // determine the DFT functional and initialize the potential object
    psi::scf::HF* scfwfn = (psi::scf::HF*)dummy_wfn.get();
    std::shared_ptr<SuperFunctional> functional = scfwfn->functional();
    //std::shared_ptr<VBase> potential = VBase::build_V(primary,functional,options_,(options_.get_str("REFERENCE") == "RKS" ? "RV" : "UV"));
    potential_ = (std::shared_ptr<VBase>)VBase::build_V(primary,functional,options_,(options_.get_str("REFERENCE") == "RKS" ? "RV" : "UV"));

    // initialize potential object
    potential_->initialize();

    // apparently compute_Vx wants me to set the density
    potential_->set_D({Da_});

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

        std::shared_ptr<DiskDFJK> myjk = (std::shared_ptr<DiskDFJK>)(new DiskDFJK(primary,auxiliary));

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

        std::shared_ptr<CDJK> myjk = (std::shared_ptr<CDJK>)(new CDJK(primary,options_.get_double("CHOLESKY_TOLERANCE")));

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



}

double PolaritonicRTDDFT::compute_energy() {

    outfile->Printf("\n");
    outfile->Printf("    No. basis functions:            %5i\n",nso_);
    outfile->Printf("    No. electrons:                  %5i\n",nalpha_ + nbeta_);

    if ( n_photon_states_ > 1 ) {
        update_cavity_terms();
    }

    // transform dipole integrals to MO basis
    dipole_[0]->transform(Ca_);
    dipole_[1]->transform(Ca_);
    dipole_[2]->transform(Ca_);

    int o = nalpha_;
    int v = nso_ - nalpha_;

    double * ea = epsilon_a_->pointer();

    double coupling_factor_x = cavity_frequency_[0] * cavity_coupling_strength_[0];
    double coupling_factor_y = cavity_frequency_[1] * cavity_coupling_strength_[1];
    double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];

    double lambda_x = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
    double lambda_y = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

    double ** dx = dipole_[0]->pointer();
    double ** dy = dipole_[1]->pointer();
    double ** dz = dipole_[2]->pointer();

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
    int N = 2*o*v + 2;

    // number of desired roots
    int M = options_.get_int("NUM_ROOTS");

    // maximum subspace size
    size_t MD = options_.get_int("MAXDIM");
    size_t ID = options_.get_int("INDIM");
    size_t maxdim = MD*M;
    size_t init_dim = ID*M;

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
    //david->solve(Hdiag,
    //             N,
    //             M,
    //             revalp,
    //             rerp,
    //             relp,
    //             [&](int N, int maxdim, int L, double **Q, double **sigmar, double **sigmal) {
    //                    build_sigma(N,maxdim,L, Q,sigmar,sigmal);
    //             }, maxdim, init_dim, residual_norm, use_residual_norm);
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
    outfile->Printf("    cQED-TDDFT energies:\n");
    outfile->Printf("\n");
    outfile->Printf("    ");
    outfile->Printf("%5s","state");
    outfile->Printf(" %5s","type");
    outfile->Printf(" %20s","energy (Eh)");
    outfile->Printf(" %20s","ex energy (Eh)");
    outfile->Printf(" %20s","photon weight");
    outfile->Printf(" %10s","mu_x");
    outfile->Printf(" %10s","mu_y");
    outfile->Printf(" %10s","mu_z");
    outfile->Printf(" %10s","f");
    outfile->Printf("\n");

    for (int state = 0; state < M; state++) {

        // normalization
        double m_weight = rerp[state][2*o*v  ] * rerp[state][2*o*v  ];
        double n_weight = rerp[state][2*o*v+1] * rerp[state][2*o*v+1];
        double photon_weight = m_weight - n_weight;
        double x_weight = 0.0;
        double y_weight = 0.0;
        for (int j = 0; j < o*v; j++) {
            x_weight += rerp[state][j    ] * rerp[state][j    ];
            y_weight += rerp[state][j+o*v] * rerp[state][j+o*v];
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
        for (int j = 0; j < o; j++) {
            for (int b = 0; b < v; b++) {
                double cr = (rerp[state][j*v+b] + rerp[state][j*v+b+o*v])/nrm;
                mu_x_r += dx[j][b+o] * cr;
                mu_y_r += dy[j][b+o] * cr;
                mu_z_r += dz[j][b+o] * cr;
            }
        }
        double w = revalp[state];//sqrt(revalp[state]);
        double f = 2.0 * 2.0/3.0*w*(mu_x_r * mu_x_r + mu_y_r * mu_y_r + mu_z_r * mu_z_r);

        outfile->Printf("    %5i %5s %20.12lf %20.12lf %20.12lf %10.6lf %10.6lf %10.6lf %10.6lf\n", state, type.c_str(),w,energy_ + w,photon_weight,mu_x_r,mu_y_r,mu_z_r,f);
    }
    outfile->Printf("%d %d %d\n",o,v,N);
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
    }


    return 0.0;
}

// cQED-TDDFT:
// 
// |  A  B  g*  g* | ( X )     |  1  0  0  0 | ( X )
// |  B  A  g*  g* | ( Y ) = W |  0 -1  0  0 | ( Y )
// |  g  g  w   0  | ( M )     |  0  0  1  0 | ( M )
// |  g  g  0   w  | ( N )     |  0  0  0 -1 | ( N )
// 
// 
// 
// this function solves the generalized eigenvalue problem S c = 1/w H c
void PolaritonicRTDDFT::build_sigma_generalized(int N, int maxdim, int L, double **Q, double **sigmah, double **sigmas){

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

    int o = nalpha_;
    int v = nso_ - nalpha_;

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

    // unpack input vectors
    for (int I = 0; I < L; I++) {
        for (int j = 0; j < o*v; j++) {
            x[I*N+j] = Q[I][j];
            y[I*N+j] = Q[I][j + o*v];
        }
        m[I] = Q[I][2*o*v];
        n[I] = Q[I][2*o*v+1];
    }

    build_Au_Bu(N, L, x, Ax, Bx);
    build_Au_Bu(N, L, y, Ay, By);
    build_gm(N, L, m, gm);
    build_gm(N, L, n, gn);
    build_sigma_m(N, L, x, y, m, sigma_m_h, sigma_m_s); // sigma_m_s need to be overwritten

    // x, y
    for (int i = 0; i < o*v; i++) {
        for (int I = 0; I < L; I++) {
            sigmah[i      ][I] =  Ax[I*N+i] + By[I*N+i] + gm[I*N+i] + gn[I*N+i];
            sigmah[i + o*v][I] =  Bx[I*N+i] + Ay[I*N+i] + gm[I*N+i] + gn[I*N+i];

            sigmas[i      ][I] =  x[I*N+i];
            sigmas[i + o*v][I] = -y[I*N+i];
        }
    }
    // m, n
    for (int I = 0; I < L; I++) {
        sigmah[2*o*v  ][I] =  sigma_m_h[I] + cavity_frequency_[2] * m[I];
        sigmas[2*o*v  ][I] =  m[I]; 

        sigmah[2*o*v+1][I] =  sigma_m_h[I] + cavity_frequency_[2] * n[I]; 
        sigmas[2*o*v+1][I] = -n[I];
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

// cQED-TDDFT:
// 
// |  A  B  g*  g* | ( X )     ( X )
// | -B -A -g* -g* | ( Y ) = W ( Y )
// |  g  g  w   0  | ( M )     ( M )
// | -g -g  0  -w  | ( N )     ( N )
// 
// 
// 
// ( X Y M N ) |  A  B  g*  g* |     W ( X Y M N )
//             | -B -A -g* -g* |  =  
//             |  g  g  w   0  |     
//             | -g -g  0  -w  |     
//           
// this function actually solves the square of the RPA problem
void PolaritonicRTDDFT::build_sigma(int N, int maxdim, int L, double **Q, double **sigmar, double **sigmal){

    for (int i = 0; i < N; i++) {
        for (int I = 0; I < maxdim; I++) {
            sigmar[i][I] = 0.0;
            sigmal[i][I] = 0.0;
        }
    }

    // TDDFT-PF:

    // sigmar_x =  Ax + By + g*m + g*n
    // sigmar_y = -Bx - Ay - g*m - g*n
    // sigmar_m =  gx + gy +  wm
    // sigmar_n = -gx - gy       -  wn

    // sigmal_x =  Ax - By + g*m - g*n
    // sigmal_y =  Bx - Ay + g*m - g*n
    // sigmal_m =  gx - gy +  wm
    // sigmal_n =  gx - gy       -  wn

    int o = nalpha_;
    int v = nso_ - nalpha_;

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
    double *sigma_m_r = (double*)malloc(L*sizeof(double));
    double *sigma_m_l = (double*)malloc(L*sizeof(double));
    double *n = (double*)malloc(L*sizeof(double));
    double *gn = (double*)malloc(N*L*sizeof(double));
    //double *sigma_n_r = (double*)malloc(L*sizeof(double));
    //double *sigma_n_l = (double*)malloc(L*sizeof(double));

    memset((void*)m,'\0',L*sizeof(double));
    memset((void*)gm,'\0',N*L*sizeof(double));
    memset((void*)sigma_m_r,'\0',L*sizeof(double));
    memset((void*)sigma_m_l,'\0',L*sizeof(double));
    memset((void*)n,'\0',L*sizeof(double));
    memset((void*)gn,'\0',N*L*sizeof(double));
    //memset((void*)sigma_n_r,'\0',L*sizeof(double));
    //memset((void*)sigma_n_l,'\0',L*sizeof(double));

    // unpack input vectors
    for (int I = 0; I < L; I++) {
        for (int j = 0; j < o*v; j++) {
            x[I*N+j] = Q[I][j];
            y[I*N+j] = Q[I][j + o*v];
        }
        m[I] = Q[I][2*o*v];
        n[I] = Q[I][2*o*v+1];
    }

    build_Au_Bu(N, L, x, Ax, Bx);
    build_Au_Bu(N, L, y, Ay, By);
    build_gm(N, L, m, gm);
    build_gm(N, L, n, gn);
    build_sigma_m(N, L, x, y, m, sigma_m_r, sigma_m_l);
    //build_sigma_m(N, L, x, y, n, sigma_n_r, sigma_n_l);

    // sigmar_x =  Ax + By + g*m + g*n
    // sigmar_y = -Bx - Ay - g*m - g*n
    // sigmar_m =  gx + gy +  wm
    // sigmar_n = -gx - gy       -  wn

    // sigmal_x =  Ax - By + g*m - g*n
    // sigmal_y =  Bx - Ay + g*m - g*n
    // sigmal_m =  gx - gy +  wm
    // sigmal_n =  gx - gy       -  wn

    // x, y
    for (int i = 0; i < o*v; i++) {
        for (int I = 0; I < L; I++) {
            sigmar[i      ][I] =  Ax[I*N+i] + By[I*N+i] + gm[I*N+i] + gn[I*N+i];
            sigmar[i + o*v][I] = -Bx[I*N+i] - Ay[I*N+i] - gm[I*N+i] - gn[I*N+i];

            sigmal[i      ][I] =  Ax[I*N+i] - By[I*N+i] + gm[I*N+i] - gn[I*N+i];
            sigmal[i + o*v][I] =  Bx[I*N+i] - Ay[I*N+i] + gm[I*N+i] - gn[I*N+i];
        }
    }
    // m, n
    for (int I = 0; I < L; I++) {
        sigmar[2*o*v  ][I] =  sigma_m_r[I] + cavity_frequency_[2] * m[I];
        sigmal[2*o*v  ][I] =  sigma_m_l[I] + cavity_frequency_[2] * m[I];
        sigmar[2*o*v+1][I] = -sigma_m_r[I] - cavity_frequency_[2] * n[I]; // note sign
        sigmal[2*o*v+1][I] =  sigma_m_l[I] - cavity_frequency_[2] * n[I];
    }

    // ok ... now square of the problem

    // unpack input vectors (right)
    for (int I = 0; I < L; I++) {
        for (int j = 0; j < o*v; j++) {
            x[I*N+j] = sigmar[j][I];
            y[I*N+j] = sigmar[j + o*v][I];
        }
        m[I] = sigmar[2*o*v  ][I];
        n[I] = sigmar[2*o*v+1][I];
    }

    build_Au_Bu(N, L, x, Ax, Bx);
    build_Au_Bu(N, L, y, Ay, By);
    build_gm(N, L, m, gm);
    build_gm(N, L, n, gn);
    build_sigma_m(N, L, x, y, m, sigma_m_r, sigma_m_l);
    //build_sigma_m(N, L, x, y, n, sigma_n_r, sigma_n_l);

    // x, y
    for (int i = 0; i < o*v; i++) {
        for (int I = 0; I < L; I++) {
            sigmar[i      ][I] =  Ax[I*N+i] + By[I*N+i] + gm[I*N+i] + gn[I*N+i];
            sigmar[i + o*v][I] = -Bx[I*N+i] - Ay[I*N+i] - gm[I*N+i] - gn[I*N+i];
        }
    }
    // m, n
    for (int I = 0; I < L; I++) {
        sigmar[2*o*v  ][I] =  sigma_m_r[I] + cavity_frequency_[2] * m[I];
        sigmar[2*o*v+1][I] = -sigma_m_r[I] - cavity_frequency_[2] * n[I]; // note sign
    }

    // unpack input vectors (left)
    for (int I = 0; I < L; I++) {
        for (int j = 0; j < o*v; j++) {
            x[I*N+j] = sigmal[j      ][I];
            y[I*N+j] = sigmal[j + o*v][I];
        }
        m[I] = sigmal[2*o*v  ][I];
        n[I] = sigmal[2*o*v+1][I];
    }
    
    build_Au_Bu(N, L, x, Ax, Bx);
    build_Au_Bu(N, L, y, Ay, By);
    build_gm(N, L, m, gm);
    build_gm(N, L, n, gn);
    build_sigma_m(N, L, x, y, m, sigma_m_r, sigma_m_l);
    //build_sigma_m(N, L, x, y, n, sigma_n_r, sigma_n_l);

    // x, y
    for (int i = 0; i < o*v; i++) {
        for (int I = 0; I < L; I++) {
            sigmal[i      ][I] =  Ax[I*N+i] - By[I*N+i] + gm[I*N+i] - gn[I*N+i];
            sigmal[i + o*v][I] =  Bx[I*N+i] - Ay[I*N+i] + gm[I*N+i] - gn[I*N+i];
        }
    }
    // m, n
    for (int I = 0; I < L; I++) {
        sigmal[2*o*v  ][I] =  sigma_m_l[I] + cavity_frequency_[2] * m[I];
        sigmal[2*o*v+1][I] =  sigma_m_l[I] - cavity_frequency_[2] * n[I];
    }

    free(x);
    free(y);
    free(Ax);
    free(Bx);
    free(Ay);
    free(By);
    free(m);
    free(sigma_m_r);
    free(sigma_m_l);
    free(gm);
    free(n);
    //free(sigma_n_r);
    //free(sigma_n_l);
    free(gn);

}

/*
// square of TDDFT problem
void PolaritonicRTDDFT::build_sigma(int N, int maxdim, int L, double **Q, double **sigmar, double **sigmal){

    for (int i = 0; i < N; i++) {
        for (int I = 0; I < maxdim; I++) {
            sigmar[i][I] = 0.0;
            sigmal[i][I] = 0.0;
        }
    }

    // TDA: 
    // sigmar = Au
    // sigmal = Au

    // RPA: 
    // sigmar = AAu + ABu - BAu - BBu
    // sigmal = AAu - ABu + BAu - BBu

    double * Au = (double*)malloc(N*L*sizeof(double));
    double * Bu = (double*)malloc(N*L*sizeof(double));
    double * u  = (double*)malloc(N*L*sizeof(double));

    memset((void*)Au,'\0',N*L*sizeof(double));
    memset((void*)Bu,'\0',N*L*sizeof(double));
    memset((void*)u, '\0',N*L*sizeof(double));

    // unpack input vectors
    for (int I = 0; I < L; I++) {
        for (int j = 0; j < N; j++) {
            u[I*N+j] = Q[I][j];
        }
    }

    double *AAu, *ABu, *BAu, *BBu;
    if ( !options_.get_bool("TDSCF_TDA") ) {
        AAu = (double*)malloc(N*L*sizeof(double));
        ABu = (double*)malloc(N*L*sizeof(double));
        BAu = (double*)malloc(N*L*sizeof(double));
        BBu = (double*)malloc(N*L*sizeof(double));
        memset((void*)AAu,'\0',N*L*sizeof(double));
        memset((void*)ABu,'\0',N*L*sizeof(double));
        memset((void*)BAu,'\0',N*L*sizeof(double));
        memset((void*)BBu,'\0',N*L*sizeof(double));
    }

    // TDA or RPA?
    if ( options_.get_bool("TDSCF_TDA") ) {

        build_Au_Bu(N, L, u, Au, Bu);

        for (int i = 0; i < N; i++) {
            for (int I = 0; I < L; I++) {
                sigmar[i][I] = Au[I*N+i];
                sigmal[i][I] = Au[I*N+i];
            }
        }

    }else {

        build_Au_Bu(N, L,  u,  Au,  Bu);
        build_Au_Bu(N, L, Au, AAu, BAu);
        build_Au_Bu(N, L, Bu, ABu, BBu);

        for (int i = 0; i < N; i++) {
            for (int I = 0; I < L; I++) {
                sigmar[i][I] = AAu[I*N+i] + ABu[I*N+i] - BAu[I*N+i] - BBu[I*N+i] ;
                sigmal[i][I] = AAu[I*N+i] - ABu[I*N+i] + BAu[I*N+i] - BBu[I*N+i] ;
            }
        }

    }

    free(u);
    free(Au);
    free(Bu);
    if ( !options_.get_bool("TDSCF_TDA") ) {
        free(AAu);
        free(ABu);
        free(BAu);
        free(BBu);
    }

}
*/

// build Au and Bu according to Yihan's formalism that includes only (1) single 
// excitations and de-excitations and (2) photon excitations and de-excitations, 
// but no coupled photon / electron excitations. Also, this formalism only works
// for n_photon_states_ = 2

// Yihan calls this TDDFT-JC
// note that Au will contain contribution from m

void PolaritonicRTDDFT::build_Au_Bu(int N, int L, double *u, double *Au, double *Bu){

    if ( n_photon_states_ > 2 ) {
        throw PsiException("plugin cqed-tddft only works for n_photon_states <= 2",__FILE__,__LINE__);
    }

    int o = nalpha_;
    int v = nso_ - nalpha_;

    double * c = (double*)malloc(N*sizeof(double));

    double ** cap = Ca_->pointer();

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

        double ** clp = cla->pointer();
        double ** crp = cra->pointer();

        cra->zero();
        cla->zero();

        for (int mu = 0; mu < nso_; mu++) {
            for (int i = 0; i < o; i++) {

                // left is plain orbitals
                clp[mu][i] = cap[mu][i];

                // right is modified orbitals
                double dum = 0.0;
                for (int a = 0; a < v; a++) {
                    int ia = i * v + a;
                    dum += cap[mu][a+o] * c[ia];
                }
                crp[mu][i] = dum;

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

    double * tmp_a = (double*)malloc(o*v*sizeof(double));
    double * tmp_b = (double*)malloc(o*o*sizeof(double));

    for (int I = 0; I < L; I++) {

        // unpack a vector
        for (int j = 0; j < N; j++) {
            c[j] = u[I*N+j];
        }

        // a <- a coulomb
        std::shared_ptr<Matrix> sa (new Matrix(jk_->J()[count]));

        // b <- b coulomb
        sa->scale(2.0);

        // xc?
        // a <- a
        if ( needs_xc_ ) {
            sa->axpy(2.0,Vx[count]);
        }

        // exact exchange?
        // a <- a
        // b <- b
        if (is_x_hybrid_) {
            sa->axpy(-x_alpha_,jk_->K()[count]);
        }

        // LRC functional?
        // a <- a
        // b <- b
        if (is_x_lrc_) {
            double beta = 1.0 - x_alpha_;
            sa->axpy(-beta,jk_->wK()[count]);
        }

        // update counter
        count += 1;

        // transform jk, e.g., j(a,i) = j(mu,nu) c(mu,a) c(nu,i)

        sa->transform(Ca_);

        double ** sap = sa->pointer();

        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                int ia = i * v + a;

                Au[I*N+ia] = sap[i][a+o];
                Bu[I*N+ia] = sap[a+o][i];
            }
        }

        // singles diagonals
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                int ia = i * v + a;
                Au[I*N+ia] += c[ia] * (ea[a+o] - ea[i]);
            }
        }

        // dipole self energy (for singles)

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

    }

    free(tmp_a);
    free(tmp_b);
    free(c);

}

void PolaritonicRTDDFT::build_gm(int N, int L, double *m, double *gm) {

    if ( n_photon_states_ > 2 ) {
        throw PsiException("plugin cqed-tddft only works for n_photon_states <= 2",__FILE__,__LINE__);
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

void PolaritonicRTDDFT::build_sigma_m(int N, int L, double *x, double *y, double *m, double *sigma_m_r, double *sigma_m_l) {

    if ( n_photon_states_ > 2 ) {
        throw PsiException("plugin cqed-tddft only works for n_photon_states <= 2",__FILE__,__LINE__);
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

// build Au and Bu, where u includes |0,0> and no de-excitation |0,-n>
/*
void PolaritonicRTDDFT::build_Au_Bu(int N, int L, double *u, double *Au, double *Bu){

    int o = nalpha_;
    int v = nso_ - nalpha_;

    double * c = (double*)malloc(N*sizeof(double));

    double ** cap = Ca_->pointer();

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

        // loop over photon states
        for (int n = 0; n < n_photon_states_; n++) {

            // alpha
            std::shared_ptr<Matrix> cra (new Matrix(Ca_) );
            std::shared_ptr<Matrix> cla (new Matrix(Ca_) );

            double ** clp = cla->pointer();
            double ** crp = cra->pointer();

            cra->zero();
            cla->zero();

            for (int mu = 0; mu < nso_; mu++) {
                for (int i = 0; i < o; i++) {

                    // left is plain orbitals
                    clp[mu][i] = cap[mu][i];

                    // right is modified orbitals
                    double dum = 0.0;
                    for (int a = 0; a < v; a++) {
                        int ian_a = (i * v + a      ) * n_photon_states_ + n;
                        dum += cap[mu][a+o] * c[ian_a];
                    }
                    crp[mu][i] = dum;

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

    double * tmp_a = (double*)malloc(o*v*sizeof(double));

    for (int I = 0; I < L; I++) {

        // unpack a vector
        for (int j = 0; j < N; j++) {
            c[j] = u[I*N+j];
        }

        for (int n = 0; n < n_photon_states_; n++) {

            // a <- a coulomb
            std::shared_ptr<Matrix> sa (new Matrix(jk_->J()[count]));

            // b <- b coulomb
            sa->scale(2.0);

            // xc?
            // a <- a
            if ( needs_xc_ ) {
                sa->axpy(2.0,Vx[count]);
                //sb->axpy(1.0,Vx[count]);

                //sa->axpy(1.0,Vx[count+1]);
                //sb->axpy(1.0,Vx[count+1]);
            }

            // exact exchange?
            // a <- a
            // b <- b
            if (is_x_hybrid_) {
                sa->axpy(-x_alpha_,jk_->K()[count]);
            }

            // LRC functional?
            // a <- a
            // b <- b
            if (is_x_lrc_) {
                double beta = 1.0 - x_alpha_;
                sa->axpy(-beta,jk_->wK()[count]);
            }

            // update counter
            count += 1;

            // transform jk, e.g., j(a,i) = j(mu,nu) c(mu,a) c(nu,i)

            sa->transform(Ca_);

            double ** sap = sa->pointer();

            for (int i = 0; i < o; i++) {
                for (int a = 0; a < v; a++) {
                    int ian_a = (i * v + a      ) * n_photon_states_ + n;

                    Au[I*N+ian_a] = sap[i][a+o];
                    Bu[I*N+ian_a] = sap[a+o][i];
                }
            }
        }

        // singles diagonals
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                for (int n = 0; n < n_photon_states_; n++) {
                    int ian_a = (i * v + a      ) * n_photon_states_ + n;

                    Au[I*N+ian_a] += c[ian_a] * (ea[a+o] - ea[i] + n * cavity_frequency_[2]);
                }
            }
        }

        // now, |0,n> diagonals
        int off = o * v * n_photon_states_;
        for (int n = 0; n < n_photon_states_; n++) {
            Au[I*N+off+n] += n * cavity_frequency_[2] * c[off+n];
        }

        // dipole self energy
        for (int n = 0; n < n_photon_states_; n++) {

            // J-like contribution from dipole self energy
            double dipole_Ja = 0.0;
            for (int i = 0; i < o; i++) {
                for (int a = 0; a < v; a++) {
                    int ian_a = (i * v + a      ) * n_photon_states_ + n;

                    dipole_Ja += c[ian_a] * dz[i][a+o];
                }
            }

            // intermediate for K-like contribution from dipole self energy
            memset((void*)tmp_a,'\0',o*v*sizeof(double));
            for (int i = 0; i < o; i++) {
                for (int a = 0; a < v; a++) {
                    double dum_a = 0.0;
                    double dum_b = 0.0;
                    for (int j = 0; j < o; j++) {
                        int jan_a = (j * v + a      ) * n_photon_states_ + n;

                        dum_a += c[jan_a] * dz[i][j];
                    }

                    tmp_a[a*o+i] = dum_a;
                }
            }

            for (int i = 0; i < o; i++) {
                for (int a = 0; a < v; a++) {

                    double dipole_Ja_ia = dipole_Ja * dz[i][a+o];

                    double dipole_Ka = 0.0;
                    for (int b = 0; b < v; b++) {
                        dipole_Ka += tmp_a[b*o+i] * dz[a+o][b+o];
                    }

                    int ian_a = (i * v + a      ) * n_photon_states_ + n;

                    Au[I*N+ian_a] += lambda_z * lambda_z * (dipole_Ja_ia + dipole_Ja_ia - dipole_Ka);
                    Bu[I*N+ian_a] += lambda_z * lambda_z * (dipole_Ja_ia + dipole_Ja_ia - dipole_Ka);
                }
            }
        }

        // coupling terms
        for (int n = 0; n < n_photon_states_; n++) {

            double sqrt_n = sqrt(n);

            for (int i = 0; i < o; i++) {
                for (int a = 0; a < v; a++) {

                    int ian_a = (i * v + a      ) * n_photon_states_ + n;

                    double dum_a = 0.0;

                    // d(n,m+1)
                    int m = n - 1;
                    if ( m > -1 ) {
                        // d(i,j)
                        for (int b = 0; b < v; b++) {
                            int ibm_a = (i * v + b      ) * n_photon_states_ + m;
                            dum_a -= sqrt_n * c[ibm_a] * dz[a+o][b+o];
                        }
                        // d(a,b)
                        for (int j = 0; j < o; j++) {
                            int jam_a = (j * v + a      ) * n_photon_states_ + m;
                            dum_a += sqrt_n * c[jam_a] * dz[i][j];
                        }
                        // d(i,j) d(a,b)
                        int iam_a = (i * v + a      ) * n_photon_states_ + m;
                        for (int k = 0; k < o; k++) {
                            dum_a -= sqrt_n * c[iam_a] * dz[k][k];
                        }
                        // more coherent-state terms ... these affect diagonals in the electronic basis
                        dum_a += sqrt_n * tot_dip_z_ * c[iam_a];
                    }

                    // d(n,m-1)
                    m = n + 1;
                    if ( m < n_photon_states_ ) {
                        double sqrt_m = sqrt(m);
                        // d(i,j)
                        for (int b = 0; b < v; b++) {
                            int ibm_a = (i * v + b      ) * n_photon_states_ + m;
                            dum_a -= sqrt_m * c[ibm_a] * dz[a+o][b+o];
                        }
                        // d(a,b)
                        for (int j = 0; j < o; j++) {
                            int jam_a = (j * v + a      ) * n_photon_states_ + m;
                            dum_a += sqrt_m * c[jam_a] * dz[i][j];
                        }
                        // d(i,j) d(a,b)
                        int iam_a = (i * v + a      ) * n_photon_states_ + m;
                        int iam_b = (i * v + a + o*v) * n_photon_states_ + m;
                        for (int k = 0; k < o; k++) {
                            dum_a -= sqrt_m * c[iam_a] * dz[k][k];
                        }
                        // more coherent-state terms ... these affect diagonals in the electronic basis
                        dum_a += sqrt_m * tot_dip_z_ * c[iam_a];
                    }

                    Au[I*N+ian_a] += coupling_factor_z * dum_a;

                }
            }

        }

        // lastly, couple |0,n+1> and |0,n-1> to |ia,n>
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                for (int n = 0; n < n_photon_states_; n++) {
                    int ian_a = (i * v + a      ) * n_photon_states_ + n;
                    // case 1
                    if ( n < n_photon_states_ - 1 ) {
                        double factor = -sqrt(2.0)*coupling_factor_z * sqrt(n+1) * dz[i][a+o];

                        double val_a   = factor * c[ian_a];
                        double val_np1 = factor * c[off+n+1];

                        Au[I*N+ian_a]   += val_np1;
                        Au[I*N+off+n+1] += val_a;

                    }

                    // case 2
                    if ( n > 0 ) {
                        double factor = -sqrt(2.0)*coupling_factor_z * sqrt(n) * dz[i][a+o];

                        double val_a   = factor * c[ian_a];
                        double val_nm1 = factor * c[off+n-1];

                        Au[I*N+ian_a]   += val_nm1;
                        Au[I*N+off+n-1] += val_a;

                    }
                }
            }
        }

    }

    free(tmp_a);
    free(c);

}
*/

/*
// build approximation to diagonal of A matrix in basis that includes |0,0> and no de-excitation |0,-n>
double * PolaritonicRTDDFT::build_hamiltonian_diagonals(){

    int o = nalpha_;
    int v = nso_ - nalpha_;

    double * H = (double*)malloc((o*v+1)*n_photon_states_*sizeof(double));; 
    memset((void*)H,'\0',(o*v+1)*n_photon_states_*sizeof(double));

    double * ea = epsilon_a_->pointer();

    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            for (int n = 0; n < n_photon_states_; n++) {
                int ian_a = (i * v + a      ) * n_photon_states_ + n;
                H[ian_a] = ea[a+o] - ea[i] + n * cavity_frequency_[2];
            }
        }
    }

    // now, |0,n> diagonals
    int off = o * v * n_photon_states_;
    for (int n = 0; n < n_photon_states_; n++) {
        H[off+n] = n * cavity_frequency_[2];
    }

    // square diagonals for TDDFT
    if ( !options_.get_bool("TDSCF_TDA") ) {
        for (int i = 0; i < (o*v+1)*n_photon_states_; i++) {
            H[i] *= H[i];
        }
    }

    return H;
}
*/

// build approximation to diagonal of A matrix in basis that excludes |0,0> and de-excitation |0,-n>
double * PolaritonicRTDDFT::build_hamiltonian_diagonals(){

    int o = nalpha_;
    int v = nso_ - nalpha_;

    int dim = 2*o*v+2;

    double * H = (double*)malloc(dim*sizeof(double));
    memset((void*)H,'\0',dim*sizeof(double));

    double * ea = epsilon_a_->pointer();

    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            int ia = i * v + a;
            H[ia]     = ea[a+o] - ea[i];
            H[o*v+ia] = H[ia];
        }
    }

    // now, |0,1> diagonals
    if ( n_photon_states_ > 2 ) {
        throw PsiException("plugin cqed-tddft only works for n_photon_states <= 2",__FILE__,__LINE__);
    }
    int off = 2 * o * v;
    H[off  ] = cavity_frequency_[2];
    H[off+1] = cavity_frequency_[2];

    //for (int i = 0; i < dim; i++) {
    //    H[i] *= H[i];
    //}

    return H;
}
double * PolaritonicRTDDFT::build_overlap_diagonals(){

    int o = nalpha_;
    int v = nso_ - nalpha_;

    int dim = 2*o*v+2;

    double * S = (double*)malloc(dim*sizeof(double));
    memset((void*)S,'\0',dim*sizeof(double));

    double * ea = epsilon_a_->pointer();

    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            int ia = i * v + a;
            S[ia]     =  1.0;
            S[o*v+ia] = -1.0;
        }
    }

    // now, |0,1> diagonals
    if ( n_photon_states_ > 2 ) {
        throw PsiException("plugin cqed-tddft only works for n_photon_states <= 2",__FILE__,__LINE__);
    }
    int off = 2 * o * v;
    S[off  ] =  1.0;
    S[off+1] = -1.0;

    return S;
}



} // End namespaces

