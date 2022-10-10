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

    free(int1_);
    free(int2_);

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

    // UTTDFT only works with TDA for now
    if ( !options_.get_bool("TDSCF_TDA") ) {
        throw PsiException("polaritonic utddft only works with TDA df for now",__FILE__,__LINE__);
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

    // determine the DFT functional and initialize the potential object
    psi::scf::HF* scfwfn = (psi::scf::HF*)dummy_wfn.get();
    std::shared_ptr<SuperFunctional> functional = scfwfn->functional();
    //std::shared_ptr<VBase> potential = VBase::build_V(primary,functional,options_,(options_.get_str("REFERENCE") == "RKS" ? "RV" : "UV"));
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

double PolaritonicUTDDFT::compute_energy() {

    outfile->Printf("\n");
    outfile->Printf("    No. basis functions:            %5i\n",nso_);
    outfile->Printf("    No. auxiliary basis functions:  %5i\n",nQ_);
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

    std::shared_ptr<Matrix> ham (new Matrix((2*o*v+1)*n_photon_states_,(2*o*v+1)*n_photon_states_));

    double ** hp = ham->pointer();
    double * ea = epsilon_a_->pointer();
    double * eb = epsilon_b_->pointer();

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

    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            for (int n = 0; n < n_photon_states_; n++) {
                int ian_a = (i * v + a      ) * n_photon_states_ + n;
                int ian_b = (i * v + a + o*v) * n_photon_states_ + n;
                for (int j = 0; j < o; j++) {
                    for (int b = 0; b < v; b++) {
                        for (int m = 0; m < n_photon_states_; m++) {
                            int jbm_a = (j * v + b      ) * n_photon_states_ + m;
                            int jbm_b = (j * v + b + o*v) * n_photon_states_ + m;

                            // cases:
                            hp[ian_a][jbm_a] = 0.0;
                            hp[ian_a][jbm_b] = 0.0;
                            hp[ian_b][jbm_a] = 0.0;
                            hp[ian_b][jbm_b] = 0.0;

                            // full diagonal 
                            if ( i == j && a == b && m == n ) {
                                hp[ian_a][jbm_a] += ea[a+o] - ea[i] + HCavity_z->pointer()[m][m];
                                hp[ian_b][jbm_b] += eb[a+o] - eb[i] + HCavity_z->pointer()[m][m];
                            }

                            // diagonal in photon states
                            if ( m == n ) {
                                int iajb = i * o * v * v + a * o * v + j * v + b;
                                int ijab = i * o * v * v + j * v * v + a * v + b;
                                hp[ian_a][jbm_a] += int1_[iajb] - int2_[ijab];
                                hp[ian_a][jbm_b] += int1_[iajb];
                                hp[ian_b][jbm_a] += int1_[iajb];
                                hp[ian_b][jbm_b] += int1_[iajb] - int2_[ijab];

                                // dipole self energy contribution
                                hp[ian_a][jbm_a] += lambda_z * lambda_z * (dz[i][a+o]*dz[j][b+o] - dz[i][j]*dz[a+o][b+o]);
                                hp[ian_a][jbm_b] += lambda_z * lambda_z * (dz[i][a+o]*dz[j][b+o]);
                                hp[ian_b][jbm_a] += lambda_z * lambda_z * (dz[i][a+o]*dz[j][b+o]);
                                hp[ian_b][jbm_b] += lambda_z * lambda_z * (dz[i][a+o]*dz[j][b+o] - dz[i][j]*dz[a+o][b+o]);

                                // coherent-state basis terms ... actually, these should be contained in orbital energies
                                //if ( i == j ) {
                                //    hp[ian][jbm] -= lambda_z * lambda_z * dz[a+o][b+o] * tot_dip_z_;
                                //}
                                //if ( a == b ) {
                                //    hp[ian][jbm] -= lambda_z * lambda_z * dz[i][j] * tot_dip_z_;
                                //}

                            }

                            if ( i == j && n == m + 1 ) {
                                hp[ian_a][jbm_a] -= coupling_factor_z * sqrt(n) * dz[a+o][b+o];
                                hp[ian_b][jbm_b] -= coupling_factor_z * sqrt(n) * dz[a+o][b+o];
                            }

                            if ( i == j && n == m - 1 ) {
                                hp[ian_a][jbm_a] -= coupling_factor_z * sqrt(m) * dz[a+o][b+o];
                                hp[ian_b][jbm_b] -= coupling_factor_z * sqrt(m) * dz[a+o][b+o];
                            }

                            if ( a == b && n == m + 1 ) {
                                hp[ian_a][jbm_a] += coupling_factor_z * sqrt(n) * dz[i][j];
                                hp[ian_b][jbm_b] += coupling_factor_z * sqrt(n) * dz[i][j];
                            }

                            if ( a == b && n == m - 1 ) {
                                hp[ian_a][jbm_a] += coupling_factor_z * sqrt(m) * dz[i][j];
                                hp[ian_b][jbm_b] += coupling_factor_z * sqrt(m) * dz[i][j];
                            }

                            if ( a == b && i == j && n == m + 1 ) {
                                for (int k = 0; k < o; k++) {
                                    hp[ian_a][jbm_a] -= coupling_factor_z * sqrt(n) * dz[k][k];
                                    hp[ian_b][jbm_b] -= coupling_factor_z * sqrt(n) * dz[k][k];
                                }
                            }

                            if ( a == b && i == j && n == m - 1 ) {
                                for (int k = 0; k < o; k++) {
                                    hp[ian_a][jbm_a] -= coupling_factor_z * sqrt(m) * dz[k][k];
                                    hp[ian_b][jbm_b] -= coupling_factor_z * sqrt(m) * dz[k][k];
                                }
                            }

                            // more coherent-state terms ... these affect diagonals in electronic basis
                            if ( i == j && a == b && n == m + 1 ) {
                                hp[ian_a][jbm_a] += coupling_factor_z * sqrt(n) * tot_dip_z_;
                                hp[ian_b][jbm_b] += coupling_factor_z * sqrt(n) * tot_dip_z_;
                            }

                            if ( i == j && a == b && n == m - 1 ) {
                                hp[ian_a][jbm_a] += coupling_factor_z * sqrt(m) * tot_dip_z_;
                                hp[ian_b][jbm_b] += coupling_factor_z * sqrt(m) * tot_dip_z_;
                            }

                        }
                    }
                }
            }
        }
    }

    // now, couple |0,n+1> and |0,n-1> to |ia,n>
    int off = 2 * o * v * n_photon_states_;
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            for (int n = 0; n < n_photon_states_; n++) {
                int ian_a = (i * v + a      ) * n_photon_states_ + n;
                int ian_b = (i * v + a + o*v) * n_photon_states_ + n;
                // case 1
                if ( n < n_photon_states_ - 1 ) {
                    hp[ian_a][off+n+1] = -coupling_factor_z * sqrt(n+1) * dz[i][a+o];
                    hp[off+n+1][ian_a] = -coupling_factor_z * sqrt(n+1) * dz[i][a+o];
                    hp[ian_b][off+n+1] = -coupling_factor_z * sqrt(n+1) * dz[i][a+o];
                    hp[off+n+1][ian_b] = -coupling_factor_z * sqrt(n+1) * dz[i][a+o];
                }

                // case 2
                if ( n > 0 ) {
                    hp[ian_a][off+n-1] = -coupling_factor_z * sqrt(n) * dz[i][a+o];
                    hp[off+n-1][ian_a] = -coupling_factor_z * sqrt(n) * dz[i][a+o];
                    hp[ian_b][off+n-1] = -coupling_factor_z * sqrt(n) * dz[i][a+o];
                    hp[off+n-1][ian_b] = -coupling_factor_z * sqrt(n) * dz[i][a+o];
                }
            }
        }
    }

    // now, couple |0,n> to all other |0,m> ... really just updating diagonals
    for (int n = 0; n < n_photon_states_; n++) {
        for (int m = 0; m < n_photon_states_; m++) {
            hp[off+n][off+m] += HCavity_z->pointer()[n][m];
        }
    }

    std::shared_ptr<Matrix> eigvec (new Matrix((2*o*v+1)*n_photon_states_,(2*o*v+1)*n_photon_states_));
    std::shared_ptr<Vector> eigval (new Vector((2*o*v+1)*n_photon_states_));
    ham->diagonalize(eigvec,eigval);

    eigval->print();

/*
    for (int i = 0; i < 50; i++) {
        double photon_weight = 0.0;
        for (int j = off + 1; j < off + n_photon_states_; j++) {
            double dum = eigvec->pointer()[j][i];
            photon_weight += dum*dum;
        }
        //printf(" %20.12lf %20.12lf",energy_ + eigval->pointer()[i],photon_weight);
        printf(" %20.12lf %20.12lf\n",energy_ + eigval->pointer()[i],photon_weight);
    }
    printf("\n");
    fflush(stdout);
*/
    
    // print orbital energies
    //epsilon_a_->print();

    // ok, try solving iteratively with davison
    std::shared_ptr<Nonsym_DavidsonSolver> david (new Nonsym_DavidsonSolver());

    // dimension of the problem
    int N = (2*o*v+1)*n_photon_states_;

    // number of desired roots
    int M = options_.get_int("NUM_ROOTS");

    // maximum subspace size
    int maxdim = 8 * M;

    // (approximate) diagonal of Hamiltonian
    double * Hdiag = build_hamiltonian_diagonals();

    std::shared_ptr<Matrix> rel ( new Matrix(M,N) );
    std::shared_ptr<Matrix> rer ( new Matrix(M,N) );
    std::shared_ptr<Vector> revals ( new Vector(M) );

    double * revalp = revals->pointer();
    double ** rerp = rer->pointer();
    double ** relp = rel->pointer();

    int init_dim = 4 * M;
    bool use_residual_norm = true;
    double residual_norm = 1e-5;

    // call eigensolver
    david->solve(Hdiag,
                 N,
                 M,
                 revalp,
                 rerp,
                 relp,
                 [&](int N, int maxdim, int L, double **Q, double **sigmar, double **sigmal) {
                        build_sigma(N,maxdim,L, Q,sigmar,sigmal);
                 }, maxdim, init_dim, residual_norm, use_residual_norm);

    free(Hdiag);


    return 0.0;
}

void PolaritonicUTDDFT::build_sigma(int N, int maxdim, int L, double **Q, double **sigmar, double **sigmal){

    for (int i = 0; i < N; i++) {
        for (int I = 0; I < maxdim; I++) {
            sigmar[i][I] = 0.0;
            sigmal[i][I] = 0.0;

        }
    }

    int o = nalpha_;
    int v = nso_ - nalpha_;

    double * c = (double*)malloc(N*sizeof(double));
    double * s = (double*)malloc(N*sizeof(double));

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

    double ** dx = dipole_[0]->pointer();
    double ** dy = dipole_[1]->pointer();
    double ** dz = dipole_[2]->pointer();

    std::vector<std::shared_ptr<Matrix> > Vx;
    std::vector<std::shared_ptr<Matrix> > Dx;

    // J/K-like contributions
    for (int I = 0; I < L; I++) {

        // unpack a vector
        for (int j = 0; j < N; j++) {
            c[j] = Q[I][j];
        }

        // push density matrices on JK object
        // we need density matrices for alpha and beta, as well as
        // alpha and beta plus a photon

        // Da(mu,nu) = ca(j,b) Ca(mu,j) Ca(nu,b) = Ca(mu,j) Ca'(nu,j)
        // Ca'(nu,j) = ca(j,b) Ca(nu,b)

        // Db(mu,nu) = cb(j,b) Cb(mu,j) Cb(nu,b) = Cb(mu,j) Cb'(nu,j)
        // Cb'(nu,j) = cb(j,b) Cb(nu,b)

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

            // beta
            std::shared_ptr<Matrix> crb (new Matrix(Cb_) );
            std::shared_ptr<Matrix> clb (new Matrix(Cb_) );

            clp = clb->pointer();
            crp = crb->pointer();

            crb->zero();
            clb->zero();

            for (int mu = 0; mu < nso_; mu++) {
                for (int i = 0; i < o; i++) {

                    // left is plain orbitals
                    clp[mu][i] = cbp[mu][i];

                    // right is modified orbitals
                    double dum = 0.0;
                    for (int a = 0; a < v; a++) {
                        int ian_b = (i * v + a + o*v) * n_photon_states_ + n;
                        dum += cbp[mu][a+o] * c[ian_b];
                    }
                    crp[mu][i] = dum;

                }
            }

            // push alpha orbitals onto JK object
            C_left.push_back(clb);
            C_right.push_back(crb);

            // build pseudo densities for xc contribution
            if (needs_xc_) {

                auto Dx_a = linalg::doublet(cla, cra, false, true);
                auto Dx_b = linalg::doublet(clb, crb, false, true);
                //std::shared_ptr<Matrix> Dx_a (new Matrix(nso_,nso_));
                //std::shared_ptr<Matrix> Dx_b (new Matrix(nso_,nso_));
                //C_DGEMM('n','t',nso_,nso_,o,1.0,&(cla->pointer()[0][0]),nso_,&(cra->pointer()[0][0]),nso_,0.0,&(Dx_a->pointer()[0][0]),nso_);
                //C_DGEMM('n','t',nso_,nso_,o,1.0,&(clb->pointer()[0][0]),nso_,&(crb->pointer()[0][0]),nso_,0.0,&(Dx_b->pointer()[0][0]),nso_);


                Vx.push_back(std::make_shared<Matrix>("Vax temp", Dx_a->rowspi(), Dx_a->colspi(), Dx_a->symmetry()));
                Vx.push_back(std::make_shared<Matrix>("Vbx temp", Dx_b->rowspi(), Dx_b->colspi(), Dx_b->symmetry()));
                Dx.push_back(Dx_a);
                Dx.push_back(Dx_b);
    
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
    double * eb = epsilon_b_->pointer();

    double * tmp_a = (double*)malloc(o*v*sizeof(double));
    double * tmp_b = (double*)malloc(o*v*sizeof(double));

    for (int I = 0; I < L; I++) {

        // unpack a vector
        for (int j = 0; j < N; j++) {
            c[j] = Q[I][j];
        }

        for (int n = 0; n < n_photon_states_; n++) {

            // a,b <- a coulomb
            std::shared_ptr<Matrix> sa (new Matrix(jk_->J()[count]));
            std::shared_ptr<Matrix> sb (new Matrix(jk_->J()[count]));

            // a,b <- b coulomb
            sa->add(jk_->J()[count+1]);
            sb->add(jk_->J()[count+1]);

            // xc?
            // a <- a
            // b <- b
            if ( needs_xc_ ) {
                sa->axpy(1.0,Vx[count]);
                //sb->axpy(1.0,Vx[count]);

                //sa->axpy(1.0,Vx[count+1]);
                sb->axpy(1.0,Vx[count+1]);
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

            for (int i = 0; i < o; i++) {
                for (int a = 0; a < v; a++) {
                    int ian_a = (i * v + a      ) * n_photon_states_ + n;
                    int ian_b = (i * v + a + o*v) * n_photon_states_ + n;

                    sigmar[ian_a][I] = sap[i][a+o];
                    sigmar[ian_b][I] = sbp[i][a+o];
                                                  
                    sigmal[ian_a][I] = sap[i][a+o];
                    sigmal[ian_b][I] = sbp[i][a+o];
                }
            }
        }

        // singles diagonals
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                for (int n = 0; n < n_photon_states_; n++) {
                    int ian_a = (i * v + a      ) * n_photon_states_ + n;
                    int ian_b = (i * v + a + o*v) * n_photon_states_ + n;

                    sigmar[ian_a][I] += c[ian_a] * (ea[a+o] - ea[i] + n * cavity_frequency_[2]);
                    sigmar[ian_b][I] += c[ian_b] * (eb[a+o] - eb[i] + n * cavity_frequency_[2]);

                    sigmal[ian_a][I] += c[ian_a] * (ea[a+o] - ea[i] + n * cavity_frequency_[2]);
                    sigmal[ian_b][I] += c[ian_b] * (eb[a+o] - eb[i] + n * cavity_frequency_[2]);
                }
            }
        }

        // now, |0,n> diagonals
        int off = 2 * o * v * n_photon_states_;
        for (int n = 0; n < n_photon_states_; n++) {
            sigmar[off+n][I] += n * cavity_frequency_[2] * c[off+n];
            sigmal[off+n][I] += n * cavity_frequency_[2] * c[off+n];
        }

        // dipole self energy
        for (int n = 0; n < n_photon_states_; n++) {

            // J-like contribution from dipole self energy
            double dipole_Ja = 0.0;
            double dipole_Jb = 0.0;
            for (int i = 0; i < o; i++) {
                for (int a = 0; a < v; a++) {
                    int ian_a = (i * v + a      ) * n_photon_states_ + n;
                    int ian_b = (i * v + a + o*v) * n_photon_states_ + n;

                    dipole_Ja += c[ian_a] * dz[i][a+o];
                    dipole_Jb += c[ian_b] * dz[i][a+o];
                }
            }

            // intermediate for K-like contribution from dipole self energy
            memset((void*)tmp_a,'\0',o*v*sizeof(double));
            memset((void*)tmp_b,'\0',o*v*sizeof(double));
            for (int i = 0; i < o; i++) {
                for (int a = 0; a < v; a++) {
                    double dum_a = 0.0;
                    double dum_b = 0.0;
                    for (int j = 0; j < o; j++) {
                        int jan_a = (j * v + a      ) * n_photon_states_ + n;
                        int jan_b = (j * v + a + o*v) * n_photon_states_ + n;

                        dum_a += c[jan_a] * dz[i][j];
                        dum_b += c[jan_b] * dz[i][j];
                    }

                    tmp_a[a*o+i] = dum_a;
                    tmp_b[a*o+i] = dum_b;
                }
            }

            for (int i = 0; i < o; i++) {
                for (int a = 0; a < v; a++) {

                    double dipole_Ja_ia = dipole_Ja * dz[i][a+o];
                    double dipole_Jb_ia = dipole_Jb * dz[i][a+o];

                    double dipole_Ka = 0.0;
                    double dipole_Kb = 0.0;
                    for (int b = 0; b < v; b++) {
                        dipole_Ka += tmp_a[b*o+i] * dz[a+o][b+o];
                        dipole_Kb += tmp_b[b*o+i] * dz[a+o][b+o];
                    }

                    int ian_a = (i * v + a      ) * n_photon_states_ + n;
                    int ian_b = (i * v + a + o*v) * n_photon_states_ + n;

                    sigmar[ian_a][I] += lambda_z * lambda_z * (dipole_Ja_ia + dipole_Jb_ia - dipole_Ka);
                    sigmar[ian_b][I] += lambda_z * lambda_z * (dipole_Ja_ia + dipole_Jb_ia - dipole_Kb);

                    sigmal[ian_a][I] += lambda_z * lambda_z * (dipole_Ja_ia + dipole_Jb_ia - dipole_Ka);
                    sigmal[ian_b][I] += lambda_z * lambda_z * (dipole_Ja_ia + dipole_Jb_ia - dipole_Kb);
                }
            }
        }

        // coupling terms
        for (int n = 0; n < n_photon_states_; n++) {

            double sqrt_n = sqrt(n);

            for (int i = 0; i < o; i++) {
                for (int a = 0; a < v; a++) {

                    int ian_a = (i * v + a      ) * n_photon_states_ + n;
                    int ian_b = (i * v + a + o*v) * n_photon_states_ + n;

                    double dum_a = 0.0;
                    double dum_b = 0.0;

                    // d(n,m+1)
                    int m = n - 1;
                    if ( m > -1 ) {
                        // d(i,j)
                        for (int b = 0; b < v; b++) {
                            int ibm_a = (i * v + b      ) * n_photon_states_ + m;
                            int ibm_b = (i * v + b + o*v) * n_photon_states_ + m;
                            dum_a -= sqrt_n * c[ibm_a] * dz[a+o][b+o];
                            dum_b -= sqrt_n * c[ibm_b] * dz[a+o][b+o];
                        }
                        // d(a,b)
                        for (int j = 0; j < o; j++) {
                            int jam_a = (j * v + a      ) * n_photon_states_ + m;
                            int jam_b = (j * v + a + o*v) * n_photon_states_ + m;
                            dum_a += sqrt_n * c[jam_a] * dz[i][j];
                            dum_b += sqrt_n * c[jam_b] * dz[i][j];
                        }
                        // d(i,j) d(a,b)
                        int iam_a = (i * v + a      ) * n_photon_states_ + m;
                        int iam_b = (i * v + a + o*v) * n_photon_states_ + m;
                        for (int k = 0; k < o; k++) {
                            dum_a -= sqrt_n * c[iam_a] * dz[k][k];
                            dum_b -= sqrt_n * c[iam_b] * dz[k][k];
                        }
                        // more coherent-state terms ... these affect diagonals in the electronic basis
                        dum_a += sqrt_n * tot_dip_z_ * c[iam_a];
                        dum_b += sqrt_n * tot_dip_z_ * c[iam_b];
                    }

                    // d(n,m-1)
                    m = n + 1;
                    if ( m < n_photon_states_ ) {
                        double sqrt_m = sqrt(m);
                        // d(i,j)
                        for (int b = 0; b < v; b++) {
                            int ibm_a = (i * v + b      ) * n_photon_states_ + m;
                            int ibm_b = (i * v + b + o*v) * n_photon_states_ + m;
                            dum_a -= sqrt_m * c[ibm_a] * dz[a+o][b+o];
                            dum_b -= sqrt_m * c[ibm_b] * dz[a+o][b+o];
                        }
                        // d(a,b)
                        for (int j = 0; j < o; j++) {
                            int jam_a = (j * v + a      ) * n_photon_states_ + m;
                            int jam_b = (j * v + a + o*v) * n_photon_states_ + m;
                            dum_a += sqrt_m * c[jam_a] * dz[i][j];
                            dum_b += sqrt_m * c[jam_b] * dz[i][j];
                        }
                        // d(i,j) d(a,b)
                        int iam_a = (i * v + a      ) * n_photon_states_ + m;
                        int iam_b = (i * v + a + o*v) * n_photon_states_ + m;
                        for (int k = 0; k < o; k++) {
                            dum_a -= sqrt_m * c[iam_a] * dz[k][k];
                            dum_b -= sqrt_m * c[iam_b] * dz[k][k];
                        }
                        // more coherent-state terms ... these affect diagonals in the electronic basis
                        dum_a += sqrt_m * tot_dip_z_ * c[iam_a];
                        dum_b += sqrt_m * tot_dip_z_ * c[iam_b];
                    }

                    sigmar[ian_a][I] += coupling_factor_z * dum_a;
                    sigmar[ian_b][I] += coupling_factor_z * dum_b;

                    sigmal[ian_a][I] += coupling_factor_z * dum_a;
                    sigmal[ian_b][I] += coupling_factor_z * dum_b;

                }
            }

        }

        // lastly, couple |0,n+1> and |0,n-1> to |ia,n>
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                for (int n = 0; n < n_photon_states_; n++) {
                    int ian_a = (i * v + a      ) * n_photon_states_ + n;
                    int ian_b = (i * v + a + o*v) * n_photon_states_ + n;
                    // case 1
                    if ( n < n_photon_states_ - 1 ) {
                        double factor = -coupling_factor_z * sqrt(n+1) * dz[i][a+o];

                        double val_a   = factor * c[ian_a];
                        double val_b   = factor * c[ian_b];
                        double val_np1 = factor * c[off+n+1];

                        sigmar[ian_a][I]   += val_np1;
                        sigmar[ian_b][I]   += val_np1;
                        sigmar[off+n+1][I] += val_a + val_b;

                        sigmal[ian_a][I]   += val_np1;
                        sigmal[ian_b][I]   += val_np1;
                        sigmal[off+n+1][I] += val_a + val_b;

                    }

                    // case 2
                    if ( n > 0 ) {
                        double factor = -coupling_factor_z * sqrt(n) * dz[i][a+o];

                        double val_a   = factor * c[ian_a];
                        double val_b   = factor * c[ian_b];
                        double val_nm1 = factor * c[off+n-1];

                        sigmar[ian_a][I]   += val_nm1;
                        sigmar[ian_b][I]   += val_nm1;
                        sigmar[off+n-1][I] += val_a + val_b;

                        sigmal[ian_a][I]   += val_nm1;
                        sigmal[ian_b][I]   += val_nm1;
                        sigmal[off+n-1][I] += val_a + val_b;
                    }
                }
            }
        }

    }


    free(tmp_a);
    free(tmp_b);

    free(c);
    free(s);

}

double * PolaritonicUTDDFT::build_hamiltonian_diagonals(){

    int o = nalpha_;
    int v = nso_ - nalpha_;

    double * H = (double*)malloc((2*o*v+1)*n_photon_states_*sizeof(double));; 
    memset((void*)H,'\0',(2*o*v+1)*n_photon_states_*sizeof(double));

    double * ea = epsilon_a_->pointer();
    double * eb = epsilon_b_->pointer();

    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            for (int n = 0; n < n_photon_states_; n++) {
                int ian_a = (i * v + a      ) * n_photon_states_ + n;
                int ian_b = (i * v + a + o*v) * n_photon_states_ + n;
                H[ian_a] = ea[a+o] - ea[i] + n * cavity_frequency_[2];
                H[ian_b] = eb[a+o] - eb[i] + n * cavity_frequency_[2];
            }
        }
    }

    // now, |0,n> diagonals
    int off = 2 * o * v * n_photon_states_;
    for (int n = 0; n < n_photon_states_; n++) {
        H[off+n] = n * cavity_frequency_[2];
    }

    return H;
}


} // End namespaces
