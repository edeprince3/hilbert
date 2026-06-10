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

#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
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
//#include <mkl_blas.h>


// jk object
#include <psi4/libfock/jk.h>

// for dft
#include "psi4/libscf_solver/hf.h"

#include <psi4/psifiles.h>
#include <psi4/libtrans/integraltransform.h>

#include "rtddft.h"
#include "rdfrpa.h"
//#include "_deps/libsdp_external-src/src/blas_helper.h"

#include <misc/blas.h>
#include <misc/hilbert_psifiles.h>
#include <misc/omp.h>
#include <misc/threeindexintegrals.h>
#include <misc/nonsym_davidson_solver.h>
//#include <mkl_blas.h>

//extern "C" {
//void dpotrf(char* uplo, int* n, double* a, int* lda, int* info);
//void dpotrs(char* uplo, int* n, int* nrhs, double* a, int* lda,
//                 double* b, int* ldb, int* info);
//}

using namespace psi;
using namespace fnocc;

namespace hilbert{ 

PolaritonicRDFRPA::PolaritonicRDFRPA(std::shared_ptr<Wavefunction> reference_wavefunction, Options& options_, std::shared_ptr<Wavefunction> dummy_wfn):
    PolaritonicHF(reference_wavefunction,options_) {
    common_init(dummy_wfn);
}

PolaritonicRDFRPA::~PolaritonicRDFRPA() {

    free(diag_);
    free(siap_);


}

void PolaritonicRDFRPA::common_init(std::shared_ptr<Wavefunction> dummy_wfn) {

    outfile->Printf("\n");
    outfile->Printf("\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    Polaritonic Restricted DFRPA                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf("\n");


    // ensure scf_type df
    if ( options_.get_str("SCF_TYPE") != "DF" && options_.get_str("SCF_TYPE") != "CD" ) {
        throw PsiException("polaritonic RPA only works with scf_type df and cd for now",__FILE__,__LINE__);
    }

    // ensure running in c1 symmetry
    if ( reference_wavefunction_->nirrep() > 1 ) {
        throw PsiException("polaritonic RPA only works with c1 symmetry for now.",__FILE__,__LINE__);
    }

    // ensure closed shell
    if ( nalpha_ != nbeta_ ) {
        throw PsiException("polaritonic RPA only works with nalpha = nbeta (for now)",__FILE__,__LINE__);
    }

    // get primary basis:
    std::shared_ptr<BasisSet> primary = reference_wavefunction_->get_basisset("ORBITAL");

    o_ = nalpha_;
    v_ = nso_ - nalpha_;

    nQ_ = 0;

    if ( options_.get_str("SCF_TYPE") == "DF" ) {

        // get auxiliary basis:
        std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");

        // total number of auxiliary basis functions
        nQ_ = auxiliary->nbf();
        siap_ = (double*)malloc(o_*v_*nQ_*sizeof(double));
        memset((void*)siap_,'\0',o_*v_*nQ_*sizeof(double));

        std::shared_ptr<DFTensor> DF (new DFTensor(primary,auxiliary,Ca_,o_,v_,o_,v_,options_));

        std::shared_ptr<Matrix> tmpoo = DF->Qoo();
        std::shared_ptr<Matrix> tmpov = DF->Qov();
        std::shared_ptr<Matrix> tmpvv = DF->Qvv();

        double ** Qoo = tmpoo->pointer();
        double ** Qov = tmpov->pointer();
        double ** Qvv = tmpvv->pointer();


        for (size_t Q = 0; Q < nQ_; Q++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t a = 0; a < v_; a++) {
                    size_t ai = i * v_ + a;
                    siap_[Q * o_ * v_ + i * v_ + a] = *&Qov[Q][ai];
                }
            }
        }
        //tmpov->print();
        tst_ = (std::shared_ptr<Matrix>)(new Matrix(nQ_,o_*v_));
        tst_ = tmpov;
        //DEBUG
        /*
        std::shared_ptr<Matrix> tmp2 = (std::shared_ptr<Matrix>)(new Matrix(nQ_,nQ_));
        tmp2->gemm(false,true,1.0,tmpov,tmpov,0.0);
        tmp2->print();
        std::shared_ptr<Matrix> evecs = (std::shared_ptr<Matrix>)(new Matrix(nQ_,nQ_));
        std::shared_ptr<Vector> evals = (std::shared_ptr<Vector>) (new Vector(nQ_));
        tmp2->diagonalize(evecs,evals);
        evals->print();
        double erpa = 0.0;
        for (size_t i = 0; i < nQ_; i++) {
            erpa += evals->pointer()[i];
        }
        outfile->Printf("    =====> RPA correlation energy <=====\n");
        outfile->Printf("%20.12lf \n",erpa);
         */


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

        double * Qoo = (double*)malloc(o_*o_*nQ_*sizeof(double));
        double * Qov = (double*)malloc(o_*v_*nQ_*sizeof(double));
        double * Qvv = (double*)malloc(v_*v_*nQ_*sizeof(double));


        for (size_t Q = 0; Q < nQ_; Q++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    Qoo[Q*o_*o_+i*o_+j] = Qmo[Q*nmo_*(nmo_+1)/2+INDEX(i,j)];
                }
            }
            for (size_t i = 0; i < o_; i++) {
                for (size_t a = 0; a < v_; a++) {
                    Qov[Q*o_*v_+i*v_+a] = Qmo[Q*nmo_*(nmo_+1)/2+INDEX(i,a+o_)];
                }
            }
            for (size_t a = 0; a < v_; a++) {
                for (size_t b = 0; b < v_; b++) {
                    Qvv[Q*v_*v_+a*v_+b] = Qmo[Q*nmo_*(nmo_+1)/2+INDEX(a+o_,b+o_)];
                }
            }
        }

        free(Qmo);

        //outfile->Printf("\n");
        //outfile->Printf("        Time for integral transformation:  %7.2lf s\n",end-start);
        //outfile->Printf("\n");

        siap_ = (double*)malloc(o_*v_*nQ_*sizeof(double));
        memset((void*)siap_,'\0',o_*v_*nQ_*sizeof(double));

        tst_ = (std::shared_ptr<Matrix>)(new Matrix(nQ_,o_*v_));
        for (size_t Q = 0; Q < nQ_; Q++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t a = 0; a < v_; a++) {
                    siap_[Q*o_*v_+i*v_+a] = Qov[Q*o_*v_+i*v_+a];
                    tst_->pointer()[Q][i*v_+a] = Qov[Q*o_*v_+i*v_+a];
                }
            }
        }

        free(Qoo);
        free(Qov);
        free(Qvv);
    }

    // orbital energies
    epsilon_a_ = std::make_shared<Vector>(nmopi_);
    epsilon_a_->copy(*reference_wavefunction_->epsilon_a().get());
    epsilon_b_ = std::make_shared<Vector>(nmopi_);
    epsilon_b_->copy(*reference_wavefunction_->epsilon_b().get());

    diag_ = (double*)malloc(o_*v_*sizeof(double));

    for (size_t i = 0; i < o_; i++) {
        for (size_t a = 0; a < v_; a++) {
            size_t ai = i * v_ + a;
            diag_[ai] = epsilon_a_->pointer()[a + o_] - epsilon_a_->pointer()[i];
        }
    }
}

double PolaritonicRDFRPA::compute_energy() {

    bool batched = options_.get_bool("BATCHED_DFRPA");

    // if batched option is set, compute RPA energy using a batched approach that avoids explicit diagonalization of Qmat_ for each grid point
    if (batched) {
         outfile->Printf("Batched DFRPA version (not fully tested yet, use with caution!)\n");
         double ecrpa =  compute_energy_batched();
          return ecrpa;
    }

    // otherwise compute RPA energy using unbatched approach that diagonalizes Qmat_ for each grid point
    return compute_energy_unbatched();

}

double PolaritonicRDFRPA::compute_energy_unbatched() {

    outfile->Printf("\n");
    outfile->Printf("    No. basis functions:            %5i\n",nso_);
    outfile->Printf("    No. electrons:                  %5i\n",nalpha_ + nbeta_);
    if ( options_.get_bool("QED_USE_RELAXED_ORBITALS") ) {
        outfile->Printf("    Relaxed orbitals are used. K terms included in g\n");
    }
    if ( options_.get_bool("RPA_EXCHANGE") ) {
        outfile->Printf("    Exchange added to A and B in RPA kernel\n");
    }

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


    // build diagonal matrix for quadrature
    std::shared_ptr<Vector> diag_matrix = build_diag_matrix();

    // build grid points and weights for quadrature
    // Curtis–Clenshaw quadrature based on diagonal approximation to RPA energy expression
    auto grid = build_grid();
    // alternative: minimax quadrature based on diagonal approximation to RPA energy expression
    // needs debugging! 
//    double tol = 1e-6;
//    auto grid = build_minimax_grid(tol);

    // print out grid points and weights for debugging
    if (n_output_ > 2)  {
        outfile->Printf("     =====> orbital energy differences <=====\n");
        for (size_t i = 0; i < o_ * v_ ; i++) {
            outfile->Printf("%20.12lf \n",diag_[i]);
        }
        outfile->Printf("     =====> (ia|ia) elements <=====\n");
        for (size_t i = 0; i < o_ * v_ ; i++) {
            outfile->Printf("%20.12lf \n",diag_matrix->pointer()[i]);
        }
        outfile->Printf("     =====> quadrature weights <=====\n");
        for (size_t i = 0; i < n_grid_points_ ; i++) {
            outfile->Printf("%20.12lf \n",grid->pointer()[i][1]);
        }
        outfile->Printf("     =====> quadrature positions <=====\n");
        for (size_t i = 0; i < n_grid_points_ ; i++) {
            outfile->Printf("%20.12lf \n",grid->pointer()[i][0]);
        }
    }

    // Static DSE contribution
    std::vector<double> mu_ai(o_*v_);

    for (size_t i = 0; i < o_; i++) {
        for (size_t a = 0; a < v_; a++) {
            size_t ai = i * v_ + a;
            mu_ai[ai] = dz[i][a+o_];
        }
    } 

    // prefactor for DSE contribution
    double omega_c = cavity_frequency_[2];
    double alpha_DSE = lambda_z * lambda_z / (2.0 * omega_c * omega_c);

    // build synthetic DSE matrix in auxiliary basis
    std::vector<double> S_DSE(o_*v_);
    for (size_t ai = 0; ai < o_*v_; ++ai)
      // factor 0.5 cancel later scaling of 3-index integrals with factor 2
      //S_DSE[ai] = std::sqrt(0.5 * alpha_DSE) * mu_ai[ai];
      S_DSE[ai] = lambda_z * mu_ai[ai];


    //extend DF tensor by one row
    auto tst_ext = std::make_shared<Matrix>(nQ_ + 1, o_*v_);
    tst_ext->zero();


    // manually copy original DF rows
    for (size_t Q = 0; Q < nQ_; ++Q)
      for (size_t ai = 0; ai < o_*v_; ++ai)
        tst_ext->set(Q, ai, tst_->get(Q, ai));

    outfile->Printf("DEBUG: tst_ext nrow = %d, ncol = %d (expected %d x %d)\n",
                tst_ext->nrow(), tst_ext->ncol(), nQ_ + 1, o_*v_);

    // add DSE row
    for (size_t ai = 0; ai < o_*v_; ++ai)
      tst_ext->set(nQ_, ai, S_DSE[ai]);


    // DEBUG
    //outfile->Printf(" Do we get here? \n");

    // after filling S_DSE and tst_ext
    //outfile->Printf("DEBUG: built S_DSE and tst_ext\n");
    //outfile->Printf("DEBUG: tst_ dims      = %d x %d\n", tst_->nrow(), tst_->ncol());
    //outfile->Printf("DEBUG: tst_ext dims   = %d x %d\n", tst_ext->nrow(), tst_ext->ncol());
    //outfile->Printf("DEBUG: expected dims  = %d x %d\n", nQ_ + 1, o_*v_);

    // allocate temporary arrays for RPA calculation
    int * int1_ = (int*)malloc(nQ_*sizeof(int));
    memset((void*)int1_,'\0',nQ_*sizeof(int));
    std::shared_ptr<Matrix> tstcopy = (std::shared_ptr<Matrix>)(new Matrix(nQ_ + 1, o_*v_));

    // compute Qmat_ for each grid point and accumulate Ecrpa
    double ecrpa = 0.0;
    // loop over grid points
    for (int i = 0; i < n_grid_points_; ++i) {

        std::shared_ptr<Matrix> tmp2 = (std::shared_ptr<Matrix>)(new Matrix(nQ_ +1 ,nQ_ + 1));
        tmp2->zero();
        // copy 3-index integrals to local variable tstcopy
        tstcopy->zero();
        tstcopy->copy(tst_ext);
        // scale 3 index integrals with e-e contribution
        for (size_t ivir = 0; ivir < v_ ; ivir++) {
            for (size_t iocc = 0; iocc < o_ ; iocc++) {
                size_t ai = iocc * v_ + ivir;
                double scal = sqrt(diag_[ai] / (diag_[ai] *diag_[ai] + grid->pointer()[i][0] * grid->pointer()[i][0]));
                if (n_output_ > 2) {
                    outfile->Printf("   =====> scaling factor for gridpoint %i, orbital %i <=====\n", i, ai);
                    outfile->Printf("%20.12lf \n", scal);
                }
                tstcopy->scale_column(0,ai,scal);
            }
        }

        if (n_output_ > 2) {
            outfile->Printf("   =====> tstcopy after scaling for gridpoint %i <=====\n");
            tstcopy->print();
        }

        // accumulate contribution on Qmat_
        // scaling factor for 3 index integrals
        double scal2 = 2.0*2.0;
        // this is the e-e contribution to Q + static DSE
        tmp2->gemm(false,true,scal2,tstcopy,tstcopy,1.0);

        if (n_output_ > 2) {
            outfile->Printf("   =====> Qmat_ for gridpoint %i <=====\n", i);
            tmp2->print();
        }

        // compute trace of Qmat_ (unscreened)
        double trace_Q = 0.0;
        for (size_t Q = 0; Q < nQ_ + 1; Q++) {
            trace_Q += tmp2->pointer()[Q][Q];
        }
        // -Tr(Q)
        double ecrpa_pnt = 0.0;
        ecrpa_pnt = -trace_Q;

        // add QED-kernel contribution to Q: static DSE term + photon coupling term
        double freq = grid->pointer()[i][0]; // omega
        double denom = std::sqrt(cavity_frequency_[2]*cavity_frequency_[2] + grid->pointer()[i][0]*grid->pointer()[i][0]);
        double Sph = freq / denom;

        // scale last row and last column
        for (size_t Q = 0; Q < nQ_ + 1; ++Q) {
          tmp2->pointer()[Q][nQ_] *= Sph;   // column
          tmp2->pointer()[nQ_][Q] *= Sph;   // row
        }
        
        // accumulate the correlation energy
        auto evecs = std::make_shared<Matrix>(nQ_ + 1, nQ_ + 1);
        auto evals = std::make_shared<Vector>(nQ_ + 1);

        // compute RPA correlation energy
        // form Q_eff + 1
        for (size_t Q = 0; Q < nQ_ + 1; Q++) {
            //ecrpa_pnt -= tmp2->pointer()[Q][Q];
            tmp2->pointer()[Q][Q] += 1.0;
        }
//        tmp2->print();
        // diagonalize Qmat_ + 1
        evals->zero();
        evecs->zero();
        tmp2->diagonalize(evecs,evals);
        //evals->print();

        // accumulate log of eigenvalues for RPA energy
        double ecln = 0.0;
        for (size_t Q = 0; Q < nQ_ + 1; Q++) {
            ecln +=  log(evals->pointer()[Q]);
        }
        // Tr(ln(Q_eff+1))
        ecrpa_pnt +=  ecln;
        // scale with quadrature weight and prefactors
        ecrpa_pnt = ecrpa_pnt/(4*M_PI) * grid->pointer()[i][1];
        // accumulate contribution to total RPA correlation energy
        ecrpa += ecrpa_pnt;

    }
    outfile->Printf("     =====> RPA correlation energy <=====\n");
    //printf("%20.12lf\n", ecrpa);
    outfile->Printf("%20.12lf \n",ecrpa);

    return ecrpa;

}

// function to build the diagonal matrix
std::shared_ptr<Vector> PolaritonicRDFRPA::build_diag_matrix() {

    std::shared_ptr<Vector> D (new Vector(o_ * v_));

    double gvals;
    for (size_t a = 0; a < v_; a++) {
        for (size_t i = 0; i < o_; i++) {
            size_t ai = i * v_ + a;
            gvals = 0.0;
            for (size_t Q = 0; Q < nQ_; Q++) {
                gvals += siap_[Q*o_*v_+ ai] * siap_[Q*o_*v_+ ai];
            }
            D->pointer()[ai] = 2.0 * gvals;
        }
    }
    return D;
}

// function to build the grid points and weights for quadrature
std::shared_ptr<Matrix> PolaritonicRDFRPA::build_grid() {

        std::shared_ptr<Matrix> grid(new Matrix(n_grid_points_, 2));

        std::shared_ptr<Vector> diag_matrix = build_diag_matrix();
        // set up quadrature using diagonal approximation to RPA energy expression

        // Evaluate e_approx
        double e_approx = 0.0;
        for (size_t a = 0; a < v_; a++) {
            for (size_t i = 0; i < o_; i++) {
                size_t ai = i * v_ + a;
                e_approx +=
                        sqrt(diag_[ai] * diag_[ai] + 2.0 * diag_[ai] * diag_matrix->pointer()[ai]) - diag_[ai] - diag_matrix->pointer()[ai];
            }
        }
        e_approx = 0.5 * e_approx;

        // setup
        double weight[n_grid_points_];
        double place[n_grid_points_];
        double t_fac = M_PI / (2.0 * n_grid_points_);
        for (size_t i = 0; i < n_grid_points_; i++) {
            double t_j = t_fac * (i + 1);
            weight[i] = M_PI / (n_grid_points_ * sin(t_j) * sin(t_j));
            place[i] = 1 / tan(t_j);
        }
        weight[n_grid_points_ - 1] = 0.5 * weight[n_grid_points_ - 1];

        // Newton iteration
        double egrid = 0.0;
        double egrid_deriv = 0.0;
        do {
            int ii = 0;
            egrid = 0.0;
            egrid_deriv = 0.0;
            for (size_t i = 0; i < n_grid_points_; i++) {
                double freq = scal_fac_ * place[i];
                double fsq = freq * freq;
                double egrid_tmp = 0.0; double egrid_deriv_tmp = 0.0;
                for (size_t j = 0; j < o_ * v_; j++) {
                    double tmp = diag_matrix->pointer()[j];
                    double tmp2 = diag_[j] / (diag_[j] * diag_[j] + fsq);
                    double tmp3 = -2.0 * tmp2/(diag_[j] * diag_[j] + fsq) * freq;
                    double tmp4 = 2.0 * tmp2 * tmp;
                    egrid_tmp +=  tmp4 - log(1+tmp4);
                    egrid_deriv_tmp +=  2.0 * tmp3 * tmp * tmp4 / (1.0 + tmp4);

                }
                egrid +=  egrid_tmp * weight[i];
                egrid_deriv +=  egrid_deriv_tmp * weight[i] * place[i];
            }
            egrid = -0.25* egrid * scal_fac_ / M_PI;
            egrid_deriv = -0.25 * egrid_deriv * scal_fac_ / M_PI;
            egrid_deriv +=  egrid / scal_fac_;

            // Newton step
            double step = (e_approx - egrid) / egrid_deriv;
            scal_fac_ = scal_fac_ + step;

        } while (abs(e_approx - egrid) > 1e-13);

        if (n_output_ > 0) {

            outfile->Printf("     =====> Number of grid points <=====\n");
            outfile->Printf("%10d \n",n_grid_points_);
            outfile->Printf("     =====> Scaling factor <=====\n");
            outfile->Printf("%20.12lf \n",scal_fac_);
            outfile->Printf("     =====> Sensitivity Parameter <=====\n");
            outfile->Printf("%20.12lf \n",egrid_deriv);
        }
        // scale grid points and weights
        for (size_t i = 0; i < n_grid_points_; i++) {
            place[i] = scal_fac_ * place[i];
            weight[i] = scal_fac_ * weight[i];
        }
        // put in reverse order to compare to turbomole
        for (size_t i = 0; i < n_grid_points_; i++) {
            grid->pointer()[n_grid_points_-1-i][0] = place[i];
            grid->pointer()[n_grid_points_-1-i][1] = weight[i];
        }

        return grid;

}

// function to build the grid points and weights for quadrature using minimax algorithm
// should be more efficient than Curtis-Clenshaw
// adapts to meet tolerance
// based on the GreenX formulation

std::shared_ptr<Matrix> PolaritonicRDFRPA::build_minimax_grid(double tol) {

    // 1. Compute Δ_max
    double Delta_max = 0.0;
    for (size_t ai = 0; ai < o_*v_; ++ai)
        Delta_max = std::max(Delta_max, diag_[ai]);

    if (Delta_max == 0.0) Delta_max = 1.0; // avoid division by zero

    // 2. Precompute diagonal-approximation reference energy
    std::shared_ptr<Vector> diag_matrix = build_diag_matrix();
    double E_ref = 0.0;
    for (size_t ai = 0; ai < o_*v_; ++ai) {
        double D = diag_[ai];
        double K = diag_matrix->pointer()[ai];
        E_ref += std::sqrt(D*D + 2.0*D*K) - D - K;
    }
    E_ref *= 0.25;

    // 3. Adaptive loop over N
    double u_max = 4.0; // max u value for grid (maps to large frequencies)
    size_t N = 4;  // starting guess
    while (true) {
        std::vector<double> omega;
        std::vector<double> weight;

        // 3a. Adaptive loop over N (uniform spacing step 'h' mapped
        // exponentially to [0,∞))
        // expand the limits symmetrically around zero to capture negative
        // frequencies as well
        double h = u_max / static_cast<int>(N); // uniform spacing in t-space
        int steps = static_cast<int>(N);
        for (int j = -steps; j <= steps; ++j) {
            double t_j = j * h;
            double sinh_t = std::sinh(t_j);
            double cosh_t = std::cosh(t_j);
            double omega_j = Delta_max * std::exp(sinh_t);
            double weight_j = h * Delta_max * std::exp(sinh_t) * cosh_t;
            omega.push_back(omega_j);
            weight.push_back(weight_j);
        }


        // 3b. Evaluate diagonal-approximation energy on this grid
        double E_test = 0.0;
        size_t grid_size = omega.size();

        for (size_t j = 0; j < grid_size; ++j) {
            double w = weight[j];
            double f2 = omega[j] * omega[j];
            double contrib = 0.0;

            for (size_t ai = 0; ai < o_*v_; ++ai) {
                double D = diag_[ai];
                double K = diag_matrix->pointer()[ai];
                double argument = (2.0 * D * K) / (D*D + f2);

                // use series expansion for small argument to avoid numerical
                // issues
                if (std::abs(argument) < 1e-4) {
                  double term2 = argument * argument;
                  double term3 = term2 * argument;
                  double term4 = term3 * argument;
                  contrib += - term2/2.0 + term3/3.0 - term4/4.0;
                } else {
                contrib += std::log(1.0 + argument) - argument;
                //outfile->Printf("D=%lf  K=%lf\n", D, K);
                }
            }
            E_test += w * contrib;
        }
        E_test *= 0.25 / M_PI;

        double current_error = std::abs(E_test - E_ref);

        outfile->Printf("Minimax Adaptive: N=%zu  E_ref=%20.12lf  E_test=%20.12lf  err=%20.12lf\n", N, E_ref, E_test, current_error);

        // 3c. Check error
        if (current_error < tol)
            break;

        // Otherwise increase N
        N += 8; // increase by 8 points (4 on each side) for next iteration
        if (N > 64)
            throw PsiException("grid construction failed to converge within reasonable N", __FILE__, __LINE__);

    }

    // 4. Build final grid matrix
    
    auto grid_matrix = std::make_shared<Matrix>("Optimized Grid", N, 2);

    // Populate grid matrix with points and weights, putting in reverse order to
    // match Curtis–Clenshaw convention
    double h = u_max /((N-1)/2); // uniform spacing in t-space
    int half_steps = (N - 1) / 2;
    int idx = 0;

    for (int j = -half_steps; j <= half_steps; ++j) {
        double t_j = j * h;
        double sinh_t = std::sinh(t_j);
        double cosh_t = std::cosh(t_j);
        double omega_j = Delta_max * std::exp(sinh_t);
        double weight_j = h * Delta_max * std::exp(sinh_t) * cosh_t;

        grid_matrix->pointer()[N-1-j][0] = omega_j;
        grid_matrix->pointer()[N-1-j][1] = weight_j;
        idx++;
    }

    return grid_matrix;
}

// function to compute batchsize based on available memory and size of system
//
int PolaritonicRDFRPA::determine_optimal_batch_size() {
    // 1. Get total memory allocated to Psi4 (in bytes)
    // df_helper_->get_memory() returns the size in bytes (e.g., 2GB = 2,147,483,648)
    size_t total_psi4_mem = df_helper_->get_memory();

    // 2. Set aside a safety margin (e.g., 25%) for static matrices and overhead
    // This prevents the batch calculation from pushing the node into out-of-memory paging
    double safety_factor = 0.75;
    size_t available_mem = static_cast<size_t>(total_psi4_mem * safety_factor);

    // 3. Compute memory footprint of a single (ia) column block
    int nQ_eff = nQ_ + 1;
    size_t bytes_per_ia = 2 * nQ_eff * sizeof(double); 

    // 4. Calculate how many ia elements can fit in that pool
    size_t calculated_batch = available_mem / bytes_per_ia;

    // 5. Impose sanity bounds
    int total_ov = o_ * v_;
    int optimal_batch_size = static_cast<int>(calculated_batch);

    // Bound: Cannot be larger than the actual total number of orbital pairs
    if (optimal_batch_size > total_ov) {
        optimal_batch_size = total_ov;
    }

    // Bound: Must be at least 1 to prevent division by zero or empty loops on tiny-memory profiles
    if (optimal_batch_size < 1) {
        optimal_batch_size = 1;
    }

    // Print a breakdown to the output file so you can audit it
    outfile->Printf("\n  ==> Dynamic Batch Size Allocation <==\n");
    outfile->Printf("  Total Psi4 Memory Allowed: %12.3f MB\n", total_psi4_mem / (1024.0 * 1024.0));
    outfile->Printf("  Available Batch Memory:    %12.3f MB\n", available_mem / (1024.0 * 1024.0));
    outfile->Printf("  Total (ia) Pairs to Process: %d\n", total_ov);
    outfile->Printf("  Calculated Optimal Batch Size: %d\n\n", optimal_batch_size);

    return optimal_batch_size;
}

// batched version of RPA
//
void PolaritonicRDFRPA::common_init_batched(std::shared_ptr<Wavefunction> dummy_wfn) {

    outfile->Printf("\n");
    outfile->Printf("\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    Polaritonic Restricted DFRPA                     *\n");
    outfile->Printf( "        *    Batched Version                                  *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf("\n");


    // ensure scf_type df
    if ( options_.get_str("SCF_TYPE") != "DF" ) {
        throw PsiException("batched polaritonic RPA only works with scf_type df for now",__FILE__,__LINE__);
    }

    // set up basis set and integrals (should be replaced with DFHelper)
    o_ = nalpha_;
    v_ = nso_ - nalpha_;
    std::shared_ptr<BasisSet> primary = reference_wavefunction_->get_basisset("ORBITAL");
    std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");
    nQ_ = auxiliary->nbf();

//    siap_ = (double*)malloc(o_*v_*nQ_*sizeof(double));
//    memset((void*)siap_,'\0',o_*v_*nQ_*sizeof(double));
//
//    std::shared_ptr<DFTensor> DF (new DFTensor(primary,auxiliary,Ca_,o_,v_,o_,v_,options_));
//    std::shared_ptr<Matrix> tmpov = DF->Qov();
//    double ** Qov = tmpov->pointer();
//
//    for (size_t Q = 0; Q < nQ_; Q++) {
//        for (size_t i = 0; i < o_; i++) {
//            for (size_t a = 0; a < v_; a++) {
//                size_t ai = i * v_ + a;
//                siap_[Q * o_ * v_ + i * v_ + a] = *&Qov[Q][ai];
//            }
//        }
//    }
//
    // initialize dfhelper instead of manual construction of 3-index integrals
    df_helper_ = std::make_shared<DFHelper>(primary, auxiliary); 
    df_helper_->initialize();

    // copy coefficients from Ca_ to occ and virt matrices for df_helper
    auto Cocc = std::make_shared<Matrix>(nso_, o_);
    auto Cvirt = std::make_shared<Matrix>(nso_, v_);

    for (int mu = 0; mu < nso_; mu++) {
        for (int i = 0; i < o_; i++) {
            Cocc->pointer()[mu][i] = Ca_->pointer()[mu][i];
        }
        for (int a = 0; a < v_; a++) {
            Cvirt->pointer()[mu][a] = Ca_->pointer()[mu][a + o_];
        }
    }

    // register spaces in DFHelper
    df_helper_->add_space("occ", Cocc);
    df_helper_->add_space("virt", Cvirt);

    // define 3-index transformation
    df_helper_->add_transformation("Qov", "occ", "virt", "Qpq");

    // transform 3-index integrals using DFHelper
    df_helper_->transform();

    //
    // initialize diagonal matrix for quadrature based on orbital energy
    // differences
    diag_ = (double*)malloc(o_*v_*sizeof(double));
    for (size_t i = 0; i < o_; i++) {
        for (size_t a = 0; a < v_; a++) {
            size_t ai = i * v_ + a;
            diag_[ai] = reference_wavefunction_->epsilon_a()->pointer()[a + o_] - reference_wavefunction_->epsilon_a()->pointer()[i];
        }
    }

}

double PolaritonicRDFRPA::compute_energy_batched() {

    common_init_batched(reference_wavefunction_);

    outfile->Printf("\n");
    outfile->Printf("    No. basis functions:            %5i\n",nso_);
    outfile->Printf("    No. electrons:                  %5i\n",nalpha_ + nbeta_);
    if ( options_.get_bool("QED_USE_RELAXED_ORBITALS") ) {
        outfile->Printf("    Relaxed orbitals are used. K terms included in g\n");
    }
    if ( options_.get_bool("RPA_EXCHANGE") ) {
        outfile->Printf("    Exchange added to A and B in RPA kernel\n");
    }

    if ( n_photon_states_ > 1 ) {
        update_cavity_terms();
    }

    // transform dipole integrals to MO basis
    dipole_[0]->transform(Ca_);
    dipole_[1]->transform(Ca_);
    dipole_[2]->transform(Ca_);

    int o = nalpha_;
    int v = nso_ - nalpha_;
    int ov_total = o * v;
    // set batchsize manually
//    int ia_batch_size = 100; // size of batch for processing (ia|ia) elements, can be tuned for performance
//  use function to determine batchsize based on available memory and system
//  size
    int ia_batch_size = determine_optimal_batch_size();

    int nQ_eff = nQ_ + 1; // size of effective Q matrix including photon contribution
    

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


    // build diagonal matrix for quadrature
    std::shared_ptr<Vector> diag_matrix = build_diag_matrix();

    // build grid points and weights for quadrature
    // Curtis–Clenshaw quadrature based on diagonal approximation to RPA energy expression
    auto grid = build_grid();
    // alternative: minimax quadrature based on diagonal approximation to RPA energy expression
    // needs debugging! 
//    double tol = 1e-8;
//   auto grid = build_minimax_grid(tol);

    // print out grid points and weights for debugging
    if (n_output_ > 2)  {
        outfile->Printf("     =====> orbital energy differences <=====\n");
        for (size_t i = 0; i < o_ * v_ ; i++) {
            outfile->Printf("%20.12lf \n",diag_[i]);
        }
        outfile->Printf("     =====> (ia|ia) elements <=====\n");
        for (size_t i = 0; i < o_ * v_ ; i++) {
            outfile->Printf("%20.12lf \n",diag_matrix->pointer()[i]);
        }
        outfile->Printf("     =====> quadrature weights <=====\n");
        for (size_t i = 0; i < n_grid_points_ ; i++) {
            outfile->Printf("%20.12lf \n",grid->pointer()[i][1]);
        }
        outfile->Printf("     =====> quadrature positions <=====\n");
        for (size_t i = 0; i < n_grid_points_ ; i++) {
            outfile->Printf("%20.12lf \n",grid->pointer()[i][0]);
        }
    }

    // Static DSE contribution
    std::vector<double> mu_ai(o_*v_);

    for (size_t i = 0; i < o_; i++) {
        for (size_t a = 0; a < v_; a++) {
            size_t ai = i * v_ + a;
            mu_ai[ai] = dz[i][a+o_];
        }
    } 

    // prefactor for DSE contribution
    double omega_c = cavity_frequency_[2];
    double alpha_DSE = lambda_z * lambda_z / (2.0 * omega_c * omega_c);

    // build synthetic DSE matrix in auxiliary basis
    std::vector<double> S_DSE(o_*v_);
    for (size_t ai = 0; ai < o_*v_; ++ai)
      // factor 0.5 cancel later scaling of 3-index integrals with factor 2
      //S_DSE[ai] = std::sqrt(0.5 * alpha_DSE) * mu_ai[ai];
      S_DSE[ai] = lambda_z * mu_ai[ai];


    // fetch the underlying MO-transformed fitting coefficients (3-index
    // integrals) from DFHelper
    auto full_Qov = df_helper_->get_tensor("Qov");
    double** Qov_pointer = full_Qov->pointer();

    // compute Qmat_ for each grid point and accumulate Ecrpa
    double ecrpa = 0.0;
    // loop over grid points
    for (int i = 0; i < n_grid_points_; ++i) {

      double freq = grid->pointer()[i][0];

      auto tmp2 = std::make_shared<Matrix>(nQ_eff, nQ_eff);
      tmp2->zero();
      // compute contribution from 3-index integrals to Qmat_ for current grid point

      // Batch over ia pairs
      for (int ia_start = 0; ia_start < ov_total; ia_start += ia_batch_size) {
          int ia_end = std::min(ia_start + ia_batch_size, ov_total);
          int batch_size = ia_end - ia_start;

          //outfile->Printf("     =====> ov_total <=====\n");
          //outfile->Printf("%10i \n",ov_total);
          if (n_output_ > 2) {
            outfile->Printf("     =====> batchsize <=====\n");
            outfile->Printf("%10i, %4i, %4i \n",batch_size, i, ia_start);
          }

          // obtain 3-index integrals for current batch 
          auto Qov_batch_scaled = std::make_shared<Matrix>(nQ_eff, batch_size);

          //populate Qov_batch_scaled with 3-index integrals from DFHelper and
          //add DSE contribution, applying frequency-dependent scaling on the
          //fly to avoid storing large unscaled 3-index integrals
          for (int q = 0; q < nQ_eff; ++q) {
              for (int idx = 0; idx < batch_size; ++idx) {
                  int ai = ia_start + idx;
                  double raw_val = 0.0;
                  
                  if (q < nQ_) {
                    raw_val = Qov_pointer[q][ai];
                  }  else {
                    // photon coupling contribution to 3-index integral
                    raw_val = S_DSE[ai];
                  }


                  // apply frequency-dependent scaling to 3-index integrals for current grid point
                  double scal = sqrt(diag_[ai] / (diag_[ai] * diag_[ai] + freq * freq));
                  Qov_batch_scaled->pointer()[q][idx] = raw_val * scal;
              }
          }

          // Accumulate contribution to Qmat_ from current batch using scaled 3-index integrals
          // tmp2 += Qov_batch_scaled * Qov_batch_scaled^T * 4.0; // factor of 4 for spin summation in closed-shell case
          tmp2 ->gemm(false, true, 4.0, Qov_batch_scaled, Qov_batch_scaled, 1.0);

          if (n_output_ > 2) {
            outfile->Printf("   =====> Qov_batch for gridpoint %i <=====\n", i);
            Qov_batch_scaled->print();
          }
          
      }

      //outfile->Printf("   =====> Qmat_ for gridpoint %i <=====\n", i);
      //tmp2->print();
      // compute trace of Qmat_ (unscreened)
      double trace_Q = 0.0;
      for (int q = 0; q < nQ_eff; ++q) {
          trace_Q += tmp2->pointer()[q][q];
      }
      // -Tr(Q)
      double ecrpa_pnt = 0.0;
      ecrpa_pnt = -trace_Q;

      // add QED-kernel contribution to Q: static DSE term + photon coupling term
      double denom = std::sqrt(cavity_frequency_[2]*cavity_frequency_[2] + grid->pointer()[i][0]*grid->pointer()[i][0]);
      double Sph = freq / denom;

      // scale last row and last column by Sph to include photon coupling contribution to Q
      for (int q = 0; q < nQ_eff; ++q) {
        tmp2->pointer()[q][nQ_] *= Sph;   // column
        tmp2->pointer()[nQ_][q] *= Sph;   // row
      }
        
      // accumulate the correlation energy
      auto evecs = std::make_shared<Matrix>(nQ_eff, nQ_eff);
      auto evals = std::make_shared<Vector>(nQ_eff);

      // compute RPA correlation energy
      // form Q_eff + 1
      for (int q = 0; q < nQ_eff; ++q) {
          //ecrpa_pnt -= tmp2->pointer()[Q][Q];
          tmp2->pointer()[q][q] += 1.0;
      }
//        tmp2->print();
      // diagonalize Qmat_ + 1
      evals->zero();
      evecs->zero();
      tmp2->diagonalize(evecs,evals);
      //evals->print();

      // accumulate log of eigenvalues for RPA energy
      double ecln = 0.0;
      for (int q = 0; q < nQ_eff; ++q) {
          ecln +=  log(evals->pointer()[q]);
      }
      // Tr(ln(Q_eff+1))
      ecrpa_pnt +=  ecln;
      // scale with quadrature weight and prefactors
      ecrpa_pnt = ecrpa_pnt/(4.0 * M_PI) * grid->pointer()[i][1];
      // accumulate contribution to total RPA correlation energy
      ecrpa += ecrpa_pnt;

    }
    outfile->Printf("     =====> RPA correlation energy <=====\n");
    //printf("%20.12lf\n", ecrpa);
    outfile->Printf("%20.12lf \n",ecrpa);

    return ecrpa;

}

} // End namespace hilbert

