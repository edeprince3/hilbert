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

    siap_ = (double*)malloc(o_*v_*nQ_*sizeof(double));
    memset((void*)siap_,'\0',o_*v_*nQ_*sizeof(double));

    std::shared_ptr<DFTensor> DF (new DFTensor(primary,auxiliary,Ca_,o_,v_,o_,v_,options_));
    std::shared_ptr<Matrix> tmpov = DF->Qov();
    double ** Qov = tmpov->pointer();

    for (size_t Q = 0; Q < nQ_; Q++) {
        for (size_t i = 0; i < o_; i++) {
            for (size_t a = 0; a < v_; a++) {
                size_t ai = i * v_ + a;
                siap_[Q * o_ * v_ + i * v_ + a] = *&Qov[Q][ai];
            }
        }
    }

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
    int ia_batch_size = 1000; // size of batch for processing (ia|ia) elements, can be tuned for performance
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

          // to be updated: should use DFHelper. 
          // compute 3-index integrals for current batch
          //auto tstcopy = compute_3index_integrals_batch(ia_start, batch_size);

          // obtain 3-index integrals for current batch from siap_
          auto Qov_batch = std::make_shared<Matrix>(nQ_eff, batch_size);
          auto Qov_batch_scaled = std::make_shared<Matrix>(nQ_eff, batch_size);

          //populate Qov_batch with 3-index integrals from siap_
          for (int q = 0; q < nQ_eff; ++q) {
              for (int idx = 0; idx < batch_size; ++idx) {
                  int ai = ia_start + idx;
                  double raw_val = 0.0;
                  if (q < nQ_) {
                    // get 3-index integral from siap_
                    raw_val = siap_[q * ov_total + ai];
                  }  else {
                    // photon coupling contribution to 3-index integral
                    raw_val = S_DSE[ai];
                  }

                  Qov_batch->pointer()[q][idx] = raw_val;

                  // apply frequency-dependent scaling to 3-index integrals for current grid point
                  double scal = sqrt(diag_[ai] / (diag_[ai] * diag_[ai] + freq * freq));
                  Qov_batch_scaled->pointer()[q][idx] = raw_val * scal;
              }
          }

          // Accumulate contribution to Qmat_ from current batch using scaled 3-index integrals
          // tmp2 += Qov_batch * Qov_batch_scaled^T * 4.0; // factor of 4 for spin summation in closed-shell case
          tmp2 ->gemm(false, true, 4.0, Qov_batch, Qov_batch_scaled, 1.0);
          
      }

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

// end namespace hilbert
}
