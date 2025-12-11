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
    free(tmp1_);
    free(Qmat_);


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
        std::shared_ptr<Matrix> tmp1 = (std::shared_ptr<Matrix>)(new Matrix(o_*v_,o_*v_));
        std::shared_ptr<Matrix> tmp2 = (std::shared_ptr<Matrix>)(new Matrix(nQ_,nQ_));
        tmp1->gemm(true,false,1.0,tmpov,tmpov,0.0);
        tmp2->gemm(false,true,1.0,tmpov,tmpov,0.0);
        tmp1->print();
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

    std::shared_ptr<Vector> diag_matrix = build_diag_matrix();
    std::shared_ptr<Matrix> grid = build_grid();
    if (n_output_ > 1)  {
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

    tmp1_ = (double*)malloc(nQ_*o_*v_*sizeof(double));
    memset((void*)tmp1_,'\0',nQ_*o_*v_*sizeof(double));
    int * int1_ = (int*)malloc(nQ_*sizeof(int));
    memset((void*)int1_,'\0',nQ_*sizeof(int));
    Qmat_ = (double*)malloc(nQ_*nQ_*sizeof(double));
    memset((void*)Qmat_,'\0',nQ_*nQ_*sizeof(double));
    std::shared_ptr<Matrix> tstcopy = (std::shared_ptr<Matrix>)(new Matrix(nQ_,o_*v_));

    // compute Qmat_ for each grid point and accumulate Ecrpa
    double ecrpa = 0.0;
    // loop over grid points
    for (int i = 0; i < n_grid_points_; ++i) {

        //std::shared_ptr<Matrix> tmp1 = (std::shared_ptr<Matrix>)(new Matrix(o_*v_,o_*v_));
        std::shared_ptr<Matrix> tmp2 = (std::shared_ptr<Matrix>)(new Matrix(nQ_,nQ_));
        //tmp1->gemm(true,false,1.0,tst_,tst_,0.0);
        //tmp2->gemm(false,true,1.0,tst_,tst_,0.0);
        //tmp1->print();
        //tmp2->print();
        //tst_->print();
        tmp2->zero();
        // copy 3-index integrals to local variable tstcopy
        tstcopy->zero();
        tstcopy->copy(tst_);
        // scale 3 index integrals
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
            outfile->Printf("   =====> tstcopy after scaling <=====\n");
            tstcopy->print();
        }

        // accumulate contribution on Qmat_
        // scaling factor for 3 index integrals
        double scal2 = 2.0*2.0;
        tmp2->gemm(false,true,scal2,tstcopy,tstcopy,1.0);

        if (n_output_ > 2) {
            outfile->Printf("   =====> Qmat_ for gridpoint %i <=====\n", i);
            tmp2->print();
        }

        std::shared_ptr<Matrix> evecs = (std::shared_ptr<Matrix>)(new Matrix(nQ_,nQ_));
        std::shared_ptr<Vector> evals = (std::shared_ptr<Vector>) (new Vector(nQ_));
        tmp2->diagonalize(evecs,evals);
        //evals->print();
        double erpa = 0.0;
        for (size_t i = 0; i < nQ_; i++) {
            erpa += evals->pointer()[i];
        }
        outfile->Printf("    =====> tr(Qmat) for gridpoint %i <=====\n", i);
        outfile->Printf("%20.12lf \n",erpa);
        // compute RPA correlation energy
        double ecrpa_pnt = 0.0;
        for (size_t Q = 0; Q < nQ_; Q++) {
            ecrpa_pnt -= tmp2->pointer()[Q][Q];
            tmp2->pointer()[Q][Q] += 1.0;
        }
//        tmp2->print();
        // diagonalize Qmat_ + 1
        evals->zero();
        evecs->zero();
        tmp2->diagonalize(evecs,evals);
        //evals->print();

        double ecln = 0.0;
        double epnt_copy = ecrpa_pnt;
        for (size_t Q = 0; Q < nQ_; Q++) {
            ecln +=  log(evals->pointer()[Q]);
        }
        ecrpa_pnt +=  ecln;
        ecrpa_pnt = ecrpa_pnt/(4*M_PI) * grid->pointer()[i][1];
        ecrpa += ecrpa_pnt;

    }
    outfile->Printf("     =====> RPA correlation energy <=====\n");
    //printf("%20.12lf\n", ecrpa);
    outfile->Printf("%20.12lf \n",ecrpa);

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

// function to compute RPA correlation energy
//std::shared_ptr<Matrix> PolaritonicRDFRPA::rpa_ecorr() {


} // End namespaces

