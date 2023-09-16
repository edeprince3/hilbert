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
#include "psi4/libscf_solver/hf.h"

#include <psi4/psifiles.h>
#include <psi4/libtrans/integraltransform.h>

#include "rtddft.h"
#include "rrpa.h"

#include <misc/blas.h>
#include <misc/hilbert_psifiles.h>
#include <misc/omp.h>
#include <misc/threeindexintegrals.h>
#include <misc/nonsym_davidson_solver.h>

using namespace psi;
using namespace fnocc;

namespace hilbert{ 

PolaritonicRRPA::PolaritonicRRPA(std::shared_ptr<Wavefunction> reference_wavefunction, Options& options_, std::shared_ptr<Wavefunction> dummy_wfn):
    PolaritonicHF(reference_wavefunction,options_) {
    common_init(dummy_wfn);
}

PolaritonicRRPA::~PolaritonicRRPA() {

    free(int1_);
    free(int2_);

}

void PolaritonicRRPA::common_init(std::shared_ptr<Wavefunction> dummy_wfn) {

    outfile->Printf("\n");
    outfile->Printf("\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    Polaritonic Restricted RPA                       *\n");
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

    o_ = nalpha_;
    v_ = nso_ - nalpha_;

    nQ_ = 0;
    if ( options_.get_str("SCF_TYPE") == "DF" ) {

        // get auxiliary basis:
        std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");

        // total number of auxiliary basis functions
        nQ_ = auxiliary->nbf();

        std::shared_ptr<DFTensor> DF (new DFTensor(primary,auxiliary,Ca_,o_,v_,o_,v_,options_));

        std::shared_ptr<Matrix> tmpoo = DF->Qoo();
        std::shared_ptr<Matrix> tmpov = DF->Qov();
        std::shared_ptr<Matrix> tmpvv = DF->Qvv();

        double ** Qoo = tmpoo->pointer();
        double ** Qov = tmpov->pointer();
        double ** Qvv = tmpvv->pointer();

        int1_ = (double*)malloc(o_*o_*v_*v_*sizeof(double));
        int2_ = (double*)malloc(o_*o_*v_*v_*sizeof(double));

        memset((void*)int1_,'\0',o_*o_*v_*v_*sizeof(double));
        memset((void*)int2_,'\0',o_*o_*v_*v_*sizeof(double));

        F_DGEMM('n','t',o_*v_,o_*v_,nQ_,1.0,&(Qov[0][0]),o_*v_,&(Qov[0][0]),o_*v_,0.0,int1_,o_*v_);
        F_DGEMM('n','t',v_*v_,o_*o_,nQ_,1.0,&(Qvv[0][0]),v_*v_,&(Qoo[0][0]),o_*o_,0.0,int2_,v_*v_);

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

        int1_ = (double*)malloc(o_*o_*v_*v_*sizeof(double));
        int2_ = (double*)malloc(o_*o_*v_*v_*sizeof(double));

        memset((void*)int1_,'\0',o_*o_*v_*v_*sizeof(double));
        memset((void*)int2_,'\0',o_*o_*v_*v_*sizeof(double));

        F_DGEMM('n','t',o_*v_,o_*v_,nQ_,1.0,Qov,o_*v_,Qov,o_*v_,0.0,int1_,o_*v_);
        F_DGEMM('n','t',v_*v_,o_*o_,nQ_,1.0,Qvv,v_*v_,Qoo,o_*o_,0.0,int2_,v_*v_);

        free(Qoo);
        free(Qov);
        free(Qvv);
    }

    // orbital energies
    epsilon_a_= SharedVector(new Vector(nirrep_, nmopi_));
    epsilon_a_->copy(reference_wavefunction_->epsilon_a().get());
    epsilon_b_= SharedVector(new Vector(nirrep_, nmopi_));
    epsilon_b_->copy(reference_wavefunction_->epsilon_b().get());

}

double PolaritonicRRPA::compute_energy() {

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

    std::shared_ptr<Matrix> H = build_rpa_matrix();
    std::shared_ptr<Matrix> eigvec (new Matrix(2 * o_ * v_ + 2, 2 * o_ * v_ + 2));
    std::shared_ptr<Vector> eigval (new Vector(2 * o_ * v_ + 2));

    H->diagonalize(eigvec, eigval);
    eigval->print();

    return 0.0;
}

// cQED-RPA:
// 
// |  A  B  g*  g* | ( X )     |  1  0  0  0 | ( X )
// |  B  A  g*  g* | ( Y ) = W |  0 -1  0  0 | ( Y )
// |  g  g  w   0  | ( M )     |  0  0  1  0 | ( M )
// |  g  g  0   w  | ( N )     |  0  0  0 -1 | ( N )
// 
// 
// 
// this function solves the generalized eigenvalue problem S c = 1/w H c
std::shared_ptr<Matrix> PolaritonicRRPA::build_rpa_matrix() {

    std::shared_ptr<Matrix> H (new Matrix(2 * o_ * v_ + 2, 2 * o_ * v_ + 2));

    for (size_t a = 0; a < v_; a++) {
        for (size_t i = 0; i < o_; i++) {
            size_t ai = a * o_ + i;
            for (size_t b = 0; b < v_; b++) {
                for (size_t j = 0; j < o_; j++) {
                    size_t bj = b * o_ + j;
                    double A_aibj = (i == j) * (a == b) * epsilon_a_->pointer()[a + o_]
                                  - (i == j) * (a == b) * epsilon_a_->pointer()[i]
                                  + 2.0 * int1_[i * o_ * v_ * v_ + a * o_ * v_ + j * v_ + b]
                                  - int2_[a * o_ * o_ * v_ + b * o_ * o_ + i * o_ + j];
                    H->pointer()[ai][bj] = A_aibj;
                }
            }
        }
    }
    

    return H;

}


} // End namespaces

