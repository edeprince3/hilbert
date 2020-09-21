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

#include "rcis.h"

#include <misc/blas.h>

using namespace psi;
using namespace fnocc;

namespace hilbert{ 

PolaritonicRCIS::PolaritonicRCIS(std::shared_ptr<Wavefunction> reference_wavefunction, Options& options_):
    PolaritonicHF(reference_wavefunction,options_) {
    common_init();
}

PolaritonicRCIS::~PolaritonicRCIS() {

    free(int1_);
    free(int2_);

}

void PolaritonicRCIS::common_init() {

    outfile->Printf("\n\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    Polaritonic RCIS                                 *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *******************************************************\n");

    // ensure scf_type df
    if ( options_.get_str("SCF_TYPE") != "DF" ) {
        throw PsiException("polaritonic rcis only works with scf_type df for now",__FILE__,__LINE__);
    }

    // ensure running in c1 symmetry
    if ( reference_wavefunction_->nirrep() > 1 ) {
        throw PsiException("polaritonic rcis only works with c1 symmetry for now.",__FILE__,__LINE__);
    }

    // ensure closed shell
    if ( nalpha_ != nbeta_ ) {
        throw PsiException("polaritonic rcis only works for closed shells",__FILE__,__LINE__);
    }

    // get primary basis:
    std::shared_ptr<BasisSet> primary = reference_wavefunction_->get_basisset("ORBITAL");

    // get auxiliary basis:
    std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");

    // total number of auxiliary basis functions
    int nQ = auxiliary->nbf();

    int o = nalpha_;
    int v = nso_ - nalpha_;
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

    F_DGEMM('n','t',o*v,o*v,nQ,1.0,&(Qov[0][0]),o*v,&(Qov[0][0]),o*v,0.0,int1_,o*v);
    F_DGEMM('n','t',v*v,o*o,nQ,1.0,&(Qvv[0][0]),v*v,&(Qoo[0][0]),o*o,0.0,int2_,v*v);

}

double PolaritonicRCIS::compute_energy() {

    // get auxiliary basis:
    std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");

    // total number of auxiliary basis functions
    int nQ = auxiliary->nbf();

    outfile->Printf("\n");
    outfile->Printf("    No. basis functions:            %5i\n",nso_);
    outfile->Printf("    No. auxiliary basis functions:  %5i\n",nQ);
    outfile->Printf("    No. electrons:                  %5i\n",nalpha_ + nbeta_);

    // print cavity properties
    //print_cavity_properties_ = true;
    if ( n_photon_states_ > 1 ) {
        build_cavity_hamiltonian();
    }
    //print_cavity_properties_ = false;

    // transform dipole integrals to MO basis
    dipole_[0]->transform(Ca_);
    dipole_[1]->transform(Ca_);
    dipole_[2]->transform(Ca_);

    int o = nalpha_;
    int v = nso_ - nalpha_;

    std::shared_ptr<Matrix> ham (new Matrix((o*v+1)*n_photon_states_,(o*v+1)*n_photon_states_));

    double ** hp = ham->pointer();
    double * ep = epsilon_a_->pointer();

    double coupling_factor_x = cavity_frequency_[0] * cavity_coupling_strength_[0];
    double coupling_factor_y = cavity_frequency_[1] * cavity_coupling_strength_[1];
    double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];

    double lambda_x = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
    double lambda_y = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

    double ** dx = dipole_[0]->pointer();
    double ** dy = dipole_[1]->pointer();
    double ** dz = dipole_[2]->pointer();

    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            for (int n = 0; n < n_photon_states_; n++) {
                int ian = i * v * n_photon_states_ + a * n_photon_states_ + n;
                for (int j = 0; j < o; j++) {
                    for (int b = 0; b < v; b++) {
                        for (int m = 0; m < n_photon_states_; m++) {
                        int jbm = j * v * n_photon_states_ + b * n_photon_states_ + m;

                            // cases:
                            hp[ian][jbm] = 0.0;

                            // full diagonal
                            if ( i == j && a == b && m == n ) {
                                hp[ian][jbm] += ep[a+o] - ep[i] + HCavity_z_->pointer()[m][m];
                            }

                            // diagonal in photon states
                            if ( m == n ) {
                                int iajb = i * o * v * v + a * o * v + j * v + b;
                                int ijab = i * o * v * v + j * v * v + a * v + b;
                                hp[ian][jbm] += 2.0 * int1_[iajb] - int2_[ijab];

                                //if ( i != j || a != b) {
                                    // dipole self energy contribution?
                                    //hp[ian][jbm] += 0.5 * lambda_z * lambda_z * (2.0 * dz[i][a+o]*dz[j][b+o] - dz[i][j]*dz[a+o][b+o]);
                                //}
                            }

                            if ( i == j && a != b && n == m + 1 ) {
                                hp[ian][jbm] -= sqrt(2.0) * coupling_factor_z * sqrt(n+1) * dz[a+o][b+o];
                            }

                            if ( i == j && a != b && n == m - 1 ) {
                                hp[ian][jbm] -= sqrt(2.0) * coupling_factor_z * sqrt(n) * dz[a+o][b+o];
                            }

                            if ( i != j && a == b && n == m + 1 ) {
                                hp[ian][jbm] -= sqrt(2.0) * coupling_factor_z * sqrt(n+1) * dz[i][j];
                            }

                            if ( i != j && a == b && n == m - 1 ) {
                                hp[ian][jbm] -= sqrt(2.0) * coupling_factor_z * sqrt(n) * dz[i][j];
                            }

                        }
                    }
                }
            }
        }
    }

    // now, couple |0,n+1> and |0,n-1> to |ia,n>
    int off = o * v * n_photon_states_;
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            for (int n = 0; n < n_photon_states_; n++) {
                int ian = i * v * n_photon_states_ + a * n_photon_states_ + n;
                // case 1
                if ( n < n_photon_states_ - 1 ) {
                    hp[ian][off+n+1] = -sqrt(2.0) * coupling_factor_z * sqrt(n+1) * dz[i][a+o];
                    hp[off+n+1][ian] = -sqrt(2.0) * coupling_factor_z * sqrt(n+1) * dz[i][a+o];
                }

                // case 2
                if ( n > 0 ) {
                    hp[ian][off+n-1] = -sqrt(2.0) * coupling_factor_z * sqrt(n) * dz[i][a+o];
                    hp[off+n-1][ian] = -sqrt(2.0) * coupling_factor_z * sqrt(n) * dz[i][a+o];
                }
            }
        }
    }

    // now, couple |0,n> to all other |0,m>
    for (int n = 0; n < n_photon_states_; n++) {
        for (int m = 0; m < n_photon_states_; m++) {
            hp[off+n][off+m] += HCavity_z_->pointer()[n][m];
        }
    }

    std::shared_ptr<Matrix> eigvec (new Matrix((o*v+1)*n_photon_states_,(o*v+1)*n_photon_states_));
    std::shared_ptr<Vector> eigval (new Vector((o*v+1)*n_photon_states_));
    ham->diagonalize(eigvec,eigval);

    eigval->print();

/*
    //for (int i = 0; i < (o*v+1)*n_photon_states_; i++) {
    for (int i = 0; i < 50; i++) {
        double photon_weight = 0.0;
        for (int j = off + 1; j < off + n_photon_states_; j++) {
            double dum = eigvec->pointer()[j][i];
            photon_weight += dum*dum;
        }
        printf(" %20.12lf %20.12lf",energy_ + eigval->pointer()[i],photon_weight);
    }
    printf("\n");
    fflush(stdout);
*/
    
    // print orbital energies
    //epsilon_a_->print();

    return 0.0;
}

} // End namespaces
