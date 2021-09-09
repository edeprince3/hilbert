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

#include <psi4/libmints/mintshelper.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/molecule.h>
#include <psi4/libmints/factory.h>
#include <psi4/libmints/integral.h>
#include <psi4/libmints/basisset.h>

#include <psi4/libqt/qt.h>

#include <psi4/libpsi4util/process.h>

#include "hf.h"

using namespace psi;

namespace hilbert{ 

PolaritonicHF::PolaritonicHF(std::shared_ptr<Wavefunction> reference_wavefunction, Options& options):
    Wavefunction(options) {
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

PolaritonicHF::~PolaritonicHF() {
}

void PolaritonicHF::common_init() {

    // from reference:

    shallow_copy(reference_wavefunction_);

    energy_   = reference_wavefunction_->energy();
    nalpha_   = reference_wavefunction_->nalpha();
    nbeta_    = reference_wavefunction_->nbeta();
    nalphapi_ = reference_wavefunction_->nalphapi();
    nbetapi_  = reference_wavefunction_->nbetapi();
    doccpi_   = reference_wavefunction_->doccpi();
    soccpi_   = reference_wavefunction_->soccpi();
    frzcpi_   = reference_wavefunction_->frzcpi();
    frzvpi_   = reference_wavefunction_->frzvpi();
    nmopi_    = reference_wavefunction_->nmopi();
    nirrep_   = reference_wavefunction_->nirrep();
    nso_      = reference_wavefunction_->nso();
    nmo_      = reference_wavefunction_->nmo();
    nsopi_    = reference_wavefunction_->nsopi();
    molecule_ = reference_wavefunction_->molecule();
    enuc_     = molecule_->nuclear_repulsion_energy({0.0,0.0,0.0});

    AO2SO_ = std::shared_ptr<Matrix>(reference_wavefunction_->aotoso());

    Ca_ = std::shared_ptr<Matrix>(reference_wavefunction_->Ca());
    Cb_ = std::shared_ptr<Matrix>(reference_wavefunction_->Cb());

    S_  = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->S()));

    H_  = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->H()));

    Fa_  = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Fa()));
    Fb_  = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Fb()));

    Da_  = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Da()));
    Db_  = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Db()));

    // Lagrangian matrix
    Lagrangian_ = std::shared_ptr<Matrix>(reference_wavefunction_->Lagrangian());

    epsilon_a_= SharedVector(new Vector(nirrep_, nmopi_));
    epsilon_a_->copy(reference_wavefunction_->epsilon_a().get());
    epsilon_b_= SharedVector(new Vector(nirrep_, nmopi_));
    epsilon_b_->copy(reference_wavefunction_->epsilon_b().get());

    // other stuff:
    gradient_     =  reference_wavefunction_->matrix_factory()->create_shared_matrix("Total gradient", molecule_->natom(), 3);
    multiplicity_ = Process::environment.molecule()->multiplicity();

    // memory is from process::environment
    memory_ = Process::environment.get_memory();

    average_electric_dipole_self_energy_ = 0.0;

    // initialize cavity parameters
    initialize_cavity();
}

std::shared_ptr<Matrix> PolaritonicHF::OrbitalGradient(std::shared_ptr<Matrix> D,
                                                       std::shared_ptr<Matrix> F,
                                                       std::shared_ptr<Matrix> Shalf) {

    std::shared_ptr<Matrix> ShalfGradShalf(new Matrix("ST^{-1/2}(FDS - SDF)S^{-1/2}", nso_, nso_));

    std::shared_ptr<Matrix> FDSmSDF(new Matrix("FDS-SDF", nso_, nso_));
    std::shared_ptr<Matrix> DS(new Matrix("DS", nso_, nso_));

    DS->gemm(false,false,1.0,D,S_,0.0);
    FDSmSDF->gemm(false,false,1.0,F,DS,0.0);

    DS.reset();

    std::shared_ptr<Matrix> SDF(FDSmSDF->transpose());
    FDSmSDF->subtract(SDF);

    SDF.reset();

    std::shared_ptr<Matrix> ShalfGrad(new Matrix("ST^{-1/2}(FDS - SDF)", nso_, nso_));
    ShalfGrad->gemm(true,false,1.0,Shalf,FDSmSDF,0.0);
    FDSmSDF.reset();

    ShalfGradShalf->gemm(false,false,1.0,ShalfGrad,Shalf,0.0);

    ShalfGrad.reset();

    return ShalfGradShalf;
}

void PolaritonicHF::initialize_cavity() {

    n_photon_states_ = options_.get_int("N_PHOTON_STATES");

    cavity_frequency_ = (double*)malloc(sizeof(double)*3);
    memset((void*)cavity_frequency_,'\0',3*sizeof(double));
    if (options_["CAVITY_FREQUENCY"].has_changed()){
       if (options_["CAVITY_FREQUENCY"].size() != 3)
          throw PsiException("The CAVITY E array has the wrong dimensions",__FILE__,__LINE__);
       for (int i = 0; i < 3; i++) cavity_frequency_[i] = options_["CAVITY_FREQUENCY"][i].to_double();
    }else{
       cavity_frequency_[0] = 0.0;
       cavity_frequency_[1] = 0.0;
       cavity_frequency_[2] = 2.042/27.21138;
    }

    // Read in the cavity coupling stregth
    cavity_coupling_strength_ = (double*)malloc(sizeof(double)*3);
    memset((void*)cavity_coupling_strength_,'\0',3*sizeof(double));
    if (options_["CAVITY_COUPLING_STRENGTH"].has_changed()){
       if (options_["CAVITY_COUPLING_STRENGTH"].size() != 3)
          throw PsiException("The CAVITY_COUPLING_STRENGTH array has the wrong dimensions",__FILE__,__LINE__);
       for (int i = 0; i < 3; i++) cavity_coupling_strength_[i] = options_["CAVITY_COUPLING_STRENGTH"][i].to_double();
    }else{
       cavity_coupling_strength_[0] = 0.0;
       cavity_coupling_strength_[1] = 0.0;
       cavity_coupling_strength_[2] = 2990.0/2.54175;
    }

    // CCSD currently won't work with any polarization other that z. 
    // throw an exception here to be safe for now. 
    // 
    // TODO: verify whether or not RHF/ROHF/UHF/CIS work correctly with non-z-polarized 
    // modes.  if so, move this exception to the CCSD code.
    // 
    if ( fabs(cavity_coupling_strength_[0]) > 1e-12 ) {
        throw PsiException("cQED codes currently only work with z-polarized modes",__FILE__,__LINE__);
    }
    if ( fabs(cavity_coupling_strength_[1]) > 1e-12 ) {
        throw PsiException("cQED codes currently only work with z-polarized modes",__FILE__,__LINE__);
    }
    if ( fabs(cavity_frequency_[0]) > 1e-12 ) {
        throw PsiException("cQED codes currently only work with z-polarized modes",__FILE__,__LINE__);
    }
    if ( fabs(cavity_frequency_[1]) > 1e-12 ) {
        throw PsiException("cQED codes currently only work with z-polarized modes",__FILE__,__LINE__);
    }

    // get nuclear contribution to the molecular total dipole moment

    nuc_dip_x_ = 0.0;
    nuc_dip_y_ = 0.0;
    nuc_dip_z_ = 0.0;
    for (int i = 0; i < molecule_->natom(); i++) {
        nuc_dip_x_ += molecule_->Z(i) * molecule_->x(i);
        nuc_dip_y_ += molecule_->Z(i) * molecule_->y(i);
        nuc_dip_z_ += molecule_->Z(i) * molecule_->z(i);
    }

    // integrals

    std::shared_ptr<MintsHelper> mints (new MintsHelper(reference_wavefunction_));

    // dipole integrals
    dipole_ = mints->so_dipole();

    // quadrupole integrals
    std::vector< std::shared_ptr<Matrix> > quadrupole = mints->so_quadrupole();

    double lambda_x = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
    double lambda_y = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

    quadrupole[0]->scale(0.5 * lambda_x * lambda_x);
    quadrupole[1]->scale(1.0 * lambda_x * lambda_y); // scaled by 2 for xy + yx
    quadrupole[2]->scale(1.0 * lambda_x * lambda_z); // scaled by 2 for xz + zx
    quadrupole[3]->scale(0.5 * lambda_y * lambda_y);
    quadrupole[4]->scale(1.0 * lambda_y * lambda_z); // scaled by 2 for yz + zy
    quadrupole[5]->scale(0.5 * lambda_z * lambda_z);

    quadrupole_scaled_sum_ = (std::shared_ptr<Matrix>)(new Matrix(nso_,nso_));
    quadrupole_scaled_sum_->zero();
    for (int i = 0; i < 6; i++) {
        quadrupole_scaled_sum_->add(quadrupole[i]);
    }

    // for two-electron part of 1/2 (lambda.d)^2 (1/2 picked up in rhf.cc)

    dipole_scaled_sum_ = (std::shared_ptr<Matrix>)(new Matrix(nso_,nso_));
    dipole_scaled_sum_->zero();
    dipole_scaled_sum_->axpy(lambda_x, dipole_[0]);
    dipole_scaled_sum_->axpy(lambda_y, dipole_[1]);
    dipole_scaled_sum_->axpy(lambda_z, dipole_[2]);

/*
    // e-(n-<d>) contribution 0.5 * 2 (lambda . de) ( lambda . (dn - <d>) )

    // initialize with <de> = 0
    e_dip_x_ = 0.0;
    e_dip_y_ = 0.0;
    e_dip_z_ = 0.0;

    tot_dip_x_ = nuc_dip_x_;
    tot_dip_y_ = nuc_dip_y_;
    tot_dip_z_ = nuc_dip_z_;

    // e contribution: lambda . de
    std::shared_ptr<Matrix> el_dipdot (new Matrix(nso_,nso_));
    el_dipdot->zero();
    el_dipdot->axpy(lambda_x,dipole_[0]);
    el_dipdot->axpy(lambda_y,dipole_[1]);
    el_dipdot->axpy(lambda_z,dipole_[2]);

    // n-<d> contribution: lambda . (dn - <d>)
    double nuc_dipdot = 0.0;
    nuc_dipdot += lambda_x * ( nuc_dip_x_ - tot_dip_x_ );
    nuc_dipdot += lambda_y * ( nuc_dip_y_ - tot_dip_y_ );
    nuc_dipdot += lambda_z * ( nuc_dip_z_ - tot_dip_z_ );

    // e-(n-<d>) contribution 0.5 * 2 (lambda . de) ( lambda . (dn - <d>) )
    scaled_e_n_dipole_squared_ = (std::shared_ptr<Matrix>)(new Matrix(nso_,nso_));
    scaled_e_n_dipole_squared_->copy(el_dipdot);
    scaled_e_n_dipole_squared_->scale(nuc_dipdot);

    // constant terms:  0.5 ( lambda . ( dn - <d> ) )^2 = 0.5 ( lambda . <de> )^2
    average_electric_dipole_self_energy_ = 0.5 * nuc_dipdot * nuc_dipdot;
*/

    update_cavity_terms();

}

/* 
 * update_cavity_terms():
 *
 * update one-electron terms defined in initialize cavity to account 
 * for changing <d>. this is just the de . ( dn - <d> ) term
 *
*/
void PolaritonicHF::update_cavity_terms(){

    int nS = n_photon_states_;

    // evaluate the electronic contribute to the molecule's dipole moment

    e_dip_x_ = C_DDOT(nso_*nso_,&(Da_->pointer())[0][0],1,&(dipole_[0]->pointer())[0][0],1); 
    e_dip_y_ = C_DDOT(nso_*nso_,&(Da_->pointer())[0][0],1,&(dipole_[1]->pointer())[0][0],1); 
    e_dip_z_ = C_DDOT(nso_*nso_,&(Da_->pointer())[0][0],1,&(dipole_[2]->pointer())[0][0],1);

    if ( same_a_b_dens_ ) {
        e_dip_x_ *= 2.0;
        e_dip_y_ *= 2.0;
        e_dip_z_ *= 2.0;
    }else {
        e_dip_x_ += C_DDOT(nso_*nso_,&(Db_->pointer())[0][0],1,&(dipole_[0]->pointer())[0][0],1); 
        e_dip_y_ += C_DDOT(nso_*nso_,&(Db_->pointer())[0][0],1,&(dipole_[1]->pointer())[0][0],1); 
        e_dip_z_ += C_DDOT(nso_*nso_,&(Db_->pointer())[0][0],1,&(dipole_[2]->pointer())[0][0],1);
    }

    // evaluate the total dipole moment:

    tot_dip_x_ = e_dip_x_ + nuc_dip_x_;
    tot_dip_y_ = e_dip_y_ + nuc_dip_y_;
    tot_dip_z_ = e_dip_z_ + nuc_dip_z_;

    // e-(n-<d>) contribution 0.5 * 2 (lambda . de) ( lambda . (dn - <d>) )

    double lambda_x = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
    double lambda_y = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

    // e contribution: lambda . de
    std::shared_ptr<Matrix> el_dipdot (new Matrix(nso_,nso_));
    el_dipdot->zero();
    el_dipdot->axpy(lambda_x,dipole_[0]);
    el_dipdot->axpy(lambda_y,dipole_[1]);
    el_dipdot->axpy(lambda_z,dipole_[2]);

    // n contribution: lambda . dn
    // 
    // or ... in the coherent-state basis:
    // 
    // n-<d> contribution: lambda . (dn - <d>)
    double nuc_dipdot = 0.0;
    if ( use_coherent_state_basis_ ) {
        nuc_dipdot += lambda_x * ( nuc_dip_x_ - tot_dip_x_ );
        nuc_dipdot += lambda_y * ( nuc_dip_y_ - tot_dip_y_ );
        nuc_dipdot += lambda_z * ( nuc_dip_z_ - tot_dip_z_ );
    }else {
        nuc_dipdot += lambda_x * nuc_dip_x_;
        nuc_dipdot += lambda_y * nuc_dip_y_;
        nuc_dipdot += lambda_z * nuc_dip_z_;
    }

    // e-(n-<d>) contribution 0.5 * 2 (lambda . de) ( lambda . (dn - <d>) )
    scaled_e_n_dipole_squared_ = (std::shared_ptr<Matrix>)(new Matrix(nso_,nso_));
    scaled_e_n_dipole_squared_->copy(el_dipdot);
    scaled_e_n_dipole_squared_->scale(nuc_dipdot);

    // constant terms:  0.5 ( lambda . dn ^2
    // 
    // or ... in the coherent state basis:
    // 
    // constant terms:  0.5 ( lambda . ( dn - <d> ) )^2
    average_electric_dipole_self_energy_ = 0.5 * nuc_dipdot * nuc_dipdot;

}

} // End namespaces
