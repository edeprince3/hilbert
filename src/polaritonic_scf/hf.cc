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

    Ca_ = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Ca()));
    Cb_ = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Cb()));

    S_  = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->S()));

    H_  = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->H()));

    Fa_  = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Fa()));
    Fb_  = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Fb()));

    Da_  = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Da()));
    Db_  = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Db()));

    // Lagrangian matrix
    Lagrangian_ = std::shared_ptr<Matrix>(reference_wavefunction_->lagrangian());

    epsilon_a_ = std::make_shared<Vector>(nmopi_);
    epsilon_a_->copy(*reference_wavefunction_->epsilon_a().get());
    epsilon_b_ = std::make_shared<Vector>(nmopi_);
    epsilon_b_->copy(*reference_wavefunction_->epsilon_b().get());

    // other stuff:
    gradient_     =  reference_wavefunction_->matrix_factory()->create_shared_matrix("Total gradient", molecule_->natom(), 3);
    multiplicity_ = Process::environment.molecule()->multiplicity();

    // memory is from process::environment
    memory_ = Process::environment.get_memory();

    // use coherent-state basis? 
    use_coherent_state_basis_ = options_.get_bool("USE_COHERENT_STATE_BASIS");

    average_electric_dipole_self_energy_ = 0.0;

    // initialize cavity parameters
    initialize_cavity();
}

void PolaritonicHF::update_charge(int charge) {
    // update electron count for anion

    molecule_->set_molecular_charge(molecule_->molecular_charge() + charge);

    // update spin count (keep nalpha_ >= nbeta_)
    nalpha_ -= charge;
    if (nalpha_ < nbeta_)
        std::swap(nalpha_, nbeta_);

    same_a_b_dens_ = nalpha_ == nbeta_;

    // calculate new multiplicity
    multiplicity_ = abs(nalpha_ - nbeta_) + 1;

    molecule_->set_multiplicity(multiplicity_);

    // update geometry
    molecule_->update_geometry();

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

    size_t x = 0, y = 1, z = 2;
    size_t xx = 0, xy = 1, xz = 2, yy = 3, yz = 4, zz = 5;

    // quadrupole integrals
    if ( options_.get_bool("USE_QUADRUPOLE_INTEGRALS") ) {
        quadrupole_ = mints->so_quadrupole();
    } else {

        // do not use integrals; compute as mu_pr * mu_rq
        auto square_dipole = [](const std::shared_ptr<Matrix>& d1, const std::shared_ptr<Matrix>& d2) {
            std::shared_ptr<Matrix> result(new Matrix(d1->nrow(),d2->ncol()));
            result->gemm(false,false,1.0,d1,d2,1.0);
            return result;
        };

        quadrupole_ = std::vector<std::shared_ptr<Matrix>>(6);

        // 0: xx
        quadrupole_[xx] = square_dipole(dipole_[x],dipole_[x]);
        // 1: xy
        quadrupole_[xy] = square_dipole(dipole_[x],dipole_[y]);
        // 2: xz
        quadrupole_[xz] = square_dipole(dipole_[x],dipole_[z]);
        // 3: yy
        quadrupole_[yy] = square_dipole(dipole_[y],dipole_[y]);
        // 4: yz
        quadrupole_[yz] = square_dipole(dipole_[y],dipole_[z]);
        // 5: zz
        quadrupole_[zz] = square_dipole(dipole_[z],dipole_[z]);
    }



    // reorder dipole / quadrupole integrals integrals
    // if polarization other than "z" is desired
    //
    axis_order_ = options_.get_str("ROTATE_POLARIZATION_AXIS");

    // copy dipoles
    std::vector< std::shared_ptr<Matrix> > old_d(3);
    for (int i = 0; i < 3; i++) 
        old_d[i] = std::make_shared<Matrix>(dipole_[i]->clone());
    
    // copy quadrupoles
    std::vector< std::shared_ptr<Matrix> > old_q(6);
    for (int i = 0; i < 6; i++) 
        old_q[i] = std::make_shared<Matrix>(quadrupole_[i]->clone()); 
    
    double old_nuc_dip_x = nuc_dip_x_;
    double old_nuc_dip_y = nuc_dip_y_;
    double old_nuc_dip_z = nuc_dip_z_;

    auto swap_quadrupoles = [this, &old_q, &old_d](size_t i, size_t j, bool transpose) {
        quadrupole_[i]->copy(old_q[j]);
        if ( transpose ) {
            quadrupole_[i]->transpose_this();
        }
    };

    if ( axis_order_ == "XYZ"){
        // do nothing
    } else if ( axis_order_ == "YZX" ) {
        // 0: x -> y
        dipole_[x]->copy(old_d[y]);
        nuc_dip_x_ = old_nuc_dip_y;
        // 1: y -> z
        dipole_[y]->copy(old_d[z]);
        nuc_dip_y_ = old_nuc_dip_z;
        // 2: z -> x
        dipole_[z]->copy(old_d[x]);
        nuc_dip_z_ = old_nuc_dip_x;

        // 0: xx -> yy
        swap_quadrupoles(xx,yy,false);
        // 1: xy -> yz
        swap_quadrupoles(xy,yz,false);
        // 2: xz -> yx <- xy
        swap_quadrupoles(xz,xy,true);
        // 3: yy -> zz
        swap_quadrupoles(yy,zz,false);
        // 4: yz -> zx <- xz
        swap_quadrupoles(yz,xz,true);
        // 5: zz -> xx
        swap_quadrupoles(zz,xx,false);
    } else if ( axis_order_ == "ZXY" ) {
        // 0: x -> z
        dipole_[x]->copy(old_d[z]);
        nuc_dip_x_ = old_nuc_dip_z;
        // 1: y -> x
        dipole_[y]->copy(old_d[x]);
        nuc_dip_y_ = old_nuc_dip_x;
        // 2: z -> y
        dipole_[z]->copy(old_d[y]);
        nuc_dip_z_ = old_nuc_dip_y;

        // 0: xx -> zz
        swap_quadrupoles(xx,zz,false);
        // 1: xy -> zx <- xz
        swap_quadrupoles(xy,xz,true);
        // 2: xz -> zy <- yz
        swap_quadrupoles(xz,yz,true);
        // 3: yy -> xx
        swap_quadrupoles(yy,xx,false);
        // 4: yz -> xy
        swap_quadrupoles(yz,xy,false);
        // 5: zz -> yy
        swap_quadrupoles(zz,yy,false);
    } else if ( axis_order_ == "XZY" ) {
        // 0: x -> x (do nothing)
        // 1: y -> z
        dipole_[y]->copy(old_d[z]);
        nuc_dip_y_ = old_nuc_dip_z;
        // 2: z -> y
        dipole_[z]->copy(old_d[y]);
        nuc_dip_z_ = old_nuc_dip_y;

        // 0: xx -> xx (do nothing)
        // 1: xy -> xz
        swap_quadrupoles(xy,xz,false);
        // 2: xz -> xy
        swap_quadrupoles(xz,xy,false);
        // 3: yy -> zz
        swap_quadrupoles(yy,zz,false);
        // 4: yz -> zy <- yz
        swap_quadrupoles(yz,yz,true);
        // 5: zz -> yy
        swap_quadrupoles(zz,yy,false);

    } else if ( axis_order_ == "ZYX" ) {
        // 0: x -> z
        dipole_[x]->copy(old_d[z]);
        nuc_dip_x_ = old_nuc_dip_z;
        // 1: y -> y (do nothing)
        // 2: z -> x
        dipole_[z]->copy(old_d[x]);
        nuc_dip_z_ = old_nuc_dip_x;

        // 0: xx -> zz
        swap_quadrupoles(xx,zz,false);
        // 1: xy -> zy <- yz
        swap_quadrupoles(xy,yz,true);
        // 2: xz -> zx <- xz
        swap_quadrupoles(xz,xz,true);
        // 3: yy -> yy (do nothing)
        // 4: yz -> yx <- xy
        swap_quadrupoles(yz,xy,true);
        // 5: zz -> xx
        swap_quadrupoles(zz,xx,false);
    } else if ( axis_order_ == "YXZ" ) {
        // 0: x -> y
        dipole_[x]->copy(old_d[y]);
        nuc_dip_x_ = old_nuc_dip_y;
        // 1: y -> x
        dipole_[y]->copy(old_d[x]);
        nuc_dip_y_ = old_nuc_dip_x;
        // 2: z -> z (do nothing)

        // 0: xx -> yy
        swap_quadrupoles(xx,yy,false);
        // 1: xy -> yx <- xy 
        swap_quadrupoles(xy,xy,true);
        // 2: xz -> yz
        swap_quadrupoles(xz,yz,false);
        // 3: yy -> xx
        swap_quadrupoles(yy,xx,false);
        // 4: yz -> xz
        swap_quadrupoles(yz,xz,false);
        // 5: zz -> zz (do nothing)
    } else {
        throw PsiException("invalid choice for ROTATE_POLARIZATION_AXIS",__FILE__,__LINE__);
    }
    

    double lambda_x = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
    double lambda_y = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

    //quadrupole[0]->scale(0.5 * lambda_x * lambda_x);
    //quadrupole[1]->scale(1.0 * lambda_x * lambda_y); // scaled by 2 for xy + yx
    //quadrupole[2]->scale(1.0 * lambda_x * lambda_z); // scaled by 2 for xz + zx
    //quadrupole[3]->scale(0.5 * lambda_y * lambda_y);
    //quadrupole[4]->scale(1.0 * lambda_y * lambda_z); // scaled by 2 for yz + zy
    //quadrupole[5]->scale(0.5 * lambda_z * lambda_z);

    quadrupole_scaled_sum_ = (std::shared_ptr<Matrix>)(new Matrix(nso_,nso_));
    quadrupole_scaled_sum_->zero();

    quadrupole_scaled_sum_->axpy(0.5 * lambda_x * lambda_x, quadrupole_[0]);
    quadrupole_scaled_sum_->axpy(1.0 * lambda_x * lambda_y, quadrupole_[1]); // scaled by 2 for xy + yx
    quadrupole_scaled_sum_->axpy(1.0 * lambda_x * lambda_z, quadrupole_[2]); // scaled by 2 for xz + zx
    quadrupole_scaled_sum_->axpy(0.5 * lambda_y * lambda_y, quadrupole_[3]);
    quadrupole_scaled_sum_->axpy(1.0 * lambda_y * lambda_z, quadrupole_[4]); // scaled by 2 for yz + zy
    quadrupole_scaled_sum_->axpy(0.5 * lambda_z * lambda_z, quadrupole_[5]);

    //for (int i = 0; i < 6; i++) {
    //    quadrupole_scaled_sum_->add(quadrupole[i]);
    //}

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
    if ( use_coherent_state_basis_ && options_.get_bool("QED_USE_RELAXED_ORBITALS") ) {
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

    // constant terms:  0.5 ( lambda . dn ) ^2 
    // 
    // or ... in the coherent state basis:
    // 
    // constant terms:  0.5 ( lambda . ( dn - <d> ) )^2
    average_electric_dipole_self_energy_ = 0.5 * nuc_dipdot * nuc_dipdot;

}

/*
 *
 * evaluate_dipole_self_energy():
 *
 * evaluate one- and two-electron contributions to the dipole self 
 * energy in the coherent-state basis: 
 *
 * < ( lambda.[mu - <mu>] )^2 >  =  lambda^2 ( <mu^2> - <mu>^2 )
 *
 * or not
 *
 * < (lambda.mu)^2 > = ( lambda . <mu> )^2
*/
void PolaritonicHF::evaluate_dipole_self_energy() {

    if ( use_coherent_state_basis_ ) {

        double one_electron = 0.0;
        double two_electron = 0.0;

        if ( n_photon_states_ < 1 ) return;

        // one-electron part if <mu^2> depends on quadrupole integrals: -Tr(D.q)
        std::shared_ptr<Matrix> oei (new Matrix(quadrupole_scaled_sum_));
        oei->scale(-1.0);

        one_electron += Da_->vector_dot(oei);
        if ( same_a_b_dens_ ) {
            one_electron *= 2.0;
        }else {
            one_electron += Db_->vector_dot(oei);
        }

        // two-electron part of <mu^2> is just the exchange contribution
        std::shared_ptr<Matrix> dipole_Ka (new Matrix(nso_,nso_));
        std::shared_ptr<Matrix> dipole_Kb (new Matrix(nso_,nso_));

        dipole_Ka->zero();
        dipole_Kb->zero();

        // Kpq += mu_pr * mu_qs * Drs
        double ** dp  = dipole_scaled_sum_->pointer();
        double ** dap = Da_->pointer();
        double ** dbp = Db_->pointer();
        double ** kap = dipole_Ka->pointer();
        double ** kbp = dipole_Kb->pointer();

        std::shared_ptr<Matrix> tmp (new Matrix(nso_,nso_));
        double ** tp = tmp->pointer();

        C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dp[0][0]),nso_,&(dap[0][0]),nso_,0.0,&(tp[0][0]),nso_);
        C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dp[0][0]),nso_,0.0,&(kap[0][0]),nso_);

        if ( !same_a_b_dens_ ) {
            C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dp[0][0]),nso_,&(dbp[0][0]),nso_,0.0,&(tp[0][0]),nso_);
            C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dp[0][0]),nso_,0.0,&(kbp[0][0]),nso_);
        }

        two_electron -= 0.5 * Da_->vector_dot(dipole_Ka);
        if ( same_a_b_dens_ ) {
            two_electron *= 2.0;
        }else {
            two_electron -= 0.5 * Db_->vector_dot(dipole_Kb);
        }

        outfile->Printf("\n");
        outfile->Printf("    ==> dipole self-energy: 1/2 lambda^2 ( <mu_e^2> - <mu_e>^2 ) <==\n");
        outfile->Printf("\n");
        outfile->Printf("    one-electron part;  %20.12lf\n",one_electron);
        outfile->Printf("    two-electron part:  %20.12lf\n",two_electron);
        outfile->Printf("    total:              %20.12lf\n",one_electron + two_electron);
        outfile->Printf("\n");

    }else{

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

        double lambda_x = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
        double lambda_y = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
        double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

        // 1/2 (lambda.<mu>)^2
        double dse = lambda_x * tot_dip_x_ + lambda_y * tot_dip_y_ + lambda_z * tot_dip_z_;
        dse *= 0.5 * dse;

        outfile->Printf("\n");
        outfile->Printf("    ==> dipole self-energy <==\n");
        outfile->Printf("\n");
        outfile->Printf("    1/2 ( lambda . <mu> )^2: %20.12lf\n",dse);
        outfile->Printf("\n");

        // add dse to scf energy. only necessary if not relaxing the orbitals and not using coherent state basis
        if ( !options_.get_bool("QED_USE_RELAXED_ORBITALS") ) {
            energy_ += dse;
        }
    }
}

/*
 *
 * evaluate_dipole_variance():
 *
 * evaluate one- and two-electron contributions to the dipole variance
 *
 * < ( [mu - <mu>] )^2 > = <mu_e>^2 - <mu_e>^2
 *
 * technically, this quantity should have 6 components:
 *
 * x.x
 * x.y + y.x
 * x.z + z.x
 * y.y
 * y.z + z.y
 * z.z
 *
*/
void PolaritonicHF::evaluate_dipole_variance() {

    // evaluate the dipole moment integrals in the molecular basis

    double dipole_mo_x = 0.0;
    double dipole_mo_y = 0.0;
    double dipole_mo_z = 0.0;

    dipole_mo_x = C_DDOT(nso_*nso_,&(Da_->pointer())[0][0],1,&(dipole_[0]->pointer())[0][0],1);
    dipole_mo_y = C_DDOT(nso_*nso_,&(Da_->pointer())[0][0],1,&(dipole_[1]->pointer())[0][0],1);
    dipole_mo_z = C_DDOT(nso_*nso_,&(Da_->pointer())[0][0],1,&(dipole_[2]->pointer())[0][0],1);
    if ( same_a_b_dens_ ) {
        dipole_mo_x *= 2.0;
        dipole_mo_y *= 2.0;
        dipole_mo_z *= 2.0;
    }else {
        dipole_mo_x += C_DDOT(nso_*nso_,&(Db_->pointer())[0][0],1,&(dipole_[0]->pointer())[0][0],1);
        dipole_mo_y += C_DDOT(nso_*nso_,&(Db_->pointer())[0][0],1,&(dipole_[1]->pointer())[0][0],1);
        dipole_mo_z += C_DDOT(nso_*nso_,&(Db_->pointer())[0][0],1,&(dipole_[2]->pointer())[0][0],1);
    }

    outfile->Printf("\n");
    outfile->Printf("    ==> dipole integrals <== \n");
    outfile->Printf("\n");
    outfile->Printf("    <mu_x> = %20.12lf\n",dipole_mo_x);
    outfile->Printf("    <mu_y> = %20.12lf\n",dipole_mo_y);
    outfile->Printf("    <mu_z> = %20.12lf\n",dipole_mo_z);
    outfile->Printf("\n");

    double one_electron_xx = 0.0;
    double one_electron_xy = 0.0;
    double one_electron_xz = 0.0;
    double one_electron_yy = 0.0;
    double one_electron_yz = 0.0;
    double one_electron_zz = 0.0;

    double two_electron_xx = 0.0;
    double two_electron_xy = 0.0;
    double two_electron_xz = 0.0;
    double two_electron_yy = 0.0;
    double two_electron_yz = 0.0;
    double two_electron_zz = 0.0;

    // one-electron part depends on quadrupole integrals: -Tr(D.q)

    one_electron_xx =       Da_->vector_dot(quadrupole_[0]);
    one_electron_xy = 2.0 * Da_->vector_dot(quadrupole_[1]);
    one_electron_xz = 2.0 * Da_->vector_dot(quadrupole_[2]);
    one_electron_yy =       Da_->vector_dot(quadrupole_[3]);
    one_electron_yz = 2.0 * Da_->vector_dot(quadrupole_[4]);
    one_electron_zz =       Da_->vector_dot(quadrupole_[5]);
    if ( same_a_b_dens_ ) {
        one_electron_xx *= 2.0;
        one_electron_xy *= 2.0;
        one_electron_xz *= 2.0;
        one_electron_yy *= 2.0;
        one_electron_yz *= 2.0;
        one_electron_zz *= 2.0;
    }else {
        one_electron_xx +=       Db_->vector_dot(quadrupole_[0]);
        one_electron_xy += 2.0 * Db_->vector_dot(quadrupole_[1]);
        one_electron_xz += 2.0 * Db_->vector_dot(quadrupole_[2]);
        one_electron_yy +=       Db_->vector_dot(quadrupole_[3]);
        one_electron_yz += 2.0 * Db_->vector_dot(quadrupole_[4]);
        one_electron_zz +=       Db_->vector_dot(quadrupole_[5]);
    }

    one_electron_xx *= -1.0;
    one_electron_xy *= -1.0;
    one_electron_xz *= -1.0;
    one_electron_yy *= -1.0;
    one_electron_yz *= -1.0;
    one_electron_zz *= -1.0;

    // two-electron part is just exchange part of <mu^2>
    std::shared_ptr<Matrix> dipole_Ka (new Matrix(nso_,nso_));
    std::shared_ptr<Matrix> dipole_Kb (new Matrix(nso_,nso_));

    dipole_Ka->zero();
    dipole_Kb->zero();

    // Kpq += mu_pr * mu_qs * Drs

    double ** dx  = dipole_[0]->pointer();
    double ** dy  = dipole_[1]->pointer();
    double ** dz  = dipole_[2]->pointer();

    double ** dap = Da_->pointer();
    double ** dbp = Db_->pointer();
    double ** kap = dipole_Ka->pointer();
    double ** kbp = dipole_Kb->pointer();

    std::shared_ptr<Matrix> tmp (new Matrix(nso_,nso_));
    double ** tp = tmp->pointer();

    // xx

    C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dx[0][0]),nso_,&(dap[0][0]),nso_,0.0,&(tp[0][0]),nso_);
    C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dx[0][0]),nso_,0.0,&(kap[0][0]),nso_);

    if ( !same_a_b_dens_ ) {
        C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dx[0][0]),nso_,&(dbp[0][0]),nso_,0.0,&(tp[0][0]),nso_);
        C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dx[0][0]),nso_,0.0,&(kbp[0][0]),nso_);
    }

    two_electron_xx -= 0.5 * Da_->vector_dot(dipole_Ka);
    if ( same_a_b_dens_ ) {
        two_electron_xx *= 2.0;
    }else {
        two_electron_xx -= 0.5 * Db_->vector_dot(dipole_Kb);
    }

    // xy

    C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dx[0][0]),nso_,&(dap[0][0]),nso_,0.0,&(tp[0][0]),nso_);
    C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dy[0][0]),nso_,0.0,&(kap[0][0]),nso_);

    if ( !same_a_b_dens_ ) {
        C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dx[0][0]),nso_,&(dbp[0][0]),nso_,0.0,&(tp[0][0]),nso_);
        C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dy[0][0]),nso_,0.0,&(kbp[0][0]),nso_);
    }

    two_electron_xy -= 0.5 * Da_->vector_dot(dipole_Ka);
    if ( same_a_b_dens_ ) {
        two_electron_xy *= 2.0;
    }else {
        two_electron_xy -= 0.5 * Db_->vector_dot(dipole_Kb);
    }

    // xz

    C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dx[0][0]),nso_,&(dap[0][0]),nso_,0.0,&(tp[0][0]),nso_);
    C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dz[0][0]),nso_,0.0,&(kap[0][0]),nso_);

    if ( !same_a_b_dens_ ) {
        C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dx[0][0]),nso_,&(dbp[0][0]),nso_,0.0,&(tp[0][0]),nso_);
        C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dz[0][0]),nso_,0.0,&(kbp[0][0]),nso_);
    }

    two_electron_xz -= 0.5 * Da_->vector_dot(dipole_Ka);
    if ( same_a_b_dens_ ) {
        two_electron_xz *= 2.0;
    }else {
        two_electron_xz -= 0.5 * Db_->vector_dot(dipole_Kb);
    }

    // yy

    C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dy[0][0]),nso_,&(dap[0][0]),nso_,0.0,&(tp[0][0]),nso_);
    C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dy[0][0]),nso_,0.0,&(kap[0][0]),nso_);

    if ( !same_a_b_dens_ ) {
        C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dy[0][0]),nso_,&(dbp[0][0]),nso_,0.0,&(tp[0][0]),nso_);
        C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dy[0][0]),nso_,0.0,&(kbp[0][0]),nso_);
    }

    two_electron_yy -= 0.5 * Da_->vector_dot(dipole_Ka);
    if ( same_a_b_dens_ ) {
        two_electron_yy *= 2.0;
    }else {
        two_electron_yy -= 0.5 * Db_->vector_dot(dipole_Kb);
    }

    // yz

    C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dy[0][0]),nso_,&(dap[0][0]),nso_,0.0,&(tp[0][0]),nso_);
    C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dz[0][0]),nso_,0.0,&(kap[0][0]),nso_);

    if ( !same_a_b_dens_ ) {
        C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dy[0][0]),nso_,&(dbp[0][0]),nso_,0.0,&(tp[0][0]),nso_);
        C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dz[0][0]),nso_,0.0,&(kbp[0][0]),nso_);
    }

    two_electron_yz -= 0.5 * Da_->vector_dot(dipole_Ka);
    if ( same_a_b_dens_ ) {
        two_electron_yz *= 2.0;
    }else {
        two_electron_yz -= 0.5 * Db_->vector_dot(dipole_Kb);
    }

    // zz

    C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dz[0][0]),nso_,&(dap[0][0]),nso_,0.0,&(tp[0][0]),nso_);
    C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dz[0][0]),nso_,0.0,&(kap[0][0]),nso_);

    if ( !same_a_b_dens_ ) {
        C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dz[0][0]),nso_,&(dbp[0][0]),nso_,0.0,&(tp[0][0]),nso_);
        C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dz[0][0]),nso_,0.0,&(kbp[0][0]),nso_);
    }

    two_electron_zz -= 0.5 * Da_->vector_dot(dipole_Ka);
    if ( same_a_b_dens_ ) {
        two_electron_zz *= 2.0;
    }else {
        two_electron_zz -= 0.5 * Db_->vector_dot(dipole_Kb);
    }

    outfile->Printf("\n");
    outfile->Printf("    ==> dipole variance: ( <mu_e^2> - <mu_e>^2 ) <== \n");
    outfile->Printf("\n");
    outfile->Printf("    one electron (xx): %20.12lf\n",one_electron_xx);
    outfile->Printf("    one electron (xy): %20.12lf\n",one_electron_xy);
    outfile->Printf("    one electron (xz): %20.12lf\n",one_electron_xz);
    outfile->Printf("    one electron (yy): %20.12lf\n",one_electron_yy);
    outfile->Printf("    one electron (yz): %20.12lf\n",one_electron_yz);
    outfile->Printf("    one electron (zz): %20.12lf\n",one_electron_zz);
    outfile->Printf("\n");
    outfile->Printf("    two electron (xx): %20.12lf\n",two_electron_xx * 2.0);
    outfile->Printf("    two electron (xy): %20.12lf\n",two_electron_xy * 2.0);
    outfile->Printf("    two electron (xz): %20.12lf\n",two_electron_xz * 2.0);
    outfile->Printf("    two electron (yy): %20.12lf\n",two_electron_yy * 2.0);
    outfile->Printf("    two electron (yz): %20.12lf\n",two_electron_yz * 2.0);
    outfile->Printf("    two electron (zz): %20.12lf\n",two_electron_zz * 2.0);
    outfile->Printf("\n");
    outfile->Printf("    total (xx):        %20.12lf\n",one_electron_xx + two_electron_xx * 2.0);
    outfile->Printf("    total (xy):        %20.12lf\n",one_electron_xy + two_electron_xy * 2.0);
    outfile->Printf("    total (xz):        %20.12lf\n",one_electron_xz + two_electron_xz * 2.0);
    outfile->Printf("    total (yy):        %20.12lf\n",one_electron_yy + two_electron_yy * 2.0);
    outfile->Printf("    total (yz):        %20.12lf\n",one_electron_yz + two_electron_yz * 2.0);
    outfile->Printf("    total (zz):        %20.12lf\n",one_electron_zz + two_electron_zz * 2.0);
    outfile->Printf("\n");

}

} // End namespaces
