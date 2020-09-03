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

    print_cavity_properties_ = false;

    // initialize cavity parameters
    initialize_cavity();
}

void PolaritonicHF::initialize_cavity() {

    n_photon_states_ = options_.get_int("N_PHOTON_STATES");

    // cavity Hamiltonian in basis of photon number states
    HCavity_x_               = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    HCavity_y_               = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    HCavity_z_               = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));

    // molecule->cavity interactoin Hamiltonian in basis of photon number states
    HCavityInteraction_x_    = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    HCavityInteraction_y_    = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    HCavityInteraction_z_    = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));

    // total cavity Hamiltonian in basis of photon number states
    HCavityTotal_x_          = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    HCavityTotal_y_          = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    HCavityTotal_z_          = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));

    // cavity dipole operator in basis of photon number states
    CavityDipole_x_          = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    CavityDipole_y_          = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    CavityDipole_z_          = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));

    // dipole potential felt by molecule due to presence of
    // cavity.  Memory will be allocated later.
    //CavityDipolePotential_x_ = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    //CavityDipolePotential_y_ = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    //CavityDipolePotential_z_ = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    // Read in the energies of each cavity mode

    cavity_e_ = (double*)malloc(sizeof(double)*3);
    memset((void*)cavity_e_,'\0',3*sizeof(double));
    if (options_["CAVITY_FREQUENCY"].has_changed()){
       if (options_["CAVITY_FREQUENCY"].size() != 3)
          throw PsiException("The CAVITY E array has the wrong dimensions",__FILE__,__LINE__);
       for (int i = 0; i < 3; i++) cavity_e_[i] = options_["CAVITY_FREQUENCY"][i].to_double();
    }else{
       cavity_e_[0] = 2.042/27.21138;  // Energy for the Au nanoparticle taken from Gray's paper
       cavity_e_[1] = 2.042/27.21138;
       cavity_e_[2] = 2.042/27.21138;
    }

    // Read in the cavity transition dipole moment

    cavity_tdm_ = (double*)malloc(sizeof(double)*3);
    memset((void*)cavity_tdm_,'\0',3*sizeof(double));
    if (options_["CAVITY_TDM"].has_changed()){
       if (options_["CAVITY_TDM"].size() != 3)
          throw PsiException("The CAVITY TDM array has the wrong dimensions",__FILE__,__LINE__);
       for (int i = 0; i < 3; i++) cavity_tdm_[i] = options_["CAVITY_TDM"][i].to_double();
    }else{
       cavity_tdm_[0] = 2990.0/2.54175;
       cavity_tdm_[1] = 2990.0/2.54175;
       cavity_tdm_[2] = 2990.0/2.54175;
    }

    // Read in the cavity coordinates

    cavity_coordinates_ = (double*)malloc(sizeof(double)*3);
    memset((void*)cavity_coordinates_,'\0',3*sizeof(double));
    if (options_["CAVITY_COORDINATES"].has_changed()){
       if (options_["CAVITY_COORDINATES"].size() != 3)
          throw PsiException("The CAVITY COORDINATES array has the wrong dimensions",__FILE__,__LINE__);
       for (int i = 0; i < 3; i++) cavity_coordinates_[i] = options_["CAVITY_COORDINATES"][i].to_double();
    }else{
       cavity_coordinates_[0] =  0.0;
       cavity_coordinates_[1] =  0.0;
       cavity_coordinates_[2] = 40.0;
    }

    // compute molecule center-of-mass

    std::shared_ptr<Molecule>mol=Process::environment.molecule();
    int natom = mol->natom();
    center_of_mass_ = (double*)malloc(sizeof(double)*3);
    memset((void*)center_of_mass_,'\0',3*sizeof(double));

    double temp_x = 0.0;
    double temp_y = 0.0;
    double temp_z = 0.0;
    double temp_m = 0.0;
    for (int i = 0; i < natom; i++){
       temp_x += mol->mass(i) * mol->x(i);
       temp_y += mol->mass(i) * mol->y(i);
       temp_z += mol->mass(i) * mol->z(i);
       temp_m += mol->mass(i);
    }
    center_of_mass_[0] = temp_x/temp_m;
    center_of_mass_[1] = temp_y/temp_m;
    center_of_mass_[2] = temp_z/temp_m;

    // Get Nuclear contribution to the molecular total dipole moment

    nuc_dip_x_ = 0.0;
    nuc_dip_y_ = 0.0;
    nuc_dip_z_ = 0.0;
    for (int i = 0; i < natom; i++) {
        nuc_dip_x_ += mol->Z(i) * mol->x(i);
        nuc_dip_y_ += mol->Z(i) * mol->y(i);
        nuc_dip_z_ += mol->Z(i) * mol->z(i);
    }

    std::shared_ptr<MintsHelper> mints (new MintsHelper(reference_wavefunction_));
    dipole_ = mints->so_dipole();

    // Build dipole potential integrals that represent the influence of the
    // cavity dipole moment on the molecule

    dipole_potential_integrals();
}

/* 
 * build_cavity_hamiltonian():
 *
 * Build the cavity Hamiltonian, the molecule->cavity interaction
 * Hamiltonian, and the cavity dipole moment operator in the
 * basis of photon number states.  Transform all three operators
 * to the basis that diagonalizes the total cavity Hamiltonian
 * (the sum of the cavity Hamiltonian and the molecule->cavity
 * Hamiltonian).
*/
void PolaritonicHF::build_cavity_hamiltonian(){

    // First, build the cavity dipole moment operator 
    // in the basis of photon number states ( n_photon_states_ )

    int nS = n_photon_states_;

    CavityDipole_x_->zero();
    CavityDipole_y_->zero();
    CavityDipole_z_->zero();

    for (int A=0; A<nS-1; A++){
        CavityDipole_x_->pointer()[A+1][A] += cavity_tdm_[0]*sqrt(A+1);
        CavityDipole_x_->pointer()[A][A+1] += cavity_tdm_[0]*sqrt(A+1);
        CavityDipole_y_->pointer()[A+1][A] += cavity_tdm_[1]*sqrt(A+1);
        CavityDipole_y_->pointer()[A][A+1] += cavity_tdm_[1]*sqrt(A+1);
        CavityDipole_z_->pointer()[A+1][A] += cavity_tdm_[2]*sqrt(A+1);
        CavityDipole_z_->pointer()[A][A+1] += cavity_tdm_[2]*sqrt(A+1);
    }

    // Build the cavity Hamiltonian operator in the basis of 
    // photon number states ( n_photon_states_ )

    HCavity_x_->zero();
    HCavity_y_->zero();
    HCavity_z_->zero();

    for (int A=0; A<nS; A++){
        HCavity_x_->pointer()[A][A] = A*cavity_e_[0];
        HCavity_y_->pointer()[A][A] = A*cavity_e_[1];
        HCavity_z_->pointer()[A][A] = A*cavity_e_[2];
    }

    // Evaluate the distance from the center of mass of the molecule to cavity

    double r = (center_of_mass_[0]-cavity_coordinates_[0])*(center_of_mass_[0]-cavity_coordinates_[0])
             + (center_of_mass_[1]-cavity_coordinates_[1])*(center_of_mass_[1]-cavity_coordinates_[1])
             + (center_of_mass_[2]-cavity_coordinates_[2])*(center_of_mass_[2]-cavity_coordinates_[2]);
    r = sqrt(r);

    double delta_x = cavity_coordinates_[0] - center_of_mass_[0];
    double delta_y = cavity_coordinates_[1] - center_of_mass_[1];
    double delta_z = cavity_coordinates_[2] - center_of_mass_[2];

    double * r_vector;

    r_vector = (double*)malloc(sizeof(double)*3);
    memset((void*)r_vector,'\0',3*sizeof(double));

    r_vector[0] = delta_x;
    r_vector[1] = delta_y;
    r_vector[2] = delta_z;

    double oer3 = 1.0 /(r*r*r);
    double coupling_strength = 1.0 * oer3;

    // Evaluate the electronic contribute to the molecule's dipole moment

    double e_dip_x = 2.0*C_DDOT(nso_*nso_,&(Da_->pointer())[0][0],1,&(dipole_[0]->pointer())[0][0],1);
    double e_dip_y = 2.0*C_DDOT(nso_*nso_,&(Da_->pointer())[0][0],1,&(dipole_[1]->pointer())[0][0],1);
    double e_dip_z = 2.0*C_DDOT(nso_*nso_,&(Da_->pointer())[0][0],1,&(dipole_[2]->pointer())[0][0],1);

    // Build the molecule->cavity interaction Hamiltonian operator 
    // in the basis of photon number states ( n_photon_states_ )

    HCavityInteraction_x_->zero();
    HCavityInteraction_y_->zero();
    HCavityInteraction_z_->zero();

    for (int A = 0; A<nS; A++) {
        if (A < nS - 1) {

            HCavityInteraction_x_->pointer()[A][A+1] += (e_dip_x + nuc_dip_x_)*cavity_tdm_[0];
            HCavityInteraction_x_->pointer()[A][A+1] -= 3.0*cavity_tdm_[0]*delta_x
                                                       * ( r_vector[0]*(nuc_dip_x_ + e_dip_x)
                                                         + r_vector[1]*(nuc_dip_y_ + e_dip_y)
                                                         + r_vector[2]*(nuc_dip_z_ + e_dip_z) ) / (r*r);
            HCavityInteraction_x_->pointer()[A][A+1] *= coupling_strength*sqrt(A+1);

            HCavityInteraction_y_->pointer()[A][A+1] += (e_dip_y + nuc_dip_y_)*cavity_tdm_[1];
            HCavityInteraction_y_->pointer()[A][A+1] -= 3.0*cavity_tdm_[1]*delta_y
                                                       * ( r_vector[0]*(nuc_dip_x_ + e_dip_x)
                                                         + r_vector[1]*(nuc_dip_y_ + e_dip_y)
                                                         + r_vector[2]*(nuc_dip_z_ + e_dip_z) ) / (r*r);
            HCavityInteraction_y_->pointer()[A][A+1] *= coupling_strength*sqrt(A+1);

            HCavityInteraction_z_->pointer()[A][A+1] += (e_dip_z + nuc_dip_z_)*cavity_tdm_[2];
            HCavityInteraction_z_->pointer()[A][A+1] -= 3.0*cavity_tdm_[2]*delta_z
                                                       * ( r_vector[0]*(nuc_dip_x_ + e_dip_x)
                                                         + r_vector[1]*(nuc_dip_y_ + e_dip_y)
                                                         + r_vector[2]*(nuc_dip_z_ + e_dip_z) ) / (r*r);
            HCavityInteraction_z_->pointer()[A][A+1] *= coupling_strength*sqrt(A+1);

        }

        if (A > 0) {

            HCavityInteraction_x_->pointer()[A][A-1] += (e_dip_x + nuc_dip_x_)*cavity_tdm_[0];
            HCavityInteraction_x_->pointer()[A][A-1] -= 3.0*cavity_tdm_[0]*delta_x
                                                       * ( r_vector[0]*(nuc_dip_x_ + e_dip_x)
                                                         + r_vector[1]*(nuc_dip_y_ + e_dip_y)
                                                         + r_vector[2]*(nuc_dip_z_ + e_dip_z) ) / (r*r);
            HCavityInteraction_x_->pointer()[A][A-1] *= coupling_strength*sqrt(A);

            HCavityInteraction_y_->pointer()[A][A-1] += (e_dip_y + nuc_dip_y_)*cavity_tdm_[1];
            HCavityInteraction_y_->pointer()[A][A-1] -= 3.0*cavity_tdm_[1]*delta_y
                                                       * ( r_vector[0]*(nuc_dip_x_ + e_dip_x)
                                                         + r_vector[1]*(nuc_dip_y_ + e_dip_y)
                                                         + r_vector[2]*(nuc_dip_z_ + e_dip_z)) / (r*r);
            HCavityInteraction_y_->pointer()[A][A-1] *= coupling_strength*sqrt(A);

            HCavityInteraction_z_->pointer()[A][A-1] += (e_dip_z + nuc_dip_z_)*cavity_tdm_[2];
            HCavityInteraction_z_->pointer()[A][A-1] -= 3.0*cavity_tdm_[2]*delta_z
                                                       * ( r_vector[0]*(nuc_dip_x_ + e_dip_x)
                                                         + r_vector[1]*(nuc_dip_y_ + e_dip_y)
                                                         + r_vector[2]*(nuc_dip_z_ + e_dip_z) ) / (r*r);
            HCavityInteraction_z_->pointer()[A][A-1] *= coupling_strength*sqrt(A);
        }
    }

    // Build the total cavity Hamiltonian

    HCavityTotal_x_->copy(HCavity_x_);
    HCavityTotal_y_->copy(HCavity_y_);
    HCavityTotal_z_->copy(HCavity_z_);
    HCavityTotal_x_->add(HCavityInteraction_x_);
    HCavityTotal_y_->add(HCavityInteraction_y_);
    HCavityTotal_z_->add(HCavityInteraction_z_);

    // Diagonalize total cavity Hamiltonian and transform cavity
    // Hamiltonian, molecule->cavity interaction Hamiltonian, and
    // cavity dipole moment operator to the basis that diagonalizes
    // the total cavity Hamiltonian
    std::shared_ptr<Matrix> eigvec (new Matrix(n_photon_states_,n_photon_states_));
    std::shared_ptr<Vector> eigval (new Vector(n_photon_states_));


    // x
    HCavityTotal_x_->diagonalize(eigvec, eigval);
    HCavity_x_->transform(eigvec);
    HCavityInteraction_x_->transform(eigvec);
    CavityDipole_x_->transform(eigvec);
    double ex = eigval->pointer()[0];

    // y
    HCavityTotal_y_->diagonalize(eigvec, eigval);
    HCavity_y_->transform(eigvec);
    HCavityInteraction_y_->transform(eigvec);
    CavityDipole_y_->transform(eigvec);
    double ey = eigval->pointer()[0];

    // z
    HCavityTotal_z_->diagonalize(eigvec, eigval);
    HCavity_z_->transform(eigvec);
    HCavityInteraction_z_->transform(eigvec);
    CavityDipole_z_->transform(eigvec);
    double ez = eigval->pointer()[0];

    if ( print_cavity_properties_ ) {
        outfile->Printf("\n");
        outfile->Printf("    ==> Cavity Properties <==\n");
        outfile->Printf("\n");

        outfile->Printf("        Energy x: %20.12lf    Dipole x: %20.12lf\n",ex,CavityDipole_x_->pointer()[0][0]);
        outfile->Printf("        Energy y: %20.12lf    Dipole y: %20.12lf\n",ey,CavityDipole_y_->pointer()[0][0]);
        outfile->Printf("        Energy z: %20.12lf    Dipole z: %20.12lf\n",ez,CavityDipole_z_->pointer()[0][0]);
        outfile->Printf("\n");
    }
}

/* 
 * dipole_potential_integrals(): 
 *
 * Evaluate dipole potential integrals that represent the potential
 * felt by the molecule due to the presence of the dipole moment
 * of the cavity
*/
void PolaritonicHF::dipole_potential_integrals() {

    std::shared_ptr<OneBodyAOInt> ints(reference_wavefunction_->integral()->ao_multipole_potential(3,0));

    int nbf = reference_wavefunction_->basisset()->nbf();
    int nao = reference_wavefunction_->basisset()->nao();

    std::vector< std::shared_ptr<Matrix> > mats;
    for(int i=0; i < 20; i++) {
        mats.push_back(std::shared_ptr<Matrix> (new Matrix(nao, nao)));
        mats[i]->zero();
    }

    // cavity dipole potential felt by the molecule
    Vector3 coords(cavity_coordinates_[0],cavity_coordinates_[1],cavity_coordinates_[2]);
    ints->set_origin(coords);
    ints->compute(mats);

    // ao/so transformation

    std::shared_ptr<Matrix> Vx = Matrix::triplet(AO2SO_,mats[1],AO2SO_,true,false,false);
    std::shared_ptr<Matrix> Vy = Matrix::triplet(AO2SO_,mats[2],AO2SO_,true,false,false);
    std::shared_ptr<Matrix> Vz = Matrix::triplet(AO2SO_,mats[3],AO2SO_,true,false,false);

    CavityDipolePotential_x_ = (std::shared_ptr<Matrix>)(new Matrix(Vx));
    CavityDipolePotential_y_ = (std::shared_ptr<Matrix>)(new Matrix(Vy));
    CavityDipolePotential_z_ = (std::shared_ptr<Matrix>)(new Matrix(Vz));
}

} // End namespaces