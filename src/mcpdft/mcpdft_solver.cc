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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <tuple>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#include "psi4/libpsi4util/libpsi4util.h"

#include "psi4/libqt/qt.h"

// jk object
#include "psi4/libfock/jk.h"

// for dft
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"

// for grid
#include "psi4/libfock/points.h"
#include "psi4/libfock/cubature.h"

#include "psi4/psi4-dec.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.hpp"

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/gshell.h"
#include "psi4/libmints/molecule.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libqt/qt.h"

// for potential object
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libscf_solver/hf.h"

// for reading 2RDM
#include "psi4/psi4-dec.h"
#include <psi4/psifiles.h>
#include <psi4/libiwl/iwl.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libtrans/integraltransform.h>

#include <psi4/libpsi4util/PsiOutStream.h>

// mcpdft 
#include "mcpdft_solver.h"

// blas
#include <misc/blas.h>
#include <misc/omp.h>

using namespace psi;
using namespace fnocc;

namespace hilbert{ 

// the MCPDFTSolver class derives from the Wavefunction class and inherits its members
MCPDFTSolver::MCPDFTSolver(std::shared_ptr<Wavefunction> reference_wavefunction,Options & options_):
    Wavefunction(options_){

    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

MCPDFTSolver::~MCPDFTSolver() {

}

// initialize members of the MCPDFTSolver class
void MCPDFTSolver::common_init() {

    outfile->Printf("\n\n");
    outfile->Printf( "        ********************************************************************\n");
    outfile->Printf( "        *                                                                  *\n");
    outfile->Printf( "        *    MCPDFT:                                                       *\n");
    outfile->Printf( "        *                                                                  *\n");
    outfile->Printf( "        *    Multiconfigurational Pair Density Functional Theory           *\n");
    outfile->Printf( "        *    M. Mostafanejad and A. E. DePrince III                        *\n");
    outfile->Printf( "        *                                                                  *\n");
    outfile->Printf( "        ********************************************************************\n");
    outfile->Printf("\n\n");

    is_low_memory_ = false;

    reference_energy_ = Process::environment.globals["V2RDM TOTAL ENERGY"];
    
    shallow_copy(reference_wavefunction_);

    // number of alpha electrons
    nalpha_   = reference_wavefunction_->nalpha();

    // number of beta electrons
    nbeta_    = reference_wavefunction_->nbeta();

    // number of alpha electrons per irrep
    nalphapi_ = reference_wavefunction_->nalphapi();

    // number of beta electrons per irrep
    nbetapi_  = reference_wavefunction_->nbetapi();

    // number of frozen core orbitals per irrep
    frzcpi_   = reference_wavefunction_->frzcpi();

    // number of frozen virtual orbitals per irrep
    frzvpi_   = reference_wavefunction_->frzvpi();

    // number of molecular orbials per irrep
    nmopi_    = reference_wavefunction_->nmopi();

    // number of symmetry orbials per irrep
    nsopi_    = reference_wavefunction_->nsopi();

    // number of irreducible representations
    nirrep_   = reference_wavefunction_->nirrep();

    // total number of symmetry orbitals
    nso_      = reference_wavefunction_->nso();

    // total number of molecular orbitals
    nmo_      = reference_wavefunction_->nmo();

    // grab the molecule from the reference wave function
    molecule_ = reference_wavefunction_->molecule();

    // SO/MO transformation matrices
    Ca_ = std::shared_ptr<Matrix>(reference_wavefunction_->Ca());
    Cb_ = std::shared_ptr<Matrix>(reference_wavefunction_->Cb());

    // overlap integrals
    S_  = std::shared_ptr<Matrix>(reference_wavefunction_->S());

    std::shared_ptr<Matrix> Sevec;
    std::shared_ptr<Vector> Seval;

    // SO-basis Fock matrices
    Fa_ = std::shared_ptr<Matrix>(reference_wavefunction_->Fa());
    Fb_ = std::shared_ptr<Matrix>(reference_wavefunction_->Fb());

    // SO-basis density matrices
    Da_ = std::shared_ptr<Matrix>(reference_wavefunction_->Da());
    Db_ = std::shared_ptr<Matrix>(reference_wavefunction_->Db());

    // orbital energies
    epsilon_a_ = std::make_shared<Vector>(nmopi_);
    epsilon_a_->copy(*reference_wavefunction_->epsilon_a());
    epsilon_b_ = std::make_shared<Vector>(nmopi_);
    epsilon_b_->copy(*reference_wavefunction_->epsilon_b());

    // orbital energies
    occupation_a_= std::shared_ptr<Vector>(new Vector(nmopi_));
    occupation_b_= std::shared_ptr<Vector>(new Vector(nmopi_));

    for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < nalphapi_[h]; i++) {
                    occupation_a_->set(h, i, 1.0);
            }
            for (int i = 0; i < nbetapi_[h]; i++) {
                    occupation_b_->set(h, i, 1.0);
            }
    }

    // occupation_a_->print();
    // set the wavefunction name
    name_ = "DFT";

    // restricted orbitals, unrestricted rdms
    same_a_b_orbs_ = true;
    same_a_b_dens_ = false;

    // symmetry of orbitals:
    symmetry_ = (int*)malloc(nmo_*sizeof(int));
    memset((void*)symmetry_,'\0',nmo_*sizeof(int));
    int count = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < nsopi_[h]; i++) {
            symmetry_[count++] = h;
        }
    }

    // pitzer offset
    pitzer_offset_ = (int*)malloc(nirrep_*sizeof(int));
    count = 0;
    for (int h = 0; h < nirrep_; h++) {
        pitzer_offset_[h] = count;
        count += nsopi_[h];
    }

    // evaluate basis function values on a grid:

    scf::HF* scfwfn = (scf::HF*)reference_wavefunction_.get();
    std::shared_ptr<SuperFunctional> functional = scfwfn->functional();
    potential_ = scfwfn->V_potential();

    // phi matrix (real-space <- AO basis mapping)
    // since grid is stored in blocks, we need to build a full phi matrix
    // from the individual blocks:

    std::shared_ptr<BasisSet> primary = reference_wavefunction_->get_basisset("ORBITAL");
    int nao = primary->nbf();
    std::shared_ptr<Matrix> Da_ao (new Matrix(nao,nao));
    std::shared_ptr<Matrix> Db_ao (new Matrix(nao,nao));

    points_func_ = potential_->properties()[0];
    points_func_->set_pointers(Da_ao,Db_ao);

    // determine number of grid points
    int nblocks = potential_->nblocks();
    phi_points_    = 0;
    max_functions_ = 0;
    max_points_    = 0;
    for (int myblock = 0; myblock < nblocks; myblock++) {
        std::shared_ptr<BlockOPoints> block = potential_->get_block(myblock);
        points_func_->compute_points(block);
        int npoints = block->npoints();
        phi_points_ += npoints;
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();
        if ( nlocal > max_functions_ ) max_functions_ = nlocal;
        if ( npoints > max_points_ )   max_points_    = npoints;
    }

    // what is the derivative level?
    deriv_   = functional->deriv();
    is_gga_  = functional->is_gga();
    is_meta_ = functional->is_meta();
    is_unpolarized_ = functional->is_unpolarized();

    std::vector < std::shared_ptr<Matrix> > JK = BuildJK();

    double caa = Da_->vector_dot(JK[0]);
    double cab = Da_->vector_dot(JK[1]);
    double cba = Db_->vector_dot(JK[0]);
    double cbb = Db_->vector_dot(JK[1]);

    coulomb_energy_ = 0.5 * ( caa + cab + cba + cbb );

    /* ==============================================
       estimating the memory requirement for MCPDFT
       ============================================== */

    outfile->Printf("\n"); 
    outfile->Printf("    ==> Memory requirements <==\n");
    outfile->Printf("\n");
    
    // memory is from process::environment
    memory_ = Process::environment.get_memory();

    // assume 500 MB is already allocated by psi4
    long int n_extra       = 500 * 1024 * 1024 / 8;
    long int n_mem         = n_extra;
    long int available_mem = 0;
    long int n_init        = 0;
    long int n_phiAO       = 0;
    long int n_grids       = 0;
    long int n_transf      = 0;
    long int n_rho         = 0;
    long int n_pi          = 0;
    long int n_R           = 0;
    long int n_transl      = 0;
    long int n_df          = 0;

    // default basic vectors and matrices initialized using common_init()
    n_init  = nso_ * nso_;                                 // S_ matrix
    n_init += 2 * nso_ * nso_;                             // Ca_ and Cb_ matrices
    n_init += 2 * nso_ * nso_;                             // Fa_ and Fb_ matrices
    n_init += 2 * nmo_;                                    // epsilon_a_ and epsilon_b_ vectors
    n_mem += n_init;

    // phiAO, phiAO_x, phiAO_y, phiAO_z
    n_phiAO = phi_points_ * nso_;
    if ( is_gga_ || is_meta_ ) { 
       n_phiAO   += 3 * phi_points_ * nso_;
    }
    n_mem += n_phiAO;

    // required memory for grids x, y, and z and weights w vectors
    n_grids = 4 * phi_points_;
    n_mem += n_grids;

    // memory needed for AO->MO transformation
    n_transf = nirrep_ * phi_points_;                     // phi_points_list vector
    n_transf += phi_points_ * nso_;                       // super_phi_ matrix
    n_transf += nso_ * nso_;                              // AO2SO matrix in TransformPhiMatrixAOMO()
    n_transf += phi_points_ * nso_;                       // temporary matrix called three times in TransformPhiMatrixAOMO()
    if ( is_gga_ || is_meta_ ) { 
        n_transf += 3 * phi_points_ * nso_;               // for super_phi_x_, _y_ and _z_ gradient matrices
    }
    n_mem += n_transf;

    // memory needed for storing rho and rho' matrices
    n_rho  = 3 * phi_points_;                             // rho_a_, rho_b_ and rho_ vectors
    if ( is_gga_ || is_meta_ ) { 
        n_rho += 3 * phi_points_;                         // rho_a_x_, rho_a_y_ and rho_a_z_ gradient vectors
        n_rho += 3 * phi_points_;                         // rho_b_x_, rho_b_y_ and rho_b_z_ gradient vectors
        n_rho += 3 * phi_points_;                         // sigma_aa_, sigma_ab_ and sigma_bb_ vectors
    }
    n_mem += n_rho;

    // memory needed for storing pi vector
    n_pi   = phi_points_;                                 // pi_ vector
    if ( is_gga_ || is_meta_ ) { 
        n_pi  += 3 * phi_points_;                         // pi_x_, pi_y_ and pi_z_ gradient vectors
    }
    n_mem += n_pi;
 
    // memory needed for storing on-top ratio vector
    n_R    = phi_points_;
    n_mem += n_R;
    
    // memory needed for (full-) translation step
    n_transl  = 2 * phi_points_;                          // (f)tr_rho_a_, (f)tr_rho_b_ vectors
    if ( is_gga_ || is_meta_ ) { 
        n_transl += 3 * phi_points_;                      // (f)tr_rho_a_x_, (f)tr_rho_a_y_ and (f)tr_rho_a_z_ gradient vectors
        n_transl += 3 * phi_points_;                      // (f)tr_rho_b_x_, (f)tr_rho_b_y_ and (f)tr_rho_b_z_ gradient vectors
        n_transl += 3 * phi_points_;                      // (f)tr_sigma_aa_, (f)tr_sigma_ab_ and (f)tr_sigma_bb_ vectors
    }
    n_mem += n_transl;

    if ( n_mem * 8.0 > (double)memory_ ) {
        outfile->Printf("\n"); 
        outfile->Printf("      <<< Warning >>> \n");
        outfile->Printf("\n");
        outfile->Printf("      Not enough memory. Switching to low-memory algorithm.\n");
        outfile->Printf("      For optimal performance, increase the available memory\n");
        outfile->Printf("      by at least %7.2lf mb.\n",(8.0 * n_mem - memory_)/1024.0/1024.0);
        outfile->Printf("\n");

        is_low_memory_ = true;
        n_mem -= n_transf;
        n_mem -= n_phiAO;
        if ( n_mem * 8.0 > (double)memory_ ) {
            outfile->Printf("\n");
            outfile->Printf("      Wow.  On second thought, there isn't even enough memory for the low-memory algorithm!\n");
            outfile->Printf("\n");
            throw PsiException("Not enough memory.",__FILE__,__LINE__); fflush(stdout);
        }
    }

    outfile->Printf("    =========================================================\n");
    outfile->Printf("    memory specified by the user:                %7.2lf MB\n",(double)memory_ / 1024.0 / 1024.0);
    outfile->Printf("    ---------------------------------------------------------\n");
    outfile->Printf("    initialization:                              %7.2lf MB\n",n_init  * sizeof(double)   / 1024.0 / 1024.0);
    outfile->Printf("    grid-points and weights:                     %7.2lf MB\n",n_grids * sizeof(double)   / 1024.0 / 1024.0);
    if ( !is_low_memory_ ) {
        outfile->Printf("    phi & phi' (AO):                             %7.2lf MB\n",n_phiAO * sizeof(double)   / 1024.0 / 1024.0);
        outfile->Printf("    AO->SO transformation:                       %7.2lf MB\n",n_transf* sizeof(double)   / 1024.0 / 1024.0);
    }
    outfile->Printf("    rho and rho':                                %7.2lf MB\n",n_rho *    sizeof(double)   / 1024.0 / 1024.0);
    outfile->Printf("    on-top pair density:                         %7.2lf MB\n",n_pi  *    sizeof(double)   / 1024.0 / 1024.0);
    outfile->Printf("    on-top ratio:                                %7.2lf MB\n",n_R   *    sizeof(double)   / 1024.0 / 1024.0);
    outfile->Printf("    translation step:                            %7.2lf MB\n",n_transl * sizeof(double)   / 1024.0 / 1024.0);
    outfile->Printf("    extra:                                       %7.2lf MB\n",n_extra *  sizeof(double)   / 1024.0 / 1024.0);
    outfile->Printf("    ---------------------------------------------------------\n");
    outfile->Printf("    total memory for MCPDFT:                     %7.2lf MB\n",n_mem *    sizeof(double)   / 1024.0 / 1024.0);
    outfile->Printf("    =========================================================\n");

    // memory available after allocating all we need for MCPDFT (deleting the temporary matrices and vectors)
    n_mem -= nirrep_ * phi_points_;                // phi_points_list vector
    n_mem -= nso_ * nso_;                          // AO2SO matrix in TransformPhiMatrixAOMO()
    if ( !is_low_memory_ ) {
        n_mem -= n_phiAO;                          // ao-basis phi matrices
        n_mem -= phi_points_ * nso_;               // temporary matrix called three times in TransformPhiMatrixAOMO()
    }

    available_memory_ = memory_ - n_mem * 8L;

    outfile->Printf("\n");
    outfile->Printf("    available memeory for building JK object:    %7.2lf MB\n",(double)available_memory_ / 1024.0 / 1024.0);
    outfile->Printf("\n");

    grid_x_      = std::shared_ptr<Vector>(new Vector("GRID X",phi_points_));
    grid_y_      = std::shared_ptr<Vector>(new Vector("GRID Y",phi_points_));
    grid_z_      = std::shared_ptr<Vector>(new Vector("GRID Z",phi_points_));
    grid_w_      = std::shared_ptr<Vector>(new Vector("GRID W",phi_points_));

    if ( is_low_memory_ ) {

        GetGridInfo();

    }else {
        outfile->Printf("\n"); 
        outfile->Printf("    ==> Build Phi and Phi' matrices ..."); fflush(stdout);

        std::shared_ptr<Matrix> super_phi_ao = std::shared_ptr<Matrix>(new Matrix("SUPER PHI (AO)",phi_points_,nso_));

        std::shared_ptr<Matrix> super_phi_x_ao;
        std::shared_ptr<Matrix> super_phi_y_ao;
        std::shared_ptr<Matrix> super_phi_z_ao;

        if ( is_gga_ || is_meta_ ) {

            super_phi_x_ao = std::shared_ptr<Matrix>(new Matrix("SUPER PHI X (AO)",phi_points_,nso_));
            super_phi_y_ao = std::shared_ptr<Matrix>(new Matrix("SUPER PHI Y (AO)",phi_points_,nso_));
            super_phi_z_ao = std::shared_ptr<Matrix>(new Matrix("SUPER PHI Z (AO)",phi_points_,nso_));
            
        }

        // build phi matrix and derivative phi matrices (ao basis)
        BuildPhiMatrixAO("PHI",super_phi_ao);

        if ( is_gga_ || is_meta_ ) {

            BuildPhiMatrixAO("PHI_X",super_phi_x_ao);
            BuildPhiMatrixAO("PHI_Y",super_phi_y_ao);
            BuildPhiMatrixAO("PHI_Z",super_phi_z_ao);

        }

        // transform the orbital label in phi matrix and derivative phi matrix to the MO basis

        Dimension phi_points_list = Dimension(nirrep_,"phi_points");
        for (int h = 0; h < nirrep_; h++) {
            phi_points_list[h] = phi_points_;
        }

        super_phi_ = std::shared_ptr<Matrix>(new Matrix("SUPER PHI",phi_points_list,nsopi_));

        TransformPhiMatrixAOMO(super_phi_ao,super_phi_);

        if ( is_gga_ || is_meta_ ) {

            super_phi_x_ = std::shared_ptr<Matrix>(new Matrix("SUPER PHI X",phi_points_list,nsopi_));
            super_phi_y_ = std::shared_ptr<Matrix>(new Matrix("SUPER PHI Y",phi_points_list,nsopi_));
            super_phi_z_ = std::shared_ptr<Matrix>(new Matrix("SUPER PHI Z",phi_points_list,nsopi_));

            TransformPhiMatrixAOMO(super_phi_x_ao,super_phi_x_);
            TransformPhiMatrixAOMO(super_phi_y_ao,super_phi_y_);
            TransformPhiMatrixAOMO(super_phi_z_ao,super_phi_z_);

        }
        outfile->Printf(" Done. <==\n"); 
    }

}// end of common_init()

void MCPDFTSolver::GetGridInfo() {

    int nblocks = potential_->nblocks();

    phi_points_ = 0;
    for (int myblock = 0; myblock < nblocks; myblock++) {
        std::shared_ptr<BlockOPoints> block = potential_->get_block(myblock);
        points_func_->compute_points(block);
        int npoints = block->npoints();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        double ** phi = points_func_->basis_value("PHI")->pointer();

        double * x = block->x();
        double * y = block->y();
        double * z = block->z();
        double * w = block->w();

        for (int p = 0; p < npoints; p++) {

            grid_x_->pointer()[phi_points_ + p] = x[p];
            grid_y_->pointer()[phi_points_ + p] = y[p];
            grid_z_->pointer()[phi_points_ + p] = z[p];
            grid_w_->pointer()[phi_points_ + p] = w[p];

        }
        phi_points_ += npoints;
    }
}//end of GetGridInfo()

void MCPDFTSolver::BuildPhiMatrixAO(std::string phi_type, std::shared_ptr<Matrix> myphi) {

    int nblocks = potential_->nblocks();

    phi_points_ = 0;
    for (int myblock = 0; myblock < nblocks; myblock++) {
        std::shared_ptr<BlockOPoints> block = potential_->get_block(myblock);
        points_func_->compute_points(block);
        int npoints = block->npoints();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        double ** phi = points_func_->basis_value(phi_type)->pointer();

        double * x = block->x();
        double * y = block->y();
        double * z = block->z();
        double * w = block->w();

        // Doing some test to see everything including Pi etc is correct on 
        // Molcas' grid points through comparison.

        // std::ifstream dataIn;
       
        // dataIn.open("H2.grids_test");
        // 
        // if (!dataIn)
        //    std::cout << "Error opening file.\n";
        // else { 
        //      int p = 0;        
        //      while (!dataIn.eof()){
        //    
        //            dataIn >> x[p];
        //            dataIn >> y[p];    
        //            dataIn >> z[p];
        //            p++;
        //      }        
        // }
        // dataIn.close(); 

        // for (int p = 0; p < npoints; p++) {

        //     outfile->Printf("\n     y[");
        //     outfile->Printf("%d",p);
        //     outfile->Printf("] = %20.7lf\n",y[p]);
        // }

        for (int p = 0; p < npoints; p++) {

            grid_x_->pointer()[phi_points_ + p] = x[p];
            grid_y_->pointer()[phi_points_ + p] = y[p];
            grid_z_->pointer()[phi_points_ + p] = z[p];
            grid_w_->pointer()[phi_points_ + p] = w[p];

            for (int nu = 0; nu < nlocal; nu++) {
                int nug = function_map[nu];
                myphi->pointer()[phi_points_ + p][nug] = phi[p][nu];
            }
        }
        phi_points_ += npoints;
    }
}

void MCPDFTSolver::TransformPhiMatrixAOMO(std::shared_ptr<Matrix> phi_in, std::shared_ptr<Matrix> phi_out) {
        
    std::shared_ptr<Matrix> phi_temp (new Matrix(phi_out));

    // grab AO->SO transformation matrix
    std::shared_ptr<Matrix> ao2so = reference_wavefunction_->aotoso();

    for (int h = 0; h < nirrep_; h++) {
        if (nsopi_[h] != 0 )
           F_DGEMM('n','n',nsopi_[h],phi_points_,nso_,1.0,ao2so->pointer(h)[0],nsopi_[h],phi_in->pointer()[0],nso_,0.0,phi_temp->pointer(h)[0],nsopi_[h]);
    }

    for (int h = 0; h < nirrep_; h++) {
        if (nsopi_[h] != 0 )
           F_DGEMM('n','n',nsopi_[h],phi_points_,nsopi_[h],1.0,Ca_->pointer(h)[0],nsopi_[h],phi_temp->pointer(h)[0],nsopi_[h],0.0,phi_out->pointer(h)[0],nsopi_[h]);
    }
}

double MCPDFTSolver::compute_energy() {
    
    /* ========================================================================================
     * (1) Calculation of the range-separated global (double-)hybrid formalis is based
     * on Eq. (12) of the reference: C. Kalai and J, Toulouse J. Chem. Phys. 148, 164105 (2018)
     * (2) Calculation of linearly scaled 1-parameter double-hybrid (LS1DH) is based upon Eq.
     * (10) of the manuscript: Toulouse J. et al., J. Chem. Phys. 135, 101102 (2011)
       ======================================================================================== */

    // read 1- and 2-RDM from disk and build rho(r), rho'(r), pi(r), and pi'(r)

    if ( options_.get_str("MCPDFT_REFERENCE") == "V2RDM" ) {

        ReadOPDM();
        ReadTPDM();

    }else if ( options_.get_str("MCPDFT_REFERENCE") == "CI" ) {

        if ( nirrep_ > 1 ) {
            throw PsiException("error, MCPDFT_REFERENCE = CI only works with symmetry c1",__FILE__,__LINE__);
        }

        // read 1-RDM
        ReadCIOPDM(Da_,"opdm_a.txt");
        ReadCIOPDM(Db_,"opdm_b.txt");

        // allocate memory 2-RDM
        double * D2ab = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
        memset((void*)D2ab,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

        // read 2-RDM
        ReadCITPDM(D2ab,"tpdm_ab.txt");

        // build alpha- and beta-spin densities and gradients (already built for MCPDFT_REFERENCE = V2RDM)
        outfile->Printf("\n"); 
        outfile->Printf("    ==> Build Rho <== \n ");

        BuildRho();

        // build on-top pair density (already built for MCPDFT_REFERENCE = V2RDM)
        outfile->Printf("\n");
        outfile->Printf("    ==> Build Pi <==\n\n");

        BuildPi(D2ab);
        free(D2ab);

    }else {
        throw PsiException("invalid MCPDFT_REFERENCE type",__FILE__,__LINE__);
    }

    // build R(r) = 4 * Pi(r) / rho(r)
    outfile->Printf("    ==> Build the on-top ratio R ...");
    Build_R();
    outfile->Printf(" Done. <==\n\n");

    // translate the alpha and beta densities and their corresponding gradients
    if ( options_.get_str("MCPDFT_TRANSLATION_TYPE") == "REGULAR") {

        outfile->Printf("    ==> Regular translation of densities and/or density gradients ...\n");
        Translate();    
        outfile->Printf("    ... Done. <==\n\n");

    }else {

        outfile->Printf("    ==> Full translation of densities and/or density gradients ...\n");
        Fully_Translate();    
        outfile->Printf("    ... Done. <==\n\n");

    }

    // calculate the on-top energy

    outfile->Printf("    ==> Evaluate the on-top energy contribution <==\n");
    outfile->Printf("\n");

    /* ===================================================================================================
       calculate the complement short-range MCPDFT XC functional energy:
       E = min(Psi->N) { <Psi| T + Wee_LR(w) + lambda * Wee_SR(w) + Vne |Psi> + E_HXC_(w,lambda)[rho,pi] }
       =================================================================================================== */

    // computing the exchange and correlation density functional energies
    double mcpdft_ex = 0.0;
    double mcpdft_ec = 0.0;
    if (options_.get_str("MCPDFT_METHOD") == "MCPDFT") {

         if ( options_.get_str("MCPDFT_FUNCTIONAL") == "SVWN" ) {

             mcpdft_ex = EX_LSDA(tr_rho_a_, tr_rho_b_);
             mcpdft_ec = EC_VWN3_RPA_III(tr_rho_a_, tr_rho_b_);

         }else if ( options_.get_str("MCPDFT_FUNCTIONAL") == "BLYP" ) {

             mcpdft_ex = EX_B88_I(tr_rho_a_, tr_rho_b_, tr_sigma_aa_, tr_sigma_bb_);
             mcpdft_ec = EC_LYP_I(tr_rho_a_, tr_rho_b_, tr_sigma_aa_, tr_sigma_ab_, tr_sigma_bb_); 

         }else if (  options_.get_str("MCPDFT_FUNCTIONAL") == "PBE" 
                  || options_.get_str("MCPDFT_FUNCTIONAL") == "REVPBE" ) {
 
                  mcpdft_ex = EX_PBE_I(tr_rho_a_, tr_rho_b_, tr_sigma_aa_, tr_sigma_bb_);
                  mcpdft_ec = EC_PBE_I(tr_rho_a_, tr_rho_b_, tr_sigma_aa_, tr_sigma_ab_, tr_sigma_bb_);

         }else if ( options_.get_str("MCPDFT_FUNCTIONAL") == "BOP" ) {
            
                  mcpdft_ex = EX_B88_I(tr_rho_a_, tr_rho_b_, tr_sigma_aa_, tr_sigma_bb_);
                  mcpdft_ec = EC_B88_OP(tr_rho_a_, tr_rho_b_, tr_sigma_aa_, tr_sigma_bb_);
         }else{
              throw PsiException("only SVWN, BLYP, PBE, REVPBE, and BOP are currecntly supported in MCPDFT",__FILE__,__LINE__);
         }
    }else {
          throw PsiException("Please choose a valid method for MCDPFT",__FILE__,__LINE__);
    }
    
    // evaluate the kinetic, potential, and coulomb energies
    outfile->Printf("    ==> Evaluate kinetic, potential, and coulomb energies <==\n");
    outfile->Printf("\n");

    // one-electron terms:
    std::shared_ptr<MintsHelper> mints(new MintsHelper(reference_wavefunction_));

    // kinetic-energy 
    SharedMatrix Ta (new Matrix(mints->so_kinetic()));
    Ta->transform(Ca_);
    
    SharedMatrix Tb (new Matrix(mints->so_kinetic()));
    Tb->transform(Cb_);
    
    double kinetic_energy = Da_->vector_dot(Ta) 
                          + Db_->vector_dot(Tb);

    // nuclear-attraction potential energy
    SharedMatrix Va (new Matrix(mints->so_potential()));
    Va->transform(Ca_);
    
    SharedMatrix Vb (new Matrix(mints->so_potential()));
    Vb->transform(Cb_);
    
    double en_potential_energy = Da_->vector_dot(Va) 
                               + Db_->vector_dot(Vb);
    
    // classical nuclear repulsion energy
    double nuclear_repulsion_energy = molecule_->nuclear_repulsion_energy({0.0,0.0,0.0});

    // print total energy and its components
    outfile->Printf("    ==> Energetics <==\n");
    outfile->Printf("\n");

    outfile->Printf("        nuclear repulsion energy =          %20.12lf\n",nuclear_repulsion_energy);
    outfile->Printf("        electron-nucleus potential energy = %20.12lf\n",en_potential_energy);
    outfile->Printf("        electron kinetic energy =           %20.12lf\n",kinetic_energy);
    outfile->Printf("        classical coulomb energy  =         %20.12lf\n",coulomb_energy_);
    outfile->Printf("        Ex                        =         %20.12lf\n",mcpdft_ex);
    outfile->Printf("        Ec                        =         %20.12lf\n",mcpdft_ec);
    outfile->Printf("\n");

    double total_energy = en_potential_energy + kinetic_energy + coulomb_energy_ + mcpdft_ex + mcpdft_ec + nuclear_repulsion_energy;

    if( options_.get_str("MCPDFT_METHOD") == "MCPDFT") {
      outfile->Printf("    * MCPDFT total energy   =     ");
    } 
    outfile->Printf("   %20.12lf\n\n",total_energy);

    Process::environment.globals["CURRENT ENERGY"]     = total_energy;
    Process::environment.globals["MCPDFT TOTAL ENERGY"] = total_energy;

    return total_energy;
}

void MCPDFTSolver::BuildPi(double * D2ab) {

    pi_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double * pi_p = pi_->pointer();

    double ** phi = super_phi_->pointer();

    for (int p = 0; p < phi_points_; p++) {

        double dum = 0.0;

        // pi(r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
        for (int mu = 0; mu < nmo_; mu++) {
            for (int nu = 0; nu < nmo_; nu++) {
                for (int lambda = 0; lambda < nmo_; lambda++) {
                    for (int sigma = 0; sigma < nmo_; sigma++) {

                        dum += phi[p][mu] * phi[p][lambda] * phi[p][sigma] * phi[p][nu] * D2ab[mu*nmo_*nmo_*nmo_ + nu*nmo_*nmo_ + lambda*nmo_ + sigma];
                    }
                }
            }
        }

        pi_p[p] = dum;
    }

    if ( is_gga_ || is_meta_ ) {

       pi_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       pi_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       pi_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       
       double * pi_xp = pi_x_->pointer();
       double * pi_yp = pi_y_->pointer();
       double * pi_zp = pi_z_->pointer();
       
       double ** phi_x = super_phi_x_->pointer();
       double ** phi_y = super_phi_y_->pointer();
       double ** phi_z = super_phi_z_->pointer();

       for (int p = 0; p < phi_points_; p++) {

           double dum_x = 0.0;
           double dum_y = 0.0;
           double dum_z = 0.0;

           // pi(r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
           for (int mu = 0; mu < nmo_; mu++) {
               for (int nu = 0; nu < nmo_; nu++) {
                   for (int lambda = 0; lambda < nmo_; lambda++) {
                       for (int sigma = 0; sigma < nmo_; sigma++) {

                           dum_x += ( phi_x[p][mu] * phi[p][lambda] * phi[p][sigma] * phi[p][nu] +
                                      phi[p][mu] * phi_x[p][lambda] * phi[p][sigma] * phi[p][nu] +
                                      phi[p][mu] * phi[p][lambda] * phi_x[p][sigma] * phi[p][nu] +
                                      phi[p][mu] * phi[p][lambda] * phi[p][sigma] * phi_x[p][nu] ) * D2ab[mu*nmo_*nmo_*nmo_ + nu*nmo_*nmo_ + lambda*nmo_ + sigma];

                           dum_y += ( phi_y[p][mu] * phi[p][lambda] * phi[p][sigma] * phi[p][nu] +
                                      phi[p][mu] * phi_y[p][lambda] * phi[p][sigma] * phi[p][nu] +
                                      phi[p][mu] * phi[p][lambda] * phi_y[p][sigma] * phi[p][nu] +
                                      phi[p][mu] * phi[p][lambda] * phi[p][sigma] * phi_y[p][nu] ) * D2ab[mu*nmo_*nmo_*nmo_ + nu*nmo_*nmo_ + lambda*nmo_ + sigma];


                           dum_z += ( phi_z[p][mu] * phi[p][lambda] * phi[p][sigma] * phi[p][nu] +
                                      phi[p][mu] * phi_z[p][lambda] * phi[p][sigma] * phi[p][nu] +
                                      phi[p][mu] * phi[p][lambda] * phi_z[p][sigma] * phi[p][nu] +
                                      phi[p][mu] * phi[p][lambda] * phi[p][sigma] * phi_z[p][nu] ) * D2ab[mu*nmo_*nmo_*nmo_ + nu*nmo_*nmo_ + lambda*nmo_ + sigma];
                       }
                   }
               }
           }

           pi_xp[p] = dum_x;
           pi_yp[p] = dum_y;
           pi_zp[p] = dum_z;

       }
    }
}

void MCPDFTSolver::BuildPiLowMemory(tpdm * D2ab, int nab) {

    // grab AO->SO transformation matrix
    std::shared_ptr<Matrix> ao2so = reference_wavefunction_->aotoso();

    int nblocks = potential_->nblocks();

    Dimension phi_points_list = Dimension(nirrep_,"phi_points");
    for (int h = 0; h < nirrep_; h++) {
        phi_points_list[h] = max_points_;
    }
    std::shared_ptr<Matrix> myPhiAO(new Matrix(max_points_,nso_));
    std::shared_ptr<Matrix> myPhiMO(new Matrix(phi_points_list,nsopi_));
    std::shared_ptr<Matrix> tempPhi(new Matrix(phi_points_list,nsopi_));

    std::shared_ptr<Matrix> myPhiMO_x;
    std::shared_ptr<Matrix> myPhiMO_y;
    std::shared_ptr<Matrix> myPhiMO_z;

    // allocate memory for Pi, Pi_x, etc.
    pi_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    double * pi_p  = pi_->pointer();

    double * pi_xp;
    double * pi_yp;
    double * pi_zp;

    if ( is_gga_ || is_meta_ ) {

        myPhiMO_x = (std::shared_ptr<Matrix>)(new Matrix(phi_points_list,nsopi_));
        myPhiMO_y = (std::shared_ptr<Matrix>)(new Matrix(phi_points_list,nsopi_));
        myPhiMO_z = (std::shared_ptr<Matrix>)(new Matrix(phi_points_list,nsopi_));

        pi_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        pi_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        pi_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        
        pi_xp = pi_x_->pointer();
        pi_yp = pi_y_->pointer();
        pi_zp = pi_z_->pointer();
    }

    // memory for rho
    rho_a_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    rho_b_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    rho_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();
    double * rho_p  = rho_->pointer();

    double * sigma_aap;
    double * sigma_bbp;
    double * sigma_abp;

    double * rho_a_xp;
    double * rho_b_xp;

    double * rho_a_yp;
    double * rho_b_yp;

    double * rho_a_zp;
    double * rho_b_zp;

    if ( is_gga_ || is_meta_ ) {

        sigma_aa_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        sigma_bb_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        sigma_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        sigma_aap = sigma_aa_->pointer();
        sigma_bbp = sigma_bb_->pointer();
        sigma_abp = sigma_ab_->pointer();

        rho_a_xp = rho_a_x_->pointer();
        rho_b_xp = rho_b_x_->pointer();

        rho_a_yp = rho_a_y_->pointer();
        rho_b_yp = rho_b_y_->pointer();

        rho_a_zp = rho_a_z_->pointer();
        rho_b_zp = rho_b_z_->pointer();

    }

    // number of nonzero elements in Da and Db:
    long int nb;
    long int na;

    std::shared_ptr<PSIO> psio (new PSIO());

    psio->open(PSIF_V2RDM_D1A,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_V2RDM_D1A,"length",(char*)&na,sizeof(long int));
    psio->close(PSIF_V2RDM_D1A,1);

    psio->open(PSIF_V2RDM_D1B,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_V2RDM_D1B,"length",(char*)&nb,sizeof(long int));
    psio->close(PSIF_V2RDM_D1B,1);
        
    phi_points_ = 0;
    for (int myblock = 0; myblock < nblocks; myblock++) {
        std::shared_ptr<BlockOPoints> block = potential_->get_block(myblock);
        points_func_->compute_points(block);
        int npoints = block->npoints();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        double ** phi   = points_func_->basis_value("PHI")->pointer();

        for (int p = 0; p < npoints; p++) {

            // grab block of AO-basis phi matrix
            myPhiAO->zero();
            for (int nu = 0; nu < nlocal; nu++) {
                int nug = function_map[nu];
                myPhiAO->pointer()[p][nug] = phi[p][nu];
            }

            // transform AO-basis block of Phi matrix to MO basis
            for (int h = 0; h < nirrep_; h++) {
                if (nsopi_[h] != 0 ) {
                    F_DGEMM('n','n',nsopi_[h],npoints,nso_,1.0,ao2so->pointer(h)[0],nsopi_[h],myPhiAO->pointer()[0],nso_,0.0,tempPhi->pointer(h)[0],nsopi_[h]);
                }
            }

            for (int h = 0; h < nirrep_; h++) {
                if (nsopi_[h] != 0 ) {
                    F_DGEMM('n','n',nsopi_[h],npoints,nsopi_[h],1.0,Ca_->pointer(h)[0],nsopi_[h],tempPhi->pointer(h)[0],nsopi_[h],0.0,myPhiMO->pointer(h)[0],nsopi_[h]);
                }
            }

            // rho_a(r)
            double duma = 0.0;
            for (int n = 0; n < na; n++) {
                
                int i = opdm_a_[n].i;
                int j = opdm_a_[n].j;
                
                int hi = symmetry_[i];
                int hj = symmetry_[j];
                
                int ii = i - pitzer_offset_[hi];
                int jj = j - pitzer_offset_[hj];
                
                duma += myPhiMO->pointer(hi)[p][ii] * 
                        myPhiMO->pointer(hj)[p][jj] * opdm_a_[n].value;
            
            }
            rho_ap[phi_points_ + p] = duma;

            // rho_b(r)
            double dumb = 0.0;
            for (int n = 0; n < nb; n++) {

                int i = opdm_b_[n].i;
                int j = opdm_b_[n].j;

                int hi = symmetry_[i];
                int hj = symmetry_[j];

                int ii = i - pitzer_offset_[hi];
                int jj = j - pitzer_offset_[hj];

                dumb += myPhiMO->pointer(hi)[p][ii] *
                        myPhiMO->pointer(hj)[p][jj] * opdm_b_[n].value;

            }
            rho_bp[phi_points_ + p] = dumb;

            rho_p[phi_points_ + p] = rho_ap[phi_points_ + p] + rho_bp[phi_points_ + p];


            // pi(r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
            double dum = 0.0;
            for (int n = 0; n < nab; n++) {

                int i = D2ab[n].i;
                int j = D2ab[n].j;
                int k = D2ab[n].k;
                int l = D2ab[n].l;

                int hi = symmetry_[i];
                int hj = symmetry_[j];
                int hk = symmetry_[k];
                int hl = symmetry_[l];

                int ii = i - pitzer_offset_[hi];
                int jj = j - pitzer_offset_[hj];
                int kk = k - pitzer_offset_[hk];
                int ll = l - pitzer_offset_[hl];

                dum += myPhiMO->pointer(hi)[p][ii] *
                       myPhiMO->pointer(hj)[p][jj] *
                       myPhiMO->pointer(hk)[p][kk] *
                       myPhiMO->pointer(hl)[p][ll] * D2ab[n].value;

            }

            pi_p[phi_points_ + p] = dum;

            if ( is_gga_ || is_meta_ ) {

                double ** phi_x = points_func_->basis_value("PHI_X")->pointer();
                double ** phi_y = points_func_->basis_value("PHI_Y")->pointer();
                double ** phi_z = points_func_->basis_value("PHI_Z")->pointer();

                // grab block of AO-basis phi_x matrix
                myPhiAO->zero();
                for (int nu = 0; nu < nlocal; nu++) {
                    int nug = function_map[nu];
                    myPhiAO->pointer()[p][nug] = phi_x[p][nu];
                }

                // transform AO-basis block of Phi_x matrix to MO basis
                for (int h = 0; h < nirrep_; h++) {
                    if (nsopi_[h] != 0 ) {
                        F_DGEMM('n','n',nsopi_[h],npoints,nso_,1.0,ao2so->pointer(h)[0],nsopi_[h],myPhiAO->pointer()[0],nso_,0.0,tempPhi->pointer(h)[0],nsopi_[h]);
                    }
                }

                for (int h = 0; h < nirrep_; h++) {
                    if (nsopi_[h] != 0 ) {
                        F_DGEMM('n','n',nsopi_[h],npoints,nsopi_[h],1.0,Ca_->pointer(h)[0],nsopi_[h],tempPhi->pointer(h)[0],nsopi_[h],0.0,myPhiMO_x->pointer(h)[0],nsopi_[h]);
                    }
                }

                // grab block of AO-basis phi_y matrix
                myPhiAO->zero();
                for (int nu = 0; nu < nlocal; nu++) {
                    int nug = function_map[nu];
                    myPhiAO->pointer()[p][nug] = phi_y[p][nu];
                }

                // transform AO-basis block of Phi_y matrix to MO basis
                for (int h = 0; h < nirrep_; h++) {
                    if (nsopi_[h] != 0 ) {
                        F_DGEMM('n','n',nsopi_[h],npoints,nso_,1.0,ao2so->pointer(h)[0],nsopi_[h],myPhiAO->pointer()[0],nso_,0.0,tempPhi->pointer(h)[0],nsopi_[h]);
                    }
                }

                for (int h = 0; h < nirrep_; h++) {
                    if (nsopi_[h] != 0 ) {
                        F_DGEMM('n','n',nsopi_[h],npoints,nsopi_[h],1.0,Ca_->pointer(h)[0],nsopi_[h],tempPhi->pointer(h)[0],nsopi_[h],0.0,myPhiMO_y->pointer(h)[0],nsopi_[h]);
                    }
                }
                // grab block of AO-basis phi_z matrix
                myPhiAO->zero();
                for (int nu = 0; nu < nlocal; nu++) {
                    int nug = function_map[nu];
                    myPhiAO->pointer()[p][nug] = phi_z[p][nu];
                }

                // transform AO-basis block of Phi_z matrix to MO basis
                for (int h = 0; h < nirrep_; h++) {
                    if (nsopi_[h] != 0 ) {
                        F_DGEMM('n','n',nsopi_[h],npoints,nso_,1.0,ao2so->pointer(h)[0],nsopi_[h],myPhiAO->pointer()[0],nso_,0.0,tempPhi->pointer(h)[0],nsopi_[h]);
                    }
                }

                for (int h = 0; h < nirrep_; h++) {
                    if (nsopi_[h] != 0 ) {
                        F_DGEMM('n','n',nsopi_[h],npoints,nsopi_[h],1.0,Ca_->pointer(h)[0],nsopi_[h],tempPhi->pointer(h)[0],nsopi_[h],0.0,myPhiMO_z->pointer(h)[0],nsopi_[h]);
                    }
                }

                // rho'_a(r)

                double duma_x = 0.0;
                double duma_y = 0.0;
                double duma_z = 0.0;

                for (int n = 0; n < na; n++) {

                    int i = opdm_a_[n].i;
                    int j = opdm_a_[n].j;

                    int hi = symmetry_[i];
                    int hj = symmetry_[j];

                    int ii = i - pitzer_offset_[hi];
                    int jj = j - pitzer_offset_[hj];

                    duma_x += ( myPhiMO_x->pointer(hi)[p][ii] * myPhiMO->pointer(hj)[p][jj]
                            +   myPhiMO->pointer(hi)[p][ii] * myPhiMO_x->pointer(hj)[p][jj] ) * opdm_a_[n].value;

                    duma_y += ( myPhiMO_y->pointer(hi)[p][ii] * myPhiMO->pointer(hj)[p][jj]
                            +   myPhiMO->pointer(hi)[p][ii] * myPhiMO_y->pointer(hj)[p][jj] ) * opdm_a_[n].value;

                    duma_z += ( myPhiMO_z->pointer(hi)[p][ii] * myPhiMO->pointer(hj)[p][jj]
                            +   myPhiMO->pointer(hi)[p][ii] * myPhiMO_z->pointer(hj)[p][jj] ) * opdm_a_[n].value;

                }

                // rho'_b(r)

                double dumb_x = 0.0;
                double dumb_y = 0.0;
                double dumb_z = 0.0;

                for (int n = 0; n < nb; n++) {

                    int i = opdm_b_[n].i;
                    int j = opdm_b_[n].j;

                    int hi = symmetry_[i];
                    int hj = symmetry_[j];

                    int ii = i - pitzer_offset_[hi];
                    int jj = j - pitzer_offset_[hj];

                    dumb_x += ( myPhiMO_x->pointer(hi)[p][ii] * myPhiMO->pointer(hj)[p][jj]
                            +   myPhiMO->pointer(hi)[p][ii] * myPhiMO_x->pointer(hj)[p][jj] ) * opdm_b_[n].value;

                    dumb_y += ( myPhiMO_y->pointer(hi)[p][ii] * myPhiMO->pointer(hj)[p][jj]
                            +   myPhiMO->pointer(hi)[p][ii] * myPhiMO_y->pointer(hj)[p][jj] ) * opdm_b_[n].value;

                    dumb_z += ( myPhiMO_z->pointer(hi)[p][ii] * myPhiMO->pointer(hj)[p][jj]
                            +   myPhiMO->pointer(hi)[p][ii] * myPhiMO_z->pointer(hj)[p][jj] ) * opdm_b_[n].value;

                }

                rho_a_xp[phi_points_ + p] = duma_x;
                rho_b_xp[phi_points_ + p] = dumb_x;

                rho_a_yp[phi_points_ + p] = duma_y;
                rho_b_yp[phi_points_ + p] = dumb_y;

                rho_a_zp[phi_points_ + p] = duma_z;
                rho_b_zp[phi_points_ + p] = dumb_z;

                sigma_aap[phi_points_ + p] = ( rho_a_xp[phi_points_ + p] * rho_a_xp[phi_points_ + p] ) 
                                           + ( rho_a_yp[phi_points_ + p] * rho_a_yp[phi_points_ + p] )
                                           + ( rho_a_zp[phi_points_ + p] * rho_a_zp[phi_points_ + p] );
                sigma_bbp[phi_points_ + p] = ( rho_b_xp[phi_points_ + p] * rho_b_xp[phi_points_ + p] ) 
                                           + ( rho_b_yp[phi_points_ + p] * rho_b_yp[phi_points_ + p] ) 
                                           + ( rho_b_zp[phi_points_ + p] * rho_b_zp[phi_points_ + p] );
                sigma_abp[phi_points_ + p] = ( rho_a_xp[phi_points_ + p] * rho_b_xp[phi_points_ + p] ) 
                                           + ( rho_a_yp[phi_points_ + p] * rho_b_yp[phi_points_ + p] ) 
                                           + ( rho_a_zp[phi_points_ + p] * rho_b_zp[phi_points_ + p] );

                double dum_x = 0.0;
                double dum_y = 0.0;
                double dum_z = 0.0;

                // pi(r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
                for (int n = 0; n < nab; n++) {

                    int i = D2ab[n].i;
                    int j = D2ab[n].j;
                    int k = D2ab[n].k;
                    int l = D2ab[n].l;
                    
                    int hi = symmetry_[i];
                    int hj = symmetry_[j];
                    int hk = symmetry_[k]; 
                    int hl = symmetry_[l];
                    
                    int ii = i - pitzer_offset_[hi];
                    int jj = j - pitzer_offset_[hj];
                    int kk = k - pitzer_offset_[hk];
                    int ll = l - pitzer_offset_[hl];

                    dum_x += ( myPhiMO_x->pointer(hi)[p][ii] *
                               myPhiMO->pointer(hj)[p][jj] *
                               myPhiMO->pointer(hk)[p][kk] * 
                               myPhiMO->pointer(hl)[p][ll] +

                               myPhiMO->pointer(hi)[p][ii] *
                               myPhiMO_x->pointer(hj)[p][jj] *
                               myPhiMO->pointer(hk)[p][kk] * 
                               myPhiMO->pointer(hl)[p][ll] +

                               myPhiMO->pointer(hi)[p][ii] *
                               myPhiMO->pointer(hj)[p][jj] *
                               myPhiMO_x->pointer(hk)[p][kk] * 
                               myPhiMO->pointer(hl)[p][ll] +

                               myPhiMO->pointer(hi)[p][ii] *
                               myPhiMO->pointer(hj)[p][jj] *
                               myPhiMO->pointer(hk)[p][kk] * 
                               myPhiMO_x->pointer(hl)[p][ll] ) * D2ab[n].value;

                    dum_y += ( myPhiMO_y->pointer(hi)[p][ii] *
                               myPhiMO->pointer(hj)[p][jj] *
                               myPhiMO->pointer(hk)[p][kk] * 
                               myPhiMO->pointer(hl)[p][ll] +

                               myPhiMO->pointer(hi)[p][ii] *
                               myPhiMO_y->pointer(hj)[p][jj] *
                               myPhiMO->pointer(hk)[p][kk] * 
                               myPhiMO->pointer(hl)[p][ll] +

                               myPhiMO->pointer(hi)[p][ii] *
                               myPhiMO->pointer(hj)[p][jj] *
                               myPhiMO_y->pointer(hk)[p][kk] * 
                               myPhiMO->pointer(hl)[p][ll] +

                               myPhiMO->pointer(hi)[p][ii] *
                               myPhiMO->pointer(hj)[p][jj] *
                               myPhiMO->pointer(hk)[p][kk] * 
                               myPhiMO_y->pointer(hl)[p][ll] ) * D2ab[n].value;

                    dum_z += ( myPhiMO_z->pointer(hi)[p][ii] *
                               myPhiMO->pointer(hj)[p][jj] *
                               myPhiMO->pointer(hk)[p][kk] * 
                               myPhiMO->pointer(hl)[p][ll] +

                               myPhiMO->pointer(hi)[p][ii] *
                               myPhiMO_z->pointer(hj)[p][jj] *
                               myPhiMO->pointer(hk)[p][kk] * 
                               myPhiMO->pointer(hl)[p][ll] +

                               myPhiMO->pointer(hi)[p][ii] *
                               myPhiMO->pointer(hj)[p][jj] *
                               myPhiMO_z->pointer(hk)[p][kk] * 
                               myPhiMO->pointer(hl)[p][ll] +

                               myPhiMO->pointer(hi)[p][ii] *
                               myPhiMO->pointer(hj)[p][jj] *
                               myPhiMO->pointer(hk)[p][kk] * 
                               myPhiMO_z->pointer(hl)[p][ll] ) * D2ab[n].value;
                }

                pi_xp[phi_points_ + p] = dum_x;
                pi_yp[phi_points_ + p] = dum_y;
                pi_zp[phi_points_ + p] = dum_z;
            }

        }
        phi_points_ += npoints;
    }

    free(opdm_a_);
    free(opdm_b_);
}

void MCPDFTSolver::BuildPiFast(tpdm * D2ab, int nab) {

    pi_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double * pi_p = pi_->pointer();

    for (int p = 0; p < phi_points_; p++) {

        double dum = 0.0;

        // pi(r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
        for (int n = 0; n < nab; n++) {

            int i = D2ab[n].i;
            int j = D2ab[n].j;
            int k = D2ab[n].k;
            int l = D2ab[n].l;

            int hi = symmetry_[i];
            int hj = symmetry_[j];
            int hk = symmetry_[k];
            int hl = symmetry_[l];

            int ii = i - pitzer_offset_[hi];
            int jj = j - pitzer_offset_[hj];
            int kk = k - pitzer_offset_[hk];
            int ll = l - pitzer_offset_[hl];

            dum += super_phi_->pointer(hi)[p][ii] * 
                   super_phi_->pointer(hj)[p][jj] * 
                   super_phi_->pointer(hk)[p][kk] * 
                   super_phi_->pointer(hl)[p][ll] * D2ab[n].value;

        }

        pi_p[p] = dum;
    }

    if ( is_gga_ || is_meta_ ) {

        pi_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        pi_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        pi_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        
        double * pi_xp = pi_x_->pointer();
        double * pi_yp = pi_y_->pointer();
        double * pi_zp = pi_z_->pointer();
        
        double ** phi_x = super_phi_x_->pointer();
        double ** phi_y = super_phi_y_->pointer();
        double ** phi_z = super_phi_z_->pointer();

        for (int p = 0; p < phi_points_; p++) {

            double dum_x = 0.0;
            double dum_y = 0.0;
            double dum_z = 0.0;

            // pi(r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
            for (int n = 0; n < nab; n++) {

                int i = D2ab[n].i;
                int j = D2ab[n].j;
                int k = D2ab[n].k;
                int l = D2ab[n].l;
                
                int hi = symmetry_[i];
                int hj = symmetry_[j];
                int hk = symmetry_[k]; 
                int hl = symmetry_[l];
                
                int ii = i - pitzer_offset_[hi];
                int jj = j - pitzer_offset_[hj];
                int kk = k - pitzer_offset_[hk];
                int ll = l - pitzer_offset_[hl];

                dum_x += ( super_phi_x_->pointer(hi)[p][ii] *
                           super_phi_->pointer(hj)[p][jj] *
                           super_phi_->pointer(hk)[p][kk] * 
                           super_phi_->pointer(hl)[p][ll] +

                           super_phi_->pointer(hi)[p][ii] *
                           super_phi_x_->pointer(hj)[p][jj] *
                           super_phi_->pointer(hk)[p][kk] * 
                           super_phi_->pointer(hl)[p][ll] +

                           super_phi_->pointer(hi)[p][ii] *
                           super_phi_->pointer(hj)[p][jj] *
                           super_phi_x_->pointer(hk)[p][kk] * 
                           super_phi_->pointer(hl)[p][ll] +

                           super_phi_->pointer(hi)[p][ii] *
                           super_phi_->pointer(hj)[p][jj] *
                           super_phi_->pointer(hk)[p][kk] * 
                           super_phi_x_->pointer(hl)[p][ll] ) * D2ab[n].value;

                dum_y += ( super_phi_y_->pointer(hi)[p][ii] *
                           super_phi_->pointer(hj)[p][jj] *
                           super_phi_->pointer(hk)[p][kk] * 
                           super_phi_->pointer(hl)[p][ll] +

                           super_phi_->pointer(hi)[p][ii] *
                           super_phi_y_->pointer(hj)[p][jj] *
                           super_phi_->pointer(hk)[p][kk] * 
                           super_phi_->pointer(hl)[p][ll] +

                           super_phi_->pointer(hi)[p][ii] *
                           super_phi_->pointer(hj)[p][jj] *
                           super_phi_y_->pointer(hk)[p][kk] * 
                           super_phi_->pointer(hl)[p][ll] +

                           super_phi_->pointer(hi)[p][ii] *
                           super_phi_->pointer(hj)[p][jj] *
                           super_phi_->pointer(hk)[p][kk] * 
                           super_phi_y_->pointer(hl)[p][ll] ) * D2ab[n].value;

                dum_z += ( super_phi_z_->pointer(hi)[p][ii] *
                           super_phi_->pointer(hj)[p][jj] *
                           super_phi_->pointer(hk)[p][kk] * 
                           super_phi_->pointer(hl)[p][ll] +

                           super_phi_->pointer(hi)[p][ii] *
                           super_phi_z_->pointer(hj)[p][jj] *
                           super_phi_->pointer(hk)[p][kk] * 
                           super_phi_->pointer(hl)[p][ll] +

                           super_phi_->pointer(hi)[p][ii] *
                           super_phi_->pointer(hj)[p][jj] *
                           super_phi_z_->pointer(hk)[p][kk] * 
                           super_phi_->pointer(hl)[p][ll] +

                           super_phi_->pointer(hi)[p][ii] *
                           super_phi_->pointer(hj)[p][jj] *
                           super_phi_->pointer(hk)[p][kk] * 
                           super_phi_z_->pointer(hl)[p][ll] ) * D2ab[n].value;

            }

            pi_xp[p] = dum_x;
            pi_yp[p] = dum_y;
            pi_zp[p] = dum_z;

        }
    }
}

void MCPDFTSolver::BuildRho() {

    rho_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    rho_b_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    rho_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double ** phi   = super_phi_->pointer();

    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();
    double * rho_p = rho_->pointer();

    double temp_tot = 0.0;
    double temp_a = 0.0;
    double temp_b = 0.0;

    double ** dap = Da_->pointer();
    double ** dbp = Db_->pointer();

    for (int p = 0; p < phi_points_; p++) {
        double duma   = 0.0;
        double dumb   = 0.0;
        for (int sigma = 0; sigma < nmo_; sigma++) {
            for (int nu = 0; nu < nmo_; nu++) {
                duma += phi[p][sigma] * phi[p][nu] * dap[sigma][nu];
                dumb += phi[p][sigma] * phi[p][nu] * dbp[sigma][nu];
            }
        }
        rho_ap[p] = duma;
        rho_bp[p] = dumb;
        rho_p[p] = rho_ap[p] + rho_bp[p];    

        temp_tot += rho_p[p] * grid_w_->pointer()[p];
        temp_a += rho_ap[p] * grid_w_->pointer()[p];
        temp_b += rho_bp[p] * grid_w_->pointer()[p];

    }
    outfile->Printf("\n");
    outfile->Printf("      Integrated total density = %20.12lf\n",temp_tot);
    outfile->Printf("      Integrated alpha density = %20.12lf\n",temp_a);
    outfile->Printf("      Integrated beta density  = %20.12lf\n",temp_b);
    outfile->Printf("\n");

    if ( is_gga_ || is_meta_ ) {

        sigma_aa_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        sigma_bb_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        sigma_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        double ** phi_x = super_phi_x_->pointer();
        double ** phi_y = super_phi_y_->pointer();
        double ** phi_z = super_phi_z_->pointer();

        double * sigma_aap = sigma_aa_->pointer();
        double * sigma_bbp = sigma_bb_->pointer();
        double * sigma_abp = sigma_ab_->pointer();
        
        double * rho_a_xp = rho_a_x_->pointer();
        double * rho_b_xp = rho_b_x_->pointer();

        double * rho_a_yp = rho_a_y_->pointer();
        double * rho_b_yp = rho_b_y_->pointer();

        double * rho_a_zp = rho_a_z_->pointer();
        double * rho_b_zp = rho_b_z_->pointer();

        for (int p = 0; p < phi_points_; p++) {
            double duma_x = 0.0;
            double dumb_x = 0.0;
            double duma_y = 0.0;
            double dumb_y = 0.0;
            double duma_z = 0.0;
            double dumb_z = 0.0;

            for (int sigma = 0; sigma < nmo_; sigma++) {
                for (int nu = 0; nu < nmo_; nu++) {
                    duma_x += ( phi_x[p][sigma] * phi[p][nu] + phi[p][sigma] * phi_x[p][nu] ) * dap[sigma][nu];
                    dumb_x += ( phi_x[p][sigma] * phi[p][nu] + phi[p][sigma] * phi_x[p][nu] ) * dbp[sigma][nu];

                    duma_y += ( phi_y[p][sigma] * phi[p][nu] + phi[p][sigma] * phi_y[p][nu] ) * dap[sigma][nu];
                    dumb_y += ( phi_y[p][sigma] * phi[p][nu] + phi[p][sigma] * phi_y[p][nu] ) * dbp[sigma][nu];

                    duma_z += ( phi_z[p][sigma] * phi[p][nu] + phi[p][sigma] * phi_z[p][nu] ) * dap[sigma][nu];
                    dumb_z += ( phi_z[p][sigma] * phi[p][nu] + phi[p][sigma] * phi_z[p][nu] ) * dbp[sigma][nu];
                }
            }
            rho_a_xp[p] = duma_x;
            rho_b_xp[p] = dumb_x;

            rho_a_yp[p] = duma_y;
            rho_b_yp[p] = dumb_y;

            rho_a_zp[p] = duma_z;
            rho_b_zp[p] = dumb_z;

            sigma_aap[p] = ( rho_a_xp[p] * rho_a_xp[p] ) +  ( rho_a_yp[p] * rho_a_yp[p] ) + ( rho_a_zp[p] * rho_a_zp[p] );
            sigma_bbp[p] = ( rho_b_xp[p] * rho_b_xp[p] ) +  ( rho_b_yp[p] * rho_b_yp[p] ) + ( rho_b_zp[p] * rho_b_zp[p] );
            sigma_abp[p] = ( rho_a_xp[p] * rho_b_xp[p] ) +  ( rho_a_yp[p] * rho_b_yp[p] ) + ( rho_a_zp[p] * rho_b_zp[p] );
        }
    }
}

void MCPDFTSolver::BuildRhoFast(int na, int nb) {

    rho_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    rho_b_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    rho_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();
    double * rho_p = rho_->pointer();

    double temp_tot = 0.0;
    double temp_a = 0.0;
    double temp_b = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        // rho_a(r)
        double duma = 0.0;
        for (int n = 0; n < na; n++) {

            int i = opdm_a_[n].i;
            int j = opdm_a_[n].j;

            int hi = symmetry_[i];
            int hj = symmetry_[j];

            int ii = i - pitzer_offset_[hi];
            int jj = j - pitzer_offset_[hj];

            duma += super_phi_->pointer(hi)[p][ii] *
                    super_phi_->pointer(hj)[p][jj] * opdm_a_[n].value;

        }
        rho_ap[p] = duma;

        // rho_b(r)
        double dumb = 0.0;
        for (int n = 0; n < nb; n++) {

            int i = opdm_b_[n].i;
            int j = opdm_b_[n].j;

            int hi = symmetry_[i];
            int hj = symmetry_[j];

            int ii = i - pitzer_offset_[hi];
            int jj = j - pitzer_offset_[hj];

            dumb += super_phi_->pointer(hi)[p][ii] *
                    super_phi_->pointer(hj)[p][jj] * opdm_b_[n].value;

        }
        rho_bp[p] = dumb;

        rho_p[p] = rho_ap[p] + rho_bp[p];    

        temp_tot += rho_p[p]  * grid_w_->pointer()[p];
        temp_a   += rho_ap[p] * grid_w_->pointer()[p];
        temp_b   += rho_bp[p] * grid_w_->pointer()[p];

    }
    outfile->Printf("\n");
    outfile->Printf("      Integrated total density = %20.12lf\n",temp_tot);
    outfile->Printf("      Integrated alpha density = %20.12lf\n",temp_a);
    outfile->Printf("      Integrated beta density  = %20.12lf\n",temp_b);
    outfile->Printf("\n");

    if ( is_gga_ || is_meta_ ) {

        sigma_aa_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        sigma_bb_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        sigma_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        tau_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        tau_b_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        tw_      = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        double ** phi_x = super_phi_x_->pointer();
        double ** phi_y = super_phi_y_->pointer();
        double ** phi_z = super_phi_z_->pointer();

        double * sigma_aap = sigma_aa_->pointer();
        double * sigma_bbp = sigma_bb_->pointer();
        double * sigma_abp = sigma_ab_->pointer();

        double * rho_a_xp = rho_a_x_->pointer();
        double * rho_b_xp = rho_b_x_->pointer();

        double * rho_a_yp = rho_a_y_->pointer();
        double * rho_b_yp = rho_b_y_->pointer();

        double * rho_a_zp = rho_a_z_->pointer();
        double * rho_b_zp = rho_b_z_->pointer();

        double * tau_ap = tau_a_->pointer();
        double * tau_bp = tau_b_->pointer();

        double * tw_p = tw_->pointer();

        for (int p = 0; p < phi_points_; p++) {

            // rho'_a(r)

            double duma_x = 0.0;
            double duma_y = 0.0;
            double duma_z = 0.0;
            double dumta = 0.0;

            for (int n = 0; n < na; n++) {

                int i = opdm_a_[n].i;
                int j = opdm_a_[n].j;
                
                int hi = symmetry_[i];
                int hj = symmetry_[j];
                
                int ii = i - pitzer_offset_[hi];
                int jj = j - pitzer_offset_[hj];
                
                duma_x += ( super_phi_x_->pointer(hi)[p][ii] * super_phi_->pointer(hj)[p][jj] 
                        +   super_phi_->pointer(hi)[p][ii] * super_phi_x_->pointer(hj)[p][jj] ) * opdm_a_[n].value;

                duma_y += ( super_phi_y_->pointer(hi)[p][ii] * super_phi_->pointer(hj)[p][jj] 
                        +   super_phi_->pointer(hi)[p][ii] * super_phi_y_->pointer(hj)[p][jj] ) * opdm_a_[n].value;

                duma_z += ( super_phi_z_->pointer(hi)[p][ii] * super_phi_->pointer(hj)[p][jj] 
                        +   super_phi_->pointer(hi)[p][ii] * super_phi_z_->pointer(hj)[p][jj] ) * opdm_a_[n].value;

                dumta  += ( super_phi_x_->pointer(hi)[p][ii] * super_phi_x_->pointer(hj)[p][jj] 
		        +   super_phi_y_->pointer(hi)[p][ii] * super_phi_y_->pointer(hj)[p][jj]
                        +   super_phi_z_->pointer(hi)[p][ii] * super_phi_z_->pointer(hj)[p][jj] ) * opdm_a_[n].value;
            }

            // rho'_b(r)

            double dumb_x = 0.0;
            double dumb_y = 0.0;
            double dumb_z = 0.0;
            double dumtb = 0.0;

            for (int n = 0; n < nb; n++) {
                
                int i = opdm_b_[n].i;
                int j = opdm_b_[n].j;
                
                int hi = symmetry_[i];
                int hj = symmetry_[j];
                
                int ii = i - pitzer_offset_[hi];
                int jj = j - pitzer_offset_[hj];
                
                dumb_x += ( super_phi_x_->pointer(hi)[p][ii] * super_phi_->pointer(hj)[p][jj] 
                        +   super_phi_->pointer(hi)[p][ii] * super_phi_x_->pointer(hj)[p][jj] ) * opdm_b_[n].value;
                
                dumb_y += ( super_phi_y_->pointer(hi)[p][ii] * super_phi_->pointer(hj)[p][jj] 
                        +   super_phi_->pointer(hi)[p][ii] * super_phi_y_->pointer(hj)[p][jj] ) * opdm_b_[n].value;
                
                dumb_z += ( super_phi_z_->pointer(hi)[p][ii] * super_phi_->pointer(hj)[p][jj] 
                        +   super_phi_->pointer(hi)[p][ii] * super_phi_z_->pointer(hj)[p][jj] ) * opdm_b_[n].value;

                dumtb  += ( super_phi_x_->pointer(hi)[p][ii] * super_phi_x_->pointer(hj)[p][jj] 
		        +   super_phi_y_->pointer(hi)[p][ii] * super_phi_y_->pointer(hj)[p][jj]
                        +   super_phi_z_->pointer(hi)[p][ii] * super_phi_z_->pointer(hj)[p][jj] ) * opdm_b_[n].value;
            }

            rho_a_xp[p] = duma_x;
            rho_b_xp[p] = dumb_x;

            rho_a_yp[p] = duma_y;
            rho_b_yp[p] = dumb_y;

            rho_a_zp[p] = duma_z;
            rho_b_zp[p] = dumb_z;

            sigma_aap[p] = ( rho_a_xp[p] * rho_a_xp[p] ) +  ( rho_a_yp[p] * rho_a_yp[p] ) + ( rho_a_zp[p] * rho_a_zp[p] );
            sigma_bbp[p] = ( rho_b_xp[p] * rho_b_xp[p] ) +  ( rho_b_yp[p] * rho_b_yp[p] ) + ( rho_b_zp[p] * rho_b_zp[p] );
            sigma_abp[p] = ( rho_a_xp[p] * rho_b_xp[p] ) +  ( rho_a_yp[p] * rho_b_yp[p] ) + ( rho_a_zp[p] * rho_b_zp[p] );

            tau_ap[p] = dumta;
            tau_bp[p] = dumtb;

	    tw_p[p] = ( sigma_aap[p] + 2.0 * sigma_abp[p] + sigma_bbp[p] ) / (8.0 * rho_p[p]);
        }
    }
}

std::vector< std::shared_ptr<Matrix> > MCPDFTSolver::BuildJK() {

    std::shared_ptr<PSIO> psio (new PSIO());

    if ( !psio->exists(PSIF_V2RDM_D1A) ) throw PsiException("No D1a on disk",__FILE__,__LINE__);
    if ( !psio->exists(PSIF_V2RDM_D1B) ) throw PsiException("No D1b on disk",__FILE__,__LINE__);

    // D1a

    psio->open(PSIF_V2RDM_D1A,PSIO_OPEN_OLD);

    long int na;
    psio->read_entry(PSIF_V2RDM_D1A,"length",(char*)&na,sizeof(long int));

    opdm_a_ = (opdm *)malloc(na * sizeof(opdm));
    psio->read_entry(PSIF_V2RDM_D1A,"D1a",(char*)opdm_a_,na * sizeof(opdm));
    psio->close(PSIF_V2RDM_D1A,1);

    for (int n = 0; n < na; n++) {

        int i = opdm_a_[n].i;
        int j = opdm_a_[n].j;

        int hi = symmetry_[i];
        int hj = symmetry_[j];

        if ( hi != hj ) {
            throw PsiException("error: something is wrong with the symmetry of the alpha OPDM",__FILE__,__LINE__);
        }

        int ii = i - pitzer_offset_[hi];
        int jj = j - pitzer_offset_[hi];

        Da_->pointer(hi)[ii][jj] = opdm_a_[n].value;

    }

    // D1b

    psio->open(PSIF_V2RDM_D1B,PSIO_OPEN_OLD);

    long int nb;
    psio->read_entry(PSIF_V2RDM_D1B,"length",(char*)&nb,sizeof(long int));

    opdm_b_ = (opdm *)malloc(nb * sizeof(opdm));
    psio->read_entry(PSIF_V2RDM_D1B,"D1b",(char*)opdm_b_,nb * sizeof(opdm));
    psio->close(PSIF_V2RDM_D1B,1);

   for (int n = 0; n < nb; n++) {

        int i = opdm_b_[n].i;
        int j = opdm_b_[n].j;

        int hi = symmetry_[i];
        int hj = symmetry_[j];

        if ( hi != hj ) {
            throw PsiException("error: something is wrong with the symmetry of the beta OPDM",__FILE__,__LINE__);
        }

        int ii = i - pitzer_offset_[hi];
        int jj = j - pitzer_offset_[hi];

        Db_->pointer(hi)[ii][jj] = opdm_b_[n].value;

    }

    free(opdm_a_);
    free(opdm_b_);

    // get primary basis:
    std::shared_ptr<BasisSet> primary = reference_wavefunction_->get_basisset("ORBITAL");

    std::vector< std::shared_ptr<Matrix> > JK;

    // JK object (note this is hard-coded to use density fitting ...)
    if (options_.get_str("MCPDFT_TYPE") == "DF") {

         // get auxiliary basis:
         std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");
        
        std::shared_ptr<DiskDFJK> jk = (std::shared_ptr<DiskDFJK>)(new DiskDFJK(primary,auxiliary,options_));
        
        // memory for jk 
        jk->set_memory(0.90 * Process::environment.get_memory());

        // integral cutoff
        jk->set_cutoff(options_.get_double("INTS_TOLERANCE"));

        jk->set_do_J(true);
        jk->set_do_K(false);
        jk->set_do_wK(false);
        jk->set_omega(false);

        jk->initialize();

        std::vector<SharedMatrix>& C_left  = jk->C_left();
        std::vector<SharedMatrix>& C_right = jk->C_right();

        // use alpha and beta C matrices
        std::shared_ptr<Matrix> myCa (new Matrix(Ca_) );
        std::shared_ptr<Matrix> myCb (new Matrix(Cb_) );

        C_left.clear();
        C_right.clear();

        // alpha first 
        myCa->zero();
        myCa->gemm('t','n',1.0,Da_,Ca_,0.0);
        myCa->transpose_this();
        C_left.push_back(myCa);
        C_right.push_back(Ca_);

        // beta second 
        myCb->zero();
        myCb->gemm('t','n',1.0,Db_,Cb_,0.0);
        myCb->transpose_this();
        C_left.push_back(myCb);
        C_right.push_back(Cb_);

        // Let jk compute for the given C_left/C_right

        jk->compute();

        std::shared_ptr<Matrix> Ja = jk->J()[0];
        std::shared_ptr<Matrix> Jb = jk->J()[1];

        Ja->transform(Ca_);
        Jb->transform(Cb_);

        JK.push_back(Ja);
        JK.push_back(Jb);

        jk.reset();

        return JK;

    }else if (options_.get_str("MCPDFT_TYPE") == "PK") { 

        std::shared_ptr<PKJK> jk = (std::shared_ptr<PKJK>)(new PKJK(primary,options_));

        // memory for jk 
        jk->set_memory(0.90 * Process::environment.get_memory());

        // integral cutoff
        jk->set_cutoff(options_.get_double("INTS_TOLERANCE"));

        jk->set_do_J(true);
        jk->set_do_K(false);
        jk->set_do_wK(false);
        jk->set_omega(false);

        jk->initialize();

        std::vector<SharedMatrix>& C_left  = jk->C_left();
        std::vector<SharedMatrix>& C_right = jk->C_right();

        // use alpha and beta C matrices

        std::shared_ptr<Matrix> myCa (new Matrix(Ca_) );
        std::shared_ptr<Matrix> myCb (new Matrix(Cb_) );

        C_left.clear();
        C_right.clear();

        // alpha first 
        myCa->zero();
        myCa->gemm('t','n',1.0,Da_,Ca_,0.0);
        myCa->transpose_this();
        C_left.push_back(myCa);
        C_right.push_back(Ca_);

        // beta second 
        myCb->zero();
        myCb->gemm('t','n',1.0,Db_,Cb_,0.0);
        myCb->transpose_this();
        C_left.push_back(myCb);
        C_right.push_back(Cb_);

        // Let jk compute for the given C_left/C_right

        jk->compute();

        std::shared_ptr<Matrix> Ja = jk->J()[0];
        std::shared_ptr<Matrix> Jb = jk->J()[1];

        Ja->transform(Ca_);
        Jb->transform(Cb_);

        JK.push_back(Ja);
        JK.push_back(Jb);

        return JK;
        
    }else {
        throw PsiException("invalid MCSCF_TYPE",__FILE__,__LINE__);
    }
}

void MCPDFTSolver::Build_R(){

    double tol = 1.0e-20;

    R_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    
    double * R_p = R_->pointer();
    double * rho_p = rho_->pointer();
    double * pi_p = pi_->pointer();
 
    for (int p = 0; p < phi_points_; p++) {
        
        R_p[p] = 4 * pi_p[p] / (rho_p[p] * rho_p[p]); 
    }
} 

void MCPDFTSolver::Translate(){

    double tol = 1.0e-20;

    tr_rho_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    tr_rho_b_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double * tr_rho_ap = tr_rho_a_->pointer();
    double * tr_rho_bp = tr_rho_b_->pointer();

    double * rho_p = rho_->pointer();
    double * pi_p = pi_->pointer();
    double * R_p = R_->pointer();

    double temp_tot = 0.0;
    double temp_a = 0.0;
    double temp_b = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double rho = rho_p[p];
        double pi = pi_p[p];
        
        double zeta = 0.0;
        double R = 0.0;

	if ( !(rho < tol) && !(pi < 0.0) ) {

            R = R_p[p];

            if ( (1.0 - R) > tol ) {
                zeta = sqrt(1.0 - R);
            }else{
                zeta = 0.0;
            }

            tr_rho_ap[p] = (1.0 + zeta) * (rho/2.0);
            tr_rho_bp[p] = (1.0 - zeta) * (rho/2.0);

        }else {
            tr_rho_ap[p] = 0.0;
            tr_rho_bp[p] = 0.0;
        }

        temp_tot += (tr_rho_ap[p] + tr_rho_bp[p]) * grid_w_->pointer()[p];
        temp_a += tr_rho_ap[p] * grid_w_->pointer()[p];
        temp_b += tr_rho_bp[p] * grid_w_->pointer()[p];
    }

    outfile->Printf("\n");
    outfile->Printf("      Integrated translated total density = %20.12lf\n",temp_tot);
    outfile->Printf("      Integrated translated alpha density = %20.12lf\n",temp_a);
    outfile->Printf("      Integrated translated beta density  = %20.12lf\n",temp_b);
    outfile->Printf("\n");

    if ( is_gga_ || is_meta_ ) {

        tr_rho_a_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        tr_rho_b_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        tr_rho_a_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        tr_rho_b_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        tr_rho_a_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        tr_rho_b_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        tr_sigma_aa_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        tr_sigma_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        tr_sigma_bb_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        double * rho_a_xp = rho_a_x_->pointer();
        double * rho_b_xp = rho_b_x_->pointer();

        double * rho_a_yp = rho_a_y_->pointer();
        double * rho_b_yp = rho_b_y_->pointer();

        double * rho_a_zp = rho_a_z_->pointer();
        double * rho_b_zp = rho_b_z_->pointer();

        double * tr_rho_a_xp = tr_rho_a_x_->pointer();
        double * tr_rho_b_xp = tr_rho_b_x_->pointer();

        double * tr_rho_a_yp = tr_rho_a_y_->pointer();
        double * tr_rho_b_yp = tr_rho_b_y_->pointer();

        double * tr_rho_a_zp = tr_rho_a_z_->pointer();
        double * tr_rho_b_zp = tr_rho_b_z_->pointer();

        double * tr_sigma_aap = tr_sigma_aa_->pointer();
        double * tr_sigma_abp = tr_sigma_ab_->pointer();
        double * tr_sigma_bbp = tr_sigma_bb_->pointer();

        for (int p = 0; p < phi_points_; p++) {

            double rho = rho_p[p];
            double pi = pi_p[p];

            double rho_x = rho_a_xp[p] + rho_b_xp[p];
            double rho_y = rho_a_yp[p] + rho_b_yp[p];
            double rho_z = rho_a_zp[p] + rho_b_zp[p];
            
            double zeta = 0.0;
            double R = 0.0;

            if ( !(rho < tol) && !(pi < 0.0) ) {

                R = R_p[p];
                // R = tanh(R);
                 
                if ( (1.0 - R) > tol )  {

                    zeta = sqrt(1.0 - R);

                }else{

                    zeta = 0.0;
                }

                tr_rho_a_xp[p] = (1.0 + zeta) * (rho_x/2.0);
                tr_rho_b_xp[p] = (1.0 - zeta) * (rho_x/2.0);

                tr_rho_a_yp[p] = (1.0 + zeta) * (rho_y/2.0);
                tr_rho_b_yp[p] = (1.0 - zeta) * (rho_y/2.0);

                tr_rho_a_zp[p] = (1.0 + zeta) * (rho_z/2.0);
                tr_rho_b_zp[p] = (1.0 - zeta) * (rho_z/2.0);

            }else {

                tr_rho_a_xp[p] = 0.0;
                tr_rho_b_xp[p] = 0.0;

                tr_rho_a_yp[p] = 0.0;
                tr_rho_b_yp[p] = 0.0;

                tr_rho_a_zp[p] = 0.0;
                tr_rho_b_zp[p] = 0.0;
            }

            tr_sigma_aap[p] = (tr_rho_a_xp[p] * tr_rho_a_xp[p]) + (tr_rho_a_yp[p] * tr_rho_a_yp[p]) + (tr_rho_a_zp[p] * tr_rho_a_zp[p]);
            tr_sigma_abp[p] = (tr_rho_a_xp[p] * tr_rho_b_xp[p]) + (tr_rho_a_yp[p] * tr_rho_b_yp[p]) + (tr_rho_a_zp[p] * tr_rho_b_zp[p]);
            tr_sigma_bbp[p] = (tr_rho_b_xp[p] * tr_rho_b_xp[p]) + (tr_rho_b_yp[p] * tr_rho_b_yp[p]) + (tr_rho_b_zp[p] * tr_rho_b_zp[p]);
        }
    }
}

void MCPDFTSolver::Fully_Translate(){
    
    double tol = 1.0e-20;

    double const R0 = 0.9;
    double const R1 = 1.15;
    double const A = -475.60656009;
    double const B = -379.47331922;
    double const C = -85.38149682;

    tr_rho_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    tr_rho_b_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    
    double * tr_rho_ap = tr_rho_a_->pointer();
    double * tr_rho_bp = tr_rho_b_->pointer();
    
    double * rho_p = rho_->pointer();
    double * pi_p = pi_->pointer();
    double * R_p = R_->pointer();

    double temp_tot = 0.0;
    double temp_a = 0.0;
    double temp_b = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double rho = rho_p[p];
        double pi = pi_p[p];
        double DelR = R_p[p] - R1;

        double zeta = 0.0;
        double R = 0.0;

	if ( !(rho < tol) && !(pi < 0.0) ) {

            R = R_p[p];

            if ( ((1.0 - R) > tol) && ( R < R0 ) ) {

                zeta = sqrt(1.0 - R);

            }else if( !(R < R0) && !(R > R1) ) {

                zeta = A * pow(DelR, 5.0) + B * pow(DelR, 4.0) + C * pow(DelR, 3.0);

            }else if( R > R1 ) {

                zeta = 0.0;

            }

            tr_rho_ap[p] = (1.0 + zeta) * (rho/2.0);
            tr_rho_bp[p] = (1.0 - zeta) * (rho/2.0);

        }else{

            tr_rho_ap[p] = 0.0;
            tr_rho_bp[p] = 0.0;
        }

        temp_a += tr_rho_ap[p] * grid_w_->pointer()[p];
        temp_b += tr_rho_bp[p] * grid_w_->pointer()[p];
        temp_tot += ( tr_rho_bp[p] + tr_rho_ap[p] ) * grid_w_->pointer()[p];

    }

    outfile->Printf("\n");
    outfile->Printf("      Integrated fully translated total density = %20.12lf\n",temp_tot);
    outfile->Printf("      Integrated fully translated alpha density = %20.12lf\n",temp_a);
    outfile->Printf("      Integrated fully translated beta density  = %20.12lf\n",temp_b);
    outfile->Printf("\n");

    if ( is_gga_ || is_meta_ ) {

        tr_rho_a_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        tr_rho_b_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        tr_rho_a_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        tr_rho_b_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        tr_rho_a_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        tr_rho_b_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        
        tr_sigma_aa_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        tr_sigma_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        tr_sigma_bb_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        double * tr_rho_a_xp = tr_rho_a_x_->pointer();
        double * tr_rho_b_xp = tr_rho_b_x_->pointer();

        double * tr_rho_a_yp = tr_rho_a_y_->pointer();
        double * tr_rho_b_yp = tr_rho_b_y_->pointer();
 
        double * tr_rho_a_zp = tr_rho_a_z_->pointer();
        double * tr_rho_b_zp = tr_rho_b_z_->pointer();
        
        double * tr_sigma_aap = tr_sigma_aa_->pointer();
        double * tr_sigma_abp = tr_sigma_ab_->pointer();
        double * tr_sigma_bbp = tr_sigma_bb_->pointer();

        double * rho_a_xp = rho_a_x_->pointer();
        double * rho_b_xp = rho_b_x_->pointer();

        double * rho_a_yp = rho_a_y_->pointer();
        double * rho_b_yp = rho_b_y_->pointer();
 
        double * rho_a_zp = rho_a_z_->pointer();
        double * rho_b_zp = rho_b_z_->pointer();

        double * pi_xp = pi_x_->pointer();
        double * pi_yp = pi_y_->pointer();
        double * pi_zp = pi_z_->pointer();

        for (int p = 0; p < phi_points_; p++) {

            double rho_x = rho_a_xp[p] + rho_b_xp[p];
            double rho_y = rho_a_yp[p] + rho_b_yp[p];
            double rho_z = rho_a_zp[p] + rho_b_zp[p];

            double rho = rho_p[p];
            double pi = pi_p[p];
            double DelR = R_p[p] - R1;

            double zeta = 0.0;
            double R = 0.0;
 
            if ( !(rho < tol) && !(pi < 0.0) ) {
            
                R = R_p[p];
            
                if ( ((1.0 - R) > tol) && ( R < R0 ) ) {
            
                    zeta = sqrt(1.0 - R);
            
                    tr_rho_a_xp[p] = (1.0 + zeta) * (rho_x/2.0) + (R * rho_x) / (2.0*zeta) - pi_xp[p] / (rho*zeta);
                    tr_rho_b_xp[p] = (1.0 - zeta) * (rho_x/2.0) - (R * rho_x) / (2.0*zeta) + pi_xp[p] / (rho*zeta);
                    
                    tr_rho_a_yp[p] = (1.0 + zeta) * (rho_y/2.0) + (R * rho_y) / (2.0*zeta) - pi_yp[p] / (rho*zeta);
                    tr_rho_b_yp[p] = (1.0 - zeta) * (rho_y/2.0) - (R * rho_y) / (2.0*zeta) + pi_yp[p] / (rho*zeta);

                    tr_rho_a_zp[p] = (1.0 + zeta) * (rho_z/2.0) + (R * rho_z) / (2.0*zeta) - pi_zp[p] / (rho*zeta);
                    tr_rho_b_zp[p] = (1.0 - zeta) * (rho_z/2.0) - (R * rho_z) / (2.0*zeta) + pi_zp[p] / (rho*zeta);

                }else if( !(R < R0) && !(R > R1) ) {
            
                    zeta = A * pow(DelR, 5.0) + B * pow(DelR, 4.0) + C * pow(DelR, 3.0);

                    tr_rho_a_xp[p] = (1.0 + zeta) * (rho_x/2.0) 
                                    + (A * pow(DelR, 4.0)) * ( (10.0 * pi_xp[p] / rho) - (5.0 * R * rho_x) )
                                    + (B * pow(DelR, 3.0)) * ( (8.0  * pi_xp[p] / rho) - (4.0 * R * rho_x) )
                                    + (C * pow(DelR, 2.0)) * ( (6.0  * pi_xp[p] / rho) - (3.0 * R * rho_x) );

                    tr_rho_b_xp[p] = (1.0 - zeta) * (rho_x/2.0) 
                                    + (A * pow(DelR, 4.0)) * (-(10.0 * pi_xp[p] / rho) + (5.0 * R * rho_x) )
                                    + (B * pow(DelR, 3.0)) * (-(8.0  * pi_xp[p] / rho) + (4.0 * R * rho_x) )
                                    + (C * pow(DelR, 2.0)) * (-(6.0  * pi_xp[p] / rho) + (3.0 * R * rho_x) );
            
                    tr_rho_a_yp[p] = (1.0 + zeta) * (rho_y/2.0) 
                                    + (A * pow(DelR, 4.0)) * ( (10.0 * pi_yp[p] / rho) - (5.0 * R * rho_y) )
                                    + (B * pow(DelR, 3.0)) * ( (8.0  * pi_yp[p] / rho) - (4.0 * R * rho_y) )
                                    + (C * pow(DelR, 2.0)) * ( (6.0  * pi_yp[p] / rho) - (3.0 * R * rho_y) );

                    tr_rho_b_yp[p] = (1.0 - zeta) * (rho_y/2.0) 
                                    + (A * pow(DelR, 4.0)) * (-(10.0 * pi_yp[p] / rho) + (5.0 * R * rho_y) )
                                    + (B * pow(DelR, 3.0)) * (-(8.0  * pi_yp[p] / rho) + (4.0 * R * rho_y) )
                                    + (C * pow(DelR, 2.0)) * (-(6.0  * pi_yp[p] / rho) + (3.0 * R * rho_y) );

                    tr_rho_a_zp[p] = (1.0 + zeta) * (rho_z/2.0) 
                                    + (A * pow(DelR, 4.0)) * ( (10.0 * pi_zp[p] / rho) - (5.0 * R * rho_z) )
                                    + (B * pow(DelR, 3.0)) * ( (8.0  * pi_zp[p] / rho) - (4.0 * R * rho_z) )
                                    + (C * pow(DelR, 2.0)) * ( (6.0  * pi_zp[p] / rho) - (3.0 * R * rho_z) );

                    tr_rho_b_zp[p] = (1.0 - zeta) * (rho_z/2.0) 
                                    + (A * pow(DelR, 4.0)) * (-(10.0 * pi_zp[p] / rho) + (5.0 * R * rho_z) )
                                    + (B * pow(DelR, 3.0)) * (-(8.0  * pi_zp[p] / rho) + (4.0 * R * rho_z) )
                                    + (C * pow(DelR, 2.0)) * (-(6.0  * pi_zp[p] / rho) + (3.0 * R * rho_z) );
                }else if( R > R1 ) {
            
                    zeta = 0.0;

                    tr_rho_a_xp[p] = (1.0 + zeta) * (rho_x/2.0);
                    tr_rho_b_xp[p] = (1.0 - zeta) * (rho_x/2.0);

                    tr_rho_a_yp[p] = (1.0 + zeta) * (rho_y/2.0);
                    tr_rho_b_yp[p] = (1.0 - zeta) * (rho_y/2.0);

                    tr_rho_a_zp[p] = (1.0 + zeta) * (rho_z/2.0);
                    tr_rho_b_zp[p] = (1.0 - zeta) * (rho_z/2.0);
                }
            
            }else{
            
                tr_rho_a_xp[p] = 0.0;
                tr_rho_b_xp[p] = 0.0;
                
                tr_rho_a_yp[p] = 0.0;
                tr_rho_b_yp[p] = 0.0;
                
                tr_rho_a_zp[p] = 0.0;
                tr_rho_b_zp[p] = 0.0;
            }
            // outfile->Printf("\n     Fully translated density gradient (x-component) =      %12.15lf\n",tr_rho_a_xp[p]);

            tr_sigma_aap[p] = (tr_rho_a_xp[p] * tr_rho_a_xp[p]) + (tr_rho_a_yp[p] * tr_rho_a_yp[p]) + (tr_rho_a_zp[p] * tr_rho_a_zp[p]);  
            tr_sigma_abp[p] = (tr_rho_a_xp[p] * tr_rho_b_xp[p]) + (tr_rho_a_yp[p] * tr_rho_b_yp[p]) + (tr_rho_a_zp[p] * tr_rho_b_zp[p]);  
            tr_sigma_bbp[p] = (tr_rho_b_xp[p] * tr_rho_b_xp[p]) + (tr_rho_b_yp[p] * tr_rho_b_yp[p]) + (tr_rho_b_zp[p] * tr_rho_b_zp[p]);  

        }
    }
}

} // end of namespaces

