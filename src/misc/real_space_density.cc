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

// psi4 stuff
#include "psi4/libpsi4util/process.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libscf_solver/hf.h"

// real_space_density 
#include "real_space_density.h"

// blas
#include <misc/blas.h>

// openmp
#include <misc/omp.h>

// tpdm and opdm structs live here
#include <v2rdm_casscf/v2rdm_solver.h>

// cs-ksdft
#include <cs_ksdft/cs_solver.h>

using namespace psi;
using namespace fnocc;

namespace hilbert{ 

// the RealSpaceDensity class derives from the Wavefunction class and inherits its members
RealSpaceDensity::RealSpaceDensity(std::shared_ptr<Wavefunction> reference_wavefunction,Options & options_):
    Wavefunction(options_){

    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

RealSpaceDensity::~RealSpaceDensity() {
}

// initialize members of the RealSpaceDensity class
void RealSpaceDensity::common_init() {

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

    // SO-basis Fock matrices
    Fa_ = std::shared_ptr<Matrix>(reference_wavefunction_->Fa());
    Fb_ = std::shared_ptr<Matrix>(reference_wavefunction_->Fb());

    // SO-basis density matrices
    Da_ = std::shared_ptr<Matrix>(reference_wavefunction_->Da());
    Db_ = std::shared_ptr<Matrix>(reference_wavefunction_->Db());

    // orbital energies
    epsilon_a_ = std::make_shared<Vector>(nmopi_);
    epsilon_a_->copy(*reference_wavefunction_->epsilon_a().get());
    epsilon_b_ = std::make_shared<Vector>(nmopi_);
    epsilon_b_->copy(*reference_wavefunction_->epsilon_b().get());

    // set the wavefunction name
    name_ = "WHOKNOWS";

    // restricted orbitals, unrestricted rdms
// TODO get these from reference wave function
    same_a_b_orbs_ = false;
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

    // estimate memory requirements

    //outfile->Printf("\n"); 
    //outfile->Printf("    ==> Memory requirements <==\n");
    //outfile->Printf("\n");
    
    // memory is from process::environment
    memory_ = Process::environment.get_memory();

    // assume 500 MB is already allocated by psi4
    long int n_mem         = 500 * 1024 * 1024 / 8;

    // default basic vectors and matrices initialized using common_init()
    n_mem +=     nso_ * nso_;                             // S_ matrix
    n_mem += 2 * nso_ * nso_;                             // Ca_ and Cb_ matrices
    n_mem += 2 * nso_ * nso_;                             // Fa_ and Fb_ matrices
    n_mem += 2 * nmo_;                                    // epsilon_a_ and epsilon_b_ vectors

    // phiAO, phiAO_x, phiAO_y, phiAO_z
    n_mem += phi_points_ * nso_;
    n_mem += 3 * phi_points_ * nso_;

    // required memory for grids x, y, and z and weights w vectors
    n_mem += 4 * phi_points_;

    // memory for AO->MO transformation
    n_mem += nirrep_ * phi_points_;    // phi_points_list vector
    n_mem += phi_points_ * nso_;       // super_phi_ matrix
    n_mem += nso_ * nso_;              // AO2SO matrix in TransformPhiMatrixAOMO()
    n_mem += phi_points_ * nso_;       // temporary matrix called three times in TransformPhiMatrixAOMO()
    n_mem += 3 * phi_points_ * nso_;   // for super_phi_x_, _y_ and _z_ gradient matrices

    // memory for rho and rho'
    n_mem += 3 * phi_points_;          // rho_a_, rho_b_ and rho_ vectors
    n_mem += 3 * phi_points_;          // rho_a_x_, rho_a_y_ and rho_a_z_ gradient vectors
    n_mem += 3 * phi_points_;          // rho_b_x_, rho_b_y_ and rho_b_z_ gradient vectors

    // memory needed for storing pi vector
    n_mem += phi_points_;              // pi_ vector
    n_mem += 3 * phi_points_;          // pi_x_, pi_y_ and pi_z_ gradient vectors
 
    if ( n_mem * 8.0 > (double)memory_ ) {
        throw PsiException("not enough memory. see plugin mcpdft for a low-memory algorithm",__FILE__,__LINE__);
    }

    grid_x_      = std::shared_ptr<Vector>(new Vector("GRID X",phi_points_));
    grid_y_      = std::shared_ptr<Vector>(new Vector("GRID Y",phi_points_));
    grid_z_      = std::shared_ptr<Vector>(new Vector("GRID Z",phi_points_));
    grid_w_      = std::shared_ptr<Vector>(new Vector("GRID W",phi_points_));

    outfile->Printf("\n"); 
    outfile->Printf("    ==> Build Phi and Phi' matrices ..."); fflush(stdout);

    std::shared_ptr<Matrix> super_phi_ao = std::shared_ptr<Matrix>(new Matrix("SUPER PHI (AO)",phi_points_,nso_));

    std::shared_ptr<Matrix> super_phi_x_ao (new Matrix("SUPER PHI X (AO)",phi_points_,nso_));
    std::shared_ptr<Matrix> super_phi_y_ao (new Matrix("SUPER PHI Y (AO)",phi_points_,nso_));
    std::shared_ptr<Matrix> super_phi_z_ao (new Matrix("SUPER PHI Z (AO)",phi_points_,nso_));

    // build phi matrix and derivative phi matrices (ao basis)
    BuildPhiMatrixAO("PHI",super_phi_ao);
    BuildPhiMatrixAO("PHI_X",super_phi_x_ao);
    BuildPhiMatrixAO("PHI_Y",super_phi_y_ao);
    BuildPhiMatrixAO("PHI_Z",super_phi_z_ao);

    // transform the orbital label in phi matrix and derivative phi matrix to the MO basis

    Dimension phi_points_list = Dimension(nirrep_,"phi_points");
    for (int h = 0; h < nirrep_; h++) {
        phi_points_list[h] = phi_points_;
    }

    super_phi_ = std::shared_ptr<Matrix>(new Matrix("SUPER PHI",phi_points_list,nsopi_));

    TransformPhiMatrixAOMO(super_phi_ao,super_phi_);

    super_phi_x_ = std::shared_ptr<Matrix>(new Matrix("SUPER PHI X",phi_points_list,nsopi_));
    super_phi_y_ = std::shared_ptr<Matrix>(new Matrix("SUPER PHI Y",phi_points_list,nsopi_));
    super_phi_z_ = std::shared_ptr<Matrix>(new Matrix("SUPER PHI Z",phi_points_list,nsopi_));

    TransformPhiMatrixAOMO(super_phi_x_ao,super_phi_x_);
    TransformPhiMatrixAOMO(super_phi_y_ao,super_phi_y_);
    TransformPhiMatrixAOMO(super_phi_z_ao,super_phi_z_);

    outfile->Printf(" Done. <==\n"); 
    outfile->Printf("\n");

}// end of common_init()

void RealSpaceDensity::GetGridInfo() {

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

void RealSpaceDensity::BuildPhiMatrixAO(std::string phi_type, std::shared_ptr<Matrix> myphi) {

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

void RealSpaceDensity::TransformPhiMatrixAOMO(std::shared_ptr<Matrix> phi_in, std::shared_ptr<Matrix> phi_out) {
        
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

void RealSpaceDensity::SetD1(std::vector<opdm> my_opdm, std::shared_ptr<Matrix> D1) {

    for (size_t n = 0; n < my_opdm.size(); n++) {

        int i = my_opdm[n].i;
        int j = my_opdm[n].j;

        int hi = symmetry_[i];
        int hj = symmetry_[j];

        if ( hi != hj ) {
            throw PsiException("error: something is wrong with the symmetry of the OPDM",__FILE__,__LINE__);
        }

        int ii = i - pitzer_offset_[hi];
        int jj = j - pitzer_offset_[hi];

        D1->pointer(hi)[ii][jj] = my_opdm[n].value;
    }
}

void RealSpaceDensity::SetOPDM(std::vector<opdm> opdm_a, std::vector<opdm> opdm_b) {

    opdm_a_ = opdm_a;
    opdm_b_ = opdm_b;

    SetD1(opdm_a_, Da_);
    SetD1(opdm_b_, Db_);
}

void RealSpaceDensity::SetTPDM(std::vector<tpdm> tpdm, std::string tpdm_type) {

    std::transform(tpdm_type.begin(), tpdm_type.end(), tpdm_type.begin(),
        [](unsigned char c){ return std::tolower(c); });

    if (tpdm_type == "ab"){
        tpdm_ab_ = tpdm;
    }else if (tpdm_type == "aa"){
        tpdm_aa_ = tpdm;
    }else if (tpdm_type == "bb"){
        tpdm_bb_ = tpdm;
    }else {
        throw PsiException("invalid tpdm_type",__FILE__,__LINE__);
    }

}

void RealSpaceDensity::ReadOPDM() {
    
    // read opdm from disk 

    std::shared_ptr<PSIO> psio (new PSIO());

    if ( !psio->exists(PSIF_V2RDM_D1A) ) throw PsiException("No D1a on disk",__FILE__,__LINE__);
    if ( !psio->exists(PSIF_V2RDM_D1B) ) throw PsiException("No D1b on disk",__FILE__,__LINE__);

    // D1a

    psio->open(PSIF_V2RDM_D1A,PSIO_OPEN_OLD);

    long int na;
    psio->read_entry(PSIF_V2RDM_D1A,"length",(char*)&na,sizeof(long int));

    opdm_a_.resize(na);
    psio->read_entry(PSIF_V2RDM_D1A,"D1a",(char*)&opdm_a_[0],na * sizeof(opdm));
    psio->close(PSIF_V2RDM_D1A,1);

    SetD1(opdm_a_, Da_);

    // D1b

    psio->open(PSIF_V2RDM_D1B,PSIO_OPEN_OLD);

    long int nb;
    psio->read_entry(PSIF_V2RDM_D1B,"length",(char*)&nb,sizeof(long int));

    opdm_b_.resize(nb);
    psio->read_entry(PSIF_V2RDM_D1B,"D1b",(char*)&opdm_b_[0],nb * sizeof(opdm));
    psio->close(PSIF_V2RDM_D1B,1);

    SetD1(opdm_b_, Db_);

}

void RealSpaceDensity::BuildExchangeCorrelationHole(size_t p) {

    outfile->Printf("\n");
    outfile->Printf("    ==> Build Exchange-Correlation Hole ...");

    if (tpdm_aa_.size() == 0){
        throw PsiException("no aa block of tpdm. build or set tpdm.",__FILE__,__LINE__);
    }else if (tpdm_ab_.size() == 0){
        throw PsiException("no ab block of tpdm. build or set tpdm.",__FILE__,__LINE__);
    }else if (tpdm_bb_.size() == 0){
        throw PsiException("no bb block of tpdm. build or set tpdm.",__FILE__,__LINE__);
    }

    double integral = 0.0;
    double * rho_p = rho_->pointer();
    double * w_p   = grid_w_->pointer();

    xc_hole_.reset();
    xc_hole_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    double * xc_hole_p = xc_hole_->pointer();
    for (int q = 0; q < phi_points_; q++) {

        // pi(r,r') = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r',nu) * phi(r,lambda) * phi(r',sigma)
        double pi = 0.0;
        for (size_t n = 0; n < tpdm_ab_.size(); n++) {

            int i = tpdm_ab_[n].i;
            int j = tpdm_ab_[n].j;
            int k = tpdm_ab_[n].k;
            int l = tpdm_ab_[n].l;

            int hi = symmetry_[i];
            int hj = symmetry_[j];
            int hk = symmetry_[k];
            int hl = symmetry_[l];

            int ii = i - pitzer_offset_[hi];
            int jj = j - pitzer_offset_[hj];
            int kk = k - pitzer_offset_[hk];
            int ll = l - pitzer_offset_[hl];

            pi += super_phi_->pointer(hi)[p][ii] * 
                  super_phi_->pointer(hj)[q][jj] * 
                  super_phi_->pointer(hk)[p][kk] * 
                  super_phi_->pointer(hl)[q][ll] * tpdm_ab_[n].value * 0.5;

            pi += super_phi_->pointer(hi)[q][ii] * 
                  super_phi_->pointer(hj)[p][jj] * 
                  super_phi_->pointer(hk)[q][kk] * 
                  super_phi_->pointer(hl)[p][ll] * tpdm_ab_[n].value * 0.5;

        }
        for (size_t n = 0; n < tpdm_aa_.size(); n++) {

            int i = tpdm_aa_[n].i;
            int j = tpdm_aa_[n].j;
            int k = tpdm_aa_[n].k;
            int l = tpdm_aa_[n].l;

            int hi = symmetry_[i];
            int hj = symmetry_[j];
            int hk = symmetry_[k];
            int hl = symmetry_[l];

            int ii = i - pitzer_offset_[hi];
            int jj = j - pitzer_offset_[hj];
            int kk = k - pitzer_offset_[hk];
            int ll = l - pitzer_offset_[hl];

            pi += super_phi_->pointer(hi)[p][ii] * 
                  super_phi_->pointer(hj)[q][jj] * 
                  super_phi_->pointer(hk)[p][kk] * 
                  super_phi_->pointer(hl)[q][ll] * tpdm_aa_[n].value * 0.5;

        }
        for (size_t n = 0; n < tpdm_bb_.size(); n++) {

            int i = tpdm_bb_[n].i;
            int j = tpdm_bb_[n].j;
            int k = tpdm_bb_[n].k;
            int l = tpdm_bb_[n].l;

            int hi = symmetry_[i];
            int hj = symmetry_[j];
            int hk = symmetry_[k];
            int hl = symmetry_[l];

            int ii = i - pitzer_offset_[hi];
            int jj = j - pitzer_offset_[hj];
            int kk = k - pitzer_offset_[hk];
            int ll = l - pitzer_offset_[hl];

            pi += super_phi_->pointer(hi)[p][ii] * 
                  super_phi_->pointer(hj)[q][jj] * 
                  super_phi_->pointer(hk)[p][kk] * 
                  super_phi_->pointer(hl)[q][ll] * tpdm_bb_[n].value * 0.5;

        }
        double nxc = (2.0 * pi - rho_p[p] * rho_p[q]) / rho_p[p];
        //if ( fabs(grid_x_->pointer()[q] < 1e-6 ) ) {
        //printf("%20.12lf %20.12lf %20.12lf %20.12lf\n",grid_x_->pointer()[q],grid_y_->pointer()[q],grid_z_->pointer()[q],nxc);
        //}

        xc_hole_p[q] = nxc;

        integral += nxc * w_p[q];

    }

    outfile->Printf(" Done. <==\n");

    //printf("integral over xc hole (should be -1): %20.12lf\n",integral);
    //exit(0);
}

// build pi. note that there is a low-memory version of this
// function in edeprince3/real_space_density.git, should we ever need it
void RealSpaceDensity::BuildPiFast(std::vector<tpdm> D2ab) {

    outfile->Printf("\n");
    outfile->Printf("    ==> Build Pi ...");

    pi_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double * pi_p = pi_->pointer();

    for (int p = 0; p < phi_points_; p++) {

        double dum = 0.0;

        // pi(r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
        for (size_t n = 0; n < D2ab.size(); n++) {

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
        for (size_t n = 0; n < D2ab.size(); n++) {

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
    outfile->Printf(" Done. <==\n\n");
}

void RealSpaceDensity::BuildRhoFast(){

    outfile->Printf("\n");
    outfile->Printf("    ==> Build Rho ...\n");

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
        for (size_t n = 0; n < opdm_a_.size(); n++) {

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
        for (size_t n = 0; n < opdm_b_.size(); n++) {

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

    rho_a_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    rho_b_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    rho_a_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    rho_b_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    rho_a_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    rho_b_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double ** phi_x = super_phi_x_->pointer();
    double ** phi_y = super_phi_y_->pointer();
    double ** phi_z = super_phi_z_->pointer();

    double * rho_a_xp = rho_a_x_->pointer();
    double * rho_b_xp = rho_b_x_->pointer();

    double * rho_a_yp = rho_a_y_->pointer();
    double * rho_b_yp = rho_b_y_->pointer();

    double * rho_a_zp = rho_a_z_->pointer();
    double * rho_b_zp = rho_b_z_->pointer();

    for (int p = 0; p < phi_points_; p++) {

        // rho'_a(r)

        double duma_x = 0.0;
        double duma_y = 0.0;
        double duma_z = 0.0;
        double dumta = 0.0;

        for (size_t n = 0; n < opdm_a_.size(); n++) {

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

        for (size_t n = 0; n < opdm_b_.size(); n++) {
            
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
    }

    outfile->Printf("    ... Done. <==\n");
    outfile->Printf("\n");
}

void RealSpaceDensity::ReadTPDM(std::string tpdm_type) {

    std::transform(tpdm_type.begin(), tpdm_type.end(), tpdm_type.begin(),
        [](unsigned char c){ return std::tolower(c); });

    std::shared_ptr<PSIO> psio (new PSIO());

    if (tpdm_type == "ab"){

        if ( !psio->exists(PSIF_V2RDM_D2AB) ) throw PsiException("No D2ab on disk",__FILE__,__LINE__);
        
        psio->open(PSIF_V2RDM_D2AB,PSIO_OPEN_OLD);

        long int n;
        psio->read_entry(PSIF_V2RDM_D2AB, "length", (char*)&n, sizeof(long int));

        tpdm_ab_.resize(n);
        memset((void*)&tpdm_ab_[0], '\0',n * sizeof(tpdm));

        psio->read_entry(PSIF_V2RDM_D2AB, "D2ab",(char*)&tpdm_ab_[0], n * sizeof(tpdm));
        psio->close(PSIF_V2RDM_D2AB,1);

    }else if (tpdm_type == "aa"){

        if ( !psio->exists(PSIF_V2RDM_D2AA) ) throw PsiException("No D2aa on disk",__FILE__,__LINE__);
        
        psio->open(PSIF_V2RDM_D2AA,PSIO_OPEN_OLD);

        long int n;
        psio->read_entry(PSIF_V2RDM_D2AA, "length", (char*)&n, sizeof(long int));

        tpdm_aa_.resize(n);
        memset((void*)&tpdm_aa_[0], '\0',n * sizeof(tpdm));

        psio->read_entry(PSIF_V2RDM_D2AA, "D2aa",(char*)&tpdm_aa_[0], n * sizeof(tpdm));
        psio->close(PSIF_V2RDM_D2AA,1);

    }else if (tpdm_type == "bb"){

        if ( !psio->exists(PSIF_V2RDM_D2BB) ) throw PsiException("No D2bb on disk",__FILE__,__LINE__);
        
        psio->open(PSIF_V2RDM_D2BB,PSIO_OPEN_OLD);

        long int n;
        psio->read_entry(PSIF_V2RDM_D2BB, "length", (char*)&n, sizeof(long int));

        tpdm_bb_.resize(n);
        memset((void*)&tpdm_bb_[0], '\0',n * sizeof(tpdm));

        psio->read_entry(PSIF_V2RDM_D2BB, "D2bb",(char*)&tpdm_bb_[0], n * sizeof(tpdm));
        psio->close(PSIF_V2RDM_D2BB,1);

    }else {
        throw PsiException("invalid tpdm_type",__FILE__,__LINE__);
    }
}

std::shared_ptr<Vector> RealSpaceDensity::xc_hole(double x, double y, double z) {

    // find appropriate grid point
    double * x_p = grid_x_->pointer();
    double * y_p = grid_y_->pointer();
    double * z_p = grid_z_->pointer();

    int p = 0;
    double min = 9e9;
    for (int myp = 0; myp < phi_points_; myp++) {
        double dx = x_p[myp] - x;
        double dy = y_p[myp] - y;
        double dz = z_p[myp] - z;
        double r2 = dx*dx + dy*dy + dz*dz;
        if ( r2 < min ) {
            min = r2;
            p = myp;
        }
    }
    BuildExchangeCorrelationHole(p);
    return xc_hole_;
}

std::shared_ptr<Vector> RealSpaceDensity::slater_potential(){

    std::shared_ptr<Vector> vs (new Vector("Slater potential", phi_points_));

    double * x_p = grid_x_->pointer();
    double * y_p = grid_y_->pointer();
    double * z_p = grid_z_->pointer();
    double * w_p = grid_w_->pointer();

    for (int p = 0; p < phi_points_; p++) {
        double x = x_p[p];
        double y = y_p[p];
        if ( fabs(x) > 1e-6 || fabs(y) > 1e-6 ){
            continue;
        }
        BuildExchangeCorrelationHole(p);
        double * xc_hole_p = xc_hole_->pointer();
        double my_vs = 0.0;
        for (int q = 0; q < phi_points_; q++) {
            if (p != q) {
                double dx = x_p[p] - x_p[q];
                double dy = y_p[p] - y_p[q];
                double dz = z_p[p] - z_p[q];

                double r2 = dx*dx + dy*dy + dz*dz;

                my_vs +=  w_p[q] * xc_hole_p[q] / sqrt(r2);
            }else {
                //my_vs += w_p[q] * xc_hole_p[q] * 2.0 * M_PI * pow(3.0 * w_p[q] / (4.0 * M_PI), 2.0 / 3.0);
            }
        }
        vs->pointer()[p] = my_vs;
    }

    return vs;
}

/// compute and return closest Kohn-Sham orbitals to a given reference density
void RealSpaceDensity::ks_orbitals(){

    if ( nirrep_ > 1 ) {
        throw PsiException("cannot determine ks orbitals for symmetry other than c1", __FILE__,__LINE__);
    }

    // grab the kinetic energy integrals from MintsHelper:
    std::shared_ptr<MintsHelper> mints (new MintsHelper(reference_wavefunction_));

    // one-electron kinetic energy integrals
    std::shared_ptr<Matrix> T = mints->so_kinetic();

    // MO basis
    T->transform(Ca_);
    Ca_->print_out();

    std::shared_ptr<cs_solver> cs (new cs_solver(options_));

    std::vector<double> T_vec(T->pointer()[0], T->pointer()[0] + nso_*nso_);
    std::vector<double> phi_vec(super_phi_->pointer()[0], super_phi_->pointer()[0] + phi_points_*nso_);
    std::vector<double> rho_vec(rho_->pointer(), rho_->pointer() + phi_points_);
    std::vector<double> x_vec(grid_x_->pointer(), grid_x_->pointer() + phi_points_);
    std::vector<double> y_vec(grid_y_->pointer(), grid_y_->pointer() + phi_points_);
    std::vector<double> z_vec(grid_z_->pointer(), grid_z_->pointer() + phi_points_);
    std::vector<double> w_vec(grid_w_->pointer(), grid_w_->pointer() + phi_points_);

    cs->solve(nalpha_, 
              nbeta_, 
              T_vec,
              phi_vec,
              rho_vec,
              {x_vec, y_vec, z_vec, w_vec});

    std::vector<std::shared_ptr<Matrix>> C;
    //C.push_back(Ca_);
    //C.push_back(Cb_);
    //return C;
}

} //end namespaces
