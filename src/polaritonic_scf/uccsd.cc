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

#include <psi4/libtrans/integraltransform.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/mintshelper.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/basisset.h>
#include <psi4/lib3index/dftensor.h>
#include <psi4/libqt/qt.h>

#include "uccsd.h"

#include <misc/blas.h>
#include <misc/hilbert_psifiles.h>

using namespace psi;
using namespace fnocc;

namespace hilbert{ 

PolaritonicUCCSD::PolaritonicUCCSD(std::shared_ptr<Wavefunction> reference_wavefunction, Options& options_):
    PolaritonicHF(reference_wavefunction,options_) {
    common_init();
}

PolaritonicUCCSD::~PolaritonicUCCSD() {

    if ( !is_hubbard_ ) {
        free(Qoo_);
        free(Qov_);
        free(Qvo_);
        free(Qvv_);
    }else {
        free(eri_abcd_);
        free(eri_abic_);
        free(eri_aibc_);
        free(eri_abij_);
    }

    free(eri_ijkl_);
    free(eri_iajb_);
    free(eri_ijab_);
    free(eri_jkia_);
    free(eri_iajk_);

    free(tmp1_);
    free(tmp2_);
    free(tmp3_);
    free(residual_);
    free(ccamps_);

}

void PolaritonicUCCSD::common_init() {

    outfile->Printf("\n\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    //outfile->Printf( "        *    Polaritonic UCCSD                                *\n");
    outfile->Printf( "        *    UCCSD                                            *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *******************************************************\n");

    is_hubbard_ = options_.get_bool("HUBBARD_HAMILTONIAN");

    outfile->Printf("\n");
    outfile->Printf("    ==> Hamiltonian type <==\n");
    outfile->Printf("\n");
    outfile->Printf("        %s\n",is_hubbard_ ? "Hubbard" : "molecular");
    outfile->Printf("\n");

    // update cavity terms once more
    if ( n_photon_states_ > 1 ) {
        update_cavity_terms();
    }

    same_a_b_orbs_ = false;
    same_a_b_dens_ = false;

    // include amplitudes for photon transitions?
    include_u0_ = options_.get_bool("POLARITONIC_CC_INCLUDE_U0");
    include_u1_ = options_.get_bool("POLARITONIC_CC_INCLUDE_U1");
    include_u2_ = options_.get_bool("POLARITONIC_CC_INCLUDE_U2");

    // molecular hamiltonian
    if ( !is_hubbard_ ) {

        initialize_with_molecular_hamiltonian();

    }else {

        initialize_with_hubbard_hamiltonian();

    }

    // allocate memory for amplitudes, residual, and temporary buffers

    ccamps_   = (double*)malloc(ccamps_dim_*sizeof(double));
    residual_ = (double*)malloc(ccamps_dim_*sizeof(double));

    size_t off = 0;
    t2_  = ccamps_   + off; 
    rt2_ = residual_ + off; off += o_*o_*v_*v_;

    t1_  = ccamps_   + off; 
    rt1_ = residual_ + off; off += o_*v_;

    if ( include_u0_ ) {
        u0_    = ccamps_   + off;
        ru0_   = residual_ + off; off++;
    }
    if ( include_u1_ ) {
        u1_    = ccamps_   + off;
        ru1_   = residual_ + off; off += o_*v_;
    }
    if ( include_u2_ ) {
        u2_    = ccamps_   + off;
        ru2_   = residual_ + off; off += o_*o_*v_*v_;
    }

    memset((void*)ccamps_,'\0',ccamps_dim_*sizeof(double));
    memset((void*)residual_,'\0',ccamps_dim_*sizeof(double));

    // initialize diis solver
    diis = (std::shared_ptr<DIIS>)(new DIIS(ccamps_dim_));

}

void PolaritonicUCCSD::initialize_with_hubbard_hamiltonian() {

    if ( ( options_.get_int("N_HUBBARD_SPINS") % 2 ) != 0 ) {
        throw PsiException("Hubbard model only works for even numbers of electrons currently.",__FILE__,__LINE__);
    }

    energy_       = 0.0;

    nalpha_       = options_.get_int("N_HUBBARD_SPINS") / 2;
    nbeta_        = options_.get_int("N_HUBBARD_SPINS") / 2;

    nso_          = options_.get_int("N_HUBBARD_SITES");
    nmo_          = options_.get_int("N_HUBBARD_SITES");

    size_t ms = (int)nalpha_ - (int)nbeta_;
    multiplicity_ = 2 * ms + 1;

    nirrep_       = 1;

    o_ = nalpha_ + nbeta_;
    v_ = (nmo_-nalpha_) + (nmo_-nbeta_);

    ccamps_dim_ = o_*o_*v_*v_ + o_*v_;

    if ( include_u0_ ) {
        ccamps_dim_++;
    }
    if ( include_u1_ ) {
        ccamps_dim_ += o_*v_;
    }
    if ( include_u2_ ) {
        ccamps_dim_ += o_*o_*v_*v_;
    }

    // temporary storage ... reduce later

    int dim = o_*o_*v_*v_;
    if ( o_ > v_ ) {
        dim = o_*o_*o_*o_;
    }

    tmp1_ = (double*)malloc(dim*sizeof(double));
    tmp2_ = (double*)malloc(dim*sizeof(double));
    tmp3_ = (double*)malloc(dim*sizeof(double));
    memset((void*)tmp1_,'\0',dim*sizeof(double));
    memset((void*)tmp2_,'\0',dim*sizeof(double));
    memset((void*)tmp3_,'\0',dim*sizeof(double));


    // alpha + beta MO transformation matrix (identity)
    C_ = (std::shared_ptr<Matrix>)(new Matrix(2L*nso_,2L*nmo_));
    double ** cp = C_->pointer();

    // solve hubbard scf problem
    std::shared_ptr<Matrix> Ca = hubbard_hartree_fock();
    double ** ca = Ca->pointer();

    for (size_t mu = 0; mu < nso_; mu++) {
        size_t count = 0;
        for (size_t i = 0; i < nalpha_; i++) {
            cp[mu][count++] = ca[mu][i];
        }
        for (size_t i = 0; i < nbeta_; i++) {
            cp[mu+nso_][count++] = ca[mu][i];
        }
        for (size_t i = nalpha_; i < nmo_; i++) {
            cp[mu][count++] = ca[mu][i];
        }
        for (size_t i = nbeta_; i < nmo_; i++) {
            cp[mu+nso_][count++] = ca[mu][i];
        }
    }

    // core hamiltonian
    std::shared_ptr<Matrix> Ha = (std::shared_ptr<Matrix>)(new Matrix(nso_,nso_));
    double ** ha_p = Ha->pointer();

    H_ = (std::shared_ptr<Matrix>)(new Matrix(2*nso_,2*nso_));
    double ** h_p = H_->pointer();
    H_->zero();

    double t = options_.get_double("HUBBARD_T");

    for (size_t i = 0; i < nso_ - 1; i++) {
        ha_p[i][i+1] = -t;
        ha_p[i+1][i] = -t;
    }

    for (size_t mu = 0; mu < nso_; mu++) {
        for (size_t nu = 0; nu < nso_; nu++) {
            h_p[mu][nu]           = ha_p[mu][nu];
            h_p[mu+nso_][nu+nso_] = ha_p[mu][nu];
        }
    }

    Dipole_x_ = (std::shared_ptr<Matrix>)(new Matrix(2*nso_,2*nso_));
    Dipole_y_ = (std::shared_ptr<Matrix>)(new Matrix(2*nso_,2*nso_));
    Dipole_z_ = (std::shared_ptr<Matrix>)(new Matrix(2*nso_,2*nso_));

    std::shared_ptr<Matrix> dip_a (new Matrix(nso_,nso_));
    if ( nso_ != 4 ) {
        throw PsiException("dipole terms hard-coded for four site Hubbard model",__FILE__,__LINE__);
    }
    double ** da_p = dip_a->pointer();
    da_p[0][0] = -1.5;
    da_p[1][1] = -0.5;
    da_p[2][2] =  0.5;
    da_p[3][3] =  1.5;

    double ** dz = Dipole_z_->pointer();
    for (size_t mu = 0; mu < nso_; mu++) {
        for (size_t nu = 0; nu < nso_; nu++) {
            dz[mu][nu]           = da_p[mu][nu];
            dz[mu+nso_][nu+nso_] = da_p[mu][nu];
        }
    }

    // construct four-index integrals (and build fock matrix)
    F_ = (std::shared_ptr<Matrix>)(new Matrix(H_));

    // orbital energies
    epsilon_ = (double*)malloc(2*nso_*sizeof(double));
    memset((void*)epsilon_,'\0',2*nso_*sizeof(double));

    // build eris
    t1_transformation_hubbard_hamiltonian(C_,C_,true);

}

std::shared_ptr<Matrix> PolaritonicUCCSD::hubbard_hartree_fock() {

    // assuming closed shell

    std::shared_ptr<Matrix> Ca = (std::shared_ptr<Matrix>)(new Matrix(nso_,nmo_));
    std::shared_ptr<Matrix> Da = (std::shared_ptr<Matrix>)(new Matrix(nso_,nso_));
    std::shared_ptr<Matrix> Fa = (std::shared_ptr<Matrix>)(new Matrix(nso_,nso_));
    std::shared_ptr<Matrix> Ha = (std::shared_ptr<Matrix>)(new Matrix(nso_,nso_));

    std::shared_ptr<Vector> eps = (std::shared_ptr<Vector>)(new Vector(nso_));

    double ** cp = Ca->pointer();
    double ** dp = Da->pointer();
    double ** fp = Fa->pointer();
    double ** hp = Ha->pointer();

    double t = options_.get_double("HUBBARD_T");

    for (size_t i = 0; i < nso_ - 1; i++) {
        hp[i][i+1] = -t;
        hp[i+1][i] = -t;
    }

    // form F = H
    Fa->copy(Ha);

    // diagonalize F' to obtain C'
    Fa->diagonalize(Ca,eps,ascending);

    // Construct density from C
    C_DGEMM('n','t',nso_,nso_,nalpha_,1.0,&(Ca->pointer()[0][0]),nso_,&(Ca->pointer()[0][0]),nso_,0.0,&(Da->pointer()[0][0]),nso_);

    double e_last    = 0.0;
    double dele      = 0.0;

    double energy = 0.0;
    energy += Da->vector_dot(Ha);
    energy += Da->vector_dot(Fa);

    outfile->Printf("    ==>  Begin Hubbard SCF Iterations <==\n");
    outfile->Printf("\n");
    outfile->Printf("    Guess energy:  %20.12lf\n",energy);
    outfile->Printf("\n");
    outfile->Printf("    ");
    outfile->Printf(" Iter ");
    outfile->Printf("              energy ");
    outfile->Printf("                  dE ");
    outfile->Printf("\n");

    // two-electron integrals
    double * eri = (double*)malloc(nso_*nso_*nso_*nso_*sizeof(double));
    memset((void*)eri,'\0',nso_*nso_*nso_*nso_*sizeof(double));

    double U = options_.get_double("HUBBARD_U");

    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

    std::shared_ptr<Matrix> dip_a (new Matrix(nso_,nso_));
    double ** da_p = dip_a->pointer();
    da_p[0][0] = -1.5;
    da_p[1][1] = -0.5;
    da_p[2][2] =  0.5;
    da_p[3][3] =  1.5;

    for (size_t p = 0; p < nso_; p++) {
        eri[p*nso_*nso_*nso_+p*nso_*nso_+p*nso_+p] = U;
    }
    // dipole self energy
    for (size_t p = 0; p < nso_; p++) {
        for (size_t q = 0; q < nso_; q++) {
            eri[p*nso_*nso_*nso_+p*nso_*nso_+q*nso_+q] += 1.0 * da_p[p][p] * da_p[q][q] * lambda_z * lambda_z;
        }
    }

    double e_convergence = options_.get_double("E_CONVERGENCE");
    size_t maxiter       = options_.get_int("MAXITER");

    int iter = 0;
    do {

        e_last = energy;

        // F = H + 2 J - K
        Fa->copy(Ha);
        for (size_t p = 0; p < nso_; p++) {
            for (size_t q = 0; q < nso_; q++) {
                double dum = 0.0;
                for (size_t r = 0; r < nso_; r++) {
                    for (size_t s = 0; s < nso_; s++) {
                        dum += 2.0 * dp[r][s] * eri[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s];
                        dum -=       dp[r][s] * eri[p*nso_*nso_*nso_+s*nso_*nso_+r*nso_+q];
                    }
                }
                // dipole self energy
                double val = 0.0;
                for (size_t r = 0; r < nso_; r++) {
                    val += da_p[p][r] * da_p[r][q];
                }
                fp[p][q] += dum + 0.5 * val * lambda_z * lambda_z;
            }
        }

        // cavity-molecule terms

        double dipole = 2.0 * Da->vector_dot(dip_a);

        // cavity 

        // bare cavity Hamiltonian
        std::shared_ptr<Matrix> cavity_H (new Matrix(n_photon_states_,n_photon_states_));
        double ** cavity_Hp = cavity_H->pointer();
        for (size_t A = 0; A < n_photon_states_; A++) {
            cavity_Hp[A][A] = A * cavity_frequency_[2];
        }

        // cavity dipole operator
        std::shared_ptr<Matrix> cavity_D (new Matrix(n_photon_states_,n_photon_states_));
        double ** cavity_Dp = cavity_D->pointer();
        for (int A=0; A<n_photon_states_-1; A++){
            cavity_Dp[A+1][A] += cavity_frequency_[2] * cavity_coupling_strength_[2] * sqrt(A+1);
            cavity_Dp[A][A+1] += cavity_frequency_[2] * cavity_coupling_strength_[2] * sqrt(A+1);
        }

        // cavity interaction
        std::shared_ptr<Matrix> cavity_I (new Matrix(n_photon_states_,n_photon_states_));
        double ** cavity_Ip = cavity_I->pointer();
        for (int A = 0; A<n_photon_states_; A++) {
            if (A < n_photon_states_ - 1) {
                cavity_Ip[A][A+1] = -dipole * cavity_frequency_[2] * cavity_coupling_strength_[2] * sqrt(A+1);
            }
            if (A > 0) {
                cavity_Ip[A][A-1] = -dipole * cavity_frequency_[2] * cavity_coupling_strength_[2] * sqrt(A);
            }
        }

        // total cavity hamiltonian
        cavity_H->add(cavity_I);

        // diagonalize total cavity hamiltonian
        std::shared_ptr<Matrix> eigvec (new Matrix(n_photon_states_,n_photon_states_));
        std::shared_ptr<Vector> eigval (new Vector(n_photon_states_));
        cavity_H->diagonalize(eigvec, eigval);

        // transform cavity dipole matrix to basis that diagonalizes total cavity hamiltonian
        cavity_D->transform(eigvec);

        // add cavity interaction term to fock matrix
        Fa->axpy(-cavity_Dp[0][0] * cavity_frequency_[2] * cavity_coupling_strength_[2],dip_a);

        // diagonalize F' to obtain C'
        Fa->diagonalize(Ca,eps,ascending);

        // Construct density from C
        C_DGEMM('n','t',nso_,nso_,nalpha_,1.0,&(Ca->pointer()[0][0]),nso_,&(Ca->pointer()[0][0]),nso_,0.0,&(Da->pointer()[0][0]),nso_);

        energy = 0.0;
        energy += Da->vector_dot(Ha);
        energy += -2.0 * cavity_Dp[0][0] * cavity_frequency_[2] * cavity_coupling_strength_[2] * Da->vector_dot(dip_a);
        energy += Da->vector_dot(Fa);

        // dele
        dele = energy - e_last;

        outfile->Printf("    %5i %20.12lf %20.12lf\n",iter,energy,dele);

        iter++;
        if ( iter > maxiter ) break;

    }while( fabs(dele) > e_convergence );

    free(eri);

    if ( iter > maxiter ) {
        throw PsiException("Maximum number of iterations exceeded!",__FILE__,__LINE__);
    }

    outfile->Printf("\n");
    outfile->Printf("    Hubbard SCF iterations converged!\n");
    outfile->Printf("\n");

    outfile->Printf("    * Hubbard SCF total energy: %20.12lf\n",energy);

    return Ca;
}

void PolaritonicUCCSD::initialize_with_molecular_hamiltonian() {

    // ensure scf_type df
    if ( options_.get_str("SCF_TYPE") != "DF" ) {
        throw PsiException("polaritonic uhf only works with scf_type df for now",__FILE__,__LINE__);
    }

    // ensure running in c1 symmetry
    if ( reference_wavefunction_->nirrep() > 1 ) {
        throw PsiException("polaritonic uccsd only works with c1 symmetry for now.",__FILE__,__LINE__);
    }

    o_ = nalpha_ + nbeta_;
    v_ = (nmo_-nalpha_) + (nmo_-nbeta_);

    ccamps_dim_ = o_*o_*v_*v_ + o_*v_;

    if ( include_u0_ ) {
        ccamps_dim_++;
    }
    if ( include_u1_ ) {
        ccamps_dim_ += o_*v_;
    }
    if ( include_u2_ ) {
        ccamps_dim_ += o_*o_*v_*v_;
    }

    // memory requirements:

    // get auxiliary basis:
    std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_CC");

    // total number of auxiliary basis functions
    nQ_ = auxiliary->nbf();

    size_t required_memory = 0;
    required_memory += ccamps_dim_ * 4L;         // amplitudes, residual, diis storage
    //required_memory += o_*v_*v_*v_;              // <ai||bc>
    required_memory += o_*o_*v_*v_;              // <ia||jb>
    required_memory += o_*o_*v_*v_;              // <ij||ab>
    required_memory += o_*o_*o_*v_;              // <jk||ia>
    required_memory += o_*o_*o_*v_;              // <ia||jk>
    required_memory += o_*o_*o_*o_;              // <ij||kl>
    required_memory += o_*o_*v_*v_ * 2L;         // tmp2_, tmp3_
    required_memory += nQ_*(o_+v_)*(o_+v_) * 2L; // three-index integrals: Qmo, Qoo, Qov, Qvo, Qvv
    //required_memory += v_*v_*v_*v_;              // <ab||cd>
    //required_memory += o_*v_*v_*v_;              // <ab||ic>
    //required_memory += o_*o_*v_*v_;              // <ab||ij>

    // tmp1 used to be ov^3, but we have eliminated eri_aibc_...
    size_t dim = o_*o_*v_*v_; 
    if ( o_ > v_ ) {
        dim = o_*o_*o_*o_;
    }
    if ( 2L*o_*v_*nQ_ > dim ) {
        dim = 2L*o_*v_*nQ_;
    }
    if ( 2L*v_*v_*v_ > dim ) {
        dim = 2L*v_*v_*v_;
    }
    size_t ns = 2L * (size_t)nso_;
    if ( nQ_*ns*ns > dim ) {
        dim = nQ_*ns*ns;
    }
    required_memory += dim;                      // tmp1_

    outfile->Printf("\n");
    outfile->Printf("    Required memory:          %8.1lf mb\n",(double)required_memory*8L/(1024.0*1024.0));
    outfile->Printf("\n");

    // available memory
    memory_ = Process::environment.get_memory();

    if ( 8L * required_memory > memory_ ) {
        throw PsiException("not enough memory",__FILE__,__LINE__);
    }

    // temporary storage ... reduce later

    tmp1_ = (double*)malloc(dim*sizeof(double));
    memset((void*)tmp1_,'\0',dim*sizeof(double));

    tmp2_ = (double*)malloc(o_*o_*v_*v_*sizeof(double));
    tmp3_ = (double*)malloc(o_*o_*v_*v_*sizeof(double));
    memset((void*)tmp2_,'\0',o_*o_*v_*v_*sizeof(double));
    memset((void*)tmp3_,'\0',o_*o_*v_*v_*sizeof(double));

    // three-index integral containers
    Qoo_ = (double*)malloc(nQ_*o_*o_*sizeof(double));
    Qov_ = (double*)malloc(nQ_*o_*v_*sizeof(double));
    Qvo_ = (double*)malloc(nQ_*v_*o_*sizeof(double));
    Qvv_ = (double*)malloc(nQ_*v_*v_*sizeof(double));

    memset((void*)Qoo_,'\0',nQ_*o_*o_*sizeof(double));
    memset((void*)Qov_,'\0',nQ_*o_*v_*sizeof(double));
    memset((void*)Qvo_,'\0',nQ_*v_*o_*sizeof(double));
    memset((void*)Qvv_,'\0',nQ_*v_*v_*sizeof(double));

    // alpha + beta MO transformation matrix
    C_ = (std::shared_ptr<Matrix>)(new Matrix(2*nso_,2*nmo_));
    C_->zero();
    double ** cp = C_->pointer();
    double ** ca = Ca_->pointer();
    double ** cb = Cb_->pointer();

    for (size_t mu = 0; mu < nso_; mu++) {
        size_t count = 0;
        for (size_t i = 0; i < nalpha_; i++) {
            cp[mu][count++] = ca[mu][i];
        }
        for (size_t i = 0; i < nbeta_; i++) {
            cp[mu+nso_][count++] = cb[mu][i];
        }
        for (size_t i = nalpha_; i < nmo_; i++) {
            cp[mu][count++] = ca[mu][i];
        }
        for (size_t i = nbeta_; i < nmo_; i++) {
            cp[mu+nso_][count++] = cb[mu][i];
        }
    }

    // core hamiltonian
    H_ = (std::shared_ptr<Matrix>)(new Matrix(2*nso_,2*nso_));
    H_->zero();
    double ** h_p = H_->pointer();

    auto coreH = reference_wavefunction_->H()->clone();
    double ** coreH_p = coreH->pointer();

    for (size_t mu = 0; mu < nso_; mu++) {
        for (size_t nu = 0; nu < nso_; nu++) {
            double dum = coreH_p[mu][nu];
            h_p[mu][nu]           = dum;
            h_p[mu+nso_][nu+nso_] = dum;
        }
    }

    // extra one-electron terms introduced by cavity
    oe_cavity_terms_ = (std::shared_ptr<Matrix>)(new Matrix(2*nso_,2*nso_));
    oe_cavity_terms_->zero();
    double ** oe_cavity_terms_p = oe_cavity_terms_->pointer();

    // e-(n-<d>) contribution to dipole self energy
    double ** scaled_e_n_dipole_squared_p = scaled_e_n_dipole_squared_->pointer();

    // one-electron part of electron-electron contribution to dipole self energy
    double ** quadrupole_scaled_sum_p = quadrupole_scaled_sum_->pointer();

    for (size_t mu = 0; mu < nso_; mu++) {
        for (size_t nu = 0; nu < nso_; nu++) {
            double dum = scaled_e_n_dipole_squared_p[mu][nu] - quadrupole_scaled_sum_p[mu][nu];
            oe_cavity_terms_p[mu][nu]           = dum;
            oe_cavity_terms_p[mu+nso_][nu+nso_] = dum;
        }
    }

    // fock matrix
    F_ = (std::shared_ptr<Matrix>)(new Matrix(2*nso_,2*nso_));
    F_->zero();
    double ** fp = F_->pointer();
    double ** fa = Fa_->pointer();
    double ** fb = Fb_->pointer();
    for (size_t mu = 0; mu < nso_; mu++) {
        for (size_t nu = 0; nu < nso_; nu++) {
            fp[mu][nu]           = fa[mu][nu];
            fp[mu+nso_][nu+nso_] = fb[mu][nu];
        }
    }
    F_->transform(C_);

    // scaled sum of dipole integrals
    Dipole_x_ = (std::shared_ptr<Matrix>)(new Matrix(2*nso_,2*nso_));
    Dipole_y_ = (std::shared_ptr<Matrix>)(new Matrix(2*nso_,2*nso_));
    Dipole_z_ = (std::shared_ptr<Matrix>)(new Matrix(2*nso_,2*nso_));
    Dipole_x_->zero();
    Dipole_y_->zero();
    Dipole_z_->zero();
    double ** dx = Dipole_x_->pointer();
    double ** dy = Dipole_y_->pointer();
    double ** dz = Dipole_z_->pointer();
    double ** dipole_x_p = dipole_[0]->pointer();
    double ** dipole_y_p = dipole_[1]->pointer();
    double ** dipole_z_p = dipole_[2]->pointer();
    for (size_t mu = 0; mu < nso_; mu++) {
        for (size_t nu = 0; nu < nso_; nu++) {

            dx[mu][nu]           = dipole_x_p[mu][nu];
            dx[mu+nso_][nu+nso_] = dipole_x_p[mu][nu];

            dy[mu][nu]           = dipole_y_p[mu][nu];
            dy[mu+nso_][nu+nso_] = dipole_y_p[mu][nu];

            dz[mu][nu]           = dipole_z_p[mu][nu];
            dz[mu+nso_][nu+nso_] = dipole_z_p[mu][nu];

        }
    }
    Dipole_x_->transform(C_);
    Dipole_y_->transform(C_);
    Dipole_z_->transform(C_);

    // orbital energies
    epsilon_ = (double*)malloc(2*nso_*sizeof(double));
    memset((void*)epsilon_,'\0',2*nso_*sizeof(double));

    for (size_t i = 0; i < 2*nso_; i++) {
        epsilon_[i] = fp[i][i];
    }

    // generate and write SO-basis three-index integrals to disk
    write_three_index_ints();

    // construct four-index integrals
    double val = t1_transformation_molecular_hamiltonian(C_,C_,true);
}

// TODO: this is so wasteful. need to
// 1. fold t1 into 3-index integrals (done!)
// 2. don't build full four-index tensor
// 3. don't build (ac|bd) or any ov^4 quantities

// write three-index integrals to disk (file PSIF_DCC_QSO)
void PolaritonicUCCSD::write_three_index_ints() {

    // get primary basis:
    std::shared_ptr<BasisSet> primary = reference_wavefunction_->get_basisset("ORBITAL");

    // get auxiliary basis:
    std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_CC");

    // total number of auxiliary basis functions
    nQ_ = auxiliary->nbf();

    // three-index integrals
    std::shared_ptr<DFTensor> DF (new DFTensor(primary,auxiliary,Ca_,nalpha_,nso_-nalpha_,nalpha_,nso_-nalpha_,options_));

    std::shared_ptr<Matrix> tmp_so = DF->Qso();
    double ** tmp_so_p = tmp_so->pointer();

    // write Qso to disk
    auto psio = std::make_shared<PSIO>();

    psio->open(PSIF_DCC_QSO, PSIO_OPEN_NEW);
    psio->write_entry(PSIF_DCC_QSO, "Qso CC", (char*)&(tmp_so_p[0][0]), nQ_ * nso_ * nso_ * sizeof(double));
    psio->close(PSIF_DCC_QSO, 1);

}

double PolaritonicUCCSD::t1_transformation() {

    // AO->MO->t1 transformation matrix
    std::shared_ptr<Matrix> CR (new Matrix(C_));
    std::shared_ptr<Matrix> CL (new Matrix(C_));

    double ** CR_p = CR->pointer();
    double ** CL_p = CL->pointer();
    double ** C_p  = C_->pointer();

    size_t n  = 2L*(size_t)nmo_;
    size_t ns = 2L*(size_t)nso_;

#pragma omp parallel for schedule(static)
    for (size_t mu = 0; mu < ns; mu++) {
        for (size_t a = 0; a < v_; a++) {
            double dum = 0.0;
            for (size_t i = 0; i < o_; i++) {
                dum += C_p[mu][i] * t1_[a * o_ + i];
            }
            CL_p[mu][a + o_] -= dum;
        }
    }
#pragma omp parallel for schedule(static)
    for (size_t mu = 0; mu < ns; mu++) {
        for (size_t i = 0; i < o_; i++) {
            double dum = 0.0;
            for (size_t a = 0; a < v_; a++) {
                dum += C_p[mu][a + o_] * t1_[a * o_ + i];
            }
            CR_p[mu][i] += dum;
        }
    }

    double ec = 0.0;

    if ( is_hubbard_ ) {

        ec = t1_transformation_hubbard_hamiltonian(CL,CR,false);

    }else {

        ec = t1_transformation_molecular_hamiltonian(CL,CR,false);

    }
    return ec;

}

double PolaritonicUCCSD::t1_transformation_hubbard_hamiltonian(std::shared_ptr<Matrix> CL, std::shared_ptr<Matrix> CR,bool do_allocate_memory) {

    size_t n  = 2L*(size_t)nmo_;

    // (pr|qs) => (p_sigma p_sigma|p_tau p_tau) 
    double * eri = (double*)malloc(n*n*n*n*sizeof(double));
    memset((void*)eri,'\0',n*n*n*n*sizeof(double));
  
    // add dipole self energy term to eris
    std::shared_ptr<Matrix> dip_a (new Matrix(nso_,nso_));
    if ( nso_ != 4 ) {
        throw PsiException("dipole terms hard-coded for four site Hubbard model",__FILE__,__LINE__);
    }

    double U = options_.get_double("HUBBARD_U");

    for (size_t pa = 0; pa < nso_; pa++) {
        size_t pb = pa + nso_;
        eri[pa*n*n*n+pa*n*n+pb*n+pb] = U;
        eri[pb*n*n*n+pb*n*n+pa*n+pa] = U;
    }

    double ** da_p = dip_a->pointer();
    da_p[0][0] = -1.5;
    da_p[1][1] = -0.5;
    da_p[2][2] =  0.5;
    da_p[3][3] =  1.5;

    double ** dz = Dipole_z_->pointer();
    Dipole_x_->zero();
    Dipole_y_->zero();
    Dipole_z_->zero();
    for (size_t mu = 0; mu < nso_; mu++) {
        for (size_t nu = 0; nu < nso_; nu++) {
            dz[mu][nu]           = da_p[mu][nu];
            dz[mu+nso_][nu+nso_] = da_p[mu][nu];
        }
    }
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

    // dipole self-energy: 1/2 d^2 lambda^2 (added in unpack_eris)
    //for (size_t pa = 0; pa < nso_; pa++) {
    //    for (size_t qa = 0; qa < nso_; qa++) {
    //        size_t pb = pa + nso_;
    //        size_t qb = qa + nso_;
    //        eri[pa*n*n*n+pa*n*n+qa*n+qa] = dz[pa][pa] * dz[qa][qa] * lambda_z * lambda_z;
    //        eri[pb*n*n*n+pb*n*n+qb*n+qb] = dz[pb][pb] * dz[qb][qb] * lambda_z * lambda_z;
    //        eri[pa*n*n*n+pa*n*n+qb*n+qb] = dz[pa][pa] * dz[qb][qb] * lambda_z * lambda_z;
    //        eri[pb*n*n*n+pb*n*n+qa*n+qa] = dz[pb][pb] * dz[qa][qa] * lambda_z * lambda_z;
    //    }
    //}
 
    // transform eris
    double * tmp_eri = (double*)malloc(n*n*n*n*sizeof(double));
    memset((void*)tmp_eri,'\0',n*n*n*n*sizeof(double));

    double ** CR_p = CR->pointer();
    double ** CL_p = CL->pointer();

/*
    for (size_t p = 0; p < n; p++) {
        for (size_t q = 0; q < n; q++) {
            for (size_t r = 0; r < n; r++) {
                for (size_t s = 0; s < n; s++) {
                    double dum = 0.0;
                    for (size_t mu = 0; mu < n; mu++) {
                        for (size_t nu = 0; nu < n; nu++) {
                            for (size_t lambda = 0; lambda < n; lambda++) {
                                for (size_t sigma = 0; sigma < n; sigma++) {
                                    dum += eri[mu*n*n*n+nu*n*n+lambda*n+sigma] 
                                         * CL_p[mu][p] 
                                         * CR_p[nu][q] 
                                         * CL_p[lambda][r] 
                                         * CR_p[sigma][s];
                                }
                            }
                        }
                    }
                    tmp_eri[p*n*n*n+q*n*n+r*n+s] = dum;
                }
            }
        }
    }
    C_DCOPY(n*n*n*n,tmp_eri,1,eri,1);
*/

    // (mu nu | lambda s)
    for (size_t mu = 0; mu < n; mu++) {
        for (size_t nu = 0; nu < n; nu++) {
            for (size_t lambda = 0; lambda < n; lambda++) {
                for (size_t s = 0; s < n; s++) {
                    double dum = 0.0;
                    for (size_t sigma = 0; sigma < n; sigma++) {
                        dum += eri[mu*n*n*n+nu*n*n+lambda*n+sigma] * CR_p[sigma][s];
                    }
                    tmp_eri[mu*n*n*n+nu*n*n+lambda*n+s] = dum;
                }
            }
        }
    }
    // (mu nu | r s)
    for (size_t mu = 0; mu < n; mu++) {
        for (size_t nu = 0; nu < n; nu++) {
            for (size_t r = 0; r < n; r++) {
                for (size_t s = 0; s < n; s++) {
                    double dum = 0.0;
                    for (size_t lambda = 0; lambda < n; lambda++) {
                        dum += tmp_eri[mu*n*n*n+nu*n*n+lambda*n+s] * CL_p[lambda][r];
                    }
                    eri[mu*n*n*n+nu*n*n+r*n+s] = dum;
                }
            }
        }
    }
    // (mu q | r s)
    for (size_t mu = 0; mu < n; mu++) {
        for (size_t q = 0; q < n; q++) {
            for (size_t r = 0; r < n; r++) {
                for (size_t s = 0; s < n; s++) {
                    double dum = 0.0;
                    for (size_t nu = 0; nu < n; nu++) {
                        dum += eri[mu*n*n*n+nu*n*n+r*n+s] * CR_p[nu][q];
                    }
                    tmp_eri[mu*n*n*n+q*n*n+r*n+s] = dum;
                }
            }
        }
    }
    // (p q | r s)
    for (size_t p = 0; p < n; p++) {
        for (size_t q = 0; q < n; q++) {
            for (size_t r = 0; r < n; r++) {
                for (size_t s = 0; s < n; s++) {
                    double dum = 0.0;
                    for (size_t mu = 0; mu < n; mu++) {
                        dum += tmp_eri[mu*n*n*n+q*n*n+r*n+s] * CL_p[mu][p];
                    }
                    eri[p*n*n*n+q*n*n+r*n+s] = dum;
                }
            }
        }
    }
    free(tmp_eri);

    // transform oeis

    // core hamiltonian
    F_->copy(H_);
    double ** hp = H_->pointer();
    double ** fp = F_->pointer();
    for (size_t mu = 0; mu < n; mu++) {
        for (size_t nu = 0; nu < n; nu++) {

            // dipole self-energy: 1/2 d^2 lambda^2
            double val = 0.0;
            for (size_t sigma = 0; sigma < n; sigma++) {
                val += dz[nu][sigma] * dz[sigma][mu];
            }
            fp[mu][nu] += 0.5 * val * lambda_z * lambda_z;
        }
    }
    t1_transform_oei(CL,CR,F_);

    // dipole integrals
    t1_transform_oei(CL,CR,Dipole_z_);

    // unpack different classes of eris
    unpack_eris(eri,do_allocate_memory,false);

    // calculate t1 contribution to correlation energy
    double ec = 0.0;
    for (size_t i = 0; i < o_; i++) {
        ec += 0.5 * fp[i][i];
    }

    // build remainder of fock matrix:

    for (size_t p = 0; p < n; p++) {
        for (size_t q = 0; q < n; q++) {
            double dum = 0.0;
            for (size_t i = 0; i < o_; i++) {
                dum += eri[p*n*n*n+q*n*n+i*n+i] + lambda_z * lambda_z * dz[p][q] * dz[i][i];
                dum -= eri[p*n*n*n+i*n*n+i*n+q] + lambda_z * lambda_z * dz[p][i] * dz[i][q];
            }
            fp[p][q] += dum;
        }
    }

    free(eri);

    for (size_t i = 0; i < o_; i++) {
        ec += 0.5 * fp[i][i];
    }

    // update orbital energies
    for (size_t i = 0; i < n; i++) {
        epsilon_[i] = fp[i][i];
    }

    // now that dipole integrals have been used for self-energy term, scale for coupling terms
    //double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];
    //Dipole_z_->scale(coupling_factor_z);

    // return t1 contribution to correlation energy
    return ec;    
}

double PolaritonicUCCSD::t1_transformation_molecular_hamiltonian(std::shared_ptr<Matrix> CL, std::shared_ptr<Matrix> CR, bool do_allocate_memory) {

    size_t n  = 2L*(size_t)nmo_;
    size_t ns = 2L*(size_t)nso_;

    // read Qso from disk

    double * Qmo = (double*)malloc(nQ_*n*ns*sizeof(double));
    memset((void*)Qmo,'\0',nQ_*n*ns*sizeof(double));

    auto psio = std::make_shared<PSIO>();

    psio->open(PSIF_DCC_QSO, PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_QSO, "Qso CC", (char*)Qmo, nQ_ * nso_ * nso_ * sizeof(double));
    psio->close(PSIF_DCC_QSO, 1);

    memset((void*)tmp1_,'\0',ns*ns*nQ_*sizeof(double));
    for (size_t Q = 0; Q < nQ_; Q++) {
        for (size_t mu = 0; mu < nso_; mu++) {
            for (size_t nu = 0; nu < nso_; nu++) {
                tmp1_[Q*ns*ns+(mu     )*ns+(nu     )] = Qmo[Q*nso_*nso_+mu*nso_+nu];
                tmp1_[Q*ns*ns+(mu+nso_)*ns+(nu+nso_)] = Qmo[Q*nso_*nso_+mu*nso_+nu];
            }
        }
    }

    double ** CR_p = CR->pointer();
    double ** CL_p = CL->pointer();

    // I(Q,mu,p) = C(nu,p) Qso(Q,mu,nu)
    F_DGEMM('n','n',n,ns*nQ_,ns,1.0,&(CL->pointer()[0][0]),n,tmp1_,ns,0.0,Qmo,n);
    for (size_t Q = 0; Q < nQ_; Q++) {
        for (size_t p = 0; p < n; p++) {
            for (size_t mu = 0; mu < ns; mu++) {
                tmp1_[Q*n*ns+p*ns+mu] = Qmo[Q*n*ns+mu*n+p];
            }
        }
    }
    // Qmo(Q,p,q) = C(mu,q) I(Q,p,mu)
    F_DGEMM('n','n',n,n*nQ_,ns,1.0,&(CR->pointer()[0][0]),n,tmp1_,ns,0.0,Qmo,n);

/*
    double * eri = (double*)malloc(n*n*n*n*sizeof(double));
    memset((void*)eri,'\0',n*n*n*n*sizeof(double));

    // (pq|rs) = Qmo(Q,rs) Qmo(Q,pq)
    F_DGEMM('n','t',n*n,n*n,nQ_,1.0,Qmo,n*n,Qmo,n*n,0.0,eri,n*n);
*/

    // transform one-electron integrals

    // core hamiltonian plus one-electron contributions from dipole self energy
    F_->copy(H_);
    F_->add(oe_cavity_terms_);
    t1_transform_oei(CL,CR,F_);

    // dipole integrals
    Dipole_x_->zero();
    Dipole_y_->zero();
    Dipole_z_->zero();
    double ** dx = Dipole_x_->pointer();
    double ** dy = Dipole_y_->pointer();
    double ** dz = Dipole_z_->pointer();
    double ** dipole_x_p = dipole_[0]->pointer();
    double ** dipole_y_p = dipole_[1]->pointer();
    double ** dipole_z_p = dipole_[2]->pointer();
    for (size_t mu = 0; mu < nso_; mu++) {
        for (size_t nu = 0; nu < nso_; nu++) {

            dx[mu][nu]           = dipole_x_p[mu][nu];
            dx[mu+nso_][nu+nso_] = dipole_x_p[mu][nu];

            dy[mu][nu]           = dipole_y_p[mu][nu];
            dy[mu+nso_][nu+nso_] = dipole_y_p[mu][nu];

            dz[mu][nu]           = dipole_z_p[mu][nu];
            dz[mu+nso_][nu+nso_] = dipole_z_p[mu][nu];
        
        }
    }

    t1_transform_oei(CL,CR,Dipole_x_);
    t1_transform_oei(CL,CR,Dipole_y_);
    t1_transform_oei(CL,CR,Dipole_z_);

    double lambda_x = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
    double lambda_y = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

    // add two-electron part of dipole self energy to eris (added in unpack_eris())
    //for (size_t p = 0; p < n; p++) {
    //    for (size_t q = 0; q < n; q++) {
    //        for (size_t r = 0; r < n; r++) {
    //            for (size_t s = 0; s < n; s++) {
    //                eri[p*n*n*n+q*n*n+r*n+s] += lambda_x * lambda_x * dx[p][q] * dx[r][s];
    //                eri[p*n*n*n+q*n*n+r*n+s] += lambda_y * lambda_y * dy[p][q] * dy[r][s];
    //                eri[p*n*n*n+q*n*n+r*n+s] += lambda_z * lambda_z * dz[p][q] * dz[r][s];
    //            }
    //        }
    //    }
    //}

    // unpack different classes of eris. 2e dipole self energy terms will be added here
    //unpack_eris(eri,do_allocate_memory,false);
    //free(eri);
    //unpack_eris(Qmo,do_allocate_memory,true);
    unpack_eris_df(Qmo,do_allocate_memory,true);

    // calculate t1 contribution to correlation energy
    double ** fp = F_->pointer();
    double ec = 0.0;
    for (size_t i = 0; i < o_; i++) {
        ec += 0.5 * fp[i][i];
    }

    // build remainder of fock matrix:

    for (size_t Q = 0; Q < nQ_; Q++) {

        // exchange:
        // sum k (q|rk) (q|ks) 
        F_DGEMM('n','n',n,n,o_,-1.0,Qmo + Q*n*n, n,Qmo + Q*n*n, n,1.0,&(fp[0][0]),n);

        // coulomb
        // sum k (q|kk) (q|rs)
        double dum = 0.0; 
        for (size_t k = 0; k < o_; k++) { 
            dum += Qmo[Q*n*n + k*n + k];
        }
        C_DAXPY(n*n, dum, Qmo + Q*n*n, 1, &(fp[0][0]), 1);
    }

    // dipole self energy contribution to fock matrix
    for (size_t p = 0; p < n; p++) {
        for (size_t q = 0; q < n; q++) {
            double dum = 0.0;
            for (size_t i = 0; i < o_; i++) {
                dum += lambda_x * lambda_x * dx[p][q] * dx[i][i];
                dum -= lambda_x * lambda_x * dx[p][i] * dx[i][q];

                dum += lambda_y * lambda_y * dy[p][q] * dy[i][i];
                dum -= lambda_y * lambda_y * dy[p][i] * dy[i][q];

                dum += lambda_z * lambda_z * dz[p][q] * dz[i][i];
                dum -= lambda_z * lambda_z * dz[p][i] * dz[i][q];
            }
            fp[p][q] += dum;
        }
    }

    // now that dipole integrals have been used for self-energy term, scale for coupling terms
    //double coupling_factor_x = cavity_frequency_[0] * cavity_coupling_strength_[0];
    //double coupling_factor_y = cavity_frequency_[1] * cavity_coupling_strength_[1];
    //double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];
    //Dipole_x_->scale(coupling_factor_x);
    //Dipole_y_->scale(coupling_factor_y);
    //Dipole_z_->scale(coupling_factor_z);

    for (size_t i = 0; i < o_; i++) {
        ec += 0.5 * fp[i][i];
    }
    // adjust ec and recall that reference energy contains both 
    // nuclear repulsion energy and 1/2 [lambda.<de>]^2
    ec += enuc_ - reference_wavefunction_->energy() + average_electric_dipole_self_energy_;

    // update orbital energies
    for (size_t i = 0; i < n; i++) {
        epsilon_[i] = fp[i][i];
    }
    free(Qmo);

    // return t1 contribution to correlation energy
    return ec;    
}

void PolaritonicUCCSD::t1_transform_oei(std::shared_ptr<Matrix> CL, std::shared_ptr<Matrix> CR, std::shared_ptr<Matrix> in) {

    size_t n  = 2L*(size_t)nmo_;
    size_t ns = 2L*(size_t)nso_;

    double ** in_p = in->pointer();

    double ** CL_p = CL->pointer();
    double ** CR_p = CR->pointer();

    double * tmp = (double*)malloc(ns*ns*sizeof(double));
    memset((void*)tmp,'\0',ns*ns*sizeof(double));

    // t1 transformation 
    for (size_t mu = 0; mu < ns; mu++) {
        for (size_t p = 0; p < n; p++) {
            double dum = 0.0;
            for (size_t nu = 0; nu < ns; nu++) {
                dum += CL_p[nu][p] * in_p[nu][mu];
            }
            tmp[p * ns + mu] = dum;
        }
    }
    for (size_t p = 0; p < n; p++) {
        for (size_t q = 0; q < n; q++) {
            double dum = 0.0;
            for (size_t nu = 0; nu < ns; nu++) {
                dum += CR_p[nu][q] * tmp[p * ns + nu];
            }
            in_p[p][q] = dum;
        }
    }
    free(tmp);
}

double PolaritonicUCCSD::compute_eri(size_t p, size_t q, size_t r, size_t s, double * eri, bool is_df) {

    size_t n  = 2L*(size_t)nmo_;

    size_t pq = p * n + q;
    size_t rs = r * n + s;

    double lambda_x = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
    double lambda_y = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);
    
    double ** dx = Dipole_x_->pointer();
    double ** dy = Dipole_y_->pointer();
    double ** dz = Dipole_z_->pointer();

    double dipole_self_energy = 0.0;

    dipole_self_energy += lambda_x * lambda_x * dx[p][q] * dx[r][s];
    dipole_self_energy += lambda_y * lambda_y * dy[p][q] * dy[r][s];
    dipole_self_energy += lambda_z * lambda_z * dz[p][q] * dz[r][s];

    if ( !is_df ) {
        return eri[pq*n*n + rs] + dipole_self_energy;
    }

    return C_DDOT(nQ_,eri + pq ,n*n,eri + rs,n*n) + dipole_self_energy;

}

void PolaritonicUCCSD::unpack_eris_df(double * Qmo, bool do_allocate_memory, bool is_df) {

    size_t n  = 2L*(size_t)nmo_;

    // unpack 3-index integrals

    memset((void*)Qoo_,'\0',nQ_*o_*o_*sizeof(double));
    memset((void*)Qov_,'\0',nQ_*o_*v_*sizeof(double));
    memset((void*)Qvo_,'\0',nQ_*v_*o_*sizeof(double));
    memset((void*)Qvv_,'\0',nQ_*v_*v_*sizeof(double));

#pragma omp parallel for schedule(static)
    for (size_t Q = 0; Q < nQ_; Q++) {
        for (size_t i = 0; i < o_; i++) {
            for (size_t j = 0; j < o_; j++) {
                Qoo_[Q*o_*o_+i*o_+j] = Qmo[Q*n*n+i*n+j];
            }
        }
    }
#pragma omp parallel for schedule(static)
    for (size_t Q = 0; Q < nQ_; Q++) {
        for (size_t i = 0; i < o_; i++) {
            for (size_t a = 0; a < v_; a++) {
                Qov_[Q*o_*v_+i*v_+a] = Qmo[Q*n*n+i*n+(a+o_)];
            }
        }
    }
#pragma omp parallel for schedule(static)
    for (size_t Q = 0; Q < nQ_; Q++) {
        for (size_t a = 0; a < v_; a++) {
            for (size_t i = 0; i < o_; i++) {
                Qvo_[Q*v_*o_+a*o_+i] = Qmo[Q*n*n+(a+o_)*n+i];
            }
        }
    }
#pragma omp parallel for schedule(static)
    for (size_t Q = 0; Q < nQ_; Q++) {
        for (size_t a = 0; a < v_; a++) {
            for (size_t b = 0; b < v_; b++) {
                Qvv_[Q*v_*v_+a*v_+b] = Qmo[Q*n*n+(a+o_)*n+(b+o_)];
            }
        }
    }

    double lambda_x = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
    double lambda_y = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

    double lx2 = lambda_x * lambda_x;
    double ly2 = lambda_y * lambda_y;
    double lz2 = lambda_z * lambda_z;

    double ** dx = Dipole_x_->pointer();
    double ** dy = Dipole_y_->pointer();
    double ** dz = Dipole_z_->pointer();

    // <ij||kl>
    if ( do_allocate_memory ) {
        eri_ijkl_ = (double*)malloc(o_*o_*o_*o_*sizeof(double));
    }
    memset((void*)eri_ijkl_,'\0',o_*o_*o_*o_*sizeof(double));

    F_DGEMM('n','t',o_*o_,o_*o_,nQ_,1.0,Qoo_,o_*o_,Qoo_,o_*o_,0.0,tmp1_,o_*o_);

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t j = 0; j < o_; j++) {
            for (size_t k = 0; k < o_; k++) {
                for (size_t l = 0; l < o_; l++) {
                    //size_t ikjl = i*n*n*n+k*n*n+j*n+l;
                    //size_t iljk = i*n*n*n+l*n*n+j*n+k;
                    //eri_ijkl_[i*o_*o_*o_+j*o_*o_+k*o_+l] = eri[ikjl] - eri[iljk];
                    //eri_ijkl_[i*o_*o_*o_+j*o_*o_+k*o_+l] = compute_eri(i,k,j,l,Qmo,is_df) - compute_eri(i,l,j,k,Qmo,is_df);

                    eri_ijkl_[i*o_*o_*o_+j*o_*o_+k*o_+l] = tmp1_[i*o_*o_*o_+k*o_*o_+j*o_+l] - tmp1_[i*o_*o_*o_+l*o_*o_+j*o_+k];

                    double dipole_self_energy = 0.0;

                    dipole_self_energy += lx2 * dx[i][k] * dx[j][l];
                    dipole_self_energy += ly2 * dy[i][k] * dy[j][l];
                    dipole_self_energy += lz2 * dz[i][k] * dz[j][l];

                    dipole_self_energy -= lx2 * dx[i][l] * dx[j][k];
                    dipole_self_energy -= ly2 * dy[i][l] * dy[j][k];
                    dipole_self_energy -= lz2 * dz[i][l] * dz[j][k];

                    eri_ijkl_[i*o_*o_*o_+j*o_*o_+k*o_+l] += dipole_self_energy;

                }
            }
        }
    }

/*
    // <ab||cd>
    if ( do_allocate_memory ) {
        eri_abcd_ = (double*)malloc(v_*v_*v_*v_*sizeof(double));
    }
    memset((void*)eri_abcd_,'\0',v_*v_*v_*v_*sizeof(double));

#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t c = 0; c < v_; c++) {
                for (size_t d = 0; d < v_; d++) {
                    //size_t acbd = (a+o_)*n*n*n+(c+o_)*n*n+(b+o_)*n+(d+o_);
                    //size_t adbc = (a+o_)*n*n*n+(d+o_)*n*n+(b+o_)*n+(c+o_);
                    //eri_abcd_[a*v_*v_*v_+b*v_*v_+c*v_+d] = eri[acbd] - eri[adbc];
                    //eri_abcd_[a*v_*v_*v_+b*v_*v_+c*v_+d] = compute_eri(a+o_,c+o_,b+o_,d+o_,Qmo,is_df) - compute_eri(a+o_,d+o_,b+o_,c+o_,Qmo,is_df);
                    eri_abcd_[a*v_*v_*v_+b*v_*v_+c*v_+d] = compute_eri(a+o_,c+o_,b+o_,d+o_,Qmo,is_df);// - compute_eri(a+o_,d+o_,b+o_,c+o_,Qmo,is_df);
                }
            }
        }
    }
*/

    // <ij||ab>
    if ( do_allocate_memory ) {
        eri_ijab_ = (double*)malloc(o_*o_*v_*v_*sizeof(double));
    }
    memset((void*)eri_ijab_,'\0',o_*o_*v_*v_*sizeof(double));

    F_DGEMM('n','t',o_*v_,o_*v_,nQ_,1.0,Qov_,o_*v_,Qov_,o_*v_,0.0,tmp1_,o_*v_);

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t j = 0; j < o_; j++) {
            for (size_t a = 0; a < v_; a++) {
                for (size_t b = 0; b < v_; b++) {
                    //size_t iajb = i*n*n*n+(a+o_)*n*n+j*n+(b+o_);
                    //size_t ibja = i*n*n*n+(b+o_)*n*n+j*n+(a+o_);
                    //eri_ijab_[i*o_*v_*v_+j*v_*v_+a*v_+b] = eri[iajb] - eri[ibja];
                    //eri_ijab_[i*o_*v_*v_+j*v_*v_+a*v_+b] = compute_eri(i,a+o_,j,b+o_,Qmo,is_df) - compute_eri(i,b+o_,j,a+o_,Qmo,is_df);

                    eri_ijab_[i*o_*v_*v_+j*v_*v_+a*v_+b] = tmp1_[i*o_*v_*v_+a*o_*v_+j*v_+b] - tmp1_[i*o_*v_*v_+b*o_*v_+j*v_+a];

                    double dipole_self_energy = 0.0;

                    dipole_self_energy += lx2 * dx[i][a+o_] * dx[j][b+o_];
                    dipole_self_energy += ly2 * dy[i][a+o_] * dy[j][b+o_];
                    dipole_self_energy += lz2 * dz[i][a+o_] * dz[j][b+o_];

                    dipole_self_energy -= lx2 * dx[i][b+o_] * dx[j][a+o_];
                    dipole_self_energy -= ly2 * dy[i][b+o_] * dy[j][a+o_];
                    dipole_self_energy -= lz2 * dz[i][b+o_] * dz[j][a+o_];

                    eri_ijab_[i*o_*v_*v_+j*v_*v_+a*v_+b] += dipole_self_energy;

                }
            }
        }
    }

/*
    // <ab||ij>
    if ( do_allocate_memory ) {
        eri_abij_ = (double*)malloc(o_*o_*v_*v_*sizeof(double));
    }
    memset((void*)eri_abij_,'\0',o_*o_*v_*v_*sizeof(double));

    F_DGEMM('n','t',o_*v_,o_*v_,nQ_,1.0,Qvo_,o_*v_,Qvo_,o_*v_,0.0,tmp1_,o_*v_);

#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    //size_t aibj = (a+o_)*n*n*n+i*n*n+(b+o_)*n+j;
                    //size_t ajbi = (a+o_)*n*n*n+j*n*n+(b+o_)*n+i;
                    //eri_abij_[a*o_*o_*v_+b*o_*o_+i*o_+j] = eri[aibj] - eri[ajbi];
                    //eri_abij_[a*o_*o_*v_+b*o_*o_+i*o_+j] = compute_eri(a+o_,i,b+o_,j,Qmo,is_df) - compute_eri(a+o_,j,b+o_,i,Qmo,is_df);

                    eri_abij_[a*o_*o_*v_+b*o_*o_+i*o_+j] = tmp1_[a*o_*o_*v_+i*o_*v_+b*o_+j] - tmp1_[a*o_*o_*v_+j*o_*v_+b*o_+i];

                    double dipole_self_energy = 0.0;

                    dipole_self_energy += lx2 * dx[a+o_][i] * dx[b+o_][j];
                    dipole_self_energy += ly2 * dy[a+o_][i] * dy[b+o_][j];
                    dipole_self_energy += lz2 * dz[a+o_][i] * dz[b+o_][j];

                    dipole_self_energy -= lx2 * dx[a+o_][j] * dx[b+o_][i];
                    dipole_self_energy -= ly2 * dy[a+o_][j] * dy[b+o_][i];
                    dipole_self_energy -= lz2 * dz[a+o_][j] * dz[b+o_][i];

                    eri_abij_[a*o_*o_*v_+b*o_*o_+i*o_+j] += dipole_self_energy;
                }
            }
        }
    }
*/

    // <ia||jb> = (ij|ab) - (ib|aj)
    if ( do_allocate_memory ) {
        eri_iajb_ = (double*)malloc(o_*o_*v_*v_*sizeof(double));
    }
    memset((void*)eri_iajb_,'\0',o_*o_*v_*v_*sizeof(double));

    // (ij|ab)
    F_DGEMM('n','t',v_*v_,o_*o_,nQ_,1.0,Qvv_,v_*v_,Qoo_,o_*o_,0.0,tmp1_,v_*v_);
    // (ib|aj)
    F_DGEMM('n','t',o_*v_,o_*v_,nQ_,1.0,Qvo_,o_*v_,Qov_,o_*v_,0.0,tmp2_,o_*v_);

#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    //size_t ijab = i*n*n*n+j*n*n+(a+o_)*n+(b+o_);
                    //size_t ibaj = i*n*n*n+(b+o_)*n*n+(a+o_)*n+j;
                    //eri_iajb_[i*o_*v_*v_+a*o_*v_+j*v_+b] = eri[ijab] - eri[ibaj];
                    //eri_iajb_[i*o_*v_*v_+a*o_*v_+j*v_+b] = compute_eri(i,j,a+o_,b+o_,Qmo,is_df) - compute_eri(i,b+o_,a+o_,j,Qmo,is_df);

                    eri_iajb_[i*o_*v_*v_+a*o_*v_+j*v_+b] = tmp1_[i*v_*v_*o_+j*v_*v_+a*v_+b] - tmp2_[i*o_*v_*v_+b*o_*v_+a*o_+j];

                    double dipole_self_energy = 0.0;

                    dipole_self_energy += lx2 * dx[i][j] * dx[a+o_][b+o_];
                    dipole_self_energy += ly2 * dy[i][j] * dy[a+o_][b+o_];
                    dipole_self_energy += lz2 * dz[i][j] * dz[a+o_][b+o_];

                    dipole_self_energy -= lx2 * dx[i][b+o_] * dx[a+o_][j];
                    dipole_self_energy -= ly2 * dy[i][b+o_] * dy[a+o_][j];
                    dipole_self_energy -= lz2 * dz[i][b+o_] * dz[a+o_][j];

                    eri_iajb_[i*o_*v_*v_+a*o_*v_+j*v_+b] += dipole_self_energy;
                }
            }
        }
    }

    // <ia||jk> = (ij|ak) - (ik|aj)
    if ( do_allocate_memory ) {
        eri_iajk_ = (double*)malloc(o_*o_*o_*v_*sizeof(double));
    }
    memset((void*)eri_iajk_,'\0',o_*o_*o_*v_*sizeof(double));

    F_DGEMM('n','t',o_*v_,o_*o_,nQ_,1.0,Qvo_,o_*v_,Qoo_,o_*o_,0.0,tmp1_,o_*v_);

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t a = 0; a < v_; a++) {
            for (size_t j = 0; j < o_; j++) {
                for (size_t k = 0; k < o_; k++) {
                    //size_t ijak = i*n*n*n+j*n*n+(a+o_)*n+k;
                    //size_t ikaj = i*n*n*n+k*n*n+(a+o_)*n+j;
                    //eri_iajk_[i*o_*o_*v_+a*o_*o_+j*o_+k] = eri[ijak] - eri[ikaj];
                    //eri_iajk_[i*o_*o_*v_+a*o_*o_+j*o_+k] = compute_eri(i,j,a+o_,k,Qmo,is_df) - compute_eri(i,k,a+o_,j,Qmo,is_df);

                    eri_iajk_[i*o_*o_*v_+a*o_*o_+j*o_+k] = tmp1_[i*o_*o_*v_+j*o_*v_+a*o_+k] - tmp1_[i*o_*o_*v_+k*o_*v_+a*o_+j];

                    double dipole_self_energy = 0.0;

                    dipole_self_energy += lx2 * dx[i][j] * dx[a+o_][k];
                    dipole_self_energy += ly2 * dy[i][j] * dy[a+o_][k];
                    dipole_self_energy += lz2 * dz[i][j] * dz[a+o_][k];

                    dipole_self_energy -= lx2 * dx[i][k] * dx[a+o_][j];
                    dipole_self_energy -= ly2 * dy[i][k] * dy[a+o_][j];
                    dipole_self_energy -= lz2 * dz[i][k] * dz[a+o_][j];

                    eri_iajk_[i*o_*o_*v_+a*o_*o_+j*o_+k] += dipole_self_energy;
                }
            }
        }
    }

    // <jk||ia> = (ji|ka) - (ja|ki)
    if ( do_allocate_memory ) {
        eri_jkia_ = (double*)malloc(o_*o_*o_*v_*sizeof(double));
    }
    memset((void*)eri_jkia_,'\0',o_*o_*o_*v_*sizeof(double));

    F_DGEMM('n','t',o_*v_,o_*o_,nQ_,1.0,Qov_,o_*v_,Qoo_,o_*o_,0.0,tmp1_,o_*v_);

#pragma omp parallel for schedule(static)
    for (size_t j = 0; j < o_; j++) {
        for (size_t k = 0; k < o_; k++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t a = 0; a < v_; a++) {
                    //size_t jika = j*n*n*n+i*n*n+k*n+(a+o_);
                    //size_t jaki = j*n*n*n+(a+o_)*n*n+k*n+i;
                    //eri_jkia_[j*o_*o_*v_+k*o_*v_+i*v_+a] = eri[jika] - eri[jaki];
                    //eri_jkia_[j*o_*o_*v_+k*o_*v_+i*v_+a] = compute_eri(j,i,k,a+o_,Qmo,is_df) - compute_eri(j,a+o_,k,i,Qmo,is_df);

                    eri_jkia_[j*o_*o_*v_+k*o_*v_+i*v_+a] = tmp1_[j*o_*o_*v_+i*o_*v_+k*v_+a] - tmp1_[k*o_*o_*v_+i*o_*v_+j*v_+a];

                    double dipole_self_energy = 0.0;

                    dipole_self_energy += lx2 * dx[j][i] * dx[k][a+o_];
                    dipole_self_energy += ly2 * dy[j][i] * dy[k][a+o_];
                    dipole_self_energy += lz2 * dz[j][i] * dz[k][a+o_];

                    dipole_self_energy -= lx2 * dx[j][a+o_] * dx[k][i];
                    dipole_self_energy -= ly2 * dy[j][a+o_] * dy[k][i];
                    dipole_self_energy -= lz2 * dz[j][a+o_] * dz[k][i];

                    eri_jkia_[j*o_*o_*v_+k*o_*v_+i*v_+a] += dipole_self_energy;
                }
            }
        }
    }

/*
    // <ai||bc> = (ab|ic) - (ac|ib)
    if ( do_allocate_memory ) {
        eri_aibc_ = (double*)malloc(o_*v_*v_*v_*sizeof(double));
    }
    memset((void*)eri_aibc_,'\0',o_*v_*v_*v_*sizeof(double));

    F_DGEMM('n','t',o_*v_,v_*v_,nQ_,1.0,Qov_,o_*v_,Qvv_,v_*v_,0.0,tmp1_,o_*v_);

#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t c = 0; c < v_; c++) {
                for (size_t i = 0; i < o_; i++) {
                    //size_t abic = (a+o_)*n*n*n+(b+o_)*n*n+i*n+(c+o_);
                    //size_t acib = (a+o_)*n*n*n+(c+o_)*n*n+i*n+(b+o_);
                    //eri_aibc_[a*o_*v_*v_+i*v_*v_+b*v_+c] = eri[abic] - eri[acib];
                    //eri_aibc_[a*o_*v_*v_+i*v_*v_+b*v_+c] = compute_eri(a+o_,b+o_,i,c+o_,Qmo,is_df) - compute_eri(a+o_,c+o_,i,b+o_,Qmo,is_df);

                    eri_aibc_[a*o_*v_*v_+i*v_*v_+b*v_+c] = tmp1_[a*o_*v_*v_+b*o_*v_+i*v_+c] - tmp1_[a*o_*v_*v_+c*o_*v_+i*v_+b];

                    double dipole_self_energy = 0.0;

                    dipole_self_energy += lx2 * dx[a+o_][b+o_] * dx[i][c+o_];
                    dipole_self_energy += ly2 * dy[a+o_][b+o_] * dy[i][c+o_];
                    dipole_self_energy += lz2 * dz[a+o_][b+o_] * dz[i][c+o_];

                    dipole_self_energy -= lx2 * dx[a+o_][c+o_] * dx[i][b+o_];
                    dipole_self_energy -= ly2 * dy[a+o_][c+o_] * dy[i][b+o_];
                    dipole_self_energy -= lz2 * dz[a+o_][c+o_] * dz[i][b+o_];

                    eri_aibc_[a*o_*v_*v_+i*v_*v_+b*v_+c] += dipole_self_energy;
                }
            }
        }
    }
*/

/*
    // <ab||ic> = (ai|bc) - (ac|bi)
    if ( do_allocate_memory ) {
        eri_abic_ = (double*)malloc(o_*v_*v_*v_*sizeof(double));
    }
    memset((void*)eri_abic_,'\0',o_*v_*v_*v_*sizeof(double));

    F_DGEMM('n','t',v_*v_,o_*v_,nQ_,1.0,Qvv_,v_*v_,Qvo_,o_*v_,0.0,tmp1_,v_*v_);

#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t c = 0; c < v_; c++) {
                for (size_t i = 0; i < o_; i++) {
                    //size_t aibc = (a+o_)*n*n*n+i*n*n+(b+o_)*n+(c+o_);
                    //size_t acbi = (a+o_)*n*n*n+(c+o_)*n*n+(b+o_)*n+i;
                    //eri_abic_[a*o_*v_*v_+b*o_*v_+c*o_+i] = eri[aibc] - eri[acbi];
                    //eri_abic_[a*o_*v_*v_+b*o_*v_+c*o_+i] = compute_eri(a+o_,i,b+o_,c+o_,Qmo,is_df) - compute_eri(a+o_,c+o_,b+o_,i,Qmo,is_df);

                    eri_abic_[a*o_*v_*v_+b*o_*v_+i*v_+c] = tmp1_[a*o_*v_*v_+i*v_*v_+b*v_+c] - tmp1_[b*o_*v_*v_+i*v_*v_+a*v_+c];

                    double dipole_self_energy = 0.0;

                    dipole_self_energy += lx2 * dx[a+o_][i] * dx[b+o_][c+o_];
                    dipole_self_energy += ly2 * dy[a+o_][i] * dy[b+o_][c+o_];
                    dipole_self_energy += lz2 * dz[a+o_][i] * dz[b+o_][c+o_];

                    dipole_self_energy -= lx2 * dx[a+o_][c+o_] * dx[b+o_][i];
                    dipole_self_energy -= ly2 * dy[a+o_][c+o_] * dy[b+o_][i];
                    dipole_self_energy -= lz2 * dz[a+o_][c+o_] * dz[b+o_][i];

                    eri_abic_[a*o_*v_*v_+b*o_*v_+i*v_+c] += dipole_self_energy;
                }
            }
        }
    }
*/

    // write Qvv and its transpose to disk for double particle ladder diagram
    auto psio = std::make_shared<PSIO>();
    psio->open(PSIF_DCC_QSO, PSIO_OPEN_OLD);
    psio->write_entry(PSIF_DCC_QSO, "Qmo", (char *)&Qmo[0], nQ_ * n * n * sizeof(double));
    psio->write_entry(PSIF_DCC_QSO, "Qvv", (char *)&Qvv_[0], nQ_ * v_ * v_ * sizeof(double));
#pragma omp parallel for schedule(static)
    for (size_t Q = 0; Q < nQ_; Q++) {
        C_DCOPY(v_ * v_, Qvv_ + Q * v_ * v_, 1, Qmo + Q, nQ_);
    }
    psio->write_entry(PSIF_DCC_QSO, "Qvv transpose", (char *)&Qmo[0], nQ_ * v_ * v_ * sizeof(double));

    // be sure to read Qmo back in because it will be used after this function
    psio->read_entry(PSIF_DCC_QSO, "Qmo", (char *)&Qmo[0], nQ_ * n * n * sizeof(double));
    psio->close(PSIF_DCC_QSO, 1);

}

void PolaritonicUCCSD::unpack_eris(double * eri, bool do_allocate_memory, bool is_df) {

    size_t n  = 2L*(size_t)nmo_;

    // <ij||kl>
    if ( do_allocate_memory ) {
        eri_ijkl_ = (double*)malloc(o_*o_*o_*o_*sizeof(double));
    }
    memset((void*)eri_ijkl_,'\0',o_*o_*o_*o_*sizeof(double));

    for (size_t i = 0; i < o_; i++) {
        for (size_t j = 0; j < o_; j++) {
            for (size_t k = 0; k < o_; k++) {
                for (size_t l = 0; l < o_; l++) {
                    //size_t ikjl = i*n*n*n+k*n*n+j*n+l;
                    //size_t iljk = i*n*n*n+l*n*n+j*n+k;
                    //eri_ijkl_[i*o_*o_*o_+j*o_*o_+k*o_+l] = eri[ikjl] - eri[iljk];
                    eri_ijkl_[i*o_*o_*o_+j*o_*o_+k*o_+l] = compute_eri(i,k,j,l,eri,is_df) - compute_eri(i,l,j,k,eri,is_df);
                }
            }
        }
    }

    // <ab||cd>
    if ( do_allocate_memory ) {
        eri_abcd_ = (double*)malloc(v_*v_*v_*v_*sizeof(double));
    }
    memset((void*)eri_abcd_,'\0',v_*v_*v_*v_*sizeof(double));

    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t c = 0; c < v_; c++) {
                for (size_t d = 0; d < v_; d++) {
                    //size_t acbd = (a+o_)*n*n*n+(c+o_)*n*n+(b+o_)*n+(d+o_);
                    //size_t adbc = (a+o_)*n*n*n+(d+o_)*n*n+(b+o_)*n+(c+o_);
                    //eri_abcd_[a*v_*v_*v_+b*v_*v_+c*v_+d] = eri[acbd] - eri[adbc];
                    //eri_abcd_[a*v_*v_*v_+b*v_*v_+c*v_+d] = compute_eri(a+o_,c+o_,b+o_,d+o_,eri,is_df) - compute_eri(a+o_,d+o_,b+o_,c+o_,eri,is_df);
                    eri_abcd_[a*v_*v_*v_+b*v_*v_+c*v_+d] = compute_eri(a+o_,c+o_,b+o_,d+o_,eri,is_df);// - compute_eri(a+o_,d+o_,b+o_,c+o_,eri,is_df);
                }
            }
        }
    }

    // <ij||ab>
    if ( do_allocate_memory ) {
        eri_ijab_ = (double*)malloc(o_*o_*v_*v_*sizeof(double));
    }
    memset((void*)eri_ijab_,'\0',o_*o_*v_*v_*sizeof(double));

    for (size_t i = 0; i < o_; i++) {
        for (size_t j = 0; j < o_; j++) {
            for (size_t a = 0; a < v_; a++) {
                for (size_t b = 0; b < v_; b++) {
                    //size_t iajb = i*n*n*n+(a+o_)*n*n+j*n+(b+o_);
                    //size_t ibja = i*n*n*n+(b+o_)*n*n+j*n+(a+o_);
                    //eri_ijab_[i*o_*v_*v_+j*v_*v_+a*v_+b] = eri[iajb] - eri[ibja];
                    eri_ijab_[i*o_*v_*v_+j*v_*v_+a*v_+b] = compute_eri(i,a+o_,j,b+o_,eri,is_df) - compute_eri(i,b+o_,j,a+o_,eri,is_df);
                }
            }
        }
    }

    // <ab||ij>
    if ( do_allocate_memory ) {
        eri_abij_ = (double*)malloc(o_*o_*v_*v_*sizeof(double));
    }
    memset((void*)eri_abij_,'\0',o_*o_*v_*v_*sizeof(double));

    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    //size_t aibj = (a+o_)*n*n*n+i*n*n+(b+o_)*n+j;
                    //size_t ajbi = (a+o_)*n*n*n+j*n*n+(b+o_)*n+i;
                    //eri_abij_[a*o_*o_*v_+b*o_*o_+i*o_+j] = eri[aibj] - eri[ajbi];
                    eri_abij_[a*o_*o_*v_+b*o_*o_+i*o_+j] = compute_eri(a+o_,i,b+o_,j,eri,is_df) - compute_eri(a+o_,j,b+o_,i,eri,is_df);
                }
            }
        }
    }

    // <ia||jb>
    if ( do_allocate_memory ) {
        eri_iajb_ = (double*)malloc(o_*o_*v_*v_*sizeof(double));
    }
    memset((void*)eri_iajb_,'\0',o_*o_*v_*v_*sizeof(double));
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    //size_t ijab = i*n*n*n+j*n*n+(a+o_)*n+(b+o_);
                    //size_t ibaj = i*n*n*n+(b+o_)*n*n+(a+o_)*n+j;
                    //eri_iajb_[i*o_*v_*v_+a*o_*v_+j*v_+b] = eri[ijab] - eri[ibaj];
                    eri_iajb_[i*o_*v_*v_+a*o_*v_+j*v_+b] = compute_eri(i,j,a+o_,b+o_,eri,is_df) - compute_eri(i,b+o_,a+o_,j,eri,is_df);
                }
            }
        }
    }

    // <ia||jk>
    if ( do_allocate_memory ) {
        eri_iajk_ = (double*)malloc(o_*o_*o_*v_*sizeof(double));
    }
    memset((void*)eri_iajk_,'\0',o_*o_*o_*v_*sizeof(double));
    for (size_t i = 0; i < o_; i++) {
        for (size_t a = 0; a < v_; a++) {
            for (size_t j = 0; j < o_; j++) {
                for (size_t k = 0; k < o_; k++) {
                    //size_t ijak = i*n*n*n+j*n*n+(a+o_)*n+k;
                    //size_t ikaj = i*n*n*n+k*n*n+(a+o_)*n+j;
                    //eri_iajk_[i*o_*o_*v_+a*o_*o_+j*o_+k] = eri[ijak] - eri[ikaj];
                    eri_iajk_[i*o_*o_*v_+a*o_*o_+j*o_+k] = compute_eri(i,j,a+o_,k,eri,is_df) - compute_eri(i,k,a+o_,j,eri,is_df);
                }
            }
        }
    }

    // <jk||ia>
    if ( do_allocate_memory ) {
        eri_jkia_ = (double*)malloc(o_*o_*o_*v_*sizeof(double));
    }
    memset((void*)eri_jkia_,'\0',o_*o_*o_*v_*sizeof(double));
    for (size_t j = 0; j < o_; j++) {
        for (size_t k = 0; k < o_; k++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t a = 0; a < v_; a++) {
                    //size_t jika = j*n*n*n+i*n*n+k*n+(a+o_);
                    //size_t jaki = j*n*n*n+(a+o_)*n*n+k*n+i;
                    //eri_jkia_[j*o_*o_*v_+k*o_*v_+i*v_+a] = eri[jika] - eri[jaki];
                    eri_jkia_[j*o_*o_*v_+k*o_*v_+i*v_+a] = compute_eri(j,i,k,a+o_,eri,is_df) - compute_eri(j,a+o_,k,i,eri,is_df);
                }
            }
        }
    }

    // <ai||bc>
    if ( do_allocate_memory ) {
        eri_aibc_ = (double*)malloc(o_*v_*v_*v_*sizeof(double));
    }
    memset((void*)eri_aibc_,'\0',o_*v_*v_*v_*sizeof(double));

    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t c = 0; c < v_; c++) {
                for (size_t i = 0; i < o_; i++) {
                    //size_t abic = (a+o_)*n*n*n+(b+o_)*n*n+i*n+(c+o_);
                    //size_t acib = (a+o_)*n*n*n+(c+o_)*n*n+i*n+(b+o_);
                    //eri_aibc_[a*o_*v_*v_+i*v_*v_+b*v_+c] = eri[abic] - eri[acib];
                    eri_aibc_[a*o_*v_*v_+i*v_*v_+b*v_+c] = compute_eri(a+o_,b+o_,i,c+o_,eri,is_df) - compute_eri(a+o_,c+o_,i,b+o_,eri,is_df);
                }
            }
        }
    }

    // <ab||ic> = (ai|bc) - (ac|bi)
    if ( do_allocate_memory ) {
        eri_abic_ = (double*)malloc(o_*v_*v_*v_*sizeof(double));
    }
    memset((void*)eri_abic_,'\0',o_*v_*v_*v_*sizeof(double));
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t c = 0; c < v_; c++) {
                for (size_t i = 0; i < o_; i++) {
                    //size_t aibc = (a+o_)*n*n*n+i*n*n+(b+o_)*n+(c+o_);
                    //size_t acbi = (a+o_)*n*n*n+(c+o_)*n*n+(b+o_)*n+i;
                    //eri_abic_[a*o_*v_*v_+b*o_*v_+c*o_+i] = eri[aibc] - eri[acbi];
                    eri_abic_[a*o_*v_*v_+b*o_*v_+i*v_+c] = compute_eri(a+o_,i,b+o_,c+o_,eri,is_df) - compute_eri(a+o_,c+o_,b+o_,i,eri,is_df);
                }
            }
        }
    }

}

double PolaritonicUCCSD::compute_energy() {

    // grab some input options_
    double e_convergence = options_.get_double("E_CONVERGENCE");
    double r_convergence = options_.get_double("R_CONVERGENCE");
    size_t maxiter          = options_.get_int("MAXITER");

    outfile->Printf("\n");
    outfile->Printf("    No. basis functions:            %5i\n",nso_);
    outfile->Printf("    No. auxiliary basis functions:  %5i\n",nQ_);
    outfile->Printf("    No. alpha electrons:            %5i\n",nalpha_);
    outfile->Printf("    No. beta electrons:             %5i\n",nbeta_);
    outfile->Printf("    e_convergence:             %10.3le\n",e_convergence);
    outfile->Printf("    r_convergence:             %10.3le\n",r_convergence);
    outfile->Printf("    maxiter:                        %5i\n",maxiter);
    outfile->Printf("\n");
    outfile->Printf("\n");

    // CCSD iterations without photon
    double ec = cc_iterations();
/*
    include_u0_ = false;
    include_u1_ = false;
    include_u2_ = false;

    ccamps_dim_ = o_*o_*v_*v_ + o_*v_;

    double ec = cc_iterations();

    // CCSD iterations with photon
    include_u0_ = options_.get_bool("POLARITONIC_CC_INCLUDE_U0");;
    include_u1_ = options_.get_bool("POLARITONIC_CC_INCLUDE_U1");;
    include_u2_ = options_.get_bool("POLARITONIC_CC_INCLUDE_U2");;

    if ( include_u0_ || include_u1_ || include_u2_ ) {

        if ( include_u0_ ) {
            ccamps_dim_++;
        }
        if ( include_u1_ ) {
            ccamps_dim_ += o_*v_;
        }
        if ( include_u2_ ) {
            ccamps_dim_ += o_*o_*v_*v_;
        }
        diis->restart();

        ec = cc_iterations();

    }
*/

    outfile->Printf("    * Polaritonic UCCSD total energy: %20.12lf\n",energy_ + ec);

    Process::environment.globals["UCCSD TOTAL ENERGY"] = energy_ + ec;
    Process::environment.globals["CURRENT ENERGY"] = energy_ + ec;

    return energy_;

}

double PolaritonicUCCSD::cc_iterations() {

    // grab some input options_
    double e_convergence = options_.get_double("E_CONVERGENCE");
    double r_convergence = options_.get_double("R_CONVERGENCE");
    size_t maxiter          = options_.get_int("MAXITER");

    double e_last  = 0.0;
    double dele    = 0.0;
    double tnorm   = 0.0;

    outfile->Printf("\n");
    outfile->Printf("    ==>  Begin CCSD Iterations <==\n");
    outfile->Printf("\n");
    outfile->Printf("    ");
    outfile->Printf(" Iter ");
    outfile->Printf("              energy ");
    outfile->Printf("                  dE ");
    outfile->Printf("                |dT| ");
    outfile->Printf("\n");

    double ec = 0.0;

    size_t iter = 0;
    do {

        e_last = energy_ + ec;

        // build residual 
        residual();

        // update amplitudes 
        tnorm = update_amplitudes();

        // t1-transformation e(-T1) H e(T1)
        ec = t1_transformation();

        // evaluate t1 contribution to correlation energy
        ec += correlation_energy();

        // dele
        dele = energy_ + ec - e_last;

        outfile->Printf("    %5i %20.12lf %20.12lf %20.12lf\n",iter,ec,dele,tnorm);

        iter++;
        if ( iter > maxiter ) break;

    }while(fabs(dele) > e_convergence || tnorm > r_convergence );

    if ( iter > maxiter ) {
        throw PsiException("Maximum number of iterations exceeded!",__FILE__,__LINE__);
    }

    outfile->Printf("\n");
    outfile->Printf("    CCSD iterations converged!\n");
    outfile->Printf("\n");

    return ec;

}

// evaluate total residual
void PolaritonicUCCSD::residual() {

    memset((void*)residual_,'\0',ccamps_dim_*sizeof(double));

    residual_t1();

    residual_t2();

    if ( include_u0_ ) {
        residual_u0();
    }

    if ( include_u1_ ) {
        residual_u1();
    }

    if ( include_u2_ ) {
        residual_u2();
    }

}

// 
// u1 residual
// 
// - 1.00000 d+(e,m) ...done
// + 1.00000 u1(e,m) w0 ...done

// - 1.00000 F(i,m) u1(e,i) ...done
// + 1.00000 F(e,a) u1(a,m) ...done

// + 1.00000 <i,e||a,m> u1(a,i) ...done

// - 1.00000 d+(i,a) t2(a,e,i,m) ... done
// + 1.00000 <i,j||a,b> t2(b,e,j,m) u1(a,i)
// + 0.50000 <i,j||a,b> t2(b,e,i,j) u1(a,m)
// + 0.50000 <i,j||a,b> t2(a,b,j,m) u1(e,i)

// + 1.00000 d-(i,m) u1(e,i) u0 ...done
// - 1.00000 d-(e,a) u1(a,m) u0 ...done

// - 1.00000 d-(i,a) u1(a,i) u1(e,m) ...done
// + 2.00000 d-(i,a) u1(a,m) u1(e,i) ...done

// + 1.00000 F(i,a) u2(a,e,i,m)  ...done
// - 1.00000 d-(i,a) u2(a,e,i,m) u0 ...done
// - 0.50000 <i,j||a,m> u2(a,e,i,j) ...done
// + 0.50000 <i,e||a,b> u2(a,b,i,m) ...done

void PolaritonicUCCSD::residual_u1() {

    memset((void*)ru1_,'\0',o_*v_*sizeof(double));

    double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];
    std::shared_ptr<Matrix> dip (new Matrix(Dipole_z_));
    dip->scale(coupling_factor_z);
    double ** dp = dip->pointer();

    double ** fp = F_->pointer();

    // - 1.00000 d+(e,m)
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t m = 0; m < o_; m++) {
            ru1_[e*o_+m] = -dp[e+o_][m];
        }
    }
    // + 1.00000 u1(e,m) w0
    // TODO: generalize for x,y,z
    double w0 = cavity_frequency_[2];
    C_DAXPY(o_*v_,w0,u1_,1,ru1_,1);

    // - 1.00000 F(i,m) u1(e,i)
    for (size_t i = 0; i < o_; i++) {
        for (size_t m = 0; m < o_; m++) {
            tmp1_[i*o_+m] = fp[i][m];
        }
    }
    F_DGEMM('n','n',o_,v_,o_,-1.0,tmp1_,o_,u1_,o_,1.0,ru1_,o_);

    // + 1.00000 u1(a,m) F(e,a) 
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t a = 0; a < v_; a++) {
            tmp1_[e*v_+a] = fp[e+o_][a+o_];
        }
    }
    F_DGEMM('n','n',o_,v_,v_,1.0,u1_,o_,tmp1_,v_,1.0,ru1_,o_);

  
    // - 1.00000 <i,e||m,a> u1(a,i)
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t m = 0; m < o_; m++) {
            double dum = 0.0;
            for (size_t a = 0; a < v_; a++) {
                for (size_t i = 0; i < o_; i++) {
                    dum += eri_iajb_[i*o_*v_*v_+e*o_*v_+m*v_+a] * u1_[a*o_+i];
                }
            }
            ru1_[e*o_+m] -= dum;
        }
    }

    // - 1.00000 d+(i,a) t2(a,e,i,m)
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t m = 0; m < o_; m++) {
            double dum = 0.0;
            for (size_t a = 0; a < v_; a++) {
                for (size_t i = 0; i < o_; i++) {
                    dum += t2_[a*o_*o_*v_+e*o_*o_+i*o_+m] * dp[i][a+o_];
                }
            }
            ru1_[e*o_+m] -= dum;
        }
    }


    // + 1.00000 <i,j||a,b> t2(b,e,j,m) u1(a,i)
#pragma omp parallel for schedule(static)
    for (size_t j = 0; j < o_; j++) {
        for (size_t b = 0; b < v_; b++) {
            double dum = 0.0;
            for (size_t a = 0; a < v_; a++) {
                for (size_t i = 0; i < o_; i++) {
                    dum += eri_ijab_[i*o_*v_*v_+j*v_*v_+a*v_+b] * u1_[a*o_+i];
                }
            }
            tmp1_[j*v_+b] = dum;
        }
    }
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t m = 0; m < o_; m++) {
            double dum = 0.0;
            for (size_t j = 0; j < o_; j++) {
                for (size_t b = 0; b < v_; b++) {
                    dum += tmp1_[j*v_+b] * t2_[b*o_*o_*v_+e*o_*o_+j*o_+m];
                }
            }
            ru1_[e*o_+m] += dum;
        }
    }

    // - 0.50000 <i,j||b,a> t2(b,e,i,j) u1(a,m)
    // I(i,j,b,m) = u1(a,m) <i,j||b,a>
    F_DGEMM('n','n',o_,o_*o_*v_,v_,1.0,u1_,o_,eri_ijab_,v_,0.0,tmp1_,o_);
    // t'(e,i,j,b) = t2(b,e,i,j)
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t i = 0; i < o_; i++) {
            for (size_t j = 0; j < o_; j++) {
                for (size_t b = 0; b < v_; b++) {
                    tmp2_[e*o_*o_*v_+i*o_*v_+j*v_+b] = t2_[b*o_*o_*v_+e*o_*o_+i*o_+j];
                }
            }
        }
    }
    // r(e,m) = I(i,j,b,m) t'(e,i,j,b)
    F_DGEMM('n','n',o_,v_,o_*o_*v_,-0.5,tmp1_,o_,tmp2_,o_*o_*v_,1.0,ru1_,o_);

    // + 0.50000 <i,j||a,b> t2(a,b,j,m) u1(e,i)
#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    tmp1_[a*o_*o_*v_+b*o_*o_+j*o_+i] = eri_ijab_[i*o_*v_*v_+j*v_*v_+a*v_+b];
                }
            }
        }
    }
    // I(m,i) = v'(a,b,j,i) t2(a,b,j,m)
    F_DGEMM('n','t',o_,o_,o_*v_*v_,1.0,tmp1_,o_,t2_,o_,0.0,tmp2_,o_);
    // r(e,m) = I(m,i) u1(e,i)
    F_DGEMM('t','n',o_,v_,o_,0.5,tmp2_,o_,u1_,o_,1.0,ru1_,o_);

    if ( include_u0_ ) {

        // + 1.00000 d-(i,m) u1(e,i) u0
        for (size_t i = 0; i < o_; i++) {
            for (size_t m = 0; m < o_; m++) {
                tmp1_[i*o_+m] = dp[i][m];
            }
        }
        F_DGEMM('n','n',o_,v_,o_,u0_[0],tmp1_,o_,u1_,o_,1.0,ru1_,o_);

        // - 1.00000 d-(e,a) u1(a,m) u0
        for (size_t e = 0; e < v_; e++) {
            for (size_t a = 0; a < v_; a++) {
                tmp1_[e*v_+a] = dp[e+o_][a+o_];
            }
        }
        F_DGEMM('n','n',o_,v_,v_,-u0_[0],u1_,o_,tmp1_,v_,1.0,ru1_,o_);

    }

    // - 1.00000 d-(i,a) u1(a,i) u1(e,m)
    double Iia = 0.0;
    for (size_t a = 0; a < v_; a++) {
        for (size_t i = 0; i < o_; i++) {
            Iia += dp[i][a+o_] * u1_[a*o_+i];
        }
    }
    C_DAXPY(o_*v_,-Iia,u1_,1,ru1_,1);

    // + 2.00000 d-(i,a) u1(a,m) u1(e,i)
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t a = 0; a < v_; a++) {
            tmp1_[i*v_+a] = dp[i][a+o_];
        }
    }
    // I(i,m) = u1_(a,m) d-(i,a)
    F_DGEMM('n','n',o_,o_,v_,1.0,u1_,o_,tmp1_,v_,0.0,tmp2_,o_);
    // r(e,m) = 2 I(i,m) u1(e,i)
    F_DGEMM('n','n',o_,v_,o_,2.0,tmp2_,o_,u1_,o_,1.0,ru1_,o_);

    if ( include_u2_ ) {

        // + 1.00000 F(i,a) u2(a,e,i,m)
        // - 1.00000 d-(i,a) u2(a,e,i,m) u0
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t m = 0; m < o_; m++) {
                double dum = 0.0;
                for (size_t a = 0; a < v_; a++) {
                    for (size_t i = 0; i < o_; i++) {
                        dum += u2_[a*o_*o_*v_+e*o_*o_+i*o_+m] * fp[i][a+o_];
                    }
                }
                ru1_[e*o_+m] += dum;

                if ( include_u0_ ) {

                    dum = 0.0;
                    for (size_t a = 0; a < v_; a++) {
                        for (size_t i = 0; i < o_; i++) {
                            dum -= u2_[a*o_*o_*v_+e*o_*o_+i*o_+m] * dp[i][a+o_];
                        }
                    }
                    ru1_[e*o_+m] += dum * u0_[0];

                }
            }
        }

// TODO: refactor o^3v integral terms

        // + 0.50000 <i,j||m,a> u2(a,e,i,j)
        // u'(e,i,j,a) = u2(a,e,i,j)
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    for (size_t a = 0; a < v_; a++) {
                        tmp2_[e*o_*o_*v_+i*o_*v_+j*v_+a] = u2_[a*o_*o_*v_+e*o_*o_+i*o_+j];
                    }
                }
            }
        }
        // v'(i,j,a,m) = <i,j||m,a>
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < o_; i++) {
            for (size_t j = 0; j < o_; j++) {
                for (size_t a = 0; a < v_; a++) {
                    for (size_t m = 0; m < o_; m++) {
                        tmp1_[i*o_*o_*v_+j*o_*v_+a*o_+m] = eri_jkia_[i*o_*o_*v_+j*o_*v_+m*v_+a];
                    }
                }
            }
        }
        // r(e,m) = v'(i,j,a,m) u'(e,i,j,a)
        F_DGEMM('n','n',o_,v_,o_*o_*v_,0.5,tmp1_,o_,tmp2_,o_*o_*v_,1.0,ru1_,o_);


        // - 0.50000 <e,i||a,b> u2(a,b,i,m)
        // u'(i,a,b,m) = u2(a,b,i,m)
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < o_; i++) {
            for (size_t a = 0; a < v_; a++) {
                for (size_t b = 0; b < v_; b++) {
                    for (size_t m = 0; m < o_; m++) {
                        tmp3_[i*o_*v_*v_+a*o_*v_+b*o_+m] = u2_[a*o_*o_*v_+b*o_*o_+i*o_+m];
                    }
                }
            }
        }

        // r(e,m) = -0.5 u'(i,a,b,m) <e,i||a,b>
	if ( is_hubbard_ ) { 
            F_DGEMM('n','n',o_,v_,o_*v_*v_,-0.5,tmp3_,o_,eri_aibc_,o_*v_*v_,1.0,ru1_,o_);
	}else {

            // four terms to avoid storing <e,i||a,b>...

            // r(e,m) = -0.5 u'(i,a,b,m) [(ea|Q)(Q|ib) - (eb|Q)(Q|ia) + d(ea)d(ib) - d(eb)d(ia) ]

            // 1: -0.5 u'(i,a,b,m) (ea|ib)
            // 2: +0.5 u'(i,a,b,m) (eb|ia)

            // u''(i,b,a,m) = u'(i,a,b,m) = u2(a,b,i,m)
#pragma omp parallel for schedule(static)
            for (size_t i = 0; i < o_; i++) {
                for (size_t b = 0; b < v_; b++) {
                    for (size_t a = 0; a < v_; a++) {
                        for (size_t m = 0; m < o_; m++) {
                            //tmp2_[i*o_*v_*v_+b*o_*v_+a*o_+m] = u2_[a*o_*o_*v_+b*o_*o_+i*o_+m];
                            tmp2_[i*o_*v_*v_+b*o_*v_+a*o_+m] = u2_[a*o_*o_*v_+b*o_*o_+i*o_+m] - u2_[b*o_*o_*v_+a*o_*o_+i*o_+m];
                        }
                    }
                }
            }
            // I(Q,a,m) = u''(i,b,a,m) (Q,i,b)
            F_DGEMM('n','n',o_*v_,nQ_,o_*v_,1.0,tmp2_,o_*v_,Qov_,o_*v_,0.0,tmp3_,o_*v_);

            // overwrite Qvv: I'(Q,a,e) = (Q|ea)
#pragma omp parallel for schedule(static)
            for (size_t q = 0; q < nQ_; q++) {
                for (size_t a = 0; a < v_; a++) {
                    for (size_t e = 0; e < v_; e++) {
                        tmp1_[q*v_*v_+a*v_+e] = Qvv_[q*v_*v_+e*v_+a];
                    }
                }
            }

            // r(e,m) = -0.5 I(Q,a,m) I'(Q,a,e)
            F_DGEMM('n','t',o_,v_,nQ_*v_,-0.5,tmp3_,o_,tmp1_,v_,1.0,ru1_,o_);

/*
            // 2: +0.5 u'(i,a,b,m) (eb|ia)

            // u'(i,a,b,m) = u2(a,b,i,m)
#pragma omp parallel for schedule(static)
            for (size_t i = 0; i < o_; i++) {
                for (size_t a = 0; a < v_; a++) {
                    for (size_t b = 0; b < v_; b++) {
                        for (size_t m = 0; m < o_; m++) {
                            tmp3_[i*o_*v_*v_+a*o_*v_+b*o_+m] = u2_[a*o_*o_*v_+b*o_*o_+i*o_+m];
                        }
                    }
                }
            }
            // I(Q,b,m) = +0.5 u'(i,a,b,m) (Q|i,a)
            F_DGEMM('n','n',o_*v_,nQ_,o_*v_,0.5,tmp3_,o_*v_,Qov_,o_*v_,0.0,tmp2_,o_*v_);

            // r(e,m) = I(Q,b,m) (Q|e,b) 
            // r(e,m) = I(Q,b,m) I'(Q,b,e) ... transposed Qvv still in tmp1
            F_DGEMM('n','t',o_,v_,nQ_*v_,1.0,tmp2_,o_,tmp1_,v_,1.0,ru1_,o_);
*/

            // 3: u'(i,a,b,m) d(e,a) d(i,b)
            double ** dz = Dipole_z_->pointer();
            double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);
            double lz2 = lambda_z * lambda_z;

           // I(a,m) = u'(i,a,b,m) d(i,b)
#pragma omp parallel for schedule(static)
           for (size_t a = 0; a < v_; a++) {
               for (size_t m = 0; m < o_; m++) {
                   double dum = 0.0;
                   for (size_t i = 0; i < o_; i++) {
                       for (size_t b = 0; b < v_; b++) {
                           //dum += tmp3_[i*o_*v_*v_+a*o_*v_+b*o_+m] * dz[i][b+o_];
                           dum += tmp2_[i*o_*v_*v_+b*o_*v_+a*o_+m] * dz[i][b+o_];
                       }
                   }
                   tmp1_[a*o_+m] = dum;
               }
           }

           // r(e,m) = -0.5 I(a,m) d(e,a) ... TODO could replace with DGEMM, but likely doesn't matter
#pragma omp parallel for schedule(static)
           for (size_t e = 0; e < v_; e++) {
               for (size_t m = 0; m < o_; m++) {
                   double dum = 0.0;
                   for (size_t a = 0; a < v_; a++) {
                       dum += tmp1_[a*o_+m] * dz[e+o_][a+o_];
                   }
                   ru1_[e*o_+m] -= 0.5 * lz2 * dum;
               }
           }

/*
           // 4: +0.5 u'(i,a,b,m) d(e,b) d(i,a)

           // I(b,m) = u'(i,a,b,m) d(i,a)
#pragma omp parallel for schedule(static)
           for (size_t b = 0; b < v_; b++) {
               for (size_t m = 0; m < o_; m++) {
                   double dum = 0.0;
                   for (size_t i = 0; i < o_; i++) {
                       for (size_t a = 0; a < v_; a++) {
                           dum += tmp3_[i*o_*v_*v_+a*o_*v_+b*o_+m] * dz[i][a+o_];
                       }
                   }
                   tmp1_[b*o_+m] = dum;
               }
           }

           // r(e,m) = +0.5 I(b,m) d(e,b) ... TODO could replace with DGEMM, but likely doesn't matter
#pragma omp parallel for schedule(static)
           for (size_t e = 0; e < v_; e++) {
               for (size_t m = 0; m < o_; m++) {
                   double dum = 0.0;
                   for (size_t b = 0; b < v_; b++) {
                       dum += tmp1_[b*o_+m] * dz[e+o_][b+o_];
                   }
                   ru1_[e*o_+m] += 0.5 * lz2 * dum;
               }
           }
*/

        }

    }

}

// 
// u2 residual
// 

// + 1.00000 u2(e,f,m,n) w0

// + 1.00000 F(i,n) u2(e,f,i,m)
// - 1.00000 F(i,m) u2(e,f,i,n)
// + 1.00000 F(e,a) u2(a,f,m,n)
// - 1.00000 F(f,a) u2(a,e,m,n)

// + 0.50000 <i,j||m,n> u2(e,f,i,j)
// + 0.50000 <e,f||a,b> u2(a,b,m,n)

// - 1.00000 <i,e||a,n> u2(a,f,i,m)
// + 1.00000 <i,e||a,m> u2(a,f,i,n)
// + 1.00000 <i,f||a,n> u2(a,e,i,m)
// - 1.00000 <i,f||a,m> u2(a,e,i,n)

// - 0.50000 <i,j||a,b> t2(e,f,i,m) u2(a,b,j,n)
// + 0.50000 <i,j||a,b> t2(e,f,i,n) u2(a,b,j,m)
// + 0.25000 <i,j||a,b> t2(e,f,i,j) u2(a,b,m,n)
// - 0.50000 <i,j||a,b> t2(b,f,m,n) u2(a,e,i,j)
// + 1.00000 <i,j||a,b> t2(b,f,i,m) u2(a,e,j,n)
// - 1.00000 <i,j||a,b> t2(b,f,i,n) u2(a,e,j,m)
// - 0.50000 <i,j||a,b> t2(b,f,i,j) u2(a,e,m,n)
// + 0.50000 <i,j||a,b> t2(b,e,m,n) u2(a,f,i,j)
// - 1.00000 <i,j||a,b> t2(b,e,i,m) u2(a,f,j,n)
// + 1.00000 <i,j||a,b> t2(b,e,i,n) u2(a,f,j,m)
// + 0.50000 <i,j||a,b> t2(b,e,i,j) u2(a,f,m,n)
// + 0.25000 <i,j||a,b> t2(a,b,m,n) u2(e,f,i,j)
// - 0.50000 <i,j||a,b> t2(a,b,i,m) u2(e,f,j,n)
// + 0.50000 <i,j||a,b> t2(a,b,i,n) u2(e,f,j,m)

// + 1.00000 <i,e||m,n> u1(f,i)
// - 1.00000 <i,f||m,n> u1(e,i)
// + 1.00000 <e,f||a,n> u1(a,m)
// - 1.00000 <e,f||a,m> u1(a,n)

// + 1.00000 F(i,a) t2(e,f,i,m) u1(a,n)
// - 1.00000 F(i,a) t2(e,f,i,n) u1(a,m)
// - 1.00000 F(i,a) t2(a,f,m,n) u1(e,i)
// + 1.00000 F(i,a) t2(a,e,m,n) u1(f,i)

// + 1.00000 <i,j||a,n> t2(e,f,j,m) u1(a,i)
// + 0.50000 <i,j||a,n> t2(e,f,i,j) u1(a,m)
// - 1.00000 <i,j||a,n> t2(a,f,j,m) u1(e,i)
// + 1.00000 <i,j||a,n> t2(a,e,j,m) u1(f,i)
// - 1.00000 <i,j||a,m> t2(e,f,j,n) u1(a,i)
// - 0.50000 <i,j||a,m> t2(e,f,i,j) u1(a,n)
// + 1.00000 <i,j||a,m> t2(a,f,j,n) u1(e,i)
// - 1.00000 <i,j||a,m> t2(a,e,j,n) u1(f,i)
// + 1.00000 <i,e||a,b> t2(b,f,m,n) u1(a,i)
// + 1.00000 <i,e||a,b> t2(b,f,i,m) u1(a,n)
// - 1.00000 <i,e||a,b> t2(b,f,i,n) u1(a,m)
// + 0.50000 <i,e||a,b> t2(a,b,m,n) u1(f,i)
// - 1.00000 <i,f||a,b> t2(b,e,m,n) u1(a,i)
// - 1.00000 <i,f||a,b> t2(b,e,i,m) u1(a,n)
// + 1.00000 <i,f||a,b> t2(b,e,i,n) u1(a,m)
// - 0.50000 <i,f||a,b> t2(a,b,m,n) u1(e,i)

// - 1.00000 d+(i,n) t2(e,f,i,m)
// + 1.00000 d+(i,m) t2(e,f,i,n)
// - 1.00000 d+(e,a) t2(a,f,m,n)
// + 1.00000 d+(f,a) t2(a,e,m,n)

// - 1.00000 d-(i,n) u1(e,i) u1(f,m)
// + 1.00000 d-(i,n) u1(e,m) u1(f,i)
// + 1.00000 d-(i,m) u1(e,i) u1(f,n)
// - 1.00000 d-(i,m) u1(e,n) u1(f,i)
// + 1.00000 d-(e,a) u1(a,n) u1(f,m)
// - 1.00000 d-(e,a) u1(a,m) u1(f,n)
// - 1.00000 d-(f,a) u1(a,n) u1(e,m)
// + 1.00000 d-(f,a) u1(a,m) u1(e,n)

// - 1.00000 d-(i,a) u1(a,i) u2(e,f,m,n)
// + 2.00000 d-(i,a) u1(a,n) u2(e,f,m,i)
// - 2.00000 d-(i,a) u1(a,m) u2(e,f,n,i)
// + 2.00000 d-(i,a) u1(e,i) u2(a,f,m,n)
// - 1.00000 d-(i,a) u1(e,n) u2(a,f,m,i)
// + 1.00000 d-(i,a) u1(e,m) u2(a,f,n,i)
// - 2.00000 d-(i,a) u1(f,i) u2(a,e,m,n)
// + 1.00000 d-(i,a) u1(f,n) u2(a,e,m,i)
// - 1.00000 d-(i,a) u1(f,m) u2(a,e,n,i)

// - 1.00000 d-(i,n) u2(e,f,i,m) u0
// + 1.00000 d-(i,m) u2(e,f,i,n) u0
// - 1.00000 d-(e,a) u2(a,f,m,n) u0
// + 1.00000 d-(f,a) u2(a,e,m,n) u0

void PolaritonicUCCSD::residual_u2() {

    memset((void*)ru2_,'\0',o_*o_*v_*v_*sizeof(double));

    // + 1.00000 u2(e,f,m,n) w0
    // TODO: generalize for x,y,z
    double w0 = cavity_frequency_[2];
    C_DAXPY(o_*o_*v_*v_,w0,u2_,1,ru2_,1);

    double ** fp = F_->pointer();

    double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];
    std::shared_ptr<Matrix> dip (new Matrix(Dipole_z_));
    dip->scale(coupling_factor_z);
    double ** dp = dip->pointer();

    // + F(i,m) u2(e,f,n,i)
    // - F(i,n) u2(e,f,m,i)
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t m = 0; m < o_; m++) {
            tmp1_[i*o_+m] = fp[i][m];
        }
    }
    // r'(e,f,n,m) =  F(i,m) u2(e,f,n,i)
    F_DGEMM('n','n',o_,o_*v_*v_,o_,1.0,tmp1_,o_,u2_,o_,0.0,tmp2_,o_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[e*o_*o_*v_+f*o_*o_+n*o_+m];
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[e*o_*o_*v_+f*o_*o_+m*o_+n];
                }
            }
        }
    }
    // + F(e,a) u2(a,f,m,n)
    // - F(f,a) u2(a,e,m,n)
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t a = 0; a < v_; a++) {
            tmp1_[e*v_+a] = fp[e+o_][a+o_];
        }
    }
    // r'(e,f,n,m) =  u2(a,f,m,n) F(e,a)
    F_DGEMM('n','n',o_*o_*v_,v_,v_,1.0,u2_,o_*o_*v_,tmp1_,v_,0.0,tmp2_,o_*o_*v_);
    C_DAXPY(o_*o_*v_*v_,1.0,tmp2_,1,ru2_,1);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[f*o_*o_*v_+e*o_*o_+m*o_+n];
                }
            }
        }
    }

    // + 0.5 <i,j||m,n> u2(e,f,i,j)
    F_DGEMM('n','n',o_*o_,v_*v_,o_*o_,0.5,eri_ijkl_,o_*o_,u2_,o_*o_,1.0,ru2_,o_*o_);

    // + 0.5 <e,f||a,b> u2(a,b,m,n) = <e,f|a,b> u2(a,b,m,n)
    if ( is_hubbard_ ) {
        F_DGEMM('n','n',o_*o_,v_*v_,v_*v_,1.0,u2_,o_*o_,eri_abcd_,v_*v_,1.0,ru2_,o_*o_);
    }else {
        double_particle_ladder_diagram(u2_,ru2_);
    }

    // - 1.00000 <i,e||n,a> u2(a,f,m,i)
    // + 1.00000 <i,e||m,a> u2(a,f,n,i)
    // + 1.00000 <i,f||n,a> u2(a,e,m,i)
    // - 1.00000 <i,f||m,a> u2(a,e,n,i)
    //
    // or
    //
    // - P(e,f) P(m,n) <i,e||n,a> u2(a,f,m,i)

    // u'(a,i,f,m) = u2(a,f,m,i)
#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v_; a++) {
        for (size_t i = 0; i < o_; i++) {
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    tmp1_[a*o_*o_*v_+i*o_*v_+f*o_+m] = u2_[a*o_*o_*v_+f*o_*o_+m*o_+i];
                }
            }
        }
    }

    // v'(a,i,e,n) = <ie||na>
#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v_; a++) {
        for (size_t i = 0; i < o_; i++) {
            for (size_t e = 0; e < v_; e++) {
                for (size_t n = 0; n < o_; n++) {
                    tmp2_[a*o_*o_*v_+i*o_*v_+e*o_+n] = eri_iajb_[i*o_*v_*v_+e*o_*v_+n*v_+a];
                }
            }
        }
    }

    // I(e,n,f,m) = u'(a,i,f,m) v'(a,i,e,n) 
    F_DGEMM('n','t',o_*v_,o_*v_,o_*v_,1.0,tmp1_,o_*v_,tmp2_,o_*v_,0.0,tmp3_,o_*v_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {

                    double dum = 0.0;
                    dum -= tmp3_[e*o_*o_*v_+n*o_*v_+f*o_+m];
                    dum += tmp3_[e*o_*o_*v_+m*o_*v_+f*o_+n];
                    dum += tmp3_[f*o_*o_*v_+n*o_*v_+e*o_+m];
                    dum -= tmp3_[f*o_*o_*v_+m*o_*v_+e*o_+n];

                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += dum;
                }
            }
        }
    }

    // - 0.5        <i,j||a,b> u2(a,b,n,j) t2(e,f,m,i)
    // + 0.5        <i,j||a,b> u2(a,b,m,j) t2(e,f,n,i)
    //
    // or
    //
    // - 0.5 P(m,n) <i,j||a,b> u2(a,b,n,j) t2(e,f,m,i) 
    //


    // u'(n,j,a,b) = u2(a,b,n,j)
#pragma omp parallel for schedule(static)
    for (size_t n = 0; n < o_; n++) {
        for (size_t j = 0; j < o_; j++) {
            for (size_t a = 0; a < v_; a++) {
                for (size_t b = 0; b < v_; b++) {
                    tmp1_[n*o_*v_*v_+j*v_*v_+a*v_+b] = u2_[a*o_*o_*v_+b*o_*o_+n*o_+j];
                }
            }
        }
    }

    // I(i,n) = u'(n,j,a,b) <i,j||a,b> 
    F_DGEMM('t','n',o_,o_,o_*v_*v_,1.0,tmp1_,o_*v_*v_,eri_ijab_,o_*v_*v_,0.0,tmp2_,o_);

    // r(e,f,m,n) = 0.5 I(i,n) t2(e,f,m,i)
    F_DGEMM('n','n',o_,o_*v_*v_,o_,0.5,tmp2_,o_,t2_,o_,0.0,tmp1_,o_);

    C_DAXPY(o_*o_*v_*v_,-1.0,tmp1_,1,ru2_,1);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp1_[e*o_*o_*v_+f*o_*o_+n*o_+m];
                }
            }
        }
    }

    // this next term is the same as the last one, swapping u2/t2

    // - 0.5        <i,j||a,b> t2(a,b,n,j) u2(e,f,m,i)
    // + 0.5        <i,j||a,b> t2(a,b,m,j) u2(e,f,n,i)
    //
    // or
    //
    // - 0.5 P(m,n) <i,j||a,b> t2(a,b,n,j) u2(e,f,m,i) 
    //


    // t'(n,j,a,b) = t2(a,b,n,j)
#pragma omp parallel for schedule(static)
    for (size_t n = 0; n < o_; n++) {
        for (size_t j = 0; j < o_; j++) {
            for (size_t a = 0; a < v_; a++) {
                for (size_t b = 0; b < v_; b++) {
                    tmp1_[n*o_*v_*v_+j*v_*v_+a*v_+b] = t2_[a*o_*o_*v_+b*o_*o_+n*o_+j];
                }
            }
        }
    }

    // I(i,n) = t'(n,j,a,b) <i,j||a,b> 
    F_DGEMM('t','n',o_,o_,o_*v_*v_,1.0,tmp1_,o_*v_*v_,eri_ijab_,o_*v_*v_,0.0,tmp2_,o_);

    // r(e,f,m,n) = 0.5 I(i,n) u2(e,f,m,i)
    F_DGEMM('n','n',o_,o_*v_*v_,o_,0.5,tmp2_,o_,u2_,o_,0.0,tmp1_,o_);

    C_DAXPY(o_*o_*v_*v_,-1.0,tmp1_,1,ru2_,1);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp1_[e*o_*o_*v_+f*o_*o_+n*o_+m];
                }
            }
        }
    }

    // + 0.25000 <i,j||a,b> t2(a,b,m,n) u2(e,f,i,j)

    // I(i,j,m,n) = t2(a,b,m,n) <i,j||a,b>
    F_DGEMM('n','n',o_*o_,o_*o_,v_*v_,1.0,t2_,o_*o_,eri_ijab_,v_*v_,0.0,tmp1_,o_*o_);

    // r(e,f,m,n) = I(i,j,m,n) u2(e,f,i,j)
    F_DGEMM('n','n',o_*o_,v_*v_,o_*o_,0.25,tmp1_,o_*o_,u2_,o_*o_,1.0,ru2_,o_*o_);

    // + 0.25000 <i,j||a,b> u2(a,b,m,n) t2(e,f,i,j)

    // I(i,j,m,n) = u2(a,b,m,n) <i,j||a,b>
    F_DGEMM('n','n',o_*o_,o_*o_,v_*v_,1.0,u2_,o_*o_,eri_ijab_,v_*v_,0.0,tmp1_,o_*o_);

    // r(e,f,m,n) = I(i,j,m,n) t2(e,f,i,j)
    F_DGEMM('n','n',o_*o_,v_*v_,o_*o_,0.25,tmp1_,o_*o_,t2_,o_*o_,1.0,ru2_,o_*o_);

    // - 0.5       <i,j||a,b> t2(a,e,m,n) u2(b,f,i,j)
    // + 0.5       <i,j||a,b> t2(a,f,m,n) u2(b,e,i,j)
    //
    // or
    //
    // - 0.5 P(e,f)<i,j||a,b> t2(a,e,m,n) u2(b,f,i,j) 

    // v'(a,b,i,j) = <i,j||a,b>
#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    tmp1_[a*o_*o_*v_+b*o_*o_+i*o_+j] = eri_ijab_[i*o_*v_*v_+j*v_*v_+a*v_+b];
                }
            }
        }
    }

    // u'(f,b,i,j) = u2(b,f,i,j)
#pragma omp parallel for schedule(static)
    for (size_t f = 0; f < v_; f++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    tmp2_[f*o_*o_*v_+b*o_*o_+i*o_+j] = u2_[b*o_*o_*v_+f*o_*o_+i*o_+j];
                }
            }
        }
    }

    // I(f,a) = v'(a,b,i,j) u'(f,b,i,j)
    F_DGEMM('t','n',v_,v_,o_*o_*v_,1.0,tmp1_,o_*o_*v_,tmp2_,o_*o_*v_,0.0,tmp3_,v_);

    // r(f,e,m,n) = 0.5 t2(a,e,m,n) I(f,a)
    F_DGEMM('n','n',o_*o_*v_,v_,v_,0.5,t2_,o_*o_*v_,tmp3_,v_,0.0,tmp1_,o_*o_*v_);

    C_DAXPY(o_*o_*v_*v_,1.0,tmp1_,1,ru2_,1);
#pragma omp parallel for schedule(static)
    for (size_t f = 0; f < v_; f++) {
        for (size_t e = 0; e < v_; e++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[f*o_*o_*v_+e*o_*o_+m*o_+n] -= tmp1_[e*o_*o_*v_+f*o_*o_+m*o_+n];
                }  
            }  
        }  
    }

    // this next term is the same as the last one, swapping u2/t2

    // - 0.5       <i,j||a,b> u2(a,e,m,n) t2(b,f,i,j) 
    // + 0.5       <i,j||a,b> u2(a,f,m,n) t2(b,e,i,j) 
    //
    // or
    //
    // - 0.5 P(e,f)<i,j||a,b> u2(a,e,m,n) t2(b,f,i,j) 

    // v'(a,b,i,j) = <i,j||a,b>
#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    tmp1_[a*o_*o_*v_+b*o_*o_+i*o_+j] = eri_ijab_[i*o_*v_*v_+j*v_*v_+a*v_+b];
                }
            }
        }
    }

    // t'(f,b,i,j) = t2(b,f,i,j)
#pragma omp parallel for schedule(static)
    for (size_t f = 0; f < v_; f++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    tmp2_[f*o_*o_*v_+b*o_*o_+i*o_+j] = t2_[b*o_*o_*v_+f*o_*o_+i*o_+j];
                }
            }
        }
    }

    // I(f,a) = v'(a,b,i,j) t'(f,b,i,j)
    F_DGEMM('t','n',v_,v_,o_*o_*v_,1.0,tmp1_,o_*o_*v_,tmp2_,o_*o_*v_,0.0,tmp3_,v_);

    // r(f,e,m,n) = 0.5 u2(a,e,m,n) I(f,a)
    F_DGEMM('n','n',o_*o_*v_,v_,v_,0.5,u2_,o_*o_*v_,tmp3_,v_,0.0,tmp1_,o_*o_*v_);

    C_DAXPY(o_*o_*v_*v_,1.0,tmp1_,1,ru2_,1);
#pragma omp parallel for schedule(static)
    for (size_t f = 0; f < v_; f++) {
        for (size_t e = 0; e < v_; e++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[f*o_*o_*v_+e*o_*o_+m*o_+n] -= tmp1_[e*o_*o_*v_+f*o_*o_+m*o_+n];
                }  
            }  
        }  
    }

    // +        <i,j||a,b> u2(a,e,n,j) t2(b,f,m,i) 
    // -        <i,j||a,b> u2(a,e,m,j) t2(b,f,n,i) 
    //
    // or
    //
    // + P(m,n) <i,j||a,b> u2(a,e,n,j) t2(b,f,m,i) 

    //I(i,b,j,a) = <i,j||a,b>
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t j = 0; j < o_; j++) {
            for (size_t a = 0; a < v_; a++) {
                for (size_t b = 0; b < v_; b++) {
                    tmp1_[i*o_*v_*v_+b*o_*v_+j*v_+a] = eri_ijab_[i*o_*v_*v_+j*v_*v_+a*v_+b];
                }
            }
        }
    }

    // t'(i,b,m,f) = t2(b,f,m,i)
#pragma omp parallel for schedule(static)
    for (size_t b = 0; b < v_; b++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t i = 0; i < o_; i++) {
                    tmp2_[i*o_*v_*v_+b*o_*v_+m*v_+f] = t2_[b*o_*o_*v_+f*o_*o_+m*o_+i];
                }
            }
        }
    }

    //I''(m,f,j,a) = I(i,b,j,a) t'(i,b,m,f) 
    F_DGEMM('n','t',o_*v_,o_*v_,o_*v_,1.0,tmp1_,o_*v_,tmp2_,o_*v_,0.0,tmp3_,o_*v_);

    // u'(n,e,j,a) = u2(a,e,n,j) 
#pragma omp parallel for schedule(static)
    for (size_t n = 0; n < o_; n++) {
        for (size_t e = 0; e < v_; e++) {
            for (size_t j = 0; j < o_; j++) {
                for (size_t a = 0; a < v_; a++) {
                    tmp1_[e*o_*o_*v_+n*o_*v_+j*v_+a] = u2_[a*o_*o_*v_+e*o_*o_+n*o_+j];
                }
            }
        }
    }

    //r'(e,n,m,f) = I''(m,f,j,a) u'(n,e,j,a)
    F_DGEMM('t','n',o_*v_,o_*v_,o_*v_,1.0,tmp3_,o_*v_,tmp1_,o_*v_,0.0,tmp2_,o_*v_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[e*o_*o_*v_+n*o_*v_+m*v_+f];
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[e*o_*o_*v_+m*o_*v_+n*v_+f];
                }
            }
        }
    }

    // this next term is the same as the last one, swapping u2/t2

    // +        <i,j||a,b> t2(a,e,n,j) u2(b,f,m,i) 
    // -        <i,j||a,b> t2(a,e,m,j) u2(b,f,n,i) 
    //
    // or
    //
    // + P(m,n) <i,j||a,b> t2(a,e,n,j) u2(b,f,m,i) 

    //I(i,b,j,a) = <i,j||a,b>
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t j = 0; j < o_; j++) {
            for (size_t a = 0; a < v_; a++) {
                for (size_t b = 0; b < v_; b++) {
                    tmp1_[i*o_*v_*v_+b*o_*v_+j*v_+a] = eri_ijab_[i*o_*v_*v_+j*v_*v_+a*v_+b];
                }
            }
        }
    }

    // u'(i,b,m,f) = u2(b,f,m,i)
#pragma omp parallel for schedule(static)
    for (size_t b = 0; b < v_; b++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t i = 0; i < o_; i++) {
                    tmp2_[i*o_*v_*v_+b*o_*v_+m*v_+f] = u2_[b*o_*o_*v_+f*o_*o_+m*o_+i];
                }
            }
        }
    }

    //I''(m,f,j,a) = I(i,b,j,a) u'(i,b,m,f) 
    F_DGEMM('n','t',o_*v_,o_*v_,o_*v_,1.0,tmp1_,o_*v_,tmp2_,o_*v_,0.0,tmp3_,o_*v_);

    // t'(n,e,j,a) = t2(a,e,n,j) 
#pragma omp parallel for schedule(static)
    for (size_t n = 0; n < o_; n++) {
        for (size_t e = 0; e < v_; e++) {
            for (size_t j = 0; j < o_; j++) {
                for (size_t a = 0; a < v_; a++) {
                    tmp1_[e*o_*o_*v_+n*o_*v_+j*v_+a] = t2_[a*o_*o_*v_+e*o_*o_+n*o_+j];
                }
            }
        }
    }

    //r'(e,n,m,f) = I''(m,f,j,a) t'(n,e,j,a)
    F_DGEMM('t','n',o_*v_,o_*v_,o_*v_,1.0,tmp3_,o_*v_,tmp1_,o_*v_,0.0,tmp2_,o_*v_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[e*o_*o_*v_+n*o_*v_+m*v_+f];
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[e*o_*o_*v_+m*o_*v_+n*v_+f];
                }
            }
        }
    }

    // - F(i,a) u1(a,n) t2(e,f,m,i) 
    // + F(i,a) u1(a,m) t2(e,f,n,i) 
    //
    // or
    //
    // - P(m,n) F(i,a) u1(a,n) t2(e,f,m,i)

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t a = 0; a < v_; a++) {
            tmp1_[i*v_+a] = fp[i][a+o_];
        }
    }
    // I(i,n) = u1(a,n) F(i,a)
    F_DGEMM('n','n',o_,o_,v_,1.0,u1_,o_,tmp1_,v_,0.0,tmp2_,o_);
    // r'(e,f,m,n) = -I(i,n) t2(e,f,m,i)
    F_DGEMM('n','n',o_,o_*v_*v_,o_,-1.0,tmp2_,o_,t2_,o_,0.0,tmp1_,o_);
    C_DAXPY(o_*o_*v_*v_,1.0,tmp1_,1,ru2_,1);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp1_[e*o_*o_*v_+f*o_*o_+n*o_+m];
                }
            }
        }
    }


    // - F(i,a) u1(e,i) t2(a,f,m,n)
    // + F(i,a) u1(f,i) t2(a,e,m,n)
    //
    // or
    //
    // - P(e,f) F(i,a) u1(e,i) t2(a,f,m,n)
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t a = 0; a < v_; a++) {
            tmp1_[i*v_+a] = fp[i][a+o_];
        }
    }
    // I(e,a) = F(i,a) u1(e,i)
    F_DGEMM('n','n',v_,v_,o_,1.0,tmp1_,v_,u1_,o_,0.0,tmp2_,v_);
    // r'(e,f,m,n) = -u2(a,f,m,n) I(e,a)
    F_DGEMM('n','n',o_*o_*v_,v_,v_,-1.0,t2_,o_*o_*v_,tmp2_,v_,0.0,tmp1_,o_*o_*v_);
    C_DAXPY(o_*o_*v_*v_,1.0,tmp1_,1,ru2_,1);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp1_[f*o_*o_*v_+e*o_*o_+m*o_+n];
                }
            }
        }
    }

// TODO: refactor o^3v integral terms

    // + 1.00000 <i,e||m,n> u1(f,i)
    // - 1.00000 <i,f||m,n> u1(e,i)
    //
    // or
    //
    // P(e,f) <i,e||m,n> u1(f,i)

    // r'(f,e,m,n) = <i,e||m,n> u1(f,i)
    F_DGEMM('n','n',o_*o_*v_,v_,o_,1.0,eri_iajk_,o_*o_*v_,u1_,o_,0.0,tmp1_,o_*o_*v_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp1_[f*o_*o_*v_+e*o_*o_+m*o_+n];
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp1_[e*o_*o_*v_+f*o_*o_+m*o_+n];
                }
            }
        }
    }

    // + 1.00000 <i,j||a,n> t2(e,f,j,m) u1(a,i)
    // - 1.00000 <i,j||a,m> t2(e,f,j,n) u1(a,i)
    //
    // or
    //
    // P(m,n) <i,j||n,a> t2(e,f,m,j) u1(a,i)

    // I(j,n) = <i,j||n,a> u1(a,i)
#pragma omp parallel for schedule(static)
    for (size_t j = 0; j < o_; j++) {
        for (size_t n = 0; n < o_; n++) {
            double dum = 0.0;
            for (size_t i = 0; i < o_; i++) {
                for (size_t a = 0; a < v_; a++) {
                    dum += eri_jkia_[i*o_*o_*v_+j*o_*v_+n*v_+a] * u1_[a*o_+i];
                }
            }
            tmp1_[j*o_+n] = dum;
        }
    }
    // r'(e,f,m,n) = I(j,n) t2(e,f,m,j)
    F_DGEMM('n','n',o_,o_*v_*v_,o_,1.0,tmp1_,o_,t2_,o_,0.0,tmp2_,o_);
    C_DAXPY(o_*o_*v_*v_,1.0,tmp2_,1,ru2_,1);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[e*o_*o_*v_+f*o_*o_+n*o_+m];
                }
            }
        }
    }

    // + 1.00000 <i,j||n,a> t2(a,e,m,j) u1(f,i)
    // - 1.00000 <i,j||m,a> t2(a,e,n,j) u1(f,i)
    // - 1.00000 <i,j||n,a> t2(a,f,m,j) u1(e,i)
    // + 1.00000 <i,j||m,a> t2(a,f,n,j) u1(e,i)
    //
    // or
    //
    // P(e,f) P(m,n) <i,j||n,a> t2(a,e,m,j) u1(f,i)

    // I(f,j,n,a) = <i,j||n,a> u1(f,i)
    F_DGEMM('n','n',o_*o_*v_,v_,o_,1.0,eri_jkia_,o_*o_*v_,u1_,o_,0.0,tmp1_,o_*o_*v_);

    // I'(a,j,f,n) = I(f,j,n,a)
#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v_; a++) {
        for (size_t j = 0; j < o_; j++) {
            for (size_t f = 0; f < v_; f++) {
                for (size_t n = 0; n < o_; n++) {
                    tmp2_[a*o_*o_*v_+j*o_*v_+f*o_+n] = tmp1_[f*o_*o_*v_+j*o_*v_+n*v_+a];
                }
            }
        }
    }
    // t'(e,m,a,j) = t2(a,e,m,j)
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t m = 0; m < o_; m++) {
            for (size_t a = 0; a < v_; a++) {
                for (size_t j = 0; j < o_; j++) {
                    tmp1_[e*o_*o_*v_+m*o_*v_+a*o_+j] = t2_[a*o_*o_*v_+e*o_*o_+m*o_+j];
                }
            }
        }
    }
    // r'(e,m,f,n) = I'(a,j,f,n) t'(e,m,a,j)
    F_DGEMM('n','n',o_*v_,o_*v_,o_*v_,1.0,tmp2_,o_*v_,tmp1_,o_*v_,0.0,tmp3_,o_*v_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp3_[e*o_*o_*v_+m*o_*v_+f*o_+n];
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp3_[e*o_*o_*v_+n*o_*v_+f*o_+m];
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp3_[f*o_*o_*v_+m*o_*v_+e*o_+n];
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp3_[f*o_*o_*v_+n*o_*v_+e*o_+m];
                }
            }
        }
    }

    // - 0.50000 <i,j||n,a> t2(e,f,i,j) u1(a,m)
    // + 0.50000 <i,j||m,a> t2(e,f,i,j) u1(a,n)
    //
    // or
    //
    // - 0.5 P(m,n) <i,j||n,a> t2(e,f,i,j) u1(a,m)

    // I(i,j,n,m) = u1(a,m) <i,j||n,a> 
    F_DGEMM('n','n',o_,o_*o_*o_,v_,1.0,u1_,o_,eri_jkia_,v_,0.0,tmp1_,o_);

    // r'(e,f,n,m) = -0.5 I(i,j,n,m) t2(e,f,i,j)
    F_DGEMM('n','n',o_*o_,v_*v_,o_*o_,-0.5,tmp1_,o_*o_,t2_,o_*o_,0.0,tmp2_,o_*o_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    // fixed sign bug here 8/11/21
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[e*o_*o_*v_+f*o_*o_+n*o_+m];
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[e*o_*o_*v_+f*o_*o_+m*o_+n];
                }
            }
        }
    }


    // - 1.00000 <e,f||n,a> u1(a,m)
    // + 1.00000 <e,f||m,a> u1(a,n)
    //
    // or
    //
    // + P(m,n) u1(a,n) <e,f||m,a> 
    //
    if ( is_hubbard_ ) {
        F_DGEMM('n','n',o_,o_*v_*v_,v_,1.0,u1_,o_,eri_abic_,v_,0.0,tmp1_,o_);
        C_DAXPY(o_*o_*v_*v_,1.0,tmp1_,1,ru2_,1);
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {
                        ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp1_[e*o_*o_*v_+f*o_*o_+n*o_+m];
                    }
                }
            }
        }
    }else {
        // P(m,n) u1(a,n) [ (em|fa) - (ea|fm) ]

        // I(Q,f,n) = u1(a,n) (Q|fa)
        F_DGEMM('n','n',o_,nQ_*v_,v_,1.0,u1_,o_,Qvv_,v_,0.0,tmp1_,o_);

        // P(m,n) [ I(Q,f,n) (Q|em) - I(Q,e,n) (Q|fm) ]
        //
        // or
        //
        // P(m,n) P(e,f) I(Q,f,n) (Q|em)

        // r'(e,m,f,n) = I(Q,f,n) (Q|em)
        F_DGEMM('n','t',o_*v_,o_*v_,nQ_,1.0,tmp1_,o_*v_,Qvo_,o_*v_,0.0,tmp2_,o_*v_);

#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {

                        double dum = 0.0;

                        dum += tmp2_[e*o_*o_*v_+m*o_*v_+f*o_+n];
                        dum -= tmp2_[e*o_*o_*v_+n*o_*v_+f*o_+m];
                        dum -= tmp2_[f*o_*o_*v_+m*o_*v_+e*o_+n];
                        dum += tmp2_[f*o_*o_*v_+n*o_*v_+e*o_+m];

                        ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += dum;
                    }
                }
            }
        }

        // and dipole terms: P(m,n) u(a,n) [ d(e,m) d(f,a) - d(e,a) d(m,f) ] = P(m,n) P(e,f) u(a,n) d(e,m) d(f,a)
        double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);
        double ** dz = Dipole_z_->pointer();

#pragma omp parallel for schedule(static)
        for (size_t f = 0; f < v_; f++) {
            for (size_t a = 0; a < v_; a++) {
                tmp2_[f*v_+a] = dz[f+o_][a+o_];
            }
        }
        // I(f,n) = u(a,n) d(f,a) * lambda^2
        F_DGEMM('n','n',o_,v_,v_,lambda_z * lambda_z,u1_,o_,tmp2_,v_,0.0,tmp1_,o_);

#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {

                        double dum = 0.0;

                        dum += tmp1_[f*o_+n] * dz[e+o_][m];
                        dum -= tmp1_[f*o_+m] * dz[e+o_][n];
                        dum -= tmp1_[e*o_+n] * dz[f+o_][m];
                        dum += tmp1_[e*o_+m] * dz[f+o_][n];

                        ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += dum;
                    }
                }
            }
        }
        

    }

    // - 1.00000 <e,i||a,b> t2(b,f,m,n) u1(a,i)
    // + 1.00000 <f,i||a,b> t2(b,e,m,n) u1(a,i)
    //
    // or
    //
    // - P(e,f) <e,i||a,b> t2(b,f,m,n) u1(a,i)

    // I(e,b) = <e,i||a,b> u1(a,i)
    if ( is_hubbard_ ) { 

#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t b = 0; b < v_; b++) {
                double dum = 0.0;
                for (size_t i = 0; i < o_; i++) {
                    for (size_t a = 0; a < v_; a++) {
                        dum += eri_aibc_[e*o_*v_*v_+i*v_*v_+a*v_+b] * u1_[a*o_+i];
                    }
                }
                tmp3_[e*v_+b] = dum;
            }
        }

    }else {

        // four terms to avoid storing <e,i||a,b>...

        // I(e,b) = [ (ea|Q)(Q|ib) - (eb|ia) + d(ea)d(ib) - d(eb)d(ia) ] u1(a,i)

        // 1: I(e,b) = u1(a,i) (Q|ea) (Q|ib)

        // I(Q,e,i) = u(a,i) (Q|e,a)
        F_DGEMM('n','n',o_,nQ_*v_,v_,1.0,u1_,o_,Qvv_,v_,0.0,tmp1_,o_);

        // I'(Q,i,e) = I(Q,e,i)
#pragma omp parallel for schedule(static)
        for (size_t q = 0; q < nQ_; q++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t e = 0; e < v_; e++) {
                    tmp1_[q*o_*v_+i*v_+e + nQ_*o_*v_] = tmp1_[q*o_*v_+e*o_+i];
                }
            }
        }

        // I(e,b) = (Q|ib) I'(Q,i,e)
        F_DGEMM('n','t',v_,v_,nQ_*o_,1.0,Qov_,v_,tmp1_+nQ_*o_*v_,v_,0.0,tmp3_,v_);

        // 2: I(e,b) = -u1(a,i) (Q|eb) (Q|ia)

        // I(Q) = u(a,i) (Q|ia)
#pragma omp parallel for schedule(static)
        for (size_t q = 0; q < nQ_; q++) {
            double dum = 0.0;
            for (size_t a = 0; a < v_; a++) {
                for (size_t i = 0; i < o_; i++) {
                    dum += u1_[a*o_+i] * Qov_[q*o_*v_+i*v_+a];
                }
            }
            tmp1_[q] = dum;
        }

        // I(e,b) = -I(Q) (Q|eb)
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t b = 0; b < v_; b++) {
                double dum = 0.0;
                for (size_t q = 0; q < nQ_; q++) {
                    dum += tmp1_[q] * Qvv_[q*v_*v_+e*v_+b];
                }
                tmp3_[e*v_+b] -= dum;
            }
        }

        // 3: I(e,b) = u(a,i) d(e,a) d(i,b)

        double ** dz = Dipole_z_->pointer();
        double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);
        double lz2 = lambda_z * lambda_z;

        // I'(e,i) = u(a,i) d(e,a)
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t i = 0; i < o_; i++) {
                double dum = 0.0;
                for (size_t a = 0; a < v_; a++) {
                    dum += u1_[a*o_+i] * dz[e+o_][a+o_];
                }
                tmp1_[e*o_+i] = dum;
            }
        }
  
        // I(e,b) += d(i,b) I(e,i)
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t b = 0; b < v_; b++) {
                double dum = 0.0;
                for (size_t i = 0; i < o_; i++) {
                    dum += tmp1_[e*o_+i] * dz[i][b+o_];
                }
                tmp3_[e*v_+b] += lz2 * dum;
            }
        }

        // 4: I(e,b) = -u(a,i) d(e,b) d(i,a)

        // I = u(a,i) d(i,a) 
        double dum = 0.0;
        for (size_t a = 0; a < v_; a++) {
            for (size_t i = 0; i < o_; i++) {
                dum += u1_[a*o_+i] * dz[i][a+o_];
            }
        }

        // r(e,m) -= I d(e,b)
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t b = 0; b < v_; b++) {
                tmp3_[e*v_+b] -= lz2 * dum * dz[e+o_][b+o_];
            }
        }


    }

    // r'(e,f,m,n) = t2(b,f,m,n) I(e,b)
    F_DGEMM('n','n',o_*o_*v_,v_,v_,-1.0,t2_,o_*o_*v_,tmp3_,v_,0.0,tmp2_,o_*o_*v_);
    C_DAXPY(o_*o_*v_*v_,1.0,tmp2_,1,ru2_,1);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[f*o_*o_*v_+e*o_*o_+m*o_+n];
                }
            }
        }
    }

    // - 1.00000 <e,i||a,b> t2(b,f,i,m) u1(a,n)
    // + 1.00000 <e,i||a,b> t2(b,f,i,n) u1(a,m)
    // + 1.00000 <f,i||a,b> t2(b,e,i,m) u1(a,n)
    // - 1.00000 <f,i||a,b> t2(b,e,i,n) u1(a,m)
    //
    // or
    //
    //   P(e,f) P(m,n) <e,i||a,b> t2(a,f,i,m) u1(b,n)

    // I(e,i,a,n) = u1(b,n) <e,i||a,b>
    if ( is_hubbard_ ) { 

        F_DGEMM('n','n',o_,o_*v_*v_,v_,1.0,u1_,o_,eri_aibc_,v_,0.0,tmp1_,o_);

        // I'(e,n,a,i) = I(e,i,a,n)
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t n = 0; n < o_; n++) {
                for (size_t a = 0; a < v_; a++) {
                    for (size_t i = 0; i < o_; i++) {
                        tmp2_[e*o_*o_*v_+n*o_*v_+a*o_+i] = tmp1_[e*o_*o_*v_+i*o_*v_+a*o_+n];
                    }
                }
            }
        }

    }else {

        // four terms to eliminate <e,i||a,b> in generating I'

        // 1: (Q|e,a) (Q|i,b) t2(a,f,i,m) u1(b,n)

        // I''(Q,i,n) = u1(b,n) (Q|i,b)
        F_DGEMM('n','n',o_,nQ_*o_,v_,1.0,u1_,o_,Qov_,v_,0.0,tmp1_,o_);

        // I'''(e,a,i,n) = I''(Q,i,n) (Q|e,a)
        F_DGEMM('n','t',o_*o_,v_*v_,nQ_,1.0,tmp1_,o_*o_,Qvv_,v_*v_,0.0,tmp3_,o_*o_);

        // I'(e,n,a,i) = I'''(e,a,i,n)
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t n = 0; n < o_; n++) {
                for (size_t a = 0; a < v_; a++) {
                    for (size_t i = 0; i < o_; i++) {
                        tmp2_[e*o_*o_*v_+n*o_*v_+a*o_+i] = tmp3_[e*o_*o_*v_+a*o_*o_+i*o_+n];
                    }
                }
            }
        }

        // 2: -(Q|e,b) (Q|i,a) t2(a,f,i,m) u1(b,n)

        // I''(Q,e,n) = u(b,n) (Q|e,b)
        F_DGEMM('n','n',o_,nQ_*v_,v_,1.0,u1_,o_,Qvv_,v_,0.0,tmp1_,o_);

        // I'''(i,a,e,n) = I''(Q,e,n) (Q|i,a)
        F_DGEMM('n','t',o_*v_,o_*v_,nQ_,1.0,tmp1_,o_*v_,Qov_,o_*v_,0.0,tmp3_,o_*v_);

        // I'(e,n,a,i) -= I'''(i,a,e,n)
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t n = 0; n < o_; n++) {
                for (size_t a = 0; a < v_; a++) {
                    for (size_t i = 0; i < o_; i++) {
                        tmp2_[e*o_*o_*v_+n*o_*v_+a*o_+i] -= tmp3_[i*o_*v_*v_+a*o_*v_+e*o_+n];
                    }
                }
            }
        }

        // 3: d(e,a) d(i,b) t2(a,f,i,m) u1(b,n)
        double ** dz = Dipole_z_->pointer();
        double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);
        double lz2 = lambda_z * lambda_z;


        // I(i,n) = u1(b,n) d(i,b)
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < o_; i++) {
            for (size_t n = 0; n < o_; n++) {
                double dum = 0.0;
                for (size_t b = 0; b < v_; b++) {
                    dum += u1_[b*o_+n] * dz[i][b+o_];
                }
                tmp1_[i*o_+n] = lz2 * dum;
            }
        }

        // I'(e,n,a,i) += I(i,n) d(e,a)
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t n = 0; n < o_; n++) {
                for (size_t a = 0; a < v_; a++) {
                    for (size_t i = 0; i < o_; i++) {
                        tmp2_[e*o_*o_*v_+n*o_*v_+a*o_+i] += tmp1_[i*o_+n] * dz[e+o_][a+o_];
                    }
                }
            }
        }

        // 4: -d(e,b) d(i,a) t2(a,f,i,m) u1(b,n)

        // I(e,n) = u1(b,n) d(e,b)
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t n = 0; n < o_; n++) {
                double dum = 0.0;
                for (size_t b = 0; b < v_; b++) {
                    dum += u1_[b*o_+n] * dz[e+o_][b+o_];
                }
                tmp1_[e*o_+n] = lz2 * dum;
            }
        }

        // I'(e,n,a,i) -= I(e,n) d(i,a)
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t n = 0; n < o_; n++) {
                for (size_t a = 0; a < v_; a++) {
                    for (size_t i = 0; i < o_; i++) {
                        tmp2_[e*o_*o_*v_+n*o_*v_+a*o_+i] -= tmp1_[e*o_+n] * dz[i][a+o_];
                    }
                }
            }
        }

    }
    // t'(a,i,f,m) = t2(a,f,i,m)
#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v_; a++) {
        for (size_t i = 0; i < o_; i++) {
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    tmp1_[a*o_*o_*v_+i*o_*v_+f*o_+m] = t2_[a*o_*o_*v_+f*o_*o_+i*o_+m];
                }
            }
        }
    }
    // r'(e,n,f,m) = t'(a,i,f,m) I'(e,n,a,i)
    F_DGEMM('n','n',o_*v_,o_*v_,o_*v_,1.0,tmp1_,o_*v_,tmp2_,o_*v_,0.0,tmp3_,o_*v_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    double dum = 0.0;
                    dum += tmp3_[e*o_*o_*v_+n*o_*v_+f*o_+m];
                    dum -= tmp3_[e*o_*o_*v_+m*o_*v_+f*o_+n];
                    dum -= tmp3_[f*o_*o_*v_+n*o_*v_+e*o_+m];
                    dum += tmp3_[f*o_*o_*v_+m*o_*v_+e*o_+n];
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += dum;
                }
            }
        }
    }

// TODO: refactor ov^3 integral terms

    // - 0.50000 <e,i||a,b> t2(a,b,m,n) u1(f,i)
    // + 0.50000 <f,i||a,b> t2(a,b,m,n) u1(e,i)
    //
    // or
    //
    // - 0.5 P(e,f) <e,i||a,b> t2(a,b,m,n) u1(f,i)

    // I(m,n,e,i) = <e,i||a,b> t2(a,b,m,n) 
    if ( is_hubbard_ ) { 

        F_DGEMM('t','t',o_*v_,o_*o_,v_*v_,1.0,eri_aibc_,v_*v_,t2_,o_*o_,0.0,tmp1_,o_*v_);

        // r'(m,n,e,f) = u1(f,i) I(m,n,e,i)
        F_DGEMM('t','n',v_,o_*o_*v_,o_,-0.5,u1_,o_,tmp1_,o_,0.0,tmp2_,v_);

#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {
                        ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[m*o_*v_*v_+n*v_*v_+e*v_+f];
                        ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[m*o_*v_*v_+n*v_*v_+f*v_+e];
                    }
                }
            }
        }
    }else {

        // four terms to avoid storing <e,i||a,b>

        // 1: +(Q|ea) (Q|ib) t2(a,b,m,n) u1(f,i)
        // 2: -(Q|eb) (Q|ia) t2(a,b,m,n) u1(f,i)

        // i think the only term that actually differs between these
        // two will be the contraction involving t2. all other 
        // intermediates are common, so let's try to do everything
        // at once.

        // need to loop over "e" to avoid storing large tensors. 
        // unfortunately, this term is still o^3v^3 cost, and much 
        // less efficient than it used to be because we're doing v o^3v^2 gemms
        for (size_t e = 0; e < v_; e++) {
           
            // I(Q,a) = (Q|e,a)
#pragma omp parallel for schedule(static)
            for (size_t q = 0; q < nQ_; q++) {
                for (size_t a = 0; a < v_; a++) {
                    tmp1_[q*v_+a] = Qvv_[q*v_*v_+e*v_+a];
                }
            }

            // I'(i,b,a) = I(Q,a) (Q|i,b)
            F_DGEMM('n','t',v_,o_*v_,nQ_,1.0,tmp1_,v_,Qov_,o_*v_,0.0,tmp2_,v_);

            // t'(b,a,m,n) = t2(a,b,m,n)
#pragma omp parallel for schedule(static)
            for (size_t b = 0; b < v_; b++) {
                for (size_t a = 0; a < v_; a++) {
                    for (size_t m = 0; m < o_; m++) {
                        for (size_t n = 0; n < o_; n++) {
                            tmp1_[b*o_*o_*v_+a*o_*o_+m*o_+n] = t2_[a*o_*o_*v_+b*o_*o_+m*o_+n] - t2_[b*o_*o_*v_+a*o_*o_+m*o_+n];
                        }
                    }
                }
            }

            // I''(i,m,n) = t'(b,a,m,n) I(i,b,a)
            F_DGEMM('n','n',o_*o_,o_,v_*v_,1.0,tmp1_,o_*o_,tmp2_,v_*v_,0.0,tmp3_,o_*o_);

            // I'''(f,m,n) = 0.5 I'''(i,m,n) u1(f,i)
            F_DGEMM('n','n',o_*o_,v_,o_,0.5,tmp3_,o_*o_,u1_,o_,0.0,tmp1_,o_*o_);

            // r(e,f,m,n) -= I'''(f,m,n) 
            // r(f,e,m,n) += I'''(f,m,n) 
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {
                        ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp1_[f*o_*o_+m*o_+n];
                        ru2_[f*o_*o_*v_+e*o_*o_+m*o_+n] += tmp1_[f*o_*o_+m*o_+n];
                    }
                }
            }

        }

/*
        // 2: -(Q|eb) (Q|ia) t2(a,b,m,n) u1(f,i) ... folded into term 1 above

        // need to loop over "e" to avoid storing large tensors. 
        // unfortunately, this term is still o^3v^3 cost, and much 
        // less efficient than it used to be because we're doing v o^3v^2 gemms
        for (size_t e = 0; e < v_; e++) {
           
            // I(Q,b) = (Q|e,b)
#pragma omp parallel for schedule(static)
            for (size_t q = 0; q < nQ_; q++) {
                for (size_t b = 0; b < v_; b++) {
                    tmp1_[q*v_+b] = Qvv_[q*v_*v_+e*v_+b]; // identical to term above
                }
            }

            // I'(i,a,b) = I(Q,b) (Q|i,a)
            F_DGEMM('n','t',v_,o_*v_,nQ_,1.0,tmp1_,v_,Qov_,o_*v_,0.0,tmp2_,v_); // identical to term above

            // I''(i,m,n) = t2(a,b,m,n) I(i,a,b)
            F_DGEMM('n','n',o_*o_,o_,v_*v_,1.0,t2_,o_*o_,tmp2_,v_*v_,0.0,tmp3_,o_*o_);

            // I'''(f,m,n) = -0.5 I'''(i,m,n) u1(f,i)
            F_DGEMM('n','n',o_*o_,v_,o_,-0.5,tmp3_,o_*o_,u1_,o_,0.0,tmp1_,o_*o_); // identical to term above (except for sign)

            // r(e,f,m,n) -= I'''(f,m,n) 
            // r(f,e,m,n) += I'''(f,m,n) 
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {
                        ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp1_[f*o_*o_+m*o_+n];
                        ru2_[f*o_*o_*v_+e*o_*o_+m*o_+n] += tmp1_[f*o_*o_+m*o_+n];
                    }
                }
            }

        }
*/

        // 3: +d(e,a) d(i,b) t2(a,b,m,n) u1(f,i)
        // 4: -d(e,b) d(i,a) t2(a,b,m,n) u1(f,i)

        double ** dz = Dipole_z_->pointer();
        double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);
        double lz2 = lambda_z * lambda_z;

        // d'(e,a) = d(e,a)
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t a = 0; a < v_; a++) {
                tmp3_[e*v_+a] = lz2 * dz[e+o_][a+o_];
            }
        }

        // t'(a,b,m,n) = t2(a,b,m,n) - t2(b,a,m,n) ... 2nd term here accounts for term 4 below
        C_DCOPY(o_*o_*v_*v_,t2_,1,tmp1_,1);
#pragma omp parallel for schedule(static)
        for (size_t b = 0; b < v_; b++) {
            for (size_t a = 0; a < v_; a++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {
                        tmp1_[a*o_*o_*v_+b*o_*o_+m*o_+n] -= t2_[b*o_*o_*v_+a*o_*o_+m*o_+n];
                    }
                }
            }
        }
        

        // I(e,b,m,n) = t2(a,b,m,n) d'(e,a)
        F_DGEMM('n','n',o_*o_*v_,v_,v_,1.0,tmp1_,o_*o_*v_,tmp3_,v_,0.0,tmp2_,o_*o_*v_);

        // I'(f,b) = d(i,b) u(f,i)
#pragma omp parallel for schedule(static)
        for (size_t f = 0; f < v_; f++) {
            for (size_t b = 0; b < v_; b++) {
                double dum = 0.0;
                for (size_t i = 0; i < o_; i++) {
                    dum += dz[i][b+o_] * u1_[f*o_+i];
                }
                tmp3_[f*v_+b] = dum;
            }
        }

        // I''(b,e,m,n) = I(e,b,m,n)
#pragma omp parallel for schedule(static)
        for (size_t b = 0; b < v_; b++) {
            for (size_t e = 0; e < v_; e++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {
                        tmp1_[b*o_*o_*v_+e*o_*o_+m*o_+n] = tmp2_[e*o_*o_*v_+b*o_*o_+m*o_+n];
                    }
                }
            }
        }


        // r'(f,e,m,n) += 0.5 I''(b,e,m,n) I'(f,b)
        F_DGEMM('n','n',o_*o_*v_,v_,v_,0.5,tmp1_,o_*o_*v_,tmp3_,v_,0.0,tmp2_,o_*o_*v_);

        // r(e,f,m,n) += -r'(f,e,m,n) + r'(e,f,m,n) 
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {
                        ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[f*o_*o_*v_+e*o_*o_+m*o_+n];
                        ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[e*o_*o_*v_+f*o_*o_+m*o_+n];
                    }
                }
            }
        }

        // 4: -d(e,b) d(i,a) t2(a,b,m,n) u1(f,i)
/*
        // t'(b,a,m,n) = t2(a,b,m,n)
#pragma omp parallel for schedule(static)
        for (size_t b = 0; b < v_; b++) {
            for (size_t a = 0; a < v_; a++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {
                        tmp1_[b*o_*o_*v_+a*o_*o_+m*o_+n] = t2_[a*o_*o_*v_+b*o_*o_+m*o_+n];
                    }
                }
            }
        }

        // d'(e,b) = d(e,b)
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t b = 0; b < v_; b++) {
                tmp3_[e*v_+b] = lz2 * dz[e+o_][b+o_];
            }
        }

        // I(e,a,m,n) = t'(b,a,m,n) d'(e,b)
        F_DGEMM('n','n',o_*o_*v_,v_,v_,1.0,tmp1_,o_*o_*v_,tmp3_,v_,0.0,tmp2_,o_*o_*v_);

        // I'(f,a) = d(i,a) u(f,i)
#pragma omp parallel for schedule(static)
        for (size_t f = 0; f < v_; f++) {
            for (size_t a = 0; a < v_; a++) {
                double dum = 0.0;
                for (size_t i = 0; i < o_; i++) {
                    dum += dz[i][a+o_] * u1_[f*o_+i];
                }
                tmp3_[f*v_+a] = dum;
            }
        }

        // I''(a,e,m,n) = I(e,a,m,n)
#pragma omp parallel for schedule(static)
        for (size_t a = 0; a < v_; a++) {
            for (size_t e = 0; e < v_; e++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {
                        tmp1_[a*o_*o_*v_+e*o_*o_+m*o_+n] = tmp2_[e*o_*o_*v_+a*o_*o_+m*o_+n];
                    }
                }
            }
        }


        // r'(f,e,m,n) -= 0.5 I''(a,e,m,n) I'(f,a)
        F_DGEMM('n','n',o_*o_*v_,v_,v_,-0.5,tmp1_,o_*o_*v_,tmp3_,v_,0.0,tmp2_,o_*o_*v_);

        // r(e,f,m,n) = -r'(f,e,m,n) + r'(e,f,m,n);
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {
                        ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[f*o_*o_*v_+e*o_*o_+m*o_+n];
                        ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[e*o_*o_*v_+f*o_*o_+m*o_+n];
                    }
                }
            }
        }
*/


    }

    // - d+(i,m) t2(e,f,n,i)
    // + d+(i,n) t2(e,f,m,i)

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t m = 0; m < o_; m++) {
            tmp1_[i*o_+m] = dp[i][m];
        }
    }
    // r'(e,f,n,m) =  d+(i,m) t2(e,f,n,i)
    F_DGEMM('n','n',o_,o_*v_*v_,o_,1.0,tmp1_,o_,t2_,o_,0.0,tmp2_,o_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[e*o_*o_*v_+f*o_*o_+n*o_+m];
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[e*o_*o_*v_+f*o_*o_+m*o_+n];
                }
            }
        }
    }

    // - d+(e,a) t2(a,f,m,n)
    // + d+(f,a) t2(a,e,m,n)
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t a = 0; a < v_; a++) {
            tmp1_[e*v_+a] = dp[e+o_][a+o_];
        }
    }
    // r'(e,f,n,m) =  t2(a,f,m,n) d+(e,a)
    F_DGEMM('n','n',o_*o_*v_,v_,v_,1.0,t2_,o_*o_*v_,tmp1_,v_,0.0,tmp2_,o_*o_*v_);
    C_DAXPY(o_*o_*v_*v_,-1.0,tmp2_,1,ru2_,1);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[f*o_*o_*v_+e*o_*o_+m*o_+n];
                }
            }
        }
    }

    // + d-(i,n) u1(f,i) u1(e,m) 
    // - d-(i,m) u1(f,i) u1(e,n) 
    // - d-(i,n) u1(e,i) u1(f,m)
    // + d-(i,m) u1(e,i) u1(f,n)
    //
    // or
    // 
    // P(e,f) P(m,n) d-(i,n) u1(f,i) u1(e,m)

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t n = 0; n < o_; n++) {
            tmp1_[i*o_+n] = dp[i][n];
        }
    }
    // I(f,n) = d-(i,n) u1(f,i)
    F_DGEMM('n','n',o_,v_,o_,1.0,tmp1_,o_,u1_,o_,0.0,tmp2_,o_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {

                    double dum = 0.0;
                    dum += u1_[e*o_+m] * tmp2_[f*o_+n];
                    dum -= u1_[e*o_+n] * tmp2_[f*o_+m];
                    dum -= u1_[f*o_+m] * tmp2_[e*o_+n];
                    dum += u1_[f*o_+n] * tmp2_[e*o_+m];

                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += dum;
                }
            }
        }
    }

    // + d-(e,a) u1(a,n) u1(f,m)
    // - d-(e,a) u1(a,m) u1(f,n)
    // - d-(f,a) u1(a,n) u1(e,m)
    // + d-(f,a) u1(a,m) u1(e,n)
    //
    // or
    //
    // P(e,f) P(m,e) d-(e,a) u1(a,n) u1(f,m)

#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t a = 0; a < v_; a++) {
            tmp1_[e*v_+a] = dp[e+o_][a+o_];
        }
    }
    // I(e,n) = u1(a,n) d-(e,a)
    F_DGEMM('n','n',o_,v_,v_,1.0,u1_,o_,tmp1_,v_,0.0,tmp2_,o_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {

                    double dum = 0.0;
                    dum += u1_[f*o_+m] * tmp2_[e*o_+n];
                    dum -= u1_[f*o_+n] * tmp2_[e*o_+m];
                    dum -= u1_[e*o_+m] * tmp2_[f*o_+n];
                    dum += u1_[e*o_+n] * tmp2_[f*o_+m];

                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += dum;
                }
            }
        }
    }

    // - 1.00000 d-(i,a) u1(a,i) u2(e,f,m,n)
    double Iia = 0.0;
    for (size_t a = 0; a < v_; a++) {
        for (size_t i = 0; i < o_; i++) {
            Iia += dp[i][a+o_] * u1_[a*o_+i];
        }
    }
    C_DAXPY(o_*o_*v_*v_,-Iia,u2_,1,ru2_,1);
    
    // + 2.0 d-(i,a) u1(a,n) u2(e,f,m,i)
    // - 2.0 d-(i,a) u1(a,m) u2(e,f,n,i)
    //
    // or
    //
    // P(m,n) 2.0 d-(i,a) u1(a,n) u2(e,f,m,i)

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t a = 0; a < v_; a++) {
            tmp1_[i*v_+a] = dp[i][a+o_];
        }
    }
    // I(i,n) = u1(a,n) d-(i,a)
    F_DGEMM('n','n',o_,o_,v_,1.0,u1_,o_,tmp1_,v_,0.0,tmp2_,o_);
    // r'(e,f,m,n) = 2.0 I(i,n) u2(e,f,m,i)
    F_DGEMM('n','n',o_,o_*v_*v_,o_,2.0,tmp2_,o_,u2_,o_,0.0,tmp1_,o_);
    C_DAXPY(o_*o_*v_*v_,1.0,tmp1_,1,ru2_,1);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp1_[e*o_*o_*v_+f*o_*o_+n*o_+m];
                }
            }
        }
    }

    // + 2.0 d-(i,a) u1(e,i) u2(a,f,m,n)
    // - 2.0 d-(i,a) u1(f,i) u2(a,e,m,n)
    //
    // or
    //
    // P(e,f) 2.0 d-(i,a) u1(e,i) u2(a,f,m,n)

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t a = 0; a < v_; a++) {
            tmp1_[i*v_+a] = dp[i][a+o_];
        }
    }
    // I(e,a) = d-(i,a) u1(e,i)
    F_DGEMM('n','n',v_,v_,o_,1.0,tmp1_,v_,u1_,o_,0.0,tmp2_,v_);
    // r'(e,f,m,n) = 2.0 u2(a,f,m,n) I(e,a)
    F_DGEMM('n','n',o_*o_*v_,v_,v_,2.0,u2_,o_*o_*v_,tmp2_,v_,0.0,tmp1_,o_*o_*v_);
    C_DAXPY(o_*o_*v_*v_,1.0,tmp1_,1,ru2_,1);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp1_[f*o_*o_*v_+e*o_*o_+m*o_+n];
                }
            }
        }
    }

    // + 1.00000 d-(i,a) u1(e,m) u2(a,f,n,i)
    // - 1.00000 d-(i,a) u1(e,n) u2(a,f,m,i)
    // - 1.00000 d-(i,a) u1(f,m) u2(a,e,n,i)
    // + 1.00000 d-(i,a) u1(f,n) u2(a,e,m,i)
    //
    // or
    //
    // P(m,n) P(e,f) d-(i,a) u2(a,f,n,i) u1(e,m)

#pragma omp parallel for schedule(static)
    for (size_t f = 0; f < v_; f++) {
        for (size_t n = 0; n < o_; n++) {
            double dum = 0.0;
            for (size_t i = 0; i < o_; i++) {
                for (size_t a = 0; a < v_; a++) {
                    dum += dp[i][a+o_] * u2_[a*o_*o_*v_+f*o_*o_+n*o_+i];
                }
            }
            tmp1_[f*o_+n] = dum;
        }
    }
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {

                    double dum = 0.0;
                    dum += u1_[e*o_+m] * tmp1_[f*o_+n];
                    dum -= u1_[e*o_+n] * tmp1_[f*o_+m];
                    dum -= u1_[f*o_+m] * tmp1_[e*o_+n];
                    dum += u1_[f*o_+n] * tmp1_[e*o_+m];

                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += dum;
                }
            }
        }
    }

    // - d-(i,m) u2(e,f,n,i) u0
    // + d-(i,n) u2(e,f,m,i) u0
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t m = 0; m < o_; m++) {
            tmp1_[i*o_+m] = dp[i][m];
        }
    }
    // r'(e,f,n,m) =  d-(i,m) u2(e,f,n,i) u0
    F_DGEMM('n','n',o_,o_*v_*v_,o_,u0_[0],tmp1_,o_,u2_,o_,0.0,tmp2_,o_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[e*o_*o_*v_+f*o_*o_+n*o_+m];
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[e*o_*o_*v_+f*o_*o_+m*o_+n];
                }
            }
        }
    }
    // - d-(e,a) u2(a,f,m,n) u0
    // + d-(f,a) u2(a,e,m,n) u0
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t a = 0; a < v_; a++) {
            tmp1_[e*v_+a] = dp[e+o_][a+o_];
        }
    }
    // r'(e,f,n,m) =  u2(a,f,m,n) d-(e,a) u0
    F_DGEMM('n','n',o_*o_*v_,v_,v_,u0_[0],u2_,o_*o_*v_,tmp1_,v_,0.0,tmp2_,o_*o_*v_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[e*o_*o_*v_+f*o_*o_+m*o_+n];
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[f*o_*o_*v_+e*o_*o_+m*o_+n];
                }
            }
        }
    }

    // additional term from revised pdaggerq ... missing in PRX 10, 041043 (2020)
    //    - P(m,n) 1.00000 d-(i,a) t2(e,f,i,m) u1(a,n) u0
    //   
    //  or
    //   
    //      P(m,n) 1.00000 d-(i,a) t2(e,f,m,i) u1(a,n) u0

    // I(i,n) = d-(i,a) u1(a,n) 
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t n = 0; n < o_; n++) {
            double dum = 0.0;
            for (size_t a = 0; a < v_; a++) {
                dum += dp[i][a+o_] * u1_[a*o_+n];
            }
            tmp1_[i*o_+n] = dum;
        }
    }
    // r'(e,f,m,n) = u0 I(i,n) t2(e,f,m,i)
    F_DGEMM('n','n',o_,o_*v_*v_,o_,u0_[0],tmp1_,o_,t2_,o_,0.0,tmp2_,o_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[e*o_*o_*v_+f*o_*o_+m*o_+n];
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[e*o_*o_*v_+f*o_*o_+n*o_+m];
                }
            }
        }
    }

    // + P(e,f) 1.00000 d-(i,a) t2(a,f,m,n) u1(e,i) u0

    // I(e,a) = d-(i,a) u1(e,i)
    // r'(e,f,m,n) = u0 t2(a,f,m,n) I(e,a)
#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v_; a++) {
        for (size_t e = 0; e < v_; e++) {
            double dum = 0.0;
            for (size_t i = 0; i < o_; i++) {
                dum += dp[i][a+o_] * u1_[e*o_+i];
            }
            tmp1_[e*v_+a] = dum;
        }
    }
    // r'(e,f,m,n) = u0 t2(a,f,m,n) I(e,a)
    F_DGEMM('n','n',o_*o_*v_,v_,v_,u0_[0],t2_,o_*o_*v_,tmp1_,v_,0.0,tmp2_,o_*o_*v_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[e*o_*o_*v_+f*o_*o_+m*o_+n];
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[f*o_*o_*v_+e*o_*o_+m*o_+n];
                }
            }
        }
    }

}

// 
// u0 residual
// 
// fully-contracted strings:
// - 1.00000 d+(i,i) 
// + 1.00000 F(i,a) u1(a,i) 
// + 0.25000 <i,j||a,b> u2(a,b,i,j) 
// + 1.00000 u0 w0 
//
// plus terms from coherent basis
//
// + 1.00000 b+ 
//


void PolaritonicUCCSD::residual_u0() {

    // TODO: generalize for x,y,z
    double w0 = cavity_frequency_[2];
    double r0 = u0_[0] * w0;

    // - d+(i,i) 

    double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];
    std::shared_ptr<Matrix> dip (new Matrix(Dipole_z_));
    dip->scale(coupling_factor_z);
    double ** dp = dip->pointer();

    for (size_t i = 0; i < o_; i++) {
        r0 -= dp[i][i];
    }

    // molecular hamiltonian treated in coherent state basis
    if ( !is_hubbard_ ) {
        // + 1.00000 b+
        r0 += coupling_factor_z * e_dip_z_;
    }

    if ( include_u1_ ) {

        // F(i,a) u1(a,i) 
        double **fp = F_->pointer();
        for (size_t i = 0; i < o_; i++) {
            for (size_t a = 0; a < v_; a++) {
                r0 += fp[i][a+o_] * u1_[a*o_+i];
            }
        }

        // additional term from revised pdaggerq ... missing in PRX 10, 041043 (2020)
        //     - 1.00000 d-(i,a) u1(a,i) u0
        for (size_t i = 0; i < o_; i++) {
            for (size_t a = 0; a < v_; a++) {
                r0 -= dp[i][a+o_] * u1_[a*o_+i] * u0_[0];
            }
        }

    }

    if ( include_u2_ ) {

        // + 0.25 <i,j||a,b> u2(a,b,i,j) 
        for (size_t a = 0; a < v_; a++) {
            for (size_t b = 0; b < v_; b++) {
                for (size_t i = 0; i < o_; i++) {
                    for (size_t j = 0; j < o_; j++) {
                        r0 += 0.25 * eri_ijab_[i*o_*v_*v_+j*v_*v_+a*v_+b] * u2_[a*o_*o_*v_+b*o_*o_+i*o_+j];
                    }
                }
            }
        }
    }

    ru0_[0] = r0;

}

// 
// t1 residual
// 
// fully-contracted strings:
//     + 1.00000 F(e,m) 
//     + 1.00000 F(i,a) t2(a,e,i,m) 
//     - 0.50000 <i,j||a,m> t2(a,e,i,j) 
//     + 0.50000 <i,e||a,b> t2(a,b,i,m) 

//     - 1.00000 d-(e,m) u0 
//     - 1.00000 d-(i,i) u1(e,m) 
//     + 1.00000 d-(i,m) u1(e,i) 
//     - 1.00000 d-(e,a) u1(a,m) 
//     - 1.00000 d-(i,a) u2(a,e,i,m) 
//     - 1.00000 d-(i,a) t2(a,e,i,m) u0
//
// plus terms from coherent basis
//
//     + 1.00000 u1(e,m) b- 
//


void PolaritonicUCCSD::residual_t1() {

    memset((void*)rt1_,'\0',o_*v_*sizeof(double));

    double ** fp = F_->pointer();

#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t m = 0; m < o_; m++) {
            rt1_[e*o_+m] = fp[(e+o_)][m];
        }
    }

    // F(i,a) t2(a,e,i,m) 
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t m = 0; m < o_; m++) {
            double dum = 0.0;
            for (size_t i = 0; i < o_; i++) {
                for (size_t a = 0; a < v_; a++) {
                    dum += fp[i][(a+o_)] * t2_[a*o_*o_*v_+e*o_*o_+i*o_+m];
                }
            }
            rt1_[e*o_+m] += dum;
        }
    }

// TODO: refactor o^3v integral terms
    
    // + 0.5 <i,j||m,a> t2(a,e,i,j)

    // v'(m,i,j,a) = <i,j||m,a>
#pragma omp parallel for schedule(static)
    for (size_t m = 0; m < o_; m++) {
        for (size_t i = 0; i < o_; i++) {
            for (size_t j = 0; j < o_; j++) {
                for (size_t a = 0; a < v_; a++) {
                    //tmp1_[m*o_*o_*v_+i*o_*v_+j*v_+a] = eri_iajk_[m*o_*o_*v_+a*o_*o_+i*o_+j];
                    tmp1_[m*o_*o_*v_+i*o_*v_+j*v_+a] = eri_jkia_[i*o_*o_*v_+j*o_*v_+m*v_+a];
                }
            }
        }
    }

    // t'(i,j,a,e) = t2(a,e,i,j)
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t j = 0; j < o_; j++) {
            for (size_t a = 0; a < v_; a++) {
                for (size_t e = 0; e < v_; e++) {
                    tmp2_[i*o_*v_*v_+j*v_*v_+a*v_+e] = t2_[a*o_*o_*v_+e*o_*o_+i*o_+j];
                }
            }
        }
    }

    // r(e,m) = v'(m,i,j,a) t'(i,j,a,e)
    F_DGEMM('t','t',o_,v_,o_*o_*v_,0.5,tmp1_,o_*o_*v_,tmp2_,v_,1.0,rt1_,o_);

    // - 0.5 <e,i||a,b> t2(a,b,i,m)

    // t'(i,a,b,m) = t2(a,b,i,m)
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t a = 0; a < v_; a++) {
            for (size_t b = 0; b < v_; b++) {
                for (size_t m = 0; m < o_; m++) {
                    tmp1_[i*o_*v_*v_+a*o_*v_+b*o_+m] = t2_[a*o_*o_*v_+b*o_*o_+i*o_+m];
                }   
            }   
        }   
    }


    // r(e,m) = -0.5 t'(i,a,b,m) <e,i||a,b>
    if ( is_hubbard_ ) { 
        F_DGEMM('n','n',o_,v_,o_*v_*v_,-0.5,tmp1_,o_,eri_aibc_,o_*v_*v_,1.0,rt1_,o_);
    }else {

        // four terms to avoid storing <e,i||a,b>...

        // r(e,m) = -0.5 t'(i,a,b,m) [(ea|Q)(Q|ib) - (eb|Q)(Q|ia) + d(ea)d(ib) - d(eb)d(ia) ]

        // 1: -0.5 t'(i,a,b,m) (ea|ib)
        // 2: +0.5 t'(i,a,b,m) (eb|ia)

        // t''(i,b,a,m) = t'(i,a,b,m) = t2(a,b,i,m) - t2(b,a,i,m)
        // ... account for terms 1 and 2 simultaneously
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < o_; i++) {
            for (size_t b = 0; b < v_; b++) {
                for (size_t a = 0; a < v_; a++) {
                    for (size_t m = 0; m < o_; m++) {
                        tmp2_[i*o_*v_*v_+b*o_*v_+a*o_+m] = t2_[a*o_*o_*v_+b*o_*o_+i*o_+m] - t2_[b*o_*o_*v_+a*o_*o_+i*o_+m];
                    }
                }
            }
        }
        // I(Q,a,m) = t''(i,b,a,m) (Q,i,b)
        F_DGEMM('n','n',o_*v_,nQ_,o_*v_,1.0,tmp2_,o_*v_,Qov_,o_*v_,0.0,tmp3_,o_*v_);

        // overwrite Qvv: I'(Q,a,e) = (Q|ea)
#pragma omp parallel for schedule(static)
        for (size_t q = 0; q < nQ_; q++) {
            for (size_t a = 0; a < v_; a++) {
                for (size_t e = 0; e < v_; e++) {
                    tmp1_[q*v_*v_+a*v_+e] = Qvv_[q*v_*v_+e*v_+a];
                }
            }
        }

        // r(e,m) = -0.5 I(Q,a,m) I'(Q,a,e)
        F_DGEMM('n','t',o_,v_,nQ_*v_,-0.5,tmp3_,o_,tmp1_,v_,1.0,rt1_,o_);

/*
        // 2: +0.5 t'(i,a,b,m) (eb|ia)

        // t'(i,a,b,m) = t2(a,b,i,m)
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < o_; i++) {
            for (size_t a = 0; a < v_; a++) {
                for (size_t b = 0; b < v_; b++) {
                    for (size_t m = 0; m < o_; m++) {
                        tmp3_[i*o_*v_*v_+a*o_*v_+b*o_+m] = t2_[a*o_*o_*v_+b*o_*o_+i*o_+m];
                    }
                }
            }
        }
        // I(Q,b,m) = +0.5 t'(i,a,b,m) (Q|i,a)
        F_DGEMM('n','n',o_*v_,nQ_,o_*v_,0.5,tmp3_,o_*v_,Qov_,o_*v_,0.0,tmp2_,o_*v_);

        // r(e,m) = I(Q,b,m) (Q|e,b) 
        // r(e,m) = I(Q,b,m) I'(Q,b,e) ... transposed Qvv still in tmp1
        F_DGEMM('n','t',o_,v_,nQ_*v_,1.0,tmp2_,o_,tmp1_,v_,1.0,rt1_,o_);
*/

        // 3: t'(i,a,b,m) d(e,a) d(i,b)
        double ** dz = Dipole_z_->pointer();
        double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);
        double lz2 = lambda_z * lambda_z;

        // I(a,m) = t'(i,a,b,m) d(i,b)
#pragma omp parallel for schedule(static)
       for (size_t a = 0; a < v_; a++) {
           for (size_t m = 0; m < o_; m++) {
               double dum = 0.0;
               for (size_t i = 0; i < o_; i++) {
                   for (size_t b = 0; b < v_; b++) {
                       //dum += tmp3_[i*o_*v_*v_+a*o_*v_+b*o_+m] * dz[i][b+o_];
                       dum += tmp2_[i*o_*v_*v_+b*o_*v_+a*o_+m] * dz[i][b+o_]; // account for terms 3 and 4 simultaneously
                   }
               }
               tmp1_[a*o_+m] = dum;
           }
       }

       // r(e,m) = -0.5 I(a,m) d(e,a) ... TODO could replace with DGEMM, but likely doesn't matter
#pragma omp parallel for schedule(static)
       for (size_t e = 0; e < v_; e++) {
           for (size_t m = 0; m < o_; m++) {
               double dum = 0.0;
               for (size_t a = 0; a < v_; a++) {
                   dum += tmp1_[a*o_+m] * dz[e+o_][a+o_];
               }
               rt1_[e*o_+m] -= 0.5 * lz2 * dum;
           }
       }

/*
       // 4: +0.5 t'(i,a,b,m) d(e,b) d(i,a)

       // I(b,m) = t'(i,a,b,m) d(i,a)
#pragma omp parallel for schedule(static)
       for (size_t b = 0; b < v_; b++) {
           for (size_t m = 0; m < o_; m++) {
               double dum = 0.0;
               for (size_t i = 0; i < o_; i++) {
                   for (size_t a = 0; a < v_; a++) {
                       dum += tmp3_[i*o_*v_*v_+a*o_*v_+b*o_+m] * dz[i][a+o_];
                   }
               }
               tmp1_[b*o_+m] = dum;
           }
       }

       // r(e,m) = +0.5 I(b,m) d(e,b) ... TODO could replace with DGEMM, but likely doesn't matter
#pragma omp parallel for schedule(static)
       for (size_t e = 0; e < v_; e++) {
           for (size_t m = 0; m < o_; m++) {
               double dum = 0.0;
               for (size_t b = 0; b < v_; b++) {
                   dum += tmp1_[b*o_+m] * dz[e+o_][b+o_];
               }
               rt1_[e*o_+m] += 0.5 * lz2 * dum;
           }
       }
*/

    }

    // cavity terms:

    double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];
    std::shared_ptr<Matrix> dip (new Matrix(Dipole_z_));
    dip->scale(coupling_factor_z);
    double ** dp = dip->pointer();

    if ( include_u0_ ) {

        // - d-(e,m) u0 
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t m = 0; m < o_; m++) {
                rt1_[e*o_+m] -= dp[e+o_][m] * u0_[0];
            }
        }


        // - d-(i,a) t2(a,e,i,m) u0
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t m = 0; m < o_; m++) {
                double dum = 0.0;
                for (size_t i = 0; i < o_; i++) {
                    for (size_t a = 0; a < v_; a++) {
                        dum += dp[i][a+o_] * t2_[a*o_*o_*v_+e*o_*o_+i*o_+m];
                    }
                }
                rt1_[e*o_+m] -= dum * u0_[0];
            }
        }

    }

    if ( include_u1_ ) {

        // - d-(i,i) u1(e,m) 
        double dii = 0.0;
        for (size_t i = 0; i < o_; i++) {
            dii -= dp[i][i];
        }

        // molecular hamiltonian treated in coherent state basis
        if ( !is_hubbard_ ) {
            // + 1.00000 u1(e,m) b- 
            dii += coupling_factor_z * e_dip_z_;
        }


        C_DAXPY(o_*v_,dii,u1_,1,rt1_,1);

        // + d-(i,m) u1(e,i) 
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < o_; i++) {
            for (size_t m = 0; m < o_; m++) {
                tmp1_[i*o_+m] = dp[i][m];
            }
        }
        F_DGEMM('n','n',o_,v_,o_,1.0,tmp1_,o_,u1_,o_,1.0,rt1_,o_);

        // - u1(a,m) d-(e,a) 
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t a = 0; a < v_; a++) {
                tmp1_[e*v_+a] = dp[e+o_][a+o_];
            }
        }
        F_DGEMM('n','n',o_,v_,v_,-1.0,u1_,o_,tmp1_,v_,1.0,rt1_,o_);

    }

    if ( include_u2_ ) {

        // - d-(i,a) u2(a,e,i,m) 
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t m = 0; m < o_; m++) {
                double dum = 0.0;
                for (size_t i = 0; i < o_; i++) {
                    for (size_t a = 0; a < v_; a++) {
                        dum -= dp[i][a+o_] * u2_[a*o_*o_*v_+e*o_*o_+i*o_+m];
                    }
                }
                rt1_[e*o_+m] += dum;
            }
        }
    }


}

// 
// t2 residual
// 
// <e,f||m,n>
// + 0.5 <i,j||m,n> t2(e,f,i,j)
// + 0.5 <e,f||a,b> t2(a,b,m,n)
// + P(e,f) P(m,n) <i,e||a,m> t2(a,f,i,n)
// + 0.25 <i,j||a,b> t2(a,b,m,n) t2(e,f,i,j)
// + 0.5 P(m,n) <i,j||a,b> t2(a,b,m,j) t2(e,f,n,i)  **
// - 0.5 P(e,f) <i,j||a,b> t2(a,e,i,j) t2(b,f,m,n)  **
// + P(m,n) <i,j||a,b> t2(a,e,n,j) t2(b,f,m,i)

// these terms get folded in to the ones marked ** above
// - F(i,m) t2(e,f,i,n)
// + F(i,n) t2(e,f,i,m)
// + F(e,a) t2(a,f,m,n)
// - F(f,a) t2(a,e,m,n)

// and terms introduced by cavity

// + d-(e,n) u1(f,m)
// - d-(e,m) u1(f,n)
// - d-(f,n) u1(e,m)
// + d-(f,m) u1(e,n)

// - d-(i,a) t2(e,f,i,m) u1(a,n)
// + d-(i,a) t2(e,f,i,n) u1(a,m)
// + d-(i,a) t2(a,f,m,n) u1(e,i)
// - d-(i,a) t2(a,e,m,n) u1(f,i)

// - d-(i,a) t2(a,e,i,m) u1(f,n)
// + d-(i,a) t2(a,e,i,n) u1(f,m)
// + d-(i,a) t2(a,f,i,m) u1(e,n)
// - d-(i,a) t2(a,f,i,n) u1(e,m)

// - d-(i,i) u2(e,f,m,n)
// - d-(i,n) u2(e,f,i,m)
// + d-(i,m) u2(e,f,i,n)
// - d-(e,a) u2(a,f,m,n)
// + d-(f,a) u2(a,e,m,n)

// - d-(i,n) t2(e,f,i,m) u0
// + d-(i,m) t2(e,f,i,n) u0
// - d-(e,a) t2(a,f,m,n) u0
// + d-(f,a) t2(a,e,m,n) u0

//
// plus terms from coherent basis
//
// + u2(e,f,m,n) b- 
//

void PolaritonicUCCSD::residual_t2() {

    memset((void*)rt2_,'\0',o_*o_*v_*v_*sizeof(double));

    double ** fp = F_->pointer();

    // <e,f||m,n> 
    if ( is_hubbard_ ) {

        C_DCOPY(o_*o_*v_*v_,eri_abij_,1,rt2_,1);

    }else {

        F_DGEMM('n','t',o_*v_,o_*v_,nQ_,1.0,Qvo_,o_*v_,Qvo_,o_*v_,0.0,tmp1_,o_*v_);

        double ** dz = Dipole_z_->pointer();
        double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);
        double lz2 = lambda_z * lambda_z;

#pragma omp parallel for schedule(static)
        for (size_t a = 0; a < v_; a++) {
            for (size_t b = 0; b < v_; b++) {
                for (size_t i = 0; i < o_; i++) {
                    for (size_t j = 0; j < o_; j++) {

                        rt2_[a*o_*o_*v_+b*o_*o_+i*o_+j] = tmp1_[a*o_*o_*v_+i*o_*v_+b*o_+j] - tmp1_[a*o_*o_*v_+j*o_*v_+b*o_+i];

                        double dipole_self_energy = 0.0;

                        dipole_self_energy += lz2 * dz[a+o_][i] * dz[b+o_][j];
                        dipole_self_energy -= lz2 * dz[a+o_][j] * dz[b+o_][i];

                        rt2_[a*o_*o_*v_+b*o_*o_+i*o_+j] += dipole_self_energy;
                    }
                }
            }
        }
    }

/*
    // + F(i,m) t2(e,f,n,i)
    // - F(i,n) t2(e,f,m,i)
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t m = 0; m < o_; m++) {
            tmp1_[i*o_+m] = fp[i][m];
        }
    }
    // r'(e,f,n,m) =  F(i,m) t2(e,f,n,i)
    F_DGEMM('n','n',o_,o_*v_*v_,o_,1.0,tmp1_,o_,t2_,o_,0.0,tmp2_,o_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[e*o_*o_*v_+f*o_*o_+n*o_+m];
                    rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[e*o_*o_*v_+f*o_*o_+m*o_+n];
                }
            }
        }
    }
    // + F(e,a) t2(a,f,m,n)
    // - F(f,a) t2(a,e,m,n)
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t a = 0; a < v_; a++) {
            tmp1_[e*v_+a] = fp[(e+o_)][(a+o_)];
        }
    }
    // r'(e,f,n,m) =  t2(a,f,m,n) F(e,a)
    F_DGEMM('n','n',o_*o_*v_,v_,v_,1.0,t2_,o_*o_*v_,tmp1_,v_,0.0,tmp2_,o_*o_*v_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[e*o_*o_*v_+f*o_*o_+m*o_+n];
                    rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[f*o_*o_*v_+e*o_*o_+m*o_+n];
                }
            }
        }
    }
*/
        
    // + 0.5 <i,j||m,n> t2(e,f,i,j)
    F_DGEMM('n','n',o_*o_,v_*v_,o_*o_,0.5,eri_ijkl_,o_*o_,t2_,o_*o_,1.0,rt2_,o_*o_);

    // + 0.5 <e,f||a,b> t2(a,b,m,n) = <e,f|a,b> t2(a,b,m,n)
    if ( is_hubbard_ ) {
        F_DGEMM('n','n',o_*o_,v_*v_,v_*v_,1.0,t2_,o_*o_,eri_abcd_,v_*v_,1.0,rt2_,o_*o_);
    }else {
        double_particle_ladder_diagram(t2_,rt2_);
    }

    // - P(e,f) P(m,n) <i,e||n,a> t2(a,f,m,i)

    // t'(a,i,f,m) = t2(a,f,m,i)
#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v_; a++) {
        for (size_t i = 0; i < o_; i++) {
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    tmp1_[a*o_*o_*v_+i*o_*v_+f*o_+m] = t2_[a*o_*o_*v_+f*o_*o_+m*o_+i];
                }
            }
        }
    }

    // v'(a,i,e,n) = <ie||na>
#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v_; a++) {
        for (size_t i = 0; i < o_; i++) {
            for (size_t e = 0; e < v_; e++) {
                for (size_t n = 0; n < o_; n++) {
                    tmp2_[a*o_*o_*v_+i*o_*v_+e*o_+n] = eri_iajb_[i*o_*v_*v_+e*o_*v_+n*v_+a];
                }
            }
        }
    }

    // I(e,n,f,m) = t'(a,i,f,m) v'(a,i,e,n) 
    F_DGEMM('n','t',o_*v_,o_*v_,o_*v_,1.0,tmp1_,o_*v_,tmp2_,o_*v_,0.0,tmp3_,o_*v_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {

                    double dum = 0.0;
                    dum -= tmp3_[e*o_*o_*v_+n*o_*v_+f*o_+m];
                    dum += tmp3_[e*o_*o_*v_+m*o_*v_+f*o_+n];
                    dum += tmp3_[f*o_*o_*v_+n*o_*v_+e*o_+m];
                    dum -= tmp3_[f*o_*o_*v_+m*o_*v_+e*o_+n];

                    rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += dum;
                }
            }
        }
    }

    // + 0.25 <i,j||a,b> t2(a,b,m,n) t2(e,f,i,j) 

    // I(i,j,m,n) = t2(a,b,m,n) <i,j||a,b>
    F_DGEMM('n','n',o_*o_,o_*o_,v_*v_,1.0,t2_,o_*o_,eri_ijab_,v_*v_,0.0,tmp1_,o_*o_);

    // r(e,f,m,n) = I(i,j,m,n) t'(e,f,i,j)
    F_DGEMM('n','n',o_*o_,v_*v_,o_*o_,0.25,tmp1_,o_*o_,t2_,o_*o_,1.0,rt2_,o_*o_);

    // - 0.5 P(m,n) <i,j||a,b> t2(a,b,n,j) t2(e,f,m,i) 
    //
    // and
    //
    // - F(i,m) t2(e,f,i,n)
    // + F(i,n) t2(e,f,i,m)

    // t'(n,j,a,b) = t2(a,b,n,j)
#pragma omp parallel for schedule(static)
    for (size_t n = 0; n < o_; n++) {
        for (size_t j = 0; j < o_; j++) {
            for (size_t a = 0; a < v_; a++) {
                for (size_t b = 0; b < v_; b++) {
                    tmp1_[n*o_*v_*v_+j*v_*v_+a*v_+b] = t2_[a*o_*o_*v_+b*o_*o_+n*o_+j];
                }
            }
        }
    }

    // I(i,n) = t'(n,j,a,b) <i,j||a,b> + 2 F(i,n)
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t n = 0; n < o_; n++) {
            tmp2_[i*o_+n] = 2.0 * fp[i][n];
        }
    }
    F_DGEMM('t','n',o_,o_,o_*v_*v_,1.0,tmp1_,o_*v_*v_,eri_ijab_,o_*v_*v_,1.0,tmp2_,o_);

    // r(e,f,m,n) = 0.5 I(i,n) t2(e,f,m,i)
    F_DGEMM('n','n',o_,o_*v_*v_,o_,0.5,tmp2_,o_,t2_,o_,0.0,tmp1_,o_);

    C_DAXPY(o_*o_*v_*v_,-1.0,tmp1_,1,rt2_,1);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp1_[e*o_*o_*v_+f*o_*o_+n*o_+m];
                }
            }
        }
    }

    // - 0.5 P(e,f)<i,j||a,b> t2(a,e,m,n) t2(b,f,i,j) 
    //
    // and 
    //
    // + F(e,a) t2(a,f,m,n)
    // - F(f,a) t2(a,e,m,n)
    
    // v'(a,b,i,j) = <i,j||a,b>
#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    tmp1_[a*o_*o_*v_+b*o_*o_+i*o_+j] = eri_ijab_[i*o_*v_*v_+j*v_*v_+a*v_+b];
                }
            }
        }
    }

    // t'(f,b,i,j) = t2(b,f,i,j)
#pragma omp parallel for schedule(static)
    for (size_t f = 0; f < v_; f++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    tmp2_[f*o_*o_*v_+b*o_*o_+i*o_+j] = t2_[b*o_*o_*v_+f*o_*o_+i*o_+j];
                }
            }
        }
    }

    // I(f,a) = v'(a,b,i,j) t'(f,b,i,j) + 2 F(f,a) 
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t a = 0; a < v_; a++) {
            tmp3_[e*v_+a] = 2.0 * fp[(e+o_)][(a+o_)];
        }
    }
    F_DGEMM('t','n',v_,v_,o_*o_*v_,1.0,tmp1_,o_*o_*v_,tmp2_,o_*o_*v_,1.0,tmp3_,v_);

    // r(f,e,m,n) = 0.5 t2(a,e,m,n) I(f,a)
    F_DGEMM('n','n',o_*o_*v_,v_,v_,0.5,t2_,o_*o_*v_,tmp3_,v_,0.0,tmp1_,o_*o_*v_);

    C_DAXPY(o_*o_*v_*v_,1.0,tmp1_,1,rt2_,1);
#pragma omp parallel for schedule(static)
    for (size_t f = 0; f < v_; f++) {
        for (size_t e = 0; e < v_; e++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    rt2_[f*o_*o_*v_+e*o_*o_+m*o_+n] -= tmp1_[e*o_*o_*v_+f*o_*o_+m*o_+n];
                }   
            }   
        }   
    }

    // + P(m,n) <i,j||a,b> t2(a,e,n,j) t2(b,f,m,i) 

    //I(i,b,j,a) = <i,j||a,b>
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t j = 0; j < o_; j++) {
            for (size_t a = 0; a < v_; a++) {
                for (size_t b = 0; b < v_; b++) {
                    tmp1_[i*o_*v_*v_+b*o_*v_+j*v_+a] = eri_ijab_[i*o_*v_*v_+j*v_*v_+a*v_+b];
                }
            }
        }
    }

    // t'(i,b,m,f) = t2(b,f,m,i)
#pragma omp parallel for schedule(static)
    for (size_t b = 0; b < v_; b++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t i = 0; i < o_; i++) {
                    tmp2_[i*o_*v_*v_+b*o_*v_+m*v_+f] = t2_[b*o_*o_*v_+f*o_*o_+m*o_+i];
                }
            }
        }
    }

    //I''(m,f,j,a) = I(i,b,j,a) t'(i,b,m,f) 
    F_DGEMM('n','t',o_*v_,o_*v_,o_*v_,1.0,tmp1_,o_*v_,tmp2_,o_*v_,0.0,tmp3_,o_*v_);

    // t''(n,e,j,a) = t2(a,e,n,j) 
#pragma omp parallel for schedule(static)
    for (size_t n = 0; n < o_; n++) {
        for (size_t e = 0; e < v_; e++) {
            for (size_t j = 0; j < o_; j++) {
                for (size_t a = 0; a < v_; a++) {
                    tmp1_[e*o_*o_*v_+n*o_*v_+j*v_+a] = t2_[a*o_*o_*v_+e*o_*o_+n*o_+j];
                }
            }
        }
    }

    //r'(e,n,m,f) = I''(m,f,j,a) t''(n,e,j,a)
    F_DGEMM('t','n',o_*v_,o_*v_,o_*v_,1.0,tmp3_,o_*v_,tmp1_,o_*v_,0.0,tmp2_,o_*v_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[e*o_*o_*v_+n*o_*v_+m*v_+f];
                    rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[e*o_*o_*v_+m*o_*v_+n*v_+f];
                }
            }
        }
    }

    // now, terms introduced by cavity

    if ( include_u1_ ) {

        double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];
        std::shared_ptr<Matrix> dip (new Matrix(Dipole_z_));
        dip->scale(coupling_factor_z);
        double ** dp = dip->pointer();

        // + d-(e,n) u1(f,m)
        // - d-(e,m) u1(f,n)
        // - d-(f,n) u1(e,m)
        // + d-(f,m) u1(e,n)
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {
                        rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += dp[e+o_][n] * u1_[f*o_+m];
                        rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= dp[e+o_][m] * u1_[f*o_+n];
                        rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= dp[f+o_][n] * u1_[e*o_+m];
                        rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += dp[f+o_][m] * u1_[e*o_+n];
                    }
                }
            }
        }

        // - d-(i,a) t2(a,e,i,m) u1(f,n)
        // + d-(i,a) t2(a,e,i,n) u1(f,m)
        // + d-(i,a) t2(a,f,i,m) u1(e,n)
        // - d-(i,a) t2(a,f,i,n) u1(e,m)
        //
        // or 
        //
        // -P(e,f) P(m,n) d-(i,a) t2(a,e,i,m) u1(f,n)
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t m = 0; m < o_; m++) {
                double dum = 0.0;
                for (size_t a = 0; a < v_; a++) {
                    for (size_t i = 0; i < o_; i++) {
                        dum += dp[i][a+o_] * t2_[a*o_*o_*v_+e*o_*o_+i*o_+m];
                    }
                }
                tmp1_[e*o_+m] = dum;
            }
        }
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {
                        rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp1_[e*o_+m] * u1_[f*o_+n];
                        rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp1_[e*o_+n] * u1_[f*o_+m];
                        rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp1_[f*o_+m] * u1_[e*o_+n];
                        rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp1_[f*o_+n] * u1_[e*o_+m];
                    }
                }
            }
        }

        // - d-(i,a) t2(e,f,i,m) u1(a,n)
        // + d-(i,a) t2(e,f,i,n) u1(a,m)
        //
        // or
        //
        // P(m,n) d-(i,a) t2(e,f,m,i) u1(a,n)

        // I(i,n) = u1(a,n) d-(i,a)
        for (size_t i = 0; i < o_; i++) {
            for (size_t a = 0; a < v_; a++) {
                tmp1_[i*v_+a] = dp[i][a+o_];
            }
        }
        F_DGEMM('n','n',o_,o_,v_,1.0,u1_,o_,tmp1_,v_,0.0,tmp2_,o_);

        // r'(e,f,m,n) = I(i,n) t2(e,f,m,i)
        F_DGEMM('n','n',o_,o_*v_*v_,o_,1.0,tmp2_,o_,t2_,o_,0.0,tmp1_,o_);

        C_DAXPY(o_*o_*v_*v_,1.0,tmp1_,1,rt2_,1);
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {
                        rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp1_[e*o_*o_*v_+f*o_*o_+n*o_+m];
                    }
                }
            }
        }

        // + d-(i,a) t2(a,f,m,n) u1(e,i)
        // - d-(i,a) t2(a,e,m,n) u1(f,i)
        //
        // or
        //
        // P(e,f) d-(i,a) t2(a,f,m,n) u1(e,i)

        // I(e,a) =  d-(i,a) u1(e,i)
        for (size_t i = 0; i < o_; i++) {
            for (size_t a = 0; a < v_; a++) {
                tmp1_[i*v_+a] = dp[i][a+o_];
            }
        }
        F_DGEMM('n','n',v_,v_,o_,1.0,tmp1_,v_,u1_,o_,0.0,tmp2_,v_);

        // r'(e,f,m,n) = t2(a,f,m,n) I(e,a)
        F_DGEMM('n','n',o_*o_*v_,v_,v_,1.0,t2_,o_*o_*v_,tmp2_,v_,0.0,tmp1_,o_*o_*v_);
        C_DAXPY(o_*o_*v_*v_,1.0,tmp1_,1,rt2_,1);
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {
                        rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp1_[f*o_*o_*v_+e*o_*o_+m*o_+n];
                    }
                }
            }
        }
        
    }

    if ( include_u2_ ) {

        double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];
        std::shared_ptr<Matrix> dip (new Matrix(Dipole_z_));
        dip->scale(coupling_factor_z);
        double ** dp = dip->pointer();

        // - d-(i,i) u2(e,f,m,n)
        double dii = 0.0;
        for (size_t i = 0; i < o_; i++) {
            dii -= dp[i][i];
        }
        // molecular hamiltonian treated in coherent state basis
        if ( !is_hubbard_ ) {
            // + u2(e,f,m,n) b- 
            dii += coupling_factor_z * e_dip_z_;
        }
        C_DAXPY(o_*o_*v_*v_,dii,u2_,1,rt2_,1);
        
        // + P(m,n) d-(i,n) u2(e,f,m,i)
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < o_; i++) {
            for (size_t n = 0; n < o_; n++) {
                tmp1_[i*o_+n] = dp[i][n];
            }
        }
        F_DGEMM('n','n',o_,o_*v_*v_,o_,1.0,tmp1_,o_,u2_,o_,0.0,tmp2_,o_);
        C_DAXPY(o_*o_*v_*v_,1.0,tmp2_,1,rt2_,1);
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {
                        rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[e*o_*o_*v_+f*o_*o_+n*o_+m];
                    }
                }
            }
        }

        // - P(e,f) u2(a,f,m,n) d-(e,a) 
        for (size_t e = 0; e < v_; e++) {
            for (size_t a = 0; a < v_; a++) {
                tmp1_[e*v_+a] = dp[e+o_][a+o_];
            }
        }
        F_DGEMM('n','n',o_*o_*v_,v_,v_,1.0,u2_,o_*o_*v_,tmp1_,v_,0.0,tmp2_,o_*o_*v_);
        C_DAXPY(o_*o_*v_*v_,-1.0,tmp2_,1,rt2_,1);
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {
                        rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[f*o_*o_*v_+e*o_*o_+m*o_+n];
                    }
                }
            }
        }

    }

    if ( include_u0_ ) {

        double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];
        std::shared_ptr<Matrix> dip (new Matrix(Dipole_z_));
        dip->scale(coupling_factor_z);
        double ** dp = dip->pointer();

        // + u0 P(m,n) d-(i,n) t2(e,f,m,i)
        for (size_t i = 0; i < o_; i++) {
            for (size_t n = 0; n < o_; n++) {
                tmp1_[i*o_+n] = dp[i][n];
            }
        }
        F_DGEMM('n','n',o_,o_*v_*v_,o_,u0_[0],tmp1_,o_,t2_,o_,0.0,tmp2_,o_);
        C_DAXPY(o_*o_*v_*v_,1.0,tmp2_,1,rt2_,1);
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {
                        rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[e*o_*o_*v_+f*o_*o_+n*o_+m];
                    }
                }
            }
        }

        // - u0 P(e,f) t2(a,f,m,n) d-(e,a)
        for (size_t e = 0; e < v_; e++) {
            for (size_t a = 0; a < v_; a++) {
                tmp1_[e*v_+a] = dp[e+o_][a+o_];
            }
        }
        F_DGEMM('n','n',o_*o_*v_,v_,v_,u0_[0],t2_,o_*o_*v_,tmp1_,v_,0.0,tmp2_,o_*o_*v_);
        C_DAXPY(o_*o_*v_*v_,-1.0,tmp2_,1,rt2_,1);
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v_; e++) {
            for (size_t f = 0; f < v_; f++) {
                for (size_t m = 0; m < o_; m++) {
                    for (size_t n = 0; n < o_; n++) {
                        rt2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[f*o_*o_*v_+e*o_*o_+m*o_+n];
                    }
                }
            }
        }

    }

}

double PolaritonicUCCSD::update_amplitudes() { 
    // dt = -residual / eps

    // t2
    for (size_t a = 0; a < v_; a++) {
        double da = epsilon_[a+o_];
        for (size_t b = 0; b < v_; b++) {
            double dab = da + epsilon_[b+o_];
            for (size_t i = 0; i < o_; i++) {
                double dabi = dab - epsilon_[i];
                for (size_t j = 0; j < o_; j++) {
                    double dabij = dabi - epsilon_[j];
                    size_t abij = a*o_*o_*v_+b*o_*o_+i*o_+j;
                    rt2_[abij] /= -dabij;
                }
            }
        }
    }

    // t1
    for (size_t a = 0; a < v_; a++) {
        double da = epsilon_[a+o_];
        for (size_t i = 0; i < o_; i++) {
            double dai = da - epsilon_[i];
            size_t ai = a*o_+i;
            rt1_[ai] /= -dai;
        }
    }

    if ( include_u2_ ) {
        // u2
        for (size_t a = 0; a < v_; a++) {
            double da = epsilon_[a+o_];
            for (size_t b = 0; b < v_; b++) {
                double dab = da + epsilon_[b+o_];
                for (size_t i = 0; i < o_; i++) {
                    double dabi = dab - epsilon_[i];
                    for (size_t j = 0; j < o_; j++) {
                        double dabij = dabi - epsilon_[j];
                        size_t abij = a*o_*o_*v_+b*o_*o_+i*o_+j;
                        ru2_[abij] /= -dabij;
                    }
                }
            }
        }
    }
    if ( include_u1_ ) {
        // u1
        double w0 = cavity_frequency_[2];
        for (size_t a = 0; a < v_; a++) {
            double da = epsilon_[a+o_];
            for (size_t i = 0; i < o_; i++) {
                double dai = da - epsilon_[i];
                size_t ai = a*o_+i;
                ru1_[ai] /= -(dai+w0);
            }
        }
    }
    if ( include_u0_ ) {
        // 0   = u0 w0 + ru0
        // u0(k+1) = u0(k) + du0(k)
        // u0(k+1) = -ru0 / w0
        // du0 = u0(k+1) - u(k) = -ru0 / w0 - u0(k)

        // TODO: generalize for x,y,z
        double w0 = cavity_frequency_[2];
        //ru0_[0] = -ru0_[0] / w0 - u0_[0];
        ru0_[0] /= -w0;

        // catch zero frequency:
        if ( fabs(w0) < 1e-12 ) {
           ru0_[0] = 0.0;
        }
    }

    // diis 
    C_DAXPY(ccamps_dim_,1.0,residual_,1,ccamps_,1);
    diis->WriteVector(ccamps_);
    diis->WriteErrorVector(residual_);
    diis->Extrapolate(ccamps_);

    return C_DNRM2(ccamps_dim_,residual_,1);

}

// fully-contracted strings:
// - 1.00000 d-(i,i) u0 
// - 1.00000 d-(i,a) u1(a,i) 
// + 0.25000 <i,j||a,b> t2(a,b,i,j) 
//
// plus t1 parts handled elsewhere:
//
// + 1.00000 h(i,i) 
// + 0.50000 <i,j||i,j> 
//
// plus terms from coherent basis
//
// + 1.00000 u0 b- 

double PolaritonicUCCSD::correlation_energy() {

    double ec = 0.0;

    // + 0.25000 <i,j||a,b> t2(a,b,i,j)
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    ec += 0.25 * eri_ijab_[i*o_*v_*v_+j*v_*v_+a*v_+b] * t2_[a*o_*o_*v_+b*o_*o_+i*o_+j];
                }
            }
        }
    }

    if ( include_u0_ ) {

        double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];
        std::shared_ptr<Matrix> dip (new Matrix(Dipole_z_));
        dip->scale(coupling_factor_z);
        double ** dp = dip->pointer();

        // - 1.00000 d-(i,i) u0 
        for (size_t i = 0; i < o_; i++) {
            ec -= dp[i][i] * u0_[0];
        }

        // molecular hamiltonian treated in coherent state basis
        if ( !is_hubbard_ ) {
            // + 1.00000 u0 b- 
            ec += u0_[0] * coupling_factor_z * e_dip_z_;
        }

    }

    if ( include_u1_ ) {

        double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];
        std::shared_ptr<Matrix> dip (new Matrix(Dipole_z_));
        dip->scale(coupling_factor_z);
        double ** dp = dip->pointer();

        // - 1.00000 d-(i,a) u1(a,i) 
        for (size_t i = 0; i < o_; i++) {
            for (size_t a = 0; a < v_; a++) {
                ec -= dp[i][a+o_] * u1_[a*o_+i];
            }
        }

    }

    return ec;
}

/**
 *  double particle ladder diagram
 */
void PolaritonicUCCSD::double_particle_ladder_diagram(double * t2, double * r2) {

    size_t oov = o_ * o_ * v_;
    size_t oo = o_ * o_;
    size_t otri = o_ * (o_ + 1L) / 2L;
    size_t vtri = v_ * (v_ + 1L) / 2L;

    double * Abij = (double *)malloc(otri * v_ * sizeof(double));
    double * Sbij = (double *)malloc(otri * v_ * sizeof(double));

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o_; i++) {
        for (size_t j = i; j < o_; j++) {
            size_t ij = INDEX(i, j);
            for (size_t a = 0; a < v_; a++) {
                for (size_t b = a; b < v_; b++) {
                    tmp2_[INDEX(a, b) * otri + ij] =
                        (t2[a * oov + b * oo + i * o_ + j] + t2[b * oov + a * oo + i * o_ + j]);
                    tmp2_[INDEX(a, b) * otri + ij + vtri * otri] =
                        (t2[a * oov + b * oo + i * o_ + j] - t2[b * oov + a * oo + i * o_ + j]);
                }
                tmp2_[INDEX(a, a) * otri + ij] = t2[a * oov + a * oo + i * o_ + j];
            }
        }
    }

    int nthreads = Process::environment.get_n_threads();

    double * Vcdb = tmp1_;
    double * Vm   = tmp1_ + v_ * v_ * v_;
    double * Vp   = Vm;

    // read transpose of Qvv from disk
    auto psio = std::make_shared<PSIO>();
    psio->open(PSIF_DCC_QSO, PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_QSO, "Qvv transpose", (char *)&Qvv_[0], nQ_ * v_ * v_ * sizeof(double));
    psio->close(PSIF_DCC_QSO, 1);

    for (size_t a = 0; a < v_; a++) {

        int nb = v_ - a;
        F_DGEMM('t', 'n', v_, v_ * nb, nQ_, 1.0, Qvv_ + a * v_ * nQ_, nQ_, Qvv_ + a * v_ * nQ_, nQ_, 0.0, Vcdb, v_);

#pragma omp parallel for schedule(static)
        for (size_t b = a; b < v_; b++) {
            size_t cd = 0;
            size_t ind1 = (b - a) * vtri;
            size_t ind2 = (b - a) * v_ * v_;
            size_t v1, v2;
            for (size_t c = 0; c < v_; c++) {
                for (size_t d = 0; d <= c; d++) {
                    Vp[ind1 + cd] = Vcdb[ind2 + d * v_ + c] + Vcdb[ind2 + c * v_ + d];
                    cd++;
                }
            }
        }

        F_DGEMM('n', 'n', otri, nb, vtri, 0.5, tmp2_, otri, Vp, vtri, 0.0, Abij, otri);
#pragma omp parallel for schedule(static)
        for (size_t b = a; b < v_; b++) {
            size_t cd = 0;
            size_t ind1 = (b - a) * vtri;
            size_t ind2 = (b - a) * v_ * v_;
            size_t v1, v2;
            for (size_t c = 0; c < v_; c++) {
                for (size_t d = 0; d <= c; d++) {
                    Vm[ind1 + cd] = Vcdb[ind2 + d * v_ + c] - Vcdb[ind2 + c * v_ + d];
                    cd++;
                }
            }
        }
        F_DGEMM('n', 'n', otri, nb, vtri, 0.5, tmp2_ + otri * vtri, otri, Vm, vtri, 0.0, Sbij, otri);

        // contribute to residual
#pragma omp parallel for schedule(static)
        for (size_t b = a; b < v_; b++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    int sg = (i > j) ? 1 : -1;
                    r2[a * oo * v_ + b * oo + i * o_ + j] +=
                        Abij[(b - a) * otri + INDEX(i, j)] + sg * Sbij[(b - a) * otri + INDEX(i, j)];
                    if (a != b) {
                        r2[b * oov + a * oo + i * o_ + j] +=
                            Abij[(b - a) * otri + INDEX(i, j)] - sg * Sbij[(b - a) * otri + INDEX(i, j)];
                    }
                }
            }
        }

    }

    // read Qvv back in
    psio->open(PSIF_DCC_QSO, PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_QSO, "Qvv", (char *)&Qvv_[0], nQ_ * v_ * v_ * sizeof(double));
    psio->close(PSIF_DCC_QSO, 1);

    free(Abij);
    free(Sbij);

    // now, cavity contribution to this term: t2(c,d,i,j) d(a,c) d(b,d) lambda^2
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);
    double ** dz = Dipole_z_->pointer();

#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v_; a++) {
        for (size_t c = 0; c < v_; c++) {
            tmp1_[a*v_+c] = dz[a+o_][c+o_];
        }
    }
    // I(a,d,i,j) = t2(c,d,i,j) d(a,c) lambda_z^2
    F_DGEMM('n','n',o_*o_*v_,v_,v_,lambda_z*lambda_z,t2,o_*o_*v_,tmp1_,v_,0.0,tmp2_,o_*o_*v_);

    // I'(d,a,i,j) = I(a,d,i,j)
#pragma omp parallel for schedule(static)
    for (size_t d = 0; d < v_; d++) {
        for (size_t a = 0; a < v_; a++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    tmp3_[d*o_*o_*v_+a*o_*o_+i*o_+j] = tmp2_[a*o_*o_*v_+d*o_*o_+i*o_+j];
                }
            }
        }
    }

    // r'(b,a,i,j) = I'(d,a,i,j) d(b,d)
    F_DGEMM('n','n',o_*o_*v_,v_,v_,1.0,tmp3_,o_*o_*v_,tmp1_,v_,0.0,tmp2_,o_*o_*v_);

#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    r2[a*o_*o_*v_+b*o_*o_+i*o_+j] += tmp2_[b*o_*o_*v_+a*o_*o_+i*o_+j];
                    //double dumx = 0.0;
                    //for (size_t c = 0; c < v_; c++) {
                    //    for (size_t d = 0; d < v_; d++) {
                    //        dumz += t2[c*o_*o_*v_+d*o_*o_+i*o_+j] * dz[a+o_][c+o_] * dz[b+o_][d+o_];
                    //    }
                    //}
                    //r2[a*o_*o_*v_+b*o_*o_+i*o_+j] += lambda_z * lambda_z * dumz;
                }
            }
        }
    }

}

} // End namespaces
