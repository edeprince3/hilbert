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
    outfile->Printf("  ==> Hamiltonian type <==\n");
    outfile->Printf("\n");
    outfile->Printf("        %s\n",is_hubbard_ ? "Hubbard" : "molecular");
    outfile->Printf("\n");

    same_a_b_orbs_ = false;
    same_a_b_dens_ = false;

    // include amplitudes for photon transitions?
    include_u0_ = false;
    include_u1_ = false;
    include_u2_ = false;

    // molecular hamiltonian
    if ( !is_hubbard_ ) {

        initialize_with_molecular_hamiltonian();

    }else {

        initialize_with_hubbard_hamiltonian();

    }

    // allocate memory for amplitudes, residual, and temporary buffers

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

    // temporary storage ... reduce later

    tmp1_ = (double*)malloc(o_*o_*v_*v_*sizeof(double));
    tmp2_ = (double*)malloc(o_*o_*v_*v_*sizeof(double));
    tmp3_ = (double*)malloc(o_*o_*v_*v_*sizeof(double));
    memset((void*)tmp1_,'\0',o_*o_*v_*v_*sizeof(double));
    memset((void*)tmp2_,'\0',o_*o_*v_*v_*sizeof(double));
    memset((void*)tmp3_,'\0',o_*o_*v_*v_*sizeof(double));

    // initialize diis solver
    diis = (std::shared_ptr<DIIS>)(new DIIS(ccamps_dim_));

}

void PolaritonicUCCSD::initialize_with_hubbard_hamiltonian() {

    throw PsiException("hubbard hamiltonian not yet implemented for UCCSD",__FILE__,__LINE__);

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

    // alpha + beta MO transformation matrix
    C_ = (std::shared_ptr<Matrix>)(new Matrix(2*nso_,2*nmo_));
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
    double ** oe_cavity_terms_p = oe_cavity_terms_->pointer();

    // electron-nucleus contribution to dipole self energy
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

    double ** eps = F_->pointer();
    for (size_t i = 0; i < 2*nso_; i++) {
        epsilon_[i] = eps[i][i];
    }

    // generate and write SO-basis three-index integrals to disk
    write_three_index_ints();

    // construct four-index integrals
    build_mo_eris();
}

// TODO: this is so wasteful. need to
// 1. fold t1 into 3-index integrals (done!)
// 2. don't build full four-index tensor
// 3. don't build (ac|bd) or any ov^4 quantities

// write three-index integrals to disk (file PSIF_DCC_QSO)
void PolaritonicUCCSD::write_three_index_ints() {

    // get primary basis:
    std::shared_ptr<BasisSet> primary = reference_wavefunction_->get_basisset("ORBITAL");

    // TODO: use DF_BASIS_CC

    // get auxiliary basis:
    std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");

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

    double ec = 0.0;

    if ( is_hubbard_ ) {

        ec = t1_transformation_hubbard_hamiltonian();

    }else {

        ec = t1_transformation_molecular_hamiltonian();

    }

    return ec;

}

double PolaritonicUCCSD::t1_transformation_hubbard_hamiltonian() {

    throw PsiException("hubbard hamiltonian not yet implemented for UCCSD",__FILE__,__LINE__);

}

double PolaritonicUCCSD::t1_transformation_molecular_hamiltonian() {

    size_t n  = 2L*(size_t)nmo_;
    size_t ns = 2L*(size_t)nso_;

    // read Qso from disk

    double * tmp = (double*)malloc(nQ_*ns*ns*sizeof(double));
    memset((void*)tmp,'\0',nQ_*n*n*sizeof(double));

    double * Qmo = (double*)malloc(nQ_*n*ns*sizeof(double));
    memset((void*)Qmo,'\0',nQ_*n*ns*sizeof(double));

    auto psio = std::make_shared<PSIO>();

    psio->open(PSIF_DCC_QSO, PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_QSO, "Qso CC", (char*)Qmo, nQ_ * nso_ * nso_ * sizeof(double));
    psio->close(PSIF_DCC_QSO, 1);

    for (size_t Q = 0; Q < nQ_; Q++) {
        for (size_t mu = 0; mu < nso_; mu++) {
            for (size_t nu = 0; nu < nso_; nu++) {
                tmp[Q*ns*ns+(mu     )*ns+(nu     )] = Qmo[Q*nso_*nso_+mu*nso_+nu];
                tmp[Q*ns*ns+(mu+nso_)*ns+(nu+nso_)] = Qmo[Q*nso_*nso_+mu*nso_+nu];
            }
        }
    }

    // AO->MO->t1 transformation

    std::shared_ptr<Matrix> CR (new Matrix(C_));
    std::shared_ptr<Matrix> CL (new Matrix(C_));

    double ** CR_p = CR->pointer();
    double ** CL_p = CL->pointer();
    double ** C_p  = C_->pointer();

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

    // I(Q,mu,p) = C(nu,p) Qso(Q,mu,nu)
    F_DGEMM('n','n',n,ns*nQ_,ns,1.0,&(CL->pointer()[0][0]),n,tmp,ns,0.0,Qmo,n);
    for (size_t Q = 0; Q < nQ_; Q++) {
        for (size_t p = 0; p < n; p++) {
            for (size_t mu = 0; mu < ns; mu++) {
                tmp[Q*n*ns+p*ns+mu] = Qmo[Q*n*ns+mu*n+p];
            }
        }
    }
    // Qmo(Q,p,q) = C(mu,q) I(Q,p,mu)
    F_DGEMM('n','n',n,n*nQ_,ns,1.0,&(CR->pointer()[0][0]),n,tmp,ns,0.0,Qmo,n);

    free(tmp);

    double * eri = (double*)malloc(n*n*n*n*sizeof(double));
    memset((void*)eri,'\0',n*n*n*n*sizeof(double));

    // (pq|rs) = Qmo(Q,rs) Qmo(Q,pq)
    F_DGEMM('n','t',n*n,n*n,nQ_,1.0,Qmo,n*n,Qmo,n*n,0.0,eri,n*n);

    // add two-electron part of dipole self energy to eris
    double * oei_tmp = (double*)malloc(ns*ns*sizeof(double));

    // transform dipole integrals
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
    std::shared_ptr<Matrix> dipole_tmp (new Matrix(ns,ns));

    // t1 transformation (dipole x)
    for (size_t mu = 0; mu < ns; mu++) {
        for (size_t p = 0; p < n; p++) {
            double dum = 0.0;
            for (size_t nu = 0; nu < ns; nu++) {
                dum += CL_p[nu][p] * dx[nu][mu];
            }
            oei_tmp[p * ns + mu] = dum;
        }
    }
    for (size_t p = 0; p < n; p++) {
        for (size_t q = 0; q < n; q++) {
            double dum = 0.0;
            for (size_t nu = 0; nu < ns; nu++) {
                dum += CR_p[nu][q] * oei_tmp[p * ns + nu];
            }
            dx[p][q] = dum;
        }
    }
    // t1 transformation (dipole y)
    for (size_t mu = 0; mu < ns; mu++) {
        for (size_t p = 0; p < n; p++) {
            double dum = 0.0;
            for (size_t nu = 0; nu < ns; nu++) {
                dum += CL_p[nu][p] * dy[nu][mu];
            }
            oei_tmp[p * ns + mu] = dum;
        }
    }
    for (size_t p = 0; p < n; p++) {
        for (size_t q = 0; q < n; q++) {
            double dum = 0.0;
            for (size_t nu = 0; nu < ns; nu++) {
                dum += CR_p[nu][q] * oei_tmp[p * ns + nu];
            }
            dy[p][q] = dum;
        }
    }
    // t1 transformation (dipole z)
    for (size_t mu = 0; mu < ns; mu++) {
        for (size_t p = 0; p < n; p++) {
            double dum = 0.0;
            for (size_t nu = 0; nu < ns; nu++) {
                dum += CL_p[nu][p] * dz[nu][mu];
            }
            oei_tmp[p * ns + mu] = dum;
        }
    }
    for (size_t p = 0; p < n; p++) {
        for (size_t q = 0; q < n; q++) {
            double dum = 0.0;
            for (size_t nu = 0; nu < ns; nu++) {
                dum += CR_p[nu][q] * oei_tmp[p * ns + nu];
            }
            dz[p][q] = dum;
        }
    }

    double lambda_x = cavity_coupling_strength_[0] * sqrt(2.0 * cavity_frequency_[0]);
    double lambda_y = cavity_coupling_strength_[1] * sqrt(2.0 * cavity_frequency_[1]);
    double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);

    // add two-electron part of dipole self energy to eris
    for (size_t p = 0; p < n; p++) {
        for (size_t q = 0; q < n; q++) {
            for (size_t r = 0; r < n; r++) {
                for (size_t s = 0; s < n; s++) {
                    eri[p*n*n*n+q*n*n+r*n+s] += lambda_x * lambda_x * dx[p][q] * dx[r][s];
                    eri[p*n*n*n+q*n*n+r*n+s] += lambda_y * lambda_y * dy[p][q] * dy[r][s];
                    eri[p*n*n*n+q*n*n+r*n+s] += lambda_z * lambda_z * dz[p][q] * dz[r][s];
                }
            }
        }
    }

    // now that dipole integrals have been used for self-energy term, scale for coupling terms
    double coupling_factor_x = cavity_frequency_[0] * cavity_coupling_strength_[0];
    double coupling_factor_y = cavity_frequency_[1] * cavity_coupling_strength_[1];
    double coupling_factor_z = cavity_frequency_[2] * cavity_coupling_strength_[2];
    Dipole_x_->scale(coupling_factor_x);
    Dipole_y_->scale(coupling_factor_y);
    Dipole_z_->scale(coupling_factor_z);

    // unpack different classes of eris

    // <ij||kl>
    memset((void*)eri_ijkl_,'\0',o_*o_*o_*o_*sizeof(double));
    for (size_t i = 0; i < o_; i++) {
        for (size_t j = 0; j < o_; j++) {
            for (size_t k = 0; k < o_; k++) {
                for (size_t l = 0; l < o_; l++) {
                    size_t ikjl = i*n*n*n+k*n*n+j*n+l;
                    size_t iljk = i*n*n*n+l*n*n+j*n+k;
                    eri_ijkl_[i*o_*o_*o_+j*o_*o_+k*o_+l]  = eri[ikjl] - eri[iljk];
                    
                }
            }
        }
    }

    // <ab||cd>
    memset((void*)eri_abcd_,'\0',v_*v_*v_*v_*sizeof(double));
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t c = 0; c < v_; c++) {
                for (size_t d = 0; d < v_; d++) {
                    size_t acbd = (a+o_)*n*n*n+(c+o_)*n*n+(b+o_)*n+(d+o_);
                    size_t adbc = (a+o_)*n*n*n+(d+o_)*n*n+(b+o_)*n+(c+o_);
                    eri_abcd_[a*v_*v_*v_+b*v_*v_+c*v_+d] = eri[acbd] - eri[adbc];
                }
            }
        }
    }

    // <ij||ab>
    memset((void*)eri_ijab_,'\0',o_*o_*v_*v_*sizeof(double));
    for (size_t i = 0; i < o_; i++) {
        for (size_t j = 0; j < o_; j++) {
            for (size_t a = 0; a < v_; a++) {
                for (size_t b = 0; b < v_; b++) {
                    size_t iajb = i*n*n*n+(a+o_)*n*n+j*n+(b+o_);
                    size_t ibja = i*n*n*n+(b+o_)*n*n+j*n+(a+o_);
                    eri_ijab_[i*o_*v_*v_+j*v_*v_+a*v_+b] = eri[iajb] - eri[ibja];
                }
            }
        }
    }

    // <ab||ij>
    memset((void*)eri_abij_,'\0',o_*o_*v_*v_*sizeof(double));
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    size_t aibj = (a+o_)*n*n*n+i*n*n+(b+o_)*n+j;
                    size_t ajbi = (a+o_)*n*n*n+j*n*n+(b+o_)*n+i;
                    eri_abij_[a*o_*o_*v_+b*o_*o_+i*o_+j] = eri[aibj] - eri[ajbi];
                }
            }
        }
    }

    // <ia||jb>
    memset((void*)eri_iajb_,'\0',o_*o_*v_*v_*sizeof(double));
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    size_t ijab = i*n*n*n+j*n*n+(a+o_)*n+(b+o_);
                    size_t ibaj = i*n*n*n+(b+o_)*n*n+(a+o_)*n+j;
                    eri_iajb_[i*o_*v_*v_+a*o_*v_+j*v_+b] = eri[ijab] - eri[ibaj];
                }
            }
        }
    }

    // <ia||jk>
    memset((void*)eri_iajk_,'\0',o_*o_*o_*v_*sizeof(double));
    for (size_t i = 0; i < o_; i++) {
        for (size_t a = 0; a < v_; a++) {
            for (size_t j = 0; j < o_; j++) {
                for (size_t k = 0; k < o_; k++) {
                    size_t ijak = i*n*n*n+j*n*n+(a+o_)*n+k;
                    size_t ikaj = i*n*n*n+k*n*n+(a+o_)*n+j;
                    eri_iajk_[i*o_*o_*v_+a*o_*o_+j*o_+k] = eri[ijak] - eri[ikaj];
                }
            }
        }
    }

    // <jk||ia>
    memset((void*)eri_jkia_,'\0',o_*o_*o_*v_*sizeof(double));
    for (size_t j = 0; j < o_; j++) {
        for (size_t k = 0; k < o_; k++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t a = 0; a < v_; a++) {
                    size_t jika = j*n*n*n+i*n*n+k*n+(a+o_);
                    size_t jaki = j*n*n*n+(a+o_)*n*n+k*n+i;
                    eri_jkia_[j*o_*o_*v_+k*o_*v_+i*v_+a] = eri[jika] - eri[jaki];
                }
            }
        }
    }

    // <ai||bc>
    memset((void*)eri_aibc_,'\0',o_*v_*v_*v_*sizeof(double));
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t c = 0; c < v_; c++) {
                for (size_t i = 0; i < o_; i++) {
                    size_t abic = (a+o_)*n*n*n+(b+o_)*n*n+i*n+(c+o_);
                    size_t acib = (a+o_)*n*n*n+(c+o_)*n*n+i*n+(b+o_);
                    eri_aibc_[a*o_*v_*v_+i*v_*v_+b*v_+c] = eri[abic] - eri[acib];
                }
            }
        }
    }

    // <ab||ic>
    memset((void*)eri_abic_,'\0',o_*v_*v_*v_*sizeof(double));
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t c = 0; c < v_; c++) {
                for (size_t i = 0; i < o_; i++) {
                    size_t aibc = (a+o_)*n*n*n+i*n*n+(b+o_)*n+(c+o_);
                    size_t acbi = (a+o_)*n*n*n+(c+o_)*n*n+(b+o_)*n+i;
                    eri_abic_[a*o_*v_*v_+b*o_*v_+c*o_+i] = eri[aibc] - eri[acbi];
                }
            }
        }
    }

    free(eri);

    // transform oeis

    // core hamiltonian
    double ** hp = H_->pointer();
    double ** fp = F_->pointer();
    for (size_t mu = 0; mu < ns; mu++) {
        for (size_t p = 0; p < n; p++) {
            double dum = 0.0;
            for (size_t nu = 0; nu < ns; nu++) {
                dum += CL_p[nu][p] * hp[nu][mu];
            }
            oei_tmp[p * ns + mu] = dum;
        }
    }

    for (size_t p = 0; p < n; p++) {
        for (size_t q = 0; q < n; q++) {
            double dum = 0.0;
            for (size_t nu = 0; nu < ns; nu++) {
                dum += CR_p[nu][q] * oei_tmp[p * ns + nu];
            }
            fp[p][q] = dum;
        }
    }


    // calculate t1 contribution to correlation energy
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
    for (size_t i = 0; i < o_; i++) {
        ec += 0.5 * fp[i][i];
    }
    ec += enuc_ - reference_wavefunction_->energy();

    // update orbital energies
    for (size_t i = 0; i < n; i++) {
        epsilon_[i] = fp[i][i];
    }
    free(Qmo);

    // now that we've copied the diagonals of the fock matrix, add extra
    // terms introduced by the cavity to F_

    // transform extra oeis introduced by cavity
    std::shared_ptr<Matrix> oe_tmp (new Matrix(ns,ns));
    double ** oe_tmp_p = oe_tmp->pointer();
    double ** oe_cavity_terms_p = oe_cavity_terms_->pointer();

    for (size_t mu = 0; mu < ns; mu++) {
        for (size_t p = 0; p < n; p++) {
            double dum = 0.0;
            for (size_t nu = 0; nu < ns; nu++) {
                dum += CL_p[nu][p] * oe_cavity_terms_p[nu][mu];
            }
            oei_tmp[p * ns + mu] = dum;
        }
    }
    for (size_t p = 0; p < n; p++) {
        for (size_t q = 0; q < n; q++) {
            double dum = 0.0;
            for (size_t nu = 0; nu < ns; nu++) {
                dum += CR_p[nu][q] * oei_tmp[p * ns + nu];
            }
            oe_tmp_p[p][q] = dum;
        }
    }
    F_->add(oe_tmp);

    free(oei_tmp);

    // return t1 contribution to correlation energy
    return ec;    
}

void PolaritonicUCCSD::build_mo_eris() {

    // read Qso from disk
    size_t n  = 2L*(size_t)nmo_;
    size_t ns = 2L*(size_t)nso_;

    double * tmp = (double*)malloc(nQ_*ns*ns*sizeof(double));
    memset((void*)tmp,'\0',nQ_*n*n*sizeof(double));

    double * Qmo = (double*)malloc(nQ_*n*ns*sizeof(double));
    memset((void*)Qmo,'\0',nQ_*n*ns*sizeof(double));

    auto psio = std::make_shared<PSIO>();

    psio->open(PSIF_DCC_QSO, PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_QSO, "Qso CC", (char*)Qmo, nQ_ * nso_ * nso_ * sizeof(double));
    psio->close(PSIF_DCC_QSO, 1);

    for (size_t Q = 0; Q < nQ_; Q++) {
        for (size_t mu = 0; mu < nso_; mu++) {
            for (size_t nu = 0; nu < nso_; nu++) {
                tmp[Q*ns*ns+(mu     )*ns+(nu     )] = Qmo[Q*nso_*nso_+mu*nso_+nu];
                tmp[Q*ns*ns+(mu+nso_)*ns+(nu+nso_)] = Qmo[Q*nso_*nso_+mu*nso_+nu];
            }
        }
    }

    // AO->MO transformation

    // I(Q,mu,p) = C(nu,p) Qso(Q,mu,nu)
    F_DGEMM('n','n',n,ns*nQ_,ns,1.0,&(C_->pointer()[0][0]),n,tmp,ns,0.0,Qmo,n);
    for (size_t Q = 0; Q < nQ_; Q++) {
        for (size_t p = 0; p < n; p++) {
            for (size_t mu = 0; mu < ns; mu++) {
                tmp[Q*n*ns+p*ns+mu] = Qmo[Q*n*ns+mu*n+p];
            }
        }
    }
    // Qmo(Q,p,q) = C(mu,q) I(Q,p,mu)
    F_DGEMM('n','n',n,n*nQ_,ns,1.0,&(C_->pointer()[0][0]),n,tmp,ns,0.0,Qmo,n);

    free(tmp);

    double * eri = (double*)malloc(n*n*n*n*sizeof(double));
    memset((void*)eri,'\0',n*n*n*n*sizeof(double));

    // (pq|rs) = Qmo(Q,rs) Qmo(Q,pq)
    F_DGEMM('n','t',n*n,n*n,nQ_,1.0,Qmo,n*n,Qmo,n*n,0.0,eri,n*n);

    free(Qmo);

    // unpack different classes of eris

    // <ij||kl>
    eri_ijkl_ = (double*)malloc(o_*o_*o_*o_*sizeof(double));
    memset((void*)eri_ijkl_,'\0',o_*o_*o_*o_*sizeof(double));

    for (size_t i = 0; i < o_; i++) {
        for (size_t j = 0; j < o_; j++) {
            for (size_t k = 0; k < o_; k++) {
                for (size_t l = 0; l < o_; l++) {
                    size_t ikjl = i*n*n*n+k*n*n+j*n+l;
                    size_t iljk = i*n*n*n+l*n*n+j*n+k;
                    eri_ijkl_[i*o_*o_*o_+j*o_*o_+k*o_+l] = eri[ikjl] - eri[iljk];
                }
            }
        }
    }

    // <ab||cd>
    eri_abcd_ = (double*)malloc(v_*v_*v_*v_*sizeof(double));
    memset((void*)eri_abcd_,'\0',v_*v_*v_*v_*sizeof(double));

    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t c = 0; c < v_; c++) {
                for (size_t d = 0; d < v_; d++) {
                    size_t acbd = (a+o_)*n*n*n+(c+o_)*n*n+(b+o_)*n+(d+o_);
                    size_t adbc = (a+o_)*n*n*n+(d+o_)*n*n+(b+o_)*n+(c+o_);
                    eri_abcd_[a*v_*v_*v_+b*v_*v_+c*v_+d] = eri[acbd] - eri[adbc];
                }
            }
        }
    }

    // <ij||ab>
    eri_ijab_ = (double*)malloc(o_*o_*v_*v_*sizeof(double));
    memset((void*)eri_ijab_,'\0',o_*o_*v_*v_*sizeof(double));

    for (size_t i = 0; i < o_; i++) {
        for (size_t j = 0; j < o_; j++) {
            for (size_t a = 0; a < v_; a++) {
                for (size_t b = 0; b < v_; b++) {
                    size_t iajb = i*n*n*n+(a+o_)*n*n+j*n+(b+o_);
                    size_t ibja = i*n*n*n+(b+o_)*n*n+j*n+(a+o_);
                    eri_ijab_[i*o_*v_*v_+j*v_*v_+a*v_+b] = eri[iajb] - eri[ibja];
                }
            }
        }
    }

    // <ab||ij>
    eri_abij_ = (double*)malloc(o_*o_*v_*v_*sizeof(double));
    memset((void*)eri_abij_,'\0',o_*o_*v_*v_*sizeof(double));

    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    size_t aibj = (a+o_)*n*n*n+i*n*n+(b+o_)*n+j;
                    size_t ajbi = (a+o_)*n*n*n+j*n*n+(b+o_)*n+i;
                    eri_abij_[a*o_*o_*v_+b*o_*o_+i*o_+j] = eri[aibj] - eri[ajbi];
                }
            }
        }
    }

    // <ia||jb>
    eri_iajb_ = (double*)malloc(o_*o_*v_*v_*sizeof(double));
    memset((void*)eri_iajb_,'\0',o_*o_*v_*v_*sizeof(double));
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    size_t ijab = i*n*n*n+j*n*n+(a+o_)*n+(b+o_);
                    size_t ibaj = i*n*n*n+(b+o_)*n*n+(a+o_)*n+j;
                    eri_iajb_[i*o_*v_*v_+a*o_*v_+j*v_+b] = eri[ijab] - eri[ibaj];
                }
            }
        }
    }

    // <ia||jk>
    eri_iajk_ = (double*)malloc(o_*o_*o_*v_*sizeof(double));
    memset((void*)eri_iajk_,'\0',o_*o_*o_*v_*sizeof(double));
    for (size_t i = 0; i < o_; i++) {
        for (size_t a = 0; a < v_; a++) {
            for (size_t j = 0; j < o_; j++) {
                for (size_t k = 0; k < o_; k++) {
                    size_t ijak = i*n*n*n+j*n*n+(a+o_)*n+k;
                    size_t ikaj = i*n*n*n+k*n*n+(a+o_)*n+j;
                    eri_iajk_[i*o_*o_*v_+a*o_*o_+j*o_+k] = eri[ijak] - eri[ikaj];
                }
            }
        }
    }

    // <jk||ia>
    eri_jkia_ = (double*)malloc(o_*o_*o_*v_*sizeof(double));
    memset((void*)eri_jkia_,'\0',o_*o_*o_*v_*sizeof(double));
    for (size_t j = 0; j < o_; j++) {
        for (size_t k = 0; k < o_; k++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t a = 0; a < v_; a++) {
                    size_t jika = j*n*n*n+i*n*n+k*n+(a+o_);
                    size_t jaki = j*n*n*n+(a+o_)*n*n+k*n+i;
                    eri_jkia_[j*o_*o_*v_+k*o_*v_+i*v_+a] = eri[jika] - eri[jaki];
                }
            }
        }
    }

    // <ai||bc>
    eri_aibc_ = (double*)malloc(o_*v_*v_*v_*sizeof(double));
    memset((void*)eri_aibc_,'\0',o_*v_*v_*v_*sizeof(double));

    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t c = 0; c < v_; c++) {
                for (size_t i = 0; i < o_; i++) {
                    size_t abic = (a+o_)*n*n*n+(b+o_)*n*n+i*n+(c+o_);
                    size_t acib = (a+o_)*n*n*n+(c+o_)*n*n+i*n+(b+o_);
                    eri_aibc_[a*o_*v_*v_+i*v_*v_+b*v_+c] = eri[abic] - eri[acib];
                }
            }
        }
    }

    // <ab||ic>
    eri_abic_ = (double*)malloc(o_*v_*v_*v_*sizeof(double));
    memset((void*)eri_abic_,'\0',o_*v_*v_*v_*sizeof(double));
    for (size_t a = 0; a < v_; a++) {
        for (size_t b = 0; b < v_; b++) {
            for (size_t c = 0; c < v_; c++) {
                for (size_t i = 0; i < o_; i++) {
                    size_t aibc = (a+o_)*n*n*n+i*n*n+(b+o_)*n+(c+o_);
                    size_t acbi = (a+o_)*n*n*n+(c+o_)*n*n+(b+o_)*n+i;
                    eri_abic_[a*o_*v_*v_+b*o_*v_+c*o_+i] = eri[aibc] - eri[acbi];
                }
            }
        }
    }

    free(eri);
}


double PolaritonicUCCSD::compute_energy() {

    if ( include_u0_ || include_u1_ || include_u2_ ) {

        throw PsiException("polaritonic UCCSD is not yet working",__FILE__,__LINE__);

    }

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
    include_u0_ = false;
    include_u1_ = false;
    include_u2_ = false;

    double ec = cc_iterations();

/*
    // CCSD iterations with photon
    include_u0_ = true;
    include_u1_ = true;
    include_u2_ = true;
    diis->restart();

    ec = cc_iterations();
*/

    //outfile->Printf("    * Polaritonic UCCSD total energy: %20.12lf\n",energy_ + ec);
    outfile->Printf("    * UCCSD total energy: %20.12lf\n",energy_ + ec);

    // print cavity properties
    //if ( n_photon_states_ > 1 ) {
    //    print_cavity_properties_ = true;
    //    build_cavity_hamiltonian();
    //    print_cavity_properties_ = false;
    //}
    
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

    double ** dp = Dipole_z_->pointer();
    double ** fp = F_->pointer();

    // - 1.00000 d+(e,m)
    for (size_t e = 0; e < v_; e++) {
        for (size_t m = 0; m < o_; m++) {
            ru1_[e*v_+m] = -dp[e+o_][m];
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
    for (size_t e = 0; e < v_; e++) {
        for (size_t a = 0; a < v_; a++) {
            tmp1_[e*v_+a] = fp[e+o_][a+o_];
        }
    }
    F_DGEMM('n','n',o_,v_,v_,1.0,u1_,o_,tmp1_,v_,1.0,ru1_,o_);
  
    // - 1.00000 <i,e||m,a> u1(a,i)
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
    for (size_t e = 0; e < v_; e++) {
        for (size_t m = 0; m < o_; m++) {
            double dum = 0.0;
            for (size_t a = 0; a < v_; a++) {
                for (size_t i = 0; i < o_; i++) {
                    dum += t2_[a*o_*o_*v_+e*o_*o_+m*o_+i] * dp[i][a+o_];
                }
            }
            ru1_[e*o_+m] -= dum;
        }
    }

    // + 1.00000 <i,j||a,b> t2(b,e,j,m) u1(a,i)
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
    for (size_t i = 0; i < o_; i++) {
        for (size_t a = 0; a < v_; a++) {
            tmp1_[i*v_+a] = dp[i][a+o_];
        }
    }
    // I(i,m) = u1_(a,m) d-(i,a)
    F_DGEMM('n','n',o_,o_,v_,1.0,u1_,o_,tmp1_,v_,0.0,tmp2_,o_);
    // r(e,m) = I(i,m) u1(e,i)
    F_DGEMM('n','n',o_,v_,o_,1.0,tmp2_,o_,u1_,o_,1.0,ru1_,o_);

    if ( include_u2_ ) {

        // + 1.00000 F(i,a) u2(a,e,i,m)
        // - 1.00000 d-(i,a) u2(a,e,i,m) u0
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
        for (size_t e = 0; e < v_; e++) {
            for (size_t i = 0; i < o_; i++) {
                for (size_t j = 0; j < o_; j++) {
                    for (size_t a = 0; a < v_; a++) {
                        tmp1_[e*o_*o_*v_+i*o_*v_+j*v_+a] = u2_[a*o_*o_*v_+e*o_*o_+i*o_+j];
                    }
                }
            }
        }
        // v'(i,j,a,m) = <i,j||m,a>
        for (size_t i = 0; i < o_; i++) {
            for (size_t j = 0; j < o_; j++) {
                for (size_t a = 0; a < v_; a++) {
                    for (size_t m = 0; m < o_; m++) {
                        tmp2_[i*o_*o_*v_+j*o_*v_+a*o_+m] = eri_jkia_[i*o_*o_*v_+j*o_*v_+m*v_+a];
                    }
                }
            }
        }
        // r(e,m) = v'(i,j,a,m) u'(e,i,j,a)
        F_DGEMM('n','n',o_,v_,o_*o_*v_,0.5,tmp2_,o_,tmp1_,o_*o_*v_,1.0,ru1_,o_);

// TODO: refactor ov^3 integral terms

        // - 0.50000 <e,i||a,b> u2(a,b,i,m)
        // u'(i,a,b,m) = u2(a,b,i,m)
        for (size_t i = 0; i < o_; i++) {
            for (size_t a = 0; a < v_; a++) {
                for (size_t b = 0; b < v_; b++) {
                    for (size_t m = 0; m < o_; m++) {
                        tmp1_[i*o_*v_*v_+a*o_*v_+b*o_+m] = u2_[a*o_*o_*v_+b*o_*o_+i*o_+m];
                    }
                }
            }
        }
        // r(e,m) = u'(i,a,b,m) <e,i||a,b>
        F_DGEMM('n','n',o_,v_,o_*v_*v_,-0.5,tmp1_,o_,eri_aibc_,o_*v_*v_,1.0,ru1_,o_);

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
    double ** dp = Dipole_z_->pointer();

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
            tmp1_[e*v_+a] = fp[(e+o_)][(a+o_)];
        }
    }
    // r'(e,f,n,m) =  u2(a,f,m,n) F(e,a)
    F_DGEMM('n','n',o_*o_*v_,v_,v_,1.0,u2_,o_*o_*v_,tmp1_,v_,0.0,tmp2_,o_*o_*v_);
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

    // + 0.5 <i,j||m,n> u2(e,f,i,j)
    F_DGEMM('n','n',o_*o_,v_*v_,o_*o_,0.5,eri_ijkl_,o_*o_,u2_,o_*o_,1.0,ru2_,o_*o_);

    // + 0.5 <e,f||a,b> u2(a,b,m,n)
    F_DGEMM('n','n',o_*o_,v_*v_,v_*v_,0.5,u2_,o_*o_,eri_abcd_,v_*v_,1.0,ru2_,o_*o_);

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
    F_DGEMM('n','n',o_*o_,o_*o_,v_*v_,1.0,t2_,o_*o_,eri_ijab_,v_*v_,0.0,tmp2_,o_*o_);

    // r(e,f,m,n) = I(i,j,m,n) u2(e,f,i,j)
    F_DGEMM('n','n',o_*o_,v_*v_,o_*o_,0.25,tmp2_,o_*o_,u2_,o_*o_,1.0,ru2_,o_*o_);

    // + 0.25000 <i,j||a,b> u2(a,b,m,n) t2(e,f,i,j)

    // I(i,j,m,n) = u2(a,b,m,n) <i,j||a,b>
    F_DGEMM('n','n',o_*o_,o_*o_,v_*v_,1.0,u2_,o_*o_,eri_ijab_,v_*v_,0.0,tmp2_,o_*o_);

    // r(e,f,m,n) = I(i,j,m,n) t2(e,f,i,j)
    F_DGEMM('n','n',o_*o_,v_*v_,o_*o_,0.25,tmp2_,o_*o_,t2_,o_*o_,1.0,ru2_,o_*o_);

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
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[e*o_*o_*v_+f*o_*o_+n*o_+m];
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[e*o_*o_*v_+f*o_*o_+m*o_+n];
                }
            }
        }
    }

// TODO: refactor ov^3 integral terms

    // - 1.00000 <e,f||n,a> u1(a,m)
    // + 1.00000 <e,f||m,a> u1(a,n)
    //
    // or
    //
    // + P(m,n) u1(a,n) <e,f||m,a> 
    //
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


// - 1.00000 <e,i||a,b> t2(b,f,m,n) u1(a,i)
// + 1.00000 <f,i||a,b> t2(b,e,m,n) u1(a,i)
//
// or
//
// - P(e,f) <e,i||a,b> t2(b,f,m,n) u1(a,i)
    // I(e,b) = <e,i||a,b> u1(a,i)
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t b = 0; b < v_; b++) {
            double dum = 0.0;
            for (size_t i = 0; i < o_; i++) {
                for (size_t a = 0; a < v_; a++) {
                    dum += eri_aibc_[e*o_*v_*v_+i*v_*v_+a*v_+b] * u1_[a*o_+i];
                }
            }
            tmp1_[e*v_+b] = dum;
        }
    }
    // r'(e,f,m,n) = t2(b,f,m,n) I(e,b)
    F_DGEMM('n','n',o_*o_*v_,v_,v_,-1.0,t2_,o_*o_*v_,tmp1_,v_,0.0,tmp2_,o_*o_*v_);
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

    // - 0.50000 <e,i||a,b> t2(a,b,m,n) u1(f,i)
    // + 0.50000 <f,i||a,b> t2(a,b,m,n) u1(e,i)
    //
    // or
    //
    // - 0.5 P(e,f) <e,i||a,b> t2(a,b,m,n) u1(f,i)

    // I(m,n,e,i) = <e,i||a,b> t2(a,b,m,n) 
    F_DGEMM('t','t',o_*v_,o_*o_,v_*v_,1.0,eri_aibc_,v_*v_,t2_,o_*o_,0.0,tmp1_,o_*v_);
    // r'(m,n,e,f) = u1(f,i) I(m,n,e,i)
    F_DGEMM('t','n',v_,o_*o_*v_,o_,-0.5,u1_,o_,tmp1_,o_,0.0,tmp2_,v_);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v_; e++) {
        for (size_t f = 0; f < v_; f++) {
            for (size_t m = 0; m < o_; m++) {
                for (size_t n = 0; n < o_; n++) {
                    double dum = 0.0;
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] += tmp2_[m*o_*v_*v_+n*v_*v_+e*v_+f];
                    ru2_[e*o_*o_*v_+f*o_*o_+m*o_+n] -= tmp2_[m*o_*v_*v_+n*v_*v_+f*v_+e];
                }
            }
        }
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
            tmp1_[e*v_+a] = dp[(e+o_)][(a+o_)];
        }
    }
    // r'(e,f,n,m) =  t2(a,f,m,n) d+(e,a)
    F_DGEMM('n','n',o_*o_*v_,v_,v_,1.0,t2_,o_*o_*v_,tmp1_,v_,0.0,tmp2_,o_*o_*v_);
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
    // + 1.00000 d-(i,a) u1(f,n) u2(a,e,m,i)
    // - 1.00000 d-(i,a) u1(f,m) u2(a,e,n,i)
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
            tmp1_[e*v_+a] = dp[(e+o_)][(a+o_)];
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

}

// 
// u0 residual
// 
// fully-contracted strings:
// - 1.00000 d+(i,i) 
// + 1.00000 F(i,a) u1(a,i) 
// + 0.25000 <i,j||a,b> u2(a,b,i,j) 
// + 1.00000 u0 w0 
void PolaritonicUCCSD::residual_u0() {

    // TODO: generalize for x,y,z
    double w0 = cavity_frequency_[2];
    double r0 = u0_[0] * w0;

    // - d+(i,i) 
    double ** dp = Dipole_z_->pointer();
    for (size_t i = 0; i < o_; i++) {
        r0 -= dp[i][i];
    }

    if ( include_u1_ ) {

        // F(i,a) u1(a,i) 
        double **fp = F_->pointer();
        for (size_t i = 0; i < o_; i++) {
            for (size_t a = 0; a < v_; a++) {
                r0 += fp[i][a+o_] * u1_[a*o_+i];
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

void PolaritonicUCCSD::residual_t1() {

    memset((void*)rt1_,'\0',o_*v_*sizeof(double));

    double ** fp = F_->pointer();

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

// TODO: refactor ov^3 integral terms

    // r(e,m) = t'(i,a,b,m) <e,i||a,b>
    F_DGEMM('n','n',o_,v_,o_*v_*v_,-0.5,tmp1_,o_,eri_aibc_,o_*v_*v_,1.0,rt1_,o_);

    // cavity terms:

    double ** dp = Dipole_z_->pointer();

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
            dii += dp[i][i];
        }
        C_DAXPY(o_*v_,-dii,u1_,1,rt1_,1);

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


void PolaritonicUCCSD::residual_t2() {

    memset((void*)rt2_,'\0',o_*o_*v_*v_*sizeof(double));

    double ** fp = F_->pointer();

    // <e,f||m,n> 
    C_DCOPY(o_*o_*v_*v_,eri_abij_,1,rt2_,1);

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

    // + 0.5 <e,f||a,b> t2(a,b,m,n)
    F_DGEMM('n','n',o_*o_,v_*v_,v_*v_,0.5,t2_,o_*o_,eri_abcd_,v_*v_,1.0,rt2_,o_*o_);

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
    F_DGEMM('n','n',o_*o_,o_*o_,v_*v_,1.0,t2_,o_*o_,eri_ijab_,v_*v_,0.0,tmp2_,o_*o_);

    // r(e,f,m,n) = I(i,j,m,n) t'(e,f,i,j)
    F_DGEMM('n','n',o_*o_,v_*v_,o_*o_,0.25,tmp2_,o_*o_,t2_,o_*o_,1.0,rt2_,o_*o_);

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

        double ** dp = Dipole_z_->pointer();

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

        double ** dp = Dipole_z_->pointer();

        // - d-(i,i) u2(e,f,m,n)
        double dii = 0.0;
        for (size_t i = 0; i < o_; i++) {
            dii += dp[i][i];
        }
        C_DAXPY(o_*o_*v_*v_,-dii,u2_,1,rt2_,1);

        
        // + P(m,n) d-(i,n) u2(e,f,m,i)
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

        double ** dp = Dipole_z_->pointer();

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
        for (size_t a = 0; a < v_; a++) {
            double da = epsilon_[a+o_];
            for (size_t i = 0; i < o_; i++) {
                double dai = da - epsilon_[i];
                size_t ai = a*o_+i;
                ru1_[ai] /= -dai;
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
    }

    // diis
    C_DAXPY(ccamps_dim_,1.0,residual_,1,ccamps_,1);
    diis->WriteVector(ccamps_);
    diis->WriteErrorVector(residual_);
    diis->Extrapolate(ccamps_);

    //if ( include_u0_ ) {
    //    printf("u0 = %20.12lf\n",u0_[0]);
    //}

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

        // - 1.00000 d-(i,i) u0 
        double ** dp = Dipole_z_->pointer();
        for (size_t i = 0; i < o_; i++) {
            ec -= dp[i][i] * u0_[0];
        }

    }

    if ( include_u1_ ) {

        // - 1.00000 d-(i,a) u1(a,i) 
        double ** dp = Dipole_z_->pointer();
        for (size_t i = 0; i < o_; i++) {
            for (size_t a = 0; a < v_; a++) {
                ec -= dp[i][a+o_] * u1_[a*o_+i];
            }
        }

    }

    return ec;
}

} // End namespaces
