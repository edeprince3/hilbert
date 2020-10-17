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

    same_a_b_orbs_ = false;
    same_a_b_dens_ = false;

    // ensure scf_type df
    if ( options_.get_str("SCF_TYPE") != "DF" ) {
        throw PsiException("polaritonic uhf only works with scf_type df for now",__FILE__,__LINE__);
    }

    // ensure running in c1 symmetry
    if ( reference_wavefunction_->nirrep() > 1 ) {
        throw PsiException("polaritonic uccsd only works with c1 symmetry for now.",__FILE__,__LINE__);
    }

    // include amplitudes for photon transitions?
    include_u0_ = false;
    include_u1_ = false;
    include_u2_ = false;

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
    Dipole_ = (std::shared_ptr<Matrix>)(new Matrix(2*nso_,2*nso_));
    double ** dp = Dipole_->pointer();
    double ** dipole_scaled_sum_p = dipole_scaled_sum_->pointer();
    for (size_t mu = 0; mu < nso_; mu++) {
        for (size_t nu = 0; nu < nso_; nu++) {
            dp[mu][nu]           = dipole_scaled_sum_p[mu][nu];
            dp[mu+nso_][nu+nso_] = dipole_scaled_sum_p[mu][nu];
        }
    }
    Dipole_->transform(C_);

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

    // allocate memory for amplitudes, residual, and temporary buffers

    size_t o = nalpha_ + nbeta_;
    size_t v = (nmo_-nalpha_) + (nmo_-nbeta_);

    ccamps_dim_ = o*o*v*v + o*v;

    if ( include_u0_ ) {
        ccamps_dim_++;
    }
    if ( include_u1_ ) {
        ccamps_dim_ += o*v;
    }
    if ( include_u2_ ) {
        ccamps_dim_ += o*o*v*v;
    }

    ccamps_   = (double*)malloc(ccamps_dim_*sizeof(double));
    residual_ = (double*)malloc(ccamps_dim_*sizeof(double));

    size_t off = 0;
    t2_  = ccamps_   + off; 
    rt2_ = residual_ + off; off += o*o*v*v;

    t1_  = ccamps_   + off; 
    rt1_ = residual_ + off; off += o*v;

    if ( include_u0_ ) {
        u0_    = ccamps_   + off;
        ru0_   = residual_ + off; off++;
    }
    if ( include_u1_ ) {
        u1_    = ccamps_   + off;
        ru1_   = residual_ + off; off += o*v;
    }
    if ( include_u2_ ) {
        u2_    = ccamps_   + off;
        ru2_   = residual_ + off; off += o*o*v*v;
    }

    memset((void*)ccamps_,'\0',ccamps_dim_*sizeof(double));
    memset((void*)residual_,'\0',ccamps_dim_*sizeof(double));

    // temporary storage ... reduce later

    tmp1_ = (double*)malloc(o*o*v*v*sizeof(double));
    tmp2_ = (double*)malloc(o*o*v*v*sizeof(double));
    tmp3_ = (double*)malloc(o*o*v*v*sizeof(double));
    memset((void*)tmp1_,'\0',o*o*v*v*sizeof(double));
    memset((void*)tmp2_,'\0',o*o*v*v*sizeof(double));
    memset((void*)tmp3_,'\0',o*o*v*v*sizeof(double));

    // initialize diis solver
    diis = (std::shared_ptr<DIIS>)(new DIIS(ccamps_dim_));

}

// TODO: this is so wasteful. need to
// 1. fold t1 into 3-index integrals (done!)
// 2. don't build full four-index tensor
// 3. don't build (ac|bd) 

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

    size_t n  = 2L*(size_t)nmo_;
    size_t ns = 2L*(size_t)nso_;
    size_t o  = nalpha_ + nbeta_;
    size_t v  = n - o;

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
        for (size_t a = 0; a < v; a++) {
            double dum = 0.0;
            for (size_t i = 0; i < o; i++) {
                dum += C_p[mu][i] * t1_[a * o + i];
            }
            CL_p[mu][a + o] -= dum;
        }
    }
#pragma omp parallel for schedule(static)
    for (size_t mu = 0; mu < ns; mu++) {
        for (size_t i = 0; i < o; i++) {
            double dum = 0.0;
            for (size_t a = 0; a < v; a++) {
                dum += C_p[mu][a + o] * t1_[a * o + i];
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
    double ** dp = Dipole_->pointer();

    std::shared_ptr<Matrix> dipole_tmp (new Matrix(ns,ns));
    double ** dipole_tmp_p = dipole_tmp->pointer();

    for (size_t mu = 0; mu < ns; mu++) {
        for (size_t p = 0; p < n; p++) {
            double dum = 0.0;
            for (size_t nu = 0; nu < ns; nu++) {
                dum += CL_p[nu][p] * dp[nu][mu];
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
            dipole_tmp_p[p][q] = dum;
        }
    }
    for (size_t p = 0; p < n; p++) {
        for (size_t q = 0; q < n; q++) {
            for (size_t r = 0; r < n; r++) {
                for (size_t s = 0; s < n; s++) {
                    eri[p*n*n*n+q*n*n+r*n+s] += dipole_tmp_p[p][q] * dipole_tmp_p[r][s];
                }
            }
        }
    }

    // unpack different classes of eris

    // <ij||kl>
    memset((void*)eri_ijkl_,'\0',o*o*o*o*sizeof(double));
    for (size_t i = 0; i < o; i++) {
        for (size_t j = 0; j < o; j++) {
            for (size_t k = 0; k < o; k++) {
                for (size_t l = 0; l < o; l++) {
                    size_t ikjl = i*n*n*n+k*n*n+j*n+l;
                    size_t iljk = i*n*n*n+l*n*n+j*n+k;
                    eri_ijkl_[i*o*o*o+j*o*o+k*o+l]  = eri[ikjl] - eri[iljk];
                    
                }
            }
        }
    }

    // <ab||cd>
    memset((void*)eri_abcd_,'\0',v*v*v*v*sizeof(double));
    for (size_t a = 0; a < v; a++) {
        for (size_t b = 0; b < v; b++) {
            for (size_t c = 0; c < v; c++) {
                for (size_t d = 0; d < v; d++) {
                    size_t acbd = (a+o)*n*n*n+(c+o)*n*n+(b+o)*n+(d+o);
                    size_t adbc = (a+o)*n*n*n+(d+o)*n*n+(b+o)*n+(c+o);
                    eri_abcd_[a*v*v*v+b*v*v+c*v+d] = eri[acbd] - eri[adbc];
                }
            }
        }
    }

    // <ij||ab>
    memset((void*)eri_ijab_,'\0',o*o*v*v*sizeof(double));
    for (size_t i = 0; i < o; i++) {
        for (size_t j = 0; j < o; j++) {
            for (size_t a = 0; a < v; a++) {
                for (size_t b = 0; b < v; b++) {
                    size_t iajb = i*n*n*n+(a+o)*n*n+j*n+(b+o);
                    size_t ibja = i*n*n*n+(b+o)*n*n+j*n+(a+o);
                    eri_ijab_[i*o*v*v+j*v*v+a*v+b] = eri[iajb] - eri[ibja];
                }
            }
        }
    }

    // <ab||ij>
    memset((void*)eri_abij_,'\0',o*o*v*v*sizeof(double));
    for (size_t a = 0; a < v; a++) {
        for (size_t b = 0; b < v; b++) {
            for (size_t i = 0; i < o; i++) {
                for (size_t j = 0; j < o; j++) {
                    size_t aibj = (a+o)*n*n*n+i*n*n+(b+o)*n+j;
                    size_t ajbi = (a+o)*n*n*n+j*n*n+(b+o)*n+i;
                    eri_abij_[a*o*o*v+b*o*o+i*o+j] = eri[aibj] - eri[ajbi];
                }
            }
        }
    }

    // <ia||jb>
    memset((void*)eri_iajb_,'\0',o*o*v*v*sizeof(double));
    for (size_t a = 0; a < v; a++) {
        for (size_t b = 0; b < v; b++) {
            for (size_t i = 0; i < o; i++) {
                for (size_t j = 0; j < o; j++) {
                    size_t ijab = i*n*n*n+j*n*n+(a+o)*n+(b+o);
                    size_t ibaj = i*n*n*n+(b+o)*n*n+(a+o)*n+j;
                    eri_iajb_[i*o*v*v+a*o*v+j*v+b] = eri[ijab] - eri[ibaj];
                }
            }
        }
    }

    // <ia||jk>
    memset((void*)eri_iajk_,'\0',o*o*o*v*sizeof(double));
    for (size_t i = 0; i < o; i++) {
        for (size_t a = 0; a < v; a++) {
            for (size_t j = 0; j < o; j++) {
                for (size_t k = 0; k < o; k++) {
                    size_t ijak = i*n*n*n+j*n*n+(a+o)*n+k;
                    size_t ikaj = i*n*n*n+k*n*n+(a+o)*n+j;
                    eri_iajk_[i*o*o*v+a*o*o+j*o+k] = eri[ijak] - eri[ikaj];
                }
            }
        }
    }

    // <jk||ia>
    memset((void*)eri_jkia_,'\0',o*o*o*v*sizeof(double));
    for (size_t j = 0; j < o; j++) {
        for (size_t k = 0; k < o; k++) {
            for (size_t i = 0; i < o; i++) {
                for (size_t a = 0; a < v; a++) {
                    size_t jika = j*n*n*n+i*n*n+k*n+(a+o);
                    size_t jaki = j*n*n*n+(a+o)*n*n+k*n+i;
                    eri_jkia_[j*o*o*v+k*o*v+i*v+a] = eri[jika] - eri[jaki];
                }
            }
        }
    }

    // <ai||bc>
    memset((void*)eri_aibc_,'\0',o*v*v*v*sizeof(double));
    for (size_t a = 0; a < v; a++) {
        for (size_t b = 0; b < v; b++) {
            for (size_t c = 0; c < v; c++) {
                for (size_t i = 0; i < o; i++) {
                    size_t abic = (a+o)*n*n*n+(b+o)*n*n+i*n+(c+o);
                    size_t acib = (a+o)*n*n*n+(c+o)*n*n+i*n+(b+o);
                    eri_aibc_[a*o*v*v+i*v*v+b*v+c] = eri[abic] - eri[acib];
                }
            }
        }
    }

    free(eri);

    // build fock matrix:

    // transform oeis
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
    for (size_t i = 0; i < o; i++) {
        ec += 0.5 * fp[i][i];
    }

    for (size_t Q = 0; Q < nQ_; Q++) {

        // exchange:
        // sum k (q|rk) (q|ks) 
        F_DGEMM('n','n',n,n,o,-1.0,Qmo + Q*n*n, n,Qmo + Q*n*n, n,1.0,&(fp[0][0]),n);

        // coulomb
        // sum k (q|kk) (q|rs)
        double dum = 0.0; 
        for (size_t k = 0; k < o; k++) { 
            dum += Qmo[Q*n*n + k*n + k];
        }
        C_DAXPY(n*n, dum, Qmo + Q*n*n, 1, &(fp[0][0]), 1);
    }
    for (size_t i = 0; i < o; i++) {
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

    size_t o = nalpha_ + nbeta_;
    size_t v = (nmo_-nalpha_) + (nmo_-nbeta_);

    // <ij||kl>
    eri_ijkl_ = (double*)malloc(o*o*o*o*sizeof(double));
    memset((void*)eri_ijkl_,'\0',o*o*o*o*sizeof(double));

    for (size_t i = 0; i < o; i++) {
        for (size_t j = 0; j < o; j++) {
            for (size_t k = 0; k < o; k++) {
                for (size_t l = 0; l < o; l++) {
                    size_t ikjl = i*n*n*n+k*n*n+j*n+l;
                    size_t iljk = i*n*n*n+l*n*n+j*n+k;
                    eri_ijkl_[i*o*o*o+j*o*o+k*o+l] = eri[ikjl] - eri[iljk];
                }
            }
        }
    }

    // <ab||cd>
    eri_abcd_ = (double*)malloc(v*v*v*v*sizeof(double));
    memset((void*)eri_abcd_,'\0',v*v*v*v*sizeof(double));

    for (size_t a = 0; a < v; a++) {
        for (size_t b = 0; b < v; b++) {
            for (size_t c = 0; c < v; c++) {
                for (size_t d = 0; d < v; d++) {
                    size_t acbd = (a+o)*n*n*n+(c+o)*n*n+(b+o)*n+(d+o);
                    size_t adbc = (a+o)*n*n*n+(d+o)*n*n+(b+o)*n+(c+o);
                    eri_abcd_[a*v*v*v+b*v*v+c*v+d] = eri[acbd] - eri[adbc];
                }
            }
        }
    }

    // <ij||ab>
    eri_ijab_ = (double*)malloc(o*o*v*v*sizeof(double));
    memset((void*)eri_ijab_,'\0',o*o*v*v*sizeof(double));

    for (size_t i = 0; i < o; i++) {
        for (size_t j = 0; j < o; j++) {
            for (size_t a = 0; a < v; a++) {
                for (size_t b = 0; b < v; b++) {
                    size_t iajb = i*n*n*n+(a+o)*n*n+j*n+(b+o);
                    size_t ibja = i*n*n*n+(b+o)*n*n+j*n+(a+o);
                    eri_ijab_[i*o*v*v+j*v*v+a*v+b] = eri[iajb] - eri[ibja];
                }
            }
        }
    }

    // <ab||ij>
    eri_abij_ = (double*)malloc(o*o*v*v*sizeof(double));
    memset((void*)eri_abij_,'\0',o*o*v*v*sizeof(double));

    for (size_t a = 0; a < v; a++) {
        for (size_t b = 0; b < v; b++) {
            for (size_t i = 0; i < o; i++) {
                for (size_t j = 0; j < o; j++) {
                    size_t aibj = (a+o)*n*n*n+i*n*n+(b+o)*n+j;
                    size_t ajbi = (a+o)*n*n*n+j*n*n+(b+o)*n+i;
                    eri_abij_[a*o*o*v+b*o*o+i*o+j] = eri[aibj] - eri[ajbi];
                }
            }
        }
    }

    // <ia||jb>
    eri_iajb_ = (double*)malloc(o*o*v*v*sizeof(double));
    memset((void*)eri_iajb_,'\0',o*o*v*v*sizeof(double));
    for (size_t a = 0; a < v; a++) {
        for (size_t b = 0; b < v; b++) {
            for (size_t i = 0; i < o; i++) {
                for (size_t j = 0; j < o; j++) {
                    size_t ijab = i*n*n*n+j*n*n+(a+o)*n+(b+o);
                    size_t ibaj = i*n*n*n+(b+o)*n*n+(a+o)*n+j;
                    eri_iajb_[i*o*v*v+a*o*v+j*v+b] = eri[ijab] - eri[ibaj];
                }
            }
        }
    }

    // <ia||jk>
    eri_iajk_ = (double*)malloc(o*o*o*v*sizeof(double));
    memset((void*)eri_iajk_,'\0',o*o*o*v*sizeof(double));
    for (size_t i = 0; i < o; i++) {
        for (size_t a = 0; a < v; a++) {
            for (size_t j = 0; j < o; j++) {
                for (size_t k = 0; k < o; k++) {
                    size_t ijak = i*n*n*n+j*n*n+(a+o)*n+k;
                    size_t ikaj = i*n*n*n+k*n*n+(a+o)*n+j;
                    eri_iajk_[i*o*o*v+a*o*o+j*o+k] = eri[ijak] - eri[ikaj];
                }
            }
        }
    }

    // <jk||ia>
    eri_jkia_ = (double*)malloc(o*o*o*v*sizeof(double));
    memset((void*)eri_jkia_,'\0',o*o*o*v*sizeof(double));
    for (size_t j = 0; j < o; j++) {
        for (size_t k = 0; k < o; k++) {
            for (size_t i = 0; i < o; i++) {
                for (size_t a = 0; a < v; a++) {
                    size_t jika = j*n*n*n+i*n*n+k*n+(a+o);
                    size_t jaki = j*n*n*n+(a+o)*n*n+k*n+i;
                    eri_jkia_[j*o*o*v+k*o*v+i*v+a] = eri[jika] - eri[jaki];
                }
            }
        }
    }

    // <ai||bc>
    eri_aibc_ = (double*)malloc(o*v*v*v*sizeof(double));
    memset((void*)eri_aibc_,'\0',o*v*v*v*sizeof(double));

    for (size_t a = 0; a < v; a++) {
        for (size_t b = 0; b < v; b++) {
            for (size_t c = 0; c < v; c++) {
                for (size_t i = 0; i < o; i++) {
                    size_t abic = (a+o)*n*n*n+(b+o)*n*n+i*n+(c+o);
                    size_t acib = (a+o)*n*n*n+(c+o)*n*n+i*n+(b+o);
                    eri_aibc_[a*o*v*v+i*v*v+b*v+c] = eri[abic] - eri[acib];
                }
            }
        }
    }

    free(eri);
}


double PolaritonicUCCSD::compute_energy() {

    // grab some input options_
    double e_convergence = options_.get_double("E_CONVERGENCE");
    double d_convergence = options_.get_double("D_CONVERGENCE");
    size_t maxiter          = options_.get_int("MAXITER");

    outfile->Printf("\n");
    outfile->Printf("    No. basis functions:            %5i\n",nso_);
    outfile->Printf("    No. auxiliary basis functions:  %5i\n",nQ_);
    outfile->Printf("    No. alpha electrons:            %5i\n",nalpha_);
    outfile->Printf("    No. beta electrons:             %5i\n",nbeta_);
    outfile->Printf("    e_convergence:             %10.3le\n",e_convergence);
    outfile->Printf("    d_convergence:             %10.3le\n",d_convergence);
    outfile->Printf("    maxiter:                        %5i\n",maxiter);
    outfile->Printf("\n");
    outfile->Printf("\n");

    // CCSD iterations

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

    }while(fabs(dele) > e_convergence || tnorm > d_convergence );

    if ( iter > maxiter ) {
        throw PsiException("Maximum number of iterations exceeded!",__FILE__,__LINE__);
    }

    outfile->Printf("\n");
    outfile->Printf("    CCSD iterations converged!\n");
    outfile->Printf("\n");

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

// evaluate total residual
void PolaritonicUCCSD::residual() {

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
void PolaritonicUCCSD::residual_u1() {
}

// 
// u2 residual
// 
void PolaritonicUCCSD::residual_u2() {
    throw PsiException("residual_u2() not yet implemented",__FILE__,__LINE__);
}

// 
// u0 residual
// 
// fully-contracted strings:
//     + 1.00000 u0 w0 
//     - 1.00000 d+(i,i) 
//     + 1.00000 F(i,a) u1(a,i) 
//     - 1.00000 d+(i,a) t1(a,i) 
//     + 0.25000 <i,j||a,b> u2(a,b,i,j) 
//     + 1.00000 <i,j||a,b> t1(b,j) u1(a,i) 
void PolaritonicUCCSD::residual_u0() {

    size_t o = nalpha_ + nbeta_;
    size_t v = 2L * nmo_ - o;

    double r0 = 0.0;

    // - d+(i,i) 
    double ** dp = Dipole_->pointer();
    for (size_t i = 0; i < o; i++) {
        r0 -= dp[i][i];
    }


    if ( include_u1_ ) {

        // F(i,a) u1(a,i) 
        double **fp = F_->pointer();
        for (size_t i = 0; i < o; i++) {
            for (size_t a = 0; a < v; a++) {
                r0 += fp[i][a+o] * u1_[a*o+i];
            }
        }

        // + <i,j||a,b> t1(b,j) u1(a,i) 
        for (size_t a = 0; a < v; a++) {
            for (size_t b = 0; b < v; b++) {
                for (size_t i = 0; i < o; i++) {
                    for (size_t j = 0; j < o; j++) {
                        r0 += 0.25 * eri_ijab_[i*o*v*v+j*v*v+a*v+b] * u1_[a*o*+i] * t1_[b*o+j];
                    }
                }
            }
        }

    }

    // - d+(i,a) t1(a,i) 
    for (size_t i = 0; i < o; i++) {
        for (size_t a = 0; a < v; a++) {
            r0 -= dp[i][a+o] * t1_[a*o+i];
        }
    }

    if ( include_u2_ ) {

        // + 0.25 <i,j||a,b> u2(a,b,i,j) 
        for (size_t a = 0; a < v; a++) {
            for (size_t b = 0; b < v; b++) {
                for (size_t i = 0; i < o; i++) {
                    for (size_t j = 0; j < o; j++) {
                        r0 += 0.25 * eri_ijab_[i*o*v*v+j*v*v+a*v+b] * u2_[a*o*o*v+b*o*o+i*o+j];
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

void PolaritonicUCCSD::residual_t1() {

    size_t o = nalpha_ + nbeta_;
    size_t v = 2L * nmo_ - o;

    memset((void*)rt1_,'\0',o*v*sizeof(double));

    double ** fp = F_->pointer();

    for (size_t e = 0; e < v; e++) {
        for (size_t m = 0; m < o; m++) {
            rt1_[e*o+m] = fp[(e+o)][m];
        }
    }

    // F(i,a) t2(a,e,i,m) 
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v; e++) {
        for (size_t m = 0; m < o; m++) {
            double dum = 0.0;
            for (size_t i = 0; i < o; i++) {
                for (size_t a = 0; a < v; a++) {
                    dum += fp[i][(a+o)] * t2_[a*o*o*v+e*o*o+i*o+m];
                }
            }
            rt1_[e*o+m] += dum;
        }
    }
    
    // + 0.5 <i,j||m,a> t2(a,e,i,j)

    // v'(m,i,j,a) = <i,j||m,a>
#pragma omp parallel for schedule(static)
    for (size_t m = 0; m < o; m++) {
        for (size_t i = 0; i < o; i++) {
            for (size_t j = 0; j < o; j++) {
                for (size_t a = 0; a < v; a++) {
                    //tmp1_[m*o*o*v+i*o*v+j*v+a] = eri_iajk_[m*o*o*v+a*o*o+i*o+j];
                    tmp1_[m*o*o*v+i*o*v+j*v+a] = eri_jkia_[i*o*o*v+j*o*v+m*v+a];
                }
            }
        }
    }

    // t'(i,j,a,e) = t2(a,e,i,j)
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o; i++) {
        for (size_t j = 0; j < o; j++) {
            for (size_t a = 0; a < v; a++) {
                for (size_t e = 0; e < v; e++) {
                    tmp2_[i*o*v*v+j*v*v+a*v+e] = t2_[a*o*o*v+e*o*o+i*o+j];
                }
            }
        }
    }

    // r(e,m) = v'(m,i,j,a) t'(i,j,a,e)
    F_DGEMM('t','t',o,v,o*o*v,0.5,tmp1_,o*o*v,tmp2_,v,1.0,rt1_,o);

    // - 0.5 <e,i||a,b> t2(a,b,i,m)

    // t'(i,a,b,m) = t2(a,b,i,m)
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o; i++) {
        for (size_t a = 0; a < v; a++) {
            for (size_t b = 0; b < v; b++) {
                for (size_t m = 0; m < o; m++) {
                    tmp1_[i*o*v*v+a*o*v+b*o+m] = t2_[a*o*o*v+b*o*o+i*o+m];
                }   
            }   
        }   
    }

    // r(e,m) = t'(i,a,b,m) <e,i||a,b>
    F_DGEMM('n','n',o,v,o*v*v,-0.5,tmp1_,o,eri_aibc_,o*v*v,1.0,rt1_,o);

    // cavity terms:

    double ** dp = Dipole_->pointer();

    if ( include_u0_ ) {

        // - d-(e,m) u0 
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v; e++) {
            for (size_t m = 0; m < o; m++) {
                rt1_[e*o+m] -= dp[e+o][m] * u0_[0];
            }
        }
    }

    if ( include_u1_ ) {

        // - d-(i,i) u1(e,m) 
        double dii = 0.0;
        for (size_t i = 0; i < o; i++) {
            dii += dp[i][i];
        }
        C_DAXPY(o*v,-dii,u1_,1,rt1_,1);

        // + d-(i,m) u1(e,i) 
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < o; i++) {
            for (size_t m = 0; m < o; m++) {
                tmp1_[i*o+m] = dp[i][m];
            }
        }
        F_DGEMM('n','n',o,v,o,1.0,tmp1_,o,u1_,o,1.0,rt1_,o);

        // - u1(a,m) d-(e,a) 
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v; e++) {
            for (size_t a = 0; a < v; a++) {
                tmp1_[e*v+a] = dp[e][a];
            }
        }
        F_DGEMM('n','n',o,v,v,-1.0,u1_,o,tmp1_,v,1.0,rt1_,o);

    }

    if ( include_u2_ ) {

        // - d-(i,a) u2(a,e,i,m) 
#pragma omp parallel for schedule(static)
        for (size_t e = 0; e < v; e++) {
            for (size_t m = 0; m < o; m++) {
                double dum = 0.0;
                for (size_t i = 0; i < o; i++) {
                    for (size_t a = 0; a < v; a++) {
                        dum -= dp[i][a+o] * u2_[a*o*o*v+e*o*o+i*o+m];
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

void PolaritonicUCCSD::residual_t2() {

    size_t o = nalpha_ + nbeta_;
    size_t v = 2L * nmo_ - o;

    memset((void*)rt2_,'\0',o*o*v*v*sizeof(double));

    double ** fp = F_->pointer();

    // <e,f||m,n> 
    C_DCOPY(o*o*v*v,eri_abij_,1,rt2_,1);

/*
    // + F(i,m) t2(e,f,n,i)
    // - F(i,n) t2(e,f,i,m)
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o; i++) {
        for (size_t m = 0; m < o; m++) {
            tmp1_[i*o+m] = fp[i][m];
        }
    }
    // r'(e,f,n,m) =  F(i,m) t2(e,f,n,i)
    F_DGEMM('n','n',o,o*v*v,o,1.0,tmp1_,o,t2_,o,0.0,tmp2_,o);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v; e++) {
        for (size_t f = 0; f < v; f++) {
            for (size_t m = 0; m < o; m++) {
                for (size_t n = 0; n < o; n++) {
                    rt2_[e*o*o*v+f*o*o+m*o+n] += tmp2_[e*o*o*v+f*o*o+n*o+m];
                    rt2_[e*o*o*v+f*o*o+m*o+n] -= tmp2_[e*o*o*v+f*o*o+m*o+n];
                }
            }
        }
    }
    // + F(e,a) t2(a,f,m,n)
    // - F(f,a) t2(a,e,m,n)
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v; e++) {
        for (size_t a = 0; a < v; a++) {
            tmp1_[e*v+a] = fp[(e+o)][(a+o)];
        }
    }
    // r'(e,f,n,m) =  t2(a,f,m,n) F(e,a)
    F_DGEMM('n','n',o*o*v,v,v,1.0,t2_,o*o*v,tmp1_,v,0.0,tmp2_,o*o*v);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v; e++) {
        for (size_t f = 0; f < v; f++) {
            for (size_t m = 0; m < o; m++) {
                for (size_t n = 0; n < o; n++) {
                    rt2_[e*o*o*v+f*o*o+m*o+n] += tmp2_[e*o*o*v+f*o*o+m*o+n];
                    rt2_[e*o*o*v+f*o*o+m*o+n] -= tmp2_[f*o*o*v+e*o*o+m*o+n];
                }
            }
        }
    }
*/
        


    // + 0.5 <i,j||m,n> t2(e,f,i,j)
    F_DGEMM('n','n',o*o,v*v,o*o,0.5,eri_ijkl_,o*o,t2_,o*o,1.0,residual_,o*o);

    // + 0.5 <e,f||a,b> t2(a,b,m,n)
    F_DGEMM('n','n',o*o,v*v,v*v,0.5,t2_,o*o,eri_abcd_,v*v,1.0,residual_,o*o);

    // - P(e,f) P(m,n) <i,e||n,a> t2(a,f,m,i)

    // t'(a,i,f,m) = t2(a,f,m,i)
#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v; a++) {
        for (size_t i = 0; i < o; i++) {
            for (size_t f = 0; f < v; f++) {
                for (size_t m = 0; m < o; m++) {
                    tmp1_[a*o*o*v+i*o*v+f*o+m] = t2_[a*o*o*v+f*o*o+m*o+i];
                }
            }
        }
    }

    // v'(a,i,e,n) = <ie||na>
#pragma omp parallel for schedule(static)
    for (size_t a = 0; a < v; a++) {
        for (size_t i = 0; i < o; i++) {
            for (size_t e = 0; e < v; e++) {
                for (size_t n = 0; n < o; n++) {
                    tmp2_[a*o*o*v+i*o*v+e*o+n] = eri_iajb_[i*o*v*v+e*o*v+n*v+a];
                }
            }
        }
    }

    // I(e,n,f,m) = t'(a,i,f,m) v'(a,i,e,n) 
    F_DGEMM('n','t',o*v,o*v,o*v,1.0,tmp1_,o*v,tmp2_,o*v,0.0,tmp3_,o*v);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v; e++) {
        for (size_t f = 0; f < v; f++) {
            for (size_t m = 0; m < o; m++) {
                for (size_t n = 0; n < o; n++) {

                    double dum = 0.0;
                    dum -= tmp3_[e*o*o*v+n*o*v+f*o+m];
                    dum += tmp3_[e*o*o*v+m*o*v+f*o+n];
                    dum += tmp3_[f*o*o*v+n*o*v+e*o+m];
                    dum -= tmp3_[f*o*o*v+m*o*v+e*o+n];

                    rt2_[e*o*o*v+f*o*o+m*o+n] += dum;
                }
            }
        }
    }

    // + 0.25 <i,j||a,b> t2(a,b,m,n) t2(e,f,i,j) 

    // I(i,j,m,n) = t2(a,b,m,n) <i,j||a,b>
    F_DGEMM('n','n',o*o,o*o,v*v,1.0,t2_,o*o,eri_ijab_,v*v,0.0,tmp2_,o*o);

    // r(e,f,m,n) = I(i,j,m,n) t'(e,f,i,j)
    F_DGEMM('n','n',o*o,v*v,o*o,0.25,tmp2_,o*o,t2_,o*o,1.0,rt2_,o*o);

    // - 0.5 P(m,n) <i,j||a,b> t2(a,b,n,j) t2(e,f,m,i) 
    //
    // and
    //
    // - F(i,m) t2(e,f,i,n)
    // + F(i,n) t2(e,f,i,m)

    // t'(n,j,a,b) = t2(a,b,n,j)
#pragma omp parallel for schedule(static)
    for (size_t n = 0; n < o; n++) {
        for (size_t j = 0; j < o; j++) {
            for (size_t a = 0; a < v; a++) {
                for (size_t b = 0; b < v; b++) {
                    tmp1_[n*o*v*v+j*v*v+a*v+b] = t2_[a*o*o*v+b*o*o+n*o+j];
                }
            }
        }
    }

    // I(i,n) = t'(n,j,a,b) <i,j||a,b> + 2 F(i,n)
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o; i++) {
        for (size_t n = 0; n < o; n++) {
            tmp2_[i*o+n] = 2.0 * fp[i][n];
        }
    }
    F_DGEMM('t','n',o,o,o*v*v,1.0,tmp1_,o*v*v,eri_ijab_,o*v*v,1.0,tmp2_,o);

    // r(e,f,m,n) = 0.5 I(i,n) t2(e,f,m,i)
    F_DGEMM('n','n',o,o*v*v,o,0.5,tmp2_,o,t2_,o,0.0,tmp1_,o);

    C_DAXPY(o*o*v*v,-1.0,tmp1_,1,rt2_,1);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v; e++) {
        for (size_t f = 0; f < v; f++) {
            for (size_t m = 0; m < o; m++) {
                for (size_t n = 0; n < o; n++) {
                    rt2_[e*o*o*v+f*o*o+m*o+n] += tmp1_[e*o*o*v+f*o*o+n*o+m];
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
    for (size_t a = 0; a < v; a++) {
        for (size_t b = 0; b < v; b++) {
            for (size_t i = 0; i < o; i++) {
                for (size_t j = 0; j < o; j++) {
                    tmp1_[a*o*o*v+b*o*o+i*o+j] = eri_ijab_[i*o*v*v+j*v*v+a*v+b];
                }
            }
        }
    }

    // t'(f,b,i,j) = t2(b,f,i,j)
#pragma omp parallel for schedule(static)
    for (size_t f = 0; f < v; f++) {
        for (size_t b = 0; b < v; b++) {
            for (size_t i = 0; i < o; i++) {
                for (size_t j = 0; j < o; j++) {
                    tmp2_[f*o*o*v+b*o*o+i*o+j] = t2_[b*o*o*v+f*o*o+i*o+j];
                }
            }
        }
    }

    // I(f,a) = v'(a,b,i,j) t'(f,b,i,j) + 2 F(f,a) 
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v; e++) {
        for (size_t a = 0; a < v; a++) {
            tmp3_[e*v+a] = 2.0 * fp[(e+o)][(a+o)];
        }
    }
    F_DGEMM('t','n',v,v,o*o*v,1.0,tmp1_,o*o*v,tmp2_,o*o*v,1.0,tmp3_,v);

    // r(f,e,m,n) = 0.5 t2(a,e,m,n) I(f,a)
    F_DGEMM('n','n',o*o*v,v,v,0.5,t2_,o*o*v,tmp3_,v,0.0,tmp1_,o*o*v);

    C_DAXPY(o*o*v*v,1.0,tmp1_,1,rt2_,1);
#pragma omp parallel for schedule(static)
    for (size_t f = 0; f < v; f++) {
        for (size_t e = 0; e < v; e++) {
            for (size_t m = 0; m < o; m++) {
                for (size_t n = 0; n < o; n++) {
                    rt2_[f*o*o*v+e*o*o+m*o+n] -= tmp1_[e*o*o*v+f*o*o+m*o+n];
                }   
            }   
        }   
    }

    // + P(m,n) <i,j||a,b> t2(a,e,n,j) t2(b,f,m,i) 

    //I(i,b,j,a) = <i,j||a,b>
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < o; i++) {
        for (size_t j = 0; j < o; j++) {
            for (size_t a = 0; a < v; a++) {
                for (size_t b = 0; b < v; b++) {
                    tmp1_[i*o*v*v+b*o*v+j*v+a] = eri_ijab_[i*o*v*v+j*v*v+a*v+b];
                }
            }
        }
    }

    // t'(i,b,m,f) = t2(b,f,m,i)
#pragma omp parallel for schedule(static)
    for (size_t b = 0; b < v; b++) {
        for (size_t f = 0; f < v; f++) {
            for (size_t m = 0; m < o; m++) {
                for (size_t i = 0; i < o; i++) {
                    tmp2_[i*o*v*v+b*o*v+m*v+f] = t2_[b*o*o*v+f*o*o+m*o+i];
                }
            }
        }
    }

    //I''(m,f,j,a) = I(i,b,j,a) t'(i,b,m,f) 
    F_DGEMM('n','t',o*v,o*v,o*v,1.0,tmp1_,o*v,tmp2_,o*v,0.0,tmp3_,o*v);

    // t''(n,e,j,a) = t2(a,e,n,j) 
#pragma omp parallel for schedule(static)
    for (size_t n = 0; n < o; n++) {
        for (size_t e = 0; e < v; e++) {
            for (size_t j = 0; j < o; j++) {
                for (size_t a = 0; a < v; a++) {
                    tmp1_[e*o*o*v+n*o*v+j*v+a] = t2_[a*o*o*v+e*o*o+n*o+j];
                }
            }
        }
    }

    //r'(e,n,m,f) = I''(m,f,j,a) t''(n,e,j,a)
    F_DGEMM('t','n',o*v,o*v,o*v,1.0,tmp3_,o*v,tmp1_,o*v,0.0,tmp2_,o*v);
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < v; e++) {
        for (size_t f = 0; f < v; f++) {
            for (size_t m = 0; m < o; m++) {
                for (size_t n = 0; n < o; n++) {
                    rt2_[e*o*o*v+f*o*o+m*o+n] += tmp2_[e*o*o*v+n*o*v+m*v+f];
                    rt2_[e*o*o*v+f*o*o+m*o+n] -= tmp2_[e*o*o*v+m*o*v+n*v+f];
                }
            }
        }
    }

}

double PolaritonicUCCSD::update_amplitudes() {

    size_t o = nalpha_ + nbeta_;
    size_t v = 2 * nmo_ - o;

    // dt = -residual / eps

    // t2
    for (size_t a = 0; a < v; a++) {
        double da = epsilon_[a+o];
        for (size_t b = 0; b < v; b++) {
            double dab = da + epsilon_[b+o];
            for (size_t i = 0; i < o; i++) {
                double dabi = dab - epsilon_[i];
                for (size_t j = 0; j < o; j++) {
                    double dabij = dabi - epsilon_[j];
                    size_t abij = a*o*o*v+b*o*o+i*o+j;
                    rt2_[abij] /= -dabij;
                }
            }
        }
    }

    // t1
    for (size_t a = 0; a < v; a++) {
        double da = epsilon_[a+o];
        for (size_t i = 0; i < o; i++) {
            double dai = da - epsilon_[i];
            size_t ai = a*o+i;
            rt1_[ai] /= -dai;
        }
    }

    if ( include_u2_ ) {
        // u2
        for (size_t a = 0; a < v; a++) {
            double da = epsilon_[a+o];
            for (size_t b = 0; b < v; b++) {
                double dab = da + epsilon_[b+o];
                for (size_t i = 0; i < o; i++) {
                    double dabi = dab - epsilon_[i];
                    for (size_t j = 0; j < o; j++) {
                        double dabij = dabi - epsilon_[j];
                        size_t abij = a*o*o*v+b*o*o+i*o+j;
                        ru2_[abij] /= -dabij;
                    }
                }
            }
        }
    }
    if ( include_u1_ ) {
        // u1
        for (size_t a = 0; a < v; a++) {
            double da = epsilon_[a+o];
            for (size_t i = 0; i < o; i++) {
                double dai = da - epsilon_[i];
                size_t ai = a*o+i;
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
        ru0_[0] = -ru0_[0] / w0 - u0_[0];
    }

    // diis
    C_DAXPY(ccamps_dim_,1.0,residual_,1,ccamps_,1);
    diis->WriteVector(ccamps_);
    diis->WriteErrorVector(residual_);
    diis->Extrapolate(ccamps_);

    return C_DNRM2(ccamps_dim_,residual_,1);

}

double PolaritonicUCCSD::correlation_energy() {

    size_t o = nalpha_ + nbeta_;
    size_t v = 2 * nmo_ - o;

    double ec = 0.0;

    // + 0.25000 <i,j||a,b> t2(a,b,i,j)
    for (size_t a = 0; a < v; a++) {
        for (size_t b = 0; b < v; b++) {
            for (size_t i = 0; i < o; i++) {
                for (size_t j = 0; j < o; j++) {
                    ec += 0.25 * eri_ijab_[i*o*v*v+j*v*v+a*v+b] * t2_[a*o*o*v+b*o*o+i*o+j];
                }
            }
        }
    }

    return ec;
}

} // End namespaces
